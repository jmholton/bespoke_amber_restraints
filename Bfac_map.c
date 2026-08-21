/* Bfac_map.c  -  store B factors as a spatial field B(x,y,z) instead of per atom
 *                                                            -James Holton 8-20-26
 *
 * Rationale
 * ---------
 * In an MD-derived model the "right" B factor is a property of a POSITION, not of
 * an atom.  A water that wanders onto an ordered site should immediately become
 * sharp; the same water out in a solvent channel should be diffuse.  Tying B to
 * atom identity (Bfac.pdb, matched by atom order) cannot do that: the per-atom B
 * can only be relaxed slowly or the refinement goes unstable.
 *
 * So: keep B in a CCP4 map.  Build it once per optimization stage by letting each
 * atom donate its B to the space around it, then have the structure-factor engine
 * look B up by position at every frame.
 *
 * The field
 * ---------
 * Normalized (Shepard-style) Gaussian inheritance with a background prior:
 *
 *      B(x) = ( SUM_i w_i(x) B_i  +  w0 farB ) / ( SUM_i w_i(x) + w0 )
 *
 *      w_i(x) = occ_i * exp( -r_i(x)^2 / (2 sigma^2) )      r = periodic distance
 *      w0     = exp( -fardist^2 / (2 sigma^2) )
 *
 * Properties:
 *   - at an atom, B -> that atom's B (blended with neighbours within ~sigma)
 *   - more than ~fardist from every atom, B -> farB (999 = "not there")
 *   - everywhere continuous, so trilinear interpolation is well behaved
 *   - w0 is exactly the weight one atom has at r = fardist, which is what makes
 *     fardist the crossover distance rather than an arbitrary knob
 *
 * The map holds FINAL B values: minB/maxB clipping and Bscale/Boffset are applied
 * to the atomic B factors at deposit time, so consumers just interpolate and use.
 * farB is deliberately NOT clipped to maxB - empty space is meant to read as 999.
 *
 * Symmetry
 * --------
 * Sampling is periodic in the map's cell, so an atom anywhere in space - several
 * cells away, or with unwrapped MD coordinates - reads the equivalent voxel by
 * lattice translation.  Space-group symmetry is handled at BUILD time instead:
 * sg=<name> deposits each atom at all of its symmetry-equivalent positions, so
 * the stored field already fills the cell and every consumer needs nothing but a
 * trilinear lookup.  Folding a sample point into an ASU at probe time would need
 * space-group machinery in every consumer, and interpolation on an ASU boundary
 * would need neighbouring voxels that belong to a different symmetry copy.
 *
 * Use sg= only when the input really is one asymmetric unit.  An MD supercell is
 * already expanded and is deliberately NOT symmetric - expanding it again would
 * average symmetry mates together and erase the differences between subcells.
 *
 * Modes
 * -----
 *   build   B-bearing PDB            -> CCP4 map of B(x,y,z)
 *   probe   PDB + CCP4 map           -> same PDB with B taken from the map
 *   stats   CCP4 map                 -> distribution of B in the map
 *
 * Compile:  gcc -O3 -o Bfac_map Bfac_map.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* unit cell                                                          */
/* ------------------------------------------------------------------ */
typedef struct {
    double a, b, c, alpha, beta, gamma;   /* Angstroms, degrees */
} CELL;

/* fractional -> Cartesian metric tensor  G_ij = a_i . a_j
   so that  r^2 = G11 dx^2 + G22 dy^2 + G33 dz^2
                + 2 ( G12 dx dy + G13 dx dz + G23 dy dz )
   for a difference expressed in FRACTIONAL coordinates */
static void cell_metric(CELL c, double *g11, double *g22, double *g33,
                                double *g12, double *g13, double *g23)
{
    double ca = cos(c.alpha * M_PI / 180.0);
    double cb = cos(c.beta  * M_PI / 180.0);
    double cg = cos(c.gamma * M_PI / 180.0);
    *g11 = c.a * c.a;      *g22 = c.b * c.b;      *g33 = c.c * c.c;
    *g12 = c.a * c.b * cg; *g13 = c.a * c.c * cb; *g23 = c.b * c.c * ca;
}

/* Cartesian -> fractional (standard PDB orientation: a along x, b in the xy plane) */
typedef struct { double f00, f01, f02, f11, f12, f22; } FRAC;

static FRAC cell_frac(CELL c)
{
    FRAC f;
    double ca = cos(c.alpha * M_PI / 180.0);
    double cb = cos(c.beta  * M_PI / 180.0);
    double cg = cos(c.gamma * M_PI / 180.0);
    double sg = sin(c.gamma * M_PI / 180.0);
    double V  = sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg);
    f.f00 = 1.0 / c.a;
    f.f01 = -cg / (c.a * sg);
    f.f02 = (ca*cg - cb) / (c.a * V * sg);
    f.f11 = 1.0 / (c.b * sg);
    f.f12 = (cb*cg - ca) / (c.b * V * sg);
    f.f22 = sg / (c.c * V);
    return f;
}

/* inverse metric diagonals: max fractional extent per Angstrom along each axis.
   Used to size the deposit box tightly for non-orthogonal cells. */
static void cell_recip_diag(CELL c, double *rfx, double *rfy, double *rfz)
{
    double g11, g22, g33, g12, g13, g23;
    cell_metric(c, &g11, &g22, &g33, &g12, &g13, &g23);
    double det_ab = g11*g22 - g12*g12;
    double det_G  = g33*det_ab - g22*g13*g13 + 2.0*g12*g13*g23 - g11*g23*g23;
    *rfx = sqrt((g22*g33 - g23*g23) / det_G);
    *rfy = sqrt((g11*g33 - g13*g13) / det_G);
    *rfz = sqrt(det_ab / det_G);
}

static int cell_is_ortho(CELL c)
{
    return (fabs(c.alpha - 90.0) < 1e-4 &&
            fabs(c.beta  - 90.0) < 1e-4 &&
            fabs(c.gamma - 90.0) < 1e-4);
}

/* ------------------------------------------------------------------ */
/* CCP4 map                                                           */
/* ------------------------------------------------------------------ */
/* Only what this program needs: mode-2 (float), full-cell coverage, any axis
 * order on read, X-fast on write.  No symmetry records (ISPG=1, NSYMBT=0):
 * a B field over an MD supercell is not a symmetric object. */
typedef struct {
    int   nx, ny, nz;        /* grid points along cell axes a,b,c            */
    int   ispg;              /* space group number to declare in the header  */
    CELL  cell;
    float *rho;              /* nx*ny*nz, x fastest                          */
} MAP;

#define MAPIDX(m, ix, iy, iz) ((size_t)(iz)*(m)->ny*(m)->nx + (size_t)(iy)*(m)->nx + (ix))

static void map_free(MAP *m) { if (m->rho) free(m->rho); m->rho = NULL; }

static int write_ccp4(const char *fname, MAP *m, const char *label)
{
    FILE *f;
    int   header[256];        /* the whole 1024-byte header: 56 words then 10x80 labels */
    char *labels = (char *)&header[56];
    size_t nvox = (size_t)m->nx * m->ny * m->nz;
    size_t i;
    double sum = 0.0, sumsq = 0.0;
    float  amin =  1e30f, amax = -1e30f;

    for (i = 0; i < nvox; ++i) {
        float v = m->rho[i];
        if (v < amin) amin = v;
        if (v > amax) amax = v;
        sum   += v;
        sumsq += (double)v * v;
    }
    double amean = sum / (double)nvox;
    double arms  = sqrt(sumsq / (double)nvox - amean * amean);

    f = fopen(fname, "wb");
    if (!f) { fprintf(stderr, "ERROR: cannot write %s\n", fname); return 0; }

    memset(header, 0, sizeof(header));
    header[0] = m->nx;                    /* NC   - columns  (X here)        */
    header[1] = m->ny;                    /* NR   - rows     (Y here)        */
    header[2] = m->nz;                    /* NS   - sections (Z here)        */
    header[3] = 2;                        /* MODE - 32-bit float             */
    header[4] = 0; header[5] = 0; header[6] = 0;      /* NCSTART/NRSTART/NSSTART */
    header[7] = m->nx;                    /* NX   - sampling along a         */
    header[8] = m->ny;                    /* NY   - sampling along b         */
    header[9] = m->nz;                    /* NZ   - sampling along c         */
    ((float *)header)[10] = (float)m->cell.a;
    ((float *)header)[11] = (float)m->cell.b;
    ((float *)header)[12] = (float)m->cell.c;
    ((float *)header)[13] = (float)m->cell.alpha;
    ((float *)header)[14] = (float)m->cell.beta;
    ((float *)header)[15] = (float)m->cell.gamma;
    header[16] = 1; header[17] = 2; header[18] = 3;   /* MAPC/MAPR/MAPS = X,Y,Z */
    ((float *)header)[19] = amin;
    ((float *)header)[20] = amax;
    ((float *)header)[21] = (float)amean;
    header[22] = m->ispg ? m->ispg : 1;   /* ISPG                            */
    header[23] = 0;                       /* NSYMBT - no symmetry records    */
    memcpy(&header[52], "MAP ", 4);       /* word 53                         */
    /* machine stamp: little-endian IEEE, same as CCP4 on x86 */
    ((unsigned char *)&header[53])[0] = 0x44;
    ((unsigned char *)&header[53])[1] = 0x41;
    ((unsigned char *)&header[53])[2] = 0x00;
    ((unsigned char *)&header[53])[3] = 0x00;
    ((float *)header)[54] = (float)arms;
    header[55] = 1;                       /* NLABL                           */

    memset(labels, ' ', 800);
    {
        char line[81];
        int n = snprintf(line, sizeof(line), "%s", label);
        if (n > 80) n = 80;
        memcpy(labels, line, n);
    }

    fwrite(header, 4, 256, f);
    fwrite(m->rho, sizeof(float), nvox, f);
    fclose(f);

    printf("wrote %s  %dx%dx%d  B: %.2f - %.2f  mean %.2f  rms %.2f\n",
           fname, m->nx, m->ny, m->nz, amin, amax, amean, arms);
    return 1;
}

static int read_ccp4(const char *fname, MAP *m)
{
    FILE *f;
    int   header[256];
    int   nc, nr, ns, mode, ncstart, nrstart, nsstart;
    int   mx, my, mz, mapc, mapr, maps, nsymbt;
    size_t nvox;
    float *raw;

    f = fopen(fname, "rb");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", fname); return 0; }
    if (fread(header, 4, 256, f) != 256) {
        fprintf(stderr, "ERROR: %s is too short to be a CCP4 map\n", fname);
        fclose(f); return 0;
    }
    if (memcmp(&header[52], "MAP ", 4) != 0)
        fprintf(stderr, "WARNING: no 'MAP ' stamp in %s - assuming CCP4 anyway\n", fname);

    nc = header[0]; nr = header[1]; ns = header[2]; mode = header[3];
    ncstart = header[4]; nrstart = header[5]; nsstart = header[6];
    mx = header[7]; my = header[8]; mz = header[9];
    m->cell.a     = ((float *)header)[10];
    m->cell.b     = ((float *)header)[11];
    m->cell.c     = ((float *)header)[12];
    m->cell.alpha = ((float *)header)[13];
    m->cell.beta  = ((float *)header)[14];
    m->cell.gamma = ((float *)header)[15];
    mapc = header[16]; mapr = header[17]; maps = header[18];
    nsymbt = header[23];

    if (mode != 2) {
        fprintf(stderr, "ERROR: %s has MODE %d; only MODE 2 (float) is supported\n",
                fname, mode);
        fclose(f); return 0;
    }
    if (mapc < 1 || mapc > 3 || mapr < 1 || mapr > 3 || maps < 1 || maps > 3) {
        fprintf(stderr, "WARNING: bad MAPC/MAPR/MAPS (%d %d %d) in %s - assuming 1 2 3\n",
                mapc, mapr, maps, fname);
        mapc = 1; mapr = 2; maps = 3;
    }

    /* the field must tile the whole cell, or position -> value is undefined */
    {
        int ext[3], start[3], samp[3];
        ext[mapc-1] = nc; ext[mapr-1] = nr; ext[maps-1] = ns;
        start[mapc-1] = ncstart; start[mapr-1] = nrstart; start[maps-1] = nsstart;
        samp[0] = mx; samp[1] = my; samp[2] = mz;
        if (ext[0] != samp[0] || ext[1] != samp[1] || ext[2] != samp[2] ||
            start[0] != 0 || start[1] != 0 || start[2] != 0) {
            fprintf(stderr,
                "ERROR: %s covers only part of the cell "
                "(extent %d %d %d of sampling %d %d %d, start %d %d %d).\n"
                "       A B-factor field must span the whole cell.  Fix with:\n"
                "         echo 'xyzlim cell' | mapmask mapin1 %s mapout full.map\n",
                fname, ext[0], ext[1], ext[2], samp[0], samp[1], samp[2],
                start[0], start[1], start[2], fname);
            fclose(f); return 0;
        }
    }

    if (nsymbt) fseek(f, nsymbt, SEEK_CUR);

    nvox = (size_t)nc * nr * ns;
    raw  = (float *)malloc(nvox * sizeof(float));
    if (!raw) { fprintf(stderr, "ERROR: out of memory reading %s\n", fname); fclose(f); return 0; }
    if (fread(raw, sizeof(float), nvox, f) != nvox) {
        fprintf(stderr, "ERROR: %s ended early (wanted %lu voxels)\n",
                fname, (unsigned long)nvox);
        free(raw); fclose(f); return 0;
    }
    fclose(f);

    m->nx = mx; m->ny = my; m->nz = mz;
    m->ispg = header[22];

    if (mapc == 1 && mapr == 2 && maps == 3) {
        m->rho = raw;                      /* already X fast, Y, Z           */
    } else {
        /* transpose whatever order the file used into X-fast */
        int ic, ir, is;
        int idx[3];
        m->rho = (float *)malloc(nvox * sizeof(float));
        if (!m->rho) { fprintf(stderr, "ERROR: out of memory\n"); free(raw); return 0; }
        for (is = 0; is < ns; ++is)
            for (ir = 0; ir < nr; ++ir)
                for (ic = 0; ic < nc; ++ic) {
                    idx[mapc-1] = ic; idx[mapr-1] = ir; idx[maps-1] = is;
                    m->rho[MAPIDX(m, idx[0], idx[1], idx[2])] =
                        raw[((size_t)is*nr + ir)*nc + ic];
                }
        free(raw);
    }
    return 1;
}

/* trilinear interpolation at a fractional coordinate, periodic in all three axes */
static double map_interp(const MAP *m, double xf, double yf, double zf)
{
    double gx, gy, gz, tx, ty, tz;
    int i0, j0, k0, i1, j1, k1;

    xf -= floor(xf); yf -= floor(yf); zf -= floor(zf);
    gx = xf * m->nx; gy = yf * m->ny; gz = zf * m->nz;
    i0 = (int)floor(gx); j0 = (int)floor(gy); k0 = (int)floor(gz);
    tx = gx - i0;       ty = gy - j0;       tz = gz - k0;
    /* floor of a value one ulp below 1.0 can still land on n */
    if (i0 >= m->nx) { i0 = 0; tx = 0.0; }
    if (j0 >= m->ny) { j0 = 0; ty = 0.0; }
    if (k0 >= m->nz) { k0 = 0; tz = 0.0; }
    i1 = (i0 + 1) % m->nx; j1 = (j0 + 1) % m->ny; k1 = (k0 + 1) % m->nz;

    return  (1-tx)*(1-ty)*(1-tz) * m->rho[MAPIDX(m, i0, j0, k0)]
          +    tx *(1-ty)*(1-tz) * m->rho[MAPIDX(m, i1, j0, k0)]
          + (1-tx)*   ty *(1-tz) * m->rho[MAPIDX(m, i0, j1, k0)]
          +    tx *   ty *(1-tz) * m->rho[MAPIDX(m, i1, j1, k0)]
          + (1-tx)*(1-ty)*   tz  * m->rho[MAPIDX(m, i0, j0, k1)]
          +    tx *(1-ty)*   tz  * m->rho[MAPIDX(m, i1, j0, k1)]
          + (1-tx)*   ty *   tz  * m->rho[MAPIDX(m, i0, j1, k1)]
          +    tx *   ty *   tz  * m->rho[MAPIDX(m, i1, j1, k1)];
}

/* ------------------------------------------------------------------ */
/* symmetry operators                                                 */
/* ------------------------------------------------------------------ */
#define MAXOPS 200

typedef struct {
    double r[3][3];          /* rotation, acting on FRACTIONAL coordinates */
    double t[3];             /* translation, fractional                   */
} SYMOP;

/* Parse one triplet, e.g. "X,Y,Z"  "1/2-X,-Y,1/2+Z"  "-Y,X-Y,Z+1/3".
   Coefficients are only 0 and +/-1 and translations are rationals, which is all
   symop.lib ever contains, but the scan is general enough not to care. */
static int parse_symop(const char *str, SYMOP *op)
{
    int row = 0;
    const char *p = str;
    double sign = 1.0;

    memset(op, 0, sizeof(*op));
    while (*p && row < 3) {
        if (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') { ++p; continue; }
        if (*p == ',') { ++row; sign = 1.0; ++p; continue; }
        if (*p == '+') { sign =  1.0; ++p; continue; }
        if (*p == '-') { sign = -1.0; ++p; continue; }
        if (*p == 'X' || *p == 'x') { op->r[row][0] += sign; ++p; continue; }
        if (*p == 'Y' || *p == 'y') { op->r[row][1] += sign; ++p; continue; }
        if (*p == 'Z' || *p == 'z') { op->r[row][2] += sign; ++p; continue; }
        if (*p >= '0' && *p <= '9') {
            char *end;
            double v = strtod(p, &end);
            p = end;
            if (*p == '/') {
                double d = strtod(p + 1, &end);
                p = end;
                if (d == 0.0) return 0;
                v /= d;
            }
            op->t[row] += sign * v;
            continue;
        }
        return 0;              /* something unexpected */
    }
    return (row == 2);         /* exactly three components */
}

/* "X,Y,Z ; -X,-Y,Z ; ..."  (';' or newline separated) */
static int read_symops_str(const char *str, SYMOP *ops, int maxops)
{
    char *buf = strdup(str), *save = buf, *tok;
    int n = 0;
    if (!buf) return 0;
    for (tok = strtok(buf, ";\n"); tok && n < maxops; tok = strtok(NULL, ";\n")) {
        int blank = 1;
        const char *q;
        for (q = tok; *q; ++q) if (*q != ' ' && *q != '\t') { blank = 0; break; }
        if (blank) continue;
        if (!parse_symop(tok, &ops[n])) {
            fprintf(stderr, "ERROR: cannot parse symmetry operator '%s'\n", tok);
            free(save); return 0;
        }
        ++n;
    }
    free(save);
    return n;
}

/* Look a space group up in CCP4's symop.lib.  Header lines are
     number nsym nsymp SHORTNAME PGname SYSTEM 'full name'
   followed by nsym operator lines.  Matches the short name or the number. */
static int read_symop_lib(const char *libpath, const char *sgname,
                          SYMOP *ops, int maxops, int *ispg)
{
    FILE *f;
    char line[512];
    int  want_num = atoi(sgname), n = 0, nsym = 0, insg = 0;

    f = fopen(libpath, "r");
    if (!f) {
        fprintf(stderr, "ERROR: cannot read %s\n"
                "       set CLIBD, or pass symoplib=<path>, or give the operators "
                "directly with symops=\n", libpath);
        return 0;
    }
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != ' ' && line[0] != '\t' && !strchr(line, ',')) {
            int num, ns, nsp;
            char name[64];
            if (insg) break;                       /* finished our block */
            if (sscanf(line, "%d %d %d %63s", &num, &ns, &nsp, name) == 4) {
                if (num >= 500) continue;          /* alternative settings */
                if (num == want_num || !strcasecmp(name, sgname)) {
                    insg = 1; nsym = ns; *ispg = num;
                }
            }
            continue;
        }
        if (!insg) continue;
        if (!strchr(line, ',')) break;
        if (n >= maxops || n >= nsym) break;
        if (!parse_symop(line, &ops[n])) {
            fprintf(stderr, "ERROR: cannot parse '%s' from %s\n", line, libpath);
            fclose(f); return 0;
        }
        ++n;
    }
    fclose(f);
    if (!insg) {
        fprintf(stderr, "ERROR: space group '%s' not found in %s\n", sgname, libpath);
        return 0;
    }
    if (n != nsym)
        fprintf(stderr, "WARNING: %s lists %d operators for %s but %d were read\n",
                libpath, nsym, sgname, n);
    return n;
}

/* ------------------------------------------------------------------ */
/* PDB                                                               */
/* ------------------------------------------------------------------ */
typedef struct {
    double x, y, z, occ, B;
} ATOM;

/* order by position, so exact duplicates end up adjacent */
static int atom_pos_cmp(const void *va, const void *vb)
{
    const ATOM *a = (const ATOM *)va, *b = (const ATOM *)vb;
    if (a->x != b->x) return (a->x < b->x) ? -1 : 1;
    if (a->y != b->y) return (a->y < b->y) ? -1 : 1;
    if (a->z != b->z) return (a->z < b->z) ? -1 : 1;
    return 0;
}

static double pdb_field(const char *line, int start, int len)
{
    char buf[32];
    int  i, n = (int)strlen(line);
    if (len > 31) len = 31;
    for (i = 0; i < len; ++i) buf[i] = (start + i < n) ? line[start + i] : ' ';
    buf[len] = '\0';
    return atof(buf);
}

/* read CRYST1 + ATOM/HETATM.  atoms==NULL just fetches the cell. */
static int read_pdb(const char *fname, CELL *cell, ATOM **atoms, int *natoms)
{
    FILE *f;
    char line[512];
    int  nmax = 65536, n = 0, got_cryst1 = 0;
    ATOM *a = NULL;

    f = fopen(fname, "r");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", fname); return 0; }
    if (atoms) {
        a = (ATOM *)malloc(nmax * sizeof(ATOM));
        if (!a) { fprintf(stderr, "ERROR: out of memory\n"); fclose(f); return 0; }
    }

    while (fgets(line, sizeof(line), f)) {
        if (!strncmp(line, "CRYST1", 6)) {
            cell->a     = pdb_field(line,  6, 9);
            cell->b     = pdb_field(line, 15, 9);
            cell->c     = pdb_field(line, 24, 9);
            cell->alpha = pdb_field(line, 33, 7);
            cell->beta  = pdb_field(line, 40, 7);
            cell->gamma = pdb_field(line, 47, 7);
            got_cryst1  = 1;
            continue;
        }
        if (strncmp(line, "ATOM  ", 6) && strncmp(line, "HETATM", 6)) continue;
        if (!atoms) { ++n; continue; }
        if (n == nmax) {
            nmax *= 2;
            a = (ATOM *)realloc(a, nmax * sizeof(ATOM));
            if (!a) { fprintf(stderr, "ERROR: out of memory at %d atoms\n", n); fclose(f); return 0; }
        }
        a[n].x   = pdb_field(line, 30, 8);
        a[n].y   = pdb_field(line, 38, 8);
        a[n].z   = pdb_field(line, 46, 8);
        a[n].occ = pdb_field(line, 54, 6);
        a[n].B   = pdb_field(line, 60, 6);
        ++n;
    }
    fclose(f);

    if (atoms) *atoms = a;
    *natoms = n;
    return got_cryst1 ? 1 : -1;   /* -1: atoms fine, but no cell in the file */
}

/* ------------------------------------------------------------------ */
/* build                                                              */
/* ------------------------------------------------------------------ */
static int good_grid_size(int n)
{
    /* next 2,3,5-smooth number >= n, so the map stays usable by FFT-based
       CCP4 tools (mapmask, map arithmetic) without regridding */
    int i2, i3, i5, best = n * 10;
    for (i2 = 1; i2 <= n * 2; i2 *= 2)
        for (i3 = i2; i3 <= n * 2; i3 *= 3)
            for (i5 = i3; i5 <= n * 2; i5 *= 5)
                if (i5 >= n && i5 < best) best = i5;
    return best;
}

static int do_build(const char *pdbfile, const char *outmap,
                    CELL cell_override, int have_cell, double super_mult[3],
                    double sigma, double gridspacing, int nxyz[3],
                    double farB, double fardist,
                    double minB, double maxB, double Bscale, double Boffset,
                    int useocc, double wtol,
                    const SYMOP *ops, int nops, int ispg)
{
    CELL   cell;
    ATOM  *atoms = NULL;
    int    natoms, i, cellstat;
    MAP    m;
    float *num = NULL, *den = NULL;
    size_t nvox, v;
    double g11, g22, g33, g12, g13, g23;
    double rfx, rfy, rfz, rcut, w0, twosig2;
    FRAC   fr;
    int    ortho;
    double Bsum = 0.0, Bmin = 1e30, Bmax = -1e30;
    int    nclip_lo = 0, nclip_hi = 0;
    time_t t0 = time(NULL);

    cellstat = read_pdb(pdbfile, &cell, &atoms, &natoms);
    if (!cellstat) return 0;
    if (natoms < 1) { fprintf(stderr, "ERROR: no atoms in %s\n", pdbfile); return 0; }

    if (have_cell) {
        cell = cell_override;
    } else {
        if (cellstat < 0) {
            fprintf(stderr, "ERROR: no CRYST1 in %s and no cell= given\n", pdbfile);
            return 0;
        }
        cell.a *= super_mult[0];
        cell.b *= super_mult[1];
        cell.c *= super_mult[2];
    }
    if (cell.a <= 0 || cell.b <= 0 || cell.c <= 0) {
        fprintf(stderr, "ERROR: bad cell %g %g %g\n", cell.a, cell.b, cell.c);
        return 0;
    }

    printf("%d atoms from %s\n", natoms, pdbfile);
    printf("cell %.3f %.3f %.3f  %.2f %.2f %.2f\n",
           cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);

    /* Catch the classic trap: supercell coordinates carrying a primitive CRYST1.
       The field would then be built on the wrong frame and every probe after it
       would read the wrong voxel.  Wrapped MD coords can legitimately spill a
       little past the box, so only shout when the span looks like a multiple. */
    {
        double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
        double edge[3];
        int ax;
        for (i = 0; i < natoms; ++i) {
            double p[3] = { atoms[i].x, atoms[i].y, atoms[i].z };
            for (ax = 0; ax < 3; ++ax) {
                if (p[ax] < lo[ax]) lo[ax] = p[ax];
                if (p[ax] > hi[ax]) hi[ax] = p[ax];
            }
        }
        edge[0] = cell.a; edge[1] = cell.b; edge[2] = cell.c;
        for (ax = 0; ax < 3; ++ax) {
            double ratio = (hi[ax] - lo[ax]) / edge[ax];
            if (ratio > 1.9)      /* ~2x is a whole extra cell, not mere spill */
                fprintf(stderr,
                    "WARNING: the coordinates span %.2f x the cell along %c.  If these are\n"
                    "         supercell coordinates with a primitive CRYST1, pass\n"
                    "         super_mult=na,nb,nc or cell=a,b,c so the field is built on\n"
                    "         the frame the atoms actually live in.\n",
                    ratio, "abc"[ax]);
        }
    }

    /* condition the atomic B factors: the map holds final values */
    for (i = 0; i < natoms; ++i) {
        double B = atoms[i].B * Bscale + Boffset;
        if (B < minB)             { B = minB; ++nclip_lo; }
        if (maxB > 0 && B > maxB) { B = maxB; ++nclip_hi; }
        atoms[i].B = B;
        Bsum += B;
        if (B < Bmin) Bmin = B;
        if (B > Bmax) Bmax = B;
    }
    printf("atomic B: %.2f - %.2f  mean %.2f   (Bscale=%g Boffset=%g minB=%g maxB=%g",
           Bmin, Bmax, Bsum / natoms, Bscale, Boffset, minB, maxB);
    if (nclip_lo || nclip_hi) printf("; %d at min, %d at max", nclip_lo, nclip_hi);
    printf(")\n");

    /* Coincident atoms are poison for a spatial field: they contribute N
       different B values to one location, so the field can only return their
       average while the per-atom model keeps them distinct.  Real production
       Bfac.pdb files carry piles of placeholder waters parked at the origin
       (padding to match orignames.pdb), which is exactly this failure. */
    {
        ATOM  *srt = (ATOM *)malloc(natoms * sizeof(ATOM));
        size_t ndup = 0, biggest = 1, run = 1;
        double bx = 0, by = 0, bz = 0;
        if (srt) {
            memcpy(srt, atoms, natoms * sizeof(ATOM));
            qsort(srt, natoms, sizeof(ATOM), atom_pos_cmp);
            for (i = 1; i < natoms; ++i) {
                if (srt[i].x == srt[i-1].x && srt[i].y == srt[i-1].y &&
                    srt[i].z == srt[i-1].z) {
                    ++run; ++ndup;
                    if (run > biggest) {
                        biggest = run;
                        bx = srt[i].x; by = srt[i].y; bz = srt[i].z;
                    }
                } else run = 1;
            }
            free(srt);
            if (biggest > 8 || ndup > (size_t)natoms / 20) {
                fprintf(stderr,
                    "WARNING: %lu of %d atoms share a position with another atom, and\n"
                    "         %lu of them are stacked at (%.3f, %.3f, %.3f).  A B field\n"
                    "         can only return one value per location, so stacked atoms\n"
                    "         get their average.  If these are placeholders, build the\n"
                    "         field from coordinates where every atom is real - e.g. a\n"
                    "         trajectory frame with the B factors merged onto it.\n",
                    (unsigned long)ndup, natoms, (unsigned long)biggest, bx, by, bz);
            }
        }
    }

    /* grid */
    if (nxyz[0] > 0) {
        m.nx = nxyz[0]; m.ny = nxyz[1]; m.nz = nxyz[2];
    } else {
        m.nx = good_grid_size((int)ceil(cell.a / gridspacing));
        m.ny = good_grid_size((int)ceil(cell.b / gridspacing));
        m.nz = good_grid_size((int)ceil(cell.c / gridspacing));
    }
    m.cell = cell;
    m.ispg = ispg;
    nvox = (size_t)m.nx * m.ny * m.nz;
    printf("grid %d x %d x %d  (%.3f x %.3f x %.3f A per voxel)  %.0f MB x2\n",
           m.nx, m.ny, m.nz,
           cell.a / m.nx, cell.b / m.ny, cell.c / m.nz,
           nvox * sizeof(float) / 1048576.0);

    cell_metric(cell, &g11, &g22, &g33, &g12, &g13, &g23);
    cell_recip_diag(cell, &rfx, &rfy, &rfz);
    fr      = cell_frac(cell);
    ortho   = cell_is_ortho(cell);
    twosig2 = 2.0 * sigma * sigma;
    w0      = exp(-fardist * fardist / twosig2);

    /* Truncate each atom's donation where its weight has fallen to wtol*w0, so
       the background already dominates there and the cut is invisible.  The
       resulting step in B is at most wtol*(farB-B) ~ 1 A^2 at the default. */
    rcut = sqrt(fardist * fardist + twosig2 * log(1.0 / wtol));
    printf("sigma %.2f A  fardist %.2f A (w0 = %.3g)  farB %.1f  rcut %.2f A\n",
           sigma, fardist, w0, farB, rcut);
    if (nops > 1) {
        printf("depositing each atom at %d symmetry-equivalent positions\n", nops);
        fprintf(stderr,
            "NOTE: symmetry expansion is for a model that is ONE asymmetric unit.\n"
            "      An MD supercell is already expanded and is deliberately not\n"
            "      symmetric - expanding it again would average symmetry mates\n"
            "      together and erase the differences between subcells.\n");
    }
    if (w0 * farB > 0.5)
        fprintf(stderr, "WARNING: w0*farB = %.2f A^2 - the background is biasing B "
                        "at atom centres.  Raise fardist or lower farB.\n", w0 * farB);

    num = (float *)calloc(nvox, sizeof(float));
    den = (float *)calloc(nvox, sizeof(float));
    if (!num || !den) {
        fprintf(stderr, "ERROR: cannot allocate 2 x %.0f MB for the grid\n",
                nvox * sizeof(float) / 1048576.0);
        free(num); free(den); return 0;
    }

    /* separable-Gaussian scratch: one row of exp() per axis, reused every atom */
    double *ex = NULL, *ey = NULL, *ez = NULL;
    if (ortho) {
        int nxb = (int)ceil(2.0*rcut*rfx*m.nx) + 3;
        int nyb = (int)ceil(2.0*rcut*rfy*m.ny) + 3;
        int nzb = (int)ceil(2.0*rcut*rfz*m.nz) + 3;
        ex = (double *)malloc(nxb * sizeof(double));
        ey = (double *)malloc(nyb * sizeof(double));
        ez = (double *)malloc(nzb * sizeof(double));
        if (!ex || !ey || !ez) {
            fprintf(stderr, "ERROR: out of memory for the deposit box\n");
            free(num); free(den); free(ex); free(ey); free(ez); return 0;
        }
    }

    for (i = 0; i < natoms; ++i) {
        double X = atoms[i].x, Y = atoms[i].y, Z = atoms[i].z;
        double xf0 = fr.f00*X + fr.f01*Y + fr.f02*Z;
        double yf0 =            fr.f11*Y + fr.f12*Z;
        double zf0 =                       fr.f22*Z;
        double Bi = atoms[i].B;
        double wi = useocc ? atoms[i].occ : 1.0;
        int iop;

        if (wi <= 0.0) continue;

      for (iop = 0; iop < nops; ++iop) {
        const SYMOP *op = &ops[iop];
        double xf = op->r[0][0]*xf0 + op->r[0][1]*yf0 + op->r[0][2]*zf0 + op->t[0];
        double yf = op->r[1][0]*xf0 + op->r[1][1]*yf0 + op->r[1][2]*zf0 + op->t[1];
        double zf = op->r[2][0]*xf0 + op->r[2][1]*yf0 + op->r[2][2]*zf0 + op->t[2];
        int ix0, ix1, iy0, iy1, iz0, iz1, kx, ky, kz;

        /* keep the image in the box; the field is periodic anyway */
        xf -= floor(xf); yf -= floor(yf); zf -= floor(zf);

        ix0 = (int)floor((xf - rcut*rfx) * m.nx); ix1 = (int)ceil((xf + rcut*rfx) * m.nx);
        iy0 = (int)floor((yf - rcut*rfy) * m.ny); iy1 = (int)ceil((yf + rcut*rfy) * m.ny);
        iz0 = (int)floor((zf - rcut*rfz) * m.nz); iz1 = (int)ceil((zf + rcut*rfz) * m.nz);

        if (ortho) {
            /* separable: exp(-(dx^2+dy^2+dz^2)/2sig^2) factorizes, so the
               exponentials cost O(box edge) instead of O(box volume) */
            int nxb = ix1 - ix0 + 1, nyb = iy1 - iy0 + 1, nzb = iz1 - iz0 + 1;
            int j;
            for (j = 0; j < nxb; ++j) {
                double d = ((double)(ix0 + j) / m.nx - xf) * cell.a;
                ex[j] = exp(-d * d / twosig2);
            }
            for (j = 0; j < nyb; ++j) {
                double d = ((double)(iy0 + j) / m.ny - yf) * cell.b;
                ey[j] = exp(-d * d / twosig2);
            }
            for (j = 0; j < nzb; ++j) {
                double d = ((double)(iz0 + j) / m.nz - zf) * cell.c;
                ez[j] = exp(-d * d / twosig2);
            }
            for (kz = iz0; kz <= iz1; ++kz) {
                int gz = ((kz % m.nz) + m.nz) % m.nz;
                double wz = wi * ez[kz - iz0];
                for (ky = iy0; ky <= iy1; ++ky) {
                    int gy = ((ky % m.ny) + m.ny) % m.ny;
                    double wzy = wz * ey[ky - iy0];
                    size_t row = (size_t)gz*m.ny*m.nx + (size_t)gy*m.nx;
                    if (wzy < 1e-30) continue;
                    for (kx = ix0; kx <= ix1; ++kx) {
                        int gx = ((kx % m.nx) + m.nx) % m.nx;
                        double w = wzy * ex[kx - ix0];
                        num[row + gx] += (float)(w * Bi);
                        den[row + gx] += (float)w;
                    }
                }
            }
        } else {
            double rcut2 = rcut * rcut;
            for (kz = iz0; kz <= iz1; ++kz) {
                double dfz = (double)kz / m.nz - zf;
                int gz = ((kz % m.nz) + m.nz) % m.nz;
                for (ky = iy0; ky <= iy1; ++ky) {
                    double dfy = (double)ky / m.ny - yf;
                    int gy = ((ky % m.ny) + m.ny) % m.ny;
                    size_t row = (size_t)gz*m.ny*m.nx + (size_t)gy*m.nx;
                    for (kx = ix0; kx <= ix1; ++kx) {
                        double dfx = (double)kx / m.nx - xf;
                        double r2  = g11*dfx*dfx + g22*dfy*dfy + g33*dfz*dfz
                                   + 2.0*(g12*dfx*dfy + g13*dfx*dfz + g23*dfy*dfz);
                        double w;
                        int gx;
                        if (r2 >= rcut2) continue;
                        gx = ((kx % m.nx) + m.nx) % m.nx;
                        w  = wi * exp(-r2 / twosig2);
                        num[row + gx] += (float)(w * Bi);
                        den[row + gx] += (float)w;
                    }
                }
            }
        }
      }
        if (natoms > 20000 && (i % 50000) == 0 && i)
            printf("  %d / %d atoms deposited\n", i, natoms);
    }

    free(ex); free(ey); free(ez);

    /* normalize, with the far-field prior */
    {
        size_t nfar = 0;
        double halfway = 0.5 * (farB + Bmin);
        for (v = 0; v < nvox; ++v) {
            double B = ((double)num[v] + w0 * farB) / ((double)den[v] + w0);
            num[v] = (float)B;
            if (B > halfway) ++nfar;
        }
        printf("%.2f%% of the cell is further than %.2f A from any atom (B -> %.0f)\n",
               100.0 * nfar / (double)nvox, fardist, farB);
    }
    free(den);
    m.rho = num;

    {
        char label[81];
        snprintf(label, sizeof(label),
                 "Bfac_map sigma=%.2f fardist=%.2f farB=%.0f from %.28s",
                 sigma, fardist, farB, pdbfile);
        if (!write_ccp4(outmap, &m, label)) { map_free(&m); free(atoms); return 0; }
    }
    printf("build took %.0f s\n", difftime(time(NULL), t0));

    map_free(&m);
    free(atoms);
    return 1;
}

/* ------------------------------------------------------------------ */
/* probe                                                              */
/* ------------------------------------------------------------------ */
/* Rewrite only columns 61-66 of each ATOM/HETATM record, so whatever the caller
 * had in the rest of the line (element, serial, trailing tags) survives. */
static int do_probe(const char *pdbfile, const char *mapfile, const char *outpdb,
                    CELL cell_override, int have_cell, double super_mult[3],
                    double minB, double maxB, int verbose)
{
    MAP  m;
    CELL cell;
    FRAC fr;
    FILE *in, *out;
    char line[512];
    int  n = 0, cellstat, nover = 0;
    double Bsum = 0.0, Bmin = 1e30, Bmax = -1e30;

    m.rho = NULL;
    if (!read_ccp4(mapfile, &m)) return 0;

    /* The coordinates and the map must live in the same frame.  Trust the map's
       own cell, but shout if the PDB disagrees: a supercell PDB carrying a
       primitive CRYST1 (which is what nc2mtz writes) is the classic trap. */
    {
        CELL pc;
        int  dummy;
        cellstat = read_pdb(pdbfile, &pc, NULL, &dummy);
        if (!cellstat) { map_free(&m); return 0; }
        if (have_cell) {
            cell = cell_override;
        } else if (cellstat > 0) {
            cell = pc;
            cell.a *= super_mult[0];
            cell.b *= super_mult[1];
            cell.c *= super_mult[2];
        } else {
            cell = m.cell;
        }
        if (fabs(cell.a - m.cell.a) > 0.02 * m.cell.a ||
            fabs(cell.b - m.cell.b) > 0.02 * m.cell.b ||
            fabs(cell.c - m.cell.c) > 0.02 * m.cell.c) {
            fprintf(stderr,
                "ERROR: coordinate frame %.2f %.2f %.2f does not match the cell of\n"
                "       %s (%.2f %.2f %.2f).  Sampling would read the wrong voxels.\n"
                "       Pass cell=a,b,c,al,be,ga or super_mult=na,nb,nc.\n",
                cell.a, cell.b, cell.c, mapfile, m.cell.a, m.cell.b, m.cell.c);
            map_free(&m); return 0;
        }
        cell = m.cell;   /* the map is the authority once they agree */
    }
    fr = cell_frac(cell);

    in = fopen(pdbfile, "r");
    if (!in) { fprintf(stderr, "ERROR: cannot open %s\n", pdbfile); map_free(&m); return 0; }
    if (!strcmp(outpdb, "-")) {
        out = stdout;
    } else {
        out = fopen(outpdb, "w");
        if (!out) { fprintf(stderr, "ERROR: cannot write %s\n", outpdb); fclose(in); map_free(&m); return 0; }
    }

    while (fgets(line, sizeof(line), in)) {
        double X, Y, Z, xf, yf, zf, B;
        int len;
        if (strncmp(line, "ATOM  ", 6) && strncmp(line, "HETATM", 6)) {
            fputs(line, out);
            continue;
        }
        len = (int)strlen(line);
        while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        /* pad short records out to the B column so the write below is in range */
        while (len < 66) line[len++] = ' ';
        line[len] = '\0';

        X = pdb_field(line, 30, 8);
        Y = pdb_field(line, 38, 8);
        Z = pdb_field(line, 46, 8);
        xf = fr.f00*X + fr.f01*Y + fr.f02*Z;
        yf =            fr.f11*Y + fr.f12*Z;
        zf =                       fr.f22*Z;
        B  = map_interp(&m, xf, yf, zf);

        if (B < minB)             B = minB;
        if (maxB > 0 && B > maxB) B = maxB;
        if (B > 999.99) { B = 999.99; ++nover; }

        {
            char Bstr[16];
            snprintf(Bstr, sizeof(Bstr), "%6.2f", B);
            memcpy(line + 60, Bstr, 6);
        }
        fprintf(out, "%s\n", line);

        ++n; Bsum += B;
        if (B < Bmin) Bmin = B;
        if (B > Bmax) Bmax = B;
    }
    fclose(in);
    if (out != stdout) fclose(out);

    if (verbose && n)
        fprintf(stderr, "probed %d atoms from %s: B %.2f - %.2f  mean %.2f%s\n",
                n, mapfile, Bmin, Bmax, Bsum / n,
                nover ? "  (some clipped to 999.99 to fit the PDB column)" : "");
    map_free(&m);
    return 1;
}

/* ------------------------------------------------------------------ */
/* stats                                                              */
/* ------------------------------------------------------------------ */
static int do_stats(const char *mapfile)
{
    MAP m;
    size_t nvox, v;
    double sum = 0.0;
    float  lo = 1e30f, hi = -1e30f;
    /* histogram edges chosen to show the sharp/diffuse split at a glance */
    static const double edge[] = { 0, 5, 10, 20, 40, 60, 80, 100, 200, 500, 1e30 };
    int nbin = sizeof(edge) / sizeof(edge[0]) - 1, b;
    size_t *cnt;

    m.rho = NULL;
    if (!read_ccp4(mapfile, &m)) return 0;
    nvox = (size_t)m.nx * m.ny * m.nz;
    cnt  = (size_t *)calloc(nbin, sizeof(size_t));

    for (v = 0; v < nvox; ++v) {
        float x = m.rho[v];
        if (x < lo) lo = x;
        if (x > hi) hi = x;
        sum += x;
        for (b = 0; b < nbin; ++b)
            if (x >= edge[b] && x < edge[b+1]) { ++cnt[b]; break; }
    }

    printf("%s\n", mapfile);
    printf("  grid %d x %d x %d   cell %.3f %.3f %.3f  %.2f %.2f %.2f\n",
           m.nx, m.ny, m.nz, m.cell.a, m.cell.b, m.cell.c,
           m.cell.alpha, m.cell.beta, m.cell.gamma);
    printf("  voxel %.3f x %.3f x %.3f A\n",
           m.cell.a/m.nx, m.cell.b/m.ny, m.cell.c/m.nz);
    printf("  B: %.2f - %.2f  mean %.2f\n", lo, hi, sum / (double)nvox);
    for (b = 0; b < nbin; ++b) {
        if (!cnt[b]) continue;
        if (edge[b+1] > 1e29)
            printf("  B >= %-6.0f      %8.4f%% of cell\n", edge[b],
                   100.0 * cnt[b] / (double)nvox);
        else
            printf("  B %6.0f - %-6.0f %8.4f%% of cell\n", edge[b], edge[b+1],
                   100.0 * cnt[b] / (double)nvox);
    }
    free(cnt);
    map_free(&m);
    return 1;
}

/* ------------------------------------------------------------------ */
static void usage(void)
{
    printf(
"Bfac_map - keep atomic B factors in a CCP4 map instead of on the atoms\n"
"\n"
"  Bfac_map build pdb=Bfac.pdb [outmap=Bfac.map] [options]\n"
"  Bfac_map probe pdb=frame.pdb map=Bfac.map [outpdb=-] [options]\n"
"  Bfac_map stats map=Bfac.map\n"
"\n"
"build - each atom donates its B to the space around it:\n"
"    B(x) = ( SUM_i w_i B_i + w0 farB ) / ( SUM_i w_i + w0 )\n"
"    w_i  = occ_i exp( -r_i^2 / 2 sigma^2 ),   w0 = exp( -fardist^2 / 2 sigma^2 )\n"
"  so B(x) is the local average of nearby atomic B, decaying to farB in space\n"
"  that has no atoms near it.\n"
"\n"
"  sigma=0.5        Gaussian donation width, A.  This is the one parameter that\n"
"                   matters: it sets how far B blends between neighbouring atoms\n"
"                   of different mobility.  Measured against the per-atom path\n"
"                   (1aho P1, dmin 1.0, grid 0.5) R on F was 0.1%% at sigma 0.25,\n"
"                   0.3%% at 0.35, 1.2%% at 0.5, 3.4%% at 1.0, 15%% at 1.5.\n"
"  fardist=4.5      distance at which one atom's vote equals the farB prior, A.\n"
"                   Note this is a threshold, not a taper: with a small sigma the\n"
"                   Gaussian tail is steep, so B goes from local to farB within a\n"
"                   few tenths of an A of fardist.\n"
"  farB=999         B assigned to space far from every atom\n"
"  grid=0.5         voxel size, A (rounded up to a 2,3,5-smooth grid).  Barely\n"
"                   matters next to sigma, but sets the file size (4 bytes per\n"
"                   voxel).  Same test: R on F 1.1%% at 0.25, 1.2%% at 0.5, 1.8%%\n"
"                   at 1.0, 2.4%% at 1.5 - so a big supercell can go coarse.\n"
"  nxyz=nx,ny,nz    exact grid instead of grid=\n"
"  minB=0 maxB=0    clip atomic B into this range before depositing (maxB=0: none)\n"
"  Bscale=1 Boffset=0   applied to atomic B before clipping\n"
"  useocc=0         weight each atom's donation by its occupancy\n"
"  sg=<name|number> deposit each atom at all of its symmetry-equivalent positions\n"
"                   too, so a model that is one asymmetric unit still fills the\n"
"                   cell.  Default P1 (no expansion).  Do NOT use on an MD\n"
"                   supercell: it is already expanded and is deliberately not\n"
"                   symmetric, so expanding it would average symmetry mates.\n"
"  symops=<list>    the operators directly, ';' separated: \"X,Y,Z;-X,-Y,Z\"\n"
"  symoplib=<path>  where to look sg= up  (default CLIBD/symop.lib)\n"
"  wtol=1e-3        donation cutoff radius, as a fraction of the farB prior\n"
"\n"
"probe - trilinear interpolation of the map at each atom, rewriting only the\n"
"  B column (columns 61-66) of each ATOM/HETATM record.  Sampling is periodic in\n"
"  the map's cell, so atoms several cells out, or with unwrapped MD coordinates,\n"
"  read the equivalent voxel.  Space-group symmetry is not applied here - it is\n"
"  baked into the field by build's sg=, so this stays a plain lookup:\n"
"  map=<file>       the B-factor map to read\n"
"  outpdb=-         output PDB ('-' is stdout)\n"
"  minB=0 maxB=0    clip the sampled B (maxB=0: none)\n"
"\n"
"both:\n"
"  cell=a,b,c,al,be,ga   coordinate frame, overriding CRYST1\n"
"  super_mult=na,nb,nc   multiply CRYST1 by this to get the frame.  Needed when\n"
"                        the PDB carries a primitive CRYST1 but supercell coords,\n"
"                        which is what nc2mtz writes.\n"
"\n"
"The map holds final B values: build applies Bscale/Boffset/minB/maxB, probe\n"
"just interpolates.  farB is never clipped by maxB - empty space reads as 999.\n");
}

int main(int argc, char **argv)
{
    const char *mode   = NULL;
    const char *pdb    = NULL;
    const char *mapf   = NULL;
    const char *outmap = "Bfac.map";
    const char *outpdb = "-";
    double sigma = 0.5, gridspacing = 0.5, farB = 999.0, fardist = 4.5;
    double minB = 0.0, maxB = 0.0, Bscale = 1.0, Boffset = 0.0, wtol = 1e-3;
    double super_mult[3] = { 1.0, 1.0, 1.0 };
    int    nxyz[3] = { 0, 0, 0 };
    int    useocc = 0, have_cell = 0, verbose = 1, i;
    CELL   cell_override;
    const char *sgname = NULL, *symops_str = NULL, *symoplib = NULL;
    SYMOP  ops[MAXOPS];
    int    nops = 0, ispg = 1;

    memset(&cell_override, 0, sizeof(cell_override));

    for (i = 1; i < argc; ++i) {
        char *arg = argv[i], *eq = strchr(arg, '=');
        if (!strcmp(arg, "-h") || !strcmp(arg, "--help") || !strcmp(arg, "-help")) {
            usage(); return 0;
        }
        if (!eq) {
            if (!mode && (!strcmp(arg, "build") || !strcmp(arg, "probe") ||
                          !strcmp(arg, "stats"))) { mode = arg; continue; }
            /* bare filenames: .pdb is coordinates, .map is a map */
            {
                size_t L = strlen(arg);
                if (L > 4 && !strcmp(arg + L - 4, ".pdb")) { pdb  = arg; continue; }
                if (L > 4 && !strcmp(arg + L - 4, ".map")) { mapf = arg; continue; }
            }
            fprintf(stderr, "WARNING: ignoring unrecognized argument '%s'\n", arg);
            continue;
        }
        *eq = '\0';
        {
            char *key = arg, *val = eq + 1;
            char *p;
            for (p = key; *p; ++p) *p = (char)tolower((unsigned char)*p);
            if      (!strcmp(key, "pdb") || !strcmp(key, "pdbfile"))  pdb    = val;
            else if (!strcmp(key, "map") || !strcmp(key, "bfacmap") ||
                     !strcmp(key, "mapin") || !strcmp(key, "inmap"))  mapf   = val;
            else if (!strcmp(key, "outmap") || !strcmp(key, "mapout"))outmap = val;
            else if (!strcmp(key, "outpdb") || !strcmp(key, "pdbout"))outpdb = val;
            else if (!strcmp(key, "sigma") || !strcmp(key, "bsigma")) sigma  = atof(val);
            else if (!strcmp(key, "grid")  || !strcmp(key, "gridspacing"))
                                                                 gridspacing = atof(val);
            else if (!strcmp(key, "farb"))                            farB   = atof(val);
            else if (!strcmp(key, "fardist"))                         fardist= atof(val);
            else if (!strcmp(key, "minb"))                            minB   = atof(val);
            else if (!strcmp(key, "maxb"))                            maxB   = atof(val);
            else if (!strcmp(key, "bscale"))                          Bscale = atof(val);
            else if (!strcmp(key, "boffset"))                         Boffset= atof(val);
            else if (!strcmp(key, "useocc"))                          useocc = atoi(val);
            else if (!strcmp(key, "sg") || !strcmp(key, "spacegroup") ||
                     !strcmp(key, "space_group") || !strcmp(key, "smallsg"))
                                                                      sgname = val;
            else if (!strcmp(key, "symops") || !strcmp(key, "symop"))  symops_str = val;
            else if (!strcmp(key, "symoplib") || !strcmp(key, "symop_lib"))
                                                                      symoplib = val;
            else if (!strcmp(key, "wtol"))                            wtol   = atof(val);
            else if (!strcmp(key, "verbose") || !strcmp(key, "debug")) verbose= atoi(val);
            else if (!strcmp(key, "nxyz") || !strcmp(key, "grid_size")) {
                if (sscanf(val, "%d,%d,%d", &nxyz[0], &nxyz[1], &nxyz[2]) != 3 &&
                    sscanf(val, "%d %d %d", &nxyz[0], &nxyz[1], &nxyz[2]) != 3) {
                    fprintf(stderr, "ERROR: nxyz must be nx,ny,nz\n"); return 9;
                }
            }
            else if (!strcmp(key, "super_mult") || !strcmp(key, "mult") ||
                     !strcmp(key, "md_mult")) {
                char buf[128], *q;
                snprintf(buf, sizeof(buf), "%s", val);
                for (q = buf; *q; ++q) if (*q == 'x' || *q == ',') *q = ' ';
                if (sscanf(buf, "%lf %lf %lf",
                           &super_mult[0], &super_mult[1], &super_mult[2]) != 3) {
                    fprintf(stderr, "ERROR: super_mult must be na,nb,nc\n"); return 9;
                }
            }
            else if (!strcmp(key, "cell")) {
                char buf[256], *q;
                snprintf(buf, sizeof(buf), "%s", val);
                for (q = buf; *q; ++q) if (*q == ',') *q = ' ';
                cell_override.alpha = cell_override.beta = cell_override.gamma = 90.0;
                if (sscanf(buf, "%lf %lf %lf %lf %lf %lf",
                           &cell_override.a, &cell_override.b, &cell_override.c,
                           &cell_override.alpha, &cell_override.beta,
                           &cell_override.gamma) < 3) {
                    fprintf(stderr, "ERROR: cell must be a,b,c[,al,be,ga]\n"); return 9;
                }
                have_cell = 1;
            }
            else fprintf(stderr, "WARNING: ignoring unrecognized option '%s'\n", key);
        }
    }

    if (!mode) { usage(); return 9; }

    if (sigma <= 0)   { fprintf(stderr, "ERROR: sigma must be > 0\n");   return 9; }
    if (fardist <= 0) { fprintf(stderr, "ERROR: fardist must be > 0\n"); return 9; }
    if (wtol <= 0 || wtol >= 1) { fprintf(stderr, "ERROR: wtol must be in (0,1)\n"); return 9; }
    if (nxyz[0] > 0 && (nxyz[1] < 1 || nxyz[2] < 1)) {
        fprintf(stderr, "ERROR: nxyz must be nx,ny,nz\n"); return 9;
    }

    if (!strcmp(mode, "build")) {
        if (!pdb) { fprintf(stderr, "ERROR: build needs pdb=<file>\n"); return 9; }
        /* resolve the symmetry operators: explicit list wins over a lookup */
        if (symops_str) {
            nops = read_symops_str(symops_str, ops, MAXOPS);
            if (!nops) return 9;
        } else if (sgname && strcasecmp(sgname, "P1") && strcasecmp(sgname, "P 1")) {
            char path[1024];
            if (symoplib) {
                snprintf(path, sizeof(path), "%s", symoplib);
            } else {
                const char *clibd = getenv("CLIBD");
                snprintf(path, sizeof(path), "%s/symop.lib",
                         clibd ? clibd : "/usr/local/lib/data");
            }
            nops = read_symop_lib(path, sgname, ops, MAXOPS, &ispg);
            if (!nops) return 9;
        }
        if (!nops) {                       /* identity only */
            memset(ops, 0, sizeof(ops[0]));
            ops[0].r[0][0] = ops[0].r[1][1] = ops[0].r[2][2] = 1.0;
            nops = 1;
            ispg = 1;
        }
        return do_build(pdb, outmap, cell_override, have_cell, super_mult,
                        sigma, gridspacing, nxyz, farB, fardist,
                        minB, maxB, Bscale, Boffset, useocc, wtol,
                        ops, nops, ispg) ? 0 : 9;
    }
    if (!strcmp(mode, "probe")) {
        if (!pdb)  { fprintf(stderr, "ERROR: probe needs pdb=<file>\n"); return 9; }
        if (!mapf) { fprintf(stderr, "ERROR: probe needs map=<file>\n"); return 9; }
        return do_probe(pdb, mapf, outpdb, cell_override, have_cell, super_mult,
                        minB, maxB, verbose) ? 0 : 9;
    }
    if (!strcmp(mode, "stats")) {
        if (!mapf) { fprintf(stderr, "ERROR: stats needs map=<file>\n"); return 9; }
        return do_stats(mapf) ? 0 : 9;
    }

    fprintf(stderr, "ERROR: unknown mode '%s'\n", mode);
    usage();
    return 9;
}
