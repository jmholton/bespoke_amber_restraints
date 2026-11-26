"""Local Fo-Fc at the His ring atoms and ring-H sites of an MD run, averaged over iterations (2026-11-25; controls for
the 0.68 M charge tests, after the refinement session's H omit maps and His54 ring-flip test).

For each iteration's ensemble render avg_N.mtz (FCavg/PHICavg, P212121): Fobs (refme.mtz FP, P1) matched at a
P212121-equivalent index, k and overall B fitted (least squares, Fo ~ k Fc exp(-B s^2/4)), and the difference
coefficients D = Fo/k - Fc exp(-B s^2/4) on model phases, in the model's electron units.  D is averaged over the
iterations in reciprocal space (= averaging the maps), in blocks for an iteration-to-iteration spread, and Fourier
transformed.  Sites are taken from the crystal model (refined.pdb, copy 1; His54 there is the M ring orientation, the
one the flip test favours), so every arm is read at the same places whatever its labels:
  ring atoms CG ND1 CD2 CE1 NE2, and the ideal ring-H sites at X-ray distances on the external bisector (N-H 0.86,
  C-H 0.93 A), for His54 and His64;
  every backbone amide H site (0.86 A, bisector of C(i-1)-N-CA): modelled in every arm, so its mean residual is the
  arm's zero point for a modelled H.
Values are the map at the site (trilinear) and its mean within 0.5 A, in e/A^3 and in units of the map's rms.
  swap index = mean(ND1, NE2 sites) - mean(CD2, CE1 sites): > 0 when the data want more electrons at the crystal N
  positions than the model has there (a ring with N and C exchanged reads about +1 e of difference per pair).
usage: ccp4-python his_fofc.py rundir first_itr last_itr [block=10] [model=build_s68_PD/refined.pdb] [fobs=refme.mtz] [free=0|1]
  free=1 uses only the free-set reflections for the map (default: all)."""
import sys, os, numpy as np, gemmi

args = [a for a in sys.argv[1:] if '=' not in a]
opt = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
rundir, i0, i1 = args[0], int(args[1]), int(args[2])
block = int(opt.get('block', 10))
model = opt.get('model', '/scratch/jamesh/projects/amber/1aho_hurry_test/build_s68_PD/refined.pdb')
fobs_fn = opt.get('fobs', os.path.join(rundir, 'refme.mtz'))

# --- Fobs, P1 -> lookup by any 222 sign combination
m = gemmi.read_mtz_file(fobs_fn)
hkl = m.make_miller_array().astype(int)
FP = m.column_with_label('FP').array.astype(float)
free = m.column_with_label('FreeR_flag').array.astype(float) if 'FreeR_flag' in [c.label for c in m.columns] else None
obs = {}
for j, (h, k, l) in enumerate(hkl):
    if np.isfinite(FP[j]):
        obs[(h, k, l)] = (FP[j], free[j] if free is not None else 1)

ref_mtz = None; base_hkl = None; Fo = None; sel_ok = None; s2 = None


def load_avg(fn):
    global ref_mtz, base_hkl, Fo, sel_ok, s2
    a = gemmi.read_mtz_file(fn)
    h = a.make_miller_array().astype(int)
    F = a.column_with_label('FCavg').array.astype(float)
    P = np.radians(a.column_with_label('PHICavg').array.astype(float))
    if ref_mtz is None:
        ref_mtz = a; base_hkl = h
        fo = np.full(len(h), np.nan); fl = np.full(len(h), -1.0)
        for j, (x, y, z) in enumerate(h):
            for e in ((x, y, z), (-x, y, z), (x, -y, z), (x, y, -z), (-x, -y, z), (-x, y, -z), (x, -y, -z), (-x, -y, -z)):
                if e in obs:
                    fo[j], fl[j] = obs[e]; break
        Fo = fo
        sel_ok = np.isfinite(fo)
        if opt.get('free', '0') == '1':
            sel_ok &= (fl == 0)
        s2 = 1 / a.make_d_array() ** 2
    elif len(h) != len(base_hkl) or (h != base_hkl).any():
        raise SystemExit('%s: reflection list differs from the first file' % fn)
    return F, P


def diff_coeffs(F, P):
    ok = sel_ok & np.isfinite(F) & np.isfinite(P)
    best = None
    for B in np.arange(-20, 20.01, 0.25):
        c = F[ok] * np.exp(-B * s2[ok] / 4)
        k = (Fo[ok] * c).sum() / (c * c).sum()
        r = ((Fo[ok] - k * c) ** 2).sum()
        if best is None or r < best[0]:
            best = (r, k, B)
    _, k, B = best
    D = np.zeros(len(F), complex)
    D[ok] = (Fo[ok] / k - F[ok] * np.exp(-B * s2[ok] / 4)) * np.exp(1j * P[ok])
    R = np.abs(Fo[ok] / k - F[ok] * np.exp(-B * s2[ok] / 4)).sum() / (Fo[ok] / k).sum()
    return D, k, B, R


def to_map(D):
    mt = gemmi.Mtz(with_base=True)
    mt.spacegroup = ref_mtz.spacegroup; mt.set_cell_for_all(ref_mtz.cell)
    mt.add_dataset('d'); mt.add_column('D', 'F'); mt.add_column('PHID', 'P')
    data = np.column_stack([base_hkl.astype(np.float32), np.abs(D).astype(np.float32),
                            np.degrees(np.angle(D)).astype(np.float32)])
    mt.set_data(data)
    g = mt.transform_f_phi_to_map('D', 'PHID', sample_rate=4)
    return g


# --- sites from the crystal model
st = gemmi.read_structure(model)
res = {}
for ch in st[0]:
    for r in ch:
        if 1 <= r.seqid.num <= 64 and gemmi.find_tabulated_residue(r.name) is not None and \
                gemmi.find_tabulated_residue(r.name).is_amino_acid():
            res.setdefault(r.seqid.num, r)
def P_(r, n):
    a = r.find_atom(n, '*')
    return None if a is None else np.array(a.pos.tolist())
def ext(atom, a, b, dist):
    u = (atom - a) / np.linalg.norm(atom - a) + (atom - b) / np.linalg.norm(atom - b)
    return atom + dist * u / np.linalg.norm(u)
sites = []   # (label, xyz)
for rn in (54, 64):
    r = res[rn]
    X = {n: P_(r, n) for n in ('CG', 'ND1', 'CD2', 'CE1', 'NE2')}
    for n in ('CG', 'ND1', 'CD2', 'CE1', 'NE2'):
        sites.append(('His%d %s' % (rn, n), X[n]))
    sites.append(('His%d HD1 site' % rn, ext(X['ND1'], X['CG'], X['CE1'], 0.86)))
    sites.append(('His%d HE2 site' % rn, ext(X['NE2'], X['CD2'], X['CE1'], 0.86)))
    sites.append(('His%d HD2 site' % rn, ext(X['CD2'], X['CG'], X['NE2'], 0.93)))
    sites.append(('His%d HE1 site' % rn, ext(X['CE1'], X['ND1'], X['NE2'], 0.93)))
amide = []
for rn in sorted(res):
    if rn < 2 or rn - 1 not in res or res[rn].name == 'PRO': continue
    N, CA, C0 = P_(res[rn], 'N'), P_(res[rn], 'CA'), P_(res[rn - 1], 'C')
    if N is None or CA is None or C0 is None: continue
    amide.append(ext(N, C0, CA, 0.86))

# small sphere for the 0.5 A mean
g1 = np.mgrid[-0.5:0.51:0.125, -0.5:0.51:0.125, -0.5:0.51:0.125].reshape(3, -1).T
g1 = g1[np.linalg.norm(g1, axis=1) <= 0.5]


def region(rn, rad=1.0, step=0.2):
    # sample points within rad of His rn's ring atoms and ring-H sites (the local fit region)
    cen = np.array([x for lab, x in sites if lab.startswith('His%d ' % rn) and lab.split()[1] != 'CG' or lab == 'His%d CG' % rn])
    lo, hi = cen.min(0) - rad, cen.max(0) + rad
    pts = np.mgrid[lo[0]:hi[0]:step, lo[1]:hi[1]:step, lo[2]:hi[2]:step].reshape(3, -1).T
    dmin = np.min(np.linalg.norm(pts[:, None, :] - cen[None, :, :], axis=2), axis=1)
    return pts[dmin <= rad]
reg = {rn: region(rn) for rn in (54, 64)}


def local_rms(g, rn):
    v = np.array([g.interpolate_value(gemmi.Position(*x)) for x in reg[rn]])
    return np.sqrt((v ** 2).mean())


def sample(g, xyz):
    pt = g.interpolate_value(gemmi.Position(*xyz))
    sph = np.mean([g.interpolate_value(gemmi.Position(*(xyz + d))) for d in g1])
    return pt, sph


# --- run through the iterations
itrs = [i for i in range(i0, i1 + 1) if os.path.exists(os.path.join(rundir, 'avg_%d.mtz' % i))]
if not itrs:
    raise SystemExit('no avg_N.mtz in %s for itr %d-%d' % (rundir, i0, i1))
blocks = {}; tot = None; ks = []; Bs = []; Rs = []
for i in itrs:
    F, P = load_avg(os.path.join(rundir, 'avg_%d.mtz' % i))
    D, k, B, R = diff_coeffs(F, P)
    ks.append(k); Bs.append(B); Rs.append(R)
    tot = D if tot is None else tot + D
    b = (i - itrs[0]) // block
    blocks[b] = blocks.get(b, 0) + D
    blocks.setdefault(('n', b), 0); blocks[('n', b)] += 1
tot /= len(itrs)
g = to_map(tot)
arr = np.array(g, copy=False); rms = arr.std()
bvals = []
nb = [b for b in blocks if not isinstance(b, tuple)]
for b in nb:
    gb = to_map(blocks[b] / blocks[('n', b)])
    bvals.append([sample(gb, x)[1] for _, x in sites] + [np.mean([sample(gb, x)[1] for x in amide])] +
                 [local_rms(gb, 54), local_rms(gb, 64), np.array(gb, copy=False).std()])
bvals = np.array(bvals)

name = os.path.basename(os.path.abspath(rundir))
print('%s: itr %d-%d (%d renders, %d blocks of %d), Fo-Fc map rms %.4f e/A^3; fitted k %.4g, B %+.2f, R %.2f%%'
      % (name, itrs[0], itrs[-1], len(itrs), len(nb), block, rms, np.mean(ks), np.mean(Bs), 100 * np.mean(Rs)))
am = [sample(g, x) for x in amide]
zp = np.mean([s[1] for s in am])
print('  backbone amide H sites (n %d, modelled): mean(point) %+.4f  mean(0.5 A) %+.4f e/A^3  = zero point' %
      (len(am), np.mean([s[0] for s in am]), zp))
print('  %-16s %9s %9s %9s %8s %9s' % ('site', 'point', '0.5A', '0.5A-zp', '/rms', 'block_sd'))
vals = {}
for j, (lab, x) in enumerate(sites):
    pt, sp = sample(g, x)
    vals[lab] = sp
    sd = bvals[:, j].std(ddof=1) if len(nb) > 1 else np.nan
    print('  %-16s %+9.4f %+9.4f %+9.4f %+8.2f %9.4f' % (lab, pt, sp, sp - zp, sp / rms, sd))
for rn in (54, 64):
    s = (vals['His%d ND1' % rn] + vals['His%d NE2' % rn] - vals['His%d CD2' % rn] - vals['His%d CE1' % rn]) / 2
    j = [i for i, t in enumerate(sites) if t[0] in ('His%d ND1' % rn, 'His%d NE2' % rn, 'His%d CD2' % rn, 'His%d CE1' % rn)]
    sb = (bvals[:, j[0]] + bvals[:, j[3]] - bvals[:, j[1]] - bvals[:, j[2]]) / 2
    print('  His%d swap index (N sites - C sites) %+.4f e/A^3 (%+.2f rms; block sd %.4f)' %
          (rn, s, s / rms, sb.std(ddof=1) if len(nb) > 1 else np.nan))
ns = len(sites) + 1
for c, rn in enumerate((54, 64)):
    lr = local_rms(g, rn)
    rel = bvals[:, ns + c] / bvals[:, ns + 2] if len(nb) > 1 else np.array([np.nan])
    print('  His%d local Fo-Fc rms within 1.0 A of ring atoms + ring-H sites (%d points): %.4f e/A^3 = %.3f x map rms '
          '(block sd %.3f)' % (rn, len(reg[rn]), lr, lr / rms, rel.std(ddof=1) if len(nb) > 1 else np.nan))
