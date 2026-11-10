# biggest mFo-DFc features of a phenix refinement: flood-fill blobs above 3 rms, ranked by integrated density,
# with the nearest model atom; plus difference-map rms inside / outside the protein mask
# usage: diffblobs.py tag model_001.pdb model_001.mtz [ntop]
import sys
import numpy as np
import gemmi
tag, pdb, mtz = sys.argv[1:4]
ntop = int(sys.argv[4]) if len(sys.argv) > 4 else 8
st = gemmi.read_structure(pdb); st.setup_entities()
m = gemmi.read_mtz_file(mtz)
grid = m.transform_f_phi_to_map("FOFCWT", "PHFOFCWT", sample_rate=3)
arr = np.array(grid, copy=False)
rms = float(np.sqrt((arr ** 2).mean()))
# protein mask (waters removed): 1 = solvent region, 0 = protein
dry = st.clone(); dry.remove_waters()
mask = gemmi.FloatGrid(grid.nu, grid.nv, grid.nw); mask.set_unit_cell(grid.unit_cell); mask.spacegroup = grid.spacegroup
masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Cctbx)
masker.put_mask_on_float_grid(mask, dry[0])
mk = np.array(mask, copy=False) > 0.5
# 2019 core water sites (deposited waters restrained in >=24 of 48 copies)
core = []
dep = gemmi.read_structure("/home2/jamesh/projects/amber/bulk_solvent/phenix_check/withwater.pdb")
cnt = np.loadtxt("/home2/jamesh/projects/amber/bulk_solvent/water_models/lineage/waters/counts_deposited.txt", usecols=1)
iw = 0
for res in dep[0][0]:
    if res.name == "HOH":
        if cnt[iw] >= 24: core.append(res[0].pos)
        iw += 1
ns = gemmi.NeighborSearch(st[0], st.cell, 6).populate()
cell = st.cell
def nearest(pos):
    best = (99, "none")
    for mk_ in ns.find_atoms(pos, "\0", radius=6):
        cra = mk_.to_cra(st[0])
        if cra.atom.is_hydrogen(): continue
        d = cell.find_nearest_image(pos, cra.atom.pos, gemmi.Asu.Any).dist()
        if d < best[0]:
            best = (d, "%s%s %s%d %s" % (cra.chain.name, "", cra.residue.name, cra.residue.seqid.num, cra.atom.name + (cra.atom.altloc if cra.atom.altloc != "\0" else "")))
    return best
def coredist(pos):
    return min([cell.find_nearest_image(pos, c, gemmi.Asu.Any).dist() for c in core]) if core else 99
def blobs(g, sign):
    b = gemmi.find_blobs_by_flood_fill(g, cutoff=3 * rms, min_volume=0.1, min_score=0, min_peak=0)
    b = sorted(b, key=lambda x: -x.score)
    return b
print("== %s   map rms %.4f   rms in solvent region %.4f, in protein %.4f  (solvent region %.0f%% of cell)" %
      (tag, rms, np.sqrt((arr[mk] ** 2).mean()), np.sqrt((arr[~mk] ** 2).mean()), 100 * mk.mean()))
for sign in (1, -1):
    g = gemmi.FloatGrid(np.ascontiguousarray(sign * arr, dtype=np.float32), grid.unit_cell, grid.spacegroup)
    bl = blobs(g, sign)
    tot = sum(x.score for x in bl)
    atcore = [x for x in bl if coredist(x.peak_pos) < 1.0]
    print("  %s blobs > 3 rms: %d, total score %.1f; on 2019 core sites (<1 A): %d, score %.1f; top %d:" %
          ("positive" if sign > 0 else "negative", len(bl), tot, len(atcore), sum(x.score for x in atcore), ntop))
    for x in bl[:ntop]:
        d, lab = nearest(x.peak_pos)
        print("   %+6.1f rms  score %6.2f  vol %5.1f A^3  nearest %-22s %4.1f A  2019-core-site %4.1f A" %
              (sign * x.peak_value / rms, x.score, x.volume, lab, d, coredist(x.peak_pos)))
