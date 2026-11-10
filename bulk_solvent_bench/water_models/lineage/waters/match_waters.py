# fold the 2019 restrained supercell waters (j5w8760/superwater.pdb, 2x2x3 P1) into the P212121 unit cell and
# match them to the waters of a refinement model; prints per-model-water counts of restrained copies within cut
import sys
import numpy as np
from cctbx import crystal, sgtbx, uctbx
uc = uctbx.unit_cell((45.9, 40.7, 30.1, 90, 90, 90))
ops = [op for op in sgtbx.space_group_info("P 21 21 21").group()]
def waters(fn, only_O=True):
    out = []
    for l in open(fn):
        if l.startswith(("ATOM", "HETATM")) and l[17:20] in ("HOH", "WAT") and (not only_O or l[12:16].strip() == "O"):
            out.append((l[17:27], float(l[30:38]), float(l[38:46]), float(l[46:54]), float(l[54:60]), float(l[60:66]), l[16]))
    return out
sup = waters("/home2/jamesh/projects/amber/1aho_test/j5w8760/superwater.pdb")
S = np.array([uc.fractionalize((x, y, z)) for _, x, y, z, *_ in sup])   # supercell coords in unit-cell fractions
model = waters(sys.argv[1])
cut = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
M = np.array([uc.fractionalize((x, y, z)) for _, x, y, z, *_ in model])
orth = np.array(uc.orthogonalization_matrix()).reshape(3, 3)
# all symmetry copies of model waters
best = np.full(len(S), 1e9); who = np.zeros(len(S), int)
for op in ops:
    R = np.array(op.r().as_double()).reshape(3, 3); t = np.array(op.t().as_double())
    Mo = M @ R.T + t
    d = S[:, None, :] - Mo[None, :, :]
    d -= np.round(d)
    dist = np.sqrt(((d @ orth.T) ** 2).sum(-1))
    j = dist.argmin(1); dm = dist[np.arange(len(S)), j]
    better = dm < best; best[better] = dm[better]; who[better] = j[better]
near = best < cut
print("%d restrained supercell waters (%.1f per protein copy); %d (%.0f%%) within %.1f A of a water in %s" %
      (len(S), len(S) / 48, near.sum(), 100 * near.mean(), cut, sys.argv[1].split("/")[-1]))
cnt = np.bincount(who[near], minlength=len(M))
print("model waters with >=24 of 48 copies restrained nearby: %d of %d;  >=1: %d" % ((cnt >= 24).sum(), len(M), (cnt >= 1).sum()))
np.savetxt(sys.argv[3] if len(sys.argv) > 3 else "counts.txt",
           np.c_[np.arange(len(M)), cnt, [m[4] for m in model], [m[5] for m in model]], fmt="%d %d %.2f %.2f",
           header="model_water_index restrained_copies_within_cut occ B")
