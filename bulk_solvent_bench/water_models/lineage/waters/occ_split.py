# refined water occupancy summed over the 2019 'core' sites (>=24 of 48 supercell copies restrained within 1 A)
# vs the other sites.  usage: occ_split.py counts.txt refined_001.pdb [more refined pdbs, averaged]
import sys
import numpy as np
cnt = np.loadtxt(sys.argv[1], usecols=1).astype(int)
core = cnt >= 24
tot = []
for fn in sys.argv[2:]:
    occ = [float(l[54:60]) for l in open(fn) if l.startswith(("ATOM", "HETATM")) and l[17:20] in ("HOH", "WAT") and l[12:16].strip() == "O"]
    occ = np.array(occ)
    if len(occ) != len(cnt):
        print("water count mismatch", fn, len(occ), len(cnt)); sys.exit(1)
    tot.append((occ[core].sum(), occ[~core].sum()))
tot = np.array(tot).mean(0)
print("core %d sites occ sum %5.1f (mean %.2f) | other %d sites occ sum %5.1f (mean %.2f)" % (core.sum(), tot[0], tot[0] / core.sum(), (~core).sum(), tot[1], tot[1] / (~core).sum()))
