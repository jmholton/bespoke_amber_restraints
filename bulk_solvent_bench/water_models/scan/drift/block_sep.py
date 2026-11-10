# drift vs sampling noise in the w1_100ns 10-ns block maps: mean CC by shell for block pairs grouped by how far
# apart they are (1 = adjacent).  Noise alone gives the same CC at every separation; drift makes it fall.
import numpy as np
import gemmi
W = "/home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns/windows/"
f = {}
for b in range(10):
    m = gemmi.read_mtz_file(W + "map_%d_%d.mtz" % (10 * b + 1, 10 * b + 10))
    f[b] = m.column_with_label("FPART").array * np.exp(1j * np.radians(m.column_with_label("PHIPART").array))
d = m.make_d_array()
shells = [(99, 8), (8, 5), (5, 3.5), (3.5, 2.5), (2.5, 1.9), (1.9, 1.4), (1.4, 0.96)]
def cc(x, y): return np.real(np.sum(x * np.conj(y))) / np.sqrt(np.sum(abs(x)**2) * np.sum(abs(y)**2))
print("%-12s" % "separation" + "".join("%9s" % ("%g-%g" % s) for s in shells))
for lo, hi, lab in [(1, 1, "10 ns"), (2, 2, "20 ns"), (3, 4, "30-40 ns"), (5, 6, "50-60 ns"), (7, 9, "70-90 ns")]:
    pairs = [(a, b) for a in range(10) for b in range(a + 1, 10) if lo <= b - a <= hi]
    row = [np.mean([cc(f[a][(d <= h) & (d > l)], f[b][(d <= h) & (d > l)]) for a, b in pairs]) for h, l in shells]
    print("%-12s" % lab + "".join("%9.4f" % r for r in row))
