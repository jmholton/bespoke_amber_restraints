# splits the early (0-20 ns) -> late (80-100 ns) change of the w1_100ns solvent map into its water and ion
# (NH4 + acetate) parts: by shell, the share of |late-early|^2 each part carries (plus the cross term), and how
# much each part's own map changed (CC late vs early)
import numpy as np
import gemmi
W = "/home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns/windows/"
def get(fn):
    m = gemmi.read_mtz_file(W + fn)
    return m, m.column_with_label("FPART").array * np.exp(1j * np.radians(m.column_with_label("PHIPART").array))
m, we = get("WAT_HOH/map_1_20.mtz")
_, wl = get("WAT_HOH/map_81_100.mtz")
_, ie = get("NH4_ACY/map_1_20.mtz")
_, il = get("NH4_ACY/map_81_100.mtz")
d = m.make_d_array()
dw, di = wl - we, il - ie
def cc(x, y): return np.real(np.sum(x * np.conj(y))) / np.sqrt(np.sum(abs(x)**2) * np.sum(abs(y)**2))
print("%-10s %8s %8s %8s %10s %10s %10s" % ("shell", "water", "ions", "cross", "CC water", "CC ions", "|ions|/|all|"))
for hi, lo in [(99, 8), (8, 5), (5, 3.5), (3.5, 2.5), (2.5, 1.9), (1.9, 1.4), (1.4, 0.96), (99, 0.96)]:
    s = (d <= hi) & (d > lo)
    tot = np.sum(abs(dw[s] + di[s])**2)
    print("%-10s %7.0f%% %7.0f%% %7.0f%% %10.4f %10.4f %10.3f" % ("%g-%g" % (hi, lo), 100 * np.sum(abs(dw[s])**2) / tot,
          100 * np.sum(abs(di[s])**2) / tot, 100 * 2 * np.real(np.sum(dw[s] * np.conj(di[s]))) / tot,
          cc(wl[s], we[s]), cc(il[s], ie[s]), np.sqrt(np.sum(abs(ie[s])**2) / np.sum(abs(ie[s] + we[s])**2))))
