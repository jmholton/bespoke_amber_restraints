# true pressure of each scan run's restrained cell from the Monte Carlo barostat (mc_scan.csh) vs what production's
# ntp=4 reported (pressure_vs_seg.txt, prod_1 and prod_2).  usage: python mc_scan.py   (in pressure_test/)
import re
S = "/home2/jamesh/projects/amber/bulk_solvent/water_models/scan/"
V0 = 56230.713
def vols(fn):
    v, t = [], []
    for line in open(fn):
        if "A V E R A G E S" in line: break
        m = re.search(r"TIME\(PS\) =\s*([0-9.]+)", line)
        if m: t.append(float(m.group(1)))
        m = re.search(r"VOLUME\s*=\s*([0-9.]+)", line)
        if m: v.append(float(m.group(1)))
    return [x for x, tt in zip(v, t) if tt - t[0] >= 200]
def mean_se(x):
    n = len(x); m = sum(x) / n
    b = [sum(x[i * n // 10:(i + 1) * n // 10]) / (n // 10) for i in range(10)]
    return m, (sum((y - m) ** 2 for y in b) / 90) ** 0.5
# w1 restrained compressibility: <V> at 1 and 1000 bar
v1, _ = mean_se(vols("release/mc_r_1.out"))
v1000, _ = mean_se(vols("release/mc_r_1000.out"))
slope = (v1 - v1000) / 999        # A^3 per bar
print("restrained cell: %.3f A^3/bar (kappa %.1f e-6/bar)" % (slope, slope / v1 * 1e6))
print("%-8s %6s %14s %16s %18s %8s" % ("run", "waters", "ntp=4 reading", "MC <V> at 1 bar", "true P (bar)", "offset"))
for t in ("pm500", "w1", "pp500", "pp1500", "pp3000"):
    v, se = mean_se(vols("mc_scan/%s.out" % t))
    p = 1 + (v - V0) / slope
    rd = [float(l.split()[1]) for l in open(S + t + "/pressure_vs_seg.txt") if l.startswith("prod_")]
    nw = [l.split()[-1] for l in open(S + t + "/pressure_vs_seg.txt") if l.startswith("hydrate_")][-1]
    print("%-8s %6s %14.0f %10.1f+/-%3.0f %12.0f +/- %3.0f %8.0f" % (t, nw, sum(rd) / len(rd), v, se, p, se / slope, sum(rd) / len(rd) - p))
