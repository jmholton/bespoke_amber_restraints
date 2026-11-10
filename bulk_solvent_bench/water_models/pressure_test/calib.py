# summarizes calib.csh.  usage: python release.py-style, in pressure_test/: python calib.py
import glob, re
D = "calib/"
kT = 1.380649e-23 * 287      # J
def series(fn, key):
    out = []
    for line in open(fn):
        if "A V E R A G E S" in line:
            break
        m = re.search(r"\b%s\s*=\s*(-?[0-9.]+)" % re.escape(key), line)
        if m:
            out.append(float(m.group(1)))
    return out
# Monte Carlo reference: <V> and compressibility from the volume fluctuations, after 50 ps
v = series(D + "water_mc.out", "VOLUME")[100:]
n = len(v)
vm = sum(v) / n
var = sum((x - vm) ** 2 for x in v) / n
kappa = var * 1e-30 / (kT * vm) * 1e5          # 1/bar
blk = [sum(v[i * n // 10:(i + 1) * n // 10]) / (n // 10) for i in range(10)]
se = (sum((b - vm) ** 2 for b in blk) / 90) ** 0.5
print("MC barostat at 1 bar: <V> = %.1f +/- %.1f A^3 (%d samples), kappa from fluctuations %.1f e-6/bar"
      % (vm, se, n, kappa * 1e6))
for p in ("p4", "p1", "p2"):
    fn = D + "water_%s.out" % p
    vol = series(fn, "VOLUME")
    pr = series(fn, "PRESS")[100:]
    if not pr:
        print(p, "no output"); continue
    V = vol[0]
    m = sum(pr) / len(pr)
    k = len(pr)
    blk = [sum(pr[i * k // 10:(i + 1) * k // 10]) / (k // 10) for i in range(10)]
    s = (sum((b - m) ** 2 for b in blk) / 90) ** 0.5
    truth = 1 + (vm - V) / (kappa * V)
    print("ntp=%s at V = %.1f: reads %7.1f +/- %4.1f bar; the MC equation of state says %7.1f +/- %.0f; error %+.0f bar"
          % (p[1], V, m, s, truth, se / (kappa * V), m - truth))
print()
print("frozen coordinates (dt 1e-5 ps), reported pressure at steps 1, 2, 10, 20:")
for fn in sorted(glob.glob(D + "frozen_*.out")):
    pr = series(fn, "PRESS")
    vir = series(fn, "VIRIAL")
    if len(pr) < 20:
        print(fn, "only", len(pr)); continue
    print("%-22s P %8.1f %8.1f %8.1f %8.1f   VIRIAL %8.1f %8.1f" % (fn.split("/")[-1][:-4], pr[0], pr[1], pr[9], pr[19], vir[0], vir[19]))
