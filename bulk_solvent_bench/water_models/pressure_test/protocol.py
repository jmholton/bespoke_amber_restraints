# summarizes unrestrained.csh (1 ns, all restraints off, ntp=4/1/2 and MC) and barometer_protocol.csh
# (optimize_weights_runme.com's barometer: 5 ps at dt 0.5 fs, ntp=4 vs ntp=1; OPC box at dt 0.5/1/2 fs).
# usage: python protocol.py   (in pressure_test/)
import re, glob
V0 = 56230.713
def recs(fn, key):
    out = []
    for line in open(fn):
        if "A V E R A G E S" in line: break
        m = re.search(r"\b%s\s*=\s*(-?[0-9.]+)" % re.escape(key), line)
        if m: out.append(float(m.group(1)))
    return out
def mean_se(x, nb=10):
    n = len(x); m = sum(x) / n
    if n < nb: return m, float("nan")
    b = [sum(x[i * n // nb:(i + 1) * n // nb]) / (n // nb) for i in range(nb)]
    return m, (sum((y - m) ** 2 for y in b) / (nb * (nb - 1))) ** 0.5
# unrestrained compressibility of the cell from release/mc_u_1 and mc_u_1000 (after 200 ps)
vu1 = mean_se(recs("release/mc_u_1.out", "VOLUME")[200:])[0]
vu1000 = mean_se(recs("release/mc_u_1000.out", "VOLUME")[200:])[0]
slope = (vu1 - vu1000) / 999
print("1 ns, every restraint off, 1aho cell (709 waters), mean after 0.1 ns:")
for i in (50, 100):
    row = []
    for v in ("u4", "u1", "u2"):
        fn = "unrestrained/prod_%d_%s.out" % (i, v)
        p = recs(fn, "PRESS")[100:] if glob.glob(fn) else []
        row.append("%s %6.0f +/- %3.0f" % (v, *mean_se(p)) if p else "%s   -" % v)
    mfn = "unrestrained/prod_%d_mc.out" % i if i == 50 else "release/mc_u_1.out"
    vm, vse = mean_se(recs(mfn, "VOLUME")[200:])
    row.append("MC true %6.0f +/- %3.0f" % (1 + (vm - V0) / slope, vse / slope))
    print("  prod_%-4d" % i + "   ".join(row))
print()
print("optimize_weights barometer (5 ps, dt 0.5 fs, mean of 10 printed pressures), 1aho cell, all restraints off:")
d = []
for i in range(10, 101, 10):
    a = recs("barometer_protocol/xtal_%d_4.out" % i, "PRESS")
    b = recs("barometer_protocol/xtal_%d_1.out" % i, "PRESS")
    if a and b:
        pa, pb = sum(a) / len(a), sum(b) / len(b)
        d.append((pa, pb))
        print("  prod_%-4d ntp=4 %7.0f   ntp=1 %7.0f   difference %+6.0f" % (i, pa, pb, pa - pb))
if d:
    n = len(d)
    md = sum(a - b for a, b in d) / n
    sd = (sum((a - b - md) ** 2 for a, b in d) / (n - 1)) ** 0.5
    print("  mean over %d frames: ntp=4 %.0f, ntp=1 %.0f, ntp=4 - ntp=1 = %+.0f +/- %.0f" % (
        n, sum(a for a, b in d) / n, sum(b for a, b in d) / n, md, sd / n ** 0.5))
print()
print("OPC box, true +103 bar at this volume, 50 ps, mean after 5 ps:")
for dt in ("0.5", "1", "2"):
    row = []
    for p in ("4", "1"):
        fn = "barometer_protocol/water_%s_%s.out" % (dt, p)
        x = recs(fn, "PRESS") if glob.glob(fn) else []
        x = x[len(x) // 10:]
        row.append("ntp=%s %6.0f +/- %3.0f" % (p, *mean_se(x)) if x else "ntp=%s -" % p)
    print("  dt %-3s fs   " % dt + "   ".join(row))
