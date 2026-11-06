# summarizes release.csh: mean reported pressure (and its parts) vs time after the restart for restrained (r4),
# released (u4) and restrained-via-ntp=1 (r1) continuations of the same 10 frames; then the Monte Carlo runs.
# usage: python3 release.py   (in pressure_test/)
import glob, re
D = "release/"
KCAL2BAR = 69476.9      # kcal/mol/A^3 -> bar

def records(fn):
    """every printed step before the averages: dict of the energy terms"""
    out, cur = [], None
    for line in open(fn):
        if "A V E R A G E S" in line:
            break
        if line.startswith(" NSTEP"):
            if cur: out.append(cur)
            cur = {}
        if cur is not None:
            for k, v in re.findall(r"([A-Za-z0-9\-\(\) ]+?)\s*=\s*(-?[0-9.]+)", line):
                cur[k.strip()] = float(v)
    if cur: out.append(cur)
    return out

bins = [(0, 0.15), (0.15, 0.55), (0.55, 1.05), (1.05, 2.05), (2.05, 5.05), (5.05, 10.05), (10.05, 20.05)]
runs = {}
for v in ("r4", "u4", "r1"):
    runs[v] = {}
    for fn in sorted(glob.glob(D + "prod_*_%s.out" % v)):
        recs = records(fn)
        if recs:
            runs[v][fn.split("_")[1]] = recs
print("frames:", {v: len(runs[v]) for v in runs})
# are r4 and r1 the same trajectory?  first step where the printed pressure differs
for seg in sorted(runs["r4"], key=int)[:3]:
    a, b = runs["r4"][seg], runs["r1"].get(seg, [])
    diff = next((x["NSTEP"] for x, y in zip(a, b) if abs(x["PRESS"] - y["PRESS"]) > 0.05), None)
    print("prod_%s: r4 and r1 first differ at step %s of %d" % (seg, diff, len(a)))
terms = ["PRESS", "VIRIAL", "EKCMT", "VDWAALS", "EELEC", "1-4 EEL", "RESTRAINT", "TEMP(K)"]
for v in ("r4", "u4", "r1"):
    print("\n%s: mean over frames, by ps after the restart" % v)
    print("%-11s" % "ps" + "".join("%11s" % t for t in terms))
    for lo, hi in bins:
        sel = [r for recs in runs[v].values() for r in recs if lo < r["TIME(PS)"] - recs[0]["TIME(PS)"] + 0.1 <= hi]
        if not sel: continue
        print("%4.1f-%-6.1f" % (lo, hi) + "".join("%11.1f" % (sum(r.get(t, 0) for r in sel) / len(sel)) for t in terms))
# Monte Carlo barostat: mean volume after 0.2 ns
print()
mc = {}
for fn in sorted(glob.glob(D + "mc_*.out")):
    recs = [r for r in records(fn) if "VOLUME" in r]
    if not recs: continue
    t0 = recs[0]["TIME(PS)"]
    late = [r["VOLUME"] for r in recs if r["TIME(PS)"] - t0 >= 200]
    n = len(late)
    if n < 10:
        print(fn, "only", len(recs), "records"); continue
    m = sum(late) / n
    # standard error from 10 blocks
    blk = [sum(late[i * n // 10:(i + 1) * n // 10]) / (n // 10) for i in range(10)]
    se = (sum((b - m) ** 2 for b in blk) / 9 / 10) ** 0.5
    tag = fn.split("/")[-1][3:-4]
    mc[tag] = (m, se)
    print("%-8s <V> %9.1f +/- %4.1f A^3 over %.0f ps (start %.1f)" % (tag, m, se, recs[-1]["TIME(PS)"] - t0 - 200, recs[0]["VOLUME"]))
V0 = 56230.713
for v in ("r", "u"):
    if v + "_1" in mc and v + "_1000" in mc:
        (v1, s1), (v2, s2) = mc[v + "_1"], mc[v + "_1000"]
        kappa = (v1 - v2) / (999 * (v1 + v2) / 2)
        p = 1 + 999 * (v1 - V0) / (v1 - v2)
        print("%s: compressibility %.1f e-6/bar; volume %.1f A^3 at 1 bar is the crystal's volume at P = %.0f bar (+/- %.0f)"
              % ({"r": "restrained", "u": "unrestrained"}[v], kappa * 1e6, V0, p, 999 * s1 / (v1 - v2)))
