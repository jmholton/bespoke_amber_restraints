# how the w1_100ns solvent map converges: complex CC of window maps (w1_100ns/windows/map_<first>_<last>.mtz,
# production ns first..last) by resolution shell.  Uses whatever maps exist so far.
import glob, os, re
from iotbx import mtz
from cctbx.array_family import flex
W = "/home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns/windows/"
def get(fn):
    for a in mtz.object(fn).as_miller_arrays():
        if a.is_complex_array(): return a
m = {}
for fn in glob.glob(W + "map_*_*.mtz"):
    f, l = map(int, re.match(r"map_(\d+)_(\d+)[.]mtz", os.path.basename(fn)).groups())
    m[(f, l)] = get(fn)
ref = list(m.values())[0]
for a in m.values(): ref = ref.common_set(a)
m = {k: a.common_set(ref).data() for k, a in m.items()}
d = ref.d_spacings().data()
shells = [(99, 8), (8, 5), (5, 3.5), (3.5, 2.5), (2.5, 1.9), (1.9, 1.4), (1.4, 0.96)]
sels = [(d <= hi) & (d > lo) for hi, lo in shells] + [d > 0]
hdr = "".join("%9s" % ("%g-%g" % s) for s in shells) + "%9s" % "all"
def cc(x, y):
    return flex.sum((x * flex.conj(y)).parts()[0]) / (flex.sum(flex.abs(x)**2) * flex.sum(flex.abs(y)**2))**0.5
def row(a, b):
    return "".join("%9.4f" % cc(m[a].select(s), m[b].select(s)) for s in sels)
def name(k): return "%d-%d ns" % (k[0] - 1, k[1])
# the longest cumulative map so far is the reference
full = max((k for k in m if k[0] == 1), key=lambda k: k[1])
last = full[1]
cum = sorted(k for k in m if k[0] == 1 and k[1] < last)
if cum:
    print("cumulative maps vs the 0-%d ns map" % last)
    print("%-12s" % "" + hdr)
    for k in cum: print("%-12s" % name(k) + row(k, full))
    print()
blocks = sorted(k for k in m if k[1] - k[0] == 9 and k[0] % 10 == 1)
if len(blocks) > 1:
    print("10-ns blocks: CC with the other blocks, all shells (drift shows as early blocks matching late ones worse)")
    print("%-12s" % "" + "".join("%9s" % name(b) for b in blocks))
    for a in blocks:
        print("%-12s" % name(a) + "".join("%9s" % ("" if a == b else "%.4f" % cc(m[a], m[b])) for b in blocks))
    print()
    print("mean CC between different 10-ns blocks, by shell; 'vs truth' = sqrt(2r/(1+r)), the CC of one 10-ns map")
    print("with the infinitely long average if the blocks are independent")
    rs = []
    for s in sels:
        v = [cc(m[a].select(s), m[b].select(s)) for i, a in enumerate(blocks) for b in blocks[i + 1:]]
        rs.append(sum(v) / len(v))
    print("%-12s" % "" + hdr)
    print("%-12s" % "r(10 ns)" + "".join("%9.4f" % r for r in rs))
    print("%-12s" % "vs truth" + "".join("%9.4f" % (2 * r / (1 + r))**0.5 for r in rs))
    print()
    short = sorted(k for k in m if k[1] - k[0] == 1)
    if short:
        print("2-ns windows: mean CC with the 10-ns blocks that do not contain them (a drift test for short runs)")
        print("%-12s" % "" + hdr)
        for k in short:
            others = [b for b in blocks if b[1] < k[0] or b[0] > k[1]]
            vals = ["%9.4f" % (sum(cc(m[k].select(s), m[b].select(s)) for b in others) / len(others)) for s in sels]
            print("%-12s" % name(k) + "".join(vals))
