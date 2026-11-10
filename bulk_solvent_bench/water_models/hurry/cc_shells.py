# complex CC by resolution shell between solvent maps (FPART/PHIPART), each listed map against each reference
# usage: phenix.python cc_shells.py ref:name=map.mtz ... name=map.mtz ...
import sys
from iotbx import mtz
from cctbx.array_family import flex
def get(fn):
    for a in mtz.object(fn).as_miller_arrays():
        if a.is_complex_array(): return a
refs, maps = [], []
for a in sys.argv[1:]:
    isref = a.startswith("ref:")
    name, fn = (a[4:] if isref else a).split("=", 1)
    (refs if isref else maps).append((name, get(fn)))
common = refs[0][1]
for _, a in refs + maps: common = common.common_set(a)
data = {n: a.common_set(common).data() for n, a in refs + maps}
d = common.d_spacings().data()
shells = [(99, 5), (5, 2.5), (2.5, 1.9), (1.9, 0.96)]
sels = [(d <= hi) & (d > lo) for hi, lo in shells]
def cc(x, y):
    return flex.sum((x * flex.conj(y)).parts()[0]) / (flex.sum(flex.abs(x)**2) * flex.sum(flex.abs(y)**2))**0.5
for rn, _ in refs:
    print("complex CC with %s, shells %s A" % (rn, "  ".join("%g-%g" % s for s in shells)))
    for n, _ in refs + maps:
        if n == rn: continue
        print("  %-12s" % n + "".join("%8.3f" % cc(data[n].select(s), data[rn].select(s)) for s in sels))
