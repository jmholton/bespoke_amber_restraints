# how the scan solvent maps (and the example, 2019 and 10-ns maps) correlate: complex CC by resolution shell
# and the amplitude ratio, then a matrix of the low-resolution CC (inf-5 A, where bulk solvent matters most)
from iotbx import mtz
from cctbx.array_family import flex
X = "/home2/jamesh/projects/amber/bulk_solvent/water_models/example_check/"
def get(fn):
    for a in mtz.object(X + fn).as_miller_arrays():
        if a.is_complex_array(): return a
tags = ["w01", "w03", "w1", "w1_r2", "w3", "w10", "flat10", "pm500", "pp500", "pp1500", "ac0", "ac0_r2", "opc_r2"]
m = {t: get("scan_%s_map.mtz" % t) for t in tags}
m["opc_10ns"] = get("opc_10ns_map.mtz")
m["example"] = get("example_FPART.mtz")
m["2019"] = get("md2019_map.mtz")
names = list(m)
ref = m["w1"]
for a in m.values(): ref = ref.common_set(a)
m = {k: a.common_set(ref).data() for k, a in m.items()}
d = ref.d_spacings().data()
shells = [(99, 8), (8, 5), (5, 3.5), (3.5, 2.5), (2.5, 1.9), (1.9, 1.4), (1.4, 0.96)]
sels = [(d <= hi) & (d > lo) for hi, lo in shells] + [d > 0]
hdr = "%-10s" % "" + "".join("%9s" % ("%g-%g" % s) for s in shells) + "%9s" % "all"
def cc(x, y):
    return flex.sum((x * flex.conj(y)).parts()[0]) / (flex.sum(flex.abs(x)**2) * flex.sum(flex.abs(y)**2))**0.5
def ratio(x, y):
    return (flex.sum(flex.abs(x)**2) / flex.sum(flex.abs(y)**2))**0.5
for base in ("w1", "example", "2019"):
    print("complex CC with %s by shell (A)" % base)
    print(hdr)
    for n in names:
        if n == base: continue
        print("%-10s" % n + "".join("%9.3f" % cc(m[n].select(s), m[base].select(s)) for s in sels))
    print()
print("rms amplitude ratio to w1 by shell")
print(hdr)
for n in names:
    if n == "w1": continue
    print("%-10s" % n + "".join("%9.3f" % ratio(m[n].select(s), m["w1"].select(s)) for s in sels))
print()
low = d > 5
print("CC matrix, inf-5 A")
print("%-8s" % "" + "".join("%7s" % n[:6] for n in names))
for a in names:
    print("%-8s" % a[:8] + "".join("%7.3f" % cc(m[a].select(low), m[b].select(low)) for b in names))
print()
high = d <= 3.5
print("CC matrix, 3.5-0.96 A")
print("%-8s" % "" + "".join("%7s" % n[:6] for n in names))
for a in names:
    print("%-8s" % a[:8] + "".join("%7.3f" % cc(m[a].select(high), m[b].select(high)) for b in names))
