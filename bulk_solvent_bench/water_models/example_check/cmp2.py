from iotbx import mtz
from cctbx.array_family import flex
def get(fn):
    for a in mtz.object(fn).as_miller_arrays():
        if a.is_complex_array(): return a
names = ["example", "opc_10ns", "opc_r2_2ns", "spce_10ns", "spceold_2ns", "spceold_20ns", "tip3p_10ns", "fb4_10ns"]
m = {n: get(n + "_map.mtz") for n in names}
ref = m["example"]
for n in names: ref = ref.common_set(m[n])
m = {n: m[n].common_set(ref) for n in names}
shells = [(99, 8), (8, 5), (5, 3.5), (3.5, 2.5), (2.5, 1.9), (1.9, 1.4), (1.4, 0.96)]
d = ref.d_spacings().data()
def cc(x, y):
    return flex.sum((x * flex.conj(y)).parts()[0]) / (flex.sum(flex.abs(x)**2) * flex.sum(flex.abs(y)**2))**0.5
print("complex CC by shell (A); first block vs example, second vs opc_10ns")
print("%-13s" % "" + "".join("%9s" % ("%g-%g" % s) for s in shells))
for base in ("example", "opc_10ns"):
    for n in names:
        if n == base: continue
        row = "%-13s" % n
        for hi, lo in shells:
            sel = (d <= hi) & (d > lo)
            row += "%9.3f" % cc(m[base].data().select(sel), m[n].data().select(sel))
        print(row)
    print()
