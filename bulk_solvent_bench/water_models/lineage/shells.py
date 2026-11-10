from iotbx import mtz
from cctbx.array_family import flex
def get(fn, lab):
    for a in mtz.object(fn).as_miller_arrays():
        if a.info().label_string() == lab: return a.as_non_anomalous_array().map_to_asu()
E = "/home/jamesh/projects/phenix_user_solvent/example/solvent_Fpart.mtz"
P = "/home2/jamesh/projects/amber/bulk_solvent/water_models/phenix/"
md19 = get("solvent_RestAll_sf.mtz", "FC,PHIC")
ex = get(E, "Fpart,PHIpart")
ours = {"opc 10 ns": get(P + "opc_10ns.mtz", "FPART,PHIPART"), "spce 10 ns": get(P + "spce_10ns.mtz", "FPART,PHIPART")}
md19 = md19.customized_copy(crystal_symmetry=ex.crystal_symmetry())
shells = [(99, 8), (8, 5), (5, 3.5), (3.5, 2.5), (2.5, 1.9), (1.9, 1.4), (1.4, 0.96)]
def row(t, x, y):
    x, y = x.common_sets(y); d = x.d_spacings().data(); s = "%-28s" % t
    for hi, lo in shells:
        sel = (d <= hi) & (d > lo); a, b = x.data().select(sel), y.data().select(sel)
        s += "%8.3f" % (flex.sum((a * flex.conj(b)).parts()[0]) / (flex.sum(flex.abs(a)**2) * flex.sum(flex.abs(b)**2))**0.5)
    print(s)
print("%-28s" % "complex CC, shell (A)" + "".join("%8s" % ("%g-%g" % s) for s in shells))
row("2019 MD vs example", md19, ex)
for t, m in ours.items():
    row(t + " vs 2019 MD", m, md19)
for t, m in ours.items():
    row(t + " vs example", m, ex)
