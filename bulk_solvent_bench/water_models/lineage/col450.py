from iotbx import mtz
from cctbx.array_family import flex
S = "/home2/jamesh/projects/amber/1aho_refine/best_Rfree/longrun/squish25/"
E = "/home/jamesh/projects/phenix_user_solvent/example/"
def load(fn): return {a.info().label_string(): a for a in mtz.object(fn).as_miller_arrays()}
ex = load(E + "solvent_Fpart.mtz")["Fpart,PHIpart"]
exr = load(E + "refme.mtz")["Fpart,PHIpart"]
def cmp(t, a):
    x, y = ex.common_sets(a)
    cc = flex.sum((x.data() * flex.conj(y.data())).parts()[0]) / (flex.sum(flex.abs(x.data())**2) * flex.sum(flex.abs(y.data())**2))**0.5
    print("%-40s CC %.5f ratio %.4f maxdiff %.3g" % (t, cc, (flex.sum(flex.abs(y.data())**2) / flex.sum(flex.abs(x.data())**2))**0.5, flex.max(flex.abs(x.data() - y.data()))))
cmp("example refme.mtz Fpart", exr)
for f in ("refme_scaled_400.mtz", "refme_scaled_450.mtz", "refme_scaled_500.mtz", "refme_500.mtz", "refme1.mtz"):
    d = load(S + f)
    for lab in ("Fpart,PHIpart", "Fpart2,PHIpart2"):
        if lab in d: cmp(f + " " + lab.split(",")[0], d[lab])
    if "Fpart2,PHIpart2" in d:
        a, b = d["Fpart,PHIpart"], d["Fpart2,PHIpart2"]
        a, b = a.common_sets(b)
        cmp(f + " 1.0078*Fpart+0.7365*Fpart2", a.customized_copy(data=a.data() * 1.0078 + b.data() * 0.7365))
