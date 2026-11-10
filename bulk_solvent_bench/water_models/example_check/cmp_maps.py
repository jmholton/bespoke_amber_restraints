# compare bulk-solvent map coefficients with the phenix_user_solvent example Fpart, by resolution shell
import sys
from iotbx import mtz
from cctbx.array_family import flex
E = "/home/jamesh/projects/phenix_user_solvent/example/"
P = "/home2/jamesh/projects/amber/bulk_solvent/water_models/phenix/"
def arrays(fn):
    return {a.info().label_string(): a for a in mtz.object(fn).as_miller_arrays()}
ref = arrays(E + "refme.mtz"); small = arrays("/home2/jamesh/projects/amber/bulk_solvent/phenix_check/refme_small.mtz")
for lab in ("FP,SIGFP", "FreeR_flag"):
    a, b = ref[lab], small[lab]
    a, b = a.common_sets(b)
    print("%-11s refme vs refme_small: %d common, max |diff| %.4g" % (lab, a.size(), flex.max(flex.abs(a.data() - b.data()))))
fr = ref["FreeR_flag"]
print("FreeR_flag values:", sorted(set(fr.data())) [:5], "... free(0) fraction %.3f" % ((fr.data() == 0).count(True) / fr.size()))
ex = arrays(E + "solvent_Fpart.mtz")["Fpart,PHIpart"]
maps = [("example Fpart", ex)] + [(t, arrays(P + f)["FPART,PHIPART"]) for t, f in
        [("opc 10 ns", "opc_10ns.mtz"), ("opc r2 2 ns", "opc_r2_2ns_map.mtz"), ("spce old 2 ns", "old_2ns_map.mtz"),
         ("spce old 20 ns", "old_20ns_map.mtz"), ("tip3p 10 ns", "tip3p_10ns.mtz")]]
nb = 8
ex.setup_binner(n_bins=nb)
print("\nvs example Fpart: complex correlation / amplitude ratio (map / example) per shell")
hdr = "%-15s" % "d (A)" + "".join("%11s" % ("%.2f-%.2f" % ex.binner().bin_d_range(i)) for i in ex.binner().range_used())
print(hdr + "   overall")
for t, m in maps:
    a, b = ex.common_sets(m)
    a.setup_binner(n_bins=nb)
    row = "%-15s" % t
    for i in a.binner().range_used():
        s = a.binner().selection(i)
        x, y = a.data().select(s), b.data().select(s)
        cc = flex.sum((x * flex.conj(y)).parts()[0]) / (flex.sum(flex.abs(x)**2) * flex.sum(flex.abs(y)**2))**0.5
        r = (flex.sum(flex.abs(y)**2) / flex.sum(flex.abs(x)**2))**0.5
        row += "%6.2f/%-4.2f" % (cc, r)
    x, y = a.data(), b.data()
    cc = flex.sum((x * flex.conj(y)).parts()[0]) / (flex.sum(flex.abs(x)**2) * flex.sum(flex.abs(y)**2))**0.5
    row += "  %5.2f/%.2f" % (cc, (flex.sum(flex.abs(y)**2) / flex.sum(flex.abs(x)**2))**0.5)
    print(row)
