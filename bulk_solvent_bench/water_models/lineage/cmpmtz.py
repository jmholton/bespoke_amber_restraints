# compare candidate mtz files with the phenix_user_solvent example: Fpart (complex CC, ratio, max diff), FP, FreeR_flag
import sys
from iotbx import mtz
from cctbx.array_family import flex
E = "/home/jamesh/projects/phenix_user_solvent/example/"
def load(fn):
    out = {}
    for a in mtz.object(fn).as_miller_arrays():
        out[a.info().label_string()] = a
    return out
ex = load(E + "solvent_Fpart.mtz")["Fpart,PHIpart"]
er = load(E + "refme.mtz")
efp, efr = er["FP,SIGFP"], er["FreeR_flag"]
for fn in sys.argv[1:]:
    try:
        arr = load(fn)
    except Exception as e:
        print(fn, "unreadable", e); continue
    out = [fn.replace("/home2/jamesh/projects/amber/", "~amber/").replace("/home/jamesh/projects/", "~p/")]
    for lab, a in arr.items():
        if not a.is_similar_symmetry(ex):
            out.append("%s: different symmetry %s" % (lab, a.unit_cell().parameters())); break
        if a.is_complex_array():
            x, y = ex.common_sets(a)
            if x.size() == 0: continue
            cc = flex.sum((x.data() * flex.conj(y.data())).parts()[0]) / (flex.sum(flex.abs(x.data())**2) * flex.sum(flex.abs(y.data())**2))**0.5
            r = (flex.sum(flex.abs(y.data())**2) / flex.sum(flex.abs(x.data())**2))**0.5
            md = flex.max(flex.abs(x.data() - y.data()))
            out.append("%s: CC %.4f ratio %.4f maxdiff %.3g n %d" % (lab, cc, r, md, x.size()))
        elif lab.startswith("FP,") or lab == "FP":
            x, y = efp.common_sets(a)
            if x.size(): out.append("%s vs ex FP maxdiff %.3g n %d" % (lab, flex.max(flex.abs(x.data() - y.data())), x.size()))
        elif "free" in lab.lower():
            x, y = efr.common_sets(a)
            if x.size():
                fx = (x.data() == 0)
                same = (y.data().select(fx) == 0).count(True)
                out.append("%s: ex-free==0 here %d/%d" % (lab, same, fx.count(True)))
    print("\n   ".join(out))
