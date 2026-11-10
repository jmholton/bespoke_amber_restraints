from iotbx import mtz
def fl(fn, lab="FreeR_flag"):
    return {x.info().label_string(): x for x in mtz.object(fn).as_miller_arrays()}[lab]
E = "/home/jamesh/projects/phenix_user_solvent/example/"
ref = fl(E + "refme.mtz")
for fn in (E + "1aho.mtz", "/home2/jamesh/projects/amber/bulk_solvent/test_1aho_long/1aho.mtz"):
    for lab in ("FreeR_flag", "FREE"):
        a = fl(fn, lab); r, b = ref.common_sets(a)
        free = (r.data() == 0)
        print(fn.split("/")[-3], lab, "n", b.size(), "values", sorted(set(b.data()))[:3], "...",
              "refme-free that are", lab, "==0: %d/%d" % ((b.data().select(free) == 0).count(True), free.count(True)),
              " ==1: %d" % ((b.data().select(free) == 1).count(True)))
