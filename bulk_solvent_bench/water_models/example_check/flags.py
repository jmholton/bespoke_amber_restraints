from iotbx import mtz
from collections import Counter
a = {x.info().label_string(): x for x in mtz.object("/home/jamesh/projects/phenix_user_solvent/example/refme.mtz").as_miller_arrays()}["FreeR_flag"]
b = {x.info().label_string(): x for x in mtz.object("/home2/jamesh/projects/amber/bulk_solvent/phenix_check/refme_small.mtz").as_miller_arrays()}["FreeR_flag"]
a, b = a.common_sets(b)
print("refme values:", Counter(a.data()).most_common(3), "n values", len(set(a.data())))
print("refme_small values:", Counter(b.data()).most_common(3), "n values", len(set(b.data())))
c = Counter(zip(a.data(), b.data()))
print("refme==0 -> refme_small:", sorted((k[1], v) for k, v in c.items() if k[0] == 0)[:6])
print("refme_small==0 -> refme:", sorted((k[0], v) for k, v in c.items() if k[1] == 0)[:6])
