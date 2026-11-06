# example Fpart = a*Fpart(refme250) + b*exp(-B s^2/4)*Fpart2(refme250), grid over B
from iotbx import mtz
import numpy as np
S = "/home2/jamesh/projects/amber/1aho_refine/best_Rfree/longrun/squish25/"
E = "/home/jamesh/projects/phenix_user_solvent/example/"
def load(fn): return {a.info().label_string(): a for a in mtz.object(fn).as_miller_arrays()}
ex = load(E + "solvent_Fpart.mtz")["Fpart,PHIpart"]
import sys
d = load(S + (sys.argv[1] if len(sys.argv) > 1 else "refme250.mtz"))
a, b = d["Fpart,PHIpart"].common_sets(d["Fpart2,PHIpart2"])
e, a = ex.common_sets(a); e, b = ex.common_sets(b)
s2 = np.array(list(a.d_star_sq().data()))
y = np.array(list(e.data())); A0 = np.array(list(a.data())); B0 = np.array(list(b.data()))
best = None
for B in np.arange(0, 40.1, 0.5):
    A = np.vstack([A0, B0 * np.exp(-B * s2 / 4)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None); r = y - A @ c
    rr = np.sqrt(np.mean(abs(r)**2))
    if best is None or rr < best[0]: best = (rr, B, c, abs(r).max())
print("best B %.1f  a %.4f  b %.4f  rms residual %.4f (of %.2f)  max %.3f" % (best[1], best[2][0].real, best[2][1].real, best[0], np.sqrt(np.mean(abs(y)**2)), best[3]))
