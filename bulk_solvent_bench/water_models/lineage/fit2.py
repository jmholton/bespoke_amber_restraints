# least-squares: example Fpart = a*Fpart(start) + b*Fpart2(increment), complex a, b
from iotbx import mtz
from cctbx.array_family import flex
import numpy as np
S = "/home2/jamesh/projects/amber/1aho_refine/best_Rfree/longrun/squish25/"
E = "/home/jamesh/projects/phenix_user_solvent/example/"
def load(fn): return {a.info().label_string(): a for a in mtz.object(fn).as_miller_arrays()}
ex = load(E + "solvent_Fpart.mtz")["Fpart,PHIpart"]
d = load(S + "refme1.mtz")
a, b = d["Fpart,PHIpart"], d["Fpart2,PHIpart2"]
a, b = a.common_sets(b); e, a = ex.common_sets(a); e, b = ex.common_sets(b)
A = np.array([list(a.data()), list(b.data())]).T; y = np.array(list(e.data()))
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
r = y - A @ coef
print("coef a = %.4f%+.4fi  b = %.4f%+.4fi" % (coef[0].real, coef[0].imag, coef[1].real, coef[1].imag))
print("rms example %.3f  rms residual %.4f  max residual %.3f  (n=%d)" % (np.sqrt(np.mean(abs(y)**2)), np.sqrt(np.mean(abs(r)**2)), abs(r).max(), len(y)))
