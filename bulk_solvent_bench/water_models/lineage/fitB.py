# is Y = k*exp(-B s^2/4)*X (+ small) ?  prints CC after best B, k, residual
import sys
from iotbx import mtz
import numpy as np
def get(fn, lab):
    for a in mtz.object(fn).as_miller_arrays():
        if a.info().label_string() == lab: return a
def fit(xf, xl, yf, yl, tag):
    x, y = get(xf, xl), get(yf, yl)
    if not x.is_similar_symmetry(y):
        x = x.customized_copy(crystal_symmetry=y.crystal_symmetry())
    x = x.as_non_anomalous_array().map_to_asu(); y = y.as_non_anomalous_array().map_to_asu()
    if x.anomalous_flag() != y.anomalous_flag(): pass
    x, y = x.common_sets(y)
    s2 = np.array(list(x.d_star_sq().data())); X = np.array(list(x.data())); Y = np.array(list(y.data()))
    best = None
    for B in np.arange(-30, 60.1, 0.25):
        Xb = X * np.exp(-B * s2 / 4)
        k = np.vdot(Xb, Y) / np.vdot(Xb, Xb); r = Y - k * Xb
        rr = np.sqrt(np.mean(abs(r)**2))
        if best is None or rr < best[0]:
            cc = np.real(np.vdot(Xb, Y)) / np.sqrt(np.vdot(Xb, Xb).real * np.vdot(Y, Y).real)
            best = (rr, B, k, cc)
    print("%-55s n %6d  B %6.2f  k %.4f  CC %.5f  rms resid %.3f of %.3f" % (tag, x.size(), best[1], abs(best[2]), best[3], best[0], np.sqrt(np.mean(abs(Y)**2))))
args = sys.argv[1:]
for i in range(0, len(args), 5):
    fit(*args[i:i+5])
