"""Complex mean of per-frame FC/PHIC mtzs (gemmi sfcalc --to-mtz or sfcalc_gpu outmtz).
usage: ccp4-python avg_frames.py out.mtz frame1.mtz frame2.mtz ...   (writes FCavg/PHICavg; no scaling)"""
import sys, numpy as np, gemmi
out, fns = sys.argv[1], sys.argv[2:]
acc = None
for fn in fns:
    m = gemmi.read_mtz_file(fn)
    hkl = m.make_miller_array(); F = m.column_with_label('FC').array; P = np.radians(m.column_with_label('PHIC').array)
    z = dict(zip(map(tuple, hkl), F * np.exp(1j * P)))
    if acc is None: acc, keys, cell, sg = z, list(z), m.cell, m.spacegroup
    else:
        for k in keys: acc[k] += z.get(k, np.nan)
keys = [k for k in keys if np.isfinite(acc[k])]
v = np.array([acc[k] for k in keys]) / len(fns)
o = gemmi.Mtz(with_base=True); o.cell = cell; o.spacegroup = sg
o.add_dataset('avg'); o.add_column('FCavg', 'F'); o.add_column('PHICavg', 'P')
o.set_data(np.column_stack([np.array(keys, float), np.abs(v), np.degrees(np.angle(v)) % 360]).astype(np.float32))
o.write_to_file(out); print('%s: mean of %d frames, %d reflections' % (out, len(fns), len(keys)))
