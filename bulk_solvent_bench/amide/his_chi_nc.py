#!/usr/bin/env python3
"""His ring chi2 (CA-CB-CG-ND1) over saved amber_N.nc of a run, every His in xtal.prmtop (protein atoms come first,
so prmtop indices hold for the smaller saved trajectories).
usage: his_chi_nc.py rundir first_itr last_itr"""
import sys, os, numpy as np
from scipy.io import netcdf_file
d, i0, i1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
txt = open(os.path.join(d, 'xtal.prmtop')).read().split('%FLAG ')
sec = {s.split('\n', 1)[0].strip(): s.split('\n', 2)[2] if s.count('\n') > 1 else '' for s in txt[1:]}
def fixed(s, w): s = s.replace('\n', ''); return [s[i:i + w].strip() for i in range(0, len(s), w)]
names = fixed(sec['ATOM_NAME'], 4); labs = fixed(sec['RESIDUE_LABEL'], 4)
ptr = [int(x) - 1 for x in sec['RESIDUE_POINTER'].split()] + [len(names)]
his = []
for r, lab in enumerate(labs):
    if lab.startswith(('HI', 'HD', 'HE', 'HP')) or lab in ('CHIE', 'CHID', 'CHIP', 'NHIE', 'NHID'):
        a = {names[i]: i for i in range(ptr[r], ptr[r + 1])}
        if all(k in a for k in ('CA', 'CB', 'CG', 'ND1')):
            his.append((r + 1, lab, [a[k] for k in ('CA', 'CB', 'CG', 'ND1')], [a.get(k) for k in ('N', 'CA', 'CB', 'CG')]))
def dih(p):
    b0, b1, b2 = p[:, 0] - p[:, 1], p[:, 2] - p[:, 1], p[:, 3] - p[:, 2]
    b1 = b1 / np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - (b0 * b1).sum(1)[:, None] * b1; w = b2 - (b2 * b1).sum(1)[:, None] * b1
    return np.degrees(np.arctan2((np.cross(b1, v) * w).sum(1), (v * w).sum(1)))
c2 = {h[0]: [] for h in his}; c1 = {h[0]: [] for h in his}; nf = 0
for i in range(i0, i1 + 1):
    f = os.path.join(d, 'amber_%d.nc' % i)
    if not os.path.exists(f): continue
    X = netcdf_file(f, 'r', mmap=False).variables['coordinates'][:].astype(float)
    nf += len(X)
    for h in his:
        c2[h[0]] += list(dih(X[:, h[2]])); c1[h[0]] += list(dih(X[:, h[3]]))
print('%s itr %d-%d, %d frames' % (os.path.basename(d.rstrip('/')), i0, i1, nf))
for r, lab, _, _ in his:
    a2 = np.array(c2[r]); a1 = np.array(c1[r])
    # chi1 is circular: report rotamer fractions (g+ 0..120, t |chi1|>120, g- -120..0) and the circular mean
    cm1 = np.degrees(np.arctan2(np.sin(np.radians(a1)).mean(), np.cos(np.radians(a1)).mean()))
    gp, tt = ((a1 >= 0) & (a1 < 120)).mean(), (np.abs(a1) >= 120).mean()
    print('  res %4d (%d in copy %d) %-4s  chi1 g+ %.2f t %.2f g- %.2f (circ mean %6.1f)  chi2: frac>0 %.2f  '
          'mean(<0) %6.1f  mean(>0) %6.1f' %
          (r, (r - 1) % 64 + 1, (r - 1) // 64 + 1, lab, gp, tt, 1 - gp - tt, cm1, (a2 > 0).mean(),
           a2[a2 < 0].mean() if (a2 < 0).any() else np.nan, a2[a2 > 0].mean() if (a2 > 0).any() else np.nan))
