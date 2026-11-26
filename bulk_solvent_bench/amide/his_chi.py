#!/usr/bin/env python3
"""chi1/chi2 of His 54/64 (every chain/copy/altloc) in PDB files.  chi2 = CA-CB-CG-ND1.
usage: his_chi.py file.pdb [...] [res=54,64]  (residue numbers taken mod 64 for supercell copies; resn HIS/HIE/HID/HIP/CHI*)"""
import sys, numpy as np
res = [54, 64]
files = []
for a in sys.argv[1:]:
    if a.startswith('res='): res = [int(x) for x in a[4:].split(',')]
    else: files.append(a)
def dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1; w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))
for f in files:
    R = {}
    for l in open(f):
        if not l.startswith(('ATOM', 'HETATM')): continue
        rn = l[17:20]
        if not (rn.startswith('HI') or rn.startswith('CH') or rn == 'NHE'): continue
        r = int(l[22:26]); rr = (r - 1) % 64 + 1
        if rr not in res: continue
        k = (l[21], r, l[16])
        R.setdefault(k, {})[l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    print('==', f)
    for k in sorted(R):
        a = R[k]
        if not all(x in a for x in ('N', 'CA', 'CB', 'CG', 'ND1')): continue
        c1 = dih(a['N'], a['CA'], a['CB'], a['CG']); c2 = dih(a['CA'], a['CB'], a['CG'], a['ND1'])
        print('  chain %s res %4d alt %s  chi1 %7.1f  chi2 %7.1f  %s' % (k[0], k[1], k[2], c1, c2, 'P' if c2 > 0 else 'M'))
