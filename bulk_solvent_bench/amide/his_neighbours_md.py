"""What sits next to the His54 / His64 ring N in the MD (2026-11-26): per run, over saved amber_N.nc, the nearest
non-ring heavy atoms to each ring N and the ring-edge atoms of copy 1-4 (minimum image), as fractions of frames
within 3.2 A, by partner (residue type + atom, residues mod 64).  Also the restrained waters whose reference point
lies within 3.2 A of a His ring N in the crystal model (restraints_for_N.pdb).
usage: python3 his_neighbours_md.py rundir first last"""
import sys, os, numpy as np
from scipy.io import netcdf_file
d, i0, i1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
txt = open(os.path.join(d, 'xtal.prmtop')).read().split('%FLAG ')
sec = {s.split('\n', 1)[0].strip(): s.split('\n', 2)[2] if s.count('\n') > 1 else '' for s in txt[1:]}
def fixed(s, w): s = s.replace('\n', ''); return [s[i:i + w].strip() for i in range(0, len(s), w)]
names = fixed(sec['ATOM_NAME'], 4); labs = fixed(sec['RESIDUE_LABEL'], 4)
ptr = [int(x) - 1 for x in sec['RESIDUE_POINTER'].split()] + [len(names)]
resof = np.zeros(len(names), int)
for r in range(len(labs)): resof[ptr[r]:ptr[r + 1]] = r
heavy = np.array([not n.startswith('H') and n != 'EPW' for n in names])
def lab(i):
    r = resof[i]; L = labs[r]
    return '%s %s%s' % (names[i], L, '' if L in ('WAT', 'HOH', 'ACY', 'NH4') else str(r % 64 + 1))
targets = []
for c in range(4):
    for rn, ats in ((54, ('ND1', 'NE2', 'CE1')), (64, ('ND1', 'NE2', 'CE1'))):
        r = rn - 1 + 64 * c
        for a in ats:
            i = [j for j in range(ptr[r], ptr[r + 1]) if names[j] == a][0]
            targets.append(('His%d %s' % (rn, a), i, r))
cnt = {}; nf = 0
for it in range(i0, i1 + 1):
    f = os.path.join(d, 'amber_%d.nc' % it)
    if not os.path.exists(f): continue
    nc = netcdf_file(f, 'r', mmap=False)
    X = nc.variables['coordinates'][:].astype(float); L = nc.variables['cell_lengths'][:].astype(float)
    n = X.shape[1]
    hv = np.where(heavy[:n])[0]
    for fr in range(len(X)):
        nf += 1
        for t, i, r in targets:
            dv = X[fr, hv] - X[fr, i]; dv -= L[fr] * np.round(dv / L[fr])
            dd = np.sqrt((dv ** 2).sum(1))
            ok = (dd < 3.2) & (resof[hv] != r)
            for j in hv[ok]:
                k = (t, lab(j).split()[0] + ' ' + (lab(j).split()[1] if labs[resof[j]] not in ('WAT', 'HOH') else 'HOH'))
                cnt[k] = cnt.get(k, 0) + 1
print('%s itr %d-%d, %d frames (x 4 copies): partners within 3.2 A, fraction of copy-frames' % (os.path.basename(d), i0, i1, nf))
for t in sorted(set(k[0] for k in cnt)):
    row = sorted([(v, k[1]) for k, v in cnt.items() if k[0] == t], reverse=True)[:6]
    print('  %-10s ' % t + '  '.join('%s %.2f' % (p, v / nf / 4) for v, p in row))
