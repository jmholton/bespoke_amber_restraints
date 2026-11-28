"""NH4+ at the HOH51 site (His64 ND1 / Tyr42 O) in the _n runs (2026-11-28): per copy, nearest NH4 N and nearest
water O to His64 ND1 and Tyr42 O in each end-of-iteration restart amber_N.rst7 (full system, ascii).
usage: nh4_site.py rundir first_itr last_itr [step]"""
import sys, os, numpy as np
d, i0, i1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); st = int(sys.argv[4]) if len(sys.argv) > 4 else 1
txt = open(os.path.join(d, 'xtal.prmtop')).read().split('%FLAG ')
def sec(name, w):
    s = [t for t in txt if t.startswith(name + ' ')][0].split('\n', 2)[2].replace('\n', '')
    return [s[i:i + w].strip() for i in range(0, len(s), w)]
names = sec('ATOM_NAME', 4); rlab = sec('RESIDUE_LABEL', 4); rptr = [int(x) - 1 for x in sec('RESIDUE_POINTER', 8)]
nat = len(names); rptr.append(nat)
res = np.zeros(nat, int)
for r in range(len(rlab)): res[rptr[r]:rptr[r + 1]] = r
def atom(r, n): return [a for a in range(rptr[r], rptr[r + 1]) if names[a] == n][0]
nd1 = [atom(63 + 64 * k, 'ND1') for k in range(4)]; o42 = [atom(41 + 64 * k, 'O') for k in range(4)]
nh4 = np.array([a for a in range(nat) if rlab[res[a]] == 'NH4' and names[a] == 'N'])
wat = np.array([a for a in range(nat) if rlab[res[a]] in ('WAT', 'HOH') and names[a] == 'O'])
def read(f):
    L = open(f).read().split('\n'); n = int(L[1].split()[0])
    v = np.array(' '.join(L[2:]).split(), float)
    return v[:3 * n].reshape(n, 3), v[-6:-3]
def mind(x, Y, box):
    dd = Y - x; dd -= box * np.round(dd / box); r = np.sqrt((dd ** 2).sum(1)); return r.min(), r.argmin()
print('itr  per copy: NH4-ND1 / NH4-O42 (which NH4 res) | nearest water O-ND1 (only waters of the final topology: approximate before the last water add)')
for i in range(i0, i1 + 1, st):
    f = os.path.join(d, 'amber_%d.rst7' % i)
    if not os.path.exists(f): continue
    x, box = read(f); row = []
    if len(x) < nh4.max() + 1: print("%4d  skip: %d atoms vs %d in prmtop" % (i, len(x), nat)); continue
    # earlier iterations had fewer waters (added later); protein and salt come first, so their indices still hold
    w = wat[wat < len(x)] if len(x) != nat else wat
    for k in range(4):
        dn, j = mind(x[nd1[k]], x[nh4], box); do = np.sqrt(((lambda v: v - box * np.round(v / box))(x[nh4[j]] - x[o42[k]]) ** 2).sum())
        dw, _ = mind(x[nd1[k]], x[w], box)
        row.append('%4.2f/%4.2f(%d)|%4.2f' % (dn, do, res[nh4[j]] + 1, dw))
    print('%4d  ' % i + '   '.join(row))
