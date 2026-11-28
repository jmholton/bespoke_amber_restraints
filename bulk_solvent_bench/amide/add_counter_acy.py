# Add nadd acetates (copies of ACY 1, C on the B-999 water farthest from protein and salt; waters within 2.6 A of
# them dropped), as setup_s68_builds.sh does.  usage: python3 add_counter_acy.py amberme.pdb nadd   (edits in place)
import sys, numpy as np
f, nadd = sys.argv[1], int(sys.argv[2])
L = open(f).read().split('\n')
cell = np.array([float(x) for x in [l for l in L if l.startswith('CRYST1')][0][6:33].split()])
isat = lambda l: l.startswith(('ATOM', 'HETATM'))
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
atoms = [l for l in L if isat(l)]
prot = np.array([xyz(l) for l in atoms if l[17:20] not in ('HOH', 'WAT', 'ACY', 'NH4') and l[76:78].strip() != 'H'])
salt = [xyz(l) for l in atoms if l[17:20] in ('ACY', 'NH4') and l[76:78].strip() != 'H']
def mind(p, q):
    d = q - p; d -= cell * np.round(d / cell)
    return np.sqrt((d ** 2).sum(1)).min()
cand = [l for l in atoms if l[17:20] in ('HOH', 'WAT') and float(l[60:66]) > 900]
acy1 = [l for l in atoms if l[17:20] == 'ACY' and int(l[22:26]) == 1]
c1 = xyz([l for l in acy1 if l[12:16] == ' C  '][0])
nacy = max(int(l[22:26]) for l in atoms if l[17:20] == 'ACY')
new = []
for k in range(nadd):
    best = max(cand, key=lambda l: min(mind(xyz(l), prot), mind(xyz(l), np.array(salt))))
    p = xyz(best); nacy += 1
    for l in acy1:
        q = xyz(l) - c1 + p
        new.append(l[:22] + '%4d' % nacy + l[26:30] + '%8.3f%8.3f%8.3f' % tuple(q) + l[54:])
        if l[76:78].strip() != 'H': salt.append(q)
    print('  ACY %d at %s: %.1f A from protein' % (nacy, np.round(p, 2), mind(p, prot)))
    cand.remove(best)
heavy = np.array([xyz(s) for s in new if s[76:78].strip() != 'H'])
lastacy = [l for l in L if isat(l) and l[17:20] == 'ACY'][-1]
out, dropped = [], 0
for l in L:
    if isat(l) and l[17:20] in ('HOH', 'WAT') and mind(xyz(l), heavy) < 2.6:
        dropped += 1; continue
    out.append(l)
    if l is lastacy:
        out += new
print('  added %d ACY, dropped %d waters; closest new ACY heavy atom to protein %.2f A, to other salt %.2f A' %
      (nadd, dropped, min(mind(h, prot) for h in heavy),
       min(mind(h, np.array([s for s in salt[:len(salt) - len(heavy)]])) for h in heavy)))
open(f, 'w').write('\n'.join(out))
