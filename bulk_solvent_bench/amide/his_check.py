# His54/His64 residue labels, atoms and charges per copy in a build's xtal.prmtop, plus the total charge.
# usage: his_check.py xtal.prmtop
import sys
lines = open(sys.argv[1]).read().split('%FLAG ')
blk = {l.split('\n')[0].strip(): l.split('\n')[2:] for l in lines}
names = [s[i:i+4].strip() for s in blk['ATOM_NAME'] for i in range(0, len(s), 4) if s[i:i+4].strip()]
labels = [s[i:i+4].strip() for s in blk['RESIDUE_LABEL'] for i in range(0, len(s), 4) if s[i:i+4].strip()]
ptr = [int(x) for s in blk['RESIDUE_POINTER'] for x in s.split()] + [len(names) + 1]
chg = [float(x) / 18.2223 for s in blk['CHARGE'] for x in s.split()]
for c in range(4):
    out = []
    for r in (54 + 64 * c, 64 + 64 * c):
        a0, a1 = ptr[r - 1] - 1, ptr[r] - 1
        ring = [n for n in names[a0:a1] if n in ('HD1', 'HE2', 'HN1', 'HN2')]
        out.append('%s%d %s q %+.3f' % (labels[r - 1], r, '/'.join(ring), sum(chg[a0:a1])))
    print('copy %d: %s' % (c + 1, '   '.join(out)))
print('total charge %+.3f, %d atoms' % (sum(chg), len(names)))
