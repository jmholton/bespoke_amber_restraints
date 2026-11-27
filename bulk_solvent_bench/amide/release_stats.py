"""Which restraints optimize_weights releases (2026-11-26).  release_worst_N.log: seeds = restraints at the weight cap
("finding challenged restraints with weight > ..." block), then symmetry mates and neighbour shells are reset to the
release weight ("moved" lines).  Per run over an iteration range: iterations with a release, seeds per iteration, the
most frequent seeds (protein atoms pooled over copies, residue mod 64; waters labelled by nearest protein atom in the
crystal model), and how many released atoms are water.
usage: python3 release_stats.py run1 [run2 ...] [itr=101-199] [top=8]"""
import sys, os, numpy as np
args = [a for a in sys.argv[1:] if '=' not in a]
opt = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
i0, i1 = [int(x) for x in opt.get('itr', '101-199').split('-')]
top = int(opt.get('top', 8))
model = '/scratch/jamesh/projects/amber/1aho_hurry_test/build_s68_PD/refined.pdb'
cell = None; prot = []; plab = []
for l in open(model):
    if l.startswith('CRYST1'): cell = np.array([float(x) for x in l[6:33].split()])
    if l.startswith(('ATOM', 'HETATM')) and l[17:20] != 'HOH' and l[76:78].strip() != 'H' and int(l[22:26]) <= 256:
        prot.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        plab.append('%s %s%d' % (l[12:16].strip(), l[17:20].strip().title(), (int(l[22:26]) - 1) % 64 + 1))
prot = np.array(prot)
def lab(l):
    rn = l[17:20].strip()
    if rn in ('HOH', 'WAT'):
        x = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        d = prot - x; d -= cell * np.round(d / cell); r = np.sqrt((d ** 2).sum(1))
        return 'HOH@%s' % plab[r.argmin()]
    if rn in ('ACY', 'NH4'): return '%s %s' % (l[12:16].strip(), rn)
    r = rn.title()
    if r in ('Hip', 'Hid', 'Hie'): r = 'His'
    return '%s %s%d' % (l[12:16].strip(), r, (int(l[22:26]) - 1) % 64 + 1)
for rd in args:
    name = os.path.basename(os.path.abspath(rd)).replace('ff19_opc_s68_', '')
    nit = nrel = 0; seeds = {}; nseed = 0; nmoved = 0; nwat = 0
    for i in range(i0, i1 + 1):
        f = os.path.join(rd, 'release_worst_%d.log' % i)
        if not os.path.exists(f): continue
        nit += 1
        L = open(f).read().split('\n'); inseed = False; s = 0
        for l in L:
            if l.startswith('finding challenged'): inseed = True; continue
            if inseed and l.startswith(('ATOM', 'HETATM')):
                k = lab(l); seeds[k] = seeds.get(k, 0) + 1; s += 1; continue
            if inseed and not l.startswith(('ATOM', 'HETATM')): inseed = False
            if ' moved ' in l:
                nmoved += 1
                if ' HOH ' in l: nwat += 1
        nseed += s; nrel += (s > 0)
    if not nit: continue
    print('%-6s itr %d-%d: release in %d of %d itr, %.1f seeds and %.0f atoms reset per itr (%.0f%% water); top seeds:' %
          (name, i0, i1, nrel, nit, nseed / nit, nmoved / nit, 100 * nwat / max(nmoved, 1)))
    print('        ' + '  '.join('%s %d' % (k, v) for k, v in sorted(seeds.items(), key=lambda x: -x[1])[:top]))
