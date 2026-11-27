"""Worst restraints of MD runs, per atom, averaged over iterations and the 4 copies (2026-11-26; user: "I'm more
interested in how protonation impacted the worst restraints. Did they move?").

Per iteration optimize_weights writes restraints_challenge_vs_weight_N.txt: for every restrained atom that has moved
off its reference point, challenge = w d^2 (restraint energy), d (A) and w (weight), sorted.  Atoms not listed have
not moved (counted as 0).  restraints_for_N.pdb holds every restraint's reference point and weight (B x pdbscale).
Protein atoms are pooled over the 4 copies (residue number mod 64).  Waters are labelled by the nearest protein heavy
atom to their reference point in the crystal model (minimum image in the P1 cell that holds the 4 copies), so they
can be matched across runs whose water numbering differs.
usage: python3 worst_restraints.py run1 [run2 ...] [itr=101-199] [top=20] [focus=1] [model=.../refined.pdb]"""
import sys, os, glob, numpy as np

args = [a for a in sys.argv[1:] if '=' not in a]
opt = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
i0, i1 = [int(x) for x in opt.get('itr', '101-199').split('-')]
top = int(opt.get('top', 20))
model = opt.get('model', '/scratch/jamesh/projects/amber/1aho_hurry_test/build_s68_PD/refined.pdb')
PDBSCALE = 0.01

# crystal model: protein heavy atoms (4 copies in the P1 cell)
cell = None; prot = []; plab = []
for l in open(model):
    if l.startswith('CRYST1'):
        cell = np.array([float(x) for x in l[6:33].split()])
    if l.startswith(('ATOM', 'HETATM')) and l[17:20] not in ('HOH', 'WAT') and l[76:78].strip() != 'H':
        r = int(l[22:26])
        if r > 256: continue
        prot.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        plab.append('%s %s%d' % (l[12:16].strip(), l[17:20].strip().title(), (r - 1) % 64 + 1))
prot = np.array(prot)


def near_label(x):
    d = prot - x; d -= cell * np.round(d / cell)
    r = np.sqrt((d ** 2).sum(1)); j = r.argmin()
    return plab[j], r[j]


def key_of(atom, resn, resi):
    if resn in ('HOH', 'WAT'):
        return None
    if resn in ('ACY', 'NH4', 'ACT'):
        return '%s %s' % (atom, resn)
    rn = resn.title()
    if rn in ('Hip', 'Hid', 'Hie', 'Chi', 'Chie', 'Chid'): rn = 'His'
    return '%s %s%d' % (atom, rn, (resi - 1) % 64 + 1)


def run_stats(rd):
    E = {}; W = {}; D = {}; nit = 0
    wcache = {}
    for i in range(i0, i1 + 1):
        f = os.path.join(rd, 'restraints_challenge_vs_weight_%d.txt' % i)
        rf = os.path.join(rd, 'restraints_for_%d.pdb' % i)
        if not os.path.exists(f) or not os.path.exists(rf):
            continue
        nit += 1
        # water reference points this iteration -> label
        wpos = {}
        for l in open(rf):
            if l.startswith(('ATOM', 'HETATM')) and l[17:20] in ('HOH', 'WAT'):
                wpos[int(l[22:26])] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            elif l.startswith(('ATOM', 'HETATM')):
                k = key_of(l[12:16].strip(), l[17:20].strip(), int(l[22:26]))
                W.setdefault(k, []).append(float(l[60:66]) * PDBSCALE)
        for l in open(f):
            p = l.split()
            if len(p) < 7: continue
            c, d, w = float(p[0]), float(p[1]), float(p[2])
            atom, resn, resi = p[3], p[4], int(p[6])
            k = key_of(atom, resn, resi)
            if k is None:
                x = wpos.get(resi)
                if x is None: continue
                ck = tuple(np.round(x, 1))
                if ck not in wcache: wcache[ck] = near_label(x)
                lab, r = wcache[ck]
                k = 'HOH@%s' % lab
                W.setdefault(k, []).append(w)
            E[k] = E.get(k, 0) + c
            D.setdefault(k, []).append(d)
    return E, W, D, nit


runs = {}
for rd in args:
    name = os.path.basename(os.path.abspath(rd)).replace('ff19_opc_s68_', '')
    E, W, D, nit = run_stats(rd)
    if nit == 0:
        print('%s: no iterations in %d-%d' % (name, i0, i1)); continue
    runs[name] = dict(E={k: v / nit / 4 for k, v in E.items()}, W={k: np.mean(v) for k, v in W.items()},
                      D={k: np.mean(v) for k, v in D.items()}, nit=nit,
                      tot=sum(E.values()) / nit)
    print('%-6s itr %d-%d (%d itr): total restraint energy %.1f per iteration (whole cell)' %
          (name, i0, i1, nit, runs[name]['tot']))

names = list(runs)
print('\nTop %d restraints per run by mean energy w d^2 per copy (waters: HOH@nearest protein atom, summed over the'
      ' waters with that label)' % top)
for n in names:
    E = runs[n]['E']
    print('  %s:' % n)
    for k in sorted(E, key=E.get, reverse=True)[:top]:
        print('    %-22s E %6.3f  w %6.2f  d %5.2f' % (k, E[k], runs[n]['W'].get(k, np.nan), runs[n]['D'].get(k, np.nan)))

# union of each run's top restraints, side by side
keys = []
for n in names:
    E = runs[n]['E']
    for k in sorted(E, key=E.get, reverse=True)[:top]:
        if k not in keys: keys.append(k)
if opt.get('focus', '1') == '1':
    for k in ['NZ Lys2', 'O His64', 'OXT His64', 'N His64', 'O Thr57', 'OH Tyr49', 'O His54',
              'CG His54', 'ND1 His54', 'CD2 His54', 'CE1 His54', 'NE2 His54',
              'CG His64', 'ND1 His64', 'CD2 His64', 'CE1 His64', 'NE2 His64',
              'HOH@ND1 His54', 'HOH@NE2 His54', 'HOH@CE1 His54', 'HOH@CD2 His54',
              'HOH@ND1 His64', 'HOH@NE2 His64', 'HOH@CE1 His64', 'HOH@CD2 His64', 'HOH@O His64', 'HOH@OXT His64',
              'HOH@NZ Lys2']:
        if k not in keys: keys.append(k)
print('\nEnergy per copy (E) and weight (w), side by side: %s' % '  '.join(names))
print('  %-22s' % 'restraint' + ''.join('%15s' % n for n in names))
for k in keys:
    row = ''
    for n in names:
        e = runs[n]['E'].get(k, 0.0); w = runs[n]['W'].get(k, np.nan)
        row += '%8.3f /%5.2f' % (e, w) if np.isfinite(w) or e else '%15s' % '-'
    print('  %-22s%s' % (k, row))
