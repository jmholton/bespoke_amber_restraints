# usage: close_waters.py context.pdb new_water.pdb mindist clashdist > waterconf_in.pdb
# context.pdb: CRYST1 and every atom a new water must be checked against (old waters renamed WAT).
# Each new water moves to its symmetry image nearest the model; one within mindist of any context atom (or of a
# new water kept before it) is dropped.  Out: "DIST i j d" lines for every pair closer than clashdist, symmetry
# mates included, numbered the way waterconf.awk counts atoms; then CRYST1, the context atoms and the kept waters.
import sys
import numpy as np
import gemmi
ctxfile, newfile = sys.argv[1:3]
mindist, clashdist = float(sys.argv[3]), float(sys.argv[4])
def atoms(fn):
    return [l.rstrip("\n") for l in open(fn) if l.startswith(("ATOM", "HETATM"))]
cryst = [l.rstrip("\n") for l in open(ctxfile) if l.startswith("CRYST1")][0]
cell = gemmi.UnitCell(*[float(cryst[i:j]) for i, j in ((6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54))])
ops = gemmi.find_spacegroup_by_name(cryst[55:66].strip()).operations()
rot = [np.array(op.rot) / op.DEN for op in ops]
tran = [np.array(op.tran) / op.DEN for op in ops]
F = np.array(cell.frac.mat.tolist())
O = np.array(cell.orth.mat.tolist())
def xyz(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines]).reshape(-1, 3)
ctx = atoms(ctxfile)
fc = xyz(ctx) @ F.T
def images(f):
    # every symmetry image of fractional position f
    return [r @ f + t for r, t in zip(rot, tran)]
def symdist(f, others):
    # shortest distance from fractional f to each of others (fractional), over symmetry and lattice
    best = np.full(len(others), 1e9)
    for g in images(f):
        d = g - others
        d -= np.round(d)
        best = np.minimum(best, np.linalg.norm(d @ O.T, axis=1))
    return best
kept, fk, dist = [], [], []
nc = len(ctx)
for l in atoms(newfile):
    f = xyz([l])[0] @ F.T
    # the image nearest any context atom, moved next to it
    bestd, bestf = 1e9, f
    for g in images(f):
        d = g - fc
        sh = np.round(d)
        r = np.linalg.norm((d - sh) @ O.T, axis=1)
        k = int(np.argmin(r))
        if r[k] < bestd:
            bestd, bestf = r[k], g - sh[k]
    dc = symdist(bestf, fc)
    dn = symdist(bestf, np.array(fk).reshape(-1, 3))
    if dc.min(initial=1e9) < mindist or dn.min(initial=1e9) < mindist:
        continue
    n = nc + len(kept) + 1
    dist += ["DIST %d %d %.3f" % (i + 1, n, d) for i, d in enumerate(dc) if d < clashdist]
    dist += ["DIST %d %d %.3f" % (nc + i + 1, n, d) for i, d in enumerate(dn) if d < clashdist]
    x = O @ bestf
    kept.append("%s%8.3f%8.3f%8.3f%s" % (l[:30], x[0], x[1], x[2], l[54:]))
    fk.append(bestf)
print("\n".join(dist + [cryst] + ctx + kept))
