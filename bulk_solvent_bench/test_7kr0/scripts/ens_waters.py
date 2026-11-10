# usage: ens_waters.py model.pdb maps.mtz out.pdb [levels=6,5,4,3,2.5] [min_new=5] [rom=2.6] [rom_mode=own]
#                      [clash=2.3] [mindist=1.0] [maxdist=4.0] [B=15] [free=1]
# Water building for a multi-conformer ensemble (7KR0 S0112 style: every protein atom has an altloc, conformers at
# equal occupancy), conservatively: peaks of the work-only mFo-DFc map (FOFCWT/PHFOFCWT, free set zeroed) are taken
# at the highest of the sigma levels that still gives at least min_new waters, so lower levels are used only once the
# higher ones have run out (the lowest level takes whatever it gives).  Peaks go tallest first, each moved to its
# symmetry image nearest the model.  A peak is skipped if it is within mindist of any atom (that atom's density).
# Conformer L can hold the water if no protein heavy atom of L or blank is within rom, no heavy atom of L or blank (old
# waters, waters placed before it) within clash, and one (protein or water) within maxdist.  rom_mode=any is the
# range-of-motion rule of 7KR0_SETTINGS.md instead: no water within rom of a protein atom of ANY conformer (on a
# dry model it drops the conformer-specific first shell: 108 waters at 6 sigma against 173).  Held by every
# conformer: one blank water at occupancy 1; otherwise one altloc per conformer that can hold it, each at that
# conformer's share (the held-water rule).  Old atoms pass through, ANISOU records included.  Prints
# "level <sigma>" (or "level none" and exit status 3 when nothing can be placed).
import sys
import numpy as np
import gemmi
from scipy import ndimage

args = [a for a in sys.argv[1:] if "=" not in a]
opt = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
pdbin, mtzin, pdbout = args
levels = sorted((float(s) for s in opt.get("levels", "6,5,4,3,2.5").split(",")), reverse=True)
min_new = int(opt.get("min_new", 5))
rom = float(opt.get("rom", 2.6))
rom_mode = opt.get("rom_mode", "own")
clash = float(opt.get("clash", 2.3))
mindist = float(opt.get("mindist", 1.0))
maxdist = float(opt.get("maxdist", 4.0))
Bnew = float(opt.get("B", 15))

lines = open(pdbin).read().splitlines()
atoms = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
st = gemmi.read_structure(pdbin)
cell = st.cell
ops = st.find_spacegroup().operations()
rot = [np.array(op.rot) / op.DEN for op in ops]
tran = [np.array(op.tran) / op.DEN for op in ops]
F = np.array(cell.frac.mat.tolist())
O = np.array(cell.orth.mat.tolist())

heavy = [l for l in atoms if l[76:78].strip() not in ("H", "D")]
fc0 = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in heavy]) @ F.T
alt0 = np.array([l[16] for l in heavy])
isprot = np.array([l[17:20] not in ("HOH", "WAT") for l in heavy])
letters = sorted(set(alt0) - {" "})
occ = {L: np.median([float(l[54:60]) for l in heavy if l[16] == L]) for L in letters}
# phenix writes occupancy to 2 decimals (1/6 -> 0.17): normalize the conformers to sum to 1
occ = {L: occ[L] / sum(occ.values()) for L in letters}
print("conformers %s, occupancy %s" % ("".join(letters), " ".join("%.4f" % occ[L] for L in letters)))


def symdist(f, others):
    best = np.full(len(others), 1e9)
    for r, t in zip(rot, tran):
        d = r @ f + t - others
        d -= np.round(d)
        best = np.minimum(best, np.linalg.norm(d @ O.T, axis=1))
    return best


def nearest_image(f):
    bestd, bestf = 1e9, f
    for r, t in zip(rot, tran):
        g = r @ f + t
        d = g - fc0
        sh = np.round(d)
        dd = np.linalg.norm((d - sh) @ O.T, axis=1)
        k = int(np.argmin(dd))
        if dd[k] < bestd:
            bestd, bestf = dd[k], g - sh[k]
    return bestf


mtz = gemmi.read_mtz_file(mtzin)
# work reflections only: zero the free set's mFo-DFc (free=0 keeps them).  The free set is the minority value of
# the flag column phenix copies into its output mtz
if opt.get("free", "1") != "0":
    data = np.array(mtz, copy=True)
    labels = [c.label for c in mtz.columns]
    fl = data[:, labels.index(opt.get("free_label", "R-free-flags"))]
    vals, cnt = np.unique(fl[~np.isnan(fl)], return_counts=True)
    test = vals[np.argmin(cnt)]
    isfree = fl == test
    data[isfree, labels.index("FOFCWT")] = 0
    mtz.set_data(data)
    print("left the %d free reflections (flag %g) out of the difference map" % (isfree.sum(), test))
grid = mtz.transform_f_phi_to_map("FOFCWT", "PHFOFCWT", sample_rate=3)
a = np.array(grid, copy=False)
rms = a.std()
mx = ndimage.maximum_filter(a, size=3, mode="wrap")
idx = np.argwhere((a == mx) & (a > levels[-1] * rms))
h = a[tuple(idx.T)]
order = np.argsort(-h)
peaks = [(h[i] / rms, nearest_image(idx[i] / np.array(a.shape))) for i in order]
print("%d grid peaks above %.1f sigma (rms %.4f); tallest %.1f sigma" %
      (len(peaks), levels[-1], rms, peaks[0][0] if peaks else 0))


def place(level):
    fc, alt = fc0.copy(), alt0.copy()
    new, newf = [], []
    counts = dict(blank=0, split=0, atom=0, rom=0, none=0, dup=0)
    for height, f in peaks:
        if height < level:
            break
        if newf and symdist(f, np.array(newf)).min() < mindist:
            counts["dup"] += 1
            continue
        d = symdist(f, fc)
        if d.min() < mindist:
            counts["atom"] += 1
            continue
        dp, ap = d[:len(fc0)][isprot], alt0[isprot]
        if rom_mode == "any" and dp.min() < rom:
            counts["rom"] += 1
            continue
        ok = [L for L in letters if clash <= d[(alt == L) | (alt == " ")].min() <= maxdist
              and dp[(ap == L) | (ap == " ")].min() >= rom]
        if not ok:
            counts["none"] += 1
            continue
        newf.append(f)
        if len(ok) == len(letters):
            counts["blank"] += 1
            confs = [(" ", 1.0)]
            fc = np.vstack([fc, f])
            alt = np.append(alt, " ")
        else:
            counts["split"] += 1
            confs = [(L, occ[L]) for L in ok]
            fc = np.vstack([fc] + [f] * len(ok))
            alt = np.append(alt, ok)
        new.append((O @ f, height, confs))
    return new, counts


for level in levels:
    new, counts = place(level)
    print("at %.1f sigma: %d shared, %d conformer-specific waters; skipped %d on atoms, %d within %.1f A of protein "
          "(rom_mode=any), %d no conformer could hold, %d duplicates" % (level, counts["blank"], counts["split"],
          counts["atom"], counts["rom"], rom, counts["none"], counts["dup"]))
    if len(new) >= min_new or level == levels[-1]:
        break
if not new:
    print("level none")
    sys.exit(3)
print("level %g" % level)

# old atoms as they are; new waters in chain W, numbered on from the highest water number there
oldw = [int(l[22:26]) for l in atoms if l[17:20] == "HOH" and l[21] == "W"]
n0 = max(oldw, default=0)
out = [l for l in lines if l.startswith(("CRYST1", "SCALE"))] + \
    [l for l in lines if l.startswith(("ATOM", "HETATM", "ANISOU", "TER"))]
serial = len(atoms)
for i, (x, height, confs) in enumerate(new):
    for L, o in confs:
        serial += 1
        out.append("HETATM%5d  O  %sHOH W%4d    %8.3f%8.3f%8.3f%6.4f%6.2f           O"
                   % (serial % 100000, L, n0 + i + 1, x[0], x[1], x[2], o, Bnew))
out.append("END")
open(pdbout, "w").write("\n".join(out) + "\n")
nalt = sum(len(c) for _, _, c in new)
print("wrote %s: %d new waters at >= %.1f sigma (%d atoms), occupancy sum %.1f" %
      (pdbout, len(new), level, nalt, sum(o for _, _, c in new for _, o in c)))
