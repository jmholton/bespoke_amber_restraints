# what changes in the w1_100ns solvent map between early and late production: (late - early) density map from
# the 10-ns block maps, its biggest peaks and the model atoms nearest them (symmetry-aware), plus the
# shell-by-shell CC and amplitude ratio late/early.  usage: drift_map.py model.pdb [npeaks]
import sys
import numpy as np
import gemmi
W = "/home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns/windows/"
early, late = ["map_1_10", "map_11_20"], ["map_81_90", "map_91_100"]
def coefs(names):
    tot = None
    for n in names:
        m = gemmi.read_mtz_file(W + n + ".mtz")
        f = m.column_with_label("FPART").array * np.exp(1j * np.radians(m.column_with_label("PHIPART").array))
        tot = f if tot is None else tot + f
    return m, tot / len(names)
me, fe = coefs(early)
ml, fl = coefs(late)
assert np.array_equal(me.make_miller_array(), ml.make_miller_array())
hkl = me.make_miller_array()
d = me.make_d_array()
print("shell       CC(late,early)  rms|F| late/early")
for hi, lo in [(99, 8), (8, 5), (5, 3.5), (3.5, 2.5), (2.5, 1.9), (1.9, 1.4), (1.4, 0.96)]:
    s = (d <= hi) & (d > lo)
    cc = np.real(np.sum(fl[s] * np.conj(fe[s]))) / np.sqrt(np.sum(abs(fl[s])**2) * np.sum(abs(fe[s])**2))
    print("%4g-%-5g %10.4f %12.3f" % (hi, lo, cc, np.sqrt(np.sum(abs(fl[s])**2) / np.sum(abs(fe[s])**2))))
def grid(f):
    a = gemmi.ComplexAsuData(me.cell, me.spacegroup, hkl, f.astype(np.complex64))
    return a.transform_f_phi_to_map(exact_size=[0, 0, 0], sample_rate=3.0)
ge, gl = grid(fe), grid(fl)
diff = np.array(gl, copy=False) - np.array(ge, copy=False)
rms = diff.std()
print("early map rms %.4f, late-early rms %.4f (%.0f%% of the early map's)" % (np.array(ge).std(), rms, 100 * rms / np.array(ge).std()))
st = gemmi.read_structure(sys.argv[1])
st.setup_entities()
ns = gemmi.NeighborSearch(st[0], st.cell, 6).populate()
npk = int(sys.argv[2]) if len(sys.argv) > 2 else 20
dg = gemmi.FloatGrid(diff.astype(np.float32), me.cell, me.spacegroup)
# local extrema on the grid above 4 rms; symmetry copies of a kept peak are skipped below
peaks = []
for sign in (1, -1):
    a = sign * diff
    m = a > 4 * rms
    for ax in range(3):
        for sh in (1, -1):
            m &= a >= np.roll(a, sh, axis=ax)
    for i, j, k in zip(*np.nonzero(m)):
        p = dg.get_position(int(i), int(j), int(k))
        peaks.append((a[i, j, k] / rms * sign, p))
seen, out = [], []
for h, p in sorted(peaks, key=lambda x: -abs(x[0])):
    if any(st.cell.find_nearest_image(q, p, gemmi.Asu.Any).dist() < 1.0 for q in seen):
        continue
    seen.append(p)
    marks = ns.find_atoms(p, '\0', radius=6)
    best = None
    for mk in marks:
        cra = mk.to_cra(st[0])
        r = st.cell.find_nearest_image(p, cra.atom.pos, gemmi.Asu.Any).dist()
        if best is None or r < best[0]:
            best = (r, "%s %s%d %s%s occ %.2f" % (cra.chain.name, cra.residue.name, cra.residue.seqid.num, cra.atom.name, cra.atom.altloc if cra.atom.altloc != '\0' else '', cra.atom.occ))
    out.append("%+6.1f rms  nearest %-28s %4.1f A" % (h, best[1] if best else "-", best[0] if best else 99))
    if len(out) >= npk: break
print("biggest late-early peaks (rms of the difference map), nearest model atom:")
print("\n".join(out))
# which residues the drift sits on: rms of the difference map within 1.5 A of each residue's side-chain heavy
# atoms (any conformer), and whether the MD left that side chain free (an atom's altlocs > 0.8 A apart)
mask_r = 1.5
print()
print("residues by rms(late-early) within %.1f A of side-chain heavy atoms, in units of the map's rms" % mask_r)
rows = []
for ch in st[0]:
    for res in ch:
        if res.name == "HOH": continue
        pos = {}
        for a in res:
            if a.element.name in ("H", "D") or a.name in ("N", "CA", "C", "O"): continue
            pos.setdefault(a.name, []).append(a.pos)
        if not pos: continue
        split = max(max(p.dist(q) for p in v for q in v) for v in pos.values())
        mk = gemmi.FloatGrid(dg.shape[0], dg.shape[1], dg.shape[2])
        mk.set_unit_cell(me.cell)
        mk.spacegroup = me.spacegroup
        for v in pos.values():
            for p in v:
                mk.set_points_around(p, mask_r, 1.0)
        mk.symmetrize_max()
        sel = np.array(mk, copy=False) > 0
        rows.append((np.sqrt(np.mean(diff[sel]**2)) / rms, "%s%d" % (res.name, res.seqid.num), split))
rows.sort(reverse=True)
free = [r for r in rows if r[2] > 0.8]
fixed = [r for r in rows if r[2] <= 0.8]
print("%d free side chains, mean %.2f; %d restrained, mean %.2f" % (len(free), np.mean([r[0] for r in free]), len(fixed), np.mean([r[0] for r in fixed])))
for r in rows[:15]:
    print("  %-8s %5.2f   altloc spread %4.1f A  %s" % (r[1], r[0], r[2], "free" if r[2] > 0.8 else "restrained"))
# fairer reference: share of the drift (sum of squared late-early density) near free side chains vs their share
# of the solvent volume (grid points farther than 1.5 A from every restrained or main-chain heavy atom)
def mask_of(atoms):
    mk = gemmi.FloatGrid(dg.shape[0], dg.shape[1], dg.shape[2])
    mk.set_unit_cell(me.cell)
    mk.spacegroup = me.spacegroup
    for p in atoms:
        mk.set_points_around(p, mask_r, 1.0)
    mk.symmetrize_max()
    return np.array(mk, copy=False) > 0
freeat, fixedat = [], []
for ch in st[0]:
    for res in ch:
        if res.name == "HOH": continue
        pos = {}
        for a in res:
            if a.element.name in ("H", "D"): continue
            pos.setdefault(a.name, []).append(a.pos)
        for name, v in pos.items():
            sc_free = name not in ("N", "CA", "C", "O") and max(p.dist(q) for p in v for q in v) > 0.8
            (freeat if sc_free else fixedat).extend(v)
fixm = mask_of(fixedat)
freem = mask_of(freeat) & ~fixm
solv = ~fixm
d2 = diff**2
print()
print("solvent volume (beyond %.1f A of restrained/main-chain atoms): %.0f%% of cell, %.0f%% of the drift" % (mask_r, 100 * solv.mean(), 100 * d2[solv].sum() / d2.sum()))
print("near free side chains: %.0f%% of that solvent volume, %.0f%% of its drift; rms there %.2f vs %.2f in the rest" % (
    100 * freem.sum() / solv.sum(), 100 * d2[freem].sum() / d2[solv].sum(),
    np.sqrt(d2[freem].mean()) / rms, np.sqrt(d2[solv & ~freem].mean()) / rms))
