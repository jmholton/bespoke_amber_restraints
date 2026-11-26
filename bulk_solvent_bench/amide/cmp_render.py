"""Compare MD-average renders (FCavg/PHICavg) from nc2mtz_gemmi.com and nc2mtz_gpu.com (2026-11-25).
usage: ccp4-python cmp_render.py ref.mtz test.mtz [more.mtz ...] [fobs=refme.mtz]
For each test file against the first: least-squares scale k (test*k ~ ref), amplitude R after k, vector R on the
complex F, the F-weighted mean phase difference and CC(F), overall and in 5 resolution shells.  With fobs=, each file's
R against Fobs (P1 FP, looked up at the same or a P212121-equivalent index) after a fitted k and overall B."""
import sys, numpy as np, gemmi

args = [a for a in sys.argv[1:] if '=' not in a]
opt = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)


def load(fn, fcol='FCavg', pcol='PHICavg'):
    m = gemmi.read_mtz_file(fn)
    labs = [c.label for c in m.columns]
    if fcol not in labs:
        fcol = [c.label for c in m.columns if c.type == 'F'][0]
        pcol = [c.label for c in m.columns if c.type == 'P'][0]
    hkl = m.make_miller_array()
    F = m.column_with_label(fcol).array.astype(float)
    P = np.radians(m.column_with_label(pcol).array.astype(float))
    d = m.make_d_array()
    ok = np.isfinite(F) & np.isfinite(P)
    return {tuple(h): (f, p, dd) for h, f, p, dd in zip(hkl[ok], F[ok], P[ok], d[ok])}


def compare(a, b, label):
    keys = sorted(set(a) & set(b))
    Fa = np.array([a[k][0] for k in keys]); Pa = np.array([a[k][1] for k in keys])
    Fb = np.array([b[k][0] for k in keys]); Pb = np.array([b[k][1] for k in keys])
    d = np.array([a[k][2] for k in keys])
    k = (Fa * Fb).sum() / (Fb * Fb).sum()
    print('%s: %d common reflections, scale k = %.5f (1/k = %.4f)' % (label, len(keys), k, 1 / k))
    edges = np.quantile(1 / d ** 2, np.linspace(0, 1, 6))
    rows = [('all', np.ones(len(d), bool))]
    for i in range(5):
        s = (1 / d ** 2 >= edges[i]) & (1 / d ** 2 <= edges[i + 1])
        rows.append(('%.2f-%.2f A' % (1 / np.sqrt(edges[i + 1]), 1 / np.sqrt(edges[i])), s))
    print('   %-14s %8s %8s %8s %9s %7s' % ('shell', 'n', 'R(F)%', 'Rvec%', 'dphi_deg', 'CC(F)'))
    for name, s in rows:
        fa, fb, pa, pb = Fa[s], k * Fb[s], Pa[s], Pb[s]
        R = np.abs(fa - fb).sum() / fa.sum() * 100
        Rv = np.abs(fa * np.exp(1j * pa) - fb * np.exp(1j * pb)).sum() / fa.sum() * 100
        dp = np.degrees(np.abs(np.angle(np.exp(1j * (pa - pb)))))
        print('   %-14s %8d %8.3f %8.3f %9.3f %7.5f' % (name, s.sum(), R, Rv, (dp * fa).sum() / fa.sum(),
                                                       np.corrcoef(fa, fb)[0, 1]))


def vs_fobs(fobs_fn, calc, label):
    m = gemmi.read_mtz_file(fobs_fn)
    hkl = m.make_miller_array(); F = m.column_with_label('FP').array.astype(float)
    obs = {tuple(h): f for h, f in zip(hkl, F) if np.isfinite(f)}
    pairs = []
    for (h, k, l), (fc, p, d) in calc.items():
        for e in ((h, k, l), (-h, k, l), (h, -k, l), (-h, -k, l), (-h, -k, -l), (h, k, -l), (-h, k, -l), (h, -k, -l)):
            if e in obs:
                pairs.append((obs[e], fc, 1 / d ** 2)); break
    fo, fc, s2 = map(np.array, zip(*pairs))
    best = None
    for B in np.arange(-20, 20.01, 0.1):
        c = fc * np.exp(-B * s2 / 4)
        k = (fo * c).sum() / (c * c).sum()
        R = np.abs(fo - k * c).sum() / fo.sum()
        if best is None or R < best[0]:
            best = (R, k, B)
    print('%-28s vs Fobs: %d refl, R %.3f%% (k %.4g, B %+.1f)' % (label, len(fo), best[0] * 100, best[1], best[2]))


ref = load(args[0])
maps = {fn: load(fn) for fn in args[1:]}
for fn, m in maps.items():
    compare(ref, m, '%s vs %s' % (fn, args[0]))
if 'fobs' in opt:
    for fn, m in [(args[0], ref)] + list(maps.items()):
        vs_fobs(opt['fobs'], m, fn)
