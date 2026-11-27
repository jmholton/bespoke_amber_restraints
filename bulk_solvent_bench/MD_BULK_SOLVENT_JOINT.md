# Which MD strategies agree best with the 1aho data: both sessions' measures side by side (2026-11-24)

Two Claude sessions score the same 1aho MD runs in different ways. This file puts their results together so both work
from the same picture. Details: `MD_BULK_SOLVENT_FINDINGS.md` and `MD_BULK_SOLVENT_STRATEGY.md` (bulk-solvent session,
here) and `WEIGHTS_AND_B.md` (restraint-weight / B-factor session, bespoke_amber_restraints repo root).

## The two measures

- **Pipeline R** (weights/B session): the R factor of each iteration's MD ensemble average (protein + solvent, every
  atom at the run's own B: uniform 5 or 10, or per-atom) against Fobs after scaleit, from the run's
  `fofc_Rplot.txt`, averaged over the iterations that make up the solvent session's map window. Iteration-to-
  iteration sd 0.5-1.8.
- **R-free with the MD solvent map** (bulk-solvent session): the run's solvent-only average, rendered at uniform B 10,
  used as phenix's bulk solvent for the example model; "none" before any water is built, "r3" after three rounds of
  water building. Noise ~0.002.

Pipeline R scores the whole model, protein included, with whatever B the run uses. The R-free test scores only the
solvent map, on a refined protein. They answer different questions, so agreement between them is informative.

"B in the render" below is the B in the file. Every pipeline render so far (`nc2mtz_gemmi.com`) also carries gemmi's
anti-aliasing blur, max(0, 7.43 - lowest B) A^2 at 0.965 A: uniform 5 renders at 7.4, and per-atom B reaching 0
adds 7.4 to every atom. A fitted overall B absorbs it, so pipeline R is unaffected; renders at uniform B 10 (the
R-free test) carry none. Details and the blur-free GPU render in WEIGHTS_AND_B.md (2026-11-25).

## Same runs, both measures

| run | what it is | pipeline R (window) | B in the render | R-free none | R-free r3 |
|---|---|---|---|---|---|
| 2019 map | all heavy atoms + 48 waters/ASU at k 10, 48 copies, 100 ns (Amber18, ff14SB) | not run through the pipeline | uniform 5 | 0.1051 | **0.1030** |
| longcalm | fb3mod, no decay, release 100 / repick 30 | 30.5% (itr 1000-2000) | per-atom | 0.1142-0.1157 | 0.1052-0.1071 |
| af_rel4 pool | 5 altfix replicates, Cys12 copy 4 relaxed | 30.1-30.7% (101-200) | per-atom | 0.1156-0.1214 | 0.1069-0.1074 |
| g_s005_hie | OPC, 0.05 M salt | 30.7% (101-200) | per-atom | - | 0.1069 |
| af_std | OPC, no decay, standard settings | 31.1% (101-200) | per-atom | 0.1205 | 0.1067 |
| b5sd1 | 2s2r, opt2 options, weight_scaledown 1 | 33.0% (101-200) | uniform 5 | 0.1186 | 0.1067-0.1076 |
| longconv | longcalm with decay 0.9 | 33.4% (1000-1790) | per-atom | 0.1201 | 0.1074 |
| b5 | 2s2r, opt2 options, weight_scaledown 0.9 | 33.5% (101-200) | uniform 5 | - | 0.1094-0.1099 |
| k1B | 2s2r, k 1 on every heavy atom + 484 water O, never re-picked | 42.0% (101-200) | per-atom | 0.1281 | 0.1125-0.1129 |
| CA-only (mult, loca, locw, ff19_opc) | 232 CA, k decayed to 0.02 | 52.5-55.9% | per-atom / field | 0.1386-0.1429 | 0.1162-0.1180 |

## Where the two measures agree

1. **Too little restraint fails.** CA-only at k 0.02 is worst by both. Its protein drifts ~1 A rms from the crystal
   positions, and the solvent organises around the wrong surface.
2. **Uniform k 1 with per-atom B (k1B) is next worst by both**, although it holds the protein closest (0.27 A).
3. **Holding restraints beats decaying them.** Within the 2s2r build, b5sd1 beats b5 by both measures (R 33.0 vs
   33.5%; r3 ~0.107 vs ~0.110). Between longcalm and longconv, pipeline R prefers no decay strongly (30.5 vs 33.4%);
   R-free prefers it before water building and ties after.
4. **The front group is density-weighted restraints on all protein heavy atoms plus density-picked waters, held**:
   longcalm, the af_rel4 pool, g_s005_hie, af_std and b5sd1 are within noise of each other by R-free after building
   (0.1067-0.1074).

## Where they differ, and why

- **B model.** Inside the front group, pipeline R separates per-atom-B runs (30-31%) from uniform-B 5 runs (33%).
  That spread is the B model, not the solvent. The R-free test renders every map at uniform B 10 and does not see
  it. Freeing per-atom B at frozen weights lowered pipeline R from 34.0 to 31.7% (WEIGHTS_AND_B.md); the run's own
  per-atom B gained the R-free test ~0.002 before water building and nothing after. So B matters for the protein
  fit, and not for the solvent map once waters are built.
- **The 2019 map is still best by R-free and has no pipeline R.** Its recipe differs in a 48-copy supercell, fixed
  k 10 on every heavy atom and on crystallographic waters, and 100 ns at constant volume. No 2026 run has tried those.
- **Chemistry is tested only by pipeline R.** 1aho's His64 is a C-terminal amide (WEIGHTS_AND_B.md). In 4 matched
  pairs the amide lowered pipeline R by 0.1-0.7 points (itr 112-268) and the restraint energy in every pair, and
  reproduced a crystal H-bond the carboxylate breaks. But those amide builds also had no salt: to neutralise the
  amide's +4 per cell the three NH4+ were removed and one acetate added (ACY x8), against 0.2 M salt in the
  carboxylate builds. The comparison is repeated at fixed 0.68 M salt below. Every map in the R-free table comes from
  a carboxylate build.

## The strategy both sessions now recommend for 1aho

ff19SB with its CMAP intact (cpptraj 7 for strips), OPC or SPC/E water, His64 as the C-terminal amide with the
crystallisation salt kept and only the counter-acetates changed (0.68 M ammonium acetate: 11 NH4 + 19 ACY per cell
with the amide, 11 + 15 with the carboxylate), His54 as HIP and His64 as HID (H omit maps; pending the charge tests,
which add 4 more ACY for HIP54), hydrated to a true ~0 bar (ntp=1 pressure readings), density-weighted restraints on all
protein heavy atoms plus density-picked waters, held without decay (`weight_scaledown=1`), new waters at B 50, the
B floor at 0, and a map averaged over 32 ns or more.

## Open question: which restraint class does the work (tests launched 2026-11-24)

The bulk-solvent session's suggested experiment, run on the amide build so it also tests the chemistry in the R-free
measure. One build (`1aho_hurry_test/build_2s2r_amide`), the same code, uniform B 5, 200 iterations, b5sd1's weight
settings (opt2 options with weight_scaledown 1) unless stated. Run directories `1aho_hurry_test/ff19_opc_rc_<arm>`:

| arm | restrained | how |
|---|---|---|
| `rc_dw_w` | density-weighted protein heavy atoms + density-picked waters | b5sd1 settings (the reference arm) |
| `rc_dw` | density-weighted protein heavy atoms, no waters | start list without HOH; `repick_hohscale=1000` so waters are never picked |
| `rc_bb_w` | backbone N, CA, C, O + density-picked waters, density-weighted | start list of those only; repick and release off |
| `rc_ca1` | CA only, fixed k 1, no waters | start list of CA at k 1; weights frozen (variant `fixed`) |
| `rc_k10` | every heavy atom + the density-sited waters, fixed k 10 (the 2019 recipe in 4 copies) | start list at k 10 (capped 9.9999); weights frozen |

Score each by both measures: pipeline R from its `fofc_Rplot.txt`, and the R-free test on itr 101-200.

**All five finished itr 200 (2026-11-25).** Pipeline R and restraint energy, mean over itr 101-200 (sd = iteration
to iteration; one run per arm, and successive iterations are correlated, so treat 1-point gaps as suggestive):

| arm | restrained | pipeline R | sd | restraint E | R-free none / r3 |
|---|---|---|---|---|---|
| `rc_bb_w` | backbone + waters | **31.44%** | 0.29 | 322 | pending |
| `rc_dw` | all heavy atoms, no waters | 32.41% | 0.53 | 122 | pending |
| `rc_dw_w` | all heavy atoms + waters (reference) | 32.66% | 0.39 | 132 | pending |
| `rc_ca1` | CA only, k 1 | 34.85% | 0.42 | 92 | pending |
| `rc_k10` | everything, fixed k 10 | 40.27% | 0.10 | 867 | pending |

By pipeline R, free side chains fit better than restrained ones (backbone-only beats all-heavy-atom by 1.2 points;
its restraint energy is higher because many backbone atoms sit at the k 10 cap), restrained waters change nothing
within noise, and CA-only at k 1 is 2.2 points worse. Fixed k 10 on everything is worst, but at uniform B 5 that is
largely the B model: held that tightly the protein renders too sharp, which pipeline R cannot separate from the
solvent. The 2019 map (the same recipe in 48 copies) is still best by R-free, so the R-free test decides `rc_k10`.

The 8 runs
`ff19_opc_{cooh,amide}_{k1,rel,rel2,b5}` (amide vs carboxylate, salt-confounded) were stopped on 2026-11-25 at itr
296-344; the 0.68 M charge tests below replace them.

The restraint-class arms run on that ACY x8, no-NH4 amide build, so they compare restraint classes at no added salt.

## Charge tests at the crystallisation salt (launched 2026-11-25)

The crystal grew from salt, and the 1aho example puts it at 0.68 M ammonium acetate (11 NH4+ / acetate pairs per cell
from `add_salt_runme.com`). For reference, the 2019 map's cell held 84 acetate + 36 ammonium in 48 copies, i.e. 7 + 3
per 4 copies, the same as our 0.2 M builds. These arms share one 0.68 M placement; only the counter-acetates follow
each arm's charge. They test the C-terminus (amide vs carboxylate) without the salt confound, and the His tautomers
from H omit maps of the phenix bxfield gb3 model, refined omit plus add-H controls: His54 is doubly protonated
(HIP: HD1 1.0, HE2 0.67 of an amide H, null sd 0.33), His64 is HID (HE2 absent, HD1 about 0.6). R-free cannot see
the ring N-H; the difference density can. Every 2026 MD build before these used HIE54 / HIE64.
All run as `rc_dw_w` (b5sd1 settings, `min_align_weight=0`, uniform B 5, 200 itr); runs `1aho_hurry_test/ff19_opc_s68_<arm>`:

| arm | C-terminus | His54 | His64 | NH4 / ACY per cell |
|---|---|---|---|---|
| `cEE` | carboxylate | HIE | HIE | 11 / 15 |
| `EE` | amide | HIE | HIE | 11 / 19 |
| `ED` | amide | HIE | HID | 11 / 19 |
| `PD` | amide | HIP | HID | 11 / 23 |

`PD` is the omit maps' answer and `EE` the old default. `DD` (HID54) and `MD` (His54 HIP/HID mixture) were built but
not run once the refined omit maps came in.

Score each by both measures over itr 101-200. Expect R to move little for His54 (one surface residue); the local
measures (Fo-Fc and restraint weights at the ring atoms, waters at each ring N) and the R-free test are more telling.

**His54's ring orientation is right in every MD build.** The refinement flip test put His54's ring in the M
orientation (χ2 −79 to −119°) and found best_so_far's χ2 > 0 altlocs mislabelled. Every 2026 MD build starts from the
1aho model's single conformer, which is M (χ2 −110°), and the ring never flips in the MD (all 4 copies, PD and EE,
itr 1-107: χ2 −99 to −104° on average, no frame > 0). His64 stays in its one orientation (χ2 about +62°).
(`amide/his_chi.py`, `amide/his_chi_nc.py`; `amide/` is `/scratch/jamesh/projects/amber/amide`, copied to
bespoke_amber_restraints `bulk_solvent_bench/amide/`.)

### Controls (launched 2026-11-25 21:43, voltron)

Run exactly as the charge-test arms (rc_dw_w settings, 200 itr); runs `1aho_hurry_test/ff19_opc_s68_<arm>`
(`amide/launch_s68_controls.sh`):

| arm | what | why |
|---|---|---|
| `PDp` | PD with His54's ring relabelled to the P orientation in all 4 copies (`amide/setup_s68_flip.sh`: ND1↔CD2 and CE1↔NE2 trade coordinates) | positive control: known-wrong ring labels. The restraint picker already treats ND1/CD2 and CE1/NE2 as one atom, so only the chemistry differs. Can the MD measures see a ring N/C label at all? |
| `DD` | HID54, HID64 amide (built 2026-11-25, not run before) | completes HID / HIE / HIP at His54 |
| `PD_r2` | replicate of PD (same start, new random seeds) | run-to-run scatter of every measure |
| `ED_r2` | replicate of ED | the same, for the arm that differs from PD only at His54 |

New local measure (`amide/his_fofc.py`): each iteration's Fo-Fc (Fobs vs the render, fitted k and B, model phases)
averaged over iterations, read at the crystal model's His ring atoms and ring-H sites, plus the rms of Fo-Fc within
1.0 Å of each ring (lower = better local fit, independent of the labels).

Early look, itr 88-107, one run per arm (block-to-block sd 0.02-0.13), His ring local Fo-Fc rms in units of the map
rms:

| arm | His54 | His64 |
|---|---|---|
| `cEE` | 1.06 | 1.25 |
| `EE` | 1.30 | 0.91 |
| `ED` | 1.33 | 1.05 |
| `PD` | 1.50 | 0.93 |

His64 fits worst as the carboxylate, as expected. His54 fits worst as HIP so far, the opposite of the omit maps. The
replicates decide whether gaps of this size are real.

### Worst restraints moved with protonation (2026-11-26)

Pipeline R is not the measure for one residue's chemistry; the restraints are. Per iteration optimize_weights logs
each restrained atom's displacement d from its reference point and its weight w (`restraints_challenge_vs_weight_N.txt`);
w d² is its restraint energy. `amide/worst_restraints.py` averages that over iterations and the 4 copies (waters
labelled by the nearest protein atom to their reference point). Energy per copy, itr 101-199 (the leaders in every arm,
about 0.4-0.6, are Arg62 CZ, Arg18 NH1, Asp8 OD2, Lys58 O, and do not depend on protonation):

| restraint | `cEE` | `EE` | `ED` | `PD` |
|---|---|---|---|---|
| His64 O | 0.37 | 0.00 | 0.02 | 0.07 |
| Thr57 O | 0.20 | 0.06 | 0.03 | 0.06 |
| water at Cys63 O | 0.52 | 0.13 | 0.08 | 0.15 |
| water at Asp8 OD2 | 0.50 | 0.14 | 0.01 | 0.03 |
| Lys2 NZ | 0.49 | 0.43 | 0.29 | 0.38 |
| His54 NE2 | 0.00 | 0.07 | 0.04 | 0.46 |
| His54 CE1 | 0.00 | 0.05 | 0.02 | 0.31 |
| His64 NE2 | 0.16 | 0.25 | 0.32 | 0.51 |
| His64 CE1 | 0.20 | 0.19 | 0.29 | 0.38 |
| total, whole cell | 151 | 135 | 137 | 135 |

Checked in all eight runs over itr 101-119 (controls included):

1. The amide fixes the C-terminal junction: His64 O ≤ 0.10 in all 7 amide runs vs 0.50 in `cEE`; Thr57 O and three
   nearby waters likewise. Lys2 NZ improves only partly (0.24-0.46, weight 3.5-5; lowest in `ED` and `ED_r2`).
2. HIP54 puts His54's outer ring edge (NE2, CE1) among the worst restraints: 0.41-0.46 in `PD`, `PD_r2` and `PDp`
   vs ≤ 0.09 for HIE or HID54; weights rise from about 0.05 to 1.1-1.9 and the edge sits about 0.5 Å off its density
   point. `PDp` shows it too, so it comes from the charge, not the ring labels. The charged ring draws more water
   (`amide/his_neighbours_md.py`: 1.5 waters within 3.2 Å of NE2 vs 1.0 for HIE; 0.5 vs 0.2 at CE1); acetate touches
   it in < 4 % of frames. So in ff19SB a fully charged His54 does not sit where the density puts it; the omit maps' HE2
   of 0.6-0.7 allows a partly protonated His54.
3. HID64 strains His64's ring edge, most with HIP or HID54: His64 NE2 0.24-0.35 with HIE54 (`ED`, `ED_r2`), 0.54
   with HID54 (`DD`), 0.67-0.73 with HIP54 (`PD`, `PD_r2`), vs 0.22-0.25 with HIE64. His64's contacts do not change
   between arms; the cause is open.

Two crystal waters sit on His64 ND1 (3.04 and 3.12 Å, full occupancy, density rank 17 and 21 of 121 waters) and are
held at weights 4-7. The refinement's partial "waters" 1.9-2.4 Å from His64 ND1 in best_so_far are not in the 1aho
model; at those distances they could be ions (JH), which would also bear on His64's HD1 signal.

**Restraint release.** These runs release every iteration (`release_itr=1`): restraints at the weight cap (> 9, up to
10 per iteration) are reset to 0.01 with their symmetry mates and neighbours (shells 2.2 + 1.8 + 1.8 + 1.8 Å), then
regrow. In practice about 3 seeds and 65 atoms per iteration, the seeds mostly disulfide SG (Cys16, 22, 46, 48), then
His64's terminal O in the amide arms and Ser33 N; water is 0-1 % of what is released (`amide/release_stats.py`).
Waters held at weights 6-8 never reach the cap, so they are never released. Water restraints change only by repick
(every iteration, within 2 Å), teleport (one water per iteration) and hydration.

### More His protonation arms (launched 2026-11-26)

Complete His54 {HIE, HID, HIP, HIP/HID mixture} × His64 {HIE, HID}, all amide at 0.68 M, rc_dw_w settings
(`amide/setup_s68_more.sh`, `amide/launch_s68_controls.sh`):

| arm | His54 | His64 | NH4 / ACY per cell |
|---|---|---|---|
| `PE` | HIP | HIE | 11 / 23 |
| `DE` | HID | HIE | 11 / 19 |
| `ME` | HIP in copies 1+3, HID in 2+4 | HIE | 11 / 21 |
| `MD` | HIP in copies 1+3, HID in 2+4 | HID | 11 / 21 |

With `EE`, `ED`, `DD`, `PD` that is every cell. The mixtures stand for the omit maps' partial HE2 at His54 (HD1 1.0,
HE2 0.6-0.7). Salt concentration is to be varied once protonation is settled (JH).
