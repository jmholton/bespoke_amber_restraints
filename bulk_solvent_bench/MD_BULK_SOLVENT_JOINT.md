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
  pairs the amide lowered pipeline R by 0.1-1.2 points and the restraint energy in every pair, and reproduced a
  crystal H-bond the carboxylate breaks. Every map in the R-free table comes from a carboxylate build.

## The strategy both sessions now recommend for 1aho

ff19SB with its CMAP intact (cpptraj 7 for strips), OPC or SPC/E water, His64 as the C-terminal amide (salt rebalanced:
ACY x8, no NH4 per cell), hydrated to a true ~0 bar (ntp=1 pressure readings), density-weighted restraints on all
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

Score each by both measures: pipeline R from its `fofc_Rplot.txt`, and the R-free test on itr 101-200. The 8 runs
`ff19_opc_{cooh,amide}_{k1,rel,rel2,b5}` (amide vs carboxylate) continue to itr 500 alongside.
