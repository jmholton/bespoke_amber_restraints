# MD bulk-solvent maps: what we learned (1aho, 2026-10-31 to 2026-11-24)

A summary of the findings. The tests, tables and methods behind every number are in `MD_BULK_SOLVENT_STRATEGY.md`
(here), in date order; section names below point into it.

## The test

- **Model and data:** the phenix_user_solvent example model (`starthere.pdb`, nearly dry: 128 partial water atoms,
  occupancy sum 21.5) against `refme.mtz`. Each MD map goes in as phenix's user bulk solvent (FPART/PHIPART).
- **Refinement:** phenix 2.2-6143, 3 macro cycles of xyz + ADP + water occupancy, from 2 jiggled starts.
- **Water building:** then 3 rounds of the `altconf` copy of add_waters.com (`water_models/addwater/rounds_close.csh`).
- **Scores:** "none" is R-free before any water is built, "r3" after round 3. Both are means of the 2 starts.
- **Reproducibility:** on one host the test is deterministic. af_std's 16 R values reproduced exactly three times,
  including from the /scratch copy of the install.
- **Noise floor ~0.002:**
  - Replicate MD runs of the same settings differ by 0.001-0.002.
  - Water picks move a map's round 3 by up to 0.003. The same frames rendered at B 5 instead of B 10 move a round by
    up to 0.001.
  - So differences under ~0.002 mean nothing.
- **Programs:** since 2026-11-24 the scripts take amber, phenix and ccp4 from `programs.sourceme` / `programs.sh`
  (`$PROGRAMS`, default `/data/home/programs`). The af_std map and its refinement reproduce exactly from that tree.

## The runs, by name

Unless stated, a 200-iteration run's map is iterations 101-200, every 5th, 5 frames each: 100 frames over its last
50 ns, rendered at uniform B 10. All the `1aho_hurry_test` runs are single-cell (4 copies) runs of the bespoke
`optimize_weights` workflow, 0.5 ns per iteration, under `/scratch/jamesh/projects/amber/1aho_hurry_test/<run>/`.

| name | what it is |
|---|---|
| 2019 map | `1aho_test/j5w8760/solvent_RestAll.map`: Amber18 ff14SB + OPC, 2x2x3 supercell (48 copies), 100 ns at constant volume, 287 K. Every protein heavy atom and 2298 crystallographic water O (48 per ASU) restrained at k = 10. Rendered at B 5. As FPART: `water_models/example_check/md2019_map.mtz` |
| example map | the phenix_user_solvent example's `solvent_Fpart.mtz`: the 2019 map, rescaled, plus one round of squish editing |
| longcalm | fb3mod water, 11 NH4 / 15 acetate, His54 HID. 2000 iterations (1 us), no weight decay, release 100 / repick 30, per-atom B. Map: iterations 1000-2000, every 10th (505 frames). `hl_lcb5` is this map rendered at B 5, `hl_lcbf` at the run's own per-atom B |
| longconv | longcalm's sibling from the same seed, with weight decay 0.9; it stalled at iteration 1790 |
| af_std | one of the other session's "altfix" tests, launched 2026-10-18 after bespoke fix 40ba58d. It ran the Aug 29 `bespoke_master` code with standard settings ("std"), seeded from iteration 200 of noalt_opc_hie54 with the altloc reference points put back and the seed's water claims stripped. That seed's line: salt02 (0.2 M NH4 acetate = 3 NH4 / 7 acetate) -> salt02_opt -> OPC rebuild with His54 HIE (Oct 11) -> noalt_opc_hie54. OPC, ~700 waters, all four Cys12 copies strained, 200 iterations (100 ns, Oct 18 20:45 to Oct 19 11:48), no weight decay, netfrc=0 |
| af_rel4 pool | five near-replicates of the altfix test seeded with copy 4's Cys12 relaxed: af_rel4, maw0, sxyz, bcsx, altB1. They differ only in alignment and claiming options. Their second halves were pooled (5 x 50 ns) |
| g_s005_hie | a fresh OPC / His HIE build at 0.05 M salt (0 NH4 / 4 acetate), 200 iterations |
| b5, b5sd1, k1B | the `ff19_opc_*_2s2r` runs (Nov 4). The salt02 system with the Cys12/Cys63 B conformers grafted onto copies 3 and 4 (2 strained + 2 relaxed), ff19SB + OPC with CMAP intact, restraints reset to k = 1. b5 = the example's opt2 options (weight_scaledown 0.9) plus render B 5. b5sd1 = b5 with weight_scaledown 1. k1B = k = 1 on every heavy atom and 484 water O, never re-picked, with per-atom B |
| CA-only | `ff19_opc_mult`, `_loca`, `_locw` (build_ff19_opc, Nov 1, ff19SB + OPC, CMAP intact). 232 CA restrained, leaving a 5 A hole around Cys12. k decays from ~1 to 0.02 by iteration 113, then holds. Nothing else restrained. The three differ only in the run's own B handling. `ff19_opc` is their parent, run to iteration 751 |
| w1, n1, ... | `bulk_solvent_runme.com` maps (OPC unless named, 2 x 1 ns, protein k = 1 x rho per sigma, cap 10, split altlocs free, no water restraints). w1 is the default. n1 is w1 hydrated to a true 0 bar (pressure branch, `ntp=1`). pp500 / pp1500 / pp3000 / pm500 are pressure targets on the `ntp=4` scale. w01-w10 and flat10 vary the restraint weight. wm_* compare 8 water models over 10 ns |

## Findings

### 1. The 2019 map is still the best solvent map

| map | none | r3 |
|---|---|---|
| 2019 map | 0.1051 | **0.1030** |
| example map (from the 2019 map) | 0.1062 | 0.1035 |
| longcalm at B 5 (hl_lcb5) | 0.1157 | 0.1052 |
| longcalm at its own B (hl_lcbf) | 0.1142 | 0.1056 |
| af_std | 0.1205 | 0.1067 |
| b5sd1 | 0.1186 | 0.1067-0.1076 |
| g_s005_hie | 0.1214 | 0.1069 |
| longcalm at B 10 | 0.1156 | 0.1071 |
| af_rel4 pool | 0.1180 | 0.1073 |
| longconv | 0.1201 | 0.1074 |
| b5 | 0.1212 | 0.1094-0.1099 |
| best `bulk_solvent_runme.com` maps (0-32 ns of w1; +500 bar) | 0.124-0.125 | 0.1077-0.1079 |
| w1 / n1 (2 ns, the defaults) | 0.1242 / 0.1250 | 0.1108 / 0.1095 |
| k1B | 0.1281 | 0.1129 |
| CA-only (4 maps) | 0.1386-0.1429 | 0.1162-0.1180 |
| flat bulk solvent | 0.1462 | 0.1305 |

- **The overnight queue found nothing better.** It water-built 76 distinct maps that had never been through water
  building (2026-11-05). None beats the front-runners. Every `bulk_solvent_runme.com` map lands at 0.1077-0.1136,
  whatever its water model, pressure, averaging window or uniform B.
- **The pure 2019 MD map is slightly better than the example map made from it.** Squish editing and rescaling did not
  help (0.1051 vs 0.1062 before building, 0.1030 vs 0.1035 after).
- **The bespoke-workflow maps agree with each other, not with the 2019 map.** At 1.9-0.96 A their complex CC with
  longcalm is 0.84-0.95, but with the 2019 map only 0.72-0.77. The `bulk_solvent_runme.com` maps reach 0.59-0.68 with
  the 2019 map.

### 2. Water building erases most of the differences between maps

- **Gaps shrink by about half or more:**
  - The bespoke maps led ours by up to 0.009 before building and by ~0.003 after.
  - CA-only trailed by 0.02 before and 0.01 after.
  - On 7KR0 the MD map led a flat mask by 0.055 on the dry model. After rebuilding the waters the two finished within
    0.002, with mixed signs (finding 8).
- **Rankings before building mislead in both directions.** Stiff protein restraints, early 2-ns windows and OPC all
  looked better before building and not after. Late windows looked worse before building and better after.
- **What the 2019 map supplies is first-shell water.**
  - With our maps, about half of the positive Fo-Fc density sat on the 2019 run's core restrained-water sites: waters
    H-bonded to Glu32 OE1, the Arg62 backbone, the Val1 N-terminus, Tyr5 OH and Lys2 NZ.
  - With the 2019 map those sites were nearly clean, and adding waters gained it nothing.
  - No protein restraint weight or pressure made those waters appear in a `bulk_solvent_runme.com` map.
- **The `altconf` add_waters rule matters.** It counts only single-conformer atoms in the 2.3 A exclusion and gives
  each new water a conformer letter no neighbour has. It beats the old rules by 0.004-0.005 on every one of our maps.

### 3. Restraints: holding the protein on the crystal positions matters, high k does not

What each run restrains in its map window, per cell of 4 copies (2020 protein heavy atoms, 256 CA, ~700 waters;
`water_models/hurry/restraint_census.py`). Hydrogens are never restrained. Weights are density-derived, 0.1-10
kcal/mol/A^2; at 287 K, k = 10 holds an atom to 0.29 A rms and k = 1 to 0.92 A.

| run | protein heavy atoms (at k >= 1) | CA | water O (at k >= 1) | r3 |
|---|---|---|---|---|
| 2019 | all 2020, k = 10 | 256 | 192 (48 per ASU), k = 10 | 0.1030 |
| longcalm | 825 (528) | 102 | 177 (129) | 0.1052-0.1071 |
| af_std | 645 (427) | 56 | 205 (147) | 0.1067 |
| g_s005_hie | 612 (383) | 56 | 240 (139) | 0.1069 |
| b5sd1 | 715 (300) | 63 | 61 (50) | 0.1067-0.1076 |
| longconv | 180 (126) | 8 | 54 (42) | 0.1074 |
| b5 | 382 (152) | 15 | 37 (24) | 0.1094-0.1099 |
| k1B | all 2020, k = 1, never re-picked | 256 | 484, k = 1 | 0.1125-0.1129 |
| CA-only | none but CA | 232, k ~ 0.02 | none | 0.1162-0.1180 |

- **Too little restraint fails.** CA-only with k decayed to 0.02 makes the worst maps tested. The longer run of the
  same protocol (iterations 651-750) is no better.
  - The cause is the protein, not the lattice. After a translation fit, its high-density atoms sit 0.92-1.14 A rms
    from their reference positions, against 0.27-0.47 A in the restrained runs.
  - The two runs that stray most (loca 1.07 A, ff19_opc 1.14 A) make the two worst maps.
  - The solvent organises around a protein surface that is not the crystal's. Its high-resolution CC with longcalm is
    0.22-0.32 (restrained runs 0.76-0.92), and its first shell is thinner and flatter.
- **High k gains nothing once water is built.**
  - Tested with `bulk_solvent_runme.com` maps: k = 10 x rho (mean 51), 3 x rho, and flat 10 on every heavy atom.
  - Before building, the stiff maps were 0.002-0.003 worse (0.1263-0.1268 vs w1 0.1242).
  - After building they finished 0.1090-0.1106, against w1 0.1108 and n1 0.1095.
  - A 100-fold range of protein weight did not make the 2019 first-shell waters appear.
  - No 2026 `optimize_weights` run raised max_weight above 9.9999.
- **Uniform restraint does not work either.** k1B holds the protein closest of any run (0.27 A) and is still the worst
  restrained map, behind even the 2-ns w1. Protein rmsd alone does not predict the map.
  - k1B changes three things at once: k = 1 on everything, weights never re-picked, and 484 frozen water restraints.
    So which of them hurts is not known.
- **Within a build, more density-weighted restraint gave the better map.**
  - b5sd1 beats b5 at every stage (~0.002 at r3). They differ only in weight_scaledown.
  - longcalm beats longconv before building (0.1156 vs 0.1201); they tie after.
  - In both pairs the run without weight decay made the better map. Every front-runner, and the 2019 run at fixed
    k = 10, holds its restraints rather than decaying them.
- **Water restraints: suggestive, not proven.** Every front-runner restrains water O at density peaks, and the 2019 run
  restrained 48 crystallographic waters per ASU at k = 10. But b5sd1, with only 15 per ASU, ties af_std, which has
  51 per ASU.

### 4. The "0.4 A per iteration drift" was thermal jitter, and frame alignment changes nothing

- **It is jitter, not drift.** `alignment_drift_vs_itr.txt` records the last frame's translation. Over all 25 frames
  of an iteration the translation is already at full size 20 ps in, flat to 500 ps, and has mean zero.
- **Its size follows from the restraints.** It is the whole lattice rattling in the well of all its restraints
  together: rms = sqrt(3kT / (2 sum k)).
  - Predicted 0.43 A for 232 CA at k = 0.02; measured 0.42 A.
  - b5sd1 0.034 A, af_std 0.022 A.
  - Re-centring between sub-runs cannot remove it.
- **Replaying the pipeline's own alignment changes no R-free beyond noise.**
  - `extract_frames.csh ... align` replays the pipeline's alignment (cpptraj `rmsd ... norotate` onto
    `align_ref.crd`, the mask from each `align_<itr>.log`), and checks every frame against the log (all 600 matched).
  - Before building, every map moved by 0.0007 or less; at r3, by 0.0011 or less.
  - The b5sd1 and af_std maps were unchanged (CC 1.000).
  - It sharpens the CA-only edge a little, closing 15-25% of the first-shell gap; the protein's own 1 A departure keeps
    the rest.

### 5. The second-order settings, after water building

Each of these is worth about 0.002 or less once water is built.

| setting | result after building |
|---|---|
| averaging length | Means of 32 ns or more give 0.1077-0.1096 (mean 0.1086), against 0.1107 for 16 ns or less and for 2-ns windows. longcalm: 500 ns vs its last 2 ns, 0.1071 vs 0.1081 (before building 0.003) |
| early vs late frames | Before building the solvent drifts in 30-50 ns to a state that fits worse: ions move most, at charged side chains. After building late windows are as good or better (50-100 ns 0.1089 vs 0-100 ns 0.1082) |
| water model | 8 models within ~0.002 at 10 ns. SPC/E is best or tied in both its runs (0.1089, 0.1084). OPC's 0.001 lead before building does not survive. fb3mod vs OPC in the bespoke protocol ties (0.1100 vs 0.1094) |
| pressure | Positive targets (+500 to +3000 bar) average 0.1085 vs 0.1097 for w1, within noise. Only the under-filled map (698 waters, true -956 bar) is worse (0.1136-0.1142). True 1 bar is ~730 waters in this cell |
| uniform render B | B 2-30 within noise. B 10 was best before building. longcalm at B 5 (0.1052) is the best non-2019 map, but neither b5 run nor af_std gained from B 5: not replicated |
| B from MD motion (`rmsd2B`) | 0.002-0.010 worse after building (0.013-0.016 before): it blurs motion the frame average already holds |
| run's own per-atom B | longcalm and k1B gain ~0.002 before building, nothing after. af_std's own B is scrambled by the old OPC extra-point bug |
| salt | no trend: 0.05 M 0.1069, 0.2 M 0.1093, 0.68 M 0.1094 (one run each) |
| altloc_copies | ~0.002 cost with altloc_copies=0 on the deposited model; not significant on the example |

### 6. Pressure readings in Amber

`PRESSURE.md` on bespoke_amber_restraints branch `pressure` has the full account.

- **`ntp=4` reads ~500 bar high.** Every production run that leap2amber writes uses `ntp=4`, and it overcounts by about
  N_water kT / V. `ntp=1` and `ntp=2` (as the 2019 runs used) and the Monte Carlo barostat agree with each other. The
  `pressure` branch hydrates with `ntp=1`.
- **Restraint forces are not in pmemd's virial.** The same coordinates give the same reading with restraints on or off.
- **The restraints do change the true pressure, by about -560 bar.** The force field's preferred geometry pushes harder
  once released: the reading rises ~650 bar within 0.1 ps.
- **So the 2019 run's "+450 bar" came from its restraints, not from over-packing.** Unrestrained, that system reads
  about -400 bar.
- **The default `bulk_solvent_runme.com` runs were under-filled.** At a true ~-600 bar they held 709 waters. R-free
  improves up to a true 0 bar and is flat beyond it.

### 7. How the 2019 map differs from ours

- **What it is:** essentially the solvent-only average of one 100-ns 2019 run, `1aho_test/j5w8760/RestAll.nc`.
  Re-running the original 2019 scripts today gives a byte-identical map.
- **Its unusual features:** a 48-copy supercell, k = 10 on every protein heavy atom, k = 10 on 2298 crystallographic
  water O, 100 ns at constant volume, ff14SB + OPC.
- **Tested here, each worth ~0.002 or less:** length, water model, pressure, and protein restraint weight.
- **Never tested in isolation:**
  - water restraints at crystallographic sites at k = 10;
  - the 48-copy supercell (every 2026 run is 4 copies);
  - all heavy atoms held at fixed k = 10 for 100 ns.

### 8. 7KR0: the MD lead on a dry model was unbuilt water

- **Setup:** PRISM's six-conformer ensemble (0.77 A), refined with `bulk_solvent_runme.com` maps.
- **With PRISM's waters:** the MD maps are within 0.002-0.003 of flat.
- **Dried, then rebuilt conformer-aware:**
  - Dry, MD led flat by 0.055.
  - The lead fell under 0.002 by round 6.
  - The final R-free of MD and flat are within ±0.002, with mixed signs.
- **Scaling bistability:** with ANISOU and an MD map, phenix's k_iso/k_aniso fit has two solutions 0.002 apart.
- **Fixed along the way:** leap2amber lost the monoclinic beta angle (fixed on branch pressure, 3f076d4).

### 9. Practical lessons

- **Each `amber_N.nc` needs its own topology.** The barometer adds and drops waters, so each nc is read with a copy of
  `padded.parm7` stripped to its atom count. `padded.parm7` itself is never touched.
- **Stripped topologies can lose CMAP.**
  - cpptraj before V7.00.0 drops CMAP when stripping a topology. longcalm, af_std and g_s005_hie ran without the
    ff19SB backbone CMAP for that reason; the 2s2r runs have it.
  - `cpptraj_strip` (`$PROGRAMS/cpptraj/7.12.6`) avoids it.
- **`CCP4_SCR` must stay short (`/tmp/jamesh`).** With a longer path, pick.com found different peaks and add_waters.com
  placed waters in a different order.
- **Free sets differ between the two test datasets.** The example's `refme.mtz` and the deposited-model
  `refme_small.mtz` carry different free sets (overlap at chance level). Compare R-free within one protocol only.
- **The `/petabyte/home/programs` phenix is a stale copy.** It is the Aug 12 install; `/data/home/programs` has the
  Aug 26 install. octamus1, voltron and bl831 all mount `/data`.

## What is not known

- **Which restraint class matters.** Every restraint comparison above crosses builds:
  - af_std vs b5sd1 differ in CMAP, His54 state, the Cys12 split and protocol, as well as restraints;
  - k1B differs in three settings at once.
  - Only two within-build pairs exist (b5/b5sd1, longcalm/longconv), and each is one MD run.
- **Whether water restraints are what the 2019 map has and we lack.** b5sd1 argues against a simple count effect.
- **Whether the supercell matters.** It has never been tried in 2026.
- **Whether longcalm's B 5 lead is real.** It needs a replicate.

## Suggested next experiment

The bespoke project's goal is the minimal restraint set. We have bracketed it: CA-only at k 0.02 is too little, and
density-weighted all-heavy-atom restraints work. Settling which part is needed takes arms from one build (2s2r,
b5sd1's settings, no decay) that differ only in what is restrained:

1. density-weighted protein heavy atoms, no water restraints;
2. the same plus density-picked water restraints;
3. backbone only, plus water;
4. CA only held at k = 1-10 (no decay), with and without water.

Each arm is ~200 iterations of GPU MD (af_std took 15 h), then the same map and water-building test. A 2019-style
arm (all heavy atoms + crystallographic waters at fixed k = 10) would show whether the 2019 recipe alone still wins
in today's pipeline.

## Where things are

- `MD_BULK_SOLVENT_STRATEGY.md`: every test, in date order, with methods and the full tables.
- `water_models/hurry/`: frame extraction (`extract_frames.csh`, with the `align` option), map rendering
  (`render_frames.py`, `make_map.csh`, `bfac_map.csh`), map sets (`sets.txt`), chains (`*_chain.sh`, `alchain.sh`),
  `rfree_table.sh`, `restraint_census.py`, `cc_shells.py`.
- `water_models/addwater/`: the add_waters.com copy and `rounds_close.csh`.
- `water_models/waterbuild/`: the 76-map overnight queue and `results.txt` (mirrored in `/petabyte/jamesh/waterbuild/`).
- `water_models/example_check/`: refinement scripts and maps (`md2019_map.mtz` = the pure 2019 map).
- `test_7kr0/`: the 7KR0 ensemble test.
- `programs.sourceme` (tcsh) and `programs.sh` (bash): the program tree the scripts use.
- MD runs: `/scratch/jamesh/projects/amber/1aho_hurry_test/<run>/` (octamus1); run lineage in
  `/home2/jamesh/projects/amber/6c2r_AMPPNP/claude/solvent_candidates.txt`.
