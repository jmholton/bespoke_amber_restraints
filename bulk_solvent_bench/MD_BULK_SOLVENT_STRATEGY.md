# MD bulk-solvent strategy

How to make an MD bulk-solvent map (FPART/PHIPART for phenix's `refinement.input.bulk_solvent_map`) that refines as
well as the phenix_user_solvent example map, tested on 1aho. It starts with where the example map came from, then
records each test in date order. Scripts and runs are under `water_models/` (here).

## Lineage of the phenix_user_solvent example bulk-solvent map

`~/projects/phenix_user_solvent/example/solvent_Fpart.mtz` (Fpart/PHIpart), traced 2026-10-31.
Every link below was checked numerically, not just by file name or date. The scripts used for the checks are in
`water_models/lineage/` (here).

### Summary

- The example map is, shell for shell, essentially the solvent-only average of one 100 ns AMBER run from July
  2019: `~/projects/amber/1aho_test/j5w8760/RestAll.nc`. The complex correlation with that 2019 map is 0.96-0.999
  in every resolution shell. Everything that happened afterwards was rescaling plus one round of squish density
  editing.
- The 2019 map (`j5w8760/solvent_RestAll.map`, 1 Aug 2019) used all 500 frames (100 ns) of `RestAll.nc`. Re-running
  the original 2019 scripts on that trajectory today gives a **byte-identical** map.
- The pure 2019 MD map, before any squish editing, refines slightly *better* than the example map: R-free 0.1051 vs
  0.1062 with the example model and water occupancies refined. It is also better than any map from
  `bulk_solvent_runme.com` in both test protocols (see the end of this document).

### The chain, newest first

| # | file | date | relation to the step above it | check |
|---|---|---|---|---|
| 0 | `phenix_user_solvent/example/solvent_Fpart.mtz` | 2026-04-22 | (the example) | |
| 1 | `~/projects/amber.old/1aho_refine/best_Rfree/phenix_usersolvent/solvent_Fpart.mtz` | 2026-04-22 17:09 | byte-identical copy; also identical to the Fpart columns of `example/refme.mtz` | md5 |
| 2 | `best_Rfree/longrun/squish25/refme_scaled_450.mtz`, column Fpart (`_500` is the same) | 2026-03-24 | the example = this column (= 1.0078 x Fpart of `squish25/refme250.mtz`) | max diff 0.005 e, identical MTZ history |
| 3 | `squish25/refme1.mtz` (Fpart = squish25 start, Fpart2 = squish increment 1) | 2026-03-04 | example = 1.0151 x start + 1.016 x exp(-14.5 s^2/4) x increment 1 | residual rms 0.0005 of 21.2 |
| 4 | `best_Rfree/longrun/refme.mtz` | 2025-12-23 (written by `longrun/notes.com`), unchanged to 2026-03-03 | squish25 start = 0.0665 x exp(-9.75 s^2/4) x Fpart of `best_Rfree/thisone_small.mtz` (`cad ... scale file 2 0.0666, 9.5757`); untouched through 1280 refmac cycles | residual 0.1% |
| 5 | `best_Rfree/thisone_small.mtz` | 2025-09-01 | `superxform.com` reindex of `thisone.mtz` -> `1aho_test/rerefine_expanded8/post_watershed2/refme.mtz` -> `../watershed2/refme.mtz` (2x2x3 P1 supercell) back to the P212121 cell (`best_Rfree/notes.com` lines 4-14) | MTZ history |
| 6 | `1aho_test/rerefine_expanded8/watershed2/refme.mtz` | 2019-08-10 12:16 | `rerefine_expanded8/setup_mtz_notes.com`: `solvent.map` -> `../j5w8760/solvent_RestAll.map` extended to the supercell, `sfall`, FC/PHIC -> Fpart/PHIpart, `cad` with `1aho.mtz` | Fpart = 12 x FT(`solvent_RestAll.map`), CC 1.00000, B 0 |
| 7 | `1aho_test/j5w8760/solvent_RestAll.map` | 2019-08-01 20:34 | solvent average of all 500 frames of `RestAll.nc` (below) | re-made byte-identical (md5 ed39b1ea7b0636576ff635a6370634c0) |

A quirk at step 2: from cycle 250 on, squish25's refmac runs used a **second** partial structure, Fpart2 (squish
increment 2, scale 0.7365 at cycle 450). The example kept only the Fpart column, so it leaves out increment 2, and it
is not quite the solvent model those refmac runs used. Increment 2 was merged into Fpart at cycle 500 (`refme_500.mtz`
and later, CC 0.9993 with the example).

### How the 1 Aug 2019 map was made

#### The MD run: `1aho_test/j5w8760`, stage RestAll (July 2019)

- **Set up by `j5w8760/setup.com`:**
  - Protein: from `1aho_test/rerefine_expanded5/refmacout_minRfree.pdb` (2019-07-11), a refmac refinement of 48
    independent copies in the P1 2x2x3 supercell (91.8 x 81.4 x 90.3 A). Its waters with occupancy >= 0.1 were kept:
    2299 of them, ~48 per copy, each copy's waters refined separately. They became `superwater.pdb`.
  - Force field: Amber18 ff14SB (`parm10.dat` + `frcmod.ff14SB`) with **OPC** water.
  - Salt: 84 acetate (ACT) and 36 ammonium (AMM).
  - Water: 8760 in total (the "w8760" in the directory name); 81,984 atoms in all.
- **Stages:** Min, Cool, Heat, then RestHeat (0.2 ns), then RestAll, from `j5w8760/runme_Rest.com`.
- **RestAll:** 50,000,000 x 2 fs = 100 ns at 287 K (Langevin, gamma 2), cut 9 A, SHAKE.
  - `ntb=2 ntp=2` but `taup=99999999999`, so the volume is effectively constant; the pressure ran at +450 to +800 bar.
  - A frame every 100,000 steps (200 ps), so `RestAll.nc` has 500 frames.
- **Restraints:** 10 kcal/mol/A^2 to the start positions (`ref.crd` -> `start.crd`) on two groups:
  - protein: GROUP input SG, `CA * M`, `* * M`, `* * S`, `* * 3`, `* * B`, `O * E`, i.e. essentially every heavy atom;
  - **water:** the O of every water whose coordinates matched a crystallographic water in `superwater.pdb`
    (`boundwater_range.txt`: residues 3193-5490, 2298 waters).

#### The map: `1aho_test/allmapme.com RestAll.nc` (written 2019-06-06, run 2019-08-01)

`allmapme_RestAll.log` records "md2map solvent on 500 pdb files", i.e. the whole trajectory.

1. **cpptraj:** `parm start.top`, `reference start.crd`, `trajin RestAll.nc`, then
   `rms reference '!:ACT,AMM,WAT&@CA' ... norotate`, a translation-only fit of each frame on the protein CA atoms.
   Then `outtraj trajectory/md.pdb pdb multi pdbv3 keepext sg "P 1"`. The OPC extra points are not written to the
   PDBs.
2. **`1aho_test/md2map_multi.com P212121 pdbfiles.txt -solventonly`:** splits the frames over CPUs, each running
   `1aho_test/md2map.com` with B = 5 (the default; `allmapme.com` does not override it). For every frame:
   - keep the solvent residues only (HOH/WAT, AMM, ACT), H atoms included;
   - write each atom by element with occupancy 1 and B = 5;
   - wrap the coordinates into the 45.9 x 40.7 x 30.1 A unit cell (supercell / `MD_mult.txt` 2 2 3);
   - run `sfall mode atmmap`, `SYMM P212121`, `FORMFAC NGAUSS 5`, `GRID 192 144 256`.
3. **Average:** sum the per-frame maps and scale by 1/(500 frames x 12 cells x 4 ASUs) = 1/24000. That is `avg.map`,
   renamed `solvent_RestAll.map`: mean 0.13581 e/A^3 (41% of bulk water), max 7.649, rms 0.253.

**Reproduced 2026-10-31.** `water_models/lineage/md2map_rerun/rerun.csh` re-ran this recipe: today's cpptraj (V6.18)
with copies of the 2019 `md2map_multi.com` and `md2map.com`. The output `avg.map` has the same md5 as
`solvent_RestAll.map`.

## Free flags

Every file in the chain carries the same FreeR_flag column, generated by FREERFLAG on 2019-03-08 in `1aho_test`.
The deposited data have no test set: the FREE column of `1aho.mtz` is all 1. The `1aho.mtz` used for the
`bulk_solvent_runme.com` tests (`phenix_check/refme_small.mtz`) carries a different, later FREERFLAG set. It matches
the example's free set only at the random-chance level (83 of 1767).

## The restrained waters, and refined occupancies at their sites

Folding the 2299 restrained supercell waters into the unit cell (`water_models/lineage/waters/match_waters.py`):

- 81% lie within 1 A of a water in the deposited model (129 waters). 45 deposited waters have at least 24 of their 48
  copies restrained nearby ("core" sites), and 87 have at least one.
- Only 42% lie within 1 A of a water in the example model `starthere.pdb`. Its 128 water atoms are partial
  conformers with an occupancy sum of 21.5, and only 5 of them sit at core sites. The example model is nearly dry,
  so the bulk-solvent map has to supply most of the ordered water.

R-free with water occupancies refined (3 macro cycles, xyz + ADP + water occupancy; patched phenix 2.2-6143 on
octamus1; mean of 2 jiggled starts, 4 for the deposited model):

| model / data | flat | example map | 2019 MD map (pure) | our opc 10 ns | our spce 10 ns |
|---|---|---|---|---|---|
| `starthere.pdb` / `refme.mtz` | 0.1462 | 0.1062 | **0.1051** | 0.1265 | 0.1256 |
| `starthere.pdb` with the 129 deposited waters instead of its own / `refme.mtz` | 0.1324 | 0.1224 | 0.1192 | 0.1235 | 0.1237 |
| deposited 1aho (`withwater.pdb`) / `refme_small.mtz` | 0.1575 | - | **0.1479** | 0.1527 | 0.1542 |

Without occupancy refinement (example protocol, xyz + ADP, 4 starts), `starthere.pdb` gives flat 0.1470, example
0.1105, and our maps 0.1262-0.1284.

Mean refined occupancy at the 45 core sites vs the 84 other deposited-water sites (deposited model):
flat 0.82 / 0.71, 2019 MD map 0.58 / 0.48, our maps 0.66 / 0.53. The 2019 map lowers both groups. It holds more
solvent structure everywhere, not only at the sites its MD restrained.

## The "+450 bar" was the restraints

After RestAll, the 2019 run went through stages with the same atoms, volume and temperature but different restraints:

| stage | restrained | mean reported pressure |
|---|---|---|
| RestAll | protein heavy atoms + 2298 water O | +451 bar |
| RestPro1-4, RestMain, RestCA | protein atoms only | -170 to -280 bar |
| RestS | SG only | -374 bar |
| RestAllUnrst | nothing | -398 bar |

Unrestrained, that system reads about -400 bar. Its solvent was full but not over-packed: about 730 waters + 10 ions
per cell, about 747 water-equivalents, which is what our hydrated cells hold too.

Our own system (spce, end of 10 ns), 1 ns at the same volume (`water_models/pressure_test/`):

| setting | restrained | unrestrained |
|---|---|---|
| `ntp=1`, taup 1e11 | -448 | +70 |
| `ntp=4` (as in production) | +38 | +462 |

Restraints shift the reported pressure by hundreds of bar in both systems, in opposite directions. Even the absolute
value for the same unrestrained state depends on the pressure-control code path. So the restrained pressure that
`bulk_solvent_runme.com`'s hydration loop steers on is not the solvent pressure. (2026-11-03, below: the `ntp=4` path
reads ~500 bar high here; `ntp=1` and `ntp=2`, which the 2019 runs used, read true. The restraint shift itself is
real, not bookkeeping: pmemd leaves restraint forces out of the virial.)

## Difference maps, and filling their peaks with waters

`water_models/diffmaps/diffblobs.py` finds mFo-DFc flood-fill blobs above 3 rms, ranked by integrated density.

- **Example model (nearly dry):**
  - With our maps, about half of the positive difference density sits on the 2019 core water sites: 19-20 blobs,
    score 16 of 36. They are first-shell waters H-bonded to Glu32 OE1, the Arg62 backbone, the Val1 N-terminus,
    Tyr5 OH and Lys2 NZ.
  - With the 2019 or example map, those sites are nearly clean: 5-6 blobs, score 1.4-2.5.
  - Difference-map rms in the solvent region: flat 0.214, ours 0.157, 2019 0.117, example 0.112.
- **Deposited model (explicit waters):** the top features are the same for every map: Gly17 O, Arg62 NH2, Asp8 OD1,
  Arg56 CZ and HOH 99. These are atomic-model errors, not bulk solvent.

Filling the peaks with a copy of `~/Develop/add_waters.com` (`water_models/addwater/`; the copy adds `conf=` and
`heavy`). Example model, water occupancies refined, R-free as the mean of 2 starts:

| map | start | round 1 (4 sigma, occ 0.5, B 15) | round 2 | round 3 (distance from heavy atoms only, 2.3 A) |
|---|---|---|---|---|
| flat | 0.1462 | 0.1363 | 0.1310 | 0.1305 |
| example | 0.1062 | 0.1041 | 0.1046 | 0.1041 |
| 2019 MD | 0.1051 | 0.1057 | 0.1049 | 0.1049 |
| our opc 10 ns | 0.1265 | 0.1190 | 0.1191 | 0.1142 |
| our spce 10 ns | 0.1256 | 0.1164 | 0.1174 | 0.1143 |

Notes on the filling:
- **Round 3 option:** with the original 1.8 A cutoff measured from riding H atoms too, every water H-bonded to an
  N-H donor is skipped. That is why the `heavy` option was added.
- **Our maps:** filling closes a little over half the gap to the 2019 map. It takes about 60 waters (17
  occupancy-equivalents), and leaves a solvent-region rms of 0.135.
- **The 2019 map:** it gains nothing from added waters.

## Restraint weights and target pressure (scan, 2026-10-31)

`water_models/scan/`: opc, 2 x 1 ns, `bulk_solvent_runme.com` working tree (altloc_copies=1, new `pressure_target`)
frozen in `water_models/repo_scan/`. Example = starthere.pdb + refme.mtz, water occupancies, 2 jiggled starts;
deposited = withwater.pdb + refme_small.mtz, 4 starts; both bl831. "Unrestrained P" = 0.5 ns from the end of
production with every restraint off (`barometer.csh`), mean of the last 0.4 ns. Core score = positive Fo-Fc blob
score on the 2019 core water sites / all positive blobs (example model, `diffblobs.py`).

| run | weights (kcal/mol/A^2) | waters | restrained P | unrestrained P | example Rfree | deposited Rfree | core / total |
|---|---|---|---|---|---|---|---|
| w01 | 0.1 x rho, mean 0.5 | 693 | -228 | 201 | 0.1256 | 0.1539 | 16.7 / 38.1 |
| w03 | 0.3 x rho, mean 1.6 | 708 | 187 | 499 | 0.1244 | 0.1541 | 16.7 / 37.5 |
| w1 (default) | 1 x rho, mean 5.1 | 709 | -11 | 658 | 0.1241 | 0.1543 | 16.1 / 36.7 |
| w1_r2 (replicate) | same | 709 | -145 | 376 | 0.1241 | 0.1543 | 15.5 / 36.6 |
| w3 | 3 x rho, mean 15 | 719 | -59 | 734 | 0.1245 | 0.1545 | 15.7 / 38.0 |
| w10 | 10 x rho, mean 51 | 730 | -84 | 1165 | 0.1265 | 0.1544 | 15.8 / 38.0 |
| flat10 | 10 on every heavy atom | 719 | 44 | 776 | 0.1270 | 0.1537 | 16.9 / 39.7 |
| pm500 | w1, target -500 bar | 698 | -448 | 128 | 0.1263 | 0.1537 | 16.0 / 37.6 |
| pp500 | w1, target +500 bar | 719 | 297 | 768 | 0.1243 | 0.1535 | 15.8 / 37.4 |
| pp1500 | w1, target +1500 bar | 753 | 1438 | 1921 | 0.1254 | 0.1534 | 16.7 / 39.4 |
| ac0, ac0_r2 | w1, altloc_copies=0 | 712, 707 | 31, -112 | 513, 355 | 0.1260, 0.1241 | 0.1532, 0.1522 | 15.5-16.1 / 37.7 |
| opc_r2 (old script, = altloc_copies 0) | w1 | 706 | -126 | | 0.1272 | 0.1519 | 16.1 / 36.9 |
| 2019 map | | | | | 0.1051 | 0.1479 | ~1 / ~10 |

- A 100-fold range of protein restraint weight and target pressures from -500 to +1500 bar move the deposited-model
  R-free only within 0.1534-0.1545 and the example-model R-free within 0.1241-0.1270. The ~16 units of Fo-Fc density
  on the 2019 core water sites are there in every run: holding the protein harder or softer does not make those
  first-shell waters. The 2019 map's advantage is not its protein restraint weight or its pressure.
- The unrestrained pressure depends on the water count alone: P = 30.4 bar/water x (N - 692), rms residual 103 bar
  over 12 runs, whatever the restraint weight (kappa_T predicts ~34 bar/water). The restraint weight only moves the
  restrained reading, by +300 to +1250 bar (more for stronger restraints), and so where hydration stops: 709 waters at
  the default weights is ~+500 bar for the liquid; 0 bar would be ~692 waters. (All `ntp=4` readings: the true
  values are ~500 bar lower, and the restraint shift is real; see 2026-11-03.)
- altloc_copies: w1 and w1_r2 reproduce each other to 0.0000 (example) and 0.0000 (deposited). With altloc_copies=0
  (ac0, ac0_r2, opc_r2) the deposited protocol gives 0.1524 +/- 0.0007, vs 0.1543 with altloc_copies=1: a ~0.002 cost.
  The example protocol goes the other way, 0.1258 +/- 0.0016 vs 0.1241, not significant.
- Occupancy groups: starthere.pdb is `example/under20_defaultsolvent_xyzoB_001.pdb` (phenix 2.1rc2, 10 cycles, flat
  solvent, xyz+ADP+occupancies, `under20_occ_groups.eff` = one constrained group per residue, 64 groups, 623
  residue-conformers, + individual=water). phenix's automatic groups in our runs are the same 623 residue-conformer
  groups; re-running with the explicit .eff gives 0.1061 vs 0.1060 (example map) and 0.1247 vs 0.1239 (w1 map).

## Waters next to split side chains, pressure, and 100 ns (2026-11-01)

Example model only from here on. The add_waters copy has an `altconf` option: the 2.3 A exclusion counts only
single-conformer heavy atoms, a peak within `mindist` (1.0 A) of any atom is dropped, and each new water starts as
conformer A and goes through `waterconf.awk` (unmodified, clashdist 2.3 A) with symmetry-aware DIST lines from
`addwater/close_waters.py`, so it gets a letter no atom within 2.3 A has. `addwater/rounds_close.csh`: 3 rounds, all
`occ=0.5 B=15 heavy dist=2.3 altconf`, each filling the j1 refinement of the round before. R-free, mean of 2 starts
(j1/j2 agree within 0.0006):

| map | no added waters | old rules, round 3 | altconf round 1 | round 2 | round 3 |
|---|---|---|---|---|---|
| w1 | 0.1241 | 0.1150 | 0.1139 | 0.1111 | 0.1096 |
| w1_r2 | 0.1241 | 0.1149 | 0.1142 | 0.1101 | 0.1090 |
| pm500 | 0.1263 | 0.1164 | 0.1149 | 0.1129 | 0.1142 |
| pp500 | 0.1243 | 0.1149 | 0.1130 | 0.1089 | 0.1082 |
| pp1500 | 0.1254 | 0.1122 | 0.1111 | 0.1074 | 0.1074 |
| pp1500_r2 | 0.1285 | 0.1147 | 0.1131 | 0.1083 | 0.1082 |
| pp3000 (789 waters) | 0.1274 | 0.1132 | 0.1124 | 0.1083 | 0.1080 |
| w1 0-100 ns | 0.1245 | | 0.1123 | 0.1100 | 0.1083 |
| w1 50-100 ns | 0.1262 | | 0.1139 | 0.1100 | 0.1095 |
| 2019 map | 0.1051 | 0.1049 (old rules, round 5) | 0.1042 | 0.1025 | 0.1038 |

- altconf beats the old rules by 0.004-0.005 on every one of our maps. The Arg62 cluster, the Gly61 peak and the
  Lys2/30/50/58 NZ peaks leave the difference map: positive blob score 17-22 -> 9-12 (2019 map 9.6 -> 7.9;
  `scan/diffblobs_after_close_waters.txt`). What remains sits on protein atoms (Glu24 OE1 F, Gly43 N, Lys28 NZ,
  Asp9) plus the negative HOH164 and Arg56 NH1 B, which added waters cannot fix.
- Pressure: before water building the ranking is noise (pp1500 0.1254 vs its replicate 0.1285). After altconf
  building, every positive target (+500 to +3000 bar) gives 0.1074-0.1083, the default w1 0.1090-0.1096, and -500 bar
  0.1129: more water helps by ~0.001 and saturates by +500 bar; fewer waters costs ~0.004.
- 100 ns (`scan/w1_100ns/`, `scan/windows_driver.csh`): before water building, every cumulative map from 0-2 to
  0-100 ns gives 0.1240-0.1248, but later 2-ns windows give 0.1256-0.1286 and 50-100 ns 0.1262. 10-ns block maps
  correlate less the farther apart they are, at every resolution (1.4-0.96 A: 0.973 adjacent, 0.961 at 20 ns, 0.956
  beyond 50 ns; `scan/drift/block_sep.txt`), levelling off after 30-50 ns: the solvent equilibrates in ~30-50 ns to
  a state that fits the crystal worse than its start. The 26 NH4/acetate move most (ion-only map CC 0-20 vs 80-100
  ns 0.30-0.74; water-only 0.97-0.99); at low resolution the ion and water changes cancel, and beyond 3.5 A water
  carries 85-90% of the change (`scan/drift/drift_split.txt`). The drift is largest at charged side chains (Arg18,
  then Asp9, Glu24, Arg56, Lys30), not at the unrestrained split side chains as such (33% of the solvent volume, 30%
  of the drift; `scan/drift/drift_early_late.txt`). After altconf water building the penalty is gone: 0-100 ns
  0.1083 and 50-100 ns 0.1095 vs 0.1090-0.1096 for 2 ns.

## What the pressure really is (2026-11-03)

`water_models/pressure_test/` (`decomp.csh`, `release.csh`, `calib.csh`, `mc_scan.csh` and their `.py` summaries),
pmemd.cuda_SPFP, Amber22, production settings (287 K, Langevin, cut 9).

- **Restraint forces are not in the virial.** The same coordinates and velocities give the same reported pressure
  to 0.1 bar with restraints on or off, through `ntp=4`, `ntp=1` and `ntp=2`. That holds on 22 frames (step 1) and on
  two frames held still for 20 steps with dt = 1e-5 ps.
- **`ntp=4`, which leap2amber writes for every production run, reads high.** In a pure OPC box (2542 waters) whose
  volume the Monte Carlo barostat (`barostat=2`, no virial used) set at 1 bar, the fixed-volume reading is:

  | path | reads | true, from the MC equation of state |
  |---|---|---|
  | `ntp=4` | +1405 | +103 |
  | `ntp=1`, Berendsen, taup 1e11 | +116 | +103 |
  | `ntp=2`, Berendsen, taup 1e11 | +105 | +103 |

  The error appears only once atoms move (the frozen frames agree). It is about N_water kT / V: 1329 bar predicted vs
  1302 measured in the box, and ~500 bar in the 1aho cell. Every pressure we reported from `ntp=4` (hydration
  targets, production, the barometer) is about 500 bar too high. The 2019 numbers above (`ntp=2`) are true readings.
- **The restraints change the true pressure, by ~-560 bar.** From the end of `w1_100ns`, the MC barostat (pmemd
  scales the restraint reference points with the box: restraint energy 313 kcal/mol at 1 and 1000 bar) puts the
  crystal's own volume at -817 +/- 44 bar restrained and -259 +/- 51 bar unrestrained. `ntp=4` shows the same gap
  (-72 restrained over 100 ns vs +376 to +658 in the barometer runs). Released, the reading rises by ~650 bar within
  0.1 ps, while the protein's 1-4 electrostatic energy falls by 90 kcal/mol: the restrained atoms are held off the
  force field's local geometry, and the free protein pushes harder on its solvent at once. Our maps come from the
  restrained MD, so the restrained true pressure is the one that describes their solvent.
- True restrained pressure of the scan runs (MC at 1 bar, 1 ns from prod_2, converted with the restrained cell's
  0.449 A^3/bar):

  | run | waters | `ntp=4` reading | true P (bar) | R-free after altconf round 3 |
  |---|---|---|---|---|
  | pm500 | 698 | -448 | -956 +/- 39 | 0.1142 (round 2 0.1129) |
  | w1 | 709 | -11 | -624 +/- 37 | 0.1096 |
  | pp500 | 719 | +297 | -466 +/- 63 | 0.1082 |
  | pp1500 | 753 | +1438 | +1115 +/- 39 | 0.1074 |
  | pp3000 | 789 | +3074 | +2576 +/- 45 | 0.1080 |

  True 1 bar is ~730 waters in this cell, `pressure_target` ~+550 on the `ntp=4` scale. The default runs (w1) were
  under-filled at about -600 bar. R-free improves up to about true 0 bar and is flat beyond it. By coincidence, the
  `ntp=4` overcount (+500) roughly cancels the restraint shift (-560), so w1's 709 waters is about the 0-bar count of
  the *unrestrained* cell.
- Releasing all restraints does not fix `ntp=4` (`unrestrained.csh`, `barometer_protocol.csh`, `protocol_summary.txt`).
  optimize_weights_runme.com's barometer protocol (5 ps at dt 0.5 fs, ntr=0) from 10 frames reads +478 through
  `ntp=4` vs -145 through `ntp=1` (+623 +/- 100), single readings +/-230 bar. Over 1 ns unrestrained, `ntp=4` gives
  +460 vs `ntp=1`/`ntp=2` -77 to +65 and MC -169 to -259. In water the `ntp=4` error is the same at dt 0.5, 1 and
  2 fs. Write-up: bespoke_amber_restraints branch `pressure`, `PRESSURE.md`.

## Hydrating to a true 0 bar, and water B factors (2026-11-03, example model)

**True-pressure hydration.** The `pressure` branch uses `ntb=2,ntp=1,barostat=1,taup=99999999999` in place of
`ntp=4`. It was run from the `water_models/repo_ntp1` snapshot via `scan/launch_ntp1.csh` and `scan/close_chain.csh`,
with otherwise the w1 settings and target 0. The first try died because the `*****` blow-up checks matched pmemd's
`taup =**********` header; the checks are fixed on the branch.

- **Hydration:** n1 and n1_r2 reached 723 and 729 waters (old default w1: 709). Production read -164 and +50 bar,
  now true pressures.
- **R-free:** mean of 2 starts after each round of altconf water building:

  | map | waters | none | round 1 | round 2 | round 3 |
  |---|---|---|---|---|---|
  | n1 | 723 | 0.1250 | 0.1138 | 0.1088 | 0.1099 |
  | n1_r2 | 729 | 0.1266 | 0.1134 | 0.1109 | 0.1105 |
  | w1 / w1_r2 (`ntp=4`) | 709 | 0.1241 / 0.1241 | 0.1139 / 0.1142 | 0.1111 / 0.1101 | 0.1096 / 0.1090 |
  | pp500 (`ntp=4`) | 719 | 0.1243 | 0.1130 | 0.1089 | 0.1082 |

- **No gain:** filling to a true 0 bar did not beat the old default. pp500 and n1 hold the same water count yet
  differ by 0.0017 after round 3.
- **The maps themselves:** the n1 maps correlate with w1, pp500 and pp1500 exactly as the w1 replicates correlate
  with each other (complex CC 0.985-0.990 below 3.5 A, 0.93-0.95 at 3.5-1.9 A, 0.91-0.94 at 1.9-0.96 A).
- **Conclusion:** the R-free differences between 709 and 789 waters (0.1074-0.1105) are at the level of
  MD-replicate noise. Only the under-filled pm500 (698 waters, 0.113-0.114) stands out. The "+0.001 from more
  water" seen earlier needs more replicates before it can be believed.

**Water B factors** (`scan/bmap.csh`). The n1 solvent trajectory was re-rendered with nc2mtz_gemmi.com from the
bulk_solvent clone (it has the water extra-point B fix), then refined against the example model with no added waters,
2 starts each (agreeing within 0.0007):

| water B | R-free |
|---|---|
| uniform 2 | 0.1261 |
| uniform 5 | 0.1261 |
| uniform 10 (default, = n1) | 0.1250 |
| uniform 30 | 0.1281 |
| `rmsd2B` per atom, 2-100 | 0.1408 |
| `rmsd2B` B field by location (`Bfac_map=auto`, `selresn=HOH`), 2-100 | 0.1414 |
| `rmsd2B` B field, 10-50 | 0.1375 |
| `rmsd2B` B field, 2-10 (sharpen bound waters only) | 0.1278 |

- **Uniform B:** phenix's per-bin k_mask absorbs most of it, but not all, and 10 is the best of these values.
- **B from MD motion:** `rmsd2B` B describes motion that the frame average already contains. Over 2 ns, 17% of water O
  have B < 30 and the median is 255. Adding it as B blurs the mobile water a second time and costs 0.013-0.016.
  Per atom or by location makes no difference.
- **Technical notes:**
  - The frames name water `HOH` (cpptraj `pdbv3`), so the field needs `selresn=HOH`.
  - Two renders started at once on a fresh checkout race to compile `Bfac_map`.

## Long bespoke runs as solvent maps (2026-11-03, example model, octamus1)

Source: candidates from `/home2/jamesh/projects/amber/6c2r_AMPPNP/claude/solvent_candidates.txt`. These are
1aho_hurry_test runs of the bespoke workflow under `/scratch/jamesh/projects/amber/1aho_hurry_test/<run>/`: 0.5-ns
`amber_N.nc` iterations of 25 frames, on the crystal frame.

**Method** (`water_models/hurry/`):
- **Extraction (`extract_frames.csh`):** per iteration, a copy of `padded.parm7` stripped to that nc's atom count.
  Only WAT, NH4 and ACY are kept, extra points dropped.
- **Rendering (`render_frames.py`):** each frame is rendered as nc2mtz_gemmi.com renders it (occupancy 1, uniform
  B 10, gemmi sfcalc in P212121 to 0.965 A, rate 2), then the frames are averaged and converted to FPART/PHIPART.
  Re-rendering the n1 frames this way reproduces n1's map exactly (CC 1.0, max difference 2e-7).
- **Frame sets:** `sets.txt`. The 500 ns sets take every 10th iteration (5 frames each); the 200-iteration runs take
  their second halves, every 5th iteration (50 ns, 100 frames).
- **Refinement (`hurry_chain.csh`):** example refinement plus 3 altconf rounds, all on octamus1. The w1, n1 and 2019
  maps were re-run there as same-host references; their starting values match bl831 to 0.0001.

R-free, mean of 2 starts:

| map | frames | none | round 1 | round 2 | round 3 |
|---|---|---|---|---|---|
| longcalm (fb3mod, 500 ns) | 505 | 0.1156 | 0.1091 | 0.1068 | 0.1071 |
| longcalm, last 2 ns only | 100 | 0.1186 | 0.1115 | 0.1101 | 0.1081 |
| longconv (fb3mod, 500 ns) | 505 | 0.1201 | 0.1116 | 0.1081 | 0.1074 |
| af_rel4 group pooled (OPC, 5 x 50 ns) | 500 | 0.1180 | 0.1100 | 0.1074 | 0.1073 |
| af_std (OPC) | 100 | 0.1205 | 0.1089 | 0.1066 | 0.1067 |
| noalt (fb3mod) | 100 | 0.1249 | 0.1142 | 0.1092 | 0.1100 |
| noalt_opc (OPC) | 100 | 0.1214 | 0.1144 | 0.1104 | 0.1094 |
| g_s005_hie (0.05 M) | 100 | 0.1214 | 0.1121 | 0.1088 | 0.1069 |
| g_s068_hie (0.68 M) | 100 | 0.1234 | 0.1114 | 0.1091 | 0.1094 |
| salt02_opt (0.2 M) | 100 | 0.1239 | 0.1133 | 0.1109 | 0.1093 |
| w1 (bulk_solvent_runme, 2 ns) | 100 | 0.1242 | 0.1137 | 0.1109 | 0.1108 |
| n1 (bulk_solvent_runme, true 0 bar, 2 ns) | 100 | 0.1250 | 0.1142 | 0.1099 | 0.1095 |
| 2019 map | 500 | 0.1051 | 0.1043 | 0.1024 | 0.1030 |

- **Map character:** the hurry maps resemble each other far more than they resemble ours or the 2019 map. Complex CC at
  1.9-0.96 A: 0.92-0.95 with longcalm (longconv 0.88), against 0.76-0.77 for w1 and n1 and 0.72-0.77 for the 2019
  map. Their mean density is 0.129-0.131 e/A^3, against 0.136-0.139 for ours.
- **Flatness:** 1aho has almost no true bulk: only 3.3% of the cell lies 3.5-5 A from protein, and essentially none
  beyond 5 A. In that band the density sd is 0.16-0.21 e/A^3 for most maps and 0.114 for longconv. longcalm's 2 ns and
  500 ns maps have the same sd (0.160) and correlate at 0.95 even at high resolution.
- **Sampling:** going from longcalm's last 2 ns to 500 ns gains 0.003 R-free before water building.
- **The model:** the setup counts for more than the sampling. Even longcalm's 2-ns map beats our 2-ns maps by 0.006.
- **Water building:** it narrows the spread. The best rounds are:

  | map | best round, R-free |
  |---|---|
  | af_std | 0.1066 |
  | longcalm | 0.1068 |
  | g_s005_hie | 0.1069 |
  | af_rel4 pool | 0.1073 |
  | longconv | 0.1074 |
  | our n1 / w1 | 0.1095 / 0.1108 |
  | 2019 map | 0.1024 |

  The best hurry maps beat ours by ~0.003 after building, where they led by up to 0.009 before it. The 2019 map is
  still ~0.004 ahead of the best of them.
- **Sampling after building:** 500 ns beats 2 ns of longcalm by only 0.001 (0.1071 vs 0.1081).
- **Water model:** noalt (fb3mod) vs noalt_opc (OPC) is 0.1100 vs 0.1094, a tie. Salt shows no clear trend (0.05 M
  0.1069, 0.2 M 0.1093, 0.68 M 0.1094).
- **Replicability:** same-host re-runs of a map diverge by ~0.001 by round 2 (n1: 0.1099 here vs 0.1088 on bl831),
  because each round fills peaks from the round before.

## Lower B, and each run's own B, on the long-run maps (2026-11-04, example model, octamus1)

longcalm (500 ns) and af_std (50 ns) re-rendered from the same frames as above:
- **Uniform B 5 and 2:** `render_frames.py` (`sets.txt` lines with `B=`).
- **The run's own per-atom B:** `bfac_map.csh`.
  - Each `amber_N.nc` is rendered by `nc2mtz_gemmi.com` (my clone, with the 01502c8 extra-point fix) with
    `premod_Bfac_<N+1>.pdb` as the B sidecar. That is the Bfac.pdb the pipeline itself rendered amber_N with.
  - The pipeline's clip [1 : 999.99] is applied, and only the solvent is rendered (`bfac_iter.csh`).
  - Every iteration's sidecar was checked against the frame atom by atom (name and residue).
  - longcalm's last iteration (2000) has no sidecar and was skipped (100 of 101 iterations).
  - The B values given to the solvent, as 10/25/50/75/90th percentiles:

    | run | 10th | 25th | 50th | 75th | 90th | water atoms B <= 5 | water atoms B > 100 |
    |---|---|---|---|---|---|---|---|
    | longcalm | 2 | 2 | 3.3 | 8.4 | 37 | 66% | 2% |
    | af_std | 2 | 2.3 | 6.0 | 16 | 52 | 45% | 4% |

    Most waters are sharper than the uniform 10, and a tail is nearly blurred out.
  - af_std is an OPC run, so its B was refined while the old extra-point bug applied water B to the wrong atoms
    (`solvent_candidates.txt`).

Flatness (`bulk_flatness.py`) is the density sd/mean beyond 5 A, a handful of voxels:

| map | B 10 | B 5 | B 2 | own B | 2019 map |
|---|---|---|---|---|---|
| longcalm | 0.38 | 0.47 | 0.49 | 0.35 | 0.32 |
| af_std | 0.58 | 0.80 | 0.86 | 0.39 | 0.32 |

Lower uniform B makes the bulk spikier. Own B is as flat as or flatter than B 10, and closest to the 2019 map.

R-free, mean of two jiggled starts (the B 10 column is the run in the section above):

| map | B | none | r1 | r2 | r3 |
|---|---|---|---|---|---|
| longcalm | 10 | 0.1156 | 0.1091 | 0.1068 | 0.1071 |
| longcalm | 5 | 0.1157 | 0.1082 | 0.1066 | **0.1052** |
| longcalm | 2 | 0.1159 | 0.1115 | 0.1075 | 0.1069 |
| longcalm | own | **0.1142** | 0.1082 | 0.1070 | 0.1056 |
| af_std | 10 | 0.1205 | 0.1089 | 0.1066 | 0.1067 |
| af_std | 5 | 0.1205 | 0.1086 | 0.1070 | 0.1069 |
| af_std | 2 | 0.1205 | 0.1089 | 0.1084 | 0.1074 |
| af_std | own | 0.1207 | 0.1092 | 0.1072 | 0.1075 |
| 2019 map | 5 | 0.1051 | 0.1043 | 0.1024 | 0.1030 |

- **longcalm B 5** ends at 0.1052, 0.002 below B 10 and the best non-2019 map so far. It leads B 10 in all three
  water rounds, but only by 0.0009, 0.0002 and 0.0019, so it needs a replicate.
- **longcalm own B** is the only variant that gains before water building (0.0014), and it ends at 0.1056.
- **B 2** is worse for both runs once waters are built.
- **af_std** gains nothing from any B. That fits its own B being scrambled by the extra-point bug.
- The 2019 map (uniform B 5) is now ~0.003 ahead of the best hurry map, down from ~0.004.

## Three overnight ff19/OPC runs (2026-11-05, example model, octamus1)

New 1aho_hurry_test runs: ff19 + OPC, 200 iterations, ~730 waters each. Their settings differ only in these:

| run | weight_scaledown | max_mult | release/repick | water B in the run |
|---|---|---|---|---|
| ff19_opc_b5_2s2r | 0.9 | 2 | on | uniform 5 (render_B) |
| ff19_opc_b5sd1_2s2r | 1 | 2 | on | uniform 5 (render_B) |
| ff19_opc_k1B_2s2r | 1 | 1 | off | per-atom, refined (Bfac_maxmod 2) |

Maps are made from the same frames as af_std: iterations 101-200, every 5th, 5 frames each, 100 frames (`sets.txt`,
chain `hurry/newmaps_chain.sh`). Each run is rendered twice:
- at uniform B 10, like every map above;
- the way the run renders: uniform B 5 for the b5 runs, and k1B's own per-atom B via `bfac_map.csh`.

af_std's B 10 map was run again alongside, as a repeat.

R-free, mean of two jiggled starts:

| map | B | none | r1 | r2 | r3 |
|---|---|---|---|---|---|
| b5 | 10 | 0.1212 | 0.1149 | 0.1099 | 0.1094 |
| b5 | 5 | 0.1212 | 0.1149 | 0.1088 | 0.1099 |
| b5sd1 | 10 | 0.1186 | 0.1094 | 0.1067 | 0.1076 |
| b5sd1 | 5 | 0.1188 | 0.1104 | 0.1076 | 0.1074 |
| k1B | 10 | 0.1281 | 0.1168 | 0.1130 | 0.1129 |
| k1B | own | 0.1260 | 0.1149 | 0.1132 | 0.1125 |
| af_std, repeat | 10 | 0.1205 | 0.1089 | 0.1066 | 0.1067 |
| longcalm (above) | 10 / 5 / own | 0.1156 / 0.1157 / 0.1142 | 0.1091 / 0.1082 / 0.1082 | 0.1068 / 0.1066 / 0.1070 | 0.1071 / 0.1052 / 0.1056 |
| w1 / n1 (above) | 10 | 0.1242 / 0.1250 | 0.1137 / 0.1142 | 0.1109 / 0.1099 | 0.1108 / 0.1095 |
| 2019 map (above) | 5 | 0.1051 | 0.1043 | 0.1024 | 0.1030 |

- **Repeat:** the af_std run reproduced every one of its 16 R values from 2026-11-03 exactly, so on one host the test
  is deterministic.
  - The noise is in the water picks. Rendering the same frames at B 5 instead of B 10 moves a round by up to 0.001.
- **b5sd1 is the best of the three:**
  - Before building it reaches 0.1186, against af_std 0.1205 and the five-run af_rel4 pool's 0.1180.
  - After building it reaches 0.1067-0.1076, a tie with af_std and with longcalm at B 10.
  - It does not reach longcalm at B 5 or own B (0.1052/0.1056), or the 2019 map (0.1024-0.1030).
- **b5 vs b5sd1:** b5 differs from b5sd1 only in weight_scaledown (0.9 instead of 1). It trails b5sd1 at every stage, by
  0.001-0.005, and by ~0.002 at round 3. That is one MD run each.
- **k1B is the worst of the long runs:**
  - It is behind even our 2-ns maps: 0.1281 before building (w1 0.1242), and 0.1125-0.1129 after (w1 0.1108, n1 0.1095).
  - Its own B gains 0.002 before building, as longcalm's own B did, and nothing after.
- **B 5 vs B 10:** neither b5 run gains from B 5, and af_std did not either, so longcalm's B 5 lead is still not
  replicated.
- **R-work:** the b5 maps reach a lower R-work (0.086 at round 3, against 0.0885 for af_std) at the same R-free.
- **Map character:**
  - Mean density is 0.1344-0.1345 e/A^3, against 0.129-0.131 for the earlier long runs and 0.136-0.139 for ours.
  - Complex CC at 1.9-0.96 A with longcalm is 0.84 for b5, 0.87 for b5sd1 and 0.76 for k1B. With the 2019 map it is
    0.74, 0.76 and 0.59.
  - Density sd 3.5-5 A from the protein, at B 10: b5 0.136, b5sd1 0.146, k1B 0.300. For comparison, longcalm 0.160,
    af_std 0.199 and the 2019 map 0.173.
  - Smoothness does not track R-free. k1B is the spikiest and the worst, but b5 is smoother than b5sd1 and also worse.

## Which atoms the runs restrain, and the CA-only end (2026-11-05, example model, octamus1)

The bespoke project's goal is the minimal set of restraints. `hurry/restraint_census.py` counts what each run
restrained in its map window, from the ntr=1 groups of its `amber_<itr>.in` (per cell of 4 copies: 2020 protein heavy
atoms, 256 CA, ~700 waters):

| run | protein heavy atoms (at k >= 1) | CA | water O (at k >= 1) |
|---|---|---|---|
| longcalm | 825 (528) | 102 | 177 (129) |
| af_std | 645 (427) | 56 | 205 (147) |
| g_s005_hie | 612 (383) | 56 | 240 (139) |
| b5sd1 | 715 (300) | 63 | 61 (50) |
| b5 | 382 (152) | 15 | 37 (24) |
| longconv | 180 (126) | 8 | 54 (42) |
| k1B | all 2020, k = 1, never re-picked | 256 | 484, k = 1 |
| ff19_opc_mult, _loca, _locw | none but CA | 232, k ~ 0.02 | none |

- No run restrains hydrogens. The 2019 run restrained every protein heavy atom and 48 deposited water O per ASU, at
  k = 10.
- Weights are density-derived, 0.1-10 kcal/mol/A^2. At 287 K, k = 10 holds an atom to 0.29 A rms, k = 1 to 0.92 A.
- Within a build, more restraint gave the better map: b5sd1 over b5, longcalm over longconv (before building).

**The CA-only runs** come from `build_ff19_opc` (Nov 1: 11 NH4 / 15 acetate, ff19SB + OPC, CMAP intact):
- 232 CA restrained, leaving a 5 A hole around Cys12. k is ramped from ~1 to 0.02 by iteration 113 (halving every
  20 iterations) and then held. Nothing else is restrained.
- mult, loca and locw differ only in the run's own B handling, which a uniform-B render ignores. ff19_opc is their
  parent, the same protocol run to iteration 751.

R-free, mean of two jiggled starts, frames as for every map above (iterations 101-200, B 10):

| map | none | r1 | r2 | r3 |
|---|---|---|---|---|
| ff19_opc_mult | 0.1386 | 0.1218 | 0.1151 | 0.1162 |
| ff19_opc_loca | 0.1399 | 0.1215 | 0.1172 | 0.1179 |
| ff19_opc_locw | 0.1393 | 0.1219 | 0.1164 | 0.1174 |
| ff19_opc, iterations 651-750 | 0.1429 | 0.1229 | 0.1168 | 0.1180 |

- These are the worst maps tested: 0.02 behind the restrained runs before building and 0.01 behind after three
  rounds, and behind k1B (0.1281 / 0.1129) and our 2-ns w1 / n1 maps. The longer run is no better.
- **Map character:** complex CC with longcalm at 1.9-0.96 A is 0.22-0.32, against 0.76-0.92 for every restrained
  run. The four CA-only maps agree with each other no better (0.51-0.67). The first shell (2.5-3.5 A from the
  protein) is thinner and flatter: mean 0.22-0.23, sd 0.14-0.16, against 0.26-0.28 and 0.26-0.44.

### Frame alignment

The hurry maps render raw `amber_<itr>.nc` frames. Before it makes its own mtz, `optimize_weights_runme.com`
translates every frame onto its alignment reference: cpptraj `rmsd ... reference norotate` onto `align_ref.crd`,
over the protein atoms with reference density >= 1 (the mask and all 25 translation vectors are in each
`align_<itr>.log`).

- **What the "drift" is:** `alignment_drift_vs_itr.txt` records the last frame's translation. The per-frame rms
  translation is already at its full size 20 ps into an iteration and stays flat to 500 ps, with mean zero:

  | run | rms translation, frame 1 | frame 25 |
  |---|---|---|
  | ff19_opc_mult | 0.40 A | 0.42 A |
  | b5sd1 | 0.033 A | 0.035 A |
  | af_std | 0.022 A | 0.022 A |

  So it is thermal jitter of the whole lattice in the well made by all its restraints together, not drift. The rms is
  sqrt(3kT / (2 sum k)): 0.43 A predicted for 232 CA at k = 0.02, 0.42 A measured. Re-centring the MD between
  sub-runs would not remove it. For the restrained runs it is far too small for a map to see.
- **`extract_frames.csh ... align`** (an `align` token on the `sets.txt` line) reproduces the pipeline's step: the same
  cpptraj command on each iteration, with that iteration's mask from its log. Every frame's translation is checked
  against the log, and the extraction stops if any differs by more than 0.001 A. All 600 frames of the six aligned
  maps matched. Tags without the token render as before.
- **Effect on the maps:** for the CA-only runs, CC between the aligned and raw maps is 0.86-0.92 at 1.9-0.96 A and
  ~1.0 below 2.5 A. Their high-resolution CC with longcalm rises only from ~0.32 to ~0.34. b5sd1 and af_std are
  unchanged (CC 1.000 in every shell).
- **R-free of the aligned maps** (`hurry/alchain.sh`), mean of two starts, raw -> aligned:

  | map | none | r3 |
  |---|---|---|
  | ff19_opc_mult | 0.1386 -> 0.1387 | 0.1162 -> 0.1155 |
  | ff19_opc_loca | 0.1399 -> 0.1406 | 0.1179 -> 0.1184 |
  | ff19_opc_locw | 0.1393 -> 0.1390 | 0.1174 -> 0.1163 |
  | ff19_opc, iterations 651-750 | 0.1429 -> 0.1425 | 0.1180 -> 0.1177 |
  | b5sd1 | 0.1186 -> 0.1185 | 0.1076 -> 0.1070 |
  | af_std | 0.1205 -> 0.1205 | 0.1067 -> 0.1071 |

  Alignment changes nothing. Before building, every map moves by 0.0007 or less. The round-3 changes (up to 0.0011)
  are the size of the water-pick noise: the aligned b5sd1 and af_std maps are identical to the raw ones (CC 1.000), yet
  their rounds moved by up to 0.0013.
- **What does separate the CA-only runs is the protein.** After the translation fit, the high-density protein atoms
  sit 0.92-1.14 A rms from their reference positions in the CA-only runs (locw 0.92, mult 0.98, loca 1.07, ff19_opc at
  iterations 651-750 1.14), against 0.27-0.47 A in the restrained runs (mean over every frame of the window, from the
  align logs). loca and ff19_opc stray most and make the two worst maps. The solvent organises around a protein
  surface that is not the crystal's, and translating the lattice back cannot fix that. Among the restrained runs,
  protein rmsd does not predict the map: k1B holds the protein closest (0.27 A) and makes the worst restrained map.

## 7KR0, PRISM's six-conformer ensemble (2026-11-04/05, octamus1)

`test_7kr0/`: Mac1, C2, 0.77 A, PRISM's ensemble (github prism-science/ag-refined-proteins, targets/7KR0): six whole
conformers at 1/6 each, R-free 0.0725.

- **MD:** `bulk_solvent_runme.com` on the ensemble (super_mult 1,3,1 = 12 copies, 2 per conformer). leap2amber lost the
  monoclinic beta angle; fixed on branch pressure (3f076d4, merged into main). Waters follow their conformer
  (altloc_waters=1).
- **Refinement recipe** (`scripts/run_ref.sh`): the agent recipe of `7KR0_SETTINGS.md` (10 cycles, wc/wu fixed,
  optimisers, N/Q/H flips, SA, ordered solvent and link residues off, 3 seeds, median) with explicit hydrogens.
  phenix 2.2 cannot switch Reduce off, and it times out on this model. `scripts/reduce_wrap/molprobity.reduce`
  replaces it for the clashscore call (no flips).
- **With PRISM's waters** (R-free, median of 3): ANISOU flat 0.0727, MD 2 ns 0.0709, MD 20 ns 0.0731; ANISOU stripped
  flat 0.0813, 0.0784, 0.0782. Restraint weight (0.3/0.1/0.03 per sigma) and the altloc merge distance (0.8 or 0.4 A)
  do not matter.
- **Scaling bistability:** with ANISOU and an MD map, phenix's k_iso/k_aniso fit has two solutions 0.002 apart in
  R-free, and any bulk-solvent step can land in either. Flat maps and isotropic models do not toggle.
- **Dry and rewet** (`scripts/rewet_track.sh`; waters deleted, then rebuilt conformer-aware by `scripts/ens_waters.py`
  from 6 sigma down, with the free set zeroed in the difference map):
  - Dry, the MD map leads by 0.055 (ANISOU flat 0.2245 vs MD 0.168).
  - The lead shrinks as water is built: 0.006-0.009 after round 1, under 0.002 by round 6.
  - Final R-free (20-cycle convergence): ANISOU flat 0.1164, md2 0.1181, md20 0.1169; isotropic 0.1113, 0.1112,
    0.1103. 388-436 waters built (occupancy sum 324-349, PRISM 545).
  - So MD and flat finish within 0.002 of each other, with mixed sign. The dry-model MD lead was water the rebuild had
    not yet placed.

## How our `bulk_solvent_runme.com` maps differ

| | 2019 RestAll | `bulk_solvent_runme.com` (1aho tests) |
|---|---|---|
| MD box | 2x2x3 supercell, 48 copies | 1x1x1 cell, 4 copies |
| length averaged | 100 ns, 500 frames | 2-20 ns |
| protein restraints | 10 kcal/mol/A^2, every heavy atom | rho x 1 kcal/mol/A^2 per sigma (cap 10); split altlocs free |
| water restraints | 2298 crystallographic water O, 10 kcal/mol/A^2 | none |
| pressure | constant volume, +450 to +800 bar | hydrated to ~0 bar, barostat on |
| water model | OPC | spce default (8 models tried) |
| temperature | 287 K | 287 K |

Complex correlation with the 2019 MD map, by resolution shell (8-5 / 5-3.5 / 3.5-2.5 / 1.9-1.4 / 1.4-0.96 A):
example 0.99 / 0.99 / 0.99 / 0.98 / 0.96, our opc 10 ns 0.96 / 0.90 / 0.81 / 0.68 / 0.59. Our maps agree with
each other at 0.92-0.97 in every shell, whatever the water model.

## Files

- `water_models/lineage/cmpmtz.py`, `fitB.py`, `fit2.py`, `fit3.py`, `col450.py`, `shells.py`: the numerical links
  in the table.
- `water_models/lineage/md2map_rerun/`: the byte-identical re-make of `solvent_RestAll.map` (intermediates deleted).
- `water_models/lineage/waters/`: restrained-water matching (`match_waters.py`, `counts_*.txt`) and the occupancy
  split (`occ_split.py`).
- `water_models/example_check/`: the example-protocol refinements (`ex_*`, `exocc_*`). `md2019_map.mtz` is the pure
  2019 MD map as FPART/PHIPART.
- 2019 sources: `~/projects/amber/1aho_test/` (`allmapme.com`, `md2map_multi.com`, `md2map.com`,
  `rerefine_expanded8/setup_mtz_notes.com`, `expand_solvent_notes.com`, `j5w8760/setup.com`,
  `j5w8760/runme_Rest.com`).
