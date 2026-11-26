# 1aho restraint weights, B factors and map probes (November 2026)

Findings from the restraint-weight and B-factor runs on 1aho, 4-24 November 2026, written for anyone building
solvent systems or bulk-solvent maps with this repo. Runs are under `/scratch/jamesh/projects/amber/1aho_hurry_test/`
(names `ff19_opc_*`), analysis tools under `/scratch/jamesh/projects/amber/tools_b5/`.

## What changed in the code that affects solvent builds

| commit | change | effect on `bulk_solvent_runme.com` and solvent maps |
|---|---|---|
| 660bced, 1ec6916 | Every density probe goes through `check_map_sites.com`: mapman, else a C spline (`Bfac_map peek`), else phenix. | `rholabel_runme.com` (Section 4 restraint weights, `w = rho x weight_per_sigma`, capped at `max_weight`) reads density at atoms **about 1.2x higher** at 1aho's 0.965 A than before, so the fixed restraints are stronger and a few more atoms pass `min_rho`. Water teleport/dehydrate/remap thresholds now mean what they say. |
| 45c802e | `leap2amber.com padwater_B` (default **50**, was 999). | Every water that hydration adds takes its padding slot's B. Any map rendered with those B values changes. |
| fab3bf0 | B floor 0 everywhere (`optimize_weights Bfac_min`, default 0; was 1 in the render and 2 in the B updates). | Nothing in `bulk_solvent_runme.com` itself; renders of per-atom B can now reach 0. |
| 660bced | `altloc_waters`: a water-letter count that does not match the water atoms now only switches the water filtering off. | Before, such a mismatch also cancelled the protein conformer handout (`altloc_copies`), with a warning that blamed the protein atoms. |

Runs started before these commits used the old behaviour (frozen code copies under `1aho_hurry_test/code_*`).

## 1aho's C-terminus is an amide, not a carboxylate

AaH II is C-terminally amidated: UniProt P01484 (SCX2_ANDAU), MOD_RES His-83 = histidine amide, which is His64 in the
mature numbering. `1aho.pdb` models a free carboxylate (O + OXT), and every 1aho topology before 2026-11-05 used CHIE
(charge -1). With the amide, each protein copy is **one unit more positive** (+4 per cell), so a neutral cell needs
four more anions. Keep the crystallisation salt and change only the counterions: the crystal grew from ammonium
acetate, 0.68 M in the example (`salt_conc`), which `add_salt_runme.com` turns into 11 NH4+ / acetate pairs per cell.
That gives 11 NH4 + 15 ACY with the carboxylate and 11 NH4 + 19 ACY with the amide.

The 2026-11-05 amide test below did not do that. Its carboxylate build had 0.2 M salt (ACY x7 + NH4 x3: 4 counter-
acetates plus 3 salt pairs). For the amide build I removed all three NH4+ and added one acetate (ACY x8, no NH4): the
neutral cell with the fewest new molecules to place. That was a convenience, not a requirement, and it took the salt
out of the amide build, so those amide-vs-carboxylate comparisons changed the ion composition as well as the
C-terminus. Salt alone moves pipeline R by about a point (0.05 vs 0.68 M: 32.5 vs 33.7% at uniform B). The 2026-11-25
charge tests (end of this section) repeat the comparison at fixed 0.68 M salt.

Evidence in the deposited model:
- The deposited His64 "O" sits 3.04 A from a symmetry mate's Thr57 carbonyl O. Two acceptors that close make sense
  only if that atom is the amide N donating an H.
- The deposited "OXT" accepts from Lys2 NZ (2.81 A) and the mate's Thr57 N-H (2.93 A), so it is the carbonyl O.
- With the carboxylate, His64 "O" had the largest Fo-Fc at its crystal position of any protein atom other than the
  disulfide sulfurs (+11 sigma), and Lys2 NZ needed one of the strongest restraints in the model. Released, Lys2 NZ
  moved 1.3 A to Asp53 O (Fo-Fc -10 sigma at the atom, +14 at its crystal site) because the -1 carboxylate pulls
  Thr57 OG1 onto itself.

Amide vs carboxylate in MD (same start, same protocol, 8 runs; amide with ACY x8 and no salt, carboxylate with 0.2 M
salt, see above). R and restraint energy over itr 112-268, the H-bond over itr 249-268 (first report, itr 92-111, in
`amide/reports/`):

| arm | R amide / cooh | restraint E | His64 "O" - Thr57' O (crystal 3.04 A), amide / cooh |
|---|---|---|---|
| weights frozen at k 1 | 38.80 / 39.16% | 383 / 391 | 3.0 A in 88-98% of frames / 3.3 A, ~50% |
| junction released | 38.96 / 39.09% | 379 / 391 | 3.0 A, 87-98% / 3.3-3.6 A, 5-51% |
| replicate | 38.91 / 39.18% | 380 / 389 | 3.0 A, 90-98% / 3.3-3.6 A, 6-53% |
| weight convergence (b5) | 32.86 / 33.53% | 80 / 90 | 3.0 A, 84-93% / 3.3-3.7 A, 12-45% |

The carboxylate broke the junction completely in one copy of each released arm (OXT - Thr57' N 4.1 A, 0% of frames);
the amide never did. In b5 the Lys2 NZ weight is 3.89 with the carboxylate (second highest of any protein atom) and
1.21 with the amide.

How the amide was built (`1aho_hurry_test/build_2s2r_amide`): `tleap_stub.in` redefines CHIE, the residue tleap uses
for a chain-final HIE, as internal ff19SB HIE plus the atoms and charges of Amber's NHE cap (N -0.4630, HN1 and HN2
+0.2315; net 0). The deposited "O" becomes the amide N and "OXT" the carbonyl O, so every atom and residue name in
the model, refpoint and restraint files still matches. The snippet is `1aho_starter_kit/his64_amide.tleap`; append it
to `tleap_stub.in` after the leaprc line. It is chemically the same as a separate NHE residue after His64, which
tleap also knows (aminoct12.lib), but that needs the deposited atoms renamed and adds a 65th residue per copy.
Gotchas: leap2amber's protonation file only knows PROTON/DEPROTON/HID/HIE/HIP (a new residue name there is
ignored), and leaprc.protein.ff19SB aliases CHIS to the library CHIE, so the snippet re-points CHIS too.

Nothing else in 1aho's restraint weights points at a protonation change at pH 6.8 (REMARK 280): the high weights sit
on backbone H-bonds (helix 21-26, turns 40-43 and 52-55, the 47-35 sheet) and on the disulfide sulfurs; Asp3 and Asp8
carry moderate weights but only receive H-bonds in the crystal. His54 and His64 touch only waters, so no contact
decides their tautomer, but the density does. H omit maps from the phenix bxfield gb3 model (another session,
2026-11-25; plain HIS with no ring N-H, re-refined without the H as a refined omit, plus add-H controls refined with
HID, HIE, HIP and His54-HIP/His64-HID; fractions of a backbone amide H, null sd 0.33):
- **His54 is doubly protonated (HIP).** Refined omit HD1 1.01, HE2 0.67. Whichever single H is modelled, the other
  returns as positive residual (+0.82, +1.11); with both modelled both residuals sit at the fit's zero point.
- **His64 is HID.** HE2 is absent (refined omit +0.13; modelled, it leaves -0.74 to -0.78). HD1 reads about 0.6 of
  an H, partial or more mobile than a riding B allows.
- R-free cannot see the ring N-H (within 0.0002); the difference density can. A second model (best_so_far, unrefined
  omit) agrees. Ring angles carry no information: the restraints hold both rings at His+ geometry.

All 2026 MD builds before these tests used HIE54 / HIE64, which the omit maps rule out for both residues.

**Charge tests at the crystallisation salt (launched 2026-11-25).** One 0.68 M ammonium acetate placement
(`add_salt_runme.com`, 11 NH4+ / acetate pairs per cell) shared by every arm; only the counter-acetates follow each
arm's charge. Build directories `1aho_hurry_test/build_s68_<arm>`, runs `ff19_opc_s68_<arm>`, all run as the
restraint-class reference arm (b5sd1 settings with `min_align_weight=0`, uniform B 5, 200 itr):

| arm | C-terminus | His54 | His64 | cell charge before counterions | NH4 / ACY |
|---|---|---|---|---|---|
| `cEE` | carboxylate | HIE | HIE | +4 | 11 / 15 |
| `EE` | amide | HIE | HIE | +8 | 11 / 19 |
| `ED` | amide | HIE | HID | +8 | 11 / 19 |
| `PD` | amide | HIP | HID | +12 | 11 / 23 |

`PD` is the omit maps' answer; `EE` the old default; `EE` vs `cEE` is amide vs carboxylate at fixed salt, `ED` vs
`EE` the His64 tautomer, `PD` vs `ED` the His54 charge. Two more builds were made but not run once the refined omit
maps came in: `DD` (HID54 / HID64, +8) and `MD` (His54 HIP in copies 1 and 3, HID in 2 and 4; +10, 11 / 21).
The extra counter-acetates sit on the placed waters farthest from protein and other salt (nested: `EE`'s 4 are `PD`'s
first 4). For His64 as HID, `tleap_stub.in` redefines CHID as the amide the same way as CHIE (HID and HIE share their
C and O template coordinates). Scripts: `amide/salt068_base.csh`, `setup_s68_builds.sh`, `build_s68.csh`,
`launch_s68_tests.sh` next to `1aho_hurry_test` (copies in `bulk_solvent_bench/amide/`). Because the MD maps include hydrogens, Fo-Fc at the modelled ring N-H
(negative where an H is wrong) is a direct check alongside R and the restraint weights.

## B factors

- **Water B 999 is not invisible.** Rendered, it is a blob about 8 A wide. Fill waters at B 999 cost R: in one
  test, 54.2% with them at 999 vs 43.9% at B 10 (62.5% with them deleted). Hence the new `padwater_B` 50.
- **Free per-atom B helps once the weights are fixed.** From b5 at itr 200 (weights converged at uniform B 5,
  R 33.9%), freeing per-atom B (`Bfac_modmode=mult Bfac_maxmod=2`, no smoothing) gave R 31.7% within 5 iterations.
  B got sharper, not softer: by itr 230 half the protein heavy atoms and over a third of the waters were on the
  old B = 2 floor, and the distribution was still splitting. That is why the floor is now 0. The MD ensemble
  supplies the motion, so rendered B near 0 is physical.
- **The disulfide sulfurs were rendered too soft.** At uniform B 5 they sat in +14 to +16 sigma Fo-Fc peaks at
  their targets and +5 sigma at the atom itself: the residual of a heavy, sharp atom drawn blurry, not a misplaced
  one. Their restraint weights climbed to the cap, but a restraint cannot sharpen an atom. With free B the residual
  at the atom dropped to about +1 sigma.
- **The default render is blurrier than its B factors; the GPU render is not** (tested 2026-11-25). gemmi's
  structure factors handle B = 0 (`gemmi sfcalc --test` on a 1aho frame, d 0.965, rate 1.5: FFT vs exact R 0.002%
  with every atom at B 0). But `nc2mtz_gemmi.com` sums gemmi's `--write-map` maps, and those keep gemmi's
  anti-aliasing blur: every atom gets max(0, 7.43 - lowest B) A^2 more at these settings (gemmi reports `B_add`:
  2.43 at uniform B 5, 7.43 once any atom reaches B 0, 0 at uniform B 10). So every uniform-B-5 run here was
  rendered at B 7.4, and the B floor of 0 renders at 7.4. `nc2mtz_gpu.com` (`sfcalc_gpu_collapse`, which removes
  its own blur in reciprocal space) has none. On rc_dw_w's itr-199 frames (25, uniform B 5), the GPU average matched
  the mean of gemmi's own de-blurred per-frame F to R 0.005% (phases 0.002 deg), and matched the gemmi render to R
  0.004% once blurred by 2.43 A^2. Against Fobs both give R 26.83%, with fitted overall B -2.3 A^2 (gemmi render)
  and +0.2 (GPU). A fitted overall B absorbs a uniform blur, so pipeline R does not change; what is off is the
  meaning of the B in the file, and the blur shifts whenever the lowest B does. Renders at uniform B 10, such as
  the solvent session's R-free test, carry none. The GPU render is now `optimize_weights render_engine=gpu`
  (draft option, default gemmi). Tested end to end on build_s68_PD (3 itr, 2026-11-25): on the same itr-0 frames
  R 49.89% (GPU) vs 50.03% (gemmi), fitted B -1.98 vs -4.54; itr 1-2 (different MD from there) 43.98 / 37.79% vs
  43.46 / 37.11%. The 25 GPU frames took 14 s per render (36 s on the first, cold).
- **Overall scale B** (scaleit's correction to F; negative = model too blurry): about -1.3 A^2 for uniform B 5
  runs (most of it gemmi's 2.4 A^2 blur, above), -5 to -6 for runs with per-atom B starting from the 999 fill
  waters.
- **B by location (`Bfac_map=refine`):** before 660bced its field update moved every voxel, so a water-only field
  was driven to the floor inside the protein (64-89% of voxels below B 5 by itr 199) and never beat per-atom B. The
  update is now masked to where the selected atoms went in the frames (5-13% of voxels); not retested since.

## Map probes

- `phenix.map_value_at_point` samples its own 0.5 A grid by default and interpolates linearly. At 1aho's 0.965 A
  that reads Fo-Fc at atom sites 27% low (slope 0.72 against a d/12 tricubic reference). It also prints its points
  in its own hierarchy order (amino acids N CA C O CB ..., not Amber's N CA CB ... C O; 1462 of 2785 atoms of a 1aho
  frame came out of order), and its output format has changed between versions. Scripts that counted its output
  lines got the right answers only because they probed water O atoms, one atom per residue.
- `check_map_sites.com atoms.pdb map|mtz` is now the one probe: one value per atom, in atom order. It tries mapman
  (`peek value ... spline`; fast, not installed everywhere), then `Bfac_map peek interp=cubic` (the same 4-point
  cubic as mapman, compiled from `Bfac_map.c`; matches mapman to 0.0005 sigma), then phenix. Given an mtz it ffts it
  on a d/4 grid (`grid_sample`) and sigma-scales it.
- Accuracy against the reference, Fo-Fc at 2785 atoms (rms 2.52 sigma):

  | sampler | grid | slope | rms error (sigma) |
  |---|---|---|---|
  | phenix, old default | 0.5 A (d/2) | 0.715 | 0.80 |
  | mapman / Bfac_map peek cubic | d/3 | 0.966 | 0.10 |
  | mapman / Bfac_map peek cubic | d/4 | 0.987 | 0.04 |
  | mapman / Bfac_map peek cubic | d/6 | 0.997 | 0.013 |
  | trilinear | d/6 | 0.961 | 0.11 |

## Restraint weights

- **Only restraints with k >= `cutoff_weight` reach the MD.** At the end of b5 (weight_scaledown 0.9 for 200
  iterations, uniform B 5, R 33.6%) that was 409 of 2258 restraints with `cutoff_weight` 0.1. Without the decay
  (`weight_scaledown=1`, b5sd1) about twice as many stay, with R 32.9%.
- **The most-restrained atoms earn their weights.** Freezing b5's weights and releasing one top atom per run in all
  four copies (Trp38 O, Ala45 CB, Lys58 CG, Pro52 O, Val55 N, Lys2 NZ, Lys50 O, Ser33 N): each drifted 0.05-0.15 A
  further off its target and Fo-Fc at the target rose 1-4 sigma, while R stayed within noise (sd 0.43) for all
  but Lys2 NZ (+1.4), which the amide explains.
- **Backbone-only restraints fit best by pipeline R** (restraint-class arms, 2026-11-25; no-salt amide build,
  HIE54/HIE64, uniform B 5, one run each, itr 101-200; table in the joint summary,
  `bulk_solvent_bench/MD_BULK_SOLVENT_JOINT.md`): backbone N/CA/C/O + waters 31.44%, every heavy atom 32.41%
  without waters and 32.66% with them, CA only at k 1 34.85%, everything at fixed k 10 40.27%. Free side chains fit
  better than restrained ones; restrained waters change nothing within noise.
- **Weights swing a lot between iterations** under decay 0.9: Ser40 O averaged k 2.4 over itr 191-200 but was 0.1
  at itr 200.
- **What it takes to hold a weight.** New weight = old x `weight_scaledown` x `max_mult`^signif, signif =
  (dFo-Fc / halfrho)^5 clipped to +-1, so a weight holds only if dFo-Fc >= halfrho x (ln(1/sd) / ln(max_mult))^(1/5).
  With 1aho's auto halfrho of 3.97 sigma: 2.72 sigma at sd 0.9 and max_mult 2, 2.37 at max_mult 4, 2.14 at 10;
  1.70 at sd 0.99 and max_mult 2. Because of the 5th power, max_mult moves it little; scaledown and halfrho move it
  more.

## Site environment after the November 2026 outage

- /home is now petabyte:/petabyte/home, and /programs (-> /home/programs) is its copy of the program tree, still being
  synced: `/programs/phenix` points at phenix-1.21, and its phenix-2.2-6143 is the Aug 12 install where the complete
  pre-outage tree, /data/home/programs, has the Aug 26 one. /data is mounted on octamus1, voltron, graphics2 and bl831.
- So take programs from /data/home/programs: source `/home2/jamesh/projects/amber/bulk_solvent/programs.sourceme`
  (tcsh) or `programs.sh` (bash), which set up amber22, phenix 2.2 and ccp4 9 from `$PROGRAMS` (default
  /data/home/programs; it works around the setup files' hard-coded /programs) and then source the site
  `compute_settings.sourceme` for the MD engines. The 1aho wrappers `run_2s2r.csh` and `run_local.csh` do.
- cpptraj 7.12.6 (`$PROGRAMS/cpptraj/7.12.6`) is needed for ff19SB topology strips; Amber22's own cpptraj drops the
  CMAP terms.
- optimize_weights stopped making progress on the 1aho runs after the outage: pmemd's timing line
  (`seconds/ns = **********`, when the last step is slow, e.g. writing to NFS) tripped its Fortran-overflow check,
  and the stage's combined .out kept the line across retries, so every retry failed the same way. Fixed: the check
  skips the timing lines and each attempt starts the file empty.

## tcsh traps met along the way

- A one-line `if ( cond ) cmd` expands every variable on the line and performs any redirect even when cond is false:
  `if ( $?X ) set y = "$X"` fails with "X: Undefined variable" when X is unset. Use if/then/endif. (An error inside a
  sourced file stops the rest of that file, not the caller.)
- `$var:` followed by a space is read as a variable modifier ("Bad : modifier"); write `${var}:`.
- ``set v = `cmd` `` splits the output into words and joins them with single spaces, so column-aligned text loses its
  spacing; compare columns with awk instead of matching formatted text.
- tcsh reads a running script by byte offset: never edit a script that is running in place; write a new file and
  `mv` it over.
