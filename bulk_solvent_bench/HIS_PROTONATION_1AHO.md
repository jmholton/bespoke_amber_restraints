# His protonation in 1aho: best guess as of 2026-11-28

Two lines of evidence, each with its own blind spots:

- **X-ray** (the refinement session, 25 Nov): H omit maps and add-H controls on two refined models (bxfield gb3 and
  best_so_far), plus a His54 ring-flip test. This evidence sees the H directly but is noisy: one H is 2-3 null sd.
- **MD** (this session, 25-27 Nov): 16 restrained ff19SB/OPC runs of the P212121 cell (4 copies), crystallisation salt
  0.68 M ammonium acetate, every His54 × His64 combination, most in duplicate; 13 more runs (launched 27 Nov) are
  at itr ~123 of 200, with interim results below. This sees the H only indirectly:
  through how far each atom (and nearby waters) must be pulled to stay in its density. It depends on the force field
  and solvent model.

## Best guess

| | ring orientation | C-terminus | tautomer | confidence |
|---|---|---|---|---|
| **His54** | M (χ2 about −100°) | – | **HD1 present; HE2 partial**: a HID / His⁺ mixture, mostly HID | orientation firm; HD1 firm; HE2 weak |
| **His64** | one orientation (χ2 about +60°) | **amide** | **HID** (HE2 absent); HD1 possibly partial | amide firm; HE2 absent firm in X-ray, contested by MD |

For MD the best single choice is **HID54 / HID64 (`DD`)**. It keeps every firm X-ray result: HD1 on both His, no
HE2 on His64, M ring. It drops only the weakest one, His54 HE2 at 1.2 null sd. It also avoids the one large MD
penalty, the strain of a fully charged His54. Its weak point is His64, where MD prefers HIE (below). `DD` now has a
replicate (`DD_r2`, interim).

**New since 27 Nov (interim, itr 101-121):** an NH4+ placed on crystal water HOH51 stays H-bonded to His64 ND1 only
when ND1 carries no H (HIE64); with HID64 it is pushed away at once. So the MD's preference for HIE64 and the
"NH4+ at HOH51" idea go together, and the X-ray HD1 peak decides between them (analysis 3 below).

## His54

**Ring orientation: M, firm.**
- The X-ray ring N/C swap test wants every P-orientation altloc relabelled. It does so by 15-24 noise sd in every H
  treatment. M sits at the baseline [measured].
- best_so_far's 8 χ2 > 0 altlocs (occupancy 0.42) are mislabelled; it is one ring, not a flip mixture [inferred].
- Every MD build started in M.
  - It stayed in M in every copy of every run, except one: copy 3 of the carboxylate control `cEE` flipped to P for
    86 % of frames (itr 101-199).
  - `PDp` was built in P on purpose and stayed there.

**HD1: present, firm.**
- In the refined omit maps (one orientation, M) HD1 reads +1.2 of an H.
- Every arm without HD1 leaves it as +0.9 of an H of residual.
- gb3 agrees: +1.0 at 3.1 null sd [measured].

**HE2: partial, weak.**
- X-ray:
  - In the M refined omit, HE2 reads +0.59 ± 0.50 [measured].
  - The HID arm leaves +0.94 at HE2 [measured].
  - The HIP arm leaves both N-H at −0.3 / −0.2, within one null sd [measured].
  - So X-ray ranks His⁺ first, but cannot separate full from partial HE2 [inferred].
- MD: a fully charged His54 does not sit where the density puts it [measured].
  - Restraint energy w·d² per copy at the ring's outer edge, itr 101-199:

    | His54 in the run | NE2 | CE1 |
    |---|---|---|
    | HIP in all copies (`PD`, `PD_r2`, `PDp`, `PE`, `PE_r2`) | 0.33-0.46 | 0.26-0.36 |
    | mixed: HIP in copies 1+3, HID in 2+4 (`MD`, `ME`, two runs each) | 0.18-0.23 | 0.15-0.18 |
    | HIE or HID (8 runs) | 0.00-0.13 | 0.00-0.10 |

  - The strain scales with the number of charged copies.
  - The weights at those atoms rise from about 0.05 to 1-2 (the density insists).
  - The charged ring draws more water: 1.5 vs 1.0 waters within 3.2 Å of NE2, and 0.5 vs 0.2 at CE1. Acetate touches
    it in under 4 % of frames.
  - The `PDp` control (wrong ring labels, same charge) shows the same strain, so it comes from the charge, not the
    labels.
- Pipeline R also groups by His54 charge [measured].
  - Values: neutral His54 32.98-33.32 %, mixed 33.50-33.86 %, all HIP 33.71-33.99 %. Replicates agree within 0.17.
  - But it is confounded: each HIP copy adds a counter-acetate (19 → 21 → 23 per cell), so this is not clean evidence.
- An early MD-render Fo-Fc look (itr 101-199, four arms) agrees: His54 ring rms 1.52 for HIP vs 1.17-1.36 for
  HIE [measured].

**Reading [inferred].**
- HD1 is certain. HE2 is at most partial: the X-ray wants some, and MD rejects a full charge on every copy.
- A His near its pKa in the crystallisation buffer, mostly HID with a minority His⁺, fits both. The actual
  crystallisation pH is not recorded here. Unbuffered ammonium acetate is near pH 7, above His's usual pKa of about 6.5.
- Partners: in the 1aho model a water sits 2.81 Å from ND1. It also bridges to Gln37 OE1 of a neighbouring copy
  (copy 4, residue 229). This suits ND1-H donating to it [inferred].
- NE2 has a water at 3.20 Å, which is nearer to Leu29 O (2.88 Å).
- Nothing in the structure demands a charge.

## His64 (C-terminus)

**C-terminal amide: firm.**
- UniProt P01484 lists His64 as amidated.
- In MD the carboxylate is the worst model of the C-terminal junction:
  - His64 O carries 0.37 as the carboxylate vs 0.00-0.07 in all 15 amide runs.
  - Thr57 O and three nearby waters behave the same way.
  - Lys2 NZ improves only partly (0.49 → 0.24-0.44).

**Tautomer: HID by X-ray, HIE by MD.** This is the open conflict.

- X-ray [measured]:
  - HE2 is absent in every test of both models: refined omit +0.13 / +0.25.
  - A modelled HE2 leaves −0.74 to −1.22, about 2.6 null sd below zero.
  - HD1 is present at 0.58 (gb3) to about 1.0 (best_so_far).
  - The carry-forward refinement `bsfflip_Mmix` has His54 M HIP and His64 HID. Its His residuals are all flat.
- MD [measured], energy per copy over itr 101-199:

  | runs | His64 NE2 | water about 2.9 Å from His64 ND1 | whole cell |
  |---|---|---|---|
  | HIE64 (`DE`, `PE`, `ME` ×2 each, `EE`, `cEE`) | 0.16-0.31 | 0.01-0.05, 0.15 in `PE_r2`, absent in 2 runs | 127-139 (`cEE` 151) |
  | HID64 with HIE54 (`ED` ×2) | 0.29-0.32 | 0.15-0.30 | 137, 148 |
  | HID64 with HD1 on His54 (`DD`, `PD` ×2, `PDp`, `MD` ×2) | 0.45-0.56 | 0.17-0.23 | 135-150 |

  - With HID64 the His64 ring edge (NE2, CE1) is strained.
  - A water restraint that repick placed about 2.9 Å from ND1 is also strained. It carries weight 6-8, and is not one of
    the two crystal waters.
  - The His64 strain is worst when His54 also carries HD1, a cross-talk whose cause is open. It is not a rotamer change:
    His64 stays χ1 trans and χ2 +56 to +73° in every frame. Nor is it a ring-label assignment.
  - But the `DD` replicate does not repeat it so far: His64 NE2 0.19 in `DD_r2` vs 0.53 in `DD` over itr 101-121
    (interim). The cross-talk may be run-to-run scatter.
  - `DE` (HID54 / HIE64) has the lowest whole-cell pair, 127 and 129.
  - Replicates scatter by up to 11, so the totals only lean.

**Reading [inferred].**
- The X-ray result is direct and replicated across two models and every arm. The MD result is indirect, and the
  energies involved (0.2-0.3 per copy) are below the strongest restraints in every run (0.4-0.6, at Arg62 CZ, Arg18
  NH1, Asp8 OD2, Lys58 O).
- So HID64 stands as the best guess.
- The MD dissent most likely points at the solvent around ND1, not at the tautomer itself:
  - Crystal HOH51 sits 3.04 Å from ND1. It is 2.79 Å from Tyr42 O and has two water neighbours: four tetrahedral
    partners.
  - That is equally consistent with a water accepting from ND1-H (HID) or with an NH4+ donating to an unprotonated ND1
    (HIE). Crystals grew from ammonium acetate, and NH4+ is isoelectronic with water.
  - The second crystal water, HOH58 (3.12 Å from ND1, 3.04 Å from Arg62 NH2), is not a cation site next to an Arg. It is
    4.15 Å from HOH51, so the pair is not an acetate.
  - best_so_far has partial waters (occupancy 0.10-0.11) 1.93 and 2.39 Å from ND1. That is too close for an H-bonded
    water, but about where the H of a donor to ND1 would sit.
- A partial HD1 (0.58 in gb3) would fit a minority of His64 with an NH4+ on ND1 instead of an H.
- **MD test of NH4+ at HOH51 (interim, itr 0-125)** [measured]: NH4+ replaced HOH51 in all 4 copies.
  - With HIE64 (`DE_n`) NH4+ holds 2.8-3.0 Å from ND1 and 2.7-2.9 Å from Tyr42 O. Over itr 100-125 it is within
    3.3 Å of ND1 in 80 % and 72 % of iterations in copies 1 and 3. Copy 2 left at itr ~30 and copy 4 at ~75 (4 %
    and 8 %).
  - With HID64 (`DD_n`) NH4+ starts already pushed to 4.0 Å from ND1 after equilibration, keeps Tyr42 O briefly, and
    is 5-7 Å away by itr 20. A water takes the ND1 site at about 2.9 Å.
  - So an NH4+ at this site requires an unprotonated ND1 [measured, in MD]. The X-ray "HD1" peak would then be the
    NH4+ H, about 2.0 Å from ND1, not 0.86 Å [inferred].
  - Cost: `DE_n` and `DD_n` score 10-25 higher whole-cell than `DE`/`DD` in the same window. Of that, the His64
    C-terminal O and amide N (OXT) add 0.2-0.4 per copy and the NH4+ restraint 0.38. The rest is spread over nearby
    waters. NH4+ on this site in MD fits the density no better than water.
- Chloride at full occupancy on these sites is unlikely. Cl⁻ has nearly twice water's electrons, and HOH51/HOH58 rank
  only 17th and 21st of 121 waters by density. Partial chloride is not excluded.

## What R-free and pipeline R say

- Refinement R-free cannot see the ring N-H: all H arms lie within 0.0005 [measured].
- Pipeline R in MD moves about 0.6 points with His54's charge, but that is confounded by the counter-acetates.
- Neither decides the tautomers; the difference density and the restraints do.

## New MD tests: interim results (itr 101-121; runs launched 2026-11-27 18:50, at itr ~123 of 200)

All runs use the same settings as above; run directories are `1aho_hurry_test/ff19_opc_<name>`. These are 21
iterations against 99 in the main window; differences under about 0.15 per copy, or about 10 whole-cell, are noise.
Full scoring over itr 101-199 follows when they reach itr 200.

The runs:
- **Replicate:** `s68_DD_r2`, a second run of the best-guess arm.
- **Same-salt controls** (`amide/setup_s68_samesalt.sh`). Every arm gets DD's salt unchanged (11 NH4 / 19 ACY). Extra
  His54 charge is left to PME's uniform neutralizing background, not countered with acetates. This removes the acetate
  confound in the HIP comparisons. Arms: `s68_QD_u` (HIP54 in 1 copy of 4, HID64; net +1), `s68_MD_u` (HIP54 in 2 of
  4; +2), `s68_PD_u` (HIP54 in all, HID64; +4), `s68_PE_u` (HIP54 in all, HIE64; +4).
- **NH4+ at His64 ND1:** crystal water HOH51 replaced by NH4+ in all 4 copies (the 4 NH4+ farthest from protein were
  moved, so the salt is unchanged). `s68_DE_n` is HIE64; `s68_DD_n` (HID64, NH4+ facing an N-H) is the control.
  Results are under His64 above.
- **Salt series** (`amide/setup_salt_series.sh`) on `DD` and `DE`: `s017_*` 0.17 M (2 NH4 / 10 ACY), `s034_*` 0.34 M
  (5 / 13), `s136_*` 1.36 M (23 / 31), beside the 0.68 M runs (11 / 19).

**Same salt vs acetate-countered** (energy per copy; whole cell per iteration; `DD` and `DE` as run / replicate):

| | `DD` | `DE` | `QD_u` | `MD` | `MD_u` | `PD` | `PD_u` | `PE` | `PE_u` |
|---|---|---|---|---|---|---|---|---|---|
| His54 NE2 | 0.00 / 0.06 | 0.00 / 0.11 | 0.09 | 0.18 | 0.15 | 0.44 | 0.16 | 0.23 | 0.30 |
| His64 NE2 | 0.53 / 0.19 | 0.16 / 0.21 | 0.31 | 0.44 | 0.38 | 0.68 | 0.74 | 0.26 | 0.38 |
| whole cell | 145 / 141 | 131 / 127 | 121 | 143 | 136 | 154 | 155 | 132 | 120 |

- His54 NE2 strain still rises with the number of charged His54 without the extra acetates: neutral, then 1 of 4
  (`QD_u`), then 2 of 4, then all. The step is smaller: `PD_u` 0.16 vs `PD` 0.44, so part of `PD`'s strain may come
  from the acetates.
- `QD_u` (one His⁺54 in four, HID64) is among the lowest whole-cell so far, level with `PE_u` and within noise of
  `DE`. A minority His⁺54 costs nothing measurable, consistent with "mostly HID".

**Salt series** (whole cell; His64 NE2 per copy):

| salt | `DD` | `DE` |
|---|---|---|
| 0.17 M | 134; 0.36 | 124; 0.16 |
| 0.34 M | 122; 0.39 | 131; 0.20 |
| 0.68 M | 145 / 141; 0.53 / 0.19 | 131 / 127; 0.16 / 0.21 |
| 1.36 M | 141; 0.38 | 147; 0.23 |

- HIE64 gives the lower His64 NE2 strain at every salt concentration.
- At 1.36 M His54 NE2 is strained even when neutral (0.18 `DD`, 0.32 `DE`; ≤ 0.11 at the other salts).
- No whole-cell trend with salt yet. Lys2 NZ stays at 0.23-0.34 throughout.

## Difference-density analyses not yet done (refinement side)

The MD Fo-Fc maps are too noisy for H, so these need the refinement models.

1. **His54 HE2 occupancy.** Scan or refine HE2 occupancy (0-1) in the M model (`bsfflip_Mmix`) and find where the
   residual crosses the zero point. This sets the HID : His⁺ ratio directly. MD can only mix by copy (1, 2 or 4 of 4).
2. **His64 HD1: partial or mobile.** Run the same scan for His64 HD1, plus a refined H B. That separates "0.6 of an H"
   from "a full H with a larger B".
3. **Which side holds the H between His64 ND1 and HOH51.** Now the deciding test: the MD says NH4+ can sit on HOH51
   only if ND1 has no H.
   - HID64 with water puts the H at the HD1 site, 0.86 Å from ND1.
   - HIE64 with NH4+ puts it about 1.0 Å from HOH51, about 2.0 Å from ND1. That is near best_so_far's partial "waters"
     at 1.93 and 2.39 Å.
   - The two sites are about 1.1 Å apart, separable at 0.96 Å. Fit both H templates jointly in the refined omit.
   - Also test the other H of a water vs an NH4+ on HOH51 (2 vs 4 H, tetrahedral).
4. **Chloride.**
   - No anomalous test is possible: the deposited data (`1aho-sf.cif`) are merged F only, measured at 0.92 Å.
   - Instead, refine HOH51, HOH58 and the other strong waters as partial Cl⁻ against water, and compare occupancy, B and
     the residual.
5. **Ring angles refined free.** Refine the His ring angles with their restraints relaxed (large sigma). At 0.96 Å
   the angle at each N (105° unprotonated vs 109° protonated) may be measurable. The 105.2° edits refined back to
   107.1°, so the data have something to say.

## Sources

- **X-ray:** refinement session report of 2026-11-25 (`his_omit.py`; `his_omit_gb3.txt`, `his_omit_bsf*.txt`,
  `flip_{P,M,X}*.txt`, `flip_Mmix.txt`; bl831 `maperror_protmodel_2026-11-25/bsfhis`, `bsfflip`). Null sd is in units
  of one H after refinement.
- **MD runs:** `/scratch/jamesh/projects/amber/1aho_hurry_test/ff19_opc_s68_<arm>`.
  - All runs: b5sd1 settings, code_391d233, `min_align_weight=0`, uniform B 5, 200 itr, scored over itr 101-199.
  - Builds: `amide/setup_s68_builds.sh`, `setup_s68_flip.sh`, `setup_s68_more.sh`, `setup_s68_samesalt.sh`,
    `setup_salt_series.sh`.
- **MD measures** (in `/scratch/jamesh/projects/amber/amide/`, mirrored in bespoke_amber_restraints
  `bulk_solvent_bench/amide/`):
  - `worst_restraints.py`; tables `worst_restraints_s68_16arms_101-199.txt`, `worst_restraints_interim_101-121.txt`;
  - `nh4_site.py` (NH4+ and water distances to His64 ND1 / Tyr42 O per iteration); table `nh4_site_s68_n.txt`;
  - `his_chi_nc.py`, `his_neighbours_md.py`, `release_stats.py`, `his_fofc.py`.
- **Joint doc:** `MD_BULK_SOLVENT_JOINT.md` has the run-by-run detail.
