# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`bespoke_amber_restraints` is a collection of command-line tools for using X-ray crystallography data to iteratively optimize harmonic restraints in AMBER molecular dynamics simulations. The goal is to find the minimum restraint weights needed to keep a structure consistent with its electron density.

## Compiling the C Programs

Three C programs are included and must be compiled before use:

```bash
gcc -o float_add float_add.c -lm
gcc -o float_func float_func.c -lm
gcc -O3 -o Bfac_map Bfac_map.c -lm
```

These have no dependencies beyond the standard C math library.
(`Bfac_map_runme.com`, `nc2mtz_gemmi.com` and `nc2mtz_gpu.com` compile `Bfac_map`
themselves if it is missing.)

## Language and Script Conventions

- All scripts use the `.com` extension and are written in **tcsh** (`#!/bin/tcsh -f`)
- AWK programs use the `.awk` extension and are self-contained utilities
- Scripts accept `key=value` command-line arguments; the argument parsing pattern checks for an `=` sign and sets a tcsh variable matching the key name
- Configuration is read from `sourceme` files (e.g., `xtal_properties.sourceme`, `settings.sourceme`) via tcsh `source`
- Intermediate files use a `$tempfile` prefix (often `/dev/shm/${USER}/temp_<script>_$$_`) and are cleaned up on exit
- Error handling uses a `$?BAD` pattern: set `BAD = "message"` and `goto exit`; the exit label checks `$?BAD` and exits with status 9

## Architecture

### Master Script

`optimize_weights_runme.com` is the main entry point. It:
1. Reads `xtal_properties.sourceme` and optionally `settings.sourceme` for crystal/run parameters
2. Iterates: runs AMBER MD → generates maps from structure factors → analyzes Fo-Fc difference map → adjusts restraint weights → repeat
3. Calls most other scripts in this repo as sub-steps

### Crystal Configuration (`xtal_properties.sourceme`)

This file defines the crystal system and must be created before running the master script:
```
set reso = 1.0          # resolution in Angstroms
set smallSG = P212121   # space group
set super_mult = 1,1,1  # supercell multiplier
set modulo = 64         # number of residues in one protein chain
```

### Required External Programs

- **AMBER** (18+): `pmemd.cuda_SPFP`, `cpptraj`, `tleap`, `AddToBox` — sourced from `/programs/amber22/amber.csh`
- **CCP4 Suite**: `mtzdmp`, `cad`, `mapmask`, `sfall`
- **Phenix Suite**: `phenix.refine`, `phenix.fetch_pdb`
- **gemmi**: for MTZ/map operations (`nc2mtz_gemmi.com`)

### Key Input Files

The master script expects these files in the working directory:
- `../refme.mtz` — structure factor MTZ file with columns FP, SIGFP, FreeR_flag
- `../centroids/centroids_in_density.pdb` — reference points in electron density
- `${template_dir}/amberme.pdb` — starting model
- `../*.mol2` and `../*.frcmod` — ligand parameter files
- `tleap_stub.in` — tleap input template

### Script Categories

**PDB manipulation** (AWK): `convert_pdb.awk`, `filter_pdb.awk`, `reformatpdb.awk`, `sequence.awk`, `jigglepdb.awk`, `build_c2n.awk`, `build_n2c.awk`

**Map/MTZ handling**: `nc2mtz_runme.com`, `addup_maps_runme.com`, `addup_mtzs_runme.com`, `map_scaleB_runme.com`

**Restraint management**: `restraintlist2amber.com`, `restraintlist_update_diffmap.com`, `delete_worst_restraints.com`, `release_worst_restraints_runme.com`, `no_new_nonbonds_runme.com`

**Water/hydration**: `hydrate_runme.com`, `dehydrate_amber_runme.com`, `add_waters.com`, `water_teleport_runme.com`, `remap_waters_runme.com`

**Geometry analysis**: `distanceify_runme.com`, `centroids_nearby_runme.com`, `measure_voids_runme.com`, `molprobify_runme.com`, `quick_geo_check_runme.com`

**B factor refinement**: `Bfac_update_diffmap.com`, `map_scaleB_runme.com`, `scaleB_search_diffmap_runme.com`, `thrubond_avgB_runme.com`, `rmsd2B`, `Bfac_map`, `Bfac_map_runme.com`

### B factors as a spatial field (`Bfac.map`)

Normally B factors live per atom in `Bfac.pdb`, applied to each trajectory frame
by atom order.  `Bfac_map` keeps them as a continuous field `B(x,y,z)` in a CCP4
map instead, so B belongs to the location rather than to the atom: a water that
lands on an ordered site is immediately sharp, one that wanders off is immediately
diffuse, and nothing has to be kept in step with atom order.

```
Bfac_map build pdb=Bfac.pdb outmap=Bfac.map [sigma=0.5] [grid=0.5] [farB=999] [fardist=4.5] [sg=P1]
Bfac_map probe pdb=frame.pdb map=Bfac.map [outpdb=-] [minB=] [maxB=]
Bfac_map stats map=Bfac.map
```

- The field is `B(x) = (SUM_i w_i B_i + w0 farB) / (SUM_i w_i + w0)` with
  `w_i = exp(-r^2/2 sigma^2)` and `w0 = exp(-fardist^2/2 sigma^2)`: the local
  average of nearby atomic B, decaying to `farB` where there is no atom within
  `fardist`.
- The map holds **final** B values — `Bfac_map build` applies
  `Bscale`/`Boffset`/`minB`/`maxB`, and `probe` just interpolates. `farB` is
  never clipped by `maxB`.
- **Build the field from coordinates where every atom is real.** A production
  `Bfac.pdb` is a *B sidecar*, applied to frames by atom order — its own
  coordinates can be placeholders. In `6c2r_37C_2x/opt84/Bfac.pdb`, **400 000 of
  801 634 atoms are stacked at the origin** (100 000 padding waters × O/H1/H2/EPW),
  a spatial field cannot represent a stack of different B values at one point,
  and building from that file gave rms ΔB 12.9 Å² and **R on F 10.2%**. Building
  from the same B values on real coordinates gave rms ΔB **0.47** and R on F
  **0.60%** (401 634 atoms, P41212 2×2×2, dmin 1.9, bimodal B 2–100). So
  `Bfac_map=auto` merges the per-atom B onto a trajectory frame first and builds
  from that — which also makes the atom-order merge a once-per-stage step instead
  of once per frame. `Bfac_map build` warns when it finds stacked atoms.
- **`sigma` is the parameter that matters.** Against the per-atom path (1aho
  expanded to P1, `dmin` 1.0) R on F was 0.1% at 0.25, 0.3% at 0.35, 1.2% at 0.5,
  3.4% at 1.0, 15% at 1.5. Voxel size matters far less — 1.1% at 0.25 Å, 1.2% at
  0.5, 1.8% at 1.0, 2.4% at 1.5 — so a large supercell can use a coarse grid to
  keep the file small (4 bytes/voxel; `grid=1.0` is 8× smaller than `grid=0.5`).
  On a bimodal model (core B=8, solvent B=90) `sigma` 0.5 reproduced F to 0.01%.
- **Cost.** The 801 634-atom production supercell takes 106 s to build at
  `sigma=0.5 grid=0.5` (159 MB map, 360×360×320); probing 800 k atoms takes 1.7 s.
  Once per stage, so this is cheap next to the MD — but `grid=1.0` is 8× smaller
  and ~8× faster for ~0.5% more error on F.
- **Symmetry.** Sampling is periodic in the map's cell, so an atom several cells
  away, or one with unwrapped MD coordinates, reads the equivalent voxel by
  lattice translation (verified exact to the last printed digit for a +5a −3b +2c
  shift). Space-group symmetry is *not* applied at probe time; it is baked into
  the field by `build sg=<name>`, which deposits each atom at all of its
  symmetry-equivalent positions. Doing it at build time keeps every consumer —
  including `sfcalc_gpu_collapse` — a plain trilinear lookup, and avoids
  interpolating across an ASU boundary where the neighbouring voxels belong to a
  different symmetry copy. Building from the 1aho ASU with `sg=P212121` gives a
  field identical to one built from the same model expanded to P1 by hand (max
  difference 6e-5 Å² over 553k voxels, i.e. float32 rounding), and drops the
  fraction of the cell reading `farB` from 57% to 0.22%. Operators come from
  `$CLIBD/symop.lib` (or `symops=`/`symoplib=`); screw axes, `X-Y` forms and
  centring translations are all handled (checked on P3121 and H3).
- **Do not use `sg=` on an MD supercell.** It is already expanded and is
  deliberately not symmetric — expanding it again would average symmetry mates
  together and erase exactly the differences between subcells that the supercell
  exists to capture. `Bfac_map_runme.com` therefore defaults `sg=P1` and does
  *not* inherit `$smallSG` from `xtal_properties.sourceme`.
- **Frame pitfall.** The field must be built on the cell the coordinates live in,
  which for MD is the *supercell*. `nc2mtz` writes per-frame PDBs carrying the
  *primitive* CRYST1 with supercell coordinates, so `probe` needs
  `cell=` or `super_mult=`; both `Bfac_map` and `sfcalc_gpu_collapse` refuse to
  sample when the two frames disagree rather than silently reading wrong voxels.
- Wiring: `Bfac_map=<file>` or `Bfac_map=auto` on `nc2mtz_gemmi.com`,
  `nc2mtz_gpu.com` and `optimize_weights_runme.com`. Empty (default) = per-atom
  path, unchanged. `auto` builds the field into the temp dir from whatever
  per-atom B the old path would have used; pass a file name to keep a field
  across iterations.
- On the GPU path the sampling happens inside `sfcalc_gpu_collapse`
  (`bfacmap=`), so no per-atom sidecar is written at all. On the gemmi path the
  field is sampled into each frame's B column, because `gemmi sfcalc` only takes
  per-atom B — that path exists as the reference implementation.

**Supercell/symmetry**: `expand2supercell_runme.com`, `xsame_runme.com`, `hsame_runme.com`, `wrap_into_cell.com`, `superxform.com`

### `float_add` and `float_func`

Low-level utilities for operating on raw binary float files (including CCP4 electron density maps). They compute statistics (mean, RMS, skewness, kurtosis, CC) and can scale/offset/combine float files. Useful for direct map arithmetic outside of CCP4 tools.

## Demo Systems and Starter Kits

Two crystal systems are used for development and demonstration:

### 1aho — scorpion toxin (reference demo)
64 residues, P212121, 0.965 Å resolution, super_mult=1,1,1 (single unit cell).
No protein ligand; salt = NH4 + ACY at 0.68 M.  All data is public (PDB deposit).
MTZ labels are already FP/SIGFP/FreeR_flag — no cad renaming needed.

Starter kit: `1aho_starter_kit/` in this repo.
Tarball also available: `https://bl831.als.lbl.gov/~jamesh/amber/1aho/example_starter4.tgz`
See `1aho_starter_kit/CLAUDE.md` for full details.

### 6c2r — Aurora Kinase A (production system)
272 residues, P41212, two datasets:
  - AMPPNP (2.4 Å, ligand LIG/AMPPNP, salt NH4+SO4) — projects 6c2r_AMPPNP and 6c2r_LIG
  - EG7 compound (1.9 Å, ligand EG7, salt NH4+SO4) — project 6c2r_37C_2x

Both use super_mult=2,2,2; data from Dirk's private models; badlinks=1 required
(phenix creates a spurious Lys-LIG covalent bond that crashes Amber).

Starter kit: `6c2r_AMPPNP/claude/starter_kit/` in the project directory.
See that kit's `CLAUDE.md` and `PROTOCOL.txt` for full details.

## Example Workflow

See `example_setup_notes.com` for a complete walkthrough starting from a PDB ID,
downloading data, preparing AMBER topology, and running the optimization loop.
The 1aho system in `1aho_starter_kit/` is the recommended starting point.
