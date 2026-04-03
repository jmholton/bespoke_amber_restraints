# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`bespoke_amber_restraints` is a collection of command-line tools for using X-ray crystallography data to iteratively optimize harmonic restraints in AMBER molecular dynamics simulations. The goal is to find the minimum restraint weights needed to keep a structure consistent with its electron density.

## Compiling the C Programs

Two C programs are included and must be compiled before use:

```bash
gcc -o float_add float_add.c -lm
gcc -o float_func float_func.c -lm
```

These have no dependencies beyond the standard C math library.

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

**B factor refinement**: `Bfac_update_diffmap.com`, `map_scaleB_runme.com`, `scaleB_search_diffmap_runme.com`, `thrubond_avgB_runme.com`, `rmsd2B`

**Supercell/symmetry**: `expand2supercell_runme.com`, `xsame_runme.com`, `hsame_runme.com`, `wrap_into_cell.com`, `superxform.com`

### `float_add` and `float_func`

Low-level utilities for operating on raw binary float files (including CCP4 electron density maps). They compute statistics (mean, RMS, skewness, kurtosis, CC) and can scale/offset/combine float files. Useful for direct map arithmetic outside of CCP4 tools.

## Example Workflow

See `example_setup_notes.com` for a complete walkthrough starting from a PDB ID, downloading data, preparing AMBER topology, and running the optimization loop. Example data for the 1aho system is available at: `https://bl831.als.lbl.gov/~jamesh/amber/1aho/example_starter4.tgz`
