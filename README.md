# bespoke_amber_restraints

Command-line tools for using X-ray crystallography data to iteratively optimize harmonic restraints in AMBER molecular dynamics simulations. The goal is to find the **minimum restraint weights** needed to keep a structure consistent with its electron density — for both Bragg and diffuse scattering data.

Example data for the 1aho system: `https://bl831.als.lbl.gov/~jamesh/amber/1aho/example_starter4.tgz`

---

## Overview

Standard AMBER force fields produce accurate local geometry, but they know nothing about the crystal. Without guidance, an MD simulation of a crystal structure will drift from its electron density. The naive fix — apply strong harmonic restraints to all atoms — overcorrects: it forces the model to stay where the crystallographic refinement put it, defeating the purpose of running MD.

**Bespoke restraints** solve this by finding the minimum per-atom restraint weight that is actually needed. Atoms that are well-served by the AMBER force field converge toward zero weight. Atoms that genuinely need to be pulled back into density (by the data) retain nonzero weights, and those weights encode how strongly the data constrains each atom.

### The feedback loop

Each iteration of the optimization loop:

1. **MD run** — AMBER propagates the structure under force-field + restraints.
2. **Map generation** — Trajectory frames are converted to electron density maps (via structure factor calculation and Fourier transform). The averaged map gives calculated structure factors FC/PHIC. If diffuse scattering is enabled, the same frames also accumulate unphased F² sums (ΣF², the diffuse signal).
3. **Fo-Fc difference map** — Observed amplitudes (Fo) phased with calculated phases give a difference map: positive peaks where the model should move *toward* density, negative peaks where it is being pushed *too hard*.
4. **Weight update** — Each atom's restraint weight is multiplied by a factor proportional to its Fo-Fc peak height. Positive peaks increase weights; negative peaks decrease them. Weights below a cutoff are dropped entirely.
5. **Repeat** — Updated restraints feed back into the next MD run. Weights converge when the structure is consistent with the data.

### Bragg vs. diffuse data

- **Bragg data** (standard MTZ with FP/SIGFP/FreeR_flag) drives the restraint weight update. The iterative loop is designed around Bragg amplitudes.
- **Diffuse scattering** data can be accumulated in parallel. Each MD frame contributes phased structure factors (whose average gives Bragg-like maps) and unphased F² (which accumulates as diffuse intensity). The resulting ΣF² can be compared directly with experimental diffuse data — measuring whether the simulated atomic disorder matches the observed correlated motion in the crystal.

### B factors and restraint weights

Restraint weights in AMBER (kcal/mol/Å²) are stored as B-factor values in the restraint PDB files using the conversion:

```
AMBER_weight = B_column × pdbscale      (default pdbscale = 0.01)
```

A B value of 100 in the restraint file → weight of 1.0 kcal/mol/Å². The optimization adjusts B values; the conversion to AMBER weights happens at each MD run setup.

### B factors by location instead of by atom

Atomic B factors are normally kept per atom, in `Bfac.pdb`, matched to each
trajectory frame by atom order. That has two costs: the matching is fragile
(atoms get stripped, reordered, added), and B can only be relaxed slowly or the
refinement destabilizes — so a water that has just settled onto an ordered site
still carries the high B it earned while it was wandering.

`Bfac_map` stores B as a continuous field `B(x,y,z)` in a CCP4 map instead. Each
atom donates its B to the space around it:

```
B(x) = ( SUM_i w_i(x) B_i + w0 farB ) / ( SUM_i w_i(x) + w0 )
w_i  = exp( -r_i(x)^2 / 2 sigma^2 )      w0 = exp( -fardist^2 / 2 sigma^2 )
```

so `B(x)` is the local average of nearby atomic B, decaying to `farB` (999, i.e.
"nothing is here") in space with no atom within `fardist`. At structure-factor
time each atom takes its B by sampling the field at its own position: B belongs
to the location, so an atom that moves gets the new location's B immediately
while the field itself can be updated as slowly as stability requires.

`sigma` is the one parameter that matters — how far B blends between neighbouring
atoms of different mobility. Measured against the per-atom path (1aho expanded to
P1, `dmin` 1.0), R on F was 0.1% at `sigma` 0.25, 1.2% at 0.5, 3.4% at 1.0 and 15%
at 1.5. The voxel size barely matters by comparison (1.1% at 0.25 Å to 2.4% at
1.5 Å), so a large supercell can use a coarse grid to keep the file small — it
costs 4 bytes per voxel, so `grid=1.0` is 8× smaller than `grid=0.5`. On a bimodal model (sharp
core at B=8, solvent shell at B=90 — the case this is meant to help) `sigma` 0.5
reproduced the per-atom structure factors to 0.01%.

The field must be built from coordinates where every atom is real. `Bfac.pdb` is
a *sidecar* — it is applied to frames by atom order, so its own coordinates may be
placeholders (in one production run, 400 000 of 801 634 atoms sit stacked at the
origin). Building from that file costs R = 10% on F; building the same B values on
a trajectory frame's coordinates costs 0.6%. `Bfac_map=auto` therefore merges the
per-atom B onto a frame first, which also reduces the fragile atom-order matching
from once per frame to once per stage.

Sampling is periodic in the map's cell, so an atom several cells away — or one
with unwrapped MD coordinates — reads the equivalent voxel by lattice
translation. Space-group symmetry is handled when the field is *built*
(`sg=<name>` deposits each atom at all of its symmetry-equivalent positions), not
when it is sampled, so every consumer stays a plain trilinear lookup. Use `sg=`
only when the source model is a single asymmetric unit: an MD supercell is
already expanded and is deliberately not symmetric, so expanding it again would
average symmetry mates together.

Turn it on with `Bfac_map=` (a field you maintain) or `Bfac_map=auto` (built from
`Bfac_file` each time) in `optimize_weights_runme.com`, `nc2mtz_gemmi.com` or
`nc2mtz_gpu.com`. Empty (the default) keeps the per-atom path.

---

## Prerequisites

| Package | Purpose | Key programs |
|---|---|---|
| **AMBER 18+** | MD engine, topology tools | `pmemd.cuda_SPFP`, `cpptraj`, `tleap`, `AddToBox` |
| **CCP4 Suite** | Map/MTZ manipulation | `mtzdump`, `cad`, `mapmask`, `sfall`, `sftools`, `reindex` |
| **Phenix Suite** | Structure refinement, B-factor optimization | `phenix.refine`, `phenix.fetch_pdb`, `phenix.cif_as_mtz` |
| **gemmi** | Fast MTZ/map conversion | `gemmi map2sf` |
| **gnuplot** | Diagnostic plots (optional) | `gnuplot` |
| **SLURM** | Parallel frame processing (optional) | `srun`, `sinfo` |

GPU-accelerated AMBER (`pmemd.cuda_SPFP`) is strongly recommended. The scripts default to SLURM for job dispatch; a local fallback is used automatically when SLURM is not available.

AMBER is sourced from `/programs/amber22/amber.csh` by default. Override with the `pmemd` variable in `settings.sourceme`.

---

## Installation and Compilation

Clone the repository and add it to your `PATH`. Three C programs must be compiled before use:

```bash
gcc -o float_add float_add.c -lm
gcc -o float_func float_func.c -lm
gcc -O3 -o Bfac_map Bfac_map.c -lm
```

These have no dependencies beyond the standard C math library.

**[float_add](docs/float_add.md)** — Add, subtract, scale, and offset raw floating-point flat files (including CCP4 `.map` files). Computes map statistics (mean, RMS, skewness, kurtosis, CC). Used for direct map arithmetic and comparison of simulated vs. observed diffuse intensities.

**[float_func](docs/float_func.md)** — Apply any C math function to one or two floating-point files. Used for more complex map transformations.

---

## Quick Start

This walkthrough starts from a PDB ID and ends with a running optimization loop. See `example_setup_notes.com` for the full annotated script this is based on.

### Step 1 — Download data and extract crystal parameters

```tcsh
set pdbid = 6c2r
getcif.com $pdbid
```

This fetches the PDB file and structure factors. Extract the relevant MTZ columns (FP/SIGFP/FreeR_flag) using CCP4's `cad`:

```tcsh
set F = `mtzdmp ${pdbid}.mtz | awk 'NF>10 && ! /2FOFCWT|FWT|DELF/ && $(NF-1)=="F"{print $NF;exit}'`
cad hklin1 ${pdbid}.mtz hklout refme_small.mtz << EOF
labin file 1 E1=$F E2=SIG$F E3=FreeR_flag
labou file 1 E1=FP E2=SIGFP E3=FreeR_flag
EOF
```

If you are using a supercell (recommended — multiple copies of the unit cell allow more atoms to sample independently), expand and reindex the MTZ:

```tcsh
set super_mult = 2,2,2
# expand symmetry, then reindex for supercell
cad hklin1 refme_small.mtz hklout expanded.mtz << EOF
labin file 1 all
outlime space 1
EOF
cad hklin1 expanded.mtz hklout refme_cell.mtz << EOF
labin file 1 all
symm 1
EOF
reindex hklin refme_cell.mtz hklout refme.mtz << EOF
reindex h2,k2,l2
EOF
```

Adjust the `reindex` command for your chosen `super_mult`.

### Step 2 — Create `xtal_properties.sourceme`

This file is sourced by almost every script in the repo. It must live in (or be symlinked into) the run directory.

```tcsh
cat << EOF >! xtal_properties.sourceme
set reso      = 1.8          # data resolution in Angstroms
set smallSG   = P212121      # space group (Hermann-Mauguin notation)
set super_mult = 2,2,2       # supercell multiplier relative to unit cell
set modulo    = 231          # residues in one protein chain
set ligands   = ( )          # three-letter codes of ligands, e.g. ( ATP MG )
set salt      = ( )          # crystallization salt species
set salt_conc = 0.15         # molar concentration
EOF
```

`smallSG` is the space group of the unit cell (not the supercell, which always has space group P1). `super_mult` can be `1,1,1` to work in the unit cell; larger values give more independent sampling but require more memory and compute.

### Step 3 — Prepare ligand parameters

For each non-standard residue:

```tcsh
mkdir ligands
# generate mol2 and frcmod for each ligand using antechamber / phenix
# place results in ligands/
```

See `example_setup_notes.com` for an automated approach using `phenix.elbow` and `antechamber`.

### Step 4 — Build the AMBER system

Starting from a PDB file that has been through normal crystallographic refinement (with correct protonation, no LINK records that AMBER cannot handle):

```tcsh
# prepare a run directory with required input files
mkdir template_dir
cp refmacout.pdb template_dir/amberme.pdb
cp tleap_stub.in template_dir/
cp ligands/*.mol2 ligands/*.frcmod template_dir/

# build AMBER topology and run initial minimization/equilibration
leap2amber.com template_dir/amberme.pdb \
  stages=Cpu,Min,Cool,Heat,Equi,EquiMin \
  leapfile=tleap_stub.in \
  pdbscale=0.01 \
  equi_ns=0.2 prod_ns=0.5
```

This produces `padded.parm7`, `xtal.prmtop`, `orignames.pdb`, and initial trajectory files. It runs minimization, cooling, heating, and equilibration stages automatically.

### Step 5 — Generate centroid reference points

Centroid reference points are positions in electron density space that atoms are restrained toward. They are typically placed on the electron density from a high-quality map:

```tcsh
mkdir centroids
# place centroids_in_density.pdb in centroids/
# format: PDB ATOM records with the last column = normalized density value (0–1)
```

The script `centroids_nearby_runme.com` then matches each AMBER atom to its nearest reference point when generating the restraint list.

### Step 6 — Run the optimization loop

```tcsh
mkdir run01
cd run01
ln -sf ../xtal_properties.sourceme .
ln -sf ../refme.mtz .

optimize_weights_runme.com
```

The script self-discovers prior restart files and loops indefinitely. To stop it cleanly:

```tcsh
touch exit      # stops after the current iteration completes
```

Key output files to monitor:

| File | Contents |
|---|---|
| `details.log` | Per-iteration R factors, weight statistics, difference map peaks |
| `restraints_for_N.pdb` | Per-atom restraint weights at iteration N (B column × pdbscale = AMBER weight) |
| `amber_N.rst7` | AMBER restart file (structure coordinates at end of run N) |
| `amber_N.nc` | AMBER trajectory (all frames from run N) |
| `avg.map` | Averaged electron density map from run N |
| `cootme.mtz` | MTZ with map coefficients for Coot visualization |

---

## Directory Layout

```
project/
  refme.mtz                        # FP SIGFP FreeR_flag (supercell)
  xtal_properties.sourceme         # crystal parameters
  centroids/
    centroids_in_density.pdb       # reference points with density values
    reference0.mtz                 # best-phased reference MTZ (optional)
  template_dir/
    amberme.pdb                    # starting AMBER-ready model
    tleap_stub.in                  # tleap input template
    orignames.pdb                  # atom name mapping
    reference.mtz                  # reference map MTZ
    *.mol2  *.frcmod               # ligand parameters
  ligands/
    *.cif                          # ligand restraint dictionaries for phenix/refmac
  run01/                           # optimization run directory
    xtal_properties.sourceme       # symlinked from ..
    settings.sourceme              # optional run-time parameter overrides
    details.log                    # per-iteration log
    restraints_for_N.pdb           # per-atom weights at iteration N
    amber_N.rst7 / amber_N.nc      # AMBER restart and trajectory
```

---

## The `settings.sourceme` File

Any variable recognized by `optimize_weights_runme.com` can be placed in `settings.sourceme` in the run directory. This file is re-read at the **start of every iteration**, so changes take effect immediately without restarting. This is the preferred way to tune the optimization while it is running.

Example `settings.sourceme` to ramp up run length and tighten the cutoff after the first 20 iterations:

```tcsh
# increase production run length
set prod_ns = 1.0

# tighten cutoff after weights have stabilized
set cutoff_weight = 0.001

# enable phenix.refine for B-factor optimization
set phenix_itr = 5
```

---

## Parameter Reference

All parameters are set with `key=value` syntax on the command line or in `settings.sourceme`. The tables below list the most important ones. Defaults come from `optimize_weights_runme.com`.

### Crystal and map generation

| Parameter | Default | Description |
|---|---|---|
| `reso` | `1.0` | Data resolution in Å. Set from `xtal_properties.sourceme`. |
| `smallSG` | `P212121` | Space group of the unit cell. |
| `super_mult` | `1,1,1` | Supercell multiplier (e.g., `2,2,1`). Larger = more independent sampling. |
| `render_reso` | `0.95` | Map grid spacing in Å for structure factor calculation. |
| `render_B` | `10` | Overall B factor applied during map calculation. Adjusted automatically. |
| `render_B_adjust` | `0.05` | Step size for automatic `render_B` adjustment per iteration. |
| `Bfac_map` | *(empty)* | CCP4 map of B(x,y,z) to sample per-atom B from instead of `Bfac_file`. A file name uses that field; `auto` builds one from `Bfac_file`; empty keeps the per-atom path. |
| `Bfac_map_sigma` | `0.5` | Å. How far B blends between neighbouring atoms when the field is built. |
| `Bfac_map_grid` | `0.5` | Å. Voxel size of the field. |
| `Bfac_map_farB` | `999` | B assigned to space with no atom within `fardist` (4.5 Å). |
| `avglast` | `1` | Number of previous MD runs to average electron density maps over. Increase to reduce noise. |
| `fft_B` | `0` | Apply this B factor (in Å²) to smooth the difference map before weight update. |
| `shan_B` | `auto` | Shannon-entropy smoothing B for the difference map. |

### MD simulation

| Parameter | Default | Description |
|---|---|---|
| `prod_ns` | `0.5` | Production MD run length in nanoseconds. |
| `equi_ns` | `0.2` | Equilibration run length in nanoseconds. |
| `dt` | `0.002` | MD time step in ps. |
| `temperature` | `287` | Simulation temperature in K. |
| `gamma_ln` | `1.0` | Langevin thermostat coupling (ps⁻¹). |
| `barostat` | `0` | Enable NPT barostat (0 = NVT, 1 = isotropic NPT). |
| `netfrc` | `0` | Net force correction (0 = enabled, prevents center-of-mass drift). **Leave at 0.** |
| `write_ps` | `20` | Write trajectory frames every N ps. |

### Restraint weights

| Parameter | Default | Description |
|---|---|---|
| `pdbscale` | `0.01` | Converts PDB B column to AMBER weight (kcal/mol/Å²). |
| `weight0` | `1` | Starting AMBER weight for newly created restraints. |
| `cutoff_weight` | `0.01` | Drop restraints whose weight falls below this (in AMBER units). |
| `cutoff_forget` | `0` | If 1, permanently remove reference points whose weight drops to cutoff. |
| `max_weight` | `9.9999` | Cap on any single restraint weight. |
| `allatom_weight` | `0` | Apply this weight to every atom, even those without explicit restraints. |
| `min_CA_weight` | `0` | Minimum weight for Cα atoms used as alignment references. |
| `min_lig_weight` | `0.01` | Minimum weight for ligand atoms. |

### Weight update behavior

| Parameter | Default | Description |
|---|---|---|
| `adjust_itr` | `1` | Update weights every N iterations. |
| `max_mult` | `2.0` | Maximum ratio by which any weight can increase in one update. |
| `weight_scaledown` | `1` | Multiply all weights by this factor each round (< 1 = gradual global decay). |
| `weight_negscaledown` | `0.95` | Extra scale-down for weights with negative Fo-Fc peaks (atoms pushed too hard). |
| `weight_power` | `1.01` | Raise small weights to this power each round (drives near-zero weights to zero). |
| `weight_power_kTmult` | `1` | Define "small" as ≤ N × kT for the `weight_power` operation. |
| `halfrho_pos` | `auto` | Sigma threshold for a significant positive Fo-Fc peak. Weights below threshold are not increased. |
| `halfrho_neg` | `auto` | Sigma threshold for a significant negative Fo-Fc peak. |
| `ambig_same_weight` | `1` | Force symmetry-equivalent atoms (e.g., OD1/OD2 on ASP) to share the same weight. |

### Diffuse scattering

Diffuse scattering accumulation is controlled within `nc2mtz_gemmi.com`. To enable it, add to `settings.sourceme` before the first iteration (or set it in the `nc2mtz_gemmi.com` call inside the master script):

```tcsh
set addmtzs = 1    # accumulate per-frame MTZ files for diffuse intensity
```

With `addmtzs=1`, after processing each trajectory frame into an MTZ file, `addup_mtzs_diffuse.com` is called to merge all frame MTZs into a single `sum.mtz` containing:
- `Fsum` / `PHIsum` — phased structure factor sum (Bragg signal, same as what drives the difference map)
- `Isum` — unphased F² sum (= ΣF², diffuse scattering signal)

The `Isum` column can be compared to experimental diffuse data using `float_add` to compute a correlation coefficient:

```bash
float_add -header 1104 simulated_diffuse.map experimental_diffuse.map
```

### Periodic maintenance operations

These operations run every N iterations (0 = disabled).

| Parameter | Default | Description |
|---|---|---|
| `repick_itr` | `10` | Re-discover nearest reference points for each atom. Use after significant structural changes. |
| `release_itr` | `50` | Release (reduce weight of) the most persistently challenged restraints. |
| `release_maxbad` | `10` | Maximum number of restraints to release per `release_itr` cycle. |
| `delete_badrest_itr` | `0` | Permanently delete the most persistently challenged restraints. |
| `randel_itr` | `0` | Randomly delete a fraction (`randel_fraction`) of restraints each cycle. |
| `hydrate_itr` | `1` | Check for voids (add waters) and over-pressure (remove waters). |
| `teleport_itr` | `1` | Move waters from negative to positive Fo-Fc density regions. |
| `teleport_waters` | `10` | Maximum number of waters to teleport per cycle. |
| `wrap_itr` | `2` | Wrap non-restrained atoms back into the supercell. |
| `align_itr` | `1` | Re-center the system on restrained atoms. |
| `Badjust_itr` | `1` | Adjust per-atom B factors from the difference map. |
| `phenix_itr` | `0` | Run `phenix.refine` to optimize B factors every N iterations. |
| `refmac_itr` | `0` | Run `refmac` to optimize B factors every N iterations. |
| `remap_itr` | `0` | Re-map water molecule names to match topology. |
| `filter_itr` | `1` | Remove water restraints that are too close to protein atoms. |
| `filter_dist` | `1.8` | Minimum distance (Å) for water restraints from protein. |

### Ramp parameters

Many parameters accept a `_ramp` variant to linearly interpolate their value over iterations. The format is `start end [nstep]`. For example:

```tcsh
set prod_ns_ramp     = "0.5 2.0"    # increase from 0.5 to 2.0 ns over maxitr iterations
set cutoff_weight_ramp = "0.1 0.01" # tighten cutoff as convergence improves
set max_mult_ramp    = "4.0 2.0"    # start aggressive, become conservative
```

Available ramps: `prod_ns_ramp`, `equi_ns_ramp`, `cutoff_weight_ramp`, `max_mult_ramp`, `avglast_ramp`, `pressure_avglast_ramp`, `pressure_scale_ramp`, `repick_maxdist_ramp`, `repick_hohscale_ramp`, `fft_B_ramp`, `halfrho_ramp`, `thrubond_avg_weight_ramp`, `thrubond_avg_B_ramp`, `temp_ramp`, `chiral_weight_ramp`, `omega_weight_ramp`.

### Hydration / pressure

| Parameter | Default | Description |
|---|---|---|
| `hydrate` | `voids` | Criterion for adding waters: `voids` (add near density peaks), `pressure` (add when pressure low). |
| `dehydrate` | `pressure` | Criterion for removing waters: `pressure` (remove when pressure high). |
| `minvoid` | `60` | Minimum void volume (Å³) to trigger water addition. |
| `pressure_deadband` | `10` | Pressure tolerance (atm) before hydration/dehydration acts. |
| `add_radius` | `2.4` | Minimum distance (Å) from existing atoms for new water placement. |
| `water_lock` | `0` | If 1, keep total water count constant (add and remove in equal numbers). |
| `barometer_cycles` | `10000` | MD steps for pressure measurement run (all restraints removed). |

### Chiral centers and peptide bonds

| Parameter | Default | Description |
|---|---|---|
| `chiral_weight` | `0` | AMBER restraint weight on flippable chiral centers (e.g., 10 to prevent chirality flips). |
| `omega_weight` | `0` | AMBER restraint weight on peptide bond dihedrals (e.g., 50 to enforce planarity). |

---

## Script Reference

All scripts use the `.com` extension (tcsh) and accept `key=value` command-line arguments. AWK utilities use the `.awk` extension.

### Master script

**`optimize_weights_runme.com`** — Main entry point. Reads `xtal_properties.sourceme` and optionally `settings.sourceme`, then runs the MD → map → difference map → weight update loop. Self-discovers restart files so it can resume a previous run. Stop by creating an `exit` file in the run directory.

**`example_setup_notes.com`** — Annotated walkthrough for setting up a new system from a PDB ID. Not meant to be run as-is; read it as a recipe.

**`leap2amber.com`** — Runs the full AMBER setup pipeline: tleap, minimization, cooling, heating, equilibration. Produces `padded.parm7`, `xtal.prmtop`, `orignames.pdb`, and initial trajectory stages.

### Restraint management

| Script | Purpose |
|---|---|
| `centroids_nearby_runme.com` | Match AMBER atoms to nearest reference centroid points; generate initial restraint PDB. Key args: `reffile`, `maxdist`, `softener`, `hohscale`. |
| `restraintlist2amber.com` | Convert restraint PDB (B column = weight × pdbscale) to AMBER group restraint format (`.in` file). |
| `restraintlist_update_diffmap.com` | Update restraint weights based on Fo-Fc difference map peaks at each reference point. Core of the feedback loop. Key args: `pdbfile`, `diffmap`, `mtzfile`, `max_mult`, `halfrho_pos`, `halfrho_neg`. |
| `restraint_bomb_detector.com` | Detect and defuse "restraint bombs" — cases where a restraint reference point is unreachably far from the restrained atom, which can cause trajectory explosions. |
| `delete_worst_restraints.com` | Identify and permanently delete restraints that are historically most challenged (high restraint energy). |
| `release_worst_restraints_runme.com` | Temporarily reduce the weight of the most challenged restraints to let atoms relax. |
| `no_new_nonbonds_runme.com` | Prevent restraints from being assigned to atoms that have non-bonded clashes with symmetry mates. |

### Map and MTZ handling

| Script | Purpose |
|---|---|
| `nc2mtz_gemmi.com` | Convert AMBER trajectory (`.nc`) to electron density maps and/or MTZ files. Calls `addup_maps_runme.com` for Bragg maps and `addup_mtzs_diffuse.com` for diffuse accumulation. Key args: `ncfile`, `smallSG`, `reso`, `B`, `addmaps`, `addmtzs`. |
| `addup_maps_runme.com` | Average a set of CCP4 electron density maps. Parallel-capable. |
| `addup_mtzs_runme.com` | Average a set of MTZ files (phased, Bragg). |
| `addup_mtzs_diffuse.com` | Pairwise-parallel merge of frame MTZ files, accumulating both phased `Fsum`/`PHIsum` (Bragg) and unphased `Isum` = ΣF² (diffuse). |
| `map_scaleB_runme.com` | Scale a map by applying an overall B factor correction. |
| `map_func.com` | Apply a mathematical function to map values (wrapper around `float_func`). |
| `nc2mtz_runme.com` | Older map generation script using `cpptraj` + `sfall`. Use `nc2mtz_gemmi.com` for new systems. |

### B-factor refinement

| Script | Purpose |
|---|---|
| `Bfac_update_diffmap.com` | Update per-atom B factors based on the Fo-Fc difference map. Positive peaks → decrease B (atom scatters too little); negative peaks → increase B. |
| `scaleB_search_diffmap_runme.com` | Grid search for the overall B factor that best matches the observed map. |
| `thrubond_avgB_runme.com` | Smooth B factors by averaging through covalent bonds to neighboring atoms. Prevents isolated outlier B values. |
| `rmsd2B` | Compute per-atom B factors from RMSD variation across trajectory frames (compiled C program). |
| `Bfac_map` | Store B factors as a spatial field: `build` a CCP4 map of B(x,y,z) from a PDB, `probe` a map to set per-atom B, `stats` to see the distribution (compiled C program). |
| `Bfac_map_runme.com` | Build `Bfac.map` from `Bfac.pdb` with pipeline defaults, then report how faithfully the field reproduces the B factors that went into it. |

### Water and hydration

| Script | Purpose |
|---|---|
| `hydrate_runme.com` | Add water molecules near electron density peaks (voids in the model). Updates topology, restart, and orignames. |
| `dehydrate_amber_runme.com` | Remove water molecules, typically in response to high pressure. |
| `water_teleport_runme.com` | Move water molecules from negative Fo-Fc density to positive Fo-Fc density. |
| `remap_waters_runme.com` | Re-match water atom names and numbering after structural changes. |
| `add_waters.com` | Low-level water addition using `AddToBox`. |
| `reorganize_waters.com` | Renumber and reorganize water molecules in a PDB file. |
| `pinchedwater_runme.com` | Detect water molecules pinched between symmetry mates. |

### Geometry analysis

| Script | Purpose |
|---|---|
| `quick_geo_check_runme.com` | Fast geometry check: bond lengths, angles, Ramachandran, rotamers. |
| `molprobify_runme.com` | Run MolProbity geometry analysis on the current structure. |
| `distanceify_runme.com` | Compute pairwise distances between atoms or groups. |
| `measure_voids_runme.com` | Measure solvent void volumes in the crystal packing. |
| `phipsichi.com` | Compute φ/ψ/χ dihedral angles for all residues. |

### Supercell and symmetry

| Script | Purpose |
|---|---|
| `expand2supercell_runme.com` | Expand an ASU PDB + MTZ into a supercell (P1) system. Requires `xtal_properties.sourceme`. |
| `xsame_runme.com` | Generate crystallographic symmetry copies within a cutoff distance. Used for restraint generation across symmetry mates. |
| `hsame_runme.com` | Generate hydrogen-bond symmetry mate contacts. |
| `wrap_into_cell.com` | Wrap atoms into the primary supercell. |
| `superxform.com` | Apply a supercell transformation to atom coordinates. |

### Structure utilities

| Script | Purpose |
|---|---|
| `rst2pdb_runme.com` | Convert AMBER restart (`.rst7`) to PDB format using `cpptraj`. Tries available topology files. |
| `update_centroid_positions_runme.com` | Update the reference coordinate file (`ref.crd`) from a current restraint PDB and structure. |
| `combine_pdbs_runme.com` | Merge two PDB files, matching atoms by name/residue and filling in from the reference. |
| `reorganize_pdb_runme.com` | Renumber and reorganize a PDB file. |
| `flip_to_target_runme.com` | Flip ambiguous sidechain conformations (e.g., ASN/GLN amide) to best match a target. |
| `align_to_reference_runme.com` | Align a structure to a reference using selected atoms. |
| `graft_atoms_runme.com` | Graft atoms from one PDB onto another. |
| `buildout_pdb_runme.com` | Build missing atoms and sidechains into a PDB. |
| `add_salt_runme.com` | Add salt ions to the simulation cell. |

### PDB manipulation AWK scripts

| Script | Purpose |
|---|---|
| `convert_pdb.awk` | General-purpose PDB format conversion and atom selection. |
| `filter_pdb.awk` | Select atoms by type (protein, ligand, water, H), residue, chain, or element. |
| `reformatpdb.awk` | Reformat PDB records; set B factors, occupancy, or other columns. |
| `sequence.awk` | Extract one-letter sequence from SEQRES or ATOM records. |
| `jigglepdb.awk` | Add small random displacements to PDB coordinates. |
| `build_c2n.awk` | Build C-terminal to N-terminal peptide backbone. |
| `build_n2c.awk` | Build N-terminal to C-terminal peptide backbone. |
| `build_side.awk` | Build sidechain atoms from backbone geometry. |
| `waterconf.awk` | Generate water molecule orientation from oxygen position. |

---

## Troubleshooting

**"WARNING: no xtal_properties.sourceme"**
Create this file in the run directory (see Step 2 of Quick Start). The script will use hardcoded 1aho defaults and prompt for confirmation if it cannot find the file.

**"leap2amber failed"**
- Check `tleap_stub.in` for correct file paths.
- Verify all `.mol2` and `.frcmod` files are present for every non-standard residue.
- Check `protonation.txt` — it must list the correct protonation state for every histidine, cysteine, and any other titratable residue.
- Inspect `leap2amber_0.log` for the specific tleap error.

**"unfixable restraint bombs detected"**
A restraint reference point is too far from its matched atom. This can happen after large structural changes. Run `restraint_bomb_detector.com` manually to inspect which atoms are affected. You may need to delete the relevant lines from `current_restraints.pdb` and let `centroids_nearby_runme.com` re-assign them.

**Map is empty or all zeros in the Fo-Fc map**
- Check `refme.mtz` column labels: they must be exactly `FP`, `SIGFP`, `FreeR_flag`. Re-run the `cad` step if they differ.
- Check that `reso` in `xtal_properties.sourceme` matches the actual data resolution.
- Check that `smallSG` matches the space group in `refme.mtz`.

**All restraint weights collapse to cutoff**
- Increase `weight0` (starting weight for new restraints).
- Decrease `cutoff_weight` to allow weights to survive at lower values.
- Inspect `details.log` for difference map sigma levels — if Fo-Fc peaks are tiny, the map may be over-smoothed (decrease `fft_B`) or the overall B may be wrong (check `render_B`).
- Decrease `weight_negscaledown` (closer to 1.0) to slow the decay.

**Another job is already running (script waits without starting)**
The script detects another instance running in the same directory and waits for it to exit. To stop the other job:
```tcsh
touch exit      # in the run directory
```
Or find and kill the blocking process manually.

**GPU not found / pmemd fails**
Override the `pmemd` variable in `settings.sourceme`:
```tcsh
set pmemd = "pmemd.cuda_SPFP"     # no SLURM, direct GPU
set pmemd = "sander"              # CPU-only fallback (slow)
```

**Water count drifts unexpectedly**
Set `water_lock=1` to enforce a constant water count. Adjust `pressure_deadband` to reduce over-sensitivity of the hydration/dehydration logic.

---

## Author

<ADDRESS><A HREF="mailto:JMHolton@lbl.gov">James Holton &lt;JMHolton@lbl.gov&gt;</A></ADDRESS>
