#!/bin/bash
# usage: run_ref.sh model.pdb flat|map.mtz prefix seed [cycles]   (on octamus1, in a directory holding 7kr0.mtz)
# One phenix.refine step with the 7KR0 six-conformer recipe of
# /home/jamesh/projects/agentic_refine/hierarchical-refinement/docs/7KR0_SETTINGS.md (the harness command at 25ec394):
# fixed weights wc 0.3 / wu 0.25 with both optimisers off, no N/Q/H flips, no ordered_solvent, SA, H addition or
# link_residues, nproc 1 (run seeds side by side), ADP types from the input (ANISOU -> anisotropic), occupancies
# held, 10 macro-cycles (20 for a convergence step).  One change from it, by the user's choice (2026-11-04):
# hydrogens refined individually (explicit), not riding as the recipe has it.  Free set: FreeR_flag 0, the
# deposited 8 712.
# map.mtz: the patched phenix's user-supplied bulk solvent (FPART/PHIPART) in place of the flat mask.
# Writes <prefix>.log, <prefix>_001.pdb/.mtz; prints phenix's Final R line.
model=$1; sol=$2; prefix=$3; seed=$4; cycles=${5:-10}
source /programs/phenix-2.2-6143/phenix_env.sh
# Reduce still runs for the clashscore (start and end statistics) and times out at 200 s on the ensemble once it has
# new full-occupancy waters, killing the run.  scripts/reduce_wrap makes it a 2 s call; with nqh_flips off it only
# changes the reported clashscore.  The phenix dispatchers put their own bin first on PATH unless
# CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES is set, so set it and do what they would have done: clear these, put
# phenix's bin on PATH (here after the wrapper, ahead of any other phenix in the login environment).
export PATH=/home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts/reduce_wrap:/programs/phenix-2.2-6143/bin:$PATH
export CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES=1
unset PYTHONHOME PYTHONPATH LD_LIBRARY_PATH
opts=(7kr0.mtz
  refinement.input.xray_data.labels=FP,SIGFP
  refinement.input.xray_data.r_free_flags.label=FreeR_flag
  data_manager.fmodel.xray_data.r_free_flags.test_flag_value=0
  data_manager.fmodel.xray_data.r_free_flags.generate=False
  data_manager.fmodel.xray_data.r_free_flags.ignore_r_free_flags=False
  refinement.main.number_of_macro_cycles=$cycles
  refinement.main.random_seed=$seed
  refinement.main.nproc=1
  refinement.refine.strategy=individual_sites+individual_adp
  refinement.main.nqh_flips=False
  refinement.hydrogens.refine=individual
  refinement.main.find_and_add_hydrogens=False
  refinement.main.ordered_solvent=False
  refinement.main.simulated_annealing=False
  refinement.target_weights.wc=0.3
  refinement.target_weights.wu=0.25
  refinement.target_weights.optimize_xyz_weight=False
  refinement.target_weights.optimize_adp_weight=False
  refinement.pdb_interpretation.automatic_linking.link_residues=False)
if [ "$sol" != flat ]; then
  opts+=(refinement.input.bulk_solvent_map.file_name=$sol
    refinement.input.bulk_solvent_map.amplitudes_label=FPART
    refinement.input.bulk_solvent_map.phases_label=PHIPART)
fi
phenix.refine "$model" "${opts[@]}" output.prefix=$prefix ${DRYRUN:+--dry-run} --overwrite > $prefix.log 2>&1
grep "Final R-work" $prefix.log
