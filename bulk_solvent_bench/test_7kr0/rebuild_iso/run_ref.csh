#! /bin/tcsh -f
# usage: run_ref.csh model_prefix flat|map.mtz tag seed   (on octamus1)
# PRISM's standing recipe for 7KR0 S0112: wc 0.3, wu 0.25, no occupancy refinement (held waters), input ADP types
# kept, riding H; 5 macro cycles; free set = FreeR_flag 0 (the deposition's)
set model = $1
set sol = $2
set tag = $3
set seed = $4
source /programs/phenix-2.2-6143/phenix_env.csh
# phenix's clashscore Reduce times out (200 s) on the ensemble with altloc waters: scripts/reduce_wrap drops -flip -allalt
set path = ( /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts/reduce_wrap $path )
# (the phenix dispatchers put their own bin first on PATH unless this is set; then clear what they would have cleared)
setenv CCTBX_CONDA_USE_ENVIRONMENT_VARIABLES 1
unsetenv PYTHONHOME PYTHONPATH LD_LIBRARY_PATH
set opts = ( 7kr0.mtz refinement.main.number_of_macro_cycles=5 refinement.main.random_seed=$seed \
  refinement.refine.strategy=individual_sites+individual_adp \
  refinement.target_weights.wxc_scale=0.3 refinement.target_weights.wxu_scale=0.25 \
  xray_data.r_free_flags.label=FreeR_flag xray_data.r_free_flags.test_flag_value=0 \
  refinement.main.nproc=4 --overwrite )
if( "$sol" != "flat" ) set opts = ( $opts refinement.input.bulk_solvent_map.file_name=$sol \
  refinement.input.bulk_solvent_map.amplitudes_label=FPART refinement.input.bulk_solvent_map.phases_label=PHIPART )
phenix.refine ${model}.pdb $opts output.prefix=ref_${model}_${tag}_s${seed} >&! ref_${model}_${tag}_s${seed}.log
grep -E "Final R-work" ref_${model}_${tag}_s${seed}.log
