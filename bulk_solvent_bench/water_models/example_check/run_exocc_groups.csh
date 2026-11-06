#! /bin/tcsh -f
# run_exocc.csh plus the example's explicit per-residue occupancy groups (starthere_occ_groups.eff)
# usage: run_exocc_groups.csh model_prefix flat|map.mtz tag j
set model = $1
set sol = $2
set tag = $3
set j = $4
source /programs/phenix-2.2-6143/phenix_env.csh
set opts = ( refme.mtz starthere_occ_groups.eff refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp+occupancies --overwrite )
if( "$sol" != "flat" ) set opts = ( $opts refinement.input.bulk_solvent_map.file_name=$sol \
  refinement.input.bulk_solvent_map.amplitudes_label=FPART refinement.input.bulk_solvent_map.phases_label=PHIPART )
phenix.refine ${model}_j${j}.pdb $opts output.prefix=exgrp_${model}_${tag}_j${j} >&! exgrp_${model}_${tag}_j${j}.log
grep -E "Final R-work" exgrp_${model}_${tag}_j${j}.log
