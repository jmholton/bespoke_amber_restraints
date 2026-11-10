#! /bin/tcsh -f
# example data (refme.mtz), 3 cycles xyz+ADP+water occupancies
# usage: run_exocc.csh model_prefix flat|map.mtz tag j   (start model ${model_prefix}_j$j.pdb)
set model = $1
set sol = $2
set tag = $3
set j = $4
source /programs/phenix-2.2-6143/phenix_env.csh
set opts = ( refme.mtz refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp+occupancies \
  refinement.refine.occupancies.individual=water --overwrite )
if( "$sol" != "flat" ) set opts = ( $opts refinement.input.bulk_solvent_map.file_name=$sol \
  refinement.input.bulk_solvent_map.amplitudes_label=FPART refinement.input.bulk_solvent_map.phases_label=PHIPART )
phenix.refine ${model}_j${j}.pdb $opts output.prefix=exocc_${model}_${tag}_j${j} >&! exocc_${model}_${tag}_j${j}.log
grep -E "Final R-work" exocc_${model}_${tag}_j${j}.log
