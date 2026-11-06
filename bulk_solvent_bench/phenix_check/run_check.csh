#! /bin/tcsh -f
# flat mask vs MD bulk solvent (FPART/PHIPART), with and without the deposited waters
# usage: run_check.csh solventmtz tag
set solmtz = $1
set tag = $2
source /programs/phenix-2.2-6143/phenix_env.csh
set common = ( refme_small.mtz refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp --overwrite )
foreach model ( withwater nowater )
  if( "$tag" == "flat" ) then
    phenix.refine ${model}.pdb $common output.prefix=${model}_flat >&! ${model}_flat.log &
  else
    phenix.refine ${model}.pdb $common output.prefix=${model}_${tag} \
      refinement.input.bulk_solvent_map.file_name=$solmtz \
      refinement.input.bulk_solvent_map.amplitudes_label=FPART \
      refinement.input.bulk_solvent_map.phases_label=PHIPART >&! ${model}_${tag}.log &
  endif
end
wait
