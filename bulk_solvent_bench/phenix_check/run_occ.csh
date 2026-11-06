#! /bin/tcsh -f
# deposited waters kept, their occupancies refined; flat mask or an MD bulk solvent map
# usage: run_occ.csh flat|solvent.mtz tag
set sol = $1
set tag = $2
source /programs/phenix-2.2-6143/phenix_env.csh
set opts = ( refme_small.mtz refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp+occupancies \
  refinement.refine.occupancies.individual=water --overwrite )
if( "$sol" != "flat" ) set opts = ( $opts refinement.input.bulk_solvent_map.file_name=$sol \
  refinement.input.bulk_solvent_map.amplitudes_label=FPART refinement.input.bulk_solvent_map.phases_label=PHIPART )
phenix.refine withwater.pdb $opts output.prefix=occ_${tag} >&! occ_${tag}.log
grep -E "Start R-work|Final R-work" occ_${tag}.log
