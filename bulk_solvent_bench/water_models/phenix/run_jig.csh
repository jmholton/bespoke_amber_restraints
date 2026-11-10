#! /bin/tcsh -f
# usage: run_jig.csh flat|map.mtz tag j   (start model withwater_j$j.pdb, xyz jiggled 0.02 A)
set sol = $1
set tag = $2
set j = $3
source /programs/phenix-2.2-6143/phenix_env.csh
set opts = ( refme_small.mtz refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp+occupancies \
  refinement.refine.occupancies.individual=water --overwrite )
if( "$sol" != "flat" ) set opts = ( $opts refinement.input.bulk_solvent_map.file_name=$sol \
  refinement.input.bulk_solvent_map.amplitudes_label=FPART refinement.input.bulk_solvent_map.phases_label=PHIPART )
phenix.refine withwater_j${j}.pdb $opts output.prefix=jig_${tag}_j${j} >&! jig_${tag}_j${j}.log
grep -E "Final R-work" jig_${tag}_j${j}.log
