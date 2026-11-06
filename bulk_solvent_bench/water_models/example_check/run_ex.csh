#! /bin/tcsh -f
# example protocol (phenix_user_solvent README): starthere.pdb + refme.mtz, 3 cycles xyz+ADP
# usage: run_ex.csh flat|map.mtz tag j [Flabel PHIlabel]   start model starthere_j$j.pdb (xyz jiggled 0.02 A)
set sol = $1
set tag = $2
set j = $3
set fl = FPART
set pl = PHIPART
if( $#argv > 4 ) then
  set fl = $4
  set pl = $5
endif
source /programs/phenix-2.2-6143/phenix_env.csh
set opts = ( refme.mtz refinement.main.number_of_macro_cycles=3 \
  refinement.refine.strategy=individual_sites+individual_adp --overwrite )
if( "$sol" != "flat" ) set opts = ( $opts refinement.input.bulk_solvent_map.file_name=$sol \
  refinement.input.bulk_solvent_map.amplitudes_label=$fl refinement.input.bulk_solvent_map.phases_label=$pl )
phenix.refine starthere_j${j}.pdb $opts output.prefix=ex_${tag}_j${j} >&! ex_${tag}_j${j}.log
grep -E "Final R-work" ex_${tag}_j${j}.log
