#! /bin/tcsh -f
# usage: one_variant.csh tag key=value ...   (on octamus1)  a 2 ns bulk_solvent_runme.com variant of the ensemble_opc
# run (PRISM 7KR0 S0112, opc, 1,3,1) in v_<tag>, then its map refined against the ensemble with PRISM's waters
# (refine/run_ref.csh, seeds 1 and 2)
set B = /home2/jamesh/projects/amber/bulk_solvent/test_7kr0
set tag = $1
shift
mkdir -p $B/v_$tag
cd $B/v_$tag
cp $B/ensemble_opc/7kr0_ens.pdb $B/ensemble_opc/7kr0.mtz $B/ensemble_opc/map_coeffs.mtz .
ln -sf ../../compute_settings.sourceme .
/home2/jamesh/projects/amber/bulk_solvent/bespoke_amber_restraints/bulk_solvent_runme.com 7kr0_ens.pdb 7kr0.mtz \
  mapmtz=map_coeffs.mtz super_mult=1,3,1 watertype=opc $argv >&! run.log
if(! -e bulk_solvent.mtz) exit 9
cp bulk_solvent.mtz $B/refine/md_$tag.mtz
cd $B/refine
./run_ref.csh ens md_$tag.mtz $tag 1 >&! ${tag}_s1.out &
./run_ref.csh ens md_$tag.mtz $tag 2 >&! ${tag}_s2.out &
wait
