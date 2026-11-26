#! /bin/tcsh -f
# One 0.68 M ammonium acetate placement for the 2026-11-25 charge tests, made exactly as the 1aho example does it
# (add_salt_runme.com, reorganize_pdb_runme.com, filter_pdb.awk), from salt02's refined.pdb (the model every 2s2r build
# starts from) at the carboxylate cell charge +4.  Arms with more positive charge add counter-acetates on top
# (setup_s68_builds.sh), so every arm shares this salt.   usage: salt068_base.csh
set H = /scratch/jamesh/projects/amber/1aho_hurry_test
set d = $H/salt068_base
source /home2/jamesh/projects/amber/bulk_solvent/programs.sourceme
set path = ( $H/code_noh $path )
mkdir -p $d ; cd $d
foreach f ( refined.pdb NH4.pdb NH4.mol2 NH4.frcmod NH4.cif ACY.pdb ACY.mol2 ACY.frcmod ACY.cif )
  cp -p ../salt02/$f .
end
ln -sf NH4.pdb cation.pdb ; ln -sf ACY.pdb anion.pdb
cp ../xtal_properties.sourceme .
echo "add_salt `date`"
add_salt_runme.com refined.pdb conc=0.68 RIP=4 RIW=3 charge=4 anion=ACY cation=NH4 >&! add_salt.log
if (! -e salty.pdb) then
  echo "add_salt FAILED"; tail -5 add_salt.log; exit 9
endif
grep "new charge should be" add_salt.log
echo "reorganize `date`"
reorganize_pdb_runme.com salty.pdb ignore_zero=0 refpdb=refined.pdb \
  outfile=reorganized.pdb phenix_bumpcheck=0 declash=1 >&! reorganize_final.log
if (! -e reorganized.pdb) then
  echo "reorganize FAILED"; tail -5 reorganize_final.log; exit 9
endif
filter_pdb.awk -v skip=water reorganized.pdb | egrep -v "^END" >! amberme.pdb
filter_pdb.awk -v skip=H -v only=water,atoms reorganized.pdb >> amberme.pdb
echo "done `date`"
