#! /bin/tcsh -f
# usage: barometer.csh tag   (on octamus1) after the run's map is written: 0.5 ns from prod_2 with every restraint
# off, same cell, same water, same pressure path as production -> barometer.out
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan
set tag = $1
cd $S/$tag
while ( 1 )
  grep -q "^ERROR" run.log
  if( ! $status ) exit 9
  grep -q "^wrote bulk_solvent.mtz" run.log
  if( ! $status ) break
  sleep 60
end
source $S/compute_settings.sourceme >& /dev/null
awk '/ ntr=/{sub(/ntr=1/,"ntr=0")} / nstlim=/{print "  nstlim=250000,";next} / ntwr=/{print "  ntwr=250000,";next}\
     / ntpr=/{print "  ntpr=2500,";next} / ntwx=/{print "  ntwx=0,";next} {print}' prod.in >! barometer.in
$pmemd -O -i barometer.in -o barometer.out -p hydrated.prmtop -c prod_2.rst7 -ref hydrated_ref.crd \
  -r barometer.rst7 -inf barometer.mdinfo
awk '/NSTEP/ && /PRESS/{n++; if(n>20){s+=$NF;m++}} END{printf("barometer (0.1-0.5 ns, unrestrained): P = %.0f bar over %d samples\n",s/m,m)}' barometer.out
