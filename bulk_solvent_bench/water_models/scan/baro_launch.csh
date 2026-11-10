#! /bin/tcsh -f
# usage: baro_launch.csh tag   (detaches barometer.csh on this host)
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/scan/$1
nohup ../barometer.csh $1 >& barometer.log < /dev/null &
