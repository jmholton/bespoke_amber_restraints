#! /bin/tcsh -f
# usage: time_reduce2.csh model.pdb   time_reduce.csh through scripts/reduce_wrap
source /programs/phenix-2.2-6143/phenix_env.csh
set path = ( /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts/reduce_wrap $path )
which molprobity.reduce
set t0 = `date +%s`
molprobity.reduce -oh -his -flip -keep -allalt -limit120 -pen9999 - < $1 |& grep -c "^ATOM"
set t1 = `date +%s`
echo "$1 `expr $t1 - $t0` s"
