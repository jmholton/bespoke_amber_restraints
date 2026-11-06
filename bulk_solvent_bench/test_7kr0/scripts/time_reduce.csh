#! /bin/tcsh -f
# usage: time_reduce.csh model.pdb   how long phenix.refine's clashscore Reduce call takes (phenix times it out at 200 s)
source /programs/phenix-2.2-6143/phenix_env.csh
set t0 = `date +%s`
molprobity.reduce -oh -his -flip -keep -allalt -limit120 -pen9999 - < $1 >& /dev/null
set t1 = `date +%s`
echo "$1 `expr $t1 - $t0` s"
