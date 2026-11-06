#! /bin/tcsh -f
# usage: track.csh dir sol first   (on octamus1)  resume one bulk solvent's rounds of a rebuild chain (chain.csh,
# chain_iso.csh) from round <first>: wait for the round before it to finish, refine it again if it left no model,
# then build (ens_waters.py, unless w<n>_<sol>.pdb is already there) and refine each round to 3 (seed 2 in round 3)
set S = /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts
set dir = $1
set sol = $2
set first = $3
cd $dir
set map = flat
if( $sol == md ) set map = md_opc.mtz
set n = $first
while ( $n <= 3 )
  set p = `expr $n - 1`
  set prev = ref_w${p}_${sol}_${sol}_s1
  # the round before: wait (at most 3 h) for its log to end, then refine it again if it left no model
  set waited = 0
  while ( -e ${prev}.log && `egrep -c "Job complete|^Sorry|Traceback" ${prev}.log` == 0 && $waited < 180 )
    sleep 60
    @ waited++
  end
  if(! -e ${prev}_001.pdb) then
    echo "$sol round ${p}: no model from $prev, refining it again"
    ./run_ref.csh w${p}_${sol} $map $sol 1 >&! w${p}_${sol}_s1.out
    if(! -e ${prev}_001.pdb) then
      echo "$sol round ${p}: refinement failed again (see ${prev}.log)"
      exit 9
    endif
  endif
  if(! -e w${n}_${sol}.pdb) then
    /programs/phenix-2.2-6143/bin/phenix.python $S/ens_waters.py ${prev}_001.pdb ${prev}_001.mtz w${n}_${sol}.pdb \
      sigma=3.0 >&! build_w${n}_${sol}.log
    echo "$sol round ${n}: `tail -2 build_w${n}_${sol}.log | tr '\n' ' '`"
  endif
  ./run_ref.csh w${n}_${sol} $map $sol 1 >&! w${n}_${sol}_s1.out &
  if( $n == 3 ) then
    ./run_ref.csh w${n}_${sol} $map $sol 2 >&! w${n}_${sol}_s2.out &
  endif
  wait
  echo "$sol round $n done: `cat w${n}_${sol}_s1.out`"
  @ n++
end
