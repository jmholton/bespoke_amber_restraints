#! /bin/tcsh -f
# usage: round.csh N prev_prefix   : add_waters (occ 0.5, B 15) on each map's j1 refinement from round N-1,
#                                    then re-refine with the same map from 2 jiggled starts -> exocc_awN_<map>_<map>_j*
set n = $1
set prev = $2
set A = /home2/jamesh/projects/amber/bulk_solvent/water_models/addwater
set X = /home2/jamesh/projects/amber/bulk_solvent/water_models/example_check
set path = ( /home/jamesh/Develop $path )
foreach t ( flat example md2019 opc_10ns spce_10ns )
  mkdir -p $A/r${n}_$t ; cd $A/r${n}_$t
  set p = ${X}/${prev}_${t}_j1_001
  if( "$prev" =~ exocc_aw* ) set p = ${X}/${prev}_${t}_${t}_j1_001
  cp ${p}.pdb model.pdb ; cp ${p}.mtz maps.mtz
  $A/add_waters.com model.pdb maps.mtz conf=_ occ=0.5 B=15 >&! add_waters.log
  echo "$t round ${n}: added `grep -c HOH new_water.pdb` waters"
  foreach j ( 1 2 )
    awk -v seed=`expr $j \* 1111` -v shift=0.02 -f ~/projects/git/bespoke_amber_restraints/jigglepdb.awk new.pdb >! $X/aw${n}_${t}_j${j}.pdb
  end
end
cd $X
foreach j ( 1 2 )
  ./run_exocc.csh aw${n}_flat flat flat $j >&! exocc_aw${n}_flat_j${j}.out &
  ./run_exocc.csh aw${n}_example example_FPART.mtz example $j >&! exocc_aw${n}_example_j${j}.out &
  foreach m ( md2019 opc_10ns spce_10ns )
    ./run_exocc.csh aw${n}_${m} ${m}_map.mtz $m $j >&! exocc_aw${n}_${m}_j${j}.out &
  end
end
wait
echo "round $n done"
