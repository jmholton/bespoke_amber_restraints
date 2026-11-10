#! /bin/tcsh -f
# usage: rounds_scan.csh tag ...   : the three add_waters rounds (rounds.csh then round_heavy.csh) on scan maps.
# Rounds 1-2: occ 0.5, B 15; round 3 adds heavy dist=2.3.  Each round fills the j1 refinement of the round before
# (round 1: exocc_starthere_scan_<tag>_j1_001) and re-refines with the same map from 2 jiggled starts
# -> example_check/exocc_sw<N>_<tag>_scan_<tag>_j*
set A = /home2/jamesh/projects/amber/bulk_solvent/water_models/addwater
set X = /home2/jamesh/projects/amber/bulk_solvent/water_models/example_check
set path = ( /home/jamesh/Develop $path )
set tags = ( $argv )
foreach n ( 1 2 3 )
  set extra = ( )
  if( $n == 3 ) set extra = ( heavy dist=2.3 )
  foreach t ( $tags )
    mkdir -p $A/s${n}_$t ; cd $A/s${n}_$t
    set p = $X/exocc_starthere_scan_${t}_j1_001
    if( $n > 1 ) set p = $X/exocc_sw`expr $n - 1`_${t}_scan_${t}_j1_001
    cp ${p}.pdb model.pdb ; cp ${p}.mtz maps.mtz
    $A/add_waters.com model.pdb maps.mtz conf=_ occ=0.5 B=15 $extra >&! add_waters.log
    echo "$t round ${n}: added `grep -c HOH new_water.pdb` waters"
    foreach j ( 1 2 )
      awk -v seed=`expr $j \* 1111` -v shift=0.02 -f ~/projects/git/bespoke_amber_restraints/jigglepdb.awk new.pdb >! $X/sw${n}_${t}_j${j}.pdb
    end
  end
  cd $X
  foreach t ( $tags )
    foreach j ( 1 2 )
      ./run_exocc.csh sw${n}_$t scan_${t}_map.mtz scan_$t $j >&! exocc_sw${n}_${t}_j${j}.out &
    end
  end
  wait
  echo "round $n done"
end
