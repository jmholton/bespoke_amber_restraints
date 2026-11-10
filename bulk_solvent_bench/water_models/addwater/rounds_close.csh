#! /bin/tcsh -f
# usage: rounds_close.csh tag ...   : rounds_scan.csh with waters allowed close to split-conformer atoms
# (add_waters.com altconf: dist 2.3 A from single-conformer heavy atoms only, mindist 1 A from any atom, each new
# water in a conformer letter no atom within 2.3 A has), all three rounds.  tag md2019 = the 2019 map; else scan maps.
# Round N fills the j1 refinement of round N-1 (round 1: the starthere refinement with that map) and re-refines
# from 2 jiggled starts -> example_check/exocc_cw<N>_<tag>_<maptag>_j*
set A = /home2/jamesh/projects/amber/bulk_solvent/water_models/addwater
set X = /home2/jamesh/projects/amber/bulk_solvent/water_models/example_check
set path = ( /home/jamesh/Develop $path )
set tags = ( $argv )
foreach n ( 1 2 3 )
  foreach t ( $tags )
    set mt = scan_$t
    set map = scan_${t}_map.mtz
    if( $t == md2019 ) then
      set mt = md2019
      set map = md2019_map.mtz
    endif
    mkdir -p $A/c${n}_$t ; cd $A/c${n}_$t
    set p = $X/exocc_starthere_${mt}_j1_001
    if( $n > 1 ) set p = $X/exocc_cw`expr $n - 1`_${t}_${mt}_j1_001
    cp ${p}.pdb model.pdb ; cp ${p}.mtz maps.mtz
    $A/add_waters.com model.pdb maps.mtz occ=0.5 B=15 heavy dist=2.3 altconf >&! add_waters.log
    echo "$t round ${n}: `grep altconf: add_waters.log`"
    foreach j ( 1 2 )
      awk -v seed=`expr $j \* 1111` -v shift=0.02 -f ~/projects/git/bespoke_amber_restraints/jigglepdb.awk new.pdb >! $X/cw${n}_${t}_j${j}.pdb
    end
  end
  cd $X
  foreach t ( $tags )
    set mt = scan_$t
    set map = scan_${t}_map.mtz
    if( $t == md2019 ) then
      set mt = md2019
      set map = md2019_map.mtz
    endif
    foreach j ( 1 2 )
      ./run_exocc.csh cw${n}_$t $map $mt $j >&! exocc_cw${n}_${t}_j${j}.out &
    end
  end
  wait
  echo "round $n done"
end
