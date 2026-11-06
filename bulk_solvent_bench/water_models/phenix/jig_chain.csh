#! /bin/tcsh -f
# 4 jiggled starts for flat and for every map the chains deliver; results in jig_*.out
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/phenix
foreach j ( 1 2 3 4 )
  ./run_jig.csh flat flat $j >&! jig_flat_j${j}.out &
end
set n = 0
while ( $n < 24 )
  foreach m ( `ls -1 *ns.mtz |& grep mtz` )
    set tag = $m:r
    if( -e jig_${tag}_j1.out ) continue
    sleep 20
    foreach j ( 1 2 3 4 )
      ./run_jig.csh $m $tag $j >&! jig_${tag}_j${j}.out &
    end
    @ n = ( $n + 1 )
  end
  sleep 60
end
wait
