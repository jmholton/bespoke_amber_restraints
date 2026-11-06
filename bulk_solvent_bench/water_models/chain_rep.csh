#! /bin/tcsh -f
# usage: chain_rep.csh rundir   (replicate run, 2 x 1 ns): map from ns 1-2, refined from 4 jiggled starts
set r = $1
set W = /home2/jamesh/projects/amber/bulk_solvent/water_models
while ( ! -e $W/$r/prod_2.rst7 )
  if( -e $W/$r/run.log ) then
    grep -q "ERROR" $W/$r/run.log
    if( ! $status ) then
      echo "$r run failed"
      exit 9
    endif
  endif
  sleep 30
end
sleep 10
$W/submap.csh $W/$r 2 >&! $W/$r/submap_2ns.log
cp $W/$r/maps_2ns/bulk_solvent_2ns.mtz $W/phenix/${r}_2ns_map.mtz
cd $W/phenix
foreach j ( 1 2 3 4 )
  ./run_jig.csh ${r}_2ns_map.mtz ${r}_2ns $j >&! jig_${r}_2ns_j${j}.out &
end
wait
