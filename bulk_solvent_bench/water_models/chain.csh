#! /bin/tcsh -f
# usage: chain.csh watertype
# as segments finish: maps from ns 1-2, 1-10 and 9-10, each refined (deposited waters, occupancies refined)
set wt = $1
set W = /home2/jamesh/projects/amber/bulk_solvent/water_models
foreach job ( 2:1 10:1 10:9 )
  set n = `echo $job | cut -d: -f1`
  set f = `echo $job | cut -d: -f2`
  set tag = ${n}ns
  if( $f != 1 ) set tag = ${f}-${n}ns
  while ( ! -e $W/$wt/prod_${n}.rst7 )
    if( -e $W/$wt/run.log ) then
      grep -q "ERROR" $W/$wt/run.log
      if( ! $status ) then
        echo "$wt run failed"
        exit 9
      endif
    endif
    sleep 30
  end
  sleep 10
  $W/submap.csh $W/$wt $n $f >&! $W/$wt/submap_${tag}.log
  cp $W/$wt/maps_${tag}/bulk_solvent_${tag}.mtz $W/phenix/${wt}_${tag}.mtz
  ( cd $W/phenix ; ./run_occ.csh ${wt}_${tag}.mtz ${wt}_${tag} >&! occ_${wt}_${tag}.out & )
end
