#! /bin/tcsh -f
# usage: close_chain.csh tag   (on bl831) : water_chain.csh with the altconf rounds: wait for scan/<tag>/bulk_solvent.mtz,
# example refinement from 2 jiggled starts, then addwater/rounds_close.csh (3 altconf add_waters rounds)
set W = /home2/jamesh/projects/amber/bulk_solvent/water_models
set tag = $1
while ( 1 )
  if( -e $W/scan/$tag/run.log ) then
    grep -q "^ERROR" $W/scan/$tag/run.log
    if( ! $status ) then
      echo "$tag run failed"
      exit 9
    endif
    grep -q "^wrote bulk_solvent.mtz" $W/scan/$tag/run.log
    if( ! $status ) break
  endif
  sleep 60
end
sleep 5
cp $W/scan/$tag/bulk_solvent.mtz $W/example_check/scan_${tag}_map.mtz
cd $W/example_check
foreach j ( 1 2 )
  ./run_exocc.csh starthere scan_${tag}_map.mtz scan_$tag $j >&! exocc_starthere_scan_${tag}_j${j}.out &
end
wait
$W/addwater/rounds_close.csh $tag
echo "$tag done"
