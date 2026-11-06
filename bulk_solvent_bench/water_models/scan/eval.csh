#! /bin/tcsh -f
# usage: eval.csh tag [map.mtz]   (default map: scan/<tag>/bulk_solvent.mtz, waited for)
# refines with the map on bl831: example model (starthere, water occupancies) from 2 jiggled starts, deposited
# model (withwater, refme_small) from 4; then the difference-map blobs of the first example refinement
set W = /home2/jamesh/projects/amber/bulk_solvent/water_models
set tag = $1
set map = $W/scan/$tag/bulk_solvent.mtz
if( $#argv > 1 ) set map = $2
if( $#argv < 2 ) then
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
endif
# no more than 12 of these phenix jobs on bl831 at once
while ( `ps -u $USER -o args | grep -c "^/home/programs/phenix-2.2-6143/bin/python.*command_line"` > 6 )
  sleep 60
end
cp $map $W/example_check/scan_${tag}_map.mtz
cp $map $W/phenix/scan_${tag}_map.mtz
cd $W/example_check
foreach j ( 1 2 )
  ./run_exocc.csh starthere scan_${tag}_map.mtz scan_$tag $j >&! exocc_starthere_scan_${tag}_j${j}.out &
end
cd $W/phenix
foreach j ( 1 2 3 4 )
  ./run_jig.csh scan_${tag}_map.mtz scan_$tag $j >&! jig_scan_${tag}_j${j}.out &
end
wait
cd $W/example_check
source /programs/phenix-2.2-6143/phenix_env.csh
phenix.python $W/diffmaps/diffblobs.py $tag exocc_starthere_scan_${tag}_j1_001.pdb exocc_starthere_scan_${tag}_j1_001.mtz 3 \
  >&! $W/scan/diffblobs_${tag}.txt
echo "$tag done"
