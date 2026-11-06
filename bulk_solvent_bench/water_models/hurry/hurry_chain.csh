#! /bin/tcsh -f
# usage: hurry_chain.csh tag ...   (on octamus1, which has the cores; compare only with runs made on this host)
# the water-addition test for each map: example model from 2 jiggled starts (run_exocc.csh), then the three altconf
# add_waters rounds (addwater/rounds_close.csh), all tags side by side.  Maps: hurry/map_<tag>.mtz for hl_* tags;
# o_w1, o_n1, o_md2019 re-run the scan w1 and n1 maps and the 2019 map here as same-host references.
set W = /home2/jamesh/projects/amber/bulk_solvent/water_models
set X = $W/example_check
set tags = ( $argv )
cd $X
foreach t ( $tags )
  set src = $W/hurry/map_${t}.mtz
  if( $t == o_w1 ) set src = $X/scan_w1_map.mtz
  if( $t == o_n1 ) set src = $X/scan_n1_map.mtz
  if( $t == o_md2019 ) set src = $X/md2019_map.mtz
  while ( ! -e $src )
    sleep 30
  end
  cp $src scan_${t}_map.mtz
  foreach j ( 1 2 )
    ./run_exocc.csh starthere scan_${t}_map.mtz scan_$t $j >&! exocc_starthere_scan_${t}_j${j}.out &
  end
end
wait
echo "starthere done"
$W/addwater/rounds_close.csh $tags
echo "hurry chain done"
