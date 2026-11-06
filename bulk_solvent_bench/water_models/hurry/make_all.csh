#! /bin/tcsh -f
# usage: make_all.csh tag ...   (on octamus1)  make_map.csh for each tag, all at once
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
foreach t ( $argv )
  $D/make_map.csh $t >& $D/make_map_${t}.log &
end
wait
echo "maps done"
