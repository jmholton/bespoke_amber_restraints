#! /bin/tcsh -f
# lower-B and own-B re-renders of longcalm and af_std (sets.txt), then the water-addition test (hurry_chain.csh)
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
cd $D
$D/make_all.csh hl_lcb5 hl_lcb2 hl_afb5 hl_afb2 &
foreach t ( hl_lcbf hl_afbf )
  $D/bfac_map.csh $t >& $D/bfac_map_${t}.log &
end
wait
cat $D/make_map_hl_lcb5.log $D/make_map_hl_lcb2.log $D/make_map_hl_afb5.log $D/make_map_hl_afb2.log | grep mtz
cat $D/bfac_map_hl_lcbf.log $D/bfac_map_hl_afbf.log
echo "maps done"
$D/hurry_chain.csh hl_lcb5 hl_lcb2 hl_afb5 hl_afb2 hl_lcbf hl_afbf
