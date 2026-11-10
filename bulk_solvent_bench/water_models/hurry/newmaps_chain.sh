#!/bin/bash
# usage: bash newmaps_chain.sh   (on octamus1, from hurry/)
# The 2026-11-05 maps of sets.txt (the overnight ff19/OPC runs at uniform B 10 and at their own rendering), af_std's
# B 10 map again as a same-host repeat (hl_afstd_rep), then the water-addition test on all of them side by side
# (hurry_chain.csh: example model from 2 jiggled starts, then 3 altconf rounds).
D=/home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
cd $D || exit 9
cp map_hl_afstd.mtz map_hl_afstd_rep.mtz
./make_all.csh hl_b5 hl_b5sd1 hl_k1B hl_b5b5 hl_b5sd1b5 &
./bfac_map.csh hl_k1Bbf > bfac_map_hl_k1Bbf.log 2>&1 &
wait
cat make_map_hl_b5.log make_map_hl_b5sd1.log make_map_hl_k1B.log make_map_hl_b5b5.log make_map_hl_b5sd1b5.log
cat bfac_map_hl_k1Bbf.log
echo "maps done $(date +%H:%M)"
./hurry_chain.csh hl_b5 hl_b5sd1 hl_k1B hl_b5b5 hl_b5sd1b5 hl_k1Bbf hl_afstd_rep
echo "newmaps chain done $(date +%H:%M)"
