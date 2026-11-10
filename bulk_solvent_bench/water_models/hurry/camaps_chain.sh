#!/bin/bash
# usage: bash camaps_chain.sh   (on octamus1, from hurry/)
# The CA-only runs of build_ff19_opc (sets.txt hl_mult, hl_loca, hl_locw, hl_ff19opc) at uniform B 10, then the
# water-addition test on them side by side (hurry_chain.csh: example model from 2 jiggled starts, then 3 altconf rounds)
D=/home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
cd $D || exit 9
tags="hl_mult hl_loca hl_locw hl_ff19opc"
./make_all.csh $tags
for t in $tags; do cat make_map_$t.log; done
echo "maps done $(date +%H:%M)"
./hurry_chain.csh $tags
bash rfree_table.sh $tags
echo "camaps chain done $(date +%H:%M)"
