#!/bin/bash
# usage: bash alchain.sh   (on octamus1, from hurry/)
# The aligned maps of sets.txt (frames translated onto each run's alignment reference first, as the run does before it
# makes its own mtz): the CA-only runs and two restrained controls, then the water-addition test on them side by side
D=/home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
cd $D || exit 9
tags="hl_mult_al hl_loca_al hl_locw_al hl_ff19opc_al hl_b5sd1_al hl_afstd_al"
./make_all.csh $tags
for t in $tags; do cat make_map_$t.log; done
echo "maps done $(date +%H:%M)"
./hurry_chain.csh $tags
bash rfree_table.sh $tags hl_mult hl_loca hl_locw hl_ff19opc hl_b5sd1 hl_afstd
echo "alchain done $(date +%H:%M)"
