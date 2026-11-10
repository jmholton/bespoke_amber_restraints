#!/bin/bash
# usage: bash highk_chain.sh   (on octamus1, from hurry/)
# the high-weight maps of the 2026-10-31 scan (w3 = 3 x rho cap 30, w10 = 10 x rho cap 99, flat10 = 10 on every heavy
# atom) through the water-addition test here (hurry_chain.csh o_<tag>: example model from 2 jiggled starts, then 3
# altconf rounds)
D=/home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
cd $D || exit 9
tags="o_w3 o_w10 o_flat10"
./hurry_chain.csh $tags
bash rfree_table.sh $tags o_w1 o_n1
echo "highk chain done $(date +%H:%M)"
