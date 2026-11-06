#!/bin/bash
# usage: bash rfree_table.sh tag ...
# R-free of the water-addition test for each map, mean of the 2 jiggled starts (each in brackets): before water
# building (starthere refinement) and after altconf rounds 1-3, from the "Final R-work" line of each refinement log
X=/home2/jamesh/projects/amber/bulk_solvent/water_models/example_check
printf "%-12s %-22s %-22s %-22s %-22s\n" tag none r1 r2 r3
for t in "$@"; do
  line=$(printf "%-12s" $t)
  for s in starthere cw1 cw2 cw3; do
    p=$X/exocc_${s}_${t}_scan_${t}
    [ $s = starthere ] && p=$X/exocc_starthere_scan_${t}
    r=$(for j in 1 2; do grep -h "Final R-work" ${p}_j${j}_001.log 2>/dev/null | tail -1 | awk '{print $NF}'; done)
    line+=$(echo $r | awk '{if(NF==2)printf " %-22s",sprintf("%.4f (%s %s)",($1+$2)/2,$1,$2); else printf " %-22s","-"}')
  done
  echo "$line"
done
