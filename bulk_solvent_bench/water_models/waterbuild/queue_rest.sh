#!/bin/bash
# usage: cd /scratch/jamesh/projects/amber/waterbuild; touch stop; setsid nohup bash queue_rest.sh > queue_rest.log 2>&1 < /dev/null &
# 2026-11-05 21:00: the power goes off tonight, before queue.sh (one batch an hour) would get past batch 3.  This
# runs every batch after 1 at once instead, one hurry_chain.csh per batch side by side (78 maps = 156 refinements of
# ~2.2 GB; octamus1 has 224 cores / 448 threads), while queue.sh finishes batch 1 and quits (the stop file).
# Every 10 min results_live.txt (rfree_table.sh: the rounds finished so far) and this log go to
# /petabyte/jamesh/waterbuild/.  At the end results.txt and done/<tag> as in queue.sh.  To resume after the
# outage: rm stop; bash queue.sh (it redoes from the start every tag without both round-3 refinements).
W=/scratch/jamesh/projects/amber/waterbuild
source $W/env.sh
cd $W || exit 9
all=$(awk '!/^#/{print $2}' maps.txt)
mirror() {
  bash rfree_table.sh $all > results_live.txt
  ( mkdir -p /petabyte/jamesh/waterbuild && cp results_live.txt queue_rest.log /petabyte/jamesh/waterbuild/ ) > /dev/null 2>&1 &
}
pids=""
for b in $(awk '!/^#/ && $1 > 1 {print $1}' maps.txt | uniq); do
  todo=""
  for t in $(awk -v b=$b '!/^#/ && $1==b {print $2}' maps.txt); do [ -e done/$t ] || todo="$todo $t"; done
  [ -z "$todo" ] && continue
  echo "batch $b start $(date +%H:%M):$todo"
  ./hurry_chain.csh $todo > chain_batch$b.log 2>&1 &
  pids="$pids $!"
done
# wait only for the chains: a background cp stuck on a dead /petabyte mount must not hold up the results
n=0
while kill -0 $pids 2> /dev/null; do
  sleep 60
  n=$((n + 1))
  [ $((n % 10)) -eq 0 ] && mirror
done
wait $pids
bash rfree_table.sh $all > results.txt
for t in $all; do awk -v t=$t '$1==t && $NF ~ /\)$/ {f=1} END{exit !f}' results.txt && touch done/$t; done
echo "batches 2-$b done $(date +%H:%M): $(ls done | wc -l) of $(echo $all | wc -w) maps through round 3"
mirror
( cp results.txt /petabyte/jamesh/waterbuild/ ) > /dev/null 2>&1 &
