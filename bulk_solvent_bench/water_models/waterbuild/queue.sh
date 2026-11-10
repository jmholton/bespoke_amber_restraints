#!/bin/bash
# usage: cd /scratch/jamesh/projects/amber/waterbuild; setsid nohup bash queue.sh > queue.log 2>&1 < /dev/null &
# Runs the batches of maps.txt in order through hurry_chain.csh (all tags of a batch side by side), from /scratch
# only.  A tag counts as done when both starts finished round 3 (its r3 column in results.txt; done/<tag>);
# re-running skips done tags.  touch stop
# to quit before the next batch.  After each batch: results.txt (rfree_table.sh over every tag), mirrored to
# /petabyte/jamesh/waterbuild/ in the background (a hung mount cannot stall the queue).
W=/scratch/jamesh/projects/amber/waterbuild
source $W/env.sh
cd $W || exit 9
for b in $(awk '!/^#/{print $1}' maps.txt | uniq); do
  if [ -e $W/stop ]; then echo "stop file: quitting before batch $b $(date +%H:%M)"; break; fi
  todo=""
  for t in $(awk -v b=$b '!/^#/ && $1==b {print $2}' maps.txt); do [ -e done/$t ] || todo="$todo $t"; done
  [ -z "$todo" ] && continue
  echo "batch $b start $(date +%H:%M):$todo"
  ./hurry_chain.csh $todo > chain_batch$b.log 2>&1
  bash rfree_table.sh $(awk '!/^#/{print $2}' maps.txt) > results.txt
  for t in $todo; do awk -v t=$t '$1==t && $NF ~ /\)$/ {f=1} END{exit !f}' results.txt && touch done/$t; done
  echo "batch $b done $(date +%H:%M): $(for t in $todo; do [ -e done/$t ] && printf "%s " $t; done)"
  ( mkdir -p /petabyte/jamesh/waterbuild && cp results.txt queue.log /petabyte/jamesh/waterbuild/ ) > /dev/null 2>&1 &
done
echo "queue done $(date +%H:%M)"
