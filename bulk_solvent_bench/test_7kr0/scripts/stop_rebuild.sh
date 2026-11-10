#!/bin/bash
# usage: stop_rebuild.sh   (on octamus1)  stop the superseded 3-sigma rebuild tracks and every process running in
# test_7kr0/rebuild or test_7kr0/rebuild_iso (the refine/ jobs are left alone)
B=/home2/jamesh/projects/amber/bulk_solvent/test_7kr0
pkill -u "$USER" -f "scripts/track.csh" && echo "stopped track.csh"
for p in $(pgrep -u "$USER"); do
  d=$(readlink /proc/$p/cwd 2>/dev/null) || continue
  case "$d" in
    $B/rebuild|$B/rebuild_iso) echo "kill $p ($(basename "$d")): $(ps -o args= -p $p | cut -c1-90)"; kill $p ;;
  esac
done
sleep 3
left=0
for p in $(pgrep -u "$USER"); do
  d=$(readlink /proc/$p/cwd 2>/dev/null) || continue
  case "$d" in $B/rebuild|$B/rebuild_iso) left=$((left+1)) ;; esac
done
echo "$left processes left in rebuild dirs"
