#!/bin/bash
# usage: count_in.sh dir   (on octamus1)  how many of my processes run under dir, and which drivers are alive
dir=$1; n=0
for p in $(pgrep -u "$USER"); do
  d=$(readlink /proc/$p/cwd 2>/dev/null) || continue
  case "$d" in "$dir"|"$dir"/*) n=$((n+1)) ;; esac
done
echo "$n processes under $dir"
ps -u "$USER" -o pid=,args= | grep -E "[r]ewet_track.sh|[b]aseline.sh" | cut -c1-110
