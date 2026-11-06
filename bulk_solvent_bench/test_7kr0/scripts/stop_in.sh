#!/bin/bash
# usage: stop_in.sh dir   (on octamus1)  stop every process of mine whose working directory is dir or below it (the
# drivers run there too, so two passes catch anything they start in between).  No pkill -f: a pattern given on the
# ssh command line matches that command line itself.
dir=$1
for try in 1 2; do
  for p in $(pgrep -u "$USER"); do
    d=$(readlink /proc/$p/cwd 2>/dev/null) || continue
    case "$d" in "$dir"|"$dir"/*) kill $p 2>/dev/null ;; esac
  done
  sleep 3
done
left=0
for p in $(pgrep -u "$USER"); do
  d=$(readlink /proc/$p/cwd 2>/dev/null) || continue
  case "$d" in "$dir"|"$dir"/*) left=$((left+1)) ;; esac
done
echo "$left processes left under $dir"
