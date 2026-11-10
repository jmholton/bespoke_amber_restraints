#!/bin/bash
# usage: rewet_track.sh dir dry.pdb flat|map.mtz [max_rounds]   (on octamus1; dir holds 7kr0.mtz, the map and dry.pdb)
# Dry-and-re-wet of the PRISM 7KR0 ensemble against one bulk solvent, with the 7KR0_SETTINGS.md recipe (run_ref.sh):
# round 0 refines the water-free model; each round after that adds waters with ens_waters.py (from 6 sigma down to
# 2.5 only as the higher levels run out; per-conformer range-of-motion 2.6 A; held occupancies) to the median-seed
# model of the round before and refines it.  Three seeds per round (1, 2, 3), the model carried on is the median
# R-free seed's.  Stops when nothing more can be placed, or after max_rounds (10), then a 20-cycle convergence step.
# Resumes: finished refinements (<prefix>_001.pdb) and builds (w<n>.pdb) are kept.  Summary lines in track.log.
S=/home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts
dir=$1; dry=$2; sol=$3; maxr=${4:-10}
cd "$dir" || exit 9
log() { echo "$(date +%H:%M) $*" >> track.log; }

refine_round() {   # tag model cycles: three seeds side by side, then one retry for any seed that left no model
  local tag=$1 model=$2 cyc=$3 s try
  for try in 1 2; do
    for s in 1 2 3; do
      [ -e r${tag}_s${s}_001.pdb ] || bash $S/run_ref.sh $model $sol r${tag}_s${s} $s $cyc > r${tag}_s${s}.out 2>&1 &
    done
    wait
  done
}
median_seed() {    # tag: the seed with the median R-free among those that finished
  local s r
  for s in 1 2 3; do
    r=$(awk '/Final R-work/{gsub(",", ""); r=$NF} END{print r}' r$1_s${s}.log 2>/dev/null)
    [ -n "$r" ] && [ -e r$1_s${s}_001.pdb ] && echo "$r $s"
  done | sort -g | awk '{s[NR]=$2} END{if(NR) print s[int((NR+1)/2)]}'
}
rfrees() { local s; for s in 1 2 3; do awk '/Final R-work/{gsub(",", ""); r=$NF} END{printf("%s ", r==""?"-":r)}' r$1_s${s}.log 2>/dev/null; done; }

log "start: $dry against $sol"
refine_round 0 $dry 10
m=$(median_seed 0)
[ -z "$m" ] && { log "round 0: no model from any seed"; exit 9; }
log "round 0 (no waters): R-free $(rfrees 0)median seed $m"
n=1
while [ $n -le $maxr ]; do
  p=$((n - 1)); mp=$(median_seed $p)
  if [ ! -e w${n}.pdb ]; then
    /programs/phenix-2.2-6143/bin/phenix.python $S/ens_waters.py r${p}_s${mp}_001.pdb r${p}_s${mp}_001.mtz \
      w${n}.pdb > build_w${n}.log 2>&1
    rc=$?
    if [ $rc -eq 3 ]; then log "round ${n}: nothing left to place at >= 2.5 sigma"; break; fi
    [ $rc -ne 0 ] && { log "round ${n}: ens_waters.py failed (build_w${n}.log)"; exit 9; }
  fi
  refine_round $n w${n}.pdb 10
  m=$(median_seed $n)
  [ -z "$m" ] && { log "round ${n}: no model from any seed"; exit 9; }
  log "round ${n}: $(grep '^wrote' build_w${n}.log | sed 's/^wrote w[0-9]*.pdb: //'); R-free $(rfrees $n)median seed $m"
  n=$((n + 1))
done
last=$((n - 1)); ml=$(median_seed $last)
refine_round c r${last}_s${ml}_001.pdb 20
log "convergence (20 cycles from round $last): R-free $(rfrees c)median seed $(median_seed c)"
