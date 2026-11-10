#!/bin/bash
# usage: time_reduce.sh model.pdb reduce-args...   time one Reduce call through scripts/reduce_wrap, as phenix makes it
source /programs/phenix-2.2-6143/phenix_env.sh
export PATH=/home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts/reduce_wrap:$PATH
pdb=$1; shift
s=$(date +%s)
n=$(molprobity.reduce "$@" < "$pdb" 2>/dev/null | grep -c '^ATOM')
echo "$(basename "$pdb") [$*]: $n ATOM lines out, $(( $(date +%s) - s )) s"
