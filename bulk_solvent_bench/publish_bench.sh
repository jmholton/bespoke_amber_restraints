#!/bin/bash
# usage: bash publish_bench.sh   (in this bench, /home2/jamesh/projects/amber/bulk_solvent)
# Copy the bench's notes, own scripts, configs and small result tables into the bespoke_amber_restraints clone as
# bulk_solvent_bench/, with the bench's layout kept so the paths in MD_BULK_SOLVENT_STRATEGY.md resolve.  It only
# copies (any old bulk_solvent_bench/ is replaced); commit and push in the clone afterwards.  Only files tracked
# here are copied, so commit here first.
# Left out: the outputs of the run directories (any directory holding a current_restraints.txt, with its
# subdirectories), copies of the repo's own scripts (same file name as a file of the repo), PRISM's files
# (test_7kr0/prism), the per-map frame lists (water_models/hurry/hl_*/), generated Amber inputs over 10 kB, and
# symlinks (the compute_settings.sourceme links point at the site file outside the bench; the repo has the template).
B=/home2/jamesh/projects/amber/bulk_solvent
R=$B/bespoke_amber_restraints
cd $B || exit 9
tmp=$(mktemp -d)
git -C $R ls-files | grep -v '^bulk_solvent_bench/' | xargs -n1 basename | sort -u > $tmp/tool
git ls-files | grep '/current_restraints.txt$' | sed 's#/current_restraints.txt$##' | sort -u > $tmp/run
git ls-files | awk -v rd=$tmp/run -v tn=$tmp/tool '
  BEGIN { while ((getline l < rd) > 0) run[l] = 1; while ((getline l < tn) > 0) tool[l] = 1 }
  { n = split($0, p, "/"); d = ""
    for (i = 1; i < n; ++i) { d = (i == 1 ? p[1] : d "/" p[i]); if (d in run) next }
    if (p[n] in tool || $0 ~ /^test_7kr0\/prism\// || $0 ~ /^water_models\/hurry\/hl_[^\/]*\//) next
    print }' |
while read -r f; do
  [ -L "$f" ] && continue
  case "$f" in *.in) [ $(stat -c %s "$f") -gt 10000 ] && continue ;; esac
  echo "$f"
done > $tmp/list
rm -rf $R/bulk_solvent_bench
mkdir -p $R/bulk_solvent_bench
rsync -a --files-from=$tmp/list ./ $R/bulk_solvent_bench/
echo "$(wc -l < $tmp/list) files, $(du -sh $R/bulk_solvent_bench | cut -f1) -> $R/bulk_solvent_bench"
rm -rf $tmp
