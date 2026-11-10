#!/bin/bash
# usage: baseline.sh   (on octamus1)  the comparisons so far, again with the 7KR0_SETTINGS.md recipe (run_ref.sh,
# 10 cycles, 3 seeds): PRISM's S0112 with its waters, anisotropic (ens) and ANISOU-stripped (ens_iso), against flat,
# the 2 ns and the 20 ns MD maps; and the 2 ns variant maps (restraint weight 0.3, 0.1, 0.03 per sigma; altloc merge
# 0.4 A) on ens.  In test_7kr0/rewet/base; results in base/summary.txt
S=/home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts
cd /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/rewet/base || exit 9
for s in 1 2 3; do
  for m in ens ens_iso; do
    for sol in flat md2 md20; do
      map=flat; [ $sol != flat ] && map=$sol.mtz
      bash $S/run_ref.sh $m.pdb $map ${m}_${sol}_s$s $s > ${m}_${sol}_s$s.out 2>&1 &
    done
  done
  for v in wps0.3 wps0.1 wps0.03 merge0.4; do
    bash $S/run_ref.sh ens.pdb v_$v.mtz ens_${v}_s$s $s > ens_${v}_s$s.out 2>&1 &
  done
done
wait
for t in ens_flat ens_md2 ens_md20 ens_iso_flat ens_iso_md2 ens_iso_md20 ens_wps0.3 ens_wps0.1 ens_wps0.03 ens_merge0.4; do
  printf "%-14s" $t
  for s in 1 2 3; do awk '/Final R-work/{gsub(",", ""); r=$NF; w=$4} END{printf("  %s/%s", w==""?"-":w, r==""?"-":r)}' ${t}_s$s.log; done
  echo
done > summary.txt
cat summary.txt
