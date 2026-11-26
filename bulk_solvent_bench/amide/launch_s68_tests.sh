#!/bin/bash
# Charge tests at the 0.68 M crystallisation salt (2026-11-25): build each system (build_s68.csh), print its His labels
# and charges (his_check.py), then run it exactly as the restraint-class reference arm ff19_opc_rc_dw_w (b5sd1,
# code_391d233, min_align_weight=0, uniform B 5, 200 itr).  Default arms cEE EE ED PD; DD and MD are built by
# setup_s68_builds.sh but not run (the refined H omit maps put His54 at HIP, not HID or a HIP/HID mixture).
# Run on octamus1:  ssh octamus1 'bash /scratch/jamesh/projects/amber/amide/launch_s68_tests.sh [cEE EE ...]'
S=/scratch/jamesh/projects/amber; H=$S/1aho_hurry_test; A=$S/amide; code=code_391d233
arms="${*:-cEE EE ED PD}"
for v in $arms; do
  b=build_s68_$v; name=ff19_opc_s68_$v
  if pgrep -f "run_local.csh $name " >/dev/null || pgrep -f "build_s68.csh $H/$b" >/dev/null; then
    echo "$name already running, skipping"; continue
  fi
  ( setsid nohup bash -c "
      if [ ! -e $H/$b/Prod.rst7 ]; then /bin/tcsh -f $A/build_s68.csh $H/$b; fi
      [ -e $H/$b/Prod.rst7 ] || { echo '$name: build failed'; exit 9; }
      /scratch/jamesh/programs/amber22/miniconda/bin/python3 $A/his_check.py $H/$b/xtal.prmtop
      cd $H; BUILD=$b /bin/tcsh -f $S/run_local.csh $name 200 b5sd1 $code min_align_weight=0
    " </dev/null >> $S/$name.log 2>&1 & )
  echo "launched $name (build $b): $(tr '\n' ' ' < $H/$b/protonation.txt)"
  sleep 5
done
