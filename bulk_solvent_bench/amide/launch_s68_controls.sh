#!/bin/bash
# Controls for the 0.68 M charge tests (2026-11-25; user: "voltron fairly idle. do more controls"), each run exactly as
# the s68 arms (b5sd1, code_391d233, min_align_weight=0, uniform B 5, 200 itr):
#   PDp    build_s68_PDp  His54 ring relabelled to the P orientation (setup_s68_flip.sh): known-wrong ring labels,
#                         a positive control for whether the MD measures see the ring N/C label
#   DD     build_s68_DD   HID54 (built by setup_s68_builds.sh, not run before): completes HID / HIE / HIP at His54
#   PD_r2  build_s68_PD   replicate of PD (new random seeds, ig=-1, same start): run-to-run scatter
#   ED_r2  build_s68_ED   replicate of ED
# Run on octamus1:  ssh octamus1 'bash /scratch/jamesh/projects/amber/amide/launch_s68_controls.sh [PDp:build_s68_PDp ...]'
S=/scratch/jamesh/projects/amber; H=$S/1aho_hurry_test; A=$S/amide; code=code_391d233
arms="${*:-PDp:build_s68_PDp DD:build_s68_DD PD_r2:build_s68_PD ED_r2:build_s68_ED}"
for spec in $arms; do
  v=${spec%%:*}; b=${spec#*:}; name=ff19_opc_s68_$v
  if pgrep -f "run_local.csh $name " >/dev/null || pgrep -f "build_s68.csh $H/$b\$" >/dev/null; then
    echo "$name already running, skipping"; continue
  fi
  [ -e $H/$name ] && { echo "$H/$name exists, skipping"; continue; }
  ( setsid nohup bash -c "
      if [ ! -e $H/$b/Prod.rst7 ]; then /bin/tcsh -f $A/build_s68.csh $H/$b; fi
      [ -e $H/$b/Prod.rst7 ] || { echo '$name: build failed'; exit 9; }
      /scratch/jamesh/programs/amber22/miniconda/bin/python3 $A/his_check.py $H/$b/xtal.prmtop
      cd $H; BUILD=$b /bin/tcsh -f $S/run_local.csh $name 200 b5sd1 $code min_align_weight=0
    " </dev/null >> $S/$name.log 2>&1 & )
  echo "launched $name (build $b): $(tr '\n' ' ' < $H/$b/protonation.txt)"
  sleep 5
done
