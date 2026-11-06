#! /bin/tcsh -f
# usage: windows_driver.csh   (on bl831)
# follows w1_100ns production: as each window's last segment finishes, its solvent map (window_map.csh on
# octamus1), then for the "r" windows the example refinement from 2 jiggled starts (as eval.csh) with tag
# scan_win<first>_<last>; "c" windows are only for map correlations (window_cc.py).  1_2 is the w1 map itself.
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan
set D = $S/w1_100ns
set X = $S/../example_check
set windows = ( 1_1_r 1_2_c 1_4_r 1_8_r 1_10_c 11_12_r 1_16_r 11_20_c 21_30_c 31_32_r 1_32_r 31_40_c 41_50_c \
  51_52_r 51_60_c 61_70_c 1_64_r 71_72_r 71_80_c 81_90_c 91_92_r 91_100_c 1_100_r 51_100_r )
foreach w ( $windows )
  set first = `echo $w | awk -F_ '{print $1}'`
  set last = `echo $w | awk -F_ '{print $2}'`
  set what = `echo $w | awk -F_ '{print $3}'`
  set tag = win${first}_${last}
  while ( 1 )
    if( -e $D/prod_${last}.out ) then
      grep -q "Total wall time" $D/prod_${last}.out
      if( ! $status ) break
    endif
    grep -q "^ERROR" $D/run.log
    if( ! $status ) then
      echo "MD failed: `grep ^ERROR $D/run.log`"
      exit 9
    endif
    sleep 120
  end
  if( ! -e $D/windows/map_${first}_${last}.mtz ) ssh octamus1 $S/window_map.csh $first $last
  # NFS can be slow to show the new file here
  set n = 0
  while ( ! -e $D/windows/map_${first}_${last}.mtz && $n < 30 )
    sleep 10
    @ n++
  end
  if( "$what" != "r" ) continue
  if( -e $X/exocc_starthere_scan_${tag}_j1.log ) continue
  while ( `ps -u $USER -o args | grep -c "^/home/programs/phenix-2.2-6143/bin/python.*command_line"` > 8 )
    sleep 60
  end
  cp $D/windows/map_${first}_${last}.mtz $X/scan_${tag}_map.mtz
  cd $X
  foreach j ( 1 2 )
    ./run_exocc.csh starthere scan_${tag}_map.mtz scan_$tag $j >&! exocc_starthere_scan_${tag}_j${j}.out &
  end
  echo "${tag}: refining"
end
wait
source /programs/phenix-2.2-6143/phenix_env.csh
phenix.python $S/window_cc.py >&! $D/windows/cc.txt
cd $X
foreach w ( $windows )
  set tag = `echo $w | awk -F_ '{print "win" $1 "_" $2}'`
  if( "$w" =~ *_c ) continue
  cat exocc_starthere_scan_${tag}_j?.log | awk -v t=$tag '/Final R-work/{w+=$4;f+=$7;++n} END{printf("%s Rwork %.4f Rfree %.4f\n",t,w/n,f/n)}'
end
echo "windows done"
