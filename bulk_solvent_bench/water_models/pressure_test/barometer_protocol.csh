#! /bin/tcsh -f
# usage: barometer_protocol.csh   (on octamus1)
# optimize_weights_runme.com's barometer: every restraint off, barometer_cycles=10000 steps at dt=0.0005 (5 ps),
# ntpr=1000, ntb=2,ntp=4, mean of the printed pressures.  Here also through ntp=1 (Berendsen, taup 1e11).
# 1) 1aho cell, w1_100ns prod_10..prod_100 -> barometer_protocol/xtal_<seg>_<4|1>.out
# 2) the calib/ OPC box (true pressure +103 bar at its volume) at dt 0.5, 1, 2 fs, 50 ps each
#    -> barometer_protocol/water_<dt fs>_<4|1>.out
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns
set W = /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/calib
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/barometer_protocol
if( ! -d $D ) mkdir $D
cd $D
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
echo "#! /bin/tcsh -f" >! inner.csh
foreach i ( 10 20 30 40 50 60 70 80 90 100 )
  foreach p ( 4 1 )
    set ntp = "ntb=2,ntp=4,"
    if( $p == 1 ) set ntp = "ntb=2,ntp=1,barostat=1,taup=99999999999,"
    awk -v ntp=$ntp '/ ntr=/{sub(/ntr=1/,"ntr=0")} / nstlim=/{print "  nstlim=10000,";next}\
      / dt=/{print "  dt=0.0005,";next} / ntpr=/{print "  ntpr=1000,";next} / ntwx=/{print "  ntwx=0,";next}\
      / ntwr=/{print "  ntwr=10000,";next} / ntb=/{print "  " ntp;next} {print}' $S/prod.in >! xtal_${i}_$p.in
    echo "pmemd.cuda_SPFP -O -i xtal_${i}_$p.in -o xtal_${i}_$p.out -p $S/xtal.prmtop -c $S/prod_${i}.rst7 -ref $S/ref.crd -r xtal_${i}_$p.rst7 -inf xtal_${i}_$p.mdinfo" >> inner.csh
  end
end
foreach dt ( 0.5 1 2 )
  set n = `echo $dt | awk '{print int(50000/$1+0.5)}'`
  set pr = `echo $dt | awk '{print int(100/$1+0.5)}'`
  foreach p ( 4 1 )
    set ntp = "ntb=2, ntp=4,"
    if( $p == 1 ) set ntp = "ntb=2, ntp=1, barostat=1, taup=99999999999,"
    awk -v ntp="$ntp" -v n=$n -v pr=$pr -v dt=$dt '/nstlim=/{sub(/nstlim=[0-9]+/,"nstlim=" n); sub(/dt=[0-9.]+/,"dt=" dt/1000)}\
      /ntb=/{$0="  ntpr=" pr ", ntwr=" n ", " ntp} {print}' \
      $W/water_p4.in >! water_${dt}_$p.in
    echo "pmemd.cuda_SPFP -O -i water_${dt}_$p.in -o water_${dt}_$p.out -p $W/water.prmtop -c $W/heat.rst7 -r water_${dt}_$p.rst7 -inf water_${dt}_$p.mdinfo" >> inner.csh
  end
end
srun --partition=gpu --gres=gpu:1 tcsh -f inner.csh >& inner.log
echo "barometer_protocol done"
