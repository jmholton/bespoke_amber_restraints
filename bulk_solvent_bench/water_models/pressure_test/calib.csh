#! /bin/tcsh -f
# usage: calib.csh   (on octamus1)
# 1) which pressure path reads right: a pure OPC box (~2000 waters, production settings: 287 K, cut 9, Langevin),
#    equilibrated by the Monte Carlo barostat at 1 bar (no virial involved), then at its last volume: 0.5 ns MC
#    (-> <V>, kappa), ntp=4 (as production), ntp=1 Berendsen taup 1e11, ntp=2 taup 1e11 (as the 2019 runs)
#    -> calib/water_*.out
# 2) the same paths on frozen coordinates: w1_100ns prod_50, prod_100 with dt = 1e-5 ps, 20 steps, every step
#    printed, restrained and not -> calib/frozen_<seg>_<r|u><ntp>.out
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/calib
if( ! -d $D ) mkdir $D
cd $D
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
cat << EOF >! tleap.in
source leaprc.water.opc
x = loadpdb $S/opc.pdb
solvatebox x OPCBOX 20.0
saveamberparm x water.prmtop water.rst7
quit
EOF
tleap -f tleap.in >&! tleap.log
# mdin pieces in production's settings
set common = "ntf=2, ntc=2, temp0=287, cut=9, ntt=3, gamma_ln=1.0, ntwx=0, ig=-1,"
cat << EOF >! min.in
min
 &cntrl
  imin=1, maxcyc=1000, ntb=1, cut=9, ntpr=200,
 /
EOF
cat << EOF >! heat.in
heat, Monte Carlo barostat 1 bar
 &cntrl
  imin=0, ntx=1, irest=0, tempi=287, nstlim=150000, dt=0.002, $common
  ntpr=1000, ntwr=150000, ntb=2, ntp=1, barostat=2, mcbarint=100, pres0=1.0,
 /
EOF
foreach v ( mc p4 p1 p2 )
  set ntp = "ntb=2, ntp=1, barostat=2, mcbarint=100, pres0=1.0,"
  if( $v == p4 ) set ntp = "ntb=2, ntp=4,"
  if( $v == p1 ) set ntp = "ntb=2, ntp=1, barostat=1, taup=99999999999,"
  if( $v == p2 ) set ntp = "ntb=2, ntp=2, barostat=1, taup=99999999999,"
  cat << EOF >! water_$v.in
water $v
 &cntrl
  imin=0, ntx=5, irest=1, nstlim=250000, dt=0.002, $common
  ntpr=500, ntwr=250000, $ntp
 /
 &ewald
  netfrc=0,
 /
EOF
end
set gpu = "srun --partition=gpu --gres=gpu:1"
( $gpu pmemd.cuda_SPFP -O -i min.in -o min.out -p water.prmtop -c water.rst7 -r min.rst7 -inf min.mdinfo && \
  $gpu pmemd.cuda_SPFP -O -i heat.in -o heat.out -p water.prmtop -c min.rst7 -r heat.rst7 -inf heat.mdinfo ) >& water_eq.log
foreach v ( mc p4 p1 p2 )
  $gpu pmemd.cuda_SPFP -O -i water_$v.in -o water_$v.out -p water.prmtop -c heat.rst7 -r water_$v.rst7 \
    -inf water_$v.mdinfo >& water_$v.log &
end
# frozen coordinates
echo "#! /bin/tcsh -f" >! inner.csh
foreach i ( 50 100 )
  foreach v ( r4 u4 r1 u1 r2 u2 )
    set ntr = 0
    if( $v =~ r* ) set ntr = 1
    set ntp = "ntb=2,ntp=4,"
    if( $v =~ *1 ) set ntp = "ntb=2,ntp=1,barostat=1,taup=99999999999,"
    if( $v =~ *2 ) set ntp = "ntb=2,ntp=2,barostat=1,taup=99999999999,"
    awk -v ntr=$ntr -v ntp=$ntp '/ ntr=/{sub(/ntr=[01]/,"ntr=" ntr)} / nstlim=/{print "  nstlim=20,";next}\
      / dt=/{print "  dt=0.00001,";next} / ntpr=/{print "  ntpr=1,";next} / ntwx=/{print "  ntwx=0,";next}\
      / ntwr=/{print "  ntwr=20,";next} / ig=/{print "  ig=71277,";next} / ntb=/{print "  " ntp;next} {print}' \
      $S/prod.in >! frozen_${i}_$v.in
    echo "pmemd.cuda_SPFP -O -i frozen_${i}_$v.in -o frozen_${i}_$v.out -p $S/xtal.prmtop -c $S/prod_${i}.rst7 -ref $S/ref.crd -r frozen_${i}_$v.rst7 -inf frozen_${i}_$v.mdinfo" >> inner.csh
  end
end
$gpu tcsh -f inner.csh >& inner.log
echo "frozen done"
wait
echo "calib done"
