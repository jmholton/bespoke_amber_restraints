#! /bin/tcsh -f
# usage: release.csh   (on octamus1)
# 1) what happens to the reported pressure in the first 20 ps after the restraints come off: w1_100ns prod_10..100,
#    20 ps each as r4 (restrained, production ntp=4), u4 (unrestrained, ntp=4) and r1 (restrained, ntp=1 Berendsen
#    taup 1e11), all with the same Langevin seed, printed every 0.1 ps -> release/<seg>_<variant>.out
# 2) virial-free pressure: Monte Carlo barostat (barostat=2) from prod_100, 1 ns at 1 and 1000 bar, restraints on
#    and off -> release/mc_<r|u>_<pres0>.out; the pressure that keeps the crystal's volume is where <V> = 56230.7
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/release
if( ! -d $D ) mkdir $D
cd $D
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
foreach p ( 1 1000 )
  foreach v ( r u )
    set ntr = 0
    if( $v == r ) set ntr = 1
    awk -v ntr=$ntr -v p=$p '/ ntr=/{sub(/ntr=[01]/,"ntr=" ntr)} / nstlim=/{print "  nstlim=500000,";next}\
      / ntpr=/{print "  ntpr=500,";next} / ntwx=/{print "  ntwx=0,";next}\
      / ntb=/{print "  ntb=2,ntp=1,barostat=2,mcbarint=100,pres0=" p ".0,";next} {print}' $S/prod.in >! mc_${v}_$p.in
    srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP -O -i mc_${v}_$p.in -o mc_${v}_$p.out -p $S/xtal.prmtop \
      -c $S/prod_100.rst7 -ref $S/ref.crd -r mc_${v}_$p.rst7 -inf mc_${v}_$p.mdinfo >& mc_${v}_$p.log &
  end
end
echo "#! /bin/tcsh -f" >! inner.csh
foreach i ( 10 20 30 40 50 60 70 80 90 100 )
  foreach v ( r4 u4 r1 )
    set ntr = 0
    if( $v =~ r* ) set ntr = 1
    set ntp = "ntb=2,ntp=4,"
    if( $v == r1 ) set ntp = "ntb=2,ntp=1,barostat=1,taup=99999999999,"
    awk -v ntr=$ntr -v ntp=$ntp '/ ntr=/{sub(/ntr=[01]/,"ntr=" ntr)} / nstlim=/{print "  nstlim=10000,";next}\
      / ntpr=/{print "  ntpr=50,";next} / ntwx=/{print "  ntwx=0,";next} / ntwr=/{print "  ntwr=10000,";next}\
      / ig=/{print "  ig=71277,";next} / ntb=/{print "  " ntp;next} {print}' $S/prod.in >! prod_${i}_$v.in
    echo "pmemd.cuda_SPFP -O -i prod_${i}_$v.in -o prod_${i}_$v.out -p $S/xtal.prmtop -c $S/prod_${i}.rst7 -ref $S/ref.crd -r prod_${i}_$v.rst7 -inf prod_${i}_$v.mdinfo" >> inner.csh
  end
end
srun --partition=gpu --gres=gpu:1 tcsh -f inner.csh >& inner.log
echo "release runs done"
wait
echo "mc runs done"
