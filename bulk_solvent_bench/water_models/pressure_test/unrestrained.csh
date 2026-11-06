#! /bin/tcsh -f
# usage: unrestrained.csh   (on octamus1)
# the release-and-read protocol: every restraint off, same cell, 1 ns, read the mean reported pressure through
# ntp=4 (production), ntp=1 and ntp=2 (Berendsen, taup 1e11), from w1_100ns prod_50 and prod_100, plus the Monte
# Carlo barostat at 1 bar from prod_50 (prod_100's is release/mc_u_1) -> unrestrained/<seg>_<variant>.out
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/unrestrained
if( ! -d $D ) mkdir $D
cd $D
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
foreach i ( 50 100 )
  foreach v ( u4 u1 u2 mc )
    if( $v == mc && $i == 100 ) continue
    set ntp = "ntb=2,ntp=4,"
    if( $v == u1 ) set ntp = "ntb=2,ntp=1,barostat=1,taup=99999999999,"
    if( $v == u2 ) set ntp = "ntb=2,ntp=2,barostat=1,taup=99999999999,"
    if( $v == mc ) set ntp = "ntb=2,ntp=1,barostat=2,mcbarint=100,pres0=1.0,"
    awk -v ntp=$ntp '/ ntr=/{sub(/ntr=1/,"ntr=0")} / nstlim=/{print "  nstlim=500000,";next}\
      / ntpr=/{print "  ntpr=500,";next} / ntwx=/{print "  ntwx=0,";next} / ntb=/{print "  " ntp;next} {print}' \
      $S/prod.in >! prod_${i}_$v.in
    srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP -O -i prod_${i}_$v.in -o prod_${i}_$v.out -p $S/xtal.prmtop \
      -c $S/prod_${i}.rst7 -ref $S/ref.crd -r prod_${i}_$v.rst7 -inf prod_${i}_$v.mdinfo >& prod_${i}_$v.log &
  end
end
wait
echo "unrestrained done"
