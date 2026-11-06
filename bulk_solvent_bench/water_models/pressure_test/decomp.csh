#! /bin/tcsh -f
# usage: decomp.csh   (on octamus1)
# the pressure pmemd reports for identical coordinates and velocities, restraints on vs off, and through the
# production code path (ntp=4) vs isotropic Berendsen with taup 1e11 (ntp=1): 2 steps from each restart, step 1
# read.  Frames: w1_100ns prod_10..prod_100 and prod_1, prod_2 of six scan runs -> decomp/<run>_<seg>_<variant>.out
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/decomp
mkdir -p $D
cd $D
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
set frames = ()
foreach i ( 10 20 30 40 50 60 70 80 90 100 )
  set frames = ( $frames w1_100ns+prod_$i )
end
foreach t ( w01 w1 w10 pm500 pp1500 pp3000 )
  set frames = ( $frames ${t}+prod_1 ${t}+prod_2 )
end
echo "#! /bin/tcsh -f" >! inner.csh
foreach f ( $frames )
  set run = `echo $f | awk -F+ '{print $1}'`
  set seg = `echo $f | awk -F+ '{print $2}'`
  foreach v ( r4 u4 r1 u1 )
    set ntr = `echo $v | awk '{print (substr($1,1,1)=="r")?1:0}'`
    set ntp = "ntb=2,ntp=4,"
    if( $v =~ *1 ) set ntp = "ntb=2,ntp=1,barostat=1,taup=99999999999,"
    awk -v ntr=$ntr -v ntp=$ntp '/ ntr=/{sub(/ntr=[01]/,"ntr=" ntr)} / nstlim=/{print "  nstlim=2,";next}\
      / ntpr=/{print "  ntpr=1,";next} / ntwx=/{print "  ntwx=0,";next} / ntwr=/{print "  ntwr=2,";next}\
      / ig=/{print "  ig=71277,";next} / ntb=/{print "  " ntp;next} {print}' $S/$run/prod.in >! ${run}_${seg}_$v.in
    echo "pmemd.cuda_SPFP -O -i ${run}_${seg}_$v.in -o ${run}_${seg}_$v.out -p $S/$run/xtal.prmtop -c $S/$run/${seg}.rst7 -ref $S/$run/ref.crd -r /dev/null -inf /dev/null" >> inner.csh
  end
end
srun --partition=gpu --gres=gpu:1 tcsh -f inner.csh
echo "decomp done: `ls *.out | wc -l` outputs"
