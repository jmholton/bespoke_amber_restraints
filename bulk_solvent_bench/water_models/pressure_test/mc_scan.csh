#! /bin/tcsh -f
# usage: mc_scan.csh   (on octamus1)  restrained Monte Carlo barostat at 1 bar, 1 ns from prod_2 of the scan runs
# (pmemd scales the restraint references with the box) -> mc_scan/<run>.out; the true pressure of each run's
# restrained cell follows from <V> and the w1 restrained compressibility (release/mc_r_1, mc_r_1000)
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test/mc_scan
if( ! -d $D ) mkdir $D
cd $D
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
foreach t ( pm500 w1 pp500 pp1500 pp3000 )
  awk '/ nstlim=/{print "  nstlim=500000,";next} / ntpr=/{print "  ntpr=500,";next} / ntwx=/{print "  ntwx=0,";next}\
    / ntb=/{print "  ntb=2,ntp=1,barostat=2,mcbarint=100,pres0=1.0,";next} {print}' $S/$t/prod.in >! $t.in
  srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP -O -i $t.in -o $t.out -p $S/$t/xtal.prmtop -c $S/$t/prod_2.rst7 \
    -ref $S/$t/ref.crd -r $t.rst7 -inf $t.mdinfo >& $t.log &
end
wait
echo "mc_scan done"
