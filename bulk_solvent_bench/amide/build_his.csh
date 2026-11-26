#! /bin/tcsh -f
# leap2amber build of one His-protonation test system (setup_his_builds.sh made the inputs).  The same leap2amber call
# and code (code_noh) as build_one.csh, which built build_2s2r_amide; software from programs.sourceme and MD through
# the SITE compute settings (SLURM), as the run wrappers now do.   usage: build_his.csh build_dir
set d = $1
source /home2/jamesh/projects/amber/bulk_solvent/programs.sourceme
setenv OMP_NUM_THREADS 4
set path = ( /scratch/jamesh/projects/amber/1aho_hurry_test/code_noh $path )
cd $d
echo "build $d on `hostname -s`: `tr '\n' ';' < protonation.txt` `date`"
leap2amber.com amberme.pdb stages=Cool,Heat,Equi,EquiMin,Prod leaprc=leaprc.protein.ff19SB \
    protons=protonation.txt watertype=opc flexwater=0 \
    refpoints=initial_restraints.pdb restraint_mult=1 \
    pdbscale=0.01 gamma_ln=1.0 barostat=1 \
    leapfile=tleap_stub.in padwater=1300 \
    cool_ns=0.001 heat_ns=0.5 equi_ns=0.5 prod_ns=0.5 \
    cool_slowdown=5 heat_slowdown=1 equi_slowdown=1 \
    restrain_omega=0 omega_weight=0 chiral_weight=0 \
    debug=0 >&! leap2amber_0.log
if(! -e Prod.rst7) then
  echo "BUILD FAILED `date`"; tail -8 leap2amber_0.log; exit 9
endif
egrep -h "parmstrip:|WARNING: cpptraj|net charge" leap2amber_0.log
set cmap = `awk '/^%FLAG CMAP_COUNT/{getline;getline;print $1;exit}' xtal.prmtop`
echo "build done: CMAP terms '$cmap' `date`"
if ( "$cmap" == "" ) then
  echo "BUILD FAILED: no CMAP terms"; mv Prod.rst7 Prod_noCMAP.rst7; exit 9
endif
