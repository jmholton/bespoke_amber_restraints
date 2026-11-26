#! /bin/tcsh -f
# Build one 0.68 M charge-test system (setup_s68_builds.sh made amberme.pdb, protonation.txt and tleap_stub.in):
# all-atom restraints at the nearest centroid at k 1 (run_local.csh's build, verbatim), then the same leap2amber call
# as every 2s2r build.  Code: code_noh (as build_2s2r_amide); software from programs.sourceme, MD through the SITE
# compute settings (SLURM).   usage: build_s68.csh build_dir
set d = $1
source /home2/jamesh/projects/amber/bulk_solvent/programs.sourceme
setenv OMP_NUM_THREADS 4
set path = ( /scratch/jamesh/projects/amber/1aho_hurry_test/code_noh $path )
cd $d
echo "build $d on `hostname -s`: `tr '\n' ';' < protonation.txt` `date`"
cp ../centroids/centroids_in_density.pdb all_possible_refpoints.pdb
centroids_nearby_runme.com amberme.pdb reffile=all_possible_refpoints.pdb \
  softener=2 weight=Bfac maxdist=1 hohscale=1 \
  outfile=initial_restraints.pdb debug=0 >&! c2r_0.log
if(! -e initial_restraints.pdb) then
  echo "BUILD FAILED at centroids_nearby"; tail -8 c2r_0.log; exit 9
endif
mv initial_restraints.pdb c2r_restraints.pdb
awk '/^ATOM|^HETAT/{printf("%s%6.2f%s\n",substr($0,1,60),100,substr($0,67));next} {print}' \
  c2r_restraints.pdb >! initial_restraints.pdb
echo "restraints: `egrep -c '^ATOM|^HETAT' initial_restraints.pdb`"
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
  echo "BUILD FAILED at leap2amber `date`"; tail -8 leap2amber_0.log; exit 9
endif
egrep -h "parmstrip:|WARNING: cpptraj|net charge" leap2amber_0.log
set cmap = `awk '/^%FLAG CMAP_COUNT/{getline;getline;print $1;exit}' xtal.prmtop`
echo "build done: CMAP terms '$cmap' `date`"
if ( "$cmap" == "" ) then
  echo "BUILD FAILED: no CMAP terms"; mv Prod.rst7 Prod_noCMAP.rst7; exit 9
endif
