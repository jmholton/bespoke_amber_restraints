#! /bin/tcsh -f
# w1 (best example R-free) continued: same directory contents, production extended to 100 x 1 ns
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan
cd $S/w1_100ns
nohup /home/jamesh/projects/git/bespoke_amber_restraints/bulk_solvent_runme.com 1aho.pdb 1aho.mtz \
  salt=NH4,ACY salt_conc=0.68 prod_ns=1 nprod=100 max_frames=500 watertype=opc >& run.log < /dev/null &
