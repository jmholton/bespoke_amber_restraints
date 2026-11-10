#! /bin/tcsh -f
# two rounds of water filling with sensible starting occupancy/B, from the original occupancy-refined runs
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/addwater
./round.csh 3 exocc_starthere
./round.csh 4 exocc_aw3
