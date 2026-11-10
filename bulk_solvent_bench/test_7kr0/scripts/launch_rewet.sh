#!/bin/bash
# usage: launch_rewet.sh   (on octamus1)  start baseline.sh and the six rewet tracks, detached
cd /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/rewet/base && nohup bash /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts/baseline.sh > baseline.log 2>&1 < /dev/null &
for c in aniso iso; do
  for s in flat md2 md20; do
    map=flat; [ $s != flat ] && map=$s.mtz
    cd /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/rewet/${c}_$s && nohup bash /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts/rewet_track.sh /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/rewet/${c}_$s dry.pdb $map > track.out 2>&1 < /dev/null &
  done
done
