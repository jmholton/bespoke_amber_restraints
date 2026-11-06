#!/bin/bash
# usage: refine_long.sh   (on octamus1)  when the 20 ns run (ensemble_opc_long) writes its map, refine the ensemble
# with PRISM's waters against it, seeds 1 and 2, the same test as the 2 ns map (refine/run_ref.csh)
B=/home2/jamesh/projects/amber/bulk_solvent/test_7kr0
L=$B/ensemble_opc_long
until [ -e $L/bulk_solvent.mtz ]; do
  grep -q ERROR $L/run.log && { echo "20 ns run failed: $(grep ERROR $L/run.log | head -1)"; exit 9; }
  sleep 60
done
cp $L/bulk_solvent.mtz $B/refine/md_20ns.mtz
cd $B/refine
./run_ref.csh ens md_20ns.mtz 20ns 1 > 20ns_s1.out 2>&1 &
./run_ref.csh ens md_20ns.mtz 20ns 2 > 20ns_s2.out 2>&1 &
wait
cat 20ns_s1.out 20ns_s2.out
