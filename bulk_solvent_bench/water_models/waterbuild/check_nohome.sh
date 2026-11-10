#!/bin/bash
# usage: bash check_nohome.sh   (on octamus1, after setup_scratch.sh and the phenix copy)
# Run the scratch add_waters.com on the inputs of an earlier round (hl_afstd round 1) and a phenix import under strace:
# the new waters must match the original run's, and no file access may leave /scratch and the local system.
W=/scratch/jamesh/projects/amber/waterbuild
source $W/env.sh
cd $W/check || exit 9
cp -p ref/model.pdb ref/maps.mtz .
strace -f -e trace=file -o trace_aw.txt $W/addwater/add_waters.com model.pdb maps.mtz occ=0.5 B=15 heavy dist=2.3 altconf > aw.log 2>&1
grep altconf: aw.log
cmp new.pdb ref/new.pdb && echo "add_waters.com from /scratch: new.pdb identical to the original run's"
strace -f -e trace=file -o trace_px.txt $P/phenix-2.2-6143/bin/phenix.python -c "import libtbx.load_env, mmtbx.f_model, iotbx.pdb, cctbx.miller; print('phenix imports ok', libtbx.env.build_path)"
echo "== file accesses outside /scratch and the local system (count, path):"
grep -h -o -E '"/(home|home2|programs|petabyte)[^"]*"' trace_aw.txt trace_px.txt | sort | uniq -c | sort -rn | head -20
echo "== end of list"
