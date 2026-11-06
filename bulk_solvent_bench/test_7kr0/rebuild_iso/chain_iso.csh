#! /bin/tcsh -f
# usage: chain_iso.csh   (on octamus1)  PRISM 7KR0 S0112 ensemble with its ANISOU records and waters deleted, rebuilt against flat and MD
# bulk solvent: round 0 refines the bare ensemble; rounds 1-3 add mFo-DFc peaks with ens_waters.py (conformer-aware,
# held occupancies) to the previous round's seed-1 model and refine (run_ref.csh: PRISM's standing recipe, 5 cycles).
# Round 3 also refines from seed 2.
set S = /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/scripts
cd /home2/jamesh/projects/amber/bulk_solvent/test_7kr0/rebuild_iso
foreach sol ( flat md )
  set map = flat
  if( $sol == md ) set map = md_opc.mtz
  cp nowat_iso.pdb w0_${sol}.pdb
  ./run_ref.csh w0_${sol} $map $sol 1 >&! w0_${sol}.out &
  ./run_ref.csh ens_iso $map $sol 1 >&! ens_iso_${sol}.out &
end
wait
echo "round 0 done"
foreach n ( 1 2 3 )
  set p = `expr $n - 1`
  foreach sol ( flat md )
    set map = flat
    if( $sol == md ) set map = md_opc.mtz
    /programs/phenix-2.2-6143/bin/phenix.python $S/ens_waters.py ref_w${p}_${sol}_${sol}_s1_001.pdb \
      ref_w${p}_${sol}_${sol}_s1_001.mtz w${n}_${sol}.pdb sigma=3.0 >&! build_w${n}_${sol}.log
    echo "$sol round ${n}: `tail -2 build_w${n}_${sol}.log | tr '\n' ' '`"
    ./run_ref.csh w${n}_${sol} $map $sol 1 >&! w${n}_${sol}_s1.out &
    if( $n == 3 ) ./run_ref.csh w${n}_${sol} $map $sol 2 >&! w${n}_${sol}_s2.out &
  end
  wait
  echo "round $n done"
end
grep Final w*.out
