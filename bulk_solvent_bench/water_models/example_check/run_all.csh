#! /bin/tcsh -f
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/example_check
foreach j ( 1 2 3 4 )
  ./run_ex.csh flat flat $j >&! ex_flat_j${j}.out &
  ./run_ex.csh example_map.mtz example $j Fpart PHIpart >&! ex_example_j${j}.out &
  foreach m ( opc_10ns tip3p_10ns opc3_10ns tip4pew_10ns fb4_10ns spce_10ns opc_r2_2ns spceold_2ns spceold_20ns )
    ./run_ex.csh ${m}_map.mtz $m $j >&! ex_${m}_j${j}.out &
  end
end
wait
