#! /bin/tcsh -f
# example model + waters added at the Fo-Fc peaks of each map's refinement, re-refined with the same map
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/example_check
foreach j ( 1 2 )
  ./run_exocc.csh aw1_flat flat flat $j >&! exocc_aw1_flat_j${j}.out &
  ./run_exocc.csh aw1_example example_FPART.mtz example $j >&! exocc_aw1_example_j${j}.out &
  foreach m ( md2019 opc_10ns spce_10ns )
    ./run_exocc.csh aw1_${m} ${m}_map.mtz $m $j >&! exocc_aw1_${m}_j${j}.out &
  end
end
wait
