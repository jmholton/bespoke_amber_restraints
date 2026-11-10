#! /bin/tcsh -f
# water occupancies refined: example model as is, example model with the 129 deposited waters instead,
# and (deposited model, phenix_check protocol) the pure 2019 MD map
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/example_check
foreach j ( 1 2 )
  foreach model ( starthere starthere_depwat )
    ./run_exocc.csh $model flat flat $j >&! exocc_${model}_flat_j${j}.out &
    ./run_exocc.csh $model example_FPART.mtz example $j >&! exocc_${model}_example_j${j}.out &
    foreach m ( md2019 opc_10ns spce_10ns )
      ./run_exocc.csh $model ${m}_map.mtz $m $j >&! exocc_${model}_${m}_j${j}.out &
    end
  end
end
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/phenix
foreach j ( 1 2 3 4 )
  ./run_jig.csh md2019_map.mtz md2019 $j >&! jig_md2019_j${j}.out &
end
wait
