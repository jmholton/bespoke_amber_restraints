#! /bin/tcsh -f
# 4 jiggled starts for the maps of the earlier 20 ns spce run (test_1aho_long)
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/phenix
foreach t ( 2ns 10ns 9-10ns 20ns 19-20ns 11-20ns )
  foreach j ( 1 2 3 4 )
    ./run_jig.csh old_${t}_map.mtz old_${t} $j >&! jig_old_${t}_j${j}.out &
  end
end
wait
