#! /bin/tcsh -f
# same atoms, same volume (barostat on with taup=1e11, as the 2019 runs, so pressure is computed), restraints on vs off
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/pressure_test
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
foreach t ( restrained4 unrestrained4 )
  srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP -O -i ${t}.in -o ${t}.out -p /home2/jamesh/projects/amber/bulk_solvent/water_models/spce/hydrated.prmtop     -c /home2/jamesh/projects/amber/bulk_solvent/water_models/spce/prod_10.rst7 -ref /home2/jamesh/projects/amber/bulk_solvent/water_models/spce/hydrated_ref.crd -r ${t}.rst7 -inf ${t}.mdinfo &
end
wait
echo finished4
