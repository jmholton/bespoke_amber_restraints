#! /bin/tcsh -f
# usage: launch.csh tag   (settings from settings.txt; 2 x 1 ns production, as the water-model replicates)
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan
set L = /home2/jamesh/projects/amber/bulk_solvent/water_models/opc_r2
set tag = $1
set args = `awk -v t=$tag '$1==t{$1="";print}' $S/settings.txt`
mkdir -p $S/$tag
cd $S/$tag
cp $L/1aho.pdb $L/1aho.mtz $L/user_settings.sourceme .
cp -r $L/ligands .
ln -sf ../compute_settings.sourceme .
echo "$args" >! scan_args.txt
nohup /home2/jamesh/projects/amber/bulk_solvent/water_models/repo_scan/bulk_solvent_runme.com 1aho.pdb 1aho.mtz \
  salt=NH4,ACY salt_conc=0.68 prod_ns=1 nprod=2 max_frames=100 $args >& run.log < /dev/null &
