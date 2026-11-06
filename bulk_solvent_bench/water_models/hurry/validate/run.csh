#! /bin/tcsh -f
# re-render scan/n1's solvent frames with render_frames.py, to compare with n1/bulk_solvent.mtz (nc2mtz_gemmi.com)
set V = /home2/jamesh/projects/amber/bulk_solvent/water_models/hurry/validate
set N = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan/n1
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
mkdir -p $V/frames
cd $V
cat << EOT >! x.in
parm $N/solvent.parm7
trajin $N/solvent.nc
strip :WAT,HOH@Y1,EPW
trajout frames/n1.pdb pdb multi pdbv3 keepext sg "P 1"
go
EOT
cpptraj < x.in >& x.log
ls $V/frames/n1.*.pdb >! frames.txt
/programs/phenix-2.2-6143/bin/phenix.python ../render_frames.py $V/n1_rerender.mtz frames.txt B=10
