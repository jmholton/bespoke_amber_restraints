#! /bin/tcsh -f
cd /home2/jamesh/projects/amber/bulk_solvent/test_dehydrate
source compute_settings.sourceme
set path = ( /home/jamesh/projects/git/bespoke_amber_restraints $path )
dehydrate_amber_runme.com xtal.prmtop hydrated.rst7 maxreject=5 mtzfile=reference0.mtz mtzlabel=Fref notthese=none >&! dehydrate.log
echo "dehydrate exit $status"; egrep "eliminating|stripping res|ERROR" dehydrate.log
head -2 drierref.crd | cut -c1-60
echo list | cpptraj -p drier.parm7 -y drierref.crd | egrep "atoms|frames|Error"
$pmemd -O -i short.in -o short.out -p drier.parm7 -c drier.rst7 -ref drierref.crd -r short.rst7
echo "pmemd exit $status"; grep -E "Final Performance|NSTEP =  *500" short.out | head -2; ls -l short.rst7
