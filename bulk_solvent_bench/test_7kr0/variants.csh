#! /bin/tcsh -f
# usage: variants.csh   (on octamus1)  smaller restraint weights and a tighter altloc merge, side by side
cd /home2/jamesh/projects/amber/bulk_solvent/test_7kr0
./one_variant.csh wps0.3 weight_per_sigma=0.3 &
./one_variant.csh wps0.1 weight_per_sigma=0.1 &
./one_variant.csh wps0.03 weight_per_sigma=0.03 &
./one_variant.csh merge0.4 altloc_merge=0.4 &
wait
