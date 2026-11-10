#!/bin/bash
# usage: bash setup_scratch.sh   (on octamus1, while /home and /home2 are up)
# Stage the water-building test (hurry_chain.csh: example model from 2 jiggled starts, then 3 altconf add_waters
# rounds) in /scratch/jamesh/projects/amber/waterbuild so that it runs from /scratch alone: the scripts are made from
# the bench's own by sed (paths only), phenix and ccp4 come from /scratch/jamesh/programs, and pick.com, waterconf.awk
# and jigglepdb.awk are copied in.  Maps: maps.txt (batch tag source) -> maps/<tag>.mtz.  Re-running it re-stages.
B=/home2/jamesh/projects/amber/bulk_solvent
WM=$B/water_models
W=/scratch/jamesh/projects/amber/waterbuild
P=/scratch/jamesh/programs
mkdir -p $W/bin $W/example_check $W/addwater $W/maps $W/tmp $W/home $W/done || exit 9
cp -p /home/jamesh/Develop/pick.com $B/bespoke_amber_restraints/waterconf.awk $B/bespoke_amber_restraints/jigglepdb.awk $W/bin/
cp -p $WM/example_check/refme.mtz $WM/example_check/starthere_j1.pdb $WM/example_check/starthere_j2.pdb $W/example_check/
cp -p $WM/addwater/close_waters.py $W/addwater/
cp -p $WM/waterbuild/maps.txt $WM/waterbuild/queue.sh $WM/waterbuild/env.sh $WM/waterbuild/check_nohome.sh $W/
sed 's#/programs/phenix-2.2-6143#'$P'/phenix-2.2-6143#g' $WM/example_check/run_exocc.csh > $W/example_check/run_exocc.csh
sed -e 's#/programs/phenix-2.2-6143#'$P'/phenix-2.2-6143#g' \
    -e 's#~/projects/git/bespoke_amber_restraints/waterconf.awk#'$W'/bin/waterconf.awk#' \
    $WM/addwater/add_waters.com > $W/addwater/add_waters.com
sed -e 's#^set A = .*#set A = '$W'/addwater#' -e 's#^set X = .*#set X = '$W'/example_check#' \
    -e '/^set path = ( \/home\/jamesh\/Develop $path )$/d' \
    -e 's#~/projects/git/bespoke_amber_restraints/jigglepdb.awk#'$W'/bin/jigglepdb.awk#' \
    $WM/addwater/rounds_close.csh > $W/addwater/rounds_close.csh
sed -e 's#^set W = .*#set W = '$W'#' -e 's#^set X = .*#set X = $W/example_check#' \
    -e 's#^  set src = $W/hurry/map_${t}.mtz$#  set src = $W/maps/${t}.mtz#' \
    -e '/^  if( $t =~ o_\* ) then$/,/^  endif$/d' -e '/^  if( $t == o_md2019 ) set src/d' \
    $WM/hurry/hurry_chain.csh > $W/hurry_chain.csh
sed 's#^X=.*#X='$W'/example_check#' $WM/hurry/rfree_table.sh > $W/rfree_table.sh
chmod +x $W/example_check/run_exocc.csh $W/addwater/add_waters.com $W/addwater/rounds_close.csh $W/hurry_chain.csh $W/bin/pick.com
# check_nohome.sh inputs: round 1 of the hl_afstd test (model, maps and the new.pdb add_waters.com made from them)
mkdir -p $W/check/ref && cp -p $WM/addwater/c1_hl_afstd/{model.pdb,maps.mtz,new.pdb} $W/check/ref/
# the maps, under their queue tags
awk '!/^#/{print $2, $3}' $W/maps.txt | while read t src; do cp -p $B/$src $W/maps/$t.mtz || echo "MISSING $t $src"; done
echo "staged $(ls $W/maps | wc -l) maps"
echo "== lines still naming /home, /programs or ~ in the scratch scripts (should be none):"
grep -n -E '/home|/programs|~/' $W/example_check/run_exocc.csh $W/addwater/add_waters.com $W/addwater/rounds_close.csh \
  $W/hurry_chain.csh $W/rfree_table.sh $W/bin/pick.com $W/addwater/close_waters.py $W/queue.sh $W/env.sh
echo "== hurry_chain.csh map source and rounds call:"; grep -n -E 'set src|rounds_close' $W/hurry_chain.csh
