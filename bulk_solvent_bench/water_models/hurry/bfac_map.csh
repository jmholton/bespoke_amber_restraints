#! /bin/tcsh -f
# usage: bfac_map.csh tag   (on octamus1)  solvent map of the iterations sets.txt lists for tag, each rendered by
# nc2mtz_gemmi.com with the run's own per-atom B (bfac_iter.csh), averaged over frames -> hurry/map_<tag>.mtz.
# Iterations without a B sidecar (the run's last) are skipped and listed.
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
set tag = $1
set specs = `awk -v t=$tag '$1==t{$1="";print}' $D/sets.txt`
if( "$specs" == "" ) then
  echo "${tag}: not in sets.txt"
  exit 9
endif
mkdir -p $D/$tag
set jobs = 0
foreach s ( $specs )
  set p = ( `echo $s | awk -F: '{print $1,$2,$3,$4}'` )
  foreach i ( `seq $p[2] $p[4] $p[3]` )
    $D/bfac_iter.csh $p[1] $i $D/$tag >& $D/$tag/iter_$p[1]_${i}.log &
    @ jobs++
    if( $jobs % 24 == 0 ) wait
  end
end
wait
cat $D/$tag/iter_*.log | grep -v "sidecar ok"
set dirs = `ls $D/$tag/*/nframes | xargs -n1 dirname`
echo "${tag}: $#dirs iterations rendered"
/programs/phenix-2.2-6143/bin/phenix.python $D/avg_maps.py $D/map_${tag}.mtz $dirs
