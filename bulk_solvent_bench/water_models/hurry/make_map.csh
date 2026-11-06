#! /bin/tcsh -f
# usage: make_map.csh tag   (on octamus1)  solvent map of the frames sets.txt lists for tag -> hurry/map_<tag>.mtz
# (extract_frames.csh per run, render_frames.py at uniform B, 10 unless the line has B=<value>; frames translated onto
# the run's alignment reference first if the line has align); the frame PDBs are removed afterwards
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/hurry
set tag = $1
set specs = `awk -v t=$tag '$1==t{for(i=2;i<=NF;++i)if($i!~/^B=/&&$i!="align")print $i}' $D/sets.txt`
set B = `awk -v t=$tag '$1==t{for(i=2;i<=NF;++i)if($i~/^B=/)print substr($i,3)}' $D/sets.txt`
if( "$B" == "" ) set B = 10
set align = `awk -v t=$tag '$1==t{for(i=2;i<=NF;++i)if($i=="align")print $i}' $D/sets.txt`
if( "$specs" == "" ) then
  echo "${tag}: not in sets.txt"
  exit 9
endif
mkdir -p $D/$tag/frames
foreach s ( $specs )
  set p = ( `echo $s | awk -F: '{print $1,$2,$3,$4,$5}'` )
  $D/extract_frames.csh $p[1] $p[2] $p[3] $p[4] $p[5] $D/$tag/frames $align >& $D/$tag/extract_$p[1].log &
end
wait
cat $D/$tag/extract_*.log
ls $D/$tag/frames/*.pdb >! $D/$tag/frames.txt
/programs/phenix-2.2-6143/bin/phenix.python $D/render_frames.py $D/map_${tag}.mtz $D/$tag/frames.txt B=$B nproc=48
if( -e $D/map_${tag}.mtz ) rm -f $D/$tag/frames/*.pdb
