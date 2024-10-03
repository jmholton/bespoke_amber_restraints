#! /bin/tcsh -f
#
# check for any restraints that will blow up amber - try to correct them
#
#
set rstfile  = ""
set parmfile  = xtal.prmtop
set orignames = orignames.pdb

set refpoints = "current_restraints.pdb"
set outfile   = "new_restraints.pdb"
set outcrd    = "new_ref.crd"

set maxdist = 10
set maxkcal = 60
set kT      = 0.6
set pdbscale = 0.01
set maxbombs = 10000
set nofix   = 0
set keepweak = 0
set keepstrong = 1

set quiet = 0
set debug = 0
set tempfile = /dev/shm/$USER/tempfile_bomb_$$_

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          if(! $quiet ) echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "output") set outfile = "$Val"
      if("$key" == "outpdb") set outfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set refpoints = $Arg
      if("$Arg" =~ *.rst7 ) set rstfile = $Arg
      if("$Arg" =~ *.rst ) set rstfile = $Arg
      if("$Arg" =~ *.crd ) set rstfile = $Arg
      if("$Arg" =~ *.parmtop ) set parmfile = $Arg
      if("$Arg" =~ *.prmtop ) set parmfile = $Arg
      if("$Arg" =~ *.parm7 ) set parmfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_bomb_
set tmpdir = `dirname $tempfile`
mkdir -p $tmpdir

set t = $tempfile


if( ! $quiet ) then
cat << EOF
parmfile = $parmfile
rstfile = $rstfile
refpoints = $refpoints

maxdist = $maxdist
maxbombs = $maxbombs
maxkcal  = $maxkcal
kT       = $kT
nofix    = $nofix

outfile = $outfile
outcrd  = $outcrd

debug   = $debug
tempfile = $t
EOF
endif

if(! -e "$rstfile") then
   set BAD = "rst7 file: $rstfile does not exist"
   goto exit
endif

echo "extracting xyz from $rstfile "
rst2pdb_runme.com $rstfile pdbfile=${t}checkme.pdb orignames=$orignames Bfactors="" >! ${t}refpoint_check.log
if( $status ) then
   set BAD = "cannot read $rstfile"
   goto exit
endif

cat ${t}checkme.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  substr($0,77,2)~/ H|XP/{next}\
  {print substr($0,1,60) "  0.00" substr($0,67)}' |\
cat >! ${t}zeroB.pdb

echo "$maxdist $maxkcal $kT $pdbscale" >! ${t}params.txt

cat $refpoints |\
awk '{id=substr($0,12,15);\
  print id "|",( ! /NOT_A_BOMB/ ),( / moved/ ),"MATCHES"}' >! ${t}matches.txt

rmsd -v debug=1 $refpoints ${t}zeroB.pdb |\
cat ${t}params.txt ${t}matches.txt - |\
awk 'NR==1{maxdist=$1;maxkcal=$2;kT=$3;pdbscale=$4}\
  $NF=="MATCHES"{id=substr($0,1,15);canbomb[id]=$(NF-2);teleported[id]=$(NF-1);next}\
  ! /moved/{next}\
    {id=substr($0,1,15);d=substr($0,25,9)+0;B=substr($0,53,9)+0;w=sqrt(B*B)*pdbscale;kcal=w*d*d;\
     bigdist=(d>maxdist);\
     bigkcal=(d>maxkcal);\
     weak=(w<kT);\
     isbomb=( bigdist && bigkcal && canbomb[id] && ! teleported[id] );\
     print isbomb,bigdist,bigkcal,weak,canbomb[id],teleported[id],w,d,kcal,"|"id"| BOMB"}' |\
sort -k9gr |\
tee ${t}all_challenges.txt |\
awk '$1==1' |\
cat >! ${t}bombers.txt

echo "top challenges:"
cat ${t}all_challenges.txt |\
awk '{split($0,k,"|");id=k[2];\
     bigdist=$2;bigkcal=$3;weak=$4;canbomb=$5;teleported=$6;w=$7;d=$8;kcal=$9;\
      key = bigdist bigkcal weak canbomb teleported;whynot=""}\
  ! bigdist{whynot=whynot" bigdist"}\
  ! bigkcal{whynot=whynot" bigkcal"}\
 # ! weak{whynot=whynot" weak"}\
  ! canbomb{whynot=whynot" canbomb"}\
  teleported{whynot=whynot" teleported"}\
  ! seen[key]{++seen[key];print kcal,"kcal",d "A w="w,id, whynot}' |\
 tee ${t}worsts.txt

set worstchallenge = `head -n 1 ${t}worsts.txt`
echo "worst challenge: $worstchallenge"
set bombers = `cat ${t}bombers.txt | wc -l`
if( ! $bombers ) then
  echo "no restraint bombs detected."
  cp $refpoints ${t}.pdb
  cp ${t}.pdb $outfile 
  goto updatecrd
endif

echo "WARNING: $bombers restraint bombs"
if( $bombers && $nofix ) then
  head ${t}bombers.txt
  set BAD = "$bombers restrains are extremely strained. will blow up amber"
  goto exit
endif
if( $bombers > $maxbombs ) then
  set BAD = "too many bomb restraints. try re-picking refponts"
  goto exit
endif

# extract named atoms from MD system
cat ${t}bombers.txt ${t}checkme.pdb |\
awk '$NF=="BOMB"{split($0,k,"|");id=k[2];++sel[id];next}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next} {id=substr($0,12,15)}\
  sel[id]{print}' >! ${t}flying.pdb

# extract named atoms from restraint list
cat ${t}bombers.txt $refpoints |\
awk '$NF=="BOMB"{split($0,k,"|");id=k[2];++sel[id];next}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next} {id=substr($0,12,15)}\
  sel[id]{print}' >! ${t}wrongrest.pdb

# gather the refpoints as close as possible to flying atoms
echo "gathering ..."
rm -f ${t}newrefloc.pdb ${t}gemmi.log >& /dev/null
touch ${t}newrefloc.pdb
foreach j ( `seq 1 $bombers` )
  awk -v j=$j '/^ATOM|^HETAT/{++n} /^CRYST/ || n==j{print}' ${t}wrongrest.pdb |\
  cat >! ${t}compme.pdb
  awk '/^ATOM|^HETAT/{print $0,"THIS"}' ${t}compme.pdb |\
  cat - ${t}flying.pdb |\
  awk '! /^ATOM|^HETAT/{next} {id=substr($0,12,17)}\
     $NF=="THIS"{++sel[id];next}\
     sel[id]{print "ATOM      1  Q   XXX z9999   ",substr($0,31)}' |\
  cat >> ${t}compme.pdb
  gemmi contact -d 100 --sort ${t}compme.pdb |\
  head -n 1 |\
  tee -a ${t}gemmi.log |\
  awk '{s=$(NF-1);print "SHIFT FRAC",5-substr(s,2,1),5-substr(s,3,1),5-substr(s,4,1);exit}' |\
  pdbset xyzin ${t}compme.pdb xyzout ${t}shifted.pdb >> ${t}refpoint_check.log
  egrep "^ATOM|^HETAT" ${t}shifted.pdb | head -n 1 |\
  cat - ${t}compme.pdb |\
  awk 'NR==1{xyz=substr($0,31,24);next}\
   ! /^ATOM|^HETAT/{next}\
    {print substr($0,1,30) xyz substr($0,55);exit}' >> ${t}newrefloc.pdb
end

combine_pdbs_runme.com ${t}newrefloc.pdb $refpoints printref=1 \
  outfile=${t}newrefpoints.pdb >> ${t}refpoint_check.log

set test = `diff $refpoints ${t}newrefpoints.pdb | tee ${t}diff.txt | egrep "^>" | wc -l`
echo "$test changes between $refpoints in $outfile"
head ${t}diff.txt
echo "..."

if( $test != $bombers ) then
  echo "failed to defuse all bombs"
  if( $nofix ) then
    set BAD = "failed to defuse all bombs"
    goto exit
  endif
endif

cp ${t}newrefpoints.pdb "$outfile"


echo "checking again..."
rmsd -v debug=1 $outfile ${t}zeroB.pdb |\
cat ${t}params.txt ${t}matches.txt - |\
awk 'NR==1{maxdist=$1;maxkcal=$2;kT=$3;pdbscale=$4}\
  $NF=="MATCHES"{id=substr($0,1,15);canbomb[id]=$(NF-2);teleported[id]=$(NF-1);next}\
  ! /moved/{next}\
    {id=substr($0,1,15);d=substr($0,25,9)+0;B=substr($0,53,9)+0;w=sqrt(B*B)*pdbscale;kcal=w*d*d;\
     bigdist=(d>maxdist);\
     bigkcal=(d>maxkcal);\
     weak=(w<kT);\
     isbomb=( bigdist && bigkcal && canbomb[id] && ! teleported[id] );\
     print isbomb,bigdist,bigkcal,weak,canbomb[id],teleported[id],w,d,kcal,"|"id"| BOMB"}' |\
sort -k9gr |\
tee ${t}postgather_challenges.txt |\
awk '$1==1' |\
cat >! ${t}postgather_bombers.txt

set bombers = `cat ${t}postgather_bombers.txt | wc -l`
echo "$bombers bomb restraints remain."
if( ! $bombers ) then
   goto updatecrd
endif

cat ${t}postgather_bombers.txt |\
awk '{split($0,k,"|");id=k[2];\
     bigdist=$2;bigkcal=$3;weak=$4;canbomb=$5;teleported=$6;w=$7;d=$8;kcal=$9;\
      why="  "}\
  bigdist{why=why" bigdist"}\
  bigkcal{why=why" bigkcal"}\
  weak{why=why" weak"}\
  canbomb{why=why" canbomb"}\
  teleported{why=why" teleported"}\
  {print kcal,"kcal",d "A w="w,id, why}' |\
head

# some bond restraints remain
if( $keepweak ) then
    set BAD = "not rejecting weak bombs"
    goto exit
endif


echo "removing weak restraints below $kT kcal/A^2"
rmsd -v debug=1 $outfile ${t}zeroB.pdb |\
cat ${t}params.txt ${t}matches.txt - |\
awk 'NR==1{maxdist=$1;maxkcal=$2;kT=$3;pdbscale=$4}\
  $NF=="MATCHES"{id=substr($0,1,15);canbomb[id]=$(NF-2);teleported[id]=$(NF-1);next}\
  ! /moved/{next}\
    {id=substr($0,1,15);d=substr($0,25,9)+0;B=substr($0,53,9)+0;w=sqrt(B*B)*pdbscale;kcal=w*d*d;\
     bigdist=(d>maxdist);\
     bigkcal=(d>maxkcal);\
     weak=(w<kT);\
     isbomb=( bigdist && bigkcal && weak && canbomb[id] && ! teleported[id] );\
     print isbomb,bigdist,bigkcal,weak,canbomb[id],teleported[id],w,d,kcal,"|"id"| BOMB"}' |\
awk '$1==1' |\
sort -k9gr |\
cat >! ${t}weak_bombers.txt


cat ${t}weak_bombers.txt $outfile |\
awk '$NF=="BOMB"{split($0,k,"|");id=k[2];++sel[id];next}\
  ! /^ATOM|^HETAT/{next} {id=substr($0,12,15)}\
  sel[id]{print}' |\
head

cat ${t}weak_bombers.txt $outfile |\
awk '$NF=="BOMB"{split($0,k,"|");id=k[2];++bad[id];next}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next} {id=substr($0,12,15)}\
  ! bad[id]{print}' >! ${t}culled.pdb


set test = `diff $outfile ${t}culled.pdb | egrep "^<" | wc -l`
echo "$test new changes between $refpoints in $outfile"

cp ${t}culled.pdb $outfile



echo "checking again..."
rmsd -v debug=1 $outfile ${t}zeroB.pdb |\
cat ${t}params.txt ${t}matches.txt - |\
awk 'NR==1{maxdist=$1;maxkcal=$2;kT=$3;pdbscale=$4}\
  $NF=="MATCHES"{id=substr($0,1,15);canbomb[id]=$(NF-2);teleported[id]=$(NF-1);next}\
  ! /moved/{next}\
    {id=substr($0,1,15);d=substr($0,25,9)+0;B=substr($0,53,9)+0;w=sqrt(B*B)*pdbscale;kcal=w*d*d;\
     bigdist=(d>maxdist);\
     bigkcal=(d>maxkcal);\
     weak=(w<kT);\
     isbomb=( bigdist && bigkcal && canbomb[id] && ! teleported[id] );\
     print isbomb,bigdist,bigkcal,weak,canbomb[id],teleported[id],w,d,kcal,"|"id"| BOMB"}' |\
sort -k9gr |\
tee ${t}postweak_challenges.txt |\
awk '$1==1' |\
cat >! ${t}strong_bombers.txt
set bombers = `cat ${t}strong_bombers.txt | wc -l`
echo "$bombers bomb restraints remain."
if( ! $bombers ) then
   goto updatecrd
endif

cat ${t}strong_bombers.txt |\
awk '{split($0,k,"|");id=k[2];\
     bigdist=$2;bigkcal=$3;weak=$4;canbomb=$5;teleported=$6;w=$7;d=$8;kcal=$9;\
      why="  "}\
  bigdist{why=why" bigdist"}\
  bigkcal{why=why" bigkcal"}\
  weak{why=why" weak"}\
  canbomb{why=why" canbomb"}\
  teleported{why=why" teleported"}\
  {print kcal,"kcal",d "A w="w,id, why}' |\
head

if( $keepstrong ) then
    set BAD = "not rejecting strong bombs"
    goto exit
endif



echo "removing them."
cat ${t}strong_bombers.txt $outfile |\
awk '$NF=="BOMB"{split($0,k,"|");id=k[2];++sel[id];next}\
  ! /^ATOM|^HETAT/{next} {id=substr($0,12,15)}\
  sel[id]{print}'

cat ${t}strong_bombers.txt $outfile |\
awk '$NF=="BOMB"{split($0,k,"|");id=k[2];++bad[id];next}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next} {id=substr($0,12,15)}\
  ! bad[id]{print}' >! ${t}strongculled.pdb

set test = `diff $outfile ${t}strongculled.pdb | egrep "^<" | wc -l`
echo "$test new changes between $refpoints in $outfile"

cp ${t}strongculled.pdb $outfile






echo "checking one last time..."
rmsd -v debug=1 $outfile ${t}zeroB.pdb |\
cat ${t}params.txt ${t}matches.txt - |\
awk 'NR==1{maxdist=$1;maxkcal=$2;kT=$3;pdbscale=$4}\
  $NF=="MATCHES"{id=substr($0,1,15);canbomb[id]=$(NF-2);teleported[id]=$(NF-1);next}\
  ! /moved/{next}\
    {id=substr($0,1,15);d=substr($0,25,9)+0;B=substr($0,53,9)+0;w=sqrt(B*B)*pdbscale;kcal=w*d*d;\
     bigdist=(d>maxdist);\
     bigkcal=(d>maxkcal);\
     weak=(w<kT);\
     isbomb=( bigdist && bigkcal && canbomb[id] && ! teleported[id] );\
     print isbomb,bigdist,bigkcal,weak,canbomb[id],teleported[id],w,d,kcal,"|"id"| BOMB"}' |\
sort -k9gr |\
tee ${t}final_challenges.txt |\
awk '$1==1' |\
cat >! ${t}final_bombers.txt
set bombers = `cat ${t}final_bombers.txt | wc -l`
echo "$bombers bomb restraints remain."

cat ${t}final_bombers.txt |\
awk '{split($0,k,"|");id=k[2];\
     bigdist=$2;bigkcal=$3;weak=$4;canbomb=$5;teleported=$6;w=$7;d=$8;kcal=$9;\
     print kcal,"kcal",d "A w="w,id, whynot}' | head


if( $bombers ) then
   set BAD = "giving up"
   goto exit
endif



updatecrd:
if(! -e ref.crd) then
  echo "no ref.crd to update."
  goto exit
endif
echo "making $outcrd"
update_centroid_positions_runme.com $outfile outfile=${outcrd} orignames=${t}checkme.pdb >> ${t}refpoint_check.log
if($status || ! -e "$outcrd") then
  set BAD = "failed to update reference crd: $outcrd"
  goto exit
endif

echo "recommend:"
if("$outfile" != "$refpoints") echo "mv $outfile $refpoints"
if("$outcrd" != "ref.crd") echo "mv ${outcrd} ref.crd"

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif

exit

###########################################################################################################
#  notes and usage
#
#



