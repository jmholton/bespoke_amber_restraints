#! /bin/tcsh -f
#
#  faster than tleap: make new ref.crd file to reflect updated centroids
#
#  needs:
#  pdb with new centroid coordinates
#  topology file
#  flying.pdb for current positions that matches names in provided pdb
#
#
set pdbfile    = restraints_to_use.pdb
set flying     = orignames.pdb
set topfile    = xtal.prmtop
set outfile    = newref.crd
set minwt      = 0
set pdbscale   = 0.01

# flag to debug things
set debug = 0
set quiet = 0

set tempfile = /dev/shm/${USER}/tempfile_ucp_$$_


foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print tolower($1)}'`
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
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "prmtop") set topfile = "$Val"
      if("$key" == "parmfile") set topfile = "$Val"
      if("$key" == "output") set outfile = "$Val"
      if("$key" == "orignames") set flying = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $debug && $tempfile =~ /dev/shm/* )  set tempfile = ./tempfile_ucp_


if(! -e $flying && -e Bfac.pdb ) then 
    set flying = Bfac.pdb
endif

cat << EOF
pdbfile    = $pdbfile
flying     = $flying
topfile    = $topfile
outfile    = $outfile

tempfile   = $tempfile
debug      = $debug
EOF

set t = $tempfile

set goalatoms = `echo list | cpptraj -p $topfile | awk '$4=="atoms,"{print $3}' | head -n 1`
echo "$goalatoms atoms in $topfile"

cat $flying |\
awk -v ga=$goalatoms '! /^ATOM|^HETAT/{next}\
  n+0<ga{++n;print}' |\
tee ${t}atomnames.pdb | wc -l >! ${t}atomcount.txt
set nnames = `cat ${t}atomcount.txt`
if( $nnames != $goalatoms ) then
  set BAD = "not enough atoms in $flying : $nnames need $goalatoms"
  goto exit
endif

# filter out too-small weights
set minB = `echo $minwt $pdbscale | awk '{print $1/$2}'`
echo "filtering out centroid atoms with B < $minB"
echo $minB | cat - $pdbfile |\
awk 'NR==1{minB=$1;next} ! /^ATOM|^HETAT/{print;next}\
  {B=substr($0,61,6)+0}\
  B>=minB{print}' |\
cat >! ${t}active_centroids.pdb 

echo "overlaying xyz of active centroids from $pdbfile"
combine_pdbs_runme.com printref=1 ${t}active_centroids.pdb ${t}atomnames.pdb \
  outfile=${t}newref.pdb >! ${t}combref.log


echo | cpptraj -p $topfile -y ${t}newref.pdb -x ${t}.rst7  >! ${t}combref.log
if( $status || ! -e ${t}.rst7 ) then
  # maybe it was too short?
  echo "WARNING: conversion failed. adding dummy atoms at the end..."
  set goalatoms = `echo list | cpptraj -p $topfile | awk '$4=="atoms,"{print $3}' | head -n 1`
  set pdbatoms = `egrep "^ATOM|^HETAT" ${t}newref.pdb | wc -l`
  echo $goalatoms |\
  cat - ${t}newref.pdb |\
  awk 'NR==1{goal=$1;next}\
    /^CRYST/{print} ! /^ATOM|^HETAT/{next}\
    {++n}\
    n<=goal{print}\
    END{while(n<goal){++n;\
      print}}' |\
  cat >! ${t}newref2.pdb

  echo | cpptraj -p $topfile -y ${t}newref2.pdb -x ${t}.rst7 >! ${t}combref2.log
  if( $status ) then
    set BAD = "conversion failed"
  endif
endif

cp ${t}.rst7 $outfile
echo "updated $outfile"

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

