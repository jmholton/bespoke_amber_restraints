#! /bin/tcsh -f
#
#  faster than tleap: edit the ref.crd file to reflect updated centroids
#
#  needs:
#  pdb with new coordinates
#  previous ref.crd
#  topology file
#  orignames.pdb that matches names in provided pdb
#
#
set pdbfile    = restraints_to_use.pdb
set orignames  = orignames.pdb
set topfile    = xtal.prmtop
set outfile    = newref.crd

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
      if("$key" == "output") set outfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $debug && $tempfile =~ /dev/shm/* )  set tempfile = ./tempfile_ucp_


if(! -e $orignames && -e Bfac.pdb ) then 
    set orignames = Bfac.pdb
endif

cat << EOF
pdbfile    = $pdbfile
orignames  = $orignames
topfile    = $topfile
outfile    = $outfile

tempfile   = $tempfile
debug      = $debug
EOF

set t = $tempfile

set goalatoms = `echo list | cpptraj -p $topfile | awk '$4=="atoms,"{print $3}' | head -n 1`
echo "$goalatoms atoms in $topfile"

cat $orignames |\
awk -v ga=$goalatoms '! /^ATOM|^HETAT/{next}\
  n+0<ga{++n;print}' |\
tee ${t}atomnames.pdb | wc -l >! ${t}atomcount.txt
set nnames = `cat ${t}atomcount.txt`
if( $nnames != $goalatoms ) then
  set BAD = "not enough atoms in $orignames : $nnames need $goalatoms"
  goto exit
endif

echo "overlaying xyz from $pdbfile"
combine_pdbs_runme.com printref=1 $pdbfile ${t}atomnames.pdb \
  outfile=${t}newref.pdb >! ${t}combref.log


echo | cpptraj -p $topfile -y ${t}newref.pdb -x ${t}.rst7  >! ${t}combref.log
if( $status || ! -e ${t}.rst7 ) then
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

