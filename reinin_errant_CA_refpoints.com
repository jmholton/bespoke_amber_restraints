#! /bin/tcsh -f
#
#  find CA refpoints in good density (indicated by last word on each line in pdb file)
#  restrain them if they have strayed
#
set allposs = all_possible_refpoints.pdb
set flying = refme.pdb

set restraints = current_restraints.pdb
set outfile = new_restraints.pdb

set thisdir = `dirname $0`
set path = ( $thisdir $path )
set tempfile = tempfile


set errant_CA_dist = 1.5
set errant_CA_mult = 2
set kT = 0.6
set pdbscale = 0.01

set max_weight = 999.99
set min_rho = 2.0

set debug = 0

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
          echo "$Key = $Val"
          continue
      endif
      # synonyms
    else
      # no equal sign
    endif
    if("$key" == "debug") set debug = "1"
end

set t = tempfile

cat << EOF
allposs = $allposs
flying  = $flying
pdbscale = $pdbscale
kT      = $kT
min_rho = $min_rho
maxdist = $errant_CA_dist
CA_mult = $errant_CA_mult
outfile = $outfile
EOF

grep CA $allposs |\
awk -v min=$min_rho '$NF>min' |\
convert_pdb.awk -v only=protein -v CONF=" " |\
cat - $flying |\
rmsd -v debug=1 |\
grep moved |\
sort -k1.25gr >! CA_rmsd.txt

echo $errant_CA_dist |\
cat -  CA_rmsd.txt |\
awk 'NR==1{maxd=$1}\
  substr($0,25)+0>maxd{print substr($0,1,18),"| ",substr($0,25)+0}' |\
 tee bad_CAs.txt
set test = `cat bad_CAs.txt | wc -l`
if( ! $test ) then
  if( "$restraints" != "$outfile" ) cp $restraints $outfile
  echo "nothign to do"
  goto exit
endif

echo "reining in $test errant CAs "
echo "$errant_CA_mult $kT $pdbscale $errant_CA_dist PARAMS" >! ${t}_params.txt
set eCAB = `awk '{print $1*$2/$3}' ${t}_params.txt | awk '$1>999.00{$1=999.99} {print $1}'`
grep CA $allposs |\
convert_pdb.awk -v CONF=" " -v BFAC=$eCAB |\
  awk '! /^ATOM|^HETAT/{next} {id=substr($0,12,17)}\
  ! seen[id]{print substr($0,1,80) "            NOT_A_BOMB";\
  ++seen[id]}' >! ${t}CA_allposs.pdb
combine_pdbs_runme.com $restraints saveXYZ=1 xor=1 printref=1 ${t}CA_allposs.pdb outfile=${t}new_CA.pdb > /dev/null

awk -F "|" '{print $0,"| BAD"}' bad_CAs.txt |\
cat ${t}_params.txt - $restraints |\
awk 'NR==1{mult=$1;kT=$2;s=$3;dmin=$4;maxB=999.99}\
 $NF=="BAD"{split($0,w,"|");id=substr($0,1,17);++sel[id];d[id]=w[2];next}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,12,17);B=substr($0,61,6);pre=substr($0,1,60);post=substr($0,67,12)}\
  {B=mult*B}\
  {B=B+kT/s*(1+d[id]-dmin)**2}\
  B>maxB{B=maxB}\
  sel[id]{#print "REMARK",B,kT,s,id,d[id],dmin;\
    printf("%s%6.2f%s         NOT_A_BOMB\n",pre,B,post)}' >! ${t}tighter_CA.pdb
combine_pdbs_runme.com $restraints ${t}tighter_CA.pdb ${t}new_CA.pdb $flying outfile=${t}new_refpoints.pdb > /dev/null
#  rmsd $restraints clip_CA.pdb | egrep "MAXD.Bfac"

cp ${t}new_refpoints.pdb $outfile

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

