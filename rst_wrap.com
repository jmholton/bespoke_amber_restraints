#! /bin/tcsh -f
#
# wrap back errant residues in rst7 file that are not near restraint reference points
#
#
set rstfile  = ""
set parmfile  = xtal.prmtop
set orignames = orignames.pdb
set outprefix = "wrapped"
set refpoints = "all_possible_refpoints.pdb current_restraints.pdb"

set maxdist = 6.1

set quiet = 0
set debug = 0
set tempfile = /dev/shm/$USER/tempfile_r2w_$$_

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
      if("$key" == "output") set outprefix = "$Val"
      if("$key" == "outfile") set outprefix = "$Val"
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

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_rw_
set tmpdir = `dirname $tempfile`
mkdir -p $tmpdir

set t = $tempfile


if("$outprefix" =~ *.rst7) set outprefix = `basename $outprefix .rst7`
if("$outprefix" == "") set outprefix = wrapped

if(! $quiet ) then
cat << EOF
parmfile = $parmfile
rstfile = $rstfile
refpoints = $refpoints

maxdist = $maxdist

outfile = ${outprefix}.rst7

debug   = $debug
tempfile = $t
EOF
endif

echo "extracting xyz from $rstfile "
rst2pdb_runme.com $rstfile pdbfile=${t}orig.pdb >! ${t}cpptraj_extract.log
set test = `grep "new parmfile: resized.parm7" ${t}cpptraj_extract.log | wc -l`
if( $test ) then
   echo "using resized.parm7 from now on"
   set parmfile = resized.parm7
endif

# expand the reference points a little
echo "padding reference points"
cat $refpoints |\
awk -v d=$maxdist '/^CRYST/{a=$2;b=$3;c=$4} ! /^ATOM|^HETAT/{next}\
  {X=substr($0,31,8)+0;Y=substr($0,39,8)+0;Z=substr($0,47,8)+0;\
   pre=substr($0,1,30);post=substr($0,55)}\
  seen[X,Y,Z]{next} {++seen[X,Y,Z]}\
  X>0 && X<d{printf("%s%8.3f%8.3f%8.3f%s\n",pre,X+a,Y,Z,post)}\
  Y>0 && Y<d{printf("%s%8.3f%8.3f%8.3f%s\n",pre,X,Y+b,Z,post)}\
  Z>0 && Z<d{printf("%s%8.3f%8.3f%8.3f%s\n",pre,X,Y,Z+c,post)}\
  X<a && X>a-d{printf("%s%8.3f%8.3f%8.3f%s\n",pre,X-a,Y,Z,post)}\
  Y<b && Y>b-d{printf("%s%8.3f%8.3f%8.3f%s\n",pre,X,Y-b,Z,post)}\
  Z<c && Z>c-d{printf("%s%8.3f%8.3f%8.3f%s\n",pre,X,Y,Z-c,post)}\
  ' |\
cat >! ${t}extrapoints.pdb


echo "extracting non-protein atoms"
convert_pdb.awk -v skip=protein,H ${t}orig.pdb >! ${t}fluff.pdb

egrep "^CRYST" ${t}fluff.pdb | head -n 1 >! ${t}compme.pdb
egrep "^ATOM|^HETAT" ${t}fluff.pdb >> ${t}compme.pdb 
cat $refpoints ${t}extrapoints.pdb |\
awk '! /^ATOM|^HETAT/{next} {xyz=substr($0,31,24)}\
   seen[xyz]{next}\
  {++seen[xyz];\
   printf("ATOM      1  O   HOH z%8d%s\n",++n,substr($0,31))}' |\
 hy36_encode.awk >> ${t}compme.pdb
echo "looking for those within $maxdist A of a potential reference point"
gemmi contact --nosym -d $maxdist --sort --ignore=3 ${t}compme.pdb >! ${t}gemmi.log
grep "HOH z" ${t}gemmi.log |\
awk 'gsub("HOH z","")>1{next}\
  {if(debug>1)print $0,"DEBUG IN";\
   id1=substr($0,12,17);\
   print id1 "|  " $NF+0,"DIST";}' |\
cat - ${t}fluff.pdb |\
awk '$NF=="DIST"{id=substr($0,1,17);++sel[id]}\
 ! /^ATOM|^HETAT/{next}\
   {id=substr($0,12,17)}\
  ! sel[id]{print $NF}' |\
cat >! ${t}dowrap.txt
set num = `cat ${t}dowrap.txt | wc -l`
echo "$num molecules can be wrapped"
if( $num == 0 ) then
   echo "nothing to do"
   cp $rstfile ${outprefix}.rst7
   goto exit
endif
cat ${t}dowrap.txt |\
 sort -u | sort -g |\
 awk 'NR==1{s=e=$1;next} $1==e+1{e=$1;next} {print s"-"e;s=e=$1} END{print s"-"e}' |\
 awk -F "-" '$1==$2{$0=$1} {print}' |\
 sort -u | sort -g >! ${t}ranges.txt 
set ranges = `cat ${t}ranges.txt `
set ranges = `echo $ranges | awk '{gsub(" ",",");print}' `

echo "imaging selected ranges"
cpptraj -p $parmfile -y $rstfile << EOF >! ${t}cpptraj_wrap.log
image center :$ranges
trajout ${outprefix}.rst7
EOF
#rst2pdb_runme.com wrapped.rst7 > /dev/null ; rmsd ref.pdb wrapped.pdb 

rst2pdb_runme.com ${outprefix}.rst7 pdbfile=${t}wrapped.pdb >! ${t}rst2pdb.log
rmsd -v debug=1 ${t}wrapped.pdb ${t}orig.pdb |\
awk '{id=substr($0,10,7)} /moved/ && substr($0,25,7)+0!=0 && ! seen[id]{\
  ++seen[id];print substr($0,1,34)}' >! wraps.txt
wc -l wraps.txt

rm -f ${t}rmsd.dat
cpptraj -p $parmfile -y $rstfile -y ${outprefix}.rst7 << EOF >> cpptraj_wrap.log
rmsd out ${t}rmsd.dat
EOF
tail -n 1 ${t}rmsd.dat | awk '{print "rmsd_wrap:",$NF}'

ls -l ${outprefix}.rst7

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


