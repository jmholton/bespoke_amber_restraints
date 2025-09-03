#! /bin/tcsh -f
#
#  make H atoms have the same B factor as the heavier atom they are bonded to
#
#
#
set infile = current_restraints.pdb
set outfile = hsame_restraints.pdb


set tempfile = /dev/shm/${USER}/tempfile_hsame_$$_

set quiet = 0
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
      if( "$Arg" =~ *.pdb ) set infile = $Arg
    endif
    if("$key" == "debug") set debug = "1"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_hsame_
set tmpdir = `dirname $tempfile`
if("$tmpdir" != "") mkdir -p "$tmpdir"

set t = $tempfile


if(! -e "$infile") then
    set BAD = "no coordinates provided"
    goto exit
endif

# use supplied parent list if available?
echo "measuring select distances in $infile"
cat $infile |\
awk '! /^ATOM|^HETAT/{next}\
 {atom=substr($0,12,5);typ=substr($0,18,3);ch=substr($0,22,1);rn=substr($0,23,8);\
  rid=ch rn;\
  Ee=substr($0,77,2);pre=substr($0,1,21);post=substr($0,23)}\
 Ee~/ H|XP| Y/{ch="z"}\
 ! seen[atom,typ]{++printing[rid];++seen[atom,typ]}\
  printing[rid]{print pre ch post,"|"atom"|"typ"|";next}' |\
cat >! ${t}distme.pdb 
gemmi contact -d 2 --sort ${t}distme.pdb |\
awk '{dist=$NF} dist==0{next}\
   {id1=substr($0,12,17);\
    id2=substr($0,12+30,17);\
    at=id2;gsub(" ","",at)}\
 at~/^H|^EPW|^Y1/{tmp=id1;id1=id2;id2=tmp}\
  {print id1 "|" id2 "|  " dist,"DIST";}' |\
cat >! ${t}_distances.txt

cat ${t}_distances.txt |\
awk -F "|" '{id1=$1;id2=$2}\
 {at1=substr(id1,1,5);typ=substr(id1,7,3);at2=substr(id2,1,5);\
  atm1=at1;atm2=at2;gsub(" ","",atm1);gsub(" ","",atm2)}\
 atm1!~/^H|^EPW|^Y1/{next} atm2~/^H|^EPW|^Y1/{next}\
 ! seen[at1 typ]{\
   print "HPARENT",typ, at1,at2;\
   ++seen[at1 typ]}' |\
cat >! ${t}hparents.txt

cat << EOF >! ${t}xsame.txt
XSAME ASP OD[12] OD
XSAME GLU OE[12] OE
XSAME PHE|TYR CD[12] CD
XSAME PHE|TYR CE[12] CE
XSAME ASN [NO]D XD
XSAME GLN [NO]E XE
XSAME HIS [NC]D XD
XSAME HIS [NC]E XE
XSAME SO4 . O
XSAME EDO . O
XSAME MLI . O
EOF

echo "extracting parent B factors from $infile"
cat ${t}hparents.txt $infile |\
awk '/^HPARENT/{++h;htyp[h]=$2;hatom[h]=$3;hparent[h]=$4;++ishparent[$4,$2];next}\
  /^ATOM|^HETAT/{typ=substr($0,18,3);atom=substr($0,12,5);gsub(" ","",atom);\
  Bfac=substr($0,61,6);res=substr($0,21,8);\
  if(ishparent[atom,typ]){\
    for(i=1;i<=h;++i){\
      if(typ!=htyp[i])continue;\
      if(atom!=hparent[i])continue;\
      print "NEWB","|"res"|",typ,hatom[i],"  ",Bfac;\
    }\
  }\
}' >! ${t}newBs.txt
wc -l ${t}newBs.txt

echo "applying them to bound H atoms"
cat ${t}newBs.txt $infile |\
awk '/^NEWB/{res=substr($0,7,8);atom=substr($0,20,5);gsub(" ","",atom);\
    newB[res,atom]=$NF;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {atom=substr($0,12,5);gsub(" ","",atom);res=substr($0,21,8)}\
  ! newB[res,atom]{print;next}\
  {pre=substr($0,1,60);post=substr($0,67);\
   printf("%s%6.2f%s\n",pre,newB[res,atom],post)}' |\
cat >! $outfile

echo "changes:"
rmsd $outfile $infile 

ls -l $outfile

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif


