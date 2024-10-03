#! /bin/tcsh -f
#
#  select atoms that are x-ray equivalent and make them have same B factor
#
#
#
set infile = current_restraints.pdb
set outfile = xsame_restraints.pdb


set tempfile = /dev/shm/${USER}/tempfile_xsame_$$_

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

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_xsame_
set tmpdir = `dirname $tempfile`
if("$tmpdir" != "") mkdir -p "$tmpdir"

set t = $tempfile


if(! -e "$infile") then
    set BAD = "no coordinates provided"
    goto exit
endif

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

echo "examining $infile"
cat ${t}xsame.txt $infile |\
awk '/^XSAME/{++x;xtyp[x]=$2;xatom[x]=$3;xdgn[x]=$4;next}\
  /^ATOM|^HETAT/{typ=substr($0,18,3);atom=substr($0,12,5);gsub(" ","",atom);\
  Bfac=substr($0,61,6);res=substr($0,21,8);\
  for(x in xtyp){\
    if(typ ~ xtyp[x] && atom ~ xatom[x]){\
      id=res" "typ" "xdgn[x];\
      ++count[id];sum[id]+=Bfac;\
      atoms[id]=atoms[id]" "atom;\
    }\
  }\
}\
END{for(id in sum){\
  res=substr(id,1,8);\
  split(atoms[id],a);\
  for( atm in a) {\
    print "NEWB",res,a[atm],"  ",sum[id]/count[id];\
  }\
 }\
}' >! ${t}newBs.txt
wc -l ${t}newBs.txt

echo "applying"
cat ${t}newBs.txt $infile |\
awk '/^NEWB/{res=substr($0,6,8);atom=substr($0,13,5);gsub(" ","",atom);\
    newB[res,atom]=$NF;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {atom=substr($0,12,5);gsub(" ","",atom);res=substr($0,21,8)}\
  ! newB[res,atom]{print;next}\
  {pre=substr($0,1,60);post=substr($0,67);\
   printf("%s%6.2f%s\n",pre,newB[res,atom],post)}' |\
cat >! $outfile

echo "changes:"
rmsd $infile $outfile

ls -l $outfile

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif


