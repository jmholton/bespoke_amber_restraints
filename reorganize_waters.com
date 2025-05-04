#! /bin/tcsh -f
#
#  sort waters by electron density height
#
#

set pdbfile = refmacout.pdb
set mtzfile = refmacout.mtz
set mtzlabel = ""
set mapfile = ""

set keep_numbers = 0
set super_mult = 1,1,1
set smallSG = ""
set smallmtz = ""

set debug = 0
set quiet = 0

set tempfile = /dev/shm/${USER}/reorg_$$_
mkdir -p /dev/shm/${USER}

set outfile = rewatered.pdb

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set key = `echo $Key | awk '{print tolower($0)}'`
    set val = `echo $Val | awk '{print tolower($0)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`

    if( $assign ) then
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          if(! $quiet ) echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "temp") set tempfile = "$Val"
      if("$key" == "parallel") set parallel = "$Val"
    else
      # no equal sign
      if("$key" =~ *.pdb ) set pdbfile = "$Arg"
      if("$key" =~ *.mtz ) set mtzfile = "$Arg"
      if("$key" =~ *.map ) set mapfile = "$Arg"
      if("$key" == "debug") set debug = 1
      if("$key" == "keep_numbers") set keep_numbers = 1
    endif
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_row_
if( $tempfile =~ /dev/shm/$USER/* ) mkdir -p /dev/shm/$USER/
set t = $tempfile

if(! -e "$pdbfile") then
    set BAD = "cannot find $pdbfile"
    goto exit
endif

if( "$super_mult" == "1,1,1" && ! -e "$smallmtz") then
    set smallmtz = "$mtzfile"
endif

if("$smallSG" == "" && ! -e "$smallmtz") then
    echo "checking for small-cell mtz..."
    set smallmtz = `ls -1rt *small*.mtz |& tail -n 1 | grep -v "No match"`
endif
if("$smallSG" == "" && -e "$smallmtz") then
    set smallSGnum = `echo head | mtzdump hklin $smallmtz | awk '/Space group =/{print $NF+0}' | tail -n 1`
    set smallSG = `awk -v num=$smallSGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
    echo "got $smallSG from $smallmtz"
endif

set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
set smallCELL = `echo $CELL $super_mult | awk -F "[ x,]" '{print $1/$7,$2/$8,$3/$9,$4,$5,$6}'`

if("$smallSG" == "") then
    set pdbSG = `awk '/^CRYST1/{SG=substr($0,56,14);if(length(SG)==14)while(gsub(/[^ ]$/,"",SG));print SG;exit}' $pdbfile | head -1`
    if("$pdbSG" == "R 32") set pdbSG = "R 3 2"
    if("$pdbSG" == "I 21") set pdbSG = "I 1 21 1"
    if("$pdbSG" == "A 1") set pdbSG = "P 1"
    if("$pdbSG" == "P 21 21 2 A") set pdbSG = "P 21 21 2"
    if("$pdbSG" == "P 1-") set pdbSG = "P -1"
    if("$pdbSG" == "R 3 2" && $CELL[6] == 120.00) set pdbSG = "H 3 2"
    if("$pdbSG" == "R 3" && $CELL[6] == 120.00) set pdbSG = "H 3"
    set smallSG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2 || pdbSG==$4{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
   echo "got $smallSG from $pdbfile"
endif



# separate protein from water and salt
awk '{typ=substr($0,18,3)} typ~/HOH|WAT|NH4|ACY|AMM|ACT/ || /^END/{next} {print}' $pdbfile >! ${t}dry.pdb
awk '{typ=substr($0,18,3)} typ~/NH4|ACY|AMM|ACT/{print}' $pdbfile >! ${t}salt.pdb
awk '{typ=substr($0,18,3)} typ~/HOH|WAT/{print}' $pdbfile >! ${t}oldwater.pdb

set nwaters = `grep " O " ${t}oldwater.pdb | wc -l`
if( $nwaters == 0 ) then
  echo "nothing to do"
  cp $pdbfile $outfile
  goto exit
endif

echo $smallCELL |\
awk '{printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f\n",$1,$2,$3,$4,$5,$6)}' |\
cat >! ${t}smallcell.pdb

echo "wrapping into cell: $smallCELL"
awk '{print substr($0,1,80)}' ${t}smallcell.pdb ${t}oldwater.pdb >! ${t}wrapme.pdb
wrap_into_cell.com ${t}wrapme.pdb outfile=${t}wrapped_oldwater.pdb >! ${t}wrap1.log

if(-e "$mapfile") goto gotmap


if(-e "$mtzfile") then
  set reidx = `echo $super_mult | awk -F "[ ,x]" '{print "h/"$1",k/"$2",l/"$3}'`
  echo "generating 2fofc from $mtzfile in $smallSG $smallCELL"
  reindex hklin $mtzfile hklout ${t}small.mtz << EOF >! ${t}reindex.log
  reindex $reidx
  symm $smallSG
EOF
  echo head | mtzdump hklin ${t}small.mtz >&! ${t}mtzdmp.txt
  set smallCELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' ${t}mtzdmp.txt`
  if( "$mtzlabel" == "" ) then
    awk '/Column Labels :/{getline;getline;for(i=4;i<=NF;++i)print $i}' ${t}reindex.log >! ${t}labels.txt
    set firstp = `egrep "^PH" ${t}labels.txt | egrep -v "F-obs|F-model" | head -n 1`
    set firstf = ""
    set si = 1
    while ( "$firstf" == "" && $si < 10 )
      set key = `echo $firstp $si | awk '{print substr($1,$2)}'`
      set firstf = `grep $key ${t}labels.txt | grep -v $firstp | head -n 1`
      @ si = ( $si + 1 )
    end
    set mtzlabel = "${firstf},${firstp}"
    echo "selected mtzlabel=$mtzlabel"
  endif
  set F = `echo $mtzlabel | awk -F "," '{print $1}'`
  set P = `echo $mtzlabel | awk -F "," '{print $2}'`
  echo "making map with $F $P"
  echo "labin F1=$F PHI=$P" |\
   fft hklin ${t}small.mtz mapout ${t}2fofc_small.map >! ${t}fft.log
  set mapfile = ${t}2fofc_small.map
  goto gotmap
endif


# generate fcalc density
echo "generating fcalc density in $smallSG with cell: $smallCELL"

sfall xyzin ${t}wrapped_oldwater.pdb mapout ${t}sfalled.map << EOF >! ${t}sfall.log
mode atmmap
CELL $smallCELL
SYMM $smallSG
SFSG 1
BRESET 5
EOF
set mapfile = ${t}sfalled.map
goto gotmap



gotmap:

check_map_sites.com $mapfile ${t}wrapped_oldwater.pdb |\
cat - ${t}oldwater.pdb |\
awk 'NF==1{++n;rho[n]=$1;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {++m;line[m]=$0;\
   res=substr($0,17,12);\
   rhores[res]+=rho[m]}\
  END{for(i=1;i<=m;++i){\
   res=substr(line[i],17,12);\
   print rhores[res],i,line[i];\
   }}' |\
sort -k1,1gr -k2,2g |\
awk '{print substr($0,match($0,/HETAT|ATOM/)),"    |RHO:",$1}' |\
awk -v k=$keep_numbers '{res=substr($0,17,12)}\
     res!=lastres{++n;lastres=res}\
   {pn=n} k==1{pn=res}\
   {pre=substr($0,1,22);post=substr($0,27);}\
   pn>9999{post=substr($0,27+(length(pn)-4))}\
   {printf("%s%4d%s\n",pre,pn,post)}' |\
cat >! ${t}newwater.pdb

awk '{print substr($0,1,21) "A" substr($0,23)}' ${t}salt.pdb |\
cat ${t}dry.pdb - ${t}newwater.pdb |\
cat >! ${outfile}

ls -l $outfile

exit:
if( $?BAD ) then
    echo "ERROR: $BAD"
    exit 9
endif

if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?debug && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    rm -f ${tempfile}*
endif


exit



