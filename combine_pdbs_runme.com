#! /bin/tcsh -f
#
#  filter/mix/cut/paste two PDB files into one
#
#
set pdbs = ""
set refpdb = ""
set outfile = "new.pdb"

set printref = 0
set xor = 0
set newocc = ""
set newB =   ""

set maxocc = 1
set minocc = 0
set maxB   = 999.99
set minB   = 0

set saveXYZ = 0
set saveOcc = 0
set saveBfac = 0
set savesuff = 0

foreach arg ( $* )
    set Key = `echo $arg | awk -F "=" '{print $1}'`
    set key = `echo $arg | awk -F "=" '{print tolower($1)}'`
    set val = `echo $arg | awk -F "=" '{print $2}'`

    if("$key" =~ *.pdb && "$val" == "") then
        set pdbs = ( $pdbs "$Key" )
        set refpdb = $Key
    endif
    if("$key" == "refpdb") set refpdb = "$val"
    if("$key" == "outfile") set outfile = "$val"
    if("$key" == "outpdb") set outfile = "$val"
    if("$key" == "printref") set printref = "$val"
    if("$key" == "keepref") set printref = "$val"
    if("$key" == "xor") set xor = "$val"
    if("$key" == "occ") set newocc = "$val"
    if("$Key" == "B") set newB = "$val"
    if("$key" == "minocc") set minocc = "$val"
    if("$key" == "maxocc") set maxocc = "$val"
    if("$Key" == "minB") set minB = "$val"
    if("$Key" == "maxB") set maxB = "$val"
    if("$key" =~ savex*) set saveXYZ = 1
    if("$key" =~ saveo*) set saveOcc = 1
    if("$key" =~ saveb*) set saveBfac = 1
    if("$key" =~ savesuf*) set savesuff = 1
end
if("$pdbs[$#pdbs]" == "$refpdb") then
   set pdbs[$#pdbs] = ""
endif
if("$pdbs" == "") then
   set name = `basename $0`
   cat << EOF
usage: $name *.pdb ref.pdb [printref=1] [xor=0] [occ=1] [B=5] [saveXYZ] [saveOcc] [saveB] [minocc=0] [maxocc=1] [minB=0] [maxB=999]  [outfile=${outfile}]"
printref - print all atoms in ref.pdb, modified or not
xor      - only print atoms not listed in *.pdb
occ      - reset all printed values
minocc/maxocc/minB/maxB - set clipping limits
saveXYZ  - use XYZ from ref.pdb
saveOcc  - use occ from ref.pdb
saveB    - use B from ref.pdb
EOF
    exit 9
endif

if( $xor ) set printref = 1


cat << EOF
splicing $pdbs into $refpdb
printref = $printref
xor = $xor
newocc = $newocc
newB = $newB
saveXYZ = $saveXYZ
saveOcc = $saveOcc
saveBfac = $saveBfac
savesuff = $savesuff

outfile  = $outfile
EOF

echo "$printref $xor ${newocc}O ${newB}B  $saveXYZ $saveOcc $saveBfac $savesuff $minocc $maxocc $minB $maxB" |\
    cat - $pdbs |\
    awk 'NR==1{print;next}\
      /^ATOM|^HETAT/{\
        print "NEW",substr($0,5)}' |\
    cat - $refpdb |\
    awk 'NR==1{printref=$1;xxor=$2;newocc=$3;newB=$4;\
               saveXYZ=$5;saveOcc=$6;saveBfac=$7;savesuff=$8;\
               minocc=$9;maxocc=$10;minB=$11;maxB=$12;next}\
    {ID=substr($0,12,19);conf=substr($0,17,1);ncID=substr($0,12,5)" "substr($0,18,10)} \
    /^NEW/{newline[ID]=substr($0,5);next}\
    /^CRYST/{if(! cryst)print;++cryst;next}\
    ! /^ATOM|^HETAT/{print;next}\
    newline[ID] && xxor{next}\
    {oldline=$0}\
    newline[ID]{$0="ATOM" newline[ID]}\
    {gsub("^ATOMTM","ATOM  ")}\
    saveXYZ{$0=substr($0,1,30) substr(oldline,31,24) substr($0,55)}\
    saveOcc{$0=substr($0,1,54) substr(oldline,55,6) substr($0,61)}\
    saveBfac{$0=substr($0,1,60) substr(oldline,61,6) substr($0,67)}\
    savesuff{$0=substr($0,1,66) substr(oldline,67)}\
    {pre=substr($0,1,54);occ=substr($0,55,6)+0;B=substr($0,61,6)+0;post=substr($0,67)}\
    occ<minocc{occ=minocc} occ>maxocc{occ=maxocc} B<minB{B=minB} B>maxB{B=maxB}\
    newocc!="O"{occ=newocc}\
      newB!="B"{B=newB}\
    {$0=sprintf("%s%6.2f%6.2f%s",pre,occ,B,post)}\
    newline[ID] || printref{print}' |\
    cat >! $outfile
exit

rmsd $refpdb $outfile

