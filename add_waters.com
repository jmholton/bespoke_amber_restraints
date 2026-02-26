#! /bin/tcsh -f
#
#        pick top difference peak far from any other and make it water
#
#   -James Holton   3-31-18
#
#
set tempfile = ${CCP4_SCR}/addwater$$
set tempfile = tempfile

set sigma = 4
set dist  = 1.8A

set pdbfile = refmacout.pdb
set mtzfile = refmacout.mtz
set mapfile = ""
set maskfile = ""

set top
set rmsrho
foreach arg ( $* )
    if("$arg" =~ sigma=*) then
        set user_sigma = `echo $arg | awk -F "=" '{print $NF+0}'`
        continue
    endif
    if("$arg" =~ chain=*) then
        set user_chain = `echo $arg | awk -F "=" '{print $NF}'`
        continue
    endif
    if("$arg" =~ occ=*) then
        set user_occ = `echo $arg | awk -F "=" '{print $NF}'`
        continue
    endif
    if("$arg" =~ B=*) then
        set user_B = `echo $arg | awk -F "=" '{print $NF}'`
        continue
    endif
    if("$arg" =~ mask=*map) then
       set maskfile = `echo $arg | awk -F "=" '{print $NF}'`
       if(! -e "$maskfile") echo "WARNING: $maskfile does not exist! "
       continue
    endif
    if("$arg" =~ map=*map) then
       set mapfile = `echo $arg | awk -F "=" '{print $NF}'`
       if(! -e "$mapfile") echo "WARNING: $mapfile does not exist! "
       continue
    endif
    if("$arg" =~ "-top"*) then
        set top = $arg
        continue
    endif
    if("$arg" =~ *A) then
        set dist = "$arg"
        continue
    endif
    if("$arg" =~ dist=*) then
        set dist = `echo $arg | awk -F "=" '{print $NF+0 "A"}'`
        continue
    endif
    if("$arg" =~ *.pdb) then
        set pdbfile = $arg
        if(! -e "$pdbfile") echo "WARNING: $pdbfile does not exist! "
        continue
    endif
    if("$arg" =~ *.mtz) then
        set mtzfile = $arg
        if(! -e "$mtzfile") echo "WARNING: $mtzfile does not exist! "
        continue
    endif
    if("$arg" =~ *.map) then
        set mapfile = $arg
        if(! -e "$mapfile") echo "WARNING: $mapfile does not exist! "
        continue
    endif
end
if($?user_sigma) then
    echo "user-specified sigma: $user_sigma"
    set sigma = $user_sigma
else
    set user_sigma = 0
endif

if(! -e "$pdbfile" || ! ( -e "$mtzfile" || -e "$mapfile" ) ) then
    echo "usage: $0 refmacout.pdb refmacout.mtz"
    exit 9
endif
echo "pdb: $pdbfile"
if(-e "$mtzfile") echo "mtz: $mtzfile"
if(-e "$mapfile") echo "map: $mapfile"

set waterchain = `awk '/^ATOM|^HETAT/ && /HOH/{print substr($0,22,1)}' $pdbfile | tail -1`
if("$waterchain" == "") then
    # is there a non-water "S" chain?
    set test = `awk '/^ATOM|^HETAT/ && ! /HOH/{print substr($0,22,1)}' $pdbfile | grep "S" | wc -l`
    if( $test == 0 ) set waterchain = "S"
endif
if("$waterchain" == "") set waterchain = "_"
if($?user_chain) then
    set waterchain = "$user_chain"
endif
echo "waters appear to be in chain: $waterchain"
set watername = `awk '/^ATOM|^HETAT/ && /HOH/{print substr($0,14,3)}' $pdbfile | tail -1`
if("$watername" == "") set watername = "O"
echo "waters appear to be named: $watername"

set minocc = `awk '/^ATOM|^HETAT/{print substr($0,55,6)}' $pdbfile | sort -gr | awk '$1>0' | tail -1`
set B      = `awk '/^ATOM|^HETAT/{print substr($0,61,6)}' $pdbfile | sort -g | awk '{++n;b[n]=$1} END{print b[int(n/2)]}'`
set conf = `awk '/^ATOM|^HETAT/{c=substr($0,17,1);++count[c]} END{for(c in count) print count[c],c}' $pdbfile | sort -k1gr -k2 | awk '{print $2;exit}'`
echo "lowest occ is $minocc"
echo "median B is $B"
echo "dominant conformer: $conf"
if("$conf" == "") set conf = "_"

set occ = "auto"
if($?user_occ) then
    set occ = $user_occ
endif
if($?user_B) then
    set B = $user_B
endif
if($?user_conf) then
    set conf = $user_conf
endif

egrep -v "HOH" $pdbfile >! ${tempfile}protein.pdb
if("$dist" == "0A") then
    echo "" >! ${tempfile}protein.pdb
endif

# make sure new map matches grid of mask
set GRID = ""
set AXIS = ""
if(-e "$maskfile") then
    # extract map organization so we can match it
    echo xyzlim asu |\
    mapmask mapin $maskfile \
            mapout ${tempfile}mask.map > /dev/null
    echo | mapdump mapin ${tempfile}mask.map >&! ${tempfile}mapdump.txt
    set GRID   = `awk '/Grid sampling/{print "GRID",$(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump.txt`
    set AXIS   = `awk '/Fast, medium, slow/{print $(NF-2), $(NF-1), $NF}' ${tempfile}mapdump.txt | awk '! /[^ XYZ]/{print "AXIS",$0}'`
endif

# make a map from a provided mtz file (refmac or phenix)
if(-e "$mtzfile" && ! -e "$mapfile") then
    # make the fofc map
     fft hklin $mtzfile mapout ${tempfile}ffted.map << EOF >& /dev/null
    labin F1=FOFCWT PHI=PHFOFCWT
    $GRID
EOF
     fft hklin $mtzfile mapout ${tempfile}ffted.map << EOF >& /dev/null
    labin F1=DELFWT PHI=PHDELWT
    $GRID
EOF
    mapmask mapin ${tempfile}ffted.map mapout ${tempfile}fofc.map << EOF > /dev/null
    $AXIS
    xyzlim asu
EOF
#    set rmsrho = `echo | mapdump mapin ${tempfile}fofc.map | awk '/Rms deviation/{print $NF}'`
    if(-e "$maskfile") then
        echo "applying $maskfile to Fo-Fc map"
        echo maps mult |\
        mapmask mapin1 ${tempfile}fofc.map mapin2 ${tempfile}mask.map \
        mapout ${tempfile}masked.map >! ${tempfile}.log
        if(! -e ${tempfile}masked.map) then
            set BAD = "mask application failed"
            goto exit
        endif
        mv ${tempfile}masked.map ${tempfile}fofc.map
    endif
endif

if(-e "$mapfile" && "$mapfile" != "${tempfile}fofc.map") then
    cp -p $mapfile ${tempfile}fofc.map
endif
if("$rmsrho" == "") then
    set rmsrho = `echo | mapdump mapin ${tempfile}fofc.map | awk '/Rms deviation/{print $NF}'`
endif

set found = 0
set finding = 0
while ( $found == 0 && $finding < 1000 )
    pick.com ${sigma} ${tempfile}fofc.map $dist ${tempfile}protein.pdb $top -fast

    echo "$occ $B $watername $waterchain $user_sigma $conf $rmsrho"
    echo "$occ $B $watername $waterchain $user_sigma $conf $rmsrho" |\
    cat - pick.pdb |\
    awk 'NR==1{occ=$1;B=$2;c=65;atom=$3;chain=$4;sigma=$5;conf=$6;rmsrho=$7;\
         pi=4*atan2(1,1);\
         if(conf=="_")conf=" "\
         if(chain=="_")chain=" ";next}\
      /^ATOM|^HETAT/{++n;X=substr($0,31,8);Y=substr($0,39,8);Z=substr($0,47,8);height=substr($0,61,6)+0;\
        rho=height*rmsrho;ele=rho/(4*pi/B)**1.5;\
        if(occ=="auto"){occ=ele/8.0;if(occ>1)occ=1;if(occ<=0.01)occ=0.01};\
       if(height>=sigma) printf("ATOM%7d  %-3s%sHOH %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n"),\
      n,atom,conf,chain,n,X,Y,Z,occ,B,"O"}' |\
     tee new_water.pdb

    set found = `cat new_water.pdb | wc -l`
    if("$found" == "0" && "$top" == "") set finding = 9999
    if("$found" == "0" && "$top" != "") then
        set top = `echo $top | awk -F "=" '{print $1 "=" $2+1}'`
        @ finding = ( $finding + 1 )
    endif
end


egrep -v "^END|^TER" $pdbfile |\
awk '/^ATOM|^HETAT/ && substr($0,55,6)+0<=0 && /HOH/{\
   bad=substr($0,12,15)} bad==substr($0,12,15){next}\
   {print}' |\
cat >! new.pdb
awk '/^ATOM|^HETAT/{print "TAKEN",substr($0,23,4)}' new.pdb >! ${tempfile}taken.txt
cat ${tempfile}taken.txt new_water.pdb |\
awk 'BEGIN{n=0} /^TAKEN/{++taken[$2];next}\
   {++n;while(taken[n])++n;\
    printf("%s%4d%s\n",substr($0,1,22),n,substr($0,27))}' |\
tee ${tempfile}.pdb >> new.pdb
mv ${tempfile}.pdb new_water.pdb
echo "END" >> new.pdb

exit:
if($?BAD) then
    cat ${tempfile}.log
    echo "ERROR: $BAD"
    exit 9
endif

# clean up
rm -f ${tempfile}fofc.map >& /dev/null
rm -f ${tempfile}protein.pdb >& /dev/null
rm -f ${tempfile}ffted.map >& /dev/null
rm -f ${tempfile}masked.map >& /dev/null
rm -f ${tempfile}mask.map >& /dev/null
rm -f ${tempfile}mapdump.txt >& /dev/null
rm -f ${tempfile}.pdb >& /dev/null
rm -f ${tempfile}taken.txt >& /dev/null

exit

cat $pdbfile |\
egrep -v "^END|^TER" |\
awk '/^ATOM|^HETAT/ && substr($0,55,6)+0<=0 && /HOH/{\
   bad=substr($0,12,15)} bad==substr($0,12,15){next}\
   {print}' |\
cat >! new.pdb


