#! /bin/csh -f
#
#	Correlate two maps, regardless of grid sampling
#
#
#
# set this to wherever your awk program is
alias nawk /usr/bin/nawk
nawk 'BEGIN{print}' >& /dev/null
if($status) alias nawk awk


set map1 = ""
set map2 = ""

if("$1" != "") set map1 = "$1"
if("$2" != "") set map2 = "$2"

set tempfile = ${CCP4_SCR}/corr_temp$$
set logfile  = correlate.log
set logfile = /dev/null

if((! -e "$map1")||(! -e "$map2")) then
    cat << EOF

usage: $0 map1.map map2.map

where map1.map and map2.map are the two maps you want to compare

EOF
    exit 2
endif

#########################################################
# check out the maps
echo "go" | mapdump mapin $map1 >! ${tempfile}mapdump1
set SG1     = `nawk '/Space-group/{print $NF; exit}' ${tempfile}mapdump1`
set CELL1   = `nawk '/Cell dimensions/{print $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump1`
set GRID1   = `nawk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump1`
set LIMITS1 = `nawk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${tempfile}mapdump1`
set AXIS1   = `nawk '/Fast, medium, slow/{print $(NF-2), $(NF-1), $NF}' ${tempfile}mapdump1 | nawk '! /[^ XYZ]/'`

echo "go" | mapdump mapin $map2 >! ${tempfile}mapdump2
set SG2     = `nawk '/Space-group/{print $NF; exit}' ${tempfile}mapdump2`
set CELL2   = `nawk '/Cell dimensions/{print $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump2`
set GRID2   = `nawk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump2`
set LIMITS2 = `nawk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${tempfile}mapdump2`
set AXIS2   = `nawk '/Fast, medium, slow/{print $(NF-2), $(NF-1), $NF}' ${tempfile}mapdump2 | nawk '! /[^ XYZ]/'`

rm -f ${tempfile}mapdump1 ${tempfile}mapdump2 >& /dev/null

if("$SG1" != $SG2) then
    echo "WARNING: maps have different space groups! "
    echo "         they are NOT going to correlate well."
endif

# re-sample the second map if grids/limits don't match
if("$GRID1" != "$GRID2") then    
    goto resample
endif

# see if only re-organization is necessary
if(("$AXIS1" != "$AXIS2")||("$LIMITS1" != "$LIMITS2")) then
    cp $map2 ${tempfile}.map
    goto reorganize
endif



again:
# actually calculate the correlation coefficient
echo "correlate section" |\
overlapmap mapin1 $map1 mapin2 $map2 | tee -a $logfile |\
grep "Total corr"

if($?KEEP_RESAMPLED && $?RESAMPLED) then
    echo "keeping resampled.map"
    cp $map2 resampled.map
endif
rm -f ${tempfile}extended.map >& /dev/null
rm -f ${tempfile}map2.map >& /dev/null



exit

resample:
set RESAMPLED
echo "re-sampling $map2"
echo "from $GRID2 to $GRID1"

mapmask mapin $map2 mapout ${tempfile}extended.map << EOF >> $logfile
xyzlim 0 1 0 1 0 1
axis $AXIS2
EOF

set hiRES1 = `echo "$GRID1 $CELL1" | nawk '{print 3*$4/$1; print 3*$5/$3; print 3*$6/$3}' | sort -n | tail -1`
set hiRES2 = `echo "$GRID2 $CELL2" | nawk '{print 3*$4/$1; print 3*$5/$3; print 3*$6/$3}' | sort -n | tail -1`

# back-transform the second map to HKLs
sfall mapin ${tempfile}extended.map \
      hklout ${tempfile}map2.mtz << EOF | tee ${tempfile}.log >> $logfile
symm $SG2
mode sfcalc mapin
RESOLUTION $hiRES2
EOF
if($status) then
    # check axis order
    set temp = `cat ${tempfile}.log | nawk '/Check Iuvw/{printf "%c %c %c", 87+$6, 87+$7, 87+$8}' | nawk '! /[^ XYZ]/'`
    if(("$temp" != "$AXIS2")&&($#temp == 3)) then
	set AXIS2 = ( $temp )
	goto resample
    endif
    echo "ERROR: unsolvable problem with SFALL."
    echo "       perhaps it would be easier to just"
    echo "       re-calculate $map2 from wherever it came from"
    echo "       using: GRID $GRID1"
    echo "       in fft"
    exit 2
endif
rm -f ${tempfile}.log >& /dev/null

# no re-make the 2nd map on the 1st map's grid
fft hklin ${tempfile}map2.mtz mapout ${tempfile}.map << EOF >> $logfile
RESOLUTION 1000 $hiRES1
GRID $GRID1
LABIN F1=FC PHI=PHIC
EOF
rm -f ${tempfile}map2.mtz >& /dev/null

reorganize:
# and extend it to the same limits
mapmask mapin ${tempfile}.map mapout ${tempfile}map2.map << EOF >> $logfile
xyzlim $LIMITS1
axis   $AXIS1
pad    0
EOF
rm -f ${tempfile}.map >& /dev/null

set map2 = ${tempfile}map2.map  

goto again

