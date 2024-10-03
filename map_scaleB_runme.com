#! /bin/tcsh -f
#
#   apply a scale and B factor to a map file
#   preserve offset
#
#

set mapfile = "solvent.map"
set outputmap  = "scaled.map"

set reso = 0.95
set scale = 1
set B     = 0
set preclip  = ""
set clip  = ""

set sfallAXIS = ( Z X Y )

set pdir = `dirname $0`
set path = ( $path $pdir )


set tempfile = /dev/shm/${USER}/solmap_$$_
mkdir -p /dev/shm/${USER}

set logfile = maprescale_details.log

foreach arg ( $* )
    set Key = `echo $arg | awk -F "=" '{print $1}'`
    set Val = `echo $arg | awk -F "=" '{print $2}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if("$Key" =~ *.map && "$val" == "") set mapfile = $Key
    if("$key" == "mapfile") set mapfile = "$Val"

    if("$key" == "output") set outputmap = "$Val"
    if("$key" == "outfile") set outputmap = "$Val"
    if("$key" == "reso") set reso = "$num"
    if("$key" == "scale") set scale = "$num"
    if("$Key" == "B") set B = "$num"
    if("$key" == "clip") set clip = "$num"
    if("$key" == "preclip") set preclip = "$num"

    if("$key" == "tempfile") then
        set tempfile = "$Val"
        setenv DEBUG
    endif
    if("$key" == "debug") then
        setenv DEBUG
    endif
end

if( "$preclip" != "" ) then
  echo "clipping at rho=$preclip "
  map_func.com -func max -param $preclip $mapfile \
  -output ${tempfile}positive.map >> $logfile
  mv ${tempfile}positive.map ${tempfile}new.map
else
  cp $mapfile ${tempfile}new.map
endif

# apply any scale factor first
echo "scaling $mapfile by $scale"
mapmask mapin ${tempfile}new.map mapout ${tempfile}scaled.map << EOF >> $logfile
scale factor $scale
EOF

if("$B" == "0") then
    cp ${tempfile}scaled.map $outputmap
    goto exit
endif

echo "applying B= $B"
echo | mapdump mapin $mapfile >! ${tempfile}mapdump.txt
set GRID = `awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump.txt`
set AXES = `awk '/Fast, medium, slow/{print $(NF-2), $(NF-1), $NF}' ${tempfile}mapdump.txt | awk '! /[^ XYZ]/'`
set CELL = `awk '/Cell dimensions /{print $4,$5,$6,$7,$8,$9;exit}' ${tempfile}mapdump.txt`
set SGnum  = `awk '/ Space-group /{print $NF;exit}' ${tempfile}mapdump.txt`
set LIMITS = `awk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${tempfile}mapdump.txt`

echo | mapdump mapin ${tempfile}scaled.map >! ${tempfile}mapdump.txt
set avgrho = `awk '/Mean density/{print $NF}' ${tempfile}mapdump.txt`

sfallagain:
mapmask mapin ${tempfile}scaled.map mapout ${tempfile}sfallme.map << EOF >> $logfile
xyzlim cell
axis $sfallAXIS
EOF

echo "sfall"
sfall mapin ${tempfile}sfallme.map hklout ${tempfile}sfalled.mtz << EOF | tee ${tempfile}sfall.log >> $logfile
mode sfcalc mapin
SFSG 1
reso $reso
EOF
if($status) then
  # check axis order
  set temp = `cat ${tempfile}sfall.log | awk '/Check Iuvw/{printf "%c %c %c", 87+$6, 87+$7, 87+$8}' | awk '! /[^ XYZ]/'`
  if(("$temp" != "$sfallAXIS")&&($#temp == 3)) then
	  set sfallAXIS = ( $temp )
	  goto sfallagain
  endif
  set BAD = "problem with SFALL"
  goto exit
endif


echo "fft"
fft hklin ${tempfile}sfalled.mtz mapout ${tempfile}ffted.map << EOF >> $logfile
labin F1=FC PHI=PHIC
scale F1 1 $B
GRID $GRID
EOF

echo "adding back offset: $avgrho"
mapmask mapin ${tempfile}ffted.map mapout ${tempfile}new.map << EOF >> $logfile
scale factor 1 $avgrho
EOF

if( "$clip" != "" ) then
  echo "clipping at rho=$clip "
  map_func.com -func max -param $clip ${tempfile}new.map \
  -output ${tempfile}positive.map >> $logfile
  mv ${tempfile}positive.map ${tempfile}new.map
endif


# make sure output header is same as input
echo axis $AXES |\
mapmask mapin ${tempfile}new.map mapout ${tempfile}axis.map >> $logfile
mapmask mapin ${tempfile}axis.map maplim $mapfile mapout ${tempfile}maplim.map << EOF >> $logfile
xyzlim match
EOF

cp ${tempfile}maplim.map $outputmap



exit:

ls -l $outputmap

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?DEBUG && "$tempfile" != "./") then
    echo "cleaning up..."
    rm -rf ${tempfile}* >& /dev/null
endif

exit


mapmask mapin solvent.map mapout sfallme.map << EOF
xyzlim cell
axis Z X Y
EOF

sfall mapin sfallme.map hklout sfalled.mtz << EOF
mode sfcalc mapin
resolution $reso
sfsg 1
EOF


echo reindex h/2,k/2,l/3 |\
 reindex hklin erefmac_Prod/refme.mtz hklout temp.mtz
cad hklin1 temp.mtz hklout this.mtz << EOF
labin file 1 E1=Fpart E2=PHIpart
labou file 1 E1=Fthis E2=PHIthis
symm P212121
EOF
cad hklin1 refmac_Prod/refme.mtz hklout that.mtz << EOF
labin file 1 E1=Fpart E2=PHIpart
labou file 1 E1=Fthat E2=PHIthat
symm P212121
EOF
cad hklin1 this.mtz hklin2 sfalled.mtz hklout phases.mtz << EOF
labin file 1 E1=PHIthis
labin file 2 E1=PHIC
EOF
mtz2txt phases.mtz 

awk '$4~/[0-9]/ && $5~/[0-9]/{print sqrt(($4-$5)^2),$0}' phases.csh  | sort -g


