#! /bin/tcsh -f
#
#   look for waters with bad non-bonds  these will probably explode amber
#   also check for vacuum-filled voids
#
set pdbfile = ""

set outfile = voids.map

set water_radius  = 2.4

set tempfile = /dev/shm/${USER}/temp_void_$$_
mkdir -p /dev/shm/${USER}
setenv CCP4_SCR ${tempfile}_dir/
mkdir -p $CCP4_SCR

# scan command line for settings
foreach Arg ( $* )
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk -F "=" '{print $2}'`
    set assign = `echo $Arg | awk '{print ( /=/ )}'`
    set arg = `echo $Arg | awk '{print tolower($0)}'`

    if( $assign ) then
        set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
        if ( $test ) then
           set $Key = $Val
           echo "$Key = $Val"
        endif
    else
        if( "$Arg" =~ *.pdb ) then
          if( "$pdbfile" == "" ) then
            set pdbfile = $Arg
          endif
        endif
        if( "$arg" == "debug" ) set debug
    endif
end

if( $?DEBUG ) set debug
if( $?debug ) set tempfile = tempfile


cat << EOF
pdbfile = $pdbfile

water_radius = $water_radius

tempfile       = $tempfile
EOF


set t = ${tempfile}

convert_pdb.awk -v skip=H $pdbfile |\
 awk '/^CRYST/ && ! c{print;++c;next}\
 ! /^ATOM|^HETAT/{next} \
 {print "ATOM  ",substr($0,8,73)}' >! ${t}_noH.pdb


echo "looking for voids"
sfall xyzin ${t}_noH.pdb mapout ${t}voidme.map << EOF >! ${t}sfall.log
MODE ATMMAP SOLVMAP
SYMM 1
grid 128 128 128
vdwrad $water_radius
EOF
if($status) then
   mv ${t}_noH.pdb sfallme.pdb
   mv ${t}sfall.log sfall.log
   goto geom
endif
float_func -func segment -header 1104 -xsize 128 -ysize 128 \
  ${t}voidme.map -outfile ${t}segments.map >! ${t}segment.log

float_add -histogram -header 1104 ${t}segments.map -outfile ${t}.bin |\
 tee ${t}segment.log |\
 awk 'NF==2 && $1+0!=0 && $2+0!=0{print $2}' |\
 sort -gr | head -n 1 >! ${t}biggest_void.txt

set void = `awk '{print $1;exit}' ${t}biggest_void.txt`
if("$void" == "") set void = 0
echo "biggest void: $void voxels"


map_scaleB_runme.com ${t}voidme.map B=50 outfile=${t}scaled.map  >> /dev/null
echo xyzlim -5 133 -5 133 -5 133 |\
 mapmask mapin ${t}scaled.map mapout ${outfile} >> /dev/null


exit:
if( $?BAD ) then
  echo "ERROR: $BAD"
  exit 9
endif

if(! $?debug) then
  echo "clearng temp files"
  rm -rf ${t}*
endif

exit


