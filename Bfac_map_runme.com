#! /bin/tcsh -f
#
#  build a B-factor field B(x,y,z) as a CCP4 map from a PDB with B factors
#                                                                -James Holton 8-20-26
#
#  The map replaces the per-atom B sidecar (Bfac.pdb) at structure-factor time:
#  each atom takes its B by sampling the map at its position, so B belongs to the
#  location rather than to the atom.  A water that lands on an ordered site is
#  immediately sharp; one that wanders into a channel is immediately diffuse.
#
#  usage:  Bfac_map_runme.com [Bfac.pdb] [outfile=Bfac.map] [sigma=0.5] ...
#
#  The heavy lifting is in Bfac_map (C).  Run "Bfac_map -h" for every option.
#
set pdbfile   = Bfac.pdb
set outfile   = Bfac.map

# Gaussian donation width, A.  The one parameter that matters: how far B blends
# between neighbouring atoms of different mobility.  0.5 reproduces the per-atom
# structure factors to ~1% on F; 1.0 costs ~3%; 1.5 is too much (15%).
set sigma     = 0.5
# voxel size, A.  Hardly matters next to sigma, but sets the file size.
set grid      = 0.5
# B assigned to space that has no atom within fardist
set farB      = 999
set fardist   = 4.5

# conditioning applied to the atomic B factors before they are deposited, so the
# map holds final values.  Same meaning as in nc2mtz.  maxB=0 disables clipping;
# farB is never clipped by it.
set minB      = 0
set maxB      = 0
set Bscale    = 1
set Boffset   = 0
# weight each atom's donation by its occupancy
set useocc    = 0

# Space group to expand the source model by, so that a model which is only one
# asymmetric unit still fills the cell.  Deliberately NOT defaulted from
# $smallSG in xtal_properties.sourceme: the usual source here is Bfac.pdb, which
# matches the amber system and is therefore already expanded - and an MD
# supercell is not symmetric, so expanding it would average symmetry mates
# together and erase the very differences between subcells that it exists to
# capture.  Set this only when the input really is a single ASU.
set sg        = P1
# or give the operators directly, ';' separated, e.g. "X,Y,Z;-X,-Y,Z"
set symops    = ""
# where to look sg= up (default $CLIBD/symop.lib)
set symoplib  = ""

# the coordinate frame of the map.  Leave cell empty to take CRYST1 from the PDB
# multiplied by super_mult.  A supercell PDB carrying a primitive CRYST1 (what
# nc2mtz writes) needs super_mult to get this right.
set cell        = ""
set super_mult  = ( 1 1 1 )
# exact grid, overriding grid= (use to match an existing map)
set nxyz      = ""

# probe the source PDB afterwards and report how faithfully the field reproduces
# the B factors that went into it
set check     = 1

set debug = 0

set pdir = `dirname $0`
set path = ( $pdir $path )

set tempfile = /dev/shm/${USER}/temp_Bmap_$$_

foreach sourceme ( xtal_properties.sourceme user_settings.sourceme )
   if(-e $sourceme ) then
      echo "sourcing $sourceme"
      source $sourceme
   endif
end

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
          set $Key = "$Val"
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "output" || "$key" == "outmap") set outfile = "$Val"
      if("$key" == "pdb") set pdbfile = "$Val"
      if("$key" == "md_mult") set super_mult = ( $Val )
      if("$key" == "bsigma") set sigma = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = "$Arg"
      if("$Arg" =~ *.map ) set outfile = "$Arg"
    endif
    if("$arg" == "debug") set debug = "1"
end

mkdir -p `dirname $tempfile` > /dev/null
if( $status ) set tempfile = ./temp_Bmap_$$_
set t = "$tempfile"

set super_mult = `echo $super_mult | awk '{gsub("[,x]"," ");print}'`
if( "$super_mult" == "" ) set super_mult = ( 1 1 1 )
while ( $#super_mult != 3 )
  set super_mult = ( $super_mult $super_mult[$#super_mult] )
end
set super_mult_csv = `echo $super_mult | awk '{print $1","$2","$3}'`

# the C program does the work
set exe = "$pdir/Bfac_map"
if( ! -x "$exe" ) then
  if( -e "$pdir/Bfac_map.c" ) then
    echo "compiling Bfac_map ..."
    gcc -O3 -o "$exe" "$pdir/Bfac_map.c" -lm
    if( $status || ! -x "$exe" ) then
      set BAD = "cannot compile $pdir/Bfac_map.c"
      goto exit
    endif
  else
    set BAD = "cannot find Bfac_map next to $0 - compile it with: gcc -O3 -o Bfac_map Bfac_map.c -lm"
    goto exit
  endif
endif

if( ! -r "$pdbfile" ) then
  set BAD = "cannot read $pdbfile"
  goto exit
endif

set cellopt = ""
if( "$cell" != "" ) then
  set cellopt = "cell=`echo $cell | awk '{gsub(" +",",");print}'`"
else
  set cellopt = "super_mult=$super_mult_csv"
endif
set gridopt = "grid=$grid"
if( "$nxyz" != "" ) set gridopt = "nxyz=`echo $nxyz | awk '{gsub("[ ,x]+",",");print}'`"

set symopt = ""
if( "$sg" != "" ) set symopt = "sg=$sg"
if( "$symops" != "" ) set symopt = "symops=$symops"
if( "$symoplib" != "" ) set symopt = "$symopt symoplib=$symoplib"

cat << EOF
pdbfile   $pdbfile
outfile   $outfile
sigma     $sigma A       (B blending width)
$gridopt
farB      $farB          beyond fardist $fardist A from any atom
B conditioning: Bscale=$Bscale Boffset=$Boffset minB=$minB maxB=$maxB
frame     $cellopt
symmetry  $symopt
EOF

$exe build pdb="$pdbfile" outmap="$outfile" \
     sigma=$sigma $gridopt farB=$farB fardist=$fardist \
     minB=$minB maxB=$maxB Bscale=$Bscale Boffset=$Boffset useocc=$useocc \
     $cellopt $symopt
if( $status || ! -e "$outfile" ) then
  set BAD = "Bfac_map build failed"
  goto exit
endif

if( $check ) then
  echo ""
  echo "checking the field against the B factors that went into it..."
  $exe probe pdb="$pdbfile" map="$outfile" outpdb=${t}probed.pdb $cellopt
  if( $status ) then
    set BAD = "Bfac_map probe failed"
    goto exit
  endif
  # same conditioning the build applied, so like is compared with like
  echo "$Bscale $Boffset $minB $maxB" |\
  cat - "$pdbfile" |\
  awk 'NR==1{Bscale=$1;Boffset=$2;minB=$3;maxB=$4;next}\
    ! /^ATOM|^HETAT/{next}\
    {B=substr($0,61,6)*Bscale+Boffset}\
    B<minB{B=minB} maxB>0 && B>maxB{B=maxB}\
    {printf("%.4f IN\n",B)}' |\
  cat - ${t}probed.pdb |\
  awk '$NF=="IN"{++i;in_[i]=$1;next}\
    ! /^ATOM|^HETAT/{next}\
    {++n;d=substr($0,61,6)-in_[n];ss+=d*d;\
     if(d*d>worst*worst){worst=d;wn=n}}\
    END{if(n<1){print "no atoms to check";exit}\
      printf("%d atoms: rms dB = %.3f   worst dB = %+.2f (atom %d)\n",\
             n,sqrt(ss/n),worst,wn)}'
endif

echo ""
$exe stats map="$outfile"

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9
endif

if( $debug ) exit
rm -f ${t}* > /dev/null

exit
