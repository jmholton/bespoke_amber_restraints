#! /bin/tcsh -f
#
#	report "rho" values of atoms in a map     - James Holton 9-3-25
#
#  usage: check_map_sites.com atoms.pdb density.map|coefficients.mtz [mtzlabel=auto] [scale=auto|sigma|volume]
#
#  prints the spline-interpolated density at each atom, one value per line in the order of the atoms.
#  Samplers, the first that gives a value for every atom wins:
#    mapman            "peek value ... spline" (fast, but not installed everywhere)
#    Bfac_map peek     the same spline in C (Bfac_map.c next to this script, compiled on demand)
#    phenix.map_value_at_point   linear, last resort (its output format and order change between versions)
#  An mtz is first made into a map with fft, grid_sample points per resolution element (mtzlabel=auto: DELFWT or
#  FOFCWT, else FWT or 2FOFCWT, else the first F; the phase is the next P column).  scale=sigma divides the map by
#  its rms over the cell (as phenix.map_value_at_point scale=sigma); scale=auto is sigma for an mtz, as-is for a map.
#

set MAPMAN =  /programs/rave/lx_mapman
if(! -e "$MAPMAN") set MAPMAN = mapman

set pdbfile = ""
set mapfile = fofc.map
set mtzfile = ""
set mtzlabel = auto
set scale = auto
# spline at d/4 reads fofc at atoms with slope 0.987 against a d/12 tricubic reference (d/3: 0.966, d/6: 0.997)
set grid_sample = 4
set pdir = `dirname $0`

set tempfile = /dev/shm/tempfile$$_

foreach arg ( $* )
    if("$arg" =~ *.pdb) set pdbfile = "$arg"
    if("$arg" =~ *.map) set mapfile = "$arg"
    if("$arg" =~ *.mtz) set mtzfile = "$arg"
    if("$arg" =~ mtzlabel=*) set mtzlabel = `echo "$arg" | awk -F "=" '{print $2}'`
    if("$arg" =~ scale=*) set scale = `echo "$arg" | awk -F "=" '{print $2}'`
    if("$arg" =~ grid_sample=*) set grid_sample = `echo "$arg" | awk -F "=" '{print $2}'`
    # mapman=none skips it (to test the other samplers)
    if("$arg" =~ mapman=*) set MAPMAN = `echo "$arg" | awk -F "=" '{print $2}'`
    if("$arg" =~ tempfile=*) then
        set tempfile = `echo "$arg" | awk -F "=" '{print $2}'`
        setenv DEBUG
    endif
end

if(! -e "$pdbfile" || ( ! -e "$mapfile" && ! -e "$mtzfile" ) ) then
    cat << EOF
usage: $0 atoms.pdb density.map|coefficients.mtz [mtzlabel=auto] [scale=auto|sigma|volume]

the spline-interpolated electron density at each atom position will be
printed out in the order of the atoms in the PDB file

EOF
    exit 9
endif

cat ${pdbfile} |\
awk '/^ATOM|^HETAT/{print substr($0,1,55),"1.00 15.00"}' |\
cat >! ${tempfile}_probe.pdb

set natoms = `cat ${tempfile}_probe.pdb | wc -l`
if( $natoms == 0 ) then
   set BAD = "no atoms in $pdbfile"
   goto exit
endif

# map coefficients: fft them to a map
if( -e "$mtzfile" ) then
    echo head | mtzdump hklin $mtzfile |\
    awk '/Column Labels :/{getline;getline;for(i=1;i<=NF;++i)label[i]=$i}\
      /Column Types :/{getline;getline;for(i=1;i<=NF;++i)print $i,label[i]}' >! ${tempfile}labels.txt
    if( "$mtzlabel" == "auto" ) then
        set mtzlabel = `awk '$2=="DELFWT"||$2=="FOFCWT"{print $2;exit}' ${tempfile}labels.txt`
        if( "$mtzlabel" == "" ) set mtzlabel = `awk '$2=="FWT"||$2=="2FOFCWT"{print $2;exit}' ${tempfile}labels.txt`
        if( "$mtzlabel" == "" ) set mtzlabel = `awk '$1=="F"{print $2;exit}' ${tempfile}labels.txt`
    endif
    set phlabel = `awk -v f="$mtzlabel" '$2==f{on=1;next} on && $1=="P"{print $2;exit}' ${tempfile}labels.txt`
    if( "$mtzlabel" == "" || "$phlabel" == "" ) then
        set BAD = "no amplitude ($mtzlabel) and phase ($phlabel) in $mtzfile"
        goto exit
    endif
    fft hklin $mtzfile mapout ${tempfile}fft.map << EOF >&! ${tempfile}fft.log
labin F1=$mtzlabel PHI=$phlabel
grid sample $grid_sample
EOF
    if( $status || ! -e ${tempfile}fft.map ) then
        set BAD = "fft of $mtzlabel $phlabel in $mtzfile failed"
        goto exit
    endif
    set mapfile = ${tempfile}fft.map
    if( "$scale" == "auto" ) set scale = sigma
endif
if( "$scale" == "sigma" ) then
    mapmask mapin $mapfile mapout ${tempfile}sigma.map << EOF >&! ${tempfile}sigma.log
xyzlim cell
scale sigma 1 0
EOF
    if( $status || ! -e ${tempfile}sigma.map ) then
        set BAD = "sigma-scaling $mapfile failed"
        goto exit
    endif
    set mapfile = ${tempfile}sigma.map
endif
touch ${tempfile}peeks ${tempfile}mapman.log ${tempfile}phenix.log

# 1: mapman
if( "$MAPMAN" != "none" ) then
  echo border 3  |\
    mapmask xyzin ${tempfile}_probe.pdb \
    mapin $mapfile mapout ${tempfile}_probeme.map >&! ${tempfile}mapmask.log
  if( $status ) then
    mapmask mapin $mapfile mapout ${tempfile}_P1.map << EOF >! ${tempfile}mapmask.log
    xyzlim cell
EOF
    if( $status ) then
      mapmask mapin $mapfile mapout ${tempfile}_P1.map << EOF >! ${tempfile}mapmask.log
      xyzlim cell
      pad 0
EOF
    endif

    mapmask xyzin ${tempfile}_probe.pdb \
      mapin ${tempfile}_P1.map mapout ${tempfile}_probeme.map << EOF >> ${tempfile}mapmask.log
    extend xtal
    symm 1
    border 3
EOF
  endif

  setenv MAPSIZE `ls -l ${tempfile}_probeme.map | awk '{printf "%d", $5/3.5}'`
  $MAPMAN -b mapsize $MAPSIZE << mapman-end >&! ${tempfile}mapman.log
read map1 ${tempfile}_probeme.map ccp4
peek value map1 ${tempfile}_probe.pdb ${tempfile}_probed.pdb spline ;
quit
mapman-end

  # maybe this version of mapman doesnt peek-a-boo
  awk '/PEEK-A-BOO/{print $NF+0}' ${tempfile}mapman.log >! ${tempfile}peeks
endif

# 2: the same spline in C
set test = `cat ${tempfile}peeks | wc -l`
if( "$test" != "$natoms" ) then
    set Bmapexe = ${pdir}/Bfac_map
    set stale = 1
    if( -x "$Bmapexe" ) set stale = `find ${pdir}/Bfac_map.c -newer $Bmapexe | wc -l`
    if( $stale && -e ${pdir}/Bfac_map.c ) then
        gcc -O3 -o ${Bmapexe}_$$ ${pdir}/Bfac_map.c -lm >&! ${tempfile}gcc.log && mv -f ${Bmapexe}_$$ $Bmapexe
    endif
    if( -x "$Bmapexe" ) then
        echo xyzlim cell |\
        mapmask mapin ${mapfile} mapout ${tempfile}cell.map >&! ${tempfile}cell.log
        $Bmapexe peek pdb=${tempfile}_probe.pdb map=${tempfile}cell.map interp=cubic >! ${tempfile}peeks
    endif
endif

# 3: phenix
set test = `cat ${tempfile}peeks | wc -l`
if( "$test" != "$natoms" ) then
    set phenixlabel = miller_array.labels.name
    set test = `phenix.version | awk '/Release tag/{print ( $NF < 5000 )}'`
    if( "$test" == "1" ) set phenixlabel = label

    echo xyzlim cell |\
    mapmask mapin ${mapfile} mapout ${tempfile}cell.map >&! ${tempfile}phenix.log
    echo symm 1 |\
    mapmask mapin ${tempfile}cell.map mapout ${tempfile}P1.map >&! ${tempfile}phenix.log
    phenix.map_to_structure_factors ${tempfile}P1.map \
       d_min=0.5 b_blur=0 scale_max=None \
       output_file_name=${tempfile}test.mtz >>& ${tempfile}phenix.log
    cad hklin1 ${tempfile}test.mtz hklout ${tempfile}onecol.mtz << EOF >>& ${tempfile}cad.log
    labin file 1 E1=F E2=PHIF
EOF
    echo head | mtzdump hklin ${tempfile}onecol.mtz |\
    awk '/Cell Dimensions :/{getline;getline;\
      a=$1;b=$2;c=$3;al=$4;be=$5;ga=$6;next}\
      /Space group =/{split($0,s,"\047");sg=s[2];next}\
      END{printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s\n",a,b,c,al,be,ga,sg}' |\
    cat - ${tempfile}_probe.pdb >! ${tempfile}_probe_c1.pdb
    # on the input map's own grid spacing, and back in atom order (phenix reorders atoms)
    set gridstep = `echo | mapdump mapin ${mapfile} | awk '/Grid sampling/{nx=$(NF-2);ny=$(NF-1);nz=$NF} /Cell dimensions/{a=$(NF-5);b=$(NF-4);c=$(NF-3)} END{s=a/nx;if(b/ny<s)s=b/ny;if(c/nz<s)s=c/nz;printf("%.3f",s)}'`
    phenix.map_value_at_point ${tempfile}onecol.mtz ${tempfile}_probe_c1.pdb \
      ${phenixlabel}=F scale=volume grid_step=$gridstep |\
    tee -a ${tempfile}phenix.log |\
    cat ${tempfile}_probe_c1.pdb - | map_values_inorder.awk |\
    awk '/Map value:/{print $NF}' >! ${tempfile}peeks
endif

set test = `cat ${tempfile}peeks | wc -l`
if("$test" != "$natoms") then
    echo "WARNING: $test peeks for $natoms atoms! "
    cat ${tempfile}mapman.log
    cat ${tempfile}phenix.log
    stat ${tempfile}_probe.pdb
    set BAD = "neither $MAPMAN, Bfac_map peek nor phenix.map_value_at_point are working."
    goto exit
endif
cat ${tempfile}peeks

if(! $?DEBUG) rm -f ${tempfile}*

exit:

if( $?BAD ) then
    echo "ERROR: $BAD"
    exit 9
endif


exit
