#! /bin/tcsh -f
#
#	report "rho" values of atoms in a map     - James Holton 9-3-25
#
#


set MAPMAN =  /programs/rave/lx_mapman
if(! -e "$MAPMAN") set MAPMAN = mapman

set pdbfile = ""
set mapfile = fofc.map

set tempfile = /dev/shm/tempfile$$_

foreach arg ( $* )
    if("$arg" =~ *.pdb) set pdbfile = "$arg"
    if("$arg" =~ *.map) set mapfile = "$arg"
    if("$arg" =~ tempfile=*) then
        set tempfile = `echo "$arg" | awk -F "=" '{print $2}'`
        setenv DEBUG
    endif
end

if(! -e "$pdbfile" || ! -e "$mapfile") then
    cat << EOF
usage: $0 atoms.pdb density.map

the spline-interpolated electron density at each atom position will be 
printed out in the order of the atoms in the PDB file

EOF
    exit 9
endif

cat ${pdbfile} |\
awk '/^ATOM|^HETAT/{print substr($0,1,55),"1.00 15.00"}' |\
cat >! ${tempfile}_probe.pdb

set test = `cat ${tempfile}_probe.pdb | wc -l`
if( $test == 0 ) then
   set BAD = "no atoms in $pdbfile"
   goto exit
endif

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


awk '/PEEK-A-BOO/{print $NF+0}' ${tempfile}mapman.log | tee ${tempfile}peeks

# maybe this version of mapman doesnt peek-a-boo 
set test = `cat ${tempfile}peeks | wc -l`
if("$test" == "0") then
    # mapman gave no peeks (missing, wrong version, or crashed) - fail over to phenix
    if( 1 ) then
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
        phenix.map_value_at_point ${tempfile}onecol.mtz ${tempfile}_probe_c1.pdb \
          ${phenixlabel}=F scale=volume |\
        tee -a ${tempfile}phenix.log |\
        awk '/Map value:/{print $NF}' | tee ${tempfile}phenixpeeks
    endif
    set test = `cat ${tempfile}phenixpeeks | wc -l`
    if( ! $test ) then
        set BAD = "neither $MAPMAN nor phenix.map_value_at_point are working."
        goto exit
    endif
    # phenix produced the peeks - hand them back in the mapman-format output file
    cp ${tempfile}phenixpeeks ${tempfile}peeks
endif
set test = `cat ${tempfile}peeks | wc -l`
if("$test" == "0") then
    echo "WARNING: no peeks! "
    cat ${tempfile}mapman.log
    cat ${tempfile}phenix.log
    stat ${tempfile}_probeme.map
    stat ${tempfile}_probe.pdb
    awk '/^ATOM|^HETAT/{print substr($0,61)}' ${tempfile}_probed.pdb
endif

if(! $?DEBUG) rm -f ${tempfile}*

exit:

if( $?BAD ) then
    echo "ERROR: $BAD"
    exit 9
endif


exit



