#! /bin/tcsh -f
#
#  label a PDB file with density values from map/mtz                     -James Holton  10-24-24
#
#
#
set pdbfile = ""
set mtzfile = "reference.mtz"
set mapfile = ""

set outfile = "rholabeled.pdb"

set tempfile = /dev/shm/${USER}/temp_rho_$$_
#set tempfile = ./tempfile_rho_
mkdir -p /dev/shm/${USER}
mkdir -p ${CCP4_SCR}



set mtzlabel = DELFWT
set phenixlabel = miller_array.labels.name
set test = `phenix.version | awk '/Release tag/{print ( $NF < 5000 )}'`
if( "$test" == "1" ) set phenixlabel = label



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
      if("$key" == "outpdb") set outfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
      if("$Arg" =~ *.mtz ) set mtzfile = $Arg
      if("$Arg" =~ *.map ) set mapfile = $Arg
    endif
    if("$key" == "debug") set debug = "$Val"
end

if(! -e "$pdbfile") then
    set BAD = "no coordinates provided"
    goto exit
endif

#if( $?debug && "$tempfile" =~ /dev/shm/* ) set tempfile = tempfile_Bud_

set t = "$tempfile"


if(! -e "$pdbfile" ) then
    set BAD = "no pdb file: $pdbfile "
    goto exit
endif
if(! -e "$mtzfile" && ! -e "$mapfile" ) then
    set BAD = "no mtz file: $mtzfile "
    goto exit
endif
if(! -e "$mtzfile" && ! -e "$mapfile" ) then
    set BAD = "no map file: $mapfile "
    goto exit
endif

#set test = `echo head | mtzdump hklin $mtzfile | egrep -i fpart | wc -l`
set mtzreso = `echo head | mtzdump hklin $mtzfile | awk '/Resolution Range/{getline;getline;print $6}'`



cat << EOF
pdbfile = $pdbfile
mtzfile = $mtzfile
outfile = $outfile

tempfile = $tempfile
debug    = $?debug
EOF

set t = $tempfile


awk '/^ATOM|^HETAT/' $pdbfile >! ${t}.pdb
awk '{x=substr($0,31,8);y=substr($0,39,8);z=substr($0,47,8);\
  printf("(%.3f,%.3f,%.3f) %d KEY\n",x,y,z,++n)}' ${t}.pdb >! ${t}key.txt
phenix.map_value_at_point $mtzfile ${t}.pdb \
          ${phenixlabel}=$mtzlabel scale=sigma |\
cat ${t}key.txt - |\
awk '$NF=="KEY"{n[$1]=$2;next}\
  /Map value:/{print "RHO",n[$3],$3,++i,$NF}' |\
cat - $pdbfile |\
awk '/^RHO/{rho[$2]=$NF;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {++n;print $0,"           ",rho[n]}' |\
cat >! $outfile





exit:

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?debug && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    echo "clearing temp files"
    rm -f ${t}*
endif


exit

