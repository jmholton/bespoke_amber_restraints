#! /bin/tcsh -f
#
#  label a PDB file with density values from map/mtz                     -James Holton 7-7-26
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



set mtzlabel = auto
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

cat << EOF
pdbfile = $pdbfile
mtzfile = $mtzfile
outfile = $outfile

tempfile = $tempfile
debug    = $?debug
EOF

set test = `tail -n 50 $pdbfile | awk '/^ATOM|^HETAT/ && substr($0,12,1)!=" "{found=1} END{print found+0}'`
if( $test ) then
    awk '/^ATOM|^HETAT/{print;exit}' $pdbfile
    tail -n 50 $pdbfile | tac |awk '/^ATOM|^HETAT/{print;exit}' 
    set BAD = "register-shifted PDB (6-digit atom numbers) please check what generated it"
    goto exit
endif

set t = $tempfile

echo head | mtzdump hklin $mtzfile |\
tee ${t}mtzdump.txt |\
awk '/Column Labels :/{getline;getline;for(i=1;i<=NF;++i)label[i]=$i}\
  /Column Types :/{getline;getline;for(i=1;i<=NF;++i)print $i,label[i]}' |\
cat >! ${t}labels.txt


if( "$mtzlabel" == "auto" ) then
  set firstlabel = `awk '$2=="DELFWT" || $2=="FOFCWT"{print $2}' ${t}labels.txt | head -n 1`
  if( "$firstlabel" != "" ) set mtzlabel = "$firstlabel"
endif
if( "$mtzlabel" == "auto" ) then
  set firstlabel = `awk '$2=="FWT" || $2=="2FOFCWT"{print $2}' ${t}labels.txt | head -n 1`
  if( "$firstlabel" != "" ) set mtzlabel = "$firstlabel"
endif
if( "$mtzlabel" == "auto" ) then
  set firstlabel = `awk '$1=="F"{print $2}' ${t}labels.txt | head -n 1`
  if( "$firstlabel" != "" ) set mtzlabel = "$firstlabel"
endif
set test = `grep $mtzlabel ${t}labels.txt | wc -l`
if( ! $test ) then
  set BAD = "no $mtzlabel in $mtzfile"
  goto exit
endif
echo "mtzlabel = $mtzlabel"


set mtzreso = `awk '/Resolution Range/{getline;getline;print $6}' ${t}mtzdump.txt`

# some new versions of phenix must have a cell
cat ${t}mtzdump.txt |\
awk '/Cell Dimensions :/{getline;getline;\
  a=$1;b=$2;c=$3;al=$4;be=$5;ga=$6;next}\
  /Space group =/{split($0,s,"\047");sg=s[2];next}\
  END{printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s\n",a,b,c,al,be,ga,sg}' \
 >! ${t}.pdb
awk '/^ATOM|^HETAT/{print substr($0,1,80)}' $pdbfile >> ${t}.pdb

# key could be off name or xyz position, hard to say
awk '{x=substr($0,31,8);y=substr($0,39,8);z=substr($0,47,8);\
  printf("(%.3f,%.3f,%.3f) %d KEY\n",x,y,z,++n)}' ${t}.pdb >! ${t}key.txt
phenix.map_value_at_point $mtzfile ${t}.pdb \
          ${phenixlabel}=$mtzlabel scale=sigma |\
tee ${t}debug_rawlog.txt |\
cat ${t}key.txt - |\
awk '$NF=="KEY"{n[$1]=$2;next}\
  /Map value:/ && $3~/^\(/{print "RHO",n[$3],$3,++i,$NF;next}\
  /Map value:/ && /^\"/{rho=$NF;point=substr($0,index($0," Point: "));\
     split(point,w);xyz="(" w[2] w[3] w[4] ")";\
     print "RHO",n[xyz],xyz,++i,$NF}' |\
tee ${t}debug_rhovalues.txt |\
cat - $pdbfile |\
awk '/^RHO/{rho[$4]=rho[$2]=$NF;next}\
  {gsub("\r","")}\
  ! /^ATOM|^HETAT/{print;next}\
  {++n;print $0,"           ",rho[n]}\
  rho[n]==""{print "ERROR: missing rho for atom",n}' |\
cat >! $outfile


set test = `egrep -l "^ERROR" $outfile | wc -l`
if( $test ) then
    egrep "^ERROR" $outfile | head
    set BAD = "density probe failed"
endif


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

