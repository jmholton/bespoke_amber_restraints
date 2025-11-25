#! /bin/tcsh -f
#
#   discard least-needed waters
#
#
set maxreject = 1000

set orignames = orignames.pdb
set Bfac_file = Bfac.pdb
set parmfile = ""
set rstfile  = ""
set notthese = current_restraints.pdb
set outprefix = drier

set mtzfile = cootme.mtz
set mtzlabel = DELFWT

set phenixlabel = miller_array.labels.name
set test = `phenix.version | awk '/Release tag/{print ( $NF < 5000 )}'`
if( "$test" == "1" ) set phenixlabel = label

set quiet = 0
set debug = 0
set tempfile = /dev/shm/${USER}/tempfile_dry$$_

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
      if("$key" == "restraints") set notthese = "$Val"
      if("$key" == "output") set outfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set notthese = $Arg
      if("$Arg" =~ *.mtz ) set mtzfile = $Arg
      if("$Arg" =~ *.rst7 ) set rstfile = $Arg
      if("$Arg" =~ *.parmtop ) set parmfile = $Arg
      if("$Arg" =~ *.prmtop ) set parmfile = $Arg
      if("$Arg" =~ *.parm7 ) set parmfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( "$parmfile" == "" ) set parmfile = xtal.prmtop

if( ! -e "$mtzfile" ) then
    set BAD = "cannot score waters without an mtz file"
    goto exit
endif

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = tempfile_dry_

if(! $quiet) then
cat << EOF
parmfile = $parmfile
rstfile  = $rstfile
mtzfile  = $mtzfile
mtzlabel = $mtzlabel

notthese = $notthese

maxreject = $maxreject

tempfile = $tempfile
debug = $debug
EOF
endif

set t = $tempfile

echo "reading $rstfile"
rm -f ${t}resized.parm7 >& /dev/null
rst2pdb_runme.com parmfile=$parmfile $rstfile outfile=${t}labeled.pdb \
  newparm=${t}resized.parm7 neworig=${t}orignames.pdb >! ${t}cpptraj.log
if( $status ) then
   cat ${t}cpptraj.log
   set BAD = "unable to read $rstfile"
   goto exit
endif

if(-e ${t}resized.parm7) set parmfile = ${t}resized.parm7
if(-e ${t}orignames.pdb) set orignames = ${t}orignames.pdb

echo "extracting water..."
egrep "O  .HOH" ${t}labeled.pdb | egrep "^ATOM|^HETAT" >! ${t}water.pdb

if(-e "$notthese") then
   echo "keeping atoms already restrained in $notthese"
   cat $notthese |\
   awk '/^ATOM|^HETAT/{print $0,"NOTME"}' |\
   cat - ${t}water.pdb |\
   awk '! /^ATOM|^HETAT/{next}\
    {id=substr($0,12,17)}\
    $NF=="NOTME"{++sel[id];next}\
    sel[id]{print}' >! ${t}restrained.pdb
   set lastslot = `tail -n 1 ${t}restrained.pdb | awk '{print $NF}'`

   cat $notthese |\
   awk '/^ATOM|^HETAT/{print $0,"NOTME"}' |\
   cat - ${t}water.pdb |\
   awk -v lastslot=$lastslot '! /^ATOM|^HETAT/{next}\
    {id=substr($0,12,17)}\
    $NF=="NOTME"{++skip[id];next}\
    skip[id] || $NF<=lastslot{next}\
    {print}' >! ${t}disposable.pdb
   set disposable = `cat ${t}disposable.pdb | wc -l`
   echo "$disposable atoms are disposable"
   if( $disposable < $maxreject ) then
      set maxreject = $disposable
      echo "WARNING: unable to reject waters before last restrained water"
   endif
else
   egrep "^ATOM|^HETAT" ${t}water.pdb >! ${t}disposable.pdb
endif

echo "probing map values"
phenix.map_value_at_point $mtzfile ${t}disposable.pdb \
          ${phenixlabel}=$mtzlabel scale=sigma >! ${t}map_values.log

cat ${t}map_values.log ${t}disposable.pdb |\
awk '/Map value:/{++n;rho[n]=$NF;next}\
  ! /^ATOM|^HETAT/{next}\
  {++i;orn=$NF;print rho[i]+0,orn}' |\
sort -g >! ${t}water_rhos.txt
head -n $maxreject ${t}water_rhos.txt |\
awk '$1<0 && NF==2' >! ${t}strip_water.txt 
# rho resnum
set rejecting = `cat ${t}strip_water.txt | wc -l`
if( $rejecting == 0 ) then
    echo "nothing to do."
    cp $orignames ${outprefix}_orignames.pdb
    cp $parmfile ${outprefix}.parm7
    cp $rstfile ${outprefix}.rst7
    cp ref.crd ${outprefix}ref.crd
    goto exit
endif

echo "eliminating lowest-density $rejecting waters"
cat ${t}strip_water.txt |\
cat - ${t}disposable.pdb |\
awk 'NF==2{++sel[$2];sigma[$2]=$1;next}\
  ! /^ATOM|^HETAT/{next}\
  sel[$NF]{print $0,sigma[$NF]}' >! ${t}_stripped.pdb
head -n 1 ${t}_stripped.pdb
tail -n 1 ${t}_stripped.pdb

if(-e "$Bfac_file") then
  echo "also collapsing out of $Bfac_file"
  cat ${t}strip_water.txt |\
  cat - $Bfac_file |\
  awk 'NF==2{++sel[$2];next}\
    ! /^ATOM|^HETAT/{next}\
    ! sel[$NF]{print ++n,substr($0,61,6),"BFAC"}' |\
  cat - $Bfac_file |\
  awk '$NF=="BFAC"{B[$1]=$2;next}\
    ! /^ATOM|^HETAT/{print;next}\
      {++n;pre=substr($0,1,60);post=substr($0,67);\
       printf("%s%6.2f%s\n",pre,B[n],post)}' |\
  cat >! ${outprefix}_Bfac.pdb
endif

touch rejected_water.pdb
cat ${t}_stripped.pdb >> rejected_water.pdb

cat ${t}strip_water.txt |\
awk '{print $NF}' |\
 sort -u | sort -g |\
awk 'NR==1{s=e=$1;next} $1==e+1{e=$1;next} {print s"-"e;s=e=$1} END{print s"-"e}' |\
awk -F "-" '$1==$2{print $1;next} {print}' |\
sort -u | sort -g >! ${t}ranges.txt 
set striprange = `cat ${t}ranges.txt `
set striprange = `echo $striprange | awk '{gsub(" ",",");print}'`
echo "stripping residues: $striprange"

echo "parmed"
rm -f ${outprefix}.parm7 ${outprefix}.rst7 >& /dev/null
parmed -p $parmfile -c $rstfile << EOF >! ${t}parmed.log
strip :$striprange
parmout ${outprefix}.parm7 ${outprefix}.rst7
EOF

# now update reference coords?
if(-e ref.crd) then

  echo "stripping ref.crd too"
  rm -f ${outprefix}ref.crd >& /dev/null
  cpptraj -p $parmfile -y ref.crd << EOF >! ${t}trajref.log
strip :$striprange
trajout ${outprefix}ref.crd
EOF

endif

if( 0 ) then
echo "assembling ${outprefix}.pdb"
egrep -v "HOH|END" ${t}labeled.pdb >! ${outprefix}.pdb

cat ${t}strip_water.txt ${t}labeled.pdb |\
awk 'NF==2{++bad[$2]}\
  ! /^ATOM|^HETAT/{next}\
  ! /HOH/{next}\
  bad[$NF]{next}\
  {print}' |\
cat >> ${outprefix}.pdb

echo "editing $orignames -> ${outprefix}_orignames.pdb"
cat ${t}strip_water.txt $orignames |\
awk 'NF==2{++bad[$2];next}\
  /^CRYST/{print;next}\
  ! /^ATOM|^HETAT/{next}\
  ! /HOH/{print;next}\
  bad[$NF]{next}\
  {print}' |\
cat >! ${outprefix}_orignames.pdb
endif
echo "leaving $orignames with extra padding"


if(-e ref.crd && ! -e ${outprefix}ref.crd) then
  # make one up?
  echo list | cpptraj -p ${outprefix}.parm7 -y ${outprefix}.rst7 >! ${t}countref.log
  
  set natoms = `awk '$4~/^atom/{print $3}' ${t}countref.log`
  echo "$natoms atoms"

  # now reconstruct ref file
  echo "    \n$natoms" >! ${outprefix}ref.crd

  cat ${outprefix}.pdb |\
  awk '! /^ATOM|^HETAT/{next} {++n;\
       X=substr($0,31,8);Y=substr($0,39,8);Z=substr($0,47,8);\
       printf("%12.7f%12.7f%12.7f",X,Y,Z)}\
       n%2==0{print ""}\
       END{if(n%2!=0)print ""}' |\
  cat >> ${outprefix}ref.crd

  awk '/^CRYST1/{printf("%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n",$2,$3,$4,$5,$6,$7);exit}' ${outprefix}.pdb |\
  tail -n 1 >> ${outprefix}ref.crd

  echo "WARNING: be sure to run update_centroid_positions_runme.com on ${outprefix}ref.crd"

endif


if(! $quiet ) then
ls -l ${outprefix}*
endif

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif

exit



egrep "^CRYST1" current_restraints.pdb >! sfallme.pdb
awk '{print substr($0,1,80)}' rejected_water.pdb |\
reformatpdb.awk -v BFAC=80 | tee -a sfallme.pdb

sfall xyzin sfallme.pdb mapout rejects.map << EOF
mode atmmap
symm 1
EOF

