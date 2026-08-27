#! /bin/tcsh -f
#
#  add ions and some water to neutralize pdb in preparation for MD
#
#
set pdbfile = ""
set outfile = salty.pdb
set charge = ""
set concentration = 0.15
set watertype = spce
set bulkwaters = auto

set cation = Na
set anion = Cl

set RIP = 4
set RIW = 3

set RP = 4
set RW = 2.8

set leaplog = "leap.log"

set logfile = debuglog.log

set tempfile = /dev/shm/${USER}/temp_salt_$$_
mkdir -p /dev/shm/${USER}

echo "command-line arguments: $* "

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set key = `echo $Key | awk '{print tolower($0)}'`
    set val = `echo $Val | awk '{print tolower($0)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`

    if( $assign ) then
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "conc") set concentration = "$Val"
    else
      # no equal sign
      if("$key" =~ *.pdb ) set pdbfile = "$Arg"
    endif
    if("$arg" == "debug") set debug = 1
end

if( $?debug ) then
    set tempfile = tempfile
endif
set t = "$tempfile"
#touch $logfile

if( ! -e "$pdbfile" ) then
  set BAD = "pdbfile $pdbfile does not exist"
  goto exit
endif

if(! $?AMBERHOME) source /programs/amber22/amber.csh 
if(! $?src) set src = ${AMBERHOME}/XtalUtilities

if( -e "$leaplog" && "$charge" == "" ) then
  echo "looking in $leaplog"
  set charge = `awk '/unperturbed charge/{gsub(/[)(]/,"");print int($7);exit}' $leaplog`
  echo "charge = $charge"
endif

if( ! -e ${src}/${watertype}.pdb) then
  set BAD = "could not find $watertype water pdb file"
  goto exit
endif





if(-e cation.pdb) then
  echo "using cation.pdb for cations"
else if(-e ${cation}.pdb) then
  cp ${cation}.pdb cation.pdb
else
  phenix.elbow --chemical_component=$cation --opt --amber_force_field_files
  if(-e ${cation}.pdb) cp ${cation}.pdb cation.pdb
endif
if(! -e cation.pdb) then
  set BAD = "cation.pdb not found for $cation — elbow RM1 fails on Phenix 2.0+ (>1000 MB virtual memory); copy ${cation}.pdb to working directory"
  goto exit
endif

if(-e anion.pdb) then
  echo "using anion.pdb for anions"
else if(-e ${anion}.pdb) then
  cp ${anion}.pdb anion.pdb
else
  phenix.elbow --chemical_component=$anion --opt --amber_force_field_files
  if(-e ${anion}.pdb) cp ${anion}.pdb anion.pdb
endif
if(! -e anion.pdb) then
  set BAD = "anion.pdb not found for $anion — elbow RM1 fails on Phenix 2.0+ (>1000 MB virtual memory); copy ${anion}.pdb to working directory"
  goto exit
endif


set anion_charge = `awk '/^ATOM|^HETAT/{sum+=substr($0,79)} END{print sum}' anion.pdb`
set cation_charge = `awk '/^ATOM|^HETAT/{sum+=substr($0,79)} END{print sum}' cation.pdb`
echo "cation charge: $cation_charge"
echo "anion charge:  $anion_charge"

echo "target concentration: $concentration M"

# measure empty space for salt and water
set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
echo | rwcontents xyzin $pdbfile >! ${t}rwcontents.log
grep "% of cell without atoms" ${t}rwcontents.log
#set waterslots = `awk '/Cell volume/{V=$NF} /% of cell without atoms/{print int($NF/100*V*55/10)}' rwcontents.log`
set waterslots = `awk '/Cell volume/{V=$NF} /% of cell without atoms/{print int(55*6.022e23/1e27*V)}' ${t}rwcontents.log`
echo "looks like room for $waterslots waters"

set saltatoms = `echo $concentration $waterslots | awk '{print int($1*$2/55.)}'`
set halfsalt = `echo $saltatoms | awk '{print int($1/2)}' `

set cations = $halfsalt
set anions  = $halfsalt

if( $charge < 0 ) then
  set cations = `echo $charge $cations $cation_charge | awk '{print ($2*$3-$1)*$3}'`
  set anions  = `echo $charge $anions $anion_charge $cations $cation_charge | awk '{print int( -($1+$2*$3+$4*$5)/$3 + $2 )}'`
endif
if( $charge > 0 ) then
  set anions = `echo $charge $anions $anion_charge | awk '{print $2+int(-$1/$3)}'`
endif
set cations = `echo $charge $anions $anion_charge $cations $cation_charge | awk '{print int( $4 - ($1+$2*$3+$4*$5)/$5 )}'`

set newcharge = `echo $charge $anions $anion_charge $cations $cation_charge | awk '{print $1+$2*$3+$4*$5}'`
echo "new charge should be: $newcharge ( $charge protein, $anions anions and $cations cations, salt: $saltatoms"

if( "$bulkwaters" == "auto" ) then
  set bulkwaters = `echo $waterslots $anions $cations | awk '{print $1-($2-$3)}'`
  set bulkwaters = 99999
endif

set P = `awk '/HOH/{print n;exit} /^ATOM|^HETAT/{++n}' $pdbfile`

AddToBox -c $pdbfile -a cation.pdb -na $cations -P $P -RP $RIP -RW $RIW \
 -o ${t}positive.pdb -V 1 |& tee ${t}addcat.log
if(! -e ${t}positive.pdb) then
  set BAD = "AddToBox failed adding cations — check output above"
  goto exit
endif

AddToBox -c ${t}positive.pdb -a anion.pdb -na $anions -P $P -RP $RIP -RW $RIW \
 -o ${t}neutral.pdb -V 1 |& tee ${t}addani.log
if(! -e ${t}neutral.pdb) then
  set BAD = "AddToBox failed adding anions — check output above"
  goto exit
endif

AddToBox -c ${t}neutral.pdb -a ${src}/${watertype}.pdb -na $bulkwaters -P $P -RP $RP -RW $RW \
 -o ${t}wet.pdb |& tee ${t}addwater.log
if(! -e ${t}wet.pdb) then
  set BAD = "AddToBox failed adding waters — check output above"
  goto exit
endif

set added_cations = `awk '/Added/{sum+=$4} END{print sum}' ${t}addcat.log`
set added_anions = `awk '/Added/{sum+=$4} END{print sum}' ${t}addani.log`
set added_waters = `awk '/Added/{sum+=$4} END{print sum}' ${t}addwater.log`
set newcharge = `echo $charge $added_anions $anion_charge $added_cations $cation_charge | awk '{print $1+$2*$3+$4*$5}'`
echo "new charge should be: $newcharge ( $charge protein, $added_anions anions and $added_cations cations, salt: $saltatoms"
echo "$added_waters new waters added"

egrep -v "^END" $pdbfile >! $outfile

echo "appending new atoms to $pdbfile in $outfile"
egrep "^ATOM|^HETAT" $pdbfile |\
wc -l |\
cat - ${t}wet.pdb |\
awk 'NR==1{n=$1;next}\
   /^ATOM|^HETAT/{++i} i>n{print}' |\
convert_pdb.awk -v only=atoms -v fixEe=1 -v OCC=1 -v BFAC=999 >> $outfile

cat << EOF
suggestions: 
reorganize_pdb_runme.com $outfile
EOF

exit:
if($?BAD) then
   echo "ERROR: $BAD"
   exit 9
endif

# clean up
if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?debug && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    rm -f ${tempfile}*
endif


exit


