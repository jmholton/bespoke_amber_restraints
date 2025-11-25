#! /bin/tcsh -f
#
#                                                           - James Holton 1-16-24
#  add water and strip-to-fit residue count
#  updating xtal.prmtop , ref.crd and orignames.pdb
#
set pdbfile = ""

set outpdb = wetter.pdb

set rstfile = toodry.rst7
set outrst = wetter.rst7

set paddedparm = padded.parm7

set outtop = xtal.prmtop
set orignames = orignames.pdb

set restraints = current_restraints.pdb


set protein_radius = 4
set water_radius = 2.8
set centroid_radius = 6

set nadd = 1000000
set force_nadd = 0

set minimize = 0
set energycheck = 0


set debug = 0

set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"

set tempfile = /dev/shm/${USER}/temp_hydrate_$$_


set itr = 0

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
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
      if("$key" == "debug") set debug = "$Val"
      if("$key" == "topout") set outtop  = "$Val"
      if("$key" == "topout") set outtop  = "$Val"
      if("$key" == "output") set outrst  = "$Val"
      if("$key" == "outfile") set outrst  = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) then
         set pdbfile = $Arg
      endif
      if("$Arg" =~ *top || "$Arg" =~ *.parm7 ) then
         set topfile = $Arg
      endif
      if("$Arg" =~ *.rst7 || "$Arg" =~ "*.crd" ) then
         set rstfile = $Arg
      endif
      if("$arg" == "debug") set debug = 1
    endif
end

# fix anything invalid
set test = `echo $nadd | awk '{print int($1)*($1>0)}'`
if( "$test" != "$nadd" ) then
  echo "WARNING: changing nadd=$nadd to nadd=$test"
  set nadd = $test
endif

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_hyd_
set tmpdir = `dirname $tempfile`
mkdir -p $tmpdir

set t = ${tempfile}

if(! -e "$pdbfile" && ! -e "$rstfile" ) then
    set BAD = "no coordinates provided"
    goto exit
endif

if( $nadd == 0 ) then
  echo "nothing to do. leaving $outtop and ref.crd unchanged."
  cp $pdbfile $outpdb
  goto exit
endif

if(! -e "$pdbfile" && -e "$rstfile" ) then
  rm -f ${t}.pdb
  rst2pdb_runme.com $rstfile ${t}.pdb >! ${t}rst2pdb.log
  if( $status || ! -e ${t}.pdb ) then
    set BAD = "could not generate pdb file"
    goto exit
  endif
  set pdbfile = ${t}.pdb
endif





if( ! $?AMBERHOME ) then
source /programs/amber22/amber.csh
endif
set src = ${AMBERHOME}/XtalUtilities
set pdir = /home/jamesh/projects/amber/1aho_refine
set path = ( $path $pdir )
set pwd = `pwd`


echo "counting protein atoms in $pdbfile"
set protein_atoms = `convert_pdb.awk -v only=protein,atoms $pdbfile | wc -l`
echo "$protein_atoms"

echo "extracting last water from $pdbfile"
tac $pdbfile |\
 awk '/^TER/{++ter} ter>1{exit}\
 ! /^ATOM|^HETAT/{next}\
 {atom=substr($0,12,5)} seen[atom]{exit}\
 {print;++seen[atom]}' |\
 tac |\
awk '/^ATOM|^HETAT/{print substr($0,1,80)}' |\
 tee ${t}onewater.pdb


echo "AddToBox"
set force = ""
if( $force_nadd ) set force = " -V 1" 
AddToBox -c $pdbfile -a ${t}onewater.pdb -na $nadd -P $protein_atoms -RP $protein_radius -RW $water_radius \
-o ${t}wet.pdb $force >&! ${t}addwater.log

set added = `awk '/Added/{sum+=$4} END{print sum}' ${t}addwater.log`
echo "added $added"

# check for corruption?

egrep "^CRYST1" $pdbfile | head -n 1 >! ${t}importme.pdb
convert_pdb.awk -v fixEe=1 -v OCC=1 \
  -v renumber=ordinal,watS,w4,chain,chainrestart -v append=ordresnum \
   ${t}wet.pdb | egrep "^ATOM|^HETAT" | tee -a ${t}importme.pdb | tail

cp ${t}importme.pdb $outpdb

echo "hoping that $orignames is long enough"
#cp $outpdb $orignames


if(! -e "$paddedparm" ) then
  set BAD = "need padded parm file to add atoms"
  goto exit
endif

echo "grafting into current system"
graft_atoms_runme.com $outpdb $rstfile paddedparm=$paddedparm \
  minimize=$minimize energycheck=$energycheck \
  outtop=$outtop outrst=$outrst debug=$debug
if( $status ) then
  set BAD = "grafting failed"
  goto exit
endif



if(! -e "$restraints") goto exit

echo "updating ref.crd"
cp -p ref.crd prevrf.crd >& /dev/null
rm -f ref.crd ref.rst7 > /dev/null
combine_pdbs_runme.com $restraints $outpdb \
      printref=1 \
      outfile=${t}refpoints_in_sys.pdb >! ${t}mkref.log
  egrep "^CRYST|^ATOM|^HETAT|^SSBO|^LINK" ${t}refpoints_in_sys.pdb |\
  awk '{print substr($0,1,80)}' >! ${t}refpoints.pdb

cpptraj -p $outtop -y ${t}refpoints.pdb << EOF >> ${t}mkref.log
trajout ${t}ref.rst7
EOF
if( $status ) then
  cat ${t}mkref.log
  set BAD = "updating ref.crd failed"
  goto exit
endif
mv ${t}ref.rst7 ref.crd



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

