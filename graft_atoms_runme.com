#! /bin/tcsh -f
#
#  Graft a PDB file into an MD simulation - potentially changing number of waters     - James Holton 1-3-24
#
#  change the number of atoms in the simulation by adding/dropping waters at the end
#  copy available velocities from previous run
#  not yet updating orignames.pdb? Bfac.pdb ?
#
#
set topfile = xtal.prmtop
set rstfile = wrapped.rst7 
set pdbfile = orignames.pdb

set paddedparm = padded.parm7

set minimize = 1
set energycheck = 1

set outtop  = grafted.prmtop
set outrst  = grafted.rst7
set outpdb  = ""

set quiet = 0
set debug = 0

set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"


set tempfile = /dev/shm/${USER}/tempfile_graft_$$_

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
      if("$key" == "debug") set debug = "$Val"
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

if( $minimize ) set debug = 1

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_ga_
set tmpdir = `dirname $tempfile`
mkdir -p $tmpdir


if(! -e "$pdbfile" ) then
    set BAD = "no pdb file provided"
    goto exit
endif
if(! -e "$rstfile" ) then
    set BAD = "need previous rst7 file to preserve velocities"
    goto exit
endif
if(! -e "$topfile") then
    set BAD = "no parmtop provided"
    goto exit
endif

if( ! $quiet ) then
cat << EOF
topfile = $topfile
rstfile = $rstfile
pdbfile = $pdbfile

paddedparm = $paddedparm

minimize = $minimize
energycheck = $energycheck

outtop  = $outtop
outrst  = $outrst

set tempfile = $tempfile
set debug = $debug
EOF
endif

set t = ${tempfile}

#if(! $debug) then
#   set pwd = `pwd`
#   mkdir -p $tempdir
#   cp $mtzfile $topfile $rstfile $refpoints $tempdir
#   cd $tempdir
#endif

if(-e "$paddedparm") then
  set maxatoms = `echo list | cpptraj -p $paddedparm | awk '$4=="atoms,"{print $3}' | head -n 1`
else
  set maxatoms = 0
endif


# make sure this is text
echo "checking $topfile with $rstfile"
cpptraj -p $topfile -y $rstfile << EOF >&! ${t}cpptraj.log
trajout ${t}old.rst7
EOF
if( $status ) then
  echo "making new parm file from $paddedparm"
  mv ${t}cpptraj.log ${t}cpptraj_error1.log
  set topfile = ${t}.parm7
  set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");for(i=NF;i>3;--i)if($i+0>0)print $i;exit}' ${t}cpptraj_error1.log | head -n 1`
  if( "$rstatoms" == "" ) then
    set BAD = "unable to count atoms in $rstfile"
    goto exit
  endif
  set stripmask = `echo $rstatoms $maxatoms | awk '{print $1+1"-"$2}'`
  cpptraj -p $paddedparm << EOF >&! ${t}strip1.log
  parmstrip @$stripmask
  parmwrite out $topfile
EOF
  cpptraj -p $topfile -y $rstfile << EOF >&! ${t}cpptraj.log
trajout ${t}old.rst7
EOF

endif
if( ! -e ${t}old.rst7 ) then
    set BAD = "$topfile and $rstfile are incompatible"
endif
set oldatoms = `awk '$4=="atoms," && $3+0>0{print $3}' ${t}cpptraj.log | head -n 1`
set oldres = `awk '$6=="res," && $5+0>0{print $5}' ${t}cpptraj.log | head -n 1`
set simtime = `head -n 2 ${t}old.rst7 | tail -n 1 | awk '{print $2}'`
echo "$oldatoms atoms, $oldres residues, $simtime ps"
if( "$oldatoms" == "" || "$oldres" == "" ) then
  set BAD = "cannot read $rstfile"
endif


# extract velocities
head -n -1 ${t}old.rst7 |\
awk 'NR==2{atoms=$1}\
  NR>2{for(i=0;i<=5;++i){f=substr($0,i*12+1,12);\
   if(length(f)==12){++n;if(n>atoms*3 && n<=atoms*6)print f}}}' |\
cat  >! ${t}vel.txt
set nvel = `cat ${t}vel.txt | wc -l | awk '{print $1/3}'`
echo "extracted $nvel velocities from $rstfile"
tail -n 1 ${t}old.rst7 >! ${t}box.txt


cat $pdbfile |\
awk '! /^ATOM|^HETAT/{next}\
  {++natoms}\
  {resid=substr($0,22,9)}\
  resid!=last{last=resid;++nres}\
  END{print natoms,nres}' >! ${t}nres.txt
set natoms = `awk '{print $1}' ${t}nres.txt`
set nres = `awk '{print $2}' ${t}nres.txt`
echo "$natoms atoms and $nres residues in $pdbfile"



if( $natoms == $oldatoms ) then
  if("$topfile" != "$outtop") cp $topfile $outtop
  cp ${t}vel.txt ${t}newvel.txt
  goto paste
endif

if( $natoms < $oldatoms ) then
  echo "stripping residues from the end"
set stripmask = `echo $nres $oldres | awk '{print $1+1"-"$2}'`
cpptraj -p $topfile << EOF >&! ${t}strip.log
parmstrip :$stripmask
parmwrite out $outtop
EOF

  # trim the velocities
  set nnum = `echo $natoms 3 | awk '{print $1*$2}'`
  head -n $nnum ${t}vel.txt >! ${t}newvel.txt

  goto paste
endif



if(! -e "$paddedparm") then
    set BAD = "no padded parmtop provided"
    goto exit
endif

echo "stripping out padding from $paddedparm"

set nresmax = `echo list | cpptraj -p $paddedparm | awk '$11=="res,"{print $10}' | head -n 1`

if( $nres > $nresmax ) then
   set BAD = "not enough padding: need $nres residues"
   goto exit
endif

set stripmask = `echo $nres $nresmax | awk '{print $1+1"-"$2}'`
cpptraj -p $paddedparm << EOF >&! ${t}strip.log
parmstrip :$stripmask
parmwrite out $outtop
EOF
set newatoms = `awk '$4=="atoms,"{print $3}' ${t}strip.log | tail -n 1`
set newres = `awk '$6=="res,"{print $5}' ${t}strip.log | tail -n 1`
echo "now there are $newatoms atoms and $newres residues in $outtop"

# append zero velocities
set newvel = `echo $newatoms $oldatoms | awk '{print ($1-$2)}'`
echo "appending $newvel zero velocities to list"
cp ${t}vel.txt ${t}newvel.txt
echo $newvel |\
awk '{for(i=1;i<=$1*3;++i)printf("%12.7g\n",0)}' |\
cat >> ${t}newvel.txt


paste:
echo "converting $pdbfile into rst7 format"
cpptraj -p $outtop -y $pdbfile << EOF >>& ${t}strip.log
trajout ${t}xyz.rst7
EOF

if( $minimize == 0 ) then
  cp ${t}xyz.rst7 ${t}Min.rst7
  goto skipmin
endif


cat << EOF >! ${t}Cpu.in
Minimize
 &cntrl
  imin=1,
  ntx=1,
  irest=0,
  ntf=1,
  ntc=1,
  maxcyc=10,
  ncyc=0,
  ntpr=1,
  ntwr=10,
  ntwx=0,
  cut=9,
  nsnb=1,
  ntr=1,
  restraintmask=':*&!@H=',
  restraint_wt=10,
 /
EOF

cat << EOF >! ${t}Min.in
Minimize
 &cntrl
  imin=1,
  ntx=1,
  irest=0,
  ntf=1,
  ntc=1,
  maxcyc=100,
  ncyc=0,
  ntpr=1,
  ntwr=100,
  ntwx=0,
  cut=9,
  nsnb=1,
  ntr=1,
  restraintmask=':*&!@H=',
  restraint_wt=10,
 /
EOF

set amber = sander.OMP
set laStage = xyz
rst2pdb_runme.com ${t}${laStage}.rst7 $outtop > /dev/null

foreach Stage ( Cpu Min )

echo "minimizing with $amber ..."
$amber -O -i ${t}${Stage}.in -o ${t}${Stage}.out \
   -p ${outtop} \
   -c ${t}${laStage}.rst7 \
   -ref ${t}${laStage}.rst7 \
   -r ${t}${Stage}.rst7 \
   -x ${t}${Stage}.nc \
   -inf ${t}${Stage}.mdinfo

  cat ${t}${Stage}.out |\
  awk '/  FINAL RESULTS/{++p;next}\
     / A V E R A G E S /{++p;next}\
     / F L U C T U A T /{p=0;next}\
     p && NF>1{print}\
     /EAMBER|---------------/{p=0}' |\
    tee ${t}junk.txt

echo "$laStage vs $Stage"
rst2pdb_runme.com ${t}${Stage}.rst7 $outtop > /dev/null
rmsd ${t}${laStage}.pdb ${t}${Stage}.pdb | egrep -v Bfac

set amber = "$pmemd"
set laStage = $Stage

end

echo "xyz vs Min"
rmsd ${t}xyz.pdb ${t}Min.pdb | grep -v Bfac
#rmsd -v debug=1 ${t}minimized.pdb ${t}moved.pdb | grep moved | sort -k1.25g | tail 

# check with gemmi
echo "overall closest contact:"
gemmi contact -d 2 --sort ${t}Min.pdb >! ${t}gemmi_final.log
head -n 3 ${t}gemmi_final.log

skipmin:

# make sure its text and ends at last atom xyz
echo "converting to $outrst"
rm -f ${t}hybrid.rst7
cpptraj -p $outtop -y ${t}Min.rst7 << EOF >! ${t}cpptraj_newxyz.log
trajout ${t}hybrid.rst7 nobox time0 $simtime
EOF
if($status || ! -e ${t}hybrid.rst7 ) then
  set BAD = "xyz conversion error"
  goto exit
endif

# preserve simulation time stamp in header
#head -n 2 ${t}old.rst7 >! ${t}hybrid.rst7
#tail -n +3 ${t}newxyz.rst7 >> ${t}hybrid.rst7
#set lines = `wc -l ${t}newxyz.rst7 | awk '{print $1+1}'`

echo "appending velocities from $rstfile"
cat ${t}newvel.txt |\
awk '{printf("%s",$0)}\
  {++n} n%6==0{print ""}\
  END{if(n%6!=0)print ""}' |\
cat >> ${t}hybrid.rst7
echo "adding box"
cat ${t}box.txt >> ${t}hybrid.rst7
#tail -n 1 ${t}hybrid.rst7

echo "testing..."
rm -f $outrst >& /dev/null
cpptraj -p $outtop -y ${t}hybrid.rst7 << EOF >! ${t}cpptraj_test.log
trajout $outrst
go
EOF
if($status || ! -e ${outrst}) then
  set BAD = "failed final test"
  goto exit
endif


echeck:
if( $energycheck == 0 ) goto exit
# final energy check
echo "energy check:"
cpptraj.OMP -p $topfile -y $rstfile << EOF >! ${t}cpptraj_final.log
energy e out ${t}energy.dat 
EOF
awk '{print "pre:",$0}' ${t}energy.dat
cpptraj.OMP -p $outtop -y $outrst << EOF >! ${t}cpptraj_final.log
energy e out ${t}energy.dat 
go
EOF
awk '{print "now:",$0}' ${t}energy.dat | tail -n 1

wc $rstfile $outrst

#set natoms = `egrep "^ATOM|^HETAT" $pdbfile | wc -l`
#set oatoms = `egrep "^ATOM|^HETAT" orignames.pdb | wc -l`
#if("$pdbfile" != "orignames.pdb" && $natoms > $oatoms ) then
#   echo "suggestions:"
#   echo "cp $pdbfile orignames.pdb"
#endif

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



