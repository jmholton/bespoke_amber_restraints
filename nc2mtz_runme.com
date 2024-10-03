#! /bin/tcsh -f
#
# convert nc file to an mtz file
#
#
set smallSG = ""
set reso = 1
set B = 10
set MD_mult = ( 1 1 1 )
set ncfile = ""
set topfile = xtal.prmtop
set paddedparm = padded.parm7

set debug = 0

set keeptraj = 0

set tempfile = /dev/shm/jamesh/temp_$$_traj/

set pdir = `dirname $0`

# read command line
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
    else
      # no equal sign
      if("$arg" =~ [pcifrh][1-6]*) set smallSG = `echo $arg | awk '{print toupper($0)}'`
      if("$key" =~ *.nc ) set ncfile = "$Arg"
    endif
    if("$key" == "debug") set debug = "$Val"
end

set t = $tempfile

set smallSGnum = `awk -v SG=$smallSG '$4 == SG && $1 < 500 {print $1}' $CLIBD/symop.lib | head -1`
set smallSG = `awk -v num=$smallSGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
if("$smallSG" == "") then
    set BAD = "bad space group."
    goto exit
endif

set MD_mult = `echo $MD_mult | awk '{gsub(","," ");print}'`
set nsymops = `awk -v SG=$smallSG '$4 == SG && $1 < 500 {print $2}' $CLIBD/symop.lib | head -1`

cat << EOF
ncfile = $ncfile
smallSG = $smallSG
MD_mult = $MD_mult
tempfile = $t
EOF


mkdir -p ${t}
rm -f trajectory > /dev/null
ln -sf ${t} trajectory
cat << EOF >! ${t}cpptraj.in
image byatom
outtraj trajectory/md.pdb pdb multi pdbv3 keepext sg "P 1"
go
EOF

cat << EOF >! ${t}sfall.in
mode sfcalc xyzin
symm $smallSG
sfsg 1
reso $reso
#breset $B
EOF

again:
echo "cpptraj..."
cpptraj -y $ncfile -p $topfile < ${t}cpptraj.in >&! ${t}cpptraj.log
if($status && -e "$paddedparm" && ! $?RESIZE) then
  set RESIZE
  echo "making new parm file from $paddedparm"
  mv ${t}cpptraj.log ${t}cpptraj_error1.log
  set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");for(i=NF;i>3;--i)if($i+0>0)print $i;exit}' ${t}cpptraj_error1.log | head -n 1`
  if( "$rstatoms" == "" ) then
   set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");print $(NF-2);exit}' ${t}cpptraj_error1.log`
  endif
  set maxatoms = `echo list | cpptraj -p $paddedparm | awk '$4=="atoms,"{print $3}' | head -n 1`
  set stripmask = `echo $rstatoms $maxatoms | awk '{print $1+1"-"$2}'`
  set topfile = ${t}resized.parm7
  echo "new parmfile: $topfile"
  cpptraj -p $paddedparm << EOF >&! ${t}strip1.log
  parmstrip @$stripmask
  parmwrite out $topfile
EOF
  goto again
endif

set pdb = trajectory/md.1.pdb
if( ! -e $pdb ) then
    set BAD = "conversion failed"
    goto exit
endif

head $pdb >! ${t}dummy.pdb
echo $MD_mult |\
cat - ${t}dummy.pdb |\
awk 'NR==1{na=$1;nb=$2;nc=$3}\
  na<1{na=1} nb<1{nb=1} nc<1{nc=1} \
  /^CRYST1/{print $2/na,$3/nb,$4/nc,$5,$6,$7;exit}' |\
cat >! ${t}cell.txt
set CELL = `cat ${t}cell.txt`

pdbset xyzin ${t}dummy.pdb xyzout ${t}cell.pdb << EOF >! ${t}pdbset.log
CELL $CELL
SPACE $smallSG
EOF

setenv MEMSIZE `echo $CELL $reso | awk '{print int($1*$2*$3/($NF**3)*100)}'`

set ns = `ls trajectory/md.*.pdb | awk '{gsub("[^0-9]","");print}' | sort -g`
echo "$#ns pdb files in $ncfile"
set mtzs = ""
echo "sfall"
foreach n ( $ns )
  set pdb = trajectory/md.${n}.pdb
  egrep "^CRYST1" ${t}cell.pdb >! ${t}pdb${n}.pdb
  cat $pdb |\
  awk -v B=$B '! /^ATOM|^HETAT/{next}\
       {RES= substr($0, 18, 3);\
        X =  substr($0, 31, 8);\
        Y =  substr($0, 39, 8);\
        Z =  substr($0, 47, 8);\
       atom=substr($0,12,5);gsub(" ","",atom)\
       Ee = $NF;\
       if(atom=="SE")Ee=atom;\
       #if(Ee!="SE" && Ee!="S" && Ee!="ZN")next;\
       rms+=X*X+Y*Y+Z*Z;\
    printf("ATOM   %4d  %-2s  %3s     1    %8.3f%8.3f%8.3f  1.00%6.2f%12s\n",++n%10000,Ee,RES,X,Y,Z,B,Ee);}\
    END{print "REMARK rms=",rms;print "END"}' |\
  cat >> ${t}pdb${n}.pdb
  set rms = `tail -n 2 ${t}pdb${n}.pdb | awk '{print $NF;exit}'`
  echo "$n $rms" | tee -a ${t}rms.log
  if($rms == 0) then
    set BAD = "all zero coordinates at $n"
    continue
  endif
  set newmtz = ${t}/${n}.mtz
  cat ${t}sfall.in | sfall xyzin ${t}pdb${n}.pdb hklout $newmtz >! ${t}/sfall.${n}.log &
  set mtzs = ( $mtzs $newmtz )
end
wait

addmtzs:
${pdir}/addup_mtzs_runme.com $mtzs


set scale = `echo $nsymops $#ns | awk '{print 1/$1/$2}'`
echo scaling by $scale for avg.mtz
cad hklin1 sum.mtz hklout avg.mtz << EOF >! ${t}scaledown.log
labin file 1 E1=Fsum E2=PHIsum
scale file 1 $scale
labou file 1 E1=FCavg E2=PHICavg
EOF



cleanup:
if( ! $keeptraj ) then
  echo "cleaning up..."
  if(! $debug  ) rm -rf ${t}
  rm -f trajectory
endif

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9
endif

exit


