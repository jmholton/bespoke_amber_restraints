#! /bin/tcsh -f
#
# quick calculator for maximum displacement from start to end of a trajectory
#
#
set ncfile  = ""
set last = ""
set parmfile  = xtal.prmtop
set paddedparm  = padded.parm7
set orignames = orignames.pdb

# will be created from paddedparm and orignames if needed
set newparm = resized.parm7

set parallel = 0
set grep = ""
set quiet = 0
set debug = 0
set tempfile = /dev/shm/${USER}/tempfile_maxdiff_$$_

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
      if("$key" == "output") set outprefix = "$Val"
      if("$key" == "outfile") set pdbfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
      if("$Arg" =~ *.rst7 ) set ncfile = $Arg
      if("$Arg" =~ *.nc ) set ncfile = $Arg
      if("$Arg" =~ *.rst ) set ncfile = $Arg
      if("$Arg" =~ *.crd ) set ncfile = $Arg
      if("$Arg" =~ *.parmtop ) set parmfile = $Arg
      if("$Arg" =~ *.prmtop ) set parmfile = $Arg
      if("$Arg" =~ *.parm7 ) set parmfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $parallel ) goto parallel


if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_maxdiff_
if( $tempfile =~ /dev/shm/$USER/* ) mkdir -p /dev/shm/$USER/

if(! -e "$ncfile") then
   set BAD = "usage: $0 amber.nc [$parmfile] [$orignames]"
   goto exit
endif

set dirname = `dirname $ncfile`
if(! -e "$parmfile") set parmfile = ${dirname}/$parmfile


if(! $quiet) then
cat << EOF
ncfile = $ncfile
parmfile = $parmfile

tempfile = $tempfile
EOF
endif

set t = $tempfile

# just in case xtal.prmtop does not exist
if(! -e "$parmfile" ) then
  echo "WARNING: missing $parmfile"
  if( -e "$newparm") cp "$newparm" "$parmfile"
  if( -e "$paddedparm") cp "$paddedparm" "$parmfile"
endif

if( "$last" == "" ) then
  set last = lastframe
else
  set last = "$last $last"
endif

if(! -e "$paddedparm") then
  ln -s $parmfile  ${t}new.parm7
  goto gotparm
endif

# will always bonk
cpptraj -p $paddedparm << EOF >&! ${t}cpptraj_error1.log
trajin $ncfile 
EOF
set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");for(i=NF;i>3;--i)if($i+0>0)print $i;exit}' ${t}cpptraj_error1.log | head -n 1`
if( "$rstatoms" == "" ) then
  set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");print $(NF-2);exit}' ${t}cpptraj_error1.log`
endif
set maxatoms = `awk '/in associated topology/{gsub("[)(]"," ");print $5;exit}' ${t}cpptraj_error1.log`
if( "$maxatoms" == "") then
  set maxatoms = `echo list | cpptraj -p $paddedparm | awk '$4=="atoms,"{print $3}' | head -n 1`
endif
set stripmask = `echo $rstatoms $maxatoms | awk '{print $1+1"-"$2}'`
set parmfile = $newparm
echo "stripping: $paddedparm"
cpptraj -p $paddedparm << EOF >&! ${t}strip1.log
parmstrip @$stripmask
parmwrite out ${t}new.parm7
EOF
gotparm:
echo "extracting non-H atoms from first and $last frames"
cpptraj -p ${t}new.parm7 << EOF >&! ${t}cpptraj.log
trajin $ncfile 1 1
trajin $ncfile $last
strip @H=
trajout ${t}first.pdb onlyframes 1
trajout ${t}last.pdb onlyframes 2
EOF
if( $status ) then
  set BAD = "conversion failed"
  goto exit
endif

set nframes = `awk '/trajectory with coordinates/{print $NF+0;exit}' ${t}cpptraj.log`

foreach name ( first last )
egrep "$grep" ${t}${name}.pdb |\
awk '/^ATOM/{xyz=substr($0,30);printf("ATOM %21d   %s\n",++n,xyz)}' |\
  cat >! ${t}${name}A.pdb
end
echo "rms, max-distance-traveled #frames atom#  atomid:"
rmsd ${t}firstA.pdb ${t}lastA.pdb | egrep "MAXD.all|RMSD.all" >! ${t}maxd.txt
set maxd = `awk '{print $2}' ${t}maxd.txt`
set n = `awk '{print $NF}' ${t}maxd.txt | tail -n 1`

echo -n "$maxd   $nframes    $n  "
awk -v n=$n '! /^ATOM/{next} {++i} i==n{print substr($0,12,45);exit}' ${t}last.pdb


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


parallel:

foreach nc ( `ls -1rt amber_*[0-9].nc` )

set basename = `basename $nc .nc`
if(-e maxd_${basename}.txt) continue
echo -n "$nc  "
srun $0 $nc last="$last" grep="$grep" >! maxd_${basename}.txt &
sleep 0.2

end
wait

rm -f opts_vs_runme.txt
set ns = `ls -1rt runme*.log | awk -F "/" '{print substr($NF,6)+0}'`
foreach n ( $ns )
set opts = `awk '/^running|^scratch/{exit} $2=="="{print $1 $2 $3}' runme${n}.log`
set itrs = `awk '/^running amber/{print substr($2,7)}' runme${n}.log | awk 'NR==1{print} END{print}'`
if( $#itrs != 2 ) continue
echo $n $itrs $opts | tee -a opts_vs_runme.txt
end


tail -n 1 maxd*[0-9].txt |\
 awk '/^==/{f=substr($2,12);getline;print f+0,$0}' |\
sort -g |\
tee maxd_vs_itr.txt

cat opts_vs_runme.txt |\
awk '{for(i=1;i<=NF;++i)if($i~/allatom_weigh/){split($i,w,"=");print $2,w[2],"opt"}}' |\
cat - nwaters_vs_itr.txt |\
awk '{s=$1-0.1;print s,$0}' |\
cat - maxd_vs_itr.txt |\
sort -g |\
awk '$4=="opt"{opt=$3;next}\
  NF==5{nwater[$2]=$3;next}\
  nwater[$1]==""{nwater[$1]=nwater[$1-1]}\
  {print $1,$2,$3,$4,opt+0,nwater[$1]}' |\
tee maxd_vs_allatomw.txt

exit


foreach o ( `seq 82 -1 1` )

  if( ! -e ../opt$o ) continue
  if(-e ../opt${o}/maxd_vs_itr.txt) continue
  cd ../opt$o
  pwd
  rm -f maxd_amber*
  traj_maxd_runme.com parallel=1 last=25 >& /dev/null
  tail maxd_vs_itr.txt

end



rm ../maxd_vs_aw.txt
foreach o ( `seq 82 -1 1` )

  if( ! -e ../opt$o ) continue
  cd ../opt$o
  pwd
foreach aw ( `awk '{print $5}' maxd_vs_allatomw.txt | sort -u` )
  set medmadW = `awk -v aw=$aw '$5==aw{print $3}' maxd_vs_allatomw.txt | median.awk`
  set medmadN = `awk -v aw=$aw '$5==aw{print $6}' maxd_vs_allatomw.txt | median.awk`
  if( $#medmadN != 3 || $#medmadW != 3 ) continue
  echo "$aw $medmadW $medmadN $o" | tee -a ../maxd_vs_aw.txt
end

end


end

