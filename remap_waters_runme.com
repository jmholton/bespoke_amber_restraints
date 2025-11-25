#! /bin/tcsh -f
#
#  reorganize all waters, moving restrained and high-density ones to top     -James Holton 8-15-24
#  stil need to move B factors to follow
#
#
set parmfile = xtal.prmtop
set rstfile = ""
set restraints = current_restraints.pdb
set orignames  = orignames.pdb
set Bfac_file  = Bfac.pdb

set remap_restraints = all
set remap_density = all

set outprefix = remapped

set mtzfile = cootme.mtz
set mtzlabel = DELFWT

set quiet = 0
set debug = 0

set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"

if(! $?AMBERHOME) then
  set BAD = "need amber set up"
  goto exit
endif

set phenixlabel = miller_array.labels.name
set test = `phenix.version | awk '/Release tag/{print ( $NF < 5000 )}'`
if( "$test" == "1" ) set phenixlabel = label

set tempfile = /dev/shm/${USER}/tempfile_remap_$$_

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
      if("$key" == "current") set restraints = "$Val"
      if("$key" == "F") set mtzlabel = "$Val"

      if("$key" == "output") set outprefix = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) then
         set restraints = $Arg
      endif
      if("$Arg" =~ *.mtz ) then
         set mtzfile = $Arg
      endif
      if("$Arg" =~ *.rst7 ) then
         set rstfile = $Arg
      endif
    endif
    if("$key" == "debug") set debug = "$Val"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_remap_

if(! -e "$rstfile" ) then
    set BAD = "no coordinates provided"
    goto exit
endif
if(! -e "$parmfile") then
    set BAD = "no parmtop provided"
    goto exit
endif

if("$outprefix" =~ *.rst7) set outprefix = `basename $outprefix .rst7`
if("$outprefix" =~ *.pdb) set outprefix = `basename $outprefix .pdb`
if("$outprefix" =~ *.parm7) set outprefix = `basename $outprefix .parm7`

set t = ${tempfile}

#if(! $debug) then
#   set pwd = `pwd`
#   mkdir -p $tempdir
#   cp $mtzfile $parmfile $rstfile $refpoints $tempdir
#   cd $tempdir
#endif
if(! $quiet) then
cat << EOF
rstfile = $rstfile
parmfile = $parmfile
orignames = $orignames
restraints = $restraints
mtzfile = $mtzfile
mtzlabel = $mtzlabel

tempfile = $tempfile
outprefix = $outprefix
EOF
endif



echo "reading xyz from $rstfile"
rst2pdb_runme.com parmfile=$parmfile orignames=$orignames \
   $rstfile outfile=${t}labeled.pdb  >! ${t}cpptraj.log
grep WARN ${t}cpptraj.log
if(! $status) then
   set BAD = "unable to read $rstfile"
   goto exit
endif

# need atom ids for residue-based remapping
cat ${t}labeled.pdb |\
awk '/^ATOM|^HETAT/{++a;atoms[$NF]=atoms[$NF]" "a}\
  END{for(r in atoms)print "ATOMS",r,atoms[r]}' |\
sort -k2g >! ${t}atommap.txt

egrep "O  .HOH" ${t}labeled.pdb |\
tee ${t}water.pdb |\
awk '/^ATOM|^HETAT/{id=substr($0,17,12);ordresnum=$NF;\
  print id "|",ordresnum,"| ORDER"}' |\
cat >! ${t}id_order.txt


# recover original names of reference points
if(-e "$restraints") then
  echo "extracting water from $restraints"
  egrep "O  .HOH" $restraints >! ${t}restrained.pdb

  # double check that restrained waters exist
  cat ${t}water.pdb |\
  awk -F "|" '{print $0 "| EXISTS"}' |\
  cat - ${t}restrained.pdb |\
  egrep "O  .HOH" |\
  awk '{id=substr($0,12,17)}\
    $NF=="EXISTS"{++exists[id];next}\
    ! /^ATOM|^HETAT/{next}\
    exists[id]{print}' |\
  cat >! ${t}exists.pdb
  diff ${t}restrained.pdb ${t}exists.pdb
  if( $status ) echo "WARNING: above restrained waters from $restraints dont exist, ignoring them."
  cp ${t}exists.pdb ${t}restrained_water.pdb
else
  touch ${t}restrained_water.pdb
endif

cat ${t}restrained_water.pdb |\
awk -F "|" '{print $0 "| TAKEN"}' |\
cat - ${t}water.pdb |\
egrep "O  .HOH" |\
awk '{id=substr($0,12,17)}\
  $NF=="TAKEN"{++taken[id];next}\
  ! /^ATOM|^HETAT/{next}\
  ! taken[id]{print}' |\
cat >! ${t}probeme.pdb

if(-e "$mtzfile") then
  echo "probing map: $mtzlabel in $mtzfile"
  phenix.map_value_at_point $mtzfile ${t}probeme.pdb \
            ${phenixlabel}=$mtzlabel scale=sigma >! ${t}map_values_ref.log

  cat ${t}map_values_ref.log ${t}probeme.pdb |\
  awk '/Map value:/{++n;rho[n]=$NF;next}\
    /^CRYST/{print;next}\
    ! /^ATOM|^HETAT/{next}\
    {++m;id=substr($0,12,17);\
     print rho[m],$0,rho[m]}' |\
  sort -gr |\
  awk '{print substr($0,index($0,$2))}' |\
  cat >! ${t}sorted_unrestrained.pdb
  echo "top and bottom:"
  head -n 1 ${t}sorted_unrestrained.pdb
  tail -n 1 ${t}sorted_unrestrained.pdb
else
  cp ${t}probeme.pdb ${t}sorted_unrestrained.pdb
endif


echo "assigning new names"
cat ${t}id_order.txt ${t}restrained_water.pdb ${t}sorted_unrestrained.pdb |\
awk -F "|" '$NF~/ORDER$/{++o;newid[o]=$1;resnum[$1]=$2+0;next}\
  ! /^ATOM|^HETAT/{next}\
  {++n;id=substr($0,17,12);\
   print id "|" newid[n] "|", resnum[id],"|",resnum[newid[n]],"| RENAME";}' |\
cat >! ${t}newnames.txt
# format: oldid newid oldordresnum newordresnum
echo "first and last:"
echo "   old id       new id      old      new"
head -n 1 ${t}newnames.txt
tail -n 1 ${t}newnames.txt

cp ${t}newnames.txt ${outprefix}_pairs.txt


if(-e "$Bfac_file") then
  echo "re-mapping water B factors in $Bfac_file"
  awk '/^ATOM|^HETAT/ && /HOH/{print $NF,"|",substr($0,61,6),"| BFAC"}' $Bfac_file |\
  cat - ${t}newnames.txt |\
  awk -F "|" '/BFAC$/{Bfac[$1+0]=$2;next}\
    /RENAME/ && Bfac[$3+0]!=""{print $4,Bfac[$3+0],"NEWB"}' |\
  cat - $Bfac_file |\
  awk '$NF=="NEWB"{newB[$1]=$2;next}\
    ! /^ATOM|^HETAT/ || ! /HOH/{print;next}\
    newB[$NF]==""{print "REMARK WARNING: no new B for",$NF}\
    newB[$NF]==""{newB[$NF]=substr($0,61,6)}\
    {pre=substr($0,1,60);post=substr($0,67);\
     printf("%s%6.2f%s\n",pre,newB[$NF],post)}' |\
  cat >! ${outprefix}_Bfac.pdb
  echo "recommend: cp ${outprefix}_Bfac.pdb $Bfac_file"
endif


if(! -e "$restraints") goto renameorig
echo "re-naming waters in $restraints -> ${outprefix}_restraints.pdb"
egrep -v "HOH|TER|END" $restraints >! ${outprefix}_restraints.pdb
cat ${t}newnames.txt ${t}restrained_water.pdb |\
awk -F "|" '/RENAME$/{id=$1;newid[id]=$2;ordres[id]=$4;next}\
  ! /HOH/{next}\
  {id=substr($0,17,12);pre=substr($0,1,16);post=substr($0,29)}\
  newid[id]==""{print "REMARK WARNING: new id for ",id,"missing"; newid[id]=id;next}\
  {print ordres[id]+0, pre newid[id] post}' |\
sort -g |\
awk '{print substr($0,index($0,$2))}' |\
cat >> ${outprefix}_restraints.pdb
grep WARN ${outprefix}_restraints.pdb
if(! $status) then
    set BAD = "missing ids in $restraints"
    goto exit
endif


renameorig:
if(1 || ! -e "$orignames") goto renamerst
set nres  = `awk '/^ATOM|^HETAT/{print substr($0,17,12)}' ${t}labeled.pdb | sort -u | wc -l`
echo "re-naming waters in $orignames -> ${outprefix}_orignames.pdb"
egrep -v "HOH|TER|END" $orignames >! ${outprefix}_orignames.pdb
grep "HOH" $orignames |\
awk -v nres=$nres '/^ATOM|^HETAT/ && $NF<=nres' |\
cat ${t}newnames.txt - |\
awk -F "|" '/RENAME$/{id=$1;newid[id]=$2;ordres[id]=$4;next}\
  ! /HOH/{print;next}\
  {id=substr($0,17,12);pre=substr($0,1,16);post=substr($0,29,80-29)}\
  newid[id]==""{print "REMARK WARNING: new id for ",id,"missing"; newid[id]=id}\
  {print ordres[id]+0, pre newid[id] post,"          ",ordres[id]+0}' |\
sort -g |\
awk '{print substr($0,index($0,$2))}' |\
cat >> ${outprefix}_orignames.pdb
grep WARN ${outprefix}_orignames.pdb
if(! $status) then
   set BAD = "missing ids in ${outprefix}_orignames.pdb"
   goto exit
endif

renamerst:
echo "moving water xyzs in $rstfile -> ${outprefix}.rst7"
set natoms = `tail -n 1 ${t}atommap.txt | awk '{print $NF}'`
cat ${t}newnames.txt |\
awk -F "|" '{print $3,$4}' |\
sort -g |\
cat ${t}atommap.txt - |\
awk -v natoms=$natoms '/^ATOM/{r=$2;na[r]=NF-2;\
  for(i=1;i+2<=NF;++i){atom[r,i]=$(i+2)};next}\
  {old=$1;new=$2;if(na[old]!=na[new]){next}\
    for(i=1;i<=na[new];++i){\
       print atom[new,i],atom[old,i];\
       ++seen[atom[new,i]]}}\
 END{for(a=1;a<=natoms;++a)if(! seen[a])print a,a}' |\
sort -g |\
awk '{print $2}' >! ${t}remap.dat

rm -f ${outprefix}.rst7 >& /dev/null
cpptraj -p $parmfile -y $rstfile << EOF >! ${t}trajremap.log
readdata ${t}remap.dat name map
remap data map parmout ${outprefix}.parm7
trajout ${outprefix}.rst7
EOF
if( $status ) then
   set BAD = "unable to remap $rstfile"
   goto exit
endif
echo "should be same: $parmfile == ${outprefix}.parm7"

echo "reading back ${outprefix}.rst7 -> ${outprefix}.pdb"
rst2pdb_runme.com ${outprefix}.rst7 orignames=$orignames  >! ${t}cpptraj.log
if( $status ) then
   set BAD = "unable to read back ${outprefix}.rst7"
   goto exit
endif
grep WARN ${t}cpptraj.log
if(! $status) then
   set BAD = "corrupted orignames"
   goto exit
endif


if(-e "$restraints") then
   echo "biggest restraint challenge:"
   rmsd ${outprefix}.pdb ${outprefix}_restraints.pdb | egrep "MAXD.all"
   echo "previous:"
   rmsd ${t}labeled.pdb $restraints | egrep "MAXD.all"
endif


goto exit


exit:

#if($debug) wc $rstfile ${outprefix}.rst7

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif

exit

#######################################################################################################
#  notes and usage
#

set itr = `ls -1rt amber_*[0-9].rst7 | tail -n 1 | awk -F "_" '{print $NF+0}'`
rst2pdb_runme.com amber_${itr}.rst7 debug=1

remap_waters_runme.com amber_${itr}.rst7 restraints_for_${itr}.pdb | tee remap_water_${itr}.log



