#! /bin/tcsh -f
#
#  move provided waters into vacant green density      -James Holton 1-5-23
#
#  can also provide current and complete refpoints list, which will try to fill vacant refpoints
#
#
set topfile = xtal.prmtop
set rstfile = ""
set pdbfile = ""
set notthese = current_restraints.pdb
set fullref = all_possible_refpoints.pdb
set teleportee = ""

set maxmoves = 10
set minrho = 0
set mindist = 2.0

set minimize = 1
set energycheck = 1

set outfile = teleported.rst7
set outref = teleported_restraints.pdb

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

set tempfile = /dev/shm/${USER}/tempfile_wtp_$$_

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
      if("$key" == "refpoints") set fullref = "$Val"
      if("$key" == "current") set notthese = "$Val"
      if("$key" == "maxmove") set maxmoves = "$Val"

      if("$key" == "output") set outfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) then
         set pdbfile = $Arg
      endif
      if("$Arg" =~ *.rst7 ) then
         set rstfile = $Arg
      endif
    endif
    if("$key" == "debug") set debug = "$Val"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./teletemp_
set tmpdir = `dirname $tempfile`
if("$tmpdir" != "") mkdir -p "$tmpdir"

if(! -e "$rstfile" || ! -e "$fullref" ) then
    set BAD = "no coordinates provided"
    goto exit
endif
if(! -e "$topfile") then
    set BAD = "no parmtop provided"
    goto exit
endif
if(! -e "$mtzfile") then
    set BAD = "no mtz provided"
    goto exit
endif
if("$teleportee" != "" && ! -e "$teleportee") then
  echo "WARNING: $teleportee does not exist"
endif


set t = ${tempfile}

#if(! $debug) then
#   set pwd = `pwd`
#   mkdir -p $tempdir
#   cp $mtzfile $topfile $rstfile $refpoints $tempdir
#   cd $tempdir
#endif

set smallSG = `echo head | mtzdump hklin $mtzfile | awk '/Space group/{print $NF+0}'`


if(-e "$notthese" && -e "$fullref" && ! -e "$pdbfile" ) then 
# look for orphan reference points
  echo "extracting goal points from $fullref that are not in $notthese"
  cat $notthese $fullref |\
   awk '! /^ATOM|^HETAT/{next} {typ=" "substr($0,18,3)" "}\
    typ !~ /HOH|EDO|NH4| CL | LI /{next}\
    {xyz=substr($0,31,24);++seen[xyz]}\
   END{for(xyz in seen)if(seen[xyz]==1){\
    printf("ATOM      1  O   HOH z%8d%s  1.00  0.00           O         \n",++n,xyz)}}' |\
  hy36_encode.awk |\
  cat >! ${t}input_refpoints.pdb
endif

if(-e "$pdbfile" && ! -e "$notthese") then
  echo "extracting goal points from $pdbfile"
  cat $pdbfile |\
   awk '! /^ATOM|^HETAT/{next}\
    {xyz=substr($0,31,24);\
    printf("ATOM      1  O   HOH z%8d%s  1.00  0.00           O         \n",++n,xyz)}' |\
  hy36_encode.awk |\
  cat >! ${t}input_refpoints.pdb
endif

if(-e "$pdbfile" && -e "$notthese") then
  echo "extracting goal points from $pdbfile that are not in $notthese"
  awk '/^ATOM|^HETAT/{print substr($0,31,24),"NOTME"}' $notthese |\
  cat - $pdbfile |\
   awk '$NF=="NOTME"{++notme[substr($0,1,24)];next}\
    ! /^ATOM|^HETAT/{next}\
    {xyz=substr($0,31,24)}\
    notme[xyz]{next}\
    {printf("ATOM      1  O   HOH z%8d%s  1.00  0.00           O         \n",++n,xyz)}' |\
  hy36_encode.awk |\
  cat >! ${t}input_refpoints.pdb
endif

set num = `cat ${t}input_refpoints.pdb | wc -l`
echo "$num points are selected"


echo "reading xyz from $rstfile"
rst2pdb_runme.com parmfile=$topfile $rstfile outfile=${t}cpptraj.pdb >! ${t}cpptraj.log
echo "encoding"
egrep "^ATOM|^HETAT" ${t}cpptraj.pdb |\
awk '{pre=substr($0,1,21);post=substr($0,23+8);resnum=substr($0,23,8);}\
     resnum!=prevres{++n;prevres=resnum}\
    {printf("%sA%8d%s\n",pre,n,post)}' |\
cat >! ${t}w8.pdb
echo "ecnoding with hy36"
cat ${t}w8.pdb |\
hy36_encode.awk |\
cat >! ${t}hy36.pdb

echo "looking for potential new clashes with selected goal points"
egrep "^CRYST" $fullref >! ${t}gemmime.pdb
if( $status ) then
   egrep "^CRYST" $pdbfile >! ${t}gemmime.pdb
endif
if( $status ) then
   egrep "^CRYST" ${t}cpptraj.pdb >! ${t}gemmime.pdb
endif
egrep "^CRYST" ${t}gemmime.pdb > /dev/null
if( $status ) then
   set BAD = "cannot find CRYST1"
   goto exit
endif
cat ${t}hy36.pdb ${t}input_refpoints.pdb |\
cat  >> ${t}gemmime.pdb
gemmi contact -d 4 --sort --ignore=3 ${t}gemmime.pdb >! ${t}gemmiclash.log
grep "HOH z" ${t}gemmiclash.log |\
awk -v debug=$debug 'debug>1{print $0,"DEBUG IN"}\
  {seg2=substr($0,31);\
   resid=substr(seg2,index(seg2,"HOH z")+5,8);\
   print resid "|  " $NF+0,"DIST";}' |\
awk '! seen[$1]{print;++seen[$1]}' |\
cat - ${t}input_refpoints.pdb |\
awk '$NF=="DIST"{resid=$1;dist[resid]=$3;next}\
 ! /^ATOM|^HETAT/{next}\
   {resid=substr($0,23,6)+0;resid=$NF}\
  dist[resid]==""{dist[resid]=999}\
  {print substr($0,1,80)," ",dist[resid]}' |\
sort -k1.81gr >! ${t}open_refpoints.pdb

echo "probing fofc map with goal points"
phenix.map_value_at_point $mtzfile ${t}open_refpoints.pdb \
          ${phenixlabel}=$mtzlabel scale=sigma >! ${t}map_values_ref.log

echo "looking for locations with rho>$minrho and >$mindist A to nearest neighbor"
cat ${t}map_values_ref.log ${t}open_refpoints.pdb |\
awk '/Map value:/{++n;rho[n]=$NF;next}\
  /^CRYST/{print;next}\
  ! /^ATOM|^HETAT/{next}\
  {++m}\
  {print $0,"  ",rho[m]}' >! ${t}rholabel_ref.pdb

echo $minrho $mindist |\
cat - ${t}rholabel_ref.pdb |\
awk 'NR==1{minrho=$1;mindist=$2;next}\
 $NF>minrho && $(NF-1)>mindist{print $NF*$(NF-1),"|"  $0}' |\
sort -gr |\
awk '{print substr($0,index($0,"|")+1)}' |\
cat >! ${t}vacant_refpoints.pdb 

set vacancies = `cat ${t}vacant_refpoints.pdb | wc -l`
if( $vacancies == 0 ) then
  echo "there are none. Biggest void:"
  sort -k1.81gr ${t}rholabel_ref.pdb | head -n 1
  goto exit
endif

echo "looking for post-teleport clashes"
egrep "^CRYST" $fullref >! ${t}candidate_refpoints.pdb
egrep "^ATOM|^HETAT" ${t}vacant_refpoints.pdb >> ${t}candidate_refpoints.pdb
 
set clashes = 1
while ( $clashes ) 

gemmi contact -d $mindist ${t}candidate_refpoints.pdb >! ${t}gemmi.log
grep "HOH z" ${t}gemmi.log |\
awk -v debug=$debug 'debug>1{print $0,"DEBUG IN"}\
  {resid1=substr($0,index($0,"HOH z")+5,8);\
   seg2=substr($0,31);\
   resid2=substr(seg2,index(seg2,"HOH z")+5,8);}\
  ! deleted[resid1]{++deleted[resid2]}\
  deleted[resid1]{next}\
  {print resid2 "|" resid1 "|  " $NF+0,"CLASH";++printed[resid1]}' |\
cat >! ${t}clashes.txt 

# decode resids
hy36_encode.awk ${t}candidate_refpoints.pdb |\
cat ${t}clashes.txt - |\
awk '$NF=="CLASH"{++bad[$1];next}\
  ! /^ATOM|^HETAT/{print;next}\
  {resid=$NF}\
  bad[resid]{next}\
  {print}' |\
cat >! ${t}declashed.pdb

set clashes = `cat ${t}clashes.txt | wc -l`
set num = `egrep "^ATOM|^HETAT" ${t}declashed.pdb | wc -l`
if( $num == 0 ) set clashes = 0
echo "$clashes clashes leaves $num targets"

mv ${t}declashed.pdb ${t}candidate_refpoints.pdb

end

set candidates = `egrep "^ATOM|^HETAT" ${t}candidate_refpoints.pdb | wc -l`
egrep "^ATOM|^HETAT" ${t}candidate_refpoints.pdb |\
head -n $maxmoves |\
awk '{printf("ATOM  %5d%s\n",++n%10000,substr($0,12))}' |\
cat >! ${t}selected_refpoints.pdb

set num = `cat ${t}selected_refpoints.pdb | wc -l`
if( $maxmoves > $num ) set maxmoves = $num

echo "                                                                                   dist    Fo-Fc"
head -n 1 ${t}selected_refpoints.pdb
if( $maxmoves > 1 ) head -n $maxmoves ${t}selected_refpoints.pdb | tail -n 1

echo "taking top $maxmoves of $candidates points"
head -n $maxmoves ${t}selected_refpoints.pdb |\
awk '{print substr($0,31,8),substr($0,39,8),substr($0,47,8),"NEWXYZ"}' |\
cat >! ${t}newxyz.txt


if(-e "$notthese") then
   echo "excluding atoms already restrained in $notthese"
   cat $notthese |\
   awk '/^ATOM|^HETAT/{print $0,"NOTME"}' |\
   cat - ${t}cpptraj.pdb |\
   awk '! /^ATOM|^HETAT/{next}\
    {id=substr($0,12,17)}\
    $NF=="NOTME"{++skip[id];next}\
    skip[id]{next}\
    /O   HOH/{print}' >! ${t}xyz.pdb
#   set natom = `egrep "^ATOM|^HETAT" ${t}cpptraj.pdb | wc -l`
#   set unrest = `egrep "^ATOM|^HETAT" ${t}unrestrained.pdb | wc -l`
#   echo "were $natom now $unrest left."
else
    cat ${t}cpptraj.pdb |\
    egrep "O   WAT" |\
    convert_pdb.awk -v skip=EP,H >! ${t}xyz.pdb
endif

if(-e "$teleportee") then
  echo "using teleportation candidates listed in $teleportee"
  cat $teleportee |\
  awk '! /^ATOM|^HETAT/{next}\
    {++m;rho=$NF;oresnum=substr($0,80)+0;\
      X=substr($0,31,8)+0;\
      Y=substr($0,39,8)+0;\
      Z=substr($0,47,8)+0;\
     print rho,oresnum,X,Y,Z,"OLDXYZ"}' |\
  head -n $maxmoves >! ${t}rho_resnum_oldxyz.txt

  set test = `cat ${t}rho_resnum_oldxyz.txt | wc -l`
  if( $test ) goto skipold
  echo "WARNING: no useful atoms in $teleportee"
endif

echo "probing fofc map at current water positions"
phenix.map_value_at_point $mtzfile ${t}xyz.pdb \
          ${phenixlabel}=$mtzlabel scale=sigma >! ${t}map_values.log

cat ${t}map_values.log ${t}xyz.pdb |\
awk '/Map value:/{rho[++n]=$NF;next}\
  ! /^ATOM|^HETAT/{next}\
  {++m;X=substr($0,31,8)+0;\
       Y=substr($0,39,8)+0;\
       Z=substr($0,47,8)+0;}\
  {print rho[m],$NF,X,Y,Z,"OLDXYZ"}' |\
sort -g |\
head -n $maxmoves >! ${t}rho_resnum_oldxyz.txt

skipold:
head -n 1 ${t}rho_resnum_oldxyz.txt
if( $maxmoves > 1 ) tail -n 1 ${t}rho_resnum_oldxyz.txt
# rho resnum X Y Z OLDXYZ


# find high slots for new restraints
cat ${t}rho_resnum_oldxyz.txt ${t}xyz.pdb |\
awk '$NF=="OLDXYZ"{++taken[$2];next}\
! /^ATOM|^HETAT/{next}\
  taken[$NF]{next}\
  {++n;X=substr($0,31,8)+0;\
       Y=substr($0,39,8)+0;\
       Z=substr($0,47,8)+0;}\
  {print $NF,$NF,X,Y,Z,"NEWSLT"}' |\
head -n $maxmoves >! ${t}slots.txt

set new = `cat ${t}newxyz.txt | wc -l`
set old = `cat ${t}rho_resnum_oldxyz.txt | wc -l`
set num = `echo $new $old | awk '$1>$2{$1=$2} {print $1}'`
echo "teleporting $num waters"
cat ${t}newxyz.txt ${t}slots.txt ${t}rho_resnum_oldxyz.txt ${t}cpptraj.pdb |\
awk '$NF=="NEWXYZ"{++n;newX[n]=$1;newY[n]=$2;newZ[n]=$3;next}\
     $NF=="NEWSLT"{++s;oresnum=$2;\
        if(newX[s]=="")next;\
        ++moving[oresnum];++dest[oresnum]\
        stayX[s]=$3;stayY[s]=$4;stayZ[s]=$5;\
        dX[oresnum]=newX[s]-stayX[s];dY[oresnum]=newY[s]-stayY[s];dZ[oresnum]=newZ[s]-stayZ[s];\
        next}\
     $NF=="OLDXYZ"{++o;oresnum=$2;if(o>n)next;\
        if(newX[o]=="")next;\
        ++moving[oresnum];\
        oldX[o]=$3;oldY[o]=$4;oldZ[o]=$5;\
        dX[oresnum]=stayX[o]-oldX[o];dY[oresnum]=stayY[o]-oldY[o];dZ[oresnum]=stayZ[o]-oldZ[o];\
        next}\
  ! /^ATOM|^HETAT/{print;next}\
  {oresnum=$NF}\
  ! moving[oresnum]{print;next}\
  {tag=""} dest[oresnum]{tag="moved"}\
  {pre=substr($0,1,30);post=substr($0,55);\
   X=substr($0,31,8)+0;Y=substr($0,39,8)+0;Z=substr($0,47,8)+0;\
   X+=dX[oresnum];Y+=dY[oresnum];Z+=dZ[oresnum];\
   printf("%s%8.3f%8.3f%8.3f%s  %s\n",pre,X,Y,Z,post,tag)}' |\
cat >! ${t}moved.pdb

grep moved ${t}moved.pdb | tee ${t}sanity.pdb | grep " O " | tee ${t}movedwater.pdb
#rmsd $fullref ${t}sanity.pdb | head
awk '/moved/{print $(NF-1)}' ${t}sanity.pdb >! teleported_residues.txt

echo "grafting into current system"
rm -f $outfile > /dev/null
graft_atoms_runme.com ${t}moved.pdb $rstfile \
  minimize=$minimize energycheck=$energycheck debug=$debug \
  outrst=$outfile >! ${t}graft.log
if( $status || ! -e "$outfile") then
  set BAD = "graft failed"
  goto exit
endif

if( $minimize ) then
  awk '/^minimizing/,/^appending/' ${t}graft.log
endif
if( $energycheck ) then
  awk '/^energy/,/^now/' ${t}graft.log
endif

# now remap those new reference points?
if(-e "$notthese") then
   echo ""
   combine_pdbs_runme.com $notthese ${t}movedwater.pdb ${t}cpptraj.pdb \
    outfile=${t}newrefpoints.pdb >! ${t}final_combine.log
    
endif


goto exit


exit:

if($debug) wc $rstfile $outfile

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif

exit


awk 'NR>2{for(i=0;i<5;++i){f=substr($0,i*12+1,12);if(length(f==12))print ++n/3,f}}' teleported.rst7 |\
cat  >! teleported.txt

rst2pdb_runme.com teleported.rst7 >& /dev/null ; rst2pdb_runme.com wrapped.rst7 > /dev/null
rmsd teleported.pdb wrapped.pdb

egrep -v "EPW|EG7" teleported.pdb >! testme.pdb
phenix.geometry_minimization macro_cycles=0 testme.pdb | tee phenixgeo.log



