#! /bin/tcsh -f
#
#  find biggest outlier restraints in a B-weighted PDB file
#  release them, their symmetry mates, and neighbors
#
#
set smallSG = ""
set mtzfile = reference.mtz

set tootight  = 9
set resetweight = 0.1
set sigma = 5
set radius = 2.2
set radii = 2.2,1.8,1.8,1.8
set keepcenter = 0
set maxbad    = 10
set pdbscale  = 0.01

set debug = 0

set restraintpdb = current_restraints.pdb
set outfile      = new_restraints.pdb

set tempfile = /dev/shm/${USER}/tempfile_rwr_$$_

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
      if("$key" == "restpdb") set restraintpdb = "$Val"
      if("$key" == "refpoints") set restraintpdb = "$Val"

      if("$key" == "multsfile") set multsfile = "$Val"
      if("$key" == "output") set outfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set restraintpdb = $Arg
    endif
    if("$key" == "debug") set debug = "$Val"
end

if(! -e "$restraintpdb") then
    set BAD = "no coordinates provided"
    goto exit
endif

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_rwr_
set tmpdir = `dirname $tempfile`

set t = "$tempfile"

# gather stats on restraints
echo "gathering restraint challenge history"
sort -k1.61gr $restraintpdb |\
 egrep  "^ATOM|^HETAT" |\
tee sorted_weights.pdb |\
awk -v pdbscale=$pdbscale '{print substr($0,61,6)*pdbscale}' |\
cat >! sorted_weights.txt

# these are on amber scale
set challenges = `ls -1rt restraints_challenge_vs_weight_*.txt | tail -n 200`
awk '{id=substr($0,length($0)-16);itr=substr(FILENAME,32);E=$1;d[id]=$2;w[id]=$3}\
  {sum[id]+=E*itr} END{for(id in sum)print sum[id],d[id],w[id],"|"id}' $challenges |\
 sort -g |\
 cat >! sorted_iwr.txt 
tail -n 3 sorted_iwr.txt 
set badB = `echo $tootight $pdbscale | awk '{print $1/$2}'`
echo "finding challenged restraints with weight > $tootight"
egrep "^CRYST" $restraintpdb >! worst.pdb
cat sorted_iwr.txt |\
cat - $restraintpdb |\
awk -F "|" -v tootight=$tootight -v pdbscale=$pdbscale '\
 NF==2{id=substr($2,1,16);score[id]=$1+0;next}\
 ! /^ATOM|^HETAT/{next}\
 /NOT_A_BOMB| moved/{next}\
 {id=substr($0,12,16);w=substr($0,61,6)*pdbscale}\
 score[id] && w>tootight{print score[id],"|"$0}' |\
sort -gr |\
awk -F "|" '{print $2}' |\
cat >! sorted_worstchallenge.pdb
head -n $maxbad sorted_worstchallenge.pdb |\
tee -a worst.pdb

goto gotworst



gotworst:

set test = `egrep "^ATOM|^HETAT" worst.pdb | wc -l`

echo "$test restraints are too tight"

if( $test == 0 ) then
    cp -p $restraintpdb $outfile
    goto exit
endif

if(-e "$mtzfile") then
    set smallSGnum = `echo head | mtzdump hklin $mtzfile | awk '/Space group =/{print $NF+0}' | tail -n 1`
    set smallSG = `awk -v num=$smallSGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
    echo "got $smallSG from $mtzfile"
endif



echo "selecting $smallSG symmetry mates of worst restraints"
neighbors_big.com $smallSG $restraintpdb worst.pdb 1A -include -outfile ${t}symmates.pdb

if( "$radii" == "" ) set radii = "$radius"
set radii = `echo $radii | awk '{gsub(","," ");print}'`
set maxradius = `echo $radii | awk '{for(i=1;i<=NF;++i)if($i>max)max=$i;print max}'`

distanceify_runme.com P1 ${maxradius}A $restraintpdb


awk '/^ATOM|^HETAT/{print substr($0,12,17) "| SELECTED"}' worst.pdb ${t}symmates.pdb >! ${t}core.txt
cp ${t}core.txt ${t}allselected.txt

echo "extention radii: $radii "
foreach rad ( $radii )
  echo "extending selection by ${rad}A"
  cat ${t}core.txt distances.txt |\
  awk -v rad=$rad -F "|" '/SELECTED/{++selected[$1];next}\
    $3+0>rad{next}\
    selected[$1]{print $2 "| SELECTED"}\
    selected[$2]{print $1 "| SELECTED"}' |\
  cat >! ${t}new.txt
  set new = `cat ${t}new.txt | wc -l`
  echo "$new new selections"
  cat ${t}new.txt >> ${t}allselected.txt
  sort -u ${t}new.txt >! ${t}core.txt
end
sort -u ${t}allselected.txt >! ${t}new.txt
set new = `cat ${t}new.txt | wc -l`
echo "$new unique selections"
mv ${t}new.txt ${t}allselected.txt 


sort -u ${t}allselected.txt |\
 cat - $restraintpdb |\
 awk -F "|" '/SELECTED/{++sel[$1];next}\
   ! /^ATOM|^HETAT/{print;next}\
   {id=substr($0,12,17)} sel[id]{print}' |\
cat >! selected.pdb
echo $resetweight $tootight $pdbscale  |\
cat - selected.pdb |\
awk 'NR==1{maxB=$1;tootight=$2/$3}\
  ! /^ATOM|^HETAT/{print;next}\
   {B=substr($0,61,6)+0;pre=substr($0,1,60);post=substr($0,67)}\
   #B<tootight{next}\
   B>maxB{B=maxB}\
   {printf("%s%6.2f%s\n",pre,B,post)}' |\
cat >! releaseme.pdb

combine_pdbs_runme.com printref=1 releaseme.pdb $restraintpdb outfile=$outfile | head -n 1
rmsd -v debug=1 $restraintpdb $outfile | grep -v "  0.00 (B)" | sort -k1.53g


exit:
ls -l $outfile

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
#   echo "GOTEHRE: cleaning up ${t}"
#   ls -l ${t}*
   rm -f ${t}* >& /dev/null
endif

exit




awk '{id=substr($0,12,16);w=substr($0,61,6);itr=substr(FILENAME,16)+0;sum[id]+=itr*w} END{for(x in sum)print sum[x]"|"x}' restraints_for_???.pdb | sort -g | tee sorted_cumw.txt


awk '{id=substr($0,length($0)-16);itr=substr(FILENAME,32);rw=$1;sum[id]+=rw*itr} END{for(id in sum)print sum[id],"|"id}' restraints_challenge_vs_weight_???.txt | sort -g | tee sorted_iwr.txt






awk '/^CRYST/ || substr($0,61,6)+0>0' current_restraints.pdb >! nonzero.pdb

gemmi contact --sort -d 2.8 current_restraints.pdb |\
tee ${t}gemmi.log |\
awk -v debug=$debug '{dist=$NF}\
   {if(debug)print $0,"DEBUG";\
   id1=substr($0,12,17);\
   id2=substr($0,12+30,17);\
      atom1=substr(id1,1,5);atom2=substr(id2,1,5);\
      typ1=substr(id1,7,3);typ2=substr(id2,7,3);\
      chain1=substr(id1,11,1);chain2=substr(id2,11,1);\
      resnum1=substr(id1,12);resnum2=substr(id2,12);\
      gsub(" ","",atom1);gsub(" ","",atom2)}\
   resnum2+0==resnum1+1 && chain1==chain2 && atom1=="C" && atom2=="N"{next}\
   atom2=="ZN" && atom1~/^[ON]/ && typ1~/KCX|HIS/{next}\
   {print id1 "|" id2 "|  " dist,"DIST";}' |\
tee ${t}contact_distances.txt |\
cat - current_restraints.pdb |\
awk -F "|" '/DIST/{dist[$1,$2]=$NF+0}\
 /DIST/ && ! mindist[$1]{mindist[$1]=$NF+0;mate[$1]=$2}\
 /DIST/ && ! mindist[$2]{mindist[$2]=$NF+0;mate[$2]=$1}\
 /DIST/{next}\
 ! /^ATOM|^HETAT/{next}\
 {id=substr($0,12,17);w[id]=substr($0,61,6)+0;}\
 END{for(id in w){\
   if(! mindist[id])continue;\
   if(w[id] <= 0)continue;\
   if(w[mate[id]] <= 0)continue;\
   print (w[id]+w[mate[id]])/mindist[id]**2,mindist[id],"|"id"|"mate[id]"|",w[id],w[mate[id]]}}\
 ' |\
sort -g |\
tee contact_vs_id.txt



  echo "looking for too-close water refpoints"
  convert_pdb.awk -v skip=protein current_restraints.pdb |\
    egrep -v "EG7| ZN " >! solvent_restraints.pdb 
  distanceify_runme.com current_restraints.pdb solvent_restraints.pdb \
     maxdist=2.4 outfile=distances.txt >! dist.log
  egrep "HOH|EDO|SO4| LI | CL |NH4|ACY" restraints_challenge_vs_weight.txt |\
  awk '{print substr($0,length($0)-16) "|",$1,"|",$2,"|",$3,"CHALLENGE"}' |\
  cat - distances.txt |\
  awk -v d2=1.8 -F "|" '/CHALLENGE/{challenge[$1]=$2+0;d[$1]=$3+0;w[$1]=$4+0;next}\
     /DIST/{dist[$1,$2]=$3+0;dist[$2,$1]=dist[$1,$2]}\
     /DIST/ && ( $1==$2 || dist[$1,$2]>2.4 ){;next}\
     elim[$1] || elim[$2]{next}\
     {sol1=sol2=polya1=polya2=0;typ1=substr($1,7,3);typ2=substr($2,7,3)}\
     typ1~/HOH|EDO|SO4| LI | CL |NH4|ACY/{++sol1}\
     typ2~/HOH|EDO|SO4| LI | CL |NH4|ACY/{++sol2}\
     typ1~/EDO|SO4|ACY/{++polya1}\
     typ2~/EDO|SO4|ACY/{++polya2}\
     # do not eliminate like-to-like polyatom restraints \
     polya1 && typ1==typ2{next}\
     ! sol1 && sol2 && dist[$1,$2]<d2{print $2 "| ELIM";++elim[$2]}\
     sol1 && ! sol2 && dist[$1,$2]<d2{print $1 "| ELIM";++elim[$1]}\
     challenge[$1]>challenge[$2] && w[$1]>700 && sol1{print $1 "| ELIM";++elim[$1]}\
     challenge[$2]>challenge[$1] && w[$2]>700 && sol2{print $2 "| ELIM";++elim[$2]}' |\
  cat - current_restraints.pdb |\
  awk '$NF=="ELIM"{++elim[substr($0,1,16)];next}\
    ! /^ATOM|^HETAT/{print;next}\
     {id=substr($0,12,16)}\
    elim[id]{next}\
     {print}' >! filtered_refpoints.pdb





append_file_date.com all_possible_refpoints.pdb
tac restraints_challenge_vs_weight.txt |\
awk '$3>9 && $2>1. || NR==1{\
  l=length($0);print substr($0,l-16)"   BAD"}' |\
cat - all_possible_refpoints.pdb |\
awk '$NF=="BAD"{id=substr($0,1,15);++bad[id];next}\
  {id=substr($0,12,15)}\
  ! bad[id]{print}' |\
cat >! culled.pdb
diff all_possible_refpoints.pdb culled.pdb
mv culled.pdb all_possible_refpoints.pdb



