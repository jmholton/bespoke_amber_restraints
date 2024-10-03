#! /bin/tcsh -f
#
#  make new restraints pdb file based on atoms near centroids             -James Holton 4-19-24
#
#  needs convert_pdb.awk scipt and ccp4
#
#
set reffile    = all_possible_refpoints.pdb
set pdbfile    = amberme.pdb
set outfile    = restraints_to_use.pdb

# max/default restraint weight
set weight = 1.0
# maximum distance to still be restrained
set maxdist = 1.0
# reference points that are within nonbond distance
set tooclose = 2.4
# scaling factor for to-water distances
set hohscale = 2
# softening factor for atoms far from centroid points.
# it is in order of magnitude: 2 drops by factor of 100 at max dist, 0 means no softening
set softener = 2
# make sure ambivalent pairs of atoms are both included
set twirl = 1
# represent all reference points as waters to gemmi
set gemmibug = 0
# flag to debug things
set debug = 0

set path = ( `dirname $0` $path )

set tempfile = /dev/shm/${USER}/tempfile_$$_
#set tempfile = ./tempfile_c2r_


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
      if("$key" == "reference") set reffile = "$Val"
      if("$key" == "ref") set reffile = "$Val"
      if("$key" == "output") set outfile = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
      if("$key" == "weight") set weight = "$Val"
      if("$key" == "maxdist") set maxdist = "$Val"
      if("$key" == "softener") set softener = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $debug && "$tempfile" =~ /dev/shm/* ) set tempfile = ./tempfile_c2r_
set tmpdir = `dirname $tempfile`
mkdir -p $tmpdir >& /dev/null

cat << EOF
pdbfile  $pdbfile
reffile  $reffile
outfile  $outfile

weight   $weight
maxdist  $maxdist
hohscale $hohscale
softener $softener

tempfile  $tempfile
EOF

if(! -r "$pdbfile") then
  set BAD = "cannot read $pdbfile"
  goto exit
endif
if(! -r "$reffile") then
  set BAD = "cannot read $reffile"
  goto exit
endif

set t = "$tempfile"

# label everything with ordinal resnums first?


echo "renaming relevant MSEs to METs ..."
awk '/SD  MET/{print substr($0,18,13),"METME"}' $pdbfile |\
cat - $reffile |\
awk '$NF=="METME"{id=substr($0,5,8);++met[id];next}\
  ! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,22,8)}\
  met[id]{gsub("SE   MSE"," SD  MET");gsub("MSE","MET");gsub("SE"," S")}\
  {print}' >! ${t}centroids_met.pdb

set changes = `diff ${t}centroids_met.pdb $reffile | grep MET | wc -l`
echo "$changes changes"

echo "encoding atom names"
egrep "^CRYST" $reffile        | head -n 1 >! ${t}temp.pdb

egrep "^ATOM|^HETAT" $pdbfile |\
  awk 'BEGIN{ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";\
       ci=1;c=substr(ABC,ci,1)}\
      substr($0,77,2)~/ H|XP/{next}\
      {origid=substr($0,12,19);xyz=substr($0,31,24);pre=substr($0,1,12);\
       resid=substr($0,22,9);\
       gsub(" ","_");\
       atom=substr($0,14,3);f=substr($0,17,1);typ=substr($0,18,4)}\
      resid!=lastid{r+=1;lastid=resid}\
       r>36**4{++ci;c=substr(ABC,ci,1);r=1}\
      debug>1{print "DEBUG" substr($0,6)}\
      {printf("%s%s%s_%s%s%8d%s  1.00  0.00           O            |%s| MODEL1 |\n",pre,f,atom,typ,c,r,xyz,origid)}' |\
  cat >> ${t}temp.pdb

egrep "^ATOM|^HETAT" ${t}centroids_met.pdb |\
awk 'BEGIN{ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";\
       ci=length(ABC);c=substr(ABC,ci,1)}\
      {origid=substr($0,12,19);xyz=substr($0,31,24);pre=substr($0,1,12);\
       resid=substr($0,22,9);\
       gsub(" ","_");\
       atom=substr($0,14,3);f=substr($0,17,1);typ=substr($0,18,4)}\
      resid!=lastid{r+=1;lastid=resid}\
       r>36**4{--ci;c=substr(ABC,ci,1);r=1}\
      debug>1{print "DEBUG" substr($0,6)}\
      {printf("%s%s%s_%s%s%8d%s  1.00  0.00           O            |%s| MODEL2 |\n",pre,f,atom,typ,c,r,xyz,origid)}' |\
cat >> ${t}temp.pdb


# final pass to check for duplicates
#convert_pdb.awk -v renumber=hy36 -v append=encresnum ${t}temp.pdb |\
#awk '{print $0, "   "}' >! ${t}compme.pdb
# final pass ... encode
hy36_encode.awk ${t}temp.pdb |\
awk '{print $0, "   "}' >! ${t}compme.pdb
#convert_pdb.awk -v renumber=w4,ordinal,chain ${t}temp.pdb |\
#cp ${t}temp.pdb ${t}compme.pdb

set errors = `grep ERROR ${t}compme.pdb | wc -l`
if( $errors ) then
    set BAD = "$errors errrors converting pdb files with convert_pdb.awk"
    goto exit
endif

echo "measuring distances..."
gemmi contact -d $maxdist --sort --nosym --ignore=3 ${t}compme.pdb >! ${t}contact.log
# sort -k1.73g

echo "decoding atom names..."
cat ${t}contact.log |\
awk -v debug=$debug '{dist=$NF}\
   {if(debug>1)print $0,"DEBUG IN";\
   id1=substr($0,12,19);\
   id2=substr($0,12+30,23);\
   while(gsub(/^ /,"",id2));\
   while(gsub(/ $/,"",id2));\
   id2=sprintf(" %-18s",id2);\
   print id1 "|" id2 "|  " dist,"DIST";}' |\
cat >! ${t}encoded_distances.txt

# change back to original names
egrep -vh DEBUG ${t}compme.pdb ${t}encoded_distances.txt |\
awk -v debug=$debug 'debug>3{print $0,"DEBUG IN"}\
   {split($0,w,"|")}\
   {m=1} / MODEL2 /{m=2}\
   /^ATOM|^HETAT/{id=substr($0,12,19);\
     # reconstruct gemmi-decoded hy36 \
     ordresnum=sprintf("%-8s",sprintf("%4d",$NF));\
     id=substr($0,12,11) ordresnum;\
     origid[id,m]=substr(w[2],1,19);\
     if(debug>1)print "DEBUG |"id"|"origid[id,m]"|";}\
   $NF!="DIST"{next}\
     {id1=w[1];id2=w[2];\
      dist=$(NF-1);\
      if(debug>1)print "DEBUG |"id1"|"id2"|---|"origid[id1,1]"|"origid[id2,2]"|"}\
    origid[id1,1]!="" && origid[id2,2]!=""{\
      print origid[id1,1] "|" origid[id2,2] "|",dist,"DIST"}' |\
cat >! ${t}decoded_distances.txt

# make water distances longer so protein gets a better chance
echo "scaling up to-water distances by $hohscale and cutting at $maxdist"
egrep -vh DEBUG  ${t}decoded_distances.txt |\
awk -v s=$hohscale -F "|" '/HOH/{id1=$1;id2=$2;dist=$3*s;\
   print id1 "|" id2 "|",dist,"DIST";next}\
  {print}' |\
awk -v d=$maxdist -F "|" '$3+0<d' |\
sort -t "|" -k3g |\
cat >! ${t}scaled_distances.txt

# also get intra-reference-point distances

echo "measuring intra-reference distances"
egrep "^CRYST1| MODEL2" ${t}compme.pdb >! ${t}refpoints.pdb
gemmi contact -d $tooclose --sort ${t}refpoints.pdb >! ${t}intracontact.log


echo "prioritizing intra-refpoint distances"
echo "decoding atom names..."
cat ${t}intracontact.log |\
awk -v debug=$debug '{dist=$NF}\
   {if(debug>1)print $0,"DEBUG IN";\
   id1=substr($0,12,19);\
   n1=sprintf("%4d",substr(id1,12,10)+0);\
   id2=substr($0,12+30+length(n1)-4,19);\
   print id1 "|" id2 "|  " dist,"DIST";}' |\
cat >! ${t}encoded_intradistances.txt

egrep -vh DEBUG ${t}refpoints.pdb ${t}encoded_intradistances.txt |\
awk -v debug=$debug 'debug>3{print $0,"DEBUG IN"}\
   {split($0,w,"|")}\
   /^ATOM|^HETAT/{id=substr($0,12,20);\
     oresnum=sprintf("%-8s",sprintf("%4d",w[4]));\
     id=substr($0,12,11) oresnum;\
     origid[id]=substr(w[2],1,19);\
     if(debug>1)print "DEBUG |"id"|"origid[id]"|";}\
   $NF!="DIST"{next}\
     {id1=w[1];id2=w[2];\
      dist=$(NF-1);\
      if(debug>3)print "DEBUG |"origid[id1]"|"origid[id2]"|";\
      if(debug>1)print "DEBUG |"id1"|"id2"|---|"origid[id1]"|"origid[id2]"|"}\
    origid[id1]!="" && origid[id2]!=""{\
      print origid[id1] "|" origid[id2] "|",dist,"DIST"}' |\
cat >! ${t}intraref_distances.txt


# now determine which atoms to restrain to what:
# protein must have same name, or solvent can also be protein restraints
echo "selecting protein..."
egrep -v DEBUG ${t}scaled_distances.txt |\
awk -v debug=$debug '$NF=="DIST"{\
      split($0,w,"|");\
      id1=w[1];id2=w[2];\
      atom1=substr(id1,1,5);atom2=substr(id2,1,5);\
      typ1=substr(id1,7,3);typ2=substr(id2,7,3);\
      resid1=substr(id1,11);resid2=substr(id2,11);\
      dist=$(NF-1);}\
  debug>1{print $0,"DEBUG"}\
  resid1==resid2 && typ1==typ2 && typ1=="ASP" && atom1~/^  OD/ && atom2~/^  OD/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="GLU" && atom1~/^  OE/ && atom2~/^  OE/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1~/PHE|TYR/ && atom1~/^  CD/ && atom2~/^  CD/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1~/PHE|TYR/ && atom1~/^  CE/ && atom2~/^  CE/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="ASN" && atom1~/^  [NO]D/ && atom2~/^  [NO]D/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="GLN" && atom1~/^  [NO]E/ && atom2~/^  [NO]E/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="HIS" && atom1~/^  [NC]D/ && atom2~/^  [NC]D/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="HIS" && atom1~/^  [NC]E/ && atom2~/^  [NC]E/{atom2=atom1}\
      {ent1=atom1" "typ1" "resid1;ent2=atom2" "typ2" "resid2;\
       atm1=atom1;atm2=atom2;gsub(" ","",atm1);gsub(" ","",atm2)}\
  seen1[ent1]{if(debug>1)print "seen1",ent1,"DEBUG";next} \
  seen2[ent2]{if(debug>1)print "seen2",ent2,"DEBUG";next} \
  {protein1=protein2=solvent1=solvent2=solvent=metal1=metal=polya=0}\
  typ1~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein1}\
  typ1~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein1}\
  typ2~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein2}\
  typ2~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein2}\
  typ1~/HOH|CL |AMM|NH4|LI /{++solvent1}\
  typ2~/HOH|CL |AMM|NH4|LI |EDO|SO4/{++solvent2}\
  typ1~/EDO|SO4|ACY/ && typ1==typ2 && atm1==atm2{++polya}\
  typ1~/ZN |EG7/{++metal1}\
  metal1 && typ1==typ2 && atm1==atm2{++metal}\
    ent1==ent2 || polya || metal || dist<2 && solvent2 && ! metal1 {\
    if(debug>1){print "DEBUG",(ent1==ent2),(polya),(metal),(dist<2 && solvent2)}\
    print;\
    ++seen1[ent1];++seen2[ent2]}' |\
cat >! ${t}protein_selections.txt


# add other member of ambivalent pairs
egrep -v DEBUG ${t}protein_selections.txt |\
awk -v debug=$debug 'BEGIN{\
  other["OD1 ASP"]="OD2";other["OD2 ASP"]="OD1";\
  other["OD1 ASN"]="ND2";other["ND2 ASN"]="OD1";\
  other["OE1 GLU"]="OE2";other["OE2 GLU"]="OE1";\
  other["OE1 GLN"]="NE2";other["NE2 GLN"]="OE1";\
  other["CD1 PHE"]="CD2";other["CD2 PHE"]="CD1";\
  other["CE1 PHE"]="CE2";other["CE2 PHE"]="CE1";\
  other["CD1 TYR"]="CD2";other["CD2 TYR"]="CD1";\
  other["CE1 TYR"]="CE2";other["CE2 TYR"]="CE1";\
  other["CE1 HIS"]="NE2";other["NE2 HIS"]="CE1";\
  other["CD2 HIS"]="ND1";other["ND1 HIS"]="CD2";\
  rigid["SO4"]="_O1 _O2 _O3 _O4 _S";\
  rigid["EDO"]="_C1 _O1 _C2 _O2";\
}\
debug>1{print $0,"DEBUG"}\
  {split($0,w,"|");id1=w[1];id2=w[2];\
      atom1=substr(id1,1,5);atom2=substr(id2,1,5);\
      typ1=substr(id1,7,3);typ2=substr(id2,7,3);\
      resid1=substr(id1,11);resid2=substr(id2,11);\
      atm1=atom1;atm2=atom2;gsub(" ","",atm1);gsub(" ","",atm2)}\
 other[atm1" "typ1]!=""{\
  if(debug>1){print "OTHER",atm1,typ1,"=",other[atm1" "typ1],"DEBUG"};\
  post1=substr(id1,6);post2=substr(id2,6);;dist=w[3];\
  oa1=other[atm1" "typ1];oa2=other[atm2" "typ2];\
  newid1="  "oa1 post1;newid2="  "oa2 post2;\
  print newid1 "|" newid2 "|" dist;\
 }\
 rigid[typ2]==""{print;next}\
 rigid[typ2]!="" && typ1==typ2{\
  if(debug>1){print "RIGID",rigid[typ2],"DEBUG"};\
  post1=substr(id1,6);post2=substr(id2,6);;dist=w[3];\
  n=split(rigid[typ2],a," ");\
  for(i=1;i<=n;++i){\
   gsub("_"," ",a[i]);\
   oa=sprintf(" %-4s",a[i]);\
   print oa post1 "|" oa post2 "|" dist;\
  }\
 }' |\
cat >! ${t}expanded_selections.txt

if(! $twirl) then
   cp ${t}protein_selections.txt ${t}expanded_selections.txt
endif

echo "selecting solvent..."
# 2nd pass to do solvent restraints
# make list of non-taken refpoints and their neighbors
awk '{print $0,"TAKEN"}' ${t}expanded_selections.txt |\
cat - ${t}decoded_distances.txt |\
egrep -v DEBUG |\
awk -v debug=$debug '{\
      split($0,w,"|");id1=w[1];id2=w[2];\
      atom1=substr(id1,1,5);atom2=substr(id2,1,5);\
      typ1=substr(id1,7,3);typ2=substr(id2,7,3);\
      resid1=substr(id1,11);resid2=substr(id2,11);\
      ent1=atom1" "typ1" "resid1;ent2=atom2" "typ2" "resid2;\
      atm1=atom1;atm2=atom2;gsub(" ","",atm1);gsub(" ","",atm2)}\
  $NF=="TAKEN"{++seen1[ent1];++seen2[ent2];next}\
  debug>2{print $0,"DEBUG"}\
  seen1[ent1]{if(debug>1)print "seen1",ent1,"DEBUG";next} \
  seen2[ent2]{if(debug>1)print "seen2",ent2,"DEBUG";next} \
  {protein1=protein2=solvent1=solvent2=metal=0}\
  typ1~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein1}\
  typ1~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein1}\
  typ2~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein2}\
  typ2~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein2}\
  typ1~/HOH|CL |AMM|NH4|LI /{++solvent1}\
  typ2~/HOH|CL |AMM|NH4|LI /{++solvent2}\
  typ2~/ZN |EG7|SO4|EDO/ {++metal}\
  ! metal && ! protein2 {\
    print;\
    ++seen1[ent1];++seen2[ent2]}' |\
cat >! ${t}nottaken_refpoints.txt

# be sure not to use solvent refpoints too close together
# first one in a cluster of points wins
awk '{print $0,"INTRA"}' ${t}intraref_distances.txt |\
cat - ${t}nottaken_refpoints.txt |\
egrep -v DEBUG |\
awk -F "|" -v debug=$debug '{id1=$1;id2=$2;resid2=substr(id2,7)}\
  /INTRA/{nearby[id1]=nearby[id1] "|" id2;\
          nearby[id2]=nearby[id2] "|" id1;next}\
  desel[id2]{next}\
  taken[resid2]{next}\
  /DIST$/{n=split(nearby[id2],id,"|");\
    for(i=1;i<=n;++i){++desel[id[i]];\
    if(debug)print "DEBUG deselecting",id[i],resid2;\
    ++taken[resid2]}\
    print}' |\
cat >! ${t}solvent_refpoints.txt


cat ${t}expanded_selections.txt ${t}solvent_refpoints.txt |\
egrep -v DEBUG |\
tee ${t}total_selections.txt |\
awk -v debug=$debug '{\
      split($0,w,"|");id1=w[1];id2=w[2];\
      atom1=substr(id1,1,5);atom2=substr(id2,1,5);\
      typ1=substr(id1,7,3);typ2=substr(id2,7,3);\
      resid1=substr(id1,11);resid2=substr(id2,11);\
      ent1=atom1" "typ1" "resid1;ent2=atom2" "typ2" "resid2;}\
  debug>2{print $0,"DEBUG"}\
  seen2[ent2]{if(debug>1)print "seen2",ent2,"DEBUG";next} \
  {print;++seen2[ent2]}' |\
cat >! ${t}unique_selections.txt


set n = `egrep -v "DEBUG" ${t}unique_selections.txt | wc -l`
echo "unique selections: $n"
if( $n == 0 ) then
    set BAD = "gemmi failed"
    goto exit
endif

egrep -v DEBUG ${t}unique_selections.txt | head -n 3
echo "..."
egrep -v DEBUG  ${t}unique_selections.txt | tail -n 3


set Zns = `grep ZN $pdbfile | egrep "^ATOM|^HETAT" | wc -l`
if( $Zns ) then
  echo "$Zns Zn atoms in $pdbfile"
  set Zns = `grep ZN ${t}unique_selections.txt | wc -l`
  echo "$Zns Zn atoms restrained"
endif


# print out only reference atoms, named after nearest floating atoms
cat ${t}unique_selections.txt ${t}centroids_met.pdb |\
egrep -v DEBUG |\
awk -v debug=$debug '$NF=="DIST"{\
      split($0,w,"|");id1=w[1];id2=w[2];\
      dist[id2]=$(NF-1);\
      newid[id2]=id1;next}\
  {id=substr($0,12,19);\
   pre=substr($0,1,11);post=substr($0,31)}\
  ! /^ATOM|^HETAT/{;next}\
  debug>1{print "DEBUG",substr($0,7),"sel |"newid[id]"|"}\
  newid[id]{\
    print pre newid[id] post, "   ",dist[id]}' |\
cat >! ${t}centroids_near_atoms.pdb


# also get original names of selected centroids - maybe do something with them
cat ${t}unique_selections.txt ${t}centroids_met.pdb |\
egrep -v DEBUG |\
awk -v debug=$debug '$NF=="DIST"{\
      split($0,w,"|");id1=w[1];id2=w[2];\
      dist[id2]=$(NF-1);\
      newid[id2]=id1;next}\
  {id=substr($0,12,19)}\
  ! /^ATOM|^HETAT/{;next}\
  debug>1{print "DEBUG",substr($0,7),"sel |"newid[id]"|"}\
  newid[id]{\
    print $0, "   ",dist[id]}' |\
cat >! ${t}selected_centroids.pdb


# awk 'substr($0,77,2)!=" H"' ${t}centroids_near_atoms.pdb $pdbfile | rmsd | less




echo "adding weights "
# soften weight for large distances for later optimization
# weight = 1xmaxweight if at center, weight = 0.1 at 1.0 A out
# option maxweight=B to use previous B factor as maxweight
# option maxweight=last to use last word on line as weight
egrep "^CRYST" $pdbfile | head -n 1 >! ${t}weighted.pdb
echo "$weight $maxdist $softener " |\
cat - ${t}centroids_near_atoms.pdb |\
awk 'NR==1{maxweight=$1;maxdist=$2;softener=$3;next}\
 ! /^ATOM|^HETAT/{next}\
  {pre=substr($0,1,60);post=substr($0,67,80-67);Bfac=substr($0,61,6);dist=$NF;}\
  {w0=maxweight;}\
  maxweight~/^B/ {w0=Bfac;}\
  maxweight~/^last/ {w0=$(NF-1);}\
  {weight= w0/10*(10**(1-softener*dist/maxdist));\
  printf("%s%6.2f%s    %s\n",pre,weight,post,dist)}' |\
cat >>  ${t}weighted.pdb

# make sure they are in the same order as input pdb file
combine_pdbs_runme.com ${t}weighted.pdb $pdbfile \
  outfile=$outfile >! ${t}combine.log


wc $outfile



exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
rm ${t}* > /dev/null

exit

