#! /bin/tcsh -f
#
#   general-purpose distance measurements
#
#

set pdbfile = ""
set probepdb = ""
set maxdist = ""
set SG = ""
set outfile = distances.txt

set tempfile = tempfile
set debug = 0

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set ARG = `echo $Arg | awk '{print toupper($0)}'`
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
      if("$key" == "output") set outfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ [0-9]*A ) set maxdist = `echo $Arg | awk '{print $1+0}'`
      if( "$ARG" =~ [PCIFRH][1-6]* ) then
        set SG = "$ARG"
      endif
      if("$Arg" =~ *.pdb ) then
          if(! -e "$Arg" ) then
             set BAD = "$Arg does not exist."
             goto exit
          endif
          if( "$pdbfile" == "" ) then
            set pdbfile = "$Arg"
          else
            set probepdb = "$Arg"
          endif
      endif
    endif
    if("$arg" == "debug") set debug = "1"
end


if(! -e "$pdbfile") set pdbfile = refmacout.pdb
if(! -e "$probepdb") set probepdb = $pdbfile
if("$maxdist" == "") set maxdist = 3
set maxdist = `echo $maxdist 0.1 | awk '$1<$2{$1=$2} {print $1}'`

# either change it or lest ncont read the header
set symm = "symm $SG"
if("$SG" == "") set symm = ""

set t = "$tempfile"

cat << EOF
pdbfile = $pdbfile
probepdb = $probepdb
maxdist = $maxdist
SG = $SG
outfile = $outfile

tempfile = $tempfile
debug = $debug
EOF

# use or edit cell
egrep "^CRYST" $pdbfile        | head -n 1 >! ${t}temp.pdb
if("$SG" != "") then
   echo "SPACE $SG" | pdbset xyzin ${t}temp.pdb xyzout ${t}cell.pdb >! ${t}cell.log
   egrep "^CRYST" ${t}cell.pdb >! ${t}temp.pdb
endif

echo "converting atom names"
awk '{print substr($0,1,80)}' $pdbfile |\
convert_pdb.awk -v only=atoms -v skip=H,EP -v fixEe=1 -v append=origid |\
awk '{f=substr($0,17,1);pre=substr($0,1,12);atyp=substr($0,13,8);\
  atyp=substr($0,13,4)" "substr($0,18,3);\
  gsub(" ","_",atyp);mid=substr($0,22,27-22);post=substr($0,28);\
  printf("%s%s %s%s%s    \n",pre,atyp,mid,f,post)}' >> ${t}temp.pdb

awk '{print substr($0,1,80)}' $probepdb |\
convert_pdb.awk -v only=atoms -v skip=H,EP -v fixEe=1 -v append=origid |\
awk '{f=substr($0,17,1);pre=substr($0,1,12);atyp=substr($0,13,8);\
  atyp=substr($0,13,4)" "substr($0,18,3);\
  gsub(" ","_",atyp);mid=substr($0,22,27-22);post=substr($0,28);\
  printf("%s%s %s%s%s    | MODEL2\n",pre,atyp,mid,f,post)}' >> ${t}temp.pdb

convert_pdb.awk -v renumber=w4,ordinal,chain ${t}temp.pdb |\
awk '{print $0, "   "}' >! ${t}compme.pdb

echo "measuring distances..."
gemmi contact -d $maxdist ${t}compme.pdb |\
sort -k1.73g >! ${t}contact.log


echo "formatting name pairs..."
cat ${t}contact.log |\
awk -v debug=$debug '{dist=$NF}\
   {if(debug)print $0,"DEBUG";\
   id1=substr($0,12,17);\
   id2=substr($0,12+30,17);\
   print id1 "|" id2 "|  " dist,"DIST";}' |\
cat >! ${t}encoded_distances.txt

# change back to original names, excluding intra-model distances
cat ${t}compme.pdb ${t}encoded_distances.txt |\
awk -v debug=$debug '/DEBUG/{next}\
   debug{print $0,"IN"}\
   {m=1} $NF=="MODEL2"{m=2}\
   /^ATOM|^HETAT/{id=substr($0,12,17);\
   origid[id,m]=substr($0,index($0,"|")+1,17);\
   if(debug)print "DEBUG |"id"|"origid[id,m]"|";}\
   $NF!="DIST"{next}\
     {id1=substr($0,1,17);id2=substr($0,19,17);\
      dist=substr($0,19+17+1,6);\
      if(debug)print "DEBUG |"id1"|"id2"|" origid[id1,1]"|"origid[id2,2]"|"}\
    origid[id1,1]!="" && origid[id2,2]!=""{\
      print origid[id1,1] "|" origid[id2,2] "|" dist,"DIST"}' |\
cat >! $outfile

wc $outfile

# now determine which atoms to restrain to what:
# protein must have same name, solvent does not
cat $outfile |\
awk '$NF=="DIST"{\
      id1=substr($0,1,15);id2=substr($0,17,15);\
      atom1=substr(id1,1,4);atom2=substr(id2,1,4);\
      typ1=substr(id1,6,3);typ2=substr(id2,6,3);\
      resid1=substr(id1,10);resid2=substr(id2,10);\
      dist=sprintf("%6.2f",substr($0,32)+0);}\
  resid1==resid2 && typ1==typ2 && typ1=="ASP" && atom1~/^ OD/ && atom2~/^ OD/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="GLU" && atom1~/^ OE/ && atom2~/^ OE/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1~/PHE|TYR/ && atom1~/^ CD/ && atom2~/^ CD/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1~/PHE|TYR/ && atom1~/^ CE/ && atom2~/^ CE/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="ASN" && atom1~/^ [NO]D/ && atom2~/^ [NO]D/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="GLN" && atom1~/^ [NO]E/ && atom2~/^ [NO]E/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="HIS" && atom1~/^ [NC]D/ && atom2~/^ [NC]D/{atom2=atom1}\
  resid1==resid2 && typ1==typ2 && typ1=="HIS" && atom1~/^ [NC]E/ && atom2~/^ [NC]E/{atom2=atom1}\
      {ent1=atom1" "typ1" "resid1;ent2=atom2" "typ2" "resid2;}\
  seen1[ent1] || seen2[ent2]{next} \
      {protein1=protein2=solvent=metal=0}\
  typ1~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein1}\
  typ1~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein1}\
  typ2~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein2}\
  typ2~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein2}\
  typ1~/EDO|HOH|CL |AMM|NH4|SO4|LI / && typ2~/EDO|HOH|CL |AMM|NH4|SO4|LI /{++solvent}\
  typ1=="ZN " && typ2=="ZN "{++metal}\
    ent1==ent2 || solvent || metal{\
    print;\
    ++seen1[ent1];++seen2[ent2]}' |\
cat >! unique_selections.txt

