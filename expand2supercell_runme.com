#! /bin/tcsh -f
#
#  use SG and super_mult to expand a PDB and mtz into a supercell
#

#defaults
set outprefix = "super"
set pdbfile = ""
set refpdb = ""
set mtzfile = ""
set mono_map = ""
set new_mono_map = new_monomer_rot_trans.txt

set logfile = details.log

set super_mult = 1,1,1
set SG = ""

set reorganize = 1
set wrap = 1
set renumber = 1
set declash = 1
set phenix_bumpcheck = 1

set pwd = `pwd`
set pdir = `dirname $0`
set path = ( $pdir $path )

set debug = 0
set tempfile = /dev/shm/${USER}/tempfile_e2s_$$_

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
      if("$key" == "fullpdb") set refpdb = "$Val"
      if("$key" == "fulllength") set refpdb = "$Val"
      if("$key" == "rtofile") set mono_map = "$Val"
      if("$key" == "rottrans") set mono_map = "$Val"
      if("$key" == "output") set outprefix = "$Val"
      if("$key" == "outfile") set outprefix = `basename $Val .pdb`
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
      if("$Arg" =~ *.mtz ) set mtzfile = $Arg
    endif
    if("$arg" == "debug") set debug = "1"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_e2s_
if( $tempfile =~ /dev/shm/$USER/* ) mkdir -p /dev/shm/$USER/

set t = ${tempfile}

# re-interpretation
if("$renumber" == "1") set renumber = ordinal,watS,w4,chainrestart
if("$renumber" == "0") set renumber = none

foreach file ( "$pdbfile" "$mtzfile" "$refpdb" "$mono_map" )
  if(! -e "$file" && "$file" != "" ) then
     echo "WARNING: $file does not exist"
  endif
end
if(! -e "$pdbfile" && ! -e "$mtzfile") then
   set BAD = "usage: $0 asu.pdb asu.mtz super_mult=2,2,2"
   goto exit
endif
set test = `echo $super_mult | awk -F "[ ,x]" '{print $1,$2,$3}'`
if( $#test != 3 ) then
   set BAD = "use comma or x-separated list for super_mult=2,2,2"
   goto exit
endif
touch $logfile

if(-e "$pdbfile") then
  echo "extract info from asu in pdb file"
  set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $pdbfile | head -1`
  set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
  if("$SG" == "") then
      set SG = `echo $pdbSG | awk '{gsub(" ","");print}'`
      set SG = `awk -v SG=$SG '$4 == SG && $1 < 500 {print $4}' $CLIBD/symop.lib | head -1`
  endif
  if("$SG" == "") then
      set SG = `echo $pdbSG | awk '{gsub(" ","");print}'`
  endif
  set SG = `echo $SG | awk '{gsub("R","H"); print}'`
  set pdbCELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
else
  set pdbCELL = ""
endif

if(-e "$mtzfile") then

  echo "extract info from asu in mtz file"
  echo head | mtzdump hklin $mtzfile >! ${t}smallmtzdump.txt
  set SGnum = `awk '/Space group =/{print $NF+0}' ${t}smallmtzdump.txt | tail -n 1`
  set smallSG = `awk -v num=$SGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
  set smallCELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' ${t}smallmtzdump.txt`
  set nsymops = `awk -v SG=$smallSG '$4==SG {print $2}' ${CLIBD}/symop.lib | head -1`

else
  set smallSG = ""
  set smallCELL = ""
endif

if( "$smallSG" == "" && "$SG" != "") set smallSG = "$SG"
if( $#smallCELL != 6 && $#pdbCELL == 6) set smallCELL = ( $pdbCELL )

cat << EOF 
smallSG = $smallSG
smallCELL = $smallCELL
super_mult = $super_mult
outprefix = $outprefix
pdbfile = $pdbfile
refpdb = $refpdb
mtzfile = $mtzfile
mono_map = $mono_map

reorganize = $reorganize
renumber = $renumber
declash = $declash
phenix_bumpcheck = $phenix_bumpcheck

tempfile = $tempfile
debug = $debug
EOF

if(-e "$mtzfile") then
echo "expanding $mtzfile to P1"
cad hklin1 $mtzfile hklout ${t}expanded.mtz << EOF >> $logfile
labin file 1 all
outlime space 1
EOF
cad hklin1 ${t}expanded.mtz hklout ${t}refme_cell.mtz << EOF >> $logfile
labin file 1 all
symm 1
EOF
rm -f ${t}expanded.mtz

echo "further expanding mtz to $super_mult supercell"
set reidx = `echo $super_mult | awk -F "[ ,x]" '{print "reindex h"$1",k"$2",l"$3}'`
echo $reidx |\
 reindex hklin ${t}refme_cell.mtz hklout ${t}super.mtz  >> $logfile
echo head | mtzdump hklin ${t}super.mtz >! ${t}mtzdump.txt
set CELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' ${t}mtzdump.txt`

# change scale of supercell mtz?

cp ${t}super.mtz ${outprefix}.mtz
echo "${outprefix}.mtz"
endif

if(! -e "$pdbfile") goto exit

if(! $?CELL ) set CELL = ( $smallCELL )
if( $#CELL != 6 && $#smallCELL == 6) set CELL = ( $smallCELL )

echo "getting PDB version of cell dimensions"
echo "" >! ${t}blank.pdb
echo "CELL $CELL\nSPACE 1" | pdbset xyzin ${t}blank.pdb xyzout ${t}set.pdb >> $logfile
egrep "^CRYST1" ${t}set.pdb >! ${t}cell.pdb
rm ${t}blank.pdb

if(-e "$refpdb") then
    echo "using CA/OXT info in $refpdb as full sequence for renumbering"
    cat $refpdb |\
    awk '! /^ATOM|^HETAT/{next} {atom=substr($0,12,5)}\
       atom~/  CA |  OXT/{print substr($0,1,53),"  0.00  0.00          He"}' |\
    cat >! ${t}temp.pdb
    set test = `egrep "ATOM ......  OXT" ${t}temp.pdb | wc -l`
    if( $test == 0 ) then
      echo "adding dummy OXT atom"
      filter_pdb.awk -v only=protein,atoms $refpdb | tail -n 1 |\
      awk '{print "ATOM      1  OXT", substr($0,18,36),"  0.00  0.00          He"}' |\
      cat >> ${t}temp.pdb
    endif
    filter_pdb.awk -v skip=water $pdbfile >> ${t}temp.pdb
    combine_pdbs_runme.com ${t}temp.pdb $refpdb outfile=${t}dry.pdb >> $logfile
    egrep "^CRYST|^ATOM|^HETAT" ${t}dry.pdb >! ${t}padded.pdb
    filter_pdb.awk -v only=water,atoms $pdbfile >> ${t}padded.pdb
    set pdbfile = ${t}padded.pdb
    set test = `egrep "  0.00  0.00          He" ${t}padded.pdb | wc -l`
    echo "$test dummy CA/OXT atoms added"
endif

awk '{print substr($0,1,80)}' $pdbfile >! ${t}monomer.pdb

if(-e "$mono_map") then
  echo "using $mono_map"
#  set wrap = "water,salt" - salt not implemented
  set wrap = "water"

  egrep "^CRYST" ${t}cell.pdb >! ${t}pile.pdb

  set asus = `awk '{print $1}' $mono_map`
  echo "$#asus ASUs found"
  filter_pdb.awk -v only=protein ${t}monomer.pdb | egrep -v "^END" >! ${t}monolabeled.pdb
  filter_pdb.awk -v only=ligand ${t}monomer.pdb | awk '{print $0,"            SALT"}' >> ${t}monolabeled.pdb
  filter_pdb.awk -v only=water ${t}monomer.pdb  | awk '{print $0,"            SOLV"}' >> ${t}monolabeled.pdb
  foreach asu ( $asus )
    head -n $asu $mono_map | tail -n 1 | tee ${t}monomer.txt
    cat ${t}monomer.txt ${t}monolabeled.pdb |\
    awk 'NR==1{a=$1;\
          xx[a]=$2;xy[a]=$3;xz[a]=$4;\
          yx[a]=$5;yy[a]=$6;yz[a]=$7;\
          zx[a]=$8;zy[a]=$9;zz[a]=$10;\
          tx[a]=$11;ty[a]=$12;tz[a]=$13;\
          oresnum[a]=$14;ch[a]=$15;delrn[a]=$16;next}\
      /^TER/{print "TER"}\
      ! /^ATOM|^HETAT/{next}\
      /^ATOM|^HETAT/{pre=substr($0,1,21);post=substr($0,55,25);typ=substr($0,18,4);\
          solv=( $NF=="SOLV" || $NF=="SALT" );\
          salt=( $NF=="SALT" );\
          chain=substr($0,22,1);\
          resnum=substr($0,23,4)+0;\
          X=substr($0,31,8)+0;\
          Y=substr($0,39,8)+0;\
          Z=substr($0,47,8)+0;\
          newresnum=resnum;\
          if(! solv){chain=ch[a];newresnum=resnum+delrn[a];}\
          newX=xx[a]*X+xy[a]*Y+xz[a]*Z+tx[a];\
          newY=yx[a]*X+yy[a]*Y+yz[a]*Z+ty[a];\
          newZ=zx[a]*X+zy[a]*Y+zz[a]*Z+tz[a];\
         printf("%s%s%4d    %8.3f%8.3f%8.3f%s\n",pre,chain,newresnum,\
            newX,newY,newZ,post)}' |\
    cat >> ${t}pile.pdb
  end
  goto havepile
endif

echo "wrapping into small cell"
wrap_into_cell.com ${t}monomer.pdb doprotein=1 bychain=1 outfile=${t}wrapped_asu.pdb >> $logfile

echo "expanding pdb into P1"
pdbset xyzin ${t}wrapped_asu.pdb xyzout ${t}P1.pdb << EOF >> $logfile
symgen $smallSG
space 1
EOF

echo "populating pdb supercell"
awk '{print substr($0,1,80)}' ${t}P1.pdb >! ${t}exme.pdb

egrep "^CRYST" ${t}cell.pdb >! ${t}pile.pdb
foreach dx ( `echo $super_mult | awk -F "[ ,x]" '{for(x=0;x<$1;++x)print x}'` )
foreach dy ( `echo $super_mult | awk -F "[ ,x]" '{for(x=0;x<$2;++x)print x}'` )
foreach dz ( `echo $super_mult | awk -F "[ ,x]" '{for(x=0;x<$3;++x)print x}'` )
echo "shift: $dx $dy $dz"
pdbset xyzin ${t}exme.pdb xyzout ${t}shifted.pdb << EOF >> $logfile
shift frac $dx $dy $dz
EOF
egrep "^ATOM|^HETAT" ${t}shifted.pdb >> ${t}pile.pdb
end
end
end

havepile:

if( ! $reorganize ) then
   cp ${t}pile.pdb ${outprefix}.pdb
   goto exit
endif

# do not use refpdb here!
echo "reorganizing  "
reorganize_pdb_runme.com ${t}pile.pdb renumber=$renumber \
  declash=$declash outfile=${t}reorg.pdb \
  wrap=$wrap phenix_bumpcheck=$phenix_bumpcheck debug=$debug |\
tee ${t}reorg.log | grep -v ${t}

if(! -e ${t}reorg.pdb) then
   set BAD = "reorganization failed"
   goto exit
endif

if( -e "$refpdb" ) then
   echo "removing placeholder CA/OXT atoms"
   cat ${t}reorg.pdb |\
   awk '! /^ATOM|^HETAT/{print;next}\
     {atom=substr($0,12,5);occB=substr($0,55,12)}\
     atom=="  CA " && occB=="  0.00  0.00"{next}\
     atom=="  OXT" && occB=="  0.00  0.00"{print "TER";next}\
     {print}' >! ${t}stripped.pdb
    set test = `diff ${t}stripped.pdb ${t}reorg.pdb | awk '/^>/' | wc -l`
    echo "$test removed from supercell"
    cp ${t}stripped.pdb ${outprefix}.pdb
else
  cp ${t}reorg.pdb ${outprefix}.pdb
endif

if( "$new_mono_map" == "" ) goto skipmap

# make the supercell trans/rot map

echo renumber 1 | pdbset xyzin ${t}monomer.pdb xyzout ${t}ref.pdb > /dev/null

convert_pdb.awk -v only=protein -v renumber=terify ${t}reorg.pdb |\
awk '/^ATOM|^HETAT/{chain=substr($0,22,1);res=substr($0,23,4)+0;ores=$NF}\
   s==""{s=res;o=ores} {e=res}\
   /^TER/{if(s!=e)print o,chain,s,e;s=""}' |\
tee ${t}resmap.txt
set monomers = `awk '{print NR}' ${t}resmap.txt`

echo "making $new_mono_map and checking alignment..."
rm -f $new_mono_map
foreach m ( $monomers )
    set monomer = `head -n $m ${t}resmap.txt | tail -n 1`
    set o = "$monomer[1]"
    set chain = "$monomer[2]"
    set s = "$monomer[3]"
    echo "$monomer" |\
    cat - ${t}reorg.pdb |\
    awk 'NR==1{o=$1;c=$2;s=$3;e=$4;next}\
       /^ATOM|^HETAT/{chain=substr($0,22,1);res=substr($0,23,4)+0;ores=$NF}\
       c==chain && res>=s && res<=e{print substr($0,1,80)}' |\
    filter_pdb.awk -v only=protein -v CHAIN=A >! ${t}temp.pdb
    pdbset xyzin ${t}temp.pdb xyzout ${t}renum.pdb << EOF >! ${t}pdbset.log
    renumber 1
    chain A
EOF
    lsqkab xyzin1 ${t}renum.pdb xyzin2 ${t}ref.pdb << EOF >! ${t}lsq.log
FIT RESIDUE main 1 to 9999 A
MATCH            1 to 9999 A
EOF
    awk '/ROTATION MATRIX:/{getline;print;getline;print;getline;print}' ${t}lsq.log |\
    awk '{for(i=1;i<=NF;++i)$i=sprintf("%.4f",$i)+0}\
      {print}' >! ${t}rotmat.txt
    set rotat = `cat ${t}rotmat.txt`
    set trans = `awk '/VECTOR IN AS/{print $5,$6,$7}' ${t}lsq.log`
    echo  "$m  $rotat $trans  $o   $chain $s" | tee -a $new_mono_map

    pdbset xyzin ${t}ref.pdb xyzout ${t}test.pdb << EOF > /dev/null
ROTA MATR $rotat
SHIFT $trans
EOF
rmsd ${t}test.pdb ${t}renum.pdb | egrep MAXD.all | tee ${t}test.txt
set test = `awk '{print ( $2 > 0.1 )}' ${t}test.txt`
if( "$test" != "0" ) then
    set BAD = "post-alignment test failed"
    goto exit
endif

end
echo "final monomer rotation-translation and residue-number-offset map in: $new_mono_map"

skipmap:


ls -l ${outprefix}.pdb ${outprefix}.mtz

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

