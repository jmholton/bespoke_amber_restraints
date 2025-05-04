#! /bin/tcsh -f
#
#   put fluffier stuff at the end of the file
#
#
set pdbfile = "$1"
set refpdb = ""
set outfile = reorganized.pdb

set renumber = ordinal,watS,w4,chainrestart
set wrap = 0
set declash = 1
set phenix_bumpcheck = 1
set ignore_zero = 1
set autorerun = 1
set debug = 0
set tempfile = /dev/shm/${USER}/tempfile_reorg_$$_
set logfile = details.log

set maxchain4resnum = ""

echo "command-line arguments: $* "

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print tolower($1)}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set key = `echo $Key | awk '{print tolower($0)}'`
    set val = `echo $Val | awk '{print tolower($0)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`

    if( $assign ) then
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "fulllength") set refpdb = "$Val"
      if("$key" == "fullpdb") set refpdb = "$Val"
    else
      # no equal sign
      if("$key" =~ *.pdb ) set pdbfile = "$Arg"
    endif
    if("$arg" == "debug") set debug = 1
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_reorg_
if( $tempfile =~ /dev/shm/$USER/* ) mkdir -p /dev/shm/$USER/

set t = "$tempfile"
touch $logfile

if( ! -e "$pdbfile" ) then
  set BAD = "pdbfile $pdbfile does not exist"
  goto exit
endif


rerun:

if(-e "$refpdb") then
    echo "using CA info in $refpdb as full sequence for renumbering"
    cat $refpdb |\
    awk 'substr($0,12,5)=="  CA "{print substr($0,1,30),"  0.000   0.000   0.000  0.00  0.00"}' |\
    cat >! ${t}temp.pdb
    filter_pdb.awk -v only=protein $pdbfile >> ${t}temp.pdb
    combine_pdbs_runme.com ${t}temp.pdb $refpdb outfile=${t}dry.pdb >> $logfile
    egrep "^CRYST|^ATOM|^HETAT" ${t}dry.pdb >! ${t}padded.pdb
    filter_pdb.awk -v only=atoms -v skip=protein $pdbfile >> ${t}padded.pdb
    set pdbfile = ${t}padded.pdb
    set test = `egrep "  0.000   0.000   0.000  0.00  0.00" ${t}padded.pdb | wc -l`
    echo "$test dummy CA atoms added"
endif


if( "$maxchain4resnum" == "" && "$renumber" != "none" ) then
  echo "determining how many chains are needed for protein"
  cat $pdbfile |\
  convert_pdb.awk -v renumber=$renumber \
   -v only=protein,atoms -v skip="H" -v append=ordresnum |\
  awk '{c=substr($0,22,1);++seen[c];ord=$NF}\
   END{for(c in seen){++n};print ord/n}' |\
  cat >! ${t}maxchain4resnum.txt
  set maxchain4resnum = `cat ${t}maxchain4resnum.txt`
  echo "determined maxchain4resnum = $maxchain4resnum"
endif


# decide what is a ligand or a salt
filter_pdb.awk -v only=ligand,salt,atoms -v skip=H $pdbfile |\
awk '{tlc=substr($0,18,3);resid=substr($0,22,5);atom=substr($0,12,5);\
    ++seen[tlc" "resid]}\
   ! aseen[atom" "tlc]{++aseen[atom" "tlc];++atoms[tlc]}\
   END{for(x in seen){\
     tlc=substr(x,1,3);\
     c=tlc;gsub(" ","",c);n=length(c);\
     print tlc,n,seen[x],atoms[tlc]}}' |\
sort -u |\
sort -k2g -k3gr -k4gr |\
awk '! seen[$1]{print $0;++seen[$1]}' |\
cat >! ${t}tlcs.txt


# always do the cell
egrep "^CRYST" $pdbfile | head -n 1 >! ${t}cell.pdb

# do SSBONDS?

# put protein first
filter_pdb.awk -v only=protein $pdbfile |\
convert_pdb.awk -v renumber=$renumber \
  -v maxchain4resnum=$maxchain4resnum -v only=protein,atoms |\
cat >! ${t}protein.pdb
echo "end of protein:"
tail -n 30 ${t}protein.pdb | convert_pdb.awk -v skip=H -v only=atoms | tail -n 3

filter_pdb.awk -v only=ligand $pdbfile >! ${t}ligands.pdb
filter_pdb.awk -v only=water $pdbfile >! ${t}water.pdb

if( "$wrap" == "1" || "$wrap" =~ "*protein*" ) then
  echo "wrapping protein into supercell..."
  cat ${t}cell.pdb ${t}protein.pdb | convert_pdb.awk -v renumber=terify >! ${t}wrapme_protein.pdb
  wrap_into_cell.com ${t}wrapme_protein.pdb byter=1 doprotein=1 outfile=${t}protein_wrap.pdb >! ${t}wrap_protein.log
  rmsd ${t}protein_wrap.pdb ${t}wrapme_protein.pdb | egrep -v "WARN|Bfac"

  cp ${t}protein.pdb ${t}protein0.pdb
  egrep -v CRYST ${t}protein_wrap.pdb >! ${t}protein.pdb

endif

if( "$wrap" == "1" || "$wrap" =~ "*nonprotein*"  || "$wrap" =~ "*notprotein*" ) then
  echo "wrapping non-protein atoms into supercell ..."
  convert_pdb.awk -v renumber=ordinal ${t}cell.pdb ${t}ligands.pdb ${t}water.pdb >! ${t}wrapme_saltwater.pdb
  wrap_into_cell.com ${t}wrapme_saltwater.pdb outfile=${t}saltwater_wrapped.pdb >! ${t}wrap_salt.log
  rmsd ${t}saltwater_wrapped.pdb ${t}wrapme_saltwater.pdb | egrep -v "WARN|Bfac"

  cp ${t}ligands.pdb ${t}ligands0.pdb
  cp ${t}water.pdb ${t}water0.pdb
  filter_pdb.awk -v only=ligand,atoms ${t}saltwater_wrapped.pdb >! ${t}ligands.pdb
  filter_pdb.awk -v only=water,atoms ${t}saltwater_wrapped.pdb >! ${t}water.pdb

endif

if( "$wrap" == "1" || "$wrap" =~ "*water*" ) then
  echo "wrapping water into supercell ..."
  convert_pdb.awk -v renumber=ordinal -v only=water ${t}cell.pdb ${t}water.pdb >! ${t}wrapme_water.pdb
  wrap_into_cell.com ${t}wrapme_water.pdb outfile=${t}water_wrapped.pdb >! ${t}wrap_water.log
  rmsd ${t}water_wrapped.pdb ${t}wrapme_water.pdb | egrep -v "WARN|Bfac"

  cp ${t}water.pdb ${t}water0.pdb
  cp ${t}water_wrapped.pdb ${t}water.pdb

endif

# put fluffy stuff at end of file
echo -n "" >! ${t}reordered.pdb
echo "ligands:"
foreach tlc ( `awk '{print NR}' ${t}tlcs.txt` )
  head -n $tlc ${t}tlcs.txt | tail -n 1 
  head -n $tlc ${t}tlcs.txt | tail -n 1 |\
  cat - ${t}ligands.pdb |\
  awk 'NR==1{sel=$1;next}\
     {tlc=substr($0,18,3);gsub(" ","",tlc)}\
      tlc==sel{print}' |\
  cat >> ${t}reordered.pdb
end
convert_pdb.awk -v only=water,atoms -v renumber=$renumber ${t}water.pdb  >> ${t}reordered.pdb

echo "adding ordinal residue numbers"
cat ${t}reordered.pdb |\
convert_pdb.awk -v renumber=$renumber -v fixEe=1 |\
cat ${t}cell.pdb ${t}protein.pdb - |\
convert_pdb.awk -v append=ordresnum >! ${t}renumbered.pdb

if( -e "$refpdb" ) then
   echo "removing placeholder CA atoms"
   cat ${t}renumbered.pdb |\
   awk '! /^ATOM|^HETAT/{print;next}\
     {atom=substr($0,12,5);xyzoccB=substr($0,30,37)}\
     atom=="  CA " && xyzoccB=="    0.000   0.000   0.000  0.00  0.00"{next}\
     {print}' >! ${t}stripped.pdb
    set test = `diff ${t}stripped.pdb ${t}renumbered.pdb | awk '/^>/' | wc -l`
    echo "$test removed"
    cp ${t}stripped.pdb ${t}renumbered.pdb
endif

if( ! $declash ) then
  cp ${t}renumbered.pdb $outfile
  goto exit
endif

if( $ignore_zero ) then
  echo "ignoring atoms with zero occ or zero B"
  cat ${t}renumbered.pdb |\
  awk '! /^ATOM|^HETAT/{print;next}\
    {occ=substr($0,55,6)+0;B=substr($0,61,6)+0}\
    occ==0 || B==0{next}\
    {print}' |\
  cat >! ${t}checkgeo.pdb
else
  cp ${t}renumbered.pdb ${t}checkgeo.pdb
endif

echo "checking for overlapping atoms..."
gemmi contact -d 1.5 --sort ${t}checkgeo.pdb >! ${t}contact.log
cat ${t}contact.log |\
awk -v debug=$debug '{dist=$NF}\
   {if(debug>1)print $0,"DEBUG IN";\
   id1=substr($0,12,19);\
   id2=substr($0,12+31,19);badid=id2}\
   substr(id1,6,3)=="HOH" && substr(id2,6,3)!="HOH"{badid=id1}\
   {print badid "|  " dist,"DELETE";}' |\
cat >! ${t}baddies.txt
set test = `cat ${t}baddies.txt | wc -l`
echo "$test clashes detected by gemmi"

if( $phenix_bumpcheck ) then
  echo "looking for clashes with phenix..."
  phenix.geometry_minimization macro_cycles=0 \
    stop_for_unknowns=false \
    output_file_name_prefix=${t}geometry ${t}checkgeo.pdb >&! ${t}phenix_geo.log
  if( $status ) then
    cat ${t}phenix_geo.log
    set BAD = geometry_minimization failed
    goto exit
  endif

  cat ${t}geometry.geo |\
  awk '/nonbonded pdb=/{key="NONBOND";split($0,w,"\"");id1=w[2];\
       getline;split($0,w,"\"");id2=w[2];\
       getline;getline;\
       obs=$1;ideal=$2;\
       if(obs+0<0.1) print id2,"DELETE"}' |\
  cat >! ${t}phenix_baddies.txt
  set test = `cat ${t}phenix_baddies.txt | wc -l`
  echo "$test clashes detected by phenix"
#  head ${t}phenix_baddies.txt ${t}baddies.txt
  cp ${t}phenix_baddies.txt ${t}baddies.txt
endif

# remove clashing waters
awk 'substr($0,6,3)=="HOH"' ${t}baddies.txt |\
cat - ${t}renumbered.pdb |\
awk '/DELETE/{id=substr($0,1,16);\
     ++del[id];next}\
  ! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,13,16)}\
  del[id]{next} {print}' |\
cat >! ${t}declashed.pdb
diff ${t}declashed.pdb ${t}renumbered.pdb | awk '/^>/' >! ${t}diff.txt
set test = `cat ${t}diff.txt | wc -l`
echo "$test clashing waters removed from supercell"

if( $test ) then
  echo "such as..."
  head ${t}diff.txt

  if( $autorerun ) then
    echo "running again on declashed result..."
    @ autorerun = ( $autorerun - 1 )
    cp ${t}declashed.pdb ${t}rerun.pdb
    set pdbfile = ${t}rerun.pdb
    goto rerun
  endif
  echo "recommend running again on $outfile"
#  cp ${t}renumbered.pdb renumbered.pdb 
else
  echo "no overlaps"
endif

cp ${t}declashed.pdb $outfile

exit:
ls -l $outfile
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if( ! $debug && "$tempfile" != "" && "$tempfile" != "./") then
   rm -f ${t}*
endif


