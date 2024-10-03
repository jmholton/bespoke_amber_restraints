#! /bin/tcsh -f
#                                                                   -James Holton 7-14-24
#   Script for standard amber rigamrol given:
#    leap input file (will be copied and edited by command-line options)
#    original names PDB
#    any .mol2 and .frcmod files you need
#
#
if(! $?AMBERHOME) source /programs/amber22/amber.csh

set src = ${AMBERHOME}/XtalUtilities

set leapinclude = ""
if($?PHENIX) then
  if(-e ${PHENIX}/conda_base/dat/leap/cmd/) then
    set leapinclude = "-I ${PHENIX}/conda_base/dat/leap/cmd/"
  endif
endif

#set tempfile = /dev/shm/${USER}/temp_l2a$$_
set tempfile = ./temp_l2a$$_
mkdir -p /dev/shm/${USER}
mkdir -p ${CCP4_SCR}

set logfile = debuglog.log


set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"

set Stages = ""

set watertype = opc

set pdbfile = refmacout.pdb
set refpointspdb = ""
set mtzfile = refmacout.mtz
set leapfile = xtal_tleap.in
set outprefix = ""
set cut    = 9           # electrostatics switch-over distance
set dt     = 0.002       # in picoseconds
set gamma_ln    = 1.0    # thermostat coupling
set min_ns      = 0.02   # in nanoseconds
set cool_ns     = 0.2    # in nanoseconds
set heat_ns     = 2.0    # in nanoseconds
set equi_ns     = 2.0    # in nanoseconds
set prod_ns     = 100    # in nanoseconds
set ncyc   = 0   # cycles of SD before CG
set ntpr   = 100 # print out every ntpr
set ntf    = 1   # 1= no shake ; 2=shake ; 3= no bond energy ; 4= no H angle
set ntc    = 1   # 1= no shake ; 2=fixed C-H ; 3=fixed bonds 

# add and then delete dummy waters
set padwater = 0

# use reduced charges - not supported anymore?
set redq = 0


# allow providing per-atom restraints in format: atom ordinalresid weight
set restraint_file = ""
# allow keyword specification to protonate specific residues
set protons_files = ""


# defaults for single-line overall restraints
set restrainedatoms = "@C=,N=,O=,S="
set restrainedatoms = "! @H= & ! @EP="
set restraint_range = "all"

set restraint_mult = 1
set restraint_wt = 1
set pdbscale = 1


# supercell stuff
set shanB = 5
set super_mult = ( 1 1 1 )
set smallSG    = ""
set nsymops = ""

# simulation temperature target
set temperature = 287
set barostat = 1

echo "command-line arguments: $* "

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
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
    endif
    if("$key" =~ *.pdb && ! $assign ) set pdbfile = "$arg"
    if("$key" == "pdbfile") set pdbfile = "$Val"
    if("$key" =~ *.mtz && ! $assign ) set mtzfile = "$arg"
    if("$key" == "mtzfile") set mtzfile = "$Val"
    if("$Val" =~ *.pdb && "$key" == "refpoints") set refpointspdb = "$Val"
    if("$Val" =~ *.pdb && "$key" =~ cation*) set cationpdb = "$Val"
    if("$Val" =~ *.pdb && "$key" =~ anion*) set anionpdb = "$Val"
    if("$key" == "leapfile") set leapfile = "$Val"
    if("$Val" =~ *.tleap) set leapfile = "$Val"
    if("$Val" =~ *.in) set leapfile = "$Val"
    if("$key" =~ restrain* && "$key" =~ *file) set restraint_file = "$Val"
    if("$key" =~ proton*) set protons_files = ( $protons_files "$Val" )


    if("$key" == "watertype") then
      if(-e ${AMBERHOME}/dat/leap/cmd/leaprc.water.${Val} ) then
        set watertype = "$Val"
      else
        echo "WARNING: $Val is not a Valid watertype"
      endif
    endif
    if("$key" == "cut") set cut = "$Val"
    if("$key" == "dt") set dt = "$Val"
    if("$key" == "ncyc") set ncyc = "$Val"
    if("$key" == "ntpr") set ntpr = "$Val"
    if("$key" == "ntc") set ntc = "$Val"
    if("$key" == "ntf") set ntf = "$Val"
    if("$key" == "min_ns") set min_ns = "$Val"
    if("$key" == "cool_ns") set cool_ns = "$Val"
    if("$key" == "heat_ns") set heat_ns = "$Val"
    if("$key" == "equi_ns") set equi_ns = "$Val"
    if("$key" == "prod_ns") set prod_ns = "$Val"

    if("$key" == "redq") set redq = "$Val"
    if("$key" == "restraint_file") set restraint_file = "$Val"

    if("$key" == "pdbfile") set pdbfile = "$Val"
    if("$key" == "outprefix") set outprefix = "$Val"
    if("$key" == "tempfile") set tempfile = "$Val"
    if("$key" == "debug") set DEBUG = 1

    if("$key" == "pmemd") set pmemd = "$Val"
    if("$key" == "stages") set Stages = ( $Stages $Csv )
end


if( "$Stages" == "" ) then
    set Stages = ( Min0 Min Cool Heat HeatMin Equi EquiMin Prod ProdMin )
endif

# use a barostat or not
set ntb = "ntb=1,ntp=0"
if( $barostat ) then
#   set ntb = "ntb=2,ntp=1,taup=9999999999"
   set ntb = "ntb=2,ntp=4"
endif

# ntb=2 is required to monitor the pressure, but also implies changing volume
# ntb=1,ntp=0             -> periodic boundary, no pressure reporting - faster
# ntp=1                   -> isotropic pressure scaling - changes volume
# ntb=2,ntp=0             -> not allowed
# ntb=2,ntp=1,taup=999999 -> pressure regulation without actually regulating
# ntb=2,ntp=2,taup=999999 -> pressure regulation without actually regulating
# ntb=2,ntp=4             -> pressure regulation to preserve volume


if( $?DEBUG ) then
    set tempfile = tempfile
endif

# shorthand
set t = "$tempfile"


cat << EOF
Stages = $Stages
pdbfile = $pdbfile
refpoints = $refpointspdb
protons = $protons_files

redq = $redq

outprefix = $outprefix
tempfile = $tempfile
EOF

if("$refpointspdb" != "" && ! -e "$refpointspdb") then
  set BAD = "refpointspdb: $refpointspdb does not exist"
endif


# get the cell
set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
echo $CELL |\
awk 'NF==6{s=atan2(1,1)/45; A=cos(s*$4); B=cos(s*$5); G=cos(s*$6); \
 skew = 1 + 2*A*B*G - A*A - B*B - G*G ; if(skew < 0) skew = -skew;\
 printf("%.3f",$1*$2*$3*sqrt(skew))}' |\
cat >! ${t}volume
set CELLvolume = `cat ${t}volume`
rm -f ${t}volume



# strip out distractions
egrep "^CRYST|^ATOM|^HETAT|^SSBO|^LINK" $pdbfile |\
awk '{print substr($0,1,80)}' >! ${t}orig.pdb

touch ${t}refpoints.pdb
if(-e "$refpointspdb") then
    echo "combining $refpointspdb with $pdbfile for reference structure"
    combine_pdbs_runme.com $refpointspdb $pdbfile \
      printref=1 \
      outfile=${t}refpoints_in_start.pdb >> $logfile
  egrep "^CRYST|^ATOM|^HETAT|^SSBO|^LINK" ${t}refpoints_in_start.pdb |\
  awk '{print substr($0,1,80)}' >! ${t}refpoints.pdb
endif



# get the disulfides from refmac output
egrep "^SSBOND" $pdbfile |\
awk '{print $5,$8}' | sort -u | sort -g >! ${t}disulfides.txt

set test = `cat ${t}disulfides.txt | wc -l`
if("$test" == "0") then
    # do something clever?
    echo "looking for S-S bonds..."
    cat $pdbfile |\
    convert_pdb.awk -v renumber=ordinal,w4,watS,chain,chainrestart -v fixEe=1 |\
    awk '/^CRYST/{print} ! /^ATOM|^HETAT/{next}\
     {Ee=substr($0,77,2);gsub(" ","",Ee)}\
    Ee=="S"{print}' |\
    cat >! ${t}S.pdb
distang xyzin ${t}S.pdb << EOF >&! ${t}distang
SYMM 1
DIST ALL
RADII S 1.2
DMIN 0.8
END
EOF
    awk '$3==$7 && $3=="SG"{print $4,$2,$8,$6}' ${t}distang |\
    sort -u |\
    awk '! seen[$1,$2,$3,$4]{print;++seen[$1,$2,$3,$4];++seen[$3,$4,$1,$2]}' |\
    sort -k2g >! ${t}disulfides.txt
endif

cat ${t}disulfides.txt |\
awk 'NF>=4{print "bond x."$2".SG","x."$4".SG"}' |\
cat >! ${t}disulfides_tleap.in

cat ${t}disulfides.txt |\
awk 'NF>=2{print "CYX",$1 $2;print "CYX",$3 $4}' |\
cat >! ${t}disulfides_protonation.txt
set protons_files = ( $protons_files ${t}disulfides_protonation.txt )



# try to predict what tleap is going to do
# amber wants N-ternimal NH3 labelled H1 H2 H3
# refmac labels N-terminal NH3 as H H2 H3
foreach suffix ( orig refpoints )
# cat $protons_files ${t}${suffix}.pdb |\
# convert_pdb.awk -v output=amber -v norenamelig=1 \
#    -v renumber=ordinal,watS,w4,chainrestart -v fixEe=1 -v append=ordresnum |\
# egrep -v "^END" >! ${t}renumbered_${suffix}.pdb
 cat $protons_files ${t}${suffix}.pdb |\
 convert_pdb.awk -v output=amber -v norenamelig=1 \
    -v fixEe=1 -v append=ordresnum |\
 egrep -v "^END" >! ${t}renumbered_${suffix}.pdb
end





set leaprc = leaprc.protein.ff19SB
if( $redq ) then
   set places = `echo $leapinclude | awk '{$1="";print}'`
   set file = `which tleap`
   set dir  = `echo $file | awk '{gsub("bin/tleap$","");print}'`
   set dirs = `ls -1d ${dir}/*/leap/ |& grep $dir`
   set places = ( `dirname $file` $places )
   if( $?MSANDERHOME ) then
      set dirs = `ls -1d ${MSANDERHOME}/*/leap/* |& grep leap`
      set places = ( $places $dirs )
   endif
   if( $?AMBERHOME ) then
      set dirs = `ls -1d ${AMBERHOME}/*/leap/* |& grep leap`
      set places = ( $places $dirs )
   endif
   foreach place ( $places )
      echo "checking $place"
      if( -r ${place}/leaprc.ff14SB.redq ) set leaprc = leaprc.ff14SB.redq
      if( -r ${place}/oldff/leaprc.ff14SB.redq ) set leaprc = oldff/leaprc.ff14SB.redq
      if( $leaprc =~ *redq ) break
   end
   if( $leaprc =~ *redq ) then
     echo "found $leaprc in $place"
   else
     set BAD = "could not find redq leaprc"
     goto exit
   endif
endif

# use tleap to set up force field
cat << EOF >! ${t}tleap.in
source $leaprc
# reduced charges
#source leaprc.ff14SB.redq
source leaprc.water.${watertype}
EOF

if(-e "$leapfile") then
  cat $leapfile |\
  awk 'tolower($0) ~ /^x = loadpdb |^save|^quit/{$0="#"$0}\
       tolower($0) ~ /^set x /{$0="#"$0}\
       tolower($0) ~ /^source/ && /water/{$0="#"$0}\
    {print}' |\
  cat >> ${t}tleap.in
endif

# remove duplicates
cat ${t}tleap.in |\
awk '/^source/ && seen[$0]{next}\
  {print}\
  {++seen[$0]}' |\
cat >! ${t}new.in
mv ${t}new.in ${t}tleap.in

echo 'x = loadPdb "'${t}'tleapme.pdb"' >> ${t}tleap.in
cat ${t}disulfides_tleap.in >> ${t}tleap.in

cat << EOF >> ${t}tleap.in
set x box { $CELL[1] $CELL[2] $CELL[3] }
set default nocenter on
saveAmberParm x ${t}tleaped.parm7 ${t}tleaped.rst7
#savePdb x ${t}tleap0.pdb
quit
EOF


egrep -v "^END" ${t}renumbered_orig.pdb >! ${t}tleapme.pdb
if( $padwater > 0 ) then

  echo $padwater |\
  awk '{for(i=1;i<=$1;++i){\
    j=i%10000;\
    printf("ATOM      1  O   WAT z%4d       0.000   0.000   0.000  1.00  0.00\n",j);;\
    print "TER";\
    }}' >> ${t}tleapme.pdb

   # come up with names for new waters that dont conflict with provided names
   egrep "WAT|HOH" ${t}renumbered_orig.pdb >! ${t}oldwater.pdb

   egrep "WAT|HOH" ${t}tleapme.pdb |\
    convert_pdb.awk -v renumber=ordinal,watS,w4,chain,chainrestart |\
   cat >! ${t}water_autonames.pdb

   set lastwater = `tail ${t}oldwater.pdb | awk '/^ATOM|^HETAT/{print $NF}' | tail -n 1`
   grep "     0.000   0.000   0.000  1.00  0.00" ${t}water_autonames.pdb >! ${t}pad_waters.pdb

   cat ${t}oldwater.pdb ${t}pad_waters.pdb |\
   convert_pdb.awk -v dedupe=1 |\
   grep "     0.000   0.000   0.000  1.00  0.00" |\
   awk -v oresnum=$lastwater '{print $0,"       ",++oresnum}' >! ${t}new_waters.pdb

#   awk '/^ATOM|^HETAT/{print substr($0,1,80)}' ${t}water_autonames.pdb >! one.pdb
#   awk '/^ATOM|^HETAT/{print substr($0,1,80)}' ${t}unique_waters.pdb >! two.pdb
#   awk '/^ATOM|^HETAT/{print substr($0,1,80)}' ${t}oldwater.pdb ${t}new_waters.pdb >! three.pdb

endif
touch ${t}new_waters.pdb


rm -f ${t}tleaped.rst7 > /dev/null
#if(! -e leap.log) ln -sf /dev/null leap.log
echo "running tleap"
tleap -f ${t}tleap.in $leapinclude >&! ${t}tleap.out 

# make sure full cell is in there?
#ChBox -c ${t}tleaped.crd -o ${t}celled.crd -X $CELL[1] -Y $CELL[2] -Z $CELL[3] \
#  -al $CELL[4] -bt $CELL[5] -gm $CELL[6] 

if(! -e ${t}tleaped.rst7 ) then
    set BAD = "tleap failed"
    goto exit
endif

set charge = `awk '/unperturbed charge/{gsub("[)(]","");print $7}' ${t}tleap.out | tail -n 1`
echo "net charge: $charge"

set test = `awk '/Added missing heavy atom/' ${t}tleap.out | wc -l`
if( $test ) then
    echo  "WARNING: tleap changed heavy atoms"
endif

set test = `awk '/Added missing heavy atom/' ${t}tleap.out | grep OXT | wc -l`
if( $test ) then
    set BAD = "ERROR: tleap added chain breaks"
    goto exit
endif

chbox:





# write inputfile for cpptraj to extract pdbs
cat << EOF >! ${t}cpptraj_stage.in
parm ${t}xtal.prmtop
trajin ${t}tleaped.rst7 1 1
outtraj ${t}tleaped.pdb include_ep sg "P 1"
go
EOF

# convert back to PDBs
cp ${t}tleaped.parm7 ${t}xtal.prmtop
awk -v stage=tleaped '{gsub("tleaped",stage);print}' ${t}cpptraj_stage.in |\
cpptraj >> $logfile

# space group does not come through this way
#echo | cpptraj -p ${t}tleaped.parm7 -y ${t}tleaped.rst7 -ya "1 1" \
#  -x ${t}tleaped.pdb -xa 'include_ep sg "P 1"' >> $logfile

# count residues in last column
convert_pdb.awk -v fixEe=1 \
  -v append=ordresnum ${t}tleaped.pdb >! ${t}tleaped_labeled.pdb


# convert back to refmac-compatible PDBs
#awk -v stage=tleaped '{gsub("tleaped",stage);print}' ${t}cpptraj_stage.in |\
#cpptraj >> $logfile
#convert_pdb.awk -v renumber=ordinal,w4,watS,chain,chainrestart -v fixEe=1 \
#  -v append=ordresnum ${t}tleaped.pdb >! ${t}tleaped_labeled.pdb


echo "mapping original atom labels onto amber coordinates in orignames.pdb"
# map all names back to original file using xyz coordinates, extrapolating any added hydrogens
egrep "^SSBOND|^LINK" ${t}orig.pdb >! ${t}tleaped_orignames.pdb
convert_pdb.awk -v append=ordresnum ${t}orig.pdb ${t}new_waters.pdb |\
awk '/^ATOM|^HETAT/{print $0,"ORIG"}' |\
cat - ${t}tleaped_labeled.pdb |\
awk '/^TER/{print "TER";next}\
  ! /^ATOM|^HETAT/{print;next}\
     {conf=substr($0,17,1);chain=substr($0,22,1);\
      resnum=substr($0,23,8);\
      occ=substr($0,55,6);B=substr($0,61,6);\
      atom=atm=substr($0,13,4);gsub(" ","",atm);\
      type=typ=substr($0,18,4);gsub(" ","",typ);\
      xyz=substr($0,31,24);\
      x=substr($0,31,8)+0;y=substr($0,39,8)+0;z=substr($0,47,8)+0;zeroxyz=0;\
      id=atm" "typ" "chain" "resnum;}\
  x==0 && y==0 && z==0 && /O   HOH/{++zeroxyz}\
  $NF=="ORIG"{oorn=$(NF-1);orig_ordresnum[id]=oorn;repid[oorn]=id;\
      orig_chain[id]=chain;orig_resnum[id]=resnum;\
      orig_atom[id]=atom;orig_typ[id]=typ;orig_type[id]=type;\
      orig_conf[id]=conf;orig_B[id]=B;\
      if(! zeroxyz){orig_id[xyz]=id;orig_id[x,y,z]=id};next}\
    {ordresnum=$NF;pre=substr($0,1,12);rest=substr($0,67,12);\
      isH=(atm ~ /^H/ || atm == "EPW");inherit=0;\
      if(orig_id[xyz]!=""){id=orig_id[xyz]}\
      if(orig_id[x,y,z]!=""){id=orig_id[x,y,z]}\
      if( zeroxyz ){oorn=ordresnum;id=repid[oorn]}\
      if(orig_conf[id]==""){print "REMARK WARNING",id,"at",x,y,z,", orig id not found"}\
      if(orig_conf[id]!=""){oorn=orig_ordresnum[id];repid[oorn]=id}\
      if(orig_conf[id]=="" && isH && typ==orig_typ[repid[oorn]]){\
        id=repid[oorn];++inherit;\
         print "REMARK inheriting from",id}\
      if(orig_conf[id]!=""){\
        conf=orig_conf[id];B=orig_B[id];oorn=orig_ordresnum[id];\
        type=orig_type[id];chain=orig_chain[id];resnum=orig_resnum[id];\
        if(! inherit){atom=orig_atom[id]};\
      }\
      printf("%s%s%s%s%s%s%s%s%s%s %10s %10s\n",pre,atom,conf,type,chain,resnum,xyz,occ,B,rest,oorn,ordresnum)}' |\
cat >> ${t}tleaped_orignames.pdb

egrep -v "^REMARK" ${t}tleaped_orignames.pdb >! orignames.pdb

egrep "WARNING" ${t}tleaped_orignames.pdb | egrep -v "EPW|WARNING H. HOH" | egrep -v "at 0 0 0"

# rmsd $pdbfile orignames.pdb 
# awk 'substr($0,77,2)!~/ H|XP/' $pdbfile orignames.pdb | rmsd


if( $padwater > 0 ) then
  echo "stripping out padding"
  cp ${t}tleaped.parm7 ${t}padded.parm7
  cp ${t}tleaped.rst7 ${t}padded.crd
  cp orignames.pdb ${t}padded_orignames.pdb
  set nres0 = `echo list | cpptraj -p ${t}padded.parm7 | awk '$11=="res,"{print $10}' | head -n 1`
  set stripmsk = `echo $nres0 $padwater | awk '{print $1-$2+1"-"$1}'`
  cpptraj -p ${t}padded.parm7 -y ${t}padded.crd << EOF >&! ${t}strip.log
strip :$stripmsk parmout ${t}xtal.prmtop
trajout ${t}start.rst7
EOF
  if( $status ) then
      set BAD = "padding strip failed"
      goto exit
  endif

  egrep "^CRYST|^ATOM|^HETAT" ${t}padded_orignames.pdb |\
  awk -v m=$stripmsk '$NF+0>=m+0{exit} {print}' |\
  cat >! ${t}start_orignames.pdb

  # using crd messes up cell box, so...
  # check
  echo | cpptraj -p ${t}xtal.prmtop -y ${t}start.rst7 -x ${t}test.pdb >> ${t}strip.log
  echo "should be same xyz:"
  egrep "O   WAT" ${t}test.pdb | tail -n 1
  egrep "O   WAT" ${t}renumbered_orig.pdb | tail -n 1
else
  cp ${t}tleaped.parm7 ${t}xtal.prmtop
  cp ${t}tleaped.rst7 ${t}start.rst7
endif




# no need to tleap reference, we can simply convert the pdb
if(-e "$refpointspdb") then
    # convert reference model
    echo "creating ref.crd"
    echo "combining $refpointspdb with orignames.pdb for reference structure"
    combine_pdbs_runme.com $refpointspdb ${t}start_orignames.pdb \
      printref=1 \
      outfile=${t}refpoints_in_start.pdb >> $logfile
    egrep "^CRYST|^ATOM|^HETAT|^SSBO|^LINK" ${t}refpoints_in_start.pdb |\
    awk '{print substr($0,1,80)}' >! ${t}refpoints.pdb
    rm -f ${t}ref.rst7 >& /dev/null
#    cpptraj -p ${t}xtal.prmtop -y ${t}refpoints.pdb -x ${t}ref.rst7 >& ${t}refpoints_gen.log
cpptraj -p ${t}xtal.prmtop << EOF >&! ${t}refpoints_gen.log
trajin ${t}refpoints.pdb 1 1
outtraj ${t}ref.rst7
EOF
    if( $status ) then
        set BAD = "refpoints_gen failed"
        goto exit
    endif

#    cpptraj -p ${t}xtal.prmtop -y ${t}ref.rst7 -x ${t}test.pdb

#   make sure it has cell in there
#   ChBox -c ${t}ref.crd -o ${t}ref.crd -X $CELL[1] -Y $CELL[2] -Z $CELL[3] \
#     -al $CELL[4] -bt $CELL[5] -gm $CELL[6] 
endif

if(! -e ${t}ref.rst7 && -e ${t}start.rst7 ) then
    echo "using starting structure as restraint reference"
    cp -p ${t}start.rst7 ${t}ref.rst7
endif

# given rst7 name, even though it has no velocities, because otherwise cell gets messed up
ln -sf ${t}ref.rst7 ${t}ref.crd
ln -sf ${t}start.rst7 ${t}start.crd

# might as well do this now?
cp ${t}start.crd ${outprefix}start.crd
cp ${t}ref.crd ${outprefix}ref.crd
cp ${t}xtal.prmtop ${outprefix}xtal.prmtop
if(-e ${t}padded.parm7) cp ${t}padded.parm7 ${outprefix}padded.parm7






# extract reference atoms, for reference
#ln -sf ${t}ref.crd ${t}ref.rst7
awk -v stage=ref '{gsub("tleaped",stage);print}' ${t}cpptraj_stage.in |\
cpptraj >> $logfile
egrep "^SSBOND|^LINK" ${t}orig.pdb >! ${t}ref_orignames.pdb
# take coordinates only, using names from starting point
awk '/^ATOM|^HETAT/{print $0,"ORIG"}' orignames.pdb |\
cat - ${t}ref.pdb |\
awk '$NF=="ORIG"{++o;pre[o]=substr($0,1,30);post[o]=substr($0,55,length($0)-55-4);next}\
  ! /^ATOM|^HETAT/{print;next}\
      {++n;\
      printf("%s%s%s    %d\n",pre[n],substr($0,31,24),post[n],n)}' |\
cat >> ${t}ref_orignames.pdb

echo "ref vs start:"
awk 'substr($0,77,2)!=" H"' ${t}ref_orignames.pdb ${t}start_orignames.pdb |\
  rmsd | grep -v Bfac
convert_pdb.awk -v dedupe=0 -v only=protein -v skip=H ${t}ref_orignames.pdb ${t}start_orignames.pdb |\
    rmsd | awk '/MAXD.all/{print $0,"protein"}'










# now that tleap has determined order, extract residue ranges
awk '/^ATOM|^HETAT/{print $NF}' ${t}start_orignames.pdb |\
 sort -u | sort -g |\
 awk 'NR==1{s=e=$1;next} $1==e+1{e=$1;next} {print s"-"e;s=e=$1} END{print s"-"e}' |\
 sort -u | sort -g >! ${t}ranges.txt 
set ranges = `cat ${t}ranges.txt `
set all_range = `echo $ranges | awk '{gsub(" ",",");print}' `

# residue number of last protein atom    
egrep "^ATOM|^HETAT" ${t}start_orignames.pdb |\
 convert_pdb.awk -v only=protein,atoms |\
 awk '/^ATOM|^HETAT/{print $NF}' | sort -u | sort -g |\
 awk 'NR==1{s=e=$1;next} $1==e+1{e=$1;next} {print s"-"e;s=e=$1} END{print s"-"e}' |\
 sort -u | sort -g >! ${t}ranges.txt 
set ranges = `cat ${t}ranges.txt `
set protein_range = `echo $ranges | awk '{gsub(" ",",");print}' `

if(-e "$refpointspdb") then
  combine_pdbs_runme.com orignames.pdb saveBfac=1 $refpointspdb outfile=${t}restme.pdb >> $logfile
  # rmsd ${t}restme.pdb $refpointspdb

  # get residue ranges
  awk '/^ATOM|^HETAT/{print $NF}' ${t}restme.pdb |\
   sort -u | sort -g |\
   awk 'NR==1{s=e=$1;next} $1==e+1{e=$1;next} {print s"-"e;s=e=$1} END{print s"-"e}' |\
   sort -u | sort -g >! ${t}ranges.txt 
  set ranges = `cat ${t}ranges.txt `
  set refpoint_range = `echo $ranges | awk '{gsub(" ",",");print}' `
  
  if(! -e "$restraint_file") then
     # use B factors in reference points PDB as restraint weights
#     awk '/^ATOM|^HETAT/{print substr($0,12,5),$NF,substr($0,61,6)}' ${t}restme.pdb |\
#     cat >! ${outprefix}restraint_list.txt
#     set restraint_file = ${outprefix}restraint_list.txt
     set restraint_file = $refpointspdb
  endif
endif


if("$restraint_range" == "protein") then
   set restraint_range = "$protein_range"
endif
if("$restraint_range" == "all") then
   set restraint_range = "$all_range"
endif
if("$restraint_range" == "refpoint") then
   set restraint_range = "$refpoint_range"
endif





# set up restraints
if(! -e "$restraint_file") then
    echo "using single restraint weight of $restraint_wt on $restrainedatoms atoms in residues $restraint_range"
endif


set op = "$outprefix"

echo "writing input files..."
set timesteps = 100
set ntwr      = `echo 0.2 $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntwr      = 1
set ntpr      = `echo $timesteps 100 $ntwr 100 | awk '{ts=int($1/$2);nt=int($3/$4)} nt<ts{ts=nt} {print ts}'`
set ntpr = 1

cat << EOF >! ${t}Cpu.in
Minimize
 &cntrl
  imin=1,
  ntx=1,
  irest=0,
  ntf=${ntf},
  ntc=${ntc},
  maxcyc=${timesteps},
  ntmin=1,
  ncyc=0,
  ntpr=${ntpr},
  ntwr=${ntwr},
  ntwx=0,
  nsnb=1,
  ntr=1,
  restraintmask=':${restraint_range}${restrainedatoms}',
  restraint_wt=${restraint_wt},
 /
EOF



set timesteps = `echo $min_ns $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntwr      = `echo $timesteps 10 | awk '{print int($1/$2+1)}'`
set ntpr      = `echo $timesteps 100 $ntwr 100 | awk '{ts=int($1/$2);nt=int($3/$4)} nt<ts{ts=nt} {print ts}'`
set ntpr = 10

cat << EOF >! ${t}Min.in
Minimize
 &cntrl
  imin=1,
  ntx=1,
  irest=0,
  ntf=${ntf},
  ntc=${ntc},
  maxcyc=${timesteps},
  ncyc=${ncyc},
  ntpr=${ntpr},
  ntwr=${ntwr},
  ntwx=0,
  ntr=1,
  restraintmask=':${restraint_range}${restrainedatoms}',
  restraint_wt=${restraint_wt},
 /
EOF



set timesteps = `echo $cool_ns $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntwr      = `echo 0.002 $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set sdt       = `echo $dt 10 | awk '{print $1/$2}'`
set ntpr      = `echo $timesteps 10 $ntwr 10 | awk '{ts=int($1/$2);nt=int($3/$4)} nt<ts{ts=nt} {print ts}'`
set ntpr = 10

cat << EOF >! ${t}Cool.in
Cool
 &cntrl
  imin=0,
  ntx=1,
  irest=0,
  nstlim=${timesteps},
  dt=$sdt,
  ntf=2,
  ntc=2,
  tempi=0.0,
  temp0=0.0,
  ntpr=${ntpr},
  ntwx=0,
  ntwr=${ntwr},
  cut=$cut,
  ${ntb},
  !taup=99999999999,
  ntt=3,
  gamma_ln=100.0,
  nmropt=1,
  nsnb=1,
  ig=-1,
  ntr=1,
  restraintmask=':${restraint_range}${restrainedatoms}',
  restraint_wt=${restraint_wt},
 /
&wt type='TEMP0', istep1=0, istep2=${timesteps}, value1=0.0, value2=0.0 /
&wt type='END' /
EOF




set timesteps = `echo $heat_ns $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set tenth     = ( $timesteps / 10 )
set ntwr      = `echo 0.02 $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntpr      = `echo $timesteps 100 $ntwr 100 | awk '{ts=int($1/$2);nt=int($3/$4)} nt<ts{ts=nt} {print ts}'`
@ step2     = ( $tenth )
@ step3     = ( $step2 + 1 )
@ step4     = ( $timesteps - $tenth )
@ step5     = ( $step4 + 1 )

cat << EOF >! ${t}Heat.in
Heat
 &cntrl
  imin=0,
  ntx=5,
  irest=0,
  nstlim=${timesteps},
  dt=$dt,
  ntf=2,
  ntc=2,
  tempi=0.0,
  temp0=${temperature},
  ntpr=${ntpr},
  ntwx=0,
  ntwr=${ntwr},
  cut=$cut,
  ${ntb},
  !taup=99999999999,
  ntt=3,
  gamma_ln=2.0,
  nmropt=1,
  ig=-1,
  ntr=1,
  restraintmask=':${restraint_range}${restrainedatoms}',
  restraint_wt=${restraint_wt},
 /
&wt type='TEMP0', istep1=0, istep2=${step2}, value1=0.0, value2=1.0 /
&wt type='TEMP0', istep1=${step3}, istep2=${step4}, value1=1.0, value2=${temperature} /
&wt type='TEMP0', istep1=${step5}, istep2=${timesteps}, value1=${temperature}, value2=${temperature} /
&wt type='END' /
EOF


set timesteps = `echo $equi_ns $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntwx      = `echo $timesteps  10 | awk '{print int($1/$2)}'`
set ntwr      = `echo 0.2 $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntpr      = `echo $timesteps 10 $ntwr 10 | awk '{ts=int($1/$2);nt=int($3/$4)} nt<ts{ts=nt} {print ts}'`

cat << EOF >! ${t}Equi.in
Equilibrate
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  nstlim=${timesteps},
  dt=$dt,
  ntf=2,
  ntc=2,
  temp0=${temperature},
  ntpr=${ntpr},
  ntwx=${ntwx},
  ntwr=${ntwr},
  cut=$cut,
  ${ntb},
  ntt=3,
  gamma_ln=${gamma_ln},
  ig=-1,
  ntr=1,
  restraintmask=':${restraint_range}${restrainedatoms}',
  restraint_wt=${restraint_wt},
 /
EOF

set timesteps = `echo $prod_ns $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntwx      = `echo 0.2 $dt | awk '{print int(($1*1e-9)/($2*1e-12))}'`
set ntpr      = `echo $timesteps 10 $ntwx | awk '{ts=int($1/$2)} $3<ts{ts=$3} {print ts}'`

cat << EOF >! ${t}Prod.in
Production
 &cntrl
  imin=0,
  ntx=5,
  irest=1,
  nstlim=${timesteps},
  dt=$dt,
  ntf=2,
  ntc=2,
  temp0=${temperature},
  ntpr=${ntpr},
  ntwx=${ntwx},
  ntwr=${ntwx},
  cut=$cut,
  ${ntb},
  ntt=3,
  gamma_ln=${gamma_ln},
  ig=-1,
  ntr=1,
  restraintmask=':${restraint_range}${restrainedatoms}',
  restraint_wt=${restraint_wt},
 /
EOF




if(-e "$restraint_file") then
 echo "applying $restraint_mult x $pdbscale x restraints from $restraint_file"
 foreach Stage ( Cpu Min Cool Heat Equi Prod )
    if("$Stage" == "Cpu") rm -f ${outprefix}restraint_suffix.in
    echo -n "$Stage  "
    cp ${t}${Stage}.in ${t}${Stage}0.in
    if(-e ${outprefix}restraint_suffix.in) then
      egrep -v "restraint[_m]" ${t}${Stage}0.in >! ${t}${Stage}.in
      cat ${outprefix}restraint_suffix.in >> ${t}${Stage}.in
      echo "done"
    else
      restraintlist2amber.com \
        list=$restraint_file \
        orignames=orignames.pdb \
        pdbscale=$pdbscale \
        scale=$restraint_mult restraint_minwt=0.001 \
        ${t}${Stage}0.in newinfile=${t}${Stage}.in |\
           tee ${tempfile}_rl2a_${Stage}.log | grep bins
    endif
  end
endif




if(! $?MSANDERHOME) goto minimize
# create stuff for X-ray term in msander





minimize:
set laStage = start
foreach Stage ( $Stages )

  if("$Stage" =~ *Min && "$Stage" != "Min") then
      echo "copying Min.in -> ${Stage}.in"
     cp ${t}Min.in ${t}${Stage}.in
  endif

  if(! -e ${t}${Stage}.rst7 && -e ${t}${Stage}.crd) then
    ln -sf ${t}${Stage}.crd ${t}${Stage}.rst7
  endif
  if(! -e ${t}${laStage}.rst7 && -e ${t}${laStage}.crd) then
    ln -sf ${t}${laStage}.crd ${t}${laStage}.rst7
  endif

  set pmemdx = "$pmemd"
  if("$Stage" == "Cpu" ) set pmemdx = sander
try2:
  echo "running $laStage -> $Stage"
  rm -f ${t}${Stage}.rst7
  $pmemdx -O -i ${t}${Stage}.in -o ${t}${Stage}.out \
   -p ${t}xtal.prmtop \
   -c ${t}${laStage}.rst7 \
   -ref ${t}ref.crd \
   -r ${t}${Stage}.rst7 \
   -x ${t}${Stage}.nc \
   -inf ${t}${Stage}.mdinfo

  if(! -e ${t}${Stage}.rst7 ) then
      echo "did that not work? ..."
      ls -lrt ${t}* > /dev/null
      ls -l ${t}${Stage}.rst7
      if(-e ${t}${Stage}.rst7 ) then
          echo "oh, yes it did."
      endif
  endif
  if(! -e ${t}${Stage}.rst7 && "$pmemdx" != "sander" && ( $Stage =~ Cool* || $Stage =~ Min*  || $Stage =~ Cpu* ) ) then
      echo "$pmemdx failed at $Stage"
      cp -p ${t}${Stage}.in ${t}${Stage}.in.failed
      cp -p ${t}${Stage}.out ${t}${Stage}.out.failed
      echo "trying sander with dt=0.0001"
      set pmemdx = sander
      cat ${t}${Stage}.in |\
      awk '/dt=/{$0="  dt=0.0001,"}\
        /nstlim=/{$0="  nstlim=1000,"}\
        /ntwx=/{print "  ntwr=10,";$0= "  ntwx=0,"}\
        /ntpr=/{$0="  ntpr=10,"}\
        /ntwr=/{$0="  ntwr=10,"}\
        {print}' >! ${t}.x
      mv ${t}.x ${t}${Stage}.in
      goto try2
  endif

  if(! -e ${t}${Stage}.rst7) then
      set BAD = "$pmemdx failed at $Stage "
      # keep going to do all the copying
  endif
  
  if(-e ${t}${Stage}.rst7) then

    awk -v stage=${Stage} '{gsub("tleaped",stage);print}' ${t}cpptraj_stage.in |\
    cpptraj >> $logfile

    egrep "^SSBOND|^LINK|^CRYST" ${t}start_orignames.pdb >! ${t}${Stage}_orignames.pdb
    # take coordinates only, using names from starting point
    awk '/^ATOM|^HETAT/{print $0,"ORIG"}' ${t}start_orignames.pdb |\
    cat - ${t}${Stage}.pdb |\
    awk '$NF=="ORIG"{++o;pre[o]=substr($0,1,30);post[o]=substr($0,55,length($0)-55-4);next}\
      ! /^ATOM|^HETAT/{print;next}\
         {++n;\
          printf("%s%s%s    %d\n",pre[n],substr($0,31,24),post[n],n)}' |\
    cat >> ${t}${Stage}_orignames.pdb

    awk 'substr($0,77,2)!~/ H|XP/' ${t}${laStage}_orignames.pdb ${t}${Stage}_orignames.pdb |\
     rmsd | grep -v Bfac
    convert_pdb.awk -v dedupe=0 -v only=protein -v skip=H ${t}${laStage}_orignames.pdb ${t}${Stage}_orignames.pdb |\
       rmsd | awk '/MAXD.all/{print $0,"protein"}'

#    wrap_into_cell.com ${t}${Stage}_orignames.pdb outfile=${t}wrapped.pdb >> $logfile
    convert_pdb.awk \
       -v append=ordresnum ${t}${Stage}_orignames.pdb >! ${outprefix}${Stage}_orignames.pdb
    cp ${t}${Stage}.rst7 ${outprefix}${Stage}.rst7
  endif

  cp -p ${t}start_orignames.pdb ${outprefix}start_orignames.pdb
  cp -p ${t}${Stage}.in ${Stage}.in
  cp -p ${t}${Stage}.out ${outprefix}${Stage}.out
  #cp -p ${t}start.crd ${outprefix}start.crd
  #cp -p ${t}ref.crd ${outprefix}ref.crd
  if(-e ${t}${Stage}.nc) mv ${t}${Stage}.nc ${outprefix}${Stage}.nc

  cat ${outprefix}${Stage}.out |\
  awk '/  FINAL RESULTS/{++p;next}\
     / A V E R A G E S /{++p;next}\
     / F L U C T U A T /{p=0;next}\
     p && NF>1{print}\
     /EAMBER|---------------/{p=0}' |\
    tee ${t}junk.txt

  set worstatom = `awk '/NUMBER/{getline;print $NF}' ${t}junk.txt`
  echo $worstatom |\
  cat - ${outprefix}${Stage}_orignames.pdb |\
  awk 'NR==1{bad=$1;next}\
     ! /^ATOM|^HETAT/{next}\
      {++n;a=substr($0,1,13);gsub("[^0-9]","",a)}\
      a==bad || n==bad{print;exit}'

  convert_pdb.awk -v only=water -v skip=H -v dedupe=0 \
     ${outprefix}${laStage}_orignames.pdb ${outprefix}${Stage}_orignames.pdb |\
   rmsd -v debug=1 |\
   awk '/moved/{print substr($0,11,1),substr($0,12,4)+0,substr($0,25,10)}' |\
   cat >! ${outprefix}${Stage}_water_wander.txt
  convert_pdb.awk -v only=protein -v dedupe=0 -v skip=H \
     ${outprefix}${laStage}_orignames.pdb ${outprefix}${Stage}_orignames.pdb |\
  rmsd -v debug=1 |\
  awk '/moved/{dx=substr($0,25,10);\
     rn=substr($0,12,4);\
     c=substr($0,11,1);\
     print rn-offset,dx}\
     /OXT/{offset=rn}' |\
   cat >! ${outprefix}${Stage}_protein_wander.txt

  if("$Stage" !~ *Min || "$Stage" == "Min") then
     set laStage = $Stage
     #echo "laStage = $laStage"
  endif
  if("$Stage" == "Prod") set didProd

  if($?BAD) break

end

# make the solvent map?
if(-e ${outprefix}Prod.nc && $?didProd ) then
    solmapme2.com ${outprefix}Prod.nc
endif

exit:

if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?DEBUG && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    rm -f ${tempfile}*
endif

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

exit

