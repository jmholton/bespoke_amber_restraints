#! /bin/tcsh -f
#
#   bulk_solvent_runme.com                                          -James Holton 10-30-26
#
#   Quick-and-dirty explicit-solvent "bulk solvent" map from Amber MD.
#
#   Takes any PDB file (and its data, if you have them), expands it to a supercell, adds
#   salt and water, and runs MD with ONLY fixed, density-derived positional restraints for
#   the whole run - the recipe optimize_weights_runme.com uses for its alignment reference
#   (align_ref.pdb): every protein (and ligand) heavy atom sitting in at least min_rho sigma
#   of 2mFo-DFc is held where it is with weight = rho(sigma) x weight_per_sigma kcal/mol/A^2.
#   Altlocs no more than altloc_merge (0.8 A) apart are merged into one reference point; every
#   other atom with altlocs is left free.  No weight optimization, no Fo-Fc
#   feedback, no B factors: the protein is only ever as good as the model you start with.
#   The water count is tuned until the mean pressure is pressure_target (zero), then the production trajectory
#   is reduced to the average electron density of the solvent alone (water + ions),
#   collapsed into the unit cell, and written as FPART/PHIPART (electron units) for
#     phenix.refine with the phenix_user_solvent patch:
#       refinement.input.bulk_solvent_map.file_name=bulk_solvent.mtz
#     refmac:  LABIN FPART1=FPART PHIP1=PHIPART  with SOLVENT NO and SCPART
#   Every water is in this map, so refine against it a model without waters, or with only
#   the waters you are happy to count twice.
#
# usage:
#   bulk_solvent_runme.com model.pdb [data.mtz] [super_mult=2,2,2] [salt=NH4,ACY salt_conc=0.68]
#
#   data.mtz: F/SIGF give a 0-cycle phenix.refine 2mFo-DFc for the restraint weights, unless
#             it already has map coefficients (2FOFCWT or FWT).  With no mtz at all the
#             weights come from Fcalc of the model, at reso= (default 1.5 A).
#   ligands:  every HETATM residue other than water needs <RES>.mol2 and <RES>.frcmod in
#             ligdir= (default ./ligands; add <RES>.cif for phenix).  Residues without them
#             are DELETED, with a warning, and MD water fills their sites.
#   salt:     cation,anion residue names, with .mol2 .frcmod and .pdb in ligdir.  Default: no
#             salt, just Na+/Cl- counter-ions to neutralize the cell.
#   HIS:      HIE unless protons= gives HID/HIE/HIP lines in the supercell numbering.
#
# Run it in an empty directory.  Re-run the same command to resume: each section skips if its
# output already exists.  MD runs as compute_settings.sourceme (here or in ..) says.
#
set pdbfile = ""
set mtzfile = ""
# map coefficients for the restraint weights (default: FWT/2FOFCWT in mtzfile, else 0-cycle phenix.refine)
set mapmtz = ""
set reso = ""
set super_mult = 1,1,1
set ligdir = ligands
set salt = ""
set salt_conc = 0
set protons = ""
# SPC/E gave the lowest R in the 1aho tests (see 1aho_example_runme.com)
set watertype = spce

# restraints: the align_ref.pdb recipe.  weight (kcal/mol/A^2) = rho x weight_per_sigma
set restrain = protein,ligand
set min_rho = 1
set weight_per_sigma = 1
set max_weight = 9.9999
# altlocs no more than this far apart (A) are merged into their occupancy-weighted mean, which becomes the reference
# point (and the start position); atoms with altlocs further apart are not restrained
set altloc_merge = 0.8
# 1: the copies in the MD cell start in different conformers of those unrestrained split altlocs, in proportion to
# their occupancies.  0: every copy starts in the first conformer
set altloc_copies = 1
set pdbscale = 0.01

# MD: leap2amber stages, then hydration segments until the mean pressure is within the deadband of pressure_target
set stages = Cool,Heat,Equi,EquiMin,Prod
set heat_ns = 0.1
set equi_ns = 0.1
set hydrate_ns = 0.1
set hydrate_max = 20
# stop hydrating after this many segments in a row inside the deadband
set hydrate_calm = 2
# bar.  Means over 0.1 ns scatter by ~150 bar in a 1aho unit cell; 100 bar is ~3 waters there
set pressure_deadband = 100
# bar: the true pressure of the restrained cell (leap2amber's ntp=1, taup 1e11 reads true; restraint forces are not
# in the virial).  The old ntp=4 read ~N_water kT/V high, ~+500 bar in a 1aho cell (PRESSURE.md)
set pressure_target = 0
# waters added or removed per bar of pressure.  auto: kappa_T x nwater_target = 4.5e-5/bar x the waters the cell
# holds when full, worked out once from the cell (not recounted as the water changes)
set pressure_scale = auto
# waters the cell should hold.  auto: (cell volume - 1.212 A^3/Da x non-water mass) / 29.9 A^3, i.e.
# partial specific volume 0.73 mL/g for everything but water.  Pressure is no guide while the cell is
# under-filled enough to cavitate, so while the pressure is negative hydration steps at least this far.
set nwater_target = auto
set add_radius = 2.1
set add_protein_radius = 3
set prod_ns = 0.5
set nprod = 4

# solvent map: residues that count as solvent (auto: water + salt + counter-ions)
set solvent_resn = auto
set solvent_B = 10
set max_frames = 200
set outfile = bulk_solvent.mtz

set pmemd = "pmemd.cuda_SPFP"
set debug = 0
set tempfile = /dev/shm/${USER}/temp_bulksol_$$_

# site settings: here, or the site's copy one directory up
if(! -e compute_settings.sourceme && -e ../compute_settings.sourceme) then
  ln -s ../compute_settings.sourceme .
  echo "using the site's compute_settings.sourceme from the parent directory"
endif
foreach sourceme ( compute_settings.sourceme xtal_properties.sourceme user_settings.sourceme )
   if(-e $sourceme ) then
      echo "sourcing $sourceme"
      source $sourceme
   endif
end

# read the command line to update variables and other settings
foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set key = `echo $Key | awk '{print tolower($1)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = "$Val"
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "conc") set salt_conc = "$Val"
      if("$key" == "pdb") set pdbfile = "$Val"
      if("$key" == "mtz") set mtzfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = "$Arg"
      if("$Arg" =~ *.mtz ) set mtzfile = "$Arg"
      if("$arg" == "debug") set debug = 1
    endif
end
set salt = `echo $salt | awk '{gsub(","," ");print}'`
set restrain = `echo $restrain | awk '{gsub(" ",",");print}'`

if( $debug && "$tempfile" =~ /dev/shm/* ) set tempfile = ./tempfile_bulksol_
if( "$tempfile" =~ /dev/shm/* ) mkdir -p /dev/shm/${USER}
set t = $tempfile

# no Amber in the environment: try /programs/amber22 (amber.csh takes AMBERHOME from $_, which is empty in a
# non-interactive script, so source it from its own directory)
if(! $?AMBERHOME && -e /programs/amber22/amber.csh) then
  set amberhere = $cwd
  cd /programs/amber22 && source amber.csh
  cd $amberhere
endif

set pwd = `pwd`
set pdir = `dirname $0`
set pdir = `cd $pdir ; pwd`
set path = ( $pdir $path )

foreach dependency ( filter_pdb.awk convert_pdb.awk rholabel_runme.com expand2supercell_runme.com \
    reorganize_pdb_runme.com add_salt_runme.com leap2amber.com hydrate_runme.com dehydrate_amber_runme.com \
    nc2mtz_gemmi.com update_centroid_positions_runme.com tleap cpptraj AddToBox parmed gemmi mtzdmp mtzdump cad pdbset rwcontents distang )
   which $dependency > /dev/null
   if( $status ) then
       set BAD = "need $dependency in "'$'"path"
       goto exit
   endif
end

if(! -e "$pdbfile") then
  set BAD = "usage: $0 model.pdb [data.mtz] [super_mult=2,2,2] [salt=NH4,ACY salt_conc=0.68]"
  goto exit
endif
if( "$mtzfile" != "" && ! -e "$mtzfile" ) then
  set BAD = "$mtzfile does not exist"
  goto exit
endif

# ions: the salt, or by default just Na+/Cl- counter-ions from the Amber ion library that every
# leaprc.water.* loads
set counterions = 0
set ions = ( $salt )
if( "$salt" == "" ) then
  set ions = ( Na+ Cl- )
  set counterions = 1
endif
# regex for filter_pdb.awk, so salt residues are neither ligands nor restrained (noglob: it holds [+])
set noglob
set saltregex = `echo $ions | awk '{for(i=1;i<=NF;++i){s=$i;gsub("[+]","[+]",s);printf("%s%s",(i>1?"|":""),s)}}'`
unset noglob


#==============================================================================
# 1. crystal properties and data
#==============================================================================
if(-e xtal_properties.sourceme) then
  echo "=== 1: crystal properties: already done, skipping ==="
  goto sec1done
endif
echo "=== 1: crystal properties and data ==="

set smallSG = ""
set smallCELL = ""
set mtzreso = ""
if(-e "$mtzfile") then
  echo head | mtzdump hklin $mtzfile >! ${t}mtzdump.txt
  set SGnum     = `awk '/Space group =/{print $NF+0}' ${t}mtzdump.txt | tail -n 1`
  set smallSG   = `awk -v num=$SGnum '$1==num && NF>5{print $4;exit}' ${CLIBD}/symop.lib`
  set smallCELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' ${t}mtzdump.txt`
  set mtzreso   = `awk '/Resolution Range :/{getline;getline;print $(NF-2)+0;exit}' ${t}mtzdump.txt`

  # F/SIGF (plus the R-free flags, if any) for a 0-cycle refinement, as in 1aho_example_runme.com
  set F    = `mtzdmp $mtzfile | awk 'NF>10 && ! /2FOFCWT|FWT|DELF/ && $(NF-1)=="F"{print $NF;exit}'`
  set SIGF = `mtzdmp $mtzfile | awk 'NF>10 && ! /2FOFCWT|FWT|DELF/ && $(NF-1)=="Q"{print $NF;exit}'`
  set FREE = `mtzdmp $mtzfile | awk 'NF>10 && $(NF-1)=="I" && /[Ff]ree/ && $3!=$4{print $NF;exit}'`
  if( "$F" != "" && "$SIGF" != "" ) then
    echo "data columns: F=$F  SIGF=$SIGF  FREE=$FREE"
    if( "$FREE" != "" ) then
      cad hklin1 $mtzfile hklout refme_small.mtz << EOF > /dev/null
labin file 1 E1=$F E2=$SIGF E3=$FREE
labou file 1 E1=FP E2=SIGFP E3=FreeR_flag
EOF
    else
      cad hklin1 $mtzfile hklout refme_small.mtz << EOF > /dev/null
labin file 1 E1=$F E2=$SIGF
labou file 1 E1=FP E2=SIGFP
EOF
    endif
    # absent or uniform flags: phenix.refine would refuse, so make a 5% test set
    set frcol = `mtzdmp refme_small.mtz | awk '$NF=="FreeR_flag" && $(NF-1)=="I"{print $3,$4;exit}'`
    set frregen = 1
    if ( $#frcol == 2 ) then
      if ( "$frcol[1]" != "$frcol[2]" ) set frregen = 0
    endif
    if ( $frregen ) then
      echo "R-free flags absent or uniform - generating a 5% test set"
      cad hklin1 refme_small.mtz hklout ${t}noflag.mtz << EOF > /dev/null
labin file 1 E1=FP E2=SIGFP
EOF
      freerflag hklin ${t}noflag.mtz hklout refme_small.mtz << EOF >&! freerflag.log
FREERFRAC 0.05
END
EOF
      rm -f ${t}noflag.mtz
    endif
  else
    echo "no F/SIGF in ${mtzfile}, so using it for the cell and any map coefficients only"
  endif
else
  echo "no mtz: cell from $pdbfile, restraint weights from Fcalc of the model"
  # quoted space group name in symop.lib, as expand2supercell_runme.com does
  set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $pdbfile | head -1`
  set smallSG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
  if("$smallSG" == "") then
    set smallSG = `echo $pdbSG | awk '{gsub(" ","");print}'`
    set smallSG = `awk -v SG=$smallSG '$4 == SG && $1 < 500 {print $4}' $CLIBD/symop.lib | head -1`
  endif
  set smallCELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7;exit}' $pdbfile`
endif
if( "$smallSG" == "" || $#smallCELL != 6 ) then
  set BAD = "cannot work out the space group and cell (SG=$smallSG cell=$smallCELL)"
  goto exit
endif
set nsymops = `awk -v SG=$smallSG '$4==SG {print $2}' ${CLIBD}/symop.lib | head -1`
if( "$reso" == "" ) set reso = "$mtzreso"
if( "$reso" == "" ) set reso = 1.5

cat << EOF >! xtal_properties.sourceme
set reso       = $reso
set smallSG    = $smallSG
set smallCELL  = ( $smallCELL )
set nsymops    = $nsymops
set super_mult = $super_mult
set salt       = ( $salt )
set salt_conc  = $salt_conc
EOF
sec1done:
echo "cell: $smallCELL   SG: $smallSG ($nsymops ops)   reso: $reso   super_mult: $super_mult"
set ncells = `echo $super_mult | awk -F "[ ,x]" '{print $1*$2*$3}'`


#==============================================================================
# 2. single-conformer, H-free, Amber-ready ASU model
#==============================================================================
if(-e asu_noalt.pdb && -e tleap_stub.in) then
  echo "=== 2: ASU model: already done, skipping ==="
  set ligands = `cat ligands.txt`
  goto sec2done
endif
echo "=== 2: ASU model ==="

# selenomethionine as methionine: ff19SB has no MSE
awk '/^ENDMDL/{exit} /^CRYST1|^SSBOND/{print;next} ! /^ATOM|^HETAT/{next}\
  {$0=substr($0,1,80)} substr($0,18,3)=="MSE"{$0="ATOM  " substr($0,7,11) "MET" substr($0,21);\
     if(substr($0,13,4)=="SE  ")$0=substr($0,1,12) " SD " substr($0,17,60) " S"} {print}' $pdbfile |\
filter_pdb.awk -v skip=H,EP >! asu_clean.pdb

# every hetero residue needs an Amber force field; drop the ones that have none
set hetres = `filter_pdb.awk -v only=ligand,atoms asu_clean.pdb | awk '{r=substr($0,18,3);gsub(" ","",r);print r}' | sort -u`
set ligands = ""
set dropres = ""
foreach res ( $hetres )
  set issalt = `echo $res $ions | awk '{for(i=2;i<=NF;++i)if($i==$1)print 1}'`
  if( "$issalt" == "1" ) continue
  if( -e ${ligdir}/${res}.mol2 && -e ${ligdir}/${res}.frcmod ) then
    set ligands = ( $ligands $res )
  else
    set dropres = ( $dropres $res )
  endif
end
if( "$dropres" != "" ) then
  echo "WARNING: no ${ligdir}/<RES>.mol2 + .frcmod for: $dropres"
  echo "WARNING: deleting them from the model; MD water will fill their sites"
endif
echo $dropres |\
cat - asu_clean.pdb |\
awk 'NR==1{for(i=1;i<=NF;++i)drop[$i]=1;next}\
  /^ATOM|^HETAT/{r=substr($0,18,3);gsub(" ","",r)} drop[r]{next} {print}' |\
cat >! asu_kept.pdb
echo "$ligands" >! ligands.txt
echo "ligands: $ligands"
echo "salt: $ions"

# force-field files go next to the run
foreach res ( $ligands $ions )
  if( $counterions && " $ions " =~ *" $res "* ) continue
  foreach ext ( mol2 frcmod pdb cif )
    if(-e ${ligdir}/${res}.$ext && ! -e ${res}.$ext) cp ${ligdir}/${res}.$ext .
  end
  if(! -e ${res}.mol2 || ! -e ${res}.frcmod) then
    set BAD = "need ${ligdir}/${res}.mol2 and ${ligdir}/${res}.frcmod"
    goto exit
  endif
end

# one conformer, occupancy 1.  Altlocs no more than altloc_merge apart are merged into their occupancy-weighted mean
# (the reference point, and the start position: at most altloc_merge/2 from either); the others keep the first
# conformer here, and section 6 hands the other conformers out to other copies.  B = the widest separation of the
# atom's altlocs (A, 0 for one conformer), so section 4 can leave split atoms unrestrained.  Where residues have
# alternate residue types, keep the first type.
awk -v merge=$altloc_merge 'FNR==NR && ! /^ATOM|^HETAT/{next}\
  /^ATOM|^HETAT/{rid=substr($0,22,6);typ=substr($0,18,3);\
    if(rtyp[FNR==NR,rid]=="")rtyp[FNR==NR,rid]=typ;\
    if(typ!=rtyp[FNR==NR,rid])next;\
    id=substr($0,13,4) substr($0,18,10)}\
  FNR==NR{k=++n[id];x[id,k]=substr($0,31,8)+0;y[id,k]=substr($0,39,8)+0;z[id,k]=substr($0,47,8)+0;\
    o[id,k]=substr($0,55,6)+0;next}\
  /^CRYST1|^SSBOND/{print;next} ! /^ATOM|^HETAT/{next}\
  seen[id]++{next}\
  {X=x[id,1];Y=y[id,1];Z=z[id,1];d=0;\
   for(i=1;i<=n[id];++i)for(j=i+1;j<=n[id];++j){\
     r=sqrt((x[id,i]-x[id,j])^2+(y[id,i]-y[id,j])^2+(z[id,i]-z[id,j])^2);if(r>d)d=r}}\
  n[id]>1 && d<=merge{sw=sx=sy=sz=0;\
   for(i=1;i<=n[id];++i){w=o[id,i];if(w<=0)w=1e-6;sw+=w;sx+=w*x[id,i];sy+=w*y[id,i];sz+=w*z[id,i]}\
   X=sx/sw;Y=sy/sw;Z=sz/sw}\
  {if(d>999.99)d=999.99;\
   printf("%s %s%8.3f%8.3f%8.3f  1.00%6.2f%s\n",substr($0,1,16),substr($0,18,13),X,Y,Z,d,substr($0,67))}' \
  asu_kept.pdb asu_kept.pdb >! asu_noalt.pdb

set nalt = `awk '/^ATOM|^HETAT/ && substr($0,61,6)+0>0' asu_noalt.pdb | wc -l`
set nmerged = `awk -v m=$altloc_merge '/^ATOM|^HETAT/{d=substr($0,61,6)+0} d>0 && d<=m' asu_noalt.pdb | wc -l`
echo "$nalt atoms have altlocs: $nmerged no more than $altloc_merge A apart merged to one reference point, the rest unrestrained"

# force fields for tleap.  leap2amber.com adds the water model and comments out loadpdb/save/quit itself
cat << EOF >! tleap_stub.in
source leaprc.protein.ff19SB
source leaprc.water.${watertype}
source leaprc.gaff2
loadAmberParams frcmod.ff19SBmodAA
loadOff mod_amino19.lib
EOF
foreach res ( $ligands $ions )
  if( $counterions && " $ions " =~ *" $res "* ) continue
  cat << EOF >> tleap_stub.in
$res = loadMol2 ${res}.mol2
loadAmberParams ${res}.frcmod
EOF
end
sec2done:


#==============================================================================
# 3. reference density: 2mFo-DFc (or Fcalc) map coefficients Fref/PHIref
#==============================================================================
if(-e reference0.mtz) then
  echo "=== 3: reference density: already done, skipping ==="
  goto sec3done
endif
echo "=== 3: reference density ==="

set refsrc = ""
foreach mtz ( $mapmtz $mtzfile )
  if(! -e "$mtz") continue
  set Fwt   = `mtzdmp $mtz | awk 'NF>5 && /2FOFCWT|FWT/ && ! / DELF/ && $(NF-1)=="F"{print $NF;exit}'`
  set PHIwt = `mtzdmp $mtz | awk 'NF>5 && /PH2FOFCWT|PHWT/ && ! /DEL/ && $(NF-1)=="P"{print $NF;exit}'`
  if( "$Fwt" != "" && "$PHIwt" != "" ) then
    set refsrc = $mtz
    break
  endif
end
if( "$refsrc" == "" && -e refme_small.mtz ) then
  echo "0-cycle phenix.refine for a 2mFo-DFc map (details: refdensity.log)"
  set cifs = `ls -1 *.cif |& egrep "\.cif"'$'`
  phenix.refine asu_clean.pdb refme_small.mtz $cifs main.number_of_macro_cycles=0 \
     output.prefix=refdensity --overwrite >&! refdensity.log
  if(! -e refdensity_001.mtz) then
    echo "WARNING: phenix.refine failed on the full model, trying without the deleted residues"
    phenix.refine asu_kept.pdb refme_small.mtz $cifs main.number_of_macro_cycles=0 \
       output.prefix=refdensity --overwrite >&! refdensity.log
  endif
  if(-e refdensity_001.mtz) then
    set refsrc = refdensity_001.mtz
    set Fwt = 2FOFCWT
    set PHIwt = PH2FOFCWT
  else
    echo "WARNING: phenix.refine failed (see refdensity.log): restraint weights will come from Fcalc instead"
  endif
endif
if( "$refsrc" == "" ) then
  echo "Fcalc of the model at $reso A"
  awk '/^CRYST1/{print;exit}' asu_noalt.pdb >! ${t}fcme.pdb
  awk '/^ATOM|^HETAT/' asu_kept.pdb >> ${t}fcme.pdb
  gemmi sfcalc --dmin=$reso --to-mtz=modelFC.mtz ${t}fcme.pdb >&! sfcalc.log
  set refsrc = modelFC.mtz
  set Fwt = FC
  set PHIwt = PHIC
endif
echo "reference density: $Fwt $PHIwt from $refsrc"
cad hklin1 $refsrc hklout reference0.mtz << EOF > /dev/null
labin file 1 E1=$Fwt E2=$PHIwt
labou file 1 E1=Fref E2=PHIref
EOF
if(! -e reference0.mtz) then
  set BAD = "could not make reference0.mtz from $refsrc"
  goto exit
endif
sec3done:


#==============================================================================
# 4. restraint weights on the ASU: B = rho x weight_per_sigma / pdbscale
#==============================================================================
if(-e asu_labeled.pdb) then
  echo "=== 4: restraint weights: already done, skipping ==="
  goto sec4done
endif
echo "=== 4: density-derived restraint weights ==="

rholabel_runme.com asu_noalt.pdb reference0.mtz mtzlabel=Fref outfile=asu_rho.pdb >&! rholabel.log
if( $status || ! -e asu_rho.pdb ) then
  set BAD = "map probe with rholabel failed (see rholabel.log)"
  goto exit
endif

# restrained: $restrain atoms (never salt or water) in >= min_rho sigma, with one conformer or merged altlocs
filter_pdb.awk -v only=${restrain},atoms -v salt="$saltregex" asu_noalt.pdb >! ${t}restrainable.pdb
awk -v min_rho=$min_rho -v wps=$weight_per_sigma -v maxw=$max_weight \
    -v merge=$altloc_merge -v s=$pdbscale \
  'FNR==NR{if(/^ATOM|^HETAT/) ok[substr($0,13,4) substr($0,18,10)]=1;next}\
   /^CRYST1|^SSBOND/{print;next} ! /^ATOM|^HETAT/{next}\
   {id=substr($0,13,4) substr($0,18,10);rho=$NF+0;gap=substr($0,61,6)+0;w=0}\
   ok[id] && gap<=merge && rho>=min_rho{w=rho*wps; if(w>maxw)w=maxw}\
   {B=w/s; if(B>999.99)B=999.99;\
    printf("%s%6.2f%s\n",substr($0,1,60),B,substr($0,67,14))}' \
  ${t}restrainable.pdb asu_rho.pdb >! asu_labeled.pdb

set nres = `awk '/^ATOM|^HETAT/' ${t}restrainable.pdb | wc -l`
awk -v s=$pdbscale -v n0=$nres '/^ATOM|^HETAT/ && substr($0,61,6)+0>0{++n;w=substr($0,61,6)*s;sum+=w;if(w>max)max=w}\
  END{printf("%d of %d candidate atoms per ASU restrained, weight mean %.2f max %.2f kcal/mol/A^2\n",\
     n,n0,(n?sum/n:0),max)}' asu_labeled.pdb
sec4done:


#==============================================================================
# 5. cell charge, from tleap on the ASU (as 1aho_example_runme.com does)
#==============================================================================
if(-e cell_charge.txt) then
  echo "=== 5: cell charge: already done, skipping ==="
  goto sec5done
endif
echo "=== 5: cell charge ==="
mkdir -p amber_asu
cd amber_asu
foreach res ( $ligands $ions )
  if(-e ../${res}.mol2) cp ../${res}.mol2 ../${res}.frcmod .
end

# disulfides as CYX, so tleap does not build a reduced thiol and then choke on the S-S bond
egrep "^SSBOND" ../asu_noalt.pdb | awk '{print $4,$5,$7,$8}' | sort -u >! disulfides.txt
if (! -s disulfides.txt) then
  convert_pdb.awk -v fixEe=1 ../asu_noalt.pdb |\
  awk '/^CRYST/{print} ! /^ATOM|^HETAT/{next} {Ee=substr($0,77,2);gsub(" ","",Ee)} Ee=="S"{print}' >! ${t}_S.pdb
  distang xyzin ${t}_S.pdb << EOF >&! ${t}_distang.log
SYMM 1
DIST ALL
RADII S 1.2
DMIN 0.8
END
EOF
  awk '$3==$7 && $3=="SG"{print $4,$2,$8,$6}' ${t}_distang.log | sort -u |\
  awk '! seen[$1,$2,$3,$4]{print;++seen[$1,$2,$3,$4];++seen[$3,$4,$1,$2]}' >! disulfides.txt
endif
awk 'NF>=4{print "CYX",$1 $2; print "CYX",$3 $4}' disulfides.txt | sort -u >! disulfides_settings.txt

cat ../tleap_stub.in >! tleap.in
echo "x = loadpdb tleapme.pdb" >> tleap.in
awk 'NF>=4{print "bond x."$2".SG x."$4".SG"}' disulfides.txt >> tleap.in
cat << EOF >> tleap.in
set default nocenter on
charge x
saveAmberParm x asu.prmtop asu.crd
quit
EOF
cat disulfides_settings.txt ../asu_noalt.pdb |\
convert_pdb.awk -v output=amber -v fixEe=1 |\
egrep -v "^LINK|^SSBOND|^CONECT" >! tleapme.pdb
tleap -f tleap.in >&! tleap.log
set charge0 = `awk '/unperturbed charge/{gsub(/[)(:]/," ");for(i=1;i<=NF;++i)if($i~/^[-+]?[0-9]+[.][0-9]*$/){q=$i+0;print int(q+(q<0?-0.5:0.5));exit}}' tleap.log | head -n 1`
cd ..
if( ! -s amber_asu/asu.prmtop || "$charge0" == "" ) then
  set BAD = "tleap failed on the ASU - no topology or charge (see amber_asu/tleap.log)"
  goto exit
endif
set charge = `echo $charge0 $nsymops $ncells | awk '{print $1*$2*$3}'`
echo "ASU charge $charge0, supercell charge $charge" | tee cell_charge.txt
sec5done:
set charge = `awk '{print $NF;exit}' cell_charge.txt`


#==============================================================================
# 6. supercell: start model, and its restrained atoms as the reference points
#==============================================================================
if(-e start_super.pdb && -e refpoints.pdb) then
  echo "=== 6: supercell: already done, skipping ==="
  goto sec6done
endif
echo "=== 6: supercell ==="
mkdir -p supercell
cd supercell
set smallmtz = ""
if(-e ../refme_small.mtz) set smallmtz = ../refme_small.mtz

# protein stays put, ligands/waters get ordinal names (Section 8 of 1aho_example_runme.com).
# The restraint weight rides along in the B column: the 2mFo-DFc map has the crystal's symmetry,
# so every copy of an atom gets the same weight.
egrep "^CRYST|^SSBON" ../asu_labeled.pdb >! renamed.pdb
filter_pdb.awk -v only=protein,atoms ../asu_labeled.pdb >> renamed.pdb
convert_pdb.awk -v only=ligand,atoms -v renumber=ordinal,watS -v fixEe=1 ../asu_labeled.pdb >> renamed.pdb
convert_pdb.awk -v only=water,atoms -v renumber=ordinal,watS ../asu_labeled.pdb >> renamed.pdb

# twice: first to discover the monomer map, then to apply it
echo "expanding to $super_mult supercell (details: supercell/expand.log)"
expand2supercell_runme.com renamed.pdb $smallmtz super_mult=$super_mult \
  outprefix=super phenix_bumpcheck=0 >&! expand0.log
if(! -e new_monomer_rot_trans.txt) then
  set BAD = "supercell expansion failed (see supercell/expand0.log)"
  cd ..
  goto exit
endif
cp new_monomer_rot_trans.txt monomer_rot_trans.txt
expand2supercell_runme.com renamed.pdb $smallmtz super_mult=$super_mult \
  outprefix=super phenix_bumpcheck=0 mono_map=monomer_rot_trans.txt >&! expand.log

# Split altlocs (more than altloc_merge apart, unrestrained) need not start as the first conformer in every copy.
# Copy c of N gets, residue by residue, the conformer whose cumulative occupancy spans (c-0.5)/N, letters taken in
# the same order everywhere, so residues refined with the same occupancies keep the same letter in a copy.
# super.pdb holds each block of renamed.pdb (protein, ligands, water) copy by copy, with copy c at R_c x + t_c
# (monomer_rot_trans.txt) give or take a lattice vector, so the swap is x + R_c (x_conf - x_first).
set ncopy = `awk 'NF>=13' monomer_rot_trans.txt | wc -l`
set nsplit = `awk -v m=$altloc_merge '/^ATOM|^HETAT/ && substr($0,61,6)+0>m' ../asu_noalt.pdb | wc -l`
if( $altloc_copies && $ncopy > 1 && $nsplit > 0 ) then
  set nonomatch
  rm -f asu_copy_*.pdb renamed_*.pdb super_altcopies.pdb
  unset nonomatch
  awk -v merge=$altloc_merge -v ncopy=$ncopy 'FNR==NR{if(! /^ATOM|^HETAT/)next;\
      rid=substr($0,22,6);typ=substr($0,18,3);if(rtyp[rid]=="")rtyp[rid]=typ;if(typ!=rtyp[rid])next;\
      id=substr($0,13,4) substr($0,18,10);k=++n[id];res[id]=rid;alt[id,k]=substr($0,17,1);\
      x[id,k]=substr($0,31,8)+0;y[id,k]=substr($0,39,8)+0;z[id,k]=substr($0,47,8)+0;o[id,k]=substr($0,55,6)+0;next}\
    FNR==1{for(id in n){if(n[id]<2)continue;d=0;\
        for(i=1;i<=n[id];++i)for(j=i+1;j<=n[id];++j){\
          r=sqrt((x[id,i]-x[id,j])^2+(y[id,i]-y[id,j])^2+(z[id,i]-z[id,j])^2);if(r>d)d=r}\
        if(d<=merge)continue;split_[id]=1;r=res[id];\
        for(i=1;i<=n[id];++i){L=alt[id,i];if(!((r,L) in w)){let[r,++nl[r]]=L;w[r,L]=0};if(o[id,i]>w[r,L])w[r,L]=o[id,i]}}\
      for(r in nl){m=nl[r];for(i=2;i<=m;++i){v=let[r,i];for(j=i-1;j>=1 && let[r,j]>v;--j)let[r,j+1]=let[r,j];let[r,j+1]=v}\
        s=0;for(i=1;i<=m;++i){ww=w[r,let[r,i]];if(ww<=0)ww=1e-6;s+=ww}\
        for(c=1;c<=ncopy;++c){u=(c-0.5)/ncopy;cum=0;use[r,c]=let[r,m];\
          for(i=1;i<=m;++i){ww=w[r,let[r,i]];if(ww<=0)ww=1e-6;cum+=ww/s;if(u<=cum){use[r,c]=let[r,i];break}}\
          if(use[r,c]!=let[r,1])++moved[c]}\
        ++nres}\
      for(c=1;c<=ncopy;++c)line=line " " moved[c]+0;\
      printf("%d residues with split conformers; per copy, residues not in the first conformer:%s\n",nres,line)}\
    {for(c=1;c<=ncopy;++c){X=substr($0,31,8);Y=substr($0,39,8);Z=substr($0,47,8);id=substr($0,13,4) substr($0,18,10);\
       if(/^ATOM|^HETAT/ && split_[id])for(k=1;k<=n[id];++k)if(alt[id,k]==use[res[id],c]){X=x[id,k];Y=y[id,k];Z=z[id,k]}\
       if(/^ATOM|^HETAT/)printf("%s%8.3f%8.3f%8.3f%s\n",substr($0,1,30),X,Y,Z,substr($0,55)) > "asu_copy_" c ".pdb";\
       else print > "asu_copy_" c ".pdb"}}' ../asu_kept.pdb ../asu_labeled.pdb | tee altloc_copies.log
  # the renamed.pdb recipe above, so the atom order matches line for line
  foreach c ( `seq 1 $ncopy` )
    egrep "^CRYST|^SSBON" asu_copy_${c}.pdb >! renamed_${c}.pdb
    filter_pdb.awk -v only=protein,atoms asu_copy_${c}.pdb >> renamed_${c}.pdb
    convert_pdb.awk -v only=ligand,atoms -v renumber=ordinal,watS -v fixEe=1 asu_copy_${c}.pdb >> renamed_${c}.pdb
    convert_pdb.awk -v only=water,atoms -v renumber=ordinal,watS asu_copy_${c}.pdb >> renamed_${c}.pdb
  end
  set np = `filter_pdb.awk -v only=protein,atoms ../asu_labeled.pdb | egrep -c "^ATOM|^HETAT"`
  set nl = `convert_pdb.awk -v only=ligand,atoms ../asu_labeled.pdb | egrep -c "^ATOM|^HETAT"`
  # check every atom sits at R_c x + t_c (mod the cell) under the right name before moving any
  awk -v np=$np -v nl=$nl 'FILENAME=="monomer_rot_trans.txt"{if(NF<13)next;c=++nc;\
        for(m=1;m<=9;++m)R[c,m]=$(m+1);t[c,1]=$11;t[c,2]=$12;t[c,3]=$13;next}\
      FILENAME=="renamed.pdb"{if(/^ATOM|^HETAT/){i=++na;x0[i]=substr($0,31,8);y0[i]=substr($0,39,8);z0[i]=substr($0,47,8);\
        nm[i]=substr($0,13,8)};next}\
      FILENAME ~ /^renamed_[0-9]+[.]pdb$/{if(/^ATOM|^HETAT/){c=substr(FILENAME,9)+0;i=++nn[c];\
        dx[c,i]=substr($0,31,8)-x0[i];dy[c,i]=substr($0,39,8)-y0[i];dz[c,i]=substr($0,47,8)-z0[i]};next}\
      /^CRYST1/{pi=atan2(1,1)*4;a=$2;b=$3;cc=$4;ca=cos($5*pi/180);cb=cos($6*pi/180);cg=cos($7*pi/180);sg=sin($7*pi/180);\
        v=sqrt(1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg)}\
      ! /^ATOM|^HETAT/{out[++no]=$0;next}\
      {s=++ns;nw=na-np-nl;\
       if(s<=nc*np){c=int((s-1)/np)+1;i=(s-1)%np+1}\
       else if(s<=nc*(np+nl)){s-=nc*np;c=int((s-1)/nl)+1;i=np+(s-1)%nl+1}\
       else{s-=nc*(np+nl);c=int((s-1)/nw)+1;i=np+nl+(s-1)%nw+1}\
       X=substr($0,31,8);Y=substr($0,39,8);Z=substr($0,47,8);\
       ex=R[c,1]*x0[i]+R[c,2]*y0[i]+R[c,3]*z0[i]+t[c,1]-X;ey=R[c,4]*x0[i]+R[c,5]*y0[i]+R[c,6]*z0[i]+t[c,2]-Y;\
       ez=R[c,7]*x0[i]+R[c,8]*y0[i]+R[c,9]*z0[i]+t[c,3]-Z;\
       fz=ez*sg/(cc*v);fy=ey/(b*sg)+ez*(cb*cg-ca)/(b*v*sg);fx=ex/a-ey*cg/(a*sg)+ez*(ca*cg-cb)/(a*v*sg);\
       fx-=int(fx+(fx<0?-0.5:0.5));fy-=int(fy+(fy<0?-0.5:0.5));fz-=int(fz+(fz<0?-0.5:0.5));\
       ex=a*fx+b*cg*fy+cc*cb*fz;ey=b*sg*fy+cc*(ca-cb*cg)/sg*fz;ez=cc*v/sg*fz;\
       if(ex*ex+ey*ey+ez*ez>0.03^2 || substr($0,13,8)!=nm[i])++bad;\
       X+=R[c,1]*dx[c,i]+R[c,2]*dy[c,i]+R[c,3]*dz[c,i];Y+=R[c,4]*dx[c,i]+R[c,5]*dy[c,i]+R[c,6]*dz[c,i];\
       Z+=R[c,7]*dx[c,i]+R[c,8]*dy[c,i]+R[c,9]*dz[c,i];\
       if(dx[c,i]!=0 || dy[c,i]!=0 || dz[c,i]!=0)++nmoved;\
       out[++no]=sprintf("%s%8.3f%8.3f%8.3f%s",substr($0,1,30),X,Y,Z,substr($0,55))}\
    END{if(bad || ns!=nc*na){printf("WARNING: %d of %d atoms in super.pdb do not match their copy: every copy keeps the first conformer\n",bad,ns);exit}\
      for(k=1;k<=no;++k)print out[k] > "super_altcopies.pdb";\
      printf("%d atoms moved to other conformers across %d copies\n",nmoved,nc)}' \
    monomer_rot_trans.txt renamed.pdb renamed_*.pdb super.pdb | tee -a altloc_copies.log
  if( -e super_altcopies.pdb ) then
    mv super.pdb super_firstconf.pdb
    mv super_altcopies.pdb super.pdb
  endif
endif

# consolidate the copies to the chain scheme the Amber model will have (hurry Section 9 of 1aho_example_runme.com)
reorganize_pdb_runme.com super.pdb outfile=consolidated.pdb \
  phenix_bumpcheck=0 declash=0 debug=0 >&! reorg.log
cd ..
if(! -e supercell/consolidated.pdb) then
  set BAD = "supercell consolidation failed (see supercell/expand.log and supercell/reorg.log)"
  goto exit
endif
cp supercell/consolidated.pdb start_super.pdb

awk '/^CRYST1/{print} /^ATOM|^HETAT/ && substr($0,61,6)+0>0' start_super.pdb >! refpoints.pdb
set nasu = `awk '/^ATOM|^HETAT/ && substr($0,61,6)+0>0' asu_labeled.pdb | wc -l`
set nref = `awk '/^ATOM|^HETAT/' refpoints.pdb | wc -l`
set expect = `echo $nasu $nsymops $ncells | awk '{print $1*$2*$3}'`
echo "$nref reference points in the supercell ( $nasu per ASU x $nsymops x $ncells = $expect )"
if( $nref == 0 ) then
  set BAD = "no reference points survived the supercell expansion"
  goto exit
endif
sec6done:


#==============================================================================
# 7. salt and water (Section 10 of 1aho_example_runme.com)
#==============================================================================
if(-e amberme.pdb && -e padwater.txt) then
  echo "=== 7: salt and water: already done, skipping ==="
  goto sec7done
endif
echo "=== 7: salt and water ==="
if( $counterions ) then
  # Amber ion library names; charge in columns 79-80 for add_salt_runme.com
  echo "HETATM    1 Na+  Na+ A   1       0.000   0.000   0.000  1.00 20.00          NA+1" >! cation.pdb
  echo "HETATM    1 Cl-  Cl- A   1       0.000   0.000   0.000  1.00 20.00          CL-1" >! anion.pdb
  set cation = Na+
  set anion = Cl-
else
  foreach ion ( $ions )
    set icharge = `awk '/TRIPOS.ATOM/,/TRIPOS.BOND/{sum+=$NF} END{printf("%.1g",sum)}' ${ion}.mol2`
    set ispos = `echo $icharge | awk '{print ($1 > 0)}'`
    echo "$ion charge $icharge"
    if ($ispos) then
      set cation = "$ion"
      ln -sf ${ion}.pdb cation.pdb
    else
      set anion = "$ion"
      ln -sf ${ion}.pdb anion.pdb
    endif
  end
  if(! $?cation || ! $?anion) then
    set BAD = "salt=$ions needs one cation and one anion"
    goto exit
  endif
endif

echo "adding $cation $anion at $salt_conc M to neutralize $charge, and water (details: add_salt.log)"
add_salt_runme.com start_super.pdb conc=$salt_conc watertype=$watertype \
  RIP=4 RIW=3 charge=$charge anion=$anion cation=$cation >&! add_salt.log
if (! -e salty.pdb) then
  set BAD = "add_salt_runme.com did not produce salty.pdb (see add_salt.log)"
  goto exit
endif
egrep "new charge|new waters" add_salt.log

reorganize_pdb_runme.com salty.pdb ignore_zero=0 refpdb=start_super.pdb \
  outfile=reorganized.pdb phenix_bumpcheck=0 declash=1 >&! reorganize_final.log
if (! -e reorganized.pdb) then
  set BAD = "reorganize failed (see reorganize_final.log)"
  goto exit
endif
filter_pdb.awk -v skip=water reorganized.pdb | egrep -v "^END" >! amberme.pdb
filter_pdb.awk -v skip=H -v only=water,atoms reorganized.pdb >> amberme.pdb

# padding waters, so hydration can add up to 150% of the cell's water capacity
egrep -v HOH amberme.pdb >! ${t}dry.pdb
echo | rwcontents xyzin ${t}dry.pdb >! ${t}rwcontents.log
set waterslots = `awk '/Cell volume/{V=$NF;n=55*6.022e23/1e27*V} /% of cell without atoms/{print int(n*$NF/100),int(n)}' ${t}rwcontents.log`
set gotwater = `grep "O   HOH" amberme.pdb | wc -l`
set padwater = `echo $waterslots $gotwater | awk '{v=$1*1.5;t=$2-$3;print 0+sprintf("%.2g",(v>t?v:t))}'`
echo "room for $waterslots[1] waters ($waterslots[2] in the whole cell), have $gotwater, padwater = $padwater"
echo $padwater >! padwater.txt
sec7done:
set padwater = `cat padwater.txt`


#==============================================================================
# 8. Amber system and the warm-up MD (leap2amber.com)
#==============================================================================
if(-e Prod.rst7 && -e Prod.in) then
  echo "=== 8: leap2amber: already done, skipping ==="
  goto sec8done
endif
echo "=== 8: leap2amber stages $stages (details: leap2amber.log) ==="
touch protonation.txt
if(-e "$protons") cp $protons protonation.txt
leap2amber.com amberme.pdb stages=$stages declash_waters=1 \
  protons=protonation.txt watertype=$watertype flexwater=0 \
  refpoints=refpoints.pdb restraint_mult=1 \
  pdbscale=$pdbscale gamma_ln=1.0 barostat=1 \
  leapfile=tleap_stub.in padwater=$padwater \
  cool_ns=0.001 heat_ns=$heat_ns equi_ns=$equi_ns prod_ns=$hydrate_ns \
  cool_slowdown=5 heat_slowdown=1 equi_slowdown=1 \
  restrain_omega=0 omega_weight=0 chiral_weight=0 \
  debug=0 >&! leap2amber.log
egrep "net charge|ERROR|clash|inverted" leap2amber.log | head -20
if (! -e Prod.rst7 || ! -e Prod.in) then
  set BAD = "leap2amber.com did not produce Prod.rst7 (see leap2amber.log)"
  goto exit
endif
set test = `awk '/^net charge:/{q=$NF+0} END{print (q*q>0.25)}' leap2amber.log`
if( $test ) echo "WARNING: the system is not neutral (protons= HIP/ASH/GLH are not in the ASU charge count)"
sec8done:

# MD inputs from leap2amber's Prod.in, which carries the restraint groups.  Water is added and
# removed at the end of the atom list, so the restrained atom numbers never change.
set dt = `awk -F "=" '/ dt=/{print $2+0;exit}' Prod.in`
foreach stage ( hydrate prod )
  set ns = $hydrate_ns
  if( $stage == prod ) set ns = $prod_ns
  set ntwx = 0
  if( $stage == prod ) set ntwx = `awk -F "=" '/ ntwx=/{print $2+0;exit}' Prod.in`
  echo $ns $dt $ntwx |\
  cat - Prod.in |\
  awk 'NR==1{nstlim=int($1/$2*1000+0.5);ntwx=$3;ntpr=int(nstlim/100);if(ntpr<1)ntpr=1;next}\
    / nstlim=/{print "  nstlim=" nstlim ",";next}\
    / ntpr=/{print "  ntpr=" ntpr ",";next}\
    / ntwx=/{print "  ntwx=" ntwx ",";next}\
    / ntwr=/{print "  ntwr=" nstlim ",";next}\
    {print}' >! ${stage}.in
end
# after waters change: fresh velocities at 1/40 time step with every-step pair lists, then 1/4 step
cat Prod.in |\
awk '/ irest=/{print "  irest=0,";next}\
  / ntx=/{print "  ntx=1,";next}\
  / dt=/{print "  dt=0.00005,";next}\
  / nstlim=/{print "  nstlim=2000,";next}\
  / ntpr=| ntwr=/{split($1,w,"=");print "  " w[1] "=2000,";next}\
  / ntwx=/{print "  ntwx=0,";next}\
  / temp0=/{print;split($1,w,"=");print "  tempi=" w[2];next}\
  / cut=/{print;print "  nsnb=1,";next}\
  {print}' >! settle.in
cat Prod.in |\
awk '/ dt=/{print "  dt=0.0005,";next}\
  / nstlim=/{print "  nstlim=10000,";next}\
  / ntpr=| ntwr=/{split($1,w,"=");print "  " w[1] "=10000,";next}\
  / ntwx=/{print "  ntwx=0,";next}\
  {print}' >! resettle.in


#==============================================================================
# 9. hydration: add or remove water until the mean pressure stays within the deadband
#==============================================================================
if(-e hydrated.rst7 && -e hydrated.prmtop) then
  echo "=== 9: hydration: already done, skipping ==="
  goto sec9done
endif
echo "=== 9: hydration (target $pressure_target bar, deadband +/- $pressure_deadband bar) ==="

# resume from the last finished segment
set h = 1
set last = Prod
while ( -e hydrate_${h}.rst7 && -e hydrate_${h}.prmtop && -e hydrate_${h}_ref.crd )
  set last = hydrate_${h}
  @ h = ( $h + 1 )
end
if( "$last" != "Prod" ) then
  echo "resuming after $last"
  cp ${last}.prmtop xtal.prmtop
  cp ${last}_ref.crd ref.crd
endif
if(! -e pressure_vs_seg.txt) touch pressure_vs_seg.txt
set test = `awk -v s=$last '$1==s' pressure_vs_seg.txt | wc -l`
if( ! $test ) then
  set press = `awk '/PRESS/ && p{n=$3;print $NF} /A V E R A/{++p} END{print n+0}' ${last}.out | tail -n 3`
  set nwater = `awk '/^%FLAG RESIDUE_LABEL/{f=1;getline;next} /^%FLAG/{f=0} f{for(i=1;i<=length($0);i+=4){r=substr($0,i,4);gsub(/ /,"",r);n+=(r=="WAT"||r=="HOH")}} END{print n+0}' xtal.prmtop`
  echo "$last $press $nwater" | tee -a pressure_vs_seg.txt
endif
# columns: segment  <P>  rms(P)  nsteps  Nwater

if( "$nwater_target" == "auto" ) then
  set cellvol = `awk '/A V E R A/{p=1} p && /VOLUME/{print $NF;exit}' Prod.out`
  set nwater_target = `awk -v V=$cellvol '/^%FLAG/{f=$2;getline;next}\
    f=="MASS"{for(i=1;i<=NF;++i)m[++na]=$i+0}\
    f=="RESIDUE_LABEL"{for(i=1;i<=length($0);i+=4){r=substr($0,i,4);gsub(/ /,"",r);if(r!="")lab[++nr]=r}}\
    f=="RESIDUE_POINTER"{for(i=1;i<=NF;++i)ptr[++np]=$i}\
    END{ptr[nr+1]=na+1;for(r=1;r<=nr;++r)if(lab[r]!="WAT" && lab[r]!="HOH")for(a=ptr[r];a<ptr[r+1];++a)M+=m[a];\
        printf("%d\n",(V-1.212*M)/29.9)}' xtal.prmtop`
  echo "cell volume $cellvol A^3 should hold about $nwater_target waters"
endif
set ps = `echo $pressure_scale $nwater_target | awk '$1+0>0{print $1+0;exit} {print 4.5e-5*$2}'`
echo "pressure_scale: $ps waters/bar"

while ( 1 )
  set press = `awk -v s=$last '$1==s{print $2}' pressure_vs_seg.txt | tail -n 1`
  set nwater = `awk -v s=$last '$1==s{print $5}' pressure_vs_seg.txt | tail -n 1`
  set calm = `awk -v d=$pressure_deadband -v p0=$pressure_target '{c=(($2-p0)^2<=d*d)?c+1:0} END{print c+0}' pressure_vs_seg.txt`
  if( $calm >= $hydrate_calm ) then
    echo "pressure inside $pressure_target +/- $pressure_deadband bar for $calm segments: hydrated"
    break
  endif
  if( $h > $hydrate_max ) then
    echo "WARNING: pressure still $press bar after $hydrate_max hydration segments - going on anyway"
    break
  endif

  # kappa_T step, or up to the target count if that is further (only while under-filled: once the
  # pressure is positive the liquid is merely compressed, and kappa_T is right)
  set dw = `echo $press $ps $pressure_deadband $nwater $nwater_target $pressure_target | awk '{$1-=$6} $1*$1<=$3*$3{print 0;exit}\
    {d=-$1*$2;t=$5-$4} $1<0 && t>d{d=t} {printf("%.0f\n",d)}'`
  echo "segment ${h}: P = $press bar with $nwater waters (target $nwater_target), $ps waters/bar -> change by $dw"

  set laStage = $last
  if( $dw > 0 ) then
    rm -f wetter.rst7 >& /dev/null
    hydrate_runme.com ${last}.rst7 outrst=wetter.rst7 nadd=$dw water_radius=$add_radius \
      protein_radius=$add_protein_radius restraints=refpoints.pdb debug=$debug >&! hydrate_${h}.log
    if( $status || ! -e wetter.rst7 ) then
      set BAD = "hydration failed (see hydrate_${h}.log)"
      goto exit
    endif
    egrep "^added " hydrate_${h}.log
    set laStage = wetter
  endif
  if( $dw < 0 ) then
    # least-needed waters: the most negative 2mFo-DFc
    set nrem = `echo $dw | awk '{print -$1}'`
    rm -f drier.parm7 drier.rst7 drierref.crd >& /dev/null
    dehydrate_amber_runme.com xtal.prmtop ${last}.rst7 maxreject=$nrem \
      mtzfile=reference0.mtz mtzlabel=Fref notthese=none >&! dehydrate_${h}.log
    if( $status || ! -e drier.parm7 || ! -e drier.rst7 ) then
      set BAD = "dehydration failed (see dehydrate_${h}.log)"
      goto exit
    endif
    egrep "eliminating|nothing to do" dehydrate_${h}.log
    mv drier.parm7 xtal.prmtop
    # rebuild ref.crd from the reference points, as optimize_weights_runme.com does
    update_centroid_positions_runme.com refpoints.pdb outfile=ref.crd >> dehydrate_${h}.log
    if( $status ) then
      set BAD = "could not rebuild ref.crd after dehydration (see dehydrate_${h}.log)"
      goto exit
    endif
    set laStage = drier
  endif
  if( "$laStage" != "$last" ) then
    foreach settle ( settle resettle )
      $pmemd -O -i ${settle}.in -o ${settle}.out -p xtal.prmtop -c ${laStage}.rst7 -ref ref.crd \
        -r ${settle}.rst7 -inf ${settle}.mdinfo
      if( $status || ! -e ${settle}.rst7 ) then
        set BAD = "MD failed at $settle after changing the water (see ${settle}.out)"
        goto exit
      endif
      set laStage = $settle
    end
  endif

  echo "running hydrate_${h} ( $hydrate_ns ns )"
  $pmemd -O -i hydrate.in -o hydrate_${h}.out -p xtal.prmtop -c ${laStage}.rst7 -ref ref.crd \
    -r hydrate_${h}.rst7 -inf hydrate.mdinfo
  set mdstatus = $status
  set test = `egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' hydrate_${h}.out | egrep -c 'NaN|\*\*\*\*\*'`
  if( $mdstatus || ! -e hydrate_${h}.rst7 || $test ) then
    set BAD = "MD failed at hydrate_${h} (see hydrate_${h}.out)"
    goto exit
  endif
  cp xtal.prmtop hydrate_${h}.prmtop
  cp ref.crd hydrate_${h}_ref.crd
  set last = hydrate_${h}
  set press = `awk '/PRESS/ && p{n=$3;print $NF} /A V E R A/{++p} END{print n+0}' ${last}.out | tail -n 3`
  set nwater = `awk '/^%FLAG RESIDUE_LABEL/{f=1;getline;next} /^%FLAG/{f=0} f{for(i=1;i<=length($0);i+=4){r=substr($0,i,4);gsub(/ /,"",r);n+=(r=="WAT"||r=="HOH")}} END{print n+0}' xtal.prmtop`
  echo "$last $press $nwater" | tee -a pressure_vs_seg.txt
  @ h = ( $h + 1 )
end
cp ${last}.rst7 hydrated.rst7
cp xtal.prmtop hydrated.prmtop
cp ref.crd hydrated_ref.crd
sec9done:


#==============================================================================
# 10. production
#==============================================================================
echo "=== 10: production: $nprod x $prod_ns ns ==="
cp hydrated.prmtop xtal.prmtop
cp hydrated_ref.crd ref.crd
set last = hydrated
foreach i ( `seq 1 $nprod` )
  if(-e prod_${i}.rst7 && -e prod_${i}.nc) then
    set last = prod_${i}
    continue
  endif
  echo "running prod_${i}"
  $pmemd -O -i prod.in -o prod_${i}.out -p xtal.prmtop -c ${last}.rst7 -ref ref.crd \
    -r prod_${i}.rst7 -x prod_${i}.nc -inf prod.mdinfo
  set mdstatus = $status
  set test = `egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' prod_${i}.out | egrep -c 'NaN|\*\*\*\*\*'`
  if( $mdstatus || ! -e prod_${i}.rst7 || ! -e prod_${i}.nc || $test ) then
    set BAD = "MD failed at prod_${i} (see prod_${i}.out)"
    goto exit
  endif
  set press = `awk '/PRESS/ && p{n=$3;print $NF} /A V E R A/{++p} END{print n+0}' prod_${i}.out | tail -n 3`
  echo "prod_${i} $press" | tee -a pressure_vs_seg.txt
  set last = prod_${i}
end


#==============================================================================
# 11. average solvent density -> FPART/PHIPART
#==============================================================================
if(-e $outfile) then
  echo "=== 11: $outfile already exists, done ==="
  goto exit
endif
echo "=== 11: solvent map ==="

if( "$solvent_resn" == "auto" ) set solvent_resn = ( WAT HOH $ions )
set solvent_resn = `echo $solvent_resn | awk '{gsub(","," ");print}'`
# cpptraj masks: + and - are operators, so match ion names like Na+ with the = wildcard
set solvmask = `echo $solvent_resn | awk '{for(i=1;i<=NF;++i){r=$i;gsub("[-+]","=",r);printf("%s%s",(i>1?",":":"),r)}}'`

set ncs = `ls -1 prod_*.nc |& awk -v n=$nprod -F "[_.]" '/^prod_[0-9]*[.]nc$/ && $2+0<=n' | sort -t_ -k2g`
if( "$ncs" == "" ) then
  set BAD = "no production trajectories"
  goto exit
endif
set nframes = `echo $#ncs $prod_ns | awk '{print int($1*$2/0.02+0.5)}'`
set stride = `echo $nframes $max_frames | awk '{s=int(($1+$2-1)/$2);print (s<1?1:s)}'`
echo "solvent = $solvent_resn ( $solvmask ), $nframes frames, using every $stride"

rm -f solvent.nc solvent.parm7 >& /dev/null
echo "parm xtal.prmtop" >! ${t}strip.in
foreach nc ( $ncs )
  echo "trajin $nc 1 last $stride" >> ${t}strip.in
end
cat << EOF >> ${t}strip.in
strip !($solvmask) parmout solvent.parm7
trajout solvent.nc
go
EOF
cpptraj < ${t}strip.in >&! strip_solvent.log
if(! -e solvent.nc || ! -e solvent.parm7) then
  set BAD = "cpptraj could not extract the solvent (see strip_solvent.log)"
  goto exit
endif

echo "rendering solvent density (details: nc2mtz_solvent.log)"
nc2mtz_gemmi.com solvent.nc topfile=solvent.parm7 paddedparm=none \
  reso=$reso B=$solvent_B outmap=solvent_avg.map outfile=solvent_avg.mtz >&! nc2mtz_solvent.log
if(! -e solvent_avg.mtz) then
  set BAD = "nc2mtz_gemmi.com failed (see nc2mtz_solvent.log)"
  goto exit
endif

cad hklin1 solvent_avg.mtz hklout $outfile << EOF > /dev/null
labin file 1 E1=FCavg E2=PHICavg
labou file 1 E1=FPART E2=PHIPART
EOF
if(! -e $outfile) then
  set BAD = "could not write $outfile"
  goto exit
endif

# the cell-average density over 0.334 e/A^3 (water) is roughly the solvent fraction
set rhomean = `echo | mapdump mapin solvent_avg.map | awk '/Mean density/{print $NF;exit}'`
set pct = `echo $rhomean | awk '{printf("%.0f",$1/0.334*100)}'`
echo "solvent map mean density: $rhomean e/A^3 ( ~ ${pct}% of the cell at bulk-water density )"
echo ""
echo "wrote $outfile ( FPART PHIPART, electron units, $smallSG unit cell, to $reso A )"
echo "phenix.refine model.pdb data.mtz refinement.input.bulk_solvent_map.file_name=$outfile \"
echo "   refinement.input.bulk_solvent_map.amplitudes_label=FPART refinement.input.bulk_solvent_map.phases_label=PHIPART"


exit:
if( ! $debug && "$t" =~ /dev/shm/* ) then
  rm -f ${t}* >& /dev/null
endif

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9
endif

exit
