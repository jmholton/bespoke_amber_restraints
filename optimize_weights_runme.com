#! /bin/tcsh -f
#
#  start with files:
#  ../refme.mtz
#  ../centroids/centroids_in_density.pdb
#  ${template_dir}/amberme.pdb
#  ../*.mol2 ../*.frcmod
#  tleap_stub.in
#  ${PHENIX}
#
#  iteratively optimize restraint weights based on fofc difference map
#
#

# defaults for all variables - can be overridden with file: settings.sourceme

# properties of the data in the mtz file
set reso = 1.85
set smallSG = P41212
# supercell parameters
set super_mult = 2,2,2
# number of residues in one protein
set modulo  = 272

# for map generation
set render_reso = 0.95
set render_B = 10
# iterative adjustment of overall B factor
set render_B_min = 10
set render_B_adjust = 0.05
# for gemmi
set render_rate = 1.5
# file for storing all atomic B factors
set Bfac_file = Bfac.pdb
set Bfac_maxmod = 2
set Bfac_modmode = add
# use rmsd variation of atoms to define B factor
#set Bfac_file = rmsd2B
set rmsd2B_range = 1-10
set rmsd2B_scale = 1

# useful when agreement is terrible
set scale_grid_search = 0

# place to look for starting files
set template_dir = ../template_dir/

# keep track of where we are
set thisdir = `dirname $0`
set path = ( $thisdir $path )
set tempfile = tempfile

set logfile = details.log
set itr = 0
set delete_traj = 0

# re-discover reference points every so often
set repick_itr = 30
set repick_itr_ramp = none
# remove water restraints that are too close
set filter_itr = 1
set filter_dist = 1.8
# periodically delete restraints that are historically too highly challenged
set release_itr = 50
# periodically force CA atoms with possible restraints to have them
set reinin_itr = 0
# periodically discard restraints that are not in density
set rhocheck_itr = 0
set min_ref_rho = 0.5
# periodically delete restraints to atoms that have alt lcs
set delete_altconf_itr = 0
# periodically just delete most-challenged restraints
set delete_badrest_itr = 0
# periodically reset all B factors
set Breset_itr = 0
# smooth the difference map
set fft_B = 15
set fft_B_ramp = none
set shan_B = auto
# criteria for statistical significance in difference peaks
set halfrho_pos = auto
set halfrho_neg = auto
set halfrho_ramp = none
# apply an overall scale factor to all weights every round
set weight_scale = 1
# apply a different scale factor to negative peaks
set weight_negscale = 0.95
# raise small weights to a power to drive them toward zero
set weight_power = 1.01
# criterion for a weak weight (multiples of kT)
set weight_power_kTmult = 1
# scale factor for converting "B factors" in restraint file to amber weights
set pdbscale = 0.01
# smallest amber weight to use
set min_weight = 0.0001
set min_weight_ramp = none
# apply a constant restraint weight to every atom, even if unspecified
set allatom_weight = 0
# starting value for newly created restraints
set weight0 = 5
# minimum weight to apply to ligands (
set min_lig_weight = 0.01
# weight to apply to CA atoms that have wandered from their reference points
set errant_CA_weight = 2
# upon release_itr, also release restraints on neighboring atoms
set release_radius = 1
# new weight to give to "released" restarints
set release_weight = 0.1
# maximum number of highly challenged restraints to release
set release_maxbad = 10
# also reset B factors of released atoms
set release_Breset = 0
set release_fftBreset = 0
# periodically delete randomly-selected restraints
set randel_itr = 0
set randel_fraction = 0.01
# default weight for re-discovered restraints
set repick_maxweight = 1
# maximum distance to look for nearest reference point
set repick_maxdist   = 6
set repick_maxdist_ramp = none
# scale distances to "HOH" in restraint file to allow protein restraints to "win"
set repick_hohscale = 2
set repick_hohscale_ramp = none
# amber restraints on flippable chiral centers
set chiral_weight = 10
set chiral_weight_ramp = none
# amber restraints on peptide bond dihedrals
set omega_weight = 50
set omega_weight_ramp = none
# periodically re-map water names
set remap_itr = 0
# periodically wrap non-restrained atoms to inside the supercell
set wrap_itr = 2
# periodically re-center the system on restrained atoms
set align_itr = 1
# periodically move water molecules from bad Fo-Fc density to good Fo-Fc density
set teleport_itr = 1
set teleport_waters = 10
set teleport_weight = 5
set teleport_mindist = 2.0
# periodically add or remove waters, depending on pressure and void size
set hydrate_itr = 1
set minvoid = 6
set hydrate = voids
set dehydrate = pressure
set pressure_scale = 0.1
# average pressure over previous runs
set pressure_avglast = 1
set void_scale = 1
# parametes for AddToBox
set add_radius = 2.4
# always make sure water count stays the same
set water_lock = 0
# run pressure measurement after normal amber run with all restraints removed
set barometer_cycles = 10000

# average electron density over this and previous amber runs
set avglast = 1
set avglast_ramp = none
# equilibrate the system for a bit before doing a produciton run
set equi_ns = 0.2
set equi_dt = 0.001
set equi_ns_ramp = none
# production run lenght in nanoseconds
set prod_ns = 0.5
set write_ps = 20
set dt = 0.002
set prod_ns_ramp = none
# simulation parameters
set temperature = 287
set temp_ramp = none
set thermostat = ""
set gamma_ln = 1.0
set barostat = 0
set kT = 0.6

# maximum allowed restraint weight
set max_weight = 999.99
# maximum scale factor for updating restraint weights
set max_mult = 2.0
set max_mult_ramp = none
# average restraint weights through bonds to neighboring atoms
set thrubond_avg_weight = 0
set thrubond_avg_weight_spread = 1
set thrubond_avg_weight_ramp = none
# make things like O1 and O2 on ASP have the same weight
set ambig_same_weight = 1

# trigger a release cycle only if weight is pegged at max value
set release_trigger = 0
set pegged_weight = 0
# trigger a randel cycle only after R factors have stabilized
set randel_trigger = 0

# keep memory of previous randel bins
set randel_bin = ""
set randel_last = 0

set nwaters = ""

set debug = 0

set scratch = /scratch/${USER}/opt_`hostname -s`_$$_
mkdir -p /scratch/${USER}

# commands to submit jobs to GPU or CPU cluster queues
set sruncpu = "srun --partition=refmac --exclude=crush18"
set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"


# read the command line to update variables and other settings
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
    else
      # no equal sign
    endif
    if("$key" == "debug") set debug = "1"
end

# shorthand for temporary file
set t = $tempfile

# find a place for scratch
foreach scratch ( $scratch /data/${USER}/scratch/ /scratch/${USER}/ /dev/shm/${USER}/ /tmp/${USER}/ )
  mkdir -p ${scratch} >& /dev/null
  if(-w ${scratch}) then
    echo "scratch = $scratch is writable"
    break
  endif
end
echo "scratch = $scratch"

echo "running as $0"


# start the ramps
if( "$fft_B_ramp" != "none" ) set fft_B = `echo $fft_B_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$min_weight_ramp" != "none" ) set min_weight = `echo $min_weight_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$repick_maxdist_ramp" != "none" ) set repick_maxdist = `echo $repick_maxdist_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$repick_hohscale_ramp" != "none" ) set repick_hohscale = `echo $repick_hohscale_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$avglast_ramp" != "none" ) set avglast = `echo $avglast_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$equi_ns_ramp" != "none" ) set equi_ns = `echo $equi_ns_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$prod_ns_ramp" != "none" ) set prod_ns = `echo $prod_ns_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$max_mult_ramp" != "none" ) set max_mult = `echo $max_mult_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$thrubond_avg_weight_ramp" != "none" ) set thrubond_avg_weight = `echo $thrubond_avg_weight_ramp | awk '$1~/^[0-9]/{print $1+0}'`

if( "$halfrho_ramp" != "none" ) set halfrho_pos = `echo $halfrho_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$halfrho_ramp" != "none" ) set halfrho_neg = `echo $halfrho_ramp | awk '$1~/^[0-9]/{print $1+0}'`

if( "$chiral_weight_ramp" != "none" ) set chiral_weight = `echo $chiral_weight_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$omega_weight_ramp" != "none" ) set omega_weight = `echo $omega_weight_ramp | awk '$1~/^[0-9]/{print $1+0}'`


touch $logfile

# find structure factor file
if(! -e refme.mtz) then
    ln -sf ../refme.mtz .

cat << EOF >! refmac_opts.txt
solvent no
blim 10 50
damp 0 0.5
weigh matrix 1
ncyc 5
make link Y
make hydr Y
make hout Y
EOF

# might want to add refinement step

endif

# in case we need to run tleap
if(! -e tleap_stub.in) then
  cp -p ${template_dir}/*.mol2 .
  cp -p ${template_dir}/*.frcmod .
#  cp -p ${PHENIX}/modules/amber_library/m/MSE.* .
  cp ${template_dir}/protonation.txt protonation.txt
  cp ${template_dir}/tleap_stub.in .
endif

if(! -e amberme.pdb) then

  cp ${template_dir}/amberme.pdb amberme.pdb
endif

# structure factors for best-phased reference map
if(! -e reference.mtz ) then
  cp  ../centroids/reference0.mtz reference.mtz
endif

# master list of points in space to use as restraint reference points
if(! -e all_possible_refpoints.pdb ) then
  set B0 = `echo $weight0 $pdbscale | awk '{print $1/$2}'`
  cat ../centroids/centroids_in_density.pdb |\
  awk -v B0=$B0 '! /^ATOM|^HETAT/{print;next}\
     {pre=substr($0,1,60);post=substr($0,67);rho=$NF;\
      B=B0*rho}\
     B>B0{B=B0}\
     {printf("%s%6.2f%s\n",pre,B,post)}' |\
  cat >! all_possible_refpoints.pdb
endif

if( $randel_itr ) then
  echo "generating randomized restraint deletion bins"
  echo $randel_fraction |\
  cat - all_possible_refpoints.pdb |\
  awk 'NR==1{f=$1;next}\
    ! /^ATOM|^HETAT/{next} {++i;\
    bin=int(rand()/f+1);B=bin/100;\
    printf("%s%6.2f%s\n",substr($0,1,60),B,substr($0,67))}' |\
  cat >! randB.pdb
  xsame_runme.com randB.pdb 
  cat xsame_restraints.pdb |\
  awk '{bin=substr($0,61,6)*100;id=substr($0,12,19);xyz=substr($0,31,24);\
     print bin,"|"xyz"|"id"| RANDEL"}' |\
  cat >! randel_selections.txt
endif



if(! -e restraints_for_${itr}.pdb) then
  # first try at restraints
  centroids_nearby_runme.com amberme.pdb reffile=all_possible_refpoints.pdb \
    softener=2 weight=Bfac maxdist=$repick_maxdist  \
    hohscale=$repick_hohscale \
    outfile=restraints_for_${itr}.pdb debug=$debug | tee c2r_${itr}.log

#  cat *_protonation.txt >! protonation.txt
   cp restraints_for_${itr}.pdb current_restraints.pdb

endif

# generate amber input files from "amberme.pdb", tleap_stub.in and any mol2 or frcmod files
if(! -e padded.parm7) then

  cp -p ${template_dir}/*.mol2 ${template_dir}/*.frcmod .
  if(! -e protonation.txt) cp ${template_dir}/protonation.txt .
  if(! -e tleap_stub.in) cp ${template_dir}/tleap_stub.in .

leap2amber.com amberme.pdb stages=Cpu,Min,Cool,Heat,Equi,EquiMin \
  protons=protonation.txt \
  refpoints=current_restraints.pdb \
  pdbscale=$pdbscale \
  gamma_ln=$gamma_ln barostat=$barostat \
  leapfile=tleap_stub.in padwater=50000 \
  cool_ns=0.01 heat_ns=0.2 equi_ns=$equi_ns prod_ns=$prod_ns \
  debug=$debug | tee leap2amber_${itr}.log
if( $status ) then
   set BAD = "leap2amber failed"
   goto exit
endif

set Stage = `awk '/^Stages /{print $(NF-1)}' leap2amber_${itr}.log`

endif


if(! -e ref.crd) then
   echo "generating new ref.crd "
   update_centroid_positions_runme.com current_restraints.pdb  outfile=ref.crd 
endif


if (! $?Stage ) then
  echo "no previous Stage defined"
  set itr = `ls -lLrt amber_*.rst7 |& egrep -v "equi|settle|unwrap" | awk -F "_" '/amber_/{print $NF+0}' | tail -n 1`
  if( "$itr" != "" ) set Stage = amber_${itr}
  #echo "no previous amber runs"
endif

if( ! $?Stage ) then
  set log = `ls -1Lrt *amber_* |& grep amber_ | grep -v equi | tail -n 1`
  set itr = `echo $log | awk -F "_" '{print $NF+0}' | tail -n 1`
  set Stage = `awk '/^Stages /{print $(NF-1)}' leap2amber_${itr}.log | tail -n 1`
  if( "$Stage" == "" ) set Stage = `awk '/^Stages /{print $(NF-1)}' leap2amber_*.log | tail -n 1`
endif


if(! -e orignames.pdb ) then 
    echo "cp ${template_dir}/orignames.pdb ."
    cp ${template_dir}/orignames.pdb .
endif
if(! -e Bfac.pdb ) then 
   echo "generating Bfac.pdb from orignames.pdb"
   cat orignames.pdb |\
   awk '! /^ATOM|^HETAT/{print;next} \
     {printf("%s%6.2f%s\n",substr($0,1,60),10,substr($0,67))}' orignames.pdb |\
   cat >! Bfac.pdb 
endif
set norig = `egrep "^ATOM|^HETAT" orignames.pdb | wc -l`
set nBfac = `egrep "^ATOM|^HETAT" Bfac.pdb | wc -l`
if( $nBfac < $norig ) then
  echo "lengthening Bfac.pdb to match orignames.pdb"
  combine_pdbs_runme.com Bfac.pdb printref=1 orignames.pdb >! Bfac_edits.log
  mv new.pdb Bfac.pdb
endif
if( $nBfac > $norig ) then
  echo "shortening Bfac.pdb to match orignames.pdb"
  combine_pdbs_runme.com Bfac.pdb orignames.pdb  >! Bfac_edits.log
  mv new.pdb Bfac.pdb
endif

echo "taking restraints from itr $itr"
cp restraints_for_${itr}.pdb current_restraints.pdb
echo "current Stage = $Stage"

if("$Stage" == "") then
  set BAD = "cannont determine Stage"
  goto exit
endif

# test and find a topfile that works with the current stage
echo "${Stage}.rst7 -> this.pdb"
rst2pdb_runme.com ${Stage}.rst7 outprefix=this >! rst2pdb_${Stage}.log
if( $status ) then
  set BAD = "no top files work with ${Stage}.rst7"
  goto exit
endif
set newparm = `awk '/new parmfile:/{print $3}' rst2pdb_${Stage}.log | tail -n 1`
if( -e "$newparm" ) then
  echo "new xtal.prmtop"
  mv $newparm xtal.prmtop
  echo "prev: "
  tail -n 4 orignames.pdb | egrep "^ATOM|^HETAT" | tail -n 1 
  echo "now:  "
  tail -n 4 this.pdb | egrep "^ATOM|^HETAT" | tail -n 1 
#  echo "new orignames.pdb"
#  mv new_orignames.pdb orignames.pdb
  update_centroid_positions_runme.com current_restraints.pdb outfile=ref.crd orignames=this.pdb |\
  tee ucp_0.log | awk 'NF==3 && $2=="="{next} {print}'
endif

if(! -e all_possible_refpoints.pdb) then
  echo "generating all_possible_refpoints.pdb with weight0 = $weight0"
  set B0 = `echo $weight0 $pdbscale | awk '{print $1/$2}'`
  cat ../centroids/centroids_in_density.pdb |\
  awk -v B0=$B0 '! /^ATOM|^HETAT/{print;next}\
      {pre=substr($0,1,60);post=substr($0,67);rho=$NF;\
      B=B0*rho}\
      B>B0{B=B0}\
      {printf("%s%6.2f%s\n",pre,B,post)}' |\
  cat >! all_possible_refpoints.pdb
endif

echo "checking for restraint bombs"
rm -f defused_restraints.pdb
restraint_bomb_detector.com ${Stage}.rst7 \
    refpoints=current_restraints.pdb \
    outfile=defused_restraints.pdb \
    outcrd=ref.crd >! restraint_bomb_check.log
if( $status ) then
  set BAD = "unfixable restraints bombs detected"
  goto exit
endif
egrep "WARNING|worst|ATOM|HETAT" restraint_bomb_check.log
egrep -v "^recommend:|^making ref.crd|^mv " restraint_bomb_check.log |\
tail -n 1 
diff current_restraints.pdb defused_restraints.pdb > /dev/null
if( $status ) then
  cat restraint_bomb_check.log
  echo "using defused_restraints.pdb"
  cp defused_restraints.pdb current_restraints.pdb
endif


while ( 1 )

# allow run-time edits to variables
if(-e settings.sourceme) then
    cat settings.sourceme
    source settings.sourceme
endif
rm -f ./exit >& /dev/null

set prev = "$itr"
@ itr = ( $itr + 1 )
set laStage = $Stage
set Stage = amber_${itr}


# see if we need a wrap - unprintable atoms
echo "extracting xyz from ${laStage}.rst7 -> this.pdb"
rst2pdb_runme.com ${laStage}.rst7 outprefix=this >! rst2pdb.log
cat this.pdb |\
awk '! /^ATOM|^HETAT/{next}\
  /*/{print;exit}\
  {X=substr($0,31,8)+0;Y=substr($0,39,8)+0;Z=substr($0,47,8)+0}\
  X<-800 || X>800 || Y<-800 || Y>800 || Z<-800 || Z>800{print;exit}' |\
cat >! temp.txt
set need_wrap_now = `cat  temp.txt | wc -l`
if( $need_wrap_now ) then
   echo "WARNING: atoms near edge of pritable box:"
   cat temp.txt
endif

if ( $repick_itr && ( $itr % $repick_itr == 0 ) ) set need_wrap_now = 1

if ( $need_wrap_now || $wrap_itr && ( $itr % $wrap_itr == 0 ) ) then
  echo "wrapping ${laStage}.rst7"

  # defaults to "all_possible_refpoints.pdb and current_restraints.pdb"
  rst_wrap.com ${laStage}.rst7 | tee wrap_${itr}.log | egrep "rmsd_wrap|maxd_wrap"

  cp ${laStage}.rst7 ${laStage}_unwrapped.rst7
  mv wrapped.rst7 ${laStage}.rst7
  if(-e wraps.txt) mv wraps.txt wraps_${itr}.txt
  # keep everything else same for laStage

  # update the this.pdb
  rst2pdb_runme.com ${laStage}.rst7 outprefix=this >! rst2pdb.log

#  if(-e ${laStage}.prmtop) cp ${laStage}.prmtop wrapped.prmtop
#  set laStage = wrapped

  echo "wrapping any errant refpoints..."
  rm -f defused_restraints.pdb
  restraint_bomb_detector.com ${laStage}.rst7 \
    refpoints=current_restraints.pdb \
    outfile=defused_restraints.pdb \
    outcrd=ref.crd >! restraint_bomb_check.log
  if( $status ) then
    set BAD = "restraints moved relative to wrapped atoms"
    goto exit
  endif
  egrep "WARNING|worst|ATOM|HETAT" restraint_bomb_check.log
  egrep -v "^recommend:|^making ref.crd|^mv " restraint_bomb_check.log |\
  tail -n 1
  diff current_restraints.pdb defused_restraints.pdb > /dev/null
  if( $status ) then
    cat restraint_bomb_check.log
    echo "using defused_restraints.pdb"
    cp defused_restraints.pdb current_restraints.pdb
  endif
endif

cat this.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  substr($0,77,2)!="XP" && substr($0,12,5)!=" EPW "{print}' |\
cat >! refme.pdb
#cp refme.pdb refme0.pdb


# check for inverted chiral centers and cis peptides
echo "checking chiral centers and peptide bonds"
awk '! /HOH/{print substr($0,1,80)}' refme.pdb >! protein.pdb
( phenix.omegalyze protein.pdb >! omegalyze_${prev}.log ) >& /dev/null &
( phenix.chiral_validation protein.pdb >! chiralyze_${prev}.log ) >& /dev/null &


# report on last itr teleportation results
if(-e teleported_waters_${prev}.pdb) then
  echo "previous teleported waters from run $prev are now:"
  awk '/^ATOM|^HETAT/{print substr($0,12,17)"| SEL"}' teleported_waters_${prev}.pdb |\
  cat - refme.pdb |\
  awk '$NF=="SEL"{id=substr($0,1,17);++sel[id];next}\
     {id=substr($0,12,17)}\
     sel[id]{print}' >! new_positions_${prev}.pdb

  rmsd -v debug=1 teleported_waters_${prev}.pdb new_positions_${prev}.pdb |\
    grep moved | grep "O   HOH" | sort -k1.25g | tail 
endif

if(0) then
echo "combining xyz with B factors from Bfac.pdb in refme.pdb"
egrep "^SSBOND|^LINK|^CISP|^CRYST" Bfac.pdb >! refme.pdb
# take coordinates only, using names from starting point
awk '/^ATOM|^HETAT/{print $0,"ORIG"}' Bfac.pdb |\
cat - this.pdb |\
awk '$NF=="ORIG"{++o;pre[o]=substr($0,1,30);post[o]=substr($0,55,length($0)-55-4);next}\
  ! /^ATOM|^HETAT/{next}\
      {++n;resid=substr($0,22,9)}\
      lastres!=resid{lastres=resid;++ordresnum}\
      pre[n]==""{print "REMARK WARNING atom",n,"missing from orig, max="o;\
       pre[n]=substr($0,1,30);post[n]=substr($0,55)}\
      {printf("%s%s%s\n",pre[n],substr($0,31,24),post[n])}' |\
egrep -v "EPW" >> refme.pdb
endif


if( "$fft_B_ramp" != "none" ) then
    set params = `echo $fft_B_ramp | awk -F "[: _-]" '{print $1,$2}'`
    set value = `echo $params $fft_B | awk '$NF<=0{print 0;exit} /x/{print $NF*$2;exit} $1>$2 && $NF>$2{print $NF-$2;exit} {print $NF-$1}'`
    set fft_B = $value
    echo "fft_B = $fft_B"
endif
if( "$avglast_ramp" != "none" ) then
    set params = `echo $avglast_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $avglast | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 1 | awk '$1<$2{$1=$2} {print $1}'`
    set avglast = $value
    echo "avglast = $avglast_avgB"
endif

if( "$repick_maxdist_ramp" != "none" ) then
    set params = `echo $repick_maxdist_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $repick_maxdist | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=$2} {print $1}'`
    set repick_maxdist = $value
    echo "repick_maxdist = $repick_maxdist"
endif

if( "$repick_hohscale_ramp" != "none" ) then
    set params = `echo $repick_hohscale_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $repick_hohscale | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=$2} {print $1}'`
    set repick_hohscale = $value
    echo "repick_hohscale = $repick_hohscale"
endif

if( "$max_mult_ramp" != "none" ) then
    set params = `echo $max_mult_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $max_mult | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=$2} {print $1}'`
    set max_mult = $value
    echo "max_mult = $max_mult"
endif

if( "$thrubond_avg_weight_ramp" != "none" ) then
    set params = `echo $thrubond_avg_weight_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $thrubond_avg_weight | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=$2} {print $1}'`
    set thrubond_avg_weight = $value
    echo "thrubond_avg_weight = $thrubond_avg_weight"
endif

if( "$min_weight_ramp" != "none" ) then
    set params = `echo $min_weight_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $min_weight | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.0001 | awk '$1<$2{$1=$2} {print $1}'`
    set min_weight = $value
    echo "min_weight = $min_weight"
endif

if( "$halfrho_ramp" != "none" ) then
    set params = `echo $halfrho_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $halfrho_pos | awk 'NF>2 && $NF>=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=$2} {print $1}'`
    set halfrho_pos = $value
    set halfrho_neg = $value
    echo "halfrho = $halfrho_pos"
endif

# set omega_weight_ramp = 50-0:0.5x
if( "$omega_weight_ramp" != "none" ) then
    set params = `echo $omega_weight_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $omega_weight | awk 'NF>2 && $3>$2 && $NF<=$2{print $2;exit} NF>2 && $3<$2 && $NF>=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=0} {print $1}'`
    set omega_weight = $value
    echo "omega_weight = $omega_weight"
endif

if( "$chiral_weight_ramp" != "none" ) then
    set params = `echo $chiral_weight_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $chiral_weight | awk 'NF>2 && $3>$2 && $NF<=$2{print $2;exit} NF>2 && $3<$2 && $NF>=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=0} {print $1}'`
    set chiral_weight = $value
    echo "chiral_weight = $chiral_weight"
endif


if ( $align_itr && ( $itr % $align_itr == 0 ) ) then
  # align rst7 and nc files to reference
  echo "aligning ${laStage} to current restraints"
  touch align_${itr}.log
  #rst2pdb_runme.com ${laStage}.rst7 this.pdb Bfactors=none >> align_${itr}.log
  #if(-e resized.parm7) mv resized.parm7 xtal.prmtop
  #update_centroid_positions_runme_new.com current_restraints.pdb orignames=this.pdb topfile=xtal.prmtop >> align_${itr}.log
  # mv newref.crd ref.crd
  cat current_restraints.pdb |\
  awk 'substr($0,61,6)+0>0{print substr($0,1,80),"     SEL"}' |\
  cat >! labeled.pdb
  combine_pdbs_runme.com labeled.pdb this.pdb printref=1 >> align_${itr}.log
  egrep "^ATOM|^HETAT" new.pdb |\
  awk -v pdbscale=$pdbscale '{++n} $NF=="SEL"{print n,substr($0,61,6)*pdbscale}' |\
  cat >! restrained_atom_weights.txt
  foreach w ( 1 med )
    if( $w == med ) then
      set h = `cat restrained_atom_weights.txt | wc -l | awk '{print int($1/2)}'`
      set w = `awk '{print $2}' restrained_atom_weights.txt | sort -g | head -n $h | tail -n 1`
    endif
    awk -v w=$w '$2>w' restrained_atom_weights.txt |\
    awk 'NR==1{s=e=$1;next}\
         $1==e+1{e=$1;next} \
         {print s"-"e;s=e=$1}\
         END{print s"-"e}' |\
    awk -F "-" '$1==$2{print $1;next} {print}' >! rest_ranges.txt 
    set rest_ranges = `cat rest_ranges.txt `
    set rest_mask = `echo $rest_ranges | awk '{gsub(" ",",");print}'`
    if( "$rest_mask" != "" ) then
      break
    endif
  end
  if( "$rest_mask" == "" ) then
    set BAD = "unable to assing restraints"
    goto exit
  endif

  cpptraj -p xtal.prmtop -y ${laStage}.rst7 -c ref.crd << EOF >> align_${itr}.log
  rmsd rmsd reference norotate @$rest_mask out rmsd.txt savevectors combined vecsout vecsout.txt
  trajout aligned.rst7
EOF
  cat rmsd.txt vecsout.txt >> align_${itr}.log

  cpptraj -p xtal.prmtop -y ${laStage}.nc -c ref.crd << EOF >> align_${itr}.log
  rmsd rmsd reference norotate @$rest_mask out rmsd.txt savevectors combined vecsout vecsout.txt
  trajout aligned.nc
EOF
  cat rmsd.txt vecsout.txt >> align_${itr}.log

  cp ${laStage}.in aligned.in
  set laStage = aligned

  echo -n "final shift: "
  tail -n 1 vecsout.txt | awk '{print $2,$3,$4}'

  # update the this and noEP pdb files
  rst2pdb_runme.com ${laStage}.rst7 outprefix=this >! rst2pdb.log
  cat this.pdb |\
  awk '! /^ATOM|^HETAT/{print;next}\
    substr($0,77,2)!="XP" && substr($0,12,5)!=" EPW "{print}' |\
  cat >! refme.pdb

endif



# convert nc file to an mtz
set nc2log = nc2mtz_${prev}.log
if(-e ${laStage}.nc && ( ! -e avg_${prev}.mtz || ! -e trajectory/md.1.pdb ) ) then
  set minB = 1
  set maxB = 999.99
  if( "$Bfac_file" == "rmsd2B" ) then
    set minB = `echo $rmsd2B_range | awk -F "-" '{print $1}'`
    set maxB = `echo $rmsd2B_range | awk -F "-" '{print $2}'`
  endif
  set nc2mtz_extraopt = ""
  if( -e avg_${prev}.mtz ) set nc2mtz_extraopt = "domaps=0"
  echo "nc2mtz gemmi ${laStage}.nc"
  nc2mtz_gemmi.com $smallSG super_mult=$super_mult ${laStage}.nc \
    reso=$render_reso B=$render_B \
    Bfac_file=$Bfac_file minB=$minB maxB=$maxB \
    keeptraj=1 wrap=1 rate=$render_rate \
    tempfile=${scratch}/nc2mtz_$$_ $nc2mtz_extraopt >>& $nc2log
  if($status) then
    set BAD = "nc2mtz gemmi failed"
    goto exit
  endif
  if("$nc2mtz_extraopt" == "") cp avg.mtz avg_${prev}.mtz
  
  # dont forget to delete the expanded trajectory when we are done with it
  set delete_traj = 1
else
  echo "no need for nc2mtz"
endif

echo "averaging last ${avglast} itrs"
ls -1rt avg_*.mtz | awk -F "[_.]" '$2!~/[a-z]/' | tail -n ${avglast} >! latest_avgs.txt
set avg_these = `cat latest_avgs.txt`
echo "$#avg_these mtz files to average"
if( ! $#avg_these ) goto skipfofc

addup_mtzs_runme.com $avg_these  >>& $nc2log
if($status) then
  set BAD = "addup_mtzs failed"
  goto exit
endif

set scale = `echo $#avg_these | awk '{print 1/$1}'`
cad hklin1 sum.mtz hklout avg_lastN.mtz << EOF >> $nc2log
labin file 1 E1=Fsum E2=PHIsum
scale file 1 $scale
labou file 1 E1=FCavg E2=PHICavg
EOF


if( $scale_grid_search ) then
echo "scaling grid search..."
scaleB_search_diffmap_runme.com avg_lastN.mtz |& tee sBs_${itr}.log | grep best
set sgstats = `awk '/best-fit scale,B/{s=$4;B=$5} /^TOTAL/{CC=$NF} /correct F:/{R=$7} END{print R,s,B,CC}' sBs_${itr}.log`
echo "$itr $sgstats" | tee -a fofc_Rplot_grid.txt
# generates cootme-grid.mtz
#cp cootme-grid.mtz cootme.mtz
else
  set sgstats = ( 99 1 0 0 )
endif

echo "generating normalized fofc map - scaleit"
diff.com reference.mtz avg_lastN.mtz >> $nc2log

cad hklin1 Fdiff.mtz \
    hklin2 avg_lastN.mtz \
    hklin3 reference.mtz \
    hklout fftme.mtz << EOF >> $nc2log
labin file 1 E1=Fref E2=Ftest
labin file 2 E1=PHICavg
labou file 2 E1=PHItest
labin file 3 E1=PHIref
resolution over_all $reso
EOF

rm -f cootme-scaleit.mtz
sftools << EOF >> $nc2log
read fftme.mtz
map correl Fref PHIref Ftest PHItest
calc ( COL DELFWT PHDELWT ) = ( COL Fref PHIref ) ( COL Ftest PHItest ) -
calc ( COL FWT PHWT ) = ( COL Fref PHIref )
write cootme-scaleit.mtz col FWT PHWT DELFWT PHDELWT
quit
EOF
set sitstats = `awk '/scale=/{s=$2;B=$4} /^TOTAL/{CC=$NF} /correct F:/{R=$7} END{print R,s,B,CC}' $nc2log`
echo "$itr $sitstats" | tee -a fofc_Rplot.txt

set stats = "$sitstats"
cp cootme-scaleit.mtz cootme.mtz

set test = `echo $sgstats[1] $sitstats[1] | awk '{print ($1<$2)}'`
if( $test ) then
  echo "using grid-search scaling results"
  set stats = "$sgstats"
  cp cootme-grid.mtz cootme.mtz
else
  echo "using scaleit results"
  set stats = "$sitstats"
  cp cootme-scaleit.mtz cootme.mtz
endif


if( "$render_B_adjust" != "0" ) then
  if("$Bfac_file" == "rmsd2B" ) then
     echo "adjusting sfall rendering B factor"
     set render_B = `echo $render_B $render_B_adjust $render_B_min $stats | awk '{newB=$1+$2*$6} newB<$3{newB=$3} {print newB}'`
     echo "next render_B will be: $render_B"
  endif
  if(-e "$Bfac_file") then
     set deltaB = `echo $stats $render_B_adjust | awk '{print $3*$NF}'`
     echo "adding $deltaB to all B factors in $Bfac_file"
     echo "$deltaB 2 999.99" |\
     cat - "$Bfac_file" |\
     awk 'NR==1{deltaB=$1;minB=$2;maxB=$3;next}\
      ! /^ATOM|^HETAT/{print;next}\
      {pre=substr($0,1,60);post=substr($0,67);B=substr($0,61,6);\
       newB=B+deltaB;}\
      newB<minB{newB=minB}\
      newB>maxB{newB=maxB}\
      {printf("%s%6.2f%s\n",pre,newB,post)}' |\
     cat >! newB.pdb
     cp ${Bfac_file} Bfac_preadjust_${itr}.pdb
     cp newB.pdb ${Bfac_file}
  endif
endif


if(! -e reference.map) then
echo "creating reference.map"
fft hklin reference.mtz mapout ffted.map << EOF >> $nc2log
labin F1=Fref PHI=PHIref
reso $reso
EOF
mapmask mapin ffted.map mapout reference.map << EOF >> $nc2log
xyzlim asu
EOF
endif




# make the Fo-Fc map for restraint updates
if(-e fofc.map) mv fofc.map prev_fofc.map
fft hklin cootme.mtz mapout ffted.map << EOF >> $nc2log
labin F1=DELFWT PHI=PHDELWT
scale F1 1 $fft_B
reso $reso
EOF
mapmask mapin ffted.map mapout fofc.map << EOF >> $nc2log
scale sigma 1
xyzlim asu
EOF
rm -f ffted.map
if(-e prev_fofc.map) then
   echo -n "CC to last fofc map: "
   correlate.com fofc.map prev_fofc.map | awk '{print $NF}'
endif

if(! -e trajectory/md.1.pdb ) then
  echo "ERROR: need trajectory."
endif

# use difference map to update restraints
echo "adjusting restraints based on fofc map"
restraintlist_update_diffmap.com fofc.map \
  refmap=reference.map \
  trajectory=trajectory/ \
  tempfile=${scratch}/rud_ \
  refme.pdb modulo=$modulo \
  overall_scale=$weight_scale \
  negative_scale=$weight_negscale \
  max_weight=$max_weight max_mult=$max_mult \
  ambig_same_weight=$ambig_same_weight \
  halfrho_pos=$halfrho_pos halfrho_neg=$halfrho_neg \
  refpointspdb=current_restraints.pdb \
  outmults=sorted_mults_${itr}.txt \
  outfile=new_restraints.pdb | tee restraint_update_${itr}.log |\
  awk 'NF<=2 || NF==3 && $2=="="{next} /^srun/{next} {print}'
if( $status ) then
  set BAD = "restraintlist_update_diffmap failed"
  goto exit
endif

cat current_restraints.pdb new_restraints.pdb |\
awk '/^ATOM|^HETAT/{id=substr($0,12,19);B=substr($0,61,6)+0}\
  lastB[id]+0>0{print id,B/lastB[id],"ACTUAL"} {lastB[id]=B}' |\
cat >! actual_mults.txt

cp new_restraints.pdb current_restraints.pdb


# now update B factors for map calculation
if(-e premod_Bfac_${itr}.pdb) then
   echo "assuming B factors in premod_Bfac_${itr}.pdb are for amber_${itr}.nc"
   cp premod_Bfac_${itr}.pdb Bfac.pdb
endif
echo "updating B factors based on difference map"
Bfac_update_diffmap.com Bfac.pdb modulo=$modulo mtzfile=cootme.mtz \
  max_mod=$Bfac_maxmod mod_mode=$Bfac_modmode \
  tempfile=${scratch}/Bud_ |&\
  tee Bfac_update_${itr}.log | awk 'NF<=2 || NF==3 && $2=="=" || /^\[/ || /\/dev\/shm/{next} {print}'

cp Bfac.pdb premod_Bfac_${itr}.pdb
cp new_Bfac.pdb Bfac.pdb
cp Bfac.pdb Bfac_${itr}.pdb
cp sorted_Bmods.txt sorted_Bmods_${itr}.txt

# re-generate this.pdb and refme.pdb with new B factors
rst2pdb_runme.com ${laStage}.rst7 outprefix=this >! rst2pdb.log
cat this.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  substr($0,77,2)!="XP" && substr($0,12,5)!=" EPW "{print}' |\
cat >! refme.pdb


if( $delete_traj && -e ./trajectory ) then
  # dont forget to clean up
  set deldir = `ls -ld ./trajectory | awk '{print $NF}'`
  if( "$deldir" != "" && -e "$deldir" ) then
    echo "deleting $deldir"
    rm -rf $deldir
  endif
  rm trajectory
endif


skipfofc:

echo "measuring challenges to previous restraints"
cat this.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  substr($0,77,2)!="XP" && substr($0,12,5)!=" EPW "{print}' |\
tee refme.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  substr($0,77,2)!=" H"{print substr($0,1,60) "  0.00" substr($0,67)}' |\
cat >! zeroB.pdb

# find restraints bigger than a kT
echo $min_weight $pdbscale $kT |\
cat - current_restraints.pdb |\
awk 'NR==1{minB=$1/$2;scale=$2;kT=$3;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {B=substr($0,61,6)+0}\
  #w=B*scale;quadsub=w**2-kT**2}\
  #quadsub<0{next}\
  {pre=substr($0,1,60);post=substr($0,67)}\
  B>minB{printf("%s%6.2f%s\n",pre,B,post)}' |\
cat >! active_restraints.pdb

# look for changes in restraints
rmsd -v debug=1 active_restraints.pdb zeroB.pdb |\
awk -v scale=$pdbscale '/moved/{print substr($0,25,9),substr($0,53,9)*scale,substr($0,1,17)}' |\
awk '{print $1*$1*$2,$0}' |\
sort -g >! restraints_challenge_vs_weight.txt
cp restraints_challenge_vs_weight.txt restraints_challenge_vs_weight_${prev}.txt
echo -n "worst challenge: "
echo -n "$prev " | tee -a worstchallenge_vs_itr.txt
awk '$4>0.5' restraints_challenge_vs_weight_${prev}.txt |\
tail -n 1 | tee -a worstchallenge_vs_itr.txt










echo "measuring voids..."
set maxvoid = `measure_voids_runme.com zeroB.pdb | awk '/biggest void:/{print $3/92}'`

# also use gemmi
gemmi mask refme.pdb void.msk
gemmi map2sf void.msk void.mtz Fbulk PHIbulk
gemmi blobs void.mtz -f Fbulk -p PHIbulk refme.pdb  >! gemmi_blobs.txt
# voids only count if more than 3 waters fit
awk '$3=="el"{print $5/92}' gemmi_blobs.txt |\
awk '$1>1' >! void_hist.txt
set nbulk_sum = `awk -v minvoid=$minvoid '$1>=minvoid{sum+=$1} END{print sum+0}' void_hist.txt`
set nbulk_max = `awk 'NR==1{print $1}' void_hist.txt`
if("$nbulk_max" == "") set nbulk_max = 0

#dont forget to update this.pdb if the rst7 file gets shuffled


# check the pressure, if any
foreach log ( barometer_${prev}.out ${laStage}.out amber_${prev}.out )
  if(! -e "$log") continue
  set press = `awk '/PRESS/ && p{n=$3;print $NF} /A V E R A/{++p} END{print P,r,n+0}' $log`
  if( "$press" == "" || "$press" == "0" ) set press = ( 0 - - )
  if( "$press[1]" != "0" ) break
end
if( "$press[1]" == "0" ) then
   set ektot = `awk '/EKtot/ && ! p{++n} /EKtot/ && p{print $6} /A V E R A/{++p} END{print n+0}' ${laStage}.out`
endif

if( "$nwaters" == "" ) then
  set nwaters = `echo list | cpptraj -p xtal.prmtop | awk '$4=="atoms,"{print $11, $5, $3}' | head -n 1`
endif
echo "$itr pressure was $press   void: $maxvoid $nbulk_sum $nbulk_max   nwater: $nwaters" |\
  tee -a pressure_vs_itr.txt

if( $pressure_avglast > 1 ) then
  set press = `tail -n $pressure_avglast pressure_vs_itr.txt | awk '{++n;sum+=$4;vsum+=$5*$5} END{print sum/n,sqrt(vsum/n),$6}'`
  echo "using pressure: $press  averaged over $pressure_avglast itrs"
endif

set changed_waters = 0
set dehydrate_thistime = 0
set nadd = 0
if ( $hydrate_itr && $itr % $hydrate_itr == 0  ) then
  echo "checking pressure: $press[1] > 10 "
  set test = `echo $press | awk '{print ( $1>10 )}'`
  if( $test && "$dehydrate" != "0" ) then
    if( "$dehydrate" == "pressure" ) then
      set dehydrate_thistime = `echo $press $pressure_scale | awk '{print int($1*$NF)}'`
      if( $dehydrate_thistime == 0 ) set dehydrate_thistime = 1
    else
      set dehydrate_thistime = "$dehydrate"
    endif
  endif

  echo "checking void size: $nbulk_max vs $minvoid"
  set bigvoid = `echo $nbulk_max $minvoid | awk '{print ( $1 > 2*$2 )}'`
  set lowpress = `echo $press -10 | awk '{print ( $1 < $2 )}'`
  set addv = `echo $nbulk_max $void_scale | awk '{print int($1*$2)}'`
  set addp = `echo $press $pressure_scale | awk '{print int(-$1*$NF)*($1<0)}'`
  echo "nadd options: $addv $addp"
  if( $bigvoid || $lowpress ) then
    echo "need to add water... "
    set nadd = `echo $addp $addv | awk '$1<$2{$1=$2} {print $1}'`
  endif

  if( $water_lock ) then
    # keep water count constant
    set dwater = `echo $nadd $dehydrate_thistime | awk '{print int(($1+$2)/2)}'`
    echo "average delta-water: $dwater"
    set nadd = $dwater
    set dehydrate_thistime = $dwater
    set revertStage = $laStage
  endif

  if( $dehydrate_thistime ) then
    echo "dropping ${dehydrate_thistime} waters"
    set changed_waters = 1

    remap_waters_runme.com ${laStage}.rst7 current_restraints.pdb | tee dehydrate_${itr}.log
    if( $status ) then
      set BAD = "remap failed"
      goto exit
    endif
    cp remapped_restraints.pdb current_restraints.pdb
    egrep -v "^REMARK" remapped_Bfac.pdb >! Bfac.pdb
    cp refme.pdb preremap_${itr}.pdb
    cp remapped_Bfac.pdb remapped_${itr}.pdb
    cp remapped_pairs.txt remapped_pairs_${itr}.txt

    # clear previous teleport records, as they are now scrambled
#    if(-e teleported_waters_${prev}.pdb) then
#       echo "previous teleport locations no longer interesting"
#       mv teleported_waters_${prev}.pdb old_teleported_waters_${prev}.pdb
#    endif

    dehydrate_amber_runme.com xtal.prmtop remapped.rst7 maxreject=${dehydrate_thistime} | tee -a dehydrate_${itr}.log

    echo "updating xtal.prmtop"
    mv drier.parm7 xtal.prmtop
#    mv orignames.pdb prev_orignames.pdb
#    echo "updating orignames.pdb"
#    mv drier_orignames.pdb orignames.pdb
    echo "updating Bfac.pdb"
    cp -p Bfac.pdb predry_Bfac_${itr}.pdb
    cp -p drier_Bfac.pdb Bfac.pdb
    cp ${laStage}.in drier.in
    set laStage = drier

    echo "re-generating ref.crd from current_restraints.pdb"
    update_centroid_positions_runme.com current_restraints.pdb  outfile=ref.crd >> dehydrate_${itr}.log
    if( $status ) then
       set BAD = "restraint update fail"
       goto exit
    endif

    echo "checking for errant refpoints..."
    rm -f defused_restraints.pdb
    restraint_bomb_detector.com ${laStage}.rst7 \
      refpoints=current_restraints.pdb kT=1.2 \
      outfile=defused_restraints.pdb \
      outcrd=ref.crd >! restraint_bomb_check.log
    if( $status ) then
      set BAD = "restraint bomb detected after dehydration"
      goto exit
    endif
    egrep "WARNING|worst|ATOM|HETAT" restraint_bomb_check.log
    egrep -v "^recommend:|^making ref.crd|^mv " restraint_bomb_check.log |\
    tail -n 1 
    diff current_restraints.pdb defused_restraints.pdb > /dev/null
    if( $status ) then
      cat restraint_bomb_check.log
      echo "using defused_restraints.pdb"
      cp defused_restraints.pdb current_restraints.pdb
    endif
  else
    echo "no need to delete waters"
  endif


  if( $nadd ) then
    echo "adding $nadd waters"

    if(-e padded.parm7) then
     set changed_waters = 1
     hydrate_runme.com ${laStage}.rst7 outrst=wetter.rst7 \
        nadd=$nadd water_radius=$add_radius force_nadd=$water_lock debug=$debug \
        restraints=current_restraints.pdb | tee hydrate_${itr}.log
      # interanlly updates xtal.prmtop, ref.crd 
      # never ever change orignames.pdb

      cp ${laStage}.in wetter.in
      set laStage = wetter

      set nadded = `awk '/^added/{print $2;exit}' hydrate_${itr}.log`
      if( "$nadded" != "$nadd" && $water_lock ) then
        echo "WARNING: failed to add $nadd waters, reverting to $revertStage"
        set laStage = $revertStage
        rst2pdb_runme.com ${laStage}.rst7 temp.pdb
        mv resized.parm7 xtal.prmtop
        update_centroid_positions_runme.com current_restraints.pdb  outfile=ref.crd >> hydrate_${itr}.log
      endif

      echo "checking for errant refpoints..."
      rm -f defused_restraints.pdb
      restraint_bomb_detector.com ${laStage}.rst7 \
        refpoints=current_restraints.pdb \
        outfile=defused_restraints.pdb \
        outcrd=ref.crd >! restraint_bomb_check.log
      if( $status ) then
        set BAD = "restraint bomb detected after hydration"
        goto exit
      endif
      egrep "WARNING|worst|ATOM|HETTAT" restraint_bomb_check.log
      egrep -v "^recommend:|^making ref.crd|^mv " restraint_bomb_check.log |\
      tail -n 1 
      diff current_restraints.pdb defused_restraints.pdb > /dev/null
      if( $status ) then
        cat restraint_bomb_check.log
        echo "using defused_restraints.pdb"
        cp defused_restraints.pdb current_restraints.pdb
      endif
    else
      echo "WARNING cannot add water without padded.parm7"
    endif
  else
    echo "no need to add more water"
  endif
endif

if( $changed_waters ) then
  # add or drop waters in the Bfac.pdb file too
  set norig = `egrep "^ATOM|^HETAT" orignames.pdb | wc -l`
  set nBfac = `egrep "^ATOM|^HETAT" Bfac.pdb | wc -l`
  if( $nBfac < $norig ) then
    echo "padding Bfac.pdb to same length as orignames.pdb"
    combine_pdbs_runme.com Bfac.pdb printref=1 orignames.pdb > /dev/null
    cp Bfac.pdb prepad_Bfac_${itr}.pdb
    mv new.pdb Bfac.pdb
  endif
  if( $nBfac > $norig ) then
    echo "truncating Bfac.pdb to same length as orignames.pdb"
    combine_pdbs_runme.com Bfac.pdb orignames.pdb  > /dev/null
    cp Bfac.pdb pretrunc_Bfac_${itr}.pdb
    mv new.pdb Bfac.pdb
  endif

  echo -n "number of waters: "
  set nwaters = `echo list | cpptraj -p xtal.prmtop | awk '$4=="atoms,"{print $11, $5, $3}' | head -n 1`
  echo "$itr $nwaters" | tee -a nwaters_vs_itr.txt

  # update the this and noEP pdb files
  rst2pdb_runme.com ${laStage}.rst7 outprefix=this >! rst2pdb.log
  cat this.pdb |\
  awk '! /^ATOM|^HETAT/{print;next}\
    substr($0,77,2)!="XP" && substr($0,12,5)!=" EPW "{print}' |\
  cat >! refme.pdb
endif

set repick_now = 0

if ( $teleport_waters && -e cootme.mtz && $teleport_itr && ( $itr % $teleport_itr == 0  ) ) then
  echo "looking for teleportable waters in ${laStage}.rst7"
else
  goto skipteleport
endif

cat restraints_challenge_vs_weight.txt |\
awk '$2<2.4{print substr($0,index($0,$5)-6),"NOTME"}' |\
cat - current_restraints.pdb |\
awk '/NOTME/{++notme[substr($0,1,17)];next}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,12,17)}\
  notme[id]{print substr($0,31,24),"NOTME"}' |\
cat - all_possible_refpoints.pdb |\
awk '/NOTME/{++notme[substr($0,1,24)];next}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next} {typ=" "substr($0,18,3)" "}\
  typ !~ /HOH|EDO|NH4| CL | LI /{next}\
  {xyz=substr($0,31,24)}\
  notme[xyz]{next}\
  {print}' >! teleport_goals.pdb
pick.com 4.5 fofc.map | head
echo "symgen $smallSG" | pdbset xyzin pick.pdb xyzout bigpick.pdb >> /dev/null
egrep "^ATOM|^HETAT" bigpick.pdb >> teleport_goals.pdb

rm -f teleported.rst7
water_teleport_runme.com ${laStage}.rst7 maxmoves=$teleport_waters \
  teleport_goals.pdb notthese=current_restraints.pdb \
  mtzfile=cootme.mtz \
  mindist=$teleport_mindist \
  debug=$debug minimize=0 energycheck=0 |\
  tee teleport_${itr}.log 
#| awk 'NF==3 && $2=="="{next} {print}'

if(! -e teleported.rst7) then
  echo "no teleportation options"
  goto skipteleport
endif

# edit restraints file
cp current_restraints.pdb blunt_restraints.pdb
#rst2pdb_runme.com teleported.rst7 >! spikeweight.log
awk '$NF=="moved"' teleport_${itr}.log >! moved.pdb
#  neighbors_big.com all_possible_refpoints.pdb moved.pdb 1.5A -include  >> spikeweight.log
#  near_atoms3.com moved.pdb neighbors.pdb 2A onepair onewrong >> spikeweight.log
set teleport_B = `echo $teleport_weight $pdbscale | awk '{print $1/$2}'`
convert_pdb.awk -v BFAC=$teleport_B moved.pdb |\
egrep "^ATOM|^HETAT" >! spikeweight.pdb
awk '/^ATOM|^HETAT/{print $0,"new"}' spikeweight.pdb |\
cat - current_restraints.pdb |\
awk '/^CRYST1/{print} ! /^ATOM|^HETAT/{next}\
  {xyz=substr($0,31,24)}\
  seen[xyz]{next}\
  {print;++seen[xyz];next}' |\
cat - spikeweight.pdb  |\
awk '/^CRYST1/{print} ! /^ATOM|^HETAT/{next}\
  {xyz=substr($0,31,24)}\
  seen[xyz]{next}\
  {print;++seen[xyz];next}' |\
cat >! filled.pdb
combine_pdbs_runme.com spikeweight.pdb filled.pdb printref=1 >! spikeweight.log
set changes = `diff new.pdb current_restraints.pdb | egrep "^<" | wc -l`
mv new.pdb current_restraints.pdb
echo "spiked $changes weights in current_restraints.pdb"

echo "updating ref.crd"
update_centroid_positions_runme.com current_restraints.pdb outfile=ref.crd >> spikeweight.log

# check for bombs ?

cp ${laStage}.in teleported.in
if(-e ${laStage}.prmtop) cp ${laStage}.prmtop teleported.prmtop
set laStage = teleported
#set repick_now = 1

# update the this and noEP pdb files
rst2pdb_runme.com ${laStage}.rst7 outprefix=this >! rst2pdb.log
cat this.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  substr($0,77,2)!="XP" && substr($0,12,5)!=" EPW "{print}' |\
cat >! refme.pdb

skipteleport:

if ( $Breset_itr && $itr % $Breset_itr == 0 ) then
#   echo "resetting all B factors to 20"
#   mv Breset.pdb refme.pdb
endif



set just_repicked = 0
if ( $repick_itr && ( $itr % $repick_itr == 0 ) || $repick_now ) then
  set just_repicked = 1

  egrep "^CRYST" all_possible_refpoints.pdb | head -n 1 >! unfiltered_refpoints.pdb
  egrep -v "HOH|^END" all_possible_refpoints.pdb >> unfiltered_refpoints.pdb
  grep HOH current_restraints.pdb                >> unfiltered_refpoints.pdb
  grep HOH all_possible_refpoints.pdb            >> unfiltered_refpoints.pdb

  egrep "^CRYST|^ATOM|^HETAT" unfiltered_refpoints.pdb |\
  awk -v debug=$debug '/^CRYST/ && ! cryst{print;++cryst;next}\
       {id=substr($0,12,17);xyz=substr($0,31,24);\
        front=substr($0,1,30);back=substr($0,31);\
        ++seen[xyz];++seen[id]}\
    ! /HOH/{print;++seen[id];++seen[xyz];next}\
    seen[xyz]==1 && seen[id]==1{if(debug)print "DEBUG passing thru";\
       print;++seen[id];next}\
    seen[xyz]>1 && seen[id]==1{if(debug)print "DEBUG banking id",id;\
       ++f;bankf[f]=front;}\
    seen[xyz]>1{next}\
    {if(debug) print "DEBUG banking xyz",xyz;\
        ++x;bankx[x]=back}\
    END{if(debug)print "DEBUG banks:",f,x;\
      for(i=1;i<=x;++i){j=i;if(j>f)j=f;print bankf[j] bankx[i],i}}' |\
  convert_pdb.awk -v renumber=dedupe >! inherited_refpoints.pdb

  # remove close overlaps?
  egrep -v "HOH|^END" inherited_refpoints.pdb >! sorted_refpoints.pdb
  egrep "HOH" inherited_refpoints.pdb |\
  sort -k1.61gr >> sorted_refpoints.pdb

  # re-pick reference points
  echo "re-discovering nearest reference points"
  centroids_nearby_runme.com refme.pdb reffile=sorted_refpoints.pdb \
    softener=1 weight=$repick_maxweight maxdist=$repick_maxdist  \
    hohscale=$repick_hohscale \
    outfile=new_reference_points.pdb debug=$debug | tee c2r_${itr}.log | awk 'NF<=2 || NF==3 && $2=="="{next} {print}'

  echo "forgetting weights for big-move waters"
  # exclude recently teleported waters?
  egrep -h HOH current_restraints.pdb new_reference_points.pdb |\
  rmsd -v debug=1 |\
    awk -v r=$repick_maxdist '/^WARN/{print substr($0,10,17),"BIGMOVE"}\
      /moved/ && /HOH/ && substr($0,25)+0>r{print substr($0,1,17),"BIGMOVE"}' |\
  tee water_bigmoves.txt | wc -l
  cat water_bigmoves.txt current_restraints.pdb |\
  awk '$NF=="BIGMOVE"{id=substr($0,1,17);++bigmove[id];next}\
   ! /^ATOM|^HETAT/{print;next}\
     {id=substr($0,12,17)}\
   ! bigmove[id]{print}' |\
  cat >! pruned_restraints.pdb

  echo "migrating previous weights into new list"
  combine_pdbs_runme.com printref=1 saveXYZ \
     pruned_restraints.pdb new_reference_points.pdb \
     outfile=new_points_old_weights.pdb >> c2r_${itr}.log

  echo "updating xyz positions in ref.crd"
  cp -p ref.crd lastref.crd
  update_centroid_positions_runme.com \
     pdbfile=new_points_old_weights.pdb \
     orignames=orignames.pdb topfile=xtal.prmtop \
     outfile=ref.crd >> c2r_${itr}.log

  egrep "^CRYST|^ATOM|^HETAT" new_points_old_weights.pdb >! current_restraints.pdb

  # make sure teleported waters "stick"
  if(-e teleport_${itr}.log) then
    set teleport_B = `echo $teleport_weight $pdbscale | awk '{print $1/$2}'`
    cp current_restraints.pdb blunt_restraints.pdb
    awk '$NF=="moved"{print $0,"SPIKEME"}' teleport_${itr}.log >! probe.pdb
    neighbors_big.com blunt_restraints.pdb probe.pdb 2A -include
    awk '/^ATOM|^HETAT/{print $0,"SPIKEME"}' neighbors.pdb |\
    cat - probe.pdb blunt_restraints.pdb |\
    awk -v w=$teleport_B '! /^ATOM|^HETAT/{print;next}\
       {xyz=substr($0,31,24)}\
      $NF=="SPIKEME"{++sel[xyz];next}\
       {pre=substr($0,1,60);post=substr($0,67)}\
       sel[xyz]{printf("%s%6.2f%s\n",pre,w,post);next}\
       {print}' >! accented_restraints.pdb
    set accented = `diff blunt_restraints.pdb accented_restraints.pdb | egrep ">" | wc -l`
    echo "applied weight=$teleport_weight to $accented teleported waters."
    cp accented_restraints.pdb current_restraints.pdb
  endif
endif

if ( $filter_itr && ( $itr % $filter_itr == 0 ) ) then
  echo "looking for too-close water refpoints"
  set dist = `echo $filter_dist 0.1 | awk '{print $1+$2}'`
  convert_pdb.awk -v skip=protein current_restraints.pdb |\
    egrep -v "EG7| ZN " >! solvent_restraints.pdb 
  distanceify_runme.com current_restraints.pdb solvent_restraints.pdb \
     maxdist=$dist outfile=distances.txt >! dist.log
  egrep "HOH|EDO|SO4| LI | CL |NH4|ACY|ACT" restraints_challenge_vs_weight.txt |\
  awk '{print substr($0,length($0)-16) "|",$1,"|",$2,"|",$3,"CHALLENGE"}' |\
  cat - distances.txt |\
  awk -v d2=$filter_dist -F "|" '/CHALLENGE/{challenge[$1]=$2+0;d[$1]=$3+0;w[$1]=$4+0;next}\
     /DIST/{dist[$1,$2]=$3+0;dist[$2,$1]=dist[$1,$2]}\
     /DIST/ && ( $1==$2 || dist[$1,$2]>2.4 ){;next}\
     elim[$1] || elim[$2]{next}\
     {sol1=sol2=polya1=polya2=0;typ1=substr($1,7,3);typ2=substr($2,7,3)}\
     typ1~/HOH|EDO|SO4| LI | CL |NH4|ACY/{++sol1}\
     typ2~/HOH|EDO|SO4| LI | CL |NH4|ACY/{++sol2}\
     typ1~/EDO|SO4|ACY|ACT/{++polya1}\
     typ2~/EDO|SO4|ACY|ACT/{++polya2}\
     # do not eliminate like-to-like polyatom restraints \
     polya1 && typ1==typ2{next}\
     ! sol1 && sol2 && dist[$1,$2]<d2{print $2 "| ELIM";++elim[$2]}\
     sol1 && ! sol2 && dist[$1,$2]<d2{print $1 "| ELIM";++elim[$1]}\
     challenge[$1]>challenge[$2] && w[$1]>7.0 && sol1{print $1 "| ELIM";++elim[$1]}\
     challenge[$2]>challenge[$1] && w[$2]>7.0 && sol2{print $2 "| ELIM";++elim[$2]}\
     sol1 && sol2 && dist[$1,$2]<d2 && w[$1]>w[$2]*2{print $2 "| ELIM";++elim[$2]}\
     sol1 && sol2 && dist[$1,$2]<d2 && w[$2]>w[$1]*2{print $1 "| ELIM";++elim[$1]}' |\
  cat - current_restraints.pdb |\
  awk '$NF=="ELIM"{++elim[substr($0,1,16)];next}\
    ! /^ATOM|^HETAT/{print;next}\
     {id=substr($0,12,16)}\
    elim[id]{next}\
     {print}' >! filtered_refpoints.pdb

  set filtered = `diff current_restraints.pdb filtered_refpoints.pdb | awk '/^</' | wc -l`
  echo "filtered out $filtered solvent anchors as too close and/or too challenged"

  egrep "^CRYST|^ATOM|^HETAT" filtered_refpoints.pdb >! current_restraints.pdb
endif

if(-e teleported_waters_${itr}.pdb) then
  echo "teleported waters are now restrained as:"
  awk '/O   HOH/{print substr($0,12,17),"SEL"}' teleported_waters_${itr}.pdb |\
  cat - current_restraints.pdb |\
  awk '$NF=="SEL"{++sel[substr($0,1,17)];next}\
     ! /^ATOM|^HETAT/{next} {id=substr($0,12,17)}\
     sel[id]{print}'
endif

if( $remap_itr &&  $itr % $remap_itr == 0 ) then

  remap_waters_runme.com ${laStage}.rst7 current_restraints.pdb | tee remap_waters_${itr}.log
  if( $status ) then
    set BAD = "remap failed"
    goto exit
  endif

  cp remapped_restraints.pdb current_restraints.pdb

  echo "updating ref.crd"
  update_centroid_positions_runme.com current_restraints.pdb outfile=ref.crd >> remap_waters_${itr}.log

  # check for bombs ?

  rst2pdb_runme.com ${laStage}.rst7 outprefix=this >! rst2pdb.log
  cat this.pdb |\
  awk '! /^ATOM|^HETAT/{print;next}\
     substr($0,77,2)!="XP" && substr($0,12,5)!=" EPW "{print}' |\
  cat >! refme.pdb

  cp ${laStage}.in remapped.in
  set laStage = remapped

endif

if( "$thrubond_avg_weight" != "0" ) then
  echo "averaging weights through bonds with avgfac=$thrubond_avg_weight"
  cp current_restraints.pdb presmooth_restraints_for_${itr}.pdb
  thrubond_avgB_runme.com current_restraints.pdb debug=$debug \
    avgfac=$thrubond_avg_weight spread=$thrubond_avg_weight_spread >> restraint_update_${itr}.log
  cp thrubondBavg_restraints.pdb current_restraints.pdb
endif

set kT = `echo $temperature 0.00198587763237322 | awk '{print $1*$2}'`

if( "$weight_power" != "1" ) then
  set small_weight = `echo $kT $weight_power_kTmult | awk '{print $1*$2}'`
  set small_B = `echo $small_weight $pdbscale | awk '{print $1/$2}'`
  echo "raising weights below $small_weight to the power of $weight_power"
  cp current_restraints.pdb prepower_restraints.pdb
  echo "$small_B $weight_power" |\
  cat - prepower_restraints.pdb |\
  awk 'NR==1{B0=$1;power=$2;next}\
    ! /^ATOM|^HETAT/{print;next}\
    {B=substr($0,61,6)+0;pre=substr($0,1,60);post=substr($0,67)}\
    B>B0{print;next}\
    {B=B0*(B/B0)**power;\
     printf("%s%6.2f%s\n",pre,B,post)}' |\
  cat >! restraints_powered.pdb
  rmsd prepower_restraints.pdb restraints_powered.pdb
  cp restraints_powered.pdb current_restraints.pdb
endif



echo "changes in restraints:"
@ prev = ( $itr - 1 )
cat restraints_for_${prev}.pdb current_restraints.pdb |\
 egrep -v "HOH|EDO" |\
 awk '! /^ATOM|^HETAT/{next} {id=substr($0,12,15);B=substr($0,61,6)+0}\
    B==0 && zeroB[id]{next} B==0{++zeroB[id]} {print}' |\
 rmsd >! rmsd.txt
head rmsd.txt | grep -v pairs
set changes = `grep WARN rmsd.txt | wc -l`
set pairs = `awk '/atom pairs found/{print $1;exit}' rmsd.txt`

egrep -h "HOH|EDO" restraints_for_${prev}.pdb current_restraints.pdb |\
 awk '! /^ATOM|^HETAT/{next} {id=substr($0,12,15);B=substr($0,61,6)+0}\
    B==0 && zeroB[id]{next} B==0{++zeroB[id]} {print}' |\
 rmsd >! rmsdHOH.txt
set HOHchanges = `grep WARN rmsdHOH.txt | wc -l`
set HOHpairs = `awk '/atom pairs found/{print $1;exit}' rmsdHOH.txt`

echo "$itr   $changes $pairs   $HOHchanges $HOHpairs changes" | tee -a changes_vs_itr.txt


# try deleting forever, tight restraints that are highly challenged 
if (  $delete_badrest_itr && $itr % $delete_badrest_itr == 0 ) then
  echo "deleting highly challenged refpoints "

  cp current_restraints.pdb predelrest_restraints_for_${itr}.pdb

  delete_worst_restraints.com outprefix="" | tee delete_worst_${itr}.log 

endif


# try releasing alt-conf atoms
if ( ( ! $just_repicked ) && $delete_altconf_itr && $itr % $delete_altconf_itr == 0 ) then
  echo "deleting all alt-conf restraints "

   cp current_restraints.pdb predelalts_restraints_for_${itr}.pdb

   awk 'substr($0,17,1)!=" "{print substr($0,31,24),"XYZ"}' all_possible_refpoints.pdb |\
   cat - current_restraints.pdb |\
   awk '$NF=="XYZ"{++altconf[substr($0,1,24)];next}\
     ! /^ATOM|^HETAT/{print;next} {xyz=substr($0,31,24)}\
     ! altconf[xyz]{print}' |\
   cat >! oneconf_restraints.pdb

  cp oneconf_restraints.pdb current_restraints.pdb

endif

set best_itr = `sort -k2g fofc_Rplot.txt | awk '{print $1;exit}'`
sort -k2g fofc_Rplot.txt |\
 awk 'NR==1{best=$1}\
  sqrt(($1-best)^2)<5{print $2}' |\
awk '{++n;v[n]=$1;sum+=$1}\
  END{if(n)avg=sum/n;\
    for(i=1;i<=n;++i){sumd+=(v[i]-avg)**2};\
    if(n)rmsd=sqrt(sumd/n);\
    if(n)print avg,rmsd,n}' >! temp.txt
set bestR = `cat temp.txt`

awk -v itr=$randel_last '$1>itr{print $2}' fofc_Rplot.txt |\
tail -n 2 |\
awk '{++n;v[n]=$1;sum+=$1}\
  END{if(n)avg=sum/n;\
    for(i=1;i<=n;++i){sumd+=(v[i]-avg)**2};\
    if(n)print avg,sqrt(sumd/n),n}' >! temp.txt
set thisR = `cat temp.txt`

set randel_trigger_now = 0
if( $randel_trigger ) then
  set randel_trigger_now = `echo $thisR $bestR | awk '{print ($1<($4+$5/2) && $6>=5 && $3>=2)}'`
  echo "randel_trigger_now = $randel_trigger_now $itr  $thisR $bestR"
endif


# try releasing random atoms
if ( $randel_itr && ( ( ( $itr - $randel_last ) > $randel_itr ) || $randel_trigger_now ) ) then
   cp current_restraints.pdb predel_restraints_for_${itr}.pdb

   if(! -e randel_selections.txt ) then
      echo $randel_fraction |\
      cat - all_possible_refpoints.pdb |\
      awk 'NR==1{f=$1;next}\
        ! /^ATOM|^HETAT/{next} {++i;\
        bin=int(rand()/f+1);B=bin/100;\
        printf("%s%6.2f%s\n",substr($0,1,60),B,substr($0,67))}' |\
      cat >! randB.pdb
      xsame_runme.com randB.pdb 
      cat xsame_restraints.pdb |\
      awk '{bin=substr($0,61,6)*100;id=substr($0,12,19);xyz=substr($0,31,24);\
         print bin,"|"xyz"|"id"| RANDEL"}' |\
      cat >! randel_selections.txt
   endif
   set bins = `awk '{print $1}' randel_selections.txt | sort -u | sort -g | wc -l`
   if(-e randel_bin.txt) then
      set randel_bin = `awk '{print $NF;exit}' randel_bin.txt`
      echo "taking randel_bin=$randel_bin from randel_bin.txt"
   endif
   if( "$randel_bin" == "" ) then
      echo "WARNING: starting at beginning of randel bin list"
      set randel_bin = 0
   endif
   @ randel_bin = ( $randel_bin + 1 )
   if( $randel_bin > $bins ) set randel_bin = 1
   echo "deleting random restraints, bin: $randel_bin / $bins at itr: $itr"

   echo $randel_bin | tee randel_bin.txt |\
   cat - randel_selections.txt current_restraints.pdb |\
   awk 'NR==1{sel=$1;print "REMARK bin="sel;next}\
     $NF=="RANDEL"{split($0,w,"|");bin[w[2]]=w[1];next}\
     ! /^ATOM|^HETAT/{print;next}\
     {xyz=substr($0,31,24)}\
     bin[xyz]==sel{next}\
     {print}' |\
   cat >! randel_restraints.pdb

   diff randel_restraints.pdb current_restraints.pdb |\
   awk '/^>/{print substr($0,3)}' |\
   cat >! randelled_restraints_${itr}.pdb
   set diff = `cat randelled_restraints_${itr}.pdb | wc -l`
   echo "$diff restraints deleted"

  cp randel_restraints.pdb current_restraints.pdb
  set randel_last = $itr
endif





# try releasing tightest restraints, see if they grow back
if ( $release_itr && $itr % $release_itr == 0 || "$release_trigger" != "0" && $pegged_weight ) then
  echo "releasing most-challenged restraints "
  set tootight = `echo 9 $pdbscale | awk '{print $1/$2}'`
  release_worst_restraints_runme.com restraintpdb=current_restraints.pdb \
    refinedpdb=refme.pdb radius=$release_radius pdbscale=$pdbscale \
    tootight=9 maxbad=$release_maxbad \
    resetweight=$release_weight sigma=5 | tee release_worst_${itr}.log

  cp current_restraints.pdb prerelease_restraints_for_${itr}.pdb
  cp new_restraints.pdb current_restraints.pdb

  if( $release_fftBreset ) then
    echo "setting fft_B factor back to $release_fftBreset"
    set fft_B = $release_fftBreset
  endif
  if( $release_Breset ) then
    echo "resetting atomic B factors"
    convert_pdb.awk -v BFAC=20 refme.pdb >! new.pdb
    mv new.pdb refme.pdb
  endif
endif

# reject reference points that are not in density
if ( $rhocheck_itr && $itr % $rhocheck_itr == 0 ) then
  echo -n "checking density > $min_ref_rho : "
  set pdbfile = current_restraints.pdb
  #set pdbfile = all_possible_refpoints.pdb
  cp $pdbfile pre_rhocheck.pdb
  if( ! -e reference_sigma.map) then
    echo scale sigma | mapmask mapin reference.map mapout reference_sigma.map > /dev/null
  endif
  check_map_sites.com reference_sigma.map $pdbfile |\
  cat - $pdbfile |\
  awk 'NF==1{++i;rho[i]=$1;next}\
    ! /^ATOM|^HETAT/{print;next}\
    {++n;print $0,"          ",rho[n]}' |\
  tee rholabeled.pdb |\
  awk -v mrr=$min_ref_rho '! /^ATOM|^HETAT/{print;next}\
     $NF>mrr{print substr($0,1,80)}' |\
  cat >! indensity.pdb

  set rejected = `rmsd indensity.pdb $pdbfile | egrep "WARNING" | wc -l`
  echo "$rejected out-of-density reference points rejected"
  cp indensity.pdb $pdbfile
endif


# keep CAs restrained, no matter what
if ( $reinin_itr && $itr % $reinin_itr == 0 ) then
  # taking all defaults
  echo -n "checking CAs : "
  cp current_restraints.pdb pre_clip_CA.pdb
  reinin_errant_CA_refpoints.com outfile=clip_CA.pdb | tee reinin.log | awk '$2!="="'

  rmsd clip_CA.pdb current_restraints.pdb | egrep "MAXD.Bfac|WARNING"
  cp clip_CA.pdb current_restraints.pdb
endif

# make sure any ligands stay put
set test = `egrep "EG7 | ZN " current_restraints.pdb | wc -l`
if( $test ) then

  echo "ensuring slight restraints on EG7|ZN "
  egrep "EG7| ZN " all_possible_refpoints.pdb >! lig.pdb
  cat current_restraints.pdb lig.pdb |\
  awk '{id=substr($0,12,5)" "substr($0,18,12)} ! seen[id]{print;++seen[id]}' |\
  awk -v mlw=$min_lig_weight -v scale=$pdbscale '! /^ATOM|^HETAT/{print;next}\
    ! /EG7| ZN /{print;next}\
    {B=substr($0,61,6)+0;pre=substr($0,1,60);post=substr($0,67,12);\
     w=B*scale}\
    w<mlw{w=mlw}\
    {printf("%s%6.2f%s    NOT_A_BOMB\n",pre,w/scale,post)}' |\
  cat >! wrongorder.pdb
  combine_pdbs_runme.com wrongorder.pdb refme.pdb outfile=clip_ligand.pdb > /dev/null
  rmsd clip_ligand.pdb current_restraints.pdb | egrep "MAXD.Bfac"
  cp clip_ligand.pdb current_restraints.pdb

endif

if( $ambig_same_weight ) then
  echo "same-weight check:"
  egrep ".D. ASN" current_restraints.pdb | head -n 2
endif

# keep track of teleports
if(-e teleport_${itr}.log) then
  awk '$NF=="moved"{print $0,"SEL"}' teleport_${itr}.log |\
  cat - current_restraints.pdb |\
  awk '/^CRYST1/{print} ! /^ATOM|^HETAT/{next}\
    {xyz=substr($0,31,24);num[xyz]=$(NF-1)}\
    $NF=="SEL"{++sel[xyz];next}\
    sel[xyz]{print substr($0,1,80);next}' |\
  cat >! teleported_waters_${itr}.pdb
  set count = `egrep "^ATOM|^HETAT" teleported_waters_${itr}.pdb | wc -l`
  echo "$count teleported waters left"
endif

# final bomb check
echo "final check for errant refpoints..."
rm -f defused_restraints.pdb
restraint_bomb_detector.com ${laStage}.rst7 \
  refpoints=current_restraints.pdb \
  outfile=defused_restraints.pdb \
  outcrd=ref.crd >! restraint_bomb_check.log
if( $status ) then
  cat restraint_bomb_check.log
  set BAD = "too many restraint bombs detected"
  goto exit
endif
egrep "WARNING|worst|^ATOM|^HETAT" restraint_bomb_check.log
egrep -v "^recommend:|^making ref.crd|^mv " restraint_bomb_check.log |\
tail -n 1 
diff current_restraints.pdb defused_restraints.pdb > /dev/null
if( $status ) then
  cat restraint_bomb_check.log
  echo "using defused_restraints.pdb"
  cp defused_restraints.pdb current_restraints.pdb
endif

# make easier to see in Coot: color range from 10 to 100
cat current_restraints.pdb |\
awk -v pdbscale=$pdbscale '! /^ATOM|^HETAT/{print;next}\
  {B=substr($0,61,6)+0;pre=substr($0,1,60);post=substr($0,67);\
   printf("%s%6.2f%s\n",pre,B*pdbscale*10,post)}' |\
cat >! cootme_restraints.pdb


# create ref.crd from most recent rst7 and current_restraints.pdb
rst2pdb_runme.com ${laStage}.rst7 now.pdb Bfactors=none >! final_centroids.log
update_centroid_positions_runme.com \
     pdbfile=current_restraints.pdb \
     orignames=now.pdb topfile=xtal.prmtop \
     outfile=ref.crd >> final_centroids.log



echo "replacing restraint list in ${laStage}.in with new list in ${Stage}.in"
restraintlist2amber.com list=current_restraints.pdb \
  pdbscale=$pdbscale \
  minwt=$min_weight \
  orignames=orignames.pdb \
  allatom_weight=$allatom_weight \
  ${laStage}.in \
  outfile=${Stage}.in >> restraint_update_${itr}.log

# set itr = `ls -1rt amber_*[0-9].rst7 | tail -n 1 | awk -F "_" '{print $NF+0}'`
# set Stage = `ls -1rt *.rst7 | egrep -v "ref.rst7|equi|unwrapped" | tail -n 1 | awk -F "[.]" '{print $1}'`
cat ${Stage}.in |\
awk '/Specific/{atom=$2;getline;w=$1}\
  /^RES/{for(n=$2;n<=$3;++n){print w,atom,n}}' |\
sort -gr |\
tee actual_weights.txt |\
cat - orignames.pdb |\
awk 'NF==3{++sel[$2,$3];next}\
  ! /^ATOM|^HETAT/{next}\
    {atom=substr($0,12,5);gsub(" ","",atom)}\
   sel[atom,$NF]{print substr($0,12,17),"SEL"}' |\
cat - current_restraints.pdb |\
awk '$NF=="SEL"{++sel[substr($0,1,17)]}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next}\
  {id=substr($0,12,17)}\
  sel[id]{print}' |\
cat >! actual_restraints.pdb

# now set in stone
cp actual_restraints.pdb current_restraints.pdb
cp current_restraints.pdb restraints_for_${itr}.pdb

echo -n "max min density: "
awk '{split(FILENAME,w,"_");f=w[3]+0} /Maximum dens/{print f,m,$NF} /Minimum dens/{m=$NF}' restraint_update_${itr}.log |\
 tee -a dens_vs_itr.txt

awk -v scale=$pdbscale '/^ATOM|^HETAT/{print substr($0,61,6)*scale,substr($0,12,15)}' current_restraints.pdb |\
   sort -gr >! sorted_weights.txt
set avg = `awk '{sum+=$1;++n} END{print sum/n}' sorted_weights.txt`
set medmad = `awk '{print $1}' sorted_weights.txt | median.awk`
set max = `awk '{print $1;exit}' sorted_weights.txt`
set worst = `awk '! t{$1="";print;++t} ! / ZN | SE /{print;exit}' sorted_weights.txt`
set mults = `awk -F ":" '/^highest|^lowest/{getline;print $2+0}' restraint_update_${itr}.log`
echo -n "worstweight: "
echo "$itr $max $avg $mults  $medmad  $worst" | tee -a worstweight_vs_itr.txt

set pegged_weight = `echo $max $max_weight | awk '{print ( $1 >= $2 )}'`
if( $pegged_weight ) echo "WARNING: max weight pegged at limit"

echo -n "worst atom: "
sort -k1.61gr current_restraints.pdb | awk '/^ATOM|^HETAT/{print;exit}' >! worstatom.pdb
 awk -v scale=$pdbscale -v itr=$itr '! /^ATOM|^HETAT/{next}\
   {id=substr($0,12,15)} NR==1{++sel[id];next} \
  sel[id]{print itr,substr($0,61,6)*scale,$0}' worstatom.pdb current_restraints.pdb |\
tee -a worstatom_weight_hist.txt
 awk '! /^ATOM|^HETAT/{next}\
   {id=substr($0,12,15)} NR==1{++sel[id];next} \
  sel[id]{j=split(FILENAME,w,"_");print w[j]+0,$0}' worstatom.pdb refme.pdb |\
sort -g |\
awk '{i=index($0,"ATOM");pre=substr($0,i,22);post=substr($0,i+26);\
  {printf("%s%4d%s\n",pre,$1,post)}}' |\
tee -a worstatom_path.pdb > /dev/null


cat ${laStage}.in |\
awk '/Specific/ && $NF+0<0.01{exit}\
  /^RES/{for(i=$2;i<=$3;++i)print i}' |\
 sort -u | sort -g |\
awk 'NR==1{s=e=$1;next} $1==e+1{e=$1;next} {print s"-"e;s=e=$1} END{print s"-"e}' |\
awk -F "-" '$1==$2{print $1;next} {print}' |\
sort -u | sort -g >! rest_ranges.txt 
set rest_ranges = `cat rest_ranges.txt `
set rest_ranges = `echo $rest_ranges | awk '{gsub(" ",",");print}'`


# update chiral and omega restraint lists
wait
set omegalogs = `ls -1rt omegalyze_*.log | tail -n 100`
cat $omegalogs |\
awk -F ":" '/ to /{;\
  dev=(($3+270)%180)-90;\
  rn1=substr($0,2,6);rn2=substr($0,17,6);\
  if(! seen[rn1]) print "BAD",rn1,rn2,dev;\
  ++seen[rn1]}' |\
cat >! bad_omega.txt
cat bad_omega.txt orignames.pdb |\
awk '/^BAD/{++bonds;\
  dev[bonds]=$NF;\
  rn1=substr($0,5,6);rn2=substr($0,12,6);\
  sel1[rn1]=sel2[rn2]=bonds;res[bonds]=rn1;next}\
 /^ATOM|^HETAT/{++a;atom=substr($0,12,5);gsub(" ","",atom);\
    rn=substr($0,22,6)}\
  sel1[rn] && atom=="CA"{b=sel1[rn];list[b] = list[b] a ",";}\
  sel1[rn] && atom=="C" {b=sel1[rn];list[b] = list[b] a ",";}\
  sel2[rn] && atom=="N" {b=sel2[rn];list[b] = list[b] a ",";}\
  sel2[rn] && atom=="CA"{b=sel2[rn];list[b] = list[b] a ",";}\
  END{for(b=1;b<=bonds;++b){\
    print "# trans-omega for",res[b],"off by",dev[b];\
    end="&end";if(b==1)end="";\
    print " &rst iat=" list[b],end;\
    if(b==1)print "  r1=150., r2=170., r3=190., r4=210., rk2 =50, rk3=50,  &end"}}' |\
cat >! bad_omega.rst


set chiralogs = `ls -1rt chiralyze_*.log | tail -n 100`
cat $chiralogs |\
awk '$2=="ideal" || $2=="delta:" || NF==0{next}\
  {rn=substr($0,4,6);atom=substr($0,15,5);\
   gsub(" ","",atom);\
   list=list" "atom;}\
     /*sigma$/{print "BAD",rn,list;++bads;list=""}' |\
sort -u >! badchir.txt
cat badchir.txt orignames.pdb |\
awk '/^BAD/{++baddies;\
   rn=substr($0,5,6);list=substr($0,12);n=split(list,a);\
   res[baddies]=rn;\
   for(i=1;i<=n;++i){\
      ++sel[rn,a[i]];\
      atname[rn,i]=a[i];\
   };\
   next}\
 /^ATOM|^HETAT/{++an;atom=substr($0,12,5);gsub(" ","",atom);\
    rn=substr($0,22,6)}\
  sel[rn,atom]{atomnum[rn,atom]=an}\
  END{for(b=1;b<=baddies;++b){\
    rn=res[b];list="";\
    for(i=1;i<=4;++i){\
      list=list atomnum[rn,atname[rn,i]]",";\
    }\
    print "# chirality for",res[b];\
    end="&end";if(b==1)end="";\
    print " &rst iat=" list,end;\
    if(b==1)print "  r1=10., r2=60.,  r3=80.,  r4=130., rk2 =10, rk3=10,  &end"}}' |\
cat >! bad_chir.rst



# adjust temp for unrestrained atoms
if ( "$thermostat" != "" ) then
  set prevtemp = `awk -F "=" '/temp0/{print $2+0}' ${laStage}.in | tail -n 1`
  echo "previous temperature setting was: $prevtemp"
  if("$prevtemp" == "0") set prevtemp = "$thermostat"

  rm -f tempW.dat tempU.dat tempUW.dat >& /dev/null
    cat << EOF >! cpptraj_temp.in
temperature Tw :WAT ntc 2 out tempW.dat
temperature Tr !:$rest_ranges ntc 2 out tempU.dat
strip :$rest_ranges
temperature Tuw :WAT ntc 2 out tempUW.dat
EOF
  cat cpptraj_temp.in |\
    cpptraj -p xtal.prmtop -y ${laStage}.rst7 >! cpptraj_temp.log

  set Tw = `awk 'NR>1{sum+=$2;++n} END{print sum/n}' tempW.dat`
  set Tu = `awk 'NR>1{sum+=$2;++n} END{print sum/n}' tempU.dat`
  set Tuw = `awk 'NR>1{sum+=$2;++n} END{print sum/n}' tempUW.dat`
  echo "water temperature: $Tw"
  echo "unrestrained residues temperature: $Tu"
  echo "unrestrained water temperature: $Tuw"
  set temperature = `echo $prevtemp $thermostat $Tuw | awk '{dT=($2-$3);print $1+0.5*(dT)}'`
  #echo "changing thermostat temperature to $temperature"
endif

if( "$temp_ramp" != "none" ) then
  set mintemp = `echo $temp_ramp | awk -F "-" 'NF>2{print $NF}'`
  if( "$mintemp" == "" ) set mintemp = 0
  set temperature = `echo $temperature $temp_ramp | awk '/x/{print $1*$2;exit} {print $1+$2}'`
  set temperature = `echo $temperature $mintemp | awk '$1<=$2{$1=$2} {print $1}'`
  echo "adjusting temperature to $temperature"
  if ( "$thermostat" != "" ) then
    set thermostat = `echo $thermostat $temp_ramp | awk '/x/{print $1*$2;exit} {print $1+$2}'`
    set thermostat = `echo $thermostat | awk '$1<=0{$1=0} {print}'`
    echo "adjusting thermostat temperature to $thermostat"
  endif
endif

if( "$equi_ns_ramp" != "none" ) then
    set equi_ns = `echo $equi_ns $equi_ns_ramp | awk '/x/{print $1*$2;exit} {print $1+$2}'`
    set equi_ns = `echo $equi_ns | awk '$1<=0.1{$1=0.1} {print}'`
    echo "adjusting amber equilibration time to $equi_ns ns"
endif
if( "$prod_ns_ramp" != "none" ) then
    set prod_ns = `echo $prod_ns $prod_ns_ramp | awk '/x/{print $1*$2;exit} {print $1+$2}'`
    set prod_ns = `echo $prod_ns | awk '$1<=0.1{$1=0.1} {print}'`
    echo "adjusting amber production time to $prod_ns ns"
endif

echo "temperature = $temperature"
echo $temperature |\
cat - ${Stage}.in |\
awk 'NR==1{temp=$1;next} \
  {space=substr($0,1,index($0,$1)-1)}\
  / temp0/{split($1,w,"=");print space w[1] "=" temp ",";next}\
  {print}' |\
cat >! tempfile.in
mv tempfile.in ${Stage}.in

if( ! $barostat ) then
 cat ${Stage}.in |\
 awk '{space=substr($0,1,index($0,$1)-1)}\
  / ntb=/{split($1,w,"=");print space "ntb=1,ntp=0,";next}\
  {print}' |\
 cat >! tempfile.in
 mv tempfile.in ${Stage}.in
endif

set nsnb = `head -n 100 ${Stage}.in | awk -F "=" '/ nsnb=/{print $2+0;exit}'`
echo "equi_ns = $equi_ns"
echo $equi_ns $equi_dt $nsnb |\
cat - ${Stage}.in |\
awk 'NR==1{ns=$1;dt=$2+0;nsnb=$3+0;if(dt<1e-6)dt=0.002;nstlim=int($1*1000/dt);next} \
  {space=substr($0,1,index($0,$1)-1)}\
  # shorter time step \
  / dt=/{split($1,w,"=");print space w[1] "=" dt ",";next}\
  # turn off shake \
#  / ntf=2/{print space "ntf=1";next}\
#  / ntc=2/{print space "ntc=1";next}\
  # nonbond updates every cycle \
#  ! nsnb && / cut=/{print space "nsnb=1,"}\
  # just a few cycles \
  / nstlim=/{split($1,w,"=");print space w[1] "=" nstlim ",";next}\
  # heavy thermostat \
#  / gamma_ln=/{split($1,w,"=");print space w[1] "=" w[2]*10 ",";next}\
  {print}' |\
cat >! settle.in

echo $equi_ns $equi_dt |\
cat - ${Stage}.in |\
awk 'NR==1{ns=$1;dt=$2+0;if(dt<1e-6)dt=0.002;nstlim=int($1*1000/dt);next} \
  {space=substr($0,1,index($0,$1)-1)}\
  # still shorter time step \
  / dt=/{split($1,w,"=");print space w[1] "=" dt ",";next}\
  # just a few cycles \
  / nstlim=/{split($1,w,"=");print space w[1] "=" nstlim ",";next}\
  {print}' |\
cat >! equi.in

echo $barometer_cycles 0.0005 |\
cat - ${Stage}.in |\
awk 'NR==1{n=$1;dt=$2;next}\
  {space=substr($0,1,index($0,$1)-1)}\
# just a few cycles \
/ nstlim=/{split($1,w,"=");print space w[1] "=" n ",";next}\
/ dt=/{split($1,w,"=");print space w[1] "=" dt ",";next}\
/ ntpr=/{split($1,w,"=");print space w[1] "=" n/10 ",";next}\
/ ntwx=| ntwr=| ntr=/{split($1,w,"=");print space w[1] "=" 0 ",";next}\
/ ntb=/{split($1,w,"=");print space "ntb=2,ntp=4,";next}\
/     Specific /{exit}\
/ restraint|wt type=|DISANG|LISTIN=POUT|dummy |nmropt=1/{next}\
{print}' |\
cat >! barometer.in

echo "prod_ns = $prod_ns"
echo $prod_ns $dt $write_ps |\
cat - ${Stage}.in |\
awk 'NR==1{prod=$1;dt=$2+0;write=$3+0;if(dt<1e-6)dt=0.002;\
     nstlim=int(prod*1000/dt);\
     ntwx=int(write/dt);next} \
  {space=substr($0,1,index($0,$1)-1)}\
  / nstlim=/{split($1,w,"=");print space w[1] "=" nstlim ",";next}\
  / ntwx=/{split($1,w,"=");print space w[1] "=" ntwx ",";next}\
  / dt=/{split($1,w,"=");print space w[1] "=" dt ",";next}\
  {print}' |\
cat >! tempfile.in
mv tempfile.in ${Stage}.in

if(-e chir_omega.rst) then
  echo "adjusting chiral/omega weights in chir_omega.rst : $chiral_weight $omega_weight"
  echo $chiral_weight $omega_weight |\
  cat - chir_omega.rst |\
  awk 'NR==1{cw=$1;ow=$2;next}\
    $1=="r1=10.," {print "   r1=10., r2=60.,  r3=80.,  r4=130., rk2 ="cw", rk3="cw",  &end";next}\
    $1=="r1=150.,"{print "  r1=150., r2=170., r3=190., r4=210., rk2 ="ow", rk3="ow",  &end";next}\
    {print}' |\
  cat >! new.txt
  mv new.txt chir_omega.rst
endif


set disang = `grep DISANG ${Stage}.in | wc -l`
if( ! $disang && ( $omega_weight != 0 ) && -e chir_omega.rst) then
  echo "applying DISANG restraints"

  cat ${Stage}.in |\
  awk 'NF==1 && $1=="/"{;\
    print "  nmropt=1, /";\
    print "&wt type=\047END\047 /";\
    print "DISANG=chir_omega.rst";\
    print " &dummy  i=1, ";}\
   {print}' |\
  cat >! tempfile.in
  mv tempfile.in ${Stage}.in
endif
if( $disang && $omega_weight == 0 && $chiral_weight == 0 ) then
  echo "removing chiral/omega restraints"

  egrep -v "dummy|DISANG|LISTIN=PO|wt type=|nmropt" ${Stage}.in >! tempfile.in
  mv tempfile.in ${Stage}.in
endif


foreach substage ( settle equi )

  rm -f ${Stage}_${substage}.rst7
  echo "${substage}-ing $Stage "
  $pmemd -O -i ${substage}.in -o ${Stage}_${substage}.out \
   -p xtal.prmtop \
   -c ${laStage}.rst7 \
   -ref ref.crd \
   -r ${Stage}_${substage}.rst7 \
   -x ${Stage}_${substage}.nc \
   -inf ${Stage}_${substage}.mdinfo
  if($status) then
    set BAD = "amber run failed at ${Stage}_${substage}"
    goto exit
  endif

  grep NaN ${Stage}_${substage}.out
  if(! $status) break
  egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' ${Stage}_${substage}.out | grep '\*\*\*\*\*'
  if(! $status) break

  set laStage = ${Stage}_${substage}

end
  grep NaN ${Stage}_${substage}.out
  if(! $status) break
  egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' ${Stage}_${substage}.out | grep '\*\*\*\*\*' 
  if(! $status) break

  rm -f ${Stage}.rst7
  echo "running $Stage "
  $pmemd -O -i ${Stage}.in -o ${Stage}.out \
   -p xtal.prmtop \
   -c ${laStage}.rst7 \
   -ref ref.crd \
   -r ${Stage}.rst7 \
   -x ${Stage}.nc \
   -inf ${Stage}.mdinfo
  if($status) then
    set BAD = "amber run failed at ${Stage}"
    goto exit
  endif

  grep NaN ${Stage}.out
  if(! $status) break
  egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' ${Stage}.out | grep '\*\*\*\*\*'
  if(! $status) break

  echo "running barometer "
  $pmemd -O -i barometer.in -o barometer_${itr}.out \
   -p xtal.prmtop \
   -c ${Stage}.rst7 \
   -ref ref.crd \
   -r deleteme.rst7 \
   -inf barometer.mdinfo

#rm -f deleteme.nc deleteme.rst7

if(-e ./exit) then
  echo "exiting because ./exit exists."
  rm -f ./exit
  break
endif

end

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

exit

#####################################################################################################
#
#  notes and post-analysis
#

set sruncpu = "srun --partition=refmac --exclude=crush18"
set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"
set itr = `ls -lrt amber_*.rst7 | grep -v equi | awk -F "_" '{print $NF+0}' | tail -n 1`
set laStage = amber_${itr}

set Stage = Min
restraintlist2amber.com list=restraints_for_${itr}.pdb \
  orignames=orignames.pdb \
  ../opt17/Min.in \
  outfile=${Stage}.in 

  $pmemd -O -i ${Stage}.in -o ${Stage}.out \
   -p xtal.prmtop \
   -c ${laStage}.rst7 \
   -ref ref.crd \
   -r ${Stage}.rst7 \
   -x ${Stage}.nc \
   -inf ${Stage}.mdinfo &

echo "outtraj xyz.pdb include_ep" |\
cpptraj -p xtal.prmtop -y ${Stage}.rst7 >> $logfile

echo "combining with B factors from Bfac.pdb"
egrep "^SSBOND|^LINK|^CISP|^CRYST" Bfac.pdb >! pinchme.pdb
# take coordinates only, using names from starting point
awk '/^ATOM|^HETAT/{print $0,"ORIG"}' Bfac.pdb |\
cat - xyz.pdb |\
awk '$NF=="ORIG"{++o;pre[o]=substr($0,1,30);post[o]=substr($0,55,length($0)-55-4);next}\
  /^CRYST/{next} ! /^ATOM|^HETAT/{print;next}\
      {++n;\
      printf("%s%s%s\n",pre[n],substr($0,31,24),post[n])}' |\
convert_pdb.awk -v skip=EP >> pinchme.pdb

pinchedwater_runme.com pinchme.pdb current_restraints.pdb energy_thresh=0.1 | tee pinchedwater_Min.log


egrep "^CRYST" tempfile_rud_refpoints_small.pdb >! sfallme.pdb
egrep -h  "^ATOM" ../opt*/badwaters.pdb | convert_pdb.awk >> sfallme.pdb 
sfall xyzin sfallme.pdb mapout sfalled.map << EOF >! sfall.log
MODE ATMMAP
SYMM 1
grid 128 128 128
EOF



  release_worst_restraints_runme.com restraintpdb=current_restraints.pdb \
    refinedpdb=refme.pdb radius=6 \
    tootight=10 resetweight=0.1 sigma=5
  mv new_restraints.pdb current_restraints.pdb




set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"
set itr = `ls -lrt amber_*.rst7 | grep -v equi | awk -F "_" '{print $NF+0}' | tail -n 1`
set laStage = amber_${itr}
set Stage = `ls -1rt amber_*.in | grep -v equi | awk -F "." '{print $1}' | tail -n 1`




# post analysis

grep PRESS amber_*.out |\
 awk '$(NF-3)>200{itr=substr($0,7)+0;++seen[itr];print itr+seen[itr]/50,$NF,$(NF-3)}' |\
 sort -g | tee press.log

rm -f nwater_vs_press.txt
foreach opt ( opt*/ )
   set nwater = `tail -n 5 ${opt}/refme.pdb | awk '/^ATOM/{print substr($0,23,5)}' | tail -n 1`
   set logs = `ls -1rt ${opt}/amber_*.out | tail -n 10`
   set press = `grep PRESS $logs | awk '{print $NF}' | avg.awk | awk '{print $1;exit}'`
   echo "$nwater $press $opt" | tee -a nwater_vs_press.txt
end


awk '{split(FILENAME,w,"_");f=w[3]+0} /Maximum dens/{print f,m,$NF} /Minimum dens/{m=$NF}' restraint_update_*.log |\
 sort -g | tee dens_vs_itr.txt

sort -k1.61gr current_restraints.pdb | awk '/^ATOM|^HETAT/{print;exit}' >! worstatom.pdb
 awk -v scale=$pdbscale '! /^ATOM|^HETAT/{next}\
   {id=substr($0,12,15)} NR==1{++sel[id];next} \
  sel[id]{split(FILENAME,w,"_");print w[3]+0,substr($0,61,6)*scale}' worstatom.pdb restraints_for_*.pdb |\
sort -g |\
 tee worstatom_weight_hist.txt
 awk '! /^ATOM|^HETAT/{next}\
   {id=substr($0,12,15)} NR==1{++sel[id];next} \
  sel[id]{j=split(FILENAME,w,"_");print w[j]+0,$0}' worstatom.pdb refmacout_*.pdb |\
sort -g |\
awk '{i=match($0,"ATOM|HETAT");pre=substr($0,i,22);post=substr($0,i+26);\
  {printf("%s%4d%s\n",pre,$1,post)}}' |\
tee worstatom_path.pdb



egrep "^CRYST" refmacout.pdb >! symme.pdb
cat worstatom.pdb >> symme.pdb
echo symgen $smallSG | pdbset xyzin symme.pdb xyzout symmates.pdb
wrap_into_cell.com doprotein=1 symmates.pdb outfile=worstatom_symmates.pdb






awk '/^ATOM|^HETAT/ && ! /HOH/{print substr($0,61,6)}' refmacout.pdb | histogram.awk | tee Bhist.txt



awk '/^ATOM|^HETAT/{id=substr($0,12,15);B=substr($0,61,6);sum[id]+=B;++count[id]} END{for(id in sum)printf("%.6g %d %s\n",sum[id],count[id],id)}' ../opt*/restraints_for*.pdb | sort -g | tee ../restraint_sums.txt

awk '/^ATOM|^HETAT/{id=substr($0,18,11);B=substr($0,61,6);sum[id]+=B;++count[id]} END{for(id in sum)printf("%.6g %d %s\n",sum[id],count[id],id)}' ../opt*/restraints_for*.pdb | sort -g | tee ../restraint_residue_sums.txt

awk '/^ATOM|^HETAT/{id=substr($0,31,24);B=substr($0,61,6);sum[id]+=B;++count[id]} END{for(id in sum)printf("%.6g %d %s\n",sum[id],count[id],id)}' ../opt*/restraints_for*.pdb | sort -g | tee ../restraint_xyz_sums.txt

cat  ../opt16/all_possible_refpoints.pdb ../restraint_xyz_sums.txt |\
awk '/^ATOM|^HETAT/{xyz=substr($0,31,24);id[xyz]=substr($0,12,15);next}\
  {xyz=substr($0,length($0)-23)} id[xyz]{print}' |\
tee ../restraint_id_sums.txt



rm -f flappers.txt
foreach itr ( `seq 1 1000` )

echo "$itr"
if(! -e restraints_for_${itr}.pdb) break

@ prev = ( $itr - 1 )
cat restraints_for_${prev}.pdb restraints_for_${itr}.pdb | egrep -v "HOH|EDO" | rmsd | grep once | tee -a flappers.txt

end
awk '{++seen[$0]} END{for(x in seen)print seen[x],x}' flappers.txt  | sort -g

awk '/^WARN/{print substr($0,10,15)}' flappers.txt |\
cat - current_restraints.pdb |\
awk 'length($0)<16{++sel[$0];next}\
  {id=substr($0,12,15)}\
  /^ATOM|^HETAT/ && sel[id]{print}' |\
sort -k1.61g



awk '$2=="="{printf("%s ",$1 $2 $3)} FNR>15{print FILENAME;nextfile}' `ls -1 ../opt*/runme1.log`

awk 'FNR==1{printf("%d ",substr(FILENAME,7)+0)} /^reffile/{print "";nextfile} $2=="="{printf("%s ",$1 $2 $3)} FNR>10{print "";nextfile}' ../opt*/runme1.log |\
 sort -g |\
tee ../optrun_options.txt

tail -n 1 ../opt*/fofc_Rplot.txt |\
awk '/^==/{o=substr($2,7)+0;getline;print o,$0}' |\
sort -g | tee ../optrun_vs_Rfac.txt

cat ../optrun_vs_Rfac.txt ../optrun_options.txt |\
awk '$3~/%/ && $2>10{R[$1]=$3;next}\
  {print R[$1],$0}'



rm opt_options.txt
foreach opt ( `ls -1d opt*/ | awk '{print substr($0,4)+0}' | sort -g` )

set options = `awk '/^reffile/{exit} $2=="="{n=split($1,w,"_");v="";for(i=1;i<=NF;++i){v=v substr(w[i],1,1)};print v $2 $3} FNR>10{exit}' opt${opt}/runme1.log | sort`

set R = `tail -n 1 opt${opt}/fofc_Rplot.txt | awk '{print $1,$2}'`

echo "$opt $R $options" | tee -a opt_options.txt

end
justify.awk opt_options.txt | sort -k2gr








foreach opt ( opt*/ )

cd $opt
pwd
if(! -s dens_vs_itr.txt) then
   set dens = `echo | mapdump mapin tempfile_rud_fofc_sigma.map | awk '/mum dens/{print $NF}'`
   echo "last $dens" | tee dens_vs_itr.txt
endif
if(! -s sorted_weights.txt) then
awk '/^ATOM|^HETAT/{print substr($0,61,6),substr($0,12,15)}' current_restraints.pdb | sort -gr >! sorted_weights.txt
endif

set infile = `ls -1rt amber_*.in | tail -n 1`
awk -F "=" '/temp/{print $NF+0 "K";exit}' $infile | tee temp.txt

set infile = `ls -1rt rest*.log | tail -n 1`
awk -F "=" '/fft/ && /B/{print $1,$NF}' $infile | tail -n 1 | tee fftB.txt

cd ..

end


rm opt_scores.txt
foreach opt ( opt*/ )

set maxw = `awk '{print $1;exit}' ${opt}/sorted_weights.txt`

set tight = `awk '$1>=10{++cnt} $1>0.06{sum+=$1} END{print sum+0,cnt+0}' ${opt}/sorted_weights.txt`

set dens = `tail -n 20 ${opt}/dens_vs_itr.txt | awk '{++n;sum1+=$2;sum2+=$3} END{print $1,sum1/n,sum2/n}'`

set R = `tail -n 1 ${opt}/refmac_Rplot.txt | awk '{print $2,$3}'`

set infile = `ls -1rt ${opt}/amber_*.in | tail -n 1`
set temp = `cat ${opt}/temp.txt`
set fftB = `awk '{print $NF}' ${opt}/fftB.txt | tail -n 1`

echo "$opt $R $maxw  $tight $dens  $temp  $fftB" | tee -a opt_scores.txt

end
justify.awk opt_scores.txt | sort -k2gr






ls -1rt amber_*.rst7 |\
head -n -1 |\
awk '{print "trajin",$0}' |\
tee these.in

echo "rmsd :WAT out rmsd.dat" | cat these.in - | cpptraj -p xtal.prmtop ; cat rmsd.dat 


rm -f restraints_dists_hists.txt
foreach rst ( `ls -1rt amber*.rst7` )

set prefix = `basename $rst .rst7`
set itr = `echo $prefix | awk -F "_" '{print $2}'`

if(! -e ${prefix}.pdb) then
rst2pdb_runme.com $rst > /dev/null
endif

cat ${prefix}.pdb |\
awk '/^ATOM|^HETAT/{print substr($0,1,60)"  0.00   "}' |\
rmsd -v debug=1 restraints_for_${itr}.pdb - |\
  awk '/moved/{print substr($0,25,10),substr($0,53,7),substr($0,1,17)}' |\
cat >! restraint_vs_deviate_${itr}.txt

cat restraint_vs_deviate_${itr}.txt |\
 histogram.awk -v bs=0.5 |\
awk -v itr=$itr '{print itr,$0}' |\
tee -a restraints_dists_hists.txt

end





rm -f changes_vs_itr.txt
set itrs = `ls -1rt amber_*.rst7 | grep -v equi | awk -F "_" '{print $NF+0}' | tail -n 1`
foreach itr ( `seq 2 $itrs` )
@ prev = ( $itr - 1 )

cat restraints_for_${prev}.pdb restraints_for_${itr}.pdb |\
 egrep -v "HOH|EDO" |\
 awk '! /^ATOM|^HETAT/{next} {id=substr($0,12,15);B=substr($0,61,6)+0}\
    B==0 && zeroB[id]{next} B==0{++zeroB[id]} {print}' |\
 rmsd >! rmsd.txt
head rmsd.txt
set changes = `grep WARN rmsd.txt | wc -l`
set pairs = `awk '/atom pairs found/{print $1;exit}' rmsd.txt`

egrep -h "HOH|EDO" restraints_for_${prev}.pdb restraints_for_${itr}.pdb |\
 awk '! /^ATOM|^HETAT/{next} {id=substr($0,12,15);B=substr($0,61,6)+0}\
    B==0 && zeroB[id]{next} B==0{++zeroB[id]} {print}' |\
 rmsd >! rmsdHOH.txt
set HOHchanges = `grep WARN rmsdHOH.txt | wc -l`
set HOHpairs = `awk '/atom pairs found/{print $1;exit}' rmsdHOH.txt`

echo "$itr   $changes $pairs   $HOHchanges $HOHpairs changes" | tee -a changes_vs_itr.txt

end

cat restraints_for_???.pdb |\
awk '/^ATOM/{xyz=substr($0,31,24);B=substr($0,61,6);\
  ++seen[xyz];sum[xyz]+=B;w[xyz,seen[xyz]]=B}\
 END{for(xyz in seen){\
   n=seen[xyz];avg=sum[xyz]/n;sumd=0;\
   for(i=1;i<=n;++i){sumd+=(w[xyz,i]-avg)^2};\
  print sqrt(sumd/n),avg,n,xyz}}' |\
cat >! restraints_rmsd.txt
 

rm -f sorted_restraints_vs_itr.txt
foreach pdb ( `ls -1rt restraints_for_*.pdb | tail -n 100 | tac` )
  set itr = `echo $pdb | awk -F "_" '{print $NF+0}' `
  awk '/^ATOM|^HETAT/{print substr($0,61,6)}' restraints_for_${itr}.pdb |\
  sort -gr |\
  awk -v itr=$itr 'NR%137==0{print itr,NR,$0}' |\
  tee -a sorted_restraints_vs_itr.txt
end


   echo 40 |\
   cat - randel_selections.txt current_restraints.pdb |\
   awk 'NR==1{sofar=$1}\
     $NF=="RANDEL"{split($0,w,"|");bin[w[2]]=w[1];next}\
     ! /^ATOM|^HETAT/{print;next}\
     {xyz=substr($0,31,24);B=substr($0,61,6)+0}\
     bin[xyz]>sofar{next}\
     B>900{print}' |\
   cat >! randel_restraints.pdb


foreach itr ( `seq 1 34` )

set laStage = amber_${itr}
set Stage = barometer

 cat ${laStage}.in |\
 awk '{space=substr($0,1,index($0,$1)-1)}\
  # just one cycle \
  / nstlim=/{split($1,w,"=");print space w[1] "=" 1 ",";next}\
  / ntpr=/{split($1,w,"=");print space w[1] "=" 1 ",";next}\
  / ntb=/{split($1,w,"=");print space "ntb=2,ntp=4,";next}\
  {print}' |\
 cat >! barometer.in

if(! -e resized_${itr}.parm7) then
rst2pdb_runme.com amber_${itr}.rst7 >! rst2pdb.log
mv resized.parm7 resized_${itr}.parm7
endif
set nwaters = `grep "O   HOH" amber_${itr}.pdb  | wc -l`
update_centroid_positions_runme.com restraints_for_${itr}.pdb topfile=resized_${itr}.parm7 > /dev/null

$pmemd -O -i ${Stage}.in -o ${Stage}.out \
   -p resized_${itr}.parm7 \
   -c ${laStage}.rst7 \
   -ref newref.crd \
   -r ${Stage}.rst7 \
   -x ${Stage}.nc \
   -inf ${Stage}.mdinfo

set press = `awk '/PRESS/ && p{n=$3;print $NF} /A V E R A/{++p} END{print P,r,n+0}' ${Stage}.out`
if( "$press" == "" || "$press" == "0" ) set press = "0 - -"

echo "$itr pressure was $press $nwaters" | tee -a pressure2_vs_itr.txt

end

