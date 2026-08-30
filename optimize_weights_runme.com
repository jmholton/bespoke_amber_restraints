#! /bin/tcsh -f
#
#  start with files:
#  ../refme.mtz
#  ../centroids/centroids_in_density.pdb
#  ${template_dir}/amberme.pdb
#  ../*.mol2 ../*.frcmod
#  tleap_stub.in
#  ${PHENIX}
#  ../xtal_properties.sourceme
#
#  iteratively optimize restraint weights based on fofc difference map
#
#   add the netfrc=0 trick by default?
#

# defaults for all variables - can be overridden with file: settings.sourceme

# properties of the data in the mtz file
set reso = 1.0
set smallSG = P212121
# supercell parameters
set super_mult = 1,1,1
# number of residues in one protein
set modulo  = 64

if(! -e xtal_properties.sourceme && -e ../xtal_properties.sourceme) then
  echo "using ../xtal_properties.sourceme"
  cp ../xtal_properties.sourceme .
endif
if(-e xtal_properties.sourceme) then
  echo "reading xtal_properties.sourceme"
  source xtal_properties.sourceme
else
  echo "WARNING: no xtal_properties.sourceme , defaulting to 1aho parameters"
  echo "confirm?"
  set in = ( $< )
endif

if(! -e xtal_properties.sourceme) then
  echo "WARNING: generating default xtal_properties.sourceme file"
  cat << EOF | tee xtal_properties.sourceme
set reso = $reso
set smallSG = $smallSG
# supercell parameters
set super_mult = $super_mult
# number of residues in one protein
set modulo  = $modulo
EOF
  echo "confirm?"
  set in = ( $< )
endif

# for map generation
set render_reso = mtz
set render_B = 10
# iterative adjustment of overall B factor
set render_B_min = 2
set render_B_adjust = 0.05
# for gemmi
set render_rate = 1.5
# file for storing all atomic B factors
set Bfac_file = Bfac.pdb
set Bfac_maxmod = 2
set Bfac_modmode = add
# B factors by LOCATION rather than by atom: a CCP4 map of B(x,y,z) that the
# structure-factor step samples at each atom's position.  Set to a file name to
# use a field you maintain yourself (Bfac_map_runme.com), or to "auto" to have
# nc2mtz build one from $Bfac_file each time.  Empty = the per-atom path.
set Bfac_map = ""
# map-field knobs, only used when Bfac_map is active (forwarded to nc2mtz).
#  sigma = Gaussian width (A) of the B field; must be wide enough to span how far
#          atoms move between frames, or displaced atoms drop into farB.
#  grid  = map sampling (A).
#  farB  = B assigned where no atom is nearby; keep modest (not 999) so coverage
#          misses degrade gracefully instead of going invisible.
set Bfac_map_sigma = 0.5
set Bfac_map_grid  = 0.5
set Bfac_map_farB  = 999
# restrict the field to these residues (default: water only).  Everything else
# keeps its per-atom B, so protein B is never smoothed into the high-B solvent.
# Set to "" to field all atoms, or to another residue name.
set Bfac_map_selresn = HOH
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
set maxitr = 999999
set delete_traj = 0

# adjust restraint weights based on difference map
set adjust_itr = 1
# adjust B factors based on difference map
set Badjust_itr = 1
# run refmac to optimize B factors
set refmac_itr = 0
# run phenix.refine to optimize B factors
set phenix_itr = 0
# re-discover reference points every so often
set repick_itr = 10
set repick_itr_ramp = none
# remove water restraints that are too close
set filter_itr = 1
set filter_dist = 1.8
# periodically delete restraints that are historically too highly challenged
set release_itr = 50
# periodically force CA atoms with possible restraints to have them
set reinin_itr = 0
# periodically force alignment reference atoms to have restraint weights
set minwt_itr = 0
# periodically discard restraints that are not in density
set rhocheck_itr = 0
set min_ref_rho = 0.5
# periodically delete restraints to atoms that have alt locs
set delete_altconf_itr = 0
# periodically just delete most-challenged restraints
set delete_badrest_itr = 0

# periodically reset all B factors
set Breset_itr = 0
# average B factors through bonds to neighboring atoms
# fraction to move into neighbors
set thrubond_avg_B = 0
# how much to favor larger values
set thrubond_avg_B_spread = 0
set thrubond_avg_B_ramp = none

# smooth the difference map
set fft_B = 0
set fft_B_ramp = none
# fft_B for the B-factor update path: an fft blur (in A^2) applied to the difference
# map that Bfac_update_diffmap.com probes, smoothing the noisy per-atom fofc signal so
# B factors stop random-walking at high Bfac_maxmod.  This is SEPARATE from fft_B above,
# which only blurs the fofc.map used for the weight/restraint update - so the B and
# weight regularizations can be tuned independently.  0 = no blur (unchanged behaviour).
set Bfac_fft_B = 0
set shan_B = auto
# criteria for statistical significance in difference peaks
set halfrho_pos = auto
set halfrho_neg = auto
set halfrho_ramp = none
# apply an overall scale factor to all weights every round
set weight_scaledown = 1
# apply a different scale factor to negative peaks
set weight_negscaledown = 0.95
# raise small weights to a power to drive them toward zero
set weight_power = 1.01
# criterion for a weak weight (multiples of kT)
set weight_power_kTmult = 1
# scale factor for converting "B factors" in restraint file to amber weights
set pdbscale = 0.01
# smallest amber weight to use
set cutoff_weight = 0.009
set cutoff_weight_ramp = none
# remember or forget refpoints that dip below cutoff weight
set cutoff_forget = 0
# apply a constant restraint weight to every atom, even if unspecified
set allatom_weight = 0
# starting value for newly created restraints
set weight0 = 1
# minimum weight to apply to CA atoms used for alignment 
set min_CA_weight = 0
# minimum weight to apply to any atoms used for alignment 
set min_align_weight = 0.01
# minimum weight to apply to ligands 
set min_lig_weight = 0.01
# weight to apply to CA atoms that have wandered from their reference points
set errant_CA_weight = 2
# upon release_itr, also release restraints on symmate atoms
set release_radius = 1
# new weight to give to "released" restarints
set release_weight = 0.01
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
#set chiral_weight = 10
set chiral_weight = 0
set chiral_weight_ramp = none
# amber restraints on peptide bond dihedrals
#set omega_weight = 50
set omega_weight = 0
set omega_weight_ramp = none
# exit with an error as soon as any cis peptide appears (rather than fighting it with omega weights)
set exit_on_cis = 0
# periodically re-map water names
set remap_itr = 0
# periodically wrap non-restrained atoms to inside the supercell
set wrap_itr = 2
# periodically re-center the system on restrained atoms
set align_itr = 1
set align_nstlim = 0
# type of atoms to use for alignment

set align_target = centroids
# seems to avoid need for re-centering
set netfrc = 0

# periodically move water molecules from bad Fo-Fc density to good Fo-Fc density
set teleport_itr = 1
set teleport_waters = 10
set teleport_weight = 5
set teleport_mindist = 2.0
set teleport_minrho = 2.0
# periodically add or remove waters, depending on pressure and void size
set hydrate_itr = 1
set minvoid = 60
set pressure_deadband = 10
set hydrate = voids
set dehydrate = pressure
set pressure_scale = auto
set pressure_scale_ramp = none
# average pressure over previous runs
set pressure_avglast = 1
set pressure_avglast_ramp = none
set void_scale = 1
# parametes for AddToBox
set add_radius = 2.4
# always make sure water count stays the same
set water_lock = 0
# always make sure water count always changes by at least one
set water_dither = 0
# run pressure measurement after normal amber run with all restraints removed
set barometer_cycles = 10000

# average electron density over this and previous amber runs
set avglast = 1
set avglast_ramp = none
# equilibrate the system for a bit before doing a produciton run
set equi_ns = 0.2
set equi_dt = 0.0005
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
set equi_gamma = same
set settle_slowdown = 10
set barostat = 0
set kT = 0.6

# maximum allowed restraint weight
set max_weight = 9.9999
# maximum scale factor for updating restraint weights
set max_mult = 2.0
set max_mult_ramp = none
# average restraint weights through bonds to neighboring atoms
# fraction to move into neighbors
set thrubond_avg_weight = 0
# how much to favor larger values
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

# need better scratch-finding logic perhaps
set scratch = /scratch/${USER}/opt_`hostname -s`_$$_
mkdir -p /scratch/${USER}

# commands to submit jobs to GPU or CPU cluster queues
set sruncpu = "srun --partition=refmac --exclude=crush18"
set pmemd = "srun --partition=gpu --gres=gpu:1 pmemd.cuda_SPFP"

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
      if("$key" == "min_weight") set cutoff_weight = "$Val"
      if("$key" == "align_cyc") set align_nstlim = "$Val"
      if("$key" == "weight_scale") set weight_scaledown = "$Val"
      if("$key" == "weight_negscale") set weight_negscaledown = "$Val"
    else
      # no equal sign
    endif
    if("$key" == "debug") set debug = "1"
end

# shorthand for temporary file
set t = $tempfile

# find a place for scratch
set base = `basename $scratch`
foreach scratch ( $scratch /data/${USER}/scratch/${base}/ /scratch/${USER}/${base}/ /dev/shm/${USER}/${base}/ /tmp/${USER}/${base}/ )
  mkdir -p ${scratch} >& /dev/null
  if(-w ${scratch}) then
    echo "scratch = $scratch is writable"
    break
  endif
end
echo "scratch = $scratch"

set md5 = `md5sum $0 | awk '{print $1}'`
echo "running as $0 ($md5)"

# is there another job running?
set pwd = `pwd`
set mypid = $$
set myname = `basename $0 | awk '/[0-9.]/{$0=substr($0,1,match($0,/[0-9.]/)-1)} {print}'`
set alone = 0
while ( ! $alone )
  
  set pidirs = `ps -fu $USER | tee ${t}debug1.log | egrep "$myname" | egrep -v " egrep -v | grep -E |^UID" | awk -v pid=$mypid '$2!=pid && ! ( / srun / && /optimize_weights/ ) {print "/proc/"$2"/cwd"}'`

  set otherpid = `ls -l $pidirs |& tee ${t}debug2.log | awk -v pwd="$pwd" '$NF==pwd{print}' | awk -F "/" '{print $3}'`
  if( "$otherpid" == "") then
    set alone = 1
    break
  endif

  set sleepids = `awk '/sleep/{print $2}' ${t}debug1.log`
  if( "$sleepids" != "" ) then
    echo "kill sleep at $sleepids ? "
    sleep `echo $$ | awk '{srand($1);print 3+rand()}'`
  endif
#  ps -flea | grep $otherpid
  echo "other job running here: $otherpid , we are $mypid "
  # try to make the other job exit
  touch exit
  rm -f ${t}debug1.log ${t}debug2.log
  sleep 60
end
rm -f ${t}debug1.log ${t}debug2.log

if( $render_reso == "mtz" || "$render_reso" == "auto" ) then
  set render_reso = $reso
endif 

set test = `echo $min_lig_weight $min_CA_weight $min_align_weight | awk '{print $1+$2+$3}'`
if( "$minwt_itr" == "0" && "$test" != "0" ) then
  echo "WARNING: setting minwt_itr = 1 because weights: lig $min_lig_weight CA $min_CA_weight align $min_align_weight"
  set minwt_itr = 1
endif
set test = `echo $max_mult | awk '{print $1+0}'`
if( "$adjust_itr" == "0" && "$test" != "1" ) then
  echo "WARNING: setting adjust_itr = 1 because max_mult = $max_mult (!= 1)"
  set adjust_itr = 1
endif
set test = `echo $Bfac_maxmod $Bfac_modmode | awk '$2=="add"{print ( $1+0 > 0 )} $2=="mult"{print ( $1*1 > 1 )}'`
if ( $test == "" ) set test = 0
if( "$Badjust_itr" == "0" && $test ) then
  echo "WARNING: setting Badjust_itr = 1 because Bfac_maxmod = $Bfac_maxmod ( $Bfac_modmode mode )"
  set Badjust_itr = 1
endif

# start the ramps
if( "$fft_B_ramp" != "none" ) set fft_B = `echo $fft_B_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$cutoff_weight_ramp" != "none" ) set cutoff_weight = `echo $cutoff_weight_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$repick_maxdist_ramp" != "none" ) set repick_maxdist = `echo $repick_maxdist_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$repick_hohscale_ramp" != "none" ) set repick_hohscale = `echo $repick_hohscale_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$avglast_ramp" != "none" ) set avglast = `echo $avglast_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$pressure_avglast_ramp" != "none" ) set pressure_avglast = `echo $pressure_avglast_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$pressure_scale_ramp" != "none" ) set pressure_scale = `echo $pressure_scale_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$equi_ns_ramp" != "none" ) set equi_ns = `echo $equi_ns_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$prod_ns_ramp" != "none" ) set prod_ns = `echo $prod_ns_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$max_mult_ramp" != "none" ) set max_mult = `echo $max_mult_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$thrubond_avg_weight_ramp" != "none" ) set thrubond_avg_weight = `echo $thrubond_avg_weight_ramp | awk '$1~/^[0-9]/{print $1+0}'`
if( "$thrubond_avg_B_ramp" != "none" ) set thrubond_avg_B = `echo $thrubond_avg_B_ramp | awk '$1~/^[0-9]/{print $1+0}'`

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
blim 2 999
damp 0 0.5
weigh matrix 1
ncyc 5
make link Y
make hydr Y
make hout Y
EOF


endif

if( "$ligands" == "auto" ) then
    set ligands = ""
    echo "WARNING: ligands=auto means ligands not set up"
  endif
endif

set ligcifs = `echo $ligands | awk '{for(i=1;i<=NF;++i) print $i ".cif"}'`
foreach cif ( $ligcifs )
  if(! -e $cif ) then
    echo "cp ../ligands/$cif ."
    cp ../ligands/$cif .
  endif
end


# in case we need to run tleap
if(! -e tleap_stub.in) then
  cp -p ${template_dir}/*.mol2 .
  cp -p ${template_dir}/*.frcmod .
#  cp -p ${PHENIX}/modules/amber_library/m/MSE.* .
  cp ${template_dir}/protonation.txt protonation.txt
  cp ${template_dir}/tleap_stub.in .
endif

if(! -e amberme.pdb) then
  echo "cp ${template_dir}/amberme.pdb amberme.pdb"
  cp ${template_dir}/amberme.pdb amberme.pdb
endif

# structure factors for best-phased reference map
if(! -e reference.mtz ) then
  echo "cp ${template_dir}/reference.mtz reference.mtz"
  cp ${template_dir}/reference.mtz .
endif
if(! -e reference.mtz ) then
  echo "cp  ../centroids/reference0.mtz reference.mtz"
  cp  ../centroids/reference0.mtz reference.mtz
endif

# master list of points in space to use as restraint reference points
if(! -e all_possible_refpoints.pdb ) then
   # sometimes gets edited
#  echo "cp ${template_dir}/reference.mtz reference.mtz"
#  cp  cp ${template_dir}/reference.mtz .
endif
# start fresh
if(! -e all_possible_refpoints.pdb ) then
  echo "generating all_possible_refpoints.pdb from ../centroids/centroids_in_density.pdb"
  set B0 = `echo $weight0 $pdbscale | awk '{print $1/$2}'`
  cat ../centroids/centroids_in_density.pdb |\
  awk -v B0=$B0 '! /^ATOM|^HETAT/{print;next}\
     {pre=substr($0,1,60);post=substr($0,67);rho=$NF;\
      B=B0*rho}\
     B>B0{B=B0}\
     {printf("%s%6.2f%s\n",pre,B,post)}' |\
  cat >! all_possible_refpoints.pdb

  if( "$min_lig_weight" != "0" ) then
      if(! -e this.pdb ) cp ../centroids/fulllength_super.pdb this.pdb
      mv all_possible_refpoints.pdb refpoints_in_density.pdb
      set liglist = `echo $ligands | awk '{gsub(" ",",");print}'`
      set Blig = `echo $min_lig_weight $pdbscale | awk '{print $1/$2}'`
      awk '/^CRYST|^ATOM|^HETAT/' this.pdb |\
      filter_pdb.awk -v only=ligand -v ligands="$liglist" -v skip=H |\
       reformatpdb.awk -v BFAC=$B0 >! lig_restraints.pdb

      combine_pdbs_runme.com lig_restraints.pdb refpoints_in_density.pdb this.pdb \
        outfile=all_possible_refpoints.pdb
  endif
endif

if( "$randel_itr" != "0" ) then
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

foreach stage ( `awk '/^Stages /{print}' leap2amber_${itr}.log` )
  if(-e ${stage}.nc) set Stage = $stage
end

endif


if(! -e ref.crd) then
   echo "generating new ref.crd "
   update_centroid_positions_runme.com current_restraints.pdb  outfile=ref.crd 
endif


if (! $?Stage ) then
  echo "no previous Stage defined"
  set itr = `ls -1Lrt amber_*.rst7 |& egrep -v "equi|settle|unwrap" | awk -F "_" '/amber_/ && NF==2{print $2+0}' | tail -n 1`
  if( "$itr" != "" ) set Stage = amber_${itr}
  #echo "no previous amber runs"
endif

if( ! $?Stage ) then
  echo "still no previous Stage defined"
  set log = `ls -1Lrt *amber_* |& grep amber_ | egrep -v "equi|settle" | tail -n 1`
  set i = `echo $log | awk -F "_" 'NF==2{print $2+0}' | tail -n 1`
  if("$i" == "") set i = 0
  if(-e amber_${i}.nc) then
    set Stage = amber_${i}
  else
    set i = 0
  endif

  if(-e leap2amber_${i}.log) then
    set itr = "$i"
    echo "checking leap2ampber_${itr}.log"
    set Stage = `awk '/^Stages /{for(i=3;i<=NF;++i)print $i}' leap2amber_${itr}.log | grep -v Min | tail -n 1`
    if( "$Stage" == "" ) set Stage = `ls -1rt *.nc |& egrep -v "equi|settle|unwrap" | awk -F "." '{print $1}' | tail -n 1`
  endif
endif
if( ! $?Stage && $?i ) then
  if( -e amber_${i}.nc && -e amber_${i}.rst7 ) then
    set Stage = amber_${i}
    echo "setting Stage = $Stage"
  endif
endif
if( ! $?Stage ) then
   set BAD = "unable to determine stage of optimization."
   goto exit
endif


if(! -e orignames.pdb ) then 
    echo "cp ${template_dir}/orignames.pdb ."
    cp ${template_dir}/orignames.pdb .
endif
if(! -e Bfac.pdb ) then 
   echo "generating Bfac.pdb from orignames.pdb"
#   cp  orignames.pdb Bfac.pdb 
   cat orignames.pdb |\
   awk '! /^ATOM|^HETAT/{print;next} \
     {TYP=substr($0,18,3)}\
     TYP!="HOH"{print;next}\
     {B=substr($0,61,6)+0}\
     B>maxB{maxB=B} B<maxB{B=maxB}\
     {printf("%s%6.2f%s\n",substr($0,1,60),B,substr($0,67))}' |\
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
if(! $?laStage) set laStage = "$Stage"

# test and find a topfile that works with the current stage
echo "${Stage}.rst7 -> this.pdb"
rst2pdb_runme.com ${Stage}.rst7 this.pdb >! rst2pdb_${Stage}.log
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

  if( "$min_lig_weight" != "0" ) then
      mv all_possible_refpoints.pdb refpoints_in_density.pdb
      set liglist = `echo $ligands | awk '{gsub(" ",",");print}'`
      set Blig = `echo $min_lig_weight $pdbscale | awk '{print $1/$2}'`
      awk '/^CRYST|^ATOM|^HETAT/' this.pdb |\
      filter_pdb.awk -v only=ligand -v ligands="$liglist" -v skip=H |\
       reformatpdb.awk -v BFAC=$B0 >! lig_restraints.pdb

      combine_pdbs_runme.com lig_restraints.pdb refpoints_in_density.pdb this.pdb \
        outfile=all_possible_refpoints.pdb
  endif
endif

echo "checking for restraint bombs"
rm -f defused_restraints.pdb
restraint_bomb_detector.com ${Stage}.rst7 \
    refpoints=current_restraints.pdb \
    outfile=defused_restraints.pdb min_weight=$cutoff_weight \
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


while ( $itr < $maxitr )

# allow run-time edits to variables
if(-e settings.sourceme) then
    cat settings.sourceme
    source settings.sourceme
endif
rm -f ./exit >& /dev/null

set prev = "$itr"
@ itr = ( $itr + 1 )
set lalaStage = "$laStage"
set laStage = $Stage
set Stage = amber_${itr}

if(! -e ${laStage}.in) then
  set BAD = "no previous input file: ${laStage}.in"
  goto exit
endif
set trajectory = ${laStage}.nc
if(! -e "$trajectory" ) then
    # maybe last run crashed early
    echo "WARNING: no trajectory: $trajectory "
    set trajectory = ${laStage}_1.nc
    echo "WARNING: trying $trajectory "
endif
if(! -e "$trajectory" ) then
  echo "WARNING: "
  echo "WARNING: no trajectory: $trajectory "
  echo "WARNING: skipping ahead"
  echo "WARNING: "
  goto wrap
endif

# measure longest diffusion 
traj_maxd_runme.com $trajectory >! traj_maxd.log
touch maxd_vs_itr.txt
echo -n "itr rms max frames atom: "
echo -n "$prev " | tee -a maxd_vs_itr.txt
tail -n 1 traj_maxd.log | tee -a maxd_vs_itr.txt


# check for out-of-date trajectory
if(-e trajectory/md.1.pdb) then
  set test = `ls -1rt $trajectory trajectory/md.1.pdb | head -n 1 | grep md.1.pdb | wc -l`
  if( $test ) then
    echo "WARNING: trajectory/md.1.pdb is older than $trajectory "
    echo "maybe: rm -f trajectory/md.1.pdb"
  endif
endif


wrap:
# see if we need a wrap - unprintable atoms
echo "extracting xyz from ${laStage}.rst7 -> this.pdb"
rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log
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

set trigger = `echo $itr $repick_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) set need_wrap_now = 1

set trigger = `echo $itr $wrap_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $need_wrap_now || $trigger ) then
  echo "wrapping ${laStage}.rst7"

  # defaults to "all_possible_refpoints.pdb and current_restraints.pdb"
  rst_wrap.com ${laStage}.rst7 | tee wrap_${itr}.log | egrep "rmsd_wrap|maxd_wrap"

  cp ${laStage}.rst7 ${laStage}_unwrapped.rst7
  mv wrapped.rst7 ${laStage}.rst7
  if(-e wraps.txt) mv wraps.txt wraps_${itr}.txt
  # keep everything else same for laStage

  # update the this.pdb
  rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log

#  if(-e ${laStage}.prmtop) cp ${laStage}.prmtop wrapped.prmtop
#  set laStage = wrapped

  echo "wrapping any errant refpoints..."
  rm -f defused_restraints.pdb
  restraint_bomb_detector.com ${laStage}.rst7 \
    refpoints=current_restraints.pdb \
    outfile=defused_restraints.pdb min_weight=$cutoff_weight \
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

filter_pdb.awk -v skip=EP this.pdb >! refme.pdb
#cp refme.pdb refme0.pdb


# check for inverted chiral centers and cis peptides
echo "checking chiral centers and peptide bonds"
set liglist = `echo $ligands | awk '{gsub(" ",",");print}'`
filter_pdb.awk -v only=protein,ligand -v ligand=$liglist -v skip=H,water refme.pdb |\
awk '{print substr($0,1,80)}' >! protein.pdb
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
    echo "avglast = $avglast"
endif
if( "$pressure_avglast_ramp" != "none" ) then
    set params = `echo $pressure_avglast_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $pressure_avglast | awk 'NF>2 && $3>$2 && $NF<=$2{print $2;exit} NF>2 && $3<$2 && $NF>=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 1 | awk '$1<$2{$1=$2} {print $1}'`
    set pressure_avglast = $value
    echo "pressure_avglast = $pressure_avglast"
endif
if( "$pressure_scale_ramp" != "none" ) then
    set params = `echo $pressure_scale_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $pressure_scale | awk 'NF>2 && $3>$2 && $NF<=$2{print $2;exit} NF>2 && $3<$2 && $NF>=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0 | awk '$1<$2{$1=$2} {print $1}'`
    set pressure_scale = $value
    echo "pressure_scale = $pressure_scale"
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

if( "$thrubond_avg_B_ramp" != "none" ) then
    set params = `echo $thrubond_avg_B_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $thrubond_avg_B | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.1 | awk '$1<$2{$1=$2} {print $1}'`
    set thrubond_avg_B = $value
    echo "thrubond_avg_B = $thrubond_avg_B"
endif

if( "$cutoff_weight_ramp" != "none" ) then
    set params = `echo $cutoff_weight_ramp | awk -F "[: _-]" '{print $3,$2,$1}'`
    set value = `echo $params $cutoff_weight | awk 'NF>2 && $NF<=$2{print $2;exit} $1~/x$/{print $NF*$1;exit} $(NF-1)~/x$/{print $NF*$(NF-1);exit} {print $NF-$1}'`
    set value = `echo $value 0.0001 | awk '$1<$2{$1=$2} {print $1}'`
    set cutoff_weight = $value
    echo "cutoff_weight = $cutoff_weight"
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


set trigger = `echo $itr $align_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger || ! -e align_ref.pdb ) then
  # align rst7 and nc files to reference
  touch align_${prev}.log

  set threshlist = ( med )
  if( "$align_target" == "restraints" ) then
    echo "aligning ${laStage} to current restraints"
    cp ref.crd align_ref.crd
    cp current_restraints.pdb align_ref.pdb
    set threshlist = ( 0.6 med )
  endif
  if( "$align_target" == "CA" ) then
    set threshlist = ( 1 med )
    echo "aligning $laStage to high-density CA atoms"
    if(-e align_ref.pdb) then
      echo "re-using align_ref.pdb"
    else
      grep "  CA  " all_possible_refpoints.pdb >! temp_ref.pdb
      rholabel_runme.com reference.mtz temp_ref.pdb mtzlabel=Fref >> align_${prev}.log
      cat rholabeled.pdb |\
      awk -v s=$pdbscale '! /^ATOM/{next} substr($0,12,6)!="  CA  "{next}\
       {pre=substr($0,1,60);B=substr($0,61,6)+0;rho=$NF/s}\
        B<0.01{next} rho<0.5{next} rho>999.99{rho=999.99}\
        {printf("%s%6.2f%20s\n",pre,rho,"")}' |\
      cat >! align_ref.pdb
    endif
  endif
  if( "$align_target" == "centroids" ) then
    set threshlist = ( 1 med )
    echo "aligning $laStage to high-density atoms"
    if(-s align_ref.pdb) then
      echo "re-using align_ref.pdb"
    else
      echo "generating new align_ref.pdb from all_possible_refpoints.pdb"
      filter_pdb.awk -v only=protein all_possible_refpoints.pdb |\
      reformatpdb.awk -v CONF=" " |\
      rmsd2B |\
      awk '/^CRYST/ || substr($0,61,6)+0<12' >! ${t}_ref.pdb
      rholabel_runme.com reference.mtz ${t}_ref.pdb mtzlabel=Fref >> align_${prev}.log
      if( $status ) then
        set BAD = "map probe with rholabel failed"
        goto exit
      endif
      cat rholabeled.pdb |\
      awk -v s=$pdbscale '! /^ATOM/{next}\
       {pre=substr($0,1,60);rho=$NF/s}\
        rho>999.99{rho=999.99}\
        {printf("%s%6.2f%20s\n",pre,rho,"")}' |\
      cat >! align_ref.pdb
    endif
  endif
  if( "$align_target" == "Bfactors" ) then
    set threshlist = ( med )
    echo "aligning $laStage to low-B-factor CA atoms"
    grep "  CA  " Bfac.pdb >! temp_B.pdb
    grep "  CA  " all_possible_refpoints.pdb >! temp_ref.pdb
    combine_pdbs_runme.com temp_B.pdb temp_ref.pdb >> align_${prev}.log
    cat new.pdb |\
    awk '! /^ATOM/{next} substr($0,12,6)!="  CA  "{next}\
     {pre=substr($0,1,60);B=substr($0,61,6);\
      invB=1000.0/(B+1e-6)} invB>999.99{invB=999.99}\
      {printf("%s%6.2f%20s\n",pre,invB,"")}' |\
    cat >! align_ref.pdb
    # needs to be re-generated
    rm align_ref.crd
  endif
  if(! -s align_ref.pdb) then
    set BAD = "failed to find/create alignment reference align_ref.pdb"
    goto exit
  endif 
  # now threshold the weights assigned to decide which atoms to use
  cat align_ref.pdb |\
  awk 'substr($0,61,6)+0>0{print substr($0,1,80),"     SEL"}' |\
  cat >! selected.pdb
  combine_pdbs_runme.com selected.pdb this.pdb printref=1 >> align_${prev}.log
  egrep "^ATOM|^HETAT" new.pdb |\
  awk -v pdbscale=$pdbscale '{++n} $NF=="SEL"{print n,substr($0,61,6)*pdbscale}' |\
  cat >! align_atom_weights.txt
  foreach w ( $threshlist )
    if( $w == med ) then
      set h = `cat align_atom_weights.txt | wc -l | awk '{print int($1/2)}'`
      set w = `awk '{print $2}' align_atom_weights.txt | sort -g | head -n $h | tail -n 1`
    endif
    awk -v w=$w '$2>=w' align_atom_weights.txt >! alignment_atnums.txt
    cat alignment_atnums.txt |\
    awk 'NR==1{s=e=$1;next}\
          $1==e+1{e=$1;next} \
          {print s"-"e;s=e=$1}\
          END{print s"-"e}' |\
    awk -F "-" '$1==$2{print $1;next} {print}' >! align_ranges.txt 
    set align_ranges = `cat align_ranges.txt `
    set align_mask = `echo $align_ranges | awk '{gsub(" ",",");print}'`
    if( "$align_mask" != "" ) then
      break
    endif
  end
  if( 1 || ! -e align_ref.crd ) then
    echo "generating align_ref.crd"
    update_centroid_positions_runme.com \
         pdbfile=align_ref.pdb \
         outfile=align_ref.crd  >> align_${prev}.log
  endif
  if( "$align_mask" == "" ) then
    set BAD = "unable to assing alignment mask"
    goto exit
  endif

  rm -f aligned.rst7 >& /dev/null
  cpptraj -p xtal.prmtop -y ${laStage}.rst7 -c align_ref.crd << EOF >> align_${prev}.log
  rmsd rmsd reference norotate @$align_mask out rmsd.txt savevectors combined vecsout vecsout.txt
  trajout aligned.rst7
EOF
  cat rmsd.txt vecsout.txt >> align_${prev}.log

  if(-e "$trajectory") then
    rm -f aligned.nc >& /dev/null
    cpptraj -p xtal.prmtop -y $trajectory -c align_ref.crd << EOF >> align_${prev}.log
    rmsd rmsd reference norotate @$align_mask out rmsd.txt savevectors combined vecsout vecsout.txt
    trajout aligned.nc
EOF
    cat rmsd.txt vecsout.txt >> align_${prev}.log

    set trajectory = aligned.nc
  endif
  cp ${laStage}.in aligned.in
  set lalaStage = $laStage
  set laStage = aligned

  echo -n "final shift: "
  set drift = `tail -n 1 vecsout.txt | awk '{print $2,$3,$4,"(",sqrt($2*$2+$3*$3+$4*$4),")"}'`
  echo "$prev $drift" | tee -a alignment_drift_vs_itr.txt

  # give a warning for big drift?
  set test = `echo $drift | awk '{print ( $(NF-1) > 0.1 ) }'`
  if( $test ) then
    echo "WARNING: drift > 0.1 A detected! consider lowering align_nstlim or increasing min_align_weight"
  endif

  # update the this and noEP pdb files
  rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log
  filter_pdb.awk -v skip=EP this.pdb >! refme.pdb

endif


if(! -e "$trajectory" ) then
  echo "WARNING: "
  echo "WARNING: no trajectory: $trajectory "
  echo "WARNING: skipping ahead"
  echo "WARNING: "
  rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log
  filter_pdb.awk -v skip=EP this.pdb >! refme.pdb
  goto skipfofc
endif


# convert nc file to an mtz
set nc2log = nc2mtz_${prev}.log
set needtraj = 0
if( ! -e avg_${prev}.mtz ) set needtraj = 1
if( ! -e trajectory/md.1.pdb && $adjust_itr != 0 ) set needtraj = 1
if( ! -e trajectory/md.1.pdb && $Badjust_itr != 0 ) set needtraj = 1

if(-e "$trajectory" && $needtraj ) then
  set minB = 1
  set maxB = 999.99
  if( "$Bfac_file" == "rmsd2B" ) then
    set minB = `echo $rmsd2B_range | awk -F "-" '{print $1}'`
    set maxB = `echo $rmsd2B_range | awk -F "-" '{print $2}'`
  endif
  set nc2mtz_extraopt = ""
  if( -e avg_${prev}.mtz ) set nc2mtz_extraopt = "domaps=0"
  echo "nc2mtz gemmi $trajectory  ${render_reso}A $Bfac_file "
  set nc2mtz_Bmap = ""
  if( "$Bfac_map" != "" ) set nc2mtz_Bmap = "Bfac_map=$Bfac_map Bfac_map_sigma=$Bfac_map_sigma Bfac_map_grid=$Bfac_map_grid Bfac_map_farB=$Bfac_map_farB Bfac_map_selresn=$Bfac_map_selresn"
  nc2mtz_gemmi.com $smallSG super_mult=$super_mult $trajectory \
    reso=$render_reso B=$render_B \
    Bfac_file=$Bfac_file minB=$minB maxB=$maxB $nc2mtz_Bmap \
    keeptraj=1 wrap=1 rate=$render_rate \
    addmtzs=1 \
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

set int_avglast = `echo $avglast | awk '{print int($1)}'`
echo "averaging last ${avglast} itrs"
ls -1rt avg_*.mtz | awk -F "[_.]" '$2!~/[a-z]/' | tail -n $int_avglast >! latest_avgs.txt
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
  echo "$prev $sgstats" | tee -a fofc_Rplot_grid.txt
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

awk '/RFAC  RF/,/TOTALS/{print}' diff_details.log |\
 tee R_table.txt |\
 awk 'NF>15 && $7+0>0{print sqrt(1/$1),$7}' |\
 cat >! ${t}Rbins.txt
set bestR = `awk '{print $2}' ${t}Rbins.txt | sort -g | head -n 1`
set Rthresh = `echo $bestR 0.20 | awk '{print sqrt($1*$1+$2*$2)}'`
set bestRrange = `awk -v t=$Rthresh '$2<1.0*t{print}' ${t}Rbins.txt | awk 'NR==1{print $1} END{print $1}'`

set sitstats = `awk '/scale=/{s=$2;B=$4} /^TOTAL/{CC=$NF} /correct F:/{R=$7} END{print R,s,B,CC}' $nc2log`
echo "$prev $sitstats $bestR $bestRrange" | tee -a fofc_Rplot.txt

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
     set head = `cat this.pdb | wc -l`
     head -n $head newB.pdb >! ${t}newB.pdb
     set pegmin = `awk '/^ATOM|^HETAT/ && substr($0,61,6)+0<=2' ${t}newB.pdb | wc -l`
     set pegmax = `awk '/^ATOM|^HETAT/ && substr($0,61,6)+0>=999.99' ${t}newB.pdb | wc -l`
     echo "adding $deltaB to all B factors in $Bfac_file ($pegmin at min, $pegmax at max)"
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

if( ! -e trajectory/md.1.pdb ) then
  echo "WARNING: no trajectory."
  goto skipfofc
endif

set trigger = `echo $itr $adjust_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
  echo "adjusting restraints based on fofc map"
else
  goto skipadjust
endif

# use difference map to update restraints
rm -f new_restraints.pdb
set max_wB = `echo $max_weight $pdbscale | awk '{print $1/$2}'`
set rud_try = 0
set rud_dbg = ""
rud_retry:
@ rud_try ++
restraintlist_update_diffmap.com fofc.map \
  refmap=reference.map \
  trajectory=trajectory/ \
  tempfile=${scratch}/rud_ \
  refme.pdb modulo=$modulo \
  overall_scale=$weight_scaledown \
  negative_scale=$weight_negscaledown \
  max_weight=$max_wB max_mult=$max_mult \
  ambig_same_weight=$ambig_same_weight \
  halfrho_pos=$halfrho_pos halfrho_neg=$halfrho_neg \
  refpointspdb=current_restraints.pdb $rud_dbg \
  outmults=sorted_mults_${itr}.txt \
  outfile=new_restraints.pdb >&! restraint_update_${itr}.log
if( $status || ! -e new_restraints.pdb) then
  # Usually a transient cluster hiccup: some of the parallel srun map-peek jobs
  # exit 9, so the peek count comes up an atom short ("ref and xyz do not match")
  # and it bails.  Re-running clears it - one retry (in debug mode) before giving
  # up; a second debug repeat rarely teaches us anything the first one didn't.
  if( $rud_try < 2 ) then
    echo "restraintlist_update_diffmap failed (likely transient srun/map-peek) - retry #$rud_try in debug mode (-p debug partition, temp files kept for post-mortem)"
    rm -f new_restraints.pdb
    set rud_dbg = "debug=1"
    sleep 10
    goto rud_retry
  endif
  set BAD = "restraintlist update diffmap failed after $rud_try tries"
  goto exit
endif
awk 'NF<=2 || NF==3 && $2=="=" || /^[\[-]/ || /\/dev\/shm|^clearing|\/scratch\//{next} /^srun/ || $4=="srun"{next} {print}' restraint_update_${itr}.log

cat current_restraints.pdb new_restraints.pdb |\
awk '/^ATOM|^HETAT/{id=substr($0,12,19);B=substr($0,61,6)+0}\
  lastB[id]+0>0{print id,B/lastB[id],"ACTUAL"} {lastB[id]=B}' |\
cat >! actual_mults.txt

cp new_restraints.pdb current_restraints.pdb
cp new_restraints.pdb rud_restraints.pdb

skipadjust:


set trigger = `echo $itr $Badjust_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
  echo "adjusting B factors based on difference map"
else
  goto skipBadjust
endif

# now update B factors for map calculation
if(-e premod_Bfac_${itr}.pdb) then
   echo "assuming B factors in premod_Bfac_${itr}.pdb are for $trajectory "
   cp premod_Bfac_${itr}.pdb Bfac.pdb
endif
rm -f new_Bfac.pdb
Bfac_update_diffmap.com Bfac.pdb modulo=$modulo mtzfile=cootme.mtz \
  max_mod=$Bfac_maxmod mod_mode=$Bfac_modmode fft_B=$Bfac_fft_B \
  tempfile=${scratch}/Bud_ >&! Bfac_update_${itr}.log
if( $status || ! -e new_Bfac.pdb) then
  set BAD = "Bfac update diffmap failed"
  goto exit
endif

awk 'NF<=2 || NF==3 && $2=="=" || /^[\[-]/ || /\/dev\/shm|^clearing|\/scratch\//{next} {print}' Bfac_update_${itr}.log

cp Bfac.pdb premod_Bfac_${itr}.pdb
cp new_Bfac.pdb Bfac.pdb
cp Bfac.pdb Bfac_${itr}.pdb
cp sorted_Bmods.txt sorted_Bmods_${itr}.txt

if( "$thrubond_avg_B" != "0" ) then
  echo "averaging B factors through bonds with avgfac=$thrubond_avg_B and spread=$thrubond_avg_B_spread"
  cp Bfac.pdb presmooth_Bfac_for_${itr}.pdb
  thrubond_avgB_runme.com Bfac.pdb debug=$debug \
    avgfac=$thrubond_avg_B spread=$thrubond_avg_B_spread \
    outfile=smoothB.pdb >> Bfac_update_${itr}.log
  cp smoothB.pdb Bfac.pdb
endif


skipBadjust:

# re-generate this.pdb and refme.pdb with new B factors
rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log
filter_pdb.awk -v skip=EP this.pdb >! refme.pdb

touch Bhist_vs_itr.txt
awk '/^ATOM|^HETAT/{print substr($0,61,6)}' refme.pdb |\
 histogram.awk -v bs=1 |\
awk -v itr=$itr '{print itr,$0}' |\
cat >> Bhist_vs_itr.txt
echo "" >> Bhist_vs_itr.txt

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


set trigger = `echo $itr $refmac_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
  echo "running refmac to refine B factors"
else
  goto skiprefmac
endif

# now use refmac to refine B factors - better than previous?
rm refmacout.pdb >& /dev/null
converge_refmac.com refme.mtz refme.pdb trials=1 $ligcifs append nosalvage >&! refmac_${itr}.log
if( $status || ! -e refmacout.pdb) then
  set BAD = "refmac failed"
  goto exit
endif

cp Bfac.pdb prerefmac_Bfac_${itr}.pdb
set lastB = `tail refmacout.pdb | awk '/^ATOM|^HETAT/{print substr($0,61,6)}' | sort -gr | head -n 1`
if( "$lastB" == "") then
  echo "WARNING: could not get last B factors from refmacout.pdb"
  set lastB = 999
endif
echo "last B was $lastB, applying to unused waters"
grep HOH Bfac.pdb >! dummywater.pdb
combine_pdbs_runme.com B=$lastB dummywater.pdb dummywater.pdb outfile=Bwater.pdb > /dev/null
combine_pdbs_runme.com refmacout.pdb Bwater.pdb printref=1 Bfac.pdb outfile=new_Bfac.pdb > /dev/null
cp new_Bfac.pdb Bfac.pdb
cp Bfac.pdb Bfac_${itr}.pdb


skiprefmac:


set trigger = `echo $itr $phenix_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
  echo "running phenix.refine to refine B factors"
else
  goto skipphenix
endif

# now use phenix to refine B factors - better than previous?
set serial = `echo $itr | awk '{printf("%03d",$1)}'`
rm -f phenixBrefine_${serial}.pdb >& /dev/null
phenix.refine refme.mtz refme.pdb $ligcifs \
  prefix=phenixBrefine serial=$serial strategy=individual_adp >&! phenixBrefine_${itr}.log
if( $status ) then
  set BAD = "phenix.refine failed"
  goto exit
endif

cp Bfac.pdb prephenix_Bfac_${itr}.pdb
set lastB = `tail phenixBrefine_${serial}.pdb | awk '/^ATOM|^HETAT/{print substr($0,61,6)}' | sort -gr | head -n 1`
if( "$lastB" == "") then
  echo "WARNING: could not get last B factors from phenixBrefine_${serial}.pdb"
  set lastB = 999
endif
grep HOH Bfac.pdb >! dummywater.pdb
combine_pdbs_runme.com B=$lastB dummywater.pdb dummywater.pdb outfile=Bwater.pdb > /dev/null
combine_pdbs_runme.com phenixBrefine_${serial}.pdb Bwater.pdb printref=1 Bfac.pdb outfile=new_Bfac.pdb > /dev/null
cp new_Bfac.pdb Bfac.pdb
cp Bfac.pdb Bfac_${itr}.pdb


skipphenix:


echo "measuring challenges to previous restraints"
filter_pdb.awk -v skip=EP this.pdb |\
tee refme.pdb |\
filter_pdb.awk -v skip=H |\
awk '! /^ATOM|^HETAT/{print;next}\
  {print substr($0,1,60) "  0.00" substr($0,67)}' |\
cat >! zeroB.pdb

# find restraints bigger than the minimum value used
echo $cutoff_weight $pdbscale |\
cat - current_restraints.pdb |\
awk 'NR==1{minB=$1/$2;scale=$2;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {B=substr($0,61,6)+0}\
  B>minB{print}' |\
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








if( $void_scale == 0 ) then
  echo "skipping void calculation"
  set maxvoid = "n/d"
  set nbulk_sum = "n/d"
  set nbulk_max = "n/d"
  goto skipvoid
endif

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

skipvoid:
#dont forget to update this.pdb if the rst7 file gets shuffled


# check the pressure, if any
set press = ( 0 - - )
set outfiles = `ls -1t *.out`
foreach log ( barometer_${prev}.out ${laStage}.out amber_${prev}.out $outfiles )
  echo "checking $log"
  if(! -e "$log") continue
  set press = `awk '/PRESS/ && p{n=$3;print $NF} /A V E R A/{++p} END{print n+0}' $log | tail -n 3`
  # avg rms N
  if( "$press" == "" || "$press" == "0" ) set press = ( 0 - - )
  if( "$press[1]" != "0" ) break
end
if( "$press[1]" == "0" && -e ${laStage}.out ) then
   set ektot = `awk '/EKtot/ && ! p{++n} /EKtot/ && p{print $6} /A V E R A/{++p} END{print n+0}' ${laStage}.out`
endif

if( "$nwaters" == "" ) then
  set nwaters = `echo list | cpptraj -p xtal.prmtop | awk '$4=="atoms,"{print $11, $5, $3}' | head -n 1`
endif
echo "$prev pressure was $press   void: $maxvoid $nbulk_sum $nbulk_max   nwater: $nwaters" |\
  tee -a pressure_vs_itr.txt
# itr "pressure was"  <P> rms(P) nsteps "void:" maxvoid sum(nbulk) max(nbulk)  "nwater:" nwaters
#  1      2      3     4     5    6      7        8       9          10          11        12

set pressure_avglast_lasttime = 1
if( $?pressure_avglast_thistime ) then
  set pressure_avglast_lasttime = $pressure_avglast_thistime
endif
set pressure_avglast_thistime = "$pressure_avglast"
if( "$pressure_avglast" =~ *auto* ) then
  set pressure_avglast_thistime = `echo $pressure_avglast | awk '{print $1+0}'` 
  if( "$pressure_avglast_thistime" == "0" ) set pressure_avglast_thistime = 1
  set npressures = 0
  if(-e pressure_vs_itr.txt) then
    set npressures = `tail -n 30 pressure_vs_itr.txt | wc -l`
    echo -n "" >! pslope.txt
    foreach tail ( `seq 5 $npressures` )
      tail -n $tail pressure_vs_itr.txt | awk '{print $4,$12}' >! ${t}press.txt
      set pslope = `tac ${t}press.txt | awk '{print ++n,$1}' | linfit.awk `
      set medmad = `awk '{print $2}' ${t}press.txt | median.awk`
      echo $medmad | cat - ${t}press.txt |\
      awk 'NR==1{med=$1;mad=$3+med*0.05;next}\
         sqrt(($2-med)^2)>2*mad{print $2}' >! ${t}rejects.txt
      set rejects = `cat ${t}rejects.txt | wc -l`
      if( $rejects ) break
      echo "$tail $pslope $medmad" >> pslope.txt
    end
    echo $press |\
    cat - pslope.txt |\
    awk 'NR==1{p0=$1;sigma=$2/sqrt($3);next}\
       $1==0 && $2==0{next}\
       {print sqrt($2*$2)/sigma,$0}' |\
    cat >! pscore.txt
    set test = `sort -g pscore.txt | awk '$1<0.5{print $2;exit}'`
    if( "$test" != "" && ( $npressures > $test ) ) then
      set pressure_avglast_thistime = "$test"
    endif
  endif
  echo "pressure avglast this time: $pressure_avglast_thistime"
endif


if( "$pressure_avglast" != "1" && "$pressure_avglast_ramp" == "none" ) then
  set int_pressure_avglast = `echo $pressure_avglast_thistime | awk '{print int($1)}' | awk '$1<1{$1=1} {print}'`
  tail -n $int_pressure_avglast pressure_vs_itr.txt |\
    tee ${t}countme.txt |\
    awk '{++n;sum+=$4;vsum+=$5*$5} END{print sum/n,sqrt(vsum/n),$6}'
  set count = `cat ${t}countme.txt | wc -l`
  set press = `tail -n $int_pressure_avglast pressure_vs_itr.txt | awk '{++n;sum+=$4;vsum+=$5*$5} END{print sum/n,sqrt(vsum/n),$6}'`
  echo "using pressure: $press  averaged over $count ($pressure_avglast_thistime) itrs at itr $itr"
endif

# try to automate pressure scale
set pressure_scale_lasttime = 0
if( $?pressure_scale_thistime ) then
  set pressure_scale_lasttime = $pressure_scale_thistime
endif
set pressure_scale_thistime = "$pressure_scale"
if( "$pressure_scale" =~ *auto* ) then
  set pressure_scale_thistime = `echo $pressure_scale | awk '{print $1+0}'` 
  if( "$pressure_scale_thistime" == "0" ) set pressure_scale_thistime = $pressure_scale_lasttime
  if(-e pressure_vs_itr.txt) then
    tac pressure_vs_itr.txt |\
    awk '{print $4,$12}' |\
    awk '{p=$1;n=$2} n0+0<=0{n0=n} {dn=sqrt((n-n0)^2)}\
      dn/n0>0{++somesignal}\
      {print $0,somesignal+0,dn}\
      somesignal>10 && NR>30{exit}' >! ${t}press.txt
    # format: pressure Nwaters somesignal dn
    sort -k2gr ${t}press.txt |\
     awk 'NR==1{posP=($1>0);maxN=$2} \
          posP && $1<0 && ! zPN{zPN=$2} {print $0,maxN,zPN}\
          posP && $2<zPN-3*(maxN-zPN){exit}' |\
    cat >! ${t}notsmall.txt 
    set variety = `awk '$3+0>0{print $2}' ${t}press.txt | sort -u | wc -l`
    set linfit = `awk '{print $1,$2}' ${t}notsmall.txt | linfit.awk `
    set test = `echo $linfit | awk '{print sqrt($1*$1)}'`
    if( $#linfit == 2 && "$test" != "0" && $variety > 5 ) then
      set pressure_scale_thistime = "$test"
    endif
  endif
  echo "pressure scale this time: $pressure_scale_thistime"
endif

set changed_waters = 0
set dehydrate_thistime = 0
set dehydrate_rounded = 0
set nadd = 0
set trigger = `echo $itr $hydrate_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
  set w = `echo $press $pressure_avglast_thistime | awk '$2+0>0{snr=$1/$2*sqrt($NF);print 1.0/sqrt(1+1/snr**2)}'`
  echo "pressure snr weight = $w"
#  if( "$w" == "" || $pressure_avglast_thistime > 1 ) set w = 1
#  echo "pressure snr weight = $w"
  echo "checking pressure: $press[1] > $pressure_deadband "
  set test = `echo $press $pressure_deadband | awk '{print ( $1>$NF )}'`
  if( $test && "$dehydrate" != "0" ) then
    if( "$dehydrate" == "pressure" ) then
      set dehydrate_thistime = `echo $press $w $pressure_scale_thistime | awk '{print $1*$NF*$(NF-1)}'`
      echo "want to drop $dehydrate_thistime waters"
      set dehydrate_rounded = `echo $dehydrate_thistime | awk '{printf("%.0f",$1)}'`
#      if( $dehydrate_thistime == 0 ) set dehydrate_thistime = 1
    else
      set dehydrate_thistime = "$dehydrate"
      set dehydrate_rounded = "$dehydrate"
    endif
  endif

  echo "checking void size: $nbulk_max vs $minvoid"
  set bigvoid = `echo $nbulk_max $minvoid | awk '{print ( $1 > 2*$2 )}'`
  set lowpress = `echo $press -$pressure_deadband | awk '{print ( $1 < $2 )}'`
  set addv = `echo $nbulk_max $void_scale $minvoid | awk '{print $1*$2*( $1>$3 )}'`
  set addp = `echo $press $w $pressure_scale_thistime | awk '$NF==0{$NF=1} {print -$1*$NF*$(NF-1)*($1<0)}'`
  echo "nadd options: $addv void $addp press"
  set addv = `echo $addv | awk '{printf("%.0f",$1)}'`
  set addp = `echo $addp | awk '{printf("%.0f",$1)}'`
  if( $bigvoid || $lowpress ) then
    echo -n "may need to add water... "
    set nadd = `echo $addp $addv | awk '$1<$2{$1=$2} {print $1}'`
    if( "$nadd" == "0" ) then
      echo "but not enough."
    else
      echo "goal is $nadd"
    endif
  endif

  if( $water_lock ) then
    # keep water count constant
    set dwater = `echo $addp $dehydrate_thistime | awk '{print int(($1+$2)/2)}'`
    echo "average delta-water: $dwater"
    set nadd = $dwater
    set dehydrate_rounded = $dwater
    set revertStage = $laStage
  endif

  if( $water_dither ) then
    # keep water count moving, so that pressure scale can be tuned
    set min_dwater = `echo $press | awk '$1>0{print -1} $1<0{print 1} $1==0{print int(rand()-0.5)}'`
    set min_dwater = `echo $min_dwater $water_dither | awk '{print $1*$2}'`
    echo "delta-water signals: $addp $dehydrate_thistime min_dwater= $min_dwater"
    if( $min_dwater > 0 ) then
      set nadd = `echo $nadd $min_dwater | awk '$1<$2 && $2>0{$1=$2} {print $1}'`
    endif
    if( $min_dwater < 0 ) then
      set dehydrate_rounded = `echo $dehydrate_rounded $min_dwater | awk '{$2=-$2} $1<$2{$1=$2} {print $1}'`
    endif
    echo "new nadd = $nadd dehydrate rounded = $dehydrate_rounded"
  endif

  if( $dehydrate_rounded ) then
    echo "dropping ${dehydrate_rounded} waters"
    set changed_waters = 1

    echo "remapping waters based on Fo-Fc map..."
    rm -f remapped_restraints.pdb
    remap_waters_runme.com ${laStage}.rst7 current_restraints.pdb >! dehydrate_${itr}.log
    if( $status || ! -e remapped_restraints.pdb ) then
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

    dehydrate_amber_runme.com xtal.prmtop remapped.rst7 maxreject=${dehydrate_rounded} >> dehydrate_${itr}.log
    if( $status || ! -e drier.parm7 ) then
      set BAD = "dehydrate failed"
      goto exit
    endif
    egrep "stripping res" dehydrate_${itr}.log

    echo "updating xtal.prmtop"
    mv drier.parm7 xtal.prmtop
#    mv orignames.pdb prev_orignames.pdb
#    echo "updating orignames.pdb"
#    mv drier_orignames.pdb orignames.pdb
    echo "updating Bfac.pdb"
    cp -p Bfac.pdb predry_Bfac_${itr}.pdb
    cp -p drier_Bfac.pdb Bfac.pdb
    cp ${laStage}.in drier.in
    set lalaStage = $laStage
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
      refpoints=current_restraints.pdb kT=1.2 min_weight=$cutoff_weight \
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
     rm -f wetter.rst7
     hydrate_runme.com ${laStage}.rst7 outrst=wetter.rst7 \
        nadd=$nadd water_radius=$add_radius force_nadd=$water_lock debug=$debug \
        restraints=current_restraints.pdb >! hydrate_${itr}.log
      # interanlly updates xtal.prmtop, ref.crd 
      # never ever change orignames.pdb
      if( $status || ! -e wetter.rst7 ) then
        set BAD = "hydration failed"
        goto exit
      endif
      egrep "^added " hydrate_${itr}.log

      cp ${laStage}.in wetter.in
      set lalaStage = $laStage
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
        outfile=defused_restraints.pdb min_weight=$cutoff_weight \
        outcrd=ref.crd >! restraint_bomb_check.log
      if( $status ) then
        set BAD = "restraint bomb detected after hydration"
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
  echo "$prev $nwaters" | tee -a nwaters_vs_itr.txt

  # update the this and noEP pdb files
  rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log
  filter_pdb.awk -v skip=EP this.pdb >! refme.pdb
endif

set repick_now = 0

set trigger = `echo $itr $teleport_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $teleport_waters && -e cootme.mtz && $trigger ) then
  echo "looking for teleportable waters in ${laStage}.rst7"
else
  goto skipteleport
endif

touch actual_restraints.pdb
cat restraints_challenge_vs_weight.txt |\
awk '$2<2.4{print substr($0,index($0,$5)-6),"NOTME"}' |\
cat - actual_restraints.pdb |\
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
pick.com 4.5 fofc.map >! pick.log
awk '/height.sigma/{print;getline;print $0,"fofc.map"}' pick.log
echo "symgen $smallSG" | pdbset xyzin pick.pdb xyzout bigpick.pdb >> /dev/null
egrep "^ATOM|^HETAT" bigpick.pdb >> teleport_goals.pdb

# retrieve restraints that are actually being used
echo $cutoff_weight $pdbscale |\
cat - current_restraints.pdb |\
awk 'NR==1{minB=$1/$2;scale=$2;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {B=substr($0,61,6)+0}\
  B>minB{print}' |\
cat >! active_restraints.pdb

awk '/^ATOM|^HETAT/{print substr($0,12,17),substr($0,61,6)/100,"RESTRAINED"}' active_restraints.pdb >! restr.txt

# get beginning of trajectory
rst2pdb_runme.com include_ep=0 $trajectory trajstart.pdb Bfactors=none >! teleport_${itr}.log
#rst2pdb_runme.com ${laStage}.rst7 end.pdb >> /dev/null

filter_pdb.awk -v skip=H -v only=water refme.pdb >! rhome.pdb
rholabel_runme.com cootme.mtz rhome.pdb >> teleport_${itr}.log

filter_pdb.awk -v skip=H -v only=water trajstart.pdb rholabeled.pdb |\
 rmsd -v debug=1 |\
cat - restr.txt rholabeled.pdb |\
awk '{id=substr($0,1,17)} /moved/{moved[id]=substr($0,25,8);next}\
     $NF=="RESTRAINED"{weight[id]=$(NF-1);next}\
  ! /^ATOM|^HETAT/{next}\
   {id=substr($0,12,17);typ=substr($0,18,3);rho=$NF}\
   moved[id]=="" || weight[id]+0>0{next}\
    moved[id]+0<2.0 && typ=="HOH" && rho<0{\
  print rho/(moved[id]+1.0),substr($0,1,90),moved[id],$NF}' |\
sort -g |\
awk '{print substr($0,index($0,$2))}' >! trapped_waters.pdb


rm -f teleported.rst7
water_teleport_runme.com ${laStage}.rst7 maxmoves=$teleport_waters \
  teleport_goals.pdb notthese=active_restraints.pdb \
  mtzfile=cootme.mtz \
  mindist=$teleport_mindist \
  minrho=$teleport_minrho \
  teleportee=trapped_waters.pdb \
  debug=$debug minimize=0 energycheck=0 >> teleport_${itr}.log 

if(! -e teleported.rst7) then
  echo "no teleportation options"
  goto skipteleport
endif

awk '/Fo-Fc/{print;getline;print} /^teleporting /' teleport_${itr}.log 

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
cp new.pdb spiked_restraints.pdb
combine_pdbs_runme.com spiked_restraints.pdb orignames.pdb >> spikeweight.log
awk '{print substr($0,1,30) substr($0,61)}' new.pdb >! ${t}this.txt
awk '{print substr($0,1,30) substr($0,61)}' current_restraints.pdb >! ${t}that.txt
diff ${t}this.txt ${t}that.txt >> spikeweight.log
set changes = `diff ${t}this.txt ${t}that.txt | egrep "^<" | egrep "ATOM|HETAT" | wc -l`
mv new.pdb current_restraints.pdb
echo "spiked $changes weights in current_restraints.pdb"

echo "updating ref.crd"
update_centroid_positions_runme.com current_restraints.pdb outfile=ref.crd >> spikeweight.log

# check for bombs ?

cp ${laStage}.in teleported.in
if(-e ${laStage}.prmtop) cp ${laStage}.prmtop teleported.prmtop
set lalaStage = $laStage
set laStage = teleported
#set repick_now = 1

# update the this and noEP pdb files
rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log
filter_pdb.awk -v skip=EP this.pdb >! refme.pdb


skipteleport:

set trigger = `echo $itr $Breset_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
   echo "resetting all B factors to 20"
    convert_pdb.awk -v BFAC=20 Bfac.pdb >! new.pdb
    mv new.pdb Bfac.pdb
endif



set just_repicked = 0
set trigger = `echo $itr $repick_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger || $repick_now ) then
  set just_repicked = 1

  egrep "^CRYST" all_possible_refpoints.pdb | head -n 1 >! unfiltered_refpoints.pdb
  if(-e fofc.map) then
    pick.com 5 fofc.map >! pick.log
    echo "symgen $smallSG" | pdbset xyzin pick.pdb xyzout bigpick.pdb >> /dev/null
    convert_pdb.awk -v only=atoms -v BFAC=999 -v CHAIN=z -v output=pdb bigpick.pdb |\
    awk '{gsub("OW","O ")}' >> unfiltered_refpoints.pdb
  endif
  filter_pdb.awk -v only=atoms -v skip=water all_possible_refpoints.pdb    >> unfiltered_refpoints.pdb
  filter_pdb.awk -v only=atoms,water current_restraints.pdb                >> unfiltered_refpoints.pdb
  filter_pdb.awk -v only=atoms,water all_possible_refpoints.pdb            >> unfiltered_refpoints.pdb

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

  if( "$cutoff_forget" != "0" ) then
    echo $cutoff_weight $pdbscale |\
    cat sorted_refpoints.pdb |\
    awk 'NR==1{minB=$1/$2;next}\
      ! /^ATOM|^HETAT/{print;next}\
      {B=substr($0,61,6)+0} B>minB{print}' |\
    cat >! bigenough_refpoints.pdb
    cp sorted_refpoints.pdb full_refpoints_${itr}.pdb
    cp bigenough_refpoints.pdb sorted_refpoints.pdb
  endif

  # re-pick reference points
  echo "re-discovering nearest reference points"
  rm -f repicked_reference_points.pdb
  centroids_nearby_runme.com refme.pdb reffile=sorted_refpoints.pdb \
    softener=1 weight=$repick_maxweight maxdist=$repick_maxdist  \
    hohscale=$repick_hohscale \
    outfile=repicked_reference_points.pdb debug=$debug >! c2r_${itr}.log
  if( $status || ! -e repicked_reference_points.pdb ) then
    set BAD = "repick failed"
    goto exit
  endif
  awk 'NF<=2 || NF==3 && $2=="="{next} {print}' c2r_${itr}.log | awk '/unique selections/,/new_ref/'

  cp repicked_reference_points.pdb preprune_restraints.pdb

  echo "forgetting weights for big-move waters"
  # exclude recently teleported waters?
  egrep -h HOH current_restraints.pdb repicked_reference_points.pdb |\
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
     pruned_restraints.pdb repicked_reference_points.pdb \
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

set trigger = `echo $itr $filter_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
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

  diff current_restraints.pdb filtered_refpoints.pdb  |\
   awk '/^< / && /ATOM|HETAT/{print substr($0,3)}' >! filtered_out.pdb 

  set filtered = `cat filtered_out.pdb | wc -l`
  echo "filtered out $filtered solvent anchors as too close and/or too challenged"

  cp current_restraints.pdb unfiltered_restraints.pdb
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

set trigger = `echo $itr $remap_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if( $trigger ) then

  remap_waters_runme.com ${laStage}.rst7 current_restraints.pdb | tee remap_waters_${itr}.log
  if( $status ) then
    set BAD = "remap failed"
    goto exit
  endif

  cp remapped_restraints.pdb current_restraints.pdb

  echo "updating ref.crd"
  update_centroid_positions_runme.com current_restraints.pdb outfile=ref.crd >> remap_waters_${itr}.log

  # check for bombs ?

  rst2pdb_runme.com ${laStage}.rst7 this.pdb >! rst2pdb.log
  filter_pdb.awk -v skip=EP this.pdb >! refme.pdb

  cp ${laStage}.in remapped.in
  set lalaStage = $laStage
  set laStage = remapped

endif

if( "$thrubond_avg_weight" != "0" ) then
  echo "averaging weights through bonds with avgfac=$thrubond_avg_weight and spread=$thrubond_avg_weight_spread"
  cp current_restraints.pdb presmooth_restraints.pdb
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
  rmsd prepower_restraints.pdb restraints_powered.pdb | grep Bfac
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
set trigger = `echo $itr $delete_badrest_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if (  $trigger ) then
  echo "deleting highly challenged refpoints "

  cp current_restraints.pdb predel_restraints.pdb

  delete_worst_restraints.com outprefix="" >! delete_worst_${itr}.log 

endif


# try releasing alt-conf atoms
set trigger = `echo $itr $delete_altconf_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( ( ! $just_repicked ) && $trigger ) then
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

if(! -e fofc_Rplot.txt) then
  echo "WARNING: skipping best-itr determination. No fofc_Rplot.txt file yet"
  set bestR = "unk"
  set R = "unk"
  goto skiprandel
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


# try releasing random atoms, but only after R factor has stabilized
if ( "$randel_itr" != "0" && ( ( ( $itr - $randel_last ) > $randel_itr ) || $randel_trigger_now ) ) then
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
skiprandel:




# try releasing tightest restraints, see if they grow back
set trigger = `echo $itr $release_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger || "$release_trigger" != "0" && $pegged_weight ) then
  echo "releasing over-challenged restraints "
  set tootight = `echo $max_weight | awk '{print ($1+0.0001)*0.9}'`
  rm -f new_restraints.pdb
  release_worst_restraints_runme.com restraintpdb=current_restraints.pdb \
    refinedpdb=refme.pdb radius=$release_radius pdbscale=$pdbscale \
    tootight=$tootight maxbad=$release_maxbad \
    resetweight=$release_weight sigma=5 \
    ignore="" >! release_worst_${itr}.log
  if( $status || ! -e new_restraints.pdb) then
    set BAD = "release worst restraints failed"
    goto exit
  endif
  grep "too tight" release_worst_${itr}.log
  tail -n 3 release_worst_${itr}.log 

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
set trigger = `echo $itr $rhocheck_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
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
set trigger = `echo $itr $reinin_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( $trigger ) then
  # taking all defaults
  echo -n "checking CAs : "
  cp current_restraints.pdb pre_clip_CA.pdb
  reinin_errant_CA_refpoints.com outfile=clip_CA.pdb | tee reinin.log | awk '$2!="="'

  rmsd clip_CA.pdb current_restraints.pdb | egrep "MAXD.Bfac|WARNING"
  cp clip_CA.pdb current_restraints.pdb
endif


set trigger = `echo $itr $minwt_itr | awk -F "[ +]" '$2+0==0{print $2+0;exit} {print ($1 % $2+0 == $3+0)}'`
if ( ! $trigger ) goto skip_minwt

# make sure any ligands stay put
set liggrep = `echo $ligands | awk '{gsub(" ","|");print}'`
set test = `egrep "$liggrep" all_possible_refpoints.pdb | wc -l`
if( $test && "$min_lig_weight" != "0" && "$ligands" != "" ) then

  echo "ensuring slight restraints on ligands: $ligands "
  egrep "$liggrep" all_possible_refpoints.pdb |\
  filter_pdb.awk -v only=ligand -v ligands="$liggrep" |\
  awk '{print $0,"MAYBELIG"}' >! lig.pdb
  cat current_restraints.pdb lig.pdb |\
  awk '{id=substr($0,12,5)" "substr($0,18,12)} ! seen[id]{print;++seen[id]}' |\
  awk -v mlw=$min_lig_weight -v scale=$pdbscale '! /^ATOM|^HETAT/{print;next}\
    $NF!="MAYBELIG"{print;next}\
    {B=substr($0,61,6)+0;pre=substr($0,1,60);post=substr($0,67,12);\
     w=B*scale}\
    w<mlw{w=mlw}\
    {printf("%s%6.2f%s    NOT_A_BOMB\n",pre,w/scale,post)}' |\
  cat >! wrongorder.pdb
  combine_pdbs_runme.com wrongorder.pdb refme.pdb outfile=clip_ligand.pdb > /dev/null
  rmsd clip_ligand.pdb current_restraints.pdb | egrep "MAXD.Bfac"
  cp clip_ligand.pdb current_restraints.pdb

endif


# make sure certain CAs stay put
if( "$min_CA_weight" != 0 && -e alignment_atnums.txt ) then

  set minB = `echo $min_CA_weight $pdbscale | awk '{print $1/$2}'`
  echo "ensuring slight restraints on high-density CA atoms used for alignment "

  # quick update in case scrambled
  if( "$align_target" == "restraints" ) cp current_restraints.pdb align_ref.pdb

  # first, put weight on atoms used in alignment
  cat alignment_atnums.txt orignames.pdb |\
  awk '$1+0>0{++sel[$1];next}\
    ! /^ATOM|^HETAT/{next} {++n}\
    ! sel[n]{next}\
    {print substr($0,12,18),"SEL";}' |\
  cat - align_ref.pdb |\
  awk -v minB=$minB '$NF=="SEL"{++sel[substr($0,1,18)];next}\
    ! /^ATOM|^HETAT/{next}\
    {id=substr($0,12,5)" "substr($0,18,12);\
     pre=substr($0,1,60);post=substr($0,67,12);}\
    sel[id]{printf("%s%6.2f%s\n",pre,minB,post)}' |\
  awk 'substr($0,12,5)=="  CA "' >! ${t}align_xyz.pdb
  # make sure they are on the restratins list
  cat current_restraints.pdb ${t}align_xyz.pdb |\
  awk '{id=substr($0,12,5)" "substr($0,18,12)} ! seen[id]{print;++seen[id]}' |\
  cat>! ${t}wrongorder.pdb
  # new master restraint list with old weights and filled in with any missing alignment atoms
  combine_pdbs_runme.com ${t}wrongorder.pdb refme.pdb outfile=${t}complete.pdb > /dev/null
  # pull out all restraint weights for alignment atoms
  combine_pdbs_runme.com ${t}align_xyz.pdb saveBfac=1 ${t}complete.pdb outfile=${t}align_w.pdb  > /dev/null
  # make sure they are all at or above the threshold
  cat ${t}align_w.pdb |\
  awk -v minB=$minB '! /^ATOM|^HETAT/{next}\
    {B=substr($0,61,6)+0;pre=substr($0,1,60);post=substr($0,67,12)}\
    B<minB{B=minB}\
    {printf("%s%6.2f%s\n",pre,B,post)}' |\
  cat >! ${t}align_neww.pdb
  # now inject updated weights into the master list
  combine_pdbs_runme.com ${t}align_neww.pdb ${t}complete.pdb printref=1 outfile=new_restraints.pdb > /dev/null
  rmsd new_restraints.pdb current_restraints.pdb | egrep "MAXD.Bfac"
  cp current_restraints.pdb prealignrest_restraints.pdb
  cp new_restraints.pdb current_restraints.pdb

endif

# make sure atoms used for alignment are restrained
if( "$min_align_weight" != 0 && ! -e alignment_atnums.txt ) then
  echo "WARNING: missing alignment_atnums.txt"
endif
if( "$min_align_weight" != 0 && -e alignment_atnums.txt ) then
  # quick update in case scrambled
  if( "$align_target" == "restraints" ) cp current_restraints.pdb align_ref.pdb

  set minB = `echo $min_align_weight $pdbscale | awk '{print $1/$2}'`
  echo "ensuring slight restraints on atoms used for alignment "
  # first, must extract atoms used in alignment
  cat alignment_atnums.txt orignames.pdb |\
  awk '$1+0>0{++sel[$1];next}\
    ! /^ATOM|^HETAT/{next} {++n}\
    ! sel[n]{next}\
    {print substr($0,12,18),"SEL";}' |\
  cat - align_ref.pdb |\
  awk -v minB=$minB '$NF=="SEL"{++sel[substr($0,1,18)];next}\
    ! /^ATOM|^HETAT/{next}\
    {id=substr($0,12,5)" "substr($0,18,12);\
     pre=substr($0,1,60);post=substr($0,67,12);}\
    sel[id]{printf("%s%6.2f%s\n",pre,minB,post)}' >! ${t}align_xyz.pdb
  # make sure they are on the restratins list
  cat current_restraints.pdb ${t}align_xyz.pdb |\
  awk '{id=substr($0,12,5)" "substr($0,18,12)} ! seen[id]{print;++seen[id]}' |\
  cat>! ${t}wrongorder.pdb
  # new master restraint list with old weights and filled in with any missing alignment atoms
  combine_pdbs_runme.com ${t}wrongorder.pdb refme.pdb outfile=${t}complete.pdb > /dev/null
  # pull out all restraint weights for alignment atoms
  combine_pdbs_runme.com ${t}align_xyz.pdb saveBfac=1 ${t}complete.pdb outfile=${t}align_w.pdb  > /dev/null
  # make sure they are all at or above the threshold
  cat ${t}align_w.pdb |\
  awk -v minB=$minB '! /^ATOM|^HETAT/{next}\
    {B=substr($0,61,6)+0;pre=substr($0,1,60);post=substr($0,67,12)}\
    B<minB{B=minB}\
    {printf("%s%6.2f%s\n",pre,B,post)}' |\
  cat >! ${t}align_newweight.pdb
  # now inject updated weights into the master list
  combine_pdbs_runme.com ${t}align_newweight.pdb ${t}complete.pdb printref=1 outfile=align_restraints.pdb > /dev/null
  rmsd align_restraints.pdb current_restraints.pdb | egrep "MAXD.Bfac"
  cp current_restraints.pdb prealignrest_restraints.pdb
  cp align_restraints.pdb current_restraints.pdb

  # check for bombs?
  

endif

skip_minwt:

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
  outfile=defused_restraints.pdb min_weight=$cutoff_weight \
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

touch weight_hist_vs_itr.txt
awk -v pdbscale=$pdbscale '/^ATOM|^HETAT/{print substr($0,61,6)*pdbscale}' current_restraints.pdb |\
 histogram.awk -v bs=1 |\
awk -v itr=$itr '{print itr,$0}' |\
cat >> weight_hist_vs_itr.txt
echo "" >> weight_hist_vs_itr.txt


# create ref.crd from most recent rst7 and current_restraints.pdb
rst2pdb_runme.com ${laStage}.rst7 flying.pdb Bfactors=none >! final_centroids.log
#cp flying.pdb start_${itr}.pdb

update_centroid_positions_runme.com \
     pdbfile=current_restraints.pdb \
     orignames=flying.pdb topfile=xtal.prmtop \
     minwt=$cutoff_weight pdbscale=$pdbscale \
     outfile=ref.crd >> final_centroids.log

update_centroid_positions_runme.com \
      pdbfile=align_ref.pdb \
      outfile=align_ref.crd  >> final_centroids.log


echo "replacing restraint list in ${laStage}.in with new list in ${Stage}.in"
restraintlist2amber.com list=current_restraints.pdb \
  pdbscale=$pdbscale \
  minwt=$cutoff_weight \
  orignames=flying.pdb \
  allatom_weight=$allatom_weight \
  ${laStage}.in \
  outfile=${Stage}.in |\
 tee -a restraint_update_${itr}.log | grep bins

#  extract reference positions as a pdb file
rst2pdb_runme.com ref.crd ref.pdb >> restraint_update_${itr}.log

# convert all restrained atoms back to pdb for a bomb check
echo "extracting full restraint list"
cat ${Stage}.in |\
awk '/Specific/{atom=$2;getline;w=$1}\
  /^RES/{for(n=$2;n<=$3;++n){print w,atom,n}}' |\
sort -gr |\
tee actual_weights.txt |\
cat - ref.pdb |\
awk -v s=$pdbscale 'NF==3{w[$2,$3]=$1;next}\
  /^CRYST/{print} ! /^ATOM|^HETAT/{next}\
    {atom=substr($0,12,5);gsub(" ","",atom);\
     pre=substr($0,1,60);post=substr($0,67)}\
   w[atom,$NF]{printf("%s%6.2f%s\n",pre,w[atom,$NF]/s,post)}' |\
cat >! actual_restraints.pdb

cp actual_restraints.pdb actual_restraints_${itr}.pdb

# bomb check
echo "checking for full-list bombs"
restraint_bomb_detector.com ${laStage}.rst7 \
  refpoints=actual_restraints.pdb  >! restraint_bomb_check_actual.log
if( $status ) then
  cat restraint_bomb_check.log
  set BAD = "too many restraint bombs detected"
  goto exit
else
  echo "none detected"
endif


# do not bring allatom_weight restraints into the current_restraints
# now set in stone
#set test = `echo $max_mult | awk '{print $1+0}'`
#if( 0 && "$test" != "1" ) cp actual_restraints.pdb current_restraints.pdb
cp current_restraints.pdb restraints_for_${itr}.pdb

echo -n "max min density: "
awk '{split(FILENAME,w,"_");f=w[3]+0} /Maximum dens/{print f,m,$NF} /Minimum dens/{m=$NF}' restraint_update_${itr}.log |\
 tee -a dens_vs_itr.txt

egrep  "^ATOM|^HETAT" actual_restraints.pdb |\
awk -v pdbscale=$pdbscale '{printf("%10s %s\n", substr($0,61,6)*pdbscale,substr($0,12,17))}' |\
sort -gr >! sorted_weights.txt
set avg = `awk '{sum+=$1;++n} END{if(n)print sum/n}' sorted_weights.txt`
set medmad = `awk '{print $1}' sorted_weights.txt | median.awk`
set max = `awk '{print $1;exit}' sorted_weights.txt`
set worst = `awk '! t{$1="";print;++t} ! / ZN | SE /{print;exit}' sorted_weights.txt`
set mults = `awk -F ":" '/^highest|^lowest/{getline;print $2+0}' restraint_update_${itr}.log`
#if("$max" == "") set max = "-"
echo -n "worstweight: "
echo "$itr $max $avg $mults  $medmad  $worst" | tee -a worstweight_vs_itr.txt

set pegged_weight = `echo $max $max_weight | awk '{print ( NR>=2 && $1 >= $2 )}'`
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




# update chiral and omega restraint lists
#goto skipchir
wait
set omegalogs = `ls -1rt omegalyze_*.log | tail -n 3`
if( "$omegalogs" == "" ) then
  echo "WARNING: no informaiton on omega"
  goto skipchir
endif
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
    if(b==1)print "  r1=134., r2=135., r3=225., r4=226., rk2 =50, rk3=50,  &end"}}' |\
cat >! bad_omega.rst


set chiralogs = `ls -1rt chiralyze_*.log | tail -n 3`
if( "$chiralogs" == "" ) then
  echo "WARNING: no informaiton on chirality"
  goto skipchir
endif
cat $chiralogs |\
awk '$2=="ideal" || $2=="delta:" || NF==0{next}\
  /All restrained atoms within/{next}\
  {rn=substr($0,4,6);atom=substr($0,15,5);\
   gsub(" ","",atom);\
   list=list" "atom;}\
     /*sigma$/{print "BAD",rn,list;++bads;list=""}' |\
sort -u >! bad_chir.txt
cat bad_chir.txt orignames.pdb |\
awk '\
   {gsub("CA N C CB","N C HA CB");\
    gsub("CB CA CG1 CG2","CA CG1 CG2 HB");\
    gsub("CB CA OG1 CG2","CA CG2 OG1 HB")}\
   /^BAD/{++baddies;\
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
    n=split(list,w,",");if(n!=5)continue;\
    print "# chirality for",res[b];\
    end="&end";if(b==1)end="";\
    print " &rst iat=" list,end;\
    if(b==1)print "  r1=50., r2=60.,  r3=80.,  r4=90., rk2 =10, rk3=10,  &end"}}' |\
cat >! bad_chir.rst

set badchir = `cat bad_chir.txt | wc -l`
set badomega = `cat bad_omega.txt | wc -l`

if( $badchir ) echo "$badchir questionable chiral centers"
if( $badomega ) echo "$badomega questionable omega twists"
if( ( $badchir || $badomega ) && ( $omega_weight != 0 || $chiral_weight != 0 ) ) then
   echo "maybe replace chir_omega.rst with bad chiral and omega lists."
   echo "cat bad_chir.rst bad_omega.rst >! chir_omega.rst"
endif

skipchir:


# adjust temp for unrestrained atoms
if ( "$thermostat" != "" ) then
  set prevtemp = `awk -F "=" '/temp0/{print $2+0}' ${laStage}.in | tail -n 1`
  echo "previous temperature setting was: $prevtemp"
  if("$prevtemp" == "0") set prevtemp = "$thermostat"

  cat ${laStage}.in |\
  awk '/^RES/{for(i=$2;i<=$3;++i)print i}' |\
   sort -u | sort -g |\
  awk 'NR==1{s=e=$1;next} $1==e+1{e=$1;next} {print s"-"e;s=e=$1} END{print s"-"e}' |\
  awk -F "-" '$1==$2{print $1;next} {print}' |\
  sort -u | sort -g >! rest_ranges.txt 
  set rest_ranges = `cat rest_ranges.txt `
  set rest_ranges = `echo $rest_ranges | awk '{gsub(" ",",");print}'`

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
echo $temperature $gamma_ln |\
cat - ${Stage}.in |\
awk 'NR==1{temp=$1;gamma_ln=$2;next} \
  {space=substr($0,1,index($0,$1)-1)}\
  / temp0/{split($1,w,"=");print space w[1] "=" temp ",";next}\
  # set thermostat \
  / gamma_ln=/{split($1,w,"=");print space w[1] "=" gamma_ln ",";next}\
  {print}' |\
cat >! tempfile.in
mv tempfile.in ${Stage}.in

if( ! $barostat ) then
 cat ${Stage}.in |\
 awk '{space=substr($0,1,index($0,$1)-1);split($1,w,"=");esc=" "$NF;if(esc!=" /")esc=""}\
  / ntb=/{print space "ntb=1,ntp=0," esc;next}\
  {print}' |\
 cat >! tempfile.in
 mv tempfile.in ${Stage}.in
endif

set nsnb = `head -n 100 ${Stage}.in | awk -F "=" '/ nsnb=/{print $2+0;exit}'`
if( "$nsnb" == "" ) set nsnb = 0
echo "equi_ns = $equi_ns"
echo $equi_ns $equi_dt $align_nstlim $nsnb $equi_gamma |\
cat - ${Stage}.in |\
awk 'NR==1{ns=$1;dt=$2;align_nstlim=$3;nsnb=$4+0;equi_gamma=$5+0;\
    if(dt<1e-6)dt=1e-6;\
    nstlim=int(ns*1000/dt);\
    if(align_nstlim>nstlim)align_nstlim=nstlim; \
    if(align_nstlim)nstlim=align_nstlim;next} \
  {space=substr($0,1,index($0,$1)-1);split($1,w,"=");esc=" "$NF;if(esc!=" /")esc=""}\
  # shorter time step \
  / dt=/{print space w[1] "=" dt "," esc;next}\
  # turn off shake \
#  / ntf=2/{print space "ntf=1 ," esc;next}\
#  / ntc=2/{print space "ntc=1 ," esc;next}\
  # nonbond updates every cycle \
#  ! nsnb && / cut=/{print space "nsnb=1," esc}\
  # just a few cycles \
  / nstlim=/{print space w[1] "=" nstlim "," esc;next}\
  / ntpr=| ntwr=/{print space w[1] "=" nstlim "," esc;next}\
  / ntwx=/{print space w[1] "=" 0 "," esc;next}\
  # heavy thermostat \
  / gamma_ln=/ && equi_gamma{print space w[1] "=" equi_gamma "," esc;next}\
  {print}' |\
cat >! equi.in

# additional settling step at 1/10 speed
echo $equi_ns $equi_dt $align_nstlim $nsnb $equi_gamma $settle_slowdown |\
cat - equi.in |\
awk 'NR==1{ns=$1;dt=$2;align_nstlim=$3;nsnb=$4+0;gamma=$5+0;slowdown=$6;\
    dt=dt/slowdown;\
    ns=ns/slowdown/slowdown;\
    gamma=gamma*slowdown;\
    if(dt<1e-6)dt=1e-6;\
    nstlim=int(ns*1000/dt);\
    if(align_nstlim>nstlim)align_nstlim=nstlim; \
    if(align_nstlim)nstlim=align_nstlim;next} \
  {space=substr($0,1,index($0,$1)-1);split($1,w,"=");esc=" "$NF;if(esc!=" /")esc=""}\
  # shorter time step \
  / dt=/{print space w[1] "=" dt "," esc;next}\
  # turn off shake \
#  / ntf=2/{print space "ntf=1 ," esc;next}\
#  / ntc=2/{print space "ntc=1 ," esc;next}\
  # nonbond updates every cycle \
  ! nsnb && / cut=/{print space "nsnb=1," esc}\
  # just a few cycles \
  / nstlim=/{print space w[1] "=" nstlim "," esc;next}\
  / ntpr=| ntwr=/{print space w[1] "=" nstlim "," esc;next}\
  / ntwx=/{print space w[1] "=" 0 "," esc;next}\
  # heavy thermostat \
  / gamma_ln=/ && equi_gamma{print space w[1] "=" equi_gamma "," esc;next}\
  {print}' |\
cat >! settle.in

echo $barometer_cycles 0.0005 |\
cat - ${Stage}.in |\
awk 'NR==1{n=$1;dt=$2;next}\
  {space=substr($0,1,index($0,$1)-1);split($1,w,"=");esc=" "$NF;if(esc!=" /")esc=""}\
# just a few cycles \
/ nstlim=/{print space w[1] "=" n "," esc;next}\
/ dt=/{print space w[1] "=" dt "," esc;next}\
/ ntpr=/{print space w[1] "=" n/10 "," esc;next}\
/ ntwx=| ntwr=| ntr=/{print space w[1] "=" 0 "," esc;next}\
/ ntb=/{print space "ntb=2,ntp=4," esc;next}\
/     Specific /{exit}\
/ restraint|wt type=|DISANG|LISTIN=POUT|dummy |nmropt=1/{\
   if(esc != "") print esc;\
   next}\
{print}' |\
cat >! barometer.in

echo "prod_ns = $prod_ns"
echo $prod_ns $dt $write_ps |\
cat - ${Stage}.in |\
awk 'NR==1{prod=$1;dt=$2+0;write=$3+0;if(dt<1e-6)dt=0.002;\
     nstlim=int(prod*1000/dt);\
     ntwx=int(write/dt);next} \
  {space=substr($0,1,index($0,$1)-1);split($1,w,"=");}\
  / nstlim=/{print space w[1] "=" nstlim ",";next}\
  / ntpr=/{print space w[1] "=" ntwx ",";next}\
  / ntwx=/{;print space w[1] "=" ntwx ",";next}\
  / ntwr=/{print space w[1] "=" ntwx ",";next}\
  / dt=/{print space w[1] "=" dt ",";next}\
  {print}' |\
cat >! tempfile.in
mv tempfile.in ${Stage}.in
set ntwx = `awk -F "=" '/ ntwx=/{print $2+0;exit}' ${Stage}.in`

if(-e chir_omega.rst) then
  echo "adjusting chiral/omega weights in chir_omega.rst : $chiral_weight $omega_weight"
  echo $chiral_weight $omega_weight |\
  cat - chir_omega.rst |\
  awk 'NR==1{cw=$1;ow=$2;next}\
    /chirality/{type="chiral"}\
    /trans-omega/{type="omega"}\
    ! ( /r1=/ && $NF=="&end" ) {print;next}\
    type=="chiral"{print "   r1=50., r2=60.,  r3=80.,  r4=90., rk2 ="cw", rk3="cw",  &end"}\
    type=="omega" {print "   r1=134., r2=135., r3=225., r4=226., rk2 ="ow", rk3="ow",  &end"}' |\
  cat >! new.txt
  if( $debug ) diff chir_omega.rst new.txt
  mv new.txt chir_omega.rst
endif


set disang = `grep DISANG ${Stage}.in | wc -l`
if( ! $disang && ( $omega_weight != 0 || $chiral_weight != 0 ) && -e chir_omega.rst) then
  echo "applying DISANG restraints"

  cat ${Stage}.in |\
  awk '$1~/^&/{section=substr($1,2)}\
    NF==1 && $1=="/" && section=="cntrl"{\
    print "  nmropt=1, /";\
    print "&wt type=\047END\047 /";\
    print "DISANG=chir_omega.rst";\
    print " &dummy  i=1, ";}\
   {print}' |\
  cat >! tempfile.in
  if( $debug ) diff ${Stage}.in tempfile.in
  mv tempfile.in ${Stage}.in
endif
if( $disang && "$omega_weight" == "0" && "$chiral_weight" == "0" ) then
  echo "removing chiral/omega restraints"

  egrep -v "dummy|DISANG|LISTIN=PO|wt type=|nmropt" ${Stage}.in >! tempfile.in
  if( $debug ) diff ${Stage}.in tempfile.in
  mv tempfile.in ${Stage}.in
endif


set netfrc_set = `grep netfrc ${Stage}.in | wc -l`
if( ! $netfrc_set && "$netfrc" != "1" ) then
  echo "applying netfrc=$netfrc to stabilize run with weak restarints"

  cat ${Stage}.in |\
  awk -v netfrc=$netfrc '$1~/^&/{section=substr($1,2)}\
    NF==1 && $1=="/" && section=="cntrl"{\
    print prev,"/";\
    print " &ewald";\
    print "  netfrc=" netfrc ",";}\
   {prev=$0;print}' |\
  cat >! tempfile.in
  if( $debug ) diff ${Stage}.in tempfile.in
  mv tempfile.in ${Stage}.in
endif



# make input file for re-alignmnent sub-run
echo $prod_ns $dt $align_nstlim |\
cat - ${Stage}.in |\
awk 'NR==1{nstlim=int($1/$2*1000);align_nstlim=$3;\
    if(align_nstlim>nstlim)align_nstlim=nstlim;\
    if(align_nstlim)nstlim=align_nstlim;next}\
  /=/{space=substr($0,1,index($0,$1)-1);split($1,w,"=");\
   value=w[2]+0;if(value>nstlim)value=nstlim}\
  / nstlim=/{print space w[1] "=" nstlim ",";next}\
  / ntpr=| ntwr=/{print space w[1] "=" value ",";next}\
  {print}' |\
cat >! ${Stage}_i.in
sync ${Stage}_i.in


foreach substage ( settle equi )

  if( "$equi_ns" == "0" ) continue

  touch drift_byframe_${substage}.txt 
  touch ${Stage}_${substage}.out 
  echo -n "" >! ${t}trajin.txt
  set nstlim0 = `echo $equi_ns $equi_dt | awk '{print int($1/$2*1000)}'`
  set subruns = `echo $nstlim0 $align_nstlim | awk '$2>0{print int($1/$2)}'`
  if( "$subruns" == "" || "$subruns" == "0" ) set subruns = 1

  foreach i ( `seq 1 $subruns` )
    echo "${substage}-ing $Stage ( $i / $subruns )"
    set ss = ${Stage}_${substage}_${i}
    $pmemd -O -i ${substage}.in -o ${ss}.out \
     -p xtal.prmtop \
     -c ${laStage}.rst7 \
     -ref ref.crd \
     -r ${ss}.rst7 \
     -x ${ss}.nc \
     -inf ${Stage}_${substage}.mdinfo
    if($status) then
      set BAD = "amber pre-run failed with status $status at $ss "
      goto exit
    endif

    cat ${ss}.out >> ${Stage}_${substage}.out
    if(-e ${ss}.nc) then
      echo "trajin ${ss}.nc" >> ${t}trajin.txt
    endif

    grep NaN ${ss}.out
    if(! $status) break
    egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' ${ss}.out | grep '\*\*\*\*\*'
    if(! $status) break

    echo -n "re-aligning "
    cpptraj -p xtal.prmtop -y ${ss}.rst7 -c align_ref.crd << EOF >! align.log
    rmsd rmsd reference norotate @$align_mask out rmsd.txt savevectors combined vecsout vecsout.txt
    trajout realigned.rst7
EOF
    cp realigned.rst7 ${ss}.rst7
    set lalaStage = $laStage
    set laStage = ${ss}

    set drift = `tail -n 1 vecsout.txt | awk '{print $2,$3,$4,"(",sqrt($2*$2+$3*$3+$4*$4),")"}'`
    echo -n "drift: "
    echo -n "$itr $i " >> drift_byframe_${substage}.txt 
    echo "$drift"  | tee -a drift_byframe_${substage}.txt 

    # give a warning for big drift?
    set test = `echo $drift | awk '{print ( $(NF-1) > 0.1 ) }'`
    if( $test ) then
      echo "WARNING: drift > 0.1 A detected! consider lowering align_nstlim or increasing min_align_weight"
  endif

  end

  cp ${ss}.rst7 ${Stage}_${substage}.rst7

  set test = `cat ${t}trajin.txt | wc -l`
  if( $test ) then
    echo "consolidating nc files into ${Stage}_${substage}.nc"
    echo "trajout ${Stage}_${substage}.nc" >> ${t}trajin.txt
    cat ${t}trajin.txt | cpptraj -p xtal.prmtop >! ${t}concat.log
  endif

end

if(-e ${Stage}_${substage}.out) then
  grep NaN ${Stage}_${substage}.out
  if(! $status) break
  egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' ${Stage}_${substage}.out | grep '\*\*\*\*\*' 
  if(! $status) break
endif

cp ${laStage}.rst7 start_${itr}.rst7
rm -f ${Stage}.rst7
touch ${Stage}.out 
set nstlim0 = `echo $prod_ns $dt | awk '{print int($1/$2*1000)}'`
set subruns = `echo $nstlim0 $align_nstlim | awk '$2>0{print int($1/$2)}'`
if( "$subruns" == "" || "$subruns" == "0" ) set subruns = 1

# see if we are aggretating ncfiles or rst7 files
set nstlimi = `head -n 1000 ${Stage}_i.in | awk -F "=" '/ nstlim=/{print $2+0}'`
set ntwxi = `head -n 1000 ${Stage}_i.in | awk -F "=" '/ ntwx=/{print $2+0}'`
set ntwri = `head -n 1000 ${Stage}_i.in | awk -F "=" '/ ntwr=/{print $2+0}'`
set ext = rst7 ; if( $ntwxi < $nstlimi ) set ext = nc
echo -n "" >! ${t}trajin.txt

foreach i ( `seq 1 $subruns` )
    echo "running $Stage ( $i / $subruns )"
    $pmemd -O -i ${Stage}_i.in -o ${Stage}_${i}.out \
     -p xtal.prmtop \
     -c ${laStage}.rst7 \
     -ref ref.crd \
     -r ${Stage}_${i}.rst7 \
     -x ${Stage}_${i}.nc \
     -inf ${Stage}.mdinfo
    if($status) then
      set BAD = "amber run failed with status $status at ${Stage} $i"
      goto exit
    endif

    grep NaN ${Stage}_${i}.out
    if(! $status) break
    egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' ${Stage}_${i}.out | grep '\*\*\*\*\*'
    if(! $status) break

    cpptraj -p xtal.prmtop -y ${Stage}_${i}.rst7 -c align_ref.crd << EOF >! align.log
    rmsd rmsd reference norotate @$align_mask out rmsd.txt savevectors combined vecsout vecsout.txt
    trajout realigned.rst7
EOF
    set drift = `tail -n 1 vecsout.txt | awk '{print $2,$3,$4,"(",sqrt($2*$2+$3*$3+$4*$4),")"}'`
    echo -n "drift: "
    echo "$itr $i $drift" | tee -a drift_byframe.txt 
    cp realigned.rst7 ${Stage}_${i}.rst7

     # give a warning for big drift?
    set test = `echo $drift | awk '{print ( $(NF-1) > 0.1 ) }'`
    if( $test ) then
      echo "WARNING: drift > 0.1 A detected! consider lowering align_nstlim or increasing min_align_weight"
    endif


    set laStage = ${Stage}_${i}
   
    # aggregate output logs
    cat ${Stage}_${i}.out >> ${Stage}.out
    rm ${Stage}_${i}.out
    # build list of files to merge at the end
    echo "trajin ${Stage}_${i}.${ext}" >> ${t}trajin.txt
  end
  # check for errors
  grep NaN ${Stage}.out
  if(! $status) break
  egrep -v 'Mask|NSTEP2=|^\*\*\*\*\*\*' ${Stage}.out | grep '\*\*\*\*\*'
  if(! $status) break

  cp ${Stage}_${i}.rst7 ${Stage}.rst7

#  ls -1rt | egrep "^${Stage}_" | awk -F "_" '$3+0>0 && /.rst7$|.nc$/' >! ${t}ls.txt
#  egrep "${ext}"'$' ${t}ls.txt |\
#    awk '{print "trajin",$1}' >! ${t}trajin.txt
  set test = `cat ${t}trajin.txt | wc -l`
  if( $test > 0 ) then
    echo "consolidating ${ext} files into ${Stage}.nc"
    echo "trajout ${Stage}.nc" >> ${t}trajin.txt
    cat ${t}trajin.txt | cpptraj -p xtal.prmtop >! ${t}concat.log
  endif

  awk '/A V E R A/{++p;n=split(FILENAME,w,"_");printf("%d ",w[n]);next}\
      p && /NSTEP/{printf("%s %s %s ",$6,$9,$12);next}\
      p && NF>2{for(i=1;i<=NF;++i)if($i!~/[A-Za-z=]/ && $i!="1-4")printf("%s ",$i)}\
      p && ( /EAMBER/ || /-----------------/ ){print "";p=0}' ${Stage}.out |\
  tail -n 1 |\
  tee -a amber_energy_vs_itr.txt
# amber_energy_vs_itr.txt columns (from the amber .out "A V E R A G E S" block):
#   1 itr        2 TIME(ps)    3 TEMP(K)   4 PRESS     5 Etot      6 EKtot
#   7 EPtot      8 BOND        9 ANGLE    10 DIHED    11 1-4 NB   12 1-4 EEL
#  13 VDWAALS   14 EELEC      15 EHBOND   16 RESTRAINT 17 EAMBER (non-restraint)
# col16 (RESTRAINT) is the total restraint energy; a flat col16 over the last N
# iterations means the weight optimization has stabilized (used by the opt2
# monitor in the example runme).  EHBOND (col15) is always 0 in these runs.

  set files = `awk '/^trajin/{print $2}' ${t}trajin.txt | awk -F "." '{print $1".rst7";print $1".nc"}'`
  if( $debug ) echo "deleting $files"
  if( $#files ) rm $files
#  rm `awk '/^trajin /{print $NF}' ${t}trajin.txt`
#  rm ${Stage}_*[0-9].out
#  rm ${Stage}_*[0-9].nc

  # reduce number of frames?
  set skip = `echo $nstlimi $ntwxi | awk '{print int($2/$1)}'`
  if( $skip > 1 ) then
    echo "skipping every $skip frames for final nc file"
    cpptraj -p xtal.prmtop << EOF >> ${t}concat.log
    trajin ${Stage}.nc 1 last $skip
    trajout smaller.nc
EOF
    mv smaller.nc ${Stage}.nc
  endif

  if(! -e ${Stage}.nc) then
    set BAD = "amber didnt work: no ${Stage}.nc file"
    goto exit
  endif

  # now quick check of geometry
  cpptraj -p xtal.prmtop -y ${Stage}.rst7 << EOF >! checks.log
checkchirality chir out chircheck.txt
strip :WAT,HOH@Y1,EPW
check reportfile geocheck.txt
multidihedral omega omega out omegaline.txt range360
EOF
  set ninv = `awk '$2!=1 && $1+0>0' chircheck.txt | wc -l`
  set worstchir = `awk '$2!=1 && $1+0>0{print $1+0;exit}' chircheck.txt `
  if("$worstchir" == "") set worstchir = "n/a"
  cat omegaline.txt |\
   awk 'NR==1{for(i=2;i<=NF;++i){split($i,w,":");oresnum[i]=w[2]};next}\
     {for(i=2;i<=NF;++i){dev=sqrt(($i-180)^2);print oresnum[i],dev,$i}}' >! omega.txt
  set worstomega = `sort -k2gr omega.txt | head -n 1`
  set ncis = `awk '$2>90{print}' omega.txt | wc -l`
  set ngeoproblems = `cat geocheck.txt | wc -l`
  set worstgeo = `awk '{print;exit}' geocheck.txt`
  echo "$ninv inverted chiral centers, $ncis cis peptides ( $worstomega[2] deg) and $ngeoproblems geometry problems"
  touch quickgeo_vs_itr.txt
  echo "$itr   $ninv $ncis $ngeoproblems    $worstchir $worstomega  $worstgeo" >> quickgeo_vs_itr.txt

  # optionally bail out (with error) the moment any cis peptide shows up
  if( $exit_on_cis && $ncis ) then
    set BAD = "$ncis cis peptides at itr $itr ( worst omega residue $worstomega[1], see quickgeo_vs_itr.txt )"
    goto exit
  endif

  if( $ninv == 0 && "$chiral_weight" != "0" ) then
    echo "turning off chiral restraints"
    set chiral_weight_ramp = none
  endif
  set rising = `echo $chiral_weight_ramp | awk -F "[: _-]" '$NF~/x/{print ( $NF+0 > 1 ) }' | awk '{print int($1+0)}'`
  if( "$rising" == "" ) set rising = 0
  if( $ninv && "$chiral_weight" == "0"  && "$chiral_weight_ramp" == "none" ) then
    set chiral_weight = `echo $chiral_weight 0.1 | awk '{print $1+$2}'`
    set chiral_weight_ramp = ${chiral_weight}-50:1.1x
    echo "turning on chiral restraint ramp: $chiral_weight_ramp"
    set ninv = 0
  endif
  set rising = `echo $omega_weight_ramp | awk -F "[: _-]" '$NF~/x/{print ( $NF+0 > 1 ) }' | awk '{print int($1+0)}'`
  if( "$rising" == "" ) set rising = 0
  if( $ncis == 0 && "$omega_weight" != "0" && $rising ) then
    echo "turning down omega restraints"
#    set omega_weight = `echo $omega_weight 0.5 | awk '{print $1*$2}'`
    set omega_weight_ramp = "${omega_weight}-0:0.5x"
  endif
  if( $ncis && ! $rising ) then
    set omega_weight = `echo $omega_weight 0.1 | awk '{print $1+$2}'`
    set omega_weight_ramp = ${omega_weight}-1000:2x
    echo "turning on omega restraint ramp: $omega_weight_ramp"
    set ncis = 0
  endif
  if( $ninv || $ncis || $ngeoproblems ) then
    echo "WARNING: should be exiting because of geometry problems."
#    break
  endif


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

