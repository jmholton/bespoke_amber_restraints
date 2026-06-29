#! /bin/tcsh -f
# 6c2r_EG7_2x_example_runme.com
# Aurora Kinase A with EG7 inhibitor — 2×2×2 supercell crystallographic Amber MD
# Starting from PDB deposit 6c2r (EG7 inhibitor, 1.9 Å, P41212).
# EG7 force-field uses Dirk's CIF (elbow.EG7_04202024.cif) with correct protonation
# state, which differs from the PDB monomer library version.
#
# See also: example_setup_notes.com (1aho, no ligand, 0.965 Å, single cell)
#
# Requirements: CCP4 Suite, Phenix Suite, Amber 22+, gnuplot
#
# Usage: copy this file to a fresh working directory and run section by section.
# This is a reference script, not a single pipeline — paste and execute each
# block interactively, reviewing outputs between sections.
#
# After leap2amber.com completes (Section 10), the opt1/opt2/... workflow is
# structurally identical for all systems; only a few parameter values differ:
#
#   System         min_lig_weight   halfrho_neg   notes
#   6c2r / EG7     0.1              3.5           this script
#   6c2r / AMPPNP  0.01             3.5           Dirk's private dataset; badlinks=1
#   1aho           0                3.5           no ligand


#==============================================================================
# Environment setup
#==============================================================================
set pdbid     = 6c2r
set super_mult = 2,2,2
set salt      = ( NH4 SO4 )
set salt_conc = 0.15
set badlinks  = 0      # EG7 does not create spurious Lys-LIG bonds
                       # Set badlinks=1 for the AMPPNP / Dirk dataset

set pdir = ~/projects/git/bespoke_amber_restraints
set skit = ~/projects/amber/6c2r_AMPPNP/claude/starter_kit   # starter kit location
set t    = tempfile
set path = ( $pdir $path )

if (! -d $skit) then
  echo "ERROR: starter kit not found at $skit"
  echo "       Edit the 'set skit = ...' line near the top of this script"
  goto exit
endif

if (-e compute_settings.sourceme) source compute_settings.sourceme

if (! $?AMBERHOME) source /programs/amber22/amber.csh
if (! $?PHENIX)    source /programs/phenix/phenix_env.csh
if (! $?CCP4)      source /programs/ccp4-9/bin/ccp4.setup-csh

if (-e user_settings.sourceme) source user_settings.sourceme

set debug = 1


#==============================================================================
# SECTION 1 — Download and sequence
#==============================================================================
# Download structure factors and coordinates from PDB (replaces deprecated getcif.com)
phenix.fetch_pdb $pdbid action=all
# Produces: 6c2r.pdb  6c2r-sf.cif

phenix.cif_as_mtz ${pdbid}-sf.cif
# Produces a .mtz file — name varies by phenix version; grab the newest one:
set rawmtz = `ls -1t ${pdbid}*.mtz | head -1`
if ($rawmtz == "") then
  echo "ERROR: phenix.cif_as_mtz produced no .mtz file"
  goto exit
endif
echo "Raw MTZ: $rawmtz"

# Extract protein sequence from SEQRES
grep SEQRES ${pdbid}.pdb |\
sequence.awk |\
awk 'NR==1{$0="> "$0} {print} NF==0{exit}' |\
tee seq.fasta

set modulo = `awk 'NR>1{printf("%s",$1)}' seq.fasta | wc -c`
echo "modulo (residues per chain) = $modulo"

# Use PDB deposit as starting ASU coordinates
cp ${pdbid}.pdb starthere_asu.pdb


#==============================================================================
# SECTION 2 — Reflection data and crystal properties
#==============================================================================
# Auto-detect structure factor amplitude column from the raw MTZ
set F    = `mtzdmp $rawmtz | awk 'NF>10 && ! /2FOFCWT|FWT|DELF/ && $(NF-1)=="F"{print $NF;exit}'`
set SIGF = `mtzdmp $rawmtz | awk 'NF>10 && ! /2FOFCWT|FWT|DELF/ && $(NF-1)=="Q"{print $NF;exit}'`
set FREE = `mtzdmp $rawmtz | awk '/FreeR|Free_R|free_R|R-free/{print $NF;exit}'`
echo "Using columns: F=$F  SIGF=$SIGF  FREE=$FREE"
if ($F == "" || $SIGF == "" || $FREE == "") then
  echo "ERROR: could not auto-detect MTZ columns — check mtzdmp output above"
  goto exit
endif

cad hklin1 $rawmtz hklout refme_small.mtz << EOF
labin file 1 E1=$F E2=$SIGF E3=$FREE
labou file 1 E1=FP E2=SIGFP E3=FreeR_flag
EOF
if ($status) then
  echo "ERROR: cad failed — check column names above"
  goto exit
endif

echo head | mtzdump hklin refme_small.mtz >! smallmtzdump.txt
set SGnum   = `awk '/Space group =/{print $NF+0}' smallmtzdump.txt | tail -n 1`
set smallSG = `awk -v num=$SGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
set smallCELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' smallmtzdump.txt`
set nsymops = `awk -v SG=$smallSG '$4==SG {print $2}' ${CLIBD}/symop.lib | head -1`
set reso    = `awk '/Resolution Range :/{getline;getline;print $(NF-2)+0;exit}' smallmtzdump.txt | head -1`

# Write crystal properties file — this will be linked into each subdirectory
# and auto-sourced by generate_alignment_reference_runme.com to pass badlinks
cat << EOF >! xtal_properties.sourceme
set reso       = $reso
set smallSG    = $smallSG
set smallCELL  = ( $smallCELL )
set nsymops    = $nsymops
set super_mult = $super_mult
set modulo     = $modulo
set ligands    = auto
set salt       = ( $salt )
set salt_conc  = $salt_conc
set badlinks   = $badlinks
EOF
source xtal_properties.sourceme

# Expand MTZ to 2×2×2 supercell
cad hklin1 refme_small.mtz hklout expanded.mtz << EOF | tee mtz_expand.log
labin file 1 all
outlime space 1
EOF
cad hklin1 expanded.mtz hklout refme_cell.mtz << EOF | tee -a mtz_expand.log
labin file 1 all
symm 1
EOF
rm -f expanded.mtz

set reidx = `echo $super_mult | awk -F "[ ,x]" '{print "reindex h"$1",k"$2",l"$3}'`
echo $reidx |\
reindex hklin refme_cell.mtz hklout refme.mtz | tee -a mtz_expand.log
if ($status) then
  echo "ERROR: MTZ expansion to supercell failed"
  goto exit
endif
rm -f refme_cell.mtz

echo head | mtzdump hklin refme.mtz >! mtzdump.txt
set CELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' mtzdump.txt`

echo "" >! blank.pdb
echo "CELL $CELL\nSPACE 1" | pdbset xyzin blank.pdb xyzout ${t}cell.pdb
egrep "^CRYST1" ${t}cell.pdb >! cell.pdb
rm -f blank.pdb ${t}cell.pdb


#==============================================================================
# SECTION 3 — Ligand force-field files
#==============================================================================
# Auto-detect ligands from PDB, excluding salt ions and water
set ligands = `filter_pdb.awk -v only=ligand,atoms starthere_asu.pdb |\
  awk '/^HETAT/{print substr($0,18,3)}' | sort -u |\
  awk -v salt="NH4 SO4 HOH WAT" \
    'BEGIN{n=split(salt,s);for(i=1;i<=n;i++)bad[s[i]]=1} {if($1 in bad)next; print}'`
echo "Ligands found: $ligands"

mkdir ligands
cd ligands

# EG7: use Dirk's CIF, which has the correct protonation state.
# The PDB monomer library (--chemical_component EG7) has a different protonation
# that does not match the experimental model.
cp ${skit}/37C_EG7/ligands/elbow.EG7_04202024.cif EG7.cif
if ($status) then
  echo "ERROR: EG7.cif not found in starter kit at $skit"
  goto exit
endif

foreach lig ( $ligands )
  if (-e ${lig}.cif) then
    phenix.elbow ${lig}.cif --id=${lig} --opt --opt_nproc=10 --amber_force_field_files
  endif
  if (-e ${lig}.mol2) continue
  phenix.elbow --chemical_component $lig --id=${lig} --amber_force_field_files \
    --opt --opt_nproc=10
end

# Salt ions — always from PDB library
foreach lig ( $salt )
  if (-e ${lig}.mol2) continue
  phenix.elbow --chemical_component $lig --id=${lig} --amber_force_field_files \
    --opt --opt_nproc=10
end

cd ..


#==============================================================================
# SECTION 4 — HIS protonation (QM-determined)
#==============================================================================
# HIS protonation is determined by QM using the Phenix quantum interface.
# Pre-computed results for this protein are in the starter kit.
#
# To regenerate from scratch (only if the protein sequence changes):
#   Working directory: ~/projects/his_flips/ (or any clean directory)
#   cp starthere_asu.pdb Dirk.pdb
#   cp ligands/EG7.cif .
#   phenix.ready_set Dirk.pdb           # adds H; produces Dirk.updated.pdb
#   echo 1 | mmtbx.quantum_interface iterate_NQH=HIS Dirk.updated.pdb \
#     | tee options_Dirk.log            # identifies HIS residues; generates .phil files
#   # For each generated Dirk.updated_A_NNN_HIS.phil:
#   mmtbx.quantum_interface Dirk.updated.pdb iterate_NQH=HIS \
#     Dirk.updated_A_NNN_HIS.phil run_qmr=True qi.nproc=6 EG7.cif |& tee philN.log
#   # Collect votes and tabulate (see ~/projects/his_flips/HIS_protonation_QM_notes.com)
#   # Output: HIS_settings.txt — copy to this directory as HIS_settings_asu.txt
#
# Format: one line per HIS — HIE (Nε), HID (Nδ), HIP (both / positively charged)
# 6c2r result: HIP176 HIP187 HIE190 HID201 HIP248 HIP254 HID280 HIE306 HID366 HIE380
cp ${skit}/HIS_settings_asu.txt HIS_settings_asu.txt
if ($status) then
  echo "ERROR: HIS_settings_asu.txt not found in starter kit at $skit"
  echo "       Regenerate using ~/projects/his_flips/HIS_protonation_QM_notes.com"
  goto exit
endif


#==============================================================================
# SECTION 5 — Build and sanitize starting model (build1)
#==============================================================================
mkdir build1
cd build1

ln -sf ../starthere_asu.pdb starthere.pdb
ln -sf ../refme_small.mtz refme.mtz
ln -sf ../xtal_properties.sourceme .
cp ../ligands/*.cif .
set ligcifs = `ls -1 *.cif |& awk '/.cif/'`

buildout_pdb_runme.com starthere.pdb $ligcifs badlinks=$badlinks >&! buildout.log
if (! -e minRfree.pdb) then
  echo "ERROR: buildout_pdb_runme.com did not produce minRfree.pdb — check buildout.log"
  goto exit
endif

# Review output and pick best model by Rfree:
grep "Final R" *.log | justify.awk | sort -k7g | head -5

ln -sf minRfree.pdb thisone.pdb

cd ..


#==============================================================================
# SECTION 6 — ASU tleap charge check (amber_asu)
#==============================================================================
mkdir amber_asu
cd amber_asu

cp ../ligands/*.mol2 .
cp ../ligands/*.frcmod .
ln -sf ../build1/thisone.pdb starthere.pdb
cp ../HIS_settings_asu.txt HIS_settings.txt

# Strip altlocs and duplicate atoms
awk '! /^ATOM|^HETAT/{print;next}\
     {print substr($0,1,16),substr($0,18)}' starthere.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
cat >! amberme.pdb

# Dynamic tleap stub — one entry per mol2 found
cat << EOF >! tleap_stub.in
source leaprc.protein.ff19SB
source leaprc.water.opc
set default FlexibleWater on
source leaprc.gaff2
loadAmberParams frcmod.ff19SBmodAA
loadOff mod_amino19.lib
EOF
foreach mol2 ( *.mol2 )
  set lig = `basename $mol2 .mol2`
  cat << EOF >> tleap_stub.in
$lig = loadMol2 $mol2
loadAmberParams ${lig}.frcmod
EOF
end
cat << EOF >> tleap_stub.in
x = loadpdb tleapme.pdb
set x box { $CELL[1] $CELL[2] $CELL[3] }
set default nocenter on
saveAmberParm x xtal.prmtop start.crd
quit
EOF

cat HIS_settings.txt amberme.pdb |\
convert_pdb.awk -v output=amber -v fixEe=1 |\
awk '/HIS|HIE|HID|HIP/ && $NF=="H"{next} {print}' |\
egrep -v "LINK" >! tleapme.pdb

tleap -f tleap_stub.in | tee tleap.log

grep "unperturbed charge" tleap.log
set charge0 = `awk '/unperturbed charge/{gsub(/[)(]/,"");print int($7);exit}' tleap.log`
set charge  = `echo $charge0 $nsymops | awk '{print $1*$2}'`
echo "ASU charge $charge0, cell charge $charge" | tee cell_charge.txt

cd ..


#==============================================================================
# SECTION 7 — Generate centroid reference points (garr1)
#==============================================================================
mkdir garr1
cd garr1

ln -sf ../build1/thisone.pdb starthere.pdb
ln -sf ../refme_small.mtz refme.mtz
cp ../ligands/*.cif .
ln -sf ../xtal_properties.sourceme .   # REQUIRED: passes badlinks to no_new_nonbonds_runme.com
set ligcifs = `ls -1 *.cif |& awk '/.cif/'`

# Note on badlinks (AMPPNP only): badlinks=1 suppresses a spurious Lys141-NZ to
# phosphate covalent bond that phenix creates by proximity. EG7 does NOT need it.
# If badlinks=1 is omitted for AMPPNP on the very first garr run, phenix_opts_unbump.eff
# will contain a NZ-LIG bond reference that crashes any subsequent restart when phenix
# can no longer find that bond. Fix: delete phenix_opts_unbump.eff AND start.geo, then
# restart from itr=1. The xtal_properties.sourceme link above prevents this.

generate_alignment_reference_runme.com starthere.pdb $ligcifs \
  repulse_nb=100 crush_nb=100 repulse_scale=0.5 >&! garr.log
if (! -e centroids_final_001.pdb) then
  echo "ERROR: garr did not produce centroids_final_001.pdb — check garr.log"
  goto exit
endif
# Re-run (continuing from last itr) if Rfree has not converged.

cd ..


#==============================================================================
# SECTION 8 — Build supercell centroid set (centroids/)
#==============================================================================
mkdir centroids
cd centroids

# Pick best map: lowest Rfree across build1 and garr1
grep "Final R" ../build1/*.log ../garr1/*.log | justify.awk | sort -k7g | tee sorted.txt | head
set minRfree = `awk '/_/{gsub(".log:"," ");print $1;exit}' sorted.txt`
ln -sf ${minRfree}.mtz minRfree.mtz

ln -sf ../garr1/centroids_final_001.pdb centroids_asu.pdb
ln -sf ../amber_asu/amberme.pdb fulllength_asu.pdb
ln -sf ../refme_small.mtz .
cp ../ligands/*.cif .
set ligcifs = `ls -1 *.cif |& awk '/.cif/'`

# Verify centroid orientation vs full-length model
flip_to_target_runme.com centroids_asu.pdb fulllength_asu.pdb > /dev/null
filter_pdb.awk -v skip=H,water fulllength_asu.pdb flipped.pdb | rmsd | head

# 2Fo-Fc reference density map
set Fwt  = `mtzdmp minRfree.mtz | awk 'NF>5 && /2FOFCWT|FWT/ && ! / DELF/ && $(NF-1)=="F"{print $NF;exit}'`
set PHIwt = `mtzdmp minRfree.mtz | awk 'NF>5 && /PH2FOFCWT|PHWT/ && ! /DEL/ && $(NF-1)=="P"{print $NF;exit}'`
cad hklin1 minRfree.mtz hklout reference0.mtz << EOF
labin file 1 E1=$Fwt E2=$PHIwt
labou file 1 E1=Fref E2=PHIref
EOF
fft hklin reference0.mtz mapout ffted.map << EOF
labin F1=Fref PHI=PHIref
EOF
mapmask mapin ffted.map mapout reference0.map << EOF
xyzlim asu
scale sigma
EOF

# Rename atoms for supercell expansion: protein stays put, ligands/waters get ordinal names
foreach prefix ( fulllength centroids )
  egrep "^CRYST|^SSBON" ${prefix}_asu.pdb >! renamed.pdb
  filter_pdb.awk -v only=protein,atoms ${prefix}_asu.pdb >> renamed.pdb
  convert_pdb.awk -v only=ligand,atoms -v renumber=ordinal,watS -v fixEe=1 \
    ${prefix}_asu.pdb >> renamed.pdb
  convert_pdb.awk -v only=water,atoms -v renumber=ordinal,watS \
    ${prefix}_asu.pdb >> renamed.pdb
  cp renamed.pdb ${prefix}_renamed.pdb
  awk '! /^ATOM|^HETAT/{print;next} {print substr($0,1,16),substr($0,18)}' \
    ${prefix}_renamed.pdb |\
  awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
  cat >! ${prefix}_asu_noalt.pdb
end
combine_pdbs_runme.com centroids_asu_noalt.pdb fulllength_asu_noalt.pdb > /dev/null
filter_pdb.awk -v skip=water new.pdb centroids_asu_noalt.pdb | rmsd
# All RMS should be zero with no warnings

# Expand to supercell — run twice: first to discover monomer map, then apply it
expand2supercell_runme.com fulllength_renamed.pdb refme_small.mtz super_mult=$super_mult \
  outprefix=fulllength_super phenix_bumpcheck=0 debug=1 | tee fulllength_expand0.log
cp fulllength_super.pdb fulllength_super0.pdb
cp new_monomer_rot_trans.txt monomer_rot_trans.txt

expand2supercell_runme.com fulllength_renamed.pdb refme_small.mtz super_mult=$super_mult \
  outprefix=fulllength_super phenix_bumpcheck=0 debug=1 \
  mono_map=monomer_rot_trans.txt | tee fulllength_expand.log

expand2supercell_runme.com centroids_renamed.pdb refme_small.mtz super_mult=$super_mult \
  refpdb=fulllength_renamed.pdb outprefix=centroids_super \
  phenix_bumpcheck=0 debug=1 \
  mono_map=monomer_rot_trans.txt | tee centroids_expand.log

filter_pdb.awk -v skip=H,water centroids_super.pdb fulllength_super.pdb | rmsd | head
filter_pdb.awk -v skip=H,water centroids_super.pdb fulllength_super.pdb | grep CA | rmsd | head

foreach prefix ( fulllength centroids )
  awk '! /^ATOM|^HETAT/{print;next} {print substr($0,1,16),substr($0,18)}' \
    ${prefix}_super.pdb |\
  awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
  cat >! ${prefix}_noalt.pdb
end
combine_pdbs_runme.com centroids_noalt.pdb fulllength_noalt.pdb
filter_pdb.awk -v skip=water new.pdb centroids_noalt.pdb | rmsd
# All RMS should be zero with no warnings

# Label centroids by their 2Fo-Fc density value; only atoms in density are kept
filter_pdb.awk -v skip=H centroids_super.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  {occ=substr($0,55,6)+0}\
  occ>0.01{print substr($0,1,80),"           |",$NF}' >! all_possible_centroids.pdb

rholabel_runme.com all_possible_centroids.pdb reference0.mtz mtzlabel=Fref

cat rholabeled.pdb |\
awk '/^CRYST/{print;next} ! /^ATOM|^HETAT/{next}\
  {rho=$NF;pre=substr($0,1,60);post=substr($0,67)}\
  rho>0{printf("%s%6.2f%s\n",pre,rho,post)}' |\
tee centroids_in_density.pdb | grep EG7
if (! -e centroids_in_density.pdb) then
  echo "ERROR: rholabel produced no centroids_in_density.pdb"
  goto exit
endif

filter_pdb.awk -v skip=water,H centroids_in_density.pdb fulllength_super.pdb | rmsd | head
filter_pdb.awk -v skip=H,water centroids_in_density.pdb fulllength_super.pdb | grep CA | rmsd | head

cd ..


#==============================================================================
# SECTION 9 — Supercell phenix refinement (super_refine1)
#==============================================================================
mkdir super_refine1
cd super_refine1

ln -sf ../refme.mtz .
ln -sf ../centroids/fulllength_super.pdb starthere.pdb
cp ../ligands/*.cif .
set ligcifs = `ls -1 *.cif |& awk '/.cif/'`

cat << EOF >! opts.eff
refinement {
  refine {
    occupancies {
      individual = water
    }
  }
  bulk_solvent_and_scale {
    apply_back_trace = False
  }
  main {
    max_number_of_iterations = 100
  }
  pdb_interpretation {
    automatic_linking {
      link_none = True
    }
  }
}
EOF

awk '{print substr($0,1,80)}' starthere.pdb >! refme.pdb
phenix.refine ../refme.mtz refme.pdb prefix=phenix opts.eff $ligcifs >&! phenix1.log
if (! -e phenix_001.pdb) then
  echo "ERROR: phenix.refine failed — check phenix1.log"
  goto exit
endif

# Select single conformer via B-weighted jiggling
awk '{print substr($0,1,80)}' phenix_001.pdb >! jiggleme.pdb
jigglepdb.awk -v seed=1 -v shift=byB -v shift_scale=0.01 \
  -v disulfide_links=1 -v independent_confsel=1 jiggleme.pdb |\
awk '/^CRYST/{print} ! /^ATOM|^HETAT/{next}\
  {occ=substr($0,55,6)+0;}\
  occ>0{print}' |\
awk '{pre=substr($0,1,16);mid=substr($0,18,38);B=substr($0,61);\
   print pre,mid " 1.00" B}' |\
cat >! confsel.pdb

awk '! /^ATOM|^HETAT/{print;next}\
  {print substr($0,1,16),substr($0,18)}' phenix_001.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
cat >! fulllength_noalt.pdb

reorganize_pdb_runme.com confsel.pdb refpdb=fulllength_noalt.pdb outfile=oneconf.pdb \
  phenix_bumpcheck=0 debug=1 autorerun=0 | tee reorg_oneconf.log

reorganize_waters.com oneconf.pdb super_mult=$super_mult \
  smallmtz=../centroids/reference0.mtz | tee rewater.log

awk '{print substr($0,1,80)}' rewatered.pdb >! reorgme.pdb
reorganize_pdb_runme.com reorgme.pdb refpdb=fulllength_noalt.pdb outfile=reorged.pdb \
  phenix_bumpcheck=0 debug=1 autorerun=0 | tee reorg_rewatered.log

awk '{print substr($0,1,80)}' reorged.pdb >! refme_noalt.pdb
phenix.geometry_minimization refme_noalt.pdb prefix=confsel_min $ligcifs \
  automatic_linking.link_none=True nonbonded_weight=500 >&! confsel_min.log
if (! -e confsel_min_001.pdb) then
  echo "ERROR: phenix.geometry_minimization failed — check confsel_min.log"
  goto exit
endif

# Fix cis-peptides using generate_omega_fix_runme.com (produces omega_fix.eff)
generate_omega_fix_runme.com confsel_min_001.pdb >&! omega_fix_gen.log
phenix.refine confsel_min_001.pdb ../refme.mtz prefix=omegafix1 opts.eff $ligcifs \
  omega_fix.eff >&! phenix_omegafix1.log

ln -sf omegafix1_001.pdb thisone.pdb

cd ..


#==============================================================================
# SECTION 10 — Initial Amber MD run (amber1)
#==============================================================================
# IMPORTANT: skip Cpu and Min stages — they reliably produce NaN coordinates
# ("black holes") on 6c2r. Start directly from Cool.
mkdir amber1
if (-e compute_settings.sourceme) ln -sf ../compute_settings.sourceme amber1/
cd amber1

cp ../ligands/*.mol2 .
cp ../ligands/*.frcmod .
cp ../ligands/*.pdb .
set ligs = `ls -1 *.mol2 | awk -F "." '{print $1}'`
ln -sf ../refme.mtz .
ln -sf ../super_refine1/thisone.pdb starthere.pdb
cp starthere.pdb refined.pdb

# Expand HIS protonation states from ASU to supercell via modulo arithmetic
cp ../HIS_settings_asu.txt HIS_settings_asu.txt
awk '{print $1,substr($2,2)}' HIS_settings_asu.txt >! unshifted_HIS_settings.txt
grep "CA  HIS" refined.pdb >! HIS.pdb
head -n 1 HIS.pdb | awk '{print substr($0,23,4)}' |\
cat - unshifted_HIS_settings.txt |\
awk 'NR==1{firstHIS=$1;next} NR==2{unshifted=$2}\
  {print $1,$2+firstHIS-unshifted}' |\
cat >! shifted_HIS_settings.txt

echo $modulo |\
cat - shifted_HIS_settings.txt HIS.pdb |\
awk 'NR==1{modulo=$1;next}\
  /^HI[PEDS]/ && NF==2{proton[$2]=$1;next}\
  ! /^ATOM|^HETAT/{next}\
  {chain=substr($0,22,1);resnum=substr($0,23,6)+0;\
   modresnum=(resnum)%modulo;\
   print proton[modresnum],chain resnum;}' |\
cat >! HIS_settings.txt

cp HIS_settings.txt protonation.txt

set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' refined.pdb`

# Dynamic tleap stub — one entry per mol2; water model handled by leap2amber.com
cat << EOF >! tleap_stub.in
source leaprc.protein.ff19SB
source leaprc.water.opc
set default FlexibleWater on
source leaprc.gaff2
loadAmberParams frcmod.ff19SBmodAA
loadOff mod_amino19.lib
EOF
foreach lig ( $ligs )
  cat << EOF >> tleap_stub.in
$lig = loadMol2 ${lig}.mol2
loadAmberParams ${lig}.frcmod
EOF
end
cat << EOF >> tleap_stub.in
x = loadpdb tleapme.pdb
set x box { $CELL[1] $CELL[2] $CELL[3] }
set default nocenter on
saveAmberParm x xtal.prmtop start.crd
quit
EOF

# Initial density restraints from centroids in density
set itr = 0
set pdbscale = 0.01
set B0 = `echo 1 $pdbscale | awk '{print $1/$2}'`
cat ../centroids/centroids_in_density.pdb |\
awk -v B0=$B0 '/^CRYST|^LINK|^SSBO/{print} ! /^ATOM|^HETAT/{next}\
    {pre=substr($0,1,60);post=substr($0,67);rho=$NF; B=B0*rho}\
    B>B0{B=B0}\
    {printf("%s%6.2f%s\n",pre,B,post)}' |\
cat >! all_possible_refpoints.pdb

centroids_nearby_runme.com refined.pdb reffile=all_possible_refpoints.pdb \
  softener=2 weight=Bfac maxdist=1 hohscale=1 \
  outfile=restraints_for_${itr}.pdb debug=$debug | tee c2r_${itr}.log
cp restraints_for_${itr}.pdb current_restraints.pdb

# Identify salt cation and anion from mol2 charges
foreach ion ( $salt )
  set icharge = `awk '/TRIPOS.ATOM/,/TRIPOS.BOND/{sum+=$NF} END{printf("%.1g",sum)}' ${ion}.mol2`
  echo "$ion charge $icharge"
  set ispos = `echo $icharge | awk '{print ($1 > 0)}'`
  if ($ispos) then
    set cation = "$ion"
    ln -sf ${ion}.pdb cation.pdb
  else
    set anion = "$ion"
    ln -sf ${ion}.pdb anion.pdb
  endif
end

# Supercell charge: unit cell charge (from amber_asu) × number of cells
set cellcharge = `awk '/cell charge/{print $NF+0}' ../amber_asu/cell_charge.txt`
set ncells     = `echo $super_mult | awk -F "[ ,x]" '{print $1*$2*$3}'`
set charge     = `echo $cellcharge $ncells | awk '{print $1*$2}'`
echo "Supercell charge = $charge  (unit cell $cellcharge × $ncells cells)"

add_salt_runme.com refined.pdb conc=$salt_conc \
  RIP=4 RIW=3 charge=$charge anion=$anion cation=$cation | tee add_salt.log

reorganize_pdb_runme.com salty.pdb ignore_zero=0 refpdb=refined.pdb \
  outfile=reorganized.pdb phenix_bumpcheck=0 declash=1 | tee reorganize_final.log

filter_pdb.awk -v skip=water reorganized.pdb | egrep -v "^END" >! amberme.pdb
filter_pdb.awk -v skip=H -v only=water,atoms reorganized.pdb >> amberme.pdb

# Estimate padding waters: fill void volume to 150% of water capacity
egrep -v HOH amberme.pdb >! dry.pdb
echo | rwcontents xyzin dry.pdb >! ${t}rwcontents.log
grep "% of cell without atoms" ${t}rwcontents.log
set waterslots = `awk '/Cell volume/{V=$NF;n=55*6.022e23/1e27*V} \
  /% of cell without atoms/{print int(n*$NF/100),int(n)}' ${t}rwcontents.log`
echo "Room for $waterslots[1] waters ($waterslots[2] total in cell)"
set gotwater = `grep "O   HOH" amberme.pdb | wc -l`
set padwater = `echo $waterslots $gotwater | awk '{print 0+sprintf("%.2g",($1-$3)*1.5)}'`
echo "padwater = $padwater"

cp restraints_for_${itr}.pdb initial_restraints.pdb

# Run MD stages: skip Cpu and Min (cause NaN / black holes on this system)
leap2amber.com amberme.pdb stages=Cool,Heat,Equi,EquiMin,Prod \
  protons=protonation.txt watertype=fb3mod flexwater=0 \
  refpoints=initial_restraints.pdb restraint_mult=1 \
  pdbscale=0.01 gamma_ln=1.0 barostat=1 \
  leapfile=tleap_stub.in padwater=$padwater \
  cool_ns=0.001 heat_ns=0.5 equi_ns=0.5 prod_ns=0.5 \
  cool_slowdown=5 heat_slowdown=1 equi_slowdown=1 \
  restrain_omega=0 omega_weight=0 chiral_weight=0 \
  debug=1 >&! leap2amber_${itr}.log
if (! -e Prod.rst7) then
  echo "ERROR: leap2amber.com did not produce Prod.rst7 — check leap2amber_${itr}.log"
  goto exit
endif
# Produces: Prod.rst7  xtal.prmtop  padded.parm7  orignames.pdb
#           chir_omega0.rst  Bfac_0.pdb  restraints_for_0.pdb

cd ..


#==============================================================================
# CONVERGENCE POINT
# Systems diverge from PDB download through leap2amber.com (different ligands,
# crystal forms, MTZ columns). After leap2amber.com, all systems use the same
# optimize_weights_runme.com structure with only these parameter differences:
#
#   System         min_lig_weight   halfrho_neg
#   6c2r / EG7     0.1              3.5           (this script)
#   6c2r / AMPPNP  0.01             3.5
#   1aho           0                3.5
#   2qpx           TBD              TBD
#==============================================================================


#==============================================================================
# SECTION 11 — Weight optimization setup (opt1)
#==============================================================================
# Set these to point to the last good iteration from amber1 (or a previous opt):
set previtr  = 0          # last good iteration number (0 = use Prod directly)
set prevdir  = amber1     # directory containing that iteration
set prevprod = Prod       # file stem of the amber rst7/in/out/nc to continue from

set o = 1
mkdir opt${o}
if (-e compute_settings.sourceme) ln -sf ../compute_settings.sourceme opt${o}/
cd opt${o}
cp ${pdir}/optimize_weights_runme.com .

cp ../${prevdir}/restraints_for_${previtr}.pdb current_restraints.pdb
cp current_restraints.pdb restraints_for_0.pdb
cp ../${prevdir}/${prevprod}.rst7 amber_0.rst7
cp ../${prevdir}/${prevprod}.in  amber_0.in
cp ../${prevdir}/${prevprod}.out amber_0.out
cp ../${prevdir}/barometer_${previtr}.out barometer_0.out
cp ../${prevdir}/leap2amber_${previtr}.log .
ln -sf ../${prevdir}/${prevprod}.nc amber_0.nc
cp ../centroids/centroids_in_density.pdb all_possible_refpoints.pdb
cp ../${prevdir}/xtal.prmtop .
cp ../${prevdir}/padded.parm7 .
cp ../${prevdir}/orignames.pdb .
cp ../${prevdir}/Bfac_${previtr}.pdb Bfac.pdb
cp ../${prevdir}/chir_omega0.rst .
cp chir_omega0.rst chir_omega.rst
cp ../xtal_properties.sourceme .

# Stage 1: hydrate and settle pressure (no weight changes yet)
optimize_weights_runme.com prod_ns=0.5 max_mult=1 Bfac_maxmod=0 weight_power=1 \
    teleport_waters=1 hydrate_itr=1 add_radius=1.8 \
    pressure_avglast=auto pressure_scale=1,auto void_scale=1 \
    release_itr=0 repick_itr=0 \
    min_lig_weight=0.1 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=1 weight_negscale=1 randel_itr=0 \
    min_align_weight=5 maxitr=20 >&! runme1.log
if ($status) then
  echo "ERROR: optimize_weights Stage 1 failed — check runme1.log"
  goto exit
endif
if (! -e fofc_Rplot.txt) then
  echo "ERROR: optimize_weights Stage 1 produced no fofc_Rplot.txt"
  goto exit
endif

# Stage 2: scale up max_mult — let pressure re-equilibrate before enabling Bfac_maxmod
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=0 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0.1 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=0.9 weight_negscale=0.5 randel_itr=0 \
    min_align_weight=1 align_target=centroids align_nstlim=0 \
    halfrho_neg=3.5 halfrho_pos=auto maxitr=20 >&! runme2.log
if ($status) then
  echo "ERROR: optimize_weights Stage 2 failed — check runme2.log"
  goto exit
endif

# Backtrack recipe — use if Rfree climbs significantly:
# set lastitr = `tail -n 1 fofc_Rplot.txt | awk '{print $1}'`
# set gooditr = `tail -n 100 fofc_Rplot.txt | awk '{print $1,$2-f;f+=0.01}' | sort -k2g | awk '{print $1}'`
# cp -p restraints_for_${gooditr}.pdb current_restraints.pdb
# cp -p restraints_for_${gooditr}.pdb restraints_for_${lastitr}.pdb
# cp amber_${gooditr}.rst7 amber_${lastitr}.rst7
# cp amber_${gooditr}.in   amber_${lastitr}.in
# cp amber_${gooditr}.out  amber_${lastitr}.out
# cp amber_${gooditr}.nc   amber_${lastitr}.nc
# cp barometer_${gooditr}.out barometer_${lastitr}.out

cd ..


#==============================================================================
# SECTION 12 — Continue optimization (opt2 template)
# Repeat this block for opt3, opt4, ... adjusting parameters as needed.
#==============================================================================
set prevdir = opt1
if (! -e ${prevdir}/fofc_Rplot.txt) then
  echo "ERROR: ${prevdir}/fofc_Rplot.txt not found — opt1 did not complete"
  goto exit
endif
set previtr = `sort -k2g ${prevdir}/fofc_Rplot.txt | awk 'NR==1{print $1}'`
echo "opt2 continuing from ${prevdir} iteration $previtr (best Rfree)"
set prevprod = amber_${previtr}

set o = 2
mkdir opt${o}
if (-e compute_settings.sourceme) ln -sf ../compute_settings.sourceme opt${o}/
cd opt${o}
cp ${pdir}/optimize_weights_runme.com .

cp ../${prevdir}/restraints_for_${previtr}.pdb current_restraints.pdb
cp current_restraints.pdb restraints_for_0.pdb
cp ../${prevdir}/${prevprod}.rst7 amber_0.rst7
cp ../${prevdir}/${prevprod}.in  amber_0.in
cp ../${prevdir}/${prevprod}.out amber_0.out
cp ../${prevdir}/barometer_${previtr}.out barometer_0.out
cp ../${prevdir}/leap2amber_0.log .
ln -sf ../${prevdir}/${prevprod}.nc amber_0.nc
cp ../centroids/centroids_in_density.pdb all_possible_refpoints.pdb
cp ../${prevdir}/xtal.prmtop .
cp ../${prevdir}/padded.parm7 .
cp ../${prevdir}/orignames.pdb .
cp ../${prevdir}/Bfac_${previtr}.pdb Bfac.pdb
cp ../${prevdir}/chir_omega0.rst .
cp chir_omega0.rst chir_omega.rst
cp ../xtal_properties.sourceme .

# Typical opt2 parameters: moderate repicking, gentle down-weighting
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0.1 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=0.95 weight_negscale=0.9 randel_itr=0 \
    min_align_weight=0.5 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto maxitr=20 >&! runme1.log
if ($status) then
  echo "ERROR: optimize_weights opt2 failed — check runme1.log"
  goto exit
endif

cd ..

# Three-stage opt workflow:
#   Stage 1 (opt1 Stage 1): max_mult=1, Bfac_maxmod=0 — settle pressure and waters
#   Stage 2 (opt1 Stage 2): max_mult=2, Bfac_maxmod=0 — scale up weights, re-equilibrate pressure
#   Stage 3 (opt2+):        max_mult=2, Bfac_maxmod=1 — enable B-factor modification
#
# Goal: minimum Rfree across all opt directories and iterations.
# Pick the best iteration:
#   cat opt*/fofc_Rplot.txt | sort -k2g | head
#
# Continue with opt3, opt4 ... using the opt2 block as a template.
# Tuning guidance:
#   - Turn off void_scale (void_scale=0) once pressure is stable
#   - Add align_nstlim=250000 if restraint weights are oscillating

exit:
