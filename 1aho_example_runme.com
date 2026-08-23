#! /bin/tcsh -f
# 1aho_example_runme.com
# Scorpion toxin II (PDB 1aho) — single unit cell (super_mult=1,1,1) crystallographic Amber MD
# The bespoke_amber_restraints reference DEMO system: 64 residues, 0.965 Å,
# P212121, no protein ligand, all data public.  The starting model is the PDB
# deposit itself (no private collaborator model).
#
# Runs standalone with NO starter kit: data + model come from the PDB download
# (Section 1), salt (NH4, acetate) force fields from elbow/antechamber, HIS by QM.
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
#   1aho           0                3.5           this script; no ligand, single cell
#   6c2r / EG7     0.1              3.5           see 6c2r_EG7_1x_example_runme.com
#   6c2r / AMPPNP  0.01             3.5           Dirk's private dataset; badlinks=1
#
# Optimization stages (each seeds from the previous stage's last iteration and
# runs in the background until a monitor detects its convergence signal):
#   opt1 s1/s2 (Sec 11): rough hydration, then pressure-stable
#   opt2       (Sec 12): weight optimization  (exit: RESTRAINT energy flat)
#   opt3       (Sec 13): subtle B-factor mods (exit: fofc_R level) - only after
#                        weights converge
#   opt4+      (Sec 14): continuation of opt3 (exit: fofc_R level)
# See the "Opt-stage workflow" header above Section 13 for the full recipe.
#
# Re-run safety: each section checks for its key output file and skips if already done.
# This makes re-running safe after a crash — only incomplete sections re-execute.


#==============================================================================
# Environment setup
#==============================================================================
set pdbid     = 1aho
set super_mult = 1,1,1
set salt      = ( NH4 ACY )
set salt_conc = 0.68
set badlinks  = 0      # 1aho has no ligand, so no spurious Lys-LIG bonds

set pwd  = `pwd`
set pdir = ${pwd}/bespoke_amber_restraints
set skit = ${pwd}/starter_kit   # starter kit location
set t    = tempfile
set path = ( $pdir $path )

# Compile utility binaries into pdir if not already in PATH
which float_func >& /dev/null
if( $status && -e ${pdir}/float_func.c ) then
  echo "Compiling float_func and float_add..."
  gcc -O -o ${pdir}/float_func ${pdir}/float_func.c -lm
  gcc -O -o ${pdir}/float_add  ${pdir}/float_add.c  -lm
endif

# The starter kit is OPTIONAL — it only caches ligand/salt parameters and the HIS
# protonation, all of which this script can regenerate.  If it is absent we run
# fully standalone (data + model come from the PDB download in Section 1).
if (! -d $skit) then
cat << EOF
NOTE: no starter kit at $skit — running fully standalone.
  Salt (NH4, ACY) force-field files will be generated on the fly (elbow if it
  works, otherwise antechamber AM1-BCC), and HIS protonation is determined by QM.
  1aho has no protein ligand, so there is nothing else to parameterize.  Setup is
  a bit slower than with a kit; to use one, symlink it and point skit at it:
    ln -sf ~/projects/git/bespoke_amber_restraints/1aho_starter_kit starter_kit
EOF
endif

# compute-environment defaults (override these in compute_settings.sourceme)
set sruncpu = "srun"                # command prefix for CPU / serial jobs
set srungpu = "srun --gres=gpu:1"   # command prefix for GPU jobs

if (-e compute_settings.sourceme) source compute_settings.sourceme

if (! $?AMBERHOME) then
  # amber.csh derives AMBERHOME from $_, which is EMPTY in a non-interactive
  # tcsh -f script, so it would set AMBERHOME=cwd and leave tleap/sander off
  # PATH.  cd into the amber dir first so its fallback invocationpath ('.')
  # resolves to the real amber root.
  set _here = $cwd
  cd /programs/amber22 && source amber.csh
  cd $_here
endif
if (! $?PHENIX)    source /programs/phenix/phenix_env.csh
if (! $?CCP4)      source /programs/ccp4-9/bin/ccp4.setup-csh

if (-e user_settings.sourceme) source user_settings.sourceme

set debug = 1


#==============================================================================
# SECTION 1 — Download and sequence
#==============================================================================
if (-e ${pdbid}.pdb && -e seq.fasta) then
  echo ""
  echo "=== Section 1: already done, skipping ==="
  set rawmtz = `ls -1t ${pdbid}*.mtz | head -1`
  set modulo = `awk 'NR>1{printf("%s",$1)}' seq.fasta | wc -c`
  goto sec1done
endif
echo ""
echo "=== Section 1: Downloading $pdbid from PDB ==="
# Download structure factors and coordinates from PDB (replaces deprecated getcif.com)
phenix.fetch_pdb $pdbid action=all >&! fetch_pdb.log
# Produces: 6c2r.pdb  6c2r-sf.cif
echo "  fetched: `ls ${pdbid}.pdb ${pdbid}-sf.cif |& grep -v 'No such file'`"

phenix.cif_as_mtz ${pdbid}-sf.cif >&! cif_as_mtz.log
# Produces a .mtz file — name varies by phenix version; grab the newest one:
set rawmtz = `ls -1t ${pdbid}*.mtz | head -1`
if ($rawmtz == "") then
  echo "ERROR: phenix.cif_as_mtz produced no .mtz file — details: cif_as_mtz.log"
  goto exit
endif
echo "  raw MTZ: $rawmtz"

# Extract protein sequence from SEQRES
grep SEQRES ${pdbid}.pdb |\
sequence.awk |\
awk 'NR==1{$0="> "$0} {print} NF==0{exit}' |\
cat >! seq.fasta

set modulo = `awk 'NR>1{printf("%s",$1)}' seq.fasta | wc -c`
echo "  modulo (residues per chain) = $modulo"

# Use PDB deposit as starting ASU coordinates
cp ${pdbid}.pdb starthere_asu.pdb
sec1done:


#==============================================================================
# SECTION 2 — Reflection data and crystal properties
#==============================================================================
if (-e refme.mtz && -e xtal_properties.sourceme) then
  echo ""
  echo "=== Section 2: already done, skipping ==="
  source xtal_properties.sourceme
  # Re-derive supercell CELL parameters needed by Section 6 tleap
  set CELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' mtzdump.txt`
  goto sec2done
endif
echo ""
echo "=== Section 2: Reflection data and crystal properties ==="
# Auto-detect structure factor amplitude column from the raw MTZ
set F    = `mtzdmp $rawmtz | awk 'NF>10 && ! /2FOFCWT|FWT|DELF/ && $(NF-1)=="F"{print $NF;exit}'`
set SIGF = `mtzdmp $rawmtz | awk 'NF>10 && ! /2FOFCWT|FWT|DELF/ && $(NF-1)=="Q"{print $NF;exit}'`
set FREE = `mtzdmp $rawmtz | awk 'NF>10 && $(NF-1)=="I" && /[Ff]ree/{print $NF;exit}'`
echo "Using columns: F=$F  SIGF=$SIGF  FREE=$FREE"
if ($F == "" || $SIGF == "") then
  echo "ERROR: could not auto-detect F/SIGF MTZ columns — check mtzdmp output above"
  goto exit
endif

# Select F/SIGF (plus the deposit's R-free flags, if it has any) into refme_small.mtz
if ( "$FREE" != "" ) then
  cad hklin1 $rawmtz hklout refme_small.mtz << EOF > /dev/null
labin file 1 E1=$F E2=$SIGF E3=$FREE
labou file 1 E1=FP E2=SIGFP E3=FreeR_flag
EOF
else
  cad hklin1 $rawmtz hklout refme_small.mtz << EOF > /dev/null
labin file 1 E1=$F E2=$SIGF
labou file 1 E1=FP E2=SIGFP
EOF
endif
if ($status) then
  echo "ERROR: cad failed — check column names above"
  goto exit
endif

# Ensure a usable R-free set.  Older deposits (e.g. 1aho) carry no R-free flags, or
# a uniform column (all one value = no test set) that phenix.refine later rejects
# ("No array of R-free flags found").  If FreeR_flag is absent or uniform, generate
# a fresh 5% test set — CCP4 freerflag if available, else phenix.  Real test sets
# are left untouched.
set frcol = `mtzdmp refme_small.mtz | awk '$NF=="FreeR_flag" && $(NF-1)=="I"{print $3,$4;exit}'`
set frregen = 1
if ( $#frcol == 2 ) then
  if ( "$frcol[1]" != "$frcol[2]" ) set frregen = 0
endif
if ( $frregen ) then
  echo "  R-free flags absent or uniform — generating a fresh 5% test set"
  cad hklin1 refme_small.mtz hklout ${t}noflag.mtz << EOF > /dev/null
labin file 1 E1=FP E2=SIGFP
EOF
  which freerflag >& /dev/null
  if ( $status == 0 ) then
    freerflag hklin ${t}noflag.mtz hklout refme_small.mtz << EOF >&! freerflag.log
FREERFRAC 0.05
END
EOF
  else
    phenix.reflection_file_converter ${t}noflag.mtz --generate-r-free-flags \
      --r-free-flags-fraction=0.05 --mtz=refme_small.mtz >&! freerflag.log
  endif
  rm -f ${t}noflag.mtz
  if ( ! -e refme_small.mtz ) then
    echo "ERROR: R-free flag generation failed — see freerflag.log"
    goto exit
  endif
endif

echo head | mtzdump hklin refme_small.mtz >! smallmtzdump.txt
set SGnum   = `awk '/Space group =/{print $NF+0}' smallmtzdump.txt | tail -n 1`
set smallSG = `awk -v num=$SGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
set smallCELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' smallmtzdump.txt`
set nsymops = `awk -v SG=$smallSG '$4==SG {print $2}' ${CLIBD}/symop.lib | head -1`
set reso    = `awk '/Resolution Range :/{getline;getline;print $(NF-2)+0;exit}' smallmtzdump.txt | head -1`

# Write crystal properties file — this will be linked into each subdirectory
# and auto-sourced by generate_alignment_reference_runme.com to pass badlinks
if(-e xtal_properties.sourceme) then
  echo "keeping existing xtal_properties.sourceme"
else
  cat << EOF >! xtal_properties.sourceme
set reso       = $reso
set smallSG    = $smallSG
set smallCELL  = ( $smallCELL )
set nsymops    = $nsymops
set super_mult = $super_mult
set modulo     = $modulo
set ligands    = ""
set salt       = ( $salt )
set salt_conc  = $salt_conc
set badlinks   = $badlinks
EOF
endif
source xtal_properties.sourceme

# computing environment stuff goes here
if(-e compute_settings.sourceme) then
  echo "keeping existing compute_settings.sourceme"
else
  cat << EOF >! compute_settings.sourceme
set pdir       = $pdir
set sruncpu    = "$sruncpu"
set srungpu    = "$srungpu"
EOF
endif
source compute_settings.sourceme

# Expand MTZ to the crystallographic cell (super_mult=1,1,1 = identity reindex, no expansion)
echo "  expanding MTZ to $super_mult supercell..."
cad hklin1 refme_small.mtz hklout expanded.mtz << EOF >&! mtz_expand.log
labin file 1 all
outlime space 1
EOF
cad hklin1 expanded.mtz hklout refme_cell.mtz << EOF >>&! mtz_expand.log
labin file 1 all
symm 1
EOF
rm -f expanded.mtz

set reidx = `echo $super_mult | awk -F "[ ,x]" '{print "reindex h"$1",k"$2",l"$3}'`
echo $reidx |\
reindex hklin refme_cell.mtz hklout refme.mtz >>&! mtz_expand.log
if ($status) then
  echo "ERROR: MTZ expansion to supercell failed — details: mtz_expand.log"
  goto exit
endif
rm -f refme_cell.mtz

echo head | mtzdump hklin refme.mtz >! mtzdump.txt
set CELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' mtzdump.txt`

echo "" >! blank.pdb
echo "CELL $CELL\nSPACE 1" | pdbset xyzin blank.pdb xyzout ${t}cell.pdb > /dev/null
egrep "^CRYST1" ${t}cell.pdb >! cell.pdb
rm -f blank.pdb ${t}cell.pdb
echo "  cell: $CELL   SG: $smallSG   reso: $reso"
sec2done:


#==============================================================================
# SECTION 3 — Ligand force-field files
#==============================================================================
# Resolve the ligand list BEFORE the skip guard so re-runs that already have the
# FF files still set $ligands (and heal a stale "auto" in xtal_properties.sourceme)
# for the downstream sections that copy ${lig}.cif/pdb/mol2/frcmod (e.g. amber_asu).
if( "$ligands" == "auto" ) then
  # Auto-detect ligands from PDB, excluding salt ions and water
  set ligands = `filter_pdb.awk -v only=ligand,atoms starthere_asu.pdb |\
    awk '/^HETAT/{print substr($0,18,3)}' | sort -u |\
    awk -v salt="$salt HOH WAT" \
      'BEGIN{n=split(salt,s);for(i=1;i<=n;i++)bad[s[i]]=1} {if($1 in bad)next; print}'`
  echo "Ligands found: $ligands"
endif
if( "$ligands" != "" ) then
  if( `grep -c "^set ligands.*auto" xtal_properties.sourceme` ) then
    sed -i "s|^set ligands.*|set ligands    = $ligands|" xtal_properties.sourceme
    echo "Updated xtal_properties.sourceme: ligands = $ligands"
  endif
endif

if (-e ligands/NH4.mol2 && -e ligands/ACY.mol2) then
  echo ""
  echo "=== Section 3: already done, skipping ==="
  goto sec3done
endif
echo ""
echo "=== Section 3: Ligand force-field files ==="
mkdir -p ligands

# 1aho has no protein ligand (ligands=""), so nothing is fetched here; only the
# crystallisation salt (NH4, ACY) needs force-field files, built by the loop below.

cd ligands


# Antechamber binary for AM1-BCC charge fallback.
# Elbow --opt fails on Phenix 2.0+ (RM1 semiempirical > 1 GB virtual memory limit).
# Elbow --amber is broken in Phenix 2.1rc2 (amber_adaptbx conda_base path is wrong;
#   fix: change '../conda_base' to '../..' in AmberPrepClass.py, or apply the nightly patch).
# Preferred: phenix-bundled antechamber ($PHENIX/bin) — uses am1bcc, no sqm needed.
# Fallback: amber18 antechamber (amber22 sqm fails: missing libgfortran.so.4).
set antechamber_cmd = ""
set parmchk2_cmd = ""
foreach ab ( ${PHENIX}/bin /home/programs/amber18/bin /programs/amber18/bin )
  if (-e $ab/antechamber) then
    set antechamber_cmd = $ab/antechamber
    set parmchk2_cmd = $ab/parmchk2
    break
  endif
end

touch elbow.log
foreach lig ( $ligands )
  if (-e ${lig}.mol2 && -e ${lig}.pdb && -e ${lig}.cif && -e ${lig}.frcmod) continue
  echo "building force field for $lig (elbow first, then antechamber; verbose output and non-fatal failures go to elbow.log)"
  # Try elbow first - it may now succeed (a newer Phenix, and --amber_force_field_files
  # needs amber on PATH, which is now set up correctly).  It runs from the starter
  # kit's correct-protonation .cif.  If every elbow tier fails, the loop falls
  # through to the antechamber AM1-BCC path below (which uses the kit .mol2/.cif).
  # 2. Elbow with --opt (full QM; ~16 GB; works on Phenix dev-5353 and earlier)
  if (-e ${lig}.cif) then
    echo "elbow $lig from cif..."
    phenix.elbow ${lig}.cif --id=${lig} --opt --opt_nproc=10 --amber_force_field_files >>& elbow.log
  endif
  if (-e ${lig}.mol2 && -e ${lig}.pdb && -e ${lig}.cif && -e ${lig}.frcmod) continue
  if (-e ${lig}.cif ) then
    echo "elbow $lig from cif with more memory..."
    phenix.elbow ${lig}.cif --id=${lig} --opt --opt_nproc=10 --amber_force_field_files --memory=128G >>& elbow.log
  endif
  echo "elbow $lig from scratch"
  phenix.elbow --chemical_component $lig --id=${lig} --opt --opt_nproc=10 \
    --amber_force_field_files >>& elbow.log
  if (-e ${lig}.mol2 && -e ${lig}.pdb && -e ${lig}.cif && -e ${lig}.frcmod) continue
  echo "elbow $lig with more memory"
  phenix.elbow --chemical_component $lig --id=${lig} --opt --opt_nproc=10 \
    --amber_force_field_files --memory=128G >>& elbow.log
  if (-e ${lig}.mol2 && -e ${lig}.pdb && -e ${lig}.cif && -e ${lig}.frcmod) continue
  if (! -e ${lig}.pdb && -e ${lig}.mol2) then
    echo "elbow to get pdb from ${lig}.mol2"
    phenix.elbow ${lig}.mol2 --id=${lig} >>& elbow.log
  endif
  if (! -e ${lig}.pdb) then
    echo "elbow to get pdb from scratch"
    phenix.elbow --chemical_component $lig --id=${lig} >>& elbow.log
  endif
  # 3. Antechamber AM1-BCC fallback (bypasses Phenix RM1 limit)
  # Elbow geometry-only (no --opt) produces .pdb and .cif without QM charges.
  if (-e ${lig}.pdb && "$antechamber_cmd" != "") then
    # Net charge is essential: sqm rejects an odd-electron system, so a wrong nc
    # (the old default 0 for a carboxylate like EG7) is a hard failure.  Read nc
    # from the local .cif; if elbow wrote none, fall back to the CCP4 monomer
    # library ($CLIBD_MON/<l>/<LIG>.cif), which carries the formal charge.
    set nccif = ${lig}.cif
    if (! -e $nccif && $?CLIBD_MON) then
      set l1 = `echo $lig | awk '{print tolower(substr($1,1,1))}'`
      if (-e ${CLIBD_MON}/${l1}/${lig}.cif) set nccif = ${CLIBD_MON}/${l1}/${lig}.cif
    endif
    set nc = ""
    if (-e "$nccif") set nc = `awk 'BEGIN{in_a=0;sum=0} /^loop_/{in_a=0} /_chem_comp_atom/{in_a=1} \
      in_a && /^[A-Z][A-Z0-9]* / && NF>=5{sum+=$5} \
      END{print int(sum<0?sum-0.5:sum+0.5)}' $nccif`
    if ("$nc" == "") set nc = 0   # no parseable .cif anywhere -> assume neutral
    echo "Generating ${lig}.mol2 via antechamber AM1-BCC (net charge $nc)"
    ${antechamber_cmd} -i ${lig}.pdb -fi pdb -o ${lig}.mol2 -fo mol2 \
      -c bcc -nc $nc -at gaff2 > ${lig}_antechamber.log
    if (-e ${lig}.mol2) ${parmchk2_cmd} -i ${lig}.mol2 -f mol2 -o ${lig}.frcmod
  endif
  # 4. Last resort: phenix geostd ships pre-built amber force fields (mol2+frcmod,
  # GAFF2, correct BCC charges) for many ligands.  Used ONLY if elbow AND
  # antechamber both failed to leave a usable mol2 + non-empty frcmod (e.g. EG7:
  # antechamber makes the mol2 but parmchk2 hits a dummy atom and writes an empty
  # frcmod).  CAUTION: geostd atom/H names can differ from the PDB deposit, and a
  # name scramble blows up the amber heating step - so this is the fallback, not
  # the default.  (For EG7 the geostd names match the deposit, verified.)
  if( (! -e ${lig}.mol2 || ! -e ${lig}.frcmod || -z ${lig}.frcmod) && $?PHENIX ) then
    set l1 = `echo $lig | awk '{print tolower(substr($1,1,1))}'`
    set gmol2 = `ls ${PHENIX}/lib/python*/site-packages/chem_data/geostd/${l1}/${lig}.mol2 ${PHENIX}/modules/chem_data/geostd/${l1}/${lig}.mol2 |& grep -v 'o such' | head -1`
    if( "$gmol2" != "" ) then
      set gdir = $gmol2:h
      if( -e ${gdir}/${lig}.frcmod && ! -z ${gdir}/${lig}.frcmod ) then
        echo "elbow+antechamber failed for $lig; taking amber force field from phenix geostd ($gdir)"
        echo "  NOTE: verify geostd atom/H names match the deposit - a mismatch can break amber heating"
        cp ${gdir}/${lig}.mol2 ${lig}.mol2
        cp ${gdir}/${lig}.frcmod ${lig}.frcmod
        if( ! -e ${lig}.cif && -e ${gdir}/data_${lig}.cif ) cp ${gdir}/data_${lig}.cif ${lig}.cif
      endif
    endif
  endif
  if (! -e ${lig}.mol2 || ! -e ${lig}.frcmod || -z ${lig}.frcmod) then
    echo "ERROR: could not build a usable force field for ${lig}"
    echo "  Tried: starter kit, elbow, antechamber AM1-BCC, phenix geostd"
    echo "  This ligand may need curated params - drop ${lig}.mol2/.frcmod into a"
    echo "  starter kit, or use an earlier phenix whose elbow works (e.g. dev-5353)."
    echo "  Check: elbow.log, ${lig}_antechamber.log"
    goto exit
  endif
end
# Per-ligand force-field outcome - runs for every ligand no matter which path it
# took above, so a noisy (but non-fatal) elbow traceback in elbow.log is never
# mistaken for failure of the ligand.  READY means amber has a usable mol2+frcmod.
foreach lig ( $ligands )
  set havefiles = ""
  foreach ext ( mol2 frcmod cif pdb )
    if (-e ${lig}.$ext) set havefiles = "$havefiles $ext"
  end
  if (-e ${lig}.mol2 && -e ${lig}.frcmod) then
    echo "  ==> $lig force field READY:$havefiles"
  else
    echo "  ==> $lig force field INCOMPLETE:$havefiles  (see elbow.log / ${lig}_antechamber.log)"
  endif
end

# Salt ions — use starter kit files (mol2/frcmod/pdb); fall back to antechamber.
# add_salt_runme.com needs ${ion}.pdb in the working directory (not ligands/).
# Elbow --opt fails on Phenix 2.0+ (RM1 semiempirical > 1 GB virtual memory limit).
foreach lig ( $salt )
  if (! -e ${lig}.mol2 && -e ${skit}/ligands/${lig}.mol2) then
    cp ${skit}/ligands/${lig}.mol2 .
    cp ${skit}/ligands/${lig}.frcmod .
    if (-e ${skit}/ligands/${lig}.pdb) cp ${skit}/ligands/${lig}.pdb .
  endif
  if (-e ${lig}.mol2) continue
  # Elbow geometry-only to get .pdb + .cif, then antechamber for charges
  if (! -e ${lig}.pdb) then
    phenix.elbow --chemical_component $lig --id=${lig} > ${lig}_elbow_geom.log
  endif
  if (-e ${lig}.pdb && "$antechamber_cmd" != "") then
    # Net charge is essential: sqm rejects an odd-electron system, so a wrong nc
    # (the old default 0 for a carboxylate like EG7) is a hard failure.  Read nc
    # from the local .cif; if elbow wrote none, fall back to the CCP4 monomer
    # library ($CLIBD_MON/<l>/<LIG>.cif), which carries the formal charge.
    set nccif = ${lig}.cif
    if (! -e $nccif && $?CLIBD_MON) then
      set l1 = `echo $lig | awk '{print tolower(substr($1,1,1))}'`
      if (-e ${CLIBD_MON}/${l1}/${lig}.cif) set nccif = ${CLIBD_MON}/${l1}/${lig}.cif
    endif
    set nc = ""
    if (-e "$nccif") set nc = `awk 'BEGIN{in_a=0;sum=0} /^loop_/{in_a=0} /_chem_comp_atom/{in_a=1} \
      in_a && /^[A-Z][A-Z0-9]* / && NF>=5{sum+=$5} \
      END{print int(sum<0?sum-0.5:sum+0.5)}' $nccif`
    if ("$nc" == "") set nc = 0   # no parseable .cif anywhere -> assume neutral
    echo "Generating ${lig}.mol2 via antechamber AM1-BCC (net charge $nc)"
    ${antechamber_cmd} -i ${lig}.pdb -fi pdb -o ${lig}.mol2 -fo mol2 \
      -c bcc -nc $nc -at gaff2 > ${lig}_antechamber.log
    if (-e ${lig}.mol2) ${parmchk2_cmd} -i ${lig}.mol2 -f mol2 -o ${lig}.frcmod
  else
    # Last resort: elbow --opt (fails on Phenix 2.0+, works on dev-5353 and earlier)
    phenix.elbow --chemical_component $lig --id=${lig} --opt --opt_nproc=10 \
      --amber_force_field_files
  endif
  if (! -e ${lig}.mol2) then
    echo "ERROR: could not generate ${lig}.mol2 for salt ion $lig"
    echo "  Copy ${lig}.mol2, ${lig}.frcmod, ${lig}.pdb from starter kit or another source"
    goto exit
  endif
end

cd ..
sec3done:


#==============================================================================
# SECTION 4 — HIS protonation (QM-determined)
#==============================================================================
if (-e HIS_settings_asu.txt) then
  echo ""
  echo "=== Section 4: already done, skipping ==="
  goto sec4done
endif
echo ""
echo "=== Section 4: HIS protonation ==="
# HIS protonation (HID/HIE/HIP) is determined by QM with the Phenix quantum
# interface.  Pre-computed results live in the starter kit as HIS_settings_asu.txt;
# if that file is absent, this section regenerates it by running
# HIS_protonation_QM_runme.com in a his_qm/ subdir (SLOW - QM refinement of every
# histidine).  To (re)generate by hand for a new system:
#   HIS_protonation_QM_runme.com pdbfile=starthere_asu.pdb nproc=6
#   # 1aho has no protein ligand, so no ligcifs=; add ligcifs=<lig>.cif for systems
#   # that do.  Writes HIS_settings.txt -> copy here as HIS_settings_asu.txt
#
# Format: one line per HIS — HIE (Nε), HID (Nδ), HIP (both / positively charged)
# 1aho result: HIE A54  HIE A64  (both epsilon-protonated)
if (-e ${skit}/HIS_settings_asu.txt) then
  cp ${skit}/HIS_settings_asu.txt HIS_settings_asu.txt
else
  # not in the kit - determine the states by QM in a subdir (SLOW)
  echo "  HIS_settings_asu.txt not in starter kit - computing HIS protonation by QM (slow)..."
  mkdir -p his_qm
  cd his_qm
  cp ../starthere_asu.pdb .
  cp ../ligands/*.cif . >& /dev/null
  HIS_protonation_QM_runme.com pdbfile=starthere_asu.pdb nproc=6 >&! his_qm.log
  cd ..
  if (-e his_qm/HIS_settings.txt && ! -z his_qm/HIS_settings.txt) then
    cp his_qm/HIS_settings.txt HIS_settings_asu.txt
    echo "  QM HIS protonation result:"
    cat HIS_settings_asu.txt
  else
    echo "ERROR: QM HIS protonation failed - see his_qm/his_qm.log"
    echo "       or drop a pre-computed HIS_settings_asu.txt into the starter kit"
    goto exit
  endif
endif
sec4done:


#==============================================================================
# SECTION 5 — Build and sanitize starting model (build1)
#==============================================================================
if (-e build1/minRfree.pdb) then
  echo ""
  echo "=== Section 5: already done, skipping ==="
  goto sec5done
endif
echo ""
echo "=== Section 5: Building starting model (build1) ==="
mkdir -p build1
cd build1

ln -sf ../starthere_asu.pdb starthere.pdb
ln -sf ../refme_small.mtz refme.mtz
ln -sf ../xtal_properties.sourceme .
foreach lig ( $ligands $salt )
  cp ../ligands/${lig}.cif .
end
set ligcifs = `echo $ligands | awk '{for(i=1;i<=NF;++i) print $i ".cif"}'`

echo "  building missing atoms (details: buildout.log)..."
buildout_pdb_runme.com starthere.pdb $ligcifs badlinks=$badlinks >&! buildout.log
if (! -e built.pdb) then
  echo "ERROR: buildout_pdb_runme.com did not produce built.pdb — details: buildout.log"
  goto exit
endif

# Iterative refinement to convergence — produces refmacout_minRfree.pdb
echo "  refmac convergence refinement (details: converge.log)..."
converge_refmac.com built.pdb refme.mtz $ligcifs >&! converge.log
if (! -e refmacout_minRfree.pdb) then
  echo "ERROR: converge_refmac.com did not produce refmacout_minRfree.pdb — details: converge.log"
  goto exit
endif
ln -sf refmacout_minRfree.pdb minRfree.pdb

grep "REMARK  FREE R VALUE" refmacout_minRfree.pdb | tail -1

ln -sf minRfree.pdb thisone.pdb

cd ..
sec5done:


#==============================================================================
# SECTION 6 — ASU tleap charge check (amber_asu)
#==============================================================================
if (-e amber_asu/cell_charge.txt) then
  echo ""
  echo "=== Section 6: already done, skipping ==="
  goto sec6done
endif
echo ""
echo "=== Section 6: ASU charge check (amber_asu) ==="
mkdir -p amber_asu
cd amber_asu

foreach lig ( $ligands $salt )
  foreach ext ( cif pdb mol2 frcmod )
    cp ../ligands/${lig}.${ext} .
  end
end
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

tleap -f tleap_stub.in >&! tleap.log
set tleap_status = $status

set charge0 = `awk '/unperturbed charge/{gsub(/[)(]/,"");print int($7);exit}' tleap.log`
set charge  = `echo $charge0 $nsymops | awk '{print $1*$2}'`
# fail loudly if tleap did not run / build a topology / report a charge -
# e.g. amber tools off PATH shows "tleap: Command not found" in tleap.log and
# leaves the ASU charge blank (do not silently continue with a bogus 0).
if ( $tleap_status || ! -s xtal.prmtop || "$charge0" == "" ) then
  echo "ERROR: tleap failed in amber_asu - no topology / ASU charge (see amber_asu/tleap.log)."
  echo "       Amber may be off PATH - check 'which tleap' and that AMBERHOME points to amber22."
  goto exit
endif
echo "  ASU charge $charge0, cell charge $charge" | tee cell_charge.txt

cd ..
sec6done:


#==============================================================================
# SECTION 7 — Generate centroid reference points (garr1)
#==============================================================================
if (-e garr1/centroids_final_001.pdb) then
  echo ""
  echo "=== Section 7: already done, skipping ==="
  goto sec7done
endif
echo ""
echo "=== Section 7: Generate centroid reference points (garr1) ==="
echo "  running generate_alignment_reference (details: garr1/garr.log)..."
mkdir -p garr1
cd garr1

ln -sf ../build1/thisone.pdb starthere.pdb
ln -sf ../refme_small.mtz refme.mtz
foreach lig ( $ligands $salt )
  foreach ext ( cif )
    cp ../ligands/${lig}.${ext} .
  end
end
ln -sf ../xtal_properties.sourceme .   # REQUIRED: passes badlinks to no_new_nonbonds_runme.com
if( ! $?ligcifs ) then
  set ligcifs = `echo $ligands | awk '{for(i=1;i<=NF;++i) print $i ".cif"}'`
  #set ligcifs = `cd ../ligands ; ls -1 *.cif |& awk '/.cif/'`
endif

# Note on badlinks (AMPPNP only): badlinks=1 suppresses a spurious Lys141-NZ to
# phosphate covalent bond that phenix creates by proximity. EG7 does NOT need it.
# If badlinks=1 is omitted for AMPPNP on the very first garr run, phenix_opts_unbump.eff
# will contain a NZ-LIG bond reference that crashes any subsequent restart when phenix
# can no longer find that bond. Fix: delete phenix_opts_unbump.eff AND start.geo, then
# restart from itr=1. The xtal_properties.sourceme link above prevents this.

generate_alignment_reference_runme.com starthere.pdb $ligcifs \
  repulse_nb=100 crush_nb=100 repulse_scale=0.5 >&! garr.log
if (! -e centroids_final_001.pdb) then
  echo "ERROR: garr did not produce centroids_final_001.pdb — details: garr.log"
  goto exit
endif
# Re-run (continuing from last itr) if Rfree has not converged.

cd ..
sec7done:


#==============================================================================
# SECTION 8 — Build supercell centroid set (centroids/)
#==============================================================================
if (-e centroids/centroids_in_density.pdb) then
  echo ""
  echo "=== Section 8: already done, skipping ==="
  goto sec8done
endif
echo ""
echo "=== Section 8: Build supercell centroid set (centroids/) ==="
mkdir -p centroids
cd centroids

# Pick best map: lowest Rfree across build1 and garr1
grep "Final R" ../build1/*.log ../garr1/*.log | justify.awk | sort -k7g | tee sorted.txt | head
set minRfree = `awk '/_/{gsub(".log:"," ");print $1;exit}' sorted.txt`
ln -sf ${minRfree}.mtz minRfree.mtz
echo "min Rfree is $minRfree"

ln -sf ../garr1/centroids_final_001.pdb centroids_asu.pdb
ln -sf ../build1/thisone.pdb fulllength_asu.pdb
ln -sf ../refme_small.mtz .
foreach lig ( $ligands $salt )
  foreach ext ( cif )
    cp ../ligands/${lig}.${ext} .
  end
end
set ligcifs = `echo $ligands | awk '{for(i=1;i<=NF;++i) print $i ".cif"}'`
#set ligcifs = `ls -1 *.cif |& awk '/.cif/'`

# Verify centroid orientation vs full-length model
echo "flip-to-target RMSD  [ flipped.pdb (centroids re-oriented to target)  vs  fulllength_asu.pdb (target model) ]:"
flip_to_target_runme.com centroids_asu.pdb fulllength_asu.pdb > /dev/null
filter_pdb.awk -v skip=H,water fulllength_asu.pdb flipped.pdb | rmsd | head

# 2Fo-Fc reference density map
set Fwt  = `mtzdmp minRfree.mtz | awk 'NF>5 && /2FOFCWT|FWT/ && ! / DELF/ && $(NF-1)=="F"{print $NF;exit}'`
set PHIwt = `mtzdmp minRfree.mtz | awk 'NF>5 && /PH2FOFCWT|PHWT/ && ! /DEL/ && $(NF-1)=="P"{print $NF;exit}'`
cad hklin1 minRfree.mtz hklout reference0.mtz << EOF > /dev/null
labin file 1 E1=$Fwt E2=$PHIwt
labou file 1 E1=Fref E2=PHIref
EOF
fft hklin reference0.mtz mapout ffted.map << EOF > /dev/null
labin F1=Fref PHI=PHIref
EOF
mapmask mapin ffted.map mapout reference0.map << EOF > /dev/null
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
echo "combine check, expect ALL-ZERO and no warnings  [ new.pdb (centroids+fulllength merged)  vs  centroids_asu_noalt.pdb ]:"
combine_pdbs_runme.com centroids_asu_noalt.pdb fulllength_asu_noalt.pdb > /dev/null
filter_pdb.awk -v skip=water new.pdb centroids_asu_noalt.pdb | rmsd | head
# All RMS should be zero with no warnings

# Expand to supercell — run twice: first to discover monomer map, then apply it
echo "  expanding full-length ASU to supercell (details: fulllength_expand.log)..."
expand2supercell_runme.com fulllength_renamed.pdb refme_small.mtz super_mult=$super_mult \
  outprefix=fulllength_super phenix_bumpcheck=0 debug=1 >&! fulllength_expand0.log
cp fulllength_super.pdb fulllength_super0.pdb
cp new_monomer_rot_trans.txt monomer_rot_trans.txt

expand2supercell_runme.com fulllength_renamed.pdb refme_small.mtz super_mult=$super_mult \
  outprefix=fulllength_super phenix_bumpcheck=0 debug=1 \
  mono_map=monomer_rot_trans.txt >&! fulllength_expand.log

echo "  expanding centroids to supercell (details: centroids_expand.log)..."
expand2supercell_runme.com centroids_renamed.pdb refme_small.mtz super_mult=$super_mult \
  refpdb=fulllength_renamed.pdb outprefix=centroids_super \
  phenix_bumpcheck=0 debug=1 \
  mono_map=monomer_rot_trans.txt >&! centroids_expand.log

echo "supercell RMSD  [ centroids_super.pdb  vs  fulllength_super.pdb ]:"
filter_pdb.awk -v skip=H,water centroids_super.pdb fulllength_super.pdb | rmsd | head
echo "  CA-only  [ centroids_super.pdb  vs  fulllength_super.pdb ]:"
filter_pdb.awk -v skip=H,water centroids_super.pdb fulllength_super.pdb | grep CA | rmsd | head

foreach prefix ( fulllength centroids )
  awk '! /^ATOM|^HETAT/{print;next} {print substr($0,1,16),substr($0,18)}' \
    ${prefix}_super.pdb |\
  awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
  cat >! ${prefix}_noalt.pdb
end
combine_pdbs_runme.com centroids_noalt.pdb fulllength_noalt.pdb | head -n 1
echo "combine check, expect ALL-ZERO  [ new.pdb (centroids+fulllength merged, supercell)  vs  centroids_noalt.pdb ]:"
filter_pdb.awk -v skip=water new.pdb centroids_noalt.pdb | rmsd | head
# All RMS should be zero 

# Label centroids by their 2Fo-Fc density value; only atoms in density are kept
filter_pdb.awk -v skip=H centroids_super.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  {occ=substr($0,55,6)+0}\
  occ>0.01{print substr($0,1,80),"           |",$NF}' >! all_possible_centroids.pdb

echo "  labeling centroids by 2Fo-Fc density..."
rholabel_runme.com all_possible_centroids.pdb reference0.mtz mtzlabel=Fref >&! rholabel.log

cat rholabeled.pdb |\
awk '/^CRYST/{print;next} ! /^ATOM|^HETAT/{next}\
  {rho=$NF;pre=substr($0,1,60);post=substr($0,67)}\
  rho>0{printf("%s%6.2f%s\n",pre,rho,post)}' |\
cat >! centroids_in_density.pdb 
if (! -e centroids_in_density.pdb) then
  echo "ERROR: rholabel produced no centroids_in_density.pdb"
  goto exit
endif

echo "in-density RMSD  [ centroids_in_density.pdb  vs  fulllength_super.pdb (protein) ]:"
filter_pdb.awk -v skip=water,H centroids_in_density.pdb fulllength_super.pdb | rmsd | head
echo "  CA-only  [ centroids_in_density.pdb  vs  fulllength_super.pdb ]:"
filter_pdb.awk -v skip=H,water centroids_in_density.pdb fulllength_super.pdb | grep CA | rmsd | head

cd ..
sec8done:


#==============================================================================
# SECTION 9 — Supercell phenix refinement (super_refine1)
#==============================================================================
if (-e super_refine1/omegafix1_001.pdb) then
  echo ""
  echo "=== Section 9: already done, skipping ==="
  goto sec9done
endif
echo ""
echo "=== Section 9: Supercell phenix refinement (super_refine1) ==="
mkdir -p super_refine1
cd super_refine1

ln -sf ../refme.mtz .
ln -sf ../centroids/fulllength_super.pdb starthere.pdb
foreach lig ( $ligands $salt )
  foreach ext ( cif )
    cp ../ligands/${lig}.${ext} .
  end
end
if( ! $?ligcifs ) then
  set ligcifs = `echo $ligands | awk '{for(i=1;i<=NF;++i) print $i ".cif"}'`
  #set ligcifs = `cd ../ligands ; ls -1 *.cif |& awk '/.cif/'`
endif

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

# Check protein residue occ sums — catches blank-altloc partials and missing alternates.
# jigglepdb treats CA occ sum < 1.0 as a partial-occupancy conformer and may drop the residue.
filter_pdb.awk -v only=protein starthere.pdb |\
  awk '/^ATOM|^HETATM/ && substr($0,13,4)==" CA " {\
    res = substr($0,22,1) " " substr($0,23,4)+0 " " substr($0,18,3)\
    occ_sum[res] += substr($0,55,6)+0\
  }\
  END { for (r in occ_sum) if (occ_sum[r] < 0.99) print r }' | sort >! badocc.txt
if (-s badocc.txt) then
  echo "ERROR: protein residues with CA occ sum < 1.0 in starthere.pdb:"
  awk '{print "  " $0}' badocc.txt
  echo "Fix occupancies before proceeding (missing alternate conformer, or phenix occ artifact)"
  goto exit
endif

if (! -e phenix_001.pdb) then
  awk '{print substr($0,1,80)}' starthere.pdb >! refme.pdb
  echo "  phenix.refine pass 1 (details: super_refine1/phenix1.log)..."
  phenix.refine ../refme.mtz refme.pdb prefix=phenix opts.eff $ligcifs >&! phenix1.log
  if (! -e phenix_001.pdb) then
    echo "ERROR: phenix.refine failed — details: phenix1.log"
    goto exit
  endif
else
  echo "  phenix.refine pass 1: already done, skipping"
endif

if (! -e confsel_min.pdb) then
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

  # check that confsel has all protein atoms; fall back to fulllength_noalt if not
  filter_pdb.awk -v only=protein fulllength_noalt.pdb confsel.pdb | rmsd >! rmsd_confsel.txt
  if (`grep -c WARN rmsd_confsel.txt`) then
    echo "  WARNING: confsel missing atoms — using fulllength_noalt instead (see rmsd_confsel.txt)"
    cp fulllength_noalt.pdb confsel.pdb
  endif

  echo "  reorganizing conformers and waters..."
  reorganize_pdb_runme.com confsel.pdb refpdb=fulllength_noalt.pdb outfile=oneconf.pdb \
    phenix_bumpcheck=0 debug=1 autorerun=0 >&! reorg_oneconf.log

  reorganize_waters.com oneconf.pdb super_mult=$super_mult \
    smallmtz=../centroids/reference0.mtz >&! rewater.log

  awk '{print substr($0,1,80)}' rewatered.pdb >! reorgme.pdb
  reorganize_pdb_runme.com reorgme.pdb refpdb=fulllength_noalt.pdb outfile=reorged.pdb \
    phenix_bumpcheck=0 debug=1 autorerun=0 >&! reorg_rewatered.log

  awk '{print substr($0,1,80)}' reorged.pdb >! refme_noalt.pdb
  echo "  geometry minimization (details: super_refine1/confsel_min.log)..."
  phenix.geometry_minimization refme_noalt.pdb prefix=confsel_min $ligcifs \
    automatic_linking.link_none=True nonbonded_weight=500 >&! confsel_min.log
  # phenix.geometry_minimization outputs prefix.pdb (not prefix_001.pdb)
  if (! -e confsel_min.pdb) then
    echo "ERROR: phenix.geometry_minimization failed — details: confsel_min.log"
    goto exit
  endif
else
  echo "  confsel + geometry minimization: already done, skipping"
endif

# Fix cis-peptides using generate_omega_fix_runme.com (produces omega_fix.eff)
generate_omega_fix_runme.com confsel_min.pdb >&! omega_fix_gen.log
echo "  phenix.refine with omega fix (details: super_refine1/phenix_omegafix1.log)..."
phenix.refine confsel_min.pdb ../refme.mtz prefix=omegafix1 opts.eff $ligcifs \
  omega_fix.eff >&! phenix_omegafix1.log
if (! -e omegafix1_001.pdb) then
  echo "ERROR: phenix.refine omegafix failed — details: phenix_omegafix1.log"
  goto exit
endif
echo "  super_refine1 done"

ln -sf omegafix1_001.pdb thisone.pdb

cd ..
sec9done:


#==============================================================================
# SECTION 10 — Initial Amber MD run (amber1)
#==============================================================================
if (-e amber1/Prod.rst7) then
  echo ""
  echo "=== Section 10: already done, skipping ==="
  goto sec10done
endif
echo ""
echo "=== Section 10: Initial Amber MD run (amber1) ==="
# IMPORTANT: skip Cpu and Min stages — they reliably produce NaN coordinates
# ("black holes") on 6c2r. Start directly from Cool.
mkdir -p amber1
if (-e compute_settings.sourceme) ln -sf ../compute_settings.sourceme amber1/
cd amber1

foreach lig ( $ligands $salt )
  foreach ext ( cif pdb mol2 frcmod )
    cp ../ligands/${lig}.${ext} .
  end
end
#set ligs = `ls -1 *.mol2 | awk -F "." '{print $1}'`
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
foreach lig ( $ligands $salt )
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
  outfile=restraints_for_${itr}.pdb debug=$debug >&! c2r_${itr}.log
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

echo "  adding salt ions and padding waters..."
add_salt_runme.com refined.pdb conc=$salt_conc \
  RIP=4 RIW=3 charge=$charge anion=$anion cation=$cation >&! add_salt.log
if (! -e salty.pdb) then
  echo "ERROR: add_salt_runme.com did not produce salty.pdb — details: add_salt.log"
  goto exit
endif

echo "  reorganizing supercell..."
reorganize_pdb_runme.com salty.pdb ignore_zero=0 refpdb=refined.pdb \
  outfile=reorganized.pdb phenix_bumpcheck=0 declash=1 >&! reorganize_final.log

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
set padwater = `echo $waterslots $gotwater | awk '{v=$1*1.5;t=$2-$3;print 0+sprintf("%.2g",(v>t?v:t))}'`
echo "padwater = $padwater"

cp restraints_for_${itr}.pdb initial_restraints.pdb

# Run MD stages: skip Cpu and Min (cause NaN / black holes on this system)
echo "  running leap2amber MD stages: Cool → Heat → Equi → EquiMin → Prod (details: leap2amber_${itr}.log)..."
leap2amber.com amberme.pdb stages=Cool,Heat,Equi,EquiMin,Prod \
  protons=protonation.txt watertype=fb3mod flexwater=0 \
  refpoints=initial_restraints.pdb restraint_mult=1 \
  pdbscale=0.01 gamma_ln=1.0 barostat=1 \
  leapfile=tleap_stub.in padwater=$padwater \
  cool_ns=0.001 heat_ns=0.5 equi_ns=0.5 prod_ns=0.5 \
  cool_slowdown=5 heat_slowdown=1 equi_slowdown=1 \
  restrain_omega=0 omega_weight=0 chiral_weight=0 \
  debug=0 >&! leap2amber_${itr}.log
if (! -e Prod.rst7) then
  echo "ERROR: leap2amber.com did not produce Prod.rst7 — details: leap2amber_${itr}.log"
  goto exit
endif
# Produces: Prod.rst7  xtal.prmtop  padded.parm7  orignames.pdb
#           chir_omega0.rst  Bfac_0.pdb  restraints_for_0.pdb

# template_dir is used by optimize_weights_runme.com to find amberme.pdb,
# tleap_stub.in, mol2/frcmod, protonation.txt, reference.mtz, orignames.pdb.
# Point it at amber1 (which has all these files); optimize_weights adds more.
if (! -e reference.mtz) ln -sf ../centroids/reference0.mtz reference.mtz
if (! -e ../template_dir) ln -sf amber1 ../template_dir

cd ..
sec10done:


#==============================================================================
# CONVERGENCE POINT
# Systems diverge from PDB download through leap2amber.com (different ligands,
# crystal forms, MTZ columns). After leap2amber.com, all systems use the same
# optimize_weights_runme.com structure with only these parameter differences:
#
#   System         min_lig_weight   halfrho_neg
#   1aho           0                3.5           (this script)
#   6c2r / EG7     0.1              3.5
#   6c2r / AMPPNP  0.01             3.5
#   2qpx           TBD              TBD
#==============================================================================


#==============================================================================
# SECTION 11 — Weight optimization setup (opt1)
#==============================================================================
if (-e opt1/runme2.log) then
  echo ""
  echo "=== Section 11: already done, skipping ==="
  goto sec11done
endif
echo ""
echo "=== Section 11: Weight optimization (opt1) ==="
set o = 1
mkdir -p opt${o}
if (-e compute_settings.sourceme) ln -sf ../compute_settings.sourceme opt${o}/
cd opt${o}

# Resume: runme1.log present means Stage 1 already ran, so skip the one-time
# setup and Stage 1 and jump straight to Stage 2.
if (-e runme1.log) goto opt1_stage2

# One-time setup: seed opt1 from the last good iteration of the previous run.
# Set these to point to the last good iteration from amber1 (or a previous opt):
set previtr  = 0          # last good iteration number (0 = use Prod directly)
set prevdir  = amber1     # directory containing that iteration
set prevprod = Prod       # file stem of the amber rst7/in/out/nc to continue from

cp ${pdir}/optimize_weights_runme.com .
cp ../${prevdir}/restraints_for_${previtr}.pdb current_restraints.pdb
cp current_restraints.pdb restraints_for_0.pdb
cp ../${prevdir}/${prevprod}.rst7 amber_0.rst7
cp ../${prevdir}/${prevprod}.in  amber_0.in
cp ../${prevdir}/${prevprod}.out amber_0.out
if (-e ../${prevdir}/barometer_${previtr}.out) then
  cp ../${prevdir}/barometer_${previtr}.out barometer_0.out
else
  cp ../${prevdir}/${prevprod}.out barometer_0.out
endif
if (-e ../${prevdir}/leap2amber_${previtr}.log) cp ../${prevdir}/leap2amber_${previtr}.log .
ln -sf ../${prevdir}/${prevprod}.nc amber_0.nc
cp ../centroids/centroids_in_density.pdb all_possible_refpoints.pdb
cp ../${prevdir}/xtal.prmtop .
cp ../${prevdir}/padded.parm7 .
cp ../${prevdir}/orignames.pdb .
if (-e ../${prevdir}/Bfac_${previtr}.pdb) cp ../${prevdir}/Bfac_${previtr}.pdb Bfac.pdb
cp ../${prevdir}/chir_omega0.rst .
cp chir_omega0.rst chir_omega.rst
cp ../xtal_properties.sourceme .

# Stage 1: rough hydration — fill the voids as fast as possible without breaking
# anything.  No weight or B changes.  dehydrate=pressure and pressure_scale=1,auto
# are the defaults and let water content self-tune.  Runs in the background; the
# monitor exits it the moment the mean pressure has gone POSITIVE (voids full).
./optimize_weights_runme.com prod_ns=0.5 max_mult=1 Bfac_maxmod=0 weight_power=1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=1,auto void_scale=1 \
    release_itr=0 repick_itr=0 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=1 weight_negscale=1 randel_itr=0 \
    min_align_weight=0.01 maxitr=20 >&! runme1.log &
set owpid = $!
while ( 1 )
  sleep 120
  ps -p $owpid >& /dev/null
  if ( $status ) break                          # optimizer already stopped (hit cap)
  if ( ! -e pressure_vs_itr.txt ) continue
  # col4 = mean pressure; voids are full once it has gone positive
  set pos = `awk '$4+0 > 0 {f=1} END{print f+0}' pressure_vs_itr.txt`
  if ( "$pos" == "1" ) then
    echo "opt1 stage 1: pressure has gone positive (voids filled) — touching ./exit"
    touch exit
    break
  endif
end
wait
set lastitr = -1
if (-e fofc_Rplot.txt) set lastitr = `tail -n 1 fofc_Rplot.txt | awk '{print $1+0}'`
if ( $lastitr < 1 ) then
  echo "ERROR: optimize_weights Stage 1 made no progress past the seed — see runme1.log"
  goto exit
endif

opt1_stage2:
# Stage 2: confirm the pressure is STABLE before any weight/B optimization.  Same
# rough-hydration regime (no weight changes, no repick) — just keep running until
# the mean pressure stops trending.
./optimize_weights_runme.com prod_ns=0.5 max_mult=1 Bfac_maxmod=0 weight_power=1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=1,auto void_scale=1 \
    release_itr=0 repick_itr=0 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=1 weight_negscale=1 randel_itr=0 \
    min_align_weight=0.01 maxitr=20 >&! runme2.log &
set owpid = $!
set pwin = 10
while ( 1 )
  sleep 120
  ps -p $owpid >& /dev/null
  if ( $status ) break
  if ( ! -e pressure_vs_itr.txt ) continue
  set nP = `wc -l < pressure_vs_itr.txt`
  if ( $nP < $pwin ) continue
  # flat when the least-squares slope of col4 over the last $pwin iters is < half a
  # standard error (pscore < 0.5, the same test optimize_weights uses internally)
  set flat = `tail -n $pwin pressure_vs_itr.txt | awk '{n++;sx+=NR;sy+=$4;sxx+=NR*NR;sxy+=NR*$4;syy+=$4*$4} END{d=(n*sxx-sx*sx);sl=(d!=0)?(n*sxy-sx*sy)/d:0;m=sy/n;sd=sqrt(syy/n-m*m);ps=(sd>0)?sqrt(sl*sl)*sqrt(n)/sd:0;print (ps<0.5)?1:0}'`
  if ( "$flat" == "1" ) then
    echo "opt1 stage 2: pressure stable over last $pwin iters — touching ./exit"
    touch exit
    break
  endif
end
wait
set lastitr = -1
if (-e fofc_Rplot.txt) set lastitr = `tail -n 1 fofc_Rplot.txt | awk '{print $1+0}'`
if ( $lastitr < 1 ) then
  echo "ERROR: optimize_weights Stage 2 made no progress past the seed — see runme2.log"
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
sec11done:


#==============================================================================
# SECTION 12 — Continue optimization (opt2 template)
# Repeat this block for opt3, opt4, ... adjusting parameters as needed.
#==============================================================================
if (-e opt2/fofc_Rplot.txt) then
  echo ""
  echo "=== Section 12: already done, skipping ==="
  goto sec12done
endif
echo ""
echo "=== Section 12: Optimization stage 2 (opt2) ==="
set prevdir = opt1
if (! -e ${prevdir}/fofc_Rplot.txt) then
  echo "ERROR: ${prevdir}/fofc_Rplot.txt not found — opt1 did not complete"
  goto exit
endif
# continue from the LAST iteration (most-evolved state), not the lowest-Rfree
# one - Rfree is noisy and an early frame can win by chance.  Override $previtr
# by hand only if there is a real reason to backtrack (e.g. the last iter blew up).
set previtr = `tail -n 1 ${prevdir}/fofc_Rplot.txt | awk '{print $1+0}'`
# refuse to seed a new stage from a predecessor that never advanced past its seed
# (last itr < 1) - that means the predecessor's amber MD failed, and continuing
# would just cascade the failure (opt5 from opt4 from a dead opt3, etc.).
if ( "$previtr" == "" ) set previtr = 0
if ( $previtr < 1 ) then
  echo "ERROR: ${prevdir} did not advance past its seed (last itr $previtr) - it likely failed."
  echo "       Not starting a new stage from a failed one; see ${prevdir}/runme1.log."
  goto exit
endif
echo "opt2 continuing from ${prevdir} iteration $previtr (last iteration)"
set prevprod = amber_${previtr}

set o = 2
mkdir -p opt${o}
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
if (-e ../${prevdir}/Bfac_${previtr}.pdb) cp ../${prevdir}/Bfac_${previtr}.pdb Bfac.pdb
cp ../${prevdir}/chir_omega0.rst .
cp chir_omega0.rst chir_omega.rst
cp ../xtal_properties.sourceme .

# opt2: weight optimization.  Pressure is stable now, so turn on max_mult=2,
# decaying weights, and repick + release.  void_scale=0 — water is self-tuning by
# pressure, no more void-based addition.  B-factor mods stay OFF until the weights
# have stabilized (opt3).  Exit when the amber RESTRAINT energy
# (amber_energy_vs_itr.txt col16) stops trending — the automatable proxy for the
# sorted_weights.txt log-log shape settling.  Runs in the background; maxitr caps.
./optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=0 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=1,auto void_scale=0 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=0.9 weight_negscale=0.5 randel_itr=0 \
    min_align_weight=0.01 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto maxitr=90 >&! runme1.log &
set owpid = $!
set pwin = 10
while ( 1 )
  sleep 120
  ps -p $owpid >& /dev/null
  if ( $status ) break
  if ( ! -e amber_energy_vs_itr.txt ) continue
  set nE = `wc -l < amber_energy_vs_itr.txt`
  if ( $nE < $pwin ) continue
  # col16 = RESTRAINT energy; flat (pscore < 0.5) = weights stabilized
  set flat = `tail -n $pwin amber_energy_vs_itr.txt | awk '{n++;sx+=NR;sy+=$16;sxx+=NR*NR;sxy+=NR*$16;syy+=$16*$16} END{d=(n*sxx-sx*sx);sl=(d!=0)?(n*sxy-sx*sy)/d:0;m=sy/n;sd=sqrt(syy/n-m*m);ps=(sd>0)?sqrt(sl*sl)*sqrt(n)/sd:0;print (ps<0.5)?1:0}'`
  if ( "$flat" == "1" ) then
    echo "opt2: restraint energy stable over last $pwin iters (weights settled) — touching ./exit"
    touch exit
    break
  endif
end
wait
set lastitr = -1
if (-e fofc_Rplot.txt) set lastitr = `tail -n 1 fofc_Rplot.txt | awk '{print $1+0}'`
if ( $lastitr < 1 ) then
  echo "ERROR: optimize_weights opt2 made no progress past the seed — see opt2/runme1.log"
  goto exit
endif

cd ..
sec12done:

# Opt-stage workflow (each stage seeds from the previous stage's LAST iteration;
# every stage runs in the background and a monitor stops it via ./exit on its own
# convergence signal, not a fixed iteration count):
#   opt1 stage 1: rough hydration     - fill voids fast; exit when pressure > 0
#   opt1 stage 2: pressure stable     - same regime; exit when pressure is flat
#   opt2:         weight optimization - max_mult=2, decaying weights, repick +
#                 release, void_scale=0; exit when the amber RESTRAINT energy is
#                 flat (amber_energy_vs_itr.txt col16 = weights converged)
#   opt3:         B-factor mods ON    - ONLY after the weights converge; keep B
#                 subtle (Bfac_maxmod=2); exit when fofc_R is level
#   opt4+:        continuation of opt3 - more iterations of the coupled B+weight
#                 optimization; exit when fofc_R is level.  Add opt5, opt6 ... the
#                 same way only if fofc_R is still improving.
#
# Goal: a low, LEVEL fofc_R with the MINIMUM restraint energy.  Key points:
#   - B and weight optimization are coupled - neither samples accurate density
#     until both converge - so fofc_R (not any single knob) is the convergence
#     signal for opt3+.  Do NOT enable more than subtle B-factor mods before the
#     weights have converged (watch RESTRAINT energy / the sorted_weights.txt
#     log-log shape stop changing).  Raise Bfac_maxmod toward 50 only once you
#     trust the density.
#   - dehydrate=pressure and pressure_scale=1,auto are the DEFAULTS: let water
#     content self-tune, don't hardcode a pressure_scale number unless it comes
#     from a long known-good run.
#   - allatom_weight is a blanket on every atom and dominates restraint energy:
#     keep it 0 (R floors ~37-38; below ~37 needs a small allatom_weight and costs
#     ~30,000 restraint-energy units per ~0.05 step).
#   - omega/chiral weights are for RECOVERY only (they auto-ramp on when chiral
#     inversions / cis-peptides appear); leave them at their default 0.
# Monitor:
#   cat opt*/fofc_Rplot.txt | sort -k2g | head                       # best R
#   grep -H . opt*/amber_energy_vs_itr.txt | awk '{print $0}' | tail  # RESTRAINT=col16
# Water model (set at leap2amber / Section 10, not here): SPC/E gives the lowest R
# at allatom_weight=0 (~37.7); OPC3 gives the lowest restraint energy.

#==============================================================================
# SECTION 13 — B-factor optimization (opt3)
# Weights are stable now (opt2); turn on B-factor modification.  B and weight
# optimization are coupled, so convergence is judged by the difference-map R.
# Seeds from opt2's last iteration.  Exit when fofc_R is level (low and flat).
#==============================================================================
if (-e opt3/fofc_Rplot.txt) then
  echo ""
  echo "=== Section 13: already done, skipping ==="
  goto sec13done
endif
echo ""
echo "=== Section 13: B-factor optimization (opt3) ==="
set prevdir = opt2
if (! -e ${prevdir}/fofc_Rplot.txt) then
  echo "ERROR: ${prevdir}/fofc_Rplot.txt not found — ${prevdir} did not complete"
  goto exit
endif
# continue from the LAST iteration (most-evolved state), not the lowest-Rfree
# one - Rfree is noisy and an early frame can win by chance.  Override $previtr
# by hand only if there is a real reason to backtrack (e.g. the last iter blew up).
set previtr = `tail -n 1 ${prevdir}/fofc_Rplot.txt | awk '{print $1+0}'`
# refuse to seed a new stage from a predecessor that never advanced past its seed
# (last itr < 1) - that means the predecessor's amber MD failed, and continuing
# would just cascade the failure (opt5 from opt4 from a dead opt3, etc.).
if ( "$previtr" == "" ) set previtr = 0
if ( $previtr < 1 ) then
  echo "ERROR: ${prevdir} did not advance past its seed (last itr $previtr) - it likely failed."
  echo "       Not starting a new stage from a failed one; see ${prevdir}/runme1.log."
  goto exit
endif
echo "opt3 continuing from ${prevdir} iteration $previtr (last iteration)"
set prevprod = amber_${previtr}

set o = 3
mkdir -p opt${o}
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
if (-e ../${prevdir}/Bfac_${previtr}.pdb) cp ../${prevdir}/Bfac_${previtr}.pdb Bfac.pdb
cp ../${prevdir}/chir_omega0.rst .
cp chir_omega0.rst chir_omega.rst
cp ../xtal_properties.sourceme .

# opt3: turn on B-factor modification, now that the weights are stable.  Raise
# Bfac_maxmod toward 50 for serious optimization once you trust the density (start
# gentle here).  B and weight optimization are coupled — neither samples accurate
# density until both converge — so the convergence signal is the difference-map R
# itself: exit when fofc_Rplot.txt col2 is LEVEL.  A low, flat fofc_R is the target
# for the production run.  Runs in the background; maxitr is a safety cap.
./optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=2 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=1,auto void_scale=0 \
    release_itr=30 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=0.9 weight_negscale=0.5 randel_itr=0 \
    min_align_weight=0.01 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto maxitr=90 >&! runme1.log &
set owpid = $!
set pwin = 10
while ( 1 )
  sleep 120
  ps -p $owpid >& /dev/null
  if ( $status ) break                          # optimizer already stopped (hit cap)
  if ( ! -e fofc_Rplot.txt ) continue
  set nR = `wc -l < fofc_Rplot.txt`
  if ( $nR < $pwin ) continue
  # col2 = difference-map R; flat (pscore < 0.5) = B + weight optimization converged
  set flat = `tail -n $pwin fofc_Rplot.txt | awk '{n++;sx+=NR;sy+=$2;sxx+=NR*NR;sxy+=NR*$2;syy+=$2*$2} END{d=(n*sxx-sx*sx);sl=(d!=0)?(n*sxy-sx*sy)/d:0;m=sy/n;sd=sqrt(syy/n-m*m);ps=(sd>0)?sqrt(sl*sl)*sqrt(n)/sd:0;print (ps<0.5)?1:0}'`
  if ( "$flat" == "1" ) then
    echo "opt3: fofc_R level over last $pwin iters (B + weight converged) — touching ./exit"
    touch exit
    break
  endif
end
wait                                             # let opt3 finish its current iteration cleanly
# ran in the background, so check it advanced PAST the seed rather than leaving a
# seed-only fofc_Rplot.txt (a failed amber MD exits at the seed and must not cascade).
set lastitr = -1
if (-e fofc_Rplot.txt) set lastitr = `tail -n 1 fofc_Rplot.txt | awk '{print $1+0}'`
if ( $lastitr < 1 ) then
  echo "ERROR: optimize_weights opt3 (B-factors) made no progress past the seed (last itr $lastitr)."
  echo "       The amber MD likely failed - check opt3/runme1.log (GPU/CUDA errors, ns/day ~0)."
  goto exit
endif

cd ..
sec13done:


#==============================================================================
# SECTION 14 — B-factor optimization continued (opt4, continuation of opt3)
# Same regime as opt3 (subtle B-factor mods; weights already converged) — just
# more iterations to let the coupled B + weight optimization settle.  Seeds from
# opt3's last iteration.  Exit when fofc_R is level.  Add more continuation
# sections (opt5, opt6, ...) the same way only if fofc_R is still improving.
#==============================================================================
if (-e opt4/fofc_Rplot.txt) then
  echo ""
  echo "=== Section 14: already done, skipping ==="
  goto sec14done
endif
echo ""
echo "=== Section 14: B-factor optimization continued (opt4) ==="
set prevdir = opt3
if (! -e ${prevdir}/fofc_Rplot.txt) then
  echo "ERROR: ${prevdir}/fofc_Rplot.txt not found — ${prevdir} did not complete"
  goto exit
endif
# continue from the LAST iteration (most-evolved state), not the lowest-Rfree
# one - Rfree is noisy and an early frame can win by chance.  Override $previtr
# by hand only if there is a real reason to backtrack (e.g. the last iter blew up).
set previtr = `tail -n 1 ${prevdir}/fofc_Rplot.txt | awk '{print $1+0}'`
# refuse to seed a new stage from a predecessor that never advanced past its seed
# (last itr < 1) - that means the predecessor's amber MD failed, and continuing
# would just cascade the failure (opt5 from opt4 from a dead opt3, etc.).
if ( "$previtr" == "" ) set previtr = 0
if ( $previtr < 1 ) then
  echo "ERROR: ${prevdir} did not advance past its seed (last itr $previtr) - it likely failed."
  echo "       Not starting a new stage from a failed one; see ${prevdir}/runme1.log."
  goto exit
endif
echo "opt4 continuing from ${prevdir} iteration $previtr (last iteration)"
set prevprod = amber_${previtr}

set o = 4
mkdir -p opt${o}
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
if (-e ../${prevdir}/Bfac_${previtr}.pdb) cp ../${prevdir}/Bfac_${previtr}.pdb Bfac.pdb
cp ../${prevdir}/chir_omega0.rst .
cp chir_omega0.rst chir_omega.rst
cp ../xtal_properties.sourceme .

# Continuation of opt3: same subtle B-factor optimization, more iterations.  Runs
# in the background; exit when fofc_R (fofc_Rplot.txt col2) is level.
./optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=2 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=1,auto void_scale=0 \
    release_itr=30 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scale=0.9 weight_negscale=0.5 randel_itr=0 \
    min_align_weight=0.01 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto maxitr=90 >&! runme1.log &
set owpid = $!
set pwin = 10
while ( 1 )
  sleep 120
  ps -p $owpid >& /dev/null
  if ( $status ) break
  if ( ! -e fofc_Rplot.txt ) continue
  set nR = `wc -l < fofc_Rplot.txt`
  if ( $nR < $pwin ) continue
  # col2 = difference-map R; flat (pscore < 0.5) = B + weight optimization converged
  set flat = `tail -n $pwin fofc_Rplot.txt | awk '{n++;sx+=NR;sy+=$2;sxx+=NR*NR;sxy+=NR*$2;syy+=$2*$2} END{d=(n*sxx-sx*sx);sl=(d!=0)?(n*sxy-sx*sy)/d:0;m=sy/n;sd=sqrt(syy/n-m*m);ps=(sd>0)?sqrt(sl*sl)*sqrt(n)/sd:0;print (ps<0.5)?1:0}'`
  if ( "$flat" == "1" ) then
    echo "opt4: fofc_R level over last $pwin iters — touching ./exit"
    touch exit
    break
  endif
end
wait
set lastitr = -1
if (-e fofc_Rplot.txt) set lastitr = `tail -n 1 fofc_Rplot.txt | awk '{print $1+0}'`
if ( $lastitr < 1 ) then
  echo "ERROR: optimize_weights opt4 made no progress past the seed (last itr $lastitr)."
  echo "       The amber MD likely failed - check opt4/runme1.log (GPU/CUDA errors, ns/day ~0)."
  goto exit
endif

cd ..
sec14done:

exit:
