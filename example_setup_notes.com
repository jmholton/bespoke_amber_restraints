#! /bin/tcsh -f
#
#   example demonstration of applying and optimzing bespoke restraints
#
#   requirements: CCP4 Suite, Phenix Suite, Amber 18 or higher, and gnuplot
#

# start from scratch
set pdbid = 6c2r

# crystal-specific stuff
# faster to do just one cell
set super_mult = 2,2,2

set ligands = ""
set salt = ( NH4 SO4 )
set salt_conc = 0.15

# local file configuration stuff
set pwd = `pwd`
set pdir = $pwd
#set pdir = ~/projects/amber/1aho_refine
set t = tempfile
set path = ( $pdir $path )
source /programs/amber22/amber.csh

set srun = "srun"
set debug = 1



# download the model
#phenix.fetch_pdb 1aho action=all
#phenix.cif_as_mtz 1aho-sf.cif 
getcif.com $pdbid
# get full sequence from PDB
grep SEQRES ${pdbid}.pdb |\
sequence.awk |\
awk 'NR==1{$0="> "$0} {print} NF==0{exit}' |\
tee seq.fasta

set modulo = `awk 'NR>1{printf("%s",$1)}' seq.fasta |  wc -c`

set pdbfile = ${pdbid}.pdb
set mtzfile = ${pdbid}.mtz

# extract only interesting columns so phenix knows which ones to use
cad hklin1 $mtzfile hklout refme_small.mtz << EOF
labin file 1 E1=F E2=SIGF E3=FreeR_flag
labou file 1 E1=FP E2=SIGFP E3=FreeR_flag
EOF

# standard nomenclature from now on
set mtzfile = refme_small.mtz
set pdbfile = starthere_asu.pdb




# extract small-cell parameters
set pdbCELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
echo head | mtzdump hklin $mtzfile >! smallmtzdump.txt
set SGnum = `awk '/Space group =/{print $NF+0}' smallmtzdump.txt | tail -n 1`
set smallSG = `awk -v num=$SGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
set smallCELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' smallmtzdump.txt`
set nsymops = `awk -v SG=$smallSG '$4==SG {print $2}' ${CLIBD}/symop.lib | head -1`
set reso = `awk '/Resolution Range :/{getline;getline;print $(NF-2)+0;exit}' smallmtzdump.txt | head -1`

# record here for many subsequent programs so they dont have to keep doing the above
cat << EOF >! xtal_properties.sourceme
set reso = $reso
set smallSG = $smallSG
set smallCELL = ( $smallCELL )
set nsymops = $nsymops
# supercell parameters
set super_mult = $super_mult
# number of residues in one protein
set modulo  = $modulo
# any ligands
set ligands = ( $ligands )
# salt used in crystallization
set salt = ( $salt )
# molar
set salt_conc = $salt_conc
EOF
source xtal_properties.sourceme









# expand mtz data to supercell
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
echo head | mtzdump hklin refme.mtz >! mtzdump.txt
set CELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' mtzdump.txt`

echo "" >! blank.pdb
echo "CELL $CELL\nSPACE 1" | pdbset xyzin blank.pdb xyzout ${t}cell.pdb
egrep "^CRYST1" ${t}cell.pdb >! cell.pdb

rm blank.pdb refme_cell.mtz






# try to generate amber input files for any ligands
mkdir ligands
if( "$ligands" == "" ) then
  set ligands = `filter_pdb.awk -v only=ligand,atoms $pdbfile | awk '/^HETAT|^ATOM/{print substr($0,18,3)}' | sort -u`
endif
cd ligands
foreach lig ( $ligands $salt )
  if(-e ${lig}.cif) then
    phenix.elbow ${lig}.cif --id=${lig} --opt --opt_nproc=10 --amber_force_field_files
  endif
  if(-e ${lig}.mol2) continue

  phenix.elbow --chemical_component $lig --id=${lig} --amber_force_field_files --opt --opt_nproc=10
end

cd ..


# use phenix quantum interface to estimate HIS protonations
# see separate repo for this
cp ~/projects/his_flips/${pdbid}/HIS_votes.txt HIS_votes.txt 
awk '! seen[$1,$NF]{string[$1]=string[$1]" "$NF;++seen[$1,$NF]}\
  END{for(r in string)print r,string[r]}' HIS_votes.txt |\
sort -g |\
awk '{print $2,"A"$1}' |\
tee HIS_settings_asu.txt



# quickly use refmac phenix and coot refinement to build missing atoms
# see sanitize_pdb_notes.com
mkdir build1
cd build1

ln -sf ../starthere_asu.pdb starthere0.pdb
cp starthere0.pdb starthere.pdb
ln -sf ../refme_small.mtz refme.mtz

# see sanitize_pdb_notes.com
#ln -sf phenix_002.pdb thisone.pdb
ln -sf bestgeo.pdb thisone.pdb



# test amber conversion for just one monomer
mkdir ../amber_asu/
cd ../amber_asu

cp ../ligands/*.mol2 .
cp ../ligands/*.frcmod .

ln -sf ../build1/bestgeo.pdb starthere.pdb

cat starthere.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
     {print substr($0,1,16),substr($0,18)}' |\
awk '! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
cat >! amberme.pdb


# these are appropriate for single ASU but also OK to not have this and take tleap defaults in amber
cp ../HIS_settings_asu.txt HIS_settings.txt

# minimalistic tleap input
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

tleap -f tleap_stub.in  | tee tleap.log

grep "unperturbed charge" tleap.log
set charge0 = `awk '/unperturbed charge/{gsub(/[)(]/,"");print int($7);exit}' ../amber_asu/tleap.log`
set charge = `echo $charge0 $nsymops | awk '{print $1*$2}'`
echo "$charge" | tee cell_charge.txt

#cd ..



# make a pdb file that only has named atoms for things that are certainly correct
# all others get to be waters
# this will serve as basis of bespoke restraint reference points
mkdir -p ../garr1/
cd ../garr1/

ln -sf ../build1/minRfree.pdb starthere0.pdb
cp starthere0.pdb starthere.pdb
ln -sf ../refme_small.mtz refme.mtz

cat << EOF >! opts.eff
refinement {
  refine {
    occupancies {
      individual = water
    }
  }
  bulk_solvent_and_scale {
    apply_back_trace=False
  }
  main {
    max_number_of_iterations=100
  }
}
EOF

generate_alignment_reference_runme.com starthere.pdb \
  repulse_nb=100 crush_nb=100 repulse_scale=0.5 >&! garr.log 



#cd - 



mkdir ../centroids/
cd ../centroids

# best map
grep "Final R" ../build1/*.log ../garr*/*.log | justify.awk  | sort -k7g | tee sorted.txt | head
set minRfree = `awk '/_/{gsub(".log:"," ");print $1;exit}' sorted.txt`
ln -sf ${minRfree}.mtz minRfree.mtz

# reference points
ln -sf ../garr1/centroids_final_001.pdb centroids_asu.pdb
#cp ../ground_truth.pdb centroids_asu.pdb

# all atoms
ln -sf ../amber_asu/amberme.pdb fulllength_asu.pdb
#cp ../ground_truth.pdb fulllength_asu.pdb


# refinement data
ln -sf ../refme_small.mtz .

# any needed ligands
cp ../ligands/*.cif .
set ligcifs = `ls -1 *.cif |& awk '/.cif/'`

# how different are they?
# maybe do something here to normalize nomenclature
flip_to_target_runme.com centroids_asu.pdb fulllength_asu.pdb > /dev/null
filter_pdb.awk -v skip=H,water fulllength_asu.pdb flipped.pdb | rmsd | head

# create ref map
set F = `mtzdmp minRfree.mtz | awk 'NF>5 && /2FOFCWT|FWT/ && ! / DELF/ && $(NF-1)=="F"{print $NF;exit}'`
set P = `mtzdmp minRfree.mtz | awk 'NF>5 && /PH2FOFCWT|PHWT/ && ! /DEL/ && $(NF-1)=="P"{print $NF;exit}'`

cad hklin1 minRfree.mtz hklout reference0.mtz << EOF
#labin file 1 E1=2FOFCWT E2=PH2FOFCWT
labin file 1 E1=$F E2=$P
labou file 1 E1=Fref E2=PHIref
EOF
fft hklin reference0.mtz mapout ffted.map << EOF
labin F1=Fref PHI=PHIref
EOF
mapmask mapin ffted.map mapout reference0.map << EOF
xyzlim asu
scale sigma
EOF


# populate supercell
foreach prefix ( fulllength centroids )
  egrep "^CRYST|^SSBON" ${prefix}_asu.pdb >! renamed.pdb
  filter_pdb.awk -v only=protein,atoms ${prefix}_asu.pdb >> renamed.pdb
  convert_pdb.awk -v only=ligand,atoms -v renumber=ordinal,watS -v fixEe=1 ${prefix}_asu.pdb >> renamed.pdb
  convert_pdb.awk -v only=water,atoms -v renumber=ordinal,watS ${prefix}_asu.pdb >> renamed.pdb
  cp renamed.pdb ${prefix}_renamed.pdb 
  awk '! /^ATOM|^HETAT/{print;next} {print substr($0,1,16),substr($0,18)}' ${prefix}_renamed.pdb  |\
  awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
  cat >! ${prefix}_asu_noalt.pdb
end
# check that centroids are subset of fulllength
combine_pdbs_runme.com centroids_asu_noalt.pdb fulllength_asu_noalt.pdb > /dev/null
filter_pdb.awk -v skip=water new.pdb centroids_asu_noalt.pdb | rmsd
# should be all zero and no warnings

expand2supercell_runme.com fulllength_renamed.pdb refme_small.mtz super_mult=$super_mult \
  outprefix=fulllength_super \
  phenix_bumpcheck=0 debug=1 | tee fulllength_expand0.log

cp fulllength_super.pdb fulllength_super0.pdb
cp new_monomer_rot_trans.txt monomer_rot_trans.txt

# 2nd time, to use mono map
expand2supercell_runme.com fulllength_renamed.pdb refme_small.mtz super_mult=$super_mult \
  outprefix=fulllength_super \
  phenix_bumpcheck=0 debug=1 \
  mono_map=monomer_rot_trans.txt | tee fulllength_expand.log

# also use mono map
expand2supercell_runme.com centroids_renamed.pdb refme_small.mtz super_mult=$super_mult \
  refpdb=fulllength_renamed.pdb outprefix=centroids_super \
  phenix_bumpcheck=0 debug=1 \
  mono_map=monomer_rot_trans.txt | tee centroids_expand.log


# should not be very far away
filter_pdb.awk -v skip=H,water centroids_super.pdb fulllength_super.pdb | rmsd | head
filter_pdb.awk -v skip=H,water centroids_super.pdb fulllength_super.pdb | grep CA | rmsd | head

foreach prefix ( fulllength centroids )
  awk '! /^ATOM|^HETAT/{print;next} {print substr($0,1,16),substr($0,18)}' ${prefix}_super.pdb  |\
  awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
  cat >! ${prefix}_noalt.pdb
end
# check that centroids are subset of fulllength
combine_pdbs_runme.com centroids_noalt.pdb fulllength_noalt.pdb 
filter_pdb.awk -v skip=water new.pdb centroids_noalt.pdb | rmsd
# should be all zero and no warnings


filter_pdb.awk -v skip=H centroids_super.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  {occ=substr($0,55,6)+0}\
  occ>0.01{print substr($0,1,80),"           |",$NF}' >! all_possible_centroids.pdb

rholabel_runme.com all_possible_centroids.pdb reference0.mtz mtzlabel=Fref

cat rholabeled.pdb |\
awk '/^CRYST/{print;next} ! /^ATOM|^HETAT/{next}\
  {rho=$NF;pre=substr($0,1,60);post=substr($0,67)}\
  rho>0{printf("%s%6.2f%s\n",pre,rho,post)}' |\
tee centroids_in_density.pdb | grep LIG

# make sure these match
filter_pdb.awk -v skip=water,H centroids_in_density.pdb fulllength_super.pdb | rmsd | head
filter_pdb.awk -v skip=H,water centroids_in_density.pdb fulllength_super.pdb | grep CA | rmsd | head
# do a test?


#cd ..



# perhaps optional - refine the supercell structure in phenix
mkdir -p ../super_refine1
cd ../super_refine1

ln -sf ../refme.mtz .
ln -sf ../centroids/fulllength_super.pdb starthere.pdb

cp ../build1/opts.eff .
cp ../ligands/???.cif .
set ligcifs = `ls -1 *.cif |& awk '/.cif/' `
set opts = `ls -1 opts.eff |& awk '/.eff$/'`

awk '{print substr($0,1,80)}' starthere.pdb >! refme.pdb

# fastest option
#cp refme.pdb thisone.pdb

phenix.refine ../refme.mtz refme.pdb prefix=phenix $opts $ciffiles >&! phenix1.log

# might be ok, but better to remove altlocs below
ln -sf phenix_001.pdb thisone.pdb



# focus on good geometry
cat << EOF >! geo.eff
refinement {
  refine {
    strategy = *individual_sites individual_sites_real_space rigid_body \
               *individual_adp group_adp tls occupancies group_anomalous den
    occupancies {
      remove_selection = All
    }
  }
  pdb_interpretation {
    nonbonded_weight = 512
  }
  target_weights {
    wxc_scale = 0.1
    wxu_scale = 1
  }
}
EOF
set opts = `ls -1 opts.eff geo.eff |& awk '/.eff$/'`

set seed = 1
#cp ../centroids/fulllength_super.pdb multiconf.pdb
ln -sf phenix_001.pdb multiconf.pdb

#set ss = 1e-6
set ss = 0.01 
  awk '{print substr($0,1,80)}' multiconf.pdb |\
  tee jiggleme.pdb |\
  jigglepdb.awk -v seed=$seed -v shift=byB -v shift_scale=$ss \
    -v disulfide_links=1 \
    -v independent_confsel=1 -v debug=1 |\
  tee jiggle_debug.pdb |\
  awk '/^CRYST/{print} ! /^ATOM|^HETAT/{next}\
    {occ=substr($0,55,6)+0;}\
    occ>0{print}' |\
  awk '{pre=substr($0,1,16);mid=substr($0,18,38);B=substr($0,61);\
     print pre,mid " 1.00" B}' |\
  cat >! confsel.pdb
  # xxdiff jiggleme.pdb confsel.pdb
end

# make sure nothing went missing
awk '! /^ATOM|^HETAT/{print;next}\
  {print substr($0,1,16),substr($0,18)}' multiconf.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
cat >! fulllength_noalt.pdb

reorganize_pdb_runme.com confsel.pdb refpdb=fulllength_noalt.pdb outfile=oneconf.pdb phenix_bumpcheck=0 \
  debug=1 autorerun=0 | tee reorg_oneconf.log

reorganize_waters.com oneconf.pdb super_mult=$super_mult smallmtz=../centroids/reference0.mtz | tee rewater.log

awk '{print substr($0,1,80)}' rewatered.pdb >! reorgme.pdb
reorganize_pdb_runme.com reorgme.pdb refpdb=fulllength_noalt.pdb outfile=reorged.pdb phenix_bumpcheck=0 \
  debug=1 autorerun=0 | tee reorg_rewatered.log

awk '{print substr($0,1,80)}' reorged.pdb >! refme_noalt.pdb
phenix.geometry_minimization refme_noalt.pdb prefix=confsel_min $ligcifs \
 automatic_linking.link_none=True \
 nonbonded_weight=500 >&! confsel_min.log

phenix.refine ../refme.mtz refme_noalt.pdb prefix=confsel $opts $ligcifs serial=002 >! confsel_refine.log

#ln -sf confsel_001.pdb thisone.pdb
ln -sf confsel_002.pdb thisone.pdb

set pdbfile = confsel_002.pdb
set n = 0

phenix.omegalyze $pdbfile |\
tee omegalyze${n}.log |\
awk -F ":" 'BEGIN{RTD=45/atan2(1,1)}\
   ! /^SUMMARY|^resid/ {om=$3/RTD;\
   n=substr($0,3,4);c=substr($0,2,1);\
   energy=(sin(om)/0.07)^2+(1+cos(om))^10;\
   print "OMEGA",energy,c,n,$0}' |\
sort -k2gr |\
tee badomega${n}.txt |\
awk 'BEGIN{print "refinement {\n  geometry_restraints.edits {"}\
  {c=$3;n=$4;\
  print "    dihedral {";\
  print "      action = change";\
  print "      atom_selection_1 = \"name CA and resseq",n,"and chain",c,"\"";\
  print "      atom_selection_2 = \"name C  and resseq",n,"and chain",c,"\"";\
  print "      atom_selection_3 = \"name N  and resseq",n+1,"and chain",c,"\"";\
  print "      atom_selection_4 = \"name CA and resseq",n+1,"and chain",c,"\"";\
  print "      angle_ideal = 180.00";\
  print "      sigma = 1";\
  print "    }";}\
END{print "  }\n}"}' >! omega_fix.eff

@ n = ( $n + 1 )
phenix.refine $pdbfile refme.mtz prefix=omegafix${n} $ciffiles \
 omega_fix.eff >&! phenix_omegafix${n}.log

set pdbfile = omegafix${n}_001.pdb

ln -sf omegafix${n}_001.pdb thisone.pdb


cat << EOF >! refmac_opts.txt
make build Y
make hydr Y
make hout Y
weight matrix 0.05
#vdwrest 10
restr tors include link TRANS  name omega value 180 sigma 2 period 0
restr tors include link PTRANS name omega value 180 sigma 2 period 0
EOF
set n = 1

converge_refmac.com ../refme.mtz $ciffiles refme.pdb trials=1 append nosalvage >&! converge${n}.log &

reoccupy.awk refmacout.pdb >! reocced.pdb

converge_refmac.com ../refme.mtz $ciffiles reocced.pdb trials=1 append nosalvage >>& converge${n}.log &








# now set up the amber MD run
mkdir -p ../amber1
cd ../amber1


cp ../ligands/*.mol2 .
cp ../ligands/*.frcmod .
cp ../ligands/*.pdb .
set ligs = `ls -1 *.mol2 | awk -F "." '{print $1}'`
foreach lig ( $ligs )
  cp ../ligands/${lig}.cif .
end

ln -sf ../refme.mtz .

#ln -sf ../super_refine1/confsel_001.pdb starthere.pdb
ln -sf ../super_refine1/confsel_min.pdb starthere.pdb
#ln -sf ../super_refine1/thisone.pdb  starthere.pdb
#ln -sf ../centroids/bestgeo_super.pdb starthere.pdb
#ln -sf ../super_refine1/refmac12.pdb starthere.pdb
cp starthere.pdb refined.pdb


set debug = 1
set itr = 0
set weight0 = 1
set pdbscale = 0.01
echo "generating all_possible_refpoints.pdb with weight0 = $weight0"
set B0 = `echo $weight0 $pdbscale | awk '{print $1/$2}'`
cat ../centroids/centroids_in_density.pdb |\
awk -v B0=$B0 '/^CRYST|^LINK|^SSBO/{print} ! /^ATOM|^HETAT/{next}\
    {pre=substr($0,1,60);post=substr($0,67);rho=$NF;\
    B=B0*rho}\
    B>B0{B=B0}\
    {printf("%s%6.2f%s\n",pre,B,post)}' |\
cat >! refpoints_in_density.pdb


# make sure all ligand atoms have potential restraints
set liglist = `echo $ligands | awk '{gsub(" ",",");print}'`

awk '/^CRYST|^ATOM|^HETAT/' refined.pdb |\
filter_pdb.awk -v only=ligand -v ligands="$liglist" -v skip=H |\
 reformatpdb.awk -v BFAC=$B0 >! lig_restraints.pdb

combine_pdbs_runme.com lig_restraints.pdb refpoints_in_density.pdb refined.pdb \
  outfile=all_possible_refpoints.pdb



# first try at restraints
centroids_nearby_runme.com refined.pdb reffile=all_possible_refpoints.pdb \
  softener=2 weight=Bfac maxdist=1  \
  hohscale=1 \
  outfile=density_restraints.pdb debug=$debug | tee c2r_${itr}.log

# make sure ligands have restraints
set liglist = `echo $ligands | awk '{gsub(" ",",");print}'`

cat amberme.pdb |\
filter_pdb.awk -v only=ligand -v ligands="$liglist" -v skip=H |\
 reformatpdb.awk -v BFAC=$B0 >! lig_restraints.pdb

combine_pdbs_runme.com density_restraints.pdb lig_restraints.pdb all_possible_refpoints.pdb \
  saveXYZ=1 outfile=initial_restraints.pdb
cp initial_restraints.pdb restraints_for_0.pdb
cp restraints_for_0.pdb current_restraints.pdb


if( 0 ) the
# alternative: restrain to starting point
set B0 = `echo $weight0 $pdbscale | awk '{print $1/$2}'`
set liglist = `echo $ligands | awk '{gsub(" ",",");print}'`
awk '/^CRYST|^LINK|^SSBO|^ATOM|^HETAT/' refined.pdb |\
filter_pdb.awk -v only=protein,ligand -v ligands="$liglist" -v skip=H |\
 reformatpdb.awk -v BFAC=$B0 >! uniform_restraints.pdb
cp uniform_restraints.pdb restraints_for_0.pdb
endif

rmsd current_restraints.pdb refined.pdb | head | grep MAXD

set avgB = `awk '/^ATOM|^HETAT/{print substr($0,61,6)}' all_possible_refpoints.pdb | avg.awk`
#set weight_scaledown = `echo $avgB | awk '{print 5/($1/100)}'`

set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' refined.pdb`

# minimalistic tleap input
cat << EOF >! tleap_stub.in
source leaprc.protein.ff19SB
source leaprc.water.opc
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
set default FlexibleWater on
set default nocenter on
saveAmberParm x xtal.prmtop start.crd
quit
EOF

# guess charge from single-asu run
set charge0 = `awk '/unperturbed charge/{gsub(/[)(]/,"");print int($7);exit}' ../amber_asu/tleap.log`
set ncells = `echo $super_mult | awk -F "[ ,x]" '{print $1*$2*$3}'`
set nmonomers = `echo $nsymops $ncells | awk '{print $1*$2}'`
set pcharge = `echo $charge0 $nmonomers | awk '{print $1*$2}'`
echo "predicted charge from asu run: $pcharge"
# set charge = $pcharge

# bring in qantum-determined HIS protonation states for the monomer
cp ../HIS_settings_asu.txt HIS_settings_asu.txt 
awk '{print $1,substr($2,2)}' HIS_settings_asu.txt |\
cat >! unshifted_HIS_settings.txt
set firstHIS = `grep ATOM refined.pdb | grep HIS | head -n 1 | awk '{print substr($0,23,4)}'`
if(! $?modulo ) then
 set modulo = `grep OXT refined.pdb | head -n 1 | awk '{print substr($0,23,4)}'`
endif

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

# cp -p ../../6c2r_37C_2x/amber12/HIS_settings.txt .

# run tleap now to get charge?
if ( 0 ) then
  cat HIS_settings.txt refined.pdb |\
   convert_pdb.awk -v fixEe=1 |\
  cat >! protonated.pdb

  egrep "HIS|HID|HIE|HIP" protonated.pdb | grep " CA "

  cat protonated.pdb | convert_pdb.awk -v output=amber -v skip=H,EP >! tleapme.pdb 

  tleap -f tleap_stub.in  | tee tleap.log

  awk '/unperturbed charge/{gsub(/[)(]/,"");print;exit}' tleap.log
  set charge = `awk '/unperturbed charge/{gsub(/[)(]/,"");print int($7);exit}' tleap.log`
  echo "charge is $charge "

  echo "energy check:"
  echo "energy e out energy.dat " |\
  cpptraj.OMP -p xtal.prmtop -y start.crd >! cpptraj_energycheck.log
  cat energy.dat
else
  # calculate from single-ASU run
  set charge = $pcharge
endif

addsalt:

foreach ion ( $salt )
#  set icharge = `awk '/^ATOM|^HETAT/{sum+=substr($0,79)} END{print sum}' ${ion}.pdb`
  set icharge = `awk '/TRIPOS.ATOM/,/TRIPOS.BOND/{sum+=$NF} END{printf("%.1g",sum)}' ${ion}.mol2`
  echo "$ion charge $icharge"
  set test = `echo $icharge | awk '{print ($1 > 0)}'`
  if( $test ) then
    echo "$ion is the salt cation"
    set cation = "$ion"
    ln -sf ${ion}.pdb cation.pdb
  else
    echo "$ion is the salt anion"
    set anion = "$ion"
    ln -sf ${ion}.pdb anion.pdb
  endif
end

# neutralize charge with specified salt species
add_salt_runme.com refined.pdb conc=$salt_conc \
  RIP=4 RIW=3 \
  charge=$charge anion=$anion cation=$cation | tee add_salt.log
# also adds some water
# creates: salty.pdb


# put fluffy stuff at end of file
reorganize_pdb_runme.com salty.pdb ignore_zero=0 refpdb=refined.pdb \
   outfile=reorganized.pdb phenix_bumpcheck=0 declash=1 | tee reorganize_final.log

# preserve SSBONDS?

# should not need to re-do centroids_nearby because top waters should keep same names
# make sure these match
filter_pdb.awk -v skip=water,H current_restraints.pdb reorganized.pdb | rmsd | head
# do a restraint bomb test?

# sometimes tleap hates the hydrogens
filter_pdb.awk -v skip=water reorganized.pdb | egrep -v "^END" >! amberme.pdb
filter_pdb.awk -v skip=H -v only=water,atoms reorganized.pdb >> amberme.pdb

# make this the protonation record
cp HIS_settings.txt protonation.txt

# use refmac for quick B factor optimization
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

converge_refmac.com refme.mtz amberme.pdb trials=3 $ligcifs append nosalvage >&! refmac_${itr}.log &


# estimate how much water could possibly fit
egrep -v HOH amberme.pdb >! dry.pdb
echo | rwcontents xyzin dry.pdb >! ${t}rwcontents.log
grep "% of cell without atoms" ${t}rwcontents.log
set waterslots = `awk '/Cell volume/{V=$NF;n=55*6.022e23/1e27*V} /% of cell without atoms/{print int(n*$NF/100),int(n)}' ${t}rwcontents.log`
echo "looks like room for $waterslots[1] waters, $waterslots[2] in cell"
set gotwater = `grep "O   HOH" amberme.pdb | wc -l`
#set padwater = `echo $waterslots  | awk '{print 0+sprintf("%.1g",($2)*1.5)}'`
set padwater = `echo $waterslots $gotwater | awk '{print 0+sprintf("%.2g",($1-$3)*1.5)}'`
#@ padwater = ( 100000 - $gotwater )

# stages=Cpu,Min,Cool,Heat,Equi,EquiMin,Prod
#cp restraints_for_${itr}.pdb initial_restraints.pdb
cp density_restraints.pdb initial_restraints.pdb
cp density_restraints.pdb restraints_for_0.pdb

leap2amber.com amberme.pdb stages=Cool,Heat,Equi,EquiMin,Prod \
  protons=protonation.txt watertype=fb3mod flexwater=0 \
  refpoints=initial_restraints.pdb restraint_mult=1 \
  pdbscale=0.01 \
  gamma_ln=1.0 barostat=1 \
  leapfile=tleap_stub.in padwater=$padwater \
  cool_ns=0.001 heat_ns=0.5 equi_ns=0.5 prod_ns=0.5 \
  cool_slowdown=1 heat_slowdown=1 equi_slowdown=1 \
  restrain_omega=0 omega_weight=0 chiral_weight=0 \
  debug=1 >&! leap2amber_${itr}.log

# reset: rm -f `ls -1rt | awk '/xtal_properties.sourceme/,""'`


# wait for refmac job to finish
wait 
set lastB = `tail refmacout.pdb | awk '/^ATOM|^HETAT/{print substr($0,61,6)}' | sort -gr | head -n 1`
if( "$lastB" == "") then
  echo "WARNING: could not get last B factors from refmacout.pdb"
  set lastB = 999
endif
grep HOH orignames.pdb >! water.pdb
combine_pdbs_runme.com B=$lastB water.pdb water.pdb outfile=Bwater.pdb > /dev/null
combine_pdbs_runme.com refmacout.pdb Bwater.pdb printref=1 orignames.pdb outfile=Bfac.pdb > /dev/null
cp Bfac.pdb Bfac_${itr}.pdb


# essential files for leap2amber setup
if ( 0 ) then
mkdir ../amber2
cd ../amber2

set prevdir = amber1
foreach file ( amberme.pdb tleap_stub.in protonation.txt restraints_for_0.pdb )
 cp ../${prevdir}/$file .
end
cp ../${prevdir}/*.mol2 ../${prevdir}/*.frcmod .

endif

set n = 1

# hydrate, and wait for pressure to settle
cp ${pdir}/optimize_weights_runme.com optimize_weights_runme${n}.com
./optimize_weights_runme${n}.com prod_ns=0.5 \
    adjust_itr=0 max_mult=1 weight_power=1 \
    Badjust_itr=0 Bfac_maxmod=0 \
    teleport_waters=0 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=0 repick_itr=0 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0.1 cutoff_weight=0.09 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0.01 equi_dt=0.001 equi_gamma=10 settle_slowdown=10 \
    netfrc=0 \
    min_align_weight=0.1 align_target=centroids align_nstlim=0 \
    >&! runme${n}.log &

# reset: rm -f `ls -1rt | awk '/avg_0.mtz/,""'`


# wait for pressure to peek over zero

@ n = ( $n + 1 )
# start teleporting waters
cp ${pdir}/optimize_weights_runme.com optimize_weights_runme${n}.com
./optimize_weights_runme${n}.com prod_ns=0.5 \
    adjust_itr=0 max_mult=1 weight_power=1 \
    Badjust_itr=0 Bfac_maxmod=0 \
    teleport_waters=100 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=0 repick_itr=0 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0.1 cutoff_weight=0.09 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    netfrc=0 \
    min_align_weight=0.1 align_target=centroids align_nstlim=0 \
    refmac_itr=1 >&! runme${n}.log &


# wait for pressure_avglast to become large
set stable = `awk '/avglast/{print $NF}' runme${n}.log | tail -n 1 | awk '{print ( $1 > 20 )}'`

# get starting pressure scale for future runs
grep "pressure scale this time" ../${prevdir}/runme?.log | tee pressure_scale.log
set pressure_scale = `tail -n 10 pressure_scale.log | tac | awk '{print NR,$NF}' | linfit.awk | awk '{print $2+0}'`
echo "default pressure scale of $pressure_scale from now on"

cd ..
ln -sf amber1 template_dir

set previtr = 112
set prevdir = template_dir
set prevprod = amber_${previtr}


set o = 1
set n = 0
mkdir ../opt${o}
cd ../opt${o}

  cp ../${prevdir}/restraints_for_${previtr}.pdb current_restraints.pdb
  cp current_restraints.pdb restraints_for_0.pdb
  cp ../${prevdir}/${prevprod}.rst7 amber_0.rst7
  cp ../${prevdir}/${prevprod}.in amber_0.in
  cp ../${prevdir}/${prevprod}.out amber_0.out
  cp ../${prevdir}/barometer_${previtr}.out barometer_0.out
  cp ../${prevdir}/leap2amber_0.log .
  ln -sf ../${prevdir}/${prevprod}.nc amber_0.nc
  cp ../centroids/centroids_in_density.pdb all_possible_refpoints.pdb
  cp ../${prevdir}/xtal.prmtop .
  cp ../${prevdir}/padded.parm7 .
  cp ../${prevdir}/orignames.pdb .
  cp ../${prevdir}/Bfac_${previtr}.pdb Bfac.pdb
  cp ../xtal_properties.sourceme .


  cp ../${prevdir}/refme.pdb .


# rough restraint opt
@ n = ( $n + 1 )
cp ${pdir}/optimize_weights_runme${n}.com .
./optimize_weights_runme${n}.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=${pressure_scale},auto void_scale=0 \
    release_itr=0 repick_itr=0 repick_maxdist=2 repick_maxweight=0.1 \
    min_lig_weight=0.1 cutoff_weight=0.09 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0.1 align_target=centroids align_nstlim=0 \
    halfrho_neg=3.5 halfrho_pos=auto |& tee runme${n}.log &

# wait for... ?

@ n = ( $n + 1 )
# start repicking
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.09 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0.1 align_target=centroids align_nstlim=0 \
    halfrho_neg=2 halfrho_pos=auto |& tee runme${n}.log &

@ n = ( $n + 1 )
# start down-weighting - will need substages or alignment weight
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.09 allatom_weight=0 \
    weight_scaledown=0.9 weight_negscaledown=0.5 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0.1 align_target=centroids align_nstlim=250000 \
    halfrho_neg=2.5 halfrho_pos=auto |& tee runme${n}.log &


@ n = ( $n + 1 )
# start down-weighting - will need substages or alignment weight
optimize_weights_runme.com prod_ns=0.5 \
    max_mult=2 \
    Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=100 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.09 allatom_weight=0 \
    weight_scaledown=0.9 weight_negscaledown=0.5 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0.1 align_target=centroids align_nstlim=0 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &


@ n = ( $n + 1 )
# try reducing align weight
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.95 weight_negscaledown=0.9 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0.3 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &


@ n = ( $n + 1 )
# zero align weight - does not work!
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.99 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &


# backtrack after something goes wrong
set lastitr = `tail -n 1 fofc_Rplot.txt | awk '{print $1;exit}'`
set gooditr = `tail -n 100 fofc_Rplot.txt | awk '{print $1,$2-f;f+=0.01}' | tee deleteme | sort -k2g | awk '{print $1;exit}'`
cp -p restraints_for_${gooditr}.pdb current_restraints.pdb
cp -p restraints_for_${gooditr}.pdb restraints_for_${lastitr}.pdb 
cp amber_${gooditr}.rst7 amber_${lastitr}.rst7
cp amber_${gooditr}.in amber_${lastitr}.in
cp amber_${gooditr}.out amber_${lastitr}.out
cp amber_${gooditr}.nc amber_${lastitr}.nc
cp barometer_${gooditr}.out barometer_${lastitr}.out 



@ n = ( $n + 1 )
# soften down-weighting - use substages
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.99 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=2000 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &



@ n = ( $n + 1 )
# longer sub-stages
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.99 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=5000 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &


@ n = ( $n + 1 )
# even longer substages - starting to get unstable
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.99 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=10000 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &



@ n = ( $n + 1 )
# repick less frequently - turn on randel
optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=10+5 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.99 weight_negscaledown=1 \
    randel_itr=10 randel_fraction=0.05 randel_trigger=1 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=2000 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &


@ n = ( $n + 1 )
# repick less frequently - turn on randel - ignore restricted worsts
../optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=5 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.99 weight_negscaledown=1 \
    randel_itr=100 randel_fraction=0.05 randel_trigger=1 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=2000 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &



@ n = ( $n + 1 )
# R factor getting high turn off randel and other down-weights
../optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=10 repick_maxdist=2 repick_maxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=2000 \
    halfrho_neg=auto halfrho_pos=auto fftB=10 |& tee runme${n}.log &


@ n = ( $n + 1 )
# increase repick distance
../optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=3 repick_maSxweight=0.11 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=2000 \
    halfrho_neg=auto halfrho_pos=auto fftB=10 |& tee runme${n}.log &

# increase repick_maxweight

# try forgetting refpoints that dip below cutoff
cp optimize_weights_runme.com optimize_weights_runme.com
@ n = ( $n + 1 )
# increase repick distance
../optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1.1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=3 repick_maxweight=2 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=2000 \
    halfrho_neg=auto halfrho_pos=auto cutoff_forget=1 |& tee runme${n}.log &



@ n = ( $n + 1 )
# reduce repick weight/dist
../optimize_weights_runme.com prod_ns=0.5 max_mult=1.5 Bfac_maxmod=1 weight_power=1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2.5 repick_maxweight=0.15 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=2000 \
    halfrho_neg=auto halfrho_pos=auto cutoff_forget=1 |& tee runme${n}.log &



@ n = ( $n + 1 )
# longer runs and less frequent releases
../optimize_weights_runme.com prod_ns=5 max_mult=1.5 Bfac_maxmod=1 weight_power=1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=100 repick_itr=5 repick_maxdist=2.5 repick_maxweight=0.15 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto cutoff_forget=1 |& tee runme${n}.log &


@ n = ( $n + 1 )
# longer runs and less frequent releases
../optimize_weights_runme.com prod_ns=5 max_mult=1.5 Bfac_maxmod=1 weight_power=1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=1 repick_itr=1 repick_maxdist=2.5 repick_maxweight=0.15 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=1 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=0 align_target=centroids align_nstlim=250000 \
    halfrho_neg=auto halfrho_pos=auto cutoff_forget=1 |& tee runme${n}.log &








######################################


@ n = ( $n + 1 )
# try to suck-in after disaster
../optimize_weights_runme.com prod_ns=0.5 max_mult=2 Bfac_maxmod=1 weight_power=1 \
    teleport_waters=1 hydrate_itr=1 add_radius=2.1 \
    pressure_avglast=auto pressure_scale=auto void_scale=1 \
    release_itr=0 repick_itr=1 repick_maxdist=6 repick_maxweight=2 \
    min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
    weight_scaledown=0.99 weight_negscaledown=1 \
    randel_itr=0 randel_fraction=0.05 randel_trigger=0 \
    equi_ns=0 \
    min_align_weight=1 align_target=centroids align_nstlim=9999999 \
    halfrho_neg=auto halfrho_pos=auto |& tee runme${n}.log &



cd ..

ln -sf amber1 template_dir

cd template_dir




# for long runs
mkdir old
foreach itr ( `seq 1 5000` )

echo $itr
  mv wraps_${itr}.txt wrap_${itr}.log Bfac_preadjust_${itr}.pdb sorted_mults_${itr}.txt Bfac_update_${itr}.log premod_Bfac_${itr}.pdb Bfac_${itr}.pdb sorted_Bmods_${itr}.txt teleport_${itr}.log c2r_${itr}.log release_worst_${itr}.log prerelease_restraints_for_${itr}.pdb teleported_waters_${itr}.pdb start_${itr}.pdb restraint_update_${itr}.log actual_restraints_${itr}.pdb restraints_for_${itr}.pdb amber_${itr}.in amber_${itr}_i.in amber_${itr}.mdinfo amber_${itr}.out amber_${itr}.nc barometer_${itr}.out amber_${itr}.rst7 amber_${itr}_unwrapped.rst7 new_positions_${itr}.pdb align_${itr}.log omegalyze_${itr}.log chiralyze_${itr}.log avg_${itr}.mtz nc2mtz_${itr}.log restraints_challenge_vs_weight_${itr}.txt hydrate_${itr}.log dehydrate_${itr}.log predry_Bfac_${itr}.pdb preremap_${itr}.pdb remapped_${itr}.pdb remapped_pairs_${itr}.txt old

end





set previtr = 1480
set prevdir = amber1
set prevprod = amber_${previtr}


set o = 1
mkdir ../opt${o}
cd ../opt${o}
  cp ${pdir}/optimize_weights_runme.com optimize_weights_runme.com

  cp ../${prevdir}/restraints_for_${previtr}.pdb current_restraints.pdb
  cp current_restraints.pdb restraints_for_0.pdb
  cp ../${prevdir}/${prevprod}.rst7 amber_0.rst7
  cp ../${prevdir}/${prevprod}.in amber_0.in
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


set n = 1
   optimize_weights_runme.com prod_ns=0.5 max_mult=1 Bfac_maxmod=0 weight_power=1 \
   teleport_waters=0 hydrate_itr=0 add_radius=1.8 pressure_avglast=auto pressure_scale=auto \
   void_scale=0 release_itr=0 repick_itr=0 repick_maxdist=2 repick_maxweight=0.11 \
   min_lig_weight=0 cutoff_weight=0.1 allatom_weight=0 \
   weight_scaledown=1 weight_negscaledown=1 \
   min_align_weight=0 equi_ns=0 \
   align_nstlim=250000 maxitr=2 |& tee runme${n}.log &





# change out water model on-the-fly
cpptraj -p xtal.prmtop << EOF
trajin amber_${itr}.rst7
strip @EPW parmout noEP.parm7
trajout noEP.rst7
EOF
foreach pdb ( orignames.pdb Bfac.pdb )
  grep -v "EPW" $pdb >! noEP.pdb
  mv noEP.pdb $pdb
end
mv noEP.parm7 xtal.prmtop
mv noEP.rst7 amber_${itr}.rst7

update_centroid_positions_runme.com current_restraints.pdb  outfile=ref.crd 

set laStage = amber_$itr

set t = tempfile
awk -v t=$t '{gsub(t,"");print}' ../$prevdir/${t}tleap.in >! tleap.in
cpptraj -p xtal.prmtop -y amber_2.rst7 -x tleapme.pdb 
tleap -f tleap.in >&! tleap.out 

cp tleaped.parm7 xtal.prmtop

set new = amber_3
cp ${laStage}.in ${new}.in

    $pmemd -O -i ${new}.in -o amber.out \
     -p xtal.prmtop \
     -c ${laStage}.rst7 \
     -ref ref.crd \
     -r ${new}.rst7 \
     -x ${new}.nc \
     -inf ${new}.mdinfo &

















