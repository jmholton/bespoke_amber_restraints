#! /bin/tcsh -f
#
#
#
set pdbfile = starthere.pdb


if(-e xtal_properties.sourceme) source xtal_properties.sourceme

if(! $?srun) set srun = "srun"
if(! $?badlinks) set badlinks = 0
if(! $?ligands) set ligands = ""

set ligcifs = `echo $ligands | awk '{for(i=1;i<=NF;++i)print $i ".cif"}'`

zerocyc:
# do zero-cycle refinement in case this is best Rfree in the end
phenix.refine $pdbfile refme.mtz $ligcifs \
  prefix=phenix0 main.number_of_macro_cycles=0 >! phenix0.log &

# regular default phenix refine
phenix.refine $pdbfile refme.mtz $ligcifs \
  prefix=phenix1 >! phenix1.log &

set phenixopts = ""
if( $badlinks ) then
   egrep -v "^LINK" $pdbfile >! nolinks.pdb
   set pdbfile = nolinks.pdb
   set phenixopts = "$phenixopts automatic_linking.link_none=True exclude_from_automatic_linking.selection_1=All exclude_from_automatic_linking.selection_2=All "
endif 

# turn up non-bonds to make it easier to move to Amber
phenix.refine $pdbfile refme.mtz $ligcifs \
  nonbonded_weight=500 wxc_scale=0.1 \
  $phenixopts \
  prefix=phenix_debump >! phenix_debump.log &

phenix.refine $pdbfile refme.mtz $ligcifs \
  nonbonded_weight=500 \
  $phenixopts \
  prefix=phenix_debump2 >! phenix_debump2.log &


# build in any protein that is missing
buildout_pdb_runme.com $pdbfile $ligcifs |& tee buildout1.log
set pdbfile = built.pdb
#set pdbfile = coot.pdb

set geo_opts = ""
if( $badlinks ) then
   egrep -v "^LINK" $pdbfile >! nolinks.pdb
   set pdbfile = nolinks.pdb
   set geo_opts = "$geo_opts link_all=False link_none=True link_ligands=False"
endif 

# make selection mask so that only loose stuff is minimized
filter_pdb.awk -v only=ligands,atoms $pdbfile >! ligands.pdb
rholabel_runme.com $pdbfile phenix0_001.mtz mtzlabel=2FOFCWT
awk '/^ATOM|^HETAT/ && $NF>1.2 && substr($0,12,5)~/  [CN][A ] /' rholabeled.pdb |\
cat - ligands.pdb |\
awk '/^ATOM|^HETAT/{print substr($0,12,5),substr($0,22,1),substr($0,23,8)}' |\
awk '! seen[$0]{++seen[$0];print}' |\
awk 'NF==3{s = "name "$1" and chain "$2" and resseq "$3;\
     printf("(%s) or ",s)}\
     NF==2{s = "name "$1" and resseq "$2;\
     printf("(%s) or ",s)}' |\
awk '{$NF="";print "selection = \"not (" $0 ")\""}' >! density_deselect.eff

phenix.geometry_minimization $pdbfile \
  density_deselect.eff \
  $ligcifs $geo_opts \
  prefix=dd_minimized >&! dd_geomin1.log 
set pdbfile = dd_minimized.pdb

cp dd_minimized.pdb bestgeo.pdb


#
#  might be good enough to stop here
#
if( $?FAST ) then
   wait 
  cp phenix0_001.mtz minRfree.mtz
  cp phenix0_001.pdb minRfree.pdb
  cp dd_minimized.pdb bestgeo.pdb

  exit
endif


echo "" >! refmac_opts.txt
converge_refmac.com $pdbfile ./refme.mtz trials=1 $ligcifs >&! refmac_debump.log 

echo "make link Y" >! refmac_opts.txt
converge_refmac.com $pdbfile ./refme.mtz trials=1 $ligcifs >&! refmac_link.log 

echo "checking for links made by refmac:"
grep LINK refmacout.pdb

echo "adding hydrogens with refmac..."
echo "make build Y" >! refmac_opts.txt
echo "make hout Y" >> refmac_opts.txt
converge_refmac.com refmacout.pdb ./refme.mtz trials=1 $ligcifs append >&! refmac_hout.log 
echo "make hydr Y" >> refmac_opts.txt

cp refmacout.pdb refmacH.pdb

echo "bringing refmac hydrogens to full occupancy"
# make protein with no conformer letter full occupancy
convert_pdb.awk -v only=protein refmacH.pdb |\
awk '/^CRYST/{print} ! /^ATOM|^HETAT/{next}\
  {conf=substr($0,17,1);occ=substr($0,55,6)+0}\
  conf==" " && occ!=1{pre=substr($0,1,54);post=substr($0,61);\
    $0=sprintf("%s%6.2f%s",pre,1,post)} \
  {print}' >! reocc.pdb
convert_pdb.awk -v only=atoms -v skip=protein refmacH.pdb >> reocc.pdb
echo "re-occupy changes:"
rmsd reocc.pdb refmacH.pdb 

# try to get full hydrogen decoration
# reduce does not want to read cif files...
echo "running phenix.reduce"
phenix.reduce reocc.pdb |\
 awk '{gsub(" H1 "," H  ");print}' |\
cat >! reduced.pdb 
echo "refining reduced.pdb in refmac"
converge_refmac.com reduced.pdb ./refme.mtz trials=1 $ligcifs append >&! refmac_reduced.log 
cp refmacout.pdb redref.pdb
echo "re-conforming refmac-refinded reduced"
reconform.com redref.pdb >! reconformed.pdb
echo "changes:"
rmsd redref.pdb reconformed.pdb
echo "refining re-conformed model in refmac"
converge_refmac.com reconformed.pdb ./refme.mtz trials=1 $ligcifs append >&! refmac_reconformed.log 

echo "now back to phenix again"
cp refmacout.pdb phenixme.pdb 
phenix.refine phenixme.pdb refme.mtz prefix=phenixH $ligcifs >! phenixH.log

echo "and back to refmac"
awk '{gsub(" H1 "," H  ");print}' phenixH_001.pdb  >! refmacme.pdb 
convert_pdb.awk -v skip=H refmacme.pdb >! occme.pdb
# tell refmac to keep the input hydrogens
echo "make hydr Y" >! refmac_opts.txt
echo "make hout Y" >> refmac_opts.txt
echo "setting up occupancy refinement in refmac"
refmac_occupancy_setup.com occme.pdb >> refmac_opts.txt 
rm -f refmacout_minRfree.pdb refmacout.pdb
echo "refining with refmac"
converge_refmac.com refmacme.pdb ./refme.mtz $ligcifs append nosalvage keep_zeroocc trials=3 >&! refmac_long.log
echo "starting long refmac refinement..."
converge_refmac.com refmacme.pdb ./refme.mtz $ligcifs append nosalvage keep_zeroocc trials=100 >>& refmac_long.log &


echo "back to phenix with all hydrogens"
cp refmacout.pdb phenixme2.pdb 


# never can remember which direction weight should go
echo "trying varied weighting schemes"
foreach wxc_scale ( 10 2 1 0.5 0.1 0.01 )
foreach wxu_scale ( 10 1 0.1 )

set wxc = `echo $wxc_scale | awk '{gsub("[^0-9]","");print}'`
set wxu = `echo $wxu_scale | awk '{gsub("[^0-9]","");print}'`
 $srun phenix.refine wxc_scale=$wxc_scale wxu_scale=$wxu_scale \
 phenixme2.pdb $ligcifs refme.mtz prefix=wxc${wxc}wxu${wxu} \
  main.number_of_macro_cycles=10 >! wxc${wxc}wxu${wxu}.log &

end
end
echo "waiting for jobs to finish"
wait

if( ! -e phenix0_001.pdb ) then
  echo "oh no. first phenix job failed."
  cat phenix0.log |\
  awk  '/Number of atoms with unknown nonbonded energy type symbols/{\
    getline;print substr($0,match($0,/ATOM|HETAT/),27),"REMOVE"}' |\
  cat - starthere.pdb |\
  awk '{id=substr($0,12,17)}\
     $NF=="REMOVE"{++remove[id];next}\
     remove[id]{next}\
     {print}' >! pruned.pdb
  diff starthere.pdb pruned.pdb
  if( $status ) then
    set pdbfile = pruned.pdb
    goto zerocyc
  else
    set BAD = "cannot refine"
    goto exit
  endif
endif
cp phenix0_001.mtz minRfree.mtz
cp phenix0_001.pdb minRfree.pdb

echo "evaluating geometry"
foreach pdb ( wxc*.pdb *_min*.pdb phenix*_001.pdb )
  set prefix = `basename $pdb .pdb`
  $srun molprobify_runme.com $pdb $ligcifs >&! molprobify_${prefix}.log &
end
wait

grep "Final R" wxc10*.log | justify.awk  | sort -k7g | tee sorted.txt
set heavyXweight = `awk -F "." '{print $1;exit}' sorted.txt`

grep "Final R" *.log | justify.awk  | sort -k7g | tee sorted.txt
set minRfree = `awk -F "." '{print $1;exit}' sorted.txt`
grep wE molprobify_* | sort -k4g | awk '{gsub("^molprobify_|.log:"," ");print}' | tee sorted_geo.txt
set bestgeo = `awk '{print $1;exit}' sorted_geo.txt`

echo "min Rfree was: $minRfree"
echo "best geometry was: $bestgeo"
ln -sf ${minRfree}.mtz minRfree.mtz
ln -sf ${minRfree}.pdb minRfree.pdb
ln -sf ${heavyXweight}.pdb heavyXweight.pdb
ln -sf ${bestgeo}.pdb bestgeo.pdb

# how different are they?
echo "difference between heavy X-ray weight and best geometry:"
flip_to_target_runme.com heavyXweight.pdb bestgeo.pdb > /dev/null
filter_pdb.awk -v skip=H bestgeo.pdb flipped.pdb | rmsd

echo "difference between min Rfree and best geometry"
flip_to_target_runme.com minRfree.pdb bestgeo.pdb > /dev/null
filter_pdb.awk -v skip=H bestgeo.pdb flipped.pdb | rmsd


set pdbfile = bestgeo.pdb
set n = 0
set badomegas = 1

echo "assessing peptide bond omegas"
while ( $badomegas )
  @ n = ( $n + 1 )

( phenix.omegalyze $pdbfile >! omegalyze.log ) >&! omegalyze_stderr.log
cat omegalyze.log |\
awk -F ":" 'BEGIN{RTD=45/atan2(1,1)}\
   ! /^SUMMARY|^resid/ {om=$3/RTD;\
   n=substr($0,3,4);c=substr($0,2,1);\
   energy=(sin(om)/0.07)^2+(1+cos(om))^10;\
   print "OMEGA",energy,c,n,$0}' |\
sort -k2gr |\
tee badomegas.txt |\
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
END{print "  }\n}"}' >! omega_fix${n}.eff

set badomegas = `grep selection badomegas.txt | wc -l`
echo "$badomegas omega angles are bad" 
if( ! $badomegas ) break

phenix.refine $pdbfile refme.mtz prefix=omegafix${n} $ligcifs \
 omega_fix${n}.eff wxc_scale=0.01 nonbonded_weight=512 >! phenix_omegafix${n}.log

set pdbfile = omegafix${n}_001.pdb

end

if("$pdbfile" != "bestgeo.pdb") ln -sf $pdbfile bestgeo.pdb



echo "checking if anything is missing"
awk '/^ATOM|^HETAT/{print substr($0,1,16),substr($0,18)}' bestgeo.pdb |\
awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
cat >! noalt.pdb
filter_pdb.awk -v only=protein,atoms -v skip=H noalt.pdb helix.pdb | grep WARN
# should print nothing
awk '/^ATOM|^HETAT/ && ! /HOH/{print substr($0,1,16),substr($0,18)}' bestgeo.pdb |\
awk '{id=substr($0,12,15);occsum[id]+=substr($0,55,6)}\
  END{for(id in occsum) if(occsum[id]<1) print id,occsum[id]}' 
# should print nothing

echo "best geometry:"
molprobify_runme.com bestgeo.pdb keepgeo |& tee molprobify_bestgeo.log 



