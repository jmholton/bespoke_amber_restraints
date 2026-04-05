#! /bin/tcsh -f
#
# prune off everything that might not be correctly placed protein
#
#
set mtzfile = ../refme_small.mtz
set pdbfile = ""

set initial_reject = 0

set clear_links = 1
set clear_H = 1
set clear_aniso = 1

set repulse_cycles = 3
set repulse_scale = 0.1
set repulse_nb = 100

set crush_cycles = 10
set crush_scale = 10
set crush_nb = 100

# re-label water confs so they dont clash
set reorg_water = 0

# dont allow ligand atoms to be rejected
set ligand_immunity = 1


# minimum density in 2FoFc
set minrho = 0.8
# only look at top outliers in each geo category
set maxgeo = 10
# large moves means atom is loose
set bigmove = 1.5
# define B factor that is too big
set maxBsig = 3
# overall max number of atoms to remove each cycle
set maxbaddies = 100
# minimum before liquefying a fragment
set min_CAs = 2
# only take so many outliers from each category: geo Bfac move
set topbad = 10

# any cif files
set ciffiles = ""

set itr = ""

set tempfile = tempfile
set debug = 0

if(-e xtal_properties.sourceme) source xtal_properties.sourceme

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
      if("$Arg" =~ *.pdb ) then
        set pdbfile = $Arg
        echo "pdbfile = $pdbfile"
        continue
      endif
      if("$Arg" =~ *.cif ) then
        set ciffiles = ( $ciffiles $Arg )
        echo "ciffiles = $ciffiles"
        continue
      endif
      if("$Arg" =~ *.mtz ) then
        set mtzfile = $Arg
        echo "mtzfile = $mtzfile"
        continue
      endif
      if("$Arg" == "$num" ) then
        set itr = $Arg
        echo "user itr = $itr"
        continue
      endif
    endif
    if("$key" == "debug") set debug = "1"
end

# shorthand for temporary file
set t = $tempfile

# make sure other scripts are in the $path
set path = ( $path `dirname $0` )

foreach dependency ( filter_pdb.awk convert_pdb.awk add_waters.com no_new_nonbonds_runme.com combine_pdbs_runme.com rholabel_runme.com molprobify_runme.com median.awk waterconf.awk flip_to_target_runme.com )
   echo -n "using: "
   which $dependency
   if( $status ) then
       set BAD = "need $dependency in "'$'"path"
       goto exit
   endif
end


if( "$itr" == "" ) then
  set itr = `ls -1 | awk -F "_" '/^centroids_/{printf("%03d\n",$2)}' | sort -g | tail -n 1`
  echo "previous last itr = $itr"
endif
if( "$itr" == "" ) set itr = 0
if( "$itr" == "0" ) then

endif

set done = 0

#if(-e ../xtal_properties.sourceme) source ../xtal_properties.sourceme

#ln -sf $pdbfile centroids_start.pdb


if(! -e opts.eff) then
  echo "WARNING: no opts.eff provided. making one"
cat << EOF >! opts.eff
refinement {
  refine {
    strategy = *individual_sites individual_sites_real_space rigid_body \
               *individual_adp group_adp tls *occupancies group_anomalous den
    occupancies {
      individual = water
    }
  }
  pdb_interpretation {
    restraints_library {
      cdl = False
      mcl = False
    }
    flip_symmetric_amino_acids = False
    correct_hydrogens = False
    allow_polymer_cross_special_position = True
    automatic_linking {
      link_none = True
    }
    exclude_from_automatic_linking {
      selection_1 = All
      selection_2 = All
    }
  }
  bulk_solvent_and_scale {
    apply_back_trace=False
  }
  main {
    nqh_flips=False
    max_number_of_iterations=100
  }
}
EOF
endif

if(! -e centroids_start0.pdb) then
  cp $pdbfile centroids_start0.pdb
endif

cp $pdbfile centroids_start.pdb

if ( $clear_H ) then
  echo "removing H atoms"
  filter_pdb.awk -v skip=H centroids_start.pdb >! ${t}new.pdb
  mv ${t}new.pdb centroids_start.pdb
endif

if ( $clear_aniso ) then
  echo "removing anisotropic B factors"
  egrep -v "^ANIS" centroids_start.pdb >! ${t}new.pdb
  mv ${t}new.pdb centroids_start.pdb
endif

set pdbfile = centroids_start.pdb

echo "initial zero-cylce refine..."
# creates phenix_opts_unbump.eff and nnnbp_001.mtz
no_new_nonbonds_runme.com centroids_start.pdb $mtzfile $ciffiles >! nnb.log
if( $status || ! -e phenix_opts_unbump.eff ) then
    tail nnb.log
    set BAD = "failed to create new non-bond deactivation list at repulse step"
    goto exit
endif

# do a probe now
echo "probing density from nnnbp_001.mtz map"
rm -f rholabeled.pdb
filter_pdb.awk -v skip=H -v only=protein centroids_start.pdb |\
 awk '{print substr($0 "       ",1,80)}' >! rhome.pdb
rholabel_runme.com rhome.pdb nnnbp_001.mtz mtzlabel=2FOFCWT outfile=rholabeled.pdb >! rhoprobe.log
if( $status || ! -e rholabeled.pdb ) then
    tail rhoprobe.log 
    set BAD = "rho-label step failed for starting model"
    goto exit
endif

sort -k1.81g rholabeled.pdb  |\
 awk -v minrho=$minrho '$NF<minrho && /^ATOM|^HETAT/{print substr($0,12,15),$NF,"BADRHO"}' |\
 cat >! badrho0.txt
cat badrho0.txt centroids_start.pdb |\
awk '$NF~/^BAD/{++bad[substr($0,1,15)];next}\
  ! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,12,15);xyz=substr($0,31,24)}\
  seen[xyz]{next} {++seen[xyz]}\
  ! bad[id]{print substr($0,1,80)}' >! indensity.pdb
set test = `cat badrho0.txt | wc -l`
echo "$test protein atoms in 2Fo-Fc density < $minrho"

if( $initial_reject ) then
  echo "eliminating them"
  set pdbfile = indensity.pdb
endif


while ( ! $done )
# liquify loop
while ( ! $done )
# add water loop
while ( ! $done )
# prune loop

set itr = `echo $itr | awk '{printf("%03d",$1+1)}'`

if( $reorg_water ) then
    # make waters non-clashy
    echo "de-clashing waters..."
    reorganize_waters.com $pdbfile >! reorg.log
    if( $status ) then
        tail reorg.log
        set BAD = "failed to reoganize waters"
        goto exit
    endif
    filter_pdb.awk -v skip=water rewatered.pdb | egrep -v "^END|^LINK" >! waterconf.pdb
    convert_pdb.awk -v only=water,atoms -v CONF="A" rewatered.pdb |\
    waterconf.awk |\
    egrep "^ATOM|^HETAT" >> waterconf.pdb

    cp waterconf.pdb refme.pdb
else
    if("$pdbfile" != "refme.pdb") cp $pdbfile refme.pdb
endif

if ( $clear_links ) then
  echo "clearing old links"
  egrep -v "^LINK" refme.pdb >! ${t}new.pdb
  mv ${t}new.pdb refme.pdb
endif

set test = `egrep HOH refme.pdb | wc -l`
if( $test == 0 ) then
   echo "model is dry. removing any references to waters from opts.eff"
   grep -v water opts.eff >! new.eff
   mv new.eff opts.eff
endif

echo "zero-cylce refine for new bonds list ..."
# creates phenix_opts_unbump.eff
no_new_nonbonds_runme.com refme.pdb $mtzfile $ciffiles >! nnb.log
if( $status || ! -e phenix_opts_unbump.eff ) then
    tail nnb.log
    set BAD = "failed to create new non-bond deactivation list at repulse step"
    goto exit
endif

echo "refining $itr repulse"
# allow removal of clashes
phenix.refine wxc_scale=$repulse_scale wxu_scale=$repulse_scale nonbonded_weight=$repulse_nb \
  refme.pdb $ciffiles \
  main.number_of_macro_cycles=$repulse_cycles \
  refine.sites.individual="not water" \
  opts.eff phenix_opts_unbump.eff \
  $mtzfile prefix=repulse serial=$itr >! repulse${itr}.log 
if( $status || ! -e repulse_${itr}.pdb ) then
    tail repulse${itr}.log 
    set BAD = "phenix.refine repulsion step failed"
    goto exit
endif

# creates phenix_opts_unbump.eff
rm phenix_opts_unbump.eff
no_new_nonbonds_runme.com repulse_${itr}.pdb $mtzfile $ciffiles >! nnb.log
if( $status || ! -e phenix_opts_unbump.eff ) then
    tail nnb.log
    set BAD = "failed to create new non-bond deactivation list at crush step"
    goto exit
endif

echo "refining $itr crush"
# allow geometry and B factors to go hog wild, go for lowest R and centroids
phenix.refine wxc_scale=$crush_scale wxu_scale=$crush_scale nonbonded_weight=$crush_nb \
 repulse_${itr}.pdb $ciffiles \
  main.number_of_macro_cycles=$crush_cycles \
  opts.eff phenix_opts_unbump.eff \
  $mtzfile prefix=centroids serial=$itr >! centroids${itr}.log 
if( $status || ! -e centroids_${itr}.pdb ) then
    tail centroids${itr}.log 
    set BAD = "phenix.refine crush step failed"
    goto exit
endif

echo "molprobify..."
molprobify_runme.com keepgeo centroids_${itr}.pdb $ciffiles >! molprobify_centroids_${itr}.log
if( $status || ! -e centroids_${itr}_fullgeo.txt ) then
    tail molprobify_centroids_${itr}.log 
    set BAD = "molprobify step failed"
    goto exit
endif
set Rstats = `tail -n 1 molprobify_centroids_${itr}.log | awk '{print $2,$3,$4,$(NF-2)}'`
echo "$itr $Rstats" | tee -a Rstats_vs_itr.txt

set best_itr = `sort -k3g Rstats_vs_itr.txt | awk '{print $1;exit}'`

echo "probing density from centroids_${best_itr} map"
rm -f rholabeled_start.pdb
filter_pdb.awk -v skip=H centroids_start0.pdb | awk '{print substr($0,1,80)}' >! rhome.pdb
rholabel_runme.com rhome.pdb centroids_${best_itr}.mtz mtzlabel=2FOFCWT outfile=rholabeled_start.pdb >! rhoprobe.log
if( $status || ! -e rholabeled_start.pdb ) then
    tail rhoprobe.log 
    set BAD = "rho-label step failed for starting model"
    goto exit
endif

rm -f rholabeled.pdb
filter_pdb.awk -v skip=H centroids_${itr}.pdb | awk '{print substr($0,1,80)}' >! rhome.pdb
rholabel_runme.com rhome.pdb centroids_${best_itr}.mtz mtzlabel=2FOFCWT outfile=rholabeled.pdb >! rhoprobe.log
if( $status || ! -s rholabeled.pdb ) then
    tail rhoprobe.log 
    set BAD = "rho-label step failed"
    goto exit
endif

cat rholabeled_start.pdb rholabeled.pdb |\
awk '! /^ATOM|^HETAT/{next}\
  {id=substr($0,12,17)}\
  startrho[id]==""{startrho[id]=$NF;next}\
  {print startrho[id],$NF,"|" id}' |\
tee rho_diffs.txt |\
awk -v minrho=$minrho '$1>0.9*minrho && $2<$1*0.5' |\
awk -F "|" '{print $2,"REVERT"}' |\
cat - rholabeled_start.pdb |\
awk '$NF=="REVERT"{id=substr($0,1,17);++revert[id];next}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,12,17)}\
  revert[id]{print}' >! revertme.pdb
set revertants = `cat revertme.pdb | wc -l`
echo "$revertants good atoms moved out of density"
if( $revertants ) then
  echo "moving them back"  
  combine_pdbs_runme.com revertme.pdb rholabeled.pdb printref=1 outfile=reverted.pdb >! revert.log

  combine_pdbs_runme.com xor=1 revertme.pdb rholabeled.pdb outfile=moveme.pdb >! xor.log
  cat revertme.pdb |\
  awk '! /^ATOM|^HETAT/{next}\
    {a=substr($0,12,5);c=substr($0,22,1);r=substr($0,23,6);f=substr($0,17,1);\
    s = "name "a" and resseq "r}\
    c!=" "{s=s" and chain "c}\
    f!=" "{s=s" and altid "f}\
    {printf("(%s) or ",s)}' |\
   cat >! holdme.txt
  awk '{$NF=")\"";print "refinement.refine.sites.individual= \"not (", $0}' holdme.txt >! sel.eff
  
  awk '{print substr($0,1,80)}' reverted.pdb >! pinchme.pdb
   #phenix.geometry_minimization sel.eff LIG.cif refme.pdb

  no_new_nonbonds_runme.com pinchme.pdb $ciffiles $mtzfile >! nnb.log
  phenix.refine wxc_scale=$crush_scale wxu_scale=$crush_scale nonbonded_weight=$crush_nb \
    pinchme.pdb $ciffiles \
    main.number_of_macro_cycles=$crush_cycles \
    opts.eff phenix_opts_unbump.eff sel.eff \
    $mtzfile prefix=reverted serial=$itr >&! revert${itr}.log 

  rm -f rholabeled.pdb
    filter_pdb.awk -v skip=H reverted_${itr}.pdb | awk '{print substr($0,1,80)}' >! rhome.pdb
    rholabel_runme.com rhome.pdb centroids_${best_itr}.mtz mtzlabel=2FOFCWT outfile=rholabeled.pdb >! rhoprobe.log
endif


sort -k1.82g rholabeled.pdb  |\
 awk -v minrho=$minrho '$NF<minrho && /^ATOM|^HETAT/{print substr($0,12,15),$NF,"BADRHO"}' |\
cat >! badrho.txt
# idstring  rho BADRHO

awk '/^BOND/' centroids_${itr}_fullgeo.txt |\
awk -F "|" '{gsub("_"," ");split($2,s,"-");\
  id1=substr(s[1],1,5) substr(s[1],7,1) substr(s[1],9,4) substr(s[1],14,1) substr(s[1],16,4);\
  id2=substr(s[2],1,5) substr(s[2],7,1) substr(s[2],9,4) substr(s[2],14,1) substr(s[2],16,4);\
  print id1 "|" id2 "| BONDED" }' |\
cat badrho.txt - |\
awk '{id=substr($0,1,15);id2=substr($0,17,15)}\
  $NF=="BADRHO"{++badrho[id];rho[id]=$(NF-1);next}\
  $NF!="BONDED"{next}\
  ! badrho[id]{++goodneighbor[id2]}\
  ! badrho[id2]{++goodneighbor[id]}\
  END{for(id in badrho){if(badrho[id]+0>0 && ! goodneighbor[id]) print id,goodneighbor[id],"|",rho[id],"BADCONN"}}' |\
cat >! badconn.txt
# idstring 0 | rho BADCONN
# bad density and no bonded neighbors have good density

awk -v maxgeo=$maxgeo '$2>maxgeo{print $0,"|",$2}' centroids_${itr}_worstgeo.txt |\
awk -F "|" '{print $2,"-",$3}' |\
awk -F "-" '{E=$NF;for(i=1;i<NF;++i)print $i,E}' |\
awk '{gsub("_"," ");\
  print substr($0,1,5) substr($0,7,1) substr($0,9,4) substr($0,14,1) substr($0,16,4),"|",$NF,"BADGEO"}' |\
cat >! badgeo.txt

awk -v maxgeo=$maxgeo '$2>maxgeo && $1=="TORSION"{print $0,"|",$2}' centroids_${itr}_worstgeo.txt |\
awk 'BEGIN{DTR=atan2(1,1)/45} cos(DTR*$4)>0' |\
awk -F "|" '{print $2,"-",$3}' |\
awk -F "-" '$1~/ CA / && $4~/ CA /{print $2,$NF;print $3,$NF}' |\
awk '{gsub("_"," ");\
  print substr($0,1,5) substr($0,7,1) substr($0,9,4) substr($0,14,1) substr($0,16,4),"|",$NF,"BADOMEGA"}' |\
cat >! badomega.txt

awk -v maxgeo=$maxgeo '$2>maxgeo && $1=="CHIR" && $4*$5<0{print $0,"|",$2}' centroids_${itr}_worstgeo.txt |\
awk -F "|" '{print $2,"-",$3}' |\
awk -F "-" '{print $1,$NF}' |\
awk '{gsub("_"," ");\
  print substr($0,1,5) substr($0,7,1) substr($0,9,4) substr($0,14,1) substr($0,16,4),"|",$NF,"BADCHIR"}' |\
cat >! badchir.txt

flip_to_target_runme.com centroids_${itr}.pdb centroids_start.pdb >! flip.log
filter_pdb.awk -v skip=H,water centroids_start.pdb flipped.pdb |\
rmsd -v debug=1  |\
  sort -k1.25gr |\
  awk -v bigmove=$bigmove 'substr($0,25)+0>bigmove && /moved/{print $0,"BIGMOVE"}' >! bigmoves.txt

set medmadB = `awk '/^ATOM|^HETAT/{print substr($0,61,6)}' centroids_${itr}.pdb | median.awk`
echo $medmadB $maxBsig |\
 cat - centroids_${itr}.pdb |\
awk 'NR==1{maxBsig=$NF;thresh=$1+maxBsig*$3;next}\
  ! /^ATOM|HETAT/{next}\
  {B=substr($0,61,6)+0;id=substr($0,12,15)}\
  B>thresh{print id,B,"BIGB"}' |\
sort -k1.16gr >! bigB.txt

head -n $topbad badgeo.txt bigmoves.txt bigB.txt |\
cat - badrho.txt badconn.txt |\
 awk '/^==/ || NF==0{next}\
   {id=substr($0,1,15)}\
   $NF=="BIGB"{++bigB[id]}\
   $NF=="BIGMOVE"{++bigmove[id]}\
   $NF=="BADGEO"{++badgeo[id]}\
   $NF=="BADCONN"{print}\
   $NF=="BADRHO" && ( badgeo[id] || bigB[id]){\
      print id," |",badgeo[id]+0,bigmove[id]+0,bigB[id]+0,"BAD"}' |\
cat >! baddies.txt

# maybe less of these
@ topbad_oc  = ( $topbad / 2 + 1 )
head -n $topbad_oc badomega.txt badchir.txt |\
 awk '/^==/ || NF==0{next}\
   {print}' >> baddies.txt

set baddies = `cat baddies.txt | wc -l`
echo "$baddies baddies"
cp baddies.txt baddies_${itr}.txt

head -n $maxbaddies baddies.txt |\
cat - rholabeled.pdb |\
awk '$NF~/^BAD/{++bad[substr($0,1,15)];next}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,12,15)}\
   bad[id]{print}' >! removed_${itr}.pdb

head -n $maxbaddies baddies.txt |\
cat - rholabeled.pdb |\
awk '$NF~/^BAD/{++bad[substr($0,1,15)];next}\
  ! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,12,15);xyz=substr($0,31,24)}\
  seen[xyz]{next} {++seen[xyz]}\
  ! bad[id]{print substr($0,1,80)}' >! survived.pdb

set pdbfile = survived.pdb

if( ! $baddies ) set done = 1

end
echo "done with pruning loop at $itr"
set done = 0

echo "adding waters..."
rm -f new.pdb
add_waters.com survived.pdb centroids_${itr}.mtz sigma=4.5 >! addwater_${itr}.log
if( $status || ! -e new.pdb ) then
    tail addwater_${itr}.log 
    set BAD = "water addition step failed"
    goto exit
endif
mv new.pdb wetter.pdb 

set test = `egrep "^ATOM|^HETAT" new_water.pdb | wc -l`
echo "$test waters found"

# preemtively check if these will get eliminated anyway
rholabel_runme.com new_water.pdb centroids_${best_itr}.mtz mtzlabel=2FOFCWT outfile=new_water_rho.pdb >! rhoprobe_water.log

set rhos = `awk '{print $NF}' new_water_rho.pdb`
echo "DEBUG 2fofc : $rhos"
cat new_water_rho.pdb  |\
 awk -v minrho=$minrho '$NF<minrho*1.2 && /^ATOM|^HETAT/{print substr($0,12,15),$NF,"BADRHO"}' |\
 cat >! badrho_water.txt
cat badrho_water.txt wetter.pdb |\
awk '$NF~/^BAD/{++bad[substr($0,1,15)];next}\
  ! /^ATOM|^HETAT/{print;next}\
  {id=substr($0,12,15);xyz=substr($0,31,24)}\
  seen[xyz]{next} {++seen[xyz]}\
  ! bad[id]{print substr($0,1,80)}' >! wetter_safe.pdb
set test = `cat badrho_water.txt | wc -l`
echo "$test new water atoms in 2Fo-Fc density < $minrho"


set pdbfile = wetter_safe.pdb

combine_pdbs_runme.com xor=1 survived.pdb wetter_safe.pdb outfile=new_water_over_thresh.pdb > /dev/null

set test = `egrep "^ATOM|^HETAT" new_water_over_thresh.pdb | wc -l`
echo "$test waters added"
if( "$test" == "0" ) set done = 1

end
echo "done with water addition loop at $itr"
set done = 0

echo "clustering bonded fragments"
awk '/^BOND/' centroids_${itr}_fullgeo.txt |\
awk -F "|" '{gsub("_"," ");split($2,s,"-");\
  id1=substr(s[1],1,5) substr(s[1],7,1) substr(s[1],9,4) substr(s[1],14,1) substr(s[1],16,4);\
  id2=substr(s[2],1,5) substr(s[2],7,1) substr(s[2],9,4) substr(s[2],14,1) substr(s[2],16,4);\
  print id1 "|" id2 }' |\
tee bonded_atoms.txt |\
awk -F "|" '\
     clust[$1] && ! clust[$2]{clust[$2]=clust[$1];conn[$2]=$1;next}\
     clust[$2] && ! clust[$1]{clust[$1]=clust[$2];conn[$1]=$2;next}\
     clust[$1] && clust[$2] && clust[$1]!=clust[$2]{old=clust[$2];\
        for(ds in clust){if(clust[ds]==old)clust[ds]=clust[$1]};conn[$2]=conn[$2]" "$1;next}\
    ! clust[$1] && ! clust[$2]{\
        ++n;clust[$1]=clust[$2]=n;conn[$1]=$2;conn[$2]=conn[$2]" "$1;next}\
      END{for(ds in clust)if(clust[ds])print clust[ds],"|"ds}' |\
sort -g |\
awk -F "|" '! seen[$1]{newn[$1]=++n;++seen[$1]} {print newn[$1],"|"$2}' |\
awk -F "|" '{list[$1]=list[$1]"|"$2} END{for(l in list)print length(list[l]),list[l]}' |\
sort -gr |\
awk -F "|" '{print NF-1,"|" substr($0,index($0,$2))}' |\
tee cluster_lists.txt | wc -l

set test = `ls -1 | awk '/^cluster_/ && /.pdb$/' | wc -l`
if( $test ) then
  rm -f cluster_*.pdb
endif
rm -f cluster_CA_counts.txt
echo "cluster CAs"
echo -n "" >! liquify.pdb
foreach cluster ( `awk '{print NR}' cluster_lists.txt` )

head -n $cluster cluster_lists.txt | tail -n 1 |\
cat - $pdbfile |\
awk -F "|" 'NR==1{for(i=2;i<=NF;++i)++sel[$i];next}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,12,15)}\
  sel[id]{print}' |\
cat >! cluster_${cluster}.pdb 

set CAs = `filter_pdb.awk -v only=protein cluster_${cluster}.pdb | awk 'substr($0,12,5)=="  CA "' | wc -l`
echo $cluster $CAs | tee -a cluster_CA_counts.txt

if( $CAs < $min_CAs ) then
  filter_pdb.awk -v only=protein cluster_${cluster}.pdb >> liquify.pdb
endif
cp cluster_CA_counts.txt cluster_CA_counts_${itr}.txt

end

echo "combining cluster pdbs..."
rm -f unbonded.pdb
combine_pdbs_runme.com cluster_*.pdb xor=1 survived.pdb outfile=unbonded.pdb >! combine.log
if( $status || ! -e unbonded.pdb ) then
    tail combine.log 
    set BAD = "cluster combination step failed"
    goto exit
endif

echo "combining solid/liquid "
combine_pdbs_runme.com liquify.pdb unbonded.pdb xor=1 printref=1 $pdbfile outfile=solid.pdb >! combine.log

cat unbonded.pdb liquify.pdb |\
awk '/^ATOM|^HETAT/{print "ATOM      1  O   HOH S   1   " substr($0,30,40) "        O"}' |\
cat >> solid.pdb

egrep "^CRYST|HOH" solid.pdb |\
convert_pdb.awk -v renumber=1 >! oldwater.pdb

set nwaters = `egrep "^ATOM|^HEATAT" oldwater.pdb | wc -l`
rm -f rewatered.pdb
reorganize_waters.com oldwater.pdb 
if( ( $status || ! -e rewatered.pdb ) && $nwaters ) then
    tail combine.log 
    set BAD = "reorganization of waters step failed"
    goto exit
endif

egrep -v "HOH|^TER|^END" solid.pdb >! centroids_new.pdb
egrep "HOH" rewatered.pdb >> centroids_new.pdb

set prev = `ls -1rt | egrep "^refragged_" | awk '/.pdb$/' | tail -n 1`
cp centroids_new.pdb refragged_${itr}.pdb

if( -e "$prev" ) then
  set test = `rmsd $prev centroids_new.pdb | grep WARN | wc -l`
  echo "$test differences from last time"
  if( $test == "" ) set test = 0
  if( $test == "0" ) set done = 1
  if( $test != "0" ) set done = 0
endif

set pdbfile = centroids_new.pdb

# fluid loop
end
echo "done with liquidation loop at $itr"

# now pitch anything in bad density?


# creates phenix_opts_unbump.eff
no_new_nonbonds_runme.com $pdbfile $mtzfile $ciffiles >! nnb.log
if( $status || ! -e phenix_opts_unbump.eff ) then
    tail nnb.log
    set BAD = "failed to create new non-bond deactivation list at final step"
    goto exit
endif

echo "final refine"
rm -f centroids_final_001.pdb
#phenix.refine wxc_scale=$crush_scale wxu_scale=$crush_scale nonbonded_weight=$crush_nb \
phenix.refine  \
 $pdbfile $ciffiles \
  opts.eff phenix_opts_unbump.eff \
  $mtzfile prefix=centroids_final >! centroids_final_refine.log 
if( $status || ! -e centroids_final_001.pdb ) then
    tail centroids_final_refine.log 
    set BAD = "phenix.refine failed at final step"
    goto exit
endif

set pdbfile = centroids_final_001.pdb

filter_pdb.awk -v skip=water $pdbfile >! ${t}geotest.pdb

gemmi contact -d 1.2 --sort $pdbfile >! bad_contacts.txt
set test = `cat bad_contacts.txt | wc -l`
echo "$test non-bond contacts < 1.2A"
head -n 1 bad_contacts.txt

set restyps = `awk '/^ATOM|^HETAT/{print substr($0,18,3)}' $pdbfile | sort -u`
mkdir -p ${t}/list/
cp ${CLIBD_MON}/*.cif ${t}/
cp ${CLIBD_MON}/*.txt ${t}/

cat << EOF >! ${t}/list/mon_lib_list.cif
data_mon_lib_list

loop_
_lib_name
EOF
foreach restyp ( $restyps )
  echo "$restyp" >> ${t}/list/mon_lib_list.cif
  set a = `echo $restyp | awk '{print tolower(substr($0,1,1))}'`
  mkdir -p ${t}/${a}/
  set ciffile = ${restyp}.cif
  if(! -e $ciffile) set ciffile = ${CLIBD_MON}/${a}/$ciffile
  cp $ciffile ${t}/${a}/
end
awk '/^data_link_list/,""' ${CLIBD_MON}/list/mon_lib_list.cif >> ${t}/list/mon_lib_list.cif

gemmi rmsz --cutoff=6 --monomers=${t} ${t}geotest.pdb >! bad_geo.txt
set badomegas = `egrep "torsion CA-C-N-CA:" bad_geo.txt | wc -l`
echo "$badomegas peptide omega outliers"
set wrongchiral = `awk '/wrong chirality:/{print $3}' bad_geo.txt`
if("$wrongchiral" == "") set wrongchiral = "unknown"
echo "$wrongchiral inverted chirals"

exit:

if("$tempfile" == "") set  tempfile = "./"
set tempbase = `basename $tempfile`
set tempdir = `dirname $tempfile`
if(! $debug && ! ( "$tempdir" == "." && "$tempbase" == "" ) ) then
    rm -rf ${tempfile}*
endif

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif


exit

################################################################################

expand2supercell_runme.com centroids_asu.pdb refme_small.mtz super_mult=$super_mult \
  refpdb=bestgeo.pdb outprefix=centroids_super | tee centroids_expand.log

cp monomer_rot_trans.txt monomer_rot_trans_centroids.txt

expand2supercell_runme.com bestgeo.pdb refme_small.mtz super_mult=$super_mult \
  outprefix=bestgeo_super | tee bestgeo_expand.log

exit



rholabel_runme.com badgeo.pdb reference.mtz mtzlabel=Fref

filter_pdb.awk -v only=protein rholabeled.pdb |\
awk '! /^ATOM|^HETAT/{next}\
   {res=substr($0,22,8);B=substr($0,61,6);rho=$NF;\
    atom=substr($0,12,5);gsub(" ","",atom);\
    conf=substr($0,17,1)}\
    conf!=" "{next}\
    ! ( atom~/^[CNO]$/ || atom~/^C[AB]$/ || atom=="OXT" ){next}\
  {print}' >! mc.pdb
cat mc.pdb |\
awk '{res=substr($0,22,8);B=substr($0,61,6);rho=$NF;\
    atom=substr($0,12,5);gsub(" ","",atom)}\
   atom=="CA"{++hasca[res]} \
   {++count[res];sum[res]+=rho;sumB[res]+=B}\
  END{for(res in sum)if(hasca[res] && count[res]>1)\
     print res"|",sum[res]/count[res],"|",sumB[res]/count[res],"RESTAT"}' |\
cat - mc.pdb |\
awk '$NF=="RESTAT"{res=substr($0,1,8);\
     split($0,w,"|");avg[res]=w[2];avgB[res]=w[3]+0;next}\
 ! /^ATOM|HETAT/{next}\
   {res=substr($0,22,8);B=substr($0,61,6);rho=$NF}\
   avg[res]<0.1{next}\
   {++count[res];sum[res]+=(rho-avg[res])^2;sumB[res]+=(B-avgB[res])^2}\
  END{for(res in sum)if(count[res]>1){\
     rms=sqrt(sum[res]/count[res]);rmsB=sqrt(sumB[res]/count[res]);\
     print avg[res],rms,avgB[res],rmsB,"|",res}}' |\
 sort -gr | tee plotme


set medmad_rmsB = ``



# look for non-bonds that didn't used to be there
molprobify_runme.com starthere0.pdb keepgeo >&! molprobify_starthere0.log



awk '{print "PREV" $0}' starthere0_fullgeo.txt |\
 cat -  centroids_003_fullgeo.txt |\
 awk -F "|" '! /NONBOND/{next} /^PREV/{++seen[$2];next} \
 ! seen[$2]{print $2}' |\
 awk '{\
     a1=$1;f1=$2;c1=$4;r1=$5;\
     a2=$7;f2=$8;c2=$10;r2=$11;\
     v=1.5;s=99;\
   print "    #",c1,c2,a1,a2,r1,r2,f1,f2,v,s;\
   print "    bond {"\
   print "      action = *add";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      distance_ideal =",v;\
   print "      sigma = 99";\
   print "      slack = 6";\
   print "    }";}' |\
awk '{gsub("and altid _","");gsub("and chain _","");print}' |\
awk 'BEGIN{print "  geometry_restraints.edits {"} {print} END{print "  }"}' |\
awk 'BEGIN{print "refinement {"} {print} END{print "}"}' |\
cat >! phenix_opts_unbump.eff



