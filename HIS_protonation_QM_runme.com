#! /bin/tcsh -f
#
#   Determine HIS protonation states (HID / HIE / HIP) by QM using the Phenix
#   quantum interface, and write them out in the format the *_example_runme.com
#   scripts consume as HIS_settings_asu.txt: one line per HIS, "STATE Achain+resid"
#   (e.g. "HIP A176").
#
#   HID = proton on Ndelta only, HIE = proton on Nepsilon only,
#   HIP = both (positively charged).
#
#   VOTING: the QM is run on the main model AND on one model per alternate
#   conformation (altloc A, B, ...).  For each HIS the winning tautomer of each
#   run casts a vote; the majority wins (ties are broken by lowest QM energy).
#   Residues whose vote is split, or whose winner is within 1 kcal/mol of the
#   next tautomer, are flagged as close calls.  (For a HIS with no nearby altloc
#   every model agrees, so it is effectively a single determination.)
#
#   SLOW: runs mmtbx.quantum_interface run_qmr=True on every histidine of every
#   conformer model - a few minutes each (scaling with 1/nproc).
#
#   usage:
#     HIS_protonation_QM_runme.com pdbfile=starthere_asu.pdb ligcifs=EG7.cif nproc=6
#
#     pdbfile   starting model, ligand present; hydrogens are added by ready_set
#     ligcifs   ligand restraint cif(s) for the QM refinement
#               (default: every *.cif already in this directory)
#     nproc     qi.nproc for the QM refinement (default 3)
#     outfile   where to write the result (default HIS_settings.txt)
#
set pdbfile = ""
set ligcifs = ""
set nproc   = 3
set outfile = HIS_settings.txt

foreach Arg ( $* )
    set assign = `echo $Arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    if( $assign ) then
      # re-set any variable declared above by its own name
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = "$Val"
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "pdb")  set pdbfile = "$Val"
      if("$key" == "cif" || "$key" == "cifs") set ligcifs = "$Val"
      if("$key" == "np")   set nproc = "$Val"
    else
      # positional
      if("$Arg" =~ *.pdb ) set pdbfile = "$Arg"
      if("$Arg" =~ *.cif ) set ligcifs = "$ligcifs $Arg"
    endif
end

if( "$pdbfile" == "" || ! -e "$pdbfile" ) then
  echo "ERROR: need an existing pdbfile=  (starting model, ligand present)"
  exit 9
endif

# default: every ligand cif already present in this directory
if( "$ligcifs" == "" ) set ligcifs = `ls -1 *.cif |& awk '/\.cif$/'`

set prefix = $pdbfile:t
set prefix = $prefix:r

echo "HIS protonation by QM (conformer voting): pdbfile=$pdbfile prefix=$prefix ligcifs=($ligcifs) nproc=$nproc"

# 1. build the conformer models: the main model plus one per altloc label.
#    Each altloc model keeps the blank-altloc atoms plus that altloc's atoms,
#    with the altloc column blanked so it is a clean single conformation.
phenix.ready_set ${prefix}.pdb
if( ! -e ${prefix}.updated.pdb ) then
  echo "ERROR: ready_set did not produce ${prefix}.updated.pdb"
  exit 7
endif
set models = ( $prefix )
set altlocs = `awk '/^ATOM|^HETATM/{a=substr($0,17,1); if(a!=" ")print a}' ${prefix}.pdb | sort -u`
foreach al ( $altlocs )
awk -v al=$al '/^ATOM|^HETATM/{a=substr($0,17,1); if(a==" "||a==al) print substr($0,1,16)" "substr($0,18); next} {print}' ${prefix}.pdb >! ${prefix}_alt${al}.pdb
phenix.ready_set ${prefix}_alt${al}.pdb
if( -e ${prefix}_alt${al}.updated.pdb ) set models = ( $models ${prefix}_alt${al} )
end
echo "conformer models to vote across: $models"

# 2. run QM for every HIS in every conformer model; append one vote per (model,HIS).
#    Winner = lowest-relative-energy tautomer ("... ~>  0.00 kcal/mol"):
#      HD1 + HE2 -> HIP,  HD1 only -> HID,  HE2 only -> HIE   ("flipped" = same state)
#    Vote line: "resid state winner_abs_energy margin_to_next model"
rm -f HIS_votes.txt
foreach m ( $models )
echo 1 | mmtbx.quantum_interface iterate_NQH=HIS ${m}.updated.pdb | tee options_${m}.log
set opts = `awk '/resname HIS/ && $2==":"{print $1}' options_${m}.log`
foreach i ( $opts )
echo $i | mmtbx.quantum_interface iterate_NQH=HIS ${m}.updated.pdb | tee iterate_${m}_${i}.log
end
foreach phil ( `ls -1 ${m}.updated_*_HIS.phil` )
set resid = `echo $phil | awk -F '_' '{print $(NF-1)}'`
echo "  run_qmr HIS $resid  (model $m) ..."
mmtbx.quantum_interface ${m}.updated.pdb iterate_NQH=HIS $phil run_qmr=True qi.nproc=$nproc $ligcifs |& tee phil_${m}_${resid}.log
# Keep only the numbered ranking lines and drop the '!!!'/'><' summary markers by
# POSITIVELY matching "  N. ..." - a literal !!! in the grep pattern would trip
# tcsh history expansion ("0: Event not found") and kill the parse after the 1st HIS.
grep 'kcal/mol ~>' phil_${m}_${resid}.log | grep -E '^ *[0-9]+\. ' |\
awk -v r=$resid -v src=$m \
  '{rel=1e30;abs=1e30;\
    for(i=1;i<=NF;i++){if($i=="kcal/mol"&&abs==1e30)abs=$(i-1)+0; if($i=="~>")rel=$(i+1)+0}\
    s="";if(/HD1, HE2/)s="HIP";else if(/HD1 only/)s="HID";else if(/HE2 only/)s="HIE";\
    if(s=="" || rel==1e30)next;\
    if(!(s in smin) || rel<smin[s]){smin[s]=rel; sabs[s]=abs}}\
   END{best="";b1=1e30;b2=1e30;\
    for(st in smin){if(smin[st]<b1){b2=b1;b1=smin[st];best=st}else if(smin[st]<b2){b2=smin[st]}}\
    if(best!="")print r,best,sabs[best],(b2-b1),src}' |\
tee -a HIS_votes.txt
end
end

# 3. tally: majority state per residue, ties broken by lowest winner energy;
#    flag split votes and close calls.  Writes "STATE Achain+resid" to $outfile.
if( ! -s HIS_votes.txt ) then
  echo "ERROR: no HIS verdicts collected - check phil_*.log"
  exit 6
endif
sort -g HIS_votes.txt |\
awk -v out=$outfile 'function emit(){\
     if(cur==""){return}\
     best="";bv=-1;ba=1e30;nd=0;ss="";\
     for(s in cnt){nd++;ss=ss" "s; if(cnt[s]>bv||(cnt[s]==bv&&mab[s]<ba)){bv=cnt[s];best=s;ba=mab[s]}}\
     print best,"A"cur > out;\
     if(nd>1) printf("  NOTE: HIS A%s split vote (%s ) -> chose %s by lowest energy\n",cur,ss,best);\
     else if(mm<1.0) printf("  NOTE: HIS A%s close call - winner within %.2f kcal/mol of next tautomer\n",cur,mm);\
     delete cnt; delete mab; mm=1e30}\
   BEGIN{cur="";mm=1e30}\
   {if(cur!="" && $1!=cur) emit(); cur=$1;\
    s=$2;abs=$3+0;marg=$4+0;\
    cnt[s]++; if(!(s in mab)||abs<mab[s])mab[s]=abs; if(marg<mm)mm=marg}\
   END{emit()}'

echo ""
echo "wrote $outfile :"
cat $outfile
