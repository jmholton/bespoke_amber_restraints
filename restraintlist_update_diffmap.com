#! /bin/tcsh -f
#
#  use x-ray difference map to update restraint weights                     -James Holton  3-29-24
#
#
#
set pdbfile = ""
set trajectory = ""
set refpointspdb = ""
set diffmap = ""
set refmap  = "reference.map"
set mtzfile = ""

set GRID = ""

set outfile = "new_restraints.pdb"
set outmults = "sorted_mults.txt"

set tempfile = /dev/shm/${USER}/temp_rud_$$_
#set tempfile = ./tempfile_rud_
mkdir -p /dev/shm/${USER}
mkdir -p ${CCP4_SCR}
set logfile = details.log

# for printing
set modulo = 10000

set overall_scale = 1
set negative_scale = 0.9
set fft_B = 0
set shan_B = auto
set halfrho_pos = auto
set halfrho_neg = auto
set thresh_refmap = 1.0
set elim_refmap = 0.0
set max_mult = 2.0
set max_weight = 999.99

set ambig_same_weight = 0

set smallSG = ""
set smallSGnum = ""
set smallCELL = ""
set GRID = ""

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
      if("$key" == "restpdb") set refpointspdb = "$Val"
      if("$key" == "refpoints") set refpointspdb = "$Val"

      if("$key" == "outfile") set outfile = "$Val"

      if("$key" == "fft_b") set fft_B = "$num"
    else
      # no equal sign
      if("$key" =~ *.in ) set infile = "$Arg"
      if("$key" =~ *.txt ) set restraint_list = "$Arg"
    
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
      if("$Arg" =~ *.mtz ) set mtzfile = $Arg
      if("$Arg" =~ *.map ) set diffmap = $Arg
    endif
    if("$key" == "debug") set debug = "$Val"
end

if(! -e "$pdbfile") then
    set BAD = "no coordinates provided"
    goto exit
endif
if(! -e "$refpointspdb") then
    set BAD = "no reference-points pdb provided"
    goto exit
endif

#if( $?debug && "$tempfile" =~ /dev/shm/* ) set tempfile = tempfile_rud_

set CPUs = `grep proc /proc/cpuinfo | wc -l | awk '{print int($1/4)}'`
if( "$CPUs" == "" ) set CPUs = 1
set thishost = `hostname -s`
# use cluster or not?
# cannot migrate hosts because of temp files
set test = `sinfo -h -n $thishost |& egrep -v "drain|n/a" | awk '$2=="up"' | wc -l`
if ( $test ) then
  echo "using slurm"
  set srun = "srun -w $thishost"
  set test = `echo $tempfile | awk '{print ( ! /\/dev\/shm/ )}'`
  if( $test ) then
    echo "full cluster"
    set srun = "srun"
  endif
else
  set srun = ""
endif

set t = "$tempfile"

touch $logfile


if(-e "$mtzfile") then
  set test = `echo head | mtzdump hklin $mtzfile | egrep -i fpart | wc -l`
  set mtzreso = `echo head | mtzdump hklin $mtzfile | awk '/Resolution Range/{getline;getline;print $6}'`

  echo "2fofc = refmap"
  fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
  labin F1=FWT PHI=PHWT
  $GRID
EOF
  # always do sigma scale?
  echo scale sigma 1 0 |\
  mapmask mapin ${t}ffted.map mapout ${t}refmap_sigma.map >> $logfile


  echo "fofc with B=$fft_B"
  fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
  labin F1=DELFWT PHI=PHDELWT
  scale F1 1 $fft_B
  reso 1
EOF
  echo scale sigma 1 0 |\
  mapmask mapin ${t}ffted.map mapout ${t}fofc_sigma.map >> $logfile
  echo | mapdump mapin ${t}fofc_sigma.map | grep "mum dens"

  echo "back-and-forth fft to purify symmetry"
  mapmask mapin1 ${t}ffted.map mapout ${t}extended.map << EOF >> $logfile
  xyzlim cell
  axis Z X Y
EOF
  sfall mapin ${t}extended.map hklout ${t}sfalled.mtz << EOF >> $logfile
  mode sfcalc mapin
  sfsg 1
  reso 1
EOF
  fft hklin ${t}sfalled.mtz mapout ${t}reffted.map << EOF >> $logfile
  labin F1=FC PHI=PHIC
  scale F1 0.5 0
  $GRID
  reso $mtzreso
EOF

  # always do sigma scale?
  echo scale sigma 1 0 |\
  mapmask mapin ${t}reffted.map mapout ${t}fofc_sigma.map >> $logfile

  set diffmap = ${t}fofc_sigma.map
endif


if(-e "$refmap") then
  echo scale sigma 1 0 |\
  mapmask mapin $refmap mapout ${t}refmap_sigma.map >> $logfile
endif

if(! -e ${t}refmap_sigma.map) then
    set BAD = "no reference (2fofc) map file: $refmap "
    goto exit
endif

if(! -e "$diffmap" ) then
    set BAD = "no map file: $diffmap "
    goto exit
endif

echo "" | mapdump mapin $diffmap >! ${t}mapdump.txt
set GRID = `awk '/Grid sampling/{print "GRID", $(NF-2), $(NF-1), $NF; exit}' ${t}mapdump.txt`
set AXIS = `awk '/Fast, medium, slow/{print "AXIS", $(NF-2), $(NF-1), $NF}' ${t}mapdump.txt | awk '$NF !~ /[^ XYZ]/'`
set LIMITS = `awk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print "XYZLIM",b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${t}mapdump.txt`
set smallCELL = `awk '/Cell dimensions /{print $4,$5,$6,$7,$8,$9;exit}' ${t}mapdump.txt`
set smallSGnum  = `awk '/ Space-group /{print $NF;exit}' ${t}mapdump.txt`
set smallSG = `awk -v n=$smallSGnum '$1==n{print $4;exit}' ${CLIBD}/symop.lib`
set nsymops = `awk -v n=$smallSGnum '$1==n{print $3;exit}' ${CLIBD}/symop.lib`

grep "mum dens" ${t}mapdump.txt


cat << EOF
pdbfile = $pdbfile
trajectory = $trajectory
refpoints = $refpointspdb
mtzfile = $mtzfile
diffmap = $diffmap
refmap  = $refmap

tempfile = $tempfile
debug    = $?debug
EOF

# ignore hydrogens here
#convert_pdb.awk -v skip=H,EP $pdbfile >! ${t}noH.pdb
awk 'substr($0,77,2) !~ / H|XP/' $pdbfile >! ${t}noH.pdb

# extract atoms that have restraints ... aka are listed in the refpoints file
# but make sure they appear in the order that they appear in the refpoints file
combine_pdbs_runme.com ${t}noH.pdb $refpointspdb \
  outfile=${t}hasref.pdb >> $logfile

# make sure order and number of restrained atoms  is identical to that of hasref file
combine_pdbs_runme.com ${t}noH.pdb $refpointspdb \
  saveXYZ \
  outfile=${t}refpoints.pdb >> $logfile

echo "hasref vs refpoints"
rmsd ${t}hasref.pdb ${t}refpoints.pdb | grep -v Bfac

echo "" >! ${t}empty.pdb
pdbset xyzin ${t}empty.pdb xyzout ${t}cell.pdb << EOF >> $logfile
CELL $smallCELL
SPACE $smallSG
EOF

egrep "^CRYST" ${t}cell.pdb >! ${t}wrapme.pdb
egrep "^ATOM|^HETAT" ${t}hasref.pdb >> ${t}wrapme.pdb

echo "wrapping $pdbfile into small cell"
wrap_into_cell.com ${t}wrapme.pdb \
  doprotein=1 breakbonds=1 \
  outfile=${t}hasref_small.pdb >> $logfile

egrep "^CRYST" ${t}cell.pdb >! ${t}wrapme.pdb
egrep "^ATOM|^HETAT" ${t}refpoints.pdb >> ${t}wrapme.pdb

echo "wrapping $refpointspdb into small cell"
wrap_into_cell.com ${t}wrapme.pdb \
  doprotein=1 breakbonds=1 \
  outfile=${t}refpoints_small.pdb >> $logfile

#echo "hasref vs refpoints"
#rmsd ${t}hasref_small.pdb ${t}refpoints_small.pdb | grep -v Bfac

# sanity check
set refatoms = `egrep "^ATOM|^HETAT" ${t}refpoints_small.pdb | wc -l`
set xyzatoms = `egrep "^ATOM|^HETAT" ${t}hasref_small.pdb | wc -l`
if( $refatoms != $xyzatoms ) then
   set BAD = "atom counts for ref and xyz do not match: $refatoms vs $xyzatoms"
   goto exit
endif

echo "peek refpoints"
check_map_sites.com ${t}refpoints_small.pdb $diffmap |\
  awk '{print ++n,$1,"FOFCref"}' >! ${t}fofc_ref.txt

if(! -e "$trajectory") then
  echo "probing xyz from $pdbfile"
  check_map_sites.com ${t}hasref_small.pdb $diffmap |\
    awk '{print ++n,$1,"FOFCxyz"}' >! ${t}fofc_xyz.txt
  goto skiptraj
endif



echo "probing density in trajectory"
echo -n "" >! ${t}fofc_traj.txt
# assume nothing other than order of atoms in each pdb is the same
awk '{print $0,"HASREF"}' ${t}hasref_small.pdb |\
cat - $pdbfile |\
awk '! /^ATOM|^HETAT/{next}\
     {id=substr($0,12,17);}\
     $NF=="HASREF"{++hasref[id];next}\
     {++n}\
     hasref[id]{print $0,n}' >! ${t}hasref_atomnums.pdb

cat << EOF >! ${t}job.csh
#! /bin/tcsh -f
set n = "\$1"
set t = ${t}_\$\$_
set pdb = ${trajectory}/md.\${n}.pdb
awk '{print \$NF}' ${t}hasref_atomnums.pdb |\
    cat - \$pdb |\
    awk 'NF==1{++sel[\$1]} ! /^ATOM|^HETAT/{next}\
      /EPW/{next}\
      {++n} sel[n]{print}' |\
    cat >! \${t}peekme.pdb
    check_map_sites.com \${t}peekme.pdb $diffmap |\
    awk '{print ++n,\$1}'
EOF
chmod a+x ${t}job.csh

foreach pdb ( ${trajectory}/md.*.pdb )
    set n = `echo $pdb | awk -F "." '{print $2}'`
    echo "peek $pdb"
    $srun ${t}job.csh $n >! ${t}peek_${n}.txt &

    if( "$srun" == "" ) then
      @ m = ( $n % $CPUs )
      if( $m == 0 ) wait
    endif
end
wait
echo "averaging results"
cat ${t}peek_*.txt |\
awk '{sum[$1]+=$2;++count[$1]}\
  END{max=$1;for(n=1;n<=max;++n)if(count[n])print n,sum[n]/count[n],count[n],"FOFCxyz"}' |\
cat >! ${t}fofc_xyz.txt

skiptraj:

touch ${t}refmap_ref.txt
if(-e ${t}refmap_sigma.map) then
 check_map_sites.com ${t}refpoints_small.pdb ${t}refmap_sigma.map |\
  awk '{print ++n,$1,"REFMAPref"}' >! ${t}refmap_ref.txt
endif

# sanity check
set refrhos = `cat ${t}fofc_ref.txt | wc -l`
set xyzrhos = `cat ${t}fofc_xyz.txt | wc -l`
if( $refrhos != $xyzrhos ) then
   set BAD = "map peek counts for ref and xyz do not match: $refrhos vs $xyzrhos"
   goto exit
endif
if( $refrhos != $refatoms ) then
   set BAD = "map peek counts do not match atom count for ref: $refrhos vs $refatoms"
   goto exit
endif
if( $xyzrhos != $xyzatoms ) then
   set BAD = "map peek counts do not match atom count for input xyz: $xyzrhos vs $xyzatoms"
   goto exit
endif

# combine into one list
cat ${t}fofc_???.txt ${t}refmap_ref.txt ${t}signif_???.txt |\
awk '$NF=="FOFCref"{fr[++n]=$(NF-1)}\
     $NF=="REFMAPref"{tfr[++m]=$(NF-1)}\
     $NF=="FOFCxyz"{fx[++p]=$2}\
  END{for(i=1;i<=n;++i)print i,fr[i]+0,fx[i]+0,tfr[i]+0}' |\
cat >! ${t}rhos.txt 
# line fofc_ref  fofc_xyz  refmap_ref 
#  1     2          3         4      


if( "$halfrho_pos" == "auto" || "$halfrho_neg" == "auto" ) then
  # number of shannon voxels expected
  if( "$shan_B" == "auto") then
    set shan_B = `echo $fft_B 2 | awk '{print $1+$2}'`
  endif
  echo $smallCELL |\
  awk 'NF==6{s=atan2(1,1)/45; A=cos(s*$4); B=cos(s*$5); G=cos(s*$6); \
    skew = 1 + 2*A*B*G - A*A - B*B - G*G ; if(skew < 0) skew = -skew;\
    printf("%.3f",$1*$2*$3*sqrt(skew))}' |\
  cat >! ${t}volume
  set smallCELLvolume = `cat ${t}volume`

  set exponent = `echo $smallCELLvolume $nsymops $shan_B | awk '{V=$1;n=$2;B=$3;pi=atan2(1,1)*4;print V/n/2/(sqrt((B+9.484)*log(2))/pi)**3}'`
  echo "exponent= $exponent"

  # ' "probability that something is there" '
  # ' prob(rho) = sign(rho)*pow(erf(abs(rho/sigma(rho))/sqrt(2)),V/2/d^3) '

  # come up with smoother function
  gnuplot << EOF >&! ${t}gnuplot.log 
  exponent=$exponent
  sigma=1.0
  # sigma level that is 50% likely to be something
  halfrho = sqrt(2)*inverf((0.5)**(1./exponent))
  print halfrho
EOF
  set halfrho = `tail -n 1 ${t}gnuplot.log | awk '{print $NF}'`
endif
if( ! $?halfrho ) set halfrho
if( "$halfrho" == "" ) set halfrho = 5
if( "$halfrho_pos" == "auto") set halfrho_pos = $halfrho
if( "$halfrho_neg" == "auto") set halfrho_neg = $halfrho


echo "half-rho: $halfrho_pos $halfrho_neg"

# convert density probes into weight multipliers
# if ref is positive and xyz is negative, increase weight
# if ref is negative and xyz is positive, decrease weight
# if ref is positive and xyz is positive, look at delta
# if ref is negative and xyz is negative, decrease weight
# if ref (2fofc) is below threshold and ref is negative, no action
echo "$halfrho_pos $halfrho_neg $thresh_refmap $max_mult $negative_scale" |\
cat - ${t}rhos.txt |\
awk 'NR==1{halfrho_pos=$1;halfrho_neg=$2;thresh_refmap=$3;max_mult=$4;negative_scale=$5;next}\
   {a=$1;fofc_ref=$2;fofc_xyz=$3;twofofc_ref=$4;\
    delta_fofc=fofc_ref-fofc_xyz;\
    # transform sigmas into significance \
    signif_ref=(fofc_ref/halfrho_pos)**5;\
    signif_xyz=(fofc_xyz/halfrho_pos)**5;\
    signif_del=(delta_fofc/halfrho_pos)**5;}\
    # allow different settings for pos and neg sidebands \
    fofc_ref<0{signif_ref=(fofc_ref/halfrho_neg)**5}\
    fofc_xyz<0{signif_xyz=(fofc_xyz/halfrho_neg)**5}\
    delta_fofc<0{signif_del=(delta_fofc/halfrho_neg)**5}\
    # use reference score as overall down-weight \
    {squelch=signif_ref}\
    squelch>0{squelch=0}\
    # clip to [-1:1] \
    signif_ref>1{signif_ref=1} signif_ref<-1{signif_ref=-1}\
    signif_del>1{signif_del=1} signif_del<-1{signif_del=-1}\
   {signif=signif_del+squelch}\
    signif>1{signif=1} signif<-1{signif=-1}\
    # compute multipliers \
   {mult=max_mult**(signif)}\
   # low 2fofc and reference is negative -> no action \
   twofofc_ref<thresh_refmap && fofc_ref<0{mult=1}\
   # reference is negative -> neg scale \
   signif_ref<0{mult*=negative_scale**(-signif_ref)}\
   {print $0,signif_xyz,signif_ref,mult,"MULT"}' |\
cat >! ${t}mults.txt
# line fofc_ref fofc_xyz   refmap signif_xyz signif_ref weight_mult MULT
#  1     2          3         4        5         6           7         8


# print out extrema - but label with rank of previous weight
awk '/^ATOM|^HETAT/{print $0,++n}' $refpointspdb |\
sort -k1.61gr |\
awk '{print $NF,$0,++n}' |\
sort -g |\
awk '{print substr($0,match($0,/ATOM|HETAT/))}' >! ${t}labeled.pdb

echo $modulo |\
cat - ${t}mults.txt ${t}labeled.pdb |\
awk 'NR==1{modulo=$1}\
  $NF=="MULT"{mult[$1]=$(NF-1);ref[$1]=$2;xyz[$1]=$3;twofofc_ref[$1]=$4;next}\
  ! /^ATOM|^HETAT/{next}\
  {++n;id=substr($0,12,11)" "substr($0,23,9);rank=$NF;\
      resnum=substr($0,23,8);\
      occ=substr($0,55,6);B=substr($0,61,6)}\
  {printf("%s %4d %6d : %.2f x %6.2f fofcs: %.2f %.2f  2fofc: %.2f rank: %d\n",id,(resnum-1)%modulo+1,n,mult[n],B,ref[n],xyz[n],twofofc_ref[n],rank)}' |\
sort -t ":" -k2gr  >! ${t}sorted.txt

echo "highest:"
awk '! seen[$5]{print;++seen[$5]}' ${t}sorted.txt | head -n 3
echo "lowest:"
tail -n 1 ${t}sorted.txt
# atomid resnum%modulo ordatom : weightmult fofc_ref fofc_xyz erfmap_ref

echo "saving weight multipliers in: $outmults"
cp ${t}sorted.txt $outmults





echo "updating weights x$overall_scale into $outfile"
echo $overall_scale $thresh_refmap $elim_refmap $max_weight  |\
cat - ${t}mults.txt $refpointspdb |\
awk 'NR==1{overall_scale=$1;thresh_refmap=$2;elim_refmap=$3;max_weight=$4;next}\
  $NF=="MULT"{mult[$1]=$(NF-1);twofofc_ref[$1]=$4;next}\
  ! /^ATOM|^HETAT/{print;next}\
  {++n;pre=substr($0,1,60);B=substr($0,61,6)+0;post=substr($0,67);\
   dw=overall_scale*mult[n];\
   if(B==0 && dw>1)B=0.01;\
   weight=B*dw;\
     newB=sprintf("%6.2f",weight)+0}\
  newB==B && dw<1{newB=newB-0.01}\
  newB==B && dw>1 && twofofc_ref[n]>=thresh_refmap{newB=B+0.01}\
  twofofc_ref[n]<elim_refmap{newB=0}\
  newB<0{newB=0}\
  newB>999.99{newB=999.99}\
  newB>max_weight{newB=max_weight}\
  {printf("%s%6.2f%s\n",pre,newB,post)}' |\
cat >! $outfile


if( $ambig_same_weight ) then
    echo "enforcing x-ray ambiguous atoms to have same restraint weight"
    xsame_runme.com infile=$outfile outfile=${t}xsame_restraints.pdb >&! ${t}xsame.log
    if($status) then
       set BAD = "xsame failed"
    endif
    egrep ".D. ASN" $outfile | head -n 2
    mv ${t}xsame_restraints.pdb $outfile
    egrep ".D. ASN" $outfile | head -n 2
endif


exit:

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?debug && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    echo "clearing temp files"
    rm -f ${t}*
endif


exit

