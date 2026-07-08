#! /bin/tcsh -f
#
#  use x-ray difference map to update B factors                     -James Holton  8-12-24
#
#
#
set Bfacfile = ""
set trajectory = "trajectory/"
set mtzfile = "cootme.mtz"

set GRID = ""

set outfile = "new_Bfac.pdb"
set outmods = "sorted_Bmods.txt"

set tempfile = /dev/shm/${USER}/temp_Bud_$$_
#set tempfile = ./tempfile_Bud_
mkdir -p /dev/shm/${USER}
mkdir -p ${CCP4_SCR}
set logfile = details.log

# for printing
set modulo = 10000

set overall_scale = 1
set overall_offset = 0
set negative_scale = 0.9
set fft_B = 0
set shan_B = auto
set Wilson_B = 2
set halfrho_pos = auto
set halfrho_neg = auto
set max_mod = 2.0
set mod_mode = add
set max_Bfac = 999.99
set min_Bfac = 2

set ambig_same_Bfac = 1

set mtzlabel = DELFWT
set phenixlabel = miller_array.labels.name
set test = `phenix.version | awk '/Release tag/{print ( $NF < 5000 )}'`
if( "$test" == "1" ) set phenixlabel = label

set sruncpu = ""

set smallSG = ""
set smallSGnum = ""
set smallCELL = ""

foreach sourceme ( compute_settings.sourceme xtal_properties.sourceme user_settings.sourceme )
    if(-e $sourceme ) then
      echo "sourcing $sourceme"
      source $sourceme
    else
      if(-e ../$sourceme ) then
        echo "sourcing ../$sourceme"
        source ../$sourceme
      endif
    endif
end

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
    
      if("$Arg" =~ *.pdb ) set Bfacfile = $Arg
      if("$Arg" =~ *.mtz ) set mtzfile = $Arg
      if("$Arg" =~ *.map ) set diffmap = $Arg
    endif
    if("$key" == "debug") set debug = "$Val"
end

if(! -e "$Bfacfile") then
    set BAD = "no coordinates provided"
    goto exit
endif

#if( $?debug && "$tempfile" =~ /dev/shm/* ) set tempfile = tempfile_Bud_

set t = "$tempfile"

touch $logfile

if(! -e "$trajectory") then
    set BAD = "need a trajectory in $trajectory"
    goto exit
endif

if(! -e "$mtzfile" ) then
    set BAD = "no mtz file: $mtzfile "
    goto exit
endif

#set test = `echo head | mtzdump hklin $mtzfile | egrep -i fpart | wc -l`
set mtzreso = `echo head | mtzdump hklin $mtzfile | awk '/Resolution Range/{getline;getline;print $6}'`

set fftreso = `echo $mtzreso 1.2 | awk '{print ($2*($1^-3))^(-1/3)}'`

echo "fofc with B=$fft_B"
fft hklin $mtzfile mapout ${t}ffted.map << EOF >> $logfile
labin F1=$mtzlabel PHI=PHDELWT
scale F1 1 $fft_B
reso $fftreso
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
cad hklin1 ${t}sfalled.mtz hklout ${t}diffmap.mtz << EOF >> $logfile
labin file 1 E1=FC E2=PHIC
labou file 1 E1=$mtzlabel E2=PHDELWT
reso over_all $mtzreso
EOF



if(! -s "${t}diffmap.mtz" ) then
    set BAD = "map reduction to $smallSG failed "
    goto exit
endif

echo "" | mapdump mapin ${t}ffted.map >! ${t}mapdump.txt
set GRID = `awk '/Grid sampling/{print "GRID", $(NF-2), $(NF-1), $NF; exit}' ${t}mapdump.txt`
set AXIS = `awk '/Fast, medium, slow/{print "AXIS", $(NF-2), $(NF-1), $NF}' ${t}mapdump.txt | awk '$NF !~ /[^ XYZ]/'`
set LIMITS = `awk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print "XYZLIM",b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${t}mapdump.txt`
set smallCELL = `awk '/Cell dimensions /{print $4,$5,$6,$7,$8,$9;exit}' ${t}mapdump.txt`
set smallSGnum  = `awk '/ Space-group /{print $NF;exit}' ${t}mapdump.txt`
set smallSG = `awk -v n=$smallSGnum '$1==n{print $4;exit}' ${CLIBD}/symop.lib`
set nsymops = `awk -v n=$smallSGnum '$1==n{print $3;exit}' ${CLIBD}/symop.lib`



cat << EOF
Bfacfile = $Bfacfile
trajectory = $trajectory
mtzfile = $mtzfile
fft_B = $fft_B

tempfile = $tempfile
debug    = $?debug
EOF

set CPUs = `grep proc /proc/cpuinfo | wc -l | awk '{print int($1/4)}'`
if( "$CPUs" == "" ) set CPUs = 1
set thishost = `hostname -s`
# use cluster or not?
# maybe cannot migrate hosts because of temp files
set test = `sinfo -h -n $thishost |& egrep -v "drain|n/a" | awk '$2=="up"' | wc -l`
if ( $test ) then
  echo "using slurm"
  if( "$sruncpu" == "" ) set sruncpu = "srun"
  set srun = "$sruncpu -w $thishost"
  set trajdir = `echo ${trajectory} | awk '{gsub("/$","");print}'`
  set trajdir = `ls -ld ${trajdir} | awk '{print $NF}'`
  if( "$trajdir" != "" && -e "$trajdir" ) then
    echo "actual trajectory: $trajdir"
    set trajtest = `echo $trajdir | awk '{print ( ! /\/dev\/shm/ )}'`
  endif
  set temptest = `echo $tempfile | awk '{print ( ! /\/dev\/shm/ )}'`
  if( $trajtest && $temptest ) then
    echo "full cluster"
    set srun = "$sruncpu"
  endif
  if( $debug ) set srun = "$srun -p debug"
else
  set srun = ""
endif

set t = $tempfile

set pdbSG = `awk -v SGnum="$smallSGnum" -F "[\047]" '$1+0==SGnum{print $2;exit}' ${CLIBD}/symop.lib`
echo $smallCELL $pdbSG |\
awk '{printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)}' |\
cat >! ${t}smallcell.pdb


cat << EOF >! ${t}peek.csh
#! /bin/tcsh -f
set pdb = "\$1"
set t = /dev/shm/${USER}/temp_\$\$_
mkdir -p /dev/shm/${USER}
if( $debug ) set t = ${t}_temp_\$\$_
#set t = tempfile_temp_
cp ${t}smallcell.pdb \${t}.pdb
awk '/EPW|Y1  HOH|Y 1  HOH/{next} \
 /^ATOM|^HETAT/ && substr(\$0,77,2) !~ / H|XP| Y/' \$pdb |\
tee -a \${t}.pdb |\
awk '{x=substr(\$0,31,8);y=substr(\$0,39,8);z=substr(\$0,47,8);\
  printf("(%.3f,%.3f,%.3f) %d KEY\n",x,y,z,++n)}' >! \${t}key.txt
phenix.map_value_at_point ${t}diffmap.mtz \${t}.pdb \
          ${phenixlabel}=$mtzlabel scale=sigma |\
cat \${t}key.txt - |\
awk '\$NF=="KEY"{n[\$1]=\$2;next}\
  /Map value:/ && \$3~/^\(/{print n[\$3],\$3,++i,\$NF;next}\
  /Map value:/ && /^"/{point=substr(\$0,index(\$0," Point: "));\
    split(point,w);xyz="(" w[2] w[3] w[4] ")";\
    print n[xyz],xyz,++i,\$NF}'
rm -f \${t}*
EOF
chmod a+x ${t}peek.csh

echo "probing density in trajectory"
# assume nothing other than order of atoms in each pdb is the same
foreach pdb ( ${trajectory}/md.*.pdb )
    set n = `echo $pdb | awk -F "." '{print $2}'`
    echo "peek $pdb"
    $srun ${t}peek.csh $pdb >! ${t}peek_${n}.txt &

    if( "$srun" == "" ) then
      @ m = ( $n % $CPUs )
      if( $m == 0 ) wait
    endif
end
wait

echo "averaging results"
cat ${t}peek_*.txt |\
awk '{sum[$1]+=$NF;++count[$1]} $1>max{max=$1+0}\
  END{for(n=1;n<=max;++n)if(count[n])print n,sum[n]/count[n],count[n],"FOFC"}' |\
cat >! ${t}avgfofc.txt
if( ! -s ${t}avgfofc.txt) then
  set BAD = "no peek results"
  goto exit
endif

# sanity checks?
set test = `awk '{print $3}' ${t}avgfofc.txt | sort -u | wc -l`
if( "$test" != "1" ) then
  set BAD = "non-equal peek counts"
  goto exit
endif


# number of shannon voxels expected
if( "$shan_B" == "auto") then
  set shan_B = `echo $fft_B $Wilson_B | awk '{print $1+$2}'`
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
if( "$halfrho" == "" ) set halfrho = 5
if( "$halfrho_pos" == "auto") set halfrho_pos = $halfrho
if( "$halfrho_neg" == "auto") set halfrho_neg = $halfrho


echo "half-rho: $halfrho_pos $halfrho_neg"

# convert density probes into B factor modifiers
# if fofc is positive decrease B factor
# if fofc is negative increase B factor
# if fofc is negative always increase B factor a bit
echo "$halfrho_pos $halfrho_neg $max_mod $mod_mode $negative_scale" |\
cat - ${t}avgfofc.txt |\
awk 'NR==1{halfrho_pos=$1;halfrho_neg=$2;max_mod=$3;mod_mode=$4;;negative_scale=$5;next}\
   {a=$1;fofc=$2;\
    # transform sigmas into significance \
    signif=(fofc/halfrho_pos)**5;}\
    # allow different settings for pos and neg sidebands \
    fofc<0{signif=(fofc/halfrho_neg)**5}\
    # clip to [-1:1] \
    signif>1{signif=1} signif<-1{signif=-1}\
    # compute modifier \
   mod_mode=="mult"{mod=max_mod**(-signif)}\
   mod_mode=="add"{mod=max_mod*-signif}\
   # reference is negative -> neg scale \
   signif<0{mod*=negative_scale**(signif)}\
   {print $1,$2,signif,mod,"MOD"}' |\
cat >! ${t}mods.txt
# line fofc  signif modifier    MOD
#  1     2    3         4        5  

# only non-hydrogen at this point
filter_pdb.awk -v skip=H -v only=atoms $Bfacfile >! ${t}B_noH.pdb

# sanity check?
# there ought to be more non-H atoms in Bfac.pdb than flying model


# print out extrema - but label with rank of previous B
awk '/^ATOM|^HETAT/{print $0,++n}' ${t}B_noH.pdb |\
sort -k1.61gr |\
awk '{print $NF,$0,++n}' |\
sort -g |\
awk '{print substr($0,match($0,/ATOM|HETAT/))}' >! ${t}labeled.pdb
# ATOM xxx order rank

echo $modulo $mod_mode |\
cat - ${t}mods.txt ${t}labeled.pdb |\
awk 'NR==1{modulo=$1;mod_mode=$2}\
  {op="x"} mod_mode=="add"{op="+"}\
  $NF=="MOD"{mod[$1]=$(NF-1);fofc[$1]=$2;next}\
  ! /^ATOM|^HETAT/{next}\
  {++n;id=substr($0,12,11)" "substr($0,23,9);rank=$NF;\
      resnum=substr($0,23,8);\
      occ=substr($0,55,6);B=substr($0,61,6)}\
  fofc[n]==""{next}\
  {printf("%s %4d %6d : %.2f %s %6.2f fofc: %.2f  rank: %d\n",id,(resnum-1)%modulo+1,n,mod[n],op,B,fofc[n],rank)}' |\
sort -t ":" -k2gr  >! ${t}sorted.txt

echo "highest:"
awk '! seen[$5]{print;++seen[$5]}' ${t}sorted.txt | head -n 3
echo "lowest:"
tail -n 1 ${t}sorted.txt
# atomid resnum%modulo ordatom : modifier fofc_ref fofc_xyz erfmap_ref

echo "saving weight modifiers in: $outmods"
cp ${t}sorted.txt $outmods





echo "updating B factors x$overall_scale +$overall_offset into $outfile"
echo $overall_scale $overall_offset $mod_mode $max_Bfac $min_Bfac |\
cat - ${t}mods.txt ${t}B_noH.pdb |\
awk 'NR==1{overall_scale=$1;overall_offset=$2;mod_mode=$3;max_Bfac=$4;min_Bfac=$5;next}\
  $NF=="MOD"{mod[$1]=$(NF-1);next}\
  ! /^ATOM|^HETAT/{print;next}\
  {++n;pre=substr($0,1,60);B=substr($0,61,6)+0;post=substr($0,67);\
   dB=overall_scale*mod[n]+overall_offset}\
   B==0 && dB>1{B=0.01}\
   mod_mode=="mult"{newB=B*dB}\
   mod_mode=="add"{newB=B+dB}\
  {newB=sprintf("%6.2f",newB)+0}\
  newB==B && dB<1{newB=newB-0.01}\
  newB==B && dB>1{newB=B+0.01}\
  newB<min_Bfac{newB=min_Bfac}\
  newB>max_Bfac{newB=max_Bfac}\
  newB>999.99{newB=999.99}\
  {printf("%s%6.2f%s\n",pre,newB,post)}' |\
cat >! ${t}newB_noH.pdb

echo "placing new B factors in $Bfacfile"
combine_pdbs_runme.com ${t}newB_noH.pdb $Bfacfile printref=1 \
  outfile=${t}heavyB.pdb >! ${t}combine.log

if( $ambig_same_Bfac ) then
    echo "enforcing x-ray ambiguous atoms to have same B factor"
    cp ${t}heavyB.pdb ${t}heavyB_prexsame.pdb
    xsame_runme.com infile=${t}heavyB.pdb outfile=${t}xsame_restraints.pdb >&! ${t}xsame.log
    if($status) then
       set BAD = "xsame failed"
       goto exit
    endif
    egrep ".D. ASN" ${t}heavyB.pdb | head -n 2
    mv ${t}xsame_restraints.pdb ${t}heavyB.pdb
    egrep ".D. ASN" ${t}heavyB.pdb | head -n 2
endif

echo "propagating B factors to hydrogens via bonds"
hsame_runme.com ${t}heavyB.pdb outfile=$outfile
if($status) then
    set BAD = "hsame failed"
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

