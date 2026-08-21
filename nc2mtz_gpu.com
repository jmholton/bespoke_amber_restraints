#! /bin/tcsh -f
#
# convert an amber nc trajectory to a time-averaged mtz (FCavg/PHICavg), using
# the GPU structure-factor engine sfcalc_gpu_collapse instead of gemmi sfcalc.
#
# Drop-in replacement for nc2mtz_gemmi.com: same arguments, same avg.mtz/avg.map.
# sfcalc_gpu_collapse takes the P1 supercell PDB + super_mult + primitive SG and
# does the FFT + symmetry collapse + (B>bmax) stripping itself, one GPU per frame.
# Shines most when the B-factor distribution is bimodal (sharp low-B atoms): its
# finer grid does not alias the narrow Gaussians the way gemmi/phenix grids do.
#
set smallSG = ""
set reso = 1
set rate = 2
set B = 10
set super_mult = ( 1 1 1 )
set ncfile = ""
set topfile = xtal.prmtop
set paddedparm = padded.parm7
set orignames = ""
set Bfac_file = ""
set auto_Bfac_file = rmsd2B.pdb

# B factors by LOCATION instead of by atom.  Give a CCP4 map of B(x,y,z) and the
# GPU sfcalc samples it at each atom's position, so nothing has to be matched up
# by atom order and a water is sharp or diffuse according to where it currently
# is.  Bfac_map=auto builds the field here, once for this trajectory, from the
# same per-atom B the old path would have used.  See Bfac_map_runme.com.
set Bfac_map = ""
set Bfac_map_sigma = 0.5
set Bfac_map_grid = 0.5
set Bfac_map_farB = 999

set outtraj = trajectory
set outmap = avg.map
set outfile = avg.mtz

set Bscale = 1
set Boffset = 0
set maxB = 100
set minB = 2

set wrap = 0
set keeptraj = 0
set domaps = 1
set addmaps = 1
set domtzs = 0
set addmtzs = 0

set debug = 0
set sruncpu = ""

set tempfile = /dev/shm/${USER}/temp_$$_traj/

set pdir = `dirname $0`

set CPUs = `grep proc /proc/cpuinfo | wc -l | awk '{print int($1/4)}'`
if( "$CPUs" == "" ) set CPUs = 1

set thishost = `hostname -s`

set startepoch = `msdate.com | awk '{print $7}'`

foreach sourceme ( compute_settings.sourceme xtal_properties.sourceme user_settings.sourceme )
   if(-e $sourceme ) then
      echo "sourcing $sourceme"
      source $sourceme
   endif
end

# read command line
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
      if( $arg == md_mult ) set super_mult = ( $Val )
    else
      # no equal sign
      if("$arg" =~ [pcifrh][1-6]*) set smallSG = `echo $arg | awk '{print toupper($0)}'`
      if("$key" =~ *.nc ) set ncfile = "$Arg"
    endif
    if("$key" == "debug") set debug = "$Val"
end

if( $domtzs ) set domaps = 1
if( $addmaps ) set domaps = 1
if( $addmtzs ) set domtzs = 1

mkdir -p ${tempfile}
if($status) then
  echo "WARNING: reverting to local tempfile"
  set tempfile = /dev/shm/${USER}/nc2mtz_$$_
  mkdir -p ${tempfile}
  if( $status ) then
    set tempfile = ./nc2mtz_$$_
    mkdir -p ${tempfile}
  endif
endif

# use cluster or not?
# cannot migrate hosts because of temp files
set test = `sinfo -h -n $thishost |& egrep -v "drain|n/a" | awk '$2=="up"' | wc -l`
if ( $test ) then
  echo "using slurm"
  if( "$sruncpu" == "" ) set sruncpu = "srun"
  set srun = "$sruncpu -w $thishost"
  set test = `echo $tempfile | awk '{print ( ! /\/dev\/shm/ )}'`
  if( $test ) then
    echo "full cluster"
    set srun = "$sruncpu"
  endif
  if( $debug ) set srun = "$srun -p debug"
else
  set srun = ""
endif

set t = $tempfile

if( "$smallSG" == "" ) set smallSG = P1
set smallSGnum = `awk -v SG=$smallSG '$4 == SG && $1 < 500 {print $1}' $CLIBD/symop.lib | head -1`
set smallSG = `awk -v num=$smallSGnum '$1==num && NF>5{print $4}' ${CLIBD}/symop.lib`
if("$smallSG" == "") then
    set BAD = "bad space group."
    goto exit
endif

set super_mult = `echo $super_mult | awk '{gsub("[,x]"," ");print}'`
if( "$super_mult" == "" ) set super_mult = ( 1 1 1 )
while ( $#super_mult != 3 )
  set super_mult = ( $super_mult $super_mult[$#super_mult] )
end
set nsymops = `awk -v SG=$smallSG '$4 == SG && $1 < 500 {print $2}' $CLIBD/symop.lib | head -1`

cat << EOF
ncfile = $ncfile
smallSG = $smallSG
super_mult = $super_mult
tempfile = $t
EOF


mkdir -p ${t}
rm -f $outtraj > /dev/null
ln -sf ${t} $outtraj


dowrap:
set image = ""
if( $wrap ) set image = "image byatom"
cat << EOF >! ${t}cpptraj.in
$image
strip :WAT,HOH@Y1,EPW
outtraj ${outtraj}/md.pdb pdb multi pdbv3 keepext sg "P 1"
go
EOF
# outtraj trajectory/md.pdb pdb multi pdbv3 keepext sg "P 1" onlyframes ${s}-${f}

again:
echo "cpptraj..."
echo list |\
cpptraj -y $ncfile -p $topfile >&! ${t}cpptraj.log
if($status && -e "$paddedparm" && ! $?RESIZE) then
  set RESIZE
  echo "making new parm file from $paddedparm"
  mv ${t}cpptraj.log ${t}cpptraj_error1.log
  set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");for(i=NF;i>3;--i)if($i+0>0)print $i;exit}' ${t}cpptraj_error1.log | head -n 1`
  if( "$rstatoms" == "" ) then
   set rstatoms = `awk '/Error: Number of atoms in /{gsub("[)(]","");print $(NF-2);exit}' ${t}cpptraj_error1.log`
  endif
  set maxatoms = `echo list | cpptraj -p $paddedparm | awk '$4=="atoms,"{print $3}' | head -n 1`
  set stripmask = `echo $rstatoms $maxatoms | awk '{print $1+1"-"$2}'`
  set topfile = ${t}resized.parm7
  echo "new parmfile: $topfile"
  cpptraj -p $paddedparm << EOF >&! ${t}strip1.log
  parmstrip @$stripmask
  parmwrite out $topfile
EOF
  goto again
endif

set nframes = `awk '/Coordinate processing will occur on/{print $6}' ${t}cpptraj.log`
echo "$nframes frames in $ncfile"
if( $nframes < 10 ) then
  set chunks = 1
else
  set chunks = $CPUs
endif
if( $chunks > ( $nframes / 2 ) ) @ chunks = ( $nframes / 2 )
if( $chunks < 1 ) set chunks = 1

#set t0 = `msdate.com | awk '{print $NF}'`
@ chunksize = ( $nframes / $chunks )
foreach chunk ( `seq 1 $chunks` )
  set s = `echo $chunk $chunksize | awk '{print ($1-1)*$2+1}'`
  set f = `echo $s $chunksize | awk '{print $1+$2-1}'`
  if( $chunk == $chunks ) set f = $nframes
  echo "chunk $chunk is $s - $f"
  cat << EOF >! ${t}cpptraj_${chunk}.in
$image
strip :WAT,HOH@Y1,EPW
outtraj ${outtraj}/md.pdb pdb multi pdbv3 keepext sg "P 1" onlyframes ${s}-${f}
go
EOF
  cat ${t}cpptraj_${chunk}.in | $srun cpptraj -y $ncfile -p $topfile >&! ${t}cpptraj_${chunk}.log &
end
wait
#set dt = `msdate.com $t0 | awk '{print $NF}'`
#echo $chunks $dt | tee -a results.txt


set pdb = ${outtraj}/md.${nframes}.pdb
if( ! -e $pdb ) then
    set BAD = "conversion failed"
    goto exit
endif

set pdbs = `ls -1 ${outtraj}/md.*.pdb`
if( $#pdbs != $nframes ) then
    set BAD = "conversion count mismatch"
    goto exit
endif

egrep "^CRYST|^ATOM|^HETAT" $pdb | head >! ${t}dummy.pdb
echo $super_mult |\
cat - ${t}dummy.pdb |\
awk 'NR==1{na=$1;nb=$2;nc=$3}\
  na<1{na=1} nb<1{nb=1} nc<1{nc=1} \
  /^CRYST1/{print $2/na,$3/nb,$4/nc,$5,$6,$7;exit}' |\
cat >! ${t}cell.txt
set CELL = `cat ${t}cell.txt`
if( $#CELL != 6 ) then
    set BAD = "bad unit cell: $CELL"
    goto exit
endif

pdbset xyzin ${t}dummy.pdb xyzout ${t}cell.pdb << EOF >! ${t}pdbset.log
CELL $CELL
SPACE $smallSG
EOF
if($status) then
    cat ${t}pdbset.log
    set BAD = "pdbset failed"
    goto exit
endif

# --- GPU structure-factor engine (sfcalc_gpu_collapse) ------------------------
# Locate the GPU sfcalc binary (ships with sfcalc_gpu.so + libcufft.so.11 in the
# same dir).  Override with sfcalc_gpu=/path/to/sfcalc_gpu_collapse on the cmdline.
if( ! $?sfcalc_gpu ) set sfcalc_gpu = ""
if( "$sfcalc_gpu" == "" ) then
  foreach d ( ${pdir}/sfcalc_gpu ${pdir} /home/jamesh/projects/fft_symmetry/claude_test )
    if( -x $d/sfcalc_gpu_collapse ) then
      set sfcalc_gpu = $d/sfcalc_gpu_collapse
      break
    endif
  end
endif
if( "$sfcalc_gpu" == "" ) then
  set BAD = "sfcalc_gpu_collapse not found - pass sfcalc_gpu=/path/to/sfcalc_gpu_collapse"
  goto exit
endif
echo "GPU sfcalc: $sfcalc_gpu"
# sfcalc_gpu wants the SUPERCELL cell in the PDB (it derives the primitive cell
# from super_mult itself), so keep a supercell-CRYST1 reference and a csv mult.
egrep "^CRYST1" ${t}dummy.pdb >! ${t}scell.pdb
set super_mult_csv = `echo $super_mult | awk '{print $1","$2","$3}'`
# the frame the trajectory coordinates actually live in, for Bfac_map
set SUPERCELL = `egrep "^CRYST1" ${t}dummy.pdb | head -n 1 | awk '{print $2","$3","$4","$5","$6","$7}'`
# one frame at a time on the GPU (a single srun --gres=gpu:1)
if( "$srun" == "" ) then
  set srungpu = ""
else
  set srungpu = "$srun --gres=gpu:1"
endif

#setenv MEMSIZE `echo $CELL $reso | awk '{print int($1*$2*$3/($NF**3)*100)}'`

touch ${t}Bfac.pdb
if(-e "$Bfac_file") then
   cat  $Bfac_file |\
   awk '/^ATOM|^HETAT/{print $0,"BFACTOR"}' |\
   cat >! ${t}Bfac.pdb
endif

set ns = `ls ${outtraj}/md.*.pdb | awk '{gsub("[^0-9]","");print}' | sort -g`
echo "$#ns md.*.pdb files in ${outtraj}/"

if(-e "$orignames") then
  echo "applying $orignames"
  foreach n ( $ns )
    set pdb = ${outtraj}/md.${n}.pdb
    egrep "^CRYST1" ${t}scell.pdb >! ${t}out.pdb
    awk '/^ATOM|^HETAT/ && ! /EPW|Y1  HOH|Y 1  HOH/{print $0,"ORIG"}' $orignames |\
    cat - $pdb |\
    awk '$NF=="ORIG"{++o;pre[o]=substr($0,1,30);post[o]=substr($0,55,length($0)-55-4);next}\
      ! /^ATOM|^HETAT/{next}\
          {++n}\
          pre[n]==""{print "REMARK WARNING atom",n,"missing from orignames.pdb";\
           pre[n]=substr($0,1,30);post[n]=substr($0,55)}\
          {printf("%s%s%s\n",pre[n],substr($0,31,24),post[n])}' |\
    cat >> ${t}out.pdb
    mv ${t}out.pdb $pdb
  end
endif


if( "$Bfac_file" == "rmsd2B" ) then
  echo "rmsd2B..."
  echo "atomicfluct out ${t}bfac.out bfactor" |\
  cpptraj.OMP -y $ncfile -p $topfile >! ${t}cpptraj_rmsd.log

  cat ${t}bfac.out $pdb |\
  awk 'NF==2{a=int($1+0);B[a]=$2;next}\
    ! /^ATOM|^HETAT/{print;next}\
    {++i;pre=substr($0,1,60);post=substr($0,67);\
     printf("%s%6.2f%s\n",pre,B[i],post)}' |\
  cat >! $auto_Bfac_file
  set Bfac_file = $auto_Bfac_file
endif

if(-e "$Bfac_file") then
   echo "using B factors clipped to [ $minB : $maxB ] from $Bfac_file"
   echo "$minB $maxB $Bscale $Boffset" |\
   cat - $Bfac_file |\
   awk 'NR==1{minB=$1;maxB=$2;Bscale=$3;Boffset=$4;next}\
     ! /^ATOM|^HETAT/{next}\
     {B=substr($0,61,6)*Bscale+Boffset;\
      pre=substr($0,1,60);post=substr($0,67)}\
     B<minB{B=minB} B>maxB{B=maxB}\
     {printf("%s%6.2f%s   BFACTOR\n",pre,B,post)}' |\
   cat >! ${t}Bfac.pdb
endif

# ---- B-factor field ---------------------------------------------------------
set bfacmapopt = ""
if( "$Bfac_map" != "" ) then
  set Bmapexe = ${pdir}/Bfac_map
  if( ! -x "$Bmapexe" ) then
    echo "compiling Bfac_map ..."
    gcc -O3 -o $Bmapexe ${pdir}/Bfac_map.c -lm
    if( $status || ! -x "$Bmapexe" ) then
      set BAD = "cannot compile ${pdir}/Bfac_map.c"
      goto exit
    endif
  endif
  if( "$Bfac_map" == "auto" || "$Bfac_map" == "1" ) then
    # ${t}Bfac.pdb already has minB/maxB/Bscale/Boffset applied, so the field
    # holds final B values and nothing needs to be re-conditioned downstream
    set test = `egrep -c "^ATOM|^HETAT" ${t}Bfac.pdb`
    if( "$test" == "0" ) then
      set BAD = "Bfac_map=auto needs a source of per-atom B: pass Bfac_file=<pdb> or Bfac_file=rmsd2B"
      goto exit
    endif
    # The field's coordinates have to be real, so build it from a trajectory
    # frame with the per-atom B merged onto it rather than from $Bfac_file's own
    # coordinates.  A production Bfac.pdb is a B sidecar: it can carry hundreds
    # of thousands of placeholder atoms parked at the origin (padding to match
    # orignames.pdb), and a spatial field cannot represent a stack of different
    # B values at one point - it can only return their average.
    # This merge is also the LAST place atom order matters, and it now happens
    # once per stage instead of once per frame.
    set refpdb = ${outtraj}/md.${nframes}.pdb
    if( ! -e "$refpdb" ) set refpdb = ${outtraj}/md.$ns[1].pdb
    echo "merging per-atom B onto $refpdb to build the field"
    cat ${t}Bfac.pdb $refpdb |\
    awk -v defB=$B '! /^ATOM|^HETAT/{next}\
      $NF=="BFACTOR"{++i;Bfac[i]=substr($0,61,6)+0;next}\
      {++n;B=defB;\
       if(Bfac[n]=="") ++miss; else B=Bfac[n];\
       printf("%s%6.2f%s\n",substr($0,1,60),B,substr($0,67))}\
      END{if(miss)printf("REMARK %d frame atoms had no B in the sidecar, given B=%g\n",miss,defB)}' |\
    cat >! ${t}Bref.pdb
    set test = `egrep -c "^ATOM|^HETAT" ${t}Bref.pdb`
    if( "$test" == "0" ) then
      set BAD = "could not merge B factors onto $refpdb"
      goto exit
    endif
    echo "building B-factor field from $Bfac_file ($test atoms) on cell $SUPERCELL"
    $Bmapexe build pdb=${t}Bref.pdb outmap=${t}Bfac.map cell=$SUPERCELL \
       sigma=$Bfac_map_sigma grid=$Bfac_map_grid farB=$Bfac_map_farB
    if( $status || ! -e ${t}Bfac.map ) then
      set BAD = "Bfac_map build failed"
      goto exit
    endif
    set Bfac_map = ${t}Bfac.map
  endif
  if( ! -e "$Bfac_map" ) then
    set BAD = "cannot read Bfac_map $Bfac_map"
    goto exit
  endif
  echo "sampling B from $Bfac_map instead of per-atom B factors"
  set bfacmapopt = "bfacmap=$Bfac_map"
endif

if( ! $domaps && ! -e "$Bfac_file" && "$Bfac_map" == "" ) goto cleanup

cat << EOF >! ${t}job.csh
#! /bin/tcsh -f
  set n = "\$1"
  set t = $t
  set B = $B
  set reso = $reso
  set rate = $rate
  set domaps = $domaps
  set pdb = ${outtraj}/md.\${n}.pdb
  egrep "^CRYST1" \${t}scell.pdb >! ${outtraj}/pdb\${n}.pdb
  cat \${t}Bfac.pdb \$pdb |\\
  awk -v B=\$B 'BEGIN{occ=1}\\
     ! /^ATOM|^HETAT/{next}\\
       \$NF=="BFACTOR"{++i;Occ[i]=substr(\$0,55,6);Bfac[i]=substr(\$0,61,6)+0;next}\\
       {++n;mid=substr(\$0, 15, 40);\\
        X =  substr(\$0, 31, 8);\\
        Y =  substr(\$0, 39, 8);\\
        Z =  substr(\$0, 47, 8);\\
       atom=substr(\$0,12,5);gsub(" ","",atom)\\
       Ee = \$NF;}\\
     mid~/\\*/{print "ERROR - corrupt coordinates";exit}\\
     atom=="SE"{Ee=atom;}\\
     Occ[n]!=""{occ=Occ[n]}\\
     Bfac[n]!=""{B=Bfac[n]}\\
     {printf("ATOM  %5d %-2s%40s%6.2f%6.2f%12s%12d\\n",n%100000,Ee,mid,occ,B,Ee,n);}\\
     END{print "END"}' |\\
  cat >> ${outtraj}/pdb\${n}.pdb
  set test = \`tail -n 2 ${outtraj}/pdb\${n}.pdb | awk '/ERROR/{print;exit}'\`
  if( "\$test" != "" ) then
      echo "ERROR: all zero coordinates at \$n"
      exit 9
  endif

  if( ! \$domaps ) exit

  set newmap = \${t}/\${n}.map
  set newmtz = \${t}/\${n}.mtz
  $sfcalc_gpu ${outtraj}/pdb\${n}.pdb sg=$smallSG super_mult=$super_mult_csv \\
     dmin=\$reso rate=\$rate bmax=$maxB $bfacmapopt outmtz=\$newmtz

EOF
chmod a+x ${t}job.csh

echo "applying B factors to ${outtraj} as pdb##.pdb"
if( $domaps ) echo "and rendering maps in $smallSG with cell $CELL"
set mtzs = ""
set maps = ""
# The GPU does one frame at a time, so run sequentially (a single srun --gres=gpu:1)
# rather than the CPU path's parallel fan-out.  ~a few seconds per frame.
foreach n ( $ns )
  set pdb = ${outtraj}/md.${n}.pdb
  set newmap = ${t}/${n}.map
  set newmtz = ${t}/${n}.mtz
  $srungpu ${t}job.csh $n >&! ${t}/job.${n}.log
  # sfcalc_gpu_collapse ignores outmap= (not implemented), so build the per-frame
  # map from its now-valid FC/PHIC mtz here on the CPU host (gemmi is on PATH here,
  # not necessarily on the GPU node).  -s $rate matches the sfcalc grid so every
  # frame shares one grid and addup_maps can sum them.
  if( -e $newmtz ) then
    gemmi sf2map -f FC -p PHIC -s $rate $newmtz $newmap >&! ${t}/sf2map.${n}.log
  endif
  set maps = ( $maps $newmap )
  set mtzs = ( $mtzs $newmtz )
end

set test = `tail -n 2 ${outtraj}/pdb*.pdb | awk '/^==/{f=$2} /ERROR/{print $0,f;exit}'`
if( "$test" != "" ) then
  if( $wrap ) then
      set BAD = "all zero coordinates at $n"
      goto exit
  else
     echo "trying with wrap"
     set wrap = 1
     goto dowrap
  endif
endif
if( $?BAD ) goto exit

if( ! $addmaps && $addmtzs ) goto addmtzs
if( ! $addmaps ) goto cleanup

rm -f sum.map >& /dev/null
${pdir}/addup_maps_runme.com $maps tempfile=${t}/addmaps/  outfile=${t}sum.map

rm -f ${t}sum.mtz >& /dev/null
gemmi map2sf -v --dmin=$reso ${t}sum.map ${t}sum.mtz Fsum PHIsum

# sfcalc_gpu_collapse already collapsed each frame to the primitive ASU (symmetry
# and supercell normalization handled internally), so the only averaging factor
# is 1/nframes.  Downstream scaleB_search absorbs residual scale; validate the
# absolute scale against nc2mtz_gemmi on a test frame before relying on it.
echo "scaling factor: 1 / $#ns frames"
set scale = `echo $#ns | awk '{print 1/$1}'`
echo "scaling by $scale for $outmap"
echo scale factor $scale |\
mapmask mapin1 ${t}sum.map mapout $outmap >! ${t}scaledown.log

rm -f ${outfile}
gemmi map2sf -v --dmin=$reso $outmap ${outfile} FCavg PHICavg


if( ! $addmtzs ) goto cleanup

addmtzs:
${pdir}/addup_mtzs_diffuse.com $mtzs

cleanup:
if( ! $keeptraj ) then
  echo "cleaning up..."
  if(! $debug ) rm -rf ${t}
  rm -f ${outtraj}
endif

exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9
endif

msdate.com $startepoch
if( $domaps ) ls -l ${outfile}
echo ""

exit


#################################
# notes and tests
#

foreach minB ( 1 5 10 20 )
foreach maxB ( 10 20 50 100 200 500 1000 )
 if( $maxB < $minB ) break

 set nc2log = nc2mtz_${minB}_${maxB}.log

 echo "B: $minB $maxB"
 nc2mtz_gemmi_new.com amber_${itr}.nc P212121 \
   Bfac_file=rmsd2B minB=$minB maxB=$maxB |& tee $nc2log | tail -n 1

 diff.com reference.mtz avg.mtz |& tee diff_${minB}_${maxB}.log | grep correct

cad hklin1 Fdiff.mtz \
    hklin2 avg.mtz \
    hklin3 reference.mtz \
    hklout fftme.mtz << EOF >> $nc2log
labin file 1 E1=Fref E2=Ftest
labin file 2 E1=PHICavg
labou file 2 E1=PHItest
labin file 3 E1=PHIref
resolution over_all $reso
EOF

rm -f cootme.mtz
sftools << EOF >> $nc2log
read fftme.mtz
map correl Fref PHIref Ftest PHItest
calc ( COL DELFWT PHDELWT ) = ( COL Fref PHIref ) ( COL Ftest PHItest ) -
calc ( COL FWT PHWT ) = ( COL Fref PHIref )
write cootme.mtz col FWT PHWT DELFWT PHDELWT
quit
EOF
fft hklin cootme.mtz mapout fofc.map << EOF >> $nc2log
labin F1=DELFWT PHI=PHDELWT
EOF

scaleB_search_diffmap_runme.com avg.mtz |& tee -a $nc2log | grep best
fft hklin cootme-grid.mtz mapout fofc-grid.map << EOF >> $nc2log
labin F1=DELFWT PHI=PHDELWT
EOF

 mv avg.mtz avg_${minB}-${maxB}.mtz
 mv avg.map avg_${minB}-${maxB}.map
 mv cootme.mtz cootme_${minB}-${maxB}.mtz
 mv fofc.map fofc_${minB}-${maxB}.map
 mv fofc-grid.map fofc-grid_${minB}-${maxB}.map

end
end

cad hklin1 avg_1-10.mtz hklin2 avg_1-1000.mtz hklout unscaled.mtz << EOF
labin file 1 E1=FCavg E2=PHICavg
labin file 2 E1=FCavg E2=PHICavg
labou file 1 E1=F1 E2=PHI1
labou file 2 E1=F2 E2=PHI2
EOF
fft hklin unscaled.mtz mapout fofc-mtz-unscaled.map << EOF
labin F1=F1 PHI=PHI1 F2=F2 PHI2=PHI2
EOF
echo scale factor -1 |\
mapmask mapin1 avg_1-1000.map mapout neg.map
echo maps add |\
mapmask mapin1 avg_1-10.map mapin2 neg.map mapout fofc-unscaled.map

awk '/^ATOM|^HETAT/ && ! /EPW/{print $0,"ORIG"}' orignames.pdb |\
cat - rmsd2B.pdb |\
awk '$NF=="ORIG"{++o;pre[o]=substr($0,1,30);post[o]=substr($0,55,length($0)-55-4);next}\
  ! /^ATOM|^HETAT/{next}\
      {++n}\
      pre[n]==""{print "REMARK WARNING atom",n,"missing from orignames.pdb";\
       pre[n]=substr($0,1,30);post[n]=substr($0,55)}\
      {printf("%s%s\n",pre[n],substr($0,31))}' |\
cat >! rmsd2B_orignames.pdb



fft hklin avg.mtz hklout this.map << EOF
labin
EOF



