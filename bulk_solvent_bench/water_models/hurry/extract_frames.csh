#! /bin/tcsh -f
# usage: extract_frames.csh run first last step frame_stride outdir [align]   (on octamus1: cpptraj)
# solvent-only frames (WAT/HOH, NH4, ACY; water extra points dropped) from 1aho_hurry_test/<run>/amber_<itr>.nc for
# itr = first, first+step, ... last, every frame_stride-th frame of each 25-frame iteration.  The atom count changes
# between iterations (the barometer adds and drops waters), so each nc is read with a copy of padded.parm7 stripped
# to its own atom count (padded.parm7 itself is never touched) -> outdir/<run>_<itr>.<k>.pdb
# align: first translate each frame onto the run's alignment reference the way optimize_weights_runme.com does before
# it makes its mtz (rmsd norotate on align_ref.crd, read with the run's xtal.prmtop, over that iteration's mask from
# align_<itr>.log).  The translations go to outdir/<run>_<itr>_tv.dat and are checked against the ones in the log.
set H = /scratch/jamesh/projects/amber/1aho_hurry_test
set run = $1
set first = $2
set last = $3
set step = $4
set fstride = $5
set outdir = $6
set align = "$7"
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
mkdir -p $outdir
set ncs = ()
foreach i ( `seq $first $step $last` )
  if( -e $H/$run/amber_${i}.nc ) set ncs = ( $ncs $H/$run/amber_${i}.nc )
end
# atom count of every nc at once
/programs/phenix-2.2-6143/bin/phenix.python `dirname $0`/natoms.py $ncs >! $outdir/natoms_$run.txt
if( `cat $outdir/natoms_$run.txt | wc -l` != $#ncs ) then
  echo "could not read atom counts for $run"
  exit 9
endif
foreach line ( "`cat $outdir/natoms_$run.txt`" )
  set nc = `echo $line | awk '{print $1}'`
  set n = `echo $line | awk '{print $2}'`
  set itr = `echo $nc:t:r | awk -F_ '{print $NF}'`
  cat << EOF >! $outdir/${run}_${itr}.in
parm $H/$run/padded.parm7
parmstrip @`expr $n + 1`-9999999
EOF
  if( "$align" == align ) then
    set alog = $H/$run/align_${itr}.log
    set mask = ""
    if( -e $alog ) then
      set mask = `grep -m1 -o 'rmsd rmsd reference norotate @[0-9,-]*' $alog | awk '{print $NF}'`
    endif
    if( "$mask" == "" ) then
      echo "${run} itr ${itr}: no alignment mask in $alog"
      exit 9
    endif
    cat << EOF >> $outdir/${run}_${itr}.in
parm $H/$run/xtal.prmtop
reference $H/$run/align_ref.crd parmindex 1 [aref]
trajin $nc $fstride last $fstride
rmsd R $mask ref [aref] norotate savevectors combined vecsout $outdir/${run}_${itr}_tv.dat
EOF
  else
    echo "trajin $nc $fstride last $fstride" >> $outdir/${run}_${itr}.in
  endif
  cat << EOF >> $outdir/${run}_${itr}.in
strip !(:WAT,HOH,NH4,ACY)
strip :WAT,HOH@Y1,EPW
trajout $outdir/${run}_${itr}.pdb pdb multi pdbv3 keepext sg "P 1"
go
EOF
  cpptraj < $outdir/${run}_${itr}.in >&! $outdir/${run}_${itr}.log
  if( $status ) then
    echo "cpptraj failed on $nc (see $outdir/${run}_${itr}.log)"
    exit 9
  endif
  if( "$align" == align ) then
    # largest difference from the translations the run itself applied to these frames: the TV block that follows
    # the log's "Reading 'amber_<itr>.nc'" (frame k of the extract is frame k*frame_stride of the nc)
    set d = `awk -v f=$fstride -v nc="'amber_${itr}.nc'" 'FNR==1{++fi} fi==1&&$1~/^[0-9]+$/{t[$1*f]=$2" "$3" "$4;++nt} fi==2&&index($0,nc){on=1} fi==2&&p&&$1!~/^[0-9]+$/{p=on=0} fi==2&&on&&/^#Frame/&&/TV/{p=1;next} fi==2&&p&&($1 in t){split(t[$1],a," ");for(k=1;k<=3;++k){x=$(k+1)-a[k];if(x<0)x=-x;if(x>m)m=x};++n} END{printf "%d %d %.4f\n",nt,n,m}' $outdir/${run}_${itr}_tv.dat $alog`
    if( $d[1] == 0 || $d[2] != $d[1] || `echo $d[3] | awk '{print ($1>0.001)}'` ) then
      echo "${run} itr ${itr}: translations differ from $alog ($d[2] of $d[1] frames compared, max diff $d[3] A)"
      exit 9
    endif
  endif
  rm -f $outdir/${run}_${itr}.in $outdir/${run}_${itr}.log
end
set nf = `ls $outdir | grep -cE "^${run}_[0-9]+\.[0-9]+\.pdb"`
echo "${run}: $nf frames from $#ncs iterations"
if( "$align" == align ) then
  cat $outdir/${run}_[0-9]*_tv.dat | awk '$1~/^[0-9]+$/{s+=$2^2+$3^2+$4^2;++n} END{printf "%s: aligned as in the run (translations match its logs); rms translation %.3f A over %d frames\n","'$run'",sqrt(s/n),n}'
endif
