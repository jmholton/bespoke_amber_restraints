#! /bin/tcsh -f
# usage: bfac_iter.csh run itr outdir   (on octamus1)
# solvent density of one iteration rendered by nc2mtz_gemmi.com with the run's OWN per-atom B, as the pipeline rendered
# it: amber_<itr>.nc was rendered in iteration itr+1 with the Bfac.pdb saved as premod_Bfac_<itr+1>.pdb, so that file
# is the B sidecar, matched to the frames by atom order.  The nc is read with a copy of padded.parm7 stripped to its own
# atom count, then cut to the solvent (WAT/HOH, NH4, ACY) with the sidecar cut the same way.  B is clipped to the
# pipeline's [1 : 999.99] (optimize_weights).  Same frames as the hurry maps (every 5th of 25, whatever the frame stride in sets.txt) -> outdir/<run>_<itr>/avg.map
set H = /scratch/jamesh/projects/amber/1aho_hurry_test
set N = /home2/jamesh/projects/amber/bulk_solvent/bespoke_amber_restraints
set run = $1
set itr = $2
set outdir = $3
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
set nc = $H/$run/amber_${itr}.nc
set Bfile = $H/$run/premod_Bfac_`expr $itr + 1`.pdb
if( ! -e $nc || ! -e $Bfile ) then
  echo "${run} ${itr}: missing $nc or $Bfile"
  exit 9
endif
set w = $outdir/${run}_${itr}
rm -rf $w
mkdir -p $w
cd $w
set n = `/programs/phenix-2.2-6143/bin/phenix.python $N/../water_models/hurry/natoms.py $nc | awk '{print $2}'`
cat << EOF >! strip.in
parm $H/$run/padded.parm7
parmstrip @`expr $n + 1`-9999999
trajin $nc 5 last 5
strip !(:WAT,HOH,NH4,ACY) parmout solvent.parm7
trajout solvent.nc
trajout check.pdb pdb onlyframes 1
go
EOF
cpptraj < strip.in >&! strip.log
if( $status || ! -e solvent.nc ) then
  echo "${run} ${itr}: cpptraj failed (see $w/strip.log)"
  exit 9
endif
# the sidecar's solvent atoms, in order (extra points kept: nc2mtz drops them from frames and sidecar alike)
awk '/^ATOM|^HETAT/{r=substr($0,18,3)} r~/^(WAT|HOH|NH4|ACY)$/' $Bfile >! Bsolvent.pdb
# atom order check: every solvent atom of the frame has the sidecar atom of the same name and residue name (extra
# points skipped, as nc2mtz skips them)
set bad = `awk '/^ATOM|^HETAT/ && ! /EPW|Y1  HOH/{nm=substr($0,13,4);rn=substr($0,18,3);gsub(" ","",nm);sub("WAT","HOH",rn);\
             if(FILENAME==ARGV[1]){++i;S[i]=rn" "nm}else{++j;if(S[j]!=rn" "nm)++bad}}\
             END{print bad+0,j,i}' Bsolvent.pdb check.pdb`
if( $bad[1] != 0 || $bad[2] == 0 ) then
  echo "${run} ${itr}: sidecar order mismatch: $bad[1] of $bad[2] frame atoms ($bad[3] in sidecar)"
  exit 9
endif
$N/nc2mtz_gemmi.com solvent.nc topfile=solvent.parm7 paddedparm=none P212121 reso=0.965 rate=2 \
  Bfac_file=Bsolvent.pdb minB=1 maxB=999.99 outmap=avg.map outfile=avg.mtz \
  tempfile=/dev/shm/${USER}/bfac_${run}_${itr}_ >&! nc2mtz.log
if( $status || ! -e avg.map ) then
  echo "${run} ${itr}: nc2mtz_gemmi.com failed (see $w/nc2mtz.log)"
  exit 9
endif
set nf = `awk '/frames in/{print $1;exit}' nc2mtz.log`
echo $nf >! nframes
rm -f solvent.nc avg.mtz check.pdb
echo "${run} ${itr}: $nf frames, $bad[2] solvent atoms, sidecar ok"
