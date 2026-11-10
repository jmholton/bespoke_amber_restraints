#! /bin/tcsh -f
# usage: window_map.csh first last [max_frames [residues]]   (on octamus1: cpptraj)
# solvent map of w1_100ns production segments first..last (1 ns, 50 frames each), at most max_frames (500) frames,
# rendered as bulk_solvent_runme.com section 11 does -> w1_100ns/windows/map_<first>_<last>.mtz
# residues (e.g. WAT,HOH or NH4,ACY): only those, -> w1_100ns/windows/<residues>/map_<first>_<last>.mtz
set D = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan/w1_100ns
set first = $1
set last = $2
set max_frames = 500
if( $#argv > 2 ) set max_frames = $3
set residues = WAT,HOH,NH4,ACY
set out = ..
if( $#argv > 3 ) then
  set residues = $4
  set out = ../`echo $residues | sed 's/,/_/g'`
endif
cd $D
source compute_settings.sourceme >& /dev/null
source xtal_properties.sourceme
set path = ( /home/jamesh/projects/git/bespoke_amber_restraints $path )
set tag = ${first}_${last}
set wdir = w$tag
if( $out != .. ) set wdir = w${tag}_$out:t
mkdir -p windows/$wdir
cd windows/$wdir
mkdir -p $out
# nc2mtz_gemmi.com takes the space group (and the rest) from these, as in the run directory
foreach f ( compute_settings xtal_properties user_settings )
  ln -sf ../../${f}.sourceme .
end
set stride = `echo $first $last $max_frames | awk '{n=($2-$1+1)*50;s=int((n+$3-1)/$3);print (s<1?1:s)}'`
echo "parm ../../xtal.prmtop" >! strip.in
foreach i ( `seq $first $last` )
  echo "trajin ../../prod_${i}.nc 1 last $stride" >> strip.in
end
cat << EOF >> strip.in
strip !(:$residues) parmout solvent.parm7
trajout solvent.nc
go
EOF
cpptraj < strip.in >&! strip.log
nc2mtz_gemmi.com solvent.nc topfile=solvent.parm7 paddedparm=none reso=$reso B=10 \
  outmap=solvent_avg.map outfile=solvent_avg.mtz >&! nc2mtz.log
cad hklin1 solvent_avg.mtz hklout $out/map_${tag}.mtz << EOF > /dev/null
labin file 1 E1=FCavg E2=PHICavg
labou file 1 E1=FPART E2=PHIPART
EOF
rm -f solvent.nc
if( -e $out/map_${tag}.mtz ) then
  echo "$out:t/map_${tag}.mtz  every ${stride}th frame"
else
  echo "$out:t/map_${tag}.mtz FAILED"
endif
