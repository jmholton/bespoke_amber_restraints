#! /bin/tcsh -f
# solvent map from production segments first..N (1 ns each) of one water-model run, as section 11 of
# bulk_solvent_runme.com does
# usage: submap.csh rundir N [first]  ->  rundir/maps_${tag}/bulk_solvent_${tag}.mtz, tag = Nns or first-Nns
set top = `cd $1 ; pwd`
set n = $2
set first = 1
if( $#argv > 2 ) set first = $3
set tag = ${n}ns
if( $first != 1 ) set tag = ${first}-${n}ns
set path = ( /home2/jamesh/projects/amber/bulk_solvent/water_models/repo_3140b19 $path )
cd $top
source compute_settings.sourceme
mkdir -p maps_${tag}
cd maps_${tag}
foreach f ( xtal_properties.sourceme user_settings.sourceme compute_settings.sourceme )
  ln -sf ../$f .
end
echo "parm ../hydrated.prmtop" >! strip.in
foreach i ( `seq $first $n` )
  if(! -e ../prod_${i}.rst7) then
    echo "ERROR: prod_${i} not finished"
    exit 9
  endif
  echo "trajin ../prod_${i}.nc" >> strip.in
end
cat << EOT >> strip.in
strip !(:WAT,HOH,NH4,ACY) parmout solvent.parm7
trajout solvent.nc
go
EOT
rm -f solvent.nc solvent.parm7
cpptraj < strip.in >&! strip.log
nc2mtz_gemmi.com solvent.nc topfile=solvent.parm7 paddedparm=none reso=0.965 B=10 \
  outmap=solvent_avg.map outfile=solvent_avg.mtz >&! nc2mtz.log
cad hklin1 solvent_avg.mtz hklout bulk_solvent_${tag}.mtz << EOT > /dev/null
labin file 1 E1=FCavg E2=PHICavg
labou file 1 E1=FPART E2=PHIPART
EOT
set nf = `awk '/frames in solvent.nc/{print $1}' nc2mtz.log`
ls -l bulk_solvent_${tag}.mtz
echo "${tag}: $nf frames"
