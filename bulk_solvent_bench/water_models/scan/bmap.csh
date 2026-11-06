#! /bin/tcsh -f
# usage: bmap.csh run variant [nc2mtz args]   (on octamus1: cpptraj)
# re-renders scan/<run>/solvent.nc (the solvent-only production trajectory bulk_solvent_runme kept) with other water
# B factors, through nc2mtz_gemmi.com from the bulk_solvent clone of bespoke_amber_restraints (has the water
# extra-point B fix) -> scan/<run>/bmaps/map_<variant>.mtz (FPART PHIPART, as bulk_solvent.mtz)
set S = /home2/jamesh/projects/amber/bulk_solvent/water_models/scan
set R = /home2/jamesh/projects/amber/bulk_solvent/bespoke_amber_restraints
set run = $1
set variant = $2
set args = ( $argv[3-] )
cd $S/$run
source compute_settings.sourceme >& /dev/null
set reso = `awk '/^set reso/{print $NF}' xtal_properties.sourceme`
if( "$reso" == "" ) set reso = 0.965
mkdir -p bmaps/$variant
cd bmaps/$variant
foreach f ( compute_settings xtal_properties user_settings )
  ln -sf ../../${f}.sourceme .
end
ln -sf ../../solvent.nc ../../solvent.parm7 .
$R/nc2mtz_gemmi.com solvent.nc topfile=solvent.parm7 paddedparm=none reso=$reso $args \
  outmap=solvent_avg.map outfile=solvent_avg.mtz >&! nc2mtz.log
cad hklin1 solvent_avg.mtz hklout ../map_${variant}.mtz << EOF > /dev/null
labin file 1 E1=FCavg E2=PHICavg
labou file 1 E1=FPART E2=PHIPART
EOF
if( -e ../map_${variant}.mtz ) then
  echo "${run} ${variant}: wrote bmaps/map_${variant}.mtz  ( $args )"
else
  echo "${run} ${variant}: FAILED (see bmaps/$variant/nc2mtz.log)"
endif
