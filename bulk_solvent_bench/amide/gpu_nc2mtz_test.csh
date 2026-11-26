#! /bin/tcsh -f
# One trajectory through nc2mtz_gemmi.com and nc2mtz_gpu.com (2026-11-25): ff19_opc_rc_dw_w's itr-199 frames
# (aligned.nc, 25 frames; the run's own avg_199.mtz came from them), with the arguments optimize_weights gives the
# 1aho runs (uniform B 5, reso 0.965, rate 1.5, P212121, 1x1x1).  Scripts from the bespoke_amber_restraints checkout.
# usage (on octamus1): gpu_nc2mtz_test.csh [engines]      default engines: gemmi gpu
source /home2/jamesh/projects/amber/bulk_solvent/programs.sourceme
set R = /scratch/jamesh/projects/amber/bespoke_amber_restraints
set path = ( $R $path )
set src = /scratch/jamesh/projects/amber/1aho_hurry_test/ff19_opc_rc_dw_w
set T = /scratch/jamesh/projects/amber/1aho_hurry_test/gpu_nc2mtz_test
set engines = ( $* )
if( $#engines == 0 ) set engines = ( gemmi gpu )
mkdir -p $T ; cd $T
foreach f ( aligned.nc xtal.prmtop padded.parm7 xtal_properties.sourceme avg_199.mtz refme.mtz reference.mtz )
  if( ! -e $f && -e $src/$f ) cp -p $src/$f .
end
# the site settings, plus the GPU queue for the structure factors
cat << EOF >! compute_settings.sourceme
source /home2/jamesh/projects/amber/compute_settings.sourceme
set srungpu = "srun --partition=gpu --gres=gpu:1"
EOF
set args = ( P212121 super_mult=1,1,1 aligned.nc reso=0.965 B=5 Bfac_file=none minB=0 maxB=999.99 \
             keeptraj=1 wrap=1 rate=1.5 addmtzs=1 )
foreach eng ( $engines )
  set s0 = `msdate.com | awk '{print $7}'`
  nc2mtz_${eng}.com $args tempfile=/scratch/jamesh/nc2mtz_test_${eng}/ \
     outfile=avg_${eng}.mtz outmap=avg_${eng}.map >&! nc2mtz_${eng}.log
  set st = $status
  set dt = `msdate.com $s0 | awk '{print $NF}'`
  set host = `hostname -s`
  echo "${eng}: exit $st, wall $dt s on $host"
  grep -E "ERROR|WARNING|GPU|scaling by" nc2mtz_${eng}.log | head -6
end
