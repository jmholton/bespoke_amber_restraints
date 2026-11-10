#! /bin/tcsh -f
# re-make j5w8760/solvent_RestAll.map as 1aho_test/allmapme.com did on 2019-08-01:
# all 500 RestAll.nc frames, translation-only fit on protein CA, solvent (HOH/WAT/AMM/ACT) atoms, B=5, GRID 192 144 256
cd /home2/jamesh/projects/amber/bulk_solvent/water_models/lineage/md2map_rerun
source /home2/jamesh/projects/amber/bulk_solvent/compute_settings.sourceme >& /dev/null
rm -rf trajectory ; mkdir -p trajectory
cpptraj << EOT >&! cpptraj.log
parm /home2/jamesh/projects/amber/1aho_test/j5w8760/start.top
reference /home2/jamesh/projects/amber/1aho_test/j5w8760/start.crd
trajin /home2/jamesh/projects/amber/1aho_test/j5w8760/RestAll.nc
rms reference '!:ACT,AMM,WAT&@CA' before norotate out rmsd_RestAll.txt time 0.2
rms reference '!:ACT,AMM,WAT&@CA' after norotate out rmsd_RestAll.txt time 0.2
strip !:WAT,HOH,AMM,ACT
outtraj trajectory/md.pdb pdb multi pdbv3 keepext sg "P 1"
go
EOT
echo "cpptraj status $status"
echo "2 2 3" >! MD_mult.txt
echo "GRID 192 144 256" >! grid.txt
ls -1 trajectory/md*.pdb | sort -k1.15g >! pdbfiles.txt
echo "md2map solvent on `cat pdbfiles.txt | wc -l` pdb files ..."
./md2map_multi.com P212121 pdbfiles.txt -solventonly >&! md2map_solvent.log
ls -l avg.map
echo done
