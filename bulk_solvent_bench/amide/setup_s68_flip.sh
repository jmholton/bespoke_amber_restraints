#!/bin/bash
# His54 ring-flip positive control (2026-11-25; user: "voltron fairly idle. do more controls", after the refinement
# flip test put His54's ring in the M orientation, chi2 -79 to -119, and found best_so_far's chi2 > 0 altlocs
# mislabelled).  Our builds already have His54 in M (chi2 -110, and it never flips in the MD).  build_s68_PDp = PD
# (HIP54, amide HID64, 0.68 M salt) with His54's ring relabelled to the P orientation in all 4 copies: ND1<->CD2 and
# CE1<->NE2 trade coordinates (names and order kept), so N and C trade places on the same ring.  The restraint
# picker already treats ND1/CD2 and CE1/NE2 as one atom (centroids_nearby "twirl"), so restraints are unchanged; only
# the chemistry (which ring atoms are N and carry the N-H) differs.  Known-wrong labels: tests whether the MD
# measures can see a ring N/C label at all.     usage: setup_s68_flip.sh   (then build_s68.csh / launch)
set -e
H=/scratch/jamesh/projects/amber/1aho_hurry_test; s=$H/build_s68_PD; d=$H/build_s68_PDp
[ -e $d/Prod.rst7 ] && { echo "$d already built"; exit 0; }
mkdir -p $d
for f in refined.pdb protonation.txt tleap_stub.in netfrc0_stub.in omegafix_stub.in NH4.cif NH4.frcmod NH4.mol2 \
         NH4.pdb anion.pdb cation.pdb HIS.pdb HIS_settings_asu.txt restraint_suffix.in Invert.in refme.mtz \
         reference.mtz ACY.cif ACY.frcmod ACY.mol2 ACY.pdb xtal_properties.sourceme salt068_amberme.pdb; do
  [ -e $s/$f ] && cp -p $s/$f $d/
done
ln -sf ../compute_settings.sourceme $d/compute_settings.sourceme
# swap coordinates (columns 31-66: xyz, occupancy, B) between ND1<->CD2 and CE1<->NE2 of every His54 copy; atom
# names and order unchanged
awk 'BEGIN{sw["ND1"]="CD2";sw["CD2"]="ND1";sw["CE1"]="NE2";sw["NE2"]="CE1"}
     function key(l){return substr(l,22,5)}
     FNR==NR{if(/^ATOM|^HETAT/ && substr($0,18,3)=="HIS" && ((substr($0,23,4)+0)-54)%64==0) xyz[key($0)" "substr($0,14,3)]=substr($0,31,36); next}
     /^ATOM|^HETAT/ && substr($0,18,3)=="HIS" && ((substr($0,23,4)+0)-54)%64==0 && (substr($0,14,3) in sw){
       $0=substr($0,1,30) xyz[key($0)" "sw[substr($0,14,3)]] substr($0,67); ++n}
     {print} END{print n" His54 ring atoms moved to their flip partner" > "/dev/stderr"}' $s/amberme.pdb $s/amberme.pdb > $d/amberme.pdb
[ $(wc -l < $d/amberme.pdb) -eq $(wc -l < $s/amberme.pdb) ]
/scratch/jamesh/programs/amber22/miniconda/bin/python3 /scratch/jamesh/projects/amber/amide/his_chi.py $s/amberme.pdb $d/amberme.pdb res=54
