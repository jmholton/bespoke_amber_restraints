#!/bin/bash
# Salt series (2026-11-27; user: "also want to look at different salt concentrations once protonation is more
# settled", then "yes, also try those next tests"): ammonium acetate at 0.17, 0.34 and 1.36 M beside the 0.68 M runs,
# on the best-guess His state DD (HID54 / HID64) and the MD-preferred alternative DE (HID54 / HIE64), both amide.
# Each concentration: one salt placement (salt_base.csh, = salt068_base.csh at that conc, cell charge +4), then per
# arm exactly as setup_s68_builds.sh: Cys12/Cys63 alt-B graft into copies 3+4, +4 counter-acetates for the amide's +8
# (add_counter_acy.py), protonation and tleap stub from build_s68_<arm>.
# Builds build_s<tag>_<arm>; runs ff19_opc_s<tag>_<arm> (launch_s68_controls.sh takes name:build, see below).
# usage: setup_salt_series.sh ["0.17:017 0.34:034 1.36:136"] ["DD DE"]
set -e
H=/scratch/jamesh/projects/amber/1aho_hurry_test; A=/scratch/jamesh/projects/amber/amide
P=/scratch/jamesh/programs/amber22/miniconda/bin/python3
for ct in ${1:-0.17:017 0.34:034 1.36:136}; do
  conc=${ct%%:*}; tag=${ct#*:}; base=$H/salt${tag}_base
  if [ ! -e $base/amberme.pdb ]; then
    /bin/tcsh -f $A/salt_base.csh $conc $tag
    [ -e $base/amberme.pdb ] || { echo "salt base $conc FAILED"; exit 9; }
  fi
  echo "salt $conc M: NH4 $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="NH4" && substr($0,13,4)==" N  "' $base/amberme.pdb | wc -l) ACY $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="ACY" && substr($0,13,4)==" C  "' $base/amberme.pdb | wc -l) at cell charge +4"
  for v in ${2:-DD DE}; do
    d=$H/build_s${tag}_$v; s=$H/build_s68_$v
    [ -e $d/Prod.rst7 ] && { echo "$d already built"; continue; }
    mkdir -p $d
    for f in refined.pdb protonation.txt tleap_stub.in netfrc0_stub.in omegafix_stub.in NH4.cif NH4.frcmod NH4.mol2 \
             NH4.pdb anion.pdb cation.pdb HIS.pdb HIS_settings_asu.txt restraint_suffix.in Invert.in refme.mtz \
             reference.mtz ACY.cif ACY.frcmod ACY.mol2 ACY.pdb xtal_properties.sourceme; do
      [ -e $s/$f ] && cp -p $s/$f $d/
    done
    ln -sf ../compute_settings.sourceme $d/compute_settings.sourceme
    cp -p $base/amberme.pdb $d/salt${tag}_amberme.pdb
    # graft: alt-B Cys12/Cys63 of copies 3 and 4 (as setup_s68_builds.sh / run_local.csh)
    awk 'FNR==NR{if(/^ATOM|^HETAT/ && substr($0,17,1)=="B" && substr($0,18,3)=="CYS"){r=substr($0,23,4)+0;\
          if(r==140||r==191||r==204||r==255) xyz[r" "substr($0,13,4)]=substr($0,31,24)}; next}\
         /^ATOM|^HETAT/ && substr($0,18,3)=="CYS"{k=(substr($0,23,4)+0)" "substr($0,13,4); if(k in xyz){\
          $0=substr($0,1,30) xyz[k] "  0.45" substr($0,61); ++n}} {print}\
         END{print n" atoms grafted" > "/dev/stderr"}' $H/centroids/centroids_in_density.pdb $d/salt${tag}_amberme.pdb \
         > $d/amberme.pdb
    $P $A/add_counter_acy.py $d/amberme.pdb 4
    echo "s${tag}_$v: $(tr '\n' ' ' < $d/protonation.txt)| NH4 $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="NH4" && substr($0,13,4)==" N  "' $d/amberme.pdb | wc -l) ACY $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="ACY" && substr($0,13,4)==" C  "' $d/amberme.pdb | wc -l) HOH $(grep -c -E ' (HOH|WAT) ' $d/amberme.pdb)"
  done
done
