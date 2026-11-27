#!/bin/bash
# More His protonation arms at 0.68 M salt (2026-11-26; user: "queue up more His protonation tests", after the worst-
# restraint analysis: HIP54 puts His54's NE2/CE1 edge among the worst restraints, HID64 strains His64's ring edge, most
# with HIP/HID54).  Completes His54 {HIE, HID, HIP, HIP/HID mix} x His64 {HIE, HID}, all amide:
#   PE  HIP54, HIE64              +12  amberme (11 NH4 / 23 ACY) from build_s68_PD
#   DE  HID54, HIE64              +8   from build_s68_DD (11 / 19)
#   ME  HIP54 copies 1+3, HID54 2+4, HIE64   +10  from build_s68_MD (11 / 21)
#   (MD = the same mix with HID64 was built by setup_s68_builds.sh)
# HIE64 is written as EE's protonation.txt does (blank name: tleap's default HIE -> the CHIE amide of tleap_stub.in).
# usage: setup_s68_more.sh   (then launch_s68_controls.sh PE:build_s68_PE ...)
set -e
H=/scratch/jamesh/projects/amber/1aho_hurry_test
declare -A SRC=( [PE]=build_s68_PD [DE]=build_s68_DD [ME]=build_s68_MD )
declare -A P54=( [PE]="HIP HIP HIP HIP" [DE]="HID HID HID HID" [ME]="HIP HID HIP HID" )
for v in PE DE ME; do
  d=$H/build_s68_$v; s=$H/${SRC[$v]}
  [ -e $d/Prod.rst7 ] && { echo "$d already built"; continue; }
  mkdir -p $d
  for f in refined.pdb tleap_stub.in netfrc0_stub.in omegafix_stub.in NH4.cif NH4.frcmod NH4.mol2 NH4.pdb anion.pdb \
           cation.pdb HIS.pdb HIS_settings_asu.txt restraint_suffix.in Invert.in refme.mtz reference.mtz ACY.cif \
           ACY.frcmod ACY.mol2 ACY.pdb xtal_properties.sourceme salt068_amberme.pdb amberme.pdb; do
    [ -e $s/$f ] && cp -p $s/$f $d/
  done
  ln -sf ../compute_settings.sourceme $d/compute_settings.sourceme
  set -- ${P54[$v]}
  printf "%s A54\n A64\n%s A118\n A128\n%s A182\n A192\n%s A246\n A256\n" $1 $2 $3 $4 > $d/protonation.txt
  echo "s68_$v: $(tr '\n' ' ' < $d/protonation.txt)| ACY $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="ACY" && substr($0,13,4)==" C  "' $d/amberme.pdb | wc -l) (from ${SRC[$v]})"
done
