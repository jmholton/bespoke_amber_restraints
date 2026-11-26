#!/bin/bash
# His protonation tests (2026-11-25; user: "queue the His protonation tests").  Five builds from build_2s2r_amide's
# inputs (His64 C-terminal amide, 2 strained + 2 relaxed Cys12, ACY x8), differing only in the His54 / His64 ring
# protonation, plus one acetate per HIP so every cell stays neutral:
#   EE  HIE54 HIE64          rebuilt control (= build_2s2r_amide's protonation, built again for build-to-build noise)
#   ED  HIE54 HID64          His64 tautomer alone
#   DD  HID54 HID64
#   PD  HIP54 HID64          +4 per cell, 4 more ACY
#   MD  His54 HIP in copies 1 and 3, HID in copies 2 and 4; HID64   +2 per cell, 2 more ACY
# Basis: H omit maps from the phenix bxfield gb3 model (other session; ring N-H never modelled, so the difference map
# omits them): His64 HD1 0.54 / HE2 0.15, His54 HD1 0.76 / HE2 0.57, in units of a backbone amide H (+-0.26 null sd).
# A chain-final HID becomes CHID in tleap (leaprc addPdbResMap), so CHID is redefined as the amide exactly like CHIE.
# usage: setup_his_builds.sh        (then launch_his_tests.sh)
set -e
H=/scratch/jamesh/projects/amber/1aho_hurry_test; O=$H/build_2s2r_amide
P=/scratch/jamesh/programs/amber22/miniconda/bin/python3
declare -A H54=( [EE]="HIE HIE HIE HIE" [ED]="HIE HIE HIE HIE" [DD]="HID HID HID HID" [PD]="HIP HIP HIP HIP" [MD]="HIP HID HIP HID" )
declare -A H64=( [EE]="-" [ED]="HID" [DD]="HID" [PD]="HID" [MD]="HID" )
for v in EE ED DD PD MD; do
  d=$H/build_2s2r_amide_h$v
  [ -e $d/Prod.rst7 ] && { echo "$d already built"; continue; }
  mkdir -p $d
  for f in refined.pdb amberme.pdb tleap_stub.in netfrc0_stub.in omegafix_stub.in NH4.cif NH4.frcmod NH4.mol2 NH4.pdb \
           anion.pdb cation.pdb HIS.pdb HIS_settings_asu.txt restraint_suffix.in Invert.in refme.mtz reference.mtz \
           ACY.cif ACY.frcmod ACY.mol2 ACY.pdb xtal_properties.sourceme initial_restraints.pdb c2r_restraints.pdb \
           all_possible_refpoints.pdb salt02_amberme.pdb amide.done; do
    [ -e $O/$f ] && cp -p $O/$f $d/
  done
  ln -sf ../compute_settings.sourceme $d/compute_settings.sourceme
  # protonation: His54 per copy, His64 per copy ("-" keeps the amide build's line: HIS -> CHIS -> the CHIE amide)
  if [ "$v" = EE ]; then
    cp -p $O/protonation.txt $d/
  else
    c=0; : > $d/protonation.txt
    for t in ${H54[$v]}; do
      echo "$t A$((54 + 64*c))" >> $d/protonation.txt
      echo "${H64[$v]} A$((64 + 64*c))" >> $d/protonation.txt
      c=$((c+1))
    done
  fi
  # CHID amide (same atoms, charges and template positions as the CHIE amide; HID and HIE share C/O coordinates)
  if ! grep -q 'CHID = copy HID' $d/tleap_stub.in; then
    awk '{print} /^CHIS = CHIE/{
      print "# CHID: the same His64 amide on the HID tautomer (ND1-H), for protonation.txt \"HID A64\" (tleap maps a"
      print "# chain-final HID to CHID).  Internal ff19SB HID + the NHE cap atoms and charges; total charge 0."
      print "CHID = copy HID"
      print "set CHID.1.O name \"OXT\""
      print "hidN = createAtom O N -0.4630"
      print "set hidN element N"
      print "set hidN position { 6.2522 3.7981 -0.0000 }"
      print "add CHID.1 hidN"
      print "bond CHID.1.C hidN"
      print "hidH1 = createAtom HN1 H 0.2315"
      print "set hidH1 element H"
      print "set hidH1 position { 7.2583 3.7092 -0.0000 }"
      print "add CHID.1 hidH1"
      print "bond hidN hidH1"
      print "hidH2 = createAtom HN2 H 0.2315"
      print "set hidH2 element H"
      print "set hidH2 position { 5.8262 4.7138 -0.0000 }"
      print "add CHID.1 hidH2"
      print "bond hidN hidH2"
      print "set CHID tail null"
      print "set CHID.1 connect1 null"
      print "set CHID.1 name \"HID\""}' $d/tleap_stub.in > $d/stub.new
    mv $d/stub.new $d/tleap_stub.in
  fi
  # one acetate per HIP: a copy of ACY 1 with its C on the placed (B 999) water farthest from protein and other
  # salt (minimum image), waters within 2.6 A of its heavy atoms dropped
  nhip=$(grep -c '^HIP' $d/protonation.txt || true)
  if [ "$nhip" -gt 0 ] && [ ! -e $d/hip_acy.done ]; then
    $P - $d/amberme.pdb $nhip << 'PYEOF'
import sys, numpy as np
f, nadd = sys.argv[1], int(sys.argv[2])
L = open(f).read().split('\n')
cell = np.array([float(x) for x in [l for l in L if l.startswith('CRYST1')][0][6:33].split()])
isat = lambda l: l.startswith(('ATOM', 'HETATM'))
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
atoms = [l for l in L if isat(l)]
prot = np.array([xyz(l) for l in atoms if l[17:20] not in ('HOH', 'WAT', 'ACY', 'NH4') and l[76:78].strip() != 'H'])
salt = [xyz(l) for l in atoms if l[17:20] in ('ACY', 'NH4') and l[76:78].strip() != 'H']
def mind(p, q):
    d = q - p; d -= cell * np.round(d / cell)
    return np.sqrt((d ** 2).sum(1)).min()
cand = [l for l in atoms if l[17:20] in ('HOH', 'WAT') and float(l[60:66]) > 900]
acy1 = [l for l in atoms if l[17:20] == 'ACY' and int(l[22:26]) == 1]
c1 = xyz([l for l in acy1 if l[12:16] == ' C  '][0])
nacy = max(int(l[22:26]) for l in atoms if l[17:20] == 'ACY')
new, picked = [], []
for k in range(nadd):
    best = max(cand, key=lambda l: min(mind(xyz(l), prot), mind(xyz(l), np.array(salt))))
    p = xyz(best); picked.append(p)
    nacy += 1
    for l in acy1:
        q = xyz(l) - c1 + p
        new.append(l[:22] + '%4d' % nacy + l[26:30] + '%8.3f%8.3f%8.3f' % tuple(q) + l[54:])
        if l[76:78].strip() != 'H': salt.append(q)
    print('ACY %d at %s: %.1f A from protein' % (nacy, np.round(p, 2), mind(p, prot)))
    cand.remove(best)
heavy = np.array([xyz(s) for s in new if s[76:78].strip() != 'H'])
lastacy = [l for l in L if isat(l) and l[17:20] == 'ACY'][-1]
out, dropped = [], 0
for l in L:
    if isat(l) and l[17:20] in ('HOH', 'WAT') and mind(xyz(l), heavy) < 2.6:
        dropped += 1; continue
    out.append(l)
    if l is lastacy:
        out += new
print('added %d ACY (%d atoms), dropped %d water atoms; min ACY heavy atom - protein %.2f A' %
      (nadd, len(new), dropped, min(mind(h, prot) for h in heavy)))
open(f, 'w').write('\n'.join(out))
PYEOF
    touch $d/hip_acy.done
  fi
  echo "h$v: $(tr '\n' ' ' < $d/protonation.txt)| ACY $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="ACY" && substr($0,13,4)==" C  "' $d/amberme.pdb | wc -l)  HOH $(grep -c -E ' (HOH|WAT) ' $d/amberme.pdb)  CHID $(grep -c 'CHID = copy' $d/tleap_stub.in)"
done
