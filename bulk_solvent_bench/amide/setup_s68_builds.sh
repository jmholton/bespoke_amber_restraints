#!/bin/bash
# Charge tests at the crystallisation salt (2026-11-25; user: "the real crystal was formed in salt ... lets look at
# charge effects at close to what we think the actual ionic strength was").  Six builds sharing one 0.68 M ammonium
# acetate placement (salt068_base.csh: 11 NH4 + 15 ACY at cell charge +4); only the counter-acetates follow the charge:
#   cEE  carboxylate His64 (CHIE), HIE54             +4   11 NH4 / 15 ACY
#   EE   amide His64 (HIE), HIE54                    +8   11 / 19
#   ED   amide His64 as HID, HIE54                   +8   11 / 19
#   DD   amide HID64, HID54                          +8   11 / 19
#   MD   amide HID64, HIP54 copies 1+3, HID54 2+4    +10  11 / 21
#   PD   amide HID64, HIP54                          +12  11 / 23
# Each amberme.pdb: the base, Cys12/Cys63 alt-B grafted into copies 3 and 4 (as run_local.csh's build, 2 strained +
# 2 relaxed), then the extra acetates: copies of ACY 1 with the C on the placed (B 999) water farthest from protein and
# other salt (minimum image; greedy, so the arms' extra acetates are nested: EE/ED/DD's 4 are MD's first 4, MD's 6 are
# PD's first 6), waters within 2.6 A of the new heavy atoms dropped.  Protonation and tleap_stub.in (CHIE/CHID amide)
# from the 2026-11-25 no-salt His builds (build_2s2r_amide_h*), the carboxylate's from build_2s2r_cooh.
# usage: setup_s68_builds.sh        (then launch_s68_tests.sh)
set -e
H=/scratch/jamesh/projects/amber/1aho_hurry_test; base=$H/salt068_base
P=/scratch/jamesh/programs/amber22/miniconda/bin/python3
declare -A SRC=( [cEE]=build_2s2r_cooh [EE]=build_2s2r_amide_hEE [ED]=build_2s2r_amide_hED [DD]=build_2s2r_amide_hDD
                 [MD]=build_2s2r_amide_hMD [PD]=build_2s2r_amide_hPD )
declare -A NADD=( [cEE]=0 [EE]=4 [ED]=4 [DD]=4 [MD]=6 [PD]=8 )
for v in cEE EE ED DD MD PD; do
  d=$H/build_s68_$v; s=$H/${SRC[$v]}
  [ -e $d/Prod.rst7 ] && { echo "$d already built"; continue; }
  mkdir -p $d
  for f in refined.pdb protonation.txt tleap_stub.in netfrc0_stub.in omegafix_stub.in NH4.cif NH4.frcmod NH4.mol2 \
           NH4.pdb anion.pdb cation.pdb HIS.pdb HIS_settings_asu.txt restraint_suffix.in Invert.in refme.mtz \
           reference.mtz ACY.cif ACY.frcmod ACY.mol2 ACY.pdb xtal_properties.sourceme; do
    [ -e $s/$f ] && cp -p $s/$f $d/
  done
  ln -sf ../compute_settings.sourceme $d/compute_settings.sourceme
  cp -p $base/amberme.pdb $d/salt068_amberme.pdb
  # graft: alt-B Cys12/Cys63 of copies 3 and 4 (run_local.csh's build, verbatim)
  awk 'FNR==NR{if(/^ATOM|^HETAT/ && substr($0,17,1)=="B" && substr($0,18,3)=="CYS"){r=substr($0,23,4)+0;\
        if(r==140||r==191||r==204||r==255) xyz[r" "substr($0,13,4)]=substr($0,31,24)}; next}\
       /^ATOM|^HETAT/ && substr($0,18,3)=="CYS"{k=(substr($0,23,4)+0)" "substr($0,13,4); if(k in xyz){\
        $0=substr($0,1,30) xyz[k] "  0.45" substr($0,61); ++n}} {print}\
       END{print n" atoms grafted" > "/dev/stderr"}' $H/centroids/centroids_in_density.pdb $d/salt068_amberme.pdb > $d/amberme.pdb
  if [ ${NADD[$v]} -gt 0 ]; then
    $P - $d/amberme.pdb ${NADD[$v]} << 'PYEOF'
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
new = []
for k in range(nadd):
    best = max(cand, key=lambda l: min(mind(xyz(l), prot), mind(xyz(l), np.array(salt))))
    p = xyz(best); nacy += 1
    for l in acy1:
        q = xyz(l) - c1 + p
        new.append(l[:22] + '%4d' % nacy + l[26:30] + '%8.3f%8.3f%8.3f' % tuple(q) + l[54:])
        if l[76:78].strip() != 'H': salt.append(q)
    print('  ACY %d at %s: %.1f A from protein' % (nacy, np.round(p, 2), mind(p, prot)))
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
print('  added %d ACY, dropped %d waters; closest new ACY heavy atom to protein %.2f A, to other salt %.2f A' %
      (nadd, dropped, min(mind(h, prot) for h in heavy),
       min(mind(h, np.array([s for s in salt[:len(salt) - len(heavy)]])) for h in heavy)))
open(f, 'w').write('\n'.join(out))
PYEOF
  fi
  echo "s68_$v: $(tr '\n' ' ' < $d/protonation.txt)| NH4 $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="NH4" && substr($0,13,4)==" N  "' $d/amberme.pdb | wc -l) ACY $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="ACY" && substr($0,13,4)==" C  "' $d/amberme.pdb | wc -l) HOH $(grep -c -E ' (HOH|WAT) ' $d/amberme.pdb)  amide stub: $(grep -c 'createAtom' $d/tleap_stub.in)"
done
