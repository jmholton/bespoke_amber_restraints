#!/bin/bash
# Same-salt His controls and the NH4+-on-His64-ND1 test at 0.68 M (2026-11-27; user: "maybe try controls with
# same-salt and the uniform charge field hack? yes, also try those next tests").
#
# Same salt (_u): every arm gets DD's amberme.pdb unchanged (11 NH4 / 19 ACY, neutral for HID/HIE).  Extra His54
# charge is NOT countered with acetates; the cell is left charged and PME's uniform neutralizing background absorbs it
# (the "uniform charge field hack").  So arms differ only in His protons; the 0.68 M runs PD/PE/MD differed also in
# 2-4 extra acetates.
#   QD_u  HIP54 in copy 1 only, HID54 in 2-4 (25 % His+), HID64    net +1
#   MD_u  HIP54 copies 1+3, HID54 2+4, HID64                        net +2
#   PD_u  HIP54, HID64                                              net +4
#   PE_u  HIP54, HIE64                                              net +4
# NH4+ test (_n): crystal water HOH51 (3.04 A from His64 ND1, 2.79 A from Tyr42 O; waters 51/172/293/414 in the four
# copies) replaced by NH4+, moving the 4 NH4+ farthest from protein there (salt composition and charge unchanged).
# NH4+ donating to ND1 needs ND1 unprotonated (HIE64): DE_n is the test, DD_n (HID64, ND1-H) the control.
#   DE_n  HID54, HIE64, NH4+ at HOH51                               net 0
#   DD_n  HID54, HID64, NH4+ at HOH51                               net 0
# usage: setup_s68_samesalt.sh [arms]   (then launch_s68_controls.sh QD_u:build_s68_QD_u ...)
set -e
H=/scratch/jamesh/projects/amber/1aho_hurry_test
P=/scratch/jamesh/programs/amber22/miniconda/bin/python3
declare -A SRC=( [QD_u]=build_s68_DD [MD_u]=build_s68_DD [PD_u]=build_s68_DD [PE_u]=build_s68_DD
                 [DE_n]=build_s68_DE [DD_n]=build_s68_DD )
declare -A P54=( [QD_u]="HIP HID HID HID" [MD_u]="HIP HID HIP HID" [PD_u]="HIP HIP HIP HIP" [PE_u]="HIP HIP HIP HIP"
                 [DE_n]="HID HID HID HID" [DD_n]="HID HID HID HID" )
declare -A P64=( [QD_u]=HID [MD_u]=HID [PD_u]=HID [PE_u]=HIE [DE_n]=HIE [DD_n]=HID )
for v in ${*:-QD_u MD_u PD_u PE_u DE_n DD_n}; do
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
  h=${P64[$v]}; [ $h = HIE ] && h=""      # HIE64 blank: tleap default -> the stub's CHIE amide (as EE, DE)
  printf "%s A54\n%s A64\n%s A118\n%s A128\n%s A182\n%s A192\n%s A246\n%s A256\n" $1 "$h" $2 "$h" $3 "$h" $4 "$h" \
    | sed 's/^ *A/ A/' > $d/protonation.txt
  if [[ $v == *_n ]]; then
    $P - $d/amberme.pdb << 'PYEOF'
import sys, numpy as np
f = sys.argv[1]
L = open(f).read().split('\n')
cell = np.array([float(x) for x in [l for l in L if l.startswith('CRYST1')][0][6:33].split()])
isat = lambda l: l.startswith(('ATOM', 'HETATM'))
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
def mind(p, q):
    dd = q - p; dd -= cell * np.round(dd / cell); return np.sqrt((dd ** 2).sum(1)).min()
atoms = [l for l in L if isat(l)]
prot = np.array([xyz(l) for l in atoms if l[17:20] not in ('HOH', 'WAT', 'ACY', 'NH4') and l[76:78].strip() != 'H'])
sites = {}
for k, r in enumerate((51, 172, 293, 414)):
    w = [l for l in atoms if l[17:20] == 'HOH' and int(l[22:26]) == r and l[12:16] == ' O  ']
    assert len(w) == 1, r
    nd1 = xyz([l for l in atoms if int(l[22:26]) == 64 + 64 * k and l[12:16] == ' ND1'][0])
    o42 = xyz([l for l in atoms if int(l[22:26]) == 42 + 64 * k and l[12:16] == ' O  '][0])
    x = xyz(w[0]); dn, do = mind(x, nd1[None]), mind(x, o42[None])
    assert abs(dn - 3.04) < 0.1 and abs(do - 2.79) < 0.1, (r, dn, do)
    sites[r] = (x, w[0][60:66])
nh4 = sorted({int(l[22:26]) for l in atoms if l[17:20] == 'NH4'},
             key=lambda r: -mind(xyz([l for l in atoms if l[17:20] == 'NH4' and int(l[22:26]) == r and l[12:16] == ' N  '][0]), prot))
move = dict(zip(nh4[:4], sites))
out = []
for l in L:
    if isat(l) and l[17:20] == 'HOH' and int(l[22:26]) in sites:
        continue
    if isat(l) and l[17:20] == 'NH4' and int(l[22:26]) in move:
        r = int(l[22:26]); n = xyz([m for m in atoms if m[17:20] == 'NH4' and int(m[22:26]) == r and m[12:16] == ' N  '][0])
        x, b = sites[move[r]]
        l = l[:30] + '%8.3f%8.3f%8.3f' % tuple(xyz(l) - n + x) + '  1.00' + b + l[66:]
    out.append(l)
open(f, 'w').write('\n'.join(out))
for r, w in move.items():
    print('  NH4 %d -> HOH %d site (B %s)' % (r, w, sites[w][1].strip()))
PYEOF
  fi
  echo "s68_$v: $(tr '\n' ' ' < $d/protonation.txt)| NH4 $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="NH4" && substr($0,13,4)==" N  "' $d/amberme.pdb | wc -l) ACY $(awk '/^ATOM|^HETAT/ && substr($0,18,3)=="ACY" && substr($0,13,4)==" C  "' $d/amberme.pdb | wc -l) HOH $(grep -c -E ' (HOH|WAT) ' $d/amberme.pdb)"
done
