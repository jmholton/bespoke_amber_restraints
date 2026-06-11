#! /bin/tcsh -f
#
#
#
set pdbfile = "$1"

set prefix = `basename $pdbfile .pdb`

phenix.omegalyze $pdbfile |\
tee omegalyze_${prefix}.log |\
awk -F ":" 'BEGIN{RTD=45/atan2(1,1)}\
   ! /^SUMMARY|^resid/ {om=$3/RTD;\
   n=substr($0,3,4);c=substr($0,2,1);\
   print "OMEGA",c,n,$0}' |\
sort -k2gr |\
awk 'BEGIN{print "refinement {\n  geometry_restraints.edits {"}\
  {c=$2;n=$3;\
  print "    dihedral {";\
  print "      action = change";\
  print "      atom_selection_1 = \"name CA and resseq",n,"and chain",c,"\"";\
  print "      atom_selection_2 = \"name C  and resseq",n,"and chain",c,"\"";\
  print "      atom_selection_3 = \"name N  and resseq",n+1,"and chain",c,"\"";\
  print "      atom_selection_4 = \"name CA and resseq",n+1,"and chain",c,"\"";\
  print "      angle_ideal = 180.00";\
  print "      sigma = 1";\
  print "    }";}\
END{print "  }\n}"}' >! omega_fix.eff

