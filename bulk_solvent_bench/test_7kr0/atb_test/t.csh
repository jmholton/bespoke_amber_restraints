#! /bin/tcsh -f
cd /programs/amber22 && source amber.csh
cd $1
AddToBox -c cell.pdb -a w.pdb -na 900 -P 0 -RW 2.6 -RP 3 -G 0.2 -o out.pdb >& atb.log
tail -3 atb.log
grep -c WAT out.pdb
head -2 out.pdb
echo help | parmed |& grep -i box
