#! /bin/tcsh -f
#
#
#
set pdbid = 6c2r
set sequence = ""

getcif.com $pdbid

set pdbfile = ${pdbid}.pdb
set firstresnum = `sequence.awk $pdbfile | awk '/^seqres start/{print $NF}'`
if( "$firstresnum" == "" ) set firstresnum = 1

if(! -e "$sequence") then
  # get full sequence from PDB
  grep SEQRES $pdbfile |\
  sequence.awk |\
  awk 'NR==1{$0="> "$0} {print} NF==0{exit}' |\
  tee seq.fasta
  set sequence = seq.fasta
  awk '/SEQRES/{print substr($0,18)}' $pdbfile |\
  awk '{for(i=1;i<=NF;++i)print $i}' >! tlc.txt
endif
set sequence = seq.fasta

sequence.awk -v tlc=1 $sequence | awk '/^[A-V][A-Y][A-Y]$/' >! tlc.txt
awk '{print "BUILD",$1,-60,-40}' tlc.txt | build_n2c.awk >! backbone.pdb
awk '{print "BUILD",$1, ++n,"? ? ? ? ?"}' tlc.txt |\
cat - backbone.pdb |\
build_side.awk >! side.pdb
echo "renumber ${firstresnum}\nchain A" | pdbset xyzin side.pdb xyzout helix.pdb

egrep "^ATOM|^HETAT" helix.pdb |\
awk '{c=substr($0,22,1)} c==" "{c="_"}\
   {print substr($0,12,5),substr($0,18,3),c,substr($0,23,6)}' |\
awk '{while(gsub("  "," "));gsub(" $","");print $0,"EXPECT"}' |\
sort -u |\
sort -k3,3 -k4g >! expected_atoms.txt



# check if anything is missing
awk '/^ATOM|^HETAT/{print substr($0,1,16),substr($0,18)}' $pdbfile |\
awk '{id=substr($0,12,15)} ! seen[id]{print;++seen[id]}' |\
cat >! noalt.pdb

egrep "^ATOM|^HETAT" noalt.pdb |\
awk '{c=substr($0,22,1)} c==" "{c="_"}\
   {print substr($0,12,5),substr($0,18,3),c,substr($0,23,6)}' |\
awk '{while(gsub("  "," "));gsub(" $","");print $0,"EXIST"}' |\
sort -u |\
sort -k3,3 -k4g >! existing_atoms.txt

cat existing_atoms.txt expected_atoms.txt |\
awk '{id=$1" "$2" "$3" "$4}\
    $NF=="EXIST"{++exist[id]}\
    $NF=="EXPECT" && ! exist[id]{print id,"MISSING"}' |\
tee missing_atoms.txt |\
awk '{print $2,$3,$4}' |\
sort -u | sort -k1.6g |\
tee missing_residues.txt
set nmissres = `cat missing_residues.txt | wc -l`

rm -f add2Nterm.txt add2Cterm.txt rebuildme.txt
touch add2Nterm.txt add2Cterm.txt rebuildme.txt
foreach line ( `seq $nmissres -1 1` )
   set info = `tail -n +$line missing_residues.txt | head -n 1`
   set typ = $info[1]
   set chain = $info[2]
   set resnum = $info[3]

   echo $info |\
   cat - existing_atoms.txt |\
   awk 'NR==1{\
      mres="CA "$2" " $3;\
      pres="CA "$2" " ($3+1);\
      nres="CA "$2" " ($3-1);\
      next}\
      {id=$1" "$3" "$4}\
      id==pres{++gotpres}\
      id==nres{++gotnres}\
    END{print gotnres+0,gotpres+0}' >! gotres.txt
   set gotres = `cat gotres.txt`

   if( "$gotres" == "1 1" ) then
     echo "$info" | tee -a rebuildme.txt
   endif
   if( "$gotres" == "0 1" ) then
     echo "$info" | tee -a add2Nterm.txt
   endif
   if( "$gotres" == "1 0" ) then
     echo "$info" | tee -a add2Cterm.txt
   endif
   echo "CA $info" >> existing_atoms.txt
end

cat << EOF >! tack.py
imol_coords = handle_read_draw_molecule("phenix0_001.pdb")
imol_map = make_and_draw_map("phenix0_001.mtz","2FOFCWT","PH2FOFCWT","",0,0)
EOF

foreach line ( `awk '{print NR}' add2Nterm.txt` )
   set info = `awk -v i=$line 'NR==i' add2Nterm.txt`
   set TYP = $info[1]
   set chain = $info[2]
   set resnum = $info[3]
   @ nres = ( $resnum + 1 )

   cat << EOF >> tack.py
set_go_to_atom_chain_residue_atom_name("${chain}",${nres}," CA ")
add_terminal_residue(imol_coords, "${chain}",${nres}, "${TYP}", 1)
with_auto_accept([sphere_refine, 3.5])
EOF
end

foreach line ( `awk '{print NR}' add2Cterm.txt` )
   set info = `tac add2Cterm.txt | awk -v i=$line 'NR==i'`
   set TYP = $info[1]
   set chain = $info[2]
   set resnum = $info[3]
   @ pres = ( $resnum - 1 )

   cat << EOF >> tack.py
set_go_to_atom_chain_residue_atom_name("${chain}",${pres}," CA ")
add_terminal_residue(imol_coords, "${chain}",${pres}, "", "${TYP}", 0)
with_auto_accept([sphere_refine, 3.5])
EOF
end

foreach line ( `awk '{print NR}' rebuildme.txt` )
   set info = `cat rebuildme.txt | awk -v i=$line 'NR==i'`
   set TYP = $info[1]
   set chain = $info[2]
   set resnum = $info[3]

   cat << EOF >> tack.py
set_go_to_atom_chain_residue_atom_name("${chain}",$resnum," CA ")
mutate_and_auto_fit(${resnum},"${chain}",imol_coords,imol_map,"${TYP}")
with_auto_accept([sphere_refine, 3.5])
EOF
end

cat << EOF >> tack.py
save_coordinates(0,"coot.pdb")
coot_no_state_real_exit(0)
EOF

coot --no-graphics --script tack.py



