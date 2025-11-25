#! /bin/awk -f
#
#  jiffy script for filtering PDB files by content
#
#
BEGIN{
   # may be: all atoms protein notprotein ligand notligand water notwater notEP zeroH
   if(! only) only = "all"
   if(! skip) skip = "none"

   # for pattern matching
   only = ","only","
   skip = ","skip","

   # allow user-specified patterns
   if(ligand) ligands = ligand
   if(salt) salts = salt
   if(! ligands) ligands = ""
   if(! salts) salts = ""

   # default
   if(only !~ /protein,|ligand|salt|water,|EP,|H,/) only = ","only",all,"
   if(only ~ /,all,/) only = ","only",protein,ligand,salt,water,EP,zeroH,"

   if(debug) print "DEBUG only =",only
}

! /^ATOM|HETAT/{
    # scrape out non-atom records if asked
    if( only ~ /,atoms,/  ) next;
    # otherwise, always pass through
    print;
    next;
}

# parse the atom line
{
    if(debug) print "DEBUG" substr($0,6);
    atom=substr($0,12,5);
      atm=atom;gsub(" ","",atm);
    typ = substr($0,18,3);
    Ee = substr($0,77,4);gsub(" ","",Ee);
#    atomEe=substr($0,13,2);gsub(" ","",atomEe);
    protein=water=ligand=salt=EP=0;
}


# classify the residue
typ~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL|MET|MSE/{++protein}
typ~/HID|HIE|HIP|HIS|ILE|LEU|LYS|KCX|PHE|PRO|SER|THR|TRP|TYR/{++protein}
typ~/HOH|WAT/{++water}
ligands!="" && typ~ligands{++ligand}
salts!="" && typ~salts{++salt}
! protein && ! water && ! salt && ligands==""{++ligand}
! protein && ! water && ! ligand && salts==""{++salt}

debug > 5 {print "DEBUG: protein",protein,"ligand",ligand,"salt",salt,"water",water,"prevprotein",prevprotein}

# skip if asked
protein && skip ~ /,protein,/{next}
water && skip ~ /,water,/{next}
ligand && skip ~ /,ligand/{next}
salt && skip ~ /,salt/{next}
Ee=="XP" && skip ~ /,EP,|,H,/{next}
atm=="EPW" && skip ~ /,EP,|,H,/{next}
Ee=="Y" && water && skip ~ /,EP,|,H,/{next}
atm=="Y1" && water && skip ~ /,EP,|,H,/{next}
Ee=="H" && occ==0 && skip ~ /,zeroH,/{next}
Ee=="H" && skip ~ /,H,/{next}


protein   && ( only !~ /,protein,/ ) {next}
water     && ( only !~ /,water,/ ) {next}
ligand    && ( only !~ /,ligand/) {next}
salt      && ( only !~ /,salt/) {next}
! protein && ! water &&  ! ligand && ( only ~ /,ligand/) {next}
protein && only ~ /notprotein|nonprotein/ {next}
ligand && only ~ /notligand|nonligand/ {next}
water && only ~ /notwater|nonwater/ {next}
Ee=="XP" && only ~ /noEP|notEP|nonEP/{next}
Ee=="H" && occ==0 && only ~ /nonzeroH|nozeroH|notzeroH/{next}

debug > 5 {print "DEBUG: not skipped"}


# actual, final printing out
{
   print;
}

