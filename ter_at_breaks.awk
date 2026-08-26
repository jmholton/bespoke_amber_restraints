#! /bin/awk -f
# Insert a TER wherever the protein backbone is actually broken: a peptide
# C(i)-N(i+1) distance beyond a cutoff (default 2.5 A), or a chain-id change with
# no TER already present.
#
# This is a safety net for models whose chain terminators would otherwise be
# found only from OXT atoms (as convert_pdb output=amber does).  A raw PDB deposit
# run without buildout (hurry mode), or supercell copies that got renumbered into
# one continuous chain, carry no OXT at the break - so tleap sees two
# consecutively-numbered residues in one chain and builds a spurious peptide bond
# straight across the gap (a 20-130 A "bond"), which then blows up the geometry.
# Passing the model through here first turns every real break into a genuine chain
# termination, so tleap caps the fragments instead of bonding them.
#
# Usage:  ter_at_breaks.awk [-v cutoff=2.5] [-v debug=1] model.pdb
BEGIN{
  if(cutoff=="") cutoff = 2.5
  # standard amino acids + common amber / protonation-state variants
  split("ALA ARG ASN ASP CYS CYX CYM GLN GLU GLY HIS HID HIE HIP ILE LEU " \
        "LYS LYN MET PHE PRO SER THR TRP TYR VAL ASH GLH HYP MSE SEC PYL", a, " ")
  for(i in a) isres[a[i]] = 1
}
/^TER/ { print; haveC = 0; prevkey = ""; next }
$0 !~ /^ATOM|^HETAT/ { print; next }
{
  atom  = substr($0,13,4); gsub(/ /,"",atom)
  resn  = substr($0,18,3); gsub(/ /,"",resn)
  chain = substr($0,22,1)
  resi  = substr($0,23,5)                 # resnum + insertion code
  x = substr($0,31,8)+0; y = substr($0,39,8)+0; z = substr($0,47,8)+0
  key = chain "|" resi
  aa  = isres[resn]
}
# arriving at a new amino-acid residue (its backbone N): is the peptide bond from
# the previous residue's C broken?  If so, terminate the chain before this residue.
aa && atom == "N" && key != prevkey && haveC {
  dd = sqrt((cx-x)^2 + (cy-y)^2 + (cz-z)^2)
  if( chain != prevchain || dd > cutoff ) {
    if(debug) printf("REMARK ter_at_breaks: TER before %s %s (C-N = %.1f A)\n", resn, key, dd)
    print "TER"
    nter++
  }
  haveC = 0
}
{ print }
aa && atom == "C" { cx = x; cy = y; cz = z; haveC = 1 }
aa && atom == "N" { prevkey = key; prevchain = chain }
END{ if(debug) printf("REMARK ter_at_breaks: inserted %d TER(s), cutoff %.1f A\n", nter+0, cutoff) > "/dev/stderr" }
