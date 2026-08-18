# 1aho Starter Kit — Bespoke Amber Restraints Demo System

Avian pancreatic polypeptide / scorpion toxin (PDB 1aho): the reference demo
system for the bespoke_amber_restraints workflow.  Ultra-high resolution (0.965 Å),
64 residues, P212121, single unit cell (super_mult=1,1,1).  All data is public.

The full worked example is also available as a tarball:
  https://bl831.als.lbl.gov/~jamesh/amber/1aho/example_starter4.tgz


## Why 1aho is the demo system

- Small (64 residues) → fast runs, easy to inspect
- Ultra-high resolution (0.965 Å) → dense, unambiguous electron density
- No ligand → simpler tleap setup, no force-field generation needed
- Public PDB deposit → no private data required
- Single unit cell → no supercell expansion step


## Files in this kit

1aho.pdb
    PDB deposit coordinates.  Used as the starting model directly.
    No private model from a collaborator — this IS the starting point.

1aho.mtz
    Reflection data from the PDB deposit, converted to MTZ.
    Column labels are already correct: FP SIGFP FreeR_flag (plus FREE).
    Symlink as refme_small.mtz — no cad renaming step needed.
    Resolution: 0.965 Å, space group P212121.

ligands/
  ACY.{mol2,frcmod,cif,pdb}  — acetate (crystallisation salt anion, 0.68 M)
  NH4.{mol2,frcmod,cif,pdb}  — ammonium (crystallisation salt cation, 0.68 M)
  No protein ligand — set ligands = "" in xtal_properties.sourceme.

HIS_settings.txt
    Protonation states for the 2 HIS residues: HIE A54, HIE A64.
    Both are epsilon-protonated (HIE) by QM vote.

HIS_votes.txt
    Raw QM votes: residues 54 and 64.
    See HIS_protonation_QM_runme.com in the repo for methodology.

xtal_properties.sourceme
    Crystal properties:
      reso=0.965, P212121, cell 45.9×40.7×30.1 Å, nsymops=4
      super_mult=1,1,1 (single unit cell — no expansion)
      modulo=64, no ligands, salt=NH4+ACY at 0.68 M

seq.fasta
    64-residue protein sequence.


## Key differences from 6c2r workflow

| Property | 1aho | 6c2r (AMPPNP) |
|----------|------|---------------|
| Residues | 64 | 272 |
| Resolution | 0.965 Å | 2.4 Å |
| Space group | P212121 | P41212 |
| super_mult | 1,1,1 | 2,2,2 |
| Protein ligand | none | LIG (AMPPNP) |
| Salt | NH4 + ACY | NH4 + SO4 |
| HIS count | 2 | 10 |
| Private model? | No (use PDB deposit) | Yes (Dirk's model) |
| MTZ labels | Already FP/SIGFP/FreeR_flag | Must rename with cad |
| badlinks | Not needed | Required (Lys-LIG bond) |


## Setup (tcsh)

    source /programs/amber22/amber.csh
    set pdir = ~/projects/git/bespoke_amber_restraints
    set path = ( $pdir $path )

    # Link MTZ — labels already correct, no cad renaming needed
    ln -sf 1aho.mtz refme_small.mtz
    # refme.mtz = refme_small.mtz for 1,1,1 supercell (reindex h1,k1,l1)

    source xtal_properties.sourceme
    cp HIS_settings.txt HIS_settings_asu.txt


## MTZ notes

1aho.mtz has columns: H K L FREE FP SIGFP FreeR_flag
The FreeR_flag column (integer 0/1) is what the workflow uses.
FREE (real) can be ignored.

For the 1,1,1 supercell the cad/reindex expansion step still runs but
is a no-op (multiplier 1 in all directions).  The result refme.mtz = refme_small.mtz.
