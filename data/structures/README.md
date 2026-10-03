# Structures

Six PDB files in two categories.

## AlphaFold DB (full-length, predicted)

| File | Protein | UniProt | Residues | Source |
|------|---------|---------|----------|--------|
| alphafold_EGFR_P00533.pdb | Human EGFR | P00533 | 1,210 | AlphaFold DB v6 |
| alphafold_HER2_P04626.pdb | Human HER2/ERBB2 | P04626 | 1,255 | AlphaFold DB v6 |

Full-length models including signal peptide and all domains.

## Experimental PDB (kinase domain)

| File | PDB | Contents | Resolution |
|------|-----|----------|------------|
| 1M17.pdb | 1M17 | EGFR kinase + erlotinib | 2.6 A |
| 3POZ.pdb | 3POZ | EGFR T790M kinase | 2.9 A |
| 5D41.pdb | 5D41 | EGFR T790M/V948R + EAI045 | 2.8 A |
| 3PP0.pdb | 3PP0 | HER2 kinase | 2.5 A |

Kinase-domain only, with co-crystallized ligands and waters where present.
Used for high-resolution context of clinical variants in the kinase domain.

## Coordinate systems

Two numbering conventions are used in this project:

1. UniProt P00533 canonical numbering -- used by AlphaFold DB and by
   the variant annotations from Project 3 (VEP protein change strings).
2. PDB structure numbering -- used by experimental crystal structures.
   These often start at a different residue number (e.g., 1M17 starts
   at residue 672 in EGFR canonical numbering).

Mapping between them is done in Phase 2.
