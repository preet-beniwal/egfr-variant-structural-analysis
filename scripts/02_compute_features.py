#!/usr/bin/env python3
"""
Phase 3.2 - Compute structural features for each variant residue.

For every unique residue in the Phase 2 variant table, compute:
    - domain assignment (based on EGFR canonical residue ranges)
    - SASA (solvent accessible surface area) on the AlphaFold model
    - distance from residue to the ATP site
    - distance from residue to the allosteric pocket from Project 1

The ATP site and allosteric pocket references are taken from published
EGFR structures and from the EAI045 co-crystal in Project 1.

Output:
    data/variants/features_structural.tsv
"""

import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

ALPHAFOLD = "data/structures/alphafold_EGFR_P00533.pdb"
VARIANTS  = "data/variants/luad_egfr_variants_with_conservation.tsv"
OUT       = "data/variants/features_structural.tsv"

# EGFR canonical domain boundaries (UniProt P00533 numbering)
DOMAINS = [
    ("signal_peptide",   1,   24),
    ("extracellular",   25,  645),
    ("transmembrane",  646,  668),
    ("juxtamembrane",  669,  706),
    ("kinase",         707,  982),
    ("c_terminal",     983, 1210),
]

# ATP site reference residues (EGFR kinase domain, canonical numbering)
ATP_SITE = [745, 790, 793, 797, 855, 856]

# Allosteric pocket residues (from Project 1's PLIF analysis)
ALLOSTERIC = [745, 856, 858]


def assign_domain(resnum):
    for name, start, end in DOMAINS:
        if start <= resnum <= end:
            return name
    return "unknown"


def get_atoms(structure, resnum, atom_name=None):
    """Return all atoms of the residue across all chains."""
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                if res.id[1] == resnum:
                    if atom_name is None:
                        atoms.extend(res.get_atoms())
                    else:
                        if atom_name in res:
                            atoms.append(res[atom_name])
        break
    return atoms


def min_distance(atoms_a, atoms_b):
    """Minimum distance between any atom in a and any atom in b."""
    if not atoms_a or not atoms_b:
        return float("nan")
    coords_a = np.array([a.coord for a in atoms_a])
    coords_b = np.array([b.coord for b in atoms_b])
    d = np.linalg.norm(coords_a[:, None, :] - coords_b[None, :, :], axis=2)
    return float(d.min())


def main():
    print("[Phase 3.2] Computing structural features\n")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("egfr", ALPHAFOLD)
    print(f"  loaded {ALPHAFOLD}")

    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    print("  computed SASA (Shrake-Rupley)")

    variants = pd.read_csv(VARIANTS, sep="\t")
    unique_residues = sorted(set(
        int(v["anchor_residue"])
        for _, v in variants.iterrows()
        if pd.notna(v["anchor_residue"])
    ))
    print(f"  {len(unique_residues)} unique variant residues\n")

    atp_atoms = []
    for r in ATP_SITE:
        atp_atoms.extend(get_atoms(structure, r))

    allo_atoms = []
    for r in ALLOSTERIC:
        allo_atoms.extend(get_atoms(structure, r))

    print(f"  ATP site reference atoms:       {len(atp_atoms)}")
    print(f"  Allosteric pocket atoms:        {len(allo_atoms)}\n")

    rows = []
    for resnum in unique_residues:
        atoms = get_atoms(structure, resnum)
        if not atoms:
            continue

        resname = None
        sasa = 0.0
        for model in structure:
            for chain in model:
                for res in chain:
                    if res.id[1] == resnum and res.id[0] == " ":
                        resname = res.resname
                        sasa = res.sasa
                        break
            break

        rows.append({
            "residue": resnum,
            "resname": resname,
            "domain": assign_domain(resnum),
            "sasa": round(sasa, 2),
            "dist_to_ATP_site": round(min_distance(atoms, atp_atoms), 2),
            "dist_to_allosteric_pocket": round(min_distance(atoms, allo_atoms), 2),
        })

    features = pd.DataFrame(rows)
    features.to_csv(OUT, sep="\t", index=False)

    print(f"  wrote {OUT}")
    print(f"  rows: {len(features)}\n")

    print("  domain distribution of variant residues:")
    for d, n in features["domain"].value_counts().items():
        print(f"    {d:<18} {n}")

    print("\n  sample features (first 8 residues):")
    print(features.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
