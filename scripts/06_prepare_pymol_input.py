#!/usr/bin/env python3
"""
Phase 5.1 - Prepare a PDB with variant scores written into the B-factor column.

PyMOL can color residues by B-factor, so this lets us map ΔΔG directly onto
the structure. Residues without ΔΔG (deletions, insertions) keep B-factor 0.

Output:
    data/structures/alphafold_EGFR_with_ddg.pdb
    data/structures/alphafold_EGFR_with_conservation.pdb
"""

import pandas as pd

INPUT  = "data/structures/alphafold_EGFR_P00533.pdb"
MASTER = "data/variants/luad_egfr_variants_master.tsv"
OUT_DDG = "data/structures/alphafold_EGFR_with_ddg.pdb"
OUT_CONS = "data/structures/alphafold_EGFR_with_conservation.pdb"


def make_score_map():
    df = pd.read_csv(MASTER, sep="\t")
    df = df.drop_duplicates("protein_change")
    # average across variants at the same residue (rare; usually 1)
    ddg_map = (df.dropna(subset=["ddg_kcal_per_mol"])
                 .groupby("anchor_residue")["ddg_kcal_per_mol"].mean().to_dict())
    cons_map = (df.dropna(subset=["conservation_score"])
                  .groupby("anchor_residue")["conservation_score"].mean().to_dict())
    return ddg_map, cons_map


def rewrite_bfactor(in_path, out_path, residue_to_score, scale=1.0):
    """Rewrite the B-factor column so non-variant residues have B=0 and
    variant residues carry their score (shifted non-negative)."""
    with open(in_path) as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    f_out.write(line)
                    continue
                if resnum in residue_to_score:
                    b = residue_to_score[resnum] * scale
                else:
                    b = 0.0
                b = max(b, 0.0)
                b = min(b, 99.99)
                line = line[:60] + f"{b:6.2f}" + line[66:]
            f_out.write(line)
def main():
    print("[Phase 5.1] Preparing PyMOL-ready structures\n")
    ddg_map, cons_map = make_score_map()
    print(f"  residues with ΔΔG:        {len(ddg_map)}")
    print(f"  residues with conservation: {len(cons_map)}\n")

    # shift ΔΔG by +3 so all values are positive (PyMOL B-factor is non-negative)
    ddg_shifted = {k: v + 3.0 for k, v in ddg_map.items()}
    rewrite_bfactor(INPUT, OUT_DDG, ddg_shifted)
    print(f"  wrote {OUT_DDG}")

    rewrite_bfactor(INPUT, OUT_CONS, cons_map)
    print(f"  wrote {OUT_CONS}")


if __name__ == "__main__":
    main()
