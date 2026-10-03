#!/usr/bin/env python3
"""
Phase 3.4 - Parse FoldX BuildModel output into a ΔΔG table and merge with
the conservation table from Phase 2.

The Dif file format in FoldX 5.1:
    EGFR_clean_Repair_N_M.pdb  <total_energy>  <term1>  <term2>  ...

    - N = 1-based index into individual_list.txt
    - M = run number (0, 1, 2) for the three replicate runs
    - total_energy = ΔΔG in kcal/mol (mutant - wildtype)

The individual_list.txt uses the format <wt><chain><resnum><mut>; e.g. LA858R.
We convert back to the project's protein-change format (L858R) and merge
with the conservation table.

Output:
    data/variants/luad_egfr_variants_complete.tsv
"""

import os
import re
import pandas as pd

FOLDX_DIR  = "data/foldx"
DDG_FILE   = f"{FOLDX_DIR}/ddg/Dif_EGFR_clean_Repair.fxout"
MUT_FILE   = f"{FOLDX_DIR}/individual_list.txt"
CONS_FILE  = "data/variants/luad_egfr_variants_with_conservation.tsv"
OUT_FILE   = "data/variants/luad_egfr_variants_complete.tsv"


def load_mutation_list(path):
    """Return the ordered list of FoldX mutation strings (e.g. LA858R)."""
    mutations = []
    with open(path) as f:
        for line in f:
            m = line.strip().rstrip(";").strip()
            if m:
                mutations.append(m)
    return mutations


def parse_dif(path, mutations):
    """Parse the Dif file, averaging ΔΔG across the 3 runs per mutation."""
    with open(path) as f:
        lines = f.readlines()

    # join wrapped rows: a new row starts when a line begins with the PDB prefix
    rows = []
    current = ""
    for line in lines:
        if re.match(r"^EGFR_clean_Repair_\d+_\d+\.pdb", line):
            if current:
                rows.append(current)
            current = line.strip()
        elif current:
            current += " " + line.strip()
    if current:
        rows.append(current)

    # aggregate per mutation index
    ddg = {}
    for row in rows:
        parts = row.split()
        if not parts or not parts[0].startswith("EGFR_clean_Repair_"):
            continue
        m = re.search(r"_(\d+)_(\d+)\.pdb", parts[0])
        if not m:
            continue
        idx = int(m.group(1))  # 1-based
        # column 2 is the total ΔΔG
        try:
            energy = float(parts[1])
        except (IndexError, ValueError):
            continue
        ddg.setdefault(idx, []).append(energy)

    records = []
    for idx, values in ddg.items():
        if not (1 <= idx <= len(mutations)):
            continue
        foldx_mut = mutations[idx - 1]
        # parse foldx_mut: <wt><chain><resnum><mut>
        m = re.match(r"^([A-Z])([A-Z])(\d+)([A-Z])$", foldx_mut)
        if not m:
            continue
        wt, chain, resnum, mut = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        records.append({
            "protein_change": f"{wt}{resnum}{mut}",
            "foldx_mutation": foldx_mut,
            "residue": resnum,
            "wt_aa": wt,
            "mut_aa": mut,
            "ddg_kcal_per_mol": round(sum(values) / len(values), 3),
            "ddg_std": round(pd.Series(values).std(), 3) if len(values) > 1 else 0.0,
            "n_runs": len(values),
        })
    return pd.DataFrame(records)


def main():
    print("[Phase 3.4] Parsing FoldX ΔΔG results\n")

    mutations = load_mutation_list(MUT_FILE)
    print(f"  mutation list: {len(mutations)} entries")

    ddg_df = parse_dif(DDG_FILE, mutations)
    print(f"  parsed ΔΔG for {len(ddg_df)} mutations\n")

    print("  ΔΔG distribution:")
    print(f"    mean:          {ddg_df['ddg_kcal_per_mol'].mean():.3f}")
    print(f"    median:        {ddg_df['ddg_kcal_per_mol'].median():.3f}")
    print(f"    std:           {ddg_df['ddg_kcal_per_mol'].std():.3f}")
    print(f"    min:           {ddg_df['ddg_kcal_per_mol'].min():.3f}")
    print(f"    max:           {ddg_df['ddg_kcal_per_mol'].max():.3f}")
    print(f"    destabilising (ΔΔG > 1): {int((ddg_df['ddg_kcal_per_mol'] > 1).sum())}")
    print(f"    stabilising  (ΔΔG < -1): {int((ddg_df['ddg_kcal_per_mol'] < -1).sum())}")

    print("\n  classic EGFR variants:")
    for m in ["L858R", "T790M", "G719S", "G719C", "L861Q", "S768I", "E746_A750del"]:
        sub = ddg_df[ddg_df["protein_change"] == m]
        if len(sub):
            v = sub.iloc[0]
            print(f"    {m:<15} residue {v['residue']:<5}  "
                  f"ΔΔG = {v['ddg_kcal_per_mol']:>7.3f}  "
                  f"(n={v['n_runs']})")

    # ---- merge with conservation ----
    cons = pd.read_csv(CONS_FILE, sep="\t")
    # Note: cons has one row per (protein_change, cohort). Merge on protein_change
    # but keep cohort rows duplicated — we want the conservation per variant,
    # not per patient.

    merged = cons.merge(
        ddg_df[["protein_change", "ddg_kcal_per_mol", "ddg_std", "n_runs"]],
        on="protein_change", how="left"
    )
    merged.to_csv(OUT_FILE, sep="\t", index=False)

    print(f"\n  wrote {OUT_FILE}")
    print(f"  total rows: {len(merged)}")
    print(f"  rows with ΔΔG: {merged['ddg_kcal_per_mol'].notna().sum()}")


if __name__ == "__main__":
    main()
