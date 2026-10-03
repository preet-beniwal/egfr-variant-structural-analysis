#!/usr/bin/env python3
"""
Phase 3.5 - Merge structural features into the master variant table.

Combines:
    data/variants/luad_egfr_variants_complete.tsv   (conservation + ddG)
    data/variants/features_structural.tsv           (SASA, distances)

Output:
    data/variants/luad_egfr_variants_master.tsv     (final Phase 3 product)
"""

import pandas as pd

COMPLETE = "data/variants/luad_egfr_variants_complete.tsv"
FEATURES = "data/variants/features_structural.tsv"
OUT      = "data/variants/luad_egfr_variants_master.tsv"


def main():
    print("[Phase 3.5] Merging structural features into master table\n")

    variants = pd.read_csv(COMPLETE, sep="\t")
    features = pd.read_csv(FEATURES, sep="\t")
    print(f"  variants: {len(variants)} rows, {len(variants.columns)} cols")
    print(f"  features: {len(features)} rows, {len(features.columns)} cols")

    master = variants.merge(
        features[["residue", "resname", "domain", "sasa",
                  "dist_to_ATP_site", "dist_to_allosteric_pocket"]],
        left_on="anchor_residue", right_on="residue", how="left",
        suffixes=("", "_feat")
    )

    # Drop the duplicate residue column introduced by the merge
    if "residue_feat" in master.columns:
        master = master.drop(columns=["residue_feat"])

    master.to_csv(OUT, sep="\t", index=False)
    print(f"\n  wrote {OUT}")
    print(f"  master shape: {master.shape}")
    print(f"  rows with ΔΔG: {master['ddg_kcal_per_mol'].notna().sum()}")
    print(f"  rows with SASA: {master['sasa'].notna().sum()}")

    # quick summary of what we have per row
    print("\n  complete-table coverage:")
    for col in ["conservation_score", "ddg_kcal_per_mol", "sasa",
                "dist_to_ATP_site"]:
        if col in master.columns:
            n = master[col].notna().sum()
            print(f"    {col:<25} {n}/{len(master)} rows populated")

    print("\n  classic variants with full data:")
    klassik = ["L858R", "T790M", "G719S", "G719C", "L861Q", "S768I"]
    for m in klassik:
        sub = master[master["protein_change"] == m]
        if len(sub):
            v = sub.iloc[0]
            d = v.get("ddg_kcal_per_mol", "NA")
            s = v.get("conservation_score", "NA")
            d_atp = v.get("dist_to_ATP_site", "NA")
            print(f"    {m:<10} cons={s}  ddG={d}  dist_ATP={d_atp}")


if __name__ == "__main__":
    main()
