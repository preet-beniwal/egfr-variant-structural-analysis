#!/usr/bin/env python3
"""
Phase 4 - Three-layer correlation analysis.

Correlates:
    1. Evolutionary conservation vs structural destabilisation (per variant)
    2. Mean ΔΔG per EGFR class vs clinical hazard ratio (class level)
    3. Conservation vs clinical HR (class level)

Inputs:
    data/variants/luad_egfr_variants_master.tsv
    ~/EGFR-NSCLC-Clinical-Survival-Analysis/results/cox_stratified_OS.tsv

Outputs:
    figures/phase4_conservation_vs_ddg.png
    figures/phase4_class_ddg_vs_hr.png
    figures/phase4_class_summary.png
    results/phase4_correlations.tsv
    results/phase4_class_table.tsv
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

MASTER = "data/variants/luad_egfr_variants_master.tsv"
COX_TCGA = os.path.expanduser(
    "~/EGFR-NSCLC-Clinical-Survival-Analysis/results/cox_stratified_OS.tsv")
OUT_FIG = "figures"
OUT_RES = "results"
os.makedirs(OUT_FIG, exist_ok=True)
os.makedirs(OUT_RES, exist_ok=True)


def spearman(x, y):
    """Spearman correlation with p-value, handling NaN."""
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 3:
        return float("nan"), float("nan"), 0
    rho, p = stats.spearmanr(x[mask], y[mask])
    return rho, p, mask.sum()


def classify_to_egfr_class(pc):
    """Map a protein change string to an EGFR clinical class."""
    if pd.isna(pc):
        return None
    pc = str(pc)
    if pc == "L858R":
        return "L858R"
    if pc == "T790M":
        return "Rare"
    if pc.startswith("E746") or pc.startswith("L747") or pc.startswith("A750") \
       or pc.startswith("T751") or pc.startswith("E749") or pc.startswith("P753"):
        return "Exon19del"
    if pc in ("G719A", "G719C", "G719D", "G719S", "L861Q", "S768I"):
        return "Rare"
    if "_" in pc and "ins" in pc:
        return "Rare"
    return "Other"


def main():
    print("[Phase 4] Three-layer correlation analysis\n")

    # ---------- load ----------
    master = pd.read_csv(MASTER, sep="\t")
    print(f"  master table: {len(master)} rows")

    subs = master[master["ddg_kcal_per_mol"].notna()].copy()
    subs = subs.drop_duplicates("protein_change")
    print(f"  substitutions with ΔΔG: {len(subs)}\n")

    # ---------- 4.1 conservation vs ΔΔG ----------
    print("  [4.1] Conservation score vs ΔΔG (per variant)")
    rho1, p1, n1 = spearman(
        subs["conservation_score"].values.astype(float),
        subs["ddg_kcal_per_mol"].values
    )
    print(f"        Spearman rho = {rho1:.3f}  p = {p1:.4f}  n = {n1}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(subs["conservation_score"], subs["ddg_kcal_per_mol"],
               s=30, alpha=0.6, color="#1f77b4", edgecolor="black", linewidth=0.5)
    # jitter x slightly since conservation is integer 1-3
    jitter = np.random.RandomState(42).uniform(-0.1, 0.1, len(subs))
    ax.scatter(subs["conservation_score"] + jitter, subs["ddg_kcal_per_mol"],
               s=30, alpha=0.6, color="#1f77b4", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Conservation score (0-3)")
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_title(f"Conservation vs structural destabilisation (n={n1})\n"
                 f"Spearman ρ = {rho1:.3f}, p = {p1:.4f}", fontsize=11)
    ax.axhline(0, color="grey", linestyle="--", alpha=0.5)
    ax.set_xticks([1, 2, 3])
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{OUT_FIG}/phase4_conservation_vs_ddg.png", dpi=300)
    plt.close()
    print(f"        wrote figures/phase4_conservation_vs_ddg.png")

    # ---------- 4.2 class-level ΔΔG vs HR ----------
    print("\n  [4.2] Mean ΔΔG per EGFR class vs clinical HR")
    subs["egfr_class"] = subs["protein_change"].apply(classify_to_egfr_class)
    class_ddg = subs.groupby("egfr_class")["ddg_kcal_per_mol"].agg(
        ["mean", "median", "std", "count"]).reset_index()
    print(f"        class ΔΔG summary:")
    print(class_ddg.to_string(index=False))

    cox = pd.read_csv(COX_TCGA, sep="\t", index_col=0)
    hr_map = {
        "Exon19del": "egfr_Exon19del",
        "L858R":     "egfr_L858R",
        "Other":     "egfr_Other",
        "Rare":      "egfr_Rare",
    }
    class_ddg["clinical_HR"] = class_ddg["egfr_class"].map(
        lambda c: cox.loc[hr_map[c], "HR"] if c in hr_map else None
    )
    class_ddg = class_ddg.dropna(subset=["clinical_HR"])
    print(f"\n        class-level comparison:")
    print(class_ddg.to_string(index=False))

    rho2, p2, n2 = spearman(class_ddg["mean"].values, class_ddg["clinical_HR"].values)
    print(f"\n        Spearman rho = {rho2:.3f}  p = {p2:.4f}  n = {n2}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(class_ddg["mean"], class_ddg["clinical_HR"],
               s=120, color="#d62728", edgecolor="black", linewidth=1.2, zorder=3)
    for _, row in class_ddg.iterrows():
        ax.annotate(row["egfr_class"],
                    (row["mean"], row["clinical_HR"]),
                    xytext=(8, 5), textcoords="offset points", fontsize=10)
    ax.axhline(1.0, color="grey", linestyle="--", alpha=0.5, label="HR = 1 (no effect)")
    ax.set_xlabel("Mean ΔΔG (kcal/mol)")
    ax.set_ylabel("Hazard ratio (overall survival)")
    ax.set_title(f"Structural destabilisation vs clinical HR\n"
                 f"Spearman ρ = {rho2:.3f}, p = {p2:.4f}, n = {n2}", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{OUT_FIG}/phase4_class_ddg_vs_hr.png", dpi=300)
    plt.close()
    print(f"        wrote figures/phase4_class_ddg_vs_hr.png")

    # ---------- 4.3 conservation vs HR per class ----------
    print("\n  [4.3] Conservation score vs clinical HR")
    class_cons = subs.groupby("egfr_class")["conservation_score"].mean().reset_index()
    class_cons = class_cons.merge(
        class_ddg[["egfr_class", "clinical_HR"]], on="egfr_class", how="inner"
    )
    print(class_cons.to_string(index=False))

    # ---------- 4.4 composite summary figure ----------
    print("\n  [4.4] Composite summary")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # panel A: conservation vs ΔΔG
    axes[0].scatter(subs["conservation_score"] + jitter,
                    subs["ddg_kcal_per_mol"],
                    s=30, alpha=0.6, color="#1f77b4",
                    edgecolor="black", linewidth=0.5)
    axes[0].set_xlabel("Conservation score")
    axes[0].set_ylabel("ΔΔG (kcal/mol)")
    axes[0].set_title(f"Conservation vs ΔΔG\nρ = {rho1:.3f}", fontsize=10)
    axes[0].set_xticks([1, 2, 3])
    axes[0].grid(True, alpha=0.3)

    # panel B: mean ΔΔG by class
    order = ["Exon19del", "L858R", "Other", "Rare"]
    order = [o for o in order if o in class_ddg["egfr_class"].values]
    class_ddg_ord = class_ddg.set_index("egfr_class").loc[order].reset_index()
    axes[1].bar(class_ddg_ord["egfr_class"], class_ddg_ord["mean"],
                color="#2ca02c", edgecolor="black", linewidth=1.2)
    axes[1].errorbar(class_ddg_ord["egfr_class"], class_ddg_ord["mean"],
                     yerr=class_ddg_ord["std"], fmt="none",
                     ecolor="black", capsize=5)
    axes[1].set_ylabel("Mean ΔΔG (kcal/mol)")
    axes[1].set_title("Mean ΔΔG per EGFR class", fontsize=10)
    axes[1].axhline(0, color="grey", linestyle="--", alpha=0.5)
    axes[1].tick_params(axis="x", rotation=20)
    axes[1].grid(True, alpha=0.3, axis="y")

    # panel C: HR by class
    axes[2].bar(class_ddg_ord["egfr_class"], class_ddg_ord["clinical_HR"],
                color="#9467bd", edgecolor="black", linewidth=1.2)
    axes[2].axhline(1.0, color="red", linestyle="--", alpha=0.6)
    axes[2].set_ylabel("Hazard ratio (OS)")
    axes[2].set_title("Clinical HR per EGFR class (TCGA)", fontsize=10)
    axes[2].tick_params(axis="x", rotation=20)
    axes[2].grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(f"{OUT_FIG}/phase4_class_summary.png", dpi=300)
    plt.close()
    print(f"        wrote figures/phase4_class_summary.png")

    # ---------- save tables ----------
    corr_table = pd.DataFrame([
        {"analysis": "conservation_vs_ddg (per variant)",
         "n": n1, "spearman_rho": round(rho1, 3), "p_value": round(p1, 4)},
        {"analysis": "class_ddg_vs_HR",
         "n": n2, "spearman_rho": round(rho2, 3), "p_value": round(p2, 4)},
    ])
    corr_table.to_csv(f"{OUT_RES}/phase4_correlations.tsv", sep="\t", index=False)
    class_ddg.to_csv(f"{OUT_RES}/phase4_class_table.tsv", sep="\t", index=False)
    print(f"\n  wrote results/phase4_correlations.tsv")
    print(f"  wrote results/phase4_class_table.tsv")

    print("\n  === Interpretation ===")
    if abs(rho1) > 0.3 and p1 < 0.05:
        print(f"  Conservation and ΔΔG are correlated (ρ={rho1:.2f}, p={p1:.4f}).")
        print("  Evolutionarily constrained positions tend to resist substitution structurally.")
    else:
        print(f"  Conservation and ΔΔG show weak correlation (ρ={rho1:.2f}, p={p1:.4f}).")
        print("  The two layers capture different aspects of variant impact.")

    if abs(rho2) > 0.5 and p2 < 0.2:
        print(f"\n  Class-level ΔΔG and HR align (ρ={rho2:.2f}, n={n2}).")
        print("  Classes with higher mean destabilisation tend to have worse prognosis.")
    else:
        print(f"\n  Class-level ΔΔG and HR show weak alignment (ρ={rho2:.2f}, n={n2}).")
        print("  Structural destabilisation does not directly predict clinical outcome.")


if __name__ == "__main__":
    main()
