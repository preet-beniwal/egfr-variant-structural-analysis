#!/usr/bin/env python3
"""
Phase 2 - Map LUAD EGFR variants from Project 3 to UniProt residue numbers
and compute evolutionary conservation using the ErbB codon alignment from
Project 2.

Inputs:
    Project 3 (TCGA):   ~/EGFR-NSCLC-Clinical-Survival-Analysis/data/mutations_EGFR.tsv
    Project 3 (OncoSG): ~/EGFR-NSCLC-Clinical-Survival-Analysis/data/oncosg/mutations_EGFR.tsv
    Project 2:          ~/erbb_project/results/alignments/erbb_codon_aligned.fasta

Outputs:
    data/variants/luad_egfr_variants_raw.tsv
    data/variants/luad_egfr_variants_mapped.tsv
    data/variants/luad_egfr_variants_with_conservation.tsv

Conservation scoring:
    Codons from the alignment are translated to amino acids, and the score
    is based on amino acid identity across the EGFR vertebrate orthologs
    (human, mouse, chicken, zebrafish). This avoids penalising synonymous
    codon differences.

    Score 3 = all vertebrates share the human amino acid (invariant)
    Score 2 = all but one share (highly conserved)
    Score 1 = at least one but not all share (partially conserved)
    Score 0 = none share (human-specific)
"""

import os
import re
import pandas as pd
from Bio import AlignIO

PROJECT3_TCGA = os.path.expanduser(
    "~/EGFR-NSCLC-Clinical-Survival-Analysis/data/mutations_EGFR.tsv")
PROJECT3_ONCOSG = os.path.expanduser(
    "~/EGFR-NSCLC-Clinical-Survival-Analysis/data/oncosg/mutations_EGFR.tsv")
PROJECT2_ALIGN = os.path.expanduser(
    "~/erbb_project/results/alignments/erbb_codon_aligned.fasta")

OUT_DIR = "data/variants"
os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------- genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(codon):
    """Translate a 3-letter codon to a 1-letter amino acid.
    Returns 'X' for gaps or unknown codons."""
    if not codon or len(codon) != 3 or '-' in codon or 'N' in codon:
        return 'X'
    return CODON_TABLE.get(codon.upper(), 'X')


# ------------------------------------------------------------- HGVS parsing
PAT_SUB = re.compile(r"^([A-Z])(\d+)([A-Z])$")                    # L858R
PAT_DEL = re.compile(r"^([A-Z])(\d+)_([A-Z])(\d+)del$")           # E746_A750del
PAT_DEL_INS = re.compile(
    r"^([A-Z])(\d+)_([A-Z])(\d+)delins([A-Z]+)$")                 # delinsAA
PAT_INS = re.compile(r"^([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$")   # insAA


def parse_protein_change(pc):
    """Return (variant_type, residues_affected, ref_aa, alt_aa)."""
    if pd.isna(pc):
        return None
    pc = str(pc).strip()

    m = PAT_SUB.match(pc)
    if m:
        ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
        return "substitution", [pos], ref, alt

    m = PAT_DEL_INS.match(pc)
    if m:
        start, end = int(m.group(2)), int(m.group(4))
        return "delins", list(range(start, end + 1)), None, m.group(5)

    m = PAT_DEL.match(pc)
    if m:
        start, end = int(m.group(2)), int(m.group(4))
        return "deletion", list(range(start, end + 1)), None, None

    m = PAT_INS.match(pc)
    if m:
        start, end = int(m.group(2)), int(m.group(4))
        return "insertion", list(range(start, end + 1)), None, m.group(5)

    # fallback for variants using alternate syntax like "A750_E758delinsP"
    m = re.match(r"^([A-Z])(\d+)_([A-Z])(\d+)(del|ins)", pc)
    if m:
        start, end = int(m.group(2)), int(m.group(4))
        kind = "deletion" if m.group(5) == "del" else "insertion"
        return kind, list(range(start, end + 1)), None, None

    return "other", None, None, None


# ---------------------------------------------------------- loading helpers
def load_variants(path, cohort_label):
    """Load EGFR mutations from a Project 3 TSV."""
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    print(f"  {cohort_label}: {len(df)} mutation records, "
          f"{df['proteinChange'].nunique()} unique protein changes")
    out = df[["proteinChange", "mutationType"]].copy()
    out["cohort"] = cohort_label
    return out.drop_duplicates(subset=["proteinChange", "cohort"])


def build_conservation_reference(alignment_path):
    """Load the ErbB alignment and build the CDS position -> column map."""
    aln = AlignIO.read(alignment_path, "fasta")
    human_id = None
    for rec in aln:
        if rec.id.startswith("EGFR_human"):
            human_id = rec.id
            break
    if human_id is None:
        raise RuntimeError("EGFR_human not found in alignment")

    seqs = {rec.id: str(rec.seq).upper() for rec in aln}
    human_seq = seqs[human_id]

    cds_to_col = {}
    cds_pos = 0
    for col, base in enumerate(human_seq):
        if base != "-":
            cds_pos += 1
            cds_to_col[cds_pos] = col

    egfr_species = [sid for sid in seqs.keys() if sid.startswith("EGFR_")]
    vertebrate = [s for s in egfr_species if "drosophila" not in s.lower()]

    return {
        "seqs": seqs,
        "human_id": human_id,
        "cds_to_col": cds_to_col,
        "vertebrate_species": vertebrate,
    }


def conservation_for_residue(residue_num, ref):
    """Return species codons, species amino acids, score, alignment columns.

    Score is based on amino acid identity across the vertebrate orthologs.
    """
    cds_start = 3 * (residue_num - 1) + 1

    cols = []
    for i in range(3):
        c = cds_start + i
        if c in ref["cds_to_col"]:
            cols.append(ref["cds_to_col"][c])

    if len(cols) < 3:
        return {}, {}, 0, []

    species_codons = {}
    species_aa = {}
    for sid in ref["vertebrate_species"]:
        codon = "".join(ref["seqs"][sid][c] for c in cols)
        species_codons[sid] = codon
        species_aa[sid] = translate(codon)

    human_aa = species_aa.get(ref["human_id"], "X")

    valid = [s for s in ref["vertebrate_species"]
             if species_aa.get(s, "X") != "X"]
    matches = sum(1 for s in valid if species_aa[s] == human_aa)

    if not valid:
        score = 0
    elif matches == len(valid):
        score = 3
    elif matches >= len(valid) - 1:
        score = 2
    elif matches >= 1:
        score = 1
    else:
        score = 0

    return species_codons, species_aa, score, cols


def short(sid):
    for key in ["human", "mouse", "chicken", "zebrafish", "drosophila"]:
        if key in sid.lower():
            return key[:5]
    return sid[:5]


# ------------------------------------------------------------------- main
def main():
    print("[Phase 2] Mapping LUAD EGFR variants to structure residues\n")

    print("Loading Project 3 variant lists:")
    tcga = load_variants(PROJECT3_TCGA, "TCGA-LUAD")
    onco = load_variants(PROJECT3_ONCOSG, "OncoSG")
    all_vars = pd.concat([tcga, onco], ignore_index=True)
    print(f"  combined: {len(all_vars)} variant records\n")
    all_vars.to_csv(f"{OUT_DIR}/luad_egfr_variants_raw.tsv",
                    sep="\t", index=False)

    print("Parsing protein change strings:")
    parsed = []
    for _, row in all_vars.iterrows():
        p = parse_protein_change(row["proteinChange"])
        if p is None:
            continue
        vtype, residues, ref_aa, alt_aa = p
        parsed.append({
            "protein_change": row["proteinChange"],
            "cohort": row["cohort"],
            "mutation_type": row["mutationType"],
            "variant_type": vtype,
            "residues": ",".join(str(r) for r in residues) if residues else "",
            "first_residue": residues[0] if residues else None,
            "last_residue": residues[-1] if residues else None,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
        })

    mapped = pd.DataFrame(parsed).drop_duplicates(
        subset=["protein_change", "cohort"])
    print(f"  parsed: {len(mapped)} variant-cohort records")
    print("  by type:")
    for t, n in mapped["variant_type"].value_counts().items():
        print(f"    {t:<15} {n}")
    mapped.to_csv(f"{OUT_DIR}/luad_egfr_variants_mapped.tsv",
                  sep="\t", index=False)

    print("\nLoading Project 2 alignment:")
    ref = build_conservation_reference(PROJECT2_ALIGN)
    print(f"  alignment: {len(ref['seqs'])} species, "
          f"{len(ref['cds_to_col'])} CDS positions")
    print(f"  vertebrate species: {ref['vertebrate_species']}")

    print("\nComputing conservation (amino-acid level):")
    results = []
    for _, row in mapped.iterrows():
        if row["variant_type"] not in (
                "substitution", "deletion", "insertion", "delins"):
            continue
        if row["first_residue"] is None:
            continue

        res_num = int(row["first_residue"])
        codons, aa, score, cols = conservation_for_residue(res_num, ref)

        codons_short = {short(k): v for k, v in codons.items()}
        aa_short = {short(k): v for k, v in aa.items()}

        results.append({
            **row.to_dict(),
            "anchor_residue": res_num,
            "human_codon": codons.get(ref["human_id"], ""),
            "human_aa": aa.get(ref["human_id"], ""),
            "codons_by_species": "|".join(
                f"{k}={v}" for k, v in sorted(codons_short.items())),
            "aa_by_species": "|".join(
                f"{k}={v}" for k, v in sorted(aa_short.items())),
            "conservation_score": score,
        })

    final = pd.DataFrame(results)
    out_path = f"{OUT_DIR}/luad_egfr_variants_with_conservation.tsv"
    final.to_csv(out_path, sep="\t", index=False)

    print(f"\n  wrote {out_path}")
    print(f"  total mapped variants: {len(final)}")
    print("\n  conservation score distribution (anchor residue):")
    for s in [3, 2, 1, 0]:
        n = int((final["conservation_score"] == s).sum())
        print(f"    score {s}: {n}")

    print("\n  classic EGFR variants:")
    for pc in ["L858R", "T790M", "E746_A750del", "L747_P753del",
               "G719S", "L861Q", "S768I"]:
        sub = final[final["protein_change"] == pc]
        if len(sub):
            row = sub.iloc[0]
            print(f"    {pc:<20} residue {row['anchor_residue']:<5} "
                  f"codon {row['human_codon']}  aa {row['human_aa']}  "
                  f"score {row['conservation_score']}")

    print("\n  top 10 most conserved variants (score 3):")
    top = final[final["conservation_score"] == 3].head(10)
    for _, row in top.iterrows():
        print(f"    {row['protein_change']:<20} residue {row['anchor_residue']:<5} "
              f"({row['cohort']})")


if __name__ == "__main__":
    main()
