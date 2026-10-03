# EGFR Variant Structural Analysis

**Structural impact of clinical EGFR variants: a three-layer analysis combining
evolutionary conservation, protein stability, and clinical outcome.**

## Motivation

Clinical interpretation of EGFR variants in lung adenocarcinoma typically relies
on population frequency and clinical annotations. Two additional layers of
evidence are available but rarely integrated:

1. Evolutionary conservation -- Project 2 (erbb-phylogenomics) produced
   conservation scores for every ErbB alignment column.
2. Clinical outcome -- Project 3 (EGFR-NSCLC-Clinical-Survival-Analysis)
   produced hazard ratios for EGFR variant classes in TCGA-LUAD and OncoSG.

Neither analysis asked where the variants physically sit in the protein, or
whether structural destabilisation correlates with either layer. This project
adds the structural dimension.

## Pipeline

| Phase | Purpose |
|-------|---------|
| 1 | Structure acquisition (AlphaFold + PDB) |
| 2 | Variant-to-residue mapping |
| 3 | In-silico mutagenesis (PyMOL + FoldX) |
| 4 | Three-layer analysis (conservation vs ddG vs HR) |
| 5 | PyMOL figure generation |
| 6 | Manuscript writeup |

Status: Phase 1 complete.

## Author

**Preet Beniwal** -- Independent bioinformatics portfolio project.
GitHub: [@preet-beniwal](https://github.com/preet-beniwal)

## License

MIT
