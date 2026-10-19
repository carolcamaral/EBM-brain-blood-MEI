# Brain–blood mobile element insertions in healthy aging and neurodegeneration

Code for the analyses in:

> Kõks S, Pfaff A, Amaral C. *Brains of healthy elderly controls show more mobile element insertion calls than matched blood.* Experimental Biology and Medicine (under revision).

Paired brain and blood whole-genome sequencing from 23 individuals (9 healthy elderly controls, 5 Alzheimer’s disease, 9 motor neuron disease). Non-reference Alu, LINE-1 and SVA insertions were called with MELT v2.1.5 in multiple-sample mode; this repository starts from the MELT VCFs.

Individual-level genomic data and patient information are **not** included (see `data/README.md`).

## Contents

| Folder | What it does | Paper |
|---|---|---|
| `01_melt_vcf_analysis/` | Notebook that reads the three MELT VCFs, applies the filters (FILTER = PASS, SR > 2, ASSESS ≥ 3), counts calls per sample, classifies each call as shared, brain-specific or blood-specific, compares filter settings and writes per-call tables. `nb_source.py` is the notebook source; `build_nb.py` rebuilds and executes it. | Tables 1–4, S1, S3, S4 |
| `02_statistics_and_figures/` | Paired Wilcoxon signed-rank tests (brain vs blood) with Bonferroni correction, group comparisons, age analyses (Spearman, linear model) and Figures 1–3. | Tables 2–3, Figures 1–3, age results |
| `03_annotation/` | Overlap with gnomAD-SV v4.1 mobile element insertions (≤ 100 bp) and L1 endonuclease motif at the insertion site (`annotate_known_and_motif.py`, run with `run_annotation.sh` on Pawsey Setonix); chance rate of the motif (`motif_chance.py`); gene and genomic location from the MELT RefSeq annotation (`gene_annotation.py`). | Table 4, S3, genomic distribution |
| `04_read_level_evidence/` | Slurm array job that downloads each BAM from the Pawsey Data Portal, counts soft-clipped and discordant reads at every MELT site (`scan_evidence.py`), cuts BAM slices for genome browser views and runs `samtools stats`. `igv/` draws the read views from the slices. | Table S2, Supplementary Figures S1–S2 |
| `05_supplementary_tables/` | How the supplementary tables were assembled. | S1–S4 |

## Running

```bash
pip install -r requirements.txt

# 1. MELT VCFs -> per-call tables and counts (expects the inputs in data/)
cd 01_melt_vcf_analysis && python build_nb.py brain_blood_MEI_original_VCFs.ipynb

# 2. statistics and figures
cd ../02_statistics_and_figures && EBM_DATA=../data python compute.py && python figures.py

# 3. annotation (on an HPC login node with the GRCh38 FASTA)
cd ../03_annotation && bash run_annotation.sh && python motif_chance.py

# 4. read-level evidence (Pawsey Setonix; edit account and paths at the top of the script)
cd ../04_read_level_evidence && sbatch scan_job.sh
```

## Software

MELT v2.1.5; samtools 1.15; Python 3 with pandas, NumPy, SciPy, statsmodels, matplotlib and pysam.

## License

MIT (see `LICENSE`).
