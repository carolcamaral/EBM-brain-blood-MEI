# Supplementary tables

The supplementary tables (S1–S4) were assembled from the outputs of the other folders:

| Table | Source |
|---|---|
| S1 – MEI calls per sample | `01_melt_vcf_analysis` (`counts_per_person.csv`) |
| S2 – Sequencing metrics | `samtools stats` output written by `04_read_level_evidence/scan_job.sh` |
| S3 – Brain-specific calls | `01_melt_vcf_analysis` (`tableS3_brain_only_calls.csv`) with `03_annotation` (gnomAD-SV overlap, L1 endonuclease motif, genes) |
| S4 – Filter comparison | `01_melt_vcf_analysis` (`filter_comparison.csv`) and `02_statistics_and_figures` |
