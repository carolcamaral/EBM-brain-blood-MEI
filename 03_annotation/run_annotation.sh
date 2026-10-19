#!/bin/bash
# Annotate MEI calls with (1) overlap with germline MEI catalogues and (2) the L1 endonuclease motif.
# No BAMs needed. Inputs: a calls table from the notebook, the GRCh38 FASTA, and catalogue files.
# Run on the Setonix LOGIN node (the downloads need internet; the annotation itself takes about a minute).
#
#   bash run_annotation.sh
set -euo pipefail

WORK=/scratch/pawsey1276/$USER/ebm_annotation          # edit if you prefer another folder
# Any GRCh38 FASTA with "chr" names works (only ~9 bp around each site are read). Override with: REF=/path/x.fa bash run_annotation.sh
REF=${REF:-/software/projects/pawsey0360/apfaff/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa}
CALLS=$WORK/sites_all.csv                               # copy this file from outputs/ in the notebook zip
mkdir -p "$WORK/catalogues" && cd "$WORK"

[[ -f "$CALLS" ]] || { echo "Copy sites_all.csv (from brain_blood_outputs.zip) to $WORK first"; exit 1; }
[[ -f "$REF" ]] || { echo "Reference FASTA not found: $REF  (set REF=/path/to/GRCh38.fa)"; exit 1; }
if [[ ! -f "$REF.fai" ]]; then                          # index it if only the .fai is missing
  echo "No $REF.fai - indexing (a few minutes)"
  module load samtools 2>/dev/null || true
  if command -v samtools >/dev/null; then samtools faidx "$REF"
  else echo "samtools not found: run 'module load samtools' and retry"; exit 1; fi
fi

# ---- 1. gnomAD-SV v4.1 (GRCh38): keep mobile element insertions only -------------------------
# ~2 GB download; check the URL on https://gnomad.broadinstitute.org/downloads#v4-structural-variants if it fails
G=catalogues/gnomad_sv_ME.bed
if [[ ! -s $G ]]; then
  curl -L -o catalogues/gnomad.v4.1.sv.sites.vcf.gz \
    https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
  zcat catalogues/gnomad.v4.1.sv.sites.vcf.gz | awk -F'\t' '!/^#/ && $5 ~ /INS:ME/ {print $1"\t"$2"\t"$2}' > $G
  echo "gnomAD-SV mobile element insertions: $(wc -l < $G)"
fi

# ---- 2. 1000 Genomes MEI sites shipped with MELT (prior_files/) -------------------------------
# Set MELT_DIR to your MELT install. These files are hg19 in older MELT releases: if the header says
# GRCh37/hg19, lift them over to GRCh38 (e.g. with CrossMap) before use, or skip them.
MELT_DIR=${MELT_DIR:-}
KNOWN=(--known gnomADSV="$G")
if [[ -n "$MELT_DIR" && -d "$MELT_DIR/prior_files" ]]; then
  for f in "$MELT_DIR"/prior_files/*.vcf; do
    n=$(basename "$f" .vcf); KNOWN+=(--known "1KGP_${n}=$f")
  done
fi

# ---- 3. (optional) dbRIP / euL1db: download from their websites, lift to GRCh38, save as BED and add:
#   KNOWN+=(--known dbRIP=catalogues/dbRIP_hg38.bed --known euL1db=catalogues/euL1db_hg38.bed)

# ---- 4. annotate -------------------------------------------------------------------------------
python3 annotate_known_and_motif.py --calls "$CALLS" --ref "$REF" "${KNOWN[@]}" \
        --window 100 --out sites_all.annotated.csv
echo "done: $WORK/sites_all.annotated.csv  (send this file back)"
