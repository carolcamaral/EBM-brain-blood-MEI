#!/bin/bash -l
#SBATCH --job-name=bam_scan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=pawsey1276
#SBATCH --partition=copy
#SBATCH --export=NONE
#SBATCH --array=1-50%4
#SBATCH --output=/scratch/pawsey1276/cdoamaral/scan/logs/scan-%A_%a.log
#SBATCH --error=/scratch/pawsey1276/cdoamaral/scan/logs/scan-%A_%a.err
#SBATCH --time=24:00:00
# One task per BAM in bams.txt (50, BWA-MEM2 alignments in BloodBrainWGS/BWAMEM2array). For each BAM:
#   1. pshell get (BAM + .bai)
#   2. read evidence at all 3,482 MELT sites (PASS + SR>2 + ASSESS>=3)  -> evidence/<bam>.evidence.tsv
#   3. +-1 kb slices for the IGV panels in loci.tsv                     -> igv_slices/
#   4. samtools stats for Table S2 (library QC)                          -> stats/<bam>.stats
#   5. delete the BAM
#   mkdir -p /scratch/pawsey1276/cdoamaral/scan/logs ; put this script, bams.txt, loci.tsv, sites_scan.tsv,
#   sites_scan.bed and scan_evidence.py there; then: cd /scratch/pawsey1276/cdoamaral/scan && sbatch scan_job.sh
#   When all 50 finish:  bash scan_job.sh pack   -> scan_results.tar.gz (a few MB; send it back)

WORK=/scratch/pawsey1276/cdoamaral/scan
PSHELL=/software/projects/pawsey0360/cdoamaral/pshell
REMOTE=/projects/MNDgenomics/BloodBrainWGS/BWAMEM2array
FLANK=1000
cd $WORK
mkdir -p evidence igv_slices stats

if [[ "$1" == "pack" ]]; then
  echo "evidence: $(ls evidence/*.tsv 2>/dev/null | wc -l)/50   stats: $(ls stats/*.stats 2>/dev/null | wc -l)/50   IGV slices: $(ls igv_slices/*.bam 2>/dev/null | wc -l)"
  grep -l "FAILED\|download failed" logs/*.log 2>/dev/null
  cp loci.tsv igv_slices/ && tar czf scan_results.tar.gz evidence stats igv_slices bams.txt && ls -lh scan_results.tar.gz
  exit 0
fi

module load pawseyenv/2025.08 2>/dev/null || true
module load samtools/1.15--h3843a85_0
module load python/3.11.6                           # load last, so pawseyenv does not unload it
samtools --version >/dev/null || { echo "samtools not working"; exit 1; }
PY=$(command -v python3); echo "python: $PY $($PY --version 2>&1)"

b=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.txt)
f=${b}.sorted.bam
ev_ok() { [[ -s evidence/$b.evidence.tsv ]] && awk -F'\t' 'NR>1{c+=$3} END{exit !(c>0)}' evidence/$b.evidence.tsv; }
if ev_ok && [[ -s stats/$b.stats ]]; then echo "skip $b (done)"; exit 0; fi
TMP=$WORK/tmp_$b; mkdir -p $TMP

for i in 1 2 3 4 5; do
    python3 $PSHELL "lcd $TMP && cd $REMOTE && get $f && get $f.bai" && break
    echo "Attempt $i failed, retrying in 30 seconds..."; sleep 30
done
[[ -s $TMP/$f && -s $TMP/$f.bai ]] || { echo "download failed for $b"; exit 1; }
touch $TMP/$f.bai                                   # index must not look older than the BAM
samtools view -H $TMP/$f | grep -m3 "^@SQ"          # check reference naming (expects chr1 ...)

# 2. evidence at all sites. samtools 1.15 can fail with "Iterator could not be created" on -M -L,
#    so try three ways and keep the first that gives non-zero coverage.
rm -f evidence/$b.evidence.tsv
for mode in M regions stream; do
  case $mode in
    M)       samtools view -@ 3 -F 0xF04 -M -L sites_scan.bed $TMP/$f ;;
    regions) samtools view -@ 3 -F 0xF04 $TMP/$f $(awk '{print $1":"$2+1"-"$3}' sites_scan.bed) ;;
    stream)  samtools view -@ 3 -h -F 0xF04 $TMP/$f | samtools view -L sites_scan.bed - ;;
  esac | $PY scan_evidence.py sites_scan.tsv $b > evidence/$b.evidence.tsv.part
  mv evidence/$b.evidence.tsv.part evidence/$b.evidence.tsv
  if ev_ok; then echo "evidence via $mode"; break; else echo "evidence via $mode gave no reads"; rm -f evidence/$b.evidence.tsv; fi
done
[[ -s evidence/$b.evidence.tsv ]] || echo "FAILED evidence $b"
echo "evidence: $(awk -F'\t' 'NR>1 && $4>=3' evidence/$b.evidence.tsv 2>/dev/null | wc -l) sites with >=3 clipped reads"

# 3. IGV slices
awk -F'\t' -v b="$b" 'NR>1{n=split($8,a,","); for(i=1;i<=n;i++) if(a[i]==b) print $1"\t"$6"\t"$7}' loci.tsv |
while IFS=$'\t' read -r panel chrom pos; do
  out=igv_slices/${panel}_${b}_${chrom}_${pos}.bam
  samtools view -b -o $out $TMP/$f "${chrom}:$((pos-FLANK))-$((pos+FLANK))" && samtools index $out
  echo "IGV $panel $b $chrom:$pos $(samtools view -c $out) reads"
done

# 4. samtools stats (whole BAM)
[[ -s stats/$b.stats ]] || samtools stats -@ 3 $TMP/$f > stats/$b.stats.part && mv stats/$b.stats.part stats/$b.stats 2>/dev/null; [[ -s stats/$b.stats ]] || echo "FAILED stats $b"
grep "^SN" stats/$b.stats | grep -E "raw total sequences|reads duplicated|insert size average|average length|error rate|pairs on different chromosomes" 

rm -rf $TMP
echo "finished $b"
