#!/bin/bash -l
#SBATCH --job-name=igv_slices
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=pawsey1276
#SBATCH --partition=copy
#SBATCH --export=NONE
#SBATCH --array=1-8%4
#SBATCH --output=/scratch/pawsey1276/cdoamaral/igv/logs/slice-%A_%a.log
#SBATCH --error=/scratch/pawsey1276/cdoamaral/igv/logs/slice-%A_%a.err
#SBATCH --time=24:00:00
# One array task per BAM (bams.txt, 8 lines): pshell get -> cut +-1 kb around each locus in loci.tsv -> delete the full BAM.
# Max 4 at a time (%4), so at most 4 full BAMs sit on /scratch.
#   mkdir -p /scratch/pawsey1276/cdoamaral/igv/logs   (put loci.tsv, bams.txt and this script there)
#   cd /scratch/pawsey1276/cdoamaral/igv && sbatch get_and_slice.sh
#   (to keep the full BAMs too: sbatch --export=NONE,KEEP_BAMS=1 get_and_slice.sh)
#   when all 8 finish:  bash get_and_slice.sh pack

WORK=/scratch/pawsey1276/cdoamaral/igv
PSHELL=/software/projects/pawsey0360/cdoamaral/pshell
KEEP_BAMS=${KEEP_BAMS:-0}                         # 1 = keep the full BAMs in $WORK/bams (needs ~1-1.5 TB)
REMOTE=/projects/MNDgenomics/BloodBrainWGS/bam     # absolute path on the Data portal
FLANK=1000
cd $WORK
OUT=$WORK/igv_slices; mkdir -p $OUT/snapshots

if [[ "$1" == "pack" ]]; then
  cp loci.tsv $OUT/ && tar czf igv_slices.tar.gz igv_slices
  echo "$(ls $OUT/*.bam | wc -l) slices (expected 20) -> $WORK/igv_slices.tar.gz ($(du -h igv_slices.tar.gz | cut -f1))"
  exit 0
fi

module load python/3.11.6
module load pawseyenv/2025.08 2>/dev/null || true
module load singularity/4.1.0-slurm 2>/dev/null || module load singularity 2>/dev/null || true
module load samtools/1.15--h3843a85_0
samtools --version >/dev/null || { echo "samtools not working on this partition"; exit 1; }

b=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bams.txt)
f=${b}.bwa.sorted.marked.bam
TMP=$WORK/tmp_$b; mkdir -p $TMP

RETRIES=5
for i in $(seq 1 $RETRIES); do
    python3 $PSHELL "lcd $TMP && cd $REMOTE && get $f && get $f.bai" && break
    echo "Attempt $i failed, retrying in 30 seconds..."
    sleep 30
done
[[ -s $TMP/$f && -s $TMP/$f.bai ]] || { echo "download failed for $b"; exit 1; }

# every locus that uses this BAM
awk -F'\t' -v b="$b" 'NR>1{n=split($8,a,","); for(i=1;i<=n;i++) if(a[i]==b) print $1"\t"$6"\t"$7}' loci.tsv |
while IFS=$'\t' read -r panel chrom pos; do
  out=$OUT/${panel}_${b}_${chrom}_${pos}.bam
  samtools view -b -o $out $TMP/$f "${chrom}:$((pos-FLANK))-$((pos+FLANK))"
  samtools index $out
  echo "$panel $b $chrom:$pos $(samtools view -c $out) reads"
done

if [[ $KEEP_BAMS == 1 ]]; then mkdir -p $WORK/bams && mv $TMP/$f $TMP/$f.bai $WORK/bams/ ; fi
rm -rf $TMP
