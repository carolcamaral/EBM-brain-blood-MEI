#!/usr/bin/env python3
"""
Add two annotations to a table of MEI calls (TSV or CSV with CHROM/POS or chrom/pos
columns, e.g. sites_all.csv from the notebook). Needs no BAMs, only the reference
FASTA and the catalogue files:

  known_catalogue : name(s) of germline MEI catalogues with a site within
                    --window bp (default 100) of the call
  EN_context      : reference sequence around the insertion site (pos-4..pos+4)
  EN_motif        : True if it matches the L1 endonuclease cleavage motif
                    (5'-TT/AAAA-3' or its reverse complement, 1 mismatch allowed)

Catalogues: any number of --known NAME=PATH, where PATH is a BED (chrom start end)
or a VCF (e.g. the 1000 Genomes MEI prior VCFs shipped in MELT's prior_files/,
dbRIP / euL1db exports, or gnomAD-SV sites filtered to INS:ME). Coordinates
must be GRCh38 (lift over hg19 files first).

Usage (Setonix):
  "$CONDA_PREFIX/bin/python" annotate_known_and_motif.py \
      --calls sites_all.csv --ref GRCh38_full_analysis_set_plus_decoy_hla.fa \
      --known 1KGP_ALU=ALU.1KGP.sites.vcf --known 1KGP_LINE1=LINE1.1KGP.sites.vcf \
      --known 1KGP_SVA=SVA.1KGP.sites.vcf --known gnomADSV=gnomad_sv_ME.bed \
      --out sites_all.annotated.csv
Pure Python: reads the FASTA through its .fai index (no pysam needed).
"""
import argparse, bisect, collections, csv, gzip

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def load_catalogue(path):
    idx = collections.defaultdict(list)
    for line in opener(path):
        if line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        if path.rstrip(".gz").endswith(".vcf"):
            s = e = int(f[1])
        else:
            s, e = int(f[1]), int(f[2])
        c = f[0] if f[0].startswith("chr") else "chr" + f[0]
        idx[c].append((s, e))
    for c in idx:
        idx[c].sort()
    return idx

def near(idx, chrom, pos, w):
    lst = idx.get(chrom, [])
    i = bisect.bisect_left(lst, (pos - w - 10**6, 0))
    for s, e in lst[max(0, i - 1):]:
        if s > pos + w:
            break
        if s - w <= pos <= e + w:
            return True
    return False

class Fasta:
    """Minimal indexed FASTA reader (uses the .fai next to the FASTA)."""
    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        for line in open(path + ".fai"):
            name, length, offset, lbases, lwidth = line.split("\t")[:5]
            self.idx[name] = (int(length), int(offset), int(lbases), int(lwidth))

    def fetch(self, chrom, start, end):          # 0-based, end-exclusive, like pysam
        length, offset, lb, lw = self.idx[chrom]
        start, end = max(0, start), min(end, length)
        if start >= end:
            return ""
        a = offset + (start // lb) * lw + start % lb
        b = offset + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        self.fh.seek(a)
        return self.fh.read(b - a).decode().replace("\n", "").replace("\r", "")


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def motif_ok(seq):
    for target in ("TTAAAA", "TTTTAA"):
        for i in range(len(seq) - 5):
            if sum(a != b for a, b in zip(seq[i:i + 6], target)) <= 1:
                return True
    return False

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--known", action="append", default=[], help="NAME=PATH")
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    fa = Fasta(a.ref)
    cats = {}
    for kv in a.known:
        name, path = kv.split("=", 1)
        cats[name] = load_catalogue(path)
    delim = "," if a.calls.endswith(".csv") else "\t"
    rows = list(csv.DictReader(open(a.calls, newline="", encoding="utf-8"), delimiter=delim))
    ck = "CHROM" if "CHROM" in rows[0] else "chrom"
    pk = "POS" if "POS" in rows[0] else "pos"
    for r in rows:
        chrom, pos = r[ck], int(r[pk])
        hits = [n for n, idx in cats.items() if near(idx, chrom, pos, a.window)]
        r["known_catalogue"] = ",".join(hits)
        try:
            seq = fa.fetch(chrom, max(0, pos - 5), pos + 4).upper()
        except (KeyError, ValueError):
            seq = ""
        r["EN_context"] = seq
        r["EN_motif"] = bool(seq) and (motif_ok(seq) or motif_ok(rc(seq)))
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter=delim)
        w.writeheader(); w.writerows(rows)
    n_known = sum(1 for r in rows if r["known_catalogue"])
    print(f"{len(rows)} calls; {n_known} near a known germline MEI; "
          f"{sum(1 for r in rows if r['EN_motif'])} with an EN motif")

if __name__ == "__main__":
    main()
