#!/usr/bin/env python3
"""Count read evidence for an insertion at every MELT site, from `samtools view` (SAM) on stdin.
Usage: samtools view -F 0xF04 -M -L sites_scan.bed X.bam | python3 scan_evidence.py sites_scan.tsv SAMPLE > SAMPLE.evidence.tsv
Per site: coverage at the breakpoint; reads soft-clipped (>=15 bp) with the clip point within 30 bp of the site;
how many of those clips are poly-A/T (>=10 nt); discordant reads (MAPQ>=20, mate on another chromosome or >1.5 kb away)
starting within 300 bp of the site. Pure Python, no pysam."""
import sys, re, bisect
sites, sample = sys.argv[1], sys.argv[2]
S = {}
for i, l in enumerate(open(sites)):
    if i == 0: continue
    sid, fam, chrom, pos, tsd = l.rstrip('\n').split('\t')
    S.setdefault(chrom, []).append((int(pos), sid))
for c in S: S[c].sort()
P = {c: [p for p, _ in v] for c, v in S.items()}
cnt = {sid: [0, 0, 0, 0] for v in S.values() for _, sid in v}   # cov, clip, clip_polyA, disc
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
PA = re.compile('A{10,}|T{10,}')
seen = set()
for line in sys.stdin:
    t = line.split('\t', 11)
    flag, chrom, start, mapq, cigar, rnext, tlen, seq = int(t[1]), t[2], int(t[3]) - 1, int(t[4]), t[5], t[6], int(t[8]), t[9]
    if chrom not in P or cigar == '*': continue
    key = (t[0], flag & 0xC0)               # -M can return a read twice for overlapping regions
    if key in seen: continue
    seen.add(key)
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    reflen = sum(n for n, o in ops if o in 'MDN=X')
    end = start + reflen
    # nearest site(s) within 350 bp
    j = bisect.bisect_left(P[chrom], start - 350)
    while j < len(P[chrom]) and P[chrom][j] <= end + 350:
        pos, sid = S[chrom][j]; c = cnt[sid]; j += 1
        if start <= pos < end: c[0] += 1
        if ops[0][1] == 'S' and ops[0][0] >= 15 and abs(start - pos) <= 30:
            c[1] += 1; c[2] += bool(PA.search(seq[:ops[0][0]]))
        if ops[-1][1] == 'S' and ops[-1][0] >= 15 and abs(end - pos) <= 30:
            c[1] += 1; c[2] += bool(PA.search(seq[-ops[-1][0]:]))
        if (flag & 1) and not (flag & 8) and mapq >= 20 and abs(start - pos) <= 300 and (rnext != '=' or abs(tlen) > 1500):
            c[3] += 1
print('site_id\tsample\tcoverage\tclipped\tclipped_polyA\tdiscordant')
for sid, (a, b, cpa, d) in cnt.items():
    print(f'{sid}\t{sample}\t{a}\t{b}\t{cpa}\t{d}')
