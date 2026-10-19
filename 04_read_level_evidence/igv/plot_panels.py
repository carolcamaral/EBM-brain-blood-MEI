import pysam, glob, os, sys, re
import pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
D, VCFD, OUT = sys.argv[1:4]
MAP = dict(l.split() for l in open(sys.argv[4])) if len(sys.argv) > 4 else {}
REV = {v: k for k, v in MAP.items()}
os.makedirs(OUT, exist_ok=True)
loci = pd.read_csv(f'{D}/loci.tsv', sep='\t', dtype=str)
ev = pd.read_csv(f'{D}/../igv_evidence.csv')
GRAY, CLIP, DISC = '#c9c8c2', '#eb6834', '#2a78d6'
# MELT genotype + GL for each sample at each locus
want = {(r.chrom, int(r.pos)) for r in loci.itertuples()}
gt = {}
for f in ['ALUall.final_comp_original.vcf', 'LINE1.all.final_comp_original.vcf', 'SVAall.final_comp_original.vcf']:
    fam = {'ALU': 'Alu', 'LIN': 'LINE-1', 'SVA': 'SVA'}[f[:3]]
    for l in open(f'{VCFD}/{f}'):
        if l.startswith('#CHROM'): S = [s.split('.')[0] for s in l.rstrip().split('\t')[9:]]; continue
        if l[0] == '#': continue
        t = l.rstrip().split('\t')
        if (t[0], int(t[1])) in want:
            for s, x in zip(S, t[9:]): gt[(fam, t[0], int(t[1]), s)] = x.split(':')[0]
            info = dict(kv.split('=', 1) for kv in t[7].split(';') if '=' in kv)
            gt[(fam, t[0], int(t[1]), 'TSD')] = info.get('TSD', '')
W = 120
for r in loci.itertuples():
    pos = int(r.pos); files = [f'{D}/{r.panel}_{b}_{r.chrom}_{pos}.bam' for b in r.bams.split(',')]
    fig, axes = plt.subplots(len(files), 1, figsize=(9, 2.3 * len(files)), sharex=True, squeeze=False)
    tsd = gt.get((r.family, r.chrom, pos, 'TSD'), '')
    tsdlen = len(tsd) if re.fullmatch('[ACGTN]+', tsd or '') else 0
    for ax, f, b in zip(axes[:, 0], files, r.bams.split(',')):
        bam = pysam.AlignmentFile(f); rows_end = []
        reads = [x for x in bam.fetch(r.chrom, pos - W, pos + W) if not (x.is_unmapped or x.is_duplicate or x.is_secondary or x.is_supplementary)]
        reads.sort(key=lambda x: x.reference_start)
        for x in reads:
            ct = x.cigartuples; lc = ct[0][1] if ct[0][0] == 4 else 0; rc = ct[-1][1] if ct[-1][0] == 4 else 0
            s, e = x.reference_start - lc, x.reference_end + rc
            row = next((i for i, en in enumerate(rows_end) if en < s - 2), None)
            if row is None: rows_end.append(e); row = len(rows_end) - 1
            else: rows_end[row] = e
            disc = x.is_paired and not x.mate_is_unmapped and (x.next_reference_name != r.chrom or abs(x.template_length) > 1500)
            ax.add_patch(Rectangle((x.reference_start, row), x.reference_end - x.reference_start, 0.8, color=DISC if disc else GRAY, lw=0))
            if lc >= 5: ax.add_patch(Rectangle((s, row), lc, 0.8, color=CLIP, lw=0))
            if rc >= 5: ax.add_patch(Rectangle((x.reference_end, row), rc, 0.8, color=CLIP, lw=0))
        ax.set_ylim(len(rows_end) + 0.5, -0.5); ax.set_yticks([])
        ax.axvline(pos, color='#333', lw=0.8, ls='--')
        if tsdlen: ax.axvspan(pos - tsdlen, pos, color='#eda100', alpha=0.25, lw=0)
        e = ev[(ev['panel'] == r.panel) & (ev['sample'] == b)].iloc[0]
        tissue = 'blood' if b.endswith('B') else ('brain' if b[-1] in 'TM' else 'PD study sample')
        g = gt.get((r.family, r.chrom, pos, REV.get(b, b)), '?')
        ax.set_title(f'{b} ({tissue})  -  MELT genotype {g}  |  {e.clipped_reads_pm30} soft-clipped reads at the breakpoint, '
                     f'{e.discordant_pm600} discordant pairs', loc='left', fontsize=9, color='#222')
        for sp in ['top', 'right', 'left']: ax.spines[sp].set_visible(False)
    axes[-1, 0].set_xlim(pos - W, pos + W)
    axes[-1, 0].set_xlabel(f'{r.chrom} position (GRCh38)', fontsize=9)
    axes[-1, 0].ticklabel_format(style='plain', axis='x', useOffset=False)
    fig.suptitle((f'Panel {r.panel}: ' if not os.environ.get('NOPANEL') else '') + f'{r.category} - {r.family}, subject {r.subject} ({r.group}), {r.chrom}:{pos}', x=0.01, ha='left', fontsize=11, fontweight='bold')
    fig.text(0.01, 0.005, 'Grey: aligned read  |  orange: soft-clipped bases  |  blue: read whose mate maps to another chromosome or >1.5 kb away\n'
             'Dashed line: MELT breakpoint; shaded: target-site duplication', fontsize=7.5, color='#555')
    fig.tight_layout(rect=(0, 0.035, 1, 0.97))
    fig.savefig(f'{OUT}/{r.panel}_{r.family}_{r.subject.replace(",", "-")}_{r.chrom}_{pos}.png', dpi=200); plt.close(fig)
print('done')
