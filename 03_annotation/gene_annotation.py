"""Gene symbol and genomic location for each MEI site, from the MELT INTERNAL field (RefSeq accession, location).
RefSeq accessions are converted to gene symbols with ANNOVAR's refGene table (columns 2 and 13), e.g.:
    cut -f2,13 humandb/hg19_refGene.txt | sort -u > refseq_to_gene.tsv
(the accession-to-symbol mapping does not depend on the genome build).

Usage: python gene_annotation.py sites_all.csv calls_brain_only.csv calls_shared.csv refseq_to_gene.tsv out_dir
"""
import sys, os
import pandas as pd

sites_csv, brain_csv, shared_csv, mapping, out = sys.argv[1:6]
FILTER = 'kept: PASS + SR>2 + ASSESS≥3'
os.makedirs(out, exist_ok=True)
sym = pd.read_csv(mapping, sep='\t', header=None, names=['acc', 'gene']).drop_duplicates('acc').set_index('acc').gene

def parse(v):
    v = str(v)
    if v in ('nan', 'null,null'):
        return [], []
    accs, regs = [], []
    for p in v.split('|'):
        a, _, r = p.partition(',')
        accs.append(a); regs.append(r)
    return accs, regs

def location(regs):
    if not regs:
        return 'intergenic'
    for key, name in [('EXON', 'exon'), ('5_UTR', "5' UTR"), ('3_UTR', "3' UTR"), ('PROMOTER', 'promoter'), ('TERMINATOR', 'terminator'), ('INTRONIC', 'intron')]:
        if any(r.startswith(key) for r in regs):
            return name
    return regs[0]

s = pd.read_csv(sites_csv, low_memory=False).set_index('site_id')
parsed = s.INTERNAL.map(parse)
s['gene'] = parsed.map(lambda p: ';'.join(sorted({sym.get(a, a) for a in p[0]})))
s['location'] = parsed.map(lambda p: location(p[1]))

summary = {}
for name, path in [('brain_specific', brain_csv), ('shared', shared_csv)]:
    calls = pd.read_csv(path, low_memory=False)
    calls = calls[calls[FILTER]]
    site_ids = calls.site_id.unique()
    s.loc[site_ids, ['family', 'CHROM', 'POS', 'gene', 'location']].to_csv(f'{out}/{name}_sites_genes.csv')
    summary[name] = s.loc[site_ids, 'location'].value_counts(normalize=True).mul(100).round(1)
print(pd.DataFrame(summary).fillna(0))
