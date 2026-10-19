import pysam, glob, os, re, sys, collections
import pandas as pd
D=sys.argv[1]
loci=pd.read_csv(f'{D}/loci.tsv',sep='\t',dtype=str)
def polyA(s): return bool(re.search('A{10,}|T{10,}',s))
rows=[]
for f in sorted(glob.glob(f'{D}/*.bam')):
    b=os.path.basename(f)[:-4]; panel=b[0]; chrom,pos=b.split('_')[-2],int(b.split('_')[-1])
    sample=b[2:].rsplit('_',2)[0]
    bam=pysam.AlignmentFile(f)
    clips=collections.Counter(); clipreads=[]; disc=[];cov=0
    for r in bam.fetch(chrom,pos-700,pos+700):
        if r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail: continue
        if r.reference_start<=pos<r.reference_end: cov+=1
        ct=r.cigartuples
        for side,(op,ln) in (('L',ct[0]),('R',ct[-1])):
            if op==4 and ln>=15:
                bp=r.reference_start if side=='L' else r.reference_end
                if abs(bp-pos)<=30:
                    seq=r.query_sequence[:ln] if side=='L' else r.query_sequence[-ln:]
                    clipreads.append((bp,side,ln,polyA(seq),r.has_tag('SA')))
        if r.is_paired and not r.mate_is_unmapped and abs(r.reference_start-pos)<=600:
            if r.next_reference_name!=chrom or abs(r.template_length)>1500:
                disc.append((r.reference_start<pos, r.is_reverse, r.mapping_quality, r.next_reference_name))
    # dominant clip cluster
    cc=collections.Counter((bp,s) for bp,s,*_ in clipreads)
    top=cc.most_common(2)
    rows.append(dict(panel=panel,sample=sample,locus=f'{chrom}:{pos}',coverage=cov,
        clipped_reads_pm30=len(clipreads),clip_polyA=sum(c[3] for c in clipreads),clip_SA=sum(c[4] for c in clipreads),
        top_clip_sites='; '.join(f'{bp-pos:+d}{s}x{n}' for (bp,s),n in top),
        discordant_pm600=len(disc),disc_left_fwd=sum(1 for l,rv,*_ in disc if l and not rv),disc_right_rev=sum(1 for l,rv,*_ in disc if not l and rv),
        disc_mate_chroms=','.join(f'{k}:{v}' for k,v in collections.Counter(d[3] for d in disc).most_common(3))))
df=pd.DataFrame(rows)
pd.set_option('display.width',250)
print(df.to_string(index=False))
df.to_csv(f'{D}/../igv_evidence.csv',index=False)
