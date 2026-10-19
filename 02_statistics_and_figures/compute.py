import os
DATA = os.environ.get('EBM_DATA', '../data')  # folder with the patient information sheet and the notebook outputs
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, spearmanr
from stats_lib import oneway_tissue, posthoc_bonf
from scipy.stats import wilcoxon
def pw(x, f):
    w = x.pivot_table(index='subject', columns='tissue', values=f)
    return wilcoxon(w.brain, w.blood).pvalue
F = 'PASS + SR>2 + ASSESS≥3'
c = pd.read_csv(os.path.join(DATA, 'outputs', 'counts_per_person.csv'))
c = c[(c['filter'] == F)]
fams = {'Alu': 'Alu', 'LINE-1': 'LINE1', 'SVA': 'SVA'}
# long sample table like the original RIPs csv
rows = []
for _, r in c[c.family != 'All TEs'].iterrows():
    for t in ['blood', 'brain']:
        rows.append(dict(subject=r.subject, dgn=r.group.lower(), tissue=t, fam=fams[r.family], n=r[t]))
L = pd.DataFrame(rows).pivot_table(index=['subject', 'dgn', 'tissue'], columns='fam', values='n').reset_index()
L['dgn2'] = L.dgn.map({'hc': 'hc', 'ad': 'nd', 'mnd': 'nd'})
L.to_csv('pass_counts_long.csv', index=False)
R = {}
# Table 2
t2 = []
for f in ['Alu', 'LINE1', 'SVA']:
    b, r = L[L.tissue == 'blood'][f], L[L.tissue == 'brain'][f]
    p = pw(L, f)
    t2.append(dict(fam=f, blood=f"{b.mean():.0f} ({b.std():.0f})", brain=f"{r.mean():.0f} ({r.std():.0f})", p=p, padj=min(1, 3 * p)))
R['table2'] = t2
# Table 3
t3 = []
for f in ['Alu', 'LINE1', 'SVA']:
    ph = posthoc_bonf(L, f, 'dgn'); ph2 = posthoc_bonf(L, f, 'dgn2')
    for g, lab, n in [('hc', 'HC', 9), ('mnd', 'MND', 9), ('ad', 'AD', 5), ('nd', 'ND', 14)]:
        col = 'dgn2' if g == 'nd' else 'dgn'
        x = L[L[col] == g]; b, r = x[x.tissue == 'blood'][f], x[x.tissue == 'brain'][f]
        p = pw(x, f)
        t3.append(dict(fam=f, dgn=f"{lab}, n={n}", blood=f"{b.mean():.0f} ({b.std():.0f})", brain=f"{r.mean():.0f} ({r.std():.0f})", p=p, padj=min(1, 3 * p)))
R['table3'] = t3
# per-subject deltas
W = L.pivot_table(index=['subject', 'dgn'], columns='tissue', values=['Alu', 'LINE1', 'SVA'])
D = pd.DataFrame({f: (W[(f, 'brain')] - W[(f, 'blood')]).astype(int) for f in ['Alu', 'LINE1', 'SVA']}).reset_index()
D['All'] = D[['Alu', 'LINE1', 'SVA']].sum(axis=1)
info = pd.read_excel(os.path.join(DATA, 'RIPs (Brain_Blood) - Patient Info.xlsx'), sheet_name='SUBJECT INFO').dropna(subset=['subject'])
info['subject'] = info.subject.astype(int)
D = D.merge(info[['subject', 'age at death', 'sex']].rename(columns={'age at death': 'age'}), on='subject')
D['HC'] = (D.dgn == 'hc').astype(int)
R['table1'] = D[['subject', 'dgn', 'Alu', 'LINE1', 'SVA']].to_dict('records')
# Results text: averages
R['avg_blood'] = {f: round(L[L.tissue == 'blood'][f].mean()) for f in ['Alu', 'LINE1', 'SVA']}
R['avg_brain'] = {f: round(L[L.tissue == 'brain'][f].mean()) for f in ['Alu', 'LINE1', 'SVA']}
# age analyses
hc, dz = D[D.HC == 1], D[D.HC == 0]
age = {}
for f in ['All', 'Alu', 'LINE1', 'SVA']:
    m = smf.ols(f'{f} ~ HC + age', D).fit()
    age[f] = dict(med_hc=float(hc[f].median()), med_dz=float(dz[f].median()), mwu=mannwhitneyu(hc[f], dz[f]).pvalue,
                  eff=m.params.HC, p_eff=m.pvalues.HC, p_age=m.pvalues.age,
                  rho_all=spearmanr(D.age, D[f])[0], p_all=spearmanr(D.age, D[f])[1])
    hcad = D[D.dgn.isin(['hc', 'ad'])]
    age[f].update(rho_hcad=spearmanr(hcad.age, hcad[f])[0], p_hcad=spearmanr(hcad.age, hcad[f])[1])
age['rho_hc_only_Alu'] = list(map(float, spearmanr(hc.age, hc.Alu)))
R['age'] = age
# HC vs AD
ad, mnd = D[D.dgn == 'ad'], D[D.dgn == 'mnd']
R['hc_vs_ad'] = {f: dict(med_hc=float(hc[f].median()), med_ad=float(ad[f].median()), p=mannwhitneyu(hc[f], ad[f]).pvalue) for f in ['Alu', 'LINE1', 'SVA', 'All']}
hcad = D[D.dgn.isin(['hc', 'ad'])].assign(AD=lambda x: (x.dgn == 'ad').astype(int))
m = smf.ols('Alu ~ AD + age', hcad).fit()
R['ad_adjusted'] = dict(coef=m.params.AD, se=m.bse.AD, p=m.pvalues.AD)
R['hc_vs_mnd_Alu_p'] = mannwhitneyu(hc.Alu, mnd.Alu).pvalue
# flowcell: different flowcell for all HC & AD, same for MND except 18614
D['diff_fc'] = (D.dgn != 'mnd') | (D.subject == 18614)
R['flowcell'] = dict(med_diff=float(D[D.diff_fc].Alu.median()), med_same=float(D[~D.diff_fc].Alu.median()),
                     p=mannwhitneyu(D[D.diff_fc].Alu, D[~D.diff_fc].Alu).pvalue, n_same=int((~D.diff_fc).sum()))
# yield (Gb) from the AGRF report
Y = {'A072_02B':168.53,'A072_02M':119.41,'A186_14M':115.26,'A216_99B':143.68,'A216_99M':129.94,'A319_95B':143.69,'A319_95M':114.95,
'0078_B':149.13,'0180_B':154.32,'0180_M':129.11,'0236_B':117.60,'0236_M':135.41,'0561_B':100.52,'0561_M':150.57,'0078_M':58.83,
'4011_B':184.03,'4011_M':142.84,'14_46_T':150.03,'10_07_T':131.33,'11_06_T':124.62,'11_07_T':153.48,'11_22_T':152.78,'14_04_T':133.74,
'15_28_T':131.87,'16_09_T':138.31,'09_24_T':149.18,'09_26_T':140.20,'14_46_B':128.29,'15_01_T':130.71,'15_11_T':145.23,'16_03_T':134.11,
'16_13_T':130.74,'10_07_B':139.41,'11_06_B':134.02,'11_07_B':142.71,'11_22_B':148.82,'14_04_B':122.11,'15_28_B':118.23,'16_09_B':172.95,
'09_24_B':111.96,'09_26_B':124.29,'15_01_B':153.98,'15_11_B':157.85,'16_03_B':142.46,'16_13_B':150.14,'A186_14B':135.43}
sheet = pd.read_excel(os.path.join(DATA, 'RIPs (Brain_Blood) - Patient Info.xlsx'), sheet_name='RIPsPerronProjectwoPD').dropna(subset=['sampleID'])
sheet.columns = [x.strip() for x in sheet.columns]
def yld(s):
    k = s if s in Y else s.replace('_M', 'M').replace('M', '_M') if False else None
    for cand in [s, s[:-1] + '_' + s[-1], s.replace('_', '', 0)]:
        if cand in Y: return Y[cand]
    raise KeyError(s)
sheet['yield'] = sheet.sampleID.map(yld)
yy = sheet.pivot_table(index='subject', columns='tissue', values='yield')
D = D.merge((yy.brain - yy.blood).rename('dyield'), left_on='subject', right_index=True)
hcad = D[D.dgn.isin(['hc', 'ad'])]
R['yield_corr_hcad'] = list(map(float, spearmanr(hcad.Alu, hcad.dyield)))
# ranges per group for the text
R['ranges'] = {g: {f: [int(D[D.dgn == g][f].min()), int(D[D.dgn == g][f].max())] for f in ['Alu', 'LINE1', 'SVA']} for g in ['hc', 'ad', 'mnd']}
D.to_csv('pass_deltas.csv', index=False)
json.dump(R, open('pass_numbers.json', 'w'), indent=1, default=float)
print(json.dumps({k: R[k] for k in ['table2', 'avg_blood', 'avg_brain', 'hc_vs_ad', 'ad_adjusted', 'hc_vs_mnd_Alu_p', 'flowcell', 'yield_corr_hcad', 'ranges']}, indent=0, default=float))
for r in R['table3']: print(r)
print(pd.DataFrame({k:v for k,v in age.items() if isinstance(v,dict)}).T.round(4)); print('hc only', age['rho_hc_only_Alu'])
print(D[['subject','dgn','Alu','LINE1','SVA','All']].to_string())
