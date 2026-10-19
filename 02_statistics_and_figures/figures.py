import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, ttest_rel, mannwhitneyu
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
RED, BLUE = '#B0281F', '#0D5DA6'
L = pd.read_csv('pass_counts_long.csv')
LAB = {'Alu': 'Alu', 'LINE1': 'LINE1', 'SVA': 'SVA'}
rng = np.random.default_rng(1)

def stars(p):
    return '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

def paired_p(x, fam):
    w = x.pivot_table(index='subject', columns='tissue', values=fam)
    return wilcoxon(w.brain, w.blood).pvalue

def box(ax, groups, colors, labels, title=None, p=None):
    bp = ax.boxplot(groups, widths=0.6, patch_artist=True, showfliers=False)
    for i, c in enumerate(colors):
        for k in ['boxes', 'medians']:
            bp[k][i].set(color=c, linewidth=4 if k == 'boxes' else 6, facecolor='white') if k == 'boxes' else bp[k][i].set(color=c, linewidth=6)
        for k in ['whiskers', 'caps']:
            for j in (2 * i, 2 * i + 1):
                bp[k][j].set(color=c, linewidth=4 if k == 'whiskers' else 0)
        ax.scatter(i + 1 + rng.uniform(-0.2, 0.2, len(groups[i])), groups[i], color=c, s=55, zorder=3)
    ax.set_xticks([1, 2]); ax.set_xticklabels(labels)
    ax.spines[['top', 'right']].set_visible(False)
    if title: ax.set_title(title, loc='left', fontsize=12)
    if p is not None:
        top = max(max(g) for g in groups)
        ax.text(1.5, top * 1.005, stars(p), ha='center', va='bottom', fontsize=16, fontweight='bold')

def fig1():
    for key, sub, lab in [('A', None, None), ('B', 'blood', 'blood'), ('C', 'brain', 'brain')]:
        fig, axes = plt.subplots(1, 3, figsize=(15, 6.5))
        for ax, f in zip(axes, ['Alu', 'LINE1', 'SVA']):
            if key == 'A':
                g = [L[L.tissue == 'blood'][f].values, L[L.tissue == 'brain'][f].values]
                box(ax, g, [RED, BLUE], ['blood', 'brain'], LAB[f], paired_p(L, f))
            else:
                x = L[L.tissue == sub]
                g = [x[x.dgn2 == 'nd'][f].values, x[x.dgn2 == 'hc'][f].values]
                p = mannwhitneyu(g[1], g[0]).pvalue
                box(ax, g, [RED, BLUE], ['nd', 'hc'], LAB[f], p if p < 0.05 else None)
        fig.text(0.01, 0.97, key, fontsize=14, fontweight='bold', va='top')
        fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(f'figs/Figure1{key}.png', dpi=200); plt.close(fig)

def fig2():
    for key, f in zip('ABC', ['Alu', 'LINE1', 'SVA']):
        fig, axes = plt.subplots(1, 3, figsize=(15, 6.5))
        for ax, g in zip(axes, ['hc', 'ad', 'mnd']):
            x = L[L.dgn == g]
            p = paired_p(x, f)
            box(ax, [x[x.tissue == 'blood'][f].values, x[x.tissue == 'brain'][f].values], [RED, BLUE], ['blood', 'brain'], None, p if p < 0.05 else None)
            ax.set_xlabel(g.upper())
        fig.text(0.01, 0.97, LAB[f], fontsize=14, fontweight='bold', va='top')
        fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(f'figs/Figure2{key}.png', dpi=200); plt.close(fig)

def fig3():
    order = [924, 926, 1106, 1107, 1122, 1404, 1446, 1501, 1528, 78, 180, 236, 561, 4011, 7202, 18614, 21699, 31995, 1007, 1511, 1603, 1609, 1613]
    for key, f in zip('ABC', ['Alu', 'LINE1', 'SVA']):
        fig, axes = plt.subplots(4, 6, figsize=(15, 8.5), sharex=True, sharey=True)
        lo, hi = L[f].min(), L[f].max(); pad = (hi - lo) * 0.08
        for ax in axes.flat: ax.set_visible(False)
        for ax, s in zip(axes.flat, order):
            ax.set_visible(True)
            x = L[L.subject == s].set_index('tissue')[f]; g = L[L.subject == s].dgn.iloc[0]
            ax.plot([0, 1], [x['blood'], x['brain']], color=BLUE, lw=2.5, marker='o', ms=9)
            ax.set_title(f'{s}{g}', fontsize=10, backgroundcolor='#EFEFEF')
            ax.set_xlim(-0.4, 1.4); ax.set_ylim(lo - pad, hi + pad); ax.grid(color='#DDDDDD', lw=0.8)
            ax.set_xticks([0, 1]); ax.set_xticklabels(['blood', 'brain'], fontsize=9); ax.tick_params(labelsize=9)
        fig.supylabel(f'{LAB[f]} count'); fig.supxlabel('tissues')
        fig.tight_layout(); fig.savefig(f'figs/Figure3{key}.png', dpi=200); plt.close(fig)

fig1(); fig2(); fig3()
# p-values used for the stars, for the legends
out = {f'1A_{f}': paired_p(L, f) for f in ['Alu', 'LINE1', 'SVA']}
for g in ['hc', 'ad', 'mnd']:
    for f in ['Alu', 'LINE1', 'SVA']: out[f'2_{g}_{f}'] = paired_p(L[L.dgn == g], f)
for t in ['blood', 'brain']:
    x = L[L.tissue == t]
    for f in ['Alu', 'LINE1', 'SVA']: out[f'1{t}_{f}'] = mannwhitneyu(x[x.dgn2 == 'hc'][f], x[x.dgn2 == 'nd'][f]).pvalue
for k, v in out.items(): print(k, f'{v:.3g}')
