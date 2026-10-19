import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
def oneway_tissue(d, fam):
    m = smf.ols(f'Q("{fam}") ~ C(tissue)', d).fit()
    return m.f_pvalue
def posthoc_bonf(d, fam, gcol):
    """DescTools::PostHocTest(aov(y ~ g * tissue), method='bonf'): pooled-MSE t tests between cells,
    Bonferroni over all pairs of cells; returns blood-vs-brain p within each level of g."""
    m = smf.ols(f'Q("{fam}") ~ C({gcol}) * C(tissue)', d).fit()
    mse, df = m.mse_resid, m.df_resid
    g = d.groupby([gcol, 'tissue'])[fam]; means, ns = g.mean(), g.size()
    ncells = len(means); npairs = ncells * (ncells - 1) // 2
    out = {}
    for lev in d[gcol].unique():
        a, b = (lev, 'blood'), (lev, 'brain')
        t = (means[b] - means[a]) / np.sqrt(mse * (1 / ns[a] + 1 / ns[b]))
        out[lev] = min(1.0, 2 * stats.t.sf(abs(t), df) * npairs)
    return out
