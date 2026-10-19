# %% [markdown]
# # Brain vs blood MEI calls — analysis from the original MELT VCFs
#
# This notebook regenerates every number and figure in the slide deck *Brain vs blood MEI calls*, the filter comparison,
# and one CSV per category of insertion, starting only from:
#
# * the three original MELT v2.1.5 VCFs (May 2019): `ALUall.final_comp_original.vcf`, `LINE1.all.final_comp_original.vcf`,
#   `SVAall.final_comp_original.vcf`
# * `RIPs (Brain_Blood) - Patient Info.xlsx` (sheet `SUBJECT INFO`: age at death, sex, brain region, diagnosis;
#   sheet `RIPsPerronProjectwoPD`: sample-to-subject map and the per-sample counts used in the submitted manuscript)
#
# **Definitions**
#
# * A sample *carries* an insertion when its genotype is non-reference (`0/1` or `1/1`).
# * For each person: **brain only** = carried in brain but not in blood; **blood only** = the reverse; **shared** = both.
# * Δ (brain − blood) = number of insertions carried in brain minus the number carried in blood.
# * Filters are applied per site (MELT `FILTER`, `SR`, `ASSESS` are site-level fields).
#
# Sections: 1 setup · 2 load · 3 check against the manuscript · 4 filter sets · 5 slide numbers (chosen filter)
# · 6 Venn diagrams · 7 filter comparison · 8 CSVs per insertion · 9 summary of outputs

# %%
from pathlib import Path
import json, math, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, mannwhitneyu
import statsmodels.formula.api as smf

# ---- paths: edit if your files live elsewhere -------------------------------------------
DATA = Path(".")                       # folder with the three VCFs and the patient-info workbook
VCFS = {"Alu": DATA / "ALUall.final_comp_original.vcf",
        "LINE-1": DATA / "LINE1.all.final_comp_original.vcf",
        "SVA": DATA / "SVAall.final_comp_original.vcf"}
PATIENT_INFO = DATA / "RIPs (Brain_Blood) - Patient Info.xlsx"
OUT = Path("outputs")                  # CSVs and figures are written here
OUT.mkdir(exist_ok=True)

# ---- filter used for the slides (see section 7 for the comparison) ----------------------
CHOSEN = "PASS + SR>2 + ASSESS≥3"

# If a file is not in DATA, look for it in this folder and in ~/Downloads (handles "... (1).xlsx" copies too)
def find(path, pattern):
    if path.exists():
        return path
    hits = sorted((f for d in [Path.cwd(), Path.home() / "Downloads"] for f in d.glob(pattern)),
                  key=lambda f: f.stat().st_mtime, reverse=True)
    if not hits:
        raise FileNotFoundError(f"{path.name} not found; set DATA (or the path) at the top of this cell")
    return hits[0]
VCFS = {"Alu": find(VCFS["Alu"], "ALUall*.vcf"), "LINE-1": find(VCFS["LINE-1"], "LINE1*all*.vcf"),
        "SVA": find(VCFS["SVA"], "SVAall*.vcf")}
PATIENT_INFO = find(PATIENT_INFO, "*Patient Info*.xlsx")

FAMILIES = list(VCFS)
pd.set_option("display.width", 200, "display.max_columns", 40)
for p in list(VCFS.values()) + [PATIENT_INFO]:
    print("using", p)

# %% [markdown]
# ## 2. Load the VCFs and the patient information

# %%
def read_melt_vcf(path, family):
    """One row per site: VCF columns, every INFO field, and the genotype (GT) and likelihoods (GL) per sample."""
    sites, samples = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = [s.split(".")[0] for s in f[9:]]       # '09_24_T.bwa.sorted.marked' -> '09_24_T'
                continue
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            fmt = f[8].split(":")
            rec = {"family": family, "CHROM": f[0], "POS": int(f[1]), "ID": f[2], "REF": f[3], "ALT": f[4],
                   "FILTER": f[6], **info}
            for s, v in zip(samples, f[9:]):
                d = dict(zip(fmt, v.split(":")))
                rec[f"GT|{s}"] = d.get("GT", "./.")
                rec[f"GL|{s}"] = d.get("GL", "")
            sites.append(rec)
    df = pd.DataFrame(sites)
    for c in ["SR", "ASSESS", "LP", "RP", "SVLEN"]:
        df[c] = pd.to_numeric(df[c])
    df["RA"] = pd.to_numeric(df["RA"], errors="coerce")
    df["site_id"] = df["family"] + ":" + df["CHROM"] + ":" + df["POS"].astype(str)
    return df, samples

sites, samples = {}, None
for fam, path in VCFS.items():
    sites[fam], s = read_melt_vcf(path, fam)
    assert samples is None or s == samples
    samples = s
    print(f"{fam:7s} {len(sites[fam]):5d} sites   source: {path.name}")
print(len(samples), "samples in the VCFs (including", sum(s.startswith("PD") for s in samples), "PD samples, not used)")

# %%
info = pd.read_excel(PATIENT_INFO, sheet_name="SUBJECT INFO").dropna(subset=["subject"])
info = info.rename(columns={"age at death": "age", "brain region": "region", "diagnosis 1": "diagnosis"})
info["subject"] = info["subject"].astype(int)
info["group"] = info["diagnosis"].map({"Healthy Control": "HC", "Alzheimer's Disease": "AD",
                                       "Motor Neurone Disease": "MND"})
info = info[["subject", "group", "age", "sex", "region", "diagnosis"]].set_index("subject")
assert info["age"].notna().all() and len(info) == 23

sheet = pd.read_excel(PATIENT_INFO, sheet_name="RIPsPerronProjectwoPD").dropna(subset=["sampleID"])
sheet.columns = [c.strip() for c in sheet.columns]
sheet["subject"] = sheet["subject"].astype(int)
sample_map = sheet.set_index("sampleID")[["subject", "tissue"]]          # e.g. '09_24_T' -> (924, 'brain')
assert set(sample_map.index) <= set(samples), set(sample_map.index) - set(samples)

# brain and blood sample name per subject
pairs = sample_map.reset_index().pivot(index="subject", columns="tissue", values="sampleID")
pairs = pairs.join(info)
GROUP_ORDER = {"HC": 0, "AD": 1, "MND": 2}
pairs = pairs.sort_values(["group", "subject"], key=lambda s: s.map(GROUP_ORDER) if s.name == "group" else s)
print(pairs[["brain", "blood", "group", "age", "sex", "region"]].to_string())

# %% [markdown]
# Genotype matrix: one boolean per site × sample (True = non-reference genotype).

# %%
def carrier_matrix(df):
    gt = df[[f"GT|{s}" for s in samples]]
    gt.columns = samples
    return ~gt.isin(["0/0", "./.", "0|0", "."])

carriers = {fam: carrier_matrix(df) for fam, df in sites.items()}

# %% [markdown]
# ## 3. Check: which filter reproduces the counts in the submitted manuscript?
#
# The per-sample counts used in the manuscript are in the sheet `RIPsPerronProjectwoPD`.

# %%
def counts_for(mask_by_family):
    return pd.DataFrame({fam: carriers[fam][mask_by_family[fam]].sum() for fam in FAMILIES})

manuscript_counts = sheet.set_index("sampleID")[["Alu", "LINE1", "SVA"]].rename(columns={"LINE1": "LINE-1"})
candidates = {
    "no filter":                    {f: pd.Series(True, index=sites[f].index) for f in FAMILIES},
    "SR>2":                         {f: sites[f]["SR"] > 2 for f in FAMILIES},
    "SR>2 (Alu, LINE-1); none (SVA)": {f: (sites[f]["SR"] > 2) if f != "SVA" else pd.Series(True, index=sites[f].index)
                                       for f in FAMILIES},
}
rows = []
for name, masks in candidates.items():
    c = counts_for(masks).loc[manuscript_counts.index]
    rows.append({"filter": name, **{f"{f}: samples matching": int((c[f] == manuscript_counts[f]).sum()) for f in FAMILIES}})
check = pd.DataFrame(rows).set_index("filter")
display(check)
assert (check.loc["SR>2 (Alu, LINE-1); none (SVA)"] == 46).all()
print("The manuscript counts = SR > 2 for Alu and LINE-1, and no filter for SVA (46/46 samples exact).")

# %% [markdown]
# ## 4. Filter sets compared
#
# MELT `FILTER` flags: `lc` low-complexity region · `rSD` LP/RP ratio > 2 SD · `hDP` excess split discordant pairs ·
# `s25` > 25% of samples without data · `ac0` no carriers. `PASS` = none of these.
# `ASSESS`: 3 = left-side TSD only, 4 = right-side TSD only, 5 = TSD decided with split reads.

# %%
def site_mask(fam, name):
    d = sites[fam]
    allsites = pd.Series(True, index=d.index)
    return {
        "Manuscript (SR>2; SVA unfiltered)": (d.SR > 2) if fam != "SVA" else allsites,
        "SR>2":                     d.SR > 2,
        "SR>2 + ASSESS≥3":          (d.SR > 2) & (d.ASSESS >= 3),
        "PASS":                     d.FILTER == "PASS",
        "PASS + SR>2 + ASSESS≥3":   (d.FILTER == "PASS") & (d.SR > 2) & (d.ASSESS >= 3),
        "no lc + SR>2 + ASSESS≥3":  (~d.FILTER.str.contains("lc")) & (d.SR > 2) & (d.ASSESS >= 3),
    }[name]

FILTERS = ["Manuscript (SR>2; SVA unfiltered)", "SR>2", "SR>2 + ASSESS≥3", "PASS",
           "PASS + SR>2 + ASSESS≥3", "no lc + SR>2 + ASSESS≥3"]
assert CHOSEN in FILTERS
pd.DataFrame({name: {fam: int(site_mask(fam, name).sum()) for fam in FAMILIES} for name in FILTERS}).T \
  .assign(total=lambda x: x.sum(axis=1)).rename_axis("sites kept")

# %% [markdown]
# Core per-person tables for any filter: counts per tissue, Δ, and brain-only / shared / blood-only.

# %%
def person_tables(filter_name):
    rows = []
    for fam in FAMILIES:
        m = carriers[fam][site_mask(fam, filter_name)]
        for subj, r in pairs.iterrows():
            b, l = m[r["brain"]], m[r["blood"]]
            rows.append({"subject": subj, "family": fam, "brain": int(b.sum()), "blood": int(l.sum()),
                         "brain_only": int((b & ~l).sum()), "blood_only": int((l & ~b).sum()),
                         "shared": int((b & l).sum())})
    t = pd.DataFrame(rows)
    allte = t.groupby("subject")[["brain", "blood", "brain_only", "blood_only", "shared"]].sum().reset_index()
    allte["family"] = "All TEs"
    t = pd.concat([allte, t], ignore_index=True)
    t["delta"] = t["brain"] - t["blood"]
    tot = t[["brain_only", "blood_only", "shared"]].sum(axis=1)
    for c in ["brain_only", "blood_only", "shared"]:
        t[c + "_pct"] = 100 * t[c] / tot
    return t.merge(pairs[["group", "age", "sex", "region"]], left_on="subject", right_index=True)

GROUPS = {"HC": ["HC"], "AD + MND": ["AD", "MND"], "AD": ["AD"], "MND": ["MND"]}

def paired_table(t):
    rows = []
    for g, members in GROUPS.items():
        for fam in ["All TEs"] + FAMILIES:
            x = t[(t.family == fam) & t.group.isin(members)]
            b, r = x.blood.mean(), x.brain.mean()
            rows.append({"group": g, "family": fam, "n": len(x), "blood mean": round(b), "brain mean": round(r),
                         "Δ (brain − blood)": round(r) - round(b), "paired Wilcoxon p": wilcoxon(x.brain, x.blood).pvalue})
    return pd.DataFrame(rows)

def venn_stats(t):
    a = t[t.family == "All TEs"].copy()
    a["HC"] = (a.group == "HC").astype(int)
    hc, dz = a[a.HC == 1], a[a.HC == 0]
    model = smf.ols("brain_only_pct ~ HC + age", data=a).fit()
    return {
        "HC brain only %": hc.brain_only_pct.mean(), "HC shared %": hc.shared_pct.mean(), "HC blood only %": hc.blood_only_pct.mean(),
        "AD+MND brain only %": dz.brain_only_pct.mean(), "AD+MND shared %": dz.shared_pct.mean(),
        "AD+MND blood only %": dz.blood_only_pct.mean(),
        "brain only %, HC vs AD+MND: Mann–Whitney p": mannwhitneyu(hc.brain_only_pct, dz.brain_only_pct).pvalue,
        "brain only %: age-adjusted HC effect (points)": model.params["HC"],
        "brain only %: age-adjusted p": model.pvalues["HC"],
        # "excess" = brain-only % minus blood-only %; ~0 when brain and blood differ only by genotyping noise
        "people with excess brain-only calls (> 5 points)": int((a.brain_only_pct - a.blood_only_pct > 5).sum()),
        "who (excess > 5 points)": ", ".join(str(s) for s in a[a.brain_only_pct - a.blood_only_pct > 5].subject),
    }

# %% [markdown]
# ## 5. Slide 1 — brain vs blood counts per group (chosen filter)

# %%
T = person_tables(CHOSEN)
slide1 = paired_table(T)
print("Filter:", CHOSEN)
display(slide1.style.format({"paired Wilcoxon p": "{:.2g}"}))

ages = pairs.assign(dz=pairs.group != "HC")
age_summary = pd.concat([ages[ages.group == g].age.describe()[["count", "mean", "std", "min", "50%", "max"]].rename(g)
                         for g in ["HC", "AD", "MND"]] +
                        [ages[ages.dz].age.describe()[["count", "mean", "std", "min", "50%", "max"]].rename("AD + MND")], axis=1).T
age_summary = age_summary.rename(columns={"50%": "median"}).round(1)
age_summary

# %% [markdown]
# ## 6. Slides 2 and 3 — brain only / shared / blood only (chosen filter)

# %%
V = venn_stats(T)
display(pd.Series(V, name=CHOSEN).to_frame())
per_person = T[T.family == "All TEs"][["subject", "group", "age", "brain_only_pct", "shared_pct", "blood_only_pct",
                                       "brain_only", "shared", "blood_only"]].round(1)
per_person.sort_values(["group", "brain_only_pct"], ascending=[True, False],
                       key=lambda s: s.map(GROUP_ORDER) if s.name == "group" else s)

# %% [markdown]
# Area-proportional Venn diagrams: circle area = % of the person's calls found in that tissue, overlap area = % shared,
# same scale in every panel.

# %%
def lens(r1, r2, d):
    if d >= r1 + r2: return 0.0
    if d <= abs(r1 - r2): return math.pi * min(r1, r2) ** 2
    a1 = r1**2 * math.acos((d*d + r1*r1 - r2*r2) / (2*d*r1))
    a2 = r2**2 * math.acos((d*d + r2*r2 - r1*r1) / (2*d*r2))
    return a1 + a2 - 0.5 * math.sqrt((-d+r1+r2) * (d+r1-r2) * (d-r1+r2) * (d+r1+r2))

def venn(ax, bo, sh, bl, cx=0.0, cy=0.0, k=1.0, fs=11, decimals=0):
    r1, r2 = k * math.sqrt((bo + sh) / 100), k * math.sqrt((bl + sh) / 100)
    target, lo, hi = math.pi * k * k * sh / 100, abs(r1 - r2), r1 + r2
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if lens(r1, r2, mid) > target else (lo, mid)
    d = (lo + hi) / 2
    off = cx - (d + r2 - r1) / 2
    x1, x2 = off, off + d
    dark = bo > 10
    ax.add_patch(plt.Circle((x1, cy), r1, fc="#5B4BB7", alpha=0.55 if dark else 0.2, ec="#3C3489", lw=1.5))
    ax.add_patch(plt.Circle((x2, cy), r2, fc="#1D9E75", alpha=0.22, ec="#0F6E56", lw=1.5))
    fmt = (lambda v: f"{v:.{decimals}f}%") if decimals else (lambda v: f"{v:.0f}")
    g = 0.08 * k
    ax.text(x1 - r1 - g, cy, fmt(bo), ha="right", va="center", fontsize=fs, fontweight="bold", color="#3C3489")
    ax.text(((x2 - r2) + (x1 + r1)) / 2, cy, fmt(sh), ha="center", va="center", fontsize=fs, color="#1E2430")
    ax.text(x2 + r2 + g, cy, fmt(bl), ha="left", va="center", fontsize=fs, fontweight="bold", color="#0F6E56")
    if decimals:   # group figure: name the two outer numbers
        ax.text(x1 - r1 - g, cy - 0.16 * k, "brain only", ha="right", va="center", fontsize=fs * 0.7, color="#3C3489")
        ax.text(x2 + r2 + g, cy - 0.16 * k, "blood only", ha="left", va="center", fontsize=fs * 0.7, color="#0F6E56")

def group_venn_figure(t, title_suffix=""):
    a = t[t.family == "All TEs"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    for ax, (label, members) in zip(axes, [("Healthy controls (n = 9)", ["HC"]), ("AD + MND (n = 14)", ["AD", "MND"])]):
        x = a[a.group.isin(members)][["brain_only_pct", "shared_pct", "blood_only_pct"]].mean()
        venn(ax, *x, k=1.0, fs=15, decimals=1)
        ax.set_title(label, fontsize=16, fontweight="bold")
        ax.set_xlim(-1.9, 1.9); ax.set_ylim(-1.2, 1.2); ax.set_aspect("equal"); ax.axis("off")
    fig.suptitle(title_suffix, fontsize=10, color="#4A5162", y=0.02)
    fig.tight_layout()
    return fig

def person_venn_figure(t):
    a = t[t.family == "All TEs"].sort_values("brain_only_pct", ascending=False)
    hc, dz = a[a.group == "HC"], a[a.group != "HC"]
    fig = plt.figure(figsize=(16.6, 7))
    for panel, (df, cols, x0, w, title, tag) in enumerate([(hc, 3, 0.0, 0.40, "Healthy controls (n = 9)", False),
                                                          (dz, 4, 0.42, 0.58, "AD + MND (n = 14)", True)]):
        rows = math.ceil(len(df) / cols)
        fig.text(x0 + w / 2, 0.97, title, ha="center", fontsize=14, fontweight="bold")
        for i, r in enumerate(df.itertuples()):
            ax = fig.add_axes([x0 + w * (i % cols) / cols, 0.92 - (i // cols + 1) * 0.92 / max(rows, 4), w / cols, 0.92 / max(rows, 4)])
            venn(ax, r.brain_only_pct, r.shared_pct, r.blood_only_pct, k=1.0, fs=10)
            ax.set_title(f"{r.subject}" + (f" · {r.group}" if tag else ""), fontsize=11, fontweight="bold", pad=2)
            ax.set_xlim(-1.9, 1.9); ax.set_ylim(-1.1, 1.1); ax.set_aspect("equal"); ax.axis("off")
    return fig

fig = group_venn_figure(T, f"Mean per-person % of all TE calls · filter: {CHOSEN}")
fig.savefig(OUT / "venn_groups.png", dpi=200, bbox_inches="tight"); plt.show()
fig = person_venn_figure(T)
fig.savefig(OUT / "venn_per_person.png", dpi=200, bbox_inches="tight"); plt.show()

# %% [markdown]
# ## 7. Filter comparison
#
# Same analyses under each filter set. What to look for: does the brain > blood excess in controls, the controls vs
# patients difference, and the same group of high brain-only people hold whatever filter is used?
#
# Note: in these VCFs every site with SR > 2 already has ASSESS ≥ 3, so "SR>2" and "SR>2 + ASSESS≥3" are identical.
# "Excess brain-only" = a person's brain-only % minus their blood-only %; blood-only calls measure the background
# mismatch between two samples of the same genome, so the excess is comparable across filters.

# %%
comp_rows, per_person_by_filter = [], {}
for name in FILTERS:
    t = person_tables(name)
    pt = paired_table(t).set_index(["group", "family"])
    v = venn_stats(t)
    a = t[t.family == "All TEs"]
    per_person_by_filter[name] = a.set_index("subject").brain_only_pct
    comp_rows.append({
        "filter": name,
        "sites kept": sum(int(site_mask(f, name).sum()) for f in FAMILIES),
        "calls per person (blood, mean)": round(a.blood.mean()),
        "HC Δ all TEs": pt.loc[("HC", "All TEs"), "Δ (brain − blood)"],
        "HC paired p": pt.loc[("HC", "All TEs"), "paired Wilcoxon p"],
        "AD+MND Δ all TEs": pt.loc[("AD + MND", "All TEs"), "Δ (brain − blood)"],
        "AD+MND paired p": pt.loc[("AD + MND", "All TEs"), "paired Wilcoxon p"],
        "HC brain only %": v["HC brain only %"], "AD+MND brain only %": v["AD+MND brain only %"],
        "HC blood only %": v["HC blood only %"], "AD+MND blood only %": v["AD+MND blood only %"],
        "brain only %, MWU p": v["brain only %, HC vs AD+MND: Mann–Whitney p"],
        "brain only %, age-adjusted p": v["brain only %: age-adjusted p"],
        "people with excess brain-only (> 5 points)": v["people with excess brain-only calls (> 5 points)"],
        "who": v["who (excess > 5 points)"],
    })
comparison = pd.DataFrame(comp_rows).set_index("filter")
comparison.to_csv(OUT / "filter_comparison.csv")
display(comparison.drop(columns="who").style.format(precision=3).format(
    {c: "{:.2g}" for c in ["HC paired p", "AD+MND paired p", "brain only %, MWU p", "brain only %, age-adjusted p"]}).format(
    {c: "{:.1f}" for c in comparison.columns if c.endswith("%")}))
comparison[["people with excess brain-only (> 5 points)", "who"]]

# %%
# Brain-only % per person under each filter (rows = people, columns = filters)
bo = pd.DataFrame(per_person_by_filter).join(pairs[["group"]]).round(1)
bo = bo.sort_values(["group", CHOSEN], ascending=[True, False],
                    key=lambda s: s.map(GROUP_ORDER) if s.name == "group" else s)
bo.to_csv(OUT / "brain_only_pct_by_filter.csv")
bo

# %%
fig, ax = plt.subplots(figsize=(11, 4.5))
for i, name in enumerate(FILTERS):
    x = per_person_by_filter[name]
    for subj, val in x.items():
        g = pairs.loc[subj, "group"]
        ax.scatter(i + {"HC": -0.15, "AD": 0.05, "MND": 0.15}[g] + np.random.default_rng(subj).uniform(-0.05, 0.05), val,
                   color={"HC": "#5B4BB7", "AD": "#D9822B", "MND": "#1D9E75"}[g], s=28, alpha=0.85,
                   label=g if (i == 0 and subj == pairs[pairs.group == g].index[0]) else None)
ax.set_xticks(range(len(FILTERS))); ax.set_xticklabels([f.replace(" (", "\n(").replace(" + ", "\n+ ", 1) for f in FILTERS], fontsize=9)
ax.set_ylabel("brain-only calls (% of person's calls)")
h, l = ax.get_legend_handles_labels(); order = [l.index(g) for g in ["HC", "AD", "MND"]]
ax.legend([h[i] for i in order], [l[i] for i in order], frameon=False); ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Brain-only % per person under each filter set")
fig.tight_layout(); fig.savefig(OUT / "filter_comparison_brain_only.png", dpi=200); plt.show()

# %% [markdown]
# Which MELT flags do the brain-only calls carry? (all sites with SR > 2 and ASSESS ≥ 3)

# %%
high_people = {int(x) for x in V["who (excess > 5 points)"].split(", ")}
flag_rows = []
for fam in FAMILIES:
    d = sites[fam]
    base = (d.SR > 2) & (d.ASSESS >= 3)
    m = carriers[fam]
    for subj, r in pairs.iterrows():
        b, l = m[r["brain"]], m[r["blood"]]
        for cat, sel in [("brain only", b & ~l), ("blood only", l & ~b), ("shared", b & l)]:
            for flt in d.loc[base & sel, "FILTER"]:
                flag_rows.append({"category": cat, "high brain-only person": subj in high_people,
                                  "PASS": flt == "PASS", "lc": "lc" in flt, "rSD": "rSD" in flt, "hDP": "hDP" in flt})
flags = pd.DataFrame(flag_rows)
(flags.groupby(["category", "high brain-only person"])[["PASS", "lc", "rSD", "hDP"]].mean() * 100).round(1) \
    .rename(columns=lambda c: f"% {c}")

# %% [markdown]
# ## 8. CSVs with every insertion
#
# * `outputs/sites_all.csv` — one row per MELT site (all filters), with every INFO field, a column per filter set
#   (TRUE if the site is kept) and the number of brain / blood samples carrying it.
# * `outputs/calls_brain_only.csv`, `outputs/calls_blood_only.csv`, `outputs/calls_shared.csv` — one row per insertion
#   per person, with subject, group, age, sex, brain region, every MELT field, genotype and likelihoods in both tissues,
#   how many other people carry it, and the filter-set columns (so any filter can be applied in Excel).
# * `outputs/counts_per_person.csv` — per person × family × filter: counts in brain and blood, Δ, brain only / shared /
#   blood only.

# %%
INFO_COLS = ["family", "CHROM", "POS", "ID", "REF", "ALT", "FILTER", "ASSESS", "SR", "LP", "RP", "RA", "TSD",
             "SVTYPE", "SVLEN", "MEINFO", "INTERNAL", "PRIOR", "DIFF"]
brain_samples, blood_samples = list(pairs["brain"]), list(pairs["blood"])
site_tables = []
for fam in FAMILIES:
    d = sites[fam].copy()
    for name in FILTERS:
        d[f"kept: {name}"] = site_mask(fam, name).values
    d["n brain samples carrying"] = carriers[fam][brain_samples].sum(axis=1)
    d["n blood samples carrying"] = carriers[fam][blood_samples].sum(axis=1)
    site_tables.append(d)
all_sites = pd.concat(site_tables, ignore_index=True)
FILTER_COLS = [f"kept: {n}" for n in FILTERS]
all_sites[["site_id"] + [c for c in INFO_COLS if c in all_sites] + FILTER_COLS +
          ["n brain samples carrying", "n blood samples carrying"]].to_csv(OUT / "sites_all.csv", index=False)

call_rows = []
for d in site_tables:
    fam = d["family"].iloc[0]
    m = carriers[fam]
    any_person = pd.DataFrame({s: m[r["brain"]] | m[r["blood"]] for s, r in pairs.iterrows()})
    n_people = any_person.sum(axis=1)
    for subj, r in pairs.iterrows():
        b, l = m[r["brain"]], m[r["blood"]]
        for cat, sel in [("brain only", b & ~l), ("blood only", l & ~b), ("shared", b & l)]:
            x = d.loc[sel, ["site_id"] + [c for c in INFO_COLS if c in d] + FILTER_COLS].copy()
            x.insert(0, "category", cat)
            x.insert(1, "subject", subj)
            for c in ["group", "age", "sex", "region"]:
                x.insert(x.columns.get_loc("subject") + 1 + ["group", "age", "sex", "region"].index(c), c, r[c])
            x["brain_sample"], x["blood_sample"] = r["brain"], r["blood"]
            x["brain_GT"], x["brain_GL"] = d.loc[sel, f"GT|{r['brain']}"], d.loc[sel, f"GL|{r['brain']}"]
            x["blood_GT"], x["blood_GL"] = d.loc[sel, f"GT|{r['blood']}"], d.loc[sel, f"GL|{r['blood']}"]
            x["n other people carrying"] = n_people[sel] - 1
            call_rows.append(x)
calls = pd.concat(call_rows, ignore_index=True)
for cat in ["brain only", "blood only", "shared"]:
    calls[calls.category == cat].drop(columns="category").to_csv(OUT / f"calls_{cat.replace(' ', '_')}.csv", index=False)

counts_long = pd.concat([person_tables(n).assign(filter=n) for n in FILTERS], ignore_index=True)
counts_long.round(2).to_csv(OUT / "counts_per_person.csv", index=False)

summary = calls.groupby(["category"]).agg(rows=("site_id", "size")).join(
    calls[calls[f"kept: {CHOSEN}"]].groupby("category").agg(**{f"rows kept by {CHOSEN}": ("site_id", "size")}))
summary

# %%
# Consistency check: the CSVs give back the counts used for the slides
k = calls[calls[f"kept: {CHOSEN}"]]
from_csv = k.groupby(["subject", "category"]).size().unstack()
direct = T[T.family == "All TEs"].set_index("subject")[["brain_only", "shared", "blood_only"]]
assert (from_csv["brain only"] == direct["brain_only"]).all() and (from_csv["blood only"] == direct["blood_only"]).all() \
    and (from_csv["shared"] == direct["shared"]).all()
print("CSV rows match the per-person counts used in the slides")

# %% [markdown]
# ## 9. Numbers used in the slides
#
# Written to `outputs/slide_numbers.json` so the deck can be checked against this notebook.

# %%
def f(x): return float(x) if isinstance(x, (np.floating, float)) else (int(x) if isinstance(x, (np.integer, int)) else x)
slide_numbers = {
    "filter": CHOSEN,
    "slide1": [{k: f(v) for k, v in r.items()} for r in slide1.to_dict("records")],
    "ages": {g: {k: f(v) for k, v in r.items()} for g, r in age_summary.to_dict("index").items()},
    "venn": {k: f(v) for k, v in V.items()},
    "per_person": [{k: f(v) for k, v in r.items()} for r in per_person.to_dict("records")],
}
(OUT / "slide_numbers.json").write_text(json.dumps(slide_numbers, indent=1, ensure_ascii=False))
print("outputs:", *sorted(p.name for p in OUT.iterdir()), sep="\n  ")

# %% [markdown]
# ## 10. Manuscript statistics (Tables 1–3, age analyses) for the chosen filter
#
# Brain vs blood: paired Wilcoxon signed-rank tests (all 23 pairs for Table 2, within each group for Table 3);
# the "p-adjusted" columns multiply by 3 for the three TE families (Bonferroni). Age analyses as in the response to Reviewer 1.

# %%
from scipy import stats as _st
from scipy.stats import spearmanr
FAMCOL = {"Alu": "Alu", "LINE-1": "LINE1", "SVA": "SVA"}
Lg = T[T.family != "All TEs"].melt(id_vars=["subject", "family", "group", "age"], value_vars=["blood", "brain"],
                                   var_name="tissue", value_name="n")
Lg = Lg.pivot_table(index=["subject", "group", "age", "tissue"], columns="family", values="n").reset_index()
Lg = Lg.rename(columns=FAMCOL); Lg["dgn"] = Lg.group.str.lower(); Lg["dgn2"] = Lg.dgn.map({"hc": "hc", "ad": "nd", "mnd": "nd"})

def posthoc_bonf(d, fam, gcol):
    m = smf.ols(f"{fam} ~ C({gcol}) * C(tissue)", d).fit()
    g = d.groupby([gcol, "tissue"])[fam]; means, ns = g.mean(), g.size()
    npairs = len(means) * (len(means) - 1) // 2
    out = {}
    for lev in d[gcol].unique():
        a, b = (lev, "blood"), (lev, "brain")
        t = (means[b] - means[a]) / np.sqrt(m.mse_resid * (1 / ns[a] + 1 / ns[b]))
        out[lev] = min(1.0, 2 * _st.t.sf(abs(t), m.df_resid) * npairs)
    return out

t2 = []
for f in ["Alu", "LINE1", "SVA"]:
    w = Lg.pivot_table(index="subject", columns="tissue", values=f); p = wilcoxon(w.brain, w.blood).pvalue
    b, r = Lg[Lg.tissue == "blood"][f], Lg[Lg.tissue == "brain"][f]
    t2.append({"TE": f, "Blood, n=23": f"{b.mean():.0f} ({b.std():.0f})", "Brain, n=23": f"{r.mean():.0f} ({r.std():.0f})",
               "p": p, "p-adjusted": min(1, 3 * p)})
table2 = pd.DataFrame(t2); display(table2)
t3 = []
for f in ["Alu", "LINE1", "SVA"]:
    ph, ph2 = posthoc_bonf(Lg, f, "dgn"), posthoc_bonf(Lg, f, "dgn2")
    for g, lab in [("hc", "HC, n=9"), ("mnd", "MND, n=9"), ("ad", "AD, n=5"), ("nd", "ND, n=14")]:
        col = "dgn2" if g == "nd" else "dgn"
        x = Lg[Lg[col] == g]; b, r = x[x.tissue == "blood"][f], x[x.tissue == "brain"][f]
        w = x.pivot_table(index="subject", columns="tissue", values=f); p = wilcoxon(w.brain, w.blood).pvalue
        t3.append({"TE": f, "Dgn": lab, "Blood": f"{b.mean():.0f} ({b.std():.0f})", "Brain": f"{r.mean():.0f} ({r.std():.0f})",
                   "p": p, "p-adjusted": min(1, 3 * p)})
table3 = pd.DataFrame(t3); display(table3)
table1 = T[T.family != "All TEs"].pivot_table(index=["subject", "group", "age", "sex", "region"], columns="family", values="delta").reset_index()
table1.to_csv(OUT / "table1_deltas.csv", index=False); table2.to_csv(OUT / "table2.csv", index=False); table3.to_csv(OUT / "table3.csv", index=False)
table1

# %%
Dd = T.pivot_table(index=["subject", "group", "age"], columns="family", values="delta").reset_index()
Dd = Dd.rename(columns={"All TEs": "All", "LINE-1": "LINE1"}); Dd["HC"] = (Dd.group == "HC").astype(int)
hc_, dz_ = Dd[Dd.HC == 1], Dd[Dd.HC == 0]; hcad_ = Dd[Dd.group.isin(["HC", "AD"])]
age_rows = []
for f in ["All", "Alu", "LINE1", "SVA"]:
    m = smf.ols(f"{f} ~ HC + age", Dd).fit()
    age_rows.append({"family": f, "median Δ controls": hc_[f].median(), "median Δ patients": dz_[f].median(),
                     "Mann–Whitney p": mannwhitneyu(hc_[f], dz_[f]).pvalue, "age-adjusted control effect": m.params.HC,
                     "p (control)": m.pvalues.HC, "p (age)": m.pvalues.age,
                     "Spearman ρ, all 23": spearmanr(Dd.age, Dd[f])[0], "p, all 23": spearmanr(Dd.age, Dd[f])[1],
                     "Spearman ρ, HC + AD": spearmanr(hcad_.age, hcad_[f])[0], "p, HC + AD": spearmanr(hcad_.age, hcad_[f])[1]})
age_table = pd.DataFrame(age_rows).set_index("family"); age_table.to_csv(OUT / "age_analyses.csv")
print("Age, patients vs controls: Mann–Whitney p =", round(mannwhitneyu(pairs[pairs.group != "HC"].age, pairs[pairs.group == "HC"].age).pvalue, 3))
age_table.round(4)
