
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import json
import numcodecs.blosc as blosc
if not hasattr(blosc, "cbuffer_sizes"): blosc.cbuffer_sizes = blosc._cbuffer_sizes
if not hasattr(blosc, "cbuffer_metainfo"): blosc.cbuffer_metainfo = blosc._cbuffer_metainfo
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path(str(project_path('ROOT')))
OUT = ROOT / "further" / "results" / "reinforcement"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(3405)

def raw_matrix(a):
    return a.raw.X.tocsr() if sparse.issparse(a.raw.X) else sparse.csr_matrix(a.raw.X)

def raw_genes(a):
    return pd.Index(a.raw.var_names)

def mean_genes(a, genes):
    idx = raw_genes(a).get_indexer([g for g in genes if g in raw_genes(a)])
    idx = idx[idx >= 0]
    if not len(idx): return np.full(a.n_obs, np.nan)
    return np.asarray(raw_matrix(a)[:, idx].mean(axis=1)).ravel()

def gene_vector(a, gene):
    j = raw_genes(a).get_loc(gene)
    return np.asarray(raw_matrix(a)[:, j].toarray()).ravel()

def cliff_delta(x, y):
    return float(np.mean(np.sign(np.asarray(x)[:, None] - np.asarray(y)[None, :])))

def mw_summary(df, metric):
    p = df.loc[df.condition == "Primary", metric].dropna().to_numpy()
    r = df.loc[df.condition == "Recurrent", metric].dropna().to_numpy()
    obs = float(np.median(r) - np.median(p))
    boot = []
    for _ in range(10000):
        boot.append(np.median(RNG.choice(r, len(r), replace=True)) - np.median(RNG.choice(p, len(p), replace=True)))
    return dict(metric=metric, n_primary=len(p), n_recurrent=len(r),
        primary_median=float(np.median(p)), recurrent_median=float(np.median(r)),
        median_difference_R_minus_P=obs, ci_low=float(np.quantile(boot, .025)),
        ci_high=float(np.quantile(boot, .975)), cliff_delta_R_vs_P=cliff_delta(r, p),
        mannwhitney_p=float(stats.mannwhitneyu(r, p, alternative="two-sided").pvalue))

# Clean myeloid cells; cluster 3 was previously identified as tumor-like contamination.
m = ad.read_h5ad(ROOT / "data" / "processed" / "GSE200984_myeloid_reclustered.h5ad")
m = m[m.obs["myeloid_leiden"].astype(str) != "3"].copy()
sig = gene_vector(m, "SIGLEC1")
m.obs["SIGLEC1_expr"] = sig
m.obs["SIGLEC1_detected"] = sig > 0

patient = m.obs[["sample", "condition"]].copy()
patient["SIGLEC1"] = sig
patient["CD68"] = gene_vector(m, "CD68")
patient = patient.groupby(["sample", "condition"], observed=True).agg(
    SIGLEC1=("SIGLEC1", "mean"), CD68=("CD68", "mean"),
    SIGLEC1_detection=("SIGLEC1", lambda x: float(np.mean(x > 0))), cells=("SIGLEC1", "size")
).reset_index()
patient["SIGLEC1_minus_CD68"] = patient.SIGLEC1 - patient.CD68
patient.to_csv(OUT / "patient_SIGLEC1_metrics.csv", index=False, encoding="utf-8-sig")

metrics = ["SIGLEC1", "SIGLEC1_detection", "SIGLEC1_minus_CD68"]
robust = pd.DataFrame([mw_summary(patient, x) for x in metrics])
robust.to_csv(OUT / "patient_effect_bootstrap.csv", index=False, encoding="utf-8-sig")

loo = []
for omitted in [None] + patient["sample"].tolist():
    d = patient if omitted is None else patient[patient["sample"] != omitted]
    for metric in metrics:
        z = mw_summary(d, metric)
        z["omitted_patient"] = "None" if omitted is None else omitted
        loo.append(z)
loo = pd.DataFrame(loo)
loo.to_csv(OUT / "leave_one_patient_out.csv", index=False, encoding="utf-8-sig")
loo_summary = loo[loo.omitted_patient != "None"].groupby("metric").agg(
    min_effect=("median_difference_R_minus_P", "min"), max_effect=("median_difference_R_minus_P", "max"),
    direction_consistency=("median_difference_R_minus_P", lambda x: float(np.mean(x > 0))),
    min_p=("mannwhitney_p", "min"), max_p=("mannwhitney_p", "max")
).reset_index()
loo_summary.to_csv(OUT / "leave_one_patient_out_summary.csv", index=False, encoding="utf-8-sig")

pathways = {
    "Interferon_response": ["ISG15","IFIT1","IFIT2","IFIT3","IFI6","IFI27","IFI44L","MX1","OAS1","STAT1","IRF7"],
    "Phagolysosome": ["LAMP1","CTSD","CTSB","FCGR1A","TYROBP","AIF1","LILRB1","CTSS"],
    "Lipid_C1Q_APOE": ["APOE","LPL","TREM2","C1QA","C1QB","C1QC","GPNMB","FABP5"],
    "Antigen_presentation": ["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","B2M"],
    "Immunoregulatory": ["IL10","TGFB1","CD274","ARG1","MRC1","MSR1","LGALS3"],
    "Immediate_inflammatory": ["CCL3","CCL4","IL1B","TNF","NFKBIA","CXCL8","PTGS2","FOS","JUN"],
}
for name, genes in pathways.items(): m.obs[name] = mean_genes(m, genes)

deltas = []
for sample, idx in m.obs.groupby("sample", observed=True).indices.items():
    o = m.obs.iloc[idx]
    pos, neg = o.SIGLEC1_detected.to_numpy(), ~o.SIGLEC1_detected.to_numpy()
    if pos.sum() < 5 or neg.sum() < 5: continue
    for pathway in pathways:
        deltas.append(dict(sample=sample, condition=o.condition.iloc[0], pathway=pathway,
            n_positive=int(pos.sum()), n_negative=int(neg.sum()),
            positive_mean=float(o.loc[pos, pathway].mean()), negative_mean=float(o.loc[neg, pathway].mean()),
            delta_positive_minus_negative=float(o.loc[pos, pathway].mean() - o.loc[neg, pathway].mean())))
deltas = pd.DataFrame(deltas)
deltas.to_csv(OUT / "SIGLEC1_positive_state_patient_deltas.csv", index=False, encoding="utf-8-sig")
tests = []
for (condition, pathway), d in deltas.groupby(["condition", "pathway"], observed=True):
    v = d.delta_positive_minus_negative.to_numpy()
    tests.append(dict(condition=condition, pathway=pathway, n_patients=len(v), median_delta=float(np.median(v)),
        positive_fraction=float(np.mean(v > 0)), wilcoxon_p=float(stats.wilcoxon(v).pvalue) if len(v) > 1 and np.any(v != 0) else 1.0))
tests = pd.DataFrame(tests)
tests["FDR"] = tests.groupby("condition")["wilcoxon_p"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
tests.to_csv(OUT / "SIGLEC1_positive_state_statistics.csv", index=False, encoding="utf-8-sig")

# Patient-level ecological association in the full atlas.
a = ad.read_h5ad(ROOT / "data" / "processed" / "GSE200984_strict_basic_atlas.h5ad")
cl = a.obs["leiden_0_6"].astype(str)
malignant = cl.isin(["0","1","2","3","5"])
tnk = cl == "12"
mes_genes = ["CD44","CHI3L1","VIM","LGALS1","TNC","SERPINE1","TIMP1","CAV1"]
exhaustion_genes = ["PDCD1","LAG3","TIGIT","HAVCR2","CTLA4","TOX"]
a.obs["MES_score"] = mean_genes(a, mes_genes)
a.obs["T_exhaustion_score"] = mean_genes(a, exhaustion_genes)
eco = []
for sample, idx in a.obs.groupby("sample", observed=True).indices.items():
    o = a.obs.iloc[idx]
    is_m = malignant.iloc[idx].to_numpy(); is_t = tnk.iloc[idx].to_numpy()
    eco.append(dict(sample=sample, condition=o.condition.iloc[0], total_cells=len(idx),
        malignant_cells=int(is_m.sum()), tnk_cells=int(is_t.sum()),
        malignant_MES=float(o.loc[is_m, "MES_score"].mean()) if is_m.sum() >= 20 else np.nan,
        T_exhaustion=float(o.loc[is_t, "T_exhaustion_score"].mean()) if is_t.sum() >= 5 else np.nan,
        malignant_fraction=float(is_m.mean()), tnk_fraction=float(is_t.mean())))
eco = pd.DataFrame(eco).merge(patient, on=["sample","condition"], how="inner")
eco.to_csv(OUT / "patient_ecology_metrics.csv", index=False, encoding="utf-8-sig")

corr_rows = []
for outcome in ["malignant_MES", "T_exhaustion", "malignant_fraction", "tnk_fraction"]:
    d = eco[["SIGLEC1", outcome]].dropna()
    if len(d) < 5: continue
    rho, pval = stats.spearmanr(d.SIGLEC1, d[outcome])
    corr_rows.append(dict(outcome=outcome, n_patients=len(d), spearman_rho=float(rho), p_value=float(pval)))
corr = pd.DataFrame(corr_rows)
corr["FDR"] = multipletests(corr.p_value, method="fdr_bh")[1]
corr.to_csv(OUT / "patient_ecology_correlations.csv", index=False, encoding="utf-8-sig")

sns.set_theme(style="whitegrid", context="talk")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
sns.stripplot(data=patient, x="condition", y="SIGLEC1", hue="condition", legend=False, ax=ax[0], size=8)
ax[0].set_title("Patient-level SIGLEC1")
plot_tests = tests.copy()
sns.barplot(data=plot_tests, x="median_delta", y="pathway", hue="condition", ax=ax[1])
ax[1].axvline(0, color="black", lw=1); ax[1].set_title("SIGLEC1+ state deltas")
sns.scatterplot(data=eco, x="SIGLEC1", y="malignant_MES", hue="condition", s=90, ax=ax[2])
ax[2].set_title("Myeloid state vs malignant MES")
fig.tight_layout(); fig.savefig(OUT / "reinforcement_summary.png", dpi=300, bbox_inches="tight"); plt.close(fig)

summary = {
    "n_clean_myeloid_cells": int(m.n_obs), "n_patients": int(patient.shape[0]),
    "patient_effects": robust.to_dict("records"), "loo": loo_summary.to_dict("records"),
    "state_tests": tests.to_dict("records"), "ecology_correlations": corr.to_dict("records")
}
(OUT / "reinforcement_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(robust.to_string(index=False)); print(loo_summary.to_string(index=False)); print(tests.to_string(index=False)); print(corr.to_string(index=False))
