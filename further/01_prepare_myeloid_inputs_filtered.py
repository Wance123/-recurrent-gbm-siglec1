
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import numcodecs.blosc as blosc
if not hasattr(blosc, "cbuffer_sizes"): blosc.cbuffer_sizes = blosc._cbuffer_sizes
if not hasattr(blosc, "cbuffer_metainfo"): blosc.cbuffer_metainfo = blosc._cbuffer_metainfo
import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

root = Path(str(project_path('ROOT')))
out = root / "further" / "input"
out.mkdir(parents=True, exist_ok=True)
adata = ad.read_h5ad(root / "data" / "processed" / "GSE200984_myeloid_reclustered.h5ad")
adata = adata[adata.obs["myeloid_leiden"].astype(str) != "3"].copy()

rng = np.random.default_rng(34)
selected = []
for (_, sample), idx in adata.obs.groupby(["condition", "sample"], observed=True).indices.items():
    idx = np.asarray(idx)
    selected.extend(rng.choice(idx, min(150, len(idx)), replace=False).tolist())
adata = adata[np.sort(selected)].copy()

targets = {"SIGLEC1", "STAT1", "IRF7", "IFNAR1", "FCGR1A", "CCL3"}
counts = sparse.csr_matrix(adata.layers["counts"])
detected = np.asarray((counts > 0).sum(axis=0)).ravel()
hvg = adata.var["highly_variable"].to_numpy(dtype=bool)
adata = adata[:, (hvg & (detected >= 25)) | adata.var_names.isin(targets)].copy()

for condition in ("Primary", "Recurrent"):
    subset = adata[adata.obs["condition"].astype(str) == condition]
    io.mmwrite(out / f"{condition.lower()}_genes_by_cells.mtx", sparse.csr_matrix(subset.layers["counts"]).T.tocoo())
    pd.Series(subset.var_names).to_csv(out / f"{condition.lower()}_genes.tsv", index=False, header=False)
    subset.obs[["sample", "condition", "myeloid_leiden"]].to_csv(out / f"{condition.lower()}_cells.tsv", sep="\t", encoding="utf-8-sig")

adata.obs.groupby(["condition", "sample"], observed=True).size().rename("cells").reset_index().to_csv(out / "balanced_cell_counts.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"genes_retained": [adata.n_vars], "minimum_detected_cells": [25], "source_hvg": ["highly_variable"]}).to_csv(out / "gene_filter_summary.csv", index=False, encoding="utf-8-sig")
