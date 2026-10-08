
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path

import numcodecs.blosc as blosc
if not hasattr(blosc, "cbuffer_sizes"):
    blosc.cbuffer_sizes = blosc._cbuffer_sizes
if not hasattr(blosc, "cbuffer_metainfo"):
    blosc.cbuffer_metainfo = blosc._cbuffer_metainfo

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

ROOT = Path(str(project_path('ROOT')))
OUT = ROOT / "further" / "input"
OUT.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(ROOT / "data" / "processed" / "GSE200984_myeloid_reclustered.h5ad")
keep = adata.obs["myeloid_leiden"].astype(str) != "3"
adata = adata[keep].copy()

# Cap abundant patients so that network construction is not dominated by TB4027.
rng = np.random.default_rng(34)
selected = []
for (_, sample), indices in adata.obs.groupby(["condition", "sample"], observed=True).indices.items():
    indices = np.asarray(indices)
    selected.extend(rng.choice(indices, min(150, len(indices)), replace=False).tolist())
adata = adata[np.sort(selected)].copy()

for condition in ["Primary", "Recurrent"]:
    subset = adata[adata.obs["condition"].astype(str) == condition]
    matrix = sparse.csr_matrix(subset.layers["counts"]).T.tocoo()
    io.mmwrite(OUT / f"{condition.lower()}_genes_by_cells.mtx", matrix)
    pd.Series(subset.var_names).to_csv(
        OUT / f"{condition.lower()}_genes.tsv", index=False, header=False
    )
    subset.obs[["sample", "condition", "myeloid_leiden"]].to_csv(
        OUT / f"{condition.lower()}_cells.tsv", sep="\t", encoding="utf-8-sig"
    )

adata.obs.groupby(["condition", "sample"], observed=True).size().rename("cells").reset_index().to_csv(
    OUT / "balanced_cell_counts.csv", index=False, encoding="utf-8-sig"
)
