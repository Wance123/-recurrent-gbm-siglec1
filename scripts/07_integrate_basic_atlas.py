from __future__ import annotations

# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path

import argparse
import json
from pathlib import Path

import numcodecs.blosc as blosc

if not hasattr(blosc, "cbuffer_sizes"):
    blosc.cbuffer_sizes = blosc._cbuffer_sizes
if not hasattr(blosc, "cbuffer_metainfo"):
    blosc.cbuffer_metainfo = blosc._cbuffer_metainfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--cluster-table-out", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    args = parser.parse_args()

    np.random.seed(34)
    atlas = sc.read_h5ad(args.input)
    atlas = atlas[~atlas.obs["predicted_doublet"].astype(bool)].copy()
    atlas.layers["counts"] = atlas.X.copy()

    sc.pp.normalize_total(atlas, target_sum=1e4)
    sc.pp.log1p(atlas)
    sc.pp.highly_variable_genes(
        atlas, flavor="seurat", n_top_genes=3000, batch_key="sample", inplace=True
    )
    sc.pp.scale(atlas, zero_center=False, max_value=10)
    sc.tl.pca(atlas, n_comps=50, use_highly_variable=True, zero_center=False, random_state=34)
    sc.external.pp.harmony_integrate(atlas, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony")
    sc.pp.neighbors(atlas, n_neighbors=20, n_pcs=40, use_rep="X_pca_harmony", random_state=34)
    sc.tl.umap(atlas, min_dist=0.35, spread=1.0, random_state=34)
    sc.tl.leiden(atlas, resolution=0.6, key_added="leiden_0_6", random_state=34, flavor="igraph", n_iterations=2)

    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.cluster_table_out.parent.mkdir(parents=True, exist_ok=True)
    args.report_out.parent.mkdir(parents=True, exist_ok=True)

    for color, name in [
        ("sample", "umap_by_sample.png"),
        ("condition", "umap_by_condition.png"),
        ("leiden_0_6", "umap_by_cluster.png"),
    ]:
        sc.pl.umap(atlas, color=color, show=False, frameon=False, legend_loc="right margin")
        plt.savefig(args.figure_dir / name, dpi=220, bbox_inches="tight")
        plt.close("all")

    composition = (
        atlas.obs.groupby(["sample", "condition", "leiden_0_6"], observed=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    totals = composition.groupby("sample")["cells"].transform("sum")
    composition["fraction"] = composition["cells"] / totals
    composition.to_csv(args.cluster_table_out, index=False, encoding="utf-8-sig")

    report = {
        "cells_after_doublet_filter": int(atlas.n_obs),
        "genes": int(atlas.n_vars),
        "highly_variable_genes": int(atlas.var["highly_variable"].sum()),
        "clusters_resolution_0_6": int(atlas.obs["leiden_0_6"].nunique()),
        "samples": int(atlas.obs["sample"].nunique()),
        "primary_cells": int((atlas.obs["condition"] == "Primary").sum()),
        "recurrent_cells": int((atlas.obs["condition"] == "Recurrent").sum()),
    }
    args.report_out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    atlas.write_h5ad(args.output, compression="gzip")
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

