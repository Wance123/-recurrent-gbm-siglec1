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

import numpy as np
import pandas as pd
import scanpy as sc


MARKERS = {
    "Myeloid": ["PTPRC", "LST1", "TYROBP", "AIF1", "C1QA", "FCER1G"],
    "T_NK": ["CD3D", "CD3E", "TRBC1", "TRBC2", "CCL5", "NKG7"],
    "B_cell": ["CD79A", "MS4A1", "CD37", "CD74", "HLA-DRA"],
    "Endothelial": ["CLDN5", "EMCN", "PECAM1", "VWF", "KDR"],
    "Pericyte": ["RGS5", "CSPG4", "PDGFRB", "MCAM", "COL4A1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MAG", "MOBP"],
    "OPC": ["PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10"],
    "Astrocyte": ["GFAP", "AQP4", "ALDOC", "SLC1A3", "SLC1A2"],
    "Neuron": ["RBFOX3", "SNAP25", "SYT1", "GRIN1", "SLC17A7"],
    "Tumor_program": ["SOX2", "EGFR", "OLIG2", "CHI3L1", "CD44", "VIM", "NES"],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--cluster-qc-out", type=Path, required=True)
    parser.add_argument("--marker-score-out", type=Path, required=True)
    parser.add_argument("--top-markers-out", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    args = parser.parse_args()

    atlas = sc.read_h5ad(args.input)
    cluster_key = "leiden_0_6"

    counts = (
        atlas.obs.groupby([cluster_key, "sample", "condition"], observed=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    cluster_totals = counts.groupby(cluster_key)["cells"].sum()
    sample_totals = counts.groupby([cluster_key, "sample"])["cells"].sum()
    cluster_qc = []
    for cluster, total in cluster_totals.items():
        by_sample = sample_totals.loc[cluster]
        recurrent = counts.loc[
            (counts[cluster_key] == cluster) & (counts["condition"] == "Recurrent"), "cells"
        ].sum()
        cluster_qc.append(
            {
                "cluster": cluster,
                "cells": int(total),
                "max_patient_fraction": float(by_sample.max() / total),
                "dominant_patient": str(by_sample.idxmax()),
                "patients_with_10_or_more_cells": int((by_sample >= 10).sum()),
                "recurrent_cell_fraction": float(recurrent / total),
                "single_patient_dominant": bool(by_sample.max() / total > 0.70),
            }
        )
    cluster_qc_df = pd.DataFrame(cluster_qc).sort_values("cluster")

    marker_rows = []
    for lineage, genes in MARKERS.items():
        present = [gene for gene in genes if gene in atlas.var_names]
        if not present:
            continue
        values = np.asarray(atlas[:, present].X.mean(axis=1)).ravel()
        atlas.obs[f"score_{lineage}"] = values
        grouped = atlas.obs.groupby(cluster_key, observed=True)[f"score_{lineage}"].mean()
        for cluster, value in grouped.items():
            marker_rows.append(
                {"cluster": cluster, "lineage": lineage, "mean_log_expression": float(value), "genes_present": ";".join(present)}
            )
    marker_df = pd.DataFrame(marker_rows)

    sc.tl.rank_genes_groups(atlas, groupby=cluster_key, method="t-test_overestim_var", n_genes=30)
    marker_tables = []
    for cluster in atlas.obs[cluster_key].cat.categories:
        table = sc.get.rank_genes_groups_df(atlas, group=cluster).head(30)
        table.insert(0, "cluster", cluster)
        marker_tables.append(table)
    top_markers = pd.concat(marker_tables, ignore_index=True)

    for path in [args.cluster_qc_out, args.marker_score_out, args.top_markers_out, args.report_out]:
        path.parent.mkdir(parents=True, exist_ok=True)
    cluster_qc_df.to_csv(args.cluster_qc_out, index=False, encoding="utf-8-sig")
    marker_df.to_csv(args.marker_score_out, index=False, encoding="utf-8-sig")
    top_markers.to_csv(args.top_markers_out, index=False, encoding="utf-8-sig")

    report = {
        "clusters": int(len(cluster_qc_df)),
        "single_patient_dominant_clusters": cluster_qc_df.loc[
            cluster_qc_df["single_patient_dominant"], "cluster"
        ].astype(str).tolist(),
        "clusters_present_in_fewer_than_3_patients": cluster_qc_df.loc[
            cluster_qc_df["patients_with_10_or_more_cells"] < 3, "cluster"
        ].astype(str).tolist(),
    }
    args.report_out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

