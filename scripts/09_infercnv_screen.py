from __future__ import annotations

# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path

import argparse
import gzip
import json
import re
from pathlib import Path

import numcodecs.blosc as blosc

if not hasattr(blosc, "cbuffer_sizes"):
    blosc.cbuffer_sizes = blosc._cbuffer_sizes
if not hasattr(blosc, "cbuffer_metainfo"):
    blosc.cbuffer_metainfo = blosc._cbuffer_metainfo

import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc


NORMAL_CLUSTERS = {"4", "9", "10", "11", "12"}
ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def read_gene_positions(gtf: Path) -> dict[str, tuple[str, int, int]]:
    positions: dict[str, tuple[str, int, int]] = {}
    with gzip.open(gtf, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            if chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}:
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            gene = attrs.get("gene_name")
            if not gene:
                continue
            start, end = int(fields[3]), int(fields[4])
            if gene not in positions:
                positions[gene] = (chrom, start, end)
            else:
                old_chrom, old_start, old_end = positions[gene]
                if old_chrom == chrom:
                    positions[gene] = (chrom, min(old_start, start), max(old_end, end))
    return positions


def median_abs_deviation(values: np.ndarray) -> float:
    median = np.median(values)
    return float(np.median(np.abs(values - median)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cluster-out", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    args = parser.parse_args()

    atlas = sc.read_h5ad(args.input)
    positions = read_gene_positions(args.gtf)
    present = atlas.var_names.intersection(pd.Index(positions))
    atlas = atlas[:, present].copy()
    atlas.var["chromosome"] = [positions[gene][0] for gene in atlas.var_names]
    atlas.var["start"] = [positions[gene][1] for gene in atlas.var_names]
    atlas.var["end"] = [positions[gene][2] for gene in atlas.var_names]
    atlas.obs["cnv_reference"] = np.where(
        atlas.obs["leiden_0_6"].astype(str).isin(NORMAL_CLUSTERS), "normal_reference", "test"
    )

    cnv.tl.infercnv(
        atlas,
        reference_key="cnv_reference",
        reference_cat="normal_reference",
        window_size=100,
        step=10,
        dynamic_threshold=1.5,
        chunksize=2000,
        n_jobs=4,
        key_added="cnv",
    )
    cnv_matrix = atlas.obsm["X_cnv"]
    atlas.obs["cnv_rms"] = np.sqrt(np.asarray(cnv_matrix.power(2).mean(axis=1)).ravel())
    atlas.obs["cnv_mean_abs"] = np.asarray(abs(cnv_matrix).mean(axis=1)).ravel()

    reference_scores = atlas.obs.loc[atlas.obs["cnv_reference"] == "normal_reference", "cnv_rms"].to_numpy()
    ref_median = float(np.median(reference_scores))
    ref_mad = median_abs_deviation(reference_scores)
    threshold = ref_median + 3 * 1.4826 * ref_mad
    atlas.obs["cnv_high"] = atlas.obs["cnv_rms"] > threshold

    cluster = (
        atlas.obs.groupby("leiden_0_6", observed=True)
        .agg(
            cells=("cnv_rms", "size"),
            median_cnv_rms=("cnv_rms", "median"),
            mean_cnv_rms=("cnv_rms", "mean"),
            cnv_high_fraction=("cnv_high", "mean"),
            reference_fraction=("cnv_reference", lambda x: float((x == "normal_reference").mean())),
        )
        .reset_index()
    )
    cluster["cnv_candidate_malignant"] = (
        (cluster["median_cnv_rms"] > threshold) | (cluster["cnv_high_fraction"] >= 0.50)
    ) & (cluster["reference_fraction"] == 0)

    for path in [args.output, args.cluster_out, args.report_out]:
        path.parent.mkdir(parents=True, exist_ok=True)
    cluster.to_csv(args.cluster_out, index=False, encoding="utf-8-sig")
    atlas.write_h5ad(args.output, compression="gzip")
    report = {
        "cells": int(atlas.n_obs),
        "genes_with_positions": int(atlas.n_vars),
        "reference_cells": int((atlas.obs["cnv_reference"] == "normal_reference").sum()),
        "reference_median_cnv_rms": ref_median,
        "reference_mad_cnv_rms": ref_mad,
        "cnv_high_threshold": float(threshold),
        "cnv_high_cells": int(atlas.obs["cnv_high"].sum()),
        "candidate_malignant_clusters": cluster.loc[
            cluster["cnv_candidate_malignant"], "leiden_0_6"
        ].astype(str).tolist(),
    }
    args.report_out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

