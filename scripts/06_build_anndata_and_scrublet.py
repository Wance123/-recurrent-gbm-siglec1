from __future__ import annotations

# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path

import argparse
import csv
import gzip
import json
import shlex
from pathlib import Path

import numcodecs.blosc as blosc

# zarr 2.18.x compatibility with numcodecs 0.16.x in the recorded D: Python environment.
if not hasattr(blosc, "cbuffer_sizes"):
    blosc.cbuffer_sizes = blosc._cbuffer_sizes
if not hasattr(blosc, "cbuffer_metainfo"):
    blosc.cbuffer_metainfo = blosc._cbuffer_metainfo

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def split_fields(line: str) -> list[str]:
    line = line.rstrip("\r\n")
    return line.split("\t") if "\t" in line else shlex.split(line)


def feature_name(line: str) -> str:
    tab = line.find("\t")
    space = line.find(" ")
    stops = [value for value in (tab, space) if value >= 0]
    if not stops:
        return line.strip().strip('"')
    return line[: min(stops)].strip().strip('"')


def canonical_gene(gene: str) -> str:
    upper = gene.upper()
    if upper.startswith("MT."):
        return "MT-" + gene[3:]
    return gene


def read_genes(path: Path) -> list[str]:
    genes: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        handle.readline()
        for line in handle:
            genes.append(canonical_gene(feature_name(line)))
    if len(genes) != len(set(genes)):
        raise ValueError(f"Duplicate canonical genes in {path.name}")
    return genes


def build_sample(
    path: Path,
    sample: str,
    condition: str,
    keep_barcodes: set[str],
    common_genes: list[str],
    clinical: dict[str, str],
) -> ad.AnnData:
    gene_to_col = {gene: index for index, gene in enumerate(common_genes)}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        all_barcodes = split_fields(handle.readline())
        selected = [index for index, barcode in enumerate(all_barcodes) if barcode in keep_barcodes]
        barcodes = [all_barcodes[index] for index in selected]
        cell_rows: list[np.ndarray] = []
        gene_cols: list[np.ndarray] = []
        count_values: list[np.ndarray] = []

        for line in handle:
            gene = canonical_gene(feature_name(line))
            gene_col = gene_to_col.get(gene)
            if gene_col is None:
                continue
            fields = split_fields(line)
            values = np.fromiter((int(fields[index + 1]) for index in selected), dtype=np.int32, count=len(selected))
            nonzero = np.flatnonzero(values)
            if nonzero.size:
                cell_rows.append(nonzero.astype(np.int32, copy=False))
                gene_cols.append(np.full(nonzero.size, gene_col, dtype=np.int32))
                count_values.append(values[nonzero])

    rows = np.concatenate(cell_rows)
    cols = np.concatenate(gene_cols)
    values = np.concatenate(count_values)
    matrix = sparse.coo_matrix(
        (values, (rows, cols)), shape=(len(barcodes), len(common_genes)), dtype=np.int32
    ).tocsr()
    matrix.sum_duplicates()

    obs = pd.DataFrame(index=[f"{sample}_{barcode}" for barcode in barcodes])
    obs["sample"] = sample
    obs["condition"] = condition
    obs["barcode"] = barcodes
    for field in ["case", "diagnosis", "IDH", "MGMT", "chemotherapy", "radiotherapy"]:
        obs[field] = clinical.get(field, "")
    var = pd.DataFrame(index=pd.Index(common_genes, name="gene"))
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    print(json.dumps({"sample": sample, "cells": result.n_obs, "genes": result.n_vars, "nnz": int(result.X.nnz)}), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cell-qc", type=Path, required=True)
    parser.add_argument("--h5ad-out", type=Path, required=True)
    parser.add_argument("--doublet-summary-out", type=Path, required=True)
    parser.add_argument("--build-report-out", type=Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest)
    manifest = manifest.loc[manifest["strict_discovery"].astype(bool)].copy()
    qc = pd.read_csv(args.cell_qc)
    qc = qc.loc[qc["keep_basic"].astype(bool)].copy()

    matrix_paths: dict[str, Path] = {}
    for path in args.matrix_dir.glob("GSM*_TB*_mat.txt.gz"):
        sample = path.name.split("_", 2)[1]
        matrix_paths[sample] = path

    gene_lists = {sample: read_genes(matrix_paths[sample]) for sample in manifest["sample"]}
    common = set(next(iter(gene_lists.values())))
    for genes in gene_lists.values():
        common.intersection_update(genes)
    first_sample = manifest.iloc[0]["sample"]
    common_genes = [gene for gene in gene_lists[first_sample] if gene in common]

    objects: list[ad.AnnData] = []
    for record in manifest.to_dict("records"):
        sample = record["sample"]
        keep = set(qc.loc[qc["sample"] == sample, "barcode"])
        objects.append(
            build_sample(
                matrix_paths[sample], sample, record["condition"], keep, common_genes, record
            )
        )

    atlas = ad.concat(objects, join="inner", merge="same", index_unique=None)
    atlas.var_names_make_unique()
    sc.pp.calculate_qc_metrics(atlas, qc_vars=[], inplace=True, percent_top=None, log1p=True)
    sc.pp.scrublet(atlas, batch_key="sample", expected_doublet_rate=0.06, random_state=34)

    args.h5ad_out.parent.mkdir(parents=True, exist_ok=True)
    args.doublet_summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.build_report_out.parent.mkdir(parents=True, exist_ok=True)
    atlas.write_h5ad(args.h5ad_out, compression="gzip")

    doublets = (
        atlas.obs.groupby(["sample", "condition"], observed=True)
        .agg(cells=("predicted_doublet", "size"), predicted_doublets=("predicted_doublet", "sum"), median_score=("doublet_score", "median"))
        .reset_index()
    )
    doublets["predicted_doublet_rate"] = doublets["predicted_doublets"] / doublets["cells"]
    doublets.to_csv(args.doublet_summary_out, index=False, encoding="utf-8-sig")

    report = {
        "samples": int(atlas.obs["sample"].nunique()),
        "cells": int(atlas.n_obs),
        "common_genes": int(atlas.n_vars),
        "nnz": int(atlas.X.nnz),
        "predicted_doublets": int(atlas.obs["predicted_doublet"].sum()),
        "predicted_doublet_rate": float(atlas.obs["predicted_doublet"].mean()),
    }
    args.build_report_out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

