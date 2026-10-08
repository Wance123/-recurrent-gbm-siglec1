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

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sample-qc", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    parser.add_argument("--table-out", type=Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest)
    qc = pd.read_csv(args.sample_qc)
    data = manifest.merge(qc, on=["sample", "condition"], suffixes=("_manifest", ""))
    strict = data.loc[data["strict_discovery"].astype(bool)].copy()
    strict["retention_rate"] = strict["basic_qc_pass"] / strict["raw_barcodes"]
    strict["log10_median_counts"] = np.log10(strict["median_counts"])

    metrics = ["raw_barcodes", "median_counts", "median_features", "median_percent_mt", "retention_rate"]
    group_tests = {}
    for metric in metrics:
        primary = strict.loc[strict["condition"] == "Primary", metric]
        recurrent = strict.loc[strict["condition"] == "Recurrent", metric]
        test = mannwhitneyu(primary, recurrent, alternative="two-sided")
        group_tests[metric] = {
            "primary_median": float(primary.median()),
            "recurrent_median": float(recurrent.median()),
            "u": float(test.statistic),
            "p_value": float(test.pvalue),
        }

    correlation = spearmanr(strict["median_counts"], strict["median_features"])
    flags: list[str] = []
    if (strict["basic_qc_pass"] < 300).any():
        flags.append("sample_below_300_cells")
    if (strict["retention_rate"] < 0.80).any():
        flags.append("sample_retention_below_80_percent")
    if (strict["median_features"] < 700).any():
        flags.append("sample_median_features_below_700")
    if (strict["median_percent_mt"] > 10).any():
        flags.append("sample_median_mt_above_10_percent")
    if group_tests["median_features"]["p_value"] < 0.05:
        flags.append("condition_associated_with_median_features")
    if group_tests["median_counts"]["p_value"] < 0.05:
        flags.append("condition_associated_with_median_counts")

    report = {
        "strict_samples": int(len(strict)),
        "primary_samples": int((strict["condition"] == "Primary").sum()),
        "recurrent_samples": int((strict["condition"] == "Recurrent").sum()),
        "strict_basic_qc_cells": int(strict["basic_qc_pass"].sum()),
        "sample_cell_range": [int(strict["basic_qc_pass"].min()), int(strict["basic_qc_pass"].max())],
        "median_feature_range": [float(strict["median_features"].min()), float(strict["median_features"].max())],
        "median_mt_range": [float(strict["median_percent_mt"].min()), float(strict["median_percent_mt"].max())],
        "retention_range": [float(strict["retention_rate"].min()), float(strict["retention_rate"].max())],
        "counts_features_spearman": {"rho": float(correlation.statistic), "p_value": float(correlation.pvalue)},
        "group_tests": group_tests,
        "flags": flags,
        "gate": "pass_with_covariate_control" if not flags else "review_flags",
    }

    args.report_out.parent.mkdir(parents=True, exist_ok=True)
    args.table_out.parent.mkdir(parents=True, exist_ok=True)
    args.report_out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    strict.to_csv(args.table_out, index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()

