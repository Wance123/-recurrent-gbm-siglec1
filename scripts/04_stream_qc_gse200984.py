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

import numpy as np


def split_fields(line: str) -> list[str]:
    line = line.rstrip("\r\n")
    return line.split("\t") if "\t" in line else shlex.split(line)


def qc_file(path: Path, sample: str, condition: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        barcodes = split_fields(handle.readline())
        n = len(barcodes)
        umi = np.zeros(n, dtype=np.int64)
        genes = np.zeros(n, dtype=np.int32)
        mito = np.zeros(n, dtype=np.int64)
        feature_rows = 0

        for line in handle:
            fields = split_fields(line)
            if len(fields) != n + 1:
                raise ValueError(f"{path.name}: expected {n + 1} fields, found {len(fields)}")
            gene = fields[0].strip('"')
            values = np.fromiter((int(value) for value in fields[1:]), dtype=np.int64, count=n)
            umi += values
            genes += values > 0
            if gene.upper().startswith("MT-"):
                mito += values
            feature_rows += 1

    mito_pct = np.divide(mito * 100.0, umi, out=np.zeros(n, dtype=float), where=umi > 0)
    keep_basic = (genes >= 400) & (mito_pct < 20) & (umi > 0)
    records = [
        {
            "sample": sample,
            "condition": condition,
            "barcode": barcode,
            "nCount_RNA": int(umi[i]),
            "nFeature_RNA": int(genes[i]),
            "percent_mt": float(mito_pct[i]),
            "keep_basic": bool(keep_basic[i]),
        }
        for i, barcode in enumerate(barcodes)
    ]
    summary = {
        "sample": sample,
        "condition": condition,
        "file": path.name,
        "feature_rows": feature_rows,
        "raw_barcodes": n,
        "nonzero_barcodes": int((umi > 0).sum()),
        "basic_qc_pass": int(keep_basic.sum()),
        "median_counts": float(np.median(umi)),
        "median_features": float(np.median(genes)),
        "median_percent_mt": float(np.median(mito_pct)),
    }
    return records, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--cells-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    with args.metadata.open("r", encoding="utf-8-sig", newline="") as handle:
        metadata = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}

    all_cells: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for path in sorted(args.matrix_dir.glob("GSM*_TB*_mat.txt.gz")):
        accession = path.name.split("_", 1)[0]
        meta = metadata[accession]
        records, summary = qc_file(path, meta["title"], meta["disease_state"])
        all_cells.extend(records)
        summaries.append(summary)
        print(json.dumps(summary, ensure_ascii=False), flush=True)

    args.cells_out.parent.mkdir(parents=True, exist_ok=True)
    with args.cells_out.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_cells[0]))
        writer.writeheader()
        writer.writerows(all_cells)

    with args.summary_out.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    overall = {
        "samples": len(summaries),
        "raw_barcodes": sum(int(row["raw_barcodes"]) for row in summaries),
        "basic_qc_pass": sum(int(row["basic_qc_pass"]) for row in summaries),
        "by_condition": {},
    }
    for condition in sorted({str(row["condition"]) for row in summaries}):
        subset = [row for row in summaries if row["condition"] == condition]
        overall["by_condition"][condition] = {
            "samples": len(subset),
            "raw_barcodes": sum(int(row["raw_barcodes"]) for row in subset),
            "basic_qc_pass": sum(int(row["basic_qc_pass"]) for row in subset),
        }
    args.json_out.write_text(json.dumps(overall, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

