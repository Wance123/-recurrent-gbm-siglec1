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
import tarfile
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tar", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    args = parser.parse_args()

    with args.metadata.open("r", encoding="utf-8-sig", newline="") as handle:
        metadata = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}

    results: list[dict[str, object]] = []
    with tarfile.open(args.tar, "r") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".gz"):
                continue
            raw = archive.extractfile(member)
            if raw is None:
                raise RuntimeError(f"Cannot read {member.name}")
            accession = member.name.split("_", 1)[0]
            with gzip.open(raw, "rt", encoding="utf-8", errors="replace") as handle:
                header = handle.readline().rstrip("\r\n").split("\t")
                first = handle.readline().rstrip("\r\n").split("\t")
                row_count = 1 + sum(1 for _ in handle)
            meta = metadata.get(accession, {})
            results.append(
                {
                    "member": member.name,
                    "compressed_bytes": member.size,
                    "accession": accession,
                    "title": meta.get("title", ""),
                    "disease_state": meta.get("disease_state", ""),
                    "molecule_type": meta.get("molecule_type", ""),
                    "header_fields": len(header),
                    "first_data_fields": len(first),
                    "data_rows": row_count,
                    "first_row_id": first[0] if first else "",
                    "estimated_observations": max(len(header), len(first)) - 1,
                    "field_count_consistent": len(header) == len(first),
                }
            )

    summary = {
        "tar": str(args.tar),
        "member_count": len(results),
        "all_members_have_metadata": all(item["title"] for item in results),
        "all_field_counts_consistent": all(item["field_count_consistent"] for item in results),
        "members": results,
    }
    args.report_out.parent.mkdir(parents=True, exist_ok=True)
    args.report_out.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

