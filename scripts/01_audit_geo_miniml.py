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
import hashlib
import json
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path


NS = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}


def text(node: ET.Element | None) -> str:
    if node is None or node.text is None:
        return ""
    return " ".join(node.text.split())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_miniml(archive: Path) -> list[dict[str, str]]:
    with tarfile.open(archive, "r:gz") as tar:
        xml_member = next(member for member in tar.getmembers() if member.name.endswith(".xml"))
        stream = tar.extractfile(xml_member)
        if stream is None:
            raise RuntimeError(f"Cannot read {xml_member.name}")
        root = ET.parse(stream).getroot()

    records: list[dict[str, str]] = []
    for sample in root.findall("m:Sample", NS):
        characteristics: dict[str, list[str]] = {}
        for item in sample.findall("m:Channel/m:Characteristics", NS):
            tag = item.attrib.get("tag", "characteristic").strip().lower().replace(" ", "_")
            characteristics.setdefault(tag, []).append(text(item))

        record = {
            "accession": text(sample.find("m:Accession", NS)),
            "title": text(sample.find("m:Title", NS)),
            "source": text(sample.find("m:Channel/m:Source", NS)),
            "description": " | ".join(text(x) for x in sample.findall("m:Description", NS) if text(x)),
        }
        for key, values in characteristics.items():
            record[key] = " | ".join(values)
        records.append(record)
    return records


def audit_matrix(path: Path) -> dict[str, object]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
        first = handle.readline().rstrip("\n\r").split("\t")
        rows = 1
        for _ in handle:
            rows += 1
    return {
        "path": str(path),
        "header_fields": len(header),
        "first_data_fields": len(first),
        "rows": rows,
        "first_gene": first[0] if first else "",
        "sample_headers": header,
        "missing_gene_header": len(first) == len(header) + 1,
    }


def write_tsv(records: list[dict[str, str]], path: Path) -> None:
    fields = ["accession", "title", "source", "description"]
    extra = sorted({key for record in records for key in record if key not in fields})
    fields.extend(extra)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--miniml", type=Path, required=True)
    parser.add_argument("--matrix", type=Path)
    parser.add_argument("--metadata-out", type=Path, required=True)
    parser.add_argument("--audit-out", type=Path, required=True)
    args = parser.parse_args()

    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    args.audit_out.parent.mkdir(parents=True, exist_ok=True)
    records = parse_miniml(args.miniml)
    write_tsv(records, args.metadata_out)

    audit: dict[str, object] = {
        "miniml": str(args.miniml),
        "miniml_sha256": sha256(args.miniml),
        "sample_count": len(records),
        "metadata_columns": sorted({key for record in records for key in record}),
    }
    if args.matrix:
        audit["matrix_sha256"] = sha256(args.matrix)
        audit["matrix"] = audit_matrix(args.matrix)
    args.audit_out.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

