
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import csv
import numpy as np
from scipy.io import mmread

base = Path(str(project_path('ROOT') / 'further'))
genes_of_interest = {
    "SIGLEC1", "STAT1", "IRF7", "IFNAR1", "FCGR1A", "CCL3",
    "KCNMB1", "FAM177B", "MBP", "PLP1", "CACNA1A", "CACNA1D",
    "CACNB4", "FAM118A", "NAV3", "IFI44L", "XIST"
}

rows = []
for condition in ("primary", "recurrent"):
    matrix = mmread(base / "input" / f"{condition}_genes_by_cells.mtx").tocsr()
    genes = [x.strip() for x in (base / "input" / f"{condition}_genes.tsv").read_text(encoding="utf-8").splitlines()]
    cells_file = base / "input" / f"{condition}_cells.tsv"
    cells = [x.strip() for x in cells_file.read_text(encoding="utf-8").splitlines()] if cells_file.exists() else [str(i) for i in range(matrix.shape[1])]
    for gene in sorted(genes_of_interest):
        if gene not in genes:
            continue
        values = matrix[genes.index(gene), :].toarray().ravel()
        rows.append({
            "condition": condition.capitalize(),
            "gene": gene,
            "n_cells": matrix.shape[1],
            "detected_cells": int((values > 0).sum()),
            "detected_fraction": float((values > 0).mean()),
            "mean_all_cells": float(values.mean()),
            "mean_detected_cells": float(values[values > 0].mean()) if np.any(values > 0) else 0.0,
        })

out = base / "results" / "sctenifoldknk" / "signal_gene_prevalence.csv"
with out.open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(out)
