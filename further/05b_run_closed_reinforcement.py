
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path

source_path = Path(str(project_path('ROOT') / 'further/05_closed_reinforcement.py'))
source = source_path.read_text(encoding="utf-8")
source = source.replace(
    'return a.raw.X.tocsr() if sparse.issparse(a.raw.X) else sparse.csr_matrix(a.raw.X)',
    'x = a.raw.X if a.raw is not None else a.X\n    return x.tocsr() if sparse.issparse(x) else sparse.csr_matrix(x)',
).replace(
    'return pd.Index(a.raw.var_names)',
    'return pd.Index(a.raw.var_names if a.raw is not None else a.var_names)',
)
exec(compile(source, str(source_path), "exec"), {"__name__": "__main__"})
