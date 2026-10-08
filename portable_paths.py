"""Explicit paths for public raw data; no network or credentials on import."""
from pathlib import Path
import os
ROOT=Path(os.environ.get('JOB34_PROJECT_ROOT',Path(__file__).resolve().parent))
DEFAULTS={
 'ROOT':ROOT,
 'GSE174554_METADATA':ROOT/'external_inputs/GSE174554',
 'GSE174554_RAW':ROOT/'external_inputs/GSE174554/raw',
 'CARE_METADATA':ROOT/'external_inputs/CARE',
 'CARE_RAW':ROOT/'external_inputs/CARE/raw',
 'CGGA693':ROOT/'external_inputs/CGGA693',
 'CGGA325':ROOT/'external_inputs/CGGA325',
 'AUTHOR_WORKSPACE':ROOT/'external_inputs/GSE200984_author_workspace',
}
def project_path(name):
    return Path(os.environ.get('JOB34_'+name,DEFAULTS[name]))
