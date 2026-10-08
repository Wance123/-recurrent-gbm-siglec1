
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from revision_analysis import ROOT,REFS
import gzip,re,pandas as pd,json
b=pd.read_csv(ROOT/'data/raw/GSE339484/GSE339484_counts_All_Samples.txt.gz',sep='\t',index_col=0)
b.index=b.index.str.split('.').str[0]
rows=[]
with gzip.open(ROOT/'data/reference/gencode.v50.basic.annotation.gtf.gz','rt') as f:
 for l in f:
  if '\tgene\t' not in l:continue
  name=re.search(r'gene_name "([^"]+)"',l);gid=re.search(r'gene_id "([^"]+)"',l)
  if name and name[1] in REFS+['SIGLEC1','CD68']:
   ident=gid[1].split('.')[0];rows.append({'gene':name[1],'ensembl':ident,'present':ident in b.index,'total_counts':float(b.loc[ident].sum()) if ident in b.index else None})
print(json.dumps(rows,indent=2));print('Canonical rows:');print(b.loc[b.index.isin(['ENSG00000204472','ENSG00000204482'])].to_string())
(ROOT/'Revision/source_audit/GSE339484_reference_gene_mapping.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
