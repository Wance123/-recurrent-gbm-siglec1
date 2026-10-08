
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from revision_analysis import ROOT,REV,TAB,REFS,paired,bh,save
import pandas as pd,numpy as np
d=pd.read_csv(TAB/'S11_external_specimen_metrics.csv');rows=[];tests=[];values=[]
for dataset,a in d.groupby('dataset'):
 for gene in REFS:
  x=a[gene].dropna() if gene in a else pd.Series(dtype=float)
  rows.append({'dataset':dataset,'gene':gene,'n_measured_specimens':len(x),'n_nonzero_specimens':int((x!=0).sum()),'median_deposited_or_normalized_expression':x.median() if len(x) else np.nan,'unit':a.unit.iloc[0]})
 # Matched two-gene sensitivity is a post-audit diagnostic, not a registered primary metric.
 a=a.copy();a['SIGLEC1_two_gene_reference_difference']=a.SIGLEC1-a[['TYROBP','CSF1R']].mean(axis=1)
 w=a.pivot(index='patient',columns='condition',values='SIGLEC1_two_gene_reference_difference').dropna()[['primary','recurrent']].reset_index()
 w['dataset']=dataset;w['delta']=w.recurrent-w.primary;values.append(w)
 tests.append({'dataset':dataset,**paired(w,'SIGLEC1 minus mean(TYROBP, CSF1R)',a.unit.iloc[0])})
save(pd.DataFrame(rows),'S11_reference_detection_diagnostic');save(pd.concat(values),'S11_two_gene_reference_values')
r=save(bh(pd.DataFrame(tests)),'S11_two_gene_reference_tests')
print(r.to_string(index=False))
