
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import sys,json
import h5py,numpy as np,pandas as pd
from revision_analysis import ROOT,REV,TAB,frame,matrix,gv,compare,bh,save

# Confirm that the 400-feature threshold was applied before taking the common-gene intersection.
rawqc=pd.read_csv(ROOT/'data/processed/GSE200984_cell_qc_corrected.csv')
with h5py.File(ROOT/'data/processed/GSE200984_strict_qc_scrublet.h5ad') as f:
    obs=frame(f['obs']).reset_index()
merged=obs.merge(rawqc[['sample','barcode','nCount_RNA','nFeature_RNA','percent_mt','keep_basic']],on=['sample','barcode'],how='left',validate='one_to_one')
assert len(merged)==29065
assert merged.keep_basic.notna().all()
assert merged.keep_basic.all() and (merged.nFeature_RNA>=400).all() and (merged.percent_mt<20).all()
q=merged.groupby('sample').agg(analyzed_pre_doublet_nuclei=('barcode','size'),predicted_doublets=('predicted_doublet','sum'),original_min_features=('nFeature_RNA','min'),common_gene_min_features=('n_genes_by_counts','min'),original_max_mito_percent=('percent_mt','max')).reset_index()
q['analyzed_post_doublet_nuclei']=q.analyzed_pre_doublet_nuclei-q.predicted_doublets
save(q,'S01_QC_provenance')
clinical=pd.read_csv(TAB/'S01_discovery_clinical.csv').rename(columns={'basic_qc_pass':'legacy_manifest_basic_qc_pass'})
for k in q.columns:
    if k!='sample' and k in clinical:clinical=clinical.drop(columns=k)
save(clinical.merge(q,on='sample',validate='one_to_one'),'S01_discovery_clinical')
corrected=rawqc[rawqc.keep_basic].groupby('sample').size()
assert all(corrected[s]==n for s,n in zip(q['sample'],q.analyzed_pre_doublet_nuclei))

# Additional threshold sensitivity on the processed common-gene counts.
with h5py.File(ROOT/'data/processed/GSE200984_myeloid_reclustered.h5ad') as f:
    o=frame(f['obs']);genes=pd.Index([v.decode() if isinstance(v,bytes) else v for v in f['var/gene'][:]]);x=matrix(f['raw/X'])
o['SIGLEC1']=gv(x,genes,'SIGLEC1');o['CD68']=gv(x,genes,'CD68');o['Detection']=o.SIGLEC1>0
d=o[(o.myeloid_leiden.astype(str)!='3')&(o.n_genes_by_counts>=400)]
a=d.groupby(['sample','condition']).agg(SIGLEC1=('SIGLEC1','mean'),CD68=('CD68','mean'),Detection=('Detection','mean'),n_cells=('SIGLEC1','size')).reset_index()
a['Log expression difference']=a.SIGLEC1-a.CD68
save(a,'S05_common_gene_400_patient_metrics')
tests=bh(pd.DataFrame([compare(a,m) for m in ['SIGLEC1','Detection','Log expression difference']]))
save(tests,'S05_common_gene_400_tests')
assert (tests.q>=.05).all()
units={'SIGLEC1':'Mean ln(1 + counts normalized to 10,000)','CD68':'Mean ln(1 + counts normalized to 10,000)','Detection':'Fraction of nuclei with a nonzero SIGLEC1 count','Log_expression_difference':'Difference of mean natural-log expression: SIGLEC1 minus CD68','myeloid_fraction':'Origin-category myeloid nuclei divided by all QC specimen nuclei'}
for stem in ['S12_external_tests','S12_external_paired_values']:
    d=pd.read_csv(TAB/(stem+'.csv'));d['unit']=d.metric.map(units);assert d.unit.notna().all();save(d,stem)
report={'discovery_pre_doublet':len(merged),'corrected_raw_QC_all_pass':True,'raw_feature_min':int(merged.nFeature_RNA.min()),'common_gene_feature_min':int(merged.n_genes_by_counts.min()),'legacy_manifest_difference':73,'difference_sample':'TB4416','explanation':'The legacy manifest predates corrected mitochondrial-gene QC; the corrected cell QC and reconstructed matrix agree exactly. Reduced feature counts reflect the intersection of genes across specimens.','common_gene_400_sensitivity_all_q_over_0_05':True,'sensitivity_myeloid_nuclei':int(a.n_cells.sum())}
(REV/'final_scientific_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2));print(tests.to_string(index=False))
