
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from revision_analysis import *
import requests, shutil
OUT=REV/'submission_package';OUT.mkdir(exist_ok=True)
SRC=REV/'source_audit';SRC.mkdir(exist_ok=True)
pr={}
for acc in ['GSE339484','GSE200984','GSE139533','GSE254873','GSE254145','GSE174554']:
 try:
  r=requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',params={'acc':acc,'targ':'self','form':'text','view':'full'},proxies=pr,timeout=45)
  if r.ok and '^SERIES' in r.text:(SRC/(acc+'_GEO_record.txt')).write_text(r.text,encoding='utf-8')
 except Exception as e:print(acc,str(e)[:180])
try:
 r=requests.get('https://raw.githubusercontent.com/adithyakan/reconvolving_gbm/main/README.md',proxies=pr,timeout=45);(SRC/'GSE200984_author_repository_README.txt').write_text(r.text,encoding='utf-8')
 r=requests.get('https://drive.google.com/drive/folders/1wyV76FgQ5whQTqkm08YtQLIpul7OfWY0',proxies=pr,timeout=45,allow_redirects=False)
 (SRC/'GSE200984_annotation_access.json').write_text(json.dumps({'http_status':r.status_code,'redirect':r.headers.get('Location',''),'cell_label_concordance':'not computed; accessible repository contains scripts, annotations are in externally hosted processed objects'},indent=2),encoding='utf-8')
except Exception as e:print('ANNOTATION',str(e)[:200])
q=pd.read_csv(ROOT/'results/tables/GSE200984_scrublet_by_sample.csv')
with h5py.File(ROOT/'data/processed/GSE200984_strict_qc_scrublet.h5ad') as f:q['automatic_threshold']=[float(f['uns/scrublet/batches/'+s+'/threshold'][()]) for s in q['sample']]
q['percent_removed']=100*q.predicted_doublet_rate;save(q,'S01_Scrublet_QC')
gs=[]
for n,genes in PROGRAMS.items():
 for g in genes:gs.append(dict(program=n,gene=g,definition='Investigator-selected literature-informed transcript score; not a validated functional assay',context_references='7-10,12; SIGLEC1 itself excluded'))
for n,genes in [('Myeloid reference',REFS),('Microglia marker proxy',MG),('Macrophage marker proxy',MDM)]:
 for g in genes:gs.append(dict(program=n,gene=g,definition='Sensitivity reference or operational marker proxy; not lineage tracing',context_references='7-10'))
save(pd.DataFrame(gs),'S07_gene_sets')
for file,name in [('GSE200984_clean_myeloid_subcluster_statistics.csv','S03_cluster_fraction_tests'),('GSE200984_clean_myeloid_subcluster_patient_composition.csv','S03_cluster_fraction_values')]:
 save(pd.read_csv(ROOT/'results/tables'/file),name)
for file,name in [('leave_one_patient_out.csv','S05_leave_one_out'),('leave_one_patient_out_summary.csv','S05_leave_one_out_summary'),('patient_ecology_correlations.csv','S08_original_correlations'),('patient_ecology_metrics.csv','S08_original_patient_metrics')]:
 save(pd.read_csv(ROOT/'further/results/reinforcement'/file),name)
# Explicitly document the earlier restricted tumor-like selection and its patient distribution.
with h5py.File(ROOT/'data/processed/GSE200984_strict_basic_atlas_cnv.h5ad') as f:a=frame(f['obs'])
a['cluster']=a.leiden_0_6.astype(str);c=a.groupby(['cluster','sample']).size().rename('cells').reset_index();c['cluster_cells']=c.groupby('cluster').cells.transform('sum');c['patient_share']=c.cells/c.cluster_cells
z=c.groupby('cluster').agg(cells=('cells','sum'),patients=('sample','nunique'),max_patient_share=('patient_share','max')).reset_index();z['included_earlier_MES_summary']=z.cluster.isin(['0','1','2','3','5']);z['patient_dominant_over_50_percent']=z.max_patient_share>.5
save(z,'S08_atlas_cluster_selection');save(c,'S08_atlas_patient_distribution')
# GSE139533 alternative patient aggregation: pool raw regional counts before CPM.
b=pd.read_csv(ROOT/'data/raw/GSE139533/GSE139533_mergede_protein_coding_counts.tab.gz',sep='\t',index_col=0).groupby(level=0).sum();man=pd.read_csv(TAB/'S02_GSE139533_sample_selection.csv');rr=[];cov=[]
for (pt,st),a in man.groupby(['patient id','stage']):
 c=b[a.matrix_column].sum(axis=1);e=np.log2(c/c.sum()*1e6+1);rr.append(dict(patient=pt,condition=st.lower(),SIGLEC1=e['SIGLEC1'],CD68=e['CD68'],SIGLEC1_CD68_log_difference=e['SIGLEC1']-e['CD68'],SIGLEC1_reference_log_difference=e['SIGLEC1']-e[REFS].mean()))
e=save(pd.DataFrame(rr),'S11_GSE139533_pooled_counts_values');tt=[]
for mt in ['SIGLEC1','CD68','SIGLEC1_CD68_log_difference','SIGLEC1_reference_log_difference']:tt.append(paired(e.pivot(index='patient',columns='condition',values=mt).dropna().reset_index(),mt,'log2(CPM + 1) after pooling raw regional counts'))
save(bh(pd.DataFrame(tt)),'S11_GSE139533_aggregation_tests')
# Preserve the exact exploratory perturbation outputs without changing its negative result.
for d in ['sctenifoldknk','stability']:
 p=ROOT/'further/results'/d
 if p.exists():
  for f in p.rglob('*.csv'):
   name='S13_'+d+'_'+f.stem
   save(pd.read_csv(f),name)
summ={
 'GSE200984':{'n_primary':5,'n_recurrent':6,'paired':False,'post_QC_nuclei':29038,'myeloid_before_exclusion':3768,'myeloid_after_exclusion':3423},
 'GSE339484':{'GEO_status':'Public on Jul 22 2026; live GEO full record overrides stale SRA private notice','n_paired':11,'excluded_IDH_mutant':['Patient02','Patient09','Patient12']},
 'GSE139533':{'n_pairs':13,'n_tissue_samples':len(man),'n_primary_tissue_samples':int((man.stage.str.lower()=='primary').sum()),'n_recurrent_tissue_samples':int((man.stage.str.lower()=='recurrent').sum())},
 'GSE174554':{'eligible_original_pairs':34,'eligible_original_specimens':68,'author_labels':'GBM-CARE-WT data/diaz_cell_md_20240820.RDS and diaz_clinical_snrna_metadata_20240106.txt'},
 'not_fully_addressed':['GSE200984 original-label concordance not computed because processed annotation host requires authentication in this environment','Public analysis code not released under author instruction','No CD169 protein or functional validation','No spatial CD64-negative controls or per-compartment area supplied in processed matrix']}
ext=pd.read_csv(TAB/'S12_external_eligibility.csv');proc=ext[ext.status=='processed'];pairedv=pd.read_csv(TAB/'S12_external_paired_values.csv');pairids=pairedv[(pairedv.origin_proxy=='All myeloid')&(pairedv.metric=='SIGLEC1')].patient.astype(str)
summ['GSE174554'].update(processed_specimens=len(proc),primary_specimens=int((proc.condition=='primary').sum()),recurrent_specimens=int((proc.condition=='recurrent').sum()),paired_patients=17,myeloid_nuclei=int(proc.myeloid_nuclei.sum()),paired_myeloid_nuclei=int(proc[proc.patient.astype(str).isin(pairids)].myeloid_nuclei.sum()),excluded_incomplete=int(ext.status.str.contains('incomplete').sum()),excluded_small=int(ext.status.str.contains('fewer').sum()))
(REV/'analysis_summary.json').write_text(json.dumps(summ,indent=2),encoding='utf-8');print(json.dumps(summ,indent=2))
