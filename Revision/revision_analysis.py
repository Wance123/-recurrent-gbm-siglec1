"""Local-only analysis for the Discover Oncology major revision."""

# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import gzip, hashlib, json, re
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
from statsmodels.stats.multitest import multipletests

ROOT = Path(str(project_path('ROOT')))
REV = ROOT/'Revision'
TAB = REV/'analysis_tables'
TAB.mkdir(exist_ok=True)
RNG = np.random.default_rng(341008)
PROGRAMS = {
 'Interferon response':['ISG15','IFIT1','IFIT2','IFIT3','IFI6','IFI27','IFI44L','MX1','OAS1','STAT1','IRF7'],
 'Phagolysosome':['LAMP1','CTSD','CTSB','FCGR1A','TYROBP','AIF1','LILRB1','CTSS'],
 'Lipid and complement':['APOE','LPL','TREM2','C1QA','C1QB','C1QC','GPNMB','FABP5'],
 'Antigen presentation':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','B2M'],
 'Immunoregulatory':['IL10','TGFB1','CD274','ARG1','MRC1','MSR1','LGALS3'],
 'Immediate inflammatory':['CCL3','CCL4','IL1B','TNF','NFKBIA','CXCL8','PTGS2','FOS','JUN']}
REFS = ['AIF1','LST1','TYROBP','CSF1R']
MG = ['P2RY12','TMEM119','CX3CR1','SALL1']
MDM = ['CD163','MRC1','MSR1','FCN1','S100A8','S100A9']

def save(d,name):
 d.to_csv(TAB/(name+'.csv'),index=False,encoding='utf-8-sig')
 return d

def decode(a):
 return np.asarray([x.decode() if isinstance(x,bytes) else x for x in a])

def vector(g):
 if isinstance(g,h5py.Dataset): return decode(g[:])
 cats=decode(g['categories'][:]); codes=g['codes'][:]
 return np.asarray([cats[i] if i>=0 else None for i in codes])

def frame(g):
 cols={k:vector(v) for k,v in g.items()}
 idx=g.attrs.get('_index','_index')
 if isinstance(idx,bytes): idx=idx.decode()
 return pd.DataFrame(cols).set_index(idx)

def matrix(g):
 return sparse.csr_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g.attrs['shape']))

def norm(c):
 lib=np.asarray(c.sum(axis=1)).ravel()
 x=c.astype(float).multiply((1e4/np.maximum(lib,1))[:,None]).tocsr()
 x.data=np.log1p(x.data)
 return x

def gv(x,genes,name):
 if name not in genes: return np.full(x.shape[0],np.nan)
 return np.asarray(x[:,genes.get_loc(name)].toarray()).ravel()

def score(x,genes,gs):
 ix=genes.get_indexer([g for g in gs if g in genes])
 return np.asarray(x[:,ix].mean(axis=1)).ravel() if len(ix) else np.full(x.shape[0],np.nan)

def wil(v):
 v=np.asarray(v,dtype=float); v=v[np.isfinite(v)]
 return float(stats.wilcoxon(v,alternative='two-sided',method='auto').pvalue) if len(v)>1 and np.any(v!=0) else 1.

def paired(d,metric,unit):
 p=d['primary'].to_numpy(float);r=d['recurrent'].to_numpy(float);v=r-p;n=len(v)
 boot=np.median(RNG.choice(v,size=(10000,n),replace=True),axis=1)
 return dict(metric=metric,unit=unit,n_pairs=n,increased=int((v>0).sum()),decreased=int((v<0).sum()),median_delta=float(np.median(v)),ci_low=float(np.quantile(boot,.025)),ci_high=float(np.quantile(boot,.975)),mean_delta=float(np.mean(v)),paired_dz=float(np.mean(v)/np.std(v,ddof=1)) if np.std(v,ddof=1)>0 else np.nan,p=wil(v))

def bh(d,group=None):
 d=d.copy()
 if group: d['q']=d.groupby(group,observed=True)['p'].transform(lambda x:multipletests(x.fillna(1),method='fdr_bh')[1])
 else: d['q']=multipletests(d.p.fillna(1),method='fdr_bh')[1]
 return d

def compare(d,metric):
 p=d.loc[d.condition=='Primary',metric].dropna().to_numpy(float);r=d.loc[d.condition=='Recurrent',metric].dropna().to_numpy(float)
 b=np.median(RNG.choice(r,size=(10000,len(r)),replace=True),axis=1)-np.median(RNG.choice(p,size=(10000,len(p)),replace=True),axis=1)
 return dict(metric=metric,n_primary=len(p),n_recurrent=len(r),primary_median=np.median(p),recurrent_median=np.median(r),delta=np.median(r)-np.median(p),ci_low=np.quantile(b,.025),ci_high=np.quantile(b,.975),cliff_delta=np.mean(np.sign(r[:,None]-p[None,:])),p=float(stats.mannwhitneyu(r,p).pvalue))

def discovery():
 f=h5py.File(ROOT/'data/processed/GSE200984_myeloid_reclustered.h5ad')
 o=frame(f['obs']);genes=pd.Index(vector(f['var/gene']));counts=matrix(f['layers/counts']);x=matrix(f['raw/X']);xy=f['obsm/X_umap'][:]; f.close()
 o['cluster']=o.myeloid_leiden.astype(str)
 o['SIGLEC1']=gv(x,genes,'SIGLEC1');o['CD68']=gv(x,genes,'CD68');o['detected']=o.SIGLEC1>0
 for label,gs in [('MG_score',MG),('MDM_score',MDM)]+list(PROGRAMS.items()):o[label]=score(x,genes,gs)
 gap=o.MG_score-o.MDM_score
 # Relative marker scores are an operational sensitivity analysis, not lineage tracing.
 o['origin_proxy']=np.where(gap>.1,'Microglia-like',np.where(gap<-.1,'Macrophage-like','Ambiguous'))
 markers=pd.read_csv(ROOT/'results/tables/GSE200984_myeloid_subcluster_markers.csv')
 save(markers,'S04_full_cluster_markers')
 composition=o.groupby(['sample','condition','cluster']).size().unstack(fill_value=0).stack().rename('cells').reset_index()
 composition['fraction_of_sample_myeloid']=composition.cells/composition.groupby('sample').cells.transform('sum')
 composition['fraction_of_cluster']=composition.cells/composition.groupby('cluster').cells.transform('sum')
 save(composition,'S03_patient_cluster_distribution')
 summary=o.groupby('cluster').agg(cells=('sample','size'),patients=('sample','nunique'),median_doublet_score=('doublet_score','median'),median_CNV_RMS=('cnv_rms','median'),CNV_high_fraction=('cnv_high','mean'),median_UMI=('total_counts','median'),median_genes=('n_genes_by_counts','median'),median_SIGLEC1=('SIGLEC1','median')).reset_index()
 save(summary,'S04_cluster_QC')
 save(o[['sample','condition','cluster','doublet_score','cnv_rms','cnv_high','total_counts','n_genes_by_counts','SIGLEC1','CD68','origin_proxy']].reset_index(),'S04_cell_QC')
 pats=[];test=[]
 for mode in ['Exclude cluster 3','Include cluster 3']:
  d=o[o.cluster!='3'] if mode.startswith('Exclude') else o
  a=d.groupby(['sample','condition']).agg(SIGLEC1=('SIGLEC1','mean'),CD68=('CD68','mean'),Detection=('detected','mean'),n_cells=('sample','size')).reset_index()
  a['Log expression difference']=a.SIGLEC1-a.CD68;a['filter']=mode;pats.append(a)
  for m in ['SIGLEC1','Detection','Log expression difference']: test.append(dict(filter=mode,**compare(a,m)))
 save(pd.concat(pats),'S05_discovery_patient_metrics');save(bh(pd.DataFrame(test),'filter'),'S05_discovery_inclusion_tests')
 d=o[o.cluster!='3'];origin=[]
 for label,g in d.groupby('origin_proxy'):
  a=g.groupby(['sample','condition']).agg(SIGLEC1=('SIGLEC1','mean'),Detection=('detected','mean'),n_cells=('sample','size')).reset_index();a=a[a.n_cells>=20];a['origin_proxy']=label;origin.append(a)
 save(pd.concat(origin),'S06_origin_patient_metrics')
 rr=[]
 for label,g in pd.concat(origin).groupby('origin_proxy'):
  for mt in ['SIGLEC1','Detection']:
   if g.condition.nunique()==2:rr.append(dict(origin_proxy=label,**compare(g,mt)))
 save(bh(pd.DataFrame(rr)),'S06_origin_tests')
 # Patient-level exploratory scores and detection-sensitivity analyses.
 rows=[];continuous=[]
 for (sample,condition),a in d.groupby(['sample','condition']):
  for name in PROGRAMS:
   if np.std(a.SIGLEC1)>0:
    rho=stats.spearmanr(a.SIGLEC1,a[name]).statistic
    z=np.column_stack([np.ones(len(a)),np.log1p(a.total_counts.astype(float)),np.log1p(a.n_genes_by_counts.astype(float))])
    xr=stats.rankdata(a.SIGLEC1);yr=stats.rankdata(a[name]);xr-=z@np.linalg.lstsq(z,xr,rcond=None)[0];yr-=z@np.linalg.lstsq(z,yr,rcond=None)[0]
    adjusted=np.corrcoef(xr,yr)[0,1]
    continuous.append(dict(sample=sample,condition=condition,program=name,rho=rho,depth_adjusted_rank_correlation=adjusted,n_cells=len(a)))
   for threshold in [0,1.]:
    pos=a.SIGLEC1>threshold;neg=~pos
    if pos.sum()<5 or neg.sum()<5:continue
    y1=a.loc[pos,name].to_numpy();y0=a.loc[neg,name].to_numpy();delta=np.mean(y1)-np.mean(y0)
    sd=np.sqrt(((len(y1)-1)*np.var(y1,ddof=1)+(len(y0)-1)*np.var(y0,ddof=1))/(len(y1)+len(y0)-2))
    rows.append(dict(sample=sample,condition=condition,program=name,threshold=threshold,n_positive=int(pos.sum()),n_negative=int(neg.sum()),delta=delta,standardized_difference=delta/sd if sd>0 else np.nan))
 ds=save(pd.DataFrame(rows),'S07_program_patient_effects');cs=save(pd.DataFrame(continuous),'S07_continuous_correlations')
 tests=[]
 for (c,n,t),a in ds.groupby(['condition','program','threshold']):tests.append(dict(condition=c,program=n,threshold=t,n_patients=len(a),positive=int((a.delta>0).sum()),median_delta=a.delta.median(),median_standardized_difference=a.standardized_difference.median(),p=wil(a.delta)))
 save(bh(pd.DataFrame(tests),['condition','threshold']),'S07_program_tests')
 tests=[]
 for (c,n),a in cs.groupby(['condition','program']):
  for mt in ['rho','depth_adjusted_rank_correlation']:tests.append(dict(condition=c,program=n,metric=mt,n_patients=len(a),median_correlation=a[mt].median(),p=wil(a[mt])))
 save(bh(pd.DataFrame(tests),['condition','metric']),'S07_continuous_tests')
 # Compositional sensitivity using all five annotated atlas lineages.
 f=h5py.File(ROOT/'data/processed/GSE200984_strict_basic_atlas_cnv.h5ad');ao=frame(f['obs']);f.close()
 cl=ao.leiden_0_6.astype(str)
 cmap={str(k):('Myeloid' if k==9 else 'T_NK' if k==12 else 'Oligodendrocyte' if k in [10,11,13] else 'Neuron' if k==4 else 'Tumor') for k in range(18)}
 ao['lineage']=cl.map(cmap)
 counts_by_patient=ao.groupby(['sample','condition','lineage']).size().unstack(fill_value=0)
 a=counts_by_patient.copy();v=a.to_numpy(float)+.5;frac=a.to_numpy()/a.sum(axis=1).to_numpy()[:,None];clr=np.log(v)-np.log(v).mean(axis=1)[:,None]
 for j,col in enumerate(a.columns):a[col+'_fraction']=frac[:,j];a[col+'_CLR']=clr[:,j]
 a=a.reset_index();p=pd.concat(pats);p=p[p['filter']=='Exclude cluster 3'];a=a.merge(p,on=['sample','condition'])
 save(a,'S08_patient_composition_CLR')
 rr=[]
 for mt in ['Myeloid_fraction','Tumor_fraction','Myeloid_CLR','Tumor_CLR','T_NK_CLR']:
  rho,pv=stats.spearmanr(a.SIGLEC1,a[mt]);rr.append(dict(outcome=mt,n_patients=len(a),rho=rho,p=pv))
 # Partial rank association conditional on myeloid fraction (descriptive with n=11).
 z=np.column_stack([np.ones(len(a)),stats.rankdata(a.Myeloid_fraction)]);xx=stats.rankdata(a.SIGLEC1);yy=stats.rankdata(a.Tumor_CLR)
 xx-=z@np.linalg.lstsq(z,xx,rcond=None)[0];yy-=z@np.linalg.lstsq(z,yy,rcond=None)[0];r=np.corrcoef(xx,yy)[0,1];pv=2*stats.t.sf(abs(r*np.sqrt((len(a)-3)/(1-r*r))),df=len(a)-3)
 rr.append(dict(outcome='Tumor_CLR_adjusted_myeloid_fraction',n_patients=len(a),rho=r,p=pv))
 save(bh(pd.DataFrame(rr)),'S08_composition_correlations')
 clinical=pd.read_csv(ROOT/'data/metadata/GSE228500_patient_metadata.tsv',sep='\t');manifest=pd.read_csv(ROOT/'data/metadata/GSE200984_discovery_manifest.csv');manifest=manifest[manifest.strict_discovery]
 clinical['sample']='TB'+clinical['TB #'].astype(str)
 clinical_out=manifest.merge(clinical,on='sample',how='left',suffixes=('','_published')).rename(columns={'basic_qc_pass':'legacy_manifest_basic_qc_pass'})
 qc=pd.read_csv(ROOT/'results/tables/GSE200984_scrublet_by_sample.csv')
 clinical_out=clinical_out.merge(qc[['sample','cells','predicted_doublets']],on='sample',how='left').rename(columns={'cells':'analyzed_pre_doublet_nuclei'})
 clinical_out['analyzed_post_doublet_nuclei']=clinical_out.analyzed_pre_doublet_nuclei-clinical_out.predicted_doublets
 save(clinical_out,'S01_discovery_clinical')
 # Model-based minimum detectable effect, normal shift under 5+6 exact MW design.
 from itertools import combinations
 combos=list(combinations(range(11),5));sums=np.array([sum(np.array(c)+1) for c in combos])-15
 probs=np.array([min(1,2*min(np.mean(sums<=u),np.mean(sums>=u))) for u in range(31)])
 out=[]
 for effect in np.arange(0,3.01,.1):
  p=RNG.normal(0,1,size=(20000,5));r=RNG.normal(effect,1,size=(20000,6));u=(r[:,:,None]>p[:,None,:]).sum(axis=(1,2));power=np.mean(probs[u]<=.05);out.append(dict(standardized_normal_shift=effect,simulations=20000,power=power,alpha=.05,test='two-sided exact Mann-Whitney; normal equal-variance shift model'))
 save(pd.DataFrame(out),'S09_power_sensitivity')
 return o,xy

def bulk():
 ids={}
 wanted=set(['SIGLEC1','CD68']+REFS)
 with gzip.open(ROOT/'data/reference/gencode.v50.basic.annotation.gtf.gz','rt') as f:
  for line in f:
   if '\tgene\t' not in line:continue
   g=re.search(r'gene_name "([^"]+)"',line);i=re.search(r'gene_id "([^"]+)"',line)
   if g and i and g[1] in wanted:ids[i[1].split('.')[0]]=g[1]
 b=pd.read_csv(ROOT/'data/raw/GSE339484/GSE339484_counts_All_Samples.txt.gz',sep='\t',index_col=0);b.index=b.index.str.split('.').str[0];b=b.rename(index=ids);b=b.groupby(level=0).sum();b=np.log2(b.div(b.sum(axis=0),axis=1)*1e6+1)
 e=b.loc[b.index.intersection(wanted)].T;e['patient']=e.index.str.split('_').str[0];e['condition']=e.index.str.split('_').str[1].str.lower();e=e[~e.patient.isin(['Patient02','Patient09','Patient12'])];e['dataset']='GSE339484';e['unit']='log2(CPM + 1)'
 b=pd.read_csv(ROOT/'data/raw/GSE139533/GSE139533_mergede_protein_coding_counts.tab.gz',sep='\t',index_col=0).groupby(level=0).sum();b=np.log2(b.div(b.sum(axis=0),axis=1)*1e6+1)
 man=pd.read_csv(ROOT/'data/metadata/GSE139533_paired_tissue_manifest.csv');man['matrix_column']=man.title.str.replace('_ZK','',regex=False);man=man[man.matrix_column.isin(b.columns)].copy()
 g=b.loc[b.index.intersection(wanted),man.matrix_column].T;g.index.name='matrix_column';g=g.reset_index().merge(man,on='matrix_column');save(g,'S02_GSE139533_sample_selection');g=g.groupby(['patient id','stage'])[list(b.index.intersection(wanted))].mean().reset_index().rename(columns={'patient id':'patient','stage':'condition'});g.condition=g.condition.str.lower();g['dataset']='GSE139533';g['unit']='mean log2(CPM + 1) across tissue samples'
 b=pd.read_csv(ROOT/'data/raw/GSE254873/GSE254873_processed.txt.gz',sep='\t',index_col=0);n=b.loc[b.index.intersection(wanted)].T;n['patient']=[re.search(r'\d+',s)[0] for s in n.index];n['condition']=['primary' if 'primary' in s or s.endswith('p') else 'recurrent' for s in n.index];n['dataset']='GSE254873';n['unit']='depositor-normalized log2 NanoString signal'
 b=pd.read_excel(ROOT/'data/raw/GSE254145/GSE254145_SpatialTranscriptomics_normalized_log2_cpm.xlsx',index_col=0);sp=b.loc[b.index.intersection(wanted)].T;sp['patient']=[s.split('_')[1] for s in sp.index];sp['condition']=[s.split('_')[2] for s in sp.index];sp['compartment']=[s.split('_')[3] for s in sp.index];sp['dataset']='GSE254145';sp['unit']='depositor-normalized log2 CPM';save(sp.reset_index(),'S10_spatial_all_compartments')
 tables=[];tests=[];pairedvals=[];coverage=[]
 for a in [e,g,n,sp[sp.compartment=='CD64'].copy()]:
  ds=a.dataset.iloc[0];unit=a.unit.iloc[0];refs=[x for x in REFS if x in a.columns];coverage.append(dict(dataset=ds,reference_genes=';'.join(refs),n_genes=len(refs)))
  a['SIGLEC1_CD68_log_difference']=a.SIGLEC1-a.CD68
  a['Myeloid_reference_score']=a[refs].mean(axis=1)
  a['SIGLEC1_reference_log_difference']=a.SIGLEC1-a.Myeloid_reference_score
  tables.append(a.reset_index(drop=True))
  for mt in ['SIGLEC1','CD68','SIGLEC1_CD68_log_difference','Myeloid_reference_score','SIGLEC1_reference_log_difference']:
   w=a.pivot(index='patient',columns='condition',values=mt).dropna();w=w[['primary','recurrent']].reset_index();w['dataset']=ds;w['metric']=mt;w['unit']=unit;w['delta']=w.recurrent-w.primary;pairedvals.append(w)
   tests.append(dict(dataset=ds,**paired(w,mt,unit)))
 save(pd.concat(tables),'S11_external_specimen_metrics');save(pd.concat(pairedvals),'S11_all_paired_values');res=bh(pd.DataFrame(tests),'dataset');save(res,'S11_external_tests');save(pd.DataFrame(coverage),'S11_reference_coverage')
 # Paired CD64/GFAP compartment contrasts within patient and timepoint, not CD64-negative controls.
 rr=[]
 for mt in ['SIGLEC1','CD68']:
  w=sp.pivot(index=['patient','condition'],columns='compartment',values=mt).dropna(subset=['CD64','GFAP']).reset_index();w['metric']=mt;w['CD64_minus_GFAP']=w.CD64-w.GFAP;rr.append(w)
 save(pd.concat(rr),'S10_spatial_compartment_contrasts')
 return res

def external_single_nucleus():
 source=Path(str(project_path('GSE174554_METADATA').parent))
 labels=pd.read_csv(source/'input/GSE174554_author_celltype_map.csv');clin=pd.read_csv(source/'input/GSE174554_original_clinical.csv')
 clin=clin[clin.Age.ge(18)&clin.IDH.eq('IDH wildtype')&clin.Diagnosis.eq('Glioblastoma')].copy();clin['condition']=clin.Stage.str.lower();clin['patient']=clin['Pair.'].astype(str)
 save(clin,'S12_external_clinical')
 raw=Path(str(project_path('GSE174554_RAW')));audit=[];rows=[];state=[];hashes=[];continuous=[];threshold_rows=[]
 for _,c in clin.iterrows():
  groups=list(labels[labels.UCSF_ID==c.ID].groupby('SampleID'));missing=[];pieces=[];read_errors=[]
  if c.ID=='SF8963':audit.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,status='excluded source matrix QC failure documented in local source audit'));continue
  for stem,a in groups:
   paths=[raw/(stem+s) for s in ['_matrix.mtx.gz','_features.tsv.gz','_barcodes.tsv.gz']]
   if not all(p.exists() for p in paths):missing.append(stem);continue
   try:
    gf=pd.read_csv(paths[1],sep='\t',header=None);genes=pd.Index(gf.iloc[:,1 if gf.shape[1]>1 else 0].astype(str));bars=pd.read_csv(paths[2],header=None)[0].astype(str)
    x=mmread(str(paths[0])).tocsr();assert x.shape==(len(genes),len(bars));keys=pd.Index(stem+'-'+bars);ii=keys.get_indexer(a.CellID);assert np.all(ii>=0)
    x=x[:,ii].T.tocsr();lib=np.asarray(x.sum(axis=1)).ravel();ng=np.asarray((x>0).sum(axis=1)).ravel();mt=np.asarray(x[:,genes.str.startswith('MT-')].sum(axis=1)).ravel()/np.maximum(lib,1)
    keep=(lib>=500)&(ng>=300)&(ng<=8000)&(mt<=.2);is_m=(a.CellType.to_numpy()=='Myeloid')&keep
    pieces.append((x[is_m],genes,int(keep.sum())))
    hashes.extend(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths)
   except Exception as exc:read_errors.append(stem+': '+str(exc))
  if not groups or missing or read_errors:
   audit.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,status='excluded incomplete library or read error',missing=';'.join(missing),error=';'.join(read_errors)));continue
  n=sum(x.shape[0] for x,_,_ in pieces);total=sum(t for _,_,t in pieces)
  if n<20:audit.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,status='excluded fewer than 20 QC myeloid nuclei',myeloid_nuclei=n));continue
  arrays=[]
  for x,genes,_ in pieces:
   y=norm(x);v=pd.DataFrame({'SIGLEC1':gv(y,genes,'SIGLEC1'),'CD68':gv(y,genes,'CD68'),'MG_score':score(y,genes,MG),'MDM_score':score(y,genes,MDM),'UMI':np.asarray(x.sum(axis=1)).ravel(),'genes':np.asarray((x>0).sum(axis=1)).ravel()})
   for name,gs in PROGRAMS.items():v[name]=score(y,genes,gs)
   arrays.append(v)
  a=pd.concat(arrays,ignore_index=True);gap=a.MG_score-a.MDM_score;a['origin_proxy']=np.where(gap>.1,'Microglia-like',np.where(gap<-.1,'Macrophage-like','Ambiguous'))
  for lab,g in [('All myeloid',a)]+list(a.groupby('origin_proxy')):
   if len(g)<20:continue
   rows.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,origin_proxy=lab,n_cells=len(g),SIGLEC1=g.SIGLEC1.mean(),CD68=g.CD68.mean(),Detection=(g.SIGLEC1>0).mean(),Log_expression_difference=g.SIGLEC1.mean()-g.CD68.mean(),myeloid_fraction=len(g)/total))
  for name in PROGRAMS:
   pos=a.SIGLEC1>0
   if pos.sum()>=5 and (~pos).sum()>=5:state.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,program=name,n_positive=int(pos.sum()),n_negative=int((~pos).sum()),delta=a.loc[pos,name].mean()-a.loc[~pos,name].mean(),rho=stats.spearmanr(a.SIGLEC1,a[name]).statistic))
   if a.SIGLEC1.std()>0:
    z=np.column_stack([np.ones(len(a)),np.log1p(a.UMI),np.log1p(a.genes)]);xr=stats.rankdata(a.SIGLEC1).astype(float);yr=stats.rankdata(a[name]).astype(float);xr-=z@np.linalg.lstsq(z,xr,rcond=None)[0];yr-=z@np.linalg.lstsq(z,yr,rcond=None)[0]
    continuous.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,program=name,rho=stats.spearmanr(a.SIGLEC1,a[name]).statistic,depth_adjusted_rank_correlation=np.corrcoef(xr,yr)[0,1]))
   pos=a.SIGLEC1>1
   if pos.sum()>=5 and (~pos).sum()>=5:threshold_rows.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,program=name,n_positive=int(pos.sum()),n_negative=int((~pos).sum()),delta=a.loc[pos,name].mean()-a.loc[~pos,name].mean()))
  audit.append(dict(sample=c.ID,patient=c.patient,condition=c.condition,status='processed',myeloid_nuclei=n,total_QC_nuclei=total,n_libraries=len(groups)));print('EXTERNAL',c.ID,n,flush=True)
 save(pd.DataFrame(audit),'S12_external_eligibility');d=save(pd.DataFrame(rows),'S12_external_patient_metrics');save(pd.DataFrame(state),'S12_external_program_deltas');(TAB/'S12_raw_file_hashes.json').write_text(json.dumps(hashes,indent=2))
 rr=[];pairs=[]
 for label,g in d.groupby('origin_proxy'):
  for mt in ['SIGLEC1','CD68','Detection','Log_expression_difference','myeloid_fraction']:
   w=g.pivot_table(index='patient',columns='condition',values=mt,aggfunc='mean').dropna();
   if len(w)<2 or not {'primary','recurrent'}<=set(w.columns):continue
   w=w.reset_index();w['origin_proxy']=label;w['metric']=mt;w['delta']=w.recurrent-w.primary;pairs.append(w);rr.append(dict(origin_proxy=label,**paired(w,mt,'mean ln(1 + counts per 10,000); fractions for detection/composition')))
 save(pd.concat(pairs),'S12_external_paired_values');res=save(bh(pd.DataFrame(rr),'origin_proxy'),'S12_external_tests')
 pr=[]
 for (condition,name),g in pd.DataFrame(state).groupby(['condition','program']):pr.append(dict(condition=condition,program=name,n_patients=len(g),median_delta=g.delta.median(),increased=int((g.delta>0).sum()),p=wil(g.delta)))
 save(bh(pd.DataFrame(pr),'condition'),'S12_external_program_tests')
 cs=save(pd.DataFrame(continuous),'S12_external_continuous_correlations');tt=[]
 for (c,n),g in cs.groupby(['condition','program']):
  for mt in ['rho','depth_adjusted_rank_correlation']:tt.append(dict(condition=c,program=n,metric=mt,n_patients=len(g),median_correlation=g[mt].median(),p=wil(g[mt])))
 save(bh(pd.DataFrame(tt),['condition','metric']),'S12_external_continuous_tests');ds=save(pd.DataFrame(threshold_rows),'S12_external_threshold_deltas');tt=[]
 for (c,n),g in ds.groupby(['condition','program']):tt.append(dict(condition=c,program=n,n_patients=len(g),median_delta=g.delta.median(),p=wil(g.delta)))
 save(bh(pd.DataFrame(tt),'condition'),'S12_external_threshold_tests')
 print(res.to_string(index=False),flush=True)

if __name__=='__main__':
 import sys
 if '--external-only' in sys.argv:external_single_nucleus()
 else:
  discovery();print(bulk().to_string(index=False));external_single_nucleus()
