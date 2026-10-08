
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread
from statsmodels.stats.multitest import multipletests
REV=Path(str(project_path('ROOT') / 'Revision'))
sys.path.insert(0,str(REV))
from revision_analysis import PROGRAMS, MG, MDM, norm, gv, score
OUT=Path(__file__).parent; TAB=OUT/'analysis_tables';TAB.mkdir(exist_ok=True)
RNG=np.random.default_rng(341008)
def save(d,name):d.to_csv(TAB/(name+'.csv'),index=False,encoding='utf-8-sig');return d
def wil(v):return float(stats.wilcoxon(v,method='auto').pvalue) if np.any(v!=0) else 1.
def bh(d,group):
 d=d.copy();d['q']=d.groupby(group)['p'].transform(lambda v:multipletests(v.fillna(1),method='fdr_bh')[1]);return d
def interactions():
 rows=[];vals=[]
 for fname,metrics,definition in [
 ('S12_external_continuous_correlations',['rho','depth_adjusted_rank_correlation'],'continuous'),
 ('S12_external_program_deltas',['delta'],'detected'),
 ('S12_external_threshold_deltas',['delta'],'threshold_gt1')]:
  a=pd.read_csv(REV/'analysis_tables'/(fname+'.csv'))
  for prog,g in a.groupby('program'):
   for metric in metrics:
    w=g.pivot_table(index='patient',columns='condition',values=metric).dropna();w=w[['primary','recurrent']].reset_index();v=(w.recurrent-w.primary).to_numpy()
    boot=np.median(RNG.choice(v,(10000,len(v))),axis=1)
    rows.append(dict(definition=definition,metric=metric,program=prog,n_pairs=len(v),primary_median=w.primary.median(),recurrent_median=w.recurrent.median(),median_change=np.median(v),ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975),increased=int((v>0).sum()),p=wil(v)))
    w['definition']=definition;w['metric']=metric;w['program']=prog;w['change']=v;vals.append(w)
 save(pd.concat(vals),'S14_paired_program_association_values');res=save(bh(pd.DataFrame(rows),['definition','metric']),'S14_stage_interaction_tests');print(res.to_string(index=False),flush=True)
def extract():
 cache=TAB/'S14_external_specimen_program_means.csv'
 if cache.exists():return pd.read_csv(cache)
 source=Path(str(project_path('GSE174554_METADATA')))
 labels=pd.read_csv(source/'GSE174554_author_celltype_map.csv');clinical=pd.read_csv(source/'GSE174554_original_clinical.csv')
 prev=pd.read_csv(REV/'analysis_tables/S12_external_patient_metrics.csv');prev=prev[prev.origin_proxy.eq('All myeloid')]
 paired=prev.pivot(index='patient',columns='condition',values='SIGLEC1').dropna().index
 eligible=prev[prev.patient.isin(paired)];raw=Path(str(project_path('GSE174554_RAW')));rows=[];coverage=[]
 for _,c in eligible.iterrows():
  arrays=[];total=0
  for stem,a in labels[labels.UCSF_ID.eq(c['sample'])].groupby('SampleID'):
   gf=pd.read_csv(raw/(stem+'_features.tsv.gz'),sep='\t',header=None);genes=pd.Index(gf.iloc[:,1 if gf.shape[1]>1 else 0].astype(str));bars=pd.read_csv(raw/(stem+'_barcodes.tsv.gz'),header=None)[0].astype(str)
   x=mmread(str(raw/(stem+'_matrix.mtx.gz'))).tocsr();keys=pd.Index(stem+'-'+bars);ii=keys.get_indexer(a.CellID);assert np.all(ii>=0)
   x=x[:,ii].T.tocsr();lib=np.asarray(x.sum(axis=1)).ravel();ng=np.asarray((x>0).sum(axis=1)).ravel();mt=np.asarray(x[:,genes.str.startswith('MT-')].sum(axis=1)).ravel()/np.maximum(lib,1)
   keep=(lib>=500)&(ng>=300)&(ng<=8000)&(mt<=.2);total+=int(keep.sum());sel=keep&a.CellType.eq('Myeloid').to_numpy();x=x[sel];y=norm(x)
   v=pd.DataFrame({'SIGLEC1':gv(y,genes,'SIGLEC1'),'CD68':gv(y,genes,'CD68'),'UMI':lib[sel],'genes':ng[sel],'MG':score(y,genes,MG),'MDM':score(y,genes,MDM)})
   for p,gs in PROGRAMS.items():
    v[p]=score(y,genes,gs);coverage.append(dict(sample=c['sample'],library=stem,program=p,present=sum(g in genes for g in gs),expected=len(gs)))
   arrays.append(v)
  a=pd.concat(arrays,ignore_index=True);assert len(a)==int(c.n_cells)
  assert np.isclose(a.SIGLEC1.mean(),c.SIGLEC1,atol=1e-10)
  gap=a.MG-a.MDM;row=dict(sample=c['sample'],patient=c.patient,condition=c.condition,n_cells=len(a),SIGLEC1=a.SIGLEC1.mean(),CD68=a.CD68.mean(),Detection=(a.SIGLEC1>0).mean(),Log_expression_difference=a.SIGLEC1.mean()-a.CD68.mean(),myeloid_fraction=len(a)/total,ln_UMI=np.log1p(a.UMI).mean(),ln_genes=np.log1p(a.genes).mean(),Macrophage_proxy_fraction=(gap<-.1).mean(),Microglia_proxy_fraction=(gap>.1).mean())
  row.update({p:a[p].mean() for p in PROGRAMS});rows.append(row);print('EXTRACT',c['sample'],len(a),flush=True)
 save(pd.DataFrame(coverage),'S14_program_coverage');return save(pd.DataFrame(rows),'S14_external_specimen_program_means')
def correlation(x,y):return float(stats.spearmanr(x,y).statistic)
def coupling(a):
 cols=['SIGLEC1','Detection','Log_expression_difference','myeloid_fraction','ln_UMI','ln_genes','Macrophage_proxy_fraction','Microglia_proxy_fraction']+list(PROGRAMS)
 delta=pd.DataFrame({col:a.pivot(index='patient',columns='condition',values=col).apply(lambda w:w.recurrent-w.primary,axis=1) for col in cols});assert len(delta)==17
 save(delta.reset_index(),'S14_patient_longitudinal_changes');rows=[];leave=[]
 ranks={col:stats.rankdata(delta[col]) for col in cols}
 n=len(delta);perm=np.array([RNG.permutation(n) for _ in range(19999)]);boot=RNG.integers(0,n,(10000,n))
 for metric in ['SIGLEC1','Detection','Log_expression_difference']:
  for prog in PROGRAMS:
   x=delta[metric].to_numpy();y=delta[prog].to_numpy();rho=correlation(x,y)
   xr=ranks[metric]-np.mean(ranks[metric]);yr=ranks[prog]-np.mean(ranks[prog]);pr=(xr[perm]@yr)/(np.linalg.norm(xr)*np.linalg.norm(yr));p=(1+np.sum(np.abs(pr)>=abs(rho)-1e-12))/20000
   ci=np.array([correlation(x[b],y[b]) for b in boot]);ci=ci[np.isfinite(ci)]
   row=dict(metric=metric,program=prog,n_pairs=n,rho=rho,p=float(p),spearman_asymptotic_p=float(stats.spearmanr(x,y).pvalue),ci_low=np.quantile(ci,.025),ci_high=np.quantile(ci,.975))
   for name,covs in [('depth',['ln_UMI','ln_genes']),('composition',['myeloid_fraction','Macrophage_proxy_fraction']),('depth_composition',['ln_UMI','ln_genes','myeloid_fraction','Macrophage_proxy_fraction'])]:
    z=np.column_stack([np.ones(n)]+[ranks[k] for k in covs]);rx=xr-z@np.linalg.lstsq(z,xr,rcond=None)[0];ry=yr-z@np.linalg.lstsq(z,yr,rcond=None)[0];r=np.corrcoef(rx,ry)[0,1];t=r*np.sqrt((n-len(covs)-2)/max(1e-15,1-r*r));row[name+'_rho']=r;row[name+'_p']=float(2*stats.t.sf(abs(t),n-len(covs)-2))
   rows.append(row)
   for i,pat in enumerate(delta.index):
    ix=np.arange(n)!=i;leave.append(dict(metric=metric,program=prog,omitted_patient=pat,rho=correlation(x[ix],y[ix])))
 res=bh(pd.DataFrame(rows),'metric')
 for name in ['depth','composition','depth_composition']:res[name+'_q']=res.groupby('metric')[name+'_p'].transform(lambda v:multipletests(v,method='fdr_bh')[1])
 save(res,'S14_longitudinal_coupling_tests');save(pd.DataFrame(leave),'S14_leave_one_patient_out');print(res.to_string(index=False),flush=True)
if __name__=='__main__':
 interactions()
 if '--interaction-only' not in sys.argv:coupling(extract())
