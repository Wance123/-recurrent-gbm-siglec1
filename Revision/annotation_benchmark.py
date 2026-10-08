"""Exact specimen/barcode benchmark against the original GSE200984 labels."""

# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import json, re, hashlib
import h5py, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from revision_analysis import ROOT, frame, vector, matrix, norm, gv, score, PROGRAMS

OUT=ROOT/'Revision/reviewer_requirements_refined_20261008/annotation_results'
OUT.mkdir(parents=True,exist_ok=True)
LABELS=Path(str(project_path('AUTHOR_WORKSPACE') / 'metadata/GSE200984_author_all_nuclei_metadata.csv'))
RNG=np.random.default_rng(341008)
def save(d,n):d.to_csv(OUT/(n+'.csv'),index=False,encoding='utf-8-sig');return d
def canonical_barcode(s):
    # Preserve the full 16-base barcode and GEM identifier. Only explicit specimen
    # prefixes and author-added specimen/library suffixes are discarded.
    m=re.search(r'([ACGTN]{16})[-.]([0-9]+)',str(s))
    return m.group(1)+'-'+m.group(2) if m else str(s)
author=pd.read_csv(LABELS)
clinical=pd.read_csv(ROOT/'data/metadata/GSE228500_patient_metadata.tsv',sep='\t')
sample_map=dict(zip(clinical['Case'], 'TB'+clinical['TB #'].astype(str)))
def map_sample(s):
    if str(s).startswith('TB'):return s
    return sample_map.get(re.sub(r'_[cm]$','',str(s)))
author['sample']=author['orig.ident'].map(map_sample)
author['barcode']=author.author_barcode.str.replace(r'^[^_]+(?:_[cm])?_','',regex=True).map(canonical_barcode)
author['key']=author['sample'].fillna('UNMAPPED')+'_'+author.barcode
author['ambiguous_key']=author['sample'].notna() & author.key.duplicated(keep=False)
au=author[~author.ambiguous_key & author['sample'].notna()].copy().set_index('key')
with h5py.File(ROOT/'data/processed/GSE200984_strict_basic_atlas_cnv.h5ad') as f:
    obs=frame(f['obs']);counts=matrix(f['layers/counts']);genes=pd.Index(vector(f['var/gene']))
with h5py.File(ROOT/'data/processed/GSE200984_myeloid_reclustered.h5ad') as f:
    my=frame(f['obs'])
obs['key']=obs['sample'].astype(str)+'_'+obs.barcode.map(canonical_barcode)
assert not obs.key.duplicated().any(), 'Reconstructed key is not unique'
obs['reconstructed_cluster']=obs.leiden_0_6.astype(str)
mapping={str(k):('Myeloid' if k==9 else 'T_NK' if k==12 else 'Oligodendrocyte' if k in [10,11,13] else 'Neuron' if k==4 else 'Other_tumor_like') for k in range(18)}
obs['reconstructed_lineage']=obs.reconstructed_cluster.map(mapping)
obs['author_cluster']=obs.key.map(au.Cluster)
obs['author_condition']=obs.key.map(au.Condition)
obs['matched']=obs.author_cluster.notna()
assert (obs.loc[obs.matched,'author_condition'].map({'P':'Primary','R':'Recurrent'}).to_numpy()==obs.loc[obs.matched,'condition'].to_numpy()).all()
mylabels={'Myel1','moTAM','mgTAM','prTAM'}
obs['author_myeloid']=obs.author_cluster.isin(mylabels)
clean=my[my.myeloid_leiden.astype(str)!='3'].index
obs['reconstructed_clean_myeloid']=obs.index.isin(clean)
coverage=[]
for sample,g in obs.groupby('sample',observed=True):
    a=author[author['sample'].eq(sample)]
    coverage.append(dict(sample=sample,condition=str(g.condition.iloc[0]),reconstructed_nuclei=len(g),matched_nuclei=int(g.matched.sum()),match_fraction=float(g.matched.mean()),author_nuclei=len(a),author_ambiguous_key_nuclei=int(a.ambiguous_key.sum()),reconstructed_clean_myeloid=int(g.reconstructed_clean_myeloid.sum()),matched_clean_myeloid=int((g.reconstructed_clean_myeloid & g.matched).sum())))
save(pd.DataFrame(coverage),'S16_match_coverage')
cross=obs[obs.matched].groupby(['reconstructed_lineage','author_cluster'],observed=True).size().rename('nuclei').reset_index()
save(cross,'S16_label_crosswalk')
save(obs[['sample','condition','key','reconstructed_cluster','reconstructed_lineage','reconstructed_clean_myeloid','author_cluster','author_myeloid','matched']].reset_index(),'S16_cell_mapping')
binary=[]
for mode,col in [('Broad reconstructed myeloid','reconstructed_lineage'),('Clean reconstructed myeloid','reconstructed_clean_myeloid')]:
    g=obs[obs.matched];pred=g[col].eq('Myeloid') if col=='reconstructed_lineage' else g[col];actual=g.author_myeloid
    tp=int((pred & actual).sum());fp=int((pred & ~actual).sum());fn=int((~pred & actual).sum());tn=int((~pred & ~actual).sum())
    binary.append(dict(definition=mode,matched_nuclei=len(g),true_positive=tp,false_positive=fp,false_negative=fn,true_negative=tn,precision=tp/(tp+fp),recall=tp/(tp+fn),f1=2*tp/(2*tp+fp+fn)))
save(pd.DataFrame(binary),'S16_myeloid_concordance')
x=norm(counts)
obs['SIGLEC1']=gv(x,genes,'SIGLEC1');obs['CD68']=gv(x,genes,'CD68');obs['Detection']=obs.SIGLEC1.gt(0)
obs['Log_expression_difference']=obs.SIGLEC1-obs.CD68
for p,gs in PROGRAMS.items():obs[p]=score(x,genes,gs)
metrics=[];tests=[];corr=[]
sets={'Reconstructed clean on matched cells':obs.matched & obs.reconstructed_clean_myeloid,
      'Original author all myeloid':obs.matched & obs.author_myeloid,
      'Original author TAM labels':obs.matched & obs.author_cluster.isin(['moTAM','mgTAM','prTAM']),
      'Original author moTAM':obs.matched & obs.author_cluster.eq('moTAM'),
      'Original author mgTAM':obs.matched & obs.author_cluster.eq('mgTAM')}
for label,mask in sets.items():
    d=obs[mask];a=d.groupby(['sample','condition'],observed=True).agg(SIGLEC1=('SIGLEC1','mean'),CD68=('CD68','mean'),Detection=('Detection','mean'),Log_expression_difference=('Log_expression_difference','mean'),n_cells=('sample','size')).reset_index();a=a[a.n_cells>=20].copy();a['definition']=label;metrics.append(a)
    for metric in ['SIGLEC1','Detection','Log_expression_difference']:
        p=a.loc[a.condition.eq('Primary'),metric].to_numpy();r=a.loc[a.condition.eq('Recurrent'),metric].to_numpy()
        if len(p)<2 or len(r)<2:
            tests.append(dict(definition=label,metric=metric,n_primary=len(p),n_recurrent=len(r),difference=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,status='insufficient specimens per stage; no inferential test'))
            continue
        bs=np.median(RNG.choice(r,(10000,len(r))),axis=1)-np.median(RNG.choice(p,(10000,len(p))),axis=1)
        tests.append(dict(definition=label,metric=metric,n_primary=len(p),n_recurrent=len(r),difference=float(np.median(r)-np.median(p)),ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),p=float(stats.mannwhitneyu(r,p,alternative='two-sided',method='auto').pvalue)))
    if label not in ['Original author all myeloid','Reconstructed clean on matched cells']:continue
    for (sample,condition),g in d.groupby(['sample','condition'],observed=True):
        if len(g)<20 or not g.SIGLEC1.std():continue
        z=np.column_stack([np.ones(len(g)),np.log1p(g.total_counts.astype(float)),np.log1p(g.n_genes_by_counts.astype(float))])
        for pr in PROGRAMS:
            xr=stats.rankdata(g.SIGLEC1).astype(float);yr=stats.rankdata(g[pr]).astype(float)
            xr-=z@np.linalg.lstsq(z,xr,rcond=None)[0];yr-=z@np.linalg.lstsq(z,yr,rcond=None)[0]
            corr.append(dict(definition=label,sample=sample,condition=str(condition),program=pr,rho=float(np.corrcoef(xr,yr)[0,1]),n_cells=len(g)))
save(pd.concat(metrics),'S16_patient_expression_values')
t=pd.DataFrame(tests);t['q']=t.groupby('definition',observed=True).p.transform(lambda p:multipletests(p.fillna(1),method='fdr_bh')[1]);t.loc[t.p.isna(),'q']=np.nan;save(t,'S16_original_label_expression_tests')
c=save(pd.DataFrame(corr),'S16_patient_program_correlations');r=[]
for (label,condition,pr),g in c.groupby(['definition','condition','program'],observed=True):
    v=g.rho.to_numpy();r.append(dict(definition=label,condition=condition,program=pr,n_specimens=len(v),median_rho=np.median(v),p=float(stats.wilcoxon(v,method='auto',alternative='two-sided').pvalue) if np.any(v) else 1.))
rt=pd.DataFrame(r);rt['q']=rt.groupby(['definition','condition'],observed=True).p.transform(lambda p:multipletests(p,method='fdr_bh')[1]);save(rt,'S16_original_label_program_tests')
summary={'author_workspace_url':'https://drive.google.com/file/d/1LmQX3vD6voSvO3dTbG9b2kY-kW3U-dhP/view','author_workspace_sha256':'3763b43b985ebacb8375a8461de87a213a08761a29cc38dfe28f19e87614c0cc','author_metadata_nuclei':len(author),'reconstructed_nuclei':len(obs),'matched_nuclei':int(obs.matched.sum()),'match_fraction':float(obs.matched.mean()),'ambiguous_author_keys_excluded':int(author.ambiguous_key.sum()),'myeloid_definition':sorted(mylabels),'matching':'Published Case to TB specimen map; exact 16-base barcode and GEM identifier after removal of explicit specimen/library prefixes and suffixes; no fuzzy barcode matching; duplicate author specimen/barcode keys excluded','concordance':binary}
(OUT/'benchmark_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2));print(t.to_string(index=False));print(rt[rt.program.eq('Interferon response')].to_string(index=False))
