
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import json, hashlib, zipfile
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).parent
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
GENES=['SIGLEC1','CD68','AIF1','LST1','TYROBP','CSF1R']
OUTCOMES=['SIGLEC1','CD68','SIGLEC1_minus_CD68','SIGLEC1_minus_four_refs','SIGLEC1_minus_two_refs']
RNG=np.random.default_rng(341008)
ledger=[]; tests=[]

def register(p):
    p=Path(p);ledger.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    return p

def save(d,name):d.to_csv(OUT/(name+'.csv'),index=False,encoding='utf-8-sig')
def zipread(p):
    with zipfile.ZipFile(register(p)) as z:
        fs=[n for n in z.namelist() if n.lower().endswith('.txt') and not n.startswith('__MACOSX/')];assert len(fs)==1
        return pd.read_csv(z.open(fs[0]),sep='\t')
def outcomes(d):
    d=d.copy();d['SIGLEC1_minus_CD68']=d.SIGLEC1-d.CD68
    d['SIGLEC1_minus_four_refs']=d.SIGLEC1-d[['AIF1','LST1','TYROBP','CSF1R']].mean(axis=1,skipna=False)
    d['SIGLEC1_minus_two_refs']=d.SIGLEC1-d[['TYROBP','CSF1R']].mean(axis=1,skipna=False)
    return d
def paired(d,cohort,compartment,unit):
    for metric in OUTCOMES:
        w=d.pivot(index='patient',columns='condition',values=metric).dropna(subset=['primary','recurrent'])
        if not {'primary','recurrent'}<=set(w.columns):
            tests.append(dict(cohort=cohort,compartment=compartment,metric=metric,design='paired',unit=unit,n_primary=0,n_recurrent=0,p=np.nan,status='unavailable gene or stage'));continue
        v=(w.recurrent-w.primary).to_numpy()
        if len(v)<2:
            tests.append(dict(cohort=cohort,compartment=compartment,metric=metric,design='paired',unit=unit,n_primary=len(v),n_recurrent=len(v),p=np.nan,status='insufficient complete pairs'));continue
        p=float(stats.wilcoxon(v,alternative='two-sided',method='auto').pvalue) if np.any(v) else 1.
        bs=np.median(RNG.choice(v,(10000,len(v))),axis=1)
        tests.append(dict(cohort=cohort,compartment=compartment,metric=metric,design='paired',unit=unit,n_primary=len(v),n_recurrent=len(v),primary_median=w.primary.median(),recurrent_median=w.recurrent.median(),effect=np.median(v),ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),increased=int((v>0).sum()),p=p))
        z=w.reset_index();z['metric']=metric;z['change']=v;save(z,f'{cohort}_{compartment}_{metric}_paired_values')

def cgga():
    for suffix,base in [('693',Path(str(project_path('CGGA693')))),('325',Path(str(project_path('CGGA325'))))]:
        cohort='CGGA'+suffix;cl=zipread(base/f'CGGA.mRNAseq_{suffix}_clinical.20200506.txt.zip')
        ex=zipread(base/f'CGGA.mRNAseq_{suffix}.RSEM-genes.20200506.txt.zip')
        assert 'Gene_Name' in ex
        ex=ex[ex.Gene_Name.isin(GENES)].set_index('Gene_Name');assert not ex.index.duplicated().any() and set(GENES)<=set(ex.index)
        assert np.isfinite(ex.to_numpy(dtype=float)).all() and (ex.to_numpy(dtype=float)>=0).all()
        elig=(cl.Histology.isin(['GBM','rGBM']) & cl.Grade.eq('WHO IV') & cl.IDH_mutation_status.eq('Wildtype') & pd.to_numeric(cl.Age,errors='coerce').ge(18) & cl.PRS_type.isin(['Primary','Recurrent']) & cl.CGGA_ID.isin(ex.columns))
        a=cl.copy();a['eligible']=elig;save(a,cohort+'_eligibility')
        d=cl[elig].set_index('CGGA_ID');assert not d.index.duplicated().any()
        d=d.join(np.log2(ex.loc[GENES,d.index].T.astype(float)+1));d['condition']=d.PRS_type.str.lower();d=outcomes(d)
        save(d.reset_index(),cohort+'_specimen_values')
        for metric in OUTCOMES:
            x=d.loc[d.condition.eq('primary'),metric].to_numpy();y=d.loc[d.condition.eq('recurrent'),metric].to_numpy()
            assert len(x)>=2 and len(y)>=2, 'insufficient comparison group; no result may be assigned'
            mw=stats.mannwhitneyu(y,x,alternative='two-sided',method='auto')
            bs=np.median(RNG.choice(y,(10000,len(y))),axis=1)-np.median(RNG.choice(x,(10000,len(x))),axis=1)
            tests.append(dict(cohort=cohort,compartment='tissue',metric=metric,design='unpaired',unit='log2(FPKM + 1)',n_primary=len(x),n_recurrent=len(y),primary_median=np.median(x),recurrent_median=np.median(y),effect=np.median(y)-np.median(x),ci_low=np.quantile(bs,.025),ci_high=np.quantile(bs,.975),increased=np.nan,cliffs_delta=2*mw.statistic/(len(x)*len(y))-1,p=float(mw.pvalue)))
        print(cohort,'eligible',d.condition.value_counts().to_dict(),flush=True)

def care():
    p=ROOT/'CARE_specimen_metrics.csv'
    if not p.exists():return
    a=pd.read_csv(register(p));cl=pd.read_csv(register(ROOT/'CARE_clinical_eligibility.csv'))
    a=a.merge(cl[['Sample ID','Patient ID','Time point','condition','eligible']],left_on='sample',right_on='Sample ID',validate='many_to_one')
    a=a[a.eligible.eq(True)].copy();a['patient']=a['Patient ID'];a['time_order']=a['Time point'].str.extract(r'(\d+)')[0].astype(int)
    # Choose the first clinically recorded recurrent sample BEFORE myeloid-QC exclusion.
    selected=cl[cl.eligible.eq(True) & cl.condition.notna()].copy();selected['time_order']=selected['Time point'].str.extract(r'(\d+)')[0].astype(int)
    selected=selected.sort_values('time_order').drop_duplicates(['Patient ID','condition'],keep='first')
    a=a[a['Sample ID'].isin(selected['Sample ID'])].copy();save(a,'CARE_all_selected_specimen_metrics')
    prog=['Interferon_response','Phagolysosome','Lipid_and_complement','Antigen_presentation','Immunoregulatory','Immediate_inflammatory']
    association=[];stage=[]
    for group in ['TAM','All_nuclei']:
        d=a[a.group.eq(group) & a.n_cells.ge(20)].copy()
        continuous=outcomes(d);save(continuous,'CARE_'+group+'_mean_cell_values')
        paired(continuous,'CARE',group+'_mean_cell','mean ln(1 + counts normalized to 10000)')
        pb=d.copy()
        for g in GENES:pb[g]=np.log2(pb['counts_'+g]/pb.library_counts*1e6+1)
        pb=outcomes(pb);save(pb,'CARE_'+group+'_pseudobulk_values')
        paired(pb,'CARE',group+'_pseudobulk','log2(CPM + 1)')
        if group=='TAM':
            for condition in ['primary','recurrent']:
                for pr in prog:
                    v=d.loc[d.condition.eq(condition),'depth_rho_'+pr].dropna().to_numpy()
                    p=float(stats.wilcoxon(v).pvalue) if len(v)>1 and np.any(v) else 1.
                    association.append(dict(condition=condition,program=pr,n_specimens=len(v),median_rho=np.median(v),p=p))
            for pr in prog:
                w=d.pivot(index='patient',columns='condition',values='depth_rho_'+pr).dropna(subset=['primary','recurrent']);v=(w.recurrent-w.primary).to_numpy()
                p=float(stats.wilcoxon(v).pvalue) if len(v)>1 and np.any(v) else 1.
                stage.append(dict(program=pr,n_pairs=len(v),median_change=np.median(v),p=p))
    b=pd.DataFrame(association);b['q']=b.groupby('condition').p.transform(lambda x:multipletests(x,method='fdr_bh')[1]);save(b,'CARE_depth_adjusted_program_associations')
    b=pd.DataFrame(stage);b['q']=multipletests(b.p,method='fdr_bh')[1];save(b,'CARE_paired_stage_association_changes')

if __name__=='__main__':
    cgga();care()
    d=pd.DataFrame(tests);d['q_within_five']=d.groupby(['cohort','compartment']).p.transform(lambda v:multipletests(v,method='fdr_bh')[1])
    primary=d[d.metric.eq('SIGLEC1') & (d.compartment.isin(['tissue','TAM_mean_cell']))].copy()
    primary['q_across_new_cohorts']=multipletests(primary.p,method='fdr_bh')[1]
    save(d,'all_fixed_outcome_tests');save(primary,'primary_SIGLEC1_tests');save(pd.DataFrame(ledger),'input_sha256')
    (ROOT/'analysis_summary.json').write_text(primary.to_json(orient='records',indent=2),encoding='utf-8')
    print(primary.to_string(index=False),flush=True)
