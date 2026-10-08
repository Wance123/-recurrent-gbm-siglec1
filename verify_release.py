"""Recompute supplied derived rank-test P values without raw matrices."""
from pathlib import Path
import json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parent
R=ROOT/'Revision/external_validation_expansion_20261008/results'
tests=pd.read_csv(R/'combined_all_fixed_outcome_tests.csv')
checked=0
for _,s in tests[tests.p.notna()].iterrows():
    if s.design=='paired':
        d=pd.read_csv(R/f'{s.cohort}_{s.compartment}_{s.metric}_paired_values.csv')
        v=(d.recurrent-d.primary).to_numpy()
        p=float(stats.wilcoxon(v,method='auto',alternative='two-sided').pvalue) if np.any(v) else 1.
    else:
        d=pd.read_csv(R/f'{s.cohort}_specimen_values.csv')
        x=d.loc[d.condition.eq('primary'),s.metric].to_numpy();y=d.loc[d.condition.eq('recurrent'),s.metric].to_numpy()
        p=float(stats.mannwhitneyu(y,x,method='auto',alternative='two-sided').pvalue)
    assert np.isclose(p,s.p,atol=1e-12,rtol=1e-10),(s.cohort,s.metric,p,s.p)
    checked+=1
for _,g in tests.groupby(['cohort','compartment']):
    q=multipletests(g.p.fillna(1),method='fdr_bh')[1]
    # Audit entries unavailable from source mapping remain missing in presentation.
    valid=g.p.notna().to_numpy()
    assert np.allclose(q[valid],g.loc[g.p.notna(),'q_within_five'].to_numpy(),rtol=1e-10,atol=1e-12)
p=pd.read_csv(R/'combined_primary_SIGLEC1_validation.csv')
assert np.allclose(multipletests(p.p,method='fdr_bh')[1],p.q_across_six_cohorts,rtol=1e-10,atol=1e-12)
A=ROOT/'Revision/reviewer_requirements_refined_20261008/annotation_results'
d=pd.read_csv(A/'S16_myeloid_concordance.csv')
for _,r in d.iterrows():
    tp=r.true_positive;fp=r.false_positive;fn=r.false_negative
    assert np.isclose(r.precision,tp/(tp+fp)) and np.isclose(r.recall,tp/(tp+fn)) and np.isclose(r.f1,2*tp/(2*tp+fp+fn))
a=pd.read_csv(A/'S16_patient_expression_values.csv');t=pd.read_csv(A/'S16_original_label_expression_tests.csv')
for _,r in t[t.p.notna()].iterrows():
    g=a[a.definition.eq(r.definition)];x=g.loc[g.condition.eq('Primary'),r.metric];y=g.loc[g.condition.eq('Recurrent'),r.metric]
    assert np.isclose(stats.mannwhitneyu(y,x,alternative='two-sided',method='auto').pvalue,r.p,atol=1e-12)
print(json.dumps({'expanded_rank_tests_recomputed':checked,'expanded_BH_families_checked':True,'six_cohort_BH_checked':True,'concordance_arithmetic_checked':True,'original_label_expression_tests_checked':int(t.p.notna().sum()),'raw_preprocessing_rerun':False},indent=2))
