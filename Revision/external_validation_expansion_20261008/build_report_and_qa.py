
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).parent;R=ROOT/'results'
alltests=pd.concat([pd.read_csv(R/'all_fixed_outcome_tests.csv'),pd.read_csv(R/'array_fixed_outcome_tests.csv')],ignore_index=True)
primary=alltests[alltests.metric.eq('SIGLEC1') & alltests.compartment.isin(['tissue','TAM_mean_cell'])].copy()
assert len(primary)==6 and primary.cohort.nunique()==6 and primary.p.notna().all()
primary['q_across_six_cohorts']=multipletests(primary.p,method='fdr_bh')[1]
urls={
 'CGGA693':'https://www.cgga.org.cn/download.jsp',
 'CGGA325':'https://www.cgga.org.cn/download.jsp',
 'CARE':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546',
 'GSE212067':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212067',
 'GSE4271':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4271',
 'GSE7696':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7696'}
primary['source_url']=primary.cohort.map(urls)
primary.to_csv(R/'combined_primary_SIGLEC1_validation.csv',index=False,encoding='utf-8-sig')
alltests.to_csv(R/'combined_all_fixed_outcome_tests.csv',index=False,encoding='utf-8-sig')

checks=[]
for _,s in alltests[alltests.design.eq('paired') & alltests.p.notna()].iterrows():
    p=R/f'{s.cohort}_{s.compartment}_{s.metric}_paired_values.csv'
    d=pd.read_csv(p);v=(d.recurrent-d.primary).to_numpy()
    pvalue=stats.wilcoxon(v,alternative='two-sided',method='auto').pvalue if np.any(v) else 1.
    assert len(v)==int(s.n_primary)==int(s.n_recurrent)
    assert np.isclose(pvalue,s.p,rtol=1e-10,atol=1e-12)
    assert np.isclose(np.median(v),s.effect,atol=1e-10)
    assert int((v>0).sum())==int(s.increased)
    checks.append(dict(cohort=s.cohort,compartment=s.compartment,metric=s.metric,status='paired statistic independently recomputed'))
for _,s in alltests[alltests.design.eq('unpaired')].iterrows():
    d=pd.read_csv(R/f'{s.cohort}_specimen_values.csv')
    x=d.loc[d.condition.eq('primary'),s.metric];y=d.loc[d.condition.eq('recurrent'),s.metric]
    assert len(x)==s.n_primary and len(y)==s.n_recurrent
    assert np.isclose(stats.mannwhitneyu(y,x,alternative='two-sided').pvalue,s.p,rtol=1e-10)
    checks.append(dict(cohort=s.cohort,metric=s.metric,status='unpaired statistic independently recomputed'))
for _,g in alltests.groupby(['cohort','compartment']):
    expected=multipletests(g.p.fillna(1),method='fdr_bh')[1]
    ix=g.p.notna();assert np.allclose(g.loc[ix,'q_within_five'],expected[ix.to_numpy()],atol=1e-12)

care=pd.read_csv(ROOT/'CARE_specimen_metrics.csv');cl=pd.read_csv(ROOT/'CARE_clinical_eligibility.csv')
assert care['sample'].nunique()==121 and not care[['sample','group']].duplicated().any()
assert cl.eligible.all() and cl.age_initial.ge(18).all() and cl['IDH mutation status'].eq('IDHwt').all()
coverage=care.filter(like='coverage_');assert coverage.ge(.8).all().all()
a=pd.read_csv(R/'CGGA693_specimen_values.csv');b=pd.read_csv(R/'CGGA325_specimen_values.csv')
overlap=sorted(set(a.CGGA_ID)&set(b.CGGA_ID))
assert not overlap
program=pd.read_csv(R/'CARE_depth_adjusted_program_associations.csv')
stage=pd.read_csv(R/'CARE_paired_stage_association_changes.csv')
assert len(program)==12 and len(stage)==6
for _,g in program.groupby('condition'):assert np.allclose(g.q,multipletests(g.p,method='fdr_bh')[1])
assert np.allclose(stage.q,multipletests(stage.p,method='fdr_bh')[1])

frozen=ROOT.parent/'recurrence_marker_reassessment'
manifest=json.loads((frozen/'submission_package_SHA256_manifest.json').read_text())
unchanged=[]
for item in manifest:
    # Accept the manifest's verified relative-file key without guessing an upload state.
    name=item.get('file',item.get('path'))
    p=frozen/'submission_package'/name
    assert p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
    unchanged.append(name)
qa=dict(numerical_checks_passed=True,completed_new_cohorts=6,fixed_outcome_rows=len(alltests),unavailable_outcomes=int(alltests.p.isna().sum()),recomputed_test_rows=len(checks),CARE_samples=121,CARE_patients=cl['Patient ID'].nunique(),CARE_QC_TAM_nuclei=int(care.loc[care.group.eq('TAM'),'n_cells'].sum()),CARE_primary_paired_myeloid_patients=48,CGGA_shared_sample_identifiers=overlap,cross_database_person_identity_deduplication='not established',candidate_revision_package_unchanged=True,unchanged_package_files=len(unchanged),external_submission=False)
(ROOT/'FINAL_ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')

lines=['SIGLEC1 additional external validation report','Date: 2026-10-08','',
'Conclusion: Computational validation in six additional cohorts/datasets is complete. Inclusion was not determined by positive results.',
'Additional support is concentrated in the tissue-level recurrence association in CGGA693 and the interferon-related myeloid-state association in CARE.',
'The additional data do not establish SIGLEC1 as a recurrence-specific marker. Negative results from the other cohorts are retained.','',
'1. Primary SIGLEC1 tests (two-sided; within-cohort BH across five fixed outcomes; separate BH across the six cohort-level SIGLEC1 tests)']
for _,s in primary.iterrows():
    size=f'{int(s.n_primary)} pairs' if s.design=='paired' else f'primary {int(s.n_primary)}, recurrent {int(s.n_recurrent)}'
    lines.append(f"{s.cohort}: {size}; recurrent-minus-primary effect={s.effect:.6g} ({s.unit}); P={s.p:.6g}; within-cohort q={s.q_within_five:.6g}; six-cohort q={s.q_across_six_cohorts:.6g}.")
lines.extend(['','2. Supported and unsupported interpretations',
'CGGA693: SIGLEC1 and CD68 both increase; the SIGLEC1-CD68 difference does not support a relative increase. P values for the four-reference/two-reference differences are approximately 0.821/0.800. This supports a tissue-level association only; changes in myeloid composition cannot be excluded.',
'CGGA325: The primary SIGLEC1 test does not replicate. Both CGGA comparisons are unpaired and cannot establish within-patient recurrence induction. Included sample IDs do not overlap, but person-level cross-cohort deduplication has not been established.',
'CARE: 121 specimens from 59 adult IDH-wild-type patients were reaggregated using the same QC. The primary myeloid test includes 48 pairs, with increases in 29, but P=0.2148 and q=0.3580 are nonsignificant. Sensitivity analyses of all-nucleus mean expression and myeloid/all-nucleus pseudobulk also do not support an adjusted overall increase.'])
for _,s in program[program.program.eq('Interferon_response')].iterrows():
    lines.append(f"CARE {s.condition}: {int(s.n_specimens)} evaluable myeloid specimens; depth-adjusted median SIGLEC1-interferon rho={s.median_rho:.6g}; within-six-program q={s.q:.6g}.")
s=stage[stage.program.eq('Interferon_response')].iloc[0]
lines.append(f"CARE paired change in interferon association: {int(s.n_pairs)} pairs; change={s.median_change:.6g}; P={s.p:.6g}; q={s.q:.6g}. Association at both stages does not establish recurrence specificity.")
lines.extend(['GSE212067: Seven clearly identified adult IDH-wild-type primary/recurrent GBM pairs do not support increased SIGLEC1.',
'GSE4271: Sixteen adult pairs with historical WHO grade IV tissue at both time points were selected using the original authors matching fields. IDH is unavailable; this is not modern IDH-wild-type GBM validation.',
'GSE7696: Pairing uses the original authors patient field, not assumptions from sample names; three evaluable pairs remain. The sample is very small and P=0.25 is nonsignificant. Bootstrap intervals should not be interpreted as reliable significance evidence.',
'Under the official unique-gene probe rule, CD68 mapping is unavailable in GSE4271/GSE7696. CD68 and SIGLEC1-CD68 outcomes remain missing; no zero imputation or selective probe choice is used.','',
'3. Resources not used as additional independent validation',
'GSE62153: Previously analyzed in the existing project; not counted as a new data source.',
'GSE228497: Shares the GSE228500 source with discovery and cannot count as independent replication.',
'GSE79671: Primarily a pre/post bevacizumab-combination treatment comparison in recurrent GBM, not an unambiguous primary/recurrent validation.',
'GSE42670/GSE42669: The former is a mixed superseries; the latter contains 58 human GBM expression profiles. Direct recurrence labels and patient pairing still require verification; no SIGLEC1 result is reported.',
'GLASS: Expression covers only some specimens. Access terms, source overlap, and IDH/time-point alignment require further review; numerical validation was not completed in this analysis.','',
'4. Data and methodological limitations',
'All additional tests are post-revision exploratory validation, not preregistered confirmatory research. BH correction does not retrospectively correct the earlier candidate-gene selection.',
'Effects measured in different expression units are not pooled directly. Pseudobulk log2(CPM+1), mean cellular natural-log expression, FPKM, and normalized microarray signals are reported separately.',
'New cohorts are selected by design, coverage, and pairing, not positive results. Relaxed thresholds, patient exclusions, and one-sided tests were not used to manufacture conclusions.',
'Person-level overlap across data sources has not been fully excluded. No direct overlap was found in CARE/GSE174554 specimen-source and tumor IDs, but this does not replace person-level deduplication.','',
'5. Suggested manuscript interpretation',
'CGGA693 can add tissue-level association evidence, while nonreplication in CGGA325 and the three paired tissue cohorts is reported in full.',
'CARE strengthens cross-cohort evidence that SIGLEC1 corresponds to an interferon-related myeloid state. The overall change in 48 pairs and stage-associated change in 47 pairs remain explicitly negative.',
'Retain the exploratory interpretation of a candidate, context-dependent marker of recurrence-associated myeloid remodeling; do not upgrade this to a validated recurrence-specific marker.',
'Stop open-ended searches for positive databases. Next, align the evidence level in the manuscript and response and include all additional tests in the supplementary materials.','',
'6. Quality checks and files',
f'There are {len(alltests)} primary and sensitivity outcome rows, including {alltests.p.isna().sum()} outcomes not estimable because genes are missing. The remaining tests and BH values were independently recomputed.',
'The original revision submission-package SHA256 hashes are unchanged; no manuscript was modified or uploaded during this analysis.',
'Statistics and patient-level data are in results; source metadata, processed matrices, and probe mappings are in input.',
'Only public metadata, processed matrices, and required platform annotations were downloaded. FASTQ and CEL files were not downloaded; large CARE and CGGA count inputs reused existing caches.','',
'Public records:']+[f'{k}: {v}' for k,v in urls.items()])
(ROOT/'Additional_external_validation_report.txt').write_text('\n'.join(lines),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False,indent=2));print(primary[['cohort','n_primary','n_recurrent','p','q_within_five','q_across_six_cohorts']].to_string(index=False))
