
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

lines=['SIGLEC1 新增外部验证报告','日期：2026-10-08','',
'结论：新增6个队列/数据集的计算验证已完成，不以是否阳性决定纳入。',
'新增支持集中在CGGA693组织层面的复发关联，以及CARE的干扰素相关髓系状态关联。',
'新增数据不确证SIGLEC1为复发特异性标志物。其他队列的阴性结果全部保留。','',
'一、SIGLEC1主要检验（双侧；组内BH按固定5个指标；另对6队列的SIGLEC1检验统一BH）']
for _,s in primary.iterrows():
    size=f'{int(s.n_primary)}对' if s.design=='paired' else f'原发{int(s.n_primary)}、复发{int(s.n_recurrent)}'
    lines.append(f"{s.cohort}：{size}；复发减原发效应={s.effect:.6g}（{s.unit}）；P={s.p:.6g}；组内q={s.q_within_five:.6g}；6队列q={s.q_across_six_cohorts:.6g}。")
lines.extend(['','二、可支持与不能支持的解释',
'CGGA693：SIGLEC1与CD68均升高；SIGLEC1-CD68并未支持相对升高；四参照/两参照差值的P分别约0.821/0.800。因此仅支持组织层面的关联，不能排除髓系组成变化。',
'CGGA325：主要SIGLEC1检验未复现。两套CGGA仅为非配对比较，不能证明同一患者复发诱导。纳入样本ID没有交集，但缺少个人级跨队列去重证明。',
'CARE：121个标本、59名成人IDHwt患者，按相同QC重汇总；髓系主要检验48对，29对升高，但P=0.2148、q=0.3580，未达到显著。全核平均表达/髓系及全核伪bulk敏感性也未支持经校正的总体升高。'])
for _,s in program[program.program.eq('Interferon_response')].iterrows():
    lines.append(f"CARE {s.condition}：{int(s.n_specimens)}个可评估髓系标本；深度调整后SIGLEC1-干扰素中位rho={s.median_rho:.6g}；6程序组内q={s.q:.6g}。")
s=stage[stage.program.eq('Interferon_response')].iloc[0]
lines.append(f"CARE配对干扰素关联变化：{int(s.n_pairs)}对；变化={s.median_change:.6g}；P={s.p:.6g}；q={s.q:.6g}。两阶段都有相关性不等于复发特异性。")
lines.extend(['GSE212067：7对明确成人IDHwt原发/复发GBM，结果不支持SIGLEC1升高。',
'GSE4271：依据原作者互相匹配字段，选16对成人且两次均为历史WHO IV的组织；IDH缺失，不能视作现代IDHwt GBM验证。',
'GSE7696：采用原作者patient字段，不按样本名猜配对；3个有效配对。样本量极小，P=0.25未显著。bootstrap区间不应被当作可靠显著性证据。',
'官方唯一基因探针规则下，GSE4271/GSE7696的CD68映射不可用；CD68及SIGLEC1-CD68结果留空，不补零、不挑探针。','',
'三、检索中未作为新增独立验证的资源',
'GSE62153：既有项目已分析过，不计为新增数据源。',
'GSE228497：与原发现来源同属GSE228500，不能当独立复现。',
'GSE79671：主要为复发GBM贝伐珠单抗联合治疗前后比较，不是清晰的原发/复发验证。',
'GSE42670/GSE42669：前者为混合超级系列，后者含58份人GBM表达，尚需直接复发标签和患者配对核查；未给出SIGLEC1结果。',
'GLASS：表达数据仅覆盖部分样本，访问条款、来源重叠及IDH/时间点对照尚需审核；未完成本轮数值验证。','',
'四、数据与方法限制',
'所有新增检验均为修稿后探索性验证，不是预注册的确认性研究；BH不追溯修正此前的候选基因选择。',
'不将不同表达单位的效应直接合并。伪bulk log2(CPM+1)、细胞平均ln表达、FPKM和芯片归一化信号分别报告。',
'新队列按设计/覆盖/配对筛选，不按阳性结果筛选。更宽松阈值、患者剔除或单侧检验未用于制造结论。',
'跨数据源个人身份重叠未完全排除。CARE与GSE174554的标本来源和肿瘤ID未发现直接重合，但这不替代个人级去重。','',
'五、建议写入稿件的定位',
'可增加CGGA693作为组织层面关联的补充证据，同时完整报告CGGA325及三个组织配对队列未复现。',
'CARE可加强“SIGLEC1与干扰素相关髓系状态相对应”的跨队列证据；同时明确48对总体变化和47对阶段关联变化为阴性。',
'保留“候选、情境依赖性复发相关髓系重塑标志物”的探索性主线，不升级为“已验证的复发特异性标志物”。',
'建议停止无边界地寻找阳性数据库，下一步同步更新手稿和回复中的证据等级，并将所有新增检验放入补充材料。','',
'六、质量核查与文件',
f'主要与敏感性结果共{len(alltests)}行，其中{alltests.p.isna().sum()}项因缺失基因而不可估计；其余检验及BH值已独立复算。',
'原修稿投稿包SHA256保持不变；本轮没有更改或上传手稿。',
'统计与患者级数据在results目录；原始元数据/处理后矩阵/探针映射在input目录。',
'本轮仅下载公开元数据、处理后矩阵和必要平台注释；未下载FASTQ或CEL；CARE和CGGA大型计数复用已有缓存。','',
'公共记录：']+[f'{k}: {v}' for k,v in urls.items()])
(ROOT/'新增外部验证报告.txt').write_text('\n'.join(lines),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False,indent=2));print(primary[['cohort','n_primary','n_recurrent','p','q_within_five','q_across_six_cohorts']].to_string(index=False))
