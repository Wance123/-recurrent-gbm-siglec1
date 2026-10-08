
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import h5py,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from revision_analysis import frame
ROOT=Path(str(project_path('ROOT')));REV=ROOT/'Revision';TAB=REV/'analysis_tables';OUT=REV/'submission_package/figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
COL={'Primary':'#367A98','Recurrent':'#C27444','primary':'#367A98','recurrent':'#C27444'}
PAL=['#658BB4','#779E72','#D69B62','#B9B9B9','#B181AA','#D2685B','#5BA7A7']
NAMES=['0 SLC1A3 / MS4A6A','1 APOE / C1Q','2 PLIN2 / hypoxia','3 SOX2 / PTPRZ1 mixed','4 NAV3 / SORL1','5 CCL3 / FOS inflammatory','6 cycling']
def read(n):return pd.read_csv(TAB/(n+'.csv'))
def save(fig,n):
 from matplotlib.text import Text
 import re
 for obj in fig.findobj(match=Text):
  if obj.get_fontsize()<9.5:obj.set_fontsize(9.5)
  t=obj.get_text()
  if '$' not in t:
   t=t.replace('Mean ln(1 + counts / 10,000)','Mean ln(1 + normalized counts)')
   t=t.replace('SIGLEC1–CD68 log expression difference','SIGLEC1–CD68\nlog expression difference')
   t=t.replace('Median difference after patient deletion','Median difference\nafter patient deletion')
   t=re.sub(r'\b(SIGLEC1|CD68|CCL3|FOS|SOX2|PTPRZ1|SLC1A3|MS4A6A|APOE|C1Q|PLIN2|NAV3|SORL1|STAT1|CD24|KCNMB1)\b',lambda m:r'$\it{'+m[0]+'}$',t)
   obj.set_text(t)
 fig.savefig(OUT/(n+'.png'),dpi=350);fig.savefig(OUT/(n+'.pdf'));fig.savefig(OUT/(n+'.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
def letter(ax,s):ax.text(-.12,1.10,s.lower(),transform=ax.transAxes,fontweight='bold',fontsize=12)
def clean(ax):ax.grid(axis='y',color='#e7e7e7',lw=.5);ax.set_axisbelow(True)
def patient(ax,d,m,ylabel):
 for j,c in enumerate(['Primary','Recurrent']):
  vals=d.loc[d.condition==c,m].to_numpy();jit=np.linspace(-.10,.10,len(vals));ax.scatter(j+jit,vals,s=35,color=COL[c],edgecolor='white',lw=.6,zorder=3);ax.plot([j-.16,j+.16],[np.median(vals)]*2,color='black',lw=1.6)
 ax.set_xticks([0,1],['Primary','Recurrent']);ax.set_ylabel(ylabel);clean(ax)
def paired(ax,d,title,unit,p=None,q=None):
 for _,r in d.iterrows():
  color='#428B79' if r.recurrent>r.primary else '#aeb3b6';ax.plot([0,1],[r.primary,r.recurrent],c=color,alpha=.68,lw=1);ax.scatter([0,1],[r.primary,r.recurrent],s=17,c=color,zorder=3)
 ax.set_xticks([0,1],['Primary','Recurrent']);ax.set_ylabel(unit);ax.set_title(title,pad=24)
 if p is not None:ax.text(.5,1.02,f'n = {len(d)} pairs   P = {p:.3g}   q = {q:.3g}',ha='center',transform=ax.transAxes,fontsize=8)
 clean(ax)
with h5py.File(ROOT/'data/processed/GSE200984_myeloid_reclustered.h5ad') as f:xy=f['obsm/X_umap'][:];o=frame(f['obs'])
with h5py.File(ROOT/'data/processed/GSE200984_strict_basic_atlas_cnv.h5ad') as f:axxy=f['obsm/X_umap'][:];ao=frame(f['obs'])
fig,axs=plt.subplots(2,2,figsize=(7.2,6.2),layout='constrained');a=axs[0,0];a.axis('off');letter(a,'A')
steps=[('Discovery','GSE200984: 5 primary, 6 recurrent'),('Independent nuclear analysis','GSE174554: 17 complete pairs'),('Tissue level comparison','Two RNA-seq cohorts; NanoString'),('Exploratory sensitivity','Spatial n = 3; negative virtual KO')]
for j,(h,t) in enumerate(steps):
 yy=.96-j*.245;a.text(.04,yy,h,fontweight='bold',transform=a.transAxes,fontsize=9);a.text(.04,yy-.10,t,transform=a.transAxes,va='top',fontsize=9.5,linespacing=1.1)
a=axs[0,1];cl=ao.leiden_0_6.astype(int);maps={9:'Myeloid',12:'T / NK',10:'Oligodendrocyte',11:'Oligodendrocyte',13:'Oligodendrocyte',4:'Neuron'};lab=cl.map(maps).fillna('Other / tumor-like')
for c,color in zip(['Other / tumor-like','Myeloid','T / NK','Oligodendrocyte','Neuron'],['#bababa','#367A98','#917CA1','#83A580','#D3A86A']):
 ix=lab==c;a.scatter(axxy[ix,0],axxy[ix,1],s=.7,c=color,rasterized=True,label=c)
a.legend(loc='upper right',frameon=True,facecolor='white',framealpha=.95,edgecolor='white',markerscale=5,fontsize=7);a.set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2',title='Marker annotated atlas');letter(a,'B')
a=axs[1,0]
for k in range(7):
 ix=o.myeloid_leiden.astype(str)==str(k);a.scatter(xy[ix,0],xy[ix,1],s=2,c=PAL[k],label=NAMES[k],rasterized=True)
a.legend(loc='upper left',bbox_to_anchor=(0,-.12),frameon=False,fontsize=7,ncol=1,markerscale=3);a.set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2',title='Seven marker defined subclusters');letter(a,'C')
a=axs[1,1];d=read('S03_cluster_fraction_values');print('cluster columns',d.columns.tolist());d=d[d.subcluster.astype(str)=='5'];patient(a,d,'fraction','Fraction of retained myeloid nuclei');a.set_title('CCL3 / FOS inflammatory cluster\nP = 0.00866; q = 0.05195',pad=15);letter(a,'D');save(fig,'Figure_1')

d=read('S05_discovery_patient_metrics');d=d[d['filter']=='Exclude cluster 3'];tt=read('S05_discovery_inclusion_tests');tt=tt[tt['filter']=='Exclude cluster 3'];loo=read('S05_leave_one_out');loo=loo[loo.omitted_patient.notna()&(loo.omitted_patient!='None')];fig,axs=plt.subplots(2,3,figsize=(7.5,5.8),layout='constrained')
for j,(mt,label,old) in enumerate([('SIGLEC1','Mean ln(1 + counts / 10,000)','SIGLEC1'),('Detection','Detected fraction','SIGLEC1_detection'),('Log expression difference','SIGLEC1–CD68 log expression difference','SIGLEC1_minus_CD68')]):
 a=axs[0,j];patient(a,d,mt,label);r=tt[tt.metric==mt].iloc[0];a.set_title(f'P = {r.p:.3g}; q = {r.q:.3g}',pad=12);letter(a,chr(65+j));a=axs[1,j];l=loo[loo.metric==old];labs=[s+' '+d.set_index('sample').condition[s][0] for s in l.omitted_patient];a.scatter(l.median_difference_R_minus_P,np.arange(len(l)),s=22,color='#646c73');a.axvline(0,c='#999',ls='--',lw=.8);a.set_yticks(np.arange(len(l)),labs,fontsize=6.5);a.set_xlabel('Median difference after patient deletion');a.set_title('Leave one patient out\nP / R denotes disease stage',fontsize=9);letter(a,chr(68+j))
save(fig,'Figure_2')

vals=read('S11_all_paired_values');tt=read('S11_external_tests');fig,axs=plt.subplots(3,3,figsize=(7.5,8.2),layout='constrained')
for j,ds in enumerate(['GSE339484','GSE139533','GSE254873']):
 for i,mt in enumerate(['SIGLEC1','CD68','SIGLEC1_CD68_log_difference']):
  a=axs[i,j];d=vals[(vals.dataset==ds)&(vals.metric==mt)];r=tt[(tt.dataset==ds)&(tt.metric==mt)].iloc[0];unit='Log2(CPM + 1)' if j<2 else 'Deposited normalized log2 CPM';
  if i==2:unit='SIGLEC1–CD68 log expression difference'
  paired(a,d,ds+'\n'+['SIGLEC1','CD68','SIGLEC1–CD68'][i],unit,r.p,r.q);letter(a,chr(65+i*3+j))
save(fig,'Figure_3')

vals=read('S12_external_paired_values');tt=read('S12_external_tests');fig,axs=plt.subplots(2,3,figsize=(7.5,5.8),layout='constrained')
for j,(mt,label) in enumerate([('SIGLEC1','Mean ln(1 + counts / 10,000)'),('CD68','Mean ln(1 + counts / 10,000)'),('Detection','Detected fraction'),('Log_expression_difference','SIGLEC1–CD68 log expression difference'),('myeloid_fraction','Fraction of all QC nuclei')]):
 a=axs.flat[j];d=vals[(vals.origin_proxy=='All myeloid')&(vals.metric==mt)];r=tt[(tt.origin_proxy=='All myeloid')&(tt.metric==mt)].iloc[0];paired(a,d,mt.replace('_',' '),label,r.p,r.q);letter(a,chr(65+j))
a=axs.flat[5];a.axis('off');letter(a,'F');a.text(0,.9,'GSE174554 eligibility',weight='bold');a.text(0,.80,'68 adult IDH wild type specimens\n34 original clinical pairs\n\n12 incomplete libraries excluded\n6 with fewer than 20 myeloid nuclei\n\n50 specimens for context\n17 pairs for recurrence tests\n14,683 paired myeloid nuclei',va='top',linespacing=1.35,fontsize=8);save(fig,'Figure_4')

fig,axs=plt.subplots(2,2,figsize=(7.5,6.8),layout='constrained');progs=['Antigen presentation','Interferon response','Immunoregulatory','Lipid and complement','Phagolysosome','Immediate inflammatory'];short=['Antigen presentation','Interferon','Immunoregulatory','Lipid / complement','Phagolysosome','Inflammatory'];
for j,stage in enumerate(['Primary','Recurrent']):
 a=axs[0,j];d=read('S07_program_patient_effects');d=d[(d.condition==stage)&(d.threshold==0)];
 for i,p in enumerate(progs):a.scatter(d[d.program==p].standardized_difference,np.linspace(i-.10,i+.10,len(d[d.program==p])),s=25,color=COL[stage])
 a.axvline(0,c='#999',lw=.8);a.set_yticks(range(6),short);a.invert_yaxis();a.set_xlabel('Within specimen standardized difference');a.set_title('Discovery '+stage.lower()+' (n = 4)\nAll six programs q > 0.05');letter(a,chr(65+j))
 a=axs[1,j];d=read('S12_external_continuous_correlations');d=d[d.condition==stage.lower()];ts=read('S12_external_continuous_tests');ts=ts[(ts.condition==stage.lower())&(ts.metric=='depth_adjusted_rank_correlation')];
 for i,p in enumerate(progs):
  v=d[d.program==p].depth_adjusted_rank_correlation;a.scatter(v,np.linspace(i-.13,i+.13,len(v)),s=12,alpha=.45,color=COL[stage]);a.plot([v.median()]*2,[i-.18,i+.18],c='black',lw=1.5);q=ts[ts.program==p].q.iloc[0];a.text(1.02,i,f'q={q:.3g}',transform=a.get_yaxis_transform(),fontsize=7)
 a.axvline(0,c='#999',lw=.8);a.set_yticks(range(6),short);a.invert_yaxis();a.set_xlim(-.25,.65);a.set_xlabel('Depth adjusted within specimen rank correlation');a.set_title(f'External {stage.lower()}\n'+('19' if stage=='Primary' else '29')+' specimens with variable SIGLEC1');letter(a,chr(67+j))
save(fig,'Figure_5')

fig,axs=plt.subplots(2,2,figsize=(7.5,6.3),layout='constrained');d=read('S03_patient_cluster_distribution');a=axs[0,0];m=d.pivot(index='sample',columns='cluster',values='fraction_of_sample_myeloid');im=a.imshow(m,aspect='auto',cmap='Blues',vmin=0,vmax=.8);a.set_xticks(range(7),range(7));a.set_yticks(range(len(m)),m.index,fontsize=7);a.set_xlabel('Subcluster');a.set_title('Fraction of each patient myeloid compartment');fig.colorbar(im,ax=a,fraction=.04);letter(a,'A')
a=axs[0,1];m=d.pivot(index='sample',columns='cluster',values='fraction_of_cluster');im=a.imshow(m,aspect='auto',cmap='Oranges',vmin=0,vmax=.5);a.set_xticks(range(7),range(7));a.set_yticks(range(len(m)),m.index,fontsize=7);a.set_xlabel('Subcluster');a.set_title('Patient contribution to each subcluster');fig.colorbar(im,ax=a,fraction=.04);letter(a,'B')
q=read('S04_cell_QC');
for a,col,l in [(axs[1,0],'doublet_score','Scrublet doublet score'),(axs[1,1],'cnv_rms','Expression derived CNV RMS')]:
 v=[q[q.cluster.astype(str)==str(k)][col] for k in range(7)];bp=a.boxplot(v,tick_labels=range(7),showfliers=False,patch_artist=True,widths=.65)
 for k,b in enumerate(bp['boxes']):b.set_facecolor(PAL[k]);b.set_alpha(.7)
 a.set_ylabel(l);a.set_xlabel('Subcluster');a.set_title('Cluster 3 is ambiguous, not CNV confirmed');clean(a)
letter(axs[1,0],'C');letter(axs[1,1],'D');save(fig,'Supplementary_Figure_S1')

fig,axs=plt.subplots(1,3,figsize=(7.5,3.5),layout='constrained');v=read('S11_all_paired_values')
for j,mt in enumerate(['SIGLEC1','CD68','SIGLEC1_CD68_log_difference']):
 a=axs[j];d=v[(v.dataset=='GSE254145')&(v.metric==mt)];paired(a,d,mt.replace('_',' ') if j<2 else 'Log expression difference','Deposited log2 CPM' if j<2 else 'Log expression difference');a.text(.5,1.02,'n = 3; descriptive only',ha='center',transform=a.transAxes,fontsize=8);letter(a,chr(65+j))
save(fig,'Supplementary_Figure_S2')

fig,axs=plt.subplots(1,3,figsize=(7.5,3.6),layout='constrained');a=axs[0];d=read('S09_power_sensitivity');a.plot(d.standardized_normal_shift,d.power,c='#367A98',lw=2);a.axhline(.8,c='#888',ls='--');a.axvline(2.2,c='#aaa',ls=':');a.set(xlabel='Standardized normal mean shift',ylabel='Simulated rejection probability',title='5 + 6 unpaired design');letter(a,'A')
z=read('S08_patient_composition_CLR')
for a,mt,l in [(axs[1],'Myeloid_fraction','Myeloid fraction'),(axs[2],'Tumor_CLR','Other / tumor-like CLR')]:
 for st,g in z.groupby('condition'):a.scatter(g.SIGLEC1,g[mt],c=COL[st],s=30,label=st)
 a.set(xlabel='Patient mean SIGLEC1',ylabel=l);a.legend(frameon=False,fontsize=7)
letter(axs[1],'B');letter(axs[2],'C');save(fig,'Supplementary_Figure_S3')

fig,axs=plt.subplots(1,2,figsize=(7.3,3.4),layout='constrained');d=read('S13_sctenifoldknk_run_summary_corrected');print('KO columns',d.columns.tolist());s=read('S13_stability_stability_runs');print('KO stability',s.columns.tolist());
axs[0].axis('off');axs[0].text(0,.92,'Exploratory network reconstruction',weight='bold');axs[0].text(0,.80,'744 genes; 1,417 balanced nuclei\nPrimary 558; recurrent 859\n5 networks × 300 nuclei per subsample\nPCR inference; no regulatory prior\n3 tensor components; 200 iterations\nAnalytical gene tests; BH q < 0.05',va='top',linespacing=1.5)
axs[1].axis('off');axs[1].text(0,.92,'Negative multi seed sensitivity',weight='bold');axs[1].text(0,.80,'SIGLEC1 and STAT1 examined in 3 seeds\nSeeds 34, 134 and 234\nRecurrent: 0 genes in at least 2 seeds\nPrimary: CD24 and KCNMB1 shared\nby both targets, not target specific\n\nTwo of three rule added after screening\nNo experimentally measured knockout',va='top',linespacing=1.4);letter(axs[0],'A');letter(axs[1],'B');save(fig,'Supplementary_Figure_S4')
print('FIGURES COMPLETE')
