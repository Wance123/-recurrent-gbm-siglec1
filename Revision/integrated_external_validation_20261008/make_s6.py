
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).parent;AN=ROOT.parent/'external_validation_expansion_20261008';OUT=ROOT/'submission_package/figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(2,3,figsize=(10.8,6.8));rng=np.random.default_rng(34);colors=['#507B9C','#C3755B']
for ax,cohort in zip(axs[0,:2],['CGGA693','CGGA325']):
 a=pd.read_csv(AN/'results'/f'{cohort}_specimen_values.csv')
 for j,stage in enumerate(['primary','recurrent']):
  v=a.loc[a.condition.eq(stage),'SIGLEC1'].to_numpy();ax.scatter(j+rng.uniform(-.14,.14,len(v)),v,s=9,alpha=.55,color=colors[j],linewidths=0);ax.plot([j-.2,j+.2],[np.median(v)]*2,color='black',lw=1.7)
 ax.set_xticks([0,1],['Primary','Recurrent']);ax.set_ylabel('SIGLEC1 log2(FPKM + 1)');ax.set_title(cohort+' tissue',loc='left',fontweight='bold');ax.text(.5,1.03,'n=108 / 79; q=0.000586' if cohort=='CGGA693' else 'n=73 / 14; q=0.917',transform=ax.transAxes,ha='center',fontsize=8)
pair=pd.read_csv(AN/'results/CARE_TAM_mean_cell_SIGLEC1_paired_values.csv');ax=axs[0,2]
for _,r in pair.iterrows():ax.plot([0,1],[r.primary,r.recurrent],color='#B9C2CB',lw=.6,zorder=1)
for j,s in enumerate(['primary','recurrent']):ax.scatter(np.full(len(pair),j),pair[s],s=12,color=colors[j],zorder=2)
ax.set_xticks([0,1],['Primary','Recurrent']);ax.set_ylabel('Mean TAM SIGLEC1 ln expression');ax.set_title('CARE paired expression',loc='left',fontweight='bold');ax.text(.5,1.03,'48 pairs; P=0.215; q=0.358',transform=ax.transAxes,ha='center',fontsize=8)
a=pd.read_csv(AN/'results/CARE_TAM_mean_cell_values.csv');pr='depth_rho_Interferon_response';ax=axs[1,0]
for j,s in enumerate(['primary','recurrent']):
 v=a.loc[a.condition.eq(s),pr].dropna();ax.scatter(j+rng.uniform(-.12,.12,len(v)),v,s=13,color=colors[j],alpha=.75);ax.plot([j-.2,j+.2],[v.median()]*2,color='black',lw=1.7)
ax.axhline(0,color='#AAB1B7',ls=':',lw=.8);ax.set_xticks([0,1],['Primary','Recurrent']);ax.set_ylabel('UMI-adjusted SIGLEC1–IFN correlation');ax.set_title('CARE interferon association',loc='left',fontweight='bold');ax.text(.5,1.03,'n=51 / 54; q=3.48e−9 / 9.75e−10',transform=ax.transAxes,ha='center',fontsize=8)
w=a.pivot(index='patient',columns='condition',values=pr).dropna();ax=axs[1,1]
for _,r in w.iterrows():ax.plot([0,1],[r.primary,r.recurrent],color='#B9C2CB',lw=.6)
for j,s in enumerate(['primary','recurrent']):ax.scatter(np.full(len(w),j),w[s],s=12,color=colors[j])
ax.axhline(0,color='#AAB1B7',ls=':',lw=.8);ax.set_xticks([0,1],['Primary','Recurrent']);ax.set_ylabel('UMI-adjusted SIGLEC1–IFN correlation');ax.set_title('CARE paired association',loc='left',fontweight='bold');ax.text(.5,1.03,'47 pairs; P=0.929; q=0.954',transform=ax.transAxes,ha='center',fontsize=8)
ax=axs[1,2];change=(w.recurrent-w.primary).to_numpy();ax.scatter(np.arange(len(change)),np.sort(change),s=15,color='#507B9C');ax.axhline(0,color='#AAB1B7',ls=':',lw=.8);ax.set_xlabel('Patient rank by paired change');ax.set_ylabel('Recurrent − primary IFN correlation');ax.set_title('CARE individual changes',loc='left',fontweight='bold');ax.text(.5,1.03,'Median change −0.006',transform=ax.transAxes,ha='center',fontsize=8)
for i,ax in enumerate(axs.flat):
 ax.text(-.14,1.24,'abcdef'[i],transform=ax.transAxes,fontweight='bold',fontsize=12)
 ax.set_title(ax.get_title(loc='left'),loc='left',fontweight='bold',y=1.12,pad=8)
 if i!=5:ax.set_xlim(-.3,1.3)
fig.subplots_adjust(left=.075,right=.965,bottom=.09,top=.84,wspace=.45,hspace=.75);fp=OUT/'Supplementary_Figure_S6'
fig.savefig(str(fp)+'.png',dpi=350,facecolor='white');fig.savefig(str(fp)+'.pdf',facecolor='white');fig.savefig(str(fp)+'.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
print('Scientific S6 plot generated using existing measured data; bundled runtime lacked matplotlib.')
