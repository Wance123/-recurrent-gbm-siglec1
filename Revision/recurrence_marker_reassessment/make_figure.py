
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
ROOT=Path(__file__).parent;OUT=ROOT/'submission_package/figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
a=pd.read_csv(ROOT/'analysis_tables/S14_paired_program_association_values.csv');a=a[(a.program=='Interferon response')&(a.metric=='depth_adjusted_rank_correlation')]
b=pd.read_csv(ROOT/'analysis_tables/S14_patient_longitudinal_changes.csv')
fig,ax=plt.subplots(1,2,figsize=(7.2,4.2));fig.subplots_adjust(left=.11,right=.98,bottom=.21,top=.77,wspace=.52)
for _,r in a.iterrows():ax[0].plot([0,1],[r.primary,r.recurrent],color='#A0A6AC',alpha=.65,lw=.8)
ax[0].scatter(np.zeros(len(a)),a.primary,color='#346E91',s=23,zorder=3);ax[0].scatter(np.ones(len(a)),a.recurrent,color='#C27457',s=23,zorder=3)
ax[0].set_xticks([0,1],['Primary','Recurrent']);ax[0].set_xlim(-.3,1.3);ax[0].axhline(0,c='gray',lw=.65,ls='--');ax[0].set_ylabel('Depth adjusted rank correlation')
ax[0].set_title('Paired interferon association\n15 pairs   P = 0.489   q = 0.733',fontsize=11,pad=12)
ax[1].scatter(b.SIGLEC1,b['Interferon response'],s=32,c='#346E91',edgecolors='white',linewidths=.4)
ax[1].axhline(0,c='gray',lw=.65,ls='--');ax[1].axvline(0,c='gray',lw=.65,ls='--');ax[1].set_xlabel('Change in mean SIGLEC1');ax[1].set_ylabel('Change in mean interferon score')
ax[1].set_title('Longitudinal coupling\n17 pairs   rho = 0.282   q = 0.544',fontsize=11,pad=12)
for lab,x in zip(['a','b'],ax):x.text(-.22,1.26,lab,transform=x.transAxes,fontweight='bold',fontsize=14)
fig.text(.5,.025,'Changes are recurrent minus primary. No corrected positive association.',ha='center',fontsize=10)
for ext,dpi in [('png',350),('pdf',350),('tif',600)]:
 kwargs={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tif' else {}
 fig.savefig(OUT/('Supplementary_Figure_S5.'+ext),dpi=dpi,**kwargs)
print('FIGURE S5 COMPLETE')
