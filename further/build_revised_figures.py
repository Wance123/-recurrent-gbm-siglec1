
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import warnings
import h5py
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

warnings.filterwarnings("ignore")
ROOT = Path(str(project_path('ROOT')))
OUT = ROOT / "further" / "figures_revised"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#277DA1"; ORANGE = "#F28E2B"; TEAL = "#43AA8B"; RED = "#D9534F"
PURPLE = "#7A5195"; GOLD = "#E3B341"; GREY = "#6B7280"; LIGHT = "#E8EEF3"
COND = {"Primary": BLUE, "Recurrent": ORANGE}
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "legend.fontsize": 7, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
})
sns.set_style("ticks")

def panel(ax, label):
    ax.text(-0.08, 1.17, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")

def clean(ax):
    sns.despine(ax=ax)
    ax.grid(axis="y", color="#D9E1E8", linewidth=0.45, alpha=0.8, zorder=0)

def save(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.tiff", dpi=300, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

def read_h5_categorical(group, key):
    obj = group[key]
    if isinstance(obj, h5py.Dataset):
        vals = obj[:]
        return np.array([x.decode() if isinstance(x, bytes) else x for x in vals])
    codes = obj["codes"][:]
    cats = obj["categories"][:]
    cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in cats])
    return np.array([cats[i] if i >= 0 else "NA" for i in codes])

def paired_panel(ax, df, left, right, title, ylabel, ptext):
    for _, r in df.iterrows():
        color = TEAL if r[right] > r[left] else GREY
        ax.plot([0,1], [r[left],r[right]], color=color, alpha=0.72, lw=0.9, zorder=1)
        ax.scatter([0,1], [r[left],r[right]], c=[BLUE,ORANGE], s=15,
                   edgecolor="white", linewidth=0.35, zorder=2)
    ax.set_xticks([0,1], ["Primary","Recurrent"])
    ax.set_title(title, fontweight="bold", pad=9, fontsize=8.5)
    ax.set_ylabel(ylabel)
    ax.margins(y=.15); ax.text(0.5, 0.95, ptext, transform=ax.transAxes, ha="center", va="top", fontsize=6.4, bbox=dict(facecolor="white", edgecolor="none", alpha=.82, pad=.5))
    clean(ax)

# Figure 1: study design and atlas
fig = plt.figure(figsize=(7.15, 6.2))
gs = fig.add_gridspec(2, 3, height_ratios=[0.75, 1.25], hspace=0.32, wspace=0.28)
ax = fig.add_subplot(gs[0,:]); ax.axis("off"); ax.text(0.0,0.98,"A",transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
steps = [
    ("Discovery atlas", "GSE200984\n5 primary + 6 recurrent"),
    ("Patient-aware tests", "Myeloid states\nSIGLEC1 metrics + LOO"),
    ("External validation", "3 paired cohorts\n+ spatial compartments"),
    ("Mechanistic stress test", "State programs + ecology\nvirtual knockout stability"),
]
xs = np.linspace(0.02,0.76,4)
for i,(title,sub) in enumerate(steps):
    box=FancyBboxPatch((xs[i],0.28),0.19,0.46,boxstyle="round,pad=0.012,rounding_size=0.02",
                       fc=["#E8F1F8","#EEF6F2","#FFF1E6","#F2ECF7"][i],ec="#AAB7C4",lw=0.8)
    ax.add_patch(box); ax.text(xs[i]+0.095,0.61,title,ha="center",va="center",fontweight="bold",fontsize=8.5)
    ax.text(xs[i]+0.095,0.43,sub,ha="center",va="center",fontsize=7.5,color="#374151")
    if i<3: ax.add_patch(FancyArrowPatch((xs[i]+0.195,0.51),(xs[i+1]-0.008,0.51),arrowstyle="-|>",mutation_scale=10,lw=0.9,color="#667788"))
ax.text(0.98,0.06,"Inference unit: patient, not cell",ha="right",va="bottom",fontsize=7.5,color="#374151",fontweight="bold")

cluster_map = {0:"Tumor",1:"Tumor",2:"Tumor",3:"Tumor",4:"Neuron",5:"Tumor",6:"Tumor",7:"Tumor",8:"Tumor",9:"Myeloid",10:"Oligodendrocyte",11:"Oligodendrocyte",12:"T/NK",13:"Oligodendrocyte",14:"Tumor",15:"Tumor",16:"Tumor",17:"Tumor"}
lineage_colors={"Tumor":"#B9B2D8","Myeloid":"#D95F59","T/NK":"#4C9F70","Oligodendrocyte":"#4C78A8","Neuron":"#E3B341"}
with h5py.File(ROOT/"data/processed/GSE200984_strict_basic_atlas_cnv.h5ad","r") as f:
    xy=f["obsm/X_umap"][:]; cl=read_h5_categorical(f["obs"],"leiden_0_6"); condition=read_h5_categorical(f["obs"],"condition")
lineage=np.array([cluster_map.get(int(float(x)),"Other") for x in cl])
ax=fig.add_subplot(gs[1,0]); panel(ax,"B")
for lab in ["Tumor","Myeloid","Oligodendrocyte","T/NK","Neuron"]:
    m=lineage==lab; ax.scatter(xy[m,0],xy[m,1],s=1.0,c=lineage_colors[lab],label=lab,alpha=0.65,rasterized=True,linewidths=0)
ax.set_title("Single-nucleus atlas",fontweight="bold"); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_xticks([]); ax.set_yticks([]); sns.despine(ax=ax,left=True,bottom=True)
ax.legend(frameon=False,markerscale=4,loc="upper left",bbox_to_anchor=(0,-0.08),ncol=2,handletextpad=0.3,columnspacing=0.7)

with h5py.File(ROOT/"data/processed/GSE200984_myeloid_reclustered.h5ad","r") as f:
    mxy=f["obsm/X_umap"][:]; state=read_h5_categorical(f["obs"],"myeloid_leiden")
ax=fig.add_subplot(gs[1,1]); panel(ax,"C")
palette=sns.color_palette("colorblind",7)
for i in range(7):
    m=state==str(i); lab="State 5: CCL3/FOS" if i==5 else f"State {i}"
    ax.scatter(mxy[m,0],mxy[m,1],s=2.2,c=[palette[i]],label=lab,alpha=0.72,rasterized=True,linewidths=0)
ax.set_title("Myeloid state map",fontweight="bold"); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_xticks([]); ax.set_yticks([]); sns.despine(ax=ax,left=True,bottom=True)
ax.legend(frameon=False,markerscale=3,loc="upper left",bbox_to_anchor=(0,-0.08),ncol=2,handletextpad=0.25,columnspacing=0.5)

comp=pd.read_csv(ROOT/"results/tables/GSE200984_clean_myeloid_subcluster_patient_composition.csv")
col_state=[c for c in comp.columns if c.lower() in {"subcluster","myeloid_leiden","cluster"}][0]
col_frac=[c for c in comp.columns if "fraction" in c.lower()][0]
sub=comp[comp[col_state].astype(str)=="5"].copy()
ax=fig.add_subplot(gs[1,2]); panel(ax,"D")
for j,cond in enumerate(["Primary","Recurrent"]):
    vals=sub.loc[sub.condition==cond,col_frac].astype(float)*100
    rng=np.random.default_rng(42+j); x=j+rng.uniform(-0.07,0.07,len(vals))
    ax.scatter(x,vals,s=28,c=COND[cond],edgecolor="white",linewidth=0.5,zorder=3,label=cond)
    med=np.median(vals); ax.plot([j-0.18,j+0.18],[med,med],color="#222",lw=1.4,zorder=4)
ax.set_xticks([0,1],["Primary","Recurrent"]); ax.set_ylabel("CCL3/FOS state (% myeloid)")
ax.set_title("Immediate-inflammatory state",fontweight="bold"); ax.text(0.5,0.97,"P=0.00866; FDR=0.052",transform=ax.transAxes,ha="center",va="top",fontsize=7)
clean(ax)
fig.suptitle("Study design and patient-aware single-nucleus atlas",fontsize=11,fontweight="bold",y=0.995)
save(fig,"Fig1_revised")

# Figure 2: discovery cohort and sensitivity
metrics=pd.read_csv(ROOT/"results/tables/GSE200984_clean_myeloid_SIGLEC1_CD68_patient_metrics.csv")
loo=pd.read_csv(ROOT/"further/results/reinforcement/leave_one_patient_out.csv")
fig,axes=plt.subplots(2,3,figsize=(7.15,5.3)); plt.subplots_adjust(wspace=0.38,hspace=0.48)
vars_=[("SIGLEC1","Mean SIGLEC1 expression"),("SIGLEC1_detection","SIGLEC1+ fraction"),("SIGLEC1_minus_CD68","SIGLEC1 - CD68")]
if "SIGLEC1_detection" not in metrics.columns:
    pm=pd.read_csv(ROOT/"further/results/reinforcement/patient_SIGLEC1_metrics.csv")
    for c in ["SIGLEC1_detection"]:
        if c in pm.columns: metrics=metrics.merge(pm[["sample",c]],on="sample",how="left")
for k,(var,ylabel) in enumerate(vars_):
    ax=axes[0,k]; panel(ax,chr(65+k))
    for j,cond in enumerate(["Primary","Recurrent"]):
        vals=metrics.loc[metrics.condition==cond,var].dropna().astype(float)
        rng=np.random.default_rng(100+k*2+j); xx=j+rng.uniform(-0.08,0.08,len(vals))
        ax.scatter(xx,vals,s=34,c=COND[cond],edgecolor="white",linewidth=0.6,zorder=3)
        ax.plot([j-0.17,j+0.17],[np.median(vals)]*2,color="#222",lw=1.5,zorder=4)
    ax.set_xticks([0,1],["Primary","Recurrent"]); ax.set_ylabel(ylabel); clean(ax)
    ax.set_title(["Expression","Detection","CD68-adjusted metric"][k],fontweight="bold")
    ax.text(.5,.97,["P=0.855","P=0.855","P=1.000"][k],transform=ax.transAxes,ha="center",va="top",fontsize=7)
for k,(var,_) in enumerate(vars_):
    ax=axes[1,k]; panel(ax,chr(68+k)); dd=loo[(loo.metric==var)&loo.omitted_patient.notna()].copy()
    dd=dd.sort_values("median_difference_R_minus_P"); y=np.arange(len(dd))
    ax.axvline(0,color="#555",lw=0.8,ls="--"); ax.scatter(dd.median_difference_R_minus_P,y,s=25,c=TEAL,edgecolor="white",linewidth=.4,zorder=3)
    ax.hlines(y,0,dd.median_difference_R_minus_P,color=TEAL,lw=1)
    ax.set_yticks(y,dd.omitted_patient.str.replace("TB","",regex=False)); ax.set_xlabel("Recurrent - primary median")
    ax.set_ylabel("Omitted patient"); ax.set_title("Leave-one-patient-out",fontweight="bold"); clean(ax)
fig.suptitle("Discovery-cohort SIGLEC1 signals are directional but uncertain",fontsize=11,fontweight="bold",y=.99)
save(fig,"Fig2_revised")

# Figure 3: external paired validation
deltas=pd.read_csv(ROOT/"results/tables/SIGLEC1_CD68_cross_cohort_patient_deltas.csv")
g339=deltas[deltas.dataset=="GSE339484"].rename(columns={"primary_log_ratio":"primary","recurrent_log_ratio":"recurrent"})
g139=deltas[deltas.dataset=="GSE139533"].rename(columns={"primary_log_ratio":"primary","recurrent_log_ratio":"recurrent"})
g254=pd.read_csv(ROOT/"results/tables/GSE254873_IFN_myeloid_paired_scores.csv").rename(columns={"primary_score":"primary","recurrent_score":"recurrent"})
sp=pd.read_csv(ROOT/"results/tables/GSE254145_spatial_candidate_paired_values.csv"); sp=sp[(sp.compartment=="CD64")&(sp.gene=="SIGLEC1")]
fig=plt.figure(figsize=(7.15,5.3)); gs=fig.add_gridspec(2,4,height_ratios=[1.1,.9],hspace=.42,wspace=.45)
for i,(df,title,ptext,ylabel) in enumerate([(g339,"GSE339484","10/11 increased; P=0.00684","log SIGLEC1/CD68"),(g139,"GSE139533","11/13 increased; P=0.00342","log SIGLEC1/CD68"),(g254,"GSE254873","11/17 increased; P=0.0129","absolute SIGLEC1"),(sp,"GSE254145 CD64+","2/3 increased; P=0.75","spatial SIGLEC1")]):
    ax=fig.add_subplot(gs[0,i]); panel(ax,chr(65+i)); paired_panel(ax,df,"primary","recurrent",title,ylabel,ptext)
ax=fig.add_subplot(gs[1,:2]); panel(ax,"E")
summ=[]
for name,df in [("GSE339484\nratio",g339),("GSE139533\nratio",g139),("GSE254873\nabsolute",g254),("GSE254145\nspatial",sp)]:
    v=(df.recurrent-df.primary).astype(float); summ.append((name,np.median(v),np.mean(v>0),len(v)))
sd=pd.DataFrame(summ,columns=["cohort","median","positive_fraction","n"]); y=np.arange(len(sd))
ax.axvline(0,color="#555",lw=.8,ls="--"); ax.scatter(sd["median"],y,s=45,c=[TEAL,TEAL,ORANGE,PURPLE],edgecolor="white",linewidth=.5,zorder=3)
ax.hlines(y,0,sd["median"],color=[TEAL,TEAL,ORANGE,PURPLE],lw=2); ax.set_yticks(y,sd.cohort); ax.invert_yaxis(); ax.set_xlabel("Median paired change (recurrent - primary)"); ax.set_title("Direction and magnitude",fontweight="bold"); clean(ax)
ax=fig.add_subplot(gs[1,2:]); panel(ax,"F")
bars=ax.barh(y,sd.positive_fraction*100,color=[TEAL,TEAL,ORANGE,PURPLE],height=.55)
for b,n in zip(bars,sd.n): ax.text(b.get_width()+2,b.get_y()+b.get_height()/2,f"n={n}",va="center",fontsize=7)
ax.axvline(50,color="#555",lw=.8,ls="--"); ax.set_yticks(y,sd.cohort); ax.invert_yaxis(); ax.set_xlim(0,110); ax.set_xlabel("Patients with increase (%)"); ax.set_title("Patient-level concordance",fontweight="bold"); clean(ax)
fig.suptitle("Independent paired cohorts validate context-dependent SIGLEC1 remodeling",fontsize=11,fontweight="bold",y=.99)
save(fig,"Fig3_revised")

# Figure 4: state programs and ecology
state=pd.read_csv(ROOT/"further/results/reinforcement/SIGLEC1_positive_state_patient_deltas.csv")
stats=pd.read_csv(ROOT/"further/results/reinforcement/SIGLEC1_positive_state_statistics.csv")
eco=pd.read_csv(ROOT/"further/results/reinforcement/patient_ecology_metrics.csv")
corr=pd.read_csv(ROOT/"further/results/reinforcement/patient_ecology_correlations.csv").set_index("outcome")
paths=["Antigen_presentation","Interferon_response","Immunoregulatory","Lipid_C1Q_APOE","Phagolysosome","Immediate_inflammatory"]
fig=plt.figure(figsize=(7.15,6.0)); gs=fig.add_gridspec(2,4,height_ratios=[1.05,1],hspace=.45,wspace=.45)
for idx,cond in enumerate(["Primary","Recurrent"]):
    ax=fig.add_subplot(gs[0,idx*2:(idx+1)*2]); panel(ax,chr(65+idx)); sub=state[state.condition==cond]
    mat=sub.pivot(index="sample",columns="pathway",values="delta_positive_minus_negative").reindex(columns=paths)
    vmax=max(.35,np.nanmax(np.abs(mat.values))); sns.heatmap(mat,cmap="vlag",center=0,vmin=-vmax,vmax=vmax,linewidths=.5,linecolor="white",cbar=idx==1,cbar_kws={"label":"SIGLEC1+ minus SIGLEC1-"},ax=ax)
    ax.set_title(f"{cond}: within-patient state differences",fontweight="bold"); ax.set_xlabel(""); ax.set_ylabel("Patient"); ax.set_xticklabels([p.replace("_","\n") for p in paths],rotation=35,ha="right")
outcomes=[("malignant_MES","Malignant MES score"),("T_exhaustion","T/NK exhaustion"),("malignant_fraction","Malignant fraction"),("tnk_fraction","T/NK fraction")]
for i,(outcome,label) in enumerate(outcomes):
    ax=fig.add_subplot(gs[1,i]); panel(ax,chr(67+i)); dd=eco[["condition","SIGLEC1",outcome]].dropna()
    for cond in ["Primary","Recurrent"]:
        s=dd[dd.condition==cond]; ax.scatter(s.SIGLEC1,s[outcome],s=28,c=COND[cond],edgecolor="white",linewidth=.5,label=cond,zorder=3)
    if len(dd)>2:
        z=np.polyfit(dd.SIGLEC1,dd[outcome],1); xx=np.linspace(dd.SIGLEC1.min(),dd.SIGLEC1.max(),50); ax.plot(xx,np.polyval(z,xx),color="#555",lw=1,ls="--")
    rr=corr.loc[outcome]; ax.text(.04,.96,f"rho={rr.spearman_rho:.2f}\nP={rr.p_value:.3f}\nFDR={rr.FDR:.3f}",transform=ax.transAxes,va="top",fontsize=6.8)
    ax.set_xlabel("Patient SIGLEC1"); ax.set_ylabel(label); clean(ax)
    if i==0: ax.legend(frameon=False,loc="lower left")
fig.suptitle("SIGLEC1-positive myeloid cells show a mixed state without a universal ecosystem axis",fontsize=11,fontweight="bold",y=.995)
save(fig,"Fig4_revised")

# Figure 5: virtual knockout and stability
run=pd.read_csv(ROOT/"further/results/sctenifoldknk/run_summary_corrected.csv")
stab=pd.read_csv(ROOT/"further/results/stability/stability_runs.csv")
rec=pd.read_csv(ROOT/"further/results/stability/gene_recurrence_across_seeds.csv")
fig=plt.figure(figsize=(7.15,4.25)); gs=fig.add_gridspec(1,3,wspace=.62)
ax=fig.add_subplot(gs[0,0]); panel(ax,"A")
pv=run.pivot(index="target",columns="condition",values="downstream_FDR_0_05").fillna(0); pv=pv.reindex(sorted(pv.index))
sns.heatmap(pv,annot=True,fmt=".0f",cmap=sns.light_palette(PURPLE,as_cmap=True),cbar_kws={"label":"FDR-positive genes"},linewidths=.6,linecolor="white",ax=ax)
ax.set_title("Initial screen",fontweight="bold",pad=8); ax.set_xlabel(""); ax.set_ylabel("Target")
ax=fig.add_subplot(gs[0,1]); panel(ax,"B")
ss=stab.copy(); ss["group"]=ss.condition+" / "+ss.target; order=["Primary / SIGLEC1","Recurrent / SIGLEC1","Primary / STAT1","Recurrent / STAT1"]
for yi,g in enumerate(order):
    q=ss[ss.group==g].sort_values("seed"); ax.plot(q.seed,q.n_sig,marker=["o","o","s","s"][yi],ms=4.5,lw=1.1,color=[BLUE,ORANGE,BLUE,ORANGE][yi],label=g.replace(" / ","-"))
ax.set_xlabel("Random seed"); ax.set_ylabel("FDR-positive genes"); ax.set_title("Multi-seed results",fontweight="bold",pad=8); ax.legend(frameon=False,fontsize=5.8,loc="upper center",bbox_to_anchor=(.5,-.16),ncol=2,handlelength=1.3,columnspacing=.7); clean(ax)
ax=fig.add_subplot(gs[0,2]); panel(ax,"C")
rr=rec[rec.seeds_detected_of_3>=2].copy(); rr["label"]=rr.target+" -> "+rr.gene
if len(rr):
    y=np.arange(len(rr)); colors=[BLUE if x=="Primary" else ORANGE for x in rr.condition]
    ax.barh(y,rr.seeds_detected_of_3,color=colors,height=.58); ax.set_yticks(y,rr.label); ax.invert_yaxis()
ax.set_xlim(0,3.25); ax.set_xticks([0,1,2,3]); ax.set_xlabel("Seeds detected (of 3)"); ax.set_title("Stable genes: primary only",fontweight="bold",pad=8); clean(ax)
ax.text(.98,.98,"No recurrent gene passed",transform=ax.transAxes,ha="right",va="top",fontsize=6.4,fontweight="bold",color=RED)
fig.subplots_adjust(top=.84,bottom=.25,left=.08,right=.98)
fig.suptitle("Virtual-knockout nominations fail the recurrent-network stability gate",fontsize=11,fontweight="bold",y=.98)
save(fig,"Fig5_revised")

print(OUT)
for p in sorted(OUT.glob("*.tiff")):
    print(p.name, p.stat().st_size)




