#!/usr/bin/env python3
"""Fig 3h - per-region genome copy-number gain/loss tracks (16 WGS samples)."""

import os, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["pdf.fonttype"]=42   # editable text in Illustrator
matplotlib.rcParams["ps.fonttype"]=42
matplotlib.rcParams["svg.fonttype"]="none"
from scipy.cluster.hierarchy import linkage, leaves_list

CNS="/path/to/atera_wgs/results/cnvkit/cohort/all"
FAI="/path/to/refdata/human/GRCh38/hg38.fa.fai"
OUT="/path/to/atera_wgs/results/cnvkit_summary"
PS=f"{OUT}/per_sample_gainloss"; os.makedirs(PS,exist_ok=True)
BIN=500_000
CHROMS=[f"chr{i}" for i in [1,2,4,5,8,11,12,17]]
YL=2.0  # fixed symmetric y-limit (log2), matches clip

sizes={l.split("\t")[0]:int(l.split("\t")[1]) for l in open(FAI)}
off={}; nb={}; tot=0
for c in CHROMS:
    n=int(np.ceil(sizes[c]/BIN)); nb[c]=n; off[c]=tot; tot+=n
x=np.arange(tot)

def key(f):
    m=re.search(r"cap_(\d+)_DL",f); return int(m.group(1)) if m else 999
files=[f for f in sorted(glob.glob(f"{CNS}/*.cns"),key=key) if not re.search(r"(call|bintest)",f)]
names=[os.path.basename(f).replace(".cns","") for f in files]

def bin_sample(df):
    v=np.full(tot,np.nan)
    for c in CHROMS:
        sub=df[df.chromosome==c]
        if sub.empty: continue
        mids=(np.arange(nb[c])+0.5)*BIN; st=sub.start.values; en=sub.end.values; lg=sub.log2.values
        idx=np.clip(np.searchsorted(st,mids,side="right")-1,0,len(st)-1)
        v[off[c]:off[c]+nb[c]]=np.where((mids>=st[idx])&(mids<en[idx]),lg[idx],np.nan)
    return np.clip(v,-YL,YL)

vecs={s:bin_sample(pd.read_csv(f,sep="\t")) for s,f in zip(names,files)}

# ---- order rows to match the WGS CNV phylogeny (ggtree tip order, top->bottom) ----
TIP="/path/to/atera_wgs/results/biological_insight/wgs_cnv_phylo_tree_ggtree_tiporder.txt"
tiporder=[l.strip() for l in open(TIP) if l.strip()]
names=[s for s in tiporder if s in names]+[s for s in names if s not in tiporder]

def draw(ax, v, label, xlabels=False, yticks=True):
    for i,c in enumerate(CHROMS):
        if i%2==1: ax.axvspan(off[c]-0.5, off[c]+nb[c]-0.5, color="0.93", zorder=0)
    ax.axhline(0,color="black",lw=0.5,zorder=4)
    ax.fill_between(x, 0, np.clip(v,0,None), step="mid", color="#c0392b", linewidth=0, zorder=3)  # gains (vector)
    ax.fill_between(x, 0, np.clip(v,None,0), step="mid", color="#2166ac", linewidth=0, zorder=3)  # losses (vector)
    for c in CHROMS[1:]: ax.axvline(off[c]-0.5, color="0.6", lw=0.4, zorder=1)
    ax.set_ylim(-YL,YL); ax.set_xlim(-0.5,tot-0.5)
    if yticks:
        ax.set_yticks([-2,0,2]); ax.tick_params(labelsize=7)
    else:
        ax.set_yticks([]); ax.tick_params(left=False)
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8)
    if xlabels:
        ax.set_xticks([off[c]+nb[c]/2 for c in CHROMS]); ax.set_xticklabels([c[3:] for c in CHROMS],fontsize=8)
    else:
        ax.set_xticks([])

# ---- faceted overview: one panel per sample ----
fig,axes=plt.subplots(len(names),1,figsize=(5.0,(1.05*len(names)+1)/3),sharex=True,
                      gridspec_kw=dict(hspace=0.3))
for i,(s,ax) in enumerate(zip(names,axes)):
    draw(ax, vecs[s], s, xlabels=(i==len(names)-1), yticks=False)
axes[0].set_title("Per-sample CNV (segmented log2, 0.5 Mb bins) — chr "+", ".join(c[3:] for c in CHROMS)
                  +"   red=gain  blue=loss   (clipped ±2)", fontsize=8)
axes[-1].set_xlabel("chromosome")
plt.savefig(f"{OUT}/cnvkit_gainloss_per_sample_faceted.png",dpi=150,bbox_inches="tight")
plt.savefig(f"{OUT}/cnvkit_gainloss_per_sample_faceted.pdf",dpi=300,bbox_inches="tight")
plt.close()

# ---- individual per-sample figures ----
for s in names:
    fig,ax=plt.subplots(figsize=(15,2.6))
    draw(ax, vecs[s], "log2", xlabels=True)
    for i,c in enumerate(CHROMS):
        yl=YL*0.85 if i%2==0 else -YL*0.85
        ax.text(off[c]+nb[c]/2, yl, c.replace("chr",""), ha="center", va="center", fontsize=8, color="0.3")
    ax.set_ylabel("log2 ratio\n(gain ↑ / loss ↓)")
    ax.set_title(f"{s} — CNV (segmented log2, 0.5 Mb bins)   red=gain  blue=loss")
    plt.tight_layout()
    plt.savefig(f"{PS}/{s}_gainloss.png",dpi=150)
    plt.savefig(f"{PS}/{s}_gainloss.pdf",dpi=300)
    plt.close()

print("samples:",len(names),"| bins:",tot)
print("faceted ->",f"{OUT}/cnvkit_gainloss_per_sample_faceted.pdf")
print("individual ->",PS)
