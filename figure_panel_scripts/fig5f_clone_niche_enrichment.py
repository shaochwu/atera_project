#!/usr/bin/env python3
"""Fig 5f - CNV clonal phylogeny + per-clone immune/stromal fold-enrichment."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
# keep PDF/PS text as editable TrueType (Type 42) glyphs, not Type 3 outlines,
# so every label stays selectable/editable in Illustrator
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.spatial import cKDTree
import anndata as ad

ap = argparse.ArgumentParser()
ap.add_argument("--loc", default="cell_annotations.csv")
ap.add_argument("--niche-csv", default="/path/to/data/atera_cell_assignment_to_niches_w_cones.csv")
ap.add_argument("--clone-h5ad",
                default="results/infercnv_run/infercnv_atera.h5ad")
ap.add_argument("--clone-col", default="clone_grp")
ap.add_argument("--radius", type=float, default=50.0)
ap.add_argument("--metric", choices=["enrichment", "share"], default="enrichment",
                help="enrichment = fold obs/exp (size-normalized); share = %% of population")
ap.add_argument("--orient", choices=["vertical", "horizontal"], default="horizontal",
                help="horizontal (published) = tree left, populations side-by-side; vertical = tree on top, stacked")
args = ap.parse_args()

D = Path("results/infercnv_run/biological_insight")

# ---- load cells, niche, clone labels -------------------------------------- #
loc = pd.read_csv(args.loc); loc["cell_id"] = loc["cell_id"].astype(str)
loc = loc.dropna(subset=["x", "y"]).reset_index(drop=True)
nm = pd.read_csv(args.niche_csv)[["cell_id", "niche"]]; nm["cell_id"] = nm["cell_id"].astype(str)
loc = loc.merge(nm, on="cell_id", how="left")
a = ad.read_h5ad(args.clone_h5ad, backed="r")
cg = pd.DataFrame({"cell_id": a.obs_names.astype(str), "clone": a.obs[args.clone_col].astype(str).values})
cg = cg[~cg["clone"].isin(["nan", "NaN", "None"])]
cg.loc[cg["clone"].isin(["3", "4"]), "clone"] = "Clone_3-4"
epi = loc.merge(cg, on="cell_id", how="inner")

# ---- propagate clone labels to ALL tumour cells, then assign each population --
#      cell to its nearest clone; share = (cells nearest clone c) / (all cells of
#      that population). Shares partition each population across clones (~100%). --
ref = epi[epi["clone"] != "Merged_Normal"]
is_tum = loc["Cell_states"].astype(str).str.startswith("Luminal DCIS")
tum = loc[is_tum].copy()
_, ii = cKDTree(ref[["x", "y"]].to_numpy()).query(tum[["x", "y"]].to_numpy(), k=1)
tum["clone"] = ref["clone"].to_numpy()[ii]
nrm = loc[loc["Cell_states"].astype(str).isin(["Normal Luminal", "Ductal Luminal"])].copy()
nrm["clone"] = "Merged_Normal"
epi = pd.concat([tum, nrm], ignore_index=True)
order = ["Clone_3-4", "0_1_5_9", "8", "2", "Merged_Normal"]
POPS = ["FOLR2 TRMs", "SPP1 TAMs", "niche 4", "mregDC", "TCF7_CD4T"]
etree = cKDTree(epi[["x", "y"]].to_numpy()); eclone = epi["clone"].to_numpy()
epi_share = epi["clone"].value_counts(normalize=True).reindex(order).fillna(0) * 100  # expected
masks = {pop: (loc["niche"] == 4) if pop == "niche 4"
              else (loc["Cell_states"] == pop) for pop in POPS}
n_clone = epi["clone"].value_counts().reindex(order).fillna(0)   # clone cell counts
tot_epi = len(epi)
prop = {c: {} for c in order}; enr = {c: {} for c in order}; npop = {}
for pop, m in masks.items():
    sub = loc[m].dropna(subset=["x", "y"]); npop[pop] = len(sub)
    _, j = etree.query(sub[["x", "y"]].to_numpy(), k=1)
    vc = pd.Series(eclone[j]).value_counts()
    for c in order:
        obs = vc.get(c, 0)
        prop[c][pop] = 100 * obs / len(sub)
        exp = len(sub) * n_clone[c] / tot_epi                   # expected from clone size
        enr[c][pop] = (obs + 0.5) / (exp + 0.5)                 # fold-enrichment (pseudocount)
P = pd.DataFrame(prop).T[POPS]; E = pd.DataFrame(enr).T[POPS]
P.round(2).to_csv(D / "fig5f_clone_niche_enrichment_share.tsv", sep="\t")
E.round(3).to_csv(D / "fig5f_clone_niche_enrichment_enrichment.tsv", sep="\t")
print("population totals:", npop)
print("SHARE %:\n", P.round(1).to_string())
print("FOLD-ENRICHMENT (obs/exp):\n", E.round(2).to_string())

# ---- tree geometry (from cnv_phylo_tree_mp_clones.newick) ------------------ #
tip_y = {"Clone_3-4": 4, "0_1_5_9": 3, "8": 2, "2": 1, "Merged_Normal": 0}
X = {"root": 0.0, "A": 2.25, "B": 2.70, "C": 4.05, "Merged_Normal": 0.0,
     "2": 4.05, "8": 2.70, "0_1_5_9": 4.50, "Clone_3-4": 4.50}
y = dict(tip_y)
y["C"] = (y["0_1_5_9"] + y["Clone_3-4"]) / 2
y["B"] = (y["C"] + y["8"]) / 2
y["A"] = (y["B"] + y["2"]) / 2
y["root"] = (y["A"] + y["Merged_Normal"]) / 2
edges = [("root", "A"), ("root", "Merged_Normal"), ("A", "B"), ("A", "2"),
         ("B", "C"), ("B", "8"), ("C", "0_1_5_9"), ("C", "Clone_3-4")]
CLONE_COL = {"2": "#9467bd", "8": "#d62728", "Merged_Normal": "#8c564b",
             "0_1_5_9": "#2ca02c", "Clone_3-4": "#ff7f0e"}
SHORT = {"Merged_Normal": "Normal", "2": "2", "8": "8", "0_1_5_9": "0/1/5/9", "Clone_3-4": "3-4"}
AMP = {"Clone_3-4": 1, "8": 1, "0_1_5_9": 1, "2": 0, "Merged_Normal": 0}

# ---- metric-driven values/formatting (shared by both orientations) --------- #
np_ = len(POPS)
if args.metric == "enrichment":
    V = E; fmt = lambda v: f"{v:.1f}×"; unit = "fold vs clone size"
    lims = {p: E.values.max() * 1.18 for p in POPS}             # shared scale, tight headroom
    ptitle = lambda p: f"{p}\n(n={npop[p]})"
    suptxt = "Enrichment per clone (fold vs clone size)"
else:
    V = P; fmt = lambda v: f"{v:.0f}%"; unit = "% of population"
    lims = {p: P[p].max() * 1.32 for p in POPS}
    ptitle = lambda p: f"% of all {p}\n(n={npop[p]})"
    suptxt = "Share of each population by clone"

# ---- VERTICAL: populations stacked on top, tree at the bottom -------------- #
if args.orient == "vertical":
    fig = plt.figure(figsize=(3.7, 1.2 + 0.72 * np_))
    gs = GridSpec(1 + np_, 1, height_ratios=[1.0] * np_ + [1.7], hspace=0.08)
    axB = [fig.add_subplot(gs[k]) for k in range(np_)]
    axT = fig.add_subplot(gs[np_])
    for k, pop in enumerate(POPS):
        ax = axB[k]
        for c in order:
            v = V.loc[c, pop]
            ax.bar(y[c], v, width=0.85, color=CLONE_COL[c], edgecolor="#222", linewidth=1.0, zorder=2)
            ax.text(y[c], v, fmt(v), ha="center", va="bottom", fontsize=7, zorder=3)
        ax.set_xlim(-0.6, 4.6); ax.set_ylim(0, lims[pop]); ax.set_xticks([])
        ax.set_ylabel(ptitle(pop), fontsize=8, rotation=0, ha="right", va="center", labelpad=2)
        ax.tick_params(axis="y", labelsize=6.5)
        ax.spines[["top", "right"]].set_visible(False)
    # transposed tree at bottom: x = tip order (= old y), y (up) = evolutionary distance; tips point UP
    for p, c in edges:
        axT.plot([y[p], y[c]], [X[p], X[p]], color="#333", lw=1.6)      # crossbar at parent depth
        axT.plot([y[c], y[c]], [X[p], X[c]], color="#333", lw=1.6)      # branch to child
    for t in tip_y:
        axT.scatter(y[t], X[t], s=95, color=CLONE_COL[t], edgecolor="#222", zorder=3)
    axT.set_xlim(-0.6, 4.6); axT.set_ylim(-0.3, 5.0)   # tips (large dist) at top, root at bottom
    axT.set_ylabel("evol.\ndistance", fontsize=8)
    axT.set_xticks([y[c] for c in order])
    axT.set_xticklabels([SHORT[c] for c in order], rotation=45, ha="right", fontsize=8)
    axT.spines[["top", "right", "bottom"]].set_visible(False)
    fig.suptitle(suptxt, fontsize=10, y=0.99)

# ---- HORIZONTAL: tree left, populations side-by-side ----------------------- #
else:
    fig = plt.figure(figsize=(2.9 + 0.85 * np_, 2.9))
    gs = GridSpec(1, 1 + np_, width_ratios=[2.1] + [0.9] * np_, wspace=0.08)
    axT = fig.add_subplot(gs[0]); axB = [fig.add_subplot(gs[1 + k]) for k in range(np_)]
    for p, c in edges:
        axT.plot([X[p], X[p]], [y[p], y[c]], color="#333", lw=1.6)
        axT.plot([X[p], X[c]], [y[c], y[c]], color="#333", lw=1.6)
    for t, yy in tip_y.items():
        axT.scatter(X[t], yy, s=100, color=CLONE_COL[t], edgecolor="#222", zorder=3)
        axT.text(X[t] + 0.18, yy, SHORT[t], va="center", fontsize=8, fontweight="bold")
    axT.set_xlim(-0.2, 6.0); axT.set_ylim(-0.6, 4.6)
    axT.set_title("CNV clonal phylogeny", fontsize=10)
    axT.set_xlabel("evolutionary distance", fontsize=8)
    axT.spines[["top", "right", "left"]].set_visible(False); axT.set_yticks([])
    xmax = max(lims.values())
    for k, pop in enumerate(POPS):
        ax = axB[k]
        for c in order:
            v = V.loc[c, pop]
            ax.barh(tip_y[c], v, color=CLONE_COL[c], edgecolor="#222", height=0.88, linewidth=1.0, zorder=2)
            ax.text(v, tip_y[c], " " + fmt(v), va="center", ha="left", fontsize=7, zorder=3)
        ax.set_ylim(-0.6, 4.6); ax.set_yticks([]); ax.set_xlim(0, xmax)
        ax.set_title(ptitle(pop), fontsize=8.5)
        ax.spines[["top", "right"]].set_visible(False)
        if k == 0:                       # x ticks/label on the first panel only (shared scale)
            ax.set_xlabel(unit, fontsize=8)
        else:
            ax.set_xticklabels([])
    fig.suptitle(suptxt, fontsize=10.5, y=1.13)

for ext in ("png", "pdf"):
    fig.savefig(D / f"fig5f_clone_niche_enrichment.{ext}", dpi=220, bbox_inches="tight")
print("saved", D / "fig5f_clone_niche_enrichment.png")
