#!/usr/bin/env python3
"""mregDC / pDC distance-to-nearest-clone-cell vs Merged_Normal (boxplot + beeswarm, MWU/BH stats)."""
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
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

ap = argparse.ArgumentParser()
ap.add_argument("--loc", default="cell_annotations.csv")
ap.add_argument("--dir", default="results/infercnv_run")
ap.add_argument("--clone-h5ad",
                default="results/infercnv_run/lum_state_anchors.h5ad")
ap.add_argument("--clone-col", default="lum_state")
ap.add_argument("--niche-col", default="Cell_states")
ap.add_argument("--types", default="mregDC,pDCs")
ap.add_argument("--all-types", action="store_true",
                help="use every TME Cell_state (grid layout) instead of --types")
ap.add_argument("--ncol", type=int, default=5, help="panels per row when many types")
ap.add_argument("--panel-w", type=float, default=4.2, help="width per panel (inches)")
ap.add_argument("--panel-h", type=float, default=5.2, help="height per panel row (inches)")
ap.add_argument("--min-anchor", type=int, default=15,
                help="skip a cell type with fewer than this many cells (unstable stats)")
ap.add_argument("--ref", default="Merged_Normal", help="reference clone every other clone is tested against")
ap.add_argument("--amp-clones", default="Luminal DCIS4,Luminal DCIS3")
ap.add_argument("--clone-merge", default="",
                help="merge clone labels, e.g. '3,4:Clone_3-4' (';'-separate multiple groups)")
ap.add_argument("--clone-order", default="",
                help="explicit left-to-right x-axis clone order (comma-separated)")
ap.add_argument("--swarm-max", type=int, default=300)
ap.add_argument("--nperm", type=int, default=5000,
                help="label-permutation reps for the clone-vs-ref null (accounts for group size)")
ap.add_argument("--min-effect", type=float, default=0.2,
                help="min |rank-biserial| (Cliff's delta) to call a comparison meaningful")
ap.add_argument("--no-effect-gate", action="store_true",
                help="star purely on permutation q (ignore the effect-size cutoff); "
                     "writes a separate '_qonly' figure")
ap.add_argument("--reverse", action="store_true",
                help="reverse direction: for each DC cell measure distance to the nearest "
                     "cell of each clone (dots = DC cells; paired stats across clones)")
ap.add_argument("--log-y", action="store_true")
ap.add_argument("--tag", default="lumstate_fine")
args = ap.parse_args()

DIR = Path(args.dir); outdir = DIR / "biological_insight"; outdir.mkdir(exist_ok=True)
TAG = ("_" + args.tag) if args.tag else ""
rng = np.random.default_rng(0)

LUMSTATE_COL = {
    "Luminal DCIS1": "#87CEFF",   # skyblue1
    "Luminal DCIS2": "#CD2626",
    "Luminal DCIS3": "#DAA520",
    "Luminal DCIS4": "#79CDCD",   # darkslategray3
    "Merged_Normal": "#9ACD32",
    "Normal Luminal": "#9ACD32",
    # CNV clone_grp groups -- tab10 colours matching the inferCNV heatmap legend
    "Clone_3-4": "#ff7f0e",       # orange
    "0_1_5_9":   "#2ca02c",       # green
    "8":         "#d62728",       # red
    "2":         "#9467bd",       # purple
}

# --------------------------------------------------------------------------- #
# cells + clone labels + niche annotation
# --------------------------------------------------------------------------- #
loc = pd.read_csv(args.loc); loc["cell_id"] = loc["cell_id"].astype(str)
loc[args.niche_col] = loc[args.niche_col].astype(str)

import anndata as ad
a = ad.read_h5ad(args.clone_h5ad, backed="r")
cl = pd.DataFrame({"cell_id": a.obs_names.astype(str),
                   "clone": a.obs[args.clone_col].astype(str).values})
cl = cl[~cl["clone"].isin(["nan", "NaN", "None"])]
# optional clone relabelling, e.g. "3,4:Clone_3-4" merges clones 3 and 4
for grp in [g for g in args.clone_merge.split(";") if ":" in g]:
    members, name = grp.split(":")
    cl.loc[cl["clone"].isin([m.strip() for m in members.split(",")]), "clone"] = name.strip()
clones_df = loc.merge(cl, on="cell_id", how="inner").dropna(subset=["x", "y", "clone"])

# canonical TME cell-type order (mirrors clone_niche_distance_boxplot.py)
ALL_TYPES = ["FBs", "Endothelial Cell", "Pericyte", "Smooth Muscle Cell",
             "FOLR2 TRMs", "CXCL9 TAMs", "SPP1 TAMs", "Classical MONOs",
             "Non-classical MONOs", "cDC1s", "cDC2s", "mregDC", "pDCs",
             "CD4T_Res", "CD8T_Eff", "TCF7_CD4T", "Tregs", "ISG_Tcells",
             "prolif_Tcells", "NK", "Bcells", "PCs"]
ann_arr = loc[args.niche_col].to_numpy()
XY_all = loc[["x", "y"]].to_numpy()
if args.all_types:
    types = [t for t in ALL_TYPES if (ann_arr == t).sum() >= args.min_anchor]
else:
    types = [t.strip() for t in args.types.split(",") if t.strip()]
    types = [t for t in types if (ann_arr == t).sum() >= args.min_anchor]

# x-axis order: explicit --clone-order, else ref first then sorted
present = set(clones_df["clone"].unique())
if args.clone_order:
    order = [c.strip() for c in args.clone_order.split(",") if c.strip() in present]
    clones = order + sorted(present - set(order))
else:
    clones = sorted(present, key=lambda c: (0, 0) if c == args.ref else (1, c))
amp = {c.strip() for c in args.amp_clones.split(",") if c.strip()}
clone_col = {c: LUMSTATE_COL.get(c, "#8C8C8C") for c in clones}

# --------------------------------------------------------------------------- #
# distance, both directions
#   forward  : anchors = clone cells   -> distance to nearest DC of type t
#              (dots = clone cells; each clone an INDEPENDENT sample)
#   reverse  : anchors = DC cells of t -> distance to nearest cell of each clone
#              (dots = DC cells; the SAME DC cells scored against every clone,
#               so clone-vs-ref comparisons are PAIRED)
# --------------------------------------------------------------------------- #
clone_xy = {c: clones_df.loc[clones_df["clone"] == c, ["x", "y"]].to_numpy()
            for c in clones}
panel_data = []   # panel_data[ti] = list over clones of distance arrays
for t in types:
    if args.reverse:
        dc_xy = XY_all[ann_arr == t]
        cols = [cKDTree(clone_xy[c]).query(dc_xy, k=1)[0] for c in clones]
    else:
        anc = clones_df[["x", "y"]].to_numpy()
        ccode = clones_df["clone"].map({c: i for i, c in enumerate(clones)}).to_numpy()
        dt = cKDTree(XY_all[ann_arr == t]).query(anc, k=1)[0]
        cols = [dt[ccode == i] for i in range(len(clones))]
    panel_data.append(cols)
PAIRED = args.reverse

def stars(p):
    return "n.s." if p >= 0.05 else "*" if p >= 1e-2 else "**" if p >= 1e-3 else \
           "***" if p >= 1e-4 else "****"

def compare(cv, refv):
    """Return (effect, perm_p). effect>0 = anchor sits CLOSER to `cv`'s clone
    than to the reference. Paired (sign-flip) or independent (label-perm)."""
    if PAIRED:                                   # same anchor cells -> paired
        diff = cv - refv
        eff = (np.sum(cv < refv) - np.sum(cv > refv)) / len(cv)   # matched rank-biserial
        obs = abs(diff.mean()); n = len(diff); ge = 0
        for _ in range(args.nperm):
            signs = rng.integers(0, 2, n) * 2 - 1
            if abs((diff * signs).mean()) >= obs - 1e-9:
                ge += 1
        return eff, (ge + 1) / (args.nperm + 1)
    U, _ = mannwhitneyu(cv, refv, alternative="two-sided")
    eff = 1.0 - 2.0 * U / (len(cv) * len(refv))
    pooled = np.concatenate([cv, refv]); N = len(pooled); n1 = len(cv)
    tot = pooled.sum(); obs = abs(cv.mean() - refv.mean()); ge = 0
    for _ in range(args.nperm):
        s = pooled[rng.permutation(N)[:n1]].sum()
        if abs(s / n1 - (tot - s) / (N - n1)) >= obs - 1e-9:
            ge += 1
    return eff, (ge + 1) / (args.nperm + 1)

# --------------------------------------------------------------------------- #
# figure
# --------------------------------------------------------------------------- #
ncol = min(args.ncol, len(types))
nrow = int(np.ceil(len(types) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(args.panel_w * ncol, args.panel_h * nrow),
                         squeeze=False)
stat_rows = []
ref_i = clones.index(args.ref)
others = [c for c in clones if c != args.ref]

for ti, t in enumerate(types):
    ax = axes[ti // ncol][ti % ncol]
    data = panel_data[ti]
    bp = ax.boxplot(data, positions=range(len(clones)), widths=0.6,
                    showfliers=False, patch_artist=True, zorder=2,
                    medianprops=dict(color="black", lw=1.4),
                    whiskerprops=dict(color="#444"), capprops=dict(color="#444"))
    for patch, c in zip(bp["boxes"], clones):
        patch.set_facecolor(clone_col[c]); patch.set_alpha(0.35)
        patch.set_edgecolor("#222"); patch.set_linewidth(2.4 if c in amp else 1.0)
    for xi, c in enumerate(clones):
        v = data[xi]
        vv = rng.choice(v, args.swarm_max, replace=False) if len(v) > args.swarm_max else v
        jit = (rng.random(len(vv)) - 0.5) * 0.5
        ax.scatter(xi + jit, vv, s=5, color=clone_col[c], alpha=0.5,
                   linewidths=0, zorder=3)

    # Each clone vs reference (effect>0 = anchor closer to that clone than to ref).
    ref_v = data[ref_i]
    pvals, effs = {}, {}
    for c in others:
        effs[c], pvals[c] = compare(data[clones.index(c)], ref_v)
    ps = np.array([pvals[c] for c in others]); o = ps.argsort()
    q = np.empty_like(ps)
    q[o] = np.clip(np.minimum.accumulate((ps[o] * len(ps) / (np.arange(len(ps)) + 1))[::-1])[::-1], 0, 1)
    qvals = dict(zip(others, q))
    for c in others:
        passed = (qvals[c] < 0.05) and (args.no_effect_gate or abs(effs[c]) >= args.min_effect)
        stat_rows.append({"niche_type": t, "clone": c, "ref": args.ref,
                          "clone_median_dist": float(np.median(data[clones.index(c)])),
                          "ref_median_dist": float(np.median(ref_v)),
                          "rank_biserial": effs[c], "perm_p": pvals[c], "BH_q": qvals[c],
                          "significant": passed})

    # brackets: show effect size; star only if q<0.05 AND |effect|>=min-effect
    ymax = max(np.percentile(d, 98) for d in data)
    ybase = ymax * 1.02
    step = ymax * 0.135
    for k, c in enumerate(others):
        j = clones.index(c)
        y = ybase + k * step
        x1, x2 = ref_i, j
        passed = (qvals[c] < 0.05) and (args.no_effect_gate or abs(effs[c]) >= args.min_effect)
        col = "#333" if passed else "#B0B0B0"
        ax.plot([x1, x1, x2, x2], [y, y + step * 0.14, y + step * 0.14, y],
                lw=1.0, color=col)
        lab = stars(qvals[c]) if passed else "ns"
        ax.text((x1 + x2) / 2, y + step * 0.16, lab,
                ha="center", va="bottom", fontsize=9, color=col)
    ax.set_ylim(0, ybase + len(others) * step + step)

    ax.set_xticks(range(len(clones)))
    ax.set_xticklabels([c.replace("Luminal ", "").replace("Merged_", "")
                        for c in clones], rotation=45, ha="right", fontsize=9)
    if args.log_y:
        ax.set_yscale("log")
    n_anchor = len(data[0]) if PAIRED else None   # reverse: same anchor cells across clones
    ax.set_title(f"{t} (n={n_anchor})" if PAIRED else t,
                 fontsize=11, fontweight="bold")
    if ti % ncol == 0:
        ax.set_ylabel("cell-type distance to nearest clone cell (µm)" if PAIRED
                      else "distance to nearest (µm)", fontsize=9)

for j in range(len(types), nrow * ncol):        # hide unused grid cells
    axes[j // ncol][j % ncol].axis("off")

if args.all_types:
    STEM = "nichedist_to_clone" if args.reverse else "clonedist_niche_vs_normal"
else:
    STEM = "DCdist_to_clone" if args.reverse else "clonedist_DC_vs_normal"
SUF = "_qonly" if args.no_effect_gate else ""
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(outdir / f"06{TAG}_{STEM}{SUF}.{ext}", dpi=220, bbox_inches="tight")
pd.DataFrame(stat_rows).round(5).to_csv(
    outdir / f"06{TAG}_{STEM}{SUF}_stats.tsv", sep="\t", index=False)
print("saved", outdir / f"06{TAG}_{STEM}{SUF}.png")
print(pd.DataFrame(stat_rows).round(4).to_string(index=False))
