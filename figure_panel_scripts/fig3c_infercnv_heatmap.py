#!/usr/bin/env python3
"""inferCNV heatmap (cells x genomic bins) for selected cell states only."""
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
# keep text as editable TrueType (fonttype 42) rather than Type-3 outlines so all
# labels remain live, selectable text in Illustrator
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="results/infercnv_run")
ap.add_argument("--state-col", default="IDENT")
ap.add_argument("--min-genes", type=int, default=0,
                help="drop cells with fewer usable genes; 0 = no filter")
ap.add_argument("--out", default=None)
ap.add_argument("--palette", default=None,
                help="TSV/CSV with columns group,color to fix block colours "
                     "(shared with the ggtree clone palette so heatmap & tree match)")
ap.add_argument("--row-order", default=None,
                help="file with one group label per line (top->bottom, e.g. a tree "
                     "tiporder); orders the row blocks by it and removes the dendrogram")
ap.add_argument("--merge", default="",
                help="comma-separated groups collapsed into one row block")
ap.add_argument("--merge-name", default="",
                help="label for the merged block (default: joined with '_')")
args = ap.parse_args()

DIR = Path(args.dir)


def build_clone_leaf(obs):
    """Per-cell clone-tree leaf: tumour clone or Merged_Normal, matching the tree
    (Basal/DCIS-Lumen-Border dropped everywhere; Luminal DCIS dropped from the
    near-diploid Merged_Normal); excluded cells -> 'nan'."""
    cl = obs["cnv_leiden_new"].astype(str).values
    ident = obs["IDENT"].astype(str).values
    merge_set = ["merged_lt260", "6", "13", "12", "9", "10", "7", "3"]
    tumor = ["0", "1", "2", "4", "5", "8", "11"]
    EXC = ["Basal", "DCIS Lumen Border"]
    MN_EXC = ["Luminal DCIS1", "Luminal DCIS2", "Luminal DCIS3", "Luminal DCIS4"]
    g = np.full(len(cl), "nan", dtype=object)
    for c in tumor:
        g[(cl == c) & ~np.isin(ident, EXC)] = c
    g[np.isin(cl, merge_set) & ~np.isin(ident, EXC + MN_EXC)] = "Merged_Normal"
    return g.astype(str)


# cell-state grouping: Normal/Ductal Luminal merged into "Merged Normal",
# shared colour template (skyblue1 / DCIS1-4); matplotlib hex equivalents of the
# R colour names (skyblue1=#87CEFF, darkslategray3=#79CDCD)
order = ["Luminal DCIS1", "Luminal DCIS2", "Luminal DCIS3",
         "Luminal DCIS4", "Merged Normal"]
pal = {"Luminal DCIS1": "#CD2626", "Luminal DCIS2": "#DAA520",
       "Luminal DCIS3": "#79CDCD", "Luminal DCIS4": "#9ACD32",
       "Merged Normal": "#87CEFF"}

print("loading h5ad ...")
a = ad.read_h5ad(DIR / "infercnv_atera.h5ad")
# optionally remove low-quality OBSERVATION cells (< min_genes usable genes across
# chromosomes); reference cells (Normal/Ductal Luminal) are ALWAYS kept so the
# Merged-Normal baseline stays intact.
if args.min_genes > 0:
    usable = np.asarray((a.X > 0).sum(1)).ravel()
    is_ref = a.obs["IDENT"].astype(str).isin(["Normal Luminal", "Ductal Luminal"]).values
    keepcell = (usable >= args.min_genes) | is_ref
    print(f"removing {int((~keepcell).sum())} observation cells with <{args.min_genes} "
          f"usable genes (references always kept; {a.n_obs} -> {int(keepcell.sum())})")
    a = a[keepcell].copy()
X = a.obsm["X_cnv"]
X = np.asarray(X.todense()) if not isinstance(X, np.ndarray) else X
if args.state_col == "clone_leaf":
    state = build_clone_leaf(a.obs)            # tree leaves (8 clones)
else:
    state = a.obs[args.state_col].astype(str).values
    if args.state_col == "IDENT":              # merge Normal/Ductal -> Merged Normal
        state = np.where(np.isin(state, ["Normal Luminal", "Ductal Luminal"]),
                         "Merged Normal", state)
# label the reference block "Merged Normal" for any grouping (e.g. cnv_leiden_w_ref)
state = np.where(state == "REFERENCE", "Merged Normal", state)
# optionally collapse several groups into a single row block (e.g. merge clones 3+4)
if args.merge:
    mset = [c.strip() for c in args.merge.split(",") if c.strip()]
    mname = args.merge_name or "_".join(mset)
    state = np.where(np.isin(state, mset), mname, state)
    print(f"merged {mset} -> '{mname}' ({int((state == mname).sum())} cells)")
chr_pos = a.uns["cnv"]["chr_pos"]
nbin = X.shape[1]

# for a non-cell-state grouping (e.g. clones), derive the group order
# (numeric first, merged/other last) and an automatic palette
if args.state_col != "IDENT":
    import matplotlib.cm as _cm
    present = [s for s in pd.unique(state) if s != "nan"]
    order = sorted(present, key=lambda s: (1, 0) if not str(s).isdigit() else (0, int(s)))
    # palette: tab10 assigned in the row/tip order so it MATCHES spatial_clone_facets
    if args.row_order:
        want = [l.strip() for l in open(args.row_order) if l.strip()]
        pal_order = [s for s in want if s in present] + [s for s in present if s not in want]
    else:
        pal_order = order
    tab10 = _cm.get_cmap("tab10").colors
    pal = {s: tab10[i % 10] for i, s in enumerate(pal_order)}
    TITLE = f"inferCNV heatmap by {args.state_col}"
else:
    TITLE = "inferCNV heatmap — Normal/Ductal Luminal + Luminal DCIS1-4"

# optional explicit palette (group -> colour) shared with the ggtree clone palette
if args.palette:
    sep = "\t" if str(args.palette).endswith((".tsv", ".txt")) else ","
    pdf = pd.read_csv(args.palette, sep=sep)
    pdf.columns = [c.lower() for c in pdf.columns]
    pal.update({str(g): c for g, c in zip(pdf["group"], pdf["color"])})

# keep only the requested states, ordered, and within each state sort by total CNV
keep_idx = []
row_state = []
for s in order:
    idx = np.where(state == s)[0]
    if len(idx) == 0:
        continue
    burden = np.abs(X[idx]).sum(1)
    idx = idx[np.argsort(burden)]          # quiet -> altered within the group
    keep_idx.append(idx)
    row_state += [s] * len(idx)
keep_idx = np.concatenate(keep_idx)
Xs = X[keep_idx]
row_state = np.array(row_state)
print(f"{Xs.shape[0]} cells x {nbin} bins")

# chromosome layout
order_chr = sorted(chr_pos, key=lambda c: int(c[3:]))
bounds = [(c, chr_pos[c], (chr_pos[order_chr[i + 1]] if i + 1 < len(order_chr) else nbin))
          for i, c in enumerate(order_chr)]

# light smoothing along the genome (within each chromosome) for clearer bands
W = 7
sm = np.empty_like(Xs)
for (_, s, e) in bounds:
    seg = Xs[:, s:e]
    k = min(W, seg.shape[1])
    csum = np.cumsum(np.insert(seg, 0, 0, axis=1), axis=1)
    run = (csum[:, k:] - csum[:, :-k]) / k
    pad = seg.shape[1] - run.shape[1]
    lo = pad // 2
    out = np.empty_like(seg)
    out[:, lo:lo + run.shape[1]] = run
    for j in range(lo):
        out[:, j] = run[:, 0]
    for j in range(lo + run.shape[1], seg.shape[1]):
        out[:, j] = run[:, -1]
    sm[:, s:e] = out
Xs = sm

# ---- order the row blocks -----------------------------------------------------
present = [s for s in order if (row_state == s).sum() > 0]
if args.row_order:
    # follow an explicit order (e.g. a phylogenetic-tree tiporder); no dendrogram
    want = [l.strip() for l in open(args.row_order) if l.strip()]
    leaf_states = [s for s in want if s in present] + [s for s in present if s not in want]
    print("row order from tree:", leaf_states)
else:
    # cluster the TUMOUR groups (dendrogram) and PIN the normal/reference group
    # (Normal/Ductal Luminal or Merged_Normal) to the bottom as an outgroup.
    from scipy.cluster.hierarchy import linkage, dendrogram
    NORMAL_SET = {"Normal Luminal", "Ductal Luminal", "Merged_Normal", "Merged Normal",
                  "REFERENCE"}
    norm = [s for s in present if s in NORMAL_SET]
    tum = [s for s in present if s not in NORMAL_SET]
    means = np.vstack([Xs[row_state == s].mean(0) for s in tum])
    Z = linkage(means, method="ward") if len(tum) > 1 else None
    tum_order = dendrogram(Z, labels=tum, no_plot=True)["ivl"] if Z is not None else tum
    leaf_states = tum_order + norm                 # normal/reference pinned at bottom
# reorder rows so blocks follow this order (keep within-state burden order)
new_idx = np.concatenate([np.where(row_state == s)[0] for s in leaf_states])
Xs = Xs[new_idx]; row_state = row_state[new_idx]

from matplotlib.colors import TwoSlopeNorm
# symmetric diverging scale: full blue at -1, white at 0, full red at +1
cnorm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
print(f"colour scale {cnorm.vmin} .. {cnorm.vmax} (white at {cnorm.vcenter})")
fig = plt.figure(figsize=(16, 9))
draw_dendro = not args.row_order
if draw_dendro:
    gs = fig.add_gridspec(1, 4, width_ratios=[0.9, 0.10, 6, 0.12], wspace=0.02)
    ann_col, heat_col, cbar_col = 1, 2, 3
    # ---- row dendrogram (tumour groups), aligned to the block centres ---------
    axd = fig.add_subplot(gs[0, 0])
    acc = 0; center = {}
    for s in leaf_states:
        nrows = int((row_state == s).sum()); center[s] = acc + nrows / 2.0; acc += nrows
    nrow = acc
    node_y = {i: center[tum[i]] for i in range(len(tum))}
    node_x = {i: 0.0 for i in range(len(tum))}
    if Z is not None:
        for i, (aa, bb, dist, _) in enumerate(Z):
            aa, bb, nid = int(aa), int(bb), len(tum) + i
            node_y[nid] = (node_y[aa] + node_y[bb]) / 2.0; node_x[nid] = dist
            axd.plot([node_x[aa], dist], [node_y[aa], node_y[aa]], color="0.3", lw=0.9)
            axd.plot([node_x[bb], dist], [node_y[bb], node_y[bb]], color="0.3", lw=0.9)
            axd.plot([dist, dist], [node_y[aa], node_y[bb]], color="0.3", lw=0.9)
        axd.set_xlim(float(Z[:, 2].max()) * 1.05, 0)   # root at left, leaves at right
    axd.set_ylim(nrow, 0)
    axd.axis("off")
else:
    gs = fig.add_gridspec(1, 3, width_ratios=[0.10, 6, 0.12], wspace=0.02)
    ann_col, heat_col, cbar_col = 0, 1, 2

# annotation bar (cell state / clone) -- drawn as VECTOR rectangles (one per block)
# so the colours are exact and editable in Illustrator. (imshow embeds a raster
# whose indexed colours get colour-managed differently by Illustrator vs Preview.)
from matplotlib.patches import Rectangle
axann = fig.add_subplot(gs[0, ann_col])
acc = 0; blk_centers = []
for s in leaf_states:
    n = int((row_state == s).sum())
    axann.add_patch(Rectangle((0, acc), 1, n, facecolor=pal[s], edgecolor="none"))
    blk_centers.append(acc + n / 2.0); acc += n
axann.set_xlim(0, 1); axann.set_ylim(acc, 0)       # 0 at top, matches the heatmap
axann.set_xticks([])
blk_labels = [f"clone {s}" if str(s).isdigit() else str(s) for s in leaf_states]
axann.set_yticks(blk_centers)
axann.set_yticklabels(blk_labels, fontsize=8)
axann.tick_params(axis="y", left=True, labelleft=True, length=0, pad=2)
for t in axann.get_yticklabels():
    t.set_ha("right")

# main heatmap
axh = fig.add_subplot(gs[0, heat_col])
im = axh.imshow(Xs, aspect="auto", cmap="bwr", norm=cnorm,
                interpolation="nearest", rasterized=True)
for (_, s, e) in bounds:
    axh.axvline(e - 0.5, color="black", lw=0.4)
# group separators between states (rows are in dendrogram order)
acc = 0
for s in leaf_states:
    acc += int((row_state == s).sum())
    axh.axhline(acc - 0.5, color="black", lw=0.6)
centers = [(s + e) / 2 for (_, s, e) in bounds]
axh.set_xticks(centers)
axh.set_xticklabels([c.replace("chr", "") for (c, _, _) in bounds], fontsize=8)
axh.set_yticks([]); axh.set_xlabel("genomic position (chromosome)")
axh.set_title(TITLE)

# colorbar
axc = fig.add_subplot(gs[0, cbar_col])
fig.colorbar(im, cax=axc, label="inferCNV signal")

handles = [Patch(color=pal[s], label=f"{s} (n={int((row_state==s).sum())})") for s in order
           if (row_state == s).sum() > 0]
axh.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.06),
           ncol=3, fontsize=8, frameon=False)

base = "luminal_states" if args.state_col == "IDENT" else args.state_col
suffix = f"_min{args.min_genes}genes" if args.min_genes > 0 else ""
out = args.out or str(DIR / f"heatmap_{base}{suffix}.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("saved", out)
