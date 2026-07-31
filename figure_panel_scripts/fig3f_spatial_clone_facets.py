#!/usr/bin/env python3
"""Spatial small-multiples: one panel per CNV clone."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="results/infercnv_run")
ap.add_argument("--loc", default="cell_locations.csv")
ap.add_argument("--order", default="", help="comma-sep clone order (else tiporder/auto)")
ap.add_argument("--ncol", type=int, default=1, help="number of panel columns")
ap.add_argument("--merge", default="", help="comma-sep clones collapsed into one")
ap.add_argument("--merge-name", default="", help="label for the merged clone")
ap.add_argument("--palette", default="",
                help="TSV/CSV group,color to fix clone colours (shared with the heatmap)")
args = ap.parse_args()

DIR = Path(args.dir); outdir = DIR / "biological_insight"

a = ad.read_h5ad(DIR / "infercnv_atera.h5ad", backed="r")
clone = a.obs["clone_grp"].astype(str)
cl = pd.DataFrame({"cell_id": a.obs_names.astype(str), "clone": clone.values})
# optionally collapse several clones into one (e.g. merge clones 3+4)
if args.merge:
    mset = [x.strip() for x in args.merge.split(",") if x.strip()]
    mname = args.merge_name or "_".join(mset)
    cl["clone"] = cl["clone"].where(~cl["clone"].isin(mset), mname)
    print(f"merged {mset} -> '{mname}'")
cl = cl[cl["clone"] != "nan"]

loc = pd.read_csv(args.loc)[["cell_id", "x", "y"]]
loc["cell_id"] = loc["cell_id"].astype(str)
BG = loc.dropna(subset=["x", "y"])[["x", "y"]].values     # full section -> tissue context
df = loc.merge(cl, on="cell_id", how="inner").dropna(subset=["x", "y", "clone"])
print(f"{len(df)} clone-labelled (>=3000-gene) cells of {len(BG)} total section cells "
      f"({len(df)/len(BG)*100:.0f}%); {df['clone'].nunique()} clones")

# clone order: explicit, else tree tiporder, else size
if args.order:
    order = [c.strip() for c in args.order.split(",") if c.strip()]
else:
    tip = outdir / "cnv_phylo_tree_ggtree_mp_clones_tiporder.txt"
    if tip.exists():
        order = [l.strip() for l in open(tip) if l.strip() and l.strip() != "Merged_Normal"]
        order += ["Merged_Normal"]
    else:
        order = list(df["clone"].value_counts().index)
order = [c for c in order if c in set(df["clone"])]

XY = df[["x", "y"]].values
clv = df["clone"].values
pal = {c: col for c, col in zip(order, plt.get_cmap("tab10").colors)}
# optional fixed palette (group,color) shared with the heatmap so figures match
if args.palette:
    sep = "\t" if str(args.palette).endswith((".tsv", ".txt")) else ","
    pdf = pd.read_csv(args.palette, sep=sep)
    pdf.columns = [c.lower() for c in pdf.columns]
    pal.update({str(g): col for g, col in zip(pdf["group"], pdf["color"])})

ncol = args.ncol
nrow = int(np.ceil((len(order) + 1) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 3.3 * nrow))
axes = axes.ravel()
for ax in axes:
    ax.axis("off")

# panel 0: all clones together
ax = axes[0]; ax.axis("on"); ax.set_xticks([]); ax.set_yticks([])
ax.scatter(BG[:, 0], BG[:, 1], s=0.4, c="#ECECEC", linewidths=0, rasterized=True)
for c in order:
    s = clv == c
    ax.scatter(XY[s, 0], XY[s, 1], s=1.4, c=[pal[c]], linewidths=0, alpha=1.0,
               rasterized=True,
               label=f"{'clone ' if c[0].isdigit() else ''}{c.replace('_',' ')}")
ax.set_title("all clones", fontsize=10)
ax.invert_yaxis(); ax.set_aspect("equal")
ax.legend(markerscale=8, fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)

# one highlighted clone per panel
for j, c in enumerate(order, start=1):
    ax = axes[j]; ax.axis("on"); ax.set_xticks([]); ax.set_yticks([])
    s = clv == c
    ax.scatter(BG[:, 0], BG[:, 1], s=0.4, c="#ECECEC", linewidths=0, rasterized=True)
    ax.scatter(XY[s, 0], XY[s, 1], s=1.6, c=[pal[c]], linewidths=0, alpha=1.0,
               rasterized=True)
    lab = f"{'clone ' if c[0].isdigit() else ''}{c.replace('_', ' ')}"
    ax.set_title(f"{lab}  (n={int(s.sum()):,})", fontsize=10)
    ax.invert_yaxis(); ax.set_aspect("equal")

fig.suptitle("Spatial distribution of CNV clones (≥ 3000 genes) — "
             "grey = all section cells (tissue context)", fontsize=13)
fig.tight_layout()
for ext in ("png", "pdf"):
    # PDF: points rasterized at 400 dpi (Illustrator-friendly), axes/text stay vector
    fig.savefig(outdir / f"spatial_clone_facets.{ext}", dpi=400)
print("saved", outdir / "spatial_clone_facets.png")
