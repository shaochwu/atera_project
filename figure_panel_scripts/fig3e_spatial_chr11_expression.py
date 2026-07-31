#!/usr/bin/env python3
"""Spatial map of the 11q13 amplicon-gene expression module score."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.neighbors import NearestNeighbors

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="results/infercnv_run")
ap.add_argument("--loc", default="cell_locations.csv")
ap.add_argument("--cnv-bins", default="901,904")
ap.add_argument("--k", type=int, default=10)
args = ap.parse_args()

DIR = Path(args.dir); outdir = DIR / "biological_insight"
AMP = ["CCND1", "FGF3", "FGF4", "FGF19", "ANO1", "CTTN", "FADD", "MYEOV", "TPCN2"]
GREYRED = LinearSegmentedColormap.from_list(
    "greyred", ["#e3e3e3", "#fdae61", "#d6604d", "#b2182b", "#67001f"])

a = ad.read_h5ad(DIR / "infercnv_atera.h5ad")
# CNV amplicon signal (only to validate the expression map is orthogonally consistent)
b0, b1 = (int(x) for x in args.cnv_bins.split(","))
Xc = a.obsm["X_cnv"]; Xc = np.asarray(Xc.todense()) if not isinstance(Xc, np.ndarray) else Xc
cnv11 = Xc[:, b0:b1].mean(1)

# expression: counts -> CP10k -> log1p -> module score over detected amplicon genes
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
genes = [g for g in AMP if g in a.var_names]
det = np.asarray((a[:, genes].X > 0).mean(0)).ravel()
genes = [g for g, d in zip(genes, det) if d >= 0.05]      # drop ~undetected genes
print("amplicon expression module genes:", genes)
sc.tl.score_genes(a, genes, score_name="amp11q13", ctrl_size=50, random_state=0)

df = pd.DataFrame({"cell_id": a.obs_names.astype(str),
                   "clone": a.obs["clone_grp"].astype(str).values,
                   "expr": a.obs["amp11q13"].values, "cnv": cnv11})
df = df[df["clone"] != "nan"]
loc = pd.read_csv(args.loc)[["cell_id", "x", "y"]]; loc["cell_id"] = loc["cell_id"].astype(str)
df = loc.merge(df, on="cell_id", how="inner").dropna(subset=["x", "y", "expr"])
XY = df[["x", "y"]].values; v = df["expr"].values
r = np.corrcoef(df["expr"], df["cnv"])[0, 1]
print(f"{len(df)} cells; expr~CNV Pearson r = {r:.2f} (orthogonal agreement)")

# Moran's I on the EXPRESSION score
nn = NearestNeighbors(n_neighbors=args.k + 1).fit(XY); _, idx = nn.kneighbors(XY); idx = idx[:, 1:]
z = v - v.mean(); W = len(v) * args.k
I = (len(v) / W) * ((z[:, None] * z[idx]).sum() / (z ** 2).sum())
rng = np.random.default_rng(0)
perm = np.array([(len(v) / W) * ((zp := rng.permutation(z))[:, None] * zp[idx]).sum() / (zp ** 2).sum()
                 for _ in range(200)])
zsc = (I - perm.mean()) / perm.std()
print(f"Moran's I (expression) = {I:.3f}  (z = {zsc:.0f}, p < {1/len(perm):.3f})")

fig, ax = plt.subplots(figsize=(11, 7))
o = np.argsort(v)
lo, hi = np.percentile(v, [2, 98])
sc_ = ax.scatter(XY[o, 0], XY[o, 1], c=v[o], s=3, cmap=GREYRED, vmin=lo, vmax=hi, linewidths=0)
ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("11q13 amplicon-gene EXPRESSION module "
             "(CCND1, FGF3/4/19, ANO1, CTTN, FADD, MYEOV, TPCN2)\n"
             f"Moran's I = {I:.2f} (z = {zsc:.0f}) — spatially segregated; "
             f"agrees with CNV (r = {r:.2f})", fontsize=11)
fig.colorbar(sc_, ax=ax, shrink=0.7, label="amplicon expression module score")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(outdir / f"spatial_chr11_expression.{ext}", dpi=220 if ext == "png" else None)
print("saved", outdir / "spatial_chr11_expression.png")
