#!/usr/bin/env python3
"""Spatial map of epithelial DCIS states (Luminal DCIS1-4 + Merged Normal)."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# keep text as editable TrueType (Type 42) in the PDF, not Type 3 outlines,
# so labels/legend/titles stay live-editable in Illustrator
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="results/infercnv_run")
ap.add_argument("--loc", default="cell_locations.csv")
args = ap.parse_args()

DIR = Path(args.dir); outdir = DIR / "biological_insight"

# state order (bottom->top of tree) + phylo-tree colours (fixed_pal, R names -> hex)
STATE_PAL = {
    "Luminal DCIS1": "#87CEFF",   # skyblue1
    "Luminal DCIS2": "#DAA520",   # (swapped with DCIS3)
    "Luminal DCIS3": "#CD2626",   # (swapped with DCIS2)
    "Luminal DCIS4": "#7D26CD",   # purple3
    "Merged Normal": "#9ACD32",
}
NORMAL_MAP = {"Ductal Luminal": "Merged Normal", "Normal Luminal": "Merged Normal"}
order = ["Luminal DCIS1", "Luminal DCIS2", "Luminal DCIS3", "Luminal DCIS4", "Merged Normal"]

loc = pd.read_csv(args.loc)
loc["cell_id"] = loc["cell_id"].astype(str)
loc = loc.dropna(subset=["x", "y"])
BG = loc[["x", "y"]].values                                  # full section = context

loc["state"] = loc["Cell_states"].astype(str).replace(NORMAL_MAP)
df = loc[loc["state"].isin(order)].copy()
n_by = df["state"].value_counts()
print(f"{len(df)} epithelial-state cells of {len(BG)} total section cells "
      f"({len(df)/len(BG)*100:.0f}%)")
for s in order:
    print(f"  {s:<15} n={int(n_by.get(s, 0)):,}")

XY = df[["x", "y"]].values
stv = df["state"].values

# ---- combined panel + one panel per state ---------------------------------
ncol = 3
nrow = int(np.ceil((len(order) + 1) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 3.3 * nrow))
axes = axes.ravel()
for ax in axes:
    ax.axis("off")

# panel 0: all states together
ax = axes[0]; ax.axis("on"); ax.set_xticks([]); ax.set_yticks([])
ax.scatter(BG[:, 0], BG[:, 1], s=0.4, c="#ECECEC", linewidths=0, rasterized=True)
for c in order:
    s = stv == c
    ax.scatter(XY[s, 0], XY[s, 1], s=1.4, c=[STATE_PAL[c]], linewidths=0,
               alpha=1.0, rasterized=True, label=f"{c}  (n={int(s.sum()):,})")
ax.set_title("all epithelial states", fontsize=10)
ax.invert_yaxis(); ax.set_aspect("equal")
ax.legend(markerscale=8, fontsize=6, loc="upper left",
          bbox_to_anchor=(1.0, 1.0), frameon=False)

# one highlighted state per panel
for j, c in enumerate(order, start=1):
    ax = axes[j]; ax.axis("on"); ax.set_xticks([]); ax.set_yticks([])
    s = stv == c
    ax.scatter(BG[:, 0], BG[:, 1], s=0.4, c="#ECECEC", linewidths=0, rasterized=True)
    ax.scatter(XY[s, 0], XY[s, 1], s=1.6, c=[STATE_PAL[c]], linewidths=0,
               alpha=1.0, rasterized=True)
    ax.set_title(f"{c}  (n={int(s.sum()):,})", fontsize=10)
    ax.invert_yaxis(); ax.set_aspect("equal")

fig.suptitle("Spatial distribution of Luminal DCIS states (≥ 3000 genes) — "
             "grey = all section cells (tissue context)", fontsize=13)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(outdir / f"spatial_lumstate_facets.{ext}", dpi=400)
print("saved", outdir / "spatial_lumstate_facets.png")

# ---- standalone single combined panel (larger) ----------------------------
fig2, ax = plt.subplots(figsize=(8, 8))
ax.set_xticks([]); ax.set_yticks([])
ax.scatter(BG[:, 0], BG[:, 1], s=0.5, c="#ECECEC", linewidths=0, rasterized=True)
for c in order:
    s = stv == c
    ax.scatter(XY[s, 0], XY[s, 1], s=2.0, c=[STATE_PAL[c]], linewidths=0,
               alpha=1.0, rasterized=True, label=f"{c}  (n={int(s.sum()):,})")
ax.set_title("Luminal DCIS states + Merged Normal (≥ 3000 genes)", fontsize=13)
ax.invert_yaxis(); ax.set_aspect("equal")
ax.legend(markerscale=6, fontsize=9, loc="upper left",
          bbox_to_anchor=(1.0, 1.0), frameon=False)
fig2.tight_layout()
for ext in ("png", "pdf"):
    fig2.savefig(outdir / f"spatial_lumstate.{ext}", dpi=400)
print("saved", outdir / "spatial_lumstate.png")
