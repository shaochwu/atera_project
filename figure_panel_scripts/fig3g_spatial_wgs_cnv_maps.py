#!/usr/bin/env python3
"""Fig 3g - laser-microdissected WGS regions coloured by copy number (Chr8 loss, Chr12&17 gain, Chr11 gain)."""
import re
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

WGS = "/path/to/atera_wgs/results"
CNS = f"{WGS}/cnvkit/cohort/all"
SVG = f"{WGS}/cnvkit_summary/wgs_region.svg"
BBOX = f"{WGS}/cnvkit_summary/wgs_region_bbox.html"
OUT = f"{WGS}/cnvkit_summary"
W, H = 960, 540           # SVG canvas
SCALE = 0.0019685         # SVG user-units -> canvas px
BIN = 500_000             # genomic bin width (bp)
CHR_SIZE = {"chr8": 145138636, "chr11": 135086622,
            "chr12": 133275309, "chr17": 83257441}

# panels: (title, chromosomes, direction, colormap, colorbar label)
PANELS = [
    ("Chr8 Loss", ["chr8"], "loss", plt.cm.Blues_r, "Mean Log2 ratio"),
    ("Average\nChr12 & 17 Gain", ["chr12", "chr17"], "gain", plt.cm.Reds, "Mean Log2 ratio"),
    ("Chr11 Gain", ["chr11"], "gain", plt.cm.Reds, "Mean Log2 ratio"),
]


# ----------------------------------------------------------- CNVkit -> scores
def key(f):
    m = re.search(r"cap_(\d+)_DL", f)
    return int(m.group(1)) if m else 999


files = [f for f in sorted(glob.glob(f"{CNS}/*.cns"), key=key)
         if not re.search(r"(call|bintest)", f)]
names = [os.path.basename(f).replace(".cns", "") for f in files]
dfs = {s: pd.read_csv(f, sep="\t") for s, f in zip(names, files)}


def binv(df, chrom):
    """Per-bin log2 for one chromosome (nan where no segment), clipped to [-2, 2]."""
    n = int(np.ceil(CHR_SIZE[chrom] / BIN))
    v = np.full(n, np.nan)
    s = df[df.chromosome == chrom]
    if not s.empty:
        mids = (np.arange(n) + 0.5) * BIN
        st, en, lg = s.start.values, s.end.values, s.log2.values
        idx = np.clip(np.searchsorted(st, mids, side="right") - 1, 0, len(st) - 1)
        v = np.where((mids >= st[idx]) & (mids < en[idx]), lg[idx], np.nan)
    return np.clip(v, -2, 2)


def chrom_score(chrom, direction):
    """Per-sample mean log2 over the recurrently-altered bins of `chrom`."""
    M = np.vstack([binv(dfs[s], chrom) for s in names])
    mean = np.nanmean(M, axis=0)
    alt = mean < -0.15 if direction == "loss" else mean > 0.15
    return {s: float(np.nanmean(M[i, alt])) for i, s in enumerate(names)}


def panel_score(chroms, direction):
    """Score per sample: mean log2 for one chr, or the equal-weight mean across chrs."""
    per = {c: chrom_score(c, direction) for c in chroms}
    return {s: float(np.mean([per[c][s] for c in chroms])) for s in names}


# ----------------------------------------------------- region polygons + labels
labels = {}
for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
                     r'yMax="([\d.]+)">(\d+)</word>', open(BBOX).read()):
    x0, y0, x1, y1, txt = m.groups()
    cy = (float(y0) + float(y1)) / 2
    if cy > 420:
        continue
    labels[int(txt)] = ((float(x0) + float(x1)) / 2, cy)

body = open(SVG).read().split("</defs>")[1]


def parse_d(d):
    toks = re.findall(r'[MLCZ]|-?\d*\.?\d+(?:e-?\d+)?', d)
    toks = [(str(float(t) * SCALE) if re.match(r'-?\d*\.?\d', t) else t) for t in toks]
    verts, codes, i, start, cur = [], [], 0, None, None
    while i < len(toks):
        t = toks[i]
        if t == "M":
            x, y = float(toks[i + 1]), float(toks[i + 2]); i += 3
            verts.append((x, y)); codes.append(MPath.MOVETO); start = cur = (x, y)
        elif t == "L":
            x, y = float(toks[i + 1]), float(toks[i + 2]); i += 3
            verts.append((x, y)); codes.append(MPath.LINETO); cur = (x, y)
        elif t == "C":
            p = [float(v) for v in toks[i + 1:i + 7]]; i += 7
            verts += [(p[0], p[1]), (p[2], p[3]), (p[4], p[5])]
            codes += [MPath.CURVE4] * 3; cur = (p[4], p[5])
        elif t == "Z":
            verts.append(start or cur); codes.append(MPath.CLOSEPOLY); i += 1
        else:
            i += 1
    return verts, codes


paths = []
for m in re.finditer(r'<path\b([^>]*?)\bd="([^"]*)"([^>]*?)/>', body):
    dash = "stroke-dasharray" in (m[1] + m[3])
    verts, codes = parse_d(m[2])
    if len(verts) < 3:
        continue
    v = np.array(verts)
    paths.append(dict(mp=MPath(verts, codes), pts=v, dash=dash,
                      bbox=(v[:, 0].min(), v[:, 1].min(), v[:, 0].max(), v[:, 1].max())))


def inpoly(pts, x, y):
    n = len(pts); ins = False; j = n - 1
    for i in range(n):
        xi, yi = pts[i]; xj, yj = pts[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi):
            ins = not ins
        j = i
    return ins


area = lambda b: (b[2] - b[0]) * (b[3] - b[1])

# match each numbered label to the smallest polygon containing it
match = {}
for num, (x, y) in labels.items():
    cand = [(i, area(p["bbox"])) for i, p in enumerate(paths) if inpoly(p["pts"], x, y)]
    if cand:
        match[min(cand, key=lambda t: t[1])[0]] = num
# regions 1 & 2 share one (lower-left) polygon with no separate number
tissue = max(range(len(paths)), key=lambda i: area(paths[i]["bbox"]))
R12_idx = None
for i, p in enumerate(paths):
    if i in match or i == tissue:
        continue
    cx, cy = p["pts"][:, 0].mean(), p["pts"][:, 1].mean()
    if cx < W / 2 and cy > H / 2:
        if R12_idx is None or area(paths[i]["bbox"]) > area(paths[R12_idx]["bbox"]):
            R12_idx = i
if R12_idx is not None:
    match[R12_idx] = "1,2"
    R12c = (paths[R12_idx]["pts"][:, 0].mean(), paths[R12_idx]["pts"][:, 1].mean())


def region_val(score, lab):
    if lab == "1,2":
        return float(np.mean([score["cap_1_DL"], score["cap_2_DL"]]))
    return score.get(f"cap_{lab}_DL")


# ------------------------------------------------------------------- figure
region_labs = list(labels) + (["1,2"] if R12_idx is not None else [])
fig, axes = plt.subplots(1, len(PANELS), figsize=(4.6 * len(PANELS), 5.2))
all_scores = {}
for ax, (title, chroms, direction, cmap, cbar_lab) in zip(np.ravel(axes), PANELS):
    score = panel_score(chroms, direction)
    all_scores[title.replace("\n", " ")] = score
    vals = [region_val(score, l) for l in region_labs]
    norm = Normalize(vmin=min(vals), vmax=max(vals))
    for i, p in enumerate(paths):
        if i in match:
            val = region_val(score, match[i])
            ax.add_patch(PathPatch(
                p["mp"], facecolor=cmap(norm(val)) if val is not None else "white",
                edgecolor="0.35", lw=0.7, ls=(0, (4, 2)) if p["dash"] else "-",
                rasterized=True))
        else:
            big = area(p["bbox"]) > 0.5 * W * H
            ax.add_patch(PathPatch(p["mp"], facecolor="none", edgecolor="black",
                                   lw=1.0 if big else 0.6,
                                   ls="-" if big else (0, (4, 2))))
    for n, (x, y) in labels.items():
        ax.text(x, y, str(n), ha="center", va="center", fontsize=8, zorder=6)
    if R12_idx is not None:
        ax.text(*R12c, "1,2", ha="center", va="center", fontsize=8,
                fontweight="bold", zorder=6)
    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=12)
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02).set_label(cbar_lab)

plt.tight_layout()
plt.savefig(f"{OUT}/fig3g_wgs_region_cnv_maps.png", dpi=200)
plt.savefig(f"{OUT}/fig3g_wgs_region_cnv_maps.pdf", dpi=300)
pd.DataFrame(all_scores).to_csv(f"{OUT}/fig3g_wgs_region_cnv_scores.tsv", sep="\t")
print("regions matched:", sorted(match.values(), key=str))
print("saved ->", OUT)
