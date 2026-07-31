#!/usr/bin/env python3
"""Build a CNV phylogenetic tree from the 16 WGS samples (NJ on CNVkit log2; arm-level Sankoff CN events)."""

import argparse
import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# hg38 chromosome sizes (autosomes) and centromere midpoints (bp) for arm splitting
CHR_SIZE = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}
CEN = {
    "chr1": 123400000, "chr2": 93900000, "chr3": 90900000, "chr4": 50000000,
    "chr5": 48800000, "chr6": 59800000, "chr7": 60100000, "chr8": 45200000,
    "chr9": 43000000, "chr10": 39800000, "chr11": 53400000, "chr12": 35500000,
    "chr13": 17700000, "chr14": 17200000, "chr15": 19000000, "chr16": 36800000,
    "chr17": 25100000, "chr18": 18500000, "chr19": 26200000, "chr20": 28100000,
    "chr21": 12000000, "chr22": 15000000,
}
CHR_ORDER = [f"chr{i}" for i in range(1, 23)]

# Copy-number regions to consider for branch labels: EXACTLY the set that appears
# in the inferCNV example tree cnv_phylo_tree_ggtree_mp_clones.pdf. Bare numbers =
# whole chromosome; NqNp = arm. (matches clone_leaf_states_50_merge34.tsv regions)
EXAMPLE_REGIONS = ["17", "12", "1q", "13", "11", "8q", "14", "8p",
                   "4", "18", "5", "15", "20"]
# Focal amplicons handled exactly like the example script's FOCAL dict: a focal
# region that washes out in a whole-chromosome average, placed ONCE at the
# carriers' MRCA and labelled on the branch by CHROMOSOME NUMBER (disp). The
# whole chromosome is represented by this focal amplicon instead of its diluted
# genome-wide average (chr11 is a net loss but its real event is the 11q13/CCND1
# gain; chr12's is the 12q14-15 HMGA2/MDM2 gain). name -> (chrom, start, end, disp)
FOCAL = {
    "11q13(CCND1)": ("chr11", 69_000_000, 70_200_000, "11"),
    "12q14-15(HMGA2/MDM2)": ("chr12", 65_000_000, 69_500_000, "12"),
}
# WGS samples used as the normal/root outgroup (near-diploid)
OUTGROUP = ["cap_1_DL", "cap_2_DL"]
# Known DCIS truncal driver gains to FORCE onto the tumour trunk even when the WGS
# amplitude is at the noise floor (per user; matches the inferCNV clone tree). The
# gain is placed at the MRCA of all non-normal samples; any non-normal sample that
# genuinely lacks it (centred log2 <= 0) gets an explicit reversion on its branch.
FORCE_TRUNCAL_GAINS = ["1q"]


ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}  # acrocentric: q arm only


def region_bounds(region):
    """Map a region label to (chrom, start, end).

    Bare number -> whole chromosome; trailing p/q -> that arm (split at centromere).
    """
    if region[-1] in "pq":
        chrom = "chr" + region[:-1]
        cen = CEN[chrom]
        return (chrom, 0, cen) if region[-1] == "p" else (chrom, cen, CHR_SIZE[chrom])
    chrom = "chr" + region
    return chrom, 0, CHR_SIZE[chrom]


def expand_to_arms(region):
    """Expand a bare-chromosome region to its arm(s); pass arm tokens through.

    Whole chromosomes are scored per ARM so arm-divergent events are not washed
    out by averaging (e.g. chr17 = 17q gain + 17p loss cancels at whole-chr level).
    """
    if region[-1] in "pq":
        return [region]
    ch = "chr" + region
    return [region + "q"] if ch in ACRO else [region + "p", region + "q"]


def weighted_mean_log2(seg, chrom, start, end):
    """Length-weighted mean log2 of segments overlapping [start,end) on chrom."""
    sub = seg[(seg.chromosome == chrom) & (seg.end > start) & (seg.start < end)]
    if sub.empty:
        return np.nan
    ov = (np.minimum(sub.end.values, end) - np.maximum(sub.start.values, start))
    ov = ov.clip(min=0).astype(float)
    if ov.sum() == 0:
        return np.nan
    return float(np.average(sub.log2.values, weights=ov))


def load_segs(cns_dir):
    """Return {sample: CNVkit .cns segment DataFrame} for arm/focal event calls."""
    files = sorted(glob.glob(str(Path(cns_dir) / "*.cns")))
    files = [f for f in files if not f.endswith((".call.cns", ".bintest.cns"))]
    return {Path(f).name.replace(".cns", ""): pd.read_csv(f, sep="\t") for f in files}


# ---------------------------------------------------------- neighbor joining
def neighbor_joining(dist, names):
    n = len(names)
    nodes = list(names)
    D = {a: {b: float(dist[i, j]) for j, b in enumerate(names)}
         for i, a in enumerate(names)}
    adj = {a: {} for a in names}
    inext = 0
    while len(nodes) > 2:
        m = len(nodes)
        r = {a: sum(D[a][b] for b in nodes if b != a) for a in nodes}
        best, bi, bj = None, None, None
        for i in range(m):
            for j in range(i + 1, m):
                a, b = nodes[i], nodes[j]
                q = (m - 2) * D[a][b] - r[a] - r[b]
                if best is None or q < best:
                    best, bi, bj = q, a, b
        a, b = bi, bj
        u = f"I{inext}"
        inext += 1
        dau = 0.5 * D[a][b] + (r[a] - r[b]) / (2 * (m - 2))
        dbu = D[a][b] - dau
        adj[u] = {}
        adj[a][u] = max(dau, 0.0)
        adj[u][a] = max(dau, 0.0)
        adj[b][u] = max(dbu, 0.0)
        adj[u][b] = max(dbu, 0.0)
        D[u] = {}
        for k in nodes:
            if k in (a, b):
                continue
            d = 0.5 * (D[a][k] + D[b][k] - D[a][b])
            D[u][k] = d
            D[k][u] = d
        D[u][u] = 0.0
        nodes = [x for x in nodes if x not in (a, b)] + [u]
    a, b = nodes
    L = max(D[a][b], 0.0)
    adj[a][b] = L
    adj[b][a] = L
    return adj


def _dist_matrix(rows):
    M = np.asarray(rows, float)
    return np.sqrt(((M[:, None, :] - M[None, :, :]) ** 2).sum(2))


def build_adjacency(prof_df, names, group):
    """NJ adjacency. If `group` (>=2 of `names`) is given, FORCE it to be a
    monophyletic clade: collapse the group to its centroid for the main NJ, then
    graft the group's own little subtree back in place of the centroid leaf.

    This is a PRESENTATION constraint, used only because the 12 tumour samples are
    effectively one clone (pairwise r 0.96-0.996) so their sub-structure is a weak
    purity gradient, not robust evolutionary signal. Everything downstream (Sankoff
    events, rooting) recomputes on the resulting topology.
    """
    group = [g for g in group if g in names]
    if len(group) < 2:
        return neighbor_joining(_dist_matrix([prof_df.loc[n].values for n in names]),
                                names)
    rep = "__GRP__"
    others = [n for n in names if n not in group]
    rows = [prof_df.loc[n].values for n in others] + [prof_df.loc[group].values.mean(0)]
    adj = neighbor_joining(_dist_matrix(rows), others + [rep])
    nb = next(iter(adj[rep]))
    L = adj[rep][nb]
    del adj[rep]
    del adj[nb][rep]

    gd = pd.DataFrame(_dist_matrix([prof_df.loc[g].values for g in group]),
                      index=group, columns=group)

    def link(u, v, w):
        adj.setdefault(u, {})[v] = max(w, 0.0)
        adj.setdefault(v, {})[u] = max(w, 0.0)

    if len(group) == 2:
        a, b = group
        link("GRP_ROOT", a, gd.loc[a, b] / 2)
        link("GRP_ROOT", b, gd.loc[a, b] / 2)
    else:
        pairs = [(group[i], group[j]) for i in range(len(group))
                 for j in range(i + 1, len(group))]
        x, y = min(pairs, key=lambda p: gd.loc[p[0], p[1]])
        rest = [g for g in group if g not in (x, y)]
        link("GRP_INNER", x, gd.loc[x, y] / 2)
        link("GRP_INNER", y, gd.loc[x, y] / 2)
        link("GRP_ROOT", "GRP_INNER",
             max(gd.loc[x, rest[0]] - gd.loc[x, y] / 2, 0.05))
        for z in rest:
            link("GRP_ROOT", z,
                 max((gd.loc[x, z] + gd.loc[y, z]) / 2 - gd.loc[x, y] / 2, 0.05))
    link("GRP_ROOT", nb, L)
    return adj


def root_on_outgroup(adj, outgroup):
    """Root on the edge whose removal splits `outgroup` (monophyletic) from the rest.

    Places the root at the midpoint of that edge, so the outgroup clade and the
    ingroup become the two children of the root.
    """
    S = set(outgroup)

    def side_tips(a, b):
        """Tips reachable from b without crossing back through a."""
        seen, stack, res = {a}, [b], set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            if len(adj[x]) == 1:
                res.add(x)
            for k in adj[x]:
                if k not in seen:
                    stack.append(k)
        return res

    edge = None
    for a in list(adj):
        for b in list(adj[a]):
            if side_tips(a, b) == S:   # b-side == outgroup
                edge = (a, b)
                break
        if edge:
            break
    if edge is None:
        raise ValueError(f"outgroup {outgroup} is not monophyletic in the NJ tree")

    a, b = edge   # root goes between a (ingroup side) and b (outgroup side)
    L = adj[a][b]
    del adj[a][b]
    del adj[b][a]
    root = "ROOT"
    adj[root] = {a: L / 2.0, b: L / 2.0}
    adj[a][root] = L / 2.0
    adj[b][root] = L / 2.0

    children, blen, parent = {}, {}, {}

    def dfs(node, par):
        parent[node] = par
        kids = [k for k in adj[node] if k != par]
        children[node] = kids
        for k in kids:
            blen[k] = adj[node][k]
            dfs(k, node)

    blen[root] = 0.0
    dfs(root, None)
    return root, children, parent, blen


# --------------------------------------------------------- Sankoff parsimony
STATES = [-1, 0, 1]


def sankoff(root, children, leaf_state):
    cost = {}

    def up(node):
        if not children[node]:
            cost[node] = {s: (0.0 if s == leaf_state[node] else float("inf"))
                          for s in STATES}
            return
        for c in children[node]:
            up(c)
        g = {}
        for s in STATES:
            tot = 0.0
            for c in children[node]:
                tot += min(abs(s - t) + cost[c][t] for t in STATES)
            g[s] = tot
        cost[node] = g

    up(root)
    assign = {}

    def down(node, s):
        assign[node] = s
        for c in children[node]:
            best_t, best_v = None, None
            for t in STATES:
                v = abs(s - t) + cost[c][t]
                if best_v is None or v < best_v:
                    best_v, best_t = v, t
            down(c, best_t)

    down(root, 0)  # diploid ancestor
    return assign


# ------------------------------------------------------------------- layout
def layout(root, children, blen):
    x = {}

    def setx(node, acc):
        acc = acc + blen.get(node, 0.0)
        x[node] = acc
        for c in children[node]:
            setx(c, acc)

    setx(root, 0.0)

    leaves = []

    def collect(node):
        if not children[node]:
            leaves.append(node)
        for c in children[node]:
            collect(c)

    collect(root)
    # place leaves top->bottom with the outgroup (root clade) at the BOTTOM,
    # matching the example tree (Merged_Normal drawn at the bottom)
    nlv = len(leaves)
    y = {lf: float(nlv - 1 - i) for i, lf in enumerate(leaves)}

    def sety(node):
        if children[node]:
            for c in children[node]:
                sety(c)
            y[node] = float(np.mean([y[c] for c in children[node]]))

    sety(root)
    return x, y, leaves


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cns-dir",
                    default="/path/to/atera_wgs/results/cnvkit/cohort/all")
    ap.add_argument("--window-profile",
                    default="/path/to/atera_wgs/results/cnvkit_summary/"
                            "cnvkit_on_infercnv_windows_log2.tsv",
                    help="WGS log2 on the 1516 inferCNV gene-windows (the heatmap "
                         "grid) — used as the tree distance metric")
    ap.add_argument("--outdir",
                    default="/path/to/atera_wgs/results/biological_insight")
    ap.add_argument("--gain-thr", type=float, default=0.13,
                    help="normal-centred log2 thr for arm gains (compressed WGS)")
    ap.add_argument("--loss-thr", type=float, default=-0.13,
                    help="normal-centred log2 thr for arm losses")
    ap.add_argument("--focal-gain-thr", type=float, default=0.5,
                    help="normal-centred log2 thr for focal-amplicon gains")
    ap.add_argument("--force-group", default="cap_11_DL,cap_12_DL,cap_14_DL",
                    help="comma-separated samples forced into one monophyletic "
                         "clade (presentation constraint; '' to disable)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- profile on the 1516 inferCNV gene-windows (distance metric) ----
    # "follow the heatmap": use the same autosomal inferCNV window grid as the
    # clone tree / heatmap, so the WGS tree is directly comparable to inferCNV.
    prof = pd.read_csv(args.window_profile, sep="\t", index_col=0).fillna(0.0)
    prof.columns = [f"w{c}" for c in prof.columns]
    segs = load_segs(args.cns_dir)          # segments -> arm/focal event calls
    samples = [s for s in prof.index if s in segs]
    prof = prof.loc[samples]
    print(f"Loaded {len(samples)} WGS samples, {prof.shape[1]} inferCNV windows")

    # tree over the real samples only; rooted on the near-diploid cap_1/cap_2
    # outgroup (the WGS "normal"), so the tree runs normal -> tumour.
    names = list(prof.index)
    missing = [s for s in OUTGROUP if s not in names]
    if missing:
        raise ValueError(f"outgroup samples not found: {missing}")

    group = [s.strip() for s in args.force_group.split(",") if s.strip()]
    adj = build_adjacency(prof, names, group)
    if group:
        print(f"forced monophyletic group (presentation): {group}")
    root, children, parent, blen = root_on_outgroup(adj, OUTGROUP)

    # ---- per-arm length-weighted log2 -> discrete calls per sample ----
    # regions restricted to the inferCNV example-tree chromosomes, but scored per
    # ARM (whole-chromosome averaging cancels arm-divergent events like chr17).
    # chromosomes with a focal amplicon (11, 12) are dropped here and represented
    # by their focal amplicon below.
    focal_chroms = {disp for _, _, _, disp in FOCAL.values()}
    arm_order = [a for r in EXAMPLE_REGIONS if r not in focal_chroms
                 for a in expand_to_arms(r)]
    arm_log2 = {}
    for s in samples:
        seg = segs[s]
        arm_log2[s] = {}
        for lab in arm_order:
            ch, st, en = region_bounds(lab)
            arm_log2[s][lab] = weighted_mean_log2(seg, ch, st, en)
    arm_log2_df = pd.DataFrame(arm_log2).T[arm_order]

    # centre on the cohort normal (cap_1/cap_2) baseline: the tumour-only
    # flat-reference CNVkit log2 is compressed and offset per arm (e.g. 8p baseline
    # -0.18, 20q +0.12), so events are called RELATIVE to the normals. This
    # surfaces consistent low-amplitude arm changes (e.g. +8q, ~+0.1 raw) that a
    # fixed absolute threshold misses, while normals stay neutral by construction.
    baseline = arm_log2_df.loc[OUTGROUP].mean()
    arm_log2_c = arm_log2_df.sub(baseline, axis=1)

    def classify(x):
        if np.isnan(x):
            return 0
        return 1 if x > args.gain_thr else (-1 if x < args.loss_thr else 0)

    leaf_disc = {n: {} for n in names}
    for n in names:
        for lab in arm_order:
            leaf_disc[n][lab] = classify(arm_log2_c.loc[n, lab])

    # ---- Sankoff per region -> ancestral discrete states for every node ----
    all_nodes = list(blen.keys())
    node_state = {n: {} for n in all_nodes}
    for lab in arm_order:
        ls = {n: leaf_disc[n][lab] for n in names}
        assign = sankoff(root, children, ls)
        for n, s in assign.items():
            node_state[n][lab] = s

    # branch events (parent -> child)
    events = {n: ([], []) for n in all_nodes}
    for child in all_nodes:
        par = parent[child]
        if par is None:
            continue
        for lab in arm_order:
            d = node_state[child][lab] - node_state[par][lab]
            if d > 0:
                events[child][0].append(lab)
            elif d < 0:
                events[child][1].append(lab)

    def descendant_tips(node):
        if not children[node]:
            return [node]
        out = []
        for c in children[node]:
            out.extend(descendant_tips(c))
        return out

    # ---- focal amplicons: placed ONCE at the MRCA of carrier samples ----
    # (surfaces gains that a whole-chromosome average hides, e.g. 11q13/CCND1
    # sitting inside the broad chr11 loss). Mirrors the example tree's FOCAL dict.
    def ancestors(node):
        path, cur = [], node
        while cur is not None:
            path.append(cur)
            cur = parent[cur]
        return path

    depth = {n: len(ancestors(n)) - 1 for n in all_nodes}

    def mrca(nodes):
        common = set(ancestors(nodes[0]))
        for n in nodes[1:]:
            common &= set(ancestors(n))
        return max(common, key=lambda c: depth[c])

    # ---- force known DCIS truncal driver gains onto the tumour trunk ----
    nonnormal = [s for s in samples if s not in OUTGROUP]
    trunk = mrca(nonnormal)
    for tok in FORCE_TRUNCAL_GAINS:
        if tok not in arm_order:
            continue
        # drop any Sankoff-placed scattered events for this token
        for nd in all_nodes:
            events[nd][0][:] = [t for t in events[nd][0] if t != tok]
            events[nd][1][:] = [t for t in events[nd][1] if t != tok]
        events[trunk][0].append(tok)                    # +tok on the trunk
        absent = [s for s in nonnormal if arm_log2_c.loc[s, tok] <= 0]
        for s in absent:                                # reversion where truly absent
            events[s][1].append(tok)
        print(f"forced truncal +{tok} at {trunk}"
              + (f"; reversion on {absent}" if absent else ""))

    focal_log2 = {}
    for label, (fch, fs, fe, disp) in FOCAL.items():
        vals = {s: weighted_mean_log2(segs[s], fch, fs, fe) for s in samples}
        focal_log2[label] = vals
        base = np.nanmean([vals[s] for s in OUTGROUP])   # centre on normals
        carriers = [s for s in samples
                    if not np.isnan(vals[s]) and (vals[s] - base) > args.focal_gain_thr]
        if carriers:
            events[mrca(carriers)][0].append(disp)   # branch label = chr number
            print(f"focal gain {label} (+{disp}): {len(carriers)} carriers -> "
                  f"MRCA {mrca(carriers)}")

    # ---- tip CNV-burden class (honest, self-contained sample annotation) ----
    # burden is counted genome-wide over ALL autosome arms (not just the 13
    # branch-label regions) so it still separates diploid-like (0 altered arms),
    # intermediate (1-3) and tumour (>=4); analogous to the cell-state colours on
    # the inferCNV tree.
    all_arms = []
    for ch in CHR_ORDER:
        num, cen = ch[3:], CEN[ch]
        if ch not in {"chr13", "chr14", "chr15", "chr21", "chr22"}:  # skip acro p
            all_arms.append((f"{num}p", ch, 0, cen))
        all_arms.append((f"{num}q", ch, cen, CHR_SIZE[ch]))
    # genome-wide arm log2; the burden count uses a robust ABSOLUTE 0.25 threshold
    # (not the sensitive centred one used for branch labels) so the tip counts stay
    # a clean summary: normals 0, cap_4 ~1, tumours ~6-11.
    aa_log2 = pd.DataFrame(
        {s: {lab: weighted_mean_log2(segs[s], c, st, en)
             for lab, c, st, en in all_arms} for s in samples}).T
    tip_class = {}   # sample -> (label, n_altered_arms genome-wide, absolute 0.25)
    for s in samples:
        n_alt = int((aa_log2.loc[s].abs() > 0.25).sum())
        if n_alt == 0:
            cls = "Diploid-like"
        elif n_alt <= 3:
            cls = "Intermediate"
        else:
            cls = "Tumour (high CNV)"
        tip_class[s] = (cls, n_alt)

    # ------------------------------------------------------------- outputs
    rows = []
    for child, (g, l) in events.items():
        if not g and not l:
            continue
        rows.append({
            "branch_to": child,
            "from": parent[child],
            "gains": ",".join(g),
            "losses": ",".join(l),
            "descendant_tips": ",".join(descendant_tips(child)),
        })
    ev_df = pd.DataFrame(rows)
    ev_df.to_csv(outdir / "wgs_cnv_phylo_branch_events.tsv", sep="\t", index=False)

    st_df = pd.DataFrame(node_state).T[arm_order]
    st_df.index.name = "node"
    st_df.to_csv(outdir / "wgs_cnv_phylo_arm_states.tsv", sep="\t")
    arm_log2_df.round(3).to_csv(outdir / "wgs_cnv_phylo_arm_log2.tsv", sep="\t")
    arm_log2_c.round(3).to_csv(outdir / "wgs_cnv_phylo_arm_log2_normcentred.tsv",
                               sep="\t")
    if focal_log2:
        pd.DataFrame(focal_log2).round(3).to_csv(
            outdir / "wgs_cnv_phylo_focal_log2.tsv", sep="\t")

    ts_df = pd.DataFrame(
        [(s, cls, n) for s, (cls, n) in tip_class.items()],
        columns=["wgs_sample", "cnv_class", "n_altered_arms"],
    ).set_index("wgs_sample")
    ts_df.to_csv(outdir / "wgs_cnv_phylo_tip_states.tsv", sep="\t")

    def newick(node):
        if not children[node]:
            nm = (node.replace(",", "_").replace(":", "_")
                  .replace("(", "").replace(")", "").replace(" ", "_"))
            return f"{nm}:{blen.get(node, 0):.5f}"
        inner = ",".join(newick(c) for c in children[node])
        return f"({inner}):{blen.get(node, 0):.5f}"

    (outdir / "wgs_cnv_phylo_tree.newick").write_text(newick(root) + ";")

    # ------------------------------------------------------------- figure
    x, y, leaves = layout(root, children, blen)
    # colour tips by CNV-burden class
    state_pal = {
        "Diploid-like": "#9ACD32", "Intermediate": "#DAA520",
        "Tumour (high CNV)": "#CD2626",
    }

    fig, ax = plt.subplots(figsize=(13, 8.5))
    for child in all_nodes:
        par = parent[child]
        if par is None:
            continue
        xp, xc, yc = x[par], x[child], y[child]
        ax.plot([xp, xp], [y[par], yc], color="0.3", lw=1.3, zorder=1)
        ax.plot([xp, xc], [yc, yc], color="0.3", lw=1.3, zorder=1)
        g, l = events.get(child, ([], []))
        parts = []
        if g:
            parts.append("+" + " +".join(g))
        if l:
            parts.append("-" + " -".join(l))
        if parts:
            gd, ld = bool(g), bool(l)
            color = "#b2182b" if gd and not ld else ("#2166ac" if ld and not gd else "#6a3d9a")
            ax.text((xp + xc) / 2.0, yc + 0.14, "  ".join(parts), fontsize=6.8,
                    color=color, ha="center", va="bottom", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=color,
                              lw=0.5, alpha=0.85))

    ax.set_xlim(-max(x.values()) * 0.02, max(x.values()) * 1.35)
    ax.set_ylim(-0.8, len(leaves) + 1.0)   # headroom so the top-left legend clears the top tip

    def tip_text(n):
        cls, na = tip_class.get(n, ("", 0))
        tag = "  [root/normal]" if n in OUTGROUP else ""
        return f"{n}  ({na} alt){tag}"

    for n in all_nodes:
        if not children[n]:
            is_out = n in OUTGROUP
            cls = tip_class.get(n, ("", 0))[0]
            col = state_pal.get(cls, "#333333")
            ax.scatter([x[n]], [y[n]], s=45 if is_out else 26, color=col,
                       edgecolors="#1b7837" if is_out else "none",
                       linewidths=1.4, zorder=4)
            ax.text(x[n] + max(x.values()) * 0.012, y[n], tip_text(n),
                    fontsize=9, va="center", ha="left",
                    fontweight="bold" if is_out else "normal")
        else:
            ax.scatter([x[n]], [y[n]], s=12, color="0.5", zorder=2)

    ax.set_title("WGS CNV phylogenetic tree of 16 samples\n"
                 "(rooted on cap_1/cap_2 normal outgroup; branch labels = "
                 "copy-number changes on inferCNV-tree chromosomes)", fontsize=12.5)
    ax.set_xlabel("CNV distance from cap_1/cap_2 (normal) root "
                  "(Euclidean on 1516 inferCNV gene-windows)")
    ax.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)

    handles = [
        mpatches.Patch(color="#b2182b", label="gain (+chr)"),
        mpatches.Patch(color="#2166ac", label="loss (-chr)"),
        mpatches.Patch(color="#6a3d9a", label="gain & loss"),
    ]
    leg1 = ax.legend(handles=handles, loc="lower right", fontsize=8,
                     frameon=False, title="branch CN change")
    ax.add_artist(leg1)
    used = [c for c in ["Diploid-like", "Intermediate", "Tumour (high CNV)"]
            if any(v[0] == c for v in tip_class.values())]
    shandles = [mpatches.Patch(color=state_pal[c], label=c) for c in used]
    if shandles:
        ax.legend(handles=shandles, loc="upper left", fontsize=8, frameon=False,
                  title="WGS CNV burden (tips)")

    fig.savefig(outdir / "wgs_cnv_phylo_tree.png", dpi=240, bbox_inches="tight")
    fig.savefig(outdir / "wgs_cnv_phylo_tree.pdf", bbox_inches="tight")
    plt.close(fig)

    print("Wrote tree to", outdir / "wgs_cnv_phylo_tree.png")
    print("Leaves:", len(leaves), "| branches with events:", len(ev_df))
    if not ev_df.empty:
        print(ev_df[["branch_to", "gains", "losses"]].to_string(index=False))


if __name__ == "__main__":
    main()
