#!/usr/bin/env python3
"""Build a CNV phylogenetic tree (NJ, diploid-rooted; Sankoff CN events) from inferCNV cluster centroids."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CHR_POS = {
    "chr1": 0, "chr2": 177, "chr3": 282, "chr4": 372, "chr5": 432, "chr6": 504,
    "chr7": 588, "chr8": 660, "chr9": 712, "chr10": 772, "chr11": 830,
    "chr12": 941, "chr13": 1026, "chr14": 1046, "chr15": 1094, "chr16": 1137,
    "chr17": 1202, "chr18": 1298, "chr19": 1314, "chr20": 1435, "chr21": 1476,
    "chr22": 1486,
}
ROOT_NAME = "Normal (diploid)"

# Focal regions (global X_cnv bin ranges, inclusive) labelled in addition to whole
# chromosomes, because focal amplicons wash out in a whole-chromosome average.
# inferCNV window_size=100, step=10. Bin ranges located from driver-gene positions.
#   11q13/CCND1        = bins 901-903
#   12q14-15/CDK4-MDM2 = bins 985-998  (strong amplicon diluted in whole-chr12 mean)
# name -> (start_bin, end_bin, chromosome-number label shown on the tree)
FOCAL = {
    "11q13(CCND1)": (901, 903, "11"),
    "12q14-15(CDK4/MDM2)": (985, 998, "12"),
}

# Optional custom regions (e.g. chromosome arms) used for branch CN labelling
# instead of whole chromosomes. Set from --regions; list of (label, start, end).
REGIONS = None
# Whole chromosomes to skip (handled instead by --extra-focal arm regions).
EXCLUDE_CHROMS = set()
# Extra point regions (e.g. specific arms) added on top of whole chromosomes,
# detected with their own threshold. name -> (start, end, display_label).
EXTRA = {}


# ---------------------------------------------------------------- data loading
def chrom_bounds(nbins):
    if REGIONS is not None:
        return REGIONS
    order = sorted(CHR_POS, key=lambda c: int(c[3:]))
    out = []
    for i, ch in enumerate(order):
        s = CHR_POS[ch]
        e = CHR_POS[order[i + 1]] if i + 1 < len(order) else nbins
        if ch not in EXCLUDE_CHROMS:          # drop after computing full bounds
            out.append((ch, s, e))
    return out


def chrom_matrix(profiles):
    bounds = chrom_bounds(profiles.shape[1])
    return pd.DataFrame(
        {ch: profiles.iloc[:, s:e].mean(axis=1) for ch, s, e in bounds},
        index=profiles.index,
    )


# ---------------------------------------------------------- neighbor joining
def neighbor_joining(dist, names):
    """Return adjacency dict {node: {neighbor: length}} for an unrooted NJ tree."""
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


def root_on_leaf(adj, outgroup):
    """Insert a root on the pendant edge of `outgroup`; return rooted tree dicts."""
    nb = next(iter(adj[outgroup]))
    L = adj[outgroup][nb]
    # detach
    del adj[outgroup][nb]
    del adj[nb][outgroup]
    root = "ROOT"
    adj[root] = {outgroup: L / 2.0, nb: L / 2.0}
    adj[outgroup][root] = L / 2.0
    adj[nb][root] = L / 2.0

    children, blen, parent = {}, {}, {}
    seen = set()

    def dfs(node, par):
        seen.add(node)
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
    """Ordered (|a-b|) parsimony; root forced to neutral. Returns {node: state}."""
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
    assign = {root: 0}  # diploid ancestor

    def down(node, s):
        assign[node] = s
        for c in children[node]:
            best_t, best_v = None, None
            for t in STATES:
                v = abs(s - t) + cost[c][t]
                if best_v is None or v < best_v:
                    best_v, best_t = v, t
            down(c, best_t)

    down(root, 0)
    return assign


# ------------------------------------------------------------------- layout
def layout(root, children, blen):
    """x = distance from root; y = leaf order. Returns x, y dicts and leaf order."""
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
    y = {}
    for i, lf in enumerate(leaves):
        y[lf] = float(i)

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
    ap.add_argument(
        "--profiles",
        default="results/infercnv_run/biological_insight/02_clone_mean_cnv_profile.tsv",
    )
    ap.add_argument(
        "--labels",
        default="results/infercnv_run/biological_insight/cluster_leaf_labels.tsv",
    )
    ap.add_argument(
        "--outdir",
        default="results/infercnv_run/biological_insight",
    )
    ap.add_argument("--gain-thr", type=float, default=0.03)
    ap.add_argument("--loss-thr", type=float, default=-0.02)
    ap.add_argument("--focal-gain-thr", type=float, default=0.05)
    ap.add_argument("--focal-loss-thr", type=float, default=-0.05)
    ap.add_argument("--merge", default="",
                    help="comma-separated clusters to merge into one leaf")
    ap.add_argument("--merge-name", default="Merged_Normal")
    ap.add_argument("--root-on", default="",
                    help="root on this real leaf instead of a synthetic diploid")
    ap.add_argument("--suffix", default="",
                    help="suffix appended to output filenames")
    ap.add_argument("--composition", default="",
                    help="clone x cell-state count table -> draw pies at the tips")
    ap.add_argument("--regions", default="",
                    help="TSV (label,start,end) of custom regions (e.g. chromosome "
                         "arms) used for branch labels instead of whole chromosomes; "
                         "disables focal regions")
    ap.add_argument("--exclude-chroms", default="",
                    help="comma-separated whole chromosomes to skip (e.g. chr1,chr8)")
    ap.add_argument("--extra-focal", default="",
                    help="TSV (name,start,end,label) of extra arm/point regions added "
                         "on top of whole chromosomes, with their own threshold")
    ap.add_argument("--extra-gain-thr", type=float, default=0.013)
    ap.add_argument("--extra-loss-thr", type=float, default=-0.013)
    ap.add_argument("--leaf-states", default="",
                    help="TSV (region x leaf, values -1/0/1) of explicit per-leaf CN "
                         "states; events placed by Sankoff parsimony")
    args = ap.parse_args()

    global REGIONS, EXCLUDE_CHROMS, EXTRA
    if args.regions:
        rdf = pd.read_csv(args.regions, sep="\t")
        REGIONS = [(str(r.label), int(r.start), int(r.end)) for r in rdf.itertuples()]
    if args.exclude_chroms:
        EXCLUDE_CHROMS = {c.strip() for c in args.exclude_chroms.split(",") if c.strip()}
    if args.extra_focal:
        edf = pd.read_csv(args.extra_focal, sep="\t")
        EXTRA = {}
        for r in edf.to_dict("records"):
            gthr = float(r["gain_thr"]) if r.get("gain_thr") == r.get("gain_thr") \
                and "gain_thr" in r else None
            lthr = float(r["loss_thr"]) if r.get("loss_thr") == r.get("loss_thr") \
                and "loss_thr" in r else None
            EXTRA[str(r["name"])] = (int(r["start"]), int(r["end"]), str(r["label"]),
                                     gthr, lthr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prof = pd.read_csv(args.profiles, sep="\t", index_col=0)
    prof.index = prof.index.astype(str)
    lab = pd.read_csv(args.labels, sep="\t", index_col=0)
    lab.index = lab.index.astype(str)

    # optional per-clone cell-state composition (counts) for tip pie charts
    comp = None
    comp_path = args.composition or str(Path(args.labels).with_name(
        "clone_state_composition.tsv"))
    if Path(comp_path).exists():
        comp = pd.read_csv(comp_path, sep="\t", index_col=0)
        comp.index = comp.index.astype(str)

    # optionally merge a set of clusters into a single (cell-weighted) leaf
    merge_set = [c.strip() for c in args.merge.split(",") if c.strip()]
    if merge_set:
        missing = [c for c in merge_set if c not in prof.index]
        if missing:
            raise ValueError(f"merge clusters not found: {missing}")
        w = lab.loc[merge_set, "n_cells"].astype(float)
        merged_profile = (prof.loc[merge_set].mul(w.values, axis=0).sum(axis=0)
                          / w.sum())
        prof = prof.drop(index=merge_set)
        prof.loc[args.merge_name] = merged_profile
        sub = lab.loc[merge_set]
        dom = sub.groupby("dominant")["n_cells"].sum().idxmax()
        lab = lab.drop(index=merge_set)
        lab.loc[args.merge_name] = {
            "n_cells": int(w.sum()),
            "dominant": dom,
            "frac": float(sub["frac"].mean()),
        }
        if comp is not None and all(c in comp.index for c in merge_set):
            merged_counts = comp.loc[merge_set].sum(axis=0)
            comp = comp.drop(index=merge_set)
            comp.loc[args.merge_name] = merged_counts
        print(f"Merged {merge_set} -> {args.merge_name} "
              f"(n={int(w.sum())} cells)")

    # Choose the outgroup the tree is rooted on:
    #   --root-on <leaf>  -> root on a real cluster (e.g. Merged_Normal)
    #   otherwise          -> add a synthetic all-neutral diploid outgroup
    global ROOT_NAME
    if args.root_on:
        if args.root_on not in prof.index:
            raise ValueError(f"--root-on leaf not found: {args.root_on}")
        ROOT_NAME = args.root_on
        names = list(prof.index)
    else:
        prof.loc[ROOT_NAME] = 0.0
        names = [c for c in prof.index if c != ROOT_NAME] + [ROOT_NAME]
    chrm = chrom_matrix(prof)

    # distance on full bin-level profiles
    P = prof.loc[names].values
    diff = P[:, None, :] - P[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=2))

    adj = neighbor_joining(dist, names)
    root, children, parent, blen = root_on_leaf(adj, ROOT_NAME)

    if args.leaf_states:
        # explicit per-leaf -1/0/1 states per region (e.g. penetrance calls);
        # events are placed by the same Sankoff parsimony as the rest.
        lsdf = pd.read_csv(args.leaf_states, sep="\t", index_col=0)
        lsdf.columns = [str(c) for c in lsdf.columns]
        chrom_order = [str(r) for r in lsdf.index]
        leaf_disc = {n: {r: (int(lsdf.loc[r, n]) if n in lsdf.columns else 0)
                         for r in chrom_order} for n in names}
        point_active = {}
        point_names = []
        char_order = chrom_order
    else:
        # per-region discretisation per leaf (whole chromosomes, or arms)
        if REGIONS is not None:
            chrom_order = [r[0] for r in REGIONS]
        else:
            chrom_order = [c for c in sorted(CHR_POS, key=lambda c: int(c[3:]))
                           if c not in EXCLUDE_CHROMS]

        def discretize(vec):
            return {ch: (1 if vec[ch] > args.gain_thr
                         else (-1 if vec[ch] < args.loss_thr else 0))
                    for ch in chrom_order}

        leaf_disc = {n: discretize(chrm.loc[n]) for n in names}

        # point regions (focal amplicons + extra arms), each with its own threshold
        point_active = {}
        if REGIONS is None:
            for k, (s, e, dlab) in FOCAL.items():
                point_active[k] = (s, e, dlab, args.focal_gain_thr, args.focal_loss_thr)
            for k, (s, e, dlab, gthr, lthr) in EXTRA.items():
                point_active[k] = (s, e, dlab,
                                   args.extra_gain_thr if gthr is None else gthr,
                                   args.extra_loss_thr if lthr is None else lthr)
        point_names = list(point_active.keys())
        for fname, (s, e, _lab, gthr, lthr) in point_active.items():
            fm = prof.iloc[:, s:e + 1].mean(axis=1)
            for n in names:
                v = float(fm[n])
                leaf_disc[n][fname] = (1 if v > gthr else (-1 if v < lthr else 0))

        char_order = chrom_order + point_names

    # Sankoff per character -> ancestral discrete states for every node
    node_state = {}
    all_nodes = list(blen.keys())
    for n in all_nodes:
        node_state[n] = {}
    for ch in char_order:
        ls = {n: leaf_disc[n][ch] for n in names}
        assign = sankoff(root, children, ls)
        for n, s in assign.items():
            node_state[n][ch] = s

    # ---- helpers to place a single focal event at the carriers' MRCA --------
    real_leaves = [n for n in all_nodes if not children[n] and n in lab.index]

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
        m = max(common, key=lambda c: depth[c])
        # if the MRCA is the (outgroup) root, descend to its tumour-side child
        if parent[m] is None:
            cset = set().union(*[set(ancestors(n)) for n in nodes])
            kids = [k for k in children[m] if k in cset]
            if kids:
                m = max(kids, key=lambda c: depth[c])
        return m

    # branch events (parent -> child): whole chromosomes only here
    events = {n: ([], []) for n in all_nodes}
    for child in all_nodes:
        par = parent[child]
        if par is None:
            continue
        for ch in chrom_order:
            d = node_state[child][ch] - node_state[par][ch]
            if d > 0:
                events[child][0].append(ch)
            elif d < 0:
                events[child][1].append(ch)

    # point-region events: place each ONCE, at the MRCA of its carrier clones,
    # labelled (e.g. +11 for CCND1, +1q / -8p for arms).
    for fname, (s, e, clab, gthr, lthr) in point_active.items():
        for direction, slot in (("gain", 0), ("loss", 1)):
            want = 1 if direction == "gain" else -1
            carriers = [lf for lf in real_leaves if leaf_disc[lf][fname] == want]
            if not carriers:
                continue
            m = mrca(carriers)
            events[m][slot].append(clab)

    # descendant tips per node (lets external tools map a branch to a clade/MRCA)
    def descendant_tips(node):
        if not children[node]:
            return [node]
        out = []
        for c in children[node]:
            out.extend(descendant_tips(c))
        return out

    # ------------------------------------------------------------- outputs
    # branch events table
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
    ev_df.to_csv(outdir / f"cnv_phylo_branch_events{args.suffix}.tsv", sep="\t", index=False)

    st_df = pd.DataFrame(node_state).T[char_order]
    st_df.index.name = "node"
    st_df.to_csv(outdir / f"cnv_phylo_chrom_states{args.suffix}.tsv", sep="\t")

    # post-merge cell-state composition (so ggtree pies match the tree's tips)
    if comp is not None:
        comp.to_csv(outdir / f"clone_state_composition{args.suffix}.tsv", sep="\t")

    # newick
    def newick(node):
        if not children[node]:
            nm = node.replace(",", "_").replace(":", "_").replace("(", "").replace(")", "")
            return f"{nm}:{blen.get(node,0):.5f}"
        inner = ",".join(newick(c) for c in children[node])
        return f"({inner}):{blen.get(node,0):.5f}"

    (outdir / f"cnv_phylo_tree{args.suffix}.newick").write_text(newick(root) + ";")

    # ------------------------------------------------------------- figure
    x, y, leaves = layout(root, children, blen)
    fig, ax = plt.subplots(figsize=(13, 8.5))

    def leaf_label(n):
        if n not in lab.index:          # synthetic diploid outgroup
            return n
        d = lab.loc[n]
        prefix = "clone " if n.isdigit() else ""
        tag = "  [root]" if n == ROOT_NAME else ""
        return f"{prefix}{n}  ({d['dominant']}, n={int(d['n_cells'])}){tag}"

    # draw edges (rectangular cladogram style)
    for child in all_nodes:
        par = parent[child]
        if par is None:
            continue
        xp, xc = x[par], x[child]
        yc = y[child]
        ax.plot([xp, xp], [y[par], yc], color="0.3", lw=1.3, zorder=1)  # vertical
        ax.plot([xp, xc], [yc, yc], color="0.3", lw=1.3, zorder=1)      # horizontal

        # branch CN-change label at midpoint of the horizontal segment.
        # display by chromosome number; de-duplicate so a focal + whole-chr hit
        # on the same branch shows the number only once.
        g, l = events.get(child, ([], []))
        def uniq(seq):
            out = []
            for c in seq:
                t = c.replace("chr", "")
                if t not in out:
                    out.append(t)
            return out
        gd, ld = uniq(g), uniq(l)
        parts = []
        if gd:
            parts.append("+" + " +".join(gd))
        if ld:
            parts.append("-" + " -".join(ld))
        if parts:
            txt = "  ".join(parts)
            color = "#b2182b" if gd and not ld else ("#2166ac" if ld and not gd else "#6a3d9a")
            xm = (xp + xc) / 2.0
            ax.text(xm, yc + 0.16, txt, fontsize=7.2, color=color,
                    ha="center", va="bottom", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white",
                              ec=color, lw=0.5, alpha=0.85))

    # finalise axis limits BEFORE placing pies (inset positions use the
    # data->axes transform, which depends on the limits)
    ax.set_xlim(-max(x.values()) * 0.02, max(x.values()) * 1.30)
    ax.set_ylim(-0.8, len(leaves) - 0.2)

    # state colour palette for the tip pies
    state_levels, state_pal = [], {}
    if comp is not None:
        state_levels = list(comp.columns)
        base = (plt.get_cmap("tab10").colors + plt.get_cmap("Set2").colors)
        state_pal = {s: base[i % len(base)] for i, s in enumerate(state_levels)}

    def tip_text(n):
        if n not in lab.index:
            return n
        prefix = "clone " if n.isdigit() else ""
        tag = "  [root]" if n == ROOT_NAME else ""
        return f"{prefix}{n}  (n={int(lab.loc[n, 'n_cells'])}){tag}"

    inv = ax.transAxes.inverted()
    PIE_W = 0.040           # pie size, axes fraction
    PIE_DX = 0.010          # gap between node and pie

    # nodes + tips (pie + short label)
    for n in all_nodes:
        if not children[n]:
            is_root = n == ROOT_NAME
            ax.scatter([x[n]], [y[n]], s=20 if not is_root else 45,
                       color="#1b7837" if is_root else "#333333", zorder=4)
            # node position in axes-fraction
            xa, ya = inv.transform(ax.transData.transform((x[n], y[n])))
            if comp is not None and n in comp.index and comp.loc[n].sum() > 0:
                frac = comp.loc[n, state_levels].astype(float).values
                axp = ax.inset_axes([xa + PIE_DX, ya - PIE_W / 2, PIE_W, PIE_W])
                axp.pie(frac, colors=[state_pal[s] for s in state_levels],
                        radius=1.0, counterclock=False, startangle=90,
                        wedgeprops=dict(linewidth=0.3, edgecolor="white"))
                axp.set_aspect("equal")
                tx = xa + PIE_DX + PIE_W + 0.006
            else:
                tx = xa + PIE_DX
            ax.text(tx, ya, tip_text(n), transform=ax.transAxes,
                    fontsize=9, va="center", ha="left",
                    fontweight="bold" if is_root else "normal")
        else:
            ax.scatter([x[n]], [y[n]], s=12, color="0.5", zorder=2)

    ax.set_title("inferCNV phylogenetic tree of CNV clones\n"
                 "(branches labelled with chromosome- and focal-level copy-number changes)",
                 fontsize=13)
    root_desc = ROOT_NAME if args.root_on else "diploid"
    ax.set_xlabel(f"CNV distance from {root_desc} root (Euclidean on inferCNV profile)")
    ax.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)

    import matplotlib.patches as mpatches
    handles = [
        mpatches.Patch(color="#b2182b", label="gain (+chr)"),
        mpatches.Patch(color="#2166ac", label="loss (-chr)"),
        mpatches.Patch(color="#6a3d9a", label="gain & loss"),
    ]
    leg1 = ax.legend(handles=handles, loc="lower right", fontsize=8,
                     frameon=False, title="branch CN change")
    ax.add_artist(leg1)
    if state_levels:
        shandles = [mpatches.Patch(color=state_pal[s], label=s)
                    for s in state_levels]
        ax.legend(handles=shandles, loc="upper left", fontsize=7.5,
                  frameon=False, title="cell state (tip pies)")

    # note: no tight_layout() — it repositions the axes and orphans the pie insets
    fig.savefig(outdir / f"cnv_phylo_tree{args.suffix}.png", dpi=240,
                bbox_inches="tight")
    fig.savefig(outdir / f"cnv_phylo_tree{args.suffix}.pdf", bbox_inches="tight")
    plt.close(fig)

    print("Wrote tree to", outdir / f"cnv_phylo_tree{args.suffix}.png")
    print("Leaves:", len(leaves), "| branches with events:", len(ev_df))
    print(ev_df.to_string(index=False))


if __name__ == "__main__":
    main()
