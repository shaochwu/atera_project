#!/usr/bin/env python3
"""Exhaustive maximum-parsimony tree from CN characters (Sankoff, outgroup-rooted)."""
import argparse
import numpy as np
import pandas as pd

STATES = np.array([-1, 0, 1])
C = np.abs(STATES[:, None] - STATES[None, :]).astype(float)
INF = 1e9


def sk(x):
    return str(x)


def enumerate_unrooted(taxa):
    taxa = list(taxa)
    adj0 = {taxa[0]: {0}, taxa[1]: {0}, taxa[2]: {0},
            0: {taxa[0], taxa[1], taxa[2]}}
    trees = [adj0]
    for k in range(3, len(taxa)):
        leaf = taxa[k]
        new = []
        for t in trees:
            w = sum(1 for x in t if isinstance(x, int))   # next internal id
            edges = sorted({frozenset((u, v)) for u in t for v in t[u]},
                           key=lambda fs: tuple(sorted(map(sk, fs))))
            for ed in edges:
                u, v = sorted(ed, key=sk)
                t2 = {x: set(t[x]) for x in t}
                t2[u].discard(v); t2[v].discard(u)
                t2[w] = {u, v, leaf}; t2[u].add(w); t2[v].add(w); t2[leaf] = {w}
                new.append(t2)
        trees = new
    return trees


def root_traversal(adj, root):
    parent, children = {root: None}, {n: [] for n in adj}
    stack, seen = [root], set()
    while stack:
        n = stack.pop()
        if n in seen:
            continue
        seen.add(n)
        for nb in sorted(adj[n], key=sk):          # deterministic
            if nb != parent.get(n):
                parent[nb] = n
                children[n].append(nb)
                stack.append(nb)
    post = []

    def dfs(n):
        for c in children[n]:
            dfs(c)
        post.append(n)
    dfs(root)
    return parent, children, post


def root_on_outgroup(adj, og):
    xnode = sorted(adj[og], key=sk)[0]
    radj = {x: set(adj[x]) for x in adj}
    radj[og].discard(xnode); radj[xnode].discard(og)
    radj["ROOT"] = {og, xnode}; radj[og].add("ROOT"); radj[xnode].add("ROOT")
    return radj


def sankoff_score(adj, leafcode, nchar):
    root = min((x for x in adj if isinstance(x, int)))
    _, children, post = root_traversal(adj, root)
    cost = {}
    for n in post:
        if not children[n]:
            g = np.full((3, nchar), INF); g[leafcode[n], np.arange(nchar)] = 0.0
        else:
            g = np.zeros((3, nchar))
            for c in children[n]:
                cc = cost[c]
                g += np.stack([(C[s][:, None] + cc).min(0) for s in range(3)])
        cost[n] = g
    return float(cost[root].min(0).sum())


def profile_treelen(adj, og, prof):
    """Total Euclidean branch length under subtree-mean reconstruction."""
    radj = root_on_outgroup(adj, og)
    parent, children, post = root_traversal(radj, "ROOT")
    npf, cnt = {}, {}
    for n in post:
        if not children[n]:
            npf[n] = prof[n].astype(float); cnt[n] = 1
        else:
            s = np.zeros_like(next(iter(prof.values())), dtype=float); c = 0
            for ch in children[n]:
                s += npf[ch] * cnt[ch]; c += cnt[ch]
            npf[n] = s / c; cnt[n] = c
    return float(sum(np.linalg.norm(npf[n] - npf[parent[n]])
                     for n in parent if parent[n] is not None))


def newick_of(adj, og, leafcode, nchar):
    radj = root_on_outgroup(adj, og)
    parent, children, post = root_traversal(radj, "ROOT")
    cost = {}
    for n in post:
        if not children[n]:
            g = np.full((3, nchar), INF); g[leafcode[n], np.arange(nchar)] = 0.0
        else:
            g = np.zeros((3, nchar))
            for c in children[n]:
                cc = cost[c]
                g += np.stack([(C[s][:, None] + cc).min(0) for s in range(3)])
        cost[n] = g
    assign = {"ROOT": cost["ROOT"].argmin(0)}
    for n in post[::-1]:
        for c in children[n]:
            opts = np.stack([C[assign[n], t] + cost[c][t] for t in range(3)])
            assign[c] = opts.argmin(0)
    blen = {"ROOT": 0.0}
    for n in radj:
        if n != "ROOT":
            blen[n] = float(np.abs(STATES[assign[n]] - STATES[assign[parent[n]]]).sum())
    depth = {"ROOT": 0.0}
    for n in post[::-1]:
        for c in children[n]:
            depth[c] = depth[n] + blen[c]
    maxd = max(depth.values())
    if maxd > 0:
        f = 4.5 / maxd
        blen = {n: v * f for n, v in blen.items()}

    def nwk(n):
        kids = sorted(children[n], key=sk)
        if not kids:
            nm = str(n).replace(",", "_").replace(":", "_").replace("(", "").replace(")", "")
            return f"{nm}:{blen[n]:.3f}"
        return "(" + ",".join(nwk(c) for c in kids) + f"):{blen.get(n,0):.3f}"

    return nwk("ROOT") + ";"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", default="")
    ap.add_argument("--leaf-states", default="")
    ap.add_argument("--tiebreak-profiles", default="",
                    help="profiles TSV (taxon x bins) to break ties by min tree length")
    ap.add_argument("--root", required=True)
    ap.add_argument("--thr", type=float, default=0.02)
    ap.add_argument("--prefer-basal", default="",
                    help="among co-optimal trees, prefer ones where this taxon is the "
                         "basal tumour tip (sister to all other non-root taxa)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.leaf_states:
        ls = pd.read_csv(args.leaf_states, sep="\t", index_col=0)
        ls.columns = [str(c) for c in ls.columns]
        taxa = list(ls.columns); disc = ls.values.T.astype(int)
    else:
        prof = pd.read_csv(args.profiles, sep="\t", index_col=0)
        prof.index = prof.index.astype(str)
        taxa = list(prof.index)
        disc = np.where(prof.values > args.thr, 1, np.where(prof.values < -args.thr, -1, 0))
    keep = np.array([len(set(disc[:, j])) > 1 for j in range(disc.shape[1])])
    disc = disc[:, keep]
    nchar = disc.shape[1]
    code = {-1: 0, 0: 1, 1: 2}
    leafcode = {taxa[i]: np.array([code[v] for v in disc[i]]) for i in range(len(taxa))}
    print(f"{len(taxa)} taxa, {nchar} informative characters")

    trees = enumerate_unrooted(taxa)
    scores = [sankoff_score(t, leafcode, nchar) for t in trees]
    best = min(scores)
    opt = [t for t, s in zip(trees, scores) if abs(s - best) < 1e-6]
    print(f"min parsimony score = {best:.0f}  ({len(opt)} optimal tree(s))")

    # among co-optimal trees, optionally keep only those where a chosen taxon is the
    # basal tumour tip (sister to all other non-root taxa) -- a biological prior, not
    # a score change (these trees are equally parsimonious)
    if args.prefer_basal:
        def basal_tip(adj):
            radj = root_on_outgroup(adj, args.root)
            _, children, post = root_traversal(radj, "ROOT")
            L = {}
            for n in post:
                L[n] = (frozenset([n]) if not children[n]
                        else frozenset().union(*[L[c] for c in children[n]]))
            tumclade = [c for c in children["ROOT"] if args.root not in L[c]][0]
            return set(s for c in children[tumclade] if len(L[c]) == 1 for s in L[c])
        filt = [t for t in opt if args.prefer_basal in basal_tip(t)]
        if filt:
            opt = filt
            print(f"restricted to {len(opt)} optimal tree(s) with "
                  f"'{args.prefer_basal}' as the basal tumour tip")
        else:
            print(f"WARNING: no optimal tree has '{args.prefer_basal}' basal; keeping all")

    tb = None
    if args.tiebreak_profiles:
        pf = pd.read_csv(args.tiebreak_profiles, sep="\t", index_col=0)
        pf.index = pf.index.astype(str)
        tb = {t: pf.loc[t].values.astype(float) for t in taxa}

    # deterministic choice among optima: (profile tree length, newick string)
    def key(t):
        nwk = newick_of(t, args.root, leafcode, nchar)
        plen = profile_treelen(t, args.root, tb) if tb is not None else 0.0
        return (round(plen, 6), nwk)
    chosen = min(opt, key=key)
    if tb is not None:
        print(f"tie broken by CNV-profile tree length "
              f"(min = {profile_treelen(chosen, args.root, tb):.3f})")

    open(args.out, "w").write(newick_of(chosen, args.root, leafcode, nchar))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
