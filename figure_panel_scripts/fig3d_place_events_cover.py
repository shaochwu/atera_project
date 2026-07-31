#!/usr/bin/env python3
"""Place CN events on a tree by maximal carrier-clade cover."""
import argparse
import pandas as pd


def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = [0]

    def node():
        children = []
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                children.append(node())
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                if s[pos[0]] == ")":
                    pos[0] += 1
                    break
        name = ""
        while pos[0] < len(s) and s[pos[0]] not in ",():":
            name += s[pos[0]]
            pos[0] += 1
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(s) and s[pos[0]] not in ",()":
                pos[0] += 1
        return {"name": name, "children": children}

    return node()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--newick", required=True)
    ap.add_argument("--leaf-states", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    root = parse_newick(open(args.newick).read())

    # collect nodes: id -> (frozenset leaves, parent_leaves)
    nodes = []  # (leaves, parent_leaves)

    def leaves_of(n):
        if not n["children"]:
            return frozenset([n["name"]])
        return frozenset().union(*[leaves_of(c) for c in n["children"]])

    def walk(n, parent_leaves):
        lv = leaves_of(n)
        nodes.append((lv, parent_leaves))
        for c in n["children"]:
            walk(c, lv)

    walk(root, frozenset())

    ls = pd.read_csv(args.leaf_states, sep="\t", index_col=0)
    ls.columns = [str(c) for c in ls.columns]

    # node -> (gains[], losses[])
    ev = {}
    for region in ls.index:
        for sign, slot in ((1, "gains"), (-1, "losses")):
            carriers = frozenset(c for c in ls.columns if ls.loc[region, c] == sign)
            if not carriers:
                continue
            for lv, par in nodes:
                if lv and lv <= carriers and not (par <= carriers):
                    key = tuple(sorted(lv))
                    ev.setdefault(key, {"gains": [], "losses": []})[slot].append(str(region))

    rows = []
    for tips, gl in ev.items():
        rows.append({"branch_to": "n", "from": "cover",
                     "gains": ",".join(gl["gains"]), "losses": ",".join(gl["losses"]),
                     "descendant_tips": ",".join(tips)})
    pd.DataFrame(rows)[["branch_to", "from", "gains", "losses", "descendant_tips"]].to_csv(
        args.out, sep="\t", index=False)
    print(f"wrote {len(rows)} branch-event rows -> {args.out}")


if __name__ == "__main__":
    main()
