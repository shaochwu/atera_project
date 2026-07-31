#!/usr/bin/env python3
"""Merge tumour clones 3 and 4 into 'Clone_3-4' and rebuild the MP clone-tree inputs."""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

D = Path("results/infercnv_run")
BI = D / "biological_insight"
THR, FRAC = 0.02, 0.50
MERGE = ["3", "4"]
NAME = "Clone_3-4"
# truncal amplicons forced present in every tumour clone (matches the original
# run's --truncal set; the base branch of the tree carries +17 +12 +1q +8q +14)
TRUNCAL = ["17", "12", "1q", "8q", "14"]

a = ad.read_h5ad(D / "infercnv_atera.h5ad")
X = a.obsm["X_cnv"]
X = np.asarray(X.todense()) if not isinstance(X, np.ndarray) else X
g = a.obs["clone_grp"].astype(str).values.copy()
ident = a.obs["IDENT"].astype(str).values

# further-merge 3,4 on top of the existing grouping (keeps 0_1_5_9 intact)
g[np.isin(g, MERGE)] = NAME
tumor = [c for c in ["2", "8", "0_1_5_9", NAME]]          # tumour leaves (2 kept basal)
leaves = tumor + ["Merged_Normal"]
print("cells per leaf:", {c: int((g == c).sum()) for c in leaves},
      "| dropped:", int((g == "nan").sum()))

# ---- mean CNV profile per leaf (MP-tree tie-break) -------------------------
pd.DataFrame({c: X[g == c].mean(0) for c in leaves}).T.rename_axis("clone").to_csv(
    BI / "02_clone_mean_cnv_profile_merge34.tsv", sep="\t")

# ---- composition + labels --------------------------------------------------
states = ["Merged Normal", "Luminal DCIS1", "Luminal DCIS2",
          "Luminal DCIS3", "Luminal DCIS4"]
NORMAL_MAP = {"Ductal Luminal": "Merged Normal", "Normal Luminal": "Merged Normal"}
comp, lab = {}, []
for c in leaves:
    m = (g == c)
    vc = pd.Series(ident[m]).replace(NORMAL_MAP).value_counts()
    comp[c] = {s: int(vc.get(s, 0)) for s in states}
    lab.append(dict(clone=c, n_cells=int(m.sum()), dominant=vc.index[0],
                    frac=round(vc.iloc[0] / m.sum(), 2)))
pd.DataFrame(comp).T[states].to_csv(BI / "clone_state_composition_merge34.tsv", sep="\t")
pd.DataFrame(lab).set_index("clone").to_csv(BI / "cluster_leaf_labels_merge34.tsv", sep="\t")

# ---- >=50% chr-event leaf-states (same reg dict/rule as build_clone_bio) ----
reg = {"17": (1202, 1298, "+"), "12": (985, 998, "+"), "1q": (93, 177, "+"),
       "13": (1026, 1046, "-"), "11": (901, 904, "+"), "8q": (678, 712, "+"),
       "14": (1046, 1062, "+"), "8p": (660, 678, "-"), "4": (372, 432, "+"),
       "18": (1298, 1314, "+"), "5": (432, 504, "+-"), "15": (1094, 1137, "+"),
       "20": (1435, 1476, "+")}
masks = {c: (g == c) for c in tumor}
M = pd.DataFrame(0, index=list(reg), columns=leaves)
for r, (s, e, d) in reg.items():
    rm = X[:, s:e].mean(1)
    for c in tumor:
        if "+" in d and (rm[masks[c]] > THR).mean() >= FRAC:
            M.loc[r, c] = 1
        elif "-" in d and (rm[masks[c]] < -THR).mean() >= FRAC:
            M.loc[r, c] = -1
# force chosen regions truncal (present in every tumour clone)
for r in TRUNCAL:
    if r in M.index:
        M.loc[r, [c for c in leaves if c != "Merged_Normal"]] = 1
M = M.loc[(M != 0).any(axis=1)]
M.index.name = "region"
M.to_csv(BI / "clone_leaf_states_50_merge34.tsv", sep="\t")
print("\n=== merged 50% leaf-states ===\n" + M.to_string())
