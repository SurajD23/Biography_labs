"""Turn results/per_condition.csv into a summary (markdown + figures) focused on metric disagreement."""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "results"
METRICS = ["pcc_all", "pcc_top50", "pcc_all_delta", "pcc_top50_delta", "deg_overlap50", "edist", "sinkhorn"]
HIGHER_BETTER = {"pcc_all": True, "pcc_top50": True, "pcc_all_delta": True, "pcc_top50_delta": True,
                 "deg_overlap50": True, "edist": False, "sinkhorn": False, "mse_all": False, "mse_top50": False}
MODELS = ["no_change", "global_mean_shift", "additive", "additive_ridge"]

df = pd.read_csv(f"{OUT}/per_condition.csv")
meta = json.load(open(f"{OUT}/meta.json"))
lines = []
P = lines.append

P("# Results: Norman 2019 doubles, mean-shift vs additive\n")
P(f"{meta['n_cells']} cells x {meta['n_genes']} HVGs; {meta['n_singles']} singles, {meta['n_doubles']} doubles, "
  f"{meta['n_test_doubles']} test doubles (both singles present), {meta['args']['n_folds']}-fold CV over doubles.\n")

# ---- 1. mean +- sd per model/metric
agg = df.groupby("model")[METRICS + ["mse_all", "mse_top50"]].agg(["mean", "std"])
P("## 1. Mean (sd) per model\n")
tbl = pd.DataFrame({m: [f"{agg.loc[mod,(m,'mean')]:.3f} ({agg.loc[mod,(m,'std')]:.3f})" for mod in MODELS]
                    for m in METRICS + ["mse_all"]}, index=MODELS)
P(tbl.to_markdown()); P("")

# ---- 2. model ranking under each metric
P("## 2. Model ranking under each metric (best -> worst, by mean)\n")
P("| metric | ranking |"); P("|---|---|")
rank_tbl = {}
for m in METRICS:
    s = agg[(m, "mean")].loc[MODELS]
    order = s.sort_values(ascending=not HIGHER_BETTER[m]).index.tolist()
    rank_tbl[m] = order
    P(f"| {m} | {' > '.join(order)} |")
P("")

# ---- 3. per-condition win rate additive vs mean shift
P("## 3. Per-condition win rate: additive beats global_mean_shift\n")
piv = {m: df.pivot(index="condition", columns="model", values=m) for m in METRICS}
P("| metric | additive wins | ties | n |"); P("|---|---|---|---|")
for m in METRICS:
    a, b = piv[m]["additive"], piv[m]["global_mean_shift"]
    better = (a > b) if HIGHER_BETTER[m] else (a < b)
    ties = np.isclose(a, b)
    P(f"| {m} | {better.mean():.1%} | {ties.mean():.1%} | {len(a)} |")
P("")

# ---- 4. how do the metrics rank *conditions*?  Spearman between metrics for the additive model
P("## 4. Do the metrics agree on which conditions are predicted well? (Spearman rho across conditions, additive model)\n")
sub = df[df.model == "additive"].set_index("condition")[METRICS].copy()
# orient everything so that larger = better before correlating
for m in METRICS:
    if not HIGHER_BETTER[m]:
        sub[m] = -sub[m]
rho = sub.corr(method="spearman")
P(rho.round(2).to_markdown()); P("")

# ---- 5. effect-size confound: pcc_all vs how strong the perturbation is
sub2 = df[df.model == "additive"].copy()
r_all = stats.spearmanr(sub2["edist_ctrl_vs_true"], sub2["pcc_all"])[0]
r_del = stats.spearmanr(sub2["edist_ctrl_vs_true"], sub2["pcc_all_delta"])[0]
r_ed = stats.spearmanr(sub2["edist_ctrl_vs_true"], sub2["edist"])[0]
P("## 5. Confounding by perturbation strength (E-distance of control vs. true cells)\n")
P(f"Spearman rho(strength, metric) for the additive model: pcc_all = {r_all:.2f}, pcc_all_delta = {r_del:.2f}, edist = {r_ed:.2f}\n")

# ---- 5b. DEG overlap relative to its ceiling
ceil = df[df.model == "additive"].set_index("condition")["deg_overlap50_ceiling"]
P("## 5b. DEG overlap is capped: ceiling for *any* mean-shift model\n")
P(f"Plugging the true mean into the predicted-DEG ranking recovers only {ceil.mean():.2f} (median {ceil.median():.2f}, "
  f"min {ceil.min():.2f}) of the Wilcoxon top-50, so deg_overlap50 values should be read against that cap. "
  f"Fraction of the ceiling reached: additive {(piv['deg_overlap50']['additive']/ceil.replace(0,np.nan)).mean():.2f}, "
  f"global_mean_shift {(piv['deg_overlap50']['global_mean_shift']/ceil.replace(0,np.nan)).mean():.2f}.\n")

# ---- 6. no_change vs real models under pcc_all
nc = df[df.model == "no_change"]
P("## 6. The no-change reference\n")
P(f"Predicting the control mean for every double gets pcc_all = {nc.pcc_all.mean():.3f} (min {nc.pcc_all.min():.3f}) "
  f"and pcc_top50 = {nc.pcc_top50.mean():.3f}, but deg_overlap50 is undefined/0 and its deltas are all zero.\n")

# ---- 7. biggest disagreements: conditions ranked best by pcc_all but worst by edist (additive)
P("## 7. Largest rank disagreements for the additive model (pcc_all rank vs edist rank; 1 = best)\n")
r1 = sub2.set_index("condition")["pcc_all"].rank(ascending=False)
r2 = sub2.set_index("condition")["edist"].rank(ascending=True)
r3 = sub2.set_index("condition")["deg_overlap50"].rank(ascending=False)
dis = pd.DataFrame({"rank_pcc_all": r1, "rank_edist": r2, "rank_deg_overlap": r3,
                    "pcc_all": sub2.set_index("condition")["pcc_all"], "pcc_all_delta": sub2.set_index("condition")["pcc_all_delta"],
                    "deg_overlap50": sub2.set_index("condition")["deg_overlap50"], "edist": sub2.set_index("condition")["edist"],
                    "edist_ctrl_vs_true": sub2.set_index("condition")["edist_ctrl_vs_true"],
                    "n_cells": sub2.set_index("condition")["n_cells"]})
dis["gap"] = (dis.rank_pcc_all - dis.rank_edist).abs()
P(dis.sort_values("gap", ascending=False).head(10).round(3).to_markdown()); P("")

# ---- 8. conditions where the ordering of the two real models flips between metrics
P("## 8. Conditions where metrics disagree on which model is better\n")
flip = pd.DataFrame({m: ((piv[m]["additive"] > piv[m]["global_mean_shift"]) if HIGHER_BETTER[m]
                         else (piv[m]["additive"] < piv[m]["global_mean_shift"])) for m in METRICS})
n_consistent = (flip.nunique(axis=1) == 1).sum()
P(f"{n_consistent}/{len(flip)} conditions get the same winner under all {len(METRICS)} metrics; "
  f"{len(flip)-n_consistent} flip depending on the metric.\n")
pairs = [("pcc_all", "edist"), ("pcc_all", "deg_overlap50"), ("pcc_top50", "edist"), ("deg_overlap50", "sinkhorn"), ("edist", "sinkhorn")]
P("| metric A | metric B | % conditions where A and B pick different winners |"); P("|---|---|---|")
for a, b in pairs:
    P(f"| {a} | {b} | {(flip[a] != flip[b]).mean():.1%} |")
P("")

open(f"{OUT}/summary.md", "w").write("\n".join(lines))
print("\n".join(lines))

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, len(METRICS), figsize=(3.0 * len(METRICS), 3.4))
for ax, m in zip(axes, METRICS):
    data = [df[df.model == mod][m].values for mod in MODELS]
    ax.boxplot(data, tick_labels=["none", "shift", "add", "add+r"], showfliers=False)
    ax.set_title(m + ("  (↑)" if HIGHER_BETTER[m] else "  (↓)"), fontsize=10)
    ax.tick_params(axis="x", labelsize=8)
fig.suptitle("Per-condition metric distributions, Norman 2019 doubles")
fig.tight_layout(); fig.savefig(f"{OUT}/metrics_by_model.png", dpi=130); plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
s = df[df.model == "additive"]
axes[0].scatter(s.pcc_all, s.edist, s=14, alpha=.7); axes[0].set_xlabel("pcc_all (↑)"); axes[0].set_ylabel("edist (↓)")
axes[1].scatter(s.pcc_all_delta, s.edist, s=14, alpha=.7); axes[1].set_xlabel("pcc_all_delta (↑)"); axes[1].set_ylabel("edist (↓)")
axes[2].scatter(s.deg_overlap50, s.sinkhorn, s=14, alpha=.7); axes[2].set_xlabel("deg_overlap50 (↑)"); axes[2].set_ylabel("sinkhorn (↓)")
for ax in axes: ax.grid(alpha=.3)
fig.suptitle("Additive model: the same conditions scored by different metrics")
fig.tight_layout(); fig.savefig(f"{OUT}/metric_disagreement_scatter.png", dpi=130); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(rho.values, vmin=-1, vmax=1, cmap="coolwarm")
ax.set_xticks(range(len(METRICS))); ax.set_xticklabels(METRICS, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(METRICS))); ax.set_yticklabels(METRICS, fontsize=8)
for i in range(len(METRICS)):
    for j in range(len(METRICS)):
        ax.text(j, i, f"{rho.values[i,j]:.2f}", ha="center", va="center", fontsize=7)
fig.colorbar(im, ax=ax, shrink=.8); ax.set_title("Spearman rho between metrics across conditions (additive)")
fig.tight_layout(); fig.savefig(f"{OUT}/metric_rank_correlation.png", dpi=130); plt.close(fig)
print("figures written")
