"""Benchmark global-mean-shift vs linear-additive on Norman 2019 double perturbations.

Protocol
--------
* Conditions with < MIN_CELLS cells are dropped.
* Test set = every double perturbation A_B whose singles A and B are both present (the "seen 2/2"
  setting; the only one where an additive model is defined). 5-fold CV over the doubles so every
  double is scored exactly once while the other folds' doubles are available for training.
* Train set for a fold = all singles + the doubles in the other folds.  Global mean shift averages
  the deltas of all training conditions; additive uses the singles' deltas; additive_ridge fits
  per-target effects on singles + training doubles.
* Metrics are in perturbbench/metrics.py.  Distribution metrics are computed in a 50-PC space fit
  on all cells, on <= N_CELLS_DIST cells per population; the mean-level prediction is turned into a
  population by shifting the control cells by the predicted delta.
"""
from __future__ import annotations

import argparse, json, time
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

from perturbbench.data import load_norman, condition_targets, CTRL
from perturbbench.models import NoChange, GlobalMeanShift, LinearAdditive
from perturbbench import metrics as M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/NormanWeissman2019_filtered.h5ad")
    ap.add_argument("--out", default="results")
    ap.add_argument("--n-hvg", type=int, default=5000)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--n-folds", type=int, default=5)
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--n-cells-dist", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    import os, anndata as ad
    cache = f"{args.out}/norman_hvg{args.n_hvg}_lognorm.h5ad"
    if os.path.exists(cache):
        adata = ad.read_h5ad(cache)
    else:
        adata = load_norman(args.data, n_hvg=args.n_hvg)
        adata.write_h5ad(cache)
    print(f"loaded {adata.shape} in {time.time()-t0:.0f}s")
    X = adata.X.tocsr()
    cond = adata.obs["condition"].values
    counts = pd.Series(cond).value_counts()
    keep_conds = set(counts[counts >= args.min_cells].index)
    idx = {c: np.where(cond == c)[0] for c in keep_conds}

    ctrl_idx = idx[CTRL]
    X_ctrl = X[ctrl_idx].toarray()
    ctrl_mean = X_ctrl.mean(0)
    ctrl_var = X_ctrl.var(0, ddof=1)

    singles = sorted(c for c in keep_conds if c != CTRL and len(condition_targets(c)) == 1)
    doubles = sorted(c for c in keep_conds if len(condition_targets(c)) == 2)
    single_set = set(singles)
    test_doubles = [d for d in doubles if all(g in single_set for g in condition_targets(d))]
    print(f"{len(singles)} singles, {len(doubles)} doubles, {len(test_doubles)} doubles with both singles present")

    cond_means = {c: np.asarray(X[idx[c]].mean(0)).ravel() for c in keep_conds if c != CTRL}

    # PCA space for distribution metrics
    sub = rng.choice(X.shape[0], size=min(60000, X.shape[0]), replace=False)
    pca = PCA(n_components=args.n_pcs, svd_solver="randomized", random_state=args.seed).fit(X[sub].toarray())
    ctrl_sub = rng.choice(len(ctrl_idx), size=min(args.n_cells_dist, len(ctrl_idx)), replace=False)
    X_ctrl_sub = X_ctrl[ctrl_sub]
    Z_ctrl_sub = pca.transform(X_ctrl_sub)
    # global Sinkhorn regularisation: 0.5 x median squared distance ctrl -> random perturbed cells
    pert_cells = rng.choice(np.where(cond != CTRL)[0], size=500, replace=False)
    from scipy.spatial.distance import cdist
    reg = 0.5 * float(np.median(cdist(Z_ctrl_sub, pca.transform(X[pert_cells].toarray()), "sqeuclidean")))
    print(f"sinkhorn reg = {reg:.3f}")

    models = [NoChange(), GlobalMeanShift(), LinearAdditive(fit_ridge=False), LinearAdditive(fit_ridge=True, alpha=1.0)]
    rows = []
    kf = KFold(n_splits=args.n_folds, shuffle=True, random_state=args.seed)
    test_arr = np.array(test_doubles)
    for fold, (tr, te) in enumerate(kf.split(test_arr)):
        held = set(test_arr[te])
        train_conds = singles + [d for d in doubles if d not in held]
        for m in models:
            m.fit(cond_means, ctrl_mean, train_conds)
        for c in test_arr[te]:
            X_true = X[idx[c]].toarray()
            true_mean = X_true.mean(0)
            deg = M.true_top_degs(X_true, X_ctrl, k=50)
            tsub = rng.choice(len(X_true), size=min(args.n_cells_dist, len(X_true)), replace=False)
            Z_true = pca.transform(X_true[tsub])
            ref = {"edist_ctrl_vs_true": M.energy_distance(Z_ctrl_sub, Z_true)}
            # attainable ceiling of deg_overlap50 for any mean-shift model: shift control by the TRUE delta
            odeg = M.pred_top_degs(true_mean, ctrl_mean, ctrl_var, 50)
            ref["deg_overlap50_ceiling"] = len(set(odeg.tolist()) & set(deg.tolist())) / 50
            for m in models:
                pred = m.predict(c)
                r = {"fold": fold, "condition": c, "model": m.name, "n_cells": len(X_true)}
                r.update(M.mean_level_metrics(pred, true_mean, ctrl_mean, ctrl_var, deg))
                P = pca.transform(M.shifted_population(X_ctrl_sub, pred, ctrl_mean))
                r.update(M.distribution_metrics(P, Z_true, reg))
                r.update(ref)
                rows.append(r)
        print(f"fold {fold} done ({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(f"{args.out}/per_condition.csv", index=False)
    meta = {"n_cells": int(X.shape[0]), "n_genes": int(X.shape[1]), "n_singles": len(singles),
            "n_doubles": len(doubles), "n_test_doubles": len(test_doubles), "sinkhorn_reg": reg,
            "args": vars(args)}
    json.dump(meta, open(f"{args.out}/meta.json", "w"), indent=2)
    print(df.groupby("model")[["pcc_all", "pcc_all_delta", "pcc_top50", "pcc_top50_delta", "deg_overlap50", "edist", "sinkhorn"]].mean().round(3))


if __name__ == "__main__":
    main()
