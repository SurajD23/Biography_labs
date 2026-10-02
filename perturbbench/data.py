"""Load and preprocess the scPerturb copy of Norman et al. 2019 (K562 CRISPRa Perturb-seq).

Download (~700 MB):
  https://zenodo.org/records/13350497/files/NormanWeissman2019_filtered.h5ad

Preprocessing is deliberately plain and done with numpy/scipy only, streaming over row chunks so
the 111k x 34k count matrix never has to be fully materialised:
  1. drop genes detected in < MIN_CELLS cells
  2. library-size normalise to 1e4 counts/cell, log1p
  3. keep the top N_HVG highly-variable genes (Seurat-style binned normalised dispersion,
     computed on the normalised (not logged) values, as scanpy flavor='seurat' does)
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import anndata as ad

CTRL = "control"
TARGET_SUM = 1e4


def _normalize(Xc: sp.csr_matrix) -> sp.csr_matrix:
    lib = np.asarray(Xc.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    return sp.csr_matrix(sp.diags((TARGET_SUM / lib).astype(np.float32)) @ Xc)


def _seurat_hvg_from_moments(mean: np.ndarray, var: np.ndarray, n_top: int, n_bins: int = 20) -> np.ndarray:
    mean = mean.copy(); mean[mean == 0] = 1e-12
    disp = var / mean
    disp[disp == 0] = np.nan
    disp = np.log(disp)
    mean = np.log1p(mean)
    bins = np.digitize(mean, np.linspace(mean.min(), mean.max() + 1e-9, n_bins + 1))
    norm = np.full_like(disp, np.nan)
    for b in np.unique(bins):
        m = bins == b
        d = disp[m]
        if np.sum(~np.isnan(d)) <= 1:
            norm[m] = 0.0
            continue
        mu, sd = np.nanmean(d), np.nanstd(d, ddof=1)
        norm[m] = (d - mu) / (sd if sd > 0 else 1.0)
    norm = np.nan_to_num(norm, nan=-np.inf)
    return np.sort(np.argsort(-norm)[:n_top])


def load_norman(path: str, n_hvg: int = 5000, min_cells: int = 20, chunk: int = 10000) -> ad.AnnData:
    backed = ad.read_h5ad(path, backed="r")
    pert = backed.obs["perturbation"].astype(str).values
    obs = backed.obs.copy()
    obs["condition"] = pert
    obs["n_targets"] = np.array([0 if p == CTRL else p.count("_") + 1 for p in pert])

    keep_genes = np.where(backed.var["ncells"].values >= min_cells)[0]
    n = backed.n_obs

    # pass 1: moments of normalised counts for HVG selection
    s1 = np.zeros(len(keep_genes)); s2 = np.zeros(len(keep_genes))
    for i in range(0, n, chunk):
        Xc = sp.csr_matrix(backed[i:i + chunk].to_memory().X)[:, keep_genes]
        Xn = _normalize(Xc)
        s1 += np.asarray(Xn.sum(axis=0)).ravel()
        s2 += np.asarray(Xn.multiply(Xn).sum(axis=0)).ravel()
    mean = s1 / n
    var = (s2 / n - mean**2) * n / (n - 1)
    hvg_local = _seurat_hvg_from_moments(mean, var, n_hvg)
    hvg = keep_genes[hvg_local]

    # pass 2: build log1p-normalised HVG matrix
    blocks = []
    for i in range(0, n, chunk):
        Xc = sp.csr_matrix(backed[i:i + chunk].to_memory().X)
        Xn = _normalize(Xc)[:, hvg]
        Xn.data = np.log1p(Xn.data)
        blocks.append(sp.csr_matrix(Xn, dtype=np.float32))
    X = sp.vstack(blocks).tocsr()
    var_df = backed.var.iloc[hvg].copy()
    out = ad.AnnData(X=X, obs=obs, var=var_df)
    out.var_names = backed.var_names[hvg]
    return out


def condition_targets(cond: str) -> tuple[str, ...]:
    return tuple(sorted(cond.split("_")))
