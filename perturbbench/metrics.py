"""Evaluation metrics for predicted post-perturbation expression.

Mean-level metrics (operate on a predicted mean profile vs. the true mean profile):
  * pcc_all          Pearson r over all genes, raw expression            (notoriously inflated)
  * pcc_all_delta    Pearson r over all genes, on (x - ctrl_mean)        (what actually discriminates)
  * pcc_top50        Pearson r on the true top-50 DEGs, raw expression
  * pcc_top50_delta  Pearson r on the true top-50 DEGs, on deltas
  * deg_overlap50    |top-50 predicted DEGs  ∩  top-50 true DEGs| / 50   (true: Wilcoxon z; predicted: standardised effect, see pred_top_degs)

Distribution-level metrics (operate on cell populations; the mean-level model is turned into a
population by shifting every control cell by the predicted delta):
  * edist            energy distance in PCA space (Peidli et al. 2024 / scPerturb E-distance)
  * sinkhorn         Sinkhorn divergence (entropic OT, debiased) in PCA space, via POT
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.spatial.distance import cdist
import ot


# ----------------------------------------------------------------------------- DEGs
def welch_t(Xa: np.ndarray, Xb: np.ndarray) -> np.ndarray:
    """Per-gene Welch t statistic.  NOT used for DEG calling: with < ~100 cells a gene that is all-zero in
    the perturbed group gets a near-zero standard error and a huge |t| (the top-50 of an 86-cell condition
    came out as pseudogenes/lncRNAs with |t| ~ 20 and an oracle overlap ceiling of 0)."""
    ma, mb = Xa.mean(0), Xb.mean(0)
    va, vb = Xa.var(0, ddof=1), Xb.var(0, ddof=1)
    se = np.sqrt(va / Xa.shape[0] + vb / Xb.shape[0])
    se[se == 0] = np.inf
    return (ma - mb) / se


def wilcoxon_z(Xa: np.ndarray, Xb: np.ndarray) -> np.ndarray:
    """Per-gene Wilcoxon rank-sum (Mann-Whitney) z with tie correction, group a vs group b.
    Equivalent to scanpy rank_genes_groups(method='wilcoxon', tie_correct=True)."""
    na, nb = Xa.shape[0], Xb.shape[0]
    N = na + nb
    pooled = np.vstack([Xa, Xb])
    ranks = stats.rankdata(pooled, axis=0)
    R1 = ranks[:na].sum(0)
    mu = na * (N + 1) / 2.0
    # tie correction: sum(t^3 - t) over tie groups, per gene
    srt = np.sort(pooled, axis=0)
    tie_term = np.zeros(pooled.shape[1])
    for g in range(pooled.shape[1]):
        col = srt[:, g]
        brk = np.flatnonzero(np.diff(col) != 0) + 1
        t = np.diff(np.concatenate([[0], brk, [N]]))
        tie_term[g] = np.sum(t.astype(float) ** 3 - t)
    var = na * nb / 12.0 * ((N + 1) - tie_term / (N * (N - 1)))
    sd = np.sqrt(np.maximum(var, 0))
    sd[sd == 0] = np.inf
    return (R1 - mu) / sd


def true_top_degs(X_true: np.ndarray, X_ctrl: np.ndarray, k: int = 50) -> np.ndarray:
    z = wilcoxon_z(X_true, X_ctrl)
    return np.argsort(-np.abs(z))[:k]


def pred_top_degs(pred_mean: np.ndarray, ctrl_mean: np.ndarray, ctrl_var: np.ndarray, k: int = 50) -> np.ndarray:
    """A mean-only model has no cells, so predicted DEGs are ranked by a standardised effect size
    |delta| / sqrt(ctrl_var + median(ctrl_var)).  The floor matters: without it, near-zero-variance genes
    dominate; with it the oracle (true mean plugged in) recovers ~0.8 of the Wilcoxon top-50 on average,
    the best of the rankings we tried (|delta| alone: 0.61; |delta|/sd: 0.72).  All-zero delta -> no DEGs.
    Running the Wilcoxon test itself on control cells + delta is degenerate: a constant shift of a
    near-point-mass gene gives perfect separation for any delta."""
    delta = pred_mean - ctrl_mean
    if not np.any(delta != 0):
        return np.array([], dtype=int)
    score = np.abs(delta) / np.sqrt(ctrl_var + np.median(ctrl_var))
    return np.argsort(-score)[:k]


# ----------------------------------------------------------------------------- Pearson
def pcc(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(stats.pearsonr(a, b)[0])


def mean_level_metrics(pred_mean, true_mean, ctrl_mean, ctrl_var, deg_idx, k=50) -> dict:
    out = {
        "pcc_all": pcc(pred_mean, true_mean),
        "pcc_all_delta": pcc(pred_mean - ctrl_mean, true_mean - ctrl_mean),
        "pcc_top50": pcc(pred_mean[deg_idx], true_mean[deg_idx]),
        "pcc_top50_delta": pcc(pred_mean[deg_idx] - ctrl_mean[deg_idx], true_mean[deg_idx] - ctrl_mean[deg_idx]),
        "mse_all": float(np.mean((pred_mean - true_mean) ** 2)),
        "mse_top50": float(np.mean((pred_mean[deg_idx] - true_mean[deg_idx]) ** 2)),
    }
    pdeg = pred_top_degs(pred_mean, ctrl_mean, ctrl_var, k)
    out["deg_overlap50"] = len(set(pdeg.tolist()) & set(deg_idx.tolist())) / k
    return out


# ----------------------------------------------------------------------------- distributional
def energy_distance(A: np.ndarray, B: np.ndarray) -> float:
    """E(A,B) = 2 E|a-b| - E|a-a'| - E|b-b'|  (Euclidean, includes self-pairs as scPerturb does)."""
    dab = cdist(A, B).mean()
    daa = cdist(A, A).mean()
    dbb = cdist(B, B).mean()
    return float(2 * dab - daa - dbb)


def _sinkhorn_cost(A: np.ndarray, B: np.ndarray, reg: float) -> float:
    M = cdist(A, B, metric="sqeuclidean")
    a = np.full(len(A), 1 / len(A)); b = np.full(len(B), 1 / len(B))
    return float(ot.sinkhorn2(a, b, M, reg, method="sinkhorn_log", numItermax=2000, stopThr=1e-7))


def sinkhorn_divergence(A: np.ndarray, B: np.ndarray, reg: float | None = None) -> float:
    """Debiased Sinkhorn divergence  S = OT_eps(A,B) - 0.5 OT_eps(A,A) - 0.5 OT_eps(B,B)
    with squared-Euclidean cost and entropic regularisation eps=reg (log-domain Sinkhorn).

    reg defaults to 0.5 x the median squared pairwise distance; with 500 cells in 50 PCs, smaller eps makes the
    divergence collapse onto the finite-sample W2 bias (two draws from the same distribution score ~70 while a
    0.3-sd shift scores ~73), 0.5 x median gives ~1.6 vs ~6."""
    if reg is None:
        reg = 0.5 * float(np.median(cdist(A[:500], B[:500], metric="sqeuclidean")))
    return _sinkhorn_cost(A, B, reg) - 0.5 * _sinkhorn_cost(A, A, reg) - 0.5 * _sinkhorn_cost(B, B, reg)


def shifted_population(X_ctrl: np.ndarray, pred_mean: np.ndarray, ctrl_mean: np.ndarray) -> np.ndarray:
    return X_ctrl + (pred_mean - ctrl_mean)[None, :]


def distribution_metrics(P: np.ndarray, T: np.ndarray, reg: float) -> dict:
    return {"edist": energy_distance(P, T), "sinkhorn": sinkhorn_divergence(P, T, reg)}
