"""Two deliberately simple predictors of the post-perturbation mean expression profile.

Both work in log1p-normalised gene space and predict a *mean* profile per condition.

(1) GlobalMeanShift   pred(c) = ctrl_mean + mean_{t in train conditions} (mu_t - ctrl_mean)
    i.e. one shared delta for every condition, regardless of what was perturbed.

(2) LinearAdditive    pred(A_B) = ctrl_mean + (mu_A - ctrl_mean) + (mu_B - ctrl_mean)
    the classic additive (no-interaction) model for combinatorial perturbations.
    `fit_ridge=True` instead fits  delta_c = sum_g beta_g * 1[g in c]  by ridge regression on
    all training conditions (singles *and* training doubles); with only singles in training the
    two are identical up to shrinkage.

(0) NoChange          pred(c) = ctrl_mean  -- a reference point, not a submission.
"""
from __future__ import annotations

import numpy as np
from .data import condition_targets


class NoChange:
    name = "no_change"

    def fit(self, cond_means: dict[str, np.ndarray], ctrl_mean: np.ndarray, train_conds: list[str]):
        self.ctrl_mean = ctrl_mean
        return self

    def predict(self, cond: str) -> np.ndarray:
        return self.ctrl_mean.copy()


class GlobalMeanShift:
    name = "global_mean_shift"

    def fit(self, cond_means, ctrl_mean, train_conds):
        self.ctrl_mean = ctrl_mean
        deltas = np.stack([cond_means[c] - ctrl_mean for c in train_conds])
        self.delta = deltas.mean(axis=0)
        return self

    def predict(self, cond: str) -> np.ndarray:
        return self.ctrl_mean + self.delta


class LinearAdditive:
    def __init__(self, fit_ridge: bool = False, alpha: float = 1.0):
        self.fit_ridge = fit_ridge
        self.alpha = alpha
        self.name = "additive_ridge" if fit_ridge else "additive"

    def fit(self, cond_means, ctrl_mean, train_conds):
        self.ctrl_mean = ctrl_mean
        targets = sorted({g for c in train_conds for g in condition_targets(c)})
        self.idx = {g: i for i, g in enumerate(targets)}
        if not self.fit_ridge:
            # closed form from single-perturbation means only
            self.beta = np.zeros((len(targets), ctrl_mean.shape[0]), dtype=np.float32)
            self.has = np.zeros(len(targets), dtype=bool)
            for c in train_conds:
                t = condition_targets(c)
                if len(t) == 1 and t[0] in self.idx:
                    self.beta[self.idx[t[0]]] = cond_means[c] - ctrl_mean
                    self.has[self.idx[t[0]]] = True
        else:
            D = np.zeros((len(train_conds), len(targets)), dtype=np.float32)
            Y = np.zeros((len(train_conds), ctrl_mean.shape[0]), dtype=np.float32)
            for r, c in enumerate(train_conds):
                for g in condition_targets(c):
                    D[r, self.idx[g]] = 1.0
                Y[r] = cond_means[c] - ctrl_mean
            A = D.T @ D + self.alpha * np.eye(D.shape[1], dtype=np.float32)
            self.beta = np.linalg.solve(A, D.T @ Y)
            self.has = np.ones(len(targets), dtype=bool)
        return self

    def can_predict(self, cond: str) -> bool:
        return all(g in self.idx and self.has[self.idx[g]] for g in condition_targets(cond))

    def predict(self, cond: str) -> np.ndarray:
        d = np.zeros_like(self.ctrl_mean)
        for g in condition_targets(cond):
            d += self.beta[self.idx[g]]
        return self.ctrl_mean + d
