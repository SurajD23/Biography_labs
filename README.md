# Perturb-seq baselines vs. metrics

A small, dependency-light benchmark asking one question: **when you evaluate perturbation-response
predictions, how much does the answer depend on the metric?**

Dataset: **Norman et al. 2019** (K562, CRISPRa, 105 single and 131 double perturbations, 111k cells),
scPerturb release `NormanWeissman2019_filtered.h5ad`. Two models that fit in a few lines each are scored on
all 131 double perturbations with seven metrics. Full numbers are in [`results/summary.md`](results/summary.md).

## Models

| name | prediction for double A+B |
|---|---|
| `no_change` (reference) | control mean |
| `global_mean_shift` | control mean + the average delta of every training condition (one shared shift) |
| `additive` | control mean + delta(A) + delta(B), deltas taken from the single-perturbation means |
| `additive_ridge` | same form, but per-target effects fit by ridge regression on singles **and** training doubles |

All models predict a *mean* log1p-normalised profile over 5,000 HVGs. For the distribution-level metrics the
prediction is turned into a cell population by shifting every control cell by the predicted delta.

## Metrics

| metric | what it compares | direction |
|---|---|---|
| `pcc_all` | Pearson r, predicted vs true mean, all 5,000 genes | higher better |
| `pcc_top50` | Pearson r on the true top-50 DEGs (Wilcoxon rank-sum vs control, tie-corrected) | higher better |
| `pcc_all_delta`, `pcc_top50_delta` | same, but on (x - control mean) | higher better |
| `deg_overlap50` | overlap between the true top-50 DEGs and the 50 genes with the largest predicted standardised effect \|delta\| / sqrt(ctrl var + median ctrl var) | higher better |
| `edist` | energy distance between predicted and true cells, 50-PC space, 500 cells each | lower better |
| `sinkhorn` | debiased Sinkhorn divergence, squared-Euclidean cost, eps = 0.5 x median squared distance, same space | lower better |

A note on the DEG metric, because the first version was wrong: with Welch's t the top-50 of an 86-cell condition
was pseudogenes and lncRNAs with |t| ~ 20 (all-zero in the perturbed group, so near-zero standard error), and even
the oracle overlap was 0. True DEGs therefore use a Wilcoxon z. A mean-level model has no cells, so its DEGs are
ranked by a standardised effect; the variance floor was chosen because it maximises the oracle ceiling (true mean
plugged in) over a sample of conditions. That ceiling is reported per condition as `deg_overlap50_ceiling`.
Running the Wilcoxon test on "control cells + predicted delta" is degenerate (a constant shift of a near-point-mass
gene separates perfectly for any delta), so it was not used.

## Protocol

* genes detected in < 20 cells dropped; library-size normalise to 1e4, log1p; top-5,000 Seurat-flavour HVGs
* conditions with < 30 cells dropped (none of the doubles are affected)
* test set = every double whose two singles are present (all 131 of them). 5-fold CV over doubles: each fold's
  doubles are held out, the other doubles plus all singles are the training conditions
* PCA (50 components) is fit on a random 60k-cell subset of all cells; distribution metrics use 500 cells per
  population

## Reproduce

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
mkdir -p data && curl -L -o data/NormanWeissman2019_filtered.h5ad \
  "https://zenodo.org/records/13350497/files/NormanWeissman2019_filtered.h5ad?download=1"
.venv/bin/python run.py        # ~4 min preprocessing (cached), ~15 min scoring
.venv/bin/python analyze.py    # writes results/summary.md and results/*.png
```

Code: [`perturbbench/data.py`](perturbbench/data.py) (streaming preprocessing),
[`perturbbench/models.py`](perturbbench/models.py), [`perturbbench/metrics.py`](perturbbench/metrics.py),
[`run.py`](run.py), [`analyze.py`](analyze.py).

## Findings

Numbers are means over the 131 doubles; per-condition values are in `results/per_condition.csv`.

| model | pcc_all | pcc_top50 | pcc_all_delta | deg_overlap50 | edist | sinkhorn |
|---|---|---|---|---|---|---|
| no_change | 0.984 | 0.916 | n/a | 0.00 | 2.45 | 40.4 |
| global_mean_shift | 0.989 | 0.942 | 0.57 | 0.24 | 1.80 | 32.1 |
| additive | 0.996 | 0.990 | 0.90 | 0.60 | 0.74 | 20.1 |
| additive_ridge | 0.997 | 0.990 | 0.91 | 0.60 | 0.58 | 18.5 |

**1. On the *model* leaderboard the metrics agree; on how much daylight there is, they do not.**
Every metric orders no_change < global_mean_shift < additive (the two additive variants swap places under
`pcc_top50` and `deg_overlap50`). But raw-expression Pearson compresses the whole field into a sliver:
predicting *the control mean for every perturbation* already scores pcc_all = 0.984 and pcc_top50 = 0.916,
and the best model is 0.997. The same gap is 0.57 vs 0.90 on delta-PCC, 0.24 vs 0.60 on DEG overlap and
a 2.4x difference in energy distance. Raw PCC mostly measures how correlated any K562 cell is with any other
K562 cell, not what the model learned.

**2. Per condition, the metrics disagree about *which conditions are predicted well* (section 4 of the summary).**
Spearman rank correlations across the 131 doubles, additive model:

* `pcc_all` vs `edist`: 0.97. Not because they measure the same thing, but because both are driven by
  perturbation strength: pcc_all correlates -0.44 with the control-vs-true energy distance, edist +0.47.
  Weak perturbations look easy under both, no matter how much of the effect was captured.
* `pcc_all_delta` vs `edist`: 0.11; vs `sinkhorn`: -0.01. The delta-based correlation is essentially
  independent of the distributional metrics, and it has the *opposite* strength confound (+0.49): strong
  perturbations give a clean signal to correlate with.
* `deg_overlap50` sits in the middle (0.4 to 0.7 with everything) and is capped: even the true mean, plugged
  into the predicted-DEG ranking, only recovers 0.78 of the Wilcoxon top-50 on average (min 0.10). The additive
  model reaches 76% of that ceiling, mean shift 31%.
* `edist` vs `sinkhorn`: 0.87. The two distributional metrics mostly agree with each other.

Concrete cases from section 7 of the summary:

* **KIF18B_KIF2C** (96 cells, almost no transcriptional effect): pcc_all = 0.998 and the 27th-best energy
  distance, but delta-PCC = 0.64 and DEG overlap = 0.04 (130th of 131). The model predicted roughly nothing,
  which was roughly right in absolute terms and useless in relative terms.
* **ZC3HAV1_CEBPA** (86 cells, strong effect): delta-PCC = 0.96 and DEG overlap 0.68 say the direction of the
  response is captured well; energy distance ranks it 122nd because the magnitude and spread of a strong
  CEBPA-driven shift are not reproduced by a mean shift of control cells.
* **ETS2_IKZF3**: 2nd-best DEG overlap in the whole set, 81st by energy distance.

**3. Which model wins a given condition also depends on the metric.**
Additive beats global mean shift on 88 to 99% of conditions depending on the metric. 18 of 131 doubles get a
different winner under at least one metric; `pcc_top50` and `edist` disagree on 11% of conditions, `pcc_all`
and `deg_overlap50` on 10%. The 15 conditions where mean shift beats additive under energy distance are
instructive: all are weaker-than-median perturbations and 9 involve CEBPE or CEBPB, whose doubles are known to
be sub-additive. On every one of them additive still has the far better delta-PCC (0.72 to 0.96 vs 0.10 to 0.76)
because the *direction* of the response is right, but it overshoots the *magnitude* (higher MSE than mean shift
in 12 of 15), and the distributional metrics charge for that while correlation, being scale-invariant, does not.

**4. Practical reading.** Report delta-space correlation or DEG overlap (with its ceiling) if the question is
"did the model learn the perturbation-specific response", and a distributional metric if the question is "are
the predicted cells where the real cells are". Raw PCC on all genes answers neither and should not be used to
rank perturbation models; a baseline that ignores the perturbation entirely gets 0.98.

Figures: `results/metrics_by_model.png`, `results/metric_disagreement_scatter.png`,
`results/metric_rank_correlation.png`.

## Caveats

* Only the Norman doubles are evaluated (the only setting where the additive model is defined). Mean shift on
  unseen *singles* would be the harder, more common benchmark setting.
* Distributional metrics use 500 cells per population in 50-PC space; conditions with fewer cells (min 54)
  use all their cells. The Sinkhorn regularisation (0.5 x median squared distance) was chosen on synthetic
  Gaussians so that two samples from one distribution score near zero; at 5% of the median the divergence
  collapsed onto the finite-sample W2 bias and could not separate identical from shifted distributions.
* Mean-level models produce no variance; the "shift the control cells" construction used for edist, sinkhorn
  and DEG calling is one choice, and it penalises any model for responses that change spread rather than mean.
* No cell-cycle or batch (gemgroup) correction was applied.
