# Results: Norman 2019 doubles, mean-shift vs additive

111445 cells x 5000 HVGs; 105 singles, 131 doubles, 131 test doubles (both singles present), 5-fold CV over doubles.

## 1. Mean (sd) per model

|                   | pcc_all       | pcc_top50     | pcc_all_delta   | pcc_top50_delta   | deg_overlap50   | edist         | sinkhorn        | mse_all       |
|:------------------|:--------------|:--------------|:----------------|:------------------|:----------------|:--------------|:----------------|:--------------|
| no_change         | 0.984 (0.012) | 0.916 (0.054) | nan (nan)       | nan (nan)         | 0.000 (0.000)   | 2.446 (1.615) | 40.442 (31.730) | 0.006 (0.004) |
| global_mean_shift | 0.989 (0.009) | 0.942 (0.042) | 0.568 (0.261)   | 0.691 (0.313)     | 0.243 (0.127)   | 1.803 (1.268) | 32.098 (25.965) | 0.004 (0.003) |
| additive          | 0.996 (0.004) | 0.990 (0.009) | 0.899 (0.074)   | 0.961 (0.039)     | 0.597 (0.156)   | 0.735 (0.677) | 20.066 (16.170) | 0.001 (0.001) |
| additive_ridge    | 0.997 (0.003) | 0.990 (0.009) | 0.908 (0.080)   | 0.967 (0.034)     | 0.597 (0.152)   | 0.575 (0.476) | 18.472 (14.251) | 0.001 (0.001) |

## 2. Model ranking under each metric (best -> worst, by mean)

| metric | ranking |
|---|---|
| pcc_all | additive_ridge > additive > global_mean_shift > no_change |
| pcc_top50 | additive > additive_ridge > global_mean_shift > no_change |
| pcc_all_delta | additive_ridge > additive > global_mean_shift > no_change |
| pcc_top50_delta | additive_ridge > additive > global_mean_shift > no_change |
| deg_overlap50 | additive > additive_ridge > global_mean_shift > no_change |
| edist | additive_ridge > additive > global_mean_shift > no_change |
| sinkhorn | additive_ridge > additive > global_mean_shift > no_change |

## 3. Per-condition win rate: additive beats global_mean_shift

| metric | additive wins | ties | n |
|---|---|---|---|
| pcc_all | 89.3% | 0.0% | 131 |
| pcc_top50 | 97.7% | 0.0% | 131 |
| pcc_all_delta | 99.2% | 0.0% | 131 |
| pcc_top50_delta | 98.5% | 0.0% | 131 |
| deg_overlap50 | 99.2% | 0.8% | 131 |
| edist | 88.5% | 0.0% | 131 |
| sinkhorn | 90.1% | 0.0% | 131 |

## 4. Do the metrics agree on which conditions are predicted well? (Spearman rho across conditions, additive model)

|                 |   pcc_all |   pcc_top50 |   pcc_all_delta |   pcc_top50_delta |   deg_overlap50 |   edist |   sinkhorn |
|:----------------|----------:|------------:|----------------:|------------------:|----------------:|--------:|-----------:|
| pcc_all         |      1    |        0.79 |            0.2  |              0.4  |            0.47 |    0.97 |       0.83 |
| pcc_top50       |      0.79 |        1    |            0.04 |              0.41 |            0.25 |    0.77 |       0.63 |
| pcc_all_delta   |      0.2  |        0.04 |            1    |              0.73 |            0.68 |    0.11 |      -0.01 |
| pcc_top50_delta |      0.4  |        0.41 |            0.73 |              1    |            0.49 |    0.35 |       0.18 |
| deg_overlap50   |      0.47 |        0.25 |            0.68 |              0.49 |            1    |    0.4  |       0.32 |
| edist           |      0.97 |        0.77 |            0.11 |              0.35 |            0.4  |    1    |       0.87 |
| sinkhorn        |      0.83 |        0.63 |           -0.01 |              0.18 |            0.32 |    0.87 |       1    |

## 5. Confounding by perturbation strength (E-distance of control vs. true cells)

Spearman rho(strength, metric) for the additive model: pcc_all = -0.44, pcc_all_delta = 0.49, edist = 0.47

## 5b. DEG overlap is capped: ceiling for *any* mean-shift model

Plugging the true mean into the predicted-DEG ranking recovers only 0.78 (median 0.82, min 0.10) of the Wilcoxon top-50, so deg_overlap50 values should be read against that cap. Fraction of the ceiling reached: additive 0.76, global_mean_shift 0.31.

## 6. The no-change reference

Predicting the control mean for every double gets pcc_all = 0.984 (min 0.926) and pcc_top50 = 0.916, but deg_overlap50 is undefined/0 and its deltas are all zero.

## 7. Largest rank disagreements for the additive model (pcc_all rank vs edist rank; 1 = best)

| condition      |   rank_pcc_all |   rank_edist |   rank_deg_overlap |   pcc_all |   pcc_all_delta |   deg_overlap50 |   edist |   edist_ctrl_vs_true |   n_cells |   gap |
|:---------------|---------------:|-------------:|-------------------:|----------:|----------------:|----------------:|--------:|---------------------:|----------:|------:|
| ZC3HAV1_CEBPA  |             90 |          122 |               44   |     0.995 |           0.963 |            0.68 |   1.821 |                5.395 |        86 |    32 |
| MAP2K6_IKZF3   |             22 |           50 |               15.5 |     0.998 |           0.963 |            0.76 |   0.409 |                3.381 |       300 |    28 |
| ETS2_IKZF3     |             54 |           81 |                1.5 |     0.997 |           0.953 |            0.84 |   0.66  |                4.023 |       445 |    27 |
| IRF1_SET       |             85 |           58 |               58   |     0.996 |           0.927 |            0.64 |   0.456 |                3.936 |       444 |    27 |
| CEBPB_OSR2     |             75 |          101 |               72.5 |     0.997 |           0.956 |            0.6  |   1.022 |                3.763 |       211 |    26 |
| FOSB_CEBPE     |             86 |          110 |               44   |     0.996 |           0.947 |            0.68 |   1.167 |                4.851 |       151 |    24 |
| KIF18B_KIF2C   |             50 |           27 |              130   |     0.998 |           0.635 |            0.04 |   0.259 |                0.333 |        96 |    23 |
| SAMD1_TGFBR2   |             60 |           38 |              125   |     0.997 |           0.867 |            0.3  |   0.34  |                1.544 |        72 |    22 |
| DUSP9_KLF1     |             44 |           23 |               51.5 |     0.998 |           0.948 |            0.66 |   0.248 |                4.001 |       366 |    21 |
| POU3F2_CBFA2T3 |             69 |           51 |               29.5 |     0.997 |           0.881 |            0.72 |   0.416 |                1.874 |       216 |    18 |

## 8. Conditions where metrics disagree on which model is better

113/131 conditions get the same winner under all 7 metrics; 18 flip depending on the metric.

| metric A | metric B | % conditions where A and B pick different winners |
|---|---|---|
| pcc_all | edist | 2.3% |
| pcc_all | deg_overlap50 | 9.9% |
| pcc_top50 | edist | 10.7% |
| deg_overlap50 | sinkhorn | 10.7% |
| edist | sinkhorn | 1.5% |
