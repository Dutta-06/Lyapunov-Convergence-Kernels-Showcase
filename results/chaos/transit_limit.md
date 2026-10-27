# Theorem D‴.53: the sliding limit (`exp_transit_limit.py`)

The free law with moving centres (m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025) from α = 0 with the masses ε, Gaussian kernels of width 0.12 on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, against its sliding path (polished); the targets: of D‴.51's 48 out-of-sample three-bump targets (seed 51) whose path is exact with simple entries and no exits, the 5 with the largest margins (≥ 0.05) for the kernels that never enter.

## A. The rates

Per target and ε: the kernels formed; the block at the midpoints between entries against S_k (max over the coefficients and centres); the never-entering kernels' largest coefficient (× ε^{−3/2}); their centres' largest displacement against the creep's εC_j (the largest ratio's distance from 1); the duals: obs = s_j(u_j − u_j^path)/ε^{1/2} against the first-order deviation — the entries' pulls Δ_kj plus each entry's shift σ_k = −dev_{j_k}/v_k times the drop of j's rate at it, by recursion over the entries — plus the creep's ε^{1/2}D_j (the largest |obs − pred|, the median ratio obs/pred where |pred| > 0.1, and the largest |obs − ΣΔ| with the pulls alone); the entries' crossing times (τ_θ − τ_k)/ε^{1/2} against s_θ + σ_k (the largest difference).

| target (margin) | ε | formed | block − S_k | max \|w_never\| (× ε^{−3/2}) | displacement / εC − 1 | duals: max \|obs − pred\| (median obs/pred) | the pulls alone | crossings: max \|obs − pred\| |
|---|---|---|---|---|---|---|---|---|
| 39 (0.151) | 1e-06 | 3 | 9.41e-8 | 3.15e-8 (31.5) | 0.0807 | 0.0048 (0.9990) | 0.3069 | 0.1507 |
| 39 (0.151) | 1e-07 | 3 | 1.25e-9 | 9.21e-10 (29.1) | 0.0246 | 0.0022 (0.9996) | 0.1435 | 0.0464 |
| 39 (0.151) | 1e-08 | 3 | 3.84e-11 | 2.84e-11 (28.4) | 0.00771 | 0.0009 (0.9999) | 0.0937 | 0.0141 |
| 39 (0.151) | 1e-09 | 3 | 1.21e-12 | 8.93e-13 (28.2) | 0.00243 | 0.0003 (1.0000) | 0.0781 | 0.0041 |
| 39 (0.151) | 1e-10 | 3 | 3.73e-14 | 2.82e-14 (28.2) | 0.000767 | 0.0001 (1.0000) | 0.0732 | 0.0013 |
| 39 (0.151) | 1e-11 | 3 | 2.52e-14 | 8.90e-16 (28.2) | 0.000244 | 0.0002 (1.0000) | 0.0715 | 0.0029 |
| 24 (0.140) | 1e-06 | 3 | 1.35e-6 | 3.85e-8 (38.5) | 0.0877 | 0.1894 (0.9994) | 0.7272 | 2.6912 |
| 24 (0.140) | 1e-07 | 3 | 4.46e-9 | 1.12e-9 (35.4) | 0.0264 | 0.0608 (0.9997) | 0.3559 | 0.9044 |
| 24 (0.140) | 1e-08 | 3 | 1.46e-10 | 3.45e-11 (34.5) | 0.00824 | 0.0188 (0.9999) | 0.2371 | 0.2829 |
| 24 (0.140) | 1e-09 | 3 | 4.69e-12 | 1.08e-12 (34.2) | 0.00259 | 0.0053 (1.0000) | 0.1994 | 0.0799 |
| 24 (0.140) | 1e-10 | 3 | 1.53e-13 | 3.41e-14 (34.1) | 0.000819 | 0.0011 (1.0000) | 0.1874 | 0.0149 |
| 24 (0.140) | 1e-11 | 3 | 4.88e-15 | 1.08e-15 (34.1) | 0.000294 | 0.0005 (1.0000) | 0.1836 | 0.0077 |
| 19 (0.090) | 1e-06 | 3 | 1.50e-7 | 1.58e-7 (158) | 0.131 | 0.0161 (0.9993) | 0.2739 | 0.1114 |
| 19 (0.090) | 1e-07 | 3 | 4.15e-9 | 4.37e-9 (138) | 0.0399 | 0.0016 (0.9997) | 0.0814 | 0.0351 |
| 19 (0.090) | 1e-08 | 3 | 1.26e-10 | 1.33e-10 (133) | 0.0125 | 0.0007 (0.9999) | 0.0247 | 0.0110 |
| 19 (0.090) | 1e-09 | 3 | 3.95e-12 | 4.14e-12 (131) | 0.00392 | 0.0003 (1.0000) | 0.0085 | 0.0047 |
| 19 (0.090) | 1e-10 | 3 | 1.30e-13 | 1.30e-13 (130) | 0.00124 | 0.0002 (1.0000) | 0.0034 | 0.0035 |
| 19 (0.090) | 1e-11 | 3 | 1.03e-14 | 4.12e-15 (130) | 0.00039 | 0.0003 (1.0000) | 0.0017 | 0.0040 |
| 6 (0.084) | 1e-06 | 3 | 0.186 | 2.03e-7 (203) | 0.18 | 0.4400 (1.0001) | 0.4101 | 0.1224 |
| 6 (0.084) | 1e-07 | 3 | 5.06e-6 | 5.54e-9 (175) | 0.0523 | 0.0036 (1.0000) | 0.1288 | 0.0305 |
| 6 (0.084) | 1e-08 | 3 | 2.06e-9 | 1.68e-10 (168) | 0.0161 | 0.0004 (0.9999) | 0.0743 | 0.0071 |
| 6 (0.084) | 1e-09 | 3 | 6.62e-11 | 5.24e-12 (166) | 0.00506 | 0.0003 (0.9999) | 0.0573 | 0.0047 |
| 6 (0.084) | 1e-10 | 3 | 2.11e-12 | 1.65e-13 (165) | 0.0016 | 0.0004 (0.9999) | 0.0520 | 0.0071 |
| 6 (0.084) | 1e-11 | 3 | 6.79e-14 | 5.21e-15 (165) | 0.000505 | 0.0005 (0.9999) | 0.0502 | 0.0089 |
| 41 (0.079) | 1e-06 | 3 | 2.08e-7 | 2.21e-7 (221) | 0.159 | 0.0481 (0.9991) | 0.4899 | 0.4724 |
| 41 (0.079) | 1e-07 | 3 | 5.65e-9 | 6.06e-9 (192) | 0.047 | 0.0058 (0.9996) | 0.2547 | 0.1481 |
| 41 (0.079) | 1e-08 | 3 | 1.71e-10 | 1.84e-10 (184) | 0.0145 | 0.0016 (0.9999) | 0.2018 | 0.0464 |
| 41 (0.079) | 1e-09 | 3 | 5.34e-12 | 5.76e-12 (182) | 0.00457 | 0.0005 (1.0000) | 0.1853 | 0.0145 |
| 41 (0.079) | 1e-10 | 3 | 1.65e-13 | 1.81e-13 (181) | 0.00144 | 0.0006 (1.0000) | 0.1802 | 0.0126 |
| 41 (0.079) | 1e-11 | 3 | 4.33e-15 | 5.73e-15 (181) | 0.000455 | 0.0008 (1.0000) | 0.1785 | 0.0166 |

Slopes in log ε (least squares over ε): target 39: block 1.370, w_never 1.508, displacement 1.006, duals' residual 0.322; target 24: block 1.632, w_never 1.509, displacement 1.006, duals' residual 0.538; target 19: block 1.452, w_never 1.514, displacement 1.009, duals' residual 0.344; target 6: block 2.366, w_never 1.515, displacement 1.012, duals' residual 0.499; target 41: block 1.529, w_never 1.515, displacement 1.011, duals' residual 0.356 (predicted 1.5, 1.5, 1 and, for the residual in units of ε^{1/2}, ≥ 0.25).

The entries (τ_k, kernel), s_θ and the pulls on the never-entering kernels (largest |Δ|): target 39: (11.1, 15) s_θ 74.53 max|Δ| 5.684, (11.75, 10) s_θ 77.25 max|Δ| 5.586, (13.58, 2) s_θ 87.32 max|Δ| 5.446; target 24: (9.676, 9) s_θ 66.82 max|Δ| 5.880, (14.31, 18) s_θ 85.55 max|Δ| 5.314, (15.14, 14) s_θ 76.38 max|Δ| 4.988; target 19: (11.31, 9) s_θ 75.62 max|Δ| 5.677, (14.38, 2) s_θ 86.91 max|Δ| 5.476, (16.43, 16) s_θ 97.53 max|Δ| 5.336; target 6: (10.79, 15) s_θ 67.78 max|Δ| 5.786, (15.48, 3) s_θ 90.93 max|Δ| 5.378, (15.7, 10) s_θ 85.37 max|Δ| 5.205; target 41: (10.62, 16) s_θ 71.79 max|Δ| 5.797, (20.05, 5) s_θ 103.4 max|Δ| 4.815, (23.65, 9) s_θ 130.1 max|Δ| 4.766.

Numbers are produced by `experiments/exp_transit_limit.py`; nothing is edited by hand.
