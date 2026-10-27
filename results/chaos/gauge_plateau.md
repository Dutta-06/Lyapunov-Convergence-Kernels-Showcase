# Proposition D‴.39: dissipativity of the dithered gauge (`exp_gauge_plateau.py`)

Proposition D‴.36's instance (three bumps, one particle per class, m = 3, the cohesion leak); energies over τ = 0.0001; the box L = (5.0, 3.0, 4.0) in (α, γ, κ), σ_B = 1: |a| ≤ 3.393, the representable masses W_k ≥ 0.2947, 0.2063, 0.2358 (the saturation's 0.2500, 0.1750, 0.2000).

## A. The plateau: the lowest configuration with one amplitude fixed far out

| masses | class | α_k | ΔJ/τ above the exact fit | fit error | a | c | s |
|---|---|---|---|---|---|---|---|
| uniform | 1 | ±20 | 7.718 | 0.0350 | 4.00, -3.73, 3.15 | 1.09, 0.34, -1.41 | 0.72, 0.74, 0.58 |
| uniform | 2 | ±20 | 7.762 | 0.0351 | 3.11, -4.00, 4.00 | -1.42, 0.32, 1.04 | 0.58, 0.70, 0.75 |
| uniform | 3 | ±20 | 7.718 | 0.0350 | 3.15, -3.73, 4.00 | -1.41, 0.34, 1.09 | 0.58, 0.74, 0.72 |
| sharp | 1 | ±20 | 7.637 | 0.0348 | 4.00, -3.49, 4.00 | -1.32, -0.02, 1.28 | 0.66, 1.17, 0.67 |
| sharp | 2 | ±20 | 7.707 | 0.0350 | 2.78, -4.00, 3.89 | -1.41, 0.31, 1.10 | 0.58, 0.72, 0.73 |
| sharp | 3 | ±20 | 7.637 | 0.0348 | 4.00, -3.49, 4.00 | -1.32, -0.02, 1.28 | 0.66, 1.17, 0.67 |
| uniform | 1 | ±6 | 7.692 | 0.0349 | 3.62, -3.38, 3.14 | 1.13, 0.29, -1.41 | 0.72, 0.72, 0.58 |
| uniform | 2 | ±6 | 7.710 | 0.0350 | 3.14, -3.62, 3.83 | -1.41, 0.32, 1.10 | 0.58, 0.73, 0.73 |
| uniform | 3 | ±6 | 7.692 | 0.0349 | 3.14, -3.38, 3.62 | -1.41, 0.29, 1.13 | 0.58, 0.72, 0.72 |
| sharp | 1 | ±6 | 7.722 | 0.0350 | 3.62, 3.52, -3.92 | 1.08, -1.41, 0.35 | 0.73, 0.58, 0.74 |
| sharp | 2 | ±6 | 7.670 | 0.0349 | 2.77, -3.62, 3.51 | -1.42, 0.26, 1.15 | 0.58, 0.70, 0.72 |
| sharp | 3 | ±6 | 7.681 | 0.0349 | 2.78, -3.67, 3.62 | -1.41, 0.28, 1.14 | 0.58, 0.72, 0.72 |
| uniform | 1 | ±L_α (box) | 7.662 | 0.0349 | 3.39, -3.17, 3.14 | 1.16, 0.26, -1.41 | 0.71, 0.71, 0.58 |
| uniform | 2 | ±L_α (box) | 7.735 | 0.0350 | 3.09, -3.39, 3.39 | -1.43, 0.26, 1.11 | 0.57, 0.68, 0.74 |
| uniform | 3 | ±L_α (box) | 7.662 | 0.0349 | 3.14, -3.17, 3.39 | -1.41, 0.26, 1.16 | 0.58, 0.71, 0.71 |
| sharp | 1 | ±L_α (box) | 8.048 | 0.0357 | 3.39, -3.18, 3.39 | -1.32, 0.01, 1.27 | 0.64, 0.92, 0.66 |
| sharp | 2 | ±L_α (box) | 7.626 | 0.0348 | 2.75, -3.39, 3.27 | -1.42, 0.23, 1.18 | 0.58, 0.69, 0.71 |
| sharp | 3 | ±L_α (box) | 7.650 | 0.0348 | 2.79, -3.39, 3.39 | -1.41, 0.24, 1.18 | 0.58, 0.72, 0.70 |

## B. The minimum-energy path from the exact fit to the plateau (string method)

Masses uniform, toward class 3's saturation: the largest ΔJ/τ along the path 7.715 (image 39 of 39); at every fourth image the profile 0.00, 1.86, 5.33, 7.17, 7.57, 7.66, 7.70, 7.71, 7.71, 7.71, the fit error 0.0000, 0.0172, 0.0291, 0.0337, 0.0346, 0.0349, 0.0349, 0.0350, 0.0350, 0.0350, the second width 0.40, 0.46, 0.54, 0.63, 0.68, 0.71, 0.74, 0.75, 0.75, 0.74.

Masses sharp, toward class 3's saturation: the largest ΔJ/τ along the path 8.211 (image 6 of 39); at every fourth image the profile 0.00, 8.04, 8.19, 8.17, 8.12, 7.90, 7.79, 7.74, 7.70, 7.68, the fit error 0.0000, 0.0357, 0.0360, 0.0360, 0.0359, 0.0354, 0.0352, 0.0350, 0.0349, 0.0349, the second width 0.40, 0.77, 0.91, 0.93, 0.96, 1.07, 1.12, 1.14, 1.15, 1.16.

## C. The escapes at fixed masses 1/3 (γ = 0; 24 seeds, at most 400000 slow units; an escape: max|α_k| > 6)

| τ | escapes | the rate (per slow unit, ± s.e.) | the median time | the rate × e^{ΔJ/τ} | the fit error at the escape |
|---|---|---|---|---|---|
| 0.00015 | 10/24 | 1.33e-6 ± 4.22e-7 | 4e+05 | 2.25e-4 | 0.0456 |
| 0.0002 | 19/24 | 4.06e-6 ± 9.32e-7 | 1.6e+05 | 1.90e-4 | 0.0465 |
| 0.0003 | 24/24 | 1.09e-5 ± 2.23e-6 | 6.65e+04 | 1.42e-4 | 0.0538 |

The slope of log(rate) against 1/τ: 6.31e-4 (the lowest J on the escape set's boundary, ΔJ = 7.69e-4); weighted by the escape counts (log(rate)'s variance ≈ 1/n): 6.25e-4 ± 1.09e-4, 0.81 ± 0.14 of ΔJ (referee round 31).

## D. The dead zone at fixed masses 1/3 (γ = 0, τ = 0.0002, 16 seeds, 400000 slow units)

A record is on the plateau when its second width exceeds 0.6 (A and B: the plateaus' configurations and the path's).

| box | seeds that reach the plateau | fraction of the records on it (s_2 > 0.6), by quarters | fraction with max|α| > 6 | the largest |α|, by quarters | the fit error, by quarters |
|---|---|---|---|---|---|
| none | 16/16 | 0.179, 0.401, 0.535, 0.568 | 0.090, 0.358, 0.461, 0.559 | 15.0, 19.8, 21.0, 27.6 | 0.0395, 0.0425, 0.0439, 0.0448 |
| L = (5.0, 3.0, 4.0) | 16/16 | 0.123, 0.086, 0.185, 0.157 | 0.000, 0.000, 0.000, 0.000 | 5.0, 5.0, 5.0, 5.0 | 0.0391, 0.0392, 0.0403, 0.0403 |

## E. The full loop with the box (γ = 0.01, τ = 0.0001, 12 seeds, 1e+06 slow units)

The box's vertices: 0.5579, 0.2063, 0.2358; 0.2947, 0.4695, 0.2358; 0.2947, 0.2063, 0.4990; W♯ = 0.3765, 0.2988, 0.3247.

| κ_F | s | mean of W_1, W_2, W_3 (± s.e.) | std |
|---|---|---|---|
| 0 | 100000 | 0.3578 ± 0.0111, 0.3272 ± 0.0105, 0.3150 ± 0.0125 | 0.0385, 0.0364, 0.0434 |
| 0 | 300000 | 0.3672 ± 0.0158, 0.3172 ± 0.0147, 0.3156 ± 0.0113 | 0.0547, 0.0509, 0.0392 |
| 0 | 1e+06 | 0.3658 ± 0.0141, 0.3171 ± 0.0216, 0.3171 ± 0.0164 | 0.0490, 0.0747, 0.0569 |
| 2 | 100000 | 0.3692 ± 0.0077, 0.3070 ± 0.0074, 0.3237 ± 0.0077 | 0.0267, 0.0255, 0.0265 |
| 2 | 300000 | 0.3695 ± 0.0103, 0.3017 ± 0.0099, 0.3288 ± 0.0068 | 0.0358, 0.0342, 0.0236 |
| 2 | 1e+06 | 0.3803 ± 0.0093, 0.2961 ± 0.0122, 0.3236 ± 0.0102 | 0.0322, 0.0424, 0.0354 |

| κ_F | the masses at the end, by seed (the nearest vertex for κ_F = 0, its distance) | the largest |α| | the fit error (second half) |
|---|---|---|---|
| 0 | 0.342, 0.389, 0.269 (v2, 0.099); 0.371, 0.221, 0.408 (v3, 0.120); 0.296, 0.347, 0.356 (v2, 0.172); 0.470, 0.224, 0.306 (v1, 0.114); 0.354, 0.292, 0.355 (v3, 0.178); 0.316, 0.348, 0.336 (v2, 0.159); 0.316, 0.375, 0.309 (v2, 0.122); 0.388, 0.227, 0.384 (v3, 0.149); 0.366, 0.399, 0.235 (v2, 0.100); 0.364, 0.395, 0.242 (v2, 0.102); 0.431, 0.222, 0.347 (v1, 0.170); 0.377, 0.367, 0.257 (v2, 0.133) | 5.03 | 0.0257 |
| 2 | 0.373, 0.308, 0.319 (0.010); 0.395, 0.260, 0.345 (0.039); 0.300, 0.375, 0.325 (0.076); 0.404, 0.263, 0.333 (0.036); 0.391, 0.311, 0.299 (0.026); 0.350, 0.290, 0.360 (0.036); 0.361, 0.275, 0.364 (0.039); 0.410, 0.282, 0.308 (0.034); 0.417, 0.290, 0.294 (0.040); 0.376, 0.378, 0.246 (0.080); 0.386, 0.242, 0.373 (0.057); 0.402, 0.279, 0.318 (0.026) | 5.01 | 0.026 |

κ_F = 0: the records coincide with exp_gauge_drift.py's Part H (the same seeds without the box) up to the first record at which a coordinate has left B: 179000, 101000, 23000, 411000, 348000, 17000, 92000, 76000, 130000, 65000, 236000, 270000.

κ_F = 2: the records coincide with exp_gauge_fixman.py's Part E (the same seeds without the box) up to the first record at which a coordinate has left B: never, 372000, 24000, 691000, 484000, 70000, 316000, 257000, 207000, 354000, 658000, 279000.

## F. The two channels to the saturation in the stored runs

| runs | class | seeds with records at |α_k| > 10 | such records | of them with W_k − |A_k|/a_max < 0.01 | > 0.05 | their mean fit error |
|---|---|---|---|---|---|---|
| no leak (exp_gauge_drift.py H) | 1 | 8/12 | 2350 | 2349 | 0 | 0.0245 |
| no leak (exp_gauge_drift.py H) | 2 | 3/12 | 138 | 138 | 0 | 0.0244 |
| no leak (exp_gauge_drift.py H) | 3 | 7/12 | 1373 | 1373 | 0 | 0.0244 |
| the leak, no correction (exp_gauge_fixman.py F) | 1 | 2/16 | 516 | 0 | 516 | 0.0427 |
| the leak, no correction (exp_gauge_fixman.py F) | 2 | 0/16 | 0 | 0 | 0 | — |
| the leak, no correction (exp_gauge_fixman.py F) | 3 | 1/16 | 13 | 0 | 13 | 0.0408 |

Numbers are produced by `experiments/exp_gauge_plateau.py`; nothing is edited by hand.
