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

Numbers are produced by `experiments/exp_gauge_plateau.py`; nothing is edited by hand.
