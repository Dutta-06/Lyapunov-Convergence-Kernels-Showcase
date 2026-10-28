# Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise (`exp_mass_dictionary.py`)

The loop at fixed positions: masses β (Σβ = m = 12, or free), amplitudes a = a_max tanh(α/a_max) (a_max = 1), w = β⊙a/m, G = 1, γ = 0.1; label noise on the measured error; the selection along the exact fits (`mass_dictionary`), in s = (m v/4)t. The paradigm's update: Euler on (α, log β), renormalised; the flat one: Euler on (α, √β), renormalised radially.

## A. Identities at random exact fits (5 × 12)

The closed-form drift against Katzenberger's drift of the scheme (`gauge_drift.limit_drift` on a chart of M), largest relative difference over three points:

| masses | the draws | the flat update | the paradigm's update (with its Itô term) |
|---|---|---|---|
| conserved | fresh, h → 0 | 1.2e-06 | 1.2e-06 |
| conserved | AR(1), τ_c = 200, h = 0.5 | 1.5e-06 | 1.5e-06 |
| free | fresh, h → 0 | 2.3e-06 | 2.3e-06 |
| free | AR(1), τ_c = 200, h = 0.5 | 2.0e-06 | 2.0e-06 |

The paradigm's drift against the flat drift without the masses' own selection (the γ-term of D_i differentiated in √β), largest relative difference: free, fresh 2.4e-15; free, AR(1) 1.0e-15; conserved, fresh 3.8e-03; conserved, AR(1) 8.1e-03.

K's spectrum against the flow's normal rates (−eigenvalues of its Jacobian in (α, √β)): largest relative difference conserved 3.7e-11, free 4.6e-11.

The excess-mass bound at the conserved law's points: | E/m | Σβa′/m | tr K | (G/(Nm))E²/Σ|φ_i|⁻² |
|---|---|---|---|
| 0.4964 | 0.6847 | 1.1702e-01 | 4.9288e-02 |
| 0.5413 | 0.7266 | 1.2444e-01 | 5.8599e-02 |
| 0.4350 | 0.6219 | 9.8324e-02 | 3.7844e-02 |

The gauge-minimal relevance D*(|w|) against (γa_max/m)|w| and the 4/3 law (the threshold γa_max³/(Gm) = 0.008333):

| |w| | D* | over (γa_max/m)|w| | over (3/2)(2G)^{1/3}(γ/m)^{2/3}|w|^{4/3} |
|---|---|---|---|
| 1e-05 | 1.6538e-08 | 0.1985 | 0.9882 |
| 0.0001 | 3.4139e-07 | 0.4097 | 0.9468 |
| 0.001 | 6.1370e-06 | 0.7364 | 0.7900 |
| 0.01 | 7.9402e-05 | 0.9528 | 0.4744 |
| 0.1 | 8.2904e-04 | 0.9948 | 0.2299 |
| 0.3 | 2.4957e-03 | 0.9983 | 0.1600 |

## B. The paradigm's noisy scheme against the reduced flow (5 × 12, the 2-sparse instance 4, step 1, σ = 0.1 for fresh draws and 0.03 for AR(1) draws, 64 paths)

| masses | the draws | steps | ‖w‖₁: the paths' mean (± s.e.) at s = 0, S/2, S | the flow's | the paths' mean w against the flow's (largest deviation, at the nine times) |
|---|---|---|---|---|---|
| conserved | fresh | 40000 | 0.7894 ± 0.0000, 0.8461 ± 0.0007, 0.8935 ± 0.0006 | 0.7894, 0.8439, 0.8954 | 0.000, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.002 |
| conserved | AR(1), τ_c = 1000 | 4444444 | 0.7894 ± 0.0000, 0.8415 ± 0.0054, 0.8734 ± 0.0054 | 0.7894, 0.8559, 0.9168 | 0.000, 0.007, 0.011, 0.015, 0.018, 0.021, 0.027, 0.033, 0.035 |
| free | fresh | 40000 | 0.8016 ± 0.0000, 0.9423 ± 0.0014, 0.9627 ± 0.0011 | 0.8016, 0.9422, 0.9625 | 0.000, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001 |
| free | AR(1), τ_c = 1000 | 4444444 | 0.8016 ± 0.0000, 0.8973 ± 0.0073, 0.9612 ± 0.0065 | 0.8016, 0.9191, 0.9645 | 0.000, 0.003, 0.008, 0.009, 0.009, 0.010, 0.010, 0.013, 0.008 |
| null kernel | fresh | 40000 | 0.7845 ± 0.0000, 0.8255 ± 0.0008, 0.8492 ± 0.0007 | 0.7845, 0.8233, 0.8490 | 0.000, 0.001, 0.001, 0.001, 0.001, 0.000, 0.001, 0.001, 0.001 |
| null kernel | AR(1), τ_c = 1000 | 4444444 | 0.7845 ± 0.0000, 0.8199 ± 0.0050, 0.8374 ± 0.0055 | 0.7845, 0.8308, 0.8608 | 0.000, 0.006, 0.010, 0.013, 0.016, 0.018, 0.022, 0.024, 0.024 |

## C. The selection's ends on Gaussian designs (5 × 12, 16 instances per k)

The paradigm's update; reduced flows to s = 1e+08; the start the noiseless flow's end from equal masses and amplitudes N(0, 0.05²). Means over the instances: ‖w‖₁/a_max (the truth's), the excess E/m (conserved law), the atoms used (|w| > 10⁻³), alive (β > 10⁻³); the fraction recovering the truth (within 10⁻⁴·max|w*|); the median distance to basis pursuit's solution; the null kernel's mass.

| k | masses | the noise | ‖w‖₁/a_max (truth's) | E/m | used | alive | the truth | median distance to basis pursuit | the null kernel's mass |
|---|---|---|---|---|---|---|---|---|---|
| 2 | — | basis pursuit | | | | | 0.81 | | |
| 2 | conserved | the start | 0.668 (0.485) | 0.332 | 12.0 | 12.0 | 0.00 | 0.172 | — |
| 2 | conserved | white | 0.988 (0.485) | 0.012 | 5.9 | 10.1 | 0.00 | 0.289 | — |
| 2 | conserved | slow, τ_c = 1e+03 | 0.990 (0.485) | 0.010 | 5.9 | 7.9 | 0.00 | 0.289 | — |
| 2 | conserved | slow, τ_c = 1e+04 | 0.994 (0.485) | 0.006 | 6.0 | 7.2 | 0.00 | 0.290 | — |
| 2 | conserved | white, the flat update | 0.975 (0.485) | 0.025 | 12.0 | 12.0 | 0.00 | 0.198 | — |
| 2 | free | the start | 0.676 (0.485) | nan | 11.9 | 12.0 | 0.00 | 0.172 | — |
| 2 | free | white | 0.468 (0.485) | nan | 2.6 | 2.7 | 0.81 | 0.000 | — |
| 2 | free | slow, τ_c = 1e+03 | 0.473 (0.485) | nan | 2.5 | 2.9 | 0.81 | 0.000 | — |
| 2 | free | slow, τ_c = 1e+04 | 0.471 (0.485) | nan | 2.4 | 2.6 | 0.88 | 0.000 | — |
| 2 | free | white, the flat update | 0.463 (0.485) | nan | 3.5 | 6.3 | 0.12 | 0.001 | — |
| 2 | null kernel | the start | 0.661 (0.485) | 0.336 | 12.0 | 12.0 | 0.00 | 0.171 | 0.77 |
| 2 | null kernel | white | 0.486 (0.485) | 0.485 | 2.8 | 2.9 | 0.75 | 0.000 | 6.11 |
| 2 | null kernel | slow, τ_c = 1e+03 | 0.610 (0.485) | 0.368 | 3.7 | 3.8 | 0.44 | 0.174 | 4.61 |
| 2 | null kernel | slow, τ_c = 1e+04 | 0.627 (0.485) | 0.353 | 3.6 | 3.9 | 0.44 | 0.226 | 4.41 |
| 2 | free | white: the flat update's end against the minimiser of Σ|φ_i|²D*(|w_i|) (median distance) | 1.1e-04 | | | | | that minimiser against basis pursuit's: 1.1e-03 | |
| 3 | — | basis pursuit | | | | | 0.38 | | |
| 3 | conserved | the start | 0.731 (0.671) | 0.269 | 11.9 | 12.0 | 0.00 | 0.140 | — |
| 3 | conserved | white | 0.988 (0.671) | 0.012 | 5.9 | 10.9 | 0.00 | 0.310 | — |
| 3 | conserved | slow, τ_c = 1e+03 | 0.992 (0.671) | 0.008 | 5.9 | 8.2 | 0.00 | 0.312 | — |
| 3 | conserved | slow, τ_c = 1e+04 | 0.995 (0.671) | 0.005 | 5.9 | 7.2 | 0.00 | 0.322 | — |
| 3 | conserved | white, the flat update | 0.974 (0.671) | 0.026 | 11.9 | 12.0 | 0.00 | 0.155 | — |
| 3 | free | the start | 0.763 (0.671) | nan | 11.9 | 12.0 | 0.00 | 0.138 | — |
| 3 | free | white | 0.581 (0.671) | nan | 4.1 | 5.0 | 0.44 | 0.000 | — |
| 3 | free | slow, τ_c = 1e+03 | 0.582 (0.671) | nan | 4.0 | 5.0 | 0.44 | 0.000 | — |
| 3 | free | slow, τ_c = 1e+04 | 0.591 (0.671) | nan | 3.9 | 4.9 | 0.50 | 0.003 | — |
| 3 | free | white, the flat update | 0.578 (0.671) | nan | 5.1 | 7.8 | 0.00 | 0.002 | — |
| 3 | null kernel | the start | 0.726 (0.671) | 0.273 | 11.9 | 12.0 | 0.00 | 0.141 | 0.66 |
| 3 | null kernel | white | 0.587 (0.671) | 0.402 | 4.2 | 5.4 | 0.38 | 0.004 | 4.87 |
| 3 | null kernel | slow, τ_c = 1e+03 | 0.655 (0.671) | 0.335 | 3.8 | 4.7 | 0.62 | 0.106 | 4.08 |
| 3 | null kernel | slow, τ_c = 1e+04 | 0.743 (0.671) | 0.251 | 4.1 | 5.0 | 0.44 | 0.179 | 3.02 |
| 3 | free | white: the flat update's end against the minimiser of Σ|φ_i|²D*(|w_i|) (median distance) | 1.3e-04 | | | | | that minimiser against basis pursuit's: 2.0e-03 | |

## D. The selection's ends on Gaussian kernel dictionaries (16 centres, 8 samples on [0, 1], width 0.08, 8 instances per k)

The paradigm's update; reduced flows to s = 1e+08; the start the noiseless flow's end from equal masses and amplitudes N(0, 0.05²). Means over the instances: ‖w‖₁/a_max (the truth's), the excess E/m (conserved law), the atoms used (|w| > 10⁻³), alive (β > 10⁻³); the fraction recovering the truth (within 10⁻⁴·max|w*|); the median distance to basis pursuit's solution; the null kernel's mass.

| k | masses | the noise | ‖w‖₁/a_max (truth's) | E/m | used | alive | the truth | median distance to basis pursuit | the null kernel's mass |
|---|---|---|---|---|---|---|---|---|---|
| 2 | — | basis pursuit | | | | | 0.88 | | |
| 2 | conserved | the start | 0.590 (0.453) | 0.410 | 13.9 | 16.0 | 0.00 | 0.118 | — |
| 2 | conserved | white | 0.981 (0.453) | 0.019 | 7.8 | 15.5 | 0.00 | 0.327 | — |
| 2 | conserved | slow, τ_c = 1e+03 | 0.983 (0.453) | 0.017 | 7.6 | 15.2 | 0.00 | 0.334 | — |
| 2 | conserved | slow, τ_c = 1e+04 | 0.986 (0.453) | 0.014 | 7.5 | 14.1 | 0.00 | 0.356 | — |
| 2 | conserved | white, the flat update | 0.966 (0.453) | 0.034 | 16.0 | 16.0 | 0.00 | 0.138 | — |
| 2 | free | the start | 0.611 (0.453) | nan | 14.1 | 16.0 | 0.00 | 0.124 | — |
| 2 | free | white | 0.428 (0.453) | nan | 2.2 | 2.8 | 0.88 | 0.000 | — |
| 2 | free | slow, τ_c = 1e+03 | 0.428 (0.453) | nan | 2.2 | 2.9 | 0.88 | 0.000 | — |
| 2 | free | slow, τ_c = 1e+04 | 0.453 (0.453) | nan | 2.0 | 2.0 | 1.00 | 0.000 | — |
| 2 | free | white, the flat update | 0.422 (0.453) | nan | 2.7 | 6.1 | 0.14 | 0.000 | — |
| 2 | null kernel | the start | 0.585 (0.453) | 0.412 | 13.9 | 16.0 | 0.00 | 0.117 | 0.50 |
| 2 | null kernel | white | 0.428 (0.453) | 0.537 | 2.2 | 2.8 | 0.88 | 0.000 | 6.82 |
| 2 | null kernel | slow, τ_c = 1e+03 | 0.453 (0.453) | 0.513 | 2.0 | 2.0 | 1.00 | 0.000 | 6.53 |
| 2 | null kernel | slow, τ_c = 1e+04 | 0.453 (0.453) | 0.513 | 2.0 | 2.0 | 1.00 | 0.000 | 6.53 |
| 2 | free | white: the flat update's end against the minimiser of Σ|φ_i|²D*(|w_i|) (median distance) | 5.1e-05 | | | | | that minimiser against basis pursuit's: 4.0e-04 | |
| 3 | — | basis pursuit | | | | | 0.38 | | |
| 3 | conserved | the start | 0.749 (0.654) | 0.251 | 15.8 | 16.0 | 0.00 | 0.119 | — |
| 3 | conserved | white | 0.983 (0.654) | 0.017 | 7.8 | 14.0 | 0.00 | 0.283 | — |
| 3 | conserved | slow, τ_c = 1e+03 | 0.986 (0.654) | 0.014 | 7.5 | 13.1 | 0.00 | 0.299 | — |
| 3 | conserved | slow, τ_c = 1e+04 | 0.990 (0.654) | 0.010 | 7.5 | 9.9 | 0.00 | 0.292 | — |
| 3 | conserved | white, the flat update | 0.966 (0.654) | 0.034 | 16.0 | 16.0 | 0.00 | 0.132 | — |
| 3 | free | the start | 0.800 (0.654) | nan | 15.8 | 16.0 | 0.00 | 0.129 | — |
| 3 | free | white | 0.563 (0.654) | nan | 4.4 | 7.0 | 0.50 | 0.001 | — |
| 3 | free | slow, τ_c = 1e+03 | 0.566 (0.654) | nan | 4.0 | 6.8 | 0.50 | 0.000 | — |
| 3 | free | slow, τ_c = 1e+04 | 0.591 (0.654) | nan | 3.8 | 6.0 | 0.62 | 0.012 | — |
| 3 | free | white, the flat update | 0.560 (0.654) | nan | 5.6 | 9.1 | 0.12 | 0.002 | — |
| 3 | null kernel | the start | 0.740 (0.654) | 0.258 | 16.0 | 16.0 | 0.00 | 0.118 | 0.43 |
| 3 | null kernel | white | 0.571 (0.654) | 0.409 | 4.5 | 6.5 | 0.50 | 0.001 | 5.05 |
| 3 | null kernel | slow, τ_c = 1e+03 | 0.599 (0.654) | 0.382 | 4.1 | 5.8 | 0.62 | 0.060 | 4.73 |
| 3 | null kernel | slow, τ_c = 1e+04 | 0.599 (0.654) | 0.382 | 4.1 | 6.4 | 0.62 | 0.060 | 4.73 |
| 3 | free | white: the flat update's end against the minimiser of Σ|φ_i|²D*(|w_i|) (median distance) | 5.3e-04 | | | | | that minimiser against basis pursuit's: 1.3e-03 | |

Flows over the 900-s limit, left out of the means (the instances' other flows kept): k = 2, free, flat: 1.

## E. The conserved law's white end against the representability boundary (5 × 12, 2-sparse truths, 12 instances; those whose basis-pursuit solution is not representable, ‖w_bp‖₁ ≥ a_max, left out)

| a_max | instances | ‖w‖₁/a_max at the start | at the end | E/m at the end | the atoms used, alive | ‖w_bp‖₁/a_max |
|---|---|---|---|---|---|---|
| 0.5 | 7 | 0.987 | 0.998 (min 0.997) | 0.0021 | 6.0, 6.9 | 0.870 |
| 1 | 12 | 0.697 | 0.988 (min 0.984) | 0.0119 | 5.9, 9.9 | 0.476 |
| 2 | 12 | 0.327 | 0.942 (min 0.930) | 0.0585 | 6.9, 12.0 | 0.238 |
| 4 | 12 | 0.161 | 0.653 (min 0.582) | 0.3469 | 11.5, 12.0 | 0.119 |

Numbers are produced by `experiments/exp_mass_dictionary.py`; nothing is edited by hand.
