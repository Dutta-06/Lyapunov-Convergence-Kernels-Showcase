# Proposition D‴.57: the transit's dead zone (`exp_dead_zone.py`)

D‴.56's setting — 21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]; three features at least 0.45 apart with m > 0.8q; y = f* + ξ, ξ i.i.d. N(0, s²) —, the gate open while the residual's rms exceeds (1 + δ)s, δ = 0.1.

## A. The sliding limit

Per level s: 10 draws on each of the 40 targets; on the draws whose noisy path enters exactly E first: the dead zone's stop (the first entry whose relaxed rms is ≤ (1 + δ)s) at the K-th entry, later or earlier; at the K-th, the energy dead zone's frozen fit against the oracle's (the K-th block relaxed): the errors' ratio, the mean square over its prediction given the stop at the K-th entry (± its standard error) and over s²(4K/N + 2δ + δ²) (the unconditional value), the oracle's over s²·2K/N, the unrelaxed tangent part |u|² over (1 + δ)²s²N − |r_K|², the mean |u|²/(s²N) over its prediction given the stop, the correlation of the oracle's error with u and |cos| of u against the target's first draw's:

| s | exact draws | stop at K, later, earlier | err_dz/err_or: median [5–95%] | mean err_dz²: / prediction given the stop; / s²(4K/N + 2δ + δ²) | mean err_or² / (s²2K/N) | \|u\|² / (1 + δ)²s²N − \|r_K\|²: median [5–95%] | mean \|u\|²/(s²N) / prediction given the stop | corr(err_or, u): mean ± s.e. | \|cos\| of u: median |
|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 400 | 386, 14, 0 | 2.154 [1.135, 3.816] | 1.023 ± 0.026; 1.071 | 0.991 | 1.0000 [1.0000, 1.0000] | 0.985 | -0.025 ± 0.021 | 1.0000 |
| 0.03 | 399 | 381, 18, 0 | 2.120 [1.121, 3.998] | 0.980 ± 0.025; 1.026 | 0.963 | 1.0000 [1.0000, 1.0000] | 0.986 | +0.007 ± 0.020 | 0.9998 |
| 0.1 | 369 | 356, 13, 0 | 2.079 [1.070, 3.980] | 1.006 ± 0.027; 1.054 | 1.021 | 1.0000 [0.9999, 1.0002] | 0.986 | -0.024 ± 0.021 | 0.9973 |

At first order |u|²/N = (1 + δ)²s² − |r_K|²/N, of mean s²(2δ + δ² + 2K/N): |u|/(s√N) ≈ 0.555 against the oracle's s(2K/N)^{1/2} = 0.314s; without the cross term the errors' rms ratio would be ((4K/N + 2δ + δ²)/(2K/N))^{1/2} = 2.033. The stop at the K-th entry is the event |r_K|² ≤ (1 + δ)²s²N, |r_K|²/s² ~ χ²_{N−2K} (probability 0.954); given it E|r_K|² = s²(N − 2K)P(χ²_{N−2K+2} ≤ c)/P(χ²_{N−2K} ≤ c), c = (1 + δ)²N, = 53.83s² against 55s², so E|u|²/N = 0.3276s² and E err_dz² = s²(2K + c − E[χ²_{N−2K} | ≤ c])/N = 0.4260s² (×1.047 the unconditional value; the rms ratio 2.081).

## B. The loop

The loop at ε = 1e-08 from α = 0 on 16 targets × 3 draws at s = 0.01, 0.03, 0.1 (the draws whose noisy path is exactly E), to τ_K + 3W: on the draws whose sliding limit stops at the K-th entry, the kernels formed (|w| > 0.001) and the fits' rms errors over the oracle's (the limit's K-th block relaxed), the energy dead zone's against the ungated loop's and the gated loops against the sliding limit's fits (the energy dead zone's frozen fit, the oracle's); on the late draws (the stop after the K-th entry) the kernels formed:

| s | draws: stop at K, late | formed at K: ungated, energy, formation (median [range]) | err / err_or (medians): ungated, energy, formation | energy against the ungated (median) | the loop against the limit: max \|err − err_dz\| (energy), max \|err − err_or\| (formation) | the gate's τ / τ_K (median) | formed, late: ungated, energy, formation |
|---|---|---|---|---|---|---|---|
| 0.01 | 46, 2 | 5 [4, 10]; 3 [3, 3]; 3 [3, 3] | 1.423, 2.236, 1.000 | 1.524 | 6.3e-11, 2.2e-12 | 1.0011 | 7, 7, 7; 4, 4, 4 |
| 0.03 | 46, 2 | 5 [4, 10]; 3 [3, 3]; 3 [3, 3] | 1.421, 2.225, 1.000 | 1.526 | 7.7e-10, 1.3e-11 | 1.0009 | 6, 6, 6; 4, 4, 4 |
| 0.1 | 44, 2 | 5 [4, 10]; 3 [3, 3]; 3 [3, 3] | 1.392, 2.244, 1.000 | 1.552 | 2.2e-09, 3.9e-09 | 1.0007 | 6, 6, 6; 4, 4, 4 |

## C. The closed loop under measurement noise

D‴.49's Part E (k_e = 10, ω = 0.2, three bumps of width 0.12; 21 kernels from α = 0 at ε = 0.0001 to τ = 600; the noise 16 tones in 10–30 rad/s of rms s, the same phases); the dead zone on the error's energy filtered over 5 time units, the drive gated smoothly between the rms 1.2s and 1.5s; formed: |w| > 1% of max|A*|; outside: |c| > 1.2; the measured error's largest value over the last period, the true error x − x_m's rms there and the fit's rms error on the reference's range [−1, 1] (every kernel):

| regressor | s | adaptation | formed | outside | the formed kernels' coefficients | their centres | max\|ẽ\|, last period | true error rms, last period | fit rms on [−1, 1] | the gate closed at τ |
|---|---|---|---|---|---|---|---|---|---|---|
| x̃ | 0 | ungated | 5 | 0 | 0.900, -0.535, -0.065, 0.031, 0.720 | -0.520, 0.070, 0.070, 0.631, 0.609 | 0.0000 | 0.0000 | 0.0001 | — |
| x̃ | 0.01 | ungated | 14 | 14 | 0.854, 0.849, 0.785, 0.901, 0.659, 0.816, -0.568, -0.639, 0.727, 0.730, 0.486, 0.659, 0.707, 0.687 | 1.394, 1.394, 1.394, 1.394, 1.394, 1.394, -1.334, -1.333, 1.394, 1.394, 1.393, 1.394, 1.394, 1.394 | 0.1011 | 0.0372 | 0.4270 | — |
| x̃ | 0.01 | the mirror leak λΛ = 0.003 | 3 | 0 | 0.434, -0.159, 0.330 | -0.472, 0.170, 0.603 | 0.0624 | 0.0186 | 0.2616 | — |
| x̃ | 0.01 | the energy dead zone | 13 | 11 | 0.874, 0.738, 0.906, 0.957, 0.934, -0.531, -0.079, 0.737, 0.611, 0.717, 0.707, 0.719, 0.617 | 1.393, 1.392, 1.393, 1.393, 1.393, 0.050, 0.068, 1.392, 1.392, 1.392, 1.392, 1.392, 1.392 | 0.1010 | 0.0339 | 0.3813 | 216.0 |
| x̃ | 0.01 | the formation dead zone | 14 | 14 | 0.855, 0.847, 0.784, 0.901, 0.659, 0.818, -0.566, -0.619, 0.669, 0.729, 0.486, 0.659, 0.706, 0.687 | 1.394, 1.394, 1.394, 1.394, 1.393, 1.394, -1.333, -1.333, 1.393, 1.394, 1.393, 1.393, 1.393, 1.393 | 0.1011 | 0.0372 | 0.4271 | — |
| x̃ | 0.02 | ungated | 16 | 16 | 0.369, 0.444, 0.444, 0.540, 0.420, 0.398, -0.398, -0.343, -0.340, -0.312, 0.239, 0.284, 0.291, 0.302, 0.270, 0.207 | 1.328, 1.330, 1.330, 1.332, 1.329, 1.329, -1.293, -1.291, -1.291, -1.290, 1.325, 1.326, 1.326, 1.327, 1.326, 1.325 | 0.1170 | 0.0380 | 0.4273 | — |
| x̃ | 0.02 | the mirror leak λΛ = 0.003 | 3 | 0 | -0.145, 0.205, 0.028 | 0.142, 1.179, 0.686 | 0.1171 | 0.0358 | 0.4045 | — |
| x̃ | 0.02 | the energy dead zone | 16 | 16 | 0.357, 0.346, 0.358, 0.476, 0.464, 0.370, -0.400, -0.337, -0.335, -0.297, 0.221, 0.247, 0.222, 0.221, 0.244, 0.208 | 1.324, 1.324, 1.324, 1.327, 1.327, 1.324, -1.292, -1.290, -1.290, -1.288, 1.321, 1.321, 1.321, 1.321, 1.321, 1.320 | 0.1172 | 0.0380 | 0.4273 | — |
| x̃ | 0.02 | the formation dead zone | 16 | 16 | 0.406, 0.353, 0.410, 0.536, 0.408, 0.470, -0.399, -0.336, -0.330, -0.299, 0.222, 0.251, 0.246, 0.247, 0.243, 0.234 | 1.328, 1.326, 1.328, 1.330, 1.328, 1.329, -1.292, -1.290, -1.289, -1.288, 1.323, 1.324, 1.324, 1.324, 1.324, 1.324 | 0.1174 | 0.0380 | 0.4273 | — |
| x_m | 0 | ungated | 3 | 0 | 0.898, -0.600, 0.750 | -0.520, 0.070, 0.610 | 0.0000 | 0.0000 | 0.0000 | — |
| x_m | 0.01 | ungated | 3 | 0 | 0.930, -0.609, 0.765 | -0.527, 0.058, 0.648 | 0.0345 | 0.0071 | 0.0584 | — |
| x_m | 0.01 | the mirror leak λΛ = 0.003 | 3 | 0 | 0.447, -0.164, 0.340 | -0.470, 0.169, 0.564 | 0.0624 | 0.0187 | 0.2620 | — |
| x_m | 0.01 | the energy dead zone | 3 | 0 | 0.923, -0.435, 0.762 | -0.524, 0.043, 0.621 | 0.0339 | 0.0068 | 0.0622 | 180.0 |
| x_m | 0.01 | the formation dead zone | 3 | 0 | 0.930, -0.610, 0.765 | -0.526, 0.058, 0.648 | 0.0345 | 0.0071 | 0.0584 | 216.0 |
| x_m | 0.02 | ungated | 12 | 11 | 1.002, 1.023, 1.014, 1.036, 1.009, 0.984, -0.645, 1.036, 0.743, 0.813, 0.843, 1.039 | -1.495, 1.557, -1.515, 1.561, 1.552, 1.361, 0.009, -1.476, 1.547, 1.553, 1.525, 1.499 | 0.1170 | 0.0348 | 0.3884 | — |
| x_m | 0.02 | the mirror leak λΛ = 0.003 | 3 | 0 | 0.452, -0.165, 0.343 | -0.473, 0.169, 0.568 | 0.0835 | 0.0203 | 0.2590 | — |
| x_m | 0.02 | the energy dead zone | 7 | 3 | 0.884, 0.976, 0.995, -0.253, 0.670, 0.843, 0.010 | -0.522, -1.492, -1.473, 0.009, 0.618, 1.467, 0.614 | 0.0686 | 0.0135 | 0.1237 | 177.0 |
| x_m | 0.02 | the formation dead zone | 12 | 11 | 1.006, 0.998, 1.014, 1.028, 1.017, 1.012, -0.645, 1.035, 1.013, 0.807, 0.583, 0.784 | 1.539, -1.508, -1.522, 1.567, -1.480, 1.363, 0.009, 1.484, 1.562, 1.559, 1.527, 1.524 | 0.1175 | 0.0348 | 0.3884 | 177.0 |

Near the fit the loop high-passes the noise in the measured error (ẽ = e + n, e's response to n through −k_eẽ: |iω/(iω + k_e)|), a floor of rms 0.881s: the energy dead zone starts to close when the slow error's rms falls to 1.21s and is shut below 0.81s.
The ungated loops against the noise-free loop with the same regressor (the formed kernels in the order of their centres): x̃, s = 0.01: 14 kernels against 5; x_m, s = 0.01: max|Δw| 3.16e-02, max|Δc| 3.79e-02; x̃, s = 0.02: 16 kernels against 5; x_m, s = 0.02: 12 kernels against 3.

Numbers are produced by `experiments/exp_dead_zone.py`; nothing is edited by hand.
