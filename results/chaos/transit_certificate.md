# Theorem D‴.55: the transit's certificate (`exp_transit_certificate.py`)

The sliding inverse scale space (the free law's limit from small masses, moving centres) on random targets of three Gaussian features of the kernels' width; the certificate η(c_j) = Σ_k sign(A_k)⟨φ(c*_k)φ(c_j)⟩/⟨φ(c*_k)φ(c_{e_k})⟩ on the starting points (the samples' inner product, e_k the nearest starting point to feature k, E their set), its margin m = 1 − max_{j∉E}|η(c_j)| and the cross-talk q = max_k Σ_{l≠k}G((Δ_kl − 2R)₊)/G(R). The path's outcome: exact with E (its active set E, the residual below 1e-08, no exits), exact with other kernels (three kernels, the residual below 1e-08), or not exact.

## A. One dimension

80 targets per minimal separation (seed 551; σ = 0.12, the grid's spacing 0.1); amplitudes ±U(0.4, 1):

| Δ ≥ | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | the disagreements' m (q) |
|---|---|---|---|---|---|
| 0.25 | 80 | 26, 17, 0 | 6, 31, 0 | 0.343 | +0.0738 (0.60), -0.1442 (0.32), +0.1132 (0.26), +0.1949 (0.85), +0.0459 (0.56), +0.1369 (0.47), +0.1346 (0.54), -0.1486 (0.22), -0.0253 (0.13), +0.0880 (0.42), +0.1341 (0.38), +0.0814 (0.34), +0.0777 (0.25), +0.0566 (0.40), +0.1429 (0.42), +0.0832 (0.51), -0.0278 (0.08), -0.0286 (0.15), +0.0279 (0.50), -0.0327 (0.17), +0.1500 (0.61), +0.0897 (0.38), +0.0643 (0.18) |
| 0.3 | 80 | 35, 13, 0 | 6, 26, 0 | 0.174 | -0.0248 (0.11), +0.0351 (0.57), +0.0243 (0.06), +0.0619 (0.55), -0.0789 (0.10), -0.1248 (0.25), -0.0017 (0.03), +0.0343 (0.17), +0.0204 (0.15), +0.0301 (0.21), +0.0737 (0.39), +0.0104 (0.09), -0.0284 (0.09), -0.0583 (0.15), +0.0903 (0.37), +0.0365 (0.31), +0.0803 (0.45), +0.0429 (0.14), +0.1167 (0.31) |
| 0.35 | 80 | 37, 5, 0 | 9, 29, 0 | 0.151 | -0.0180 (0.10), -0.0094 (0.08), +0.0352 (0.12), -0.0078 (0.20), -0.0359 (0.11), -0.1760 (0.36), -0.0751 (0.19), -0.1770 (0.31), +0.0131 (0.11), +0.0233 (0.21), -0.0575 (0.13), -0.0328 (0.10), +0.1083 (0.24), +0.0487 (0.21) |
| 0.45 | 80 | 60, 4, 0 | 5, 11, 0 | 0.048 | +0.0150 (0.11), -0.0003 (0.06), +0.0074 (0.14), +0.0248 (0.09), -0.0077 (0.07), +0.0127 (0.08), -0.0015 (0.06), -0.0214 (0.13), -0.0015 (0.12) |
| 0.6 | 80 | 76, 2, 0 | 2, 0, 0 | 0.008 | -0.0005 (0.01), +0.0016 (0.01), -0.0007 (0.02), +0.0014 (0.01) |

Over all 400 targets the certificate's sign predicts exactness with E on 331; the disagreements' margins |m| ≤ 0.1949, at most 0.761 times their cross-talk q. Outside the band |m| ≤ Cq the sign decides on C = 0.5: 208 of 211; C = 0.8: 135 of 135; C = 1: 122 of 122. Every path ends at an exact fit with 3 kernels on 400 of 400 targets.

## B. The amplitudes' range

40 configurations (seed 552; positions at least 0.35 apart, signs) with |A_k| = 10^{−r·u_k}, r = 0, 1, 2, 3 (the amplitudes' range up to 10^3); the certificate depends on the positions and signs only (21 of 40 with m > 0). The outcome is the same at every range on 29 of 40 configurations, while the order of the entries changes with the range on 36; it changes on 11: configuration 5 (m = -0.1052): other, other, other, nearest; configuration 7 (m = -0.1488): other, other, other, nearest; configuration 11 (m = -0.1783): other, nearest, nearest, nearest; configuration 12 (m = +0.0997): other, nearest, nearest, nearest; configuration 14 (m = -0.1265): nearest, nearest, other, nearest; configuration 15 (m = -0.1138): other, nearest, nearest, nearest; configuration 16 (m = -0.0197): nearest, nearest, other, nearest; configuration 23 (m = -0.0277): other, other, nearest, nearest; configuration 24 (m = -0.0007): other, nearest, nearest, nearest; configuration 28 (m = +0.0240): other, nearest, other, other; configuration 33 (m = -0.0828): other, other, nearest, nearest. The configurations whose outcome changes have |m|/q ≤ 0.67; those with |m| > 0.8q (6) keep their outcome at every range.

## C. The plane

40 targets per minimal separation (seed 553; σ = 0.15, the lattice's spacing 0.125; centres in [−0.7, 0.7]²), amplitudes ±U(0.4, 1):

| Δ ≥ | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | the disagreements' m (q) |
|---|---|---|---|---|---|
| 0.3 | 40 | 14, 12, 0 | 5, 9, 0 | 0.215 | +0.0380 (0.31), -0.0111 (0.11), +0.0573 (0.22), +0.0736 (0.39), -0.0524 (0.29), +0.0185 (0.61), -0.0199 (0.16), +0.1001 (0.41), -0.0226 (0.24), +0.1289 (0.64), +0.0145 (0.68), +0.1268 (0.31), +0.1368 (0.52), -0.0225 (0.17), +0.0392 (0.46), +0.0192 (0.90), +0.0128 (0.19) |
| 0.45 | 40 | 26, 6, 0 | 4, 4, 0 | 0.107 | -0.0297 (0.10), -0.0907 (0.27), +0.0841 (0.48), +0.0124 (0.06), +0.0150 (0.10), -0.0317 (0.19), +0.0116 (0.13), +0.0142 (0.20), +0.0393 (0.30), -0.1624 (0.35) |
| 0.6 | 40 | 33, 3, 0 | 1, 3, 0 | 0.028 | +0.0189 (0.08), +0.0075 (0.07), +0.0064 (0.07), -0.0049 (0.07) |

Over all 120 targets the certificate's sign predicts exactness with E on 89. Outside the band |m| ≤ Cq the sign decides on C = 0.5: 52 of 52; C = 0.8: 41 of 41; C = 1: 34 of 34. Every path ends at an exact fit with 3 kernels on 120 of 120 targets.

Numbers are produced by `experiments/exp_transit_certificate.py`; nothing is edited by hand.
