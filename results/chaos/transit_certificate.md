# Theorem D‴.55: the transit's certificate (`exp_transit_certificate.py`)

The sliding inverse scale space (the free law's limit from small masses, moving centres) on random targets of three Gaussian features of the kernels' width; the certificate η(c_j) = Σ_k sign(A_k)⟨φ(c*_k)φ(c_j)⟩/⟨φ(c*_k)φ(c_{e_k})⟩ on the starting points (the samples' inner product, e_k the starting point of largest correlation with feature k — the nearest in the continuum —, E their set), its margin m = 1 − max_{j∉E}|η(c_j)| and the cross-talk q = max_k Σ_{l≠k}G((Δ_kl − 2R)₊)/G(R). The path's outcome: exact with E (its active set E, the residual below 1e-08, no exits), exact with other kernels (three kernels, the residual below 1e-08), or not exact. The paths by `sliding_iss_nd` (in one dimension on the 1-D samples, with the exact Jacobian), the block polished by its least-squares fit unless a kernel below 10⁻² of the block's largest coefficient would land beyond one kernel width of both its start and its correlation's ascent's end (a rejected polish; the bound does not check single-kernel blocks or the block's largest kernel); the skeletons recomputed from the targets.

## A. One dimension

80 targets per minimal separation (seed 551; σ = 0.12, the grid's spacing 0.1); amplitudes ±U(0.4, 1):

| Δ ≥ | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | the disagreements' m (q) |
|---|---|---|---|---|---|
| 0.25 | 80 | 26, 17, 0 | 6, 31, 0 | 0.343 | +0.0738 (0.60), -0.1442 (0.32), +0.1132 (0.26), +0.1949 (0.85), +0.0459 (0.56), +0.1369 (0.47), +0.1346 (0.54), -0.1486 (0.22), -0.0253 (0.13), +0.0880 (0.42), +0.1341 (0.38), +0.0814 (0.34), +0.0777 (0.25), +0.0566 (0.40), +0.1429 (0.42), +0.0832 (0.51), -0.0278 (0.08), -0.0286 (0.15), +0.0279 (0.50), -0.0327 (0.17), +0.1500 (0.61), +0.0897 (0.38), +0.0643 (0.18) |
| 0.3 | 80 | 35, 13, 0 | 6, 26, 0 | 0.174 | -0.0248 (0.11), +0.0351 (0.57), +0.0243 (0.06), +0.0619 (0.55), -0.0789 (0.10), -0.1248 (0.25), -0.0017 (0.03), +0.0343 (0.17), +0.0204 (0.15), +0.0301 (0.21), +0.0737 (0.39), +0.0104 (0.09), -0.0284 (0.09), -0.0583 (0.15), +0.0903 (0.37), +0.0365 (0.31), +0.0803 (0.45), +0.0429 (0.14), +0.1167 (0.31) |
| 0.35 | 80 | 37, 5, 0 | 9, 29, 0 | 0.151 | -0.0180 (0.10), -0.0094 (0.08), +0.0352 (0.12), -0.0078 (0.20), -0.0359 (0.11), -0.1760 (0.36), -0.0751 (0.19), -0.1770 (0.31), +0.0131 (0.11), +0.0233 (0.21), -0.0575 (0.13), -0.0328 (0.10), +0.1083 (0.24), +0.0487 (0.21) |
| 0.45 | 80 | 60, 4, 0 | 5, 11, 0 | 0.048 | +0.0150 (0.11), -0.0003 (0.06), +0.0074 (0.14), +0.0248 (0.09), -0.0077 (0.07), +0.0127 (0.08), -0.0015 (0.06), -0.0214 (0.13), -0.0015 (0.12) |
| 0.6 | 80 | 76, 2, 0 | 2, 0, 0 | 0.008 | -0.0005 (0.01), +0.0016 (0.01), -0.0007 (0.02), +0.0014 (0.01) |

Over all 400 targets the certificate's sign predicts exactness with E on 331; the disagreements' margins |m| ≤ 0.1949, at most 0.7608 times their cross-talk q. Outside the band |m| ≤ Cq the sign decides on C = 0.5: 208 of 211; C = 0.8: 135 of 135; C = 1: 122 of 122. Every path ends at an exact fit with 3 kernels on 400 of 400 targets.

Exactness with E predicted by the sign of each margin (the relative margin min_k s_kη(c_{e_k}) − max_{j∉E}|η(c_j)|; the single-feature margin min_k(1 − max_{j≠e_k}⟨φ_kφ_j⟩/⟨φ_kφ_{e_k}⟩)), of 400: the certificate, 1 − max off E 331, the relative margin 344, the single feature (no signs or sums) 262, the atoms subtracted (no tangent directions) 390, the path's skeleton 397, the fitted skeleton 400, always E 262, E iff the first entrant is in E 327 (the first entrant off E — the outcome decided at the first entry — on 65; the fitted skeleton decided at its first phase on 101); the unpolished flow's outcome the polished path's on 400 of 400 (its blocks relaxed on 400).

## B. The amplitudes' range

40 configurations (seed 552; positions at least 0.35 apart, signs) with |A_k| = 10^{−r·u_k}, r = 0, 1, 2, 3 (the amplitudes' range up to 10^3); the certificate depends on the positions and signs only (21 of 40 with m > 0). The outcome is the same at every range on 31 of 40 configurations, while the order of the entries changes with the range on 35; it changes on 9: configuration 5 (m = -0.1052): other, other, other, nearest; configuration 7 (m = -0.1488): other, other, other, nearest; configuration 11 (m = -0.1783): other, nearest, nearest, nearest; configuration 12 (m = +0.0997): other, nearest, nearest, nearest; configuration 15 (m = -0.1138): other, nearest, nearest, nearest; configuration 23 (m = -0.0277): other, other, nearest, nearest; configuration 24 (m = -0.0007): other, nearest, nearest, nearest; configuration 28 (m = +0.0240): other, nearest, nearest, nearest; configuration 33 (m = -0.0828): other, other, nearest, nearest. The configurations whose outcome changes have |m|/q ≤ 0.67; of those with |m| > 0.8q, 6 of 6 keep their outcome at every range. The ranges realised, max|A|/min|A|, median [max]: r = 0: 1 [1], r = 1: 3.22 [9.08], r = 2: 10.4 [82.4], r = 3: 33.9 [749]; the amplitudes' order is the same for every r > 0 (u_k fixed), so this Part varies the range, not the order (D(c) varies the order).

Exactness with E predicted by the margins' signs, of 160: the certificate, 1 − max off E 119, the relative margin 107, the single feature (no signs or sums) 113, the atoms subtracted (no tangent directions) 140, the path's skeleton 160, the fitted skeleton 160, always E 113, E iff the first entrant is in E 125 (the first entrant off E — the outcome decided at the first entry — on 12; the fitted skeleton decided at its first phase on 34).

## C. The plane

40 targets per minimal separation (seed 553; σ = 0.15, the lattice's spacing 0.125; centres in [−0.7, 0.7]²), amplitudes ±U(0.4, 1):

| Δ ≥ | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | the disagreements' m (q) |
|---|---|---|---|---|---|
| 0.3 | 40 | 14, 12, 0 | 5, 9, 0 | 0.215 | +0.0380 (0.31), -0.0111 (0.11), +0.0573 (0.22), +0.0736 (0.39), -0.0524 (0.29), +0.0185 (0.61), -0.0199 (0.16), +0.1001 (0.41), -0.0226 (0.24), +0.1289 (0.64), +0.0145 (0.68), +0.1268 (0.31), +0.1368 (0.52), -0.0225 (0.17), +0.0392 (0.46), +0.0192 (0.90), +0.0128 (0.19) |
| 0.45 | 40 | 26, 6, 0 | 4, 4, 0 | 0.107 | -0.0297 (0.10), -0.0907 (0.27), +0.0841 (0.48), +0.0124 (0.06), +0.0150 (0.10), -0.0317 (0.19), +0.0116 (0.13), +0.0142 (0.20), +0.0393 (0.30), -0.1624 (0.35) |
| 0.6 | 40 | 33, 3, 0 | 1, 3, 0 | 0.028 | +0.0189 (0.08), +0.0075 (0.07), +0.0064 (0.07), -0.0049 (0.07) |

Over all 120 targets the certificate's sign predicts exactness with E on 89. Outside the band |m| ≤ Cq the sign decides on C = 0.5: 52 of 52; C = 0.8: 41 of 41; C = 1: 34 of 34. Every path ends at an exact fit with 3 kernels on 120 of 120 targets.

Exactness with E predicted by the margins' signs, of 120: the certificate, 1 − max off E 89, the relative margin 102, the single feature (no signs or sums) 83, the atoms subtracted (no tangent directions) 117, the path's skeleton 120, the fitted skeleton 120, always E 83, E iff the first entrant is in E 98 (the first entrant off E — the outcome decided at the first entry — on 15; the fitted skeleton decided at its first phase on 23).

## D. Designed families (referee round 43)

(a) Two features on starting points (1e-06 off) at −0.1 and −0.1 + D, the first's amplitude 1, the second's ±ρ_A: the margins (amplitude-free) and the outcome at each ρ_A (E exact with E, O exact with others, × not exact, t over 600 s):

| D | signs | m | m/q | relative | single-feature | 0.001 | 0.02 | 0.1 | 0.333 | 1 | 3 | 10 | 50 | 1e+03 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.3 | + + | -0.3400 | -1.622 | -0.1304 | +0.1594 | E | E | O | O | O | O | O | E | E |
| 0.3 | + − | +0.2216 | +1.057 | +0.0119 | +0.1594 | O | O | O | O | E | O | O | O | O |
| 0.4 | + + | -0.0502 | -0.808 | +0.0119 | +0.1594 | E | E | E | E | E | E | E | E | E |
| 0.4 | + − | +0.1724 | +2.773 | +0.1102 | +0.1594 | E | E | E | E | E | E | E | E | E |
| 0.5 | + + | +0.0972 | +7.458 | +0.1102 | +0.1594 | E | E | E | E | E | E | E | E | E |
| 0.5 | + − | +0.1613 | +12.376 | +0.1483 | +0.1594 | E | E | E | E | E | E | E | E | E |

(b) Two features at s·σ, equal amplitudes, the midpoint uniform over a grid cell (20 per row): the outcomes and how often each margin's sign is right:

| s | signs | E, O, ×, t | m > 0 | the certificate right | relative right | single-feature right | atoms subtracted right | the skeleton right | the fitted skeleton right | always E right |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.5 | + + | 5, 15, 0, 0 | 0 | 15 | 20 | 5 | 20 | 20 | 20 | 5 |
| 1.5 | + − | 0, 20, 0, 0 | 20 | 0 | 20 | 0 | 20 | 20 | 20 | 0 |
| 2 | + + | 0, 20, 0, 0 | 0 | 20 | 20 | 0 | 20 | 20 | 20 | 0 |
| 2 | + − | 2, 18, 0, 0 | 20 | 2 | 20 | 2 | 20 | 20 | 20 | 2 |
| 2.5 | + + | 0, 20, 0, 0 | 0 | 20 | 20 | 0 | 20 | 20 | 20 | 0 |
| 2.5 | + − | 1, 19, 0, 0 | 20 | 1 | 20 | 1 | 20 | 20 | 20 | 1 |
| 3 | + + | 0, 20, 0, 0 | 0 | 20 | 20 | 0 | 20 | 20 | 20 | 0 |
| 3 | + − | 15, 5, 0, 0 | 20 | 15 | 20 | 15 | 20 | 20 | 20 | 15 |
| 3.5 | + + | 12, 8, 0, 0 | 0 | 8 | 19 | 12 | 20 | 20 | 20 | 12 |
| 3.5 | + − | 16, 4, 0, 0 | 20 | 16 | 20 | 16 | 19 | 20 | 20 | 16 |
| 4 | + + | 17, 3, 0, 0 | 15 | 18 | 20 | 17 | 20 | 20 | 20 | 17 |
| 4 | + − | 16, 4, 0, 0 | 20 | 16 | 20 | 16 | 20 | 20 | 20 | 16 |

(c) 30 of A's targets with |m| > 0.8q, the amplitude magnitudes 1, 10^{−1/2}, 10^{−1} in all six orders (the signs kept): the outcome changes with the order on 0 of them.

Over D, exactness with E predicted by the margins' signs: (a) of 54: the certificate, 1 − max off E 33, the relative margin 42, the single feature (no signs or sums) 41, the atoms subtracted (no tangent directions) 34, the path's skeleton 52, the fitted skeleton 54, always E 41, E iff the first entrant is in E 42 (the first entrant off E — the outcome decided at the first entry — on 1; the fitted skeleton decided at its first phase on 1); (b) of 240: the certificate, 1 − max off E 151, the relative margin 239, the single feature (no signs or sums) 84, the atoms subtracted (no tangent directions) 239, the path's skeleton 240, the fitted skeleton 240, always E 84, E iff the first entrant is in E 205 (the first entrant off E — the outcome decided at the first entry — on 121; the fitted skeleton decided at its first phase on 125); the unpolished flow's outcome the polished path's on 240 of 240 (its blocks relaxed on 240); (c) of 180: the certificate, 1 − max off E 180, the relative margin 180, the single feature (no signs or sums) 180, the atoms subtracted (no tangent directions) 180, the path's skeleton 180, the fitted skeleton 180, always E 180, E iff the first entrant is in E 180 (the first entrant off E — the outcome decided at the first entry — on 0; the fitted skeleton decided at its first phase on 29). Paths over 600 s: 0; polishes rejected: 1 on 1 paths.

D(b)'s failures that are not the second feature's own near ties (m₀ > 0.05), per row: 1.5σ + + 11, 1.5σ + − 11, 2σ + + 5, 2σ + − 9, 2.5σ + + 16, 2.5σ + − 15, 3σ + + 9, 3σ + − 0, 3.5σ + + 0, 3.5σ + − 3 (m₀ = 0.054, 0.056, 0.055), 4σ + + 0, 4σ + − 0.

Outside the band, the certificate's sign on the designed targets: |m| > 0.5q: wrong on 47 of 399 (|m|/q up to 1.622, 34 of them with m > 0); |m| > 0.8q: wrong on 25 of 333 (|m|/q up to 1.622, 12 of them with m > 0); |m| > 1q: wrong on 12 of 271 (|m|/q up to 1.622, 8 of them with m > 0).

## E. The two-feature map

The first feature (amplitude 1) on the starting point −0.1 (1e-06 off), the second at the distance D = 0.12 … 0.6 (step 0.01; σ = 0.12) with the amplitude ±10^ℓ, ℓ = -3 … 3 (step 0.125): the fitted skeleton's outcome at every point (E: exact with the nearest kernels, ·: not), shown at every 0.06 in D and every ½ in ℓ:

Same signs:

| D/σ \ ℓ | -3 | -2.5 | -2 | -1.5 | -1 | -0.5 | 0 | 0.5 | 1 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 1.50 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 2.00 | E | E | E | E | E | · | · | · | · | · | · | · | · |
| 2.50 | E | E | E | · | · | · | · | · | · | · | E | E | E |
| 3.00 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 3.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 5.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |

Opposite signs:

| D/σ \ ℓ | -3 | -2.5 | -2 | -1.5 | -1 | -0.5 | 0 | 0.5 | 1 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.00 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 1.50 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 2.00 | · | · | · | · | · | · | · | · | · | · | · | · | · |
| 2.50 | · | · | · | · | · | · | E | · | · | · | · | · | · |
| 3.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 3.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 4.50 | E | E | E | E | E | E | E | E | E | E | E | E | E |
| 5.00 | E | E | E | E | E | E | E | E | E | E | E | E | E |

The fitted skeleton's threshold D*(ℓ) on this grid (the first feature on a starting point, |ℓ| ≤ 3: not (iii)'s Δ₀) — every D of the grid at or beyond it exact with E —, in units of σ: same signs ℓ = -3: 4.67, ℓ = -2: 4.67, ℓ = -1: 4.67, ℓ = -0.5: 4.67, ℓ = 0: 4.67, ℓ = 0.5: 4.67, ℓ = 1: 4.67, ℓ = 2: 4.67, ℓ = 3: 4.67; opposite signs ℓ = -3: 3.75, ℓ = -2: 3.75, ℓ = -1: 3.75, ℓ = -0.5: 3.75, ℓ = 0: 3.75, ℓ = 0.5: 2.92, ℓ = 1: 2.75, ℓ = 2: 2.75, ℓ = 3: 2.83; over all ℓ at most 4.67σ. Beyond 3.25σ every failure is the second feature's own near tie (its single-feature margin m₀ ≤ 0.05, near a midpoint of the grid: D = 0.44, 0.45, 0.46, 0.55); without them D*(ℓ) is same signs ℓ = -3: 3.17, ℓ = -2: 3.25, ℓ = -1: 3.25, ℓ = -0.5: 3.25, ℓ = 0: 3.33, ℓ = 0.5: 3.17, ℓ = 1: 2.67, ℓ = 2: 2.50, ℓ = 3: 2.50; opposite signs ℓ = -3: 2.83, ℓ = -2: 2.83, ℓ = -1: 2.83, ℓ = -0.5: 2.83, ℓ = 0: 2.83, ℓ = 0.5: 2.67, ℓ = 1: 2.75, ℓ = 2: 2.75, ℓ = 3: 2.83, at most 3.33σ.

The two ends of the ratio ((v)): same signs: constant for ℓ ≤ −2 at 48 and for ℓ ≥ 2 at 49 of 49 separations, the two ends' outcomes different at 7; opposite signs: constant for ℓ ≤ −2 at 49 and for ℓ ≥ 2 at 48 of 49 separations, the two ends' outcomes different at 5.

The margins' signs against the fitted skeleton on the map's 4802 points: the certificate, 1 − max off E 3128, the relative margin 4227, the single feature (no signs or sums) 2754, the atoms subtracted (no tangent directions) 3991, the path's skeleton 4719; exact with E on 2754; the linear skeleton's |m_sk| > 0.01 where its sign differs from the fitted skeleton's at 33 points (|m_sk| up to 0.134).

(ii)'s order: |m_sk − m_fit| against the ratio of the smaller feature's amplitude to the larger's, its log–log slope over 10⁻³ … 10⁻¹ (1: first order in the tail, 2: second): 2.50σ + + second small 1.00, 2.50σ + + first small 1.00, 3.50σ + + second small 0.99, 3.50σ + + first small 0.99, 4.00σ + + second small 0.99, 4.00σ + + first small 0.98, 3.50σ + − second small 0.98, 4.00σ + − second small 0.98.

The path at 240 sampled points (|ℓ| ≤ 2.5; under 300 s each): it agrees with the fitted skeleton at 240 of 240 (0 over the limit); exact with E at 134; polishes rejected 14 on 3 paths.

## F. A wide range

Three features at −s, 0, s (+0.013; s = 0.6, 0.36: 5σ, 3σ), the amplitudes 1, 0.001, 1e-06 in all six orders, the signs + + + and + − +; exits only below 10⁻¹² of the largest coefficient: the outcome per order (E exact with E, O with others, × not, t over the limit), the polishes rejected, the skeletons' signs:

| s | signs | outcomes (orders 012, 021, 102, 120, 201, 210) | polishes rejected | the skeleton right | the fitted skeleton right |
|---|---|---|---|---|---|
| 5σ | + + + | EEEEEE | 0 | 6 of 6 | 6 of 6 |
| 5σ | + − + | EEEEEE | 0 | 6 of 6 | 6 of 6 |
| 3σ | + + + | O×OOOO | 0 | 6 of 6 | 6 of 6 |
| 3σ | + − + | OEOEOO | 0 | 6 of 6 | 6 of 6 |

Numbers are produced by `experiments/exp_transit_certificate.py`; nothing is edited by hand.
