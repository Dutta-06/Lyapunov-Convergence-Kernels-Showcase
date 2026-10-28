# Proposition D‴.54: the transit in ℝⁿ (`exp_transit_nd.py`)

The free law with moving centres (m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025), Gaussian kernels of width 0.15 in the plane on a square lattice of spacing 0.125 (H: a hexagonal one), the samples on a grid of spacing 0.05; masses ε; formed: |w| > 0.01|A|; κ = πρ(2c)^{1/2} = 0.7025.

## A. The k-fold near tie

One feature near an edge (the kernels at (0, 0), (0.125, 0)) or a vertex (and (0.125, 0), (0, 0.125), (0.125, 0.125)) of the lattice — 8 × 8 kernels around it, the samples symmetric about the vertex —, its margins ε^{1/2}x at the entry; the loop to 1.32 τ_entry. The gaps x (the runner-ups', as placed; at the vertex the diagonal's as it follows), the shares of the tie's kernels in the loop against the inner problem's at the gaps of the smallest ε (the feature at the tie; extrapolated in its start from xi0 = 800, 1600, the change between those in the last column), the softmax at those gaps and, at the vertex, the product of the two edges' inner problems posed in the vertex's geometry (the feature H/2 off each edge's axis):

| tie | A | x | shares: loop (ε = 1e-06, 1e-07, 1e-08, 1e-09, 1e-10) | inner problem | softmax | the edges' product | the others' largest \|w\| (ε = 1e-06, 1e-07, 1e-08, 1e-09, 1e-10) | the loop's successive changes | xi0's change |
|---|---|---|---|---|---|---|---|---|---|
| edge | 0.3 | 0.000 | 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5 | 0.5, 0.5 | 0.5, 0.5 | — | 2.50e-08, 7.63e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 1.4e-13, 1.6e-14, 9.5e-13, 5.2e-12 | 1.1e-12 |
| edge | 0.3 | 0.500 | 0.5886, 0.4114; 0.59, 0.41; 0.5905, 0.4095; 0.5906, 0.4094; 0.5906, 0.4094 | 0.5907, 0.4093 | 0.5869, 0.4131 | — | 2.51e-08, 7.64e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 1.4e-03, 4.5e-04, 1.4e-04, 4.6e-05 | 1.6e-06 |
| edge | 0.3 | 1.000 | 0.6718, 0.3282; 0.6744, 0.3256; 0.6752, 0.3248; 0.6754, 0.3246; 0.6755, 0.3245 | 0.6756, 0.3244 | 0.6687, 0.3313 | — | 2.52e-08, 7.64e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 2.5e-03, 8.1e-04, 2.6e-04, 8.2e-05 | 2.8e-06 |
| edge | 0.3 | 2.000 | 0.8074, 0.1926; 0.8109, 0.1891; 0.812, 0.188; 0.8123, 0.1877; 0.8125, 0.1875 | 0.8125, 0.1875 | 0.803, 0.197 | — | 2.53e-08, 7.65e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 3.5e-03, 1.1e-03, 3.5e-04, 1.1e-04 | 3.9e-06 |
| edge | 0.3 | 4.000 | 0.946, 0.05402; 0.9481, 0.0519; 0.9488, 0.05123; 0.949, 0.05102; 0.949, 0.05095 | 0.9491, 0.05092 | 0.9432, 0.05679 | — | 2.54e-08, 7.66e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 2.1e-03, 6.7e-04, 2.1e-04, 6.7e-05 | 2.5e-06 |
| edge | 0.3 | 8.000 | 0.9958, 0.004171; 0.9961, 0.003929; 0.9961, 0.003855; 0.9962, 0.003832; 0.9962, 0.003824 | 0.9962, 0.003821 | 0.9964, 0.003612 | — | 2.54e-08, 7.66e-10, 2.39e-11, 7.52e-13, 2.38e-14 | 2.4e-04, 7.4e-05, 2.3e-05, 7.4e-06 | 3.0e-07 |
| edge | 1 | 0.000 | 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5; 0.5, 0.5 | 0.5, 0.5 | 0.5, 0.5 | — | 2.57e-08, 7.69e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 8.0e-15, 6.6e-14, 2.1e-12, 4.3e-13 | 9.3e-13 |
| edge | 1 | 0.500 | 0.5886, 0.4114; 0.59, 0.41; 0.5905, 0.4095; 0.5906, 0.4094; 0.5907, 0.4093 | 0.5907, 0.4093 | 0.5869, 0.4131 | — | 2.58e-08, 7.70e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 1.4e-03, 4.5e-04, 1.4e-04, 4.6e-05 | 1.6e-06 |
| edge | 1 | 1.000 | 0.6719, 0.3281; 0.6744, 0.3256; 0.6752, 0.3248; 0.6755, 0.3245; 0.6756, 0.3244 | 0.6756, 0.3244 | 0.6687, 0.3313 | — | 2.58e-08, 7.71e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 2.5e-03, 8.1e-04, 2.6e-04, 8.2e-05 | 2.8e-06 |
| edge | 1 | 2.000 | 0.8075, 0.1925; 0.811, 0.189; 0.8121, 0.1879; 0.8125, 0.1875; 0.8126, 0.1874 | 0.8126, 0.1874 | 0.803, 0.197 | — | 2.60e-08, 7.72e-10, 2.39e-11, 7.53e-13, 2.38e-14 | 3.5e-03, 1.1e-03, 3.5e-04, 1.1e-04 | 3.9e-06 |
| edge | 1 | 4.000 | 0.9464, 0.05364; 0.9485, 0.0515; 0.9492, 0.05083; 0.9494, 0.05062; 0.9494, 0.05055 | 0.9495, 0.05052 | 0.9432, 0.05679 | — | 2.61e-08, 7.72e-10, 2.40e-11, 7.53e-13, 2.38e-14 | 2.1e-03, 6.7e-04, 2.1e-04, 6.7e-05 | 2.5e-06 |
| edge | 1 | 8.000 | 0.9966, 0.003442; 0.9968, 0.003211; 0.9969, 0.003141; 0.9969, 0.003118; 0.9969, 0.003111 | 0.9969, 0.003108 | 0.9964, 0.003612 | — | 2.61e-08, 7.73e-10, 2.40e-11, 7.53e-13, 2.38e-14 | 2.3e-04, 7.1e-05, 2.2e-05, 7.0e-06 | 2.8e-07 |
| vertex | 0.3 | 0.000, 0.000, 0.000 | 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 3.83e-09, 1.20e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 1.9e-13, 4.1e-13, 9.0e-13, 2.2e-12 | 7.8e-13 |
| vertex | 0.3 | 0.500, 0.500, 1.000 | 0.3509, 0.2415, 0.2415, 0.1662; 0.3526, 0.2412, 0.2412, 0.165; 0.3532, 0.2411, 0.2411, 0.1646; 0.3534, 0.2411, 0.2411, 0.1645; 0.3534, 0.2411, 0.2411, 0.1645 | 0.3535, 0.241, 0.241, 0.1644 | 0.3445, 0.2424, 0.2424, 0.1706 | 0.3535, 0.2411, 0.2411, 0.1644 | 3.86e-09, 1.21e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 1.8e-03, 5.7e-04, 1.8e-04, 5.8e-05 | 3.9e-06 |
| vertex | 0.3 | 1.000, 1.000, 2.000 | 0.4605, 0.2181, 0.2181, 0.1034; 0.4641, 0.2171, 0.2171, 0.1017; 0.4653, 0.2168, 0.2168, 0.1012; 0.4657, 0.2167, 0.2167, 0.101; 0.4658, 0.2166, 0.2166, 0.1009 | 0.4658, 0.2166, 0.2166, 0.1009 | 0.4472, 0.2215, 0.2215, 0.1097 | 0.466, 0.2166, 0.2166, 0.1007 | 3.88e-09, 1.21e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 3.7e-03, 1.2e-03, 3.7e-04, 1.2e-04 | 8.0e-06 |
| vertex | 0.3 | 1.000, 3.000, 4.000 | 0.613, 0.2905, 0.06527, 0.0313; 0.6175, 0.2889, 0.06347, 0.03012; 0.619, 0.2884, 0.06289, 0.02975; 0.6195, 0.2882, 0.06271, 0.02963; 0.6196, 0.2881, 0.06266, 0.0296 | 0.6197, 0.2881, 0.06263, 0.02958 | 0.5963, 0.2954, 0.07247, 0.0359 | 0.6202, 0.2883, 0.06247, 0.02904 | 3.95e-09, 1.22e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 4.6e-03, 1.5e-03, 4.6e-04, 1.5e-04 | 9.9e-06 |
| vertex | 0.3 | 2.000, 2.000, 4.000 | 0.6668, 0.1496, 0.1496, 0.03399; 0.6728, 0.1473, 0.1473, 0.03274; 0.6747, 0.1465, 0.1465, 0.03234; 0.6753, 0.1463, 0.1463, 0.03221; 0.6755, 0.1462, 0.1462, 0.03218 | 0.6756, 0.1461, 0.1461, 0.03216 | 0.6448, 0.1582, 0.1582, 0.03882 | 0.676, 0.1462, 0.1462, 0.03161 | 3.92e-09, 1.21e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 6.0e-03, 1.9e-03, 6.1e-04, 1.9e-04 | 1.3e-05 |
| vertex | 0.3 | 2.000, 5.000, 7.000 | 0.7962, 0.179, 0.01962, 0.005226; 0.801, 0.1754, 0.01868, 0.004926; 0.8026, 0.1742, 0.01838, 0.004834; 0.803, 0.1739, 0.01829, 0.004805; 0.8032, 0.1738, 0.01826, 0.004796 | 0.8033, 0.1737, 0.01825, 0.004791 | 0.7797, 0.1913, 0.02326, 0.005707 | 0.8042, 0.1739, 0.01807, 0.003907 | 4.02e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.8e-03, 1.5e-03, 4.9e-04, 1.5e-04 | 1.0e-05 |
| vertex | 0.3 | 4.000, 4.000, 8.000 | 0.9049, 0.04594, 0.04594, 0.003234; 0.9091, 0.04395, 0.04395, 0.003035; 0.9104, 0.04333, 0.04333, 0.002974; 0.9108, 0.04313, 0.04313, 0.002955; 0.9109, 0.04307, 0.04307, 0.002949 | 0.911, 0.04304, 0.04304, 0.002946 | 0.8896, 0.05357, 0.05357, 0.003225 | 0.9119, 0.04304, 0.04304, 0.002032 | 4.00e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.2e-03, 1.3e-03, 4.1e-04, 1.3e-04 | 8.7e-06 |
| vertex | 1 | 0.000, 0.000, 0.000 | 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25; 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 0.25, 0.25, 0.25, 0.25 | 3.87e-09, 1.21e-10, 3.80e-12, 1.20e-13, 3.80e-15 | 1.3e-13, 5.7e-13, 8.1e-13, 4.4e-12 | 9.6e-13 |
| vertex | 1 | 0.500, 0.500, 1.000 | 0.3509, 0.2415, 0.2415, 0.1661; 0.3527, 0.2412, 0.2412, 0.1649; 0.3533, 0.2411, 0.2411, 0.1645; 0.3535, 0.2411, 0.2411, 0.1644; 0.3535, 0.2411, 0.2411, 0.1644 | 0.3536, 0.241, 0.241, 0.1643 | 0.3445, 0.2424, 0.2424, 0.1706 | 0.3536, 0.241, 0.241, 0.1643 | 3.90e-09, 1.21e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 1.8e-03, 5.7e-04, 1.8e-04, 5.8e-05 | 3.9e-06 |
| vertex | 1 | 1.000, 1.000, 2.000 | 0.4607, 0.2181, 0.2181, 0.1032; 0.4643, 0.2171, 0.2171, 0.1015; 0.4655, 0.2168, 0.2168, 0.101; 0.4659, 0.2167, 0.2167, 0.1008; 0.466, 0.2166, 0.2166, 0.1007 | 0.4661, 0.2166, 0.2166, 0.1007 | 0.4472, 0.2215, 0.2215, 0.1097 | 0.4661, 0.2166, 0.2166, 0.1007 | 3.93e-09, 1.21e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 3.7e-03, 1.2e-03, 3.7e-04, 1.2e-04 | 8.0e-06 |
| vertex | 1 | 1.000, 3.000, 4.000 | 0.6136, 0.2907, 0.06498, 0.03078; 0.6182, 0.2891, 0.06317, 0.02959; 0.6196, 0.2886, 0.06259, 0.02921; 0.6201, 0.2884, 0.06241, 0.02909; 0.6202, 0.2883, 0.06235, 0.02905 | 0.6203, 0.2883, 0.06233, 0.02904 | 0.5963, 0.2954, 0.07247, 0.0359 | 0.6204, 0.2883, 0.06231, 0.02896 | 4.00e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.6e-03, 1.5e-03, 4.7e-04, 1.5e-04 | 9.9e-06 |
| vertex | 1 | 2.000, 2.000, 4.000 | 0.6674, 0.1496, 0.1496, 0.03349; 0.6734, 0.1472, 0.1472, 0.03222; 0.6753, 0.1464, 0.1464, 0.03182; 0.6759, 0.1462, 0.1462, 0.0317; 0.6761, 0.1461, 0.1461, 0.03166 | 0.6762, 0.1461, 0.1461, 0.03164 | 0.6448, 0.1582, 0.1582, 0.03882 | 0.6763, 0.1461, 0.1461, 0.03156 | 3.97e-09, 1.22e-10, 3.81e-12, 1.20e-13, 3.80e-15 | 6.0e-03, 1.9e-03, 6.1e-04, 1.9e-04 | 1.3e-05 |
| vertex | 1 | 2.000, 5.000, 7.000 | 0.7975, 0.1791, 0.01897, 0.004466; 0.8023, 0.1755, 0.01802, 0.004174; 0.8038, 0.1744, 0.01772, 0.004085; 0.8043, 0.174, 0.01762, 0.004056; 0.8045, 0.1739, 0.01759, 0.004047 | 0.8046, 0.1738, 0.01758, 0.004043 | 0.7797, 0.1913, 0.02326, 0.005707 | 0.8048, 0.1739, 0.01755, 0.003791 | 4.07e-09, 1.23e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.8e-03, 1.5e-03, 4.9e-04, 1.5e-04 | 1.0e-05 |
| vertex | 1 | 4.000, 4.000, 8.000 | 0.9063, 0.04556, 0.04556, 0.002559; 0.9105, 0.04356, 0.04356, 0.002374; 0.9118, 0.04293, 0.04293, 0.002318; 0.9122, 0.04274, 0.04274, 0.0023; 0.9124, 0.04267, 0.04267, 0.002295 | 0.9124, 0.04264, 0.04264, 0.002292 | 0.8896, 0.05357, 0.05357, 0.003225 | 0.9127, 0.04265, 0.04265, 0.001993 | 4.04e-09, 1.22e-10, 3.82e-12, 1.20e-13, 3.80e-15 | 4.2e-03, 1.3e-03, 4.1e-04, 1.3e-04 | 8.8e-06 |

The largest difference of the shares from the inner problem's (at the smallest ε): edge — the loop 5.31e-05; vertex — the loop 9.57e-05; vertex — the edges’ product 9.14e-04; edge — the softmax 9.66e-03; vertex — the softmax 3.14e-02. The loop's successive changes fall by 3.15 per decade of ε (median; √10 = 3.16), their last 1.9e-04 at most; the inner problem's change between its starts at most 1.3e-05. The edges' product against the 4-kernel problem, by amplitude (all shares; the edges' runner-ups only — the rest is the diagonal's frozen tail, which does not factorise, and the entrant's balance): A = 0.3: 9.1e-04 (2.0e-04); A = 1: 3.0e-04 (4.9e-05).

## T. The edge's split along the face

The 2-kernel inner problem (A = 0.8) with the feature on the face at the offset t from the edge's axis: the runner-up's share p(x; t) at x = 2, 4, 6, 8 and the crossing x_θ(t) of 0.01:

| t/H | p(2) | p(4) | p(6) | p(8) | x_θ(t) |
|---|---|---|---|---|---|
| 0 | 0.18740 | 0.05056 | 0.01234 | 0.00319 | 6.300 |
| 0.125 | 0.18679 | 0.05019 | 0.01220 | 0.00315 | 6.284 |
| 0.25 | 0.18498 | 0.04906 | 0.01179 | 0.00302 | 6.233 |
| 0.375 | 0.18195 | 0.04721 | 0.01112 | 0.00282 | 6.149 |
| 0.5 | 0.17771 | 0.04470 | 0.01024 | 0.00256 | 6.034 |

x_θ averaged over the face (t uniform in [−H/2, H/2]): 6.208 (at t = 0: 6.300; at the vertex, t = H/2: 6.034). For a feature uniform in a Voronoi cell the offsets from a face are (2σ²/H)ε^{1/2}x at leading order, so P(N ≥ 2) = 4·(2σ²/H²)·x̄_θ·ε^{1/2} = 71.5·ε^{1/2} (the four faces' layers, less O(ε) where they overlap).

Beyond H/2 (T2, by bisection): x_θ = 5.888 (t = 0.625H), 5.715 (t = 0.75H), 5.517 (t = 0.875H), 5.298 (t = 1H). In the continuum x_θ depends on |t| only, so on hℤⁿ (σ/H as here) a face is the (n − 1)-cube of offsets and E[N − 1] = 2n(2σ²/H²)x̄_θ⁽ⁿ⁾ε^{1/2}: n = 2: x̄_θ = 6.208, 71.5·ε^{1/2}; n = 3: x̄_θ = 6.120, 105.8·ε^{1/2}; n = 4: x̄_θ = 6.033, 139.0·ε^{1/2}; n = 5: x̄_θ = 5.947, 171.3·ε^{1/2} (the n = 2 average by the midpoint rule; the trapezoid rule above).

## B. The formation map at a vertex

A = 0.8; the gaps (x₁₀, x₀₁) of the edge runner-ups on a grid, x₁₀ ≤ x₀₁, the diagonal's x₁₁ = x₁₀ + x₀₁ − ε^{1/2}x₁₀x₀₁ exactly (the tensor grid's correlations multiply); the number of kernels formed by the loop at ε = 1e-06, 1e-08, 1e-10, by the 4-kernel inner problem (B4's), by the independent rule #{j: x_j < x_θ} + 1 with x_θ = 6.034 (the edge's inner problem at the vertex's offset t = H/2 (T)), and by the product rule with the edges' factors at the vertex's offset t = H/2 (T):

| x₁₀ | x₀₁ | x₁₁ | loop | inner problem | the independent rule | the product rule (T) |
|---|---|---|---|---|---|---|
| 0 | 0 | 0.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 1 | 1.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 2 | 2.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 3 | 3.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 4 | 4.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 5 | 5.00 | 4, 4, 4 | 4 | 4 | 4 |
| 0 | 6 | 6.00 | 2, 2, 2 | 2 | 4 | 2 |
| 0 | 7 | 7.00 | 2, 2, 2 | 2 | 2 | 2 |
| 0 | 8 | 8.00 | 2, 2, 2 | 2 | 2 | 2 |
| 0 | 9 | 9.00 | 2, 2, 2 | 2 | 2 | 2 |
| 0 | 10 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 1 | 2.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 2 | 3.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 3 | 4.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 4 | 5.00 | 4, 4, 4 | 4 | 4 | 4 |
| 1 | 5 | 6.00 | 3, 3, 3 | 3 | 4 | 3 |
| 1 | 6 | 7.00 | 2, 2, 2 | 2 | 3 | 2 |
| 1 | 7 | 8.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 8 | 9.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 9 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 1 | 10 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 2 | 4.00 | 4, 4, 4 | 4 | 4 | 4 |
| 2 | 3 | 5.00 | 4, 4, 4 | 4 | 4 | 4 |
| 2 | 4 | 6.00 | 3, 3, 3 | 3 | 4 | 3 |
| 2 | 5 | 7.00 | 3, 3, 3 | 3 | 3 | 3 |
| 2 | 6 | 8.00 | 2, 2, 2 | 2 | 3 | 2 |
| 2 | 7 | 9.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 8 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 9 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 2 | 10 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 3 | 6.00 | 3, 3, 3 | 3 | 4 | 3 |
| 3 | 4 | 7.00 | 3, 3, 3 | 3 | 3 | 3 |
| 3 | 5 | 8.00 | 3, 3, 3 | 3 | 3 | 3 |
| 3 | 6 | 9.00 | 3, 2, 2 | 2 | 3 | 2 |
| 3 | 7 | 10.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 8 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 9 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 3 | 10 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 4 | 8.00 | 3, 3, 3 | 3 | 3 | 3 |
| 4 | 5 | 9.00 | 3, 3, 3 | 3 | 3 | 3 |
| 4 | 6 | 10.00 | 3, 2, 2 | 2 | 3 | 2 |
| 4 | 7 | 11.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 8 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 9 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 4 | 10 | 14.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 5 | 10.00 | 3, 3, 3 | 3 | 3 | 3 |
| 5 | 6 | 11.00 | 3, 3, 3 | 3 | 3 | 3 |
| 5 | 7 | 12.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 8 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 9 | 14.00 | 2, 2, 2 | 2 | 2 | 2 |
| 5 | 10 | 15.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 6 | 12.00 | 3, 3, 3 | 3 | 3 | 3 |
| 6 | 7 | 13.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 8 | 14.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 9 | 15.00 | 2, 2, 2 | 2 | 2 | 2 |
| 6 | 10 | 16.00 | 2, 2, 2 | 2 | 2 | 2 |
| 7 | 7 | 14.00 | 1, 1, 1 | 1 | 1 | 1 |
| 7 | 8 | 15.00 | 1, 1, 1 | 1 | 1 | 1 |
| 7 | 9 | 16.00 | 1, 1, 1 | 1 | 1 | 1 |
| 7 | 10 | 17.00 | 1, 1, 1 | 1 | 1 | 1 |
| 8 | 8 | 16.00 | 1, 1, 1 | 1 | 1 | 1 |
| 8 | 9 | 17.00 | 1, 1, 1 | 1 | 1 | 1 |
| 8 | 10 | 18.00 | 1, 1, 1 | 1 | 1 | 1 |
| 9 | 9 | 18.00 | 1, 1, 1 | 1 | 1 | 1 |
| 9 | 10 | 19.00 | 1, 1, 1 | 1 | 1 | 1 |
| 10 | 10 | 20.00 | 1, 1, 1 | 1 | 1 | 1 |

The loop's count equals the inner problem's at 64 of 66 (ε = 1e-06), 66 of 66 (ε = 1e-08), 66 of 66 (ε = 1e-10); the independent rule at 58 of 66, the product rule at 66 of 66.

## B4. The areas at a vertex

The 4-kernel inner problem on the grid of gaps 0 … 10 in steps of 0.25 (861 problems, x₁₀ ≤ x₀₁): the region N ≥ 3 has the area 33.5 and N = 4 15.5 in the gaps' quarter-plane (trapezoid weights; 648 of 861 problems settled to the shares' tolerance 1e-07; the 213 others rerun with the stop's reason recorded: 209 at the doubling cap (the last drift 1.0e-07 … 3.7e-07; the change still to come, from the drifts' fall, at most 1.6e-05 on 101; the drifts not yet falling on 108), 1 at the clock (the last drift 2.0e-07; the drifts not yet falling), 3 settled on rerun, their counts unchanged at 213 of 213, their shares moved by at most 1.1e-06; the unsettled shares at least 2.0e-04 from the threshold, so at their last drifts no count changes for at least 788 more doublings); the product rule with the factors at t = H/2 on the same grid: 33.5 and 15.1 (on a grid of step 0.05: 33.5 and 15.3, the coarse grid's error), its count equal to the inner problem's at 1675 of 1681 grid points. Each vertex's sector in a cell maps to the gaps by (2σ²ε^{1/2}/H)·(x₁₀, x₀₁) and a cell has four: P(N ≥ 3) = 4(2σ²/H²)²·33.5·ε = 1112·ε and P(N = 4) = 513·ε.

## C. Several features in the plane

Random three-feature targets (seed 54; centres in [−0.7, 0.7]² at least 3σ apart, amplitudes ±U(0.4, 1)) on the lattice over [−1, 1]² (17² kernels; the samples' grid 41²); the sliding path against the loop to 1.3 times the path's last entry; the never-entering kernels' displacement vectors against the creep's εC_j integrated along the path (`creep_constants_nd`; where |εC_j| is above 10⁻³ of the largest: the largest relative error and the smallest cosine):

| target | the path (entries) | ε | formed = the path's | the never-entering kernels' largest \|w\|·ε^{−3/2} | their largest displacement/ε | creep: relative error, smallest cosine | the block: max \|w − w_path\|, max \|c − c_path\| | rms residual |
|---|---|---|---|---|---|---|---|---|
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-06 | yes | 2862 | 4838 | 4.57e-01, 0.979100 | 3.01e-06, 2.14e-07 | 4.61e-08 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-07 | yes | 2039 | 3849 | 1.21e-01, 0.998353 | 6.94e-08, 4.99e-09 | 1.11e-09 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-08 | yes | 1851 | 3594 | 3.64e-02, 0.999843 | 2.01e-09, 1.43e-10 | 3.23e-11 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-09 | yes | 1797 | 3518 | 1.13e-02, 0.999985 | 6.18e-11, 4.40e-12 | 9.94e-13 |
| 0 | [146, 174, 234] (68.88, 69.53, 95.23) | 1e-10 | yes | 1780 | 3495 | 3.56e-03, 0.999998 | 1.94e-12, 1.38e-13 | 3.12e-14 |
| 0 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.508 | 1.006 | | 1.508, 1.509 | |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-06 | yes | 1.369e+05 | 4.749e+04 | 6.02e-01, 1.000000 | 1.36e-04, 3.63e-06 | 1.27e-07 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-07 | yes | 4.676e+04 | 3.729e+04 | 2.58e-01, 1.000000 | 1.41e-06, 1.16e-07 | 1.23e-08 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-08 | yes | 3.576e+04 | 3.193e+04 | 7.71e-02, 1.000000 | 3.41e-08, 2.94e-09 | 3.29e-10 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-09 | yes | 3.305e+04 | 3.035e+04 | 2.38e-02, 1.000000 | 9.95e-10, 8.63e-11 | 9.70e-12 |
| 1 | [57, 166, 226] (79.32, 91.30, 147.80) | 1e-10 | yes | 3.225e+04 | 2.987e+04 | 7.48e-03, 1.000000 | 3.07e-11, 2.66e-12 | 3.00e-13 |
| 1 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.522 | 1.014 | | 1.522, 1.522 | |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-06 | yes | 6041 | 5.122e+04 | 1.67e+00, 0.999876 | 4.89e-06, 9.40e-07 | 2.04e-07 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-07 | yes | 3088 | 2.318e+04 | 2.06e-01, 0.999988 | 9.18e-08, 1.17e-08 | 1.76e-09 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-08 | yes | 2736 | 2.027e+04 | 5.49e-02, 0.999999 | 2.60e-09, 3.22e-10 | 4.81e-11 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-09 | yes | 2642 | 1.953e+04 | 1.65e-02, 1.000000 | 7.94e-11, 9.82e-12 | 1.47e-12 |
| 2 | [91, 113, 214] (87.45, 94.62, 120.73) | 1e-10 | yes | 2613 | 1.932e+04 | 5.14e-03, 1.000000 | 2.48e-12, 3.07e-13 | 4.60e-14 |
| 2 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.510 | 1.010 | | 1.510, 1.510 | |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-06 | yes | 6.347e+04 | 3.56e+04 | 5.70e-01, 0.999984 | 6.27e-05, 3.88e-06 | 1.65e-07 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-07 | yes | 2.763e+04 | 2.731e+04 | 2.04e-01, 0.999999 | 8.32e-07, 9.93e-08 | 8.09e-09 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-08 | yes | 2.232e+04 | 2.408e+04 | 6.18e-02, 1.000000 | 2.12e-08, 2.62e-09 | 2.22e-10 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-09 | yes | 2.095e+04 | 2.312e+04 | 1.92e-02, 1.000000 | 6.29e-10, 7.79e-11 | 6.63e-12 |
| 3 | [58, 97, 224] (71.06, 74.86, 109.32) | 1e-10 | yes | 2.054e+04 | 2.282e+04 | 6.04e-03, 1.000000 | 1.95e-11, 2.42e-12 | 2.06e-13 |
| 3 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.518 | 1.012 | | 1.518, 1.517 | |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-06 | yes | 3.851e+06 | 8.48e+04 | 4.98e-01, 0.999098 | 3.87e-03, 5.69e-06 | 1.91e-07 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-07 | yes | 6.017e+05 | 2.339e+05 | 3.87e-01, 0.999963 | 1.87e-05, 1.51e-06 | 1.49e-07 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-08 | yes | 3.632e+05 | 1.894e+05 | 1.24e-01, 0.999997 | 3.49e-07, 3.69e-08 | 5.02e-09 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-09 | yes | 3.187e+05 | 1.75e+05 | 3.80e-02, 1.000000 | 9.68e-09, 1.04e-09 | 1.45e-10 |
| 4 | [98, 173, 198] (66.87, 111.72, 138.55) | 1e-10 | yes | 3.066e+05 | 1.706e+05 | 1.19e-02, 1.000000 | 2.95e-10, 3.17e-11 | 4.43e-12 |
| 4 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.537 | 1.023 | | 1.537, 1.533 | |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-06 | yes | 8.397e+04 | 4.497e+04 | 5.28e-01, 0.999917 | 8.38e-05, 3.77e-06 | 2.17e-07 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-07 | yes | 3.498e+04 | 3.532e+04 | 2.00e-01, 0.999992 | 1.06e-06, 1.03e-07 | 1.28e-08 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-08 | yes | 2.82e+04 | 3.121e+04 | 6.07e-02, 0.999999 | 2.70e-08, 2.73e-09 | 3.54e-10 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-09 | yes | 2.648e+04 | 2.998e+04 | 1.89e-02, 1.000000 | 8.03e-10, 8.13e-11 | 1.06e-11 |
| 5 | [115, 125, 192] (66.63, 97.11, 121.78) | 1e-10 | yes | 2.597e+04 | 2.96e+04 | 5.93e-03, 1.000000 | 2.49e-11, 2.52e-12 | 3.29e-13 |
| 5 | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | 1.518 | 1.012 | | 1.518, 1.517 | |

## D. The probabilities over uniformly placed features

A = 0.8 at 600 uniform positions in the Voronoi cell of the kernel at (0, 0) (seed 541; the same positions at both masses, so the two rows are not independent); the loop's count N against the product rule at each position's gaps to its two nearest faces and its offsets along them (T), and the fractions against the rule's at the sampled positions and the asymptotics (P(N ≥ 2) = 4w/H − 4w²/H², w = 2σ²x̄_θε^{1/2}/H with the face average x̄_θ = 6.208; P(N ≥ 3) from B4's area):

| ε | N = 1, 2, 3, 4 (loop) | the product rule agrees | P(N ≥ 2): loop (± s.e.) | the rule at the positions | asymptotic | P(N ≥ 3): loop | the rule at the positions | asymptotic |
|---|---|---|---|---|---|---|---|---|
| 1e-05 | 458, 133, 2, 7 | 592 of 600 | 0.2367 ± 0.0174 | 0.2267 | 0.2134 | 0.0150 | 0.0133 | 0.0111 |
| 1e-06 | 548, 49, 1, 2 | 597 of 600 | 0.0867 ± 0.0115 | 0.0833 | 0.0702 | 0.0050 | 0.0033 | 0.0011 |

The misses — ε: the loop's count/the rule's at the gaps (x₁₀, x₀₁) [x_θ at their offsets]: 1e-05: 2/1 at (26.49, 6.28) [6.09, 6.23]; 1e-05: 2/1 at (6.52, 13.35) [6.15, 6.09]; 1e-05: 4/3 at (0.59, 4.85) [6.08, 6.04]; 1e-05: 2/1 at (6.36, 8.41) [6.11, 6.09]; 1e-05: 2/1 at (42.90, 6.49) [6.09, 6.29]; 1e-05: 3/2 at (6.12, 3.45) [6.06, 6.09]; 1e-05: 2/1 at (6.33, 7.32) [6.10, 6.09]; 1e-05: 2/1 at (6.30, 46.00) [6.29, 6.09]; 1e-06: 2/1 at (6.11, 8.25) [6.06, 6.05]; 1e-06: 3/2 at (3.92, 6.03) [6.05, 6.04]; 1e-06: 2/1 at (6.29, 76.61) [6.22, 6.05].

The offset-blind rule (every factor at t = 0, x_θ(0)) agrees with the loop at 596 of 600 (ε = 1e-05), 600 of 600 (ε = 1e-06) — D does not resolve x_θ(t) (referee round 47): the positions near a face, where the offset matters, are too few.


## H. A hexagonal lattice

Kernels on a hexagonal lattice of spacing 0.125, the feature (A = 0.8) at a Voronoi vertex — three cells meet there, the kernels (0, 0), (H, 0), (H/2, √3H/2) at H/√3 —, placed where the three correlate equally with it on the samples. The 3-kernel inner problem at the exact tie: shares 0.33333, 0.33333, 0.33333; on the gaps' grid 0 … 12 (step 0.5) the region where all three form, N = 3, has the area 36.6 (the square vertex's N ≥ 3: B4). A sector of a hexagonal vertex maps to the gaps with the Jacobian (2σ²ε^{1/2})²/((√3/2)H²) and a cell has six: P(N = 3) = 32σ⁴/H⁴·36.6·ε = 2426·ε. The loop at the exact tie:

| ε | shares | formed | the others' largest \|w\| | the gaps of the runner-ups and of the fourth kernel (units of ε^{1/2}) |
|---|---|---|---|---|
| 1e-08 | 0.33333, 0.33333, 0.33333 | 3 | 2.39e-11 | 5.551e-12, 5.551e-12, 1594 |
| 1e-10 | 0.33333, 0.33333, 0.33333 | 3 | 2.38e-14 | 5.551e-11, 5.551e-11, 1.594e+04 |

The fourth kernel's margin at the vertex is ε-free, 1 − e^{−H²/(4σ²)} in the continuum, so its gap grows like ε^{−1/2}: a fourth kernel never forms at a simple vertex (k ≤ n + 1).

## C2. Designed targets in the plane

Part C's lattice; 'near tie': feature 0 (A = 0.8) near the face between (0, 0.25) and (0.125, 0.25), on their axis, its runner-up's margin 0.0003 when it enters (first), two features far from it (−0.6 at (−0.5, −0.5), 0.7 at (0.5, −0.55)); 'window edge': feature 0 (0.8) at (0.96, 0) with −0.7 at (−0.4, 0.5) and 0.9 at (−0.3, −0.5). The loop to 1.3 times the path's last entry; feature 0's entrant and runner-up (the target's correlations at the start) and its nearest kernel; the runner-up's share in the loop (of the kernels started nearest feature 0) against the edge's 2-kernel inner problem p(x; t) at its gap (T), and its coefficient at the first output after the entry against the later features' unbinding level max A′(d/σ)e^{−d²/(4σ²)} (A′, d their amplitudes and distances from feature 0); the kernels formed for feature 0 (|w| > 1% of its amplitude; by their final positions) against the rule 1 + [x < x_θ(t)], and the formed kernels that ended at another feature than the one they started nearest (captured):

| target | ε | entrant, runner-up, nearest | margin μ, gap x = μ/ε^{1/2} | runner-up's share: loop, p(x; t) | the runner-up after the entry: \|w\| (τ), unbinding level | formed for feature 0: loop, rule | captured (kernel: start's feature → end's) | formed = the path's | max\|w − w_path\| | rms residual |
|---|---|---|---|---|---|---|---|---|---|---|
| near tie | 1e-06 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 0.3 | 0.4464, 0.4452 | 0.357 (79.3), 0.00041 | 2, 2 | none | no: [72, 146, 163, 208] | 3.57e-01 | 9.87e-08 |
| near tie | 1e-07 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 0.949 | 0.3338, 0.3328 | 0.267 (79.3), 0.00041 | 2, 2 | none | no: [72, 146, 163, 208] | 2.67e-01 | 3.47e-08 |
| near tie | 1e-08 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 3 | 0.1002, 0.09973 | 0.0802 (79.3), 0.00041 | 2, 2 | none | no: [72, 146, 163, 208] | 8.01e-02 | 4.02e-08 |
| near tie | 1e-09 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 9.49 | — (captured), 0.001322 | 0.00106 (79.3), 0.00041 | 1, 1 | (0.125, 0.250): 0 → 2 | no: [72, 146, 163] | 7.00e-01 | 1.95e-13 |
| near tie | 1e-10 | (0.000, 0.250), (0.125, 0.250), (0.000, 0.250) | 3.000e-04, 30 | 8.205e-06, — | 6.48e-06 (79.3), 0.00041 | 1, 1 | none | yes | 6.33e-06 | 4.28e-08 |
| window edge | 1e-06 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 159 | 3.295e-08, — | 1.89e-09 (71.3), 1.1e-08 | 1, 1 | none | yes | 4.15e-06 | 3.94e-08 |
| window edge | 1e-07 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 504 | 9.693e-10, — | 5.98e-11 (71.3), 1.1e-08 | 1, 1 | none | yes | 8.72e-08 | 1.03e-09 |
| window edge | 1e-08 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 1.59e+03 | 2.998e-11, — | 1.89e-12 (71.3), 1.1e-08 | 1, 1 | none | yes | 2.46e-09 | 2.95e-11 |
| window edge | 1e-09 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 5.04e+03 | 9.416e-13, — | 5.97e-14 (71.3), 1.1e-08 | 1, 1 | none | yes | 7.53e-11 | 9.04e-13 |
| window edge | 1e-10 | (0.875, 0.000), (0.875, -0.125), (1.000, 0.000) | 1.594e-01, 1.59e+04 | 2.971e-14, — | 1.89e-15 (71.3), 1.1e-08 | 1, 1 | none | yes | 2.35e-12 | 2.83e-14 |

## V. A runner-up held at its feature: the static law of (v)

The feature (A = 0.8) at the origin and an unfitted one (A′ = 0.7) at the distance d (the continuum: a trapezoid of step σ/20 in one dimension, σ/8 in the plane along the diagonal); the entrant (w_e, c_e) and the runner-up's centre free, its coefficient w fixed: the local minimum's displacement Δ of the runner-up toward the unfitted feature (by continuation in w) against (v)'s law Δ² = σ²A′((d/σ)² − 2)e^{−d²/(4σ²)}·w_e/(w(w_e + w)); the fold — below it no minimum within 2σ — by bisection; the level ℓ(d) = A′(d/σ)²e^{−d²/(4σ²)}, the fold over it c(d), and the law's Δ at the fold:

| d/σ | Δ/σ (the law) at w = 0.3, 0.1, 0.03, 0.01, 0.003, 0.001 | the fold w_f | ℓ(d) | c(d) = w_f/ℓ | the law's Δ/σ at w_f |
|---|---|---|---|---|---|
| 5 | 0.3164 (0.2556), —, —, —, —, — | 0.1424 | 0.0338 | 4.21 | 0.424 |
| 5.5 | 0.1641 (0.1465), 0.4503 (0.3120), —, —, —, — | 0.05956 | 0.011 | 5.41 | 0.400 |
| 6 | 0.0837 (0.0783), 0.1983 (0.1667), 0.4823 (0.3122), —, —, — | 0.02066 | 0.00311 | 6.64 | 0.372 |
| 6.5 | 0.0407 (0.0390), 0.0912 (0.0831), 0.1862 (0.1555), 0.3708 (0.2596), —, — | 0.006018 | 0.000765 | 7.87 | 0.347 |
| 7 | 0.0186 (0.0181), 0.0407 (0.0386), 0.0791 (0.0723), 0.1397 (0.1206), 0.3325 (0.2363), — | 0.001495 | 0.000164 | 9.11 | 0.324 |

The plane, d = 6σ along the diagonal: w = 0.3: along 0.0837, across 2.3e-08 (the law 0.0783); w = 0.079: along 0.2224, across 2.2e-07 (the law 0.1831); w = 0.0208: along 0.8015, across 7.2e-08 (the law 0.3709); the fold between w = 0.02066 and 0.02066.

## C3. The capture window

C2's near tie (feature 2, 0.7, at 6.08σ from feature 0) and its variant with feature 2 at 5.00σ (−60°), at 15 masses, the loop at the tolerance 1e-12 (captures sit on a saddle): the runner-up's gap x, its coefficient after the entry against the fold at feature 2's distance (V: c(d)ℓ(d), interpolated in d), the kernels formed per feature by their final positions, the captured kernels, and the outcome — bound (feature 0 two kernels, the runner-up formed there), captured (feature 0 one kernel, feature 2 fitted by the migrant), a capture tie (feature 2 two kernels), stranded (one kernel per feature, the fit not exact at the end: the tail still on its way), one kernel (the runner-up not formed, the fit exact):

| target | ε | gap x | the runner-up after the entry | the fold c(d)ℓ(d) | kernels per feature (final) | captured | rms | outcome |
|---|---|---|---|---|---|---|---|---|
| near tie | 1e-08 | 3 | 0.0801 | 0.0172 (d = 6.08σ) | 2, 1, 1 | none | 4.0e-08 | bound |
| near tie | 7.05e-09 | 3.57 | 0.0545 | 0.0172 (d = 6.08σ) | 2, 1, 1 | none | 5.1e-08 | bound |
| near tie | 4.97e-09 | 4.26 | 0.034 | 0.0172 (d = 6.08σ) | 2, 1, 1 | none | 7.1e-08 | bound |
| near tie | 3.5e-09 | 5.07 | 0.0191 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 5.1e-10 | captured |
| near tie | 3e-09 | 5.48 | 0.0143 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 8.6e-13 | captured |
| near tie | 2.47e-09 | 6.04 | 0.00966 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 7.3e-13 | captured |
| near tie | 1.74e-09 | 7.19 | 0.00435 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 4.2e-13 | captured |
| near tie | 1.22e-09 | 8.59 | 0.00178 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.5e-13 | captured |
| near tie | 8.63e-10 | 10.2 | 0.000722 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.6e-13 | captured |
| near tie | 6.08e-10 | 12.2 | 0.000296 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.3e-13 | captured |
| near tie | 4.29e-10 | 14.5 | 0.000128 | 0.0172 (d = 6.08σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.4e-13 | captured |
| near tie | 3.13e-10 | 17 | 6.26e-05 | 0.0172 (d = 6.08σ) | 1, 1, 2 | (0.125, 0.250): 0 → 2 | 1.4e-11 | capture tie |
| near tie | 3.02e-10 | 17.3 | 5.79e-05 | 0.0172 (d = 6.08σ) | 1, 1, 1 | none | 2.9e-06 | stranded |
| near tie | 2.13e-10 | 20.6 | 2.78e-05 | 0.0172 (d = 6.08σ) | 1, 1, 1 | none | 2.8e-07 | stranded |
| near tie | 1.5e-10 | 24.5 | 1.39e-05 | 0.0172 (d = 6.08σ) | 1, 1, 1 | none | 9.8e-08 | stranded |
| near tie 5σ | 1e-08 | 3 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 4.9e-12 | captured |
| near tie 5σ | 7.05e-09 | 3.57 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.9e-12 | captured |
| near tie 5σ | 4.97e-09 | 4.26 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.7e-12 | captured |
| near tie 5σ | 3.5e-09 | 5.07 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.0e-12 | captured |
| near tie 5σ | 3e-09 | 5.48 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 8.0e-13 | captured |
| near tie 5σ | 2.47e-09 | 6.04 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 6.0e-13 | captured |
| near tie 5σ | 1.74e-09 | 7.19 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 3.5e-13 | captured |
| near tie 5σ | 1.22e-09 | 8.59 | 0.7 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.1e-13 | captured |
| near tie 5σ | 8.63e-10 | 10.2 | 0.339 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.2e-13 | captured |
| near tie 5σ | 6.08e-10 | 12.2 | 0.000303 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 7.3e-14 | captured |
| near tie 5σ | 4.29e-10 | 14.5 | 0.00013 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 4.3e-14 | captured |
| near tie 5σ | 3.13e-10 | 17 | 6.35e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.7e-14 | captured |
| near tie 5σ | 3.02e-10 | 17.3 | 5.87e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 2.6e-14 | captured |
| near tie 5σ | 2.13e-10 | 20.6 | 2.82e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 1.5e-14 | captured |
| near tie 5σ | 1.5e-10 | 24.5 | 1.41e-05 | 0.142 (d = 5.00σ) | 1, 1, 1 | (0.125, 0.250): 0 → 2 | 9.3e-15 | captured |

## T3. The product rule at a vertex of hℤ³

The 8 kernels at the corners of the cube {0, H}³ (samples on a tensor grid of spacing 0.05), the feature (A = 0.8) at the cube's centre — each edge's axis at the transverse offset H/√2 ((n − 1)^{1/2}H/2; H/2 only in the plane) —; the corner 1_S behind the entrant by Σ_{i∈S}x_i: the 8-kernel inner problem against the product of the edges' 2-kernel shares at the offset H/√2 (the planar edge problem: in the continuum it depends on |t| only) and at H/2 (the planar offset):

| gaps (x₁, x₂, x₃) | 8-kernel: formed (stop) | the product at H/√2: formed, max \|Δ share\| | at H/2: formed, max \|Δ share\| |
|---|---|---|---|
| (1, 1, 1) | 8 (converged) | 8, 8.68e-05 | 8, 1.03e-02 |
| (2, 2, 2) | 7 (converged) | 7, 3.58e-04 | 7, 1.95e-02 |
| (1, 2, 3) | 6 (converged) | 6, 3.54e-04 | 7, 1.58e-02 |
| (0.5, 3, 5) | 5 (converged) | 4, 5.56e-04 | 5, 9.59e-03 |
| (3, 3, 3) | 4 (converged) | 4, 9.28e-04 | 4, 1.93e-02 |
| (2, 4, 6) | 3 (cap) | 3, 8.79e-04 | 3, 1.42e-02 |
| (4, 4, 4) | 4 (converged) | 4, 1.16e-03 | 4, 1.38e-02 |

## W. Where the selection is the Voronoi cell

W1. The samples' window (Part C's grid: the samples on [−1, 1]², the lattice up to ±1). The faces between neighbours on the axis x₂ = 0 (where their correlations with a feature are equal on the samples; the continuum's: the midpoints):

| pair | the midpoint | the samples' face |
|---|---|---|
| 0 \| 0.125 | 0.0625 | 0.0625 |
| 0.125 \| 0.25 | 0.1875 | 0.1875 |
| 0.25 \| 0.375 | 0.3125 | 0.3125 |
| 0.375 \| 0.5 | 0.4375 | 0.4375 |
| 0.5 \| 0.625 | 0.5625 | 0.5625 |
| 0.625 \| 0.75 | 0.6875 | 0.6881 |
| 0.75 \| 0.875 | 0.8125 | 0.8260 |
| 0.875 \| 1 | 0.9375 | 1.1028 |

The largest error of the nearest kernel's neighbours' margins against the continuum's, features on the axis at x₁: 0: 3.3e-16, 0.3: 3.0e-11, 0.5: 4.3e-06, 0.55: 1.4e-05, 0.6: 3.1e-04, 0.65: 7.9e-04, 0.7: 7.2e-03, 0.75: 1.4e-02, 0.8: 2.8e-02, 0.85: 9.2e-02, 0.9: 1.4e-01, 0.95: 2.5e-01; it first exceeds ε^{1/2}x̄_θ (6.208) at x₁ = 0.64 for ε = 1e-08 (2.40σ from the window's edge), 0.57 for ε = 1e-10 (2.87σ from the window's edge); the samples' entrant is not the nearest kernel at 6 of 99 positions, from x₁ = 0.82.

W2. Unequal masses: the kernel at (H, 0) of mass βε, the rest ε (the local lattice, ε = 1e-08); the entrant is the kernel of the least entry time ρ_j/v_j (ρ = b^{−1/2}): a power diagram, the face between (0, 0) and (H, 0) at H/2 − σ²ln β/H in the continuum. Features (A = 0.8) on the axis at the Voronoi face H/2, at the power face p (the samples') and 2 layer widths (w = 2σ²x̄_θε^{1/2}/H = 2.23e-04) to either side:

| β | p: samples, continuum | x₁ | entrant: power, Voronoi | runner-up, its gap | formed | shares of (0, 0), (H, 0) |
|---|---|---|---|---|---|---|
| 1 | 0.06250, 0.06250 | 0.06205 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9997, 0.0003 |
| 1 | 0.06250, 0.06250 | 0.06250 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 1.11e-12 | (+0, +0)H, (+1, +0)H | 0.5000, 0.5000 |
| 1 | 0.06250, 0.06250 | 0.06295 | (+1, +0)H, (+1, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0003, 0.9997 |
| 1.001 | 0.06232, 0.06232 | 0.06187 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9997, 0.0003 |
| 1.001 | 0.06232, 0.06232 | 0.06232 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 2.22e-12 | (+0, +0)H, (+1, +0)H | 0.4978, 0.5022 |
| 1.001 | 0.06232, 0.06232 | 0.06250 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 5 | (+0, +0)H, (+1, +0)H | 0.0251, 0.9749 |
| 1.001 | 0.06232, 0.06232 | 0.06277 | (+1, +0)H, (+1, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0003, 0.9997 |
| 1.01 | 0.06071, 0.06071 | 0.06026 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9996, 0.0004 |
| 1.01 | 0.06071, 0.06071 | 0.06071 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 0 | (+0, +0)H, (+1, +0)H | 0.4782, 0.5218 |
| 1.01 | 0.06071, 0.06071 | 0.06116 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0003, 0.9997 |
| 1.01 | 0.06071, 0.06071 | 0.06250 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 49.6 | (+1, +0)H | 0.0000, 1.0000 |
| 1.5 | -0.01048, -0.01048 | -0.01093 | (+0, +0)H, (+0, +0)H | (+1, +0)H, 12.4 | (+0, +0)H | 0.9919, 0.0081 |
| 1.5 | -0.01048, -0.01048 | -0.01048 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 5.55e-12 | (+0, +0)H, (+1, +0)H | 0.0205, 0.9795 |
| 1.5 | -0.01048, -0.01048 | -0.01004 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 12.4 | (+1, +0)H | 0.0001, 0.9999 |
| 1.5 | -0.01048, -0.01048 | 0.06250 | (+1, +0)H, (+0, +0)H | (+0, +0)H, 1.84e+03 | (+1, +0)H | 0.0000, 1.0000 |

W3. Coarse samples along x₂ (spacing ds, the grid offset by 0.3ds; along x₁ the spacing 0.05): the face between (0, 0) and (0, H) for features at (0, x₂), its shift from H/2 on the samples and from the Poisson sum (harmonics k ≤ 3 of q = e^{−π²σ²/ds²}), the margin's slope there against the continuum's H/(2σ²), and the shift in layer widths at ε = 1e-08:

| ds | ds/σ | q | shift/H: samples, Poisson | slope ratio | shift/w |
|---|---|---|---|---|---|
| 0.05 | 0.333 | 2.650e-39 | -7.772e-16, +0.000e+00 | 1.0000 | 4.35e-13 |
| 0.1 | 0.667 | 2.269e-10 | +2.151e-09, +2.151e-09 | 1.0000 | 1.2e-06 |
| 0.15 | 1 | 5.172e-05 | +3.855e-04, +3.855e-04 | 0.9989 | 0.216 |
| 0.2 | 1.33 | 3.881e-03 | +3.131e-03, +3.131e-03 | 0.9276 | 1.75 |
| 0.3 | 2 | 8.480e-02 | -5.560e-01, -5.560e-01 | 0.7092 | 311 |

W4. Two kernels coinciding at the samples' centre of symmetry, the runner-up x_r behind, the feature (A = 0.8) δ off their centre: the inner problem's runner-up share (extrapolated in xi0) against the profile split's (iii)(b) — off the centre the kernels slide onto the feature at speeds ∝ their coefficients, apart on the way —, and the kernels' largest distance from the feature at the end:

| x_r | δ/σ | inner problem | profile split | difference | at the end (σ) |
|---|---|---|---|---|---|
| 2 | 0 | 0.197048 | 0.197048 | -1.53e-08 | 8.5e-13 |
| 2 | 0.07 | 0.196776 | 0.197048 | -2.72e-04 | 3.5e-13 |
| 2 | 0.2 | 0.194823 | 0.197048 | -2.22e-03 | 4.0e-13 |
| 2 | 0.4 | 0.188132 | 0.197048 | -8.92e-03 | 1.7e-12 |
| 4 | 0 | 0.056862 | 0.056862 | -1.06e-08 | 8.9e-13 |
| 4 | 0.07 | 0.056678 | 0.056862 | -1.84e-04 | 1.8e-12 |
| 4 | 0.2 | 0.055370 | 0.056862 | -1.49e-03 | 3.3e-12 |
| 4 | 0.4 | 0.051028 | 0.056862 | -5.83e-03 | 1.6e-12 |

W5. The 2-fold profile split's runner-up share p(x) against the softmax 1/(1 + e^{κx}) at the formation threshold: p's crossing of 0.01 (x_θ of coinciding kernels at the feature's centre) and the softmax's, the softmax's relative error at p's crossing and at x = 8, and e^{−2α/a} at the runner-up's coefficient (its mirror coordinate α; the exponential branch's correction, D‴.51(i)):

| A | x_θ: profile split, softmax | the softmax's error: at the crossing, at x = 8 | α/a | e^{−2α/a} |
|---|---|---|---|---|
| 0.3 | 6.688, 6.541 | -0.097, -0.221 | 0.886 | 0.170 |
| 0.8 | 6.579, 6.541 | -0.026, -0.091 | 1.204 | 0.090 |
| 1 | 6.568, 6.541 | -0.018, -0.071 | 1.287 | 0.076 |

Numbers are produced by `experiments/exp_transit_nd.py`; nothing is edited by hand.
