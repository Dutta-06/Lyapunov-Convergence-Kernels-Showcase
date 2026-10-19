# Numerical results for PARADIGM.md (transport adaptive control)

Every number here is produced by the named script from raw files, and no table is edited by hand.

- `tests/test_transport.py`, `tests/test_langevin.py`, `tests/test_averaging.py`, `tests/test_excitation.py`,
  `tests/test_nonreciprocity.py`, `tests/test_periodic.py` and `tests/test_barycentric.py` check the inequalities at every
  sample of simulated trajectories, or on sampled states.
- Run them with `venv\Scripts\python -m pytest tests`. The current count is in `README.md`.

## Theorem B — label-free Lyapunov function (`exp_transport_permutation.py` → `transport_permutation/`)

Each network starts as a random relabelling of the ideal configuration plus noise 0.05, with mean-field gains and 3
seeds per m. "Excess" is |(1/m)Σ Taylor remainders|, the part that actually enters V̇.

| m | labelled error (t=0) | OT error W (t=0) | RMS ‖e‖ | actual excess | bound, matched error | same bound, labelled error |
|---|---|---|---|---|---|---|
| 8 | 2.675 | 0.107 | 0.011 | 2.0e-05 | 8.3e-02 | 2.9 |
| 16 | 2.532 | 0.104 | 0.007 | 9.6e-06 | 4.4e-02 | 1.9 |
| 32 | 2.684 | 0.112 | 0.008 | 2.7e-05 | 7.1e-02 | 2.1 |
| 64 | 2.764 | 0.114 | 0.009 | 1.2e-05 | 6.8e-02 | 2.4 |

- The bound is at least the actual excess at every sample of every run.
- The matched error is about 1/25 of the labelled one.
- The bound is 30–45× smaller when evaluated with the matched error.
- Nothing grows with m.
- The bounds remain far above the actual excess, because the Hessian supremum is global.

**Corollary B.5 along (L_σ)** (test `test_corollary_B5_fails_along_lsigma_and_its_leak_variant_holds`). The test
asserts every hypothesis of C(iv) and measures sup W ≤ η̄. There, the plain B.5 bound is exceeded by more than 3×,
while the leak variant with σ_m b²/4 holds.

## Theorem C — OT σ-modification (`exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py` → `ot_sigma_mod/`)

The two starts share plant, prior, geometry and x₀. They differ only in particle labels.

| start | mechanism | RMS ‖e‖ (20–40 s) | W_end | labelled error at end |
|---|---|---|---|---|
| aligned | projection | 0.0051 ± 0.0010 | 0.295 | 0.295 |
| aligned | labelled σ-mod | 0.0091 ± 0.0018 | 0.304 | 0.304 |
| aligned | OT σ-mod | 0.0091 ± 0.0018 | 0.304 | 0.304 |
| relabelled | projection | 0.0051 ± 0.0010 | 0.295 | 0.295 |
| relabelled | labelled σ-mod | 0.0091 ± 0.0018 | 0.304 | **2.059** |
| relabelled | OT σ-mod | 0.0091 ± 0.0018 | 0.304 | 0.304 |

Reshuffling transient from the relabelled start (5 seeds):

| leak | RMS ‖e‖ [0,5] s | [5,10] s | [10,20] s | [20,40] s | particle path length |
|---|---|---|---|---|---|
| labelled | 0.0771 | 0.0086 | 0.0093 | 0.0091 | 2.26 |
| OT | 0.0743 | 0.0072 | 0.0093 | 0.0091 | 0.39 |

- The labelled leak's longer path costs only about 4% of early RMS‖e‖. Theorem F does not price leak-driven transport
  (PARADIGM.md Proposition C.5).
- The C(iv) instance in the tests satisfies every hypothesis, asserted. V_OT stays below the envelope
  e^{−αt}V(0) + β/α at every sample and contracts by more than 5×.
- The C(ii) test asserts C ≤ b² + 4ηR, √C ≤ b + 2η, and C = b² under both gap conditions, on 300 configurations.

## Theorem D — convergence in m (`exp_mean_field_m.py` → `mean_field_m/`)

| m | sup_t ‖e^m − e^{2m}‖ | W(μ^m_0, μ^{2m}_0) | W(μ^m_T, μ^{2m}_T) | ratio |
|---|---|---|---|---|
| 16 vs 32 | 4.6e-3 | 0.218 | 0.366 | 1.68 |
| 32 vs 64 | 9e-4 | 0.109 | 0.179 | 1.64 |
| 64 vs 128 | 1e-4 | 0.055 | 0.092 | 1.69 |
| 128 vs 256 | <1e-4 | 0.027 | 0.046 | 1.68 |

This is a consistency check in m, not a check of the constant.

- The guaranteed factor e^{ΛT} is at least e^{496} on this box.
- The Van der Pol plant fails the tube condition of Theorem D(b), so the run is covered by the localised D(d).

## Theorem E — label-free convergence under excitation (`exp_pe_convergence.py` → `pe_convergence/`)

| reference | seed | λ_min, first 60 s window | min over sliding 60 s windows | W(T)/W(0) | fitted rate of W, t ≥ 100 s (1/s) | labelled(T)/labelled(0) |
|---|---|---|---|---|---|---|
| rich (6 tones) | 0 | 8.6e-3 | 4.1e-3 | 0.544 | 2.3e-3 | 0.992 |
| rich (6 tones) | 1 | 5.5e-3 | 1.0e-3 | 0.409 | 2.3e-3 | 0.962 |
| rich (6 tones) | 2 | 1.1e-3 | 5.9e-4 | 0.369 | 3.1e-3 | 1.007 |
| single sine | 0 | 2.5e-2 | 2.4e-2 | 0.274 | 4.1e-3 | 0.988 |
| single sine | 1 | 3.0e-2 | 2.7e-2 | 0.314 | 4.3e-3 | 0.967 |
| single sine | 2 | 6.7e-3 | 6.6e-3 | 0.226 | 2.9e-3 | 1.008 |

- The matched error decays in all six runs. The labelled error stays within 4% of its initial value.
- The rank correlation between the uniform excitation level and the fitted rate is 0.60 over the six runs. The PE
  level is a hypothesis of Theorem E, not a rate predictor (Remark E.1).
- The proof's linearisation is checked exactly: the Jacobian equals A(t), and QA + AᵀQ = diag(−k_eI, 0).

## Theorem F — transport cost (tests)

- Along (L), the particle path length (which is at least the Wasserstein length) stays below λ_G^{1/2}L₁∫‖e‖ at
  every sample.
- Along (L_σ) with the OT leak, W(μ_{t₁},μ_{t₀}) ≤ ∫(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰)) on a grid of pairs.
- Lemma P.1 holds on 2000 adversarial states, and a tight instance of F(i) reaches ratio ≥ 0.99.

## Proposition G — slow geometry (`exp_variable_projection.py`, `exp_slow_manifold_slope.py` → `variable_projection/`)

Ranges are along the run (5 samples). "a / b" means first-quarter / last-quarter medians. Velocities are centred
differences of period means.

| configuration | (H2a) ω_ref / joint fastest | (H2b) slowest weight / geometry | (H2c) S | G(a) direct check (cosine, G⁻¹ metric) | G(b) cosine | G(b) size ratio | L* | L̄ at period-mean weights |
|---|---|---|---|---|---|---|---|---|
| violating (m = 4) | 15–26 | 0.10–0.23 | 5.4–25 | 0.9994 / 0.9985 | −0.12 / 0.39 | 4.2 / 0.47 | **increases** (3.8e-3 → 5.9e-3) | non-increasing (0.219 → 0.008) |
| satisfying (m = 2) | 6.5–7.3 | 79–87 | 0.006 | 0.985 / 0.962 (geometry block 0.999 / 0.998) | 0.999 / 0.999 | 0.97 / 0.97 | decreases (1.20e-2 → 4.11e-3) | non-increasing |

Slope counterexample (`table_slope.md`). Here m = 1, k_e = 20, gains (g_a, δ, δ), and (H2b) holds in every row. The
closed loop and the joint flow start on the slow manifold of (b′).

| case | (H2b) ratio | S | closed / VarPro (b): size, cosine | closed / induced metric (b′): size, cosine | (a) / (b′) | closed / (a) |
|---|---|---|---|---|---|---|
| δ = 0.47, g_a = 20 | 81–178 | 2.0–2.8 | 0.42, 0.930 (centre moves the opposite way) | 1.12, 0.9986 | 1.015 | 1.10 |
| δ = 0.47/4, g_a = 5 (all gains ÷ 4) | 81–178 | 2.0–2.8 | 0.40, 0.916 | 1.065, 0.9999 | 1.015 | 1.05 |
| δ = 0.047, g_a = 20 | 654–1779 | 0.20–0.32 | 0.85, 0.9986 | 1.05, 1.0000 | 1.003 | 1.04 |
| δ = 0.0047, g_a = 20 | 13162–17790 | 0.020–0.023 | 1.02, 1.0000 | 1.04, 1.0000 | 1.000 | 1.04 |

- (a)/(b′) isolates the O(1/R) error of the induced-metric flow.
- closed/(a) isolates the averaging error, which shrinks with the gains, plus a 4.4% O(1/k_e) error-channel factor.

What the tables show:

- G(b) holds in the satisfying configuration, where (H2b) and S ≪ 1 both hold. The violating configuration fails both.
- The slope rows show that (H2b) alone is not sufficient.
- The induced-metric flow G(b′) tracks the closed loop at every S.
- In the satisfying configuration, the G(a) weight block is dominated by the O(1/k_e) shift of the weight equilibrium:
  cosine −0.21 / 0.50, carrying 92–95% of the residual. The geometry block matches.

## Theorem H — Langevin adaptive control (`exp_langevin_global.py`, `exp_langevin_tradeoff.py`)

The setting is the real closed loop with m = 32, a scalar plant and a periodic reference.

- The unregularised offline optimum of the occupation loss is about 7.6·10⁻⁶ (best of 5 L-BFGS restarts, an upper
  bound). One run gave 7.0e-6, because of floating-point summation order.
- The regularised optimum L̄(μ*_τ) is 0.11118 at (σ, τ) = (0.02, 0.004) and 0.00850 at (0.002, 5·10⁻⁴), grid-converged
  (`exp_gibbs_reference.py` → `langevin_global/gibbs_reference.md`).
- In the scaling σ = 0.2/k_e, τ = 0.04/k_e of Conjecture H.6, L̄(μ*_τ) = 0.11118 for every k_e. The exactly averaged
  stationary state gives L̄ = 0.11415, 0.11148 and 0.11101 at k_e = 10, 20 and 40. The L²(ρ) distance between the two
  learned functions is 6.8e-3, 1.4e-3 and 9.7e-4.

Speed–bias trade-off (σ = 0.002, 2 seeds, `langevin_tradeoff/table.md`):

| τ | escape of the collapsed start (s) | loss, collapsed | loss, spread |
|---|---|---|---|
| 0 | never | 0.166 | 0.0070 |
| 0.0005 | ∞, 1050 | 0.025 | 0.011 |
| 0.001 | 1480, 1005 | 0.022 | 0.015 |
| 0.002 | 1115, 800 | 0.029 | 0.023 |
| 0.004 | 970, 430 | 0.039 | 0.036 |
| 0.008 | 390, 295 | 0.060 | 0.065 |
| annealed (hot start) | 990, 530 | 0.032 | 0.032 |

- **Fit.** The seed-averaged escape time fits A + B·log(1/τ) with a 4.5% residual, against 25% for an Arrhenius fit.
  On the eight per-seed values both laws fit poorly (27% and 33%). The data do not identify the escape mechanism, and
  version 2's "saddle" reading is withdrawn (PARADIGM.md Remark H.2).
- **Caveats.** The threshold grows with τ, and the censored τ = 5·10⁻⁴ point is excluded from the fit.
- **Guaranteed rate.** The Holley–Stroock rate 2τρ_τ is at least e^{−8401} s⁻¹ at (0.002, 5·10⁻⁴) and at least
  e^{−1053} s⁻¹ at (0.02, 0.004).

Two fixed settings (3 seeds, `langevin_global/table_*.md`):

- **σ = 0.02, τ = 0.004.** Both starts equalise: the loss ratio is 1.03, and the final measures are 4.2× closer in W
  than without dither. The common loss is about 0.11, which matches the grid value L̄(μ*_τ) = 0.111. It is dominated by
  the regularisation bias.
- **σ = 0.002, τ = 5·10⁻⁴.** Dither frees the collapsed start (0.166 → 0.048), but the starts do not equalise within
  900 s (loss ratio 5.2).

## Remark H.0 and Proposition H.7 — the exactly averaged closed loop (tests)

Relative error of the frozen-μ averaged drift against (𝔄):

| instance | k_e = 10 | 20 | 40 | 80 | k_e × error |
|---|---|---|---|---|---|
| n = 1 (`test_langevin.py`) | 4.2% | 2.0% | 0.91% | 0.48% | ≈ 0.4 |
| n = 2, ellipse (`test_averaging.py`) | 43% | 21% | 10% | 5.1% | ≈ 4.2 |

Structure of the first correction (`test_averaging.py`):

| check | n = 1 | n = 2 (ellipse) |
|---|---|---|
| expansion error of k_e²(P_μ − P₀) against U_even + U_odd, k_e = 40 → 80 | 0.54 → 0.25 | 0.038 → 0.019 |
| circulation kernel K | ≡ 0 (asserted at 1e-12) | antisymmetric; equals the enclosed-area integral to 7e-15 (asserted at 1e-6) |
| orbit reversal (same ρ) | U_odd ≡ 0 | U_even unchanged, U_odd flipped, to 2e-15 (asserted at 1e-10) |
| k_e²·cyclic sum of the antisymmetric second variation, exact (k_e = 80) against predicted (asserted at 3%; the full matrix at 5%) | −0.0711 against −0.0704 | −0.0471 against −0.0464 |

## Remark H.7a — the orientation of the reference orbit (`exp_orbit_orientation.py` → `orbit_orientation/`)

The setup is deterministic: n = 2, m = 6, G = diag(4,4,1,1,1), σ = 2/k_e, no dither. The same elliptic reference
is traversed both ways, so the occupation measure is the same, and so is the initial configuration. The runs
converged to machine precision: the residual drift over the last 10 periods is about 1e-16. The six kernels converge
to a common configuration (kernel spread 1.05 at the start, at most 5.6e-17 at the end of every run), so the
learned network is effectively a single kernel and the test probes the effect at a collapsed equilibrium.
Predictions are linear responses at the (+) equilibrium.

| k_e | measured ‖Θ₊ − Θ₋‖ | frozen-measure prediction (H.7): size ratio, cosine | full prediction (H.7 + oscillation term): size ratio, cosine |
|---|---|---|---|
| 10 | 5.73e-3 | 2.67, 0.758 | 0.937, 0.931 |
| 20 | 2.84e-3 | 2.41, 0.830 | 0.965, 0.981 |
| 40 | 1.40e-3 | 2.29, 0.863 | 0.969, 0.994 |

- The vector error ‖measured − predicted‖/‖predicted‖ of the full prediction is 0.36, 0.20 and 0.11 at k_e = 10, 20
  and 40.
- The measured difference scales like k_e^{−1.01}, against a predicted −1.
- The full first-order theory accounts for it. The frozen-measure terms alone account for less than half.
- A first attempt with σ = 0.5/k_e and half the horizon had not converged (residual drift about a third of the
  difference). It was discarded and rerun.


## Lemma E.0 and Theorem E with explicit constants (`exp_explicit_rate.py` → `explicit_rate/`)

The explicit rate of Lemma E.0 against the exact Floquet rate of the linearised loop. It lies below the exact rate in
every row (a numerical ordering: Ψ_M, Ψ_D are grid maxima and ε a grid optimum; margins ≥ 17×). Full tables in
`results/explicit_rate/table.md`.

| instance | regime | exact / explicit |
|---|---|---|
| Remark E.1 (m = n = 1, G = I, k_e = 3), ω = 1 … 64 | Ψ_M²/k_e ≈ 1 | 954 → 17.6 (decreasing in ω) |
| the ideal of the Theorem E tests (m = 3, n = 2, k_e = 6), two-tone reference of period 20π, G = g·diag(2,2,1), g = 30 … 0.01 | fast → slow | 6.56·10¹⁰ → 4.22·10⁵ (optimised ε; 8.74·10¹⁰ → 5.63·10⁵ with the closed-form ε) |

- In Remark E.1's instance the exact rate falls like ω⁻²; the explicit rate approaches that law (log-log slope −1.56
  over ω ≥ 4, −1.92 over the last octave ω = 32 → 64). Lemma E.0(iii) derives the ω⁻² law from Ψ_D ∝ ω.
- In the slow-adaptation limit (TΨ_M² ≪ k_e) the ratio decreases towards a constant of order κ_Ψ(1 + Ψ_D/(k_eΨ_M))² for
  the optimised ε (3κ_Ψ(…)² for the closed-form bound (ii)), with κ_Ψ the condition number of the averaged Gramian
  (1.27·10⁵ for this ideal; formula values 4.9·10⁵ and 1.47·10⁶; measured 4.22·10⁵ with the optimised ε and 5.63·10⁵ with
  the closed-form ε at g = 0.01; on this instance the exact rate is 0.66 of λ_min⟨ΨᵀΨ⟩/k_e because k_e = 6 is not large
  against the reference's bandwidth, and the T²Ψ_M³ term is still 13% of c₁; the PE level is evaluated with 800
  Gauss–Legendre nodes, converged over the 20π window). The exact rate pays for ill-conditioning once, the bound twice
  (Remark E.0a).
- Lemma E.0′ (the slow regime; table 4 of the same file): on the g = 0.01 ideal (κ(M̄) = 7.8·10³; peak factors ρ = 69.9,
  ρ_D = 3.7·10³, ρ_DD = 3.7·10⁵; threshold k₀ = 226 for δ = 0.2, 321 for δ = 0.1) the guaranteed rate is within a factor
  1.26 (δ = 0.2) and 1.12 (δ = 0.1) of the exact Floquet rate at k = 400 and 600, where Lemma E.0's optimised rate is
  1.3·10⁵ times smaller; the certificate sup V̇/V ≤ −γ holds on the grid wherever the conditions hold; exact/(λ_min/k) =
  0.9885, 0.9971, 0.9993, 0.9997, 0.9998, 0.9999 at k = 50 … 600 (the k⁻² law of N₂ = 0).
- Corollary E.0′.2 (Theorem E in the slow regime; table 5): on the single-kernel instance at k_e = 40 the nonlinear rate from
  Lemma E.0′'s V is 4.53·10⁻⁴ (δ = 0.2) against 2.63·10⁻⁵ from Lemma E.0 and an exact linear rate of 1.35·10⁻³; the explicit
  ball radius is 1.53·10⁻⁶ against 6.15·10⁻⁶ (the anisotropy of the quadratic form), the overshoot 11.8 against 1.43.
- Tests (`tests/test_excitation.py`): the pointwise inequalities of Lemma E.0 on random states for four regressors, the
  envelope along trajectories, the ordering explicit ≤ exact, the ω⁻² law, the no-growth of the ratio as the gains
  shrink, and the explicit nonlinear envelope of Theorem E on the implemented closed loop (single-kernel instance with
  G = 0.3·I, where r₀ ≈ 10⁻⁵ in the w-coordinates and the explicit nonlinear rate is 3.7·10⁻⁵: the inequality holds at
  every sample over 8000 s, about 0.3 e-folds of the explicit envelope, while the state itself decays by far more).
  Lemma E.0′: the certificate and the sharpness on Ψ = (sin t, cos t) (k = 20, 50, 200; ratio → 1/(1 − δ)) and on the
  slow-adaptation ideal at k = 400 (ratio ≤ 1.5; more than 10⁴ above Lemma E.0's rate).
  Corollary E.0′.2: the explicit envelope on the implemented single-kernel loop at k_e = 40 over 1.5 e-folds (dt = 0.005,
  below the fixed-step discretisation floor).
- Corollary E.0′.3 (table 6): the explicit bound on the stationary parameter spread of the dithered (LTV), with Σ = (2·10⁻³)^{1/2}I,
  against the exact periodic stationary covariance: on Ψ = (sin t, cos t) the ratio bound/exact is 2.69, 1.58, 1.31 (δ = 0.2) and
  2.52, 1.41, 1.17 (δ = 0.1) at k = 20, 50, 200, falling toward 1/(1 − δ); on the Theorem E ideal (κ = 7.8·10³) it is 15.4–17.5,
  the price of q = 15 directions charged at the worst excitation; the exact value equals the averaged Langevin covariance k·tr X
  within 3·10⁻⁴ (ideal) and 1.25/k² (sinusoid; Σ isotropic, so ΣΣᵀ commutes with M̄ — otherwise the deviation is O(1/k)), and
  E‖e‖² = ‖Σ‖_F²/(2k) = tr(M̄X)/k to five digits (round-7 review).
- Lemma E.0″ (`exp_pe_general.py` → `explicit_rate/general_pe.md`): the exponential-average construction on two non-periodic
  regressors. Quasi-periodic Ψ = (sin t, cos √2 t): λ_T = 0.297, 0.391, 0.444 for windows 2π, 4π, 8π; conditions from k = 20;
  exact/guaranteed (δ = 0.2 / 0.1) 2.49/2.13 → 1.63/1.45 at T = 4π and 2.10/1.81 → 1.43/1.27 at T = 8π from k = 20 to 200; the
  exact rate is 0.5/k to within 1%. Chaotic Ψ = (x₁/20, (x₃ − 25)/20) along a Lorenz orbit traversed five times slower than standard
  time: λ_T = 0.015, 0.024, 0.051 for windows 5, 10, 20 against an effective excitation 0.15 (exact rates 0.00759/k·20 …);
  conditions from k = 100; exact/guaranteed 5.0/4.3 (k = 100) → 4.2/3.7 (k = 200) at T = 20, dominated by the window loss
  λ_∞/λ_T ≈ 3; the certificate holds at every admissible row.

## Theorem J — periodic orbits of the finite-m closed loop (`exp_periodic_orbit.py` → `periodic_orbit/`)

Newton shooting finds the exact periodic orbit near the critical point Θ* of J_m (residual < 10⁻¹¹); its period mean and
monodromy matrix are compared with the averaged system. Full tables in `results/periodic_orbit/table.md`.

| instance | k_e | ‖Θ̄_p − Θ*‖ | vector error of the first-order location (full / frozen only) | parameter Floquet exponents: max, median relative error |
|---|---|---|---|---|
| A (n = 2, m = 3, collapsed minimiser) | 10 / 20 / 40 / 80 / 160 | 4.09e-2 / 1.91e-2 / 9.22e-3 / 4.53e-3 / 2.24e-3 | 0.169 / 0.083 / 0.041 / 0.020 / 0.010 (frozen: 0.40 → 0.34) | 8.3e-2, 1.2e-2 → 1.3e-2, 1.1e-3 |
| B (n = 1, m = 2, two distinct kernels, fast sweep) | 10 / 20 / 40 / 80 / 160 | 6.44e-2 / 2.22e-2 / 8.37e-3 / 3.51e-3 / 1.59e-3 | 1.995 / 1.078 / 0.551 / 0.277 / 0.139 (b_osc = 0 by symmetry) | 9.5e-2, 4.0e-2 → 2.0e-3, 9.5e-4 |

- Fitted exponent of the offset against k_e: −1.05 (A), −1.33 (B; not yet asymptotic below k_e ≈ 80).
- Error channel: Liouville's trace identity along the orbit gives a mean error-channel exponent k_e − L with L = 0.16 → 0.07 (A) and
  L = −0.05 → −0.09 (B) for k_e = 10 → 160, as Theorem J(b) states (the finite-difference monodromy resolves these multipliers only as noise).
- The vector error of the first-order location halves with every doubling of k_e in both instances; in A the
  frozen-measure terms alone leave a 35% error at every k_e, so the oscillation term of Remark H.7a is an O(1) fraction.
- Tests (`tests/test_periodic.py`): ∇J_m against finite differences; shooting residual < 1e-10; the offset, the vector
  error and the Floquet error all fall by the required factors from k_e = 20 to 40 on instance B.

## Proposition K — collapsed configurations and the orientation of the orbit (`exp_orientation_stability.py` → `orientation_stability/`)

Instance A of Theorem J (n = 2, m = 3, elliptic reference traversed both ways). Floquet exponents of the collapsed
periodic orbit by Newton shooting, restricted to the separation subspace Σu_i = 0; the second-order prediction is the
linearisation of (1/k_e)⟨H⟩ + b_frozen + b_osc at its own equilibrium Θ̄₂. Full table in `results/orientation_stability/table.md`.

| | k_e = 10 | 20 | 40 | 80 |
|---|---|---|---|---|
| softest separation exponent, orientation (+) / (−), σ₀ = 0.5 | 3.563e-3 / 2.660e-3 | 1.800e-3 / 1.614e-3 | 9.277e-4 / 8.863e-4 | 4.737e-4 / 4.640e-4 |
| difference (+) − (−) measured / predicted at Θ̄₂ / at Θ* | 9.03e-4 / 7.77e-4 / 9.01e-4 | 1.86e-4 / 1.69e-4 / 2.08e-4 | 4.14e-5 / 3.93e-5 / 4.98e-5 | 9.71e-6 / 9.44e-6 / 1.22e-5 |
| k_e² × difference, measured / predicted at Θ̄₂ | 0.090 / 0.078 | 0.075 / 0.068 | 0.066 / 0.063 | 0.062 / 0.060 |

- Fitted exponent of the difference against k_e: −2.18 (Proposition K(iii): −2); the prediction at Θ̄₂ converges to the
  measured k_e²-scaled limit, the prediction at Θ* (third-derivative term omitted) converges to a value 30% too large.
- At σ_crit = 0.4606 both orientations are unstable (−4.1e-4 and −1.27e-3 at k_e = 10); the even part dominates.
- Designed σ₀ (pred₊ + pred₋ = 0): σ₀ − σ_crit = 0.0125 (k_e = 10), 0.0061 (k_e = 20); measured softest exponents
  +8.5e-4 / −2.5e-5 and +1.4e-4 / −3.7e-5: one orientation stable, the other unstable.
- Test (`tests/test_periodic.py::test_proposition_K_hessian_splitting_and_critical_leak`): the Hessian splitting at the
  collapsed point (1e-5), σ_crit's gap equation (1e-8), the two-fold zero separation mode at σ_crit.

## Theorem D′ — propagation of chaos for the dithered closed loop (`exp_chaos.py` → `chaos/`)

Synchronous coupling of the m-particle loop (n = 1, k_e = 10, σ = 0.02, T = 20) to an M = 4096 reference that shares its
first m particles and Brownian motions; sup_t of the coupled distance D_m = ‖x^m − x^M‖ + (1/m)Σ‖ϑ^m_i − ϑ^M_i‖, 24 seeds.
Full table in `results/chaos/table.md` (raw per-seed values in `raw_sup_D.npz`).

| | τ = 0 | τ = 2·10⁻³ |
|---|---|---|
| mean sup_t D_m, m = 16 → 1024 | 1.49e-1 → 2.17e-2 | 1.52e-1 → 2.32e-2 |
| slope of the mean against m (m ≤ 256), bootstrap s.e. | −0.47 ± 0.06 | −0.46 ± 0.05 |
| slope of the median | −0.46 | −0.45 |
| slope against (1/m − 1/M)^{1/2} (predicted 1) | 0.94 ± 0.08 | 0.92 ± 0.08 |

- The m^{−1/2} rate of Theorem D′(b) holds within one standard error at both temperatures; the shared-particle coupling
  predicts the mean ∝ (1/m − 1/M)^{1/2}, which the data follow up to m = M/4.
- Single seeds differ by more than an order of magnitude at fixed m; a three-seed version of the experiment gave slopes
  −0.25 and −0.21 and was withdrawn (review round 5, D2). The constant of the theorem (T K_Φ e^{ΛT}) is not checked.
- No test: the experiment is the check.
- Theorem D″ (`exp_chaos_uniform.py` → `chaos/uniform_table.md`; n = 1, k_e = 25, G = 0.05·diag(4,1,1), a two-bump plant that
  is not odd, M = 2048): exact Lipschitz constants L_x = L_ϑ = 8.087, L₂ = L_{xϑ} = 44.4 and the a-priori tube E_X = 0.228 give
  σ_cert = 4.04, against Proposition K's critical leak σ_K = 0.00144 of the collapsed critical point (ratio 2800). At 1.2σ_cert
  with dither τ = 0.02 (80 s, 8 seeds) the coupled distance is 6.7·10⁻⁴ → 7.9·10⁻⁴ (m = 16) and 1.6·10⁻⁴ → 1.9·10⁻⁴ (m = 256)
  between the first and the last window, slope −0.51 ± 0.02, bound 0.55 and 0.14; the network stays at its prior (fit error
  1.000, rms spread 0.033 = the dither cloud). Without dither the reference collapses (spread 10⁻¹⁴) and the distance is zero;
  at σ_cert/8 it decays to 10⁻¹¹. In the learning regime (σ_K/3, 1000 s, 6 seeds) the network stays spread (rms 1.53) and fits
  the plant (fit error 0.632), and the distance follows m^{−1/2} (slope −0.62 ± 0.16) but grows threefold between the windows
  (100, 200] and (500, 1000]: no empirical uniformity there (the first draft's "spread regime" run was a slow collapse; review
  round 6).

## Theorem D‴ — uniform-in-time propagation of chaos in the learning regime (`exp_chaos_learning.py` → `chaos/learning_table.md`)

- Instance: n = 1, G = diag(1, 0.25, 0.25), three-bump plant (1, −1.5, 0.5), (−0.7, 0, 0.4), (0.8, 1.5, 0.5), x_m = 2 sin t, σ_crit = 0.04716
  (Proposition K), σ₀ = 0.15σ_crit = 0.00707; L_ϑ^G = 4.044, L₂^G = 11.11, K_Φ = 4, F_ρ = 0.999. Averaged slow loop (𝔄_m), m = 16, 64, 256
  against M = 2048, 900 slow units, 4 seeds per start; real loop at k_e = 25 over 2000 time units, 3 seeds.
- The hypothesis at the three-atom stationary measure (fit error 0.244): κ_R at the atoms 0.00004, −0.00227, 0.00411; over spheres of
  radius 0.05, 0.1, 0.2, 0.3: 0.00461, 0.00512, 0.00619, 0.00730 — (H) holds with ϱ′ = 0.1 and margin 0.00195 (1/κ = 512), the window
  κ_R[μ*](N_{0.1}) = 0.00512 < σ₀ = 0.00707 < σ_crit = 0.04716 is realised; at 0.3σ_crit the three-atom measure fails splitting
  (κ_R = 0.01696 > 0.01415). Theorem constants: C₁ = 2.2·10⁹, m₀ = 1.2·10²³ (vacuous envelopes).
- Along the warm run (jittered mixture around the atoms) the margin a(s) stays in [+0.00139, +0.00167]; the own-residual statistic
  exceeds the r*-one by ≤ 0.00013 ((0.2) allows 0.214); rms drift 6·10⁻³ → 3·10⁻⁶. Along the cold run (connected cloud) the margin is
  negative throughout (−0.196 → −0.0025) while the fit converges to 0.244 as well.
- Warm start: window maxima of D_m 0.396 → 0.448 (m = 16), 0.122 → 0.138 (64), 0.0119 → 0.0137 (256) over (0, 112.5] → (450, 900]
  (last/first 1.13–1.15, last/previous 1.01–1.02: saturated); slope −1.28 ± 0.30 (mean square; theorem −1), learned function −0.69 ± 0.16;
  per-particle maxima 1.36, 0.62, 0.21 ≫ 3ϱ′/4 = 0.075 (Remark D‴.6).
- Proposition D‴.5: ‖S_k‖ = 8.38, 4.66, 8.32; E D_∞ predicted 0.823, 0.201, 0.0454 against means 0.446 ± 0.28, 0.138 ± 0.051,
  0.0136 ± 0.0055 (four seeds, heavy-tailed); per seed the linearised prediction is within 25% of the measured D in 9 of 12 pairs
  (worst ratio 0.57 at |w^m − w^M| = 0.275) and the exact stationary distance of the replicated weighted configuration within 2.3% in
  all 12; splitting margins positive in 11 of 12 (+0.0006 to +0.0050), negative (−0.0022) for m = 16, seed 3 (weights (9,4,3)/16: a
  saddle, left later by particle migration — review round 9); per-cluster displacements along S_k(w^m − w^M) with cosines 0.92–1.000
  and norms within 0.8–20%.
- Cold start: D grows 1.81 → 5.33 (16), 0.505 → 1.79 (64), 0.0744 → 0.583 (256) (last/first 2.95, 3.54, 7.84), slope −0.87 ± 0.14;
  bimodal deviations (medians 0.40, 0.14, 0.067; maxima 6.2–6.7); fraction of pairs > 1 whitened unit apart 0.203 ± 0.047,
  0.082 ± 0.021, 0.037 ± 0.010, decaying like m^{−0.64±0.08} (medians m^{−0.62±0.07}): the floor of Remark D‴.2.
- Real loop (k_e = 25, warm start): window maxima 0.328 → 0.590 (16), 0.0605 → 0.0909 (64), 0.0188 → 0.0294 (256) over (0, 500] →
  (1000, 2000] (slow time ≤ 80 = 0.16/κ), slope −1.13 ± 0.15; on the same initial clouds the averaged flow gives D_real/D_averaged
  between 0.993 and 1.004 at every window end (across different clouds the last-window ratios are 1.57, 0.78, 2.58): the real loop
  follows its averaged slow dynamics to within 0.7% (Corollary D‴.7); two of the three m = 16 clouds have inadmissible or marginal
  weights (margins −0.0021, −0.00001).

## Proposition D‴.9 — measure-valued priors and the design of the window (`exp_prior_design.py` → `chaos/prior_design.md`)

- Same instance and leak as Theorem D‴'s numerics (σ₀ = 0.15σ_crit = 0.00707), nearest-atom leak toward three prior atoms. Matching prior
  (exact representation): fit error 0, κ_R = 0, margin σ₀ = 0.00707 against the σ-modification's 0.00196 (the theorem's relaxation time 1/κ 141 against 511 slow units; observed late-time e-folding of the coupled distance 48–84 against 107–306); ‖S_k‖ = 13.01, 6.18, 8.91 against 8.38, 4.66, 8.32; E D_∞(m = 64) 0.323 against 0.201 (1.61×); function-level floor 9.33·10⁻⁴
  against 4.29·10⁻⁴ (2.18×).
- Misplaced priors: centres +0.15/+0.3/+0.45/+0.6 → margins +0.0055/+0.0037/+0.0017/−0.0006; amplitudes ×0.5/×0.7/×1.2 → +0.0074/+0.0073/−0.0317;
  widths ×0.7/×1.5/×2 → −0.0261/+0.0058/+0.0040; the shrunk priors (×0.5, ×0.7) have the largest margins and lower floors (E D_∞ 0.199, 0.230;
  function level 7.1·10⁻⁴, 7.7·10⁻⁴) at fit errors 0.09, 0.06.
- Admissible weights (splitting margin of Θ*(w)): (7,5,4)/16: +0.0026 / −0.0054 / +0.0073; (9,4,3)/16: −0.0022 / −0.0058 / +0.0057;
  (10,3,3)/16: −0.0038 / −0.0042 / +0.0054; (0.5,0.3,0.2): −0.0006 / −0.0065 / +0.0066 (σ-modification / matching / amplitudes ×0.7): the
  shrunk prior keeps all four admissible (above σ₀ only at (7,5,4)/16; its over-weighted atom overshoots by 1–3% at the other three), and
  the matching prior's set is smallest toward the largest bump (an excess 0.04 on atom 1 alone is inadmissible, 0.08 on atom 2 or 3 alone is not).
- Matching-prior coupled runs (m = 16, 64, 256 vs M = 2048, 300 slow units, 3 seeds): D saturates (last/previous 1.04–1.05), slope −1.23 ± 0.40;
  per-seed linearised-law ratios 0.81–1.11 (parameter level) and 0.61–1.69 (function level), exact law 0.987–0.998 in all nine pairs; two of
  nine realised weight vectors are saddles (margins −0.0005, −0.0030), left by a local split inside the cell with the coupled distance on its
  plateau (round 10); observed e-folding times 48–84 slow units (median 59) against 107–306 (median 195) for the σ-modification.
- Function-level law for the σ-modification against Theorem D‴'s warm runs: predicted RMS 0.0419, 0.0207, 0.0098 vs measured means
  0.0307 ± 0.0097, 0.0154 ± 0.0033, 0.0047 ± 0.0012 (ratios 0.73, 0.74, 0.48; RMS over seeds/predicted 0.84, 0.79, 0.52; per-seed law from the
  realised weights 0.60–1.57).

## Proposition D‴.12 — the OT leak's mass constraint against the weight fluctuation (`exp_ot_leak.py` → `chaos/ot_leak.md`)

- Same instance and matching prior as Proposition D‴.9's numerics (σ₀ = 0.15σ_crit, q = (1/3, 1/3, 1/3), Δ_ν = 6.26); the nearest-atom
  and the barycentric OT leak on the same warm clouds, 600 slow units (4.2/σ₀), 4 seeds; m = 15, 60, 240 vs M = 510 (commensurate) and
  m = 16, 64, 256 vs M = 512 (incommensurate).
- Cluster weights: conserved exactly under the nearest-atom leak (max deviation 0 in 32 systems); under the OT leak the mass in transit
  ½|w(0) − q|₁ = 0.117, 0.058, 0.027, 0.024 (seed means; 1–21 particles per system) reaches zero and the weights equal q from s_eq = 93–164
  slow units (ln 2/σ₀ = 98). Fit error at the end (seed means): 0.067, 0.029, 0.017, 0.012 (nearest) → 2.8·10⁻⁴, 1.2·10⁻⁴, 1.2·10⁻⁴,
  0.8·10⁻⁴ (OT; the migrated particles still 0.06–0.16 from their atoms, the OT finals not yet stationary).
- Commensurate sizes, last window: function-level distance² 2.5·10⁻³, 5.2·10⁻⁴, 5.9·10⁻⁵ (nearest; law predicts 4.0·10⁻³, 9.1·10⁻⁴,
  1.4·10⁻⁴) vs 4.9·10⁻⁷, 4.9·10⁻⁸, 7.5·10⁻⁹ (OT, still falling ×0.02–0.035 per 150 slow units); W² 5.1, 2.7, 0.81 (nearest; cross-cluster
  transport ½|w^m − w^M|₁Δ_ν² = 4.9, 2.5, 0.70; slope −0.67 ± 0.13) vs 7.9·10⁻³, 2.3·10⁻³, 6.9·10⁻⁴ (OT, ×0.12 per 150 slow units);
  synchronous D 0.87, 0.25, 0.017 (nearest) vs 5.8, 3.3, 0.94 (OT: the reassigned particles; at 0.95–0.96 of the mismatch plateau and
  still rising). Exact K × K transport between the weighted stationary measures at the realised weights: 0.997–0.999 of the measured W² at
  the seven admissible pairs.
- Incommensurate sizes: one or two split (bridge) particles per system at 0.28–0.58 Δ_ν from the nearest atom; against ν (the mean-field
  stationary measure): W²(μ_m, ν) = 1.21, 0.33, 0.079, 0.038, m·W² = 19.4–21.0 (slope −0.98 ± 0.03), fit error 0.028, 0.0093, 0.0023, 0.0011,
  m·fit = 0.45–0.60 (squared slope −1.81 ± 0.05); against the M = 512 reference (two bridge particles of its own) OT 3.4·10⁻⁴, 3.1·10⁻⁵,
  1.1·10⁻⁶ (function level; m = 16 still relaxing, plateau 2.8·10⁻⁴ in round 11's continued run) and 1.16, 0.29, 0.041 (W²) vs nearest
  2.5·10⁻³, 5.0·10⁻⁴, 4.4·10⁻⁵ and 4.8, 2.7, 0.69; nearest fit errors 0.061, 0.035, 0.013, 0.012.
- Masses as a design input: q_δ = (1/3 + δ, 1/3 − δ/2, 1/3 − δ/2) gives fit error 0.0069–0.049 for δ = 0.01–0.1, splitting margin +0.0048,
  +0.0028, +0.0011, −0.0015, −0.0049 (inadmissible between δ = 0.03 and 0.05), diagonal plan unique at every row (ε ≤ 1.04 < 3.13).

Tests (`tests/test_ot_leak.py`): the plan is a vertex with the prescribed masses (integral when commensurate, ≤ K₀ − 1 split rows
otherwise) and the single-atom prior is the σ-modification; the leak is the G-gradient flow of J_m^{OT} against finite differences; the
replicated weighted configuration is stationary for three commensurate sizes with the diagonal plan and the same learned function (W² = 0
between sizes), for matching and mismatched masses; conservation against constraint on a short run; the forest bound; the exact W² by
replication against the transportation LP.

## Corollary D‴.14 — the real loop with the OT leak (`exp_ot_real_loop.py` → `chaos/ot_real_loop.md`)

- Real loop at k_e = 25 (σ = σ₀/k_e, RK4 plant with the control and the reference evaluated at the stages, Euler parameters, dt = 10⁻³,
  100 time units = 16 periods) from the contracted clouds of Proposition D‴.12's averaged runs (seeds 0, 1; m = 15, 60, M = 510) and from
  the collapsed lifts; last five periods on the reference orbit.
- OT leak: ‖f̂_m − f̂_M‖²_{L²(orbit)} = 7.5·10⁻⁹–2.7·10⁻⁸, against the 3-atom lift 5.2·10⁻¹⁰–3.6·10⁻⁸; fit error along the orbit
  3.5·10⁻⁵–2.9·10⁻⁴; rms tracking error 0.9·10⁻⁶–7.3·10⁻⁶ (the lift 8·10⁻¹²); W²(μ_m, μ_M) 3.3·10⁻⁴–1.4·10⁻³; within-class spreads 0.013–0.033.
- Nearest-atom leak: ‖f̂_m − f̂_M‖² = 1.1·10⁻⁴–2.5·10⁻³; fit errors 0.012–0.075; rms tracking errors 3.2·10⁻⁴–1.9·10⁻³; W²(μ_m, μ_M) 1.6–5.5.
  Three of six systems have inadmissible realised weights (margins −0.0035, −0.0011, −0.0022): their lifts are saddle orbits and the systems
  have split (spread 0.04–0.11; lift-to-lift ratios 0.960, 0.998, 0.961); at the one pair with both weights admissible the disagreement is
  the lifted orbits' distance to 1.0003, the admissible systems within 2.5·10⁻¹¹ of their lifts (spread 0.0006–0.0013).

Tests (`tests/test_ot_leak.py`): the collapsed 15-system and the 3-atom loop integrate identically under the OT leak at masses q (10⁻¹²
over 300 steps); the collapsed system with weights (6, 5, 4)/15 stays collapsed under the nearest-atom leak and loses a particle under the
OT leak; two commensurate lifts keep W² = 0 and the same learned function; the lifted orbit's class separation subspaces are invariant and
carry the ghost spectrum (Proposition D‴.15, m = 6).

## Corollary D‴.14 in two dimensions (`exp_ot_2d.py` → `chaos/ot_2d.md`)

- Planar plant (four Gaussian bumps with vector amplitudes), Lissajous reference, p = 5 parameters per kernel, matching prior (masses ¼),
  k_e = 25, σ₀ = 0.05, 2500 time units = 100 slow units, sizes 8, 32 vs 128 from the same jittered clouds, two seeds.
- OT leak: weights equalised at 13–14 slow units (ln 2/σ₀ = 14); ‖f̂_m − f̂_M‖² = 5.0·10⁻⁸–1.6·10⁻⁷; fit error along the orbit
  2.1·10⁻⁴–8.7·10⁻⁴; rms tracking error 3.9·10⁻⁶–1.6·10⁻⁵; W²(μ_m, μ_M) 2.0·10⁻⁴–4.4·10⁻⁴.
- Nearest-atom leak: weights conserved; ‖f̂_m − f̂_M‖² = 5.5·10⁻³–2.3·10⁻²; fit errors 0.11–0.39; rms tracking errors 1.9·10⁻³–6.7·10⁻³;
  W²(μ_m, μ_M) 4.7–11.8.

## Proposition D‴.15 — the transverse spectrum of a lifted orbit (`exp_lifted_orbit.py` → `chaos/lifted_orbit.md`)

- Real loop at k_e = 25; monodromies by central differences (DOP853, rtol 10⁻¹¹). OT leak at the matching masses: the 3-atom orbit (the exact
  representation, zero tracking error) lifted to m = 6, 9, 15 — return residuals ≤ 1.7·10⁻¹⁵; separation exponents σ₀/k_e = 2.830·10⁻⁴ in
  every class (= ghosts; invariance defects ≤ 1.8·10⁻⁸); collective exponents 2.935·10⁻⁴ … 7.934·10⁻³ identical across m to four digits,
  Theorem J(b)'s prediction 2.952·10⁻⁴ … 8.091·10⁻³ (0.6–3.2%).
- σ-modification at (7, 5, 4)/16 (margin +0.00257) and (9, 4, 3)/16 (−0.00219), lifts to m = 16 (residuals ≤ 4.9·10⁻¹⁴): separation exponents =
  ghosts to 2·10⁻⁵ relative, = (σ₀ − λ_j(R^G_k))/k_e within 0.2–4.1% (defects ≤ 5.7·10⁻⁸); the saddle's class 0 exponent −9.13·10⁻⁵ (predicted
  −8.77·10⁻⁵, multiplicity 8); period-start offsets from Θ*(w) 1.1·10⁻², 1.6·10⁻² (Theorem J(a)).

## Theorem D‴.18 — Theorem J uniformly in m (`exp_uniform_basin.py` → `chaos/uniform_basin.md`)

- Real loop at k_e = 25 (σ₀ = 0.15σ_crit, dt = 5·10⁻³), 8000 time units, many systems integrated together.
- Matching prior at the matching masses (constant lift), sizes 3, 15, 60, 240, per-particle random displacements ρ: inside the cells
  (ρ = 0.3, 1.0, 2.5) both leaks identical, fitted decay rates of the per-particle distance 2.69–2.97·10⁻⁴ for m ≥ 15 (predicted
  σ₀/k_e = 2.83·10⁻⁴; ratios 0.95–1.05), 2.95–4.06·10⁻⁴ at m = 3 (no separation modes; collective 2.95·10⁻⁴), constants C = 1.00–1.29.
- Across a cell boundary (ρ = 4.5): nearest-atom leak ends at counts (2, 0, 1), (6, 5, 4), (23, 19, 18), (87, 70, 83), W² 12.3–1.9; OT leak
  at q for every m, W² 0.010–0.23 (still falling).
- Moving orbit (OT leak, prior centres +0.15): rates 2.18–2.63·10⁻⁴ against 2.27·10⁻⁴ (ratios 0.96–1.16), constants 1.04–1.23.

Tests (`tests/test_uniform_basin.py`): the drift's second-order remainder ≤ C_N|||u|||² with C_N ≤ 0.04 for random per-particle displacements at
m = 6, 30, 150, 600 and exactly 0.064 for a collective one; a one-particle separation moves nothing else at first order (response −σ₀u at
m = 6, 60, 300); the averaged loop contracts |||·||| at 0.9σ₀ with one constant for m = 6, 30, 150.

## Remark D‴.16 — dither and the mass constraint (`exp_dither_ot.py` → `chaos/dither_ot.md`)

- Dithered averaged slow loop (Langevin in slow time, Gibbs law exp(−J_m/τ₀)), m = 15 from the counts (6, 5, 4), matching prior,
  σ₀ = 0.15σ_crit; τ₀ = 0.005, 0.01, 0.02 over 50 000, 200 000, 50 000 slow units; two seeds.
- Nearest-atom leak: hop rates 1.0·10⁻⁵, 5.1·10⁻⁴, 5.6·10⁻³ per slow unit (Arrhenius barriers 0.039–0.048; ridge estimate σ₀Δ_ν²/8 = 0.035);
  time at the balanced counts 0.00, 0.54, 0.25 against the leading-order Gibbs law 0.86, 0.58, 0.33 (the Laplace determinant factors
  over-weight soft separation modes and are not used).
- OT leak: balanced counts 99.6%, 99.2%, 88% of the time; label changes 1–364 per run (swaps followed by transits); the plan's masses stay q.

## Theorem C′ — the barycentric OT leak (`exp_barycentric_leak.py` → `barycentric_leak/`)

| experiment | result |
|---|---|
| pure leak, m = 6 toward a 48-atom prior, 3 starts | W² non-increasing at every step; end point centroidal (V_bar² ≤ 2.2e-25, W² = Var_γ = 1.164, 1.213, 1.203); a Lloyd fixed point; best of 20 Lloyd runs 1.171 (between the end points: a local minimiser, as for Lloyd) |
| m₀ = m | exact decay of Theorem C(i) to 5e-16 |
| closed loop, representable plant, m = 4, k_e = 4, σ_m = 0.5: (A) OT leak to the ideal / (B) barycentric leak to a 32-atom cloud / (C) labelled leak to a random 4-subset | W(μ_T, μ*) = 0 / 0.167, 0.097, 0.157 (= b_bar = 0.168, 0.097, 0.159, against b = W(μ*, cloud) = 0.45, 0.33, 0.40) / 0.61, 1.02, 0.35 (≈ W(μ*, subset)); RMS‖e‖ on [20, 40]: 0 / 0.002–0.007 / 0.03–0.04 |

Tests (`tests/test_barycentric.py`): the plan reduces to the assignment for m₀ = m and W² = V_bar² + Var_γ; the pure leak
dissipates at every step with D⁺Z ≤ −σ_mV_bar² and ends at a Lloyd fixed point; √C_bar ≤ b + 2η on 90 configurations with
C_bar = b_bar² at small η; the Lyapunov inequality (iii) along a closed loop with a 24-atom prior.

## Theorem J′ — the stationary law of the dithered single-kernel loop (`exp_dither_stationary.py` → `dither/`)

Single-kernel instance (n = 1, m = 1, ϑ = (α, γ), frozen width, G = diag(2, 2), σ₀ = 0.4, τ₀ = 0.1), reference traversed both
ways; 2·10⁴ independent replicas per (k_e, orientation) started from the first-order law, k_e = 10, 20, 40, 80; the Gibbs
measure, the first-order correction g, the current j₁ and the circulation coefficient on a 201² grid. Full tables in
`results/dither/table.md` (raw in `raw.npz`).

| | k_e = 10 | 20 | 40 | 80 |
|---|---|---|---|---|
| (m_ε(α) − π(α))·k_e, prediction −0.255 | −0.28 | −0.36 | −0.36 | −0.37 |
| (m_ε(γ) − π(γ))·k_e, prediction −0.426 | −0.44 | −0.54 | −0.51 | −0.67 |
| residual m_ε(α) − first order (s.e. 0.0013) | −0.0021 | −0.0051 | −0.0025 | −0.0015 |
| residual m_ε(γ) − first order (s.e. 0.0018) | −0.0015 | −0.0056 | −0.0020 | −0.0030 |
| orientation-odd shift of γ, measured/predicted | 0.93 ± 0.04 | 1.10 ± 0.07 | 1.15 ± 0.15 | 1.38 ± 0.29 |
| k_e²E e²/(2π(L̄)), orientation +1 | 1.054 | 1.033 | 1.017 | 1.009 |
| circulation Λ_ε (+1), measured ± 0.0005 / predicted | +0.0037 / −0.0005 | +0.0006 / −0.0003 | −0.0004 / −0.0001 | +0.0003 / −0.0001 |

- The Gibbs discrepancy decays like k_e^{−0.87} (α), k_e^{−0.83} (γ), k_e^{−1.09} (αγ), and the first-order prediction removes
  93–97% of it at k_e = 10 and 70–85% at k_e = 20–40; the residual is at the Monte Carlo floor at k_e = 80 and its second-order
  decay is not resolved.
- The orientation-odd shift (the dither version of Remark H.7a) is predicted to 7–8% at k_e = 10 and within the 10–30%
  resolution at larger k_e.
- The tracking error is the Gibbs-averaged occupation loss over k_e²/2 with an O(1/k_e) correction of 0.5–0.7/k_e.
- The first-order circulation is below the Monte Carlo resolution for k_e ≥ 20 and dominated by higher orders at k_e = 10; the
  current machinery of Proposition J′.2 is checked on the first-order averaged SDE on the grid: ∫Θ∧(j_ε − j₀)/(ε∫Θ∧j₁) = 1.47,
  1.21, 1.09, 1.03 at ε = 0.2, 0.1, 0.05, 0.025 and ‖m_ε − m₀(1 + εg)‖₁/ε² = 0.31 constant.
- Referee round 6 (`review/REVIEW_round6_Jprime.md`): the first run's batch standard errors were biased low by 1.4–1.9 and its
  burn-in too short (JP-4–JP-6); the run reported here is the corrected one.
- Corollary J′.3 (`exp_gibbs_m.py` → `dither/gibbs_m.md`): μ*_τ on the grid (∫α dμ* = 0.4963, ∫γ dμ* = 0.1873, ∫α² dμ* = 0.9353,
  L̄(μ*) = 0.1394); the Gibbs measures π_m for m = 1–32 by MALA (4000 replicas): E[L̄(μ_Θ)] = 0.1960, 0.1681, 0.1540, 0.1470, 0.1432,
  0.1413 with (E[L̄] − L̄(μ*))·m = 0.057–0.061 (slope −0.99), E[γ₁] slope −1.01, E[α₁²] from 0.536 to 0.922 with rescaled deviation
  −0.41 ± 0.01; the m = 1 row matches the quadrature of π₁ within one standard error.

## Theorem I — reciprocity at order k_e⁻² (tests)

No new experiment: every statement is checked in `tests/test_nonreciprocity.py` against the generic second-variation
code of `src/core/averaging.py`.

| check | result |
|---|---|
| A(ϑ,ϑ') = aᵀB(ξ,ξ')a' against the antisymmetric part of S_even + K, n = 1 and 2 | agreement to 1e-10 |
| A(ϑ, ϑ_z) = 0 for the zero feature; cyclic sum over (ϑ₀, ϑ₁, ϑ_z) equals A(ϑ₀, ϑ₁) | 1e-14 |
| n = 1, f := f̂_μ (perfect fit on the orbit) | max |A| < 1e-13 over 20 random pairs |
| n = 1, non-zero residual | |A| > 1e-3 found; equals the closed form −aa'((c−c')/s²)∫rφ_cφ_c'dρ to 1e-12 |
| κ(ξ,ξ') against −(1/(2s²))e^{−|c−c'|²/(4s²)}(c−c')·F_s(m), n = 2 | 1e-10 |
| F₂((0.5, 0)) > 0 on the ellipse; κ ≡ 0 on a segment traversed back and forth | F₂ > 1e-3; |κ| < 1e-12 |
| n = 2 with f ≡ f̂_μ on ℝ² (residual and its Jacobian vanish on the orbit): B = 2κI | 1e-13; A ≠ 0 on the ellipse, A ≡ 0 (< 1e-13) on the segment |
| n = 2 with f = f̂_μ + he₁, h zero on the ellipse (residual zero on the orbit, Jacobian not symmetric): B − 2κI = −∫φφ'∂₂h dρ·(e₁e₂ᵀ − e₂e₁ᵀ) | 1e-12, and not small (the round-5 correction of Theorem I(iv)) |
| positivity witness f ≡ 0, μ = δ_ϑ̄ | integral > 0; generic code equals the closed form to 1e-12 |
