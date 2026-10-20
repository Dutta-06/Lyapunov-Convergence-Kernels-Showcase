# Transport adaptive control: adaptive approximators as evolving measures

**Version 4**: version 3 after review rounds 1–4, plus the results of the research continuation reviewed in
`review/REVIEW_LOG.md` round 5 (Lemma E.0, Theorem E with explicit constants, Theorems C′, D′, I and J, Proposition K) and round 6
(Lemma E.0′ with Corollary E.0′.2, Theorem D″, Theorem J′ with Proposition J′.2 and Corollary J′.3). Version 2 is
archived as `review/PARADIGM_v2.md`; the review reports of rounds 1–3 cite its line numbers.

**Companion files.**

- **Code.**
  - `src/core/transport.py` implements laws (L) and (L_σ) and the oracle Lyapunov report.
  - `src/core/langevin.py` implements Theorem H.
  - `src/core/averaging.py` implements the frozen-measure averaging of Proposition H.7.
  - `src/core/excitation.py` implements the strict Lyapunov function of Lemma E.0 and the explicit constants of Theorem E.
  - `src/core/nonreciprocity.py` implements the reciprocity criterion of Theorem I.
  - `src/core/periodic.py` implements the periodic orbits of the finite-m closed loop and the averaged objective J_m of Theorem J, and the collapsed-configuration analysis of Proposition K.
  - `src/core/transport.py` also implements the optimal plan, the barycentric projection and the leak (L_σ,bar) of Theorem C′.
  - `src/core/dither.py` implements the single-kernel instance, the averaged data (J, ∇J, b, ⟨u∧u̇⟩), the Gibbs grid with its Poisson,
    current and stationary Fokker–Planck solvers, and the vectorised dithered loop of Theorem J′.
  - `src/core/chaos_learning.py` implements the averaged slow loop (𝔄_m), the three terms of its coupling identity, the
    residual-weighted curvature κ_R and the constants of Theorem D‴; `src/core/prior_design.py` the nearest-atom leak toward a
    measure-valued prior and the two weight-sensitivity laws of Proposition D‴.9; `src/core/ot_leak.py` the barycentric OT leak toward
    a prior measure (the vertex plan, the Laguerre/diagonal-plan checks, the exact W² between empirical measures) of Proposition D‴.12;
    `src/core/lifted_orbit.py` the real loop with per-class leak targets, the lift of a K₀-atom orbit, the class separation bases and the
    ghost monodromy of Corollary D‴.14 and Proposition D‴.15; `src/core/entropic_leak.py` the entropic plan, the entropic leak and its
    separation blocks of Proposition D‴.20; `src/core/ot_basin.py` the identities along the plan segments, the plan-segment curvatures and the
    Polyak–Łojasiewicz constants of Theorem D‴.21; `src/core/adaptive_prior.py` the adaptive prior (the joint drift and energy, the two
    gates, the Gauss–Newton block and the pair rates; the nearest-atom option of Proposition D‴.24; the projected prior of Proposition D‴.27) of Theorem D‴.23 (Proposition D‴.28's gauge in the tests and `exp_gauge_floor.py`); `src/core/floor_bias.py`
    the exact weighted representation, the bias form of the floor law and its decomposition by modes of Proposition D‴.24;
    `src/core/planar_loop.py` the averaged slow loop on ℝⁿ (n ≥ 1) with the nearest-atom leak toward fixed and adaptive priors, the weighted
    stationary point and the Gauss–Newton block, for the planar check of Proposition D‴.24; `src/core/adaptive_floquet.py` the real loop with
    the adaptive prior in the state, the lift of perfect tracking, its monodromy and the predicted exponents of Corollary D‴.26, and the
    variational equation on the collapsed subspace with its monodromy, the error channel by Liouville's formula and the averaged generator to
    second order in 1/k_e (the coefficients c_j); `src/core/sliding.py` the second-order plant under the sliding-variable controller of
    Proposition D‴.29 (the phase-plane network, the averaged slow loop along the regressor's orbit, the real servo, the filtered error);
    `src/core/measurement_noise.py` the noise matrices, the floors and the second-order force of Proposition D‴.30, the exact periodic
    and transient covariances of the linearised loop, the real loop with an Ornstein–Uhlenbeck measurement noise, and the separations' second-order drift at the lift and along the
    second-order transient for each regressor, its cross term, the transient carried to the stationary point by one-period maps
    (`long_transient`) and a separation's variance term (Proposition D‴.30(i), Corollary D‴.33); and Proposition D‴.31's desired-regressor variants of the real loop
    and of the second-order force, the washout block of the prior family, its zero-frequency operators and the response bounds; `src/core/cohesion_leak.py` the cohesion leak of Proposition D‴.32 (the kernel and
    its energy, the averaged loop, the clusters, the predicted spectrum at a replicated configuration, the real loop with the lift); `src/core/mass_transfer.py` the replicator weight law
    of Proposition D‴.34 (the weighted loop and its energy, the run on (Θ, log β), the class masses, the real loop at fixed positions with
    the relative-entropy Lyapunov function).
- **Tests.** `tests/test_transport.py`, `tests/test_langevin.py`, `tests/test_averaging.py`, `tests/test_excitation.py`,
  `tests/test_nonreciprocity.py`, `tests/test_periodic.py`, `tests/test_barycentric.py`, `tests/test_dither.py`, `tests/test_chaos_learning.py`, `tests/test_prior_design.py`, `tests/test_ot_leak.py`, `tests/test_uniform_basin.py`, `tests/test_entropic_leak.py`, `tests/test_ot_basin.py`, `tests/test_adaptive_prior.py`, `tests/test_adaptive_floor.py`, `tests/test_planar_loop.py`, `tests/test_adaptive_floquet.py`, `tests/test_gauge_floor.py`, `tests/test_sliding.py`, `tests/test_measurement_noise.py` and `tests/test_eiv_remedy.py`, `tests/test_cohesion_leak.py`, `tests/test_mass_transfer.py`, `tests/test_split_bifurcation.py`. Run them with
  `venv\Scripts\python -m pytest tests`.
- **Experiments.** The scripts are `experiments/exp_transport_*.py`, `exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py`, `exp_mean_field_m.py`,
  `exp_pe_convergence.py`, `exp_variable_projection.py`, `exp_slow_manifold_slope.py`, `exp_langevin_global.py`,
  `exp_langevin_tradeoff.py`, `exp_gibbs_reference.py`, `exp_orbit_orientation.py`, `exp_explicit_rate.py` and
  `exp_periodic_orbit.py`, `exp_orientation_stability.py`, `exp_chaos.py`, `exp_chaos_uniform.py`, `exp_chaos_learning.py`, `exp_prior_design.py`, `exp_ot_leak.py`, `exp_ot_real_loop.py`, `exp_lifted_orbit.py`, `exp_ot_2d.py`, `exp_dither_ot.py`, `exp_uniform_basin.py`, `exp_entropic_leak.py`, `exp_ot_basin.py`, `exp_adaptive_prior.py`, `exp_adaptive_floor.py`, `exp_adaptive_floor_2d.py`, `exp_adaptive_floquet.py`, `exp_gauge_floor.py`, `exp_servo_cogging.py`, `exp_measurement_noise.py`, `exp_eiv_remedy.py`, `exp_cohesion_leak.py`, `exp_mass_transfer.py`, `exp_split_bifurcation.py`, `exp_barycentric_leak.py`,
  `exp_dither_stationary.py` and `exp_pe_general.py`. Their numbers are collected in `results/RESULTS_PARADIGM.md`.
- **Literature audit.** `LITERATURE_TRANSPORT.md` gives a novelty verdict for each claim. Reference numbers [1]–[89]
  point to `review/LITERATURE_refs.md`, and [S1]–[S51] point to the audit file.
- **Review log.** `review/REVIEW_LOG.md` lists every confirmed issue and its fix. The raw reports are in
  `review/REVIEW_round1.md`, `review/REVIEW_round2.md`, `review/REVIEW_round3.md`, for round 5 the five files
  `review/REVIEW_round5_*.md` (E0, I, J, DC, consistency), for round 6 `review/REVIEW_round6_E0prime.md`,
  `review/REVIEW_round6_Dpp.md` and `review/REVIEW_round6_Jprime.md`, for round 7 `review/REVIEW_round7_E0p3_E0dp.md` and for
  round 8 `review/REVIEW_round8_Dppp.md`, for round 9 `review/REVIEW_round9_D3b.md`, for round 10 `review/REVIEW_round10_D3c.md`, for
  round 11 `review/REVIEW_round11_D3d.md`, for round 12 `review/REVIEW_round12_D3e.md`, for round 13 `review/REVIEW_round13_D3f.md`, for
  round 14 `review/REVIEW_round14_D3g.md`, for round 15 `review/REVIEW_round15_D3h.md`, for round 16 `review/REVIEW_round16_D3i.md`, for
  round 17 `review/REVIEW_round17_D3j.md` and for round 18 `review/REVIEW_round18_D3k.md` and `review/REVIEW_round18_D3l.md`.

---

## 0. What this is, what is new, and what is not

**The idea.** An adaptive approximator with nonlinearly parametrised features can be read as an *empirical measure*
μ = (1/m)Σδ_{ϑ_i} on a lifted parameter space. Here the features are Gaussian RBF kernels with adaptive weights,
centres and widths.

- The adaptive law transports this measure, and the control output f̂ = ∫Φ dμ is linear in μ.
- Stability, robustness modifications, excitation, the many-kernel limit and the slow evolution of kernel geometry
  all acquire formulations on Wasserstein space. So does the global optimality of kernel placement.
- These formulations do not depend on how the kernels are labelled.

**Why the classical theory needs it.** The classical labelled analysis has two structural defects. This holds for
version 1 of the manuscript `lyapunov.pdf` (audited in `ASSESSMENT.md`) and for the literature it builds on.

- The classical Lyapunov analyses depend on labels. A network that is a relabelling of the ideal has zero function
  error but a large parameter error.
- There is no object that describes where the kernels go, or whether their placement is good.

A third apparent defect was that constants grow with the number of kernels. It turned out to be an artefact of the
normalisation. The mean-field normalisation removes it for labelled and label-free analyses alike (Remark A.1).

**What the closed loop actually is (new in version 3).**

- At infinite feedback gain, the averaged closed loop is a Wasserstein gradient flow. With dither it is a mean-field
  Langevin dynamics (Theorem H).
- At finite gain it is in general *not* a gradient flow of any free energy in the W_{G⁻¹} metric. This already holds
  for n = 1, so it is not a two-dimensional effect (Proposition H.7). The first correction has two parts at a frozen
  measure, and a third along the real loop.
  - The first part is a *semi-gradient* term. The adaptive law differentiates the loss on the closed-loop state
    distribution while ignoring that this distribution depends on the parameters, as in performative prediction. Like
    semi-gradient TD, the resulting update is not the gradient of any objective. Its performative defect is
    present for every n; for n = 1 it is non-reciprocal exactly when the residual is non-zero on the orbit
    (Theorem I(iii)), and it is non-reciprocal in the checked n = 2 instance.
  - The second part is an antisymmetric interaction between kernels (n ≥ 2). Its kernel is an integral of
    Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ')] over the region that the reference orbit encloses (Stokes).
  - Along the real loop, a third part comes from the parameters' oscillation within each period. It is a second-order
    averaging term (Remark H.7a).
- Together these predict a new closed-loop effect: traversing the reference orbit backwards changes the learned
  kernels by O(1/k_e), although the occupation measure is unchanged. At k_e = 40 the prediction matches to 3% in size
  and 11% as a vector (cosine 0.994). The vector error halves with each doubling of k_e (`exp_orbit_orientation.py`).
- The gradient structure survives as an input-to-state Lyapunov function (Proposition H.8). This reduces the global
  conjecture for the real loop to an averaging estimate.

**Contributions and their honest status.**

- "Proved" means proved here. "Formal" means a derivation without error control.
- "Checked" means a test evaluates the statement, or every inequality in it, along simulated trajectories or on
  sampled states.
- The novelty column follows the verdicts in `LITERATURE_TRANSPORT.md`.

| result | statement | status | novelty (audit) |
|---|---|---|---|
| A | The adaptive law is the particle form of the Wasserstein gradient flow of a μ-linear potential. Projected rates coincide only for a per-kernel projection (m ≥ 2). | proved and checked | interpretive; ingredients known [29], [41], [46], [25] |
| B | V_OT = ½‖e‖² + ½W² satisfies a Dini inequality with a two-regime matched remainder. It gives conditional RMS bounds along (L), and along (L_σ) with a leak term. | proved and checked, including a moving ideal, ε ≠ 0, switches, an active projection, and the failure of the plain bound along (L_σ) | packaging new; ingredients known [62], [64], [20], [49] |
| C | OT σ-modification: exact dissipation for selection solutions; a gap-protected strict inequality with √C ≤ b + 2η; a local exponential estimate, m-uniform for a single-atom prior and O(m^{−1/2}) otherwise; reshuffling by the labelled leak | proved and checked; the C(iv) instance satisfies every hypothesis | novel as a control construction; C(i) and the obstruction are known [71], [73] |
| C′ | The OT leak toward a prior with any number of atoms, through the barycentric projection of the optimal plan: dissipation of W(μ, ν) up to the within-cell variance of the plan, fixed points of the pure leak = equal-mass centroidal configurations of the prior (Lloyd's algorithm in continuous time, the G-gradient flow of the quantisation energy), induced-coupling bound √C_bar ≤ b + 2η, Lyapunov inequality and local exponential estimate with the bias b_bar := distance of the ideal from the centroids of its prior cells (≤ b) | proved for vertex-selection solutions (existence assumed) and for generalised solutions (which exist); checked (dissipation at every step, convergence to a Lloyd fixed point, the coupling bound and plan stability, the Lyapunov inequality along the closed loop) | the law is the continuous-time Lloyd descent of coverage control [S18] with equal-mass (capacity-constrained) OT plans [S19], used as a σ-modification leak; new here are the induced-coupling bound, the Lyapunov inequality with the bias b_bar and the transfer reading; the quantisation side is Lloyd/CVT theory [S20]; it reduces to Theorem C for m₀ = m |
| D | Local well-posedness of the mean-field loop; a Dobrushin estimate; a tube bound and global existence under a self-consistent growth condition; the mean-field limit on the existence interval | proved; consistency in m checked numerically (the constant is not checked) | technique anticipated [77]–[79] |
| E | Label-free local exponential convergence under tangential PE, in a uniform tube, at fixed m, with explicit radius, overshoot and rate (Lemma E.0); exponential merging of solutions. In the slow-adaptation regime a second strict Lyapunov function (Lemma E.0′: second-order slaving and a second-order averaging correction, with the peak factors of the regressor as data) has the averaged rate λ_min⟨ΨᵀΨ⟩/k up to a factor 1 − δ, and under parameter dither bounds the stationary parameter spread by (1 + O(ε))k‖Σ‖_F²/(2(1 − δ)λ_min) (Corollary E.0′.3); for non-periodic PE regressors the same construction with an exponential average over a window gives the rate (1 − δ − δ₂)λ_T/k, λ_T the weighted PE level of the window (Lemma E.0″) | proved, self-contained; linearisation checked; the strict Lyapunov inequalities checked pointwise (E.0′ by a certificate over all states) and the explicit envelope checked on the linearised and on the implemented loop; convergence observed. The rate of E.0 is conservative in a quantified way (Remark E.0a); the rate of E.0′ is within 12% of the exact rate above its threshold | partially anticipated [49], [54]; the strict Lyapunov function is of the Mazenc–Malisoff type |
| F | Along (L), the Wasserstein length of learning is ≤ λ_G^{1/2}L₁∫‖e‖; leak versions | proved and checked, including the leak clause, a tight leak instance with e ≡ 0 and a tight velocity instance | elementary |
| G | Averaged joint descent on the occupation loss. The slow geometry follows the gradient flow of the profile loss in the metric induced on the slow manifold, which is plain variable projection only when that manifold is flat (S ≪ 1). | formal; checked directly, including a counterexample to plain VarPro at S = O(1) | a fixed-m, closed-loop analogue of the two-timescale/VarPro training of [S1], [S2] |
| H | The averaged mean-field loop with saturated features, σ-modification and dither is mean-field Langevin dynamics. It converges globally and exponentially to the regularised optimum, with the log-Sobolev constant given in control quantities. | proved for (𝔄); the reduction to (𝔄) is formal (Remark H.0); real loop: Conjecture H.6 | a bridge to [S3], [S4]; new in adaptive control |
| H.7 | At a frozen measure, the exactly averaged loop is (𝔄) plus an explicit O(1/k_e) correction: a semi-gradient term, whose performative defect is present for every n, and an orbit-area circulation kernel (n ≥ 2). Along the real loop, a second-order averaging term is added (Remark H.7a). | Frozen-measure part proved. Stokes identity and parity hold to machine precision. The expansion is checked through its O(1/k_e) convergence. Non-reciprocity of the exact second variation is checked to 3% (asserted; 1–1.5% measured) and characterised by Theorem I. The full first-order theory predicts the orientation dependence of the learned kernels in the real closed loop to 3% in size and 11% as a vector (cosine 0.994, k_e = 40). | no precedent found (targeted search, round 3); ingredients known (higher-order averaging; performative prediction [S14]; semi-gradient TD [S10], [S11]) |
| I | Reciprocity at order k_e⁻² is a pairwise property: the antisymmetric kernel A_even + 2K is a coboundary iff it vanishes. For n = 1, A_μ ≡ 0 iff the residual vanishes on the orbit, and the performative non-reciprocity is bounded by the residual. For n ≥ 2 at perfect fit, A_μ ≡ 0 iff the oriented orbit carries no circulation (n = 2: zero winding function) and the residual's Jacobian is symmetric along the orbit; A_μ = 2K when the residual vanishes near the orbit. A_even is the curl of the residual paired with the features plus the residual paired with the Lie bracket of the feature fields. No open set of measures is reciprocal, for any n: the averaged loop is nowhere a gradient flow at order k_e⁻². | proved; every formula checked against the generic second-variation code and, by the round-5 referee, against the exact second variation; the characterisations checked on perfect-fit, curl and zero-current instances | new; the analytic tools (injectivity of Gaussian convolution; Stokes with winding numbers) are classical |
| D′ | Finite-horizon propagation of chaos for the dithered closed loop (L_τ): by synchronous coupling, E sup_{t≤T}(‖x^m − x̄‖ + (1/m)Σ‖ϑ_i − ϑ̄_i‖) ≤ T K_Φ e^{ΛT} m^{−1/2}, with Λ independent of m and τ, and a per-particle sup bound of the same order; the plant and the learned function converge at the rate m^{−1/2}, whose exponent is dimension-free, and only W₁(μ^m, μ̄) carries m^{−1/p}; under a dominant leak (small-gain condition κ_eκ_ϑ > L_b′L_ϑ) the bound is uniform in time (Theorem D″), a leak-dominated regime in which the network does not learn (at τ = 0 it collapses to one kernel; at τ > 0 it stays in the dither cloud of the prior) | proved; checked: the coupled distance to a 4096-particle reference decays like m^{−0.47±0.06} (24 seeds, τ = 0; −0.46 ± 0.05 with dither) and like (1/m − 1/M)^{0.94±0.08}, the law of the shared-particle coupling | the technique is Sznitman's, and the dimension-free exponent for the state is its standard output; new in that the plant is in the loop, with constants in control quantities |
| D‴ | Uniform-in-time propagation of chaos *in the learning regime* for the averaged slow loop (𝔄_m) at τ = 0 (Theorem J′'s leading order): under the synchronous coupling with the mean-field characteristic flow, the interaction term of the coupling identity is −‖f̂_Θ − f̂_Θ̄‖²_{L²(ρ)} + O(D²) — flat convexity of the loss makes it contract whatever the gains — and only the residual-weighted feature curvature R[μ](ϑ) = ⟨r_μ∇²_ϑΦ̃⟩_ρ (Proposition K's R₂) can expand; if along the mean-field flow σ₀I − G^{1/2}RG^{1/2} ⪰ κI on a neighbourhood of the support from a time s₁ on — the strict second-order optimality condition of the fitted configuration, compatible with σ₀ < σ_crit — then sup_s E(1/m)Σ‖ϑ_i − ϑ̄_i‖²_{G⁻¹} ≤ 4e^{2Λ₁}(L_ϑ^G)²K_Φ²/(κ²m) plus a tail exponentially small in m (uniformly in time under exponential convergence of the flow, over horizons e^{cm} without it), no particle changes basin, and the learned function converges at m^{−1/2} uniformly in time while the network learns; the window κ_R[μ_∞] < σ₀ < σ_crit is where learning and uniformity coexist; at finite m the stationary coupled distance follows the weight-sensitivity law E D_∞ = (1/m − 1/M)Σ_k w̄_k tr(S̃_kC̃S̃_kᵀ), S = ∂Θ*/∂w the sensitivity of the stationary atoms to the cluster weights (Proposition D‴.5, the sharp form of the theorem's constant); the leak's prior is a design input of the window: with a measure-valued prior (nearest-atom leak) placed on the plant's structure the stationary residual vanishes and the margin is maximal, κ = σ₀, while the floors at parameter and function level rise (1.6–2.2× on the instance, with the amplitudes), and the sign of the residual at a cluster sets its admissible weights — overshoot splits, undershoot stabilises, so a shrunk prior enlarges the admissible set (Proposition D‴.9, Remark D‴.11); the barycentric OT leak of Theorem C′ toward the same prior constrains the cluster masses to the prior's — the nearest-atom leak conserves them, and the floor is the memory of the initial cloud that this conservation law carries — so on the sizes commensurate with the masses the finite-m stationary measure is the mean-field one exactly, both floors vanish and the fit loses the weight bias, while on other sizes at most K₀ − 1 split particles bridge the clusters, with floors O(1/m²) at function level and O(1/m) in W² against the nearest-atom leak's 1/m and m^{−1/2} (Proposition D‴.12); in the real loop the collapsed-by-class subspace is invariant and carries the weighted K₀-atom loop, so with the OT leak the lift of one K₀-atom periodic orbit is a periodic orbit of every commensurate size — the same orbit measure, learned function and tracking error exactly, the real loop's finite-size floor being zero — while with the nearest-atom leak it is the lift of the orbit of the realised weights, whose distance between sizes is Corollary D‴.7's floor (Corollary D‴.14); the Floquet spectrum of a lifted orbit is the K₀-atom orbit's together with each class's ghost spectrum of multiplicity n_k − 1, independent of m (Proposition D‴.15); with dither the nearest-atom leak equalises the masses only by Arrhenius crossings of the Voronoi ridges (barrier ≈ 0.05 on the instance: frozen for 350/σ₀ at τ₀ = 0.005), while under the OT leak a crossing particle stays assigned to its own atom and returns, so the counts are balanced 88–99.6% of the time at the same crossing barrier (Remark D‴.16); every trajectory of the OT dynamics converges to a Clarke-critical point of its energy (Kurdyka–Łojasiewicz for the definable, semiconcave-plus-smooth J_m^{OT}), so the stationary statement is reached from every initial condition in the basin (Proposition D‴.17) | proved for the averaged slow loop (synchronous coupling, a bootstrap on the mean square and the per-particle distance, McDiarmid for the reference's fluctuation); the hypothesis is a computable property of the mean-field flow; it fails for a connected cold start whose mass splits between basins (the support threads the ridges), where the coupling carries a floor; the real loop at finite k_e: settled on a warm start by Corollary D‴.7 through Theorem J (both loops converge to their periodic orbits, whose distance is Proposition D‴.5's up to O(1/k_e)), open for a cold start (Remark D‴.4); dither: the entropy method of [S5] (imported); the theorem's regime m ≥ m₀ is asymptotic (m₀ ~ 10²³ on the instance) and Proposition D‴.5 is formal; numerics: on a three-bump instance at σ₀ = 0.15σ_crit the window is realised (κ_R[μ*](N_{0.1}) = 0.0051 < σ₀ = 0.0071 < σ_crit = 0.047) and (H) holds along a warm-start run; the coupled distance saturates (window maxima within 15% over 900 slow units) with slope −1.28 ± 0.30 in m, and the per-seed weight-sensitivity law predicts it within 25% in 9 of 12 cases with per-cluster cosines 0.92–1.00; the exact stationary distance of the replicated weighted configuration reproduces the measured one within 2.3% in all pairs, one m = 16 seed with inadmissible weights sitting on a saddle; a connected cold start grows and carries a floor of pairs in different basins whose fraction 0.20, 0.08, 0.04 (m = 16, 64, 256) decays like m^{−0.64±0.08}; the real loop at k_e = 25 follows its averaged slow dynamics on the same clouds to within 0.7%; for the OT leak the frozen-plan equivalence and the stationary statement are proved (total unimodularity, c-cyclical monotonicity, the Hessian of D‴.5), the quantisation bounds formal, and Theorem D‴'s coupling applies only to a balanced sample (probability ≈ 0.83/m under an i.i.d. start, and even there formal: Step 8 is written for i.i.d. particles) — for a generic sample the shared-particle coupling keeps the transit's mismatch for ever and a coupling theorem for the OT dynamics is open; numerics on the same clouds under the two leaks: cluster weights conserved exactly against equalised by s ≈ 100 slow units (ln 2/σ₀ = 98), function-level and W² distances at 10⁻⁴ and 10⁻³ of the nearest-atom leak's after 4.2/σ₀ and still falling on the commensurate sizes while the synchronous D sits at 0.95 of its mismatch plateau, the quantisation laws m·W²(μ_m, ν) ≈ 20 (slope −0.98 ± 0.03) and m·fit ≈ 0.55 (squared slope −1.81 ± 0.05) on the incommensurate ones, the reference's fit error 0.012 → 8·10⁻⁵; the real loop at k_e = 25 on the same clouds: learned functions of m = 15, 60, 510 within squared distances 7.5·10⁻⁹–2.7·10⁻⁸ of each other under the OT leak against 10⁻⁴–10⁻³ under the nearest-atom leak (the lifted orbits' distance to 1.0003 at the one pair with admissible weights; three of the six nearest-atom systems sit at saddle weights and have split), rms tracking error 0.9–7·10⁻⁶ against 3·10⁻⁴–1.9·10⁻³ (on a planar plant with four kernels: 10⁻⁷ against 10⁻² and 4·10⁻⁶–1.6·10⁻⁵ against 2–7·10⁻³); the lifted orbit's transverse spectrum equals the ghosts' and the collective block the K₀-atom loop's to 10⁻⁴ relative (invariance defects ≤ 10⁻⁷) | the two ingredients — the sign of the interaction from flat convexity, and the curvature at the actual rather than the worst-case residual — are new in propagation of chaos for adaptive control; the decomposition Gram + residual×curvature is the one of the interpolation literature ([S39]); the window between Proposition K's collapsed instability and the fitted configuration's stability is new; the finite-size floor as a conservation law of the leak, broken by the plan's mass constraint so that the particle system's stationary measure is exactly the mean-field one, is new; the lift of a weighted K₀-atom orbit as the exact periodic orbit of every commensurate size, with an m-independent Floquet spectrum (the orbit-level form of Proposition K's block structure), is new in adaptive control, where the size of the approximator is usually a nuisance parameter |
| D‴.18 | Theorem J uniformly in m at a replicated configuration: for weights with margins (collective block and separation blocks ≥ λ₀) the real loop's orbit near Θ*_m(w) is the lift of the weighted K₀-atom orbit, unique in the per-particle ball of radius ρ/2 (inside the |||·|||-ball of radius ρ where the period map contracts) and attracting every solution started in it at the rate (1 − κ)min{λ_coll, λ_sep}/k_e, with k₀, ρ, C independent of m; in the unnormalised norm of Theorem J the basin is m-independent too but shrinks per particle like m^{−1/2} and the threshold k₀ grows like √m, and for the nearest-atom leak no mean-square (W₂-type) basin is possible — one particle in another cell changes the weights (Remark D‴.19) | proved as an adaptation of Theorem J's proof in the per-particle norm (three structural facts: the plant sees the class means up to second order in the separations; each particle's drift depends on the others only through the plant; the linearisation is block diagonal with an m-independent spectrum); numerics: the Floquet exponents within 0.56–3.22 % of the averaged prediction and the ghosts' to 3·10⁻⁴; on the matching prior fitted rates 0.95–1.05 of the prediction for m ≥ 15 and constants 1.03–1.29, on a moving orbit 0.96–1.16 and 1.04–1.23, all from radii 31–261 times the proved one; across a cell boundary the nearest-atom leak ends at other weights (two of them inadmissible) while the OT leak restores the masses q | uniform-in-m averaging for an adaptive law with a mean-field interaction: the size of the approximator drops out of the stability constants — new in adaptive control, where Lyapunov and averaging constants usually grow with the number of parameters |
| D‴.20 | The entropic OT leak (Sinkhorn plan, regularisation ε) is a smooth gradient flow of mL̄ + (σ₀m/2)OT_ε that interpolates between the OT leak (ε → 0) and the σ-modification toward the barycentre of the prior (ε → ∞); at a collapsed configuration the separation block is σ₀(I − (2/ε)Cov^G_k) − R^G_k with Cov^G_k the conditional covariance of the atoms, so where A_k := σ₀I − R^G_k ≻ 0 an over-weighted class is a saddle exactly below the critical regularisation ε_c = 2σ₀λ_max(A_k^{−1/2}Cov^G_kA_k^{−1/2}) (implicit, along the branch Θ*_ε(w)), and a split moves its excess toward an atom its conditional plan charges, which need not reduce it; at fixed ε the masses are equalised only down to a critical excess δ_c(ε) ≈ w_kp_c(ε), and at large m the floor returns as Proposition D‴.5's law at the entropic stationary point | (i)–(iii) proved, the separation block checked entrywise against the finite-difference Hessian; (iv) proved under A_k ≻ 0 (which fails at large excesses and small ε), the conditional plan charging both under-weighted atoms (Cov^G_k of rank two, the rank-one formula 1.22–1.31 times ε_c); (v) formal; numerics: ε_c = 5.1, 13.5, 23.1, 65.0 for excesses 0.02–0.2; at ε = 1 every size equalises; at ε = 4 every system ends, after up to 6000 slow units, at stable weights with residual excesses up to 0.017, of the order of δ_c(4) = 0.012–0.016; at ε = 16 and 64 the stable systems keep their sampled weights, with floors 0.23–0.91 of the nearest-atom leak's | the critical temperature of deterministic annealing is the precedent ([S46]), with the same roles (particles as codevectors, atoms as data), and the ε → 0/∞ interpolation is that of the Sinkhorn divergences ([S45]); new are its role as the threshold of mass equalisation in an adaptive law, where the conditional covariance of a collapsed class is an effect of its excess, and the consequence that on the commensurate sizes only the hard constraint is floor-free for every m |
| D‴.21 | A measure-level basin for the OT leak: at commensurate sizes two exact identities along the plan segments — the energy, and its gradient paired with the displacement toward the stationary atom of one's own class — give a Polyak–Łojasiewicz inequality ‖∇J‖² ≥ λ·m(𝒥 − 𝒥*) with λ set by the residual-weighted curvature along the segments and independent of m (sharp, 2σ₀, at a matching prior); the excess energy decays exponentially from an m-independent sublevel set (tiny on the instance: none of its configurations has a particle in another cell below m ≈ 2·10⁵), and a posteriori wherever the plan-segment margin is positive, which covers configurations with a fixed fraction of the particles in other cells; under a curvature hypothesis at the stationary residual the replicated configuration is the unique global minimum among the empirical measures for every m, while the nearest-atom energy has a spurious minimum at every admissible count vector (the energy-landscape form of the floor) | (i), (ii), (iv), (v), (vi) proved, (iii) proved for vertex-selection solutions (their existence assumed); the mean-field form and the transfer to the real loop formal (Remark D‴.22); numerics: the identities to 6·10⁻¹³ along 64 runs, the margin positive from s ≤ 300, the observed rate at least λ throughout and λ within 1.5 % of the observed 2σ₀ at the end (matching prior); the nearest-atom runs at admissible counts stop at the energy of the replicated configuration of their sampled weights, at inadmissible ones (5 of 8 per prior at m ≤ 60) they split below it | the Polyak–Łojasiewicz inequality ([S47]), the generalised geodesics of AGS §9.2 and PL inequalities through optimal maps for functionals that are not geodesically convex ([S48]) are known; new is the PL inequality along the generalised geodesic based at the prior, with the loss linearised at the stationary residual and the counts making the constant independent of m |
| D‴.23 | Adaptive priors: the OT leak toward a prior whose atoms move toward the class barycentres at a rate η > 0 is the gradient flow of L̄ + (σ₀/2)W²(μ_Θ, ν_Z) in a weighted metric; its replicated stationary points carry no leak bias (the critical points of the unregularised K₀-atom loss), its non-replicated ones the K₀-cluster bias, and boundedness becomes a hypothesis (the leak no longer confines the parameters); at an exact representation it linearises, for every m, into the separations (rate σ₀) and 2 × 2 collective pairs whose slower rates satisfy ηh_j/(h_j + σ₀ + η) < λ_j^− < min{h_j, η} (h_j, η > 0), h_j the measure-space excitation spectrum — a learned prior removes the fixed prior's bias σ₀(H_L + σ₀I)⁻¹x⁰ mode by mode, at the rate of the excitation; as η → ∞ a capacity-constrained K₀-means regulariser; the transit drags the prior, so the adaptation should start from the fixed-prior stationary point | (i)–(iii) proved ((i) for bounded solutions), (iv) formal, (v) sketched (the radius stated, its constants not derived), (vi) the monotonicity of 𝒥 proved and the drag numerical; numerics: the K₀-atom modes decay at λ_j^−; from separated starts the m-particle systems depart from the K₀-atom run at second order in the separations, consistent with the same expected size at m = 15 and 510 (five draws); started from the fixed-prior stationary point every size m = 15–510 follows the K₀-atom run (per-seed ratios 0.995–1.000 at η = σ₀, 0.968–1.001 at η = ∞), with fit errors 1.6·10⁻⁴ and 3.4–3.5·10⁻⁵ against the fixed prior's 9.78·10⁻³ | proximal point / iterated Tikhonov regularisation (its per-mode filter factors), elastic averaging (its quadratic analysis) and soft weight-sharing are the precedents ([S49]); new are their transfer to the mass-constrained particle system uniformly in m (H_L as the measure-space excitation spectrum, the separations decoupled at σ₀), and the transit's drag with its gate; the limits (the lost a-priori bound, the K₀-cluster bias) are the precedents' own, specific here in their particle form and in the run-off from an exact fit |
| D‴.24 | The finite-size floor is the leak's bias at the sampled weights: under the nearest-atom leak toward a fixed matching prior the function-level weight-sensitivity is F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰, the weight-derivative of the bias, and the floor law decomposes on the excitation spectrum with the band-pass weight h_jσ₀²/(h_j + σ₀)² ≤ σ₀/4, centred at h = σ₀ — at the samples where the fixed prior's replicated point is stable (D‴.5-admissible); the nearest-atom leak toward an adaptive prior is the OT leak with the masses minimised out (the cell masses), and its stationary replicated configurations are the critical points of the unregularised weighted loss, the exact representation of the sampled weights among them — no function-level floor at any representable sample, every separation relaxing at σ₀ (no splitting), perfect tracking an invariant set of the real loop — while the coupled distance of the amplitude–mass gauge grows (the sample's weight fluctuation moved from the bias into the parameters) | (i) F exact at q, the bias formula first order, the law (B) formal (the delta method, as D‴.9), (ii) proved inside the cells and across crossings for vertex-selection solutions, (iii) proved for the stationary points, (iv) the boundary case formal, (vi) the invariance proved; numerics: at the stable samples the exact bias at the sampled weights reproduces the observed floor (ratio 1.0000), while 5 of the 12 samples are fixed-prior saddles and split; the law is asymptotic (the exact floor's expectation 2.8, 1.5, 1.1 times it at m = 15, 60, 240); under the adaptive prior the function-level distance decreases without a floor (at s = 6000 at most 2.0·10⁻⁶ at η = σ₀ over the representable samples) and the coupled distance is 1.7–4.8 times the fixed prior's, while W² is the cross-class transport of order ‖δw‖₁ under both priors; the samples outside W have probability 0.36, 0.074, 0.0021 at m = 15, 60, 240; in the real loop the tracking error follows the averaged prediction (ratios 0.987–1.001) | the 1/m fluctuations of mean-field networks and their vanishing at an interpolating limit ([30], [S51]), the delta method and the filter-factor bound of ridge regression are classical; new are the floor as the bias's weight-slope, mode by mode, a leak toward a learned prior as the way to restore the interpolation property under regularisation, and the design table of Remark D‴.25 |
| D‴.26 | Perfect tracking with an adaptive prior in the real loop: the lift of the exact representation of representable sampled weights, with the prior at its atoms and the plant on the reference, is a T₀-periodic solution whose Floquet exponents are σ₀/k_e exactly on the (m − K₀)p separations — every separation stable, also where the fixed prior's replicated orbit is a saddle — and the weighted K₀-atom loop's on the collapsed subspace for every m; where ∂_xf̂ = ∂_xf on the orbit (an exact representation) these are (1/k_e)[λ_j^±(w) + c_j^±(w)/k_e² + O(k_e⁻⁴)] — in slow time the odd orders vanish (the variational equation's symmetry with the slow block at −ε, the time reversal of the reversed reference's: the slow generator's spectrum is even; the first-order correction is skew), and c_j, explicit, contains the error filter's attenuation (≤ 0), the slaving's part (≥ 0) and averaging terms, its sign not fixed; so perfect tracking is exponentially stable (locally, for each m) at the rate (1/k_e)min{σ₀, λ_1^−(w) + c_1^−(w)/k_e² + …} | (a), (b) proved, exactly; (c) proved for x_m ∈ C⁴ at simple eigenvalues, the stability where H_L(w) ≻ 0 (Lemma D‴.26a: the slaving exact, the O(k_e⁻⁴) remainder bounded, the parity an identity), c_j explicit; the basin uniform in m not proved; numerics: at (6, 5, 4)/15, m = 15 and 30, k_e = 25 and 50, η = σ₀ and ∞: the separations at σ₀/k_e to ≤ 1.7·10⁻⁶, the collective exponents the K₀-atom loop's to ≤ 6.4·10⁻⁴; by the variational equation to k_e = 400 on two references (one with a non-zero first-order correction): the first-order deviations falling by → 4 per doubling (every exponent below its prediction on the instance; 13.6% in the slowest mode at k_e = 25), the second-order residuals by → 16, c₁/λ₁ = −100.6 of which the filter's −102.1 = −0.989ω₁²; with the gain scaled by 30 two c_j > 0, those exponents above their prediction by c_j/(λ_jk_e²) | Floquet theory, averaging and the Magnus expansion's commutator term, the averaging analysis of adaptive loops (the PE level), and leaks compatible with exact tracking (e-modification; a centre that follows the parameters) are classical, and (a)–(b) are Corollary D‴.14(ii) and Proposition D‴.15 at zero residual; new are the zero-residual lift at every representable sample with its exactly m-independent spectrum, and the explicit second-order coefficient of this loop (the even expansion is the classical consequence of time symmetry, [S56]) |
| D‴.26a | The exact slow subspace of the perfect-tracking linearisation: Chang's periodic decoupling P (a contraction) makes the slaving of D‴.26 an exact block-triangularisation, so the slow Floquet multipliers are exactly the slow block's; its generator log Φ(ε)/(εT₀) has an expansion to O(ε⁴) with a bounded remainder (the Magnus series convergent for εT₀‖A‖ < π); the slow block at −ε, defined by the same contraction, is the time reversal of the reversed reference's at ε, so with the variational equation's symmetry the slow generator's whole spectrum is even, eig Ā(−ε) = eig Ā(ε), and the rates are exactly even; and Chang's second transformation decouples the slow coordinates' drift from the fast ones exactly, their noise εbν, b = −(K̇ + aK)ᵀ − ε(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ + O(ε²), shared with the fast block — at leading order the high-passed −ε(K̇ + aK)ᵀν —, which gives D‴.30(ii)'s floor as a period average | proved (the constants not explicit; ε₀ also set by the spectral gap and Ā₀'s margin): D‴.26(c) with a bounded O(k_e⁻⁴) remainder, D‴.30(ii) for the linearised loop under white noise as a period average with the relative remainder O(1/k_e) (at η = 0 with K = Ĵ, Λ = σ₀I); numerics (η = 0, σ₀, ∞; k_e = 12.5–1600): the invariance residual ≤ 3.8·10⁻¹⁵, P − P₂ = O(ε³) and P − P₍₃₎ = O(ε⁴), the slow block's rates the full equation's to 1.0·10⁻⁸, the rates at ±ε equal to 5.9·10⁻¹² of the largest (1.1·10⁻⁹ per mode) while the slowest mode's first-order deviation reaches 29.6%, a mismatch breaking it at O(1/k_e); b − b₀ − εb₁ = O(ε²); the weakest mode's covariance deviation the slow noise's attenuation 2.4–2.5 times the rate's correction, of opposite sign; the floor's period average at the formula to 3·10⁻⁵ at k_e = 1600, pointwise ±1% (the jitter) | Chang's decoupling — with its time-varying and stochastic versions (Kokotović–Khalil–O'Reilly Ch. 4–5 and App. A, Haddad–Kokotović, Yang–Zhu; Berglund–Gentz's slow manifolds under noise) —, the Magnus series' convergence and the reversal argument (the adjoint relation of symmetric methods) are classical ([S59], [S56]); new is their application here: the exact slaving, the parity as an identity at negative ε, the slow noise identified exactly |
| D‴.27 | A projected adaptive prior restores the leak's a-priori bound: with the prior's atoms started in and confined to a compact convex set 𝒵 (the projected flow; at η = ∞ the projected class barycentres), every particle's distance from 𝒵 is at most B/σ₀ + (d(0) − B/σ₀)₊e^{−σ₀s}, B the worst-case learning force, for every m, η and target; 𝒥 is non-increasing (at η = ∞ in its selected form, along a Lloyd-consistent plan; the reduced energy need not decrease); inside 𝒵 nothing changes (Theorem D‴.23's bias removal and linearisation), and where the limit has an atom on the boundary — necessarily when no stationary point has its prior inside — a bounded bias remains; at η = ∞ on collapsed classes the leak is the gradient of the exterior penalty (σ₀/2)dist²(ϑ, 𝒵), a switching σ-modification toward 𝒵 | (i), (iii) proved ((iii)'s local convergence as Theorem D‴.23(v), sketched); (ii) the monotonicity proved, the convergence argued at η = ∞ and formal at finite η; (iv) the identification at η = ∞ proved, the boundary atoms and the bias outside numerical; numerics: beyond the saturation both η converge to the same bounded point (fit error 0.1415; the projected midpoint rule consistent, the trajectory error 6.3·10⁻⁷–3.8·10⁻⁸ at dt = 0.25–0.0625), particles within 0.92 of 𝒵 (bound ≤ 3686); the box's amplitude bound 0.75–0.99 of a_max trading the fit error 0.1676–0.1101 between the fixed prior's 0.1750 and the unprojected limit 0.1050, tighter boxes (which exclude the exact atoms) worse than the fixed prior; at η = ∞ the reduced energy rising by 26% along a Lloyd-consistent plan while the selected energy falls; on the four-bump target the box stops the width's run-off at its face; at the misplaced prior the box never approached | parameter projection, the exterior quadratic penalty and the switching σ-modification (which keeps the ideal properties inside its set) are classical, and at η = ∞ on collapsed classes the construction is the latter; new is its transfer to the learned centre of a mass-constrained K₀-atom prior at finite η, compatible with Theorem D‴.23's bias removal inside 𝒵 |
| D‴.28 | The adaptive prior's floor removal is a gauge property: the masses enter the network only through w_ka(α_k), so for any nondegenerate critical point ϑ^b(q) of the unregularised loss — exact or not — its gauge transform G_wϑ^b(q) is a critical point at every sampled weight vector of W_g with the same function: at a strict local minimum no function-level floor at the samples of W_g with positive margins, whatever the approximation error; the fixed matching prior's floor is σ₀Ĵ(Ã + σ₀I)⁻¹g̃ with the full Hessian Ã = H_L − R̃; for a scalar amplitude per kernel the residual curvature has a zero amplitude row at a critical point and scales exactly by q_k/w_k under the gauge, so the adaptive prior's splitting margins are σ₀ − (q_k/w_k)λ⁺_k | (i), (iii) proved exactly, about the stationary points; (ii) as D‴.24(i); the dynamics' convergence, and with it the floor's vanishing, numerical; numerics: with a fourth bump that three kernels cannot represent (fit errors 0.4–2.5%), the critical point continued in w from ϑ^b(q) is the gauge transform to 2·10⁻¹¹ with f̂ unchanged to 6·10⁻¹³, the adaptive prior's sensitivity 10⁻¹¹–3·10⁻¹⁰ of the fixed prior's, the fixed prior's equal to the full-Hessian formula to 2·10⁻⁹ (the Gauss–Newton form off by 7·10⁻⁴–5.6·10⁻³), the margin formula exact to 10⁻¹⁴ and the dynamics' separation eigenvalues to 6·10⁻¹²; P(w^m ∈ W_g) = 0.52 at m = 15; the dynamics: under the fixed prior the floor stays and the fit error above the approximant's; under the adaptive priors the floor falls by orders of magnitude and the fit error tends to the approximant's at the samples inside W_g; outside it the systems approach the box minimum on the branch (1.11–1.59 times the approximant's fit error) | the homogeneity of the conic view ([49], [55]) and the rescaling symmetries of networks ([S55]) are classical; new are the gauge as the mechanism of the floor removal, robust to approximation error, and the exact margin formula for a scalar amplitude |
| D‴.29 | Second-order plants: with the sliding variable s = ė + λe and u = ẍ_m − λė − k_s s − f̂_Θ(z), the error channel is ṡ = −k_s s + r_Θ(z), ė = −λe + s — Theorem J's with s in place of e — so as k_s → ∞ the parameters follow the averaged slow loop along the regressor's orbit z_m(t), every averaged-loop statement of §5d transfers with the occupation measure of z_m, and the position error is the residual filtered by (d/dt + λ)⁻¹(d/dt + k_s)⁻¹; kernels on the coordinate f depends on see the orbit as an interval, while isotropic phase-plane kernels on a circular orbit (ω₀ = 1 in isotropic units) are over-parametrised (von Mises restrictions: in general position one zero Gauss–Newton eigenvalue per kernel; on other orbits weak excitation) | (i), (iv) exact; (ii), (iii) formal (first-order averaging with the fast pair (s, e); Theorem J's steps carry over for the smooth fixed-prior leak once the error-channel steps are redone; switching targets beyond the first-order theory); numerics: a position servo with a cogging torque whose four-kernel local minimum ϑ^b(q) (fit error 1.46%; the loss's infimum lies at the saturation, ≤ 0.64%) is Proposition D‴.28's regime: from clouds around ϑ^b(q) the fixed prior's floor 1.85·10⁻³–4.49·10⁻⁵ (m = 16–128), the adaptive prior's 38–780 times lower, the OT leak at the approximant, σ-modification toward 0 biased (0.328), and no leak as good (0.0180–0.0146); from an engineer's model (fit error 0.200) the fixed designs stay biased (0.048–0.10) and the adaptive prior removes most of the bias, slowly (0.022–0.027 at s = 2000; the K₀-atom flow 0.0155 at s = 10⁵); the real servo's position error at the exactly filtered prediction to 1.0001–1.0029 (the slaved one 0.931–0.996), no design drifting; noise-free, so the leaks' robustness is not exercised | the filtered error and Gaussian networks in adaptive control are classical (Sanner–Slotine), and so is the partial excitation of Gaussian networks along an orbit; new are the transfer of the transport leaks, their floors and the adaptive and projected priors to second-order plants, and the exact fibre of isotropic phase-plane kernels on a circular orbit |
| D‴.30 | Measurement noise at perfect tracking: the noise is coherent — a replicated start stays replicated under every noise path, so the floor does not average out with m (invariance, not attraction: near the lift the separations grow above a noise level, which drifts along the trajectory with the developing errors-in-variables deviation); to first order the loop high-passes the measured error's noise, ((d/dt − a)/(d/dt + k_e))ν, and the function-level floor is (N₀/2k_e)[Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) + tr H²], Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩ — the prior's mobility enters the Tikhonov filter factor as h_j + η, and the adaptive prior's floor is (N₀/2k_e)tr Q whatever the spectrum; to second order the noise inside the regressor exerts a mean force (the errors in variables) whose stationary response is (H + σ₀I)⁻¹ under the fixed prior and H⁻¹ under the adaptive prior, for every η > 0 | (i) exact as invariance, its threshold formal at the lift, the long-time threshold open; (ii) proved for the linearised loop under white noise as a period average (Lemma D‴.26a), (iii) formal, their closed forms exact for the averaged system; numerics: the floor per unit N₀ 0.1597 (fixed) to 0.1757 (η = ∞), the weakest mode's variance ×24.1 (×28.0 in the exact linearised loop at k_e = 25); the exact periodic covariance of the linearised loop at the prediction to 0.976–1.00001 (k_e = 25–800: the relative error O(1/k_e), its first-order term from the slope a, the slope-free control numerically even in 1/k_e); the errors-in-variables bias 491 (fixed) and 2397 (adaptive) per (Eν²)² averaged, 504 and 2653 in the exact linearised loop, the parameters' fluctuation negligible; the real loop (12 common noise paths, two noise levels) within 1.75 standard errors of the prediction, the mean part 89% at std 0.01; the lift m = 15 identical to the K₀-atom run to rounding; the separations' threshold std 0.031 at the lift; on antithetic pairs the noise-induced separation rate λ_eff = 0.327, 0.318, 0.294 at std 0.01, 0.03, 0.05 against the lift's 0.298 and the transient's 0.326–0.328, the small-noise limit 0.329 ± 0.003, the rate vanishing at std 0.0299 over 300 time units | measurement noise, the errors-in-variables bias, the diffusion approximation, the LMS misadjustment, leaky LMS, common noise and the lack-of-excitation drift are classical; new are the floor's interpolation in the prior's mobility (h_j + η in the filter factor), the driving matrix Q of the sensitivity-filtered measurement, and their combination with the coherence of the K₀-atom identification and its noise threshold at the lift |
| D‴.31 | The adaptive prior is a washout leak: the prior lags the class barycentre by its exponentially forgotten increments, so the leak acts in full on the separations and, for integral plans, on the barycentre through the high-pass iω/(iω + η) — mode j's stiffness h_j + σ₀iω/(iω + η) — which is its bias removal and its lost robustness to persistent forces on the barycentres; along every affine time-invariant prior the stationary response R and the bias B toward a misplaced anchor satisfy RH + B = I (per mode h_jr_jj + b_jj = 1), so no prior gives both; a force's function-level response is at most the norm of its smallest equivalent target perturbation, under every leak with S + Sᵀ ⪰ 0; and the desired regressor (the adaptive law's regressor at the reference) keeps the lift, the linearisation and the first-order noise response while turning the noise's own persistent force into the least-squares force of g = k_eE[e], its response at most ¼‖∂²_xf‖²_ρ(Eν²)² uniformly in k_e, B and the spectrum — the parameters' fluctuation adding a force of the same order that is not the learning force of a mean tracking error, so no bound in the noise's moments applies | (i)–(iii) proved; (iv)(a)–(b) proved at the lift and for the linearised equations, (iv)(c) formal, its noise part's consequences exact for the averaged system, the parameters' part numerical; numerics: the adaptive prior's bias 2397, 0.471, 0.0526 per (Eν²)² (noise part; MR, DR, DR*) at k_e = 25, with the parameters' fluctuation 2398, 0.486, 0.0481 (DR +56% at k_e = B = 100, where its noise part nearly cancels); MR over DR 3.0·10³–2.1·10⁸ over k_e = 25–200; R and B from the priors' dynamic blocks, RH + B = I to 1.4·10⁻¹⁴, also for non-commuting leaks, where HR + B − I reaches 0.47–0.79; the real loop on common noise paths: E‖f̂ − f‖² 8.7 times lower under DR (predicted 9.2), the MR − DR mean difference at its prediction to 0.5–5.3% (0.4–1.8 sampling errors) | washout filters, low-frequency learning, the desired regressor, the Tikhonov/Kalman split and the minimum-norm bound are classical; new are their transfer to the transport leaks (the washout reading of the adaptive prior, uniform in m) and the quantification; the parameters' own errors in variables is the classical second-order bias of a constant-gain stochastic approximation ([S57]: Benveniste–Métivier–Priouret, Dieuleveut–Durmus–Bach, Box), new only the observation that the desired regressor does not remove it |
| D‴.32 | The cohesion leak: the particles attract each other with a finite-range or Gaussian kernel — the leak is σ₀ × the kernel mass × the mean-shift vector of the network's own parameter distribution, no prior, no atom count, no plan — and the averaged loop is the gradient flow of mL̄ + (σ₀m/2)∫∫W dμdμ (in the mean field an aggregation equation in the fit's external potential with two interactions, the fit's and the cohesion's), its empirical measure the mean-field solution exactly; with the finite range, at a clustered configuration whose atoms are R apart the leak vanishes, so every critical point of the unregularised weighted loss is stationary, for every number of classes and every set of masses, and it acts exactly on the unexcited subspace: the separations of class k at σ₀w_k, the class means the unregularised flow's; a class can split iff σ₀w_k² < q_kλ⁺_k at a gauge transform, so σ₀ sets the resolution and an exact representation is never over-split, and a split sheds single particles (the unstable separations degenerate, the kernel concave); the Gaussian kernel pays an exponentially small bias, and twins with no third class within reach contract at σ₀ψ(Δα²/(G_αR²))(w_ka′_k + w_la′_l)²/(w_ka′_k² + w_la′_l²) — a third class within reach can win; in the real loop the lift is periodic with the separation exponents σ₀w_k/k_e and the η = ∞ adaptive prior's collective spectrum | (i)–(iii) proved; (iv) proved for the stationary points and their linearisation, the split's pattern explained; (v) the twins' rate to leading order, their merger only without a third class within reach (a counterexample), the dead class's motion and time scale stated, its merger not proved; (vi) the lift exact, the collapsed exponents by Corollary D‴.26 and Lemma D‴.26a; not proved: the transit, and with it any selection of the number of kernels — the count is the transit's and the schedule's; numerics: at the 15 count vectors of m = 15 the lift of ϑ^e(w) stationary (2.8·10⁻¹⁷) with (iii)'s spectrum to 2.8·10⁻¹¹ (the finite differences' error), the margins the Jacobian's to 2·10⁻¹¹; the twins' rate to 0.04%, opposite-sign twins with a third class within reach not merging (4·10⁵ slow units, a biased four-class fit); the floor removed from warm starts as by the adaptive prior; from broad starts the annealed Gaussian kernel at the target's three classes in 9 of 9 runs (≤ 2.3·10⁻⁶), at a fixed scale not (4–6 classes at R = 1, three biased at R = 2); the real loop's 45 slow Floquet exponents to 5.5·10⁻⁶; the dither's within-class spread 0.063–0.068 against √((1 − 1/n_k)τ/(σ₀w_k)) = 0.0635; σ₀ lowered: no class with a positive margin split and all 8 below −5·10⁻⁵ did, by shedding single particles (8–10 classes at σ₀/64 for four bumps); five and three heavy bumps biased (the transit's masses), two bumps closer than the broad scale merged into one | mean shift, aggregation equations, bounded-confidence and polarized consensus-based clustering, clustering in self-attention dynamics, Fukumizu–Amari's lines of critical points, splitting steepest descent, convex clustering and sum-of-norms identification, deterministic annealing and self-organising RBF networks are classical ([S60]); new is their combination as a leak in mean-field adaptive control: the exact blocks at clustered configurations with the mass factor, the gauge-scaled splitting margin, the twins' rate, the real loop's exponents |
| D‴.33 | The desired regressor keeps the replicated configuration attracting under measurement noise: with the adaptive law's regressor at the reference the regressor's own term ⟨E[ε_m²]G^{1/2}∂_x∇²_ϑΦ̃G^{1/2}⟩ — which carries the measured regressor's drift and does not decrease with k_e — leaves the separations' second-order drift; what remains, the mean error's term ⟨E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}⟩ and, once the parameters fluctuate, the cross term, is O(1/k_e) like the separations' stiffness, so the noise threshold is bounded below uniformly in k_e and constant when the noise's bandwidth scales with the gain (at a fixed bandwidth DR's sign reverses at k_e ≈ B); at the adaptive prior's stationary point the mean error carries the parameters' part of the force, not only the residual (I − P)g; under DR the compensation's rectification of the noise makes the drift at the lift negative in the centre directions of separated kernels, the whole matrix on the instance — not for same-sign, wide or edge kernels —, and along the trajectory the drift turns indefinite | (i) and (iv)'s formula formal at D‴.30(i)'s order, (ii)'s bound proved, (iii)'s balance formal with the exact periodic mean, (iv)'s sign proved in the centre directions of separated Gaussian kernels at the lift where the saturation's term is dominated, the whole matrix's numerical; numerics: at the lift DR's drift negative definite, DR*'s largest eigenvalue 60 times below MR's (the threshold std 0.238 against 0.031); k_e·λ_max constant for DR and DR* at B = 4k_e while MR's threshold falls like k_e^{−1/2}, at B = 100 DR positive definite from k_e = 200, the floor std ≥ 0.069; at the stationary point rms E[e] 6.7·10⁻⁴ (DR), its parameters' part 6.3·10⁻⁴, the drift indefinite (≤ +3.7·10⁻⁵); along DR's transient the drift first positive at t = 2714 and at most +1.05·10⁻⁴ (the formal threshold std 1.64); the real loop on common noise paths: MR's separations growing, DR's decaying faster than σ₀/k_e at the transient's prediction (0.97–0.98 of its noise part, the lift's 0.81–0.82), along DR's weakest directions by 0.14–0.36% of σ₀/k_e, DR*'s slower | the desired compensation law's insensitivity to measurement noise and the synchronisation of identical systems by a common noise are classical ([S53], [S61]); new is the removal of the regressor's term in the K₀-atom identification — a threshold bounded below uniformly in k_e — and the instance's quantification (the sign at the lift, its loss along the trajectory, the stationary point) |
| D‴.34 | Mass transfer: the particles carry weights adapted by the replicator law dβ_i/ds = γβ_i(⟨rΦ̃_i⟩ − λ) (in real time the exponentiated-gradient law β̇_i = γβ_ie(Φ̃_i − f̂)) beside the transport law with the cohesion leak; the pair is the gradient flow of the fit plus the cohesion energy in the metric Euclidean on the positions and Fisher–Rao on the weights (gradient descent on (√β, Θ)), Σβ and the saturation's bound kept without a projection; at fixed positions the relative entropy of the weights is a Lyapunov function of the real loop (V̇ = −k_e|e|²), the negative-entropy instance of a mirror-descent adaptation law; at clustered configurations under the finite range the stationary points are Proposition D‴.32's for every vector of masses, under the Gaussian kernel the replicator selects the masses and moves the bias; the saturated transients — a class short of mass creeping to a_max with its residual — are left at the replicator's rate, so the transit's masses no longer decide whether a target is fitted: at fixed positions bounded trajectories reach the exact fit iff Σ_k|A_k|/a_max < 1 | (i), (ii), (iii)'s stationarity at clustered configurations under the finite range, (v)'s zero modes proved; (iv) proved for bounded trajectories (boundedness not proved); the Gaussian kernel's mass selection and the transit numerical; numerics: at fixed positions the exact fit reached with the replicator (3.0·10⁻¹⁴ at γ = 10⁻²) where the transport alone stalls at 0.388 (the least error at those masses 0.387), and not when Σ|A_k|/a_max = 1.09 (0.071, creeping to the ℓ¹ ball's 0.0684); from broad starts the five-bump target at five classes in 3 of 3 runs (6.0·10⁻⁴–1.4·10⁻³, levelling off at 2.8·10⁻⁴ near the Gaussian kernel's stationary point; 6–9 biased classes without), the three heavy bumps at three (≤ 1.6·10⁻⁷; 3–5 biased without), the instance unchanged; two bumps closer than the broad scale merged as without, and the shedding at σ₀/64 unchanged (the fitness vanishes at critical points); at the class level the Gaussian kernel's stationary fit error 3.85·10⁻⁴ with the masses held and 2.63·10⁻⁴ with them free, the compact kernel's exact; the real loop's V monotone with V̇ = −k_ee² to 1.1·10⁻¹⁶; at the lift 15 zero modes in (Θ, log β) — 12 redistributions, 3 gauges | Wasserstein–Fisher–Rao distances and gradient flows, conic particle gradient descent and the √β lift, neuron birth–death, the exponentiated-gradient law and EM's multiplicative weights, mirror-descent adaptation laws with Bregman Lyapunov functions and the replicator equation's Shahshahani gradient structure are classical ([S62]); new is the replicator as the mass transfer of the cohesion leak in mean-field adaptive control: the stationary set unchanged at clustered configurations under the finite range, the saturation trap removed, the representability criterion |
| D‴.35 | The split of a class under the cohesion leak is an S_n-symmetric transcritical bifurcation: on the class's separations along the splitting direction v the flow reduces to dξ_i/ds = μξ_i + (c₃/2)(ξ_i² − Σξ²/n), μ the negative margin, c₃ the residual's third derivative along v (the cohesion, even, enters at third order); every small split state (q particles against n − q) is a saddle with the eigenvalues ∓μn/(n − 2q) and −μ, the S_{n/2} × S_{n/2} states cubic, so the split is a jump; the particles keep their order along v (equivariance), and once its excursion exceeds 2μ/\|c₃\| the extreme particle of sign(c₃)ξ runs away on the side sign(c₃)v — alone when its lead is large, with a near-tied neighbour as a bound pair —, at a time the reduced equation gives; after a departure the margin is σ₀(w_k − k/m) − λ⁺(r′, ϑ′), the stiffness's drop and the curvature's gauge rise both favouring further departures, the residual's change stopping them; the mass-free leak (the normalised mean shift) has the thresholds λ⁺_k, a factor w_k lower, and delays the cascade without stopping it | (i)–(ii) proved to leading order in μ (the centre manifold, the S_n normal form), (iii) the order by equivariance, the escape time in the quadratic truncation, (iv) a margin statement, (v) the blocks proved; the cascade's extent numerical; numerics: the split states found by Newton within 0.3%–3.0% of the prediction at μ = 5·10⁻⁷, their eigenvalues converging to ±μn/(n − 2q) and −μ (errors falling fourfold with μ); the first excursion's particle the reduced quadratic equation's in 12 of 12 runs (the linear equation's in 6, the flipped one's in 2), its time to 1.1% (the linear equation's to 22.4%); bound pairs in 3 of 6 runs at σ₀/16; under Part G's schedule 8–10 classes (the cohesion leak), 6–7 with the replicator at γ = 0.1, 5–8 at γ = 1, 9, 8, 9 classes (the fit error 0.0028, the cohesion leak's end state) under the mass-free leak continued to σ₀/256 | S_N-equivariant bifurcation theory — the standard representation's quadratic equivariant, transcritical branches with no stable small branch, the jump of sympatric speciation, the permutation-symmetric critical points of networks — is classical, as are order preservation and the Riccati escape ([S63]); new is the identification for the cohesion leak: the quadratic coefficient as the residual's third derivative, the cohesion at third order, the margin after a departure, and the instance's quantification |
| D‴.36 | The gauge under dither: at an exact representation with one particle per class the fits through the configuration form a manifold M of dimension K − 1 parametrised by the masses (f̂ constant on it); with the replicator and the positions dithered the masses diffuse along M (Katzenberger's limit) with the covariance D = (TᵀgT)⁻¹TᵀgΣgT(TᵀgT)⁻¹ in the flow's metric, for two classes 2τT_α/(T_α + T_W)², vanishing as γ → 0 and at the saturation, where the amplitude's excursion is unbounded; the entropic leak (ε/γ)Σβ log β on the masses is a leak on the unexcited subspace: it keeps the exact representation at W = w stationary (no bias where w is representable) and makes the gauge an Ornstein–Uhlenbeck direction with stiffness κ = ε(TᵀgT)⁻¹T_W, its spread solving κC + Cκᵀ = D (∝ ε^{−1/2}), the dither's drift term shifting it by O(τ/ε) | (i) proved, (ii) the diffusion coefficient by Katzenberger's theorem (the drift term not computed), (iii) the stationarity, the stiffness and the linearised spread proved; numerics: the masses' diffusion measured over 64 seeds at 0.95, 1.36, 0.90 times the prediction (the displacement's variance; within two standard errors); without the leak the spread 0.05 by 2·10⁵ slow units and the first amplitude at |a|/a_max = 0.9994; with it the spreads 0.0044, 0.0032, 0.0034 at ε = 10⁻³ against the predicted 0.0044, 0.0031, 0.0033 (4 seeds; at ε = τ inconclusive) | Diffusions forced onto a manifold by a large drift (Katzenberger) and SGD's noise-driven motion along the manifold of minimisers (Blanc et al., Li–Wang–Arora) are classical ([S64]); new: the amplitude–mass gauge of Proposition D‴.28 as such a manifold for the replicator of D‴.34, the diffusion in the flow's metric with its closed form and its excursion into the saturation, and the entropic leak on the masses as the pin |
| J | At finite m, without dither, in the scaling σ = σ₀/k_e: every hyperbolic critical point of the m-particle regularised loss J_m has a unique nearby T₀-periodic orbit of the real closed loop, at distance O(1/k_e), with the same stability type; the parameter Floquet exponents are (1/k_e)λ(G_m^{1/2}∇²J_mG_m^{1/2}) + O(k_e⁻²); the basin of a minimum is uniform in k_e; the period mean is Θ* − J⁻¹(b_frozen + b_osc) + O(k_e⁻²) | proved (tube and slaving, periodic Riccati block-triangularisation, second-order averaging, Poincaré map); checked by Newton shooting for the exact orbit, its monodromy matrix and its period mean against the averaged Hessian and the first-order prediction, at k_e = 10–80 | the first step of §8 branch 1; the tools (two-time-scale averaging under attraction) are classical [Khalil Ch. 10–11; Sanders–Verhulst–Murdock Ch. 5; Sastry–Bodson Ch. 4] |
| K | At a collapsed critical point (all kernels equal) the Hessian of J_m splits by permutation symmetry into a separation block σ₀G⁻¹ − R₂, with R₂ the Hessian of the first-variation potential at the atom (independent of m), and a collective block; the collapsed configuration is a strict local minimum iff the atom minimises the single-kernel objective and σ₀ > σ_crit := λ_max(G^{1/2}R₂G^{1/2}). The orientation-odd k_e⁻² terms of the real loop (U_odd, b_osc and the third derivative of J_m against the first-order offset) shift the separation Floquet exponents of the collapsed orbit by O(k_e⁻²), identically for the (m − 1)-fold multiplet, and in a window of width O(1/k_e) around σ_crit the direction in which the reference is traversed decides whether the kernels separate | proved ((iii)–(iv) given the second-order normal form, §6e); checked by shooting on instance A: the orientation difference of the softest separation exponent scales like k_e^{−2.18} and matches the prediction from the averaged data to 3% at k_e = 80; at a leak gain computed from that prediction, one traversal leaves the collapsed configuration stable and the other makes it unstable (k_e = 10 and 20) | new as a design statement; the splitting is the isotypic decomposition of the permutation action, and the stability criterion is the local Chizat–Bach condition at a Dirac mass |
| J′ | At finite m with dither τ = τ₀/k_e and leak σ = σ₀/k_e: the closed loop has an invariant law (unique under a smoothness hypothesis) whose parameter marginal is the Gibbs measure exp(−J_m/τ₀) up to O(1/k_e), with the explicit first-order correction ε∫φ g dπ driven by the first-order drift b of Theorem J as a vector field on parameter space (g solves τ₀∇·(πG∇g) = ∇·(πb)); the stationary tracking error is (2/k_e²)π(L̄_m) + O(k_e⁻³); the stationary parameter process is not reversible, with a circulation ε[∫Θ∧j₁ + π⟨u∧u̇⟩] + O(ε²), where j₁ = π(b − τ₀G∇g) is the first-order stationary current, nonzero exactly when b is not a G-gradient (Proposition J′.2: the finite-m form of Theorem I's reciprocity defect plus a Hessian commutator from second-order averaging; never a gradient for m ≥ 3) | proved (the Poisson-equation expansion of the invariant measure with the fast error and phase eliminated; three standard theorems imported and named); checked on a single-kernel instance: the first-order prediction removes the O(1/k_e) discrepancy between the real loop's stationary means and the Gibbs measure, the orientation-odd part of the means (the dither version of Remark H.7a) matches its prediction, and the current machinery is checked on the first-order averaged SDE | the technique is the Pavliotis–Stuart / Mattingly–Stuart–Tretyakov expansion; new in its object (a hypoelliptic learning loop with a plant) and in the identification of the first-order non-equilibrium current with the reciprocity defect of Theorem I |
| H.8 | (𝔄) is input-to-state stable under drift perturbations, with gain 1/(τρ_τ) | proved, assuming the chain rule for the perturbed flow | elementary given [S3] |

**Where the new mathematics lies.** The new mathematics this viewpoint requires lies in the open problems of §8.

1. Infinite-horizon averaging for McKean–Vlasov closed loops. Propositions H.7 and H.8 reduce it to controlling a
   near-identity averaging transformation along the moving measure. Theorem J settles the deterministic finite-m case.
2. Uniform-in-time propagation of chaos with a plant in the loop.
3. The theory of the non-reciprocal, semi-gradient mean-field dynamics that the closed loop actually produces.
   Theorem I characterises its non-reciprocity (completely for n = 1 and at perfect fit; by an explicit bilinear form
   otherwise, and nowhere reciprocal on an open set); its long-time behaviour is open.
   Proposition K turns one of its effects into a design rule: near the critical leak, the orientation of the reference
   orbit decides whether a collapsed network separates.
4. Global convergence without dither.

---

## 1. Setting

**Plant, reference, control, error.**

    ẋ = f(x,t) + u,    e = x − x_m(t),    u = ẋ_m − k_e e − f̂(x).

**Lifted parameter.** ϑ = (a, c, s) ∈ ℝⁿ × ℝⁿ × ℝ_{>0}, and p := 2n+1. The feature is Φ(x; ϑ) := a φ(x; c, s) ∈ ℝⁿ,
where φ(x; c, s) = exp(−‖x−c‖²/2s²).

**Network as a measure.** f̂(x) = (1/m)Σ_i Φ(x; ϑ_i) = ∫Φ(x;ϑ) dμ(ϑ), with μ := (1/m)Σ_i δ_{ϑ_i}. The map μ ↦ f̂ is
linear.

**Metric.** G = diag(g_a I_n, g_c I_n, g_s) ≻ 0 is the gain matrix. Distances are measured in ‖v‖²_{G⁻¹} = vᵀG⁻¹v,
and λ_G := λ_max(G).

**Admissible box and projection.**

- **The box.** Θ_K := {‖a‖ ≤ a_max} × {‖c‖ ≤ c_max} × [s_min, s_max].
- **Layer and inner set.** ε_K is the width of the boundary layer. Shrinking each bound of Θ_K by ε_K gives the inner
  set Θ_K^{ε}, with interior int Θ_K^{ε}.
- **The operator.** Proj_K(ϑ, y) is the per-particle, factor-wise projection (`projection.py`). It uses Pomet–Praly
  norm-ball operators for a and c, with a quadratic ramp h ∈ [0,1], and a box operator for s, with a linear ramp.
- **Regularity.** The operator is locally Lipschitz in (ϑ, y). The norm-ball factors are C¹ at the inner edge of
  their layers. The box factor, and every factor at the outer edge, is only Lipschitz. The proofs use only Lipschitz
  continuity and ramp values in [0,1].

**Property (P).** For ϑ ∈ Θ_K, ϑ' ∈ Θ_K^{ε} and all y,

    ⟨ϑ − ϑ', G⁻¹ Proj_K(ϑ, G y)⟩ ≤ ⟨ϑ − ϑ', y⟩.

On ∂Θ_K the projected velocity has no outward component, so Θ_K is forward invariant (Nagumo; Pomet–Praly 1992).

**Lemma P.1 (Proj_K is non-expansive).** For every ϑ = (a,c,s) ∈ Θ_K and y = (y_a, y_c, y_s),

    ‖Proj_K(ϑ, Gy)‖²_{G⁻¹} = ‖y‖²_G − Σ_{b∈{a,c,s}} h_b(2−h_b) [(n_bᵀG_b y_b)₊]² / (n_bᵀG_b n_b)  ≤  ‖y‖²_G.

The terms are as follows.

- For the balls, n_a = a and n_c = c, with the quadratic ramps h_a = h(‖a‖; a_max, ε_K) and h_c = h(‖c‖; c_max, ε_K).
- For s < s_min + ε_K, n_s = −1 and h_s = 1 − clip((s − s_min)/ε_K, 0, 1).
- For s > s_max − ε_K, n_s = +1 and h_s = 1 − clip((s_max − s)/ε_K, 0, 1).
- Otherwise h_s = 0.

*Proof.* G is block-diagonal, so both squared norms split over the blocks. On an active ball block,
Proj_b = G_b q_b with q_b = y_b − h_b n_b (n_bᵀG_b y_b)/(n_bᵀG_b n_b). Expanding q_bᵀG_b q_b gives the subtracted
term. On an inactive block (h_b = 0 or n_bᵀG_b y_b ≤ 0), Proj_b = G_b y_b. The s-block is the one-dimensional case.
Finally h_b ∈ [0,1], so h_b(2 − h_b) ∈ [0,1]. ∎

**Potential.** g_t(ϑ) := e(t)ᵀ Φ(x(t); ϑ), and F_t(μ) := −∫ g_t dμ = −e(t)ᵀ f̂(x(t)).

**Adaptive law (the algorithm).**

    ϑ̇_i = Proj_K(ϑ_i, G ∇_ϑ g_t(ϑ_i)),   i = 1..m.                                                  (L)

**Ideal measure.** For j = 1..m, ϑ*_j(t) ∈ Θ_K^{ε} is C¹ in t, and μ*_t := (1/m)Σδ_{ϑ*_j(t)}. Write
ε(x,t) := f(x,t) − ∫Φ(x;·)dμ*_t, with ‖ε‖ ≤ ε_max on the operating set X, and set

    V_*(t) := ((1/m)Σ_j ‖ϑ̇*_j‖²_{G⁻¹})^{1/2},   V_* := sup_t V_*(t).

The proofs never use optimality of μ*. It is an analytic device, and any admissible C¹ curve may be used.

**Optimal-transport distance.** For uniform atomic μ = (1/m)Σδ_{ϑ_i} and ν = (1/m)Σδ_{ϑ'_j},

    W(μ, ν)² := min_{π ∈ S_m} (1/m) Σ_i ‖ϑ_i − ϑ'_{π(i)}‖²_{G⁻¹}.

This is the Kantorovich 2-Wasserstein distance with ground cost ‖·‖²_{G⁻¹}. For equal masses the optimum is attained
at a permutation (Birkhoff). It is computed by linear assignment in O(m³). For general measures on Θ_K, W denotes the
same Kantorovich distance.

**Label-free Lyapunov function.**

    V_OT(t) := ½‖e(t)‖² + ½ W(μ_t, μ*_t)² = min_π V_π(t),
    V_π := ½‖e‖² + (1/2m) Σ_i ‖ϑ_i − ϑ*_{π(i)}‖²_{G⁻¹}.

**Standing hypotheses (S).**

- **(S0)** k_e > 0; G = diag(g_a I_n, g_c I_n, g_s) ≻ 0; s_min > 0; and 0 < ε_K < min{a_max, c_max, (s_max − s_min)/2}
  (`projection.py` enforces the last condition).
- **(S1)** f is continuous in (x,t) and locally Lipschitz in x, uniformly in t on bounded time intervals.
- **(S2)** x_m ∈ C¹, and x_m and ẋ_m are bounded.
- **(S3)** ϑ_i(0) ∈ Θ_K.
- **(S4)** The ϑ*_j are C¹ with values in Θ_K^{ε}, and V_* < ∞.
- **(S5)** x(t) ∈ X on the interval considered. For bounded f, take X := ∪_t{x : ‖x − x_m(t)‖ ≤ r(F)} with
  F = sup‖f‖ (Theorem D(b)). Then ε_max is the supremum over this X, and it depends on ‖e(0)‖.

Under (S0)–(S3) the right-hand side of (L) is continuous in t and locally Lipschitz in (x, ϑ). Solutions are
classical (C¹) and unique on their maximal interval. Under (S4) as well, each V_π is C¹.

**Notation.**

- Times are written t₀ ≤ t₁ (or T₁ < T₂), and r is the integration variable in time.
- RMS_{[T₁,T₂]}‖e‖ := ((T₂−T₁)⁻¹∫_{T₁}^{T₂}‖e‖²dr)^{1/2}.

**Constants of the Gaussian** (closed forms; `transport.gaussian_constants`). With u := ‖x−c‖²/s² and
φ = e^{−u/2}, on s ≥ s_min:

    ‖∇_c φ‖ ≤ J_c := e^{−1/2}/s_min,   |∂_s φ| ≤ J_s := 2e^{−1}/s_min,   J := (J_c² + J_s²)^{1/2},
    ‖∇²_{(c,s)} φ‖ ≤ H_φ := (max(h_cc, h_ss) + h_cs)/s_min²,
    h_cc = 1,  h_ss = 2e^{−1/2},  h_cs = √u e^{−u/2}|u−2| evaluated at u = (5−√17)/2.

J and H_φ are upper bounds. The true supremum of ‖∇_{(c,s)}φ‖ is about 0.92/s_min, and that of ‖∇²_{(c,s)}φ‖ is
about 1.57/s_min². The Hessian of the lifted potential is the block matrix [[0, u∇φᵀ],[∇φuᵀ, (uᵀa)∇²φ]]. Hence

    sup_{x ∈ ℝⁿ, ϑ ∈ Θ_K, ‖u‖=1} ‖∇²_ϑ (uᵀΦ(x;ϑ))‖ ≤ H_K := J + a_max H_φ,     K := ½H_K.          (1)

A brute-force check over the full lifted Hessian gives sup ≈ 0.76·H_K (review round 1).

---

## 2. Theorem A — identification

**Theorem A.** Assume (S0).

(i) Assume scalar classical gains Γ_W = γ_W I, Γ_c = γ_c I and Γ_σ = γ_σ I, and use the dictionary W_i = a_i/m,
γ_W = g_a/m, γ_c = m g_c, γ_σ = m g_s.

- The *raw* (unprojected) rates of (L) and of the classical law of `THEORY_REVISED.md` §2 coincide.
- The *projected* rates coincide when the classical joint projections are replaced by per-kernel ones. The classical
  law projects W onto one Frobenius ball and vec(c) onto one Euclidean ball. The per-kernel version projects:
  - W_i onto ‖W_i‖ ≤ a_max/m, with layer ε_K/m;
  - c_i onto ‖c_i‖ ≤ c_max, with layer ε_K;
  - σ_i onto [s_min, s_max], with layer ε_K.
- For m ≥ 2, no joint ball, of any radius W_max and layer ε_W, reproduces the per-kernel field. For m = 1 the two
  coincide, with W_max = a_max and layer ε_K.

(ii) Fix a closed-loop trajectory (x(·), e(·)). Then:

- The field v_t(ϑ) := Proj_K(ϑ, G∇g_t(ϑ)) is Lipschitz in ϑ on Θ_K, and Θ_K is forward invariant.
- For every μ_0 ∈ P(Θ_K), μ_t := (X_t)#μ_0 is a distributional solution of ∂_tμ + ∇·(μ v_t) = 0.
- On int Θ_K^{ε}, v_t = −G∇(δF_t/δμ). So at frozen time, μ_t is the Wasserstein gradient flow (metric G⁻¹) of the
  linear functional F_t, and (L) is its particle form.
- F_t is generated by the closed loop, so the coupled system is the one studied in Theorem D.

(iii) **Dictionary**, all in the G⁻¹ metric:

| classical device | measure-level object |
|---|---|
| projection onto Θ_K | A velocity-level (Pomet–Praly) regularisation of the support constraint μ_t ∈ P(Θ_K). For ε_K > 0 it is not the gradient flow of any functional of μ. Formally, as ε_K → 0, the field tends to the G⁻¹-projection of G∇g_t onto the tangent cone of Θ_K, which is the velocity of the W_{G⁻¹}-gradient flow of F_t + ι_{P(Θ_K)} ([46] Prop. 3.1 for the exact W₂ projection). |
| σ-modification in the G⁻¹ metric (leak outside the gain, ϑ̇ = G∇g − σ(ϑ − ϑ⁰), as in (L_σ) and Theorem H) | the potential energy (σ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ, which is geodesically σ-convex |
| classical σ-modification Γ(y − σθ) (leak inside the gain; `src/baselines/sigma_mod.py` implements its W-block, with c and σ frozen) | For the full-vector leak, under (i), the object is (σ/2)∫(ϑ−ϑ⁰)ᵀD(ϑ−ϑ⁰)dμ with D = diag(I_n/m, mI_n, m). In the G⁻¹ metric it is λ-convex with λ = σ·min(g_a/m, m g_c, m g_s), which depends on m. For the W-block alone it is (σ/2)∫(a − a⁰)ᵀ(I_n/m)(a − a⁰)dμ. |
| σ-modification toward a labelled (ϑ⁰_i) | not a functional of μ; it is a functional on the labelled lift ℝ^{mp} |
| OT σ-modification (Theorem C) | (σ/2)W(μ,μ⁰)². For atomic μ with a unique optimal assignment, its W-gradient at atom i is σ(ϑ_i − ϑ⁰_{τ(i)}). At ties, (L_σ) uses a measurable selection from the superdifferential {σ(ϑ_i − ϑ⁰_{τ(i)}) : τ optimal}. |
| e-modification | the same functionals with coefficient σ‖e‖ |
| parameter dither (Theorem H) | the entropy τ∫μ log μ |

*Proof.*

(i) **Gain map.** We have ∇_a g = φe, ∇_c g = (aᵀe)∇_cφ and ∂_s g = (aᵀe)∂_sφ. Matching these with the classical rates
gives the gain map.

**Scaling identity.** (1/m)Proj_{r,ε}(a, y) = Proj_{r/m, ε/m}(a/m, y/m). Three facts give it:

- h_{r/m,ε/m}(‖a‖/m) = h_{r,ε}(‖a‖);
- the direction a(aᵀy)/‖a‖² is homogeneous of degree 0 in a and of degree 1 in y;
- the activation aᵀy > 0 is invariant under scaling.

If the layer is not rescaled, the ramps differ.

*Only if (m ≥ 2).* Let r_in := (a_max − ε_K)/m, and let the joint ball have radius W_max, layer ε_W and
ρ := W_max − ε_W. In both cases below the driving term points outward.

- If ρ < √2 r_in, take two rows of equal norm r with √2 r ∈ (ρ, min(√2 r_in, W_max)], parallel to e. The per-kernel
  field is inactive and the joint field is active.
- If ρ ≥ √2 r_in, take one row with norm in (r_in, min(a_max/m, ρ)) and the others zero. The per-kernel field is
  active and the joint field is inactive.

The same argument applies to vec(c). For m = 1 the two operators coincide.

`test_identification_projected_rates_per_kernel` checks three things: the identity with active projections; that the
layer must be rescaled; and the only-if construction against `ProjectionOperator.project_matrix` on five joint balls. Four
are active against an inactive per-kernel field, and one is inactive against an active per-kernel field.

(ii) Proj_K is Lipschitz in ϑ. The ramp vanishes for ‖ϑ‖ ≤ θ_max − ε_K, so the direction ϑ/‖ϑ‖ is evaluated only
away from 0. Also, ∇g_t is C¹ on s ≥ s_min.

- **Invariance.** Forward invariance follows from Nagumo's theorem, since n_bᵀProj_b ≤ 0 on the boundary, with the outward normals n_b of Lemma P.1.
- **Continuity equation.** The push-forward of μ_0 by a Lipschitz flow solves the continuity equation (Golse [82];
  AGS §8.1).
- **Gradient form.** On int Θ_K^{ε} the projection is inactive, so v_t = G∇g_t = −G∇(δF_t/δμ).

(iii) **Potential-energy rows.** These follow from the computation in the metric G⁻¹: for any V, ∇_{G⁻¹}V = G∇V.
For example, V = (σ/2)‖ϑ−ϑ⁰‖²_{G⁻¹} gives G∇V = σ(ϑ−ϑ⁰).

**Classical-leak row.** Applying the gain map of (i) to Γ(y − σθ) gives the velocity −σGD(ϑ−ϑ⁰).

**Projection row.** The G⁻¹-rescaled correction field has an asymmetric Jacobian. In the example m = 1, n = 1, G = I,
x = 0.7, e = 1, ϑ = (0.93, 0.2, 1.5), a_max = 1, ε_K = 0.15 (so h(0.93) = 0.284) and c_max = 2, with a in the
layer:

- the c-derivative of the a-component is −0.060;
- the a-derivative of the c-component is 0 (review round 2).

So the field is not a gradient. ∎

*Remark A.1 (normalisation).* The algorithm's raw rates are unchanged. Two things do change: the admissible set becomes
a product of per-kernel balls, and a_max, not ‖W‖_F, is held fixed in m.

- **m-free constants.** In this normalisation the structural constants k_e, K, L₀, L₁, λ_G and Δ are free of m for
  fixed G and a_max. So are the constants of Theorems B, C(i)–(iii), D and F.
- **m-dependent data.** ε_max and V_* are data of the chosen m-atom ideal. The local constants of C(iv) and of Theorem
  E depend on m.
- **Classical convention.** Fixed classical γ's correspond to g_a = mγ_W ∝ m and g_c, g_s ∝ 1/m. Weights speed up,
  centres and widths freeze, λ_G ∝ m and Δ ∝ √m.

---

## 3. Theorem B — the optimal-assignment Lyapunov function

**Lemma B.1 (Dini derivative of a minimum).** Suppose each V_π is C¹ and V_OT = min_π V_π. Then V_OT is locally
Lipschitz, and D⁺V_OT(t) ≤ V̇_π(t) for every π active at t.

*Proof.* For h > 0, V_OT(t+h) ≤ V_π(t+h), and V_OT(t) = V_π(t). ∎

**Lemma B.2 (matched remainder in two regimes; gradient bound).** Let L₀ := 2a_max and L₁ := (1 + a_max²J²)^{1/2}.
For ϑ, ϑ' ∈ Θ_K:

(a) the remainder R := g_t(ϑ') − g_t(ϑ) − ⟨∇g_t(ϑ), ϑ' − ϑ⟩ satisfies
|R| ≤ ‖e(t)‖ · min{ K‖ϑ − ϑ'‖², L₀ + L₁‖ϑ − ϑ'‖ };

(b) ‖∇_ϑ g_t(ϑ)‖ ≤ L₁‖e(t)‖.

*Proof.* (b) We have ∇_a g = φe and ∇_{(c,s)}g = (aᵀe)∇_{(c,s)}φ, with φ ≤ 1 and ‖∇_{(c,s)}φ‖ ≤ J on s ≥ s_min. So
‖∇g‖² ≤ ‖e‖²(1 + a_max²J²).

(a) If e(t) = 0, both sides vanish. Otherwise set h := ϑ' − ϑ.

- *Quadratic branch.* Θ_K is convex, so the Lagrange remainder of τ ↦ g_t(ϑ + τh) gives R = ½hᵀ∇²g_t(ϑ + τ₀h)h for
  some τ₀ ∈ (0,1). Since g_t = ‖e‖·(uᵀΦ) with u = e/‖e‖, bound (1) gives this branch.
- *Linear branch.* It follows from |g_t| ≤ ‖e‖a_max and (b). ∎

**Theorem B.** Under (S), along every solution of (L), for every t, with W := W(μ_t, μ*_t):

    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖ ( ε_max + min{ Kλ_G W², L₀ + L₁ λ_G^{1/2} W } ) + V_*(t) W.        (B)

*Proof.* Fix t and an optimal π, and write ϑ̃_i := ϑ_i − ϑ*_{π(i)}. Then

    V̇_π = eᵀė + (1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇_i⟩ − (1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇*_{π(i)}⟩.

(a) eᵀė = −k_e‖e‖² + eᵀε + (1/m)Σ_i [g_t(ϑ*_{π(i)}) − g_t(ϑ_i)].

(b) By (P) with ϑ' = ϑ*_{π(i)}, (1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇_i⟩ ≤ (1/m)Σ⟨ϑ̃_i, ∇g_t(ϑ_i)⟩.

(c) −(1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇*_{π(i)}⟩ ≤ ((1/m)Σ‖ϑ̃_i‖²_{G⁻¹})^{1/2}((1/m)Σ_i‖ϑ̇*_{π(i)}‖²_{G⁻¹})^{1/2} = W·V_*(t),
since π is a bijection and is optimal.

Summing gives V̇_π ≤ −k_e‖e‖² + ‖e‖ε_max + (1/m)Σ_i R_i + W V_*(t), where each R_i is a matched remainder. Average
Lemma B.2(a), using three facts:

- (1/m)Σmin{A_i, B_i} ≤ min{avg A, avg B};
- (1/m)Σ‖ϑ̃_i‖² ≤ λ_G W²;
- (1/m)Σ‖ϑ̃_i‖ ≤ λ_G^{1/2}W.

Lemma B.1 then completes the proof. ∎

**Corollary B.3 (ultimate bounds).** Let Δ := diam_{G⁻¹}(Θ_K), so W ≤ Δ. Define

    D^{OT} := ε_max + min{Kλ_GΔ², L₀ + L₁λ_G^{1/2}Δ},    C^{OT} := V_*Δ,
    R^{OT} := D^{OT}/(2k_e) + (D^{OT 2}/(4k_e²) + C^{OT}/k_e)^{1/2}.

Then:

- D⁺V_OT < 0 whenever ‖e‖ > R^{OT}.
- Under (S5) on [0, ∞), limsup ‖e‖ ≤ (R^{OT 2} + Δ²)^{1/2}. This uses sup ½W² ≤ Δ²/2.
- RMS_{[T₁,T₂]}‖e‖ ≤ D^{OT}/(2k_e) + (D^{OT 2}/(4k_e²) + C^{OT}/k_e + V_OT(T₁)/(k_e(T₂−T₁)))^{1/2}.

*Proof.* As in `THEORY_REVISED.md` Theorem 1. The Dini comparison lemma, applied with the continuous, V-independent
right-hand side, gives the integrated form. ∎

*Remark B.4 (what OT buys and what it does not).*

- **Globally, the OT constant is not smaller than the classical one.** Take the test box: a_max = 8, c_max = 4,
  s_min = 0.3, s_max = 3, G = diag(2, 1, 0.5) and F_max := a_max. Compare constants on the same per-kernel admissible
  set.

  | constant | value |
  |---|---|
  | classical Lemma 2, F_max + 3a_max + 2a_max c_max J_c + a_max(s_max − s_min)J_s | 214, for every m |
  | OT linear branch, crude form L₀ + L₁λ_G^{1/2}Δ | 533 |
  | OT linear branch, block split through W | 382 |
  | OT linear branch, bounded per particle by ‖e‖(2a_max + ‖ã_i‖ + a_max(J_c‖c̃_i‖ + J_s|s̃_i|)) | 214 |

  Bounded per particle, the OT constant equals the classical one. Version 2's comparison figure of "101–138" used a
  joint centre ball of radius c_max. That radius is √m times smaller than the radius √m·c_max of the ball that
  circumscribes the per-kernel product of Theorem A.
- **Locally, the excess disturbance is quadratic in the matched error** (Corollary B.5).
- **Labels.** A labelled theory can make the same local statement only for a labelling fixed in advance, and only
  while the optimal assignment does not switch. (B) holds uniformly, with the minimum over π taken at each t.
- **m-uniformity** is a property of the normalisation (Remark A.1), not of OT.

**Corollary B.5 (conditional bound).** Along a solution of (L), suppose sup_{[T₁,T₂]} W(μ_t, μ*_t) ≤ η̄. This quantity
is available in oracle simulations, since it requires μ*_t, or it can be controlled a priori by Theorem E. With
D̄ := ε_max + min{Kλ_Gη̄², L₀ + L₁λ_G^{1/2}η̄},

    RMS_{[T₁,T₂]}(‖e‖) ≤ D̄/(2k_e) + ( D̄²/(4k_e²) + V_*η̄/k_e + V_OT(T₁)/(k_e(T₂−T₁)) )^{1/2}.

As η̄ → 0 and T₂ → ∞, this tends to ε_max/k_e, the value of the frozen-parameter error channel.

*Variant for (L_σ).* Take a solution of (L_σ) in the sense of §4 with sup_{[T₁,T₂]} W(μ_t, μ*_t) ≤ η̄ and
4η̄(b + η̄) < Δ_gap. For example, η̄ can come from Theorem C(iv). Then, with D̄ := ε_max + Kλ_Gη̄²,

    RMS_{[T₁,T₂]}(‖e‖) ≤ D̄/(2k_e) + ( D̄²/(4k_e²) + (V_*η̄ + σ_m b²/4)/k_e + V_OT(T₁)/(k_e(T₂−T₁)) )^{1/2}.

*Proof.* Apply Theorem C(iii) with C = b(t)² ≤ b², and use −σ_mη² + σ_mηb ≤ σ_m b²/4. Without the gap condition,
replace b² by (b + 2η̄)² (from √C ≤ b + 2η in C(ii)), or by b² + 4η̄R. ∎

As η̄ → 0 and T₂ → ∞ the variant tends formally to ε_max/(2k_e) + (ε_max²/(4k_e²) + σ_m b²/(4k_e))^{1/2}. Along (L_σ)
with b > 0, the measured η̄ stays of order b (about 0.9b in the test instance), so the bound is first order in b. The plain bound (B.5) is false along (L_σ). In the test instance, every hypothesis of C(iv) is asserted and
sup W ≤ η̄ is measured. There, B.5 is exceeded by more than 3×, while the variant holds
(`test_corollary_B5_fails_along_lsigma_and_its_leak_variant_holds`). The round-2 reviewer found violations of 3.8–8.3×.

**Corollary B.6 (approximation error as quantisation error).** Suppose f = ∫Φ dμ*_∞ + ε_∞ for a probability measure
μ*_∞ on Θ_K^{ε}, with ‖ε_∞‖ ≤ ε̄_∞. Then for every m-atom μ*_m supported in Θ_K^{ε},

    ε_max(m) ≤ ε̄_∞ + L₁λ_G^{1/2} W(μ*_m, μ*_∞).

*Proof.* For any coupling γ,

    ‖∫Φ d(μ*_m − μ*_∞)‖ ≤ ∫‖Φ(ϑ) − Φ(ϑ')‖dγ ≤ L₁λ_G^{1/2}∫‖ϑ − ϑ'‖_{G⁻¹}dγ ≤ L₁λ_G^{1/2}(∫‖ϑ − ϑ'‖²_{G⁻¹}dγ)^{1/2}.

Minimise over γ. The block constant (g_a + a_max²(J_c²g_c + J_s²g_s))^{1/2} may replace L₁λ_G^{1/2}. ∎

*Rates.*

(a) **First order.** Let μ*_∞ be compactly supported and p ≥ 3. Then there are equal-weight μ*_m, with atoms (cell
barycentres) in Θ_K^{ε}, such that W(μ*_m, μ*_∞) = O(m^{−1/p}). This follows from an equal-mass dyadic decomposition
(Chevallier 2018 [S13]). For p > 4, i.i.d. atoms also achieve the rate in expectation (Fournier–Guillin 2015). If
μ*_∞ has a non-trivial absolutely continuous part, the order cannot be improved, even with free weights (Zador;
Graf–Luschgy [S9] Thm 6.2 with Lemma 3.4).

(b) **Second order.** Suppose μ*_m is a centroid quantiser for an optimal plan; a best equal-weight approximation is
one. Then the first-order term cancels, and by (1) and convexity of Θ_K^{ε},

    ε_max(m) ≤ ε̄_∞ + min{Kλ_G W², L₁λ_G^{1/2}W} = ε̄_∞ + O(m^{−2/p}).

(c) **Free of p.** I.i.d. atoms drawn from μ*_∞ give, with positive probability, an m-atom ideal with
ε_max(m) ≤ ε̄_∞ + O(a_max(n log(m·diam X)/m)^{1/2} + a_max/(s_min√m)). This is Maurey's empirical method; the
constants are not optimised.

The exponent 1/p is therefore a property of controlling ε_max through W, as Theorem D(e) requires. It is not a
property of the approximation error itself.

*Remark B.7 (continuum form, via couplings).* Let μ_0 ∈ P(Θ_K) and μ_t = (X_t)#μ_0, and let μ*_t be an absolutely
continuous ideal curve in (P(Θ_K^{ε}), W) with velocity field w_t. Set V_* := ess sup_t ‖w_t‖_{L²(μ*_t;G⁻¹)}. Its
smallest admissible value is the metric derivative |μ̇*_t|_W (AGS Thm 8.3.1). Then (B) holds for a.e. t.

*Sketch.*

1. Push an optimal plan γ_t forward by (X_{t,t+h}, Y_{t,t+h}). Here Y is the flow of w, or the superposition measure
   of AGS Thm 8.2.1 when w is not Lipschitz.
2. This couples (μ_{t+h}, μ*_{t+h}), so D⁺(½W²) ≤ ∫⟨ϑ − ϑ', v_t(ϑ) − w_t(ϑ')⟩_{G⁻¹}dγ_t.
3. Integrate (P) and Lemma B.2 against γ_t. ∎

---

## 4. Theorem C — optimal-transport σ-modification

Let μ⁰ = (1/m)Σδ_{ϑ⁰_k} be a prior configuration with ϑ⁰_k ∈ Θ_K^{ε}. The **OT σ-modified law** is

    ϑ̇_i = Proj_K( ϑ_i, G∇g_t(ϑ_i) − σ_m(ϑ_i − ϑ⁰_{τ(i)}) ),   τ(t) ∈ argmin_ρ (1/m)Σ_i‖ϑ_i − ϑ⁰_{ρ(i)}‖²_{G⁻¹}.        (L_σ)

**Solution concept.** The right-hand side of (L_σ) is discontinuous at assignment ties.

- **Selection solutions.** A solution is an absolutely continuous (x, ϑ) such that, for a.e. t,
  ϑ̇_i = Proj_K(ϑ_i, G∇g_t(ϑ_i) − σ_m(ϑ_i − ϑ⁰_{τ(t)(i)})) for some τ(t) optimal at ϑ(t). This is what `transport.py`
  integrates.
- **Filippov solutions.** These exist by Filippov's theorem, since the right-hand side is measurable and locally
  bounded. But they include sliding motions on repelling ties, along which (i) below fails. For example, take
  ideal = prior = {a = 0, a = 1}, both particles at a = ½ and e ≡ 0: this is a stationary Filippov solution.
- **What holds for which class.**
  - Parts (ii)–(iv) hold for every Filippov solution, with C taken as the maximum over optimal τ.
  - Part (i) and the exact decay hold for selection solutions.
  - For Filippov solutions, (i) weakens as stated below.
- **Regularity and existence.** V_OT and Z are absolutely continuous along solutions, and every inequality holds for
  a.e. t. Existence of selection solutions is proved for the pure leak (Lemma C.3). In general it is assumed: two-way
  ties are never attracting, by monotonicity of Proj_K, but multi-way ties, the degenerate ramp h = 1 and accumulating
  switch times are not handled.

**Notation.**

- c_ρ(ϑ) := (1/m)Σ‖ϑ_i − ϑ⁰_{ρ(i)}‖²_{G⁻¹}.
- η := W(μ_t, μ*_t), and R := Δ = diam_{G⁻¹}(Θ_K).
- b(t) := W(μ*_t, μ⁰), and b := sup_t b(t).
- Δ_gap(t) := min{c_ρ(ϑ*(t)) − b(t)² : c_ρ(ϑ*(t)) > b(t)²}, with min ∅ := +∞; and Δ_gap := inf_t Δ_gap(t).
- A_{ττ'} := (1/m)Σ‖ϑ⁰_{τ(i)} − ϑ⁰_{τ'(i)}‖²_{G⁻¹}.

**Theorem C.** Under (S), along every solution of (L_σ) in the sense above:

(i) *Dissipation toward the prior (selection solutions).* Let Z := ½W(μ_t, μ⁰)² and
Γ := ((1/m)Σ‖G∇g_t(ϑ_i)‖²_{G⁻¹})^{1/2} ≤ λ_G^{1/2}L₁‖e‖. Then:

- D⁺Z ≤ −2σ_mZ + (2Z)^{1/2}Γ(t);
- hence W(μ_t, μ⁰) ≤ max{W(μ_0, μ⁰), sup Γ/σ_m};
- for the pure-leak flow (Γ ≡ 0), the decay W(μ_t, μ⁰) = e^{−σ_m t}W(μ_0, μ⁰) is exact for every G ≻ 0 (Lemma C.3).

For a Filippov solution with weights θ_k on tied optimal assignments τ_k, the bound weakens to

    D⁺Z ≤ −2σ_mZ + (σ_m/2)min_τ Σ_kθ_kA_{ττ_k} + (2Z)^{1/2}Γ,

with the minimum over assignments τ optimal at t.

(ii) *Induced coupling.* Let π* be optimal from μ to μ*, and let τ be the optimal assignment from μ to μ⁰ that the law
uses at t (for Filippov solutions, any optimal τ). Set C := (1/m)Σ‖ϑ*_{π*(i)} − ϑ⁰_{τ(i)}‖²_{G⁻¹}. Then

    C ≤ b(t)² + 4ηR   and   √C ≤ b(t) + 2η,

and C = b(t)² whenever 4η(b(t) + η) < Δ_gap.

(iii) *Lyapunov inequality.* With C as in (ii) (for Filippov solutions, its maximum over optimal τ),

    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖(ε_max + Kλ_Gη²) + V_*η − σ_mη² + σ_mηC^{1/2}.

Under 4η(b + η) < Δ_gap this becomes

    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖(ε_max + Kλ_Gη²) + V_*η − (σ_m/2)η² + (σ_m/2)b².

(iv) *Local exponential estimate.* Suppose 4η̄(b + η̄) < Δ_gap and 8K²λ_G²η̄² ≤ k_eσ_m. Set α := min(k_e, σ_m/2) and
β := ε_max²/k_e + (σ_m/2)b² + 2V_*²/σ_m, and suppose max{V_OT(0), β/α} < η̄²/2. Then:

- {V_OT ≤ c} is forward invariant for every c ∈ [max{V_OT(0), β/α}, η̄²/2);
- V_OT(t) ≤ e^{−αt}V_OT(0) + β/α;
- on this set ‖e‖ < η̄, so (S5) holds with X := ∪_t{x : ‖x − x_m(t)‖ ≤ η̄}.

(v) *The labelled leak.* Define b_lab(π)(t) := ((1/m)Σ‖ϑ*_{π(i)}(t) − ϑ⁰_i‖²_{G⁻¹})^{1/2} and
b_lab(π) := sup_t b_lab(π)(t).

- For any fixed π, the labelled σ-modification (leak toward ϑ⁰_i) satisfies the second display of (iii) at every
  state, with no gap condition. Here b is replaced by b_lab(π), V_OT by V_π, and η by the π-labelled error.
- Choosing π to minimise b_lab(π) gives the bias min_π sup_t b_lab(π)(t) ≥ b. Equality holds when one assignment is
  optimal for (μ*_t, μ⁰) at every t, for example a static ideal, or whenever Δ_gap > 0.
- The labelled leak does **not** satisfy (iii) for V_OT (Proposition C.5).

**Scope.**

- **Finite gap.** If Δ_gap < ∞, that is, not all ideal–prior assignments tie, then Δ_gap ≤ 2R²/m always. Indeed, if
  every transposition of an optimal ρ₀ were cost-neutral, the matrix S(i,j) := B(i,j) − (B(i,i) + B(j,j))/2 built
  from B(i,j) := ⟨ϑ*_i, ϑ⁰_{ρ₀(j)}⟩_{G⁻¹} would be antisymmetric, and every assignment would tie.
  - Hence (iv) needs η̄ < √Δ_gap/2 ≤ R/√(2m).
  - β/α < η̄²/2 forces b < η̄/√2.
  - So (iv) is a near-prior local result of size O(m^{−1/2}) in W, and O(1/m) in V_OT.
- **Infinite gap.** If Δ_gap = +∞, for example for a single-atom prior (the classical σ-modification), the gap
  condition is void. Then (iv) is m-uniform, with η̄ fixed by 8K²λ_G²η̄² ≤ k_eσ_m alone.
- **Switching.** Δ_gap > 0 forces the set of optimal ideal–prior assignments to be constant in t. If that set switches
  along a moving ideal, then Δ_gap = 0, and the second clause of (ii), the second display of (iii) and all of (iv)
  are vacuous.

*Proof.*

(i) For a selection solution, at a.e. t the law uses an optimal τ(t), and Z ≤ Z_τ := ½c_τ with equality at t. By
Lemma B.1 it suffices to bound Ż_τ. Apply (P) with ϑ' = ϑ⁰_{τ(i)} and Cauchy–Schwarz in the G⁻¹ inner product:
Ż_τ ≤ (2Z)^{1/2}Γ − 2σ_mZ. Write y := (2Z)^{1/2}.

- Where Z > 0, D⁺y ≤ −σ_m y + Γ.
- At Z = 0, D⁺y ≤ ((1/m)Σ‖ϑ̇_i‖²_{G⁻¹})^{1/2} ≤ Γ, by Lemma P.1 (the leak term vanishes there).
- Γ ≤ λ_G^{1/2}L₁‖e‖ by Lemma B.2(b).

The pure-leak case is Lemma C.3, and the projection is inactive there. In the layer,
aᵀ(a⁰ − a) ≤ ‖a‖(‖a⁰‖ − ‖a‖) ≤ 0, and similarly for c. For s, the leak points away from the nearby face.

For Filippov velocities, the leak part averaged over the tie is
−σ_m(1/m)Σ_iΣ_kθ_k⟨ϑ_i − ϑ⁰_{τ(i)}, ϑ_i − ϑ⁰_{τ_k(i)}⟩_{G⁻¹}. By polarisation, and since c_τ = c_{τ_k} = 2Z for tied
optimal assignments, this equals −2σ_mZ + (σ_m/2)Σ_kθ_kA_{ττ_k}. (P) holds termwise for the convex combination. The
bound is tight at the tie equilibrium of the solution-concept example (0 = −σ_m/4 + σ_m/4).

(ii) Relabel so that (1/m)Σ‖ϑ_i − ϑ*_i‖²_{G⁻¹} = η². For any ρ:

- |c_ρ(ϑ) − c_ρ(ϑ*)| ≤ 2ηR;
- |√c_ρ(ϑ) − √c_ρ(ϑ*)| ≤ η, by the reverse triangle inequality for the (1/m)-weighted G⁻¹ norm.

Let ρ₀ be optimal for (ϑ*(t), ϑ⁰). Then

    C = c_τ(ϑ*) ≤ c_τ(ϑ) + 2ηR ≤ c_{ρ₀}(ϑ) + 2ηR ≤ b(t)² + 4ηR,
    √C ≤ √c_τ(ϑ) + η ≤ √c_{ρ₀}(ϑ) + η ≤ b(t) + 2η.

So C − b(t)² ≤ 4η(b(t) + η) < Δ_gap ≤ Δ_gap(t). No assignment cost at time t lies strictly between b(t)² and
b(t)² + Δ_gap(t), so C = b(t)².

(iii) Apply Theorem B with y = ∇g_t − σ_mG⁻¹(ϑ_i − ϑ⁰_{τ(i)}) in (P). The extra term is
−σ_m(1/m)Σ⟨ϑ̃_i, ϑ̃_i + ϑ*_{π*(i)} − ϑ⁰_{τ(i)}⟩_{G⁻¹} ≤ −σ_mη² + σ_mηC^{1/2}. Under the gap condition, (ii) gives
C^{1/2} = b(t) ≤ b, and σ_mηb ≤ (σ_m/2)(η² + b²).

(iv) Apply Young's inequality three times:

- ‖e‖ε_max ≤ (k_e/4)‖e‖² + ε_max²/k_e;
- ‖e‖Kλ_Gη² ≤ (k_e/4)‖e‖² + (σ_m/8)η² on {η ≤ η̄};
- V_*η ≤ (σ_m/8)η² + 2V_*²/σ_m.

This gives D⁺V_OT ≤ −αV_OT + β. Each set {V_OT ≤ c} with β/α < c < η̄²/2 lies inside {η ≤ η̄} and has a negative
boundary derivative. Intersect over c > β/α for the endpoint, and apply the comparison lemma.

(v) Steps (a)–(c) of Theorem B hold for any fixed π, so use the classical completion of squares for V_π.

- *Equality case.* If Δ_gap > 0, the costs minus b(t)² take values in {0} ∪ [Δ_gap, ∞) and are continuous in t, so
  the optimal set is constant.
- *Strict case.* If the optimal ideal–prior assignment changes along the ideal, the inequality can be strict. For
  example, a*_1(t) = t and a*_2(t) = 1 − t on [0,1], with prior a ∈ {0, 1}, give b = 0.5 against 1.0.
- *Not always strict.* For a*(t) = (u(t), 1 − u(t)) with u(t) = −0.025 + 0.575 sin t, the optimal assignment alternates
  between the identity and the swap, yet min_π sup_t b_lab(π)(t) = sup_t b_lab(id)(t) = 0.6 = b (review round 3).

The negative statement is Proposition C.5. ∎

**Lemma C.3 (no switching along the pure leak; exact decay for any G ≻ 0).** Let τ be optimal for (μ_init, μ⁰) and
ϑ̇_i = −σ_m(ϑ_i − ϑ⁰_{τ(i)}). Then τ stays optimal and W(μ_t, μ⁰) = e^{−σ_m t}W(μ_init, μ⁰).

*Proof.* Put d_i := ϑ_i(0) − ϑ⁰_{τ(i)}, so ϑ_i(t) − ϑ⁰_{τ(i)} = λd_i with λ = e^{−σ_m t}. For any ρ, with
q_i := ϑ⁰_{τ(i)} − ϑ⁰_{ρ(i)},

    c_ρ(t) − c_τ(t) = (2λ/m)Σ⟨d_i, q_i⟩_{G⁻¹} + (1/m)Σ‖q_i‖²_{G⁻¹}.

This is affine in λ. It is non-negative at λ = 1 by optimality and at λ = 0 trivially, hence on all of [0,1].

At an exact tie at t = 0, the chosen assignment becomes strictly optimal immediately against every ρ with
A_ρ := (1/m)Σ‖q_i‖² > 0. Assignments with A_ρ = 0 induce the same velocity. So no selection solution of the pure leak
switches. Filippov sliding on the tie is a different solution, excluded by the solution concept. ∎

Checked: zero switches, and exact decay to integrator accuracy.

**Proposition C.5 (the labelled leak reshuffles).** There are configurations with b = b_lab(id) = 0 and e = 0 at
which the labelled leak gives D⁺V_OT > 0.

*Example.* Take n = 1, G = I and a_max − ε_K ≥ 1. Let ideal = prior = {a = 0, a = 1}, with the other coordinates equal
and in the inner set, and let the network weights be (1−δ, δ) with 0 < δ < ½.

- The optimal assignment to the ideal is the swap.
- The labelled leak moves each particle toward the *other* ideal atom.
- So D⁺V_OT = σ_m(δ − δ²) > 0.

From a relabelled start, the labelled leak moves every kernel in a non-trivial cycle of the relabelling to its prior
position. In the reshuffling experiment (`results/ot_sigma_mod/transient.md`):

- the steady states coincide;
- the labelled leak's particle path is 5.8× longer (2.26 against 0.39);
- the early RMS‖e‖ is only about 4% larger (0.0771 against 0.0743 on [0,5]).

Theorem F does not price this motion. With a leak, F bounds transport from above by tracking error plus a leak term,
and the leak can pay for all of it. With δ = 0 above, e ≡ 0 while every particle crosses the domain.

*Remark C.6 (the obstruction is known).*

- **Non-convexity.** For some ν, for example a two-atom ν in dimension ≥ 2, ½W²(·,ν) is not λ-convex along geodesics
  for any λ (AGS Example 9.1.5; Santambrogio [71] §4.5 of the journal version). It is 1-semiconcave in every dimension (AGS Thm 7.3.2). For
  ν = δ_{ϑ⁰} it is a potential energy, and 1-convex.
- **Generalised geodesics.** ½W²(·,ν) is 1-convex along generalised geodesics with base ν (AGS Lemma 9.2.1;
  Carlen–Craig [73]).
- **Label-free leaks.**
  - A label-free leak that is λ-convex (λ > 0) on the lift ℝ^{mp} has a unique, permutation-invariant minimiser, at
    which all particles coincide. So it cannot vanish at a μ⁰ with distinct atoms.
  - 0-convex leaks that vanish at μ⁰ exist, for example toward conv supp μ⁰. They do not single out μ⁰.

So the trade-off reflects this geometry. On one side is a global inequality for a labelled V_π, which is large at
relabelled starts (hence the reshuffling of Proposition C.5). On the other is a local, gap-protected inequality for
the label-free V_OT.

### 4b. Theorem C′ — the OT leak toward a prior with a different number of atoms: continuous-time equal-mass Lloyd

Theorem C needs a prior with as many atoms as the network. Here the prior is any ν ∈ P₂(Θ_K^{ε}), for example the
empirical measure of a network of a different size (the transfer problem of §8, branch 8), and the leak points at the
barycentre of the prior mass that the optimal plan assigns to each particle. The dissipation toward the prior survives
up to the quantisation error, the fixed points of the pure leak are the equal-mass centroidal configurations of ν (the
leak is Lloyd's algorithm in continuous time), and the Lyapunov bias is measured by the distance of the ideal from the
*centroids* of its prior cells, not by W(μ*, ν).

**The law.** Let ν ∈ P₂(Θ_K^{ε}), for instance ν = (1/m₀)Σ_{j=1}^{m₀}δ_{ϑ⁰_j} with any m₀ ≥ 1. For
μ = (1/m)Σδ_{ϑ_i}, let γ be a vertex optimal plan from μ to ν for the cost ‖·‖²_{G⁻¹} (for atomic ν a vertex of the
transportation polytope; for general ν any optimal plan), and let

    T̄_i := m∫ϑ' dγ(ϑ_i, ϑ')    (the barycentric projection of γ at ϑ_i; for atomic ν, T̄_i = mΣ_jγ_{ij}ϑ⁰_j),
    V_bar² := (1/m)Σ_i‖ϑ_i − T̄_i‖²_{G⁻¹},    Var_γ := ∫‖ϑ' − T̄(ϑ)‖²_{G⁻¹}dγ(ϑ,ϑ'),    so that    W(μ,ν)² = V_bar² + Var_γ.

The **barycentric OT leak** is

    ϑ̇_i = Proj_K( ϑ_i, G∇g_t(ϑ_i) − σ_m(ϑ_i − T̄_i) ).                                                        (L_σ,bar)

For atomic ν the plan is a vertex of the transportation polytope (row sums 1/m, column sums 1/m₀), computed by linear
programming (`transport.optimal_plan`, `transport.barycentric_projection`); for m₀ = m it is a permutation and (L_σ,bar)
is (L_σ). The right-hand side is discontinuous where the optimal plan is not unique, and here, unlike §4 (Lemma C.3),
the plan does switch along solutions (numerics below), so two solution concepts are used. *Vertex-selection solutions*
are absolutely continuous ϑ(t) with ϑ̇_i = Proj_K(ϑ_i, G∇g_t − σ_m(ϑ_i − T̄_i^{γ(t)})) for a.e. t and γ(t) a vertex optimal
plan at ϑ(t); the time-stepped implementation produces these (switches at step times), and their existence for all time
is assumed. *Generalised solutions* allow γ(t) to be any measurable selection of optimal plans. The optimal face F(ϑ) of
the transportation LP is a compact convex polytope, upper semicontinuous in ϑ (finitely many vertices, continuous costs),
and γ ↦ T̄^γ is linear, so the set of generalised velocities is compact, convex and upper semicontinuous in ϑ: it contains
the Filippov regularisation of the vertex law, and generalised solutions exist on [0, ∞) by the standard theorem for such
differential inclusions (the projection is handled as in §4). Every statement below holds for both concepts, with T̄ and
V_bar computed from the plan the solution uses. Since Θ_K^{ε} is convex, T̄_i ∈ Θ_K^{ε}, so property (P) applies with
ϑ' = T̄_i.

**Notation.** η := W(μ_t, μ*_t), b := sup_t W(μ*_t, ν), and

    b_bar(t)² := (1/m)Σ_j‖ϑ*_j(t) − T̄*_j(t)‖²_{G⁻¹} = W(μ*_t, ν)² − Var_{γ*(t)},    b_bar := sup_t b_bar(t),

the *centroid defect of the ideal*, computed with an optimal plan γ*(t) from μ*_t to ν. Δ_gap(t) is the gap between the
optimal and the next vertex cost of the transportation polytope at μ*_t (min ∅ := +∞), and Δ_gap := inf_t Δ_gap(t).

**Theorem C′.** Under (S), along every vertex-selection or generalised solution of (L_σ,bar):

(i) *Dissipation toward the centroidal configurations.* With Z := ½W(μ_t, ν)² and Γ as in Theorem C(i),

    D⁺Z ≤ −σ_m V_bar² + V_bar·Γ ≤ −2σ_m(Z − ½Var_γ) + (2Z)^{1/2}Γ    for a.e. t.

For the pure leak (Γ ≡ 0), W(μ_t, ν) is non-increasing, strictly decreasing wherever the plan in use is not centroidal
(V_bar > 0), and at every ω-limit point ϑ_∞ of a solution on the compact Θ_K some optimal plan γ_∞ at ϑ_∞ is centroidal:
ϑ_{∞,i} = T̄_i^{γ_∞} for all i — a vertex plan for a vertex-selection solution. Such configurations are the fixed points of
the equal-mass Lloyd iteration ϑ ← T̄^γ(ϑ) with an optimal (vertex) plan; where the optimal plan is unique they are the
critical points of the quantisation energy ϑ ↦ ½W²((1/m)Σδ_{ϑ_i}, ν), which is differentiable there, and the pure leak
is its G-gradient flow with the factor m: ϑ̇_i = −σ_m m G∂_{ϑ_i}[½W²]. They include degenerate configurations with
coincident particles. For m₀ = m, Var_γ = 0 and (i) is Theorem C(i).

(ii) *Induced coupling.* Let π* be optimal from μ to μ*, and let γ be the plan the law uses at t. Put
C_bar := (1/m)Σ_i‖ϑ*_{π*(i)} − T̄_i‖²_{G⁻¹}. Then

    √C_bar ≤ b(t) + 2η    always,    and, for atomic ν,    C_bar = b_bar(t)² whenever 4η(b(t) + η) < Δ_gap(t) and γ*(t) is unique.

(iii) *Lyapunov inequality.* With C_bar as in (ii),

    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖(ε_max + Kλ_Gη²) + V_*η − σ_mη² + σ_mη√C_bar,

and under the gap condition of (ii) the leak terms are at most −(σ_m/2)η² + (σ_m/2)b_bar².

(iv) *Local exponential estimate.* For atomic ν, Theorem C(iv) holds verbatim with b replaced by b_bar: under 4η̄(b + η̄) < Δ_gap,
uniqueness of γ*(t), and 8K²λ_G²η̄² ≤ k_eσ_m, with α := min(k_e, σ_m/2) and β := ε_max²/k_e + (σ_m/2)b_bar² + 2V_*²/σ_m,
if max{V_OT(0), β/α} < η̄²/2 then V_OT(t) ≤ e^{−αt}V_OT(0) + β/α and {V_OT ≤ c} is forward invariant for
c ∈ [max{V_OT(0), β/α}, η̄²/2).

(v) *What the bias measures.* b_bar ≤ b, with b_bar = 0 exactly when the ideal is itself an equal-mass centroidal
configuration of ν. So a prior with many atoms whose centroidal m-quantisation coincides with the ideal anchors the
network without bias, although W(μ*, ν) > 0; a prior with the right number of atoms but the wrong positions has
b_bar = b.

*Scope of (ii)–(iv).* The gap condition 4η(b + η) < Δ_gap(t) asks the network to be closer to the ideal than the optimal
vertex of the LP is to the next one. For a prior with many atoms, Δ_gap is at most the cost of exchanging two atoms
between two cells, which is O(1/m₀), so the condition holds only in a small neighbourhood of the ideal; in the
closed-loop experiment below 4η(0)(b + η(0)) exceeds that bound by factors 83–175, and (iv) applies to no tested
instance (the table reports the exchange bound on Δ_gap, η(0) and 4η(0)(b + η(0)) for each run). What the experiment
confirms is the fixed-point structure, (i) and (v): the loop settles at the centroidal configuration of its prior, whose
distance from the ideal is b_bar.

*Proof.* (i) At a.e. t the law uses an optimal plan γ, and for h > 0 the same plan couples μ_{t+h} and ν, so
W(μ_{t+h}, ν)² ≤ Σ_{ij}γ_{ij}‖ϑ_i(t+h) − ϑ⁰_j‖²_{G⁻¹} (for general ν, ∫‖ϑ_i(t+h) − ϑ'‖²dγ). Differentiating the
right-hand side at h = 0,

    D⁺Z ≤ Σ_{ij}γ_{ij}⟨ϑ_i − ϑ⁰_j, G⁻¹ϑ̇_i⟩ = (1/m)Σ_i⟨ϑ_i − T̄_i, G⁻¹ϑ̇_i⟩,

because Σ_jγ_{ij}(ϑ_i − ϑ⁰_j) = (1/m)(ϑ_i − T̄_i). Apply (P) with ϑ' = T̄_i ∈ Θ_K^{ε} and y = ∇g_t − σ_mG⁻¹(ϑ_i − T̄_i):
⟨ϑ_i − T̄_i, G⁻¹ϑ̇_i⟩ ≤ ⟨ϑ_i − T̄_i, ∇g_t(ϑ_i)⟩ − σ_m‖ϑ_i − T̄_i‖²_{G⁻¹}. Average and use Cauchy–Schwarz for the first term:
(1/m)Σ⟨ϑ_i − T̄_i, ∇g_t⟩ ≤ V_bar·Γ ≤ (2Z)^{1/2}Γ. The Pythagorean identity W² = V_bar² + Var_γ follows from
Σ_jγ_{ij}(ϑ⁰_j − T̄_i) = 0: for each i, Σ_jγ_{ij}‖ϑ_i − ϑ⁰_j‖² = Σ_jγ_{ij}(‖ϑ_i − T̄_i‖² + ‖T̄_i − ϑ⁰_j‖²). For the pure
leak, Z is non-increasing and D⁺Z ≤ −σ_mV_bar(t)², with V_bar computed from the plan in use (the coupling inequality
above holds for every optimal plan, vertex or not), so Z converges and ∫₀^∞V_bar²dt ≤ Z(0)/σ_m < ∞. V_bar itself is
neither monotone nor continuous — it jumps at plan switches, upward as well as downward (numerics below) — so the
ω-limit statement rests on the stability of the LP instead. Let ϑ_∞ be an ω-limit point and F_∞ the set of optimal
plans at ϑ_∞ (for a vertex-selection solution, its finitely many vertices). For a fixed plan γ the barycentres T̄^γ_i do
not depend on ϑ, so V_bar^γ(ϑ) := ((1/m)Σ‖ϑ_i − T̄^γ_i‖²_{G⁻¹})^{1/2} is 1-Lipschitz in ϑ for the norm
((1/m)Σ‖·‖²_{G⁻¹})^{1/2}; and the optimal face is upper semicontinuous — the plan costs are continuous in ϑ and the
polytope has finitely many vertices — so every plan optimal at a ϑ close enough to ϑ_∞ lies in F_∞. Suppose no plan of
F_∞ is centroidal; then δ := min_{γ∈F_∞}V_bar^γ(ϑ_∞) > 0 (a minimum over a compact set), and for ‖ϑ − ϑ_∞‖ < δ/2,
V_bar(ϑ) ≥ δ/2 whatever optimal plan is in use. The speed is bounded by σ_m·diam_{G⁻¹}(Θ_K) (Lemma P.1), so the solution,
which enters the ball of radius δ/4 around ϑ_∞ at arbitrarily large times, either stays in the ball of radius δ/2 from
some time on or spends a time ≥ δ/(4σ_m diam) in it on each of infinitely many visits; either way ∫V_bar²dt = ∞, a
contradiction. Hence some plan of F_∞ — a vertex for vertex-selection solutions — is centroidal at ϑ_∞. (A two-way tie
between vertex plans γ and ρ is not a sliding mode of the vertex law unless T̄^γ = T̄^ρ: the rate of change of c_ρ − c_γ
when following γ, plus that of c_γ − c_ρ when following ρ, equals (2σ_m/m)Σ_i‖T̄^γ_i − T̄^ρ_i‖²_{G⁻¹} ≥ 0, so at least one of
the two plans stays optimal when followed. Ties between three or more plans are not analysed.) The gradient identity is the envelope theorem for the
transportation LP at a configuration with a unique optimal plan: ∂_{ϑ_i}[½Σγ_{ij}‖ϑ_i − ϑ⁰_j‖²_{G⁻¹}] = Σ_jγ_{ij}G⁻¹(ϑ_i − ϑ⁰_j)
= (1/m)G⁻¹(ϑ_i − T̄_i), so −σ_m(ϑ_i − T̄_i) = −σ_m m G∂_{ϑ_i}[½W²]. (Where the plan is not unique, ½W²(·, ν) is not differentiable; the generalised velocities are then −σ_m m G ξ with ξ
ranging over the barycentric projections of the optimal plans — the elements of the (sub)differential of ½W²(·, ν) at
an atomic μ, cf. [AGS, §10.4].)

(ii) Relabel so that (1/m)Σ‖ϑ_i − ϑ*_i‖²_{G⁻¹} = η². For any plan ρ with row sums 1/m, the plan cost
c_ρ(ϑ) := Σρ_{ij}‖ϑ_i − ϑ⁰_j‖²_{G⁻¹} is the squared norm of (ϑ_i − ϑ⁰_j) in L²(ρ), so |√c_ρ(ϑ) − √c_ρ(ϑ*)| ≤ η by the
triangle inequality (the labels of ϑ and ϑ* are matched). The barycentric part is dominated by the plan cost:
C_bar = (1/m)Σ‖ϑ*_i − T̄_i‖² ≤ c_γ(ϑ*) because Σ_jγ_{ij}‖ϑ*_i − ϑ⁰_j‖² = (1/m)‖ϑ*_i − T̄_i‖² + Σ_jγ_{ij}‖T̄_i − ϑ⁰_j‖². Let
γ* be optimal for (ϑ*, ν). Then √C_bar ≤ √c_γ(ϑ*) ≤ √c_γ(ϑ) + η ≤ √c_{γ*}(ϑ) + η ≤ √c_{γ*}(ϑ*) + 2η = b(t) + 2η. For the
equality case: c_γ(ϑ*) − b(t)² ≤ (b + 2η)² − b² = 4η(b + η) < Δ_gap(t); the vertex γ (a vertex by the choice of the law) has cost at ϑ* strictly below the next vertex cost, so it is optimal at
ϑ*; if γ*(t) is unique, γ = γ*(t) and C_bar = (1/m)Σ‖ϑ*_i − T̄*_i‖² = b_bar(t)². [For a generalised solution whose plan
at t is not a vertex, the same chain applies to each vertex of the optimal face at ϑ, so every such vertex is optimal at
ϑ*, and uniqueness of γ*(t) again gives γ = γ*(t).]

(iii) Apply Theorem B with y = ∇g_t − σ_mG⁻¹(ϑ_i − T̄_i) in (P). The extra term is
−σ_m(1/m)Σ⟨ϑ̃_i, ϑ̃_i + ϑ*_{π*(i)} − T̄_i⟩_{G⁻¹} ≤ −σ_mη² + σ_mη√C_bar. Under the gap condition √C_bar = b_bar(t) ≤ b_bar and
σ_mηb_bar ≤ (σ_m/2)(η² + b_bar²).

(iv) is the proof of Theorem C(iv) with b_bar in place of b.

(v) b_bar(t)² = W(μ*_t, ν)² − Var_{γ*(t)} ≤ W(μ*_t, ν)², with equality iff Var_{γ*} = 0, that is, iff every cell of γ*
consists of a single atom of ν (coincident atoms counted once) — the case m₀ = m, but also m₀ | m and priors with repeated
atoms — and b_bar = 0 iff ϑ*_j = T̄*_j for all j (at a tie, with the plan γ*(t) fixed by the selection). ∎

*Remark C′.1 (reading).* The barycentric leak replaces "pull each kernel to its prior kernel" by "pull each kernel to
the centre of mass of the prior it is responsible for". Three consequences.

- *Transfer across sizes.* A network of m kernels can be anchored to a learned measure of m₀ ≠ m kernels. The bias in
  the Lyapunov inequality is b_bar, the distance of the ideal from the centroids of its cells, not the quantisation error
  W(μ*, ν) itself.
- *Quantisation on the fly.* With e ≡ 0 the law is Lloyd's algorithm for the equal-mass quantisation of ν in the G⁻¹
  metric, in continuous time — the continuous-time Lloyd descent of coverage control [S18], with the equal-mass
  constraint of capacity-constrained quantisation [S19]; so the controller's leak and the quantisation theory of
  Corollary B.6 (§8, branch 9) are the same object. Like Lloyd's algorithm it converges to a critical point of the
  quantisation energy, not necessarily the global minimiser.
- *Theorem F.* The leak clauses (i′)–(iii′) of Theorem F hold with W(μ_r, ν) replaced by V_bar(r) ≤ W(μ_r, ν), since
  ‖ϑ̇_i‖_{G⁻¹} ≤ ‖∇g_t‖_G + σ_m‖ϑ_i − T̄_i‖_{G⁻¹} by Lemma P.1 and Minkowski.

*Checked* (`tests/test_barycentric.py`): the plan reduces to the assignment for m₀ = m (1e-10) and the Pythagorean
split W² = V_bar² + Var_γ holds (1e-10); along the pure leak W² is non-increasing at every step, the Dini bound
D⁺Z ≤ −σ_mV_bar² holds at every state (within 5% of the first-order Euler quotient), the flow reaches V_bar² < 1e-8, and
the end point is a fixed point of the equal-mass Lloyd iteration with the same W²; √C_bar ≤ b + 2η on 90 random
configurations, with C_bar = b_bar² at small η; the Lyapunov inequality (iii) holds at every recorded state of a closed
loop with the barycentric leak toward a 24-atom prior. The numbers are in `results/barycentric_leak/table.md`
(`exp_barycentric_leak.py`).

*Numerics* (`results/barycentric_leak/table.md`).

- *Pure leak toward a 48-atom prior* (m = 6, lifted dimension 3, G = diag(2, 1, 0.5), σ_m = 1, three random starts). W²(μ_t, ν)
  is non-increasing at every one of the 3000 steps; the flow stops where the centroid defect is at the rounding level
  (V_bar² ≤ 2.2·10⁻²⁵) and W² = Var_γ (1.164, 1.213, 1.203 from 2.55, 2.30, 2.58), after 2–6 plan switches, at which
  V_bar jumps (upward as well as downward). The end points are fixed points of the equal-mass Lloyd iteration; the best
  of 20 Lloyd runs from random starts reaches 1.171, below two of the three end points and above the third: like Lloyd's
  algorithm, the leak finds a local minimiser of the quantisation energy, not always the global one. With m₀ = m the same code reproduces the exact decay of Theorem C(i) to 5·10⁻¹⁶.
- *Closed loop* (representable plant with an m* = 4 ideal, m = 4, k_e = 4, σ_m = 0.5, three seeds). Three leaks of
  the same gain: (A) the OT leak toward the ideal itself, (B) the barycentric leak toward a 32-atom cloud of jittered
  copies of the ideal atoms, (C) the labelled leak toward a random 4-subset of the cloud. The network converges, in
  each case, to the configuration its leak anchors: (A) to the ideal (W(μ_T, μ*) and RMS‖e‖ below the printed
  precision); (B) to the
  equal-mass centroids of the cloud, at W(μ_T, μ*) = 0.167, 0.097, 0.157, which equal the centroid defects b_bar of the
  ideal (0.168, 0.097, 0.159) although W(μ*, ν) = 0.45, 0.33, 0.40; (C) to the subset, at W(μ_T, μ*) = 0.61, 1.02, 0.35
  against W(μ*, subset) = 0.65, 1.03, 0.40. The steady tracking RMS is 0.002–0.007 for (B) and 0.03–0.04 for (C). So
  the loop settles at the centroidal configuration of its prior, at distance b_bar from the ideal, as the fixed-point
  statements (i) and (v) describe; the Lyapunov bound (iv) does not apply to these runs (scope paragraph: the exchange
  bound on Δ_gap and 4η(0)(b + η(0)) are in the table). The comparison with (C) shows the cost of a prior of the right
  size at the wrong positions, not an advantage of barycentric over labelled leaks as such: a labelled leak toward the
  cloud's own equal-mass m-quantiser would have the same bias b_bar by (v).

---

## 5. Theorem D — the mean-field closed loop

The state is (x, μ) ∈ ℝⁿ × P(Θ_K), with

    ẋ = f(x,t) + ẋ_m − k_e(x − x_m) − ∫Φ(x;ϑ)dμ_t,    ∂_tμ_t + ∇·(μ_t v_t) = 0,    v_t = Proj_K(ϑ, G∇_ϑ[(x−x_m)ᵀΦ(x;ϑ)]).

Solutions are Lagrangian: μ_t = (Ξ_t)#μ_0 for the flow Ξ_t of v_t. This is law (L). The OT leak is multivalued at ties
and is excluded here.

**Theorem D.** Assume (S0)–(S2).

(a) *Local well-posedness.* For every (x_0, μ_0) ∈ ℝⁿ × P(Θ_K) there is a unique Lagrangian solution. It continues for
as long as x stays bounded.

(b) *Adaptation-free tube (m-uniform, μ*-free).* For every μ ∈ P(Θ_K), ‖∫Φ(x;·)dμ‖ ≤ a_max. For F ≥ 0 let

    r(F) := ‖e(0)‖ + (F + a_max)/k_e,    S(F) := sup{‖f(x,t)‖ : t ≥ 0, ‖x − x_m(t)‖ ≤ r(F)}.

Suppose some F ≥ 0 satisfies S(F) ≤ F. Then every solution, for every m and every μ_0:

- exists for all t ≥ 0;
- satisfies ‖x(t) − x_m(t)‖ ≤ r(F);
- satisfies ‖e(t)‖ ≤ e^{−k_e t}‖e(0)‖ + (1 − e^{−k_e t})(F + a_max)/k_e.

Examples:

- For bounded f, F = sup‖f‖ works for every e(0).
- If ‖f(x,t)‖ ≤ α + β‖x‖ with β < k_e, then F = [α + β(sup‖x_m‖ + ‖e(0)‖ + a_max/k_e)]/(1 − β/k_e) works.
- The condition can fail for superlinear f. For the Van der Pol plant of the numerics, S(F) − F ≥ 51.6 for every F.

(c) *Dobrushin estimate.* Let two solutions of the same loop (same f, x_m, k_e, G, Θ_K and ε_K) remain in a compact X
on [0,T]. Then

    ‖x_t − x'_t‖ + W(μ_t,μ'_t) ≤ e^{Λt}(‖x_0 − x'_0‖ + W(μ_0,μ'_0)),    Λ := max(L_ϑ + L_2, L_x + M_x),

and Λ is independent of m. The constants are:

- **L_ϑ.** A Lipschitz constant of ϑ ↦ v in the G⁻¹ norm, uniform over x ∈ X and t ≤ T. It is proportional to
  Ē := sup_{X×[0,T]}‖x − x_m(t)‖, and it is O(Ē/ε_K) because of the ramp (h' ≤ 2/ε_K).
- **L_x.** A Lipschitz constant, Euclidean to G⁻¹, of x ↦ v(ϑ; x, x − x_m(t)), uniform over ϑ ∈ Θ_K, t ≤ T and
  x ∈ X.
- **M_x.** M_x := sup_{t≤T}Lip_X f(·,t) + k_e + a_max J_c. Here a_max J_c bounds ‖∇_x∫Φ(x;·)dμ‖ over x ∈ ℝⁿ and
  μ ∈ P(Θ_K).
- **L_2.** L_2 := λ_G^{1/2}L₁, with L₁ from Lemma B.2, is a W-Lipschitz constant of μ ↦ ∫Φ(x;·)dμ. For any coupling,
  ‖∫Φ(x;·)d(μ − μ')‖ ≤ L_2(∫‖ϑ − ϑ'‖²_{G⁻¹}dγ)^{1/2}.

(d) *Mean-field limit.* Let (x, μ) be the solution from (x_0, μ_0), and let T be smaller than its existence time. Put
X_1 := {y : dist(y, x([0,T])) ≤ 1}, and let Λ_1 be the Λ of (c) computed on X_1 × [0,T]. Suppose x^m_0 → x_0 and
μ^m_0 → μ_0 in W, and set δ_m := ‖x^m_0 − x_0‖ + W(μ^m_0, μ_0). Then, for all large m:

- the m-particle loop exists on [0,T] and stays in X_1;
- sup_{t≤T}(‖x^m_t − x_t‖ + W(μ^m_t, μ_t)) ≤ e^{Λ_1T}δ_m → 0.

Under the hypothesis of (b), for example for bounded f, the limit is global and T is arbitrary.

(e) *Passage of (B) to the limit.* Work in the setting of (d), with X := X_1. Let μ*^m_t = (1/m)Σ_jδ_{ϑ*^m_j(t)}
satisfy (S4), and let ε^m_max ≥ sup_{X×[0,T]}‖f − ∫Φdμ*^m_t‖. Assume μ*^m_t → μ*_t uniformly in W on [0,T],
sup_m ε^m_max ≤ ε_max and sup_m V^m_* ≤ V_*. Only tails matter, so sup_m may be replaced by limsup_m. Then:

- μ*_t is V_*-Lipschitz in W, with values in P(Θ_K^{ε});
- the limit satisfies (B) in integrated form, and in Dini form for every t ∈ [0,T), with constants
  (k_e, ε_max, K, L₀, L₁, λ_G, V_*).

V_* depends on the approximating ideals; the intrinsic version is Remark B.7.

*Proof.*

(a) Write the loop as an ODE for (x, Ξ) ∈ ℝⁿ × C(Θ_K; ℝ^p), with the sup norm (and G⁻¹ on ℝ^p):

    ẋ = f(x,t) + ẋ_m − k_e(x − x_m) − ∫Φ(x; Ξ(ϑ))dμ_0(ϑ),    ∂_tΞ(ϑ) = v(Ξ(ϑ); x, t),    Ξ_0 = id,

so that μ_t = Ξ_t#μ_0. Extend v Lipschitz to a neighbourhood of Θ_K inside {s > s_min/2} (McShane). On {‖x‖ ≤ R},
with the constants of (c), the right-hand side is Lipschitz:

- ‖v(Ξ(·);x,t) − v(Ξ'(·);x',t)‖_∞ ≤ L_ϑ‖Ξ − Ξ'‖_∞ + L_x‖x − x'‖;
- the x-field changes by at most M_x‖x − x'‖ + L_2‖Ξ − Ξ'‖_∞, because W(Ξ#μ_0, Ξ'#μ_0) ≤ ‖Ξ − Ξ'‖_{L²(μ_0;G⁻¹)}.

Picard–Lindelöf gives a unique local solution. Each characteristic stays in Θ_K (Nagumo, as in A(ii)), so the solution
continues while x is bounded. Golse [82] Thms 1.3.2 and 1.4.3, and [77]–[79], are structural analogues (pairwise C¹
kernels, W₁), not the proof.

(b) While ‖e‖ ≤ r(F), ė = −k_e e + (f − f̂) with ‖f − f̂‖ ≤ F + a_max. Variation of constants gives
‖e(t)‖ ≤ e^{−k_e t}‖e(0)‖ + (1 − e^{−k_e t})(F + a_max)/k_e < r(F) for every finite t, since a_max > 0. So x never
reaches the boundary of the tube (first exit), stays bounded, and (a) continues the solution.

(c) Let γ₀ be an optimal plan for (μ_0, μ'_0), let Ξ_t and Ξ'_t be the two flows, and set
D(t)² := ∫‖Ξ_t(ϑ) − Ξ'_t(ϑ')‖²_{G⁻¹}dγ₀. Then (Ξ_t, Ξ'_t)#γ₀ couples μ_t and μ'_t, so W(μ_t, μ'_t) ≤ D(t). Split

    v_t(Ξ) − v'_t(Ξ') = [v_t(Ξ) − v_t(Ξ')] + [v_t(Ξ') − v'_t(Ξ')]

and apply Cauchy–Schwarz:

    d/dt D² ≤ 2L_ϑD² + 2L_x‖x − x'‖D,    D⁺‖x − x'‖ ≤ M_x‖x − x'‖ + L_2D.

Apply Grönwall to D_δ + ‖x − x'‖, with D_δ := (D² + δ²)^{1/2}, and let δ ↓ 0.

(d) Let τ_m ≤ T be the first time x^m leaves X_1, with τ_m = T if it never leaves. On [0, τ_m] both solutions lie in
X_1, so (c) gives ‖x^m_t − x_t‖ ≤ e^{Λ_1T}δ_m < 1 once m is large. Then x^m(τ_m) is interior to X_1, so τ_m = T. By
(a), the m-solution continues while it is bounded.

(e) Theorem B for each m gives V^m(t₂) − V^m(t₁) ≤ ∫R^m ≤ ∫R̄^m. Here R̄^m is R^m with ε^m_max and V^m_*(t) replaced by
ε_max and V_*. By (d) and the triangle inequality for W, e^m → e and W(μ^m_t, μ*^m_t) → W(μ_t, μ*_t) uniformly. So
V^m → V and R̄^m → R̄ uniformly, and V(t₂) − V(t₁) ≤ ∫R̄. R̄ is continuous, which gives the Dini form. The identity
coupling gives W(μ*^m_t, μ*^m_s) ≤ V_*|t − s|. This passes to the uniform limit, and Θ_K^{ε} is closed. ∎

*Remark D.1 (a one-sided constant free of ε_K).* The coupling step uses v only through
⟨Ξ − Ξ', G⁻¹(v(Ξ) − v(Ξ'))⟩. So the one-sided constant

    L_ϑ⁺ := sup ⟨ϑ − ϑ', G⁻¹(v(ϑ) − v(ϑ'))⟩/‖ϑ − ϑ'‖²_{G⁻¹}

may replace L_ϑ. Let N(θ,y) := h(‖θ‖)(θ̂ᵀy)₊θ̂ be the norm-ball correction, with r = ‖θ‖, c = θ̂ᵀθ̂' and
α = (θ̂ᵀy)₊. Then

    ⟨θ − θ', N(θ,y) − N(θ',y)⟩ ≥ −‖y‖‖θ − θ'‖²/(2(θ_max − ε_K)).

This follows from three facts:

- h is monotone;
- |α − α'| ≤ ‖y‖√(2(1 − c));
- ‖θ − θ'‖² ≥ 2√(rr')|r − r'|√(2(1 − c)).

The box correction is monotone. Hence L_ϑ⁺ stays bounded as ε_K → 0, and so does Λ := max(L_ϑ⁺ + L_2, L_x + M_x).
On the experiment's box, the round-2 reviewer measured L_ϑ⁺ ≈ 88–98 for ε_K from 0.2 to 0.005. Over the same range
the two-sided constant grows like 30/ε_K. This measurement is not repeated in the repository's tests.

*Numerics* (`results/mean_field_m`). Halving W(μ^m_0, μ^{2m}_0) halves W(μ^m_T, μ^{2m}_T), with a nearly constant
factor (1.64–1.69 over the four refinements). sup_t‖e^m − e^{2m}‖ falls from 4.6e-3 (16 against 32 particles) to below 1e-4 (128 against 256).

This is a consistency check only.

- On the experiment's box the guaranteed factor e^{ΛT} is at least e^{496}, from L_2 = 33.08 alone over T = 15.
- The condition of (b) fails for this Van der Pol plant, so the run is covered by the localised (d) through its
  observed boundedness, not by (b).

---

## 5b. Theorem F — the transport cost of adaptation

**Theorem F.** Assume (S0)–(S3), and take t₀ ≤ t₁ in the maximal interval of existence (all t ≥ 0 under the
hypothesis of Theorem D(b)). μ*, (S4) and (S5) are not needed. Along (L):

- (i) ((1/m)Σ‖ϑ̇_i‖²_{G⁻¹})^{1/2} ≤ λ_G^{1/2}L₁‖e(t)‖;
- (ii) W(μ_{t₁},μ_{t₀}) ≤ λ_G^{1/2}L₁∫_{t₀}^{t₁}‖e‖dr. Hence the W-length of r ↦ μ_r on [t₀,t₁] is at most
  λ_G^{1/2}L₁∫_{t₀}^{t₁}‖e‖dr;
- (iii) W(μ_{T₂},μ_{T₁}) ≤ λ_G^{1/2}L₁(T₂ − T₁)·RMS_{[T₁,T₂]}‖e‖.

Along (L_σ) with the OT leak, for a.e. t:

- (i′) ((1/m)Σ‖ϑ̇_i‖²_{G⁻¹})^{1/2} ≤ λ_G^{1/2}L₁‖e‖ + σ_mW(μ_t,μ⁰);
- (ii′) W(μ_{t₁},μ_{t₀}) ≤ ∫_{t₀}^{t₁}(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰))dr;
- (iii′) W(μ_{T₂},μ_{T₁}) ≤ λ_G^{1/2}L₁(T₂−T₁)·RMS_{[T₁,T₂]}‖e‖ + σ_m∫_{T₁}^{T₂}W(μ_r,μ⁰)dr.

For the labelled leak, replace W(μ_r,μ⁰) by d_lab(r) := ((1/m)Σ‖ϑ_i(r) − ϑ⁰_i‖²_{G⁻¹})^{1/2} ≥ W(μ_r,μ⁰). The term
W(μ_r,μ⁰) does not suffice there.

All bounds also hold with L_G := (g_a + a_max²(g_cJ_c² + g_sJ_s²))^{1/2} ≤ λ_G^{1/2}L₁ in place of λ_G^{1/2}L₁. Since
max_i‖ϑ̇_i‖_{G⁻¹} ≤ L_G‖e‖ along (L), also W_∞(μ_{t₁},μ_{t₀}) ≤ L_G∫_{t₀}^{t₁}‖e‖dr.

*Proof.*

(i) By Lemma P.1 and Lemma B.2(b),

    ‖ϑ̇_i‖_{G⁻¹} ≤ ‖∇g_t(ϑ_i)‖_G ≤ λ_G^{1/2}‖∇g_t(ϑ_i)‖ ≤ λ_G^{1/2}L₁‖e(t)‖.

Square and average over i.

(ii) Use the identity coupling (ϑ_i(t₀), ϑ_i(t₁)) and Minkowski's inequality in L²((1/m)Σ; G⁻¹). For general
measures, use the Lagrangian coupling (X_{t₀}, X_{t₁})#μ_0, or AGS Thm 8.3.1 after the isometry ϑ ↦ G^{−1/2}ϑ. Summing
over a partition gives the length statement.

(iii) Cauchy–Schwarz.

(i′) Apply Lemma P.1 with y = ∇g_t − σ_mG⁻¹(ϑ_i − ϑ⁰_{τ(i)}). Then Minkowski's inequality gives
‖y‖_G ≤ ‖∇g_t‖_G + σ_m‖ϑ_i − ϑ⁰_{τ(i)}‖_{G⁻¹}. Average over i with τ optimal. At ties, Filippov velocities are convex
combinations of fields with τ optimal at the current point. Each such field satisfies (i′), and the norm is convex.

(ii′) and (iii′) follow as (ii) and (iii).

For L_G: ‖∇g‖²_G = g_a‖φe‖² + (aᵀe)²(g_c‖∇_cφ‖² + g_s(∂_sφ)²) ≤ ‖e‖²L_G². ∎

*Reading.* Along (L), moving the kernel measure by W costs at least W/(λ_G^{1/2}L₁) of integrated tracking error. This
is an exchange rate between learning and tracking.

Along a leaky law the proved statement is ∫_{t₀}^{t₁}(λ_G^{1/2}L₁‖e‖ + σ_m d_r)dr ≥ W(μ_{t₁},μ_{t₀}), with
d_r = W(μ_r,μ⁰) for (L_σ) and d_r = d_lab(r) for the labelled leak. This gives no lower bound on ∫‖e‖, because the
leak can pay for the whole motion. Under (L_σ) there are runs with e ≡ 0 and W(μ_10, μ_0) ≈ 1.05 (review round 2).

Checked at every sample (with the particle path length ≥ the Wasserstein length):

- (i) and Lemma P.1, including a tight instance with ratio ≥ 0.99;
- (ii) along (L);
- (ii′) along (L_σ), on a grid of pairs, including an instance with e ≡ 0 where it is tight (ratio ≥ 0.99) and the
  (L) reading fails.

## 5c. Theorem D′ — propagation of chaos for the dithered closed loop on finite horizons

Theorem D covers the deterministic projected law (L). The route to Conjecture H.6 (§6c) lists the mean-field limit of
the dithered law (L_τ) as one of its three missing pieces. This section proves it on finite horizons, with an explicit
m^{−1/2} rate that does not depend on the temperature, by a synchronous coupling. The uniform-in-time statement is
§8, branch 3.

**Setting.** The design of Theorem H: saturated features Φ̃, ‖Φ̃‖ ≤ K_Φ, with all ϑ-derivatives bounded uniformly in x,
and the m-particle closed loop

    ẋ = f(x) + ẋ_m − k_e e − (1/m)Σ_jΦ̃(x;ϑ_j),    e = x − x_m(t),
    dϑ_i = b(x, e, ϑ_i)dt + √(2τ)G^{1/2}dB_i,    b(x, e, ϑ) := G∇_ϑ[eᵀΦ̃(x;ϑ)] − σ(ϑ − ϑ⁰),    i = 1..m,          (L_τ)

with independent Brownian motions B_i, i.i.d. initial particles ϑ_i(0) ∼ μ̄_0 ∈ P₂(ℝ^p) independent of the B_i, and a
deterministic x(0) = x̄_0. f is C¹ and x_m ∈ C¹ with ẋ_m bounded; no autonomy or periodicity is needed. The
**McKean–Vlasov closed loop** is

    ẋ̄ = f(x̄) + ẋ_m − k_e ē − ∫Φ̃(x̄;ϑ)dμ̄_t(ϑ),    ē = x̄ − x_m(t),
    dϑ̄ = b(x̄_t, ē_t, ϑ̄)dt + √(2τ)G^{1/2}dB,    μ̄_t := Law(ϑ̄_t),                                                 (MV)

from the same x̄_0 and ϑ̄(0) ∼ μ̄_0. Here x̄_t is deterministic, and the particle ϑ̄ interacts with its own law only
through x̄_t.

**Standing bound.** Fix T > 0 and suppose both x^m (for every m and every realisation) and x̄ remain in a compact set
X on [0,T]. Since ‖(1/m)ΣΦ̃‖ ≤ K_Φ and ‖∫Φ̃dμ̄‖ ≤ K_Φ whatever the parameters, the tube bound of Theorem D(b) (with
a_max replaced by K_Φ) gives such an X for every T when its self-consistent growth condition holds, for example for
bounded f; otherwise the statement is local in time, as in Theorem D(d). Put E_X := sup_{t≤T, x∈X}‖x − x_m(t)‖ and let

    L_f := Lip_X f,    L_x := sup‖∂_xΦ̃‖,    L_ϑ := sup‖∇_ϑΦ̃‖,
    L_b := λ_G( L_ϑ + E_X sup‖∂_x∇_ϑΦ̃‖ + E_X sup‖∇²_ϑΦ̃‖ ) + σ    (a Lipschitz constant of (x, e, ϑ) ↦ b on X × {‖e‖ ≤ E_X} × ℝ^p
                                                                      for the norm ‖x‖ + ‖e‖ + ‖ϑ‖),
    Λ := max{ L_f + L_x + 2L_b, L_ϑ + L_b }.

All of these are finite by saturation and independent of m and τ; the suprema of the Φ̃-derivatives are over all x ∈ ℝⁿ
(the saturated features and their derivatives are globally bounded), so the segments between x and x̄ used below may
leave X. The constants grow with the output dimension n (K_Φ = a_max√n), not with the lifted dimension p beyond that.

**Theorem D′.** Under the standing bound:

(a) *Well-posedness.* (MV) has a unique solution (x̄, Law(ϑ̄)) on [0,T], and x̄ is deterministic.

(b) *Synchronous coupling.* Let ϑ̄_i, i = 1..m, be the solutions of (MV)'s particle equation driven by the same B_i and
started at the same ϑ_i(0) as the m-particle system. Then, pathwise,

    ‖x_t − x̄_t‖ + (1/m)Σ_i‖ϑ_i(t) − ϑ̄_i(t)‖ ≤ e^{Λt}∫₀^t ζ_s ds,    ζ_s := ‖(1/m)Σ_jΦ̃(x̄_s;ϑ̄_j(s)) − ∫Φ̃(x̄_s;·)dμ̄_s‖,

and since the ϑ̄_j(s) are i.i.d. with law μ̄_s and x̄_s is deterministic, E ζ_s² ≤ K_Φ²/m. Hence

    E sup_{t≤T}( ‖x_t − x̄_t‖ + (1/m)Σ_i‖ϑ_i(t) − ϑ̄_i(t)‖ ) ≤ T K_Φ e^{ΛT} m^{−1/2},
    E sup_{t≤T}( ‖x_t − x̄_t‖ + (1/m)Σ_i‖ϑ_i(t) − ϑ̄_i(t)‖ )² ≤ T² K_Φ² e^{2ΛT} m^{−1},

with constants independent of τ ≥ 0 and of μ̄_0. The exponent −1/2 does not depend on the dimension; the constants do,
through K_Φ and the Lipschitz bounds of Φ̃ (standing bound).

(c) *Consequences.*
- Propagation of chaos: for each fixed k, the law of (ϑ_1, …, ϑ_k) on [0,T] converges to the k-fold product of the
  path law of ϑ̄, with the per-particle bound E sup_{t≤T}‖ϑ_i(t) − ϑ̄_i(t)‖ ≤ 2(e^{L_bT} − 1)T K_Φ e^{ΛT}m^{−1/2} for
  every i. (Exchangeability alone gives only the fixed-t bound E‖ϑ_i(t) − ϑ̄_i(t)‖ = EΔ̄(t) ≤ T K_Φ e^{ΛT}m^{−1/2}, since
  sup_t of a mean is not the mean of the sup_t; the per-particle sup bound is the pathwise Grönwall estimate of the proof.)
- The learned function converges at the rate m^{−1/2}: at each fixed t,
  E‖f̂^m(x_t) − ∫Φ̃(x̄_t;·)dμ̄_t‖ ≤ (L_x + L_ϑ)T K_Φ e^{ΛT}m^{−1/2} + E ζ_t ≤ ((L_x + L_ϑ)T e^{ΛT} + 1)K_Φ m^{−1/2}
  (the supremum over t would need a chaining argument for ζ and is not pursued).
- The parameter measure converges in W₁ (Euclidean cost on ℝ^p; the G⁻¹ metric used elsewhere changes it by a factor
  λ_G^{±1/2}) at the i.i.d. rate: E W₁(μ^m_t, μ̄_t) ≤ T K_Φ e^{ΛT}m^{−1/2} + E W₁(μ̄^{emp,m}_t, μ̄_t), where μ̄^{emp,m}_t is
  the empirical measure of m i.i.d. samples of the time-t marginal μ̄_t; the last term is O(m^{−1/p}) for p > 2
  (Fournier–Guillin; applicable because μ̄_t has a bounded second moment on [0,T], the drift being linearly bounded). Only
  this term carries the dimension.

(d) *τ = 0.* The statement holds for the deterministic saturated law with i.i.d. initial particles, where it gives, for
this law (not the projected law of Theorem D(d)), the explicit m^{−1/2} rate in place of a rate δ_m → 0 for the state and
the paired particles.

*Proof.* (a) Let X′ be a compact neighbourhood of X on which f is Lipschitz. For a given curve x̄ ∈ C([0,T₁]; X′), the
particle equation of (MV) is an SDE with drift Lipschitz in ϑ (bounded ∇²_ϑΦ̃ and the linear leak; the derivative bounds
of Φ̃ are global in x) and constant diffusion; it has a unique strong solution, whose law μ̄_t[x̄] depends
Lipschitz-continuously on x̄: by synchronous coupling of two solutions driven by x̄ and x̄′,
d‖ϑ̄ − ϑ̄′‖/dt ≤ L_b(2‖x̄ − x̄′‖ + ‖ϑ̄ − ϑ̄′‖), so W₁(μ̄_t[x̄], μ̄_t[x̄′]) ≤ 2L_b∫₀^t e^{L_b(t−s)}‖x̄_s − x̄′_s‖ds. The plant
input ẋ_m − k_eē − ∫Φ̃dμ̄ is bounded by sup‖ẋ_m‖ + k_eE_{X′} + K_Φ whatever the measure, so for T₁ small the plant solutions
from x̄_0 stay in X′, and the map x̄ ↦ (solution of the plant ODE with input ∫Φ̃(·;ϑ)dμ̄_t[x̄]) sends C([0,T₁]; X′) into
itself and is a contraction there, since μ ↦ ∫Φ̃(x;·)dμ is W₁-Lipschitz with constant L_ϑ and f is Lipschitz on X′. This
gives a unique solution on [0,T₁] in the class of curves with values in X′; the standing bound keeps it in X, so the
local solutions concatenate to [0,T], and uniqueness holds among solutions with x̄ ∈ C([0,T]; X). This is the Picard
argument of Theorem D(a) with the Dobrushin estimate replaced by the coupling estimate. (Theorem D(b)'s tube bound, used
in the standing bound for a law it does not formally cover, rests only on ‖f̂‖ ≤ a_max — here K_Φ — and on the pathwise
error ODE, so it applies.)

(b) Set Δ_i := ϑ_i − ϑ̄_i and δ := x − x̄. The noise cancels in dΔ_i, so Δ_i is absolutely continuous with
d‖Δ_i‖/dt ≤ ‖b(x,e,ϑ_i) − b(x̄,ē,ϑ̄_i)‖ ≤ L_b(‖x − x̄‖ + ‖e − ē‖ + ‖Δ_i‖) = L_b(2‖δ‖ + ‖Δ_i‖), using e − ē = δ (same
x_m). For δ,

    δ̇ = f(x) − f(x̄) − k_eδ − [(1/m)Σ_jΦ̃(x;ϑ_j) − ∫Φ̃(x̄;·)dμ̄_t],

and the bracket splits as (1/m)Σ_j[Φ̃(x;ϑ_j) − Φ̃(x̄;ϑ_j)] + (1/m)Σ_j[Φ̃(x̄;ϑ_j) − Φ̃(x̄;ϑ̄_j)] + [(1/m)Σ_jΦ̃(x̄;ϑ̄_j) − ∫Φ̃(x̄;·)dμ̄_t],
of norms at most L_x‖δ‖, L_ϑ(1/m)Σ‖Δ_j‖ and ζ_t. Dropping the dissipative −k_eδ,
d‖δ‖/dt ≤ (L_f + L_x)‖δ‖ + L_ϑΔ̄ + ζ_t with Δ̄ := (1/m)Σ‖Δ_j‖, and dΔ̄/dt ≤ L_b(2‖δ‖ + Δ̄). For D := ‖δ‖ + Δ̄ this gives
Ḋ ≤ ΛD + ζ_t with D(0) = 0, and Grönwall gives the pathwise bound, whose right-hand side is non-decreasing in t. For
the expectation, x̄_s is deterministic and the Φ̃(x̄_s;ϑ̄_j(s)) are i.i.d. bounded by K_Φ, so
E ζ_s² = (1/m)E‖Φ̃(x̄_s;ϑ̄_1(s)) − EΦ̃(x̄_s;ϑ̄_1(s))‖² ≤ K_Φ²/m, and E∫₀^Tζ_s ds ≤ T K_Φ m^{−1/2}; for the second moment,
E(∫₀^Tζ)² ≤ T∫₀^T Eζ² ≤ T²K_Φ²/m.

(c) For the first item, d‖Δ_i‖/dt ≤ L_b(2‖δ‖ + ‖Δ_i‖) with ‖δ(s)‖ ≤ D(s) ≤ e^{Λs}∫₀^sζ and Δ_i(0) = 0 give, by Grönwall
and Λ ≥ L_b, ‖Δ_i(t)‖ ≤ 2L_b∫₀^t e^{L_b(t−s)}e^{Λs}(∫₀^sζ)ds ≤ 2(e^{L_bt} − 1)e^{Λt}∫₀^tζ, a pathwise bound whose
right-hand side is non-decreasing in t; take sup_t and expectations. Convergence of the k-particle law follows from the
pathwise coupling bound. For the second, ‖f̂^m(x_t) − ∫Φ̃(x̄_t;·)dμ̄_t‖ is bounded by the same three
terms as the bracket in (b), the first two by (L_x + L_ϑ)D(t) and the third by ζ_t. For the third, couple μ^m_t with
the empirical measure of the ϑ̄_j(t) through the identity pairing (cost (1/m)Σ‖Δ_j‖) and use the triangle inequality.

(d) Nothing in the argument uses τ > 0: the noise cancels in the coupling and only the i.i.d. sampling of the initial
particles enters through ζ. ∎

*Remark D′.1 (what carries the dimension).* The plant and the learned function see the particle measure only through
∫Φ̃dμ, the mean of a bounded function. Their fluctuation is K_Φ/√m: the exponent does not see the dimension p of the
lifted parameter, and the constant K_Φ = a_max√n grows with the output dimension only.
The curse of dimensionality of Corollary B.6 (rates m^{−1/p}) concerns the Wasserstein distance between the parameter
measures, which the closed loop never evaluates. Theorem F, by contrast, bounds the transport of μ by the tracking error
and is dimension-free as well.

*Remark D′.2 (what is still missing for Conjecture H.6).* The constant e^{ΛT} grows with the horizon; uniformity in
time would follow from a contraction of the coupled system, which the plant coupling and the non-reciprocal averaged
interaction (Proposition H.7) obstruct for the general law. For the averaged system, Remark H.5 points to the
uniform-in-time result of [S5]; with the plant in the loop the question is §8, branch 3.

*Numerics* (`exp_chaos.py` → `results/chaos/table.md`). An M = 4096 particle system stands in for (MV); for
m = 16, 32, …, 1024 the m-particle system is driven by the first m initial particles and Brownian motions of the
reference — the coupling of (b) — with n = 1, k_e = 10, σ = 0.02, T = 20, at τ = 0 and τ = 2·10⁻³, and
sup_t D_m, D_m := ‖x^m − x^M‖ + (1/m)Σ_{i≤m}‖ϑ^m_i − ϑ^M_i‖, is averaged over 24 seeds. Because the m-system shares its
particles with the reference, the sampling fluctuation that drives D_m has variance proportional to 1/m − 1/M (the shared
particles contribute (1/m − 1/M)Σ_{i≤m}, the others (1/M)Σ_{i>m}), so the predicted law for the mean is a constant times
(1/m − 1/M)^{1/2}, and m^{−1/2} for m ≪ M.

| | τ = 0 | τ = 2·10⁻³ |
|---|---|---|
| mean sup_t D_m, m = 16 → 1024 | 1.49·10⁻¹ → 2.17·10⁻² | 1.52·10⁻¹ → 2.32·10⁻² |
| log-log slope of the mean against m over m ≤ 256 (bootstrap s.e. over seeds) | −0.47 ± 0.06 | −0.46 ± 0.05 |
| slope of the median | −0.46 | −0.45 |
| slope of the mean against (1/m − 1/M)^{1/2}, all m (predicted 1) | 0.94 ± 0.08 | 0.92 ± 0.08 |
| mean / (1/m − 1/M)^{1/2}, m = 32 … 1024 | 0.70–0.80 | 0.67–0.86 |

The rate is the m^{−1/2} of (b) within one standard error, the same at both temperatures (the noise cancels in the
coupling, (d)), and the two-scale law in (1/m − 1/M)^{1/2} holds up to m = M/4. Two features of the data are worth
recording. First, single seeds differ by more than an order of magnitude at fixed m (the table lists the minimum and
maximum over seeds): sup_t D_m is governed by the projection of the initial sampling fluctuation on a few slowly
amplified directions of the closed loop, so the expectation in (b) needs tens of seeds — a first version of this
experiment with three seeds gave slopes −0.25 and −0.21 and was withdrawn (`review/REVIEW_LOG.md`, round 5, D2). Second,
the measured ratio sup_t D_m/(1/m − 1/M)^{1/2} ≈ 0.7 is astronomically smaller than the constant T K_Φ e^{ΛT} of (b): on
this horizon the theorem is a rate statement, not a quantitative bound.

*Checked* (`exp_chaos.py`): the experiment only; there is no test for Theorem D′ (its hypotheses are structural and the
coupling is the experiment itself).

**Theorem D″ (uniform-in-time propagation of chaos in the leak-dominated regime).** In the setting and with the constants of
Theorem D′, assume the standing bound for all t ≥ 0 and let L₂ := sup‖∇²_ϑΦ̃‖, L_{xϑ} := sup‖∂_x∇_ϑΦ̃‖ (over all x and ϑ; for
n > 1 the operator norms are contracted with unit vectors as in §1 (1)), and

    κ_e := k_e − L_f − L_x,    κ_ϑ := σ − λ_G E_X L₂,    L_b′ := λ_G(L_ϑ + E_X L_{xϑ}),    det A := κ_eκ_ϑ − L_b′L_ϑ.

If κ_e > 0 and det A > 0 (then κ_ϑ > 0 as well) — the leak dominates the loop gain of the adaptive law — then, under the
synchronous coupling of Theorem D′(b), for every t ≥ 0 and every τ ≥ 0,

    E‖x_t − x̄_t‖ ≤ K_Φ κ_ϑ m^{−1/2}/det A,    E(1/m)Σ_i‖ϑ_i(t) − ϑ̄_i(t)‖ = E‖ϑ_1(t) − ϑ̄_1(t)‖ ≤ K_Φ L_b′ m^{−1/2}/det A,

with constants independent of t: propagation of chaos is uniform in time, at the rate m^{−1/2}. Moreover
sup_t E W₁(μ^m_t, μ̄_t) ≤ K_Φ L_b′ m^{−1/2}/det A + sup_t E W₁(μ̄^{emp,m}_t, μ̄_t), and the second term is finite: sup_t∫‖ϑ − ϑ⁰‖²dμ̄_t ≤
max{∫‖ϑ − ϑ⁰‖²dμ̄_0, ((λ_G E_X L_ϑ)²/σ + 2τ tr G)/σ}.

*Status.* Proved; checked numerically on one slow-adaptation instance (`exp_chaos_uniform.py`), where it also shows its
limit. The contraction condition κ_ϑ > 0 makes the particle map a contraction uniformly in t, so in the certified regime the
particles of a single system stay together: at τ = 0 they synchronise exponentially and the network collapses to one kernel;
at τ > 0 they form the dither cloud of the prior, N(ϑ⁰, (τ/σ)G) to leading order, and the network's output stays at the
level of its prior (Remark D″.1). In both cases the leak dominates the learning: the certified regime is the one in which
the network does not learn. What the theorem adds to Theorem D′ is the removal of the horizon from the bound in that regime,
with explicit constants (non-trivially at τ > 0, where the coupled distance is of order m^{−1/2} and not zero), and a precise
statement of where the contraction method stops. Referee round 6 (`review/REVIEW_round6_Dpp.md`) confirmed the proof,
corrected a double-counted factor 2 in the interconnection gain L_b′ (sharp constants below), and found that the run first
reported as "learning regime" was a slow collapse to the dead kernel of an odd plant; the instance was changed, and the new
long-horizon run in the learning regime shows the coupled distance growing in time — the uniform statement is supported only
where the theorem proves it.

*Proof.* Keep the contracting terms that Theorem D′'s proof dropped. For δ := x − x̄, the term −k_eδ gives
D⁺‖δ‖ ≤ −(k_e − L_f − L_x)‖δ‖ + L_ϑΔ̄ + ζ_t with the same decomposition of the bracket as in D′(b). For Δ_i := ϑ_i − ϑ̄_i,
split b(x,e,ϑ_i) − b(x̄,ē,ϑ̄_i) = [b(x,e,ϑ_i) − b(x̄,ē,ϑ_i)] + [b(x̄,ē,ϑ_i) − b(x̄,ē,ϑ̄_i)]. Since e − ē = δ, the first bracket is
G∇_ϑΦ̃(x;ϑ_i)ᵀ(e − ē) + G[∇_ϑΦ̃(x;ϑ_i) − ∇_ϑΦ̃(x̄;ϑ_i)]ᵀē, of norm at most λ_G L_ϑ‖δ‖ + λ_G L_{xϑ}‖δ‖E_X = L_b′‖δ‖ on the tube.
The second bracket is G[∇_ϑΦ̃(x̄;ϑ_i) − ∇_ϑΦ̃(x̄;ϑ̄_i)]ᵀē − σΔ_i, and ⟨Δ_i, ·⟩ ≤ λ_G‖ē‖L₂‖Δ_i‖² − σ‖Δ_i‖² ≤ −κ_ϑ‖Δ_i‖², so
D⁺‖Δ_i‖ ≤ −κ_ϑ‖Δ_i‖ + L_b′‖δ‖, and the same for the average Δ̄. The pair u := (‖δ‖, Δ̄) thus satisfies
D⁺u ≤ Au + (ζ_t, 0)ᵀ componentwise, u(0) = 0, with the Metzler matrix A := [[−κ_e, L_ϑ], [L_b′, −κ_ϑ]] (non-negative
off-diagonal entries), whose flow e^{As} is entrywise non-negative. The comparison principle for cooperative systems gives
u(t) ≤ ∫₀^t e^{A(t−s)}(ζ_s, 0)ᵀds entrywise. A is Hurwitz exactly when its trace is negative and its determinant positive;
with κ_e > 0 and det A = κ_eκ_ϑ − L_b′L_ϑ > 0 one has κ_ϑ > L_b′L_ϑ/κ_e ≥ 0, so both hold, and then
∫₀^∞e^{As}ds = −A⁻¹ = (det A)⁻¹[[κ_ϑ, L_ϑ], [L_b′, κ_e]], entrywise non-negative. Taking expectations, with Eζ_s ≤
(Eζ_s²)^{1/2} ≤ K_Φ m^{−1/2} for every s (the ϑ̄_j(s) are i.i.d. with law μ̄_s and x̄_s is deterministic, as in D′(b)),
E u(t) ≤ K_Φ m^{−1/2}(−A⁻¹)e₁ = K_Φ m^{−1/2}(κ_ϑ, L_b′)ᵀ/det A for every t. The per-particle statement is exchangeability
at fixed t (valid here, unlike for sup_t: Theorem D′(c)); the W₁ statement is the triangle inequality as in D′(c). For the
second moment of μ̄_t: z(t) := E‖ϑ̄_t − ϑ⁰‖² satisfies, by Itô's formula with ‖G∇_ϑ(ēᵀΦ̃)‖ ≤ λ_G E_X L_ϑ =: β,
ż ≤ −2σz + 2β√z + 2τ tr G ≤ −σz + β²/σ + 2τ tr G, whence z(t) ≤ max{z(0), (β²/σ + 2τ tr G)/σ}. The noise cancels in Δ_i, so
the coupling estimate does not depend on τ. ∎

*Remark D″.1 (reading).*

- *The certified regime is the leak-dominated one.* The condition κ_ϑ = σ − λ_G E_X L₂ > 0 says that ϑ ↦ b(x̄, ē, ϑ) is a
  contraction uniformly in t. At τ = 0 the same estimate applied to two particles of the *same* system gives
  D⁺‖ϑ_i − ϑ_j‖ ≤ −κ_ϑ‖ϑ_i − ϑ_j‖: all particles synchronise exponentially, the network collapses to a single kernel, and
  the m-system and the M-system — which share their first m particles — become identical, so the coupled distance is
  not just bounded but zero (the numerics: particle spread and coupled distance at the rounding level). At τ > 0 the two
  Brownian increments do not cancel in ϑ_i − ϑ_j; the particles form a cloud around the common mean motion whose
  second moment is bounded as in the proof, N(ϑ⁰, (τ/σ)G) to leading order when the leak dominates (measured: rms spread
  within 3% of (τ tr G/σ)^{1/2} in the referee's runs), the coupled distance is of order m^{−1/2} and the bound is not
  trivial — but the network's output stays at the level of its prior: sup|f̂| is a few percent of sup|f| in the certified
  runs. Either way, propagation of chaos is uniform in time because the leak has stopped the learning. This is the price
  of the contraction method, not of the closed loop: the learning regime, in which kernels separate and the network fits
  the plant (σ below the critical leak of Proposition K at the collapsed critical point, the leak below the separating
  force of the residual), is where Theorems H, I and K live, where the averaged interaction is non-reciprocal, and where
  uniformity in time remains the open part of §8, branch 3.
- *What the conditions say in the scalings of the program.* On the a-priori tube of Theorem D(b) with e(0) = 0,
  k_eE_X = F + K_Φ, so the certified leak is σ_cert := λ_G E_X L₂ + L_b′L_ϑ/κ_e = [λ_G(F + K_Φ)L₂ + λ_G L_ϑ²]/k_e + O(k_e⁻²):
  every fixed leak gain σ > 0 certifies uniform-in-time propagation of chaos for k_e large, provided the tracking error
  starts at zero (or ‖e(0)‖ = O(1/k_e); otherwise E_X ≥ ‖e(0)‖ and the condition does not improve with k_e) — and stops
  the learning. In the scaling σ = σ₀/k_e of Theorem J the condition reads σ₀ > λ_G(F + K_Φ)L₂ + λ_G L_ϑ², which implies
  that J_m is strongly convex in the G⁻¹ metric (its Hessian is bounded below by σ₀ − λ_G sup_Θ‖r_Θ‖L₂ ≥ σ₀ − λ_G(F + K_Φ)L₂
  there) with a unique, collapsed minimiser.
- *Why the constants are conservative.* L_ϑ, L₂, L_{xϑ} are suprema over the parameter box and the tube; for sharp
  kernels (small s_min) they are large, so σ_cert can exceed the adaptation gains by an order of magnitude; the a-priori
  E_X is itself several times the measured sup‖ē‖ (0.23 against 0.07 in the numerics below). The gap to the learning regime
  is not a matter of constants, however: in the numerics the coupled distance grows over the horizon once the kernels
  separate.
- *Technique.* The comparison system is the standard contraction argument for McKean–Vlasov systems with a confining
  potential, with the plant as a second contracting channel (rate κ_e) and the two interconnection gains L_ϑ (particles →
  plant) and L_b′ (plant → particles); det A > 0 is the small-gain condition of that interconnection.

*Numerics* (`exp_chaos_uniform.py` → `results/chaos/uniform_table.md`). A slow-adaptation instance (n = 1, k_e = 25,
G = 0.05·diag(4,1,1), saturated features with a_max = 4, s ≥ 0.3, prior ϑ⁰ = 0), a two-bump plant that is not odd,
f(x) = 1.6e^{−(x−1)²/0.4} − 1.2e^{−(x+1)²/0.4} + 0.1x, on x_m = 1.5 sin t + 0.5 sin 2t (so that the dead kernel is not a
critical point; referee round 6), m = 16, 64, 256 coupled to an M = 2048 reference. The constants are exact: the Lipschitz
constants of the saturated features in closed form, L_x = L_ϑ = a_max e^{−1/2}/s_min = 8.087 and L₂ = L_{xϑ} = a_max/s_min² =
44.4 (a numerical search over 1.2·10⁵ parameter samples with Nelder–Mead polish reproduces them to 10⁻⁸), L_f = 2.29, and the
a-priori tube radius of Theorem D(b) with e(0) = 0, E_X = (F + K_Φ)/k_e = 0.228 (F = 1.70; the measured sup‖ē‖ is 0.070):
κ_e = 14.6, L_b′ = 3.64, σ_cert = 4.04. Proposition K gives the critical leak of the collapsed critical point (the best single
kernel, a = 0.136, c = 0.068, s = 1.15): σ_K = σ₀_crit/k_e = 0.00144 — the certified leak is 2800 times the leak below which
the kernels separate. Four runs:

| σ, τ | regime | T, seeds | reference network at T: rms spread in (a, c, s) (from 1.53), fit error ‖f − f̂‖/‖f‖ (from 1.00) | max D_m, first → last window, m = 16 and 256 | slope against m (per-seed, ± s.e.) | bound sup_t E D_m, m = 16 / 256 |
|---|---|---|---|---|---|---|
| 4.85 (1.2σ_cert), 0.02 | certified, with dither | 80 s, 8 | 0.033 (the dither cloud: (τ tr G/σ)^{1/2} = 0.035), 1.000 | 6.7·10⁻⁴ → 7.9·10⁻⁴ and 1.6·10⁻⁴ → 1.9·10⁻⁴ | −0.51 ± 0.02 | 0.55 / 0.14 |
| 4.85, 0 | certified, no dither | 80 s, 8 | 1.7·10⁻¹⁴, 1.000 | 2.5·10⁻²⁵ → 0 and 1.2·10⁻²⁶ → 0 | — | 0.55 / 0.14 |
| 0.505 (σ_cert/8), 0 | leak-dominated, below the threshold | 80 s, 8 | 3.2·10⁻¹⁴, 1.000 | 7.2·10⁻⁵ → 4.3·10⁻¹¹ and 3.1·10⁻⁵ → 1.8·10⁻¹¹ | −0.33 ± 0.20 | — |
| 0.00048 (σ_K/3), 0 | learning regime | 1000 s, 6 | 1.53, 0.632 | 6.1·10⁻² → 2.0·10⁻¹ and 1.3·10⁻² → 3.9·10⁻² | −0.62 ± 0.16 | — |

(Windows (10, 20], (20, 40], (40, 80] for the 80 s runs and (100, 200], (200, 500], (500, 1000] for the long run; the table
reports E max_t D_m over a window, which exceeds the theorem's sup_t E D_m.) At the certified leak with dither the
particles stay in the cloud of the prior, the network's output stays at the level of the empty network (fit error 1.000),
and the coupled distance is of order m^{−1/2} (slope −0.51 ± 0.02), bounded in time (the last window, twice as long, has
maxima 15–18% above the first) and seven hundred times below the bound: the theorem's content at τ > 0. Without dither
the reference collapses (spread 10⁻¹⁴) and the coupled distance is zero to rounding; at σ_cert/8 the particles still
synchronise and the distance decays with them. In the learning regime (σ_K/3, the collapsed configuration unstable) the
network stays spread, fits the plant (fit error 1.00 → 0.63 at T = 1000 s with these slow gains) and the coupled distance
follows the m^{−1/2} law at each time but *grows* over the horizon — by a factor 3 between the windows (100, 200] and
(500, 1000] for every m: on this horizon propagation of chaos in the learning regime is the finite-horizon statement of
Theorem D′, not a uniform one. Whether the growth saturates on longer horizons (the leak's own time scale is 1/σ = 2000 s)
is open; the earlier claim that the distance stays bounded in a "spread" regime (version 4, first draft) rested on a run
that was in fact a slow collapse, and is withdrawn (referee round 6, DPP-3).

*Checked* (`exp_chaos_uniform.py`): the experiment only.

---

## 5d. Theorem D‴ — uniform-in-time propagation of chaos in the learning regime (the averaged slow loop)

Theorem D″ ends where the leak stops dominating the loop gain: in the learning regime its small-gain argument has nothing
to say, and the long-horizon run of `exp_chaos_uniform.py` showed the coupled distance growing over the horizon. This
section isolates the structure that the small-gain argument throws away. Under the synchronous coupling the interaction
between the particles enters through the fitted function f̂_Θ, which is linear in the empirical measure; the loss L̄ is a
convex quadratic of it (flat convexity, Theorem H(i)); and the interaction term of the coupling identity is therefore the
negative of a squared L²(ρ)-norm at linear order — it contracts, whatever the gains. What remains and can expand is the
self term, whose curvature is the Hessian of the first-variation potential, σ₀G⁻¹ − R[μ](ϑ) with R the *residual-weighted*
curvature of the features — the same object as Proposition K's R₂ at a collapsed configuration. Proposition K says that
the collapsed configuration is unstable when this curvature exceeds the leak *with the collapsed residual*; at a fitted
configuration the residual is the approximation error, and the second-order optimality condition of the limit measure
says that the curvature is *below* the leak at every atom. The theorem below turns this into a uniform-in-time statement
for the averaged slow loop at τ = 0, with explicit constants, in the regime in which the network learns.

**Setting.** The averaged slow dynamics of §6f at τ₀ = 0 (Theorem J′'s leading order L̄₀ without dither; slow time s = t/k_e,
leak σ₀ = k_eσ, f = f(x) continuous, ρ the occupation measure of the T₀-periodic reference, F_ρ := sup_{supp ρ}‖f‖,
saturated features, K_Φ := sup‖Φ̃‖), written as the G-preconditioned gradient flow of J_m(Θ) = mL̄(μ_Θ) + (σ₀/2)Σ_i‖ϑ_i − ϑ⁰‖²_{G⁻¹}:

    dϑ_i/ds = −G∇_{ϑ_i}J_m(Θ) = −G∇_ϑ𝒱[μ_Θ](ϑ_i),    𝒱[μ](ϑ) := V_L[μ](ϑ) + (σ₀/2)‖ϑ − ϑ⁰‖²_{G⁻¹},    i = 1..m,          (𝔄_m)

with V_L[μ](ϑ) = −⟨r_μᵀΦ̃(·;ϑ)⟩_ρ and r_μ = f − f̂[μ], f̂[μ] := ∫Φ̃(·;ϑ′)dμ(ϑ′), as in §6c. Its mean-field equation is the
Wasserstein gradient flow (𝔄) of §6c at τ = 0 in slow time, ∂_sμ_s = ∇·(μ_sG∇𝒱[μ_s]), solved by the characteristic flow

    dX_s(ϑ)/ds = −G∇𝒱[μ_s](X_s(ϑ)),    X_0 = id,    μ_s = (X_s)_#μ₀                                                     (𝔄₀)

(a McKean–Vlasov equation without noise whose vector field has linear growth through the leak, is Lipschitz in ϑ and, through
r_μ, Lipschitz in μ for W₁; existence, uniqueness and the representation by characteristics follow from Dobrushin's argument
[S35] — the synchronous-coupling proof of Theorem D′(b) with the plant removed — and the velocity is bounded on the invariant
ball of Step 0). The constants are taken in the G⁻¹ metric on parameters ("whitened"):

    L_ϑ^G := sup_{x,ϑ}‖∇_ϑΦ̃(x;ϑ)G^{1/2}‖,    L₂^G := sup_{x,ϑ,‖u‖=1}‖G^{1/2}∇²_ϑ(uᵀΦ̃(x;ϑ))G^{1/2}‖,

so that ‖Φ̃(x;ϑ) − Φ̃(x;ϑ′)‖ ≤ L_ϑ^G‖ϑ − ϑ′‖_{G⁻¹}, ‖∇_ϑΦ̃(x;ϑ)v‖ ≤ L_ϑ^G‖v‖_{G⁻¹} and ‖(∇_ϑΦ̃(x;ϑ) − ∇_ϑΦ̃(x;ϑ′))v‖ ≤
L₂^G‖ϑ − ϑ′‖_{G⁻¹}‖v‖_{G⁻¹} (L_ϑ^G ≤ λ_G^{1/2}L_ϑ and L₂^G ≤ λ_GL₂ in the constants of §5c; `chaos_learning.whitened_lipschitz_constants`
computes the suprema by search, with closed-form bounds). The *residual-weighted curvature* of μ at ϑ is the symmetric matrix

    R[μ](ϑ) := ⟨Σ_k r_μ^k(x)∇²_ϑΦ̃^k(x;ϑ)⟩_ρ = −∇²_ϑV_L[μ](ϑ),    R^G := G^{1/2}RG^{1/2},    κ_R[μ](A) := sup_{ϑ∈A}λ_max(R^G[μ](ϑ)),

so that G^{1/2}∇²_ϑ𝒱[μ](ϑ)G^{1/2} = σ₀I − R^G[μ](ϑ); for μ = δ_ϑ̄ and ϑ = ϑ̄ it is Proposition K's
R₂ and σ_crit(ϑ̄) = κ_R[δ_ϑ̄]({ϑ̄}). It satisfies ‖R^G[μ](ϑ)‖ ≤ L₂^G⟨‖r_μ‖⟩_ρ ≤ L₂^G(F_ρ + K_Φ) and ‖R^G[μ](ϑ) − R^G[ν](ϑ)‖ ≤
L₂^G‖f̂[μ] − f̂[ν]‖_{L¹(ρ)}. N_ϱ(A) := {ϑ : dist_{G⁻¹}(ϑ, A) ≤ ϱ}.

*The coupling.* ϑ_i(0), i = 1..m, i.i.d. with law μ₀ of compact support, R₀ := sup_{supp μ₀}‖ϑ − ϑ⁰‖_{G⁻¹}; Θ(s) solves (𝔄_m) and
the reference particles are ϑ̄_i(s) := X_s(ϑ_i(0)), i.i.d. with law μ_s and ϑ̄_i(s) ∈ supp μ_s. Let δ_i := ϑ_i − ϑ̄_i,

    D(s) := (1/m)Σ_i‖δ_i(s)‖²_{G⁻¹},    Y(s) := ‖f̂_{Θ̄(s)} − f̂[μ_s]‖_{L²(ρ)},

Y being the sampling fluctuation of the reference's fitted function (E Y(s)² = (1/m)⟨Var_{μ_s}Φ̃(x;·)⟩_ρ ≤ K_Φ²/m). A priori
(Step 0 below) every particle of either system stays in the ball ‖ϑ − ϑ⁰‖_{G⁻¹} ≤ R_Θ := max{R₀, L_ϑ^G(F_ρ + K_Φ)/σ₀}, so D ≤
D_max := 4R_Θ², and the characteristics move at speed at most V_max := L_ϑ^G(F_ρ + K_Φ) + σ₀R_Θ.

**Theorem D‴ (uniform-in-time propagation of chaos in the learning regime, τ = 0).** Assume

(H) *(non-degenerate stability along the mean-field flow)* there are s₁ ≥ 0, κ > 0 and ϱ′ > 0 such that

        κ_R[μ_s](N_{ϱ′}(supp μ_s)) ≤ σ₀ − κ    for every s ≥ s₁,

and put a(s) := σ₀ − κ_R[μ_s](N_{ϱ′}(supp μ_s)) (so a ≥ κ on [s₁, ∞) and a ≥ σ₀ − L₂^G(F_ρ + K_Φ) everywhere) and
Λ₁ := ∫₀^{s₁}(κ − a(s))₊ds < ∞. Define

    ϱ := min{(2κ)^{1/2}/L₂^G, ϱ′κe^{−Λ₁}/(4(L_ϑ^G)²)},    θ := min{κ/(8L₂^G), ϱκe^{−Λ₁}/(4L_ϑ^G), ϱ′κe^{−Λ₁}/(4L_ϑ^G)},
    Δ := θ/(4L_ϑ^GV_max),    C₁ := 4e^{2Λ₁}(L_ϑ^G)²K_Φ²/κ²,    m₀ := 16K_Φ²/θ².

Then for every m ≥ m₀:

(a) *(exponentially long horizons)* for every S > 0 and every s ∈ [0, S],

        E D(s) ≤ C₁/m + 4R_Θ²(S/Δ + 2)exp(−mθ²/(32K_Φ²));

    in particular sup_{s ≤ S}E D(s) ≤ (C₁ + 12R_Θ²)/m for all S ≤ Δe^{mθ²/(64K_Φ²)} and m ≥ m₀ with e^{mθ²/(64K_Φ²)} ≥ m.

(b) *(uniform in time)* if moreover (H′) the characteristic flow converges exponentially on the support of μ₀,
    sup_{ϑ∈supp μ₀}‖X_s(ϑ) − X_∞(ϑ)‖_{G⁻¹} ≤ C₂e^{−κ₂s} for all s ≥ 0, then for every s ≥ 0

        E D(s) ≤ C₁/m + 4R_Θ²(R_eff/Δ + 3)exp(−mθ²/(32K_Φ²)),    R_eff := κ₂⁻¹(log(8L_ϑ^GC₂/θ))₊.

(c) *(the learned function and the particles)* in either case, with the same right-hand side B(s) for E D(s),

        E‖f̂_{Θ(s)} − f̂[μ_s]‖²_{L²(ρ)} ≤ 2(L_ϑ^G)²B(s) + 2K_Φ²/m,    E W₂(μ_{Θ(s)}, μ_s) ≤ B(s)^{1/2} + E W₂(μ_{Θ̄(s)}, μ_s),

    and on the event 𝒢_S := {sup_{s≤S}Y(s) ≤ θ} (case (a); 𝒢_∞ in case (b)), whose complement is the tail above, every
    particle stays within 3ϱ′/4 of its reference up to time S: max_i sup_{s≤S}‖δ_i(s)‖_{G⁻¹} ≤ 3ϱ′/4 — no particle changes basin.

The constants depend on (κ, ϱ′, Λ₁, L_ϑ^G, L₂^G, K_Φ, F_ρ, σ₀, R₀) and not on m, s or S.

*Status.* Proved (synchronous coupling with the three-term identity of Step 1, a bootstrap on the mean square and the
per-particle distance, and McDiarmid's inequality [S36] for the reference's sampling fluctuation; the finite-horizon part of
the transient is Theorem D′'s mechanism and appears in Λ₁). Hypothesis (H) is a property of the deterministic mean-field
flow and can be computed; it is the strict second-order optimality condition of the limit configuration with a margin and
a neighbourhood (Remark D‴.1), and it is compatible with the learning regime σ₀ < σ_crit of Proposition K. It requires the
support of μ_s to stay away from the ridges of 𝒱[μ_s] for s ≥ s₁: for an initial cloud whose support is a union of pieces
each flowing into one basin (a warm start, or a cluster per bump) this holds; for a connected cold start whose mass splits
between several basins it fails for every s (Remark D‴.2), and the numerics below show what the coupling does there.
Numerics: `exp_chaos_learning.py` → `results/chaos/learning_table.md` (below). The theorem is about the averaged slow
dynamics; its transfer to the real loop at finite k_e is open (Remark D‴.4). Its constants are worst-case and its regime is
asymptotic: on the instance below m₀ is of order 10²³ (Remark D‴.6), and the finite-m behaviour is described sharply by
Proposition D‴.5. Referee round 8 (`review/REVIEW_round8_Dppp.md`) re-derived every constant of Steps 0–9, found no counterexample
to (2.1), (3.1), (4.1) on an adversarial search, and corrected the wording of Remarks D‴.1–D‴.2, the good event per horizon and the
well-posedness sentence (all applied). Referee round 9 (`review/REVIEW_round9_D3b.md`) reviewed Proposition D‴.5, Corollary D‴.7 and the
numerics: all hold with fixes, the main one being that the finite-m replicated configuration is a minimiser only for admissible weights
(Proposition D‴.5(iv)) — one of the four m = 16 warm seeds sits on a saddle, which the numerics paragraph now says.

*Proof.* Step 0 (a priori bounds). For either system, d/ds‖ϑ − ϑ⁰‖²_{G⁻¹} = 2⟨ϑ − ϑ⁰, ⟨∇_ϑ(rᵀΦ̃(·;ϑ))⟩_ρ⟩ − 2σ₀‖ϑ − ϑ⁰‖²_{G⁻¹} with
|⟨ϑ − ϑ⁰, ⟨∇_ϑ(rᵀΦ̃)⟩_ρ⟩| = |⟨rᵀ∇_ϑΦ̃(ϑ − ϑ⁰)⟩_ρ| ≤ ⟨‖r‖⟩_ρL_ϑ^G‖ϑ − ϑ⁰‖_{G⁻¹} ≤ (F_ρ + K_Φ)L_ϑ^G‖ϑ − ϑ⁰‖_{G⁻¹}, so
D⁺‖ϑ − ϑ⁰‖_{G⁻¹} ≤ L_ϑ^G(F_ρ + K_Φ) − σ₀‖ϑ − ϑ⁰‖_{G⁻¹} and ‖ϑ(s) − ϑ⁰‖_{G⁻¹} ≤ R_Θ; hence ‖δ_i‖_{G⁻¹} ≤ 2R_Θ, D ≤ D_max, and
‖dϑ/ds‖_{G⁻¹} = ‖G^{1/2}⟨∇_ϑ(rᵀΦ̃)⟩_ρ − σ₀G^{−1/2}(ϑ − ϑ⁰)‖ ≤ V_max. Two Lipschitz facts: for Θ, Θ̄ with rows ϑ_j, ϑ̄_j and
A_j(x) := ∫₀¹∇_ϑΦ̃(x;ϑ̄_j + tδ_j)dt,

    f̂_Θ − f̂_Θ̄ = (1/m)Σ_jA_jδ_j =: h,    ‖h‖_{L²(ρ)} ≤ L_ϑ^G(1/m)Σ_j‖δ_j‖_{G⁻¹} ≤ L_ϑ^GD^{1/2},                        (0.1)
    ‖R^G[μ_Θ̄](ϑ) − R^G[μ_s](ϑ)‖ ≤ L₂^G‖r_Θ̄ − r_{μ_s}‖_{L¹(ρ)} ≤ L₂^GY,    ‖(∇_ϑΦ̃(x;ϑ_i) − A_i(x))v‖ ≤ (L₂^G/2)‖δ_i‖_{G⁻¹}‖v‖_{G⁻¹}.   (0.2)

Step 1 (the identity). Since dδ_i/ds = −G[∇𝒱[μ_Θ](ϑ_i) − ∇𝒱[μ_s](ϑ̄_i)] and ⟨δ, Gv⟩_{G⁻¹} = ⟨δ, v⟩,

    ½dD/ds = −(1/m)Σ_i⟨δ_i, ∇𝒱[μ_Θ](ϑ_i) − ∇𝒱[μ_s](ϑ̄_i)⟩ = T_A + T_B + T_C,
    T_A := −(1/m)Σ_i⟨δ_i, ∇𝒱[μ_Θ̄](ϑ_i) − ∇𝒱[μ_Θ̄](ϑ̄_i)⟩,
    T_B := −(1/m)Σ_i⟨δ_i, ∇V_L[μ_Θ](ϑ_i) − ∇V_L[μ_Θ̄](ϑ_i)⟩,
    T_C := −(1/m)Σ_i⟨δ_i, ∇V_L[μ_Θ̄](ϑ̄_i) − ∇V_L[μ_s](ϑ̄_i)⟩

(the leak is linear and sits in T_A; the three brackets telescope to ∇𝒱[μ_Θ](ϑ_i) − ∇𝒱[μ_s](ϑ̄_i)). The self term T_A compares
the *reference's* potential at the two points; the interaction term T_B compares the two empirical measures at the m-system's
point; T_C is the reference's sampling error. `chaos_learning.AveragedSlowLoop.decomposition` evaluates the three terms and
checks the identity to rounding.

Step 2 (the self term is local convexity). ∇𝒱[μ_Θ̄](ϑ_i) − ∇𝒱[μ_Θ̄](ϑ̄_i) = ∫₀¹∇²𝒱[μ_Θ̄](ϑ̄_i + tδ_i)dt·δ_i = σ₀G⁻¹δ_i − R̄_iδ_i with
R̄_i := ∫₀¹R[μ_Θ̄](ϑ̄_i + tδ_i)dt, so T_A = −σ₀D + (1/m)Σ_i⟨δ_i, R̄_iδ_i⟩ and ⟨δ, Rδ⟩ ≤ λ_max(R^G)‖δ‖²_{G⁻¹}. If ‖δ_i‖_{G⁻¹} ≤ ϱ′ the
segment lies in N_{ϱ′}(supp μ_s) (ϑ̄_i ∈ supp μ_s), and by (0.2) λ_max(R^G[μ_Θ̄]) ≤ κ_R[μ_s](N_{ϱ′}(supp μ_s)) + L₂^GY there. Hence,
as long as max_i‖δ_i‖_{G⁻¹} ≤ ϱ′,

    T_A ≤ −(a(s) − L₂^GY(s))D.                                                                                           (2.1)

Step 3 (the interaction term is a negative square at linear order). ∇V_L[μ_Θ](ϑ_i) − ∇V_L[μ_Θ̄](ϑ_i) = ⟨∇_ϑΦ̃(·;ϑ_i)ᵀ(f̂_Θ − f̂_Θ̄)⟩_ρ =
⟨∇_ϑΦ̃(·;ϑ_i)ᵀh⟩_ρ, so, with the L²(ρ;ℝⁿ) inner product ⟨⟨·,·⟩⟩,

    T_B = −⟨⟨(1/m)Σ_i∇_ϑΦ̃(·;ϑ_i)δ_i, h⟩⟩ = −‖h‖² − ⟨⟨(1/m)Σ_i(∇_ϑΦ̃(·;ϑ_i) − A_i)δ_i, h⟩⟩ ≤ −‖h‖² + ‖h‖(L₂^G/2)D ≤ (L₂^G)²D²/16     (3.1)

by (0.2) and max_x(−x² + cx) = c²/4. The leading part −‖h‖² = −‖f̂_Θ − f̂_Θ̄‖²_{L²(ρ)} is the flat convexity of L̄: the
interaction contracts, and only the Taylor remainder of the features, quadratic in D, is charged. This is the step that
Theorem D″'s small-gain bound replaces by +L_ϑ²D-type terms.

Step 4 (the sampling term). ∇V_L[μ_Θ̄](ϑ̄_i) − ∇V_L[μ_s](ϑ̄_i) = ⟨∇_ϑΦ̃(·;ϑ̄_i)ᵀ(f̂_Θ̄ − f̂[μ_s])⟩_ρ, so |⟨δ_i, ·⟩| ≤ L_ϑ^G‖δ_i‖_{G⁻¹}Y and

    T_C ≤ L_ϑ^GD^{1/2}Y.                                                                                                 (4.1)

Step 5 (the differential inequalities). While max_i‖δ_i‖_{G⁻¹} ≤ ϱ′,

    ½dD/ds ≤ −(a − L₂^GY)D + (L₂^G)²D²/16 + L_ϑ^GD^{1/2}Y,                                                              (5.1)

and per particle, from the same three brackets with |T_B,i| ≤ ‖δ_i‖L_ϑ^G‖h‖ ≤ ‖δ_i‖(L_ϑ^G)²D^{1/2} by (0.1) and |T_C,i| ≤ ‖δ_i‖L_ϑ^GY,

    D⁺‖δ_i‖_{G⁻¹} ≤ −(a − L₂^GY)‖δ_i‖_{G⁻¹} + (L_ϑ^G)²D^{1/2} + L_ϑ^GY.                                                    (5.2)

Step 6 (the bootstrap on the good event). Fix a horizon S ≤ ∞, let 𝒢_S := {sup_{s≤S}Y(s) ≤ θ} and T* := inf{s ≤ S : D(s) ≥ ϱ²
or max_i‖δ_i(s)‖_{G⁻¹} ≥ ϱ′} (T* > 0 since δ(0) = 0; everything below is on [0, S]). On 𝒢_S and [0, T*): L₂^GY ≤ κ/8, (L₂^G)²D/16 ≤ (L₂^G)²ϱ²/16 ≤ κ/8, and L_ϑ^GD^{1/2}Y ≤ (κ/8)D + (2(L_ϑ^G)²/κ)Y²
(Young), so (5.1) gives dD/ds ≤ −2(a − 3κ/8)D + (4(L_ϑ^G)²/κ)Y². Since a − 3κ/8 ≥ κ/2 − (κ − a)₊, ∫_r^s(a − 3κ/8) ≥ (κ/2)(s − r) − Λ₁
for r ≤ s, and Grönwall gives

    D(s) ≤ Z(s) := 4e^{2Λ₁}(L_ϑ^G)²κ⁻¹∫₀^s e^{−κ(s−r)}Y(r)²dr ≤ 4e^{2Λ₁}(L_ϑ^G)²κ⁻²θ² ≤ ϱ²/4    on 𝒢_S ∩ [0, T*)             (6.1)

by the second constraint on θ. Per particle, (5.2) with D^{1/2} ≤ ϱ/2 and Y ≤ θ, and ∫_r^s(a − κ/8) ≥ (κ/2)(s − r) − Λ₁, gives
‖δ_i(s)‖_{G⁻¹} ≤ 2e^{Λ₁}κ⁻¹[(L_ϑ^G)²ϱ/2 + L_ϑ^Gθ] ≤ ϱ′/4 + ϱ′/2 = 3ϱ′/4 by the second constraint on ϱ and the third on θ. Both
thresholds are therefore strict on [0, T*), so T* = S on 𝒢_S by continuity: D ≤ Z and max_i‖δ_i‖ ≤ 3ϱ′/4 on [0, S] on 𝒢_S. Z does
not involve the m-system.

Step 7 (expectation). For s ≤ S, E D(s) ≤ E[Z(s)] + D_maxP(𝒢_S^c) ≤ 4e^{2Λ₁}(L_ϑ^G)²κ⁻¹∫₀^s e^{−κ(s−r)}E Y(r)²dr + D_maxP(𝒢_S^c) ≤
C₁/m + 4R_Θ²P(𝒢_S^c), using E Y(r)² ≤ K_Φ²/m.

Step 8 (the good event). Y(s) is a function of the i.i.d. sample (ϑ_j(0))_{j≤m}: Y(s) = ‖(1/m)Σ_j[Φ̃(·;X_s(ϑ_j(0))) − ∫Φ̃(·;X_s(ϑ))dμ₀]‖_{L²(ρ)}.
Changing one ϑ_j(0) changes it by at most 2K_Φ/m (bounded differences), and E Y(s) ≤ (E Y²)^{1/2} ≤ K_Φ/m^{1/2}; McDiarmid's
inequality [S36] gives P(Y(s) ≥ K_Φ/m^{1/2} + t) ≤ exp(−mt²/(2K_Φ²)), so for m ≥ m₀ = 16K_Φ²/θ² (K_Φ/m^{1/2} ≤ θ/4, t = θ/4)

    P(Y(s) ≥ θ/2) ≤ exp(−mθ²/(32K_Φ²))    for every fixed s.                                                                (8.1)

In time, |Y(s) − Y(s′)| ≤ 2L_ϑ^GV_max|s − s′| (each X_s(ϑ) moves at speed ≤ V_max and Φ̃ is L_ϑ^G-Lipschitz, for the sample and for
the mean), so on a grid of step Δ = θ/(4L_ϑ^GV_max), Y ≤ θ on [s_k, s_k + Δ] whenever Y(s_k) ≤ θ/2. (a) The grid {kΔ : kΔ ≤ S} ∪ {S}
has at most S/Δ + 2 points, and a union bound with (8.1) gives P(𝒢_S^c) = P(sup_{s≤S}Y > θ) ≤ (S/Δ + 2)exp(−mθ²/(32K_Φ²)); with
Step 7 this is (a). For S ≤ Δe^{mθ²/(64K_Φ²)}, (S/Δ + 2)e^{−mθ²/(32K_Φ²)} ≤ 3e^{−mθ²/(64K_Φ²)} ≤ 3/m when e^{mθ²/(64K_Φ²)} ≥ m. (b) Under (H′),
for s ≥ R_eff every X_s(ϑ) is within θ/(8L_ϑ^G) of X_∞(ϑ) on supp μ₀, so Y(s) ≤ Y(∞) + θ/4 with Y(∞) := ‖(1/m)Σ_j[Φ̃(·;X_∞(ϑ_j(0))) −
∫Φ̃(·;X_∞)dμ₀]‖, another bounded-differences statistic with E Y(∞) ≤ K_Φ/m^{1/2}; the grid on [0, R_eff] plus the point at infinity
has at most R_eff/Δ + 3 points, and P(𝒢_∞^c) = P(sup_{s≥0}Y > θ) ≤ (R_eff/Δ + 3)exp(−mθ²/(32K_Φ²)).

Step 9 (the learned function and the measures). ‖f̂_Θ − f̂[μ_s]‖ ≤ ‖h‖ + Y ≤ L_ϑ^GD^{1/2} + Y by (0.1), so E‖·‖² ≤ 2(L_ϑ^G)²E D + 2K_Φ²/m.
The coupling (ϑ_i, ϑ̄_i) is a coupling of μ_Θ and μ_Θ̄, so W₂²(μ_Θ, μ_Θ̄) ≤ D (G⁻¹ metric), and the triangle inequality gives (c).
The per-particle statement is Step 6. ∎

*Remark D‴.1 (what the hypothesis says, and the window).* (H) at s = ∞ reads: on a neighbourhood of every atom of the limit
measure, G^{1/2}∇²𝒱[μ_∞]G^{1/2} = σ₀I − R^G[μ_∞] ⪰ κI — each atom is a non-degenerate local minimiser of the first-variation
potential 𝒱[μ_∞], with margin κ: the strict form of the second-order optimality condition for measures (the local Chizat–Bach
condition, as in Proposition K(ii)). Its non-strict form, σ₀G⁻¹ − R[μ_∞](ϑ*) ⪰ 0 at every atom ϑ* of positive mass, is necessary
for an atomic stationary measure to be a local minimiser of 𝒢 in the Wasserstein sense: splitting an atom of mass w along a
direction of negative curvature of 𝒱[μ_∞] lowers 𝒢 at second order (Proposition K(i)–(ii) with the atom's own residual). Along
the flow, (H) asks for this with the
*current* residual r_{μ_s} from the time s₁ on. It does not mention the collapsed configuration: σ_crit = κ_R[δ_ϑ̄]({ϑ̄}) is
the same curvature *with the collapsed residual*, which is the residual of the best single kernel and is large when the
plant needs several. The learning regime σ₀ < σ_crit (the kernels separate) and the uniform regime (H) (κ_R[μ_∞] ≤ σ₀ − κ at
the fitted configuration) therefore coexist exactly when the fitted configuration is non-degenerately stable while the
collapsed one is not — the window

    κ_R[μ_∞](N_{ϱ′}(supp μ_∞)) < σ₀ < σ_crit(ϑ̄)

of leaks in which the network learns *and* propagation of chaos is uniform in time. Theorem D″'s condition, by contrast,
asks the leak to dominate the worst-case residual curvature λ_GL₂(F_ρ + K_Φ) *and* the interaction's Lipschitz constant
L_ϑ²-type terms, which forces σ₀ ≥ σ_crit. The two ingredients that remove the horizon here are (i) the sign of the
interaction, which is the flat convexity of L̄ (Step 3), and (ii) the residual-weighted curvature evaluated at the actual
residual rather than at the worst case (Step 2). Both are properties of the lifted linear parametrisation f̂ = ∫Φ̃dμ, not of
the particular features.

*Remark D‴.2 (the transient, the support, and the floor of the coupling).* Two things limit the theorem, and both are visible in
the numerics. (i) The constant carries e^{2Λ₁}, the finite-horizon amplification of Theorem D′ over the time the flow needs
to enter the non-degenerate regime: for a warm start Λ₁ = 0 and C₁ = 4(L_ϑ^G)²K_Φ²/κ²; for a cold start the expanding
phase (a(s) < κ while the residual is large) lasts a time of order 1/σ₀ and Λ₁ is of order σ₀⁻¹·sup(κ_R − σ₀) — the bound is
uniform but quantitatively void, as D′'s e^{ΛT} is. (ii) (H) is a statement about the ϱ′-neighbourhood of the *support* of
μ_s. The characteristic flow is a homeomorphism, so a connected support stays connected: if the mass of a connected cold
start splits between several basins of 𝒱[μ_∞], the support of μ_s threads through the ridges between them for every s,
and the filament's points converge to the saddles of 𝒱[μ_∞] on those ridges, where 𝒱 has a direction of negative curvature
and κ_R > σ₀ — (H) fails for every s₁ (asymptotically in s; pointwise the filament may well lie in a locally convex region for
a while, referee round 8, D3-11), however little mass the filament carries. The theorem then says nothing, and the coupling can
indeed carry a floor: a particle whose reference lies near a
separatrix can be pushed to the other side by the O(m^{−1/2}) difference between the two drifts, after which the pair sits
in different basins for ever, at a distance of the order of the basin separation. A fraction p_m of such pairs contributes
p_m·O(R_Θ²) to D, so D follows the law of p_m (the μ₀-mass of a boundary layer of width O(m^{−1/2}) around the separatrices'
preimage, hence p_m = O(m^{−1/2}) when μ₀ has a bounded density there), not 1/m; the learned function is less sensitive,
‖f̂_Θ − f̂_Θ̄‖ ≤ L_ϑ^G(1/m)Σ‖δ_i‖_{G⁻¹} ≤ L_ϑ^G(p_m·2R_Θ + D_good^{1/2}), which is still O(m^{−1/2}). The measure-weighted form of
(H) — a bad set of small μ_s-mass, with the per-particle bootstrap run on the good particles only and the interaction
split as in Step 3 — is the natural extension and is listed in §8; what the present theorem shows is where the clean
statement holds (supports inside basins) and why the cold start of Theorem D″'s experiment was outside it on two counts
(the transient and the connected support). (iii) The floor has a second source at small m, on a warm start: when the realised
cluster weights fall outside the admissible set of Proposition D‴.5(iv) the replicated configuration is a saddle of J_m, the m-system
lingers on it and then splits a cluster, so that pairs end in different basins although the mean-field flow satisfies (H); its
probability is exponentially small in m but appreciable at m = 16 on the instance below.

*Remark D‴.3 (dither).* At τ₀ > 0 the synchronous coupling cancels the noise and the identity of Step 1 holds pathwise with
the same three terms, but the reference law has full support and the particles visit the ridges of 𝒱 with the Gibbs weight
of the ridge, so (H) fails in the sense of (ii) above and the coupling carries a floor of that order. Uniform-in-time
propagation of chaos at τ₀ > 0 is the province of the entropy method: with dither, (𝔄_m) is the m-particle mean-field
Langevin dynamics of the free energy F_τ of Theorem H in the G⁻¹ metric (Corollary J′.3 is its stationary statement), and
[S5] proves uniform-in-time propagation of chaos for mean-field Langevin dynamics under flat convexity, regularity of the
first variation and a uniform log-Sobolev inequality for the proximal Gibbs measures — the three hypotheses that Theorem
H(i)–(ii) verifies for F_τ, with the constant ρ_τ of H(ii). That result is imported here, not re-derived; its constants are
those of Theorem H and are exponentially small in 1/τ₀, where the constants of Theorem D‴ are those of the leak and the
residual curvature. The leak-dominated regime is covered for the real loop and every τ ≥ 0 by Theorem D″.

*Remark D‴.4 (the real loop).* Theorem D‴ is a statement about the averaged slow dynamics (𝔄_m). Theorems J and J′ identify
(𝔄_m) as the leading order of the real loop's slow dynamics, with the corrections ε·b (Theorem J′'s b, not a gradient:
Proposition J′.2) and the fast error eliminated at O(ε²); the transfer of a uniform-in-time statement through these
corrections — a near-identity change of variables that is uniform in time — is the open item of §8, branch 1; on a warm start
it is settled by Corollary D‴.7 below through Theorem J. What the theorem settles is the mechanism: in the learning regime the uniformity rests on the fitted configuration's non-degenerate
stability and on the flat convexity of the loss, not on a dominant leak; the real-loop run on a warm start (numerics,
table 5) tracks the averaged slow dynamics over the horizon run, at the same k_e at which Theorem D″'s cold start grew.

**Proposition D‴.5 (the stationary coupled distance at finite m: the weight-sensitivity law).** Let μ* = Σ_{k=1}^K w̄_kδ_{ϑ*_k} be a
K-atom stationary measure of (𝔄₀) (each ϑ*_k a critical point of 𝒱[μ*]), and for weights w ∈ (0, ∞)^K let 𝒢_w(Θ) := ½⟨‖f −
Σ_kw_kΦ̃(·;ϑ_k)‖²⟩_ρ + (σ₀/2)Σ_kw_k‖ϑ_k − ϑ⁰‖²_{G⁻¹} be the mean-field energy of μ_w := Σ_kw_kδ_{ϑ_k}, so that ∇_{ϑ_k}𝒢_w(Θ) =
w_k∇_ϑ𝒱[μ_w](ϑ_k). Assume the K-atom Hessian ∇²_Θ𝒢_{w̄}(Θ*) is positive definite. Then:

(i) *(the weighted stationary point and its replication)* there is a smooth map w ↦ Θ*(w) on a neighbourhood of w̄ with Θ*(w̄) = Θ*
and ∇_ϑ𝒱[μ_{Θ*(w)}](ϑ*_k(w)) = 0 for every k (implicit function theorem); S(w) := ∂Θ*/∂w (S_k := ∂ϑ*_k/∂w ∈ ℝ^{p×K}) is given by
S = −[∇²_Θ𝒢_w]⁻¹∂_w∇_Θ𝒢_w at (Θ*(w), w). For m particles with n_k of them in cluster k (w^m_k := n_k/m), the configuration Θ*_m(w^m)
with n_k copies of ϑ*_k(w^m) is a critical point of J_m (μ_Θ = μ_{w^m} there). Its whitened Hessian splits (Proposition K(i), applied
cluster by cluster) into the separation blocks σ₀I − R^G[μ_{w^m}](ϑ*_k(w^m)), each with multiplicity n_k − 1, and the collective block
diag(w^m)^{−1/2}·G_K^{1/2}∇²_Θ𝒢_{w^m}G_K^{1/2}·diag(w^m)^{−1/2} (the K-atom Hessian, whitened and weighted). Hence Θ*_m(w^m) is a strict
local minimiser of J_m exactly when w^m lies in the *admissible set*

    W_adm := {w : σ₀G⁻¹ − R[μ_{Θ*(w)}](ϑ*_k(w)) ≻ 0 for every k and ∇²𝒢_w(Θ*(w)) ≻ 0},

a neighbourhood of w̄ (both conditions are open and hold at w̄ by (H) and the assumption), and a saddle when the splitting
condition fails at some cluster: a cluster can then split — by Proposition K's local separation inside its basin (sub-clusters along
the unstable direction, the cluster weights unchanged) or by particles migrating to other basins; the numerics show both (§5d, referee
rounds 9 and 10).

(ii) *(the stationary distance)* Under the shared-particle coupling of Theorem D‴ with i.i.d. cluster labels of law w̄ (the m
particles are the first m of the M reference particles), suppose w^m, w^M ∈ W_adm and that the m-system converges to Θ*_m(w^m) and
the reference to Θ*_M(w^M). Then the coupling distance converges to

    D_∞ = Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖²_{G⁻¹} = Σ_kw^m_k‖S_k(w^M)(w^m − w^M)‖²_{G⁻¹} + O(‖w^m − w^M‖³),

and, with S̃_k := G^{−1/2}S_k(w̄) (replacing S(w^M) by S(w̄) costs a relative O(‖w^M − w̄‖)) and E[(w^m − w^M)(w^m − w^M)ᵀ] =
(1/m − 1/M)C̃, C̃ := diag(w̄) − w̄w̄ᵀ (multinomial counts of the shared labels), conditionally on the admissibility of both weights,

    E[D_∞ | w^m, w^M ∈ W_adm] = (1/m − 1/M)·Σ_kw̄_k tr(S̃_kC̃S̃_kᵀ) + O(m^{−3/2}) + O(m^{−1}M^{−1/2}):

the 1/m law with an explicit constant, and the (1/m − 1/M) law of the shared-particle coupling observed in Theorem D′'s numerics.

(iii) *(the fluctuations)* As M → ∞ (w^M → w̄) and then m → ∞, √m(w^m − w̄) → N(0, C̃) and √m·(ϑ_i − ϑ̄_i) → S_{k(i)}ξ with ξ ~ N(0, C̃)
for every particle i of cluster k(i): the fluctuations of the coupling are Gaussian, carried entirely by the cluster weights, and
all particles of a cluster are displaced by the same vector.

(iv) *(admissibility)* Let δ_w := sup{δ : ‖w − w̄‖_∞ ≤ δ ⇒ w ∈ W_adm} > 0. By Hoeffding's inequality for each multinomial coordinate,
P(w^m ∉ W_adm) ≤ 2K exp(−2mδ_w²): the finite-m configuration is a strict local minimiser with probability exponentially close to
one, the exceptional event being the one on which the realised weights make the replicated configuration a saddle: a cluster of the
m-system then splits, locally or by migration to other basins — in the latter case with Remark D‴.2's floor on a warm start.

*Status.* Formal: (i) and the algebra of (ii)–(iv) are exact (the block splitting is Proposition K(i) cluster by cluster, confirmed
on J₁₀ to 8·10⁻¹² in referee round 9; implicit function theorem; the multinomial covariance of shared labels; Hoeffding), the
convergence of both systems to the stated configurations is assumed (it is what Theorem D‴ proves for m ≥ m₀ and what the numerics
show for m = 16–256 at admissible weights; at inadmissible weights the m-system sits on the saddle for a long time and then leaves
it, numerics below). Checked: `weighted_stationary_atoms`, `weight_sensitivity`, `predicted_stationary_D` in
`chaos_learning.py`; on the instance below the replicated weighted stationary point is a critical point of J₁₆ to 10⁻¹⁰, the
linearisation ϑ*(w̄ + dw) − ϑ* = S·dw has a second-order remainder, and the per-cluster displacements of the coupled run point
along S_k(w^m − w^M) (test and table below). The proposition is the sharp form of Theorem D‴'s constant: S = −H_w⁻¹∂_w∇𝒢_w
has the structure (stiffness)⁻¹ × (feature correlations), of which C₁ = 4e^{2Λ₁}(L_ϑ^G)²K_Φ²/κ² is the worst-case envelope
(κ⁻¹ for the stiffness, L_ϑ^GK_Φ for the correlations).

*Remark D‴.6 (the regime of the theorem, and what the numerics test).* The sample size m₀ = 16K_Φ²/θ² from which the theorem's
tail is in force is astronomically large on realistic instances (10²³ below): θ carries ϱ′κ²e^{−2Λ₁}/(L_ϑ^G)³ from the per-particle
bootstrap of Step 6, because the interaction's effect on one particle is bounded by (L_ϑ^G)²D^{1/2} and divided by the stiffness κ.
This is not an artefact of the proof's constants only: Proposition D‴.5 says the particles of a cluster are displaced by
S_k(w^m − w^M), of order ‖S‖m^{−1/2} with ‖S‖ ~ (stiffness)⁻¹, which exceeds any fixed ϱ′ for m below (‖S‖/ϱ′)² — on the
instance below the m = 16 particles sit 0.4–1.4 whitened units from their reference copies, outside the neighbourhood N_ϱ′ in
which 𝒱[μ_Θ̄] is convex, at the replicated stationary configuration Θ*₁₆(w^m) of their own J₁₆ — a strict local minimiser when
w^m ∈ W_adm and a saddle otherwise (Proposition D‴.5(i), (iv); one of the four m = 16 seeds below). Hypothesis (H) is nevertheless the right object: it is the
condition under which the reference configuration is stable to splitting, and the m-system's configuration is the same
configuration with perturbed weights. What the numerics test is therefore (a) hypothesis (H) at the stationary measure and
along the flow, (b) the mechanism — the coupled distance saturates (no growth) and follows the 1/m law — and (c) the sharp
constant of Proposition D‴.5; the theorem's own constant C₁ is reported for the record.

**Corollary D‴.7 (the real loop on a warm start: uniform in time through Theorem J).** Let μ* = Σ_kw̄_kδ_{ϑ*_k} be a K-atom
stationary measure as in Proposition D‴.5, satisfying the splitting condition σ₀G⁻¹ − R[μ*](ϑ*_k) ≻ 0 at every atom and the
non-degeneracy of the K-atom Hessian, and let w ↦ Θ*(w) be the weighted stationary point of D‴.5(i). Fix m ≤ M and cluster
weights w^m, w^M ∈ W_adm (Proposition D‴.5(i)), and let Θ*_m :=
the configuration with n_k = mw^m_k copies of ϑ*_k(w^m) (likewise Θ*_M). Both are non-degenerate local minimisers of J_m and J_M
(Proposition D‴.5(i)), hence hyperbolic critical points, and Theorem J applies to each loop: there are k₀, ρ₀, C (depending on
Θ*_m and Θ*_M, hence on m, M and the weights) such that for k_e ≥ k₀ the real closed loop (CL) with m kernels started at
‖Θ(0) − Θ*_m‖_{G⁻¹} ≤ ρ₀, ‖e(0)‖ ≤ R converges to its T₀-periodic orbit Θ_p^m, and likewise the M-kernel loop to Θ_p^M. Under the
shared-particle coupling (the m kernels are the first m of the M, the plants are driven by the same reference) the coupled distance
D(t) := (1/m)Σ_i‖ϑ_i(t) − ϑ̄_i(t)‖²_{G⁻¹} satisfies, for all t ≥ 0,

    D(t)^{1/2} ≤ D_orb^{1/2} + Ce^{−μt}(‖Θ(0) − Θ_p^m(0)‖_{G⁻¹} + ‖x(0) − x_p^m(0)‖) + Ce^{−μt}(‖Θ̄(0) − Θ_p^M(0)‖_{G⁻¹} + ‖x̄(0) − x_p^M(0)‖),

with μ = (λ_*/k_e)(1 − κ_J) the rate of Theorem J(c) (κ_J its small parameter; the smaller of the two loops' rates) and

    D_orb := sup_t (1/m)Σ_i‖Θ_p^m(t)_i − Θ_p^M(t)_i‖²_{G⁻¹} = Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖²_{G⁻¹} + O(1/k_e)
           = Σ_kw^m_k‖S_k(w^m − w^M)‖²_{G⁻¹} + O(‖w^m − w^M‖³) + O(1/k_e):

the coupled distance of the real loop is bounded uniformly in time, converges to the distance between the two periodic orbits,
and that distance is Proposition D‴.5's stationary distance up to the O(1/k_e) offsets of Theorem J(a),(d) (which differ between
the two loops by O(k_e⁻¹‖w^m − w^M‖) + O(k_e⁻²), the first-order offset of a replicated configuration being a smooth function of the
K-atom data (Θ*(w), w) by permutation equivariance). In particular, for cluster labels i.i.d. of law w̄ and conditionally on both
weights being admissible and both starts lying in their basins, E D_orb = (1/m − 1/M)Σ_kw̄_k tr(S̃_kC̃S̃_kᵀ) + O(m^{−3/2}) + O(1/k_e).

*Status.* Proved as stated (Theorem J(a),(c),(d) applied to each loop, the triangle inequality in the G⁻¹ metric, and Proposition
D‴.5(ii) for the distance between the critical points), under Theorem J's hypotheses for both configurations: the constants k₀,
ρ₀, C, μ are those of Theorem J at Θ*_m and Θ*_M, and the corollary does not assert their uniformity in m or M (Remark D‴.8). It
is the statement of Remark D‴.4 for a warm start: the real loop's coupled distance saturates, at a value set by the cluster weights.
Consistent with the instance of §5d's numerics (table 5: the real loop at k_e = 25 from the jittered mixture over 80 slow units, a sixth
of the relaxation time 1/κ, follows its averaged slow dynamics on the same initial clouds to within a few percent at every window end;
the saturation at the orbit distance lies beyond the horizon run and is not checked).

*Proof.* Theorem J(c) for the m-kernel loop gives ‖Θ(t) − Θ_p^m(t)‖_{G⁻¹} ≤ Ce^{−μt}(‖Θ(0) − Θ_p^m(0)‖_{G⁻¹} + ‖x(0) − x_p^m(0)‖) and the
same for the M-kernel loop; the first m rows of the M-loop are the reference particles ϑ̄_i. For the (m×p)-row norm
(1/m)Σ_i‖·‖²_{G⁻¹} ≤ ‖·‖²_{G⁻¹}/m ≤ ‖·‖²_{G⁻¹} (the full configuration norm dominates the mean square over the first m rows), so by
the triangle inequality D(t)^{1/2} ≤ [(1/m)Σ_i‖Θ_p^m(t)_i − Θ_p^M(t)_i‖²]^{1/2} + (the two transients). Theorem J(a) places each
orbit within C/k_e of its critical point, uniformly in t, so the first term is Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖²_{G⁻¹} up to O(1/k_e)
(the two offsets are O(1/k_e) each, with O(k_e⁻²) remainders; at a replicated configuration the first-order offset is the same vector
for all particles of a cluster and depends smoothly on the K-atom data (Θ*(w), w), by the permutation equivariance of (CL) and the block
splitting of Proposition D‴.5(i), so the two loops' offsets differ by O(k_e⁻¹‖w^m − w^M‖) + O(k_e⁻²) — the two loops live on ℝ^{mp} and
ℝ^{Mp}, and it is this reduction to the K-atom data that makes them comparable), and Proposition D‴.5(ii) expands the distance between the
critical points in w^m − w^M. The expectation statement is D‴.5(ii)'s. ∎

*Remark D‴.8 (what the corollary does and does not say).* It transfers the uniform-in-time statement to the real loop in the regime
where Theorem J applies to both loops — a warm start in the basins of two hyperbolic minimisers at a gain k_e above both thresholds.
Theorem J's constants at a replicated configuration are controlled by the K-atom data (the Hessian of J_m there is the collective
block ∇²𝒢_{w^m} ⊕ the separation blocks σ₀G⁻¹ − R_k of Proposition K(i), both independent of m beyond the weights), so the
threshold and the basin are expected to be uniform in m for configurations with the same weights; the corollary does not prove
this, and the cold start — where the coupled distance of the real loop grew in Theorem D″'s experiment — is outside it on the
same two counts as in Remark D‴.2.

**Proposition D‴.9 (the leak toward a measure-valued prior, and the design of the window).** Let ν = Σ_{j=1}^{K₀}q_jδ_{ϑ⁰_j} be a prior
measure with distinct atoms and let the leak be the G-gradient of the distance-to-support potential

    V_ν(ϑ) := ½ min_j‖ϑ − ϑ⁰_j‖²_{G⁻¹},    −G∇V_ν(ϑ) = −(ϑ − ϑ⁰_{j(ϑ)}),    j(ϑ) := argmin_j‖ϑ − ϑ⁰_j‖_{G⁻¹}

(the *nearest-atom leak*: inside each Voronoi cell of ν it is the σ-modification toward that cell's atom; the barycentric OT leak
of Theorem C′ is its mass-constrained relative). Write (𝔄_m^ν) and (𝔄₀^ν) for (𝔄_m) and (𝔄₀) with σ₀(ϑ − ϑ⁰) replaced by
σ₀(ϑ − ϑ⁰_{j(ϑ)}), and 𝒱^ν[μ] := V_L[μ] + σ₀V_ν for the first-variation potential.

(i) *(the theorems carry over)* Assume, in addition to the hypotheses of Theorem D‴ for (𝔄₀^ν), that the assignment is locally
constant along the flow: for every s ≥ 0 the neighbourhood N_{ϱ′}(supp μ_s) meets no Voronoi boundary of ν (in particular when the
gap between the distances to the two nearest atoms exceeds 2ϱ′ on supp μ_s, the G⁻¹-distance to the boundary being at least half the
gap). Then Theorem D‴, Proposition D‴.5 and Corollary D‴.7 hold for (𝔄_m^ν) with the same constants, R_Θ replaced by
max{R₀, L_ϑ^G(F_ρ + K_Φ)/σ₀} + max_j‖ϑ⁰_j − ϑ⁰_1‖_{G⁻¹} and ϑ⁰ by ϑ⁰_{j(·)} in the a priori bounds, and Theorem J's basin radius ρ₀ at the
replicated configuration taken below min_k dist_{G⁻¹}(ϑ*_k, ∂cell_k). The masses q_j enter the nearest-atom leak only through supp ν;
they matter in (ii), and they are what the mass constraint of the OT leak (Theorem C′, Remark D‴.10) acts on.

(ii) *(the stationary residual and the window)* If ν is placed on an exact representation of the plant on the orbit — K₀ = K atoms
ϑ⁰_k with Σ_kq_kΦ̃(·;ϑ⁰_k) = f on supp ρ — then μ* := Σ_kq_kδ_{ϑ⁰_k} is a stationary measure of (𝔄₀^ν) with r_{μ*} = 0, R[μ*] ≡ 0,
and the window of Remark D‴.1 has its maximal margin: κ_R[μ*](A) = 0 for every set A (R[μ*] ≡ 0) and κ = σ₀, the only restriction
on ϱ′ being (i)'s Voronoi hypothesis, against κ = σ₀ − κ_R[μ*_σ](N_{ϱ′}) for the σ-modification's stationary measure μ*_σ, whose
residual is the leak's bias. For
a prior ν_δ displaced from the exact representation the residual of the stationary measure is O(δ) and κ_R[μ*_{ν_δ}] = O(δ): the
margin decreases with the misplacement, and the window closes when the displacement makes the residual-weighted curvature at an
atom reach σ₀.

(iii) *(the floors)* Under the shared-particle coupling the stationary distances are the weight-sensitivity laws of Proposition
D‴.5 at parameter level, E[D_∞ | w^m, w^M ∈ W_adm] = (1/m − 1/M)Σ_kw̄_k tr(S̃_kC̃S̃_kᵀ) + O(m^{−3/2}) + O(m⁻¹M^{−1/2}), and at
function level,

    E‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖²_{L²(ρ)} = (1/m − 1/M)·tr(FC̃Fᵀ)/N + O(m^{−3/2}),    F_{·j} := ∂f̂_w/∂w_j = Φ̃(·;ϑ*_j) + Σ_kw_k∇_ϑΦ̃(·;ϑ*_k)S_{k·j},

with S the sensitivity of the weighted stationary point of (𝔄₀^ν) to the weights. Neither floor vanishes for the matching prior:
away from w̄ the weighted stationary point does not stay on the zero-residual manifold (the leak pulls the amplitudes back to the
prior while the fit asks for a_k ∝ 1/w_k); on the instance both sensitivities rise with the amplitudes of the stationary atoms along
the amplitude-scaled priors — a numerical trend, not a law (at a fixed amplitude the width changes the sensitivity, table 1).

*Status.* (i) proved (the proofs of Theorem D‴, Proposition D‴.5 and Corollary D‴.7 use the leak only through its curvature σ₀G⁻¹
and its a priori bound — the latter goes through with dist(·, supp ν) in place of ‖· − ϑ⁰‖, the upper Dini derivative of the minimum
of quadratics along the flow being bounded by the active branch's. The local constancy of the assignment makes every segment [ϑ̄_i, ϑ_i]
of Step 2 lie in one cell, and it is needed for every s ≥ 0, not only from s₁ on, because (2.1) and (5.2) are used on [0, s₁) as well:
a pair straddling a boundary adds (σ₀/m)⟨δ_i, ϑ⁰_{j(ϑ_i)} − ϑ⁰_{j(ϑ̄_i)}⟩_{G⁻¹} ≥ 0 to ½dD/ds (the sum of the two cell inequalities; V_ν is
semiconcave and its leak expands across the boundaries), a forcing of size up to σ₀ diam(supp ν)√D that none of the theorem's constants
absorbs — referee round 10 found `decomposition`'s same-cell convention and the true ½dD/ds to differ by exactly this term on a
straddling pair. The same hypothesis gives the well-posedness of (𝔄₀^ν), whose vector field is discontinuous on the boundaries, and
Theorem J applies to J_m^ν, which is C² away from the Voronoi boundaries, on a basin that meets none.) (ii) proved (r_{μ*} = 0 and
∇V_ν(ϑ⁰_k) = 0 give stationarity; R[μ*] = ⟨r_{μ*}∇²Φ̃⟩_ρ = 0; the O(δ) statements are the smooth dependence of the weighted stationary
point on the prior, by the implicit function theorem at μ*, whose K-atom Hessian Gram + σ₀diag(q)G⁻¹ is positive definite — on the
instance its whitened, weighted eigenvalues are 0.0074–0.20, all at or above σ₀, referee round 10). (iii) formal as Proposition D‴.5 (the
function-level law is the delta method applied to w ↦ f̂_{Θ*(w)}); the sign of the comparison with the σ-modification is a
numerical finding (`exp_prior_design.py` → `results/chaos/prior_design.md`, below), not a theorem. Checked: `tests/test_prior_design.py`
(the single-atom prior reproduces (𝔄_m); the gradient of J_m^ν against finite differences; the coupling identity inside the cells;
the matching prior's exact representation stationary with zero residual and zero curvature, margin σ₀ against the σ-modification's;
the parameter- and function-level sensitivities, the latter against finite differences of the fitted function of the weighted
stationary point).

*Remark D‴.10 (design).* The leak's prior is a design input of the window. The σ-modification's single ϑ⁰ biases every kernel toward
one point and leaves a residual at the fitted configuration — the residual-weighted curvature that bounds the window from below
and slows the relaxation (κ = 0.00196 on the instance of §5d, the theorem's relaxation time 1/κ = 511 slow units). A prior measure
whose atoms sit on the plant's structure removes the bias and gives the maximal margin κ = σ₀ (1/κ = 141 slow units, 3.6 times
shorter; the observed e-folding times of the coupled distance are 3.3 times shorter in the median, numerics below),
at the price of a larger floor at both levels (the stationary atoms have the full amplitudes, to which the weight sensitivity is
proportional): the classical trade-off between bias and variance, here between the window's margin and the coupling's floor
(Proposition D‴.24 resolves it at the function level: the floor is the weight-slope of the bias, and an adaptive prior removes both, at the
price of the collective rate, of the leak's a-priori bound — restored by projecting the prior, Proposition D‴.27 — and of a larger
parameter-level floor). A
misplaced prior does worse than no structured prior (centres off by 0.6, amplitudes ×1.2 or widths ×0.7 close the window, which the
σ-modification keeps open with margin 0.0020), and a shrunk prior does better than the matching one (amplitudes ×0.5–0.7: margin above
σ₀ at the uniform weights, both floors lower, all four tested unbalanced weights admissible — with splitting margins above σ₀ at one of
them and slightly below at the other three; numerics below and Remark D‴.11). The barycentric OT leak of Theorem C′ adds a mass constraint
to the same geometry — each prior atom must receive its mass — and so acts on the cluster weights themselves, the source of the
floor; its interplay with the weight fluctuation (an equal-mass Lloyd dynamics against the fit) is the natural next step (§8,
branch 8).

*Numerics* (`exp_prior_design.py` → `results/chaos/prior_design.md`). The three-bump instance of Theorem D‴'s numerics at the same leak
σ₀ = 0.15σ_crit = 0.00707, with the nearest-atom leak toward a prior of three atoms placed at the exact representation (3A_k, C_k, S_k) and
misplaced from it, against the σ-modification (ϑ⁰ = 0); ϱ′ = 0.1; the floors predicted for m = 64 against M = 2048.

*The design table (table 1).* Stationary three-atom measure μ*_ν under each prior — fit error, κ_R at the atoms and over N_{0.1}, margin
κ = σ₀ − κ_R(N_{0.1}), whitened ‖S_k‖, predicted E D_∞ and E‖f̂₆₄ − f̂_M‖²:

| prior | fit error | κ_R at the atoms (max) | κ_R over N_{0.1} | margin κ | ‖S_k‖ | E D_∞ (m = 64) | E‖f̂₆₄ − f̂_M‖² |
|---|---|---|---|---|---|---|---|
| σ-modification, ϑ⁰ = 0 | 0.2439 | +0.00411 | +0.00512 | +0.00196 | 8.38, 4.66, 8.32 | 2.010·10⁻¹ | 4.287·10⁻⁴ |
| matching prior (exact representation) | 0.0000 | +0.00000 | +0.00000 | +0.00707 | 13.01, 6.18, 8.91 | 3.227·10⁻¹ | 9.334·10⁻⁴ |
| centres shifted by +0.15 / +0.3 / +0.45 / +0.6 | 0.0098 / 0.0200 / 0.0311 / 0.0436 | +0.00140 / +0.00298 / +0.00480 / +0.00698 | +0.00157 / +0.00333 / +0.00534 / +0.00770 | +0.00550 / +0.00374 / +0.00173 / −0.00063 | 12.9–13.0, 6.18–6.19, 8.92 | 3.21–3.22·10⁻¹ | 9.40–9.51·10⁻⁴ |
| amplitudes ×0.5 / ×0.7 / ×1.2 | 0.0921 / 0.0624 / 0.0861 | −0.00039 / −0.00028 / +0.03684 | −0.00033 / −0.00023 / +0.03875 | +0.00740 / +0.00731 / −0.03167 | 9.76, 5.01, 7.47 / 10.59, 5.40, 7.92 / 17.80, 6.96, 10.11 | 1.991 / 2.303 / 5.242·10⁻¹ | 7.106 / 7.657 / 18.02·10⁻⁴ |
| widths ×0.7 / ×1.5 / ×2 | 0.2195 / 0.0581 / 0.1190 | +0.03086 / +0.00113 / +0.00294 | +0.03316 / +0.00126 / +0.00306 | −0.02609 / +0.00581 / +0.00401 | 13.00, 2.82, 9.77 / 11.13, 7.19, 8.50 / 10.39, 7.27, 8.71 | 3.062 / 2.932 / 3.008·10⁻¹ | 17.16 / 7.459 / 6.148·10⁻⁴ |

The matching prior makes the exact representation stationary with zero residual (fit error 1.2·10⁻¹⁶), κ_R = 0 and the maximal margin
κ = σ₀ = 0.00707, 3.6 times the σ-modification's 0.00196 (the theorem's relaxation time 1/κ: 141 against 511 slow units); the σ-modification's bias (amplitudes
2.45, −2.29, 2.05 against 3A_k = 3.0, −2.1, 2.4; widths 0.83, 0.70, 0.88 against 0.5, 0.4, 0.5) is removed. The floors move the other way:
‖S_k‖ = 13.01, 6.18, 8.91 against 8.38, 4.66, 8.32, E D_∞ 1.61 times and the function-level floor 2.18 times the σ-modification's. Misplacing
the prior brings the residual curvature back and closes the window at centres shifted by +0.6 (margin −0.00063), at amplitudes ×1.2
(−0.03167) and at widths ×0.7 (−0.02609), while shrinking the amplitudes (×0.5, ×0.7) gives margins above σ₀ (+0.00740, +0.00731: κ_R is
negative at every atom) with floors below the matching prior's (E D_∞ 0.199, 0.230; function level 7.1·10⁻⁴, 7.7·10⁻⁴) at fit errors 0.09, 0.06.
The cell margins (the gap between the two nearest prior atoms along the configuration) are 1.2–6.3, above 2ϱ′ = 0.2 everywhere.

*Admissible weights (table 2).* The splitting margin of the replicated weighted stationary configuration Θ*(w) under the three priors:

| weights w | σ-modification | matching prior | amplitudes ×0.7 |
|---|---|---|---|
| (7, 5, 4)/16 | +0.00257 | −0.00544 | +0.00733 |
| (9, 4, 3)/16 | −0.00219 | −0.00579 | +0.00570 |
| (10, 3, 3)/16 | −0.00377 | −0.00415 | +0.00536 |
| (0.5, 0.3, 0.2) | −0.00062 | −0.00653 | +0.00663 |

The matching prior's admissible set is the smallest: any excess weight overshoots the plant (κ_R at the over-weighted atom +0.0112 to +0.0136),
so even (7, 5, 4)/16 is a saddle, and the set is smallest in the direction of the largest bump (an excess ε on the first atom alone, w = w̄ +
ε(1.5e₁ − 0.5), has margins +0.00279, −0.00031, −0.00392 at ε = 0.02, 0.04, 0.08, while ε = 0.08 on the second or third atom alone keeps
+0.00208, +0.00316); the σ-modification undershoots (f̂/f = 0.77–0.78 at the first atom) and tolerates (7, 5, 4)/16 but not (9, 4, 3)/16
(Theorem D‴'s inadmissible seed); the shrunk prior keeps all four admissible, with margins +0.00536 to +0.00733 — above σ₀ at (7, 5, 4)/16 and
slightly below it at the other three, where its over-weighted first atom overshoots the plant at the centre by 1–3% (f̂/f = 1.020, 1.028,
1.007; κ_R there +0.0014, +0.0017, +0.0004; Remark D‴.11). At the realised weights of the matching-prior coupled runs
below, two of the nine (seed, m) pairs are saddles (m = 16, seed 0 with w = (6, 5, 5)/16, margin −0.00052; m = 64, seed 1 with w = (26, 17, 21)/64,
−0.00296), the other seven admissible (+0.00056 to +0.00609).

*Coupled runs with the matching prior (table 3).* m = 16, 64, 256 against M = 2048, warm start the mixture of Gaussians of width 0.02 around the
exact representation, 300 slow units (2.1 relaxation times), three seeds. Along the run κ_R[μ*] over the spheres is 0 (r* = 0) and with the
reference's own residual at most 0.00674 < σ₀ (per-seed maxima 0.00634, 0.00674, 0.00155, each at the first record; the seed-mean curve peaks
at 0.00488); the reference's fit error falls from 0.029 to 0.008 and its rms drift from 3.6·10⁻³ to 3.0·10⁻⁵ (seed means; per seed the fit
error 0.041, 0.036, 0.011 → 0.012, 0.009, 0.003).
Window maxima of D_m over (0, 37.5], (37.5, 75], (75, 150], (150, 300]: 0.114 → 0.193 → 0.239 → 0.250 ± 0.017 (m = 16), 0.0954 → 0.157 → 0.190 →
0.197 ± 0.12 (64), 0.00933 → 0.0140 → 0.0166 → 0.0173 ± 0.010 (256): last-to-previous 1.05, 1.04, 1.04 — saturated within 2.1 of the theorem's relaxation
times 1/κ; the observed late-time e-folding times of D_m (from the geometric decay of its increments over (150, 300]) are 48–84 slow units,
median 59, against 107–306 (median 195) for the σ-modification's warm runs of Theorem D‴'s numerics — a factor 3.3 in the medians, against
3.6 in 1/κ; slope −1.23 ± 0.40 in m (the function-level distance −0.63 ± 0.21). Per seed, the linearised law from the realised weights
predicts the measured D_m at s = 300 with ratios 0.81–1.11 at parameter level (nine pairs) and 0.61–1.69 at function level, and the exact
form of Proposition D‴.5(ii), Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖²_{G⁻¹}, with ratios 0.987–0.998 in all nine pairs: the linearised law's
deviations are its second-order remainder (its two lowest ratios, 0.81 and 0.82, sit at |w^m − w^M|_∞ = 0.084 and 0.033), not the saddles'
doing. The two inadmissible pairs (ratios 0.91 and 0.81; exact 0.995 and 0.994) sit on their replicated saddles at s = 300; continued
autonomously to s = 6000 and 9000 (referee round 10, `review/REVIEW_round10_D3c.md`), both m-systems leave the saddle by Proposition K's
local split — the over-weighted cluster separates into two sub-clusters along the centre direction inside its own cell, the cluster weights
unchanged — and the coupled distance stays on its plateau (−1.3% at m = 64), unlike the σ-modification's seed of round 9, which migrated
between basins.

*The function-level law for the σ-modification (table 4).* Against the warm runs of Theorem D‴'s numerics (M = 2048, last window): the predicted
root-mean-square (1/m − 1/M)^{1/2}(tr FC̃Fᵀ/N)^{1/2} is 0.0419, 0.0207, 0.00984 against the measured means 0.0307 ± 0.0097, 0.0154 ± 0.0033,
0.00467 ± 0.0012 (ratios 0.73, 0.74, 0.48; a mean of norms over four seeds lies below the root-mean-square).

*Reading.* The prior is a design input on three counts. Its placement sets the stationary residual, hence the window's lower end and the
relaxation rate (the matching prior: margin σ₀, the theorem's relaxation time 3.6 times shorter, the observed one 3.3 times); the
amplitudes it induces set the floors (which rise with them along the amplitude-scaled priors, so removing the bias raises the floors by
1.6–2.2 times); and the sign of the residual it leaves at each cluster sets the admissible set of
weights (overshoot splits, undershoot stabilises). A shrunk prior — amplitudes 0.5–0.7 of the plant's — wins on all three at a fit error
of 0.06–0.09; a misplaced one (centres off by more than a bump width, amplitudes too large, widths too narrow) closes the window.

*Remark D‴.11 (undershoot stabilises, overshoot splits).* For a Gaussian kernel, ∂²_{cc}Φ̃(x;ϑ) = −(a/s²)(1 − (x − c)²/s²)φ is negative at the
centre, so the (c, c) entry of R[μ](ϑ_k) = ⟨r_μ∇²Φ̃⟩_ρ has the sign of −a_k·r_μ near the centre: a cluster whose kernel overshoots the plant
there (r of the sign opposite to a_k) gets R ≻ 0 in the centre direction — the separation block σ₀G⁻¹ − R loses margin, and the cluster
splits once κ_R exceeds σ₀ — while an undershooting cluster gets R ≺ 0 and a splitting margin above σ₀. For the matching prior this is
exactly what table 2 shows: at every tested weight the over-weighted atom overshoots at its centre (f̂/f = 1.07–1.10), its R_cc is positive
and the top eigenvector of R^G is the centre direction, (0.01, 1.00, −0.02) in (a, c, s), while the under-weighted atoms have R_cc < 0; the
shrunk prior undershoots at the uniform weights and at (7, 5, 4)/16, and at the three more unbalanced weights its over-weighted atom
overshoots slightly (f̂/f = 1.007–1.028), taking the margin just below σ₀. The two other instabilities of §5d and §6e have a different
anatomy (referee round 10): the σ-modification's inadmissible seed (weight 9/16 on the first bump) still undershoots at its centre
(R_cc < 0), its unstable direction being a mixed amplitude–width direction — the broad kernel (width 0.88 against the bump's 0.5) overshoots
in its tails, through the cross-curvature ⟨r φ(x − c)²/s³⟩_ρ; and Proposition K's collapsed kernel owes its R_cc = +0.055 to the overshoot
over the negative bump and to the undershoot of the outer bumps beyond its inflection points, where 1 − (x − c)²/s² < 0, not to the region
between the bumps. It is an observation about the Gaussian kernel's second derivative, not a theorem: the width and amplitude directions
carry other signs, and the residual's profile matters; the splitting margin itself is the quantity to compute.

**Proposition D‴.12 (the mass constraint removes the weight-fluctuation floor).** Let ν = Σ_{j=1}^{K₀}q_jδ_{ϑ⁰_j} be a prior measure
with distinct atoms and masses q_j > 0, Δ_ν := min_{j≠k}‖ϑ⁰_j − ϑ⁰_k‖_{G⁻¹}, and let the leak be the barycentric OT leak of Theorem C′
toward ν for the cost ‖·‖²_{G⁻¹}: with γ = γ(Θ) a vertex optimal plan from μ_Θ to ν (row sums 1/m, column sums q_j) and
T̄_i := mΣ_jγ_{ij}ϑ⁰_j its barycentric projection,

    dϑ_i/ds = G⟨∇_ϑ(r_ΘΦ̃(·;ϑ_i))⟩_ρ − σ₀(ϑ_i − T̄_i),                                                                  (𝔄_m^{OT})

the G-gradient flow of J_m^{OT}(Θ) := mL̄(μ_Θ) + (σ₀m/2)W²_{G⁻¹}(μ_Θ, ν) wherever the optimal plan is unique (Danskin's theorem, as in
Theorem C′(i)), with the vertex-selection solutions of §4b where it is not. Its mean-field form is the *Laguerre-cell leak*: for μ giving
no mass to their boundaries, the optimal plan from μ to the atomic ν is induced by the map T_μ(ϑ) = ϑ⁰_j on the Laguerre cell
L_j[μ] := {ϑ : ‖ϑ − ϑ⁰_j‖²_{G⁻¹} − ψ_j ≤ ‖ϑ − ϑ⁰_k‖²_{G⁻¹} − ψ_k for all k}, the dual weights ψ fixed by μ(L_j[μ]) = q_j (semi-discrete optimal
transport [S42]; for a disconnected support the boundaries in the gaps are not unique, but T_μ is unique on supp μ), and (𝔄₀^{OT}) is (𝔄₀)
with the leak −σ₀(ϑ − T_{μ_s}(ϑ)): the nearest-atom leak of Proposition D‴.9 with the Voronoi cells replaced by cells of prescribed
masses. Call m *commensurate* with ν if mq_j ∈ ℕ for every j.

(i) *(conservation against constraint)* Along (𝔄_m^ν) on the good event 𝒢_S of Theorem D‴ (Proposition D‴.9(i)) no particle changes
cell, so the cluster masses w_k(Θ_s) := (1/m)#{i : j(ϑ_i(s)) = k} are conserved: w^m(s) = w^m(0), the multinomial sample of the initial
cloud, and the stationary floor of Proposition D‴.5 is the memory of the initial condition that this conservation law carries. Along
(𝔄_m^{OT}) the masses assigned to the atoms are constrained instead: Σ_iγ_{ij}(s) = q_j for every s.

(ii) *(the frozen plan)* If m is commensurate, every vertex of the transportation polytope is integral — γ_{ij} ∈ {0, 1/m}: each particle
is assigned to one atom π(i) and atom j receives exactly mq_j particles — so T̄_i = ϑ⁰_{π(i)} and (𝔄_m^{OT}) is (𝔄_m^ν) with the Voronoi
assignment j(·) replaced by the plan's assignment π, locally constant in Θ wherever the optimal plan is unique; there J_m^{OT} is C² with
the leak's Hessian σ₀G⁻¹ ⊗ I. Theorem D‴'s synchronous coupling applies to (𝔄_m^{OT}) only under the hypothesis of Proposition D‴.9(i)
read with the plan — for every s ≥ 0 the plan assigns each particle to the atom of its reference's Laguerre cell L_·[μ_s], and
N_{ϱ′}(supp μ_s) meets no Laguerre boundary — and at s = 0 this requires the sample to be *balanced*, w^m(0) = q (the mean-field flow of
the mixture start has masses q from the start): an event of probability of order m^{−(K₀−1)/2} under an i.i.d. start — 0.053, 0.014,
0.0034, 0.0016 at m = 15, 60, 240, 510 on the instance, ≈ 0.83/m. On that event Theorem D‴, Proposition D‴.5 and Corollary D‴.7 go through
with the same constants except for Step 8, whose McDiarmid bound for the reference's sampling term is written for i.i.d. particles and
would have to be redone for the stratified sample (mq_k particles i.i.d. within each class). For a generic sample the plan reassigns at
least ½‖w^m(0) − q‖₁ of the mass at s = 0 and the shared-particle coupling carries that mismatch for ever ((v)); what is proved about the
OT dynamics in general is the local equivalence above and the stationary statement (iii).

(iii) *(the stationary measure does not depend on m)* Let Θ*(q) = (ϑ*_1, …, ϑ*_{K₀}) be the weighted stationary point of 𝒢_q with the
prior atoms (Proposition D‴.5 with w = q, each class k leaking toward ϑ⁰_k as in Proposition D‴.9) and ε := max_k‖ϑ*_k − ϑ⁰_k‖_{G⁻¹}. If
ε < Δ_ν/2 then, for every commensurate m, the replicated configuration Θ*_m(q) (mq_k particles at ϑ*_k) has the diagonal plan (class k →
atom k) as its unique optimal plan, is a critical point of J_m^{OT}, and is a strict local minimiser iff q ∈ W_adm (its Hessian is
Proposition D‴.5's, with separation blocks σ₀I − R^G_k of multiplicities mq_k − 1). The stationary measure μ*(q) := Σ_kq_kδ_{ϑ*_k} and the
learned function f̂*(q) = Σ_kq_kΦ̃(·;ϑ*_k) are therefore the same for every commensurate m, and the stationary floors vanish:

    W_{G⁻¹}(μ*_m, μ*_M) = 0,    ‖f̂*_m − f̂*_M‖_{L²(ρ)} = 0    for all commensurate m, M,

against E[D_∞ | adm] = (1/m − 1/M)Σ_kw̄_k tr(S̃_kC̃S̃_kᵀ) and (1/m − 1/M)tr(FC̃Fᵀ)/N under the nearest-atom leak (Propositions D‴.5 and
D‴.9(iii)), whose stationary measures are moreover at squared W₂-distance ½‖w^m − w^M‖₁·O(Δ_ν²) = O(Δ_ν²m^{−1/2}) from each other (W₂
itself O(Δ_νm^{−1/4})) — the cross-cluster transport of the weight mismatch, which the shared-particle distance D does not see. The price is the bias of the masses: the
stationary residual r*(q) = f − Σ_kq_kΦ̃(·;ϑ*_k) vanishes iff ν carries an exact representation of the plant *with its masses*,
Σ_kq_kΦ̃(·;ϑ⁰_k) = f on supp ρ; for masses q_δ displaced from those of the exact representation the residual is O(δ) and the window's margin
decreases with δ (Proposition D‴.9(ii) with the masses as the parameter), and the designed masses must themselves be admissible, q ∈ W_adm:
past the admissible set the over-weighted equal-mass cluster splits (Remark D‴.11).

(iv) *(incommensurate m: the quantisation remainder)* For any m and any configuration a vertex optimal plan splits at most K₀ − 1 particles
between two or more atoms (a vertex of the transportation polytope is supported on a forest of the bipartite graph with m + K₀ nodes, so on
at most m + K₀ − 1 edges, every particle carrying at least one). At a stationary configuration whose unsplit particles are collapsed by
class (the minimisers of Proposition D‴.5's type) a split particle sits where the fit force balances the pull toward its barycentric target
T̄_i ∈ conv{ϑ⁰_j : γ_{ij} > 0} — a *bridge particle* between the clusters — so the stationary measure differs from μ*(q) by at most
(K₀ − 1)/m of mass at the bridge positions and by the O(K₀/m) shift of the class weights:

    W²_{G⁻¹}(μ*_m, μ*(q)) ≤ 4(K₀ − 1)R_Θ²/m + O(K₀²‖S‖²/m²),    ‖f̂*_m − f̂*(q)‖_{L²(ρ)} ≤ 2(K₀ − 1)K_Φ/m + O(K₀(K_Φ + L_ϑ^G‖S‖)/m),

R_Θ the a priori radius of Step 0 (with the hull of supp ν). The orders are the content, not the constants — on the instance R_Θ = 2865
makes the first bound 4·10⁶ against a measured 1.0–1.4 at m = 16, and the second 1.00 against 0.018: the function-level floor is O(1/m²)
in square, one order in m below the (1/m − 1/M) law, and the W² floor is O(1/m) against the nearest-atom leak's O(m^{−1/2}) of (iii) —
half an order gained in W², a full order at function level — with deterministic constants set by the plan's split pattern (measured
m·W²(μ*_m, ν) ≈ 20 = 0.5Δ_ν² and m·‖f̂*_m − f‖/‖f‖ ≈ 0.55).

(v) *(the transit)* From a cloud with cluster masses w^m ≠ q the plan assigns at least ½‖w^m − q‖₁ of the mass (up to the split
particles) of the over-weighted clusters to under-weighted atoms, and these particles migrate between the clusters against the fit force
at the leak's rate — for the pure leak a migrating particle crosses the Voronoi midpoint at s = ln 2/σ₀. During the transit the
synchronous coupling of Theorem D‴ sees the mass migrating in one system and not in the other at a distance of order Δ_ν,
D = O(Δ_ν²‖w^m − w^M‖₁) — it is not a coupling the OT leak respects — W²_{G⁻¹}(μ_m, μ_M) is of the same order O(Δ_ν²m^{−1/2}), the
sampling term of the cluster masses, and the function-level distance is O(K_Φ‖w^m − w^M‖₁) = O(m^{−1/2}). After the transit the two
measure-level distances vanish on the commensurate sizes ((iii)), while the synchronous D converges from below to the plateau
Σ_{k≠k′}p_{kk′}Δ²_{kk′} of the mass reassigned in one system and not in the other (0.95–0.96 of it at s = 600 on the instance): the
memory of the initial masses, erased from the stationary measure, is kept by the shared-particle coupling. A mean-field flow whose
clusters start with masses w ≠ q has the same transit — the Laguerre boundaries cut the over-weighted clusters and the mass ½‖w − q‖₁
migrates — but the mixture start of the numerics has mean-field masses q, so there the M-references transit only as finite samples and
are poor stand-ins for the mean-field flow during the transit.

*Status.* (i) proved on 𝒢_S (Theorem D‴'s per-particle bound keeps every particle within 3ϱ′/4 of its reference, whose cell has margin
ϱ′). (ii) the local equivalence proved (the transportation constraints are the incidence matrix of a bipartite graph, totally unimodular, so
with integral right-hand sides every vertex is integral — Hoffman–Kruskal [S40]; a unique optimal vertex of a linear program stays optimal
under small perturbations of its cost; the Laguerre form of the mean-field leak is semi-discrete optimal transport [S42]); the coupling
statement formal — the proofs of Theorem D‴ and Proposition D‴.5 use the leak only through its Hessian and its a priori bound, as in
Proposition D‴.9(i), the constancy of the plan being needed for every s ≥ 0 (referee round 10), which at s = 0 is the balanced-sample
event of probability ≈ 0.83/m, and on it Step 8 would have to be redone for the stratified sample (referee round 11). (iii) proved (c-cyclical monotonicity [S41]: for ε < Δ_ν/2
every off-diagonal pair costs at least ½(Δ_ν − ε)² > ½ε² ≥ every diagonal pair, so along any cycle the diagonal plan is strictly cheaper;
the diagonal plan is therefore optimal, and every other plan, which differs from it by a nonnegative combination of cycle flows, costs
strictly more; the Hessian is Proposition D‴.5(i)'s block splitting with the plan frozen; the bias statements are Proposition D‴.9(ii)'s
implicit-function argument with the masses as the parameter). (iv) the count proved (the vertices of the transportation polytope are the
plans supported on spanning forests, Klee–Witzgall [S40]); the two bounds formal (they assume the stationary configuration's structure).
(v) formal; the transit and the three distances are numerical (`exp_ot_leak.py` → `results/chaos/ot_leak.md`, below). Open: a
coupling theorem for the OT dynamics from an i.i.d. start — Theorem D‴'s synchronous coupling needs a balanced sample and otherwise
carries the transit's mismatch for ever (the plan reassigns O(√m) particles at s = 0, which the shared-particle pairing never re-matches);
the natural candidate follows the plan, pairing a particle with its assigned atom's mean-field mass — so what is proved about the OT
dynamics is the local equivalence (ii) and the stationary statement (iii). Checked: `tests/test_ot_leak.py` (the plan is a vertex with the prescribed masses,
integral when m is commensurate and with at most K₀ − 1 split rows otherwise, and the single-atom prior gives the σ-modification; the leak is
the G-gradient flow of J_m^{OT} against finite differences, for commensurate and incommensurate m; the replicated weighted configuration is
stationary for three commensurate sizes with the diagonal plan and the same learned function — for the matching masses and for mismatched
ones — with W² = 0 between the sizes and the splitting margin of Proposition D‴.5 (σ₀ at the matching masses); the conservation of the
cluster masses under the nearest-atom leak against their equalisation under the OT leak on a short run, with the mass in transit equal to
½‖w(0) − q‖₁ and the fit error lower at the end; the forest bound on random clouds; the exact W² between empirical measures by replication
against the transportation LP).

*Remark D‴.13 (conservation, constraint, and what the mean-field limit forgets).* Theorem D‴ bounds the coupled distance uniformly in time,
and Proposition D‴.5 says what the bound's constant is made of: the m-system's memory of its initial cluster masses, which the nearest-atom
leak conserves. The OT leak replaces the conservation law by a constraint and the memory is erased from the stationary measure — the finite-m
stationary measure *is* the mean-field one on the commensurate sizes, and within (K₀ − 1)/m of mass of it otherwise — though not from
the shared-particle coupling, whose synchronous distance keeps the plateau of the reassigned particles ((v)). This is a statement about the design of the
leak, not about m: the same clouds under the nearest-atom leak keep their floors. Three consequences. (a) The masses join the atoms as
design inputs of the leak's prior — both are read from the plant's structure (Proposition D‴.9(ii)) — and the admissible set W_adm of
Proposition D‴.5, a probability statement about the realised weights under the nearest-atom leak, becomes a condition on the design,
q ∈ W_adm, under the OT leak. (b) The mean-field form of Theorem C′'s leak for an atomic prior is the Laguerre-cell leak; its finite-m
fluctuation is not a random displacement of the cells but the split particles — at most K₀ − 1 — and the migrating ones during the transit.
(c) A soft mass constraint (an entropic or penalised plan) would interpolate between the two leaks, with a floor set by the stiffness of
the penalty; open (§8, branch 8), as is the transfer of a learned measure to a plant whose representation weights are unknown, where q
itself would have to be learned.

*Numerics* (`exp_ot_leak.py` → `results/chaos/ot_leak.md`). The instance of Theorem D‴'s numerics at σ₀ = 0.15σ_crit = 0.00707 with the
matching prior of Proposition D‴.9's numerics — three atoms at the exact representation, masses q = (1/3, 1/3, 1/3), the plant's
representation weights, Δ_ν = 6.26 — under the nearest-atom leak and under the OT leak (plan recomputed every 0.2 slow units) on the same
warm clouds (the jittered mixture of width 0.02), 600 slow units = 4.2/σ₀, four seeds; commensurate sizes m = 15, 60, 240 against M = 510
and incommensurate sizes m = 16, 64, 256 against M = 512.

*Conservation against constraint (table 1).* Under the nearest-atom leak no cluster weight changes along any run (max_s|w(s) − w(0)|_∞ = 0
in all 32 systems). Under the OT leak the mass in transit at s = 0 is ½‖w(0) − q‖₁ = 0.117, 0.058, 0.027, 0.024 (seed means for m = 15, 60, 240,
510; per system 1–3, 1–5, 2–13 and 3–21 particles), it reaches zero, and the weights by nearest atom equal q exactly from s_eq = 93–164
slow units on (median 104, against the pure leak's midpoint crossing ln 2/σ₀ = 98); no plan is split. The fit error at the end falls
from 0.067, 0.029, 0.017, 0.012 under the nearest-atom leak — the bias of the realised weights — to 2.8·10⁻⁴, 1.2·10⁻⁴, 1.2·10⁻⁴,
0.8·10⁻⁴ (seed means); the migrated particles are still 0.06–0.16 whitened units from their new atoms at s = 600 — the pure leak's
remaining approach Δ_{kk′}e^{−σ₀·600} = 0.09–0.10 — and the fit force cancels their first-order effect on f̂, so the OT finals are not yet
stationary (maximal whitened drift 5–8·10⁻⁴, against ≤ 10⁻³ and mostly ≤ 10⁻⁴ under the nearest-atom leak).

*The three distances on the commensurate sizes (table 2, last window (450, 600]).* Function level ‖f̂_m − f̂_M‖²_{L²(ρ)}: nearest-atom
2.5·10⁻³, 5.2·10⁻⁴, 5.9·10⁻⁵ (stationary, last/previous window 0.97–1.00; the law (1/m − 1/M)tr(FC̃Fᵀ)/N predicts 4.0·10⁻³, 9.1·10⁻⁴,
1.4·10⁻⁴) against OT 4.9·10⁻⁷, 4.9·10⁻⁸, 7.5·10⁻⁹ — ratios 2·10⁻⁴, 9·10⁻⁵, 1.3·10⁻⁴ — still falling by a factor 0.02–0.035 per 150 slow
units (a rate of 3.0–3.8σ₀: the displaced particles' first-order effect cancels, the class means agreeing to first order) toward the zero
of (iii). W²_{G⁻¹}(μ_m, μ_M): nearest-atom 5.1, 2.7, 0.81 — the cross-cluster transport of the weight mismatch: the exact K × K transport
between the weighted stationary measures at the realised weights reproduces it to 0.997–0.999 at the seven admissible pairs (0.85–0.99 at
the five pairs with saddle weights, where the m-system is not at the replicated configuration), the formula ½‖w^m − w^M‖₁Δ_ν² = 4.9, 2.5,
0.70 to −7% to +28%; slope −0.67 ± 0.13 in m, the sampling term of the cluster masses — against OT 7.9·10⁻³, 2.3·10⁻³, 6.9·10⁻⁴ (ratios
1.5·10⁻³, 0.9·10⁻³, 0.9·10⁻³), falling by 0.12 per 150 slow units (e^{−2σ₀·150} = 0.12, the square of the spread's contraction at the
margin σ₀). The synchronous D is the one distance the OT leak does not improve: 5.8, 3.3, 0.94 against 0.87, 0.25, 0.017, still rising
(last/previous 1.08–1.10) toward the plateau Σ_{mismatched}Δ²_{kk′}/m of the pairs reassigned in one system and not in the other
(0.95–0.96 of it at s = 600), as (v) says: the coupling keeps the memory that the measure loses.

*The incommensurate sizes (table 3).* Every plan has one or two split particles (≤ K₀ − 1 = 2); the weights by nearest atom differ from q
by the split particle's share, 0.042, 0.010, 0.0026, 0.0013 = (2/3)/m; the unsplit particles in transit at s = 0 (0–2 at m = 16, 1–4 at 64, 2–11 at 256, 1–19 at 512) all
arrive by the end. At the end the split particles sit between the clusters, at 0.28–0.58 Δ_ν from the nearest atom and within 0.03–0.6
(G⁻¹ metric) of their barycentric target T̄_i — the bridge particles of (iv). Against the M = 512 reference, which carries its own two
bridge particles: function level OT 3.4·10⁻⁴, 3.1·10⁻⁵, 1.1·10⁻⁶ (the m = 16 value still relaxing, 0.82 per window — continued to
s = 2400 by referee round 11 it settles at 2.8·10⁻⁴, and m = 64 at 2.9·10⁻⁵, with slope −2.00) against nearest-atom 2.5·10⁻³, 5.0·10⁻⁴,
4.4·10⁻⁵, ratios 0.11–0.14, 0.063, 0.025; W² OT 1.16, 0.29, 0.041 (stationary) against nearest-atom 4.8, 2.7, 0.69, ratios 0.24, 0.10,
0.059. The quantisation laws of (iv) are read against the mean-field stationary measure μ*(q), here ν itself: W²(μ_m, ν) = 1.21, 0.33,
0.079, 0.038 for m = 16, 64, 256, 512, i.e. m·W² = 19.4, 21.0, 20.1, 19.5 ≈ 0.5Δ_ν² with slope −0.98 ± 0.03, and the fit error 0.028,
0.0093, 0.0023, 0.0011, i.e. m·fit = 0.45–0.60 with squared slope −1.81 ± 0.05 — the orders 1/m and 1/m² of (iv), with constants set by
the plan's split pattern (seeds with the same pattern reach identical plateaus, referee round 11) far below the worst case (K₀ − 1)Δ_ν²
= 78; against the nearest-atom leak's fit errors 0.061, 0.035, 0.013, 0.012.

*The masses as a design input (table 4).* For q_δ = (1/3 + δ, 1/3 − δ/2, 1/3 − δ/2) the weighted stationary point has fit error 0.0069,
0.013, 0.019, 0.029, 0.049 at δ = 0.01, 0.02, 0.03, 0.05, 0.1 — the bias of the mass constraint, 0.69δ to 0.49δ, decreasing with δ — κ_R at the over-weighted atom
+0.0023, +0.0043, +0.0060, +0.0086, +0.0119, margin κ +0.0047, +0.0026, +0.0008, −0.0019, −0.0054 and splitting margin +0.0048, +0.0028,
+0.0011, −0.0015, −0.0049: the designed masses leave the admissible set between δ = 0.03 and 0.05 — the asymmetry of the matching prior's
admissible set in Proposition D‴.9's numerics (an excess of 0.04 on the largest bump is inadmissible) — while the diagonal plan stays the
unique optimal plan at every row (ε ≤ 1.04 < Δ_ν/2 = 3.13).

*Reading.* The two leaks start from the same clouds and end in different worlds. The nearest-atom leak keeps the clusters' sampled masses
for ever, and its stationary distances are what Propositions D‴.5 and D‴.9 say: a function-level floor of order 1/m, and a W² floor of
order m^{−1/2} set by the cross-cluster transport of the weight mismatch. The OT leak moves one to 21 particles across in about
ln 2/σ₀ slow units and then settles on the mean-field stationary measure: on the commensurate sizes every distance but the synchronous one
decays to zero (three to four orders of magnitude below the nearest-atom leak's floors after 4.2/σ₀), and the fit loses the weight bias;
on the incommensurate sizes the floors are the one or two bridge particles — O(1/m²) at function level, O(1/m) in W² — whose position
between the clusters is the price of a mass constraint the size cannot meet. The masses are a design input with a bias and an
admissibility constraint of their own.

**Corollary D‴.14 (the real loop with the OT leak: one orbit for every commensurate size).** Consider the real closed loop (CL) of
Theorem J with the leak σ(ϑ_i − T_i), σ = σ₀/k_e, where T_i is the nearest-atom target ϑ⁰_{j(ϑ_i)} (Proposition D‴.9) or the barycentric
target T̄_i of the vertex plan (Proposition D‴.12). For a weight vector w with mw ∈ ℕ^{K₀} and a class map π₀ : {1, …, m} → {1, …, K₀} with
|π₀⁻¹(k)| = mw_k, let 𝒞_m(π₀) := {Θ : ϑ_i = ϑ_{π₀(i)}} be the *collapsed subspace*, on which f̂_Θ = Σ_kw_kΦ̃(·;ϑ_k), and let the *weighted
K₀-atom loop* (CL_w) be (CL) with f̂ = Σ_kw_kΦ̃(·;ϑ_k), the K₀ parameters ϑ_k and the leaks σ(ϑ_k − ϑ⁰_k).

(i) *(invariance)* 𝒞_m(π₀) is invariant under (CL) with the nearest-atom leak as long as every class atom ϑ_k stays in the Voronoi cell
of ϑ⁰_k, and under (CL) with the OT leak as long as w = q and max_k‖ϑ_k − ϑ⁰_k‖_{G⁻¹} < Δ_ν/2; on it (CL) is (CL_w). With the OT leak and
w ≠ q the subspace is not invariant: the plan assigns at least ½‖w − q‖₁ of the mass of the over-weighted classes to other atoms, and those
particles leave their class.

(ii) *(one orbit for every commensurate size)* Let q ∈ W_adm and let Γ_q be the T₀-periodic orbit of (CL_q) near Θ*(q) given by Theorem J,
whose proof applies to (CL_w) verbatim (its averaged drift is −(1/k_e)(diag(w)⁻¹ ⊗ G)∇𝒢_w, the gradient of 𝒢_w in the metric diag(w) ⊗ G⁻¹,
which is all its Steps 4–6 use), and assume Γ_q stays inside the cells, max_k‖ϑ^Γ_k(t) − ϑ⁰_k‖_{G⁻¹} < Δ_ν/2 — the cell condition of (i),
which W_adm does not contain. For every commensurate m the lift of Γ_q to 𝒞_m(π₀) is a T₀-periodic orbit of (CL) with the OT leak. It is
exponentially stable iff Γ_q is and the ghost monodromies of the classes contract (Proposition D‴.15 below) — at the matching masses, where
the residual vanishes, the ghosts' exponent is exactly σ₀/k_e > 0 — and it is Theorem J's orbit Γ_m at Θ*_m(q), with the basin (c) of that
theorem, whenever it lies in Theorem J's uniqueness ball ‖Θ − Θ*_m(q)‖_{G_m⁻¹} ≤ ρ₀(m): with ‖Θ‖²_w := Σ_kq_k‖ϑ_k‖²_{G⁻¹} the lift's
G_m⁻¹-distance to Θ*_m(q) is √m times the weighted distance, so the condition is √m·sup_t‖Θ^q_p(t) − Θ*(q)‖_w ≤ ρ₀(m), satisfied for
k_e ≥ max{k₀(m), √m C_q/ρ₀(m)} by Theorem J(a)'s sup_t‖Θ^q_p(t) − Θ*(q)‖_w ≤ C_q/k_e. Along the lift the orbit measure
(1/m)Σ_iδ_{ϑ_i(t)} = Σ_kq_kδ_{ϑ_k(t)}, the learned function Σ_kq_kΦ̃(·;ϑ_k(t)), the plant's orbit x_p(t) and the tracking error e(t) are the
same for every commensurate m at every t — exactly, not up to O(1/k_e): the finite-size floor of the real loop (Corollary D‴.7) is zero.

(iii) *(Corollary D‴.7 sharpened)* With the nearest-atom or the σ-modification leak and both realised weight vectors admissible,
w^m, w^M ∈ W_adm (Corollary D‴.7's hypothesis: the replicated configurations are hyperbolic minima), the same argument lifts Γ_{w^m}, the
orbit of (CL_{w^m}) with the realised weights: Theorem J's orbit near Θ*_m(w^m) is the lift of Γ_{w^m} under the same proviso, so the real
loop's coupled distance between two sizes is exactly the distance between two weighted K₀-atom orbits,

    D_orb(t) = Σ_kw^m_k‖ϑ^Γ_k(t; w^m) − ϑ^Γ_k(t; w^M)‖²_{G⁻¹},

Corollary D‴.7's O(1/k_e) terms being the k_e-dependence of w ↦ Γ_w and nothing else; and the tracking error along the orbit is that of
(CL_{w^m}) — the weight bias of the fit is a tracking error, which the mass constraint removes for every size at once. When w^m ∉ W_adm the
lift is a saddle orbit (a ghost exponent is negative, Proposition D‴.15) and the m-system leaves it by splitting its over-weighted class, as
three of the six nearest-atom systems of the numerics below do.

*Status.* (i) proved (the particle law depends on i only through ϑ_i and its target; on 𝒞_m(π₀) the targets agree within a class — by the
cell condition for the nearest-atom leak, by the lemma of Proposition D‴.12(iii) for the OT leak — so the vector field is tangent to
𝒞_m(π₀) and coincides there with (CL_w)'s, and solutions are unique away from the cell boundaries and the plan switches). (ii), (iii) proved
from (i), Proposition D‴.15 and Theorem J under the stated proviso and hypotheses (the cell condition along Γ_q; admissible weights in (iii)); the m-uniformity of Theorem J's thresholds k₀(m), ρ₀(m) at the
replicated configuration is the open point of Corollary D‴.7 — its spectral part is settled by Proposition D‴.15 (the Floquet spectrum of
the lift does not depend on m), the basin radius in the mean norm remains formal; without the proviso the lift is still a periodic orbit of
(CL) for every commensurate m, with the stability Proposition D‴.15 gives it, and only its identification with Theorem J's orbit is
conditional. Checked: `tests/test_ot_leak.py` (the collapsed 15-system and the 3-atom loop integrate to the same configuration and plant
state to 10⁻¹² over 300 steps of the real loop with the OT leak at masses q, the atoms moving; the collapsed 15-system with weights
(6, 5, 4)/15 stays collapsed under the nearest-atom leak, while under the OT leak one of its particles is assigned away from its class; two
commensurate lifts, m = 15 and 30, keep W² = 0 and the same learned function along the real loop; the lifts of two weight vectors under the
nearest-atom leak are different measures).

*Numerics* (`exp_ot_real_loop.py` → `results/chaos/ot_real_loop.md`). The real loop at k_e = 25 (σ = σ₀/k_e; RK4 plant with the control and
the reference evaluated at the stages and the parameters frozen over the step, Euler parameters, dt = 10⁻³ — an earlier zero-order hold of
the control had set a tracking floor of 4·10⁻⁵, referee round 12) for 100 time units — 16 periods, 4 slow units — from the
clouds at the end of Proposition D‴.12's averaged runs (contracted where the realised weights are admissible) (seeds 0 and 1; m = 15, 60 and M = 510; the OT runs' clouds carry the masses q, the
nearest-atom runs' their conserved sampled weights) and from the collapsed lifts (the exact representation for the OT leak — the 3-atom
loop — and the replicated Θ*(w^m) for the nearest-atom leak); statistics over the last five periods on the reference orbit. Under the OT
leak the learned functions of the three sizes agree to ‖f̂_m − f̂_M‖²_{L²(orbit)} = 7.5·10⁻⁹–2.7·10⁻⁸ and with the 3-atom lift to
5.2·10⁻¹⁰–3.6·10⁻⁸ — the clouds' within-class spread (0.013–0.033, the migrated particles still approaching their atoms) at second order —
with fit errors along the orbit 3.5·10⁻⁵–2.9·10⁻⁴ and rms tracking errors 0.9·10⁻⁶–7.3·10⁻⁶ (the fit error over k_e; the collapsed lift
tracks to 8·10⁻¹², the scheme's floor), and W²(μ_m, μ_M) = 3.3·10⁻⁴–1.4·10⁻³ at the end. Under the nearest-atom leak they disagree by
1.1·10⁻⁴–2.5·10⁻³, with fit errors 0.012–0.075 and rms tracking errors 3.2·10⁻⁴–1.9·10⁻³ — the weight bias as a tracking error, 130 to 550
times the OT leak's per system — and W²(μ_m, μ_M) = 1.6–5.5, the cross-cluster transport of the weight mismatch. Three of the six nearest-atom systems
have inadmissible realised weights (seed 0, m = 15: (6, 5, 4)/15, splitting margin −0.0035; seed 1, m = 15: (6, 3, 6)/15, −0.0011; seed 1,
m = 60: (24, 15, 21)/60, −0.0022): their lifts are saddle orbits and the m-systems are not on them — the over-weighted class has split
(within-class spread 0.04–0.11 against 0.0006–0.0013 for the admissible systems), the distance to the lift, 4·10⁻⁸–3·10⁻⁶, being the split's
second-order effect on f̂, and the lift-to-lift ratio at those pairs 0.960, 0.998, 0.961. At the one pair whose weights are both admissible
(seed 0, m = 60 against M = 510) the disagreement is the distance between the lifted orbits to a ratio 1.0003, and the admissible m-systems
sit within 2.5·10⁻¹¹ of their lifts: (iii) as stated, with its admissibility hypothesis visibly necessary.

*Numerics in two dimensions* (`exp_ot_2d.py` → `results/chaos/ot_2d.md`). The mechanism is dimension-free, and the real loop shows it on a
planar plant: ẋ = f(x) + u on ℝ² with f the sum of four Gaussian bumps with vector amplitudes, the Lissajous reference x_m(t) = (2 sin t,
1.5 sin 2t), a saturated network with p = 5 parameters per kernel, G = diag(1, 1, ¼, ¼, ¼), the matching prior (four atoms at the exact
representation, masses ¼), k_e = 25 and σ₀ = 0.05 (the leak is a design input; ln 2/σ₀ = 14 slow units), 2500 time units = 100 slow units =
5/σ₀ from the same jittered warm clouds with random class labels, sizes m = 8, 32 against M = 128, two seeds, the plant integrated as above.
Under the nearest-atom leak the sampled cluster weights — (3, 2, 2, 1)/8, (7, 4, 12, 9)/32, … — never change, the fit errors along the orbit
are 0.11–0.39, the sizes disagree by ‖f̂_m − f̂_M‖²_{L²(orbit)} = 5.5·10⁻³–2.3·10⁻² and the rms tracking error is 1.9·10⁻³–6.7·10⁻³. Under
the OT leak the weights reach (¼, ¼, ¼, ¼) at 13–14 slow units in all six systems (ln 2/σ₀ = 14), the fit errors fall to 2.1·10⁻⁴–8.7·10⁻⁴
(the within-class spread 0.012–0.018 still contracting), the sizes agree to 5.0·10⁻⁸–1.6·10⁻⁷ at the stopping time 5/σ₀ —
five orders of magnitude below the nearest-atom leak, and still falling by about 3.5 times per 10 slow units (1.8·10⁻¹² at 10/σ₀ in referee
round 13's continued run) — and the rms tracking error is 3.9·10⁻⁶–1.6·10⁻⁵, 290 to 480 times smaller per system: the weight bias is a
tracking error in any dimension, and the mass constraint removes it for every size at once.

**Proposition D‴.15 (the transverse spectrum of a lifted orbit).** Let Γ = (x_p, Θ_p) be a T₀-periodic orbit of (CL) — with any leak
whose targets are locally constant near Γ and equal within each class — lying in the collapsed subspace 𝒞_m(π₀) with classes of sizes n_k, and let
Γ_K = (x_p, (ϑ_{p,k})_k) be the K₀-atom orbit it lifts. For each class k let

    S_k := {δ ∈ ℝ^{n+mp} : δx = 0, δ_i = 0 for i ∉ π₀⁻¹(k), Σ_{i∈π₀⁻¹(k)}δ_i = 0}    (dimension (n_k − 1)p),

and let M_k be the monodromy over one period of the *ghost equation* δ̇ = A_k(t)δ, A_k(t) := ∂_ϑ[G∇_ϑ(e_p(t)ᵀΦ̃(x_p(t);ϑ))]|_{ϑ=ϑ_{p,k}(t)} − σI
— the linearisation of one particle of class k driven by the orbit's plant state, with weight zero in f̂. Then the monodromy matrix 𝕄 of Γ
leaves every S_k invariant, with 𝕄|_{S_k} = M_k ⊗ I_{n_k−1}, and on the complementary invariant subspace — the plant coordinates and the
class means — it is the monodromy 𝕄_K of Γ_K as an orbit of the weighted K₀-atom loop (CL_w). Hence the Floquet spectrum of Γ is that of Γ_K
together with the spectrum of each M_k with multiplicity n_k − 1, independent of m; Γ is exponentially stable iff Γ_K is and every M_k is
contracting; and when Γ_K is Theorem J's orbit of (CL_w) at Θ*(w), the exponents of M_k are (1/k_e)λ_j(σ₀I − R^G_k) + O(k_e⁻²), j = 1..p, with
R^G_k = G^{1/2}R[μ_{Θ*(w)}](ϑ*_k)G^{1/2} — Proposition K's separation condition and the splitting blocks of Proposition D‴.5(i), read along
the orbit.

*Proof.* Write the particle law as ϑ̇_i = F(x, e, ϑ_i) − σ(ϑ_i − T_i), F(x, e, ϑ) := G∇_ϑ(eᵀΦ̃(x;ϑ)), and note that the plant's equation
depends on Θ only through f̂_Θ(x) = (1/m)Σ_iΦ̃(x;ϑ_i). Along Γ the linearised equations are δẋ = (∂_xf − k_eI − ∂_xf̂)δx − (1/m)Σ_i∂_ϑΦ̃(x_p;ϑ_i)δ_i
and δ̇_i = (∂_xF + ∂_eF)δx + (∂_ϑF − σI)|_{ϑ_i}δ_i (e = x − x_m gives δe = δx). For δ ∈ S_k, δx = 0 and all ϑ_i with i ∈ π₀⁻¹(k) equal
ϑ_{p,k}(t), so Σ_i∂_ϑΦ̃(x_p;ϑ_i)δ_i = ∂_ϑΦ̃(x_p;ϑ_{p,k})Σ_{i∈π₀⁻¹(k)}δ_i = 0: the plant is not perturbed to first order, each δ_i obeys
δ̇_i = A_k(t)δ_i with the same A_k, and Σ_{i∈π₀⁻¹(k)}δ̇_i = A_kΣδ_i = 0. So S_k is invariant and the linearised flow acts on it as the ghost
equation on every component, whence 𝕄|_{S_k} = M_k ⊗ I_{n_k−1}. The class-constant perturbations {δ_i = u_{π₀(i)}} together with δx form the
complementary subspace; it is invariant because particles of one class have identical right-hand sides when their perturbations coincide,
and there the variation of f̂ is Σ_k(n_k/m)∂_ϑΦ̃(x_p;ϑ_{p,k})u_k = Σ_kw_k∂_ϑΦ̃u_k — the linearisation of (CL_w) along Γ_K. The spectral
statements follow. For the last one, along Theorem J's orbit of (CL_w) the tracking error is e_p(t) = r_{Θ*(w)}(x_m(t))/k_e + O(k_e⁻²)
and ϑ_{p,k}(t) = ϑ*_k(w) + O(1/k_e) (Theorem J(a) with the error equation), so the ghost matrix is A_k(t) = (1/k_e)[G∂²_ϑ(r_{Θ*(w)}(x_m(t))
Φ̃(x_m(t);ϑ*_k)) − σ₀I] + O(k_e⁻²): a T₀-periodic matrix of size O(1/k_e) whose period mean is (1/k_e)(GR_k − σ₀I) + O(k_e⁻²), with
R_k = R[μ_{Θ*(w)}](ϑ*_k) the residual-weighted curvature of Theorem D‴. First-order averaging of this linear periodic equation gives
M_k = exp(T₀(1/k_e)(GR_k − σ₀I)) + O(k_e⁻²), and GR_k − σ₀I is similar to G^{1/2}R_kG^{1/2} − σ₀I = −(σ₀I − R^G_k): the exponents are
(1/k_e)λ_j(σ₀I − R^G_k) + O(k_e⁻²), with no condition on m. ∎

*Status.* Proved (linear algebra along the orbit; first-order averaging of the ghost equation for the last claim, which needs only that
Γ_K is Theorem J's orbit of (CL_w) — referee round 12 confirmed the finite-difference ghost monodromy against the matrix ODE to 3·10⁻⁹ and
against exp(T₀Ā_k) to 7·10⁻⁹, the prediction's remaining error being R at Θ*(w) against R along Θ_p(t), falling 6–12 times from k_e = 25
to 50). The m-uniformity of the Floquet spectrum settles the spectral part
of Corollary D‴.7's open point: the stability type and the rates of Theorem J(b)–(c) at a replicated configuration do not depend on m.
Checked: `tests/test_ot_leak.py` (the lift of the 3-atom orbit to m = 6 at the matching masses is a fixed point of the Poincaré map to
10⁻¹⁰; each class's separation subspace is invariant under its monodromy to 10⁻⁶ and carries the ghost's exponents, σ₀/k_e = 2.83·10⁻⁴, to
10⁻⁶; the complement carries the 3-atom loop's exponents to 10⁻⁶).

*Numerics* (`exp_lifted_orbit.py` → `results/chaos/lifted_orbit.md`). The real loop at k_e = 25, σ = σ₀/k_e; monodromy matrices by
central differences of the Poincaré map (DOP853, rtol 10⁻¹¹). (A) The OT leak at the matching masses: the 3-atom orbit is the exact
representation with zero tracking error, and its lifts to m = 6, 9, 15 are fixed points of the Poincaré map to 10⁻¹⁵; every class's
separation subspace is invariant to 2·10⁻⁸ and carries the exponent σ₀/k_e = 2.830·10⁻⁴ with multiplicity n_k − 1 — the ghosts' exponent to
4·10⁻⁵ relative — and the collective exponents 2.935·10⁻⁴, 3.558·10⁻⁴, 3.999·10⁻⁴, 8.73·10⁻⁴, 1.109·10⁻³, 1.573·10⁻³, 3.704·10⁻³, 5.523·10⁻³,
7.934·10⁻³ are the same for m = 3, 6, 9 and 15 to 10⁻⁴ relative, within 0.6–3.2% of Theorem J(b)'s λ_j(H_q)/k_e = 2.952·10⁻⁴, …,
8.091·10⁻³, H_q the whitened, weighted collective block of Proposition D‴.5(i). (B) The σ-modification at w = (7, 5, 4)/16 (admissible, splitting margin +0.00257) and (9, 4, 3)/16 (a saddle, −0.00219): the
weighted 3-atom orbits by shooting (residuals 10⁻¹⁴, period-start offsets from Θ*(w) of 1.1·10⁻² and 1.6·10⁻² — Theorem J(a)'s O(1/k_e)), and
their lifts to m = 16 fixed points to 5·10⁻¹⁴. The separation exponents of each class equal the ghost's to 2·10⁻⁵ relative (invariance defects
≤ 6·10⁻⁸) and the prediction (σ₀ − λ_j(R^G_k))/k_e within 0.2–4.1% (class 0 at (7, 5, 4)/16: 9.989·10⁻⁵, 4.553·10⁻⁴, 1.267·10⁻³ against
1.029·10⁻⁴, 4.572·10⁻⁴, 1.275·10⁻³); at (9, 4, 3)/16 the over-weighted class carries the negative exponent −9.13·10⁻⁵ (predicted −8.77·10⁻⁵)
with multiplicity 8 — the saddle orbit of Theorem D‴'s inadmissible seed, Proposition K's split at the orbit level — while the collective
exponents (2.131·10⁻⁴, …, 9.985·10⁻³) are the weighted 3-atom orbit's to 2·10⁻⁵ relative.

*Remark D‴.16 (dither and the mass constraint).* With dither (Theorem J′: τ = τ₀/k_e, σ = σ₀/k_e) the leading-order stationary law of
the slow variables is the Gibbs measure π ∝ exp(−J_m/τ₀), and the cluster masses are no longer conserved under the nearest-atom leak: a
particle crosses a Voronoi ridge at the Arrhenius rate exp(−ΔV/τ₀), ΔV the single-particle barrier out of the configuration — at least the
leak's part σ₀Δ_ν²/8 = (σ₀/2)(Δ_ν/2)² at the lowest ridge point, plus the fit's variation — and the counts sample the Gibbs law over the count
vectors n: in the Laplace approximation multinomial(n)·exp(−m𝒢*_{n/m}/τ₀) times the determinants of the separation and collective blocks of
Proposition D‴.5(i), for n/m ∈ W_adm (at a saddle weight the split configuration lies lower, and its soft and unstable separation modes carry
entropy that the quadratic approximation misses). Dither thus does on the Arrhenius time scale, and only in distribution, what the OT leak
does deterministically in ln 2/σ₀ slow units at τ₀ = 0. Under the OT leak the Gibbs measure is π^{OT} ∝ exp(−J_m^{OT}/τ₀) (J_m^{OT} is
locally Lipschitz, piecewise smooth, with Gaussian tails), the plan's masses are q by construction, and a particle that crosses a Voronoi
ridge stays assigned to its own atom: with integral plans an exchange of the assignments of particles a ∈ class j and b ∈ class k lowers the
cost only once ⟨ϑ_a − ϑ_b, ϑ⁰_k − ϑ⁰_j⟩_{G⁻¹} > 0, that is, once they have passed each other. So the counts by nearest atom leave q at the
nearest-atom leak's ridge barrier, but only for a reversible excursion — beyond the ridge the leak toward the particle's own atom keeps
pulling it back — of stationary probability of order exp(−ΔV_ridge/τ₀); exchanges of assignments (whose barrier is at most the lowest tie
point of two particles) leave the masses at q as well. The mass constraint does not raise the barrier of a crossing; it makes the crossing
reversible. Formal (the Gibbs law of a locally Lipschitz potential with Gaussian tails; the Arrhenius estimate; the Laplace asymptotics).
*Numerics* (`exp_dither_ot.py` → `results/chaos/dither_ot.md`): the dithered averaged slow loop — Langevin in slow time (Euler–Maruyama,
ds = 0.05), with the Gibbs law exp(−J_m/τ₀) — for m = 15 particles started at the counts (6, 5, 4) under the matching prior at σ₀ = 0.15σ_crit,
at τ₀ = 0.005, 0.01, 0.02 over 50 000, 200 000 and 50 000 slow units, two seeds. On the instance the single-particle barriers out of the
balanced configuration — the lowest points of the Voronoi ridges with the fourteen other particles relaxed — are 0.040–0.074 (the leak's
part 0.035), and the lowest exchange barrier between two classes is 0.083 (referee round 13). Under the nearest-atom leak the record-level
label changes (every 10 slow units; about half of them undone at the next change, recrossings) occur at 1.0·10⁻⁵, 5.1·10⁻⁴, 5.6·10⁻³ per slow
unit — an Arrhenius barrier of 0.048 (0.043–0.052 at 95%) from τ₀ = 0.01 and 0.02, among the ridge minima, while τ₀ = 0.005 has one event in
100 000 slow units — so at τ₀ = 0.005 the counts (6, 5, 4) stand for 350/σ₀, and at τ₀ = 0.01 and 0.02 the counts sample a distribution whose
weight on (5, 5, 5) is 0.54 and 0.25, against 0.58 and 0.33 for the energy–multinomial law multinomial(n)·exp(−m𝒢*_n/τ₀) — a comparison that
cannot discriminate at τ₀ = 0.01 (an integrated autocorrelation time of 7 400–23 000 slow units leaves 9–27 effective samples per seed, a
standard error of about 0.1) and puts the law 2.1 standard errors above the data at τ₀ = 0.02; the law under-weights the saddle count
vectors (6, 4, 5) and (6, 5, 4) by 3.2 and 2.7 times (total-variation distance 0.20, against 0.07 between the two seeds): their soft and
unstable separation modes carry the excess, the very factors the law omits. Under the OT leak the counts by nearest atom are balanced in
every record after the initial transit at τ₀ = 0.005 (no excursion in 50 000 slow units; 99.6% of the whole run), and in 99.4% and 88% of the
records after it at τ₀ = 0.01 and 0.02; the excursions are returns of the crossing particle — all 121 at τ₀ = 0.01, 298 of 336 at 0.02, the other 38
exchanges — the label changes are about as frequent as the nearest-atom leak's (ratios 1.18–1.30), and the off-balance fraction
(6.3·10⁻³, 0.117) has an effective barrier of 0.058, inside the ridge range.

**Proposition D‴.17 (convergence of the OT dynamics).** Let J_m^{OT} be the energy of Proposition D‴.12 for the saturated Gaussian
features, with ρ the finite quadrature of the code or the period average of an analytic reference, on the a priori ball of Step 0.

(i) *(a Lyapunov function with a semiconcave part)* J_m^{OT} = mL̄(μ_Θ) + (σ₀m/2)min_γΣ_{ij}γ_{ij}‖ϑ_i − ϑ⁰_j‖²_{G⁻¹} is the minimum over the
finitely many vertices γ of the transportation polytope of the smooth energies J^γ_m, so it is locally Lipschitz, the sum of a smooth
function and a semiconcave one (a quadratic in Θ minus a maximum of linear functions), and its Clarke subdifferential at Θ is the convex
hull of the gradients ∇J^γ_m(Θ) of the active vertices. The vertex-selection solutions of (𝔄_m^{OT}) solve the differential inclusion
dΘ/ds ∈ −G_m∂J_m^{OT}(Θ); along them J_m^{OT} is non-increasing, with dJ_m^{OT}/ds = −‖∇J^γ_m(Θ)‖²_{G_m} wherever the plan γ is unique, and
at a tie the directional derivative along either active plan's drift is at most −‖∇J^γ_m‖²_{G_m}.

(ii) *(convergence)* J_m^{OT} is definable in the o-minimal structure ℝ_{an,exp} on the a priori ball — finitely many compositions of
polynomials, exponentials and the saturations tanh and sig, a finite quadrature sum (or an analytic function of Θ for the period average
of an analytic reference), and a minimum over finitely many such functions; for the finite quadrature it is definable in ℝ_exp on all of
ℝ^{mp}, and for the period average of an analytic reference in ℝ_an on the ball — so it satisfies the Kurdyka–Łojasiewicz inequality in the
nonsmooth sense [S43]. The change of variables Y = G_m^{−1/2}Θ maps dΘ/ds ∈ −G_m∂J into dY/ds ∈ −∂(J∘G_m^{1/2})(Y), the form in which the
theorems of [S43] are stated; definable locally Lipschitz functions obey the chain rule along absolutely continuous curves, so along a
subgradient trajectory dJ/ds = −‖dΘ/ds‖²_{G_m⁻¹} [S43], and with the KL inequality every bounded subgradient trajectory has finite length
and converges to a Clarke-critical point. Every vertex-selection and every generalised solution of (𝔄_m^{OT}) is such a trajectory, bounded
by Step 0's a priori estimate with conv supp ν; hence it converges, and so do the solutions of (𝔄_m^ν) (a minimum over the Voronoi
assignments) and of (𝔄_m) (smooth: Łojasiewicz's gradient inequality for real-analytic functions, in Kurdyka's form for definable ones
[S44]).

(iii) *(the limit)* At a Clarke-critical point whose optimal plan is unique, the plan is a vertex and the point is a critical point of the
frozen energy J^γ_m; for commensurate m the plan is integral (an assignment π) and the limit is a critical point of the frozen energy
J^π_m, every particle leaking toward its assigned atom — Θ*_m(q), the replication of another critical point of 𝒢_q, or a split
configuration. (iii) does not classify these limits: on the instance 𝒢_q has a second local minimum (fit error 0.37), whose replication to
m = 15 has the diagonal plan as its unique optimal plan and is a saddle of J^{OT}_15 — one of its separation blocks is indefinite although
q ∈ W_adm, the admissible set being defined at the branch Θ*(q) (referee round 13). A trajectory that enters the basin of Θ*_m(q) converges to
it: the stationary statement of Proposition D‴.12(iii) is reached dynamically from every initial condition in that basin, with the limit's
masses the plan's masses q.

*Status.* (i) proved (Danskin's theorem and the structure of a minimum of smooth functions; the energy decrease is checked along a run with
reassigned particles, and the two-sided decrease at a constructed tie, in `tests/test_ot_leak.py`). (ii) proved given the imported
Kurdyka–Łojasiewicz theorem for Clarke subgradient trajectories of definable functions [S43] and the definability of restricted analytic
and exponential functions [S44]; the KL exponent is not quantified, so the statement carries no rate — the uniform-in-m rate is Theorem
D‴'s question, which this does not answer. (iii) proved. Under the pure leak an integral plan (commensurate m) never switches once
optimal: Θ(s) = (1 − e^{−σ₀s})T + e^{−σ₀s}Θ(0), and the cost difference between two plans is affine in Θ (the quadratic terms cancel), so it is
the convex combination of its value at Θ(0), nonnegative, and at the targets T, positive because an integral plan costs nothing there and
every other plan costs something for distinct atoms. A split plan can switch even under the pure leak — its barycentric targets are not
atoms, and it need not be optimal at them (referee round 13: switches in 216 of 400 random incommensurate instances, none in 200
commensurate ones). With the fit, cold starts on the instance switch the plan in 2 of 6 runs at m = 15 and in all 6 at m = 16 at the
experiments' σ₀, and in none at ten times that leak — the test's run, which is therefore a smooth gradient flow; the nonsmooth part of (i)
is tested at the static tie.

**Theorem D‴.18 (Theorem J uniformly in m at a replicated configuration).** Let (CL) be the real closed loop of Theorem J with a leak
whose targets are class-constant and locally constant near the configurations considered (the σ-modification; the nearest-atom leak with
every class atom inside its Voronoi cell; the OT leak at w = q with ε < Δ_ν/2, Proposition D‴.12(iii)), and let 𝒲 be a compact set of
weight vectors with min_k w_k ≥ w_min > 0 such that for every w ∈ 𝒲 the weighted stationary point Θ*(w) of 𝒢_w is a non-degenerate
minimum with margins

    λ_coll(w) := λ_min(H_w) ≥ λ₀ > 0,    λ_sep(w) := min_k λ_min(σ₀I − R^G_k(w)) ≥ λ₀ > 0,

H_w the whitened, weighted collective block and σ₀I − R^G_k(w) the separation blocks of Proposition D‴.5(i) (so 𝒲 ⊂ W_adm). For a
configuration Θ = (ϑ_1, …, ϑ_m) with a class map π₀ whose class sizes are n_k = mw_k, w ∈ 𝒲, write u_i := ϑ_i − ϑ*_{π₀(i)}(w), ū_k for the
class means of u, u_i^⊥ := u_i − ū_{π₀(i)}, and

    |||u||| := max{ (Σ_kw_k‖ū_k‖²_{G⁻¹})^{1/2}, max_i‖u_i^⊥‖_{G⁻¹} },    so that    max_i‖u_i‖_{G⁻¹} ≤ (1 + w_min^{−1/2})|||u|||    and    |||u||| ≤ 2max_i‖u_i‖_{G⁻¹}

(the second constant is sharp up to the factor 1 − 1/max_kn_k). There are k₀, ρ, C, κ depending on the plant, the reference, the features,
σ₀, G, 𝒲, λ₀ and the leak's prior (through its cell and plan margins) — but **not on m** — such that for every m,
every w ∈ 𝒲 with mw ∈ ℕ^{K₀}, every class map π₀ with those class sizes and every k_e ≥ k₀:

(a) the lift of the weighted K₀-atom orbit Γ_w (Corollary D‴.14) is the unique T₀-periodic orbit of (CL) whose parameters satisfy
sup_t max_i‖ϑ_i(t) − ϑ*_{π₀(i)}(w)‖_{G⁻¹} ≤ ρ/2 and whose plant state stays in Theorem J's tube;

(b) every solution with max_i‖ϑ_i(0) − ϑ*_{π₀(i)}(w)‖_{G⁻¹} ≤ ρ/2 and ‖x(0) − x_m(0)‖ ≤ R/2 converges to it, and with u(t) the deviation
from the lifted orbit at the times t = jT₀,

    |||u(jT₀)||| ≤ C e^{−(1 − κ)min{λ_coll, λ_sep} jT₀/k_e} (|||u(0)||| + ‖x(0) − x_p(0)‖/k_e);

(c) its Floquet spectrum is Proposition D‴.15's — the K₀-atom orbit's and each class's ghost exponents with multiplicity n_k − 1 —
and the exponents are (1/k_e)λ_j(H_w) and (1/k_e)λ_j(σ₀I − R^G_k(w)) up to O(k_e⁻²), uniformly in m.

In particular, for k_e ≥ k₀ Corollary D‴.14(ii)'s proviso is replaced by a statement uniform in m: for every commensurate m the lift of
Γ_q is the unique periodic orbit in the per-particle ball of radius ρ/2 around Θ*_m(q), and it attracts that ball. (Theorem J's own orbit,
unique in its unnormalised ball, coincides with the lift wherever Theorem J applies, but its threshold in that norm grows like √m, Remark
D‴.19.) Corollary D‴.7's open point — the m-uniformity of Theorem J's thresholds — is settled in the per-particle norm.

*Proof.* Theorem J's proof with the norm |||·||| in place of the unnormalised ‖·‖_{G_m⁻¹}; what has to be checked is that every constant
of its Steps 1–6 is independent of m in the new norm. Three structural facts do it.

(1) *The plant sees the configuration through the class means, up to second order in the separations.* f̂_Θ = Σ_kw_k n_k⁻¹Σ_{i∈k}Φ̃(·;ϑ_i)
and Σ_{i∈k}u_i^⊥ = 0, so f̂_Θ = Σ_kw_kΦ̃(·;ϑ̄_k) + ½Σ_kw_k n_k⁻¹Σ_{i∈k}∇²Φ̃(·;ϑ̄_k)[u_i^⊥, u_i^⊥] + O(max_i‖u_i^⊥‖³): the first-order variation
of the residual r_Θ, hence of the slaved error e_Θ = εr_Θ(x_m) + O(ε²) of Steps 1–2, is a function of the class means alone, and the
separations enter at second order with a coefficient bounded by sup‖∇²Φ̃‖, whatever m. Steps 1–2 (the tube, the slaved error, the fast
variable ζ) involve Θ through r_Θ and — in the ζ-forcing −∂_Θe_Θ·Θ̇ and the remainder ∂_eF·ζ — through the particles' drifts, which are O(ε)
per particle; by (1) and (2) their constants are uniform in m.

(2) *Each particle's drift depends on the others only through the plant.* The particle law is ϑ̇_i = F(x, e, ϑ_i) − σ(ϑ_i − T_i): the right-hand
side of particle i is a function of (x, e, ϑ_i) and of its target. So the frozen-Θ averaged drift ⟨H⟩_i(Θ) = −G∇𝒱[μ_Θ](ϑ_i) and the
antiderivative u_i(Θ, t) of Step 3 are functions of (ϑ_i, f̂_Θ), bounded and C² in ϑ_i with m-independent bounds, and their derivatives in
ϑ_j (j ≠ i) are O(1/m) and pass through f̂ — that is, through the class means at first order. Writing the second-order remainder of the
drift at the replicated configuration per particle, N_i(u) = O(‖u_i‖²) + O(‖u_i‖·(1/m)Σ_j‖u_j‖) + O((1/m)Σ_j‖u_j‖²) (the second derivatives
∂²_{ϑ_iϑ_i} are O(1), ∂²_{ϑ_iϑ_j} and ∂²_{ϑ_jϑ_j} for j ≠ i are O(1/m); ∂²_{ϑ_jϑ_l} for distinct j, l ≠ i vanishes in ⟨H⟩, whose residual
is linear in the Φ̃'s, and is O(1/m²) in the second-order term H₂, which contains products of f̂), so the remainder is bounded in |||·||| by
an m-independent multiple of |||u|||²; the same holds for the near-identity transformations of Step 3 and for the O(ε³) remainders of (N),
which are built from the same functions.

(3) *The linearisation at the replicated configuration is block diagonal with an m-independent spectrum, and contracts in |||·|||.* By
Proposition D‴.5(i) and (1)–(2), D⟨H⟩(Θ*_m(w)) acts on each separation u_i^⊥ of class k separately, as the matrix GR_k − σ₀I — the same for
every particle of the class, similar to −(σ₀I − R^G_k(w)) whose eigenvalues are ≤ −λ_sep — and on the class means as minus the collective
block, self-adjoint in the weighted norm (Σ_kw_k‖ū_k‖²_{G⁻¹})^{1/2} with eigenvalues ≥ λ_coll; the two families do not mix.
Hence I + εT₀D⟨H⟩(Θ*_m(w)) has |||·|||-operator norm 1 − εT₀min{λ_coll, λ_sep} once εT₀ times the largest eigenvalue is below 1, the
same bound for every m; by (2), D⟨H⟩ is Lipschitz on the |||·|||-ball around Θ*_m(w) in that operator norm with an m-independent constant
C_N, the constant used below (the remainder of (2) is then at most (C_N/2)|||u|||², by integration).

With (1)–(3), Step 4 runs as written in the norm |||·||| for P_Y (P_ζ is unchanged), giving the Lipschitz constant
1 − εT₀min{λ_coll, λ_sep}(1 − κ) of P_Y on the |||·|||-ball of radius ρ, κ = O(ρ + ε), with ρ := min{λ₀/(4C_N), the cell and plan
margins of the leak divided by C(1 + w_min^{−1/2})} — the |||·|||-ball of radius ρ reaches the per-particle distance (1 + w_min^{−1/2})ρ, the
solutions of (b) reach C times that, and the leak's targets must stay constant there; the coupling constants Lip_ζP_Y = O(ε) and
Lip_YP_ζ = O(ε³) are m-independent by (1) and (2). Step 6's 2 × 2 Lipschitz matrix
then makes P a contraction in the adapted product norm on that ball, uniformly in m: its fixed point is unique there, and it is the lift
of Γ_w, which is a fixed point by Corollary D‴.14(i) and lies in the ball for k_e ≥ k₀ because sup_t|||ϑ^Γ(t) − ϑ*(w)||| = O(1/k_e) with an
m-independent constant (it is the weighted K₀-atom distance, Theorem J(a) for (CL_w)). This is (a) and (b): the per-particle ball of
radius ρ/2 lies in the |||·|||-ball of radius ρ (since |||u||| ≤ 2max_i‖u_i‖), which gives the uniqueness set of (a) and the start set of
(b), and the |||·|||-ball of radius ρ lies in the per-particle ball of radius (1 + w_min^{−1/2})ρ, which is what the margins have to cover.
(c) is Proposition D‴.15 with its averaging argument, the ghost matrices A_k being m-independent. ∎

*Remark D‴.19 (which norm).* In the unnormalised norm ‖·‖_{G_m⁻¹} of Theorem J the basin radius ρ₀ is m-independent too — the Hessian of
J_m has the m-independent spectrum of (3) and its Lipschitz constant in that norm is O(1) — but a ball of fixed radius in that norm
contains only configurations whose particles are within ρ₀/√m of their atoms on average: it shrinks per particle. And Theorem J's
threshold is not m-independent in that norm: the lifted orbit sits √m·O(1/k_e) away from Θ*_m(w) there (√m times the weighted K₀-atom
distance), and so does the first averaging correction, so the implicit-function step needs k_e ≥ C√m/ρ₀ — Corollary D‴.14(ii)'s proviso,
and the actual reason for the per-particle norm. That norm is the right one because the dynamics is a mean-field one — each particle feels
its own position and the class means — and it cannot be weakened to a mean-square (W₂-type) ball for the nearest-atom leak: a single
particle pushed into another Voronoi cell changes the cluster weights, and the solution converges to the lift of different weights however
small its mean-square displacement — for m large enough that the shifted weights are admissible. At small m they need not be (at m = 15 on
the matching prior, one particle moved out of a class of five gives the weights (6, 5, 4)/15, whose lift is a saddle, numerics below); the
solution then splits instead, and either way it does not converge to Γ_q. Under the OT leak the plan restores the masses (a transit,
Proposition D‴.12(v)), and the basin at the level of the empirical measure is of a different shape: for the averaged loop Theorem D‴.21
below proves an m-independent energy basin, and its a-posteriori form covers configurations with a fixed fraction of the particles in
other cells (the two sets are not nested, and D‴.18 concerns the real loop).

*Status.* Proved as an adaptation of Theorem J's proof, at Theorem J's level of detail: Steps 1–3 and 5 carry over with m-independent
constants by the structural facts (1)–(2), and Steps 4 and 6 are rerun in the norm |||·||| with (3); the m-uniformity rests on the two-sided
norm equivalence and on the mean-field structure of every map built in Steps 1–4. The inputs are those of Theorem J, Corollary D‴.14(i),
Proposition D‴.15, the margins λ₀ of 𝒲 and the leak's cell and plan margins. It replaces Corollary D‴.14(ii)'s proviso by an m-uniform
statement and settles the spectral and the basin part of Corollary D‴.7's uniformity question in the per-particle norm. Checked:
`tests/test_uniform_basin.py`, on the averaged loop at the prior with centres shifted by +0.15 at the matching masses, where the residual
is not zero and λ_sep = 0.80σ₀ (so the per-particle terms that vanish at the matching prior are exercised): the second-order remainder of
the drift is bounded by C_N|||u|||² with C_N < 0.1 for random per-particle displacements at m = 6, 30, 150, 600, not growing with m; for a
one-particle displacement the displaced particle's remainder is exactly N_a + N_b/m with N_a the same for every m — non-zero (above 10⁻⁵)
here and zero (below 10⁻¹¹) at the matching prior — and every other particle's remainder is exactly O(1/m); a separation perturbation (two
particles of one class moved by ±v) moves only those particles at first order, with the response (GR_k − σ₀I)v, not the leak's −σ₀v, at
m = 6, 60, 300; and |||u(s)||| ≤ Ce^{−μs}|||u(0)||| at μ = min{λ_coll, λ_sep} itself, with one constant C ≤ 1.05 for m = 6, 30, 150.

*Numerics* (`exp_uniform_basin.py` → `results/chaos/uniform_basin.md`). The real loop at k_e = 25 (σ₀ = 0.15σ_crit; RK4 plant with the
control and the reference at the stages, Euler parameters, dt = 5·10⁻³) for 8000 time units — 1273 periods, 320 slow units — with all
systems integrated together. *The rates.* The exponents that test (b) and (c) are Floquet exponents (Part F, by shooting). At the matching
prior the 3-atom orbit's parameter exponents are 2.9355·10⁻⁴, 3.5579·10⁻⁴, 3.9986·10⁻⁴, … against λ_j(H_q)/k_e = 2.9521·10⁻⁴, 3.6155·10⁻⁴,
4.0755·10⁻⁴, … (0.56–3.22 %), and every ghost exponent is 2.830·10⁻⁴ = σ₀/k_e. At the prior with centres shifted by +0.15 the parameter
exponents are within 0.57–3.15 % of λ_j(H_q)/k_e and the ghost exponents equal (σ₀ − λ_j(R^G_k))/k_e to 3·10⁻⁴ relative, the slowest
2.2701·10⁻⁴ (class 2) against λ_sep/k_e = 2.2696·10⁻⁴. *The proof's radius* is small: a sampled lower bound C_N ≥ 0.092 on the
|||·|||-Lipschitz constant of D⟨H⟩ gives ρ ≤ 0.019 at the matching prior and ρ ≤ 0.015 at the shifted one, so the runs below, at per-particle
radii 0.3–2.5, start 31–261 times outside the proved start radius ρ/2: the clean contraction there is evidence that the basin extends to
the cell margins, which the theorem does not claim. (A) The matching prior at the matching masses, where the lift is the exact
representation, a constant orbit: sizes m = 3, 15, 60, 240, every particle displaced from its class atom by ρ in a random direction of the
G⁻¹-unit sphere, the same displacements under both leaks. Inside the cells (ρ = 0.3, 1.0, 2.5 < Δ_ν/2 = 3.13) no particle leaves its cell
and the two leaks give identical runs; d_∞ = max_i min_k‖ϑ_i − ϑ⁰_k‖_{G⁻¹} decays at the fitted rates 2.69·10⁻⁴–2.97·10⁻⁴ for m ≥ 15 against
σ₀/k_e = 2.83·10⁻⁴ (ratios 0.95–1.05), with C(T) := max_{t≤T}d_∞(t)e^{μt}/d_∞(0) = 1.03–1.29 at T = 8000 — still growing at the horizon in
the m = 15 runs at ρ = 0.3 and 1 (1.16 and 1.27), constant from T = 2000 on in the others. The m = 3 runs (one particle per class, no
separation modes) are mixtures of collective modes: at ρ = 1 the end of the run is carried by the slowest one and the fitted rate is
2.948·10⁻⁴ (λ_coll/k_e = 2.952·10⁻⁴), at ρ = 2.5 by the third (amplitudes 0.029, 0.024, 0.086; fitted 4.06·10⁻⁴ against 4.08·10⁻⁴), at
ρ = 0.3 by all three (3.48·10⁻⁴). At ρ = 4.5 some particles start in another cell (1, 3, 13 and 55 for m = 3, 15, 60, 240). The nearest-atom
leak ends at other weights — the counts (2, 0, 1), (6, 5, 4), (23, 19, 18), (87, 70, 83) — and not only through the initial misplacement:
at m = 60 and 240 the counts change during the transient (in 3 and 12 records, the last at t = 420 and 285), and m = 240 starts at exactly q.
Its clouds have not collapsed at the horizon (within-class rms spreads 0.24–0.66 for m ≥ 15), but W² = 12.3, 2.81, 2.30, 1.92 is already
that of the lifts of the realised weights for m = 15, 60, 240 (3.21, 2.35, 1.93): the cross-cluster transport of the weight mismatch. Two of
these weight vectors are inadmissible — (6, 5, 4)/15 and (23, 19, 18)/60 have splitting margins −0.00348 and −0.00161, so their lifts are
saddles — while (87, 70, 83)/240 is admissible (+0.00167) and (2, 0, 1)/3 empties a class (Newton finds no collapsed stationary point near
the prior atoms). The OT leak ends at the counts q at every m (the last count change at t = 2555, 2695, 1610, 1445), with W² = 0.010–0.23
still falling at the leak's rate 2σ₀ per slow unit for m ≥ 15 (log-slopes −0.01411, −0.01408, −0.01417 against 2σ₀ = 0.01415): the transit
is over, and what remains is the leak's contraction. The plan relabels — 4, 15 and 53 particles end outside the cell of their initial
class — so the OT leak restores the masses q, not the class map π₀. (C) Remark D‴.19's mechanism: the lift of q with one particle of class
0 moved 0.6 of the way to atom 1, a mean-square displacement of 0.311 (m = 60) and 0.078 (m = 240). Under the nearest-atom leak W² grows to
0.79 and 0.20, toward the lifts of the admissible weights (19, 21, 20)/60 and (79, 81, 80)/240 (W² 0.8655 and 0.2163); under the OT leak it
falls to 3.0·10⁻³ and 8.2·10⁻⁴. (B) A moving orbit: the OT leak toward the prior with centres shifted by +0.15 (a non-zero residual, so Γ_q
moves and the tracking error along it is 2.55·10⁻⁴), masses q; sizes 15, 60, 240 started at the lift of Θ*(q) displaced by ρ = 0.3, 1.0
and measured against the lift of the 3-atom loop's trajectory. The second-half fits 2.18·10⁻⁴–2.63·10⁻⁴ (ratios 0.96–1.16 to 2.27·10⁻⁴)
are mixtures of separation modes whose gaps exceed the fit window; they are consistent with (1 − κ)min{λ_coll, λ_sep}/k_e, and the Floquet
exponents above are the test of the rate. C(8000) = 1.04–1.23 (still growing at the horizon at m = 15, ρ = 0.3), and every m-system's
tracking error is 2.55–2.56·10⁻⁴ at the end.

**Proposition D‴.20 (the entropic OT leak: a soft mass constraint with a critical excess).** Let ν = Σ_jq_jδ_{ϑ⁰_j}, ε > 0, let π^ε(Θ) be
the entropic plan from μ_Θ to ν for the cost ‖·‖²_{G⁻¹} — the unique minimiser of ⟨C, π⟩ + ε KL(π | μ_Θ ⊗ ν) over the plans with marginals
(1/m, q) — and T^ε_i := mΣ_jπ^ε_{ij}ϑ⁰_j its barycentric projection, and consider

    dϑ_i/ds = G⟨∇_ϑ(r_ΘΦ̃(·;ϑ_i))⟩_ρ − σ₀(ϑ_i − T^ε_i).                                                                (𝔄_m^ε)

(i) *(a smooth gradient flow between two leaks)* (𝔄_m^ε) is the G-gradient flow of J^ε_m := mL̄(μ_Θ) + (σ₀m/2)OT_ε(μ_Θ, ν) — by the
envelope theorem ∂OT_ε/∂ϑ_i = (2/m)G⁻¹(ϑ_i − T^ε_i) — and J^ε_m is real-analytic in Θ (the dual potentials are analytic in the cost: the
semi-dual is strictly concave once the gauge g₁ = 0 is fixed, and the implicit function theorem applies), so there are no plan switches and
Łojasiewicz's theorem gives the convergence of every trajectory — all are bounded, T^ε_i lying in conv supp ν and the fit's force being
bounded — without the nonsmooth machinery of Proposition D‴.17. As ε → 0, T^ε → T̄ wherever the OT plan is
unique (Proposition D‴.12); as ε → ∞, π^ε → μ_Θ ⊗ ν and T^ε_i → Σ_jq_jϑ⁰_j: the leak becomes the σ-modification toward the barycentre of ν.

(ii) *(collapsed configurations)* At a configuration collapsed by class with weights w the entropic plan is class-symmetric,
π^ε_{ij} = π^K_{π₀(i)j}/n_{π₀(i)} with π^K the K × K₀ entropic plan from Σ_kw_kδ_{ϑ_k} to ν, so every particle of class k has the class target
T_k = Σ_j(π^K_{kj}/w_k)ϑ⁰_j, and the replicated configuration is a critical point of J^ε_m iff its atoms are a critical point of the weighted
K-atom energy with the leak (σ₀/2)OT_ε(Σ_kw_kδ_{ϑ_k}, ν).

(iii) *(the separation block)* The Hessian of J^ε_m at such a configuration splits as in Proposition D‴.5(i), with separation blocks

    S^ε_k = σ₀(I − (2/ε)Cov^G_k) − R^G_k,    Cov^G_k := G^{−1/2}Cov_{π^K(·|k)}(ϑ⁰)G^{−1/2},

each of multiplicity n_k − 1: a zero-sum perturbation of a class leaves the column marginals, hence the dual potentials, stationary to first
order, and the conditional plan of a particle responds to its own displacement through its softmax, which moves its barycentric target by
(2/ε)Cov_{π^K(·|k)}(ϑ⁰)G⁻¹ — a positive feedback against the leak's restoring σ₀I.

(iv) *(the critical excess)* When w = q and ε ≪ Δ_ν², the conditional plans charge the other atoms with exponentially small mass and S^ε_k is
the OT leak's σ₀I − R^G_k up to exponentially small terms. Let class k carry an excess, and assume A_k := σ₀I − R^G_k ≻ 0 at Θ*_ε(w) — the
splitting condition of Proposition D‴.5, which fails at large excesses and small ε, where the class is a saddle for the fit's own reason.
Then S^ε_k = A_k − (2σ₀/ε)Cov^G_k has a negative eigenvalue iff

    ε < 2σ₀λ_max(A_k^{−1/2}Cov^G_kA_k^{−1/2}),

twice σ₀ times the largest generalised eigenvalue of the conditional covariance relative to A_k. A_k and Cov^G_k depend on ε through Θ*_ε(w)
and its plan, so the critical regularisation ε_c is the root of λ_min(S^ε_k) = 0 along the branch Θ*_ε(w). When the conditional plan sends
a fraction p to one atom j and 1 − p to its own, Cov^G_k = p(1 − p)Δ̂Δ̂ᵀ with Δ̂ := G^{−1/2}(ϑ⁰_j − ϑ⁰_k), and the criterion is the rank-one
ε < 2σ₀p(1 − p)Δ̂ᵀA_k⁻¹Δ̂. Below ε_c the over-weighted collapsed class is a saddle and splits; its split-off particles leave toward an atom its
conditional plan charges, which need not be an under-weighted one, so a split moves the excess but need not reduce it. Above ε_c the class's
separation block is positive, and the replicated configuration is a local minimum when the collective block and the other classes' blocks
are positive as well; the excess is then absorbed by a bias of the class atom toward T_k.

(v) *(the floor at large m)* For fixed ε the critical excess δ_c(ε) is positive. At small ε and R ≈ 0 the critical conditional fraction is
p_c(ε) ≈ ε/(2‖Δ̂‖²), and the plan then sends away exactly the excess, so δ_c ≈ w_kp_c. The sampled excess of an i.i.d. warm start has the
standard deviation (q_k(1 − q_k)/m)^{1/2}, so the entropic leak equalises the cluster masses only down to the critical excess, and for
m ≫ ((1 − q_k)/q_k)(2‖Δ̂‖²/ε)² it does not equalise them at all. The masses are then conserved, as under the nearest-atom leak, and partly
absorbed by the class biases: the stationary measure depends on the sample again, and the floor returns as Proposition D‴.5's law at the
entropic stationary point. Only the hard constraint ε → 0 removes the floor — on the commensurate sizes; on the others its K₀ − 1 bridge
particles leave a floor of order m⁻² (Proposition D‴.12's numerics).

*Status.* (i) proved (the envelope theorem; the analyticity of the entropic potentials; the limits are the convergence of the entropic plan
to an optimal plan as ε → 0 — the optimal plan of maximal entropy where it is not unique [S45] — and the product-measure limit as ε → ∞).
(ii) proved (uniqueness of the entropic plan and the class symmetry). (iii) proved (the first-order stationarity of the potentials under
zero-sum class perturbations and the derivative of the conditional softmax); checked against the finite-difference Hessian of J^ε_m at
m = 15 in `tests/test_entropic_leak.py`. (iv) proved under the hypothesis A_k ≻ 0 (the generalised-eigenvalue criterion and its rank-one
form); the reading of a split as a transport of the excess is not part of it. (v) formal (the scaling of the sampled excess and the
identification of the dynamics' end state); numerical below. The structure is that of the critical temperature of deterministic annealing
[S46], with the same roles — the particles as the codevectors, the prior's atoms as the data, ε as the temperature — and twice the largest
eigenvalue of the conditional covariance when R = 0. New are the two-sided mass constraint, under which the conditional covariance of a
collapsed class is an effect of its mass excess rather than of a data spread, and the fit's curvature in A_k. Checked:
`tests/test_entropic_leak.py` (the semi-dual Newton plan's marginals and its agreement with log-domain Sinkhorn; the gradient of J^ε_m
against finite differences; the two limits at their rates; the replicated configuration critical with a class-symmetric plan, and the
finite-difference Hessian leaving each class's separation subspace invariant and acting there as I ⊗ (σ₀(I − (2/ε)Cov^G_k) − R^G_k),
entrywise; the over-weighted class of (6, 5, 4)/15, with A₁ ≻ 0 checked, unstable at ε = 16, 17 and stable at ε = 18, 20, 64; every
weighted stationary point to a Newton residual of 10⁻⁹).

*Numerics* (`exp_entropic_leak.py` → `results/chaos/entropic_leak.md`). The instance with the matching prior (masses 1/3) at σ₀ = 0.15σ_crit;
the entropic plan by Newton's method on the semi-dual; every weighted stationary point continued from a converged neighbour and checked to a
residual of 10⁻⁹. *The critical excess.* For w_δ = q + δ(1.5e₁ − 0.5) (the excess of class 1 taken from both other atoms) the branch
Θ*_ε(w_δ), continued in ε from 256 down to 1 (residuals ≤ 8.7·10⁻¹³), changes the sign of the over-weighted class's margin once, at
ε_c = 5.11, 13.5, 23.1, 65.0 for δ = 0.02, 0.05, 0.1, 0.2. There A₁ ≻ 0 (λ_min 0.0054, 0.0050, 0.0051, 0.0020), the collective block is
positive (0.0073–0.0078) and so are the other classes' margins: above ε_c the replicated configuration is a local minimum. At ε_c the class
sends the fraction p = 0.057, 0.15, 0.33, 0.57 of its mass away, about equally to the two under-weighted atoms (conditional plans
(0.9434, 0.0277, 0.0289), (0.8485, 0.0764, 0.0751), (0.6655, 0.1850, 0.1496), (0.4257, 0.3074, 0.2669)), so Cov^G₁ has rank two; the exact threshold
2σ₀λ_max(A₁^{−1/2}Cov^G₁A₁^{−1/2}) equals ε_c to the printed digits, and the rank-one formula with all of p sent to one atom gives 1.22–1.31
times ε_c with atom 2 and 0.70–1.22 with atom 3. With the excess taken from atom 3 alone (w′_δ = q + δ(e₁ − e₃)) ε_c = 6.99, 14.6, 24.3 for
δ = 0.02, 0.05, 0.1. At fixed ε the critical excess is δ_c = 0.0042, 0.0083, 0.016, 0.030, 0.062, 0.14, 0.20 at ε = 1, 2, 4, 8, 16, 32, 64
(w_δ; 0.0035–0.18 for w′_δ). Up to ε = 8 it equals w₁p_c, with p_c within 12% (w_δ) and 26% (w′_δ) of the small-ε scale ε/(2‖Δ̂‖²), as
(v) says; from ε = 16 on the plan no longer sends away exactly the excess.
*Dynamics.* m = 15, 60, 240 against M = 510 from the warm clouds of Proposition D‴.12's numerics, four seeds, ε = 1, 4, 16, 64, 600 slow
units, and the ε = 4, 16, 64 runs continued to 6000 slow units. At ε = 1 every system equalises its masses (all sixteen at q by s = 600; fit
errors 1.7–3.0·10⁻⁴; function-level distances to the reference 3.7·10⁻⁹–1.6·10⁻⁶ per seed on (450, 600], values at the stopping time, still
falling): the OT leak's behaviour. At ε = 4 the masses move toward q slowly: at s = 600 five of the twelve larger systems sit at unstable
weights (margins −0.00035 to −0.00146 at Θ*_4(w), e-folding times 1/|margin| of 700–2900 slow units, longer than the run). Continued, every
system ends at stable weights (positive margins at its s = 6000 weights), the larger ones with residual excesses up to 1/60, 3/240 and 7/510
(‖w − q‖_∞ = 0.0167, 0.0125, 0.0137), of the order of the critical excess at ε = 4 (0.016 for an excess taken from both other atoms, 0.012
from one). At ε = 16 the over-weighted (8, 4, 3)/15 reaches (6, 6, 3) within 600 slow units, and in the continuation the three systems at
unstable weights move: (6, 5, 4)/15 to (5, 6, 4) at s ≈ 3700 — the split-off particle joins atom 2, which was at q, so the excess is moved,
not reduced — (6, 6, 3)/15 to (5, 6, 4) as well, and (25, 17, 18)/60 to (23, 19, 18); every other system keeps its sampled weights, which
are stable. At ε = 64 only (8, 4, 3)/15 moves (to (7, 5, 3), at s ≈ 1400). The floors return as Proposition D‴.5's law at the entropic
stationary point Θ*_ε(q), 0.55 (ε = 16) and 0.74 (ε = 64) of the nearest-atom leak's law; the observed floors on the same clouds are
0.23–0.84 (ε = 16) and 0.47–0.91 (ε = 64) of the nearest-atom leak's, per seed; at ε = 64 the bias toward the barycentre of ν dominates the
fit (fit errors 0.12–0.14). The soft constraint thus interpolates in two ways at once: in ε between the OT leak and the σ-modification toward
the barycentre, and in m between equalisation (small m, large sampled excesses) and conserved sampled masses with a floor of the
nearest-atom kind, smaller by the factor of the law (large m).

**Theorem D‴.21 (a measure-level basin for the OT leak: a Polyak–Łojasiewicz inequality uniform in m).** Let ν = Σ_jq_jδ_{ϑ⁰_j} and
let μ* = Σ_jq_jδ_{ϑ*_j} be a replicated stationary measure of (𝔄_m^{OT}) — every atom critical for its class potential,
∇V_L[μ*](ϑ*_j) + σ₀G⁻¹(ϑ*_j − ϑ⁰_j) = 0, and the diagonal plan optimal (ε < Δ_ν/2 suffices, Proposition D‴.12(iii)) — with residual r*.
Put 𝒥(μ) := L̄(μ) + (σ₀/2)W²_{G⁻¹}(μ, ν), so that J_m^{OT}(Θ) = m𝒥(μ_Θ), and 𝒥* := 𝒥(μ*). Let m be commensurate; every vertex optimal
plan is then an assignment π with |π⁻¹(j)| = mq_j (Proposition D‴.12(ii)), J_m^{OT} = min_πJ^π_m over these assignments with J^π_m(Θ) :=
mL̄(μ_Θ) + (σ₀/2)Σ_i‖ϑ_i − ϑ⁰_{π(i)}‖²_{G⁻¹},
and for a configuration Θ and such a π write

    u_i := ϑ_i − ϑ*_{π(i)},    D_π := (1/m)Σ_i‖u_i‖²_{G⁻¹},    Δf̂ := f̂_Θ − f̂[μ*],    E := 𝒥(μ_Θ) − 𝒥*,

and, along the plan segments S_i := {ϑ*_{π(i)} + tu_i : t ∈ [0, 1]},

    κ* := max_i sup_{S_i}λ_max R^G[μ*],    κ*_− := max_i sup_{S_i}(−λ_min R^G[μ*])₊,    κ^Θ := max_i sup_{S_i}λ_max R^G[μ_Θ].

(i) *(two identities along the plan segments)* For every assignment π with the counts mq, optimal or not,

    J^π_m(Θ)/m − 𝒥* = ½‖Δf̂‖² + (σ₀/2)D_π − (1/m)Σ_i∫₀¹(1 − t)u_iᵀR[μ*](ϑ*_{π(i)} + tu_i)u_i dt,                               (E)

    (1/m)Σ_i⟨∇_{ϑ_i}J^π_m(Θ), u_i⟩ = ‖Δf̂‖² + σ₀D_π − (1/m)Σ_i∫₀¹u_iᵀ[(1 − t)R[μ*] + tR[μ_Θ]](ϑ*_{π(i)} + tu_i)u_i dt,             (P)

the norms in L²(ρ). The counts enter through Δf̂ = (1/m)Σ_i[Φ̃(·;ϑ_i) − Φ̃(·;ϑ*_{π(i)})], which holds because |π⁻¹(j)| = mq_j, and through
the constant term of the leak. (E) holds for every coupling of μ_Θ with the column masses q, fractional plans included ((1/m)Σ_i replaced by
Σ_{ij}γ_{ij}); (P), with one test direction per particle, needs an assignment. With counts N ≠ mq both acquire a defect —
Σ_k(N_k/m − q_k)(−⟨r*, Φ̃(·;ϑ*_k)⟩ + (σ₀/2)‖ϑ*_k − ϑ⁰_k‖²_{G⁻¹}) in (E) and −⟨Δf̂, f̂[μ_N] − f̂[μ*]⟩ in (P), μ_N := Σ_k(N_k/m)δ_{ϑ*_k} — and at
incommensurate m no PL inequality relative to 𝒥* can hold: the stationary configurations of Proposition D‴.12(iv), with their bridge
particles, are critical at a positive excess energy.

(ii) *(the Polyak–Łojasiewicz inequality)* For every assignment π with the counts mq, ½‖Δf̂‖² + ((σ₀ − κ*)/2)D_π ≤ J^π_m(Θ)/m − 𝒥* ≤
½‖Δf̂‖² + ((σ₀ + κ*_−)/2)D_π (the lower bound is informative when κ* ≤ σ₀). If π is optimal at Θ, so that the middle term is E, and the
plan-segment margin α := σ₀ − (κ* + κ^Θ)/2 is positive, then

    ‖∇J^π_m(Θ)‖²_{G_m} ≥ λ·mE,    λ := α·min{2, 2α/(σ₀ + κ*_−)}.

(iii) *(a-posteriori decay)* Along every vertex-selection solution of (𝔄_m^{OT}) (§4b, Proposition D‴.17(i): an absolutely continuous
solution that follows an optimal integral plan π(s) for almost every s; their existence is assumed, as in §4b), E is non-increasing and

    E(s) ≤ E(0)exp(−∫₀^s λ(σ)₊dσ),

with λ(σ) computed at Θ(σ) with its plan (λ₊ := 0 where α ≤ 0); where κ* < σ₀, moreover W²(μ_{Θ(s)}, μ*) ≤ D_{π(s)} ≤ 2E(s)/(σ₀ − κ*) and
‖f̂_{Θ(s)} − f̂[μ*]‖² ≤ 2E(s).

(iv) *(an a-priori basin, uniform in m)* Assume (H*): κ_R[μ*](B) ≤ σ₀ − κ for some κ > 0 on a convex set B, forward-invariant for
(𝔄_m^{OT}) and containing the stationary atoms — for instance the whole parameter space, or the a-priori invariant set N_{R_Θ}(conv supp ν)
of Proposition D‴.17 (T̄_i ∈ conv supp ν and the fit's force is bounded). Then for every commensurate m, J_m^{OT} ≥ m𝒥* on B^m, with equality
exactly at the replicated configurations (μ* is the unique global minimiser among the m-point empirical measures in B, up to the
labelling); and every vertex-selection solution from B with E(0) ≤ E₀ := κ²/(2(L₂^G)²) satisfies E(s) ≤ E(0)e^{−λ₀s},
λ₀ := (κ/2)min{2, κ/(σ₀ + sup_B‖R^G[μ*]‖)}. The constants depend on κ, σ₀, L₂^G and r* — not on m.

(v) *(the matching prior; sharpness)* At the matching prior (r* = 0 on supp ρ, μ* = ν), R[μ*] ≡ 0, (H*) holds on the whole parameter space
with κ = σ₀, E = L̄(μ_Θ) + (σ₀/2)W²(μ_Θ, ν), α = σ₀ − κ^Θ/2 and λ = α·min{2, 2α/σ₀}, which tends to 2σ₀ as κ^Θ → 0. The constant 2σ₀ is attained by the
separations (each particle's displacement within its class decays at σ₀, the separation blocks of Theorem D‴.18), so (ii) is sharp there.

(vi) *(what the mass constraint removes)* For the nearest-atom leak the assignment is the Voronoi one, j(·), whose counts N are conserved
inside the cells (Proposition D‴.12(i)) and differ from mq in general. The derivation of (E)–(P) breaks at the counts:
(1/m)Σ_iΦ̃(·;ϑ*_{j(ϑ_i)}) = f̂[μ_N] with μ_N := Σ_k(N_k/m)δ_{ϑ*_k}, not f̂[μ*]. And the nearest-atom energy has a critical point at the
replicated configuration of the weights N/m — a strict local minimum when they are admissible, a saddle when a separation margin is
negative — for every count vector N in the domain of the weighted stationary point w ↦ Θ*(w) (a neighbourhood of q, Proposition D‴.5(i))
whose atoms lie in their own cells, at a positive excess energy when N ≠ mq (under (H*), by (iv)). These spurious minima carry the floor
of Proposition D‴.5. Under the mass constraint they are not critical points; under (H*) μ* is the unique global minimiser in B and the only
critical point where α > 0 (in particular in {E ≤ E₀}). Local minima of J_m^{OT} at a positive excess energy are not excluded: by (ii) they
need κ^Θ ≥ 2σ₀ − κ* on their segments.

*Proof.* (i) The stationarity of the atoms reads σ₀G⁻¹(ϑ*_j − ϑ⁰_j) = ⟨∇_ϑΦ̃(·;ϑ*_j)ᵀr*⟩_ρ. In J^π_m/m − 𝒥* the loss contributes
−⟨r*, Δf̂⟩ + ½‖Δf̂‖² and the leak (σ₀/2)D_π + (σ₀/m)Σ_i⟨u_i, ϑ*_{π(i)} − ϑ⁰_{π(i)}⟩_{G⁻¹} = (σ₀/2)D_π + ⟨r*, (1/m)Σ_i∇_ϑΦ̃(·;ϑ*_{π(i)})u_i⟩.
With the counts, Δf̂ = (1/m)Σ_i[Φ̃(ϑ_i) − Φ̃(ϑ*_{π(i)})], and Taylor's formula with integral remainder along S_i turns the two r*-terms into
−(1/m)Σ_i∫₀¹(1 − t)u_iᵀR[μ*](ϑ*_{π(i)} + tu_i)u_i dt: this is (E). For (P), ∇_{ϑ_i}J^π_m = −⟨∇_ϑΦ̃(·;ϑ_i)ᵀr_Θ⟩ + σ₀G⁻¹(ϑ_i − ϑ⁰_{π(i)})
with r_Θ = r* − Δf̂, and the stationarity gives ⟨∇_{ϑ_i}J^π_m, u_i⟩ = σ₀‖u_i‖²_{G⁻¹} − ⟨r*ᵀ(∇Φ̃(ϑ_i) − ∇Φ̃(ϑ*_{π(i)}))u_i⟩ +
⟨Δf̂ᵀ∇Φ̃(ϑ_i)u_i⟩. The second term is −∫₀¹u_iᵀR[μ*](ϑ*_{π(i)} + tu_i)u_i dt. Averaged over i, the third is ‖Δf̂‖² plus, by Taylor's formula
for Φ̃(ϑ*_{π(i)}) around ϑ_i, (1/m)Σ_i∫₀¹t·u_iᵀR_{Δf̂}(ϑ*_{π(i)} + tu_i)u_i dt, where R_g := ⟨Σ_kg^k∇²Φ̃^k⟩ and R_{Δf̂} = R[μ*] − R[μ_Θ]
because R is linear in the residual. Collecting the R[μ*]-terms gives (P).
(ii) The t-weights 1 − t and t integrate to ½ each. By (E), E = J^π_m/m − 𝒥* lies between ½‖Δf̂‖² + ((σ₀ − κ*)/2)D_π and
½‖Δf̂‖² + ((σ₀ + κ*_−)/2)D_π, and by (P), Q := (1/m)Σ_i⟨∇_{ϑ_i}J^π_m, u_i⟩ ≥ ‖Δf̂‖² + αD_π. Hence Q ≥ cE with c := min{2, 2α/(σ₀ + κ*_−)},
and Q ≥ αD_π. Cauchy–Schwarz in the duality of G and G⁻¹, mQ ≤ ‖∇J^π_m‖_{G_m}(mD_π)^{1/2}, gives ‖∇J^π_m‖²_{G_m} ≥ mQ²/D_π ≥ mαQ ≥ mαcE.
(iii) J_m^{OT} = min_π′J^{π′}_m ≤ J^{π(s)}_m with equality at Θ(s). At every s where Θ(·) and s ↦ J_m^{OT}(Θ(s)) are differentiable — almost
every s, J_m^{OT} being locally Lipschitz — dJ_m^{OT}(Θ(s))/ds ≤ ⟨∇J^{π(s)}_m, Θ̇(s)⟩ = −‖∇J^{π(s)}_m‖²_{G_m} ≤ −mλ(s)₊E(s), and Grönwall's
lemma applies. The distance bounds are (E)'s lower bound and W² ≤ D_π, the assignment being a coupling of μ_Θ and μ*.
(iv) The segments lie in B by convexity, so κ* ≤ σ₀ − κ, and (E) gives E ≥ ½‖Δf̂‖² + (κ/2)D_π ≥ 0, with equality iff D_π = 0. Since
‖R^G[μ_Θ](ϑ) − R^G[μ*](ϑ)‖ ≤ L₂^G‖Δf̂‖_{L¹(ρ)} ≤ L₂^G(2E)^{1/2}, κ^Θ ≤ κ* + L₂^G(2E)^{1/2}, and α ≥ κ − L₂^G(E/2)^{1/2} ≥ κ/2 when E ≤ E₀; E is
non-increasing and the solution stays in B (forward invariance), so this persists, and (iii) gives the rate.
(v) R[μ*] = ⟨r*∇²Φ̃⟩ = 0. For a separation u at μ* (zero class means), Δf̂ = O(‖u‖²), ∇_{ϑ_i}J^π_m = σ₀G⁻¹u_i + O(‖u‖²) and
E = (σ₀/2)D_π + O(‖u‖³), so ‖∇J^π_m‖²_{G_m}/(mE) → 2σ₀.
(vi) With counts N the identity behind Δf̂ becomes (1/m)Σ_iΦ̃(ϑ*_{j(ϑ_i)}) = f̂[μ_N]; the replicated configuration of the weights N/m is
critical for the nearest-atom energy by Proposition D‴.9, a strict local minimum when the weights are admissible and a saddle when a
separation margin is negative (Proposition D‴.5), and its excess energy is positive unless N = mq. Under (H*) a critical point with E > 0
has α ≤ 0 by (ii), that is κ* + κ^Θ ≥ 2σ₀. ∎

*Remark D‴.22 (the mean-field form, propagation of chaos, and the real loop).* The identities and their proof hold verbatim for a
probability measure μ whose optimal plan to ν is induced by a map π (for instance an absolutely continuous μ: the Laguerre map of
Proposition D‴.12), with (1/m)Σ_i replaced by ∫dμ. So formally the Laguerre-cell flow (𝔄₀^{OT}) has the same PL inequality, and in the
basin both the particle system and the mean-field flow converge to μ* exponentially in W₂, at the same m-independent rate, so
W₂(μ_{Θ(s)}, μ_s) ≤ W₂(μ_{Θ(s)}, μ*) + W₂(μ_s, μ*) → 0: for every commensurate m the particle system and the mean-field flow have the same
limit — no floor, where Theorem D‴'s synchronous coupling kept the transit's mismatch for ever (Proposition D‴.12(v)). This is not
propagation of chaos uniform in time, which needs sup_sEW₂(μ_{Θ(s)}, μ_s) → 0 as m → ∞: the bound is not small in m at any fixed time (an
i.i.d. warm start lies in the a-priori basin only at astronomically large m, numerics below), and smallness in m during the transit is not
claimed; uniform-in-time statements of that kind are known for entropy-regularised mean-field dynamics under a uniform log-Sobolev
inequality [S50]. The existence of the mean-field OT flow and its chain rule are not established here. For the real loop: once E(s) ≤
(σ₀ − κ*)(ρ/2)²/(2m), every particle is within ρ/2 of its atom (max_i‖u_i‖² ≤ mD_π), and Theorem D‴.18 takes over; on the averaged
dynamics this happens after a time O(λ⁻¹log m). The transit of the real loop is covered only on finite horizons, by averaging with plan
switches (formal).

*Status.* (i), (ii), (iv), (v) proved. (iii) proved for vertex-selection solutions, whose existence is assumed (§4b); for this leak a
two-particle tie is never attracting — with h := ⟨ϑ_i − ϑ_{i′}, ϑ⁰_k − ϑ⁰_j⟩_{G⁻¹}, J^{π₁} − J^{π₂} = σ₀h, and on the two sides of the tie
dh/ds = ℓ − σ₀‖ϑ⁰_k − ϑ⁰_j‖²_{G⁻¹} and ℓ + σ₀‖ϑ⁰_k − ϑ⁰_j‖²_{G⁻¹} with the same loss term ℓ, so no solution slides along it — while multi-way
ties are not analysed; a Clarke solution in the sense of Proposition D‴.17 may stop at a kink-type Clarke-critical point where the active
plans' gradients cancel (two particles at one point between two atoms), which a selection solution leaves. (vi) proved for the critical
point; its minimality is Proposition D‴.5's admissibility. Remark D‴.22 is formal. The classical case is the leak alone: the metric slope of
½W²(·, ν) is W(·, ν) at every μ, so it satisfies a PL inequality with constant 2 (AGS §10.4, [S47]); a particle velocity field realises the
slope where the optimal plan is a map. The plan segments are the generalised geodesic from μ* to μ_Θ based at ν (AGS §9.2): the diagonal plan
ν → μ* and the assignment ν → μ_Θ make (1 − t)ϑ*_{π(i)} + tϑ_i that curve, along which ½W²(·, ν) is 1-convex — the (σ₀/2)D_π of (E) — and
the passage from convexity along a curve ending at a point to a slope bound is AGS's. PL inequalities in Wasserstein space through optimal
maps are known for functionals that are not geodesically convex (barycentre functionals; the distance of uniform point clouds to a density)
[S48]. Here 𝒥 is semiconvex along generalised geodesics only with a negative constant (the loss contributes −L₂^G sup⟨|r|⟩), so a global
Łojasiewicz argument, as for the λ-convex functionals of [S48], is unavailable. What is specific is the loss linearised exactly in f̂ with
the curvature taken at the stationary residual, the explicit local constants, and the counts, which make the generalised geodesic a particle
interpolation and the constant independent of m. (ii) is sharp at a matching prior and not in general: it bounds E with σ₀ + κ*_− and Q with
(κ* + κ^Θ)/2 separately (numerics below). Checked: `tests/test_ot_basin.py` ((E) and (P) to 10⁻¹² relative for random assignments with the
counts mq at both priors, and with other counts the two defects exactly; the PL chain ‖∇J‖²/(mE) ≥ Q²/(D_πE) ≥ λ; the constant 2σ₀ attained by
a separation to 10⁻⁴; from the lift with one particle moved into another cell the margin positive, the particle kept in its class, E below
the a-posteriori bound and decaying at 2σ₀ at the end, while the nearest-atom leak keeps the floor; the nearest-atom spurious minimum at
(11, 9, 10)/30, admissible — its separation margins and its K-atom Hessian positive — critical for that leak, not critical for the OT leak,
with the PL inequality holding there).

*Numerics* (`exp_ot_basin.py` → `results/chaos/ot_basin.md`). The averaged slow loop at σ₀ = 0.15σ_crit from the warm clouds of Proposition
D‴.12's numerics: sizes m = 15, 60, 240, 510, four seeds, 600 slow units, under both leaks, at the matching prior and at the prior with
centres shifted by +0.15. The identities (E) and (P) hold along every run to 6·10⁻¹³ relative. Under the OT leak the initial optimal plan
is never changed (no particle is relabelled between records), and the margin α is positive from s₊ ≤ 300 on — from s ≤ 25 in 19 of the
32 runs. From then on the observed ratio ‖∇J‖²/(mE) is at least λ at every record (the smallest ratio/λ is 1.00 at the matching prior and
1.61–1.72 at the shifted one, Q²/(D_πE)/λ 1.58–1.69), and E falls from 4.1·10⁻³–4.0·10⁻² to 8.3·10⁻⁷–8.1·10⁻⁶ (means over the seeds;
W²(μ_Θ, μ*) from 0.99–5.25 to 2.3·10⁻⁴–2.3·10⁻³), below the a-posteriori bound in every run. At the end the observed ratio is
0.01416–0.01417 at the matching prior, against 2σ₀ = 0.01415 and the theorem's λ = 0.0140–0.0141: the inequality is sharp. At the shifted
prior λ = 0.0077–0.0078 against an observed 0.0137–0.0139. The local PL constant there is 2min{λ_sep, λ_coll} = 0.01135 (λ_sep = 0.00567),
so the curvature itself costs a factor 0.802 against 2σ₀ and the inequality (ii) another 0.683 (its slack: the theorem's late λ is 0.548
of 2σ₀); the window still mixes the separation modes. *The nearest-atom contrast.* Under the nearest-atom leak no particle changes cell.
At the 22 admissible count vectors (every run at m = 240 and 510, and (5, 4, 6)/15, (19, 21, 20)/60, (19, 18, 23)/60 at both priors) the
runs sit at the replicated configuration of the sampled weights — E(600) is 0.02–0.24 % below its energy and still rising toward it, the
classes collapsed (spreads ≤ 0.006) — and that configuration is a critical point of the nearest-atom energy (drift ≤ 6·10⁻¹⁴) and not of
the OT energy (drift 2.3–4.4·10⁻²): the floor is a spurious minimum. At the other 10, five per prior at m ≤ 60 — (6, 5, 4), (6, 3, 6),
(8, 4, 3)/15 and (24, 15, 21), (25, 17, 18)/60 — the weights are inadmissible (class-0 separation margins −0.0004 to −0.0061): the
replicated configuration is a saddle, class 0 splits (spreads 0.04–0.62), and E ends 0.5–16.5 % below it, still falling. *The a-priori
basin.* At the matching prior κ = σ₀ on the whole space and E₀ = 2.0·10⁻⁷. At the shifted prior the supremum of λ_max R^G[μ*] over the whole
parameter space is 6.60·10⁻³ < σ₀, approached at a saturation corner (a → a_max, s → s_min): (H*) holds on the whole space with
κ = 4.7·10⁻⁴ — μ* is the unique global minimiser among all m-point empirical measures for every commensurate m — and E₀ = 9.1·10⁻¹⁰ (over the
box of the runs' end configurations, which is not forward-invariant, the supremum is 4.39·10⁻³ and E₀ = 2.9·10⁻⁸). These basins are tiny: a
particle outside its cell is at G⁻¹-distance at least 3.13 (matching) and 2.95 (shifted) from its class atom, so {E ≤ E₀} holds a fraction
p ≤ 5.8·10⁻⁶ (matching) and 4.4·10⁻⁷ (shifted) of the particles in other cells — none below m ≈ 1.7·10⁵ and 2.3·10⁶. The runs start at
E(0) = 4.1·10⁻³–4.0·10⁻², far outside, reach the region α > 0 by themselves, and are covered from there a posteriori by (iii). *The
measure-level basin.* With a fraction p of all the particles, taken from class 0, moved 0.6 of the way into cell 1, the plan keeps them in
their class, and κ^Θ = 0.0096, 0.019, 0.029, 0.048 for p = 1/30, 1/15, 1/10, 1/6 (segment suprema on 129 points): the margin is +0.0023 at
p = 1/30 and negative from p = 1/15 on, the same at m = 30 and m = 300, which are the same measures. With one particle moved, E decays from
3.3·10⁻³ to 2.5·10⁻⁵ (m = 30) and from 2.3·10⁻⁴ to 3.1·10⁻⁶ (m = 300), at the rate 2σ₀ at the end and below E(0)exp(−∫λ₊) with λ evaluated at
every step (9.9·10⁻⁵ and 3.7·10⁻⁶; λ ≥ 0.00148 throughout), and the particle returns to its own cell. Under the nearest-atom leak it stays in
cell 1, and E rises toward the spurious minimum of the counts (9, 11, 10)/30 and (99, 101, 100)/300, both admissible: 5.72·10⁻³ and
5.56·10⁻⁴ at s = 300, 6.304·10⁻³ and 6.128·10⁻⁴ at s = 1000, and the minimum's 6.3126·10⁻³ and 6.1328·10⁻⁴ at s = 3000.

**Theorem D‴.23 (adaptive priors: a learned prior removes the leak's bias, at the rate of the measure-space excitation).** Let q be masses,
m commensurate, η ≥ 0, and let the prior's atoms Z = (z_1, …, z_{K₀}) move toward the barycentres of the particles that the plan assigns to
them while the particles follow the OT leak toward ν_Z := Σ_jq_jδ_{z_j}:

    dϑ_i/ds = G⟨∇_ϑ(r_ΘΦ̃(·;ϑ_i))⟩_ρ − σ₀(ϑ_i − z_{π(i)}),    dz_j/ds = −η(z_j − T̂_j),    T̂_j := (1/(mq_j))Σ_{i∈π⁻¹(j)}ϑ_i,      (𝔄_m^{OT,η})

π an optimal assignment from μ_Θ to ν_Z (η = 0 is the fixed prior of Proposition D‴.12).

(i) *(a gradient flow whose replicated stationary points carry no leak bias)* For η > 0, (𝔄_m^{OT,η}) is the vertex-selection gradient
flow of 𝒥(μ_Θ, Z) := L̄(μ_Θ) + (σ₀/2)W²_{G⁻¹}(μ_Θ, ν_Z) in the metric that weighs ϑ_i by (1/m)G⁻¹ and z_j by (σ₀q_j/η)G⁻¹ — by the envelope
theorem ∇_{z_j}W² = 2q_jG⁻¹(z_j − T̂_j) — so 𝒥 is non-increasing (under any schedule η(s) ≥ 0: 𝒥 does not involve η), and bounded
solutions converge to Clarke-critical points (Proposition D‴.17's argument: 𝒥 is a minimum over the vertices of definable smooth functions
of (Θ, Z)). Boundedness is now a hypothesis (Proposition D‴.27 restores it by projecting the prior onto a compact convex set). Proposition D‴.17's a-priori bound rested on the fixed prior's support, while the leak toward an
adaptive prior is an internal force between a class and its atom — with the plan held, the leak cancels from
d/ds[(1/m)Σ_{i∈π⁻¹(j)}ϑ_i + (σ₀q_j/η)z_j] = (1/m)Σ_{i∈π⁻¹(j)}G⟨r∇Φ̃(ϑ_i)⟩ — and solutions can leave every compact set (numerics, Part G(a) and (b)), even from a fit that the network attains exactly.
At a stationary point z_j = T̂_j, and each particle satisfies G⟨r∇Φ̃(ϑ_i)⟩ = σ₀(ϑ_i − T̂_{π(i)}): the learning forces cancel over each
class, but a particle may be held off its class barycentre by the leak. A replicated configuration — the particles of each class π⁻¹(j) at
one point ϑ_j — is therefore stationary with Z iff z_j = ϑ_j and ∇V_L[Σ_lq_lδ_{ϑ_l}](ϑ_j) = 0 for every j: the leak force vanishes, and the
stationary replicated configurations are the critical points of the unregularised weighted K₀-atom loss. (The classes are the plan's:
labelled by the atoms of a given configuration this reads z_{σ(k)} = ϑ_k for the plan's class-to-atom bijection σ, and at an exact
representation a prior permuted among atoms of equal mass is stationary; a permutation between atoms of different masses is not, since a
replicated configuration with W² = 0 needs q_{σ(k)} = q_k.) Non-replicated stationary points — a class held apart against the loss by the
leak toward its barycentre — exist (a saddle on the instance), and carry the K₀-cluster bias of (iv)'s regulariser (numerics, Part G(a)). The stationary equations do not involve the
prior's initial position: an exact representation with the weights q, with the prior at its atoms, is stationary and a global minimum
(𝒥 = 0 there); which stationary point is reached does depend on it ((vi)).

(ii) *(the linearisation at an exact representation)* At an exact representation μ* = Σ_kq_kδ_{ϑ*_k} (r* = 0) with Z = (ϑ*_k), the
linearisation of (𝔄_m^{OT,η}) splits, for every m, into the separations of each class, which decay at σ₀, and the collective pairs (ū, ζ) of
class means and prior offsets. In the whitened weighted coordinates y_k := √q_kG^{−1/2}ū_k, x_k := √q_kG^{−1/2}ζ_k it reads ẏ = −H_Ly − σ₀(y − x),
ẋ = −η(x − y), with the Gauss–Newton block H_L := [√(q_kq_l)G^{1/2}⟨∇_ϑΦ̃(·;ϑ*_k)ᵀ∇_ϑΦ̃(·;ϑ*_l)⟩_ρG^{1/2}]_{kl} ⪰ 0 (eigenvalues
h_1 ≤ … ≤ h_{pK₀}). On its eigenvectors every pair evolves by [[−(h_j + σ₀), σ₀], [η, −η]], whose rates λ_j^± are the roots of
λ² − (h_j + σ₀ + η)λ + ηh_j = 0; for h_j, η > 0 the slower one satisfies

    ηh_j/(h_j + σ₀ + η) < λ_j^− < min{h_j, η},

the lower bound being tight as η → 0; λ_j^− = 0 if h_jη = 0 and λ_j^− = h_j at η = ∞; and λ^− increases with h
(dλ^−/dh = (η − λ^−)/(λ^+ − λ^−) > 0), so the slowest collective rate is λ_1^−. So the adaptive prior's convergence in mode j is never faster
than the loss's own curvature h_j, while a fixed prior's is h_j + σ₀; in an unexcited mode (h_j = 0) λ_j^− = 0, and the network and the prior
drift together, as an adaptive law without a leak does without persistent excitation. The pair matrix is self-adjoint in the metric of (i) (in the coordinates (y_j, (σ₀/η)^{1/2}x_j) it is symmetric, with
off-diagonal (σ₀η)^{1/2}).

(iii) *(the bias it removes)* With a fixed prior near an exact representation, the stationary point is, to first order in the prior's
error x⁰, y* = σ₀(H_L + σ₀I)⁻¹x⁰: mode j keeps the fraction σ₀/(h_j + σ₀) of the error. The leak buys the rate h_j + σ₀ at the price of that
fraction, and the adaptive prior removes it at the rate λ_j^−. Weakly excited modes (h_j ≪ σ₀) keep almost all of the prior's error under a
fixed prior and lose it slowly under an adaptive one; the fit, ½Σ_jh_jy_j² to second order, weighs them by h_j, so the fit is restored faster
than the parameters. The spectrum {h_j} is the measure-space excitation level of the reference at μ* (§8, branch 7).

(iv) *(infinitely fast priors: a capacity-constrained K₀-means regulariser)* Formally, as η → ∞ the prior follows the class barycentres,
z_j = T̂_j, and (𝔄_m^{OT,η}) becomes the gradient flow of L̄(μ_Θ) + (σ₀/2)𝒬_q(μ_Θ), with 𝒬_q(μ) := min_Z W²(μ, ν_Z) the capacity-constrained
K₀-means distortion of μ — along the plan the flow carries, which is optimal for its own barycentres (Lloyd-consistent) and need not be
jointly optimal: the flow follows the selected piece of 𝒬_q, and 𝒬_q itself need not decrease (Proposition D‴.27(ii)) —: a regulariser that pulls the network toward K₀ clusters of the prescribed masses without biasing where they sit.
The separations keep the rate σ₀, and mode j of the class means relaxes at h_j. It does bias how many there are: a target that the network
fits exactly only with a class split is fitted with K₀ clusters (numerics, Part G), so K₀ remains a design input.

(v) *(uniformity in m)* The blocks of (ii) and their rates do not depend on m. In the norm

    |||(u, ζ)|||_η := max{ (Σ_kq_k(‖ū_k‖²_{G⁻¹} + (σ₀/η)‖ζ_k‖²_{G⁻¹}))^{1/2}, max_i‖u_i^⊥‖_{G⁻¹} },

the per-particle norm of Theorem D‴.18 with the prior's offsets added in the metric of (i), the linearisation acts on the separations as
−σ₀I and is self-adjoint on the collective part with the eigenvalues −λ_j^±; the prior is driven by the class barycentres, which are linear
in the class means while the plan is diagonal, so Theorem D‴.18's facts (1)–(2) hold for the joint drift with m-independent constants, and
for 0 < η < ∞ the averaged loop converges, for every commensurate m, in the |||·|||_η-ball of radius

    ρ(η) := min{ λ₀(η)/(4C_N), margin/(C[1 + q_min^{−1/2} + (η/(σ₀q_min))^{1/2}]) },    λ₀(η) := min{σ₀, λ_1^−(η)},

at the rate (1 − κ)λ₀(η) — C_N the Lipschitz constant of the joint drift's Jacobian in the |||·|||_η-operator norm (the prior's part of the
drift is linear, so it is the particles'), margin the plan margin of Proposition D‴.12(iii) at the exact representation, C as in Theorem
D‴.18: the ball reaches the per-particle distance (1 + q_min^{−1/2})ρ and the prior's offsets ‖ζ_k‖_{G⁻¹} ≤ ρ(η/(σ₀q_k))^{1/2}, and the plan
stays diagonal while both are within the margin. At η = ∞ the prior is the class barycentres (ζ = ū), the norm is Theorem D‴.18's,
λ₀ = min{σ₀, h_1}, and the margin term is margin/(C(1 + 2q_min^{−1/2})); for large finite η the bound degenerates like η^{−1/2}, because this
norm does not use the prior's relaxation onto the class means at the rate η (a singular perturbation). As η → 0, λ_1^− < η and the rate
and the radius vanish with it; at η = 0 the prior is a parameter, and the statement is Theorem D‴.18's contraction for the fixed prior.

(vi) *(the transit: when the prior may start)* From a sample whose cluster masses differ from q, the plan's classes contain particles of
other clusters until the transit of Proposition D‴.12 ends, and their barycentres T̂_j are pulled toward those clusters: an adaptive prior is
dragged along, and afterwards relaxes back only at the rates λ_j^−. At η = ∞ the collective modes keep no leak curvature at all, and the
drag can carry the class atoms onto a plateau of the unregularised loss, which the flow then crosses at that loss's own curvature. Freezing
the prior while mass is in transit — η = 0 as long as some particle is assigned to an atom that is not its nearest — keeps the prior out of
the drag; 𝒥 stays non-increasing, as under any schedule of η, and at η = ∞ the opening is a jump of Z to the class barycentres, which lowers
𝒥. But the adaptive prior forgets a perturbation of mode j only at λ_j^− < h_j, so whatever the transient leaves in the weakly excited modes
when the prior starts is kept for a time of order 1/λ_1^−. When the parameters, not only the fit, are wanted at a given horizon, the
adaptation should therefore start from the fixed-prior stationary point — after the transit, and after the fixed-prior dynamics has
settled, which it does at the fixed prior's rates (of order σ₀ and faster) — and then (ii) and (iii) apply with the prior's error as the
only initial displacement. The wait is a delay, which the fit error does not always repay (numerics).

*Proof.* (i) W²(μ_Θ, ν_Z) = min_γΣ_{ij}γ_{ij}‖ϑ_i − z_j‖²_{G⁻¹}; where the optimal plan is unique, ∇_{z_j}W² = Σ_iγ_{ij}2G⁻¹(z_j − ϑ_i) =
2q_jG⁻¹(z_j − T̂_j), so dz_j/ds = −(η/(σ₀q_j))G∇_{z_j}𝒥 and dϑ_i/ds = −mG∇_{ϑ_i}𝒥 (Proposition D‴.12); the vertex-selection argument and the
convergence are Proposition D‴.17's with (Θ, Z) in place of Θ, for solutions that stay bounded. With the plan held,
(1/m)Σ_{i∈π⁻¹(j)}σ₀(ϑ_i − z_j) = σ₀q_j(T̂_j − z_j) = (σ₀q_j/η)ż_j, which is the conservation law. ż = 0 forces z_j = T̂_j for every j, and then
ϑ̇_i = 0 reads G⟨r∇Φ̃(ϑ_i)⟩ = σ₀(ϑ_i − T̂_{π(i)}). At a replicated configuration T̂_j is the class point ϑ_j, so the leak term vanishes and
stationarity is z_j = ϑ_j and ∇V_L(ϑ_j) = 0. At an exact representation L̄ = 0 and W²(μ*, ν_{Z*}) = 0. (ii) At r* = 0, −G∇V_L[μ](ϑ_i) = −G⟨Δf̂ᵀ∇_ϑΦ̃(·;ϑ_i)⟩ with Δf̂ = Σ_lq_l∇_ϑΦ̃(·;ϑ*_l)ū_l + O(‖u‖²) (the separations enter
at second order, Theorem D‴.18's fact (1)), so a separation evolves by −σ₀u^⊥ and the class means by ū̇_k = −GΣ_l⟨∇Φ̃_kᵀ∇Φ̃_l⟩q_lū_l −
σ₀(ū_k − ζ_k), with ζ̇_k = −η(ζ_k − ū_k) because T̂_k − ϑ*_k = ū_k; the coordinates y, x symmetrise it. The eigenvalues of the pair matrix are
−λ with λ² − (h + σ₀ + η)λ + ηh = 0. Their product is ηh and their sum h + σ₀ + η; for h, η > 0 both roots are positive, so λ⁺ < h + σ₀ + η and
λ⁻ = ηh/λ⁺ > ηh/(h + σ₀ + η); the polynomial is −σ₀h < 0 at λ = h and −σ₀η < 0 at λ = η, so both lie strictly between the roots; and
differentiating the polynomial in h at λ⁻ gives dλ⁻/dh = (η − λ⁻)/(λ⁺ − λ⁻) > 0. (iii) The stationary equations of the fixed prior, linearised, read 0 = −H_Ly − σ₀(y − x⁰);
L̄ = ½‖Δf̂‖² = ½yᵀH_Ly + O(‖y‖³). (iv) For a fixed assignment the minimiser of W² over Z is the class barycentres; by the envelope theorem the
particles' drift of L̄ + (σ₀/2)𝒬_q is −G∇V_L − σ₀(ϑ_i − T̂_{π(i)}), the η → ∞ limit of (𝔄_m^{OT,η}) (formal: a singular perturbation in which Z is
the fast variable). (v) (sketch) The symmetrisation of (ii) makes the collective part self-adjoint in |||·|||_η, with the eigenvalues −λ_j^±,
and the separations decouple at −σ₀; the rest is Theorem D‴.18's contraction argument (its fact (3) with this norm) for the averaged loop,
the diagonal plan being locally constant at commensurate m (Proposition D‴.12(iii)) inside the ball, whose per-particle and prior extents are
as stated. The constants κ, C, C_N and the margin are not derived. (vi) 𝒥 does not involve η, so it is non-increasing for any schedule
η(s) ≥ 0, the gates' included; at η = ∞ the opening sets Z to the class barycentres of the current plan, the minimiser of 𝒥 over Z for that
plan, so 𝒥 drops there. ∎

*Status.* (i)–(iii) proved, (i) for the solutions that stay bounded (a hypothesis now); (iv) formal; (v) sketched, as an adaptation of
Theorem D‴.18's contraction argument for the averaged loop in the norm |||·|||_η — the structure (the decoupled separations, the self-adjoint
collective pairs, the m-independent blocks) verified, the radius stated, its constants not derived; (vi) the monotonicity of 𝒥 under the
gates proved (it holds under any schedule of η), the drag, the memory of the weak modes and the remedy numerical (below). The prior that
tracks the network is the continuous-time form of the proximal point method, or iterated Tikhonov regularisation — the leak toward a moving
centre — and the coupled flow of the particles and a centre drawn to their mean is elastic averaging; the fast limit is the K-means form of
soft weight-sharing, a regulariser of network compression [S49], here with the masses fixed (capacity-constrained, [S19]). The per-mode
formulas are those of the quadratic analyses: σ₀/(h_j + σ₀) is the complement 1 − h_j/(h_j + σ₀) of the Tikhonov filter factor in the
Gauss–Newton eigenbasis (the residual factor: the fraction of the prior's error that the regularised solution keeps), iterated Tikhonov
contracts mode j by it per step (in continuous time, the gradient flow of the Moreau envelope, whose curvature h_jσ₀/(h_j + σ₀) gives the
slow-prior limit ηh_j/(h_j + σ₀) of λ_j^−), and the pair matrix is the continuous-time form of elastic averaging's worker–centre map on a
quadratic [S49]. New are their transfer to the mass-constrained particle system — H_L identified as the measure-space excitation spectrum,
and the separations decoupled at σ₀, which make the trade-off between the fixed prior's bias and the adaptive prior's rate hold uniformly in
m — and the transit's drag with its gate. The limits of the bias removal are the precedents' own properties in this form: the lost a-priori
bound is the proximal point method's behaviour without a minimiser (its iterates are bounded iff a zero exists [S49], cited from memory),
and the K₀-cluster bias is the quantisation error of soft weight-sharing's hard limit [S49]; specific here are their particle and
mass-constrained form, and the run-off from a fit that the network attains exactly (Part G(a)). Checked: `tests/test_adaptive_prior.py` (the joint drift against finite differences of 𝒥 in the stated metric; the exact
representation stationary with the prior at its atoms and the fixed prior's stationary point not; the spectrum of the linearisation at
m = 6 equal to the separation rate σ₀ and the pair rates λ_j^± for two η, and the linearisation self-adjoint in the gradient-flow metric;
the fixed prior's bias at a smaller shift to the first-order formula; the transit gate frozen while mass is in transit and the settled
gate opening only when the fixed-prior drift has fallen below its tolerance, with 𝒥 non-increasing in both; the limits of (i): a
non-replicated stationary saddle held by the leak, a prior permuted among atoms of equal mass stationary at the exact representation, and
the adaptive prior's run-off — the force on the amplitude tending to the positive Karush–Kuhn–Tucker multiplier of the loss's minimum on
the face a = a_max, where a bounded limit would need it to vanish, and the run-off law (a_max − a)s → a_max²(σ₀ + η)/(4ηG_αF_∞) of the
conservation law — against the fixed prior's convergence).

*Numerics* (`exp_adaptive_prior.py` → `results/chaos/adaptive_prior.md`). The averaged slow loop at σ₀ = 0.15σ_crit with the OT leak toward a
misplaced prior — the matching prior's centres shifted by +0.15, masses 1/3. *The linear theory.* At the exact representation the whitened
weighted Gauss–Newton spectrum is h_j = 3.06·10⁻⁴ … 0.195. A fixed prior keeps the fraction σ₀/(h_j + σ₀) of its error, 0.959 in the weakest
mode and 0.035 in the strongest, with collective rates h_j + σ₀ ≥ 7.38·10⁻³; the adaptive prior's slowest rate is λ_1^− = 2.68·10⁻⁵,
1.51·10⁻⁴, 2.78·10⁻⁴, 3.06·10⁻⁴ for η = σ₀/10, σ₀, 10σ₀, ∞. *(v)'s radius.* A sampled lower bound on the Jacobian-Lipschitz constant of
the collective block is C_N ≥ 0.341, so the proof's radius is at most 2.0·10⁻⁵, 1.1·10⁻⁴, 2.0·10⁻⁴, 2.2·10⁻⁴ for the four η — λ₀/σ₀ =
0.0038–0.043 times Theorem D‴.18's at the same C_N — while the K₀-atom runs below start at |||·|||_η = 1.17, 0.37, 0.12, 0.021, 94 to 6·10⁴ times
outside, the warm clouds of the m-particle systems 2.2·10⁴–3.6·10⁵ times and the separated starts of Part F 131–3344 times: they probe the
dynamics beyond the proved ball. *The fixed prior's bias.* Θ*(q) has the fit error 9.78·10⁻³, and the first-order
formula y* = σ₀(H_L + σ₀I)⁻¹x⁰ describes it to 4.8% in norm; per mode the weak modes 1 and 3 hold 0.85 and 0.73 of the prediction, and the
error is first order in the shift (halving the shift from 0.15 to 0.01875 takes the norm error through 4.8, 2.3, 1.2, 0.58% and modes 1 and 3
to 0.98 and 0.97 of the prediction). The weakest mode keeps 0.814 of its error, against the predicted 0.959. *The K₀-atom loop* from Θ*(q),
with the prior at the shifted centres, for the four η over 6000 slow units: every mode amplitude above 10⁻⁶ — thirteen of the 36 — decays
at its predicted rate (log-slopes on [3000, 6000] between 0.981 and 1.003 times λ_j^−), at a ratio to the exact pair solution between 0.82
and 1.12 that changes by at most 1.04% from s = 1000 on for η ≥ σ₀ (6.0% at σ₀/10; the offset is set in the first, nonlinear phase); below
10⁻⁶ the faster modes stop at the level of the second-order remainder. The fit error falls from 9.78·10⁻³ to 3.6·10⁻⁴, 1.8·10⁻⁵, 8.0·10⁻⁶,
6.7·10⁻⁶ at s = 6000, while the weakest mode still holds 1.45·10⁻³, 6.7·10⁻⁴, 3.0·10⁻⁴, 2.5·10⁻⁴ of its 1.58·10⁻³: the fit is restored
long before the parameters, as (iii) says. *m-particle systems:* m = 15, 60, 240 and M = 510 from the warm clouds of Proposition D‴.12's
numerics, four seeds, 3000 slow units. Under the fixed prior every system ends at the lift of Θ*(q) (‖ȳ‖ = 2.10·10⁻², W²(μ_Θ, μ_exact) =
4.4·10⁻⁴, fit error 9.78·10⁻³). Ungated, the adaptive prior is dragged by the transit: at s = 500 ‖ȳ‖ is 0.032–0.32 (η = σ₀) and 0.048–0.71
(η = ∞) in seed means, against the K₀-atom run's 7.4·10⁻³ and 4.3·10⁻³, and at s = 3000 still 0.0094–0.10 and 0.012–0.59, the largest at
m = 15. At η = ∞ one m = 15 system is carried onto a plateau of the unregularised loss: its split class re-collapses by s = 7000, onto atoms
where the loss's Hessian has the eigenvalues −1.8·10⁻⁵ and 7.2·10⁻⁵ (fit error 3.46·10⁻² at s = 23000, the drift 1.4·10⁻⁵). The
collapsed η = ∞ dynamics from there, the unregularised K₀-atom gradient flow, falls below 0.9 of the plateau's loss 2.52·10⁻⁴ after 4.16·10⁴
more slow units, halves it at 4.78·10⁴ and falls below 10⁻¹² at 7.97·10⁴ (the crossings located as events), at the exact representation:
the plateau is crossed on the time scale 1/1.84·10⁻⁵ = 5.4·10⁴ of its negative curvature. The transit gate opens at s = 88–149 and removes the plateau, but the transient's
residue in the weak modes stays: ‖ȳ‖ = 4.9·10⁻³–3.7·10⁻² (η = σ₀) at s = 3000, against the K₀-atom run's 1.3·10⁻³, decreasing with m. The
settled gate opens at s = 1211–1304, with the systems replicated (rms separations 1.2·10⁻⁴–7.9·10⁻⁴) and their class means within
6.7·10⁻⁵ of Θ*(q); each is then the K₀-atom system started from its own class means, and at s = 3000 ‖ȳ‖ and the fit error of every
system equal the K₀-atom run's at s − t_open — per-seed ratios 0.995–1.000 (η = σ₀) and 0.968–1.001 (η = ∞) for every m from 15 to 510, as
the invariance of the replicated configurations requires, the departures from 1 being the opening offset in the weakest mode — with fit
errors 1.6·10⁻⁴ (η = σ₀) and 3.4–3.5·10⁻⁵ (η = ∞) against the fixed prior's 9.78·10⁻³, W²(μ_Θ, μ_exact) = 6.7–7.2·10⁻⁶ and 9.3–9.9·10⁻⁷
against 4.4·10⁻⁴, and the counts at q in every run. The wait costs 1211–1304 of the 3000 slow units: in the seed means the settled gate's ‖ȳ‖ and W²
are the lower everywhere, and its fit error at η = ∞ and at m ≤ 240 (at η = σ₀ and m = 510 the transit gate's fit error is the lower,
1.55·10⁻⁴ against 1.64·10⁻⁴); system by system the transit gate ends lower in all three in 6 of the 32 systems. *The structure behind (v), beyond its ball* (Part F): from the lift of Θ*(q) plus zero-mean per-particle
separations of rms 0.02 and 0.1 (m = 15 and 510, ungated, five draws of their direction), ‖ȳ‖ at s = 1750 departs from the K₀-atom run's,
in the draw means, by +0.13% and +3.1% at η = σ₀ and by +0.77% and +18.5% at η = ∞ for m = 15 (+0.09%, +2.2%, +0.63% and +15.5% for
m = 510): second order in the separations draw by draw (factors 23.8–25.0 for the factor 5), with a spread across the draws 3.3–6.0 times
larger at m = 15, and consistent with the same expected size at both m (the draw means, higher at m = 15 in all four columns, differ by
+0.9 to +1.7 combined standard errors; five draws) — the prediction of
Theorem D‴.18's facts (1)–(2) at second order, the departure being linear in the separations' class covariances, whose expectation is the
same at every m. The starts lie 131–3344 times outside (v)'s ball, so (v) itself asserts nothing there. *The limits of (i)* (Part G). On a four-bump target that six kernels fit exactly with one class split
(K₀ = 3, m = 6), the joint flow has a stationary saddle with the class split — the pair at G⁻¹-distance 1.64, against the exact fit's
3.29 — held by a leak force 5.8·10⁻³ against an unregularised gradient of the same size (fit error 4.2·10⁻²); from the exact fit the flow
collapses the class within 8000 slow units (at s = 4·10⁴, η = ∞: the merged width 1.05, fit error 2.23·10⁻²), and the merged atom then
widens toward s_max — continued, its width coordinate runs off like ln s and the width tends to s_max like 1/s, as in (b): the solution is
unbounded — while its physical image approaches the minimiser of the loss on the face s = s_max (fit error 1.943·10⁻², a Karush–Kuhn–Tucker
point with ∂L/∂s = −1.3·10⁻⁶, not a critical point), reached only as s → ∞: the K₀-cluster bias drives the solution from an exact fit to
infinity while 𝒥 decreases. With a target bump of weighted amplitude 4.8 > a_max = 4, the fixed prior's
K₀-atom solution converges (fit error 0.175), while under the adaptive prior (η = σ₀ and ∞) the amplitude coordinate grows like 2 ln s (its
increments per decade 1.99 ln 10 by s = 10⁶) and the amplitude tends to a_max like 1/s: the solution leaves every compact set. At the exact
representation (equal masses) a cyclically permuted prior is stationary (drift 1.1·10⁻¹⁷), with z_σ(k) = ϑ_k.

**Proposition D‴.24 (the finite-size floor is the leak's bias at the sampled weights; an adaptive prior removes it without a mass
constraint).** Let the target admit, for the weights w in a neighbourhood W of q, an exact weighted K₀-atom representation ϑ^e(w) —
f = Σ_kw_kΦ̃(·;ϑ^e_k(w)), a critical point of the unregularised weighted loss whose whitened weighted Gauss–Newton block H_L(w) is
invertible, smooth in w — and let the prior be matching, ϑ⁰ = ϑ^e(q). (On the instance the network's amplitude–mass degeneracy gives
ϑ^e(w) explicitly: the amplitudes A_k/w_k with the bumps' centres and widths, on the branch of ϑ^e(q) while |A_k|/w_k < a_max; these
weights are the *representable* ones, W. They are not Proposition D‴.5's admissible set W_adm of the fixed prior, where its replicated
stationary point is a strict local minimiser; a sample in W_adm is called *D‴.5-admissible*.)

(i) *(the floor is the bias)* Under the nearest-atom leak toward the fixed prior, the weighted stationary point Θ*(w) of Proposition D‴.5
satisfies, with Theorem D‴.23's whitened weighted coordinates at the weights w,

    f̂_{Θ*(w)} − f = Ĵ(w)σ₀(H_L(w) + σ₀I)⁻¹x⁰(w) + O(‖x⁰(w)‖²),    x⁰_k(w) := √w_kG^{−1/2}(ϑ⁰_k − ϑ^e_k(w)),

Ĵ(w) := [√w_k∇_ϑΦ̃(·;ϑ^e_k(w))G^{1/2}]_k (so that ĴᵀĴ = H_L): the leak's bias at the sampled weights, zero at w = q. The function-level
weight-sensitivity of Propositions D‴.5 and D‴.9 is its derivative, F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ at q — exactly, the remainder's derivative
vanishing with x⁰(q) — and on the eigenvectors v_j of H_L(q) (eigenvalues h_j) the function-level floor law of Proposition D‴.9 reads

    E‖f̂_m − f̂_M‖² = (1/m − 1/M)·Σ_j (h_jσ₀²/(h_j + σ₀)²)·z_jᵀC̃z_j + o(1/m),    z := Vᵀ∂_wx⁰,                                       (B)

∂_w taken along the zero-sum directions e_j − e_{K₀} and C̃ the multinomial covariance of (w₁, …, w_{K₀−1}) — Proposition D‴.9's law
written on the eigenvectors, and like it a law for samples in W_adm: at a saddle the replicated stationary point is left by a split
(Proposition D‴.5), and the floor is the split configuration's. Each mode is weighted by h_jσ₀²/(h_j + σ₀)² ≤ σ₀/4, with equality at
h_j = σ₀: a band-pass centred at the level of the leak, at least half its peak for h_j/σ₀ ∈ [3 − 2√2, 3 + 2√2] (a factor 5.8 either side).
The law is therefore at most (1/m − 1/M)(σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ), whatever the excitation spectrum. A weakly excited mode keeps its bias
but hardly shows it in the function; a stiff one shows little bias; which modes carry the floor is set first by the gauge's displacement
∂_wx⁰ and then re-weighted toward h ≈ σ₀.

(ii) *(the nearest-atom leak toward an adaptive prior: the masses minimised out)* Let the prior's atoms follow the means of their
Voronoi cells, dz_j/ds = −η(z_j − T̂_j) (η > 0) with T̂_j the mean of the particles nearest to z_j, while the particles follow the
nearest-atom leak toward Z. Inside the cells this is the gradient flow of

    𝒥_∘(Θ, Z) := L̄(μ_Θ) + (σ₀/2m)Σ_i min_j‖ϑ_i − z_j‖²_{G⁻¹} = min_{q′∈Δ} 𝒥(μ_Θ, Σ_jq′_jδ_{z_j})

in the metric (1/m)G⁻¹ per particle and (σ₀w_j/η)G⁻¹ per atom, w the cell masses: the OT leak of Theorem D‴.23 with the prior's masses
minimised out (semi-relaxed transport — the infinitely fast limit of a dynamics of the masses, not a learning law for them), which sends
every particle to its nearest atom. The counts change only when a particle crosses into another cell; across a crossing the flow is the
vertex-selection flow of 𝒥_∘ — a minimum over the assignments of smooth functions, in a metric that depends on the active assignment — and
𝒥_∘ stays non-increasing, as in Theorem D‴.23(i). A cell without particles has a degenerate metric: its atom is frozen, and the cell is
refilled only if a particle's own dynamics carries it in (on the instance an empty cell stays empty; tests).

(iii) *(no function-level floor at the representable samples)* At fixed counts N, w = N/m, the stationary replicated configurations are
the critical points of the unregularised weighted K₀-atom loss at w (Theorem D‴.23(i) with w in place of q), and for w ∈ W the exact
representation ϑ^e(w) with the prior at its atoms is one, with f̂ = f. At representable sampled weights the function-level sensitivity
therefore vanishes identically, and the function-level distance between two systems at their exact representations is zero. There r = 0,
so R = 0 and Proposition D‴.5's splitting margins are the maximal σ₀ for every w ∈ W: the adaptive prior also removes the fixed prior's
splitting instability, whose replicated point is a saddle wherever a cluster overshoots (on the instance at w₁ ≳ 0.4). Outside W (on the
instance some |A_k|/w_k ≥ a_max) no exact representation exists on the branch of ϑ^e(q) — a permutation of the classes over the bumps can
be representable, but the dynamics from the warm clouds follows the branch — and the adaptive prior removes the leak's bias without
reaching f: the amplitude coordinates run off (Theorem D‴.23(i)), and the sample keeps a function-level error. The probability of such
samples decays exponentially in m (on the instance at most Σ_k exp(−2m(q_k − |A_k|/a_max)²), by Hoeffding's inequality, vacuous at
m = 15), so if every representable sample converges to its exact representation the expected limiting floor is exponentially small in m
rather than of order 1/m; on the instance this asymptotic statement is not yet in force at m ≤ 60 (numerics, Part E). The parameter-level
stationary measure Σ_kw_kδ_{ϑ^e_k(w)} still carries the sample's weights — the counts fix the amplitude–mass gauge — and its coupled
distance is the larger one, conditionally on representability (outside W the parameters run off and the coupled distance diverges):
Proposition D‴.5's coupled distance has the sensitivity S^e := ∂_wϑ^e of the exact representation, where the fixed prior's is, in the
whitened weighted coordinates and to first order, H_L(H_L + σ₀I)⁻¹S^e. Per mode, the weight fluctuation of the exact representation
divides into the fraction h_j/(h_j + σ₀) carried by the parameters and the fraction σ₀/(h_j + σ₀) left as the bias, which the function
sees with the weight h_j; the adaptive prior moves all of it into the parameters. Since S^e = −∂_wx⁰, the bound of (i) is σ₀/4 times the
adaptive prior's coupled-distance law: the fixed prior's function-level law is at most σ₀/4 times the adaptive prior's parameter-level
law. At the measure level the two priors hardly differ: between the stationary measures of two samples, W² is the transport of the mass
½‖w − w′‖₁ across the inter-cluster distances, of order ‖w − w′‖₁ ~ m^{−1/2} under both priors, plus a within-class part of order
‖w − w′‖², the larger under the adaptive prior.

(iv) *(the rates)* At ϑ^e(w) the linearisation is Theorem D‴.23(ii) with w in place of q: the collective pairs relax at
λ_j^−(w) < min{h_j(w), η}, and the fit, weighted by h_j, recovers faster than the parameters. On the boundary of W (an amplitude
A_k/w_k equal to the saturation a_max) the exact representation sits at infinity in the unconstrained coordinates and the convergence is
algebraic: with a = a_max tanh(ã/a_max), a_max − a ≈ 2a_maxe^{−2ã/a_max}, the loss along the saturating coordinate is of order
e^{−4ã/a_max}, its gradient flow gives ã ≈ (a_max/4)log s, and the function-level error decays like a_max²w_k(1 + σ₀/η)/(8G_αs), the loss's
curvature cancelling (formal).

(v) *(the design table)* Fixed prior, nearest-atom leak: the function-level floor (B) at the D‴.5-admissible samples (elsewhere a cluster
splits) and, if the prior is misplaced, its bias. Fixed prior, OT leak (Proposition D‴.12): on the commensurate sizes no function- or
measure-level floor, and the prior's bias. Adaptive prior, nearest-atom leak: at the representable samples no function-level floor and no
bias, and the larger coupled distance of the gauge. Adaptive prior, OT leak (Theorem D‴.23): neither function-level floor nor bias, the
masses fixed by design. Under both OT leaks the coupled distance keeps the transit's mismatch (Proposition D‴.12's *Open*). The mass
constraint fixes the gauge, and the adaptive prior removes the bias — with it the function-level floor, which is its slope in w. This
resolves Remark D‴.10's trade-off between the window's margin and the coupling's floor at the function level: at an exact representation
the separation margin is the maximal σ₀ and both the bias and the function-level floor vanish. The price is the collective rate, λ_j^− in
place of h_j + σ₀ (Theorem D‴.23(ii)), the leak's a-priori bound (restored by projecting the prior, Proposition D‴.27), and a larger
coupled distance.

(vi) *(the real loop: perfect tracking)* With the adaptive prior, the lift of the exact representation of representable sampled weights,
with the prior at its atoms and the plant on the reference (e = 0), is an invariant set of the real closed loop of Theorem J: f̂ = f on the
reference orbit — the support of the occupation measure ρ, which is all the error equation sees at e = 0 — so e ≡ 0 solves the error
equation there (its linear stability also uses ∂_xf̂ = ∂_xf on the orbit, which the exact representation, an identity of functions,
supplies: Corollary D‴.26(c)), and the adaptive law's gradient term (proportional to e), the leak (toward each particle's own position) and the prior's motion
(toward its own atom) all vanish: perfect tracking at every representable sample. With a fixed prior the tracking error keeps Theorem J's
slaved value e = (1/k_e)r_{Θ*(w)}(x_m) + O(k_e⁻²): the leak's bias, divided by the gain.

*Proof.* (i) The derivation of Theorem D‴.23(iii) at the weights w: to first order in u_k := ϑ*_k(w) − ϑ^e_k(w), the stationarity of Θ*(w)
reads H_L(w)y + σ₀(y − x⁰(w)) = 0 in y_k := √w_kG^{−1/2}u_k, and f̂_{Θ*(w)} − f = Σ_kw_k∇_ϑΦ̃(·;ϑ^e_k)u_k + O(‖u‖²) = Ĵy + O(‖u‖²). Since
x⁰(q) = 0, the derivative at q keeps only ∂_wx⁰ (the derivatives of Ĵ, of H_L, of the factors √w_k and of the O(‖x⁰‖²) remainder multiply
x⁰(q) or u(q) = 0), so F is exact at q; the vectors Ĵv_j are orthogonal with ‖Ĵv_j‖² = h_j, which with Proposition D‴.9's law gives (B);
and hσ₀²/(h + σ₀)² ≤ σ₀/4 is (h + σ₀)² ≥ 4hσ₀, with hσ₀²/(h + σ₀)² ≥ σ₀/8 iff h² − 6hσ₀ + σ₀² ≤ 0. (ii) For fixed Z and μ every coupling of μ
with some Σ_jq′_jδ_{z_j} costs at least ∫min_j‖ϑ − z_j‖²dμ, with equality for the nearest-atom map, whose image masses are the cell masses;
inside the cells ∇_{z_j}(σ₀/2m)Σ_{i:j(i)=j}‖ϑ_i − z_j‖²_{G⁻¹} = σ₀w_jG⁻¹(z_j − T̂_j), and ∇_{ϑ_i}𝒥_∘ is Proposition D‴.9's. Across a crossing,
𝒥_∘ is the minimum over the assignments of smooth functions; a vertex-selection solution follows the gradient flow of an active
assignment, whose derivative along it is minus the squared norm of the velocity in that assignment's metric and bounds the upper Dini
derivative of the minimum, so 𝒥_∘ is non-increasing. An empty cell's atom feels no force (w_j = 0) and does not move. (iii) As Theorem
D‴.23(i); at ϑ^e(w) the residual vanishes, so R[μ](ϑ) = ⟨r∇²_ϑΦ̃⟩_ρ = 0 and the separation blocks σ₀I − R^G of Proposition D‴.5 are σ₀I. The
parameter level: two replicated measures with the weights w and w′ and atoms O(‖w − w′‖) apart differ by the mass ½‖w − w′‖₁ that a
coupling must move across inter-cluster distances. The split: ϑ*(w) = ϑ^e(w) + u(w) with ∂_wy = −σ₀(H_L + σ₀I)⁻¹S^e_y at q (S^e_y the
whitened weighted S^e), so the fixed prior's S_y = S^e_y + ∂_wy = H_L(H_L + σ₀I)⁻¹S^e_y; and S^e_y = −∂_wx⁰ by the definition of x⁰. (iv) As
Theorem D‴.23(ii); the boundary case: with ε := a_max − a and α − z_α slaved, the conservation law d/ds[α + (σ₀/η)z_α] = G_α⟨rφ_k⟩∂a/∂α
gives dε/ds = −(η/(σ₀ + η))G_α⟨rφ_k⟩(∂a/∂α)², with ⟨rφ_k⟩ = w_kε⟨φ_k²⟩ and ∂a/∂α ≈ 2ε/a_max, so ε² ≈ a_max²(σ₀ + η)/(8ηG_αw_k⟨φ_k²⟩s) and
⟨r²⟩ = w_k²ε²⟨φ_k²⟩. (vi) Substitution into the closed loop; the slaved error is Theorem J's Step 2. ∎

*Status.* (i): F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ is exact at q, the bias formula at w ≠ q is first order in x⁰(w), and the law (B) is formal — the
delta method, as Proposition D‴.9 — and holds for D‴.5-admissible samples; its value is the asymptotic one (on the instance the exact
floor's expectation is 2.8, 1.5 and 1.1 times it at m = 15, 60, 240). (ii) proved inside the cells, and across crossings for
vertex-selection solutions; (iii) proved for the stationary points, the probability bound by Hoeffding's inequality, the convergence to
the exact representations numerical; (iv) as Theorem D‴.23, the boundary case formal; (v) a summary of (i)–(iii) with Proposition D‴.12 and
Theorem D‴.23; (vi) the invariance proved, its attractivity in the real loop not proved here (it would follow from (iv) and Theorem
D‴.18's argument for the joint system, as sketched in Theorem D‴.23(v)) and checked numerically. The amplitude–mass degeneracy is the
homogeneity that the conic (Wasserstein–Fisher–Rao) view of mean-field networks exploits [49], [55]; the finite-size fluctuations of
trained mean-field networks around their mean-field limit are known to scale like 1/m, with central limit theorems for the training
dynamics ([30], [S51]), and the asymptotic deviation is known to vanish when the mean-field dynamics converges to a measure that
interpolates the data (Chen–Rotskoff–Bruna–Vanden-Eijnden, [S51]). That is the mechanism by which the adaptive prior removes the
function-level floor: its replicated stationary points are unregularised interpolants, and an exact fit is weight-invariant, so the
absence of a function-level floor in (iii) follows at once. (i) is the chain rule applied to Theorem D‴.23(iii)'s first-order bias at
x⁰(q) = 0, and its weight is the squared complement of the Tikhonov filter factor times the mode's visibility h_j, bounded by σ₀/4 by the
AM–GM bound familiar from ridge regression's filter factors. New are the identification of the function-level floor with the
weight-derivative of the regulariser's bias at the sampled weights, with its mode-by-mode law and the bound's identity with the adaptive
prior's coupled-distance law; the split of the weight fluctuation between the parameters and the bias (iii); and a leak toward a learned
prior as the way to restore the interpolation property under regularisation, which also removes the fixed prior's splitting instability.
Checked: `tests/test_adaptive_floor.py` (the nearest-atom adaptive flow against finite differences of 𝒥_∘ with the cell masses; 𝒥_∘ the
minimum of the OT energy over the masses, attained at the cell masses and strictly below the energies with masses moved by ±1/m; the bias
formula against the finite-difference sensitivity, the mode terms summing to the law, the law below the bound, and the weight moving the
shares toward h ≈ σ₀ against the unweighted ones at a small and a large σ₀; at the D‴.5-admissible sample (5, 4, 6)/15 the exact
representation, with the amplitudes A_k/w_k, stationary for the adaptive prior while the fixed prior's stationary point is a strict local
minimiser and biased; at the saddle sample (6, 5, 4)/15 the fixed prior's splitting margin negative and the exact representation's
margins σ₀; the parameter split; 𝒥_∘ decreasing across a cell crossing at η = σ₀ and ∞, and an empty cell staying empty with its atom
frozen).

*Numerics* (`exp_adaptive_floor.py` → `results/chaos/adaptive_floor.md`). The three-bump instance at σ₀ = 0.15σ_crit (and 0.05σ_crit,
0.5σ_crit for the law), the matching prior unless stated, the warm clouds of Proposition D‴.12's numerics (width 0.02, four seeds), m = 15,
60, 240 against M = 510; the checks at the sampled weights from the stored final states. *The floor is the bias* (Part A). The
function-level sensitivity F of the fixed prior's stationary point, by central differences in w, equals the bias form
Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ to the finite-difference error (1.1·10⁻⁶, 1.9·10⁻⁷, 1.2·10⁻⁸ at σ₀/σ_crit = 0.05, 0.15, 0.5, with the step 10⁻⁴); the
law (B) per unit of 1/m − 1/M is 1.356·10⁻², 6.166·10⁻², 0.2291, equal to the sum of its mode terms, at 0.496, 0.753, 0.839 of the bound,
and the bound is σ₀/4 times the adaptive prior's coupled-distance law 46.32. Modes 2–6 carry 99.4%, 99.0%, 97.9% of the law, and 97.2% of
the unweighted sum Σ_jz_jᵀC̃z_j — the concentration is mostly the gauge's — while the weight moves the shares toward h ≈ σ₀: the modes with
h_j/σ₀ ∈ [1/5, 5] carry 42%, 99%, 90% of the law against 21%, 97%, 78% unweighted. The law's value at m = 15, 60, 240 is 3.99·10⁻³,
9.07·10⁻⁴, 1.36·10⁻⁴. The parameter split S = H_L(H_L + σ₀I)⁻¹S^e holds to 1.4·10⁻⁷–2.9·10⁻⁸, and the adaptive prior's coupled-distance
law is 1.44, 2.17, 5.27 times the fixed prior's. Every sampled weight vector is representable (margins min_k(w_k − |A_k|/a_max) =
0.025–0.10) except seed 3's at m = 15, (8, 4, 3)/15, which lies on the boundary (w₃ = |A₃|/a_max). But 5 of the 12 sampled m-systems are
not D‴.5-admissible for the fixed prior: (6, 5, 4)/15, (6, 3, 6)/15, (24, 15, 21)/60, (8, 4, 3)/15 and (25, 17, 18)/60 have splitting
margins −3.48, −1.08, −2.16, −6.14, −3.98·10⁻³, all at class 1 with w₁ ≥ 0.4; the other seven and the four references have margins
≥ 3.78·10⁻³. *Fixed against adaptive priors* (Part B, the nearest-atom leak, 6000 slow units). Under the fixed prior ‖f̂_m − f̂_M‖² stays
at its stationary value, 3.0·10⁻⁴–5.5·10⁻³ at m = 15, 1.1·10⁻⁴–9.9·10⁻⁴ at m = 60 and 4.8·10⁻⁶–8.2·10⁻⁵ at m = 240 over the seeds. At
the seven D‴.5-admissible samples it equals the exact bias at the sampled weights, ‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖², to four digits (ratio
1.0000), while the delta method's ‖F(w^m − w^M)‖² misses it by the factors 0.86–1.72; the five saddle systems have split (class diameters
0.64–1.50 at the end, class 1 broken into sub-clusters) and sit at 0.81–1.01 of the replicated point's floor. Under the adaptive prior the
function-level distance decreases without a floor: at s = 6000 to 4.6·10⁻⁹–2.0·10⁻⁶ (η = σ₀) and 1.9·10⁻¹⁰–5.4·10⁻⁸ (η = ∞) at m = 15,
to 5.4·10⁻¹⁰–6.2·10⁻⁹ and 2.4·10⁻¹⁰–6.2·10⁻¹⁰ at m = 60, and to 2.5·10⁻¹²–1.9·10⁻⁸ and 7.1·10⁻¹⁴–2.2·10⁻⁹ at m = 240, every adaptive
system replicated at the end and the class means converging to the exact representation of the sampled weights. *The test of (iii)* at
the parameter level: the coupled distance D of the shared particles is 0.944–1.000 of the replicated Θ*(w)'s under the fixed prior and
0.951–1.000 of the exact representations' under the adaptive ones (both samples representable), and the adaptive prior's is 1.72–4.81
times the fixed prior's, sample by sample. W²(μ_m, μ_M) is dominated, under both priors, by the cross-class transport of the mass
½‖w^m − w^M‖₁ — 0.82–0.99 of W² under the fixed prior, 0.67–0.99 under the adaptive ones — so the two priors' W² differ by the factors
1.00–1.32 only, while their within-class parts differ by 1.71–4.84. On the boundary sample the function-level distance decays
algebraically (6000·‖f̂_m − f̂_M‖² = 0.88 at η = σ₀ and 0.42 at η = ∞, against the predicted 0.8 and 0.4), and its parameters run off
(D = 6.4 and 7.8 at s = 6000). The counts are conserved in every run. *The design table* (Part C, the misplaced prior, 3000 slow units,
seed means at m = 15, 60, 240). Nearest-atom leak: the fixed prior's fit error 6.8·10⁻², 3.1·10⁻², 2.0·10⁻² and floor 2.5·10⁻³, 5.1·10⁻⁴,
5.9·10⁻⁵ — the same five samples are saddles at the misplaced prior (margins −2.87, −0.38, −1.52, −5.76, −3.41·10⁻³), four split by
s = 3000 and the fifth partly — against the adaptive prior's (η = σ₀) 8.9·10⁻³, 5.5·10⁻⁴, 6.0·10⁻⁴ and 8.0·10⁻⁵, 1.6·10⁻⁷, 1.1·10⁻⁷ (at
m = 15 one sample is on the boundary). OT leak: the fixed prior's 9.78·10⁻³, the bias, with no floor (below 3·10⁻²¹), against the
adaptive prior's 1.6·10⁻⁴ with the floor below 5·10⁻¹¹ (settled gate). At the coupling level the OT designs end with the largest D,
6.0, 3.5, 1.0 at m = 15, 60, 240 — the transit moves 7, 14, 26 of the m-systems' particles into another cell (over the four seeds) —
against the nearest-atom leak's 0.86, 0.25, 0.017
(fixed) and 2.3, 0.48, 0.046 (adaptive), while their W² is at most 2.1·10⁻⁸. The ranking of the fit errors depends on the horizon: until the
settled gate opens (s = 1211–1304) the OT design toward the adaptive prior is the fixed one, and at s = 1250 the ungated nearest-atom
adaptive prior leads at m = 60 and 240 (3.8·10⁻³ and 2.8·10⁻³ against 9.2–9.8·10⁻³), not at m = 15 (1.75·10⁻² in the seed mean); from
s ≈ 2000 the OT design leads (4.0–4.2·10⁻⁴ against
1.2·10⁻², 1.6·10⁻³, 1.4·10⁻³). *The real loop* (Part D, k_e = 25, 3000 time units, seed 0's sampled weights). From the fixed prior's
stationary point the rms tracking error stays at 1.17–1.18·10⁻³ (m = 15) and 3.43–3.44·10⁻⁴ (m = 60), at the end 1.001 and 0.996 times
the averaged loop's (1/k_e)‖f − f̂‖_ρ at s = t/k_e — the m = 15 run on a saddle orbit (Proposition D‴.15; margin −3.48·10⁻³), held there by
the exact replication of its start; with the adaptive prior (η = ∞) it falls to 3.27·10⁻⁴ and 1.52·10⁻⁴ by t = 3000 (s = 120), along the
averaged loop (ratios 0.987 and 0.992); from the exact representation of the sampled weights with the adaptive prior at its atoms it stays
at the integrator's level (the five-period rms 9.0·10⁻¹² at the first window, falling to 1.4–1.5·10⁻¹²): perfect tracking is invariant in
the real loop. *The instance's probabilities* (Part E, exact multinomial sums). Outside W: 0.361, 0.074, 0.0021 at m = 15, 60, 240, and
on its boundary 0.118, 0.053, 0.0012, against Hoeffding's bound 1.87, 0.60, 0.036; a permutation of the classes would be representable
for 31%, 87%, 100% of the mass outside W. The fixed prior's replicated point is a saddle with probability 0.369, 0.260, 0.109 — a failure
set comparable to the adaptive prior's (outside W or on its boundary: 0.48, 0.13, 0.003) at m = 15 and much larger at m ≥ 60. Taking the
constrained minimum of the loss on the branch of ϑ^e(q) as the limit outside W (the dynamics from the warm clouds follows the branch; a
local minimum, the global infimum being 0 wherever a permutation is representable), the adaptive prior's expected
limiting floor is at least 7.3·10⁻³, 2.2·10⁻⁴, 6.6·10⁻⁷: 1.8, 0.24, 0.005 times the fixed prior's law and 0.66, 0.16, 0.004 times its
Monte Carlo expectation. *The law's accuracy* (Part F, 600 nested samples, the exact replicated stationary points). The fixed prior's
E‖f̂_m − f̂_M‖² is 2.78 ± 0.20, 1.53 ± 0.11, 1.11 ± 0.07 times the law (2.89 ± 0.26, 1.56 ± 0.15, 1.07 ± 0.07 for D‴.5-admissible pairs):
the o(1/m) remainder dominates at m = 15, where the gauge's derivative varies fast near the saturation (the matching amplitude 3 against
a_max = 4). The fixed prior's coupled distance follows its law more closely (E D/law 0.68, 0.93, 0.95); the adaptive prior's
(0.53, 1.17, 1.26) has a far heavier tail — its gauge diverges at the boundary of W; at m = 240 one sample lies at 55.7 times the law, and
without it the mean is 1.17 — so its means are not reliable estimates at this sample size; it is 2.9–3.7 times the fixed prior's.

*The planar check* (`exp_adaptive_floor_2d.py` → `results/chaos/adaptive_floor_2d.md`). The planar plant of `exp_ot_2d.py` — four bumps
with vector amplitudes, the Lissajous reference, p = 5, G = diag(1, 1, 0.25, 0.25, 0.25), σ₀ = 0.05 — with the averaged slow loop on ℝ² of
`src/core/planar_loop.py`, the matching prior, m = 16, 32, 64 against M = 128, six seeds, 6000 slow units (η = ∞: the atoms reset to the
cell means at every refresh of the cells, every slow unit, and frozen in between — a lagged prior). The excitation spectrum is
h_j/σ₀ = 0.004–1.41 (twenty modes, all below 1.5σ₀); F equals the bias form to the finite-difference error (6.0·10⁻⁷), and the law per unit
of 1/m − 1/M is 0.310, at 0.215 of the bound, spread over the modes (the largest share 0.18). Representability needs
w > (0.200, 0.150, 0.175, 0.150), whose sum 0.675 leaves little room: 7 of the 18 m-systems are representable, with the six references.
Every sample is D‴.5-admissible for the fixed prior (splitting margins +0.016 to +0.047), and the fixed prior's floor equals the exact bias
at the sampled weights in all 18 systems (ratio 1.0000): 3.2·10⁻⁴–4.2·10⁻³ at the representable systems, 3.0·10⁻³–2.1·10⁻² at the others.
At the representable systems the adaptive prior's falls without a floor, to 1.0·10⁻⁷–4.1·10⁻⁵ (η = σ₀) and 1.2·10⁻⁹–1.0·10⁻⁵ (η = ∞) at
s = 6000 — at most 0.019 and 0.006 of the fixed prior's, still falling in every system — and at the others to 0.023–0.64 of the fixed
prior's, still falling, toward the constrained minimum of the loss at their sampled weights — a floor of their own (Proposition
D‴.24(iii)) — not toward zero. The counts are conserved in every run. In the real loop (k_e = 25, dt = 2·10⁻³, 2500
time units, seed 0, m = 64) the fixed prior's rms tracking error stays at the averaged prediction (ratios 0.993–1.009) and the adaptive
prior's (η = ∞) follows the averaged decay (0.979–1.009); from the exact representation of the sampled weights with the adaptive prior at its atoms
it stays at the integrator's level, 2.0–4.7·10⁻¹⁰.

*Remark D‴.25 (the design table of the leak).* Propositions D‴.9, D‴.12, D‴.20, D‴.24 and Theorem D‴.23 choose, for the averaged loop near a
replicated configuration, what the leak fixes — the prior's atoms and its masses — and what each choice costs. With a matching prior unless
a misplaced one is named, at commensurate sizes for the OT leaks, and with the weights w sampled from an i.i.d. warm start (δw := w^m − w^M):

| leak, prior | function-level floor | coupled distance D (Proposition D‴.5) | W² (measure level) | bias of a misplaced prior | transit | collective rates | a-priori bound |
|---|---|---|---|---|---|---|---|
| nearest atom, fixed (D‴.9) | (1/m − 1/M)Σ_j h_jσ₀²/(h_j + σ₀)² z_jᵀC̃z_j at D‴.5-admissible w (D‴.24(i)) | D‴.5's law with H_L(H_L + σ₀I)⁻¹S^e | the transport of ½‖δw‖₁ between the clusters, plus O(‖δw‖²) | σ₀(H_L + σ₀I)⁻¹x⁰ | none while no particle crosses a cell boundary (D‴.9(i)'s hypothesis) | h_j + σ₀ | yes |
| nearest atom, adaptive (D‴.24) | 0 at representable w | D‴.5's law with S^e (larger), at representable w | likewise, a larger within-class part | 0 | none (likewise) | λ_j^−(w) < min{h_j(w), η} | no (yes when projected, D‴.27) |
| OT, fixed (D‴.12) | 0 | the transit's mismatch | 0 | σ₀(H_L + σ₀I)⁻¹x⁰ | the excess moved by the plan | h_j + σ₀ | yes |
| OT, adaptive (D‴.23) | 0 | the transit's mismatch | 0 | 0 | gated: adapt after the fixed-prior relaxation | λ_j^−(q) | no (yes when projected, D‴.27) |
| entropic OT, fixed (D‴.20) | D‴.9's function-level law at Θ*_ε(q) for sampled excesses below δ_c(ε) | D‴.5's law there | — | that of Θ*_ε(q) | equalises the masses down to δ_c(ε) | — | yes |

Reading. The finite-size floor has three faces: the function level (the sample's weights acting through the leak's bias), the coupling
level (the shared particles displaced by the amplitude–mass gauge, or relabelled by the transit), and the measure level (the mass
½‖δw‖₁ that has to move between the clusters). The mass constraint removes the function-level and the measure-level floors by fixing the
weights; at the coupling level it replaces the gauge's displacement by the transit's mismatch, which relabels different particles in the
two systems. The adaptive prior removes the function-level floor by removing the bias at every representable weight, enlarges the
coupled distance (the
gauge's whole displacement goes into the parameters) and leaves the measure level's leading term unchanged; it also removes the fixed
prior's splitting instability, every separation relaxing at σ₀ at an exact representation. The bias of a misplaced prior is removed only
by an adaptive prior, and every adaptive prior pays with the rate — λ_j^− < h_j, so a weakly excited collective mode converges at the
loss's own curvature rather than at the leak's — and with the leak's a-priori bound on the parameters (Theorem D‴.23(i); restored by
projecting the prior onto a compact convex set, Proposition D‴.27, at the price of a bias where the limit reaches its boundary). The design
without a function-level floor or bias is the OT leak toward an adaptive prior started from the fixed-prior stationary point; the
floor-free design when the prior is trusted is the OT leak toward a fixed one; and when the masses cannot be prescribed (the sample
decides them), the nearest-atom leak toward an adaptive prior keeps the fitted function exact at the representable samples and leaves the
fluctuation in the parameters.

**Corollary D‴.26 (perfect tracking with an adaptive prior in the real loop: the Floquet spectrum, uniformly in m).** Let (CL) be the
real closed loop of Theorem J with the nearest-atom leak toward an adaptive prior — the particles leak at σ = σ₀/k_e toward the atom of
their class, and the atoms follow the class means at the rate η/k_e, ż_j = −(η/k_e)(z_j − ϑ̄_j), 0 < η ≤ ∞ (η = ∞: z_j = ϑ̄_j) — and let
w = N/m be representable sampled weights (w ∈ W, Proposition D‴.24) with the exact representation ϑ^e(w). Let Γ_w be the lift of ϑ^e(w)
by the classes, with Z = ϑ^e(w) and the plant on the reference: a T₀-periodic solution with e ≡ 0 (Proposition D‴.24(vi)). Its Floquet
exponents (the rates −log|λ|/T₀ of its monodromy) are:

(a) σ₀/k_e exactly, with multiplicity (m − K₀)p, on the separation subspaces of the classes;

(b) on the collapsed subspace, the exponents of the weighted K₀-atom loop's periodic solution at ϑ^e(w) — the same for every m with these
weights;

(c) if moreover ∂_xf̂_{ϑ^e(w)} = ∂_xf on the reference orbit (f̂ = f on a neighbourhood of it — on the instances everywhere, the exact
representation being an identity of functions) and x_m ∈ C⁴: n exponents k_e + O(1/k_e) of the error channel, summing to nk_e − tr H_L(w)/k_e +
O(k_e⁻³), and for the K₀p collective pairs

    (1/k_e)[λ_j^±(w) + c_j^±(w)/k_e² + O(k_e⁻⁴)],

λ_j^±(w) the rates of Theorem D‴.23(ii) at the weights w (η = ∞: h_j(w)) and c_j^±(w) the explicit coefficients of the proof (for
simple eigenvalues; at a degenerate one only the rates are even). The odd orders vanish because the slow block at −ε (ε := 1/k_e) is the time reversal
of the reversed reference's at ε, so with the variational equation's symmetry the slow generator's spectrum is even in ε — an identity
(Lemma D‴.26a(iv)) that a derivative mismatch a ≢ 0 breaks; term by term, every term of Ā_k has the parity (−1)^k under reversal, the
first order in particular because the first-order averaging correction is skew. c_j is the sum of the error filter's
attenuation −(dλ_j^±/dh)v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j ≤ 0 (at η = ∞ it is −h_jω_j², ω_j² = Σ_nn²|Ĵ_nv_j|²/Σ_n|Ĵ_nv_j|² the mode's mean squared frequency
in the regressor, by Parseval, independent of η: the averaged rate of a mode holds when k_e ≫ ω_j), of the slaving's part from
−(J̄ᵀJ̄)², which speeds the modes up, and of the averaging terms, of either sign (the level repulsion through Ā₁ is ≥ 0 for the slowest
mode): its sign is not fixed — on the instances every c_j < 0, the filter's part dominating, and with the adaptation gain scaled up
some c_j > 0.

Hence, where H_L(w) ≻ 0 (so λ_1^−(w) > 0) and for k_e above a threshold — which depends on the spectral gaps of the collective block
through Lemma D‴.26a's ε₀ —, perfect tracking is exponentially stable (locally, for each m) at every representable sample, with the
rate (1/k_e)min{σ₀, λ_1^−(w) + c_1^−(w)/k_e² + O(k_e⁻⁴)}; the same holds for the OT leak toward an adaptive prior (Theorem D‴.23) at a
commensurate m, with w = q and the diagonal plan (locally constant at the lift). Without the derivative condition the mismatch
a(t) := ∂_x(f − f̂)(x_m(t)) enters the error equation and moves the collective rates at first order in 1/k_e (the proof).

*Proof.* Along Γ_w the linearisation of (CL) reads δė = (a(t) − k_e)δe − Df̂[δΘ](x_m(t)), δϑ̇_i = G∇_ϑΦ̃(x_m(t); ϑ^e_{π(i)})ᵀδe −
(σ₀/k_e)(δϑ_i − δz_{π(i)}), δż_j = −(η/k_e)(δz_j − δϑ̄_j), with a ≡ 0 under the derivative condition: the second-derivative term
G∇²_ϑ(eᵀΦ̃)δϑ_i and the x-derivative of the learning law vanish because e ≡ 0. (a) A separation (zero class means, δz = 0, δe = 0) does
not move f̂ at first order — Df̂[δΘ] = Σ_k(1/m)∇_ϑΦ̃(·; ϑ^e_k)Σ_{i∈π⁻¹(k)}δϑ_i = 0 — so it stays a separation with δϑ̇_i = −(σ₀/k_e)δϑ_i: the
subspace is invariant and the exponent is σ₀/k_e, whatever a. (b) The collapsed subspace is invariant (Corollary D‴.14(i), with the
prior's equation: the class atoms and their prior atoms follow the weighted K₀-atom equations), and the loop on it is the weighted K₀-atom
loop with its prior. (c) In the whitened weighted coordinates y_k = √w_k G^{−1/2}δϑ_k (and ζ_k likewise for δz_k) the collapsed
variational equation is the adjoint pair δė = (a − k_e)δe − Ĵ(t)y, ẏ = Ĵ(t)ᵀδe − (σ₀/k_e)(y − ζ), ζ̇ = −(η/k_e)(ζ − y), with
Ĵ(t) = [√w_k(G^{1/2}∇_ϑΦ̃(x_m(t); ϑ^e_k))ᵀ]_k (η = ∞: ζ = y, and the leak vanishes on the collapsed subspace). Let ε = 1/k_e, Y = (y, ζ),
J̄ = Ĵ on the y block and C the leak and the prior's constant block in slow time. The error channel is slaved: the linear periodic system
has the invariant subspace δe = L(t)Y, L = εL₁ + ε²L₂ + ε³L₃ + …, with L₁ = −J̄, L₂ = J̄̇ − aJ̄ and, for a ≡ 0, L₃ = −J̄̈ + J̄C − J̄J̄ᵀJ̄
(the Chang transformation, [S9]: Kokotović–Khalil–O'Reilly). This leaves, in slow time s = εt, the slow block
dY/ds = (B₀ + εB₁ + ε²B₂ + …)Y with B₀ = C − J̄ᵀJ̄, B₁ = J̄ᵀ(J̄̇ − aJ̄) and, for a ≡ 0, B₂ = −J̄ᵀJ̄̈ + J̄ᵀJ̄C − (J̄ᵀJ̄)², a system of fast
period εT₀. Its near-identity averaging gives the averaged generator Ā = Ā₀ + εĀ₁ + ε²Ā₂ + …, whose eigenvalues are the slow-time
exponents: Ā₀ = ⟨B₀⟩ = C − ⟨J̄ᵀJ̄⟩, the averaged loop's linearisation at the replicated configuration (Theorem D‴.23(ii) at the weights w,
Proposition D‴.24(iv)), self-adjoint in the metric diag(I, (σ₀/η)I) of Theorem D‴.23(i); Ā₁ = ⟨B₁⟩ + ⟨B̃₀U₁⟩ with U₁′ = B̃₀ (zero means);
Ā₂ = ⟨B₂⟩ + ⟨B₁U₁⟩ + ⟨B₀U₂⟩ with U₂′ = B₁ + B₀U₁ − U₁Ā₀ − Ā₁. With a ≡ 0, Ā₁ is skew: the self-adjoint parts of ⟨J̄ᵀJ̄̇⟩ and ⟨B̃₀U₁⟩ are
½⟨(J̄ᵀJ̄)˙⟩ = 0 and ½⟨(U₁²)′⟩ = 0, and a skew perturbation moves no simple eigenvalue of a self-adjoint matrix at first order (vᵀKv = 0;
on a degenerate eigenspace it splits the eigenvalue into a pair of unchanged real part). All odd orders vanish: the variational matrix
M(t) satisfies M(t)ᵀ = DM(t)D⁻¹ with D = diag(−I_n, I, (σ₀/η)I) (at η = ∞ on (δe, y), D = diag(−I_n, I)), so the system with its
coefficients reversed in time, M(−t), has the
same Floquet multipliers (its monodromy is conjugate to the inverse of the adjoint system's, which is the inverse transpose of the
original's); reversing time changes the sign of every fast-time derivative and antiderivative, and a term of Ā_k carries a number of
them of the parity of k (ε is a time, Ĵ and d/dt are rates, C a rate in slow time), so Ā_k ↦ (−1)^kĀ_k, and the exponents, being
unchanged, are even in ε — for simple eigenvalues of Ā₀; at a degenerate one the first-order splitting ±iεκ is imaginary and only the
rates are even. The symmetry M(t)ᵀ = DM(t)D⁻¹ alone does not suffice: it persists with a ≢ 0, which enters only the (δe, δe) entry, while
a breaks the parity count. (The mechanism is that of the even error expansions of symmetric integrators and of the time symmetry of the
Magnus exponent, [S56], [S52].) Perturbation theory gives c_j = −ℓ_jĀ₂r_j − Σ_{i≠j}(ℓ_jĀ₁r_i)(ℓ_iĀ₁r_j)/(λ_i − λ_j) (r_j, ℓ_j the right and
left eigenvectors of −Ā₀, ℓ_jr_j = 1), and ⟨B₂⟩ contains −⟨J̄ᵀJ̄̈⟩ = ⟨J̄̇ᵀJ̄̇⟩: the error filter k_e/(k_e + iω) attenuates the
regressor's harmonics, H_L − ε²⟨Ĵ̇ᵀĴ̇⟩ + O(ε⁴) in place of H_L; ⟨B₂⟩ also contains −⟨(J̄ᵀJ̄)²⟩, from the slaving's L₃, which enters c_j
with the opposite sign, and the level repulsion −Σ_{i≠j}(ℓ_jĀ₁r_i)(ℓ_iĀ₁r_j)/(λ_i − λ_j) = Σ_{i≠j}|ℓ_jĀ₁r_i|²/(λ_i − λ_j) (Ā₁ skew) is
≥ 0 for the slowest mode, so the sign of c_j is not fixed. The error channel's exponents are k_e + O(ε), and by Liouville's
formula all the exponents sum to −⟨tr M⟩ = nk_e + (σ₀ + η)K₀p/k_e (η = ∞: nk_e), the collective ones to
(1/k_e)tr(−Ā₀) + O(k_e⁻³). With a ≢ 0 the term −⟨J̄ᵀaJ̄⟩ of Ā₁ is self-adjoint and moves the rates at first order, and a — a rate that
does not change sign under time reversal — breaks the parity. ∎

*Status.* (a) and (b) proved, exactly (no averaging); (c) proved by Lemma D‴.26a for x_m ∈ C⁴ — the slaving is an exact
block-triangularisation, the expansion has a bounded O(k_e⁻⁴) remainder at simple eigenvalues (the threshold and the constant depending on
the spectral gaps) and the parity is the exact identity eig Ā(−ε) = eig Ā(ε), the stability statement under H_L(w) ≻ 0 — with the structure
(the skewness of Ā₁, c_j explicit) derived here; the stability consequence local, for each m (the nearest-atom cells are constant near the lift, the particles lying inside
their cells, so the field is smooth there and the multipliers decide). (a) and (b) are Corollary D‴.14(ii) and Proposition D‴.15 at zero residual with the
prior's coordinates added; what the adaptive prior brings is the zero-residual lift at every representable w, not only at w = q under the
OT leak. The basin uniform in m in the per-particle norm would follow from Theorem D‴.18's contraction argument for the joint system, as
sketched in Theorem D‴.23(v), and is not proved here. With a fixed prior the lift of ϑ^e(w) is a solution only when the prior sits at
ϑ^e(w) — at a matching prior, w = q; at the other sampled weights the leak's bias keeps e = O(1/k_e) (Theorem J) — so perfect tracking at
every representable sample, with this spectrum, is specific to the adaptive prior. (a) does not depend on the fixed prior's splitting
margin at w: at (6, 5, 4)/15, where the fixed prior's replicated orbit is a saddle orbit (Proposition D‴.24's numerics), every separation
of perfect tracking decays at σ₀/k_e. Checked: `tests/test_adaptive_floquet.py` (at the matching weights, η = ∞, m = 6: the return of one
period, the nine separation exponents at σ₀/k_e, the collective exponents the K₀-atom loop's, within the averaging correction of the
prediction; at the weights (2, 1, 1)/4 with the prior in the state, η = σ₀, m = 4: the real loop's Jacobian at the lift mapping the
separations to −(σ₀/k_e) times themselves and the collapsed subspace into itself as the collapsed variational equation, the monodromy's
3 separation exponents at σ₀/k_e and its 18 collective ones the variational ones, within the second-order prediction and below the
first-order one; on an asymmetric reference: Ā₁ non-zero and skew, the residual of the second-order prediction falling like k_e⁻⁴, the
measured second-order coefficient the explicit c₁, the error channel by Liouville's formula).

*Numerics* (`exp_adaptive_floquet.py` → `results/chaos/adaptive_floquet.md`). The three-bump instance at σ₀ = 0.15σ_crit, x_m = 2 sin t,
the sampled weights w = (6, 5, 4)/15 (the amplitudes A_k/w_k = 2.5, −2.1, 3.0; representable, and a saddle of the fixed prior, splitting
margin −3.48·10⁻³), η = σ₀ and ∞, k_e = 25 and 50; the monodromy by central differences of the period map (step 10⁻⁵ relative, DOP853 at
rtol 10⁻¹³). Every lift returns after one period to ≤ 1.7·10⁻¹⁵. (a) At m = 15 and 30 the (m − 3)·3 = 36 and 81 separation exponents equal
σ₀/k_e to ≤ 1.7·10⁻⁶ (relative). (b) The other 18 (η = σ₀) or 9 (η = ∞) exponents equal the weighted K₀-atom loop's to 1.4·10⁻⁵–6.4·10⁻⁴
(relative, the resolution of the differences for the slowest of them). (c) Part C integrates the variational equation on the collapsed
subspace over a period (its generator equals the central differences of the real loop's field to ≤ 5.8·10⁻¹¹ relative, its exponents the
differences' to ≤ 4.1·10⁻⁵) at k_e = 12.5, 25, …, 400. At η = σ₀ the slowest exponent lies 37.5%, 13.6%, 3.84%, 0.99%, 0.25% and 0.063%
below λ_1^−(w)/k_e — every exponent lies below its first-order prediction at every k_e, a property of the instance — and the largest
deviation falls by 2.77, 3.53,
3.86, 3.96, 3.99 per doubling of k_e (→ 4), while the residual of the second-order prediction falls by 11.1, 14.0, 15.4, 15.8, 16.0
(→ 16: O(k_e⁻⁴); 4.9·10⁻⁷ at k_e = 400). Every c_j is negative here (c_j/λ_j from −0.47 to −100.6); the slowest mode's c₁/λ₁ = −100.6
(−101.7 at η = ∞), of which the error filter's part is −102.1 = −0.989ω₁² (−103.3 = −ω₁² at η = ∞), ω₁ = 10.16 the slowest mode's
root-mean-square frequency in the regressor — the weakest-excited mode is carried by the regressor's high harmonics, which the error filter attenuates unless k_e ≫ ω₁ —
and the measured (k_eμ₁ − λ₁)k_e²/λ₁ approaches it (−99.4, −100.3, −100.5 at k_e = 100, 200, 400; μ₁ the slowest exponent). On x_m = 2 sin t the first-order correction vanishes identically (in one dimension
⟨ĴᵀĴ̇⟩ is the integral of an exact differential, and the commutator term vanishes by the symmetry t → π − t), so this reference alone
could not test the mechanism; on the asymmetric x_m = 2 sin t + 0.8 sin 2t it does not vanish (‖Ā₁‖ = 5.6·10⁻³, its self-adjoint part
1.4·10⁻¹⁶), and the deviations still fall by → 4 per doubling and the second-order residuals by → 16 (15.5 and 15.9 at the last two
doublings), with c₁/λ₁ = −108.2 (the filter's part −109.7 = −0.984ω₁², ω₁ = 10.56). With the adaptation gain scaled by 30 two
coefficients are positive (c₈/λ₈ = +4.17, c₉/λ₉ = +12.65) and those exponents lie above their first-order predictions by
c_j/(λ_jk_e²) to three digits (+4.17·10⁻⁴ and +1.27·10⁻³ at k_e = 100), the slowest still below. At η = 10σ₀, where the metric
diag(I, (σ₀/η)I) is not the identity, Ā₀ is self-adjoint in it (to 8.7·10⁻¹⁸; in I its asymmetry is 6.4·10⁻²) and Ā₁ skew (to 1.4·10⁻¹⁶).
The error channel's exponent by Liouville's formula is within O(k_e⁻³) of k_e − tr H_L(w)/k_e (24.98084 against 24.98037 at k_e = 25).

**Lemma D‴.26a (the exact slow subspace of the perfect-tracking linearisation: Corollary D‴.26(c) and Proposition D‴.30(ii) proved for
the linearised loop).** Let K: ℝ → ℝ^{n×q} be T₀-periodic and C⁴, Λ ∈ ℝ^{q×q}, a: ℝ → ℝ T₀-periodic and C¹, N₀ ≥ 0, and for real ε ≠ 0
consider

    ε dw = −(w + K(t)Y) dt − (1 + εa(t))ε⁻¹√N₀ dW,    dY = ε(K(t)ᵀw − ΛY) dt + K(t)ᵀ√N₀ dW    (★)

(W a standard n-dimensional Wiener process). At ε = 1/k_e > 0 and N₀ = 0 this is the collapsed variational equation of Corollary D‴.26
under its derivative condition — δe = εw, K = [Ĵ, 0], Y = (y, ζ), Λ = [[σ₀I, −σ₀I], [−ηI, ηI]] for 0 < η < ∞ (at η = ∞: K = Ĵ, Y = y,
Λ = 0; at η = 0, the fixed prior at the atoms: K = Ĵ, Y = y, Λ = σ₀I — with the pair's Λ the prior's block would have a zero rate there) —
and at N₀ > 0 Proposition D‴.30(ii)'s linearised loop under a white measurement noise of density N₀ (the error equation's input
−(k_e + a)ν, a = ∂_xf̂ on the orbit, the law's Ĵᵀν). K ∈ C⁴ holds when x_m ∈ C⁴. There is ε₀ > 0, depending on the C⁴ norms of K, on |Λ|,
‖a‖_{C¹} and T₀ — and, for (iv), on the distance of −λ_j from the rest of the spectrum of Ā₀ and on the condition of its eigenvector, for
(v) on Ā₀'s stability margin, as do the constants of the remainders there — such that for 0 < |ε| ≤ ε₀:

(i) *(the slow subspace)* The equation εṖ + P = −K − ε²P(KᵀP − Λ) has a unique T₀-periodic solution P(·; ε) with ‖P + K‖_∞ = O(ε), and
P = Σ_{k≤3}ε^kP_k + O(ε⁴) uniformly in t, P₀ = −K, P₁ = K̇, P₂ = −(K̈ + KKᵀK + KΛ), P_k = −Ṗ_{k−1} − Σ_{i+j=k−2}P_iKᵀP_j + P_{k−2}Λ. In
the coordinates v := w − P(t)Y the deterministic part of (★) is block-triangular,

    εv̇ = −(I + ε²PKᵀ)v,    Ẏ = εA(t)Y + εKᵀv,    A := KᵀP − Λ = −KᵀK − Λ + εKᵀK̇ − ε²Kᵀ(K̈ + KKᵀK + KΛ) + O(ε³).

(ii) *(the slow multipliers, exactly)* For ε > 0, q Floquet multipliers of (★) are those of the slow block Ẏ = εA(t)Y, of modulus
≥ e^{−εT₀‖A‖_∞}, and the other n those of the fast block, of modulus ≤ e^{−T₀/(2ε)}; the fast exponents sum to
n/ε + ε⟨tr PKᵀ⟩ = nk_e − ⟨tr KᵀK⟩/k_e + O(k_e⁻³).

(iii) *(the slow generator)* Let Φ(ε) be the slow block's monodromy and Ā(ε) := log Φ(ε)/(εT₀) (the principal logarithm). Then
Ā(ε) = Ā₀ + εĀ₁ + ε²Ā₂ + ε³Ā₃ + O(ε⁴) uniformly in |ε| ≤ ε₀, with Ā₀ = −⟨KᵀK⟩ − Λ, Ā₁ = ⟨KᵀK̇⟩ + (1/(2T₀))∫₀^{T₀}∫₀^{t₁}[A₀(t₁),
A₀(t₂)]dt₂dt₁ (A₀ := −KᵀK − Λ) and the higher coefficients those of the Magnus series of εA₍₃₎ collected by powers of ε (A₍₃₎ := KᵀP₍₃₎ − Λ,
P₍₃₎ := Σ_{k≤3}ε^kP_k, a polynomial in ε; A itself has only an asymptotic expansion); the full Magnus series of εA converges where
εT₀‖A‖_∞ < π, a sufficient form of Moan–Niesen's ∫₀^{T₀}‖εA‖₂ < π ([S59]). The averaged generator of Corollary D‴.26(c)'s near-identity
averaging, carried there to Ā₂, is conjugate to Ā(ε) up to O(ε³) — a conjugate generator, with the same eigenvalues to that order; its
coefficients differ from Ā(ε)'s Taylor coefficients by commutators with Ā₀.

(iv) *(the rates, and their parity proved)* At a simple real eigenvalue −λ_j of Ā₀ the slow-time rate μ_j(ε) := −(εT₀)⁻¹log|ν_j(ε)| of
the continued multiplier ν_j is μ_j(ε) = λ_j + εμ_j¹ + ε²μ_j² + ε³μ_j³ + O(ε⁴), the coefficients those of perturbation theory for Ā(ε)
(μ_j¹ = −ℓ_jĀ₁r_j and μ_j² = c_j of Corollary D‴.26(c)). If the full matrix M_ε(t) of (★) satisfies M_ε(t)ᵀ = DM_ε(t)D⁻¹ for an invertible
D — Corollary D‴.26's symmetry, D = diag(−ε²I, I, (σ₀/η)I) in (w, y, ζ) (diag(−ε²I, I) at η = ∞ and at η = 0), which persists with a
mismatch in the error channel — and a mismatch does not enter the error channel (it breaks the parity, below), then the spectrum of Ā(ε) is
even: eig Ā(−ε) = eig Ā(ε) as multisets of complex numbers, the slow block at −ε defined by (i). In particular μ_j(−ε) = μ_j(ε) exactly,
μ_j¹ = μ_j³ = 0 and μ_j(ε) = λ_j + c_jε² + O(ε⁴), with a bounded remainder, and at a degenerate eigenvalue the rates are even.

(v) *(the slow noise)* For 0 < ε ≤ ε₀ the anti-causal equation S − εṠ = εKᵀ − ε²(SPKᵀ + AS) has a unique T₀-periodic solution with
S = εKᵀ + ε²K̇ᵀ + O(ε³), and Y′ := Y + εS(t)v solves exactly

    dY′ = εA(t)Y′ dt + εb(t)√N₀ dW,    b := ε⁻¹[Kᵀ − ε⁻¹(1 + εa)S − εSPKᵀ] = −(K̇ + aK)ᵀ − ε(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ + O(ε²):

the slow coordinates' drift is decoupled from the fast ones, and their noise is εb(t)ν — at leading order the high-passed −ε(K̇ + aK)ᵀν;
the same Wiener process drives v, so E[vY′ᵀ] ≠ 0 — with the diffusion ε²N₀bbᵀ per unit time (in the slow time εN₀⟨bbᵀ⟩ = (N₀/k_e)Q +
O(k_e⁻²), Q := ⟨(K̇ + aK)ᵀ(K̇ + aK)⟩). If Ā₀ is Hurwitz, the periodic stationary second moments are Cov(v) = (N₀/(2ε³))(I + O(ε)),
E[vY′ᵀ] = −N₀bᵀ + O(εN₀) and Cov(Y′) = Σ̄ + O(ε²N₀), Σ̄ the solution of Ā₀Σ̄ + Σ̄Ā₀ᵀ + εN₀Q̄ = 0 (Q̄ = Q on the y block, Σ̄ = O(εN₀)); hence
Cov(Y)(t) − Cov(Y′)(t) = (εN₀/2)K(t)ᵀK(t) − (ε²N₀/2)(d(KᵀK)/dt + 2aKᵀK)(t) + O(ε³N₀) and Cov(Y)(t) = Σ̄ + (εN₀/2)K(t)ᵀK(t) + O(ε²N₀) — the
slow part plus the fast jitter — and Proposition D‴.30(ii)'s floor, (N₀/2k_e)[Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) + tr H²], holds for the
linearised loop under white noise as a period average, with a relative remainder O(1/k_e) (at η = 0 with Λ = σ₀I, whose filter factor is
h_j/(h_j + σ₀)): at time t the fast part is (N₀/2k_e)tr(HĴ(t)ᵀĴ(t)), whose oscillation over the period is of the floor's own order.

*Proof.* (i) On the Banach space of continuous T₀-periodic matrix functions (sup norm) the periodic solution of εu̇ + u = g is
𝒢g(t) = ∫₀^∞e^{−σ}g(t − εσ)dσ — for ε < 0 the anti-causal ∫₀^∞e^{−σ}g(t + |ε|σ)dσ —, so ‖𝒢g‖ ≤ ‖g‖, and 𝒢g − g = −ε𝒢ġ. The equation
is P = 𝒢F(P), F(P) := −K − ε²P(KᵀP − Λ). With r := 2|ε|(‖K̇‖ + β), β := (‖K‖ + 1)(‖K‖(‖K‖ + 1) + |Λ|), the ball ‖P + K‖ ≤ r ≤ 1 is
mapped into itself (‖𝒢F(P) + K‖ ≤ ‖ε𝒢K̇‖ + ε²β) and 𝒢F is a contraction on it with constant ε²L < 1, L := |Λ| + 2‖K‖(‖K‖ + 1): Banach's theorem
gives P, unique in the ball. The truncation P_{(3)} := Σ_{k≤3}ε^kP_k satisfies εṖ_{(3)} + P_{(3)} = F(P_{(3)}) + ρ with ‖ρ‖ = O(ε⁴) (the
recursion cancels the orders ε⁰–ε³; ρ involves K up to its fourth derivative), so P_{(3)} = 𝒢(F(P_{(3)}) + ρ) and
‖P − P_{(3)}‖ ≤ ‖ρ‖/(1 − ε²L) = O(ε⁴). With v = w − PY, εv̇ = −(v + PY) − KY − εṖY − ε²P(Kᵀ(v + PY) − ΛY) = −(I + ε²PKᵀ)v, the
coefficient of Y vanishing by P's equation; and Ẏ = ε(Kᵀ(v + PY) − ΛY). (ii) The change of coordinates is periodic, so it conjugates
the monodromies; in (v, Y) the monodromy is block-triangular, its eigenvalues those of the two diagonal blocks. The fast block's monodromy
has norm ≤ exp(−T₀/ε + εT₀‖P‖‖K‖) ≤ e^{−T₀/(2ε)} (a logarithmic-norm estimate), the slow one's eigenvalues modulus ≥ e^{−εT₀‖A‖}; by
Liouville's formula the fast exponents sum to ε⁻¹⟨tr(I + ε²PKᵀ)⟩, and ⟨tr PKᵀ⟩ = −⟨tr KᵀK⟩ + ε⟨tr K̇Kᵀ⟩ + O(ε²) with
⟨tr K̇Kᵀ⟩ = ½⟨(tr KKᵀ)˙⟩ = 0. (iii) A − A_{(3)} = Kᵀ(P − P_{(3)}) = O(ε⁴) for A_{(3)} := KᵀP_{(3)} − Λ, a polynomial in ε, so by Grönwall's
inequality the monodromies of εA and εA_{(3)} differ by O(ε⁵); the latter, Φ_{(3)}(ε), is analytic in ε with Φ_{(3)}(0) = I, and the
principal logarithm is analytic near I, so log Φ_{(3)}(ε)/(εT₀) is analytic and log Φ(ε)/(εT₀) differs from it by O(ε⁴). Its Taylor
coefficients are obtained by any formal expansion of the logarithm of the monodromy — the Magnus series of εA_{(3)} collected by powers of
ε (Ω₁ = εT₀⟨A⟩, Ω₂ = ½∫∫[εA(t₁), εA(t₂)] gives Ā₁'s commutator) or the near-identity averaging: a periodic transformation I + εU₁ + …
conjugates the system to an autonomous one up to O(ε⁴), so the monodromy is conjugate, by I + εU(0) + …, to exp(εT₀(Ā₀ + … + ε³Ā₃) +
O(ε⁵)). The convergence of the full Magnus series for ∫₀^{T₀}‖εA‖₂ < π is Moan and Niesen's theorem ([S59]). (iv) A simple real eigenvalue
of Ā₀ continues to a simple real eigenvalue of the analytic family log Φ_{(3)}(ε)/(εT₀), analytic in ε (Kato, [S59]), and simple
eigenvalues are Lipschitz under the O(ε⁴) perturbation to Ā(ε): the expansion of μ_j with a bounded remainder, the coefficients by
standard perturbation theory. Parity. Let K̃(t) := K(−t), ã(t) := a(−t). The full matrix of (★) for (K̃, ã) is M_ε(−t), and the system
with the coefficients M_ε(−t) has the multipliers of M_ε(t): its monodromy is conjugate to the transpose of M_ε's, by the symmetry
(ż = −M(s)z is D-conjugate to the adjoint system, whose monodromy is the inverse transpose). The fast multipliers being exponentially
small and the slow ones near 1, the slow multipliers of K̃ at ε are those of K at ε. On the other hand P[K](−t; −ε) solves the equation
of (i) for K̃ at ε (the substitution t ↦ −t, ε ↦ −ε leaves εṖ + P = F(P) invariant when a does not enter it), so by uniqueness
A[K̃](t; ε) = A[K](−t; −ε): the slow block of K̃ at ε is the time reversal of the slow block of K at −ε, and its monodromy is the inverse,
Φ_{K̃}(ε) = Φ_K(−ε)⁻¹, whence Ā[K̃](ε) = Ā[K](−ε). Together, the slow-time exponents of K at −ε are those at ε, and the continued simple
eigenvalue is even: μ_j(−ε) = μ_j(ε) — more, Φ_K̃(ε) = Φ_K(−ε)⁻¹ as matrices and K̃'s slow multipliers are K's, so eig Ā(−ε) = eig Ā(ε)
as multisets (a shorter route needs no P: the full system at −ε is the time reversal of K̃'s at ε, so with the D-symmetry its slow
multipliers are the inverses of K's at ε); the expansion of (iv), valid for both signs with the same coefficients (the recursion of (i) does
not depend on the sign), has no odd terms. A mismatch a in the error channel (εẇ = −(1 − εa)w − KY) adds εaP to F(P), which the substitution
changes in sign: it breaks the parity, at first order. (v) In the coordinates (v, Y), dv = −ε⁻¹(I + ε²PKᵀ)v dt − (ε⁻²(1 + εa) +
PKᵀ)√N₀dW and dY = (εAY + εKᵀv)dt + Kᵀ√N₀dW; for Y′ = Y + εSv the drift's coefficient of v is εKᵀ + εṠ − S − ε²(SPKᵀ + AS), zero by S's
equation, whose periodic solution is the anti-causal contraction as in (i); the noise of Y′ is Kᵀ − εS(ε⁻²(1 + εa) + PKᵀ) =
Kᵀ − (Kᵀ + εK̇ᵀ)(1 + εa) + O(ε²) = −ε(K̇ + aK)ᵀ + O(ε²) times √N₀dW. The second moments solve periodic Lyapunov equations: v's,
V̇ = −ε⁻¹(ΓV + VΓᵀ) + N₀ggᵀ with Γ = I + ε²PKᵀ, g = ε⁻²(1 + εa) + PKᵀ, quasi-static at the fast rate (V = N₀(1 + εa)²/(2ε³) + O(N₀/ε));
the cross moment C = E[vY′ᵀ], Ċ = −ε⁻¹ΓC + εCAᵀ − εN₀gbᵀ, likewise C = −N₀bᵀ + O(εN₀); and Y′'s, Σ̇′ = ε(AΣ′ + Σ′Aᵀ) + ε²N₀bbᵀ, a linear
periodic equation of the slow rate ε: with Σ̄ as stated, Z := Σ′ − Σ̄ is driven by ε[(A − Ā₀)Σ̄ + Σ̄(A − Ā₀)ᵀ] + ε²N₀(bbᵀ − Q̄), of period
average O(ε³N₀) and oscillating part O(ε²N₀); the oscillating part's zero-mean antiderivative is O(ε²N₀), and the slow system's periodic
response to the rest is O(ε⁻¹·ε³N₀) since its Lyapunov operator is invertible at rate ε (Ā₀ Hurwitz): Z = O(ε²N₀). Then
Cov(Y) = Cov(Y′) − ε(SC + CᵀSᵀ) + ε²SVSᵀ, where −ε(SC + CᵀSᵀ) = −ε²N₀(d(KᵀK)/dt + 2aKᵀK) + O(ε³N₀) and, with S = εKᵀ + ε²K̇ᵀ + O(ε³) and
V = N₀(1 + εa)²/(2ε³) + O(N₀/ε), ε²SVSᵀ = (εN₀/2)KᵀK + (ε²N₀/2)(d(KᵀK)/dt + 2aKᵀK) + O(ε³N₀): Cov(Y) = Σ̄ + (εN₀/2)KᵀK + O(ε²N₀). The floor
proved is the period average of tr(H Cov(y)): the slow part by D‴.30(ii)'s pair blocks (at η = 0 by −(H + σ₀I)) and the fast part
(εN₀/2)tr(H⟨ĴᵀĴ⟩) = (N₀/2k_e)tr H²; at time t the fast part is (εN₀/2)tr(HĴ(t)ᵀĴ(t)). ∎

*Status.* Proved (the constants not made explicit). It upgrades Corollary D‴.26(c) from a formal to a proved expansion, with a bounded
O(k_e⁻⁴) remainder and the odd orders vanishing exactly (for x_m ∈ C⁴ at simple eigenvalues; the stability where H_L(w) ≻ 0), and
Proposition D‴.30(ii) to proved for the linearised loop under white noise, for the floor's period average, with the relative remainder
O(1/k_e) — the order of round 21's observation, whose first-order term comes from the slope a in the noise input (with a = 0 in the noise
the numerics show an even expansion; (v) does not prove it) — at η = 0 with the identification K = Ĵ, Λ = σ₀I. The real loop's noise at
second order (D‴.30(iii)), the OU noise of the numerics and the real nonlinear loop stay as stated there. Chang's decoupling transformation
of linear singularly perturbed systems — the pair (P, S), here periodic —, its time-varying and stochastic versions (Kokotović, Khalil and
O'Reilly's Ch. 4–5 and App. A, with the slow subsystem's transformed noise and its covariance; Haddad and Kokotović; Yang and Zhu's
contraction for slowly time-varying systems), the slow manifolds of stochastic slow–fast systems (Berglund and Gentz) and the convergence
of the Magnus series are classical ([S59], [S9]), and so is the reversal argument for even expansions ([S56]); what the lemma adds is
their application to the perfect-tracking linearisation of the transport-leak loop: the slaving of Corollary D‴.26 made an exact
block-triangularisation, the parity made an exact identity eig Ā(−ε) = eig Ā(ε) through the slow block's definition at negative ε, and the
slow noise identified exactly as εb(t)ν, b = −(K̇ + aK)ᵀ − ε(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ + O(ε²) — the high-passed noise at leading order, the
drift decoupled, the noise shared with the fast block.
Checked by `exp_adaptive_floquet.py` Part D: the asymmetric reference at η = 0, σ₀ and ∞, k_e = 12.5–1600 — P's invariance residual, the
truncations against ε³ and ε⁴, the slow block's rates against the full monodromy's, the rates at ±ε per mode, a mismatch breaking their
equality at first order, the Magnus parameter; under white noise the decoupling residual, b against b₀ and b₀ + εb₁, the exact periodic
covariance in the coordinates (v, Y′) with the weakest mode's deviation split into its two sources, the jitter with its next term, and the
floor over the period. Checked separately by `tests/test_adaptive_floquet.py`: the slow block's rates against the monodromy's and P − P₂
against ε³ on the instance; on a synthetic block with an O(1), non-symmetric Λ, P − P₂ against ε³ and P − P₍₃₎ against ε⁴, with a flipped
Λ and a dropped cubic term caught; the parity and its breaking by a mismatch at η = σ₀ and ∞; and b₀ against K̇ by the chain rule and the
plant's own slope by differences, b₁ against ε², and the covariance's orders.

*Numerics* (`exp_adaptive_floquet.py` Part D → `results/chaos/adaptive_floquet.md`). Corollary D‴.26's instance (w = (6, 5, 4)/15) on the
asymmetric reference x_m = 2 sin t + 0.8 sin 2t, η = 0 (K = Ĵ, Λ = σ₀I), σ₀ and ∞, k_e = 12.5–1600, a 512-point Fourier grid. The
contraction for P converges in 3–9 iterations to an invariance residual ≤ 3.8·10⁻¹⁵ (at η = 0 and k_e = 12.5 it stalls at rounding: 499
iterations, the residual 3.3·10⁻¹⁵). |P − P₂|/ε³ tends to 3.34·10³ and |P − P₍₃₎|/ε⁴ to 5.8·10⁴ (bounded: the O(ε³) of (i) and the third
term of (iii), their constants set by the regressor's harmonics). The slow block's slow-time rates equal the full variational equation's
to ≤ 1.0·10⁻⁸ (the monodromies' accuracy). The rates at ±ε agree to ≤ 5.9·10⁻¹² of the largest rate and ≤ 1.1·10⁻⁹ per mode, while
their deviation from the first-order prediction is 9.8·10⁻²–8.9·10⁻⁶ of the largest rate and, per mode, up to 29.6% — the slowest mode's
at η = ∞: −29.6%, −12.5%, −4.0%, −1.1% at k_e = 12.5, 25, 50, 100; at η = 0 the slowest mode's is −1.7% and the largest 12.0% at k_e =
12.5. The parity of (iv) is an identity, not an asymptotic statement. With the plant's slope as a mismatch in the error channel the rates
at ±ε differ by 1.1·10⁻²–7.8·10⁻⁵ of the largest rate (per mode 6.1·10⁻²–3.5·10⁻⁴), halving per doubling of k_e: the parity broken at
first order. The Magnus parameter is 0.48 at k_e = 12.5, ≤ 0.25 from k_e = 25 and 0.004 at 1600, so the full Magnus series converges
throughout. Under white measurement noise the second transformation decouples Y′ to a residual ≤ 3·10⁻¹⁶. The slow noise b differs from
b₀ = −(K̇ + aK)ᵀ by 63%, 40%, 22%, 12%, 5.9%, 3.0%, 1.5% and 0.74% (O(ε), with the large constant of the regressor's harmonics: εω₁ ≈ 0.4
at k_e = 25), and |b − b₀ − εb₁|/ε² tends to 3.34·10³ (bounded: b₁ of (v)). At t = 0 the exact periodic covariance gives Cov(v)·2ε³ − 1
from −9.6·10⁻² to −2.4·10⁻⁵ and ‖E[vY′ᵀ] + bᵀ‖/‖b‖ from 0.15 to 1.2·10⁻³ (O(ε)). ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ falls from 0.35 to 9.4·10⁻⁵ (η =
∞; to 1.7·10⁻⁴ at η = σ₀ and 2.4·10⁻⁴ at η = 0), by 3.6–4 per doubling from k_e = 50 up to 200–400 and by 2–2.9 beyond: the O(ε) of (v), with a
small constant. In the weakest mode the deviation is O(ε²) (at η = ∞: −34%, −17%, −6.0%, −1.7% at k_e = 12.5–100, −6.6·10⁻⁵ at 1600). It
is the sum of two parts of opposite sign: the slow noise's attenuation, ⟨bbᵀ⟩ below Q̄ (Σ with ⟨bbᵀ⟩ in place of Q̄: −56%, …,
−1.1·10⁻⁴), and the rate's correction (Σ with the exact slow generator in place of Ā₀: +49%, …, +4.6·10⁻⁵). The attenuation is 2.4–2.5
times the correction from k_e = 100 at η = ∞ (1.2 times at k_e = 12.5) and 12–23 times at η = 0. Cov(Y) − Cov(Y′) equals the jitter (ε/2)K(0)ᵀK(0) to a relative 0.32, 0.16,
0.081, 0.041 and 0.020 at k_e = 100–1600 (O(ε)), and with the next term −(ε²/2)(d(KᵀK)/dt + 2aKᵀK)(0) to 0.30, 0.079, 0.020, 0.005 and
0.001 (O(ε²)). The floor at t = 0 is 0.763, 0.905, 0.971, 0.992, 0.998, 0.9995, 0.9999 and 1.0000 of the lemma's slow-plus-fast value
(η = ∞). Over the period the floor's average is 0.766, 0.907, 0.972, 0.9925, 0.9981, 0.99951, 0.99987 and 0.99997 of Proposition
D‴.30(ii)'s formula. Pointwise it ranges over [0.9952, 1.0091] of the formula at k_e = 1600, which is the range of tr(H(Σ̄ +
(ε/2)Ĵ(t)ᵀĴ(t))): the fast block's jitter, ±1%. So the formula holds as the period average (v), not pointwise.

**Proposition D‴.27 (a projected adaptive prior restores the leak's a-priori bound).** Let 𝒵 ⊂ ℝ^p be a compact convex set — on the
instance a bounded box in the unconstrained coordinates — and let the prior's atoms of Theorem D‴.23 (or of Proposition D‴.24(ii)) start
in it, Z(0) ∈ 𝒵^{K₀}, and be confined to it: dz_j/ds = Π_{T_𝒵(z_j)}(−η(z_j − T̂_j)), the velocity projected onto the tangent cone of 𝒵 at
z_j in the metric G⁻¹ (for a box and diagonal G: the outward components set to zero on the active faces); at η = ∞, z_j = Π_𝒵(T̂_j), the
G⁻¹-projection of the class barycentre. Let B := sup_{Θ,ϑ}‖G⟨r_Θ∇_ϑΦ̃(·;ϑ)⟩_ρ‖_{G⁻¹} ≤ (max_ρ|f| + a_max)·L_ϑ^G, with
L_ϑ^G := sup_{x,ϑ}‖G^{1/2}∇_ϑΦ̃(x;ϑ)‖ (|f̂| ≤ a_max for the saturated network).

(i) *(the a-priori bound)* Every particle satisfies dist_{G⁻¹}(ϑ_i(s), 𝒵) ≤ B/σ₀ + (dist_{G⁻¹}(ϑ_i(0), 𝒵) − B/σ₀)₊e^{−σ₀s}, at most
max{dist_{G⁻¹}(ϑ_i(0), 𝒵), B/σ₀}, for all s ≥ 0, and the prior stays in 𝒵; 𝒵 being bounded, every solution is bounded, for every m,
every η ∈ (0, ∞] and every target. (For an unbounded 𝒵 the distance bound still holds but bounds nothing: 𝒵 = ℝ^p is Theorem D‴.23's
prior, which runs off.)

(ii) *(the energy)* 𝒥 = L̄ + (σ₀/2)W² — for the nearest-atom prior of Proposition D‴.24(ii) its energy 𝒥_∘, with the cell masses — is
non-increasing along the projected flow, which is the vertex-selection gradient flow of 𝒥 + ι_{𝒵^{K₀}}(Z) in the metric of Theorem
D‴.23(i) (for 𝒥_∘ with the cell masses) — at η = ∞ in its selected form 𝒥(Θ, Π_𝒵T̂(γ)), γ the plan the flow carries, optimal for its own
projected barycentres (Lloyd-consistent) but not necessarily jointly optimal, so that the reduced energy min_{Z∈𝒵^{K₀}}𝒥 itself need not
decrease. At η = ∞ the bounded solutions are expected to converge, by a Kurdyka–Łojasiewicz argument on the finitely many selected pieces,
to a critical point of a Lloyd-consistent piece, which need not be Clarke-critical for the reduced energy; at finite η the solutions are
expected to converge to Clarke-critical points of 𝒥 + ι. At a
critical point with a unique plan the prior is a stationary point of Theorem D‴.23(i) where it is interior, and each boundary atom
satisfies the Karush–Kuhn–Tucker condition z_j = Π_𝒵(T̂_j), T̂_j − z_j in the normal cone of 𝒵 at z_j; at a tie of plans the condition
holds for a convex combination of the active pieces, as in Proposition D‴.17(iii).

(iii) *(nothing changes inside)* As long as the prior stays in the interior of 𝒵 the projected flow is Theorem D‴.23's. In particular, at
an exact representation whose atoms lie in int 𝒵 the stationary point and its linearisation are Theorem D‴.23(i)–(ii)'s, and the bias
removal of (iii) and the local convergence of (v) hold unchanged for starts close enough that the prior stays in int 𝒵.

(iv) *(the price outside)* If the limit has an atom on ∂𝒵 with its class barycentre outside 𝒵 — necessarily so when no stationary point of
Theorem D‴.23(i) has its prior in 𝒵 — the particles of that class leak toward the face, and the limit keeps a bias, bounded by (i), in
place of the unprojected prior's run-off (an atom on ∂𝒵 with T̂_j = z_j is a stationary point of Theorem D‴.23(i) whose prior lies on the
boundary, with no bias). At η = ∞ on a collapsed class the leak is −σ₀(ϑ − Π_𝒵ϑ), the gradient of the exterior quadratic penalty
(σ₀/2)dist²_{G⁻¹}(ϑ, 𝒵) — zero on 𝒵 and growing linearly outside: a switching σ-modification toward 𝒵, the leak off on 𝒵 ([4], [S52]).

*Proof.* (i) z_{π(i)} ∈ 𝒵 at all times, and so is the leak target mΣ_jγ_{ij}z_j of a split row (incommensurate m), a convex combination.
With Π := Π_𝒵 and d := dist_{G⁻¹}(ϑ_i, 𝒵), ½d² is C¹ with gradient ϑ − Πϑ in the metric G⁻¹, so
(d/ds)½d² = ⟨ϑ_i − Πϑ_i, G⟨r∇Φ̃(ϑ_i)⟩⟩_{G⁻¹} − σ₀⟨ϑ_i − Πϑ_i, ϑ_i − z⟩_{G⁻¹}, z the leak target. The projection inequality
⟨ϑ − Πϑ, z − Πϑ⟩_{G⁻¹} ≤ 0 for z ∈ 𝒵 gives ⟨ϑ_i − Πϑ_i, ϑ_i − z⟩_{G⁻¹} ≥ d², so (d/ds)½d² ≤ Bd − σ₀d², i.e. d′ ≤ B − σ₀d where d > 0,
and the comparison principle gives the bound; boundedness follows from that of 𝒵. The bound on B: ‖G⟨r∇Φ̃⟩‖_{G⁻¹} = ‖⟨rG^{1/2}∇Φ̃⟩‖ ≤
sup|r|·L_ϑ^G and |r| ≤ |f| + |f̂|. (ii) The particles' part of d𝒥/ds is Theorem D‴.23(i)'s, minus the squared norm of their velocity; the
prior's is ⟨∇_Z𝒥, Π_T(−c∇_Z𝒥)⟩ = −‖Π_T(−c∇_Z𝒥)‖²/c (c the metric's factor; Moreau's decomposition of −c∇_Z𝒥 into its tangent-cone and
normal-cone parts, which are orthogonal), so 𝒥 is non-increasing. At η = ∞ the prior minimises the W² part over 𝒵 for the current plan:
for each plan γ, min_Z 𝒥 = L̄ + (σ₀/2)Σ_j[Σ_iγ_{ij}‖ϑ_i − T̂_j‖²_{G⁻¹} + q_j dist²_{G⁻¹}(T̂_j, 𝒵)] (minimising
Σ_iγ_{ij}‖ϑ_i − z‖² = q_j‖z − T̂_j‖² + const over z ∈ 𝒵), each piece C^{1,1} and definable (the squared distance to a box is piecewise
quadratic), and its gradient in ϑ_i is the learning term and the leak toward Π_𝒵T̂_j (the derivative of q_j dist²(T̂_j, 𝒵) adds
γ_{ij}(T̂_j − Π_𝒵T̂_j) to the within-class term). Between plan refreshes the particles follow the gradient flow of the piece E_γ(Θ) :=
min_Z𝒥(Θ, Z; γ) of the plan γ they carry, and at a refresh the new plan is optimal for the current prior and the prior then re-minimises,
so the selected energy does not increase there either. But γ is Lloyd-consistent — optimal for ν_Z with Z = Π_𝒵T̂(γ) — and need not
attain min_γE_γ = min_Z𝒥: the flow is then not a subgradient trajectory of the reduced energy, and Proposition D‴.17's argument, whose
fixed prior makes every selected piece the active one, does not apply as it stands (Part H(e): the reduced energy rises by 26% while the
selected energy falls). A Kurdyka–Łojasiewicz argument on the finitely many definable pieces that share the limiting value should give
the bounded solutions finite length, with a critical point of a Lloyd-consistent piece as the limit; it is not carried out here. At finite η the joint energy carries the indicator of 𝒵^{K₀}, and the Kurdyka–Łojasiewicz argument
with it is not carried out here. At a critical point the projected velocity of every atom vanishes, which is the stated condition.
(iii) The projection is inactive in int 𝒵, and Theorem D‴.23(v)'s ball keeps the prior within ρ(η/(σ₀q_k))^{1/2} of the atoms.
(iv) If no stationary point of Theorem D‴.23(i) has its prior in int 𝒵, a critical point of (ii) has an atom on ∂𝒵. At η = ∞ with a class
collapsed at ϑ, T̂ = ϑ and z = Π_𝒵ϑ. ∎

*Status.* (i) proved; (ii) the monotonicity of the selected energy and the characterisation of the critical points proved, the convergence
argued at η = ∞ (a Kurdyka–Łojasiewicz argument on the selected pieces, not carried out) and formal at finite η; (iii) proved, its local convergence as Theorem D‴.23(v) (sketched); (iv) the
identification at η = ∞ on collapsed classes proved, the boundary atoms and the bias outside numerical. The bound of (i) is crude (it uses
the worst-case force; L_ϑ^G ≤ 4.658 is the closed-form bound, 4.044 by search), and its content is qualitative: boundedness is no longer a
hypothesis. Parameter projection onto a convex set, the exterior quadratic penalty and the switching σ-modification are classical
robustness modifications of adaptive laws, and the switching σ-modification keeps the ideal (unbiased) properties of the adaptive law
when the true parameter lies in its set ([4], [S52]); on collapsed classes at η = ∞ the construction is that law. New is its transfer to
the learned centre of a mass-constrained K₀-atom prior at finite η — a projected adaptive law for the prior's atoms, whose velocity jumps
at the faces while the particles' law keeps its smooth gradient form — and its compatibility with Theorem D‴.23's bias removal wherever
the target is representable inside 𝒵. Checked: `tests/test_adaptive_prior.py` (beyond the saturation, the projected prior at η = ∞
(LSODA) keeps the particles within 2 of the box and converges, while without the box the prior's coordinate passes 11 by s = 10⁴; at
η = σ₀ the projected midpoint rule at dt = 0.125 keeps the prior in the box, brings its amplitude atom onto the face, keeps the particles
bounded and 𝒥 non-increasing; at the misplaced prior near the exact representation the projected and unprojected runs are identical).

*Numerics* (`exp_adaptive_prior.py` Part H → `results/chaos/adaptive_prior.md`). The box 𝒵 holds the amplitudes and centres to 0.9 of
a_max and c_max (α ∈ [−5.889, 5.889], γ ∈ [−4.417, 4.417]) and the width coordinate in [−6, 6]; the instance's exact atoms are inside.
(a) Part G(b)'s saturated target (a weighted amplitude 4.8 > a_max), the K₀-atom system from the exact atoms with the prior there: under
the projected prior both η = σ₀ (the projected midpoint rule — positions projected inside the step — at dt = 0.25 and 0.125, to s = 10⁵)
and η = ∞ (LSODA, rtol 10⁻¹⁰) converge to the same bounded stationary point (largest drift 1.8·10⁻¹⁵, 3.5·10⁻¹⁵ and 1.2·10⁻¹⁶; the two
steps' end states within 3.4·10⁻¹²), the first prior atom's amplitude coordinate on the face of 𝒵 (z_a = 5.889) and the particle's at
6.728 (amplitude 3.7325), with the fit error 0.1415 — the two steps' agreement shows their shared fixed point; the scheme's
consistency: at s = 300 the trajectory lies within 6.3·10⁻⁷, 1.9·10⁻⁷ and 3.8·10⁻⁸ of a reference at dt = 1/64 for dt = 0.25, 0.125,
0.0625. The particles' largest G⁻¹-distance from 𝒵 along the run is 0.92 (η = σ₀) and 0.89
(η = ∞), against the bound B/σ₀ ≤ 3686 (B ≤ 26.08 with the closed-form L_ϑ^G). (b) Part C's K₀-atom loop at the misplaced prior (η = σ₀,
6000 slow units) is identical with and without 𝒵 (largest difference 0.0): the box is never approached — the prior and the particles keep
at least 16.9% of the box's width to every face — so the two runs coincide by construction. (c) The box trades the fit error against the
parameters' bound: with the amplitude bound 0.75, 0.8, 0.9, 0.95, 0.99 of a_max — boxes that contain the instance's exact atoms — the
limit's fit error is 0.1676, 0.1611, 0.1415, 0.1268, 0.1101, decreasing between the fixed prior's 0.1750 and the unprojected adaptive
prior's 0.1050 at s = 10⁶ (which that prior approaches only as its amplitude coordinate runs off), the prior's amplitude atom on the face in
every case; the bound 0.7 (0.1736) already excludes the exact atoms (the start is projected at once), and the tighter boxes 0.65, 0.6,
0.5 are worse than the fixed prior (0.1792, 0.1845, 0.1949): a box should contain the representation. (e) At η = ∞ on six particles in
a Lloyd-consistent, non-optimal pairing (a seeded search; zero amplitudes, a small three-bump target, the box) the selected energy falls
from 4.52·10⁻³ to 4.95·10⁻⁵ (its largest increment −3.4·10⁻⁶) while the reduced energy min_Z𝒥 rises from 1.251·10⁻³ to 1.574·10⁻³
before the plan becomes optimal at s = 80. (d) On Part G(a)'s
four-bump target, where without the box the merged atom's width coordinate grows like ln s and the solution leaves every compact set, the
box stops it: from the exact split fit at η = ∞ (the joint flow for 40000 slow units, then the collapsed dynamics) κ passes 2.159, 3.843,
5.529 at s = 10⁶, 10⁷, 10⁸ as before and stays at the face κ = 6 from s = 10⁹ (the atom within 6.6·10⁻⁷ of it, drift 5.3·10⁻¹⁷), with
the fit error 1.943376·10⁻², within 9.3·10⁻⁷ of the minimiser's on the face s = s_max that the unprojected solution approaches only as
s → ∞.

**Proposition D‴.28 (the adaptive prior's floor removal is a gauge property: it survives an approximation error).** Let ϑ^b(q) be a
nondegenerate critical point of the unregularised weighted loss L_q(Θ) := ½⟨|f − Σ_kq_kΦ̃(·;ϑ_k)|²⟩_ρ — not necessarily an exact
representation: r_b := f − f̂_b may be non-zero — and let the masses enter the network only through the products w_ka(α_k), as they do for the
saturated network. The amplitude–mass gauge G_w: α_k ↦ a⁻¹(q_ka(α_k)/w_k), the other coordinates unchanged, is defined on
W_g := {w : q_k|a^b_k| < w_ka_max for every k} and satisfies f̂_w∘G_w = f̂_q, hence L_w∘G_w = L_q.

(i) *(the stationary points carry no function-level floor)* For every w ∈ W_g, ϑ^b(w) := G_wϑ^b(q) is a nondegenerate critical point of
L_w with the same function f̂_b, the same residual r_b and the same Hessian signature. With the prior at its atoms it is a stationary
replicated configuration of the nearest-atom leak toward an adaptive prior (Proposition D‴.24(iii)), so two systems at their samples'
gauge transforms carry the same function whatever the approximation error, with the fit error ‖r_b‖_ρ. The function-level floor — the
limit of E‖f̂_m − f̂_M‖² along the dynamics — therefore vanishes at the samples of W_g with positive margins in (iii) when ϑ^b(q) is a
strict local minimum and the systems converge to the gauge transforms (numerical); W_g shrinks as the approximation error grows (its
thresholds q_k|a^b_k|/a_max rise with the amplitudes). (At an exact representation W_g is Proposition D‴.24's W.)

(ii) *(the fixed prior's floor)* At the matching prior ϑ⁰ = ϑ^b(q) the function-level weight-sensitivity of the fixed prior's stationary
point is

    F = σ₀Ĵ(Ã + σ₀I)⁻¹g̃,    Ã := H_L − R̃,

Ã the whitened weighted Hessian of L_q at ϑ^b(q) (R̃ the whitened weighted residual curvature) and g̃ the whitened gauge displacement, the
direct effect of the weights on f̂ being Ĵg̃: Proposition D‴.24(i) with the full Hessian in place of the Gauss–Newton block. Along the same
derivation the adaptive prior's sensitivity is F(σ) at σ = 0, which vanishes, as (i) says.

(iii) *(the splitting margins under the adaptive prior)* For a scalar amplitude per kernel (one output), at a critical point the residual
curvature R_k = ⟨r∇²_ϑΦ̃(·;ϑ_k)⟩_ρ has a zero amplitude row and column, and its block in the centre and width coordinates is a(α_k) times
a function of those coordinates alone. Hence
R_k(ϑ^b(w)) = (q_k/w_k)R_k(ϑ^b(q)) exactly, and Proposition D‴.5's splitting margins at the adaptive prior's replicated stationary point are

    σ₀ − (q_k/w_k)λ⁺_k,    λ⁺_k := λ_max(G^{1/2}R_k(ϑ^b(q))G^{1/2}) ≥ 0:

a class sampled below q_kλ⁺_k/σ₀ can split. At an exact representation λ⁺ = 0 and every margin is σ₀ (Proposition D‴.24(iii)); with an
approximation error, and ϑ^b(q) a strict local minimum, the adaptive prior's admissible set is {w ∈ W_g : w_k > q_kλ⁺_k/σ₀}. For vector
amplitudes (n ≥ 2 outputs, a(α) ∈ ℝⁿ saturated componentwise) only the centre–width block scales by q_k/w_k: stationarity in the
amplitudes gives Ψ_r = 0 for every output r but stationarity in the centres only Σ_ra_r∂_γΨ_r = 0, so the amplitude–centre block
a′_r(α_r)∂_γΨ_r survives, is gauge-invariant in the physical amplitudes, and no closed margin formula follows.

*Proof.* (i) G_w is a diffeomorphism from {α : |q_ka(α)/w_k| < a_max} onto its image with f̂_w(G_wΘ) = f̂_q(Θ), so L_w∘G_w = L_q, and it maps
critical points to critical points with the same Hessian signature; the function and the residual are unchanged by construction. The
leak toward an adaptive prior at its atoms vanishes (Theorem D‴.23(i), Proposition D‴.24(ii)). (ii) Let u := ∂_w(G_wϑ^b(q)) at w = q, D :=
diag(q_kG⁻¹), A := ∇²L_q(ϑ^b), B := ∂_w∇_ΘL_w at fixed Θ, J := ∂_Θf̂, Φ := ∂_wf̂ at fixed Θ. Differentiating f̂_w(G_wΘ) = f̂_q(Θ) gives
Φ = −Ju, and differentiating ∇L_w(G_wϑ^b) = 0 gives B = −Au. The fixed prior's stationarity ∇L_w(Θ) + σ₀D(Θ − ϑ⁰) = 0, linearised at
(q, ϑ^b), reads (A + σ₀D)δΘ = −Bδw = Auδw, so δf̂ = Φδw + JδΘ = −J(A + σ₀D)⁻¹[(A + σ₀D) − A]uδw = −σ₀J(A + σ₀D)⁻¹Duδw; in the whitened
coordinates (Ĵ = JD^{−1/2}, Ã = D^{−1/2}AD^{−1/2}, g̃ := −D^{1/2}u) this is σ₀Ĵ(Ã + σ₀I)⁻¹g̃. With σ₀ = 0, δΘ = uδw and δf̂ = 0. At an exact
representation R̃ = 0, Ã = H_L and g̃ = ∂_wx⁰. (iii) With a scalar amplitude and Ψ(γ, κ) := ⟨rφ(·; c(γ), s(κ))⟩_ρ, ⟨rΦ̃(·;ϑ)⟩_ρ = a(α)Ψ, whose gradient is
(a′Ψ, aΨ_γ, aΨ_κ) and whose Hessian is [[a″Ψ, a′Ψ_γ, a′Ψ_κ], [a′Ψ_γ, aΨ_γγ, aΨ_γκ], [a′Ψ_κ, aΨ_γκ, aΨ_κκ]]. At a critical point the gradient
vanishes, so Ψ = Ψ_γ = Ψ_κ = 0 (a′ > 0, a ≠ 0) and only the block aΨ_{(γ,κ)} remains; r and Ψ are unchanged by the gauge and a scales by
q_k/w_k. The α-direction contributes the eigenvalue 0, so λ⁺_k ≥ 0. ∎

*Status.* (i) and (iii) proved, exactly, as statements about the stationary points and their linearisation (no linearisation in w: the
gauge maps critical points to critical points); (ii) as Proposition D‴.24(i) — exact at q, first order in w; the convergence of the
dynamics to the gauge transform numerical, and with it the floor's vanishing. The amplitude–mass degeneracy is the homogeneity of the
conic view [49], [55], and the rescaling symmetries of networks, which map critical points to critical points with the same function
while changing their curvature, are classical ([S55]); what is new is that this symmetry, not exact representability, is the mechanism of
Proposition D‴.24's floor removal — so at a strict local minimum in the interior the adaptive prior keeps the fitted function independent
of the sampled weights without an exact representation — and, for a scalar amplitude per kernel, the exact scaling of the residual
curvature under the gauge, which gives the adaptive prior's splitting margins in closed form. Its hypotheses are a strict local minimum in
the interior and samples in W_g with positive margins: on the instance P(w^m ∈ W_g) is only 0.52 at m = 15 (0.87 at m = 60), and the branch
ends at the saturation near ε = 0.68; targets whose best approximant lies on the saturation boundary (a sharp friction step fitted by a
dipole of two kernels; the cogging servo of Proposition D‴.29, globally; Proposition D‴.27's beyond-the-saturation example) are outside
it. Checked: `tests/test_gauge_floor.py` (at ε = 0.3, a critical point with a 1.7% fit error: the gauge transform a critical point at two
sampled weight vectors, the critical point continued in w from ϑ^b(q) landing on it; the adaptive prior's sensitivity, from critical
points reached by Newton from ϑ^b(q), below 10⁻⁸ of the fixed prior's, the fixed prior's equal to σ₀Ĵ(Ã + σ₀I)⁻¹g̃ to 10⁻⁶ and not to
the Gauss–Newton form; the residual curvature's zero amplitude row, its scaling (q_k/w_k), the margins σ₀ − (q_k/w_k)λ⁺_k and the
dynamics' separation eigenvalues −σ₀ + (q_k/w_k)eig(G^{1/2}R_kG^{1/2}) to 10⁻⁹; adaptive-prior runs at η = ∞ and σ₀ converging to the gauge
transform, 1.2 away from ϑ^b(q), with the approximant's function).

*Numerics* (`exp_gauge_floor.py` → `results/chaos/gauge_floor.md`). The three-bump instance with a fourth bump ε·exp(−(x − 0.8)²/(2·0.3²))
that K₀ = 3 kernels cannot represent, ε = 0–0.5, the masses q = 1/3, σ₀ = 0.15σ_crit; ϑ^b(q) the weighted approximant, a critical
point of the unregularised loss continued in ε from the exact representation (relative fit error 3.6·10⁻³, 6.8·10⁻³, 1.2·10⁻², 1.7·10⁻², 2.5·10⁻² at ε = 0.05, 0.1, 0.2,
0.3, 0.5; gradient residual ≤ 1.6·10⁻¹⁵; nondegenerate, the smallest eigenvalue of Ã 4.6·10⁻⁵–3.1·10⁻⁴). (i) At six weight vectors
0.01 inside W_g the unregularised gradient at G_wϑ^b(q) is ≤ 1.6·10⁻¹⁵, and Newton's method continued in w from ϑ^b(q) (20 steps) lands
on G_wϑ^b(q) to ≤ 1.7·10⁻¹¹ with f̂ unchanged to ≤ 5.5·10⁻¹³ — plain Newton from ϑ^b(q) at w does so in 3–5 of the 6 cases, elsewhere
reaching other critical points of L_w; the adaptive prior's function-level sensitivity, from the critical points reached from ϑ^b(q),
is 1.4·10⁻¹¹–2.7·10⁻¹⁰ of the fixed prior's. The samples fall in W_g with the probability 0.52, 0.87–0.88, 0.997 at m = 15, 60, 240
(ε = 0.05–0.3; 0.26, 0.71, 0.98 at ε = 0.5), and the branch continued to ε = 0.68 reaches max_k|a^b_k|/a_max = 0.999, the smallest
eigenvalue of Ã 1.2·10⁻⁸ and W_g's threshold 0.333 = q₂ for its class: there w = q leaves W_g. (ii) The fixed matching prior's
sensitivity equals σ₀Ĵ(Ã + σ₀I)⁻¹g̃ to ≤ 2.2·10⁻⁹ and misses Proposition D‴.24's Gauss–Newton form by 7.0·10⁻⁴–5.6·10⁻³, growing with ε as
‖R̃‖/‖H_L‖ does (5.3·10⁻⁴–2.8·10⁻³); its floor law per unit of 1/m − 1/M rises from 6.17·10⁻² (ε = 0) to 8.52·10⁻² (ε = 0.5). (iii) The
residual curvature's amplitude row is ≤ 7.4·10⁻¹⁴ (its largest entry 2.6·10⁻⁴–2.9·10⁻³), λ⁺_k ≤ 7.2·10⁻⁴ (at most 0.10σ₀), and over every
count vector of m = 15 and 60 in W_g the adaptive prior's splitting margins equal σ₀ − (q_k/w_k)λ⁺_k to ≤ 2.1·10⁻¹⁴, none negative (the
smallest +6.23·10⁻³); the separation eigenvalues of the dynamics at a replicated configuration (counts (6, 4, 5), the Jacobian
of the averaged drift at η = σ₀ and ∞) equal −σ₀ + (q_k/w_k)eig(G^{1/2}R_kG^{1/2}) to ≤ 5.8·10⁻¹². *The dynamics* (Part B: warm clouds of width 0.02 around ϑ^b(q), m = 15, 60, 240 against M = 510, four seeds, 6000
slow units, ε = 0.1 and 0.3, the nearest-atom leak). Under the fixed prior ‖f̂_m − f̂_M‖² stays within 8% of its value at s = 600 (the ratio at s = 6000 0.917–1.003; at ε = 0.1:
3.1·10⁻⁴–5.9·10⁻³, 1.1·10⁻⁴–1.0·10⁻³ and 5.0·10⁻⁶–8.2·10⁻⁵ at m = 15, 60, 240) and the fit error above the approximant's
(9.0·10⁻³–1.2·10⁻¹ against 6.84·10⁻³), the classes splitting at the samples where the fixed prior's point is a saddle (diameters up to
1.6). Under the adaptive priors the classes stay replicated, their means approach the gauge transform G_wϑ^b(q), ‖f̂_m − f̂_M‖² falls (at
m = 60, ε = 0.1, from 6.9·10⁻⁶–1.2·10⁻⁴ at s = 600 to 1.0·10⁻⁹–1.4·10⁻⁸ at η = σ₀ and to 3.2·10⁻¹¹–3.8·10⁻⁹ at η = ∞ by s = 6000), and
the fit error tends to the approximant's: at the samples at least 0.02 inside W_g it exceeds it at s = 6000 by a median 1.2·10⁻⁴,
1.2·10⁻⁴, 1.6·10⁻⁵ (relative; at most 1.1·10⁻¹, 1.2·10⁻², 1.3·10⁻³) at m = 15, 60, 240. At the m = 15 sample (0.533, 0.267, 0.200),
outside W_g (margins −0.0016 and −0.0097 at ε = 0.1 and 0.3), ϑ^b(q)'s branch has no gauge transform: the adaptive priors' fit errors at
s = 6000 (2.2·10⁻² and 1.6·10⁻² at ε = 0.1, 3.7·10⁻² and 3.3·10⁻² at ε = 0.3, for η = σ₀ and ∞) are still falling toward the minimum of
L_w over the closed amplitude box on the branch, 1.107 and 1.586 times the approximant's fit error (one amplitude at the bound, its
coordinate running off), and a gauge transform with two classes exchanged exists there, which the warm-cloud dynamics does not reach; at (0.400, 0.200, 0.400), 6·10⁻⁴ inside W_g at ε = 0.3, the gauge transform's
amplitude lies near the saturation and the approach is slow (2.8–3.1 from it at s = 6000), as at the boundary of Proposition D‴.24's W.

**Proposition D‴.29 (second-order plants: the sliding-variable reduction).** Let the plant be ẍ = f(x, ẋ) + u, x_m a T₀-periodic C³
reference, e := x − x_m, s := ė + λe (λ > 0), and let the controller and the adaptive law be

    u = ẍ_m − λė − k_s s − f̂_Θ(z),    ϑ̇_i = G∇_ϑ(sΦ̃(z; ϑ_i)) − (σ₀/k_s)(ϑ_i − T_i),

z the coordinates of (x, ẋ) on which f depends and T_i any of the leak targets of §5d (a fixed prior, the nearest-atom or the OT plan, an
adaptive or a projected prior).

(i) *(the error channel)* ṡ = −k_s s + r_Θ(z) and ė = −λe + s, with z = z_m + (e, ė) the regressor at the plant's state: Theorem J's error
channel with s in place of e.

(ii) *(the averaged slow loop)* As k_s → ∞ at fixed λ, s and e are slaved — s = (d/dt + k_s)⁻¹[r_Θ(z_m)] + O(k_s⁻²) = r_Θ(z_m)/k_s +
O(k_s⁻²) and e = (d/dt + λ)⁻¹s, the periodic solutions of the two filters — and in the slow time τ = t/k_s the parameters follow Theorem
J's averaged slow loop along the regressor's orbit z_m(t),

    dϑ_i/dτ = G⟨r_Θ∇_ϑΦ̃(·; ϑ_i)⟩_ρ − σ₀(ϑ_i − T_i),    ρ the occupation measure of z_m over one period,

up to O(1/k_s), the prior's law at finite η being dz_j/dt = −(η/k_s)(z_j − T̂_j) in real time. Both forms of the slaved error — the exact
filter (d/dt + k_s)⁻¹ and r/k_s — are accurate to O(k_s⁻²) in s, i.e. to O(1/k_s) relatively, the parameters' within-period motion and
the nonlinearity r(z) − r(z_m) entering at the same order; the exact filter keeps the error filter's attenuation of the residual's
harmonics, which r/k_s drops at the relative order n²ω₀²/(2k_s²) in the rms of the n-th harmonic, and is the better prediction where that
term dominates (a small residual of high frequency).

(iii) *(the transfer)* Every statement of §5d about the averaged slow loop holds for the second-order plant with this ρ — Theorems D‴, D‴.21
and D‴.23, Propositions D‴.5, D‴.9, D‴.12, D‴.17, D‴.24(i)–(v), D‴.27 and D‴.28 — and the tracking error of the position is the filtered
residual of (ii); the real-loop statements (Proposition D‴.24(vi), Corollaries D‴.14 and D‴.26) would need the variational equation of (i).

(iv) *(where the kernels live)* When f depends on one coordinate (a cogging torque on the position, a friction law on the velocity), kernels
on that coordinate see the orbit as an interval, and the instance is §5d's one-dimensional one with another target. Isotropic Gaussian
kernels on the phase plane over a circular orbit (x_m = R sin ω₀t with ω₀ = 1 in the units where the kernel is isotropic in (x, ẋ)) are
over-parametrised: on the circle z = R(sin θ, cos θ) a kernel restricts to a·e^{−(R² + |c|²)/(2s²)}·e^{(R|c|/s²)cos(θ − θ_c)}, a von Mises
function of three effective parameters for four, so in general position the Gauss–Newton block has a zero eigenvalue for every kernel
(more at degenerate atoms, as for kernels centred on the orbit's centre) — a direction that the fixed prior's leak pins at the rate σ₀,
that an adaptive prior leaves neutral (λ^− = 0, Theorem D‴.23(ii)), and that a projected prior keeps bounded (Proposition D‴.27). On
other orbits the same directions are weakly excited rather than unexcited: the classical partial excitation of Gaussian networks along an
orbit ([6]–[8]).

*Proof.* (i) ë = f + u − ẍ_m = r_Θ(z) − λė − k_s s, so ṡ = ë + λė = −k_s s + r_Θ(z), and ė = s − λe by the definition of s. (ii) The error
subsystem (s, e) is linear, exponentially stable (eigenvalues −k_s and −λ) and driven by r_Θ(z), while the parameters move at the rate
O(1/k_s) (their velocity is G∇_ϑ(sΦ̃) with s = O(1/k_s)). Theorem J's argument — for the smooth fixed-prior leak near a hyperbolic
critical point, using the error channel's exponential stability, the regressor's smoothness and its Step 3's normal form — carries over
with the fast pair (s, e) once its error-channel steps are redone for it (s slaved at k_s, e filtered at the fixed λ: the error channel's
multipliers are e^{−k_sT₀} and e^{−λT₀}), the separation being between λ and the parameters' rate; for switching targets (plan and cell
switches, adaptive and projected priors) not even the first-order loop has more than Theorem D‴.18's locally constant case. At the slaved state z = z_m + O(1/k_s), and the period average of G∇_ϑ(sΦ̃(z_m; ϑ_i)) with s = r_Θ(z_m)/k_s is (1/k_s) times the averaged
learning term. For the n-th harmonic of r, (d/dt + k_s)⁻¹ has the gain (k_s² + n²ω₀²)^{−1/2} = (1/k_s)(1 − n²ω₀²/(2k_s²) + …). (iii) The
statements of §5d use the averaged loop only through ρ, the target on its support and the network. (iv) |z − c|² = R² + |c|² −
2R|c|cos(θ − θ_c) on the circle; the map (a, |c|, θ_c, s) ↦ (a e^{−(R²+|c|²)/(2s²)}, R|c|/s², θ_c) has a one-dimensional fibre, along which the
kernel's restriction, and hence f̂ on ρ's support, is constant. ∎

*Status.* (i) and (iv) exact; (ii) and (iii) formal — first-order averaging with the fast pair (s, e), Theorem J's steps carrying over for
the smooth fixed-prior leak near a hyperbolic critical point once the error-channel steps are redone, and switching targets beyond even the
first-order theory; (iii) covers the statements about the averaged loop, not the real-loop ones (Proposition D‴.24(vi) and Corollaries
D‴.14 and D‴.26 would need the variational equation of (i)). The filtered error is the classical construction of adaptive control for higher-order plants, used with Gaussian networks by
Sanner and Slotine ([1]); what is added is the transfer of the transport leaks, their floors and the adaptive and projected priors to such
plants, and the exact fibre of isotropic phase-plane kernels on a circular orbit (the weak excitation of the same directions on other
orbits being the classical partial excitation [6]–[8]). Checked: `tests/test_sliding.py` (the phase-plane loop as the gradient flow of L̄;
on a circular orbit one zero Gauss–Newton eigenvalue per kernel in general position, the next above 10⁻¹¹ of the largest; the servo's parameter displacement over ten periods within 3% of the
averaged loop's at k_s = 25 and closer at k_s = 50, and its position error within 1% of both predictions of (ii) at k_s = 50).

*Numerics* (`exp_servo_cogging.py` → `results/servo/servo_cogging.md`). The servo ẍ = −τ_cog(x) + u with the cogging torque τ_cog(x) = 0.4
sin(3x + 0.3) + 0.1 sin 7x, x_m = sin t, λ = 5; K₀ = 4 saturated Gaussian kernels on the position (a_max = 4, c_max = 1.6, widths in [0.05,
1]), masses 1/4, σ₀ = 0.01. (A) ϑ^b(q), the nondegenerate local minimum of the unregularised weighted loss that BFGS reaches from a fixed
start, has the relative fit error 1.46·10⁻², and Proposition D‴.28 is its local analysis: the whitened unregularised Hessian Ã spans
2.9·10⁻⁴σ₀–30σ₀ (four modes below 0.1σ₀: the adaptive prior's slowest rates 1.4·10⁻⁴σ₀ at η = σ₀ and 2.9·10⁻⁴σ₀ at η = ∞, against the fixed
prior's ≈ σ₀), the residual curvature's λ⁺_k ≤ 0.0098σ₀ (the adaptive prior's margins positive on all of W_g = {w > (0.061, 0.085, 0.119,
0.008)}), the fixed matching prior's floor law 1.29·10⁻² per unit of 1/m − 1/M and its margin +9.9·10⁻³ at w = q. Globally the loss has its
infimum at the saturation: of 120 BFGS runs from random starts one ends at ϑ^b(q) and 33 below its fit error, all with saturated amplitudes;
the best reaches 6.39·10⁻³, and with its two saturated amplitudes held at |α| = 6, 8, 10, 14 (|a| = 3.62–3.99) and the rest re-optimised the
fit error is 7.14·10⁻³, 6.65·10⁻³, 6.49·10⁻³, 6.40·10⁻³ — approached, not attained. Globally the servo is Proposition D‴.27's case (critical
points at infinity), and the runs below are statements about the basin of ϑ^b(q) and about an engineer's model. A misplaced prior — the
approximant of the cogging model with the phase 0.5 in place of 0.3, the engineer's model, its own fit error 0.200 — has a replicated
stationary point (the same from its atoms and from ϑ^b(q)) with the fit error 5.45·10⁻² and the splitting margin −3.97·10⁻³: a D‴.5 saddle.
(B) The averaged loop (m = 16, 64, 128 against M = 256, three seeds, 2000 slow units). From warm clouds around ϑ^b(q) — computed from the
true torque — with the matching prior ϑ^b(q): the fixed prior's nearest-atom leak keeps a finite-size floor (‖f̂_m − f̂_M‖² 1.85·10⁻³,
2.24·10⁻⁴, 4.49·10⁻⁵; fit error 0.174, 0.066, 0.038 against the approximant's 0.0146); the adaptive prior's floor is 38–780 times lower
(η = ∞: 4.85·10⁻⁵, 2.88·10⁻⁷, 6.36·10⁻⁸), with the fit error 0.0298 (the classes not yet replicated at m = 16 by s = 2000), 0.0148, 0.0147;
the OT leak toward the matching prior sits at the approximant (0.01462, no floor); the σ-modification toward ϑ⁰ = 0 at the same σ₀ keeps a
large bias (0.328) — and no leak at all (σ₀ = 0) does as well as the best of them (0.0180, 0.0147, 0.0146; floors 7.7·10⁻⁶, 2.4·10⁻⁷,
6.6·10⁻⁸). In this noise-free simulation the floors and biases that the adaptive and OT leaks remove are those a fixed leak creates; the
leaks' purpose — robustness to noise and disturbances, Propositions D‴.27, D‴.30, D‴.31 — is not exercised here. With the misplaced prior
from that start, at η = ∞ the prior never acts (its initial position enters only through the first cells: the runs are the matching ones),
at η = σ₀ its bias decays (0.035, 0.020, 0.019 at s = 2000), and the fixed designs keep a bias — their m-systems leave the saddle by
splitting their classes (diameters 0.32–1.07) and end at 0.048 (OT) and 0.10, 0.060, 0.048 (the nearest atom; above the replicated point's
0.0545 at m = 16 and 64). From warm clouds around the engineer's model, with the prior at it — the start an engineer has — the fixed designs
stay biased (0.102, 0.060, 0.048; OT 0.048), the adaptive prior removes most of the bias but slowly (η = ∞: 0.027, 0.025, 0.022 at s = 2000;
η = σ₀: 0.032, 0.028, 0.028), and no leak reaches the basin's fit (0.0152, 0.0147, 0.0147): its m particles spread and explore directions
that the K₀-atom structure the leaks impose leaves weakly excited. The collapsed K₀-atom flows give the time scale: from the engineer's
model at η = ∞ — where the collapsed flow is the unregularised one — the fit error is 0.0297, 0.0213, 0.0169, 0.0155 at s = 600, 2000, 10⁴,
10⁵, still 0.55 from ϑ^b(q) in the parameters (10⁵ slow units are 2.5·10⁶ time units of the real servo at k_s = 25), at η = σ₀ 0.0345,
0.0259, 0.0181, 0.0160; from ϑ^b(q) toward the misplaced prior at η = σ₀ 0.0240, 0.0183, 0.0154, 0.0149. (C) The real servo (k_s = 25, RK4
at dt = 2·10⁻³, 1000 time units, m = 16 and 64, seed 0) from the oracle-start designs' end states, with each design's leak (its targets
refreshed, and the error recorded, every 0.02 time units; the finite-η prior relaxed in real time): the position error's rms is 1.7–2.1·10⁻⁵
for the adaptive designs, the matching OT design and no leak, 7.8·10⁻⁵ for the misplaced OT design, 9.1·10⁻⁵–1.8·10⁻⁴ for the fixed
nearest-atom prior and 5.4–5.5·10⁻⁴ for the σ-modification; it equals the prediction of (ii) with the exact filter to 1.0001–1.0029 at
t = 500 and 1000 in every design (0.990–1.005 over the first five periods, where the parameters still move), and with the slaved s = r/k_s
to 0.931–0.996, the deficit largest for the adaptive designs, whose residual is the approximation error's broad fine scale; over these 40
slow units no design drifts, no leak included. (D) Phase-plane kernels (five isotropic kernels in general position, σ₀ = 0.01): on the
circle 5 exact zero Gauss–Newton eigenvalues (one per kernel) and 8 below 10⁻⁶σ₀, the (K + 1)-th eigenvalue 2.4·10⁻¹⁰ of the largest; on
ellipses with the velocity axis scaled by 1.1 and 2, 8 and 6 below 10⁻⁶σ₀ and none exact; on the orbit of 1.2 sin t + 0.3 sin 2t 3; with the
five centres at the origin 16 exact zeros of 20.

**Proposition D‴.30 (measurement noise at perfect tracking: a coherent, high-passed noise, and the errors in variables).** Let (CL) be the
real loop of Corollary D‴.26 (n = 1) at the lift Γ_w of an exact representation of representable weights w, with the prior at its atoms —
fixed there (η = 0) or adaptive (0 < η ≤ ∞) — and ∂_xf̂ = ∂_xf, ∂²_xf̂ = ∂²_xf on the orbit (an exact representation, an identity of
functions, supplies both), and let the controller and the adaptive law use the measurement
y = x + ν in place of x: u = ẋ_m − k_e(y − x_m) − f̂_Θ(y), ϑ̇_i = G∇_ϑ((y − x_m)ᵀΦ̃(y; ϑ_i)) − (σ₀/k_e)(ϑ_i − z_{π(i)}), ν a stationary
zero-mean noise, independent of the rest, with the spectral density N₀ at the frequencies of the regressor and of the error channel. Let
H := H_L(w) = ⟨ĴᵀĴ⟩ (the whitened weighted Gauss–Newton matrix, Corollary D‴.26's Ĵ) with the eigenpairs (h_j, v_j), a(t) := ∂_xf̂(x_m(t)),
Q := ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩ and Q_jj := v_jᵀQv_j.

(i) *(coherence)* The particles of a class see the same y, the same error and the same prior atom, so a replicated configuration stays
replicated under every path of ν: on the replicated subspace the m-system is the weighted K₀-atom loop driven by the same noise, for every
m with these weights, at every order in ν. Measurement noise does not average out with m, as the per-particle dither of Theorem J′ does.
This is invariance, not attraction. A separation u = G^{−1/2}δϑ_i of class k obeys u̇ = (G^{1/2}ε_m∇²_ϑΦ̃(y; ϑ_k)G^{1/2} − σ₀/k_e)u
along the noisy class trajectory, and to second order in ν its mean drift is E[ε_m²]G^{1/2}∂_x∇²_ϑΦ̃G^{1/2} + E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}
(the measured regressor; the parameters' cross term negligible), E[e] the mean tracking error. At the lift — the parameters at their
atoms, E[e] the noise's own mean error at fixed parameters — its period average S_k gives the decay rate σ₀/k_e − Eν²λ_max(S_k) to
leading order, so near the lift the replicated configuration attracts — and an m-system started near it follows the K₀-atom loop — only
below the noise level Eν² = (σ₀/k_e)/max_kλ_max(S_k). Along the noisy trajectory E[e] also carries the mean tracking error −Ĵȳ(t)/k_e of
the developing errors-in-variables deviation ȳ(t) of (iii), of the same order in Eν², so the drift S_k(t) changes on ȳ's time scale and
to this order the separation follows u̇ = (Eν²S_k(t) − σ₀/k_e)u (the first-order multiplicative noise, high-passed by the loop, adds a
negligible cumulant to its exponent); the long-time threshold — where ȳ is O(1) at these noise levels and the expansion no longer
applies — is open.

(ii) *(the floor at first order)* To first order in ν the measured error's noise is ((d/dt − a)/(d/dt + k_e))ν — the loop high-passes it —
and in the slow time the collective coordinates diffuse with the matrix (N₀/k_e)Q; their stationary covariance solves
Ā₀Σ + ΣĀ₀ᵀ + (N₀/k_e)diag(Q, 0) = 0 (Ā₀ Corollary D‴.26's averaged collective block), and the function-level floor, averaged over a
period, is

    E‖f̂ − f‖²_ρ = (N₀/2k_e)[Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) + tr H²] + O(k_e⁻²)

(at time t its fast part is (N₀/2k_e)tr(HĴ(t)ᵀĴ(t)), whose oscillation over the period is of the floor's own order):

the parameters' part, each mode's share multiplied by (h_j + η)/(h_j + η + σ₀) — the fixed prior's Tikhonov filter factor h_j/(h_j + σ₀)
at η = 0, the factor 1 at η = ∞, where the part is (N₀/2k_e)tr Q whatever the spectrum — and the jitter of the frequencies that the loop
passes, independent of the prior (the compensation error at the plant's current point has ⟨|Ĵ|⁴⟩ in place of tr H²). Mode j's parameter
variance is (N₀/2k_e)Q_jj(h_j + η)/(h_j(h_j + η + σ₀)) — from Q_jj/(2(h_j + σ₀)) to Q_jj/(2h_j) in units of N₀/k_e — plus the fast part
(N₀/2k_e)h_j.

(iii) *(the errors in variables)* To second order the noise inside the regressor Φ̃(y) exerts a mean force on the collective
coordinates. Its part carried by the noise's own second moments, at fixed parameters, is k_eF in the slow time with
F := ⟨Ĵ′ᵀE[ε_m²] + ĴᵀE[e]⟩ = O(Eν²) (Ĵ′ := ∂_xĴ, ε_m := y − x_m the measured error; for an Ornstein–Uhlenbeck noise of rate B,
E[ε_m²] = Eν²(a² + k_eB)/(k_e(k_e + B)) → Eν² as B/k_e → ∞); the parameters' first-order fluctuation adds a force of the same order in N₀
— the compensation's curvature and slope–parameter terms in the error equation, the correlation ∇²_ϑΦ̃E[ε_mδϑ] in the law — negligible
against F for the measured regressor, not for the desired one (Proposition D‴.31(iv)). The stationary response is k_e(H + σ₀I)⁻¹F under
the fixed prior and k_eH⁻¹F under the adaptive prior, for every η > 0: the adaptive prior gives up the leak's bound on the response to a
persistent force in the weakly excited modes — the robustness half of the classical trade-off, here against the noise, and the
lack-of-excitation drift of adaptive control ([S57]) — and its function-level bias k_e²FᵀH⁻¹F grows like k_e²(Eν²)² while the noise is
white well beyond the loop's bandwidth (k_e ≪ B), and saturates for k_e ≫ B, where E[ε_m²] → Eν²B/k_e. The part of F that the adaptive
prior amplifies is the regressor's own, ⟨Ĵ′ᵀE[ε_m²]⟩; evaluating the regressor at the reference removes it (Proposition D‴.31).

Hence η buys the bias-removal rate λ_j^−(η) at the price of the factor (h_j + η)/(h_j + η + σ₀) in the floor and of the H⁻¹ response to
the errors in variables, and k_e lowers the floor like 1/k_e while it raises the errors-in-variables bias like k_e² (k_e ≪ B).

*Proof.* (i) The equations of the particles of a class are identical — the same y, the same measured error and the same prior atom (the
class mean at η = ∞, the same z_k otherwise) enter them — and so are their initial values; the field is locally Lipschitz and the cells
are constant near the lift, so the trajectories coincide, and the class means follow the weighted K₀-atom loop with the same noise.
(ii) Linearised at Γ_w in Corollary D‴.26's whitened coordinates: δė = −k_eδe − Ĵy − (k_e + a)ν (from f(x) − f̂(x + ν) = −Df̂[δΘ] − aν),
ẏ = Ĵᵀ(δe + ν) − (σ₀/k_e)(y − ζ), ζ̇ = −(η/k_e)(ζ − y) (the regressor's own noise is of second order). The measured error ε_m = δe + ν
obeys ε̇_m = −k_eε_m − Ĵy + ν̇ − aν, so ε_m = −(d/dt + k_e)⁻¹[Ĵy] + n, n := (d/dt + k_e)⁻¹[ν̇ − aν] (written ((d/dt − a)/(d/dt + k_e))ν):
the first part is Corollary D‴.26's slaving, the second drives ẏ = Ĵᵀn. With ξ := (d/dt + k_e)⁻¹ν, of variance N₀/(2k_e) for white noise,
n = ξ̇ − aξ + (d/dt + k_e)⁻¹[ȧξ] and Ĵᵀn = d(Ĵᵀξ)/dt − (Ĵ̇ + aĴ)ᵀξ + O(ξ/k_e): y = Ĵᵀξ + u, the fast jitter Ĵᵀξ and the slow part u driven
by −(Ĵ̇ + aĴ)ᵀξ, white at the slow scale with the intensity (N₀/k_e²)Q per unit time, (N₀/k_e)Q in the slow time; the averaged slow
equation dY = Ā₀Y ds + (N₀/k_e)^{1/2}(Q^{1/2}dW, 0) (stochastic averaging, [S53]; exactly, for the linearised loop, Lemma D‴.26a(v)) has the
stated stationary covariance. In the eigenbasis
of H the block Ā₀ is a direct sum of the pair blocks [[−(h + σ₀), σ₀], [η, −η]] (η = 0: −(h + σ₀); η = ∞: −h), whose Lyapunov
equation with the noise q on y gives Σ_ζζ = Σ_yζ = ηΣ_yy/(h + η) and Σ_yy = q(h + η)/(2h(h + η + σ₀)); Q's off-diagonal entries in that
basis enter only the covariances of different modes, which tr(HΣ_yy) = Σ_jh_jΣ_yy,jj does not see. The fast jitter gives
E[(Ĵᵀξ)ᵀH(Ĵᵀξ)] = (N₀/2k_e)tr H² over a period, and E[(Ĵ(t)Ĵ(t)ᵀξ)²] = (N₀/2k_e)⟨|Ĵ|⁴⟩. (iii) At fixed parameters the second-order part
of E[ε_m∇_ϑΦ̃(y; ϑ)] at the lift is ∂_x∇_ϑΦ̃·E[ε_m²] + ∇_ϑΦ̃·E[e], E[e] from the error equation's second-order terms ½∂²_xf·e² −
½∂²_xf̂·(e + ν)² — one second derivative for f and f̂, by the hypothesis (otherwise ½(∂²_xf − ∂²_xf̂)Var(e) remains, of the same order);
the parameters' fluctuation adds ∇²_ϑΦ̃·E[ε_mδϑ] in the law and −∂_xDf̂[E[(δe + ν)δΘ]] − ½E[D²_Θf̂[δΘ, δΘ]] in the error equation, of the
same order in N₀ (computed in the numerics from the periodic covariance: the full second-order mean). For the Ornstein–Uhlenbeck noise,
quasi-statically along the orbit, E[νe] = −(k_e + a)Eν²/(k_e + B) and
Var(e) = (k_e + a)²Eν²/(k_e(k_e + B)), which give E[ε_m²] = Eν² + 2E[νe] + Var(e) as stated and E[e] =
−(∂²_xf̂/(2k_e))Eν²(B − k_e − 2a)/(k_e + B). A constant slow force F_s on the pair block of a mode has the stationary point ζ = y and
(h + σ₀)y − σ₀ζ = F_s, i.e. y = F_s/h for every η > 0, and y = F_s/(h + σ₀) under the fixed prior. ∎

*Status.* (i) proved, exactly, as invariance; the threshold formal at the lift (the separations' averaged second-order drift there), the
drift along the noisy trajectory computed from the second-order transient and checked against the real loop's separation rates on
antithetic noise pairs, the long-time threshold open; (ii) proved for the linearised loop under white noise, for the floor's period
average (Lemma D‴.26a(v): the
slow coordinates' drift decoupled exactly, their noise εbν, at leading order the high-passed one; at η = 0 with K = Ĵ, Λ = σ₀I), with the
relative remainder O(1/k_e), the closed form exact for the averaged system; the OU
noise of the numerics and the real loop formal — and checked against the exact periodic covariance of the linearised loop: the deficit 2.4%, 0.61%,
0.14%, 0.028% at k_e = 25–200 (η = 0) is dominated there by a term in k_e⁻², but the relative error is O(1/k_e) — the slope a in the noise
input −(k_e + a)ν, a rate that does not change sign under time reversal, breaks the parity of Corollary D‴.26(c) (with the noise input −k_eν
the deviation is numerically even in 1/k_e, to k_e = 800) — and the ratio crosses 1 near k_e = 500; (iii) formal — a quasi-static
second-moment computation for the
noise's own moments, the parameters' part computed numerically — its stationary response exact for the averaged system, the exact
linearised loop's 10.7% (adaptive) and 2.6% (fixed) higher at k_e = 25, an O(k_e⁻²) excess; the real nonlinear loop checked by Monte Carlo
against the sum of the two. Classical: measurement noise in adaptive control and the errors-in-variables bias of estimates from noisy
regressors (bias-compensated least squares and instrumental variables in identification; bounded disturbances and robust modifications,
[4]); using the reference in the regressor, as the desired compensation adaptation law does ([S53]); the Lyapunov covariance of a
constant-gain adaptive algorithm as its diffusion approximation (Benveniste–Métivier–Priouret), the LMS misadjustment's dependence on the
input's spectrum through a trace only (Widrow et al.: the η = ∞ floor's tr Q), the leaky LMS's steady state (Mayyas–Aboulnasr: the fixed
prior's factor h_j/(h_j + σ₀)), the common noise of mean-field theory (Carmona–Delarue–Lacker), of which the coherence is an instance, and
the drift of the weakly excited directions under a persistent disturbance (Anderson), of which the H⁻¹ response is ([S57]). New, in the
transport-leak setting: the interpolation (h_j + η)/(h_j + η + σ₀) of the floor in the prior's mobility, the driving matrix
Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩ — the sensitivity-filtered measurement the law sees, built on Corollary D‴.26's derivative Gram — and their
combination with the coherence of the K₀-atom identification and its noise threshold at the lift. Checked: `tests/test_measurement_noise.py` (the
closed form against the Lyapunov solution, the trace identity at η = ∞ and the monotonicity in η; the exact periodic floors of the
linearised loop at k_e = 100 for η = 0, σ₀ and ∞; a lift m = 6 under a noise path staying exactly replicated and following the K₀-atom loop
driven by the same path; a separation of the weakest class decaying at the noise std 0.01 and growing at 0.05 — a sign test, the rates
within 10⁻⁴ of the lift's prediction —; the separations' drift against an independent quadrature of its two moments, and along the
second-order transient the mean parameters moved by the averaged force and the mean error carrying the tracking error of the developing
deviation; the OU noise's stationary variance and correlation, and the law seeing the measured error at the measured point; the responses
(H + σ₀I)⁻¹ and H⁻¹ to a constant force, the quasi-static force against the full second-order mean's noise part to 2·10⁻³, the
parameters' part below 10⁻⁴ of the bias; the exact periodic mean against the filtered-Gram response under both priors, the fixed prior's
branch and its leak in the bias; and the mean function's curvature term against an independent contraction).


*Numerics* (`exp_measurement_noise.py` → `results/chaos/measurement_noise.md`). The three-bump instance (σ₀ = 0.15σ_crit, x_m = 2 sin t,
k_e = 25) at the lift of the exact representation of the matching weights q = 1/3, the prior at its atoms. (A) The linear theory: H's nine
eigenvalues span 3.06·10⁻⁴–0.195; tr Q = 8.784, of which 8.105 is the derivative Gram ⟨Ĵ̇ᵀĴ̇⟩ and 0.679 the slope term ⟨a²ĴᵀĴ⟩ (the cross
term vanishes, as on every one-dimensional reference: an exact differential); Q_jj/h_j ranges from 104 (mode 1) to 5.5 (mode 6). The
parameters' floor per unit N₀ is 0.1597 under the fixed prior, 0.1627 at η = σ₀, 0.1690 at 10σ₀ and 0.1757 = tr Q/(2k_e) at η = ∞ — the
adaptive prior costs 10% here, the strong modes carrying the floor — while the weakest mode's parameter variance grows from 0.0865 to 2.09,
×24.1 in the averaged system (×28.0 in the exact linearised loop at k_e = 25, Part B). The fast parts are tr H² = 0.0657 and ⟨|Ĵ|⁴⟩ = 0.299.
The second-order force of the noise's own moments (OU noise, B = 100) has the averaged function-level stationary bias 491 (fixed prior) and
2397 (adaptive, ×4.88) per unit (Eν²)², equal to the floor of the parameters at the noise std 0.0026 and 0.0012 respectively and dominant
above; the full second-order mean (the periodic covariance contracted with the true field's second derivatives) reproduces the force to
1.7·10⁻⁴, the parameters' fluctuation changing the bias by less than 10⁻⁵ (2397.66 → 2397.68; the mean function's 2397.65), and the exact
periodic second-order mean of the linearised loop gives 504 and 2653 (×5.27). The averaged bias grows in k_e with the local exponent
1.93–1.99 for B ≥ 20k_e but 1.69 and 1.47 at B = 8k_e and 4k_e, and saturates beyond. At the lift the separations' second-order drift
has λ_max(S_k) = 0.0148, 0.0101 and 0.298 per unit Eν² (real time) in the three classes: there the third class's separations stop
decaying at the noise std 0.031 (1.5% of the reference's amplitude), and at Part C's stds 0.003 and 0.01 the margin σ₀/k_e is used to
0.9% and 10.5%. (B) The
exact periodic stationary covariance of the linearised loop (the monodromy and a discrete Lyapunov equation) against the prediction, at
η = 0, σ₀, ∞ and k_e = 25–800: the function-level floor's ratio 0.976–0.980 at k_e = 25, 0.994–0.995 at 50, 0.9986–0.9988 at 100,
0.9997–0.9998 at 200, and above 1 at 800 (+1.2·10⁻⁵ at η = 0, +1.3·10⁻⁵ at η = ∞); the deviation's successive ratios per doubling are 3.87,
4.38, 5.26, 11.4 and a change of sign (η = ∞) — a k_e⁻² term dominating to k_e = 200 over a first-order one of the opposite sign — while
with the noise input −k_eν (the slope removed, the control) they are 3.67, 3.90, 3.98, 3.99, 4.00: without the slope the expansion is
numerically even (to k_e = 800; not proved).
The weakest mode's variance is 0.845 (fixed prior) and 0.980 (adaptive) of its slow-plus-fast prediction at k_e = 25, tending to 1. The
exact stationary second-order mean of the linearised loop lies 10.7%, 2.65%, 0.66%, 0.17% above the averaged bias under the adaptive prior
and 2.58%, 0.64%, 0.16%, 0.04% under the fixed one at k_e = 25, 50, 100, 200 — the weak modes' stiffness being the filtered Gram
H − k_e⁻²⟨Ĵ̇ᵀĴ̇⟩ of Corollary D‴.26(c). (C) The real loop (RK4, dt = 10⁻³, an OU noise exact at the steps, 12 noise paths of 300 time units
common to the four cases) against the transient covariance of the linearised loop plus the second-order mean deviation: at the noise std
0.003 within one standard error at every record (z from −0.99 to 0.99), the mean part 42% of the prediction at t = 300; at std 0.01 within
1.75 standard errors (z from −0.3 to 1.75, the records of a run correlated), the mean part 89% at t = 300 and the covariance part alone
8.8–8.9 times below the Monte Carlo — the second-order force is the dominant effect of the noise there. The adaptive and the fixed prior
drift apart, by 6.0% and 6.4% at t = 300 (predicted 5.3% and 6.1%): the weak modes, where the two responses part, equilibrate over k_e/h₁ ≈
8·10⁴ time units, so over 300 time units the transient checks the force (the weak modes still grow linearly), and the stationary response is
checked against the exact linearised loop (B). The lift m = 15 (classes of five) driven by the same noise path as the K₀-atom run gives the
same f̂ to rounding (identical in the stored float32, 5.6·10⁻¹⁶ apart in float64): the noise is coherent, as (i) states. (D) But invariance
is not attraction: from a separation of size 10⁻⁶ along S₃'s top eigenvector, on eight antithetic pairs of noise paths of 300 time units
per std (ν and −ν, whose pair means cancel the first-order part of the log-growth), the lift's third class separates at
−2.50·10⁻⁴ ± 3.8·10⁻⁷ per time unit at the std 0.01, +2.74·10⁻⁶ ± 3.3·10⁻⁶ at 0.03 and +4.53·10⁻⁴ ± 8.4·10⁻⁶ at 0.05. The noise-induced
part λ_eff := (rate + σ₀/k_e)/Eν² is 0.3271 ± 0.0038, 0.3175 ± 0.0037 and 0.2943 ± 0.0034, against the lift's drift λ_max(S₃) = 0.2981 and
the drift along the second-order transient from the lift, 0.3264, 0.3270 and 0.3281 (the same estimator applied to u̇ = (Eν²S₃(t) − σ₀/k_e)u). Along
the transient the mean error's rms per unit Eν² grows from 0.121 to 0.602 over the six 50-unit windows to t = 300 and λ_max(S₃(t)) from 0.307
to 0.349: the tracking error of the developing errors-in-variables deviation. The weighted fit λ_eff = c₀ + c₁Eν² has the small-noise
limit c₀ = 0.3291 ± 0.0033 — 10.4% above the lift's drift, within its standard error of the transient's 0.3264 — and the next order
c₁ = −13.8 ± 2.0; over this window the measured rate
vanishes at the std 0.0299, against 0.0308 by the lift's drift and 0.0294 by the transient's. Above it the replicated configuration
repels, its within-class spread staying exactly 0.0 from a replicated start only because each particle's update is the same arithmetic.

**Proposition D‴.31 (the adaptive prior is a washout leak; the desired regressor makes the noise's persistent force a least-squares
force).** Let the leak be Theorem D‴.23's toward an adaptive prior, dz_j/ds = −η(z_j − T̂_j) with 0 < η < ∞ (the nearest-atom form of
Proposition D‴.24(ii) with the class barycentres likewise), in the averaged loop or in the real loop (CL) of Corollary D‴.26 (s = t/k_e);
at an exact representation of representable weights w let H, Ĵ, a and the mean ⟨·⟩ over a period be Proposition D‴.30's,
‖g‖²_ρ := ⟨g²⟩ for a function g along the orbit, and P the L²_ρ-projection onto the span of Ĵ's components.

(i) *(the washout)* Along every solution and every path of the plan with T̂_j of bounded variation, T̂_j(s) − z_j(s) =
∫₀ˢe^{−η(s−τ)}dT̂_j(τ) + e^{−ηs}(T̂_j(0) − z_j(0)) — and for every measurable T̂_j (a chattering plan included) in the integrated form
z_j(s) = e^{−ηs}z_j(0) + η∫₀ˢe^{−η(s−τ)}T̂_j(τ)dτ: the prior lags the barycentre by the barycentre's exponentially forgotten increments.
For integral plans (the nearest-atom cells, the OT plan at a commensurate m) the leak on a particle, −σ₀(ϑ_i − z), is therefore the full
σ-modification of its offset from the barycentre plus, on the barycentre, the σ-modification passed through the washout filter
(d/ds)(d/ds + η)⁻¹ — zero at zero frequency, the full leak above η; with split rows (an incommensurate m) the barycentres also feel
−σ₀(I − M)z, M_jl := (m/q_j)Σ_iγ_ijγ_il, which is not washed out. Linearised at the exact representation, mode j of the collective block
has at the slow frequency ω the stiffness h_j + σ₀iω/(iω + η): the fixed prior's h_j + σ₀ for ω ≫ η, the bare h_j at ω = 0. This is
Theorem D‴.23(iii)'s bias removal — feedback through a washout filter leaves the equilibria where they are ([S54]), here the replicated
ones — and its price: no leak opposes a persistent force on the class barycentres (the separations keep the full leak).

(ii) *(the stationary identity)* Let the prior be an affine time-invariant filter of the barycentre and of a fixed anchor z_a — linearised,
ζ = 𝒦[y] + (I − 𝒦)[δ] in the whitened weighted coordinates, δ the anchor's misplacement from the representation; for example
dz/ds = −η(z − T̂) − κ(z − z_a), κ ≥ 0, with 𝒦(0) = η/(η + κ) (κ = 0 the adaptive prior, κ = ∞ the fixed prior at z_a) — and
S₀ := σ₀(I − 𝒦(0)) the leak's zero-frequency gain (the anchor entering through I − 𝒦: the filter returns the anchor when the barycentre
sits there). A persistent force F on the collective coordinates and the misplacement δ give the stationary deviation y = RF + Bδ with
R = (H + S₀)⁻¹ and B = (H + S₀)⁻¹S₀, so

    RH + B = I,    i.e. R = (I − B)H⁻¹:

the response is the unregularised one times the share not taken by the anchor (HR + B = I as well where S₀ commutes with H, as for the
family; not for a prior fixed on some coordinates or a projected prior with an active face). On every mode, whatever S₀,
h_jr_jj + b_jj = 1 (r_jj := v_jᵀRv_j, b_jj := v_jᵀBv_j): a vanishing bias forces the response 1/h_j. So no such prior both removes the
leak's bias and keeps its bound on the response to persistent forces — the trade binding on the modes with h_j < σ₀, where 1/h_j exceeds
the fixed prior's 1/σ₀ — and a remedy for the robustness loss of Proposition D‴.30(iii) has to change the force.

(iii) *(least-squares forces)* For H ≻ 0 and every S with S + Sᵀ ⪰ 0, Fᵀ(H + S)⁻ᵀH(H + S)⁻¹F ≤ FᵀH⁻¹F = min{‖g‖²_ρ : ⟨Ĵᵀg⟩ = F},
with equality at S = 0: the function-level stationary response to a persistent force is at most the squared norm of the smallest target
perturbation that exerts it (⟨Ĵᵀg⟩ is the learning force of the target f + g), and the adaptive prior attains the bound. A least-squares
force ⟨Ĵᵀg⟩ therefore moves the learned function by at most ‖g‖_ρ under every prior — under the adaptive prior by exactly Pg, so that
the stationary point's projected residual P(g − Ĵy) vanishes.

(iv) *(the desired regressor)* Let the adaptive law evaluate its regressor at the reference, ϑ̇_i = G∇_ϑ((y − x_m)ᵀΦ̃(x_m; ϑ_i)) −
(σ₀/k_e)(ϑ_i − z_{π(i)}) (DR), with the compensation f̂_Θ(y) as before or also at the reference, u = ẋ_m − k_e(y − x_m) − f̂_Θ(x_m) (DR*,
the desired compensation adaptation law [S53]). Then
(a) without noise the lift Γ_w is a solution and the linearisation at it is Corollary D‴.26's — for DR exactly (the regressor's variation
multiplies e ≡ 0), for DR* with −(k_e − a(t)) in place of −k_e in the error channel — and the averaged loop is unchanged at leading order
in 1/k_e (the regressor moves by x − x_m = O(1/k_e));
(b) to first order in ν, DR is Proposition D‴.30(ii)'s linearised loop driven by the same noise, path by path; DR* has the noise input
−k_eν in place of −(k_e + a)ν, and the same floor up to a relative O(1/k_e);
(c) to second order the part of the force carried by the noise's own second moments (at fixed parameters) is the least-squares force
⟨Ĵᵀg⟩ (slow time) of g := k_eE[e], the rescaled mean tracking error at the lift; for an Ornstein–Uhlenbeck noise of rate B,
quasi-statically along the orbit, g = −½∂²_xf̂·Eν²(B − k_e − 2a)/(k_e + B) for DR and g = ½∂²_xf·Eν²k_e³/((k_e − a)²(k_e − a + B)) for
DR* (k_e > max a), so |g| ≤ ½|∂²_xf|Eν²(1 + O(|a|/k_e)) for all k_e and B. By (iii) its function-level response is at most
‖g‖²_ρ ≤ ¼‖∂²_xf‖²_ρ(Eν²)²(1 + O(1/k_e)) under every prior — uniformly in k_e, in B and in the spectrum of H — and under the adaptive
prior it leaves no projected mean tracking error. The parameters' first-order fluctuation adds a second-order force of the same order in
N₀: in the error equation the compensation's curvature ½D²_Θf̂[Σ_ΘΘ] and, for MR and DR, its slope–parameter term (DR*'s compensation,
at x_m, has none) — a least-squares force again, of a perturbation set by the parameters' covariance; the curvature also shifts the mean
function, and the loop removes its tangent part —, in the learning law the correlation ∇²_ϑΦ̃E[ε_mδϑ] of the regressor's parameter
dependence with the measured error — an errors-in-variables term in the parameters, which the reference in the regressor does not remove
and which is not the learning force of a mean tracking error, so no bound in the noise's moments applies to it. Its size is set by
the parameters' stationary covariance (mode j's variance (N₀/2k_e)Q_jj(h_j + η)/(h_j(h_j + η + σ₀)), unbounded as h_j → 0 at η > 0), so
the bound ‖g‖²_ρ is the noise part's, not the total's. With the measured regressor the noise part is F₁ + ⟨Ĵᵀg⟩, F₁ := k_e⟨Ĵ′ᵀE[ε_m²]⟩ —
the regressor's errors in variables proper, the learning force of no mean tracking error — and the adaptive prior's stationary point keeps
the projected mean tracking error P E[e] = −ĴH⁻¹F₁/k_e: what grows like k_e² while k_e ≪ B is k_e times it, of squared norm F₁ᵀH⁻¹F₁, the
function-level bias; the mean tracking error itself, of squared norm F₁ᵀH⁻¹F₁/k_e² ∝ (B/(k_e + B))², stays bounded.

Hence with the desired regressor the first-order floor still falls like 1/k_e while the noise part of the errors-in-variables bias stays
bounded, and on the instance the parameters' part's own response falls like k_e⁻² — a one-dimensional reference, where the law's leading
term, the period average of a function of x_m times ẋ_m (a line integral ∮M∂_xĴ dx along the orbit), vanishes; on a curve in ℝⁿ, n > 1,
it need not: the trade-off in k_e of Proposition D‴.30 disappears, and what is left
is the plant's state jitter E[e²] = Eν²(k_e + a)²/(k_e(k_e + B)) (DR; DR*: k_e²Eν²/((k_e − a)(k_e − a + B))) — the classical price of a
high gain under measurement noise.

*Proof.* (i) The integrated form is the variation of constants for dz_j/ds = −η(z_j − T̂_j); for T̂_j of bounded variation (it jumps
where the plan does, finitely often on bounded intervals), u := T̂_j − z_j satisfies du = dT̂_j − ηu ds, whence
u(s) = e^{−ηs}u(0) + ∫₀ˢe^{−η(s−τ)}dT̂_j(τ); and ϑ_i − z = (ϑ_i − T̂) + u for an integral plan (for split rows the targets are mΣ_jγ_ijz_j,
and averaging the leak over a class gives −σ₀(T̂ − Mz)). Linearised, ζ = η(iω + η)⁻¹y on each mode (the prior and the leak are scalars in
the whitened weighted coordinates, so they commute with H), and σ₀(y − ζ) = σ₀iω(iω + η)⁻¹y. (ii) At zero frequency
ζ = 𝒦(0)y + (I − 𝒦(0))δ, so the leak on y is −σ₀(y − ζ) = −S₀(y − δ), and the stationary equation Hy + S₀(y − δ) = F gives
y = (H + S₀)⁻¹(F + S₀δ); then RH + B = (H + S₀)⁻¹(H + S₀) = I, and on the diagonal in H's eigenbasis h_jr_jj + b_jj = 1. For the family,
the prior's stationary equation η(y − ζ) = κ(ζ − δ) gives 𝒦(0) = η/(η + κ). (iii) By congruence with H + S, (H + S)⁻ᵀH(H + S)⁻¹ ≤ H⁻¹ ⇔
H ≤ (H + S)ᵀH⁻¹(H + S) = H + S + Sᵀ + SᵀH⁻¹S, true for S + Sᵀ ⪰ 0. The map g ↦ ⟨Ĵᵀg⟩ has the adjoint F ↦ ĴF, so the smallest g with ⟨Ĵᵀg⟩ = F is ĴH⁻¹F, of squared norm
FᵀH⁻¹F; for F = ⟨Ĵᵀg⟩, FᵀH⁻¹F = ‖Pg‖²_ρ ≤ ‖g‖²_ρ, and at S = 0, Ĵy = ĴH⁻¹⟨Ĵᵀg⟩ = Pg. (iv) (a) At e = 0 both laws vanish, and
their variations, δe∇_ϑΦ̃(x; ϑ) + e(⋯), coincide at x = x_m; for DR*, f(x) − f̂(x_m) = aδe − Ĵy to first order. Off the lift e is
slaved, e = r(x_m)/k_e + O(k_e⁻²), so the two regressors differ by O(1/k_e). (b) The regressor's dependence on y enters the law
multiplied by ε_m, at second order; for DR* the compensation no longer carries ν. (c) At fixed parameters the mean of ε_m∇_ϑΦ̃(x_m; ϑ)
is ∇_ϑΦ̃(x_m)E[e], since Eν = 0; the parameters' fluctuation adds ∇²_ϑΦ̃(x_m)E[ε_mδϑ] in the law and the compensation's curvature and
slope–parameter terms in the error equation (Proposition D‴.30's proof). For DR, E[e] is Proposition D‴.30's (the same compensation); for DR*,
ė = (a − k_e)e + ½∂²_xf e² − k_eν + O(3), whose linear part has Var(e) = k_e²Eν²/((k_e − a)(k_e − a + B)) for the OU noise, and the
mean balance is (k_e − a)E[e] = ½∂²_xf Var(e). The bounds: |B − k_e − 2a| ≤ (k_e + B)(1 + 2|a|/(k_e + B)) and
k_e³/((k_e − a)²(k_e − a + B)) ≤ (k_e/(k_e − a))³. With the measured regressor, at the adaptive prior's stationary point
⟨Ĵᵀ(g − Ĵy)⟩ + F₁ = 0 and the mean of the error equation gives k_eE[e] = g − Ĵy at second order, so k_eP E[e] = −ĴH⁻¹F₁. ∎

*Status.* (i) proved, exactly (every solution of the real or the averaged loop; the differential form for T̂ of bounded variation, the
integrated one always; the leak's decomposition for integral plans); (ii) proved as RH + B = I and per mode — linear algebra on the
averaged linearised block; HR + B = I where S₀ commutes with H; (iii) proved (S + Sᵀ ⪰ 0, H ≻ 0); (iv)(a) proved at the lift and for its
linearisation, the leading-order equality of the averaged loops by Theorem J's averaging argument (not rewritten); (iv)(b) the first-order
identity proved for the linearised equations, the floors themselves with Proposition D‴.30(ii)'s status; (iv)(c) formal and split — the
noise part by Proposition D‴.30(iii)'s quasi-static second-moment computation, its consequences exact for the averaged system; the
parameters' part computed numerically (the full second-order mean), not in closed form. Classical: reference regressors and the desired
compensation adaptation law ([S53]; DR* without Sadegh and Horowitz's nonlinear damping is a local statement), washout filters ([S54]),
the adaptive law kept near a low-pass filtered copy of itself — Yucelen and Haddad's low-frequency learning, which has the washout
structure ([S58]) — and elastic averaging's moving-average centre ([S49]); the identity, per mode the Tikhonov filter factor plus the
residual factor (the Kalman-gain split of a regularised estimate); the bound (iii), the minimum-norm solution with the contraction of the
regularised fit; and the split of the noise's force into the regressor's errors in variables proper and the output's, of which only the
mean biases; and the parameters' own errors in variables, the classical second-order bias of a constant-gain stochastic approximation
from its own fluctuation — the mean field's curvature in the parameters contracted with their stationary covariance, plus the
noise–parameter correlation (Benveniste–Métivier–Priouret; for constant-step stochastic gradient descent Dieuleveut–Durmus–Bach, who remove
it by Richardson–Romberg extrapolation; the curvature bias of nonlinear least squares, Box's, the analogue of ½D²_Θf̂[Σ_ΘΘ]: [S57]). New
are their transfer to the transport leaks — the washout reading of the adaptive prior in the mass-constrained, mean-field setting,
uniform in m — and the instance's quantification, with the observation that the desired regressor does not remove the parameters'
part. Checked: `tests/test_eiv_remedy.py` (the washout identity along a real-loop solution at finite η, the
barycentre's velocity integrated alongside, and on the averaged loop across a plan switch, the jump's share large; the collective block's
transfer function against the washout formula, RH + B = I along the prior family and, R and B from their dynamic blocks, for priors whose
zero-frequency leak does not commute with H — fixed on the amplitudes, one active face, a matrix anchor gain, a random one with
K + Kᵀ ⪰ 0 — where HR + B = I fails; the response bound under random leaks S ≥ 0, the smallest equivalent perturbation and
the projection; the lift exact under DR and DR*, on a common noise path DR's run within the second order of MR's, and the regressor's
point in the law; DR's and DR*'s noise forces least-squares, their pointwise bounds on g, MR's F₁ three orders of magnitude above; the
quasi-static forces against the full second-order mean's noise part, the parameters' part above 1% for DR and DR* and below 10⁻⁵ for MR).


*Numerics* (`exp_eiv_remedy.py` → `results/chaos/eiv_remedy.md`). The instance of Proposition D‴.30 with an Ornstein–Uhlenbeck noise of rate
B; MR the measured regressor and compensation, DR the regressor at the reference, DR* both at the reference. (A) The averaged system's
biases per unit (Eν²)². From the noise's own moments, the adaptive prior's function-level bias at k_e = 25, B = 100 is 2397 (MR), 0.471 (DR)
and 0.0526 (DR*); with the parameters' fluctuation (the full second-order mean) 2398, 0.486 (+3.3%) and 0.0481 (−8.7%), the mean function's
2398, 0.476 and 0.0511. Over k_e = 25–200 and B = 100, 1000 the parameters' part changes MR's bias by at most 0.03%, DR's by −0.75% to +3.3%
except at k_e = B = 100 — where DR's noise part nearly cancels (4.5·10⁻⁵) and the parameters' part raises it by 56% — and DR*'s by −8.7% to
−0.15%; its own response is 3.9·10⁻⁴, 1.1·10⁻⁴, 2.9·10⁻⁵ and 7.4·10⁻⁶ at k_e = 25, 50, 100, 200 (B = 100), like k_e⁻² — the cancellation
of the law's leading term on the one-dimensional reference —, while its share of DR's bias does not fall monotonically. MR over DR (full)
is 4.9·10³ at k_e = 25, B = 100,
and 3.0·10³–2.1·10⁸ over the grid. MR's bias grows like (k_eB/(k_e + B))² (×46.6 from k_e = 25 to 200 at B = 1000, the law ×46.7, k_e² ×64),
while DR's and DR*'s stay at most 1.18, below ¼‖∂²_xf‖²_ρ = 1.31. DR's noise part attains its bound ‖g‖²_ρ to 0.01% at k_e = 25 — g lies
almost in the tangent space, the second derivative of a Gaussian bump being in the span of its amplitude and width derivatives. The rescaled
projected mean tracking error k_e²‖P E[e]‖²_ρ under the adaptive prior is 2394 for MR (0 for DR and DR* by construction), the mean tracking
error itself 3.83 at k_e = 25 and 0.66 at 200 (B = 100). The noise level at which the adaptive prior's full bias reaches its floor rises
from the std 0.0012 (MR) to 0.085 (DR) and 0.27 (DR*) at B = 100 (stds above 0.1 lie outside the small-noise regime). (B) Along the prior
family at η = σ₀ (κ/η = 0, 0.1, 1, 10 and the fixed prior; S₀ from 0 to σ₀) the weakest mode's response falls from 3269 to 136 while its
bias factor rises from 0 to 0.959, with h₁r₁ + b₁ = 1 and max|RH + B − I| ≤ 1.1·10⁻¹⁴; the trade binds on the 3 of 9 modes with h_j < σ₀;
MR's errors-in-variables bias falls only from 2397 to 491, and DR's from 0.471 to 0.211 — MR over DR 2.3·10³–5.1·10³ along the family. For
priors whose zero-frequency leak does not commute with H, R and B again from their dynamic blocks — the prior fixed on the amplitudes or
on the centres and adaptive at η = σ₀ elsewhere, one active face of the projected adaptive prior (η = ∞), the anchor gain K = ηP_amp, a
random anchor gain with K + Kᵀ ⪰ 0 — the block's R and B equal (H + S₀)⁻¹ and (H + S₀)⁻¹S₀, S₀ by elimination of the prior, to
1.9·10⁻¹⁵: RH + B = I holds to 1.4·10⁻¹⁴ and per mode to 8.9·10⁻¹⁵, while max|HR + B − I| is 0.47–0.79. (C) The real loop
(k_e = 25, B = 100, std 0.01, η = ∞, RK4 dt = 10⁻³, 12 runs of 300 time units, the three laws on common noise paths): at t = 300,
E‖f̂ − f‖²_ρ is 2.43·10⁻⁶ under MR against 2.80·10⁻⁷ under DR and 2.75·10⁻⁷ under DR*, a ratio of 8.7 (predicted 9.2); DR and DR* sit at
their covariance part, their mean parts 0.07% and 0.008% of it. Against the prediction the z-scores run from −2.04 to −0.04 (MR; Proposition
D‴.30's runs of the same loop, on other paths, lay 7% above it) and from −1.8 to 0.82 (DR, DR*). The run mean of f̂_MR − f̂_DR, in which the
first-order responses cancel path by path, equals the predicted Ĵ(ȳ_MR − ȳ_DR) to a relative L²_ρ error of 0.5–5.3% at every record,
0.43–1.83 times the run mean's relative sampling error; the averaged prediction carries its own O(1/k_e) error, and the test checks F₁
(through MR − DR), not DR's own force, whose mean part is 7·10⁻⁴ of its covariance part.

**Proposition D‴.32 (the cohesion leak: a leak on the unexcited subspace, without a prior).** Let ψ: [0, ∞) → [0, 1] be C², non-increasing,
ψ(0) = 1 — of finite range, ψ(z) = (1 − z)³₊ (the range R), or Gaussian, ψ(z) = e^{−z/2} (the scale R) — and let the particles attract each
other instead of leaking toward a prior:

    dϑ_i/ds = G⟨∇_ϑ(r_ΘΦ̃(·;ϑ_i))⟩_ρ − (σ₀/m)Σ_jψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²)(ϑ_i − ϑ_j)                                          (𝔄_m^coh)

(in the real loop (CL), σ₀/k_e in place of σ₀ and the learning law in place of the averaged one).

(i) *(an energy; the network's own modes)* (𝔄_m^coh) is the G-gradient flow of J^coh_m(Θ) := mL̄(μ_Θ) + (σ₀m/2)∫∫W dμ_Θdμ_Θ, W(u) :=
R²ω(‖u‖²_{G⁻¹}/R²) with ω′ = ψ/2 (so ∇²W(0) = G⁻¹), which decreases along every solution; its mean-field form is the Wasserstein gradient
flow (in the metric G⁻¹) of L̄(μ) + (σ₀/2)∫∫W dμdμ — an aggregation equation in the external potential −⟨fΦ̃(·;ϑ)⟩_ρ with two interactions,
the fit's ⟨Φ̃(·;ϑ)Φ̃(·;ϑ′)⟩_ρ and the cohesion W. The leak force on ϑ_i is σ₀ψ̄_i(m_i − ϑ_i), ψ̄_i := (1/m)Σ_jψ_ij the kernel mass of the
cloud at ϑ_i and m_i := Σ_jψ_ijϑ_j/Σ_jψ_ij the kernel-weighted mean: the mean-shift step of the cloud's own kernel density, so the leak
pulls the network toward the modes of its own parameter distribution — no prior, no atom count, no plan, no cells. The self-interaction
vanishes (ψ(0)·0), so the empirical measure of (𝔄_m^coh) solves the mean-field equation exactly, and Theorem D's finite-horizon propagation
of chaos holds (the velocity field is Lipschitz in ϑ and in μ for W₁, Dobrushin's estimate). With the Gaussian kernel J^coh_m is
real-analytic (the saturated network is), so every bounded solution converges to a single stationary point (Łojasiewicz's theorem, as for
Proposition D‴.20). The forces sum to zero: the leak is internal and gives no a-priori bound.

(ii) *(no bias at any clustered configuration)* Let Θ be replicated by K classes with atoms ϑ_1, …, ϑ_K and masses w_k = n_k/m. With the
finite range and the atoms pairwise at least R apart, the leak vanishes on Θ, and Θ is stationary iff (ϑ_k) is a critical point of the
unregularised weighted loss L_w := ½⟨|f − Σ_kw_kΦ̃(·;ϑ_k)|²⟩_ρ — for every K and every w, with no prior and no K₀ fixed in advance; in
particular the lift of an exact representation ϑ^e(w) is stationary at every representable sample (w ∈ W, Proposition D‴.24), as under the
adaptive prior (Theorem D‴.23(iii)). With the Gaussian kernel the classes attract each other with the forces −σ₀w_lψ(d²_kl/R²)(ϑ_k − ϑ_l)
(d_kl the atoms' distance), and at a nondegenerate critical point (Ã_w ≻ 0, the whitened weighted Hessian of L_w) the stationary point
continues to one at the whitened distance Ã_w⁻¹b + O(|b|²) from it, b the whitened inter-class force: a bias at most of the order
σ₀‖Ã_w⁻¹‖max_kΣ_lw_ld_kle^{−d²_kl/(2R²)}, exponentially small in (min d_kl/R)². No bias is the finite range's property; the Gaussian kernel
pays this small bias for the merging of (v).

(iii) *(the leak acts exactly on the unexcited subspace)* At a configuration as in (ii) (the finite range), in the whitened weighted
coordinates, the linearisation of (𝔄_m^coh) is block diagonal: on the separations of class k (Σ_{i∈k}u_i = 0) it is −(σ₀w_kI − R^G_k), with
multiplicity n_k − 1, R^G_k := G^{1/2}⟨r∇²_ϑΦ̃(·;ϑ_k)⟩_ρG^{1/2} the residual curvature; on the class means it is −Ã_w, the linearised
unregularised weighted K-atom flow, the leak absent. The separations are the directions in which the plant does not see the network (Df̂ = 0
on them), and where H_L(w) ≻ 0 they are all of them: the leak acts exactly on the unexcited subspace. At an exact representation the
separations of class k decay at σ₀w_k and the collective modes at h_j(w) (the eigenvalues of H_L(w)), for every w ∈ W; the class means'
stationary response to a persistent force F is Ã_w⁻¹F — the adaptive prior's (Proposition D‴.31(i)–(iii)): no leak opposes it. The restoring
force is linear only near the class: under the compact kernel the pair force ψ(d²/R²)d is largest at d = R/√7 and vanishes at R, so a
separation pushed beyond R is not restored.

(iv) *(splitting, and no over-splitting)* The splitting margin of class k is σ₀w_k − λ_max(R^G_k). For a scalar amplitude per kernel at the
gauge transform ϑ^b(w) = G_wϑ^b(q) of a nondegenerate critical point ϑ^b(q) (Proposition D‴.28) it is σ₀w_k − (q_k/w_k)λ⁺_k, so a class can
split iff σ₀w_k² < q_kλ⁺_k — the leak's strength sets the resolution at which a class becomes unstable where the residual's curvature asks
for a kernel, light classes first; at an exact representation λ⁺ = 0 and no class splits, whatever σ₀ > 0: an exact representation is never
over-split (as linear stability). The unstable separations of a class are (n_k − 1)-fold degenerate, so the linear theory does not choose
which particles leave; the kernel's concavity (ω″ = ψ′/2 < 0) makes shedding single particles the cheapest split (Numerics (G)).

(v) *(twins and dead classes under the Gaussian kernel)* At a degenerate exact representation — twins (two classes with the same centre and
width) or a dead class (a_k = 0) — the fit is neutral along a curve (the transfer of amplitude between the twins at fixed w_ka_k + w_la_l;
the dead class's centre and width), and along it only the cohesion acts. For twins with no third class within the kernel's reach, the fit
relaxed quasi-statically along the curve, the twins' amplitude difference Δα contracts at the rate

    σ₀ψ(Δα²/(G_αR²))·(w_ka′_k + w_la′_l)²/(w_ka′_k² + w_la′_l²),   a′ := da/dα,

so they merge, for unsaturated twins on the time scale G_αR²e^{Δα²/(2G_αR²)}/(σ₀(w_k + w_l)Δα²) (Δα² ≫ G_αR²). A third class within reach of
one twin can win: twins need not merge, and a degenerate exact representation can be a saddle of J^coh_m that the flow leaves toward a
biased configuration (Numerics (I)). A dead class moves under the attraction of the live ones, a Gaussian mixture in its coordinates whose
weights carry the live classes' amplitude coordinates (w_le^{−α_l²/(2G_αR²)}); when the components are farther apart than a few R its local
minima lie near them, where it becomes a twin of a live class — on the time scale R²e^{d²/(2R²)}/(σ₀wd²), d ≫ R its distance from the
nearest live class and w that class's mass. With the finite range, dead and duplicated classes farther than R from the others persist.

(vi) *(the real loop)* In (CL) with the cohesion leak (finite range) the lift of an exact representation ϑ^e(w) by classes pairwise at least
R apart is a T₀-periodic solution with e ≡ 0 at every representable sample; its Floquet exponents are σ₀w_k/k_e exactly on the separations
of class k (multiplicity (n_k − 1)p) and, on the collapsed subspace, those of the unregularised weighted K-atom loop, which are Corollary
D‴.26's at η = ∞ — (1/k_e)[h_j(w) + c_j(w)/k_e² + O(k_e⁻⁴)] under the derivative condition (Lemma D‴.26a). So perfect tracking is
exponentially stable (locally, for each m) at the rate (1/k_e)min{σ₀min_kw_k, h_1(w) + c_1(w)/k_e² + …}. The replicated configuration is
invariant under a measurement noise common to the particles, and Proposition D‴.30(ii)–(iii) hold at η = ∞ (the floor (N₀/2k_e)[tr Q + tr
H²], the errors-in-variables response H⁻¹), the separations' noise threshold scaled by w_k in the variance (√w_k in the standard deviation)
— the softest of the leaks; with the desired regressor the regressor's own term leaves their noise drift, whose threshold is then bounded
below uniformly in k_e (Corollary D‴.33). Under a dither τ (the averaged loop with √(2τG) per particle) the separations of class k are an
Ornstein–Uhlenbeck process of stiffness σ₀w_k: the particles' stationary variance about the class mean is (1 − 1/n_k)τ/(σ₀w_k) per whitened
coordinate, against (1 − 1/n_k)τ/σ₀ under a leak toward a prior and the free diffusion 2τs without a leak.

*Proof.* (i) W is even, so ∂_{ϑ_i}[(σ₀/(2m))Σ_{j,l}W(ϑ_j − ϑ_l)] = (σ₀/m)Σ_j∇W(ϑ_i − ϑ_j), and G∇W(u) = R²ω′(z)·2u/R² = ψ(z)u (z =
‖u‖²_{G⁻¹}/R²); with ∂_{ϑ_i}(mL̄) = −⟨r∇_ϑΦ̃(·;ϑ_i)⟩_ρ the right-hand side of (𝔄_m^coh) is −G∇_{ϑ_i}J^coh_m, and dJ^coh_m/ds =
−Σ_i‖∇_{ϑ_i}J^coh_m‖²_G ≤ 0. Σ_jψ_ij(ϑ_i − ϑ_j) = mψ̄_i(ϑ_i − m_i), and Σ_iΣ_jψ_ij(ϑ_i − ϑ_j) = 0 by antisymmetry. The mean-field form
replaces (1/m)Σ_j by the integral against μ: L̄(μ) = ½⟨f²⟩_ρ − ∫⟨fΦ̃⟩_ρdμ + ½∫∫⟨Φ̃Φ̃′⟩_ρdμdμ, the velocity is v[μ](ϑ) = G∇_ϑ⟨r_μΦ̃(·;ϑ)⟩_ρ −
σ₀∫ψ(‖ϑ − ϑ′‖²_{G⁻¹}/R²)(ϑ − ϑ′)dμ(ϑ′) = −G∇_ϑ(δℱ/δμ), ℱ := L̄ + (σ₀/2)∫∫W dμdμ, and on μ_Θ it is the right-hand side of (𝔄_m^coh) at the
particles, the term j = i vanishing. (ii) Inside a class ϑ_i − ϑ_j = 0; across classes ψ = 0 at distances ≥ R. The learning drift of a
particle of class k is G⟨r_w∇_ϑΦ̃(·;ϑ_k)⟩_ρ = −Gw_k⁻¹∇_{ϑ_k}L_w, which vanishes for all k iff (ϑ_k) is critical for L_w. With the Gaussian
kernel the class drift is −Gw_k⁻¹∇_{ϑ_k}L_w + b_k, b_k := −σ₀Σ_lw_lψ(d²_kl/R²)(ϑ_k − ϑ_l); in the whitened weighted coordinates its
linearisation at the critical point is −Ã_w (plus the derivative of b, of the order of b/R), and the implicit function theorem gives the
continuation with the stated first-order displacement and the bound. (iii) A separation u (Σ_{i∈k}u_i = 0, the other classes fixed) leaves
f̂ unchanged at first order — Df̂[u] = (w_k/n_k)Σ_{i∈k}∇_ϑΦ̃(·;ϑ_k)u_i = 0 — so the residual is unchanged and the learning drift varies by
GR_ku_i; the leak's variation is −(σ₀/m)Σ_{j∈k}ψ(0)(u_i − u_j) = −σ₀w_ku_i (ψ′ terms vanish at coincident particles, and across classes ψ
and ψ′ vanish at distances ≥ R). A perturbation of the class means moves the leak not at all (inside a class all particles move together;
across classes ψ = ψ′ = 0), and the learning drift's variation is the unregularised K-atom flow's; the two subspaces do not mix (the
learning drift's cross terms vanish as in Proposition D‴.5(i)). Whitening, y = √w_kG^{−1/2}·, gives the blocks; Df̂ restricted to the class
means is the weighted Jacobian whose Gram matrix is H_L(w), so where H_L(w) ≻ 0 the kernel of Df̂ is exactly the separations. At an exact
representation r = 0, R_k = 0 and Ã_w = H_L(w). The stationary response follows from the linearisation; d/dd[(1 − d²/R²)³d] = (1 − d²/R²)²(1
− 7d²/R²). (iv) Proposition D‴.28(iii): R_k(ϑ^b(w)) = (q_k/w_k)R_k(ϑ^b(q)); the separations of class k are u_i = c_iv with Σc_i = 0 along
each eigenvector v, all with the same linear rate. (v) Twins: the network depends on their amplitudes only through p = w_ka(α_k) +
w_la(α_l), and their other coordinates coincide; with u := ⟨rφ⟩_ρ common to both, α̇_k = G_αa′_ku − σ₀w_lψΔα and α̇_l = G_αa′_lu + σ₀w_kψΔα
(ψ := ψ(Δα²/(G_αR²)); no third class within reach), so ṗ = G_αu(w_ka′_k² + w_la′_l²) − σ₀ψΔαw_kw_l(a′_k − a′_l). The fit relaxes p at the
rate G_α⟨φ²⟩_ρ(w_ka′_k² + w_la′_l²), fast against σ₀ψ, so to leading order ṗ = 0 fixes u, and d(Δα)/ds = G_αu(a′_k − a′_l) − σ₀ψ(w_k +
w_l)Δα = −σ₀ψΔα(w_ka′_k + w_la′_l)²/(w_ka′_k² + w_la′_l²); for unsaturated twins z := Δα²/(2G_αR²) obeys dz/ds = −2σ₀(w_k + w_l)ze^{−z}, so
the time from z₀ ≫ 1 is e^{z₀}/(2σ₀(w_k + w_l)z₀)(1 + O(1/z₀)). A third class at the distance d from one twin adds to that twin a force of
the order σ₀w ψ(d²/R²)d, which wins when d is below the twins' own distance. A dead class: L_w does not depend on its centre and width, so
J^coh_m varies there only through its attraction to the live classes, a Gaussian mixture in those coordinates (its amplitude coordinate
fixed, the live classes' entering the weights); a mixture's modes away from its components need overlapping components, so for components
farther apart than a few R its local minima lie within O(R e^{−d²/(2R²)}) of the components, and the descent toward the nearest has the rate
σ₀wψ(d²/R²): z := d²/(2R²) obeys dz/ds ≈ −2σ₀wze^{−z}, whence the time scale as for twins. (vi) As Corollary D‴.26: along the lift e ≡ 0,
the separations do not move f̂ at first order and the second-derivative term of the learning law vanishes with e, so a separation stays a
separation with δϑ̇_i = −(σ₀w_k/k_e)δϑ_i; the collapsed subspace is invariant and on it the leak vanishes, which is the η = ∞ adaptive
prior's collapsed variational equation; Lemma D‴.26a applies to it. The noise: Proposition D‴.30(i) needs only that every particle of a
class sees the same signals, and its (ii)–(iii) only the collapsed linearisation; the separations' drift there is opposed by σ₀w_k instead
of σ₀. The dither: the offsets from the class mean are forced by (1 − 1/n_k)·2τ per whitened coordinate. ∎

*Status.* (i)–(iii) proved, exactly; (iv) proved for the stationary points and their linearisation, the split's pattern explained (the
degeneracy and the kernel's concavity), its dynamics numerical; (v) the twins' contraction rate proved to leading order in σ₀ψ against the
fit's relaxation (checked to 0.04%, Numerics (I)), their merger only without a third class within reach (a counterexample, Numerics (I));
the dead class's motion and time scale stated, its merger not proved; (vi) the lift and its separation exponents exact, the collapsed
exponents by Corollary D‴.26 and Lemma D‴.26a, the noise statements as in Proposition D‴.30 (formal where it is). Not proved: the transit
from a spread cloud — which classes form and with which masses — and with it any selection of the number of kernels: the count is the
transit's and the schedule's. Under an annealed Gaussian scale it ended at the target's count on the instance and on well-separated targets
whose transit masses were representable, not for bumps closer than the broad scale or for unrepresentable transit masses, and at a fixed
scale not at all (Numerics (C), (H), (J), (K)); the replicator weight law of Proposition D‴.34 removes the mass obstruction. Nor is the basin uniform in m proved. Checked: `tests/test_cohesion_leak.py` (the drift equal to −G∇J^coh_m by central differences and J^coh_m decreasing along the
flow, for both kernels; the forces summing to zero and equal to σ₀ × the kernel mass × the mean-shift vector; at the exact representation
for the counts (3, 2, 2) and (4, 3, 3) the drift zero and the Jacobian's spectrum the predicted one — the separations at −σ₀w_k, the class
means the unregularised K-atom flow's — to 10⁻⁹; at the best one-atom fit the top separation eigenvalue λ⁺ − σ₀w, stable at σ₀ and split at
λ⁺/2; twins contracting at the predicted rate to 1% under the Gaussian kernel and persisting beyond the compact range; the kernels' range
and scale (a pair at 0.99R attracted, at 1.01R not; the Gaussian force at R e^{−1/2} of its linear part); the clusters of known
configurations; the real loop's lift T₀-periodic and a separation contracting by e^{−σ₀w_kT₀/k_e} per period to 10⁻⁶).

*Numerics* (`exp_cohesion_leak.py` → `results/chaos/cohesion_leak.md`). The three-bump instance (σ₀ = 0.15σ_crit), the averaged slow loop,
the midpoint rule (dt = 0.25), single-linkage clusters at 0.25. (A) At each of the 15 count vectors of m = 15 in W the lift of the exact
representation ϑ^e(w) is stationary — the largest drift 2.8·10⁻¹⁷ (compact kernel, R = 1) and 6.2·10⁻¹⁰ (Gaussian, R = 1) — and the
Jacobian's spectrum is (iii)'s to 2.8·10⁻¹¹ (the finite differences' error) and to 4.6·10⁻⁹ (the Gaussian kernel's cross-class terms, which
(iii) omits): the separations at σ₀w_k (1.41·10⁻³–3.77·10⁻³), the collective modes the unregularised flow's (the slowest 2.2·10⁻⁴–4.2·10⁻⁴),
the atoms at least 5.80 apart. The Gaussian kernel's bias (Newton on the class drift, counts (5, 5, 5)) is 2.1·10⁻⁸ (whitened; the fit error
6.6·10⁻¹⁰) at R = 1, 1.3·10⁻³ (3.9·10⁻⁵) at R = 1.5 and 0.068 (2.1·10⁻³) at R = 2, the first order Ã_w⁻¹b to 6·10⁻⁸, 9·10⁻⁴ and 1.0% — and
56% off at R = 3, where the classes overlap. With the fourth bump (ε = 0.3; the three-atom approximant ϑ^b(q), the fit error 0.017, λ⁺_k =
1.00·10⁻⁴, 2.86·10⁻⁴, 1.08·10⁻⁴) the predicted margins σ₀w_k − (q_k/w_k)λ⁺_k are the Jacobian's to 2·10⁻¹¹ (the finite differences' error)
at three count vectors. At σ₀/8 the counts (6, 4, 5) have one negative margin, −1.22·10⁻⁴, and that class splits, shedding one particle: its
diameter grows from 2.2·10⁻³ to 1.0 over 10⁵ slow units, the clusters 3 → 4 (the masses 0.200 and 0.067), the fit error 0.017 → 0.0148. At
(5, 5, 5) the smallest margin, +8.7·10⁻⁶, makes that class's diameter decay at 8.8·10⁻⁶ (6.9·10⁻⁴ → 3.4·10⁻⁴ from s = 2·10⁴ to 10⁵); at σ₀
every class collapses to rounding. Twins (bump 3 carried by classes of masses 3/15 and 2/15, the weighted amplitude shared 3:1, 2.31 apart)
merge under the Gaussian kernel (R = 1: their distance 6.1·10⁻³ at s = 4000, 0 by s = 18000, the fit error 4.5·10⁻¹⁰), on (v)'s time scale
(~10³), and persist beyond the compact range (2.31 throughout). A dead class (amplitude 0, mass 2/15, at 5.19 from the nearest live class, of mass 4/15) joins it under the Gaussian kernel at R = 2 (to rounding by s = 16000; the fit error then 1.6·10⁻³, the inter-class bias at that
scale) and at R = 1.5 (at s ≈ 2.8·10⁴, (v)'s time scale 1.8·10⁴; the fit error 1.7·10⁻⁵); at R = 1, where the live classes do not overlap
and the time scale is 1.4·10⁷, it does not move in 4·10⁴ slow units, nor under the compact kernel. (I) The twins' contraction: along (A)'s twins at R = 1 over s ∈ [200, 2000], ln(Δα(200)/Δα(2000)) = 1.2949 against (v)'s rate integrated along the run, 1.2943 (0.04%). Twins
far apart do not merge on the runs' time scale: with bump 3's amplitude shared 0.9 : 0.1 (5.28 apart) and 0.97 : 0.03 (8.19) their distance
does not change in 4·10⁴ slow units ((v)'s rate at their distance 1.3·10⁻⁹ and 3·10⁻¹⁸; 1.45·10⁻⁴ for (A)'s twins). A third class within
reach wins: opposite-sign twins of bump 3 (counts (5, 3, 5, 2), amplitudes +3.8 and −3.5, 12.7 apart), the light twin at the amplitude
coordinate of the class of bump 2 (3.63 from it), jittered by 10⁻³, do not merge in 4·10⁵ slow units (their distance 11.5–12.9); the light
twin is drawn toward the class of bump 2 (its centre 1.5 → 0.88) and the fit error grows from 3.2·10⁻⁴ to 4.5·10⁻³, a biased four-class
configuration; the annealed schedule does the same (11.7 apart at the end, the fit error 5.0·10⁻³). (B) From Proposition D‴.24's warm clouds
(m = 15, 60, 240 against M = 510, four seeds, 6000 slow units) both kernels remove the floor as the adaptive prior does, without a prior:
‖f̂_m − f̂_M‖² at s = 6000 is 1.9·10⁻¹⁰–5.4·10⁻⁸ at m = 15 over the three representable samples, under both kernels and the adaptive prior
alike (the fixed prior's floor 3.0·10⁻⁴–2.4·10⁻³), and 2.4·10⁻¹⁰–7.0·10⁻¹⁰ and 5.1·10⁻¹⁴–2.2·10⁻⁹ at m = 60 and 240 (the adaptive prior's
2.4·10⁻¹⁰–6.2·10⁻¹⁰, 7.2·10⁻¹⁴–2.2·10⁻⁹; the fixed prior's floor 1.1·10⁻⁴–9.9·10⁻⁴, 4.9·10⁻⁶–8.2·10⁻⁵); the class means lie within 0.012,
0.0010 and 0.0045 (whitened) of ϑ^e(w^m). The sample of seed 3 at m = 15, the counts (8, 4, 3), lies on the boundary of W and is not
representable: there every leak but the fixed prior ends at 7.0·10⁻⁵–7.1·10⁻⁵ (the fit error 0.013; the fixed prior 5.5·10⁻³ and 0.119) — no
leak repairs an unrepresentable sample —, and the compact kernel's run ends with a fourth class (a class diameter of 0.81). (C) From broad
starts (m = 60, six seeds; m = 120, three) the annealed Gaussian cohesion leak (R = 2 for 20000 slow units, then R = 1 for 20000) ends at
three classes with the exact representation's centres and widths in all nine runs, the fit error 1.8·10⁻⁷–2.3·10⁻⁶, the masses the transit's
(0.217–0.500). The nearest-atom adaptive prior told the right count (K = 3) does as well (4.2·10⁻⁹–2.1·10⁻⁶); told K = 2 it ends at three
biased classes (0.131–0.191), told K = 4 at four, one bump split (5.4·10⁻⁴–5.6·10⁻³); the σ-modification ends at three biased classes (one
run four; 0.24), no leak at 40–48 classes (8.2·10⁻⁴–1.6·10⁻³), and the compact kernel at R = 1 without the broad stage at 17–22
(6.0·10⁻⁴–1.5·10⁻³): the transit's fragments lie farther apart than R, where the leak no longer acts. (J) The schedule, not the kernel,
gives the count: at a fixed Gaussian scale the same broad starts end at 4–6 classes (R = 1; the three bumps' classes and light fragments of
0.017–0.083, the fit error 4.8·10⁻⁴–1.6·10⁻³, six seeds) or at three biased classes (R = 2; 7.8·10⁻⁴–1.4·10⁻³, the order of (A)'s bias at
that scale, three seeds); the annealed schedules R = 1.5 → 1, 2 → 0.5 and 3 → 1 (seed 0) end at three classes, 2.4·10⁻⁸–5.9·10⁻⁶. (D) The
real loop (k_e = 25): the lift of ϑ^e(w), w = (6, 5, 4)/15, returns after a period to 7.8·10⁻¹⁶, and its 45 slow Floquet exponents — 36
separations at σ₀w_k/k_e = 7.55·10⁻⁵, 9.43·10⁻⁵, 1.13·10⁻⁴ and 9 collective — are the predicted ones to 5.5·10⁻⁶ (relative; median 1.5·10⁻⁷;
the monodromy's finite differences' error). From the warm clouds (seed 0, 3000 time units) the rms tracking error falls to 5.8·10⁻⁴ (m = 15)
and 2.4·10⁻⁴ (m = 60) under the cohesion leak, against 1.3·10⁻³ and 3.9·10⁻⁴ under the fixed prior. (E) Under the dither τ = 10⁻⁵ (m = 60
from the lift of ϑ^e(q), three seeds, 20000 slow units) the within-class spread per whitened coordinate, averaged over the second half, is
0.068 (compact) and 0.063 (Gaussian) against (vi)'s √((1 − 1/n_k)τ/(σ₀w_k)) = 0.0635 (n_k = 20) — the compact kernel softening at that
spread, ψ = (1 − z)³ ≈ 1 − 3z —, 0.036 under the fixed prior (√((1 − 1/n_k)τ/σ₀) = 0.0366) and 0.030 under the σ-modification, while without
a leak it diffuses (0.44–0.56 at s = 20000); the class means' whitened deviation is 0.020 (the fixed prior's 0.0085), the fit error 1.9·10⁻³
(1.5·10⁻³). (F) On the cogging servo of Proposition D‴.29 (σ₀ = 0.01, m = 64, broad starts, three seeds; the Gaussian kernel at R = 0.7,
then 0.35), whose torque no finite mixture represents exactly, the cohesion leak settles on seven or eight classes inside the saturation
(|a|/a_max ≤ 0.84) with the fit error 0.14%–0.35% after 80000 slow units — below the four-kernel designs of Proposition D‴.29 (1.46% inside,
0.64% at the saturation), above the best equal-weight six- and eight-kernel fits (9.4·10⁻⁶ and 3.3·10⁻⁷, both saturated or nearly): the
count is the transit's, not the best configuration of that count, and its smallest classes are single particles. (G) The resolution (the
fourth bump, m = 30 from the jittered lift of ϑ^b(q) with the counts (10, 10, 10), the compact kernel at R = 1, σ halved every 10⁵ slow
units from σ₀ to σ₀/64 with a jitter 10⁻³ renewed at each stage, three seeds). At w = q, (iv) puts the first split in the middle class at σ
= 3λ⁺₂ = 0.121σ₀: down to σ₀/8 the three classes stay, the fit error 0.0170, the middle class's margin there −2.0·10⁻⁶ to +9.3·10⁻⁷
(8.7·10⁻⁶ without the jitter, (A)); at σ₀/16 it is −1.24·10⁻⁴ to −1.40·10⁻⁴ and the class splits in every seed. Over the 21 stages no class
with a positive margin split, every class of two or more particles whose margin was below −5·10⁻⁵ split within its stage (8 of 8), and of
those between −5·10⁻⁵ and 0 — where the jitter grows by e^{|margin|·10⁵} = 1.2–140 over a stage — 3 of 7; a class counts as split when a
piece passes the clustering threshold 0.25 from the jitter 10⁻³, a growth by 250 within 10⁵ slow units, i.e. a margin below −ln(250)/10⁵ =
−5.5·10⁻⁵, which is where the tallies change. The split is not a halving: the unstable class sheds one to three single particles — (iv)'s
degeneracy leaves the pattern to the nonlinear terms, and the kernel's concavity makes shedding the cheapest: for a class of ten split by 'k
particles at +d, 10 − k at −kd/(10 − k)' at a fixed Σ‖u_i‖² = 0.5 the cohesion energy is 0.43, 0.62, 0.70, 0.74 of its quadratic part for k
= 1, 2, 3, 5 (compact; Gaussian 0.87, 0.93, 0.94, 0.95), the same for all k at 10⁻⁴. A particle of mass 1/30 carries at most a_max/30 = 0.13
of amplitude, less than the fourth bump's 0.3, so the singletons do not represent it; the middle and third classes move toward it (centres
0.19–0.22 and 1.07–1.16 at σ₀/64) and the singletons fill the flanks they leave (centres −0.29 to 0.02 and 1.69–1.76, |a|/a_max up to 0.91).
The fit error falls — 0.0151–0.0152 at σ₀/16, 0.0141 at σ₀/32, 0.0028–0.0031 at σ₀/64 — but the count overshoots: 8–10 classes at σ₀/64,
five to seven of them single particles, for a target of four bumps. σ₀ sets the resolution as (iv) says; which classes form, and with which
masses, the shedding decides — the transit, which is not proved. (H) Other targets (the annealed Gaussian schedule of (C), broad starts, m =
60, three seeds each): two bumps (Σ_k|A_k|/a_max = 0.45) end at two classes with the fit error 8.9·10⁻⁹–1.6·10⁻⁸ in all three runs; the four-bump f_ε of (A) and (G) at three classes with the fit error 0.0174 in all three, near the three-atom approximant's 0.0170 — the fourth
bump lies below the resolution σ₀, and its width is the bound s_min = 0.3, so f_ε has no exact representation; five bumps (s_min = 0.2;
|A_k|/a_max = 0.225, 0.15, 0.175, 0.125, 0.2, Σ = 0.875) at 6–9 classes with the fit error 0.045–0.118. An exact representation needs w_k >
|A_k|/a_max for every bump, which leaves 0.125 of the mass free, and the transit's masses miss it: in every run two or three of the three
inner bumps get less mass than they need (0.067–0.167 against 0.125–0.175), their classes saturated (|a|/a_max 0.96–1.00), and the outer
bumps are carried by pairs of classes (one pair in seed 0, two in seeds 1 and 2). (K) Where the count fails and where it holds (the same
schedule): three equal bumps — three classes in 3 of 3, the fit error 9.4·10⁻⁵–1.8·10⁻⁴, the Gaussian kernel's bias between same-sign
classes 1.5 apart in the centre: (L) it stays at 9.39·10⁻⁵ for 2·10⁵ more slow units and falls to 2.5·10⁻¹³ under the compact kernel; four
alternating bumps (Σ_k|A_k|/a_max = 0.46) — four classes in 3 of 3, 9.6·10⁻³–1.1·10⁻², a slow transient, not a bias: (L) the fit error falls
geometrically, to 6.1·10⁻⁴ after 2·10⁵ more slow units (3.4·10⁻⁴ after the compact stage), the classes' weighted amplitudes the bumps' to
10⁻⁴; three heavy bumps (Σ_k|A_k|/a_max = 0.95, 0.05 of the mass free) — 3–5 classes, 0.058–0.20, the masses short again; two bumps closer
than the broad scale (±0.45, the widths 0.4) — the broad stage merges them into one class of width 0.67–0.72 that the fine stage keeps (its
margin positive, 6.5·10⁻³–7.1·10⁻³: a stable biased configuration), with one or two stray particles, 1–3 classes, 0.098–0.13. So the count
is the target's when the bumps are farther apart than the broad scale and the transit's masses are representable, and the approximant's
below the resolution; the scale and the masses, the memory of the transit, decide the rest — the masses are Proposition D‴.34's remedy, the
scale is not.

**Corollary D‴.33 (the desired regressor keeps the replicated configuration attracting under measurement noise).** In the setting of
Proposition D‴.30 let the adaptive law's regressor be evaluated at the reference — the desired regressor (DR) or the desired compensation
adaptation law (DR*) of Proposition D‴.31(iv) — and let s_k be the stiffness of class k's separations: σ₀/k_e under the fixed and the
adaptive prior, σ₀w_k/k_e under the cohesion leak (Proposition D‴.32(vi)). Write M_k := G^{1/2}∇²_ϑΦ̃(x_m; ϑ_k)G^{1/2}.

(i) *(the drift)* With the regressor at x_m the law's Jacobian on a separation of class k is ε_mG^{1/2}∇²_ϑΦ̃(x_m; ϑ_k(t))G^{1/2}: the
regressor's own term ⟨E[ε_m²]G^{1/2}∂_x∇²_ϑΦ̃G^{1/2}⟩ of Proposition D‴.30(i), which does not decrease with k_e and carries the measured
regressor's drift, is absent. To second order in ν the separation's mean drift is S_k + C_k: S_k = ⟨E[e]M_k⟩ (at the lift E[e] is MR's,
the regressor entering only the parameters' equation) and the cross term C_k = ⟨G^{1/2}Σ_lE[ε_mδϑ_kl]∂_{ϑ_l}∇²_ϑΦ̃G^{1/2}⟩, zero at the
lift, which builds up within a period as the class atom fluctuates — negligible against the measured regressor's drift, of the order of
DR's and DR*'s. A separation's exponent along a unit θ adds the nonnegative variance term E‖P_⊥v₁(t)‖²/t, v₁(t) = ∫₀ᵗε_mM_kθ and
P_⊥ = I − θθᵀ (to second order E log‖u(t)‖/‖u(0)‖ = −s_kt + ∫₀ᵗθᵀ(S_k + C_k)θ + E‖P_⊥v₁(t)‖², ε_mM_k being symmetric).

(ii) *(a threshold bounded uniformly in k_e)* At the lift k_eE[e] = g, |g| ≤ ½|∂²_xf|Eν²(1 + O(|a|/k_e)) (Proposition D‴.31(iv)(c)), so
‖S_k‖ ≤ (Eν²/(2k_e))‖∂²_xf‖_∞⟨‖M_k‖⟩(1 + O(1/k_e)): the drift and the stiffness are both O(1/k_e), and the noise level up to which the
replicated configuration attracts at the lift is bounded below uniformly in k_e — Eν² ≥ 2σ₀/(‖∂²_xf‖_∞max_k⟨‖M_k‖⟩) up to (1 + O(1/k_e)) —,
and constant when the noise's bandwidth is scaled with the gain (B/k_e fixed, where k_eS_k depends on k_e only through a/k_e). For a sensor
of fixed bandwidth the factor (B − k_e − 2a)/(k_e + B) of (iv) changes sign at k_e ≈ B, and DR's drift with it; the thresholds then fall
toward the bound's level. Under the measured regressor the threshold falls like 1/k_e in the variance while k_e ≪ B, and levels off for
k_e ≫ B, where E[ε_m²] → Eν²B/k_e.

(iii) *(the stationary point)* At the adaptive prior's stationary point the slow force balance carries the parameters' part of the
second-order force (Proposition D‴.31(iv)(c)): ⟨Ĵᵀ(g + g_par − Ĵȳ)⟩ + k_e⟨f_{y,par}⟩ = 0, g_par and f_{y,par} the parameters' parts of the
error equation's and of the learning law's forcing, so (averaged)

    k_eE[e] = (I − P)(g + g_par) − ĴH⁻¹k_e⟨f_{y,par}⟩,

P the projection onto the regressor's span: the residual of the noise's own target outside the tangent space, plus the parameters' part,
whose learning-law term lies in the tangent space and is not removed. The drift there is ⟨E[e]M_k⟩ + C_k with the stationary covariance's
cross term. Both parameters' terms are of the order of the noise's residual and need not be smaller: E[e] does not vanish at the stationary
point when g lies in the tangent space.

(iv) *(the sign under DR)* With the compensation f̂(y) the mean error at the lift is the compensation's rectification of the noise, E[e] =
−(∂²_xf/(2k_e))Eν²(B − k_e − 2a)/(k_e + B). A Gaussian kernel depends on x and its centre c through x − c, so ∂²_cΦ̃ = ∂²_xΦ̃, and the
drift's entry along the centre coordinate γ, c = c(γ), with the metric weight G_c := G_γc′(γ)², is −(Eν²G_c/(2k_e))⟨∂²_xf·∂²_xΦ̃_k·(B − k_e
− 2a)/(k_e + B)⟩ − c″(γ)G_γ⟨E[e]∂_xΦ̃_k⟩; for B ≫ k_e its first term is −(Eν²G_c/(2k_e))[w_k⟨(∂²_xΦ̃_k)²⟩_ρ +
Σ_{l≠k}w_l⟨∂²_xΦ̃_l∂²_xΦ̃_k⟩_ρ] (f = Σ_lw_lΦ̃_l at the exact representation, Φ̃_l := Φ̃(·; ϑ_l)), negative when the kernels are separated,
while the saturation's term, zero at c = 0, has no sign. So at the lift a noise fast against the loop (B > k_e + 2 max|a|) stiffens the
centres of separated kernels where that term is dominated, and a slower one (B < k_e − 2 max|a|) reverses the factor and the sign. The other
directions have no sign: the whole matrix is negative definite for separated kernels of alternating sign well inside the orbit, as on the
instance, and indefinite for same-sign, wide or edge kernels (their widths and amplitudes soften). Along the trajectory E[e] relaxes with
the collective modes (k_eE[e] = g − Ĵȳ(t)), the fast modes first, so it passes through the part of g carried by the slow modes, which has no
sign relation to ∂²_xf: the sign holds near the lift, not along the trajectory. With the compensation at the reference (DR*) the mean error
is the plant's rectification of its state jitter, +½∂²_xf·Var(e)/(k_e − a), and the sign at the lift is reversed.

*Proof.* (i) With the regressor at x_m the law is G∇_ϑ(ε_mΦ̃(x_m; ϑ_i)) minus the leak, and its Jacobian on a separation is
ε_mG^{1/2}∇²_ϑΦ̃(x_m; ϑ_k)G^{1/2}: the measurement enters only through ε_m, whose mean is E[e] (Eν = 0), and the class atom's fluctuation
through the Hessian, whose correlation with ε_m is the cross term. A separation does not enter f̂ at first order, so neither the
compensation nor the leak adds a noise term. The exponent: with A₁ = ε_m^{(1)}M_k, w = u e^{s_kt}, log‖w‖ = ⟨θ, w₁⟩ + ⟨θ, w₂⟩ + ½‖w₁‖² −
⟨θ, w₁⟩² + O(ν³), and E⟨θ, w₂⟩ = ∫θᵀE[A₂]θ + ½E‖w₁‖² for a symmetric A₁, whence the variance term. (ii) The mean of the error equation at
second order gives k_eE[e] = g − Ĵȳ (Proposition D‴.31's proof), at the lift ȳ = 0, and the bound on g is Proposition D‴.31(iv)(c)'s;
with B = βk_e the factor of (iv) is (β − 1 − 2a/k_e)/(1 + β). For k_e ≫ B, E[ε_m²] = Eν²(a² + k_eB)/(k_e(k_e + B)) → Eν²B/k_e. (iii) The
mean of the learning law adds the parameters' part k_e⟨f_{y,par}⟩ (the law's ∇²_ϑΦ̃E[ε_mδϑ]) to the least-squares force, and the error
equation's mean forcing is g + g_par (Proposition D‴.31(iv)(c)); stationarity gives Hȳ = ⟨Ĵᵀ(g + g_par)⟩ + k_e⟨f_{y,par}⟩, and Ĵȳ = P(g +
g_par) + ĴH⁻¹k_e⟨f_{y,par}⟩. (iv) Φ̃ depends on x and c through x − c, so ∂_cΦ̃ = −∂_xΦ̃ and ∂²_cΦ̃ = ∂²_xΦ̃; for c = c(γ),
∂²_γΦ̃ = c′²∂²_cΦ̃ + c″∂_cΦ̃, and the quasi-static E[e] of Proposition D‴.30's proof with ∂²_xf̂ = ∂²_xf gives the entry; DR*'s E[e] is in
Proposition D‴.31's proof. ∎

*Status.* (i) and (iv)'s formula exact consequences of the second-order expansion — formal at the order of Proposition D‴.30(i)'s drift,
with its caveats —; (ii)'s bound proved given (i) and D‴.31(iv)(c)'s bound on g; (iii)'s balance formal (averaged), its numbers from the
exact periodic second-order mean; (iv)'s sign proved for the centre direction of separated Gaussian kernels where the saturation's term is
dominated, at the lift; the whole matrix's sign numerical and specific to the target, and lost along the trajectory. Classical: the desired
compensation adaptation law's insensitivity to measurement noise ([S53]); the attraction of identical uncoupled systems driven by a common
noise — noise-induced synchronisation, decided by the transverse Lyapunov exponent at second order in the noise, which strong noise can make
positive (Teramae–Tanaka, Goldobin–Pikovsky, Nakao–Arai–Kawamura, [S61]) — is the frame of Proposition D‴.30(i)'s threshold and of this
corollary; the second-order drift of a regressor without the measured state is the one-line consequence of a y-independent regressor. New:
that consequence for the K₀-atom identification — the separations' drift O(1/k_e) like their stiffness, a threshold bounded uniformly in
k_e, which makes the desired regressor the natural partner of the cohesion leak, whose separations are the softest (σ₀w_k) — and the
instance's quantification: the sign at the lift, its loss along the trajectory, the stationary point. Checked: `tests/test_eiv_remedy.py`
(DR's drift MR's mean-error term alone at the lift, negative definite in every class, MR's 0.30 in class 3; DR*'s O(1/k_e), k_e times its
largest eigenvalue the same at k_e = 25 and 50, and DR*'s drift against an independent quadrature to 10⁻⁴; at the adaptive prior's
stationary point the exact periodic mean's noise and parameters' parts summing to it, the parameters' part dominant and largely in the
tangent space, the drift indefinite; the period maps continuing `second_order_transient` to 10⁻⁷ and DR's drift turning indefinite).

*Numerics* (`exp_eiv_remedy.py` Part D → `results/chaos/eiv_remedy.md`). Proposition D‴.30's instance (k_e = 25, the OU noise B = 100), per
unit Eν², real time. (D1) At the lift the separations' drift has the extreme eigenvalues [−3.7·10⁻¹, +1.5·10⁻²]; [−2.6·10⁻², +1.0·10⁻²];
[−3.5·10⁻², +3.0·10⁻¹] in the three classes under MR, [−1.1·10⁻², −5.0·10⁻⁴]; [−1.5·10⁻², −1.5·10⁻⁴]; [−7.2·10⁻³, −3.4·10⁻⁴] under DR and
[+1.6·10⁻⁴, +3.8·10⁻³]; [+5.1·10⁻⁵, +5.0·10⁻³]; [+1.1·10⁻⁴, +2.4·10⁻³] under DR*: DR's is negative definite, DR*'s largest eigenvalue,
+4.98·10⁻³, is 60 times below MR's 0.2982 — the threshold std 0.238 under the adaptive prior against MR's 0.031, and 0.138 against 0.018
under the cohesion leak's stiffness σ₀w_k (w_k = 1/3). (D2) With B = 4k_e MR's largest eigenvalue is 0.2982 → 0.3045 over k_e = 25 → 800,
its threshold std 0.031 → 0.0054 (∝ k_e^{−1/2}), while k_e times DR's is −0.0038 and k_e times DR*'s 0.1245–0.1242 at every k_e (DR*'s
threshold 0.238–0.239). At the fixed B = 100 DR's drift is indefinite at k_e = 100 and positive definite from k_e = 200, its threshold std
0.185, 0.138, 0.121 at k_e = 200, 400, 800; DR*'s falls 0.238 → 0.113 and MR's 0.031 → 0.014 (its variance by 4.6 over a factor 32 in
k_e). ‖S_k‖ is at most 0.37 of (ii)'s bound in these cases (DR 0.253, DR* 0.084 at B = 4k_e), which gives the threshold the floor std ≥ 0.069
for DR and DR*. (D3) For other exactly represented targets at the lift DR's centre entries are negative for all six; its whole drift is
negative definite for the instance and for two opposite bumps, and indefinite for three separated same-sign bumps (the largest eigenvalue
+2.9·10⁻⁴), wide bumps (S = 1; +5.6·10⁻⁴), bumps at the orbit's edge (+2.7·10⁻⁴) and two close same-sign bumps (+5.0·10⁻⁴). (D4) At the
adaptive prior's stationary point (the exact periodic second-order mean) the mean error's rms is 6.70·10⁻⁴ under DR and 6.44·10⁻⁴ under
DR* — 2.4% and 7.0% of the lift's —, its parameters' part 6.32·10⁻⁴ and 6.27·10⁻⁴ and its noise part 2.11·10⁻⁴ and 1.37·10⁻⁴ (the round-22
averaged noise-only balance (I − P)g/k_e gave 2.34·10⁻⁴ and 2.60·10⁻⁴); its tangent part, 5.80·10⁻⁴ and 5.81·10⁻⁴, is the parameters'
part's. The drift there is indefinite in every class — DR's [−4.35, +3.68]·10⁻⁵; [−0.13, +0.30]·10⁻⁵; [−2.82, +2.44]·10⁻⁵ —, the parameters'
part of E[e] contributing up to +2.6·10⁻⁴ and the cross term down to −2.6·10⁻⁴ in class 2, nearly cancelling. Under the fixed prior the
mean error keeps 41% of the lift's (rms 1.13·10⁻², DR), DR's drift stays negative definite (the largest eigenvalue −3.8·10⁻⁵) and DR*'s
positive definite (+1.50·10⁻³, the threshold std 0.435). (D5) Along DR's exact second-order transient from the lift (82 000 periods,
5.2·10⁵ time units) the drift's largest eigenvalue, −4.7·10⁻⁴, −1.6·10⁻⁴, −3.1·10⁻⁴ in the first period, turns positive at t = 4103, 2947,
2714 in classes 1, 2, 3 — with the noise part of E[e] alone at t = 4072, 2721, 2689: (iv)'s mechanism, not the parameters' terms — and stays
positive, reaching +1.05·10⁻⁴ (class 1, t ≈ 3.1·10⁴), formally the threshold std 1.64, far outside the expansion, before settling at the
stationary point's (the period maps' limit within 0.3% of the exact periodic mean's). DR*'s largest eigenvalue is positive throughout
(3.6·10⁻³ in the first period, 3.4·10⁻⁵ at the stationary point). Over the first 300 time units the cross term moves DR's class-2 window
drift by −8% ([0, 50): −1.61·10⁻⁴, without it −1.49·10⁻⁴) and DR*'s by −17% (+3.69·10⁻³ against +4.45·10⁻³), so DR*'s threshold std over
[0, 50) is 0.277 (0.238 at the lift's quasi-static moments). (D6) The real loop (the lift m = 15, one class separated by 10⁻⁶ along a
direction θ, four antithetic pairs per family on paths common to the laws, 300 time units, the rate fitted over t ≥ 50): along MR's top
eigenvector u of S₃ DR's separations decay at −2.9651·10⁻⁴ ± 1.3·10⁻⁷ (std 0.05) and −3.3650·10⁻⁴ ± 5.3·10⁻⁷ (std 0.1), faster than
σ₀/k_e = 2.8298·10⁻⁴; their noise-induced part is 0.98 ± 0.01 and 0.97 ± 0.01 of the transient's prediction — the mean error's term along
the transient, the cross term and the variance term — (2.3 and 3.5 sampling errors), against 0.82 and 0.81 of the lift's. Along DR's own
top eigenvectors, its weakest directions, the separations decay faster than σ₀/k_e by only 0.36%, 0.14%, 0.26% (classes 1–3, std 0.05;
1.4%, 0.54%, 1.0% at 0.1), at the transient's prediction to 0.9, 1.9 and 0.05 sampling errors (0.8, 1.8, 0.1); along DR*'s top
eigenvector of S₂ they decay slower than σ₀/k_e by 3.2% and 12.6%, the noise-induced part 0.97 ± 0.10 and 0.96 ± 0.10 of the transient's
prediction (0.72 of the lift's). The pair means of DR decay faster than σ₀/k_e in every family; no DR or DR* run grows. MR's separations
grow at +4.69·10⁻⁴ ± 1.3·10⁻⁵ and +1.25·10⁻³ ± 3.1·10⁻⁵ (8 of 8 runs each), the latter less than half the second-order prediction (the next
order of Proposition D‴.30). On the same noise paths the measured regressor's replicas separate and the desired regressor's stay together.

**Proposition D‴.34 (mass transfer: the replicator weight law removes the saturation trap of the transit's masses).** Let the particles of
Proposition D‴.32 carry weights β_i > 0 with Σ_iβ_i = m — the network f̂ = (1/m)Σ_iβ_iΦ̃(·; ϑ_i), its measure μ_{Θ,β} = (1/m)Σ_iβ_iδ_{ϑ_i},
the class masses W_k = Σ_{i∈k}β_i/m — and let the weights adapt by the replicator law

    dβ_i/ds = γβ_i(u_i − λ),   u_i := ⟨r Φ̃(·; ϑ_i)⟩_ρ,   λ := Σ_jβ_ju_j/Σ_jβ_j,                                        (𝔉)

beside dϑ_i/ds = β_iG⟨∇_ϑ(r Φ̃(·; ϑ_i))⟩_ρ − (σ₀/m)Σ_jψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²)(ϑ_i − ϑ_j) (r = f − f̂; in the real loop, on Σβ = m,
β̇_i = γβ_ieᵀ(Φ̃(x; ϑ_i) − f̂(x)), the replicator dynamics with the fitness eᵀΦ̃_i — an exponentiated-gradient law).

(i) *(a gradient flow; the bound kept)* The pair is the gradient flow of J(Θ, β) := mL̄(μ_{Θ,β}) + (σ₀m/2)∫∫W dμ_Θdμ_Θ (the cohesion energy
of the unweighted cloud, so that the weights see the fit only) in the metric Σ_i‖δϑ_i‖²_{G⁻¹} + Σ_iδβ_i²/(γβ_i) on {Σ_iβ_i = m} — Euclidean
on the positions, Fisher–Rao (Shahshahani) on the weights; with ϱ_i := √β_i, gradient descent on (ϱ, Θ), with the step γ/4 in ϱ, on the
sphere Σ_iϱ_i² = m, ϱ_i²/m the output weights of a two-layer network —, so J decreases; β > 0 and Σβ are invariant (from every start), and
on Σβ = m, |f̂| ≤ (1/m)Σβ_i|a_i| ≤ a_max as for equal weights: no projection is needed. With β ≡ 1 and γ = 0 it is the loop of Proposition
D‴.32. (A Wasserstein–Fisher–Rao metric would weigh the positions by their masses, Σ_iβ_i‖δϑ_i‖²_{G⁻¹}: its drift lacks the factor β_i on the
fit and divides the cohesion force by β_i; and J is not a function of μ_{Θ,β} alone, its cohesion being the unweighted cloud's.)

(ii) *(the real loop at fixed positions)* With the positions at atoms where f(x) = (1/m)Σ_iβ*_iΦ̃(x; ϑ_i) at every state the loop visits
(for n = 1 and analytic features, on the orbit's interval), Σβ* = m, the real loop with (𝔉) has V̇ = −k_e|e|² for V := ½|e|² +
(1/(γm))Σ_iβ*_i log(β*_i/β_i): the relative entropy of the weights takes the place of the quadratic parameter error — the negative-entropy
instance of the mirror-descent adaptation laws with Bregman Lyapunov functions ([S62]). The identity needs only Σβ* = m and the pointwise
representation; V ≥ 0 needs Σβ = m, which the law keeps (it does not restore it). β → β* is not claimed: it needs the persistent
excitation of (Φ̃(x_m(t); ϑ_i))_i.

(iii) *(the stationary points; the saturated transients)* (Θ, β) is stationary iff β_iG⟨r∇_ϑΦ̃(·; ϑ_i)⟩_ρ = (σ₀/m)Σ_jψ_ij(ϑ_i − ϑ_j) for
every i and u_i = λ on the support: the fit's force carries the factor β_i, the cohesion does not. At a clustered configuration under the
finite range (the classes R apart) the cohesion vanishes, the condition is ⟨r∇_ϑΦ̃(·; ϑ_k)⟩ = 0, whose amplitude component (1 −
a_k²/a_max²)⟨rφ_k⟩ = 0 forces ⟨rφ_k⟩ = 0, so u_k = a_k⟨rφ_k⟩ = 0 = λ: there the stationary points are Proposition D‴.32's, for every vector
of masses. Under the Gaussian kernel the inter-class forces do not vanish and the fit's factor β_i weighs against them, so u_k = λ fixes the
masses: the replicator selects the masses of each class configuration and moves the Gaussian kernel's bias (Numerics (G)); and a class
whose particles carry unequal weights is not exactly replicated at a biased stationary point. The difference is the saturated transient:
when a class is short of mass (|A_k| > W_ka_max) the equal-weight loop at fixed positions stays biased (it has no critical point) while
|a_k| → a_max with ⟨rφ_k⟩ ↛ 0, its amplitude gradient exponentially small in α_k (heuristically α_k grows like (a_max/2)log s); under (𝔉) the
starved class's mass grows at the rate γ(u_k − λ), u_k ≈ a_max|⟨rφ_k⟩| — dW_k/ds = γW_k(1 − W_k)u_k > 0 exactly when the class's particles
coincide and u_l = 0 for the other classes — and the others take up the gauge (Proposition D‴.28: their p_l = W_la_l kept by |a_l|
increasing).

(iv) *(fixed positions: representability)* With one particle at each of K distinct atoms whose centres and widths are fixed and whose
features φ_k are linearly independent in L²_ρ, and f = Σ_kA_kφ_k in their span (otherwise A is the least-squares coefficient vector), the
flow of (W, α) is a gradient flow of the least-squares loss in p_k := W_ka(α_k), convex in p, on the open simplex times ℝ^K; its critical
points are the exact fits p = A, which exist iff Σ_k|A_k|/a_max < 1 (the representable set of Proposition D‴.24 is not empty); every
trajectory bounded in (log W, α) converges to one (L̄ is real-analytic there: Łojasiewicz), so bounded trajectories reach the exact fit iff
the target is representable (boundedness is not proved); when Σ_k|A_k|/a_max ≥ 1 every trajectory leaves every compact set, each class
saturating or losing its mass, and the fit error stays above the least-squares error over the ℓ¹ ball Σ_k|p_k| ≤ a_max. The transit's
masses no longer decide whether the target is fitted.

(v) *(the linearisation at an exact representation)* At the lift of an exact representation (β ≡ 1) the weights add, on {Σβ = m}, Σ_k(n_k −
1) neutral directions (the redistribution of mass inside a class, which f̂ does not see) and K − 1 gauge directions (mass moved between
classes against their amplitudes, W_ka_k fixed), and no unstable one — in the coordinates (Θ, log β), the constraint not imposed, the common
rescaling of the weights against the amplitudes is a K-th gauge; the separations' block is Proposition D‴.32's (−σ₀w_k, w_k the particle
fraction), the collective block changes by the replicator's metric on the masses. In the real loop the lift is a periodic orbit in a
manifold of periodic orbits (the gauge and the redistributions), neutral along it.

*Proof.* (i) ∂_{ϑ_i}(mL̄) = −β_i⟨r∇_ϑΦ̃(·; ϑ_i)⟩_ρ and ∂_{β_i}(mL̄) = −u_i (f̂ is linear in β); the cohesion energy does not depend on β;
the Fisher–Rao gradient on {Σβ = const} is −γβ_i(∂_{β_i}J − λ′), λ′ = Σ_jβ_j∂_{β_j}J/Σ_jβ_j, which is (𝔉), and dJ/ds = −Σ_i‖dϑ_i/ds‖²_{G⁻¹}
− Σ_i(dβ_i/ds)²/(γβ_i) ≤ 0; d log β_i/ds is bounded on bounded sets, so β stays positive, and Σ_idβ_i/ds = γΣ_iβ_i(u_i − λ) = 0. With β_i =
ϱ_i², δβ_i²/(γβ_i) = 4δϱ_i²/γ. (ii) With Σβ = m, λ(t) = eᵀf̂(x); d/dt Σ_iβ*_i log(β*_i/β_i) = −γeᵀΣ_iβ*_i(Φ̃_i − f̂) = −γm eᵀ(f − f̂) (f =
(1/m)Σβ*_iΦ̃_i, Σβ* = m), and ė = f − f̂ − k_ee; Σ_iβ*_i log(β*_i/β_i) = m·KL(β*/m ‖ β/m) ≥ 0 when Σβ = m. (iii) The stationarity of the
positions and of the weights read off (i); at a clustered configuration under the finite range ψ = 0 across classes and ϑ_i − ϑ_j = 0 inside
them; the amplitude's derivative da/dα = 1 − a²/a_max² = 4e^{−2|α|/a_max}(1 + o(1)); dW_k/ds = (1/m)Σ_{i∈k}γβ_i(u_k − λ) = γW_k(u_k − λ)
and λ = Σ_lW_lu_l = W_ku_k when u_l = 0 for l ≠ k. (iv) In the variables (W, α) the flow is the restriction of (i)'s gradient flow; L̄
depends on them through p, ∂L̄/∂p_k = −⟨rφ_k⟩ =: −g_k, and at a critical point with W in the open simplex the amplitude component gives g_k =
0 for every k, i.e. p is the least-squares solution A; |a(α)| < a_max and ΣW_k = 1 make p = A attainable iff Σ_k|A_k|/a_max < 1. The flow is
the gradient flow of a real-analytic function for an analytic metric in (log W, α), so a bounded trajectory has a single limit, a critical
point (Łojasiewicz's theorem). The dissipation ∫‖∇J‖² ds ≤ J(0) makes every ω-limit point critical, so without a critical point a
trajectory has no ω-limit point and eventually leaves every compact set; and Σ_k|p_k| = Σ_kW_k|a_k| < a_max always, which bounds the fit
error below by the ℓ¹ ball's. (v) At r = 0 the Hessian of mL̄ is the Gram matrix of the derivatives of f̂, and J's Hessian is the
linearisation's up to the metric: inside a class ∂f̂/∂β_i = Φ̃_k/m is the same for every particle (the redistributions), and ∂f̂/∂β_i =
(a_k/a′_k)∂f̂/∂α_i at β_i = 1 (the gauge, one per class; K − 1 combinations of them keep Σδβ = 0); the cohesion energy acts on the
separations only and does not see the weights; J's Hessian is positive semidefinite there (the Gram matrix plus the cohesion's ∇²W(0) = G⁻¹
on the separations), so no direction is unstable; along the lift e ≡ 0, so dβ/dt = 0. ∎

*Status.* (i), (ii), (iii)'s stationarity at clustered configurations under the finite range and (v)'s zero modes proved; under the
Gaussian kernel the masses' selection numerical (G); (iii)'s transient: the rate exact under its conditions, the bias's duration heuristic;
(iv) proved for bounded trajectories, boundedness not proved; the transit with the cohesion leak numerical, and with it any selection of
the number of kernels: the replicator does not act on the separations at a critical point (u_k = 0 = λ there), so the shedding of
Proposition D‴.32(iv) is untouched (Numerics (C)). Classical: Wasserstein–Fisher–Rao (Hellinger–Kantorovich) distances and gradient flows
of measures, conic particle gradient descent and the lift to (√β, Θ) with its global convergence for sparse optimisation on measures,
neuron birth–death, the exponentiated-gradient law and EM's multiplicative update of mixture weights, mirror-descent adaptation laws with
Bregman Lyapunov functions — of which (ii) is the negative-entropy instance —, the replicator equation as the Shahshahani gradient flow with
its relative-entropy Lyapunov function ([S62]). New: the replicator as the mass transfer of the cohesion leak in mean-field adaptive control
— the stationary set unchanged at clustered configurations under the finite range, the masses selected under the Gaussian kernel, the
saturation trap of the transit's masses removed —, the representability criterion at fixed positions for bounded trajectories, and the
numerics.

*Numerics* (`exp_mass_transfer.py` → `results/chaos/mass_transfer.md`; Proposition D‴.32's instances, the midpoint rule on (Θ, log β)).
(A) Fixed positions — the five-bump target's atoms, one particle per class, the masses (0.15, 0.30, 0.20, 0.25, 0.10), two of them short
of |A_k|/a_max = (0.225, 0.150, 0.175, 0.125, 0.200), the amplitudes and the weights adapted: without the replicator the fit error stays at
0.388 from s = 4000 on — the least fit error at those masses is 0.387 (G) —, the two short classes at |a|/a_max = 0.999; with γ = 10⁻³ it
falls to 2.9·10⁻⁴ at s = 40000, with γ = 10⁻² to 3.0·10⁻¹⁴ (geometrically, by about 6 per 2000 slow units), every mass above its
|A_k|/a_max at the end (0.234, 0.199, 0.191, 0.168, 0.208). With the amplitudes scaled by 1.25 (Σ_k|A_k|/a_max = 1.094) the fit error falls
to 0.071 at both rates and still creeps down toward the ℓ¹ ball's least error 0.0684 (G), every class at |a|/a_max ≥ 0.992, the masses
already the minimiser's (0.271, 0.164, 0.194, 0.132, 0.239 against 0.270, 0.164, 0.194, 0.132, 0.239): (iv)'s dichotomy. (B) The five-bump
target from broad starts (m = 60, the annealed Gaussian schedule): without the replicator 6–9 classes, the fit error 0.045–0.118
(Proposition D‴.32, Numerics (H)); with γ = 10⁻³ 5–7 classes, 0.0076–0.0139; with γ = 10⁻² five classes in 3 of 3 runs, their centres the
bumps' to 0.01 (−1.600, −0.798, +0.004, +0.798, +1.597 for seed 1), the fit error 6.0·10⁻⁴–1.4·10⁻³, falling at the end by 2–4% per 1000
slow units. Continued (G, seed 1) it levels off — 5.99·10⁻⁴, 3.62·10⁻⁴, 3.01·10⁻⁴, 2.87·10⁻⁴, 2.82·10⁻⁴, 2.80·10⁻⁴ at s = 0, 2, …, 10·10⁴
under the Gaussian kernel, near its class-level stationary point (2.63·10⁻⁴), and 2.85·10⁻⁴ after 4·10⁴ slow units under the compact
kernel, whose class-level stationary point is exact (2.4·10⁻¹⁶): the nearly saturated classes (|a|/a_max up to 0.96) are slow. (D) The
instance (D‴.32 (C)) with γ = 10⁻²: three classes in 3 of 3 runs, the fit error 2.4·10⁻⁷–9.7·10⁻⁷ (without the replicator, the same seeds
0–2, 1.9·10⁻⁷–9.0·10⁻⁷). (C) The resolution's finest stage (the fourth bump, compact R = 1, σ₀/64, from the jittered lift of ϑ^b(q), 10⁵
slow units): with γ = 10⁻² the same 8–10 classes as without, from the same splits at nearly the same times (3 of 51 records differ by one
class; the fit error 0.0055–0.0075 against 0.0059–0.0078): at a critical point of L_w the fitness vanishes (u_k = a_k⟨rφ_k⟩ = 0 = λ), so the
replicator does not act on the separations' instability, and the shedding of D‴.32 (iv) is untouched. (F) The targets where D‴.32's count
failed (its Numerics (K); the annealed schedule, γ = 10⁻²): three heavy bumps (Σ_k|A_k|/a_max = 0.95) end at three classes at the bumps in 3
of 3 runs, the fit error 7.8·10⁻⁸–1.6·10⁻⁷, every mass above its |A_k|/a_max (0.411–0.423, 0.308–0.313, 0.265–0.280 against 0.400, 0.300,
0.250), where without the replicator 3–5 biased classes ended at 0.058–0.20; two bumps closer than the broad scale (±0.45) end as without
it, merged (1–3 classes, 0.068–0.13, one class of mass 0.88–1.00 at the centre): the scale's obstruction is not the masses', and no mass
transfer separates what the broad stage merged. (G) The stationary points at the class level (a run's classes replicated with equal
weights inside each, the conditions of (iii) solved to rounding, R = 1): for B's seed 1 (the counts 6, 6, 20, 10, 18) under the Gaussian
kernel the fit error is 3.85·10⁻⁴ with the masses held and 2.63·10⁻⁴ with u_k = λ, the masses moving (0.1573, 0.1330, 0.2618, 0.2060, 0.2419 →
0.1579, 0.1318, 0.2661, 0.2076, 0.2366) and the classes of each sign then at a common |a|/a_max (0.9578, 0.8464); under the compact kernel
both are exact (2.4·10⁻¹⁶–2.5·10⁻¹⁶), the masses free to 5·10⁻⁴. For D‴.32 (L)'s three equal bumps (28, 16, 16) the Gaussian kernel's
stationary point has 9.39·10⁻⁵ with the masses held and 9.17·10⁻⁵ with them free, the masses moving from (0.467, 0.267, 0.267) toward
(0.567, 0.217, 0.217) (the residual of the conditions there 1.7·10⁻¹⁰: the outer classes approach the saturation), the compact kernel's
exact either way: (iii)'s difference between the kernels. (E) The real loop at the instance's atoms (m = 3, k_e = 25, γ = 2, β(0) = (1.6,
0.7, 0.7), 200 time units): V from 0.0418 to 2.1·10⁻⁷, monotone; V̇ from the vector field equal to −k_ee² to 1.1·10⁻¹⁶ of max k_ee²; Σβ = m
to 2.4·10⁻¹⁴; β at the end within 1.3·10⁻³ of β* = 1. The averaged loop's linearisation in (Θ, log β) at the lift of ϑ^e(w), counts (6, 5,
4), compact R = 1: with γ = 10⁻² 15 eigenvalues below 10⁻⁹ — the 12 redistributions and the 3 gauges of (v) — and the 45 others in
[−0.1834, −3.081·10⁻⁴]; with γ = 0 the 45 are Proposition D‴.32's predicted spectrum to 2.1·10⁻¹¹ (the finite differences' error).

**Proposition D‴.35 (the split of a class is an S_n-symmetric transcritical bifurcation: shedding, its order, its cascade).** In the setting
of Proposition D‴.32 (the finite range), let class k of n ≥ 3 particles sit at an atom ϑ_k of a nondegenerate critical point of L_w, the
classes at least R apart, and let λ⁺ = λ_max(R^G_k) be simple with the whitened unit eigenvector v. Write the separations along v as u_i =
ξ_iG^{1/2}v (Σ_iξ_i = 0), μ := λ⁺ − σ₀w_k (the negative splitting margin), and

    c₃ := ∂³_t⟨r Φ̃(·; ϑ_k + tG^{1/2}v)⟩_ρ |_{t=0}   (r the residual at the critical point, frozen),

the residual's third derivative along the splitting direction.

(i) *(the reduced equation)* For |μ| small against the rates of the other modes, the flow has a centre manifold over the separations along v
— the standard representation of the symmetric group S_n, which permutes the particles — and on it

    dξ_i/ds = μξ_i + (c₃/2)(ξ_i² − (1/n)Σ_jξ_j²) + O(|ξ|³ + |μ||ξ|²):

the quadratic coefficient is the residual's alone, since the cohesion energy is even in the differences (W(−u) = W(u)) and the mean term is
common to the particles. The splitting margin of Proposition D‴.32(iv) is the linear coefficient; the quadratic one identifies the
equivariant normal form's coefficient for this flow. At third order the cohesion's concavity, the residual's own cubic term and the
slaving of the class mean all enter.

(ii) *(no small split is stable)* If c₃ ≠ 0 the bifurcation at μ = 0 is transcritical. Its equilibria near 0 off the S_{n/2} × S_{n/2}
symmetry are the split states with isotropy S_q × S_{n−q} — q particles at (n − q)t*, n − q at −qt*, t* = −2μ/(c₃(n − 2q)), 1 ≤ q < n/2 —
and on the centre manifold each has the eigenvalues −μn/(n − 2q) (multiplicity q − 1, inside the q-group), +μn/(n − 2q) (multiplicity
n − q − 1, inside the other group) and −μ (along the branch), all to leading order in μ: every one is a saddle, for either sign of μ. For
even n the quadratic term vanishes on the S_{n/2} × S_{n/2} states, which are cubic (at |ξ| ~ √|μ|, on one side of the threshold only).
Past the threshold (μ > 0) the class leaves a fixed neighbourhood of the replicated state from almost every initial condition: the split
is a jump, not a gradual separation.

(iii) *(who leaves first)* The particles keep their order along v (Fix of a transposition is invariant under any S_n-equivariant flow, so
two particles never cross), so the particle with the largest sign(c₃)ξ_i is ahead on that side at all times; in the quadratic truncation
its excursion ζ = sign(c₃)ξ obeys dζ/ds = μζ + (|c₃|/2)ζ² up to the common term, with the Riccati solution s = (1/μ)ln[ζ₁(μ +
|c₃|ζ₀/2)/(ζ₀(μ + |c₃|ζ₁/2))] — the leading behaviour of the departure, which is measured against the class's median. The quadratic term
decides once ζ exceeds 2μ/|c₃|: the top-ranked particle then runs away along sign(c₃)v, alone when its lead is large, with its neighbour
as a bound pair when their ζ₀ are nearly tied (within about 11% on the instance). A split sheds the extreme particle, alone or with a
near-tied neighbour, because the quadratic term — the residual's, not the kernel's — makes it run away first.

(iv) *(after a departure)* After k particles have left the kernel's range the class has n − k, and its margin is σ₀(w_k − k/m) − λ⁺(r′, ϑ′),
the curvature at the new residual and the remaining particles' raised amplitude (the gauge of Proposition D‴.28: by n/(n − k) when the
residual is unchanged): both the stiffness's drop and the curvature's gauge rise favour further departures, and the class stops shedding
only if the residual's change lowers λ⁺ by more than μ + kσ₀/m plus that rise. On the instance the residual's change dominates and the
margin is positive by the end of each stage; the cascade resumes when σ is lowered (Numerics (C), (D)).

(v) *(the mass-free leak)* The mass-free leak σ(m_i − ϑ_i) — the kernel-normalised mean shift, σ/(σ₀ψ̄_i) times the cohesion force —
vanishes at clustered configurations under the finite range as the cohesion leak does and has its blocks with σ in place of σw_k on the
separations (the class means free), so its thresholds are σ = λ⁺_k, lower than the cohesion leak's by the factor w_k, and (i)–(iv) hold for
it with μ = λ⁺ − σ. It is not a gradient flow, and its forces do not sum to zero. It does not stop the cascade: the curvature's change after
a departure is the same.

*Proof sketch.* (i) The centre manifold theorem for the S_n-equivariant flow; the expansion of the learning drift along v at the replicated
state (the mean term (S₂/2m)⟨Φ″Φ′⟩ is common to the particles and leaves with the projection on Σξ = 0); W even, so ∇³W(0) = 0. (ii)
Substituting the isotropy ansatz gives μ(n − q)t + (c₃/2)(n − q)(n − 2q)t² = 0; the Jacobian J_ij = (μ + c₃ξ_i)δ_ij − (c₃/n)ξ_j on {Σξ = 0}
has the stated eigenvalues on the zero-sum perturbations inside each group and along the branch; on Fix(S_{n/2} × S_{n/2}) Σξ_i² is
constant across the groups' values and the quadratic term is the common one. (iii) Equivariance under the transposition (i j) makes
{ξ_i = ξ_j} invariant; the Riccati solution. (iv) Proposition D‴.32(iv) at the new configuration with Proposition D‴.28's gauge scaling of
R_k. (v) Inside a class in range m_i is the class mean, so the force is σ(ϑ̄_k − ϑ_i), −σ on the separations, and nothing on the class mean;
across classes ψ = 0. ∎

*Status.* (i)–(ii) proved to leading order in μ (the centre manifold, the S_n-equivariant normal form; checked: the branches by Newton on
the full drift, their eigenvalues converging to the predicted ones as μ → 0); at σ₀/16 three stable modes are slower than μ, so there the
reduced equation's agreement (C) reflects the weak coupling of those modes, not the centre manifold; (iii) the order by equivariance, the
escape time in the quadratic truncation (checked); (iv) a statement of the margin, the cascade's extent numerical; (v) the blocks proved,
its cascade numerical. Classical: S_n-equivariant bifurcation theory — the standard representation's quadratic equivariant makes generic
symmetry-breaking from the fully symmetric state transcritical, with no stable small branch (Elmhirst; Golubitsky–Stewart; Arjevani–Field),
the jump of sympatric speciation (Stewart–Elmhirst–Cohen, Golubitsky–Stewart), the permutation-symmetric critical points of neural
networks (Fukumizu–Amari, Brea et al., Şimşek et al., Arjevani–Field), a cluster's split at a critical parameter in deterministic
annealing (Rose), order preservation (equivariance) and the Riccati escape time. New: the identification for the cohesion leak — the
quadratic coefficient the residual's third derivative along the splitting direction, the cohesion and the mean term entering only at third
order or in the common direction —, the margin after a departure, and the instance's quantification.

*Numerics* (`exp_split_bifurcation.py` → `results/chaos/split_bifurcation.md`; Proposition D‴.32's resolution instance: the lift of ϑ^b(q)
with the counts (10, 10, 10), m = 30, the compact kernel at R = 1). (A) At the lift λ⁺ = 1.00·10⁻⁴, 2.86·10⁻⁴, 1.08·10⁻⁴ for the three
classes, v the centre against the width ((0, +0.75, −0.66), (0, +0.72, −0.70), (0, +0.78, −0.62), whitened), c₃ = 1.42·10⁻⁵, −1.257·10⁻³,
1.20·10⁻⁴ (two difference steps agreeing to 2·10⁻⁵ relative); the middle class splits first, at σ = λ⁺/w_k = 0.1213σ₀, and leaves on the
side −v (its centre lower, its width larger). (B) At μ = 2·10⁻⁶ and 5·10⁻⁷ Newton on the full particle drift (90 coordinates, the drift to
10⁻¹⁷) finds the four split states q = 1, …, 4 from (ii)'s t*: their positions within 1.3%, 1.8%, 3.5%, 15% of the prediction at μ = 2·10⁻⁶
and within 0.3%, 0.4%, 0.8%, 3.0% at 5·10⁻⁷; on the separations along v the eigenvalues are +1.2544μ (×8) and −1.094μ for q = 1 at 2·10⁻⁶
and +1.2511μ (×8) and −1.023μ at 5·10⁻⁷, against +1.25μ and −μ, and likewise for q = 2, 3, 4 (the unstable ones +1.6708μ, +2.5165μ, +5.1444μ
at 5·10⁻⁷ against 5μ/3, 5μ/2, 5μ), the errors falling about fourfold with μ, as the next order says; every split state is a saddle, with n −
q − 1 = 8, 7, 6, 5 unstable directions of the full Jacobian. (C) From the lift with a jitter 10⁻³ (six seeds) at μ = 5·10⁻⁵ and 1.39·10⁻⁴
(σ₀/16, Part G's first split), an excursion being a particle 0.1 from the class's median (far inside the range R = 1): the reduced quadratic
equation integrated from the measured ξ(0) names the first excursion's particle in 12 of 12 runs and gives its time to 1.1% (μ = 5·10⁻⁵:
6.61·10⁴–8.84·10⁴ slow units; σ₀/16: 2.68·10⁴–3.35·10⁴); the linear equation (c₃ = 0) names it in 6 of 12 and misses the time by up to
22.4%, the equation with c₃'s sign flipped names it in 2 of 12 — the ranking, which every equivariant law preserves, does not name the
particle by itself; the quadratic term does. The first excursion is on the side sign(c₃)v in 10 of 12 runs; the two on the other side (σ₀/16, seeds 2 and 4) are the
linear growth of a particle whose excursion stays below 2μ/|c₃| = 0.22, and neither leaves the range. After 10⁵ slow units the pieces
beyond the kernel's range are, at σ₀/16, one particle in 3 of 6 runs and a bound pair in 3 of 6 — the top two of sign(c₃)ξ(0) within 1.2%,
4.4% and 11% of each other; where the lead was 22% or more the leader went alone —, at 1.5–1.9 (single) and 1.1–1.2 (pairs) from the main
piece; at μ = 5·10⁻⁵ one particle in 3 of 6 runs and none yet in the other 3. The main piece at the clustering threshold (0.25) then has 8 or
9 particles and a positive margin (7.5·10⁻⁶–1.6·10⁻⁴): the cascade stops within the stage, and resumes at the next lower σ. (D) Part G's schedule (σ halved every 10⁵ slow units, the jitter renewed, seeds 0–2): the cohesion leak ends at 8–10
classes at σ₀/64 (five to seven single particles, the fit error 0.0028–0.0031, Proposition D‴.32 (G)); with the replicator (γ = 0.1) at
6, 6 and 7, the first escapee's mass growing to 0.105–0.106 (three particles' worth) and the later ones to 0.029–0.060, the fit error
0.0031–0.0032. (E) With the replicator at γ = 1 the same schedule ends at 5, 6 and 8 classes, the first escapee's mass 0.134–0.136 (four
particles' worth), the fit error 0.0029–0.0030: with the replicator the cascade ends at 5–8 classes against 8–10 without, and the rate
makes no difference in three seeds; the jitter's seeds are not paired with Part G's. (F) The mass-free leak σ(m_i − ϑ_i) — the
separations' stiffness σ instead of σw_k — splits later, its threshold λ⁺ = 0.040σ₀ instead of λ⁺/w_k = 0.121σ₀: at σ₀/64, where Part G's
schedule stops, it has 4 classes in 3 of 3 runs (the fit error 0.0149); continued to σ₀/256 it sheds as the cohesion leak does, to 9,
10, 9 classes at σ₀/128 and 9, 8, 9 at σ₀/256 (four to six single particles, a pair in one run), the fit error 0.0028 — the cohesion
leak's end state: its 4 classes were where the schedule stopped, not a count the leak selects.

**Proposition D‴.36 (the gauge under dither: the replicator's masses diffuse along the manifold of exact fits; a leak on the masses pins
them).** Let the loop of Proposition D‴.34 sit at a representable exact representation with one particle per class — K classes, the masses W, the amplitude coordinates α_k with W_ka(α_k) = A_k —, and dither the positions by √(2τG) per particle
(Proposition D‴.32(vi)).

(i) *(the gauge manifold)* The exact fits through the configuration form a manifold M = {W_ka(α_k) = A_k, the centres and widths fixed} of
dimension K − 1 inside the representable set, parametrised by the masses: J is constant on M, and where H_L ≻ 0 the flow attracts M
exponentially, its tangent T spanned by the gauges δα_k = −A_kδW_k/(W_k²a′_k), ΣδW = 0 (a′ = 1 − (a/a_max)²). With n_k > 1 particles in a
class the redistributions of its mass among them are further zero modes (Proposition D‴.34's count), M has dimension m − 1, and T_α
below carries the factor n_k: the statements here are for one particle per class.

(ii) *(the masses diffuse)* As τ → 0, on the time scale 1/τ, the masses perform a diffusion on M (a Katzenberger limit) whose
covariance per unit slow time is

    D(W) = (TᵀgT)⁻¹Tᵀ g Σ g T (TᵀgT)⁻¹,   g = diag(G⁻¹ on the positions, (m/γ)/W_k on the masses),   Σ = 2τG on the positions,

the g-orthogonal projection of the noise onto M (the linearised gradient flow's limit map), plus a drift of the same order in τ from the
second derivative of that map (computed in Proposition D‴.37). For two classes D_11 = 2τT_α/(T_α + T_W)², T_α := Σ_k(A_k/(W_k²a′_k))²/G_α, T_W := (m/γ)Σ_k1/W_k: there is no
restoring force in D; D → 0 as γ → 0 (with equal weights the amplitudes are pinned by the fit) and D → 2τ/T_α as γ → ∞; at the boundary of
the representable set a′_k → 0, so the saturating class's row and column of D vanish quadratically in W_k − |A_k|/a_max (the other
masses keep diffusing), while in the amplitude coordinate α_k ≈ −(a_max/2) log(W_k − |A_k|/a_max) the diffusion tends to 2τG_α exactly:
heuristically the class's excursion into the saturation is unbounded in α_k — outside the theorem, which holds on compact parts of M,
and subject to the drift of Proposition D‴.37. The function f̂
does not move: the gauge's motion changes the representation's conditioning, not the fit.

(iii) *(a leak on the masses)* The entropic leak on the masses — the energy (ε/γ)Σ_iβ_i log β_i added to J, so that d log β_i/ds gains
−ε(log β_i − Σ_jβ_j log β_j/Σ_jβ_j) — is a leak on the unexcited subspace: it does not change f̂ along M, so at an exact representation with
W = w (equal weights) the configuration stays stationary — no bias where w is representable — and the gauge becomes an Ornstein–Uhlenbeck
direction: linearised at W = w the masses' limit on M is an Ornstein–Uhlenbeck process with the stiffness κ = ε(TᵀgT)⁻¹T_W (T_W the
masses' share of the metric on M, Tᵀdiag(0, (m/γ)/W)T), so their stationary covariance C solves κC + Cκᵀ = D — for two classes C =
D/(2κ), κ = εT_W/(T_α + T_W) — and the spread falls as ε^{−1/2}, while the drift term of (ii) shifts the pinned state by O(τ/ε); where w
is not representable the leak pulls the masses back toward the saturation trap that the replicator removed (Proposition D‴.34(iii)), a
bias of order ε/γ — the leak's energy (ε/γ)Σβ log β against the fit's —, not ε.

*Proof.* (i) On M each amplitude is a function of its class's mass, α_k = a⁻¹(A_k/W_k), so M is a graph over the representable masses;
the derivatives of f̂ along the complement of T are those of the fit at fixed masses, whose Gram matrix H_L is ≻ 0 by assumption. (ii) Near an attracting manifold of equilibria of a gradient flow ẋ = −g⁻¹∇J, the limit map Φ (the point of
M the flow converges to) has derivative ∂Φ = the g-orthogonal projection onto T (the linearisation g⁻¹∇²J is g-self-adjoint with kernel T);
Katzenberger's theorem gives the limit dY = ∂Φ(Y)√(2τ)B dW + τ∂²Φ(Y)[BBᵀ]ds, BBᵀ = G on the positions (so 2τBBᵀ = Σ); in the masses' coordinates, the coefficient of δW along
T is (TᵀgT)⁻¹Tᵀgξ. For K = 2 T is one vector with ‖t‖²_g = T_α + T_W and Tᵀgξ = Σ_k t_{α,k}ξ_k/G_α of variance 2τT_α. Near the boundary
a′_k ≈ 2(1 − |a_k|/a_max) ∝ W_k − |A_k|/a_max. (iii) The leak's term is the Fisher–Rao gradient of (ε/γ)Σβ log β on Σβ = m; at β ≡ 1 it
vanishes; along M f̂ is constant, so the leak's energy is the only one that varies there, its Hessian on T being diag(εm/(γW)); the
reduced flow on M in the metric TᵀgT is then −(TᵀgT)⁻¹H(W − w) to first order, and the linear SDE's stationary covariance is the Lyapunov
equation's solution. ∎

*Status.* (i) proved for one particle per class; (ii) the diffusion coefficient by Katzenberger's theorem on compact parts of M (the drift
term in Proposition D‴.37), the vanishing of D at the boundary proved, the excursion in α_k heuristic; (iii) the stationarity, the restored
stiffness and the linearised spread proved, the O(τ/ε) shift and the O(ε/γ) bias sketched. The time-scale separation is thin on the
instance: the slowest normal rate (3.06·10⁻⁴, the second bump's width near its lower bound) is only about 3τ, and at ε = 10⁻³ the leak's
stiffness exceeds it; the reduction holds because that mode barely couples to the masses. Classical: diffusions forced onto a
manifold by a large drift (Katzenberger), the noise-driven motion along manifolds of minima of overparameterised models and its implicit
regularisation (Blanc et al., Li–Wang–Arora) ([S64]). New: the gauge of the amplitude–mass representation (Proposition D‴.28) as such a
manifold for the replicator of Proposition D‴.34, the masses' diffusion coefficient computed in the flow's metric, the excursion into the saturation,
and the entropic leak on the masses as a leak on the unexcited subspace that pins the gauge without bias where the particle fractions are
representable.

*Numerics* (`exp_gauge_dither.py` → `results/chaos/gauge_dither.md`; Proposition D‴.32's three-bump instance with one particle per class at
the exact representation ϑ^e(1/3), m = 3, A = (1, −0.7, 0.8), a_max = 4 — the representable masses W_k ≥ |A_k|/a_max = 0.25, 0.175, 0.2 —,
the replicator at γ, the positions dithered at τ = 10⁻⁴, Euler–Maruyama at dt = 0.25). (A) At W = 1/3 the predicted diffusion of the masses
along M is D_kk = 2.85·10⁻⁸, 1.63·10⁻⁸, 1.81·10⁻⁸ per unit slow time at γ = 10⁻², rising with γ (4.85·10⁻¹⁰ for the first class at
γ = 10⁻³, 4.01·10⁻⁷ at γ = 1, 4.23·10⁻⁷ at γ = 10: saturating as γ → ∞, as for two classes). (B) Measured over 64 seeds at γ = 10⁻² for 5·10⁴ slow units, the
variance of the masses' displacement (the mean square less the drift's square) grows linearly (1.5·10⁻⁴ at 5·10³, 1.4·10⁻³ at 5·10⁴ for
the first class), its slope over [10⁴, 5·10⁴] (2.69 ± 0.58)·10⁻⁸, (2.21 ± 0.47)·10⁻⁸, (1.62 ± 0.38)·10⁻⁸ (bootstrap errors) against the
predicted 2.85, 1.63, 1.81·10⁻⁸: the ratios 0.95, 1.36, 0.90, within two standard errors (D's variation along the path not included). (C) Over 2·10⁵ slow units (10 seeds) the masses spread to standard deviations 0.057, 0.053, 0.052 over the second half (the free
diffusion's √(Ds) at s = 1.5·10⁵: 0.065, 0.049, 0.052, comparable within the 10 seeds' scatter), and the amplitude coordinates run into the saturation:
α_1 up to 16.1 (|a|/a_max = 0.9994), against [3.58, 4.43] and |a|/a_max ≤ 0.803 at γ = 0 (4 seeds, the masses fixed); the fit error over
the second half is the dither's in both (0.0253 with the replicator, 0.0269 without): the gauge's motion costs no fit — it changes only
the representation's conditioning. (D) The entropic leak at ε = 10⁻⁴ and 10⁻³ (4 seeds each, 2·10⁵ slow units) pins the masses: their
second-half standard deviations 0.0152, 0.0103, 0.0117 and 0.0044, 0.0032, 0.0034 against the Ornstein–Uhlenbeck prediction of (iii)
at W = 1/3, 0.0139, 0.0099, 0.0105 and 0.0044, 0.0031, 0.0033; four seeds scatter by 12–23%, so only the ε = 10τ agreement is
informative, and at ε = τ the linearisation is not expected to hold (the shift and the spread's variation along M are comparable) — the
16-seed test is Proposition D‴.37 (G). The amplitudes stay inside the range in these four seeds (|a|/a_max ≤ 0.89 and 0.82, against 0.803
at γ = 0), the maximum growing with the number of seeds; the means (0.3252, 0.3397, 0.3351) and (0.3329, 0.3338, 0.3333) are within
1.5 standard errors of 1/3, too few seeds to resolve the O(τ/ε) shift (Proposition D‴.37 (G)). The fit error is the dither's throughout
(0.0264, 0.0265). The drift itself is visible in (B): after 5·10⁴ slow units the masses' means are 0.3212 ± 0.0046, 0.3445 ± 0.0041,
0.3342 ± 0.0038 (from 1/3), the first two 2.6 and 2.7 standard errors away — the drift Proposition D‴.37 computes, which the diffusion's
√(Ds) (0.038 for the first class) dominates on this horizon.

*Numerics* (`exp_chaos_learning.py` → `results/chaos/learning_table.md`). The instance: n = 1, saturated features with a_max = 4,
c_max = 3, s ∈ [0.3, 2], prior ϑ⁰ = 0, G = diag(1, 0.25, 0.25); a plant of three Gaussian bumps, (A, C, S) = (1, −1.5, 0.5),
(−0.7, 0, 0.4), (0.8, 1.5, 0.5), which three kernels of weights 3A_k represent exactly; reference x_m = 2 sin t, ρ on a grid of 256
times; slow time, midpoint rule with step 0.02; m ∈ {16, 64, 256} sharing their particles with an M = 2048 reference (the stand-in for
the mean-field flow, as in Theorems D′ and D″); horizon 900 slow units; four seeds per start. Proposition K gives σ_crit = 0.04716
(collapsed atom a = 0.331, c = −0.057, s = 1.218); the runs use σ₀ = 0.15σ_crit = 0.00707 (leak time 141 slow units). Whitened
constants by search: L_ϑ^G = 4.044 (closed form 4.658), L₂^G = 11.11 (closed form λ_Ga_max/s_min² = 44.44); K_Φ = 4, F_ρ = 0.999.

*The hypothesis at the stationary measure, and the window (table 1).* The critical point of J₃ reached from the exact representation,
atoms (a, c, s) = (2.4504, −1.4161, 0.8287), (−2.2879, 0.0120, 0.7009), (2.0478, 1.3427, 0.8810), has ‖∇J₃‖_∞ = 4·10⁻¹⁷ (a stationary
measure μ* of the mean-field flow), fit error ‖f − f̂‖/‖f‖ = 0.244 on the orbit (the leak's bias: amplitudes 2.0–2.5 against 3A_k =
2.1–3, widths 0.70–0.88 against 0.4–0.5), and a positive whitened Hessian (eigenvalues 0.0057–0.2003). κ_R at the atoms is 0.00004,
−0.00227, 0.00411; over the spheres of radius ϱ′ around them (400 directions, radii ϱ′ and ϱ′/2) it is 0.00461 (ϱ′ = 0.05), 0.00512
(0.1), 0.00619 (0.2) and 0.00730 (0.3), so (H) holds at μ* with ϱ′ = 0.1 and margin κ = σ₀ − 0.00512 = 0.00195 (relaxation time 1/κ =
512) and fails from ϱ′ ≈ 0.3 on. The window of Remark D‴.1 is realised: κ_R[μ*](N_{0.1}) = 0.00512 < σ₀ = 0.00707 < σ_crit = 0.04716.
At σ₀ = 0.3σ_crit, by contrast, the three-atom critical point (fit error 0.406) fails the splitting condition at its third atom
(κ_R = 0.01696 > σ₀ = 0.01415): there the mean-field limit has more atoms, and the hypothesis has to be checked at that configuration.
The theorem's constants for the warm start (Λ₁ = 0, κ := ½min{κ(ϱ′), min_s a(s)} = 0.00070, the factor ½ leaving room for the sampling of
the supremum by finitely many directions): R_Θ = 2857, ϱ = 1.1·10⁻⁶, θ = 4.6·10⁻¹¹,
C₁ = 2.2·10⁹, m₀ = 1.2·10²³ — the worst-case envelope of Remark D‴.6.

*The hypothesis along the flow (table 2).* On the warm start (μ₀ the mixture of Gaussians of width 0.02 around the atoms) the margin
a(s) = σ₀ − κ_R[μ*](spheres of radius 0.1 around 256 reference particles, 24 directions) stays between +0.00139 and +0.00167 over the
whole run, so (H) holds with s₁ = 0; the same statistic with the reference's own residual is larger by at most 0.00013 (its
sampling fluctuation; the proof's bound (0.2) allows L₂^G‖r_ref − r*‖_{L¹(ρ)} ≤ 0.214), the reference's rms drift decays from 6·10⁻³
to 3·10⁻⁶ and its fit error stays at 0.245–0.248: the reference is stationary. On the cold start (a ~ N(0, 0.8²), c ~ U(−2, 2), s ~
U(0.4, 1), a connected cloud) the margin is negative over the whole run, −0.196 at the first record (−0.179 at s = 1) and −0.0025 at
s = 900, while the fit error
falls from 0.948 to 0.244 — the same fit as the warm limit — and the rms drift from 0.15 to 4.7·10⁻⁴: the cold reference converges
to a configuration of the same quality whose support reaches the ridges (Remark D‴.2).

*Warm start (table 3).* Window maxima of D_m over (0, 112.5], (112.5, 225], (225, 450], (450, 900]:

| m | (0, 112.5] | (112.5, 225] | (225, 450] | (450, 900] | mean on (450, 900] | C₁/m | mean ‖f̂_Θ − f̂_ref‖, last window | max_i‖δ_i‖_{G⁻¹} at s = 900 |
|---|---|---|---|---|---|---|---|---|
| 16 | 0.396 ± 0.24 | 0.423 ± 0.26 | 0.441 ± 0.28 | 0.448 ± 0.28 | 0.446 ± 0.28 | 1.4·10⁸ | 0.0307 ± 0.0097 | 1.36 |
| 64 | 0.122 ± 0.046 | 0.131 ± 0.049 | 0.136 ± 0.050 | 0.138 ± 0.050 | 0.138 ± 0.051 | 3.4·10⁷ | 0.0154 ± 0.0033 | 0.62 |
| 256 | 0.0119 ± 0.0050 | 0.0127 ± 0.0052 | 0.0134 ± 0.0055 | 0.0137 ± 0.0056 | 0.0136 ± 0.0055 | 8.4·10⁶ | 0.00467 ± 0.0012 | 0.21 |

The distance saturates: last-to-first window ratios 1.13, 1.13, 1.15 and last-to-previous 1.02, 1.01, 1.02. Its m-law has slope
−1.28 ± 0.30 (per-seed fits; the theorem's mean square: −1), the learned function's −0.69 ± 0.16 (m^{−1/2}: −0.5). The per-particle
maxima 1.36, 0.62, 0.21 are far above 3ϱ′/4 = 0.075: the localisation of part (c) is not in force at these m (Remark D‴.6), and the
values lie below the vacuous C₁/m by 8.4–8.8 orders of magnitude.

*Proposition D‴.5 (the same table).* The sensitivities of the three atoms to the weights have whitened norms ‖S_k‖ = 8.38, 4.66, 8.32.
The multinomial expectation E D_∞ = (1/m − 1/M)Σ_kw̄_k tr(S̃_kC̃S̃_kᵀ) is 0.823, 0.201, 0.0454 against the four-seed means 0.446 ± 0.28,
0.138 ± 0.051, 0.0136 ± 0.0055 (ratios 0.54, 0.69, 0.30: a mean over four seeds of a heavy-tailed quantity, E|w^m − w^M|² being carried by
the rare large fluctuations). The per-seed law is the sharp test: Σ_kw^m_k‖S̃_k(w^m − w^M)‖² computed from the realised cluster weights
predicts the measured D_m at s = 900 within 25% in 9 of the 12 (seed, m) pairs (for instance 0.123 → 0.112, 0.187 → 0.186, 0.206 →
0.205 at m = 16; 0.0303 → 0.0261, 0.00087 → 0.00085, 0.0081 → 0.0087, 0.0187 → 0.0190 at m = 256); the largest deviation, a ratio 0.57
at m = 16, is the seed with the largest weight fluctuation |w^m − w^M| = 0.275, where the linearisation's remainder and the saturation
a ≤ a_max enter. The exact first equality of (ii) — Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖² with the weighted stationary points themselves —
reproduces the measured D_m within 2.3% in all 12 pairs (1.292 → 1.292 for the seed just mentioned): every m-system sits at the
replicated weighted stationary configuration. The splitting margin of Θ*(w^m) is positive in 11 of the 12 pairs (+0.0006 to +0.0050);
it is negative, −0.0022, for m = 16, seed 3, whose weights (9, 4, 3)/16 lie outside the admissible set: there the replicated
configuration is a saddle of J₁₆ and the plateau at s = 900 is the saddle's — referee round 9 continued that run to s = 3600 and saw
D grow from 1.29 to 5.13 as two particles changed basins and the weights re-sorted to (7, 5, 4)/16 (Proposition D‴.5(iv); the
root-mean-square weight fluctuation is 0.204, 0.102, 0.051 at m = 16, 64, 256 against admissible margins of 0.06–0.36 in ‖·‖_∞). Per
cluster (seed 0), the measured mean displacements of the m-system's particles point along S_k(w^m − w^M) with cosines 0.92–1.000 and
norms within 0.8–20% (m = 16: 0.434 against 0.455, 0.130 against 0.129, 0.338 against 0.358): the particles of a cluster move together,
by the predicted vector. The ratios 0.54, 0.69, 0.30 of the four-seed means to the multinomial expectation are the sampling of a
heavy-tailed quantity at m = 16 and 64 and, at m = 256, four seeds whose realised fluctuations happen to be small.

*Cold start (table 4).* Window maxima of D_m: 1.81 ± 0.28 → 3.13 → 4.50 → 5.33 ± 0.16 (m = 16), 0.505 → 0.900 → 1.42 → 1.79 ± 0.93 (64),
0.0744 → 0.162 → 0.347 → 0.583 ± 0.19 (256): the distance grows over the whole horizon (last-to-first 2.95, 3.54, 7.84; last-to-previous
1.19, 1.26, 1.68), with slope −0.87 ± 0.14 in m (the learned function −0.66 ± 0.06, steeper than the −0.5 of the warm start). The per-particle deviations are bimodal: medians
0.40, 0.14, 0.067 (the m^{−1/2} law of the particles that stayed with their references) against maxima 6.7, 6.3, 6.2 for every m (the
separation between basins), and the fraction of pairs more than one whitened unit apart is 0.203 ± 0.047, 0.082 ± 0.021, 0.037 ± 0.010 —
a factor 2.5 and 2.2 per factor 4 in m, a decay like m^{−0.64±0.08} (per-seed fits; the medians decay like m^{−0.62±0.07}) — steeper than
the m^{−1/2} of Remark D‴.2's boundary-layer count, which is a lower-order prediction (the flips contribute p_m·O(R²) ≈ 0.2·38 = 7.8, 3.1,
1.4 against the measured 5.3, 1.8, 0.58). This is Theorem D″'s negative finding reproduced and diagnosed on the averaged
system: a connected cold start that splits between basins is outside (H), and its coupling distance is dominated by the pairs that
chose different basins.

*The real loop (table 5).* The closed loop at k_e = 25 (σ = σ₀/k_e = 2.8·10⁻⁴, τ = 0, dt = 10⁻³, three seeds) from the same warm start:
window maxima over (0, 500], (500, 1000], (1000, 2000] real-time units are 0.328 ± 0.11 → 0.467 → 0.590 ± 0.17 (m = 16), 0.0605 → 0.0801 →
0.0909 ± 0.045 (64), 0.0188 → 0.0253 → 0.0294 ± 0.011 (256), slope −1.13 ± 0.15 in m; the reference's fit error stays at 0.245–0.246. The
horizon is 80 slow units, a sixth of the relaxation time 1/κ = 512, and the growth (last-to-first 1.80, 1.50, 1.56) is the early
transient of the averaged dynamics: the averaged flow started from the real loop's own initial clouds (reproduced from their seeds)
gives D_real(t)/D_averaged(t/k_e) between 0.993 and 1.004 at every window end, for every seed and every m — the real loop at k_e = 25
follows its averaged slow dynamics to within 0.7% (across different clouds the last-window maxima are 1.57, 0.78 and 2.58 times the
averaged runs' 0.377, 0.117, 0.0114, within the seed-to-seed spread). Two of the three m = 16 clouds have inadmissible or marginal
weights (splitting margins −0.0021, −0.00001; the third +0.0026), so their eventual configuration is a saddle's; at m = 64 and 256 all
margins are positive (+0.0029 to +0.0049). The saturation at the orbit distance of Corollary D‴.7 lies beyond the horizon run.

*Reading.* In the window κ_R[μ*] < σ₀ < σ_crit the network learns (the collapsed configuration is unstable, the three kernels fit the
plant to 0.244) and propagation of chaos is uniform in time on a warm start with admissible weights: the coupled distance saturates with
the 1/m law at the stationary distance of the replicated weighted configuration, which the exact form of Proposition D‴.5(ii) gives
within 2.3% — the level is the weight-sensitivity law, not the theorem's constant. At m = 16 the multinomial weights leave the admissible
set with appreciable probability; the one such seed sits on a saddle and later splits a cluster. The cold start of a connected cloud
shows what the hypothesis excludes: growth over the horizon and a floor of pairs in different basins whose fraction decays like
m^{−0.64}, while the learned function still converges at a rate between m^{−1/2} and m^{−2/3} in both cases. The real loop at k_e = 25
follows its averaged slow dynamics to within 0.7% over the horizon run.

*Checked* (`tests/test_chaos_learning.py`): (𝔄_m)'s drift equals −G∇J_m of `periodic.ClosedLoop` to machine precision
(`test_averaged_slow_loop_is_the_G_gradient_flow_of_J_m`); the three-term identity of Step 1 holds to rounding on random
configurations with a mean-field stand-in, the interaction term equals −‖f̂_Θ − f̂_Θ̄‖² up to a relative O(D^{1/2}) correction (ratio
0.98 → 0.9999 as δ is scaled by 1 → 0.01) and obeys (3.1) at every scale (`test_coupling_identity_three_terms_and_the_sign_of_the_interaction`);
κ_R at the atoms of a collapsed configuration with the collapsed residual equals Proposition K's λ_max(G^{1/2}R₂G^{1/2}), which
equals σ₀ at σ₀ = σ_crit, the collapsed critical point has zero drift, and the deterministic neighbourhood search over the sphere of
radius ϱ reproduces the first-order expansion λ_max + ϱ‖G^{1/2}∇λ_max‖ with a relative error O(ϱ) (`test_residual_curvature_…`);
the searched whitened constants reproduce the document's values and lie below their closed forms, `theorem_constants` satisfies the
four constraints the bootstrap needs, and Λ₁ integrates exactly to s₁ (`test_whitened_constants_…`); at σ₀ = 0.15σ_crit the
three-atom critical point of J₃ has zero drift and (H) holds at it with ϱ′ = 0.1 inside the learning regime, the replicated weighted
stationary point is a critical point of J₁₆ to 10⁻¹⁰, the linearisation in the weights has a second-order remainder, and on a coupled
run from the jittered mixture (m = 16, 64 against M = 512, 60 slow units) κ_R[μ*] over the neighbourhood of the reference's support
stays below σ₀, the two curvature statistics differ by at most the bound (0.2), the per-cluster displacements point along
S_k(w^m − w^M) (cosine ≥ 0.8) with D within a factor 2 of Σ_kw^m_k‖S_k(w^m − w^M)‖², ‖f̂_Θ − f̂_Θ̄‖ ≤ L_ϑ^G·mean‖δ_i‖, D₆₄ < D₁₆
and the (vacuous) bound C₁/m holds (`test_hypothesis_H_at_the_three_atom_stationary_measure_and_the_weight_sensitivity_law`).

---

## 6. Theorem E — label-free local convergence under tangential excitation

**Assumptions.**

- (S0) holds, in particular k_e > 0.
- **Uniform tube.** There is ρ > 0 such that ε(x,t) = 0 whenever ‖x − x_m(t)‖ < ρ, for all t ≥ 0.
- **Interior ideal.** μ* is constant, with ϑ*_j ∈ int Θ_K^{ε}. Put δ := min_j dist(ϑ*_j, ∂Θ_K^{ε}) > 0.
- x_m ∈ C¹, and x_m and ẋ_m are bounded.
- Define Ψ(t) := [∇_ϑΦ(x_m(t);ϑ*_1), …, ∇_ϑΦ(x_m(t);ϑ*_m)] ∈ ℝ^{n×mp}.
- **(PE)** There are T, α > 0 with ∫_t^{t+T}ΨᵀΨ ⪰ αI_{mp} for all t ≥ 0.

(PE) forces a*_j ≠ 0 and pairwise distinct pairs (c*_j, s*_j).

- If a*_j = 0, the (c,s)-columns of block j vanish.
- If (c*_j, s*_j) = (c*_k, s*_k), the kernel of ∫ΨᵀΨ contains the vector with v in the a-slot of block j, −v in the
  a-slot of block k, and zeros elsewhere.

In particular the ϑ*_j are pairwise distinct, and

    d² := min_{π≠π'}(1/m)Σ‖ϑ*_{π(i)} − ϑ*_{π'(i)}‖²_{G⁻¹} = (2/m)·min_{i≠j}‖ϑ*_i − ϑ*_j‖²_{G⁻¹} > 0,

with d² := +∞ for m = 1. The closed form holds because a non-identity permutation moves at least two indices, and the
transposition of the G⁻¹-closest pair attains the bound.

**Lemma E.0 (explicit strict Lyapunov function for the skew-symmetric mixed class).** Let Ψ: [0,∞) → ℝ^{n×q} be C¹
with ‖Ψ(t)‖ ≤ Ψ_M and ‖Ψ̇(t)‖ ≤ Ψ_D for all t, and persistently exciting: there are T, α > 0 with
∫_t^{t+T}Ψ(r)ᵀΨ(r)dr ⪰ αI_q for all t ≥ 0. Let k > 0 and consider

    ė = −k e − Ψ(t)ζ,    ζ̇ = Ψ(t)ᵀe,    w := (e, ζ) ∈ ℝ^{n+q}.                                       (LTV)

Put M := ΨᵀΨ and

    Q(t) := ∫_t^{t+T}(t + T − r)M(r)dr,    Q_M := T²Ψ_M²/2,    c₁ := T²Ψ_M³ + TkΨ_M + TΨ_D,
    V(t, e, ζ) := ½‖e‖² + ½‖ζ‖² − εζᵀQ(t)ζ + εT eᵀΨ(t)ζ.

(i) *Every admissible ε.* For every ε > 0 with ε(Q_M + TΨ_M/2) < ½, V is C¹ and, along (LTV),

    p_lo‖w‖² ≤ V ≤ p_hi‖w‖²,    p_lo/hi := ½ ∓ ε(Q_M + TΨ_M/2),
    V̇ ≤ −(k − εTΨ_M²)‖e‖² − εα‖ζ‖² + εc₁‖e‖‖ζ‖.

Hence V̇ ≤ −γ(ε)‖w‖² with γ(ε) := λ_min[[a, −c/2], [−c/2, b]], a := k − εTΨ_M², b := εα, c := εc₁, whenever this is
positive, and then ‖w(t)‖ ≤ (p_hi/p_lo)^{1/2}e^{−(γ/(2p_hi))(t−t₀)}‖w(t₀)‖ for all t ≥ t₀ ≥ 0.

(ii) *Closed form.* With ε := min{1/(4(Q_M + TΨ_M/2)), k/(4TΨ_M²), kα/c₁²}, one has p_lo ≥ ¼, p_hi ≤ ¾ and
V̇ ≤ −γ‖w‖² with γ := min{k/4, εα/2}. Hence

    ‖w(t)‖ ≤ √3 · e^{−(2γ/3)(t − t₀)} · ‖w(t₀)‖    for all t ≥ t₀ ≥ 0.

(iii) *Dependence on Ψ_D and on time rescaling.* The rate tends to 0 as Ψ_D → ∞, and it must: Ψ = sin t² is bounded
and PE with Ψ_D = ∞, and (LTV) is then not uniformly exponentially stable (proof below). Under the time rescaling
Ψ_ω(t) := Ψ(ωt), with Ψ_D > 0, the data transform as (α, T, Ψ_M, Ψ_D) ↦ (α/ω, T/ω, Ψ_M, ωΨ_D), so
c₁ = TΨ_D·(1 + O(1/ω)) and the closed-form rate is kα²/(3ω²T²Ψ_D²)·(1 + O(1/ω)) as ω → ∞. This is the ω⁻² law of
Remark E.1.

*Proof.* (a) *The Gramian weight.* 0 ⪯ Q(t) ⪯ ∫_t^{t+T}(t + T − r)Ψ_M²dr·I = Q_M I. Since M is continuous, Q is C¹
with Q̇(t) = ∫_t^{t+T}M(r)dr − TM(t): the upper limit contributes (t + T − (t + T))M(t + T) = 0, the lower limit
contributes −TM(t), and differentiating the weight under the integral gives ∫_t^{t+T}M.

(b) *Bounds on V.* |ζᵀQζ| ≤ Q_M‖ζ‖² and |eᵀΨζ| ≤ Ψ_M‖e‖‖ζ‖ ≤ (Ψ_M/2)(‖e‖² + ‖ζ‖²).

(c) *Three derivatives along (LTV).*

- d/dt ½(‖e‖² + ‖ζ‖²) = −k‖e‖² − eᵀΨζ + ζᵀΨᵀe = −k‖e‖².
- d/dt(−ζᵀQζ) = −ζᵀQ̇ζ − 2ζᵀQζ̇ = −ζᵀ(∫_t^{t+T}M)ζ + T‖Ψζ‖² − 2ζᵀQΨᵀe ≤ −α‖ζ‖² + T‖Ψζ‖² − 2ζᵀQΨᵀe.
- d/dt(eᵀΨζ) = ėᵀΨζ + eᵀΨ̇ζ + eᵀΨζ̇ = −k eᵀΨζ − ‖Ψζ‖² + eᵀΨ̇ζ + ‖Ψᵀe‖².

Adding them with weights 1, ε and εT, the ‖Ψζ‖² terms cancel:

    V̇ ≤ −k‖e‖² − εα‖ζ‖² − 2εζᵀQΨᵀe − εTk eᵀΨζ + εT eᵀΨ̇ζ + εT‖Ψᵀe‖².

The last four terms are at most 2εQ_MΨ_M‖e‖‖ζ‖, εTkΨ_M‖e‖‖ζ‖, εTΨ_D‖e‖‖ζ‖ and εTΨ_M²‖e‖². This is (i); the
quadratic form −a‖e‖² − b‖ζ‖² + c‖e‖‖ζ‖ is at most −λ_min[[a, −c/2], [−c/2, b]]·‖w‖², and V ≤ p_hi‖w‖², ‖w‖² ≤ V/p_lo
give the estimate.

(d) *The closed form.* Young's inequality gives εc₁‖e‖‖ζ‖ ≤ (εα/2)‖ζ‖² + (εc₁²/(2α))‖e‖². With ε ≤ kα/c₁² the second
coefficient is at most k/2, and with ε ≤ k/(4TΨ_M²) the coefficient of ‖e‖² in (i) is at least k − k/4. So
V̇ ≤ −(k/4)‖e‖² − (εα/2)‖ζ‖² ≤ −γ‖w‖². Since V ≤ (3/4)‖w‖², V̇ ≤ −(4γ/3)V, and ‖w‖² ≤ 4V gives
‖w(t)‖² ≤ 4V(t) ≤ 4e^{−(4γ/3)(t−t₀)}V(t₀) ≤ 3e^{−(4γ/3)(t−t₀)}‖w(t₀)‖².

(e) *Rescaling.* Ψ_ω has sup norm Ψ_M, derivative bound ωΨ_D, and ∫_t^{t+T/ω}Ψ_ωᵀΨ_ω dr = (1/ω)∫_{ωt}^{ωt+T}ΨᵀΨ ⪰ (α/ω)I.
Substitute into c₁ and the closed form.

(f) *Ψ = sin t² is not uniformly exponentially stabilising* (n = q = 1, k > 0). Let F(s) := ∫_s^∞ sin r² dr, so that
sin s² = −F'(s) and |F(s)| ≤ C/s (Fresnel). From w(t₀) = (0, ζ₀), e(t) = −∫_{t₀}^t e^{−k(t−s)}sin(s²)ζ(s)ds; integrating
by parts and using |ζ| ≤ ‖w(t₀)‖ (the energy ½‖w‖² is non-increasing) and |ζ̇| = |Ψe| ≤ |e|,

    |e(t)| ≤ (C/t₀)(2 + 1)‖w(t₀)‖ + (C/t₀)∫_{t₀}^t e^{−k(t−s)}|e(s)|ds,

so sup_{t≥t₀}|e| ≤ 3C‖w(t₀)‖/(t₀(1 − C/(kt₀))) for t₀ > C/k. Hence ½‖w(t₀ + Δ)‖² ≥ ½‖w(t₀)‖² − kΔ(3C/t₀)²‖w(t₀)‖²(1 + o(1)):
the decay over any fixed Δ tends to 0 as t₀ → ∞, which excludes a uniform exponential rate. The three ratios quoted in
step 3 of the proof of Theorem E (0.9715, 0.9942, 0.9980 for t₀ = 5, 20, 50 and Δ = 200) match 1 − (1/8)(1/t₀ − 1/(t₀ + 200))
to three digits. ∎

*Checked* (`tests/test_excitation.py`): the two bounds on V and V̇ ≤ −γ‖w‖² at random (t, w) for four regressors of
different shapes and gains, with the closed-form ε and with the ε maximising the rate (`excitation.optimise_eps`; every
admissible ε gives a valid Lyapunov function, so the maximum over the grid is a valid bound for the given PE data; the
data Ψ_M, Ψ_D, α are themselves grid maxima and a refinement-checked quadrature, so the ordering against the exact rate
is numerical); the envelope along trajectories; the ordering explicit rate ≤ exact Floquet rate; the ω⁻² law; the
slow-adaptation limit; and the explicit nonlinear envelope of Theorem E on the implemented closed loop. The round-5
referee verified the pointwise form λ_max(N(t)) ≤ −γ of V̇ ≤ −γ‖w‖² over all w, with an independent quadrature, for
ten regressor families (q > n, q < n, k from 10⁻³ to 10⁴, Ψ_D/Ψ_M up to 100, non-periodic PE, near-degenerate PE) and
eight choices of ε, without finding a violation.

*Remark E.0a (how conservative the explicit rate is; `exp_explicit_rate.py` → `results/explicit_rate/table.md`).* The
construction is of the Mazenc–Malisoff type, a backward–forward integral of the PE Gramian as a cross term; its content
here is the explicit, self-contained statement for the mixed class (LTV), which removes the dependence of Theorem E on
[S8]. The price is conservativeness, which the table quantifies against the exact Floquet rate of the linearised loop.

- *Slow adaptation* (TΨ_M² ≪ k). To leading order in the gains, c₁ ≈ T(kΨ_M + Ψ_D) and α = Tλ_min⟨ΨᵀΨ⟩ for a window
  equal to the period, so the closed-form bound (2/3)γ_cf of (ii) is kα²/(3c₁²), and the rate with the optimal ε of
  the same quadratic form is kα²/c₁² (maximise εα − ε²c₁²/(4k) over ε). The exact rate is λ_min⟨ΨᵀΨ⟩/k when k is also
  large compared with the bandwidth of the reference (Remark E.1). Hence, in that doubly asymptotic regime,

      exact / explicit ≈ 3κ_Ψ(1 + Ψ_D/(kΨ_M))²  for (ii),    ≈ κ_Ψ(1 + Ψ_D/(kΨ_M))²  for the optimised ε,

  with κ_Ψ := Ψ_M²/λ_min⟨ΨᵀΨ⟩ the condition number of the averaged Gramian. The exact rate pays for ill-conditioning
  once, through λ_min; the bound pays twice, because the Young step that absorbs the cross term εc₁‖e‖‖ζ‖ compares
  c₁ ∝ Ψ_M with α ∝ λ_min (Lemma E.0′ below removes this by slaving the error to the slow manifold before the Young step). On the ideal of the Theorem E tests (m = 3, n = 2, k_e = 6), driven by a two-tone reference
  of period 20π, with G = g·diag(2,2,1): κ_Ψ = 1.27·10⁵, (1 + Ψ_D/(k_eΨ_M))² = 3.84, so the two formula values are
  1.47·10⁶ and 4.9·10⁵. The measured ratio with the optimised ε falls from 6.56·10¹⁰ at g = 30 to 4.22·10⁵ at g = 0.01
  (5.63·10⁵ with the closed-form ε, which the table evaluates with the exact 2 × 2 form and the exact p_hi rather than
  the bound (ii)). Two corrections account for the difference from 4.9·10⁵: on this instance the exact rate is 0.66
  of λ_min⟨ΨᵀΨ⟩/k at every g (k_e = 6 is not large against the reference's bandwidth; the two last columns of the table),
  and the T²Ψ_M³ term is still 13% of c₁ at g = 0.01 ((c₁/c₁^∞)² = 1.31); their product, 0.66 × 1.31 = 0.86, brings
  4.9·10⁵ to 4.2·10⁵. The ratio is still decreasing at g = 0.01 and tends to about 0.66·κ_Ψ(1 + Ψ_D/(k_eΨ_M))² ≈ 3.2·10⁵
  as g → 0. In the mean-field normalisation Ψ_G ∝ m^{−1/2}, so large networks sit in this regime at fixed geometry.
- *Fast adaptation* (TΨ_M² ≳ k). The ratio grows further with the excitation strength, through the T²Ψ_M³ term of c₁
  and the Q_M term of ε.
- *Remark E.1's instance* (m = n = 1, G = I, k_e = 3, x_m = 1.5 sin ωt). The ratio exact/explicit decreases from 954 at
  ω = 1 to 17.6 at ω = 64. Over ω ≥ 4 the log-log slopes are −1.97 (exact) and −1.56 (explicit); over the last octave,
  ω = 32 → 64, the explicit slope is −1.92, approaching the asymptotic ω⁻² of (iii).
- *The window.* Any window length T is admissible with its own α(T). On the test ideal the best window is half a period.
- *Answered below.* Lemma E.0′ is a strict Lyapunov function whose rate is λ_min⟨ΨᵀΨ⟩/k up to a factor 1 − δ, above a
  threshold in k set by the peak factors of the regressor: in the slow regime it pays for the ill-conditioning once.

**Lemma E.0′ (the slow-adaptation regime: a strict Lyapunov function whose rate is the averaged rate).** Let Ψ be
T-periodic and C², M := ΨᵀΨ, M̄ := (1/T)∫₀^T M dt ≻ 0, ‖v‖_{M̄} := (vᵀM̄v)^{1/2}, λ_min := λ_min(M̄), and

    Ψ_M := sup‖Ψ‖,  Ψ_D := sup‖Ψ̇‖,  Ψ_DD := sup‖Ψ̈‖,
    ρ := sup_{t, v≠0} ‖Ψ(t)v‖²/‖v‖²_{M̄},   ρ_D := sup ‖Ψ̇(t)v‖²/‖v‖²_{M̄},   ρ_DD := sup ‖Ψ̈(t)v‖²/‖v‖²_{M̄},

the *peak factors*: the largest ratio, over directions, between the instantaneous excitation of a direction (by Ψ, Ψ̇, Ψ̈)
and its average excitation; ρ ≥ 1, and all three are finite because M̄ ≻ 0. Write ε := 1/k and let

    P₁(t) := −(2/T)∫_t^{t+T}(t+T−r)(M(r) − M̄)dr          (the zero-mean antiderivative of 2(M − M̄)),
    F := P₁M + MP₁,   S := the zero-mean antiderivative of F,   P₂ := −(M − M̄) + S,
    e₂(ζ,t) := −εΨζ + ε²Ψ̇ζ − ε³ΨMζ   (the slow manifold of the error, expanded through order ε³),   w := e − e₂(ζ,t),
    V(e,ζ,t) := ½ζᵀ(I + εP₁(t) + ε²P₂(t))ζ + (β/2)‖w‖².

F has zero mean: ⟨P₁M + MP₁⟩ = ⟨P₁⟩M̄ + M̄⟨P₁⟩ + ½⟨d/dt(P₁²)⟩ = 0, so S exists and the second-order averaged
correction N₂ := ½⟨F⟩ vanishes identically. Put

    π₁ := (T/2)Ψ_M(√ρ + 1),   π₂ := Ψ_M(√ρ + 1) + (T²/8)Ψ_M³(2√ρ + 1),   p₁ := TΨ_M²/2,   p₂ := Ψ_M² + (T²/4)Ψ_M⁴,
    c₁ := √ρ + εΨ_Mπ₁ + ε²Ψ_Mπ₂,
    c₂ := ε²√ρ_DD + ε³(3Ψ_DΨ_M√ρ + 2Ψ_M²√ρ_D) + ε⁴(2Ψ_M⁴√ρ + Ψ_DΨ_M√ρ_D) + ε⁵(Ψ_M⁴√ρ_D + Ψ_DΨ_M³√ρ) + ε⁶Ψ_M⁶√ρ,
    d_w := εΨ_M² + ε²Ψ_DΨ_M + ε³Ψ_M⁴,
    r_ζ := ε³(π₂Ψ_M√ρ + π₁Ψ_M√ρ_D + Ψ_M²ρ) + ε⁴(π₂Ψ_M√ρ_D + π₁Ψ_M³√ρ) + ε⁵π₂Ψ_M³√ρ,
    β := 2c₁²/δ,   c_ζ := 1 − δ − r_ζ/ε.

Fix δ ∈ (0, 1). If

    (A) c₁c₂ ≤ δ/2,   (B) d_w ≤ k/4,   (C) c_ζ > 0,   (D) εp₁ + ε²p₂ < 1,

then along (LTV), for all (e, ζ, t),

    V̇ ≤ −ε c_ζ ζᵀM̄ζ − (βk/4)‖w‖² ≤ −γV,    γ := min{ 2c_ζλ_min/(k(1 + εp₁ + ε²p₂)), k/2 },

and ‖(e,ζ)(t)‖ ≤ C e^{−(γ/2)(t−t₀)}‖(e,ζ)(t₀)‖ with C := (1 + εΨ_M + ε²Ψ_D + ε³Ψ_M³)²·(max{1 + εp₁ + ε²p₂, β}/min{1 − εp₁ − ε²p₂, β})^{1/2}.

In the slow-adaptation limit (k → ∞ at fixed Ψ) the conditions hold for every δ ∈ (0, 1) once k ≥ k₀(δ), with
k₀ ≈ (4ρρ_DD/δ²)^{1/4} when (A) binds (small δ, or TΨ_M² ≲ 1; otherwise (C) sets k₀ through r_ζ, as on the sinusoidal
regressor below), and the rate is γ/2 = (1 − δ)(1 − εp₁ + O(k⁻²))λ_min/k: the averaged rate λ_min⟨ΨᵀΨ⟩/k of Remark E.1 up
to the factor 1 − δ and the first-order deficit εp₁ = TΨ_M²/(2k) of the quadratic form. The bound pays for the
ill-conditioning of M̄ once, through λ_min, as the exact rate does; the peak factors enter only the threshold k₀ and the
overshoot C = O((ρ/δ)^{1/2}), not the rate.

*Status.* Proved (below); checked pointwise by a certificate — the supremum of V̇/V over all states and over a time
grid is a generalised eigenvalue problem — on the sinusoidal regressor and on the slow-adaptation ideal of the Theorem E
tests, where the guaranteed rate is within 12% (δ = 0.1) and 26% (δ = 0.2) of the exact rate at k = 400 while Lemma E.0's
is five orders of magnitude below it (`tests/test_excitation.py`, `exp_explicit_rate.py`, table 4). The round-6 referee
re-derived every step and certified the inequality independently on five regressor families without a failure
(`review/REVIEW_round6_E0prime.md`).

*Proof.* Write M, Ψ, … for their values at t. Two facts about zero-mean antiderivatives are used throughout: if f is
T-periodic with zero mean and g is its zero-mean antiderivative, then ‖g‖_∞ ≤ (T/4)‖f‖_∞ — for the maximiser t* of g,
g(t*) = (1/T)∫₀^T(g(t*) − g(s))ds ≤ (1/T)∫₀^T dist(t*, s)‖f‖_∞ ds = (T/4)‖f‖_∞ with dist the distance on the circle, and
the same for −g — and the same bound holds for every bilinear coordinate uᵀg v, hence ‖g(t)v‖ ≤ (T/4)sup_r‖f(r)v‖. As
‖M − M̄‖ ≤ Ψ_M² (two positive semidefinite matrices of norm ≤ Ψ_M²), this gives ‖P₁‖ ≤ TΨ_M²/2 = p₁ and
‖P₁v‖ ≤ (T/2)sup‖(M − M̄)v‖, and ‖S‖ ≤ (T/4)‖F‖_∞, ‖Sv‖ ≤ (T/4)sup‖Fv‖.

(a) *The fast variable.* From ė = −ke − Ψζ and ke₂ = −Ψζ + εΨ̇ζ − ε²ΨMζ, ė = −kw − εΨ̇ζ + ε²ΨMζ. Differentiating e₂ along
the flow, with ζ̇ = Ψᵀe = Ψᵀw − εMζ + ε²ΨᵀΨ̇ζ − ε³M²ζ and Ṁ = Ψ̇ᵀΨ + ΨᵀΨ̇, the terms of order ε and ε² in ζ cancel and

    ẇ = −kw + εΨΨᵀw − ε²Ψ̇Ψᵀw + ε³ΨMΨᵀw − ε²Ψ̈ζ + ε³(2Ψ̇M + ΨΨ̇ᵀΨ + 2ΨΨᵀΨ̇)ζ − ε⁴(2ΨM² + Ψ̇ΨᵀΨ̇)ζ + ε⁵(ΨMΨᵀΨ̇ + Ψ̇M²)ζ − ε⁶ΨM³ζ.

(b) *The slow part of V̇.* With H := I + εP₁ + ε²P₂ (symmetric), d/dt[½ζᵀHζ] = ζᵀHζ̇ + ½ζᵀḢζ and Ḣ = 2ε(M − M̄) +
ε²(F − Ṁ). Insert ζ̇ and collect by powers of ε. At order ε: −ζᵀMζ + ζᵀ(M − M̄)ζ = −ζᵀM̄ζ. At order ε²:
−ζᵀP₁Mζ + ζᵀΨᵀΨ̇ζ + ½ζᵀ(F − Ṁ)ζ = −½ζᵀ(P₁M + MP₁)ζ + ½ζᵀṀζ + ½ζᵀ(P₁M + MP₁ − Ṁ)ζ = 0, using ζᵀΨᵀΨ̇ζ = ½ζᵀṀζ: the
choice of P₂ removes the whole second-order term, and this is where N₂ = 0 enters (otherwise −ζᵀN₂ζ would remain).
The remaining terms are R_ζ := ε³(−ζᵀP₂Mζ + ζᵀP₁ΨᵀΨ̇ζ − ζᵀM²ζ) + ε⁴(ζᵀP₂ΨᵀΨ̇ζ − ζᵀP₁M²ζ) − ε⁵ζᵀP₂M²ζ, plus the coupling
ζᵀHΨᵀw, so d/dt[½ζᵀHζ] = −εζᵀM̄ζ + R_ζ + ζᵀHΨᵀw.

(c) *Per-direction bounds.* For every v and t: ‖Ψv‖ ≤ √ρ‖v‖_{M̄}, ‖Ψ̇v‖ ≤ √ρ_D‖v‖_{M̄}, ‖Ψ̈v‖ ≤ √ρ_DD‖v‖_{M̄};
‖Mv‖ ≤ Ψ_M‖Ψv‖ ≤ Ψ_M√ρ‖v‖_{M̄}; ‖M̄v‖ ≤ λ_max(M̄)^{1/2}‖v‖_{M̄} ≤ Ψ_M‖v‖_{M̄}; ‖ΨᵀΨ̇v‖ ≤ Ψ_M√ρ_D‖v‖_{M̄};
‖P₁v‖ ≤ (T/2)sup_r‖(M(r) − M̄)v‖ ≤ π₁‖v‖_{M̄}; ‖Fv‖ ≤ ‖P₁Mv‖ + ‖MP₁v‖ ≤ (T/2)Ψ_M³√ρ‖v‖_{M̄} + Ψ_M²π₁‖v‖_{M̄} =
TΨ_M³(√ρ + ½)‖v‖_{M̄}, so ‖Sv‖ ≤ (T/4)sup‖Fv‖ = (T²/8)Ψ_M³(2√ρ + 1)‖v‖_{M̄} and ‖P₂v‖ ≤ ‖(M − M̄)v‖ + ‖Sv‖ ≤ π₂‖v‖_{M̄}; also
‖P₂‖ ≤ Ψ_M² + (T/4)‖F‖_∞ ≤ Ψ_M² + (T/4)·2·(TΨ_M²/2)Ψ_M² = p₂. Hence |R_ζ| ≤ r_ζ‖ζ‖²_{M̄}, |ζᵀHΨᵀw| = |⟨ΨHζ, w⟩| ≤ c₁‖ζ‖_{M̄}‖w‖, and from (a),
βwᵀẇ ≤ β[−(k − d_w)‖w‖² + c₂‖ζ‖_{M̄}‖w‖]: each ζ-term of ẇ is bounded through ‖Ψζ‖, ‖Ψ̇ζ‖, ‖Ψ̈ζ‖, ‖Mζ‖ or ‖M²ζ‖ (the terms
of c₂, in the order of (a); the last one is ‖ΨM³ζ‖ ≤ Ψ_M·Ψ_M⁴‖Mζ‖ ≤ Ψ_M⁶√ρ‖ζ‖_{M̄}), the w-terms through operator norms (d_w).

(d) *Assembly.* V̇ ≤ −ε‖ζ‖²_{M̄} + r_ζ‖ζ‖²_{M̄} + (c₁ + βc₂)‖ζ‖_{M̄}‖w‖ − β(k − d_w)‖w‖². By Young,
(c₁ + βc₂)‖ζ‖_{M̄}‖w‖ ≤ δε‖ζ‖²_{M̄} + (c₁ + βc₂)²‖w‖²/(4δε), and (c₁ + βc₂)² ≤ 2c₁² + 2β²c₂². With β = 2c₁²/δ the term
2c₁²/(4δε) equals βk/4, and (A) gives 2β²c₂²/(4δε) ≤ βk/4, so with (B) the coefficient of ‖w‖² is at most
−β(k − d_w − k/2) ≤ −βk/4; the coefficient of ‖ζ‖²_{M̄} is −ε(1 − δ − r_ζ/ε) = −εc_ζ. This is the first inequality. Since
½(1 − εp₁ − ε²p₂)‖ζ‖² + (β/2)‖w‖² ≤ V ≤ ½(1 + εp₁ + ε²p₂)‖ζ‖² + (β/2)‖w‖² by (D), and ζᵀM̄ζ ≥ λ_min‖ζ‖², V̇ ≤ −γV follows,
and the envelope with C from ‖e‖ ≤ ‖w‖ + (εΨ_M + ε²Ψ_D + ε³Ψ_M³)‖ζ‖ and its converse. ∎

*Remark E.0′.1 (reading).*

- *What is sharp and what is not.* The exact rate of (LTV) is λ_min/k·(1 + O(k⁻²)): the first-order averaging of
  ζ̇ = −(1/k)Mζ after the slaving gives λ_min/k, and the second-order correction to the symmetric part of the averaged
  generator is N₂ = 0 (table 4 measures 1 − exact/(λ_min/k) ≈ 29/k² on the ideal; on Ψ = (sin t, cos t) the deviation is
  of order k⁻⁴). Lemma E.0′ recovers this rate up to the factor c_ζ/(1 + εp₁ + ε²p₂) → 1 − δ, and δ can be taken → 0 as
  k → ∞ (condition (A) allows δ ≳ 2√(ρρ_DD)/k²). The overshoot constant C = O((ρ/δ)^{1/2}) is an artefact of the
  quadratic form, not of the dynamics.
- *Two regimes, two lemmas.* Lemma E.0 needs only (α, T, Ψ_M, Ψ_D) and gives the ω⁻² law of the fast regime; its rate
  is quadratic in λ_min. Lemma E.0′ needs k above the threshold k₀ ≈ (4ρρ_DD/δ²)^{1/4}, set by the peak factors of Ψ
  and its two derivatives, and gives the sharp rate. On the regressors of Theorem E the derivative peak factors are
  large — the directions that are poorly excited on average are excited in short bursts, so ρ_DD/ρ is of the order of
  the fourth power of the effective bandwidth — and k₀ is in the hundreds there, while on Ψ = (sin t, cos t)
  (ρ = ρ_D = ρ_DD = 2) it is about 10, set by (C) rather than (A). A third-order slaving would replace ρ_DD by the peak factor of Ψ⃛ at one order
  higher in ε and lower k₀ further; this is not pursued.
- *Hypotheses.* Only the periodicity of Ψ is used (through the zero-mean antiderivatives; a non-periodic PE regressor
  would need a windowed M̄(t) and its derivative). Theorem E's x_m is periodic, so Ψ_G is. Theorem E's steps 4–5 go
  through with this V in place of Lemma E.0's (Corollary E.0′.2 below): the nonlinear rate becomes γ/4 ≈ (1 − δ)λ_min/(2k),
  while the explicit ball radius carries the anisotropy p′_lo/p′_hi ≈ 1/β of the quadratic form and can be smaller than
  Lemma E.0's.

*Numerics* (`exp_explicit_rate.py` → `results/explicit_rate/table.md`, table 4). The ideal of table 2 at g = 0.01
(T = 20π): λ_min(M̄) = 2.16·10⁻⁷, κ(M̄) = 7.8·10³, Ψ_M = 0.168, Ψ_D = 0.957, Ψ_DD = 8.39, TΨ_M² = 1.78; peak factors
ρ = 69.9, ρ_D = 3.7·10³, ρ_DD = 3.7·10⁵ (the poorly excited directions are excited in bursts: ρ_D/ρ ≈ 53 and ρ_DD/ρ_D ≈ 99
are effective squared bandwidths); ‖N₂‖ = 3·10⁻²¹. The threshold of (A)–(D) is k₀ = 226 for δ = 0.2 and 321 for δ = 0.1,
set by (A): c₁c₂ = 2.07, 0.51, 0.13, 0.056, 0.032, 0.014 at k = 50, 100, 200, 300, 400, 600, the k⁻² law of ε²√(ρρ_DD).

| k | exact Floquet rate | exact / (λ_min/k) | guaranteed rate γ/2, δ = 0.2 (δ = 0.1) | exact / guaranteed | exact / Lemma E.0 (optimised ε) | sup V̇/V on the grid, −γ (δ = 0.1) |
|---|---|---|---|---|---|---|
| 300 | 7.21·10⁻¹⁰ | 0.9997 | 5.70·10⁻¹⁰ (below threshold) | 1.26 (—) | 1.33·10⁵ | — |
| 400 | 5.41·10⁻¹⁰ | 0.9998 | 4.30·10⁻¹⁰ (4.84·10⁻¹⁰) | 1.26 (1.12) | 1.32·10⁵ | −1.079·10⁻⁹, −9.68·10⁻¹⁰ |
| 600 | 3.61·10⁻¹⁰ | 0.9999 | 2.88·10⁻¹⁰ (3.24·10⁻¹⁰) | 1.25 (1.12) | 1.30·10⁵ | −7.20·10⁻¹⁰, −6.47·10⁻¹⁰ |

Above the threshold the guaranteed rate is within 26% (δ = 0.2) and 12% (δ = 0.1) of the exact rate — the factor
1/(1 − δ) of the lemma, with c_ζ = 0.80 and 0.90 — while Lemma E.0's optimised rate is 1.3·10⁵ times smaller; the
certificate sup V̇/V ≤ −γ holds on the 128-point grid in every row where the conditions hold, with sup V̇/V = −1.078·10⁻⁹
≈ −2λ_min/k at k = 400 (attained in the slow direction at w = 0). The exact rate deviates from λ_min/k by 1.15%, 0.29%,
0.07% at k = 50, 100, 200: the k⁻² law that N₂ = 0 predicts. The overshoot constant is C = √β·(1 + O(1/k)) = 27 (δ = 0.2)
and 38 (δ = 0.1). Below k₀ the lemma makes no claim; Lemma E.0 does (table 2).

*Checked* (`tests/test_excitation.py`): on Ψ = (sin t, cos t) the computed λ_min, ρ, ρ_D, ρ_DD equal 1/2, 2, 2, 2, N₂ = 0 to
10⁻¹², the threshold of (A)–(D) is below 12 for δ = 0.2, the certificate sup V̇/V ≤ −γ holds at k = 20, 50, 200, and the
guaranteed rate is below the exact Floquet rate and within a factor 2 of it, the ratio decreasing towards 1/(1 − δ)
(`test_lemma_E0prime_sinusoidal_regressor_certificate_and_sharpness`); on the slow-adaptation ideal of the Theorem E tests
at g = 0.01 and k = 400, δ = 0.2: κ(M̄) > 10³, the conditions hold, the certificate holds on a 128-point grid, the
guaranteed rate is within a factor 1.5 of the exact rate and exceeds Lemma E.0's optimised rate by more than 10⁴
(`test_lemma_E0prime_pays_for_ill_conditioning_once_on_the_theorem_E_ideal`).

**Corollary E.0′.2 (Theorem E in the slow regime).** Under the hypotheses of Theorem E, suppose the regressor Ψ_G of its
step 3 satisfies the conditions (A)–(D) of Lemma E.0′ at k = k_e for some δ, with constants (β, γ, p₁, p₂) and
L_M := εΨ_{G,M} + ε²Ψ_{G,D} + ε³Ψ_{G,M}³ (ε = 1/k_e), and that βL_M² < 1 − εp₁ − ε²p₂. Put

    p′_lo := ½min{β/2, 1 − εp₁ − ε²p₂ − βL_M²},    p′_hi := ½max{β(1 + L_M), 1 + εp₁ + ε²p₂ + βL_M(1 + L_M)},
    c′₀ := min{ [γp′_lo/(4(p′_hi)^{1/2}C̃_R)]², p′_lo r²/‖S⁻¹‖² },    r′₀ := (c′₀/p′_hi)^{1/2},    λ′ := γ/4,    C′ := (2p′_hi/p′_lo)^{1/2},

with r, C̃_R, S as in Theorem E. Then the conclusions of Theorem E hold with (c₀, λ, C) replaced by (r′₀²/2, λ′, C′). In the
slow regime λ′ = (1 − δ)(1 − εp₁ + O(k_e⁻²))λ_min(M̄_G)/(2k_e), M̄_G := ⟨Ψ_GᵀΨ_G⟩: half the exact rate of the linearised
loop, up to the factor 1 − δ and the first-order deficit εp₁.

*Proof.* Write V := ½wᵀ𝐇(t)w for Lemma E.0′'s V in the coordinates w = (e, ζ) of step 3: 𝐇 = [[βI, −βL], [−βLᵀ, H_ζ + βLᵀL]]
with e₂ = L(t)ζ, ‖L‖ ≤ L_M and H_ζ := I + εP₁ + ε²P₂. Then wᵀ𝐇w = β‖e − Lζ‖² + ζᵀH_ζζ, and ‖e − Lζ‖² ≤ (1 + L_M)(‖e‖² +
L_M‖ζ‖²), ‖e − Lζ‖² ≥ ½‖e‖² − L_M²‖ζ‖² give 2p′_lo‖w‖² ≤ wᵀ𝐇w ≤ 2p′_hi‖w‖². Along the nonlinear loop ẇ = [A₀ + B(t)]w +
R̃(t, w) with ‖R̃‖ ≤ C̃_R‖w‖² on {‖w‖ < r/‖S⁻¹‖} (step 4), so with Lemma E.0′,

    V̇ ≤ −γV + ⟨𝐇w, R̃⟩ ≤ −γV + ‖𝐇^{1/2}w‖‖𝐇^{1/2}‖C̃_R‖w‖² ≤ −γV + 2(p′_hi)^{1/2}C̃_R V^{3/2}/p′_lo,

using ‖𝐇^{1/2}w‖ = (2V)^{1/2}, ‖𝐇^{1/2}‖ = (2p′_hi)^{1/2} and ‖w‖² ≤ V/p′_lo. Hence V̇ ≤ −(γ/2)V wherever
V^{1/2} ≤ γp′_lo/(4(p′_hi)^{1/2}C̃_R), and the sublevel set {V ≤ c′₀} — inside the ball where R̃ is bounded, by the second
term of c′₀ — is forward invariant with V(t) ≤ e^{−(γ/2)(t−t₀)}V(t₀) on it. It contains {‖w‖ ≤ r′₀} because V ≤ p′_hi‖w‖²,
and ‖w(t)‖ ≤ (p′_hi/p′_lo)^{1/2}e^{−(γ/4)(t−t₀)}‖w(t₀)‖ follows. Step 5 of Theorem E is unchanged: ‖w‖² = 2V_π, so
V_OT(t₀) ≤ r′₀²/2 puts the state in the ball, ‖e‖ + W ≤ √2‖w‖ gives C′, and uniqueness and the limits follow as there. The
asymptotics of λ′ are those of γ in Lemma E.0′. ∎

*Numerics* (`exp_explicit_rate.py` → `results/explicit_rate/table.md`, table 5). Remark E.1's single-kernel instance (m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = 0.3·I, x_m = 1.5 sin t) at
k_e = 40: λ_min(M̄_G) = 5.41·10⁻², κ(M̄_G) = 6.7, peak factors ρ = 4.8, ρ_D = 30.3, ρ_DD = 461, threshold k₀ = 24.4 (δ = 0.2)
and 33.5 (δ = 0.1); L_J = 4.62 estimated from samples; exact linear rate 1.35·10⁻³.

| Lyapunov function | explicit ball radius r₀ (w-coordinates) | nonlinear rate λ | overshoot C | exact linear rate / λ |
|---|---|---|---|---|
| Lemma E.0, optimised ε (Theorem E) | 6.15·10⁻⁶ | 2.63·10⁻⁵ | 1.43 | 51 |
| Lemma E.0′, δ = 0.2 (this corollary) | 1.53·10⁻⁶ | 4.53·10⁻⁴ | 11.8 | 3.0 |
| Lemma E.0′, δ = 0.1 (this corollary) | 8.32·10⁻⁷ | 5.16·10⁻⁴ | 17.1 | 2.6 |

On this well-conditioned regressor (κ ≈ 7) Lemma E.0 is only an order of magnitude conservative, and the corollary trades
the certified ball for the rate: λ′ is within a factor 3 of the exact linear rate (the factor 2 of the Lyapunov argument
times 1/(1 − δ)), while r′₀ is four times smaller than Theorem E's r₀, carrying the anisotropy p′_lo/p′_hi ≈ 1/β of the
quadratic form. On ill-conditioned regressors the rate improves by the factor of Remark E.0a (table 4);
the radius comparison was not computed there.

*Checked* (`tests/test_excitation.py`): on the single-kernel instance at k_e = 40, δ = 0.2, from a start with ‖w(0)‖ ≤ r′₀,
the implemented closed loop satisfies ‖w(t)‖ ≤ C′e^{−λ′t}‖w(0)‖ at every sample over 1.5 e-folds of the explicit envelope (dt = 0.005: the fixed-step loop has a discretisation floor ≈ 6·10⁻⁶ at dt = 0.02, above the explicit ball),
with λ′ more than ten times Theorem E's rate from Lemma E.0 and within a factor 4 of the exact linear rate
(`test_corollary_E0prime_2_envelope_on_the_implemented_loop`).

**Corollary E.0′.3 (the stationary parameter spread of the dithered linearised loop).** Under the hypotheses of Lemma E.0′
(conditions (A)–(D) at (k, δ)), consider (LTV) with parameter dither,

    ė = −k e − Ψ(t)ζ,    dζ = Ψ(t)ᵀe dt + Σ dB_t,                                                        (LTV_τ)

with Σ ∈ ℝ^{q×q′} constant and B a standard Brownian motion. Let V be the strict Lyapunov function of Lemma E.0′ and
𝐇(t) its quadratic form (V = ½xᵀ𝐇(t)x, x = (e, ζ)), L(t) := −εΨ + ε²Ψ̇ − ε³ΨM the slaving map (e₂ = L(t)ζ), L_M := εΨ_M +
ε²Ψ_D + ε³Ψ_M³ ≥ sup‖L‖, and

    ν(t) := ½ tr(Σᵀ𝐇_ζζ(t)Σ) = ½ tr(Σᵀ(I + εP₁ + ε²P₂ + βLᵀL)Σ),    ν̄ := sup_t ν(t) ≤ ½‖Σ‖_F²(1 + εp₁ + ε²p₂ + βL_M²).

Then, for every initial law with finite second moment,

    limsup_{t→∞} E V(t, x_t) ≤ ν̄/γ,    limsup_{t→∞} E‖ζ_t‖² ≤ 2ν̄/(γ(1 − εp₁ − ε²p₂)),    limsup_{t→∞} E‖w_t‖² ≤ 2ν̄/(γβ),

with γ the rate of Lemma E.0′; in the slow regime γ = 2c_ζλ_min(M̄)/(k(1 + εp₁ + ε²p₂)), so that

    limsup E‖ζ_t‖² ≤ k ‖Σ‖_F² (1 + εp₁ + ε²p₂)(1 + εp₁ + ε²p₂ + βL_M²) / (2c_ζ(1 − εp₁ − ε²p₂)λ_min(M̄)) = (1 + O(ε)) k‖Σ‖_F²/(2(1 − δ)λ_min(M̄)).

For comparison, the averaged Langevin dynamics ζ̇ = −(1/k)M̄ζ + ΣḂ has the exact stationary covariance C̄ = kX with
M̄X + XM̄ = ΣΣᵀ, so E‖ζ‖² = k tr X, equal to (k/2)tr(M̄⁻¹ΣΣᵀ) when ΣΣᵀ commutes with M̄: the bound is within (1 + O(ε))/(1 − δ)
of it when M̄ is isotropic, and pays the condition number of M̄ in the trace otherwise (it replaces every 1/λ_i by 1/λ_min),
as a scalar Lyapunov function must.

*Status.* Proved (Itô's formula and the pointwise inequality of Lemma E.0′). Checked (`exp_explicit_rate.py`, table 6;
`tests/test_excitation.py`) against the exact periodic stationary covariance of (LTV_τ), computed from the periodic Lyapunov
equation.

*Proof.* V(t, x) = ½xᵀ𝐇(t)x is C¹ in t and quadratic in x, so by Itô's formula along (LTV_τ), with Σ̃ := [0; Σ] ∈ ℝ^{(n+q)×q′},

    dV = [∂_tV + ∇_xV·A(t)x] dt + ∇_xV·Σ̃ dB_t + ½ tr(Σ̃ᵀ𝐇(t)Σ̃) dt ≤ [−γV + ν(t)] dt + ∇_xV·Σ̃ dB_t,

the inequality being Lemma E.0′'s V̇ ≤ −γV pointwise in (t, x), and Σ̃ᵀ𝐇Σ̃ = Σᵀ𝐇_ζζΣ with 𝐇_ζζ = I + εP₁ + ε²P₂ + βLᵀL from the
block form of 𝐇 (Corollary E.0′.2). The stochastic integral is a martingale because its integrand is linear in x and the linear
SDE has bounded second moments (Gaussian with a bounded covariance), so d/dt E V ≤ −γ E V + ν̄ and E V(t) ≤ e^{−γt}E V(0) + ν̄/γ. The bounds on E‖ζ‖² and
E‖w‖² follow from V ≥ ½(1 − εp₁ − ε²p₂)‖ζ‖² and V ≥ (β/2)‖w‖²; the estimate of ν̄ uses ‖P₁‖ ≤ p₁, ‖P₂‖ ≤ p₂, ‖L‖ ≤ L_M. The
averaged covariance is the stationary solution of the Lyapunov equation of the time-invariant averaged system. ∎

*Numerics* (`exp_explicit_rate.py` → `results/explicit_rate/table.md`, table 6). Σ = (2·10⁻³)^{1/2}I; the exact stationary
covariance from the periodic Lyapunov equation (RK4 with k·dt ≤ ¼; the period mean and the maximum of E‖ζ‖² coincide to four
digits), the averaged value k·tr X, and the bound:

| regressor (κ(M̄)) | k | δ | exact E‖ζ‖² | averaged k·tr X | bound | bound/exact | 1/(1 − δ) |
|---|---|---|---|---|---|---|---|
| Theorem E ideal of table 4 (7.8·10³) | 300 | 0.2 | 1.509·10⁶ | 1.509·10⁶ | 2.64·10⁷ | 17.5 | 1.25 |
| | 400 | 0.2 / 0.1 | 2.012·10⁶ | 2.012·10⁶ | 3.51·10⁷ / 3.11·10⁷ | 17.4 / 15.5 | 1.25 / 1.11 |
| | 600 | 0.2 / 0.1 | 3.018·10⁶ | 3.018·10⁶ | 5.23·10⁷ / 4.65·10⁷ | 17.3 / 15.4 | 1.25 / 1.11 |
| Ψ = (sin t, cos t) (1) | 20 | 0.2 / 0.1 | 8.025·10⁻² | 8.000·10⁻² | 0.216 / 0.202 | 2.69 / 2.52 | 1.25 / 1.11 |
| | 50 | 0.2 / 0.1 | 2.001·10⁻¹ | 2.000·10⁻¹ | 0.315 / 0.282 | 1.58 / 1.41 | 1.25 / 1.11 |
| | 200 | 0.2 / 0.1 | 8.000·10⁻¹ | 8.000·10⁻¹ | 1.051 / 0.935 | 1.31 / 1.17 | 1.25 / 1.11 |

The exact stationary spread is the averaged Langevin value up to 3·10⁻⁴ (ideal) and 1.25/k² (sinusoid): on these instances,
where Σ is isotropic so that ΣΣᵀ commutes with M̄, the averaging prediction of the dithered loop's stationary covariance is
exact to O(k⁻²). When ΣΣᵀ does not commute with M̄ the deviation is O(1/k) (referee round 7 measured ±0.25/k on
two-dimensional sinusoids with anisotropic Σ, the bound still holding with ratios 2.7–5.1). On the isotropic sinusoid
the bound exceeds the exact value by a factor that falls toward 1/(1 − δ) (2.7, 1.6, 1.3 at k = 20, 50, 200 for δ = 0.2;
the excess is the O(1/k) of the factors (1 ± εp₁)… and βL_M²). On the ideal, where one direction is excited 7.8·10³ times
less than the others, the bound is 15–18 times the exact value: it charges every one of the q = 15 directions the price of
the worst, (q/λ_min)/Σ_i(1/λ_i) ≈ q = 15 (up to 1/(1 − δ)), the trace-level price of a scalar Lyapunov function. The exact
E‖e‖² is ‖Σ‖_F²/(2k) = tr(M̄X)/k to five digits on both instances (the error rides on its slow manifold e ≈ Ψζ/k, and
tr(M̄X) = ½‖Σ‖_F² by the trace of the Lyapunov equation), while the bound on E‖w‖² is coarse (factor 1/β): on the sinusoid the
exact stationary E‖w‖² = E‖e − Lζ‖², computed from the periodic covariance and the slaving map L = −εΨ + ε²Ψ̇ − ε³ΨM, is
more than three orders of magnitude below the bound (asserted in the test below).


*Checked* (`tests/test_excitation.py`): on Ψ = (sin t, cos t) with Σ = (2·10⁻³)^{1/2}I, the explicit bound dominates the exact
period mean and the exact maximum of E‖ζ‖² at k = 20, 50, 200 for δ = 0.1 and 0.2, its ratio to the exact value decreases
toward 1/(1 − δ) (within 6% at k = 200), the grid value of ν̄ lies below its explicit bound, and the exact value agrees with
the averaged k·tr X within 2/k² (`test_corollary_E0prime_3_dither_covariance_bound_against_the_periodic_lyapunov_solution`).

**Lemma E.0″ (the slow regime for non-periodic persistently exciting regressors).** Let Ψ: [0, ∞) → ℝ^{n×q} be C² with
‖Ψ‖ ≤ Ψ_M, ‖Ψ̇‖ ≤ Ψ_D, ‖Ψ̈‖ ≤ Ψ_DD, M := ΨᵀΨ, and fix a window T > 0. For a bounded continuous X let

    ⟨X⟩_T(t) := (1/T)∫₀^∞ e^{−s/T} X(t + s) ds    (the forward exponential average; T d/dt⟨X⟩_T = ⟨X⟩_T − X),

and M̄ := ⟨M⟩_T. Assume the exponentially weighted PE condition

    M̄(t) ⪰ λ I_q for all t ≥ 0, with λ > 0,                                                                   (PE_T)

which is equivalent to classical PE (∫_t^{t+T}M ⪰ αI for all t gives (PE_T) with λ ≥ e^{−1}α/T, and (PE_T) gives PE over
windows T′ > T log(Λ/λ) with level T(λ − e^{−T′/T}Λ); proof below). Let Λ := sup_t λ_max(M̄(t)) ≤ Ψ_M². Fix symmetric matrices
M̂, F̂ (any; the long-time means of M̄ and F̄ minimise the time-averaged squares of ‖M̄ − M̂‖, ‖F̄ − F̂‖, not necessarily the
suprema below) and define

    P₁ := −2T(M̄ − M̂)   (so Ṗ₁ = 2(M − M̄)),    F := P₁M + MP₁,    F̄ := ⟨F⟩_T,    S := −T(F̄ − F̂)   (Ṡ = F − F̄),
    P₂ := −(M − M̄) + S,    N̄₂ := ½(F̄ − Ṁ̄) = ½(F̄ − (M̄ − M)/T),

and e₂, w, V exactly as in Lemma E.0′ (V = ½ζᵀ(I + εP₁ + ε²P₂)ζ + (β/2)‖w‖², ε = 1/k). The constants are those of Lemma E.0′ with
the period mean replaced by the exponential average, all of them computable from the regressor: the time-local peak factors

    ρ := sup_{t,v} ‖Ψ(t)v‖²/vᵀM̄(t)v,   ρ_D := sup ‖Ψ̇(t)v‖²/vᵀM̄(t)v,   ρ_DD := sup ‖Ψ̈(t)v‖²/vᵀM̄(t)v   (all ≤ Ψ_·²/λ),

the operator and per-direction norms

    p₁ := sup_t‖P₁(t)‖ = 2T sup‖M̄ − M̂‖,   p₂ := sup_t‖P₂(t)‖ ≤ Ψ_M² + T sup‖F̄ − F̂‖,
    π₁ := sup_{t,v} ‖P₁(t)v‖/‖v‖_{M̄(t)} ≤ 2T(√Λ + ‖M̂‖/√λ),   π₂ := sup_{t,v} ‖P₂(t)v‖/‖v‖_{M̄(t)} ≤ Ψ_M√ρ + √Λ + T(2p₁Ψ_M² + ‖F̂‖)/√λ,
    n₂ := sup_t‖N̄₂(t)‖ ≤ ½(sup‖F‖ + 2Ψ_M²/T) ≤ p₁Ψ_M² + Ψ_M²/T,

and, with these π₁, π₂, p₁, p₂, the quantities c₁, c₂, d_w, r_ζ and β := 2c₁²/δ of Lemma E.0′ (same formulas), together with

    δ₂ := ε n₂/λ,    c_ζ := 1 − δ − δ₂ − r_ζ/ε.

Fix δ ∈ (0, 1). If (A) c₁c₂ ≤ δ/2, (B) d_w ≤ k/4, (C) c_ζ > 0 and (D) εp₁ + ε²p₂ < 1, then along (LTV), for all (e, ζ, t),

    V̇ ≤ −ε c_ζ ζᵀM̄(t)ζ − (βk/4)‖w‖² ≤ −γV,    γ := min{2c_ζλ/(k(1 + εp₁ + ε²p₂)), k/2},

and ‖(e(t), ζ(t))‖ ≤ C e^{−γ(t−t₀)/2}‖(e(t₀), ζ(t₀))‖ with the C of Lemma E.0′. In the slow regime the guaranteed rate is

    γ/2 = (1 − δ − δ₂)(1 − εp₁ + O(k⁻²)) λ/k,    λ = inf_t λ_min(⟨M⟩_T(t)):

the exponentially weighted PE level of the window divided by k, up to the factor 1 − δ − δ₂, where δ₂ = O(1/k) is the price of
the residual second-order term (zero for a periodic Ψ with M̂ = ⟨M⟩ in the limit T → ∞, where Lemma E.0′ applies directly).

*Status.* Proved (the proof of Lemma E.0′ with the changes listed below). Checked (`exp_pe_general.py`,
`tests/test_excitation.py`) on a quasi-periodic regressor and on a regressor driven by a Lorenz orbit, where no periodic
construction applies: the conditions hold from k = 20 on for the quasi-periodic regressor and from k = 100 on for the
chaotic one (the derivative peak factors of a chaotic orbit are large), the certificate sup V̇/V ≤ −γ holds on the time grid
wherever they hold, and the guaranteed rate is within a factor 1.3–2.5 of the exact decay rate on the quasi-periodic
regressor (k = 50–200) and, with the window T = 20, within a factor 3.5–5 on the chaotic one (shorter windows pay more:
7.5–11 at T = 10 and 14–138 at T = 5), where the worst exponential window is three times less exciting than the long-time
average. Referee round 7 (`review/REVIEW_round7_E0p3_E0dp.md`) confirmed the proof and the code and found no counterexample
on an attack set of barely exciting, bursty, near-resonant and chaotic regressors.

*Proof.* Follow the proof of Lemma E.0′ with M̄, P₁, S, P₂ as defined here. (a) is unchanged (it does not involve the
averages). In (b), Ḣ = 2ε(M − M̄) + ε²(F − F̄ − Ṁ + Ṁ̄) exactly, by Ṗ₁ = 2(M − M̄) and Ṗ₂ = −Ṁ + Ṃ̄ + Ṡ = −Ṁ + Ṁ̄ + F − F̄;
the order-ε terms give −ζᵀM̄ζ as before, and the order-ε² terms give −ζᵀP₁Mζ + ½ζᵀṀζ + ½ζᵀ(F − F̄ − Ṁ + Ṁ̄)ζ = ½ζᵀ(Ṁ̄ − F̄)ζ =
−ζᵀN̄₂ζ, since ζᵀP₁Mζ = ½ζᵀFζ: the periodic cancellation N₂ = 0 is replaced by the residual −ε²ζᵀN̄₂ζ, and R_ζ has the same
form, so d/dt[½ζᵀHζ] = −εζᵀ(M̄ + εN̄₂)ζ + R_ζ + ζᵀHΨᵀw. In (c), every per-direction bound is taken in the time-local norm
‖v‖_{M̄(t)}: ‖Ψ(t)v‖ ≤ √ρ‖v‖_{M̄(t)} and likewise for Ψ̇, Ψ̈ by the definition of the peak factors; ‖M(t)v‖ ≤ Ψ_M√ρ‖v‖_{M̄(t)};
‖M̄(t)v‖ ≤ √Λ‖v‖_{M̄(t)}; ‖P₁(t)v‖ = 2T‖(M̄(t) − M̂)v‖ ≤ 2T(√Λ + ‖M̂‖/√λ)‖v‖_{M̄(t)}, using √λ‖v‖ ≤ ‖v‖_{M̄(t)} (this is where the
exponential average pays off: P₁(t) is a function of M̄(t) alone, with no integral over other times); for F̄(t)v and S(t)v,
which do involve other times r ≥ t, ‖F(r)v‖ ≤ ‖P₁(r)‖‖M(r)v‖ + ‖M(r)‖‖P₁(r)v‖ ≤ 2p₁Ψ_M²‖v‖ ≤ (2p₁Ψ_M²/√λ)‖v‖_{M̄(t)}, an average
of such bounds bounds ‖F̄(t)v‖, and ‖S(t)v‖ ≤ T(‖F̄(t)v‖ + ‖F̂v‖), whence the bound on π₂; the stated bounds on p₂ and n₂ are
‖P₂‖ ≤ ‖M − M̄‖ + ‖S‖ and ‖N̄₂‖ ≤ ½(‖F̄‖ + ‖Ṁ̄‖) with ‖Ṁ̄‖ = ‖M̄ − M‖/T ≤ 2Ψ_M²/T and ‖F̄‖ ≤ sup‖F‖ ≤ 2p₁Ψ_M². The bounds |R_ζ| ≤
r_ζ‖ζ‖²_{M̄}, |ζᵀHΨᵀw| ≤ c₁‖ζ‖_{M̄}‖w‖ and βwᵀẇ ≤ β[−(k − d_w)‖w‖² + c₂‖ζ‖_{M̄}‖w‖] follow with the same formulas (the suprema
over v at the current time t of the norms that appear). In (d), the extra term is bounded by ε²|ζᵀN̄₂ζ| ≤ ε²n₂‖ζ‖² ≤
(ε²n₂/λ)ζᵀM̄ζ = εδ₂ζᵀM̄ζ, which turns 1 − δ − r_ζ/ε into c_ζ = 1 − δ − δ₂ − r_ζ/ε; the rest of the assembly, the sandwich
½(1 − εp₁ − ε²p₂)‖ζ‖² + (β/2)‖w‖² ≤ V ≤ ½(1 + εp₁ + ε²p₂)‖ζ‖² + (β/2)‖w‖² and ζᵀM̄(t)ζ ≥ λ‖ζ‖² give V̇ ≤ −γV. The equivalence
of (PE_T) with classical PE: if ∫_t^{t+T}M ⪰ αI then M̄(t) ⪰ (1/T)∫₀^T e^{−s/T}M(t+s)ds ⪰ (e^{−1}/T)∫_t^{t+T}M ⪰ (e^{−1}α/T)I;
conversely (PE_T) gives, for every T′ > 0, ∫_t^{t+T′}M ⪰ ∫_t^{t+T′}M e^{−(r−t)/T}dr ⪰ T(M̄(t) − e^{−T′/T}M̄(t + T′)) ⪰ T(λ − e^{−T′/T}Λ)I,
positive for T′ > T log(Λ/λ) (the identity ∫_t^{t+T′}e^{−(r−t)/T}M(r)dr = T[M̄(t) − e^{−T′/T}M̄(t+T′)] follows from the
definition of the exponential average). ∎

*Remark E.0″.1 (reading).*

- *What changes and what does not.* The construction, the conditions and the rate formula are those of Lemma E.0′; the
  period mean becomes the exponential average over a window, the peak factors become time-local, and one new term appears:
  the residual N̄₂, whose price δ₂ = εn₂/λ vanishes like 1/k. For a periodic regressor the exponential average with a finite
  window keeps residual oscillations of relative size O(period/T), so Lemma E.0′ is sharper there; Lemma E.0″ is the version
  for references that are not periodic — quasi-periodic (several incommensurate frequencies) or chaotic (a reference along a
  strange attractor, as in the Lorenz-tracking experiments of Part 1), where only PE is available.
- *The window is a design parameter.* The lemma holds for every T > 0; the guaranteed rate is (1 − δ − δ₂)λ_T/k with
  λ_T = inf λ_min(⟨M⟩_T), which increases toward the long-time excitation level as T grows, against p₁ ∝ T, n₂ ∝ T and the
  threshold in k, which grow with T. The numerics show the optimum at a window of a few characteristic times of the regressor.
- *Exact rate.* For a non-periodic regressor the exact rate is the top Lyapunov exponent of (LTV), which the certificate
  does not need: the lemma is a pointwise inequality in (t, e, ζ), checked on a grid in t by the same generalised eigenvalue
  problem as Lemma E.0′.

*Numerics* (`exp_pe_general.py` → `results/explicit_rate/general_pe.md`). Two regressors that are not periodic, on [0, 300]
(grid step 0.005; the exponential averages formed from an extended grid; the certificate sup V̇/V over 400 times and all
states; the exact rate the top Lyapunov exponent of (LTV) over 200 time units), with the measured constants (the grid
suprema of ‖P₁‖, ‖P₂‖, the per-direction norms and ‖N̄₂‖, computable constants of the regressor like the peak factors).

*Quasi-periodic Ψ = (sin t, cos √2 t)* (q = 2; the long-time mean of M is ½I, and the exact decay rate is λ_∞/k = 0.5/k to
within 1%, 0.6% at k = 20 and 0.8% at k = 400): the window average has λ_T = 0.297, 0.391, 0.444 for T = 2π, 4π, 8π (ρ = 4.0, ρ_D ≈ 6.1–6.4, ρ_DD ≈ 10.0, p₁ ≈ 2.9,
p₂ = 8.0–9.2, π₁ = 4.3–4.7, π₂ = 11.6–13.0, ‖N̄₂‖ = 0.71, 0.35, 0.17); the conditions hold from k = 20 on for both δ, the
certificate holds at every row, and

| T | k | δ₂ | guaranteed γ/2 (δ = 0.2 / 0.1) | exact | exact/guaranteed | λ_T/k |
|---|---|---|---|---|---|---|
| 4π | 20 | 0.045 | 0.0100 / 0.0117 | 0.0249 | 2.49 / 2.13 | 0.0195 |
| 4π | 50 | 0.018 | 0.00557 / 0.00631 | 0.00995 | 1.78 / 1.58 | 0.00781 |
| 4π | 200 | 0.004 | 0.00153 / 0.00172 | 0.00249 | 1.63 / 1.45 | 0.00195 |
| 8π | 20 | 0.020 | 0.0118 / 0.0137 | 0.0249 | 2.10 / 1.81 | 0.0222 |
| 8π | 200 | 0.002 | 0.00174 / 0.00196 | 0.00249 | 1.43 / 1.27 | 0.00222 |

(T = 2π: 3.69 → 2.16 for δ = 0.2 and 3.10 → 1.92 for δ = 0.1 from k = 20 to 200.) The ratio exact/guaranteed is the product of the window loss
λ_∞/λ_T (1.69, 1.28, 1.13 for the three windows), the factor 1/(1 − δ − δ₂) and the O(1/k) factors; the explicit closed-form
constants satisfy the conditions from k = 50 on and give rates below the measured-constant ones by 22–79% at k = 50 (the
longer the window, the larger the deficit), 5–19% at k = 100, at most 5% at k = 200 and about 1% at k = 400. All constants
are suprema over [0, 300] (the exponential averages are formed on an extended grid): on the quasi-periodic regressor they are
identical to four digits on [0, 10⁴] (closed form); on the Lorenz regressor below, the suprema over [0, 1500] are larger by
up to 14% (p₁ 5.31 → 5.67, p₂ 19.4 → 20.7, π₂ 49.4 → 51.9, ‖N̄₂‖ 0.42 → 0.48 at T = 20), which lowers the guaranteed rate
at k = 100 by 3% and keeps the certificate (referee round 7): the constants of a chaotic regressor are those of the stretch
of the attractor that the run visits.

*Chaotic Ψ = (x₁/20, (x₃ − 25)/20) along a Lorenz orbit* (σ = 10, r = 28, b = 8/3, traversed five times slower than standard
Lorenz time, Ψ̇ by the chain rule through the vector field; q = 2): the exact decay rates are 0.00759, 0.00304, 0.00152,
0.00076, 0.00038 at k = 20, 50, 100, 200, 400 (an effective excitation level 0.15), while the worst exponential window has
λ_T = 0.015, 0.024, 0.051 for T = 5, 10, 20 (Λ = 0.55, 0.36, 0.24; ρ ≈ 20–26, ρ_D ≈ 110–160, ρ_DD ≈ 2000–3100; ‖N̄₂‖ = 1.05,
0.34, 0.42: the residual does not follow a monotone rule in T). The conditions hold from k = 100 on (condition (A) fails below: the derivative peak factors of a chaotic orbit
are large), the certificate holds at every admissible row, and at T = 20 the guaranteed rate is 0.00030 / 0.00035 (k = 100),
0.00018 / 0.00021 (k = 200) against the exact 0.00152 and 0.00076: exact/guaranteed = 5.0 / 4.3 and 4.2 / 3.7, of which the
window loss λ_∞/λ_T ≈ 3 is the dominant factor. For a chaotic reference the uniform-in-time construction pays the least
exciting stretch of the attractor; a longer window (T = 20 against the slowed Lorenz loop time ≈ 3.75) recovers part of it
(λ_T = 0.015 → 0.051 from T = 5 to 20) at the price of a larger p₂ and π₂ (9.9 → 19.4 and 13.7 → 49.4). An exploratory
regressor with components tanh(x_i/8) traversed at standard Lorenz speed (not in the table) had ρ_DD of order 10⁶–10⁷ and a
threshold beyond k = 400: the slow regime requires the
reference to be slow relative to the loop, and the time scale of a chaotic reference is the designer's choice.


*Checked* (`tests/test_excitation.py`): on Ψ = (sin t, cos √2 t) with T = 4π, the differential identities TṀ̄ = M̄ − M and Ṗ₁ =
2(M − M̄) on the grid (10⁻⁴), the conditions at k = 20, 50, 200 for δ = 0.1, 0.2 with δ₂ ≤ 5%, the certificate sup V̇/V ≤ −γ on
400 times, and the guaranteed rate below the exact decay rate (top Lyapunov exponent over 200 time units) within a factor 2.6,
the ratio improving with k (`test_lemma_E0dprime_quasi_periodic_regressor_certificate_and_rate`).

**Theorem E (with explicit constants).** Let Ψ_G := m^{−1/2}Ψ(I_m ⊗ G^{1/2}) ∈ ℝ^{n×mp}, and let (α_G, T, Ψ_{G,M}, Ψ_{G,D})
be PE data for it: ∫_t^{t+T}Ψ_GᵀΨ_G ⪰ α_G I, ‖Ψ_G‖ ≤ Ψ_{G,M}, ‖Ψ̇_G‖ ≤ Ψ_{G,D}. One may take α_G := αλ_min(G)/m,
Ψ_{G,M} := m^{−1/2}λ_G^{1/2}sup‖Ψ‖ and Ψ_{G,D} := m^{−1/2}λ_G^{1/2}sup‖Ψ̇‖; Ψ̇ is bounded because x_m ∈ C¹ with bounded
derivative and the Gaussian's mixed derivatives are bounded for s ≥ s_min. Let (ε, γ, p_lo, p_hi, g_V) be the constants of
Lemma E.0 for the data (k_e, α_G, T, Ψ_{G,M}, Ψ_{G,D}), with any admissible ε (for example the closed form of E.0(ii)), where
g_V := 1 + ε(2Q_{G,M} + TΨ_{G,M}) bounds ‖∂V/∂w‖/‖w‖ (the Hessian of V is I plus a matrix of norm at most
εTΨ_{G,M} + 2εQ_{G,M}) and Q_{G,M} := T²Ψ_{G,M}²/2. Let r := min(ρ, δ), let L_J be a Lipschitz constant, uniform in t, of the
Jacobian of the closed-loop field F of step 2 on the ball D, and let S := diag(I_n, m^{−1/2}(I_m ⊗ G^{−1/2})) be the
similarity to the coordinates w := (e, ζ) of step 3. Put

    C̃_R := ‖S‖(L_J/2)‖S⁻¹‖²,    r_w := min{ r/‖S⁻¹‖, γ/(2g_V C̃_R) },    r₀ := r_w (p_lo/p_hi)^{1/2},
    λ := γ/(4p_hi),    c₀ := r₀²/2,    C := (2p_hi/p_lo)^{1/2}  (≤ √6 for the closed-form ε).

If V_OT(t₀) ≤ c₀, then:

- the optimal assignment π at t₀ is unique, and it is the unique optimal assignment for all large t;
- ‖e(t)‖ + W(μ_t,μ*) ≤ Ce^{−λ(t−t₀)}(‖e(t₀)‖ + W(μ_{t₀},μ*)) for t ≥ t₀;
- ϑ_i(t) → ϑ*_{π(i)}. The labelled error E_lab(t) := ((1/m)Σ‖ϑ_i(t) − ϑ*_i‖²_{G⁻¹})^{1/2} tends to
  ((1/m)Σ‖ϑ*_{π(i)} − ϑ*_i‖²_{G⁻¹})^{1/2}, which is non-zero for π ≠ id.

Every constant is an explicit function of (k_e, α_G, T, Ψ_{G,M}, Ψ_{G,D}, m, G, r, L_J), hence uniform over the class
of regressors with given PE data (T, α, sup‖Ψ‖, sup‖Ψ̇‖). With the closed-form ε, λ ≥ (1/3)·min{k_e/4, εα_G/2}. The
constants are conservative, in a quantified way (Remark E.0a); their role is to make the statement self-contained. L_J
is defined as a supremum; a closed-form bound follows from the third-derivative suprema of the Gaussian and is not
computed here (`excitation.jacobian_lipschitz_estimate` estimates it from samples for the tests). Version 3 proved this
theorem with existence constants, conditional on Morgan–Narendra [S8, Thm 1]; Lemma E.0 removes that citation.

*Proof.*

1. *Symmetry.* Relabelling the network, (e, ϑ) ↦ (e, ϑ∘σ), maps solutions to solutions. It fixes e, is an isometry of
   the matched metric, and conjugates the dynamics near (0, ϑ*∘π) to those near (0, ϑ*). So it suffices to treat
   π = id.

2. *Domain and linearisation.* Let r := min(ρ, δ) and D := {z = (e, ϑ̃) ∈ ℝ^{n+mp} : ‖z‖₂ < r}.

   - For z ∈ D we have ‖e‖ < ρ, so f = f̂*.
   - Also ‖ϑ̃_i‖ < δ for every i, so ϑ_i ∈ int Θ_K^{ε} and Proj_K = id.

   On [0,∞) × D the closed loop is F(t,z) = (−k_e e + f̂*(x) − f̂(x), (G∇_ϑΦ(x;ϑ_i)ᵀe)_i), with x = x_m(t) + e.

   - It is C¹ in (t,z), since x_m ∈ C¹, and C^∞ in z.
   - Its derivatives through order two are bounded uniformly in t, because the Gaussian's derivatives through order
     three are bounded uniformly in x for s ≥ s_min.
   - Its Jacobian at 0 is A(t) = [[−k_eI, −(1/m)Ψ], [(I_m⊗G)Ψᵀ, 0]], and the remainder is O(‖z‖²).
   - The Jacobian and the identity QA + AᵀQ = diag(−k_eI, 0), with Q = diag(I/2, (I⊗G⁻¹)/(2m)), are checked
     numerically by `test_theorem_E_linearisation_and_lyapunov_identity`.

3. *LTV stability, with explicit constants.* Put w := (e, ζ) with ζ := m^{−1/2}(I_m⊗G^{−1/2})ϑ̃, that is w = Sz, and
   Ψ_G := m^{−1/2}Ψ(I⊗G^{1/2}).

   - The system becomes ė = −k_ee − Ψ_Gζ and ζ̇ = Ψ_Gᵀe: this is (LTV) of Lemma E.0 with k = k_e and Ψ = Ψ_G, that is,
     ẇ = [A₀ + B(t)]w with A₀ = diag(−k_eI, 0) and B skew-symmetric, bounded, with bounded derivative.
   - The excitation transfers: ∫Ψ_GᵀΨ_G ⪰ (αλ_min(G)/m)I = α_G I, and ‖Ψ_G‖ ≤ Ψ_{G,M}, ‖Ψ̇_G‖ ≤ Ψ_{G,D}.
   - Lemma E.0 gives the strict Lyapunov function V(t, w) with p_lo‖w‖² ≤ V ≤ p_hi‖w‖² and V̇ ≤ −γ‖w‖² along the linear
     system, hence ‖w(t)‖ ≤ (p_hi/p_lo)^{1/2}e^{−(γ/(2p_hi))(t−t₀)}‖w(t₀)‖: uniform exponential stability with explicit
     overshoot and rate.
   - The w-coordinates are the coordinates of the matched Lyapunov function: ‖ζ‖² = (1/m)Σ‖ϑ̃_i‖²_{G⁻¹} = W_π², so
     ‖w‖² = ‖e‖² + W_π² = 2V_π.
   - The bound on Ψ̇ cannot be dropped. Ψ = sin t² is bounded and PE, but for ė = −e − Ψζ, ζ̇ = Ψe, from w(t₀) = (e, ζ) = (0, 1), the ratio
     ‖w(t₀+200)‖/‖w(t₀)‖ is 0.9715, 0.9942 and 0.9980 for t₀ = 5, 20 and 50 (review round 2), and Lemma E.0's proof
     step (f) shows that the decay over a fixed horizon tends to 0 as t₀ → ∞. So that system is not UES. Consistently,
     the rate of Lemma E.0 tends to 0 as Ψ_D → ∞.

4. *Nonlinear system, with an explicit ball.* On the closure of D (to which the Lipschitz bound extends) the field is
   F(t,z) = A(t)z + R(t,z) with ‖R(t,z)‖ ≤ (L_J/2)‖z‖² (Taylor with the Lipschitz Jacobian). In w-coordinates, ẇ = [A₀ + B(t)]w + R̃(t,w) with R̃(t,w) := S R(t, S⁻¹w), so
   ‖R̃(t,w)‖ ≤ C̃_R‖w‖² on {‖w‖ < r/‖S⁻¹‖} ⊂ S(D). Along the nonlinear system, since ‖∂V/∂w‖ ≤ g_V‖w‖,

       V̇ ≤ −γ‖w‖² + g_V C̃_R ‖w‖³ ≤ −(γ/2)‖w‖²    whenever ‖w‖ ≤ r_w.

   The sublevel set {V ≤ p_lo r_w²} is contained in {‖w‖ ≤ r_w}, so it is forward invariant, and it contains the ball
   {‖w‖ ≤ r₀}, since V ≤ p_hi‖w‖². On it V̇ ≤ −(γ/(2p_hi))V, so for ‖w(t₀)‖ ≤ r₀,

       ‖w(t)‖ ≤ (p_hi/p_lo)^{1/2} e^{−λ(t−t₀)} ‖w(t₀)‖,    λ = γ/(4p_hi).

   These trajectories stay in D, where the closed loop equals F. (Khalil Thm 4.13 gives the same conclusion with
   existence constants; the explicit form is what Lemma E.0 adds.)

5. *Transfer to (‖e‖, W).* Let π be optimal at t₀. Then ‖w(t₀)‖² = 2V_π(t₀) = 2V_OT(t₀) ≤ 2c₀ = r₀², so step 4 applies.
   Moreover ‖e‖ + W ≤ ‖e‖ + W_π = ‖e‖ + ‖ζ‖ ≤ √2‖w‖ at every t, and ‖w(t₀)‖ ≤ ‖e(t₀)‖ + ‖ζ(t₀)‖ = ‖e(t₀)‖ + W(t₀). This
   gives the estimate with C = (2p_hi/p_lo)^{1/2}.

   *Uniqueness.* Suppose V_π, V_π' ≤ c₀ with π ≠ π'. Then the state lies in the invariant balls of step 4 around both
   (0, ϑ*∘π) and (0, ϑ*∘π'). The solution would converge to two distinct equilibria, which is impossible.

   *Limits.* Hence z(t) → 0, that is, ϑ_i(t) → ϑ*_{π(i)}. For π' ≠ π,
   V_π'(t) → (1/2m)Σ‖ϑ*_{π(i)} − ϑ*_{π'(i)}‖²_{G⁻¹} ≥ d²/2 > 0, while V_π(t) → 0. So π is the unique optimal assignment
   for all large t, and E_lab converges as stated. ∎

*Remark E.1 (the rate).*

- **What converges.** "The parameters converge" is false for the labelled error and true for the matched one.
- **The rate is not set by the PE level.** λ > 0 exists for each admissible Ψ, but it depends on Ψ(·), not only on
  its PE level. No positive function of αλ_min(G)/m, k_e, T and sup‖Ψ_G‖ alone bounds it below.
  - For m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = I, k_e = 3 and x_m = 1.5 sin ωt (T = 2π), these four quantities do not
    depend on ω. The exact Floquet rate falls from 4.7e-2 (ω = 1) to 2.9e-5 (ω = 64), like ω⁻² (review round 3).
  - A lower bound must also involve sup‖Ψ̇_G‖ (step 3). Lemma E.0 gives one; it reproduces the ω⁻² law (Remark E.0a).
  - Exact Floquet rates (review round 2) are also not monotone in G or in k_e, and they decrease for large G. When k_e is large relative to the reference bandwidth, they
  approach λ_min⟨Ψ_GᵀΨ_G⟩/k_e.
- **The numbers** (`results/pe_convergence`; m = 3, n = 2, six runs):

  | quantity | range |
  |---|---|
  | min over sliding 60 s windows of λ_min(∫_t^{t+60}ΨᵀΨ) | 5.9e-4 to 2.7e-2 |
  | fitted decay rate of the matched error, second half of the run | 2.3e-3 to 4.3e-3 s⁻¹ |
  | rank correlation between the two | 0.60 |

The convergence test is therefore qualitative. The machine-checked part of the proof is the linearisation.

**Corollary E.2 (uniform-in-time exponential merging under PE).** Take two solutions of the same closed loop with
V_OT(0) ≤ c₀ and V'_OT(0) ≤ c₀. Then

    ‖x_t − x'_t‖ + W(μ_t, μ'_t) ≤ Ce^{−λt}R₀,    R₀ := ‖e(0)‖ + ‖e'(0)‖ + W(μ_0,μ*) + W(μ'_0,μ*).

Let δ₀ := ‖x_0 − x'_0‖ + W(μ_0, μ'_0). Combine this with the coupling estimate of Theorem D(c) on [0, T₀], with
T₀ := (Λ + λ)^{−1}log(CR₀/δ₀).

- Both solutions stay in D (proof step 4), so ‖x − x_m(t)‖ < r ≤ ρ.
- The segment between x_t and x'_t lies in B(x_m(t), ρ), where f = f̂*, and Proj_K = id.
- Hence (c) holds with M_x := k_e + 2a_max J_c, with L_x and L_ϑ computed for Ē = r and an inactive projection, and
  with L_2 = λ_G^{1/2}L₁. This Λ depends neither on T₀ nor on f outside the tube.

This gives uniform-in-time Hölder continuity in the basin:

    sup_{t≥0}(‖x_t − x'_t‖ + W(μ_t, μ'_t)) ≤ (CR₀)^{Λ/(Λ+λ)} δ₀^{λ/(Λ+λ)}.

*Proof.* ‖x_t − x'_t‖ = ‖e − e'‖ ≤ ‖e‖ + ‖e'‖, and W(μ_t, μ'_t) ≤ W(μ_t, μ*) + W(μ'_t, μ*). Bound each bracket by
Theorem E. On [0, T₀] use e^{ΛT₀}δ₀; on [T₀, ∞) use CR₀e^{−λT₀}. Note that δ₀ ≤ R₀ ≤ CR₀. ∎

*Remark E.3 (the collapse point).* Take the G⁻¹-closest pair (j,k). Set e = 0 and ϑ_j = ϑ_k = (ϑ*_j + ϑ*_k)/2, with
all other particles ideal.

- At this point V_OT = d²/8, and there are two optimal assignments.
- By exchange symmetry, ϑ_j ≡ ϑ_k for all t, so W(μ_t, μ*) ≥ d/2 forever.

The uniqueness argument of step 5 shows c₀ < d²/8 automatically, and d²/8 is sharp along this direction.

---

## 6b. Proposition G — the slow geometry: joint descent, and variable projection in the induced metric

**Setting.**

- f = f(x) is time-invariant and C².
- x_m is T₀-periodic, with ω_ref := 2π/T₀ and occupation measure ρ := (1/T₀)∫₀^{T₀}δ_{x_m(t)}dt.
- Law (L) acts in the interior, with gains (g_a, δγ_c, δγ_s). Put Γ := diag(γ_c I_n, γ_s), acting per kernel, and
  G_ξ := δΓ. Write the geometry as ξ := (c, s).
- L̄(a,c,s) := ½∫‖f − (1/m)Σa_iφ(·;c_i,s_i)‖²dρ, and M(c,s) := (1/m²)∫φφᵀdρ ⊗ I_n is the Hessian of L̄ in a.
- 𝒦 is the set of geometries ξ on which (H1) holds, and 𝒦′ ⋐ 𝒦 is compact. (𝒦 is not the functional 𝒢(μ) of §6c.)
- a*(ξ) := argmin_a L̄, and L*(ξ) := L̄(a*(ξ), ξ) is the profile (variable-projection) loss.
- Three dimensionless parameters:
  - the rate ratio R := g_a m₀/(δ·sup_{𝒦′}‖Γ^{1/2}∇²L*Γ^{1/2}‖);
  - the slope S := (δ/g_a)·sup_{𝒦′}‖D_ξa*·Γ^{1/2}‖²;
  - the averaging parameter ε_avg := (m/k_e)λ_max(G^{1/2}∇²L̄G^{1/2})/ω_eff, the second ratio in (H2a′). The tables
    report its inverse, with ω_ref in place of ω_eff.

**Hypotheses.**

- **(H0)** On the interval considered, the closed-loop and reduced geometries stay in 𝒦′. The parameters stay in
  int Θ_K^{ε}, and sup_{𝒦′} max_i‖a*_i‖ < a_max − ε_K.
- **(H1)** M ⪰ m₀I on 𝒦.
- **(H2a′)** *Two scales.* k_e ≫ ω_eff and ε_avg ≪ 1, where ω_eff := max{ω_ref, sup_t‖ẋ_m‖/s_min, sup_X‖∂_x(f − f̂)‖,
  sup_t‖(d/dt)r(x_m(t))‖/sup‖r‖}. The last entry is the bandwidth of the residual along the reference. Also
  sup‖r‖/k_e ≪ s_min.
- **(H2b)** *Rate separation.* R ≫ 1: the slowest weight rate g_a m₀ m/k_e is much larger than the geometry rate
  (mδ/k_e)‖Γ^{1/2}∇²L*Γ^{1/2}‖.
- **(H2c)** *Flat slow manifold.* S ≪ 1. Since Da* = −M⁻¹∂²_{a,ξ}L̄, a checkable sufficient condition is
  δ·max(γ_c, γ_s)·‖∂²_{a,ξ}L̄‖² ≪ g_a m₀².

Three remarks on the hypotheses.

- **(H2b) is about λ_min.** It is a condition on λ_min of the feature Gram matrix. Since 0 < φ ≤ 1, m₀ ≤ 1/m². So
  the slowest weight rate is at most g_a/(k_e m), and overlapping kernels make (H2b) fail.
- **The ratios are gain choices.** R and S contain neither k_e nor ω_ref. The separation between weights and geometry
  is therefore set by the gains.
- **(H2b) does not imply (H2c).** The product R·S = m₀·sup‖Da*Γ^{1/2}‖²/sup‖Γ^{1/2}∇²L*Γ^{1/2}‖ does not depend on g_a,
  δ or k_e. It is fixed by the problem and by the shape of Γ: in the slope configuration below it ranges from 166 to
  1149 over γ_s/γ_c ∈ [10⁻⁴, 10⁴] (review round 3). It is large whenever refitting the weights largely compensates a
  change of geometry.

**Proposition G (formal reduction).**

(a) *Two scales, under (H0) and (H2a′).* To leading order in ω_eff/k_e and ε_avg, the period-averaged closed loop is
the joint gradient flow ϑ̇_i = −(m/k_e)G∇_{ϑ_i}L̄(ϑ_1..ϑ_m). Its measure form is −(1/k_e)G∇(δL̄/δμ).

- This is the drift of (𝔄) at σ = τ = 0, with Φ in place of Φ̃ and (a,c,s) in place of (α,γ,κ), by the same
  averaging computation as Theorem H(i).
- It is not literally a special case of H(i): the saturated coordinates there induce a different, state-dependent
  metric.
- L̄ is non-increasing along the averaged flow. Along the closed loop, L̄ at period-averaged parameters decreases
  while the averaged dissipation (m/k_e)‖∇L̄‖²_G dominates the O(ω_eff/k_e) corrections of Proposition H.7 and
  Remark H.7a. Near equilibrium it can increase by O((ω_eff/k_e)²), as the weights settle at their shifted
  equilibrium (review round 3: +8·10⁻⁷ over 35 periods at k_e = 20, and 3.7·10⁻⁸, about 22 times less, at
  k_e = 80).

(b′) *Three scales, for the averaged system (a), under (H0), (H1), (H2a′) and (H2b).* After a boundary layer, the
weights follow the slow manifold:

    a − a*(ξ) = −(k_e/(mg_a))M⁻¹Da*·ξ̇ + h.o.t. = (δ/g_a)M⁻¹Da*(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*(ξ) + h.o.t.

This reduces to (δ/g_a)M⁻¹Da*Γ∇L*(ξ) when S ≪ 1. The geometry follows, to leading order in 1/R and for every S,

    ξ̇ = −(mδ/k_e)(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*(ξ).

- This is the gradient flow of (m/k_e)L* in the metric that the G⁻¹ metric induces on the graph {(a*(ξ), ξ)}.
- Its metric tensor is G_ξ⁻¹ + g_a⁻¹Da*ᵀDa*, and it couples all kernels through Da*.
- In the closed loop, the weight equilibrium is further shifted by O(ω_eff/k_e) (Proposition H.7). This shift can
  exceed the offset displayed above: in the `satisfying` run it is 3–6·10⁻³, against an offset of 4·10⁻⁵–3·10⁻⁴
  (review round 3).

(b) *Variable projection, under (H2c) in addition.* Let (c̄, s̄) solve ξ̄' = −(mδ/k_e)Γ∇L*(ξ̄), that is,
c̄_i' = −(mδγ_c/k_e)∇_{c_i}L* and s̄_i' = −(mδγ_s/k_e)∂_{s_i}L*, with (c̄,s̄)(0) = (c,s)(0). Then

    sup_{t ≤ τ₁k_e/(mδ)} ‖(c,s)(t) − (c̄,s̄)(t)‖ → 0   as (ω_eff/k_e, ε_avg, S, 1/R) → 0, with τ₁ fixed.

- Formally the error is O(ω_eff/k_e + ε_avg + S + 1/R).
- L* is non-increasing along ξ̄, and the weights track a*.
- This flow is the S → 0 limit of (b′).

*No mean-field form.* For fixed m, on configurations of m distinct geometries in 𝒦, the flow (b) is the gradient flow
of L*/k_e in the metric (1/m)Σ_i‖ξ̇_i‖²_{G_ξ⁻¹}. Locally, that metric is the Wasserstein metric restricted to uniform
m-atom measures.

- L* is the unregularised profile loss. It is the λ → 0 limit of λℒ^λ of [S2] (time rescaled by λ). It is not [S2]'s
  λ → 0 functional ℒ⁰, which is a min-norm interpolation functional whose flow is an ultra-fast diffusion.
- As a functional of the geometry marginal, L* depends only on the support, because the free weights absorb the masses.
  So it is not W-continuous, and the statement has no mean-field form. This is consistent with m₀ ≤ 1/m².

*Derivation (formal).* Each step is a standard theorem on its own time scale. Their composition on the O(k_e/(mδ))
horizon is formal.

1. *Error channel (exact identity).* For w := e − r(x,ϑ)/k_e, ẇ = −k_e w − ṙ/k_e. So
   ‖w(t)‖ ≤ e^{−k_et}‖w(0)‖ + sup‖ṙ‖/k_e², and e = r(x_m)/k_e + O(sup‖ṙ‖/k_e² + sup‖∂_x r‖·sup‖r‖/k_e²). This is
   O(ω_eff·sup‖r‖/k_e²), with the bandwidth term in ω_eff.
2. *Period averaging of the joint system* (Khalil Thm 10.4; Sastry–Bodson Ch. 4). This gives (a), on [0, b/ε_avg].
3. *Weights against geometry.* Use Tikhonov in slow time (Kokotović–Khalil–O'Reilly) with the uniformly exponentially
   stable averaged weights on 𝒦′ (H1).
   - On the slow manifold a = h(ξ), ȧ = Dh ξ̇ = −(mg_a/k_e)M(h − a*), so h − a* = −(k_e/(mg_a))M⁻¹Dh ξ̇. This must
     be evaluated with the (b′) velocity, not the VarPro velocity.
   - Since L̄ is quadratic in a and ∂²_{ξa}L̄ = −Da*ᵀM, we get ∇_ξL̄(h,ξ) = ∇L* − Da*ᵀM(h − a*) + O((δ/g_a)²).
   - Substituting into ξ̇ = −(mδ/k_e)Γ∇_ξL̄(h,ξ), with Dh = Da* + O(1/R), gives (b′). The resummed form solves the
     invariance equation exactly once Dh − Da* and the quadratic term are dropped, so the remaining corrections are
     in 1/R (they involve D²a* and ∂_ξM), not in S.
4. *Envelope theorem.* ∇_ξL̄(a*(ξ), ξ) = ∇L*(ξ). When S ≪ 1, the metric correction in (b′) drops out, which gives (b).

The composition requires the averaged system to remain valid up to t ~ τ₁k_e/(mδ), which is much longer than b/ε_avg.
Step 2 does not certify this. A proof would combine:

- a near-identity averaging transform;
- uniform exponential stability of the averaged weights on 𝒦′;
- Tikhonov in slow time with an O(ε_avg) perturbation;
- Grönwall in ξ (Sastry–Bodson §4.4; Sanders–Verhulst–Murdock Ch. 5).

*Numerics* (`results/variable_projection`):

- The first two rows use centred differences of period means.
- The slope rows use the displacement between the first and last period means. They start on the slow manifold of
  (b′), so no weight boundary layer enters the comparison.
- (b′) was evaluated only in the slope rows.
- k_e/ω_eff is 2.0 in the first two rows with the box's s_min, and 2.6–4.8 with the actual widths, so (H2a′) holds
  only marginally there. In the `satisfying` run the residual bandwidth, about 7.7, then dominates.
- k_e/ω_eff is about 16 in the slope rows.

| configuration | (H2b) ratio | S | G(a): closed loop against the joint flow (a) | G(b) and the induced-metric flow (b′) |
|---|---|---|---|---|
| `violating` (m = 4, overlapping) | 0.10–0.23 | 5.4–25 | joint cosine 0.9994 / 0.9985, size ratio 1.00 / 0.94 | (b) fails: cosine −0.12 / 0.39, and L* *increases* from 3.8e-3 to 5.9e-3 |
| `satisfying` (m = 2) | 79–87 | 0.006 | joint cosine 0.985 / 0.962, size 0.96 / 0.94; geometry block 0.999 / 0.998; the weight block fails (see below) | (b) holds: cosine 0.999, ratio 0.97; L* decreases from 1.20e-2 to 4.11e-3 |
| `slope`, δ = 0.47 (m = 1) | 81–178 | 2.0–2.8 | cosine 0.9992, size 1.10 | (b) fails: ratio 0.42, cosine 0.93, and the centre moves the opposite way; (b′): cosine 0.9986, ratio 1.12 |
| `slope`, δ = 0.47, all gains ÷ 4 | 81–178 | 2.0–2.8 | cosine 1.0000, size 1.05 | (b) fails: ratio 0.40, cosine 0.92; (b′): cosine 0.9999, ratio 1.065 |
| `slope`, δ = 0.047 | 654–1779 | 0.20–0.32 | cosine 1.0000, size 1.04 | (b): ratio 0.85, cosine 0.9986; (b′): ratio 1.05, cosine 1.0000 |
| `slope`, δ = 0.0047 | 13162–17790 | 0.020–0.023 | cosine 1.0000, size 1.04 | (b): ratio 1.02; (b′): ratio 1.04; both cosines 1.0000 |

In the first two rows, a / b means first quarter / last quarter (medians).

- **Joint descent.** In all runs, L̄ at the period-mean parameters is non-increasing, sampled once per period.
- **(H2b) is not enough.** The `slope` rows satisfy (H2b) by a wide margin, yet (b) fails until S is small. Meanwhile
  (b′) tracks the closed loop at every S. So (H2c) is the missing hypothesis of version 2's Proposition G(b), and (b′)
  is the correct leading-order slow flow. The configuration is the round-2 reviewer's counterexample, made
  reproducible in `exp_slow_manifold_slope.py`.
- **Error decomposition** (`slope`, δ = 0.47, started on the slow manifold):
  - ‖(a)‖/‖(b′)‖ = 1.015. This is the intrinsic O(1/R) error of (b′).
  - ‖closed‖/‖(a)‖ = 1.10. It falls to 1.05 when all gains are divided by 4, which leaves R and S unchanged and
    divides ε_avg ≈ 0.08 by 4. So this part is the averaging error.
  - The remaining ≈ 4.4% is the O(1/k_e) error-channel factor at k_e = 20. It is the same at δ = 0.047 and 0.0047.
- **The size ratio in the `satisfying` run.** The centred ratio of 0.97 is the error-channel factor times a +1.4%
  averaging correction from the weights' oscillation within each period, the analogue of Remark H.7a.
  - The frozen-parameter ratio is 0.951–0.956 along the run, and 0.996 at k_e = 80.
  - The fastest weight rate over ω_ref is about 0.14.
  - With the weights running at frozen geometry the ratio is 0.969–0.970; with g_a/4 it is 0.956–0.960 (review
    round 3).
  - The slope correction is below 0.01% here, because ∇L* is nearly orthogonal to the range of Da*ᵀ. Version 2's
    0.94–0.95 contained a 2.5% forward-difference bias.
- **The weight block.** In the `satisfying` run the weights sit at a* + O(1/k_e), within 0.1% of a*.
  - Their actual drift is tiny: 0.001–0.008 of the geometry speed in the G⁻¹ norm.
  - The prediction −(m/k_e)g_a∇_aL̄ at the period-mean weights is 0.17–0.30 of the geometry speed, 40–200 times the
    actual drift. It is the restoring velocity of the O(1/k_e) offset, which (a) does not contain.
  - So the weight block fails (cosine −0.21 / 0.50, carrying 92–95% of the residual), while the geometry block
    matches.

*Status.*

- **What the loop supplies.** The error-channel and averaging scales, and the data distribution (the occupation
  measure of the reference), come from the closed loop.
- **What the designer supplies.** The separation between weights and geometry is a gain choice, as in [S1], [S2]:
  δ ≪ g_a·min(m₀/‖Γ^{1/2}∇²L*Γ^{1/2}‖, 1/‖Da*Γ^{1/2}‖²). Equivalently, ‖G_ξ^{1/2}∇²L*G_ξ^{1/2}‖ ≪ g_a m₀ and
  ‖Da*G_ξ^{1/2}‖² ≪ g_a.
- **What is new here.**
  - The learner's data distribution is its own excitation, and at finite k_e it moves with the parameters
    (Proposition H.7(c)).
  - The slow flow is variable projection in the *induced* metric, valid for every S. The slope S measures the
    departure from plain VarPro.

---

## 6c. Theorem H — Langevin adaptive control and global convergence to the regularised optimum

**Why something new is needed.** Theorems B–G are local, or describe the dynamics. A global statement about kernel
placement needs either the Chizat–Bach theory or entropy.

- **Boundary conditions.** The Chizat–Bach theory for noiseless Wasserstein flows (partially 1-homogeneous case, [S6])
  requires boundary conditions: either Neumann conditions on a bounded convex domain, or C¹-uniform behaviour at
  infinity. Gaussian kernels with free centres and widths satisfy neither.
- **Separation.** Its Theorem 3.5 also needs the support of μ_0 to separate {−r₀}×Θ from {r₀}×Θ. A Dirac start
  violates this, and that is an obstruction distinct from the boundary conditions.
- **The collapsed configuration is invariant.** It remains a moving Dirac mass and never separates without noise, by
  permutation symmetry and uniqueness of solutions. `test_noiseless_symmetric_initialisation_never_separates_but_dither_does`
  checks this on the m = 6 real closed loop. The rows stay bit-identical at every snapshot, the common particle moves,
  and dither separates the rows within 0.1 s.

Entropy removes the obstruction: Theorem H(ii) covers Dirac starts after any positive time.

**Design.** The particles are unconstrained, ϑ = (α, γ, κ) ∈ ℝ^p, and the features are saturated coordinate by
coordinate:

    a = a_max tanh(α/a_max),    c = c_max tanh(γ/c_max),    s = s_min + (s_max − s_min)·sig(κ).

Here a_max bounds each coordinate of a, unlike §1 where it bounds ‖a‖. So Φ̃(x;ϑ) = aφ(x;c,s) satisfies
‖Φ̃‖ ≤ K_Φ := a_max√n, and every ϑ-derivative of Φ̃ is bounded uniformly in x. The law is

    dϑ_i = [ G∇_ϑg_t(ϑ_i) − σ(ϑ_i − ϑ⁰) ] dt + √(2τ) G^{1/2} dB_i,    g_t(ϑ) := e(t)ᵀΦ̃(x(t);ϑ).           (L_τ)

It has three parts: the classical gradient law, σ-modification written in the G⁻¹ metric, and parameter dither of
temperature τ.

**Averaged mean-field system.** Assume f = f(x) is continuous and time-invariant, as in Proposition G (the code's
`occupation_loss` and `first_variation` assume this). Let ρ be the occupation measure of a T₀-periodic reference, and
F_ρ := sup_{supp ρ}‖f‖. Define

    r_μ := f − ∫Φ̃(·;ϑ')dμ(ϑ'),    L̄(μ) := ½∫‖r_μ‖²dρ,    V_L[μ](ϑ) := −∫r_μ(x)ᵀΦ̃(x;ϑ)dρ(x),
    V[μ] := (1/k_e)V_L[μ] + (σ/2)‖ϑ − ϑ⁰‖²_{G⁻¹},

and the averaged system

    ∂_tμ = ∇·(μG∇V[μ]) + τ∇·(G∇μ).                                                                  (𝔄)

*Remark H.0 (from the closed loop to (𝔄): a formal reduction with a measured error).*

- **An algebraic identity.** Substitute e := r_μ(x_m)/k_e and x := x_m in the drift of (L_τ), and average over one
  reference period at frozen μ. The result is exactly −G∇V[μ].
  `test_averaged_drift_is_the_mean_field_langevin_drift` checks this identity to 1e-7. It does not test averaging.
- **The formal reduction.** (𝔄) is the leading-order description of the mean-field closed loop under two conditions:
  - *Tikhonov:* k_e ≫ max(ω_ref, Lip(f) + K_Φ/s_min, (F_ρ + K_Φ)/s_min), so that the shift x − x_m = O(1/k_e) is
    small on the kernel scale.
  - *Period averaging:* T₀·sup_ϑ‖(1/k_e)G∇V_L[μ](ϑ)‖ ≪ ℓ, σT₀ ≪ 1 and τλ_GT₀ ≪ ℓ², where ℓ is the ϑ-scale of Φ̃.
    The left-hand side is at most (T₀/k_e)λ_G(F_ρ + K_Φ)·sup‖∇_ϑΦ̃‖. The confinement drift does not depend on time and needs no
    averaging.
- **The first corrections.**
  - At frozen μ (and frozen ϑ) they are O(1/k_e) relative to the drift, do not depend on the gains, and are given
    exactly by Proposition H.7.
  - Along the loop, the within-period oscillation of the parameters adds the second-order averaging term b_osc of
    Remark H.7a. Its size relative to the frozen terms is O(|G|T₀), so it does depend on the gains.
  - In the formal limit below, both are O(k_e⁻²).

  Measured relative errors of the frozen-μ averaged drift against (𝔄), with b_osc excluded:

  | instance | k_e = 10 | 20 | 40 | 80 | k_e × error | test |
  |---|---|---|---|---|---|---|
  | n = 1 | 4.2% | 2.0% | 0.91% | 0.48% | ≈ 0.4 | `test_averaging_error_of_the_formal_reduction_decreases_like_one_over_ke` |
  | n = 2, ellipse | 43% | 21% | 10% | 5.1% | ≈ 4.2 | `test_averaged_drift_error_decays_like_one_over_ke_for_an_orbit_enclosing_area` (k_e = 10, 20, 40, 80) |

  - k_e = 10 is outside the averaging regime for the second instance.
  - Its larger error does not come from the enclosed area. H.7 attributes it to the orientation-even state-shift term:
    it predicts k_e × error ≈ 4.05, of which U_odd alone would give 1.5, and the reversed orbit gives 4.0 (review
    round 3).
  - At the gains of these instances b_osc is not small. k_e × its relative size is 0.92 (n = 1, G = diag(4,1,1)) and
    4.4 (n = 2, G = I), against 0.42 and 4.05 for the frozen terms (review round 3).
- **The exact formal limit** is k_e → ∞ with G, k_eσ and k_eτ fixed, in slow time t/k_e (Conjecture H.6).
- **Not proved here.** The mean-field limit of (L_τ) itself is not proved; Theorem D covers the deterministic
  projected law (L).

**Theorem H.** Assume the setting above, and σ > 0, τ > 0 where stated.

(i) *Identification.* V[μ] is the first variation of 𝒢(μ) := (1/k_e)L̄(μ) + (σ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ. Hence (𝔄) is the
mean-field Langevin dynamics, in the G⁻¹ metric, of the free energy F_τ := 𝒢 + τH, with H(μ) = ∫μ log μ.
Equivalently, with γ_τ := N(ϑ⁰, (τ/σ)G),

    F_τ = (1/k_e)L̄ + τH(·|γ_τ) + const.

F_τ is well defined on all of P(ℝ^p), and finite exactly when H(μ|γ_τ) < ∞.

(ii) *Global exponential convergence (σ, τ > 0).* F_τ has a unique minimiser, μ*_τ ∝ exp(−V[μ*_τ]/τ). For every
μ_0 ∈ P₂(ℝ^p) with F_τ(μ_0) < ∞,

    F_τ(μ_t) − F_τ(μ*_τ) ≤ e^{−2τρ_τ t}(F_τ(μ_0) − F_τ(μ*_τ)),
    W_{G⁻¹}(μ_t, μ*_τ)² ≤ (2/(τρ_τ))e^{−2τρ_τ t}(F_τ(μ_0) − F_τ(μ*_τ)).

Here ρ_τ is a log-Sobolev constant, for the whitened measures, of every proximal Gibbs measure μ̂ ∝ e^{−V[μ]/τ}. In
ϑ-coordinates this means Fisher information in the G metric. It satisfies

    ρ_τ ≥ (σ/τ) · exp( −2K_Φ(F_ρ + K_Φ)/(k_e τ) ).

Suppose μ_0 is only in P₂(ℝ^p), for example a Dirac mass or an empirical measure, so that H(μ_0) = +∞. Then
F_τ(μ_{t₀}) < ∞ for every t₀ > 0, and the estimates hold with (μ_0, t) replaced by (μ_{t₀}, t − t₀).

(iii) *Loss and tracking.* Let ē_t := r_{μ_t}/k_e. Then ⟨‖ē_t‖²⟩_ρ = 2L̄(μ_t)/k_e², and

    |L̄(μ_t) − L̄(μ*_τ)| ≤ (F_ρ + K_Φ)L_Φ·W_{G⁻¹}(μ_t, μ*_τ) ≤ (F_ρ + K_Φ)L_Φ(2/(τρ_τ))^{1/2}e^{−τρ_τ t}(F_τ(μ_0) − F_τ(μ*_τ))^{1/2},

where L_Φ := sup_{x,ϑ,‖u‖=1}‖G^{1/2}∇_ϑ(uᵀΦ̃)‖. The regularisation bias satisfies

    L̄(μ*_τ) ≤ inf_ν { L̄(ν) + k_eτH(ν|γ_τ) }.

The link e ≈ ē between (𝔄) and the closed loop is formal (Remark H.0).

(iv) *Annealing.* Let μ_0 ∈ P₂ with F_{τ_init}(μ_0) < ∞. Let (τ_t) be smooth and decreasing from τ_init, with
τ_t = α_a/log t for t ≥ t₁ > 1, where α_a > α* := 2K_Φ(F_ρ + K_Φ)/k_e. Then 𝒢(μ_t) − inf 𝒢 ≤ C'·log log t/log t, where inf 𝒢 is the
σ-regularised optimum. This is [S3] Thm 4.1. Its hypotheses hold:

- 𝒢 ≥ 0;
- Assumptions 1–3 hold for every τ > 0;
- for τ ∈ (0, τ_init], (ii) gives ρ_τ ≥ C₀e^{−α*/τ} with C₀ = σ/τ_init, which is all its proof uses.

*Proof.*

(i) d/dε L̄((1−ε)μ + εν)|₀ = −∫r_μᵀ∫Φ̃d(ν−μ)dρ = ∫V_L d(ν−μ). This is exact, because L̄ is quadratic in μ. The
confinement is linear in μ. In the G⁻¹ metric, ∇·(μG∇(V + τ log μ)) = ∇·(μG∇V) + τ∇·(G∇μ). Since
log γ_τ = −(σ/(2τ))‖ϑ−ϑ⁰‖²_{G⁻¹} + const, we get τH(μ|γ_τ) = τH(μ) + (σ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ + const.

(ii) Whiten with ζ := G^{−1/2}ϑ. The metric becomes Euclidean, the confinement becomes (σ/2)‖ζ − ζ⁰‖², W_{G⁻¹}
becomes W₂, and H shifts by ½log det G. (𝔄) becomes [S3] eq. (12). We verify Assumptions 1–3 of [S3].

- *Assumption 1.* ∇V[μ](ζ) is Lipschitz in (ζ, μ). The constants do not depend on μ, because
  ‖r_μ‖ ≤ F_ρ + K_Φ on supp ρ for every μ ∈ P(ℝ^p).
  - In ζ, this holds because ∇²Φ̃ is bounded (saturation, and s ≥ s_min).
  - In μ, it holds because ∫Φ̃dμ is W₁-Lipschitz, W₁ ≤ W₂, and ∇Φ̃ is bounded.
- *Assumption 2.* 𝒢 is convex: a linear map followed by a convex quadratic, plus a linear term. F_τ has a minimiser by
  confinement and lower semicontinuity ([S3] Prop. 5.1).
- *Assumption 3.* e^{−V[μ]/τ} = e^{−(σ/2τ)‖ζ−ζ⁰‖²}·e^{−V_L[μ]/(k_eτ)}.
  - The Gaussian factor satisfies LSI(σ/τ) (Bakry–Émery).
  - |V_L[μ]| ≤ ∫‖r_μ‖‖Φ̃‖dρ ≤ (F_ρ + K_Φ)K_Φ, so the perturbation has oscillation at most 2K_Φ(F_ρ + K_Φ)/(k_eτ).
    Here F_ρ may be replaced by ∫‖f‖dρ.
  - Holley–Stroock gives ρ_τ.

[S3] Theorem 3.2 and Corollary 3.3 (the latter with Otto–Villani) give the two inequalities. Uniqueness follows from
strict convexity of H.

*P₂ starts.* Work in ζ-coordinates and fix the solution μ_t. Compare the particle SDE with the Ornstein–Uhlenbeck
process dY = −σ(Y − ζ⁰)dt + √(2τ)dB started from the same law.

- The two drifts differ by b_t = −(1/k_e)G^{1/2}∇V_L[μ_t], with ‖b_t‖ ≤ B := (1/k_e)λ_G^{1/2}(F_ρ + K_Φ)·sup‖∇_ϑΦ̃‖.
- Girsanov and data processing give H(μ_{t₀}|OU_{t₀}) ≤ B²t₀/(4τ).
- OU_{t₀} has a density at most (2π(τ/σ)(1 − e^{−2σt₀}))^{−p/2}, so H(μ_{t₀}) < ∞.
- μ_{t₀} ∈ P₂, because the drift grows linearly.

Hence F_τ(μ_{t₀}) < ∞.

(iii) L̄(μ) − L̄(ν) = ½∫(r_μ − r_ν)ᵀ(r_μ + r_ν)dρ. Here ‖r_μ + r_ν‖ ≤ 2(F_ρ + K_Φ), and
‖r_μ − r_ν‖ = ‖∫Φ̃d(ν − μ)‖ ≤ L_Φ·W_{1,G⁻¹} ≤ L_Φ·W_{G⁻¹}. Conclude with (ii). For the bias:
(1/k_e)L̄(μ*_τ) + τH(μ*_τ|γ_τ) ≤ (1/k_e)L̄(ν) + τH(ν|γ_τ) for every ν, and H(μ*_τ|γ_τ) ≥ 0.

(iv) This is [S3] Thm 4.1 with the stated hypotheses. ∎

*Remark H.1 (dither adds no output noise in the mean-field limit; its cost is the entropic bias).*

- In the limit, f̂ = ∫Φ̃dμ is deterministic by the law of large numbers.
- The only effect of dither is the entropic bias, which Theorem H(iii) quantifies.
- In the experiment (σ = 0.02, τ = 0.004, G = diag(4,1,1), k_e = 10), a grid computation gives L̄(μ*_τ) = 0.11118.
  It is grid-converged (`exp_gibbs_reference.py` → `results/langevin_global/gibbs_reference.md`).
- At (σ, τ) = (0.002, 5·10⁻⁴) the same computation gives 0.0085.
- The unregularised offline optimum over 32 kernels is about 7.6·10⁻⁶ (best of 5 L-BFGS restarts, an upper bound), and
  L̄(0) = 0.311.
- At finite m, the fluctuation of f̂ is O(K_Φ/√m) at each fixed t (finite-horizon propagation of chaos). It is
  uniform in time only under Remark H.5.

*Remark H.2 (the control knobs, and what the numbers say).*

- **The knobs.** σ is the strong-convexity constant of the confinement, τ is the temperature, and a_max sets K_Φ.
- **The role of k_e.** Since k_eF_τ = L̄ + (k_eσ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ + k_eτH, μ*_τ and ρ_τ depend on (σ, τ, k_e)
  only through (k_eσ, k_eτ). Raising k_e at fixed (σ, τ) has three effects:
  - it improves ρ_τ and the guaranteed rate 2σe^{−α*/τ}, since α* ∝ 1/k_e;
  - it moves μ*_τ toward the prior N(ϑ⁰, (τ/σ)G), which increases the regularisation bias;
  - it shrinks the loss-driven drift like 1/k_e.
- **The guaranteed rate** is exponentially small in 1/τ, a known feature of Holley–Stroock bounds:
  - 2τρ_τ ≥ e^{−8401} s⁻¹ at (σ, τ) = (0.002, 5·10⁻⁴);
  - 2τρ_τ ≥ e^{−1053} s⁻¹ at (0.02, 0.004).

  The experiments do not test this rate.
- **The temperature sweep** (`results/langevin_tradeoff`; m = 32, σ = 0.002, k_e = 10, T = 1500 s, 2 seeds):

  | τ | escape time of the collapsed start (s), per seed | loss floor of the spread start |
  |---|---|---|
  | 0 | never, never | 7.0e-3 |
  | 5·10⁻⁴ | never, 1050 | 1.1e-2 |
  | 1·10⁻³ | 1480, 1005 | 1.5e-2 |
  | 2·10⁻³ | 1115, 800 | 2.3e-2 |
  | 4·10⁻³ | 970, 430 | 3.6e-2 |
  | 8·10⁻³ | 390, 295 | 6.5e-2 |
  | annealed, hot start | 990, 530 | 3.2e-2 |

- **How escape is measured.** Escape is the first sample at which the collapsed run's loss falls to at most twice the
  spread run's final loss, at the same τ and seed. Two caveats follow:
  - the threshold grows with τ, from 0.022 to 0.13, and the escape time includes convergence after escape;
  - the τ = 5·10⁻⁴ point is censored and is excluded from the fit.
- **The fit.** The seed-averaged escape time fits t_esc ≈ −1693 + 427·log(1/τ) (4.5% relative residual) better than an
  Arrhenius law (25%). On the eight per-seed values the two laws fit equally poorly (27% and 33%;
  `langevin_tradeoff/table.md`).
- **No saddle explanation.** The noiseless collapsed end point is a saddle of the averaged landscape, but a weak one
  (review round 3). Its unstable rate is about 5·10⁻⁵ s⁻¹ for (𝔄), and about 2·10⁻⁴ s⁻¹ for the exactly averaged loop
  at k_e = 10. A saddle-escape law would therefore need τ ≲ 10⁻⁵, and would have a log(1/τ) slope of 10³–10⁴ s.
- **What the data show.** At the tested temperatures escape is governed by diffusion, by the approach transient and by
  the τ-dependent threshold. The data do not identify the escape mechanism. Version 2's "saddle, not barrier" reading
  is withdrawn.
- **The annealed runs** use α_a between about 2·10⁻³ and 2·10⁻², far below α* = 4.2. They are not instances of H(iv).

*Remark H.3 (why saturation and not projection).* The mean-field Langevin results are on ℝ^d. Projecting the drift
with the Pomet–Praly operator, while keeping additive noise, leaves the noise unconfined. Reflected (Skorokhod)
dynamics on the convex body Θ_K would confine it.

- Such dynamics satisfy a log-Sobolev inequality with constant of order 1/diam², even for σ = 0.
- The argument is a diameter-based log-Sobolev inequality for the uniform measure on a convex body, followed by
  Holley–Stroock. The diameter-based inequality comes from CD(0,∞) with a convex Neumann boundary plus a diameter bound
  (F.-Y. Wang 1997), or from Frieze–Kannan (1999) for log-concave measures on convex bodies. Together these give
  LSI(c·diam⁻²·e^{−2K_Φ(F_ρ+K_Φ)/(k_eτ)}).
- Bakry–Émery alone gives nothing at σ = 0.
- That route connects to [25] and is left open.

*Remark H.4 (the role of σ-modification).* In the classical theory, σ-modification is a robustness device with a bias.
In the saturated design on ℝ^p it supplies the confinement that makes e^{−V/τ} integrable and log-Sobolev, so the
global result rests on it. Its bias is the regularisation term of the variational problem that μ*_τ solves. Under
reflection on Θ_K, compactness would play that role instead.

*Remark H.5 (finite m; conditional).* The m-particle averaged system is

    dϑ_i = −G∇V[μ^m_t](ϑ_i)dt + √(2τ)G^{1/2}dB_i,    μ^m_t := (1/m)Σδ_{ϑ_i}.

Uniform-in-time propagation of chaos would follow from [S5], under its hypotheses.

- *Available:* convexity of 𝒢; a uniform LSI, from (ii); and a bounded, smooth second variation
  (1/k_e)∫Φ̃(x;ϑ)ᵀΦ̃(x;ϑ')dρ.
- *Still to check:* the remaining regularity conditions of [S5], verbatim (`LITERATURE_TRANSPORT.md` open question 4).

**Conjecture H.6 (the real closed loop, in the scaling where (𝔄) is the limit).** Let f = f(x) satisfy the tube
condition of Theorem D(b), with a_max replaced by K_Φ. Fix G ≻ 0 and σ₀, τ₀ > 0, and set σ := σ₀/k_e and τ := τ₀/k_e.

- Then μ*_τ, the minimiser of L̄ + (σ₀/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ + τ₀H, does not depend on k_e.
- In slow time s = t/k_e, (𝔄) is the formal limit of the closed loop.

*Conjecture:* for the mean-field limit of the closed loop under (L_τ), for every μ_0 ∈ P₂ and every e(0), there are
S(μ_0) and r(k_e) → 0 such that W_{G⁻¹}(μ_t, μ*_τ) ≤ r(k_e) for all t ≥ k_e·S.

At fixed k_e the analogous statement is false. The limit of small gains is the exactly averaged system, with the true
periodic error e_μ and x = x_m + e_μ.

- Its stationary state is offset from μ*_τ: L̄ = 0.11415 against L̄(μ*_τ) = 0.11118 at k_e = 10
  (`results/langevin_global/gibbs_reference.md`).
- In the scaling of H.6, the L²(ρ) distance between the two learned functions is 6.8·10⁻³, 1.4·10⁻³ and 9.7·10⁻⁴ at
  k_e = 10, 20 and 40. It is small and decreasing, but a clean 1/k_e law is not resolved in this range.
- By Proposition H.7(e) it is not the W_{G⁻¹}-gradient flow of any free energy. The cyclic sums of its antisymmetric
  second variation are non-zero in the checked instances, including the network of the experiments (review round 3);
  by Theorem I(iii) they are non-zero at every n = 1 measure whose residual is non-zero on the orbit, and by Theorem
  I(vi) on no open set of measures do they all vanish.

**Proposition H.7 (frozen-measure structure of the exactly averaged closed loop).** Assume the following:

- f ∈ C²(ℝⁿ) is autonomous;
- x_m ∈ C² is T₀-periodic, with oriented orbit 𝒪;
- the features are the saturated features of the design;
- μ ∈ P(ℝ^p) is frozen.

Then for every R > 0 there is k₀(R), independent of μ, such that for k_e ≥ k₀(R) the following hold.

(a) *Periodic error.* The frozen-μ error equation ė = r_μ(x_m + e) − k_e e has exactly one T₀-periodic solution e_μ
with sup_t‖e_μ‖ ≤ R. It attracts exponentially every solution with ‖e(0)‖ ≤ R, and

    e_μ = r_μ(x_m)/k_e + (∇r_μ(x_m) r_μ(x_m) − (d/dt) r_μ(x_m(t)))/k_e² + O(k_e⁻³),

uniformly in t and μ.

- Here ∇r_μ is the Jacobian (∂r_k/∂x_l), and ∂_xΦ̃ = (∂Φ̃_k/∂x_l).
- If f is one-sided Lipschitz, ⟨f(x) − f(y), x − y⟩ ≤ L_f‖x − y‖², then for k_e > L_f + a_max√n·e^{−1/2}/s_min the
  periodic solution is unique and attracts every solution.
- Without such a condition uniqueness is only local. For f(x) = x² there is a second, repelling periodic solution
  near e ≈ k_e, and solutions above it blow up (review round 3).

(b) *Expansion.* The period average of the drift of (L_τ) along e_μ is G∇_ϑP_μ(ϑ) − σ(ϑ − ϑ⁰), where
P_μ(ϑ) := ⟨e_μᵀΦ̃(x_m + e_μ;ϑ)⟩_t satisfies

    P_μ(ϑ) = −(1/k_e)V_L[μ](ϑ) + (1/k_e²)(U_even[μ](ϑ) + U_odd[μ](ϑ)) + O(k_e⁻³),
    U_even[μ](ϑ) := ∫(r_μ·∇_x)(r_μᵀΦ̃(·;ϑ))dρ,    U_odd[μ](ϑ) := −(1/T₀)∮_𝒪 Φ̃(·;ϑ)ᵀdr_μ,

uniformly in ϑ and μ, together with ϑ-gradients.

(c) *Semi-gradient form.* −(1/k_e)V_L[μ] + (1/k_e²)U_even[μ] = (1/k_e)∫r_μᵀΦ̃(·;ϑ)dρ^cl_μ + O(k_e⁻³). Here
ρ^cl_μ is the occupation measure of the closed-loop state x_m + e_μ.

- So, apart from the lag term U_odd, the averaged law follows the gradient of the closed-loop loss ½∫‖r_ν‖²dρ^cl_μ at
  ν = μ, with the closed-loop state distribution ρ^cl_μ held fixed.
- This is a semi-gradient: the update freezes part of its own parameter dependence. In TD learning the frozen part is
  the bootstrapped target, and the data distribution does not depend on the parameters [S10], [S11].
- Here the frozen part is the data distribution ρ^cl_μ. This is the structure of performative prediction [S14]: steps
  taken with the distribution held at its current value converge to performatively stable points, not to optima.
- **The performative defect carries the even non-reciprocity.** With J^cl(μ) := ½∫‖r_μ‖²dρ^cl_μ,

      −(1/k_e)V_L[μ] + (1/k_e²)U_even[μ] = −(1/k_e)δJ^cl/δμ − (1/k_e²)∫(Φ̃_ϑ·∇)(½‖r_μ‖²)dρ + O(k_e⁻³).

  The first term is a true first variation, hence reciprocal. The second is the defect: it has the same order as
  U_even, and it carries all of A_even in (e).

(d) *Circulation kernel.* U_odd[μ](ϑ) = −(1/T₀)∮_𝒪Φ̃(·;ϑ)ᵀdf + ∫K(ϑ,ϑ')dμ(ϑ'), with

    K(ϑ,ϑ') := (1/T₀)∮_𝒪 Φ̃(·;ϑ)ᵀ dΦ̃(·;ϑ') = ⟨Φ̃(x_m;ϑ)ᵀ (d/dt)Φ̃(x_m;ϑ')⟩_t .

- K is antisymmetric.
- By Stokes, K(ϑ,ϑ') = (1/T₀)∫_Σ Σ_k dΦ̃_k(·;ϑ) ∧ dΦ̃_k(·;ϑ') for any 2-chain Σ with ∂Σ = 𝒪.
- For n = 2 and a simple, positively oriented orbit bounding D,
  K(ϑ,ϑ') = (1/T₀)Σ_k∫_D det[∇_xΦ̃_k(·;ϑ), ∇_xΦ̃_k(·;ϑ')] dA.
- K is an integral over the enclosed region, not the area itself. For Gaussian features,
  Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ')] = (aᵀa')·det[∇φ, ∇φ'], so K = 0 whenever aᵀa' = 0 or c = c', whatever the area.
- Hence U_odd ≡ 0 for n = 1.
- Reversing the orbit preserves ρ, flips U_odd, and leaves U_even unchanged.

(e) *Non-reciprocity.* Write S(ϑ,ϑ') := δP_μ(ϑ)/δμ(ϑ') for the derivative in the direction δ_ϑ', and call
S(ϑ,ϑ') − S(ϑ',ϑ) its antisymmetric part (without the factor ½). Adding a function of ϑ alone, or a constant c(μ), to
P_μ changes the antisymmetric part by a coboundary, which the criterion below allows.

- At order 1/k_e, S is symmetric: it equals −(1/k_e)∫Φ̃(·;ϑ)ᵀΦ̃(·;ϑ')dρ.
- At order 1/k_e² it is not. The second variation of the k_e⁻² coefficient U_even + U_odd has antisymmetric part
  A_even + 2K, with

    A_even(ϑ,ϑ') := −∫[ (Φ̃_ϑ'·∇)(r_μᵀΦ̃_ϑ) − (Φ̃_ϑ·∇)(r_μᵀΦ̃_ϑ') ] dρ .

- P_μ, up to O(k_e⁻³), can be a first variation of a functional of μ only if A_even + 2K is a coboundary
  a(ϑ) − a(ϑ'), that is, only if all its cyclic sums vanish.
- They do not vanish in the checked instances, for n = 1 (from A_even alone) and for n = 2. Theorem I shows that the
  cyclic condition is equivalent to A_even + 2K ≡ 0, and characterises when that happens: for n = 1 exactly when the
  residual vanishes on the orbit; for n ≥ 2 at perfect fit exactly when the orbit carries no circulation and the
  residual's Jacobian is symmetric along the orbit. It also shows (Theorem I(vi)) that no open set of measures satisfies
  the condition, for any n.
- So at finite k_e, the exactly averaged loop is a non-reciprocal McKean–Vlasov dynamics: it is not the
  W_{G⁻¹}-gradient flow of any free energy on any open set of measures. A different, μ-dependent metric is not excluded
  by this argument.

*Proof.*

(a) Write e = r_μ(x_m)/k_e + w. Then w solves ẇ = −k_e w + g(t, w) with g = ∇r_μ(x_m)e − ṙ_μ/k_e + O(|e|²).

- The tube {‖e‖ ≤ R} is forward invariant once k_eR > sup_{‖e‖≤R}‖r_μ(x_m + e)‖.
- Inside the tube, solutions contract at rate k_e − Lip_R(r_μ), which gives existence, uniqueness and attraction there.
  A contraction on T₀-periodic functions gives the expansion.
- The bounds on r_μ, ∇r_μ and ∇²r_μ on the R-neighbourhood of the orbit are uniform in μ, by the saturation and
  f ∈ C².
- The global statement uses L_f in place of the local Lipschitz constant of f.

(b) Taylor-expand eᵀΦ̃(x_m + e;ϑ) = eᵀΦ̃(x_m;ϑ) + eᵀ∂_xΦ̃(x_m;ϑ)e + O(|e|³), and insert (a).

- The period averages of the time-symmetric terms are integrals against ρ.
- (∇r r)ᵀΦ̃ + rᵀ(∂_xΦ̃)r = (r·∇)(rᵀΦ̃).
- ϑ-derivatives commute with the average.
- The remainders are uniform, because the x- and ϑ-derivatives of Φ̃ of the orders needed are bounded.

(c) ∫r_μᵀΦ̃_ϑ dρ^cl_μ = ⟨r_μ(x_m + e_μ)ᵀΦ̃(x_m + e_μ;ϑ)⟩ = ∫r_μᵀΦ̃_ϑ dρ + ⟨∇_x(r_μᵀΦ̃_ϑ)·e_μ⟩ + O(k_e⁻²). By (a)
this equals ∫r_μᵀΦ̃_ϑ dρ + (1/k_e)U_even + O(k_e⁻²).

(d) (d/dt)r_μ(x_m) = (d/dt)f(x_m) − ∫(d/dt)Φ̃(x_m;ϑ')dμ(ϑ').

- Antisymmetry: K(ϑ,ϑ') + K(ϑ',ϑ) = (1/T₀)∮d(Φ̃_ϑᵀΦ̃_ϑ') = 0.
- The area form is Stokes, or Green's theorem for n = 2.
- For n = 1, every 1-form on ℝ is exact, so its loop integrals vanish.
- Reversal changes the sign of oriented loop integrals and preserves ρ.

(e) The second variation of U_even is S_even(ϑ,ϑ') = −∫[(Φ̃_ϑ'·∇)(r_μᵀΦ̃_ϑ) + (r_μ·∇)(Φ̃_ϑ'ᵀΦ̃_ϑ)]dρ, whose second term
is symmetric. The second variation of U_odd is K. A kernel's antisymmetric part A is a coboundary if and only if
A(ϑ₀,ϑ₁) + A(ϑ₁,ϑ₂) + A(ϑ₂,ϑ₀) = 0 for all triples.

The remainder of (b) is not differentiated in μ. Instead the argument uses a circulation around a small triangle of
measures.

1. Fix ϑ₀, ϑ₁, ϑ₂ and ε ∈ (0,1). For (a,b) in T_ε := {a, b ≥ 0, a + b ≤ ε}, let
   μ_{a,b} := (1−ε)μ + (ε−a−b)δ_{ϑ₀} + aδ_{ϑ₁} + bδ_{ϑ₂}.
2. The circulation ∮_{∂T_ε}⟨Q_{μ_{a,b}}, dμ_{a,b}⟩ vanishes when Q_μ is a first variation plus a constant c(μ),
   because dμ_{a,b} has zero mass. In particular it vanishes for −(1/k_e)V_L.
3. U_even + U_odd is quadratic in μ, so its circulation equals −(ε²/2)·Σ_cyc(A_even + 2K)(ϑ_i, ϑ_{i+1}) at the
   centroid of T_ε.
4. The O(k_e⁻³) remainder of (b) contributes at most 6ε·sup|R| = O(εk_e⁻³).
5. Suppose that, for all large k_e, P_μ were a first variation up to O(k_e⁻³), with constants uniform in k_e. Then
   (ε²/(2k_e²))|cyc| ≤ Cεk_e⁻³, that is, |cyc| ≤ 2C/(εk_e). The cyclic sum does not depend on k_e, so letting
   k_e → ∞ at fixed ε gives cyc = 0 at the centroid measure of the triangle; letting ε → 0, since A_even + 2K is affine
   in the measure, gives cyc = 0 at μ.

This is checked: for ε = 0.05 and n = 2 the circulation equals the centroid formula to all printed digits (review
round 3). ∎

*Checked* (`tests/test_averaging.py`, using `src/core/averaging.py`; seven tests: four for H.7, one for Remark H.0, two
for Remark H.7a):

| check | result |
|---|---|
| expansion: relative error of k_e²(P_μ − P₀) against U_even + U_odd | halves when k_e doubles. At k_e = 40 → 80 it goes 0.54 → 0.25 for n = 1, 0.038 → 0.019 for n = 2 counter-clockwise, and 0.022 → 0.011 clockwise. A Richardson step at k_e = 80 and 160 brings n = 1 within 3% (asserted) |
| Stokes: loop-integral K against the area integral over the enclosed ellipse | agreement to 7e-15 relative (asserted at 1e-6); K antisymmetric; K ≡ 0 for n = 1 |
| parity under orbit reversal (same ρ) | U_even unchanged and U_odd flipped, to 2e-15 relative (asserted at 1e-10) |
| non-reciprocity: the *exact* second variation S (k_e = 80) against the prediction | the full matrix k_e²S + k_eS₀ is within 5% of S_even + K. The cyclic sums of the antisymmetric part agree to 3%: n = 1 gives −0.0711 against −0.0704, all from the state shift; n = 2 gives −0.0471 against −0.0464 (state shift −0.089, circulation +0.043) |

*Remark H.7a (the full first-order correction, and a test in the real closed loop).* Proposition H.7 freezes the
measure over a period. In the closed loop the parameters also oscillate within each period, by O(1/k_e).

*The parameter equation.* In the scaling σ = σ₀/k_e, write

    Θ̇ = (1/k_e)H(Θ,t) + (1/k_e²)H₂(Θ,t) + O(k_e⁻³),
    H_i(Θ,t) := G∇_ϑ[r_Θ(x_m(t))ᵀΦ̃(x_m(t);ϑ_i)] − σ₀(ϑ_i − ϑ⁰).

*Second-order averaging* (Sanders–Verhulst–Murdock §2.9) gives the dynamics of the period mean Y:

    Ẏ = (1/k_e)⟨H⟩(Y) + (1/k_e²)[ ⟨H₂⟩(Y) + ⟨D_ΘH(Y,t)·u(Y,t)⟩ ] + O(k_e⁻³).

- u is the zero-mean antiderivative in t of H − ⟨H⟩.
- ⟨H₂⟩ is the frozen-measure correction of Proposition H.7.
- This is formal. The averaging theorem gives it on finite slow-time horizons.
- `averaging.oscillation_correction` computes the oscillation term b_osc := ⟨D_ΘH·u⟩.

Three properties follow.

- **Parity and value.** b_osc is odd under reversal of the orbit
  (`test_oscillation_correction_is_odd_under_orbit_reversal`). Its value matches an independent one-period
  near-identity flow at k_e = 80: the test asserts 3%, and 0.7% is measured
  (`test_oscillation_correction_matches_a_one_period_flow`).
- **Size.** Relative to the frozen-measure terms it scales like |G|T₀, the parameter displacement per period times
  k_e. It is negligible when adaptation is slow within a period, which is the small-gain limit in the fixed-k_e
  discussion of Conjecture H.6. It is not negligible as k_e → ∞ at fixed G. Both terms are O(k_e⁻²) and uniformly
  bounded, so the order of the route to Conjecture H.6 is unchanged.
- **Non-reciprocity survives.** The only orientation-even part of the full correction is U_even, and its antisymmetric
  second variation is not a coboundary (Proposition H.7(e)). Suppose the full first-order dynamics were reciprocal for
  both orientations. Then the average of the two antisymmetric parts, A_even, would be a coboundary, which it is not.
  So the full first-order averaged dynamics is non-reciprocal for at least one orientation.

*A prediction tested in the real closed loop* (`exp_orbit_orientation.py` → `results/orbit_orientation/`).

- **The prediction.** Traverse the same elliptic reference both ways. The occupation measure is the same, so (𝔄),
  and any functional of ρ, cannot distinguish the two runs. The first-order theory predicts that the learned kernels
  differ by Θ₊ − Θ₋ ≈ −2J⁻¹(b_odd + b_osc) = O(1/k_e). Here J is the Jacobian of the leading averaged drift at the
  equilibrium.
- **The setting.** n = 2, m = 6, G = diag(4,4,1,1,1), σ = 2/k_e, no dither, and the same initial parameters for both
  runs. The runs converged to machine precision. The six kernels start spread out (spread 1.05) and converge to a
  common configuration (spread at most 5.6e-17 at the end of every run), so the learned network is effectively a
  single kernel: the test probes the orientation effect at a collapsed equilibrium of the m = 6 system.

| k_e | measured ‖Θ₊ − Θ₋‖ | frozen measure only (H.7): size ratio, cosine | full (H.7 + oscillation): size ratio, cosine | full: vector error |
|---|---|---|---|---|
| 10 | 5.7e-3 | 2.67, 0.76 | 0.94, 0.93 | 0.36 |
| 20 | 2.8e-3 | 2.41, 0.83 | 0.97, 0.98 | 0.20 |
| 40 | 1.4e-3 | 2.29, 0.86 | 0.97, 0.994 | 0.11 |

The size ratio is measured over predicted.

- The measured difference scales like k_e^{−1.01}.
- At k_e = 40 the full first-order theory predicts it to within 3% in size and 11% as a vector (cosine 0.994).
- The vector error halves with each doubling of k_e: 36%, 20% and 11% at k_e = 10, 20 and 40.
- The frozen-measure terms alone account for less than half of it.
- A first attempt with a weaker leak (σ = 0.5/k_e) and a shorter horizon had not converged: its residual drift was a
  third of the difference. It is not used.

The direction in which the reference orbit is traversed is therefore a property of the excitation that the learned
kernels retain, at order 1/k_e. The first-order theory predicts its effect.

**Theorem I (reciprocity at order k_e⁻²: a pairwise criterion and its characterisation).** Assume the setting of
Proposition H.7: saturated Gaussian features Φ̃(x;ϑ) = aφ_ξ(x) with ξ := (c, s), f ∈ C²(ℝⁿ) autonomous, x_m ∈ C²
T₀-periodic with oriented orbit 𝒪 and occupation measure ρ, and a frozen μ ∈ P(ℝ^p). Let S be the k_e⁻² coefficient of
the second variation of P_μ, and A_μ := S − Sᵀ = A_even + 2K its antisymmetric part, as in H.7(e). Call μ *reciprocal
at order k_e⁻²* if A_μ ≡ 0 (the 1-form that the k_e⁻² drift defines is closed at μ). Note that A_μ is affine in μ, because
r_μ is and K does not depend on μ.

(i) *Pairwise criterion.* A_μ is a coboundary a(ϑ) − a(ϑ') if and only if A_μ ≡ 0. Consequently, by Proposition H.7(e),
whose hypothesis concerns a neighbourhood of μ (its triangle of measures), P_ν can be a first variation of a functional
of ν up to O(k_e⁻³), with constants uniform in k_e and in ν near μ, only if A_ν ≡ 0 for every ν in a neighbourhood of μ.

(ii) *Structure.* With Dr the Jacobian of r_μ, (Dr)_{kl} = ∂_l r_k, and [X, Y] := (X·∇)Y − (Y·∇)X the Lie bracket of
vector fields on ℝⁿ,

    A_even(ϑ,ϑ') = −∫ { Φ̃_ϑᵀ(Dr_μ − Dr_μᵀ)Φ̃_ϑ' + r_μᵀ[Φ̃_ϑ', Φ̃_ϑ] } dρ.

For the Gaussian features A_μ is a bilinear form in the weights: A_μ(ϑ,ϑ') = aᵀB_μ(ξ,ξ')a' with

    B_μ(ξ,ξ') := −∫(Dr_μ − Dr_μᵀ)φ_ξφ_ξ' dρ − ∫ r_μ(∇φ_ξ)ᵀφ_ξ' dρ + ∫(∇φ_ξ')r_μᵀφ_ξ dρ + 2κ(ξ,ξ')I_n,
    κ(ξ,ξ') := ⟨φ_ξ(x_m(t))·(d/dt)φ_ξ'(x_m(t))⟩_t,

so A_μ ≡ 0 if and only if B_μ(ξ,ξ') = 0 for all admissible ξ, ξ'. The first term is the curl of the residual paired
with the features, the next two are the residual paired with the Lie bracket of the feature fields, and the last is the
circulation kernel K(ϑ,ϑ') = (aᵀa')κ(ξ,ξ'). The performative part is bounded by the residual and by the antisymmetric
part of its Jacobian along the orbit:

    |A_even(ϑ,ϑ')| ≤ |a||a'| [ ∫‖Dr_μ − Dr_μᵀ‖dρ + (2e^{−1/2}/s_min)∫|r_μ|dρ ].

(iii) *n = 1.* A_μ ≡ 0 if and only if r_μ = 0 ρ-a.e. on the orbit. So the one-dimensional loop is reciprocal at
order k_e⁻², at a frozen μ, exactly when μ interpolates f along the orbit, and its non-reciprocity is bounded by the
residual: |A_μ(ϑ,ϑ')| ≤ (2e^{−1/2}/s_min)|a||a'|∫|r_μ|dρ.

(iv) *n ≥ 2 at perfect fit.* Suppose r_μ = 0 ρ-a.e. Then the Lie-bracket terms vanish and

    B_μ(ξ,ξ') = −∫(Dr_μ − Dr_μᵀ)φ_ξφ_ξ' dρ + 2κ(ξ,ξ')I_n,

an antisymmetric matrix plus a multiple of the identity. Hence A_μ ≡ 0 if and only if both

    (a) K ≡ 0, which holds if and only if the 1-current of the oriented orbit vanishes,
        T_𝒪 := ∫₀^{T₀} δ_{x_m(t)} ẋ_m(t) dt = 0  (an ℝⁿ-valued measure; ∫v·dT_𝒪 = ∮_𝒪 v·dx); for n = 2, iff the
        winding-number function w_𝒪 of the orbit vanishes a.e.; and
    (b) Dr_μ = Dr_μᵀ ρ-a.e.: the Jacobian of the residual has no antisymmetric part along the orbit.

Vanishing of r_μ on the orbit constrains only its tangential derivative, so (b) is a genuine extra condition: for
f := f̂_μ + h e₁ with h = 0 on the orbit and ∂_⊥h ≠ 0, A_μ − 2K = −∫φ_ξφ_ξ'∂_⊥h dρ·(e₁e_⊥ᵀ − e_⊥e₁ᵀ) ≠ 0. When r_μ
vanishes on a neighbourhood of the orbit, Dr_μ = 0 there and A_μ = 2K: the semi-gradient defect has vanished and the
loop is non-reciprocal exactly when the orbit carries circulation, through the orbit-area kernel alone.

(v) *General n ≥ 2.* In any connected finite-dimensional real-analytic family of measures, for example atomic measures
with J atoms parametrised by their weights and lifted positions, the set of reciprocal measures is either the whole
family or a closed real-analytic subvariety of measure zero. The instances of `tests/test_averaging.py` (n = 1 and the
n = 2 ellipse) have A_μ ≢ 0, with cyclic sums −0.0704 and −0.0464 of the predicted kernel.

(vi) *Nowhere a gradient flow.* For every n, every (f, x_m) and every μ, every neighbourhood of μ contains measures ν
with A_ν ≢ 0: along the segments ν_δ := (1 − δ)μ + δδ_{ϑ'} either the zero feature ϑ' = ϑ_z or a feature with a' ≠ 0
gives A_{ν_δ} ≢ 0 for all small δ > 0. Hence, by (i), there is no open set of measures on which P_ν is a first
variation up to O(k_e⁻³) with uniform constants: the exactly averaged closed loop is nowhere the W_{G⁻¹}-gradient flow
of a free energy at order k_e⁻². The reciprocal measures of (iii) and (iv) are coincidences without interior.

*Proof.* (i) Let ϑ_z := (a = 0, ξ) be a zero-feature particle (the unconstrained weight coordinate equal to 0, any geometry),
an admissible point of ℝ^p, so Φ̃(·;ϑ_z) ≡ 0 as a function of x; its x-derivatives vanish too. Every term of A_even(ϑ,ϑ_z) contains Φ̃_{ϑ_z} or (Φ̃_{ϑ_z}·∇), and
K(ϑ,ϑ_z) = ⟨Φ̃_ϑᵀ(d/dt)Φ̃_{ϑ_z}⟩ = 0; so A_μ(ϑ,ϑ_z) = 0 for every ϑ. If A_μ(ϑ,ϑ') = a(ϑ) − a(ϑ'), then a(ϑ) = a(ϑ_z) for
every ϑ, so A_μ ≡ 0. The converse is trivial. For the consequence, apply H.7(e) at any ν near μ with ϑ₂ := ϑ_z in its
triangle of measures: the cyclic sum over (ϑ₀, ϑ₁, ϑ_z) is A(ϑ₀,ϑ₁), evaluated at the centroid measure of the triangle,
which tends to ν as the triangle's size ε → 0; A is affine in the measure, so A_ν(ϑ₀,ϑ₁) = 0.

(ii) Expand (Φ̃'·∇)(rᵀΦ̃) = Σ_{k,l}Φ̃'_l(∂_l r_k)Φ̃_k + Σ_{k,l}Φ̃'_l r_k∂_lΦ̃_k = Φ̃ᵀ(Dr)Φ̃' + rᵀ(Φ̃'·∇)Φ̃, and subtract the
expression with ϑ and ϑ' exchanged. For Φ̃ = aφ, (Φ̃'·∇)Φ̃ = a(a'ᵀ∇φ)φ', so rᵀ[Φ̃',Φ̃] = (rᵀa)(a'ᵀ∇φ)φ' − (rᵀa')(aᵀ∇φ')φ,
and K(ϑ,ϑ') = (aᵀa')κ(ξ,ξ'). Collect the bilinear form in (a, a'). Since a ranges over the open cube (−a_max, a_max)ⁿ,
A_μ ≡ 0 iff B_μ ≡ 0. The bound uses φ ≤ 1 and |∇φ| ≤ e^{−1/2}/s_min.

(iii) For n = 1, Dr − Drᵀ = 0, and for equal widths s, φ_c'Dφ_c − φ_cDφ_c' = φ_cφ_c'(c − c')/s², so
B_μ(ξ,ξ') = −((c − c')/s²)∫r_μφ_cφ_c'dρ. Now φ_cφ_c' = e^{−(c−c')²/(4s²)}·e^{−(x−m)²/s²} with m := (c + c')/2. Put
ν := r_μρ, a finite signed measure with compact support. The function m ↦ ∫e^{−(x−m)²/s²}dν(x) extends to an entire
function. If B_μ ≡ 0 on admissible pairs, it vanishes for all m ∈ (−c_max, c_max) (take c = m − δ, c' = m + δ with
0 < δ small), hence everywhere; its Fourier transform is the product of the Gaussian's, which has no zeros, with ν̂, so
ν = 0, that is r_μ = 0 ρ-a.e. Conversely, if ν = 0 then A_even(ϑ,ϑ') = −∫(φ'Dφ − φDφ')aa' dν = 0 for all widths.

(iv) If r_μ = 0 ρ-a.e. the two Lie-bracket integrals vanish, which gives the displayed B_μ; the first term is an
antisymmetric matrix and the second a multiple of I_n, so B_μ ≡ 0 iff both vanish (take the trace, which is 2nκ, and
the antisymmetric part). The first vanishes identically iff the matrix-valued measure (Dr_μ − Dr_μᵀ)ρ is zero, by the
argument of (iii) applied entrywise with equal widths (the entries of ∫(Dr − Drᵀ)φ_cφ_c'dρ are Gaussian convolutions of
the entries of that measure). For (a): with equal widths, φ_ξφ_ξ' = e^{−|c−c'|²/(4s²)}e^{−|x−m|²/s²} and
∇φ_ξ' = −((x − c')/s²)φ_ξ', so

    κ(ξ,ξ') = −(1/s²)e^{−|c−c'|²/(4s²)}⟨e^{−|x_m−m|²/s²}(x_m − c')·ẋ_m⟩_t.

Write x_m − c' = (x_m − m) + (c − c')/2. The first part integrates to zero, being the time derivative of
−(s²/2)e^{−|x_m−m|²/s²} along a period. Hence

    κ(ξ,ξ') = −(1/(2s²))e^{−|c−c'|²/(4s²)}(c − c')·F_s(m),    F_s(m) := ⟨e^{−|x_m(t)−m|²/s²}ẋ_m(t)⟩_t = (1/T₀)∫e^{−|x−m|²/s²}dT_𝒪(x),

the Gaussian convolution of the vector measure T_𝒪. If κ ≡ 0 on admissible pairs, then F_s(m)·v = 0 for every m with
‖m‖ < c_max and every direction v (take c = m − δv, c' = m + δv), so F_s ≡ 0 on that ball; F_s is entire, so F_s ≡ 0 on
ℝⁿ, and the Fourier argument gives T_𝒪 = 0. Conversely, κ(ξ,ξ') = (1/T₀)∫φ_ξ∇φ_ξ'·dT_𝒪 = 0 for all widths when
T_𝒪 = 0. For n = 2, ∮_𝒪 v·dx = ∫w_𝒪 curl v dA for every smooth v (Stokes with winding numbers), and every smooth
compactly supported g is curl v for some smooth v, so T_𝒪 = 0 iff w_𝒪 = 0 a.e. (T_𝒪 is compactly supported, so v need
not be). For the example, r = he₁ with h = 0 on 𝒪 gives (Dr − Drᵀ)_{1l} = ∂_l h and the only non-zero entries are the
normal derivative; the Lie-bracket terms vanish because r = 0 on 𝒪.

(v) The saturated features are real-analytic in the lifted coordinates and extend holomorphically to a complex
neighbourhood of ℝ^p (tanh and the logistic function are holomorphic on the strips |Im| < πa_max/2, πc_max/2, π, where
s stays away from 0), uniformly for x on the compact orbit; integrals over the orbit of such integrands are real-analytic
in the parameters, and r_μ is affine in the weights. The zero set of a real-analytic function on a connected real-analytic
manifold is either the whole manifold or a closed set of Lebesgue measure zero. The instances are those of
`test_exactly_averaged_dynamics_is_non_reciprocal_for_n1_and_n2`.

(vi) Since r_ν is affine in ν, B_{ν_δ} = B_μ + δB^{even}[g] with g := f̂_μ − Φ̃_ϑ', where B^{even}[g] denotes the three
residual terms of B with g in place of r_μ. The trace of B^{even}[g] at equal widths is −(1/s²)(c − c')·∫gφ_cφ_c'dρ (the
curl term is traceless). If A_{ν_δ} ≡ 0 for all small δ, then B^{even}[g] ≡ 0, so (c − c')·F[gρ]((c + c')/2) = 0 for all
admissible c, c', where F[gρ] is the Gaussian convolution of the vector measure gρ; as in (iv), gρ = 0, that is g = 0
ρ-a.e. For ϑ' = ϑ_z this says f̂_μ = 0 on the orbit; for ϑ' with a' ≠ 0 it then says Φ̃_ϑ' = 0 on the orbit, which is
false since a Gaussian has no zeros. So one of the two segments has A_{ν_δ} ≢ 0 for all small δ > 0. ∎

*Checked* (`tests/test_nonreciprocity.py`, using `src/core/nonreciprocity.py` against `src/core/averaging.py`): the bilinear
form aᵀB_μa' reproduces the antisymmetric part of the generic second-variation coefficient (1e-10, n = 1 and 2);
A_μ(ϑ,ϑ_z) = 0 (1e-14) and the cyclic sum over (ϑ₀,ϑ₁,ϑ_z) equals A_μ(ϑ₀,ϑ₁); for n = 1 with f := f̂_μ (perfect fit) A_μ ≡ 0
(below 1e-13) while a non-zero residual gives A_μ ≠ 0, equal to the closed form; the closed form of κ (1e-10), the second
component of F_s(m) at m = (0.5, 0), s = 0.7 is positive on the ellipse, and κ ≡ 0 on a segment traversed back and forth
(below 1e-12); with f ≡ f̂_μ on ℝ² (so Dr_μ ≡ 0), B_μ = 2κI (1e-13), A_μ ≠ 0 on the ellipse and A_μ ≡ 0 on the segment;
with f := f̂_μ + he₁, h = (x₁/1.2)² + (x₂/0.8)² − 1 (zero on the ellipse, non-zero normal derivative), B_μ − 2κI equals
−∫φφ'∂₂h dρ·(e₁e₂ᵀ − e₂e₁ᵀ) (1e-12) and is not small; the positivity witness f ≡ 0, μ = δ_ϑ̄, for which
A_μ(ϑ₀,ϑ₁) = a₀a₁ā((c₀ − c₁)/s²)∫φ_c̄φ_c₀φ_c₁dρ has a strictly positive integral. The round-5 referee re-derived every
formula and checked aᵀB_μa' against the exact second variation of the averaged potential (Richardson in k_e) to 2e-4–9e-4
relative, including parameters at the saturation limits.

*Remark I.1 (what Theorem I changes in the reading of Proposition H.7).*

- The non-reciprocity of H.7(e) is no longer an instance-by-instance check. Pointwise, for n = 1 it holds at every
  frozen μ that does not interpolate f along the orbit; for n ≥ 2 at perfect fit it holds unless the orbit carries no
  circulation and the residual's Jacobian is symmetric along the orbit; in general it fails only on a null set of any
  connected analytic family that contains one non-reciprocal member. Globally, (vi) shows that no open set of measures
  is reciprocal, for any n: the averaged loop is nowhere a gradient flow at order k_e⁻².
- The two defects scale differently. The performative (semi-gradient) defect is bounded by the residual and the
  antisymmetric part of its Jacobian along the orbit, so it is small at a well-fitting equilibrium whose residual is
  smooth across the orbit; the circulation defect does not depend on the fit at all.
- For the n = 1 Langevin experiments and Conjecture H.6: the non-gradient correction that falsifies H.6 at fixed k_e
  is, by (iii), bounded by the residual, which at the regularised optimum μ*_τ is the regularisation bias. Both the
  leading drift and the correction are proportional to the residual there, so the correction is O(1/k_e) relative to
  the leading drift and O(k_e⁻²) absolute, as in Remark H.0.
- The pairwise criterion also simplifies the test of reciprocity: a single pair (ϑ, ϑ') with B_μ(ξ,ξ') ≠ 0 suffices.

## 6d. Theorem J — the finite-m closed loop as a perturbed gradient flow

This section proves the first step of §8, branch 1: at finite m, without dither, in the scaling σ = σ₀/k_e, the real
closed loop has, near every hyperbolic critical point of the m-particle regularised occupation loss, a unique periodic
orbit with the same stability type, at distance O(1/k_e), whose first-order location is the one Proposition H.7 and
Remark H.7a predict. It turns Remark H.7a's "formal" second-order averaging into a theorem for this system.

**Setting.** The closed loop of Theorem H's design with τ = 0, at finite m, in the scaling σ = σ₀/k_e:

    ẋ = f(x) + ẋ_m − k_e e − f̂_Θ(x),    e := x − x_m(t),    f̂_Θ := (1/m)Σ_j Φ̃(·;ϑ_j),
    ϑ̇_i = G∇_ϑ[eᵀΦ̃(x;ϑ_i)] − (σ₀/k_e)(ϑ_i − ϑ⁰),    i = 1..m,                                           (CL)

with Θ := (ϑ_1, …, ϑ_m) ∈ ℝ^{mp}, G_m := I_m ⊗ G, saturated features, x_m ∈ C² T₀-periodic with occupation measure ρ,
and f ∈ C² on a neighbourhood of the tube 𝒯_R := {x : dist(x, x_m([0,T₀])) ≤ R} for some R > 0. Put

    r_Θ(x) := f(x) − f̂_Θ(x),    J_m(Θ) := mL̄_m(Θ) + (σ₀/2)Σ_i‖ϑ_i − ϑ⁰‖²_{G⁻¹},    L̄_m(Θ) := ½∫‖r_Θ‖²dρ,
    H_i(Θ,t) := G∇_ϑ[r_Θ(x_m(t))ᵀΦ̃(x_m(t);ϑ)]|_{ϑ=ϑ_i} − σ₀(ϑ_i − ϑ⁰),    ⟨H⟩(Θ) = −G_m∇J_m(Θ).

The averaged system of (CL) is Θ̇ = (1/k_e)⟨H⟩(Θ) = −(1/k_e)G_m∇J_m(Θ): the gradient flow of J_m in the G⁻¹ metric,
which is (𝔄) at τ = 0 for the m-atom measure (Remark H.0); in slow time s := t/k_e it does not depend on k_e. Its
equilibria are the critical points of J_m. A critical point Θ* is *hyperbolic* if ∇²J_m(Θ*) is non-singular; its
index is the number of negative eigenvalues; it is a *non-degenerate minimum* if ∇²J_m(Θ*) ≻ 0. Write ν_1 ≤ … ≤ ν_{mp}
for the eigenvalues of G_m^{1/2}∇²J_m(Θ*)G_m^{1/2} (those of G_m∇²J_m(Θ*)), and λ_* := ν_1 when Θ* is a minimum.

Let u(Θ,t) be the zero-mean antiderivative in t of H(Θ,t) − ⟨H⟩(Θ) (T₀-periodic), and define, as in Proposition H.7 and
Remark H.7a, at the atomic measure μ* := (1/m)Σ_jδ_{ϑ*_j}:

    J := D_Θ[(1/k_e)⟨H⟩](Θ*) = −(1/k_e)G_m∇²J_m(Θ*),
    b_frozen,i := (1/k_e²)G∇_ϑ(U_even + U_odd)[μ*](ϑ)|_{ϑ=ϑ*_i},    b_osc := (1/k_e²)⟨D_ΘH(Θ*,t)·u(Θ*,t)⟩_t.

**Theorem J.** Let Θ* be a hyperbolic critical point of J_m. There are k₀, ρ₀, C > 0, depending on (f, x_m, G, σ₀, R,
the feature bounds, Θ*) but not on k_e, such that for every k_e ≥ k₀:

(a) *Existence, uniqueness, location.* (CL) has a T₀-periodic solution (x_p, Θ_p), unique among the T₀-periodic
solutions with sup_t‖Θ_p(t) − Θ*‖_{G⁻¹} ≤ ρ₀ and sup_t‖x_p(t) − x_m(t)‖ ≤ R. It satisfies

    sup_t‖x_p(t) − x_m(t)‖ ≤ C/k_e,    sup_t‖Θ_p(t) − Θ*‖_{G⁻¹} ≤ C/k_e,    Θ_p(t) − Θ̄_p = (1/k_e)u(Θ*,t) + O(k_e⁻²),

with Θ̄_p the period mean of Θ_p.

(b) *Stability type.* The monodromy matrix of (x_p, Θ_p) has n eigenvalues of modulus at most e^{−(k_e − L)T₀}, with
L independent of k_e (the error channel), and mp eigenvalues

    λ_j = exp(−(T₀/k_e)ν_j + O(k_e⁻²)),    j = 1..mp.

Hence the orbit is hyperbolic with the index of Θ*: exponentially stable if ∇²J_m(Θ*) ≻ 0, with parameter Floquet
exponents (1/k_e)ν_j + O(k_e⁻²), and unstable otherwise.

(c) *Basin (non-degenerate minimum).* If ∇²J_m(Θ*) ≻ 0, every solution of (CL) with ‖Θ(0) − Θ*‖_{G⁻¹} ≤ ρ₀ and
‖x(0) − x_m(0)‖ ≤ R converges to the periodic orbit, and

    ‖Θ(t) − Θ_p(t)‖_{G⁻¹} + ‖x(t) − x_p(t)‖ ≤ C e^{−μt}(‖Θ(0) − Θ_p(0)‖_{G⁻¹} + ‖x(0) − x_p(0)‖),    μ := (λ_*/k_e)(1 − κ),

with κ = κ(ρ₀, 1/k_e) → 0 as (ρ₀, 1/k_e) → 0. The basin radius ρ₀ does not depend on k_e; the rate is λ_*/k_e to
leading order, that is, λ_* in slow time.

(d) *First-order location.* The period mean satisfies

    Θ̄_p = Θ* − J⁻¹(b_frozen + b_osc) + O(k_e⁻²) = Θ* + (1/k_e)(G_m∇²J_m(Θ*))⁻¹·k_e²(b_frozen + b_osc) + O(k_e⁻²).

The offset is O(1/k_e). Its orientation-odd part is the quantity measured in Remark H.7a, so the prediction
Θ₊ − Θ₋ = −2J⁻¹(b_odd + b_osc) + O(k_e⁻²) of that remark is now a statement with error control.

*Proof.* Write ε := 1/k_e. All constants below are uniform in ε ≤ 1/k₀ and in Θ in the ball 𝒦 := {‖Θ − Θ*‖ ≤ 2ρ₀},
because the saturated features and their derivatives are bounded uniformly in Θ and f ∈ C² on the compact tube 𝒯_R.

*Step 1 (the tube and the slaved error).* Let M_r := sup{‖r_Θ(x)‖ : x ∈ 𝒯_R, Θ ∈ 𝒦} and L_r := sup‖∂_x r_Θ‖ on the
same set. For k_e > M_r/R the tube {‖e‖ ≤ R} is forward invariant: on ‖e‖ = R, d/dt ½‖e‖² = eᵀr_Θ(x) − k_e‖e‖² < 0.
Put w := e − εr_Θ(x) with x the actual state. Since r_Θ(x) − k_e e = −k_e w and ẋ = ẋ_m − k_e w,

    ẇ = −(k_e I − ∂_x r_Θ(x))w − ε∂_x r_Θ(x)ẋ_m − ε∂_Θ r_Θ(x)·Θ̇,

and Θ̇ = F(e,Θ,t) := (G∇_ϑ[eᵀΦ̃(x;ϑ_i)] − εσ₀(ϑ_i − ϑ⁰))_i is O(‖e‖ + ε) = O(‖w‖ + ε). So for k_e ≥ 2L_r + C,
d/dt‖w‖ ≤ −(k_e/2)‖w‖ + C₁ε on {‖w‖ ≤ 1} ∩ {‖e‖ ≤ R}, whence ‖w(t)‖ ≤ e^{−k_e t/2}‖w(0)‖ + 2C₁ε². The set
W_ε := {‖w‖ ≤ 3C₁ε²} is forward invariant and absorbs every solution of the tube within t₁ = O(ε log(1/ε)), during
which Θ moves by O(ε‖e(0)‖ + ε² log(1/ε)) (Θ̇ = O(‖e‖ + ε) and ‖e‖ decays at rate k_e/2). On W_ε, e = εr_Θ(x_m(t)) + ε[r_Θ(x_m + e) − r_Θ(x_m)] + w = εr_Θ(x_m(t)) + O(ε²).

*Step 2 (the slow equation and the fast variable).* Let e_Θ(t) be the frozen-Θ periodic error of Proposition H.7(a):
the unique T₀-periodic solution of ė = r_Θ(x_m(t) + e) − k_e e, exponentially stable with rate at least k_e/2, with
e_Θ = εr_Θ(x_m) + ε²(∂_x r_Θ r_Θ − ṙ_Θ)(x_m) + O(ε³) in sup norm (Taylor expansion of r_Θ in x to second order; f ∈ C²).
Because Θ enters r_Θ = f − f̂_Θ only through the smooth f̂_Θ, e_Θ(t) is C² in Θ as soon as f ∈ C¹, and differentiating
the frozen equation (gain 1/k_e) gives ∂_Θ e_Θ = ε∂_Θ r_Θ(x_m) + O(ε²) and ∂²_Θ e_Θ = O(ε), uniformly on 𝒦. Let
ζ := e − e_{Θ(t)}(t) be the fast variable. Then

    ζ̇ = [r_Θ(x_m + e) − r_Θ(x_m + e_Θ)] − k_eζ − ∂_Θ e_Θ(t)·Θ̇,

where the bracket is Lipschitz in ζ with constant L_r and the last term is O(ε)·O(ε) = O(ε²) (Θ̇ = F = O(‖e‖ + ε)). So
for k_e ≥ 2L_r + C the set Z_ε := {‖ζ‖ ≤ C₂ε³} is forward invariant and absorbs every solution of the tube within
O(ε log(1/ε)); on W_ε of Step 1 one has ‖ζ‖ = O(ε²) at once and O(ε³) after this further transient. On Z_ε, since F is
C¹ in e with ∂_eF = O(1),

    Θ̇ = F(e_Θ(t), Θ, t) + R(Θ, ζ, t) = εH(Θ,t) + ε²H₂(Θ,t) + O(ε³),    ‖R‖ ≤ C‖ζ‖ = O(ε³),                      (S)

where εH + ε²H₂ is the expansion of the frozen-Θ drift along the frozen periodic error, and by Proposition H.7(b),
⟨H₂⟩(Θ)_i = G∇_ϑ(U_even + U_odd)[μ_Θ](ϑ_i), that is, ε²⟨H₂⟩(Θ*) = b_frozen. H and H₂ are T₀-periodic and C² in Θ on 𝒦,
with bounds uniform in ε.

*Regularity of the remainders (used in Steps 4 and 6).* The second-order terms of (S), and of (N) below, are used in sup
norm only, for the location (d). The implicit-function and contraction arguments of Steps 4 and 6 need first-order C¹
closeness, which holds: F(e_Θ(t), Θ, t) is C² in Θ with ∂_Θ[F(e_Θ(t), Θ, t)] = ε∂_ΘH(Θ, t) + O(ε²) (from
∂_Θ e_Θ = ε∂_Θ r_Θ(x_m) + O(ε²) and the smooth dependence of F on (e, Θ) through Φ̃); the remainder R is Lipschitz in ζ with
constant sup‖∂_eF‖ = O(1) and Lipschitz in Θ with constant O(‖ζ‖) = O(ε³); and the forcing −∂_Θ e_Θ·Θ̇ of the ζ-equation
is Lipschitz in Θ with constant O(ε²) (from ∂²_Θ e_Θ = O(ε) and ∂_ΘF = O(ε)). None of this uses more than f ∈ C².

*Step 3 (second-order averaging normal form).* Let u be the zero-mean antiderivative of H − ⟨H⟩ and Θ = Y + εu(Y,t).
Then Ẏ = ε⟨H⟩(Y) + ε²[H₂(Y,t) + D_ΘH(Y,t)u(Y,t) − D_Yu(Y,t)⟨H⟩(Y)] + O(ε³). Let u₂ be the zero-mean antiderivative of the
fluctuating part of the bracket and Y = Y₂ + ε²u₂(Y₂,t). Since ⟨D_Yu⟩ = D_Y⟨u⟩ = 0,

    Ẏ₂ = ε⟨H⟩(Y₂) + ε²Φ₂(Y₂) + O(ε³),    Φ₂ := ⟨H₂⟩ + ⟨D_ΘH·u⟩,                                              (N)

the second-order averaged equation of Remark H.7a (Sanders–Verhulst–Murdock §2.9), now with a remainder that is
uniformly O(ε³) on Z_ε × 𝒦 because (S) holds there. The transformation Θ ↦ Y₂ is a T₀-periodic near-identity
diffeomorphism from 𝒦 onto its image, O(ε)-close to the identity with O(ε) derivative, and the period mean of Θ along any curve
Y₂(t) = Y₂⁰ + O(ε²) equals the period mean of Y₂ up to O(ε³), because u and u₂ have zero mean.

*Step 4 (the Poincaré map in the slow coordinates).* Let P = (P_Y, P_ζ) be the time-T₀ map of (CL) from t = 0, written
in the coordinates (Y₂, ζ) on 𝒦 × Z_ε. By (N) and Step 2,

    P_Y(Y₂, ζ) = Y₂ + εT₀ Ã(Y₂) + O(ε³),    Ã := ⟨H⟩ + ε[Φ₂ + (T₀/2)D⟨H⟩·⟨H⟩],    ‖P_ζ(Y₂, ζ)‖ ≤ e^{−k_eT₀/2}‖ζ‖ + O(ε³),

the ε²(T₀²/2)D⟨H⟩⟨H⟩ term being the second-order term of the time-T₀ flow of (N). The Lipschitz constants of the coupling
are Lip_ζ P_Y = O(ε) (a perturbation of ζ acts on Θ̇ through ∂_eF = O(1) while it decays at rate k_e/2, so its integrated
effect over a period is O(ε)‖ζ‖) and Lip_Y P_ζ = O(ε³) (the forcing of the ζ-equation is O(ε²)-Lipschitz in Θ and is
integrated against the fast decay). Hence ζ ↦ P_ζ(Y₂, ζ) has a unique fixed point ζ*(Y₂) = O(ε³), Lipschitz in Y₂ with
constant O(ε³), and the fixed points of P are the zeros of

    Ã_ε(Y) := (P_Y(Y, ζ*(Y)) − Y)/(εT₀) = Ã(Y) + O(ε²).

By the regularity paragraph, D_Y P_Y = I + εT₀D⟨H⟩(Y) + O(ε²) (the variational equation of the slow flow over one period,
whose coefficient is ε∂_ΘH + O(ε²)), and the ζ*-dependence contributes Lip_ζ P_Y·Lip_Y ζ* = O(ε⁴), so Ã_ε is Lipschitz with
Lip(Ã_ε − ⟨H⟩) = O(ε) on 𝒦. Since D⟨H⟩(Θ*) = −G_m∇²J_m(Θ*) is invertible, Y ↦ Y − D⟨H⟩(Θ*)⁻¹Ã_ε(Y) maps a ball of
ε-independent radius ρ₀ around Θ* into itself and is a contraction there for ε small (the Lipschitz form of the implicit
function theorem), which gives a unique zero Y₂* = Θ* + O(ε), hence a unique fixed point of P with Θ-component within
O(ε) of Θ*, that is, a unique T₀-periodic solution in the stated set; the bound on x_p − x_m is Step 1. Expanding Ã_ε(Y₂*) = 0 to first order,
D⟨H⟩(Θ*)(Y₂* − Θ*) + εΦ₂(Θ*) = O(ε²) (the term (T₀/2)D⟨H⟩⟨H⟩(Y₂*) is O(ε) because ⟨H⟩(Y₂*) = O(ε)), so

    Y₂* = Θ* − εD⟨H⟩(Θ*)⁻¹Φ₂(Θ*) + O(ε²) = Θ* − J⁻¹(b_frozen + b_osc) + O(ε²),

using J = εD⟨H⟩(Θ*) and ε²Φ₂(Θ*) = b_frozen + b_osc. Along the periodic orbit Y₂(t) = Y₂* + O(ε²) (since
⟨H⟩(Y₂*) = O(ε)), so by Step 3 the period mean of Θ_p is Y₂* + O(ε²): this is (d), and Θ_p(t) − Θ̄_p = εu(Θ*,t) + O(ε²)
is (a).

*Step 5 (the monodromy matrix).* Linearise (CL) along the periodic orbit: with δe and δΘ,

    δė = (∂_x r − k_e I)δe + ∂_Θ r·δΘ,    δΘ̇ = ∂_eF·δe + ∂_ΘF·δΘ,

all coefficients T₀-periodic and evaluated along the orbit; ∂_x r, ∂_Θ r, ∂_eF are O(1) and ∂_ΘF = O(ε) (F is linear
in e with coefficients smooth in Θ, plus the leak −εσ₀(ϑ_i − ϑ⁰), so ∂_ΘF = O(‖e‖ + ε), and e = O(ε) on the orbit). Let S(t) be a T₀-periodic solution of the Riccati equation

    Ṡ = (∂_x r − k_e I)S + ∂_Θ r − S∂_eF S − S∂_ΘF.

For k_e large it exists, is unique with ‖S‖ ≤ C/k_e, and depends smoothly on the data (contraction of the periodic
solution operator of the dominant part Ṡ = −k_eS + g; cf. Proposition H.7(a)). In the coordinates η := δe − SδΘ the
linearisation is block triangular:

    η̇ = [(∂_x r − k_e I) − S∂_eF]η,    δΘ̇ = ∂_eF·η + [∂_eF S + ∂_ΘF]δΘ,

so the monodromy matrix is similar to a block-triangular matrix whose diagonal blocks are the monodromies of the two
diagonal equations. The η-block has norm at most e^{−(k_e − L)T₀} with L := L_r + sup‖S∂_eF‖ = O(1) (Grönwall): these
are the n error-channel multipliers. The δΘ-block is the linearisation along the orbit of the slow equation (S) in
Θ-coordinates, εD_ΘH(Θ_p(t),t) + O(ε²) = εD_ΘH(Θ*,t) + O(ε²), whose monodromy is I + εT₀D⟨H⟩(Θ*) + O(ε²) =
I − εT₀G_m∇²J_m(Θ*) + O(ε²). Its eigenvalues are 1 − εT₀ν_j + O(ε²) = exp(−εT₀ν_j + O(ε²)), those of G_m∇²J_m(Θ*)
being the ν_j by similarity with G_m^{1/2}∇²J_m(Θ*)G_m^{1/2}; the O(ε²) perturbation moves them by O(ε²) because
I − εT₀G_m∇²J_m(Θ*) is diagonalisable (similar to a symmetric matrix through G_m^{1/2}), so the Bauer–Fike theorem applies
with a condition number independent of ε — this covers the repeated eigenvalues of a collapsed configuration. This is (b).

*Step 6 (the basin).* Let ∇²J_m(Θ*) ≻ 0. In the norm ‖v‖_{G⁻¹} := ‖G_m^{−1/2}v‖, the matrix I − εT₀G_m∇²J_m(Θ*) is
similar to a symmetric contraction with norm exactly 1 − εT₀λ_* once εT₀ν_{mp} < 1. By Step 4,
D_Y P_Y = I − εT₀G_m∇²J_m(Y) + O(ε²) on 𝒦, and ∇²J_m is Lipschitz, so on the ball {‖Y₂ − Y₂*‖_{G⁻¹} ≤ ρ₀} the Lipschitz
constant of P_Y in Y₂ is at most 1 − εT₀λ_*(1 − κ) with κ = O(ρ₀ + ε); the coupling constants are Lip_ζ P_Y = O(ε) and
Lip_Y P_ζ = O(ε³), and P_ζ contracts ζ by e^{−k_eT₀/2} (Step 4). The 2 × 2 matrix of Lipschitz constants
[[1 − εT₀λ_*(1 − κ), Cε], [Cε³, e^{−k_eT₀/2}]] has spectral radius 1 − εT₀λ_*(1 − κ′) with κ′ = O(ρ₀ + ε), so on
{‖Y₂ − Y₂*‖_{G⁻¹} ≤ ρ₀} × Z_ε, with ρ₀ small but independent of ε, P is a contraction in an adapted product norm with
that factor: its iterates converge to the fixed point geometrically, uniformly in ε. Solutions started in the tube with
‖e(0) − e_p(0)‖ ≤ 1 enter Z_ε within O(ε log(1/ε)), during which Θ moves by O(ε)‖e(0) − e_p(0)‖ (Step 1), which is absorbed
by halving ρ₀. Between the sampling times the solutions stay within a factor C of each other
(Grönwall on one period, with the fast variable absorbing ‖e(0) − e_p(0)‖ into an O(ε)‖e(0) − e_p(0)‖ displacement of
Θ, Step 1). Transporting back through the near-identity transformations gives (c). ∎

*Numerics* (`exp_periodic_orbit.py` → `results/periodic_orbit/table.md`; `tests/test_periodic.py`). Two instances, k_e = 10, 20,
40, 80, 160. For each k_e the periodic orbit of (CL) is found by Newton shooting on the Poincaré map (DOP853 at relative
tolerance 1e-11; residual below 10⁻¹¹ in every case), its monodromy matrix is the finite-difference Jacobian of the map, and
its period mean is compared with Θ* and with the first-order prediction Θ* − J⁻¹(b_frozen + b_osc) of (d). Distances and
cosines are in the labelled G⁻¹ metric.

- **Instance A** (n = 2, m = 3, the plant and elliptic reference of Remark H.7a, σ₀ = 0.5, G = diag(4,4,1,1,1)). The
  minimiser of J_m is a collapsed configuration (kernel spread 4·10⁻¹⁴) with λ_min(∇²J_m) = 0.0148: a non-degenerate
  minimum with permutation symmetry.

  | k_e | ‖Θ̄_p − Θ*‖ | full prediction: size ratio, cosine, vector error | frozen-measure part only: vector error | parameter Floquet exponents against (1/k_e)ν_j: max, median relative error |
  |---|---|---|---|---|
  | 10 | 4.09e-2 | 1.149, 0.9972, 0.169 | 0.398 | 8.3e-2, 1.2e-2 |
  | 20 | 1.91e-2 | 1.074, 0.9993, 0.083 | 0.360 | 7.3e-2, 7.3e-3 |
  | 40 | 9.22e-3 | 1.036, 0.9998, 0.041 | 0.348 | 4.5e-2, 4.0e-3 |
  | 80 | 4.53e-3 | 1.018, 1.0000, 0.020 | 0.343 | 2.5e-2, 2.1e-3 |
  | 160 | 2.24e-3 | 1.009, 1.0000, 0.010 | 0.342 | 1.3e-2, 1.1e-3 |

  The offset scales like k_e^{−1.05} (fitted exponent; (a) predicts −1), the vector error of the full prediction
  halves with every doubling of k_e (the O(k_e⁻²) remainder of (d) against an O(1/k_e) offset), and the frozen-measure
  part alone has a 35% vector error at every k_e: the oscillation term b_osc of Remark H.7a is an O(1) fraction of the
  first-order location. The median relative error of the Floquet exponents is O(1/k_e), as (b) states; the largest
  relative error sits on the (m − 1)-fold degenerate separation pair with the smallest exponent (0.0018 at k_e = 20),
  where an O(k_e⁻²) absolute error is a large fraction of a small exponent. The two error-channel multipliers are
  below the resolution of the finite-difference monodromy (exponents of order k_e); their product is recovered from
  Liouville's identity Σ log|λ| = ∫₀^{T₀}tr ∂F/∂z dt along the orbit, which gives a mean error-channel exponent of
  k_e − 0.16, k_e − 0.11, k_e − 0.086, k_e − 0.074, k_e − 0.068 at k_e = 10 … 160 (instance B: k_e + 0.05 … k_e + 0.09):
  the multipliers are e^{−(k_e − L)T₀} with |L| < 0.2, as (b) states.

- **Instance B** (n = 1, m = 2, a two-bump plant f(x) = 1.6e^{−(x−1)²/0.4} − 1.2e^{−(x+1)²/0.4} + 0.1x on x_m = 2 sin t,
  σ₀ = 0.05, G = diag(2,1,1)). The minimiser has two distinct kernels of opposite sign (spread 2.1, λ_min(∇²J_m) = 0.066).
  Here b_osc = 0 exactly: for a one-dimensional sinusoidal sweep the reflection t ↦ π − t maps x_m to itself and u to
  −u, so ⟨D_ΘH·u⟩ = 0, and the full and frozen predictions coincide. The reference is fast on the kernel scale
  (sup‖ẋ_m‖/s is 2/s: between 1.3 and 10 over the admissible widths, about 2.6 at the minimiser's width; the residual bandwidth is of the same order), so k_e = 10–20 is outside the averaging regime:

  | k_e | ‖Θ̄_p − Θ*‖ | prediction: size ratio, cosine, vector error | Floquet: max, median relative error |
  |---|---|---|---|
  | 10 | 6.44e-2 | 2.818, 0.880, 1.995 | 9.5e-2, 4.0e-2 |
  | 20 | 2.22e-2 | 1.942, 0.929, 1.078 | 3.2e-2, 1.4e-2 |
  | 40 | 8.37e-3 | 1.464, 0.970, 0.551 | 1.2e-2, 5.4e-3 |
  | 80 | 3.51e-3 | 1.227, 0.990, 0.277 | 4.6e-3, 2.2e-3 |
  | 160 | 1.59e-3 | 1.111, 0.9969, 0.139 | 2.0e-3, 9.5e-4 |

  The vector error and the Floquet error halve with each doubling of k_e, as the theorem's O(k_e⁻²) remainders
  predict, but the constants are large: the asymptotic regime starts around k_e ≈ 80 here. This instance shows what
  the theorem does not claim: the O(1/k_e) location is sharp only once k_e exceeds the effective bandwidth ω_eff of
  Proposition G's (H2a′).

*Checked* (`tests/test_periodic.py`): ∇J_m against finite differences (1e-7) and the non-degenerate two-kernel minimum of
instance B; at k_e = 20 and 40 the shooting residual is below 1e-10, ‖Θ̄_p − Θ*‖ and max‖e‖ fall by more than 35% from 20
to 40, the vector error of the first-order location falls by more than 40% and is below 0.6 at k_e = 40, and the maximal
relative Floquet error falls by more than 40% and is below 2% at k_e = 40.

*Remark J.1 (what Theorem J says about the program of §8).*

- It is the deterministic, finite-m first step of branch 1: near every hyperbolic critical point of the m-particle
  regularised loss k_e·𝒢 at τ = 0 the real closed loop has exactly one stationary state (a periodic orbit), an
  O(1/k_e)-perturbation of it with the same stability type, and the convergence rate is the one of the averaged gradient
  flow in slow time. (The theorem does not say that every periodic orbit of (CL) arises this way.) The two missing
  pieces of Conjecture H.6 remain: dither (τ > 0, where the stationary state is a measure) and uniformity in m.
- It makes the stationary (period-mean) prediction of Remark H.7a a theorem for this system: the first-order location
  of the orbit, hence the orientation difference of the period-mean kernels, is Θ* − J⁻¹(b_frozen + b_osc) with an
  O(k_e⁻²) error. When the leading offset is itself small, the O(k_e⁻²) remainder can dominate its *direction* at
  moderate k_e: a collapsed saddle instance examined in the round-5 review (`review/REVIEW_round5_J.md`) has a
  first-order location nearly orthogonal to the true offset at k_e = 20 and aligned with it at k_e = 80, while the sizes
  obey the theorem's orders throughout.
- The leak is not needed for the structure, only for the hyperbolicity of Θ*. At a collapsed critical point the
  separation Hessian is σ₀G⁻¹ − R₂ (Proposition K); for n = 1 and σ₀ = 0 the amplitude row of R₂ vanishes at a leak-free
  critical point (⟨rφ⟩ = ∂_c⟨rφ⟩ = ∂_s⟨rφ⟩ = 0 there), so the m − 1 amplitude-redistribution directions are exactly
  flat and the theorem does not apply, while for n ≥ 2 only the amplitude block of R₂ vanishes and a collapsed critical
  point is generically non-degenerate even without a leak. Where ∇²J_m is singular the theorem does not apply.
- At a collapsed critical point (all ϑ*_i equal) the permutation symmetry makes the separation modes (m − 1)-fold
  degenerate; the Floquet exponents come in such multiplets (`results/periodic_orbit/table.md`); Proposition K (§6e) works out this case.

---

## 6e. Proposition K — collapsed configurations, the critical leak, and the orientation of the orbit as a design input

Branch 2 of §8 asks whether the circulation of the reference orbit can be designed. Theorem J makes the question
precise for the collapsed configuration, the trap that Theorem H's dither was introduced to escape: its stability is
decided by the eigenvalues of a Hessian that the orientation shifts at order k_e⁻².

**Setting.** As in Theorem J. A configuration is *collapsed* if ϑ_i = ϑ̄ for all i. Then f̂_Θ = Φ̃(·;ϑ̄), so Θ is a critical
point of J_m if and only if ϑ̄ is a critical point of the single-kernel objective J_1(ϑ) := ½∫‖f − Φ̃(·;ϑ)‖²dρ +
(σ₀/2)‖ϑ − ϑ⁰‖²_{G⁻¹}. Write r̄ := f − Φ̃(·;ϑ̄) and

    M₂ := ∫∇_ϑΦ̃(·;ϑ̄)ᵀ∇_ϑΦ̃(·;ϑ̄)dρ ⪰ 0,    R₂ := ∫Σ_k r̄_k∇²_ϑΦ̃_k(·;ϑ̄)dρ = ∇²_ϑ⟨r̄ᵀΦ̃(·;ϑ)⟩_ρ|_{ϑ̄}  (residual frozen).

R₂ is the Hessian of −V_L[μ̄](ϑ) at the atom: the second-order cost of moving a small amount of mass off the atom at
frozen residual.

**Proposition K.** Let Θ* be a collapsed critical point of J_m with atom ϑ̄.

(i) *Symmetry splitting.* ∇²J_m(Θ*) = I_m ⊗ (σ₀G⁻¹ − R₂) + (1/m)(1 1ᵀ) ⊗ M₂. On the separation subspace
{(u_1..u_m) : Σ_i u_i = 0}, of dimension (m − 1)p, it acts as σ₀G⁻¹ − R₂; on the collective subspace {u_i all equal} as
∇²J_1(ϑ̄) = σ₀G⁻¹ − R₂ + M₂. The same splitting holds for the Jacobian of any permutation-equivariant particle drift
at a collapsed configuration, in particular for the second-order averaged drift of Remark H.7a: it has the form
I_m ⊗ D + (1/m)(1 1ᵀ) ⊗ E.

(ii) *Separation stability at leading order.* Θ* is a non-degenerate local minimum of J_m iff ∇²J_1(ϑ̄) ≻ 0 and
σ₀G⁻¹ ≻ R₂, that is, iff ϑ̄ is a non-degenerate minimiser of J_1 and σ₀ > σ_crit(ϑ̄) := λ_max(G^{1/2}R₂G^{1/2}). The
separation Floquet exponents of the collapsed periodic orbit of (CL) are then, by Theorem J(b),
(1/k_e)λ_j(σ₀I − G^{1/2}R₂G^{1/2}) + O(k_e⁻²), each with multiplicity m − 1, and they do not depend on m. The condition
σ₀G⁻¹ ≻ R₂ is the statement that the atom is a strict local minimiser of the first-variation potential V_L[μ̄] + (σ₀/2)‖·‖²_{G⁻¹}
on its support, the local form of the Chizat–Bach optimality condition for measures.

(iii) *The orientation-odd shift.* Let 𝒪 be traversed in either direction (the same ρ, hence the same Θ* and the
same leading exponents). Write the second-order normal form of Theorem J's proof as

    (N)    Θ̇ = −(1/k_e)G_m∇J_m(Θ) + b(Θ) + O(k_e⁻³),    b := b_frozen + b_osc = O(k_e⁻²),

with b_frozen = (1/k_e²)G_m∇(U_even + U_odd) (Proposition H.7) and b_osc of Remark H.7a. (N) has a unique equilibrium Θ̄₂
near Θ*, which is collapsed, with Θ̄₂ − Θ* = −J₁⁻¹b(Θ*) + O(k_e⁻²) for J₁ := −(1/k_e)G_m∇²J_m(Θ*) (the first-order location
of Theorem J(d)); let J₂ be the Jacobian of (N) at Θ̄₂. The parameter Floquet exponents of the collapsed orbit are the
eigenvalues of −J₂ up to O(k_e⁻³), and

    J₂ = J₁ + ∂b(Θ*) + ∂J₁(Θ*)[Θ̄₂ − Θ*] + O(k_e⁻³),

where the third term — the third derivative of J_m contracted with the first-order offset — is O(k_e⁻²), of the same order
as ∂b, and must be kept: the Jacobian of (N) at Θ* is not the right linearisation. Under orbit reversal ρ, Θ*, J₁ and
U_even are unchanged while U_odd and b_osc change sign (Proposition H.7(d), Remark H.7a), so with b_even := (1/k_e²)G_m∇U_even
and b_odd := (1/k_e²)G_m∇U_odd + b_osc,

    J₂ = J₁ + E ± O,    E := ∂b_even(Θ*) − ∂J₁(Θ*)[J₁⁻¹b_even(Θ*)],    O := ∂b_odd(Θ*) − ∂J₁(Θ*)[J₁⁻¹b_odd(Θ*)],

for the two orientations. J₁, E and O are permutation-equivariant and preserve the separation subspace; by (i) their
restrictions to one separation block are p × p matrices, written J₁|_{sep} = −(1/k_e)(σ₀I − GR₂), D_even, D_odd. The
separation Floquet exponents of the two orientations differ by

    λ_j^{(+)} − λ_j^{(−)} = −2 wⱼᵀD_odd vⱼ / (wⱼᵀvⱼ) + O(k_e⁻³)    (simple eigenvalue λ_j of −J₁|_{sep}, right and left eigenvectors vⱼ, wⱼ),

an O(k_e⁻²) quantity that scales like k_e⁻² at fixed (G, σ₀), while the exponents themselves are O(1/k_e). The
difference is a property of the excitation's orientation alone: the occupation measure, the averaged system (𝔄) and
all of Theorems B–H are blind to it.

(iv) *The critical leak as a design point.* At σ₀ = σ_crit(ϑ̄) the softest separation exponent vanishes at leading order
and equals the k_e⁻² value −(wᵀD_even v)/(wᵀv) ∓ (wᵀD_odd v)/(wᵀv) + O(k_e⁻³). If |wᵀD_odd v| > |wᵀD_even v|, the
collapsed configuration is stable for one orientation and unstable for the other: the direction in which the reference
is traversed decides whether the kernels separate. More generally, for σ₀ in a window of width O(1/k_e) around σ_crit —
the window where the leading exponent (σ₀ − σ_crit)(1 + o(1))/k_e lies between the two k_e⁻² values — the orientation
decides the sign, and the centre of the window, σ₀ = σ_crit + k_e(wᵀD_even v)/(wᵀv) + O(k_e⁻¹), can be computed from the
averaged data alone.

*Status.* (i) and (ii) are proved (linear algebra at the collapsed point and Theorem J(b)). (iii) and (iv) are proved
given the second-order normal form (N) with an O(k_e⁻³) remainder, which is Theorem J's proof carried one order further
(periodic averaging to second order [Sanders–Verhulst–Murdock Ch. 2–3]; the first-order step is written out in the proof
of Theorem J, the second-order step is Remark H.7a's derivation and is not written out here). The orientation-odd shift,
its k_e⁻² scaling, the prediction and the design window are checked numerically in one instance
(`exp_orientation_stability.py`).

*Proof.* (i) ∂²J_m/∂ϑ_i∂ϑ_j = (1/m)∫∇Φ̃_iᵀ∇Φ̃_j dρ − δ_ij∫Σ_k r_k∇²Φ̃_{i,k} dρ + δ_ijσ₀G⁻¹ (the first term from the
quadratic dependence of ½∫‖r‖²dρ on the m features with weights 1/m, the second from the ϑ_i-curvature of the
residual term: with r = f − (1/m)ΣΦ̃_j, ∂_{ϑ_i}(mL̄_m) = −∫rᵀ∇Φ̃_i dρ and ∂_{ϑ_j}∂_{ϑ_i}(mL̄_m) = (1/m)∫∇Φ̃_jᵀ∇Φ̃_i dρ − δ_ij∫Σ_k r_k∇²Φ̃_{i,k} dρ).
At a collapsed point all blocks coincide, giving the Kronecker form; the two subspaces are the isotypic components of
the permutation action, on which (1/m)(1 1ᵀ) acts as 0 and 1 respectively. The second-order averaged drift is
permutation-equivariant (each particle moves by a function of its own position and of the measure), so its Jacobian at
a collapsed point has the same form with D the self-derivative and E the derivative through the measure.

(ii) Positive definiteness of a block-diagonal form; the exponents are Theorem J(b) applied to the separation block,
whose eigenvalues of G_m^{1/2}∇²J_mG_m^{1/2} restricted to the separation subspace are those of σ₀I − G^{1/2}R₂G^{1/2}, each
repeated m − 1 times. For the last sentence, the first variation of 𝒢 := L̄ + (σ₀/2)∫‖·−ϑ⁰‖²_{G⁻¹}dμ at μ̄ = δ_ϑ̄ is
V_L[μ̄](ϑ) + (σ₀/2)‖ϑ − ϑ⁰‖²_{G⁻¹} = −⟨r̄ᵀΦ̃(·;ϑ)⟩_ρ + (σ₀/2)‖ϑ − ϑ⁰‖²_{G⁻¹}, whose Hessian at ϑ̄ is −R₂ + σ₀G⁻¹.

(iii) The proof of Theorem J carries one order further: the near-identity change of variables of its Step 3 removes the
oscillating terms through order k_e⁻², which turns (CL) into (N) with an O(k_e⁻³) remainder (this is the normal form of
Remark H.7a's derivation, with b_osc the term that the frozen-measure expansion misses). (N) is a regular perturbation
of the averaged system near the hyperbolic point Θ*, so it has a unique equilibrium Θ̄₂ within O(1/k_e) of Θ*, obtained
by the implicit function theorem; by the first-order expansion, Θ̄₂ − Θ* = −J₁⁻¹b(Θ*) + O(k_e⁻²) (J₁⁻¹ = O(k_e), b = O(k_e⁻²)).
(N) is permutation-equivariant and Θ* is collapsed, so the equivariant Newton iterates stay collapsed and Θ̄₂ is
collapsed. Theorem J(b) applied to (N) (same proof, remainder O(k_e⁻³)) gives the Floquet exponents as the eigenvalues
of −J₂ + O(k_e⁻³). Expanding the Jacobian of (N) about Θ*: J₂ = J₁ + ∂b(Θ*) + ∂J₁(Θ*)[Θ̄₂ − Θ*] + O(‖Θ̄₂ − Θ*‖²/k_e) +
O(k_e⁻²‖Θ̄₂ − Θ*‖), and both error terms are O(k_e⁻³); substituting the offset gives E and O, which are the even and odd
parts under reversal because ∂J₁ is even and b_even, b_odd have the parities named. All three operators are
permutation-equivariant, hence preserve the separation subspace, and by the Kronecker structure of (i) each acts on it
as I_{m−1} ⊗ (a p × p block). First-order eigenvalue perturbation of the simple eigenvalue λ_j of −J₁|_{sep} gives
λ_j^{(±)} = λ_j − wⱼᵀD_even vⱼ/(wⱼᵀvⱼ) ∓ wⱼᵀD_odd vⱼ/(wⱼᵀvⱼ) + O(k_e⁻³), and the difference as stated. (For the
(m − 1)-fold degenerate separation modes the degeneracy is the trivial one of the Kronecker structure, which the
restriction to one p-dimensional block removes; a further degeneracy within the block would call for degenerate
perturbation theory.)

(iv) is (iii) at a vanishing leading eigenvalue; the window follows because the leading exponent is
(σ₀ − σ_crit)/k_e·(1 + o(1)) near σ_crit (ϑ̄, R₂ and the eigenvectors depend smoothly on σ₀) while the shifts are
O(k_e⁻²); its centre is where the leading exponent cancels the even shift. ∎

*Remark K.1 (reading).* Three things.

- The collapsed configuration's stability is a single-kernel question at leading order (the Hessian σ₀G⁻¹ − R₂ does
  not depend on m), which explains why the m = 6 and m = 32 experiments of §6c collapsed alike: a leak above σ_crit
  traps every network size.
- Dither (Theorem H) and the orientation (this proposition) are two different ways to leave the trap: the first adds
  entropy, the second tilts the Hessian by O(k_e⁻²). The second is free but small; it matters only near σ_crit, that is,
  when the leak is tuned to the edge of trapping, or for nearly flat separation modes such as the weak saddle of
  Remark H.2.
- The sign of the shift depends on the data (f, 𝒪, ϑ̄) through D_odd and has no general sign: the design rule is to
  compute D_odd (from the averaged data: U_odd, b_osc and the third derivative of J_m) and choose the orientation with the
  smaller softest separation exponent (or the negative one); the window of (iv) is where the choice matters.

*Numerics* (`exp_orientation_stability.py` → `results/orientation_stability/table.md`). Instance A of Theorem J's
numerics (n = 2, m = 3, the elliptic reference traversed in both directions, G = diag(4,4,1,1,1)). For each orientation
the collapsed periodic orbit is found by Newton shooting (residual ≤ 4·10⁻¹³), the parameter block of its monodromy
matrix is restricted to the separation and collective subspaces, and the exponents are compared with the leading order
(σ₀I − G^{1/2}R₂G^{1/2})/k_e and with the second-order prediction, the eigenvalues of −J₂ at the equilibrium Θ̄₂ of (N)
(Newton from Θ*; every Θ̄₂ found is collapsed to 10⁻⁸). The second-order equilibrium is ‖Θ̄₂ − Θ*‖_{G⁻¹} = 3.7·10⁻²,
1.8·10⁻², 9.0·10⁻³, 4.5·10⁻³ away from Θ* at k_e = 10, 20, 40, 80: the O(1/k_e) offset of Theorem J(d).

- *The shift and its prediction* (σ₀ = 0.5, a stable collapse). The softest separation exponents of the two orientations
  differ by 9.0·10⁻⁴, 1.9·10⁻⁴, 4.1·10⁻⁵, 9.7·10⁻⁶ at k_e = 10, 20, 40, 80: fitted exponent −2.18 against k_e ((iii):
  −2). Scaled by k_e², the measured corrections to the leading order are −0.032, −0.057, −0.070, −0.077 for (+) and
  −0.123, −0.132, −0.136, −0.139 for (−), approaching the predicted −0.083 and −0.143 with a gap that halves at every
  doubling of k_e (the O(k_e⁻³) remainder); the predicted difference at Θ̄₂ — 7.8·10⁻⁴, 1.7·10⁻⁴, 3.9·10⁻⁵, 9.4·10⁻⁶ — is
  within 14%, 9%, 5%, 3% of the measured one. The prediction with the Jacobian taken at Θ* instead of Θ̄₂ (the
  third-derivative term omitted) reads 9.0·10⁻⁴, 2.1·10⁻⁴, 5.0·10⁻⁵, 1.2·10⁻⁵: it agrees at k_e = 10 by accident and its
  k_e²-scaled value converges to 0.078 while the measured one converges to 0.060, the 30% error that the numerics
  exposed (`review/REVIEW_LOG.md`, round 5, K-1). The separation exponents come in exact pairs (multiplicity m − 1 = 2)
  and the collective exponents differ between the orientations by 2% at k_e = 10 and 0.2% at k_e = 80, as (i)–(iii)
  describe.
- *At the critical leak* (σ₀ = σ_crit = 0.4606). The leading exponent vanishes and the k_e⁻² terms decide: the softest
  separation exponent is −4.1·10⁻⁴ for (+) and −1.27·10⁻³ for (−) at k_e = 10, −4.8·10⁻⁵ and −8.7·10⁻⁵ at k_e = 40 —
  unstable for both orientations in this instance, the even part dominating the odd one by a factor 2–3; the predictions
  at Θ̄₂ are within 14% and 6% of the measured values at k_e = 40.
- *The orientation decides* (σ₀ designed from the averaged data). Solving pred₊(σ₀) + pred₋(σ₀) = 0 for the predicted
  softest exponents gives σ₀ − σ_crit = 0.0125 at k_e = 10 and 0.0061 at k_e = 20, the centre of the window of (iv); there
  the measured softest separation exponents are +8.5·10⁻⁴ and −2.5·10⁻⁵ (k_e = 10), +1.4·10⁻⁴ and −3.7·10⁻⁵ (k_e = 20):
  the same reference orbit, traversed one way, keeps the collapsed kernels together, and traversed the other way makes
  them separate. The window is narrow — the leading exponent at these σ₀ is 1.2·10⁻³ and 3.0·10⁻⁴ — and the design lands
  inside it although the k_e⁻³ terms move both measured exponents by about +4.7·10⁻⁴ and +6·10⁻⁵ against the prediction.

*Checked* (`tests/test_periodic.py`): the Hessian splitting (i) against finite differences (10⁻⁵ relative) at the
collapsed minimiser of instance A, the isotypic eigenvalue structure, the gap equation of σ_crit (10⁻⁸), and the two-fold
zero of the separation block at σ_crit with a positive definite collective block
(`test_proposition_K_hessian_splitting_and_critical_leak`); the experiment above for (iii)–(iv).

---

## 6f. Theorem J′ — the stationary law of the dithered finite-m closed loop

This section proves the second step of §8, branch 1 — the finite-m closed loop *with dither* — and answers the first
"next step" of branch 2 (the orientation effect when the stationary state is a measure). It also opens a new angle: the
stationary state of the dithered loop is a non-equilibrium steady state, with a stationary probability current of order
1/k_e that vanishes exactly when the first-order averaged drift is reciprocal (Theorem I), and the stationary parameter
process is not reversible in time.

**Setting.** The closed loop of Theorem H's design at finite m, in the scaling σ = σ₀/k_e, τ = τ₀/k_e (the formal limit of
Remark H.0: k_e → ∞ with G, k_eσ, k_eτ fixed), ε := 1/k_e:

    ẋ = f(x) + ẋ_m − k_e e − f̂_Θ(x),    e := x − x_m(t),
    dϑ_i = [G∇_ϑ(eᵀΦ̃(x;ϑ_i)) − (σ₀/k_e)(ϑ_i − ϑ⁰)]dt + √(2τ₀/k_e) G^{1/2} dB_i,    i = 1..m,                     (CL_τ)

with the notation of §6d: saturated features, Θ ∈ ℝ^{mp}, G_m := I_m ⊗ G, x_m ∈ C^∞ T₀-periodic with occupation measure ρ,
f ∈ C^∞ on a neighbourhood of the tube 𝒯_R, F_R := sup_{𝒯_R}‖f‖, and r_Θ, J_m, H, u, ⟨H⟩ = −G_m∇J_m as there. Assume
k_e ≥ k_T := (F_R + K_Φ)/R, so that the ball {‖e‖ ≤ R} is forward invariant for every parameter path (‖ė‖ ≤ F_R + K_Φ −
k_e‖e‖ on it): the plant never sees f outside the tube, and f may be modified outside a neighbourhood of 𝒯_R to be smooth
with bounded derivatives on ℝⁿ without changing the law of the loop.

*Slow time and the three scales.* Put s := εt, θ := t mod T₀ ∈ 𝕋 := ℝ/T₀ℤ, and η := e/ε = k_e e. Then dW_s := √ε dB_t is a
Brownian motion in slow time and (CL_τ) reads

    dη/ds = ε⁻²[r_Θ(x_m(θ) + εη) − η],    dϑ_i = [G∇_ϑ(ηᵀΦ̃(x_m(θ) + εη; ϑ_i)) − σ₀(ϑ_i − ϑ⁰)]ds + √(2τ₀) G^{1/2} dW_i,    dθ/ds = ε⁻¹:

the error relaxes at rate ε⁻² to the slaved value η*(Θ,θ) := r_Θ(x_m(θ)), the phase turns at rate ε⁻¹, and the parameters
move at rate one. The generator of the Markov process Z := (Θ, η, θ) on 𝒮 := ℝ^{mp} × {‖η‖ ≤ F_R + K_Φ} × 𝕋 is

    L_ε = ε⁻²A_ε + ε⁻¹∂_θ + B_ε,    A_ε := (r_Θ(x_m(θ) + εη) − η)·∇_η,
    B_ε := Σ_i[G∇_ϑ(ηᵀΦ̃(x_m(θ) + εη; ϑ_i)) − σ₀(ϑ_i − ϑ⁰)]·∇_{ϑ_i} + τ₀Σ_i tr(G∇²_{ϑ_i}),

with A_ε = A₀ + εA₁ + …, A₀ = (η* − η)·∇_η, A₁ = (D_xr_Θ(x_m)η)·∇_η, and B_ε = B₀ + εB₁ + …, where B₀ has the drift
Σ_i[G∇_ϑ(ηᵀΦ̃(x_m;ϑ_i)) − σ₀(ϑ_i − ϑ⁰)] (equal to H(Θ,θ) at η = η*) and B₁ = Σ_iG∇_ϑ(ηᵀD_xΦ̃(x_m;ϑ_i)η)·∇_{ϑ_i}.

*The averaged objects.* L̄₀ := −G_m∇J_m·∇ + τ₀tr(G_m∇²) is the Langevin generator in the G⁻¹ metric with invariant measure
the Gibbs measure π := Z⁻¹exp(−J_m/τ₀) (Theorem H(ii) at finite m; π has Gaussian tails because J_m ≥ (σ₀/2)Σ‖ϑ_i − ϑ⁰‖²_{G⁻¹} −
C). The first-order drift is the vector field on ℝ^{mp}

    b(Θ) := b_fr(Θ) + b_osc(Θ),    b_fr,i := G∇_ϑ(U_even + U_odd)[μ_Θ](ϑ)|_{ϑ=ϑ_i},    b_osc := ⟨(D_ΘH)(Θ,θ)·u(Θ,θ)⟩_θ,

with μ_Θ := (1/m)Σ_jδ_{ϑ_j} and U_even, U_odd the functionals of Proposition H.7 (so b = k_e²(b_frozen + b_osc) of Theorem J,
now at every Θ and not only at critical points; `periodic.second_order_drift` and `dither.SingleKernelInstance.averaged`
compute it). The *first-order correction* g: ℝ^{mp} → ℝ is the solution with π(g) = 0 of the weighted Poisson equation

    τ₀∇·(πG_m∇g) = ∇·(πb)    (equivalently L̄₀g = −τ₀π⁻¹∇·(πb) = ∇J_m·b − τ₀∇·b, which has zero π-mean),

the *first-order current* is j₁ := π(b − τ₀G_m∇g) = πb_⊥ with b_⊥ := b − τ₀G_m∇g, so that ∇·j₁ = 0 and b = G_m∇(τ₀g) + b_⊥
is the Hodge decomposition of b in L²(π; G_m⁻¹) into a G_m-gradient and a π-divergence-free field. For an observable
φ ∈ C^{3,α}(ℝ^{mp}) with polynomially bounded derivatives, ψ_φ denotes the solution of L̄₀ψ_φ = φ − π(φ), π(ψ_φ) = 0.

**Theorem J′.** Assume the setting above, with G diagonal (as in every instance of this document; for a general symmetric G
replace G_l∂_{ϑ_l} by the contraction with G below). There are ε₀ > 0 and, for every observable φ ∈ C^{3,α}(ℝ^{mp}) (three
derivatives, the third Hölder continuous) with polynomially bounded derivatives up to order three, a constant C_φ
(depending on φ, f, x_m, G, σ₀, τ₀, R, m, the feature bounds, and not on ε) such that for every ε = 1/k_e ∈ (0, ε₀]:

(a) *Invariant law.* The process Z has an invariant probability measure m_ε on 𝒮, every invariant measure is supported in
{‖η‖ ≤ F_R + K_Φ} and has moments of all orders in Θ bounded uniformly in ε, and its θ-marginal is uniform. If moreover
a_max > F_R componentwise and the tube lies in the box of the saturated centres, the invariant measure is unique, so the
law of Z_t converges to m_ε and time averages along a single trajectory converge to m_ε-expectations.

(b) *The Gibbs measure at leading order.* For every invariant measure m_ε, |m_ε(φ) − π(φ)| ≤ C_φ ε.

(c) *First-order correction.* m_ε(φ) = π(φ) − ε π(b·∇ψ_φ) + O(ε²) = π(φ) + ε∫φ g dπ + O(ε²), with |O(ε²)| ≤ C_φ ε².
In particular the stationary parameter means are π(Θ) + ε∫Θ g dπ + O(ε²).

(d) *Observables of the full state.* For φ(Θ, η, θ) smooth, T₀-periodic in θ, with polynomially bounded derivatives in Θ
(φ̄ below in the class of (b)), m_ε(φ) = π(φ̄) + O(ε), where φ̄(Θ) := ⟨φ(Θ, r_Θ(x_m(θ)), θ)⟩_θ. In particular the stationary mean-square tracking error is

    E_{m_ε}‖e‖² = (2/k_e²) π(L̄_m) + O(k_e⁻³),    L̄_m(Θ) = ½∫‖r_Θ‖²dρ:

the Gibbs average of the occupation loss, divided by k_e²/2.

(e) *Circulation and irreversibility.* Let V(Θ,η,θ) := k_e·(drift of Θ in physical time) = Σ_i[G∇_ϑ(ηᵀΦ̃(x_m(θ) + εη;ϑ_i))
− σ₀(ϑ_i − ϑ⁰)] (so Θ̇ = εV + noise). For every pair of coordinates (k, l) of ℝ^{mp},

    Λ_ε^{kl} := E_{m_ε}[Θ_kV_l − Θ_lV_k] = ε[ ∫(Θ_k j₁,l − Θ_l j₁,k) dΘ + π(⟨u_k u̇_l − u_l u̇_k⟩_θ) ] + O(ε²),

where u̇ := ∂_θu = H − ⟨H⟩. If the bracket is non-zero for some (k, l), then for ε small the stationary parameter process
(Θ_t)_{t≥0} under m_ε is not reversible in time: E[Θ_{k,t}Θ_{l,t+h}] ≠ E[Θ_{l,t}Θ_{k,t+h}] for small h > 0.

*Status.* Proved. The proof imports three standard theorems, named where they are used: the existence, growth and
regularity theory of the Poisson equation for the Langevin generator with a confining potential (Pardoux–Veretennikov
[S23], with interior Schauder estimates [S24] for the higher derivatives), and, for the uniqueness in (a) only, Hörmander's
theorem for the smoothness of the transition densities of a hypoelliptic time-periodic diffusion [S25], [S26], the
Stroock–Varadhan support theorem [S27] and the Doob–Khasminskii and Harris theorems [S28], [S32], [S33]. Parts (b)–(e) hold for every invariant measure and do not use the uniqueness.
The first-order drift b is the one of Theorem J, so (c) is the thermal version of Theorem J(d): at low temperature it
reduces to the O(1/k_e) offset of the periodic orbit (Remark J′.1). The orientation effect of Remark H.7a and Proposition
K appears in (c) as the orientation-odd part of ∫Θ g dπ, and in the numerics it is measured to within 7–8% at k_e = 10 and
within the Monte Carlo resolution (10–30%) at k_e = 20–80.
Checked numerically on the single-kernel instance below (`exp_dither_stationary.py`, `tests/test_dither.py`), where the
first-order prediction removes the O(1/k_e) discrepancy between the real loop and the Gibbs measure; the circulation of (e)
is below the Monte Carlo resolution of the real loop at the affordable cost, and the current machinery behind it is
checked on the first-order averaged SDE instead (Proposition J′.2).

*Proof.* Throughout, C denotes a constant independent of ε, and "polynomially bounded" means bounded by C(1 + ‖Θ‖^q) on
𝒮 for some q.

*Step 0 (invariance, moments, existence).* By the choice k_e ≥ k_T, {‖e‖ ≤ R} is forward invariant and ‖e(t)‖ ≤ e^{−k_et}‖e(0)‖
+ (F_R + K_Φ)/k_e, so the set {‖η‖ > F_R + K_Φ + δ} is left in finite time and never re-entered: every invariant measure
is supported in {‖η‖ ≤ F_R + K_Φ}, and on this set the Θ-drift of (CL_τ) is bounded by λ_G(F_R + K_Φ)L_ϑ + σ‖Θ − Θ⁰‖.
For W := ‖Θ − Θ⁰‖² (slow time) L_εW ≤ −2σ₀W + 2λ_G(F_R + K_Φ)L_ϑ√(mW) + 2τ₀ tr G_m ≤ −σ₀W + C, and likewise L_εW^q ≤ −c_qW^q
+ C_q for every q ≥ 1 (Itô's formula and Young's inequality). Hence E_zW^q(Z_s) ≤ e^{−c_qs}W^q(z) + C_q/c_q: the Cesàro
averages of the laws of Z_s are tight (W^q is coercive in Θ, η and θ live in compact sets), the semigroup is Feller (the
coefficients are smooth with bounded derivatives in η and θ and of linear growth in Θ), and Krylov–Bogoliubov gives an
invariant probability measure m_ε. For any invariant measure, 0 = m_ε(L_εW^q) ≤ −c_q m_ε(W^q) + C_q gives the uniform
moment bounds; the θ-marginal is uniform because θ moves at constant speed. Dynkin's formula with these moment bounds
shows m_ε(L_εψ) = 0 for every ψ ∈ C²(𝒮), T₀-periodic in θ, with polynomially bounded derivatives up to order two: the
local martingale in Itô's formula is a martingale, and the stationary expectation of ψ(Z_s) − ψ(Z_0) vanishes.

*Step 1 (uniqueness).* The noise acts on all of Θ; the x-equation carries no noise. Write the SDE for (x, Θ) in physical
time with the vector fields X₀ (drift, time-periodic) and X_{i,k} := √(2τ₀/k_e)(G^{1/2})_{·k}·∂_{ϑ_i}. The bracket
[X_{i,α_k}, X₀] contains, in the x-direction, −(1/m)(∂_{α_k}Φ̃)(x;ϑ_i)·∂_{x_k} = −(1/m)(da_k/dα_k)φ(x;c_i,s_i)∂_{x_k} with
da_k/dα_k > 0 and φ > 0: together with the X_{i,k} these brackets span ℝⁿ × ℝ^{mp} at every point and every time, so the
parabolic Hörmander condition holds everywhere. By Hörmander's theorem for time-dependent hypoelliptic diffusions with
smooth coefficients of bounded derivatives ([S25]; the non-homogeneous case is [S26]), the transition probabilities
P_{s,t}(z, ·) have densities that are smooth in z, so the period map P := P_{0,T₀} is strong Feller. Accessibility: fix
Θ₁ and let z* := (e*, Θ₁) be the state at time T₀ of the deterministic skeleton that follows, on [0, T₀/2], a C¹ parameter
path along which f̂_{Θ(t)}(x(t)) = f(x(t)) − w(t) with w(t) := λ·β(t) for a fixed bump β supported in (0, T₀/2) and λ ∈ ℝ
chosen so that e(T₀/2) = 0 (possible since ė = −k_ee + w is controllable and the required λ is O(e^{−k_eT₀/2}); such
parameter paths exist because a_max > F_R: m kernels centred at x(t) with amplitude f(x(t)) − w(t) realise f̂ = f − w, and
the saturated centre box contains the tube), and on [T₀/2, T₀] any C¹ path from that configuration to Θ₁. From every
initial (e₀, Θ₀) in the invariant set the same construction (with λ depending on e₀) reaches e(T₀/2) = 0 at the same
parameter configuration, hence z* at time T₀. By the Stroock–Varadhan support theorem [S27] (the skeleton paths are
realised by controls h ∈ L²: the Θ-noise is non-degenerate and the skeleton's Θ-path is C¹), z* ∈ supp P(z, ·) for every
z. The supports of two distinct ergodic invariant measures of a strong Feller kernel are disjoint, and both contain z*
([S28, Prop. 4.1.1 and Thm 4.2.1]; [S32, Cor. 1.14–Thm 1.16]), so P has at most one invariant measure. The periodic family
P_{0,t}^*m_ε is then the unique invariant law of the augmented process; the strong Feller property with the accessible
point gives a minorisation of P² on the sublevel sets of the Lyapunov function W ([S32, Lemma 1.13]), and Harris'
theorem ([S33]) gives the convergence of the law of Z_{nT₀} to m_ε in a W-weighted total variation norm, hence the
convergence of time averages along a single trajectory.

*Step 2 (the fast operator).* At frozen (Θ, θ), A₀ is the generator of the flow η̇ = η* − η, which contracts to η* at
rate one. For g smooth with g(Θ, η*, θ) = 0 the function

    (A₀⁻¹g)(Θ, η, θ) := −∫₀^∞ g(Θ, η* + e^{−s}(η − η*), θ) ds

is smooth, satisfies A₀(A₀⁻¹g) = g, and its derivatives in (Θ, θ) of any order are bounded by those of g; if g is linear in
η − η*, g = (η − η*)ᵀw(Θ,θ), then A₀⁻¹g = −(η − η*)ᵀw. Conversely A₀ψ = 0 with ψ bounded in η forces ψ = ψ(Θ, θ).

*Step 3 (the expansion).* Let φ ∈ C⁵(ℝ^{mp}) with polynomially bounded derivatives up to order five. Define

    ψ₁ := −u·∇φ,    w_φ := Σ_{i,l}G_l ∂_{ϑ_l}Φ̃(x_m(θ);ϑ_i) ∂_{(i,l)}φ ∈ ℝⁿ    (for a symmetric G: (w_φ)_k = Σ_i(∇_ϑΦ̃_k(x_m;ϑ_i))ᵀG∇_{ϑ_i}φ),
    ψ₂ := (η − η*)ᵀw_φ + ψ₂′(Θ,θ),

where ψ₂′ is the zero-mean antiderivative in θ of b·∇φ − [−ṙᵀw_φ + (D_xr_Θ r)ᵀw_φ + (B₀ψ₁)|_{η*} + (B₁φ)|_{η*}] with
ṙ := ∂_θ[r_Θ(x_m(θ))] and r := r_Θ(x_m(θ)), and ψ₃ := A₀⁻¹ of the η-dependent remainder defined below. Then
ψ_ε := φ + εψ₁ + ε²ψ₂ + ε³ψ₃ satisfies, order by order:

- ε⁻²: A₀φ = 0.
- ε⁻¹: A₀ψ₁ + ∂_θφ + A₁φ = 0, since ψ₁ does not depend on η.
- ε⁰: A₀ψ₂ + ∂_θψ₁ + A₁ψ₁ + (B₀ + A₂)φ. Here A₁ψ₁ = A₂φ = 0 and A₀ψ₂ = −(η − η*)ᵀw_φ = (B₀φ)|_{η*} − B₀φ, so the sum
  equals ∂_θψ₁ + (B₀φ)|_{η*} = −u̇·∇φ + H·∇φ + τ₀tr(G_m∇²φ) = ⟨H⟩·∇φ + τ₀tr(G_m∇²φ) = L̄₀φ, because u̇ = H − ⟨H⟩. This is
  Remark H.0's algebraic identity: the slaved drift averages to the Langevin drift.
- ε¹: A₀ψ₃ + ∂_θψ₂ + A₁ψ₂ + (B₀ + A₂)ψ₁ + (B₁ + A₃)φ. At η = η*: ∂_θψ₂ = −ṙᵀw_φ + ∂_θψ₂′ (the lag of the error behind its
  slaved value), A₁ψ₂ = (D_xr_Θ(x_m)r)ᵀw_φ (the state shift in the error equation), B₀ψ₁ = −H·∇(u·∇φ) − τ₀tr(G_m∇²(u·∇φ)),
  A₂ψ₁ = A₃φ = 0, and B₁φ = Σ_iG∇_ϑ(rᵀD_xΦ̃(x_m;ϑ_i)r)·∇_{ϑ_i}φ (the state shift in the feature). By the choice of ψ₂′
  the sum at η = η* equals b·∇φ pointwise in θ, once the θ-averages are identified: ⟨−ṙᵀw_φ⟩ = G∇_ϑ⟨−ṙᵀΦ̃(x_m;ϑ)⟩·∇φ
  = (G_m∇U_odd)·∇φ at the particles, ⟨(D_xr r)ᵀw_φ + (B₁φ)|_{η*}⟩ = ⟨(r·∇_x)(rᵀΦ̃)⟩-gradient·∇φ = (G_m∇U_even)·∇φ (the two
  state shifts combine into U_even), so these three give b_fr·∇φ; and ⟨B₀ψ₁|_{η*}⟩ = −⟨(D_Θu)H⟩·∇φ − ⟨H⊗u⟩:∇²φ −
  τ₀tr(G_m∇²(⟨u⟩·∇φ)) = ⟨(D_ΘH)u⟩·∇φ = b_osc·∇φ, because ⟨u⟩ = 0, ⟨(D_Θu)u̇⟩ = −⟨(D_Θu̇)u⟩ = −⟨(D_ΘH)u⟩ (integration by
  parts in θ, ⟨D_Θu⟩⟨H⟩ = 0), and ⟨H_ju_k⟩ + ⟨H_ku_j⟩ = ⟨∂_θ(u_ju_k)⟩ = 0 so the symmetric contraction with ∇²φ vanishes
  (the same cancellation as N₂ = 0 in Lemma E.0′). Writing T₁(Θ,η,θ) for the sum of the order-ε terms other than A₀ψ₃, one
  has T₁(Θ,η*,θ) = b·∇φ, so ψ₃ := A₀⁻¹(b·∇φ − T₁) is well defined (Step 2) and A₀ψ₃ + T₁ = b·∇φ: the order-ε term is
  exactly L̄₁φ := b·∇φ.
- ε²: R_ε := ε⁻²[L_εψ_ε − L̄₀φ − εL̄₁φ] is a finite sum of products of derivatives of f on the tube (orders ≤ 3), of Φ̃ in
  (x, ϑ) (bounded), of x_m, of u, w_φ, ψ₂′, ψ₃ and of derivatives of φ: ψ₂′ contains ∇³φ (through τ₀tr(G_m∇²(u·∇φ)) in
  (B₀ψ₁)|_{η*}) and ψ₃ contains ∇²φ, so B₀ψ₂ ∋ τ₀tr(G_m∇²ψ₂′) carries ∇⁵φ and B₀ψ₃ carries ∇⁴φ; using the exact operators
  A_ε, B_ε and Taylor's theorem in ε rather than their formal series, |R_ε| ≤ C(1 + ‖Θ‖^q) on 𝒮 uniformly in ε ∈ (0, ε₀]
  for φ ∈ C⁵ with polynomially bounded derivatives up to order five. (For the zeroth-order statement (b) alone ψ₂′ can be
  dropped — it only serves to make the order-ε term exactly b·∇φ — and C³ suffices.)

*Step 4 (the Poisson equation).* L̄₀ is uniformly elliptic with smooth coefficients, ∇J_m has linear growth with bounded
derivatives of all orders, and ⟨−G_m∇J_m(Θ), Θ − Θ⁰⟩ ≤ −σ₀‖Θ − Θ⁰‖² + C‖Θ − Θ⁰‖: the recurrence condition of
Pardoux–Veretennikov [S23, Thm 1] holds with every polynomial Lyapunov function. For φ ∈ C^{3,α}(ℝ^{mp}) with polynomially
bounded derivatives up to order three there is a unique ψ_φ ∈ C²(ℝ^{mp}) with L̄₀ψ_φ = φ − π(φ), π(ψ_φ) = 0, and ψ_φ, ∇ψ_φ are
polynomially bounded ([S23, Thm 1 and 2]); interior Schauder estimates on unit balls [S24, Thm 6.2 and 6.17], which gain two
derivatives and whose constants depend polynomially on ‖Θ‖ through the local Hölder norms of the coefficients, give
ψ_φ ∈ C^{5,α} with polynomially bounded derivatives up to order five.

*Step 5 (b) and (c).* Apply Step 3 to ψ := ψ_φ and Step 0's identity to ψ_ε: 0 = m_ε(L_εψ_ε) = m_ε(L̄₀ψ) + εm_ε(b·∇ψ) +
ε²m_ε(R_ε) = m_ε(φ) − π(φ) + εm_ε(b·∇ψ) + O(ε²), the O(ε²) by the moment bounds of Step 0. Since b·∇ψ is polynomially
bounded, |m_ε(φ) − π(φ)| ≤ Cε, which is (b). Now b·∇ψ ∈ C⁴(ℝ^{mp}) with polynomially bounded derivatives (b is smooth with
bounded derivatives in Θ, ψ ∈ C^{5,α}), so (b) applies to it (its Poisson solution is in C³, which the zeroth-order
expansion needs): m_ε(b·∇ψ) = π(b·∇ψ) + O(ε), and
m_ε(φ) = π(φ) − επ(b·∇ψ_φ) + O(ε²). The second form of (c) is Green's identity for the weighted Poisson equation:
π(b·∇ψ) = ∫b·∇ψ π = −∫ψ∇·(πb) = −τ₀∫ψ∇·(πG_m∇g) = −τ₀∫g∇·(πG_m∇ψ) = −∫gπL̄₀ψ = −∫g(φ − π(φ))π = −∫φg dπ, using
τ₀π⁻¹∇·(πG_m∇ψ) = L̄₀ψ and π(g) = 0; the boundary terms vanish by the Gaussian tails of π and the polynomial growth of
g, ψ and their gradients (the same regularity for g as for ψ_φ, with φ replaced by the polynomially bounded ∇J_m·b − τ₀∇·b).

*Step 6 (d).* Write φ = φ̄ + φ_θ + φ_η with φ_θ(Θ,θ) := φ(Θ,η*,θ) − φ̄(Θ) (zero θ-mean) and φ_η := φ(Θ,η,θ) − φ(Θ,η*,θ)
(zero at η*). With χ₁ := −(zero-mean antiderivative in θ of φ_θ) and χ₂ := −A₀⁻¹φ_η (Step 2), L_ε(εχ₁ + ε²χ₂) = ∂_θχ₁ +
A₀χ₂ + ε[B₀χ₁ + (∂_θ + A₁)χ₂] + O(ε²) = −(φ_θ + φ_η) + ε[…] + O(ε²), all terms polynomially bounded. Step 0's identity gives
m_ε(φ) = m_ε(φ̄) + εm_ε(B₀χ₁ + (∂_θ + A₁)χ₂) + O(ε²) = π(φ̄) + O(ε) by (b). For φ = ‖η‖², φ̄ = ⟨‖r_Θ(x_m(θ))‖²⟩_θ = 2L̄_m(Θ),
and E‖e‖² = ε²E‖η‖².

*Step 7 (e).* Expand V = V⁰ + εV¹ + O(ε²) with V⁰ the B₀-drift and V¹ the B₁-drift, and apply Step 6 to φ := Θ_kV_l − Θ_lV_k
keeping the order-ε terms. For φ⁰ := Θ∧V⁰ (∧ on the pair (k,l)): φ̄⁰ = Θ∧⟨H⟩ = Θ∧(−G_m∇J_m), and π(φ̄⁰) = 0 because
−πG_m∇J_m = τ₀G_m∇π and ∫Θ∧(G_m∇π) = 0 for a symmetric G_m (integration by parts: the integrand is a divergence plus
π·(G_{lk} − G_{kl})). Next, φ⁰_θ = Θ∧(H − ⟨H⟩) = Θ∧u̇, so χ₁ = −Θ∧u and ⟨(B₀χ₁)|_{η*}⟩ = −⟨H·∇_Θ(Θ∧u)⟩ − τ₀⟨tr G_m∇²(Θ∧u)⟩
= −⟨H∧u⟩ − ⟨Θ∧(D_Θu·H)⟩ = −⟨u̇∧u⟩ + Θ∧⟨(D_ΘH)u⟩ = ⟨u∧u̇⟩ + Θ∧b_osc (the τ₀ term and ⟨H⟩∧⟨u⟩ vanish since ⟨u⟩ = 0, and
⟨(D_Θu)H⟩ = −⟨(D_ΘH)u⟩ as in Step 3). Further, φ⁰_η = Θ∧(V⁰(η) − V⁰(η*)) is linear in η − η*, so χ₂ = φ⁰_η and
⟨((∂_θ + A₁)χ₂)|_{η*}⟩ = Θ∧⟨Σ_iG∇_ϑ((−ṙ + D_xr r)ᵀΦ̃(x_m;ϑ_i))⟩ = Θ∧(lag and error-shift parts of b_fr); with φ̄¹ = Θ∧⟨V¹|_{η*}⟩
= Θ∧(feature-shift part of b_fr) these complete Θ∧b_fr. Collecting, and using (c) for φ̄⁰ (m_ε(φ̄⁰) = ε∫Θ∧(−G_m∇J_m)g dπ
+ O(ε²)) and (b) for the order-ε observables,

    Λ_ε = ε[ ∫Θ∧(b − gG_m∇J_m) dπ + π(⟨u∧u̇⟩) ] + O(ε²),

and b − gG_m∇J_m = (b − τ₀G_m∇g) + (τ₀G_m∇g − gG_m∇J_m) = b_⊥ + τ₀π⁻¹G_m∇(gπ), whose second term integrates to zero
against Θ∧ as above. This is the formula of (e). For the irreversibility: under m_ε the Θ-process is stationary;
reversibility would give E[Θ_{k,t}Θ_{l,t+h}] = E[Θ_{l,t}Θ_{k,t+h}] for all h, i.e. E[Θ_t∧Θ_{t+h}] = 0, whereas
E[Θ_t∧(Θ_{t+h} − Θ_t)] = E[Θ_t∧∫_t^{t+h}εV ds] + 0 (the martingale term has mean zero) = hεΛ_ε + o(h) by stationarity and
the continuity of the integrand; so Λ_ε ≠ 0 excludes reversibility. ∎

*Remark J′.1 (reading).*

- *From the periodic orbit to the thermal state.* As τ₀ → 0 with Θ* the non-degenerate global minimiser of J_m, π → δ_{Θ*},
  ψ_{Θ_k} → the linear function v_kᵀ(Θ − Θ*) with −(G_m∇²J_m)ᵀv_k = e_k, and (c) gives m_ε(Θ) → Θ* + ε(G_m∇²J_m(Θ*))⁻¹b(Θ*):
  Theorem J(d)'s period mean of the periodic orbit. Theorem J′(c) is the finite-temperature version: the correction is a
  Gibbs average of the same first-order drift, weighted by the Poisson solution instead of the inverse Hessian.
- *The orientation effect with dither (branch 2).* Reversing the direction of traversal of the reference flips U_odd and
  b_osc and leaves U_even, J_m and π unchanged, so the stationary means of the two traversals differ at order ε by
  2ε∫Θ g_odd dπ with g_odd the correction driven by the odd part of b: the quantity measured in Remark H.7a for the
  periodic orbit is now a statement about stationary measures, and it is checked below (to within 7–8% at k_e = 10, and
  within the resolution of the runs at larger k_e) in the single-kernel instance (where, in one state dimension, U_odd ≡ 0 —
  a 1-form integrated over a path that retraces itself
  — and the whole odd effect comes from b_osc; a reference with the time-reversal symmetry x_m(θ₀ + θ) = x_m(θ₀ − θ), such
  as a sinusoid, has no odd effect at all, since its reversal is a time shift). As a design rule: a
  time-symmetric reference (one whose Fourier series carries only cosines after a time shift, such as x_m = X[cos θ + r₂cos 2θ])
  removes b_osc and U_odd identically, so its first-order correction is the performative term G_m∇U_even alone, which is
  reciprocal exactly when Theorem I's A_even ≡ 0; a reference traversed one way carries both, and the odd part can be used
  (Proposition K) or cancelled by reversing the traversal. Checked: for the cosine reference the code's b_osc and ⟨u∧u̇⟩ are
  10⁻¹⁷ at every point tested.
- *Non-equilibrium.* (e) says that the dithered learning loop, at stationarity, is not in detailed balance even to first
  order: the periodic forcing alone (the term π⟨u∧u̇⟩, the signed area swept by the fast parameter loop) breaks the
  time-reversal symmetry, and the slow dynamics carries in addition the current j₁ of Proposition J′.2, nonzero exactly
  when the first-order drift is not reciprocal. The first is a property of any periodically forced system; the second is
  the signature of closed-loop learning (Theorem I).
- *Hypotheses.* Smoothness of f on the tube and of the features is used through the expansion (three x-derivatives of f,
  five Θ-derivatives of the Poisson solution, i.e. three Hölder derivatives of the observable) and, for the uniqueness only,
  through Hörmander's theorem; k_e ≥ k_T keeps the
  plant in the tube for every parameter path. The constants are not explicit: the result is an asymptotic expansion of the
  invariant measure, not a quantitative bound at a given k_e, and the numerics show where it applies.

**Proposition J′.2 (the first-order current, entropy production, and reciprocity).** Consider the first-order averaged
dynamics dΘ = (−G_m∇J_m + εb)ds + √(2τ₀)G_m^{1/2}dW on ℝ^{mp}, with b as above.

(i) Its stationary law m^{(1)}_ε satisfies m^{(1)}_ε(φ) = π(φ) + ε∫φg dπ + O(ε²) for every φ as in Theorem J′ (the
expansion m^{(1)}_ε = π(1 + εg) + O(ε²) in the weak sense; on the grid it holds in L¹ with an ε² remainder, Numerics
below), and its stationary current j := m^{(1)}_ε(−G_m∇J_m + εb) − τ₀G_m∇m^{(1)}_ε equals εj₁ + O(ε²) in the same sense, with
j₁ = πb_⊥; its entropy production rate (in slow time) is

    EPR_ε := ∫‖j‖²_{G_m⁻¹}/(τ₀m^{(1)}_ε) dΘ = ε²τ₀⁻¹π(‖b_⊥‖²_{G_m⁻¹}) + O(ε³)

(formal in the remainder: the ε² term uses only the first-order current, the O(ε³) presumes the expansion of j in L²(1/π)).

(ii) The following are equivalent: j₁ ≡ 0; b_⊥ ≡ 0; b is a G_m-gradient field, b = −G_m∇J₁; D_Θ(G_m⁻¹b) is a symmetric
matrix at every Θ. In that case the first-order averaged dynamics *is* the Langevin dynamics of J_m + εJ₁, with stationary
density exactly ∝ exp(−(J_m + εJ₁)/τ₀), in detailed balance.

(iii) The antisymmetric part of D_Θ(G_m⁻¹b) is, block by block ((i,j) blocks of size p × p),

    ½[D(G_m⁻¹b) − D(G_m⁻¹b)ᵀ]_{ij} = (1/(2m))∇_ϑ∇_{ϑ'}A_{μ_Θ}(ϑ_i, ϑ_j) + ½⟨[∇²_Θ𝒥, ∇²_Θ𝒰]_{G_m}⟩_{θ,ij},

where A_μ = S − Sᵀ is the antisymmetric kernel of Theorem I at the atomic measure μ_Θ (derivatives in the kernel arguments,
μ_Θ frozen), 𝒥(Θ,θ) := (m/2)‖r_Θ(x_m(θ))‖² + (σ₀/2)Σ_i‖ϑ_i − ϑ⁰‖²_{G⁻¹} is the instantaneous loss (so H = −G_m∇𝒥, ⟨𝒥⟩ = J_m),
𝒰 its zero-mean antiderivative in θ, and [X, Y]_G := XGY − YGX. The first term is the finite-m form of the reciprocity
defect of Theorem I; the second is the failure of the instantaneous Hessians at different phases to commute, the
non-gradient part of second-order averaging of a time-periodic gradient flow.

(iv) In the setting of Theorem I (saturated Gaussian features, f ∈ C², x_m ∈ C² periodic), for m ≥ 3, b_fr is never a
G_m-gradient on ℝ^{mp} (b_fr does not depend on σ₀ and depends on G through the factor G only); consequently, along the gain
scale t ↦ tG,
the full drift b = b_fr + b_osc is a (tG)_m-gradient for at most one value of t > 0. For m = 1, 2 the condition of (ii) is
the explicit vanishing of (iii), checked to fail in the single-kernel instance below (non-gradient fraction
π(‖b_⊥‖²)/π(‖b‖²) ≈ 0.3–0.5).

*Proof.* (i) Stationary Fokker–Planck: ∇·(mG_m∇J_m) − ε∇·(mb) + τ₀∇·(G_m∇m) = 0; at order ε⁰ the Gibbs measure, at order
ε¹ with π₁ = gπ: ∇·(gπG_m∇J_m) + τ₀∇·(G_m∇(gπ)) = τ₀∇·(πG_m∇g) (using τ₀∇π = −π∇J_m), so the equation is τ₀∇·(πG_m∇g) =
∇·(πb), the defining equation of g; the expansion is justified by the same Poisson-equation argument as in Theorem J′
(the generator of this SDE is L̄₀ + εb·∇, and Steps 4–5 apply verbatim without the fast variables). The current:
j = π(−G_m∇J_m) − τ₀G_m∇π + ε[πb − gπG_m∇J_m − τ₀G_m∇(gπ)] + O(ε²) = ε π(b − τ₀G_m∇g) + O(ε²). The entropy production rate
of a diffusion with mobility G_m and temperature τ₀ in a stationary state m with current j is ∫jᵀ(τ₀G_m)⁻¹j/m [S29],
which gives the formula.
(ii) j₁ ≡ 0 iff b = τ₀G_m∇g iff b is a G_m-gradient (then τ₀g = −J₁ + const and the equation for g is satisfied); on the
simply connected ℝ^{mp} a smooth field G_m⁻¹b is a gradient iff its Jacobian is symmetric (Poincaré lemma). When b = −G_m∇J₁ the
drift −G_m∇(J_m + εJ₁) is a gradient drift, so the SDE is the Langevin dynamics of J_m + εJ₁ with the exact Gibbs law
∝ exp(−(J_m + εJ₁)/τ₀) and j ≡ 0; consistently, g = −J₁/τ₀ + const solves the defining equation.
(iii) G⁻¹b_fr,i = ∇_ϑP₂[μ_Θ](ϑ)|_{ϑ_i} with P₂ := U_even + U_odd, a functional of μ (quadratic in μ through r_μ) and of the
evaluation point. Differentiating with respect to ϑ_j moves the atom ϑ_j of μ_Θ, which gives (1/m)∇_ϑ∇_{ϑ'}S(ϑ_i, ϑ_j) with
S(ϑ,ϑ') := δP₂[μ](ϑ)/δμ(ϑ') (Theorem I's k_e⁻² second-variation kernel), and, for j = i, additionally the symmetric
∇²_ϑP₂(ϑ_i). The antisymmetric part of the (i,j) block is therefore (1/(2m))∇∇'[S(ϑ_i,ϑ_j) − S(ϑ_j,ϑ_i)] = (1/(2m))∇∇'A_μ(ϑ_i,ϑ_j).
For b_osc: H = −G_m∇𝒥 and u = −G_m∇𝒰, so G_m⁻¹b_osc = ⟨∇²𝒥G_m∇𝒰⟩ and ∂_j(G_m⁻¹b_osc)_k = ⟨∂_j∂_k∂_l𝒥G_l∂_l𝒰⟩ +
⟨(∇²𝒥G_m∇²𝒰)_{kj}⟩; the first is symmetric in (j,k), the antisymmetric part of the second is ½⟨∇²𝒥G_m∇²𝒰 − ∇²𝒰G_m∇²𝒥⟩.
(iv) Suppose b_fr is a G_m-gradient for some m ≥ 3. By (iii), ∇_ϑ∇_{ϑ'}A_{μ_Θ}(ϑ_i,ϑ_j) = 0 for all Θ and all i, j. A_μ is
affine in μ (Theorem I): A_μ = A_f + B(∫Φ̃dμ) with B linear in the function ∫Φ̃(·;ϑ)dμ(ϑ), so for Θ = (ϑ, ϑ', ϑ₃, …),
∇_ϑ∇_{ϑ'}[A_f + (1/m)Σ_jB(Φ̃(·;ϑ_j))](ϑ,ϑ') = 0 with the atom positions inside B frozen. Varying ϑ₃ alone shows that
∇_ϑ∇_{ϑ'}B(Φ̃(·;ϑ₃))(ϑ,ϑ') does not depend on ϑ₃; at the zero feature ϑ₃ = ϑ_z it vanishes (B(0) = 0), so it vanishes for all
ϑ₃, and then ∇∇'A_f(ϑ,ϑ') = 0 as well: ∇_ϑ∇_{ϑ'}A_μ(ϑ,ϑ') = 0 for every μ and all (ϑ,ϑ'). A smooth antisymmetric kernel
with vanishing mixed derivative on the connected ℝ^p × ℝ^p has the form a(ϑ) − a(ϑ'), a coboundary, hence A_μ ≡ 0 by
Theorem I(i) — for every μ, contradicting Theorem I(vi). For the gain scale: under G → tG, b_fr → t b_fr and b_osc → t²b_osc,
so the antisymmetric part of D((tG)⁻¹b) is a_fr + t a_osc with a_fr ≢ 0; it vanishes identically in Θ for at most one t. ∎

**Corollary J′.3 (the stationary law in the iterated limit k_e → ∞, m → ∞).** Let π_m be the Gibbs measure of Theorem J′ on
ℝ^{mp}, and μ*_τ the unique minimiser of Theorem H's free energy F_τ at σ = σ₀/k_e, τ = τ₀/k_e; since k_eF_τ = L̄ + (σ₀/2)∫‖ϑ −
ϑ⁰‖²_{G⁻¹}dμ + τ₀H(μ) = L̄ + τ₀H(·|γ_τ) + const with γ_τ = N(ϑ⁰, (τ₀/σ₀)G), μ*_τ does not depend on k_e. Then:

(i) *Chaoticity of the Gibbs measures.* Under π_m the empirical measure μ_Θ = (1/m)Σ_iδ_{ϑ_i} converges in probability to μ*_τ
as m → ∞, and for every k the k-particle marginal of π_m converges weakly to (μ*_τ)^{⊗k}.

(ii) *The iterated limit of the real loop.* For every φ ∈ C²(ℝ^p) bounded with bounded derivatives up to order two,

    lim_{m→∞} lim_{k_e→∞} E_{m_ε}[(1/m)Σ_i φ(ϑ_i)] = ∫φ dμ*_τ,

and the inner limit is attained at the rate O(1/k_e) of Theorem J′(b), with a constant that depends on m.

*Status.* Proved (Sanov's theorem and Varadhan's lemma, [S34, Thms 6.2.10 and 4.3.1]; Theorem H(i)–(ii); Theorem J′(b)).
The joint limit, with k_e and m related, is Conjecture H.6 and remains open; the corollary places the stationary law of the
finite loop between Theorem J′ and Theorem H. Numerics: `exp_gibbs_m.py` → `results/dither/gibbs_m.md` (below).

*Proof.* J_m(Θ) = m[L̄(μ_Θ) + (σ₀/2)∫‖ϑ − ϑ⁰‖²_{G⁻¹}dμ_Θ], so π_m(dΘ) ∝ exp(−mL̄(μ_Θ)/τ₀) γ_τ^{⊗m}(dΘ). The map μ ↦ L̄(μ) =
½∫‖f − ∫Φ̃(·;ϑ)dμ(ϑ)‖²dρ is bounded and continuous for the weak topology of P(ℝ^p) (the features are bounded and continuous
in ϑ, and ρ is a probability measure). Under γ_τ^{⊗m} the empirical measure satisfies Sanov's large-deviation principle with
rate function H(·|γ_τ); by Varadhan's lemma (tilting by the bounded continuous −mL̄/τ₀), under π_m it satisfies a
large-deviation principle with rate function I(μ) = L̄(μ)/τ₀ + H(μ|γ_τ) − inf = (k_e/τ₀)(F_τ(μ) − F_τ(μ*_τ)), by Theorem H(i).
F_τ is strictly convex (L̄ is convex, the entropy strictly convex) with the unique minimiser μ*_τ (Theorem H(ii)), so I has
the unique zero μ*_τ and μ_Θ → μ*_τ in probability. Exchangeability of π_m then gives the convergence of the k-marginals to
(μ*_τ)^{⊗k} (Sznitman's equivalence of chaoticity and the convergence of the empirical measure). For (ii), Theorem J′(b)
applied to the observable (1/m)Σ_iφ(ϑ_i) ∈ C²(ℝ^{mp}) gives lim_{k_e→∞}E_{m_ε}[(1/m)Σφ(ϑ_i)] = π_m((1/m)Σφ(ϑ_i)) = π_m^{(1)}(φ),
the one-particle marginal, which converges to ∫φ dμ*_τ by (i). ∎

*Remark J′.4 (what the corollary says about Conjecture H.6).* The conjecture asks for the limit of the real loop's law when
k_e → ∞ with m fixed or growing, uniformly in time. Theorem J′ settles the stationary part at fixed m, and the corollary
shows that the resulting Gibbs measures are chaotic with limit μ*_τ, the stationary state of (𝔄): the stationary law of the
real dithered loop is the minimiser of Theorem H's free energy in the iterated limit. What is missing for the conjecture is
the uniformity in m of Theorem J′'s O(1/k_e) (its constant grows with m through the Poisson equation on ℝ^{mp}) and the
time-dependent statement.

*Numerics* (`exp_dither_stationary.py` → `results/dither/table.md`). The single-kernel instance: n = 1, m = 1, plant
f(x) = e^{−(x−0.5)²/0.5} + 0.4x, reference x_m(t) = 1.5[sin t′ + ⅓ sin 2t′] with t′ = ±t (the two orientations), saturated
amplitude and centre ϑ = (α, γ) (a_max = 3, c_max = 2.5), frozen width s = 0.6, G = diag(2, 2), σ₀ = 0.4, ϑ⁰ = 0, τ₀ = 0.1,
k_e ∈ {10, 20, 40, 80}. The averaged data on a 201² grid: π, g, j₁ (rms b_fr = 0.143, rms b_osc = 0.225; non-gradient fraction
π(‖b_⊥‖²)/π(‖b‖²) = 0.32 and 0.46 for the two orientations; EPR coefficient τ₀⁻¹π(‖b_⊥‖²_{G⁻¹}) = 0.144 and 0.118; circulation
coefficient ∫Θ∧j₁ + π⟨u∧u̇⟩ = +0.0117 − 0.0167 = −0.0051 and −0.0103 + 0.0167 = +0.0065). The real loop: 2·10⁴ independent
replicas per (k_e, orientation), started from the first-order law π(1 + εg), Strang splitting with k_eh = ½, burn-in ≥ 10 and
averaging ≥ 40 slow-time units (whole numbers of reference periods), standard errors from the spread of the per-replica time
averages (the slowest mode of L̄₀ relaxes in 3.5 slow-time units on this instance, so batch means in time would be
correlated; referee round 6).

*The first-order law (b), (c).* The stationary means against the Gibbs value and the first-order prediction, orientation +1
(ε-coefficients ∫φ g dπ: −0.255 for α, −0.426 for γ):

| k_e | m_ε(α) − π(α) (s.e. 0.0013) | ×k_e | m_ε(α) − first order | m_ε(γ) − π(γ) (s.e. 0.0018) | ×k_e | m_ε(γ) − first order |
|---|---|---|---|---|---|---|
| 10 | −0.0276 | −0.28 | −0.0021 | −0.0441 | −0.44 | −0.0015 |
| 20 | −0.0179 | −0.36 | −0.0051 | −0.0269 | −0.54 | −0.0056 |
| 40 | −0.0089 | −0.36 | −0.0025 | −0.0127 | −0.51 | −0.0020 |
| 80 | −0.0046 | −0.37 | −0.0015 | −0.0083 | −0.67 | −0.0030 |

The discrepancy between the real loop and the Gibbs measure falls like k_e^{−0.87} (α), k_e^{−0.83} (γ), k_e^{−0.72} (α²) and
k_e^{−1.09} (αγ) (γ² has a first-order coefficient 0.03, below resolution), and the first-order prediction removes most of it:
93% and 97% of the α and γ discrepancies at k_e = 10, 70–85% at k_e = 20–40; at k_e = 80 the residuals (−0.0015 ± 0.0013
and −0.0030 ± 0.0018) are at the Monte Carlo floor. The second-order decay of the residual is *not* resolved: the residuals at
k_e = 20–80 are 1–3 standard errors with a common sign, compatible with an O(ε²) coefficient of order −1 to −2 and with noise;
at k_e = 10 the even part of the α² and γ² shifts exceeds its prediction by a factor 1.5–1.8 (an O(ε²) term of the size seen
in Theorem J's instance A, where c₂/c₁ ≈ 1.4).

*The orientation-odd part (Remark J′.1).* Measured (m₊ − m₋)/2 over predicted (p₁₊ − p₁₋)/(2k_e):

| k_e | α | γ | αγ |
|---|---|---|---|
| 10 | 0.92 ± 0.04 | 0.93 ± 0.04 | 0.89 ± 0.06 |
| 20 | 1.23 ± 0.09 | 1.10 ± 0.07 | 1.28 ± 0.13 |
| 40 | 1.28 ± 0.17 | 1.15 ± 0.15 | 1.03 ± 0.26 |
| 80 | 1.25 ± 0.34 | 1.38 ± 0.29 | 1.11 ± 0.51 |

The dither version of Remark H.7a's test: the direction in which the reference is traversed shifts the stationary kernel by
2ε∫Θ g_odd dπ (here entirely from b_osc), measured to within 7–8% at k_e = 10 (where the odd shift is 0.02–0.03 against standard
errors of 0.001) and consistent with the prediction within the 10–30% resolution of the runs at k_e = 20–80 (the ratios sit
10–30% above one with error bars of the same size: a second-order odd term of relative size O(ε) is not excluded).

*The tracking error (d).* k_e²E‖e‖²/(2π(L̄)) = 1.054, 1.033, 1.017, 1.009 (orientation +1) and 1.023, 1.014, 1.008, 1.004 (−1)
at k_e = 10, 20, 40, 80, with standard errors 0.001: the deviation from one is 0.54–0.74/k_e and 0.23–0.35/k_e, the O(1/k_e)
of (d).

*The circulation (e) and Proposition J′.2.* The predicted first-order circulation ε[∫Θ∧j₁ + π⟨u∧u̇⟩] is −5·10⁻⁴ (+1) and
+6.5·10⁻⁴ (−1) at k_e = 10 and scales like 1/k_e, against a Monte Carlo standard error of 5·10⁻⁴ at every k_e: the real loop's
circulation is not resolved at first order for k_e ≥ 20 (measured −4·10⁻⁴ … +1.2·10⁻³, all within 1.7 standard errors of the
prediction), and at k_e = 10 the measurement (+3.7·10⁻³ ± 0.5·10⁻³ for +1, −0.6·10⁻³ ± 0.5·10⁻³ for −1) is dominated by
higher orders (an even O(ε²) term of coefficient ≈ +0.2 and an odd one of the opposite sign to the first-order odd part). The
current machinery is therefore checked where the signal is clean, on the first-order averaged SDE itself: its stationary
density m_ε and current j_ε, computed on the grid (the discrete current at ε = 0, −6.6·10⁻⁵, a residual of the central scheme,
subtracted), give ∫Θ∧(j_ε − j₀)/(ε∫Θ∧j₁) = 1.47, 1.21, 1.09, 1.03 at ε = 0.2, 0.1, 0.05, 0.025 — one plus a term linear in ε —
and ‖m_ε − m₀(1 + εg)‖₁/ε² = 0.307, 0.308, 0.307, 0.305: the first-order expansion of Proposition J′.2(i) with a clean ε²
remainder.

*Corollary J′.3* (`exp_gibbs_m.py` → `results/dither/gibbs_m.md`). On the same instance, μ*_τ is computed on the 201² grid as
the fixed point of the first-variation condition (59 damped iterations, residual 10⁻¹¹): ∫α dμ* = 0.4963, ∫γ dμ* = 0.1873,
∫α² dμ* = 0.9353, L̄(μ*_τ) = 0.1394. The Gibbs measures π_m for m = 1, 2, 4, 8, 16, 32 are sampled by MALA (4000 replicas,
acceptance 0.98–0.996, standard errors from the independent replicas; the m = 1 row agrees with the quadrature of π₁ within
one standard error):

| m | E_{π_m}[α₁²] | E_{π_m}[γ₁] | E_{π_m}[L̄(μ_Θ)] | (E[α₁²] − ∫α²dμ*)·m | (E[γ₁] − ∫γdμ*)·m | (E[L̄] − L̄(μ*))·m |
|---|---|---|---|---|---|---|
| 1 | 0.536 | 0.251 | 0.1960 | −0.40 | +0.064 | +0.057 |
| 2 | 0.729 | 0.225 | 0.1681 | −0.41 | +0.075 | +0.057 |
| 4 | 0.832 | 0.204 | 0.1540 | −0.41 | +0.067 | +0.058 |
| 8 | 0.881 | 0.194 | 0.1470 | −0.43 | +0.050 | +0.061 |
| 16 | 0.909 | 0.191 | 0.1432 | −0.42 | +0.051 | +0.060 |
| 32 | 0.922 | 0.190 | 0.1413 | −0.42 | +0.084 | +0.059 |

The one-particle marginals converge to μ*_τ with a clean 1/m law (slopes −1.01 for γ and −0.99 for the Gibbs-averaged
occupation loss, whose rescaled deviation is 0.058 ± 0.002 at every m; the second moment of α moves from 0.54 at m = 1 to 0.92
at m = 32 against 0.935, with rescaled deviation −0.41 ± 0.01 throughout: the single kernel sits in a narrow well, while in
the mean field each particle sees the wide one-body potential V[μ*_τ]); the mean of α is within 0.01 of its limit already at
m = 1 (a small 1/m coefficient, its slope not resolved). The stationary law of the finite dithered loop thus approaches the
minimiser of Theorem H's free energy in the iterated limit at the rate 1/m, as the corollary states.

*Checked* (`tests/test_dither.py`): the single-kernel averaged data (J, ∇J, b_fr, b_osc) against the general implementation
of Proposition H.7 and Remark H.7a (`averaging.expansion_terms`, `averaging.oscillation_correction`, `periodic.ClosedLoop`)
at three parameter points (2·10⁻⁵ relative); the weighted Poisson solver against the Ornstein–Uhlenbeck closed forms (g =
G⁻¹b·Θ/τ₀, mean shift (Gλ)⁻¹b, the Green identity ∫φg dπ = −π(b·∇ψ_φ), zero current and entropy production for a gradient
drift); the stationary Fokker–Planck solver of the first-order SDE (OU exact; on the instance the current ratio → 1 and
‖m − π(1 + εg)‖₁ ∝ ε²); the real loop at k_e = 20 (8000 replicas: the stationary means of α and γ closer to the first-order
prediction than to π, the γ shift resolved at ≥ 3 standard errors, the tracking error within 10% of 2π(L̄)); and the
Corollary J′.3 numerics (the mean-field fixed point to 10⁻⁹, the Gibbs sampler against the quadrature of π₁)
(`test_single_kernel_averaged_data_matches_the_general_H7_implementation`,
`test_poisson_solver_reproduces_the_ornstein_uhlenbeck_closed_forms`,
`test_first_order_stationary_density_and_current_of_the_averaged_sde`,
`test_real_loop_means_are_closer_to_the_first_order_prediction_than_to_gibbs`,
`test_corollary_J3_mean_field_minimiser_and_gibbs_sampler_consistency`).

---

## 7. Verification map

Test files are abbreviated: `T` is `tests/test_transport.py`, `Lg` is `tests/test_langevin.py`, `Av` is
`tests/test_averaging.py`, `Ex` is `tests/test_excitation.py`, `NR` is `tests/test_nonreciprocity.py`, `Pe` is
`tests/test_periodic.py` and `Ba` is `tests/test_barycentric.py`.

| statement | what is checked | where |
|---|---|---|
| constants (1), closed forms | ∇²_{(c,s)}φ ≤ H_φ by brute force | `T::test_gaussian_constants_are_suprema` |
| Lemma P.1 | ‖Proj_K(ϑ,Gy)‖_{G⁻¹} ≤ ‖y‖_G on 2000 adversarial states, with equality when inactive; a tight instance of F(i) with ratio ≥ 0.99 | `T::test_lemma_P1_projection_is_nonexpansive_and_tight_velocity_bound` |
| A(i) raw rates | lifted raw rates equal classical raw rates under the gain map (1e-10) | `T::test_identification_with_classical_law` |
| A(i) projected rates, only-if | per-kernel projection with rescaled layer matches to 1e-10 with active projections; an unrescaled layer does not; five joint balls are each distinguished (two-case construction: four active against an inactive per-kernel field, one inactive against an active one) | `T::test_identification_projected_rates_per_kernel` |
| B with constant ideal and ε = 0 | the identity to fp precision; cross term ≤ 0; remainder ≤ Lemma B.2; V̇_π ≤ (B) at every sample | `T::test_theorem_B_inequality_along_trajectories` |
| B with a moving ideal, ε ≠ 0, switches and active projection | the same plus V_*W; at least one switch; projection active at ≥ 10 samples; Dini quotient ≤ (B) | `T::test_theorem_B_moving_ideal_nonrepresentable_switches_and_active_projection` |
| B, Dini form | forward-difference quotient of V_OT ≤ (B) | `T::test_dini_derivative_of_V_OT_bounded_by_active_branch` |
| B against labelling | matched error ≪ labelled error; the bound is 30–45× smaller | `T::test_permuted_ideal_labelled_vs_OT`, `exp_transport_permutation.py` |
| B.5: failure of the plain bound along (L_σ); the (L_σ) variant holds | with every C(iv) hypothesis asserted, the plain bound fails by > 3× and the leak variant holds | `T::test_corollary_B5_fails_along_lsigma_and_its_leak_variant_holds` |
| C(i), Lemma C.3 (unprojected pure leak) | zero switches and exact e^{−σt} decay (explicit Euler, pure leak) | `T::test_theorem_C_i_exact_dissipation_pure_leak` |
| C(ii) | C ≤ b² + 4ηR; √C ≤ b + 2η; C = b² under both gap conditions (300 configurations) | `T::test_lemma_C_ii_induced_coupling_bound` |
| C(iv) | an instance satisfying **every** hypothesis, asserted, including version 3's gap condition 4η̄(b + η̄) < Δ_gap; the envelope holds at every sample, with actual contraction | `T::test_theorem_C_iv_exponential_estimate` |
| C(v), Prop. C.5 | label invariance of the OT leak; the reshuffling transient of the labelled leak | `exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py` |
| C′(i) | the plan reduces to the assignment for m₀ = m (1e-10) and W² = V_bar² + Var_γ (1e-10); along the pure leak W² is non-increasing at every step and D⁺Z ≤ −σ_mV_bar² holds at every state; the flow reaches V_bar² < 1e-8 and ends at a fixed point of the equal-mass Lloyd iteration | `Ba::test_plan_reduces_to_assignment_and_pythagorean_split`, `Ba::test_pure_barycentric_leak_dissipates_to_a_centroidal_configuration`, `exp_barycentric_leak.py` → `results/barycentric_leak/table.md` |
| C′(ii) | √C_bar ≤ b + 2η on 90 random configurations; C_bar = b_bar² at small η (plan stability) | `Ba::test_induced_coupling_bound_and_plan_stability` |
| C′(iii) | the Lyapunov inequality at every recorded state of a closed loop with the barycentric leak toward a 24-atom prior | `Ba::test_lyapunov_inequality_along_the_barycentric_leak_closed_loop` |
| D | consistency in m (nearly constant amplification 1.64–1.69; not the constant e^{ΛT}) | `exp_mean_field_m.py` |
| D′(b), (d) | the m^{−1/2} rate of the synchronous coupling to a 4096-particle reference: slopes −0.47 ± 0.06 (τ = 0) and −0.46 ± 0.05 (τ = 2·10⁻³) over 24 seeds, 0.94 ± 0.08 and 0.92 ± 0.08 against (1/m − 1/M)^{1/2}; not the constant | `exp_chaos.py` → `results/chaos/table.md` |
| D″ | with exact Lipschitz constants and the a-priori tube, σ_cert = 4.04 against Proposition K's σ_K = 0.00144; at 1.2σ_cert with dither the coupled distance is of order m^{−1/2} (slope −0.51 ± 0.02), bounded in time and 700× below the bound, with the network at its prior (fit error 1.000); without dither the reference collapses (spread 10⁻¹⁴) and the distance is zero; at σ_cert/8 it decays to 10⁻¹¹; in the learning regime (σ_K/3, 1000 s) the network stays spread and fits the plant (fit error 0.63) and the distance follows m^{−1/2} (−0.62 ± 0.16) but grows threefold between the windows (100, 200] and (500, 1000] | `exp_chaos_uniform.py` → `results/chaos/uniform_table.md` |
| D‴, D‴.5, D‴.7 | the window realised at σ₀ = 0.15σ_crit (κ_R[μ*](N_{0.1}) = 0.0051 < 0.0071 < 0.047), (H) along the warm run (margin ≥ 0.0014), saturation of the coupled distance (last/first window ratios 1.13–1.15 over 900 slow units), slope −1.28 ± 0.30 in m; the weight-sensitivity law per seed within 25% in 9 of 12 pairs, per-cluster cosines 0.92–1.00; the exact stationary distance of Proposition D‴.5(ii) within 2.3% in all 12 pairs, splitting margins positive in 11 (one m = 16 seed a saddle); cold start: growth (2.9–7.8×) and the floor's fraction 0.20, 0.08, 0.04 decaying like m^{−0.64±0.08}; real loop at k_e = 25 within 0.7% of its averaged slow dynamics on the same clouds; not the constants C₁, m₀ (vacuous: 2·10⁹, 10²³) nor the saturation of Corollary D‴.7 | `exp_chaos_learning.py` → `results/chaos/learning_table.md` |
| D‴.9, D‴.11 | the matching prior: zero residual, κ_R = 0, margin σ₀ = 0.00707 against the σ-modification's 0.00196, floors 1.61× and 2.18× higher; misplacement closes the window (centres +0.6, amplitudes ×1.2, widths ×0.7); shrunk amplitudes ×0.5–0.7 give margins above σ₀ at the uniform weights and lower floors; admissible-weight table under three priors (overshoot splits, undershoot stabilises; the shrunk prior keeps all four tested weights admissible, +0.0054 to +0.0073); matching-prior coupled runs saturate within 2.1 of the theorem's relaxation times (observed e-folding 48–84 slow units against 107–306 for the σ-modification) with per-seed linearised-law ratios 0.81–1.11 (parameter) and 0.61–1.69 (function) and the exact law within 0.987–0.998 in all nine pairs; the two saddles among them are left by a local split with the coupled distance on its plateau; the function-level law for the σ-modification within ratios 0.48–0.74 of the measured means (RMS over seeds 0.52–0.84; per-seed law from the realised weights 0.60–1.57) | `exp_prior_design.py` → `results/chaos/prior_design.md` |
| D‴.12 | the same clouds under the nearest-atom and the OT leak: cluster weights conserved (max deviation 0 in 32 systems) against equalised (s_eq = 93–164 slow units; ln 2/σ₀ = 98; 1–21 particles moved); commensurate sizes: function-level distance² 4.9·10⁻⁷–7.5·10⁻⁹ against 2.5·10⁻³–5.9·10⁻⁵ and W² 7.9·10⁻³–6.9·10⁻⁴ against 5.1–0.81 (the exact K × K transport between the weighted stationary measures within 0.997–0.999 at admissible weights; the formula ½‖w^m − w^M‖₁Δ_ν² within −7% to +28%), both OT distances still falling (0.02–0.035 and 0.12 per 150 slow units) while the synchronous D rises to 0.95–0.96 of its mismatch plateau; incommensurate sizes: one or two bridge particles at 0.28–0.58 Δ_ν from the atoms; against ν: m·W² = 19.4–21.0 (slope −0.98 ± 0.03) and m·fit = 0.45–0.60 (squared slope −1.81 ± 0.05); against the M = 512 reference slopes −2.07 ± 0.05 and −1.21 ± 0.04 (the reference carries its own bridge particles; the m = 16 value still relaxing); fit error of the reference 0.012 → 8·10⁻⁵; the designed masses inadmissible between δ = 0.03 and 0.05 | `exp_ot_leak.py` → `results/chaos/ot_leak.md`; `tests/test_ot_leak.py` (the vertex plan, the gradient flow, the m-independent stationary configuration, conservation against constraint, the forest bound, the exact W²) |
| D‴.14 | the real loop at k_e = 25 from the clouds at the end of `exp_ot_leak.py` (contracted where the weights are admissible) and from the collapsed lifts: under the OT leak ‖f̂_m − f̂_M‖²_{L²(orbit)} = 7.5·10⁻⁹–2.7·10⁻⁸ (m = 15, 60 against 510) and 5.2·10⁻¹⁰–3.6·10⁻⁸ against the 3-atom lift, rms tracking error 0.9·10⁻⁶–7.3·10⁻⁶ (the lift 8·10⁻¹²); under the nearest-atom leak 1.1·10⁻⁴–2.5·10⁻³ with tracking errors 3.2·10⁻⁴–1.9·10⁻³; at the one pair with both weights admissible the disagreement equals the lifted orbits' distance to 1.0003 and the admissible systems sit within 2.5·10⁻¹¹ of their lifts; three systems at saddle weights (margins −0.0035, −0.0011, −0.0022) have split (spread 0.04–0.11) and sit off their lifts | `exp_ot_real_loop.py` → `results/chaos/ot_real_loop.md`; `tests/test_ot_leak.py` (the collapsed subspace's invariance, the lifts' equality) |
| D‴.15 | lifts of the 3-atom orbit to m = 6, 9, 15 (OT leak, matching masses) and to m = 16 (σ-modification at (7, 5, 4)/16 and (9, 4, 3)/16): fixed points of the Poincaré map to 10⁻¹⁴–10⁻¹⁵; class separation subspaces invariant to ≤ 6·10⁻⁸; separation exponents equal to the ghosts' to 2·10⁻⁵ relative and to (σ₀ − λ_j(R^G_k))/k_e within 0.2–4.1%; collective exponents equal to the K₀-atom orbit's to 10⁻⁴ relative for every m and to Theorem J(b)'s prediction within 0.6–3.2%; the saddle's exponent −9.13·10⁻⁵ (multiplicity 8) at (9, 4, 3)/16 | `exp_lifted_orbit.py` → `results/chaos/lifted_orbit.md`; `tests/test_ot_leak.py` (the m = 6 lift against the ghosts and the 3-atom loop) |
| D‴.16 | the dithered averaged slow loop, m = 15 from the counts (6, 5, 4): nearest-atom record-level label changes 1.0·10⁻⁵, 5.1·10⁻⁴, 5.6·10⁻³ per slow unit at τ₀ = 0.005, 0.01, 0.02 (Arrhenius barrier 0.048, among the ridge minima 0.040–0.074); time at the balanced counts 0.00, 0.54, 0.25 against the energy–multinomial law 0.86, 0.58, 0.33 (standard errors ≈ 0.1 at τ₀ = 0.01) and against the OT leak 0.996, 0.992, 0.88, whose excursions are returns of the crossing particle (121 of 121, 298 of 336) with the effective barrier 0.058 | `exp_dither_ot.py` → `results/chaos/dither_ot.md` |
| D‴.14 (n = 2) | a planar plant (four vector-valued bumps, a Lissajous reference, p = 5), k_e = 25, σ₀ = 0.05, sizes 8, 32 against 128 from the same clouds: the OT leak equalises the weights at 13–14 slow units (ln 2/σ₀ = 14), the sizes agree to 5.0·10⁻⁸–1.6·10⁻⁷ with rms tracking errors 3.9·10⁻⁶–1.6·10⁻⁵, against 5.5·10⁻³–2.3·10⁻² and 1.9·10⁻³–6.7·10⁻³ under the nearest-atom leak (fit errors along the orbit 2–9·10⁻⁴ against 0.11–0.39) | `exp_ot_2d.py` → `results/chaos/ot_2d.md` |
| D‴.17 | J_m^{OT} non-increasing along 300 steps of a run with reassigned particles (strong leak; the plan does not switch there, so this is a smooth gradient flow), dJ/ds = −‖drift‖²_{G⁻¹} where the plan is unique, and at a constructed tie between two optimal plans (mirror-image particles) J^{OT} equals both frozen energies and decreases along either plan's drift | `tests/test_ot_leak.py` |
| D‴.18 | the real loop at k_e = 25, 8000 time units, sizes 3–240: the 3-atom orbit's Floquet exponents within 0.56–3.22 % of λ_j(H_q)/k_e and the ghosts' equal to (σ₀ − λ_j(R^G_k))/k_e to 3·10⁻⁴ (matching and shifted priors); inside the cells the per-particle distance to the lift decays at 2.69–2.97·10⁻⁴ for m ≥ 15 (prediction 2.83·10⁻⁴) with constants 1.03–1.29 under both leaks, the m = 3 runs mixtures of collective modes; across a cell boundary (ρ = 4.5) the nearest-atom leak ends at other counts (W² 1.9–12.3, that of the lifts of the realised weights; two of them inadmissible) and the OT leak at q (W² 0.010–0.23, falling at 2σ₀); one particle across a boundary: W² grows to 0.79, 0.20 (nearest) and falls to 3.0·10⁻³, 8.2·10⁻⁴ (OT); on a moving orbit (shifted prior) 2.18–2.63·10⁻⁴ against 2.27·10⁻⁴, constants 1.04–1.23; the proof's radius at most 0.019 | `exp_uniform_basin.py` → `results/chaos/uniform_basin.md`; `tests/test_uniform_basin.py` (at the shifted prior: the remainder constant and the exact N_a + N_b/m structure of the per-particle term, the separation response (GR_k − σ₀I)v, the contraction at μ = min{λ_coll, λ_sep} with one constant, the same for every m) |
| D‴.20 | the critical excess by continuation in ε (Newton residuals ≤ 8.7·10⁻¹³): ε_c = 5.11, 13.5, 23.1, 65.0 for δ = 0.02–0.2 (the exact threshold equal to ε_c; the rank-one formula 1.22–1.31 and 0.70–1.22 times ε_c with the two atoms), 6.99, 14.6, 24.3 for an excess taken from one atom; δ_c(ε) = 0.0042–0.20 for ε = 1–64; dynamics for m = 15–510 at ε = 1, 4, 16, 64 over 600 slow units, continued to 6000: equalised at ε = 1; at ε = 4 stable end weights with residual excesses up to 1/60, 3/240, 7/510; at ε = 16 and 64 the stable systems keep their sampled weights and the four unstable ones move (one moving its excess to an atom that was at q); floors 0.55 and 0.74 of the nearest-atom law, 0.23–0.84 and 0.47–0.91 of the observed nearest-atom floors; fit errors 0.12–0.14 at ε = 64 | `exp_entropic_leak.py` → `results/chaos/entropic_leak.md`; `tests/test_entropic_leak.py` (the solvers, the gradient flow, the limits at their rates, the separation block entrywise with the invariance of the subspace, the critical regularisation by continuation in ε with A ≻ 0 checked and every residual below 10⁻⁹) |
| D‴.21 | the averaged loop from the warm clouds (m = 15–510, four seeds, matching and shifted priors, both leaks, 600 slow units): (E), (P) to 6·10⁻¹³; under the OT leak no relabelling, α > 0 from s ≤ 300, ratio/λ ≥ 1.00 (matching) and ≥ 1.61 (shifted), E from 4.1·10⁻³–4.0·10⁻² to 8.3·10⁻⁷–8.1·10⁻⁶ below the a-posteriori bound, the end rate 0.01416–0.01417 against 2σ₀ = 0.01415 and λ = 0.0140–0.0141, at the shifted prior the curvature's cost 0.802 and (ii)'s slack 0.683; the nearest-atom E at the spurious minima's energy at the 22 admissible count vectors (within 0.24 %), split and still falling 0.5–16.5 % below it at the 10 inadmissible ones (m ≤ 60); (H*) on the whole space at the shifted prior (κ = 4.7·10⁻⁴, E₀ = 9.1·10⁻¹⁰), the a-priori basins holding no particle in another cell below m ≈ 1.7·10⁵ (matching) and 2.3·10⁶ (shifted); the margin +0.0023 with 1/30 of the particles moved into another cell, negative from 1/15 | `exp_ot_basin.py` → `results/chaos/ot_basin.md`; `tests/test_ot_basin.py` (the identities for every assignment with the counts and their exact defects at other counts, the PL chain and its sharp constant, the decay from a particle in another cell below the a-posteriori bound with a positive margin, the admissible nearest-atom spurious minimum) |
| D‴.23 | the linear theory at the exact representation (h_j = 3.06·10⁻⁴–0.195; λ_1^− = 2.68·10⁻⁵–3.06·10⁻⁴ for η = σ₀/10 … ∞; (v)'s radius at most 2.0·10⁻⁵–2.2·10⁻⁴, 94–6·10⁴ times below the K₀-atom runs' starts, 2.2·10⁴–3.6·10⁵ below the warm clouds'); the fixed prior's bias within 4.8% of the first-order formula in norm (0.85 and 0.73 of it in modes 1 and 3, the error first order in the shift); the K₀-atom loop: every mode above 10⁻⁶ decays at λ_j^− (log-slopes 0.981–1.003 of it), fit error 9.78·10⁻³ → 1.8·10⁻⁵ (η = σ₀) and 6.7·10⁻⁶ (η = ∞) by s = 6000; separated starts (five draws): departures from the K₀-atom run second order in the separations, consistent with the same expected size at m = 15 and 510; m = 15–510 from the warm clouds: ungated the prior is dragged by the transit (‖ȳ‖ 0.032–0.71 at s = 500) and one m = 15 system at η = ∞ is carried onto a plateau of the unregularised loss (Hessian eigenvalues −1.8·10⁻⁵, 7.2·10⁻⁵; still there at s = 23000, crossed by the gradient flow, below 10⁻¹² after 7.97·10⁴ more slow units); the transit gate leaves a residue in the weak modes (‖ȳ‖ 4.9·10⁻³–3.7·10⁻² against 1.3·10⁻³, η = σ₀); the settled gate (t_open = 1211–1304) opens replicated and follows the K₀-atom run for every m (per-seed ratios 0.995–1.000 at η = σ₀, 0.968–1.001 at η = ∞), fit errors 1.6·10⁻⁴ and 3.4–3.5·10⁻⁵ against 9.78·10⁻³; the limits of (i): a non-replicated saddle held by the leak, and unbounded adaptive solutions, beyond the saturation and from an exact six-kernel fit (the width coordinate running off toward a Karush–Kuhn–Tucker point of the box-constrained loss) | `exp_adaptive_prior.py` → `results/chaos/adaptive_prior.md`; `tests/test_adaptive_prior.py` (the gradient structure, stationarity, the linearisation's rates and its self-adjointness in the gradient-flow metric, the first-order bias, the two gates, the limits of (i)) |
| D‴.24 | the floor as the bias: F equal to the bias form to the finite-difference error (1.1·10⁻⁶–1.2·10⁻⁸); the law (B) 1.36·10⁻², 6.17·10⁻², 0.229 per unit of 1/m − 1/M at σ₀/σ_crit = 0.05, 0.15, 0.5 (0.50–0.84 of the bound, which is σ₀/4 times the adaptive prior's coupled-distance law 46.3), modes 2–6 carrying 98–99% of it and 97% unweighted, the weight moving the shares toward h ≈ σ₀; the parameter split to 1.4·10⁻⁷, the adaptive prior's coupled-distance law 1.4–5.3 times the fixed prior's; 5 of the 12 sampled m-systems are fixed-prior saddles (class 1, w₁ ≥ 0.4) and split (floors 0.81–1.01 of the replicated point's), the exact bias at the sampled weights reproduces the floor at the other seven (ratio 1.0000), the delta method misses it by 0.86–1.72; the nearest-atom leak over 6000 slow units: the fixed prior's floor stays, the adaptive prior's decreases without a floor (at most 2.0·10⁻⁶ at η = σ₀ and 5.4·10⁻⁸ at η = ∞ over the representable samples at s = 6000); D 0.94–1.00 of the replicated stationary points' under both priors, the adaptive prior's 1.72–4.81 times the fixed prior's; W² the cross-class transport (0.67–0.99 of it), the two priors' W² within 1.00–1.32; the boundary sample (8, 4, 3)/15 decays like 1/s (6000·‖f̂_m − f̂_M‖² = 0.88 and 0.42 against 0.8 and 0.4); the design table at the misplaced prior (the fit error and the floor of the four designs; the OT designs' D the largest, 6.0–1.0, the transit's mismatch; the ranking at s = 1250 and 2000); the real loop: the fixed prior at the averaged floor (ratios 0.996–1.001; the m = 15 run on a saddle orbit), the adaptive prior along the averaged decay (0.987–0.992), the exact representation invariant (five-period rms ≤ 9.0·10⁻¹²); the instance's probabilities (outside W 0.361, 0.074, 0.0021; Hoeffding 1.87, 0.60, 0.036; the fixed prior's saddles 0.369, 0.260, 0.109) and the adaptive prior's expected limiting floor (≥ 7.3·10⁻³, 2.2·10⁻⁴, 6.6·10⁻⁷: 1.8, 0.24, 0.005 times the fixed prior's law); Monte Carlo (600 nested samples): the exact floor 2.78, 1.53, 1.11 times the law; in two dimensions (four bumps, p = 5, σ₀ = 0.05, six seeds): every sample D‴.5-admissible and the fixed prior's floor the exact bias at the sampled weights (ratio 1.0000 in all 18 systems), the adaptive prior's at most 0.019 (η = σ₀) and 0.006 (η = ∞, the atoms reset every slow unit) of it at the 7 representable systems and still falling, the real loop (η = ∞) along the averaged prediction (0.979–1.009) and invariant from the exact representation (≤ 4.7·10⁻¹⁰) | `exp_adaptive_floor.py` → `results/chaos/adaptive_floor.md`; `tests/test_adaptive_floor.py` (the nearest-atom adaptive flow as a gradient flow, the free-mass energy as the minimum of the OT energy over the masses, also against masses moved by ±1/m, the floor law as the bias derivative with the weight's band-pass against the unweighted shares, the exact representation of a D‴.5-admissible sample stationary under the adaptive prior and the margins σ₀ at a fixed-prior saddle, the parameter split, the energy across a cell crossing and an empty cell) |
| D‴.26 | the three-bump instance, w = (6, 5, 4)/15 (representable; a fixed-prior saddle), the real loop with the adaptive prior in the state, η = σ₀ and ∞, k_e = 25 and 50, the monodromy by central differences (step 10⁻⁵, DOP853 at rtol 10⁻¹³): return residuals ≤ 1.7·10⁻¹⁵; at m = 15 and 30 the 36 and 81 separation exponents at σ₀/k_e to ≤ 1.7·10⁻⁶ (relative), the other 18 (η = σ₀) or 9 (η = ∞) equal to the weighted K₀-atom loop's to 1.4·10⁻⁵–6.4·10⁻⁴; the collapsed variational equation (its generator the real loop's linearisation to 5.8·10⁻¹¹, its exponents the differences' to 4.1·10⁻⁵) at k_e = 12.5–400 for x_m = 2 sin t and 2 sin t + 0.8 sin 2t (‖Ā₁‖ = 0 and 5.6·10⁻³, skew): the first-order deviations falling by → 4 per doubling (the slowest mode 13.6% and 3.84% below at k_e = 25, 50; every mode below on the instance), the second-order residuals by → 16 (4.9·10⁻⁷ at k_e = 400), c₁/λ₁ = −100.6 (the filter's part −102.1 = −0.989ω₁², ω₁ = 10.16), measured −100.5 at k_e = 400; the gain scaled by 30: c₈/λ₈ = +4.17, c₉/λ₉ = +12.65, those exponents above their prediction by c_j/(λ_jk_e²) to three digits; at η = 10σ₀ Ā₀ self-adjoint in diag(I, (σ₀/η)I) to 8.7·10⁻¹⁸ (6.4·10⁻² in I); the error channel by Liouville's formula within O(k_e⁻³) of k_e − tr H_L(w)/k_e | `exp_adaptive_floquet.py` (Part C: `--only-C`) → `results/chaos/adaptive_floquet.md`; `tests/test_adaptive_floquet.py` (at w = q, η = ∞, m = 6: the return of one period, the nine separation exponents, the K₀-atom loop's spectrum within the averaging correction; at (2, 1, 1)/4, η = σ₀, m = 4, the prior in the state: the real loop's Jacobian splitting into the separations and the collapsed variational equation, the monodromy's 3 separation and 18 collective exponents, the second-order prediction; on an asymmetric reference: Ā₁ skew, ℓ_jĀ₁r_j = 0, the residual O(k_e⁻⁴), c₁ measured, Liouville's formula; at η = 10σ₀ Ā₀ self-adjoint in the metric and not in I, Ā₁ skew in it) |
| D‴.26a | Corollary D‴.26's instance on x_m = 2 sin t + 0.8 sin 2t, η = 0 (K = Ĵ, Λ = σ₀I), σ₀ and ∞, k_e = 12.5–1600, a 512-point Fourier grid: the contraction for P in 3–9 iterations (499 at η = 0, k_e = 12.5: stalled at rounding), invariance residual ≤ 3.8·10⁻¹⁵, |P − P₂|/ε³ → 3.34·10³, |P − P₍₃₎|/ε⁴ → 5.8·10⁴; the slow block's rates the full variational equation's to ≤ 1.0·10⁻⁸; the rates at ±ε equal to ≤ 5.9·10⁻¹² of the largest, ≤ 1.1·10⁻⁹ per mode (first-order deviation 9.8·10⁻²–8.9·10⁻⁶ of the largest, per mode up to 29.6%: the slowest mode −29.6%, −12.5%, −4.0%, −1.1% at k_e = 12.5–100, η = ∞), with a mismatch 1.1·10⁻²–7.8·10⁻⁵ (per mode 6.1·10⁻²–3.5·10⁻⁴); the Magnus parameter 0.48–0.004; the noise: decoupling residual ≤ 3·10⁻¹⁶, ‖b − b₀‖/‖b₀‖ 0.63–0.0074, |b − b₀ − εb₁|/ε² → 3.34·10³, Cov(v)·2ε³ − 1 from −9.6·10⁻² to −2.4·10⁻⁵, ‖E[vY′ᵀ] + bᵀ‖/‖b‖ 0.15–1.2·10⁻³, ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ 0.35–9.4·10⁻⁵ (η = ∞; O(ε²) to k_e = 200–400, O(ε) beyond), the weakest mode's deviation −34% → −6.6·10⁻⁵ the sum of the noise's attenuation (−56% → −1.1·10⁻⁴) and the rate's correction (+49% → +4.6·10⁻⁵), the jitter with its next term at O(ε²) (0.30 → 0.001 at k_e = 100–1600), the floor at t = 0 at 0.763–1.0000 of the prediction, its period average at 0.766–0.99997 of D‴.30(ii)'s formula and pointwise in [0.9952, 1.0091] at k_e = 1600, the jitter's predicted range | `exp_adaptive_floquet.py` (`--add-D`) → `results/chaos/adaptive_floquet.md`; `tests/test_adaptive_floquet.py` (Chang's P and the exact slow multipliers, the ε³ truncation on the instance; on a synthetic block with an O(1), non-symmetric Λ, P − P₂ at ε³ and P − P₍₃₎ at ε⁴, a flipped Λ and a dropped cubic term caught; the rates even in ε at η = σ₀ and ∞, a mismatch breaking it; the slow noise decoupled, b₀ against K̇ by the chain rule and the plant's own slope by differences, b₁ at ε², the covariances' orders) |
| D‴.27 | the box 𝒵 (amplitudes and centres to 0.9 of a_max, c_max; the width coordinate in [−6, 6]): beyond the saturation (Part G(b)'s target) the projected prior converges at η = σ₀ (the projected midpoint rule, positions projected, dt = 0.25 and 0.125, to s = 10⁵; end states within 3.4·10⁻¹², the trajectory at s = 300 within 6.3·10⁻⁷, 1.9·10⁻⁷, 3.8·10⁻⁸ of a reference for dt = 0.25, 0.125, 0.0625) and η = ∞ (LSODA) to the same bounded stationary point (drift ≤ 3.5·10⁻¹⁵), the prior's amplitude atom on the face, fit error 0.1415; the particles' largest distance from 𝒵 0.92 and 0.89 against B/σ₀ ≤ 3686; the amplitude bound 0.5–0.99 of a_max: fit errors 0.1949–0.1101, between the fixed prior's 0.1750 and the unprojected 0.1050 for 0.7–0.99, above the fixed prior's for 0.5–0.65, the exact atoms outside the boxes up to 0.7; a Lloyd-consistent, non-optimal plan at η = ∞ along which the reduced energy rises from 1.251·10⁻³ to 1.574·10⁻³ while the selected energy falls (largest increment −3.4·10⁻⁶); on Part G(a)'s four-bump target the width coordinate stopped at the face κ = 6 from s = 10⁹ (fit error within 9.3·10⁻⁷ of the face minimiser's); the misplaced-prior run identical with and without 𝒵, the box never approached (margins ≥ 16.9% of its width) | `exp_adaptive_prior.py` Part H → `results/chaos/adaptive_prior.md`; `tests/test_adaptive_prior.py` (the projected prior bounded and convergent at η = ∞, while the unprojected coordinate passes 11 by s = 10⁴; at η = σ₀ and dt = 0.125 the prior in its box with its amplitude atom on the face, the particles bounded and 𝒥 non-increasing; inactive near the exact representation) |
| D‴.28 | the three-bump instance plus ε·exp(−(x − 0.8)²/(2·0.3²)), ε = 0–0.5: the approximant's fit error 3.6·10⁻³–2.5·10⁻² (ε = 0.05–0.5), the gauge check (the gradient at the gauge transform ≤ 1.6·10⁻¹⁵; Newton continued in w from ϑ^b(q) at it to ≤ 1.7·10⁻¹¹, f̂ to ≤ 5.5·10⁻¹³; plain Newton there in 3–5 of 6 cases), the adaptive prior's sensitivity 1.4·10⁻¹¹–2.7·10⁻¹⁰ of the fixed prior's, P(w^m ∈ W_g) = 0.52–0.64, 0.87–0.88, 0.997 at m = 15, 60, 240 (0.26, 0.71, 0.98 at ε = 0.5), the branch at the saturation near ε = 0.68, the fixed prior's within 2·10⁻⁹ of σ₀Ĵ(Ã + σ₀I)⁻¹g̃ and 7.0·10⁻⁴–5.6·10⁻³ of the Gauss–Newton form, the margins σ₀ − (q_k/w_k)λ⁺_k exact to 10⁻¹⁴ over all count vectors of m = 15, 60 (none negative: λ⁺ ≤ 0.04σ₀ at ε ≤ 0.3, ≤ 0.10σ₀ at ε = 0.5), the dynamics' separation eigenvalues equal to them to ≤ 5.8·10⁻¹²; the dynamics (m = 15, 60, 240 against 510, four seeds, ε = 0.1 and 0.3, 6000 slow units): the fixed prior's floor within 8% of its value at s = 600 (3.1·10⁻⁴–5.9·10⁻³ at m = 15, ε = 0.1) and its fit error above the approximant's, the adaptive priors' floor falling (m = 60, ε = 0.1: to 1.0·10⁻⁹–1.4·10⁻⁸ at η = σ₀, 3.2·10⁻¹¹–3.8·10⁻⁹ at η = ∞) and the fit error to the approximant's (median excess 1.2·10⁻⁴, 1.2·10⁻⁴, 1.6·10⁻⁵ at m = 15, 60, 240 for the samples 0.02 inside W_g); at the sample outside W_g the box minimum on the branch 1.107 and 1.586 times the approximant's fit error, a permuted gauge transform existing | `exp_gauge_floor.py` (`--only-A` recomputes Part A, the branch's end and the outside samples) → `results/chaos/gauge_floor.md`; `tests/test_gauge_floor.py` (the gauge's critical points and function, the continuation from ϑ^b(q), the two sensitivities, the residual curvature's structure, the margins and the dynamics' separation eigenvalues, the adaptive prior's convergence to the gauge transform at η = ∞ and σ₀) |
| D‴.29 | the servo ẍ = −τ_cog(x) + u, τ_cog = 0.4 sin(3x + 0.3) + 0.1 sin 7x, x_m = sin t, λ = 5, K₀ = 4 kernels on the position, σ₀ = 0.01: ϑ^b(q) (the local minimum BFGS reaches from a fixed start) fit error 1.46·10⁻², Ã spanning 2.9·10⁻⁴σ₀–30σ₀ (four modes below 0.1σ₀), λ⁺_k ≤ 0.0098σ₀, W_g = {w > (0.061, 0.085, 0.119, 0.008)}, the floor law 1.29·10⁻²; 120 random starts: 1 at ϑ^b(q), 33 below it, all saturated, the best 6.39·10⁻³, with its saturated amplitudes held at |α| = 6–14 7.14·10⁻³–6.40·10⁻³; the misplaced prior's (the engineer's model, fit error 0.200) replicated point a saddle (margin −3.97·10⁻³, fit error 0.0545); the averaged loop (m = 16, 64, 128 against 256, three seeds, 2000 slow units) from clouds around ϑ^b(q): fixed nearest atom floor 1.85·10⁻³, 2.24·10⁻⁴, 4.49·10⁻⁵ and fit error 0.174–0.038, adaptive (η = ∞) 4.85·10⁻⁵, 2.88·10⁻⁷, 6.36·10⁻⁸ and 0.0298–0.0147, OT 0.01462 with no floor, σ-modification 0.328, no leak 0.0180, 0.0147, 0.0146 (floors 7.7·10⁻⁶–6.6·10⁻⁸); misplaced: the fixed designs split and end at 0.048–0.10, the adaptive prior unaffected at η = ∞, 0.035–0.019 at η = σ₀; from the engineer's model: fixed 0.102, 0.060, 0.048, OT 0.048, adaptive 0.027–0.022 (η = ∞) and 0.032–0.028 (η = σ₀), no leak 0.0152–0.0147, the K₀-atom flow at η = ∞ 0.0297, 0.0213, 0.0169, 0.0155 at s = 600, 2000, 10⁴, 10⁵; phase-plane kernels: the circle 5 exact zeros and 8 below 10⁻⁶σ₀, ellipses ×1.1 and ×2 8 and 6, a two-harmonic orbit 3, centred kernels 16 exact of 20; the real servo (k_s = 25, 1000 time units, records and targets every 0.02): position-error rms 1.7–2.1·10⁻⁵ (adaptive, matching OT, no leak), 7.8·10⁻⁵ (misplaced OT), 9.1·10⁻⁵–1.8·10⁻⁴ (fixed), 5.4–5.5·10⁻⁴ (σ-modification), at the exactly filtered prediction to 1.0001–1.0029 at t ≥ 500 (the slaved 0.931–0.996) | `exp_servo_cogging.py` (`--only-A` recomputes Parts A and D, `--add-B` adds the missing Part B jobs, `--only-C [designs]` reruns Part C) → `results/servo/servo_cogging.md`; `tests/test_sliding.py` (the phase-plane loop as the gradient flow of L̄, the zero Gauss–Newton eigenvalue per kernel on a circular orbit in general position, the servo following the averaged loop with the position error at both predictions) |
| D‴.30 | the three-bump instance at the lift of the exact representation of q = 1/3, k_e = 25: h_j = 3.06·10⁻⁴–0.195, tr Q = 8.784 = 8.105 (⟨Ĵ̇ᵀĴ̇⟩) + 0.679 (⟨a²ĴᵀĴ⟩), tr H² = 0.0657, ⟨|Ĵ|⁴⟩ = 0.299; the parameters' floor per unit N₀ 0.15967, 0.16271, 0.16897, 0.17567 at η = 0, σ₀, 10σ₀, ∞, the weakest mode's variance 0.0865–2.09 (×24.1); the second-order force (OU, B = 100) and its averaged stationary biases 491 and 2397 per (Eν²)², equal to the floor at std 0.0026 and 0.0012, the full second-order mean's 2397.68, the exact periodic mean's 504 and 2653; the bias's local exponents in k_e 1.93–1.99 (B ≥ 20k_e), 1.69, 1.47 (B = 8k_e, 4k_e); the separations' drift at the lift λ_max(S_k) = 0.0148, 0.0101, 0.298 per Eν², the threshold std 0.031 there; along the second-order transient the mean error's rms 0.121 → 0.602 and λ_max(S₃) 0.307 → 0.349 over 300 time units; the exact periodic covariance of the linearised loop (k_e = 25, …, 800; η = 0, σ₀, ∞) at the prediction to 0.9762–1.00001, the control (noise input −k_eν) numerically even in 1/k_e (successive ratios → 4.00), the weakest mode's variance ×28.0 at k_e = 25, the exact stationary mean 10.7%–0.17% above the averaged (adaptive, k_e = 25–200); the real loop (RK4 dt = 10⁻³, OU B = 100, 12 common noise paths of 300 time units, std 0.003 and 0.01, η = 0 and ∞): z-scores −0.99 to 1.75 against the transient covariance plus the mean part, the mean part 42% and 89% at t = 300, the priors 6.0–6.4% apart (predicted 5.3–6.1%); the lift m = 15 on a common noise path identical to rounding (5.6·10⁻¹⁶ in float64); the separations of the lift's third class on eight antithetic pairs at −2.50·10⁻⁴, +2.74·10⁻⁶, +4.53·10⁻⁴ per time unit at std 0.01, 0.03, 0.05 (λ_eff 0.3271, 0.3175, 0.2943; the lift's drift 0.2981, the transient's 0.3264, 0.3270, 0.3281), the fit's small-noise limit 0.3291 ± 0.0033, the rate vanishing at std 0.0299 (0.0308 by the lift's drift, 0.0294 by the transient's) | `exp_measurement_noise.py` (`--only-A`, `--only-B`, `--only-D` recompute a part) → `results/chaos/measurement_noise.md`; `tests/test_measurement_noise.py` (the closed form against the Lyapunov solution, the trace identity at η = ∞ and the monotonicity in η; the exact periodic floors at k_e = 100, η = 0, σ₀, ∞; a lift m = 6 staying replicated under a noise path and following the K₀-atom loop; a separation decaying at std 0.01 and growing at 0.05 (a sign test); the separations' drift against an independent quadrature; along the second-order transient the mean parameters at the averaged force and the mean error carrying the developing deviation's tracking error; the exact periodic mean against the filtered Gram, the fixed prior's branch and its leak in the bias; the mean function's curvature term; the OU noise's variance and correlation, the law seeing the measurement; the responses (H + σ₀I)⁻¹ and H⁻¹ to a constant force, the quasi-static force against the full second-order mean's noise part) |
| D‴.31 | the instance of D‴.30 with an OU noise: the adaptive prior's bias per (Eν²)² for MR, DR, DR* at k_e = 25, 50, 100, 200 and B = 100, 1000, from the noise's own moments and with the parameters' fluctuation (MR 2397–1.66·10⁵, its parameters' part ≤ 0.03%; DR and DR* ≤ 1.18 < ¼‖∂²_xf‖²_ρ = 1.31, their parameters' part −8.7% to +3.3% except +56% for DR at k_e = B = 100, its own response 3.9·10⁻⁴, 1.1·10⁻⁴, 2.9·10⁻⁵, 7.4·10⁻⁶ at k_e = 25–200 (∝ k_e⁻² on the one-dimensional reference); DR's noise part at its bound ‖g‖²_ρ to 0.01% at k_e = 25; k_e²‖P E[e]‖²_ρ 2394 for MR, 0 for DR and DR* by construction, ‖P E[e]‖²_ρ 3.83 → 0.66; the bias reaching the floor at the std 0.0012, 0.085, 0.27); the prior family at η = σ₀, κ/η = 0, 0.1, 1, 10 and fixed (the weakest mode's response 3269–136, its bias factor 0–0.959, RH + B = I to 1.1·10⁻¹⁴; MR's bias 2397–491, DR's 0.471–0.211) and five priors with non-commuting leaks, R and B from their dynamic blocks (the definitions with S₀ by elimination to 1.9·10⁻¹⁵, RH + B = I to 1.4·10⁻¹⁴, HR + B − I 0.47–0.79); the real loop (RK4 dt = 10⁻³, B = 100, std 0.01, η = ∞, 12 runs of 300 time units, common noise paths): E‖f̂ − f‖² at t = 300 2.43·10⁻⁶ (MR), 2.80·10⁻⁷ (DR), 2.75·10⁻⁷ (DR*), z-scores −2.04 to 0.82; the MR − DR run mean at Ĵ(ȳ_MR − ȳ_DR) to 0.5–5.3%, 0.43–1.83 sampling errors | `exp_eiv_remedy.py` (`--only-AB` recomputes Parts A and B, `--only-B` Part B) → `results/chaos/eiv_remedy.md`; `tests/test_eiv_remedy.py` (the washout identity along a real-loop solution and across a plan switch of the averaged loop; the block's transfer function, RH + B = I along the prior family and, from their dynamic blocks, for non-commuting leaks; the response bound under random leaks, the smallest equivalent perturbation and the projection; the lift exact under DR and DR*, DR's run within the second order of MR's on a common path, the regressor's point; DR's and DR*'s noise forces least-squares with their pointwise bounds, MR's F₁ three orders of magnitude above; the quasi-static forces against the full second-order mean) |
| D‴.32 | the three-bump instance (σ₀ = 0.15σ_crit), the averaged loop: (A) the lift of ϑ^e(w) at the 15 count vectors of m = 15: the largest drift 2.8·10⁻¹⁷ (compact, R = 1) and 6.2·10⁻¹⁰ (Gaussian), the spectrum (iii)'s to 2.8·10⁻¹¹ (the finite differences' error) and 4.6·10⁻⁹ (the Gaussian kernel's cross-class terms) (the separations 1.41·10⁻³–3.77·10⁻³, the slowest collective mode 2.2·10⁻⁴–4.2·10⁻⁴); the Gaussian kernel's bias 2.1·10⁻⁸, 1.3·10⁻³, 0.068 (whitened) at R = 1, 1.5, 2, its first order to 6·10⁻⁸–1.0% (56% off at R = 3); with the fourth bump the margins σ₀w_k − (q_k/w_k)λ⁺_k the Jacobian's to 2·10⁻¹¹ at three count vectors, the class with the margin −1.22·10⁻⁴ at σ₀/8 splitting (3 → 4 classes, the fit error 0.017 → 0.0148); twins merging under the Gaussian kernel, persisting under the compact one; the dead class joining a live one at R = 2 and 1.5 (s ≈ 2.8·10⁴ against (v)'s 1.8·10⁴), not moving at R = 1 (time scale 1.4·10⁷); (I) the twins' contraction 1.2949 against (v)'s 1.2943; far twins (5.28 and 8.19 apart) not moving in 4·10⁴ slow units; opposite-sign twins with a third class within reach not merging in 4·10⁵ (the fit error 3.2·10⁻⁴ → 4.5·10⁻³); (B) warm clouds, m = 15, 60, 240 against M = 510: ‖f̂_m − f̂_M‖² at s = 6000 1.9·10⁻¹⁰–5.4·10⁻⁸ (m = 15, the representable samples), 2.4·10⁻¹⁰–7.0·10⁻¹⁰, 5.1·10⁻¹⁴–2.2·10⁻⁹, as under the adaptive prior, against the fixed prior's floor 3.0·10⁻⁴–2.4·10⁻³, 1.1·10⁻⁴–9.9·10⁻⁴, 4.9·10⁻⁶–8.2·10⁻⁵; seed 3's sample (8, 4, 3) on the boundary of W at 7.0·10⁻⁵–7.1·10⁻⁵ under every leak but the fixed prior; (C) broad starts: the annealed Gaussian kernel at three classes in 9 of 9 runs, the fit error 1.8·10⁻⁷–2.3·10⁻⁶ (the adaptive prior told K₀ = 3: 4.2·10⁻⁹–2.1·10⁻⁶, K₀ = 2: 0.13–0.19, K₀ = 4: 5.4·10⁻⁴–5.6·10⁻³; the σ-modification 0.24; no leak 40–48 classes, the compact kernel without the broad stage 17–22); (J) at a fixed Gaussian scale 4–6 classes (R = 1, 4.8·10⁻⁴–1.6·10⁻³) or three biased (R = 2, 7.8·10⁻⁴–1.4·10⁻³), the schedules 1.5 → 1, 2 → 0.5, 3 → 1 at three classes (2.4·10⁻⁸–5.9·10⁻⁶); (D) the real loop, k_e = 25: the lift's return 7.8·10⁻¹⁶, its 45 slow exponents (36 separations at σ₀w_k/k_e) to 5.5·10⁻⁶ (median 1.5·10⁻⁷), the rms tracking error at t ≈ 2900 5.8·10⁻⁴ and 2.4·10⁻⁴ (m = 15, 60) against the fixed prior's 1.3·10⁻³ and 3.9·10⁻⁴; (E) the dither τ = 10⁻⁵: the within-class spread over the second half 0.068 (compact) and 0.063 (Gaussian) against √((1 − 1/n_k)τ/(σ₀w_k)) = 0.0635, 0.036 under the fixed prior (0.0366), 0.44–0.56 without a leak; (F) the cogging servo: 7–8 classes, the fit error 0.14%–0.35% (the four-kernel designs 1.46% and 0.64%, the best equal-weight six- and eight-kernel fits 9.4·10⁻⁶ and 3.3·10⁻⁷); (G) σ halved from σ₀ to σ₀/64 (the fourth bump, m = 30, three seeds): no split down to σ₀/8 (the predicted threshold 0.121σ₀), the middle class splitting at σ₀/16 in every seed; no class with a positive margin split, all 8 with a margin below −5·10⁻⁵ (= −ln 250/10⁵, the detection threshold) did, 3 of 7 between −5·10⁻⁵ and 0; a split sheds one to three single particles (the cohesion energy of shedding 0.43 of its quadratic part against 0.74 for a halving, compact): 8–10 classes at σ₀/64, five to seven of one particle, the fit error 0.0028–0.0031; (H) two bumps: two classes, the fit error ≤ 1.6·10⁻⁸; the four-bump f_ε: three classes, 0.0174; five bumps (Σ_k|A_k|/a_max = 0.875): 6–9 classes, 0.045–0.118, two or three inner bumps short of mass in every run; (K) three equal bumps three classes (9.4·10⁻⁵–1.8·10⁻⁴, the Gaussian kernel's bias: (L) 2.5·10⁻¹³ under the compact kernel from there), four alternating four (9.6·10⁻³–1.1·10⁻², a slow transient: (L) 6.1·10⁻⁴ after 2·10⁵ more slow units), three heavy (0.95) 3–5 biased (0.058–0.20), two close bumps (±0.45) merged into one stable class (1–3 classes, 0.098–0.13) | `exp_cohesion_leak.py` (`--only A,…,L` recomputes parts, `--merge` merges parts computed elsewhere) → `results/chaos/cohesion_leak.md`; `tests/test_cohesion_leak.py` (the drift −G∇J^coh_m and J^coh_m decreasing; the forces internal and the mean-shift identity; the exact representation stationary with the predicted spectrum at two count vectors; the splitting margin at the best one-atom fit; the twins' contraction against (v)'s rate to 1% under the Gaussian kernel and their persistence beyond the compact range; the kernels' range and scale; the clusters of known configurations; the real loop's lift periodic with the separations' exact contraction) |
| D‴.33 | the instance of D‴.30 (k_e = 25, B = 100, OU noise): the separations' drift at the lift per class — MR [−3.7·10⁻¹, +1.5·10⁻²]; [−2.6·10⁻², +1.0·10⁻²]; [−3.5·10⁻², +3.0·10⁻¹], DR [−1.1·10⁻², −5.0·10⁻⁴]; [−1.5·10⁻², −1.5·10⁻⁴]; [−7.2·10⁻³, −3.4·10⁻⁴], DR* [+1.6·10⁻⁴, +3.8·10⁻³]; [+5.1·10⁻⁵, +5.0·10⁻³]; [+1.1·10⁻⁴, +2.4·10⁻³] —, the threshold std 0.031, —, 0.238 (adaptive prior) and 0.018, —, 0.138 (cohesion leak); against k_e at B = 4k_e (25–800) MR's λ_max 0.2982 → 0.3045, k_e times DR's −0.0038 and DR*'s 0.1245–0.1242, at B = 100 DR's threshold std 0.185, 0.138, 0.121 (k_e = 200, 400, 800), DR*'s 0.238 → 0.113, MR's 0.031 → 0.014, ‖S_k‖ ≤ 0.37 of (ii)'s bound, the floor 0.069; six targets at the lift (DR's centre entries negative in all, the whole drift indefinite for same-sign, wide, edge and close bumps, up to +5.6·10⁻⁴); the stationary point from the exact periodic mean — rms E[e] 6.70·10⁻⁴ (DR), 6.44·10⁻⁴ (DR*), the parameters' part 6.32·10⁻⁴, 6.27·10⁻⁴, the drift indefinite ([−4.35, +3.68]·10⁻⁵ in class 1), under the fixed prior DR's negative definite and DR*'s positive definite (the threshold std 0.435); DR's transient over 5.2·10⁵ time units (the largest eigenvalue first positive at t = 4103, 2947, 2714, at most +1.05·10⁻⁴), the cross term −8% (DR) and −17% (DR*) of the class-2 window drift over [0, 50), DR*'s threshold std there 0.277; the real loop (four antithetic pairs per family, 300 time units): MR +4.69·10⁻⁴, +1.25·10⁻³ (std 0.05, 0.1; all runs growing), DR along u −2.9651·10⁻⁴, −3.3650·10⁻⁴ (the noise part 0.98, 0.97 of the transient's prediction), along DR's top eigenvectors faster than σ₀/k_e by 0.36%, 0.14%, 0.26% (std 0.05), DR* along its top eigenvector slower by 3.2% and 12.6% (0.97, 0.96 of the prediction) | `exp_eiv_remedy.py` (`--only-Dpred` recomputes Part D's predictions, `--only-D2` also the runs along the laws' own top eigenvectors) → `results/chaos/eiv_remedy.md`; `tests/test_eiv_remedy.py` (DR's drift MR's mean-error term alone at the lift, negative definite in every class, MR's 0.30 in class 3; DR*'s O(1/k_e) and against an independent quadrature to 10⁻⁴; at the stationary point the exact periodic mean's two parts, the parameters' part dominant, the drift indefinite; the period maps against `second_order_transient` to 10⁻⁷ and DR's drift turning indefinite) |
| D‴.34 | the instances of D‴.32: (A) fixed positions (the five bumps' atoms, the masses (0.15, 0.30, 0.20, 0.25, 0.10)): the fit error at s = 4·10⁴ 0.388 (γ = 0; the least error at those masses 0.387), 2.9·10⁻⁴ (γ = 10⁻³), 3.0·10⁻¹⁴ (γ = 10⁻²); the amplitudes scaled by 1.25 (Σ|A_k|/a_max = 1.094): 0.071 at both rates (the ℓ¹ ball's least error 0.0684), every |a|/a_max ≥ 0.992; (B) the five bumps from broad starts, m = 60: γ = 0 6–9 classes, 0.045–0.118; γ = 10⁻³ 5–7, 0.0076–0.0139; γ = 10⁻² five in 3 of 3, the centres the bumps' to 0.01, 6.0·10⁻⁴–1.4·10⁻³ (falling by 2–4% per 1000 slow units; seed 1 continued: 2.80·10⁻⁴ after 10⁵ more under the Gaussian kernel, 2.85·10⁻⁴ after 4·10⁴ under the compact one); (D) the instance, γ = 10⁻²: three classes, 2.4·10⁻⁷–9.7·10⁻⁷ (without, seeds 0–2, 1.9·10⁻⁷–9.0·10⁻⁷); (C) σ₀/64 from the lift of ϑ^b(q): 8–10 classes with and without the replicator (0.0055–0.0075 against 0.0059–0.0078); (F) three heavy bumps: three classes in 3 of 3, 7.8·10⁻⁸–1.6·10⁻⁷; two close bumps: 1–3 classes, 0.068–0.13; (G) the class-level stationary points at R = 1: B's seed 1 under the Gaussian kernel 3.85·10⁻⁴ (masses held) and 2.63·10⁻⁴ (u_k = λ), under the compact one exact (2.4·10⁻¹⁶); three equal bumps 9.39·10⁻⁵ and 9.17·10⁻⁵, the masses moving toward (0.567, 0.217, 0.217), exact under the compact kernel; (E) the real loop at the instance's atoms (k_e = 25, γ = 2): V 0.0418 → 2.1·10⁻⁷ monotone, V̇ + k_ee² ≤ 1.1·10⁻¹⁶ max k_ee², Σβ − m ≤ 2.4·10⁻¹⁴; the averaged linearisation at the lift of ϑ^e(w), (6, 5, 4): 15 eigenvalues below 10⁻⁹, the others in [−0.1834, −3.08·10⁻⁴], at γ = 0 D‴.32's spectrum to 2.1·10⁻¹¹ | `exp_mass_transfer.py` (`--only A,…,G`) → `results/chaos/mass_transfer.md`; `tests/test_mass_transfer.py` (the gradient-flow identity by central differences of J, γ = 0 Proposition D‴.32's drift; the weights' drift against an independent quadrature of the fitness; the real loop's relative-entropy V̇ = −k_ee² pointwise from the vector field and Σβ conserved; at fixed positions Σβ = m kept, the exact fit reached with the replicator and not without it, and above the ℓ¹ ball's least error when Σ|A_k|/a_max > 1; at the lift the zero modes the predicted vectors — the within-class redistributions and one gauge per class —, a separation at −σ₀w_k with γ > 0, the rest D‴.32's spectrum) |
| D‴.35 | Proposition D‴.32's resolution instance (the lift of ϑ^b(q), (10, 10, 10), compact R = 1): (A) λ⁺_k = 1.00·10⁻⁴, 2.86·10⁻⁴, 1.08·10⁻⁴, c₃ = 1.42·10⁻⁵, −1.257·10⁻³, 1.20·10⁻⁴, the middle class first at 0.1213σ₀; (B) the split states q = 1–4 by Newton within 1.3%–15% (μ = 2·10⁻⁶) and 0.3%–3.0% (5·10⁻⁷) of the prediction, the eigenvalues on the separations along v +1.2511μ (×8) and −1.023μ for q = 1 at 5·10⁻⁷ against +1.25μ and −μ, 8, 7, 6, 5 unstable directions; (C) six seeds at μ = 5·10⁻⁵ and at σ₀/16: the first excursion's particle the quadratic equation's in 12 of 12 (linear 6, flipped c₃ 2), its time to 1.1% (linear 22.4%), on the side sign(c₃)v in 10 of 12 (the other two linear growth below 2μ/\|c₃\|); beyond the range at σ₀/16 a single particle in 3 of 6 and a bound pair in 3 of 6 (the top two within 1.2%–11%), the main piece's margin positive at the end; (D) Part G's schedule, seeds 0–2: 8–10 classes (the cohesion leak), 6–7 (replicator γ = 0.1, the first escapee's mass 0.105–0.106); (E) the replicator at γ = 1: 5, 6, 8 classes, the first escapee's mass 0.134–0.136, the fit error 0.0029–0.0030; (F) the mass-free leak continued to σ₀/256: 4 classes at σ₀/64 in 3 of 3 (the fit error 0.0149), then 9, 10, 9 at σ₀/128 and 9, 8, 9 at σ₀/256, the fit error 0.0028 — the cohesion leak's end state | `exp_split_bifurcation.py` (`--only A,…,F`) → `results/chaos/split_bifurcation.md`; `tests/test_split_bifurcation.py` (the full drift along v against the reduced quadratic equation with c₃, the remainder cubic; the single-particle split state a saddle with n − 2 directions at +μn/(n − 2) and the radial one at −μ; the predicted particle leaving first on the side sign(c₃)v, at μ = 10⁻⁴ with the other classes' margins positive; the reduced crossing against the linear closed form and an independent integration; the mass-free leak's blocks) |
| D‴.36 | Proposition D‴.32's three-bump instance, one particle per class at ϑ^e(1/3), τ = 10⁻⁴: (A) the predicted D_kk = 2.85·10⁻⁸, 1.63·10⁻⁸, 1.81·10⁻⁸ at γ = 10⁻² (4.85·10⁻¹⁰ to 4.23·10⁻⁷ for the first class over γ = 10⁻³–10); (B) 64 seeds, 5·10⁴ slow units: the variance's slopes (2.69 ± 0.58, 2.21 ± 0.47, 1.62 ± 0.38)·10⁻⁸, the ratios 0.95, 1.36, 0.90; the means 0.3212 ± 0.0046, 0.3445 ± 0.0041, 0.3342 ± 0.0038 (the drift term); (C) 10 seeds, 2·10⁵: the spread 0.057, 0.053, 0.052 (√(Ds) 0.065, 0.049, 0.052), α_1 up to 16.1 (|a|/a_max = 0.9994, against 0.80 at γ = 0), the fit error the dither's (0.0253 against 0.0269); (D) the leak at ε = 10⁻⁴, 10⁻³: the spreads 0.0152, 0.0103, 0.0117 and 0.0044, 0.0032, 0.0034 against the predicted 0.0139, 0.0099, 0.0105 and 0.0044, 0.0031, 0.0033, the means within 1.5 standard errors of 1/3 (4 seeds; the 16-seed test in D‴.37 (G)) | `exp_gauge_dither.py` → `results/chaos/gauge_dither.md`; `tests/test_gauge_dither.py` (the entropic leak a Fisher–Rao gradient; the two-class diffusion's closed form; the leak keeping the exact representation stationary and making the gauge stable; the leak's spread against the two-class Ornstein–Uhlenbeck closed form; γ = 0 leaving the weights at 1; D against the projector of the code's own drift Jacobian; run_dithered's noise 2τG dt) |
| F (ii) along (L) | particle path length (≥ Wasserstein length) ≤ λ_G^{1/2}L₁∫‖e‖ at every sample | `T::test_theorem_F_wasserstein_length_bounded_by_integrated_error` |
| F (ii′) along (L_σ) | W(μ_{t₁},μ_{t₀}) ≤ ∫(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰)) on a grid of pairs; a tight instance with e ≡ 0 (ratio ≥ 0.99), where the motion is W ≈ 1.05 | `T::test_theorem_F_leak_clause_along_ot_sigma_mod`, `T::test_theorem_F_leak_clause_is_tight_when_the_leak_pays_for_all_motion` |
| E linearisation | Jacobian of the implemented loop equals A(t); QA + AᵀQ = diag(−k_eI, 0) | `T::test_theorem_E_linearisation_and_lyapunov_identity` |
| E convergence | the matched error decays and the labelled error does not (qualitative); sliding-window excitation and fitted rates | `T::test_theorem_E_matched_error_converges_labelled_does_not`, `exp_pe_convergence.py` |
| Lemma E.0 | p_lo‖w‖² ≤ V ≤ p_hi‖w‖² and V̇ ≤ −γ‖w‖² at random (t, w) for four regressors, closed-form and optimised ε; the envelope along trajectories; explicit rate ≤ exact Floquet rate | `Ex::test_lemma_E0_pointwise_inequalities_hold_on_random_states`, `Ex::test_lemma_E0_envelope_and_floquet_ordering` |
| Lemma E.0(iii), Remark E.0a | the ω⁻² law of the explicit rate and of the exact rate; the saturation of exact/explicit in the slow-adaptation limit; the tables | `Ex::test_remark_E1_omega_scaling_of_the_explicit_rate`, `Ex::test_slow_adaptation_limit_ratio_is_bounded`, `exp_explicit_rate.py` → `results/explicit_rate/table.md` |
| E, explicit constants | on the implemented loop, from a start inside the explicit ball, ‖w(t)‖ ≤ (p_hi/p_lo)^{1/2}e^{−λt}‖w(0)‖ at every sample (L_J estimated from samples) | `Ex::test_theorem_E_explicit_nonlinear_envelope_on_the_implemented_loop` |
| Corollary E.0′.2 | the explicit envelope of Theorem E with Lemma E.0′'s constants on the implemented single-kernel loop at k_e = 40: 1.5 e-folds inside the envelope from a start in the explicit ball (dt = 0.005, below the discretisation floor); λ′ > 10 × Theorem E's rate from Lemma E.0 and < 4 × the exact linear rate; the constants table | `Ex::test_corollary_E0prime_2_envelope_on_the_implemented_loop`, `exp_explicit_rate.py` (table 5) |
| Corollary E.0′.3 | the explicit bound on the stationary parameter spread of the dithered (LTV) dominates the exact periodic covariance (sinusoid: ratios 2.7 → 1.3 (δ = 0.2) and 2.5 → 1.2 (δ = 0.1) at k = 20 → 200, falling toward 1/(1 − δ); the Theorem E ideal: 15–18, the price of q = 15 directions at the worst excitation); the exact value equals the averaged k·tr X within 2/k² | `Ex::test_corollary_E0prime_3_dither_covariance_bound_against_the_periodic_lyapunov_solution`, `exp_explicit_rate.py` → `results/explicit_rate/table.md` (table 6) |
| Lemma E.0″ | on Ψ = (sin t, cos √2 t) (not periodic) with window 4π: the exponential-average identities on the grid (10⁻⁴), the conditions at k = 20, 50, 200 with δ₂ ≤ 5%, the certificate on 400 times, the guaranteed rate below the exact decay rate (top Lyapunov exponent) within a factor 2.6 and improving with k; the experiment: exact/guaranteed 2.5 → 1.3 on the quasi-periodic regressor and 5.0 → 3.7 on a slowed Lorenz-driven regressor (k = 100–400), the certificate holding at every admissible row | `Ex::test_lemma_E0dprime_quasi_periodic_regressor_certificate_and_rate`, `exp_pe_general.py` → `results/explicit_rate/general_pe.md` |
| Lemma E.0′ | the certificate sup_{t,x}V̇/V ≤ −γ (a generalised eigenvalue problem on a time grid) on Ψ = (sin t, cos t) at k = 20, 50, 200 and on the slow-adaptation ideal at g = 0.01, k = 400; N₂ = 0 numerically; the guaranteed rate below the exact Floquet rate and within a factor 2 (sinusoid) and 1.5 (ideal) of it, above Lemma E.0's by more than 10⁴ on the ideal; the table of k and δ | `Ex::test_lemma_E0prime_sinusoidal_regressor_certificate_and_sharpness`, `Ex::test_lemma_E0prime_pays_for_ill_conditioning_once_on_the_theorem_E_ideal`, `exp_explicit_rate.py` → `results/explicit_rate/table.md` (table 4) |
| G(a), two scales | direct check (cosine ≥ 0.998 in the G⁻¹ metric, `violating`; the weight block fails in `satisfying`, as explained in §6b); L̄ at period-mean parameters non-increasing | `exp_variable_projection.py` → `table_violating.md`, `table_satisfying.md` |
| G(b), (b′), (H2c) | VarPro holds when (H2b) and S ≪ 1 hold; it fails when S = O(1) although (H2b) holds; the induced-metric flow tracks the closed loop at every S (started on the slow manifold; gains ÷ 4 isolates the averaging error) | `exp_variable_projection.py` → `table_satisfying.md`; `exp_slow_manifold_slope.py` → `table_slope.md` |
| H(i), algebraic identity | with e := r/k_e and x := x_m, the averaged drift equals −G∇V (code gradients against finite differences, 1e-7); does not test averaging | `Lg::test_averaged_drift_is_the_mean_field_langevin_drift` |
| Remark H.0 | the frozen-μ averaging error decays like 1/k_e, n = 1 and n = 2 | `Lg::test_averaging_error_of_the_formal_reduction_decreases_like_one_over_ke`, `Av::test_averaged_drift_error_decays_like_one_over_ke_for_an_orbit_enclosing_area` |
| H, need for entropy | finite-m real loop: a noiseless symmetric start stays bit-identical (and moves); dither separates it | `Lg::test_noiseless_symmetric_initialisation_never_separates_but_dither_does` |
| H(ii), oscillation bound | |V_L[μ](ϑ)| ≤ K_Φ(F_ρ + K_Φ) for random and adversarial μ, ϑ | `Lg::test_first_variation_bound_used_in_holley_stroock` |
| H(ii) constant | formula transcription only | `Lg::test_lsi_lower_bound_formula` |
| H, global behaviour | collapsed against spread starts, with and without dither, two (σ, τ) settings | `exp_langevin_global.py` |
| H, speed–bias trade-off | escape time and loss floor against τ; log(1/τ) and Arrhenius fits, seed-averaged and per seed (neither identifies the mechanism) | `exp_langevin_tradeoff.py` |
| H, reference values | L̄(μ*_τ) by a concave dual on a grid; the stationary state of the exactly averaged loop; their offset in the scaling of H.6 | `exp_gibbs_reference.py` → `results/langevin_global/gibbs_reference.md` |
| H.7 | expansion, Stokes area identity, parity under reversal, non-reciprocity of the exact second variation | `Av::test_frozen_measure_expansion_converges_at_rate_one_over_ke`, `Av::test_circulation_kernel_antisymmetric_vanishes_in_1d_and_equals_the_enclosed_area_integral`, `Av::test_shift_correction_is_even_and_lag_correction_is_odd_under_orbit_reversal`, `Av::test_exactly_averaged_dynamics_is_non_reciprocal_for_n1_and_n2` |
| H.7a, oscillation term | value against an independent one-period near-identity flow (within 3% at k_e = 80); odd under orbit reversal (1e-6) | `Av::test_oscillation_correction_matches_a_one_period_flow`, `Av::test_oscillation_correction_is_odd_under_orbit_reversal` |
| I(i), (ii) | aᵀB(ξ,ξ')a' equals the antisymmetric part of S_even + K (1e-10); A(ϑ,ϑ_z) = 0 and the cyclic sum over (ϑ₀,ϑ₁,ϑ_z) equals A(ϑ₀,ϑ₁) (1e-14) | `NR::test_bilinear_form_B_reproduces_the_antisymmetric_second_variation`, `NR::test_zero_feature_collapses_the_cyclic_criterion_to_a_pairwise_one` |
| I(iii) | n = 1: A ≡ 0 at perfect fit (< 1e-13); A ≠ 0 with a residual, equal to the closed form; the positivity witness | `NR::test_n1_reciprocal_iff_residual_vanishes_on_the_orbit`, `NR::test_positivity_witness_n1` |
| I(iv) | the closed form of κ (1e-10); the second component of F_s((0.5, 0)) is positive on the ellipse; κ ≡ 0 on a zero-current orbit; with f ≡ f̂_μ on ℝ² (Dr_μ ≡ 0) B = 2κI, A ≠ 0 on the ellipse and A ≡ 0 on the segment; with f = f̂_μ + he₁, h zero on the ellipse, the curl term survives and equals −∫φφ'∂₂h dρ·(e₁e₂ᵀ − e₂e₁ᵀ) | `NR::test_circulation_kernel_closed_form_and_vanishing_for_a_zero_current_orbit`, `NR::test_perfect_fit_n2_leaves_only_the_circulation_kernel`, `NR::test_perfect_fit_on_the_orbit_alone_keeps_the_curl_term` |
| J, averaged objective | ∇J_m against finite differences (1e-7); the critical point is a non-degenerate minimum with two distinct kernels | `Pe::test_gradient_of_J_m_matches_finite_differences_and_critical_point_is_nondegenerate` |
| J(a), (b), (d) | Newton shooting finds the periodic orbit (residual < 1e-10); ‖Θ̄_p − Θ*‖ and max‖e‖ halve from k_e = 20 to 40; the vector error of the first-order prediction of the offset halves as well (it is 0.55 at k_e = 40 in the fast-sweep test instance, and 0.041 at k_e = 40 in instance A of the experiment); the parameter Floquet exponents match (1/k_e)λ(G_m^{1/2}∇²J_mG_m^{1/2}) with a maximal relative error below 2% at k_e = 40 that halves with k_e; both instances at k_e = 10–160 in the experiment; the error-channel exponent recovered from Liouville's trace identity is k_e − L with |L| < 0.2 | `Pe::test_theorem_J_periodic_orbit_offset_and_floquet_exponents`, `exp_periodic_orbit.py` → `results/periodic_orbit/table.md` |
| K(i), (ii) | the Hessian of J_m at a collapsed configuration equals I_m ⊗ (σ₀G⁻¹ − R₂) + (1/m)(11ᵀ) ⊗ M₂ (finite differences, 1e-5 relative); σ_crit solves σ = λ_max(G^{1/2}R₂(ϑ̄(σ))G^{1/2}) (1e-8) and the separation block is singular there; the separation Floquet exponents of the collapsed orbit come in (m − 1)-fold multiplets | `Pe::test_proposition_K_hessian_splitting_and_critical_leak`, `exp_orientation_stability.py` |
| K(iii), (iv) | the orientation difference of the softest separation exponent at σ₀ = 0.5: slope −2.18 against k_e, within 14% → 3% of the prediction at Θ̄₂ for k_e = 10 → 80, while the prediction at Θ* converges to a value 30% off; at σ_crit both orientations are unstable, predictions within 14% at k_e = 40; at the designed σ₀ the signs differ for k_e = 10 and 20 | `exp_orientation_stability.py` → `results/orientation_stability/table.md` |
| J′ averaged data, solver | the single-kernel J, ∇J, b_fr, b_osc against the general H.7/H.7a code (2·10⁻⁵); the weighted Poisson solver against the Ornstein–Uhlenbeck closed forms (g, the mean shift (Gλ)⁻¹b, the Green identity, zero current for a gradient drift); the stationary Fokker–Planck solver of the first-order SDE: current ratio → 1 linearly in ε, ‖m − π(1 + εg)‖₁ ∝ ε² | `tests/test_dither.py` (3 tests) |
| J′(b), (c), (d) | real loop, 2·10⁴ replicas, k_e = 10–80, both orientations: the stationary means differ from the Gibbs values by k_e^{−0.87} (α), k_e^{−0.83} (γ), k_e^{−1.09} (αγ); the first-order prediction removes 93–97% of the discrepancy at k_e = 10 and 70–85% at k_e = 20–40, the residual at the Monte Carlo floor at k_e = 80 (second order not resolved); the orientation-odd shift matches its prediction to 7–8% at k_e = 10 and within the 10–30% resolution at k_e = 20–80; k_e²E e²/(2π(L̄)) − 1 = 0.5–0.7/k_e | `tests/test_dither.py::test_real_loop_means_are_closer_to_the_first_order_prediction_than_to_gibbs` (k_e = 20), `exp_dither_stationary.py` → `results/dither/table.md` |
| J′(e), J′.2 | the circulation of the real loop is below the Monte Carlo resolution at first order (k_e ≥ 20) and dominated by higher orders at k_e = 10; the current and density expansion of the first-order averaged SDE checked on the grid (ratios 1.47 → 1.03, ε² remainder 0.31); the antisymmetric part of D(G⁻¹b) against the formula of J′.2(iii) to six digits (referee round 6) | `exp_dither_stationary.py` → `results/dither/table.md`; `review/REVIEW_round6_Jprime.md` |
| Corollary J′.3 | the Gibbs measures π_m (MALA, 4000 replicas, m = 1–32) converge to the mean-field minimiser μ*_τ (fixed point to 10⁻¹¹) with a 1/m law: slopes −1.01 (E[γ₁]) and −0.99 (E[L̄(μ_Θ)], rescaled deviation 0.058 ± 0.002 at every m), E[α₁²] from 0.54 to 0.92 against 0.935 with rescaled deviation −0.41 ± 0.01; the m = 1 row agrees with the quadrature of π₁ | `exp_gibbs_m.py` → `results/dither/gibbs_m.md`; `tests/test_dither.py::test_corollary_J3_mean_field_minimiser_and_gibbs_sampler_consistency` |
| H.7 + H.7a in the real closed loop | the orientation dependence of the learned kernels (same ρ) is predicted to within 3% in size and 11% as a vector, cosine 0.994 at k_e = 40; the measured difference scales like k_e^{−1.01} | `exp_orbit_orientation.py` → `results/orbit_orientation/table_sigma2_slow16.md` |

Not checked in the repository:

- B.5 along (L), and B.6 and B.7;
- the Filippov form of C(i), the inequality of C(v), and Proposition C.5's D⁺V_OT = σ_m(δ − δ²);
- D(b), D's constant, and the one-sided constant of Remark D.1;
- the sharpness of the explicit rates of E beyond the tables of Remark E.0a and Lemma E.0′; Lemma E.0′ below its threshold k₀ (where it makes no claim) and for non-periodic PE; Corollary E.2 and Remark E.3;
- Theorem I(v) beyond the two instances (the analytic dichotomy is a proof only);
- Theorem J(c) (the basin; a single basin run is reported in `review/REVIEW_round5_J.md`);
- Proposition K(iii)–(iv) beyond the instance of `exp_orientation_stability.py`, and the second-order normal form it rests on;
- Theorem C′(iv) (its hypotheses are not met on any tested instance) and the convergence of the pure leak to a critical
  point beyond the tested instances; Theorem D′'s constant, its per-particle bound and its k-particle convergence (only
  the rate of (b) is checked);
- H.8, which is a proof only;
- Conjecture H.6.

---

## 8. Research program

Each branch states what is proved, what is open, and the first concrete step.

1. **Infinite-horizon averaging for McKean–Vlasov closed loops (Conjecture H.6).**
   - *Proved:* the averaged system converges globally and exponentially (H(ii)); the exact first correction at frozen
     measure (H.7); input-to-state stability of the averaged system (H.8).
   - *Open:* a near-identity averaging transformation for the stochastic McKean–Vlasov loop that keeps the
     frozen-measure correction valid along the moving measure, uniformly in time.
   - *Done (Theorem J):* the deterministic case (τ = 0, σ = σ₀/k_e) at fixed m: near every hyperbolic critical point
     of J_m the real loop has a unique periodic orbit with the same stability type, k_e-uniform basin, and the
     first-order location of Proposition H.7 and Remark H.7a.
   - *Done (Theorem J′):* the dithered case at fixed m. The invariant law of the loop has the Gibbs measure exp(−J_m/τ₀)
     as its parameter marginal up to O(1/k_e), with the first-order correction driven by the same drift b as Theorem J; the
     fast error and the phase are eliminated by a Poisson-equation expansion of the generator rather than a transformation
     of the Fokker–Planck equation (§6f).
   - *Next steps:* uniformity in m — the mean-field limit of the stationary law (the Gibbs measures of J_m converge to the
     minimiser of the mean-field free energy F_τ of Theorem H, a propagation-of-chaos question for Gibbs measures which is
     where Conjecture H.6 now lives), the second-order correction, and quantitative constants in place of the O(ε²).

2. **Non-reciprocal mean-field dynamics from closed-loop learning (new; Proposition H.7).**
   - *Proved at a frozen measure:* the averaged adaptive dynamics is a semi-gradient flow plus a lag term
     (Proposition H.7). Along the real loop a second-order averaging term is added (Remark H.7a).
   - Its interaction has two antisymmetric parts:
     - the performative defect A_even of the state shift, which is present for every n and is the only part for n = 1;
     - the circulation kernel K, an integral over the region that the reference orbit encloses (n ≥ 2).
   - *Proved (Theorem I):* reciprocity at order k_e⁻² is equivalent to A_μ ≡ 0; for n = 1 it holds exactly at perfect
     fit on the orbit, and for n ≥ 2 at perfect fit exactly when the orbit carries no circulation and the residual's
     Jacobian is symmetric along the orbit; no open set of measures is reciprocal, for any n. The performative defect
     is bounded by the residual and the antisymmetric part of its Jacobian; the circulation defect depends on neither.
   - *Open questions:*
     - Are stationary states unique beyond the perturbative regime?
     - Can circulation be *designed*, choosing the orbit's orientation and area, to break symmetric traps or to speed
       escape from the collapsed configuration? (Partly answered by Proposition K, below: the orientation decides the
       stability of the collapsed configuration in a window around the critical leak.)
     - Is there a Lyapunov functional other than F_τ for the non-reciprocal flow?
   - *Context:* the closest analogues are performative prediction [S14], semi-gradient TD learning [S10], [S11] and
     non-reciprocal phase transitions [S12].
   - *Done (Remark H.7a):* deterministic n = 2 closed-loop runs with the reference traversed both ways (same ρ). The
     full first-order theory predicts the O(1/k_e) difference of the learned kernels to 3% in size and 11% as a
     vector.
   - *Done (Proposition K):* at a collapsed configuration the orientation-odd k_e⁻² terms shift the separation Floquet
     exponents of the real loop by O(k_e⁻²), with a closed-form prediction from the averaged data; at the critical leak
     σ_crit the direction of traversal decides whether the kernels separate, and the σ₀ at which it does can be computed
     in advance (§6e).
   - *Done (Theorem J′):* the test with dither. The orientation-odd part of the stationary parameter means is 2ε∫Θ g_odd dπ,
     with g_odd driven by the odd part of b; measured on the single-kernel instance of §6f (where, in one state dimension,
     the whole odd effect is b_osc).
   - *New (Proposition J′.2): the stationary state is a non-equilibrium steady state.* The first-order averaged dynamics
     carries a stationary current j₁ = π(b − τ₀G∇g) exactly when b is not a G-gradient, which for m ≥ 3 is always the case;
     its entropy production rate is ε²τ₀⁻¹π(‖b_⊥‖²_{G⁻¹}); and the real loop's stationary parameter process is not
     reversible. The Lyapunov-functional question for the non-reciprocal flow is therefore a question about a
     non-equilibrium steady state (its free energy is not a Lyapunov functional of the real loop at order ε).
   - *Next steps:* the design of the orbit's area (not only its orientation) as a control input; the orientation effect at
     non-collapsed configurations with several kernels; measuring the current of the real loop (below the Monte Carlo
     resolution in §6f's instance; a variance-reduced estimator or a stronger instance is needed); the mean-field limit of
     the current.

3. **Uniform-in-time propagation of chaos with a plant in the loop.**
   - Theorems D and D′ are finite-horizon; D′ gives the m^{−1/2} rate for the dithered loop by synchronous coupling, with
     constants independent of the temperature and an exponent independent of the dimension for the state and the learned
     function.
     Uniform-in-time results exist for mean-field Langevin dynamics [S5] and for shallow networks [36].
   - *Done (Theorem D″):* uniform in time when the leak dominates the loop gain of the adaptive law (a small-gain
     condition between the plant channel and the particle channel, explicit in Lipschitz constants). The certified regime
     is leak-dominated: the network does not learn (without dither it collapses to one kernel; with dither it stays in the
     cloud of the prior, where the bound has content: the coupled distance is of order m^{−1/2} and bounded in time). The
     result marks where the contraction method stops rather than settling the question.
   - *Open, and now with a negative indication:* the learning regime (kernels separated, leak below the critical value of
     Proposition K), where the plant couples all particles through e and (by H.7) the averaged interaction is not
     reciprocal. The long-horizon run of `exp_chaos_uniform.py` (σ = σ_K/3, 1000 s) shows the coupled distance growing
     threefold between the windows (100, 200] and (500, 1000] while keeping the m^{−1/2} law at each time: on that horizon
     only the finite-horizon Theorem D′ is supported. Whether the growth saturates (the leak's time scale is 1/σ = 2000 s),
     and whether a uniform statement holds at all in the learning regime, is the open question; the first draft's claim of
     empirical uniformity was based on a run that was in fact a slow collapse and is withdrawn (review round 6).
   - *Done (Theorem D‴, §5d):* for the averaged slow loop at τ = 0 the question is settled where the fitted configuration is
     non-degenerately stable: the interaction term of the synchronous coupling is a negative square at linear order (flat
     convexity), the self term carries the residual-weighted curvature at the actual residual, and uniform-in-time propagation
     of chaos at the rate 1/m (mean square) holds under the strict second-order optimality condition along the mean-field flow,
     with explicit constants carrying the transient's e^{2Λ₁}; the window κ_R[μ_∞] < σ₀ < σ_crit is where learning and
     uniformity coexist. The cold start of Theorem D″'s experiment was outside it on two counts (the expanding transient, and
     a connected support threading the ridges, where the coupling has a floor of pairs in different basins). At finite m the
     stationary coupled distance is the weight-sensitivity law of Proposition D‴.5, and Corollary D‴.7 carries the uniform statement
     to the real loop on a warm start through Theorem J.
   - *Done (Corollary D‴.14, Proposition D‴.15):* in the real loop the collapsed-by-class subspace is invariant and carries the weighted
     K₀-atom loop, so Theorem J's orbit at a replicated configuration is, under its uniqueness-ball proviso, the lift of one K₀-atom orbit — with the OT leak the same orbit
     for every commensurate size (the real loop's finite-size floor is zero), with the nearest-atom leak the orbit of the realised weights
     (Corollary D‴.7's floor is exactly the distance between two weighted K₀-atom orbits); the Floquet spectrum of a lifted orbit is the
     K₀-atom orbit's together with each class's single-particle (ghost) spectrum of multiplicity n_k − 1, independent of m — the spectral
     part of Corollary D‴.7's uniformity question.
   - *Done (Theorem D‴.18):* Theorem J uniformly in m at a replicated configuration with admissible weights (margins λ₀ on a compact set 𝒲), in the per-particle norm —
     one threshold k₀, one basin radius ρ and one rate for every size; for k_e ≥ k₀ the proviso of Corollary D‴.14(ii) is replaced by
     this m-uniform statement. In the mean-square norm no m-independent basin exists for the nearest-atom leak (Remark D‴.19).
   - *Done (Theorem D‴.21):* the measure-level basin of the OT leak, for the averaged loop. The OT energy satisfies a Polyak–Łojasiewicz
     inequality whose constant is set by the residual-weighted curvature along the plan segments and does not depend on m (sharp, 2σ₀, at a
     matching prior); the excess energy decays exponentially from an m-independent sublevel set (tiny on the instance) and a posteriori
     wherever the plan-segment margin is positive, which covers configurations with a fixed fraction of the particles in other cells; under
     a curvature hypothesis at the stationary residual the replicated configuration is the unique global minimum among the empirical
     measures, while the nearest-atom energy has a spurious minimum at every admissible count vector.
     *Open:* the transit of the real loop (averaging with plan switches) and the existence of the mean-field OT flow (Remark D‴.22);
     the dithered real loop with the OT leak (Remark D‴.16, checked on the averaged loop: under the OT leak a ridge crossing is reversible
     — the particle stays assigned to its own atom — while under the nearest-atom leak dither equalises the masses only on the Arrhenius
     time scale) — Theorem J′'s numerics with the two leaks.
   - *Open:* the measure-weighted form of hypothesis (H) (a bad set of small μ_s-mass, the bootstrap on the good particles, the
     flip fraction p_m and its m-law); the transfer to the real loop at finite k_e through Theorem J′'s O(ε) corrections,
     uniformly in time; dither (the entropy method of [S5], whose hypotheses Theorem H verifies, against the coupling's floor).

4. **Geometry of two-time-scale adaptation (Proposition G(b′)).**
   - *Proved formally, and checked:* the slow geometry follows the gradient flow of the profile loss in the metric
     induced on the slow manifold, and the slope S measures the departure from plain variable projection.
   - *Open:* a rigorous long-horizon statement; the O(1/R) corrections, which involve D²a* and ∂_ξM; gain designs that
     shape the induced metric.

5. **Noiseless global convergence.**
   - *Done (Proposition K):* the local criterion for the collapsed trap: a collapsed configuration is a strict local minimum
     of the regularised loss iff σ₀ > σ_crit(ϑ̄) = λ_max(G^{1/2}R₂G^{1/2}), where R₂ is the Hessian of the first-variation
     potential at the atom; below σ_crit the noiseless flow leaves the trap along the softest separation mode, without dither,
     and the residual at the atom is what drives the separation (R₂ is linear in it).
   - *Open:* a closed-loop separation condition under which the noiseless flow avoids spurious stationary points,
     despite the failure of the Chizat–Bach boundary conditions. Candidates are a compactified width coordinate and
     reflecting boundaries (Remark H.3).

6. **Unbalanced transport, meaning online structure adaptation.**
   - Treat the weights as masses (Hellinger–Kantorovich / Wasserstein–Fisher–Rao) to obtain birth and death of kernels
     with Lyapunov guarantees. The conic particle descent of [49] and [55] is the static analogue.

7. **Excitation theory in measure space.**
   - *Done (Lemma E.0):* a quantitative, uniform UES estimate for the skew-symmetric class, explicit in
     (k, α, T, Ψ_M, Ψ_D), which reproduces the ω⁻² law of Remark E.1.
   - *Done (Lemma E.0′):* a strict Lyapunov function whose rate is the averaged rate λ_min⟨ΨᵀΨ⟩/k up to a factor 1 − δ,
     for k above a threshold set by the peak factors of Ψ, Ψ̇, Ψ̈ (second-order slaving of the error to its slow manifold
     and a second-order averaging correction, whose mean part vanishes identically); the ill-conditioning of the averaged
     Gramian is paid once, as in the exact rate.
   - *Done (Corollary E.0′.3):* the same Lyapunov function bounds the stationary parameter spread of the dithered
     linearised loop, limsup E‖ζ‖² ≤ (1 + O(ε))k‖Σ‖_F²/(2(1 − δ)λ_min), within 1/(1 − δ) of the exact (averaged Langevin)
     value on isotropic regressors and paying the condition number in the trace otherwise: an explicit design formula for
     how much dither a given tracking budget tolerates.
   - *Done (Lemma E.0″):* the non-periodic case. The period mean is replaced by a forward exponential average over a
     window T (whose derivative identity TṀ̄ = M̄ − M gives the antiderivatives in closed form), the peak factors become
     time-local, and one residual second-order term N̄₂ appears, costing δ₂ = O(1/k) in the rate; the guaranteed rate is
     (1 − δ − δ₂)λ_T/k with λ_T the exponentially weighted PE level of the window. Checked on a quasi-periodic regressor
     and on a regressor driven by a Lorenz orbit, where no periodic construction applies.
   - *Done (Theorem D‴.23(ii)):* at an exact representation the measure-space excitation spectrum h_j — the whitened weighted
     Gauss–Newton spectrum of the K₀-atom loss — sets both a fixed prior's bias, the fraction σ₀/(h_j + σ₀) of its error kept in mode j,
     and an adaptive prior's rates λ_j^− < min{h_j, η} (h_j, η > 0): removing the leak's bias costs the excitation of every collective
     mode, and the leak's a-priori bound on the parameters (the leak toward an adaptive prior is an internal force, and the solutions can
     run off where no exact representation exists; projecting the prior onto a compact convex set restores the bound, Proposition D‴.27) — the measure-space form of the classical trade-off between a leak's robustness and
     persistent excitation (an unexcited mode, h_j = 0, drifts with the prior).
   - *Done (Proposition D‴.24(i)):* the same spectrum sets the finite-size floor of a fixed prior, each mode weighted by the band-pass
     h_jσ₀²/(h_j + σ₀)² ≤ σ₀/4, centred at h = σ₀ — the floor sits in the modes that carry the gauge's displacement, re-weighted toward the
     level of the leak.
   - *Open:* the continuum form of (PE); the trajectory kernel K_T(ϑ,ϑ') = ∫Φ(x(s);ϑ)Φ(x(s);ϑ')ds; identifiability of
     measures up to symmetry; the relation to [16]. Remark E.1 shows that the PE level alone does not predict the rate.
     For Lemmas E.0′ and E.0″: lowering the threshold k₀ (a third-order slaving trades ρ_DD for the peak factor of Ψ⃛ one
     order higher); the choice of the window T as an optimisation; the nonlinear Theorem E with the constants of E.0″
     (Corollary E.0′.2's steps carry over verbatim, with λ_T in place of λ_min(M̄)).

8. **Transport-regularised adaptive laws.**
   - Wasserstein proximal (JKO) discrete-time adaptive laws.
   - *Done (Theorem C′):* the OT leak toward priors with a different number of atoms, using optimal plans with
     barycentric projection; it is Lloyd's algorithm in continuous time, and its bias is the centroid defect of the ideal.
   - *Open (in part):* its mean-field form — done for an atomic prior (the Laguerre-cell leak, Proposition D‴.12), open for a general
     prior ν and a general μ (where the barycentric projection is the superdifferential of ½W²(·, ν); AGS Thm 10.4.12) — and adaptive
     priors (ν itself learned: its atoms, Theorem D‴.23; its masses, Proposition D‴.24).
   - *Done (Proposition D‴.9):* the nearest-atom leak toward a prior measure inherits Theorem D‴, Proposition D‴.5 and Corollary D‴.7
     inside the Voronoi cells; the prior is a design input of the window (maximal margin κ = σ₀ at a matching prior, the window
     closing as the prior is misplaced) and of the floors (which grow with the stationary amplitudes).
   - *Done (Proposition D‴.12):* the barycentric OT leak's mass constraint against the weight fluctuation. The nearest-atom leak conserves
     the cluster masses, and the finite-m floor of Proposition D‴.5 is the memory of the initial cloud that this conservation carries; the
     OT leak constrains the masses to the prior's, so on the sizes commensurate with the masses the finite-m stationary measure is the
     mean-field one exactly (the replicated weighted configuration is a critical point for every such m, by total unimodularity and
     c-cyclical monotonicity), and otherwise differs from it by at most K₀ − 1 split particles; for an atomic prior the mean-field form of
     the barycentric leak is the Laguerre-cell leak. *Open:* a coupling theorem for the OT dynamics from an i.i.d. start (Theorem D‴'s synchronous coupling needs a
     balanced sample, of probability ≈ 0.83/m, and otherwise keeps the transit's mismatch for ever; the convergence of every trajectory to a
     Clarke-critical point is Proposition D‴.17, without a rate; Theorem D‴.21 gives the rate in an m-independent energy basin, with the plan as
     the test direction, and formally (Remark D‴.22) propagation of chaos with no floor there; smallness in m during the transit from an
     i.i.d. start remains open); a penalised
     (unbalanced) mass constraint — the entropic one is Proposition D‴.20, which equalises the masses only down to a critical excess and
     conserves smaller ones, so that its floor returns at large m (Proposition D‴.5's law at the entropic stationary point); the
     mean-field form for a general prior ν (the barycentric projection as the superdifferential of ½W²(·, ν), AGS Thm 10.4.12) and adaptive
     priors with learned masses beyond the cell masses (Proposition D‴.24).
   - *Done (Theorem D‴.23):* adaptive priors with the masses fixed. The prior's atoms move toward the class barycentres of the plan; the
     joint flow is the gradient flow of L̄ + (σ₀/2)W²(μ_Θ, ν_Z), its replicated stationary points are the unregularised K₀-atom critical
     points (no leak bias; the non-replicated ones carry the K₀-cluster bias), and at an exact representation every collective mode relaxes
     at λ_j^− < min{h_j, η}, for every commensurate m; the
     fast limit is a capacity-constrained K₀-means regulariser. The transit drags the prior and its weak modes remember the drag, so the
     adaptation starts from the fixed-prior stationary point (a gate). *Open:* the real loop with an adaptive prior over the adaptation's own
     time scale beyond the linearisation (Corollary D‴.26); the plateaus of the unregularised loss that the fast limit exposes; the
     constants of (v).
   - *Done (Corollary D‴.26):* perfect tracking with an adaptive prior is an exponentially stable periodic solution of the real loop
     (locally, for each m) at every representable sample: its Floquet exponents are σ₀/k_e exactly on the separations and the weighted
     K₀-atom loop's on the collapsed subspace, for every m, and the latter are (1/k_e)[λ_j^±(w) + c_j^±(w)/k_e² + O(k_e⁻⁴)] — the odd orders
     vanish (the time-reversal parity of the averaged terms), c_j explicit — the error filter slowing, the slaving speeding up, its sign
     not fixed. *Open:* its basin uniformly in m (Theorem D‴.18's contraction for the joint system).
   - *Done (Lemma D‴.26a):* the slow Floquet expansion proved for the linearised loop. Chang's periodic decoupling makes the slaving exact,
     so the slow multipliers are exactly the slow block's, with an expansion to O(k_e⁻⁴) whose remainder is bounded; defining the slow block
     at negative ε by the same contraction, the time-reversal symmetry becomes the identity eig Ā(−ε) = eig Ā(ε) — the odd orders vanish
     exactly; and Chang's second transformation decouples the slow coordinates' drift exactly, their noise εbν with b the high-passed
     −(K̇ + aK)ᵀ at leading order and its next term explicit, which proves Proposition D‴.30(ii)'s floor for the linearised loop under white
     noise as a period average (pointwise the fast jitter oscillates by ±1% on the instance). *Open:* the constants explicit; the evenness
     of the noise floor without the slope (observed to k_e = 800, not proved); the real loop's nonlinear noise response.
   - *Done (Proposition D‴.28):* the adaptive prior's floor removal is a property of the amplitude–mass gauge, not of exact
     representability: at any nondegenerate critical point of the unregularised loss the gauge transform is a critical point at the
     sampled weights, with the same function, so at a strict local minimum the function-level floor vanishes at the samples of the gauge
     set W_g with positive margins, whatever the approximation error (W_g shrinking as the error grows); the fixed prior's floor is the
     bias formula with the full Hessian; for a scalar amplitude the residual curvature scales exactly by q_k/w_k under the gauge, so the
     adaptive prior's splitting margins are σ₀ − (q_k/w_k)λ⁺_k. *Open:* critical points at infinity (the best approximant on the
     saturation boundary, as for a sharp friction step, a cogging function that two kernels fit as a dipole, or Proposition D‴.29's
     servo globally), where only the projected prior keeps the parameters bounded; vector amplitudes, where the margins have no closed
     form.
   - *Done (Proposition D‴.36):* the replicator's gauge under dither. At an exact representation the fits through the configuration
     form a manifold parametrised by the masses, on which f̂ does not move; dithering the positions makes the masses diffuse along it
     (Katzenberger's limit, the noise projected in the flow's own metric), measured at the predicted rate, and their amplitudes run into
     the saturation (|a|/a_max = 0.9994) at no cost in fit — the replicator removed the masses' obstruction (Proposition D‴.34) and with
     it their stiffness. An entropic leak on the masses restores it without biasing the representable target: the spread is the
     Ornstein–Uhlenbeck one at ε = 10τ (to 3%; at ε = τ four seeds do not
     resolve it). *Open:* the limit's drift term (visible in the means, a selection of the
     masses by the dither); the gauge with several particles per class; the real loop's dither.
   - *Done (Proposition D‴.35):* how a class splits. The class's separations along the splitting direction carry the standard
     representation of S_n, and the reduced flow has the quadratic equivariant with the residual's third derivative as its coefficient
     (the cohesion, even, enters at third order): the bifurcation is transcritical, every small split state a saddle, the split a jump.
     The particles keep their order, and once its excursion passes 2μ/|c₃| the extreme particle runs away on the side the quadratic term
     picks, alone or with a near-tied neighbour — the quadratic equation names the first excursion's particle in 12 of 12 runs and its
     time to 1.1%, where the linear law names it in 6 and misses the time by 22% —; after a departure the stiffness's drop and the
     curvature's gauge rise favour more, the residual's change stops them. A mass-free leak has thresholds a factor w_k lower and only
     delays the cascade: continued to σ₀/256 it ends at 8–9 classes and the fit error 0.0028, as the cohesion leak
     does; its 4 classes at σ₀/64 were the schedule's. *Open:* steering the escape (the sign of c₃ against the residual's
     location); a schedule-free count; the transit beyond one split.
   - *Done (Proposition D‴.34):* mass transfer for the cohesion leak. The particles carry weights adapted by the replicator law (in
     real time an exponentiated-gradient law with the relative entropy of the weights as its Lyapunov function, the negative-entropy
     instance of a mirror-descent adaptation law); with the transport law the pair is a gradient flow (Euclidean on the positions,
     Fisher–Rao on the weights: gradient descent on (√β, Θ)), Σβ and the saturation's bound kept. At clustered configurations under the
     finite range its stationary points are D‴.32's; under the Gaussian kernel it selects the masses and moves the bias. The saturated
     transients that froze the transit's masses are left: at fixed positions bounded trajectories reach the exact fit iff the target is
     representable (Σ_k|A_k|/a_max < 1), and from broad starts the five-bump and the three-heavy-bump targets of D‴.32 (H), (K) are found
     with their kernels. It does not separate bumps closer than the broad scale, nor act on the shedding (the fitness vanishes at a
     critical point). *Open:* boundedness at fixed positions; the transit with the replicator (which classes form; the shedding of D‴.32
     (G)); the neutral directions under noise (a common noise random-walks the gauge; the redistributions move only when the positions are
     dithered); the real loop with the transport.
   - *Done (Corollary D‴.33):* the desired regressor keeps the replicated configuration attracting under measurement noise. With the
     regressor at the reference the regressor's own term leaves the separations' second-order drift; what remains — the mean error's
     term and the cross term — is O(1/k_e) like the separations' stiffness, so the noise threshold is bounded below uniformly in k_e
     (constant when the noise's bandwidth scales with the gain; at a fixed bandwidth DR's sign reverses at k_e ≈ B). Under DR near the
     lift the compensation's rectification of the noise stiffens the separations of separated alternating-sign kernels (the instance;
     same-sign, wide or edge kernels soften); along the trajectory the drift turns indefinite as the mean error relaxes with the
     collective modes, and at the stationary point the parameters' part of the force dominates the mean error. It makes the desired
     regressor the partner of the cohesion leak, whose separations are the softest (σ₀w_k). *Open:* the long-time threshold (the drift
     indefinite along the trajectory and at the stationary point, formally at noise levels outside the expansion); the higher orders;
     the real loop beyond 300 time units; references in ℝⁿ.
   - *Done (Proposition D‴.32):* the cohesion leak, a prior-free leak on the unexcited subspace. The particles attract each other
     toward the modes of their own distribution (the mean shift); the averaged loop is the gradient flow of the fit plus an aggregation
     energy; with a finite-range kernel, at a clustered configuration whose atoms are R apart the leak vanishes, so every critical point
     of the unregularised loss is stationary, for every count and every set of masses, and it acts exactly on the unexcited subspace (the
     separations at σ₀w_k, the class means free); a class can split iff σ₀w_k² < q_kλ⁺_k, so σ₀ sets the resolution and an exact
     representation is never over-split, and a split sheds single particles. A Gaussian kernel pays an exponentially small bias, and
     twins with no third class within reach contract at an explicit rate (checked to 0.04%) — a third class within reach can win. The
     number of kernels is not selected by the leak: an annealed Gaussian scale found the instance's three kernels from every broad start,
     a fixed scale did not, and on other targets the count failed where the transit's masses were not representable or the bumps were
     closer than the broad scale. *Open:* the transit (which classes form, with which masses), and with it any selection; the shedding
     below the resolution's threshold (Part G); the scale R chosen from the data (Parts J, K); the a-priori bound (the leak is internal:
     a projection). The masses' obstruction is removed by Proposition D‴.34.
   - *Done (Proposition D‴.31):* the adaptive prior is a washout leak, and the remedy for the errors in variables. The prior lags the
     class barycentre by its exponentially forgotten increments, so the leak acts in full on the separations and on the barycentre through
     the high-pass iω/(iω + η) — the bias removal and the lost robustness to persistent forces in one transfer function; along every affine
     time-invariant prior the stationary response R and the bias B satisfy RH + B = I, so no prior gives both; the response to a force is
     at most the norm of its smallest equivalent target perturbation, under every leak; and the desired regressor (the adaptive law's
     regressor at the reference) keeps the lift, the linearisation and the first-order floor while turning the noise's own persistent force
     into a least-squares one, bounded uniformly in k_e and in the spectrum — the parameters' fluctuation adding a force that is not the
     learning force of a mean tracking error (the classical second-order bias of a constant-gain stochastic approximation): small in
     absolute terms on the instance (≤ 3.9·10⁻⁴ per (Eν²)²), −8.7% to +3.3% of DR's and DR*'s biases and +56% of DR's at k_e = B = 100,
     its own response falling like k_e⁻² on a one-dimensional reference, where the law's leading term averages out — which on a curve in
     ℝⁿ it need not. *Open:* the second-order computation made rigorous, the parameters' part in closed form, and in ℝⁿ; the
     stationary biases in the real loop (the weak modes need ~10⁵ time units); a sensor model with quantisation and delay.
   - *Done (Proposition D‴.30):* measurement noise at perfect tracking. The noise is coherent — a replicated start stays replicated, so
     its floor does not average out with m (invariance, not attraction: near the lift the separations grow above a noise level, which
     drifts along the trajectory with the developing errors-in-variables deviation); at first order the
     loop high-passes the measured error's noise and the function-level floor is (N₀/2k_e)[Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) + tr H²], the
     adaptive prior costing at most the filter factor; at second order the noise inside the regressor exerts a mean force (the errors in
     variables) whose stationary response is (H + σ₀I)⁻¹ under the fixed prior and H⁻¹ under the adaptive prior for every η > 0 — the
     robustness half of the leak, lost in the weakly excited modes. *Open:* the real loop's noise response beyond the linearisation
     (the linearised loop under white noise: Lemma D‴.26a); a projected or partially fixed prior on the weak modes; the optimal k_e and η
     for a given noise level.
   - *Done (Proposition D‴.29):* second-order plants. With the sliding variable s = ė + λe the error channel is Theorem J's with s in
     place of e, so the averaged slow loop — and with it every averaged-loop statement of §5d — holds (formally) along the regressor's
     orbit, the position error being the residual filtered by (d/dt + λ)⁻¹(d/dt + k_s)⁻¹. On a cogging-torque servo whose torque four
     kernels cannot represent — the loss's infimum lies at the saturation, and its local minimum ϑ^b(q) is Proposition D‴.28's regime — the
     adaptive prior removes the fixed prior's floor and, from an engineer's misplaced model, most of its bias, at the rate of the weakest
     excitation; the real servo follows the exactly filtered prediction to 0.3%. Noise-free, no leak at all does as well: the floors and
     biases the leaks remove there are those a fixed leak creates. *Open:* the real-loop statements (Corollaries D‴.14, D‴.26) for
     second-order plants; the servo with measurement noise and disturbances, where the leaks' robustness is the point, against the
     classical baselines and a Fourier regressor of the known harmonics; velocity-dependent friction with the projected prior in the real
     loop.
   - *Done (Proposition D‴.27):* the a-priori bound restored. Confining the adaptive prior's atoms to a compact convex set keeps every
     particle within B/σ₀ of it (eventually) whatever the target and the energy non-increasing (at η = ∞ in its selected form, the plan
     being Lloyd-consistent, not jointly optimal); inside the set the bias removal of Theorem D‴.23 is unchanged, and where the limit has an
     atom on the boundary a bounded bias remains — the box trading the fit error against the parameters' bound, as long as it contains
     the representation; at η = ∞ on collapsed classes it is a switching σ-modification toward the set. *Open:* the convergence (at η = ∞ the
     Kurdyka–Łojasiewicz argument on the selected pieces, at finite η with the indicator), the choice of the set as a design input (how
     large a box keeps every representable target inside), and a sharper bound than the worst-case force B.
   - *Done (Proposition D‴.24):* the masses minimised out of the leak's own energy are the cell masses (min_q W²(μ, ν_{Z,q}) is the
     nearest-atom distortion), so the nearest-atom leak toward an adaptive prior is the OT leak with free masses — and the function-level
     floor does not return at the representable samples (those outside W, of probability exponentially small in m but 0.36 at m = 15 on the
     instance, keep an error): it is the weight-derivative of the leak's bias at the sampled weights, which the adaptive atoms remove. What
     remains is the coupled distance of the amplitude–mass gauge, conditional on representability, which the mass constraint replaces by
     the transit's mismatch, and at the measure level the cross-class transport common to both priors; the adaptive prior also removes the
     fixed prior's splitting instability. With an adaptive prior the real loop has perfect tracking as an invariant set at every
     representable sample. *Open:* the gauge as a design input (the masses set the excitation spectrum
     H_L(q)); the floor law and its removal for a misplaced prior beyond first order.
   - Transfer: transporting a learned measure to a new plant.

9. **Quantisation and the curse of dimensionality (Corollary B.6).**
   - *Done (Theorem C′):* an adaptive law whose leak tracks the equal-mass centroidal quantisers of a prior measure
     (centroidal Voronoi tessellations [S20]; coverage control [S18], [89]), and Theorem D′'s observation that the state and the learned
     function converge at the rate m^{−1/2}, whose exponent is dimension-free: the curse of dimensionality of Corollary B.6 lives in the
     W-distance of the parameter measures, which the loop never evaluates.
   - *Open:* second-order (centroid) quantisation rates for the ideal measure on the lifted space; approximation rates
     that do not go through W (Maurey-type).

10. **Symmetric systems beyond kernels.** The OT quotient applies verbatim to swarms of identical agents, mixtures and
    multi-head models ([88] and [89] are the control-side precedents).

11. **Deep approximators.** Hierarchical (multi-layer) mean-field measures, with a layer-wise OT Lyapunov function.

---

## 9. References

See `LITERATURE_TRANSPORT.md`: [1]–[89] are in `review/LITERATURE_refs.md`; [S1]–[S39] are in the audit file.

Key mathematical references used in proofs and remarks:

- **Transport and gradient flows:** Ambrosio–Gigli–Savaré (2008); Ambrosio–Savaré (2007); Golse (2013).
- **Control and dynamical systems:** Khalil (2002), Thms 4.11, 4.13, 10.4, 11.1; Kokotović–Khalil–O'Reilly (1986);
  Morgan–Narendra (1977), Thm 1; Pomet–Praly (1992); Sastry–Bodson (1989); Sanders–Verhulst–Murdock (2007).
- **Mean-field Langevin dynamics and functional inequalities:** Chizat (2022, TMLR), Thm 3.2, Cor. 3.3, Lemmas 2.2 and
  3.4, Thm 4.1; Nitanda–Wu–Suzuki (2022); Chizat–Bach (2018); Holley–Stroock (1987); Bakry–Émery (1985);
  Otto–Villani (2000).
- **Quantisation:** Graf–Luschgy (2000); Chevallier (2018); Fournier–Guillin (2015).
- **Semi-gradient learning, performativity and non-reciprocity:** Tsitsiklis–Van Roy (1997);
  Zhang–Cai–Yang–Chen–Wang (2020); Perdomo–Zrnic–Mendler-Dünner–Hardt (2020); Fruchart–Hanai–Littlewood–Vitelli (2021).
- **Log-Sobolev inequalities on convex bodies (Remark H.3):** F.-Y. Wang (1997); Frieze–Kannan (1999).
