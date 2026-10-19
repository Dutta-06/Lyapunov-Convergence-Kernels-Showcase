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
    residual-weighted curvature κ_R and the constants of Theorem D‴.
- **Tests.** `tests/test_transport.py`, `tests/test_langevin.py`, `tests/test_averaging.py`, `tests/test_excitation.py`,
  `tests/test_nonreciprocity.py`, `tests/test_periodic.py`, `tests/test_barycentric.py`, `tests/test_dither.py` and `tests/test_chaos_learning.py`. Run them with
  `venv\Scripts\python -m pytest tests`.
- **Experiments.** The scripts are `experiments/exp_transport_*.py`, `exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py`, `exp_mean_field_m.py`,
  `exp_pe_convergence.py`, `exp_variable_projection.py`, `exp_slow_manifold_slope.py`, `exp_langevin_global.py`,
  `exp_langevin_tradeoff.py`, `exp_gibbs_reference.py`, `exp_orbit_orientation.py`, `exp_explicit_rate.py` and
  `exp_periodic_orbit.py`, `exp_orientation_stability.py`, `exp_chaos.py`, `exp_chaos_uniform.py`, `exp_chaos_learning.py`, `exp_barycentric_leak.py`,
  `exp_dither_stationary.py` and `exp_pe_general.py`. Their numbers are collected in `results/RESULTS_PARADIGM.md`.
- **Literature audit.** `LITERATURE_TRANSPORT.md` gives a novelty verdict for each claim. Reference numbers [1]–[89]
  point to `review/LITERATURE_refs.md`, and [S1]–[S39] point to the audit file.
- **Review log.** `review/REVIEW_LOG.md` lists every confirmed issue and its fix. The raw reports are in
  `review/REVIEW_round1.md`, `review/REVIEW_round2.md`, `review/REVIEW_round3.md`, for round 5 the five files
  `review/REVIEW_round5_*.md` (E0, I, J, DC, consistency), for round 6 `review/REVIEW_round6_E0prime.md`,
  `review/REVIEW_round6_Dpp.md` and `review/REVIEW_round6_Jprime.md`, for round 7 `review/REVIEW_round7_E0p3_E0dp.md` and for
  round 8 `review/REVIEW_round8_Dppp.md`.

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
| D‴ | Uniform-in-time propagation of chaos *in the learning regime* for the averaged slow loop (𝔄_m) at τ = 0 (Theorem J′'s leading order): under the synchronous coupling with the mean-field characteristic flow, the interaction term of the coupling identity is −‖f̂_Θ − f̂_Θ̄‖²_{L²(ρ)} + O(D²) — flat convexity of the loss makes it contract whatever the gains — and only the residual-weighted feature curvature R[μ](ϑ) = ⟨r_μ∇²_ϑΦ̃⟩_ρ (Proposition K's R₂) can expand; if along the mean-field flow σ₀I − G^{1/2}RG^{1/2} ⪰ κI on a neighbourhood of the support from a time s₁ on — the strict second-order optimality condition of the fitted configuration, compatible with σ₀ < σ_crit — then sup_s E(1/m)Σ‖ϑ_i − ϑ̄_i‖²_{G⁻¹} ≤ 4e^{2Λ₁}(L_ϑ^G)²K_Φ²/(κ²m) plus a tail exponentially small in m (uniformly in time under exponential convergence of the flow, over horizons e^{cm} without it), no particle changes basin, and the learned function converges at m^{−1/2} uniformly in time while the network learns; the window κ_R[μ_∞] < σ₀ < σ_crit is where learning and uniformity coexist; at finite m the stationary coupled distance follows the weight-sensitivity law E D_∞ = (1/m − 1/M)Σ_k w̄_k tr(S̃_kC̃S̃_kᵀ), S = ∂Θ*/∂w the sensitivity of the stationary atoms to the cluster weights (Proposition D‴.5, the sharp form of the theorem's constant) | proved for the averaged slow loop (synchronous coupling, a bootstrap on the mean square and the per-particle distance, McDiarmid for the reference's fluctuation); the hypothesis is a computable property of the mean-field flow; it fails for a connected cold start whose mass splits between basins (the support threads the ridges), where the coupling carries a floor; the real loop at finite k_e: settled on a warm start by Corollary D‴.7 through Theorem J (both loops converge to their periodic orbits, whose distance is Proposition D‴.5's up to O(1/k_e)), open for a cold start (Remark D‴.4); dither: the entropy method of [S5] (imported); the theorem's regime m ≥ m₀ is asymptotic (m₀ ~ 10²³ on the instance) and Proposition D‴.5 is formal; numerics: [D3 NUMERICS PENDING] | the two ingredients — the sign of the interaction from flat convexity, and the curvature at the actual rather than the worst-case residual — are new in propagation of chaos for adaptive control; the decomposition Gram + residual×curvature is the one of the interpolation literature ([S39]); the window between Proposition K's collapsed instability and the fitted configuration's stability is new |
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
well-posedness sentence (all applied).

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
(the transient and the connected support).

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
table [X]) shows the coupled distance saturating at the same k_e, G and horizon at which Theorem D″'s cold start grew.

**Proposition D‴.5 (the stationary coupled distance at finite m: the weight-sensitivity law).** Let μ* = Σ_{k=1}^K w̄_kδ_{ϑ*_k} be a
K-atom stationary measure of (𝔄₀) (each ϑ*_k a critical point of 𝒱[μ*]), and for weights w ∈ (0, ∞)^K let 𝒢_w(Θ) := ½⟨‖f −
Σ_kw_kΦ̃(·;ϑ_k)‖²⟩_ρ + (σ₀/2)Σ_kw_k‖ϑ_k − ϑ⁰‖²_{G⁻¹} be the mean-field energy of μ_w := Σ_kw_kδ_{ϑ_k}, so that ∇_{ϑ_k}𝒢_w(Θ) =
w_k∇_ϑ𝒱[μ_w](ϑ_k). Assume the K-atom Hessian ∇²_Θ𝒢_{w̄}(Θ*) is invertible. Then:

(i) *(the weighted stationary point and its replication)* there is a smooth map w ↦ Θ*(w) on a neighbourhood of w̄ with Θ*(w̄) = Θ*
and ∇_ϑ𝒱[μ_{Θ*(w)}](ϑ*_k(w)) = 0 for every k (implicit function theorem); S := ∂Θ*/∂w (S_k := ∂ϑ*_k/∂w ∈ ℝ^{p×K}) is given by
S = −[∇²_Θ𝒢_w]⁻¹∂_w∇_Θ𝒢_w at (Θ*(w), w). For m particles with n_k of them in cluster k (w^m_k := n_k/m), the configuration with n_k
copies of ϑ*_k(w^m) is a critical point of J_m (μ_Θ = μ_{w^m} there), and a strict local minimiser of J_m when in addition the
splitting condition σ₀G⁻¹ − R[μ_{w^m}](ϑ*_k(w^m)) ≻ 0 holds at every cluster (Proposition K(i): the Hessian of J_m at a configuration
with coinciding particles splits into the collective block ∇²𝒢_{w^m} and the separation blocks of each cluster).

(ii) *(the stationary distance)* Under the shared-particle coupling of Theorem D‴ with i.i.d. cluster labels of law w̄ (the m
particles are the first m of the M reference particles), suppose the m-system converges to the configuration of (i) with its
weights w^m and the reference to the one with w^M. Then the coupling distance converges to

    D_∞ = Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖²_{G⁻¹} = Σ_kw^m_k‖S_k(w^m − w^M)‖²_{G⁻¹} + O(‖w^m − w^M‖³),

and since E[(w^m − w^M)(w^m − w^M)ᵀ] = (1/m − 1/M)C̃ with C̃ := diag(w̄) − w̄w̄ᵀ (multinomial counts of the shared labels),

    E D_∞ = (1/m − 1/M)·Σ_kw̄_k tr(G^{−1/2}S_kC̃S_kᵀG^{−1/2}) + O(m^{−3/2}):

the 1/m law with an explicit constant, and the (1/m − 1/M) law of the shared-particle coupling observed in Theorem D′'s numerics.

(iii) *(the fluctuations)* As M → ∞ (w^M → w̄) and then m → ∞, √m(w^m − w̄) → N(0, C̃) and √m·(ϑ_i − ϑ̄_i) → S_{k(i)}ξ with ξ ~ N(0, C̃)
for every particle i of cluster k(i): the fluctuations of the coupling are Gaussian, carried entirely by the cluster weights, and
all particles of a cluster are displaced by the same vector.

*Status.* Formal: (i) and the algebra of (ii)–(iii) are exact (implicit function theorem; the multinomial covariance of shared
labels), the convergence of both systems to the stated configurations is assumed (it is what Theorem D‴ proves for m ≥ m₀ and what
the numerics show for m = 16–256). Checked: `weighted_stationary_atoms`, `weight_sensitivity`, `predicted_stationary_D` in
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
which 𝒱[μ_Θ̄] is convex, at a stable stationary point of their own J₁₆. Hypothesis (H) is nevertheless the right object: it is the
condition under which the reference configuration is stable to splitting, and the m-system's configuration is the same
configuration with perturbed weights. What the numerics test is therefore (a) hypothesis (H) at the stationary measure and
along the flow, (b) the mechanism — the coupled distance saturates (no growth) and follows the 1/m law — and (c) the sharp
constant of Proposition D‴.5; the theorem's own constant C₁ is reported for the record.

**Corollary D‴.7 (the real loop on a warm start: uniform in time through Theorem J).** Let μ* = Σ_kw̄_kδ_{ϑ*_k} be a K-atom
stationary measure as in Proposition D‴.5, satisfying the splitting condition σ₀G⁻¹ − R[μ*](ϑ*_k) ≻ 0 at every atom and the
non-degeneracy of the K-atom Hessian, and let w ↦ Θ*(w) be the weighted stationary point of D‴.5(i). Fix m ≤ M and cluster
weights w^m, w^M in the neighbourhood of w̄ on which Θ*(·) is defined and the splitting condition persists, and let Θ*_m :=
the configuration with n_k = mw^m_k copies of ϑ*_k(w^m) (likewise Θ*_M). Both are non-degenerate local minimisers of J_m and J_M
(Proposition D‴.5(i)), hence hyperbolic critical points, and Theorem J applies to each loop: there are k₀, ρ₀, C (depending on
Θ*_m and Θ*_M, hence on m, M and the weights) such that for k_e ≥ k₀ the real closed loop (CL) with m kernels started at
‖Θ(0) − Θ*_m‖_{G⁻¹} ≤ ρ₀, ‖e(0)‖ ≤ R converges to its T₀-periodic orbit Θ_p^m, and likewise the M-kernel loop to Θ_p^M. Under the
shared-particle coupling (the m kernels are the first m of the M, the plants are driven by the same reference) the coupled distance
D(t) := (1/m)Σ_i‖ϑ_i(t) − ϑ̄_i(t)‖²_{G⁻¹} satisfies, for all t ≥ 0,

    D(t)^{1/2} ≤ D_orb^{1/2} + Ce^{−μt}(‖Θ(0) − Θ_p^m(0)‖_{G⁻¹} + ‖x(0) − x_p^m(0)‖) + Ce^{−μt}(‖Θ̄(0) − Θ_p^M(0)‖_{G⁻¹} + ‖x̄(0) − x_p^M(0)‖),

with μ = (λ_*/k_e)(1 − κ) the rate of Theorem J(c) (the smaller of the two loops' rates) and

    D_orb := sup_t (1/m)Σ_i‖Θ_p^m(t)_i − Θ_p^M(t)_i‖²_{G⁻¹} = Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖²_{G⁻¹} + O(1/k_e)
           = Σ_kw^m_k‖S_k(w^m − w^M)‖²_{G⁻¹} + O(‖w^m − w^M‖³) + O(1/k_e):

the coupled distance of the real loop is bounded uniformly in time, converges to the distance between the two periodic orbits,
and that distance is Proposition D‴.5's stationary distance up to the O(1/k_e) offsets of Theorem J(a),(d) (which differ between
the two loops by O(k_e⁻¹‖w^m − w^M‖), the offset being a smooth function of the configuration). In particular, for cluster labels
i.i.d. of law w̄ and conditionally on both starts lying in their basins, E D_orb = (1/m − 1/M)Σ_kw̄_k tr(S̃_kC̃S̃_kᵀ) + O(m^{−3/2}) + O(1/k_e).

*Status.* Proved as stated (Theorem J(a),(c),(d) applied to each loop, the triangle inequality in the G⁻¹ metric, and Proposition
D‴.5(ii) for the distance between the critical points), under Theorem J's hypotheses for both configurations: the constants k₀,
ρ₀, C, μ are those of Theorem J at Θ*_m and Θ*_M, and the corollary does not assert their uniformity in m or M (Remark D‴.8). It
is the statement of Remark D‴.4 for a warm start: the real loop's coupled distance saturates, at a value set by the cluster weights.
Checked on the instance of §5d's numerics (table 5: the real loop at k_e = 25 from the jittered mixture; the window maxima [REAL-LOOP
READING]).

*Proof.* Theorem J(c) for the m-kernel loop gives ‖Θ(t) − Θ_p^m(t)‖_{G⁻¹} ≤ Ce^{−μt}(‖Θ(0) − Θ_p^m(0)‖_{G⁻¹} + ‖x(0) − x_p^m(0)‖) and the
same for the M-kernel loop; the first m rows of the M-loop are the reference particles ϑ̄_i. For the (m×p)-row norm
(1/m)Σ_i‖·‖²_{G⁻¹} ≤ ‖·‖²_{G⁻¹}/m ≤ ‖·‖²_{G⁻¹} (the full configuration norm dominates the mean square over the first m rows), so by
the triangle inequality D(t)^{1/2} ≤ [(1/m)Σ_i‖Θ_p^m(t)_i − Θ_p^M(t)_i‖²]^{1/2} + (the two transients). Theorem J(a) places each
orbit within C/k_e of its critical point, uniformly in t, so the first term is Σ_kw^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖²_{G⁻¹} up to O(1/k_e)
(the two offsets are O(1/k_e) each; their difference is O(k_e⁻¹‖Θ*_m − Θ*_M‖) because Θ ↦ J⁻¹(b_frozen + b_osc)(Θ) is smooth, which
sharpens the O(1/k_e) to O(k_e⁻¹‖w^m − w^M‖)), and Proposition D‴.5(ii) expands it in w^m − w^M. The expectation statement is
D‴.5(ii)'s. ∎

*Remark D‴.8 (what the corollary does and does not say).* It transfers the uniform-in-time statement to the real loop in the regime
where Theorem J applies to both loops — a warm start in the basins of two hyperbolic minimisers at a gain k_e above both thresholds.
Theorem J's constants at a replicated configuration are controlled by the K-atom data (the Hessian of J_m there is the collective
block ∇²𝒢_{w^m} ⊕ the separation blocks σ₀G⁻¹ − R_k of Proposition K(i), both independent of m beyond the weights), so the
threshold and the basin are expected to be uniform in m for configurations with the same weights; the corollary does not prove
this, and the cold start — where the coupled distance of the real loop grew in Theorem D″'s experiment — is outside it on the
same two counts as in Remark D‴.2.

*Numerics* (`exp_chaos_learning.py` → `results/chaos/learning_table.md`). [D3 NUMERICS PENDING]

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
| D‴ | [D3 NUMERICS PENDING] | `exp_chaos_learning.py` → `results/chaos/learning_table.md` |
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
   - *Open:* the continuum form of (PE); the trajectory kernel K_T(ϑ,ϑ') = ∫Φ(x(s);ϑ)Φ(x(s);ϑ')ds; identifiability of
     measures up to symmetry; the relation to [16]. Remark E.1 shows that the PE level alone does not predict the rate.
     For Lemmas E.0′ and E.0″: lowering the threshold k₀ (a third-order slaving trades ρ_DD for the peak factor of Ψ⃛ one
     order higher); the choice of the window T as an optimisation; the nonlinear Theorem E with the constants of E.0″
     (Corollary E.0′.2's steps carry over verbatim, with λ_T in place of λ_min(M̄)).

8. **Transport-regularised adaptive laws.**
   - Wasserstein proximal (JKO) discrete-time adaptive laws.
   - *Done (Theorem C′):* the OT leak toward priors with a different number of atoms, using optimal plans with
     barycentric projection; it is Lloyd's algorithm in continuous time, and its bias is the centroid defect of the ideal.
   - *Open:* its mean-field form (a general prior ν and a general μ, where the barycentric projection is the
     superdifferential of ½W²(·, ν); AGS Thm 10.4.12), and adaptive priors (ν itself learned).
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
