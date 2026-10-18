# Transport adaptive control: adaptive approximators as evolving measures

**Version 3** (after adversarial review rounds 1–3 and the literature audit). Version 2 is archived as
`review/PARADIGM_v2.md`; the review reports cite its line numbers.

**Companion files.**

- **Code.**
  - `src/core/transport.py` implements laws (L) and (L_σ) and the oracle Lyapunov report.
  - `src/core/langevin.py` implements Theorem H.
  - `src/core/averaging.py` implements the frozen-measure averaging of Proposition H.7.
- **Tests.** `tests/test_transport.py`, `tests/test_langevin.py` and `tests/test_averaging.py`. Run them with `venv\Scripts\python -m pytest tests`.
- **Experiments.** The scripts are `experiments/exp_transport_*.py`, `exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py`, `exp_mean_field_m.py`,
  `exp_pe_convergence.py`, `exp_variable_projection.py`, `exp_slow_manifold_slope.py`, `exp_langevin_global.py`,
  `exp_langevin_tradeoff.py`, `exp_gibbs_reference.py` and `exp_orbit_orientation.py`. Their numbers are collected in `results/RESULTS_PARADIGM.md`.
- **Literature audit.** `LITERATURE_TRANSPORT.md` gives a novelty verdict for each claim. Reference numbers [1]–[89]
  point to `review/LITERATURE_refs.md`, and [S1]–[S14] point to the audit file.
- **Review log.** `review/REVIEW_LOG.md` lists every confirmed issue and its fix. The raw reports are in
  `review/REVIEW_round1.md`, `review/REVIEW_round2.md` and `review/REVIEW_round3.md`.

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
    present for every n, and it is non-reciprocal in the checked instances (n = 1 and n = 2).
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
| D | Local well-posedness of the mean-field loop; a Dobrushin estimate; a tube bound and global existence under a self-consistent growth condition; the mean-field limit on the existence interval | proved; consistency in m checked numerically (the constant is not checked) | technique anticipated [77]–[79] |
| E | Label-free local exponential convergence under tangential PE, in a uniform tube, at fixed m; exponential merging of solutions | proved, conditional on [S8, Thm 1] as stated in its abstract (its hypotheses were not checked against the full text); linearisation checked; convergence observed | partially anticipated [49], [54] |
| F | Along (L), the Wasserstein length of learning is ≤ λ_G^{1/2}L₁∫‖e‖; leak versions | proved and checked, including the leak clause, a tight leak instance with e ≡ 0 and a tight velocity instance | elementary |
| G | Averaged joint descent on the occupation loss. The slow geometry follows the gradient flow of the profile loss in the metric induced on the slow manifold, which is plain variable projection only when that manifold is flat (S ≪ 1). | formal; checked directly, including a counterexample to plain VarPro at S = O(1) | a fixed-m, closed-loop analogue of the two-timescale/VarPro training of [S1], [S2] |
| H | The averaged mean-field loop with saturated features, σ-modification and dither is mean-field Langevin dynamics. It converges globally and exponentially to the regularised optimum, with the log-Sobolev constant given in control quantities. | proved for (𝔄); the reduction to (𝔄) is formal (Remark H.0); real loop: Conjecture H.6 | a bridge to [S3], [S4]; new in adaptive control |
| H.7 | At a frozen measure, the exactly averaged loop is (𝔄) plus an explicit O(1/k_e) correction: a semi-gradient term, whose performative defect is present for every n and non-reciprocal in the checked instances (n = 1, 2), and an orbit-area circulation kernel (n ≥ 2). Along the real loop, a second-order averaging term is added (Remark H.7a). | Frozen-measure part proved. Stokes identity and parity hold to machine precision. The expansion is checked through its O(1/k_e) convergence. Non-reciprocity of the exact second variation is checked to 3% (asserted; 1–1.5% measured). The full first-order theory predicts the orientation dependence of the learned kernels in the real closed loop to 3% in size and 11% as a vector (cosine 0.994, k_e = 40). | no precedent found (targeted search, round 3); ingredients known (higher-order averaging; performative prediction [S14]; semi-gradient TD [S10], [S11]) |
| H.8 | (𝔄) is input-to-state stable under drift perturbations, with gain 1/(τρ_τ) | proved, assuming the chain rule for the perturbed flow | elementary given [S3] |

**Where the new mathematics lies.** The new mathematics this viewpoint requires lies in the open problems of §8.

1. Infinite-horizon averaging for McKean–Vlasov closed loops. Propositions H.7 and H.8 reduce it to controlling a
   near-identity averaging transformation along the moving measure.
2. Uniform-in-time propagation of chaos with a plant in the loop.
3. The theory of the non-reciprocal, semi-gradient mean-field dynamics that the closed loop actually produces.
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

**Theorem E.** There are r₀, k, λ > 0 such that the following holds.

- k and λ come from the converse Lyapunov function of the linearisation (proof steps 3–4).
- r₀ depends in addition on r := min(ρ, δ), and on the Lipschitz constant of the Jacobian on the ball D of step 2.

Put c₀ := r₀²/(2(1 + mλ_G)) and C := √2·k·max{1,(mλ_G)^{1/2}}·max{1,(mλ_min(G))^{−1/2}}. If V_OT(t₀) ≤ c₀, then:

- the optimal assignment π at t₀ is unique, and it is the unique optimal assignment for all large t;
- ‖e(t)‖ + W(μ_t,μ*) ≤ Ce^{−λ(t−t₀)}(‖e(t₀)‖ + W(μ_{t₀},μ*)) for t ≥ t₀;
- ϑ_i(t) → ϑ*_{π(i)}. The labelled error E_lab(t) := ((1/m)Σ‖ϑ_i(t) − ϑ*_i‖²_{G⁻¹})^{1/2} tends to
  ((1/m)Σ‖ϑ*_{π(i)} − ϑ*_i‖²_{G⁻¹})^{1/2}, which is non-zero for π ≠ id.

k and λ depend on Ψ(·) itself, through the converse Lyapunov function of its linearisation, and on m, G and k_e. r₀
and c₀ depend in addition on a_max, s_min, ρ and δ. The cited theorems are of existence type, so it is not proved
here that k and λ can be chosen uniformly over the class of Ψ with given (T, α, sup‖Ψ‖, sup‖Ψ̇‖).

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

3. *LTV stability.* Put w := (e, ζ) with ζ := m^{−1/2}(I_m⊗G^{−1/2})ϑ̃, and Ψ_G := m^{−1/2}Ψ(I⊗G^{1/2}).

   - The system becomes ė = −k_ee − Ψ_Gζ and ζ̇ = Ψ_Gᵀe, that is, ẇ = [A₀ + B(t)]w with A₀ = diag(−k_eI, 0) and B
     skew-symmetric, bounded, with bounded derivative.
   - The excitation transfers: ∫Ψ_GᵀΨ_G ⪰ (αλ_min(G)/m)I.
   - Morgan–Narendra [S8, Thm 1] characterises uniform asymptotic stability of [[A, −B(t)], [B(t)ᵀ, 0]] with
     A + Aᵀ ≺ 0 by richness of B. Apply it with A = −k_eI and B = Ψ_G (bounded, with bounded derivative): the origin of
     the w-system is uniformly asymptotically stable.
   - (PE) supplies the richness, since ∫_t^{t+T}‖Ψ_G w‖ ≥ αλ_min(G)/(m·sup‖Ψ_G‖) for every unit w.
   - For linear systems UAS is UES (Khalil Thm 4.11). UES transfers to (e, ϑ̃) through the constant similarity
     diag(I, m^{1/2}(I⊗G^{1/2})).
   - The bound on Ψ̇ cannot be dropped. Ψ = sin t² is bounded and PE, but for ė = −e − Ψζ, ζ̇ = Ψe, from w(t₀) = (e, ζ) = (0, 1), the ratio
     ‖w(t₀+200)‖/‖w(t₀)‖ is 0.9715, 0.9942 and 0.9980 for t₀ = 5, 20 and 50 (review round 2). So that system is not
     UES.
   - We could not check the exact richness definition and regularity class of [S8] against its full text. Its
     abstract's characterisation matches this system.

4. *Nonlinear system.* Apply Khalil Thm 4.13 in z = (e, ϑ̃) on [0,∞) × D. The Jacobian is bounded and Lipschitz on D
   uniformly in t, and the linearisation is UES. This gives r₀ ≤ r, k and λ, with ‖z(t)‖ ≤ k‖z(t₀)‖e^{−λ(t−t₀)}
   whenever ‖z(t₀)‖ < r₀. These trajectories stay in D, where the closed loop equals F.

5. *Transfer to (‖e‖, W).* If V_π(t₀) ≤ c₀, then

       ‖z(t₀)‖² ≤ 2max(1, mλ_G)V_π(t₀) ≤ 2max(1, mλ_G)c₀ < 2(1 + mλ_G)c₀ = r₀².

   The bound W ≤ (mλ_min(G))^{−1/2}‖ϑ̃‖ and the optimality of π at t₀ give C.

   *Uniqueness.* Suppose V_π, V_π' ≤ c₀ with π ≠ π'. Then the state lies in Khalil balls around both (0, ϑ*∘π) and
   (0, ϑ*∘π'). The solution would converge to two distinct equilibria, which is impossible.

   *Limits.* Hence z(t) → 0, that is, ϑ_i(t) → ϑ*_{π(i)}. For π' ≠ π,
   V_π'(t) → (1/2m)Σ‖ϑ*_{π(i)} − ϑ*_{π'(i)}‖²_{G⁻¹} ≥ d²/2 > 0, while V_π(t) → 0. So π is the unique optimal assignment
   for all large t, and E_lab converges as stated. ∎

*Remark E.1 (the rate).*

- **What converges.** "The parameters converge" is false for the labelled error and true for the matched one.
- **The rate is not set by the PE level.** λ > 0 exists for each admissible Ψ, but it depends on Ψ(·), not only on
  its PE level. No positive function of αλ_min(G)/m, k_e, T and sup‖Ψ_G‖ alone bounds it below.
  - For m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = I, k_e = 3 and x_m = 1.5 sin ωt (T = 2π), these four quantities do not
    depend on ω. The exact Floquet rate falls from 4.7e-2 (ω = 1) to 2.9e-5 (ω = 64), like ω⁻² (review round 3).
  - A lower bound must also involve sup‖Ψ̇_G‖ (step 3), and none is proved here.
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
  second variation are non-zero in the checked instances, including the network of the experiments (review round 3).

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
- They do not vanish in the checked instances, for n = 1 (from A_even alone) and for n = 2.
- So in these instances, at finite k_e, the exactly averaged loop is a non-reciprocal McKean–Vlasov dynamics. It is not
  the W_{G⁻¹}-gradient flow of any free energy. A different, μ-dependent metric is not excluded by this argument.

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
   k_e → ∞ at fixed ε gives cyc = 0.

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

**Proposition H.8 (input-to-state stability of (𝔄)).** Assume σ, τ > 0 and the setting of Theorem H. Let μ_t solve
the perturbed equation

    ∂_tμ = ∇·(μ(G∇V[μ] − b_t)) + τ∇·(G∇μ)

from μ_0 ∈ P₂ with F_τ(μ_0) < ∞. Here b_t is a velocity perturbation with β(t)² := ∫‖b_t‖²_{G⁻¹}dμ_t ≤ β̄². Assume the
chain rule of [S3] Lemma 2.2 holds along this solution: d/dt F_τ(μ_t) = ∫⟨∇(δF_τ/δμ)(μ_t), v_t⟩dμ_t, in whitened
coordinates, with v_t the velocity. It is expected for smooth bounded perturbations (parabolic regularity and the
entropy-dissipation estimate), but it is not verified here. Then

    d/dt (F_τ(μ_t) − F*) ≤ −τρ_τ (F_τ(μ_t) − F*) + ½β(t)²,
    F_τ(μ_t) − F* ≤ e^{−τρ_τ t}(F_τ(μ_0) − F*) + β̄²/(2τρ_τ),
    limsup_t W_{G⁻¹}(μ_t, μ*_τ) ≤ β̄/(τρ_τ).

*Proof.* Work in ζ = G^{−1/2}ϑ, and let ξ_t := ∇(δF_τ/δμ)(μ_t) = ∇Ṽ[μ_t] + τ∇log μ_t. The velocity is −ξ_t + b̃_t,
with ‖b̃_t‖ = ‖b_t‖_{G⁻¹}.

1. By the chain rule and Young's inequality,
   d/dt F_τ(μ_t) = −∫‖ξ_t‖²dμ_t + ∫⟨ξ_t, b̃_t⟩dμ_t ≤ −½∫‖ξ_t‖²dμ_t + ½β².
2. With μ̂_t ∝ e^{−Ṽ[μ_t]/τ},
   ∫‖ξ_t‖²dμ_t = τ²I(μ_t|μ̂_t) ≥ 2ρ_ττ²H(μ_t|μ̂_t) ≥ 2ρ_ττ(F_τ(μ_t) − F*),
   by the LSI of Theorem H(ii) and the entropy sandwich ([S3] Lemma 3.4).
3. Grönwall gives the first two lines.
4. For the third, F_τ(μ) − F* ≥ τH(μ|μ*_τ) ([S3] Lemma 3.4). μ*_τ satisfies LSI(ρ_τ), hence Talagrand's inequality
   W² ≤ (2/ρ_τ)H(·|μ*_τ) (Otto–Villani). So W² ≤ (2/(τρ_τ))(F_τ − F*), and the limsup follows from the second line.
   ∎

The gain 1/(τρ_τ) ≤ e^{α*/τ}/σ is exponentially large in 1/τ, like the inverse of the rate.

*Route to Conjecture H.6.* In the scaling of H.6:

- τρ_τ ≥ (σ₀/k_e)·e^{−2K_Φ(F_ρ+K_Φ)/τ₀}, so 1/(τρ_τ) = O(k_e).
  - In this scaling the whitened proximal Gibbs measures do not depend on k_e, so τρ_τ = Θ(1/k_e).
  - At the experiment's τ₀ = 0.04 the constant is e^{1049}.
- Define b := (1/k_e²)[G∇(U_even + U_odd) + b_osc] + O(k_e⁻³). The first part is the frozen-measure correction of
  Proposition H.7; the second is the oscillation term of Remark H.7a. Both are bounded uniformly in μ and ϑ, so
  β̄ = O(k_e⁻²).
- Proposition H.8 then confines every solution of (𝔄) + b to an asymptotic W-ball around μ*_τ of radius
  O(1/k_e)·e^{2K_Φ(F_ρ+K_Φ)/τ₀}/σ₀. This needs F_τ(μ_0) < ∞, or a P₂ start, after the Girsanov regularisation of
  H(ii), which applies because b is bounded.

Three things are still missing for H.6:

1. control of the frozen-measure approximation along the moving μ_t, that is, a near-identity averaging transformation
   for the stochastic McKean–Vlasov loop, uniformly in time;
2. the mean-field limit of (L_τ), since Theorem D covers only (L);
3. the chain rule that H.8 assumes for the perturbed flow.

**Evidence** (the real closed loop with m = 32 and n = 1, not the averaged system):

- **Dither is necessary.** Without dither the collapsed start is trapped forever. Its loss stays at 0.166, against
  0.007 from the spread start (`results/langevin_tradeoff`, 1500 s).
- **Equalisation at moderate temperature.** At σ = 0.02 and τ = 0.004 (3 seeds,
  `results/langevin_global/table_sigma0.02_tau0.004.md`):
  - collapsed and spread starts reach the same loss (ratio 1.03) and tracking RMS (0.0465 against 0.0452);
  - their final 32-particle configurations, driven by the same noise realisation, are 4.2× closer in W₂ over the
    physical coordinates (a, c, s) than without dither (0.40 against 1.69);
  - two independent 32-point samples of μ*_τ are typically about 0.6 apart, so 0.40 is at the sampling floor
    (review round 3);
  - the common loss of about 0.11 matches the grid value L̄(μ*_τ) = 0.11118. It is dominated by the regularisation
    bias.
- **Metastability at low temperature.** At σ = 0.002 and τ = 5·10⁻⁴ (`table_sigma0.002_tau0.0005.md`), dither frees
  the collapsed start (0.166 → 0.048). But the starts have not equalised within 900 s (loss ratio 5.2). The spread
  start's 9.2·10⁻³ is close to the regularised optimum L̄(μ*_τ) = 8.5·10⁻³.
- **Limits of the evidence.**
  - The Langevin runs have n = 1. By Proposition H.7, the averaged dynamics is already non-reciprocal for n = 1,
    through the state shift.
  - The only n ≥ 2 closed-loop runs are the deterministic orientation runs of Remark H.7a. They confirm the first-order
    theory quantitatively. No n ≥ 2 run with dither has been made.
  - The runs use fixed k_e, not the scaling of H.6.
  - All Langevin runs use k_e = 10, below the Tikhonov scale of Remark H.0 for this network
    ((F_ρ + K_Φ)/s_min ≈ 35). The measured frozen-μ drift error at k_e = 10 is 4.2% (n = 1), so the condition is
    conservative, but the agreement with (𝔄) is empirical.
  - So the evidence supports Conjecture H.6 qualitatively, but the equalisation time can exceed practical horizons at
    low τ.

---

## 7. Verification map

Test files are abbreviated: `T` is `tests/test_transport.py`, `Lg` is `tests/test_langevin.py`, and `Av` is
`tests/test_averaging.py`.

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
| D | consistency in m (nearly constant amplification 1.64–1.69; not the constant e^{ΛT}) | `exp_mean_field_m.py` |
| F (ii) along (L) | particle path length (≥ Wasserstein length) ≤ λ_G^{1/2}L₁∫‖e‖ at every sample | `T::test_theorem_F_wasserstein_length_bounded_by_integrated_error` |
| F (ii′) along (L_σ) | W(μ_{t₁},μ_{t₀}) ≤ ∫(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰)) on a grid of pairs; a tight instance with e ≡ 0 (ratio ≥ 0.99), where the motion is W ≈ 1.05 | `T::test_theorem_F_leak_clause_along_ot_sigma_mod`, `T::test_theorem_F_leak_clause_is_tight_when_the_leak_pays_for_all_motion` |
| E linearisation | Jacobian of the implemented loop equals A(t); QA + AᵀQ = diag(−k_eI, 0) | `T::test_theorem_E_linearisation_and_lyapunov_identity` |
| E convergence | the matched error decays and the labelled error does not (qualitative); sliding-window excitation and fitted rates | `T::test_theorem_E_matched_error_converges_labelled_does_not`, `exp_pe_convergence.py` |
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
| H.7 + H.7a in the real closed loop | the orientation dependence of the learned kernels (same ρ) is predicted to within 3% in size and 11% as a vector, cosine 0.994 at k_e = 40; the measured difference scales like k_e^{−1.01} | `exp_orbit_orientation.py` → `results/orbit_orientation/table_sigma2_slow16.md` |

Not checked in the repository:

- B.5 along (L), and B.6 and B.7;
- the Filippov form of C(i), the inequality of C(v), and Proposition C.5's D⁺V_OT = σ_m(δ − δ²);
- D(b), D's constant, and the one-sided constant of Remark D.1;
- the rate of E, Corollary E.2 and Remark E.3;
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
   - *First step:* the deterministic case (τ = 0, σ > 0) at fixed m, where the needed averaging theorem is classical
     (Sanders–Verhulst–Murdock Ch. 5, "averaging under attraction").

2. **Non-reciprocal mean-field dynamics from closed-loop learning (new; Proposition H.7).**
   - *Proved at a frozen measure:* the averaged adaptive dynamics is a semi-gradient flow plus a lag term
     (Proposition H.7). Along the real loop a second-order averaging term is added (Remark H.7a).
   - Its interaction has two antisymmetric parts:
     - the performative defect A_even of the state shift, which is present for every n and is the only part for n = 1;
     - the circulation kernel K, an integral over the region that the reference orbit encloses (n ≥ 2).
   - Non-reciprocity is checked in the instances tested, not proved in general.
   - *Open questions:*
     - Are stationary states unique beyond the perturbative regime?
     - Can circulation be *designed*, choosing the orbit's orientation and area, to break symmetric traps or to speed
       escape from the collapsed configuration?
     - Is there a Lyapunov functional other than F_τ for the non-reciprocal flow?
   - *Context:* the closest analogues are performative prediction [S14], semi-gradient TD learning [S10], [S11] and
     non-reciprocal phase transitions [S12].
   - *Done (Remark H.7a):* deterministic n = 2 closed-loop runs with the reference traversed both ways (same ρ). The
     full first-order theory predicts the O(1/k_e) difference of the learned kernels to 3% in size and 11% as a
     vector.
   - *Next steps:* the same test with dither, where the stationary state is a measure; and the design of the orbit's
     orientation and area as a control input.

3. **Uniform-in-time propagation of chaos with a plant in the loop.**
   - Theorem D is finite-horizon. Uniform-in-time results exist for mean-field Langevin dynamics [S5] and for shallow
     networks [36].
   - *Open:* the closed loop, where the plant couples all particles through e, and where (by H.7) the averaged
     interaction is not reciprocal.

4. **Geometry of two-time-scale adaptation (Proposition G(b′)).**
   - *Proved formally, and checked:* the slow geometry follows the gradient flow of the profile loss in the metric
     induced on the slow manifold, and the slope S measures the departure from plain variable projection.
   - *Open:* a rigorous long-horizon statement; the O(1/R) corrections, which involve D²a* and ∂_ξM; gain designs that
     shape the induced metric.

5. **Noiseless global convergence.**
   - *Open:* a closed-loop separation condition under which the noiseless flow avoids spurious stationary points,
     despite the failure of the Chizat–Bach boundary conditions. Candidates are a compactified width coordinate and
     reflecting boundaries (Remark H.3).

6. **Unbalanced transport, meaning online structure adaptation.**
   - Treat the weights as masses (Hellinger–Kantorovich / Wasserstein–Fisher–Rao) to obtain birth and death of kernels
     with Lyapunov guarantees. The conic particle descent of [49] and [55] is the static analogue.

7. **Excitation theory in measure space.**
   - *Open:* the continuum form of (PE); the trajectory kernel K_T(ϑ,ϑ') = ∫Φ(x(s);ϑ)Φ(x(s);ϑ')ds; identifiability of
     measures up to symmetry; the relation to [16]. Remark E.1 shows that the PE level alone does not predict the rate.
     A quantitative, uniform UES estimate for the skew-symmetric class is a first target.

8. **Transport-regularised adaptive laws.**
   - Wasserstein proximal (JKO) discrete-time adaptive laws.
   - The OT leak toward priors with a different number of atoms, using optimal plans with barycentric projection.
   - Transfer: transporting a learned measure to a new plant.

9. **Quantisation and the curse of dimensionality (Corollary B.6).**
   - *Open:* second-order (centroid) quantisation rates for the ideal measure on the lifted space; approximation rates
     that do not go through W (Maurey-type); adaptive laws that track optimal quantisers (the centroidal Voronoi
     tessellation connection [14]).

10. **Symmetric systems beyond kernels.** The OT quotient applies verbatim to swarms of identical agents, mixtures and
    multi-head models ([88] and [89] are the control-side precedents).

11. **Deep approximators.** Hierarchical (multi-layer) mean-field measures, with a layer-wise OT Lyapunov function.

---

## 9. References

See `LITERATURE_TRANSPORT.md`: [1]–[89] are in `review/LITERATURE_refs.md`; [S1]–[S14] are in the audit file.

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
