# Literature audit for the transport paradigm (PARADIGM.md)

**Method.** Six independent primary-source sweeps, each a different search angle: optimal transport in
control, mean-field networks, neuro-adaptive control, leakage and permutation symmetry, propagation of
chaos, and Wasserstein geometry. A synthesis agent deduplicated their sources into
`review/LITERATURE_refs.md` (89 entries with URLs; 79 unsourced entries dropped). The per-sweep verdicts
are in `review/LITERATURE_sweeps.md`. The final merge agent hit a session limit, so this document is the
merge, written from the sweep verdicts. It adds sources I verified directly for Theorems G and H, marked
[S1]–[S9], and, in round 3, sources for Propositions B.6, H.7 and H.8, marked [S10]–[S14]. Numbers [1]–[89]
refer to `review/LITERATURE_refs.md`.

The audit was adversarial toward novelty. A verdict of "novel" means only that none of the searches
found the construction.

---

## Supplementary sources verified for Theorems G and H

- [S1] P. Marion, R. Berthier. *Leveraging the two-timescale regime to demonstrate convergence of neural networks*. NeurIPS 2023; arXiv:2304.09576.
- [S2] R. Barboni, G. Peyré, F.-X. Vialard. *Ultra-fast feature learning for the training of two-layer neural networks in the two-timescale regime*. arXiv:2504.18208 (2025).
- [S3] L. Chizat. *Mean-field Langevin dynamics: exponential convergence and annealing*. TMLR 2022; arXiv:2202.01009. Assumptions 1–3, Theorem 3.2, Corollary 3.3, Theorem 4.1 and Proposition 5.1 were read in the PDF.
- [S4] A. Nitanda, D. Wu, T. Suzuki. *Convex analysis of the mean field Langevin dynamics*. AISTATS 2022, PMLR 151:9741–9757; arXiv:2201.10469.
- [S5] F. Chen, Z. Ren, S. Wang. *Uniform-in-time propagation of chaos for mean field Langevin dynamics*. arXiv:2212.03050 (v4, 2025).
- [S6] L. Chizat, F. Bach. NeurIPS 2018 (= [29]). Assumptions 3.4(iii), the boundary conditions for the partially 1-homogeneous case, were read in the PDF. They require Neumann conditions on a bounded convex domain, or C¹-uniform behaviour at infinity. Gaussian kernels with free centres and widths satisfy neither (see the PARADIGM.md §6c preamble). Its Theorem 3.5 also needs the support of μ_0 to separate {−r₀}×Θ from {r₀}×Θ; a Dirac start violates this.
- [S7] G. H. Golub, V. Pereyra. *Separable nonlinear least squares: the variable projection method and its applications*. Inverse Problems 19 (2003).
- [S8] A. P. Morgan, K. S. Narendra. *On the stability of nonautonomous differential equations ẋ = [A + B(t)]x with skew symmetric matrix B(t)*. SIAM J. Control Optim. 15(1):163–176 (1977). The round-1 Theorem E judge identified this as the correct citation for the mixed error model. B. D. O. Anderson, IEEE TAC 22 (1977), covers the identifier equations.
- [S9] Standard references: Ambrosio–Gigli–Savaré, *Gradient Flows* (2008); Khalil, *Nonlinear Systems* (3rd ed.); Kokotović–Khalil–O'Reilly, *Singular Perturbation Methods in Control* (1986); Sastry–Bodson, *Adaptive Control* (1989); Graf–Luschgy, *Foundations of Quantization for Probability Distributions* (2000); Holley–Stroock (1987); Bakry–Émery (1985); Otto–Villani, J. Funct. Anal. 173 (2000) (LSI implies Talagrand); Ambrosio–Savaré (2007), Handbook of Differential Equations: Evolutionary Equations 3; Sanders–Verhulst–Murdock, *Averaging Methods in Nonlinear Dynamical Systems* (2nd ed., 2007).
- [S10] Y. Zhang, Q. Cai, Z. Yang, Y. Chen, Z. Wang. *Can temporal-difference and Q-learning learn representation? A mean-field theory*. NeurIPS 2020; arXiv:2006.04761. A mean-field analysis of TD learning with two-layer networks, whose update is not the gradient of any objective.
- [S11] J. N. Tsitsiklis, B. Van Roy. *An analysis of temporal-difference learning with function approximation*. IEEE TAC 42(5):674–690 (1997). TD is a semi-gradient method: it holds the bootstrapped target fixed. The data distribution is that of a fixed policy and does not depend on the parameters.
- [S12] M. Fruchart, R. Hanai, P. B. Littlewood, V. Vitelli. *Non-reciprocal phase transitions*. Nature 592:363–369 (2021). Non-reciprocal interactions in many-body systems.
- [S14] J. C. Perdomo, T. Zrnic, C. Mendler-Dünner, M. Hardt. *Performative prediction*. ICML 2020, PMLR 119; arXiv:2002.06673. Here predictions shift the data distribution, and retraining with the distribution held at its current value converges to performatively stable points rather than optimal ones. This is the closest structural analogue of PARADIGM.md H.7(c).
- [S13] J. Chevallier. *Uniform decomposition of probability measures: quantization, clustering and rate of convergence*. J. Appl. Probab. 55(4):1037–1045 (2018). Equal-weight approximation rates in W₂ (used in PARADIGM.md Corollary B.6), with Fournier–Guillin, Probab. Theory Relat. Fields 162 (2015), for i.i.d. atoms.

---

## Claim-by-claim verdicts

### A. The adaptive law is a Wasserstein gradient flow of a linear potential

- **Closest prior.**
  - Particle gradient flows coincide with Wasserstein gradient flows from atomic data: [29] Prop. B.1 and [28]. The potential-energy flow transports particles along −∇V: [71] §4.5.
  - RBF units as Wasserstein particles in confined training: [41]. Projection as a support constraint: [46] Prop. 3.1, in a control venue.
  - Time-dependent linear potential with confinement and exogenous data: [38], [48].
  - Adaptive laws read as gradient or Onsager flows of the law of one parameter vector: [25], [17], [18]. An ensemble of estimate particles in Wasserstein space: [26], whose full text could not be read.
- **Distinction.** No source states that the classical centre/width laws with 1/m gains are the characteristics of the continuity equation for a μ-linear potential whose coefficient e(t) comes from plant feedback.
- **Strongest objection.** For a μ-linear potential the "gradient flow" is the continuity equation with velocity −∇V. Its characteristics are the particle ODEs by definition. The identification is therefore essentially tautological, and its ingredients are published.
- **Verdict: partially anticipated.** The contribution is interpretive. It is a change of viewpoint that licenses Theorems B–H, not a theorem with new mathematics. PARADIGM.md says so.

### B. Optimal-assignment Lyapunov function V_OT = ½‖e‖² + ½W²

- **Closest prior.**
  - The minimum over permutations of matched parameter error equals W₂ between neuron measures: [62], [63], [64], [66].
  - Labelled Taylor remainders with explicit curvature constants in neuro-adaptive control: [2], [19], [20], [21], [23].
  - Label-free local analyses with width-independent constants in machine learning: [49], [50], [54].
  - OT costs inside Lyapunov arguments for swarms, where the measure is the plant: [88], [89].
- **Distinction.** No closed-loop Lyapunov function couples the tracking error to W²(μ_t, μ*). None handles assignment switching through Dini derivatives. None has the two-regime remainder with a drift term V_*W.
- **Strongest objection.** The labelled inequality of [20] holds for every relabelling, so minimising over permutations gives the matched bound almost for free. The derivative of W² along transport flows is standard (AGS Thm 8.4.7). m-uniformity comes from the normalisation, not from OT; the round-1 review confirmed this.
- **Verdict: partially anticipated.** The packaging is new; every ingredient is known. The genuine content is uniformity over labellings and across assignment switches, which fixed labellings cannot provide.

### C. OT σ-modification

- **Closest prior.**
  - Labelled leakage: [3], [4], [5].
  - Offline permutation-aware regularisation toward a reference network: [61], [64], [65], [67], [68].
  - Proximal OT to the previous iterate: [45].
  - Curvature obstruction and assignment switching: [71] §4.4 two-atom example, [72], [73] (AGS Example 9.1.5, Lemma 9.2.1, Thm 7.3.2).
  - Local convexification near a target measure: [76].
- **Distinction.** No source found has online-recomputed assignment leakage inside a Lyapunov-based adaptive law, a strict σ-modification inequality with bias W(μ*, μ⁰), an assignment-gap condition, or a reshuffling result for the labelled leak (PARADIGM.md C(v)).
- **Strongest objection.** C(i) and the obstruction are textbook. The gap condition is LP sensitivity of the optimal assignment. Three probes in this direction were blocked.
- **Verdict: novel as a control construction, medium confidence.** C(i) and Remark C.2 are cited as known.

### D. Mean-field closed loop and Dobrushin estimate

- **Closest prior.**
  - ODE coupled to a continuity equation, with a Dobrushin-type estimate and mean-field limit: [77] Prop. 3.2 and Thm 3.3, [78] Thm 4.1, [79], [80].
  - Dobrushin's estimate and exact empirical solutions: [82], [83].
  - Particle systems that act as a controller: [85]. Closed-loop mean-field learners: [43], [44].
  - Uniform-in-time propagation of chaos: [36], [37], [S5].
  - Limit-interchange caveat: [81], [86].
- **Distinction.** None of these couples a tracking plant to the continuity equation of the controller's own parameter measure.
- **Strongest objection.** D is an instance of [77]–[79], with W₂ in place of W₁. Its constant e^{ΛT} is bounded-horizon only.
- **Verdict: anticipated in technique.** It is stated as an application.

### E. Label-free local parameter convergence under tangential excitation

- **Closest prior.**
  - Static analogues: [49] (non-degeneracy, local exponential convergence, label-free), [50], [52], [54] (lifted Gaussian weight–centre–width space), [57].
  - Control-side Jacobian-excitation convergence for a labelled vector: [12], [20], [22].
  - Identifiability only up to permutation: [58], [59], [60].
- **Distinction.** Here PE is a condition on the time-varying tangent regressor along the closed-loop reference. The result is stated for the matched error jointly with e, and the labelled error provably does not converge.
- **Strongest objection.** This is [49] Assumption 5 transplanted into a closed loop, with PE replacing non-degeneracy.
- **Verdict: partially anticipated.**

### F. Transport cost of adaptation (Wasserstein length ≤ λ_G^{1/2} L₁ ∫‖e‖)

- **Prior.** The length bound for absolutely continuous curves is AGS Thm 8.3.1. No sweep targeted F.
- **Verdict: elementary.** It is presented as a corollary with a control interpretation, not as a contribution.

### G. Variable-projection reduction of the slow geometry

- **Closest prior.** Two-timescale (variable-projection) training of two-layer networks: [S1] (global convergence in a univariate setting) and [S2] (mean-field variable projection with rates). Variable projection itself: [S7]. Averaging and singular perturbation in adaptive control: [S9].
- **Distinction (corrected in round 2).**
  - What the closed loop supplies: the error-channel and averaging scales, and the data distribution (the occupation measure of the reference). The loss is the occupation-measure loss, and the controller learns on its own excitation.
  - What the designer supplies: the separation between weights and geometry is a gain choice, δ ≪ g_a·min(m₀/‖Γ^{1/2}∇²L*Γ^{1/2}‖, 1/‖Da*Γ^{1/2}‖²) with (g_c, g_s) = δ(γ_c, γ_s), exactly as in [S1], [S2]. Version 2's claim that the loop produces this separation was wrong.
  - What is new: the leading-order slow flow, valid for every slope S, is the gradient flow of the profile loss in the metric induced on the slow manifold (PARADIGM.md G(b′)). It reduces to plain variable projection only when the slope S = (δ/g_a)‖Da*Γ^{1/2}‖² is small. `exp_slow_manifold_slope.py` gives a reproducible counterexample to plain VarPro under the rate condition alone.
  - L* is the *unregularised* profile loss: the λ → 0 limit of λℒ^λ of [S2] (time rescaled). It is not [S2]'s λ → 0 functional ℒ⁰, which is min-norm interpolation with an ultra-fast-diffusion flow. Version 2's statement that "the geometry marginal follows the Wasserstein gradient flow of the variable-projection functional of [S2]" was false: as a functional of the geometry marginal, L* depends only on the support, so it is not W-continuous and has no mean-field form.
- **Strongest objection.** Each reduction step is a standard theorem on its own time scale, and their composition on the O(k_e/(mδ)) horizon is formal. The induced-metric correction is a standard geometric singular-perturbation computation.
- **Verdict: a fixed-m, closed-loop analogue of the two-timescale/VarPro training of [S1], [S2].** The induced-metric slow flow and the slope condition (H2c), with its counterexample, are the content added here. My searches were standard, not an exhaustive sweep, so confidence is medium-low.

### H. Langevin adaptive control: the averaged mean-field loop is mean-field Langevin dynamics

- **Closest prior.**
  - [25] (Patil, Oct 2026): projected adaptation as an Onsager flow, σ-modification as confinement, and reflected Langevin as a constrained free-energy flow, for linearly parametrised models and the law of one parameter vector. The functional is linear in the law, and no global result is claimed.
  - [24] (LyLA-Therm 2025): Langevin parameter updates in a labelled deep-network controller, with ultimate boundedness in probability.
  - [38]: online mean-field learning with a time-varying linear potential and uniform-in-time propagation of chaos; the data are exogenous and there is no plant.
  - Mean-field Langevin theory: [S3], [S4], [S5].
- **Distinction.** The measure here is the empirical measure of kernels with nonlinear mean-field coupling, and the functional is convex but not linear in μ. The loss is the occupation-measure loss of the closed loop. Three control devices map onto the hypotheses of [S3]: saturated lifted features give bounded features, σ-modification in the G⁻¹ metric gives Bakry–Émery confinement, and parameter dither gives entropy. The log-Sobolev constant comes out in control quantities (σ, τ, k_e, a_max, F_ρ).
- **Strongest objection.** The global convergence is [S3]/[S4] applied verbatim. The step from the real closed loop to the averaged system is formal (PARADIGM.md Remark H.0), and the statement for the real loop is Conjecture H.6, restated in round 2 in the scaling σ = σ₀/k_e, τ = τ₀/k_e.
- **Verdict: novel as a bridge.** None of the searches found it in adaptive control.

### H.7 Structure of the exactly averaged closed loop: semi-gradient and orbit-area circulation

- **Closest prior.**
  - Second-order averaging of periodically forced systems [S9] (Sanders–Verhulst–Murdock).
  - Performative prediction [S14]: gradient steps that hold the model-dependent data distribution fixed. This is the structural analogue of H.7(c).
  - Semi-gradient TD [S11], which holds the bootstrapped target fixed, and its mean-field analysis with two-layer networks [S10].
  - Non-reciprocal interactions and the failure of free-energy descriptions in many-body systems [S12].
- **Distinction.** No source found states three things together.
  - The first correction to the averaged adaptive closed loop is a semi-gradient of the closed-loop loss on the closed-loop state distribution, and its performative defect is non-reciprocal.
  - It includes an antisymmetric kernel–kernel interaction for n ≥ 2. Its kernel is an integral of Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ')] over the region that the reference orbit encloses (Stokes).
  - So at finite feedback gain the averaged adaptive dynamics is, in general, not the W_{G⁻¹}-gradient flow of any free energy. This is checked for n = 1 and n = 2; the state-shift (semi-gradient) part needs no enclosed area.
- **Search record (round 3).** Three targeted searches found no combined precedent:
  - "non-reciprocal mean-field interaction adaptive control averaging closed-loop neural network parameters not gradient flow";
  - "adaptive control as semi-gradient learning closed-loop state distribution depends on parameters two-timescale averaging non-gradient correction";
  - "averaging adaptive systems periodic reference second-order correction antisymmetric term area enclosed orbit geometric phase parameter drift".

  A fourth search confirmed [S10].
- **Strongest objection.** The expansion is routine second-order averaging at a frozen measure. The term from the intra-period parameter oscillation (PARADIGM.md Remark H.7a) is the classical second-order averaging term. The semi-gradient reading is an interpretation, and the circulation kernel is Stokes' theorem.
- **What is specific here.**
  - The closed-loop consequence: the learned kernels depend on the direction in which the reference orbit is traversed, at order 1/k_e, although the occupation measure, and with it the averaged mean-field system, is unchanged.
  - The first-order theory predicts this effect quantitatively in the real closed loop: to 3% in size, with cosine 0.994, at k_e = 40.
- **Verdict: no precedent found; the ingredients are known.** Confidence is low to medium, because these were three targeted searches, not a sweep.

### E.0 Explicit strict Lyapunov function for the skew-symmetric mixed error model

- **Closest prior.** Strict Lyapunov functions for persistently excited systems built from a backward–forward integral
  of the PE Gramian: Mazenc–Malisoff, *Constructions of Strict Lyapunov Functions* (Springer 2009), and Mazenc, de Queiroz
  and Malisoff, IEEE TAC 54 (2009), for adaptively controlled systems. Explicit UES bounds for the identifier
  ζ̇ = −ΨᵀΨζ: Anderson (1977), Loría–Panteley (2002). The qualitative statement for the mixed class is Morgan–Narendra
  [S8], cited by version 3.
- **Distinction.** The lemma gives, for the mixed class ė = −ke − Ψζ, ζ̇ = Ψᵀe, a closed-form rate and overshoot in
  (k, α, T, Ψ_M, Ψ_D), with the cancellation of the ‖Ψζ‖² terms between the Gramian cross term and the eᵀΨζ cross term
  as its only idea, and it quantifies its own conservativeness (Remark E.0a: a factor 3κ(1 + Ψ_D/(kΨ_M))² in the
  slow-adaptation regime, κ the condition number of the averaged Gramian).
- **Verdict: a known construction made explicit for this class.** Its role is rigour, not novelty: Theorem E no longer
  depends on an unverified citation. The ω⁻² law of Remark E.1 now has a proof.

### I. Reciprocity criterion for the exactly averaged loop

- **Closest prior.** The coboundary (cyclic-sum) criterion for reciprocal pairwise interactions is standard in the
  non-reciprocal many-body literature [S12]; the injectivity of the Gaussian convolution (Weierstrass transform) and
  Stokes' theorem with winding numbers are classical analysis.
- **Distinction.** No source found states that for a feature family containing the zero feature the criterion collapses
  to the vanishing of the antisymmetric kernel, nor the two characterisations: for n = 1 the averaged adaptive loop is
  reciprocal at order k_e⁻² iff the residual vanishes on the orbit; for n ≥ 2 at perfect fit iff the orbit's 1-current
  vanishes. The decomposition of the performative defect into the curl of the residual and a Lie bracket of the feature
  fields was not found either.
- **Verdict: new, with classical tools; the content is specific to the closed-loop learner.** The result settles
  README open problem 3 of version 3 (a general non-reciprocity theorem beyond checked instances) for n = 1 and for the
  perfect-fit case, and reduces the general case to an explicit bilinear form.

### H.8 Input-to-state stability of the averaged system

- **Prior.** The proof is the energy-dissipation argument of [S3] Thm 3.2, with the perturbation absorbed by Young's inequality and the limit transferred by Talagrand (Otto–Villani).
- **Verdict: elementary given [S3].** Its role is to reduce Conjecture H.6 to an averaging estimate.

---

## What this means for the claims in PARADIGM.md

The paradigm's contribution is a viewpoint with a set of correct, machine-checked theorems that connect
adaptive control of nonlinearly parametrised approximators to optimal transport and mean-field learning
theory. Within it, the constructions that appear new are:

- the closed-loop OT Lyapunov function (B, its packaging);
- OT σ-modification with its reshuffling result (C);
- the induced-metric form of the occupation-measure variable-projection reduction (G);
- the mean-field Langevin bridge with its global convergence statement (H);
- the frozen-measure structure of the exactly averaged loop (H.7), and the characterisation of its non-reciprocity (I).

The tools are known. The open problems in PARADIGM.md §8 are where new mathematics is required.

## Open questions this audit could not settle

1. [26] (ICML 2026) could not be read in full. Its relation to A and B is unknown.
2. Assignment-free formation control is a likely source of min-over-permutation Lyapunov functions with Dini derivatives. Candidates are Zavlanos–Pappas and Turpin–Michael–Kumar (CAPT). The probes were blocked. This could reduce B's novelty further.
3. The averaging and singular-perturbation literature on neuro-adaptive control was not swept exhaustively for G.
4. Whether [S5]'s hypotheses hold verbatim for the averaged closed loop of H (the finite-m statement, now the conditional Remark H.5).
5. Whether the averaging literature on periodically forced gradient systems already contains the orbit-area (Stokes) term of Proposition H.7 in some form, for example as geometric-phase or Hannay-type corrections. This literature was not swept.
