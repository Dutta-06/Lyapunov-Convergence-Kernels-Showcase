# Literature audit for the transport paradigm (PARADIGM.md)

**Method.** Six independent primary-source sweeps, each a different search angle: optimal transport in
control, mean-field networks, neuro-adaptive control, leakage and permutation symmetry, propagation of
chaos, and Wasserstein geometry. A synthesis agent deduplicated their sources into
`review/LITERATURE_refs.md` (89 entries with URLs; 79 unsourced entries dropped). The per-sweep verdicts
are in `review/LITERATURE_sweeps.md`. The final merge agent hit a session limit, so this document is the
merge, written from the sweep verdicts. It adds sources I verified directly for Theorems G and H, marked
[S1]–[S9], and, in round 3, sources for Propositions B.6, H.7 and H.8, marked [S10]–[S14]; round 5 adds [S15]–[S17] for
Lemma E.0, and round 6 adds [S23]–[S34] for Theorem J′ (bibliographic details from memory, not re-verified offline). Numbers [1]–[89]
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
- [S15] F. Mazenc, M. Malisoff. *Constructions of Strict Lyapunov Functions*. Springer, Communications and Control Engineering, 2009. Strict Lyapunov functions for persistently excited systems built from a backward–forward integral of the PE Gramian. (Bibliographic details from memory; not re-verified against the publisher's record in this session.)
- [S16] F. Mazenc, M. de Queiroz, M. Malisoff. *Uniform global asymptotic stability of a class of adaptively controlled nonlinear systems*. IEEE Trans. Automat. Control 54(5), 2009. (Details from memory; not re-verified here.)
- [S17] A. Loría, E. Panteley. *Uniform exponential stability of linear time-varying systems: revisited*. Systems & Control Letters 47(1):13–24, 2002; and B. D. O. Anderson, *Exponential stability of linear equations arising in adaptive identification*, IEEE Trans. Automat. Control 22(1):83–88, 1977. Explicit UES bounds for persistently excited identifiers. (Details from memory; not re-verified here.)
- [S18] J. Cortés, S. Martínez, T. Karatas, F. Bullo. *Coverage control for mobile sensing networks*. IEEE Trans. Robotics and Automation 20(2):243–255, 2004. Continuous-time Lloyd descent ṗ_i = −k(p_i − C_{V_i}) toward the centroids of Voronoi cells, with its LaSalle convergence proof. (Details from memory; not re-verified here.)
- [S19] M. Balzer, T. Schlömer, O. Deussen. *Capacity-constrained point distributions*. ACM SIGGRAPH 2009; Q. Mérigot, *A multiscale approach to optimal transport*, Computer Graphics Forum 30(5), 2011; F. de Goes, K. Breeden, V. Ostromoukhov, M. Desbrun, *Blue noise through optimal transport*, ACM Trans. Graphics 31(6), 2012. Equal-mass (capacity-constrained) quantisation through optimal transport plans and their barycentres. (Details from memory; not re-verified here.)
- [S20] S. P. Lloyd. *Least squares quantization in PCM*. IEEE Trans. Inform. Theory 28(2):129–137, 1982; Q. Du, V. Faber, M. Gunzburger, *Centroidal Voronoi tessellations: applications and algorithms*, SIAM Review 41(4):637–676, 1999. Lloyd's algorithm and centroidal Voronoi tessellations. (Details from memory; not re-verified here.)
- [S21] A. Saberi, H. Khalil. *Quadratic-type Lyapunov functions for singularly perturbed systems*. IEEE Trans. Automat. Control 29(6):542–550, 1984; and H. K. Khalil, *Nonlinear Systems*, 3rd ed., Prentice Hall, 2002, Ch. 11 (composite Lyapunov functions for two-time-scale systems). (Details from memory; not re-verified here.)
- [S22] D. Aeyels, J. Peuteman. *A new asymptotic stability criterion for nonlinear time-variant differential equations*. IEEE Trans. Automat. Control 43(7):968–971, 1998; J. Peuteman, D. Aeyels, *Exponential stability of slowly time-varying nonlinear systems*, Math. Control Signals Systems 15:202–228, 2002. Stability of slowly varying and averaged systems through the averaged vector field. (Details from memory; not re-verified here.)
- [S13] J. Chevallier. *Uniform decomposition of probability measures: quantization, clustering and rate of convergence*. J. Appl. Probab. 55(4):1037–1045 (2018). Equal-weight approximation rates in W₂ (used in PARADIGM.md Corollary B.6), with Fournier–Guillin, Probab. Theory Relat. Fields 162 (2015), for i.i.d. atoms.

---
- [S23] E. Pardoux, A. Yu. Veretennikov. *On the Poisson equation and diffusion approximation. I*. Ann. Probab. 29(3):1061–1085, 2001 (and *II*, Ann. Probab. 31(3):1166–1192, 2003). Existence, uniqueness and polynomial growth bounds for solutions of the Poisson equation Lu = f − π(f) of a uniformly elliptic diffusion with a recurrence (confining) condition; the tool behind the invariant-measure expansion of Theorem J′. (Details from memory; not re-verified here.)
- [S24] D. Gilbarg, N. S. Trudinger. *Elliptic Partial Differential Equations of Second Order*. Springer, 2nd ed., 1983 (reprinted 2001), Theorems 6.2 and 6.17: interior Schauder estimates and higher regularity. (Details from memory; not re-verified here.)
- [S25] L. Hörmander. *Hypoelliptic second order differential equations*. Acta Math. 119:147–171, 1967; and M. Hairer, *On Malliavin's proof of Hörmander's theorem*, Bull. Sci. Math. 135(6–7):650–666, 2011 (a self-contained probabilistic proof for diffusions with smooth coefficients of bounded derivatives). (Details from memory; not re-verified here.)
- [S26] P. Cattiaux, L. Mesnager. *Hypoelliptic non-homogeneous diffusions*. Probab. Theory Relat. Fields 123:453–483, 2002. Smooth transition densities for time-inhomogeneous diffusions under the parabolic Hörmander condition. (Details from memory; not re-verified here.)
- [S27] D. W. Stroock, S. R. S. Varadhan. *On the support of diffusion processes with applications to the strong maximum principle*. Proc. Sixth Berkeley Symp. Math. Statist. Prob. III:333–359, 1972. The support of the law of a diffusion is the closure of the skeleton (controlled-ODE) trajectories. (Details from memory; not re-verified here.)
- [S28] G. Da Prato, J. Zabczyk. *Ergodicity for Infinite Dimensional Systems*. Cambridge University Press, LMS Lecture Notes 229, 1996, Theorem 4.2.1 (Doob–Khasminskii: strong Feller + irreducible ⇒ at most one invariant measure, and convergence of the laws) and §4.3. (Details from memory; not re-verified here.)
- [S29] U. Seifert. *Stochastic thermodynamics, fluctuation theorems and molecular machines*. Rep. Prog. Phys. 75:126001, 2012; and C. Maes, K. Netočný, *Time-reversal and entropy*, J. Stat. Phys. 110:269–310, 2003. Entropy production rate of a stationary diffusion, ∫jᵀD⁻¹j/ρ, and its vanishing exactly in detailed balance. (Details from memory; not re-verified here.)
- [S30] G. A. Pavliotis, A. M. Stuart. *Multiscale Methods: Averaging and Homogenization*. Springer, Texts in Applied Mathematics 53, 2008, Chapters 10–11 (averaging for SDEs with fast deterministic or stochastic variables via the perturbation expansion of the backward Kolmogorov equation; the solvability conditions used in Theorem J′'s Step 3). (Details from memory; not re-verified here.)
- [S31] J. C. Mattingly, A. M. Stuart, M. V. Tretyakov. *Convergence of numerical time-averaging and stationary measures via Poisson equations*. SIAM J. Numer. Anal. 48(2):552–577, 2010. Perturbation of invariant measures through the Poisson equation of the limiting generator; the structure of Theorem J′'s Step 5 (there for numerical discretisations, here for the averaging of a hypoelliptic closed loop). (Details from memory; not re-verified here.)
- [S32] M. Hairer. *Ergodic properties of Markov processes*. Lecture notes, University of Warwick, 2006 (Lemma 1.13, Corollary 1.14, Theorem 1.16: strong Feller + a point in the support of every transition probability ⇒ uniqueness of the invariant measure; minorisation on compacts). (Details from memory; not re-verified here.)
- [S33] M. Hairer, J. C. Mattingly. *Yet another look at Harris' ergodic theorem for Markov chains*. In: Seminar on Stochastic Analysis, Random Fields and Applications VI, Progr. Probab. 63:109–117, Birkhäuser, 2011. Geometric ergodicity in a weighted total variation norm from a Lyapunov function and a minorisation on its sublevel sets. (Details from memory; not re-verified here.)
- [S34] A. Dembo, O. Zeitouni. *Large Deviations Techniques and Applications*. Springer, 2nd ed., 1998, Theorem 6.2.10 (Sanov's theorem) and Theorem 4.3.1 (Varadhan's lemma). The large-deviation principle for the empirical measure under a mean-field Gibbs measure with a bounded continuous energy, used in Corollary J′.3. (Details from memory; not re-verified here.)
- [S35] R. L. Dobrushin. *Vlasov equations*. Functional Analysis and Its Applications 13:115–123, 1979. Existence, uniqueness and the representation by characteristics of the Vlasov equation with a bounded Lipschitz vector field depending Lipschitz-continuously (bounded-Lipschitz / W₁ metric) on the measure; stability in the measure. Used for the mean-field equation (𝔄₀) of Theorem D‴. (Details from memory; not re-verified here.)
- [S36] C. McDiarmid. *On the method of bounded differences*. In: Surveys in Combinatorics 1989, London Math. Soc. Lecture Note Ser. 141:148–188, Cambridge Univ. Press, 1989. The bounded-differences (McDiarmid) inequality: a function of m independent variables that changes by at most c_j when the j-th is changed satisfies P(f − Ef ≥ t) ≤ exp(−2t²/Σc_j²). Used for the reference's sampling fluctuation in Theorem D‴, Step 8. (Details from memory; not re-verified here.)
- [S37] F. Malrieu. *Logarithmic Sobolev inequalities for some nonlinear PDE's*. Stochastic Processes and their Applications 95(1):109–132, 2001. Uniform-in-time propagation of chaos for McKean–Vlasov diffusions with a uniformly convex confinement and a convex interaction, by synchronous coupling and functional inequalities. The convexity there is global; Theorem D‴ needs it on a neighbourhood of the support of the mean-field flow and gets the interaction's sign from flat convexity. (Details from memory; not re-verified here.)
- [S38] A. Durmus, A. Eberle, A. Guillin, R. Zimmer. *An elementary approach to uniform in time propagation of chaos*. Proceedings of the AMS 148:5387–5398, 2020. Uniform-in-time propagation of chaos by reflection coupling for McKean–Vlasov diffusions whose confinement is convex outside a ball and whose interaction is Lipschitz with a small constant. The interaction of the fitted network is not small in the learning regime; its sign (Theorem D‴, Step 3) replaces the smallness. (Details from memory; not re-verified here.)
- [S39] C. Liu, L. Zhu, M. Belkin. *Loss landscapes and optimization in over-parameterized non-linear systems and neural networks*. Applied and Computational Harmonic Analysis 59:85–116, 2022. The decomposition of the Hessian of a least-squares loss into the Gram (tangent-kernel) term and the residual-weighted second-derivative term, the latter vanishing at interpolation; the PL condition near interpolation. The same decomposition, per particle, is the self term of Theorem D‴'s coupling identity. (Details from memory; not re-verified here.)

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
  of the PE Gramian: [S15], and [S16] for adaptively controlled systems. Explicit UES bounds for the identifier
  ζ̇ = −ΨᵀΨζ: [S17]. The qualitative statement for the mixed class is Morgan–Narendra [S8], cited by version 3.
- **Distinction.** The lemma gives, for the mixed class ė = −ke − Ψζ, ζ̇ = Ψᵀe, a closed-form rate and overshoot in
  (k, α, T, Ψ_M, Ψ_D), with the cancellation of the ‖Ψζ‖² terms between the Gramian cross term and the eᵀΨζ cross term
  as its only idea, and it quantifies its own conservativeness (Remark E.0a: a factor 3κ(1 + Ψ_D/(kΨ_M))² in the
  slow-adaptation regime, κ the condition number of the averaged Gramian).
- **Verdict: a known construction made explicit for this class.** Its role is rigour, not novelty: Theorem E no longer
  depends on an unverified citation. The ω⁻² law of Remark E.1 now has a proof.

### E.0′ Sharp explicit rate in the slow-adaptation regime

- **Closest prior.** Composite (quadratic-type) Lyapunov functions for singularly perturbed systems [S21], where a fast
  variable is slaved to a slow manifold and the interconnection terms are absorbed with sup-norm constants; averaging
  Lyapunov functions for slowly time-varying systems [S22], where a zero-mean antiderivative of the fluctuation of the
  vector field is added as a correction term. Lemma E.0′ combines the two for the mixed class (LTV): a second-order
  slaving of the error (e₂ = −εΨζ + ε²Ψ̇ζ − ε³ΨMζ) and two averaging corrections (P₁, P₂) of the slow quadratic form.
- **Distinction.** The interconnection and remainder terms are bounded not in sup norm but per direction against the
  average excitation ‖v‖_{M̄} — the *peak factors* ρ, ρ_D, ρ_DD — which is what makes the rate λ_min⟨ΨᵀΨ⟩/k (the exact
  averaged rate) rather than quadratic in λ_min; the identity N₂ = ½⟨P₁M + MP₁⟩ = 0 shows that the second-order averaged
  correction to the symmetric part vanishes for this class, which explains the measured k⁻² relative deviation of the
  exact rate from λ_min/k. No source stating an explicit rate for the mixed PE class that is sharp in the slow limit, or
  using peak factors of the regressor's derivatives as the threshold data, was found (the averaging and PE literatures
  were not swept exhaustively).
- **Verdict: a combination of two classical constructions with a new choice of norms; the sharp explicit rate for the
  mixed class appears new.** Its scope is the slow regime above an explicit threshold; the fast regime stays with
  Lemma E.0.

### E.0″ The slow regime for non-periodic PE regressors

- **Claim.** Lemma E.0′'s construction with the period mean replaced by a forward exponential average over a window T,
  whose derivative identity TṀ̄ = M̄ − M gives the antiderivatives P₁ = −2T(M̄ − M̂) and S = −T(F̄ − F̂) in closed form; the
  periodic cancellation N₂ = 0 is replaced by a residual N̄₂ whose price δ₂ = εn₂/λ is O(1/k); the rate is (1 − δ − δ₂)λ_T/k with
  λ_T the exponentially weighted PE level of the window.
- **Prior.** Strict Lyapunov functions for PE systems from windowed or exponentially weighted integrals of the regressor are
  the Mazenc–Malisoff construction [S15], [S16]; averaged (auxiliary) systems with sliding or weighted windows appear in the
  stability theory of slowly varying systems [S22]; the exponentially weighted PE condition is equivalent to classical PE
  (elementary). Verdict: the ingredients are standard; what is new is the same as for Lemma E.0′ — the rate is the averaged
  rate λ_T/k up to 1 − δ − δ₂, with all constants computable from the regressor — now without periodicity, which is what the
  Lorenz-tracking experiments of Part 1 need.

### I. Reciprocity criterion for the exactly averaged loop

- **Closest prior.** The coboundary (cyclic-sum) criterion for reciprocal pairwise interactions is standard in the
  non-reciprocal many-body literature [S12]; the injectivity of the Gaussian convolution (Weierstrass transform) and
  Stokes' theorem with winding numbers are classical analysis.
- **Distinction.** No source found states that for a feature family containing the zero feature the criterion collapses
  to the vanishing of the antisymmetric kernel, nor the characterisations: for n = 1 the averaged adaptive loop is
  reciprocal at order k_e⁻² iff the residual vanishes on the orbit; for n ≥ 2 at perfect fit iff the orbit's 1-current
  vanishes and the residual's Jacobian is symmetric along the orbit; and no open set of measures is reciprocal. The
  decomposition of the performative defect into the curl of the residual and a Lie bracket of the feature fields was
  not found either.
- **Verdict: new, with classical tools; the content is specific to the closed-loop learner.** The result settles
  README open problem 3 of version 3 (a general non-reciprocity theorem beyond checked instances) for n = 1 and for the
  perfect-fit case, and reduces the general case to an explicit bilinear form.

### C′. The barycentric OT leak toward a prior with a different number of atoms

- **Closest prior.** The pure law ϑ̇_i = −σ_m(ϑ_i − T̄_i) is the continuous-time Lloyd descent of coverage control
  [S18] with the Voronoi cell replaced by the cell of an equal-mass optimal plan, that is, capacity-constrained
  quantisation through optimal transport [S19]; OT-based coverage control is [89]; Lloyd's algorithm and centroidal
  Voronoi tessellations are [S20]. The repository's earlier pointer [14] is an application of CVT to adaptive
  estimation, not the quantisation theory itself.
- **Distinction.** The law is used as a σ-modification leak inside an adaptive controller, and the analysis is the
  closed-loop one: the dissipation inequality up to the within-cell variance, the induced-coupling bound √C_bar ≤ b + 2η,
  the Lyapunov inequality with the bias b_bar (the ideal's centroid defect), and the transfer reading (anchoring m kernels
  to a learned measure of m₀ ≠ m kernels). The ω-limit statement for the pure leak is proved by LP stability for
  generalised solutions (round 5, C1), not by the LaSalle argument of [S18], whose Voronoi cells vary continuously.
- **Verdict: the construction is known in coverage control and quantisation; new as a leak in adaptive control, with
  the induced-coupling and Lyapunov analysis.** The §0 and README rows say so.

### D′. Propagation of chaos for the dithered closed loop

- **Prior.** Sznitman's synchronous coupling for McKean–Vlasov systems with Lipschitz drifts, whose standard output is
  the m^{−1/2} rate for the particles and for averages of bounded functions — a dimension-free exponent — with the
  dimension entering only through the empirical-measure term (Fournier–Guillin, [S13]).
- **Distinction.** The plant is in the loop: the mean field enters through a deterministic state x̄ driven by the
  empirical mean of the features, and the constants are control quantities (K_Φ, the Lipschitz bounds of the saturated
  features, k_e enters only through the tube). The per-particle sup bound is the pathwise Grönwall estimate (round 5, D1).
- **Verdict: the technique and the dimension-free exponent are standard; the closed-loop statement with its constants
  is new in this setting.** Uniformity in time is open (PARADIGM.md §8, branch 3).
- **Theorem D″ (uniform in time in the leak-dominated regime).** The contraction argument with a confining potential is the
  standard route to uniform-in-time propagation of chaos for McKean–Vlasov systems (e.g. the small-interaction regime of
  Malrieu 2001/2003, and [S5] for mean-field Langevin dynamics); the closed-loop version adds the plant as a second
  contracting channel and a small-gain condition. Its certified regime is the one in which the leak stops the learning
  (round 6, DPP-2, DPP-3): the result is a statement about where the method stops. Verdict: standard technique, new
    object; the learning regime, where the long-horizon numerics show the coupled distance growing, is open.

### D‴. Uniform-in-time propagation of chaos in the learning regime (the averaged slow loop)

- **Claim.** For the G-preconditioned gradient flow of J_m = mL̄(μ_Θ) + leak (the averaged slow loop at τ = 0), synchronously
  coupled to the characteristic flow of its mean-field equation, the mean-square coupling distance D obeys ½Ḋ ≤ −(σ₀ − κ_R)D +
  (L₂^G)²D²/16 + L_ϑ^G√D·Y (+ L₂^GYD), the interaction term being −‖f̂_Θ − f̂_Θ̄‖² at linear order; under the strict second-order
  optimality condition along the mean-field flow (a neighbourhood of the support, margin κ, from time s₁), sup_s E D ≤
  4e^{2Λ₁}(L_ϑ^G)²K_Φ²/(κ²m) + an exponentially small tail (PARADIGM.md §5d).
- **Prior.** Uniform-in-time propagation of chaos by synchronous coupling under *global* convexity of the confinement and
  convexity of the interaction is Malrieu's [S37]; non-convex confinements with *small* Lipschitz interactions are handled by
  reflection coupling [S38]; flat-convex energies with dither by the entropy method [S5]. The Hessian decomposition
  "Gram + residual-weighted feature curvature" is the standard tool of the interpolation / PL-condition literature [S39] and
  of lazy-training analyses, where the second term vanishes at interpolation.
- **Distinction.** (i) The interaction of the fitted network is neither small nor convex in the particle coordinates; its
  sign comes from the flat convexity of L̄ (f̂ linear in the measure), which makes −‖f̂_Θ − f̂_Θ̄‖² the leading interaction
  term under the coupling — the finite-m, two-measure version of [S39]'s decomposition. (ii) The convexity needed is local:
  σ₀I − G^{1/2}R[μ_s]G^{1/2} ⪰ κI on a neighbourhood of the support of the mean-field flow, i.e. the strict second-order
  optimality condition of the fitted configuration, evaluated at the actual residual — not a global convexity (which would
  force the collapsed configuration to be stable and the network not to learn, Theorem D″). (iii) The window κ_R[μ_∞] <
  σ₀ < σ_crit between the collapsed instability (Proposition K) and the fitted stability is a design statement about
  the leak that has no counterpart in the cited work. (iv) The limits are stated with the result: a connected initial
  support that splits between basins violates the hypothesis for all times (the characteristic flow is a homeomorphism),
  and the coupling then has a floor of pairs in different basins; the real loop at finite k_e is not covered.
- **Verdict: the technique (synchronous coupling + bootstrap + bounded differences) is standard; the sign structure of the
  interaction and the local, residual-weighted hypothesis — and with them a uniform-in-time statement in a regime where the
  network learns — are new in propagation of chaos for adaptive control.** Numerics: `exp_chaos_learning.py`.

### J. Periodic orbits of the finite-m closed loop near the critical points of the averaged loss

- **Prior.** Two-time-scale averaging under attraction (Khalil Ch. 10–11, Sanders–Verhulst–Murdock Ch. 5,
  Sastry–Bodson Ch. 4): hyperbolic equilibria of the averaged system correspond to periodic orbits of the original
  system with the same stability type. The periodic Riccati block-triangularisation of the monodromy matrix is the
  standard slow-manifold reduction for singularly perturbed linear periodic systems.
- **Distinction.** The theorem identifies the averaged system with the gradient flow of the m-particle regularised
  occupation loss J_m (the finite-m form of Theorem H's functional at τ = 0, in the scaling σ = σ₀/k_e), gives the
  O(1/k_e) location of the orbit including the second-order oscillation term of Remark H.7a, and the k_e-uniform basin.
  No source applying this to an adaptive controller read as a measure was found.
- **Verdict: classical tools, applied to a new object; the first rigorous real-loop ↔ averaged-system statement of the
  program.** The remainder regularity was supplied in round 5 (J-1).

### K. Collapsed configurations and the orientation of the orbit

- **Prior.** The isotypic (symmetry-adapted) splitting of a Hessian at a permutation-symmetric point is standard
  equivariant linear algebra; the stability condition σ₀G⁻¹ ≻ R₂ is the local form of the Chizat–Bach first-variation
  optimality condition at a Dirac mass [S3]-type analyses; higher-order averaging corrections to Floquet exponents are
  classical (Sanders–Verhulst–Murdock).
- **Distinction.** The reading of the orientation of the reference orbit as a design input — an orientation-odd O(k_e⁻²)
  shift of the separation exponents that decides the stability of the collapsed configuration in a window around the
  critical leak, computable from the averaged data — was not found. (The orientation effect itself is Remark H.7a's; the
  round-3 search found no precedent.)
- **Verdict: new as a design statement; the mathematics is classical.** Checked in one instance.

### J′. The stationary law of the dithered finite-m closed loop

- **Claim.** With dither τ = τ₀/k_e and leak σ = σ₀/k_e, the finite-m closed loop has an invariant law whose
  parameter marginal is the Gibbs measure exp(−J_m/τ₀) up to O(1/k_e), with an explicit first-order correction driven by
  the first-order averaged drift of Theorem J (Proposition H.7 + Remark H.7a) as a vector field on parameter space; the
  stationary state carries a probability current of order 1/k_e that vanishes exactly when that drift is reciprocal, and the
  stationary parameter process is not time-reversible.
- **Ingredients, all known:** averaging of SDEs with fast variables through the backward-Kolmogorov expansion [S30];
  perturbation of invariant measures through the Poisson equation of the limiting generator [S31], [S23]; hypoellipticity
  and uniqueness [S25]–[S28]; the thermodynamics of non-equilibrium steady states (currents, entropy production,
  detailed balance) [S29]; second-order averaging [Sanders–Verhulst–Murdock].
- **Verdict:** no source found that applies the invariant-measure expansion to an adaptive control loop (a hypoelliptic
  system with a deterministic fast error channel, a periodic reference and a slow Langevin parameter process), nor one that
  identifies the first-order non-equilibrium current of a learning loop with a reciprocity defect of the averaged dynamics.
  The closest objects are the stationary measures of performative/strategic learning dynamics [S14] and the non-reciprocal
  many-body steady states of [S12]; neither has a plant in the loop or a 1/k_e expansion. The proof technique is standard
  (Pavliotis–Stuart; Mattingly–Stuart–Tretyakov); the result is new in its object and in the identification of the
  correction with Theorem I's defect. As with Theorem J, the non-commutativity term ½⟨[∇²𝒥, ∇²𝒰]_G⟩ in the antisymmetric part
  of the second-order averaged drift of a periodic gradient flow is presumably known in the averaging literature; no source
  was located.
- **Corollary J′.3** (the Gibbs measures π_m are chaotic with limit the minimiser of Theorem H's free energy) is Sanov's
  theorem with Varadhan's lemma [S34] applied to a bounded continuous energy; the statement is classical in form
  (Messer–Spohn 1982 for mean-field Gibbs measures); new only as the bridge between Theorems J′ and H.
- **Not checked numerically:** the uniqueness of the invariant measure (Step 1); the circulation of the real loop (below the
  Monte Carlo resolution at the affordable cost; checked on the first-order averaged SDE).

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
- the frozen-measure structure of the exactly averaged loop (H.7), and the characterisation of its non-reciprocity (I);
- the finite-m averaging theorem for the real closed loop (J) and the orientation as a design input for the collapsed
  configuration (K);
- the explicit strict Lyapunov functions for the mixed PE class (E.0, and E.0′ with the sharp slow-regime rate through
  peak factors);
- the closed-loop propagation-of-chaos statement (D′) and the barycentric leak's Lyapunov analysis (C′), both built on
  standard constructions (Sznitman's coupling; the Lloyd descent of coverage control).

The tools are known. The open problems in PARADIGM.md §8 are where new mathematics is required.

## Open questions this audit could not settle

1. [26] (ICML 2026) could not be read in full. Its relation to A and B is unknown.
2. Assignment-free formation control is a likely source of min-over-permutation Lyapunov functions with Dini derivatives. Candidates are Zavlanos–Pappas and Turpin–Michael–Kumar (CAPT). The probes were blocked. This could reduce B's novelty further.
3. The averaging and singular-perturbation literature on neuro-adaptive control was not swept exhaustively for G.
4. Whether [S5]'s hypotheses hold verbatim for the averaged closed loop of H (the finite-m statement, now the conditional Remark H.5).
5. Whether the averaging literature on periodically forced gradient systems already contains the orbit-area (Stokes) term of Proposition H.7 in some form, for example as geometric-phase or Hannay-type corrections. This literature was not swept.
