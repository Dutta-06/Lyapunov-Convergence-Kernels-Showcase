# Paper roadmap

The results of `PARADIGM.md`, grouped into papers that can be written from this repository. For each paper: the claim, the results it rests
on with their status as `PARADIGM.md` labels them (proved / formal / sketched / numerical), the novelty as `LITERATURE_TRANSPORT.md` judges
it, what is still missing before submission, and a natural venue. The statuses and the verdicts are abbreviated from those files, not upgraded:
where a result is formal, sketched or numerical, the paper inherits that, and the exact wording is the one in `PARADIGM.md`. Every number a paper would quote is produced by a script in `experiments/` and recorded in
`results/RESULTS_PARADIGM.md`; every result has been through the adversarial review logged in `review/REVIEW_LOG.md`.

## Core papers

### 1. Adaptive kernel laws as Wasserstein gradient flows, and the transport Lyapunov function

- **Claim.** The gradient adaptive law for Gaussian kernels is the particle form of a Wasserstein gradient flow; ½‖e‖² + ½W² is a
  Lyapunov function uniform over labellings and assignment switches; the OT σ-modification and its Lloyd-type form toward priors with any
  number of atoms.
- **Rests on.** Theorems A, B, C (proved and checked), C′ (proved for vertex-selection solutions, their existence assumed, and for
  generalised solutions), Theorem F (proved and checked).
- **Novelty (audit).** A: partially anticipated (interpretive). B: packaging new, ingredients known; the content is the uniformity over
  labellings. C: novel as a control construction (medium confidence). C′: known in coverage control and quantisation, new as a leak in
  adaptive control.
- **Missing.** A comparison with the classical projection and e-modification laws on a common benchmark beyond the repository's
  instances.
- **Venue.** SIAM J. Control Optim.; Math. Control Signals Systems.

### 2. Label-free convergence under measure-space persistent excitation

- **Claim.** Local exponential convergence without labels under tangential persistent excitation, with explicit radius, overshoot and rate;
  the sharp slow-adaptation rate; non-periodic regressors by exponential averaging.
- **Rests on.** Theorem E with Lemma E.0 (proved), Lemma E.0′ and Corollaries E.0′.2–3 (proved), Lemma E.0″ (proved); Remark E.1 (the PE
  level alone does not predict the rate).
- **Novelty (audit).** E: a known construction made explicit for this class; E.0′/E.0″: classical constructions with a new choice of norms
  and a sharp explicit rate.
- **Missing.** The continuum form of (PE) and the trajectory kernel (§8, branch 7, open).
- **Venue.** Automatica; IEEE Trans. Automatic Control.

### 3. Propagation of chaos for adaptive control loops

- **Claim.** The mean-field closed loop is well posed; finite-horizon propagation of chaos for the dithered loop; uniform in time under a
  dominant leak; uniform in time in the learning regime via the flat convexity of the loss.
- **Rests on.** Theorem D (proved: local well-posedness, global existence under a self-consistent growth condition), D′ (proved), D″
  (proved), D‴ (proved for the averaged slow loop; its hypothesis (H) a computable property, checked numerically on the instance),
  Proposition D‴.5 (formal: the algebra exact, the law by the delta method), Corollary D‴.7.
- **Novelty (audit).** D′: standard technique, new closed-loop statement with constants; D‴: standard technique (synchronous coupling,
  bootstrap), new sign structure of the coupling identity.
- **Missing.** The constants of D‴ are a worst-case envelope (Remark D‴.6); an instance where (H) is proved rather than checked.
- **Venue.** Annals of Applied Probability; SIAM J. Math. Anal.; Stochastic Processes Appl.

### 4. Langevin adaptive control: global convergence, non-reciprocity and the stationary law

- **Claim.** With saturated features, σ-modification and dither the averaged loop is mean-field Langevin dynamics and converges globally to
  the regularised optimum; the exactly averaged loop is generically non-reciprocal; the finite-m dithered loop has a stationary law whose
  parameter marginal is the Gibbs measure up to O(1/k_e), with a non-equilibrium current.
- **Rests on.** Theorem H (proved for the averaged loop; the reduction to it formal, the real loop Conjecture H.6), H.7, H.8, Theorem I
  (proved), Theorem J′ and Proposition J′.2 (proved, with three standard theorems imported; the invariant law unique under a smoothness
  hypothesis), Corollary J′.3.
- **Novelty (audit).** H: novel as a bridge (none found in adaptive control); I: no precedent found; J′: no source applies the
  invariant-measure expansion to an adaptive loop.
- **Missing.** Conjecture H.6 (the real loop's global statement) reduced to an averaging estimate but open.
- **Venue.** SIAM J. Control Optim.; Annals of Applied Probability.

### 5. The finite-m closed loop: periodic orbits, collapse and the critical leak

- **Claim.** Every hyperbolic critical point of the regularised finite-m loss carries a periodic orbit of the real loop at distance
  O(1/k_e); collapsed configurations and the critical leak σ_crit; the lifted K₀-atom orbit and its m-independent Floquet spectrum; Theorem
  J uniformly in m in the per-particle norm.
- **Rests on.** Theorem J (proved), Proposition K (proved; (iii)–(iv) given the second-order normal form), Corollary D‴.14 and
  Proposition D‴.15 (proved; the exponents to O(k_e⁻²) by first-order averaging), Theorem D‴.18 (proved as an adaptation of Theorem J's
  proof, at its level of detail).
- **Novelty (audit).** J: classical tools on a new object (the first rigorous real-loop ↔ averaged-system statement for this learner);
  K: new as a design statement.
- **Missing.** Nothing structural; the constants are not explicit.
- **Venue.** SIAM J. Applied Dynamical Systems; Nonlinearity.

### 6. Measure-valued leaks: the finite-size floor, the mass constraint and its basin

- **Claim.** A leak toward a measure-valued prior is a design input of the convergence window; the nearest-atom leak conserves the cluster
  masses and leaves a finite-size floor whose law is the weight-derivative of the leak's bias; the barycentric OT leak's mass constraint
  removes it at commensurate sizes, with an m-independent Polyak–Łojasiewicz basin; the entropic constraint equalises the masses only down
  to a critical excess. On the instance the floor at the stable samples is the exact bias at the sampled weights to four digits, while the
  delta-method law is asymptotic (2.8, 1.5, 1.1 times below the exact floor at m = 15, 60, 240) and a third of the small samples are
  saddles of the fixed prior, where a cluster splits.
- **Rests on.** Propositions D‴.9 ((i) proved), D‴.12 ((i)–(ii) proved, the coupling statement formal), D‴.17 (proved; no rate), D‴.20
  ((i)–(iii) proved, (iv) under a positivity hypothesis that fails at large excesses and small ε), Theorem D‴.21 (proved, (iii) for
  vertex-selection solutions; the a-priori basin small on the instance; the mean-field form and the real loop formal), Proposition
  D‴.24(i) (F exact at q, the bias formula first order, the law (B) formal), Remark D‴.25 (the floor's three levels: function, coupling,
  measure).
- **Novelty (audit).** D‴.12: textbook ingredients, new mechanism (a conservation law as the source of the floor); D‴.21: classical
  inequality and geometry, new mechanism for a new energy; D‴.24: classical ingredients (the delta method, the filter-factor bound, the
  vanishing deviation at an interpolating limit), new mechanism (the floor as the bias's slope, mode by mode).
- **Missing.** A coupling theorem from an i.i.d. start through the transit (§8, branch 8, open); incommensurate sizes beyond the K₀ − 1
  split particles.
- **Venue.** SIAM J. Control Optim.; J. Machine Learning Research (the mean-field-network side).

### 7. Adaptive priors: removing the leak's bias without losing the leak

- **Claim.** A prior whose atoms are learned removes a misplaced prior's bias mode by mode at the rate of the measure-space excitation,
  uniformly in m; with the masses minimised out (the nearest-atom leak toward an adaptive prior) it removes the function-level floor at
  every representable sample, and the fixed prior's splitting instability with it — a property of the amplitude–mass gauge, so it
  survives an approximation error; perfect tracking is then an exponentially stable periodic solution of the real loop with an exact,
  m-independent separation exponent and collective exponents whose expansion in 1/k_e is even, with an explicit second-order term; the
  transit's drag and its gate; the limits (the K₀-cluster bias, the lost a-priori bound) and the remedy for the second (a projected
  prior, convergent at η = ∞).
- **Rests on.** Theorem D‴.23 ((i)–(iii) proved, (iv) formal, (v) sketched, (vi) partly numerical), Proposition D‴.24 ((i) exact at q
  with the law formal, (ii)–(iii), (vi) proved; the convergence numerical; checked in one and two dimensions), Corollary D‴.26 ((a)–(b)
  exact, (c) a formal expansion with its structure proved; the basin open), Proposition D‴.27 ((i), (iii) proved, (ii)'s convergence proved
  at η = ∞ and formal at finite η), Proposition D‴.28 ((i), (iii) exact, (ii) as D‴.24(i); the dynamics numerical).
- **Novelty (audit).** D‴.23: standard per-mode formulas (proximal point, iterated Tikhonov, elastic averaging), new transfer to the
  mass-constrained particle system uniformly in m, new failure mode with its remedy; D‴.24: see paper 6, and a learned prior as the way to
  restore the interpolation property under regularisation; D‴.26: classical tools (averaging, the Magnus commutator, e-modification), a
  new structural statement (the exact m-independent spectrum at every representable sample, the even expansion); D‴.27: classical
  ingredients (projection, the switching σ-modification, which it is at η = ∞), transferred to the learned centre of the mass-constrained
  prior at finite η; D‴.28: a classical symmetry (the conic view's homogeneity), a new identification of the mechanism and an exact margin
  formula.
- **Missing.** The constants of D‴.23(v); the nonlinear basin of D‴.26 and a bound on its O(k_e⁻⁴) remainder; the convergence of D‴.27 at
  finite η; critical points at infinity for D‴.28 (best approximants on the saturation boundary); the choice of the projection's set as a
  design input; the
  instance's samples without an exact representation are not rare at small m (0.36 at m = 15), so a paper should treat them (the
  expected limiting floor is then of the fixed prior's order).
- **Venue.** IEEE Trans. Automatic Control; Automatica.

## Application papers

### 8. Neural-network tracking control with transport leaks

- **Claim.** On plants with unknown smooth dynamics, Gaussian-kernel adaptive controllers with an OT or adaptive-prior leak track with no
  finite-size floor and no leak bias, in one and two dimensions, with the averaged predictions holding in the real loop.
- **Rests on.** The real-loop numerics of `exp_ot_real_loop.py`, `exp_ot_2d.py`, `exp_lifted_orbit.py`, `exp_adaptive_floor.py` (Part D),
  `exp_adaptive_floor_2d.py` (Part C), `exp_adaptive_floquet.py`; the theorems of papers 5–7, with Proposition D‴.27 for the bounded
  parameters a controller needs.
- **Missing.** Plants beyond sums of Gaussian bumps (e.g. a pendulum, a quadrotor's attitude), measurement noise, input saturation, and
  comparisons with σ-/e-modification and projection on identical geometry (the ASSESSMENT.md protocol).
- **Venue.** IEEE Trans. Control Systems Technology; Int. J. Adaptive Control Signal Processing.

### 9. Learned mixture priors for compressing kernel networks

- **Claim.** The fast limit of the adaptive prior is a capacity-constrained K₀-means regulariser (soft weight-sharing with prescribed
  masses); it pulls the network to K₀ clusters without biasing where they sit, but biases how many there are.
- **Rests on.** Theorem D‴.23(iv) (formal) and Part G of `exp_adaptive_prior.py` (numerical).
- **Missing.** The number of atoms as a learned quantity (split/merge); benchmarks against soft weight-sharing.
- **Venue.** A machine-learning venue (e.g. TMLR).

### 10. A design guide for leaks in mean-field adaptive control (tutorial)

- **Claim.** What each leak fixes — the prior's atoms, its masses — and what it costs: the window σ₀ ∈ (κ_R, σ_crit), the floor, the bias,
  the transit, the rate, the a-priori bound.
- **Rests on.** Remark D‴.1 (the window), Remark D‴.10, Remark D‴.25 (the design table), Theorem D‴.23(vi) (the gate), Proposition D‴.27
  (the projected prior).
- **Missing.** Worked examples on a standard benchmark.
- **Venue.** IEEE Control Systems Magazine; Annual Reviews in Control.

### 11. Adaptive cogging-torque compensation in servos with transport leaks

- **Claim.** For a position servo with an unknown position-dependent disturbance (cogging torque), the sliding-variable controller with
  adaptive Gaussian kernels on the position reduces to the averaged slow loop of §5d along the position orbit; without an exact
  representation (the regime of Proposition D‴.28), the adaptive prior removes the finite-size floor of the learned compensation and a
  misplaced prior's bias, the OT leak removes the floor at fixed masses, the classical σ-modification keeps a large bias, and the real
  servo's position error follows the exactly filtered prediction (d/dt + λ)⁻¹(d/dt + k_s)⁻¹r to 0.3%.
- **Rests on.** Proposition D‴.29 ((i), (iv) exact; (ii)–(iii) formal), Propositions D‴.27 and D‴.28, `exp_servo_cogging.py`.
- **Novelty (audit).** The sliding-variable reduction and Gaussian networks in adaptive control are classical ([1]); new are the transport
  leaks' floors and priors on a standard servo problem, and the over-parametrisation of phase-plane kernels on circular orbits.
- **Missing.** A measured cogging profile (hardware or a high-fidelity simulation), measurement noise, input saturation, and
  velocity-dependent friction (a sharp Stribeck step, whose best approximant lies at infinity) with the projected prior in the real loop.
- **Venue.** IEEE/ASME Trans. Mechatronics; Control Engineering Practice.
