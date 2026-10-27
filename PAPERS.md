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
  exact, (c) proved for the linearised loop by Lemma D‴.26a, x_m ∈ C⁴; the basin open), Proposition D‴.27 ((i), (iii) proved, (ii)'s convergence proved
  at η = ∞ and formal at finite η), Proposition D‴.28 ((i), (iii) exact, (ii) as D‴.24(i); the dynamics numerical).
- **Novelty (audit).** D‴.23: standard per-mode formulas (proximal point, iterated Tikhonov, elastic averaging), new transfer to the
  mass-constrained particle system uniformly in m, new failure mode with its remedy; D‴.24: see paper 6, and a learned prior as the way to
  restore the interpolation property under regularisation; D‴.26: classical tools (averaging, the Magnus commutator, e-modification), a
  new structural statement (the exact m-independent spectrum at every representable sample, the even expansion in slow time, proved by
  Chang's decoupling with its periodic and stochastic precedents — Kokotović–Khalil–O'Reilly, Haddad–Kokotović, Yang–Zhu); D‴.27: classical
  ingredients (projection, the switching σ-modification, which it is at η = ∞), transferred to the learned centre of the mass-constrained
  prior at finite η; D‴.28: a classical symmetry (the conic view's homogeneity), a new identification of the mechanism and an exact margin
  formula.
- **Missing.** The constants of D‴.23(v); the nonlinear basin of D‴.26 (its O(k_e⁻⁴) remainder bounded by Lemma D‴.26a, the constants not explicit); the convergence of D‴.27 at
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
  adaptive Gaussian kernels on the position reduces (formally) to the averaged slow loop of §5d along the position orbit. On a cogging
  torque that four kernels cannot represent — the loss has a nondegenerate local minimum ϑ^b(q) (Proposition D‴.28's regime) and its
  infimum at the saturation — the adaptive prior removes the finite-size floor that a fixed prior's leak creates and, from an engineer's
  misplaced model, most of its bias, at the rate of the weakest excitation; the OT leak removes the floor at fixed masses; the real servo's
  position error follows the exactly filtered prediction (d/dt + λ)⁻¹(d/dt + k_s)⁻¹r to 0.3%. Noise-free, no leak at all fits as well:
  the paper's case for the leaks rests on their robustness to noise and disturbances, which it must exercise.
- **Rests on.** Proposition D‴.29 ((i), (iv) exact; (ii)–(iii) formal), Propositions D‴.27, D‴.28 and D‴.30, `exp_servo_cogging.py`.
- **Novelty (audit).** The sliding-variable reduction and Gaussian networks in adaptive control are classical ([1]), and so is the partial
  excitation of Gaussian networks along an orbit ([6]–[8]); new are the transport leaks' floors and priors on a standard servo problem,
  and the exact fibre of isotropic phase-plane kernels on a circular orbit.
- **Missing.** The classical baselines (σ- and e-modification toward a nominal model, projection, a Fourier regressor of the known
  harmonics, which represents τ_cog exactly); measurement noise and disturbances, where the leaks' robustness is the point (Propositions
  D‴.30–D‴.31 for first-order plants); a measured cogging profile (hardware or a high-fidelity simulation); input saturation;
  velocity-dependent friction with the projected prior in the real loop.
- **Venue.** IEEE/ASME Trans. Mechatronics; Control Engineering Practice.

### 12. Measurement noise in mean-field adaptive control: coherent floors, the washout prior and the desired regressor

- **Claim.** In kernel adaptive control with transport leaks, measurement noise is coherent across the kernels (a replicated
  configuration stays replicated, so the floor does not average out with the network's size — invariance, not attraction: above a noise
  level the replicas separate), the closed loop high-passes it (a function-level floor ∝ N₀/k_e in closed form, the prior's mobility η
  entering the Tikhonov filter factor as h_j + η), and its errors-in-variables part is a persistent force to which an adaptive prior
  answers with the inverse excitation spectrum. The adaptive prior is a washout leak — no affine time-invariant prior removes the leak's
  bias and keeps its bound on persistent forces (RH + B = I) — so the remedy acts on the force: with the regressor at the reference, the
  noise's own force becomes a least-squares one, whose response is bounded under every prior, and the trade-off in k_e disappears; the
  parameters' fluctuation adds a force that is not the learning force of a mean tracking error (the classical second-order bias of a
  constant-gain stochastic approximation): small in absolute terms on the instance (≤ 3.9·10⁻⁴ per (Eν²)²), −8.7% to +3.3% of DR's and
  DR*'s biases and +56% of DR's at k_e = B = 100, its own response falling like k_e⁻² on a one-dimensional reference, where the law's
  leading term averages out — on a curve in ℝⁿ it need not, which matters for the second-order plants. With the regressor at the
  reference the regressor's own term leaves the separations' noise drift, which becomes O(1/k_e) like their stiffness: the noise threshold
  is bounded below uniformly in k_e (constant when the noise's bandwidth scales with the gain), and on the instance the noise stiffens the
  separations near the lift — not along the whole trajectory (Corollary D‴.33).
- **Rests on.** Proposition D‴.30 ((i) exact as invariance, its threshold formal at the lift — along the noisy trajectory the separations'
  drift grows with the developing errors-in-variables deviation, checked on antithetic noise pairs, the long-time threshold open —; (ii)
  proved for the linearised loop under white noise as a period average (Lemma D‴.26a), (iii) formal, their closed forms exact for the
  averaged system; (ii) checked against the exact linearised loop, (iii)'s stationary response against its exact periodic mean, the
  transient in the real loop), Proposition D‴.31 ((i)–(iii) proved, (ii) with R and B from the priors' blocks; (iv)(a)–(b) proved for the
  linearised equations, (iv)(c) formal,
  the parameters' part numerical; the remedy checked in the real loop on common noise paths through MR − DR), Corollary D‴.33 ((i) and
  (iv)'s formula formal at D‴.30(i)'s order, (ii)'s bound proved, (iii)'s balance formal with the exact periodic mean; checked along the
  exact second-order transient to the stationary point and on common noise paths: MR's separations growing, DR's decaying at the
  transient's prediction),
  Corollary D‴.26,
  `exp_measurement_noise.py`, `exp_eiv_remedy.py`.
- **Novelty (audit).** Classical ingredients (the errors-in-variables bias, stochastic averaging and the diffusion approximation, the LMS
  misadjustment and leaky LMS, common noise, the lack-of-excitation drift, washout filters and low-frequency learning, desired
  compensation, the Tikhonov/Kalman split); new, for the mean-field transport-leak setting: the floor's interpolation in the prior's
  mobility, the driving matrix of the sensitivity-filtered measurement, the coherence of the K₀-atom identification with its noise
  threshold at the lift, the washout reading of the adaptive prior uniform in m, and the quantification; the parameters' own errors in
  variables is classical (Benveniste–Métivier–Priouret, Dieuleveut–Durmus–Bach, Box), new only the observation that the desired
  regressor does not remove it; and the desired regressor's effect on the coherence (D‴.33): in the classical frame of noise-induced
  synchronisation ([S61]), the removal of the regressor's term, whence a threshold bounded below uniformly in k_e.
- **Missing.** The second-order computation proved (the parameters' part in closed form; the first order is, for the linearised loop under
  white noise, Lemma D‴.26a); the long-time noise threshold of the separations (the second-order expansion fails where the mean
  deviation is O(1); under DR the drift turns indefinite along the trajectory, formally at noise levels outside the expansion); the stationary biases in the real loop (the weak modes equilibrate over ~10⁵ time units); a remedy for the
  parameters' errors in variables (Richardson–Romberg extrapolation in k_e or the gain, as for constant-step SGD, is a candidate); second-order
  plants (Proposition D‴.29) and a servo with a realistic sensor (quantisation, delay).
- **Venue.** IEEE Trans. Automatic Control or Automatica (a technical note); Systems & Control Letters.

### 13. A cohesion leak and mass transfer for mean-field adaptive networks

- **Claim.** A leak that pulls each kernel toward the mean shift of the network's own parameter distribution — the particles attract
  each other, no prior, no atom count, no plan — is the gradient flow of the fit plus an aggregation energy. With a finite-range kernel it
  acts exactly on the unexcited subspace at clustered configurations and has no bias there, so every critical point of the unregularised
  weighted loss is stationary, for every number of classes and every set of masses; it holds the separations at σ₀w_k and leaves the
  collective modes to the fit; a class can split iff σ₀w_k² < q_kλ⁺_k (an exact representation is never over-split), and a split sheds
  single particles (the kernel's concavity). A Gaussian kernel pays an exponentially small bias for merging twins at an explicit rate
  when no third class is within reach. In the real loop the lift of an exact representation is periodic with the separation exponents
  σ₀w_k/k_e and the η = ∞ adaptive prior's collective spectrum, and with the desired regressor (Corollary D‴.33) the classes' noise
  threshold is bounded below uniformly in the gain. The number of kernels is not selected by the leak: it is the transit's and the
  schedule's. A replicator law on the kernels' weights (Proposition D‴.34) — the Fisher–Rao part of a Wasserstein–Fisher–Rao gradient
  flow, in real time an exponentiated-gradient law with the relative entropy as its Lyapunov function — keeps the stationary set at clustered configurations under the finite range (under the Gaussian kernel it selects the masses) and the saturation's bound, and removes the transit's mass obstruction: at fixed positions bounded trajectories reach the exact fit iff the target is representable. A split is an S_n-symmetric transcritical bifurcation whose quadratic coefficient is the residual's third derivative along the splitting direction (Proposition D‴.35): the extreme particle of the jitter runs away first, alone or with a near-tied neighbour, at the reduced equation's time, and after a departure the stiffness's drop and the curvature's gauge rise favour more. Under dither the replicator's masses diffuse along the manifold of exact fits — Katzenberger's limit, at the predicted rate — and their amplitudes run into the saturation at no cost in fit; an entropic leak on the masses pins them, without bias where the target's particle fractions are representable, with the Ornstein–Uhlenbeck spread ∝ ε^{−1/2} for ε ≫ τ (Proposition D‴.36). Under dither the masses also drift, at the order of their diffusion: in the slow-replicator limit they follow the replicator flow of the entropic free energy (τ/2) log det of the positions' Hessian, toward the saturated, flattest fits; a boundary layer of width O(√γ) slows them but, the amplitude diffusing without drift inside it, does not stop them (Proposition D‴.37).
- **Rests on.** Proposition D‴.32 ((i)–(iii) proved; (iv) proved for the stationary points and their linearisation; (v) the twins' rate
  to leading order, the merger only without a third class within reach; (vi) the lift exact, the collapsed exponents by Corollary D‴.26
  and Lemma D‴.26a), Corollary D‴.33, Propositions D‴.34–D‴.37, `exp_cohesion_leak.py` (the linear theory, warm and broad starts, the real loop,
  dither, the cogging servo, the σ₀ sweep, other targets, the scale window, the counterexample to merging), `exp_mass_transfer.py`, `exp_split_bifurcation.py`, `exp_gauge_dither.py`, `exp_gauge_drift.py`.
- **Novelty (audit).** Mean shift, aggregation equations, bounded-confidence and polarized consensus-based clustering, clustering in
  self-attention dynamics, splitting steepest descent and Fukumizu–Amari's critical lines, convex clustering and sum-of-norms
  identification, deterministic annealing and self-organising RBF networks are classical ([S60]), as are Wasserstein–Fisher–Rao flows,
  conic particle gradient descent, neuron birth–death and exponentiated-gradient laws ([S62]); new is their combination as a leak in
  mean-field adaptive control: the exact blocks at clustered configurations with the mass factor, the gauge-scaled splitting margin, the
  twins' rate, the real loop's exponents, and the replicator's removal of the saturation trap with the representability criterion (for bounded trajectories).
- **Missing.** The transit from a spread cloud (which classes form, with which masses; from small masses Proposition D‴.48's transit forms
  one kernel per feature) and any selection of the number of kernels by the leak itself; the
  shedding of single particles below the resolution's threshold (the replicator does not remove it); a scale R chosen from the data (at a
  fixed scale the broad starts do not end at the target's count, Part J); the basin uniform in m; the a-priori bound under the transport
  (the leak is internal — a projection, Proposition D‴.27, is the remedy); the replicator's neutral directions under noise beyond papers 14–15 (the redistributions within a class), boundedness at fixed positions, and the real loop with transport and replicator together.
- **Venue.** SIAM J. Control Optim. or IEEE TAC (theory); NeurIPS/ICML-style venues for the mean-field network regulariser and the mass
  transfer.

### 14. Noise-induced selection on the gauge of a kernel network: Katzenberger's drift, its Fixman correction, general noise and the dither's confinement

- **Claim.** With a replicator (Fisher–Rao) law on the kernels' weights and dither on their positions, the weights are a gauge of the
  exact fit: the fits through a configuration form a manifold parametrised by the masses, on which the function does not move. Dither
  makes the masses diffuse along it — Katzenberger's limit, the noise projected in the flow's own metric, at the predicted rate
  (Proposition D‴.36) — and drift at the same order: in the slow-replicator limit they follow the replicator flow of the entropic free
  energy (τ/2) log det of the positions' Hessian, which with one particle per class has the closed form 2Σ log(W_k a′_k) and diverges at
  every saturation, so the noise drives the masses into the saturated, flattest fits; a boundary layer of width O(√γ), where the
  replicator outruns the saturating amplitude, slows but does not stop them (Proposition D‴.37). A Fixman correction — the free energy's
  derivative added to the replicator's fitness with a gain κ_F, in closed form on the masses — makes the corrected loop a gradient flow
  with a log-barrier on the representable masses: it cancels the drift to leading order (κ_F = 1), or reverses it toward the barrier's
  maximiser, the unique D-optimal exact fit (κ_F > 1, robust to an underestimated dither level; in the full loop within 1.1 standard
  errors of it), and, unfloored, makes the saturation non-attracting, the density there δ^{κ_F − 1} (a floor that keeps the simulation
  finite restores δ⁻¹ within its cap); with a leak the pinned state is predicted to second order, the fluctuations acting through the
  leak's curvature; evaluated at the current amplitudes instead the correction turns the positions' excursions into a mass runaway
  (Proposition D‴.38). Under the cohesion leak, which has no prior, the energy is bounded on the unconstrained coordinates, so the
  dithered loop has no stationary law: the positions leave the exact fit's well, after a time of order e^{ΔJ/τ}, for a plateau of
  approximate fits with saturated amplitudes; a dead zone on the unconstrained coordinates (the Moreau–Yosida envelope of a box's
  indicator) restores the stationary law — without bias for the exact fits inside its box, though the law still weighs the plateau where
  it meets the box —, turns the gauge's saturation into a non-degenerate wall at which the entropic drift ends, and keeps the Fixman
  correction's selection (Proposition D‴.39; the escapes at fixed masses with an Arrhenius slope of (0.81 ± 0.14) of the barrier at the
  instance's moderate barriers, the dead zone holding the plateau's occupation at 0.09–0.19). With several particles per class, together
  at the class's exact-fit position, the positions' Hessian factors into the class masses' and the cohesion's within-class block: under
  the isotropic dither the selection, its correction and a frozen split act on the class masses as with one particle per class, slowed
  by the number of particles, while an anisotropic dither makes the split felt and moves the class masses at unequal splits
  (Proposition D‴.42). In preparation (Proposition D‴.40, in review): for any noise covariance on the positions the gauge identity makes
  the masses' fitness the generator's image of a function of one amplitude, so that at leading order only each amplitude's own noise
  intensity selects — the replicator flow of a weighted log-determinant Σ_kτ_k log(W_k a′_k) whatever the noise —, dithering only the
  centres and widths freezes the masses, and the correction at the classes' levels, read from the amplitudes' quadratic variation,
  cancels the selection.
- **Rests on.** Propositions D‴.36 ((i) proved for one particle per class, (ii) by Katzenberger's theorem on compact parts of M, (iii)
  proved with the O(τ/ε) shift sketched), D‴.37 ((i) proved given Katzenberger, (ii) the averaging and Laplace steps sketched and checked
  to 0.3%, (iii) the divergence and the layer's scale, the non-normalisability; (iv) first order), D‴.38 ((i) the slow-limit cancellation
  proved, the remainder numerical; (ii) proved; (iii) the projected drift proved on M, the exponents checked in Katzenberger's limit;
  (iv) checked on paired seeds; (v) numerical), D‴.39 ((i) proved; (ii) the plateau numerical, the exit time's order Freidlin–Wentzell's;
  (iii) (a), (b) proved ((b) for the deterministic objects), (c) sketched at leading order), D‴.42 ((i), (iii) proved, (ii) proved for the
  isotropic dither, its anisotropic failure's mechanism proved and its size numerical); `src/core/gauge_drift.py`,
  `exp_gauge_dither.py`, `exp_gauge_drift.py`, `exp_gauge_fixman.py`, `exp_gauge_plateau.py`, `exp_gauge_classes.py`; Proposition
  D‴.40 and `exp_gauge_noise.py` once their review and runs are complete.
- **Novelty (audit).** Diffusions forced onto a manifold of equilibria (Katzenberger), the noise-driven drift along manifolds of minima
  and its log-pseudo-determinant form for isotropic noise (Blanc et al., Li–Wang–Arora), averaging and the potential of mean force
  (Khasminskii, Kirkwood), flat minima (Hochreiter–Schmidhuber, Entropy-SGD), the Fixman correction of stiff constraints (Fixman, Hinch,
  Ciccotti–Lelièvre–Vanden-Eijnden), the adaptive biasing force and well-tempered metadynamics, the max-det selection of a log-barrier,
  dissipativity, metastable exit times and parameter projection (Mattingly–Stuart–Higham, Raginsky–Rakhlin–Telgarsky, Freidlin–Wentzell,
  Bovier et al., Pomet–Praly), Langevin sampling with a Moreau–Yosida envelope (MYULA), Jacobi coordinates are classical ([S64]–[S67]).
  New: the amplitude–mass gauge of a replicator-weighted kernel network as the slow manifold, with noise on the positions only and a slow
  mass channel; the log-determinant's closed form; the consequence that the noise brings back the saturation trap the replicator had
  removed, and the √γ layer that does not stop it; the online Fixman correction on the masses, its gain selecting the D-optimal fit; the
  loss of the dithered loop's stationary law with the prior, and the dead zone that restores it without touching the selection; the class
  masses' laws with several particles per class exactly when the dither is isotropic.
- **Missing.** Proposition D‴.40's review and runs (the general noise); separated particles of one class (the split's onset under
  dither); the fractions' order-γ²τ drift; the real closed loop (the masses' time scale 1/(γτ) is out of direct simulation's reach on the
  instance); an estimate of τ when the noise is measurement noise (Proposition D‴.41 for the error's noise); the averaging step of
  D‴.37(ii) with constants; the plateau's occupation in the box's stationary law and the box's term at small τ.
- **Venue.** SIAM J. Applied Dynamical Systems or Multiscale Model. Simul. (the limit, its correction and the confinement); a
  learning-theory venue (COLT/ALT) for the implicit bias of noisy replicator dynamics.

### 15. The colour of label noise is a sparsity prior: no selection in continuous time, the trace for white noise, the log-determinant for a slow Markov one

- **Claim.** For a gradient flow driven by an error signal — supervised learning on the square loss, an adaptive control law, the
  replicator–transport loop of an adaptive approximator — noise in the error enters along the descent's own vector fields, which are
  orthogonal to the manifold of exact fits in the flow's metric. The limit map's invariance then makes the bilinear form
  ⟨∇r ⊗ ∇(∂Φ·U)⟩ antisymmetric, and the Wong–Zakai term cancels Katzenberger's Itô drift exactly: in continuous time, for a noise white on
  the data's resolution and of any temporal spectrum, the parameters neither drift nor diffuse along the manifold at leading order — a
  loop adapting in continuous time does not select through its measurement noise. Only a discrete update selects, through the second
  moment of the whole increment: at linear response each pair of normal modes is weighted by a filter f_il of the step and the noise's
  correlation (the cross term the correlation adds cancelling by the same antisymmetry). On a flat-metric network the selection descends
  tr Ψ(∇²L), Ψ′ = f: for fresh draws Ψ = −(2/h)log(1 − hλ/2), Damian–Ma–Lee's finite-step regulariser; for a finely sampled Markov noise
  of correlation time τ_c Ψ = log(1 + λτ_c)/τ_c — (1/τ_c)log det(I + τ_c∇²L) —; for every coloured noise in the continuum Ψ is concave.
  On factorised coefficients — w = u⊙u − v⊙v, the quadratic parametrisation of a diagonal network or of a kernel expansion's two adaptive
  masses per atom — white label noise selects, as the step and the noise vanish, the weighted ℓ¹ norm (basis pursuit), a slow Markov noise
  the log-determinant of Fazel's heuristic (sparse Bayesian learning's log-det term, not its penalty), which by Cauchy–Binet counts the
  support: a k-sparse exact fit is held once τ_c passes a coherence-dependent threshold (a weighted dual certificate), and the slow
  selection recovers sparse expansions where weighted basis pursuit fails (on 10 × 40 Gaussian designs 0.75 and 0.60 of the 3- and
  4-sparse truths against basis pursuit's 0.38 and 0.17; label-noise descent with correlated draws at small σ where its reduced flow
  does) — the mechanism, not a better sparse solver (sparse Bayesian learning's EM recovers 0.88 and 0.55) —; the coloured selection is
  multistable, the concave penalty's local minima — on a replicator's masses, three stable targets where white noise has one.
- **Rests on.** Propositions D‴.41 ((i) proved; (ii) proved at order v given Katzenberger's and the Wong–Zakai theorems; (iii) proved at
  leading order; (iv) proved at order v and linear response for any step h < 2/λ_max and correlation, (c)'s concavity in the continuum
  limit, (e)'s finite-step slow limit; the examples exact per step) and D‴.43 ((i) proved at order v given D‴.41(iv) where rank K = N,
  the pruning phase numerical; (ii) the potential proved, the end at basis pursuit's solution sketched; (iii), (iv) proved; (v) proved for
  the reduced flow, the certificate sufficient; the recovery rates numerical); `src/core/label_noise.py`, `src/core/colour_selection.py`,
  `exp_gauge_label.py`, `exp_colour_selection.py`.
- **Novelty (audit).** Katzenberger's limit along manifolds of minimisers (Li–Wang–Arora), label-noise SGD's flat-minimum and sparsity
  selections and the finite-step regulariser (Blanc et al., Damian–Ma–Lee, HaoChen et al., Pillaud-Vivien–Reygner–Flammarion,
  Shalova–Schlichting–Peletier), the Itô–Stratonovich distinction and the coloured-noise crossover (Wong–Zakai, van Kampen,
  Kupferman–Pavliotis–Stuart), the active Ornstein–Uhlenbeck filter and the coloured-noise approximation's determinant (Jung–Hänggi,
  Maggi et al.), correlated and anti-correlated noise (Orvieto et al., Lucchi et al., Kühn–Rosenow, Beneventano), the noise equilibria of
  parameter symmetries (Ziyin et al.), sparse Bayesian learning, its non-separable reweighting and the log-det heuristic (Tipping,
  Wipf–Rao, Wipf–Nagarajan, Fazel–Hindi–Boyd), reweighted ℓ¹ (Candès–Wakin–Boyd), stochastic collapse (Chen et al.) and the dual
  certificate (Fuchs, Tropp) are classical ([S69], [S70], [S72]). New (not found stated): the exact continuous-time cancellation for noise
  entering through the error; the coloured filter of label noise's drift along the exact fits; that the noise's temporal correlation sets
  the sparsity prior of the explicit update's selection — the trace for white noise, a support-counting log-determinant for a slow Markov
  one, a concave penalty for every colour —, its stability law and its recovery beyond basis pursuit; the closed forms and the
  multistability on the replicator's masses.
- **Missing.** The order-v² motion along the manifold; the singular limit near sparse exact fits (the normal space degenerating, the
  pruning phase); minibatch noise; trained networks with temporally correlated label noise (augmentation repeated over epochs, sensor noise
  in online learning) at scale; the paradigm's own replicator mass law on a dictionary (Proposition D‴.44, in preparation); a real closed
  loop with sampled-data adaptation.
- **Venue.** A learning-theory venue (COLT, NeurIPS) or SIAM J. Math. Data Sci. for the filter and the sparsity prior; IEEE TAC or
  Automatica for the adaptive-control reading (continuous-time adaptation does not select through its measurement noise; a sampled one
  selects by its noise's spectrum, and a slow sensor noise prunes the kernels).

### 16. The free mass–amplitude law of self-adaptive kernels: a mirror flow, its implicit bias for every start, a leak in its mirror coordinates

- **Claim.** Free the kernels' masses from their constraint and each kernel's mass and saturating amplitude are driven by one scalar
  signal: Q = β − (γa_max²/4G)cosh(2α/a_max) is conserved for every drive, so the law only reparametrises the amplitude and its coefficient
  w = βa/m follows the mirror flow of an explicit separable potential — a mass-weighted ℓ² core, a soft wall at the saturation, an entropy
  tail; for small starting masses an ℓ¹ norm of weight ∝ β⁰^{−1/2} (Proposition D‴.45). On a dictionary the loop converges from every
  start to the Bregman projection of the start onto the sign-constrained least-squares fits — the exact fits when the frozen signs allow
  it —, the kernels whose sign cannot help dying out; starts at nonzero amplitude and small mass put the kernels on branches of their
  curves, where the signs freeze and the law is an exponentiated gradient with the temperature 1/(γa₀) (Proposition D‴.47): the starting
  masses and amplitudes select between the minimum-norm fit, weighted basis pursuit and sign-constrained basis pursuit, and on positive
  kernel dictionaries the frozen signs make the target's kernels the limit at every scale. Noise with smooth paths (Stratonovich) never
  moves the selected fit; a sampled-data update does, through a closed-form change of Q per step, which a mirror update on the curves
  removes. In the real loop ½e² plus the Bregman divergences is a Lyapunov function of the joint mass–amplitude adaptation; a leak written
  in the mirror coordinates keeps it under disturbances (an ultimate bound the free masses alone lack) and makes the loop's bias an
  explicit Bregman-regularised fit — a ridge from large starting masses, a LASSO from small ones that keeps exactly the target's kernels
  under the irrepresentable condition, the averaged loop's equilibrium written with the error filter's Gram (Proposition D‴.46).
- **Rests on.** Propositions D‴.45 ((i), (ii), (iv)(a), (v) exact; (iii) proved for Q > −c, the rich limit when weighted basis pursuit is
  unique; (iv)(b) the formula proved, the selection's end numerical), D‴.46 ((i) proved; (ii) for the averaged loop, the real loop's mean
  by classical averaging; (iii) the limits proved, the support recovery by the LASSO's theory), D‴.47 ((ii) for every start; (iii) the
  limit and the off-support rate; (v) with a certified uniqueness); `src/core/mirror_free.py`, `exp_mirror_free.py`, `exp_mirror_leak.py`,
  `exp_mirror_branch.py`; the sampled-data selection's ends under the paradigm's update (Proposition D‴.44, `exp_mass_dictionary.py`) once
  its runs are complete.
- **Novelty (audit).** Mirror flows of commuting reparametrisations and their conservation laws (Li–Wang–Lee–Arora; Marcotte–Gribonval–
  Peyré; Azulay et al.), the kernel and rich regimes and the hyperbolic entropy (Woodworth et al.; Ghai–Hazan–Singer), exponentiated
  gradient and Hessian–Riemannian flows on constrained domains (Kivinen–Warmuth; Alvarez–Bolte–Brahic), nonnegative least squares by
  overparametrisation (Chou–Maly–Verdun et al.), SGD's Itô noise and the effective initialisation (Pesme et al.; Even et al.; HaoChen et
  al.), mirror-descent adaptation laws (Boffi–Slotine), discounted mirror descent and penalty-regulated dynamics (Gao–Pavel; Romano–Pavel;
  Coucheney–Gaujal–Mertikopoulos), the σ-modification and averaging of adaptive systems (Ioannou–Kokotović; Anderson et al.; Sastry–
  Bodson) and the LASSO's irrepresentable condition (Zhao–Yu; Wainwright) are classical ([S71]–[S74]); new is the paradigm's saturating
  unit made explicit — the conserved Q, the mirror coordinate and potential in every regime, the weights 1/√β⁰, the branches with their
  temperature —, the implicit bias for every start with the extinction of wrong-signed kernels, the sampled-data change of Q and the
  mirror update that removes it, and the dual-coordinate leak realised on (α, β) with its bound and its starting-mass-dependent LASSO.
- **Missing.** Moving centres and widths (Q stays conserved while the cohesion's α-component vanishes; the mirror structure does not
  extend; the small-mass transit with moving centres is Proposition D‴.48's); the rich equilibria's slow approach (times ∝ 1/β⁰² without, ∝
  β⁰^{−1/2} with the leak) and the extinction's rate without strict complementarity; the real loop's bias beyond averaging; the
  sampled-data drift of the projection in closed form.
- **Venue.** IEEE TAC or SIAM J. Control Optim. (the adaptive law, its Lyapunov function and the leak); COLT/NeurIPS-style venues for
  the implicit bias of the saturating unit and its branches.
