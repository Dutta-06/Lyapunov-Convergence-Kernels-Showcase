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
  flow, in real time an exponentiated-gradient law with the relative entropy as its Lyapunov function — keeps the stationary set at clustered configurations under the finite range (under the Gaussian kernel it selects the masses) and the saturation's bound, and removes the transit's mass obstruction: at fixed positions bounded trajectories reach the exact fit iff the target is representable. A split is an S_n-symmetric transcritical bifurcation whose quadratic coefficient is the residual's third derivative along the splitting direction (Proposition D‴.35): the extreme particle of the jitter runs away first, alone or with a near-tied neighbour, at the reduced equation's time, and after a departure the stiffness's drop and the curvature's gauge rise favour more. Under dither the replicator's masses diffuse along the manifold of exact fits — Katzenberger's limit, at the predicted rate — and their amplitudes run into the saturation at no cost in fit; an entropic leak on the masses pins them, without bias where the target's particle fractions are representable, with the Ornstein–Uhlenbeck spread ∝ ε^{−1/2} (Proposition D‴.36).
- **Rests on.** Proposition D‴.32 ((i)–(iii) proved; (iv) proved for the stationary points and their linearisation; (v) the twins' rate
  to leading order, the merger only without a third class within reach; (vi) the lift exact, the collapsed exponents by Corollary D‴.26
  and Lemma D‴.26a), Corollary D‴.33, Propositions D‴.34–D‴.36, `exp_cohesion_leak.py` (the linear theory, warm and broad starts, the real loop,
  dither, the cogging servo, the σ₀ sweep, other targets, the scale window, the counterexample to merging), `exp_mass_transfer.py`, `exp_split_bifurcation.py`, `exp_gauge_dither.py`.
- **Novelty (audit).** Mean shift, aggregation equations, bounded-confidence and polarized consensus-based clustering, clustering in
  self-attention dynamics, splitting steepest descent and Fukumizu–Amari's critical lines, convex clustering and sum-of-norms
  identification, deterministic annealing and self-organising RBF networks are classical ([S60]), as are Wasserstein–Fisher–Rao flows,
  conic particle gradient descent, neuron birth–death and exponentiated-gradient laws ([S62]); new is their combination as a leak in
  mean-field adaptive control: the exact blocks at clustered configurations with the mass factor, the gauge-scaled splitting margin, the
  twins' rate, the real loop's exponents, and the replicator's removal of the saturation trap with the representability criterion (for bounded trajectories).
- **Missing.** The transit from a spread cloud (which classes form, with which masses) and any selection of the number of kernels; the
  shedding of single particles below the resolution's threshold (the replicator does not remove it); a scale R chosen from the data (at a
  fixed scale the broad starts do not end at the target's count, Part J); the basin uniform in m; the a-priori bound under the transport
  (the leak is internal — a projection, Proposition D‴.27, is the remedy); the replicator's neutral directions under noise (the gauge and the redistributions), boundedness at fixed positions, and the real loop with transport and replicator together.
- **Venue.** SIAM J. Control Optim. or IEEE TAC (theory); NeurIPS/ICML-style venues for the mean-field network regulariser and the mass
  transfer.
