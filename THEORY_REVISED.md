# Revised theory: projection-based adaptive RBF control with adaptive centres and widths

This document replaces the stability argument of version 1 of `lyapunov.pdf` ("Mathematical
Foundations of Self-Adaptive Kernels Achieving Lyapunov Convergence") with a
statement whose assumptions are explicit, whose algorithm is exactly the one in
`src/core/`, and whose every step is either proved here or verified numerically
by `tests/test_theorem_to_code.py`.  Section 6 lists the errata to the
manuscript; Section 7 positions the result against prior work.

Notation: `‖·‖` is the Euclidean norm for vectors and the Frobenius norm for
matrices; `vec(·)` is the row-major vectorisation used throughout the code;
`m` kernels, `n` states.

---

## 1. Setting and standing assumptions

Plant, reference and error:

    ẋ = f(x, t) + u,      e = x − x_m(t),      u = ẋ_m − k_e e − f̂(x),
    f̂(x) = Wᵀ φ(x, c, σ),   φ_i(x) = exp(−‖x − c_i‖² / 2σ_i²).

**A1 (operation on a compact set).** There is a compact set X ⊂ ℝⁿ with
x(t) ∈ X for all t ≥ 0, and F_max := sup_{x∈X, t≥0} ‖f(x, t)‖ < ∞.
*(The manuscript assumes "a compact spatial domain X" in the same way.  For
globally bounded f the assumption is discharged by Theorem 1(iv); for the
Lorenz plant it is a genuine assumption.)*

**A2 (reference).** x_m is C¹ with x_m, ẋ_m bounded, and ẋ_m is available to the controller.

**A3 (admissible ideal parameters).** There is a curve θ*(t) = (W*(t), c*(t), σ*(t))
lying for all t in the *inner* projection sets

    ‖W*‖ ≤ W_max − ε,     ‖c*‖ ≤ c_max − ε,     σ_min + ε ≤ σ*_i ≤ σ_max − ε.

**A4 (ideal-parameter regularity).** θ*(·) is C¹ with
‖Ẇ*‖ ≤ V_W, ‖ċ*‖ ≤ V_c, ‖σ̇*‖ ≤ V_σ.  (For a time-invariant plant one may take θ* constant, V_W = V_c = V_σ = 0.)

**A5 (gains and geometry).** k_e > 0; Γ_W, Γ_c symmetric positive definite on ℝ^{mn}
(a zero Γ freezes that group, see Remark 3); Γ_σ diagonal positive;
0 < ε < min(W_max, c_max, (σ_max − σ_min)/2); σ_min > 0.

Remarks on the assumptions.

* **The proof never uses optimality of θ*.**  θ* is an analytic device; any
  admissible C¹ curve works, and the approximation quality of the network
  enters the result only through the *actual* size of the disturbance d(t)
  defined below, not through the bound.  The manuscript's Tikhonov-regularised
  definition of θ* and its implicit-function-theorem (IFT) argument for A4 are
  therefore unnecessary — and, as written, not a proof (Section 6, E3).
* **No centre-separation constraint is needed.**  The set
  {c : ‖c_i − c_j‖ ≥ δ ∀ i ≠ j} is not convex, so no projection onto it has
  property (P) below; the manuscript's claim that separation is "maintained
  via the projection set construction" is withdrawn (Section 6, E1).  The
  bounded-Jacobian property that the proof needs follows from σ ≥ σ_min alone.

---

## 2. The algorithm (identical to `src/core/controller.py`)

    Ẇ = Proj_W(W, Γ_W φ eᵀ)
    ċ = Proj_c(c, Γ_c φ'_cᵀ W e)
    σ̇ = Proj_σ(σ, Γ_σ φ'_σᵀ W e)

with φ'_c = ∂φ/∂vec(c) ∈ ℝ^{m×mn} (row i has the single non-zero block
φ_i (x − c_i)ᵀ/σ_i²) and φ'_σ = diag(φ_i ‖x − c_i‖²/σ_i³).

Norm-ball projection (W and c), boundary layer of width ε:

    Proj(θ, Γy) = Γy − h(θ) · Γθ · (θᵀΓy)/(θᵀΓθ)    if h(θ) > 0 and θᵀΓy > 0,
                = Γy                                 otherwise,
    h(θ) = 0 for ‖θ‖ ≤ θ_max − ε,  ((‖θ‖ − θ_max + ε)/ε)² in between,  1 for ‖θ‖ ≥ θ_max.

Box projection (σ, component-wise): the update component is multiplied by
clip((σ_i − σ_min)/ε, 0, 1) when it points below the lower face and by
clip((σ_max − σ_i)/ε, 0, 1) when it points above the upper face.

Discrete time: the simulations integrate with fixed-step RK4 and re-clip the
parameters to their sets after every stage (`integrator.py`).  The theorem is
a continuous-time statement; the clip is a safety net that is inactive if the
step is small enough (property (P) already removes the outward component on
the boundary).

---

## 3. Lemmas

**Lemma 1 (projection property (P)).** Under A3 and A5, for every θ in its set and every y,

    θ̃ᵀ Γ⁻¹ Proj(θ, Γy) ≤ θ̃ᵀ y,        θ̃ := θ − θ*.

*Proof.* If the projection is inactive the two sides are equal.  Otherwise
θ̃ᵀΓ⁻¹(Proj − Γy) = −h (θᵀΓy)/(θᵀΓθ) · θ̃ᵀθ, with h > 0 and θᵀΓy > 0, and
θ̃ᵀθ = ‖θ‖² − θ*ᵀθ ≥ ‖θ‖(‖θ‖ − ‖θ*‖) ≥ 0 because h > 0 forces ‖θ‖ ≥ θ_max − ε ≥ ‖θ*‖.
For the box, component i: σ̃_i (s_i − 1) y_i ≤ 0 because s_i ≤ 1, and on the
lower face y_i < 0 while σ̃_i = σ_i − σ*_i ≤ σ_min + ε − σ*_i ≤ 0 (symmetrically
on the upper face).  ∎

Numerical check: `tests/test_projection.py` (randomised, isotropic and
non-isotropic Γ, vector/matrix/box) and the counterexample showing that (P)
fails when θ* is outside the inner set.

**Lemma 2 (lumped disturbance).** Define

    d := f − f̂ + W̃ᵀφ + Wᵀφ'_c c̃ + Wᵀφ'_σ σ̃ .

Then identically  ė = −k_e e − W̃ᵀφ − Wᵀφ'_c c̃ − Wᵀφ'_σ σ̃ + d,  and under A1, A3,

    ‖d‖ ≤ D_max := F_max + 3 W_max √m + W_max ( 2 c_max J_c + √m (σ_max − σ_min) J_σ ),
    J_c = e^{−1/2}/σ_min,    J_σ = 2e^{−1}/σ_min .

*Proof.* The identity is the definition of d substituted into ė = −k_e e + f − f̂.
For the bound: ‖f‖ ≤ F_max; ‖f̂‖ ≤ ‖W‖‖φ‖ ≤ W_max √m (0 < φ_i ≤ 1); ‖W̃‖ ≤ 2W_max;
‖c̃‖ ≤ 2c_max; ‖σ̃‖ ≤ √m (σ_max − σ_min); rows of φ'_c have disjoint supports so
‖φ'_c‖₂ = max_i φ_i‖x − c_i‖/σ_i² ≤ sup_{r≥0} r e^{−r²/2σ²}/σ² = e^{−1/2}/σ ≤ J_c,
and ‖φ'_σ‖₂ = max_i φ_i‖x − c_i‖²/σ_i³ ≤ 2e^{−1}/σ ≤ J_σ.  ∎
(`tests/test_rbf_kernel.py` checks the two suprema by brute force.)

*Remark 1 (local form).*  Writing φ* = φ(x; c*, σ*) and expanding to first order
in (c̃, σ̃) gives d = ε + W*ᵀO_ho + W̃ᵀφ'_c c̃ + W̃ᵀφ'_σ σ̃ with ε := f − W*ᵀφ* and
O_ho the Taylor remainder — the manuscript's decomposition.  Near convergence
(θ̃ → 0) this shows ‖d‖ ≲ ‖ε‖ + O(‖θ̃‖²), which is why the controller performs far
better than D_max suggests; but it does not yield a smaller global constant
than Lemma 2, which needs no Hessian and no Taylor expansion.

**Lemma 3 (Lyapunov derivative).** Let

    V = ½ eᵀe + ½ vec(W̃)ᵀΓ_W⁻¹vec(W̃) + ½ c̃ᵀΓ_c⁻¹c̃ + ½ σ̃ᵀΓ_σ⁻¹σ̃ .

Along trajectories of the algorithm,

    V̇ = −k_e‖e‖² + eᵀd + cross − Ω*,
    cross = ⟨W̃, Γ_W⁻¹Ẇ − φeᵀ⟩ + c̃ᵀ(Γ_c⁻¹ċ − φ'_cᵀWe) + σ̃ᵀ(Γ_σ⁻¹σ̇ − φ'_σᵀWe) ≤ 0,
    Ω* = ⟨W̃, Γ_W⁻¹Ẇ*⟩ + c̃ᵀΓ_c⁻¹ċ* + σ̃ᵀΓ_σ⁻¹σ̇*,   |Ω*| ≤ C_Ω,
    C_Ω := 2W_max‖Γ_W⁻¹‖₂V_W + 2c_max‖Γ_c⁻¹‖₂V_c + √m(σ_max − σ_min)‖Γ_σ⁻¹‖₂V_σ .

*Proof.* Differentiate V, substitute Lemma 2's identity for ė, and use
eᵀW̃ᵀφ = ⟨W̃, φeᵀ⟩, eᵀWᵀφ'_c c̃ = c̃ᵀφ'_cᵀWe, eᵀWᵀφ'_σ σ̃ = σ̃ᵀφ'_σᵀWe.  cross ≤ 0 is
Lemma 1 applied to the three groups; the bound on Ω* is Cauchy–Schwarz with A3, A4.  ∎

The identity (to machine precision) and cross ≤ 0 (at every sample) are
verified along simulated trajectories in `tests/test_theorem_to_code.py`, on a
plant that is exactly representable with known constant θ*.

---

## 4. Theorem 1 (uniform ultimate boundedness)

Under A1–A5, with R_UUB := D_max/(2k_e) + √(D_max²/(4k_e²) + C_Ω/k_e) and
V_max,Θ := sup over the projection sets of the parameter part of V:

**(i) Descent outside a ball.**  V̇ ≤ −k_e‖e‖² + D_max‖e‖ + C_Ω, hence V̇ < 0 whenever ‖e‖ > R_UUB.

**(ii) Pointwise ultimate bound.**  Every solution satisfies
V(t) ≤ max{V(0), c_inf} with c_inf := ½R_UUB² + V_max,Θ, every level c > c_inf
is reached in finite time, and

    limsup_{t→∞} ‖e(t)‖ ≤ √( R_UUB² + 2 V_max,Θ ).

**(iii) RMS bound (the quantitatively meaningful one).**  For every T₁ < T₂,

    RMS_{[T₁,T₂]}(‖e‖) ≤ D_max/(2k_e) + √( D_max²/(4k_e²) + C_Ω/k_e + V(T₁)/(k_e (T₂ − T₁)) ),

so limsup_{T₂→∞} RMS_{[T₁,T₂]}(‖e‖) ≤ R_UUB — *without* the parameter-energy term of (ii).
The same inequality holds with D_max and C_Ω replaced by the trajectory-wise
suprema of ‖d(t)‖ and |Ω*(t)| on [T₁, T₂].

**(iv) Reference bound (adaptation-free).**  Under A1 and the hard parameter
bounds alone, ė = −k_e e + (f − f̂) with ‖f − f̂‖ ≤ F_max + W_max√m =: F̃, so

    ‖e(t)‖ ≤ e^{−k_e t}‖e(0)‖ + F̃/k_e      for all t, whatever the adaptation does.

*Proof.* (i) Lemma 3 with Cauchy–Schwarz on eᵀd.  (ii) If V > c_inf then, since the
parameter part of V is at most V_max,Θ, ½‖e‖² > ½R_UUB², so ‖e‖ > R_UUB and V̇ < 0; more
precisely, with Δ := V − c_inf > 0 one has ‖e‖² ≥ R_UUB² + 2Δ and therefore
V̇ ≤ −g(√(R_UUB² + 2Δ)) =: −ψ(Δ) with g(r) = k_e r² − D_max r − C_Ω increasing for
r ≥ R_UUB and ψ of class K.  The comparison lemma gives V(t) ≤ max{V(0), c_inf},
finite-time entry into every sublevel set {V ≤ c} with c > c_inf, and limsup V ≤ c_inf, i.e. the
stated bound on ‖e‖ (‖e‖² ≤ 2V).  (iii) Integrate (i): k_e∫‖e‖² ≤ V(T₁) − V(T₂) +
D_max∫‖e‖ + C_Ω(T₂ − T₁) ≤ V(T₁) + D_max √(T₂−T₁) √∫‖e‖² + C_Ω(T₂ − T₁) by
Cauchy–Schwarz; with S := RMS this reads k_e S² − D_max S − (C_Ω + V(T₁)/(T₂−T₁)) ≤ 0,
and S is below the positive root.  (iv) Variation of constants.  ∎

**What (ii)–(iv) say together.**  (ii) is the manuscript's conclusion; it is the
textbook level-set argument (Khalil, *Nonlinear Systems*, Thm 4.18; Lewis et al.
1999) and its radius contains the parameter-energy reservoir V_max,Θ — a bias-like
term — so it does **not** establish that projection yields a tighter tracking
bound than σ-modification.  (iii) is tighter and is the statement that can be
checked without unknown constants (replace D_max by the measured sup ‖d‖, see
`experiments/exp_oracle_ruub.py`: satisfied in 48/48 runs).  (iv) shows that mere
boundedness of e is already a consequence of projection plus a stable error
channel, independent of the adaptation law; the content of adaptive control is
in making the *actual* d(t) small, which no a-priori constant captures.  In the
oracle experiments the a-priori R_UUB exceeds the observed RMS error by factors
of 915× to 10 954×.

*Remark 2 (no LaSalle / Barbalat).*  Correct, and standard: UUB statements never
need them.  Nothing here claims asymptotic convergence of e or of θ̃.

*Remark 3 (frozen groups).*  If Γ_c = 0 (fixed centres) drop the c̃ term from V;
Lemma 2's identity still holds with −Wᵀφ'_c c̃ absorbed into d, with the same
bound.  Likewise for σ.

*Remark 4 (what is not proved).*  Parameter convergence (requires persistent
excitation of a nonlinear-in-parameter regressor — not addressed); exact
tracking; any statement about the size of the ideal approximation error ε.

---

## 5. The "no bias" claim, stated precisely

With σ-modification, V̇ ≤ −k_e‖e‖² − ½σ_mod‖θ̃‖²_{Γ} + ½σ_mod‖θ*‖²_{Γ} + D‖e‖ + C_Ω: the
leak adds a positive constant and removes the parameter-energy reservoir, which
gives exponential convergence of V to a ball whose radius grows with σ_mod‖θ*‖².
With projection, V̇ ≤ −k_e‖e‖² + D‖e‖ + C_Ω: no positive constant is added, but
nothing dissipates parameter energy, which is why the pointwise bound (ii) keeps
V_max,Θ.  Projection therefore removes the *Lyapunov-derivative* bias of leakage
(true, and known since Pomet–Praly 1992 / Ioannou–Sun 1996), not the ultimate
*tracking* bound's dependence on the parameter set (false as stated in the
manuscript's Sections 7.3 and 10.3).  Empirically the two mechanisms are
compared on identical geometry in `experiments/exp_bias_factorial.py`.

---

## 6. Errata to version 1 of `lyapunov.pdf`

| # | Section | Statement in the manuscript | Status |
|---|---|---|---|
| E1 | 1.2, 9 | Centre separation ‖c_i − c_j‖ ≥ δ "maintained via the projection set construction" with "convex projection boundaries" | **Wrong.** The separation set is non-convex; no projection onto it has property (P). The implemented repulsion term violated (P) in 30% of samples and did not keep the separation (min distance 0.07 vs δ = 0.5). Removed from the algorithm; not needed by the proof. |
| E2 | 1.3 | Offline Tikhonov regularisation of the ideal target: "This strict separation eliminates the circularity" | **Irrelevant to the proof.** Optimality of θ* is never used (Section 1 remark). |
| E3 | 1.3 | IFT gives finite velocity bounds V_W, V_c, V_σ | **Not a proof.** The Hessian of a Tikhonov-regularised non-convex loss at a minimiser is PSD, not necessarily PD; the IFT hypothesis is asserted, not established. Replaced by Assumption A4. |
| E4 | 2.1 | D_max from the Taylor remainder, "a finite geometric scalar" | **Valid but not explicit.** Replaced by the explicit, Hessian-free Lemma 2. The code's former `estimate_D_max = ε(1 + 2W_max c_max)` had no derivation and is deleted. |
| E5 | 4 | tr(W̃ᵀΓ_W⁻¹W̃) with Γ_W "positive-definite learning rate matrices" | **Type error** for Γ_W ∈ ℝ^{mn×mn}; the correct form is vec(W̃)ᵀΓ_W⁻¹vec(W̃). Harmless for scalar Γ. |
| E6 | 5 | Projection inequality "holds as an absolute mathematical certainty — the cancellation gap does not exist" | **True only under A3** (θ* in the inner set) and only for the Γ-weighted operator; the previous code applied the isotropic formula to any Γ. Fixed and tested. |
| E7 | 7.2 | "A critical and non-trivial consequence of the projection operator"; its hard energy cap "converts what would be a boundary-only argument … into a global interior argument" | **Standard** level-set / comparison argument (Khalil Thm 4.18). Correct, not new; also "strict monotone descent through every intermediate level set until the infimal set is penetrated" overstates: entry into every level above c_inf is finite-time, entry into c_inf itself only asymptotic. |
| E8 | abstract, 7.3, 10.3 | Projection gives UUB "without the persistent steady-state bias" (abstract); "no modification bias is introduced" (7.3); "zero bias into the Lyapunov derivative", "not a marginal improvement but a structural result" (10.3) | **Overclaimed.** The ultimate bound √(R_UUB² + 2V_max,Θ) contains the parameter-energy term; see Section 5. The tighter, bias-free statement is the RMS bound (iii). |
| E9 | 8 | "a standard Persistent Excitation (PE) condition on the regressor φ is insufficient to guarantee Θ̃(t) → 0"; convergence "requires alternative mechanisms", such as ICL | Parameter convergence is not treated; the remark is speculation and the code's former "Concurrent Learning" hook used only the current residual (not recorded data) and is removed. |
| E10 | README | "100% of the independent simulation trials converged strictly below the calculated R_UUB threshold" | **False.** The threshold 0.5 was a hand-chosen constant (`tier1_monte_carlo.py`: "A chosen arbitrary bound"). Replaced by `exp_monte_carlo_honest.py` / `exp_oracle_ruub.py`. |
| E11 | README | "O(mn) time" | **Was false** with the O(m²n) pairwise repulsion; true now that it is off by default. |

---

## 7. Prior art and the honest novelty statement

* Adaptive Gaussian-network control with robust modification: Sanner & Slotine, *IEEE TNN* 3(6), 1992.
* Adaptation of inner-layer parameters via first-order Taylor expansion of the feature map with bounded remainder — the exact device of Sections 2–3: Lewis, Yesildirek & Liu, *IEEE TNN* 7(2), 1996; Polycarpou, *IEEE TAC* 41(3), 1996; Ge & Wang, *IEEE TNN* 13(1), 2002; dynamic-structure RBF networks with tuned centres/widths: Fabri & Kadirkamanathan, *IEEE TNN* 7(5), 1996.
* Smooth projection with a boundary layer and the inequality (P), including the "no bias when θ* ∈ Ω" observation: Pomet & Praly, *IEEE TAC* 37(6), 1992; Ioannou & Sun, *Robust Adaptive Control*, 1996 (§4.4); Lavretsky & Gibson, arXiv:1112.4232, 2011; Lavretsky & Wise, *Robust and Adaptive Control*, 2013 (ch. 11).
* σ- and e-modification and their bias: Ioannou & Kokotovic 1983; Narendra & Annaswamy, *IEEE TAC* 32(2), 1987.
* UUB via level sets: Khalil, *Nonlinear Systems* (3rd ed.), Thm 4.18; Lewis, Jagannathan & Yesildirek, *Neural Network Control of Robot Manipulators and Nonlinear Systems*, 1999.
* Concurrent learning (recorded-data stacks, finite excitation): Chowdhary & Johnson, CDC 2010 / *J. Guid. Control Dyn.* 2011.

**What is new here.** Nothing in the stability theorem: Theorem 1(i)–(ii) is the
composition of known pieces (Taylor-linearised inner-layer adaptation + boundary-
layer projection + level-set UUB), applied to Gaussian RBFs with simultaneously
adaptive weights, centres and widths.  The contribution of this repository, after
revision, is (a) a theorem-to-code map in which the implemented and the proved
algorithm coincide and are machine-checked (Lemmas 1–3 verified numerically
along trajectories), (b) the explicit Hessian-free disturbance bound of Lemma 2
and the RMS bound (iii), which is the form of the result that can actually be
confronted with data, and (c) a controlled empirical separation of
"projection vs leakage" from "adaptive vs fixed geometry" (Section 5;
`results/factorial/`).  Claims of a "structural" advantage over leaky
modifications are withdrawn.
