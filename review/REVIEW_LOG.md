# Review log — PARADIGM.md

## Round 1 (adversarial, 2026-10-14)

**Design.**  Each of Theorems A–E got three independent refuters with distinct lenses: proof gaps,
counterexamples, and hypotheses/constants. A judge then adjudicated each theorem, and a cross-cutting
reviewer read the whole document. The raw reports are in `REVIEW_round1.md`.

**Coverage.**  Judges ran for A, B, D and E. The judge for C and the cross-cut reviewer did not run,
because the session hit its subagent limit. For C, all three refuter reports are available, and every
issue they agreed on was applied below.

**Verdicts.**  A: holds with fixes (1 major). B: holds with fixes (no major). C: per the refuters, holds
with fixes (2 major). D: holds with fixes (2 major). E: holds with fixes (no major). **No theorem was
broken.**

| id | sev. | issue | fix (where) |
|---|---|---|---|
| A.1 | major | A(i) said (L) "coincides" with the classical law; true only for raw rates, because the classical law projects jointly | A(i) restated (raw rates coincide; projected rates coincide iff per-kernel projection with W_i ≤ a_max/m and layer ε_K/m); scaling identity in the proof; new test `test_identification_projected_rates_per_kernel` checks it with active projections and checks sharpness |
| A.2 | minor | σ-modification row used the Euclidean norm, not the G⁻¹ norm | dictionary rewritten in the G⁻¹ metric |
| A.3–A.7 | minor | OT-leak gradient set-valued at ties; frozen-time reading; fixed-grid row unproved; support constraint only ε-regularised; wrong citation for push-forward; Nagumo invariance uncited | A(ii) and (iii) restated; fixed-grid row withdrawn; [46] cited; Nagumo cited |
| A.8–A.12 | cosmetic | scalar-gain hypothesis; status line cited missing files; garbled constants line | fixed |
| B.1 | minor | hypotheses used but not stated (ϑ_i(0) ∈ Θ_K, x ∈ X, regularity, C¹ ideal) | standing hypotheses (S1)–(S5) in §1 |
| B.2 | minor | x ∈ X discharged circularly | adaptation-free tube bound, Theorem D(b) |
| B.3 | minor | B.4 bound conditional, limits imprecise | Corollary B.5 marked conditional and a-posteriori; min form; both limits stated |
| B.4 | minor | "the labelled theory cannot state this" overclaimed; m-uniformity attributed to OT | Remark B.4: labelled theory can state it for a fixed labelling without switches; m-uniformity attributed to the normalisation (Remark A.1); §0 rewritten |
| B.5 | minor | continuum Remark B.6 cited AGS outside its hypotheses | now Remark B.7, a sketch with explicit hypotheses and a two-curve formula flagged as needing a reference |
| B.6 | minor | "globally the same order" unproved; OT constant is 3.9–5.3× the classical one | Remark B.4 gives both numbers and the block-structured tightening (533 → 382). **Erratum (round 2, AB.2):** the "3.9–5.3×" compared against a classical constant computed with a joint centre radius c_max. On the per-kernel set of Theorem A, the classical constant is 214 for every m, equal to the OT constant bounded per particle (214). OT is not smaller globally. |
| B.7 | minor | tests never exercised a moving ideal, ε ≠ 0, switches, or an active projection; `bound_B` omitted V_*W | `lyapunov_report` adds V_*(t)·W; new test `test_theorem_B_moving_ideal_nonrepresentable_switches_and_active_projection` |
| B.8–B.15 | cosmetic | "=" should be "≤" in step (c); e = 0 case; scalar Lagrange remainder; B.3 conclusions explicit; EVI wording; constants (closed forms, J definition, H_φ an upper bound); ε clash; status files | all fixed; `gaussian_constants` now uses closed forms; layer width renamed ε_K |
| C.1 | major | C(i) exactness argued by a false homothety claim | Lemma C.3: affine-in-λ argument, valid for any G ≻ 0, no switching; the refuters found the same proof independently |
| C.2 | major | C(v) false for V_OT (counterexample D⁺V_OT = +0.0475) | C(v) restated for V_π; Proposition C.5 (labelled leak reshuffles) with the counterexample; transient experiment |
| C.3 | minor | no solution concept for the discontinuous (L_σ) | Filippov/Carathéodory paragraph before Theorem C |
| C.4 | minor | b, Δ_gap time-dependent | defined as sup and inf over t |
| C.5 | minor | (iv) not m-uniform; Δ_gap = O(1/m); b = O(1/m) forced | "Scope" paragraph in Theorem C |
| C.6 | major (test) | C(iv) test bypassed the gain condition (`or True`) and violated both hypotheses | test rebuilt: an instance satisfying every hypothesis, asserted; envelope checked; contraction checked |
| C.7–C.10 | minor/cosmetic | endpoint of forward invariance; Z = 0 case; Remark C.2 citation and the universal claim; σ used for both a permutation and the leak gain | fixed; Remark C.6 cites AGS Thm 7.3.2 and Lemma 9.2.1, Santambrogio §4.4 and Carlen–Craig; universal claim withdrawn |
| D.1 | major | "(B) passes to the limit" unproved and ill-defined | Theorem D(e): explicit hypotheses on the ideals; proof by integrated form |
| D.2 | major | compact-set hypothesis discharged circularly | Theorem D(b): m-uniform, μ*-free tube bound ‖f − f̂‖ ≤ F_max + a_max with a first-exit argument |
| D.3 | minor | Lipschitz in x not uniform in t; counterexample verified | (S1) requires uniformity in t |
| D.4–D.6 | minor | constants undefined; Λ = O(1/ε_K) unstated; KR/W₁ bookkeeping; Dobrushin citation | constants defined in D(c); O(1/ε_K) stated; Golse [82] cited |
| D.7–D.11 | cosmetic | "smooth" projection; law (L) only; labelled-constant remark; verification row; status files | fixed |
| E.1 | minor | "projection inactive (since ϑ* ∈ Θ_K^ε)" false on the boundary of a closed set | hypothesis ϑ*_j ∈ int Θ_K^{ε}; `admissible_ideal(strict=True)` |
| E.2 | minor | Khalil Thm 4.13 hypothesis misquoted | proof step 2: C^∞ field with bounded second derivatives on D |
| E.3 | minor | N&A Thm 2.16 is the −uuᵀ theorem; change of variables missing | proof step 3: change of variables to the skew-symmetric form; Morgan–Narendra [S8] and Anderson cited |
| E.4 | minor | c₀, C, uniqueness of π, norm equivalence missing | explicit c₀, C and d²; uniqueness proved |
| E.5 | minor | convergence test qualitative, status overstated | verification map marks it qualitative; new test checks the proof's linearisation and Lyapunov identity |
| E.6–E.11 | cosmetic | symmetry misidentified (network vs ideal); labelled-limit norm; m-dependence of constants | fixed |

## Round 2 (adversarial, 2026-10-14)

**Design.** There were seven targets: H, G, C, D, E, A with B, and F. Each got two independent refuters and a judge.
The raw reports are in `REVIEW_round2.md`, and version 2 as reviewed is archived as `PARADIGM_v2.md`.

**Verdicts.** Every target held with fixes. Nothing was fatal.

| target | major | minor | cosmetic |
|---|---|---|---|
| H | 2 | 7 | 4 |
| G | 2 | 7 | 3 |
| C | 1 | 4 | 5 |
| A, B | 1 | 7 | 6 |
| E | 0 | 6 | 5 |
| D | 0 | 4 | 5 |
| F | 0 | 3 | 6 |

**Fixes, all applied in PARADIGM.md version 3.** The ids refer to `REVIEW_round2.md`.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| H.1 | major | the closed-loop-to-(𝔄) step sat inside a "proved" theorem; the cited test checks only an algebraic identity | Remark H.0, formal, with measured O(1/k_e) errors. New tests: n = 1 in `test_langevin.py`, n = 2 in `test_averaging.py`. §0 and §7 rows corrected |
| H.2 | major | Conjecture H.6 ill-posed; at fixed k_e it is false (stationary offset 0.1141 vs 0.1112; non-gradient correction) | H.6 restated in the scaling σ = σ₀/k_e, τ = τ₀/k_e, with the fixed-k_e counterexample. New Proposition H.7 (frozen-measure structure, non-reciprocity) and Proposition H.8 (ISS), and a route to H.6 |
| H.3–H.7 | minor | μ_0 ∈ P₂, σ, τ > 0, Dirac starts; autonomous f; H(iii) formal and missing continuity; H(iv) hypotheses; H(v) conditional | H(ii) restated, with the Girsanov regularisation of P₂ starts; autonomy assumed explicitly (and documented in `langevin.py`); H(iii) as a continuity bound with L_Φ plus a bias bound; H(iv) restated; H(v) moved to Remark H.5 |
| H.8 | minor | Remark H.2 overclaims; omitted τ rows; "slows learning" | log-rates e^{−8401} and e^{−1053}; all six τ rows; escape-threshold and censoring caveats; the (k_eσ, k_eτ) dependence. `exp_langevin_global.py` prints log(2τρ_τ), and both tables were regenerated |
| H.9 | minor | the titles "global optimality" and "costs nothing" | retitled; bias bound L̄(μ*_τ) ≤ inf_ν{L̄(ν) + k_eτH(ν\|γ)} added |
| H.10–H.13 | cosmetic | "stationary" should be "invariant"; tests weaker than claimed; Remarks H.3 and H.4; bookkeeping | preamble rewritten; symmetry test asserts exact row equality; new tests of the V_L bound and of the averaging error; Skorokhod reflection in H.3; H.4 rescoped; `expit`; τ guard; [S6] pointer |
| G.1 | major | (H2b) is not sufficient: slow-manifold slope | (H2c) S ≪ 1 added; the first-order induced-metric flow is G(b′); the counterexample is reproduced in `exp_slow_manifold_slope.py` → `table_slope.md` |
| G.2 | major | the measure-form/[S2] sentence was false | replaced by the "No mean-field form" paragraph; `LITERATURE_TRANSPORT.md` §G corrected |
| G.3–G.9 | minor | composition horizon; (H2a) mis-scaled; relation to H(i); (H0); the separation is a gain choice; the o(1) limit; Γ in (H2b) | Derivation marked formal, with the exact error identity; (H2a′) with ω_eff; G(a) reworded; (H0) added; Status rewritten; precise limit with slow time k_e/(mδ); Γ defined and used throughout |
| G.10–G.12 | cosmetic | L̄ monotone only once per period; explanation of the 0.94 ratio; table labels | reworded; centred differences (the ratio is now 0.97, and the 2.5% bias is explained); `exp_variable_projection.py` reports S, the joint fastest rate, ratios along the run, a direct G(a) check with a block split, and first/last labels |
| C.1 | major | solution concept: C(i) is false for Filippov sliding | selection solutions defined; Filippov version of (i) with the A_{ττ'} term; (ii)–(iv) hold for every Filippov solution |
| C.2–C.5 | minor | scope (Δ_gap = ∞, Δ_gap = 0 at switches); the R-free bound √C ≤ b + 2η; the (v) bias claim; "as Theorem F prices it" | Scope rewritten (m-uniform for a single-atom prior, O(m^{−1/2}) otherwise); (ii)–(iv) use 4η(b+η) < Δ_gap, and the test asserts √C ≤ b + 2η; (v) restated with b_lab(π) := sup_t; C.5 reworded |
| C.6–C.10 | cosmetic | example details; notation; hypotheses line; strictness; Remark C.6 wording | fixed |
| AB.1 | major | B.5 applied along (L_σ) via C(iv) is false (3.8–8.3×) | B.5 restricted to (L), with an (L_σ) variant carrying σ_m b²/4. New test `test_corollary_B5_fails_along_lsigma_and_its_leak_variant_holds` |
| AB.2–AB.8 | minor | B.4 numbers; A(i) "only if"; projection row not a gradient; σ-mod rows; Remark A.1; B.6 citation and rates | Remark B.4 table (214 = 214); two-case only-if proof, tested against `project_matrix`; dictionary rows corrected; Remark A.1 rewritten; B.6 now cites Chevallier [S13], Zador, and the second-order and Maurey rates |
| AB.9–AB.14 | cosmetic | (S0)/(S5) bookkeeping; ≤ Δ²/2; B.5 title; Kantorovich–Rubinstein; OT row; B.7 | fixed; B.7 via a coupling argument |
| D.1–D.4 | minor | (d) localisation; (b) as an inequality; (a) proof; L_1 vs L₁ | (d) proved by a bootstrap on X₁; (b) with S(F) ≤ F and examples; (a) by Picard on (x, Ξ); constants M_x and L_2 |
| D.5–D.9 | cosmetic | one-sided L_ϑ⁺; (e) wording; (c) proof; (S0); §0 row and numerics caveat | Remark D.1; (e) rewritten; (c) via D_δ; (S0) added; consistency-only caveat (e^{ΛT} ≥ e^{496}) |
| E.1–E.6 | minor | uniform tube and Khalil ball; constants; c₀ and uniqueness; Morgan–Narendra step; Corollary E.2; Remark E.1 | uniform tube ρ and interior distance δ; dependency list; c₀ := r₀²/(2(1+mλ_G)) with uniqueness via distinct equilibria; Morgan–Narendra Thm 1, then UAS ⇒ UES (Khalil 4.11), with the sin t² counterexample; E.2 is now exponential merging plus Hölder continuity; Remark E.1 rewritten with sliding-window values (`exp_pe_convergence.py`) |
| E.7–E.11 | cosmetic | coordinates z and w; PE wording; Remark A.1 conflict; k_e > 0; d² closed form | fixed; collapse point added as Remark E.3 |
| F.1 | minor | the contraction identity lacked a positive part and was not block-wise | Lemma P.1, cited in C(i) and F; new test `test_lemma_P1_projection_is_nonexpansive_and_tight_velocity_bound` |
| F.2–F.3 | minor | the Reading holds only along (L); the labelled-leak clause | Reading restricted to (L); d_lab clause; Proposition C.5 reworded |
| F.4–F.9 | cosmetic | test slack; box not C¹; hypotheses; leak versions; notation; block constant | tight-instance and leak-clause tests; docstrings in §1 and `projection.py`; (S0)–(S3); (i′)–(iii′); t₀, t₁, r and RMS defined; L_G and W_∞ |

**Errata to round 1.**

- Row B.6, the "3.9–5.3×" comparison, is corrected above.
- Round 1's fix C.1 and the proof of C(i) used the identity ‖Proj(ϑ,Gy)‖²_{G⁻¹} = ‖y‖²_G − h(2−h)(ϑᵀGy)²/(ϑᵀGϑ). It lacks a positive part and is not block-wise. The correct statement is Lemma P.1 (round 2, F.1). `REVIEW_round1.md` is left verbatim.

## Round 3 (adversarial, 2026-10-15): new and changed material of version 3

**Design.** Five independent refuters each took one block:

- the mathematics of Proposition H.7;
- the code and tests behind H.7 and Remark H.0;
- Theorem H as rewritten, Remark H.0, Conjecture H.6 and Proposition H.8;
- Proposition G as rewritten;
- a regression check of every round-2 fix in §1–§6.

A judge then adjudicated their reports. The raw reports are in `REVIEW_round3.md`.

**Verdict.** The judge's list has 51 confirmed issues: 1 major, 24 minor and 26 cosmetic. Nothing is fatal.

- **Rejected (2).** The claim that the review log lacked a round-2 section was stale. The H.0 gain-independence issue
  was downgraded from major to minor.
- **Uncertain (5).** These are the citations of Remark H.3, an alternative self-contained proof of Theorem E step 3,
  one numeric claim in the H.7(e) circulation proof, the concurrent edits, and two cosmetic G numbers.
- **Independent confirmation.** The refuters and the judge confirmed independently that the implementation of
  Proposition H.7 is correct. Fourteen formula mutations all fail the tests. The oscillation term matches an
  independent one-period flow to 0.7%, and every quoted number reproduces.

**Fixes, all applied.** The ids refer to `REVIEW_round3.md`.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| G1-offset | major | Proposition G(b′)'s slow-manifold offset used the VarPro velocity. It is 1 + S times too large, and the truncated flow ascends L* | (b′) offset is now (δ/g_a)M⁻¹Da*(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*, reducing to the old formula when S ≪ 1. Step 3 states which velocity to use |
| G2–G4 | minor | closed-loop equilibrium shift exceeds the offset; the satisfying-run ratio was misattributed; the slope runs included the boundary layer | caveat added to (b′). The 0.97 ratio is now attributed to the error-channel factor plus a +1.4% intra-period averaging correction. `exp_slow_manifold_slope.py` starts on the corrected slow manifold, adds a "gains ÷ 4" row, and the numerics give the error decomposition |
| G5–G8 | minor | L̄ monotonicity near equilibrium; error-channel bound; the gain condition counted Γ twice; numerics reporting | text replaced; the bandwidth term added to ω_eff; Status uses δ ≪ g_a·min(·); headings and k_e/ω_eff values corrected (2.6–4.8 and about 16) |
| G9–G11 | cosmetic | "first order"; R·S depends on Γ; undefined symbols | "(b′) is leading order in 1/R, for every S"; ε_avg, 𝒦 and R, S as suprema defined |
| H7a-global-uniqueness | minor | the periodic error is unique only locally (f(x) = x² gives a second, repelling orbit) | H.7(a) restated on balls of radius R, plus a global version under a one-sided Lipschitz condition; `averaging.py` docstrings |
| H7e-remainder-derivative | minor | (e) differentiated an O(k_e⁻³) remainder in μ | (e) restated for the k_e⁻² coefficient and proved by circulation around a triangle of measures (step 5 corrected after the judge: the limit is k_e → ∞) |
| H0-gains | minor | Remark H.0 said the first corrections do not depend on the gains, contradicting H.7a | H.0 now separates the frozen-measure corrections from b_osc, which is O(\|G\|T₀) relative; `exp_orbit_orientation.py` docstring |
| H7-area-claims, H7c-TD-analogy | minor | the kernel "is the area"; the TD analogy; the even non-reciprocity attributed to area | K is an integral of Σ_k det[∇Φ̃_k, ∇Φ̃'_k] over the enclosed region. The analogy is now performative prediction [S14], with TD as related. H.7(c) gives the performative-defect decomposition, which carries all of A_even. §0, §8 and the literature audit updated |
| H7a-3pct-claim, H7-status-and-checked-table | minor | "3%" was a size ratio only; tolerances were reported as measurements | 3% in size and 11% as a vector (cosine 0.994), with a vector-error column in the experiment; measured values (7e-15, 2e-15) separated from tolerances |
| H7-tests-weak | minor | loose thresholds; no full-matrix or value test | full-matrix second-variation check (5%), cyclic sums at 3%, a Richardson check for n = 1, per-doubling ratio checks, and the value test `test_oscillation_correction_matches_a_one_period_flow` |
| H0-period-averaging, H0-tikhonov-vs-evidence | minor | the averaging condition was unbounded through the σ-term; the runs sit below the Tikhonov scale | condition restricted to the loss-driven drift; Evidence limits added |
| H2-fit-saddle | minor | the "saddle, not barrier" reading is not supported per seed (27% and 33%), and the saddle is weak | reading withdrawn; `exp_langevin_tradeoff.py` reports per-seed fits (table regenerated from `raw.npz`) |
| H3-bakry-emery | minor | the 1/diam² LSI is not Bakry–Émery at σ = 0 | citations changed to F.-Y. Wang 1997 and Frieze–Kannan 1999 (marked uncertain by the judge) |
| Cv-strict-case | minor | "strict whenever the assignment changes" is false | a counterexample in which equality holds |
| E1-rate-not-bounded-by-PE | minor | no function of the PE level, k_e, T and sup‖Ψ_G‖ bounds the rate (it falls like ω⁻² under rescaling) | Remark E.1, the constants sentence, and a test comment |
| E2-coupling-constants | minor | Corollary E.2 used constants of f outside the tube | constants computed inside the Khalil ball |
| E-status, F-leak-test-power | minor | E conditional on [S8]'s abstract; the F leak test could not detect a missing leak term | §0 status; new test `test_theorem_F_leak_clause_is_tight_when_the_leak_pays_for_all_motion` (e ≡ 0, ratio ≥ 0.99) |
| H1-reference-values | cosmetic | numbers came from a review scratchpad | new `exp_gibbs_reference.py` → `results/langevin_global/gibbs_reference.md` (L̄(μ*_τ) = 0.11118; exactly averaged state 0.11415; function-space offsets) |
| others | cosmetic | notation clashes (γ, α, τ₀); definitions; Route bullets; the H.6 hedge; H.8 wording; tables; the A(iii) rows; D(b), D(c), C(i) and B.4 wording; the (S5) tube symbol; map rows; `__main__` blocks | all applied as given in `REVIEW_round3.md` |

**New results added while answering round 3** (and themselves reviewed by its judge):

- Remark H.7a: the second-order averaging term of the intra-period oscillation.
- The orbit-orientation experiment: the full first-order theory predicts the orientation dependence of the learned
  kernels in the real closed loop to 3% in size and 11% as a vector at k_e = 40.
- The reproducible slope counterexample, started on the slow manifold.
- The Gibbs reference values.

## Round 4 (manuscript conversion and review, 2026-10-18)

**Design.** The theory, the audit and the results were typeset as the manuscript `lyapunov.pdf` (version 2, source in
`manuscript/`). Twelve agents each converted one block. Eleven independent reviewers then compared each section with its
source, statement by statement and formula by formula, and checked the rendered pages. A twelfth checked the
bibliography against the published records. The reports are not archived; the confirmed issues are listed here.

**Verdict.** Two critical findings, both in the manuscript text and both fixed. No reviewer found a transcription error
in any theorem, proof, formula or table cell. The reviews also found the following errors in the sources themselves.

**Fixes to the sources.**

| where | issue | fix |
|---|---|---|
| PARADIGM §6c, Remark H.0 (period averaging) | "The sup is at most (T₀/k_e)λ_G(…)": the bound is for T₀ times the sup | "The left-hand side is at most …" |
| PARADIGM §2, proof of A(ii) | "ϑᵀProj ≤ 0 on the boundary" is wrong for the box factor at s = s_min (outward normal −1) | n_bᵀProj_b ≤ 0 with the outward normals n_b of Lemma P.1 |
| PARADIGM §5 numerics and §7 | "constant factor 1.68": `mean_field_m/table.md` gives 1.64–1.69 | "nearly constant factor (1.64–1.69)" |
| PARADIGM §5b | "W(μ_10, μ_0) = 1.05": the test asserts > 1.05 | "≈ 1.05" |
| PARADIGM §0 | the two structural defects were attributed to `ASSESSMENT.md`, which does not contain them | stated as defects of the classical labelled analysis |
| Remark H.7a experiment | not disclosed: the six kernels converge to a common configuration (spread 1.05 → at most 5.6e-17), so the effect is measured at a collapsed equilibrium | disclosed in PARADIGM, RESULTS_PARADIGM, the README and the manuscript; `exp_orbit_orientation.py` now reports the spread (table regenerated from `raw_sigma2_slow16.npz`; every other number unchanged) |
| README, ASSESSMENT §4, Simulation_Report | "0.1–0.6 ms per step against 13–1900 ms per update": the proposed method takes 4.11 ms at n = 3; NTK updates every step at n = 2, 3, 6 (not every 50–100) and is cheaper than the proposed method at n = 3 | exact per-step and per-update times; the cost edge is over every GP baseline at every n and over NTK except at n = 3 |
| README, ASSESSMENT §4, Simulation_Report, `plot_comparison_hd_complete.py` | not disclosed: every kernel-regression baseline (including the fair n = 1 re-run) regresses on noise-free samples of the true f(x, t), i.e. is given the state derivative; the adaptive law uses only the tracking error | disclosed; the table caption regenerated (numbers unchanged) |
| ASSESSMENT §4 | "NTK is first … by 10–25×": the factors are 25× (n = 2) and 3.4× (n = 6); "3–4× more accurate in every seed": per-seed ratios 2.79, 3.69, 3.48; "20× lower cost": 22.6× | corrected |
| ASSESSMENT, THEORY_REVISED, Simulation_Report | "three to four orders of magnitude (900×–11 000×)": the oracle table gives 915×–10954× | the measured factors |
| THEORY_REVISED §6 (errata) | several quotations of version 1 were not verbatim; E10's "100% of trials …" is from the old README, not version 1 | verbatim quotations with correct attribution |
| THEORY_REVISED §4 | "superlevel c > c_inf" (sublevel set meant); level-set argument credited to Lewis et al. 1996 (1999 meant) | corrected |
| README | the withdrawn statements were not quoted verbatim | verbatim quotations of the old README |

**Fixes to the manuscript text only.** The abstract's "not a gradient flow" sentence had dropped "of any free energy in
the W_{G⁻¹} metric" and "in the checked instances" (critical). The status of Corollary E.2 had dropped "as stated in its
abstract (its hypotheses were not checked against the full text)" (critical). Other fixes: missing verification
pointers and novelty verdicts in status lines; an over-general sentence in the verification appendix; caption wording;
cross-references to the appendix; layout.

**Noted, not changed.** Overloaded symbols (u, ε, σ, X, R, B, J, e) and a few symbols used before or without
definition (f̂*, W_∞, X_t, θ_max, P₀, S₀, b_odd, F*, R^m) are listed in the agents' reports for a notation pass.
"Within 3% in size" at k_e = 40 is a ratio of 0.969.

---

## Round 5 (research continuation and adversarial review, 2026-10-18)

**Design.** Version 4 adds Lemma E.0 with the explicit form of Theorem E, Theorem I, Theorem J, Theorem D′, Theorem C′
and Proposition K (PARADIGM.md §§4b, 5c, 6, 6c, 6d, 6e), with their code, tests and experiments. Five independent
reviewers were run, each told to refute: one on Lemma E.0 and Theorem E (proof gaps, counterexamples, hypotheses and
constants), one on Theorem I, one on Theorem J, one on Theorems D′ and C′, and a transcription checker for every number,
test name and cross-reference of the first batch. Their reports are kept as `review/REVIEW_round5_E0.md`,
`review/REVIEW_round5_I.md`, `review/REVIEW_round5_J.md`, `review/REVIEW_round5_DC.md` and
`review/REVIEW_round5_consistency.md`. Proposition K was reviewed by its own numerics (K-1 below). The confirmed issues
and their fixes follow, grouped by result.

**Verdicts.** Every result of version 4 holds, several after fixes the referees found. Theorem I(iv) was false as first
stated (the curl term survives at perfect fit on the orbit) and was corrected and strengthened. The proof of Theorem J had
a regularity gap in its implicit-function and contraction steps, now closed with the fast variable ζ. Theorem D′(c)'s
per-particle bound did not follow from exchangeability and has a pathwise proof now, and its numerical check was
unsupported until the experiment was redone with 24 seeds. The ω-limit argument of Theorem C′(i) was invalid and is
replaced by an LP-stability argument for generalised solutions, and its novelty label was an overclaim (the law is the
Lloyd descent of coverage control). Proposition K(iii)'s first draft used the wrong linearisation, which its own numerics
exposed. Lemma E.0 and Theorem E hold as stated, with two major fixes to the quoted numbers (a catastrophic cancellation
in the code and a misattributed accounting of the conservativeness factor). Thirteen items are marked major (statements,
proofs, quoted numbers, claims); some fifty are minor or cosmetic.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| I-1 | major | Theorem I(iv) claimed A = 2K at perfect fit (r_μ = 0 on the orbit). False for n ≥ 2: the curl term −∫Φ̃ᵀ(Dr − Drᵀ)Φ̃' dρ survives, because r_μ = 0 on a curve constrains only its tangential derivative. Counterexample f := f̂_μ + he₁ with h = (x₁/1.2)² + (x₂/0.8)² − 1 on the test ellipse: sup_𝒪|r_μ| = 4e-16, yet B − 2κI = −∫φφ'∂₂h dρ·(e₁e₂ᵀ − e₂e₁ᵀ) with entries 0.02–0.12 against |2κ| = 0.05–0.13; the exact second variation differs from 2K by 38% | (iv) restated: at perfect fit B = −∫(Dr−Drᵀ)φφ' dρ + 2κI, and A ≡ 0 iff T_𝒪 = 0 and Dr_μ = Dr_μᵀ ρ-a.e.; A = 2K when r_μ vanishes near the orbit. Propagated to §0, the H.7(e) bullets, Remark I.1, §8, README, RESULTS and the literature audit. New test `test_perfect_fit_on_the_orbit_alone_keeps_the_curl_term` (the closed form of the surviving term holds to 1e-12) |
| I-2 | minor | "reciprocal at μ" was undefined, and (i) was read as a pointwise statement, while H.7(e)'s hypothesis concerns a neighbourhood of μ (the triangle of measures) | defined: reciprocal at order k_e⁻² := A_μ ≡ 0; (i) states the consequence as "only if A_ν ≡ 0 for all ν near μ"; new part (vi): since A_ν is affine in ν, no open set of measures is reciprocal, for any n, so the averaged loop is nowhere a gradient flow at order k_e⁻² (proof along the segments toward a Dirac mass) |
| I-3 | minor | "the performative defect is proportional to the residual" is false for n ≥ 2 (the curl term) and an upper bound only for n = 1 | the bound |A_even| ≤ |a||a'|[∫‖Dr−Drᵀ‖dρ + (2e^{−1/2}/s_min)∫|r|dρ] in (ii); "bounded by" throughout |
| I-4 | minor | Remark I.1 said the non-gradient correction is O(k_e⁻²) relative to the leading drift; it is O(k_e⁻¹) relative (both are proportional to the residual at μ*_τ), O(k_e⁻²) absolute | corrected |
| I-5 | minor | (v) omitted "connected" | added |
| I-6 | minor | §0 said Theorem I "settles when" the loop is non-reciprocal: an overclaim for general n ≥ 2 | rephrased |
| I-7–I-10 | cosmetic | "proportional" → "bounded by"; H.7(e) step 5 should pass from the centroid measure to μ (ε → 0, A affine); "F₂((0.5,0))" ambiguous; ϑ_z's admissibility unstated | all fixed |

**Independent confirmation (Theorem I).** The referee re-derived every formula and checked aᵀB_μa' against the exact
second variation of the averaged potential (finite differences, Richardson in k_e = 80/160), including parameters at the
saturation limits (a = ±2.985, c = ±1.9995, s = s_min): agreement 6e-15 with `averaging.py` and 2e-4–9e-4 relative with
the exact second variation; the zero-feature pairwise criterion confirmed at the exact P by the triangle-of-measures
circulation with ϑ_z as a vertex (ratio 0.993 at k_e = 160); the quoted cyclic sums are those of `test_averaging.py`.
Figure-eight orbits have zero signed area but non-zero circulation kernel (max |K| = 0.049), consistent with the
"winding" wording; the segment and a forwards-then-backwards circle give |K| ≤ 1e-16.

**Consistency review (transcription and cross-references).** The reviewer compared every number, test name and
cross-reference of the version-4 prose with the code, the tests and the generated tables, and ran the eleven new tests.
The formula-to-code side was found clean: every constant of Lemma E.0 and Theorem E and every term of Theorem I(ii)–(iv)
matches `excitation.py` and `nonreciprocity.py`, all cited test names exist and assert the quoted tolerances.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| M1 | major | PARADIGM.md and the README cited a "round 5" of the review log that did not exist yet | this section |
| M2 | major | RESULTS said the nonlinear envelope test covers "one e-fold" of the explicit envelope; with the test's constants (explicit nonlinear rate 3.7·10⁻⁵) the 8000 s horizon covers about 0.3 e-folds | RESULTS and the test's docstring and comment corrected to "about 0.3 e-folds of the explicit envelope"; the inequality is still checked at every sample |
| M3 | major | README open problem 3 still listed "a general non-reciprocity theorem beyond the checked instances" although Theorem I settles n = 1 and the perfect-fit case | reworded to what remains open (n ≥ 2 with a non-zero residual in closed form; the orientation experiment with dither) |
| M4 | major (numbers) | the PE level α of the slow-adaptation table was under-resolved (200 Gauss–Legendre nodes over the 20π window; the value varied with t although a full-period Gramian is t-independent): α = 1.07·10⁻⁵ against the converged 1.36·10⁻⁵ at g = 0.01, so the quoted ratios exact/explicit were inflated by about 1.6. The bound itself stays rigorous (α under-estimated) | `pe_level` evaluated with 800 nodes in the experiment and in the test; the table regenerated and every quoted number re-quoted (Remark E.0a, RESULTS) |
| m1 | minor | "the optimised ε accounts for the remaining factor of about 2" between the formula value and the measured ratio mis-attributed the slack: optimising ε gives ×1.4; the Young step and p_hi ≤ ¾ give the rest | the closed-form-ε column added to the table; sentence rewritten |
| m2 | minor | the "−1.92 over the last octave" slope was not printed in the table | printed by the script |
| m3 | minor | "the ratio saturates" was not shown by the table (still decreasing at g = 0.01) and the test only checks no growth | "decreases towards a constant of order 3κ_Ψ(1 + Ψ_D/(kΨ_M))²"; the test described as a no-growth check |
| m4–m7 | minor | README H.7 row, PARADIGM §0 and the H.6 route still said "non-reciprocal in the checked instances"; RESULTS and the README layout row did not list the new modules and tests | updated |
| m8 | minor | the manuscript lagged the theory | rebuilt as version 3 (this round) |
| m9 | minor | the module docstring of `excitation.py` stated a different closed-form ε and γ than the code and Lemma E.0(ii) | aligned |
| m10 | minor | the new sources of the E.0 audit section had no [S]-numbers | [S15]–[S17], marked as quoted from memory |
| c1–c7 | cosmetic | an ambiguous "tables of … Corollary E.2"; three meanings of κ in §6 and a clash of (α, γ, κ) with the PE data in the proof of Theorem I(i); a garbled version line; Q_M, Ψ_M where the G-data were meant in Theorem E; 6.8·10⁵ quoted at 2 s.f.; the slow-adaptation ideal described without its reference; "fall like ω⁻²" for a fitted slope of −1.56 | the condition number renamed κ_Ψ; ϑ_z written as (a = 0, ξ); Q_{G,M}, Ψ_{G,M}; 3 s.f.; "driven by a two-tone reference of period 20π"; "approaches that law" |

**Lemma E.0 and Theorem E (refuter).** Verdict: both hold. The referee re-derived every step of the proof and verified
the pointwise form λ_max(N(t)) ≤ −γ of V̇ ≤ −γ‖w‖² over all w, with an independent quadrature and Hessian, for ten
regressor families (q > n, q < n, n = q = 1; k from 10⁻³ to 10⁴; Ψ_D/Ψ_M up to 100; non-periodic quasi-periodic PE;
near-degenerate PE with α ≈ 5·10⁻⁶; slow adaptation over 20π) and eight choices of ε, without a violation (worst margin
−7.5·10⁻³ relative, typically ≤ −0.3); the Jacobian, the similarity to w-coordinates, the PE-data transfer, the Floquet
code (against a DOP853 monodromy and a 240-period fit) and the ω⁻² asymptotics were confirmed numerically.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| E0-1 | major (code → numbers) | γ = λ_min of the 2 × 2 form was computed as ½(a + b) − √(¼(a − b)² + ¼c²), which cancels catastrophically when b, c ≪ a (every slow-adaptation row: a ≈ k_e, γ ~ 10⁻¹⁴–10⁻¹⁶ against a floor a·ε_mach ≈ 10⁻¹⁵); the g = 30 entry of the table was exactly 2⁻⁵⁰, +133% off the true value, and for ill-conditioned PE γ underflowed to 0 | γ := (ab − c²/4)/λ_max (relative error 10⁻¹⁶); the table regenerated and the quoted numbers re-quoted |
| E0-2 | major (text) | the slow-adaptation accounting of Remark E.0a: the formula 3κ(1 + Ψ_D/(kΨ_M))² is the leading-order value of the stated bound (2/3)γ_cf, while the table tabulates the optimised rate (asymptotically kα²/c₁², a factor 3 larger); "the exact rate is λ_min⟨ΨᵀΨ⟩/k" is off by a g-independent factor 0.66 on this instance (k_e = 6 is not large against the reference's bandwidth; a filtered averaging reproduces the Floquet rate to 5%); the PE level was under-resolved (M4); the T²Ψ_M³ term is still 13% of c₁ at g = 0.01; the regime condition is TΨ_M² ≪ k, not Ψ_M² ≪ k. The product of these factors is the "factor of about 2" that the text attributed to optimising ε | the paragraph rewritten with both formula values (1.47·10⁶ for (ii), 4.9·10⁵ for the optimised ε), the 0.66 bandwidth factor read off the table, the (c₁/c₁^∞)² = 1.31 correction, the limit ≈ 3.2·10⁵, and the regime stated as TΨ_M² ≪ k; the closed-form column added to the table |
| E0-3 | minor (code) | `optimise_eps` returned None for near-degenerate PE (the whole range of ε with γ > 0 below the grid floor), and every caller indexes the result | the closed-form ε is always a candidate, the grid starts at or below 10⁻³ of it, and the closed-form constants are returned when nothing beats them |
| E0-4 | minor (numerics) | the module's default 200-node Gauss–Legendre quadrature is under-resolved when a window holds many feature widths, so `pe_level`, `Q_matrix` and `StrictLyapunov` carried O(10%) errors in such cases (21% in table 2's α) | `pe_level` repeats the quadrature with twice the nodes at the minimising t and raises if the values differ (default check); the experiment and the tests use 800 nodes over the 20π window |
| E0-5 | minor (overclaim) | "a rigorous lower bound of the exact rate in every row", "the maximum is a rigorous bound": the inputs Ψ_M, Ψ_D are grid maxima, α a quadrature, ε a finite grid, and γ was subject to E0-1 | "numerical", with the margins (≥ 17×) stated; the Checked paragraph says what is grid data |
| E0-6 | minor (text) | "Ψ = sin t² … is not uniformly exponentially stable" rested on three simulated ratios | proof step (f) added to Lemma E.0: integration by parts with the Fresnel tail F(s) = ∫_s^∞ sin r² dr = O(1/s) shows that the decay over a fixed horizon tends to 0 as t₀ → ∞; the three ratios match 1 − (1/8)(1/t₀ − 1/(t₀ + 200)) to three digits |
| E0-7 | cosmetic | g_V := 1 + 2ε(Q_M + TΨ_M) was loose (measured ‖H‖ = 1.0017 against g_V = 1.505 on the test instance), and the statement of Theorem E wrote Q_M, Ψ_M for the G-data | g_V := 1 + ε(2Q_{G,M} + TΨ_{G,M}) in the text and the code |
| E0-8 | cosmetic | step 4 allowed ‖z‖ = r on the boundary of the open ball D | D taken closed |
| E0-9–E0-11 | cosmetic | the module docstring's different closed form (fixed before the report); "sharp up to a factor" in a test docstring; (iii) needs Ψ_D > 0 and its o(1) is O(1/ω) | fixed |

**Proposition K (collapsed configurations and the orientation; found by the numerics).**

| id | sev. | issue | fix (where) |
|---|---|---|---|
| K-1 | major (statement) | the first draft of K(iii) defined the orientation-odd shift through the Jacobian of the second-order drift at the collapsed critical point Θ*. The experiment showed the predicted k_e² × difference converging to 0.078 while the measured one converged to 0.060 (slope −2.18, both O(k_e⁻²)): the prediction missed the term ∂J₁(Θ*)[Θ̄₂ − Θ*], the third derivative of J_m contracted with the O(1/k_e) offset of Theorem J(d), which is O(k_e⁻²) like ∂b | K(iii) restated with the Jacobian of the normal form (N) at its own equilibrium Θ̄₂ and the decomposition E ± O with the third-derivative terms; the proof of (iii) rewritten; the experiment computes Θ̄₂ by Newton and tabulates both predictions (`exp_orientation_stability.py`): the corrected prediction is within 3% of the measured difference at k_e = 80 and its k_e²-scaled value converges to the measured limit |

**Theorem J (refuter).** Verdict: holds as stated, all four parts, under the stated hypotheses; no counterexample. The
referee verified on instance B (k_e = 20, 40, 80) the slaving constant (max|e − e_{Θ_p(t)}(t)| = 0.29–0.30·k_e⁻³), the
oscillation amplitude (max|Θ_p − Θ̄_p − εu(Θ*,t)| = 1.1·k_e⁻²), the convergence of k_e²‖Θ̄_p − Θ* − prediction‖ (A: 0.60 → 0.57;
B: 4.6 → 5.1), the O(k_e⁻²) convergence of every parameter Floquet exponent, the real slow invariant subspace of the
monodromy as a Riccati graph δe = S(0)δΘ with S = ε∂_Θr + O(ε²) (relative error 0.067 → 0.016), b_osc = 0 to 7·10⁻¹⁵ on
instance B, and the untested unstable case: an index-1 collapsed saddle of J_2 has a collapsed periodic orbit with
exactly one multiplier above 1 (1.00168 against the predicted 1.00185 at k_e = 20; 1.00089 against 1.00092 at 40). The
report is `review/REVIEW_round5_J.md`.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| J-1 | major (proof) | Steps 4 and 6 used additive O(ε³) remainder bounds where Lipschitz/C¹ control is needed: "Ã_ε = ⟨H⟩ + O(ε)" does not give DÃ_ε = D⟨H⟩ + O(ε) for the implicit function theorem on an ε-independent ball, and an additive O(ε³) in Step 6 yields convergence into an O(ε²) ball, not a contraction; as written, (a)'s uniqueness and (c) were not established | Step 2 rewritten with the fast variable ζ := e − e_Θ(t) and its invariant set Z_ε = {‖ζ‖ ≤ Cε³}; a regularity paragraph: e_Θ is C² in Θ (f enters r_Θ only through the smooth f̂_Θ) with ∂_Θe_Θ = ε∂_Θr_Θ(x_m) + O(ε²), the remainder is O(1)-Lipschitz in ζ and O(ε³)-Lipschitz in Θ, the ζ-forcing O(ε²)-Lipschitz in Θ; Step 4 solves for ζ*(Y) first and applies the Lipschitz form of the IFT to Ã_ε with Lip(Ã_ε − ⟨H⟩) = O(ε); Step 6 contracts with the 2 × 2 Lipschitz matrix [[1 − εT₀λ_*(1 − κ), Cε], [Cε³, e^{−k_eT₀/2}]]; f ∈ C² suffices |
| J-2 | minor | Step 1: "Θ moves by O(ε² log(1/ε))" during the transient is false for ‖e(0)‖ = O(1); it is O(ε)‖e(0)‖ (measured 1.7–1.9·‖e(0)‖/k_e) | O(ε‖e(0)‖ + ε² log(1/ε)) in Step 1; Step 6 absorbs it by halving ρ₀ |
| J-3 | minor | (S)/(N) were claimed uniform on W_ε = {‖w‖ ≤ 3C₁ε²}, but the O(ε³) slaving needs the thinner set ‖ζ‖ ≤ Cε³ | Z_ε replaces W_ε in Steps 2–4 and 6 |
| J-4 | minor | the eigenvalues of I + εA + O(ε²) are 1 − εT₀ν_j + O(ε²) only because A is diagonalisable (needed for the repeated eigenvalues of a collapsed configuration) | Bauer–Fike sentence in Step 5 |
| J-5 | minor (numbers) | "residual ≤ 3·10⁻¹³ in every case" is false: 8.6·10⁻¹² at k_e = 160 (A) | "below 10⁻¹¹" |
| J-6 | minor (text) | Remark J.1 bullet 3 named the wrong mechanism: at a collapsed point the separation Hessian is σ₀G⁻¹ − R₂; for n = 1 and σ₀ = 0 the amplitude row of R₂ vanishes exactly at a leak-free critical point, so only the m − 1 amplitude-redistribution directions are flat (verified: eig ∇²J_2 = −0.0308, 2.6·10⁻¹³, 0.038, …), while for n ≥ 2 only the amplitude block vanishes and the collapsed point is generically non-degenerate | bullet rewritten with the mechanism (and a pointer to Proposition K) |
| J-7 | minor (verification) | the error-channel part of (b) was unverified (the FD monodromy's fast multipliers are noise at 10⁻⁸–10⁻⁹); (c) untested | Liouville's trace identity Σlog|λ| = ∫tr ∂F/∂z dt recovers the product of the fast multipliers from the resolved slow ones: `exp_periodic_orbit.py` now reports the mean error-channel exponent against k_e (the referee measured −(k_e − L_eff)T₀ with L_eff ≈ −0.07…−0.09); a basin run from 25× the offset converged with per-period factors inside the predicted band (not added as a test: runtime) |
| J-8 | cosmetic | "∂_ΘF = O(ε) because F = O(‖e‖ + ε)" (right conclusion, wrong reason); ρ₀ in ‖·‖ in (a) but ‖·‖_{G⁻¹} in (c); {‖w‖ ≤ 1} not intersected with the tube; "diffeomorphism of 𝒦"; "sup‖ẋ_m‖/s ≈ 3–10" (2.6 at the minimiser's width); Remark J.1 bullet 1 stated an unproved converse; "makes Remark H.7a a theorem" holds for the period-mean prediction only, and when the leading offset is small the O(k_e⁻²) remainder can dominate its direction (the saddle instance: cosine 0.006 at k_e = 20, 0.84 at 80) | all reworded; the saddle caution added to Remark J.1 |

**Theorems D′ and C′ (refuter).** Verdict: both hold, with fixes. For D′ the referee confirmed the coupling argument
(exact noise cancellation, the assembly of L_b, the i.i.d. property of the ϑ̄_j, the Grönwall step, the Fournier–Guillin
exponent) and diagnosed the failed numerics: the closed loop is contracting (sup_t D_m/∫ζ ≈ 1/k_e), D_m plateaus by
t ≈ 5, and the driving fluctuation is a persistent sample fluctuation of the first m particles that does not
self-average in time, so sup_t D_m has O(1) scatter between seeds; with 12 seeds (T = 5, M = 2048) the referee measured
the slope −0.54 ± 0.07 and mean/(1/m − 1/M)^{1/2} flat at 0.28–0.33. For C′ the referee confirmed (i)–(v) as inequalities,
the code (LP constraints, the dropped redundant constraint, the variance formula, the Pythagorean split) and the
numerics up to one rounding, and found two proof-level problems. The report is `review/REVIEW_round5_DC.md`.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D1 | major (proof) | (c): E sup_t‖ϑ_i − ϑ̄_i‖ ≤ TK_Φe^{ΛT}m^{−1/2} "by exchangeability" does not follow: exchangeability gives E‖Δ_i(t)‖ = EΔ̄(t) at fixed t only, and sup_t Δ̄ ≤ (1/m)Σ sup_t‖Δ_i‖ is the wrong direction | the pathwise per-particle bound ‖Δ_i(t)‖ ≤ 2(e^{L_bt} − 1)e^{Λt}∫₀^tζ (Grönwall for Δ_i driven by ‖δ‖ ≤ D), hence E sup_t‖Δ_i‖ ≤ 2(e^{L_bT} − 1)TK_Φe^{ΛT}m^{−1/2}; the fixed-t statement kept as the exchangeability consequence |
| D2 | major (numbers) | "proved and checked" was unsupported: the three-seed table had slopes −0.25 and −0.21 against −0.5, with per-seed scatter up to 9× at fixed m, and the numerics paragraph was a placeholder | diagnosis (above and in the diagnostic run: no particle deviating by more than 0.17, no basin switching, D_m growing smoothly or plateauing); `exp_chaos.py` rewritten with 24 seeds, the shared-particle law (1/m − 1/M)^{1/2}, mean ± s.e., median, bootstrap s.e. of the slope and the raw per-seed values saved; result −0.47 ± 0.06 (τ = 0), −0.46 ± 0.05 (τ = 2·10⁻³), 0.94 ± 0.08 and 0.92 ± 0.08 against (1/m − 1/M)^{1/2}; the numerics paragraph written from it, the three-seed table withdrawn |
| D3 | minor (text) | "constants independent of the dimension": K_Φ = a_max√n and the Lipschitz bounds grow with n; only the exponent is dimension-free | reworded in (b), (c), Remark D′.1, the §0 row, §8 and the README |
| D4 | minor (proof) | (a): the Picard map need not be a self-map of C([0,T₁]; X), and the standing bound presupposes existence | compact neighbourhood X′, bounded input, local solution and continuation through the tube bound; the uniqueness class stated |
| D5 | minor | Theorem D(b)'s tube bound was invoked for a law it does not cover | a sentence: its proof uses only ‖f̂‖ ≤ K_Φ and the pathwise error ODE |
| D6 | minor (experiment text) | "the M^{−1/2} floor": there is none (D_M ≡ 0); the forcing variance is ∝ 1/m − 1/M; the bound quoted for the M-system's first m particles had the gap of D1 | docstring and table text rewritten (the shared-particle law; the per-particle bound) |
| D7 | cosmetic | W₁ Euclidean against the G⁻¹ metric elsewhere; μ̄ used for the path law and the marginal; "(d) sharpens D(d)" compares different laws; "the dimension-free rate for the state is new" is the standard output of Sznitman's coupling; the P₂ moment condition for Fournier–Guillin was not stated | all reworded; the novelty column says "standard output … new in that the plant is in the loop" |
| D8 | cosmetic | the product norm for L_b and the global-in-x derivative bounds were not stated | stated in the standing bound |
| C1 | major (proof) | (i): "∫V_bar² < ∞ and V_bar continuous except at finitely many switches ⇒ ω-limit points have V_bar = 0" is invalid: switches may accumulate, V_bar jumps upward at switches (seed 0: 1.074 → 1.100, 0.428 → 0.451, 0.247 → 0.298), and Lemma C.3's no-switching does not extend to m₀ ≠ m (the pure barycentric leak genuinely switches), so the existence of selection solutions was assumed, not proved "as in §4" | two solution concepts (vertex-selection, existence assumed; generalised, which exist as solutions of an upper semicontinuous convex-valued differential inclusion); the ω-limit statement proved by LP stability (finitely many vertices; the optimal face is upper semicontinuous; V_bar^γ is 1-Lipschitz for a fixed plan; bounded speed ⇒ ∫V_bar² = ∞ unless some optimal plan at ϑ_∞ is centroidal); the two-way tie computed (the sum of the two switching rates is (2σ_m/m)Σ‖T̄^γ_i − T̄^ρ_i‖² ≥ 0, so two-way ties are not sliding modes); multi-way ties left open |
| C2 | major (novelty) | "new as a control construction": ṗ_i = −k(p_i − C_{V_i}) is the continuous-time Lloyd descent of coverage control (Cortés–Martínez–Karatas–Bullo 2004); OT-based coverage control is [89]; equal-mass Lloyd with OT plans is Balzer–Schlömer–Deussen 2009, Mérigot 2011, de Goes et al. 2012; [14] is an application of CVT, not CVT theory | [S18]–[S20] added to the audit with a C′ section; the §0 row, §4b, §8 and the README now say "the Lloyd descent of coverage control with equal-mass OT plans, used as a σ-modification leak; new here are the induced-coupling bound, the Lyapunov inequality with the bias b_bar and the transfer reading" |
| C3 | minor | no scope paragraph for (ii)–(iv), and "the bias is b_bar, as (ii)–(v) state" overreached: on the experiment's clouds Δ_gap is O(1/m₀) (one-atom exchange) and 4η(0)(b + η(0)) exceeds it by an order of magnitude, so (iv) is vacuous on every tested instance | scope paragraph added; the experiment reports the exchange bound on Δ_gap at μ*, η(0) and 4η(0)(b + η(0)) with the failure factors; the numerics text attributes the observed bias to the fixed-point statements (i), (v) |
| C4 | minor | the law was not defined with a vertex plan, "optimal, hence a vertex" is a non sequitur, and Δ_gap only makes sense for atomic ν | "vertex optimal plan" in the law; (ii)-equality and (iv) restricted to atomic ν; the generalised-solution case argued vertex by vertex |
| C5 | minor | the baseline (C) is a random 4-subset: the fair right-size baseline (a leak toward the cloud's equal-mass 4-quantiser) has the same bias b_bar | the "anchors better" claim qualified in §4b and the table text |
| C6 | cosmetic | "V_bar² ≤ 2·10⁻²⁵" against 2.2·10⁻²⁵; "= 0" is 0.000 at printed precision; the leak's seed 0 (1.1642) beats the best of 20 Lloyd runs (1.1708), so "compare the last two columns" did not show what it claimed | corrected |
| C7 | cosmetic | (v): Var_{γ*} = 0 also when m₀ divides m or atoms coincide; b_bar, T̄* and "centroidal" are plan-dependent at ties; degenerate coincident-particle critical points are fixed points | stated |
| C8 | cosmetic (code) | `optimal_plan` used HiGHS without a vertex guarantee, entries could be −10⁻⁹ | `method='highs-ds'` (dual simplex) and clipping; tests pass unchanged |
| C9 | cosmetic | "AGS Thm 10.4.12" cited for the superdifferential of ½W²(·, ν) with an unchecked sign convention | replaced by a pointer to [AGS, §10.4] with the generalised velocities described directly |

## Round 6 (research continuation and adversarial review, 2026-10-18)

**Design.** The second research increment adds Lemma E.0′ (the sharp explicit rate of the slow-adaptation regime),
Corollary E.0′.2 (Theorem E with its constants) and Theorem D″ (uniform-in-time propagation of chaos under a dominant
leak) to PARADIGM.md §§5c and 6, with code, tests and experiments. Three referees were run: on Lemma E.0′ (proof,
constants, code and certificate; `review/REVIEW_round6_E0prime.md`), on Theorem D″ (proof, constants, experiment and
reading; `review/REVIEW_round6_Dpp.md`) and on Theorem J′ with Proposition J′.2 and `dither.py` (proof, code, numerics;
`review/REVIEW_round6_Jprime.md`). The confirmed issues and their fixes follow.

**Verdicts.** Lemma E.0′ holds: the referee re-derived the fast-variable equation, the order-ε and order-ε² cancellations
(including N₂ = 0), the list of remainder terms and every per-direction bound, and certified the inequality V̇ ≤ −γV
independently (its own quadratic forms, Gauss–Legendre antiderivatives instead of the FFT) on five regressor families —
the sinusoid, two random multitones (n = 2, q = 4), a hand-made bursty regressor with κ(M̄) = 3.3·10³ and the Theorem E
ideal — over 416 admissible (k, δ) pairs, on full grids, at off-grid times and states and along trajectories, without a
failure; no counterexample. One coefficient had the wrong dimension, two constants carried a factor-2 slack, and the
stated asymptotics missed a first-order term (all fixed). Corollary E.0′.2's proof is right; its first test measured the
integrator's discretisation floor and was repaired (E0P-4 = E0P-f). Theorem D″ holds as a theorem; its reading did not: the
run presented as the learning regime was a slow collapse (DPP-3), the constants carried a factor 2 (DPP-1), and the
synchronisation claim holds only without dither (DPP-2); the instance, the constants and every text were redone, and the
new learning-regime run shows the coupled distance growing in time. Theorem J′ holds with fixes (a derivative count in the
regularity bookkeeping; the Monte Carlo error bars and burn-in of the first run), Proposition J′.2 holds, and the first-order
numerics were rerun with the corrected design.

**Lemma E.0′ (refuter).** Selected confirmations: the stated ẇ equals the chain-rule derivative to 5·10⁻¹⁴ on four
regressors with εΨ_M up to 3 (all orders through ε⁶ exercised) and the fourth-order finite difference of w along DOP853
trajectories to 10⁻¹⁰; Ṗ₁ = 2(M − M̄), Ṡ = F, ⟨F⟩ = 0 to 4·10⁻¹⁵; the ε² coefficient vanishes to 10⁻¹³; the remainder
R_ζ is exactly the listed one (10⁻¹²); the (T/4) antiderivative bound is sharp (ratio 0.9991 on tanh(500 sin t)); on the
ideal at k = 400, sup V̇/V = −1.0755·10⁻⁹ against −γ = −8.594·10⁻¹⁰ (δ = 0.2) and −1.0775·10⁻⁹ against −9.673·10⁻¹⁰
(δ = 0.1); exact/(λ_min/k) = 0.9885 … 0.99995 at k = 50 … 800; E.0′/E.0 = 1.07·10⁵; the sinusoid's exact/guaranteed →
1.2512 at k = 3200 (the predicted 1/(1 − δ) = 1.25). The report is `review/REVIEW_round6_E0prime.md`.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| E0P-1 | minor (statement, code) | the ε⁶ coefficient of c₂ was Ψ_M⁵√ρ; the bound of ‖ΨM³ζ‖ is Ψ_M·Ψ_M⁴‖Mζ‖ ≤ Ψ_M⁶√ρ‖ζ‖_{M̄} (wrong dimension; violated by the factor Ψ_M on 3(sin t, cos t); numerically irrelevant) | Ψ_M⁶√ρ in the statement, the proof and `constants()` |
| E0P-2 | minor (text) | "γ/2 = (1 − δ − O(k⁻²))λ_min/k" (and the corollary's λ′) missed the first-order deficit of the denominator 1 + εp₁ + ε²p₂: the ratio is (1 − δ)(1 − εp₁ + O(k⁻²)) (sinusoid: 0.638, 0.900, 0.982, 0.9968 at k = 20, 50, 200, 1000 against 1 − εp₁ = 0.843, 0.937, 0.984, 0.9969) | asymptotics rewritten in the lemma and the corollary |
| E0P-3 | minor (text) | "k₀ ≈ (4ρρ_DD/δ²)^{1/4} from (A)" holds only when (A) binds (ideal: 226 against 225; bursty: 63 against 61); on the sinusoid (C) binds (9.85 against 4.5) and on a multitone with TΨ_M² = 134 the threshold is 269 against 17.7 and not monotone in δ | qualified ("when (A) binds …; otherwise (C) sets k₀") |
| E0P-4 | major (test) | the corollary's first envelope test failed (worst 2.81 against ≤ 1): the fixed-step RK4 loop drifts to ‖w‖ ≈ 5.8·10⁻⁶ from an exact start at the ideal (3.0·10⁻⁷ at dt = 0.01, 1.7·10⁻⁸ at dt = 0.005), four times the explicit ball r′₀ = 1.42·10⁻⁶; the linearised and the step-2 fields obey the envelope (ratio 0.084) | = E0P-f: dt = 0.005 in the test (passes, 1.5 e-folds) |
| E0P-5 | cosmetic (constants) | π₂ and p₂ carried a factor 2 of slack inherited from the draft's F − 2N₂: the proof's chain gives (T²/8)Ψ_M³(2√ρ + 1) and Ψ_M² + (T²/4)Ψ_M⁴ | tightened in the statement, the proof and the code; tables 4 and 5 regenerated |
| E0P-6 | cosmetic | "single-tone regressor (ρ = ρ_D = ρ_DD = 2)" holds at unit frequency only (ρ_D = 2ω², ρ_DD = 2ω⁴) | "Ψ = (sin t, cos t)" |
| E0P-7 | cosmetic | "within 12% of the exact rate above its threshold" is the δ = 0.1 figure; δ = 0.2 gives 26% | qualified |
| E0P-8 | cosmetic | "1 − exact/(λ_min/k) ∝ k⁻²" is regressor-dependent (≈ 29/k² on the ideal, slope −1.98; ∝ k⁻⁴ on the sinusoid); "slow manifold to second order" for an expansion through ε³ | "O(k⁻²), attained on the ideal"; "expanded through order ε³" |
| E0P-9 | cosmetic (table) | table 4's "exact/E.0 = 1.3·10⁵" rests on the 80-point grid maxima of Ψ_M, Ψ_D (Ψ_D = 0.94–0.96 for 60–200 points); the converged value is 1.19–1.36·10⁵ | the table's header already declares grid maxima; the text says "five orders of magnitude" |

**Found by the author during the derivation (logged for the record).**

| id | sev. | issue | fix (where) |
|---|---|---|---|
| E0P-a | major (first draft) | the first construction (one-step slaving w = e + εΨζ, V = ½ζᵀ(I + εP₁)ζ + β‖w‖²/2) needs k ≥ √(ρρ_D)/δ, which is ≈ 2400 on the Theorem E regressor (ρ_D = 3.7·10³): vacuous | second-order slaving e₂ = −εΨζ + ε²Ψ̇ζ − ε³ΨMζ and a second averaging correction P₂; the threshold becomes (4ρρ_DD/δ²)^{1/4} ≈ 320 |
| E0P-b | major (first draft) | the first version kept a second-order averaged matrix N₂ = ½⟨P₁M + MP₁⟩ in the rate and bounded it crudely; it vanishes identically (⟨P₁⟩ = 0 and P₁M + MP₁ = P₁M̄ + M̄P₁ + ½d/dt P₁²) | N₂ = 0 proved in the statement; the rate is λ_min(M̄)/k exactly at second order, which the measured k⁻² deviation of the exact rate confirms |
| E0P-c | minor (code) | the first certificate compared λ_max(H^{-1/2}N H^{-1/2}) with −γ, missing the factor 2 of V̇/V = 2xᵀNx/xᵀHx, and reported spurious violations | fixed before the code entered the repository |
| E0P-d | minor (constants) | the zero-mean antiderivative was bounded by (T/2)‖f‖_∞; the sharp bound is (T/4)‖f‖_∞ (circle-distance argument) | adopted in the proof and the code (p₁, π₁, π₂, p₂) |
| E0P-e | minor (text) | Remark E.0′.1 claimed that Theorem E's explicit ball radius improves by the same factor as the rate; the radius carries the anisotropy p′_lo/p′_hi ≈ 1/β of the quadratic form and shrinks on well-conditioned regressors (table 5: 1.5·10⁻⁶ against 6.1·10⁻⁶ at k_e = 40) | Corollary E.0′.2 stated with its own radius; the remark corrected |
| E0P-f | minor (test) | the fixed-step RK4 closed loop (`rk4_closed_loop`, dt = 0.02) has a discretisation floor ‖w‖ ≈ 6·10⁻⁶ from an exact start at the ideal (∝ dt⁴: 4·10⁻¹⁰ at dt = 0.002), above the corollary's explicit ball r′₀ ≈ 1.4·10⁻⁶; the first envelope test failed for this reason (the state hovered at the floor instead of decaying) | the corollary test integrates with dt = 0.005 (floor ≈ 2·10⁻⁸); the Lemma E.0 envelope test (r₀ ≈ 10⁻⁵ at dt = 0.02) sits above the floor and is unchanged, with the caveat recorded here |

**Theorem D″ (refuter).** Verdict: the theorem holds as stated and the proof is correct (bracket decomposition, the
contraction of the second bracket with the G-free leak and ‖ē‖ ≤ E_X, the Metzler comparison with Dini derivatives, the
Hurwitz criterion, the −A⁻¹ formula, Eζ_s ≤ K_Φm^{−1/2}, exchangeability at fixed t, the W₁ triangle inequality); the referee
also verified the comparison inequality pathwise on six m-against-M couplings with τ > 0 (ratio ≤ 0.71 on the state channel,
≤ 0.056 on the particle channel). The *reading* did not hold in one major and several minor places. The report is
`review/REVIEW_round6_Dpp.md`.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| DPP-1 | minor (sharpness) | the first bracket was bounded by 2L_b′‖δ‖, double-counting ‖e − ē‖ = ‖δ‖; the sharp bound is L_b′‖δ‖, so A = [[−κ_e, L_ϑ],[L_b′, −κ_ϑ]], det A = κ_eκ_ϑ − L_b′L_ϑ and the bounds carry L_b′, not 2L_b′ (confirmed pathwise with factor 1) | statement, proof, §0 row, experiment, README |
| DPP-2 | minor (overstatement) | "the certified regime is the collapsed one, in which the statement is true and trivial": the synchronisation D⁺‖ϑ_i − ϑ_j‖ ≤ −κ_ϑ‖ϑ_i − ϑ_j‖ holds only at τ = 0; at τ > 0 the particles form the dither cloud N(ϑ⁰, (τ/σ)G) (rms spread within 3% of (τ tr G/σ)^{1/2}) and the coupled distance is of order m^{−1/2}, not zero — the theorem has content there, though the network is switched off (sup‖f̂‖ a few percent of sup‖f‖) | Status, Remark D″.1, §0 row, §8, README rewritten: "leak-dominated regime, in which the network does not learn"; the certified run now also made with dither |
| DPP-3 | **major** (reading of numerics) | the σ = 0.02 run was labelled "spread / learning regime": Proposition K applies to the instance, whose unique collapsed critical point is the dead kernel a = 0 (f odd, ρ symmetric, ϑ⁰ = 0) with σ_K = σ₀_crit/k_e = 0.00117, seventeen times below 0.02; the spread decays exponentially at every σ ∈ [0.02, 0.4] with rate σ − 0.0012 (half-life 37 s at σ = 0.02), the network's output stays below 0.15 against sup‖f‖ = 1.25, and the reported spread 2.57 at T = 80 was the leak-decay of saturated initial outliers measured in the unconstrained coordinates: a transient of a slow collapse, not a learning regime | the instance changed to a two-bump plant that is not odd; σ_K computed in the script by `periodic.sigma_crit`; the learning-regime run at σ = σ_K/3 over 1000 s with the particle spread in physical coordinates and the fit of the plant reported; all texts rewritten from the new table |
| DPP-4 | minor | the high-gain bullet "σ > 2λ_GL_ϑ²/k_e + O(k_e⁻²)" dropped the term λ_GE_XL₂ = λ_G(F + K_Φ)L₂/k_e of the same order (and larger here) | bullet rewritten with both terms |
| DPP-5 | minor | "every fixed σ > 0 certifies for k_e large" needs e(0) = 0 or ‖e(0)‖ = O(1/k_e), since E_X ≥ ‖e(0)‖ | hypothesis stated |
| DPP-6 | minor | "the regime in which J_m is strongly convex" follows only for the a-priori tube with e(0) = 0 (k_eE_X = F + K_Φ ≥ sup‖r_Θ‖), not for the measured E_X | stated for the a-priori tube, with the Hessian bound |
| DPP-7 | minor (experiment) | "certified" rested on an a-posteriori sample E_X = 0.052 (every 20th step, 8 seeds); with the a-priori tube radius (F + K_Φ)/k_e = 0.2236 the old constants did not certify the run (det A < 0), the sharp ones do | the script computes the a-priori E_X and certifies from it (σ_cert = 4.0, run at 1.2σ_cert); the measured E_X reported alongside |
| DPP-8 | minor (experiment) | `lipschitz_constants` under-estimated the suprema (L_x = 7.98, L_ϑ = 7.70, L₂ = 31.6, L_{xϑ} = 41.1 against the exact 8.087, 8.087, 44.4, 44.4) | closed forms for n = 1 (a_max e^{−1/2}/s_min, a_max/s_min²), confirmed in the script by a 1.2·10⁵-sample search with Nelder–Mead polish to 10⁻⁸ |
| DPP-9 | cosmetic (positive) | the second-moment hypothesis for the W₁ statement is redundant: Itô gives sup_t E‖ϑ̄ − ϑ⁰‖² ≤ max{z(0), ((λ_GE_XL_ϑ)²/σ + 2τ tr G)/σ} | hypothesis removed, the bound proved in the statement |
| DPP-10 | cosmetic | "u̇" for a Dini derivative; the Hurwitz condition listed κ_ϑ > 0 redundantly; for n > 1 the norms L₂, L_{xϑ} must be contracted with unit vectors as in §1 (1) | D⁺u; "κ_e > 0 and det A > 0 (then κ_ϑ > 0)"; the norm convention stated |
| DPP-11 | cosmetic | the slope −0.41 had no standard error; the bound column compares sup_t E D with E max_t D (right direction, unsaid) | per-seed slopes with s.e.; the comparison stated in the table |

**Theorem J′, Proposition J′.2 and `dither.py` (refuter).** Verdict: the theorem holds with fixes (every sign of Steps 3, 5,
6, 7 re-derived by hand and confirmed on exactly solvable fast–slow toys, including the sign of b_osc, which the referee
established independently of the text's integration by parts: with the sign flipped the remainder grows like 1/ε), the
proposition holds ((iii) confirmed to six digits by finite differences on the single-kernel instance; (iv) sound), the
code holds with fixes (every formula of `averaged`, the solver, the current, the circulation prediction and the splitting
confirmed independently; the Monte-Carlo error bars and the burn-in had quantifiable biases). The report is
`review/REVIEW_round6_Jprime.md`.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| JP-1 | minor (regularity) | the expansion needs five Θ-derivatives of the Poisson solution, not four: ψ₂′ contains ∇³φ through τ₀tr(G∇²(u·∇φ)), so B₀ψ₂ carries ∇⁵φ; and interior Schauder theory gives ψ_φ ∈ C^{5,α} only for φ ∈ C^{3,α}, not for φ ∈ C² as stated | the class of observables is C^{3,α} with polynomially bounded derivatives up to order three (statement, Remark J′.1, Step 4); Step 3 counts the derivatives (C⁵), with the remark that (b) alone needs C³ |
| JP-2 | cosmetic | w_φ written for a diagonal G while Step 7 allows a symmetric G | G diagonal stated in the theorem; the general contraction given in Step 3 |
| JP-3 | cosmetic | Proposition J′.2(i) stated the density expansion pointwise where the proof gives it weakly; (ii) said "up to O(ε²)" where the Langevin identification is exact | (i) stated weakly, with the L¹ check on the grid; the EPR remainder labelled formal; (ii) stated exactly |
| JP-4 | minor (numerics) | the batch-mean standard errors of `simulate_stationary` were biased low by 1.4–1.9: the spectral gap of L̄₀ on the instance is 0.286 per slow unit (relaxation 3.5 units), longer than the 2-unit batches | standard errors from the spread of the per-replica time averages (independent replicas); the experiment and the test rerun |
| JP-5 | minor (numerics) | the averaging window was not a whole number of reference periods (63.66, 127.33, … periods), so the phase was not averaged uniformly; negligible for the means, but of the size of the predicted circulation for Θ∧V | burn-in and window are whole numbers of periods |
| JP-6 | minor (numerics) | burn-in of 4 slow units ≈ 1.1 relaxation times from a start at π: ≈2.8% of the O(ε) shift survived in the window, biasing the measured shift toward π by the size of the claimed agreement | replicas start from the first-order law π(1 + εg) and burn in for 10 units (transient O(ε²)e^{−2.9}) |
| JP-7 | cosmetic (table) | the "resolved" criterion for the second-order slope used the biased standard errors | criterion kept with the corrected errors; the table says where the residual is not resolved |
| JP-8 | cosmetic | (iv) quantified over σ₀ and G idly and did not recall Theorem I's setting | reworded |
| JP-9 | minor (reading) | the first run's text claimed the orientation-odd part "to within a few percent"; the data give 7–8% at k_e = 10 and agreement within a 10–30% resolution at k_e = 20–80; the k_e⁻² residual of (c) is at the noise floor for k_e ≥ 20 and the even part exceeds its first-order prediction by 1.5–1.8× at k_e = 10 (an O(ε²) coefficient of order −0.4, the ratio c₂/c₁ ≈ 1.4 of Theorem J's instance A); Table 4 (circulation) is dominated by the O(ε²) term and, in the first run, by the phase bias of JP-5 | the numerics paragraph written from the corrected run with these qualifications; the second-order residual reported as not resolved; Table 4's caption says what it does and does not test |

## Round 7 (research continuation and adversarial review, 2026-10-18)

**Design.** The third research increment adds Corollary E.0′.3 (the explicit stationary parameter spread of the dithered
linearised loop, by Itô's formula on Lemma E.0′'s Lyapunov function, against the exact periodic Lyapunov solution) and
Lemma E.0″ (the slow-regime strict Lyapunov function for non-periodic persistently exciting regressors, with the period mean
replaced by a forward exponential average over a window T) to PARADIGM.md §6, with `excitation.dither_covariance_bound`,
`excitation.exact_stationary_covariance`, `excitation.ExpAverageLyapunov`, `excitation.lyapunov_exponent_rate`, two tests,
table 6 of `exp_explicit_rate.py` and the new `exp_pe_general.py`. One referee was run on both results (proof, constants,
code, numerics and reading; `review/REVIEW_round7_E0p3_E0dp.md`).

**Verdicts.** Corollary E.0′.3 holds: the Itô correction ½tr(Σ̃ᵀ𝐇Σ̃) with 𝐇_ζζ = I + εP₁ + ε²P₂ + βLᵀL was re-derived from the
block form of V, the martingale step, Grönwall, the three bounds and the averaged covariance kX were confirmed, and the exact
covariance routine was checked against an independent integrator (10⁻¹¹) and against Euler–Maruyama Monte Carlo (0.1–0.3%);
the Itô correction was pinned by stationarity (mean E V̇_det = −mean ν to seven digits in 24 cases) and the bounds dominate
the exact values for isotropic, diagonal, rank-one and random Σ and for anisotropic M̄. Lemma E.0″ holds: the derivative
identities, the signs of the residual N̄₂ = ½(F̄ − Ṁ̄), the per-direction bounds in the time-local norm and both directions
of the PE equivalence were confirmed, the referee's own Lyapunov function reproduced the certificate at off-grid states, and
an attack set (barely exciting, bursty, near-resonant, constant-component, n = 2, q = 3 and chaotic regressors; ≈250
admissible rows) produced no violation — in every row the slack between the certified and the true rate is exactly the Young
step. No major issue; fifteen minor and cosmetic ones, all fixed.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| R7-1 | minor | "the exact E‖e‖² is ‖Σ‖_F²/k²" (PARADIGM, `exp_explicit_rate.py`) — the numbers are ‖Σ‖_F²/(2k) = tr(M̄X)/k (the error rides on its slow manifold; tr(M̄X) = ½‖Σ‖_F² by the trace of the Lyapunov equation), ratio 1.00000 in 14 cases | formula corrected in both; asserted to 10⁻⁴ in the test |
| R7-2 | minor | "averaging exact to O(k⁻²)" holds when ΣΣᵀ commutes with M̄ (isotropic Σ in both instances); for anisotropic Σ the deviation is ±0.25/k on two sinusoids, the bound still dominating | sentence qualified in PARADIGM and in the table's reading |
| R7-3 | cosmetic | "finite moments of every order" in the martingale step; second moments are what is used | reworded |
| R7-4 | minor (test) | the check bound_w ≥ ½E‖e‖² was vacuous (E‖w‖² = 1.3·10⁻⁷ ≪ E‖e‖² = 10⁻⁴ at k = 20) | `exact_stationary_covariance` returns E‖w‖² = E‖e − Lζ‖² from the covariance trajectory and the slaving map L; the test asserts bound_w ≥ max_t E‖w‖² and records the looseness (> 10³) |
| R7-5 | minor | the Status said "conditions hold from k = 20 on" (true for the quasi-periodic regressor only; from k = 100 on the Lorenz one) and "within a factor 4–5 on the chaotic one" (only at T = 20; 7.5–11 at T = 10, 14–138 at T = 5) | Status rewritten |
| R7-6 | minor | "explicit constants give rates 10–25% below the measured ones": the deficits are 22–79% at k = 50, 5–19% at k = 100, ≤ 5% at k = 200, ≈ 1% at k = 400 | numbers corrected |
| R7-7 | minor | Lorenz ‖N̄₂‖ quoted as 1.05, 0.30, 0.17; the table has 1.054, 0.341, 0.422 (not monotone in T), and the generated reading said p₁, p₂, ‖N̄₂‖ "grow with T" (‖N̄₂‖ decreases with T on the quasi-periodic regressor) | values corrected; the reading in `exp_pe_general.py` and its table rewritten |
| R7-8 | minor | T = 2π ratios quoted as 3.7 → 1.9 and 3.1 → 1.9; the table has 3.69 → 2.16 (δ = 0.2) and 3.10 → 1.92 (δ = 0.1) | corrected |
| R7-9 | cosmetic | "0.5/k to 0.5%": 0.6% at k = 20, 0.8% at k = 400 | "within 1%" with both values |
| R7-10 | cosmetic | "the long-time means minimise the norms below": they minimise time-averaged squares, not suprema | reworded |
| R7-11 | minor | constants are suprema over [0, 300]; on the Lorenz regressor the suprema over [0, 1500] are larger by up to 14% (p₁ 5.31 → 5.67, p₂ 19.4 → 20.7, π₂ 49.4 → 51.9, ‖N̄₂‖ 0.42 → 0.48 at T = 20), the rate at k = 100 lower by 3%, the certificate intact; unchanged to four digits on the quasi-periodic regressor | stated in the numerics paragraph |
| R7-12 | cosmetic | "the regressor traversed at standard Lorenz speed has ρ_DD ~ 10⁶–10⁷" referred to an exploratory instance with other components (tanh(x_i/8)) | attributed as such |
| R7-13 | minor (code) | `ExpAverageLyapunov.constants` docstring stated p₂ = ‖M − M̄‖ + T‖F̄ − F̂‖ for explicit=True; the code uses the looser Ψ_M² + T(sup‖F̄‖ + ‖F̂‖) (valid) | docstring corrected |
| R7-14 | cosmetic (code) | `exact_stationary_covariance`: stray "Ẅ := Ẇ" in the docstring; the returned C0 is the end-of-period covariance | docstring fixed; the key C_end added (C0 kept) |
| R7-15 | minor (experiment; fixed during the review) | `exp_pe_general.py` as first read integrated the Lorenz orbit only to Lorenz time 101 while the padded grid for T = 20 reaches 108, so the spline extrapolated and the T = 20 constants were garbage (Λ = 2233, p₂ = 9·10¹⁴, all conditions failing); fixed before the published table was produced (its T = 20 row reproduces exactly) | `CubicSpline(..., extrapolate=False)` and a range assertion in `Psi`, `Psid` |

## Round 8 (research continuation and adversarial review, 2026-10-19)

**Design.** The fourth research increment adds Theorem D‴ (uniform-in-time propagation of chaos in the learning regime for the
averaged slow loop at τ = 0, §5d) with Proposition D‴.5 (the weight-sensitivity law of the stationary coupled distance at
finite m), `src/core/chaos_learning.py`, `tests/test_chaos_learning.py` and `experiments/exp_chaos_learning.py`. One referee was
run on the theorem's statement, proof and remarks, the code and the tests before the numerics existed
(`review/REVIEW_round8_Dppp.md`). The confirmed issues and their fixes follow; the numerics were designed after the report.

**Verdicts.** Theorem D‴ holds: every constant of Steps 0–9 was re-derived (the two Lipschitz facts with the ½ Taylor factor,
the telescoping identity, (2.1) with the segment inside the neighbourhood and Weyl's inequality, (3.1) with c²/4, (4.1), the
κ/8 split and the Young constant, the Grönwall exponents, the two thresholds, McDiarmid with c_j = 2K_Φ/m for the vector-valued
L²(ρ) statistic, the grid counts), an adversarial search over 3000 configurations found no violation of (2.1), (3.1), (4.1)
(the interaction term reaches 2% of its bound), and an independent brute-force implementation (m = 2, one quadrature point)
confirmed the three-term identity to 10⁻¹¹. The code holds (drift = −G∇J_m to 5·10⁻¹¹; the residual-weighted curvature is the
frozen-residual Hessian, differing from the full Hessian block by exactly M₂/m; the searched constants L_ϑ^G = 4.0435, L₂^G =
11.1111 are the true suprema on a 10⁸-point grid). Three test assertions were vacuous and the fifth test's premise was wrong
(a transient, not a stationary state; the hypothesis violated at its own ϱ′): the tests were rewritten and the numerics were
redesigned around a verified instance (a three-atom stationary measure at which (H) holds with ϱ′ = 0.1).

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3-1 | minor (premise) | (H) fails for the ϱ′ = 0.3 of the first warm-start test: the supremum of κ_R over N_ϱ′ exceeds σ₀ (0.0145 > 0.014); (H) holds there only for ϱ′ ≲ 0.27 | instance changed to σ₀ = 0.15σ_crit, where (H) holds at the three-atom stationary measure with ϱ′ = 0.1 and margin 0.0020 |
| D3-2 | minor (code) | `kappa_R_neighbourhood` was a Monte Carlo lower estimate (0.0125 against the dense 0.0145), unable to certify (H); the supremum sits on the boundary sphere | deterministic Fibonacci-sphere search on the boundary (and an inner shell), labelled a lower estimate; the test checks it against the first-order expansion λ_max + ϱ‖∇λ_max‖ |
| D3-3 | minor (tests) | `D.max() < C₁/m` compared 0.22 with 3·10⁶ (vacuous: C₁ = 5·10⁷, m₀ = 8·10¹⁸); D was increasing over the whole run; the (c)-check used K_Φ²(1/m + 1/M) where the shared-particle law is ⟨Var⟩(1/m − 1/M) | the vacuous assertions are kept only as sanity checks and labelled as such; the saturation and 1/m statements moved to the experiment (horizon 900 ≫ 1/κ); the (1/m − 1/M) law is now Proposition D‴.5(ii) |
| D3-4, D3-5 | cosmetic (tests) | `kappa_R_neighbourhood(...) ≥ λ − 10⁻¹²` true by construction; the C₁ check a tautology | replaced by the first-order-expansion check of the search and by the four constraints of the bootstrap (constants) |
| D3-6 | minor (code, anti-conservative) | `Lambda1_from_profile` integrated only to the last grid point ≤ s₁ (0.050 for 0.055); Λ₁ enters as e^{2Λ₁} | the integral runs to s₁ with a linearly interpolated value; tested |
| D3-7 | cosmetic | the closed form L₂^G ≤ λ_G a_max/s_min² is loose by λ_G/G_c (44.4 against the true 11.11 = G_c a_max/s_min²) | stated in the docstring and the table |
| D3-8 | cosmetic | the module docstring's T_A bound omitted the L₂^GY term | fixed |
| D3-9 | cosmetic (block) | the good event was defined with sup_{s≥0} while (a) bounds a horizon S; (c)'s "for all time" | 𝒢_S per horizon in Steps 6–8 and in (c); (R_eff)₊ |
| D3-10 | cosmetic (block) | "bounded vector field" is false (linear leak); a garbled fragment in the definition of R | well-posedness sentence rewritten (linear growth, Lipschitz; bounded on the invariant ball); fragment fixed |
| D3-11 | minor (Remark D‴.2) | "on a ridge 𝒱 has negative curvature, so (H) fails for every s" is false pointwise (counterexample 𝒱 = (x² − 1)² + y² + 10x²y²: the invariant separatrix {x = 0} is convex for |y| > 0.45); the filament's points converge to saddles, so (H) fails for every s₁ asymptotically | remark rewritten |
| D3-12 | cosmetic (Remark D‴.1) | the "non-strict form for every Wasserstein-stable stationary measure" needs atoms of positive mass and yields non-minimality of 𝒢; the term was undefined | stated for atomic measures with the splitting argument (Proposition K) |
| D3-13 | minor (premise) | the first test's warm start was a transient (fit error 0.037 → 0.26 over s ≤ 40): "κ_R < σ₀ along the run" said nothing about (H) for s ≥ s₁; no instance with (H) verified at a near-stationary state existed | the three-atom stationary measure of J₃ (zero drift) is computed, (H) is checked there by the sphere search for ϱ′ ∈ {0.05, 0.1, 0.2, 0.3}, and the warm start is the jittered mixture around its atoms; the hypothesis is tracked along the run with the stationary residual (the proof's bound (0.2) separates the reference's sampling fluctuation, which is asserted in the test) |
| D3-14 | cosmetic (test) | `T_B < 0` is configuration-specific (T_B = +1.67 found at D = 3.9, m = 1, consistent with (3.1)) | the test asserts T_B/(−‖h‖²) → 1 as δ → 0 and (3.1) at every scale |

**A finding of the redesign (recorded here because the referee's D3-3 pointed at it).** At m = 16–256 the particles of the
m-system sit 0.4–1.4 whitened units from their reference copies, outside the ϱ′-neighbourhood in which 𝒱[μ_Θ̄] is convex,
at a stable stationary point of their own J_m: the theorem's per-particle localisation is in force only for m ≥ m₀ = 16K_Φ²/θ²
(10²³ on the instance), and the finite-m distance is set by the clusters' weight fluctuation times the atoms' sensitivity to
the weights — Proposition D‴.5, whose prediction the per-cluster displacements of the run follow (cosine ≥ 0.8, norms within a
factor 2 at s = 60 ≪ 1/κ). The theorem's constants are reported as worst-case envelopes; what the numerics test is the hypothesis,
the mechanism (saturation, 1/m) and the sharp constant.

## Round 9 (research continuation and adversarial review, 2026-10-19)

**Design.** After round 8 the Theorem D‴ block gained Proposition D‴.5 (the weight-sensitivity law of the stationary coupled distance
at finite m), Remark D‴.6 (the theorem's regime), Corollary D‴.7 (the real loop on a warm start through Theorem J), Remark D‴.8, the
numerics of `exp_chaos_learning.py` and the readings. One referee was run on these (proofs, the experiment script, the new functions
of `chaos_learning.py`, the fifth test, every number and every reading against the raw data; `review/REVIEW_round9_D3b.md`). The
confirmed issues and their fixes follow.

**Verdicts.** Proposition D‴.5 holds with fixes: the algebra was confirmed independently (∇_{ϑ_k}𝒢_w = w_k∇𝒱 to 2·10⁻¹¹; the replicated
weighted stationary point a critical point of J_m to 10⁻¹⁶, also at weights far from uniform; the Hessian block splitting — separation
blocks σ₀I − R^G_k with multiplicity n_k − 1 and the whitened, weighted K-atom collective block — reproduces all 30 eigenvalues of J₁₀
to 8·10⁻¹²; S by the implicit function theorem against central differences to 3·10⁻⁷; the multinomial covariance (1/m − 1/M)C̃ by
simulation), and the exact first equality of (ii) holds on the data to 0.1–2.3% in all 12 (seed, m) pairs. Corollary D‴.7 and Remark
D‴.8 hold with fixes (Theorem J used correctly; the two plant states handled correctly). Numerics and readings hold with fixes for
m = 64 and 256; at m = 16 one of four warm seeds sits on a saddle (the major issue below). The experiment script and the new code hold
with minor fixes; the fifth test passes and tests what it says, with two near-vacuous assertions and one tautology.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3b-1 | **major** (reading) | at seed 3, m = 16 the realised weights (9, 4, 3)/16 make the replicated configuration a saddle of J₁₆ (the splitting condition fails at atom 1: κ_R = 0.00927 > σ₀ = 0.00707; eigenvalue −0.00219 with multiplicity 8); the m-system sits on it at s = 900 (D = 1.2917 against the saddle's exact stationary distance 1.2916) and leaves it later (the referee continued the autonomous run: D = 1.95, 4.70, 5.13 at s = 1800, 2700, 3600 as two particles migrated and the weights re-sorted to (7, 5, 4)/16); the reported m = 16 "saturation" is, in that seed, a transient plateau; Remark D‴.6's "at a stable stationary point of their own J₁₆" was false for it; generically the splitting condition fails at ‖w − w̄‖ ≈ 0.06–0.36 against a root-mean-square weight fluctuation 0.20 (m = 16) and 0.10 (m = 64): admissible in ≈ 52% of multinomial draws at m = 16 and 76% at m = 64 | Proposition D‴.5(i) states the admissible set W_adm and the saddle alternative, (ii) is conditional on admissibility, new (iv) bounds P(w^m ∉ W_adm) by Hoeffding; Remark D‴.6 corrected; Remark D‴.2(iii) names the saddle mechanism; the experiment reports the splitting margin of Θ*(w^m) and the exact stationary distance per seed (new columns), and the reading says which seed is inadmissible; the test checks the margins of (7,5,4)/16 and (9,4,3)/16 |
| D3b-2 | minor | "invertible" K-atom Hessian; the minimiser claims need positive definite | stated |
| D3b-3 | minor | the O(‖dw‖³) remainder holds with S at w^M; with S at w̄ (as coded) there is an O(‖dw‖²‖w^M − w̄‖) term | (ii) written with S(w^M), the replacement by S(w̄) priced |
| D3b-4 | minor | the corollary's smoothness argument compared maps on ℝ^{mp} and ℝ^{Mp}; the O(k_e⁻²) remainders were dropped; κ clashed with Theorem J's | argument via permutation equivariance and the block splitting (the offset depends on the K-atom data); remainders kept; κ_J |
| D3b-5 | minor | the expectation statement needs conditioning on admissible weights (probability ≈ ½ at m = 16) | conditioned |
| D3b-6 | cosmetic | "Checked on the instance" overstated: the corollary's saturation is not checked | "Consistent with …; not checked" |
| D3b-7 | cosmetic | "−0.196 at s = 1" is the first record (s = 0.02); at s = 1 the margin is −0.179 | "from the first record" |
| D3b-8 | cosmetic | "seven to eight orders of magnitude" → 8.4–8.8 | corrected |
| D3b-9 | minor (reading) | the flip fraction and the medians decay like m^{−0.64±0.08} (threshold 3: −0.73), not "the m^{−1/2} law"; the cold learned-function slope −0.66 ± 0.06 is 2.3 s.e. from −0.5 | the experiment fits the exponents; the texts say "steeper than m^{−1/2}" |
| D3b-10 | minor (support) | the real-loop claim compared different clouds (ratios 0.78–2.58); on the same clouds D_real(t)/D_avg(t/k_e) ∈ [0.985, 1.016] for all seeds and m | the experiment reproduces the real loop's clouds, runs the averaged flow on them and reports the per-seed ratios |
| D3b-11 | minor (code) | `weighted_stationary_atoms` without continuation jumps branches at weights far from w̄ (seed 3's weights gave D = 19 instead of 1.29; 77 of 300 plain solves in a multinomial sample) | the margin and the exact distance use continuation from w̄ (`splitting_margin`, `stationary_distance_exact`); `weight_sensitivity` (small steps) is unaffected |
| D3b-12 | cosmetic | κ halved without a stated reason; S̃_k undefined; "within 5–20%" → 0.8–20%; "0.245 throughout" → 0.245–0.248; the heavy-tail explanation does not fit m = 256 | fixed in the texts |
| D3b-13 | cosmetic (test) | the (0.2) and (0.1) assertions have slack ≥ 10²; `sigma0 < SIGMA_CRIT` a tautology; a docstring mismatch | tautology removed, slack stated, docstring corrected, a splitting-margin check added |
| D3b-14 | minor (code) | the recorded κ_R at the m-particles used the reference's residual and cannot see the m-system's own instability (0.0053 against 0.0101 with its own residual) | `coupled_run` records both (future runs) |

## Round 10 (adversarial review of Proposition D‴.9, 2026-10-19)

**Design.** One referee on the material new since round 9: Proposition D‴.9 (i)–(iii) and its Status, Remark D‴.10, the D‴.9 numerics
block (tables 1–4, Reading), Remark D‴.11, `src/core/prior_design.py`, `experiments/exp_prior_design.py`, `tests/test_prior_design.py`
(`review/REVIEW_round10_D3c.md`). Every number of tables 1–4 and of the prose was recomputed from the raw data (`prior_design_raw.npz`,
`learning_raw.npz`); the regeneration by `--from-raw` was byte-identical; the two inadmissible m-systems were re-simulated autonomously
to s = 6000 and 9000.

**Verdicts.** (i) holds with one major fix (a hidden hypothesis); (ii) holds; (iii) holds with fixes; Remark D‴.10, the numerics and Remark
D‴.11 hold with fixes (two attributions wrong on the instance); the code holds; the experiment holds with minor fixes; the tests pass with
one vacuous assertion. Independent confirmations: grad_J against finite differences (5.6·10⁻¹⁰), drift = −G∇J, the single-atom reduction
(8.9·10⁻¹⁶), the weighted energy gradient (6·10⁻¹²), the function sensitivity, all rows of table 1, table 2 from scratch, table 3's
windows and slopes, the per-seed laws, the K-atom Hessian at the exact representation (eigenvalues 0.0074–0.20), the exact law in all nine
pairs (0.987–0.998, not yet reported), and the saddles' fate.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3c-1 | **major** (statement/Status (i)) | the no-Voronoi-boundary hypothesis was imposed for s ≥ s₁ only, while Step 2's (2.1) and (5.2) are used on [0, s₁); a pair straddling a boundary adds (σ₀/m)⟨δ_i, ϑ⁰_{j(ϑ_i)} − ϑ⁰_{j(ϑ̄_i)}⟩_{G⁻¹} ≥ 0 to ½dD/ds (V_ν semiconcave), up to σ₀·diam(supp ν)·D^{1/2}, confirmed numerically to all digits; also needed for the well-posedness of (𝔄₀^ν) | hypothesis required for every s ≥ 0; the mechanism stated in Status (i); Theorem J's ρ₀ below the cell distance |
| D3c-2 | minor (texts ×4, results md) | "the shrunk prior keeps all four admissible with margins at or above σ₀" false for three of four (+0.0054 to +0.0066 < σ₀); "undershoots everywhere" false (the over-weighted atom overshoots by 1–3%) | `report()` computes f̂/f at the centres and the margins' position relative to σ₀; texts corrected |
| D3c-3 | minor (Remark D‴.11) | the centre-overshoot mechanism is exact for the matching prior only: the σ-modification's inadmissible seed undershoots at its centre (unstable direction amplitude–width, the tails), Proposition K's R₂ is carried by the overshoot over the negative bump and the undershoot beyond the inflection points, not "between the bumps"; "loses definiteness" → "loses margin"; the admissible set is asymmetric | attributions corrected; the asymmetry computed by `report()` (excess on one atom alone) |
| D3c-4 | minor (Numerics) | "plateaus being saddles that a longer horizon would see left": the saddles are left by Proposition K's local split inside the cell, weights unchanged, D on its plateau (−1.3%); the ratios 0.91/0.81 are the linearisation remainder (exact law 0.987–0.998) | `report()` computes the exact law for the nine pairs; the saddles' fate described; Proposition D‴.5(i)/(iv) wording ("migrating particles") widened to both mechanisms |
| D3c-5 | minor | "κ_R(r_ref) at most 0.00488" was the seed-mean curve's maximum (per seed 0.00634, 0.00674, 0.00155); fit and drift endpoints were seed means | `report()` prints per-seed maxima and labels the seed means |
| D3c-6 | minor ((iii)) | "both sensitivities grow with the amplitudes" is a trend of the amplitude rows, not a law | moved to the numerical findings, qualified |
| D3c-7 | cosmetic | "equivalently the gap exceeds 2ϱ′" is sufficient, not equivalent (distance to the boundary ≥ gap/2) | "in particular when"; `cell_margin` docstring and table header |
| D3c-8 | cosmetic (code) | module docstring said the matching prior "lowers S" | corrected (it raises S) |
| D3c-9 | cosmetic (tests) | test 2's T_B bound vacuous (T_B < 0 against a positive bound); test 3's fm < 4fs loose | test 2 asserts the content of (3.1), \|T_B + ‖h‖²\| ≤ (L₂^G/2)‖h‖D; test 3 brackets the ratios (1.4–1.9, 1.9–2.5) |
| D3c-10 | cosmetic | "3.6 times faster", "relaxation 141 against 511/512" are the theorem's 1/κ, not observed times | `report()` measures the late-time e-folding times (48–84 against 107–306 slow units, medians 59 and 195, factor 3.3); texts say "the theorem's 1/κ"; 511 |
| D3c-11 | cosmetic | q_j plays no role in the nearest-atom leak; κ_R[μ*] ≡ 0 for every set; IFT non-degeneracy unstated | stated in (i), (ii) and Status |
| D3c-12 | cosmetic | E D_∞ quoted as an equality | conditional form with the remainders |
| D3c-13 | cosmetic (table 4) | the Jensen gloss explains only part of the ratio 0.48 | `report()` adds the RMS over seeds (0.84, 0.79, 0.52) and the per-seed law from the realised weights (0.60–1.57) |

**Lesson.** A hypothesis imposed "from s₁ on" must be checked against every step of the proof that uses it, including the transient
estimates feeding Λ₁; and the convention of a diagnostic (`decomposition`'s same-cell leak) is not the identity it stands in for once its
hypothesis fails — the true derivative should be recorded alongside.

## Round 11 (adversarial review of Proposition D‴.12, 2026-10-19)

**Design.** One referee on Proposition D‴.12 (i)–(v), its Status, Remark D‴.13, the D‴.12 numerics block, the §0/§7/§8 clauses, the RESULTS
and LITERATURE entries, `src/core/ot_leak.py`, `experiments/exp_ot_leak.py` → `results/chaos/ot_leak.md`, `tests/test_ot_leak.py`
(`review/REVIEW_round11_D3d.md`). Every number recomputed from the raw data (regeneration byte-identical); the design table recomputed from
scratch; the plan checked against the transportation LP (equal to 5·10⁻¹⁵, integral when commensurate, a forest otherwise); Danskin to
3·10⁻¹⁰; the Hessian of J₁₅^{OT} equal to the nearest-atom J₁₅^ν's to 1.7·10⁻¹⁸ with the predicted spectrum; the incommensurate runs
reproduced from their seeds and continued to s = 2400.

**Verdicts.** (i) holds; (ii) holds as a conditional statement whose hypothesis is incompatible with Theorem D‴'s i.i.d. start except on an
event of probability ≈ 0.83/m (major); (iii) holds; (iv) holds as stated (the count proved, the bounds formal with vacuous constants; the
orders confirmed by the right estimators); (v) holds except one false clause; Status, Remark D‴.13, the numerics and the readings hold with
fixes; the code holds; the experiment script holds with minor fixes; the tests pass and test what they say, with one loose tolerance and
one tautology.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3d-1 | **major** ((ii), Status, §0, §8) | the frozen-plan hypothesis ("the plan assigns each particle to the atom of its reference's Laguerre cell for every s ≥ 0") forces w^m(0) = q: under the i.i.d. mixture start P(balanced) = 0.053, 0.014, 0.0034, 0.0016 (m = 15, 60, 240, 510) → 0; in all 32 OT systems the plan at s = 0 assigns 1–21 particles away from their nearest atom and the shared-particle coupling never re-matches them (D at 0.97–0.98 of the mismatch plateau at s = 600) | (ii) restated for a balanced sample, with the probability; Step 8's McDiarmid noted as written for i.i.d. particles; the coupling statement labelled formal; the generic case stated as open (a plan-following coupling as the candidate); §0, Status, §8 corrected |
| D3d-2 | minor ((v), Remark D‴.13) | "after the transit all three vanish" false for the synchronous D, which converges from below to the permanent-mismatch plateau Σp_{kk′}Δ²_{kk′} (0.95–0.96 of it at s = 600) | (v) and D‴.13 corrected: the memory moves from the measure to the coupling; `report()` prints the plateau ratio |
| D3d-3 | cosmetic | "one to twelve particles" → 1–21 per system; fit errors are seed means | per-system counts and labels |
| D3d-4 | minor | the fit-error explanation by the jitter's spread 0.02e^{−σ₀600} ≈ 3·10⁻⁴ is off by two orders: the migrated particles sit at Δ_{kk′}e^{−σ₀600} ≈ 0.09 from their atoms and the finals are not stationary (drift 5.7·10⁻⁴–1.1·10⁻³); h² decays at 3.0–3.8σ₀ | `report()` prints the migrated particles' distances and the final drifts; text corrected |
| D3d-5 | minor | the incommensurate m-laws (−2.07 ± 0.05, −1.21 ± 0.04) were measured against an M = 512 reference carrying its own bridge particles; against μ*(q) = ν the laws are m·W² = 19.3–20.9 (slope −0.98 ± 0.03) and m·fit = 0.45–0.60 (squared slope −1.81 ± 0.05) | `report()` adds the distances to ν; texts quote them as the test of (iv) |
| D3d-6 | minor | the m = 16 function-level value (0.82 per window) was called stationary; continued runs give the plateau 2.82·10⁻⁴ (21% below), ratio OT/nearest 0.114 | text says "still relaxing" and quotes the continued run |
| D3d-7 | cosmetic | W₂ written for the squared distance; (iv) compared the W² floor with the wrong law (the nearest-atom W² floor is O(m^{−1/2})) | W² throughout; "half an order in W², a full order at function level" |
| D3d-8 | cosmetic (code) | `transit_mask` undefined on a tied split row (the centroid bridge particle counted as in transit for ever); header "= ½|w(0) − q|₁" exact on commensurate sizes only | split rows excluded; the table recomputes the unsplit particles in transit from the clouds; the recorded s_tr kept for the commensurate sizes only |
| D3d-9 | cosmetic | "at m = 16 the two or three particles to move are the split ones" (seed 3 migrates two unsplit ones) | corrected |
| D3d-10 | cosmetic | "about 0.5δ" → 0.69δ–0.49δ | `report()` prints fit/δ |
| D3d-11 | cosmetic | "within 5–14%" → −7% to +28% of the formula; the exact K × K transport between the weighted stationary measures gives 0.997–0.999 at admissible weights | `report()` prints both; texts corrected |
| D3d-12 | cosmetic | (iv)'s constants vacuous on the instance (4·10⁶ against 1.0–1.4; 1.00 against 0.018); measured m·W² ≈ 20 = 0.51Δ_ν², m·fit ≈ 0.55 | stated in (iv) |
| D3d-13 | cosmetic (tests) | test 2's tolerance slack 1.4·10⁴; test 1's W² identity tautological; dead code | tolerance 10⁻⁸; the plan's cost checked against the LP; dead code removed |
| D3d-14 | cosmetic | the w̄-law column is 1.6–2.3× the seed means; the exact per-seed function-level law gives 0.999–1.009 | `report()` prints the exact per-seed law |
| D3d-15 | cosmetic | Laguerre cells not unique for a disconnected support; the M-reference is a poor stand-in for the mean-field flow during the transit | both stated |

**Lesson.** A coupling hypothesis must be checked at s = 0 against the law of the initial data, not only along the flow: the plan's mass
constraint and an i.i.d. sample are incompatible with probability 1 − O(1/m), and a stationary statement about the particle system does
not transfer to a dynamical one without a coupling that respects the constraint.

## Round 12 (adversarial review of Corollary D‴.14 and Proposition D‴.15, 2026-10-19)

**Design.** One referee on Corollary D‴.14 (i)–(iii), Proposition D‴.15 and their Status and numerics, the rows and bullets mentioning them,
`src/core/lifted_orbit.py`, `experiments/exp_ot_real_loop.py` → `results/chaos/ot_real_loop.md`, `experiments/exp_lifted_orbit.py` →
`results/chaos/lifted_orbit.md`, and the three real-loop tests (`review/REVIEW_round12_D3e.md`). Every number of both md files reproduced
(regeneration byte-identical); the block structure of Proposition D‴.15 confirmed at the field level (leakage 2·10⁻¹⁶) and at the monodromy
level (‖BᵀDP B − M_k‖ ≤ 2·10⁻⁷); the ghost monodromy against the matrix ODE (3·10⁻⁹) and against exp(T₀Ā_k) (7·10⁻⁹); the k_e-scaling of the
prediction's error checked at k_e = 50.

**Verdicts.** D‴.14(i) holds; (ii) and (iii) hold with fixes (a garbled proviso, hidden hypotheses); the Status is honest about the
uniqueness ball; the D‴.14 numerics reproduce but their reading fails in part (major: three saddle weights unmentioned; the tracking-error
comparison against a scheme floor); D‴.15 holds, with its last claim true but incompletely proved as written; the D‴.15 numerics reproduce
with two ranges wrong; the code holds; the tests pass with two weak assertions.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3e-1 | **major** (D‴.14 numerics, §7, RESULTS) | three of the six nearest-atom real-loop systems have inadmissible realised weights (seed 0 m = 15, margin −0.00348; seed 1 m = 15, −0.00108; seed 1 m = 60, −0.00216): their lifts are saddle orbits and the m-systems are not on them (the over-weighted class splits along the centre direction, spread growing); "every size on the orbit of its own weights" and "within 3.3·10⁻⁶ of its own lift" were false for them, and the two ratios 0.960 were the splits | `report()` prints the splitting margin and the within-class spread of every system and names the saddles; (iii) requires w^m, w^M ∈ W_adm; the text separates admissible pairs from saddles |
| D3e-2 | minor ((ii)) | the proviso "k_e ≥ √m·sup‖Θ^q_p − Θ*‖_w/ρ₀(m)" was off by a factor k_e; ‖·‖_w undefined | restated as √m·sup‖·‖_w ≤ ρ₀(m), i.e. k_e ≥ √m C_q/ρ₀(m) with Theorem J(a)'s C_q/k_e; ‖·‖_w defined |
| D3e-3 | minor ((ii), (iii), D‴.15) | hidden hypotheses: (i)'s cell condition along Γ_q; admissibility of both weights in (iii); "targets constant along Γ" → "locally constant near Γ" | stated |
| D3e-4 | minor (D‴.15 proof) | the identification of M_k's exponents through Theorem J(b)'s multiset does not attribute eigenvalues to classes | proof replaced by first-order averaging of the ghost equation (e_p = r/k_e + O(k_e⁻²)), which needs only that Γ_K is Theorem J's orbit of (CL_w) — no condition on m; the Status's conditional removed |
| D3e-5 | minor | "within 0.6–2.4%" is 0.56–3.22%; "within 0.3–3%" is 0.23–4.07% | `report()` prints the relative errors; texts corrected |
| D3e-6 | minor | "m = 60 and 510 within 3·10⁻¹¹" false for seed 1, m = 60 (a saddle) | text corrected |
| D3e-7 | minor | the OT loop's tracking error 3.9·10⁻⁵ and fit error 2.1·10⁻⁴ were the zero-order-hold floor of the integration scheme (the continuous orbit has e ≡ 0; the floor ∝ dt) | the plant's RK4 re-evaluates the control at the stages (Θ frozen over the step) in both real-loop experiments, which are rerun |
| D3e-8 | cosmetic | the collective prediction labelled λ(G^{1/2}∇²𝒢_qG^{1/2})/k_e is the weight-whitened block | label H_q, the collective block of Proposition D‴.5(i) |
| D3e-9 | cosmetic | "to four digits" is 10⁻⁴ relative; a round-before-format artefact | relative errors printed; grouped display |
| D3e-10 | cosmetic ((i)) | "½‖w − q‖₁ … leave" is a lower bound (a chain of reassignments can move more) | "at least" |
| D3e-11 | cosmetic (§8) | "Theorem J's orbit is the lift" stated unconditionally | "under its uniqueness-ball proviso" |
| D3e-12 | cosmetic (tests) | test 1's "atoms moved" compared particles with particle 0; test 2's W² > 10⁻³ held at t = 0 | motion measured against the start; the two lifts' displacements compared |
| D3e-13 | cosmetic | the ghost exponent at the matching masses is exactly σ₀/k_e; "−G∇𝒢_w" is the gradient in the metric diag(w) ⊗ G⁻¹; "within 10⁻⁸" are squared distances | corrected |

**Lesson.** A comparison of tracking errors must first establish the integrator's floor (here the zero-order hold of the control over the
step), and a run that is meant to sit on an orbit must check that the orbit is stable at the realised weights — the admissible set of
Proposition D‴.5 is where it decides.

## Round 13 (adversarial review of Remark D‴.16, the two-dimensional numerics, Proposition D‴.17, and the fixes of rounds 11–12, 2026-10-19)

**Design.** One referee on Remark D‴.16 and its dithered numerics (`exp_dither_ot.py`), the paragraph on the planar plant (`exp_ot_2d.py`),
Proposition D‴.17 with its test and references [S43], [S44], and a verification pass on every fix of rounds 11 and 12
(`review/REVIEW_round13_D3f.md`). Every number recomputed from the raw data (both md files regenerate byte-identically; single runs reproduced
from their seeds); the ridge and exchange barriers of the dithered instance computed by constrained minimisation; the excursions of the OT
runs classified particle by particle; the planar OT run continued to 10/σ₀; the pure-leak plan stability tested on random instances.

**Verdicts.** Remark D‴.16: the nearest-atom half holds (formal) with fixes, the OT half fails (major). The dithered numerics reproduce; several
readings fail. The planar numerics hold with cosmetic fixes. Proposition D‴.17 (i)–(ii) hold; (iii) holds as stated with a false gloss; its
Status fails in one sentence. The test passes and discriminates (a control 10⁻³ off the tie fails). The round-11/12 fixes are in and correct
except D3e-12, which was logged as applied but was not.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3f-1 | **major** (Remark D‴.16, its numerics, `dither_ot.md`, RESULTS, §8, docstring) | the OT half — "J^{OT} charges a misplaced particle ≈ σ₀Δ_ν²/2 = 0.14, so the counts are q up to exp(−σ₀Δ_ν²/(2τ₀)); a hop is answered by a swap of assignments and a transit" — is false: the plan keeps a particle that crosses a ridge assigned to its own atom (an exchange needs two particles to pass each other), so the counts change at the nearest-atom ridge barrier (0.040 at the lowest ridge, J^{OT} = J^ν there); the off-balance fractions 6.3·10⁻³, 0.117 are 6700 and 121 times the claimed order (effective barrier 0.058); the excursions are returns of the crossing particle (121 of 121 at τ₀ = 0.01, 298 of 336 at 0.02), and the plan never changes at τ₀ = 0.01 | the mechanism rewritten (a crossing is reversible, the barrier unchanged); `report()` classifies the excursions (returns/exchanges) and computes the off-balance fraction and its barrier; RESULTS, §0, §7, §8, docstring |
| D3f-2 | minor | the "leading-order Gibbs law multinomial × energy" is not an order (the Laplace factors are of the same order); the agreement 0.54/0.58 cannot discriminate (N_eff 9–27 per seed, SE ≈ 0.1); at τ₀ = 0.02 the law is 2.1 SE above; the law under-weights the saddle count vectors 2.7–3.2 times (TV 0.20 against 0.07 between seeds); "min J^ν = m𝒢*_w" needs w ∈ W_adm | "energy–multinomial law", SEs and TV quoted, the saddle excess attributed to the soft and unstable modes, "for n/m ∈ W_adm" |
| D3f-3 | minor | the barrier 0.039 rests on one event; the "hop rates" count record-level label changes, about half undone at the next change | "record-level label changes, about half undone"; 0.048 (0.043–0.052) among the ridge minima 0.040–0.074 |
| D3f-4 | minor (D‴.17 Status, §7) | "under the pure leak a plan never switches once optimal" holds for integral plans only (split plans switch: 216 of 400 random incommensurate instances); with the fit, cold starts at σ₀ switch the plan (2 of 6 at m = 15, 6 of 6 at m = 16); the test's run never switches | restricted to commensurate m; the switches reported; the test's run described as a smooth gradient flow |
| D3f-5 | minor (D‴.17(iii)) | the gloss reads as a classification of the limits; 𝒢_q has a second local minimum whose replication is a saddle of J^{OT}_15 although q ∈ W_adm; "J_m" should be the frozen energy J^π_m | reworded; the example stated |
| D3f-6 | minor (tests, log) | the round-12 test fix D3e-12 was logged but not applied (test 1 compared particles with particle 0; test 2's W² > 10⁻³ held at t = 0; docstrings cited (a)/(b)) | applied: motion against the start, the two lifts' class-atom displacements compared, docstrings (i)/(ii) |
| D3f-7 | cosmetic | "contracted clouds" for the split saddle clouds; `ot_real_loop.md` §2 title; "to four digits" in RESULTS; "from below" next to a ratio 1.116 in `ot_leak.md` | corrected (the md files regenerated) |
| D3f-8 | cosmetic | "three hundred to a thousand times" is 285–480 per system; "40 to 2000" is 126–548 per system; the planar OT distance at 5/σ₀ is a stopping-time value (1.8·10⁻¹² at 10/σ₀); the planar docstring said "control held over the step"; numpy reprs in two md headers | corrected; `report()` prints the per-system ratios |
| D3f-9 | cosmetic (D‴.17) | the change of variables Y = G_m^{−1/2}Θ unstated; ℝ_exp (or ℝ_an) suffices; Łojasiewicz's inequality attributed to [S44]; "Lipschitz" for a locally Lipschitz potential; the test docstring's "Clarke-subgradient inclusion" | stated and corrected |

**Lesson.** A barrier argument must follow the particle that actually moves: under the mass constraint the plan does not follow a single
crossing particle, so the relevant barrier is the crossing's, and the constraint acts on what happens after the crossing (reversibility), not on
its probability. And a fix logged in the review log must be checked in the repository before the round is closed — a patch script that stops
at its first failing anchor silently skips everything after it.

## Round 14 (adversarial review of Theorem D‴.18, Remark D‴.19 and their numerics, 2026-10-19)

**Design.** One referee on Theorem D‴.18 (statement, proof, Status), Remark D‴.19, the numerics paragraph, the §0 and §7 rows, the §8 bullet,
the RESULTS section, `exp_uniform_basin.py` with its raw data and `tests/test_uniform_basin.py` (`review/REVIEW_round14_D3g.md`). Every
number recomputed from the raw records (the md regenerates byte-identically); the batched integrator compared with `exp_ot_real_loop.real_run`
on single systems (4.4·10⁻¹⁶ over 10⁴ steps); the initial configurations rebuilt from the seeds; the three structural facts checked on the
averaged loop at three replicated configurations (m ≤ 160) and on the real loop's one-period map (m = 6, 30, 150); the 3-atom orbits shot
and their ghost exponents computed.

**Verdicts.** Theorem D‴.18 holds with a radius fix in (a); (b) and (c) hold; the m-uniformity is sound (the norm is equivalent to the
per-particle max-norm with constants 2 and 1 + w_min^{−1/2}, and every map of Steps 1–4 is of mean-field type). Remark D‴.19 holds in part.
The numbers hold; the reading of the runs across a cell boundary fails. The tests pass and are not vacuous, but are weaker than their names.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3g-1 | **major** (D‴.18 numerics, RESULTS, §0, §7) | "at ρ = 4.5 … the nearest-atom leak keeps them there … the lifts of other weights, as Remark D‴.19 says" is false on three counts: the counts change during the transient at m = 60, 240 (m = 240 starts at exactly q), no cloud has collapsed and two of the end weights are inadmissible (their lifts are saddles; (2, 0, 1)/3 empties a class), and the runs displace every particle, not one | the reading rewritten from a new cell-bookkeeping table (`report()` prints per system the particles outside their class cell at t = 0, the count changes and the last one, the within-class spreads, the end weights' stationary point, residual and splitting margin, the distance from and W² of their lift, the relabelled particles and the late W² slope); a new Part C moves one particle across a boundary (m = 60, 240) to show Remark D‴.19's mechanism |
| D3g-2 | minor (D‴.18 (a), proof) | the contraction is shown on the |||·|||-ball of radius ρ but uniqueness is claimed in the per-particle ball of radius ρ (which needs radius ≈ 2ρ in |||·|||); the margins must bound C(1 + w_min^{−1/2})ρ; C_N names both the remainder constant and the Jacobian-Lipschitz constant (2–5 times larger on the instance) | uniqueness in the per-particle ball of radius ρ/2; both norm equivalences stated; C_N := the |||·|||-Lipschitz constant of D⟨H⟩ (the remainder then (C_N/2)|||u|||²); ρ := min{λ₀/(4C_N), margins/(C(1 + w_min^{−1/2}))} |
| D3g-3 | minor (Remark D‴.19, §0) | Theorem J's threshold k₀ grows like √m in the unnormalised norm (the real reason for the per-particle norm), omitted; the shifted weights must be admissible (at m = 15 they are not); "the basin at the empirical-measure level is larger" unproved | added; "for m large enough that the shifted weights are admissible" with the m = 15 example; "appears larger in the numerics (open, §8)" |
| D3g-4 | minor (numerics, docstring) | the runs start 12–640 times outside the proved ball, unsaid | Part F computes a sampled lower bound C_N ≥ 0.092 and the proof's ρ ≤ 0.019 (0.015 shifted); the runs are presented as probing beyond it (31–261 times the start radius) |
| D3g-5 | minor (numerics, md) | the m = 3 rates are mode mixtures (the 4.06·10⁻⁴ is the third collective mode); C grows at the horizon for three m = 15 runs; Part B's fits are pre-asymptotic mixtures; the md compared m = 3 with σ₀/k_e | Part F's Floquet and ghost exponents are the test of the rate; the m = 3 final deviations projected on the collective modes; C(T) at four horizons with the growing runs flagged; m = 3 against λ_coll/k_e; the fits "consistent with (1 − κ)μ" |
| D3g-6 | cosmetic | "W² … still falling (the transit)": the transit has ended (last count change t ≤ 2695) and W² falls at exactly 2σ₀; the plan relabels particles | "falling at the leak's rate 2σ₀ (log-slopes −0.01411 … against 0.01415); the transit is over"; the relabelled particles counted; "restores the masses q, not the class map π₀" |
| D3g-7 | cosmetic (tests) | at zero residual the per-particle remainder vanishes and test 2's "same at every m" holds by construction; test 3 would pass at 0.8σ₀; test 2 perturbs two particles | the tests run at the shifted prior (R ≠ 0, λ_sep = 0.80σ₀): the per-particle remainder is exactly N_a + N_b/m with N_a m-independent and non-zero (zero at the matching prior), the separation response is (GR_k − σ₀I)v, the contraction is tested at μ = min{λ_coll, λ_sep} with C ≤ 1.05; names and docstrings corrected |
| D3g-8 | cosmetic | fact (1)'s "only through r_Θ"; "∂²_{ϑ_jϑ_l} vanishes" holds for ⟨H⟩ only; the dependence list omits the prior's margins; "the proviso is void"; the §0 row drops ρ/2 and (1 − κ); the §7 constants for m ≥ 15 are 1.03–1.29; the md prints a rate and C for distances that do not tend to zero | corrected; the md prints "—" for runs that do not end at q |

**Lesson.** A reading of runs started outside a theorem's ball must be checked particle by particle before it is attributed to a mechanism:
the same end weights can come from the initial misplacement or from crossings during the transient, and a lift is a limit only at
admissible weights. And a rate fitted over a window shorter than the inverse spectral gap measures a mixture: the Floquet exponents are the
test.

## Round 15 (adversarial review of Proposition D‴.20 and its numerics, 2026-10-19)

**Design.** One referee on Proposition D‴.20 (statement, Status, numerics), the §0 and §7 rows and the §8 clause, the RESULTS section,
`src/core/entropic_leak.py`, `exp_entropic_leak.py` with its raw data, `tests/test_entropic_leak.py`, [S45], [S46] and the audit bullet
(`review/REVIEW_round15_D3h.md`). Every number recomputed from the raw records (the md regenerates byte-identically); an independent
implementation of J^ε_m, its gradient, the entropic plan and the separation blocks; the critical excess recomputed by continuation in ε with
the Newton residual checked; four runs reproduced from their seeds bit for bit and continued to 3000–6000 slow units; the floors compared with
the nearest-atom leak's on the same clouds.

**Verdicts.** (i)–(iii) hold (the gradient of J^ε_m to 2–3·10⁻¹¹; the separation blocks entrywise to 1.4·10⁻¹², at a size and a quadrature the
tests do not use). (iv)'s rank-one algebra holds, but σ₀I − R^G_k ≻ 0 was a hidden hypothesis, and the reading of a split as a transport of
the excess is contradicted on the instance. (v) holds as a scaling, with loose wording. The critical-excess table holds for δ ≤ 0.1 and fails
at δ = 0.2 (unconverged Newton solves). The dynamics numbers hold, but three readings were 600-unit snapshots of systems still moving.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3h-1 | **major** (numerics, §7, RESULTS, md, `critical_row`, `weighted_stationary_atoms_entropic`) | the δ = 0.2 row (ε_c = 136, p = 0.62, rank-one 88) came from unconverged Newton solves (residuals up to 5·10⁻²), restarted from the exact representation at every ε: the bisection found Newton's convergence boundary, not a sign change | Θ*_ε(w) continued in ε from 256 down (four steps per octave), the bisection started from the nearer converged end, every solve checked to 10⁻⁹ (a RuntimeError otherwise; the solver warns on non-convergence), the branch residuals printed (≤ 8.7·10⁻¹³): ε_c = 65.0, p = 0.574, rank-one 81.9 at δ = 0.2 (5.11, 13.5, 23.1 unchanged) |
| D3h-2 | **major** (numerics, §0, §7, RESULTS, md) | "at ε = 4 the larger sizes stop at the residual excesses 0.017, 0.029, 0.028 — the critical excess at that ε": two of them are unstable weights still splitting at s = 600, the start values came from other seeds, and the critical excess at ε = 4 is 0.016 (0.012 for an excess taken from one atom) | δ_c(ε) by continuation in δ (a new table); the ε = 4, 16, 64 runs continued to 6000 slow units (a new Part 3: per system the end counts, the count changes, the final drift and the margins of Θ*_ε(w) at the weights of s = 600 and s = 6000); the reading rewritten from them: at s = 600 five of the twelve larger systems sat at unstable weights; continued, every system ends at stable weights, the larger ones with residual excesses up to 1/60, 3/240, 7/510, of the order of δ_c(4) = 0.012–0.016 |
| D3h-3 | minor (numerics, §0, md) | "the conditional plan charges mainly two atoms" is false: the outflow splits about equally between the two under-weighted atoms, Cov^G₁ has rank two, and the 22–31% depends on which atom the rank-one formula uses | the md prints the conditional plan, the eigenvalues of Cov^G₁, the exact threshold 2σ₀λ_max(A^{−1/2}Cov^G₁A^{−1/2}) (equal to ε_c) and the rank-one values for both atoms (1.22–1.31 and 0.70–1.22 times ε_c); a second family w′_δ = q + δ(e₁ − e₃) |
| D3h-4 | minor ((iv), numerics, §0, §7, RESULTS) | "the split-off part carries the excess to atom j" is contradicted ((6, 5, 4)/15 at ε = 16 ends at (5, 6, 4)/15: the excess moved, not reduced); "keep their sampled weights" (ε = 16) and "no excess is equalised" (ε = 64) fail on continued runs | (iv): the split-off particles leave toward an atom the conditional plan charges, which need not be an under-weighted one, so a split moves the excess but need not reduce it; the readings from the continued runs: at ε = 16 the three unstable systems move ((6, 5, 4)/15 → (5, 6, 4) at s ≈ 3700, (6, 6, 3)/15 → (5, 6, 4), (25, 17, 18)/60 → (23, 19, 18)) and at ε = 64 (8, 4, 3)/15 → (7, 5, 3); the stable ones keep their sampled weights |
| D3h-5 | minor ((iv), §0) | σ₀I − R^G_k ≻ 0 asserted, not assumed: it fails at (6, 5, 4)/15 and (8, 4, 3)/15 for ε ≤ 4, where the class is a saddle for the fit's own reason; ε_c is implicit; "stable" needs the collective block | the hypothesis stated, with the examples where it fails; ε_c the root of λ_min(S^ε_k) = 0 along the branch (one root per branch on the instance); a local minimum when the collective block and the other classes' blocks are positive too (printed in the md) |
| D3h-6 | minor ((v), §0, §8, audit bullet) | "only the hard constraint is floor-free for every m": on the incommensurate sizes the hard leak keeps an O(m⁻²) floor; "the floor is the nearest-atom leak's": it is 0.23–0.84 of it per seed; "the excess absorbed instead of conserved": it is conserved | (v) and the rows: floor-free on the commensurate sizes only; the floor is Proposition D‴.5's law at the entropic stationary point (a new Part 4: the law at Θ*_ε(q) against the nearest-atom law and the observed floors per seed): the law is 0.55 (ε = 16) and 0.74 (ε = 64) of the nearest-atom law, the observed floors 0.23–0.84 and 0.47–0.91 of the nearest-atom leak's; the excess conserved and partly absorbed by the class biases |
| D3h-7 | cosmetic ((v)) | p used both for the critical conditional fraction and for the sampled mass excess | p_c and δ_c = w_kp_c distinguished; the threshold m ≫ ((1 − q_k)/q_k)(2‖Δ̂‖²/ε)² |
| D3h-8 | cosmetic (Status, §0, [S45], [S46], audit bullet) | "with the roles exchanged" is inaccurate (the roles are the same: the particles as the codevectors, the atoms as the data); "twice the largest eigenvalue" holds only at R = 0; the ε → 0/∞ interpolation and Sinkhorn equipartition have precedents not cited | "with the same roles"; the generalised eigenvalue relative to A_k; [S45] extended (Ramdas et al. 2017, Genevay et al. 2018, Feydy et al. 2019, Asano et al. 2020); the audit bullet rewritten |
| D3h-9 | cosmetic (md, docstrings, (i)) | the md and the module named the wrong solver; "the strictly concave semi-dual" (modulo constants); boundedness is automatic; the ε = 1 distances are stopping-time values; max-over-seeds columns mixed seeds | corrected; the report prints per seed |
| D3h-10 | cosmetic (tests) | slacks of 10³–10¹⁷; test 3 compared sorted spectra; test 4 ignored the Newton residuals, and its ε = 4 point violated A ≻ 0 | the ε → 0 limit tested at its rate (at an ε where the error is near the tolerance) and the uniqueness of the OT plan checked (exchange gaps > 1); the separation blocks compared entrywise, with the invariance of the subspace; the threshold test by continuation in ε at points where A ≻ 0 (checked), every residual below 10⁻⁹ |
| D3h-11 | cosmetic (`LITERATURE_TRANSPORT.md`) | the D‴.21 commit had put the D‴.20 Verdict line under the D‴.21 bullet | moved back and amended per D3h-6 |

**Lesson.** A bisection on a quantity computed by an iterative solver locates the solver's failure boundary as readily as a sign change:
every solve inside a root-finder must be continued from a converged neighbour and checked. And a stopping-time snapshot of a slow saddle
(e-folding times of 700–2900 slow units at margins of −0.0004 to −0.0015) is not an end state: before a run's end weights are read as a
threshold, the margins at those weights must be computed and the run continued until they are positive.

## Round 16 (adversarial review of Theorem D‴.21, Remark D‴.22 and their numerics, 2026-10-19)

**Design.** One referee on Theorem D‴.21 (statement, proof, Status, numerics), Remark D‴.22, the §0 and §7 rows, the §8 bullet, Remark
D‴.19's last sentence, the RESULTS section, `src/core/ot_basin.py`, `exp_ot_basin.py` with its raw data, `tests/test_ot_basin.py`, [S47],
[S48] and the audit bullet (`review/REVIEW_round16_D3i.md`). An independent implementation with analytic derivatives (no import of `src/`):
(E) and (P) with genuine segment integrals at 50 configurations (optimal and non-optimal assignments, both priors); the chain of (ii) with true
segment suprema at 58 configurations; 34 single systems re-integrated with the experiment's integrator, bitwise equal to the stored records;
the a-priori constants, (H*) over the whole parameter space, the sublevel set, the nearest-atom end states and the admissibility of their
counts.

**Verdicts.** (i), (ii), (v) hold; (iii) holds for vertex-selection solutions (their existence assumed); (iv) holds for solutions that stay in
B, which must be forward-invariant; (vi) holds as stated with two overstatements; Remark D‴.22 overstates. The numbers hold; the reading of
the nearest-atom runs fails in 10 of 32 runs. The tests pass and are not vacuous, but three are weaker than their names.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3i-1 | **major** (numerics, §0, §7, RESULTS, md) | "the floor is a spurious minimum … the nearest-atom runs stop at the energy of the replicated configuration of their sampled weights" fails in 10 of 32 runs: at m ≤ 60, 5 of 8 sampled count vectors per prior are inadmissible (class-0 margins −0.0004 to −0.0061), their replicated configuration is a saddle, class 0 splits and E ends 0.5–16.5 % below it, still falling; the one-particle runs quote E at s = 300 (0.91 of the spurious minimum) | `exp_ot_basin.py` prints, per run, the admissibility of the sampled counts (the separation margins and the whitened weighted K-atom Hessian), E at s = 300, 450, 600, the end state's class spreads and its distance from Θ*(w); the reading split: at the 22 admissible count vectors E sits 0.02–0.24 % below the spurious minimum and rises toward it, at the 10 inadmissible ones the class splits; the one-particle nearest-atom runs continued to s = 3000, where E equals the spurious minimum's 6.3126·10⁻³ and 6.1328·10⁻⁴ |
| D3i-2 | minor (§0, §8, Remark D‴.19, numerics) | "the basin contains configurations with a fixed fraction of the particles in other cells" is vacuous at the proven E₀ ({E ≤ E₀} holds a fraction ≤ 6·10⁻⁶, so one misplaced particle needs m ≥ 1.7·10⁵); the 1/30 example is covered only a posteriori; "larger" is not an inclusion | the script computes the cell distance and the largest misplaced fraction p ≤ 2E₀/(κd²) (5.8·10⁻⁶ matching, 4.4·10⁻⁷ shifted on the whole space; m ≥ 1.7·10⁵ and 2.3·10⁶); the rows say "tiny on the instance" and "a posteriori wherever the margin is positive"; Remark D‴.19: "of a different shape (not nested)" |
| D3i-3 | minor ((iv), proof, numerics, docstring) | (iv) needs B forward-invariant and the numerics' box is not; the box supremum was 1 % low (random points only); (H*) actually holds on the whole space at the shifted prior | (iv): "a convex set B, forward-invariant, … for instance the whole parameter space"; the script searches the whole space (a grid over the saturation ranges, corner scans, Nelder–Mead): sup = 6.60·10⁻³ < σ₀ at a saturation corner, κ = 4.7·10⁻⁴, E₀ = 9.1·10⁻¹⁰, and the box with L-BFGS-B polish (4.39·10⁻³); the docstring names the box correctly |
| D3i-4 | minor ((iv), (vi), §0, §8) | "under (H*) no minimum other than μ* is left" is false for local minima (a strict local minimum of J^{OT}_3 with E = 0.107 at the matching prior, the referee's computation); "for every count vector N" holds only near q (at m = 15, 35 of 90 vectors are saddles, one has no stationary point found) | "the unique global minimiser"; "the only critical point where α > 0 (in particular in {E ≤ E₀}); local minima at E > 0 need κ^Θ ≥ 2σ₀ − κ*"; "for every N in the domain of w ↦ Θ*(w) whose atoms lie in their cells … a saddle when a separation margin is negative" |
| D3i-5 | minor (Remark D‴.22, §8) | "propagation of chaos uniform in time with no floor" is overstated: the bound shows that the long-time limits coincide, not smallness in m at a fixed time | "for every commensurate m the particle system and the mean-field flow have the same limit — no floor; this is not uniform-in-time propagation of chaos"; [S50] (Chen–Ren–Wang 2022; Suzuki–Nitanda–Wu 2023) as the comparison |
| D3i-6 | minor (Status, §0 novelty, [S48], audit bullet) | the plan segments are AGS's generalised geodesic based at ν; PL inequalities through optimal maps have precedents (Chewi et al. 2020; Mérigot–Santambrogio–Sarrazin 2021); the reason given for [S48] not applying is wrong (it covers λ-convex functionals); the metric-slope qualifier misplaced; the double verdict | Status and the audit bullet rewritten (the generalised geodesic, the semiconvexity with a negative constant as the reason, the specific content: the loss linearised at the stationary residual and the counts); [S48] extended; the slope W(·, ν) at every μ, a velocity field realising it where the plan is a map; the stray verdict moved back (round 15) |
| D3i-7 | cosmetic (numerics, §7, RESULTS, md, script) | "1.58–1.72" mixes two columns (ratio/λ is 1.61–1.72); W² 5.2 is 5.25; α at p = 1/30 is +0.0023 with true suprema; the factor 0.55 is 0.80 (curvature) × 0.69 (slack of (ii)); "identically at m = 30 and 300" holds by construction; two docstring/comment mismatches; ∫λ₊ over 5-unit steps | corrected; the measure-level table uses 129 points per segment (α = +0.00229); the report prints the local PL constant 2min{λ_sep, λ_coll} = 0.01135 and the split 0.683 × 0.802; the one-particle OT run evaluates λ at every step (bounds 9.9·10⁻⁵ and 3.7·10⁻⁶); "the same measures"; the comments fixed |
| D3i-8 | cosmetic (statement, proof) | the lower bound on E in (ii) needs no α > 0; "every optimal plan" should be "every vertex optimal plan"; (E) holds for every coupling, only (P) needs an assignment, and the counts' defects are explicit; the existence of vertex-selection solutions is assumed | (ii) restated with the two-sided bound for every assignment; "vertex optimal plan"; (i) states the coupling form, the two defects at counts N ≠ mq and the failure at incommensurate m; (iii) and the Status: existence assumed, a two-particle tie never attracting (the sides' dh/ds = ℓ ∓ σ₀Δ²), multi-way ties not analysed |
| D3i-9 | cosmetic (tests) | test 1 trivial at the matching prior and never exercises the counts; test 3 does not assert its margin; test 4 does not check the K-atom Hessian | test 1 at the shifted prior with counts N ≠ mq against the two defect formulas (to 10⁻¹² relative; the (P) defect is 7 times Q); test 3 asserts α > 0 at the start; test 4 checks the whitened weighted K-atom Hessian (smallest eigenvalue 0.0073); `BasinDiagnostics.terms(check_counts=False)` and `counts_mismatch` |

**Lesson.** The admissibility of the realised weights decides whether a replicated configuration is a minimum or a saddle, and round 14's
D3g-1 had already shown it: every reading of runs at sampled weights must print the margins at those weights. And an a-priori basin must be
measured in the units it is used in: an energy level that holds no misplaced particle below m = 10⁵ does not support "a fixed fraction in
other cells" — the a-posteriori form does, along the run.

## Round 17 (adversarial review of Theorem D‴.23 and its numerics, 2026-10-19)

**Design.** One referee on Theorem D‴.23 (statement (i)–(vi), proof, Status, numerics), the §0 and §7 rows, the §8 bullets of branches 7
and 8, the RESULTS section, `src/core/adaptive_prior.py`, `exp_adaptive_prior.py` with its raw data, `tests/test_adaptive_prior.py`, [S49]
and the audit bullet, README's sentence (`review/REVIEW_round17_D3j.md`). Every number of the md, the paragraph, the rows and RESULTS
recomputed from the raw records independently of `report()`; single systems re-integrated bitwise; the joint Jacobian at m = 12 and 30; a
non-replicated stationary point and an unbounded solution constructed; the plateau's Hessian by Richardson extrapolation and the
unregularised gradient flow from it; the radius of (v) bounded on the instance.

**Verdicts.** (i) holds for replicated configurations, its converse up to relabelling, with boundedness a hypothesis; (ii) holds; (iii)
holds to first order; (iv) holds as stated (formal); (v) holds as a sketch, "proved" overstated; (vi) holds in substance. Every number
reproduces; three readings fail or overreach (the uniformity in m of (v), "every mode", "to 4.8%"). The tests pass and are not vacuous; one
name says more than its tolerance checks.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3j-1 | minor ((i) heading and body, proof, §0, §8 branch 8, docstring) | only *replicated* stationary points are characterised: non-replicated stationary points with a non-zero leak force exist (a four-bump target, m = 6, one class split, a saddle); the converse holds up to the plan's relabelling (a permuted prior is stationary at the exact representation); "the stationary set does not depend on the initial prior" is a tautology | (i) "whose replicated stationary points carry no leak bias", the stationary equations G⟨r∇Φ̃(ϑ_i)⟩ = σ₀(ϑ_i − T̂_π(i)) with z = T̂, the classes labelled by the plan (z_σ(k) = ϑ_k), the non-replicated points with the K₀-cluster bias, "the stationary equations do not involve the prior's initial position; which point is reached does"; (iv): K₀ remains a design input; new Part G of the experiment (the saddle, the flow from the exact split fit, which collapses the class and approaches a boundary critical point with the merged width at s_max; the permuted prior) and a test |
| D3j-2 | minor ((i), proof, §8 branch 7) | boundedness is no longer automatic: the leak toward an adaptive prior is an internal force (the conservation law d/ds[(1/m)Σ_{i∈k}ϑ_i + (σ₀q_k/η)z_k] = (1/m)Σ_{i∈k}G⟨r∇Φ̃⟩), and beyond the saturation the amplitude coordinate grows like 2 ln s | (i) states boundedness as a hypothesis with the conservation law, the proof derives it; Part G(b) (LSODA to s = 10⁶: the fixed prior converges, η = σ₀ and ∞ run off, increments per decade 1.99 ln 10, a → a_max like 1/s) and a test; §8 branch 7: the adaptive prior gives up the leak's a-priori bound, the robustness half of the trade-off |
| D3j-3 | minor ((v), proof, Status, §0) | (v) is a two-sentence pointer: the radius bookkeeping is missing (the prior's offsets reach ρ(η/(σ₀q_k))^{1/2}), the limits η → 0, ∞ untreated, and the runs start 94–6·10⁴ times outside the ball | (v) states ρ(η) = min{λ₀(η)/(4C_N), margin/(C[1 + q_min^{−1/2} + (η/(σ₀q_min))^{1/2}])}, the η = ∞ norm (ζ = ū) and the degeneration as η → ∞ (singular perturbation) and η → 0; Status "sketched"; the script computes a sampled lower bound C_N ≥ 0.3406 (the full collective Jacobian), the radius bound per η and the starts' norms (the factors 94–59507) |
| D3j-4 | minor (numerics, §0, §7, RESULTS, md) | the settled gate's agreement with the K₀-atom run does not test (v): the systems open replicated, so it is the invariance of the replicated configurations | the report prints the opening state (rms separations 1.2·10⁻⁴–7.9·10⁻⁴, class means within 6.7·10⁻⁵ of Θ*(q), 32 runs) and the reading; new Part F tests (v) from separated starts (rms 0.02, 0.1; m = 15, 510; η = σ₀, ∞; 1750 slow units): departures +0.07–0.08%, +1.8–2.1% (σ₀) and +0.44–0.59%, +10.9–14.5% (∞), growing by 24.6–25.0 for a factor 5, alike at both m |
| D3j-5 | minor (Status, §0 novelty, [S49], audit bullet) | the per-mode formulas are the classical quadratic analyses (Tikhonov filter factors, iterated Tikhonov / Hanke–Groetsch, EASGD's worker–centre map) | Status and the audit bullet cite them as the source of the per-mode formulas (the Moreau-envelope flow as the slow-prior limit); "new" restricted to the transfer to the mass-constrained particle system uniformly in m, the drag with its gate and the limits of (i); Hanke–Groetsch added to [S49], EASGD's quadratic analysis described |
| D3j-6 | cosmetic ((ii), proof, §0, §8, (i)) | the strict upper bound fails when h_jη = 0; the lower bound is strict for h_j, η > 0; the monotonicity of λ⁻ in h is used but not stated; (i) needs η > 0 | (ii) "for h_j, η > 0: ηh_j/(h_j + σ₀ + η) < λ_j^− < min{h_j, η}", the limits h_jη = 0 and η = ∞, dλ⁻/dh = (η − λ⁻)/(λ⁺ − λ⁻) > 0 with its proof; (i) "for η > 0" |
| D3j-7 | cosmetic (numerics, md, §7, RESULTS) | `report()` computed slopes and ratios for five of nine modes; the §7 row's residue lacks "(η = σ₀)" | all nine modes printed; the report computes the ranges over every mode above 10⁻⁶ (13 of 36: log-slopes 0.981–1.003, ratios 0.820–1.115, the ratios' change from s = 1000 at most 1.04% for η ≥ σ₀, 5.96% at σ₀/10); "(η = σ₀)" added |
| D3j-8 | cosmetic (numerics, md, docstring) | "describes it to 4.8%" is a norm dominated by the strong modes (modes 1 and 3 at 0.85 and 0.73 of the prediction); y⁰ for x⁰ | the md prints actual/predicted per mode and the error against the shift (4.80, 2.34, 1.16, 0.58% in norm; modes 1 and 3 at 0.977 and 0.969 at the smallest shift); the paragraph quotes both; x⁰ throughout |
| D3j-9 | cosmetic ((vi), proof, Status, §0) | 𝒥 is non-increasing under every schedule of η, gated or not; at η = ∞ the opening is a jump | (i), (vi) and the proof: 𝒥 does not involve η; the opening at η = ∞ sets Z to the minimiser over Z for the current plan and lowers 𝒥; the gate is there against the drag; Status "the monotonicity of 𝒥" |
| D3j-10 | cosmetic (test 7) | the name says the gate opens at the fixed-prior stationary point; at the test's tolerance it opens with rms separation 0.042 | renamed `test_the_settled_gate_opens_when_the_fixed_prior_drift_first_falls_below_its_tolerance_after_the_transit` |
| D3j-11 | cosmetic ((vi), §0, README) | "must" too strong: the settled gate waits ≈ 1200 slow units, and at η = σ₀, m = 510 the transit gate's fit error is lower at s = 3000 (1.547·10⁻⁴ against 1.638·10⁻⁴); the README's remedy is the transit gate's | (vi) "when the parameters, not only the fit, are wanted at a given horizon … should"; the wait's cost quoted in the numerics and RESULTS; §0 "should"; README "after the transit and the fixed-prior relaxation" |
| D3j-12 | cosmetic (md Part E, numerics) | "the loss descends from there to the exact representation" rested on BFGS | Part E integrates the unregularised K₀-atom gradient flow (LSODA): above 0.9 of the plateau's loss for 4.18·10⁴ more slow units, halved by 4.87·10⁴, below 10⁻¹² by 8.00·10⁴, at the exact representation; the time scale 1/1.84·10⁻⁵ = 5.43·10⁴ |

**Lesson.** A characterisation proved for one class of configurations must say so in its heading: the stationary points of a joint flow
include the ones the derivation did not look at, and a regulariser that is an internal force (the leak toward a learned centre) gives up
the confinement that made bounded solutions automatic. And an agreement between two runs tests a theorem only if the runs could have
disagreed without it: systems that open replicated agree with the K₀-atom run by invariance, whatever (v) says — separations must be put
in at the start.

## Round 18 (adversarial review of Proposition D‴.24, Remark D‴.25 and their numerics, and of Theorem D‴.23 as revised in round 17, 2026-10-19)

**Design.** Two referees in parallel on a snapshot of d2092a5 (`git archive`), neither modifying the repository. D3k: Proposition D‴.24
((i)–(vi), proof, Status, numerics), Remark D‴.25, the §0 and §7 rows, the §8 bullets, the pointer in Remark D‴.10, RESULTS,
`src/core/floor_bias.py`, the nearest-atom option of `src/core/adaptive_prior.py`, `exp_adaptive_floor.py` with its raw data,
`tests/test_adaptive_floor.py`, [S51] and the audit bullet, README (`review/REVIEW_round18_D3k.md`): an independent implementation of the
instance with analytic derivatives; the D‴.5 margins at all 16 sampled weight vectors; D and W² from the stored finals; exact multinomial
probabilities, Hoeffding's bound and permuted representations; a Monte Carlo of 600 nested samples with the exact stationary points; the
boundary case integrated to s = 10⁶; crossings and an empty cell. D3l: Theorem D‴.23 after round 17 (`review/REVIEW_round18_D3l.md`): Part
G(a) re-run bitwise and continued to 10¹⁰ slow units, its end state's face minimiser; a counterexample with unequal masses; five further
separation draws for Part F; the start factors of Parts D and F; per-system gate comparisons; event location for Part E; a variant of
test 9 at η = σ₀/10.

**Verdicts (D3k).** (i) holds — F is exact at q, the bias formula first order, (B) is D‴.9's law on the eigenvectors (identical in three
forms) — but it holds only on D‴.5-admissible samples, and its value is the asymptotic one; (ii) holds inside the cells, and 𝒥_∘ stays
non-increasing across crossings; (iii) holds for the stationary points, the parameter-level claim (D larger under the adaptive prior)
confirmed dynamically, the probability statement overstated on the instance; (iv) holds with its constant (0.40/s and 0.80/s); (v) and
Remark D‴.25 hold in part (the parameter-level column mixes D and W²); (vi) holds. Every number of the md, the rows and RESULTS reproduces;
`--from-raw` is byte-identical. **Verdicts (D3l).** Round 17's fixes hold except Part G(a)'s end state, which is a point at infinity, not a
boundary critical point; the rest are wording.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3k-1 | major (D‴.24 numerics B–D, Status, test 4, RESULTS, the rows; "admissible" throughout) | 5 of the 12 fixed-prior m-systems sit at weights where Θ*(w) is a D‴.5 saddle ((6, 5, 4)/15, (6, 3, 6)/15, (24, 15, 21)/60, (8, 4, 3)/15, (25, 17, 18)/60; class 1, w₁ ≥ 0.4) and split, with floors 0.81–1.01 of the replicated point's; the W² example compared with a split configuration; Part D's m = 15 fixed run and test 4 use a saddle sample; "admissible" used for two sets; the fixed prior's saddles are as likely as the adaptive prior's non-representability (0.37, 0.26, 0.11) | the experiment's new checks stage prints per sample the D‴.5 margin of Θ*(w) (matching and misplaced prior), the end class diameters and sub-clusters, the exact bias at the sampled weights against the observed floor (ratio 1.0000 at the seven stable samples) and the delta method (0.86–1.72); "representable (w ∈ W)" and "D‴.5-admissible (w ∈ W_adm)" distinguished in D‴.24, Remark D‴.25, the rows, RESULTS and §8; (B) stated for D‴.5-admissible samples; (iii): at an exact representation every splitting margin is σ₀ (R = 0), so the adaptive prior also removes the splitting instability; Part D's saddle orbit noted; test 4 at the stable sample (5, 4, 6)/15 with its margins, and the margins σ₀ at the exact representation of the saddle sample (6, 5, 4)/15 |
| D3k-2 | major (Remark D‴.25's column and Reading, (v), numerics, §0, RESULTS) | the column "parameter-level floor" mixes the coupled distance D (nearest-atom rows) and W² (OT rows); under D the OT designs are not 0 (the transit's mismatch, the largest at m = 240), under W² the nearest-atom rows are the cross-class transport ½‖δw‖₁ × inter-cluster distances, common to both priors | the column split into D (D‴.5) and W² (measure level); the Reading rewritten (three faces of the floor: function, coupling, measure; the mass constraint removes the function- and measure-level floors and replaces the gauge's displacement by the transit's mismatch); the checks print D (observed against the replicated stationary points' 0.944–1.000, adaptive/fixed 1.72–4.81) and W² split into cross- and within-class parts (cross 0.67–0.99; the priors' W² within 1.00–1.32) for Parts B and C |
| D3k-3 | minor (numerics, (i), md) | "the floor's expectation" is the law's value; the law is far from the expectation at m ≤ 60; "the floor is at most" bounds the law | "the law's value", "the law is at most"; Part F of the experiment: a Monte Carlo of 600 nested samples (the exact floor 2.78 ± 0.20, 1.53 ± 0.11, 1.11 ± 0.07 times the law; D/law 0.68–0.95 fixed, 0.53–1.26 adaptive); Status: the law is asymptotic |
| D3k-4 | minor ((iii), (vi), §0, §8 branch 8, audit) | Hoeffding's bound vacuous at m = 15; the expected limiting floor 1.8 times the fixed law at m = 15; "no exact representation" false under permutations; the parameters run off outside W; "admissible m" | Part E of the experiment: exact P(outside W) 0.361, 0.074, 0.0021, P(boundary) 0.118, 0.053, 0.0012, Hoeffding 1.87, 0.60, 0.036, the permutation share 31%, 87%, 100%, the fixed prior's saddles, the adaptive prior's expected limiting floor ≥ 7.3·10⁻³, 2.2·10⁻⁴, 6.6·10⁻⁷ (constrained minima); (iii) "on the branch of ϑ^e(q)", "if every representable sample converges … exponentially small … not yet in force at m ≤ 60", "conditionally on representability (outside W the parameters run off)"; "at every representable sample" |
| D3k-5 | minor ((i), §0, §8 branch 7, numerics, test 3) | the concentration in modes 2–6 is mostly the gauge's (97.2% unweighted); the weight is a broad band-pass; test 3's share assertion passes without the weight | "a band-pass centred at h = σ₀, at least half its peak for h/σ₀ ∈ [3 − 2√2, 3 + 2√2]; the floor sits in the modes of the gauge's displacement, re-weighted toward h ≈ σ₀"; the report prints weighted and unweighted shares (the band [1/5, 5]: 42%, 99%, 90% against 21%, 97%, 78%); test 3 compares them at 0.05 and 0.5σ_crit |
| D3k-6 | minor (Status, §0, [S51], audit) | Chen–Rotskoff–Bruna–Vanden-Eijnden's vanishing deviation at an interpolating limit is the mechanism of (iii); the weight's bound is ridge regression's AM–GM bound; F is the chain rule on D‴.23(iii) | [S51] annotated; Status and the audit bullet cite all three as classical; "new" restricted to the identification of (i) with its law and the bound's identity with the adaptive prior's coupled-distance law (σ₀/4 times it, found by the referee), the split (iii), and a learned prior restoring the interpolation property under regularisation |
| D3k-7 | cosmetic (Status, §0, numerics) | the exactness was placed on a trivial identity; F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ is exact at q; the quoted discrepancies are finite-difference errors; (B) is formal | "F exact at q, the bias formula first order, the law (B) formal (the delta method, as D‴.9)"; "to the finite-difference error" |
| D3k-8 | cosmetic ((ii), proof) | "masses learned by the same energy" is the minimum over the masses; a tautology on crossings; the empty cell's degenerate metric; the crossings | "minimised out (semi-relaxed transport, the infinitely fast limit of a mass dynamics)"; across crossings the vertex-selection flow with 𝒥_∘ non-increasing (proof by the Dini derivative); the empty cell's atom frozen; a new test (a crossing at η = σ₀ and ∞ with 𝒥_∘ decreasing at every refresh; an empty cell staying empty, its atom frozen) |
| D3k-9 | cosmetic ((v), the D‴.10 pointer, Part C) | the trade-off's price understated (D larger, the a-priori bound, the collective rate 49× and 24× below the fixed prior's); Part C's ranking depends on the horizon | "at the price of the collective rate, of the leak's a-priori bound and of a larger parameter-level floor"; the report prints the fit errors at s = 1250 and 2000 and the gate's opening (1211–1304) |
| D3k-10 | cosmetic (rows, Remark D‴.25, (vi), `floor_bias.py`) | "10⁻⁷" for 1.4·10⁻⁷; the e_rms wording; "the stationary replicated points are the exact representations"; the entropic row's law; "none (counts conserved)" is D‴.9(i)'s hypothesis; min{h_j(w), η}; (vi) needs f̂ = f on supp ρ; `bias_floor_law(directions=…)` correct only for orthonormal directions | all reworded; `bias_floor_law` uses Vd⁺C̃Vd⁺ᵀ (the default, the explicit e_j − e_K and an orthonormal basis now give 0.06166) |
| D3k-11 | cosmetic (tests) | test 2 tests random masses only; test 3's tolerances loose and its share check weight-independent; test 4 at a saddle; no test beyond the cells | test 2 perturbs the cell masses by ±1/m (strictly larger energies); test 3's tolerances 2·10⁻⁶ and 10⁻⁶ and the weighted/unweighted comparison; test 4 as under D3k-1; the crossing and empty-cell test |
| D3l-1 | minor (md Part G(a), `boundary_job`, numerics, RESULTS, round 17's D3j-1 row) | the end state of the flow from the exact split fit is a point at infinity: the width coordinate runs off like ln s, the width tends to s_max like 1/s, and the physical image approaches the loss's minimiser on the face s = s_max, a Karush–Kuhn–Tucker point (∂L/∂s = −1.3·10⁻⁶), not a critical point; the "refinement" by least squares stopped on a slope | `boundary_job` continues the collapsed flow to 10¹⁰ slow units (κ = 2.16 … 9.75, increments per decade → 1) and computes the face minimiser (fit error 1.943283·10⁻², ∂L/∂s = −1.326·10⁻⁶); md, numerics and RESULTS: an unbounded solution from an exact fit; (i) cites Part G(a) and (b). *Correction to round 17:* D3j-1's "approaches a boundary critical point with the merged width at s_max" is wrong in the same way |
| D3l-2 | cosmetic ((i)) | a permuted prior is stationary only among atoms of equal mass (a counterexample with q = (1/2, 1/3, 1/6)) | "a prior permuted among atoms of equal mass … a permutation between atoms of different masses is not, since W² = 0 needs q_σ(k) = q_k" |
| D3l-3 | cosmetic (Part F, numerics, rows, RESULTS) | Part F is not a test of (v) (its starts lie 133–3344 times outside the ball), and its m-comparison rested on one draw per m | "the structure behind (v), beyond its ball"; five draws per (m, η, rms), means, ranges and standard errors printed (the means at m = 15 and 510 within 0.9–1.7 combined standard errors, the spread 3.3–6.0 times larger at m = 15), with the reason the expectation is m-independent at second order |
| D3l-4 | cosmetic (md Part A, numerics, RESULTS) | the start factors were Part C's; Parts D and F start farther out | the report prints them (Part D's warm clouds 2.2·10⁴–3.6·10⁵, Part F 131–3344 times outside) |
| D3l-5 | cosmetic (numerics, README) | the settled gate is lower in the seed means, not system by system (the transit gate lower in all three in 6 of 32); README made the gate a condition | "in the seed means; … 6 of the 32"; README "is best started after the transit and the fixed-prior relaxation" |
| D3l-6 | cosmetic (Status, audit bullet, [S49]) | the Tikhonov filter factor is h/(h + σ); σ/(h + σ) is its complement; the limits of the bias removal are largely the precedents' own properties | "the complement 1 − h_j/(h_j + σ₀) of the Tikhonov filter factor (the residual factor)" in all three places; the limits attributed (the proximal point method without a minimiser, [S49]; soft weight-sharing's quantisation error), specific only in their particle form and the run-off from an exact fit |
| D3l-7 | cosmetic (md Part E, numerics, rows, RESULTS) | the crossing times were grid points | located as events: 4.156·10⁴, 4.780·10⁴, 7.971·10⁴ |
| D3l-8 | cosmetic (test 9) | the test asserted finite-time numbers calibrated to η = σ₀, which a bounded solution could pass | the test asserts the mechanism: the force ⟨rφ₁⟩ tends to the positive Karush–Kuhn–Tucker multiplier of the loss's minimum on the face a₁ = a_max (to 0.3% at s = 10⁵), and the run-off law (a_max − a)s → a_max²(σ₀ + η)/(4ηG_αF_∞) holds (to 1.4%), derived from the conservation law with α − z_α slaved |

**Lesson.** The rounds-14/16 lesson again, in a new place: every reading of runs at sampled weights must print the stability margins at
those weights — and a word ("admissible") that names two sets must be split before it is used for both. And a distance is a choice: the
coupled distance of shared particles and the measure-level W² can rank the same designs in opposite orders, so a design table needs a
column for each. A law obtained by the delta method is a statement about m → ∞; at the sizes of the numerics it must be checked against
the exact quantity it linearises (here the exact bias at the sampled weights, which reproduces the floor to four digits, against the
law's factor 2.8 at m = 15). Finally, a point reached "on the boundary" in unconstrained coordinates is at infinity: the flow can only
approach a Karush–Kuhn–Tucker point of the closed problem, and the solution is unbounded.

## Round 19 (adversarial review of Corollary D‴.26, the planar check of Proposition D‴.24, Proposition D‴.27, and the application of round 18, 2026-10-19)

**Design.** Two referees in parallel on a snapshot of 442cd6c (`git archive`), neither modifying the repository. D3m: Corollary D‴.26
((a)–(c), the stability consequence, proof, Status, numerics), its rows, the §8 bullet, RESULTS and the audit bullet,
`src/core/adaptive_floquet.py`, `tests/test_adaptive_floquet.py`, `exp_adaptive_floquet.py`; the planar check (`src/core/planar_loop.py`,
`tests/test_planar_loop.py`, `exp_adaptive_floor_2d.py` with its checks stage) (`review/REVIEW_round19_D3m.md`): an independent
implementation of the instance with the variational equation along the lift (k_e = 12.5–400), the finite differences' resolution at five
steps and tolerances, an asymmetric reference, weights approaching ∂W, the OT plan under perturbations, the planar loop against
independent code, three non-representable systems continued to s = 6·10⁴. D3n: Proposition D‴.27 and the round-18 fixes
(`review/REVIEW_round19_D3n.md`): twelve stress runs of the bound, an unbounded convex set, the projected midpoint rule at six steps against
an event-driven reference, the KKT condition, interior critical points, the four-bump target with the box to 10¹⁰, the closed-form
constant against a search; for round 18 an independent model (margins, exact bias, D, W², the laws), exact probabilities and a Monte Carlo
with its own seeds.

**Verdicts (D3m).** (a) and (b) hold exactly, for finite and infinite η (the separations at σ₀/k_e to 10⁻¹² by the variational equation;
the md's 10⁻⁶ is the differences' resolution); (c) holds, in slow time O(k_e⁻²) to k_e = 400, but the numerics could not test its
mechanism (the first-order correction vanishes identically on the reference) and it needs ∂_xf̂ = ∂_xf on the orbit; the stability
consequence, the OT extension and the fixed-prior claim hold; the planar check holds in every number. **Verdicts (D3n).** (i) holds as a
distance bound, but boundedness needs a compact 𝒵; (ii) holds, its convergence provable at η = ∞; (iii) holds trivially; (iv)'s
identification is exact, "stops on the boundary" and "between" are not theorems; the projected midpoint rule can stall the prior inside a
face (Part H's numbers right only because dt = 0.25 reached the face); every round-18 fix is in place, three readings in the rewritten
text wrong.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3m-1 | minor (D‴.26 numerics, md, rows, RESULTS) | on x_m = 2 sin t the first-order averaging correction vanishes identically (⟨ĴᵀĴ̇⟩ = 0 in one dimension, the commutator term by t → π − t), so the observed second order could not test (c)'s mechanism; two values of k_e; the remainder unnamed | `CollapsedVariational` in `adaptive_floquet.py` (the variational equation on the collapsed subspace in whitened coordinates, its monodromy without differences, Liouville's formula, the averaged generator to second order) and Part C of the experiment (k_e = 12.5–400, the asymmetric reference 2 sin t + 0.8 sin 2t, where ‖Ā₁‖ = 5.6·10⁻³, skew): the first-order deviations fall by → 4 per doubling and the second-order residuals by → 16 on both references — the expansion is even in 1/k_e, proved by the time-reversal symmetry M(t)ᵀ = DM(t)D⁻¹ of the adjoint pair; (c) restated as (1/k_e)[λ_j^± + c_j^±/k_e² + O(k_e⁻⁴)] with c_j explicit (c₁/λ₁ = −100.6, the error filter's attenuation −102.1 of it; every c_j < 0 here) |
| D3m-2 | minor ((c), proof, D‴.24(vi)'s gloss) | the linearised error equation drops ∂_x(f − f̂)δe, which needs ∂_xf̂ = ∂_xf on the orbit, not supplied by "f̂ = f on the orbit"; without it the correction is first order | (c) assumes ∂_xf̂ = ∂_xf on the orbit (the exact representation supplies it, an identity of functions); the proof carries a(t) := ∂_x(f − f̂)(x_m(t)) and shows its first-order, self-adjoint term −⟨J̄ᵀaJ̄⟩ and why it breaks the parity (a rate even under time reversal); the gloss of (vi) amended |
| D3m-3 | minor (§0 last column, audit bullet) | (a)–(b) are D‴.14(ii) and D‴.15 at zero residual; a leak compatible with exact tracking is classical (e-modification; a centre that follows the parameters); the skew first correction is the Magnus commutator; "[S49]'s trade-off" is Ioannou–Kokotović's | Status, the §0 row and the audit bullet cite D‴.14(ii)–D‴.15, the e-modification, the proximal point, the Magnus expansion and the averaging analysis of adaptive loops ([S52], new: Narendra–Annaswamy 1987, Magnus 1954, Blanes et al. 2009, Kosut–Anderson–Mareels 1987; bibliographic details checked by a web search); "new" restricted to the zero-residual lift at every representable sample with its exactly m-independent spectrum and the even expansion with its explicit coefficient; [3], [4] for the trade-off |
| D3m-4 | cosmetic (proof) | M(t), N(t) lack W⁻¹; the citations of D‴.14 and D‴.15 loose | the proof rewritten in the whitened weighted coordinates, where the learning law is exactly the adjoint of the error injection; (b) cites D‴.14(i) with the prior's equation, (c) the slaving (the Chang transformation, [S9]) and second-order averaging |
| D3m-5 | cosmetic (§8, RESULTS, rows, docstring) | §8 said O(k_e⁻²) where (c) says O(k_e⁻³); 14%/3.9% against 13.6%/3.8%; "linearly stable / nonlinear basin not proved"; the remainder's sign; the error channel unmeasured; a stale docstring | §8 and the rows in the new form; one set of percentages (13.6%, 3.84% at η = σ₀); "exponentially stable (locally, for each m)", the basin uniform in m open; the sign (every exponent below its first-order prediction) printed; the error channel by Liouville's formula (24.98084 against k_e − tr H_L/k_e = 24.98037); the docstring points to `monodromy` |
| D3m-6 | cosmetic (`tests/test_adaptive_floquet.py`, `tests/test_planar_loop.py`) | the Floquet test runs at η = ∞ and w = q only (the prior in the state never executed; a fixed-prior mutant passes three of four assertions); "for every m" is m = 6; the comment's "O(1/k_e)"; the planar loop's finite-η branch untested | three new tests at (2, 1, 1)/4 with the prior in the state: the real loop's Jacobian at the lift splitting into the separations (−σ₀/k_e) and the collapsed variational equation; the m = 4 lift's monodromy (3 separations, 18 collective exponents — the count fails a fixed prior —, the second-order prediction); on the asymmetric reference Ā₁ skew, the residual O(k_e⁻⁴), c₁ measured, Liouville; the old test renamed "… at m = 6", its comment O(k_e⁻²); a planar test of the finite-η branch against DOP853 (errors falling by 4 per halving) |
| D3m-7 | cosmetic (`exp_adaptive_floquet.py`, `tests/test_planar_loop.py`, RESULTS) | "admissible" for representable weights after round 18 | "representable" in the four places |
| D3m-8 | cosmetic (the planar check, the D‴.24 §7 clause, RESULTS, md Part B) | Part B's η = ∞ is a lagged prior (reset every slow unit); Part C's adaptive prior is η = ∞ only; the non-representable systems have a floor; "four kernels" | the md, the paragraph and RESULTS state the reset and "(η = ∞)"; "still falling, toward the constrained minimum at their sampled weights — a floor of their own — not toward zero"; "four bumps"; "function-level distance" |
| D3n-1 | minor (D‴.27, rows, §8, branch 7, D‴.23(i), audit) | "closed convex" 𝒵: the distance bound holds but boundedness needs a bounded 𝒵 (an unbounded slab lets the amplitude run off, 22.2 at s = 10⁶); Z(0) ∈ 𝒵 implicit | "a compact convex set, Z(0) ∈ 𝒵^{K₀}" everywhere; the sharper bound B/σ₀ + (d(0) − B/σ₀)₊e^{−σ₀s}; the unbounded case remarked |
| D3n-2 | minor (`AdaptivePriorLoop.run`, Part H, test 10) | the projected midpoint rule projects the velocity onto the tangent cone at the midpoint: when the midpoint lies past a face the step returns the prior unchanged, a spurious fixed point 1.2–2.4·10⁻⁴ inside 𝒵 (dt = 0.2, 0.125, 0.1, 0.0625; drift 5.9·10⁻³) | `run` projects positions, not velocities (Z_mid = Π(Z + (dt/2)v), Z_new = Π(Z + dt·v(Z_mid))), unchanged without a box; all six steps now reach the face and agree to ten digits; Part H at dt = 0.25 and 0.125 (end states within 3.4·10⁻¹², drift asserted ≤ 10⁻¹²; the md's rows unchanged); test 10 at dt = 0.125 asserts the atom on the face (the old step fails it, 1.17·10⁻⁴ inside) |
| D3n-3 | minor ((iv), rows, §8, numerics, md) | "stops on the boundary" does not follow (Part G(a)'s saddle has an interior prior); "between" is a property of one run (an m = 15 system ends at 0.0534, below the unprojected limit) | (iv) conditional: "if the limit has an atom on ∂𝒵 — necessarily when no stationary point has its prior inside"; "between" an observation of the K₀-atom run, with Part H(c): the amplitude bound 0.7–0.99 of a_max giving 0.1736–0.1101; Part H(d): on Part G(a)'s four-bump target the box stops the width coordinate at its face from s = 10⁹ |
| D3n-4 | minor (Status, §0, audit, (iv)) | the switching σ-modification already keeps the ideal properties inside its set, and at η = ∞ on collapsed classes the construction is that law (the exterior quadratic penalty); "a continuous form" and "dead zone" mislead | "new" restricted to the transfer to the learned centre of a mass-constrained K₀-atom prior at finite η and its compatibility with the bias removal; the switching σ's retention of the ideal properties cited ([4], [S52]: Ioannou 1986, Ioannou–Sun 1988); "a switching σ-modification toward 𝒵, the leak off on 𝒵", "the exterior quadratic penalty" |
| D3n-5 | cosmetic (numerics, rows, RESULTS, md) | L_ϑ^G = 4.658 and B = 26.08 are upper bounds (the supremum 4.0435) | "B ≤ 26.08 with the closed-form bound 4.658 of L_ϑ^G; 4.044 by search"; B a supremum over Θ and ϑ |
| D3n-6 | cosmetic (numerics (b), md, test 10, Status, RESULTS) | Part H(b)'s identity holds by construction (the box never approached); "passes 11 by s = 10⁴" is an unasserted comment | the md prints the margins (≥ 17% of the box's width) and says "by construction"; test 10 asserts the unprojected coordinate passing 11 by s = 10⁴ |
| D3n-7 | cosmetic (Status, (ii), proof) | (iii)'s local convergence inherits D‴.23(v)'s "sketched"; (ii) should name 𝒥_∘ and cover ties; the convexity step for split rows; the convergence at η = ∞ provable | "(iii) proved, its local convergence as D‴.23(v) (sketched)"; 𝒥_∘ with the cell masses and the tie case in (ii); the split rows' targets in 𝒵 by convexity; at η = ∞ the reduced energy is a minimum over the plans of C^{1,1} definable pieces (the squared distance to a box), so Proposition D‴.17's argument applies as it stands: the convergence proved at η = ∞, formal at finite η |
| D3n-8 | cosmetic (Remark D‴.25's Reading, D‴.24(v), the D‴.10 pointer, two audit bullets) | five sentences still said the adaptive prior gives up the a-priori bound; "at every weight" | "(restored by projecting the prior, Proposition D‴.27)" in all five; "at every representable weight" |
| D3n-9 | cosmetic (D‴.24 numerics, RESULTS) | "at s = 1250 the nearest-atom adaptive prior leads" is false at m = 15 | "leads at m = 60 and 240, not at m = 15 (1.75·10⁻² in the seed mean)" |
| D3n-10 | cosmetic (`constrained_floor`, md Part E, numerics, RESULTS) | the "constrained minimum" is the branch's local minimum, not the infimum (0 wherever a permutation is representable) | "the constrained minimum on the branch of ϑ^e(q), taken as the limit (the dynamics from the warm clouds follows the branch); a local minimum, the global infimum 0 where a permutation is representable" in the docstring, the md (regenerated), the text and RESULTS |
| D3n-11 | cosmetic (D‴.24 and D‴.23 numerics, §0, RESULTS) | the adaptive prior's coupled distance is heavy-tailed (one sample at 55.7 times the law at m = 240); "a failure set as large as" understates m ≥ 60; "of the same expected size" should be "consistent with" | the md prints the tails (the largest sample's D/law: fixed 3.36, 5.64, 7.47; adaptive 1.9, 8.5, 55.7; the adaptive mean without it 1.167 at m = 240), the text "a far heavier tail … not reliable"; "comparable at m = 15 and much larger at m ≥ 60"; "consistent with the same expected size (the means higher at m = 15, five draws)" |

**Lesson.** A numerical test of a mechanism must be run where the mechanism has something to do: on a reference with a symmetry that
kills the first-order term, "second order observed" confirms nothing about why. Moving the test to an asymmetric reference and to
k_e = 400 found more than was claimed — the expansion is even in 1/k_e, a consequence of the time-reversal symmetry of the adjoint pair —
and showed which hypothesis that needs (∂_xf̂ = ∂_xf on the orbit, the only rate even under time reversal). And for a projected flow the
integrator must project positions, not velocities: a velocity cut at a midpoint past a face creates a fixed point that is not the flow's,
and a single step size that happens to land on the face hides it — run a projected scheme at two steps and assert the end drift.

## Round 20 (adversarial review of Proposition D‴.28, Proposition D‴.29 with its servo application, and Corollary D‴.26(c) and Proposition D‴.27 as revised in round 19, 2026-10-20)

**Design.** Two referees in parallel on a snapshot of ee2ed27 (`git archive`), neither modifying the repository. D3o: Proposition D‴.28 and
Corollary D‴.26(c) as revised (`review/REVIEW_round20_D3o.md`): an independent re-implementation of the instance (analytic gradients and
Hessians, Newton, continuation, the averaged drift), the gauge checked by continuation from ϑ^b(q), 200 multi-starts per ε, the gauge
set's probability, vector amplitudes on a planar orbit, the separation eigenvalues of the dynamics, five Part B systems re-run bitwise to
s = 6·10⁴; for (c) an independent variational equation, Hill's method, a fourth-order Magnus integrator, the averaged generators by its own
FFT code, a derivative mismatch, the adaptation gain scaled. D3p: Proposition D‴.29 and the revised D‴.27 (`review/REVIEW_round20_D3p.md`):
120 BFGS multi-starts on the servo's loss and the saturated fits, the collapsed flows to s = 10⁸, the engineer's-model start, a no-leak
baseline, the exact-filter prediction decomposed, phase-plane kernels on ellipses and with centred atoms, a geometric search for
Lloyd-consistent non-optimal plans at η = ∞, the box scanned below 0.7, the projected step against a reference.

**Verdicts (D3o).** D‴.28: (i) holds exactly as a statement about the stationary points, but its numerical checks were tautological
(Newton started at the gauge transform) and "no floor whatever the approximation error" omitted a minimum, positive margins and the gauge
set's probability; (ii) holds (re-derived, Ã ≻ 0); (iii) holds for a scalar amplitude and fails for vector amplitudes; every number
reproduces. D‴.26(c): holds as a formal expansion whose structure is right; the sign of c_j is not fixed, the filter's part was mislabelled
at finite η, the parity's summary named the symmetry without the parity count. **Verdicts (D3p).** D‴.29: (i) and (iv) exact, (ii) holds
numerically and is rightly formal; the application's approximant is a local minimum, its runs started at the answer, and in the
noise-free simulation no leak fits as well as any transport design. D‴.27: (i) and the projected step hold; the convergence at η = ∞ was
not proved (the plan is Lloyd-consistent, and the reduced energy can rise); "between" holds only for boxes that contain the representation.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3o-1 | minor (D‴.28 numerics, §0, RESULTS, `exp_gauge_floor.py`, test 2) | the gauge check and the adaptive sensitivity were tautological: Newton started at the gauge transform takes no step ("distance 0.0" by construction), and f̂_w∘G_w = f̂_q holds for every Θ | Part A prints the gradient at G_wϑ^b(q) (≤ 1.6·10⁻¹⁵), continues the critical point in w from ϑ^b(q) (20 steps, refined on a fold: at G_wϑ^b(q) to ≤ 1.7·10⁻¹¹, f̂ to ≤ 5.5·10⁻¹³), counts plain Newton from ϑ^b(q) (at the gauge transform in 3–5 of 6 cases), and computes the adaptive sensitivity from critical points reached from ϑ^b(q) (1.4·10⁻¹¹–2.7·10⁻¹⁰ of the fixed prior's); test 1 asserts the continuation, test 2 starts Newton at ϑ^b(q) |
| D3o-2 | minor ((iii), proof, Status, rows, §8, audit) | (iii) needs a scalar amplitude: for vector amplitudes the amplitude–centre block a′_r∂_γΨ_r survives and does not scale (a planar counterexample misses the formula by 12–17%) | (iii) "for a scalar amplitude per kernel (one output)", with the vector case stated (only the centre–width block scales; no closed margin formula); the proof, Status, rows, §8 and audit likewise |
| D3o-3 | minor ((i), (iii), §0, §8, Status) | "no floor whatever the approximation error" omits a minimum, positive margins and the samples' probability of W_g; W_g shrinks, and the branch ends at the saturation | (i) about the stationary points, the floor vanishing at the samples of W_g with positive margins for a strict local minimum when the systems converge (numerical); (iii)'s admissible set for a strict local minimum; Part A prints P(w^m ∈ W_g) (0.52, 0.87–0.88, 0.997 at m = 15, 60, 240 for ε ≤ 0.3; 0.26, 0.71, 0.98 at ε = 0.5) and the branch continued to ε = 0.68 (max|a|/a_max = 0.999, min eig Ã 1.2·10⁻⁸, the threshold 0.333 = q₂); quoted in the Status |
| D3o-4 | minor (numerics, §7, RESULTS) | at the sample outside W_g the quoted "kept errors" were transients still falling at s = 6000; a gauge transform with two classes exchanged exists | a new `outside_job`: the minimum of L_w over the closed box from the clipped gauge transform (1.107 and 1.586 times the approximant's fit error at ε = 0.1, 0.3; one amplitude at the bound) and the permutations with a gauge transform ((1, 3, 2), (2, 3, 1)); the text: the systems approach that minimum, the permuted transform exists but is not reached |
| D3o-5 | minor ((c), §0, §8, numerics, RESULTS, test 3) | the sign of c_j is not fixed (the slaving's −(J̄ᵀJ̄)² and the level repulsion are ≥ 0); with the gain scaled by 30 two c_j > 0 | (c) and the proof: c_j the filter's part (≤ 0) plus the slaving's (≥ 0) plus averaging terms, sign not fixed; Part C adds the gain scaled by 30 (c₈/λ₈ = +4.17, c₉/λ₉ = +12.65; those exponents above their first-order prediction by c_j/(λ_jk_e²) to three digits); "every exponent below" an instance property in the md, text, rows, RESULTS and test 3's comment |
| D3o-6 | cosmetic (numerics, md, (c)) | "−102.1 = −ω₁²" at η = σ₀ is mislabelled (ω₁ does not depend on η; the filter's part is −0.989ω₁² there) | ω₁ from the η-independent Rayleigh quotient (10.1625, 10.5578), the filter's part −κ₁ω₁² with κ₁ printed; (c) "−(dλ_j^±/dh)v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j" |
| D3o-7 | cosmetic ((c), proof, §0, §8, audit) | time-reversal symmetry alone does not make the expansion even (it persists with a ≢ 0); D singular at η = ∞; degenerate eigenvalues | (c) and the proof: the same multipliers and the parity (−1)^k of every term of Ā_k, which a ≢ 0 breaks; D = diag(−I_n, I) at η = ∞; "for simple eigenvalues; at a degenerate one only the rates are even" |
| D3o-8 | cosmetic (numerics, rows, RESULTS, docstring, §8) | 0.11σ₀ (0.10); the rows' 7·10⁻⁴–3·10⁻³ (to 5.6·10⁻³); "constant from s = 600" (within 8%); the margins table checks a curvature identity; g̃'s sign in a docstring; "the critical point" | 0.10σ₀; 7.0·10⁻⁴–5.6·10⁻³; the md prints the ratio h2(6000)/h2(600) (0.917–1.003) and the text says within 8%; Part A adds the dynamics' separation eigenvalues (to ≤ 5.8·10⁻¹²); the docstring g̃ := −D^{1/2}∂_w(G_wϑ^b); "a critical point" |
| D3o-9 | cosmetic (`tests/test_gauge_floor.py`, `tests/test_adaptive_floquet.py`) | tautological assertions, a fragile amplitude-row tolerance, η = ∞ only; an Ā₁ = 0 mutant passes test 4, the metric never tested, `var < p1` an instance property | test 1 asserts the continuation, test 2 Newton from ϑ^b(q), test 3 a tolerance with margin and the dynamics' separation eigenvalues (a mutant without q_k/w_k misses by 5·10⁻⁵), test 4 also at η = σ₀ (to s = 1.2·10⁵); the Floquet test 4 asserts ℓ_jĀ₁r_j = 0 and, at η = 10σ₀, Ā₀ self-adjoint in the metric and not in I, Ā₁ skew in it; `var < p1` commented as an instance property |
| D3o-10 | cosmetic (Status, rows, audits) | the gauge's precedents (rescaling symmetries) uncited; the even expansion is the classical time-symmetry argument | [S55] (Dinh et al. 2017, Neyshabur et al. 2015, Şimşek et al. 2021) and [S56] (Hairer–Lubich–Wanner 2006) cited; "new" restricted to the mechanism, the scalar-amplitude margin formula, the lift at every representable sample and the explicit c_j |
| D3p-1 | major (D‴.29 numerics, §0, §7, §8, RESULTS, PAPERS paper 11, `exp_servo_cogging.py`) | ϑ^b(q) is a local minimum, not the best approximant: the best fits saturate amplitudes and reach 0.64%; globally the servo is D‴.27's case, and the friction/cogging contrast in §8 is wrong | ϑ^b(q) "the nondegenerate local minimum that BFGS reaches from START" throughout; Part A adds the global picture (120 random starts: 1 at ϑ^b(q), 33 below its fit error, all saturated; the best 6.39·10⁻³; with its two saturated amplitudes held at |α| = 6–14 the fit error 7.14·10⁻³–6.40·10⁻³: the infimum approached at the saturation); the text: globally Proposition D‴.27's case (critical points at infinity), Parts B–C statements about the basin and about an engineer's model; the §8 contrast dropped |
| D3p-2 | major (numerics (B), §0, §8, RESULTS, PAPERS) | every design started from clouds around ϑ^b(q), computed from the true torque; at η = ∞ the misplaced prior never acts; from the engineer's model the bias decays only slowly | Part B adds the engineer's start (clouds around the engineer's model, the prior at it) for every design and the collapsed K₀-atom flows to s = 10⁵ from both starts: from the engineer's model the adaptive prior's fit error is 0.027, 0.025, 0.022 at s = 2000 (the fixed designs 0.10, 0.060, 0.048; OT 0.048), and the K₀-atom flow at η = ∞ needs s ≈ 10⁵ (0.0297, 0.0213, 0.0169, 0.0155 at s = 600, 2000, 10⁴, 10⁵); the text says the matching prior is ϑ^b(q), that at η = ∞ the misplaced runs are the matching ones, and that the bias removal from the engineer's model is real but slow |
| D3p-3 | major (D‴.27(ii), proof, Status, §0, §8, LITERATURE; D‴.23(iv)) | at η = ∞ the flow's plan is Lloyd-consistent (optimal for its own barycentres), not jointly optimal: the reduced energy min_Z𝒥 can rise (+30%) while the code's energy falls, so the convergence "by Proposition D‴.17's argument" was not proved | (ii): the selected energy 𝒥(Θ, Π_𝒵T̂(γ)) is non-increasing (also at plan switches), the reduced energy need not be; the convergence at η = ∞ argued (a Kurdyka–Łojasiewicz argument on the selected pieces, not carried out), its limit a critical point of a Lloyd-consistent piece; D‴.23(iv)'s formal limit qualified likewise; Part H(e) (a seeded search, the code's own run: the selected energy falling, largest increment −3.4·10⁻⁶, the reduced energy rising from 1.251·10⁻³ to 1.574·10⁻³ before the plan becomes optimal) and a test of both |
| D3p-4 | minor (Part B, §0, RESULTS, PAPERS) | no no-leak baseline: with σ₀ = 0 the fits and floors are at or below the adaptive prior's — the leaks' robustness is never exercised; the σ-modification compared pulls toward 0; a Fourier regressor of the harmonics represents τ_cog exactly | Part B adds "no leak (σ₀ = 0)" from both starts (oracle: 0.0180, 0.0147, 0.0146, floors 7.7·10⁻⁶, 2.4·10⁻⁷, 6.6·10⁻⁸; engineer: 0.0152, 0.0147, 0.0147 — the m-particle flow spreads and learns faster than the K₀-atom structure the leaks impose); the text says so — in the noise-free simulation the floors and biases the leaks remove are those a fixed leak creates, their purpose (robustness to noise and disturbances, Propositions D‴.27, D‴.30, D‴.31) not exercised; "the σ-modification toward ϑ⁰ = 0 at the transport leaks' σ₀"; paper 11's "Missing": the classical baselines, a disturbance or noise, the Fourier regressor of the known harmonics |
| D3p-5 | minor (D‴.29 (ii), (iii), proof, Status, §0, Part C) | "Theorem J's argument transferred" overstates (Theorem J is local, needs a smooth leak); the prior's real-time law unstated and frozen in Part C; (iii)'s list includes D‴.24(vi) | the proof and Status: first-order averaging with the fast pair (s, e), Theorem J's steps carrying over for the smooth fixed-prior leak near a hyperbolic critical point once the error-channel steps are redone, switching targets beyond even the first-order theory; (ii) states dz/dt = −(η/k_s)(z − T̂) and Part C integrates it; (iii) lists D‴.24(i)–(v) and names D‴.24(vi) among the real-loop statements |
| D3p-6 | minor ((ii), Part C, RESULTS, §0, LITERATURE) | the exact filter is not uniformly the better prediction (both O(1/k_s) relatively, the parameters' motion at the same order); Part C's records 0.2 apart; the residual's spectrum broad, not sin 7x | (ii): both forms accurate to O(1/k_s) relatively, the exact filter better where the dropped harmonics dominate; Part C records every 0.02 time units and prints both predictions at every window: the exact-filter ratios 1.0001–1.0029 at t = 500 and 1000 in every design (0.990–1.005 over the first five periods), the slaved 0.931–0.996; the residual described as the approximation error's broad fine scale; the filter dropped from the audit's "specific" list |
| D3p-7 | minor ((iv), Status, §0, LITERATURE, PAPERS) | the exact zero needs ω₀ = 1 in isotropic units; the weak excitation is generic on other orbits; "one per kernel" only in general position (centred kernels: more) | (iv): ω₀ = 1, "in general position (more at degenerate atoms)", and on other orbits weak rather than zero excitation, the classical partial excitation [6]–[8]; Part D tabulates it (five kernels: the circle 5 exact zeros and 8 below 10⁻⁶σ₀; ellipses ×1.1 and ×2: 8 and 6 below 10⁻⁶σ₀, none exact; a two-harmonic orbit 3; centred kernels 16 exact zeros of 20; the (K + 1)-th eigenvalue 2.4·10⁻¹⁰ of the largest); "new" restricted to the exact fibre |
| D3p-8 | minor (D‴.27 numerics (c), rows, RESULTS, §8) | "between" holds only on the scanned bounds; tighter boxes are worse than the fixed prior; the 0.7 box excludes the exact atom | Part H(c) scans 0.5–0.99 and prints whether the exact atoms are inside: 0.75–0.99 between the fixed prior's 0.1750 and the unprojected 0.1050, 0.7 excluding the exact atoms (0.1736), 0.65, 0.6, 0.5 worse than the fixed prior (0.1792, 0.1845, 0.1949); "a box should contain the representation" |
| D3p-9 | cosmetic (D‴.29 numerics, rows) | "below the replicated point's 0.0545" false for the nearest atom at m = 16, 64; "fit error to the approximant's" fails at m = 16; three modes far below σ₀, not two; the §7 row omits 7.8·10⁻⁵; n²ω₀² | the numbers corrected (the nearest atom 0.10, 0.060, 0.048, above 0.0545 at m = 16 and 64; the m = 16 adaptive prior not yet replicated; four modes below 0.1σ₀); the §7 row completed; n²ω₀²/(2k_s²) |
| D3p-10 | cosmetic (D‴.27) | (iv) needs T̂_j ∉ 𝒵; "at least 17%" is 16.9%; the two steps' agreement shows the fixed point, not consistency | (iv) "with its class barycentre outside 𝒵", the boundary-stationary case named; 16.9%; Part H(a) prints the trajectory error at s = 300 against dt = 1/64 (6.3·10⁻⁷, 1.9·10⁻⁷, 3.8·10⁻⁸ at dt = 0.25, 0.125, 0.0625) |

**Lesson.** A check that starts at the answer tests nothing. The gauge check started Newton at the gauge transform, the servo's runs at
the approximant computed from the true torque, and at η = ∞ the misplaced prior never acted. Each fix replaced the answer by what one
would actually have — a continuation from ϑ^b(q), an engineer's model, random starts — and the picture that results is less flattering and
more useful: the gauge result holds at a strict local minimum; the servo's approximant is a local minimum of a loss whose infimum lies at
the saturation; the adaptive prior removes a misplaced prior's bias, but at the rate of the weakest excitation; and in a noise-free
simulation no leak at all fits as well, because what the transport leaks buy — structure and robustness — needs noise or disturbances to
show. Likewise a Lyapunov function must be checked along the flow the code integrates: at η = ∞ the code's plan is Lloyd-consistent, so
its own energy decreases while the reduced energy the proof named can rise. And a baseline must be checked to be what it is named: the first
no-leak run of the real servo inherited the global σ₀ — a fixed leak in all but name — and its error crept toward the fixed design's
(2.0·10⁻⁵ → 1.4·10⁻⁴ at m = 16), a "drift" that the corrected run (the design's σ₀ = 0: 1.88–1.90·10⁻⁵, at its prediction to
1.0018–1.0020) does not have.

## Round 21 (adversarial review of Propositions D‴.30 and D‴.31, measurement noise at perfect tracking and its remedy, 2026-10-20)

**Design.** Two referees in parallel on a snapshot of 637a051 (`git archive`), neither modifying the repository. D3q: Proposition D‴.30
(`review/REVIEW_round21_D3q.md`): an independent instance (analytic features, complex-step derivatives) and linear-periodic machinery (RK4,
a Kronecker discrete Lyapunov solve), Part B to k_e = 800 at two resolutions and a control with the slope removed from the noise input,
the exact stationary second-order mean of the linearised loop with every term, the separations' second-order drift and their finite-time
exponents along noisy trajectories, an independent batched real loop reproducing the stored runs, antithetic noise pairs, test mutants.
D3r: Proposition D‴.31 (`review/REVIEW_round21_D3r.md`): an independent instance, the washout identity along a path with jumps and for an
OT plan with split rows, (ii) beyond the commuting family, (iii) for non-symmetric S, the exact second-order small-noise mean of the real
loop at the lift (every quadratic term, split into the noise and the parameter blocks of the covariance), the exact second-order
transient of Part C, its seed reproduced bitwise, test mutants.

**Verdicts (D3q).** D‴.30: (i) holds exactly, as invariance — not as attraction; (ii) holds as a leading-order formula with an absolute
O(k_e⁻²) remainder, its relative error O(1/k_e), not O(k_e⁻²); (iii) holds for the averaged system, the force right in sign and size
(confirmed in the real loop mode by mode by antithetic pairs); every number reproduces. **Verdicts (D3r).** D‴.31: (i) holds (the variation
of constants; the leak's decomposition for integral plans); (ii) holds as RH + B = I — HR + B = I only where S₀ commutes with H; (iii)
holds, and is classical; (iv)(a)–(b) hold; (iv)(c) holds for the truncated force and omits the parameters' own fluctuation, of the same
order — the one major issue.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3r-1 | major (D‴.31(iv)(c), proof, Status, §0, §8, numerics, RESULTS, PAPERS 12, README, test 5) | the second-order force omits the parameters' own fluctuation — ∇²_ϑΦ̃E[ε_mδϑ] in the law, the compensation's curvature and slope–parameter terms in the error equation: DR's bias +3.3% at k_e = 25 and +56% at k_e = B = 100, above ‖g‖²_ρ; "the only part … that the adaptive prior amplifies", "uniformly in the spectrum" and the extreme 3.3·10⁸ are properties of the truncation | `second_order_mean` (the periodic covariance of the linearised loop with the OU noise as a state, contracted with the true field's second derivatives; the noise and the parameters' parts; the exact periodic mean) and `mean_biases` (the biases of the mean parameters and of the mean function): Part A prints both parts at every (k_e, B, law) under both priors — DR +3.3% (mean function +1.2%) at k_e = 25, +56% at k_e = B = 100, DR* −8.7% to −0.15%, MR ≤ 0.03%, the parameters' part's own response 3.9·10⁻⁴ → 2.9·10⁻⁵ from k_e = 25 to 100; (iv)(c) restated for the noise part, the parameters' part named (not least-squares, set by the parameters' covariance, unbounded as h_j → 0), the closing paragraph and the Status ("formal and split"); MR over DR 2.1·10⁸ at the extreme; a test (the noise parts against the closed forms; the parameters' part above 1% for DR and DR*, below 10⁻⁵ for MR; DR's full bias three orders below MR's) |
| D3q-1 | minor (D‴.30 Status, numerics, §0, md B) | "the deficit falling like k_e⁻²" holds only on k_e = 25–200: the relative error has a +0.03/k_e term from the slope a in the noise input −(k_e + a)ν, which breaks the time-reversal parity, and the ratio crosses 1 near k_e = 500 | Part B to k_e = 400, 800 (the ratio above 1 at 800: +1.2·10⁻⁵ and +1.3·10⁻⁵ at η = 0, ∞; the deviation's successive ratios 3.87, 4.38, 5.26, 11.4, then a sign change) and the control with the noise input −k_eν (successive ratios 3.67, 3.90, 3.98, 3.99, 4.00: even in 1/k_e); the Status: the deficit dominated by a k_e⁻² term to k_e = 200, the relative error O(1/k_e), its first-order term the slope's |
| D3q-2 | minor (numerics, rows, RESULTS, PAPERS 12) | the averaged system's values quoted as the instance's at k_e = 25: the weakest mode's ×24 is ×28.0 in the exact linearised loop, the biases 491 and 2397 are 504 and 2653 (×5.27); "which the transient checks" holds for the force, not the stationary response | Part B prints the weakest mode's variance against its slow-plus-fast prediction (0.845 fixed, 0.980 adaptive at k_e = 25; adaptive over fixed 28.0, 25.1, 24.4, 24.2, 24.1, 24.1 at k_e = 25–800) and the exact stationary second-order mean of the linearised loop (503.7 and 2653.0 at k_e = 25, the excess 10.7%, 2.65%, 0.66%, 0.17% under the adaptive prior and 2.58%–0.04% under the fixed one: O(k_e⁻²)); the text quotes both, "the transient checks the force, the stationary response is checked against the exact linearised loop" |
| D3q-3 | minor ((i), Status, numerics, md, §0, §8, RESULTS) | invariance is not attraction: under noise a separation decays at σ₀/k_e − Eν²λ_max(S_k), class 3 losing its margin at std ≈ 0.03; the lift's "0.0" is the float32 record and identical arithmetic per particle | `separation_drift` (S_k: λ_max 0.0148, 0.0101, 0.298 per Eν²; the threshold std 0.031; at Part C's stds the margin used to 0.9% and 10.5%); a new Part D — the lift against the K₀-atom run in float64 (5.6·10⁻¹⁶; the within-class spread 0.0 by identical arithmetic) and the lift's third class separated by 10⁻⁶ along S₃'s top eigenvector, four noise paths of 300 time units: −2.55·10⁻⁴ ± 0.04·10⁻⁴, −0.9·10⁻⁵ ± 1.4·10⁻⁵, +4.37·10⁻⁴ ± 0.26·10⁻⁴ per time unit at std 0.01, 0.03, 0.05 (predicted −2.53·10⁻⁴, −1.5·10⁻⁵, +4.62·10⁻⁴); (i) states the threshold; a test (a separation decaying at 0.01, growing at 0.05, each at its predicted rate to 10⁻⁴) |
| D3q-4 | minor (D‴.30(iii), proof) | the parameters' fluctuation's second-order terms omitted (D3r-1); ∂²_xf̂ = ∂²_xf needed and unstated; "grows like k_e²" not realised at B = 4k_e | the hypothesis adds ∂²_xf̂ = ∂²_xf; (iii) and the proof name the parameters' terms, computed by `second_order_mean` (MR: 2397.66 → 2397.68); Part A prints the bias's local exponents in k_e at B = 100, 10³, 10⁴ (1.69, 1.47, 1.17, 0.83, 0.53 at B = 100; 1.93–1.99 for B ≥ 20k_e); "while k_e ≪ B, saturating for k_e ≫ B" |
| D3q-5 | minor (tests, Status) | the force untested (a sign flip, ×10, no E[e], E[ε_m²] = Eν² pass); test 3 passes without the noise in the law and with a wrong OU variance; test 2 at η = ∞ only | test 4 compares the quasi-static force with the full second-order mean's noise part — the periodic covariance contracted with the true field's second derivatives, independent of the closed forms — to 2·10⁻³ (1.7·10⁻⁴ here; the mutants miss by 25% or more); a test of the OU noise's variance and lag-one correlation (a stub recording ν) and of the law at the measured error and point; test 2 at η = 0, σ₀, ∞; "linear in the variance" dropped from the Status |
| D3q-6 | minor (Status, §0, audit, PAPERS 12) | the novelty overstated: the LMS misadjustment, leaky LMS, the diffusion approximation, common noise and the lack-of-excitation drift omitted | [S57] (Widrow et al. 1976, Mayyas–Aboulnasr 1997, Benveniste–Métivier–Priouret 1990, Carmona–Delarue–Lacker 2016, Anderson 1985; bibliographic details checked by a web search); "new" restricted to the floor's interpolation in η, the driving matrix Q of the sensitivity-filtered measurement and their combination with the coherence and its noise threshold |
| D3r-2 | minor (D‴.31(ii), proof, §0, §8, PAPERS, `stationary_operators`, test 2, md B) | HR + B = I is false for a zero-frequency leak not commuting with H (residual 0.21–0.79); the identity is RH + B = I; the trade binds only on h_j < σ₀; the anchor's normalisation unstated | (ii) "RH + B = I, i.e. R = (I − B)H⁻¹" (HR + B = I where S₀ commutes with H), per mode h_jr_jj + b_jj = 1, "the trade binding on the modes with h_j < σ₀" (3 of 9), the anchor entering through I − 𝒦; `zero_frequency_operators`; Part B adds four non-commuting leaks (RH + B = I to 1.4·10⁻¹⁴, per mode to 8.5·10⁻¹⁵, HR + B − I 0.36–0.79); test 2 asserts both for the active face and the fixed amplitudes |
| D3r-3 | minor (Status, §0, audit, PAPERS) | the novelty overstated: Yucelen–Haddad's low-frequency learning and elastic averaging's moving centre anticipate the washout reading; the identity is the Tikhonov/Kalman split; (iii) and the F₁ split classical | [S58] (Yucelen–Haddad 2013, with [S49]'s elastic averaging); the identity, the bound and the split labelled classical; "new" restricted to their transfer to the transport leaks (uniform in m) and the quantification; DR* without Sadegh–Horowitz's nonlinear damping a local statement |
| D3r-4 | minor ((iv)(c), numerics, RESULTS, §7, md A) | "the projected mean tracking error −ĴH⁻¹F₁ grows like k_e²" is k_e times it; the md's zeros for DR and DR* hold by construction | (iv)(c): P E[e] = −ĴH⁻¹F₁/k_e, what grows like k_e² is k_e times it, the mean tracking error itself bounded (∝ (B/(k_e + B))²); the md labels k_e²‖P E[e]‖²_ρ (zero by construction) and prints ‖P E[e]‖²_ρ (3.83 → 0.66 from k_e = 25 to 200, MR, B = 100) |
| D3q-7 | cosmetic (numerics, proof, §8, md) | the cross term vanishes on any 1-D reference; "falls from 104 … to 5.5" not monotone; "at most 10%" (10.02%); "5–6%"; the four Part C cases share their noise paths; a stale §8 "Open"; the proof's notation, n = 1; a mode's variance without its fast part | "an exact differential"; "ranges from 104 (mode 1) to 5.5 (mode 6)"; "10%"; the md prints the priors' drift (6.0%, 6.4% at t = 300; predicted 5.3%, 6.1%) and the common paths; the §8 "Open" updated; n = (d/dt + k_e)⁻¹[ν̇ − aν] and y = Ĵᵀξ + u in the proof; n = 1 in the hypothesis; the fast part (N₀/2k_e)h_j |
| D3r-5 | cosmetic (D‴.31(i), proof) | bounded variation (a chattering plan); split rows (−σ₀(I − M)z not washed out); "at a stationary point"; "no leak opposes a persistent force"; the washout analogy only for replicated equilibria | (i): the integrated form for every measurable T̂, the differential one for bounded variation; integral plans, split rows' −σ₀(I − M)z named; "at the exact representation"; "on the class barycentres (the separations keep the full leak)"; "here the replicated ones" |
| D3r-6 | cosmetic ((iii), docstrings) | symmetry of S not needed (S + Sᵀ ⪰ 0 suffices); H must be invertible | (iii) "for H ≻ 0 and every S with S + Sᵀ ⪰ 0", the proof with (H + S)ᵀH⁻¹(H + S); `response_size`'s and `equivalent_perturbation`'s docstrings |
| D3r-7 | cosmetic (numerics, md, rows, RESULTS, PAPERS, docstrings) | averaged values; "grows like k_e²" is ×46.6; the md's "(1 + O(1/k_e))"; the crossover std 4.41; "three to four orders"; the reading of the sampling error; `mean_error`'s docstring; DR*'s jitter; PAPERS' "(ii) proved", "remedy checked" | "like (k_eB/(k_e + B))²" (×46.6, the law ×46.7); the md's line without "(1 + O(1/k_e))"; crossovers with the full bias, "stds above 0.1 outside the small-noise regime"; "MR over DR 2.3·10³–5.1·10³ along the family"; the md prints the error over the relative sampling error (0.43–1.83), and the text says the averaged prediction carries an O(1/k_e) error and the test checks F₁, not DR's own force; `mean_error`'s docstring; DR*'s jitter k_e²Eν²/((k_e − a)(k_e − a + B)), k_e > max a; PAPERS 12 rewritten |
| D3r-8 | cosmetic (tests) | test 1 never exercises a plan jump; test 2 only the commuting family; test 4 passes a regressor at the true state; test 5 checks the force against its own definition | a test of the washout identity on the averaged loop across a plan switch (T̂ jumps by 0.49; the identity to 2.4·10⁻⁶ with the jump's share large); test 2 with the active face; test 4 checks the regressor's point (DR's law at x ≠ x_m is G∇_ϑ(eΦ̃(x_m; ϑ))); the forces against the full second-order mean (D3r-1) |

**Lesson.** A bound proved for a truncated force is a bound on the truncation. The second-order force had been computed from the noise's
own moments at fixed parameters, and every statement about it — the least-squares form, the bound by ‖g‖²_ρ, the extreme ratio where that
force cancels — was exact for the truncation; the parameters' own fluctuation, of the same order in the noise, is a second
errors-in-variables term, in the parameters, which putting the reference in the regressor does not remove. The remedy was to compute the
second-order mean in full and generically — the periodic covariance contracted with the true field's second derivatives — which both
checks the hand-derived formulas (the noise part to 1.7·10⁻⁴) and measures what they leave out. The same round showed twice that the
observed behaviour on a range is not the asymptotic one: an order read off k_e = 25–200 was a large k_e⁻² term masking a first-order term
of the opposite sign, which a control without the parity breaker exposed; and a subspace invariant under every noise path is not
attracting — above a noise level the replicas separate, which a lift started exactly on the subspace could never show.

## Round 22 (adversarial review of Lemma D‴.26a and of the round-21 revisions of Propositions D‴.30 and D‴.31, 2026-10-20)

**Design.** Two referees in parallel on a snapshot of 3e8f8e5 (`git archive`), neither modifying the repository. D3s: Lemma D‴.26a and
the places it upgrades — Corollary D‴.26(c), Proposition D‴.30(ii) (`review/REVIEW_round22_D3s.md`): an independent instance (Ĵ and its
x-derivatives to order four in closed form, K(t) and its t-derivatives by Faà di Bruno), order-8 Gauss–Legendre propagation of the linear
periodic systems, P and S from the right and left slow subspaces propagated forward and backward (no contraction), the periodic covariance
by Van Loan's block exponential and a Kronecker discrete Lyapunov solve, the parity's two halves as matrix identities and on cases of its
own (n = 2, complex pairs, an indefinite D, a control without the symmetry), Ā₁ by its own quadrature, the noise at η = 0, σ₀ and ∞ with
the covariance's deviation decomposed, test mutants. D3t: the round-21 revisions of D‴.30 and D‴.31 (`review/REVIEW_round22_D3t.md`; its
issues ID-1–ID-6 are D3t-1–D3t-6 here): hyper-dual Hessians and an independent periodic covariance for the full second-order mean (36
cases), the separation drift's parts and the exact second-order transient from the lift, the K₀-atom loop with the linearised separation
on 24 antithetic noise pairs per level, the dynamic block of a non-commuting anchor gain, the parameters' part's scaling in k_e, 13 code
mutants against the tests.

**Verdicts (D3s).** Lemma D‴.26a: (i)–(iii) hold; (iv) holds as an identity, and for the whole spectrum of the slow generator, not only the
simple real rates; (v) holds asymptotically, the floor as a period average; every number reproduces with independent code; no major
issue. **Verdicts (D3t).** D‴.30: (i) holds exactly as invariance, but the threshold formula holds at the lift only — the one major
issue; (ii) holds, its order supported to k_e = 1600; (iii)'s computation is right (36 cases independently). D‴.31: (i)–(iii) hold, (ii)
as RH + B = I, the non-commuting check an identity of the definitions; (iv)(c) holds as restated. Both mds regenerate byte-identically.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3t-1 | major (D‴.30(i), Status, *Numerics* (A), (D), §0, §7, §8, RESULTS, PAPERS 12, LITERATURE, `separation_drift`, md D) | the separations' drift S_k evaluates E[e] as the noise's own mean error at fixed parameters; along the noisy trajectory E[e] also carries the mean tracking error −Ĵȳ(t)/k_e of the developing errors-in-variables deviation, of the same order: the stated threshold is the lift's, the small-noise λ_eff measured 0.320 against the lift's 0.298 (+7.5%, 9 s.e.), λ_max(S₃(t)) rising 0.306 → 0.343 over 250 time units; Part D's agreement at std 0.05 a cancellation against the next order; the long-time threshold outside the expansion (the stationary mean deviation O(1) at these noise levels) | (i) restated — the decay rate σ₀/k_e − Eν²λ_max(S_k) and the threshold at the lift, along the trajectory u̇ = (Eν²S_k(t) − σ₀/k_e)u with E[e] carrying −Ĵȳ(t)/k_e, the long-time threshold open; `second_order_transient` (the second-order transient from the lift: the covariance of the augmented linearised loop, the stencil forcing of the true field, the mean) and `separation_drift_series`; Part D on eight antithetic pairs per std with λ_eff and its standard error against the lift's and the transient's predictions: the separations' rates −2.50·10⁻⁴, +2.74·10⁻⁶, +4.53·10⁻⁴ at std 0.01, 0.03, 0.05, λ_eff = 0.3271 ± 0.0038, 0.3175 ± 0.0037, 0.2943 ± 0.0034 against the lift's 0.2981 and the transient's 0.3264, 0.3270, 0.3281 (along the transient the mean error's rms 0.121 → 0.602 per Eν² and λ_max(S₃) 0.307 → 0.349 over 300 time units), the fit λ_eff = c₀ + c₁Eν² with c₀ = 0.3291 ± 0.0033 — the transient's to within its standard error, 10.4% above the lift's — and c₁ = −13.8 ± 2.0, the rate vanishing at std 0.0299 (0.0308 by the lift's drift, 0.0294 by the transient's); the Status "formal at the lift … the long-time threshold open", "at their predictions" replaced by the λ_eff figures, the novelty "its noise threshold at the lift"; a test (along the transient from the lift the mean parameters moved by the averaged force to 1.5% after one period, E[ε_m²] quasi-static to 0.5%, and in the third period the mean error the noise's own plus −Ĵȳ/k_e, the tracking part 2.6 times the noise's own) |
| D3s-1 | minor (lemma *Numerics*, §0, §7, RESULTS) | "‖Cov(Y′) − Σ̄‖ dominated by the weakest mode's rate correction, O(ε²)": the slow noise's attenuation ⟨bbᵀ⟩ < Q̄ dominates, 2.5 times the rate's correction and of the opposite sign, and beyond k_e = 400 the deviation is O(ε) | Part D to k_e = 1600 with the weakest mode's deviation split into "noise only" (Σ with ⟨bbᵀ⟩) and "rate only" (Σ with the exact slow generator): −56% and +49% at k_e = 12.5, 2.4–2.5 to 1 from k_e = 100 at η = ∞ (12–23 to 1 at η = 0), their sum the exact deviation (O(ε²) in that mode); ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ falling by 3.6–4 per doubling to k_e = 200–400 and by 2–2.9 beyond, "the O(ε) of (v), with a small constant" |
| D3s-2 | minor ((v), D‴.30(ii), Status, *Numerics*, md D) | the floor holds as a period average (pointwise ±1%: the fast jitter's oscillation is of the floor's order); the jitter untested; the remainder's next term | (v) and D‴.30(ii) "averaged over a period", the pointwise fast part (N₀/2k_e)tr(HĴ(t)ᵀĴ(t)) named; the next term of Cov(Y) − Cov(Y′), −(ε²N₀/2)(d(KᵀK)/dt + 2aKᵀK), in (v); `decoupled_covariance(period_average=True)`; Part D prints the floor's period average (0.766 → 0.99997 of the formula, η = ∞, k_e = 12.5–1600), its range over the period against the jitter's predicted range ([0.9952, 1.0091] both at k_e = 1600) and the jitter against Cov(Y) − Cov(Y′) without and with the next term (0.32 → 0.020 and 0.30 → 0.001 at k_e = 100–1600) |
| D3s-3 | minor (the identification, D‴.30's Status, `decoupled_covariance`) | η = 0 is outside the stated identification (the pair's Λ has a zero rate there, Ā₀ not Hurwitz); the fixed prior needs K = Ĵ, Λ = σ₀I | the identification adds η = 0 (K = Ĵ, Y = y, Λ = σ₀I); `CollapsedVariational` treats η = 0 as the fixed prior (the state (δe, y); the generator, the exponents, the averaged and second-order blocks, Chang's grid); Part D at η = 0 (the floor's period average 0.90325, 0.99216, 0.99949 at k_e = 25, 100, 400, the referee's to 3·10⁻⁵) |
| D3s-4 | minor (the lemma's ε₀, D‴.26(c), its Status, §0) | ε₀ and the constants also depend on the spectral gap and the eigenvector's condition ((iv)) and on Ā₀'s margin ((v)); D‴.26(c) needs x_m ∈ C⁴ and, for the stability, H_L(w) ≻ 0 | the ε₀ clause lists both; D‴.26(c), its Status and the §0 row "for x_m ∈ C⁴ at simple eigenvalues, the stability where H_L(w) ≻ 0" |
| D3s-5 | minor ((v), Status, §0, §8, LITERATURE, D‴.30's Status) | "the slow noise identified exactly as the high-passed −ε(K̇ + aK)ᵀν" overstated: the exact slow noise is εb(t)ν, b = −(K̇ + aK)ᵀ − ε(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ + O(ε²), the drift decoupled, the noise shared | (v) states b with its two terms and "the drift decoupled, the noise shared" (E[vY′ᵀ] ≠ 0); the Status, the rows, the §8 bullet, the audit; `noise_decoupling(with_b1=True)`; Part D prints |b − b₀ − εb₁|/ε² (1842 → 3341: bounded); a test (b₁ at ε²) |
| D3s-6 | minor (tests, Status) | the lemma's tests partly vacuous: P₂ with the leak's sign flipped or the cubic term dropped passes (on the instance K̈ dominates both by orders of magnitude), and so do Σ̄ without the slope and the slope dropped everywhere | a synthetic block (n = 1, q = 3, low harmonics, an O(1) non-symmetric Λ) where every term is visible: P − P₂ at ε³ (ratios 6.5–9.5), P − P₍₃₎ at ε⁴ (13–19), the two mutants at ε² (3.5–4.6); the noise test with b₀ against K̇ by the chain rule and the plant's own slope by differences of f, Q computed independently, b₁ at ε², the covariances' orders |
| D3s-10 | minor (Status, §0, audit, [S59]) | the closest precedents omitted: Chang's transformation for stochastic and time-varying singularly perturbed systems (Kokotović–Khalil–O'Reilly Ch. 4–5 and App. A; Haddad–Kokotović 1977; Yang–Zhu 2010, 2012), stochastic slow manifolds (Berglund–Gentz 2003, 2006), the adjoint relation of symmetric methods | [S59] extended; the Status, the §0 row and the audit bullet name them; "new" kept to the application: the exact slaving, the parity as an identity at negative ε, the slow noise identified |
| D3s-7 | cosmetic ((iv), its proof, D‴.26(c)) | "which holds when a enters only the noise" misleading (the D-symmetry persists with a mismatch; what a breaks is P's invariance under (t, ε) ↦ (−t, −ε)); the proof gives the whole spectrum even; a shorter route through the time reversal | (iv): "the spectrum of Ā(ε) is even: eig Ā(−ε) = eig Ā(ε)", the symmetry "persists with a mismatch in the error channel", the parity needing a ≡ 0 there; the proof by the reversed reference's system, without P |
| D3s-8 | cosmetic (*Numerics*, §0, §7, RESULTS, md D) | the parity and the deviations normalised by the largest rate; per mode the slowest mode deviates 29% at k_e = 12.5 | Part D prints the per-mode columns and the slowest mode's signed deviation (−29.6%, −12.5%, −4.0%, −1.1% at k_e = 12.5–100, η = ∞; the parity per mode ≤ 1.1·10⁻⁹); the text quotes both normalisations |
| D3s-9 | cosmetic ((iii), proof) | the Magnus series "of εA collected by powers of ε" (A has only an asymptotic expansion; the proof uses εA₍₃₎); "conjugate up to O(ε⁴)" (D‴.26's averaging is carried to Ā₂) | (iii): the Magnus series of εA₍₃₎ (A₍₃₎ := KᵀP₍₃₎ − Λ, a polynomial in ε), "conjugate to Ā(ε) up to O(ε³)", the coefficients differing by commutators with Ā₀ |
| D3s-11 | cosmetic (README, `adaptive_floquet.py`, D‴.26(c)) | README: "exponentially stable" without "locally", "even powers" without "in slow time"; the module docstring's "eigenvalues are the exponents" (their real parts are minus the rates); D‴.26(c) still explains the parity by a term count | README "locally exponentially stable … in slow time … proved for the linearised loop"; the docstring; D‴.26(c) by the identity Ā[K̃](ε) = Ā[K](−ε) |
| D3t-2 | minor (PAPERS 12, §8, LITERATURE, RESULTS, D‴.31 "Hence", *Numerics* (A)) | "small on the instance and falling with k_e" overstated: relative to DR's bias the parameters' part is +56% at k_e = B = 100 and −8.7% for DR* at k_e = 25; its fall with k_e has a one-dimensional cause (the law's leading term is a line integral along the orbit, which vanishes in ℝ) | "small in absolute terms (≤ 3.9·10⁻⁴ per (Eν²)²), −8.7% to +3.3% of DR's and DR*'s biases and +56% of DR's at k_e = B = 100"; its own response 3.9·10⁻⁴, 1.1·10⁻⁴, 2.9·10⁻⁵, 7.4·10⁻⁶ at k_e = 25–200, ∝ k_e⁻² "on a one-dimensional reference, where the law's leading term averages out — on a curve in ℝⁿ it need not" (§8's Open adds ℝⁿ) |
| D3t-3 | minor (D‴.31 Status, §0, audit, PAPERS 12, [S57]) | the parameters' own errors in variables, listed as new, is the classical second-order bias of a constant-gain stochastic approximation from its own fluctuation (Benveniste–Métivier–Priouret; Dieuleveut–Durmus–Bach for constant-step SGD; Box's curvature bias of nonlinear least squares) | [S57] extended with Dieuleveut–Durmus–Bach (Ann. Statist. 48(3):1348–1382, 2020) and Box (JRSS B 33(2):171–190, 1971), the details checked by a web search; "new" restricted to the quantification and the observation that the desired regressor does not remove it; PAPERS 12 names Richardson–Romberg extrapolation as a candidate remedy |
| D3t-4 | minor (tests, the Status lists, RESULTS) | the round-21 code paths behind the numbers untested: `separation_drift` ×0.85 and without its E[e] term pass, `mean_biases` ignoring S passes; never executed: the exact periodic mean, the fixed-prior branch, `mean_biases` with S and its mean-function term; the separation test's tolerance 35% of σ₀/k_e | four tests: the separation drift against an independent quadrature of its two moments (the Hessian by second differences of f̂, the OU moments typed independently) to 10⁻⁴, the E[e] term ≥ 2% of every class's drift; the transient (D3t-1); the exact periodic mean at k_e = 50 against the filtered-Gram response under both priors (0.17% and 0.04%; the plain Gram 2.6% and 0.6% away), the fixed-prior branch's force against the quasi-static one, `mean_biases` with S = σ₀ against (H + σ₀I)⁻¹ and not H⁻¹; the curvature term c̄ against an independent contraction to 10⁻⁴, Ĵȳ against the network's directional derivative, the mean function's excess over the mean parameters' bias to 1%; the separation test described as a sign test |
| D3t-5 | cosmetic (D‴.31 *Numerics* (B), Status, §0, §7, RESULTS, `exp_eiv_remedy.py`, md B, test) | for the non-commuting rows R and B were the definitions (H + S₀)⁻¹ and (H + S₀)⁻¹S₀, so RH + B = I was an identity of the definitions | `anchored_block`, `anchored_operators`, `anchored_leak`: a prior held at its anchor on some coordinates (fixed there, or a projected prior's active face) and moving elsewhere with a matrix anchor gain K; Part B's five non-commuting priors from their dynamic blocks — fixed amplitudes or centres with η = σ₀ elsewhere, one active face at η = ∞, K = ηP_amp, a random K with K + Kᵀ ⪰ 0 —: the blocks' R and B equal the definitions with S₀ by elimination to 1.9·10⁻¹⁵, RH + B = I to 1.4·10⁻¹⁴, HR + B − I 0.47–0.79; `--only-B`; the test |
| D3t-6 | cosmetic (wording) | the slope-free control's evenness stated as fact; "adaptive (every η > 0)" labelling columns computed at η = ∞; the module docstring's "not excited"; (iv)(c)'s "its slope–parameter term" (DR* has none) and "not a least-squares force" | "numerically even (to k_e = 800; not proved)"; "adaptive (η = ∞; the noise part for every η > 0)" in both mds; "not excited at first order"; "for MR and DR, its slope–parameter term (DR*'s compensation, at x_m, has none)", "not the learning force of a mean tracking error, so no bound in the noise's moments applies to it" |

**Lesson.** A threshold read off a linearisation at a point is the threshold at that point. The separations' drift was evaluated at the
lift, with the noise's own mean error; but the noise also moves the point — the errors-in-variables deviation develops on the slow time
scale, and its tracking error enters the drift at the same order — so along the trajectory the drift grows, and the long-time threshold,
where the deviation is O(1), lies outside the expansion. Agreement at one noise level was a cancellation against the next order; antithetic
pairs at several levels separate the orders. The lemma's review made the same point about asymptotic statements on the other side: a
formula proved for the period average holds as a period average (pointwise the fast jitter oscillates at the floor's own order), and which
term dominates a deviation is decided by computing the terms separately, not by its order.

## Round 23 (adversarial review of Proposition D‴.32, the cohesion leak, and of Corollary D‴.33, the desired regressor and the coherence under noise, 2026-10-20)

**Design.** Two referees, neither modifying the repository. D3v: Corollary D‴.33 on a snapshot of bfdbc0f (`review/REVIEW_round23_D3v.md`;
its issues ID-1–ID-10 are D3v-1–D3v-10 here): an independent instance with jets to order three, the moments re-derived, an independent
real-loop simulator (20 of Part D's 32 runs to all printed digits), the exact second-order theory of a separation's finite-time growth
(the cross term and the variance term), the exact periodic covariance and second-order mean, the second-order transient over 5.1·10⁵ time
units by period maps, eleven targets, (ii)'s bound in 36 cases, eleven code mutants. D3u: Proposition D‴.32 on a snapshot of 64cb25e
(`review/REVIEW_round23_D3u.md`; its issues are D3u-1–D3u-10, applied in the next commit): an independent implementation with analytic
Jacobians, every Part of the numerics rerun from its seeds, twin configurations and dead classes beyond the md's, other targets and
schedules, twelve mutants against the tests.

**Verdicts (D3v).** (i) holds as the mean Jacobian at the lift (the regressor's term is absent; the cross term and the variance term are
not negligible for DR and DR*); (ii)'s bound holds (‖S_k‖ ≤ 0.40 of it in 36 cases), its uniformity in k_e needs B ∝ k_e; (iii) does not
hold at its order; (iv)'s formula is exact in all eleven targets, the whole matrix's sign the instance's and the lift's; every printed
number reproduces, the md regenerates byte-identically.

**Verdicts (D3u).** (i)–(iii) hold, (iii) exactly (the blocks with the mass factor, ker Df̂ the separations where H_L(w) ≻ 0); (ii)'s bias
order is an upper estimate (7.6–8.3 times the actual on the instance); (iv) holds for the linearisation (the margin exact at gauge
transforms, Part A3's split at the predicted rate, Part G's tallies); (v) is not proved — the twins' argument is wrong (a third class within
reach wins: a counterexample on the instance; far twins do not merge in 4·10⁴ slow units) and the dead class's holds only for
well-separated live classes, on a time scale e^{d²/(2R²)}/(σ₀w); (vi) holds (the 45 exponents to 5.2·10⁻¹¹ by a variational monodromy);
every printed number reproduces and the md regenerates byte-identically, but the readings overstate the selection of the count.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3v-1 | major (D‴.33(iii), proof, Status, *Numerics*, §0, §7, §8, RESULTS, LITERATURE, `exp_eiv_remedy.py` Part D, the test) | the stationary point's balance ⟨Ĵᵀ(g − Ĵȳ)⟩ = 0 omits the parameters' part of the second-order force (D‴.31(iv)(c): of the same order, not least-squares) and the drift there omits the cross term; with both the exact stationary mean error has rms 6.70·10⁻⁴ per unit Eν² (DR) against the stated 2.34·10⁻⁴, the drift is indefinite (+3.68·10⁻⁵), and "zero when g lies in the tangent space" and "the drift vanishing at the stationary point" are false | (iii) restated: ⟨Ĵᵀ(g + g_par − Ĵȳ)⟩ + k_e⟨f_{y,par}⟩ = 0, k_eE[e] = (I − P)(g + g_par) − ĴH⁻¹k_e⟨f_{y,par}⟩, the drift ⟨E[e]M_k⟩ + C_k; `second_order_mean(periodic=True)` returns the noise and parameters' parts of the periodic mean and E[ε_mδϑ_k]; Part D's stationary table from the exact periodic mean under the adaptive and the fixed prior, the drift split into the noise part of E[e], its parameters' part, the cross term and the total (DR: rms E[e] 6.70·10⁻⁴ against the stated 2.34·10⁻⁴, its parameters' part 6.32·10⁻⁴ and tangent part 5.80·10⁻⁴; the drift [−4.35, +3.68]·10⁻⁵ in class 1; DR*: 6.44·10⁻⁴; under the fixed prior DR negative definite, DR* positive definite with the threshold std 0.435); the novelty item dropped; a test (the parameters' part more than twice the noise part's, largely in the tangent space, the drift indefinite) |
| D3v-2 | major (D‴.33(iv), §0, §8, *Numerics*, §7, RESULTS, PAPERS 12) | "the noise stiffens the separations" holds near the lift only: along DR's exact second-order transient the drift turns indefinite at t ≈ 2.7·10³–4.1·10³ (nine times the md's window), its largest eigenvalue up to +1.05·10⁻⁴ per unit Eν² (t ≈ 3·10⁴), in the corollary's own mechanism (the noise part alone does the same) | `long_transient`: the second-order transient from the lift carried by exact one-period maps (the covariance through the fundamental matrix, the forcing through the field's second-derivative tensor, the mean by the grid's RK4 maps; the forcing split into its noise and parameters' parts, the cross term from the covariance), 504 grid points per period (0.3% from the exact periodic mean at the stationary point, second order in the grid), checked against `second_order_transient` to 4·10⁻¹⁰ over the first periods; Part D's long-transient table for DR and DR* over 82 000 periods (DR's largest eigenvalue first positive at t = 4103, 2947, 2714, with the noise part alone at 4072, 2721, 2689, at most +1.05·10⁻⁴ at t ≈ 3.1·10⁴, the formal threshold std 1.64; DR*'s positive throughout); (iv) restated "near the lift … along the trajectory the sign is not kept", the mechanism (the fast collective modes relax first, so E[e] passes through the part of g carried by the slow modes, which has no sign relation to ∂²_xf) |
| D3v-3 | minor (*Numerics*, §0, (iv), PAPERS, RESULTS) | the whole matrix's negative definiteness is the instance's: for other targets DR's drift at the lift is indefinite (three separated same-sign bumps +2.9·10⁻⁴); the saturation's term in the centre entry has no sign | (iv): "negative definite for separated kernels of alternating sign well inside the orbit, as on the instance, and indefinite for same-sign, wide or edge kernels", the saturation's term "has no sign", the sign "where that term is dominated"; Part D3: six targets at the lift (the centre entries negative in all six; the whole drift indefinite for three same-sign bumps +2.9·10⁻⁴, wide bumps +5.6·10⁻⁴, edge bumps +2.7·10⁻⁴, two close same-sign bumps +5.0·10⁻⁴) |
| D3v-4 | minor ((ii), §0, Status, *Numerics*, §8, RESULTS, PAPERS, LITERATURE, README, md) | the threshold's uniformity in k_e needs B ∝ k_e: at the fixed B = 100 DR's drift turns positive definite from k_e = 200 and the thresholds fall toward the bound's level; MR's falls like 1/k_e only while k_e ≪ B | (ii): "bounded below uniformly in k_e — Eν² ≥ 2σ₀/(‖∂²_xf‖_∞max_k⟨‖M_k‖⟩) up to (1 + O(1/k_e)) — and constant when B/k_e is fixed", the fixed bandwidth's sign change at k_e ≈ B, MR's threshold levelling off for k_e ≫ B; Part D2: the scans at B = 4k_e and B = 100 for k_e = 25–800 with ‖S_k‖ against (ii)'s bound (at B = 100 DR positive definite from k_e = 200, the threshold std 0.185, 0.138, 0.121 at k_e = 200, 400, 800, DR*'s 0.238 → 0.113, MR's 0.031 → 0.014; ‖S_k‖ ≤ 0.37 of the bound, the floor 0.069) |
| D3v-5 | minor ((i), proof, §0, the thresholds, `separation_drift_series`, Part D) | the cross term and the variance term are negligible against MR's drift, not against DR's and DR*'s (60–2000 times smaller) | (i) states the cross term C_k (zero at the lift, building up within a period) and the variance term E‖P_⊥v₁(t)‖²/t of a separation's exponent (to second order E log‖u(t)‖ = −s_kt + ∫θᵀ(S_k + C_k)θ + E‖P_⊥v₁‖², ε_mM_k symmetric), the proof the expansion; `separation_cross_term`, `separation_drift_series(Eeth=…)`, `second_order_transient` returning E[ε_mδϑ_k], `separation_variance` (the extended covariance with v̇₁ = M_kθ(δe + ν)); DR*'s threshold near the lift with the cross term (the cross term −8% of DR's and −17% of DR*'s class-2 drift over [0, 50), DR*'s threshold std there 0.277) |
| D3v-6 | minor (*Numerics*, §0, RESULTS, Part D, md) | Part D compared DR with the lift's drift (its noise-induced part off by 18–19% at 22–24 s.e., explained by the transient); the start direction was DR's most favourable one; "to 3.6%" measured the whole rate | Part D6: the prediction from the transient along each start direction (the mean error's term, the cross term, the variance term), the noise-induced part measured against it with its sampling error; DR started along its own top eigenvectors of S₁, S₂, S₃ and DR* along its top eigenvector of S₂, four antithetic pairs at std 0.05 and 0.1 (DR along u: the noise part 0.98 ± 0.01 and 0.97 ± 0.01 of the transient's prediction against 0.82 and 0.81 of the lift's; along DR's own top eigenvectors faster than σ₀/k_e by 0.36%, 0.14%, 0.26% at std 0.05, at the prediction to 0.05–1.9 sampling errors; DR* slower than σ₀/k_e by 3.2% and 12.6%, 0.97 and 0.96 of the prediction); "the pair means decay faster than σ₀/k_e; no run grows", the runs slower than σ₀/k_e counted |
| D3v-7 | minor (tests, Status, RESULTS, §7) | DR*'s drift doubled passed every test; the stationary assertion tested the averaged noise part only; nothing tested the transient | three tests: DR*'s drift against an independent quadrature (the Hessian by second differences of f̂, DR*'s mean error typed independently) to 10⁻⁴; the stationary mean error from `second_order_mean(periodic=True)` (the parts' sum, the parameters' part dominant and in the tangent space, the drift indefinite); `long_transient` against `second_order_transient` to 10⁻⁷ over the first periods and DR's drift turning indefinite within 800 periods; the old stationary assertion described as the averaged noise-only residual's |
| D3v-8 | minor (Status, §0, PAPERS 12, LITERATURE) | the classical frame — noise-induced synchronisation of identical systems by a common noise, decided by the transverse Lyapunov exponent at second order (Teramae–Tanaka, Goldobin–Pikovsky, Nakao–Arai–Kawamura) — uncited; two of the four "new" items do not hold | [S61] (the details checked by a web search); the Status, the §0 row, PAPERS 12 and the audit cite it; "new" restricted to the removal of the regressor's term in the K₀-atom identification (a threshold bounded uniformly in k_e) and the instance's quantification |
| D3v-9 | minor (the setting, Status, §8, LITERATURE) | "Proposition D‴.32(vi)" was not yet in the document | Proposition D‴.32 written (commit 64cb25e, round 23's D3u); its (vi) gives σ₀w_k/k_e; the transfer's thresholds (0.1376, 0.0178) kept |
| D3v-10 | cosmetic | "+5.0·10⁻⁵" (DR*'s smallest eigenvalue in class 2) is 5.1·10⁻⁵; G_c undefined; `separation_drift_series` takes no variant; "the noise white over the loop's bandwidth" at B = 4k_e | "+5.1·10⁻⁵"; G_c := G_γc′(γ)² in (iv); the docstring: pass E[ε_m²] = 0 for DR and DR*; the md: "with the noise's bandwidth scaled with the gain (B = 4k_e)" |
| D3u-1 | major ((v), proof, Status, §0, §8, RESULTS, PAPERS, README, the audit) | "the twins' attraction decreases J^coh_m along the curve all the way to their merger … the number of classes is the representation's": both the twins' pull and a third class's are exponentially small, and a third class nearer one twin wins — on the instance opposite-sign twins of bump 3, the light one at the amplitude coordinate of the class of bump 2, form a saddle the flow leaves toward a biased four-class fit; far twins do not merge on the runs' time scale | (v) restated as a rate: for twins with no third class within reach, the fit relaxed quasi-statically, Δα contracts at σ₀ψ(Δα²/(G_αR²))(w_ka′_k + w_la′_l)²/(w_ka′_k² + w_la′_l²) (proved to leading order in σ₀ψ against the fit's relaxation), the merger on the time scale G_αR²e^{Δα²/(2G_αR²)}/(σ₀(w_k + w_l)Δα²); "a third class within reach of one twin can win", "the number of classes is the representation's" dropped; Part I (the twins' contraction 1.2949 against the integrated rate 1.2943; twins 5.28 and 8.19 apart unmoved in 4·10⁴ slow units; the referee's opposite-sign twins not merging in 4·10⁵, the fit error 3.2·10⁻⁴ → 4.5·10⁻³, the light twin drawn toward the class of bump 2; the same under the annealed schedule); a test of the rate (to 1% at 10σ₀) |
| D3u-2 | major (the title, §0, §8, README, PAPERS 13, *Numerics* (C), (H)) | "selects the number of kernels": nothing selection-related is proved; the annealed schedule selects; at a fixed Gaussian scale the broad starts do not end at three classes; the schedule fails on bumps closer than the broad scale and on unrepresentable transit masses; README gives one leak the two kernels' properties | the title "a leak on the unexcited subspace, without a prior"; the Status "not proved: the transit … and with it any selection of the number of kernels: the count is the transit's and the schedule's"; Part J (at a fixed Gaussian scale 4–6 classes at R = 1, three biased at R = 2; the schedules 1.5 → 1, 2 → 0.5, 3 → 1 at three classes) and Part K (three equal bumps three classes, the Gaussian kernel's bias; four alternating four, a slow transient — Part L continues both; three heavy (Σ|A_k|/a_max = 0.95) 3–5 biased; two bumps at ±0.45 merged into one stable class); README and PAPERS 13 state which kernel has which property and what the schedule did where; the replicator of Proposition D‴.34 removes the mass obstruction |
| D3u-3 | minor ((v), proof, Status, *Numerics* (A4)) | the dead class's merger is shown only at R = 2, where the live classes overlap; at R = 1, where the hypothesis holds, it does not move; the mixture's weights carry the live classes' amplitude coordinates | (v): the weights w_le^{−α_l²/(2G_αR²)}, the time scale R²e^{d²/(2R²)}/(σ₀wd²), "its merger not proved"; Part I (R = 1: unmoved in 4·10⁴ slow units, the time scale 1.4·10⁷; R = 1.5: joined at s ≈ 2.8·10⁴ against 1.8·10⁴) |
| D3u-4 | minor (*Numerics* (B), §7, RESULTS) | the m = 15 maximum 7.1·10⁻⁵ and the compact kernel's fourth class come from seed 3's sample (8, 4, 3), on the boundary of W, which no leak repairs | (B) reports the representable samples' range and seed 3's sample separately (the representable samples 1.9·10⁻¹⁰–5.4·10⁻⁸ under both kernels and the adaptive prior; seed 3 at 7.0·10⁻⁵–7.1·10⁻⁵) |
| D3u-5 | minor (tests, Status) | consistent mutants of the kernel's range or scale pass all tests; the twins test is loose; `clusters` untested | three tests: the kernels' range and scale (under the compact kernel a pair at 0.99R attracted and at 1.01R not, the Gaussian pair force at R e^{−1/2} of its linear part, both in the G⁻¹ metric); the twins' contraction against the rate formula to 1%; `clusters` on known configurations |
| D3u-6 | cosmetic (*Numerics* (G), §0, §7, RESULTS) | the −5·10⁻⁵ criterion and the shedding unexplained; σ₀w_k², not σ₀, sets the first split | (iv): "a class can split iff σ₀w_k² < q_kλ⁺_k", the (n_k − 1)-fold degeneracy and the kernel's concavity (ω″ = ψ′/2 < 0); Part G prints the growth estimate −ln(250)/10⁵ = −5.5·10⁻⁵ and the cohesion energy of the split patterns (for a class of ten at Σ‖u‖² = 0.5 the cohesion energy of shedding one particle 0.43 of the quadratic part against 0.74 for a halving, compact; 0.87 against 0.95, Gaussian) |
| D3u-7 | cosmetic (*Numerics* (A), (D), §0, §7, RESULTS) | finite-difference precision quoted as agreement | "(the finite differences' error)" at the spectra, the margins and the monodromy's exponents; the Gaussian kernel's 4.6·10⁻⁹ named its cross-class terms |
| D3u-8 | cosmetic ((i), proof, §0, §8, audit) | "the fit as the external field": L̄ is quadratic in μ, an external potential plus a second interaction; the metric is G⁻¹'s | (i): "the Wasserstein gradient flow (in the metric G⁻¹) … an aggregation equation in the external potential −⟨fΦ̃⟩_ρ with two interactions, the fit's ⟨Φ̃Φ̃′⟩_ρ and the cohesion W"; the proof writes L̄(μ) out |
| D3u-9 | cosmetic ((vi), *Numerics* (E), §0, §7, RESULTS) | τ/(σ₀w_k) is a variance, the compared numbers rms; the (1 − 1/n_k) factor missing | (vi): "the particles' stationary variance about the class mean is (1 − 1/n_k)τ/(σ₀w_k)"; Part E prints the second-half time averages against √((1 − 1/n_k)τ/(σ₀w_k)) = 0.0635 (0.068 compact, 0.063 Gaussian, the fixed prior 0.036 against 0.0366) |
| D3u-10 | minor (the audit, [S60], §0, PAPERS 13) | precedents omitted: polarized consensus-based dynamics, clustering in self-attention dynamics, Fukumizu–Amari's critical lines, convex clustering and sum-of-norms identification, deterministic annealing, stabilized forgetting | [S60] extended with Bungert–Roith–Wacker, Geshkovski–Letrouit–Polyanskiy–Rigollet, Fukumizu–Amari, Hocking–Joulin–Bach–Vert, Ohlsson–Ljung–Boyd, Rose (the details checked by a web search); the audit and the §0 row cite them; "new" restricted to the combination in mean-field adaptive control: the exact blocks at clustered configurations with the mass factor, the gauge-scaled margin, the twins' rate, the real loop's exponents |

**Lesson.** A sign found at a point is the sign at that point. Round 22 learned it for the threshold of the measured regressor; here it held
for the desired regressor's sign and for the stationary point: the collective modes relax over 10³–10⁵ time units, fast ones first, and on
the way the mean error passes through combinations that the lift's formula does not see; and a balance written for the noise's own force
omits the parameters' part, which at the stationary point dominates. Exact one-period maps make the whole second-order transient cheap —
about a minute for 5·10⁵ time units — so the trajectory, not only its end points, can be computed; and a prediction compared with the whole
rate, 95% of it the noise-free decay, hides a 20% miss in the part that is tested. A count found is not a count selected. The cohesion
leak's stationary set holds every clustered critical point, so what it ends at is the transit's and the schedule's; the selection claimed
was a property of one target and one schedule, and two more targets and a fixed scale undid it. An energy argument along a neutral curve
must look at every class within reach, not only at the pair that defines the curve. And a numerical check placed outside the hypothesis's
regime (the dead class at R = 2, where the live classes overlap) checks something else.

## Round 24 (adversarial review of Proposition D‴.34, the replicator weight law, and of Proposition D‴.32's round-23 revisions, 2026-10-20)

**Design.** One referee, not modifying the repository. D3w: Proposition D‴.34 and the round-23 revisions of D‴.32 on a snapshot of 1ddff8c
(`review/REVIEW_round24_D3w.md`; its issues ID-1–ID-9 are D3w-1–D3w-9 here): an independent implementation of the weighted loop, its energy
and the real loop; the gradient-flow identity in the Euclidean–Fisher–Rao metric and the √β form; the relative-entropy identity off its
stated hypotheses; the class-level stationary points under both kernels with the masses held and freed; Parts A, B (seed 1), C (seed 0),
E and F (three heavy bumps, seed 0) reproduced to the printed digits; a continuation of Part B; the zero modes against the predicted
vectors; nine code mutants against the tests; D‴.32's twins' rate, time scales and Parts I–L.

**Verdicts (D3w).** (i) holds (the identity to 9·10⁻¹⁰, the √β form with the step γ/4); (ii) holds (V̇ = −k_e|e|² to 1.1·10⁻¹⁶), with
other hypotheses than stated; (iii) holds at clustered configurations under the finite range, not under the Gaussian kernel, where the
stationarity condition's factor β_i and the condition u_k = λ fix the masses; (iv) holds for bounded trajectories; (v) holds exactly (the 15
zero modes are the 12 redistributions and 3 gauges, 14 on Σβ = m); every reproduced number matches; the md regenerates byte-identically;
D‴.32's round-23 fixes are applied correctly where checked.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3w-1 | major ((iii), Status, §0, §8, PAPERS, README, the audit) | "the stationary points are D‴.32's for every vector of masses", "the stationary set unchanged" hold only at clustered configurations under the finite range: the positions' condition carries β_i on the fit and not on the cohesion, and under the Gaussian kernel u_k = λ fixes the masses (five bumps: the masses move and the fit error 2.63·10⁻⁴; three equal bumps: (0.267, 0.267, 0.467) → (0.22, 0.22, 0.56), 9.17·10⁻⁵) | (iii) restated with the factor β_i; "there the stationary points are D‴.32's" restricted to the finite range at classes R apart; under the Gaussian kernel the replicator selects the masses and moves the bias; Part G: class-level stationary points under both kernels with the masses held and freed (B's seed 1: 3.85·10⁻⁴ with the masses held, 2.63·10⁻⁴ freed; three equal bumps 9.39·10⁻⁵ and 9.17·10⁻⁵, the masses toward (0.567, 0.217, 0.217); the compact kernel exact either way); the §0 row, the Status, the §8 bullet, PAPERS 13, README and the audit qualified |
| D3w-2 | minor (*Numerics* (B), §0, RESULTS) | "still falling at the end (by 3–4% per 1000 slow units)": the rates are 2–4%, and the fit levels off near the Gaussian stationary point; under the compact kernel a slow mode of the nearly saturated classes holds it | (B) the rates per seed and Part G's continuation of seed 1 (2.80·10⁻⁴ after 10⁵ more slow units under the Gaussian kernel, 2.85·10⁻⁴ after 4·10⁴ under the compact one) |
| D3w-3 | minor ((iv), proof, §0, §8, PAPERS, README, RESULTS, *Numerics* (A)) | "the exact fit is reached iff representable" unproved in the "if" direction (boundedness numerical); f in the span and one particle per atom assumed; the unbounded case needs the ω-limit argument; classes may lose their mass instead of saturating; Part A's 0.071 creeps to the ℓ¹ infimum 0.0684 | (iv): "bounded trajectories reach the exact fit iff …", the hypotheses stated, the ω-limit argument in the proof, "each class saturating or losing its mass", the ℓ¹ bound; Part G prints the infima (0.0684 with the masses free, 0.387 and 0.494 at the held masses); the summaries "for bounded trajectories"; the test asserts the fit above the infimum |
| D3w-4 | cosmetic ((i)) | the WFR remark omits that WFR divides the cohesion by β_i and that J is not a functional of μ_{Θ,β}; r_i clashes with the residual; the output weights are r_i²/m; the step γ/4 | (i) rewritten: ϱ_i := √β_i, the step γ/4, ϱ_i²/m the output weights, the WFR metric's two differences |
| D3w-5 | cosmetic ((ii), proof, (𝔉)) | (ii)'s hypotheses are not the needed ones (Σβ* = m and the pointwise representation suffice; V ≥ 0 needs Σβ = m, invariant, not attracting); (𝔉)'s λ and the code's agree on Σβ = m only; β → β* needs persistent excitation | (𝔉) with λ := Σβu/Σβ (the code's: Σβ conserved from every start); (ii) restated with its hypotheses, V ≥ 0 by the relative entropy on Σβ = m, "β → β* is not claimed" |
| D3w-6 | minor (tests, RESULTS, §7) | every drift mutant caught only by the gradient identity (the code against its own J); the fixed-position test passes with λ = 0; the zero-mode test passes every mutant | the fixed-position test asserts Σβ = m and the fit above the ℓ¹ infimum; the zero-mode test checks the predicted null vectors (the redistributions, one gauge per class) and a separation's −σ₀w_k with γ > 0; a new test checks the weights' drift against an independent quadrature of the fitness |
| D3w-7 | minor ([S62], the audit, Status) | Kondratyev–Monsaingeon–Vorotnikov and Chizat–Peyré–Schmitzer–Vialard (WFR), Chizat–Bach (the √β lift), Shahshahani's metric and EM's multiplicative weights omitted; (ii) is the negative-entropy instance of Boffi–Slotine's mirror-descent laws | [S62] extended (KMV, CPSV and Shahshahani checked by a web search); the Status and the audit name (ii) as that instance |
| D3w-8 | cosmetic ((iii), *Numerics* (C)) | "whose amplitudes are not saturated" vacuous; the bias at fixed positions is permanent, not "for times exponential"; the rate exact only with the class's particles coincident and the others critical; (C)'s "the same splits at the same times" differs at 3 of 51 records | (iii) rewritten; (C) "at nearly the same times (3 of 51 records differ by one class)" |
| D3w-9 | cosmetic (D‴.32 *Numerics* (A); D‴.34 (D)) | the dead class's mass is 2/15 (4/15 is the nearest live class's, the w of the time scale); (D) compared with six seeds instead of seeds 0–2 | "(amplitude 0, mass 2/15, at 5.19 from the nearest live class, of mass 4/15)"; (D) against seeds 0–2 (1.9·10⁻⁷–9.0·10⁻⁷) |

**Lesson.** A stationary set is a property of the whole law. The weights entered the positions' drift through the fit and not through the
cohesion, so the clean reduction held only where the cohesion vanishes — and every broad run used the kernel under which it does not.
And a claim of reaching a target needs the trajectory bounded; a theorem about bounded trajectories says nothing about the ones that run
to the boundary, where the saturated classes live.

## Round 25 (adversarial review of Proposition D‴.35, the split as an S_n-symmetric transcritical bifurcation, 2026-10-20)

**Design.** One referee, not modifying the repository. D3x: Proposition D‴.35 on a snapshot of 44fe0f9 (`review/REVIEW_round25_D3x.md`; its
issues ID-1–ID-11 are D3x-1–D3x-11 here): an independent instance (round 23's `inst.py`) with the quadratic coefficient re-derived and
measured, the slaving at third order, the full Jacobian's spectrum at the lift, the branches including the S₅ × S₅ states, all twelve Part C
runs and the mass-free run of seed 0 reproduced from their seeds, the mass-free schedule extended to σ₀/256, the margin after a departure
decomposed, mutants against the tests.

**Verdicts (D3x).** (i) holds to second order (the coefficient c₃/2 to 3·10⁻⁸ on a scale of 5·10⁻⁴; no even part in either kernel or the
mass-free leak), the centre-manifold regime not covering σ₀/16; (ii) holds for 1 ≤ q < n/2 and both signs of μ, the S_{n/2} × S_{n/2} states
missing; (iii) the order and the escape time hold (the quadratic equation's particle in 12 of 12 runs, its time to 1.2%; the linear one's in
6 of 12, to 22%), "alone" does not; (iv) the margin's accounting is wrong; (v) the blocks hold, "sheds once per crossing" is refuted; every
reproduced number matches; the md regenerates byte-identically.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3x-1 | major ((v), *Numerics* (D), Status, §0, §7, §8, RESULTS, PAPERS, README, the audit) | the mass-free leak's thresholds are a factor w_k lower, so Part G's schedule reaches them later; extended to σ₀/256 (seed 0) it ends at 9 classes and 0.0028, as the cohesion leak does: the 4 classes were where the schedule stopped | (v) restated: the thresholds σ = λ⁺_k, the cascade not stopped; Part F extends the mass-free schedule to σ₀/256 for seeds 0–2 (9, 10, 9 classes at σ₀/128 and 9, 8, 9 at σ₀/256, the fit error 0.0028, in all three seeds); "sheds once per crossing", "the target's count" and "the mass factor drives the cascade" removed everywhere |
| D3x-2 | major ((iv), proof, §0, §8, PAPERS, README) | the margin after a departure omits the curvature's gauge rise (λ⁺/(n − 1), 1.35–2.4 times the stiffness's drop) and the condition is off by μ; the residual's change dominates both | (iv): the margin σ₀(w_k − k/m) − λ⁺(r′, ϑ′) with the gauge factor n/(n − k), "stops only if the residual's change lowers λ⁺ by more than μ + kσ₀/m plus that rise", the residual's change dominating on the instance; the proof cites D‴.28's scaling |
| D3x-3 | major ((iii), *Numerics* (C), §8, PAPERS, README) | "alone" fails: at σ₀/16 bound pairs leave in 3 of 6 runs (the top two ζ₀ within 1–11%); the 0.1-from-the-median criterion is far inside the range | (iii): "alone when its lead is large, with its neighbour as a bound pair when nearly tied"; Part C reports the pieces beyond the range at the end (at σ₀/16 a single particle in 3 of 6 runs and a bound pair in 3; at μ = 5·10⁻⁵ a single particle in 3, none yet in 3); "sheds the extreme particle, alone or with a near-tied neighbour" |
| D3x-4 | minor (*Numerics* (C)) | the "+ side" first excursions are not the common drift (it cancels against the median) but linear growth below 2μ/\|c₃\| = 0.22; those particles do not leave | (C) rewritten; Part C prints the linear and flipped-c₃ predictions |
| D3x-5 | minor ((C), §0, §8, RESULTS) | "first on the side sign(c₃)v the top-ranked, 12 of 12" is order preservation (the linear law gives it too); the first excursion is on that side in 10 of 12 | the discriminating statistics reported: the quadratic equation's particle in 12 of 12 (linear 6, flipped c₃ 2), its time to 1.1% (linear 22%); the side in 10 of 12 |
| D3x-6 | minor ((ii)) | for even n the S_{n/2} × S_{n/2} states are cubic and missing; "leaves every neighbourhood" too strong | (ii): the cubic states named; "a fixed neighbourhood from almost every initial condition" |
| D3x-7 | minor ((i), (iii), Status) | at σ₀/16 three stable modes are slower than μ; the third order is mostly the class mean's slaving and the residual's own cubic, nearly cancelling, the kernel's concavity the smaller term | (i): "for \|μ\| small against the other modes' rates", the third-order terms named; the Status: at σ₀/16 the agreement reflects the weak coupling |
| D3x-8 | minor (tests) | test 3 ran at a negative cohesion strength (a repulsion); `reduced_crossing` untested; `MassFree` untested | test 3 at μ = 10⁻⁴ with σ > 0 and the other classes' margins positive asserted; a test of `reduced_crossing` against the linear closed form and an independent integration; a test of the mass-free leak's blocks (−σ on a separation, 0 on a class mean) |
| D3x-9 | minor ((E), RESULTS, §8) | "a faster mass transfer shortens the cascade": γ = 1 (5, 6, 8) and γ = 0.1 (6, 6, 7) have the same mean; the seeds are not paired with Part G's | (E): "with the replicator the cascade ends at 5–8 classes against 8–10 without; the rate makes no difference in three seeds", the seeds' jitter unpaired with Part G's said |
| D3x-10 | minor (Status, §0, [S63], the audit) | order preservation is equivariance and the escape time the Riccati solution; Stewart–Elmhirst–Cohen, Brea et al., Şimşek et al. and Rose missing | "new" restricted to the identification of the coefficient, the margin after a departure and the quantification; [S63] extended (Stewart–Elmhirst–Cohen checked by the referee's web search; Brea et al., Şimşek et al. from memory); Rose in the audit |
| D3x-11 | cosmetic | σ/ψ̄_i should be σ/(σ₀ψ̄_i); Δs undefined; "side" missing; the closed form is not the time against the median; (D)'s "delayed" sentence | fixed in (v), (iii), (C), (D) |

**Lesson.** A count that comes out right at the end of a schedule is a statement about the schedule. The mass-free leak "stopped" the
cascade only because its thresholds lay one stage lower; continued, it shed like the original. And a statistic that any order-preserving
law reproduces tests the order, not the term it was meant to test: the discriminating number was the time, and the linear law missed it by
22%.

## Round 26 (adversarial review of Proposition D‴.36, the gauge under dither, 2026-10-20)

**Design.** One referee, not modifying the repository. D3y: Proposition D‴.36 on a snapshot of 15dcd79 (`review/REVIEW_round26_D3y.md`; its
issues ID-1–ID-12 are D3y-1–D3y-12 here): an independent model matching the drift to 10⁻¹⁷, the spectrum and the projector, the limit map by
long LSODA integrations and its drift term, the full linear Ornstein–Uhlenbeck problem with the leak, batched simulations with 256 seeds
(Part B, at τ = 10⁻⁴ and 2.5·10⁻⁵) and 64 seeds per ε (Part D), mutants against the tests.

**Verdicts (D3y).** (i) correct for one particle per class (M normally attracting, the slowest normal rate independent of γ); (ii) correct —
the metric, the spectral projector equal to the g-orthogonal one to 2·10⁻⁹, Table A to every digit, the closed form and the γ-limits;
the uncomputed drift (−2.055, +1.580, +0.475)·10⁻⁷ explains Part B's means; (iii) mostly correct — the leak a Fisher–Rao gradient, κ and the
Lyapunov equation right (0.5% and 4% against the full linear covariance), the O(τ/ε) shift confirmed; the numbers reproduce bit for bit.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3y-1 | major ((i), §0 row) | "(or with its classes collapsed)": with n_k > 1 the within-class redistributions are zero modes too (dimension m − 1), and T_α needs n_k | restricted to one particle per class; the n_k > 1 case stated |
| D3y-2 | major ((iii)) | where w is not representable the bias is O(ε/γ), not O(ε) (the minimiser depends on ε/γ only) | "a bias of order ε/γ" |
| D3y-3 | major (*Numerics* (D), §0, §8, RESULTS, PAPERS) | at ε = τ the 4-seed agreement is chance: 64 seeds give spreads 14–23% above the prediction, and the ε^{−1/2} ratio 3.7–4.0 | (D): only ε = 10τ informative (to 3%), ε = τ unresolved at 4 seeds, the 16-seed test in Proposition D‴.37 (G); "to 11%" removed |
| D3y-4 | minor ((D), (B)) | the drift is computable and predicts Part B's means; "0.008 vs 0.0005, as O(τ/ε) says" is 4-seed noise | the drift computed in Proposition D‴.37; the claim removed |
| D3y-5 | minor ((B)) | the slope used the mean square, which includes the drift's square | the displacement's variance: ratios 0.95, 1.36, 0.90 (`exp_gauge_dither.py`) |
| D3y-6 | minor (Status) | the time-scale separation is thin (slowest normal rate ≈ 3τ; at ε = 10⁻³ the leak's stiffness exceeds it) | said in the Status |
| D3y-7 | minor (proof) | "BBᵀ = Σ" should be BBᵀ = G | fixed |
| D3y-8 | minor ((ii), Status) | for K = 3 only the saturating class's row and column of D vanish; the α-diffusion is exactly 2τG_α; the unbounded excursion is outside Katzenberger's theorem | (ii) and the Status restated |
| D3y-9 | minor (tests) | nothing ties the metric to the code's drift or the noise to τ | two tests: D against the projector of the code's drift Jacobian; `run_dithered`'s noise 2τG dt |
| D3y-10 | minor ([S64], New) | Parsons–Rogers' first author T. L. Parsons, J. Phys. A 2017; Constable–McKane–Rogers missing; "in the flow's metric" computed, not new | [S64] fixed and extended; the audit's wording |
| D3y-11 | minor ((D)) | "amplitudes back inside" depends on the number of seeds | said |
| D3y-12 | cosmetic | 0.80 vs 0.8033; (C)'s comparison without error bar; (i)'s proof circular | fixed |

**Lesson.** Four seeds can agree with a prediction to 10% by chance when their own scatter is 20%: an agreement is only as good as the
sample's error bar, and the referee's 64 seeds turned the best-looking number into the one that failed.

## Round 27 (adversarial review of Proposition D‴.37, the dither's drift along the gauge, 2026-10-20)

**Design.** One referee, not modifying the repository. D3z: Proposition D‴.37 on a snapshot of 391b29f (`review/REVIEW_round27_D3z.md`;
ID-1–ID-10 are D3z-1–D3z-10 here): the identities re-derived and implemented without the eigenbasis (a Lyapunov solve and the group
inverse), an independent brute-force limit map (BDF, Richardson) at two new points, one inside the layer; the slow limit to γ = 10⁻⁷; the
split invariance under compact and Gaussian kernels; a 2556-point independent zero search; the full loop's records; mutants.

**Verdicts (D3z).** (i) correct (≤ 6·10⁻⁶ against the independent implementation, 5.9·10⁻⁵ against the brute force inside the layer);
(ii) correct, with the closed form ℓ = 2Σ log(W_k − W_min,k²/W_k) + const, the O(γ) convergence and D = O(γ²τ) confirmed; (iii) the
divergence, the rates, δ* and the inside law right, "stops it / the masses settle there" wrong; (iv) the formula fine, its validity
generous; Part E's seven zeros confirmed; the md regenerates byte-identically.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3z-1 | major ((iii), title, §0, §8, RESULTS, PAPERS, README, audit) | inside the layer α diffuses without drift, so the zero-flux density in W ∝ 1/(W − W_min) is not normalisable; the full loop's excursions grow and 11–12% of the last quarter is beyond the saturation | (iii) restated: the layer slows but does not stop; E's nodes "attractors of the mean drift, not settling states"; H reports the excursions by quarters and the fraction beyond the saturation (`exp_gauge_drift.py`, recomputed from the stored records) |
| D3z-2 | major (Status, audit) | Li–Wang–Arora's Corollary 5.1 already gives the log-pseudo-determinant drift for isotropic noise; Fixman, Bornemann–Schütte, Damian–Ma–Lee, HaoChen et al., noisy replicators, Rose missing | "New:" narrowed; the audit names them (not seen) |
| D3z-3 | minor ((A)) | the quoted differences were stale (2.2·10⁻⁵, 1.1·10⁻⁴ after the steps' refinement) | fixed |
| D3z-4 | minor ((H)) | "means within errors at all three times" false (1.55, 2.23 s.e.); the SDE cannot test the saturation | restated |
| D3z-5 | minor ((iv), (G)) | W_2 above κ⁻¹b at every ε (2.8, 2.1, 1.7 s.e.; the same seeds); the window τ ≪ ε ≪ 3·10⁻⁴ barely met | (iv) and (G) restated |
| D3z-6 | minor ((ii)) | the closed form of ℓ; "converge toward its minimisers" empty (no minimiser) | the closed form stated, the sentence replaced; a test |
| D3z-7 | minor ((ii)) | the Gibbs law is the local well's quasi-stationary law | said |
| D3z-8 | minor (tests) | no test of the Itô term; none at γ = 10⁻² | a nonlinear-observable test; the instance against the brute force |
| D3z-9 | minor ((iii)) | within ≈ 0.1δ* the transverse mass fluctuations equal the distance | said |
| D3z-10 | cosmetic | 5.35 → 5.34; "10⁴ times" → 10³; 0.57 still decreasing | fixed |

**Lesson.** A zero of the drift is not a resting place when the diffusion vanishes there too: in the coordinate where the noise is
regular the "restoring" drift was only an Itô term, and the masses kept going. The drift field alone does not decide recurrence.

## Round 28 (adversarial review of Proposition D‴.38, the Fixman correction, 2026-10-26)

**Design.** One referee, not modifying the repository. D3aa: Proposition D‴.38 at 2b56ae1 (`review/REVIEW_round28_D3aa.md`): the closed
form ℓ against an independent Gauss–Newton Hessian at six points; the conditioning along M; the second-order stationary mean of the pinned
state (the linearised limit SDE with the leak's projected force); the floor's cap and the zero-flux density across it; the finite-γ
correction drift; equipartition with ⟨e²⟩; the reproducibility of every part by the current code; 25 mutants.

**Verdicts (D3aa).** (i)'s cancellation in the slow limit and the closed form right (log det − ℓ constant to 10⁻⁸); (ii) right, with
"best-conditioned" wrong (W♯ is D-optimal); (iii) the exponents right in Katzenberger's limit, "normalisable for every κ_F > 0" false for
the implemented (floored) correction; (iv) the pinned state's first order incomplete — the paired excess is the second-order fluctuation
term; Part D not reproducible by the current code (run before the floor); 11 of 25 mutants surviving.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3aa-1 | major ((iv), Part G) | the pinned state omits the fluctuation term −J⁻¹½∂²f[C], of the same order through the leak's curved projected force; Part G's "unexplained" excess is that term | the second-order stationary mean in (iv) and `exp_gauge_fixman.second_order`; the paired differences within 0.93 s.e. of it, D‴.37's Part G within 1.7 s.e. |
| D3aa-2 | major ((iii), title) | the floored correction's capped c_k makes the density p ∝ δ⁻¹ below δ_cap = 0.00126: not normalisable | (iii)'s law stated for the unfloored correction and the cap δ_cap with it (Part H: the exponent across the cap falls back toward −1, as predicted to 10⁻³); the title no longer claims the saturation non-attracting |
| D3aa-3 | minor | W♯ maximises det ∇²_ΘJ (D-optimal), not the conditioning | "the D-optimal exact fit" throughout; the test renamed |
| D3aa-4 | minor | "no excursion into the saturation" false at κ_F = 1 (|α_3| to 13.5) | (D), (E): the excursions along M at κ_F = 1 stated, none beyond a saturation; |α| ≤ 8.1 at κ_F = 2 |
| D3aa-5 | minor | κ_F = 1's "within 1.3 s.e. of W♯" tests nothing (the zero's relaxation time 2·10⁶) | compared with the deterministic mean path; κ_F = 2 within 1.0–1.1 s.e. |
| D3aa-6 | minor | the over-compensation from γ ≈ 2·10⁻³; the finite-γ correction drift P·F_c = (TᵀgT)⁻¹∇_WΦ missing | (i) and (A): the formula, the reversal from γ ≈ 2·10⁻³, the effective gain |
| D3aa-7 | minor | Part C tests only the τ → 0 algebra | (C) called Katzenberger's limit; the finite-τ claim stated as untested there |
| D3aa-8 | minor | equipartition with the mean of e, not the rms | ⟨e²⟩, τ̂ = m⟨F²⟩⟨e²⟩/p, p_eff reported (8.92–9.04 corrected, 8.06 without) |
| D3aa-9 | minor | Part D run before the floor came in | Part D re-run with the current code; the run cache keyed by the code's version (`runcache.VERSION`) |
| D3aa-10 | minor | 11 of 25 mutants survive (the floor, the tracked average, the analysis) | tests of the floor, the tracked average, κ's factor m, the analysis on synthetic records |
| D3aa-11 | minor | forward references carry unproved claims (⟨u_k⟩ = −τ⟨g′_k⟩ "exactly") | the identity by one integration by parts, in the positions' quasi-stationary law (the global Gibbs law has none) |
| D3aa-12 | minor | adaptive biasing force, well-tempered metadynamics, max-det/analytic centres, Helfand's stiff-limit term missing; "Fixman" needs a qualification | [S66] extended (search-verified or flagged); the audit's "New" narrowed |
| D3aa-13 | cosmetic | 0.0045 → 0.0046; τ̂ 9.4·10⁻⁵ was D's; B's eigenvalues 16–19γ²τ and 29–30γ²τ; the md header | fixed |

**Lesson.** A stationary mean has to be predicted to the order at which the fluctuations act: when the restoring force is curved on the
scale of the distance to a boundary, its fluctuation term is of the same order as the drift balance, and a paired design resolves exactly
that term. A limit statement must be made for the object the code implements, cap included; and a name ("best-conditioned") is checked
by computing the quantity it names.

## Round 29 (adversarial review of Propositions D‴.41 and D‴.42, measurement noise on the gauge and several particles per class, 2026-10-26)

**Design.** One referee, not modifying the repository. D3ab: the drafts of Propositions D‴.41 and D‴.42 at 74e955e
(`review/REVIEW_round29_D3ab.md`): the per-step expansion of the correlated-draw selection term by term; the scalar and replicator
selection against the per-step law at λτ_c up to 4; node covariances other than white; the coordinates of the update; D‴.42 under masked
dither and measurement noise per component; `limit_drift` at degenerate eigenvalues; 21 mutants.

**Verdicts (D3ab).** D‴.41: (i)–(iii) right at leading order; (iv)'s per-step-variance law false beyond λ_maxτ_c ≪ 1 (the cross term
h²∂²Φ[F, ⟨Uξ⟩] dropped); the heteroscedastic and orbit-correlated caveats O(γ), not "needs homoscedasticity". D‴.42: the frozen split
and the per-class law false for a non-isotropic dither; the measurement-noise residual 0.7–9.6% per component at unequal splits.
`limit_drift` silently wrong at degenerate eigenvalues; 11 of 21 mutants surviving.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ab-1 | major ((iv)) | the cross term dropped; the per-step law only for hλ_max ≪ 1 − ρ | the linear-response filter f_il with the cross term (vanishing in mean by the antisymmetry), `label_noise.lr_factor`, `filtered_covariance`; Part C's prediction column, D–F's ρ to 0.999 |
| D3ab-2 | major (D‴.42 (ii), (iii)) | the frozen split and the per-class law fail for masked dither | (ii), (iii) restricted to the isotropic dither, the anisotropic mechanism (three couplings, a non-product law); Part F (masks) |
| D3ab-3 | major (D‴.42 (iv)) | the residual normalised by the fractions' drift; per component up to 9.6%, not vanishing as γ → 0 | Part D per component; the pieces corrected |
| D3ab-4 | major | 11 of 21 mutants survive (Heun's shared draw, the AR(1) normalisation, the weights) | tests of the draws, a Heun step, the weighted and covariance Wong–Zakai terms, the classes' masses form |
| D3ab-5 | major | `limit_drift` garbage at degenerate eigenvalues | the Schur path (spectral projection without eigenvectors), automatic at cond(V) > 10⁸; n = 3 equal-split rows; a test |
| D3ab-6–7 | minor | the node-covariance caveats' γ-scaling; "any spectrum" overstated | (ii): white on the orbit's resolution in the hypothesis, the O(γ) residual and its second mechanism; Part A2 |
| D3ab-8 | minor | the Itô selection depends on the update's coordinates | (iii): the a-coordinate factor (1 − q)(1 + 3q) and its target |
| D3ab-9–10 | minor | D‴.42's Itô potential and the equal split's m-independence | the full trace; the invariant-subspace argument |
| D3ab-11 | minor | literature missing (Ziyin et al., Pillaud-Vivien et al., Kupferman–Pavliotis–Stuart, Kühn–Rosenow, Woodworth et al.) | [S69] and the audit, the novelty narrowed |
| D3ab-12–18 | minor/cosmetic | W_ℓ is the slow limit's target; numerics only where trivial; the saturation's qualification; numbers, notation; Heun's predictor unnormalised | fixed (the γ = 10⁻² zero; ρ to 0.999; (iii), (v); the predictor renormalised) |

**Lesson.** Each major failure promoted a leading-order identity to a universal quantifier without testing it where its assumptions are
weakest — a term dropped as O(h²|F|²) that is O(h²|F||ξ|), a Gibbs free energy borrowed for a non-isotropic noise, a caveat resting on
one profile at one γ; and the tests pinned only the degenerate reductions of the code paths that carry the novelty.

## Round 30 (adversarial review of the revised Propositions D‴.41 and D‴.42, 2026-10-27)

**Design.** One referee, not modifying the repository. D3ac: the drafts revised after round 29 and the code at 9143f28
(`review/REVIEW_round30_D3ac.md`): the per-step expansion's completeness at O(v); the filter against the full loop at ρ = 0.9, 0.95 and
at hλ_max = 1.56 (8 seeds); the scalar model to hλ = 1.78; the second mechanism's γ-scaling; the slow-noise selection's zero set and
basins (multi-start, flows), smooth against AR(1) noise; D‴.42's equal-split law in the finite differences' step; the Schur path; 45
mutants.

**Verdicts (D3ac).** The filter complete at O(v) and confirmed everywhere tested (hλ_max = 1.56 included, with a path-averaged prediction
and enough seeds); the second mechanism right at leading order. Overstated: (iv)(e) and Part G — the slow-noise selection is multistable
(three stable zeros, W♯ an unstable node) and at Part C's γ its targets are interior; (c) and (e) hold for Markov (AR(1), block) noise
only — a smooth noise has f ≈ (λτ_c)⁻² and a bounded Ψ; D‴.42's "internal differences' contribution" is the finite differences'
O(h2²) error. 38 of 45 mutants killed.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ac-1 | major ((iv)(e), Part G) | the selection multistable; "the target moves continuously" one branch of several | Part G: the whole zero set (fsolve and flows from 21 starts), stability, the basins of W₀; (e) restricted, the multistability stated ((c)'s concavity) |
| D3ac-2 | major ((c), (e)) | the log-determinant and the dither identity AR(1)-specific | (c): the Markov noise's; f = ∫ŝω²/(λ² + ω²) decreasing for every colour, Ψ concave, a smooth noise's bounded; `corr_sums`' Gaussian correlations; G2's smooth row |
| D3ac-3 | major (D‴.42 (ii)) | the equal split's "internal differences' contribution" a truncation error | Part F Richardson-extrapolated in h2 with its resolution (ratios 1.0000–1.0004 at γ ≥ 10⁻⁵); the claim removed, a proof open |
| D3ac-4 | minor | (e)'s limit only as h → 0 (2/(2 − hλ)) and away from the saturation | the finite-step limit X_h − BX_hBᵀ = 2hS, proved and in G2 (7·10⁻⁴); λ_soft·τ_c reported along the branch |
| D3ac-5 | minor | G1's guard stopped the branch | the guard 1 + 10⁻⁶ and a finer τ_c grid (to 10⁵) |
| D3ac-6 | minor | (e)'s effect ∝ vh/τ_c², unreachable | stated |
| D3ac-7 | minor | Heun's zero only as h → 0 | (iv)(d): its selection (hλ)²/4 of Euler's, from the scalar identity |
| D3ac-8 | minor | the second mechanism's constant ∝ 1/a′; O(γ) only for γ ≲ 10⁻³ | (ii) and the numerics |
| D3ac-9 | minor | Part F cannot resolve ρ ≥ 0.99 | the floor printed, such rows marked; (vi) "checked for ρ ≤ 0.9" |
| D3ac-10 | minor | the prediction at W₀ only; no large-step test with enough seeds | the path-averaged prediction (`--reanalyse-C`); C2 at h = 8, 8 seeds (0.979–1.049) |
| D3ac-11 | minor | mutants A1, A2, H9 survive | tests of the amplitude coordinates' target and of Heun with dither; Gaussian sums; the finite-step slow limit |
| D3ac-12 | minor | literature (Lucchi et al., Beneventano, Kim et al., AOUP/UCNA, Li–Wang–Arora's log-det) | [S69], the audit; "new" narrowed to the coloured filter of the drift along M |
| D3ac-13–14 | minor/cosmetic | "a few percent" is 9.6%; "any step" is h < 2/λ_max; Part D's digits | fixed |

**Lesson.** The derivations held; the failures came from reading a global statement off a few points — a target named from one
continued branch, a limit's direction from three cosines, a sub-percent residual given a mechanism without varying the step. Map the whole
zero set and its basins before naming a target, vary the noise's spectrum and not only its correlation time, and extrapolate any small
numerical effect in the step of the computation that produced it.

## Round 31 (adversarial review of Propositions D‴.39 and D‴.40, the dithered gauge's dissipativity and the gauge under general noise, and of the round-30 revisions of D‴.41 and D‴.42, 2026-10-27)

**Design.** One referee, not modifying the repository. D3ad: Proposition D‴.39 (inserted at 1cefb5f), the draft of D‴.40 and the
round-30 revisions of D‴.41 and D‴.42, with their code at 7adc47e, Part E of `exp_gauge_noise.py` still running
(`review/REVIEW_round31_D3ad.md`): the generator identity term by term in the boxed stationary law (simulations at τ = 10⁻³, 10⁻² with
every term), the leading-order claim under correlated noise, Part E's 32 finished series (the correlated noise's effective levels), Part
G1b's zero search, the plateau's escape statistics, the box's regularity, the τ̂ estimator; 26 mutants.

**Verdicts (D3ad).** D‴.40(i)'s identity right only with the box's and the cohesion's terms (without the box no stationary law exists;
with it the box's term shifts E[u_k] by up to 5.4·10⁻² from the closed form at the fit); "only each amplitude's own intensity" exact for
the prefactor only — the O(τ²) term depends on all of Σ; Part E ran without the box; D‴.41's Part G1b missed stable zeros at γ = 10⁻⁵.
D‴.39's dead zone without bias only for the deterministic objects; its Arrhenius slope 0.81 ± 0.14 of the barrier, within 1.3 standard
errors of Freidlin–Wentzell's limit (the "82%" had no error bar).

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ad-1 | major (D‴.40 (i)) | the identity omits E[g_k∂_{α_k}V]; "exactly in any stationary law" vacuous without the box | (i) with the box's and the cohesion's terms, the boxed stationary law (exact for the diffusion, the Euler chain's O(Δt) bias measured), the box-free form in the well up to exit terms; Part D every term |
| D3ad-2 | major (D‴.40) | "only each amplitude's own intensity" only for the prefactor | (i): at leading order; the O(τ²) term E[g′_k] depends on all of Σ (Part D's amplitudes-only against isotropic) |
| D3ad-3 | major (Part E) | run without the box; the correlated noise heats the centres 4–12× | Part E inside the box with a correlated Σ at the amplitude levels (2, 1, ½)τ, the classes' levels variant: done — the centres and widths dithered only, the masses within 2·10⁻⁴ of their start; with the barrier at W♯ within 8·10⁻⁴; with the class levels toward W♯_τ (0.022 from it over the second half of 10⁶ slow units) |
| D3ad-4 | major (D‴.41 G1b) | three stable zeros missed at γ = 10⁻⁵; W₀'s flow into a corner | the zero search with reference levels and a capped step (`g_ref`, `g_zero`, `g_flow`): at γ = 10⁻⁵ the three stable zeros from τ_c = 100 on and the unstable centre from τ_c = 1000 (a8074e0); at τ_c = 10⁵ the flow from W₀ stops on its step budget short of its zero (exponents −2.2, −1.9) |
| D3ad-5–6 | minor | D‴.40(ii)'s zero drift fails in the boxed law; τ̂ biased by the drift | (ii) qualified, "no block coupling"; (iii): the noise increments Δα − Δt·b, or the O(Δt) bias stated |
| D3ad-7 | minor | Part D's χ² tests E[Lψ] only; the Euler bias 0.9–1.4%; the errors underestimated | Part D: the identity's residual and z, the chain-exact Euler estimator, martingale standard errors: done (the Part D rerun) — the box-free form fails by up to 58 errors, its residual the box's term; the whole identity χ² = 24.8 for the 18 values with amplitude noise (|z| ≤ 1.71; one seed for every run, so the diagonal-noise runs' twelve values repeat three streams), the boundary term (ψ_k(end) − ψ_k(start))/T recorded (below 0.03 errors); the centres-and-widths run (τ_k = 0, no martingale error) out of the χ², its identity the boundary term to 7·10⁻¹⁹ |
| D3ad-8–9 | minor | numbers ("3·10⁻²" is 5.4·10⁻², from the box); the stationary Stein identity's precedents (configurational temperature) | numerics, docstrings; [S68] (Rugh; Jepps–Ayton–Evans; Stein's method), the novelty narrowed to the gauge identity |
| D3ad-10–12 | minor (D‴.39) | "without bias" for the stationary averages; the slope's error bar; the C^{1,1} wall and MYULA | (iii)(b) the deterministic objects only; (0.81 ± 0.14)ΔJ; the layer flagged, [S67] |
| D3ad-13 | minor | mutants of the correlated-noise path (`noise_L`), the levels, the per-class form, the replicator flow survive; no test of the generator identity | `tests/test_gauge_noise.py` (+3 tests) |
| D3ad-14–16 | minor/cosmetic | D‴.41 (iv)(c)'s Ψ(∞) off by 2π; slips; typos | fixed |

**Lesson.** An identity that is exact for one generator is not exact for the law one actually samples: every force the gauge does not
cover — a box, a cohesion, a time step — enters the stationary average, and has to be measured, not assumed away. Name the stationary
law, list its forces, and test the identity term by term where the law reaches its walls.

## Round 32 (adversarial review of Proposition D‴.43, the colour of label noise selects the kernels, 2026-10-27)

**Design.** One referee, not modifying the repository. D3ae: the draft of Proposition D‴.43 and its code at 883425a
(`review/REVIEW_round32_D3ae.md`): the coordinate change and the reduced flow against 64-path means with the w-drift scaled; integration
horizons; success thresholds, magnitude laws, a 40 × 160 design, a fit floor; smooth and anti-correlated noise; R(truth) against R(end);
certificates against τ_c on kernel dictionaries and the flow's actual return; SBL, reweighted ℓ¹ (ε, start), IRLS-ℓp, OMP on the same
instances; Euler against Heun at the same physical noise; seeds 4–7 at B's setting and 0–7 at σ = 0.3, 0.5; 23 mutants.

**Verdicts (D3ae).** The reduced flow's formula right and followed by 64-path means at σ = 0.1; its recovery advantage over basis pursuit
robust at the level of the flow (threshold, magnitude law, 40 × 160) — but SBL matches or beats it; F's "no advantage from α = 2" a
horizon artefact; B's 4/4 seeds luck (4 more: 0, 0, 0, 2/16); the descent's agreement breaks by σ = 0.3; at h = 0.5 the white comparator
is not basis pursuit; the effect an O(h) term of the explicit step (Heun 25× smaller); (i)'s "proved" outside its spectral-gap
hypothesis in the pruning phase; 10 of 23 mutants surviving.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ae-1 | major | F's α = 2 claim a horizon artefact | F run to convergence (`flow_to_convergence`): from α = 2 the slow flow recovers 0.35, 0.30 (k = 3, 4), at basis pursuit's rate (0.45, 0.25) |
| D3ae-2–3 | major | the noisy outcome depends on σ and on the seed | B over 8 seeds × 16 paths, the σ sweep per seed: the descent recovers where its reduced flow does (4 of 8 seeds, 16/16 paths), at σ = 0.3 on three of them, at σ = 0.5 on one; (i): the flow as the small-σ limit |
| D3ae-4 | major | SBL missing; reweighted ℓ¹ unfair; Wipf–Nagarajan's theorem | G with SBL and reweighted ℓ¹ at three ε from two starts: SBL 0.88, 0.55 against the slow flow's 0.75, 0.60 (k = 3, 4), reweighted ℓ¹ 0.38–0.53, 0.17–0.25; the novelty narrowed to the mechanism; [S70] |
| D3ae-5 | major | white = basis pursuit only as h, σ → 0 | (ii) restated; the finite-step comparator labelled; Damian–Ma–Lee's Ψ |
| D3ae-6 | major | an O(h) effect of the explicit step | the setting states it (rate ∝ hσ², none in continuous time); Part J, Euler at two steps against Heun: Δlog P off the support −2.39 (Euler, h = 0.5), −0.278 (Euler, h = 0.125), −0.097 (Heun, h = 0.5) |
| D3ae-7 | major | (i) proved only with a spectral gap | Status: proved where rank K = N and λ_min(K) ≫ v_h; the pruning phase and (v) the reduced flow's |
| D3ae-8–10 | minor | the certificate's τ_c threshold; Fuchs's, not the exact recovery condition; "spark" | (v) restated; H on kernel dictionaries with the flow's return and basis pursuit's success: certified in 0.85–1.00 of the instances at τ_c = 10⁵ (white 0–0.60), the flow returning in 0.90–1.00 |
| D3ae-11–12 | minor | concavity only in the continuum (slope ∝ 1 − 3ρ); the SBL analogy | (iv) in the continuum, the finite-step expansion; Fazel's log-det heuristic, SBL's log\|Σ_y\| term |
| D3ae-13–14 | minor | the truth not always the global minimiser; D's instance 13 misread | C: R(end) against R(truth); D: the flow to convergence beside the noisy runs |
| D3ae-15 | minor | mutants survive (the w-drift's factor, uv_of, the projection's truncation, the experiment's helpers) | three tests: the drifts against R's gradient in its metric (a doubled w-drift fails), `uv_of`'s round trip with the singular projection, the experiment's `ok` and `certificate` by hand with SBL, the convergence criterion and Heun's step |
| D3ae-16–19 | minor | a fit floor; a smooth noise's slow-down; unconverged white flows; "recovers" for the interpolant | Part I; E at ten times the horizon; C's white flow to convergence; "identifies the support" |
| D3ae-20–21 | cosmetic | slips; [S69]/[S70] attributions | fixed |

**Lesson.** A striking effect has to beat the classical method that minimises the same kind of penalty on the same instances, and every
"all" or "never" needs a convergence criterion, more seeds and a sweep of the noise's size before it is written; a selection read off a
discrete scheme's Itô drift is a property of the scheme, and the README and the audit must say so as plainly as the statement.

## Round 33 (adversarial review of Proposition D‴.45, the free mass–amplitude law as a mirror flow, 2026-10-27)

**Design.** One referee, not modifying the repository. D3af: the draft of Proposition D‴.45 and its code at 06247cb, with the results
then uncommitted (committed as 23b4b80 during the review, the code unchanged) (`review/REVIEW_round33_D3af.md`): the identities by hand;
η against quadrature and a 60-digit reference to q = 10¹⁰, both branches, the potential's shape and its rich remainder; the exact Euler
means by Gauss–Hermite at h = 0.5, 0.05, 0.005; the seven instances with unequal-mass rich distances above 10⁻³ re-solved by an
independent Newton method on the dual; the loop from three branch starts to t = 10⁷ and its time to the projection at s = 10⁻² … 10⁻⁴;
D‴.44's reduced flow to s = 270, and to s = 10⁸ where basis pursuit is not the truth; C(3)'s end values of Q against (iv)(b)'s zero; 36
mutants of `mirror_free.py` against the tests and a reduced copy of Parts A–E.

**Verdicts (D3af).** The identities right: Q conserved for every drive, η's closed forms on both branches (1.5·10⁻¹⁴ against quadrature),
dw/dη = mD and ψ*, the dual flow's convergence and the lazy limit, the Euler step's ΔQ — exactly D‴.44's drift (C(3)'s end values of Q
to 7 digits) —, V̇ = −k_ee². Wrong or overstated: the rich remainder (−(3/2)K|w|^{2/3}, not O(|w|(1 + |log|w||))); the rich limit needing
uniqueness, reached at scales set by the vertex gap and in a time ∝ 1/β⁰² (t ≈ 10⁶, 10⁷, 10⁹ at s = 10⁻², 10⁻³, 10⁻⁴); Part B's medians
hiding silent solver failures (|w| ≈ 2.2·10²⁹²) and a bimodal branch median; C(3)'s "basis pursuit" only where basis pursuit is the
truth; "never moves" true for noise with smooth paths (Stratonovich) only; C(2)'s 14–25 standard errors the O(h³) term; (ii)–(iii) an
instance of Li–Wang–Lee–Arora's theorem; 9 of 36 mutants passing the tests.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3af-1 | major ((iv)(b), C(3)) | "reaches the rich limit (basis pursuit)" holds only because all six C(3) instances have basis pursuit = the truth; elsewhere D‴.44's reduced flow stops on other sparse fits (s to 10⁸), and nothing backs basis pursuit under the paradigm's update (D‴.44(ii)'s D* is the flat scheme's) | (iv)(b): the state first moves away from basis pursuit, then settles on a sparse exact fit — basis pursuit's where it is the truth, elsewhere its support or another sparse fit of larger ℓ¹ norm (numerical); C(3) adds instances 8, 9, 13, 17: the six end on basis pursuit (2.7–11.6·10⁻⁹ at s = 10⁸), of the four two on its support (2.6·10⁻⁵, 6.5·10⁻⁴) and two on other fits (‖w‖₁ 0.4228 against 0.4156, 0.4714 against 0.4516) |
| D3af-2 | major (B, code) | Part B's medians include silent solver failures: `bregman_projection` returned no residual or flag, inherited `alpha_of_eta`'s clip at \|α\| = 340, and diverged to \|w\| ≈ 2.2·10²⁹² (rich, s = 10⁻⁷) and on two branch starts | the projection by Newton on the convex dual with Armijo backtracking, an adaptive radius, rejected trial points and the best iterate's true residual (`dual_newton`); `alpha_of_eta` raises instead of clipping; every aggregate counts failures and leaves them out — none in the rerun (8feae75, 7421bcf; 7421bcf adds a decimal-arithmetic Newton step where the float step stalls at extreme masses) |
| D3af-3 | major ((iii), B) | the branches misreported: the k = 2 median 6.6·10⁻⁴ the middle of a bimodal set (ten instances ≤ 1.2·10⁻⁴, eight at 0.05–0.37), "signs kept" true by construction, "an entropy against the others' ℓ¹ weight" wrong (costs 340–480 against 420–1220, λ₁ = 702), exact fits lost as s → 0 | (iii)'s branch sentences replaced by Proposition D‴.47 (the Bregman projection onto the sign-constrained least-squares fits; an exact fit may cease to exist); B reports the distribution (k = 2: 10 of 20 within 10⁻³, the largest 0.37), the instances with a truth atom on a wrong-signed branch (13 and 18 of 20 at s = 10⁻⁶) and those whose frozen signs still fit y (18, 19); "signs kept" dropped |
| D3af-4 | major ((iv)(a), pieces) | ROW0, BULLET and README: continuous-time noise of any colour "never moves that projection" — false for white noise read in Itô's sense (Pesme et al.'s stochastic gradient flow), whose effect is (iv)(b)'s | "noise with smooth paths (Stratonovich)" in the setting, (iv)(a), ROW0, BULLET and README; Itô's white noise adds (b)'s drift; C(1) discriminates as h → 0: Euler–Maruyama's mean change of Q over its Itô correction's integral 1.0106, 1.0037, 1.0011 at h = 0.1, 0.03, 0.01, Heun's 1.4·10⁻⁵ → 1.4·10⁻⁶ (∝ h), the mirror update's 2.2·10⁻¹⁶ (by construction) |
| D3af-5 | major (Status, [S72]) | novelty: (ii)–(iii) Li–Wang–Lee–Arora's theorem (with Marcotte–Gribonval–Peyré, Azulay et al.), (iv)(a) stochastic mirror descent's invariance, (iv)(b) Pesme et al.'s and Even et al.'s mechanism, (v) Boffi–Slotine's, the shape the hypentropy | Status, ROW0 and the audit recast these as classical, with Gunasekar et al., Azizan–Hassibi, HaoChen et al., Ghai–Hazan–Singer, Woodworth et al. ([S72]); new narrowed to the unit's explicit potential (Q, the wall, the entropy tail, the 1/√β⁰ weights and their correction, the branches), the identification with D‴.44's selection, the mirror update |
| D3af-6–7 | major/minor ((v), (i)) | the pair claim and Part D lie inside the cohesion's range (it acts on α_i − α_j and breaks Q); "also when centres and widths move" needs a qualifier | (v): "two particles on one atom without cohesion between them"; (i) and the setting: "as long as the cohesion's α-component vanishes on the particle" |
| D3af-8 | minor ((iv)(b)) | the selection's early phase moves away from basis pursuit (‖w‖₁ 0.706 → 0.779 over s ≤ 270); C(3) sampled only s ≥ 10⁵ | (iv)(b): "first moves away from basis pursuit (‖w‖₁ rising)"; C(3) samples s = 10², 10³, 10⁴ (instance 0: 0.706 → 0.781 by s = 10³) |
| D3af-9 | minor ((ii)) | the rich remainder is −(3/2)K\|w\|^{2/3}, not O(\|w\|(1 + \|log\|w\|\|)) | (ii): ψ = λ₁\|w\| − (3/2)K\|w\|^{2/3}(1 + o(1)), λ₁ exact; the proof's remainder O(\|w\|^{2/3}/λ₁) |
| D3af-10–11 | minor ((iii)) | the rich limit needs uniqueness and is reached at scales set by the vertex gap (k = 3, instance 13 still 0.17 away at s = 10⁻⁷); the time ∝ 1/β⁰² unstated | (iii): "when it is unique"; the approach a power at an N-atom vertex — the off-support coefficients on the core's tangent profile, a small slack delaying it —, the transit's time ∝ β⁰^{−1/2} and the last digits' ∝ 1/β⁰² (fewer than N atoms); B reports medians and largest values and the loop at s = 10⁻³ (within 10⁻⁸ of its projection from t ≈ 8.5·10⁶); Part F (7421bcf): the leading order to the rounding floor, the distance falling by 1068 and 984 from ε = 10⁻⁷ to 10⁻⁹ (ε^{3/2}), instance 13 (μ_min = 2.7·10⁻⁴) 0.17 away at 10⁻⁷ and 3.2·10⁻⁸ at 10⁻¹¹ |
| D3af-12 | minor (code) | η NaN at q = −c and not cancellation-free for \|t\| ≥ 0.9; `alpha_of_eta` clips silently and its branch bracket can fail; `mirror_step` unchecked | η without cancellation for every t and level (4·10⁻¹⁶ against a 60-digit reference to q = 10¹⁴), the critical curve, the branches from log sinh of the distance to the branch point; `alpha_of_eta` and `alpha_of_w` raise and solve on a branch in log δ; `mirror_step` checks its convergence (8feae75) |
| D3af-13 | minor (numerics) | checks true by construction ("signs kept", Q under the mirror update, Part D's split, A's shift, C(1)) | A's shift integrated along the trajectory (1.5–4.4·10⁻⁷); C(1) Euler–Maruyama, Heun and the mirror update as h → 0; the mirror update's Q marked "by construction"; D: the pair's common mirror increments |
| D3af-14–15 | minor ((iv)(a), Status) | the stationary mean under persistent noise displaced at O(v); the Status's "(iv) proved", "(vi) proved", (iv)(c) | (iv)(a): one clause; Status: (iv)(b) the formula proved, the selection's end numerical, (iv)(c) in expectation, (vi)'s identity proved ("not a mirror flow" dropped) |
| D3af-16 | minor (D) | Part D's "sensor noise" reaches only the adaptation; "parameters converge slowly" not shown | D: a smooth noise in the adaptation only and a sensor noise in the feedback, the regressor and the adaptation (Q kept to 6.8·10⁻¹² and 4.9·10⁻¹³); after round 36 the noiseless parameters' distance at t = 2000 (0.048, 0.027) in place of "converge slowly" |
| D3af-17–19 | minor | [S72] credits "isotropic noise does not" to Pesme et al. (it is HaoChen–Wei–Lee–Ma's); no link to D‴.41(ii)'s C ∝ I; the lost bound \|f̂\| ≤ a_max | [S72] fixed; (iv)(a): C ∝ I is the conserved law's centring; (v): the bound lost, D‴.46's leak the remedy |
| D3af-20 | minor (tests) | 9 of 36 mutants survive (M04, M09–M13, M22, M24, M36); the V̇ test blind to the mass law; no rich, branch or range test | `tests/test_mirror_free.py` 8 → 15 tests: η against a 60-digit reference, round trips and range errors, rich projections at 10⁻¹⁰, `real_lyapunov` reading the state's mass, branch mirror steps (round 36: M10n, M11n, M24n killed, M12n only slower) |
| D3af-21–22 | cosmetic | "falling 0.14, 0.16" (it rises); "as often as" (0.75 against 0.80); C(2) by Monte Carlo; the proof's dropped (c − Q) | fixed (the first two after round 36); C(2) by Gauss–Hermite (80 nodes) |

**Lesson.** The exact identities were right and cheap to check; what failed passed around them unexamined — medians over a solver whose
convergence was never recorded, a claim generalised from instances on which it holds automatically, checks true by construction, and
shortened statements that dropped the one qualifier (Stratonovich) the theorem needs. Every aggregate carries its failure count, a
numerical claim is tested where it can fail, and every summary keeps its hypotheses.

## Round 35 (adversarial review of Proposition D‴.46, a leak in the mirror coordinates, 2026-10-27)

**Design.** One referee, not modifying the repository. D3ah: the draft of Proposition D‴.46 and its code and results at 67c08b1
(`review/REVIEW_round35_D3ah.md`): the three-point identity by quadrature, ψ rebuilt independently of the module; the averaged loop with
the error filter L = 1/(s + k_e) (the periodic steady state by FFT on the 20π-periodic orbit) integrated to 2·10⁵; the real σ = 0 loop
continued to t = 18000; Part C's equilibrium three ways (the idealised H/k_e, the filtered H_L, the exact frozen-parameter error); the
irrepresentable condition over all 1760 three-sparse patterns of Part D's Gram; the rich start's averaged dynamics over 10⁸ time units;
an instrumented `leak_equilibrium`; 25 mutants of the leak's code and 8 of the experiment's report code.

**Verdicts (D3ah).** The algebra right: the Lyapunov identity (1.3·10⁻¹⁵ by quadrature, 1.2·10⁻¹⁰ in the module) and its ultimate bound,
B(w, 0) = ψ(w), the ridge limit, the LASSO limit at a fixed penalty σλ₁. What the numerics were said to show did not hold: Part B's
"drift" is slow convergence to the biased fit (τ ≈ 1.5·10⁴ against T = 3000); (ii) is the k_e → ∞ limit — the averaged loop has the
error filter's Gram H_L and d̄_L, whose lag moves Part C's lazy equilibrium by 8.2·10⁻⁴, more than the reported gap, and the intermediate
run lasted four time constants; "prunes to the target's kernels" one instance (the irrepresentable condition holds for 23.6% of the
patterns); w_σ → w* + H⁺d̄, not w*; the dual-coordinate leak established in learning in games; 8 of 25 mutants and all 8 report mutants
surviving.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ah-1 | major (B) | Part B's "drift" is slow convergence: on a persistently exciting orbit the averaged σ = 0 loop converges to the biased fit w* + H_L⁻¹d̄_L (max\|w\| 0.337, τ ≈ 1.5·10⁴, against T = 3000), and the real loop continued to t = 18000 approaches it | B(1) reworded ("creeps toward the biased fit"), the σ = 0 run to 2·10⁴ (1.3τ; the quarters' largest \|w\| 0.309, 0.330, 0.337, 0.341, the orbit's oscillation included); B(2), a drift without persistent excitation (φ = (1 + t)^{−1/2}, a bounded d built for it): w grows like (1 + t)^{1/4} (23.3 at t = 10⁶), σ = 0.01 and 0.1 hold it at 0.46 |
| D3ah-2 | major ((ii)) | (ii) is the k_e → ∞ limit, not the slow-adaptation limit: the averaged loop has H_L = k_e⟨φL[φ]ᵀ⟩ and d̄_L (0.9% from H at k_e = 20, symmetric on a scalar periodic orbit) | (ii) restated with H_L and d̄_L: the unique, globally attracting root of a monotone equation (Lyapunov function B(w_σ, w)), the mirror flow of F_σ to O(k_e⁻²), the real loop's mean by the averaging theorem (cited); `filtered_gram` |
| D3ah-3 | major (C) | Part C too short (T = 4/σ for the intermediate start) and its drift the lazy start's within-period oscillation | C started at w_σ^L, with period means (20π): the last within 1.6·10⁻⁴ (318 periods) and 1.7·10⁻⁴ (954) of w_σ^L, 4.4 and 2.6 times closer than to the idealised w_σ (lags 8.2·10⁻⁴, 5.1·10⁻⁴; the leak's bias 0.038, 0.12) |
| D3ah-4 | major ((iii), D, pieces) | "the leak prunes to the target's kernels" one instance: on Part D's Gram the irrepresentable condition holds for 416 of 1760 three-sparse patterns; where it fails the equilibrium keeps spurious atoms at every μ | (iii): exact support recovery under the irrepresentable condition and no sign flips (weighted after round 36); D(2) over all 1760 patterns at μ = 10⁻³: the condition 416, the prediction 416, the equilibrium exact on 415, agreement on 1745; README, BULLET, ROW0 and the audit: "under the irrepresentable condition" |
| D3ah-5 | major ((ii)) | w_σ → w* + H⁺d̄, not w*, under a biased disturbance; the d̄ path never exercised (its sign mutant survives) | (ii): w_σ → w* + H_L⁻¹d̄_L; D(3): the distance to w* + H⁻¹d̄ linear in σ (3.7·10⁻⁸ at σ = 10⁻⁸), 0.043 from w*; a test with k_e = 3 and a bias |
| D3ah-6 | major (novelty) | a dual-coordinate leak, its regulariser-perturbed rest points, Bregman (Fenchel-coupling) Lyapunov functions and the vanishing-leak selection are established (Gao–Pavel; Romano–Pavel; Coucheney–Gaujal–Mertikopoulos; Pérolat et al.; Abe et al.; Azizan–Lale–Hassibi; Alvarez–Bolte–Brahic) | [S73] and the audit with these; new narrowed to the (α, β) implementation keeping Q, the bound in the adaptive loop, the weights λ₁ ∝ β⁰^{−1/2}, the filtered averaged equilibrium as a monotone equation |
| D3ah-7–9 | minor ((iii)) | the remainder cited from D‴.45 (O(\|w\| log\|w\|)); λ₁ only to leading order (10.6% off at β⁰ = 0.4c); "soft thresholding" | (iii) with D‴.45's −(3/2)K\|w\|^{2/3}, the exact λ₁ = πa_maxκ′/(4Gc), the LASSO's KKT condition |
| D3ah-10 | minor | the rich (LASSO) regime never run dynamically; its time ∝ k_eβ⁰^{−1/2} | D(4): the averaged loop from β⁰ = 10⁻⁶ at its equilibrium at t ≈ 3.2·10⁵ (2·10⁶ where the condition fails, five atoms kept); (iii) states the escape time |
| D3ah-11 | minor (B) | the bound check without power: loose, not asserted, five single mutants of its formula surviving | the comparison bound asserted along B's runs (the printed largest ratio, 1.0000, still the t = 0 value); a test under a worst-case disturbance in which V comes within a factor 2 of C/λ |
| D3ah-12–13 | minor (tests, code) | 8 of 25 mutants survive (the noise path, the k_e scalings, d̄'s sign, the tolerance, the damping); `leak_equilibrium` fails silently (it_max at β⁰ = 10⁻⁸, rejected trial points returned, clamping) | tests with a measurement noise, k_e ≠ 1, a bias and the solver's residual (round 36: R6, R7, K1, K2, P1, T1 killed; `lasso` tested at k_e = 5 after it, 2773f33; D2: fix not traced); `leak_equilibrium` returns (w, α, iterations, residual), caps its trial steps in η, rejects points beyond the curves, keeps its best iterate |
| D3ah-14–15 | minor ((i), code) | (i)'s hypotheses incomplete (w* in each curve's range, V proper, the split of a non-representable f); η from the nominal q, not Q from the state | (i): the range, properness, the split with the bound minimised over it, "along the loop's own trajectory"; `real_rhs_leak` computes Q from the state (continuous time only; no Euler run) |
| D3ah-16–17 | minor/cosmetic | Part A's ridge rows test only linear response (the ridge moving w by 10⁻⁹–4·10⁻⁵); the ridge "support" a threshold artefact; absolute and relative distances unlabelled; "V rising from 0.18" quotes V(0) | A: ρ = 0.5 (w moved by 31–55%), and β⁰ = 3 with the wall inside the coefficients' range (0.15 off the ridge); D's lazy rows "dense (no coefficient vanishes)", the distances labelled; B: V 0.40 → 0.45 |

**Lesson.** Every agreement the proposition claimed sat at the size of an effect its own setting neglected — the error filter's lag, the
slowest averaged time constant, the irrepresentable condition. Before calling a number confirmation, compute the neglected quantity for
the same instance: here each took seconds, and each changed what the paragraph may claim.

## Round 36 (adversarial review of Proposition D‴.47, the free law's branches, and of the round-33 and round-35 revisions of Propositions D‴.45 and D‴.46, 2026-10-27)

**Design.** One referee, not modifying the repository. D3ai: the draft of Proposition D‴.47 with the round-33 and round-35 revisions of
D‴.45 and D‴.46, and their code at 7421bcf (`review/REVIEW_round36_D3ai.md`): `exp_mirror_branch.py` rerun (byte-identical) and every
number of the three numerics paragraphs checked; the identities at a_max = 0.3 and 2.5; mixed starts whose branch signs cannot fit y;
rank-deficient and duplicated designs (`sign_ls_set` against nnls); inexact fits whose least-squares set is not a point; degenerate
vertices; the nonnegative certificate's margin over every three-atom support, with 80-digit certificates; the loop to t = 10⁹ on typical
targets; where the float Newton fails; 68 mutants (54 of `mirror_free.py` against the 18 tests and a reduced experiment).

**Verdicts (D3ai).** D‴.47's proofs sound — (ii)'s zero-sum argument airtight (the limit confirmed on a rank-4 design, duplicated atoms
with an inexact fit and mixed starts, to ≤ 5.4·10⁻⁹), (iii)'s rate right at N-atom vertices, (i), (iv)'s conjugate, (v)'s proof and (vi)
correct —; every failure in what the numerics were chosen to show: every inexact instance had a single-point least-squares set (the
Bregman selection with extinction untested, the branch speed's mutants surviving); the mixed start checked only where the branch atoms'
signs fit y; Part E's targets from the top decile of certificate margins (over all three-atom supports the medians 4.0·10⁻⁶ and
3.8·10⁻⁹, 25% and 45% numerically zero); `sign_ls_set` silently wrong on rank-deficient and duplicated designs; the decimal fallback
unbudgeted. D‴.45's and D‴.46's round-33 and round-35 items fixed but for leftovers (D3ai-17, 18, 20, 21). 28 of 68 mutants surviving
the tests.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ai-1 | major ((v), E) | "converges to w* from every positive-branch start, at every mass scale" rests on top-decile targets: over all three-atom supports median margins 4.0·10⁻⁶ and 3.8·10⁻⁹ (25%, 45% numerically zero; exact margins 10⁻³²–10⁻³⁴); typical targets 0.45–0.63 away at t = 10⁹, with 0.6–1.2 of mass off the target | (v) restated: the selection's speed is the margin (Σ_off w ≤ (\|v\|/μ)(2NB/t)^{1/2}), uniqueness alone does not make the target reachable, stability constants (Eftekhari et al.) govern it; E: the margins' quantiles on the three grids (25%, 45%, 4% at most 10⁻⁹) and the loop on supports at the quantiles — reached by t = 10⁹ above a margin of 10⁻³, still 0.37–0.83 away below 10⁻⁴ (ebfc9d9); BULLET, ROW0 |
| D3ai-2 | major ((iii), C) | "the free particles keeping O(ε² log(1/ε))" false when the branch atoms' frozen signs cannot fit y (1.36 and 1.10 at ε = 10⁻⁴, 10⁻⁸); Part C's check true by construction | (iii) with the hypothesis, and otherwise the limit basis pursuit over the free atoms with the branch atoms sign-constrained and free of charge; C: an infeasible-sign row — the free particles carry 2.27 at ε = 10⁻⁴ and 10⁻⁸, 4.1·10⁻⁴ and 2.5·10⁻¹¹ from that limit; at 10⁻¹² the projection fails and is reported as failed |
| D3ai-3 | major (code) | `sign_ls_set` never read lsq_linear's status: a wrong ŷ on rank-deficient (residual 0.103 for 0.520) and duplicated designs (\|Z\| 0 for 14) | Lawson–Hanson's active set (exact at any rank), its KKT conditions checked, the linear programs' status read (17d8396); D‴.45 Part B's ŷ (to 4.6·10⁻¹⁴) and Z unchanged on all 120 instances (2773f33) |
| D3ai-4 | major (tests, B) | the Bregman selection with extinction never tested (S a point in every inexact instance); the branch-speed mutants S1, S3, S4 pass all 18 tests | a test (duplicated atoms, Z ≠ ∅, negative branches alive: the flow at the projection to 10⁻⁹, not at the equal split); B: 5 × 4 designs with every atom carried by two particles — the flow at the Bregman projection to 6.4·10⁻¹⁰ on 8 seeds, the projection 1.6·10⁻⁴–0.20 from the equal split |
| D3ai-5, 13 | minor (C, D, code) | the experiment neither asserts nor counts projection residuals (a 40-digit mutant tabulates failures silently); the decimal fallback can run for hours (10 s per 3000-digit step) | C and D count and flag failed projections (C's one, at ε = 10⁻¹²); `dual_newton` budgets its exact steps (count and wall time) and returns the residual (2773f33) |
| D3ai-6, 7 | minor ((iii), (ii)(c)) | degenerate vertices (\|S\| < N: ε^{0.13}, ε^{0.26} observed) silently excluded; "no other coefficient vanishes" false for free particles; the critical curve's extinction algebraic | (iii): at a vertex with fewer than N atoms the approach ε^{μ*} stated (μ* not computed); (ii)(c): "no other branch coefficient", the decay on the critical curve about t⁻³ |
| D3ai-8, 9 | minor | (iv) drops "without cohesion between them" (a D3af-6 regression); BULLET's unconditional "exponentially"; ROW0's rate wording | restored; BULLET and ROW0: exponentially when strictly complementary; ROW0's rate ε^{μ} |
| D3ai-10, 11 | minor (literature) | lit_S74 and lit_S72 put D‴.45's α = 0 rate under the exponential penalty at ε^{1/2} (it is the power ε^{3/2}); novelty overstated; Chou–Maly–Stöger missing | corrected (D‴.47(iii): ε^{3/2} from the core's tangent, not an entropy); Status: (ii) Alvarez–Bolte–Brahic's convergence with the range argument, (iv) the u² − v² parametrisation's, (v)'s identifiability Eftekhari et al.'s; Chou–Maly–Stöger cited |
| D3ai-12 | minor | a free sign is not enough for D‴.45(iii): η(0) = (2/π)arcsin θ·λ₁, the weights λ₁(1 ∓ (2/π)arcsin θ) | stated in D‴.47's setting and in D‴.45(iii): the symmetric limit needs θ → 0 |
| D3ai-14, 15 | minor/cosmetic | the rows where the leading order fails omitted (it holds only from ε ≈ 10⁻¹⁶, u ≈ 10⁻¹²⁸); E's time labels; "≈ 1/t" (t^{−0.88}) | C and D report every row; E prints the sampled times; t^{−0.9}–t^{−1} |
| D3ai-16, 19 | minor (tests) | untested helpers and solver rules (FG1, IR1, NC1, R2, T3, T4, X4, X5, L4, E4, Z2 survive the tests; M12n only slower) | tests of the filtered Gram's causality, the weighted irrepresentable value, the nonnegative certificate's margin and the LASSO at k_e = 5 (2773f33); the solver's rules and an iteration or time bound for M12n: (open) |
| D3ai-17, 18 | minor (D‴.45) | numerics leftovers ("falling", "as often as", the bimodal median, "converge slowly", "atoms die", Q kept by construction reported as a result); the closed form's large constants (113–138% at ε = 10⁻⁷) and the μ⁻³ hidden in O(β⁰^{3/2}); the 1/β⁰² time unqualified in ROW0 and BULLET | D‴.45's numerics edited (B's spread, "nearly as often", "no atom yet extinct", "(by construction)", F's largest values and ranges); (iii): O(β⁰^{3/2}μ_min⁻³), the relative error β⁰^{1/2}/μ with a large constant; ROW0, BULLET: "where the support has fewer than N atoms" |
| D3ai-20, 21 | minor/cosmetic (D‴.46) | the irrepresentable condition unweighted for unequal masses, "iff < 1"; leftovers (the audit's "sparsity-selecting", the ridge "support", the σ = 0 run's 0.341 above 0.337, "V/bound 1.0000" its t = 0 value, D2 and L1 surviving) | the weighted condition, < 1 sufficient and ≤ 1 necessary (`irrepresentable` takes the weights); the audit's and D's wording, 0.341 attributed to the orbit's oscillation, `lasso` tested at k_e = 5; the printed bound ratio and D2: (open) |

**Lesson.** A confirmation is only as good as the instance it was obtained on: whether the solution set is a point, how the margin is
distributed over targets, whether the hypothesis holds by construction, whether the reference solver reported success. Ask what a wrong
answer would look like on that very instance before calling a number a check.

## Round 37 (adversarial review of Proposition D‴.48, the rich transit, 2026-10-27)

**Design.** One referee, not modifying the repository. D3aj: the draft of Proposition D‴.48 and its code at 5995696
(`review/REVIEW_round37_D3aj.md`): `exp_mirror_transit.py` rerun (the md byte-identical) and every number checked against the md and the
raw file; `iss_path`'s defining properties on three designs; the loop against the path to ε = 10⁻¹² with a Q-preserving integrator; the
creep's instantaneous speed against both formulas; ties, near ties, close and off-grid bumps, coarse grids, edges, amplitudes; the
Wasserstein–Fisher–Rao limit along trajectories with unequal masses and the formed-kernel count to ε = 10⁻¹⁰; one unit's transit time
for α(0) = θε^p; `sliding_iss`'s sensitivity to its seed and fast horizon; 53 mutants.

**Verdicts (D3aj).** (o)'s expansions and (i)'s soft-analysis limit (Mosco, then Brezis–Pazy) right, with slips (the energy identity's
factor N, the composition step, the inactive margins); the ε^{1/2} rate cleaner than reported (design 5's 0.471 an artefact of the
(α, log β) integration drifting off Q; 0.4998 on the curves); (ii) right but practically invisible (at ε = 10⁻⁴ the dual 40–80% of the
weights away); (o)'s dichotomy a continuum (the scale set by α(0)/ε^{1/2}); (iii)(a)'s coefficient bound right, its creep speed too large
by π³/8 and missing b^{3/2} and the sign; (iii)(b) failing at ties (a bump midway between two kernels formed by both, 0.40 each, at every
ε); (iv)'s condensation real, but "several kernels per feature" not a limit property (15, 12, 8, 4, 3, 3, 3 kernels at ε = 10⁻⁴ … 10⁻¹⁰;
exactly three on 16 of 16 random targets at 10⁻¹⁰), so no dichotomy of counts; `sliding_iss` neither faithful nor robust; 27 of 53
mutants surviving (about 20 not equivalent).

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3aj-1 | major ((iv), title, New) | "several kernels share each feature" and the metric dichotomy false as ε → 0 | (iv), the title, New, ROW0, BULLET and the audit restated: the massless kernels' positions condense, the count reaches one per feature only at much smaller masses; C2 and C4 to ε = 10⁻¹⁰: 15, 12, 8, 3, 3 kernels at ε = 10⁻⁴, 10⁻⁵, 10⁻⁶, 10⁻⁸, 10⁻¹⁰, always in three clusters at the bumps; on the random targets exactly three at the bumps on 0, 10, 16 of 16 at 10⁻⁶, 10⁻⁸, 10⁻¹⁰ (c1b3334) |
| D3aj-2 | major ((iii)(b)) | ties: a bump midway between two kernels is formed by both at every ε; near ties non-uniform; `sliding_iss`'s argmin tie-break not the limit | (iii)(b) with ties (shared, half the amplitude each) and genericity; `sliding_iss` enters ties together; C5: the tie shared (0.4, 0.4) at every ε, a near tie formed by two kernels at ε = 10⁻⁶ and by one from 10⁻⁸; bumps 0.15 and 0.25 apart off the grid resolved by two kernels, the second entering earlier than on the path (τ = 16.2, 29.6, 33.1 against 33.5; 12.1, 16.1, 17.0 against 17.1) — the creep; the near-tie law later Proposition D‴.51 |
| D3aj-3 | major ((iii)(a), C6) | the creep speed (πa_max²/(4c))εμ^{−3}⟨r∂_cφ⟩ too large by π³/8, without b^{3/2} and the sign; C6 compared a displacement since τ = 0 with speed × window | (iii)(a): the core's profile (πa_max²/(4c))εb^{3/2}cos(πμ/2)/sin³(πμ/2)·s⟨r∂_cφ⟩; C6 the instantaneous dc/dτ — the profile's to 0.990–0.996 (ε = 10⁻⁸) and 0.9999–1.0000 (10⁻¹⁰), 0.55 at μ = 4.5·10⁻⁴ where the kernel is entering, 0.255–0.258 of the μ^{−3} form (8/π³) —; a test (to 2·10⁻³) |
| D3aj-4 | major (code) | `sliding_iss` neither faithful nor robust: exits without an inward check, a fixed fast horizon, the path dependent on the seed (t_fast = 10⁴: 0 of 17 targets) | ties entered together, the active block relaxed until stationary (the elapsed time doubled), exits only with an inward correlation (17d8396); a test: seeds 10⁻³ and 10⁻¹¹ and t_fast = 10³ give the same path; C4: 16 of 16 random targets exact, no exits |
| D3aj-5, 14 | major/minor (tests) | 27 of 53 mutants survive (`iss_path`'s exit rule, NNLS and hitting time; `moving_rhs`'s m and a′; the WFR limit's k, tan and 1/ρ; `sliding_iss`'s parameters; ParticleSet's free drive); the test design has no exits; the experiment untested | tests: an independent forward-stepping ISS integrator on designs with exits, a binding sign constraint, the drifts at a random state, the creep profile, the WFR limit along the trajectory with unequal masses, ties and seeds, the free particles' drive (17d8396); tests of the experiment's own code: (open) |
| D3aj-6, 7 | minor ((o), (iii)) | (o)'s dichotomy a continuum (crossover at α(0) ~ ε^{1/2}); the creep generic on close bumps, "well-separated" unquantified | (o): the scale set by θ = \|α(0)\|/(β⁰/(2c))^{1/2}, α(0) ∝ ε^p giving ε^{−p}log(1/ε), a branch start with α(0) ≲ ε^{1/2} on the power scale; (iii)(a): frozen only while ε ≪ μ³; C5 off-grid close bumps |
| D3aj-8 | minor (A, C) | the (α, log β) integration drifts off Q at small masses (design 5's slope 0.471) | each atom on its curve, β = εb + 2c sinh²(α/a_max) (`critical_rhs(…, b0)`, `FreeLaw.rhs_curve`); A to ε = 10⁻¹² and C rerun: slopes 0.500 on the three designs between 10⁻⁶ and 10⁻⁸, 0.503, 0.501, 0.484 between 10⁻⁸ and 10⁻¹² (c1b3334) |
| D3aj-9, 10 | minor (C2, C4) | the path's kernels entered on 15 of 16, not 14; the WFR limit checked only at its attractor, with ρ = 1 | net entries compared (15 of 16 at 10⁻⁸); C2 and a test: the limit along the trajectory with unequal masses (1.1·10⁻⁷ with ρ_i = b_i^{−1/2}, 0.45 with the masses ignored) |
| D3aj-11 | minor (proof of (i)) | the energy identity's factor (N∫\|ζ̇\|²), the composition's Rockafellar–Wets step, the superfluous −O(1), the inactive margins unstated | fixed |
| D3aj-12, 13 | minor | novelty: (ii) Berthier–Pillaud-Vivien's entropic case, (i) their route, the dichotomy refuted; Kunin et al., Labarrière et al., Maennel–Bousquet–Gelly missing; (ii) practically invisible (the ± pairs' first ratio 2.65) | New and the audit narrowed, the three cited; (ii) qualified as asymptotic (at ε = 10⁻⁴ the dual 40–80% of the weights away) |
| D3aj-15 | cosmetic | "(≈ ε)", "17–50" (18.8–40.5), the bound called a speed, c_j(0) where the code slides c, the margin's normalisation | "∝ ε", the speed, "c_j where it is" fixed; "17–50" restated as the amplitudes' range, 1/(0.1·0.6) ≈ 17 to 50 (round 41); the margin's normalisation: (open) |

**Lesson.** The soft-analysis limit held; the claims resting on finite-ε pictures and on a discretisation tuned to the showcase targets did
not. Sweep the small parameter decades past the experiment, measure instantaneous rates rather than accumulated displacements, feed
symmetric inputs, and check that each test can tell the claimed formula from a wrong one.

## Round 38 (adversarial review of Proposition D‴.49, the rich transit in the closed loop, 2026-10-27)

**Design.** One referee, not modifying the repository. D3ak: the draft of Proposition D‴.49 and its code at b59721f
(`review/REVIEW_round38_D3ak.md`): Parts A and C reproduced from the committed raw file; the real loop re-run with an independent C
integrator (DOPRI5, Q kept exactly) against the repository's `real_loop` and every real-loop row of the md; the exact averaged loop
along the frozen periodic orbit (an FFT fixed point) against the leading order `OrbitBumps`; sweeps of the loop's own knobs (k_e from 2
to 40, ω, the amplitude, a quasi-periodic reference), bumps moved toward the range, close bumps, 1% measurement noise, Part D to
ε = 10⁻¹⁰; 30 mutants.

**Verdicts (D3ak).** The orbit filter and its 1/k_e factor right; but (i)'s averaged loop was the leading-order truncation (the orbit
taken to be the reference), which misplaces the first entry by 7% and gets the identity or count wrong in 9 of 10 settings, where the
exact frozen-orbit average reproduces the real loop; the averaging hypothesis fails on the active block (its centres relax 5–13 times
faster than the reference); the count of formed kernels depended on the knob (six at k_e = 2–3 with three expelled, five at k_e = 20,
four at ω = 1); "never formed" false within 1–2σ of the range; 1% noise expels 12 kernels; 24 of 30 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ak-1 | major ((i)) | the averaged loop is the leading-order truncation along the reference | the exact averaged loop along the frozen periodic orbit (`FrozenOrbit`, the contraction e ↦ L[r(x_m + e)] for k_e > Lip r) and a C integrator of the real loop (5baef23); (i) restated, the leading order kept as such with its additive error; Part A rerun (e7ed133): the same entrants and coefficients, the entries within 0.4 in τ |
| D3ak-2 | major ((ii), proof) | averaging fails on the active block | (ii): averaging a statement about the inactive kernels' duals, the active block's agreement numerical |
| D3ak-3 | major ((ii), (iii)) | low gain expels formed kernels | Part B sweeps k_e = 2 … 40 (e7ed133); (iii)(a): below k_e ≈ 4 six kernels, three stranded at \|x\| ≈ 1.5 with \|w\| = 0.51–0.82 |
| D3ak-4 | major ((ii)) | the count depends on where the near ties fall | (ii): one kernel per feature at masses small against the entries' margins (D‴.51); Part B: five at k_e = 20, four at ω = 1 and with a quasi-periodic reference, in both loops |
| D3ak-5 | major ((i), (iii)) | "never formed" false near the range | Part C: bumps 0.92σ, 1.6σ, 2.9σ, 3.9σ below the range; (iii)(c) |
| D3ak-6 | major (draft) | measurement noise without a leak | Part E (noise with and without D‴.46's mirror leak) and E2 (the leak's strength, e7ed133); (iii)(b): the leak holds the kernels at 1% only at a strength that halves the amplitudes — its hold and its bias tied —, and every leak fails at 2% |
| D3ak-7 | major (tests) | 24 of 30 mutants survive | `tests/test_closed_transit_loop.py` (5baef23): the C loop against the Python loop, the leak, the frozen orbit as the periodic solution, the exact averaged loop against the real loop |
| D3ak-8 | minor ((iii)) | transit time ∝ k_e only for large k_e | (iv): τ₁/k_e = 4.9 at k_e = 2, 12.7–13.9 from 4 to 40 |
| D3ak-9 | minor (D) | Part D's ε and horizon | Part D at 10⁻⁶ and 10⁻⁸: Wasserstein–Fisher–Rao seven kernels in three clusters, then three; fixed centres six |
| D3ak-10, 11 | minor | the price in tracking; "keep their starting positions" | (ii), (iv) restated; the inactive displacements reported (0.119 at 10⁻⁴) |
| D3ak-12, 13 | minor | near ties flip with the evaluation point; the steady state's statement and bound | (i): the frozen orbit's choice, not the reference's; the contraction and its factor Lip(r)/k_e |
| D3ak-14 | minor (novelty) | Wang–Hill, Boffi–Slotine, Fabri–Kadirkamanathan, RAN/MRAN/GAP-RBF uncited | Status: cited |
| D3ak-15 | minor ((ii)) | close bumps at 0.83σ merge | (open; the open loop's resolution limit is Proposition D‴.52) |
| D3ak-16 | cosmetic | tables, entry times beyond their grid | entries interpolated in the md |

**Lesson.** The draft was tested at one point of the controller's parameter space, where the leading-order picture happened to give
the right features. Compute the averaged system the theorem delivers — the frozen nonlinear orbit — not its truncation; sweep the loop's
own knobs; move the "unvisited" feature toward the range; and add the smallest realistic perturbation to a law with no leak.


## Round 39 (adversarial review of Proposition D‴.51, the transit at a finite mass, 2026-10-27)

**Design.** One referee, not modifying the repository. D3al: the draft of Proposition D‴.51 and its code at aab77ca
(`review/REVIEW_round39_D3al.md`): every number of the draft's numerics against the md and raw file (Part D's 48 predictions
reproduced bit for bit); (i) re-derived; the loop's own thresholds bisected on all 48 of Part D's targets, 16 of Part F's, single-bump
near ties, D‴.48's target 10 and Part D's target 2 with its middle bump moved; controlled loops (the runner-up's centre frozen, one or
both centres held fixed until the entry, exact ties made by unequal masses); the Part E paths' margins; 34 mutants.

**Verdicts (D3al).** (i) right, at the rate ε^{1/2} (not ε), the power law 6.7% off at ξ = 15 (not "more than 20% below 15"); (ii)'s
leading order right at a symmetric tie, but its finite-margin correction is the partners' differential approach creep (the loop's
x_θ = 6.39, 7.10, 8.19 at μ* = 0.0035, 0.034, 0.099; the Gibbs value with both centres held fixed), exact ties made by unequal masses do
not split ½, a runner-up on the entrant's near side is passed over, and a near-tie loser is mobile and migrates to an unfitted feature;
(iv) false as stated — the loop fails below the predicted mass on most targets (median ε_loop/ε_c ≈ 0.65 where an entry binds), with
outliers from ×0.13 to ×10⁵ —, which Part D's decade grid cannot see (a predictor without its entry threshold still scores 183/192);
(v) the predictor's statistics, not the loop's; 16 of 34 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3al-1 | major ((iv)) | (iv) false as stated; the decade grid hides it | (iv) restated as a leading-order estimate with its measured error; Part G: the loop's threshold bisected on D's 48 targets (ε_loop/ε_c median 0.771, quartiles 0.545, 0.917, 36 of 48 within ×2; entry-bound 0.917, 0.771, 0.545 as ε_c grows past 2·10⁻⁶ and 10⁻⁴); G2: F's five outliers bisected (×9.1 … 1.15·10⁵) |
| D3al-2 | major ((ii)) | the finite-margin correction is the partners' differential approach creep; unequal-mass ties | (ii)(a), (b); H(1) free / runner-up fixed / both fixed; HI: the inner problem's own x_θ = 6.40, 6.30, 6.27 (below the Gibbs 6.54); the loop's rise over it 0.53–0.60 times the approach-creep difference (a²b/(π²c))((g′_r/g_r)² − (g′_e/g_e)²)/x_θ; H(2): 0.885 and 0.997; tests |
| D3al-3 | major ((ii), (iii)) | a near-tie loser migrates | (ii)(c), the cubic displacement scale in (iii); H(3): the loop's threshold ∝ μ*^{3.52}, ×7–240 below (μ*/x_θ)²; the estimate's restart from the pulled margin not implemented (stated) |
| D3al-4 | major (code) | `entry_pull` integrates past the next entry (D27, ×43) | stated in (iv) and the numerics; the first version's 186/192 reported; the cap not implemented |
| D3al-5 | major (code) | the creep criterion's false alarms and misses | stated in (iv) and the numerics (G, G2); the criterion unchanged |
| D3al-6 | major ((ii)) | a near-side runner-up is passed over (E117) | (ii)(d) and its premise in (ii); H(4) with the Gibbs split's 41% at x = 0.52; test |
| D3al-7 | minor | "the gap shrinking during the slide" | (ii): the runner-up falls further behind, its share below the Gibbs split |
| D3al-8 | minor | the creep formula's mass factor b, not b^{3/2}; ρ ≠ 1 untested | Part B's formula and the status; a test of 𝒲's dependence on ρξ alone |
| D3al-9 | minor ((i)) | "20%"; the rate | (i): 6.7% at ξ = 15, 19% at 8.3, 65% at 3, the next orders 1 − 8(κξ)⁻² + 64(κξ)⁻⁴ (derived; test); the rate ε^{1/2} |
| D3al-10 | minor (B) | the reference's error; no frozen control; the 1/x² term | Part B: starts to xi0 = 3200, a three-term extrapolation; the runner-up's centre held fixed (within 0.003 of the inner problem); free − fixed against the formula 1.00–1.01 at x ≈ 347, the excess a 1/x² term growing with A |
| D3al-11 | minor (proof) | the sublevel set's zero coefficients; the plateau; the D‴.53 citation | proof of (ii): Łojasiewicz where both coefficients stay nonzero (the slide onto the feature — the threshold's regime), the large-x plateau and its recession on s ~ x⁴; the loop's convergence on bounded s-intervals |
| D3al-12 | minor ((iii)) | the first-order scale is not a threshold; the core law at ξ = 3 | (iii) restated as scales, the core law's overstatement noted; Part C's bracket (0.72, 1.43)ε_c reported as such |
| D3al-13 | minor ((v)) | the predictor's fractions; h⁴ untested | (v) restated as the estimate's statistics; h⁴ dropped |
| D3al-14 | minor | the out-of-sample history | disclosed in the numerics |
| D3al-15 | minor | κ_ij an underived heuristic moving x the wrong way | (iv): named a heuristic; HP: the estimate's own x = 6.25, 6.15, 5.85 against the loop's 6.39, 7.10, 8.19 |
| D3al-16 | minor (tests) | 16 of 34 mutants survive | six tests: the profile's next orders and ρξ scaling, the inner threshold below the Gibbs value, the approach creep (free against both fixed), unequal masses, the near side, the estimate's bias on D's target 7 (pins ε_c and x_θ's scale) |
| D3al-17 | minor (novelty) | Gissin et al., exact regularisation, the linearised Bregman ISS, AGF | the Classical paragraph and the literature note S76 |
| D3al-18 | minor | equal masses only in the code | (iv) and `recovery_threshold`'s docstring: equal masses |
| D3al-19 | cosmetic | the temperature; docstrings; Part E's wording; D39 | a/(4Gc); `recovery_threshold` cites (iv); Part B's docstring; Part E's columns named as not exclusive and the predictor's; D39's ε_c noted as outside any asymptotic regime |

**Lesson.** A threshold predictor checked at masses a decade apart agreed on 96% of the runs while biased by a factor of about two on
most targets and wrong by 10–10⁵ on others, and the same grid scored a predictor with half its mechanism removed at 183/192. The cheap,
decisive check is to bisect the loop's own threshold and report ε_loop/ε_c, beside ablations that must score worse; and a leading-order
law at fixed x is to be validated where its predictions are used, the next-order mechanism isolated by a controlled run (here: freezing
the centres before the entry) before it is folded into a predictor.


## Round 40 (adversarial review of Proposition D‴.50, the rich transit under label noise, 2026-10-27)

**Design.** One referee, not modifying the repository. D3am: the draft of Proposition D‴.50, its code at 284bde6 and its results at 0cf2e9b
(`review/REVIEW_round40_D3am.md`): the transit's response split with the scheme's tangent-linear model; the face Hessian of the designs;
design 4's small coefficient; mutants of `tangent_projectors`, `linear_response` and `gibbs_rwm`; independent noise paths per σ; Part D's
window against the face's crossing time; the convergence and the scheme's stability; Part A's small-temperature claims.

**Verdicts (D3am).** The in-window failure of the first-order response is linear and finite-ε, not noise-nonlinear (its error the same
at every σ); (ii) leaned on a Laplace expansion that is not uniform in ε where the support is smaller than the samples (design 3's face
stiffnesses 2·10⁻¹²–6·10⁻¹⁰); (ii) lacked the hypothesis T·D_i ≪ w_i² (design 4's third coefficient); `tangent_projectors` and
`linear_response` untested (8 of 8 mutants survived; `gibbs_rwm` pinned T to ±20% only); Part D a small fraction of the face's crossing,
so no test of (iii); one noise path reused across σ overstated the precision; the convergence "from every start" unproved; the scheme ULA
on an exponentially growing drift (transient in the strict sense); scope, novelty, labels and numbers overclaimed in places.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3am-1 | major | the windows' failure misread as noise-nonlinear | (iv): two windows — the noise's entry shift O(κ/\|u̇_j\|) and the entrant's finite-ε absorption O(N/(Λ·D_j\|P_{T⁻}φ_j\|²)) ∝ ε^{1/2}, an error linear in the noise; C reports both widths, the error outside the ×1, ×2, ×4 windows and inside, σ-independent as predicted |
| D3am-2 | major | (ii) from a non-uniform Laplace expansion | (i)(b): the exact identity Cov_π(Φw) = T·ΦE_π[D]Φᵀ by integration by parts at every noise level; (i)(c) Laplace with its range, not uniform in the masses; (ii) derived from the exact identity, with the remainder R_off; A: the identity on the scheme's runs (2.1–2.2·10⁻²) against Laplace's 0.55 |
| D3am-3 | major | (ii)'s hypothesis T·D_i ≪ w_i² | stated, with E_π[D] as the variance; B: E[D]/D(w_BP) and (T·D)^{1/2}/\|w\| reported (design 4's third coefficient 0.59 at σ = 0.3, where D(w_BP) misses) |
| D3am-4 | major | untested response code | `tests/test_dual_langevin.py`: the response on an entries-only path equals κP_{T(τ)}W(τ) in closed form; the Gibbs identity pins the temperature (T at 1.1T fails) |
| D3am-5 | major | D cannot test (iii) | D relabelled the local quasi-equilibrium with the window over the crossing time; (iii) formal and numerically untested (status) |
| D3am-6 | major | same-path replication | `b_job`/`d_job`: four independent paths per design and σ, means ± spread across paths |
| D3am-7 | minor | convergence from every start | (i): inf D_i > 0 on free curves ⇒ strong convexity ⇒ W₂ contraction (Bakry–Émery, synchronous coupling) at a rate ≥ min inf D_i·λ_min(ΦΦᵀ)/N, small at small masses |
| D3am-8 | minor | the scheme's stability | (i): the scheme samples π up to O(h) over practical horizons, its drift grows exponentially (Roberts–Tweedie; Hutzenthaler–Jentzen–Kloeden) |
| D3am-9 | minor | scope and novelty | free curves only; the Itô (Euler) reading not covered (Q drifts); (i) a corollary of D‴.45(iii) and (iv)(a); mirrored Langevin (Hsieh–Kavis–Rolland–Cevher) cited against mirror-Langevin diffusions |
| D3am-10 | minor | A's "O(T) shift" | the label dropped; a small-T run (σ = 0.1): Cov(w) against T·D^{1/2}PD^{1/2} to 0.016 with P of rank 4 in ℝ⁸, the mean within 0.04 standard deviations |
| D3am-11 | minor | A's Heun and Itô readings | the Heun drift ∝ h√t (no stationary law at a fixed step); the Itô rate left to the unit test |
| D3am-12 | minor | C's residual mechanism | the active block's Ornstein–Uhlenbeck component κ²Σ_i1/(2Λλ_i); C: the residual over it 0.94–1.01 |
| D3am-13 | minor | C's slope from different paths | every ε driven by one Brownian path on a common τ-grid; the slope per path 0.252 ± 0.006 |
| D3am-14 | minor | numbers | the numerics rewritten from the rerun; B's relaxation times of the support's normal modes, the windows ≥ 30 of them |
| D3am-15 | minor | exits dropped silently | `linear_response` raises on a path with exits; (iv) names the exits' first-order jump as not included |
| D3am-16 | minor | (iii)'s wording | the limit ε → 0 at fixed σ, not uniform; the O(T^{3/2}) on the small-α branch only |
| D3am-17 | cosmetic | +O(T²), dimensions, order of entries | fixed: (i)(c)'s O(T²); the entries keep their order while κ is small against their gaps |

**Lesson.** A relative error that does not change with σ cannot be a nonlinearity in the noise: split a gap between a first-order theory
and a simulation with the scheme's own tangent-linear model, size the exclusion windows by each mechanism's predicted width, and prefer
the exact identities — here integration by parts, E_π[∇F∇Fᵀ] = T·E_π[∇²F] — to leading-order expansions whose uniformity fails exactly in
the rich regime the proposition is about.

## Round 41 (adversarial review of Theorem D‴.53, the sliding limit, and Lemma D‴.53a, 2026-10-27)

**Design.** One referee, not modifying the repository. D3an: the draft of Theorem D‴.53 and its code at 7f5d1e4
(`review/REVIEW_round41_D3an.md`): the experiment rerun byte for byte; every number of the numerics against the md and raw file; the
gradient-flow identity re-derived; the single-kernel fit's critical points on the real grid; the exact profile W^ε against 𝒲 for
ε = 10⁻⁴ … 10⁻⁸; an independent integration of the loop; `entry_pull`'s convergence in xi0 and tolerance; conclusion (b) at four points
of the window; the constants of all 48 targets and the full comparison on seven more; 27 mutants.

**Verdicts (D3an).** The identity exact and the isolated entrant in the continuum correct; but the constants are not uniform — the
inner problem's far-field creep tail ∝ g′²/v³ ranges from 2.5·10² to 7·10⁶ over the targets, and the five chosen ones (largest margins)
happen to have fast entrants that barely slide —; Step 3's "bounded fundamental matrix" false (the entrant's dual depends on its centre at
O(1)), though its conclusion survives a weighted estimate; the lemma's coupled case not proved; the predictions not converged (Richardson
from xi0 = 400/800 off by up to 0.014 in s_θ); statement (a)'s integrals printed with the rate instead of ⟨r∂_cφ⟩; Step 1's relative error
O(ε^{1/2}), not O(ε); 5 of 27 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3an-1 | major | non-uniform constants hidden by the targets | statement: the constants' dependence on μ₀, λ_k, v_k, g′_k, the slides and the gaps; *The regime* (ε^{1/2}\|s_θ + σ_k\| ≪ the gap, (1 + g′_k²/v_k³)ε^{1/4} ≪ 1); Part B: the constants of all 48 targets; Part A: eight targets of the twelve eligible |
| D3an-2 | major | Step 3's bounded fundamental matrix false | Step 3: a weighted Gronwall estimate — the far field's couplings carry w_e = O(ξ_e⁻³) or 𝒲′ = O(ξ_e⁻⁴), integrable —, the creep tail's O(g′²/(v³X)) shift of the arrival, X = ε^{−1/4} |
| D3an-3 | major | the lemma's coupled case | Lemma D‴.53a(iii): small overlaps under (H2), the block's Lyapunov estimate, the entrant's energy bound, LaSalle and the implicit function theorem; (H3) checked numerically at every entry of the 48 (Part B: the end within 3.1·10⁻¹¹ and 6.4·10⁻¹³) |
| D3an-4 | major | unconverged predictions | `entry_pull`'s deviation form (ζ = ξ − the outer solution), Radau at rtol 10⁻¹⁰, second-order Richardson from xi0 = 1600, 3200, 6400 (spread ≤ 2.8·10⁻⁵ in s_θ on the eligible twelve); the loop at rtol 10⁻¹¹–10⁻¹³ |
| D3an-5 | major | (a)'s integrals | ⟨r_{k(τ′)}∂_cφ(c_j⁰)⟩ in C_j and D_j |
| D3an-6 | minor | Step 1's relative error | O(ε^{1/2}), a shift of ξ by O(ε^{1/2}) (O(ε) in u) |
| D3an-7 | minor | (d)'s εD_j and the O(ε) remainder | (d): ε^{1/2}dev + O(ε^{3/4}) proved, the remainder O(ε) observed (falling ≈ √10 per decade in units of ε^{1/2}), its parts named |
| D3an-8 | minor | a weak test | `tests/test_transit_limit.py`: target 41 at 10⁻⁹ — the sign per entry and the shifts (a per-kernel sign fails), the converged predictions' spread, (H3), the block and the residual at tight thresholds |
| D3an-9 | minor | (H3) not recorded | Part B's (H3) column (the inner problem's end against S_k) and its stationarity |
| D3an-10 | minor | implicit hypotheses | (H1) no exits stated; (H4) θ a regular value; v_k with its ρ; equal masses in the code |
| D3an-11 | minor | Step 5 and (b)'s constant | (b): C ≥ 3/(2λ_k); Step 5 carries Step 3's error in the neutral directions |
| D3an-12 | minor | (b) tested at the midpoints only | Part A: the window's start τ_k + 30ε^{1/2}log(1/ε) and end τ_{k+1} − ε^{1/4} |
| D3an-13 | minor | numerics details | the paragraph rewritten from the rerun |
| D3an-14 | minor | novelty | AGF (Kunin et al.), the sliding Frank–Wolfe and ADCG, Krupa–Szmolyan's blow-up cited |
| D3an-15 | cosmetic | code nits | `near_tie_inner`'s docstring (the 1/xi0 tail); `loop_with_duals` raises on a failed segment |

**Lesson.** A matched-asymptotics theorem is only as useful as its untracked constants are small, and the numerics decide what the reader
believes about them: the targets with the largest margins also had fast entrants that barely slide, while elsewhere the inner problem's
creep tail ∝ g′²/v³ makes the first-order correction 10²–10³ times larger. Print the constants of every target, converge the predictions
past the effect being measured, and state the regime in which the expansion is informative.

## Round 42 (adversarial review of Proposition D‴.54, the transit in ℝⁿ, 2026-10-27)

**Design.** One referee, not modifying the repository. D3ao: the draft of Proposition D‴.54 and its code at 72c4fe5
(`review/REVIEW_round42_D3ao.md`): the inner problem with the feature along the face; the edge factors in the vertex's geometry; the
selection on the samples used, with unequal masses, coarse and generic sampling; the hexagonal lattice and the ℤⁿ law's genericity; the
loop's references extrapolated in xi0; Part D's statistics; 30 mutants.

**Verdicts (D3ao).** In n dimensions the one-dimensional threshold x_θ is a function on the face: the split depends on the feature's
offset t along it (x_θ = 6.30 at t = 0, 6.03 at t = h/2), so the formation rule, E[N − 1] and P(N ≥ 2) were mis-specified; the edge
factors were posed at t = 0 and used at a vertex (t = h/2), making the reported "1.6·10⁻² failure of factorisation", three phantom boundary
points and P(N ≥ 3), P(N = 4) 8–9% high artefacts of the reference; the Voronoi selection fails at the samples' window, with unequal masses
and on coarse or generic samples; the codimension law holds only within ranges, and the ℤⁿ law is non-generic; "the proofs carry over"
overstated at a near tie (H1, H2 fail); the softmax's error bound false; Part A's reference biased by its finite start; the diagonal's gap
identity misstated; `inner_k`'s stationarity test decorative; the n-D creep law in (i) correct (verified by the referee, not by the
draft).

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ao-1 | major | the split depends on the offset along the face | statement (iii)/(iv): 𝒮(x; y), x_θ(t), E[N − 1] with the face integral; Part T: p(x; t) and x_θ(t) to t = H/2, T2 beyond it to H; the face averages x̄_θ⁽ⁿ⁾ = 6.208, 6.120, 6.033, 5.947 (n = 2 … 5), P(N ≥ 2) = 71.5ε^{1/2}; D's product rule at each position's offsets |
| D3ao-2 | major | edge factors in the wrong geometry | the edges' product and the product rule at t = H/2 (A: within 9.1·10⁻⁴ of the 4-kernel problem, 2.0·10⁻⁴ for the edges' runner-ups; B: 66 of 66); B4: the areas from the 4-kernel inner problem itself (33.5, 15.5: P(N ≥ 3) = 1112ε, P(N = 4) = 513ε) |
| D3ao-3 | major | the selection outside the continuum's interior | (ii)(a) equal masses in the continuum; (ii)(b) the power diagram (the face shifted by σ²ln(b_j/b_e)/\|c_j − c_e\|; W2: the loop forms the power rule's entrant at every position); (ii)(c) the samples — aliasing by the Poisson sum (W3: the shift to the digits shown), the window ½erfc(d/σ) (W1: the zone 2.4σ at 10⁻⁸, 2.9σ at 10⁻¹⁰; C2: a feature at (0.96, 0) formed by (0.875, 0)), separation as cross-talk ≪ x_θε^{1/2} |
| D3ao-4 | major | the codimension law's ranges | c(k) = k − 1 for k ≤ n + 1 on lattices with simple Voronoi faces, zero beyond (H: the hexagonal vertex splits 1/3 each, P(N = 3) = 2426ε, the fourth kernel's gap ∝ ε^{−1/2}); ⌈log₂k⌉ on hℤⁿ with tensor samples for 2^{−c(k)} > θ (n ≥ 7: none formed as defined); the degenerate faces split on generic samples; the regions stated as the inner problem's, non-monotone |
| D3ao-5 | major | "the proofs carry over" | (i) under (H1)–(H2); the near tie's reduction to the inner problem formal; Part C: the creep's vector law εC_j checked (to 3.6·10⁻³ … 1.2·10⁻² at 10⁻¹⁰, falling as ε^{1/2}, cosines ≥ 0.999998); C2: designed targets — a near tie, the window's edge |
| D3ao-6 | minor | the softmax's error bound | (iii)(b): the error per share O(e^{−2α_j/a}), not small at the threshold (W5: −0.097, −0.026, −0.018 at A = 0.3, 0.8, 1 against e^{−2α/a} = 0.170, 0.090, 0.076) |
| D3ao-7 | minor | coinciding kernels only at the centre | (iii)(b): the profile split exact for kernels coinciding at the feature's centre with samples symmetric about it; off the centre W4 (−2.7·10⁻⁴ … −8.9·10⁻³ at 0.07σ … 0.4σ); D‴.51's wording to follow in the documents |
| D3ao-8 | minor | Part A's reference biased | `inner_k` extrapolated in its start (xi0 = 800, 1600; the change ≤ 1.3·10⁻⁵); the loop's successive changes reported (falling 3.15 per decade) |
| D3ao-9 | minor | the diagonal's gap | x₁₁ = x₁₀ + x₀₁ − ε^{1/2}x₁₀x₀₁ exactly on tensor samples (proof (iii)(d); a test) |
| D3ao-10 | minor | Łojasiewicz gaps; decorative stationarity test | status: the inner problem's convergence numerical, the argument's gaps for k ≥ 2 listed; `inner_k`: a test on the shares' drift (share_tol) and a wall-clock budget checked inside the integration (B4: 648 of 861 settled within 300 s, the counts of the rest at the budget) |
| D3ao-11 | minor | (iii)(c)'s mechanism | rephrased: the runner-ups' gaps grow while the entrant grows; below the softmax near the tie, above it far behind at small amplitudes; the slide hardly matters |
| D3ao-12 | minor | Part D's statistics | the rule evaluated at the sampled positions; the same positions at both masses stated; the misses listed with their gaps and thresholds |
| D3ao-13 | minor | the mass for exact formation | (iv): labelled as the single-feature near ties' component only (not the creep, the cross-talk, the amplitudes' x_θ or the losers' later dynamics) |
| D3ao-14 | minor | text and md mismatches | the numerics rewritten from the rerun |
| D3ao-15 | minor | thin tests | `tests/test_transit_nd.py`: the diagonal's exact gap; the product rule and the areas on synthetic tables; the power-diagram and aliasing formulas; x_θ falling along the face and its n-D averages; test 5 at the vertex's offset |
| D3ao-16 | minor | novelty | Classical: Voronoi tessellations, power diagrams, Duval–Peyré's thin-grid Lasso (the static edge split), multiplicative weights, independent Gibbs measures, tube volumes; New: the finite-mass dynamic law |
| D3ao-17 | cosmetic | exact ties' offsets; C's threshold | `c_job`: formed against the nearest feature's amplitude |

**Beyond the report.** C2's near tie at ε = 10⁻⁹ exposed a new mechanism, (v): a near tie's runner-up left at a small tail coefficient is
an active kernel held at its feature only at O(w²) and pulled by an unfitted feature at O(w); below an unbinding level ≈ A′(d/σ)e^{−d²/(4σ²)}
it leaves and fits that feature before the feature's own entrant enters (the loop's kernels then differ from the path's, the fit and the
count per feature by position do not) — in a window of masses: bound above it, too slow below it (C2 at 10⁻⁸ and 10⁻¹⁰).

**Lesson.** A one-dimensional constant becomes a function in n dimensions: the threshold x_θ depends on where the feature sits on the
face, and a reference posed in the wrong geometry turns its own error into a "finding" (the 1.6% failure of factorisation). Pose every
reference in the geometry where it is used, converge it in its own parameters, and let designed targets (a near tie, the window's edge)
probe the claims — the designed near tie found a dynamic the random targets never showed.

## Round 43 (adversarial review of Theorem D‴.55, the transit's certificate, 2026-10-28)

**Design.** One referee, not modifying the repository. D3ap: the draft of Theorem D‴.55 and its code at 6666e70
(`review/REVIEW_round43_D3ap.md`): the md regenerated byte for byte from the raw file; every number of the numerics against it; Parts A
and C rerun with the polish off, with t_fast = 10³ and with loosened polish tolerances; baselines and alternative certificates on all 520
targets; designed families — same- and opposite-sign pairs in one and two dimensions, pairs with amplitude ratios, a lattice scan, lattice
pairs, every ordering of the amplitudes at ranges 10³ and 10⁶, K = 4, 5 and edge features; 13 mutants.

**Verdicts (D3ap).** The band |m| ≤ 0.8q false both ways (sufficiency fails at m/q = 1.057, necessity at −1.62, −0.90, −0.81); the
amplitudes' range enters outside the band (the same certificate exact at one amplitude ratio and not at another); (ii)'s accounting omits
the screening by earlier-fitted blocks; (iii) a sketch with q₀, C never given, its trapping region degenerate at R = 0; (iv)'s
"opposite signs recovered closer together" contradicted (exact with E on 0, 14, 2 of 40 at 1.5σ, 2σ, 2.5σ); the evidence outside the band
trivial ("always nearest" scores 131/135); the band constant fitted after the fact; "every path exact" dependent on the polish
(`polish_block` carries an entrant still at its seed into another feature's basin); large ranges and lattice points not computable;
the finite-mass claim amplitude-dependent; 8 of 13 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ap-1 | major | The band \|m\| ≤ 0.8q false both ways | No band is claimed: at a finite separation the outcome is the skeleton's ((ii), (iv)), the certificate its zeroth order; `part_D`'s designed families (lattice pairs with amplitude ratios 10⁻³ … 10³, pairs at 1.5σ … 4σ over a grid cell, every order of the amplitudes): the certificate's sign wrong on 24 of 332 designed targets outside 0.8q (\|m\|/q up to 1.62), 11 of 270 outside q |
| D3ap-2 | major | The amplitudes' range enters outside the band | Title and (iv): the skeleton depends on the amplitudes through the order of the entries and the phases; D(a): at 2.5σ the outcome changes with the ratio; (iii) the separated regime uniform in the amplitudes (the order by \|A_k\|G(R_k) bounds the unfitted amplitudes) |
| D3ap-3 | major | (ii) omits the screening; T degenerate at R = 0; q₀, C not given | The skeleton is the screening: r_k = f* − Π_kf*, Π_k the projection on the fitted atoms and their tangent directions ((ii), `skeleton()`); T's radius max(2R, σ/2); Δ₀ and the O(δ²) band not computed (status) |
| D3ap-4 | major | "Opposite signs recovered closer together" contradicted | Deleted; D(b): opposite-sign pairs at 1.5σ … 2.5σ exact with E on 3 of 60 where the certificate predicts all 60; the skeleton right on 240 of 240 |
| D3ap-5 | major | The margin compares \|η\| with 1 while η(e_k) = 1 + O(q) | `certificate()` adds the relative margin min_k s_kη(e_k) − max_{j∉E}\|η(c_j)\| (344/400, 102/120, 239/240 pairs), reported beside the certificate; the theorem's predictor is the skeleton's margin |
| D3ap-6 | major | The out-of-band evidence trivial; no baselines | `report()`: the single-feature margin (= always the nearest kernels: 262/400, 113/160, 83/120) and the relative margin as baselines for every part; designed targets with negative margins (D(a), D(b)) |
| D3ap-7 | major | The band constant fitted after the fact | No constant is claimed; A and C still count the old bands for the record, and D tests them on targets designed to break them |
| D3ap-8 | major | "Every path exact" depends on the polish | `sliding_iss`, `sliding_iss_nd`: a polish is accepted only if every kernel at a seed-level coefficient (below 10⁻² of the largest) lands within one kernel width of its start or of its correlation's ascent in the current residual (`ascent_target`), else the flow continues and the polish is retried after each doubling (`polish_rejected` returned); over the 1168 stored paths of D‴.51, D‴.55 and D‴.56, 4 outcomes change (all in B, artefacts of the jump, now exact with E); target 267 with t_fast = 10³ now right; tests |
| D3ap-9 | major | Large ranges and lattice points not computable | Partly: D(a) runs lattice pairs to the ratio 10³ under a 600 s limit per path, 1 of 54 over it — the opposite-sign pair at 0.3 with the ratio 50, whose second block (the large feature's kernel and a kernel sliding from −0.2 onto the small feature) is stiff at the default t_fast = 10⁵; with t_fast = 10³ it ends in 0.2 s exact with other kernels, as both skeletons predict —; the tolerances are still absolute and the range 10⁶ is not attempted |
| D3ap-10 | major | The finite-mass claim depends on the amplitudes | (iii): the mass small against (μ/x_θ)² for the path's margins at the entries, x_θ = ln(\|A_k\|/θ − 1)/κ growing with the amplitudes' range, below the creep's threshold — the mass, unlike the path, depends on the amplitudes |
| D3ap-11 | minor | (i) false near the samples' edge; "ẇ has A's sign at w = 0" | (i) in the continuum; on samples the entrant is the kernel of largest correlation; the wording: {w = 0} invariant, the entrant starts with A's sign |
| D3ap-12 | minor | The 3.3σ threshold is a midpoint-only statement | (iv): the midpoint algebra kept as where η fails, "not a bound on where the path fails" |
| D3ap-13 | minor | Part B weak | B reports the realised ranges (median 33.9, at most 749 at r = 3) and says it varies the range, not the order; D(c) varies the order (all six orders of 1, 10^{−1/2}, 10^{−1} on 30 targets: no change) |
| D3ap-14 | minor | A hard-coded sentence | Computed |
| D3ap-15 | minor | Tests weak | 8 tests in `test_transit_certificate.py` (the skeleton for one feature; the skeleton following the amplitude ratio on the lattice pair where the certificate cannot; the fitted skeleton resolving the linear one's near tie; the relative margin; `outcome()` with exits; the bounded polish on target 267; `band()`) and `test_ascent_target_climbs_to_the_residuals_nearest_peak`; the mutants not rerun |
| D3ap-16 | minor | `outcome()`'s tolerances and classes | Not changed: exact with E requires no exits by definition (a path ending on E after exits is 'other'); the designed features sit 10⁻⁶ off the starting points, so the nearest point is never tied |
| D3ap-17 | minor | The generator avoids edges and wide ranges | D's designed families (lattice points, ratios 10⁻³ … 10³, positions over a cell); (i) restricted to the continuum for the edges |
| D3ap-18 | minor | Novelty overstated | Classical: Fuchs, the irrepresentable condition, the inverse scale space's sign consistency, the homotopy k-step property, sliding Frank–Wolfe's finite termination, η_V; New: the sliding inverse scale space's separated-regime exactness uniformly in the amplitudes and its skeleton |
| D3ap-19 | minor | "Proved up to constants" and "iff" | Title without "iff"; status: (iii) proved with Δ₀ not computed, (ii) first order with the band's constant not computed, (iv) numerical |
| D3ap-20 | cosmetic | Wording; `q` reused in `part_B` | Fixed |

**Beyond the report.** The revision's predictor is the path's skeleton — the inverse scale space on the starting kernels with the fitted atoms moved onto their features and their tangent directions projected out —, right on 1148 of 1153 targets (A 397/400, B 160/160, C 120/120, D 471/473), its misses four near ties and one block beyond its linear regime; with each block at its own least-squares fit instead (the fitted skeleton) it is right on all 1153: the path's outcome is K nonlinear least-squares fits and piecewise linear duals, no flow, wherever the blocks' flows end near their features.

**Lesson.** A certificate claimed up to an unspecified O(q) band and "confirmed" with a constant fitted to the same random targets
cannot be falsified, and random targets put the out-of-band mass where "the nearest kernels win" is already right; designed families —
features on starting points, opposite-sign neighbours, every ordering of the amplitudes — break the band and the amplitude-free claim in
minutes. A shortcut in the solver (the polish) can decide the outcome it is meant to compute: audit it against the flow and bound it. And
an amplitude-free certificate cannot predict an outcome that depends on the order of the entries: the path's own skeleton — the inverse
scale space with the fitted atoms moved onto their features — can, from projections alone.

## Round 46 (adversarial review of Theorem D‴.55 as revised after round 43, 2026-10-28)

**Design.** One referee, not modifying the repository. D3as: the round-43 revision of Theorem D‴.55 and its code at 010f24c
(`review/REVIEW_round46_D3as.md`): the md regenerated byte for byte from the raw file and every number of the numerics against it;
(ii)'s order measured against the blocks' displacement δ; a single entrant integrated on the critical curve; ξ_θ against x_θ; the E map's
linear-against-fitted disagreements and their scaling, D(b)'s m₀, the decision phase of all 1153 targets, ablations of the skeleton; the
eight E timeouts and D(a)'s traced through LSODA, LSODA with the exact Jacobian and BDF/Radau; a flow-only audit (polish off, BDF with the
exact Jacobian) on 635 distinct targets; designed targets in one dimension (K = 4, 5, edge features, near ties, the amplitudes 1, 10⁻³,
10⁻⁶ in all six orders, a global scale 10⁻³, same-sign chains) and in the plane; 32 mutants; two web checks.

**Verdicts (D3as).** The revision's main empirical claim survives: with the polish off and a robust integrator the flow agrees with the
fitted skeleton on 635 of 635 targets. Wrong: (ii)'s order (the linear skeleton errs at first order in the tail, slope 0.97–1.12, not
O(δ²)); (iii)'s mechanism (the entrant slides before it grows, not after) and its cross-talk bound (a kernel between two features, or beside
e_k, receives more than q); the solver (LSODA stalls on easy targets at ratios ≥ 7.5, reported as "over the limit"); the E map's "Δ₀ for
K = 2" (the fitted skeleton's threshold for one geometry, contradicted by D(b) on the same grid); 17 of 32 mutants survive.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3as-1 | major | (ii)'s order wrong | (ii), proof, status: the block's residual is r_k + O(\|A_f\|δ²), integrated over a phase ∝ 1/\|A_next\|, so the linear skeleton's duals err by O(δ²\|A_f\|/\|A_next\|) = O(δ·t), first order in the tail; its sign decides except where \|m_sk\| ≲ Cδt; the fitted skeleton is the predictor; E measures the order: \|m_sk − m_fit\| against the amplitude ratio, log–log slope 0.98–1.00 on eight rows at 2.5σ … 4σ |
| D3as-2 | major | (iii)'s proof: the mechanism backwards, the cross-talk at the competitors' distances, no packing bound | Proof (iii) "a sketch": the seed climbs its correlation before it grows (dc/dα = α⟨r∇φ⟩/⟨rφ⟩), to the ascent's end within O(q′σ) of its feature; q′ the cross-talk at the competitors' distances (a point between two features, one beside e_k); a packing bound makes q′ independent of K; status: not proved — the ascent's end and its time against the growth, T's invariance from seed level, the uniqueness, with constants |
| D3as-3 | major | The solver stalls on easy targets | `relax_chunk` (e3b77b3): a chunk of a block's fast flow that LSODA has not finished within 2·10⁵ evaluations is redone from its start by Radau — bit-identical below the budget —, in `sliding_iss` and `sliding_iss_nd`; the full rerun has no path over its limit (D: 0 of 474; E's 240 sampled paths all within 300 s and all as the fitted skeleton); Part F: the amplitudes 1, 10⁻³, 10⁻⁶ in all six orders at 5σ and 3σ, exits only below 10⁻¹² of the largest coefficient |
| D3as-4 | major | "Δ₀ for K = 2 on this grid" is not one | (E): "the fitted skeleton's threshold D*(ℓ) on this grid … not (iii)'s Δ₀"; the near-tie threshold m₀ ≤ 0.05 stated as lying between the map's m₀ values 0.034 and 0.067; D(b)'s failures with m₀ > 0.05 reported per row (beyond 3σ only three opposite-sign pairs at 3.5σ, m₀ = 0.054–0.056); the range 10⁶ (F) with the fixed solver |
| D3as-5 | minor | (iv)'s failure domain too narrow | (iv): the linear skeleton decides except within its error scale \|m_sk\| ≲ Cδt — near ties, same-sign pairs near 2σ and at ≈ 3.8σ for ratios 0.04–0.1, opposite-sign pairs —; E: \|m_sk\| > 0.01 at 33 of its disagreements (up to 0.134) |
| D3as-6 | minor | Circular evidence for the fitted skeleton | `job_1d` reruns each path with the polish off: the flow ends as the polished path on A's 400 and D(b)'s 240, its blocks relaxed on all; the flow-free basin check (each entrant's ascent target in r_{k−1}) not added |
| D3as-7 | minor | (v)'s window claim unproved | Dropped: "between them it may change more than once (E records where)" |
| D3as-8 | minor | The safeguard bounds less than stated | md and numerics: the bound does not check single-kernel blocks or the block's largest kernel; the per-kernel scale not adopted (the unpolished flow agrees with the polished path wherever computed) |
| D3as-9 | minor | The finite-mass sentence omits ξ_θ | (iii): ε small against (μ/max(x_θ, ξ_θ))², ξ_θ ∝ (θ/\|A\|)^{−1/3} — 9.3 and 45.7 at θ/\|A\| = 10⁻⁴ and 10⁻⁶ against x_θ = 13.1 and 19.7 —, so ε ≲ μ²(min\|A\|/max\|A\|)^{2/3} at wide ranges |
| D3as-10 | minor | 17 of 32 mutants survive; the stored m_skel hides regressions | `fill_margins`: both skeletons and the atoms baseline recomputed from the targets at report time; tests: the fitted skeleton on negative features (K7), the atoms baseline, the entrants near the window's edge, the ascent's end, the outcome's tolerances; round 46's 17 surviving mutants rerun against the revised tests: 13 killed — K2, K6 and K7 on stored targets each gets wrong (`test_skeletons_on_the_targets_that_break_their_mutants`), O3 and O4 (`test_outcome_counts_kernels_and_its_tolerance_is_absolute`), the safeguard's P3, P5, P6, P7, P9 and P11 in both `sliding_iss` and `sliding_iss_nd` (`test_the_polish_safeguard_on_stored_paths`), P4 and A1 by round 46's tests —; K3 and K11 no longer apply (the code they mutate was rewritten); O6 is equivalent on consistent paths and P8 (polish_wfrac 0.5) on every stored path |
| D3as-11 | minor | Edges: title, body and code disagree | E is the starting kernels of largest correlation in the title, the statement, the code (`entrants`) and the md; (i): on samples the entrant need not be the nearest kernel near the window's edge; a test (a feature at 0.97) |
| D3as-12 | minor | Trivial cases and baselines unreported | `predictors()`: "E iff the first entrant is in E" and the atoms subtracted beside the margins in every part; numerics: over 1154 targets the fitted skeleton 1154, the linear skeleton 1149, the atoms subtracted 1100, the relative margin 1014, the first entrant 977, the certificate 903, always E 763; 214 decided at the first entry |
| D3as-13 | minor | Novelty: the adaptive inverse scale space and saddle-to-saddle not cited | "Classical": Burger–Möller–Benning–Osher's adaptive inverse scale space, whose exact recursion of low-dimensional nonnegative least-squares problems is the fixed-dictionary skeleton; Pesme–Flammarion's saddle-to-saddle limit of vanishing initialisation |
| D3as-14 | cosmetic | Numerics omissions | Polishes rejected per part (A 2 on one path, B 8 on 4, D 1 on one, E 14 on 3 of its sampled paths, none in C or F); `job_2d` records `rej` |
| D3as-15 | cosmetic | Stale cross-reference in D56 | D56 cites D‴.55's E and the fitted skeleton, not the withdrawn band |

**Beyond the report.** The rerun found a header bug of its own: D(b)'s table had gained the atoms column in its rows but not in its
header (`report()` now builds the header from the predictors). Part F's one path that is neither exact with E nor with other kernels ends
at an exact fit with four kernels (3σ, + + +, the middle feature 10⁻⁶ and the right 10⁻³), which both skeletons count as not E.

**Lesson.** A claimed order of accuracy is checked by its scaling, not by sign agreement on targets most of which are easy: the skeleton
agreed with the path on 1148 of 1153 targets while its error was first order, and 214 of those targets were decided at the first entry,
where any predictor that names the first entrant is right. A solver's shortcut is audited against a robust integrator, not against the
least-squares solve it performs itself — done so, the polish was vindicated, and the "hard" paths were the stiff integrator's stalls, not
the dynamics'. And a proof's mechanism ("grows before it moves") is checked with a two-line simulation before it carries a theorem.


## Round 47 (adversarial review of Proposition D‴.54 as revised after round 42, 2026-10-28)

**Design.** One referee, not modifying the repository. D3at: the round-42 revision of Proposition D‴.54 and its code
(`review/REVIEW_round47_D3at.md`): the md regenerated byte for byte from the raw file, C2 rerun at two masses, every number of the
numerics against the md and the raw; a static test of (v)'s binding in one dimension at a held tail coefficient (d = 4–7σ); C2's near tie
at 20 further masses and at rtol 10⁻¹¹ … 10⁻¹³; two variants with the unfitted feature at 5σ; D‴.51's H(3) target in one dimension over
30 masses and five tolerances; 105 of B4's 213 unsettled problems rerun; Part D's rule with offset-blind thresholds; the product rule's
areas on grids down to step 0.02; 24 mutants.

**Verdicts (D3at).** Most of round 42's fixes hold in the plane (x_θ(t) and the face averages, the edges posed at t = H/2, the power
diagram, the aliasing and window formulas, the hexagonal check, the softmax error per share). (v), new in the revision, fails: its balance
is inconsistent (the entrant's shift that makes the binding quartic also cancels the linear pull; the surviving pull is a curvature term,
so the unbinding level was one to two orders of magnitude too low); runner-ups that form are captured too, so (iv)'s count fails near a
same-sign unfitted feature; "the count per feature by final positions does not differ" is false (a capture tie); and (v) is D‴.51(ii)(c)'s
migrating loser, uncited. The n-D product rule used the planar offset for n ≥ 3; B4's unsettled problems ran out of doublings, not time;
the tests killed 4 of 24 mutants.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3at-1 | major | (v)'s balance wrong | (v) and its proof: with the entrant's compensating shift the first-order terms cancel and the energy depends on M₂ = w(1 + w/w_e)\|δ\|² — the binding (1/8)M₂²‖(∂²_eφ)_⊥‖² against the curvature pull −½M₂⟨r′∂²_eφ⟩ —, so the runner-up rests at \|δ\|² = σ²A′((d/σ)² − 2)e^{−d²/(4σ²)}w_e/(w(w_e + w)), dimension-free, and unbinds in a fold at w_f = c(d)ℓ(d), ℓ = A′(d/σ)²e^{−d²/(4σ²)}; Part V (the continuum, one dimension and the plane): c(d) = 4.21, 5.41, 6.64, 7.87, 9.11 at d = 5 … 7σ (A = 0.8, A′ = 0.7), the law low by 3–19% at w = 0.3, the plane's fold the line's |
| D3at-2 | major | Formed runner-ups are captured; (iv)'s count fails near same-sign unfitted features | (iv) restricted to features without a same-sign unfitted feature within capture range (every runner-up's coefficient above the fold); (ii)(c): separation at the entries is not enough; "below θ" dropped from (v); Part C3 (the near tie at 15 masses, rtol 10⁻¹²): bound at the gaps x ≤ 4.26, captured at 5.07 ≤ x ≤ 14.5, a capture tie at 17.0, stranded beyond — the bound share A·p(x; t) crosses the fold at x = 5.22, the loop capturing from 1.1 times the fold —; with the feature at 5σ (the fold above every bound share) captured at all 15 masses |
| D3at-3 | major | A capture tie | (v): arriving as the feature's entrant enters, the migrant shares it (two kernels there; C3 at x = 17.0); arriving later it is stranded, its coefficient left to decay |
| D3at-4 | minor | (v) on one capture, one mass, one tolerance | C3 scans 15 masses on two geometries at rtol 10⁻¹² (`RTOL_C3`); V the static law at five distances and in the plane; the amplitude A′ not scanned |
| D3at-5 | minor | (v) not new | (v) titled "a migrating runner-up: D‴.51(ii)(c) in ℝⁿ"; the status cites D‴.51(ii)(c), (iii) and its observations (G, target 2; H(3)) and claims the threshold (the fold) and its reach to formed runner-ups and capture ties; round 42's "a new mechanism" corrected here |
| D3at-6 | major | The n-D product rule uses the planar offset | (iii)(d), (iv): the factors p(x_i; \|t_i\|), \|t_i\| the feature's distance from edge i's axis, (n − 1)^{1/2}h/2 at a vertex of hℤⁿ; Part T3: the 8-kernel inner problem at a vertex of hℤ³ against the product at H/√2 — the count on 6 of 7 gap triples, the shares to 1.2·10⁻³ — and at H/2 (0.01–0.02 off) |
| D3at-7 | minor | (iii)–(iv) for equal masses only, unstated | (iii)–(iv) say "equal masses" (unequal ones move the layers and x_θ: W2) |
| D3at-8 | minor | B4's unsettled problems' cause misstated | B4's 213 unsettled problems rerun with the stop's reason recorded (`part_B4u`): 209 at the doubling cap, 1 at the clock, 3 settled; the counts unchanged at 213 of 213, the shares moved by at most 1.1·10⁻⁶; every unsettled share at least 2.0·10⁻⁴ from the threshold, so at their last drifts (1.0–3.7·10⁻⁷ per doubling) no count changes for at least 788 doublings — the caption (`b4u_note`) and the numerics say so; the product rule's areas on a grid of step 0.05 (33.5 and 15.3, against 15.1 at 0.25: B4's grid good to about 1.5% for N = 4) |
| D3at-9 | minor | `inner_k`'s "stationary" a cap and clock flag | `inner_k` returns why it stopped ('converged', 'cap', 'clock'), the drifts per doubling and the remaining change estimated from them, drift²/(drift_prev − drift); `test_inner_k_says_why_it_stopped` (an edge converges; a vertex's slow tail creeps, its drift growing with the doubled chunk, until the cap; the clock; each drift against separate runs) |
| D3at-10 | minor | (i) still asserts D‴.51(iv) and D‴.53(c)–(d) in ℝⁿ | (i) and the status restricted to D‴.51(i), (iii) and D‴.53(a)–(b) |
| D3at-11 | minor | Part D does not resolve x_θ(t) | report(): the offset-blind rule (x_θ(0) on every face) agrees at 596 and 600 of 600 against the offset rule's 592 and 597 — D does not resolve x_θ(t) —; the numerics name the misses' two exceptions |
| D3at-12 | minor | 20 of 24 mutants survive | report()'s formulas moved into functions (`crossing`, `vertex_counts`, `p_layers`, `e_extra`, `p_vertices`, `w_unbind`, `formed_by_feature`, `c3_outcome`), the md byte-identical from the old raw; seven tests (the product rule at unequal gaps and offsets, x_θ(t) and the face averages against a quadratic, the asymptotic probabilities and the symmetrised areas through report(), the bisection and the unbinding level, the bookkeeping by final features, inner_k's stop reasons, C3's outcomes): 22 mutants killed, the report's named survivors among them |
| D3at-13 | cosmetic | Text and md mismatches | C2's md reads "— (captured)" where the runner-up was captured; C2 and C3 judge "formed" by the start's feature's amplitude and count by final positions (stated); the vertices' tail coefficients (3.8–4.1·10⁻⁹ … 3.8·10⁻¹⁵) added; "α_j its amplitude coordinate"; A = 0.8 stated with x_θ, x̄_θ and E[N − 1]; the (4, 6)/(5, 6) non-monotonicity dropped (its shares within 0.3% of θ) |

**Beyond the report.** C3's first classification called C2's bound loops "stranded" because their residual at the end (4–7·10⁻⁸, the
finite mass's relaxation) exceeded 10⁻⁸; the outcome is now decided by the final counts first (`c3_outcome`), and the caption's distance of
the same-sign feature is computed (6.08σ, not the 5.9σ written by hand).

**Lesson.** A mechanism found by balancing two terms must expand both under the same constraint: the compensation that made the binding
quartic also cancelled the pull it was weighed against, the level came out one to two orders of magnitude too low, and the one observed
capture looked consistent with it. A threshold law is validated by scanning what moves the threshold (the mass, the distance) at a
tolerance where the saddle-borne outcomes have converged, and by reading the earlier proposition that had already recorded the phenomenon
— D‴.51's H(3) table held the counterexample to the count claim. Hand-written numbers in captions (5.9σ) drift from the computed ones.

## Round 44 (adversarial review of Proposition D‴.49 as revised after round 38, 2026-10-28)

**Design.** One referee, not modifying the repository. D3aq: the round-38 revision of Proposition D‴.49 and its code at e7ed133
(`review/REVIEW_round44_D3aq.md`). The referee:
- checked every number of (ii)–(iv) and of the numerics against the md and the raw file;
- traced the low-gain runs with dense output (k_e = 2 at ε = 10⁻⁴ … 10⁻⁷), the state's range per period and each expelled kernel's path;
- ran the threshold with f × 2 and f × ½;
- compared the real loop's period averages with `FrozenOrbit`'s drive;
- computed the active block's Jacobian;
- ran the frozen orbit below Lip f, the leading order's error from k_e = 5 to 160, the count at ε = 10⁻⁷ and 10⁻⁸;
- ran Part C with every |w| listed, the in-range bump's sign flipped and the horizon at 2·10⁴;
- ran fixed centres to τ = 3·10⁴ against the grid's least-squares fit, a feature of width 0.2;
- ran noise with the desired regressor, the mirror leak at 5·10⁻³ and 10⁻², a filtered adaptation, a dead zone and a second seed;
- ran 15 mutants and three web checks.

**Verdicts (D3aq).** The revision computes the right averaged object and sweeps the knobs. Four sentences fall, to one cheap check each:
- the low-gain mechanism: the state never leaves the range, while the kernels are ejected by the fast active block within a period;
- the tie between the leak's hold and its bias: the repository's own desired regressor breaks it;
- the near-range thresholds: a doubled horizon and a flipped sign change them;
- "a smaller mass forms fewer kernels": false at low gain, where a smaller ε makes it worse.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3aq-1 | major | (iii)(a)'s mechanism false | (iii)(a): an averaging failure of the fast active block — the formed centres swing within a period by 0.42 (k_e = 2), 0.18 (3), 0.05 (10), and a kernel at an edge feature is carried out within one period while the state stays in [−1.0275, 1.0025]; Part B2(b) traces it (kernel 17: formed at 0.70, 0.60 → 1.24 in the period before it leaves, \|w\| = 0.82, the state in [−0.995, 0.995]); B(b): the state's range and the transit-time swing for every knob; the averaged loop expels nothing |
| D3aq-2 | major | The averaging's small parameter is not ε | Proof sketch: the inactive kernels' averaging is classical; the active block needs the fast subsystem (e with the active block) to settle within a period — its rates below k_e, the rest condition, numerical. B(b): the active block's rates (the fastest 4.47k_e, 2.41k_e, 1.49k_e, 1.0k_e, 0.27k_e at k_e = 2, 3, 4, 5, 10); B2(a): 1, 3, 3, 6 kernels expelled at ε = 10⁻⁴ … 10⁻⁷ (k_e = 2), f × 2 expels at k_e = 4, f × ½ clean at 2; (iii)(a): "in this instance clean for k_e ≥ 4" |
| D3aq-3 | major | "Tied" refuted by the desired regressor | (iii)(b): the desired regressor (D‴.31, D‴.33) keeps the three kernels at 1% (the fit's rms error 0.059, the true error's 0.007) — "tied" restricted to the mirror leak —; at 2% every remedy fails (E3: the desired regressor 12 kernels, 11 out; the mirror leak at 3·10⁻³, 5·10⁻³, 10⁻²: three at the wrong places, one, none; D‴.57's dead zones 7 and 12) |
| D3aq-4 | major | (iii)(c)'s thresholds horizon- and sign-specific | (iii)(c): stated in terms of the drive along the orbit, the horizon and the neighbours' signs; Part C to τ = 2·10⁴ with the bump at 0.07 as given and flipped, the orbit's error at τ = 10³, 5·10³, 10⁴, 2·10⁴: at 1.6σ the feature forms late (5.3·10³) through a compensating cluster that settles (7.6·10⁻⁷), cleanly with the flip; beyond 3σ never by 2·10⁴ |
| D3aq-5 | major | Surviving mutants | Tests (with the rerun's code): the C loop with noise, the desired regressor and a quasi-periodic reference against the Python loop; the error bookkeeping (e_first, e_last, the true error) and the entries' interpolation; frozen_orbit_rhs's WFR and fixed branches; the tones' rms |
| D3aq-6 | minor | Part A at unequal τ; the coefficients at 10⁻⁴ | Part A: every loop to τ = 600; (ii) and the numerics: agreement claimed for the sums over coinciding kernels (−0.600, 0.750), the split a gauge the loops fix differently |
| D3aq-7 | minor | The B table | B(a): both loops to the same τ, the averaged entries, centres, coefficients and the inactive displacements; the entries' agreement over the knobs within 2.5 (k_e = 20), stated; the numerics quote the md only |
| D3aq-8 | minor | (ii) needs kernel-width features | The setting: f a sum of Gaussians of the kernels' width in (ii)–(iv) (a wider feature: D‴.58) |
| D3aq-9 | minor | "A smaller mass forms fewer kernels" fails at k_e = 2 | (iv): above the gain where the active block is slow; below it a smaller mass expels more |
| D3aq-10 | minor | "The grid's residual" a plateau | (iv), D: fixed centres to τ = 3·10⁴ form eleven kernels, the orbit's error 9.2·10⁻⁴ — the averaged drive's zero on those eleven (7.4·10⁻⁶ on the whole grid) |
| D3aq-11 | minor | The averaged loop below Lip f | B: the knobs below (i)'s hypothesis marked (k_e ≤ 4.56); (ii) restricted to k_e > Lip f |
| D3aq-12 | minor | Noise numbers misstated | E: the true error (max, rms) and the in-range fit's rms error beside the measured error; the numerics restated from them |
| D3aq-13 | minor | Swing and rates not in the md | B(b): the transit-time swing (24 windows) and the end's, the active rates |
| D3aq-14 | minor | Counts depend on the threshold | Every table: the largest sub-threshold coefficient (at most 0.0086 against the threshold 0.009) |
| D3aq-15 | minor | Novelty | Status: the frozen-system averaging is classical (Anderson et al.); Sadegh–Horowitz's desired compensation and Peterson–Narendra's dead zone cited |
| D3aq-16 | minor | Part C averaged only; no error history | C: the orbit's error at four horizons |
| D3aq-17 | cosmetic | τ₁, 0.0032, "undefined", "alike" | τ_{e,1} for the first entry; 0.0032; the quasi-periodic reference has no frozen periodic orbit (the averaged loop not computed); "alike where both exist" |

**Lesson.** A mechanism inferred from end states — "the state leaves the range, so the kernels follow" — is checked with a trace of the
quantity it names; here the state never left, and the kernels were carried out by their own fast dynamics within a period. An averaging
claim has to name its small parameter; when the failure worsens as ε shrinks, ε is not it, and the quantity to compare is the active
block's rate against the error loop's. A remedy's "tied" trade-off is a claim about every remedy, and the repository's own desired
regressor refuted it at the first try.

## Round 49 (adversarial review of Proposition D‴.54 as revised after round 47, 2026-10-28)

**Design.** One referee, not modifying the repository. D3av: the round-47 revision of Proposition D‴.54 and its code and results at
abd9f36 (`review/REVIEW_round49_D3av.md`). The referee:
- freed the held runner-up's coefficient at V's rests and ran C3's near tie at fixed gaps across the masses;
- re-derived (v)'s second-moment law and took its A′ → 0 limit;
- scanned the fold's constant against A′/A;
- re-ran B4's unsettled problems to 24 doublings;
- compared Part D's loop with the offset-blind rule;
- ran 26 single-replacement mutants and checked the citations.

**Verdicts (D3av).** The second-moment mechanism is right, but the fold is not the capture threshold. With w free the held runner-up's
coefficient decays to the fold on an ε-free inner time, so capture is a race against the unfitted feature's own entry, whose clock runs
ε^{−1/2} faster in τ. C3's window in x is an artefact of the fixed margin: at fixed x the outcome changes with ε. The leading-order law
missed the pull's projection off φ (2–4% in |δ|, exact in the A′ → 0 limit once projected). The fold's rest is 0.94–0.76σ from the
feature, not 0.3–0.4σ, and c depends on A′/A. "Stranded" is an end state at a horizon. C2 printed round 42's refuted level. B4's
788-doubling bound was invalid while the drifts grew. Part D measured the finite-ε shift, not x_θ(t). T3's count did not discriminate the
offsets. (i) cited the wrong parts of D‴.53 and omitted (H3). 18 of 26 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3av-1 | major | The fold is not the capture threshold: the held coefficient decays (the hold), and capture is a race | (v) rewritten as a race: with w held the rest |δ|² = σ²ℓ(d)w_e/(w(w_e + w)) disappears at the fold w_f = c(d, A′/A)ℓ(d); with w free the rests' energy rises with w (V: the secant dE*/dw > 0 along every branch), so the coefficient decays on the inner time s = (τ − τ_e)/ε^{1/2} to the fold and the runner-up departs at an ε-free s_d; captured iff it arrives before the unfitted feature's own entry s₂ = (τ₂ − τ_e)ε^{−1/2}; bound, capture tie, stranded at the horizon otherwise; w_f is the threshold only in the zero-hold limit. (iv): restricted to features with no same-sign feature unfitted at their entry within the race's reach, which grows as ε → 0. New Part C4: the gap held at x = 3, 4.9, 7, 10, 17.3, 24.5, ε = 10⁻⁸ … 10⁻¹²; C4b: the same at rtol 10⁻¹³ (29 of 30 outcomes agree; x = 4.9, ε = 10⁻¹¹ at rtol 10⁻¹² sits on the saddle δ = 0 and stays, at 10⁻¹¹ and 10⁻¹³ it departs at s = 3.3·10⁵ and is captured). C3 relabelled as the window at a fixed margin; "captures from 1.1 times the fold" dropped |
| D3av-2 | major | The leading-order law lacked the pull's projection off φ | (v) and its proof: P = A′π^{n/2}σ^{n−2}(d/σ)²e^{−d²/(4σ²)}/4 (the projected pull), |δ|² = σ²ℓ(d)w_e/(w(w_e + w)), no transverse pull; `v_law` and V's caption; test `test_v_law_is_the_static_rest_as_the_unfitted_feature_vanishes` (the A′ → 0 limit); V: the law low by 0.5–16% at w = 0.3 |
| D3av-3 | minor | The fold's displacement misstated | (v), V: the rest at the fold 0.94σ … 0.76σ from the feature (the leading order's 0.44σ … 0.33σ), the binding softening |
| D3av-4 | minor | c depends on A′/A | V3: c(d, A′/A) at A′/A = 0.44, 0.88, 1.75, 3.5 (4.50 → 3.76 at 5σ, 6.73 → 6.18 at 6σ, 9.12 → 9.05 at 7σ; at 3.5 and 5σ the entrant leaves for A′, off the branch); `v_job` continues from 0.99A when the tail is unbound at the grid's top and stops when the entrant leaves (V3 had crashed there) |
| D3av-5 | minor | "Stranded on the way" | (v), C3's caption, `c3_outcome`: stranded at the horizon (one kernel per feature, the fit not exact, the tail away from both features then — it can still be recalled to its own) |
| D3av-6 | minor | Wrong md columns | C3's caption: on the 5σ variant the coefficient read is the captured kernel's; C2: the fold c(d)ℓ(d) of the same-sign later features (0.017 at 6.08σ for the near tie) instead of round 42's unbinding level |
| D3av-7 | minor | The 788-doubling argument | B4x: the 213 unsettled problems continued to 24 doublings — 211 settled, 2 stopped by the clock, the counts unchanged, the shares moved by at most 1.6·10⁻⁶ after the first cap, the drifts peaking at doublings 1–13; `b4u_note` cites the continuation, `remaining` an estimate |
| D3av-8 | minor | Part D's reason | D2: the edge loop's threshold at a finite mass, shifted by ≈ 130ε^{1/2} (+0.407, +0.130, +0.013 at 10⁻⁵, 10⁻⁶, 10⁻⁸, the same at two offsets): D measures the shift, more than x_θ varies along a face; at 10⁻⁸ the loop resolves x_θ(t) |
| D3av-9 | minor | T3: the count does not discriminate; the far corners' miss is systematic | T3: three triples with a corner between x_θ(H/√2) and x_θ(H/2) (the product at H/√2 9 of 10, at H/2 6 of 10); (iii)(d): the product exact to 0.2% for the near corners, below the far ones — the frozen tails — by 0.3–7% (20–28% at shares ≈ 10⁻³) |
| D3av-10 | minor | (i)'s citations | (i): Theorem D‴.53 (a)–(b) the limit, (c)–(d) the inner layers and the rate; (H1)–(H4), (H3) proved only for small overlaps and checked numerically here |
| D3av-11 | minor | Log and draft disagree on (ii)(c) | (ii)(c): the selection of the entrant; the count needs more than separation at the entries ((iv), (v)) |
| D3av-12 | minor | Tests miss the new code | tests: `b4u_note` on a synthetic B4u with the continuation; `c3_outcome` at its residual cut; `v_law` against the A′ → 0 limit; the C3 and V report formulas on synthetic raws; `inner_k` with a clock and a zero drift; `formed_by_feature` with a kernel changing feature (round-49 commit 928e4ae) |
| D3av-13 | minor | "Its scale D‴.51(iii)'s" asserted | dropped; the status names the mechanism D‴.51(ii)(c)'s migrating loser, observed there, and claims the hold and the race |
| D3av-14 | cosmetic | Notation and wording | the codimension written q(k); "depends on δ only through M₂" along e; "its own feature keeps one kernel fewer" for formed runner-ups; the aliasing bound doubled for oblique neighbours; W5's caption |

**Lesson.** A static fold computed with one variable held is a threshold only if that variable is slaved. Here the coefficient's own
gradient is positive at every held rest, so the "rest" drifts to the fold on an inner time that does not depend on ε, and the outcome is a
race whose clock runs ε^{−1/2} faster in τ. A scan along one line of the parameter plane (a fixed margin) cannot tell a threshold in x
from such a race; holding x fixed and varying ε can, and an outcome decided near a saddle has to be checked at two tolerances. For the
second round running, two balanced terms had been expanded under different constraints. The cheapest test was the limit in which the
leading order must be exact (A′ → 0).

## Round 48 (adversarial review of Proposition D‴.58's first draft, 2026-10-28)

**Design.** One referee, not modifying the repository. D3au: Proposition D‴.58's first draft and its code and results at 44b8aa3
(`review/REVIEW_round48_D3au.md`). The referee:
- checked (i) from scratch in 80–160-digit arithmetic (the nodes by Newton, the rule's error in closed form, the defects from the nodes,
  every constant of the proof) at τ/σ = ¼ … 3, K ≤ 20;
- computed the scaled Gauss–Hermite rules' errors in closed form to K = 60 (moment matching, the half-squeeze, the natural Fock and the
  Mercer-zero scalings, the best variance);
- continued the frame to σ* < σ;
- computed P_j(K) exactly for j ≤ 10, K ≤ 40 and extended Part D at t = 0.8 to K = 1280;
- checked C(c) on ℝ after the continuum optima;
- read Part B from the raw (gradients and Hessians at the first non-optimal blocks, the objective along the segment to the optimum, LM
  from the blocks, the fast flow continued with solve_ivp's status) and reran it with the polish step, gtol and the grid's offset changed
  one at a time;
- read Part C per draw;
- ran 42 mutants and a literature search (snippets only).

**Verdicts (D3au).** (i)'s algebra is right to 30–46 digits, and the ₂F₁ conjecture is better supported than the status said. The rest
fails where it interprets numerics at the edge of double precision:
- K* and "a non-optimal stationary point beyond K*" are numerical: a floor near e_opt ≈ 10⁻⁷, the polish's acceptance radius, the
  starting grid and LSODA's silent failures in `relax_chunk` account for every K*;
- the enter/exit cycles are an exit rule that does not re-relax;
- (iv)'s oracle reading is contradicted by its own data: the tolerance's bias dominates;
- (v) is wrong: the frame continues exactly to narrower features, Derevianko–Plonka's class;
- the novelty paragraph misplaced the closest precedents;
- 19 of 42 mutants survived, 14 of 17 in the experiment's analysis layer.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3au-1 | major | K* is not a property of the dynamics | (iii), (B), the status: K*'s values and law dropped — the transit forms the optimum wherever double-precision relaxation resolves it; Part B in three variants (base, the polish step 2σ, the grid shifted by 0.05), per block whether the relaxation became stationary, per path the chunks redone, the polish's rejections and the exits: the blocks optimal to K = 4, 5, 7, 9, 12, 12 at τ/σ = ¼ … 2 in every variant but one block (τ/σ = 1, the grid shifted, K = 9, its optimum at 3.0·10⁻⁷); beyond, every block (69 of 69) passes the stationarity test but stays at 1–3·10⁻⁷ above optima below 3·10⁻⁷ (round-48 commit 2fe6a28) |
| D3au-2 | major | `relax_chunk` ignored solve_ivp's failures | `rich_transit.relax_chunk`: a failed LSODA solve is redone by Radau, then as two half chunks to depth RELAX_SPLITS; the return says how the chunk ran; `sliding_iss` and `sliding_iss_nd` count `solver_redone` and `solver_failed`, a failed chunk ending its relaxation unrelaxed; a test reproduces LSODA's failure at 3·10¹¹ from a relaxed block (1c5a996). B's paths: chunks redone 0–3 times after a stall and 0–7 after a failed solve, none failed |
| D3au-3 | major | "Beyond K* a non-optimal stationary point; the path wanders through exits and re-entries" | the sentences dropped; the KKT exit rule (`exit_kkt=True`: a kernel leaves only if its correlation stays inward with its own contribution removed) in B and C — 0–2 exits per path, no cycles (1c5a996, 2fe6a28) |
| D3au-4 | major | (iv)'s oracle reading contradicted; the norm | (iv): the stopped error ≈ (e_K²‖f*‖²_N + 2Ks²/N)^{1/2} ≤ s((1 + δ)² − 1 + 4K/N)^{1/2}, a bias up to the tolerance beside the oracle's variance, the oracle reached only when δ = O(K/N); ‖·‖_N the sample rms; Part C's paths confined to ±2σ of the window, thresholds on max\|y\|, the KKT exit rule: the stops at the predicted 6, 4, 3 in 26 of 30 draws (the others one later), the error 0.92, 0.99, 1.02 times the bias–variance prediction (median), at most 0.87 of the bound, 1.1–2.7 times the oracle |
| D3au-5 | major | (v) wrong | (v) rewritten: for σ* < σ the nodes are imaginary in conjugate pairs and the K-term sum is a Gaussian of the kernels' width times a cosine sum at Hermite-zero frequencies (Derevianko–Plonka's class), its error the series of (i) at \|t\|; `narrow_rule`, `narrow_error`, Part A4: at σ*/σ = 0.5 and 0.8 the error is the series to 3·10⁻¹⁴ and 1·10⁻¹⁰ (K ≤ 10), its ratio per kernel 2.92 and 8.77 against 3 and 9; a test (532ef9e) |
| D3au-6 | major | Novelty and literature | the status: Derevianko–Plonka as (v)'s construction, Plonka–Riebe–Cuyt 2026, Kuo–Woźniakowski 2012, Kuo–Sloan–Woźniakowski 2017, Tang 1993, Steinwart–Hush–Scovel, Fasshauer–McCourt, Janszky–Vinogradov 1990, Szabó et al. 1996, Morozov's discrepancy principle for the noise trade-off; (ii) and Part A3 (`scaled_rule`, `scalings`): the natural Fock scaling falls at λ/(2 − λ) and the Mercer-zero scaling at q per kernel, both slower than t (63.50, 17.08, 5.034, 2.021 and 34.02, 9.916, 3.751, 2.021 against the half-squeeze rule's 63.58, 17.49, 5.680, 2.551 at τ/σ = ¼, ½, 1, 2); "stalls beyond K*" dropped (cd1efb7) |
| D3au-7 | major | The analysis layer untested | tests (36b4cf1): `optimal_fit` at odd K and unsymmetric; `ls_amplitudes` against the Gram solve where a 10⁻³ cut-off changes everything; `gh_defects` at K = 25–40; `sigma_star`, `blocks`, `sample_fit`; the stop's 2K/N; Part D's term count; `lead_coeffs`' degree check on a fake table; the report on synthetic results (A2's coefficients, B's threshold, C's stop at (1 + δ)s, D's Catalan verdict): 15 of the 17 experiment mutants killed (the two survivors in `c_job`, which a unit test does not run) |
| D3au-8 | minor | δ_t and the nodes misquoted; δ_t treated as constant | (iii), (A): δ_t(K) 0.078% → 0.040%, 0.73% → 0.40%, 3.5% → 2.1% at τ/σ = ½, 1, 2 from K = 1 to 8–10; the optimum's nodes within 0.008% (¼) – 5.2% (2) of the rule's; A3's best single variance (`best_scale`) recovers 96–100% of the gap |
| D3au-9 | minor | Proof gaps | proof of (i): the termwise integration dominated for small t (E\|He_n(√2Z)\| ≤ Cⁿ(n!)^{1/2}) and extended by analytic continuation; Uspensky's theorem for the product Gauss–Hermite rule; the rate exact as a K-th root (status) |
| D3au-10 | minor | (ii) unconditional; "2^K" | (ii): the order of the rates proved for λ and λ/(2 − λ), the rates numerical (A3); the gain ≈ 0.6K^{−1/4}(λ/t)^K, a factor up to 2 per kernel |
| D3au-11 | minor | (iii)'s wording | C ≥ 0 stated after an optimal block; "moving out with K" and "alternating flanks" dropped — the entrants come mostly from the flanks (126 of 128 within the optimal ranges; Beyond the report) |
| D3au-12 | minor | "Optimal" not certified | "the best local optimum found" ((iii), (B)) |
| D3au-13 | cosmetic | Wording and numbers | integer-valued polynomials; 10⁻¹⁵ … 3·10⁻⁵; the trapezoid–series difference ~10⁻¹⁷ (10⁻¹³ the cut-off); the status's P_j evidence (j ≤ 10, K ≤ 40, the j-th differences −(2j − 3)!!) and the referee's extension of Part D (K = 640, 1280 at t = 0.8) |

**Beyond the report.** Checking (iii)'s flank rule against the rerun: within the optimal ranges 126 of 128 entrants start beyond the
block's outermost node; the two others, at τ/σ = 2 with the grid shifted, start between two nodes (wide nodes, a kernel of the grid
between them first to its wall), and the blocks they join are still optimal. (iii) now says "mostly from the flanks". Part C's stop
happens at the predicted K in 26 of 30 draws, the others one kernel later.

**Lesson.** A number reported as a property of a dynamics — here the first block the transit misses — has to be invariant under every
numerical setting before it is interpreted: the precision, the tolerances, the acceptance rules, the solver's success flags and the
starting grid each moved it. A solver loop that reuses its last state must check the solve's status, or a failure looks like a converged
chunk. And a novelty paragraph should first try the result's one-line continuation (σ* < σ), which here found the closest precedent.

## Round 53 (adversarial review of Proposition D‴.54 as revised after round 49, 2026-10-28)

**Design.** One referee, not modifying the repository. D3az: the round-49 revision of Proposition D‴.54 — (v) rewritten as a race —
and its code and results at 886793f (`review/REVIEW_round53_D3az.md`). The referee:
- solved the ε-free reduced problem (C4's two kernels after feature 0's entry, features 1 and 2 in the residual) and compared it with
  C3's and C4's outcomes;
- mirrored C2's near tie across the face and ran the reduced problem and the loop on it;
- computed the far branch of V's static rests and the fold in C4's own geometry;
- traced the rtol-10⁻¹² run that C4b disagreed with, and scanned the side against the tolerance;
- ran 25 single-replacement mutants and checked the literature.

**Verdicts (D3az).** The projected second-moment law is right and the race is real: the ε-free problem reproduces C4's departures and
predicts C3's and C4's outcomes. But the hold is one-sided. At the second-moment order the rests' energy does not depend on the
coefficient at all; its w-dependence enters at the third moment, odd in the displacement, so the branch away from the unfitted feature
is the mirror image of the branch toward it: there the coefficient grows and the rests never fold. The inner layer collapses the pair
and keeps its direction, so the side is set by the start — but for formed runner-ups the collapse passes below double precision and the
integration's error picks the side; C4b's one disagreement was such a run, held on the far side, not a saddle. C4's "stranded" tails
never left their feature; V3's (5σ, 3.5) cell was an artefact of the continuation's start; the fold divided by was not the
configuration's own; several numbers were hand-computed; 13 of 25 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3az-1 | major | The hold is one-sided; the race misses the far side | (v) rewritten with both sides: the rests on either side of the feature, E* − E₀ = −P²/(2B) at the second order (the same on both sides, no w), the hold odd — the third moment M₃ = w(1 − w²/w_e²)δ³ adds E*₃ = ∓(1/6)P₃(σ²ℓ)^{3/2}(1 − w/w_e)(w(1 + w/w_e))^{−1/2}, P₃ = A′√π(d/σ)³e^{−d²/(4σ²)}/(8σ²) — toward the feature the coefficient decays to the fold and races; away from it the rests continue without a fold and the coefficient grows; the side (c_r⁰ − c_e⁰)·e kept by the inner layer; on the far side a tail is formed iff its coefficient passes θ\|A\| before s₂ (s_θ < s₂: one kernel more than (iv)'s count). Proof: the third moment, its projection off span{φ, ∇φ}, the envelope theorem, the gradient flow on the rests. (iv): the count moves by one either way (−1 a formed runner-up captured, +1 at the destination for a capture tie, +1 a far-side tail). D‴.51(ii)(c) in PARADIGM.md qualified (the loser migrates when its offset from the entrant points toward the feature; otherwise it is held and grows). Code: `e3_law`, `e3_ratios` (V's and VF's secants against the law: toward 1.042, away 0.963 at 7σ, mean 1.0025), VF (the far branches: no fold to w = 10⁻⁴, secants < 0 at 48 of 48 per d, 12 of 12 in the plane), FC (the folds in C4's own geometry), R (the ε-free race on both targets), C4M (the loop on the mirrored target: at (7, 10⁻¹³) the tail held far grows to 0.0224 and feature 0 ends with two kernels; at (7, 10⁻¹⁰) it stays below θ\|A\|) |
| D3az-2 | major | Rounding decides the side for formed runner-ups; C4b misdiagnosed | `race_job` and `c4_job` record the side (where the runner-up first passes SIDE_R = 0.05σ: a formed runner-up rests 0.13–0.20σ out), the pair's closest approach in the layer and the displacement along e at the end and at τ₂; C4's table prints the side for both tolerances (at rtol 10⁻¹² every run went toward but (4.9, 10⁻¹¹), far — C4b's disagreement —; at 10⁻¹³ all 30 toward); R's scan of the side against rtol (not toward at 6, 6, 3, 1, 0 of 8 gaps at 10⁻⁹ … 10⁻¹³, two merged to the last bit); (v) and the numerics: where the collapse passes below the centres' resolution the integration's error picks the side; C4b's disagreement is the far side, not the saddle; the departure of a formed runner-up includes the saddle exit from the offset, λ⁻¹ln(rest/offset) |
| D3az-3 | minor | "Stranded … away from both features" wrong for C4 | `c3_outcome`: by the recorded distance — not departed (never 2σ from feature 0), stranded (left, at neither feature), unsettled (no distance recorded: C3, run before the record); `c2_job` records the runner-up's largest and final distance; (v): a tail not departed by s₂ stays near its feature; C2's ε = 10⁻¹⁰ tail: not departed by R's s_d (no hand claim) |
| D3az-4 | minor | V3's (5σ, 3.5) cell an artefact | `v_job` continues from a bound state found by fresh starts (the entrant from −0.5σ … σ): the cell is 3.135 (fold at w = 0.42369), the test corrected |
| D3az-5 | minor | The fold divided by is not the configuration's own | FC: the fold in C4's geometry, 0.01668 (3.0% below the line's 0.0172, feature 1's binding); C3 and C4 divide by it (`fc_fold`); the departures' coefficients 0.97–1.03 of it |
| D3az-6 | minor | The race's statement underspecified | (v): the capture set {x: s_d(x) < s₂}, an interval about the fastest departure (x ≈ 5.5) widening on both sides as ε → 0; "the fold would be its edge only if the coefficient did not decay before s₂"; R tabulates s_d(x) from the ε-free problem at 23 gaps (rtol 10⁻¹³) and checks the rule on 74 loop outcomes (73 agree, one capture tie) |
| D3az-7 | minor | Hand-computed or unrecorded numbers | the numerics quote the md only: T3's table prints the formed corners' and the frozen tails' relative errors (the tails 6.5–50% below the product), C4's ratios to the fold with three significant digits (no "0.00"), the secant count over distinct branches (223 steps, V3's repeated column counted once: `toward_branches`); "(and 10⁻¹¹)" dropped; W4's range −1.8·10⁻⁴ … −8.9·10⁻³; B4's grid error named the product rule's coarse-grid error |
| D3az-8 | minor | (iv)'s sentence; the destination's count | (iv) rewritten as two sentences: the count's domain, then (v)'s three exceptions with the destination's |
| D3az-9 | minor | Tests miss the round-49 code | tests: `t3_product` with the complements; D2's shift at the offset; C3's fold column; the secant line over V and V3 and each branch once; the third-moment law; `race_job` on both sides and for a formed runner-up (the side at SIDE_R); `race_rule` and R's report; `c4_job` recording the side on both targets; `inner_k`'s clock after a chunk and a zero drift (monkeypatched); T3's frozen-tail column |
| D3az-10 | minor | "The binding softens" unsupported | dropped: the leading order's error is attributed to the higher moments; the third moment carries the hold |
| D3az-11 | minor | The novelty's neighbours | Status: Chizat's conic particle gradient descent (a Wasserstein–Fisher–Rao particle flow for sparse deconvolution) and Chizat–Bach; Boursier–Flammarion's early alignment and spurious stationary points from small initialisation, Pesme–Flammarion's and Kunin et al.'s dormant units; Poon–Peyré's multi-dimensional super-resolution; "New" limited to the one-sided race |
| D3az-12 | cosmetic | Dead code; formats | `w_unbind` removed with its test; A′/A printed 0.437, 0.875; the C4 table's side column makes "max w held" on the far side the grown coefficient |

**Lesson.** A sign computed on one branch does not settle a drift whose sign comes from a broken symmetry: at the second-moment order the
held runner-up's energy does not depend on its coefficient, so the hold lives at the third moment, odd in the displacement, and the
mirrored geometry reverses it — a held tail that grows instead of leaving. The same symmetry sets a numerical trap: the inner layer
collapses the pair below double precision before the instability chooses the side, so a side computed near that floor is a draw of the
integration's error. Record the quantity that decides (the side, the pair's closest approach), compare tolerances where that quantity
is resolved, and solve the ε-free problem, which reproduces the loop's departures in seconds, before scanning masses.


## Round 55 (adversarial review of Proposition D‴.52, close features, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bb: the first draft of Proposition D‴.52 with its code and results at fa25a4f
(`review/REVIEW_round55_D3bb.md`). The referee:
- reran cells of A, B, C and D and compared them with the raw;
- rebuilt (i)'s prediction from the exact one-bump residual and to the next order (T₃);
- re-derived (ii)'s determinant and spectra, and scaled the amplitude mobilities;
- computed the projected noise per kernel and identified which kernel enters second under H0 in all 8000 draws;
- reran C's (Δ, s) = (0.025, 0.01) cell at five centroids within a grid cell;
- computed the Neyman–Pearson bound per cell and the Cramér–Rao bound per placement cell;
- bootstrapped Δ₅₀ and its slope;
- ran 41 single-replacement mutants and checked the literature.

**Verdicts (D3bb).** The deterministic core holds and is stronger than the draft shows. The path is the exact-residual law, the third
moment explains the opposite-sign pair's slow convergence, the counting statistic's power sits at the Neyman–Pearson bound, and the
placement errors are at the Cramér–Rao bound. But (iii) was false as written. The noise reaches the duals through the projection off
the first fit's tangent span, not as s²⟨φ_j²⟩/N. The H0 delay is set over the whole dictionary, mostly by far kernels. And the test's
threshold was the H0 quantile at the true offset — an oracle null whose quantile swings 13× within one grid cell. Elsewhere: (i)'s
error term missed the (Δ/σ)² part, (i) had no no-tie hypothesis, r_sep's "(1 + O(Δ))" fails below Δ ≈ 0.003, the saturation in D was
mis-described, several numbers were misquoted, the knot-test literature was missing, and 24 of 41 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bb-1 | major | The test calibrated at the true offset; the null depends on the grid offset | Part E (new): five centroids across a grid cell × three cells near Δ₅₀, 1000 H0 and 400 H1 draws each; the H0 quantile against the law's at the truth; a test calibrated per draw by bootstrapping the law at the draw's own first fit (size and power at every offset), the oracle-calibrated test, the fixed threshold of the old offset, the Neyman–Pearson bound; C's caption names its oracle-calibrated simple null. Outcome: the H0 quantile × s 7.4–8.2 at c̄ ≤ 0.032 and 0.46–0.81 at 0.049, the law 0.891–1.080 of the draws'; the fixed threshold's false alarms 0.214–0.454 beyond its offset, the calibrated test's 0.019–0.062 in all 15 cells, its power within 0.05 of the oracle's in 12 of 15 |
| D3bb-2 | major | (iii)'s noise law and H0 mechanism false | `knot_of`/`knot_law`: the second knot to first order in the noise — the first entry by the duals, the rates ⟨Rξ, φ_j⟩ to the walls over the whole dictionary (the path's knot to 10⁻⁸ at s = 10⁻³, error ∝ s); (iii) rewritten with the projected rates, the law, its offset dependence and the matched-filter reason for the near-Neyman–Pearson power; C prints the H0 second entrants next to the first (16.6% overall, 70.2% of each cell's earliest 5%) |
| D3bb-3 | major | Knot tests prior art; Liu–Zhang misattributed; references missing | Status: Azaïs–De Castro–Mourareau, the covariance and spacing tests, Liu–Zhang (Liu–He–Ammari for positive sources), Smith 2005, Demanet–Nguyen, Li–Liao, the confluent-Vandermonde references, sliding Frank–Wolfe and ADCG; "New" limited to the delay law with its offset dependence, the first-order null law and its calibration, the near-Neyman–Pearson power, the mobility-weighted rates |
| D3bb-4 | minor | (i)'s error term | (1 + O(ϑ) + O((Δ/σ)²(1 + \|A₁A₂\|/A²))); `predicted_delay(order=3)` with the T₃ term — A's new column: (0.8, −0.4) 1.053, 1.012, 1.003 |
| D3bb-5 | minor | No no-tie hypothesis; the constant's offset dependence | (i) assumes a unique first entrant and a margin bounded away from a tie; E prints the noiseless Δτ₂Δ² at each offset |
| D3bb-6 | minor | r_sep's "(1 + O(Δ))" wrong | `sep_regimes`: the separation's rate from the exact −D·JᵀJ over the frozen and the relaxed formulas at Δ = 0.025 … 0.0005; (ii) states both regimes |
| D3bb-7 | minor | D's saturation mis-described; the Cramér–Rao constants unstated | D's table: the errors over their Cramér–Rao bounds, the blocks with both kernels near their features, the adjacent second entrants; the constants 0.0070, 0.0068, the split lost at s ≈ 14.8Δ³; the outlier's coefficients printed |
| D3bb-8 | minor | Misquoted numbers | 1.4–3.2%, 1.84–1.93, 2.3–4.7%, 0.006–0.012; the outlier w = (−698.6, 700); A's ratio path/prediction; Δ = 0.2's different kernels |
| D3bb-9 | minor | 24 of 41 mutants survive | tests: the opposite-sign entrant at the other wall with the T₃ term, the split's centre components, `detection` on a synthetic table, the knot law against the path, `one_fit`'s stationarity, the delay runs (bounded, noisy, distinct hypotheses, Δτ₂ = τ₂ − τ₁), placed blocks sorted and relaxed, the report's cells, the separation's two regimes |
| D3bb-10 | minor | No statistical precision | C: bootstrap 95% intervals for Δ₅₀ and its slope (2000 resamples); E: 1000 H0 draws per cell |
| D3bb-11 | cosmetic | Notation; ℓ; classical facts; eigenvector components | the proof's θ = w₂/W, W, δ defined; ℓ = Nh (the constants to 10⁻⁴); the orders called classical; the end at the exact fit marked numerical; B prints the split eigenvector's centre components (Δ/(2w_k) to 0.14%); the `__main__` block removed |

**Lesson.** A statistic's null must be calibrated over the nuisance parameters the experiment held fixed — here the feature's offset
from the dictionary's grid. The cure was the law itself: the second knot to first order in the noise is a minimum over the dictionary of
margins over projected Gaussian rates, which reproduces the path's knot to 10⁻⁸ and gives a per-draw calibration by bootstrapping it at
the data's own first fit. The same exercise showed the statistic's strength: its power at the true offset sits at the Neyman–Pearson
bound, because the runner-up's projected kernel is nearly the matched filter of the pair's deviation from one bump.

## Round 54 (adversarial review of Proposition D‴.58 as revised after round 48, 2026-10-28)

**Design.** One referee, not modifying the repository. D3ba: Proposition D‴.58's round-48 revision with its code and results at 265a133
(`review/REVIEW_round54_D3ba.md`). The referee:
- checked (i) and (v) from scratch at 80 digits against closed-form Gaussian integrals (to 10⁻⁵²–10⁻⁸⁰), and (ii)'s rates to K = 40;
- regenerated the md (byte for byte) and checked every quoted number of (A)–(D);
- instrumented the polish-2σ paths (the polish's status and evaluations, the flow continued, variable projection and LM from each stalled
  block) and reran B with the polish's budget raised;
- rebuilt the duals from the raw and read every entry;
- exercised `sliding_iss`'s failure branch by monkeypatching the solver;
- ran 56 mutants and searched the literature (snippets only).

**Verdicts (D3ba).** (i), (v) and (ii)'s rates hold; C(c) ≥ 0 with double zeros holds after all 147 optimal blocks. But the account of
where the transit's optimal ranges end was refuted: they end where the block's polish spends its 5000 evaluations, not at a
double-precision floor — with 5·10⁴ the ranges grow by two kernels. The flank-entry mechanism was misattributed (the entrant is the
kernel of the least time to the wall, never the largest current rate). The round-48 failure branch of `sliding_iss` and
`sliding_iss_nd` crashed or reused stale correlations, 38 of 56 mutants survived, and several numbers were misquoted.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ba-1 | major | The ranges' end is the polish's budget, not double precision | `sliding_iss(polish_nfev=, repolish=)`; the polish's status, evaluations and acceptance recorded per block (`polish_log`, also in `sliding_iss_nd`); B rerun in five variants (base; polish 2σ; grid +0.05; 5·10⁴ evaluations; 5·10⁴ repeated up to three times); B's table prints the first non-optimal block's polishes and its block/optimum; Part B2: variable projection from every non-optimal block. Outcome: the ranges K ≤ 4, 7, 9, 11 (5·10⁴) and 4, 7, 9, 12 (repeated) at τ/σ ≤ 1, optima down to 1.6·10⁻⁹; variable projection from the first non-optimal block of every path reaches the optimum (19 of 19) — the repeated polish stops on its gradient test 1.012, 2.677, 1.121 times above it at τ/σ = ¼, ½, ¾ —, from 13 of the 83 later blocks. (iii), (B), status: "double precision" and "ε_mach" dropped; the block relaxes to the optimum as far as its relaxation is carried |
| D3ba-2 | major | The entry mechanism | `sliding_iss` records each entry (the entrant's rate and the largest, its dual at the block's start, the time to the wall); B prints the entries (the largest current rate in 68 of 390; 339 of 360 entrants beyond the block's outermost node, a median 1.40σ; 21 between its nodes with earlier-charged duals); (iii): the next entrant is the kernel of the least time to the wall, (1 − \|u_j\|)/\|C(c_j)\| |
| D3ba-3 | major | The failure branch of `sliding_iss` / `sliding_iss_nd` | `relax_chunk`: a failed chunk is halved (LSODA, then Radau) down to a quarter before it counts as failed; the iteration's correlations are computed after the chunk loop in every case; the counts (stalls, redone, failed) recorded; tests with a fake solver (halving, Radau, never, partial failure), the counts and the exit rule on a relaxed block |
| D3ba-4 | major | 38 of 56 mutants survive | tests: the Fock and Mercer-zero closed forms and the best variance, A3/A4 and B's and C's cells on synthetic results, the stop's law, the entries line, B2's line and job; the exit rule wired into the path |
| D3ba-5 | minor | "96–100%" | 95.7–100%; A3 flags an optimum the retuned rule beats by more than 10⁻⁹ ("optimum unconverged": τ/σ = ½, K = 8) |
| D3ba-6 | minor | (iv) deterministic; e_K misattributed | `stop_law`: the rms's expectation s(b_K² + 1 − 2K/N)^{1/2} and scatter, P(stop at K); (iv) states the stop's law; C prints K(s), P(stop at K) and the expected against the observed stops |
| D3ba-7 | minor | (C)'s gap and entries | C prints where the largest gap falls against the stop; "to 11 entries" |
| D3ba-8 | minor | "q numerically" | q > t proved (the increasing (x − 1)/(x + 1) at β² against σ*/σ) |
| D3ba-9 | minor | Proof gaps | the atoms' bound E He_n(√2Z)² ≤ 4ⁿn!((3/4)ⁿ + (1/4)ⁿ)/2; Cramér's \|He_n(x)\| ≤ κ(n!)^{1/2}e^{x²/4} with the two-dimensional Mehler kernel at a = ϑ/(1 + ϑ) < ½; (v)'s measure on the imaginary axis (Mehler at ρ = \|t\|) |
| D3ba-10 | minor | Literature | Derevianko–Plonka (open: whether their frequencies are (v)'s), Plonka–Riebe–Cuyt, Kuo–Sloan–Woźniakowski, Karvonen–Oates–Girolami; the classical-rule discrepancy stated |
| D3ba-11 | cosmetic | Captions and wording | A's caption; "passed the relaxation's test" for "stationary"; the trapezoid-to-series difference attributed to rounding of the nodes, weights or residual |

**Lesson.** A plateau in a relaxed error is a claim about the arithmetic only after the solver's own limits are ruled out. Here every
range ended where the polish spent its budget or stopped on its gradient test, and a better-conditioned solver started from the very
same block — variable projection, which eliminates the amplitudes — reached the optimum every time. Record the inner solver's status
next to the result it produced, and restart a stall with a different method before naming its cause.

## Round 57 (adversarial review of Proposition D‴.59's first draft, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bd: Proposition D‴.59's first draft with its code and results at d0ad068
(`review/REVIEW_round57_D3bd.md`). The referee:
- regenerated the md (byte for byte) and reran cells of A, D and G (bit-identical);
- re-derived (i) and the dipole asymptotics, measured the moment, |Δc|, |Δσ|, the residual's excess and the coefficient's mobility on
  all six stalled blocks, and closed the runaway law's prefactor;
- instrumented C's paths at their budget stops and followed G's loops to τ = 1600 and at a second mass;
- built one-sign close targets and perturbations of them;
- refitted C's and D's pairs with equal and split widths;
- ran 55 mutants and searched the literature (snippets).

**Verdicts (D3bd).** (i) holds in the continuum; the dipole is the six budget stops, and its asymptotics are stronger than claimed
(v³ = 3μ₀K·s to 4–5 digits). But "dipoles need opposite-sign kernels: a target of one sign cannot make them" is false (two close
positive features give negative kernels); the loop's dipole does not freeze into an exact cancelling pair — it migrates, landing at an
ε-dependent time and feature; (ii)(a)'s bookkeeping misses the width split that the data use; (iii) omits its exceptions; the law's
prefactor was off by 2; (i)'s scope was overstated; C's stalls were reported at the wrong state; 25 of 55 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bd-1 | major | One-sign targets do make negative kernels; "separated" undefined | Part H: two positive features at ±0.08 (σ* 0.04, 0.06) and six perturbations — cancellations at σ* = 0.06 (6–7 negative kernels at the end, the fit exact to 4–10·10⁻¹⁰), a coinciding pair at 0.04, exact with the second amplitude 0.6; (ii) states the outcomes as observed on C's distribution (Δ ≥ 0.4 or 0.6, widths in [0.05, 0.3]) and adds the cancellation (b) |
| D3bd-2 | major | The loop's dipole is transient, not a frozen exact pair; PAIR_TOL ten times too wide | G's loop at ε = 10⁻⁶, 10⁻⁸ to τ = 400, 1600, 6400 with the closest pair's distance, coefficients and offset from the nearest feature; PAIR_TOL = 0.01; (ii)(c): the dipole migrates, separation and residual falling, and lands on an ε-dependent feature (target 9 at both masses, 27 at 10⁻⁸; the others in transit at 6400) |
| D3bd-3 | major | (ii)(a)'s bookkeeping | the width split σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ) cancels the third moment (`width_split`; the proof via ∂³_cφ = σ⁻¹∂_c∂_σφ − σ⁻²∂_cφ); C and D print Δc and the split against the prediction (to 4–5 digits); a test of the order (O(Δc³) → O(Δc⁴)) |
| D3bd-4 | major | (iii)'s exceptions | E: the per-draw oracle fit, the rms of the errors, the draws within 10⁻⁴ of the oracle fit (18 of 19), the unstopped dipole; F: the loop's closest pair and residual at 1.3 and 6 times the last entry — target 2's pair stays apart; (iii) restated |
| D3bd-5 | minor | The law's prefactor and κ | `dipole_law`: E_inf by the optimal dipole fit, K = v²(E − E_inf); G prints K, Δv³/Δs against 3μ₀K, the exponents of |w|, the distance, |Δc|, |Δσ| and the residual's excess, and the moment's σ-share; proof (c): dv/ds = −¼(μ_a + μ_b)E′(v), K = −⟨r_dip, ∇³φ:p³⟩/24 |
| D3bd-6 | minor | The pairs' mechanism; "a gauge" | (ii)(a): a kernel of an earlier block collapsing onto a fitted feature (C), the near tie only on the symmetric target (D, d = 0); "a gauge at coincidence" |
| D3bd-7 | minor | (i)'s scope | title and (i): one Gaussian feature, in the continuum; on samples J_N's limits fall below ½⟨f*²⟩, the samples checked through the paths; B labelled the continuum's |
| D3bd-8 | minor | Proof gaps and wording | (i): E ≥ J, the sublevel set, |w| bounded away from 0 and ∞; the paradigm's metric stated correctly (Gm against ∂_θE); "where the flow approaches that infimum it runs off"; the higher confluent terms (de Boor) |
| D3bd-9 | minor | C's stalls at the wrong state; circular classification | `sliding_iss_nd` records the stalled block (`stalled`) and the polishes (`polish_log`); G reports the stalled block's |w|, distance, residual against the dipole infimum and its polishes; `outcome` classifies by the stalled block (dipole run-off, stopped by the budget, cancellation, not fitted) |
| D3bd-10 | minor | 25 of 55 mutants survive | tests: `n_formed`, `outcome` at its thresholds, `pairs_of` at PAIR_TOL ± ε, `g_job`'s report, the width split's order, `dipole_law` on a synthetic dipole |
| D3bd-11 | minor | Literature and novelty | de Boor, Kammler, Petersen–Raslan–Voigtlaender, Marr–Hildreth, HyperBF, CL-OMPR, ADCG, Chizat, the mixture degeneracies; Duval–Peyré re-cited for the pair; "New" narrowed to the runaway law, the width split, the cancellations and the ε-dependent landing |
| D3bd-12 | cosmetic | Misquotes and docstrings | "(planned)" removed; the convergence from a start of the feature's sign; "formed kernels" in the md's header; the budget comment; 1.8–6.0, 0.035–0.072, "the 4 kernels of the entry cap" |

**Lesson.** A runaway that the block's flow follows to the end of its budget is only half the story when the loop around it has a
finite mass: there the pair keeps moving, and where it lands — and when — depends on ε. The same holds for a symmetry argument about
signs: a one-sign target is not a one-sign fit once features sit closer than the kernels' starting width. Follow the configuration past
the horizon chosen for the expected outcome, and test the sign rule on the targets built to break it.

## Round 58 (adversarial review of Proposition D‴.54 as revised after round 53, 2026-10-28)

**Design.** One referee, not modifying the repository. D3be: Proposition D‴.54's round-53 revision with its code and results at edb5dd3
(`review/REVIEW_round58_D3be.md`). The referee:
- regenerated the md (byte for byte) and reran four cells against the raw (bitwise);
- re-derived (v)'s moments, projections and third-moment law and checked them by quadrature at 5, 6, 7σ;
- continued V's far branch upward through w_e, R's mirrored problem to s = 10¹¹, and the loop on the mirror at ε = 10⁻¹⁷;
- traced c_r − c_e through the collapse at two tolerances on C2's target, its mirror and a y-reflection, with the full variational
  equation;
- ran c4_job at the far rule's boundary;
- ran 29 mutants and searched the literature.

**Verdicts (D3be).** The third-moment law, the reduction on the rests and the far rule at its boundary hold. Two statements of (v) were
wrong, and both fed (iv)'s count:
- **The far side ends in a swap.** E*₃ is monotone in w, so the far runner-up grows through w_e and outgrows the entrant, which then
  races and is captured at small enough ε.
- **Formed runner-ups' side and s_d.** For these the side and s_d were attributed to the inner layer, but at rtol 10⁻¹³ they go toward
  on the target and its mirror alike, pinned at about 10⁻¹⁴σ. The referee read this as a rounding floor.

Minor issues: the proof's linear pull, far tails described as collapsed, s_d quantised to its output grid, hand-computed numbers, a
literature neighbour, 16 of 29 mutants surviving, and dead outputs.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3be-1 | major | The far side ends in a swap | R continues the mirror to s = 3·10⁹ (x = 7, 10): the swap at 8.63·10⁸, 8.73·10⁸, the entrant's departure at 1.38·10⁹, 1.39·10⁹, the entrant fitting feature 2; C4M adds 10⁻¹⁷ ([1, 1, 1], the entrant captured). (v): ∂_wE*₃ ≠ 0 on (0, A), the swap and the entrant's race; the far +1 iff s_θ < s₂ < s_{d,e}. (iv): the exception per runner-up and the windows. The title and PARADIGM D‴.51(ii)(c) (two places): "grows until it outgrows the entrant, which then migrates in its turn" |
| D3be-2 | major | Formed runner-ups' side and s_d: "a rounding floor" | The offset is not rounding: the pair's relative mode obeys exactly dD/ds = (2Gm ww_e/W)⟨r, ∇²φ(c̄)⟩D + (2(w − w_e)/W)dc̄/ds (the lighter kernel trails a moving block), and on samples the block's rest moves with W near the window's edge, ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩ = 1.14·10⁻¹²σ along y, 5σ below the top edge. Part RP (`PairSums`, `rp_job`; the sums in float64 or the 80-bit long double). The pin law (no fitted constant) gives 0.902, 0.988, 1.005 of the measured pins at x = 3, 4.9, 5.48. The pin is the same in the long double and on the mirror. With the edge nearer (4.67σ) it grows by ∂c̄*/∂W's factor and s_d falls; with the edge 8.33σ away it vanishes, the pairs collapse to identical centres, and the side becomes the arithmetic's: x = 4.9 goes to the far side and stays. (v): a tail's side is the layer's; a formed runner-up's is the block drift's, so its side and s_d belong to the sampled problem. R records the closest approach's share along e: +0.87 … +0.89 formed (the pin, −y), +0.58 … +0.48 tails (the layer, +0.48). The scan's looser tolerances do not resolve the pin. The relative-coordinate formulation was tried: in float64 the block force's rounding enters D's equation at the same level, so the precision test uses the long double instead |
| D3be-3 | minor | The linear pull | proof: δ from the block's centre c_b, c_b and the coefficients slaved by the projection alone; "the linear pull moves the block" |
| D3be-4 | minor | The reduction's justification; far tails | proof: the reduction holds while the rest's rate exceeds the coefficient's drift, failing at the fold; (v): "collapsed … held" for x ≲ 8.6, far tails (0.056–0.36σ) without a rest, waiting by their mobility ∝ w |
| D3be-5 | minor | C4M only far from its boundary; labels | C4M adds the four boundary runs (x = 7 at 2·10⁻¹², 4.5·10⁻¹²; x = 10 at 2.5·10⁻¹³, 6·10⁻¹³) and 10⁻¹⁷; the held unformed tails labelled by their side (`c4m_outcome`) |
| D3be-6 | minor | s_d quantised | `race_job`: the departure, the swap and the entrant's departure bisected on the dense output; the grid point kept (`s_d_grid`) |
| D3be-7 | minor | Numbers | the md prints the two sides' mean (1.00242 at 7σ) and R at x = 30 (s_d 3.21·10⁶, the closest approach 0.36σ); C4's tolerance and T3's range corrected; the proof cites no referee's trace; the floating-point-resolution sentence dropped |
| D3be-8 | minor | 16 of 29 mutants survive | tests: the golden report (the committed md from the committed raw), the mirrored race's s_θ ≈ 5.7·10⁶ and none for a formed runner-up, the bisected departure within its grid step, `c3_outcome` at 1.5σ and 2.5σ, `c4m_outcome`, `PairSums` against the critical right-hand side, the pin law on two windows |
| D3be-9 | minor | Literature | the over-specified mixtures' moments (Chen 1995; Ho–Nguyen 2016; Dwivedi et al. 2020) as the hold's static neighbour; Pesme–Flammarion's "successively activated", Kunin et al.'s "dormant"; "New" drops "the side kept by the inner layer" for formed runner-ups and "adds a kernel on the far side" |
| D3be-10 | cosmetic | Dead outputs | c4_job's e_t2, w_t2 dropped; C3's stale w_unbind stripped from the raw; R's s_θ only before the departure |

**Lesson.** The referee saw an offset that did not move with the tolerance or with a mirror image, and called it a rounding floor. The
law of the relative mode explains it instead. A lighter kernel trails a moving block, so whatever moves the block's rest (here the
sampling window's edge, through ⟨φ, ∇φ⟩ ≠ 0) pins a formed runner-up's offset deterministically. The test that tells physics from
arithmetic is to change one and keep the other: the long double kept the pin, and moving the window's edge removed it. An unexplained
floor is a hypothesis about the arithmetic only until a law for it is written down and tried.


## Round 34 (adversarial review of Proposition D‴.44, the replicator's own masses on a kernel dictionary under measurement noise, 2026-10-28)

**Design.** One referee, not modifying the repository. D3ag: Proposition D‴.44 and its code at 99263ba (`review/REVIEW_round34_D3ag.md`).
The referee:
- re-derived (i)–(iv) by hand and regenerated the md from the raw file;
- compared the closed-form drift with Katzenberger's at Part B's own step and correlation, for all three laws;
- projected the paradigm's AR(1) paths onto M by the noiseless flow on an early and a late window (σ = 0.03 and 0.015, h = 1 and 0.5);
- measured the Itô term and the centring's cross term along the paradigm's white flows, and took loop integrals of the drifts on M;
- checked the excess bound with the conserved K at 246 points, R_null − R_free under white and slow Ψ, D*'s convexity and its minimiser's
  support;
- re-ran three off-basis-pursuit ends with solve_ivp's status recorded and Part D's timed-out flow;
- ran 34 single-replacement mutants and web searches on the literature.

**Verdicts (D3ag).** The identities hold — the closed form equals Katzenberger's drift at Part B's own h and τ_c, the free masses'
cancellation is exact, the bound is true —, but the proposition proves its minimiser statements for the flat (√β) update and runs its
numbers on the paradigm's (log β) update, which is not a gradient on M; its medians hid the instances where the two differ:
- the paradigm's drift is not a gradient (loop integrals ∝ ε², 2.1·10⁻⁸ at ε = 10⁻²); its gauge settles at Gm|w|(1 − t⁴) = γa³t³, not
  at D*'s split; its free ends leave basis pursuit on most instances where basis pursuit fails;
- the Itô term is not O(γ/G) at the ends (0.95–1.28 of the flat drift), and the conserved law's cross term dominates there;
- (iv) is false for coloured noise and for the paradigm's update;
- "every atom used near saturation" holds only for the flat update: the paradigm's conserved ends are (N + 1)-sparse;
- Part B's AR(1) lag is the order-v expansion's limit at σ²τ_c/h ≈ 1, not a bug; its "dev" mixes in an O(σ²) offset of E[w(x)];
- the tests killed 11 of 34 mutants.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ag-1 | major | Minimiser statements describe the flat update; the paradigm's is non-gradient | (i)(c): the paradigm's drift −P_gP(∇R − ν), dν = Σ(4γω_i/(Nm))ϱ_ia_ia′_i dα_i∧dϱ_i ≠ 0 on M — "minimisers" for the flat update only; (i)(d): for free masses the paradigm's selection is each atom's Euler lag −(2ω_i/(Nm))dv_i/dζ_i, and Q (D‴.45) drifts by (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)), κ = 2 (paradigm) or 1 (flat), for every colour, without quadratic variation; (ii)(b): its ends' conditions — every alive atom at the gauge equilibrium Gβ(1 − t⁴) = γa² (D*'s split with G halved) and ω_iζ(|w_i|)sign w_i ∈ range Φᵀ, ζ = (2γa_max/(Nm))t(1 − t²)/(1 + t²): a logarithmic penalty above the threshold, so sparse, path-dependent ends. A: loop integrals (flat |∮| ≤ 1.03·10⁻¹⁹, paradigm 2.40·10⁻⁹ and 6.25·10⁻⁹ at ε = 10⁻², ratios 4.00), Q's drift against Katzenberger's to 8.6·10⁻⁷–1.6·10⁻⁶, its quadratic variation ≤ 8.5·10⁻¹⁷; C: the gauge ratio at the free ends (paradigm 1.0000 [1.0000, 1.0000] for every colour, flat 0.5000 [0.5000, 0.5000]); `mass_dictionary.q_drift`, `euler_lag`, `gauge_t`, `gauge_ratio`, `mirror_rate_eq` |
| D3ag-2 | major | "O(γ/G) relative" false at the ends | (i)(b): γa_i²/(2Gβ_ia′_i² + γa_i²) of the mass drift, O(1) at the ends; A: along Part B's white flows the Itô term 0.037 → 1.000 of the flat drift (at most 1.304; free 0.028 → 1.000) |
| D3ag-3 | major | The cross term dominates at the ends | (i)(b): the cross fitness −(γ²/(2N))a_iφ_iᵀf(K)y; (iii)(b): it favours atoms opposing the target, and its push on ‖w‖₁ balances the rest's at the end; A: cross/net 0.024 → 8.2·10¹² and the pushes (cross, rest) 5.8·10⁻⁶, 4.8·10⁻⁴ at s = 0, 3.72·10⁻⁷, −3.72·10⁻⁷ at 10⁷; `mass_dictionary.cross_term` |
| D3ag-4 | major | (iv) false for coloured noise and the paradigm's update | Title, (iv), Status: white noise and the flat update; R_null − R_free = (1/τ)log(1 − τ(γ/N)yᵀ(I + τK_free)⁻¹y) for slow Ψ (A: white −1.082422·10⁻³ = −(γ/N)|y|² at six points, slow −2.11·10⁻⁵ … −9.60·10⁻⁶); C, D: the null kernel against free masses per instance (ends differ on 4–11 of 16 on the Gaussian designs, 0–2 of 8 on the kernel dictionaries; the null kernel misses the truth where both find it on 0–7 instances, 0 on the dictionaries); test of the white constant and the slow lemma |
| D3ag-5 | major | Per-instance agreement with basis pursuit hidden by medians | C, D: every row over one instance set (n), the truth and "at basis pursuit (10⁻³)" split by bp_ok, the largest distance; numerics restated (free masses within 10⁻³ of basis pursuit where it is the truth on 12–13 of 13, 6/6, 7/7, 3/3; elsewhere on 1–2 of 3, 2–3 of 10, 0–1 of 1, 1–2 of 5) |
| D3ag-6 | major | "Every kernel used near saturation" only for the flat update | (iii)(a)–(b): the flat ends use every atom; the paradigm's conserved ends are near-vertices of {Φw = y, ‖w‖₁ ≤ a_max} with at most N + 1 atoms used (C, D: used 5.9–6.0 of 12, at most 6; 7.5–7.8 of 16, at most 9) |
| D3ag-7 | major | Part B's AR(1) lag: the order-v expansion's limit | B rebuilt: from the flow's state at s = 1500, every path projected onto M by the noiseless flow (`settle`) at the window's ends, the projected advance over the flow's ± s.e., the raw mean's offset (z); AR(1) at σ²τ_c = 0.9, 0.225, 0.05 (h = 1) and 0.45, 0.1 (h = 0.5): the projected advance over the flow's 0.008, 0.113, 0.808, 0.818, 0.970 (conserved, in order of σ²τ_c/h = 0.9, 0.9, 0.225, 0.2, 0.05), free 0.857 → 0.991, null kernel −0.340 → 0.904, fresh 0.999–1.001; (iii), Status: slow-noise claims are the reduced flows' (σ²τ_c/h → 0) |
| D3ag-8 | major | Tests kill 11 of 34 mutants | Tests: D* at atol = 0 and its minimiser equation (G = 1, 2); the bound's equality case with unequal columns and a_max = 2; one noiseless step against the law; one noisy step's covariance, mean and the AR(1) stationarity and lag correlation; Q's drift (closed form, Katzenberger's, no quadratic variation); the Euler lag; the gauge equilibria and the mirror rate against η's closed form; loop integrals; the cross term; R_null − R_free; `summary` on hand-built states; `filt_of`; the report on a synthetic result; `settle`; the floor; FlowFailed and FlowTimeout's last state — 30 of 30 single-replacement mutants killed (library 20, experiment 10) |
| D3ag-9 | minor | "Weighted basis pursuit above the threshold" | (ii)(a): a dense ℓ^{4/3}-perturbation of weighted basis pursuit (D* convex, D*′(0) = 0); C: the flat end against the D* minimiser (median, largest), that minimiser against basis pursuit, the flat ends on the truth (the flat ends 5.9·10⁻⁵–5.7·10⁻⁴ from the D* minimiser (median), that minimiser 4.0·10⁻⁴–2.0·10⁻³ from basis pursuit; the flat ends on the truth 2/16, 0/16, 2/8, 0/7) |
| D3ag-10 | minor | The slow count's leading term | (ii)(c): min(k, N) log τ_c/τ_c, K's eigenvalues above 1/τ_c, the O(1/τ_c) term not uniform, gradient selections only |
| D3ag-11 | minor | The conserved-K inequality's proof; Part A and the test checked tr K_free; boundary fits are limits | (iii): both Cauchy–Schwarz steps, the boundary value tr K → (γa_max/N)Σ|φ_i|²|w_i| − (γ/N)|y|², the boundary a limit; A's table: tr K and tr K_free; the test on the conserved K with the equality case |
| D3ag-12 | minor | The E bound vacuous at a_max ≥ 3.16, loose; no saturation at a_max = 4 | (iii): the bound's range and the threshold γa_max³/(Gm); E: the bound, max|a|/a_max and the threshold against the median used |w_i| (E/m 0.0021–0.3469 against the bound 0.158–1.265; at a_max = 4 the median used coefficient 0.038 below the threshold 0.5333, max|a|/a_max 0.862) |
| D3ag-13 | minor | Part D's timed-out row mixed denominators | `flow`'s FlowTimeout and FlowFailed carry the last state; C, D: one instance set per row with n, a stopped flow listed with its last state (one: k = 3, instance 1, free, flat, over the limit at s = 2.83·10⁷, on the truth) |
| D3ag-14 | minor | Part B fresh: "within 0.002" hides an O(v) offset | B: the projected comparison (fresh rows), the raw offset in s.e. |
| D3ag-15 | minor | `flow()` ignores sol.status; the null kernel's E/m included \|w_0\| | `flow`: a failed piece redone with Radau, else FlowFailed; `summary`: E/m = 1 − ‖w‖₁/a_max over the real atoms |
| D3ag-16 | minor | Novelty overstated | Status: the cancellation a corollary of the log-coordinate Itô correction (exponential learning, stochastic mirror descent), D*'s branches the quadratic and depth-3 penalties, support counting (Tipping; Wipf–Rao; stochastic collapse); new: the Euler lag, Q's drift, the gauge equilibrium and the logarithmic selection, the inflation and its cross term, the null kernel's white equivalence |
| D3ag-17 | cosmetic | Renormalisation argument; M for free masses; representability; the null kernel's comparison; "nan" | (i)'s proof: the radial drifts P-normal, removed by the projection; M the law's own; "basis pursuit's solution representable"; the null kernel compared with free masses per instance; "—" for undefined E/m |

**Lesson.** A potential describes one update; before a minimiser theorem is put beside numbers, a one-second loop integral tells whether the
dynamics measured is a gradient flow at all. Here it was not, and the right description of the paradigm's update was not a modified
potential but a conservation law: Q, conserved by the noiseless free law for every target, has no noise of its own, so the selection's
drift on it is a one-line Itô computation — and its zero is exactly the gauge equilibrium the referee measured. Medians over instances where
every method agrees hide what distinguishes the methods; the rows have to be split by the instances where the candidate answers differ. A
reduced flow is validated in the limit it claims — σ²τ_c/h → 0 here — and against the projected state, not the raw mean.

## Round 56 (adversarial review of Proposition D‴.60's first draft, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bc: Proposition D‴.60's first draft with its code and results at 17d6217
(`review/REVIEW_round56_D3bc.md`). The referee:
- regenerated the md (byte for byte) and reran an A row and a full B3 cell (bit-identical);
- re-derived (i), `second_order` term by term, the per-kernel √K_max,i = (Gm)^{1/2}|w_i|/(σe^{1/2}) and the design bound (m ≲ 4.83);
- ran controls on the loop: nine half-integer tones at m = 3 and 12, k_e = 5–40 at 13.1 and 41.3, the commensurate 21.4 against 21.3
  over 60 periods, 300 periods at 21.37 and with the 16 tones, the 16 tones one at a time and in pairs;
- split the exact second order into its response to the forcing's constant part and to its 2ω part;
- traced an ejection, reran C's thresholds under three more seeds, with single tones, a band above √K_max and the WFR metric;
- ran 42 mutants and searched the literature (snippets).

**Verdicts (D3bc).** (i) is correct and classical, `second_order` is right term by term, the even part converges to s²⟨θ₂⟩ with an
error ∝ s², and the drift's resonance is real (⟨w₂⟩ at m = 12 over m = 3: 20 at 13.1, 1.1 at 41.3). But (iii)'s causal chain is
unsupported: a kernel leaves its feature by a sudden ejection that single resonant tones at twice the rms never cause, and a band moved
above √K_max still ejects at m = 12 while m = 3 never does; C's thresholds were single-seed and not monotone; the averaging's failure is
the commensurability of the tones' combinations with the reference's harmonics, not slow beats (and Krylov–Bogoliubov does not apply);
the odd part alternates or grows secularly instead of being a start offset; ⟨w₂⟩ saturates after about 130 periods; the four sources
are D‴.31(iv)(c)'s and D‴.33(i)'s; the WFR threshold claim had no data; several numbers were hand-computed; 19 of 42 mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bc-1 | major | (iii)'s mechanism and design rule unsupported | (iii) numerical: ejection a first-passage event of multi-tone noise, its level falling with Gmw and with resonance and rising with k_e; the rule (√K_max below the band) removes the resonance, not the ejection. C2: one tone 13.1, 15.7, 21.3 at s = 0.02, 0.04 never ejects (0.069); the band 25–75: m = 3 never to s = 0.08, m = 12 from 0.06; m = 4.8, 6, 8 from 0.08, 0.06, 0.03; tone pairs (18.9/22.5, 18.9/22.4, four seeds' phases) never at 0.02 |
| D3bc-2 | major | Thresholds single-seed, horizon-bound, non-monotone | C over the seeds 7, 1, 2, 3: per kernel the least level a seed left at (how many), the level from which every seed left at every level above, '*' where a seed leaves at one level and stays at a higher, the median period of leaving (118–256); "within 300 periods" |
| D3bc-3 | major | The averaging-failure condition; KB | (ii) and its proof: the averaged loop is the response to the quadratic forcing's constant part; the 2ω_k and ω_j ± ω_k parts give period-mean increments rotating by ω_kT — bounded unless commensurate with ω_ref, of the averaged part's order where K(t) crosses ω_k²; "slow beats" and KB dropped; B3: WFR single tones 0.9996–1.0010 against the 16 tones' 0.874–1.022 |
| D3bc-4 | major | The odd part is not a start offset | (ii): rotating Floquet increments, bounded by about 1/\|sin(ω_kT/2)\|, secular at commensurability; B4 prints the first kernel's centre odd part per period (+0.338, −0.312, … at 13.1; −0.032, +0.030, … at 21.3; rotating at 21.37) and the commensurate 21.4's coefficient odd part per s growing to −0.467 |
| D3bc-5 | major | "⟨θ₂⟩ grows with time" | B4 (300 periods, s = 0.00125, 0.0025): the first kernel's ⟨w₂⟩ +221, +591, +874, +1.01·10³ at the periods 30, 100, 200, 300 (13.1), +225 to +234 from period 200 (16 tones); (ii): a second-order offset reached with the slow coefficient mode's time constant |
| D3bc-6 | major | Novelty overstated | Status: the sources D‴.31(iv)(c), D‴.33(i), the k_e sign D‴.33(iv); Widrow et al., Glover, Gibson–Annaswamy–Lavretsky, Riedle–Kokotović, Sanders–Verhulst–Murdock, Neishtadt; New: the exact non-averaged second order with its commensurability condition, the mobility Gmw as the stiffness with √K_i and the bound on m, the ejection data |
| D3bc-7 | minor | k_e as damping only | (iii): the oscillator's damping and D‴.33(iv)'s rectification, whose sign k_e reverses |
| D3bc-8 | minor | No WFR threshold data | C2: WFR at k_e = 5 and 10, s = 0.02, 0.04, 0.08: never (largest distance 0.023) |
| D3bc-9 | minor | (ii) cannot reach (iii)'s levels | Status: the even part 0.8122 of s²⟨w₂⟩ at s = 0.005 near the resonance; (iii)'s levels start at 0.0075 |
| D3bc-10 | minor | Per-kernel resonances; the jitter's ratio | C's caption per kernel, each law's own fitted state (a bug fixed: m = 3 had used m = 12's; now 7.88, 5.25, 6.57 against 15.76, 10.51, 13.13); (i): √K_i = (Gm)^{1/2}\|w_i\|/(σe^{1/2}); A: each kernel's jitter peaks near its own √K_i, far above them the paradigm's jitter is mostly the centres' mobility |
| D3bc-11 | minor | 19 of 42 mutants survive | tests: `b3_job` on two periods (NP_B, MARKS_B, S_B3 monkeypatched), `b2_jobs`' labels, the envelopes' ē coupling against central differences of the law, the blur at x_m + ē and the correlation at x_m against direct tone averages, C's caption at each law's fitted state — 10 of 10 targeted mutants killed; 24 tests pass |
| D3bc-12 | minor | Hand-computed numbers | A prints the ranges per metric, B3 the centres' ratios at every s (0.9796–1.0013, 0.9201–1.0048, 0.6994–1.0121); 1.0%, 0.58%, 0.12%, 248, 8% and 3% replaced by printed values |
| D3bc-13 | minor | B's "averaged / loop" includes the odd part | the column dropped; B compares the averaged loop with the loop's even part only |
| D3bc-14 | cosmetic | Empty terms; the bound; the 16 tones' odd part | "three sources carry the drift" (B2); C2 prints the bound m ≲ 4.83; the odd part's departure for the 16 tones (2.78·10⁻⁴ against s⟨w₁⟩ = 4.80·10⁻⁴) |

**Lesson.** A mechanism confirmed on the quantity it predicts — here the drift's size, resonant as claimed — is not thereby the cause of
the event it is used to explain. The controls that separate the two cost one run each: the resonant tone alone at twice the level, the
band moved above the resonance, intermediate m. They turned a design rule into a condition that removes the resonance but not the
ejection. And a period mean over a periodic orbit is a Floquet sum: each tone's contribution rotates by ω_kT, so commensurability, not the
slowness of a beat, decides whether it averages out.

## Round 50 (adversarial review of Proposition D‴.57, the transit's dead zone, 2026-10-28)

**Design.** One referee, not modifying the repository. D3aw: Proposition D‴.57 with its code and results at 502c522
(`review/REVIEW_round50_D3aw.md`). The referee:
- regenerated the md (byte for byte) and checked every χ² number with mpmath (P(late) = 0.046085, m_c = 53.82663, the ratio 2.08100);
- reran all 1200 draws of Part A with the entry order and the vectors u, f̂_K, f̂_c recorded, decomposed u on the fit's tangent basis and
  linearised the fast flow at each noiseless block;
- reran B's two worst draws at rtol 10⁻⁹ and 10⁻¹¹, comparing the fits themselves;
- reran C with a dense record (the filtered energy, the entry times, the trajectories), at another seed, other thresholds, with the
  control's f̂ at the true state and the adaptation driven by the true error, and at s = 0.01 to τ = 800;
- ran three mutants and searched the literature (memory where pages would not load).

**Verdicts (D3aw).** (i) holds in the sliding limit on fixed data, close to tautologically, and Part A does not simulate gate (b);
(ii)'s identity is exact at first order and its law right for a fixed entry order, but the cross term's mean is not zero when the order
swaps (correlation −0.32 to −0.41 on those draws), the slowest-mode step is not the mechanism (u is the last entrant's amplitude
deficit, δw ≈ 1.76s), and "about twice the oracle" is a coincidence of N/K ≈ 20 and δ = 0.1 — the ratio grows like (Nδ/K)^{1/2}; (iii)
holds for sharp gates on fixed data; (iv)'s failure is real but misdescribed — the gates chatter, the twelve "formed kernels" are an
eject-and-replace cycle on two features, nothing supports "fits the noise at the edge", and the negative result was overgeneralised
from one level, one regressor and one gate setting; C's numbers were single snapshots; three mutants survived.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3aw-1 | major | "u's direction is one per target, falling with s" | Part A compares u only within one last entrant (the reference left out): \|cos\| ≥ 0.9996, 0.9954, 0.9388 (medians 1.0000, 0.9998, 0.9973); the falling median was order swaps |
| D3aw-2 | major | The cross term assumes the order fixed | (ii) and Status: given the noiseless entry order; Part A splits kept and swapped orders (`classify`): kept 1.016, 0.978, 0.957 of the prediction, swapped 1.53, 1.05, 1.24 with correlations −0.32, −0.03, −0.17; the o(s²) not uniform in the order's margins |
| D3aw-3 | major | "About twice the oracle" | title and (ii): the ratio (2 + N(2δ + δ²)/(2K))^{1/2}; A3 on N = 61, 181, 541, 1621: 2.08, 2.85, 4.76, 7.71 against 2.08, 2.91, 4.58, 7.66; the energy dead zone inconsistent in N, the shortfall δw = 1.80–1.43s not shrinking |
| D3aw-4 | major | C's "gate closed at τ" a coarse sample | dense C: the gate's shut and fully open shares, first closing and last reopening (x_m, s = 0.02: formation gate shut 0.11, fully open 0.72; energy gate 0.54, 0.30) |
| D3aw-5 | major | The mechanism misdescribed | (iv): an eject-and-replace cycle on the features at −0.52 and 0.61 (C's timelines); the variants: the true error in the adaptation 8 formed, 6 ejected; the centres boxed: every kernel inside at the end; Part D from the fitted state: kernels leave only at λ ≥ 1 within 300 periods, none under WFR to s = 0.16; the ejection is D‴.60's first passage |
| D3aw-6 | major | The negative result overgeneralised | title and (iv) restricted: 16 tones of rms 0.02 at x_m, 0.01 and 0.02 at x̃, the paradigm's metric, the gates at 1.2s/1.5s; at x_m, s = 0.01 every adaptation keeps 3 kernels; the energy gate at 2.0s/2.5s forms 1 kernel; C2: one tone, two tones and the WFR metric eject none |
| D3aw-7 | minor | C's snapshots | the fit's median [5–95%] over the last half (x_m, s = 0.01: 0.0254 ungated, 0.0613 energy gate) |
| D3aw-8 | minor | The slowest-mode step | the proof: u = δw·P⊥φ_last from the relaxed block and the overlap; u depends on ξ through the crossing |
| D3aw-9 | minor | The gates as implemented | "As implemented": A stops the flow at the crossing and does not simulate (b), B a one-time mask, C smooth gates with a partial drive for \|w\| ∈ (θ/2, θ) |
| D3aw-10 | minor | Three mutants survive | tests: the formation gate drives a formed kernel while the others rest (M1); the stop law against quadrature (M3); `gate_stats` and `kernel_events` on hand-made series and the report on synthetic raws (M8) — 16 of 16 single-replacement mutants killed |
| D3aw-11 | minor | Literature | Egardt, Peterson–Narendra, Samson, Ioannou–Sun §8.5; Blanchard–Hoffmann–Reiß, Celisse–Wahl, Blanchard–Mathé; Burger–Gilboa–Osher–Xu, Burger–Möller–Benning–Osher; New narrowed to the frozen discrepancy principle's O(s) amplitude shortfall with its first-order law given the stop |
| D3aw-12 | minor | Part A's exactness checks the entry set | `classify` and `fit_after` (one entry per iteration no longer assumed) |
| D3aw-13 | minor | B reuses the noise across levels | `b_noise`: fresh noise per level |
| D3aw-14 | minor | The x̃ rows confounded | (iv): at x̃ the noise also enters the regressor and the noise-free loop already forms 5 kernels |
| D3aw-15 | cosmetic | Hand-computed medians against an rms | A prints the rms of \|u\|/(s√N) (0.567, 0.568, 0.566 order kept) against the conditional prediction 0.572 |
| D3aw-16 | cosmetic | "To 10⁻⁴"; "to 2.2·10⁻⁹" | the identity to 7.5·10⁻⁴ at worst (5–95% [0.9999, 1.0002] at s = 0.1); B compares the fits: 2–4·10⁻¹⁰ (formation), 3–10·10⁻¹⁰ (energy) |
| D3aw-17 | cosmetic | Slips | the proof of (i); the monotone decrease from the flow; the noise floor at the desired regressor (0.881s) |

**Lesson.** A dead zone tuned and tested on fixed data says little about a loop whose data move. Here the dead zone's sliding limit had
an exact law, while the closed loop's failure — kernels pushed off their features by the noise — happened behind the gate, through
the centres' mobility. A gate on the error's energy cannot see that. The negative result needed its controls: a lower level, other
thresholds, the true error in the adaptation, bounded centres, another metric. Run before writing the claim, they show its scope:
one level, one regressor, one metric. And when a law holds "given the order", the draws where the order changes must be reported
apart; a pooled mean hides them.

## Round 59 (adversarial review of Proposition D‴.52 as revised after round 55, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bf: Proposition D‴.52 as revised after round 55, with its code and results
at 03d01bb (`review/REVIEW_round59_D3bf.md`). The referee:
- regenerated the md byte for byte and reran A, B, C/D and E cells bit-identically, E consuming the main rng as `e_job` does;
- checked (i) on two pairs it was not tested on, and the null law at the true fit in all 15 E cells;
- computed the score test at E's first fits, the runner-up's filter and its cosine with f₁ − f₀, the effective critical values;
- ran the opposite-sign pair (0.8, −0.4) and the near-tie offset c̄ = 0.049 over separations (`detect.py`);
- identified the transit's knots with the LARS knots on a fixed dictionary (τ₂λ₂ = 1 to 5·10⁻¹⁵ over 200 noisy pairs);
- ran 36 single-point mutants and searched the literature.

**Verdicts (D3bf).** The deterministic core, the new null law and the placement at the Cramér–Rao bound hold. But the knot's
"near-Neyman–Pearson" power held at one offset only: across E's 15 cells the bound exceeds the knot's power by up to 0.31, and the
matched-filter explanation predicts the opposite pattern (the runner-up's filter is best aligned where the power is worst). The losses
come from what the knot adds to the filter — the minimum over the dictionary, and near a tie the margin's own noise —, and a score
test computed from the same first fit sits at the bound in every cell. For an opposite-sign pair the knot test is biased (power below
α). The second knot is the LARS second knot, so Δτ₂ is a spacing statistic with known conditional nulls, and the calibration is the
parametric bootstrap at the fitted null. Minors: the resolution limit's offset dependence, the bootstrap conservative near a tie, D's
jump variable, biased split cells, C's 200-draw thresholds, misquotes, E's data path untested (17 of 36 mutants surviving).

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bf-1 | major | "Near-Neyman–Pearson power" fails across the cell; the matched filter does not explain it | (iii): the knot falls below the bound by what it adds to the runner-up's filter — the dictionary's minimum (far kernels share the H0 tail), the margin's noise near a tie (it swaps the first entrant); (iv) new: the score z = ⟨r̂₁, ∂²_cφ(ĉ)⟩/(s‖R̂∂²_cφ(ĉ)‖) at the first fit, N(0, 1) under H0 and with the bound's mean under H1 to first order; Part S scores C's and E's draws: z's power −0.057 … +0.052 from the bound and +0.022 … +0.332 above the calibrated knot in all 15 E cells, its pooled size 0.0493 [0.0459, 0.0529] (t-version 0.0506) |
| D3bf-2 | major | The knot test biased for opposite signs; only (0.8, 0.6) tested | Part F: (0.8, 0.6) and (0.8, −0.4) at 8 offsets × 3 levels × 6 separations; the knot's Δ₅₀ 1.846–3.344 times the bound's for (0.8, −0.4), its least power 0.000–0.025 in every cell, not monotone in Δ at c̄ = 0.032; both walls' knots at α/2 each (`knot_walls`, `walls_at_fit`) 1.707–1.840; z one-sided 0.965–1.058 for both signs; (iii) states the bias and the walls' remedy |
| D3bf-3 | major | Novelty: the LARS knot, the spacing test, the parametric bootstrap, the resolution literature | (iii): τ_k = 1/λ_k on a fixed dictionary (u_j(τ) = τc_j(1/τ)), Δτ₂ a spacing of the sliding LARS path; Status: the spacing test (Tibshirani–Taylor–Lockhart–Tibshirani 2016), its power (Azaïs–De Castro–Mourareau 2018, 2017), the covariance test, McLachlan 1987, Harris, Helstrom, Shahram–Milanfar, Liu–Nehorai (Smith's criterion), Kulaitis–Munk–Werner, Neyman–Scott's C(α), Osher–Ruan–Xiong–Yao–Yin, Amari–Park–Ozeki, Wei–Amari, Dwivedi et al., Heinrich–Kahn; New narrowed to the delay law and its branches, the knot's losses measured, the score's two readings in the paradigm (the moving width's drive, the relative mode's rate) and the mobility-weighted rates |
| D3bf-4 | minor | The limit depends on the offset; the no-tie condition | F gives Δ₅₀ per offset over both halves of the cell for six tests; (iii): near a tie, where the runner-up's margin is comparable with its noise (sd ≈ 0.16s), the quantile falls faster than 1/s |
| D3bf-5 | minor | The bootstrap conservative near a tie; s known | (iii): "valid, conservatively near a tie; s known"; E prints Clopper–Pearson intervals (0.019 [0.011, 0.030], 0.022 [0.014, 0.033]; pooled 0.0457 [0.0424, 0.0492]) |
| D3bf-6 | minor | D's jump variable | (v) and D: the positions jump at s ≈ 12Δ², where the pair stops being detected ((iii)'s order); D prints the position error against s/Δ² (≥ 0.2 in the 9 cells with s/Δ² ≥ 12, ≤ 0.152 in the 31 with s/Δ² ≤ 9.8) |
| D3bf-7 | minor | The split's error below the bound in biased cells | (v): at the bound up to s ≈ 24Δ³, below it beyond, where the bounded split is no longer estimated without bias; D prints the split's error over its bound against s/Δ³ |
| D3bf-8 | minor | C's thresholds from 200 H0 draws | C′: C's draws against the knot's law at the truth (20000 samples, `knot_walls`) and the score, with the threshold's false alarms and Δ₅₀'s bootstrap intervals |
| D3bf-9 | minor | Misquotes and a floating-point count | B prints the flow's separation rate over its second-smallest (17.8, 21.4, 11.9 at Δ = 0.0125); E's within-0.05 count with a tolerance (13 of 15); the fixed threshold's out-of-sample range c̄ < 0.032 (0.027–0.044); E's std-ratio column names j₁ and j₁ ∓ 1; the near-dipole's own split error 500.7 |
| D3bf-10 | minor | E's data path untested; tests on the code's own formulas | tests: `e_job` on a few draws (p against a direct `knot_law`, H0 from the one-bump fit, distinct seeds), `crb` on an orthogonal design and `np_power` against Φ(‖f₁ − f₀‖/s − z₀.₉₅), `d50_bootstrap` on a deterministic table, `knot_of` over several noise vectors (error ∝ s), the score statistic, the knot per wall, `f_job`/`s_job` on a few draws, F's report on synthetic jobs (24 tests) |
| D3bf-11 | cosmetic | README; E's column; the quantile ratio | the drafts README updated; the realized size at the law's quantile printed with its interval beside the ratio |

**Lesson.** A detection statistic is judged against the efficient test on the same data, not against its own calibration. The knot
looked near-optimal at the offset where it was tuned; the score test, computed from the first fit the path already makes, sits at
the Neyman–Pearson bound everywhere, and its two readings — the moving width's drive and the split mode's rate — are statistics the
paradigm computes anyway. A test with a sign-dependent drift must be run for both signs before "counts two features" is written: for
an opposite-sign pair the drift empties exactly the channel through which the knot comes early.

## Round 60 (adversarial review of Proposition D‴.58 as revised after round 54, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bg: Proposition D‴.58 as revised after round 54, with its code and results at
fef9c70 (`review/REVIEW_round60_D3bg.md`). The referee:
- regenerated the md byte for byte, ran the 45 tests, and reran three B paths, seven B2 jobs and four A rows bit-identically;
- re-derived (i), (ii) and (v), including the three proof gaps round 54 asked to close, and spot-checked untested ratios;
- reran the first non-optimal blocks of 9 paths with TRF as polished, TRF without its gradient test, LM and VP-LM, recording the
  scaled gradient, its rounding and the Jacobians' condition numbers, and reran B's repeated variant with the gradient test disabled;
- categorised B2's later misses by their reference optimum against the experiment's floor E_FLOOR;
- split the entry statistics by regime and checked C(c) ≥ 0 after all 258 optimal blocks;
- ran 36 single-point mutants and searched the literature.

**Verdicts (D3bg).** (i), (ii), (iv) and (v) hold; the numbers match. But the new explanation of where the optimal ranges end was
refuted: the repeated polish's stops were `polish_block`'s absolute gradient tolerance 10⁻¹⁵ — met along a slow mode at a fit error of
10⁻⁹–10⁻¹⁰ — not a loss of descent of the joint least squares; without it the same solver goes on (ranges 5 and 10 instead of 4 and 9
at τ/σ = 1/4, 3/4), and what separates variable projection is the conditioning. The later blocks' "basin" count used optima below the
experiment's own floor. 21 of 36 mutants of the new code survived, including every one in `sliding_iss`'s repolish, budget and entry
log. Minors: an unfixed range (95.7%), "not the stop", pooled entry statistics, a vacuous stationarity test, wording, literature.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bg-1 | major | The ranges' end: an absolute gradient tolerance, not the joint least squares | `polish_block(gtol=…)` (None: no gradient test) and `sliding_iss(polish_gtol=…)`, defaults unchanged; B's sixth variant, the repeated polish without the test: optimal to K ≤ 5, 7, 10, 12 at τ/σ ≤ 1 (one kernel further at 1/4 and 3/4; at 1/2, K = 8 1.233 times above the optimum, four polishes on the budget); (iii): the ranges are the solver's — its budget, acceptance radius or gradient test —; B2 prints the condition numbers at the first non-optimal blocks: 5.4·10³–3.6·10⁵ for variable projection's centres, 1.3·10⁶–5·10⁸ for the joint (w, c) |
| D3bg-2 | major | Later blocks judged against sub-floor optima | `b2_line` counts only optima ≥ E_FLOOR: 33 of 38 reached (the 20 first non-optimal blocks among them), 1 at 1.001 times (in the basin, unconverged), 4 at 16.36–24.74 times with nodes 2.14–3.99σ from the optimum's (`b2_job`'s dc); the 78 below the floor not counted |
| D3bg-3 | major | 21 of 36 mutants survive | tests: the polish's budget and gradient test passed to every polish, a repeated polish accepted only if it improves, at most `repolish` times, the chunks resetting to the elapsed time, a non-improving repeat changing nothing, the log's rejections the path's count, `entry_info` against the duals rebuilt from the path (three paths: a wide feature, three bumps with a negative dual, a loosely relaxed block), `polish_block`'s gtol, a failed chunk counted as neither stall nor redo; `b2_line` on asymmetric data, `entry_lines` by regime, `stop_law` against its closed form with a Monte Carlo check of the moments, `cs_vacuous` at its threshold, `merge_B`, `b_job` reproducing the raw on two paths; 14 of 15 targeted mutants of `sliding_iss`/`polish_block` killed, the survivor equivalent outside exact ties |
| D3bg-4 | minor | "95.7–100%" | 93.0–100% (A3) |
| D3bg-5 | minor | Pooled entry statistics | B's table by regime; (iii) quotes the optimal blocks': the largest rate in 5 of 316, its share a median 0.21, the dual a median 0.959, 314 beyond the outermost node (a median 1.22σ), 2 between (0.989–0.992) |
| D3bg-6 | minor | The stationarity test vacuous | `cs_vacuous`: 53 of the 116 non-optimal blocks pass it at any state by Cauchy–Schwarz (errors up to 4.0·10⁻¹⁰); (iii) says so |
| D3bg-7 | minor | Numerics wording | the gradient test's stop (status 1) at τ/σ = 1/4; within 10⁻⁶σ and C ≥ 0 to 10⁻⁶ for the base variant; the largest gap at the stop in 0, 4, 1 draws; 12 the path's cap (censored) |
| D3bg-8 | minor | Literature | Gnewuch–Hefter–Hinrichs–Ritter 2022, Irrgeher–Kritzer–Leobacher–Pillichshammer 2015, Gnewuch–Ritter–Rüßmann; the Mercer expansion's ratio ε²/(α² + δ² + ε²) = (β² − 1)/(β² + 1) = q; the classical rule's error below √2 (Kuo–Woźniakowski's divergence concerns a reweighted fixed-weight rule); Derevianko–Plonka's weighted problem stated, article 6 of volume 24; New narrowed |
| D3bg-9 | cosmetic | K(s) twice; "to within 1 − δ_t(K)"; formats | the root K*(s), the predicted stop ⌈K*(s)⌉ (md caption and table too); "its error the rule's times 1 − δ_t(K)"; `b2_line` at four digits |

**Lesson.** A stopping test is calibrated for a scale. An absolute gradient tolerance certifies nothing along a slow mode of a fit whose
residual is 10⁻¹⁰; the check that settles it — rerun the same solver from the same state with only the suspect setting changed — took
a variant of B, and moved two of the three ranges. A miss measured against an optimum below the experiment's own floor is evidence about
the solver, not about basins; and a statistic pooled over regimes is not to be quoted for one of them.

## Round 45 (adversarial review of Proposition D‴.56's first draft, 2026-10-28)

**Design.** One referee, not modifying the repository. D3ar: Proposition D‴.56, the transit under measurement noise, with its code and
results (`review/REVIEW_round45_D3ar.md`). The referee regenerated the md byte for byte and reproduced Part B bit for bit on 1200 paths;
recomputed all 40 noiseless paths and their responses (entry gaps, margins, R_K's trace, the union bound per target, a KS test of
|R_Kξ|²/s² against χ²₅₅, the PIT of W·s under each target's limit law); checked `path_response` against central finite differences on six
targets of their own (near ties of 0.86% and 1.1%, a same-sign and an opposite-sign close pair, amplitudes 1 : 0.1); ran 19 mutants;
audited Part E's raw paths (time-outs, short paths, Radau paths, exits after the (K + 1)-th entry, zero gaps) and reproduced three
anomalous ones with and without the polish; ran the two stopping rules on amplitudes (1, 0.3, 0.1), on K = 5 and to K + 8 entries; fitted
the oracle from the truth on 400 draws; searched the literature.

**Verdicts (D3ar).** The first-order response is right (second order in h on every target), and so are the window's law for small s and
the oracle at the K-th fit. Wrong: (ii) applied along the noiseless order where noisy draws swap near-tied entries (the swaps cause most
of the "exact, criterion ≥ 1" disagreements), and a probability "+O(s²)" from a local expansion; past the first noise kernel the solver
returned non-fits flagged as relaxed — an edge kernel leaving the window with a coefficient up to 2.65·10⁹, thresholds relative to max|w|,
exit/re-entry bursts at one τ —; the largest relative gap is not a stopping rule without s (it needs the path past the stop over a horizon
chosen knowing K, and fails with the amplitude range, with K = 5 and with the horizon); the tests did not protect the published numbers
(the criterion's intermediate rows, `b_job`, `e_job` and the report untested). Minors: the limit law untested, a circular χ² check of the
late stops, novelty, Part E's selection and a wall-clock budget, (D)'s wording, D‴.53's hypotheses, index bookkeeping.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ar-1 | major | (ii) linearised along the noiseless order where it breaks | (i): the ball's radius — r₁ = min(1 − \|u⁰_j(τ_k)\|)/\|U_{k,j}\| over the kernels not yet entered, the later entrants included (a near tie is the later entrant's small margin), and the blocks' nondegeneracy; (ii): the order criterion over every entry's not-yet-entered kernels (`U_all`); B: the order kept, the swaps 0, 2, 26, 64, 81, the order criterion agreeing per draw on 400, 400, 398, 392, 372, every disagreement within the level's largest linear-dual error of the wall; B2: 20 targets of small radius (r₁ 0.0026 … 0.019) at s = 0.001, 0.003, 0.01, the order kept on 400, 389, 344 against 399.8 ± 0.4, 388.3 ± 3.1, 337.5 ± 6.7 |
| D3ar-2 | major | The probability "+O(s²)" unproved; the union bound never computed | (ii) two-sided on the ball: the order kept if crit < 1 − C\|ξ\|², not if crit > 1 + C\|ξ\|²; P(crit < 1) the first-order probability, off by at most P(\|crit − 1\| ≤ C\|ξ\|²) + P(\|ξ\| > r₀) — of order s where failures are frequent; B: the linear model's expected counts with z (+0.2, +1.0, −2.0, −0.6, +1.3) and the union bound 0.1, 7.7, 40.0, 219.6, 400 (capped) against 0, 2, 27, 103, 273 failures |
| D3ar-3 | major | Non-fits past the first noise kernel flagged as relaxed | `sliding_iss` (4fd8082): `box` keeps the centres in the window ± 2σ and stops the path with 'runaway' when a relaxation carries one out; `w_scale` thresholds on the data's scale; `max_rhs` a deterministic budget; `ev_iter` maps entries to iterations; `relax_chunk` redoes failed solves and the opt-in exit rule (`exit_kkt`) lets a kernel leave only if its correlation stays inward without its own contribution (1c5a996); (iv): past τ_K + W the fit can have no minimiser; E rerun: the runaway paths reported, no exit after the (K + 1)-th entry on any path |
| D3ar-4 | major | The largest relative gap not a stopping rule without s | (iv): a selection over a fixed horizon, not a stopping time, needing the path past the stop; it selects the K-th entry iff W/τ_K exceeds every other relative gap of the horizon — the features' own (growing with the amplitude range and K) and the later noise kernels' (not bounded uniformly); Status: the rule presumes features of the kernels' width; E's four families: on the exact paths the rule selects the K-th exactly where W/τ_K exceeds both kinds of gaps — failing from s = 0.1 by the noise kernels' gaps (base), from s = 0.03 by the amplitude range (wide: the features' gaps 2.38–3.3 against W/τ_K 11.9, 4.87, 0.502) |
| D3ar-5 | major | The tests do not protect the numbers | seven new tests (finite differences of the duals at every entry against `U_all[k]`, the criteria's rows, radius and union bound, classification and both rules, a seeded `b_draw`, Part E's report on hand-built paths, the box, `w_scale` and the budget on the runaway draw): 13 of 13 mutants killed, the referee's survivors among them |
| D3ar-6 | minor | The limit law of sW untested | C on the first order's regime (the order kept, crit < 1); the PIT of W·s under each target's law, KS p = 0.28, 0.88, 0.03, 5.5·10⁻⁵, 4.5·10⁻¹⁴: W ∝ 1/s and its limit law for s ≤ 0.01; the draws keeping the order with crit ≥ 1 (0, 0, 0, 2, 7) and the swaps counted apart, with W_pred < 0 |
| D3ar-7 | minor | The χ² check of the late stops circular | (D): the path's own late stops 15, 12, 20, 18, 8 against the binomial law 18.4 ± 4.2, …, 9.6 ± 3.0 (p = 0.48, 0.15, 0.63, 0.71, 0.74) |
| D3ar-8 | minor | Novelty overstated | Osher–Ruan–Xiong–Yao–Yin, Vaškevičius–Kanade–Rebeschini, Wainwright's primal–dual witness, Duval–Peyré, the sliding Frank–Wolfe and ADCG, the covariance test, Ahn–Horenstein, Morozov, Yao–Rosasco–Caponnetto, Burger–Resmerita–He cited; New restricted to the off-grid, moving-centre first-order response, the order criterion, its two-sided form and the limit law of sW; "unbiased at first order" |
| D3ar-9 | minor | Part E scored on exact paths against the path's own fit; a wall-clock limit | both rules scored on every finished path against the oracle polished from the truth; a budget of 3·10⁶ right-hand-side evaluations |
| D3ar-10 | minor | (D) overstates the gap rule, mislabels the error | "over the first K + 1 entries"; "the K-th entry's fit" |
| D3ar-11 | minor | D‴.53 applied to any data | the loop follows the noisy path on [0, τ_K + W), on draws whose path satisfies (H1)–(H4) |
| D3ar-12 | minor | Index mismatches | `ev_iter` maps entries to fits; paths with later exits flagged (none in the rerun) |
| D3ar-13 | cosmetic | The recursion's start and the noise model unstated | (i): τ₀ = 0, δτ₀ = 0, δu(0) = 0, R₀ = I; the noise i.i.d. Gaussian of variance s² in (ii)–(iv) |

**Beyond the report.** E's rerun under the box finished in 10 h on two workers: of the 100 base paths per level 94, 94, 81, 69 finished
within the budget, and a runaway kernel stopped 19–26 of them before the horizon — the fit without a minimiser is the common case past
the first noise kernel, not an anomaly. The largest-gap rule's success is exactly its condition (W/τ_K above every other gap) on all 11
rows, so its failures are attributed per row to the noise kernels' or the features' own gaps.

**Lesson.** A first-order statement is applied only where the order it linearises along holds: near ties make the later entrants'
margins the binding ones, and a draw that swaps two entries is outside the expansion however small its dual error looks. A path past
its last justified event needs its own existence checked before its fits are scored — here the fit had no minimiser, and a "relaxed"
flag computed on the solution's own scale certified a kernel at c = −1.86 with w = 2.65·10⁹.

## Round 62 (adversarial review of Proposition D‴.54 as revised after round 58, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bi: Proposition D‴.54 as revised after round 58 — the far side's swap and the
relative mode's pin law (`review/REVIEW_round62_D3bi.md`). The referee regenerated the md byte for byte, ran the 38 tests, and reran R
(three gaps, the tolerance scan's cell, the mirror to 3·10⁸ and 3·10⁹), all 14 RP rows and a C4M loop bitwise; re-derived the relative
mode's identity and the rest's shift with W (checked against a long-double Newton solve for the rest); moved the window's edges (a right
edge, two near edges, two sample spacings) to test the pin law; continued the relative mode in log-scaled form below the floating-point
floor on a window with every edge 11.7σ away; compared the law's rate with the full Jacobian's eigenvalues; rescaled s_d by the windows'
sample counts; ran the far rule at both boundaries with two new loops; ran 36 mutants; searched the literature.

**Verdicts (D3bi).** The pin law is right about the pin — the window's edge moves the block's rest with W, and the pins follow it at
every edge and spacing tried —, but (v) drew three wrong conclusions from it. In the continuum the near tie does decide a formed
runner-up's side (the tails' rule; the offset falls to 10^{−19} … 10^{−55}σ before it grows, so Part C's pin overrides it by 5 to 40
decades); the law's rate was about twice the true one (the midpoint is not the block's rest: with the centroid slaved the rate is κ =
W²/(2(w² + w_e²)) ≈ 0.51–0.53 times the coefficient), so RP's law failed at x = 6.04 and 7 for that reason, not the layer's; and the
window comparisons offered for s_d were mostly a rescaling of time (⟨·⟩ a mean over N samples), with s_d of held runner-ups the hold's.
Minors: the nearest edge decides formed runner-ups' outcomes; (iv)'s exception; collapsed against held; the pins' quantisation; numbers
and text (2.36·10⁷, Λ₀, a garbled proof ending, s_θ's digits); C4M at 10⁻¹⁷ and the far rule's boundary; 21 of 36 mutants surviving;
the literature (Chizat's WFR flow has no trailing); stale raw fields.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bi-1 | major | In the continuum the near tie decides a formed runner-up's side | (v): the continuum's forcing is only the centroid's relaxation, the mode keeps the layer's sign as a tail's — its offset falling to 10^{−19}σ (x = 5.48), 10^{−25}–10^{−28}σ (4.9), 10^{−45}–10^{−104}σ (3) —, so whatever moves the rest with W by more decides (on samples the nearest edge, without one the rounding); RP: a window with every edge 11.67σ away (`WIN_FAR`) and the law continued below the floor (forced by the rest, unforced, forced by the centroid's measured relaxation): unforced toward on C2's target and far on the mirror at every gap, 10^{−103.5} … 10^{−7.8}σ; at x = 6.04 and 7 the float64 runs resolve the offset and take the same sides; R: s_d at 5.48 the pin's |
| D3bi-2 | major | The relative mode's rate mislabelled; RP's law integrated with it | (v) and its proof: the slaved-centroid law dD/ds = κ[Λ₀(2Gm ww_e/W)⟨r, ∇²_cφ⟩D + (2(w − w_e)/W)(ρ′D + dc_b/ds)], κ = W²/(2(w² + w_e²)) — the rest the centroid's, F(c_b) = 0, not the midpoint's —, the pin unchanged (κ cancels); `relative_law` integrates it (`slaved`); RP: the Jacobian's eigenvalues against κλ_D·eig: −0.1351, −0.1346 against −0.1373, −0.1368 at x = 4.9, s = 449 (round 58's −0.260, −0.259), the late growth to the digits shown; the law at x = 6.04 and 7 now at 0.904–0.966 and 0.877–0.988 of the run (round 58's 0.004–0.865 and 0.000–0.053); test (`test_the_relative_mode_with_the_centroid_slaved`) |
| D3bi-3 | major | s_d's attribution: a time rescaling and the hold | RP prints N and s_d × N₀/N (N₀ = 1681): at x = 3 the top edges' change is the rescaling alone, −1.3% and −9.0% beyond it at 4.9 and 6.04; (v): s_d is the hold's where there is one (x ≲ 5.07, pin-independent) and the climb's where there is none (x ≈ 5.5–6); the continuation's climb in Part C's time — at x = 5.48 1.7·10⁴ later from the continuum's offset than from Part C's pin |
| D3bi-4 | minor | The nearest edge decides formed runner-ups' outcomes | RP's right-edge rows (far on C2's target, toward on the mirror); C4W: C4's loop at (4.9, 10⁻¹⁰) on that window ends bound on the far side (C2's target) and captured (mirror); (v), (iv) and the capture set say so |
| D3bi-5 | minor | (iv)'s exception | the side's decider named; the formed far-side runner-up's −1 once s₂ > s_{d,e}; "closes" replaced by the interval ((τ₂ − τ_e)/s_{d,e})² < ε < ((τ₂ − τ_e)/s_θ)² — 4.9·10⁻¹⁷ < ε < 2.9·10⁻¹² at x = 7 (R prints it) |
| D3bi-6 | minor | "Held wherever the layer collapses" | (v): collapsed for x ≲ 8.6, held — a rest above the fold — only for x ≲ 5.07 |
| D3bi-7 | minor | The pins' resolution | RP prints the runner-up's c_y's last place (1.85·10⁻¹⁶σ): x = 3's pin 23 steps, the law within 2 |
| D3bi-8 | minor | Numbers and text | R's s_d 2.346·10⁷; Λ₀ in the identity and its proof; the garbled proof ending removed; s_θ to two digits ("read on the output grid, within a factor 1.04") |
| D3bi-9 | minor | C4M at 10⁻¹⁷; the far rule's lower boundary | C4M's outputs from s = −100; two loops at the lower boundary — [2, 1, 1] at 6.5·10⁻¹⁷ (s₂ = 0.87 s_{d,e}), [1, 1, 1] at 3.7·10⁻¹⁷ (1.15 s_{d,e}) —; and the cause of the 10⁻¹⁷ cell (beyond the report): the loops at ε ≤ 6.5·10⁻¹⁷ ran at the absolute tolerance 10⁻¹⁴, above the masses, and entered thousands of s early (the entrant fully formed by s = −2000); `c4_atol` sets 10⁻⁵ε below ε = 10⁻¹⁴ (the earlier rows unchanged), and the rerun's closest approach is R's 1.8·10⁻⁸σ at every mass, the outcomes unchanged |
| D3bi-10 | minor | 21 of 36 mutants survive | golden job tests: `race_job` at the near tie's x = 4.9 and the mirror's long run at x = 7 reproduce the raw (s_d, the share, s_swap, s_{d,e}, the destination); `rp_job` at 4.9 and 6.04 (the law, the Jacobian's rates, κ); the continuum's sides (`test_in_the_continuum_the_near_tie_decides_the_side`); `centroid_velocity`'s extrapolation; `c4m_outcome` with the side 'toward'; the windows' edges; the continuation from above its climb; `c4_atol` |
| D3bi-11 | minor | Literature | Chizat's conic particle gradient descent and Yan–Wang–Rigollet's WFR flow: positions' mobility not ∝ weight, so no trailing and no pin; Denoyelle–Duval–Peyré (clustered positive spikes) and Batenkov–Goldman–Yomdin (near-colliding sources) beside the hold's moments; New recast as a consequence of the paradigm's metric |
| D3bi-12 | cosmetic | C4/C4b raw rows carried e_t2/w_t2 | stripped on load and save |

**Beyond the report.** The continuation as first written stopped at its start wherever |D(400)| exceeded the climb's end — on the far
window, whose layer runs at N₀/N of Part C's speed, at x = 5.48, 6.04 and 7 —; it now stops only after |D| has been below the end, and
RP was rerun. The early entry at 10⁻¹⁷ (D3bi-9) was an integration artefact of the loop's absolute tolerance, found with a probe of
the gap through the layer at 10⁻¹⁷ against 10⁻¹³ and confirmed at atol 10⁻²².

**Lesson.** A law validated at a fixed point says nothing about its rate: the pin agreed to 1% while the coefficient was twice the true
one, because forcing and rate scale alike there. A selection decided at 10⁻¹⁹σ is decided by whatever is larger — here a window edge
— and the clean statement is the continuum's, which only a continuation below the arithmetic can make. And a solver's absolute tolerance
is a scale: at masses below it the state is not integrated, whatever the relative tolerance says.

## Round 64 (adversarial review of Proposition D‴.52 as revised after round 59, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bk: Proposition D‴.52 as revised after round 59, with its code and results
at 8c826c6 (`review/REVIEW_round64_D3bk.md`). The referee:
- regenerated the md byte for byte and reran all of S, 18 F jobs and two A and B rows bit-identically;
- identified the ISS's knots with the LARS and lasso knots on 300 noisy targets, and computed the separation's Cramér–Rao bound with the
  amplitudes unknown and known;
- checked z's sign and H1 mean (1 − cos = ϑ²/12), both readings of (iv) by finite differences, and z's null against s/|ŵ|;
- analysed F's raw per wall (the one-sided p_A, p_B), the α allocation, the walls' power at small Δ, the moving null and the sizes at
  α = 0.01, and decomposed E's knot loss through the runner-up's own knot;
- probed (0.8, −0.6) and (0.8, 0.2) with new seeds, ran 69 single-point mutants and searched the literature.

**Verdicts (D3bk).** The deterministic core, the null law, the score at the bound and the placement at the Cramér–Rao bound hold, and
(iii)'s loss mechanisms are right. But the knot is a de facto one-sided test: the grid's margins put 75–100% of its null rejections at
the wall of the fit's sign, which explains both why it rivals the sign-agnostic |z| for same-sign pairs and why it is biased for
opposite signs. The draft compared the two-sided walls test with the one-sided score, and its title claimed the bound without the sign
condition. Testing both walls at α/2 does not remove the bias: their power stays below α at small Δ in all 24 opposite-sign cells.
Minors: Smith's criterion with unknown unequal amplitudes sits at s ∝ Δ³; τ_k = 1/λ_k for k ≤ 2 only; z's mean has V's sign and the
brackets clash; the non-monotone power is the null moving across the tie; the calibrated knot is conservative more broadly; the lower
quantiles' wording and the unmeasured decomposition; 25 non-equivalent mutants surviving; (iv)(b) is the classical splitting criterion.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bk-1 | major | The knot is de facto one-sided; the comparisons with the score ignore it; the title's bound needs the sign | Title: "at the Neyman–Pearson bound with the pair's sign known, and at 1.04–1.18 times the bound's separation without it"; (iii): the α allocation and like-for-like comparisons — the fit's-sign wall one-sided 1.056 [0.987, 1.289] times the bound's Δ₅₀ against z one-sided 1.003, the full knot 1.064 against \|z\| 1.085, for (0.8, −0.4) the other wall one-sided 1.680 [1.621, 1.749]; F prints the per-wall share of the H0 rejections (390/430 and 328/364, 75–100% per offset) and the one-sided per-wall tests; New reworded |
| D3bk-2 | major | Both walls at α/2 stay biased | (iii) and F: "reaches power ½ sooner but stays biased; only the score is unbiased for both signs" — least power 0.010–0.045, below α in 24 of 24 cells, pooled 0.036 [0.031, 0.042] and 0.039 [0.034, 0.045] at the two smallest separations, with the mechanism (the fit's-sign wall's one-sided power 0.012 and 0.005 there, the other wall's 0.057 and 0.080); F prints each test's least power and its pooled power with Clopper–Pearson intervals (`f_alloc_lines`) |
| D3bk-3 | minor | Smith's criterion mis-derived | Status: with unknown unequal amplitudes the separation's bound ∝ s/Δ² (the split moves the separation by (Δ/2)(1/w₂ − 1/w₁) per unit of weight), the criterion at s ≈ 474Δ³ and 56Δ³, (v)'s order; with equal or known amplitudes ∝ s/Δ (49.4Δ², 51.7Δ², 109Δ²), (iii)'s; B prints `smith_line` (`crb_sep`) |
| D3bk-4 | minor | τ_k = 1/λ_k for k ≤ 2 only | (iii) and the proof: the first two knots; later ones while the block's least-squares fit keeps the duals' signs and the lasso drops no atom (the identity carried from piece to piece) |
| D3bk-5 | minor | z's sign; the brackets; D‴.54's time | (iv): the mean has V's sign, ŵz has sign(A₁A₂); proof: cos = (1 + ϑ²/6)^{−1/2} = 1 − ϑ²/12 + O(ϑ⁴) in the continuum; (·, ·) and ‖·‖ Euclidean, ⟨·⟩ the mean; D‴.54's time t; tests: z against an independent projection and the pair (−0.8, −0.6) |
| D3bk-6 | minor | The non-monotone power is the null moving across the tie | (iii): the null's fit at ĉ − c̄ = −T₃/(4Aσ²) + … ∝ Δ³ crossing the tie; F prints per separation at c̄ = 0.032 the null's ĉ − c̄, the law's quantile × s and the H1 knot's median × s (s = 0.01: +0.0097, +0.0217 at Δ = 0.05, 0.07; 4.83 → 1.97; 4.08 → 2.05), and the nulls' offsets at the H0 separations |
| D3bk-7 | minor | The calibrated knot conservative beyond the tie | (iii): "conservatively, the more so at α = 0.01 and near a tie"; F prints the sizes at α = 0.01 (the calibrated knot 0.0076 [0.0064, 0.0089]; z one-sided 0.0086 and 0.0103; \|z\| 0.0090) |
| D3bk-8 | minor | The lower quantiles offset-dependent; the decomposition unmeasured | (iii) worded by offset; E prints the adjacent share among the earliest 5% (0.34–0.52 at c̄ ≤ 0.02, 0.76–0.80 at 0.032, 1.00 at 0.045–0.049); Part U: the runner-up's own knot against its own law (size 0.039–0.068; z₂'s power to −0.015 … +0.008 at c̄ ≤ 0.032, the dictionary's minimum 0.012–0.098 below it; at c̄ = 0.049 0.665, 0.247, 0.165 against z₂'s 0.705, 0.400, 0.530); E′ prints the margin over its sd (6.94, 2.12, 0.89 at c̄ = 0.049) |
| D3bk-9 | minor | 25 non-equivalent mutants survive; report sections untested | tests (33): the walls for both signs, `s_job` on an E cell (labels, stored fits, z₂ toward either wall, swap, margin, cosine) and on a C null (the law's seed and noise), `part_S`/`part_U` seeds and fits, `f_job` under H0, z against an independent formula, z_σ ≈ sign(ŵ)z, `e_job` calibrated at the path's fit, `knot_of_kernel` and `u_job`, E′'s U column, the separation's bound, E's earliest-5% column and S and U passed on to E′ |
| D3bk-10 | minor | (iv)(b) is the classical splitting criterion; Shahram–Milanfar's test; Neyman–Scott 1965 | (iv): Fukumizu–Amari 2000, Liu–Wu–Wang 2019's splitting matrix, Shahram–Milanfar's locally optimal test; (b) needs a merged block, which the transit of (i)–(iii) does not form; Status: Neyman–Scott 1965; New: "the classical splitting criterion in the paradigm's metric" |
| D3bk-11 | cosmetic | Misquotes; A's round-off; z's range | "four times the bound's 0.01001"; the margin's noise from E′'s new column; Part A rerun with the current solvers (equal to round 59's but for the rms at 10⁻¹⁶); z's sizes in F's range printed |

A plumbing fault found while applying: `report` did not hand U to E′, so the runner-up's-knot column printed "—" from the raw; fixed,
and a test now spies on the call.

**Lesson.** Where a test spends its α is part of what it is. The knot never chose its allocation — the grid's margins put it at the
wall of the fit's sign —, so it is a one-sided test, as good as the sign-agnostic score for same-sign pairs and blind to opposite
signs. A remedy is judged by the defect it was meant to cure (power at least α at small Δ), not by a different figure of merit
(Δ₅₀); and comparisons are fair only like for like: one-sided against one-sided, sign known or not.

## Round 61 (adversarial review of Proposition D‴.59 as revised after round 57, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bh: Proposition D‴.59 as revised after round 57, with its code and results at
a1f5d85 (`review/REVIEW_round61_D3bh.md`). The referee:
- regenerated the md byte for byte and reran A, D, H and E cells field by field;
- dumped the raw in full precision (the ends' sub-sample and negative widths), followed C's pair targets block by block and H's
  equal-amplitude path at the tolerances 10⁻⁹ and 10⁻¹³;
- re-derived (i), (ii)(a) and (ii)(c), checked the width split's order on same-sign, cancelling and the loop's large-|w| pairs, and
  re-fitted the optimal dipole;
- scanned the one-sign separation at σ* = 0.06 and ran the real loop on G's six targets at ε = 10⁻⁵ … 10⁻¹⁰ (the dipole phase, the
  landing, the ε-scaling, the collapse) and on H's targets to τ = 4.1·10⁵;
- ran 49 single-point mutants and searched the literature.

**Verdicts (D3bh).** (i), the dipole's runaway law and the width split hold. But three of the draft's words were thresholds: (ii)(b)'s
"exact" was the path's stopping tolerance (with 10⁻¹³ the path keeps entering kernels and one feature is approximated by an eight-kernel
cluster, with unreported sub-sample and negative widths), σ* = 0.04's "coinciding pair" a null triple beside an exact fit, and the
loop's "migrates … lands" a midpoint and PAIR_TOL — the pair's centre of mass already sits on a feature; landed, it collapses
algebraically; and the landed half-difference follows a law the draft never stated, v ∝ ε^{−1/6}, the block's runaway in the loop's
fast time. The cancellation's condition "closer than σ₀" was contradicted by its own evidence (1.33σ₀) and the loop does not keep it.
Minors: (ii)(a)'s check selective and its constant non-uniform; the pairs' mechanism; F's pair drifts; (iii)'s noise kernel; 17
non-equivalent surviving mutants; raw provenance; (i)'s proof's order; σ unconstrained; literature and "New".

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bh-1 | major | (ii)(b)'s "the fit exact" is the path's tolerance; H's ends misdescribed | (ii)(b): the path reaches its tolerance with one feature fitted by one or two kernels and the other by a cluster of seven or eight with the feature's mass, and at a lower tolerance keeps entering kernels; H2: the six σ* = 0.06 paths to 10⁻¹³ (20–32 entries, the last at τ = 1.16·10¹² … 2.34·10⁶⁷, the residual creeping to 1.5·10⁻¹¹ … 2.0·10⁻¹⁰, the masses 0.9998–1.0002); H's columns for the sub-sample and negative widths over all active kernels (3–6 and 0–2); `outcome`: a pair counts only if its net coefficient and centre of mass fit a feature (any feature: D's two share a centre), else a null cluster, and 'exact, null cluster' (σ* = 0.04) |
| D3bh-2 | major | (ii)(c)'s real-loop sentence | (ii)(c): the dipole as a phase — v³ = 3μ₀KΛ(ε)(τ − τ_d) until the next kernel forms, the landing (centre of mass and W at a feature, within 1%), the algebraic collapse —, v ∝ ε^{−1/6} derived from the law and Λ(ε) ∝ ε^{−1/2}; G at ε = 10⁻⁵ … 10⁻¹⁰ (`phases`, `g_phase_summary`): the slope over the law 0.997–1.112 at 10⁻¹⁰, the landed v 0.71–1.09 of the law's, the last decade's ratio 1.4078–1.5002 against 10^{1/6}, landed at τ = 6400 on 0, 3 and 6 of 6 at 10⁻⁵, 10⁻⁶ and ≤ 10⁻⁷, the collapse τ^{−0.235 … −0.173} (distance) and τ^{−0.906 … −0.689} (residual), 6 cells frozen; the centre of mass's offset and W − A in place of the midpoint, the pair column at report time |
| D3bh-3 | major | The cancellation's condition contradicted; (b) the path's | (ii)(b): at 1.33σ₀ with the amplitudes within 10%; H3 (the separation scan): 0.83σ₀ stalls with 10 negative kernels, 1.67σ₀ and 2σ₀ exact with kernels of one sign, 2.5σ₀ and 3.33σ₀ a coinciding pair; HL (the loop on H's equal targets): σ* = 0.04 ends with two one-sign kernels at τ = 6400, σ* = 0.06 still has 6–12 negative kernels there |
| D3bh-4 | minor | (ii)(a)'s check selective; the constant not uniform | (ii)(a) and the proof: the expansion in \|w_i − w_j\|Δc/(\|W\|σ) (the k-th moment about c̄ written out), valid while \|w_i − w_j\|Δc ≪ \|W\|σ; the check on the three pairs a little apart and, new, on the loop's 72 landed rows (0.9828–1.0458 of the prediction, the parameter ≤ 0.284); at coincidence (C's Δ ≥ 0.6 target 2, Δc 2.78·10⁻⁸) unresolved |
| D3bh-5 | minor | The pairs' mechanism | (ii)(a): two early kernels approximating the widest feature collapse onto it together, or the last entrant, of the other sign, cancels an early kernel's excess |
| D3bh-6 | minor | F's pair drifts apart; (iii)'s noise kernel | F to 1.3, 6, 30 and 150 times the last entry at 10⁻⁶ and 10⁻⁸: the pair drifts apart at 10⁻⁶ (7.0·10⁻³ → 1.6·10⁻², the residual at the idle kernels' floor 5.0–5.4·10⁻⁶), stays 2.6·10⁻³ apart at 10⁻⁸; (iii): the stop's fourth kernel fits noise, no pair forms |
| D3bh-7 | minor | 17 non-equivalent surviving mutants; the job functions untested | tests (17): `outcome` at its thresholds (the stalled branch included), the pair rule with formed kernels and null clusters, a pair fitting any feature sharing its centre, `closest_pair` with unequal widths, `dipole_law` (the nearest opposite kernel, E_inf's ½), `h_stats`, the stalled block's summary, `phases`, `oracle_error` started at the truth |
| D3bh-8 | minor | Raw provenance; C's table | A–D regenerated with the current code; `oc` is `outcome` on C's own end, the stalled block printed; the polished \|w\| marked solver-dependent; G's pair column computed in the report |
| D3bh-9 | minor | (i)'s proof | E = J + ½‖φ‖²(w − w*)² and E ≤ ½‖f*‖² − η give (w − w*)² ≤ w*² − 2η/‖φ‖²: w keeps w*'s sign, away from 0 and below 2\|w*\|; "on the critical curve" |
| D3bh-10 | minor | σ unconstrained | stated in the setting; H and H2 count the widths below the samples' spacing and of negative sign |
| D3bh-11 | minor | Literature and "New" | Status: Gallon–Jentzen–Lindner 2022, Kranz–Gallon–Dereich–Jentzen 2025; the exponent from the mobility (v^{a+2−b} ∝ s; Chizat 2022); Pearson's moment algebra and Runnalls' merge for the width split; Schiebinger–Robeva–Recht 2018 for the cancellations; the sliding Frank–Wolfe algorithm; New narrowed to the loop's three phases and the ε^{−1/6} landing, the split and the cancellations as observed |
| D3bh-12 | cosmetic | Misquotes | 2.5·10⁻⁴ on target 12; B's nine rows one computation (the grid scales with σ*); E's agreement 7·10⁻⁸ (printed) |

**Lesson.** A name given to a configuration has to survive the configuration's invariants. "Exact" was a tolerance, "coinciding pair"
a threshold on a null triple, "lands" a threshold on a distance still falling as τ^{−0.2}, "migrates" a midpoint that unequal
coefficients displace. Measured instead — the centre of mass and net coefficient against the feature, a lower tolerance, a longer
horizon, a scan of the separation —, the loop's dipole turned out to follow the block's own runaway law in the loop's time, which
predicts the landed coefficient's ε^{−1/6}.

## Round 63 (adversarial review of Proposition D‴.60 as revised after round 56, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bj: Proposition D‴.60 as revised after round 56, with its code and results at
ada0e1c (`review/REVIEW_round63_D3bj.md`). The referee:
- regenerated the md byte for byte and reran B4 and C/C2 cells bit-identically;
- re-derived (i) and the proof of (ii), computed the linearised loop's monodromy (the Floquet multipliers) and checked the second
  order's first order against B4's odd part;
- ran the loop to 1200 periods at 21.4 and with the 16 tones, at further tones, and the averaged loop at 13.05 and 13.15 beside 13.1;
- ran 149 ejection controls (round 56's ejecting pair, C2's pairs under four seeds, single tones, integer-multiple tones, WFR to
  s = 0.32, eight more seeds, m between 3.5 and 5.5) and decomposed the band control into roll-off and resonance;
- ran 51 single-point mutants and searched the literature.

**Verdicts (D3bj).** The first order, the oscillator law and the exact second order hold. But (ii)'s "bounded unless commensurate, where
secular" was false — the linearised loop is exponentially stable (multipliers 0.99400, 0.99245, 0.99237: every coherent sum saturates
with the slow modes' τ) —, the averaging failure at 13.1 is its commensurate 2ω part (13.05 and 13.15 cross the same resonance and
average well), the 16 tones' saturation was a 300-period snapshot, and "single tones and tone pairs do not eject" fails: round 56's
pair still ejects. The resonance attribution was not isolated (most of the band effect is roll-off) and the bound on m is no
threshold. Minors: sources and "New"; single-seed controls; "after hundreds of periods" and "first passage"; the WFR margin; the 0.3
cut counting slow excursions; B3's asymptotics; 21 behaviour-changing mutants surviving; the peaks below √K_i; two wordings.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bj-1 | major | Single tones and pairs do eject; "multi-tone" unsupported | (iii): single tones drift but never jump (the largest distance ≤ 0.363 at s = 0.08), pairs and bands jump — a pair one ω_ref apart from s = 0.04, D‴.49's own pair with its draw's phases at s = 0.02; C2: single tones at 10, 13.1, 15, 15.7, 16, 21.3 and pairs 17.5, 17.75, 18ω_ref and one ω_ref apart under four seeds, and D‴.49's own pair, at s = 0.02–0.08; jumps (beyond 0.6) reported apart from leaves |
| D3bj-2 | major | "Bounded unless commensurate, where secular" false | (ii) and the proof: the period means as Floquet sums, \|m_j\| = \|δ\|\|ρ^j − e^{iΩjT}\|/\|ρ − e^{iΩT}\|, at commensurability an offset \|δ\|/(1 − ρ) reached with τ; F: the monodromy's multipliers against m and the metric (0.99400, 0.99245, 0.99237 at m = 12; τ = 166, 132, 131 periods), the Floquet recursion against the loop's odd part (3.38·10⁻³ … 2.40·10⁻² of its largest, ∝ s² relative), the commensurate tones' saturation (21.4: −0.5230 against the loop's −0.521 at period 1200) |
| D3bj-3 | major | The averaging domain misattributed | (ii): the averaged loop's error is the dropped parts' Floquet sums, large where 2ω_k or ω_j ± ω_k is commensurate; B5: at 13.05, 13.1, 13.15 the constant part (the half-sum at φ and φ + π/2) is +154, +153.6, +153.1, the 2ω part −1.02, +69.84, −1.584, the averaged loop 0.959–1.085 of the constant part at all three |
| D3bj-4 | major | The 16 tones' saturation a snapshot | B4 to 1200 periods: the single tones saturate (13.1's first coefficient +1.13·10³ from period 900), the 16 tones stay within [+235.7, +290.1], [−93.63, −87.06], [+123.4, +139.1] over the periods 601–1200 |
| D3bj-5 | major | The resonance attribution not isolated; the bound no threshold | A2: the resonance factor (the jitter over its roll-off alone) 1.088–1.495 from m = 2 to 12, 1.225 at the bound; the band moved to 25–75 cuts the jitter 2.92, 2.61, 2.79-fold, of which the roll-off 2.24, 2.21, 2.22 and the resonance 1.31, 1.18, 1.25; (i): the bound "reduces" the resonance; C2's m scan under four seeds: the first jump at s = 0.08, 0.06, 0.06, 0.06, 0.04, 0.04, 0.03 for m = 3.5 … 8 — smooth |
| D3bj-6 | minor | Novelty and sources | Status: the Floquet response, averaging with periodic signals and k·ω = 0 classical; Riedle–Cyr–Kokotović 1984, Kokotović–Riedle–Praly 1985, Fu–Sastry 1987, LMS's non-Wiener shifts; New limited to the instance (Gmw as the stiffness, the exact second order with centres as Floquet sums, the ejection data) |
| D3bj-7 | minor | Single-seed controls | C2 under the seeds 7, 1, 2, 3 (the draw's pair with its own phases) |
| D3bj-8 | minor | "After hundreds of periods", "first passage" | (iii): a drift on the slow modes, then a jump; C3 (eleven seeds at the threshold): the first kernel leaves at the periods 114–274 (τ = 131–166) with its coefficient risen 1.071–1.159-fold — (ii)'s s²⟨w₂⟩ at saturation, +0.094 … +0.116 against +0.064 … +0.143 — its centre still at the feature, then jumps beyond 0.6 within 0–1 periods; C prints the leave periods per level against τ |
| D3bj-9 | minor | The WFR margin | C2 to s = 0.32: jumps from 0.12 at k_e = 5, none to 0.24 at k_e = 10 (one at 0.32) |
| D3bj-10 | minor | The 0.3 cut counts slow excursions | C, C2 and C3 report the jumps (the largest distance beyond 0.6) apart from the leaves |
| D3bj-11 | minor | B3's asymptotics | B3 prints the remainders per doubling of s (cubic and quartic but on the 16 tones' first kernel); the caption "mostly ∝ s²" |
| D3bj-12 | minor | 21 behaviour-changing mutants survive | tests (32): `b4_job` on two periods (even_c included), `c_noise`/`c_job` options reaching `simulate` (seed, band, freqs, metric, k_e), C2's labels, `monodromy` against the loop, the Floquet means, `b5_job`, `c3_job`'s statistics |
| D3bj-13 | cosmetic | "Peaks near √K_i" | (i) and A2: about 10% below √K_i (0.898, 0.871, 0.887 at m = 12) |
| D3bj-14 | cosmetic | Wordings | the rectification's sign reversal dropped; the slow modes in the plural |

**Lesson.** A dichotomy ("bounded unless commensurate, then secular") is checked by the one computation that decides it — here the
monodromy's spectrum, which makes every coherent sum saturate with the slow modes' time constant. A control is informative where the
effect was seen: rerunning the pair that had ejected showed that it still ejects, and the special tone's ordinary neighbours placed the
averaging failure in its commensurate 2ω part. And a leave is an event with an anatomy: measured, it is the second-order offset built
up on the slow modes, then a jump.

## Round 65 (adversarial review of Proposition D‴.58 as revised after round 60, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bl: Proposition D‴.58 as revised after round 60, with its code and results at
3eb9326 (`review/REVIEW_round65_D3bl.md`). The referee:
- regenerated the md byte for byte, reran four full B paths, the sixth variant's paths to their first non-optimal blocks and all 14 B2
  jobs bit-identically;
- restarted the solve from the range ends with one setting changed at a time (budget, gradient test, step test, Jacobian scaling, start);
- computed the analytic variable-projection Jacobian, the condition numbers in several units and the smallest singular values;
- recomputed the Mercer frame's defects and the first tail ratio in every Mehler frame;
- deduplicated the entries across the variants and measured the stationarity test's effective threshold;
- ran 52 single-point mutants.

**Verdicts (D3bl).** (i), (ii), (iv), (v) and every number in (iii) and the numerics but the ones below check out. But the sixth variant's
two new range ends (τ/σ = 1/4, K = 6; 3/4, K = 11) lie below the arithmetic floor, where the joint solve stops on its step test 1.5% and
0.4% above the optimum whatever the budget — not on the budget, radius or gradient test the draft named —, and 19 single-point mutants
survived both the tests and `--from-raw`, among them `b_job` dropping the sixth variant's gradient setting (X01). Minors: B2's far misses
are variable projection's stalls at coalescing nodes; the condition numbers depend on the unit and do not explain variable projection's
success (the smallest singular values agree); the Mercer sentence names the wrong rule and the wrong measure; the entries were counted
across the variants (316 entries are 107 events); the stationarity test is weak beyond the vacuous blocks; the status-1 wording (the
"exception" block's polishes reach the optimum and the radius rejects them); a floor qualifier in the status.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bl-1 | major | The sixth variant's two range ends lie below the floor | B3 (new): each path's first non-optimal block restarted by the polish's own solve — with the path's budget and gradient test (and the path's seed rule on its end), and without the test at 5·10⁴ and 5·10⁵ evaluations; above 10⁻¹⁰ (20 ends) the solve reaches the optimum from all 20 with 5·10⁵ (1530–227797 evaluations, the step test; 16 with 5·10⁴), below it (2 ends) it stops on the step test 1.0147 and 1.0037 times above it; (iii): above the floor the budget, the radius or the gradient test binds, the two ends below it censored by the floor; status |
| D3bl-2 | major | 19 surviving mutants (X01 the sixth variant's gtol; X02 dc; X04, X05, X07 the condition numbers; R01 the polish's default; X28–X30 the `--variants` CLI; boundary mutants) | Tests: `b_job` passes each variant's polish settings (monkeypatched `sliding_iss`); `polish_block`'s defaults pinned; `b2_job` takes dc from the better start; the analytic VP Jacobian against finite differences and the condition numbers on a 3-kernel case; `main --variants` merges the new paths (stub `b_job`); B2's line and the entries at their boundaries; B3's restarts and seed rule; `b3_lines` on synthetic data; the frames and the first tail ratio |
| D3bl-3 | minor | B2's far misses are VP stalls | B2's line: the far misses' closest nodes (0.000–0.021σ apart), the largest residual–Jacobian cosine at their end (2.6·10⁻⁴–4.7·10⁻³) and their stops; (iii): "stalls with two nodes coalescing at an end that is not stationary" |
| D3bl-4 | minor | Condition numbers: units, accuracy, not the mechanism | `vp_jacobian` (Golub–Pereyra, analytic); the centres in units of σ: joint 4.5·10⁶–2.3·10⁹, VP 5.4·10³–3.8·10⁵, the smallest singular values' ratio 1–31.7 (median 1.01); (iii): the two problems differ at the stiff end, not (at most blocks) in their slow modes, and the polish's own solve reaches the optimum given the evaluations (B3) — "far better conditioned" dropped |
| D3bl-5 | minor | The Mercer sentence | A5 (new): the Mehler frames ρ = t, q, λ — the frame's own rule N(0, v_ρ), v_ρ = 2σ²τ²/(2σ² + (1 − ρ)τ²), v′ = 2, 1.029–1.5, 1; the first tail ratio −v′K(K + 2) + (v′ − 1)(K + 1)(2K + 1) exact in 60 of 60 (`gh_defect_v`), of degree 1 only at v′ = 2; the Mercer frame's rule N(0, 2τ²/(β² + 1)) at 63.49, 17.17, 5.321, 2.332 per kernel; status rewritten (the frames, the Mercer frame's rule, the Mercer-zero rule's measure); the √2 bound in the normalised space, (2σ*/σ)^{1/2} in relative L²; Kuo–Woźniakowski's condition unresolved |
| D3bl-6 | minor | Entries counted across variants | `entry_lines`: one event per width, grid, iteration and kernel whose duals agree within 10⁻⁴ (their spread 6.8·10⁻⁶), events (entries) per row, the 6 events in two rows noted; (iii): 1 of 107 events of the largest rate, 105 of 107 beyond the outermost node (median 1.17σ) |
| D3bl-7 | minor | The stationarity test weak | B's summary: the threshold over the Cauchy–Schwarz bound for the 63 non-vacuous blocks, 3.5·10⁻⁴–0.97 (median 7.7·10⁻³); (iii) |
| D3bl-8 | minor | Status-1 wording | (iii) and numerics: the gradient test "below 10⁻¹⁵"; the once-variant at 1/4 and the third repeated polish at 1/2, 3/4; at τ/σ = 1, the grid shifted, K = 9, the polish stops on the test at the optimum and the radius rejects it (B3: 1.27σ); B's table prints the repeated polishes' ends over the optimum in brackets |
| D3bl-9 | minor | Floor qualifier | Status: "the first non-optimal one whose optimum is at least 10⁻¹⁰" |
| D3bl-10 | cosmetic | "Nodes up to 2.14–3.99σ" | Sentence replaced (D3bl-3) |

**Lesson.** A range's end is a property of the solver only where the solver could have gone further: below the arithmetic floor the
solve stops on its step test wherever its budget, gradient test or radius stand, so naming those as the binding settings there was a
claim about the floor, not the solver. Restarting the solve from each end with one setting changed at a time — the cheapest experiment
on the table — separates the two, and a mutation pass on the variant plumbing (here the gradient setting a refactor dropped silently)
belongs with every new variant.

## Round 67 (adversarial review of Proposition D‴.52 as revised after round 64, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bn: Proposition D‴.52 as revised after round 64, with its code and results at
88933bf (`review/REVIEW_round67_D3bn.md`). The referee:
- regenerated the md byte for byte and reran all 15 U jobs, all 95 S jobs, five F jobs and three A rows bit-identically, and wrote U
  again from scratch;
- checked the ISS against LARS and the lasso on 400 noisy targets, re-derived (iv)'s cosine and the moving null;
- compared the tests like for like with paired bootstraps, decomposed the knot's loss with each knot against its own law, froze the
  runner-up's margin, recomputed E′'s margin uncensored and Smith's criterion by exact roots;
- ran 71 single-point mutants and ported round 64's surviving ones.

**Verdicts (D3bn).** Reproducibility, the LARS identity, (iv), the de facto one-sided knot, the walls' bias, Part U and the moving null
hold. But the two-sided half of the like-for-like comparison was not like for like: the full knot, de facto one-sided, was set against
the two-sided |z|; both walls at α/2 trail |z| in 23 of 24 cells, and a score spending α the knot's way reaches ½ at 1.02 × NP against
the knot's 1.06. Minors: the minimum's part changes sign at the tie; the Status equated Liu–Nehorai's GLRT limit with Smith's criterion
where they part; E′'s margin column was censored at the tie; 8 of round 64's non-equivalent mutants and two new ones in U survived;
leading orders without their range; misquotes.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bn-1 | major | The two-sided comparison not like for like | `rejections`' 'alloc': the score spending α as the calibrated knot does (`knot_allocation`: 90.7% and 90.1% toward the fit's sign: ŵz > 1.692 or −ŵz > 2.601); F prints its Δ₅₀ (1.019 and 1.257 × NP) and `paired` comparisons per pair — the calibrated knot over the allocated score 1.045 (smaller in 2 of 24), both walls over |z| 1.026 (1 of 24), the fit's-sign wall over the one-sided z 1.043 (0 of 24), the knot over |z| 0.987 (17 of 24); (iii) and New rewritten: one-sided against one-sided, two-sided against two-sided, and with the knot's own allocation; New quantified (4–6% in Δ₅₀, about 20% in samples) |
| D3bn-2 | minor | The minimum's part changes sign at the tie | E′: the full knot against its own law at the truth (size, power) and the minimum's part per cell in all 15 — −0.1075 … −0.0300 at c̄ ≤ 0.032, −0.0075 … +0.0175 at c̄ ≥ 0.045 (positive in 5 of 6); (iii): a loss away from the tie, a small gain at it; the margin's noise enters the knot's numerator μ̂ |
| D3bn-3 | minor | Liu–Nehorai = Smith only with equal or known amplitudes; local-fit constants | Status: the identification restricted; otherwise Smith's criterion an estimation threshold, not a detection limit; `smith_line` prints the constants at the integer powers (501Δ³, 119Δ²) beside the fitted exponents |
| D3bn-4 | minor | E′'s margin censored | `s_job`'s 'mu_u' (the runner-up's dual over the noiseless first entrant's correlation, not the largest; S rerun, every other field bit-identical); E′'s margin column from it: sd(μ̂)/s 0.148–0.170, the tie cell's margin 0.71 sd |
| D3bn-5 | minor | Surviving mutants | Tests: `draw_stats`' walls in `walls_at_fit`'s order and `knot_walls` against an independent per-kernel definition (D01, K01); `f_job`'s qA at α/2 and the bound at the draws' s (F03, F04); `f_jobs`' H0 separation (F07); E′'s t size and z_σ power and the uncensored margin (Z02, D02); `width_z`'s escape window (S03); `u_job` at a swapping tie cell and at c̄ = 0 (U09, U01); U's first entrants read against S's swaps (U10); the allocated score, `knot_allocation`, `paired` |
| D3bn-6 | cosmetic | Leading orders without their range | (iv) and the proof: cos + O(ϑ⁴ + (Δ/σ)⁴); the moving null "to leading order"; the title "≈ 1.09 (1.04–1.18)" |
| D3bn-7 | cosmetic | Misquotes, loose ends | the opposite-sign knot over 23 cells; the walls mechanism with the walls test's own α/2 components (0.0034 and 0.0345 under H1 against 0.0223 and 0.0235 under H0); the calibrated knot's size per offset (0.0512 at c̄ = 0); z's tested range s/|ŵ| = 0.0021–0.0735; z₂'s size in E′; Neyman–Scott 1966 |

**Lesson.** Like for like has two halves, and a fix applied to one of them is not the rule. Round 64 matched the one-sided comparison and
left the two-sided one comparing a test that spends α on one side with one that splits it; the measure that settles it — a score spending
α exactly as the knot does — took one line of code. A decomposition is a claim about every cell, signs included: the part quoted only
where it was a loss turned out to be a gain where it was not quoted.

## Round 66 (adversarial review of Proposition D‴.54 as revised after round 62, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bm: Proposition D‴.54 as revised after round 62, with its code and results at
6d48964 (`review/REVIEW_round66_D3bm.md`). The referee:
- regenerated the md byte for byte and reran RP's and C4's rows bitwise;
- continued the relative mode in exact arithmetic about a merged block on the far window and injected it into the full race;
- reran C4's (4.9, 10⁻¹¹) at three absolute tolerances;
- ran 34 single-point mutants.

**Verdicts (D3bm).** Two majors. A held runner-up's s_d is not pin-independent: the hold starts only when the climb ends, so s_d shifts
one for one with the climb's end (in the continuum +3.8%, +12.2%, +23.4%, +38% at x = 4.5, 4.9, 5.07, 5.48), and C3's capture at
x = 5.07 belongs to the window. RP's continuation below the floor had the right sides but artefactual numbers: it was linearised about the
run's own base after the run had left the merged state, and its "measured relaxation" was the pair's own feedback, extrapolated (the least
|D| off by up to 54 decades, the climbs by up to 4·10⁴). Minors: C4's (4.9, 10⁻¹¹) an absolute-tolerance artefact (74 of 74 at
atol ≤ 10⁻¹⁵); the c4_atol test vacuous (`np.isclose`'s default atol) and 17 of 34 mutants surviving; the rate check printed only where
best, κ up to 0.62, the slaving condition unstated; "in the continuum" meaning the ε-free problem in exact arithmetic; RP's labels and a
NaN; the overlap-singularity literature; slips.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bm-1 | major | s_d = the climb's end + the hold, not pin-independent | (v): s_d is the climb's end (\|D\| = 0.01σ) plus the hold, the hold starting when the climb ends (the merged pair's split frozen); the climbs from Part C's pin 9.9·10³ … 5.09·10⁴ against the continuum's 6.97–7.49·10⁴; the continuum's s_d 2.352·10⁷, 1.257·10⁶, 3.58·10⁵, 1.854·10⁵, 8.43·10⁴ (x = 3 … 5.48); the capture set's ends the window's; C3 and R: x = 5.07's capture the window's, bound in the continuum; RC prints the climb's and s_d's shifts, equal to 10⁻⁴ |
| D3bm-2 | major | RP's continuation: the sides right, the numbers artefacts | New Part RC: R's race to s = 400 on the far window, then the merged block run on (`merged_rhs`) and the centres' perturbation evolved by the race's exact Jacobian about it (`PairSums`, the centre block of `critical_jac_nd`, midpoint matrix exponentials, the norm in logarithms), injected into the float64 race at 10⁻⁹σ after regrowth: the least \|D\| 10^{−98.7}, 10^{−36.0}, 10^{−27.6}, 10^{−24.7}, 10^{−19.0}, 10^{−13.4}, 10^{−7.7}σ (x = 3 … 7), the sides the layer's, the climbs, s_d; checks (16 substeps, the block's offset 0, the injection at 10⁻⁷ and 10⁻¹¹σ) and the validation at x = 7 (3.2725·10⁵ against 3.2727·10⁵); RP's continuation and `centroid_velocity` ("meas") removed |
| D3bm-3 | minor | C4's (4.9, 10⁻¹¹) an absolute-tolerance artefact | `c4_atol` = min(10⁻¹⁴, 10⁻⁵ε); C4, C4b, C4M and C4W rerun: the two tolerances agree at 30 of 30 points, R's rule at 74 of 74; the artefact's 7.0·10⁻¹⁵σ dropped; C4's and C4M's captions give the tolerance |
| D3bm-4 | minor | tests | `c4_atol` tested exactly; `PairSums`' force and block against central differences (float64, long double); `merged_rhs` against the race's right-hand side; `growth_rate`'s window; κ and the slaved law's error on the raw; RC on the raw (sides, x = 3's least and climb, the checks, the shift identity) and a golden `rc_job` (x = 5.48) |
| D3bm-5 | minor | the rate check selective; κ; the slaving condition | RP prints the Jacobian against the slaved law at every float64 row (κ 0.506–0.620; error 0.0–10.1% of the largest rate, 0.4–10.1% at s ≈ 450, the most at x = 3); (v), the proof and the Status state the slaving condition and separate the identity (exact to first order) from the slaved law |
| D3bm-6 | minor | "in the continuum" | (iv), (v), the proof and New: the ε-free inner problem in exact arithmetic; no float64 loop resolves the offset below x ≈ 6 |
| D3bm-7 | minor | RP's presentation (unforced rows on pinned windows, an exact zero, a NaN) | the continuation rows leave RP (RC replaces them); RP's float64 sides are labelled per window, "None" where the pair stays merged |
| D3bm-8 | minor | literature; New overstated | Status: the overlap singularities (Fukumizu–Amari 2000; Amari–Park–Ozeki 2006; Wei–Zhang–Cousseau–Ozeki–Amari 2008; Wei–Amari 2008), Saad–Solla's plateau (1995), duplicated units' critical sets (Şimşek et al. 2021; Zhang et al. 2021), balancedness (Du–Hu–Lee 2018) [memory]; New: the relative mode's law in this setting, and the window's pin and climb a caveat of the sampled setup, not a law |
| D3bm-9 | cosmetic | slips | "within 2.3 steps"; "15–216 e-folds"; C4's artefact quote dropped |

**Lesson.** A continuation below the floating-point floor is only as good as the state it is linearised about: continued over the run's
own trajectory after that trajectory had left the merged state, a law of the right form gave the right sides and wrong numbers. Linearising
the exact Jacobian about a merged base and injecting the result back into the full race gives the offset, the side, the climb and s_d with
a validation built in — and showed that the hold begins only when the climb ends, so a departure time inherits the pin one for one. A
tolerance "below the masses" is not enough when an outcome is decided at 10⁻¹⁴σ: the absolute tolerance sets the step on every
component, the dormant ones included.

## Round 71 (adversarial review of Proposition D‴.52 as revised after round 67, 2026-10-28)

**Design.** One referee, not modifying the repository. D3br: Proposition D‴.52 as revised after round 67, with its code and results at
42aa6a3 (`review/REVIEW_round71_D3br.md`). The referee:
- regenerated the md byte for byte, compared the raw with round 67's and reran all 15 U jobs, all 95 S jobs and two F jobs bit-identically;
- re-implemented the rejections, the crossing and the pairing independently and tried the knot's allocation per offset, per cell, from
  the law-calibrated knot and at the realised per-wall rates, with paired bootstraps;
- drew U's law from E's own samples and from ten further shared seeds;
- checked (iv)'s cosine in the continuum and Smith's constants by finite differences;
- ran 70 single-point mutants, round 67's survivors among them.

**Verdicts (D3br).** Reproducibility, the like-for-like comparison (robust to every allocation tried), F's new lines, the uncensored margin,
the minimum's loss off the tie, (iv) and Smith's constants hold. One major: the minimum's "small gain at the tie" was law-sampling noise — the
full knot and the runner-up's own knot were calibrated by two independently sampled 20000-sample laws, whose α quantiles differ by up to 3%;
drawn from the same samples they coincide at the tie (the pair sets the whole tail; the full knot is never later than the pair's own) and the
part is exactly 0. Minors: the tests credited for Z02, D02 and the uncensored margin never ran the mutated code (Z02, D02, SM01, SM05
survive and would change E′); the one-sided half compared a conservative knot at the nominal level (sizes 0.0446 against 0.0478; at equal
size the gap is ≈3%); the Liu–Nehorai/Smith sentence; loose ends.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3br-1 | major | The minimum's gain at the tie a Monte Carlo artefact | `u_job` draws the own knot's law from E's law samples (rng [seed, 9]) and records the full knot's quantile from them (q_full, equal to E's q_law in all 15 cells) and the null's margin; U rerun: the part −0.1225 … −0.0250 at c̄ ≤ 0.032 and exactly 0 at c̄ ≥ 0.045 (6 of 6); E′'s summary prints the shared samples' check; (iii), New and the Numerics: "none at the tie" |
| D3br-2 | minor | Tests not reaching the mutated code | Tests: `draw_stats`' t against z·s/(‖r‖/√(N − 2)) and z_σ's sign against sign(ŵ)·z for a pair of each sign (Z02, D02); `s_job`'s mu_u against 1 − wall·g_{j₂}/\|g_{j₁}\| and mu_u < 0 exactly on the swaps at the swapping tie cell (SM01, SM05); `u_job` against the raw and q_full = E's q_law in every cell; `f_alloc_lines` on synthetic rows (FA01–FA09); `smith_line` with a closed-form bound (SL01–SL05); `equal_size_line` on synthetic rows; the own-knot test on E's samples |
| D3br-3 | minor | The one-sided half at unequal sizes | `equal_size_line`: the runner-up's wall (the fit's sign for a same-sign pair, the other for an opposite-sign one) at the level α′ whose pooled H0 size is nearest the one-sided z's — (0.8, 0.6): α′ = 0.055, size 0.0476 against 0.0478, paired 1.032 [0.973, 1.239], ahead in 4 of 24; (0.8, −0.4): α′ = 0.051, 1.652, ahead in 0 of 24; (iii) and New: ≈3% in Δ₅₀ (≈13% in samples) at equal size, 4–6% at the nominal level |
| D3br-4 | cosmetic | Liu–Nehorai and Smith | Status: a limit ∝ √CRB, like Smith's criterion, has the detection order Δ² ∝ s only with equal or known amplitudes; with unknown unequal amplitudes both sit at s ∝ Δ³, estimation thresholds |
| D3br-5 | cosmetic | Loose ends | sign(ŵ)z in the allocation's thresholds; the calibrated knot nominal at c̄ = 0 and 0.02, below nominal from 0.032; z's range over the draws' own fits (0.0021–0.0825); the two figures for the minimum's part labelled (the oracle/calibrated difference against the same-law part); E′'s margin column the null's (0.73 at the swapping tie cell); the proof's intermediate remainder O(ϑ²(Δ/σ)² + (Δ/σ)⁴); the Status's scope (the decomposition one pair over half a cell). Round 64's entry (D3bk-11) should read "equal to round 59's but for round-off: the rms by ≤ 6·10⁻¹⁶, the end coefficients by ≤ 6.5·10⁻¹² relative" |

**Lesson.** A difference between two tests, each calibrated by its own Monte Carlo law, is only as accurate as the two calibrations agree:
two independent 20000-sample laws that share their whole lower tail in distribution still differ in their sampled quantiles, and that
difference produced the "gain" in 5 of 6 cells. Drawing both laws from the same samples is a one-seed change that makes the part exactly
zero, as the inequality bounding it (the full knot never later than the pair's own) requires [**erratum, round 75 (D3bv-7):** the inequality makes the part ≥ 0 at a common threshold — the share of H1 draws whose two knots straddle it —, not zero; it is zero here because no draw's knots straddle the common threshold, which round 73 counted]. And a test written to kill a mutant must
execute the mutated function: checking a report on hand-made rows cannot see a broken score, width statistic or margin upstream.

## Round 68 (adversarial review of Proposition D‴.59 as revised after round 61, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bo: Proposition D‴.59 as revised after round 61, with its code and results at
2391685 (`review/REVIEW_round68_D3bo.md`). The referee:
- regenerated the md, reran G, H, H2, H3 and HL jobs and A–D;
- ran the loop on H's asymmetric targets to τ = 10⁵, the path without its polish, and the landing on a fine time grid down to ε = 10⁻¹²;
- predicted τ₃ from the path's duals, measured E against Δc at frozen v and ran target 20 to τ = 10⁷ and target 9 to 5·10⁹;
- traced H2's ends (stop reasons, a re-entering kernel of negative width, invisible spikes, capped polishes, a converged re-polish);
- ran 73 single-point mutants.

**Verdicts (D3bo).** (i), (ii)(a), the dipole law and v ∝ ε^{−1/6} hold. Three majors: "the cancellation is the path's", "which the loop
does not keep" were contradicted — on asymmetric pairs the loop and the polish-free path keep it, and on the symmetric pair it is the
polish's (accepted at its evaluation budget); the pair does not land at τ₃ — an intermediate phase (three formed kernels, the pair on no
feature) lasts to the next formation τ₄, the landed v differs from the law's by an O(1) factor κ, the grid's τ₃ was late, and "landed" is a
threshold p ≲ 0.28; H2's "exploding times while the residual creeps down" were artefacts (the dictionary exhausted, a negative-width
re-entry at τ = 2.3·10⁶⁷, invisible spikes, capped polishes). Minors: the collapse law (E ∝ Δc⁸: τ^{−1/6}, τ^{−2/3}; "frozen" the same law
later); 13 non-equivalent surviving mutants and untested job code; (i)'s "from any start"; literature; C's table; cosmetic.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bo-1 | major | The cancellation's attribution | HL: the asymmetric targets at ε = 10⁻⁶ and 10⁻⁸ to τ = 10⁵ (11–12 negative kernels and \|w\| 2.19–2.69, still growing, at 10⁻⁸; 4–6 at 10⁻⁶; c₂ +0.01 exact there); HP: the path without its polish (the equal pair \|w\| 0.929 against 3.24, the asymmetric ones keep the cancellation); `stop_anatomy` flags the polishes stopped at their budget; (ii)(b), Status, New: the path's and the loop's on asymmetric pairs, the polish's on the symmetric one |
| D3bo-2 | major | The landing | G on a fine grid (τ = 10 … 140 by 0.1) at ε = 10⁻⁵ … 10⁻¹² (`fine_phases`): τ₃, τ₄, the dipole phase's slopes, v(τ₃⁻) over the law, the intermediate phase, the landing's τ, κ and p; GT: τ₃ from the path's duals (`t3_from_duals`: the duals advanced to τ_d — the path's history stops at the previous entry — and moved by the block's residual at s = 10¹¹), within 0.073 of the loop's at 10⁻¹²; (ii)(c) and the proof: four stages, v_land = κv(τ₃⁻) (κ 0.752–1.116, ε-free to 0.4%), landed ⇔ p ≲ 0.28 from (a)'s moments (34067 of 36535 states agree) |
| D3bo-3 | major | H2's ends artefacts | H2 prints why each path ended (active kernels of 21, polishes stopped at their budget, negative widths, kernels invisible on the samples); (ii)(b), Status and Numerics: the ends the dictionary's and the polish's, "exploding times" dropped |
| D3bo-4 | minor | The collapse law | (ii)(c) and the proof: Δc⁻⁶ = Δc_L⁻⁶ + 48aGm(W²/(w_a² + w_b²))Λ(ε)(τ − τ_L) from E ≈ aΔc⁸ ((a)'s O(Δc⁴)) and the metric restricted to a pinned centre of mass; GL: target 20 at 10⁻¹⁰ to τ = 10⁷, the exponents −0.1669, −0.3320, −0.6667; "frozen" dropped (landed at a small p, the same law later) |
| D3bo-5 | minor | Mutants; untested job code | tests: `pair_fits` and `pair_index` at their boundaries; `h_stats`' 1% threshold; `stop_anatomy`; `fine_phases` with an intermediate phase; `t3_from_duals` against GT and the loop; `hl_job` against `loop_nd`; `h_job` with a tolerance (the last entry's τ, the last blocks) and `h3_targets`; the oracle against a direct fit |
| D3bo-6 | minor | (i)'s start | (i), the proof and the Status: from a start of energy below ½‖f*‖² (a large coefficient drives the width down: ∂_σE = ½π^{1/2}w(w − A)) |
| D3bo-7 | minor | Literature | Status: Xu–Du (COLT 2023) and Dwivedi et al. (Ann. Statist. 2020) [memory]; the unverified NeurIPS 2025 entry dropped; New narrowed to the stages, τ₃ from the path and the derived collapse law |
| D3bo-8 | minor | C's table | C's table counts a budget-stopped path at its stalled block (`moving_state`): formed kernels 2–4, the largest residual 0.12 and 0.14 |
| D3bo-9 | cosmetic | Wording | the null cluster no longer "idle"; "the five cells with a phase" replaced by the fine grid's six |

**Lesson.** A conclusion of the form "X is the path's", "the pair lands" or "the cell is frozen" needs the control that could make it
false: the loop on the asymmetric pairs and the path without its polish, the landing timed on a grid finer than the event, the collapse
followed past the window that made it look frozen. Each control here overturned the conclusion it tested and strengthened the result
around it: the dipole phase obeys the block's law, τ₃ follows from the path alone, and the collapse has a law.

## Round 73 (adversarial review of Proposition D‴.52 as revised after round 71, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bt: Proposition D‴.52 as revised after round 71, with its code and results at
1f51d48 (`review/REVIEW_round73_D3bt.md`). The referee:
- regenerated the md byte for byte, compared the raw with round 71's and reran two U jobs, one E S job and two F jobs bit-identically;
- re-implemented the equal-size comparison and tried four other ways of matching the size, with a paired bootstrap that re-selects α′;
- computed the minimum's part in the population (400 000 law samples and 40 000 fresh H1 draws per cell) and its Monte Carlo errors;
- checked the inequality between the full knot and the pair's own draw by draw, the null's margin against the swap rates, and (iv)'s
  remainder in the continuum;
- ran 69 single-point mutants, round 71's among them.

**Verdicts (D3bt).** No major issue. Reproducibility, the zero at the tie (also in the population), the inequality (proved; the paths agree
to round-off), the equal-size comparison (robust to every matching tried), the null's margin, (iv)'s remainder and the loss's first two
parts hold. Minors: the zero was argued from an inequality that only bounds the part from below — the minimum cannot lose; that it gains
nothing is the empirical absence of H1 draws between the two knots, which no script counted —; the part off the tie quoted to four decimals
without its Monte Carlo error; tests that run the mutated code on data where the mutant changes nothing; the size wording reading one noisy
offset as a pattern (c̄ = 0.08 is 0.02 from the next kernel, the mirror of 0.02); loose ends.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bt-1 | minor | The tie's zero argued from a one-sided bound | E′ prints the full knot against the own knot draw by draw (the path's exceeds it in 3293 of 21000 draws, by at most 6.5·10⁻¹² relative: round-off) and the H1 draws between the two knots at the tie's common threshold (0 of 2400); (iii) and the Numerics: the minimum cannot lose at a common threshold (its part is the share of draws between the knots) and gains nothing because none falls between them |
| D3bt-2 | minor | The part off the tie without its error | E′ prints each cell's standard error over the paired H1 draws (0.013–0.021; the laws' own Monte Carlo error not included); (iii): "other kernels' knots", the loss 0.025–0.12 with those errors |
| D3bt-3 | minor | Tests on data the mutants do not reach | Tests: `s_job`'s uncensored margin for a pair of each sign; `u_job`'s null margin for a negative pair and the null's own first entrant equal to f₁'s on every E cell; E′ on the raw (each cell's part and standard error, the inequality's excess, the draws between the knots, the margin '—' without U); `f_alloc_lines` pooling mirrored offsets; `equal_size_line` with α′ below α, a cell that never crosses ½ and a cell with a separation missing |
| D3bt-4 | minor | The size wording at a mirrored offset | F prints the calibrated knot's size by the distance to the nearest kernel, an offset and its mirror pooled: 0.0512, 0.0456, 0.0398, 0.0379, 0.0329 at 0, 0.02, 0.032, 0.045, 0.049; (iii) and the Numerics restated by distance |
| D3bt-5 | cosmetic | Loose ends | the Status's repeated "α allocation"; round 67's −0.0075 … +0.0175 replaced by what the script prints; `u_job`'s null margin over the null's own first entrant argmax\|Φᵀf₀\| (the same kernel in all 15 cells, so U is unchanged — asserted by a test); E′'s margin '—' without U; the score as the limit of the score process of the one-against-two likelihood ratio test (Dacunha-Castelle–Gassiat 1999; Azaïs–Gassiat–Mercadier 2006), Amar–Weiss 2008, El Korso et al. Round 71's entry (D3br-1 and its Lesson) should read: the full knot never later than the pair's own bounds the part from below; the zero is the absence of draws between the two knots |

**Lesson.** An inequality settles one side. "Never later" with equal quantiles shows the minimum cannot lose at the tie; that it gains
nothing is a separate fact about the draws, and a script should count it rather than print it as a consequence. A Monte Carlo difference
needs its standard error before its extremes are read, and on a grid with mirror symmetry an offset and its mirror are one data point
seen twice: pooled, the size falls smoothly with the distance from the nearest kernel instead of jumping between two offsets [**erratum, round 75 (D3bv-1):** the fall is the weak null's ((0.8, −0.4), |ŵ| ≈ 0.4); the strong null (|ŵ| ≈ 1.4) holds its level until next to the tie, and the weak null's shifted centres are binned by their own distance] [**erratum, round 78 (D3by-1):** both nulls fall toward the tie — the strong null, E's and F's 24 600 draws pooled, 0.0505 → 0.0360 (trend z = −3.28), the weak 0.0508 → 0.0296–0.0325 (z = −3.24); their slopes do not differ (likelihood-ratio p = 0.16); F's strong null alone (z = −1.45) lacked the power to show it].

## Round 69 (adversarial review of Proposition D‴.60 as revised after round 63, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bp: Proposition D‴.60 as revised after round 63, with its code and results at
e65d6a0 (`review/REVIEW_round69_D3bp.md`). The referee:
- regenerated the md byte for byte and reran cells of A2, F, B3, B5, C, C2 and C3 bit-identically;
- re-derived (i) and (ii), linearised the C loop by central differences (the full monodromy) and split the period means into mode sums
  and the period's own response;
- ran each C3 draw's own second order and its −n twin, single tones to s = 0.16, the noise-free loop from a raised coefficient against m,
  k_e and the metric, pre-raised band and tone runs, and the band without single pairs;
- ran 65 single-point mutants, 12 of them replicas of round 63's survivors.

**Verdicts (D3bp).** Reproducibility, (i)'s law and √K_i, F's monodromy (independently), (ii)'s Floquet algebra, the averaged loop as the
constant part's response at O(s²), C3's chronology and the jump-level trends hold. Two majors: (iii)'s "by about (ii)'s even part at its
saturation" compared ten other draws' leaves with one draw's saturated second order — like for like the rise is the loop's even part at
the level, about a quarter above the draw's own second order and before saturation —; and (iii) had no criterion for the jump and could not
explain single tones, which drift further and never jump: the missing piece is the fit's basin edge in the coefficient (ejection from about
+0.075 at m = 12, k_e = 10), which the noise moves. Minors: the Floquet statement described the mode coordinate (H missing) and the
stability holds for the formed block; the "resonance factor" is the loop's sensitivity (> 1 above √(K/2) whatever k_e); the peak shift is
the quasi-static law's at m = 12; B5 at s = 0.005 mixed in the averaged loop's O(s⁴); the increments are not monotone in ω; pooled quotes
and an untested first-passage dichotomy; B3's exceptions; 18 behaviour-changing surviving mutants; sources; cosmetics.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bp-1 | major | Not like for like | C4: the −n twins of C3's band runs (the half-sum over the rise 0.981–1.019 at s = 0.02, 0.972–1.030 at 0.03); W2: each draw's own second order over 300 periods (±0.00125·n) per configuration — the rise 1.024–1.593 (median 1.254) times s²⟨w₂⟩ at the leave period, which is 0.64–1.01 of its period-300 value; (iii): the rise is the loop's even part at the level, above the second order, before saturation |
| D3bp-2 | major | No criterion; single tones unexplained | EDGE: the noise-free edge against m, k_e and the metric (0.0731 at m = 12, k_e = 10; exponents 0.537 and 0.555; Δ_e(v₁/k_e)^{1/2} 0.0756–0.1017 over all rows, v₁ the centre mobility; WFR's edges 3.6–4.4 times larger), the ejection's anatomy (the orbit point repels an over-amplitude kernel and carries it out within one downward passage), s* = (Δ_e/⟨w₂⟩)^{1/2} against the first jumps with every case classified (33 of 52 consistent, 11 later, 8 earlier, 17 without a jump); C6: pre-raised band and tone runs; C7: the edge under noise — single tones raise it by E·s² (E ≈ 350–1540 by tone) while their drift saturates, the band moves it window by window and its least window edge decides the jump (3 of 3 where the draw jumped, 10 of 11 where it did not; the exception WFR at s = 0.24); the title, (iii), the Status and New |
| D3bp-3 | minor | The Floquet statement | (ii) and the proof: ⟨z⟩_j = Σ_v(Pv)a_j + Re(e^{iΩjT}H) and the bound with \|H\|; the formed kernels' linearisation stable, the unformed kernels neutral at O(ε); the averaged loop's error "at most" the sum of moduli |
| D3bp-4 | minor | The sensitivity | (i), the proof and A2: per tone at a frozen K the factor is \|S\| = \|iω(iω + k_e)\|/\|K − ω² + ik_eω\| (the plant's slope cancels), \|S\|² − 1 = K(2ω² − K)/\|K − ω² + ik_eω\|², \|S(√K)\| = (1 + K/k_e²)^{1/2}, Bode's integral (0 by quadrature to 3·10⁻¹²); at the bound \|S\| = 1.41, 1.38, 1.24, 1.11 at ω = 10, 15, 20, 30; the band's factor per kernel 1.495, 1.248, 1.370 at m = 12 and 1.225, 1.095, 1.152 at the bound |
| D3bp-5 | minor | The peak shift | (i): the quasi-static response at m = 12, the sweep broadened by k_e; A2 marks maxima at the grid's lower end |
| D3bp-6 | minor | B5 at s = 0.005 | B5's averaged loop also at s = 0.00125 (over the constant part 0.998–1.005 at 13.05, 13.1, 13.15) and 29.3 added (2ω = 293ω_ref, its increment small: 1.0000–1.0028 of the exact second order); "located" restricted to 13.1, the 16 tones' error not located |
| D3bp-7 | minor | Increments and frequencies | (ii): "falling with ω on average, with interference fringes"; the products' frequencies modulo nω_ref, the first order's transient included |
| D3bp-8 | minor | Pooled quotes; the first passage | C3 prints its summary per kernel (kernel 1's offset before −0.016 … −0.011, its jumps ending beyond ±1 in 11 of 11 at s = 0.02; kernel 3's ending in the vacated first feature in 3 cases); the dichotomy dropped, the leaves placed against the edge instead; C2's caption "out of its feature" |
| D3bp-9 | minor | B3's exceptions | Numerics: 13.15's and 13.05's first kernel (not polynomial by s = 0.005 near √K₃ = 13.13), 21.37's third and WFR 29.3's first (noise level) |
| D3bp-10 | minor | Surviving mutants | Tests (47): the Floquet means at a rotating tone (21.37), `f_job`'s sums and saturations, `b5_job`'s phases and both levels, `c3_job` with a fake loop and a raised start, `leave_stats`' window, `c_noise`'s negation and shift, `edge_run`, `edge_job`, `bracket_edge`, `c7_job` and `c7_lines` (other configurations too), `w2_job`, the C4, EDGE and C6 reports, the sensitivity and Bode's integral, the edge's scaling and `star_verdict`, C2's single-tone rises on the raw |
| D3bp-11 | minor | Sources | Status: Bode and Freudenberg–Looze; Egardt, Ioannou–Kokotović; bursting (Anderson 1985) for the phenomenology, kept out of New; Kokotović–Riedle–Praly as slow adaptation's stability criterion; Hoppensteadt and Khalil for a uniform-in-time expansion on the formed block (not attempted) |
| D3bp-12 | cosmetic | Wordings | the factors per kernel; C2's caption; round 63's entry credited tests to "c3_job's statistics" — they covered `leave_stats` and `c3_lines` only (`c3_job` is tested from this round) |

**Lesson.** A mechanism written as a sequence — a drift, then a jump — is tested where the first step occurs without the second: single
tones drift further than bands and never jump. The missing piece was a threshold, the fit's basin edge, and measured with and without
noise it explains both: the noise moves the edge, a tone by E·s² while its own drift saturates, a band window by window, so that the band's
least window decides when its drift suffices. A comparison of theory with events must also fix the draw, the time and the level: like for
like, the rise at the leave is the even part at the level, a quarter above the second order, before saturation.

## Round 70 (adversarial review of Proposition D‴.58 as revised after round 65, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bq: Proposition D‴.58 as revised after round 65, with its code and results at
409ec98 (`review/REVIEW_round70_D3bq.md`). The referee:
- regenerated the md byte for byte, ran the 66 tests, and reran cells of B3, A5, B2 and three paths to their first non-optimal block
  bit-identically;
- re-derived the Mehler-frame paragraph (ℓ², v_ρ, the first tail ratio, the Mercer frame);
- reran paths with every polish start recorded and restarted the polish from those starts;
- solved the two below-floor blocks by Levenberg–Marquardt with the residual in double and in long double;
- compared the slow modes in three unit conventions and recounted the entries;
- ran 91 single-point mutants plus a control.

**Verdicts (D3bq).** Reproduction, every printed number, (i), (ii), (iv), (v), the proof, the Mehler-frame algebra, the seed rule's
agreement with `sliding_iss` and the entries hold, and the floor is the arithmetic's (long double passes it). Two majors: the round's new
code was weakly tested (45 of 91 mutants survived, 25 of them changing what a rerun publishes); and B3 restarted from the block's final
state, not from the polish's start, with two settings changed together and the ends unjudged — at two ends its own "path" column showed the
path's settings succeeding. Minors: nothing in the repository showed the floor, the ¾ reference was 0.16% high and the floor soft; the
slow-mode comparison depended on the units; events merged entries after the paths parted; the Mercer frame's rule is slower, not "nearly as
fast"; cosmetics (the repolish's improvement test; numbers no script prints).

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bq-1 | major | Weak tests of the new code | Tests (35): `b3_job` with a variant that has a gradient test (the four restarts' budgets and tests in order), `path_polishes` (records, stops at a call, restores), `first_nonoptimal` and `b3_jobs`, `b3_binds` and `b3_lines`, `seed_rule` with two seed kernels (one decided by its ascent target; its 10⁻² pinned to `sliding_iss`'s), `part_A5` on its own numbers (the Mercer frame's rule, not the Mercer-zero one), `b2_job` where variable projection leaves the block and `vp_fit`'s order, `run_B` with variants merging, `lm_fit` in double and long double, `b4_jobs`, `b4_job` and `b4_lines` |
| D3bq-2 | major | B3 not what the path did | B3 rebuilt: the path rerun with every polish start recorded (stopped at its first non-optimal block's first polish, located from B's stored polish counts), that polish restarted from its start one setting at a time (the path's; 5·10⁵ evaluations; no gradient test; both), each end judged by the path's acceptance (`seed_rule`); above 10⁻¹⁰ the budget binds at 7 ends, the budget and the radius at 5, the budget and the gradient test at 5, the test at 2, the radius at 1; with both changed the solve reaches the optimum from all 20 (0.99998–1.00022 times it); (iii) and the Numerics end by end |
| D3bq-3 | minor | The floor | B4: Levenberg–Marquardt (`lm_fit`) with the parameters and the residual in double and in long double from the blocks below 2·10⁻⁹, from the path's block and from the stored optimum: in long double both starts end together (within 8.2·10⁻⁹), in double from the blocks 1.9·10⁻⁸ … 2.3·10⁻³ above that end (0.2% below 10⁻¹⁰) — the floor the arithmetic's and soft, crossing 0.1% near 10⁻¹⁰; the stored optimum 1.6·10⁻³ above the long-double end at ¾, K = 11; B3 prints the restarts' excess by optimum; "whatever the budget" dropped (a single solve; at ¾, K = 11 the 5·10⁵ budget is spent) |
| D3bq-4 | minor | The slow modes' units | B2 prints the condition numbers in the polish's own units (1.3·10⁶–5·10⁸) and the ratio of the smallest singular values in three conventions (median 1.01 with the centres in units of σ, 1.54 in the polish's, 570 with unit columns — 49–2790); (iii) states the dependence and drops "not in their slow modes" |
| D3bq-5 | minor | Events after the paths part | (iii) and the md's caption: an event is the same entrant at the same iteration, whether or not the variants' paths have parted |
| D3bq-6 | minor | The Mercer frame's rule | A5 prints its ratio per kernel at the last three K — level at 5.32 (τ/σ = 1) and 2.33 (τ/σ = 2), below 1/t = 5.83, 2.62 —; the Status: it falls more slowly (its first term at qv′/2 in its own frame, the tail setting the rate), and "a frame's orthogonality measure is its own Gaussian only at ρ = λ" |
| D3bq-7 | cosmetic | The repolish's improvement test | B3 restarts the first polish only, which has no improvement test; the "path" column is the path's settings from the recorded start |
| D3bq-8 | cosmetic | Numbers no script prints | Status: the j ≤ 10, K ≤ 40 check, the K = 640, 1280 gaps and the 80-digit check marked as referee rounds' computations that the repository's scripts do not reproduce |

**Lesson.** A restart tests the run it restarts only from that run's own state: from the block's final state the "path" column was a second
polish, and changing two settings at once could not say which binds; recorded starts, one setting at a time, and the path's own acceptance
answer it end by end. A floor claimed for the arithmetic is shown by changing the arithmetic: the same solve with its residual in long
double passes it from the block and from the optimum alike, and the gap by optimum shows the floor soft rather than a threshold. A
comparison of conditioning is a statement about units unless it is scale-free.

## Round 75 (adversarial review of Proposition D‴.52 as revised after round 73, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bv: Proposition D‴.52 as revised after round 73, with its code and results at
3dc9443 (`review/REVIEW_round75_D3bv.md`). The referee:
- regenerated the md byte for byte, compared the raw with round 71's and reran all 15 U jobs bit-identically;
- recomputed the minimum's part in the population (10⁶ shared law samples and 10⁵ fresh H1 draws per cell, the vectorised fits checked
  against E's 21000 draws);
- re-derived the part's decomposition at a common threshold and off the tie, and checked the excess statistic draw by draw;
- split the size by distance by null and by the null's own centre;
- ran 44 single-point mutants of round 73's code and round 73's surviving ones again.

**Verdicts (D3bv).** No major issue. Reproducibility, the numbers, the inequality and the part's decomposition, the excess statistic
(round-off, balanced in sign), the paired standard errors, the tie's zero (also in the population: 0 of 600 000 draws between the knots)
and the size by distance as computed hold. Minors: the size's fall with the distance is the weak null's (|ŵ| ≈ 0.4, trend z = −3.83); the
strong null (|ŵ| ≈ 1.4) holds its level, and the weak null's centres at s = 0.03 sit 0.0097 off their label c̄, so pooling by c̄ did not
pool mirrors; (iii) quoted the extremes of nine noisy estimates as the off-tie loss (0.034–0.087 in the population); three sources cited
for results they do not contain (Azaïs–Gassiat–Mercadier for N^{−1/4}; Dacunha-Castelle–Gassiat as if for the regression model; Amar–Weiss
tagged [snippets]); older reporting code caught only by the md diff. Cosmetic: zero standard errors printed at the tie; infinite knots
dropped by a filter; wording and log loose ends.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bv-1 | minor | The size by distance pools two nulls | F bins each H0 null by its own centre's distance to the nearest kernel (`null_distance`) and prints the line per pair (`size_by_distance`): the strong null (0.8, 0.6) 0.0517, 0.0471, 0.0421, 0.0450, 0.0383 and the weak (0.8, −0.4) 0.0508, 0.0429, 0.0379, 0.0296, 0.0325 at 0, 0.02, 0.032, 0.045, 0.049; (iii) and the Numerics: the knot holds its level, the fall with the distance is the weak null's; an erratum in place at round 73's Lesson |
| D3bv-2 | minor | The off-tie loss quoted by its extremes | E′ pools the part per offset over its three cells with its standard error (−0.104 ± 0.012, −0.085 ± 0.011, −0.029 ± 0.008 at c̄ = 0, 0.02, 0.032); (iii) quotes these |
| D3bv-3 | minor | Sources | Status: Dacunha-Castelle–Gassiat 1997, 1999 as the mixture analogue (population mixtures, positive weights) and the regression analogue in additive noise (Davies 1977, 1987; Liu–Shao 2003; Gassiat 2002 [memory]); the n^{−1/4} order from Chen 1995 and Chen–Chen–Kalbfleisch 2001 [memory]; Azaïs–Gassiat–Mercadier 2006 for the cost of searching a large parameter set, the analogue of the dictionary's minimum; Amar–Weiss [memory]; El Korso et al. beside Liu–Nehorai, with their HAL records |
| D3bv-4 | minor | Older reporting code untested | Tests: `f_alloc_lines` per pair and offset (H1 rows left out, qA ≠ qB, a `q_law`, a ŵ < 0 row, the α = 0.01 figures); `smith_line` with a non-integer power; `f_lines`' range and equal-size lines; E′'s margin with mu0 ≠ mu_true; `size_by_distance` with H1 rows and a p exactly at α |
| D3bv-5 | cosmetic | Zero standard errors at the tie | E′ prints the discordant draws and their 95% bound instead (0 of 400: at most 0.0075 per cell, 0.0012 over the 2400 draws) |
| D3bv-6 | cosmetic | Infinite knots dropped | E′ counts an infinite path knot beside a finite own knot as an exceedance and prints the draws with a finite own knot (21000 of 21000) |
| D3bv-7 | cosmetic | Wording and log loose ends | "holds its level"; the threshold falls between the two knots of no H1 draw; "3293 of 21000" over all cells; E′'s header normalises by the null's own first entrant, and U prints that f₁'s and f₀'s first entrants coincide (15 of 15 cells); the Status's nested dashes; an erratum in place at round 71's Lesson |

**Lesson.** A trend in pooled cells can belong to one of its parts. The size's fall with the distance is the weak null's, and the strong
null holds its level until the tie [**erratum, round 78 (D3by-1):** not so — E's 15 000 draws of the strong null fall too (z = −2.97), and
the two nulls' slopes do not differ (p = 0.16); a failure to detect a trend in F's 9600 is not its absence]: before reading a pattern, split by what sets the statistic's scale (here s/|ŵ|) and bin by the
quantity the pattern is about (the null's own centre, not its label). And the range of nine noisy estimates describes their noise: a
pooled value with its error is the statement.

## Round 72 (adversarial review of Proposition D‴.54 as revised after round 66, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bs: Proposition D‴.54 as revised after round 66, with its code and results at
37918e1 (`review/REVIEW_round72_D3bs.md`). The referee:
- regenerated the md byte for byte and reran RC, RP, C4 and C4M jobs bit-identically;
- wrote the ε-free race from scratch with closed-form integrals (no samples) in 130-digit arithmetic, fully nonlinear, and reproduced RC's
  least |D| (10^{−98.716} at x = 3);
- computed the block's rest's drift in float64, the long double and exact sums on the committed and on an ideal grid, and added it to RC's
  linear system as a forcing;
- checked s_d's identity against the climb's threshold, the coefficient's track, the proof's y and the κ-law's next order;
- ran 62 single-point mutants; eleven literature searches.

**Verdicts (D3bs).** No major issue; RC is right about the continuum. Minors: on samples the far window's drift is not rounding but the
samples' own float64 coordinates (3.2·10⁻¹⁶σ per unit W, the same in the long double), deciding x ≤ 5.48, and on an ideal grid the
aliasing (10⁻³⁷σ) decides x = 3; s_d's shift is an identity but its split into a climb and a hold a convention, and the coefficient's
decay begins only at the rest, after a rise; the proof's y is the slaved centroid's first-order response, not a forcing, and the κ-law's
error is ε_b at next order; 16 behaviour-changing surviving mutants (RC's checks vacuous to the tests); literature (Guo et al. 2018,
Cousseau–Ozeki–Amari 2008, Biehl–Riegler–Wöhler 1996; Du–Hu–Lee's homogeneity) and a New that claimed the generic trailing. Slips.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bs-1 | minor | What moves the rest on samples | RD: ∂c_b*/∂W at RC's merged block in float64, the long double and exact (80-digit) sums on the committed samples ((+3.151·10⁻¹⁶, −8.227·10⁻¹⁷)σ per unit W; the long double's to 3.2·10⁻⁴), exact on the ideal grid ((+1.249·10⁻³⁷, −6.96·10⁻⁴⁰)σ, the aliasing's first harmonic to 2.4·10⁻¹⁴) and on a dyadic grid (exact coordinates: 5.6·10⁻⁴³σ, each arithmetic's rounding beside it); RC: the drift as a forcing κk(∂c_b*/∂W)dW/ds continued with the relative mode — the coordinates' forced D exceeds the continuum's least at x ≤ 5.48, the aliasing's at x = 3, both toward the far side on C2's target —; (iv), (v), the proof, Status, New, the RP/RC captions: the continuum in exact arithmetic, on samples the edge's, else the aliasing's or the coordinates' |
| D3bs-2 | minor | s_d's split | RC and RP: the climb's end at θ = 10⁻³, 10⁻², 10⁻¹σ — s_d − s_climb(θ) equal on Part C's window and in the continuum to 7.5·10⁻¹² … 2·10⁻⁴ for every θ, the hold moving by ln 10/λ per decade —; the runner-up's coefficient along the injected race (rising ∝ \|D\|², peaking as D reaches its rest: +9.6·10⁻⁴ at 0.21σ, x = 3; 1% below from 6.56·10⁵); (v): the shift an identity, the split a convention, the hold beginning at the rest |
| D3bs-3 | minor | The proof's y; the κ-law's error | RP: the (D, y) system on the block's moments, the centroid unslaved — within 5·10⁻⁴ of the Jacobian's rates where the κ-law errs by up to 10.1% — and ε_b = (w − w_e)²ww_e/(w² + w_e²)²·R⁻¹Q (0.4–8.8%); the proof: y first order in D, the κ-law dropping dy/ds, its error ε_b; in the continuum nothing forces D (y = 0 moves the least by 0.016 dex) |
| D3bs-4 | minor | Surviving mutants | Tests: golden `rc_job` rows for the mirror (x = 4.9), the tail (x = 7) and every check row, each check differing from the main row where its argument acts; the new fields in the golden comparison; `part_RC`'s job list; `rp_job` on the far window at x = 3 and 5.07 (the law, D, the growth rate, the climbs); the unslaved rates and ε_b on the raw; `exact_drift`, `alias_drift`, `sample_coords` and RD on the raw (55 pass) |
| D3bs-5 | minor | Literature and New | Status: Guo et al. (JMLR 2018), Cousseau–Ozeki–Amari (IEEE TNN 2008), Biehl–Riegler–Wöhler (J. Phys. A 1996); Fukumizu–Amari's line of critical points for the frozen split; Du–Hu–Lee restricted to positively homogeneous units; New: the race's consequences in this paradigm (the side, the pin, the climb, s_d), the relative mode's law generic to (a, θ) gradient flows |
| D3bs-6 | cosmetic | Slips | 0.02 dex (y = 0: 0.016); the formed range 27–216 e-folds and 10^{−13.4} … 10^{−98.7}σ (x = 7 a tail); 6.97–7.49·10⁴ (the mirror's to 7.50·10⁴); R's x = 5.4772 beside RP's 6.10·10⁴ at 5.48; the κ-law's ranges at s ≈ 450 and 5000; the top edge at 1.5 (the right edge 6.25σ away; float64's and the long double's drifts); the mirror's −0.64; x = 6.04's validation; `c4_atol`'s docstring; C4's pin at (4.9, 10⁻⁹) resolved to a factor 2 |

**Lesson.** "The continuum" must be said in full: the integral over ℝⁿ in exact arithmetic. A sampled problem has two deterministic
sources that move what the continuum keeps fixed — the samples' aliasing and their float64 coordinates — and they act in any arithmetic,
so a long-double run cannot tell them from the continuum, while exact sums on an ideal and on a dyadic grid separate all three. And a
law's next order is the cheapest check of its error: here the unslaved two-variable system reproduces the Jacobian's rate to 5·10⁻⁴,
where the slaved law's 10% had been read as a validity limit.

## Round 74 (adversarial review of Proposition D‴.59 as revised after round 68, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bu: Proposition D‴.59 as revised after round 68, with its code and results at
18bcfa0 (`review/REVIEW_round74_D3bu.md`). The referee:
- regenerated the md byte for byte and reran G's fine phases at 10⁻¹², HL, HP and GL bit-identically, and `t3_from_duals` on the six paths;
- fitted the intermediate phase's v³ and derived K′, τ₄ and κ from the path alone (a closure fit with the τ₃ entrant, K′ from the third
  order, the duals moved on);
- computed every formed kernel's largest φ on the samples in H, HP, H2 and HL, and ran the symmetric pair and A₂ 0.99 to τ = 10⁷;
- re-derived (i)'s ∂_σE, (a)'s moments, the dipole law, the collapse's metric with the slaved coefficients, and the duals' rate, and ran
  target 9 to τ = 10¹⁰;
- ran one kernel from (c*, σ*) with w = 1.5A … 10A, and 71 single-point mutants.

**Verdicts (D3bu).** Reproducibility, every number, τ₃ from the path, (a)'s landing algebra, K's third-order formula and the collapse
law's exponents hold. Two majors: the intermediate phase is a second dipole phase — v³ linear in τ with the three-kernel block's K′ — and
τ₄, κ and v_land follow from the path (κ not "ε-free": 0.729–2.905 over the cells, converging as ε → 0); (ii)(b)'s counts included dead
spikes invisible on the samples (4–6 of the loop's 11–12 negatives at 10⁻⁸), which were the "widths below the samples' spacing", and the
symmetric loop also keeps kernels of both signs. Minors: the collapse law's constant holds only as p → 0 (the slaved coefficients' drag,
10% at p = 0.15; a undefined; τ_L misplaced); the landing statistic for p from Δc alone (a σ-dipole); (i)'s "a large one is not enough"
(the condition is sufficient, not sharp); 13 behaviour-changing surviving mutants; the τ₃ proof's residual limit and the grid offset;
Ho–Nguyen and a narrower New; cosmetic.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bu-1 | major | The intermediate phase; κ | GT (`kappa_from_path`): K from the third order (within 1.4·10⁻⁴ of G's), τ₃ from the duals, K′/K from the three-kernel block's optimal dipole fit (+0.0461, −0.0833, −0.3563, 0, +0.2278, −0.0904), τ₄ from the duals moved on (124.3, 97.89, 81.64, 51.29, 46.97 against the loop's 124.2, 97.9, 81.7, 51.3, 47), κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3} (within 0.6% of the loop's at 10⁻¹²), v_land by the law, v frozen after τ₄ (2.7·10⁻⁵); G's fine table: the second phase's slope over 3μ₀K′Λ (0.9997–1.0012 at 10⁻¹²); (ii)(c)(2)–(3), the proof, Status, New: a second dipole phase, τ₄ and κ from the path, κ converging as ε → 0 |
| D3bu-2 | major | Dead spikes counted; the symmetric loop | `visibility`: H, HP, H2 and HL count only kernels whose largest φ on the samples is at least 10⁻¹⁰ (6–7 visible negatives at the ends, every formed visible kernel ≥ 0.049 wide on the paths, ≥ 0.036 in the loop); HP prints its stop reasons (dictionary, budget); HLL: the symmetric pair at 10⁻⁶ and A₂ 0.99 at 10⁻⁸ to τ = 10⁷ (8 visible negatives, \|w\| 0.705 → 0.839; 6 visible, \|w\| 2.38 → 2.55, undecided); (ii)(b) and the Numerics: "widths below the samples' spacing" dropped, the symmetric loop fits with both signs at \|w\| < A, the polish's only the inflated cluster |
| D3bu-3 | minor | The collapse law's constant | GL (`gl_constant`): the measured rate over the law, 0.903 → 0.999 as p falls 0.146 → 0.028, and over the law divided by 1 + C, C = GmW²(W − A)²(1/\|w_a\| + 1/\|w_b\|)/(μ₀Δc²(w_a² + w_b²)) the slaved coefficients' share of the metric, 1.018 → 1.004 — the remainder ∝ p, E's next order (a = E/Δc⁸ varying as 1 + 1.05p) —; target 9 added (0.991–0.998 at p = 0.015); (c)(4) and the proof: a defined, the law as p → 0 with C, τ_L after the relaxation, the induced-metric argument |
| D3bu-4 | minor | The landing statistic | `pair_p` uses \|Δθ\| (a σ-dipole has Δc = 0): 29177 of 29425 states agree (22800 of 22978 past τ₄); the threshold called derived |
| D3bu-5 | minor | (i)'s large coefficient | A2: from (c*, σ*) w = 2.5A, 3A, 5A converge exactly at σ*/σ₀ = 0.25, 1, 4 (the width dipping to 0.0058σ*), 10A ends as a spike; (i) and the Status: the condition sufficient, not sharp |
| D3bu-6 | minor | Surviving mutants | Tests: `fine_phases` with a σ-dipole, 0.28 < p < 1 off the features, the slope sampling and dW's sign; `part_GT` on a status-0 entry; `h_job`'s polish and tolerance reaching `sliding_iss` (monkeypatched) and clusters at unequal amplitudes; `gl_job` and `gl_constant` (and the remainder's fit); `visibility`, `path_stop`, a faint kernel in `stop_anatomy`; `d3phi`/`k_third`, `first_wall` and `kappa_from_path` against the raw; A2; `loop_row` and `hll_job` with a kernel between 0.1% and 1% |
| D3bu-7 | minor | The τ₃ proof; the grid offset | the proof: the block's state runs off but its residual converges, r = r_dip + O(s^{−2/3}), the duals' error O(ε^{1/3}); `fine_phases` evaluates the law at the time of v(τ₃⁻) (0.9995–1.0130 at 10⁻¹²) |
| D3bu-8 | minor | Literature and New | Ho–Nguyen (Ann. Statist. 44(6), 2016) beside Pearson for (a)'s O(Δc⁴); Xu–Du's potential on pairwise distances [snippet]; the exponents generic; New narrowed to the stage sequence from the path and the law's constant with C |
| D3bu-9 | cosmetic | Wording | "idle" dropped from the null cluster; "72 rows, 59 distinct"; G's stale t3_pred gone with G's rerun |

**Lesson.** A phase that a script labels "numerical" may be the same law again: the intermediate phase was the dipole law with the
three-kernel block's constant, and once it was seen as such every stage — τ₃, K′, τ₄, κ and the landing — came out of the path alone.
Counting kernels needs the criterion that the samples can see: spikes invisible on every sample had been read as a fit narrower than
the samples. And a law's constant checked against a run shows what the derivation leaves out — here the slaved coefficients' share of
the metric, and after it E's next order.

## Round 77 (adversarial review of Proposition D‴.58 as revised after round 70, 2026-10-28)

**Design.** One referee, not modifying the repository. D3bx: Proposition D‴.58 as revised after round 70, with its code and results at
39a8820 (`review/REVIEW_round77_D3bx.md`). The referee:
- regenerated the md byte for byte, reran all 12 B4 cells, B3 at two rows and B2 at three blocks bit-identically, and compared all 22 B3
  'path' restarts with B's stored first polishes;
- at every B4 start ran `lm_fit` in both precisions, evaluated the long-double end in both, measured the objective's rounding noise,
  scanned the damping at the double stall, restarted `lm_fit` from its end, and ran the same solve in pure double with its residual about
  an anchor (as is, re-anchored once, every 2000 steps); at the two blocks below 10⁻¹⁰ also compensated summation, identity scaling,
  μ₀ = 10⁻¹⁰ and the polish's TRF from the final blocks;
- tested `b3_binds` on nine synthetic cases, the Schur bound at the 20 first blocks under 22 centre scalings, and A5's first-term rates
  and the tail's onset;
- ran 99 single-point mutants and a control.

**Verdicts (D3bx).** Reproducibility, B3's restarts (22 of 22 the stored first polish), the extended-precision core and every number hold.
Three majors: the floor is not the arithmetic's — the same solve in pure double with its residual evaluated about an anchor reaches the
long-double end from all 12 starts, and at the double stalls the objective computed in double still resolves the gap (the rounding
destroys the step, not the objective); "censored by the floor", 'neither' and the soft floor's 0.1% crossing do not follow — at 3/4,
K = 11 the restart spent its budget, B4's double solve from that block passes the experiment's own 0.1% test, the polish's TRF from the
final blocks stops 4–8 times above the double solve's excess, and the excesses depend on the start and the solver; 11 surviving mutants
change B3 or B4 on a rerun. Minors: `b3_binds`' labels; the rerun guard and n0; B3 restarts the first polish only; B2's "whether … depends
on the units" against the Schur bound; A5's slower fall shown at two widths; cosmetic.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3bx-1 | major | The floor is not the arithmetic's | B4 rerun with `lm_fit(anchored=True)`: double, the residual r = r₀ − Σ_kφ₀_k[(w_k − w₀_k) + w_k·expm1(Δ_k)] about an anchor, re-anchored every 2000 accepted steps and at a stall, stopping flat, stall or cap, its end also after 1, 2 and 3·10⁵ iterations: from all 12 starts −1.4·10⁻⁸ … +7·10⁻¹² from the long-double end (the cap from 7, each end within 3.1·10⁻¹⁰ over its last 2·10⁵ iterations); the optimum-start rows in B4's summary (the double stall 2.4·10⁻⁸ … 2.3·10⁻⁵ above); B4's legend, `lm_fit`'s docstring, the test and (iii), Status and Numerics: the floor of the residual as these solvers evaluate it (f − Φ(c)w by direct subtraction), not the arithmetic's; "no step … as double resolves it" deleted; round 70's "the floor is the arithmetic's" withdrawn |
| D3bx-2 | major | Censored, 'neither', the soft floor | `b3_binds` carries the stop ("neither (the budget spent)", "neither (the solve's stop, status 3)"); B4's TRF column: the polish's own solve with 5·10⁵ evaluations and no gradient test stops on its step test 1.5·10⁻⁷ … 1.5·10⁻² above the long-double end from the final blocks (1.5·10⁻² and 5.3·10⁻³ at the two below 10⁻¹⁰); B3's excesses against B4's long-double end below 2·10⁻⁹ (1.0235 and 1.0484 at the two); (iii): the two ends not reached by the polish's solve (its step test at 1/4, its budget at 3/4), the double solve at 3/4, K = 11 1.0007 times the stored optimum, the stall's excess the solver's (non-monotone; at most 2.3·10⁻⁵ from the stored optima), "censored by the floor" and the 0.1% crossing dropped |
| D3bx-3 | major | Surviving mutants | Tests: `b3_restart` judges the seed from the start's amplitudes with the path's step (a 1.5σ move accepted at 2σ, rejected at 1σ); `b3_job` on the variant without the gradient test; `first_nonoptimal` with unequal polish counts and e_opt ≠ 1, and the guard on the stored first polish; `path_polishes` passes each call's settings to the solve; `b4_jobs` with two variants at one width (the row's variant, not the first) and an optimum in [10⁻⁹, 2·10⁻⁹); `part_B4` and `main --only B4` with a fake pool; a long `lm_fit` descent from the (1/4, K = 5) block (10.85 times the optimum to 1 + 5·10⁻⁹ in 7996 iterations); the anchored solve (cap, marks, the flat stop, each step's decrease when re-anchored after every step, long double refused); 24 of 24 targeted mutants of the round's code killed |
| D3bx-4 | minor | `b3_binds`' labels | labels from the ends that reach the optimum and are accepted; "(and the radius)" only where every reaching end is rejected; "nothing" an inconsistency; below 2·10⁻⁹ judged against B4's long-double end |
| D3bx-5 | minor | The rerun guard and n0 | `b3_job` refuses a 'path' restart whose status or evaluations differ from the stored first polish (`first_nonoptimal` carries it); B3: "the 'path' restart reproduces the first polish B stored — its stop and evaluations — at 22 of 22 blocks"; n0 counts every call, a block without kernels having no polish — the guard checks it |
| D3bx-6 | minor | The first polish only | (iii): the first polish's settings bind; at (1, grid +0.05) all 24 polishes stop on the gradient test and are rejected, their ends not logged; only the first polish is restarted — at four of the six radius-bound blocks the path made 22–24 rejected polishes of the block, the later ones from flowed states |
| D3bx-7 | minor | B2's units | `b2_line` and (iii): never faster in a common scaling of the centres (a Schur complement, for Kaufman's form; the Golub–Pereyra form used differs by a term ∝ the residual), how much slower unit-dependent; unit norm not a common scaling; the condition numbers with the centres in x, the polish's parameters, and TRF's Coleman–Li v^{1/2} |
| D3bx-8 | minor | A5's slower fall | `a5_lines`: the Mercer frame's tail overtakes its first term from K ≈ 2/(q(2 − v′)) = 70.0, 22.0, 10.2, 8.0, which the K tested pass only at τ/σ = 2; Status: the slower fall shown at τ/σ = 1, 2 only, at the tail's onset, its rate not shown |
| D3bx-9 | cosmetic | Wording | "its duals within 10⁻⁴"; A5's ratios .4f; variable projection's condition number printed once; B4's ratios quoted against the long-double end |

**Lesson.** The control for "the arithmetic's floor" is the best method in that arithmetic, not the same method in a wider one. An
anchored residual — its rounding fixed once, the change from the anchor computed to its own accuracy — lets a pure-double solve pass the
long-double end, and the floor becomes a property of how the solvers evaluate the residual, not of the problem or the arithmetic. A
verdict such as "neither" must carry the solve's stop, and excesses measured from one kind of start against an unconverged reference are
the solver's, not a law.

## Round 78 (adversarial review of Proposition D‴.52 as revised after round 75, 2026-10-29)

**Design.** One referee, not modifying the repository. D3by: Proposition D‴.52 as revised after round 75, with its code and results at
6764e41 (`review/REVIEW_round78_D3by.md`). The referee:
- regenerated the md byte for byte, compared the raw part by part with round 73's, and reran all 15 U jobs bit-identically;
- checked every number on the 32 lines round 75 changed, and all 598 decimal numbers of D52.md against the md;
- recomputed the sizes by distance per H0 job, with trend tests and logistic models (pair × distance), and E's calibrated sizes;
- ran the minimum's part in the population (10⁶ shared law samples and 10⁵ fresh H1 draws per cell), the N-scaling of Δ₅₀ (N = 61 … 481)
  and Smith's constants by an independent finite-difference bound;
- ran 63 single-point mutants (35 of round 75's code, 28 named by D3bv-4) and the identity.

**Verdicts (D3by).** No major issue. Reproducibility, every number, the tie's zero in the population, the pooled bound (by AM–GM), the
round-75 code, the N^{−1/4} law (slopes −0.247 … −0.250) and Smith's constants hold. Minors: "the size's fall with the distance is the
weak null's" is contradicted by E's strong null (15 000 draws, trend z = −2.97) and the two nulls' slopes do not differ (p = 0.15); the new
tests miss most of what they are credited for (the FA mutants, the bins by label, the pooled part and bound); three sources do not carry
what they are cited for (Chen 1995 an estimation rate corrected by Heinrich–Kahn 2018; Liu–Shao and Gassiat mixture results; the squared
positive part). Cosmetic: bins labelled by points; the pooled part at c̄ = 0 2.1 SE above the population; one-sided bounds; loose ends.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3by-1 | minor | The size's fall is not the weak null's | `size_by_distance` (with `size_nulls`, `trend_z`, `logit_fit`): F's and E's H0 nulls by their own distance — the strong null in F, in E and pooled (24 600 draws: 0.0505, 0.0502, 0.0459, 0.0433, 0.0360 [0.0305, 0.0420], trend z = −3.28), the weak (0.0508 … 0.0325, z = −3.24) —, and the two nulls in one logistic model (the slopes' difference z = +1.42, likelihood-ratio p = 0.157; the common slope z = −4.37); (iii) and the Numerics: "grows conservative toward the tie for both nulls"; the errata in place at rounds 73 and 75 (its Lesson) |
| D3by-2 | minor | The tests' strength | Tests: the per-offset size line and the walls test's components (qA ≠ qB; H1 at the two smallest separations only) asserted; `size_by_distance` with a null binned by its own centre into another bin than its label's, with E's nulls, the trend and the logistic line; `trend_z` and `logit_fit` against closed forms; `crb_sep` against a finite-difference bound; E′'s pooled parts, their contrast and the pooled bound on the raw; an infinite own knot left out of the finite count |
| D3by-3 | minor | Sources | Status: the score process's squared positive part; Liu–Shao and Gassiat with Dacunha-Castelle–Gassiat (mixtures), Davies alone as the regression analogue; the n^{−1/4} order as the contiguous alternatives' in testing homogeneity, Chen 1995 as the (corrected) estimation rate, Heinrich–Kahn 2018; Neyman–Scott's test and the score process in one clause; El Korso et al.'s year marked unchecked |
| D3by-4 | cosmetic | Bins labelled by points | `size_by_distance` prints each bin as its interval of distances with its nulls' own distances; (iii) "within 0.01 of a kernel", "at 0.041–0.050" |
| D3by-5 | cosmetic | The pooled part at c̄ = 0 | E′ prints the two lowest offsets' contrast (z = −1.19); (iii): the first two not distinguishable, the referee's population values (≈ 0.08, 0.08, 0.035) quoted as not reproduced; the laws' Monte Carlo error's size (0.003–0.006 per cell, the referee's estimate) |
| D3by-6 | cosmetic | One-sided bounds | `part_se` and E′: "≤ 0.0075 one-sided, 0.0092 two-sided at 95%" per cell and "0.0012 (one-sided 95%; 0.0015 two-sided) for the cells' mean share" |
| D3by-7 | cosmetic | Loose ends | Smith's fitted constants labelled by pair and amplitude case; the per-offset size marked as pooled over both nulls by the label; the first-entrant check said to be trivial below half the spacing; the dangling dash gone |

**Lesson.** A failure to detect a trend is not evidence of its absence: F's strong null alone (p = 0.15) had less power than E's draws of
the same null, which the repository already held. Before attributing a trend to one of two subgroups, test the difference between them —
here the two slopes do not differ — and use every draw of the same null and calibration.

## Round 80 (adversarial review of Proposition D‴.59 as revised after round 74, 2026-10-29)

**Design.** One referee, not modifying the repository. D3ca: Proposition D‴.59 as revised after round 74, with its code and results at
74e57ea (`review/REVIEW_round80_D3ca.md`). The referee:
- regenerated the md byte for byte, and reran GT through `kappa_from_path` on all six paths, GL, HLL and four A2 cells bit-identically;
- ran a closure fit from every stored loop state at 10⁻¹⁰ in (τ₃, τ₄), and derived the dissolution times and the loop's κ without the
  fine grid's offset;
- computed target 20's slow manifold and every coordinate's share of its induced metric, also for the pair alone and on a wide window;
- ran `hll_job` to τ = 10⁹ on the symmetric pair and two asymmetric ones, and `a2_job` to s = 10¹³;
- derived (a)'s fourth-order coefficient for signed pairs, and v²(E − E_inf)/K − 1 against v;
- ran 84 single-point mutants.

**Verdicts (D3ca).** Reproducibility, every number (one rounding), the second phase's law for both signs of K′, the O(ε^{1/3}) rate as a
bound, the collapse law and C (re-derived independently), (a) for signed pairs, the visibility counts and the rest hold. Three majors: for
K′ < 0 the second dipole phase is no approach to an infimum — the block lies below its optimal dipole fit by |K′|/v² and the pair
retracts; without τ₄ the dipole would dissolve at τ = 145.0, 154.1, 61.5 —; on the symmetric pair the loop keeps no cancellation — to
τ = 10⁹ the extra kernels decay at 10⁻⁶ and form a null cluster beside an exact one-sign fit at 10⁻⁸, HP's end the same state —; A2's 10A
"spike end" was the horizon's — the spike is J_N's σ → 0 plateau, escaped at s ≈ 2–6·10¹⁰. Minors: the collapse's finite-p correction
(the mean width's share p²/36, c̄ moving, E's p⁹ term from the neighbouring kernels); the loop's κ carries the fine grid's offset and two
targets switch sequence at larger ε; target 20's GT row is roundoff and nothing checks κ³'s sign; O(v⁻⁴), not O(v⁻³); 10
behaviour-changing surviving mutants. Cosmetic: literature, wording.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3ca-1 | major | K′ < 0 | `kappa_from_path`: κ only where K > 0 and the bracket is positive; for K′ < 0 the dissolution time τ₃ + K(τ₃ − τ_d)/\|K′\| (GT: 145, 154.1, 61.44 for targets 18, 12, 29, after the path's τ₄ 97.89, 81.64, 46.97); (ii)(c)(2), the proof and the Status: for K′ > 0 an approach to the dipole infimum from above, for K′ < 0 a retraction below the optimal dipole fit, the law holding while K(τ₃ − τ_d) + K′(τ − τ₃) > 0; "the residual of that block's optimal dipole fit" for "its asymptotic residual" |
| D3ca-2 | major | The symmetric pair | `fit_split` (HLL; at report time H, HP and HL): per feature the main kernel's error, and the others' Σ\|w\|, sum on the samples, widths and centres — a cancellation told from a null cluster; HLL to τ = 10⁹: the symmetric pair at 10⁻⁶ the others' Σ\|w\| 0.676 → 0.185, at 10⁻⁸ the main kernels' error 7.1·10⁻² → 1.4·10⁻⁵ with the others leaving the window (their sum 8.6·10⁻² → 2.2·10⁻⁵); A₂ 0.99 and 0.95 keep their cancellations (undecided); (ii)(b), the Status and New: on the symmetric pair the cancellation is the polish's |
| D3ca-3 | major | A2's 10A | S_A2 to 10¹³ and 20A added; `a2_spike`: the plateau's outputs, its E and the escape — 3A to s = 3·10⁵, 5A to 3·10⁷, 10A to 1–3·10¹⁰, then exact; 20A still on it at 10¹³; E 0.644, 0.906, 0.976 of ½⟨f*²⟩ whatever w —; the test asserts 10A's escape; (i), the Status and the Numerics |
| D3ca-4 | minor | The collapse's finite-p correction | GL: the local slope (d ln E/d ln Δc − 8)/p of the formed kernels' energy (0.951–0.963), the force's factor (d ln E/d ln Δc)/8, and with it and C the ratio 1 within −4.7·10⁻⁵ … +3.2·10⁻⁴; the mean width's share p²/36 (at most 0.00067); the proof: c̄ moves on the manifold (c* + M₃/(6Wσ²)) and the shares' motion cancels it, the arithmetic mean width moving a third as far as the weighted one, E even in Δc for a pair alone (the reflection about its feature), so the O(p) term is the neighbouring kernels' |
| D3ca-5 | minor | The loop's κ | `kappa_loop`: κ′ with v carried by the law to the path's τ₃ (the grid's offset, +0.10 … +0.58% at 10⁻¹², removed; within −0.27 … +0.21% of the path's κ), and † where the loop's intermediate slope lies beyond 1 ± 0.5 times the path's law (another sequence: target 27 at 10⁻⁷ … 10⁻⁹, 18 at 10⁻⁷); (3): κ converges once the loop follows the path's sequence, O(ε^{1/3}) a bound |
| D3ca-6 | minor | Target 20's GT row | `kappa_from_path`: a closure whose rms is below CLOSURE_EXACT = 10⁻¹² is exact — K′ = 0, no τ₄, κ = 1 (target 20: 4.4·10⁻¹⁷); GT prints every closure's rms and "— (exact)" |
| D3ca-7 | minor | O(v⁻⁴) | the proof and `k_third`'s docstring: E − E_inf = K/v² + O(v⁻⁴), the pair's sum even in 1/v about the midpoint; the test asserts the relative error ∝ 1/v² |
| D3ca-8 | minor | Surviving mutants | Tests: `gt_job` against GT's raw in every field (targets 9 and 20); `kappa_from_path` with its fits stubbed (the dissolution, no κ past it, an exact closure); `fine_phases` with a landing before τ₃, an exit before τ₄ and the midpoint state; `gl_job` at unequal widths; `h_stats` with an unformed invisible kernel; `hl_job`'s 1% threshold and ε reaching `loop_nd`, `hll_job`'s ε and split; `a2_job` at σ*/σ₀ = 0.25, its tolerances and 10A's escape; `fit_split` (each kernel main once) and `closure_fit`'s \|σ\|; of 21 single-point mutants of this round's code the tests kill 20, the last equivalent |
| D3ca-9 | cosmetic | Literature | Ho–Nguyen's r̄(2) = 4 for one extra component [snippet, secondary], their weights positive (for signed pairs, (a)); "the pair one Gaussian to third order, its fourth-order term not vanishing" |
| D3ca-10 | cosmetic | Wording | 0.7520–1.1145; `closure_fit` returns \|σ\|; GL's p (with \|Δc\|) and the landing's (with \|Δθ\|) named, differing by p²/72; "past τ₄, which target 20 lacks" |

**Lesson.** A claim about where a flow ends is a claim about its future: A2's 10A spike and the symmetric loop's cancellation were states
at their horizons, and two to four more decades showed the escape to the exact fit and the decay into a null cluster leaving the window.
A trend that has slowed has not stopped. Naming needs the same care: one v³ law covers both signs of K′, but only K′ > 0 approaches an
infimum, and for K′ < 0 the four stages survive only because τ₄ comes before the dipole dissolves. And an isolated pair is a cheap control
for an attribution — its energy is even in Δc, so the odd term that explains GL's remainder belongs to the neighbouring kernels.

## Round 81 (adversarial review of Proposition D‴.58 as revised after round 77, 2026-10-29)

**Design.** One referee, not modifying the repository. D3cb: Proposition D‴.58 as revised after round 77, with its code and results at
dd88e96 (`review/REVIEW_round81_D3cb.md`). The referee:
- regenerated the md byte for byte and reran all 12 B4 cells bit-identically, every field;
- evaluated all 48 B4 ends and 17 marks exactly (60-digit decimal residuals, long doubles as exact hi + lo) against `lm_rms`;
- instrumented the anchored solve (its re-anchors, the period's decrease, the jump at re-anchoring; without stall re-anchoring; with
  anchor_every = 500; the anchored residual in long double from its end) and predicted its floor from the frozen anchor's rounding;
- checked B3's labels and references, B2's Schur bound for Kaufman's and Golub–Pereyra's forms, and A5's frame series in 110-digit
  decimal to K = 40;
- ran 64 single-point mutants and the ten re-created round-77 survivors.

**Verdicts (D3cb).** Reproducibility, the anchored residual's identity, the central claim (stronger in exact arithmetic: the anchored ends
lie −1.385·10⁻⁸ … +2.2·10⁻¹⁶ off the long-double end and agree from both starts to 2·10⁻¹³), B3, B2 and the round-77 survivors (all killed
now) hold. Two majors: B4's numbers below ~10⁻¹⁰ measured the long-double rms — biased low at all 12 long-double ends, the function that
solve minimised — not the solves; nine behaviour-changing mutants of B4's new code survive behind the test's fakes. Minors: the solve
re-anchored thousands of times at stalls and never periodically, its 'flat' stop dead and "to its own relative accuracy" false along the
slow modes; the accurate solves' own floors unstated; A5's onset and its inference wrong, the rate computable; `b3_lines`' reference keyed
across variants; two references for the optimum below 2·10⁻⁹. Cosmetic: the trapezoid floor's digits, wording, `b3_job`'s guard.

| id | sev. | issue | fix (where) |
|---|---|---|---|
| D3cb-1 | major | B4 quoted below its instrument | `exact_resid`/`exact_rms` (the residual in 60-digit decimal, every double and long double taken exactly; `_dec` exact in any context): every B4 end and mark evaluated exactly, against its block's best end (`b4_best`); `lm_rms`'s error printed (−2.4·10⁻¹⁰ … +1.7·10⁻¹⁰, below the exact at 12 of 12 long-double ends); B4 rerun: double +1.9·10⁻⁸ … +2.3·10⁻³, long double up to +2.2·10⁻⁸, the anchored double within 1.3·10⁻¹² of the best end from all 12 starts, a block's two starts within 2.7·10⁻¹³; (iii), the Status and the Numerics restated; the test pins the evaluator against 90 digits and `lm_rms` at the 2.55·10⁻¹¹ optimum |
| D3cb-2 | major | Surviving mutants behind fakes | Tests: the real `b4_job` on the two fastest cells against the raw in every field; the fakes return distinct values per call and moved amplitudes (the long double's anchored solve started from the anchored end); the period counted across stall re-anchors (the exact rms evaluated once per period); `b3_lines`' 'path' line with a mismatching row |
| D3cb-3 | minor | The solve run, not the solve described | `lm_fit`: the 'flat' stop on the exact rms over anchor_every accepted steps counted across stall re-anchors (10 of 12 anchored solves now stop flat, after 1.5·10⁴–1.6·10⁵ iterations); `stall_reanchor=False` ('nostall', up to 9.8·10⁻⁹ above the best end: the anchored residual, not the restarts, passes the floor); the docstring and (iii): stall re-anchoring with the damping reset, the change's error proportional to the displacement |
| D3cb-4 | minor | The accurate solves' floors | `anchor_floor`: ‖P_Je₀‖²/(2‖r‖²) from the frozen anchor's rounding (1.8·10⁻¹⁶ … 8.8·10⁻¹³) against the long double's anchored solve from the anchored end (`lm_fit` in long double about an anchor), which gains 2.2·10⁻¹⁶ … 1.3·10⁻¹²; (iii) and the Status: the direct subtraction's floor linear in its rounding, the anchored one quadratic |
| D3cb-5 | minor | A5's onset and rate | `frame_series`, `gh_nodes`, `a5_series` (110 digits, K ≤ 40) and `a5_far` (140 digits, K = 60 … 120): the tail above the first term from K = 19, 8, 6 at τ/σ = 1/2, 1, 2 (not by 40 at 1/4); the Mercer frame's rule's ratio per kernel converging as K⁻² to 17.165, 5.3094, 2.33333 below 1/t; the series against the actual errors (to 6·10⁻¹⁴ … 3·10⁻⁵); the leading-order onset dropped; the Status: the slower fall shown at τ/σ = 1/2, 1, 2, its closed form open |
| D3cb-6 | minor | `b3_lines`' reference across variants | the reference B4's best end, absolute, over each row's own optimum (`b4_best`); a test with two variants' optima at one (τ/σ, K) |
| D3cb-7 | minor | Two references below 2·10⁻⁹ | one reference, the best end; (iii): the stored optimum 1.6·10⁻³ above it at 3/4, K = 11, the 0.1% test admitting ends up to 2.6·10⁻³ above it there |
| D3cb-8 | cosmetic | Loose ends | the trapezoid floor's digits at 1/4 and 1/2 marked; "from the path's blocks at those K"; md's "above B4_MAX a negative value …"; `b3_job`'s guard inside `b3_restart`, right after the 'path' restart and with its acceptance (a test: the larger restarts never run on a mismatch) |

**Lesson.** A claim about a floor needs an instrument finer than the floor. The long-double rms was the very function the long-double solve
minimised, so it measured that solve as better than it was, and every difference B4 printed below 10⁻¹⁰ was the instrument's. Evaluating
the residual exactly costs a tenth of a second an end and turns the claim into a measurement: the anchored residual reaches the optimum to
10⁻¹², its own floor quadratic in its frozen rounding, and it is the residual, not the restarts, that passes the direct subtraction's
floor. Describe the solver that ran: here thousands of stall re-anchors, not the periodic ones the text named.
