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
