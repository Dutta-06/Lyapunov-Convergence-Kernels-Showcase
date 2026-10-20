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
