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
