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

