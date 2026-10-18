# Assessment of "Self-Adaptive Kernels Achieving Lyapunov Convergence"

Scope: the root workspace (`lyapunov.pdf`, `src/`, `experiments/`, `tests/`,
`results/`, `public_showcase/`).  The AAAI and AISTATS folders were excluded on
request.  Every numerical statement below was regenerated in this revision;
the scripts and raw files are named next to each claim.

## 1. Verdict in one paragraph

The theory is a correct but standard composition of known results
(Taylor-linearised inner-layer adaptation, boundary-layer projection, level-set
UUB), presented with language that overstates its novelty and its conclusion.
Two of its claims are wrong as written (non-convex "projection" for centre
separation; an IFT argument that is not a proof) and one was not implemented
(the Γ-weighted projection).  The experiments contained a hand-chosen
"theoretical" threshold, a confounded headline comparison against baselines
that failed representationally, GP baselines handicapped by stale predictions
and a mis-placed inducing grid, invalid significance tests on autocorrelated
samples, and the public figure omitted the one benchmark (n=3) where the method
is weak and the one baseline (NTK) that beats it.  After revision, the
implemented algorithm is exactly the proved one and is machine-checked; the
honest experiments show a real, consistent advantage of projection over σ- and
e-modification on identical geometry (5/5 seeds in every cell) and a modest but
consistent gain from adaptive geometry, while the a-priori bound is shown to be
three to four orders of magnitude loose.

## 2. Novelty

| Ingredient | Prior source | Status of the manuscript's claim |
|---|---|---|
| Adaptive RBF control with projection-type robustification | Sanner & Slotine 1992; Ioannou & Sun 1996 | known |
| Adapting centres/widths by first-order Taylor expansion with bounded remainder | Lewis, Yesildirek & Liu 1996; Polycarpou 1996; Fabri & Kadirkamanathan 1996; Ge & Wang 2002 | known; the manuscript cites none of these |
| Smooth projection with ε boundary layer and inequality (P) | Pomet & Praly 1992; Lavretsky & Gibson 2011; Lavretsky & Wise 2013 | known; the showcase README calls it "the core theoretical contribution" |
| "Projection introduces no bias in V̇ when θ* ∈ Ω" | Ioannou & Sun 1996 §4.4 | known |
| Level-set ("global interior") UUB argument | Khalil Thm 4.18; Lewis et al. 1999 | known; called "critical and non-trivial" |
| Simultaneous adaptation of W, c, σ for Gaussian RBFs with all three projected | — | the only combinational novelty; minor |

Net: no new theorem.  What this repository can honestly claim after revision is
in `THEORY_REVISED.md` §7: a machine-checked theorem-to-code map, an explicit
Hessian-free disturbance bound, the RMS form of the bound, and a controlled
mechanism-vs-geometry experiment.

## 3. Integrity findings (each verified by running code)

| # | Finding | Evidence | Fix |
|---|---|---|---|
| I1 | Implemented centre update ≠ proved update: a pairwise repulsion term broke the projection inequality the proof needs | `tests/test_theorem_to_code.py` (oracle plant): inequality violated at 6054 of 20000 samples with the legacy term, 0 of 20000 without; the term also failed its own purpose (min separation 0.07 vs δ=0.5) | removed from the algorithm (`enforce_separation=False` default); tracking unchanged on the headline run (mean SS error 0.0767 → 0.0774, max 0.46 → 0.40) |
| I2 | Projection applied the isotropic formula to any Γ; the proof needs the Γ-weighted operator | `tests/test_controller_paths.py` | Γ-weighted projection implemented; fast/dense paths agree to 1e-12 |
| I3 | "R_UUB threshold" in the Monte Carlo figure was the constant 0.5 | `tier1_monte_carlo.py` line 69: "A chosen arbitrary bound" | `exp_monte_carlo_honest.py`: only exactly checkable statements (200/200 hard-bound and ISS checks) |
| I4 | "Theoretical R_UUB" in the k_e sweep used invented constants and an underived formula | `tier2_ruub_sweep.py`, former `estimate_D_max` | `exp_oracle_ruub.py`: all constants computable; a-priori bound 900×–11 000× loose; RMS bound with measured disturbance holds 48/48 |
| I5 | Headline figure: fixed-grid baselines with a grid not covering the reference; their error (≈15) exceeded open loop, so the figure showed representational failure, not "bias" | `results/tier1_lorenz_tracking.png` vs ablation figure (projection W-only fails identically) | `exp_bias_factorial.py`: 3 mechanisms × 3 geometries, covering grid, 5 seeds |
| I6 | Wilcoxon tests on ~12 000 time samples with lag-1 autocorrelation 0.9999; README's "p > 0.3 vs Exact GP" is contradicted by its own data (p = 1.6e-11) | recomputed from `results/*.npy` | seed-level statistics only |
| I7 | n=3 benchmark and NTK baseline removed from the public figure (commits "exclude n=3 and NTK") | `results/n3/*.npy`, `results/*/ntk_*.npy`; `plot_comparison_hd_complete.py` | complete table regenerated: proposed ranks 3/5, 2/5, 3/5, 2/5 |
| I8 | GP baselines held a stale prediction for 50–100 steps; sparse GP inducing points on a diagonal line over x ∈ [−3,3] while the reference spans ±20 | `sim1`–`sim3` | `kernel_baselines_fair.py`: per-step prediction, 9×3 covering grid, seeds |
| I9 | "O(mn) per step" claim false with the O(m²n) repulsion | `projection.py` | true after I1 |
| I10 | "Concurrent Learning" hook used the current residual only | former `controller.py` | removed |
| I11 | Lyapunov monitor requires θ*; never used by any experiment, so "Lyapunov descent verification" was never run | `integrator.py` | oracle module + tests perform it on a representable plant |

## 4. How it stacks up (regenerated numbers)

Mechanism × geometry, Lorenz tracking, 5 seeds, mean steady-state ‖e‖
(`results/factorial/table.md`):

| geometry | projection | σ-mod | e-mod |
|---|---|---|---|
| fixed grid (W) | 0.158 | 0.733 | 0.562 |
| adaptive (W, c) | 0.147 | 0.400 | 0.233 |
| adaptive (W, c, σ) | 0.112 | 0.386 | 0.239 |

Projection beats both leaky baselines in every cell, 5/5 seeds, paired CIs
excluding zero.  Within projection, adaptive centres improve on the fixed grid
by 0.012 (5/5) and adaptive widths by a further 0.035 (5/5).  So the method's
real advantage is the robustness mechanism, with geometry adaptation a smaller
second-order gain — close to the reverse of the original figure's story.

Against kernel-regression baselines (saved single runs, nothing excluded,
`results/comparison_complete_table.md`): the proposed method is never first;
NTK is first at n=2 and n=6 by 10–25×, sparse GP is first at n=3 where the
proposed method's error is 7.7.  Its genuine edge is cost: 0.1–0.6 ms per step
against 13–1900 ms per update for the baselines.

Fair re-run at n=1 (per-step prediction for all methods, covering inducing
grid, 3 seeds; `results/fair_n1/table.md`): Exact GP 0.059, Sparse GP 0.187,
proposed 0.198, Deep GP 0.309, NTK 0.470 (mean steady-state ‖e‖); median
online step 6.3 / 1.35 / 0.06 / 4.5 / 0.21 ms.  The exact GP is 3–4× more
accurate than the proposed method in every seed once it is allowed to predict
every step; the proposed method is on par with a fixed-grid sparse GP at 20×
lower cost.  This is the version to quote.

Oracle bound study (`results/oracle_ruub/table.md`): the Lyapunov identity
holds to 3.6e-15, property (P) is never violated, the RMS bound with the
measured disturbance holds in 48/48 runs and tracks 1/k_e, and the a-priori
R_UUB is 26–261 against observed RMS errors of 0.01–0.09.

## 5. What changed in this revision

* `src/core/projection.py`, `controller.py`, `rbf_kernel.py`: proved algorithm only; Γ-weighted projection; explicit bound; validation of Γ; `k_cl` removed; legacy repulsion opt-in.
* `src/core/oracle.py`: representable plant, exact V, V̇, d, cross; RMS bound.
* `src/baselines/leaky_geometry.py`: σ-/e-modification on the same adaptive architecture.
* `tests/`: 20 tests (projection inequality incl. non-isotropic Γ, Jacobian suprema, path agreement, theorem-to-code identity and property (P) along trajectories, RMS bound, hard-bound invariants).  Gate: `venv\Scripts\python -m pytest tests`.
* `experiments/`: `exp_bias_factorial.py`, `exp_monte_carlo_honest.py`, `exp_oracle_ruub.py`, `kernel_baselines_fair.py` + `run_fair_n1.py` + `plot_fair_n1.py`, `plot_comparison_hd_complete.py`.  Old scripts kept untouched for reproduction of the original figures.
* `THEORY_REVISED.md`: assumptions, lemmas, theorem (i)–(iv), errata E1–E11, prior art.
* `public_showcase/README.md`: claims corrected locally (not committed, not pushed).
* Originals backed up in the session scratchpad (`backup_original/`).

## 6. Recommendation

Do not submit the manuscript in its current form.  Re-frame as an
empirical-plus-verified-implementation paper: the factorial result is a clean,
reproducible finding; the theory section should be `THEORY_REVISED.md`
Sections 1–5 with the prior work cited.  Drop the "zenith" language
("absolute mathematical certainty", "curse of dimensionality solved",
"structural result"); none of it survives the data.
