# Transport Adaptive Control: Adaptive Kernel Approximators as Evolving Measures

**Copyright (c) 2026 Odwitiyo Dutta. All Rights Reserved.** See [LICENSE](LICENSE).

This repository holds a research program on continuous-time adaptive control with Gaussian RBF networks, whose
weights, centres and widths are all adapted online. It has two parts.

1. **The controller and its audit.** A Lyapunov-based adaptive RBF controller with a boundary-layer projection
   operator, together with an honest re-evaluation of its claims ([ASSESSMENT.md](ASSESSMENT.md),
   [THEORY_REVISED.md](THEORY_REVISED.md)).
2. **Transport adaptive control.** A new theory that treats the network as a probability measure on a lifted
   parameter space and the adaptive law as transport of that measure ([PARADIGM.md](PARADIGM.md)). Stability,
   robustness modifications, excitation, the many-kernel limit, the slow evolution of kernel geometry and the global
   placement of kernels all become statements on Wasserstein space, independent of how the kernels are labelled.

**Status.** This is research in progress and has not been peer reviewed.

- Every result in [PARADIGM.md](PARADIGM.md) is labelled *proved*, *formal* (a derivation without error control)
  or *conjectural*.
- Every inequality the tests can evaluate is checked numerically along simulated trajectories, and every number in
  the documents is produced by a script in [experiments/](experiments).
- The theory went through three rounds of automated multi-agent adversarial review. Every confirmed issue and its fix
  is logged in [review/REVIEW_LOG.md](review/REVIEW_LOG.md), and the raw reports are kept beside it.
- The claim-by-claim novelty audit is in [LITERATURE_TRANSPORT.md](LITERATURE_TRANSPORT.md).

---

## Part 2 (new): transport adaptive control

The network f̂ = (1/m)Σ a_i φ(x; c_i, s_i) is read as the empirical measure μ = (1/m)Σ δ_{ϑ_i} of its lifted
particles ϑ_i = (a_i, c_i, s_i). The output is linear in μ, and the adaptive law moves μ.

| result | statement | status |
|---|---|---|
| A | The adaptive law is the particle form of a Wasserstein gradient flow of a μ-linear potential | proved; interpretive |
| B | A label-free Lyapunov function ½‖e‖² + ½W(μ, μ*)², with a two-regime matched remainder | proved and checked |
| C | Optimal-transport σ-modification: exact dissipation, a gap-protected strict inequality, a local exponential estimate; the labelled leak reshuffles kernels | proved and checked |
| D | The mean-field closed loop: well-posedness, a Dobrushin estimate, the many-kernel limit | proved |
| E | Label-free local exponential convergence under tangential persistent excitation | proved, conditional on one cited theorem |
| F | The Wasserstein length of learning is bounded by the integrated tracking error | proved and checked |
| G | The slow kernel geometry follows variable projection in the metric induced on the slow manifold; plain variable projection needs a flatness condition | formal; checked, with a counterexample |
| H | With saturated features, σ-modification and parameter dither, the averaged loop is mean-field Langevin dynamics and converges globally to the regularised optimum | proved for the averaged system |
| H.7 | At finite feedback gain, the averaged loop is *not* a gradient flow: a semi-gradient (performative) correction, and an antisymmetric kernel interaction given by the area the reference orbit encloses | frozen-measure part proved; non-reciprocity checked |
| H.8 | The averaged system is input-to-state stable, which reduces the global conjecture for the real loop to an averaging estimate | proved under a stated regularity assumption |

**A new closed-loop effect.** Traversing the same reference orbit in the opposite direction leaves the visited
states unchanged, yet the learned kernels shift by an amount of order 1/k_e. The first-order theory, H.7 plus the
second-order averaging term of Remark H.7a, predicts that shift
([results/orbit_orientation/](results/orbit_orientation)):

| k_e | measured shift, size relative to prediction | cosine with prediction | vector error |
|---|---|---|---|
| 10 | 0.94 | 0.93 | 36% |
| 20 | 0.97 | 0.98 | 20% |
| 40 | 0.97 | 0.994 | 11% |

**Dither and global placement.** Without dither, a collapsed kernel configuration never separates. With dither,
collapsed and spread initialisations reach the same regularised loss, at the price of an entropic bias.

![Speed–bias trade-off of parameter dither](results/langevin_tradeoff/figure.png)

All numbers behind Part 2 are collected in [results/RESULTS_PARADIGM.md](results/RESULTS_PARADIGM.md).

---

## Part 1: the controller and what its experiments show

**What is claimed.**

- **Stability.** With the parameters confined to convex sets by a boundary-layer projection operator (Pomet and
  Praly 1992; Lavretsky and Wise 2013), the tracking error is uniformly ultimately bounded. Its long-run RMS is bounded
  by a radius that does not contain the leakage bias of σ- or e-modification. The argument combines known ingredients,
  so it is not a new theorem, and nothing is claimed about parameter convergence.
- **Implementation equals theory.** The implemented update laws are exactly the analysed ones. The Lyapunov identity
  and the projection inequality are checked along trajectories in the tests.

**Mechanism against geometry** (Lorenz tracking, m = 27, 5 seeds, same covering grid). Mean steady-state ‖e‖:

| geometry | projection | σ-mod | e-mod |
|---|---|---|---|
| fixed grid (W) | 0.158 | 0.733 | 0.562 |
| adaptive (W, c) | 0.147 | 0.400 | 0.233 |
| adaptive (W, c, σ) | 0.112 | 0.386 | 0.239 |

Projection beats both leaky modifications in every cell and every seed. Adaptive centres and widths add smaller,
consistent gains.

![Factorial](results/factorial_figure.png)

**Kernel-regression baselines** (exact GP, sparse GP, deep GP, NTK; n = 1, 2, 3, 6).

- The proposed method is never the most accurate. NTK is more accurate at n = 2 and n = 6, and a sparse GP at n = 3.
- Its advantage is cost: 0.1–0.6 ms per step, against 13–1900 ms per update for the baselines.
- In a fair re-run at n = 1, where every baseline predicts at every step, the exact GP is 3–4 times more accurate at
  100 times the per-step cost.

![All benchmarks](results/comparison_higher_dims_complete.png)

Earlier versions of this README stated that the method was "statistically indistinguishable from the Exact GP", that
all Monte Carlo trials satisfied "the calculated R_UUB", and that it "solved the curse of dimensionality". Those
statements are withdrawn: the significance test was invalid, the threshold was an arbitrary constant, and the omitted
n = 3 and NTK results contradict them.

---

## Repository layout

| path | content |
|---|---|
| [PARADIGM.md](PARADIGM.md) | the transport theory: results, proofs, verification map (§7), research program (§8) |
| [LITERATURE_TRANSPORT.md](LITERATURE_TRANSPORT.md) | claim-by-claim literature and novelty audit |
| [ASSESSMENT.md](ASSESSMENT.md), [THEORY_REVISED.md](THEORY_REVISED.md) | the audit of the original controller and its corrected theorem |
| [review/](review) | review log, raw adversarial-review reports for rounds 1–3, archived version 2 of the theory |
| `src/core/` | `transport.py` (laws (L), (L_σ), OT assignment), `langevin.py` (Theorem H), `averaging.py` (Proposition H.7), and the original controller, projection, kernel and plant code |
| `src/baselines/` | σ- and e-modification baselines |
| [tests/](tests) | 50 tests that evaluate the theorems' inequalities on simulated trajectories and sampled states |
| [experiments/](experiments) | one script per table or figure; outputs go to `results/<name>/` |
| [results/](results) | tables, figures and saved raw runs |

## Running

```
python -m venv venv
venv/bin/pip install -r requirements.txt        # Windows: venv\Scripts\pip install -r requirements.txt
python -m pytest tests                          # 50 tests, about 5–7 minutes
python experiments/exp_transport_permutation.py # each experiment writes results/<name>/
```

Set `OPENBLAS_NUM_THREADS=1` for the experiments that use a multiprocessing pool. Several long experiments save their
raw runs and can regenerate tables without re-simulating: `VP_REUSE=1` for `exp_variable_projection.py`,
`OO_REUSE=1` for `exp_orbit_orientation.py`, and `TR_REFIT=1` for `exp_langevin_tradeoff.py`.

## Continuing the research

The working conventions are part of the method:

- **Status labels.** Every claim is labelled proved, formal or conjectural, and the label changes only with a proof
  or a disproof.
- **Numbers come from scripts.** No table is edited by hand.
- **Test before changing.** A change to a theorem, a constant or an algorithm goes in only after a test or experiment
  supports it.
- **Review and log.** New results go through an adversarial review, and every confirmed issue goes into
  [review/REVIEW_LOG.md](review/REVIEW_LOG.md).

The open problems are listed in [PARADIGM.md](PARADIGM.md) §8. The most important are:

1. Infinite-horizon averaging for the McKean–Vlasov closed loop, which would settle Conjecture H.6, the real-loop
   counterpart of Theorem H.
2. The mean-field limit of the dithered law, and uniform-in-time propagation of chaos with a plant in the loop.
3. A general non-reciprocity theorem for H.7, beyond the checked instances, and the orientation experiment with
   dither.
4. Points the third review left uncertain: the exact hypotheses of the cited Morgan–Narendra theorem behind
   Theorem E, and the citations of Remark H.3.
