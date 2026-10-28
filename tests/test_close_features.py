"""Proposition D‴.52: close features — the second entry's delay ∝ Δ⁻² (a dipole's ∝ Δ⁻¹, the transit of its derivative shape), the
two-kernel block's stiffnesses ∝ Δ² and Δ⁶ in closed form, and the fast flow's rates along those modes."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import exp_close_features as CF
import exp_mirror_transit as X
from src.core import rich_transit as RT


def test_block_stiffnesses_are_the_closed_forms():
    """At the exact fit the Gauss–Newton Hessian's two small eigenvalues are λ_split = Δ⁶⟨(P⊥e₃)²⟩/8 — the same for every pair of
    amplitudes — and λ_sep = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²), to O(Δ)."""
    L, K, c0 = CF.setup()
    D = 0.0125
    split = []
    for A1, A2 in ((0.8, 0.6), (0.8, -0.4), (0.8, 0.8)):
        r = CF.b_job((A1, A2, D))
        l1, l2 = r['pred']
        assert abs(r['hess'][0] / l1 - 1) < 5e-3
        assert abs(r['hess'][1] / l2 - 1) < 5e-2
        assert r['hess'][2] > 10 * r['hess'][1]                           # the other two are of order one
        split.append(r['hess'][0] / D ** 6)
    assert np.ptp(split) / np.mean(split) < 1e-3                         # amplitude-free
    vs = CF.taylor_basis(K, CF.CBAR)                                     # the continuum's √π/(48σ⁵ℓ) and √π/(8σ³ℓ), ℓ = Nh (round 55)
    s, ell = K.sig, K.xs.size * (K.xs[1] - K.xs[0])
    assert abs(CF.perp_norm2(vs, 3) / (np.sqrt(np.pi) / (48 * s ** 5 * ell)) - 1) < 1e-4
    assert abs(CF.perp_norm2(vs, 2) / (np.sqrt(np.pi) / (8 * s ** 3 * ell)) - 1) < 1e-4


def test_flow_relaxes_the_split_at_the_mobility_weighted_stiffness():
    """The fast flow's slowest mode is the split, at λ_split·2/(1/μ₁ + 1/μ₂) (μ_k the amplitudes' mobility); the separation relaxes at
    Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) while the amplitudes are frozen on its scale; both modes identified by their eigenvectors."""
    r = CF.b_job((0.8, 0.6, 0.0125))
    (rs, os_), (rp, op) = r['split'], r['sep']
    ps, pp = r['pred_flow']
    assert os_ > 0.999 and op > 0.999
    assert abs(rs / ps - 1) < 1e-2
    assert abs(rp / pp - 1) < 5e-2
    assert rs == min(r['rates'])
    L = X.law()
    mu_branch = 4 * L.G * L.c * np.array([0.8, 0.6]) / L.A                 # μ ≈ 4Gc|w|/a on the exponential branch
    assert np.allclose(r['mu'], mu_branch, rtol=1e-3)


def test_second_entry_delay_is_the_leading_order():
    """A same-sign pair Δ = 0.025 apart: the second entry comes after Δτ₂ = min_j (s_j − u_j)/v_j, v_j = ½V⟨ψφ_j⟩, to 0.5%; the path ends at
    the exact two-kernel fit."""
    r = CF.a_job((0.8, 0.6, 0.025))
    (t1, j1), (t2, j2) = r['entries'][:2]
    dt, j, wall = r['pred']
    assert j2 == j and wall == 1
    assert abs((t2 - t1) / dt - 1) < 5e-3
    assert r['n_act'] == 2 and np.allclose(r['c_end'], r['cs'], atol=1e-6) and r['rms'] < 1e-10


def test_dipole_is_its_derivative_shape_slowed_by_the_separation():
    """A dipole (A₂ = −A₁) Δ = 0.025 apart: τ₁ and Δτ₂ are the derivative shape −A₁∂_cφ's, divided by Δ, to 0.5% and 1.5%, the same
    kernels entering."""
    r = CF.a_job((0.8, -0.8, 0.025))
    (t1, j1), (t2, j2) = r['entries'][:2]
    sh = r['shape']
    assert [j1, j2] == sh['kernels']
    assert abs(t1 / sh['tau1'] - 1) < 5e-3
    assert abs((t2 - t1) / sh['dtau2'] - 1) < 1.5e-2


def test_predicted_delay_scales_as_the_inverse_square_separation():
    """The leading-order delay is ∝ Δ⁻²: halving Δ multiplies it by 4 (to O(Δ))."""
    L, K, c0 = CF.setup()
    d = []
    for D in (0.02, 0.01):
        f = X.bumps_target(K, CF.pair_centres(0.8, 0.6, D), (0.8, 0.6))
        d.append(CF.predicted_delay(K, f, c0, 0.8, 0.6, D)[0])
    assert abs(d[1] / d[0] / 4 - 1) < 2e-2



def test_c_job_stops_at_the_second_entry_and_places_within_a_budget(monkeypatch):
    """C's draws stop at the second entry (the delay only); the first R_D under H1 are relaxed again after it within RHS_D right-hand-side
    evaluations — placed (two kernels, their centres sorted) or out of budget, never unbounded (round 50's runs: a near-coalescent
    block's relaxation exhausted the memory)."""
    monkeypatch.setattr(CF, 'R_C', 3)
    monkeypatch.setattr(CF, 'R_D', 2)
    r = CF.c_job((0.3, 0.003, 1, 7))
    rows = r['rows']
    assert all(np.isfinite(q['dt']) for q in rows) and [('oob' in q) for q in rows] == [True, True, False]
    assert all('c' in q and q['c'][0] < q['c'][1] and not q['oob'] for q in rows[:2])
    assert all(abs(c_ - t) < 0.01 for q in rows[:2] for c_, t in zip(q['c'], CF.pair_centres(*CF.PAIR_C, 0.3)))
    monkeypatch.setattr(CF, 'RHS_D', 50)
    rows = CF.c_job((0.3, 0.003, 1, 7))['rows']
    assert all(q['oob'] and 'c' not in q for q in rows[:2]) and all(np.isfinite(q['dt']) for q in rows)
    h0 = CF.c_job((0.3, 0.003, 0, 7))['rows']
    assert not any('oob' in q for q in h0)


def test_the_opposite_sign_pair_enters_at_the_other_wall():
    """Referee round 55: for (0.8, −0.4) at Δ = 0.0125 the predicted second entrant is kernel 13 at the wall −1 (the flank, the
    minority's sign), the path's too, its delay 0.907 of the leading order (the O(ϑ) correction, ϑ = 25Δ)."""
    r = CF.a_job((0.8, -0.4, 0.0125))
    (t1, j1), (t2, j2) = r['entries'][:2]
    dt, j, wall = r['pred']
    assert (j, wall) == (13, -1) and j2 == 13 and 0.9 < (t2 - t1) / dt < 0.92
    L, K, c0 = CF.setup()                                                # with the T₃ term: within 1.3%
    f = X.bumps_target(K, CF.pair_centres(0.8, -0.4, 0.0125), (0.8, -0.4))
    d3 = CF.predicted_delay(K, f, c0, 0.8, -0.4, 0.0125, order=3)
    assert d3[1:] == (13, -1) and abs((t2 - t1) / d3[0] - 1) < 0.013


def test_the_split_moves_the_centres_by_half_the_separation_per_weight():
    """Referee round 55: the flow's split eigenvector moves each centre by Δ/(2w_k) per unit of weight moved, to 0.2% at Δ = 0.0125."""
    for A in ((0.8, 0.6), (0.8, -0.4)):
        r = CF.b_job(A + (0.0125,))
        assert r['split'][1] > 0.999 and np.allclose(r['split_c'], 1, atol=2e-3)


def test_the_separation_has_two_regimes():
    """Referee round 55: the separation's exact rate (−D·JᵀJ) is the frozen-amplitude formula's at Δ = 0.025 and tends to the
    relaxed-amplitude rate Gm·λ_sep — 2/3 of the frozen one — at Δ = 0.0005."""
    for A in ((0.8, 0.6), (0.8, -0.4)):
        fr, rl, m = CF.sep_regimes(*A, 0.025)
        assert abs(fr - 1) < 5e-3
        fr, rl, m = CF.sep_regimes(*A, 0.0005)
        assert abs(rl - 1) < 1e-2 and abs(fr - 2 / 3) < 1e-2


def test_detection_on_a_synthetic_table():
    """Referee round 55: `detection`'s H0 quantile, power and Δ₅₀ (log-linear interpolation of the power at ½) on synthetic draws."""
    C, pw = [], {}
    Ds = sorted(CF.DELTAS_C)
    for k, D in enumerate(Ds):
        for s in CF.SIGMAS_C:
            p = min(1.0, 0.1 * k * (1 + CF.SIGMAS_C.index(s) % 2))
            n1 = int(round(p * 20))
            C += [dict(D=D, s=s, hyp=0, rows=[dict(dt=float(t)) for t in range(1, 201)]),
                  dict(D=D, s=s, hyp=1, rows=[dict(dt=5.0)] * n1 + [dict(dt=500.0)] * (20 - n1))]
            pw[(D, s)] = n1 / 20
    det, d50 = CF.detection(C)
    q = np.quantile(np.arange(1, 201), CF.ALPHA_C)
    assert all(abs(det[k_]['q0'] - q) < 1e-12 and det[k_]['power'] == pw[k_] for k_ in pw)
    for s in CF.SIGMAS_C:
        p = [pw[(D, s)] for D in Ds]
        k = next(i for i in range(len(Ds) - 1) if p[i] < 0.5 <= p[i + 1])
        t = (0.5 - p[k]) / (p[k + 1] - p[k])
        assert abs(d50[s] - np.exp(np.log(Ds[k]) + t * np.log(Ds[k + 1] / Ds[k]))) < 1e-12


def test_the_knot_law_is_the_paths_first_order():
    """Referee round 55: under one bump with noise sξ, `knot_of` — the first entry by the duals, then their rates ⟨Rξ, φ_j⟩ to the walls —
    is the path's second knot to first order: the relative difference 10⁻⁸ at s = 10⁻³, falling ∝ s; its law's quantile ∝ 1/s up to the
    margins' O(s) noise."""
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, 0.025), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    xi = np.random.default_rng(3).standard_normal(K.xs.size)
    err = []
    for s in (1e-3, 1e-4):
        R = RT.sliding_iss(L, K, w1 * K.phi(np.array([c1]))[:, 0] + s * xi, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [t for t, j, k, c in R['events'] if k == 'enter']
        err.append(abs((ent[1] - ent[0]) / CF.knot_of(K, c0, c1, w1, s * xi[None, :])[0] - 1))
    assert err[0] < 1e-7 and 5 < err[0] / err[1] < 20
    q = [np.quantile(CF.knot_law(K, c0, c1, w1, s, 4000, np.random.default_rng(1)), 0.05) * s for s in (1e-5, 1e-4, 1e-3)]
    d1, d2 = q[1] / q[0] - 1, q[2] / q[0] - 1                    # ∝ 1/s up to the margins' noise, O(s)
    assert abs(d1) < 1e-3 and 5 < d2 / d1 < 15


def test_one_fit_is_stationary():
    """Referee round 55: `one_fit` returns the single-kernel least-squares fit (the gradient in (w, c) vanishes), not one Newton step."""
    L, K, c0 = CF.setup()
    for D in (0.3, 0.025):
        f = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D), CF.PAIR_C)
        c, w = CF.one_fit(K, f)
        ph, dp = K.phi(np.array([c]))[:, 0], K.dphi(np.array([c]))[:, 0]
        r = f - w * ph
        assert abs(ph @ r) < 1e-12 and abs(w * (dp @ r)) < 1e-12


def test_the_delay_draws_are_noisy_bounded_and_the_hypotheses_differ(monkeypatch):
    """Referee round 55: the delay runs stop at the second entry (max_entries = 2); the draws differ (noise); under H0 the delay is far
    longer than under H1 at Δ = 0.3; Δτ₂ is the gap τ₂ − τ₁."""
    monkeypatch.setattr(CF, 'R_C', 4)
    monkeypatch.setattr(CF, 'R_D', 0)
    kws, real = [], RT.sliding_iss

    def spy(*a, **kw):
        kws.append(kw)
        return real(*a, **kw)
    monkeypatch.setattr(RT, 'sliding_iss', spy)
    d1 = [q['dt'] for q in CF.c_job((0.3, 0.003, 1, 7))['rows']]
    d0 = [q['dt'] for q in CF.c_job((0.3, 0.003, 0, 7))['rows']]
    assert all(k.get('max_entries') == 2 for k in kws) and len(kws) == 8
    assert len(set(d1)) == 4 and len(set(d0)) == 4 and max(d1) < 10 < 100 < min(d0)
    fake = dict(events=[(5.0, 10, 'enter', 0.0), (7.5, 11, 'enter', 0.1)], hist=[], act=np.zeros(21, bool), relaxed=[True])
    monkeypatch.setattr(RT, 'sliding_iss', lambda *a, **kw: fake)
    assert CF.c_job((0.3, 0.003, 0, 7))['rows'][0]['dt'] == 2.5


def test_placed_blocks_are_sorted_and_relaxed(monkeypatch):
    """Referee round 55: a placement run's block is recorded with its coefficients in its centres' order, and only if it relaxed."""
    monkeypatch.setattr(CF, 'R_C', 1)
    monkeypatch.setattr(CF, 'R_D', 1)
    act = np.zeros(21, bool); act[[3, 7]] = True
    c = np.zeros(21); c[3], c[7] = 0.5, -0.5
    w = np.zeros(21); w[3], w[7] = 1.0, 2.0
    delay = dict(events=[(5.0, 10, 'enter', 0.0), (7.5, 11, 'enter', 0.1)], hist=[], act=act, relaxed=[True])
    for rel, placed in (([True, True], True), ([True, False], False)):
        place = dict(events=delay['events'], hist=[], act=act, c=c, w=w, relaxed=rel, out_of_budget=False)
        monkeypatch.setattr(RT, 'sliding_iss', lambda *a, **kw: place if kw.get('max_rhs') else delay)
        row = CF.c_job((0.3, 0.003, 1, 7))['rows'][0]
        assert ('c' in row) == placed and not row['oob']
        if placed:
            assert row['c'] == [-0.5, 0.5] and row['w'] == [2.0, 1.0]


def test_the_report_on_synthetic_results():
    """Referee round 55: the report's cells — A's ratio and scaled delays, B's slopes and the split's centre components, C's powers with
    the Neyman–Pearson bound, D's errors over the Cramér–Rao bounds, the blocks near their features and the adjacent second entrants."""
    A = [dict(A=(0.8, 0.6), D=0.025, entries=[(2.0, 10), (1002.0, 11)], pred=(800.0, 11, 1), n_act=2, cs=[0.0, 0.04],
              c_end=[0.0, 0.04], w_end=[0.8, 0.6], rms=1e-14)]
    md = CF.report(A=A)
    L, K, c0 = CF.setup()
    p3 = CF.predicted_delay(K, X.bumps_target(K, CF.pair_centres(0.8, 0.6, 0.025), (0.8, 0.6)), c0, 0.8, 0.6, 0.025, order=3)
    assert f'| 1000 | 800 (11, +1) | 1.2500 | {1000 / p3[0]:.4f} | 0.625 | 25 |' in md and p3[1:] == (11, 1)
    B = [dict(A=(0.8, 0.6), D=D, hess=[1e-6 * (D / 0.05) ** 6, 1e-2 * (D / 0.05) ** 2, 1.0, 2.0], rates=[1.0, 2.0, 3.0, 4.0],
              pred=(1e-6, 1e-2), split=(3e-6 * (D / 0.05) ** 6, 0.9999), sep=(5.0 * (D / 0.05) ** 2, 0.999), pred_flow=(1.0, 1.0),
              split_c=[1.0012, 1.0011], mu=[0.1, 0.1]) for D in (0.05, 0.025)]
    md = CF.report(B=B)
    assert md.count('| 6.000, 2.000 |') == 2 and '; centres 1.0012, 1.0011)' in md
    rows = lambda dts, **kw: [dict(dt=float(t), j=[10, 11], **kw) for t in dts]
    C = []
    for D in CF.DELTAS_C:
        for s in CF.SIGMAS_C:
            C.append(dict(D=D, s=s, hyp=0, rows=rows(range(100, 300))))
            cs = CF.pair_centres(*CF.PAIR_C, D)
            C.append(dict(D=D, s=s, hyp=1, rows=rows([50] * 200, c=list(cs + [0.001, -0.001]), w=[0.8, 0.6], oob=False)))
    md = CF.report(C=C)
    L, K, c0 = CF.setup()
    for D in (0.3, 0.0125):
        cs = CF.pair_centres(*CF.PAIR_C, D)
        f1 = X.bumps_target(K, cs, CF.PAIR_C)
        c1, w1 = CF.one_fit(K, f1)
        npw = CF.np_power(K, f1, w1 * K.phi(np.array([c1]))[:, 0], 0.001)
        assert f"| {D:g} | 1.00 [{npw:.2f}] (50, 200) |" in md
        cp, cf_ = CF.crb(K, cs, np.array(CF.PAIR_C), 0.001)
        assert (f"| {D:g} | 0.001 | 200, 0 (200) | 1.000e-03 | {1e-3 * D ** 2 / 1e-3:.3g} | {1e-3 / cp:.3g} | 0.000e+00 | 0 | 0 | "
                f"200 | 40 |") in md
    assert "Under H0 the second entrant is the first's neighbour in 100.0% of the 8000 draws, in 100.0% of the earliest 5%" in md
    assert 'The largest coefficient of a placed block: Δ = 0.3, s = 0.001, w = (0.8, 0.6) (their sum 1.4)' in md
