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


def test_the_e_report_on_synthetic_draws():
    """E's rows and summary (rounds 55, 59): the draws' H0 quantile against the law's, the threshold fixed at c̄ = CBAR applied out of
    sample (below CBAR) and beyond, the calibrated test's size with its Clopper–Pearson interval and the powers, on hand-made draws (H0
    delays uniform on [a, 2a], the bootstrap p-values i/100)."""
    rng = np.random.default_rng(3)
    D, s = CF.CELLS_E[0]
    E_ = []
    for cb, a in ((0.0, 1.0), (CF.CBAR, 1.0), (0.049, 0.1)):
        d0 = a * (1 + rng.random(100))
        d1 = a * (0.5 + rng.random(100))
        E_.append(dict(D=D, s=s, cbar=cb, rows0=[dict(dt=float(x), p=i / 100, j=[5, 6]) for i, x in enumerate(d0)],
                       rows1=[dict(dt=float(x), p=0.0, j=[5, 7]) for x in d1], q_law=1.05 * float(np.quantile(d0, CF.ALPHA_C)),
                       ratio=np.full(21, 0.1), c1=cb, np=0.9))
    md = CF.report(E_=E_)
    lo, hi = CF.cp_int(5, 100)
    ref = E_[1]
    qref = np.quantile([q['dt'] for q in ref['rows0']], CF.ALPHA_C)
    fa0 = np.mean(np.array([q['dt'] for q in E_[0]['rows0']]) < qref)
    assert md.count(f'| {D:g}, {s:g} | ') == 3 and f'| 0.050, 0.050 [{lo:.3f}, {hi:.3f}] |' in md and f'| 1.000, 0.050 [{lo:.3f}, {hi:.3f}] |' in md
    assert (f"the false alarms of the threshold fixed at c̄ = {CF.CBAR:g} {fa0:.3f}–{fa0:.3f} at c̄ < {CF.CBAR:g} (out of sample) and "
            "1.000–1.000 beyond, of the calibrated test 0.050–0.050 (nominal 0.05)") in md
    pw = [np.mean(np.array([q['dt'] for q in e['rows1']]) < np.quantile([q['dt'] for q in e['rows0']], CF.ALPHA_C)) for e in E_]
    assert f"power minus the oracle's {1 - max(pw):+.3f} … {1 - min(pw):+.3f}, within 0.05 in 0 of 3" in md
    assert '| 0.100; 0.100, 0.100 |' in md                                       # the fit's kernel j₁, then j₁ ∓ 1


def test_e_job_on_a_few_draws(monkeypatch):
    """Referee round 59: `e_job` draws H0 from the one-bump fit and H1 from the pair, one rng per cell; each draw's p is its knot against
    the law bootstrapped at its own first fit (`knot_law` at (ĉ, ŵ), the noise s, the draw's own rng); the law at the truth and the
    std ratio as defined; the cells' seeds distinct."""
    monkeypatch.setattr(CF, 'R0_E', 2)
    monkeypatch.setattr(CF, 'R1_E', 2)
    monkeypatch.setattr(CF, 'M_E', 64)
    monkeypatch.setattr(CF, 'M_LAW', 400)
    seen = []
    real = RT.sliding_iss
    monkeypatch.setattr(RT, 'sliding_iss', lambda L, K, y, *a, **k: (seen.append(y.copy()), real(L, K, y, *a, **k))[1])
    D, s, cb, seed = CF.e_jobs()[7]
    r = CF.e_job((D, s, cb, seed))
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, D, cb), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    ys = [f0 + s * rng.standard_normal(K.xs.size) for _ in range(2)] + [f1 + s * rng.standard_normal(K.xs.size) for _ in range(2)]
    assert len(seen) == 4 and all(np.array_equal(a, b) for a, b in zip(seen, ys))
    for hyp, rows in ((0, r['rows0']), (1, r['rows1'])):
        for k, q in enumerate(rows):
            law = CF.knot_law(K, c0, q['c1'], q['w1'], s, 64, np.random.default_rng([seed, hyp, k]))
            assert q['p'] == np.mean(law <= q['dt'])
    assert r['q_law'] == np.quantile(CF.knot_law(K, c0, c1, w1, s, 400, np.random.default_rng([seed, 9])), CF.ALPHA_C)
    Rm, Phi = CF.tangent_projector(K, c1, w1), K.phi(c0)
    assert np.allclose(r['ratio'], np.sqrt(np.mean((Rm @ Phi) ** 2, axis=0) / np.mean(Phi ** 2, axis=0)))
    assert len({j[3] for j in CF.e_jobs()}) == len(CF.e_jobs()) and len({j[3] for j in CF.c_jobs()}) == len(CF.c_jobs())
    assert len({CF.f_seed(*j[:5]) for j in CF.f_jobs()}) == len(CF.f_jobs())


def test_crb_and_np_power_closed_forms():
    """Referee round 59: `crb` on an orthogonal design (JᵀJ diagonal) and `np_power` against Φ(‖f₁ − f₀‖/s − z₀.₉₅)."""
    from scipy.stats import norm

    class Ortho:
        xs = np.zeros(4)

        def phi(self, cs):
            return np.array([[2.0, 0.0], [0.0, 3.0], [0.0, 0.0], [0.0, 0.0]])

        def dphi(self, cs):
            return np.array([[0.0, 0.0], [0.0, 0.0], [5.0, 0.0], [0.0, 7.0]])
    w, s = np.array([0.8, 0.6]), 0.3
    pos, split = CF.crb(Ortho(), None, w, s)
    W = w.sum()
    assert np.isclose(pos, s * np.sqrt((1 / (25 * w[0] ** 2) + 1 / (49 * w[1] ** 2)) / 2), rtol=1e-12)
    assert np.isclose(split, s * np.sqrt(w[1] ** 2 / W ** 4 / 4 + w[0] ** 2 / W ** 4 / 9), rtol=1e-12)
    L, K, c0 = CF.setup()
    d = np.zeros(K.xs.size)
    d[3] = 1.0
    for a in (0.0, 1.0, 2.5):
        f1 = (norm.ppf(0.95) + a) * 0.01 * d
        assert np.isclose(CF.np_power(K, f1, 0 * d, 0.01), norm.cdf(a), rtol=1e-12)


def test_d50_bootstrap_on_a_deterministic_table():
    """Referee round 59: a power that steps from 0 to 1 between two grid separations crosses ½ at their geometric mean, in every
    resample; `crossing` interpolates in log Δ."""
    Ds = sorted(CF.DELTAS_C)
    C = []
    for s in CF.SIGMAS_C:
        step = Ds[3] if s < 0.01 else Ds[5]
        for D in Ds:
            C.append(dict(D=D, s=s, hyp=0, rows=[dict(dt=float(v)) for v in np.linspace(10.0, 20.0, 20)]))
            C.append(dict(D=D, s=s, hyp=1, rows=[dict(dt=1.0 if D >= step else 50.0) for _ in range(20)]))
    ci, sl = CF.d50_bootstrap(C, n=20)
    for s in CF.SIGMAS_C:
        k = 3 if s < 0.01 else 5
        g = np.sqrt(Ds[k - 1] * Ds[k])
        assert np.isclose(ci[s][0], g) and np.isclose(ci[s][1], g) and ci[s][2] == 0
    assert np.isclose(CF.crossing([1.0, 2.0, 4.0], [0.2, 0.4, 0.6]), 2 * np.sqrt(2))
    assert np.isnan(CF.crossing([1.0, 2.0], [0.1, 0.3]))


def test_the_knot_law_on_several_noise_vectors():
    """Referee round 59: `knot_of` is the path's second knot to first order over several noise vectors — the error ∝ s."""
    L, K, c0 = CF.setup()
    f1 = X.bumps_target(K, CF.pair_centres(*CF.PAIR_C, 0.025), CF.PAIR_C)
    c1, w1 = CF.one_fit(K, f1)
    y0 = w1 * K.phi(np.array([c1]))[:, 0]
    err = {}
    for s in (1e-3, 1e-4):
        e = []
        for k in range(6):
            xi = np.random.default_rng(100 + k).standard_normal(K.xs.size)
            R = RT.sliding_iss(L, K, y0 + s * xi, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
            ent = [t for t, j, kk, c in R['events'] if kk == 'enter']
            e.append(abs((ent[1] - ent[0]) / CF.knot_of(K, c0, c1, w1, s * xi[None, :])[0] - 1))
        err[s] = np.array(e)
    assert err[1e-3].max() < 1e-2 and 4 < np.median(err[1e-3]) / np.median(err[1e-4]) < 25


def test_the_score_statistic():
    """Referee round 59: z = ⟨r, R∂²_cφ⟩/(s‖R∂²_cφ‖) at a one-bump fit — zero on a one-bump target; on a close pair's one-bump fit the
    signal ‖f₁ − f₀‖ (r ∥ R∂²_cφ) with the pair's sign sign(A₁A₂); N(0, 1) under one bump with noise; the width's z_σ equal to it to first
    order; σ̂² − σ² → V/A on a close pair (the pair's variance added)."""
    L, K, c0 = CF.setup()
    y = 0.9 * K.phi(np.array([0.03]))[:, 0]
    assert abs(CF.score_stats(K, y, 0.03, 0.9, 0.01)[0]) < 1e-9
    for A in ((0.8, 0.6), (0.8, -0.4)):
        D = 0.01                                                       # ϑ = 0.12, 0.25: r ∥ R∂²φ to O(ϑ²)
        f1 = X.bumps_target(K, CF.pair_centres(*A, D), A)
        c1, w1 = CF.one_fit(K, f1)
        z = CF.score_stats(K, f1, c1, w1, 1.0)[0]
        assert np.sign(z) == np.sign(A[0] * A[1]) and abs(z) / np.linalg.norm(f1 - w1 * K.phi(np.array([c1]))[:, 0]) > 0.99
        sg, sd = CF.width_fit(K, f1, c1, w1)
        V_A = A[0] * A[1] * D ** 2 / sum(A) ** 2
        assert abs((sg ** 2 - K.sig ** 2) / V_A - 1) < 0.05
    rng = np.random.default_rng(8)
    f0 = 1.4 * K.phi(np.array([0.032]))[:, 0]
    zs, ws = [], []
    for _ in range(600):
        yy = f0 + 0.01 * rng.standard_normal(K.xs.size)
        c, w = CF.one_fit(K, yy)
        zs.append(CF.score_stats(K, yy, c, w, 0.01)[0])
        ws.append(CF.width_z(K, yy, c, w, 0.01))
    zs, ws = np.array(zs), np.array(ws)
    assert abs(zs.mean()) < 0.15 and 0.9 < zs.std() < 1.1
    assert np.abs(ws - zs).max() < 0.05


def test_the_knot_per_wall():
    """Referee round 59: `knot_walls`' minimum is `knot_of`; `walls_at_fit` at the path's own first fit is the path's second knot."""
    L, K, c0 = CF.setup()
    Xi = 0.01 * np.random.default_rng(1).standard_normal((300, K.xs.size))
    a, b = CF.knot_walls(K, c0, 0.032, 1.4, Xi)
    assert np.array_equal(np.minimum(a, b), CF.knot_of(K, c0, 0.032, 1.4, Xi))
    f1 = X.bumps_target(K, CF.pair_centres(0.8, -0.4, 0.025), (0.8, -0.4))
    for k in range(3):
        y = f1 + 0.01 * np.random.default_rng(20 + k).standard_normal(K.xs.size)
        R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [t for t, j, kk, c in R['events'] if kk == 'enter']
        h = R['hist'][1]
        i = int(np.flatnonzero(h[2])[0])
        tA, tB = CF.walls_at_fit(K, c0, y, float(h[1][i]), float(h[2][i]))
        assert abs(min(tA, tB) / (ent[1] - ent[0]) - 1) < 1e-9


def test_f_and_s_jobs_on_a_few_draws(monkeypatch):
    """Referee round 59: `f_job`'s rows — the knot per wall at the path's first fit, the p-values against `knot_walls` at that fit with the
    draw's rng, the law at the truth with its own; `s_job` regenerates C's labels in the job's order (its scores are those of `c_job`'s
    labels at their one-bump fits)."""
    monkeypatch.setattr(CF, 'M_E', 64)
    monkeypatch.setattr(CF, 'M_LAW', 400)
    A, cb, s, D = (0.8, -0.4), 0.032, 0.01, 0.025
    r = CF.f_job((A, cb, s, D, 1, 3))
    L, K, c0 = CF.setup()
    seed = CF.f_seed(A, cb, s, D, 1)
    N = K.xs.size
    for k, q in enumerate(r['rows']):
        assert abs(min(q['tA'], q['tB']) / q['dt'] - 1) < 1e-9
        la, lb = CF.knot_walls(K, c0, q['c1'], q['w1'], s * np.random.default_rng([seed, 1, k]).standard_normal((64, N)))
        assert q['p'] == np.mean(np.minimum(la, lb) <= q['dt']) and q['pA'] == np.mean(la <= q['tA']) and q['pB'] == np.mean(lb <= q['tB'])
    f1 = X.bumps_target(K, CF.pair_centres(*A, D, cb), A)
    c1, w1 = CF.one_fit(K, f1)
    la, lb = CF.knot_walls(K, c0, c1, w1, s * np.random.default_rng([seed, 0]).standard_normal((400, N)))
    assert r['q_law'] == np.quantile(np.minimum(la, lb), CF.ALPHA_C) and r['qB'] == np.quantile(lb, CF.ALPHA_C / 2)
    seen = []
    monkeypatch.setattr(RT, 'sliding_iss', lambda L_, K_, y, *a, **kw: (seen.append(y.copy()), dict(events=[]))[1])
    Dc, sc, hc, sdc = CF.c_jobs()[3]
    monkeypatch.setattr(CF, 'R_C', 3)
    CF.c_job((Dc, sc, hc, sdc))
    out = CF.s_job(('C', Dc, sc, CF.CBAR, hc, sdc, (3,), None))
    for y, q in zip(seen, out[f'rows{hc}']):
        c, w = CF.one_fit(K, y)
        assert q['z'] == CF.score_stats(K, y, c, w, sc)[0]


def test_the_f_report_on_synthetic_jobs():
    """Referee round 59: F's report — Δ₅₀ per test from the rejections, '<'/'>' at the grid's ends, the false alarms at H0."""
    s, cb = CF.SIGMAS_F[1], CF.CBAR
    Ds = CF.DELTAS_F[s]
    F = []
    for A in CF.PAIRS_F:
        sg = np.sign(A[0] * A[1])
        for D in Ds:
            k = Ds.index(D)
            rows = [dict(dt=1.0 if i < 50 * k else 9.0, p=0.5, tA=9.0, tB=9.0, z=sg * (3.0 if i < 40 * k else 0.0),
                         zs=sg * (3.0 if i < 40 * k else 0.0), w1=1.0) for i in range(200)]
            F.append(dict(A=A, cbar=cb, s=s, D=D, hyp=1, rows=rows, q_law=2.0, qA=1.0, qB=1.0, np=min(1.0, 0.2 * (k + 1))))
        rows0 = [dict(dt=1.0 if i < 20 else 9.0, p=0.01 if i < 10 else 0.5, tA=9.0, tB=9.0, z=0.0, zs=0.0, w1=1.0) for i in range(400)]
        F.append(dict(A=A, cbar=cb, s=s, D=Ds[CF.H0_F], hyp=0, rows=rows0, q_law=2.0, qA=1.0, qB=1.0, np=0.0))
    md = '\n'.join(CF.f_lines(F))
    knot = CF.crossing(list(Ds), [min(1.0, 0.25 * j) for j in range(len(Ds))])
    z = CF.crossing(list(Ds), [min(1.0, 0.2 * j) for j in range(len(Ds))])
    npd = CF.crossing(list(Ds), [min(1.0, 0.2 * (j + 1)) for j in range(len(Ds))])
    row = f"| {s:g} | {cb:g} | {knot:.4g}, > {Ds[-1]:g}, > {Ds[-1]:g}, {z:.4g}, {z:.4g}, {z:.4g}, {npd:.4g} | 0.050, 0.025, 0.000, 0.000, 0.000, 0.000 |"
    assert md.count(row) == 2
    assert '## F. Resolution across the grid cell, for both signs' in md
