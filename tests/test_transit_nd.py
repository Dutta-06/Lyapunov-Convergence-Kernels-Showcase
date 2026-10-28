"""Proposition D‴.54: the transit in ℝⁿ (`src/core/rich_transit_nd.py`, `experiments/exp_transit_nd.py`)."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

from src.core import rich_transit as RT
from src.core import rich_transit_nd as RN
import exp_cohesion_leak as E
import exp_mirror_transit as X
import exp_transit_mass as TM
import exp_transit_nd as TN


def _grid():
    xs1, xs2 = np.linspace(-1, 1, 17), np.linspace(-1, 1, 15)
    return RN.BumpsGrid2(xs1, xs2, 0.25), RN.BumpsND(RN.BumpsGrid2(xs1, xs2, 0.25).xs, 0.25)


def test_separable_grid_equals_the_generic_dictionary():
    KG, KN = _grid()
    rng = np.random.default_rng(1)
    c, w, r = rng.uniform(-0.7, 0.7, (5, 2)), rng.normal(size=5), rng.normal(size=KG.N)
    assert np.abs(KG.synth(c, w) - KN.synth(c, w)).max() < 1e-14
    assert np.abs(KG.dphi(c) - KN.dphi(c)).max() < 1e-14
    for a, b in zip(KG.corr(c, r, second=True), KN.corr(c, r, second=True)):
        assert np.abs(a - b).max() < 1e-13
    for a, b in zip(KG.grams(c), KN.grams(c)):
        assert np.abs(a - b).max() < 1e-14
    # the second derivative against finite differences of the first
    g, h, Q = KN.corr(c, r, second=True)
    for k in range(2):
        e = np.zeros(2); e[k] = 1e-6
        dh = (KN.corr(c + e, r)[1] - KN.corr(c - e, r)[1]) / 2e-6
        assert np.abs(dh - Q[:, :, k]).max() < 1e-7 * max(1.0, np.abs(Q).max())


def test_the_jacobian_is_exact():
    L = X.law()
    KG, _ = _grid()
    rng = np.random.default_rng(2)
    c = rng.uniform(-0.6, 0.6, (4, 2))
    f = KG.synth(np.array([[0.1, -0.2], [-0.4, 0.3]]), np.array([0.8, -0.5]))
    for b0 in (0.0, 1e-3):
        fr, jr = RN.critical_rhs_nd(L, KG, f, b0), RN.critical_jac_nd(L, KG, f, b0)
        z = np.r_[rng.normal(size=4) * 1.5, c.ravel()]
        J = jr(0, z)
        Jfd = np.stack([(fr(0, z + e) - fr(0, z - e)) / 2e-6 for e in 1e-6 * np.eye(z.size)], 1)
        assert np.abs(J - Jfd).max() < 1e-7 * np.abs(J).max()


def test_one_dimension_is_the_one_dimensional_library():
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, [-0.43, 0.12, 0.61], [0.8, -0.6, 0.9])
    K1 = RN.BumpsND(K.xs[:, None], K.sig)
    R1 = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    R2 = RN.sliding_iss_nd(L, K1, f, c0[:, None], np.ones(c0.size), polish=True)
    assert [(round(e[0], 8), e[1], e[2]) for e in R1['events']] == [(round(e[0], 8), e[1], e[2]) for e in R2['events']]
    assert np.abs(R1['w'] - R2['w']).max() < 1e-10 and np.abs(R1['c'] - R2['c'][:, 0]).max() < 1e-10
    taus = np.linspace(0, 25, 26)
    _, W1, C1, _ = X.run_moving(L, K, f, c0, 1e-6, 25, 'paradigm', taus=taus)
    _, W2, C2, _, _ = RN.loop_nd(L, K1, f, c0[:, None], 1e-6, taus)
    assert np.abs(W1 - W2).max() < 1e-7 and np.abs(C1 - C2[:, :, 0]).max() < 1e-7


def test_coinciding_kernels_split_by_the_profile_and_ties_equally():
    L = X.law()
    K, c0 = TN.local()
    fb = 1.0 * K.phi(np.array([[TN.H / 2, TN.H / 2]]))[:, 0]
    for x in ([0.0, 1.0], [0.0, 0.5, 2.0], [0.0, 0.3, 1.1, 3.0]):
        out = RN.inner_k(L, K, fb, np.tile([TN.H / 2, TN.H / 2], (len(x), 1)), np.array(x))
        ps = RN.profile_split(L, 1.0, x)
        assert out['stationary'] and np.abs(out['shares'] - ps).max() < 1e-6
        # the exponential branch: the softmax, each share to O(e^{−2α_j/a}) — small for these shares, not at the threshold (W5)
        assert np.abs(ps - TN.softmax(x, L)).max() < 1e-4
    # a symmetric tie of distinct kernels: equal shares
    for kind, k in (('edge', 2), ('vertex', 4)):
        pts, T = TN.TIES[kind]
        out = RN.inner_k(L, K, K.phi(np.atleast_2d(T))[:, 0], np.array(pts, float), np.zeros(k))
        assert np.abs(out['shares'] - 1 / k).max() < 1e-8


def test_the_vertex_gaps_are_additive_and_the_loop_reaches_the_inner_problem():
    L = X.law()
    K, c0 = TN.local()
    cs, idx = TN.place(K, c0, 1.0, 'vertex', (1.0, 3.0), 1e-10)
    mu, _ = TN.margins(K, c0, 1.0, cs, idx)
    x = mu / 1e-5
    assert abs(x[1] - 1) < 1e-8 and abs(x[2] - 3) < 1e-8 and abs(x[3] - 4) < 1e-3
    r = TN.a_job(('vertex', 1.0, (1.0, 3.0), 1e-9))
    ri = TN.inner_job(('vertex', 1.0, tuple(np.r_[0.0, r['x_act'][1:]])))
    assert ri['stationary'] and np.abs(r['shares'] - ri['shares']).max() < 2e-4
    # not the softmax; the edges' product in the vertex's geometry (the feature H/2 off each edge's axis) to 3·10⁻⁴, on the edges'
    # axes (round 42's artefact) further off
    assert np.abs(TN.softmax(np.r_[0.0, r['x_act'][1:]], L) - ri['shares']).max() > 1e-2
    prod = {}
    for t in (TN.H / 2, 0.0):
        p1, p2 = (TN.inner_job(('edge', 1.0, (0.0, r['x_act'][i]), t))['shares'][1] for i in (1, 2))
        prod[t] = np.array([(1 - p1) * (1 - p2), p1 * (1 - p2), (1 - p1) * p2, p1 * p2])
    assert np.abs(prod[TN.H / 2] - ri['shares']).max() < 3e-4 < np.abs(prod[0.0] - ri['shares']).max()


def test_the_diagonal_gap_is_exactly_the_product_of_the_edges_correlations():
    """On tensor samples 1 − μ₁₁ = (1 − μ₁₀)(1 − μ₀₁): x₁₁ = x₁₀ + x₀₁ − ε^{1/2}x₁₀x₀₁ exactly (round 42: not 'to 3·10⁻⁴')."""
    K, c0 = TN.local()
    for eps in (1e-6, 1e-10):
        for x1, x2 in ((1.0, 3.0), (4.0, 4.0), (2.0, 5.0)):
            xd = TN.margins_diag(K, c0, x1, x2, eps)
            assert abs(xd - (x1 + x2 - np.sqrt(eps) * x1 * x2)) < 1e-6


def test_product_rule_and_areas_on_synthetic_tables():
    """`product_count` multiplies the edges' shares at the offsets given (interpolated in |t|); `areas` integrates N ≥ 3 and N = 4 by
    the trapezoid rule in the gaps' quarter-plane."""
    T = dict(x=np.array([0.0, 4.0, 8.0]), t=np.array([0.0, TN.H / 2]), pt=np.array([[0.5, 0.12, 0.0], [0.5, 0.05, 0.0]]))
    T['p'] = T['pt'][0]
    assert TN.product_count(T, 0.0, 0.0) == 4                       # four shares ¼
    assert TN.product_count(T, 0.0, 8.0) == 2                       # one edge tied, the other far
    assert TN.product_count(T, 4.0, 4.0, 0.0, 0.0) == 4             # on the edges' axes the diagonal's share 0.12² is above 1% …
    assert TN.product_count(T, 4.0, 4.0) == 3                       # … at the vertex's offset 0.05² is not
    assert abs(float(TN.p_edge(T, 2.0, TN.H / 4)) - 0.5 * (0.31 + 0.275)) < 1e-12
    g = np.arange(0.0, 4.01, 1.0)
    N = np.where((g[:, None] <= 2) & (g[None, :] <= 2), 4, np.where((g[:, None] <= 3) & (g[None, :] <= 3), 3, 1))
    a3, a4 = TN.areas(N, g)
    assert abs(a4 - 2.5 ** 2) < 1e-12 and abs(a3 - 3.5 ** 2) < 1e-12


def test_power_diagram_and_aliasing_formulas():
    """W2: with the kernel at (H, 0) of mass βε the entry times tie at H/2 − σ²ln β/H (the samples fine, so the continuum's value);
    W3: on samples of spacing ds along x₂ the face moves as the Poisson sum says (harmonics q^{k²}, q = e^{−π²σ²/ds²})."""
    for beta in (1.01, 1.5, 2.0):
        assert abs(TN.power_face(beta) - (TN.H / 2 - TN.SIG ** 2 * np.log(beta) / TN.H)) < 1e-9
    for ds in (0.05, 0.2, 0.3):
        r = TN.w_alias(ds)
        assert abs(r['shift'] - r['shift_p']) < 1e-9
    assert abs(TN.w_alias(0.05)['shift']) < 1e-12 and abs(TN.w_alias(0.3)['shift']) > 0.1


def test_x_theta_falls_along_the_face_and_its_average_with_the_dimension():
    """x_θ(t) from a table, and the face averages over the (n − 1)-cube of offsets decreasing with n."""
    tt = np.array([0.0, TN.H / 4, TN.H / 2, TN.H])
    xt = np.array([6.3, 6.2, 6.0, 5.3])
    fa = [TN.face_average(tt, xt, n) for n in (2, 3, 4, 5)]
    assert all(a > b for a, b in zip(fa, fa[1:])) and 6.0 < fa[0] < 6.3


def test_product_rule_at_unequal_gaps_and_offsets():
    """Round 47 (D3at-12): with x₁ ≠ x₂ the offsets t₁, t₂ belong to their own edges — swapping them changes the count; the offset
    enters by |t|; the diagonal's share decides between 3 and 4."""
    T = dict(x=np.array([0.0, 2.0, 4.0, 8.0]), t=np.array([0.0, TN.H / 2]),
             pt=np.array([[0.5, 0.2, 0.04, 0.0], [0.5, 0.1, 0.005, 0.0]]))
    assert TN.product_count(T, 2.0, 4.0, 0.0, TN.H / 2) == 2         # p = 0.2, 0.005: the second edge's share below 1%
    assert TN.product_count(T, 2.0, 4.0, TN.H / 2, 0.0) == 3         # p = 0.1, 0.04: both edges, the diagonal 0.004 not
    assert TN.product_count(T, 2.0, 4.0, 0.0, -TN.H / 2) == 2        # the offset's sign does not matter
    assert TN.product_count(T, 2.0, 2.0, 0.0, 0.0) == 4              # the diagonal 0.2² = 0.04 above 1% …
    assert TN.product_count(T, 2.0, 4.0, 0.0, 0.0) == 3              # … 0.2·0.04 = 0.008 not


def test_x_theta_t_and_its_face_averages():
    """x_θ(t) is each offset's crossing of TH_REL (exact on rows linear in x); x̄_θ the trapezoid rule over the offsets (not their
    mean); the face average over the (n − 1)-cube of a quadratic x_θ = a − b|t|² is a − b(n − 1)H²(1 − 1/m²)/12 (the midpoint rule
    with the Euclidean |t|)."""
    xc = np.array([6.3, 6.2, 5.5])
    tt = np.array([0.0, TN.H / 8, TN.H / 2])
    x = np.linspace(0.0, 10.0, 11)
    T = dict(x=x, t=tt, pt=np.array([TN.TH_REL * (1 + 0.1 * (c - x)) for c in xc]))
    assert np.abs(TN.x_theta_t(T) - xc).max() < 1e-12
    assert abs(TN.x_theta_avg(T) - 5.95) < 1e-12                       # the mean of the three would be 6.0
    a, b, m = 6.3, 20.0, 24
    tq = np.linspace(0.0, TN.H, 2001)
    for n in (2, 3, 4, 5):
        exact = a - b * (n - 1) * TN.H ** 2 * (1 - 1 / m ** 2) / 12
        assert abs(TN.face_average(tq, a - b * tq ** 2, n, m) - exact) < 1e-7
    assert TN.face_average(tq, a - b * tq ** 2, 1) == a


def test_the_asymptotic_probabilities_and_the_vertex_areas():
    """P(N ≥ 2) = 1 − (1 − 2w/H)² (the complement of the cell's inner square, w = 2σ²x̄_θε^{1/2}/H); E[N − 1] = 2n·w/H on hℤⁿ;
    P(N ≥ 3) = 4(2σ²ε^{1/2}/H)²a₃/H² (four vertices); B4's counts symmetrised before the areas, and report() uses them."""
    xbar, eps = 6.0, 1e-6
    w = 2 * TN.SIG ** 2 * xbar * np.sqrt(eps) / TN.H
    assert abs(TN.p_layers(xbar, eps) - (1 - (1 - 2 * w / TN.H) ** 2)) < 1e-15
    for n in (2, 3, 4, 5):
        assert abs(TN.e_extra(xbar, n) * np.sqrt(eps) - 2 * n * w / TN.H) < 1e-15
    assert abs(TN.p_vertices(2.25) * eps - 4 * (2 * TN.SIG ** 2 * np.sqrt(eps) / TN.H) ** 2 * 2.25 / TN.H ** 2) < 1e-15
    rows = [dict(x=(0.0, 0.0), n=4, stationary=True), dict(x=(0.0, 1.0), n=3, stationary=True), dict(x=(0.0, 2.0), n=1, stationary=True),
            dict(x=(1.0, 1.0), n=3, stationary=True), dict(x=(1.0, 2.0), n=1, stationary=True), dict(x=(2.0, 2.0), n=1, stationary=True)]
    B4 = dict(grid=np.array([0.0, 1.0, 2.0]), rows=rows)
    g, N = TN.vertex_counts(B4)
    assert (N == N.T).all() and N[1, 0] == 3
    a3, a4 = TN.areas(N, g)
    assert abs(a3 - 2.25) < 1e-15 and abs(a4 - 0.25) < 1e-15             # unsymmetrised: 1.75
    assert f"P(N ≥ 3) = {TN.p_vertices(2.25):.0f}·ε and P(N = 4) = {TN.p_vertices(0.25):.0f}·ε" in TN.report(B4=B4)


def test_the_bisection_and_the_default_mass():
    """`crossing` finds where a falling share crosses TH_REL (T2, x_θ); `margins_diag`'s default ε is B's and B4's 10⁻¹⁰."""
    assert abs(TN.crossing(lambda x: 0.5 * np.exp(-x), 0.0, 14.0, 50) - np.log(50)) < 1e-12
    K, c0 = TN.local()
    assert TN.margins_diag(K, c0, 1.0, 3.0) == TN.margins_diag(K, c0, 1.0, 3.0, 1e-10)


def test_formed_kernels_are_booked_by_their_final_feature():
    """C2's bookkeeping: a formed kernel that ends at another feature is captured and counted there; the formed test uses the
    start's feature's amplitude."""
    c0 = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    cs, As = np.array([[0.1, 0.0], [2.9, 0.0]]), np.array([1.0, -0.5])
    c = np.array([[0.1, 0.0], [2.8, 0.0], [2.0, 0.0], [2.9, 0.0]])
    w = np.array([0.9, -0.2, 0.004, -0.3])                             # kernel 2: 0.004 < 1% of 0.5, not formed
    bk = TN.formed_by_feature(c0, cs, As, w, c)
    assert list(bk['near']) == [0, 0, 1, 1] and list(bk['fin']) == [0, 1, 1, 1] and list(bk['big']) == [True, True, False, True]
    assert bk['n_start'] == [2, 1] and bk['n_fin'] == [1, 2] and bk['captured'] == [(1, 0, 1)]


def test_inner_k_says_why_it_stopped():
    """Round 47 (D3at-9): 'converged' (an edge), 'cap' (a vertex's slow tail: the drift grows with the doubled chunk, so a drift
    within share_tol is not convergence; the remaining change from the doubling sequence), 'clock'; each drift is the change over its
    own doubling (against separate runs of 0 … 3 doublings)."""
    L = X.law()
    K, _ = TN.local()
    pts, T = TN.TIES['edge']
    f = TN.A_B * K.phi(np.atleast_2d(T + np.array([0.0, TN.H / 2])))[:, 0]
    r = RN.inner_k(L, K, f, np.array(pts, float), np.array([0.0, 4.0]), share_tol=1e-7)
    assert r['stop'] == 'converged' and r['stationary'] and r['drift'] <= 1e-7
    r = RN.inner_k(L, K, f, np.array(pts, float), np.array([0.0, 4.0]), share_tol=1e-7, max_secs=0.0)
    assert r['stop'] == 'clock' and not r['stationary'] and r['s'] == 0.0 and r['drifts'] == []
    pts, T = TN.TIES['vertex']
    f = TN.A_B * K.phi(np.atleast_2d(T))[:, 0]
    x = np.array([0.0, 3.0, 10.0, 13.0])                               # one of B4's unsettled problems
    rs = [RN.inner_k(L, K, f, np.array(pts, float), x, tol=1.0, share_tol=0.0, max_doublings=k) for k in range(4)]
    assert all(q['stop'] == 'cap' and len(q['drifts']) == k for k, q in enumerate(rs))
    d = rs[-1]['drifts']
    assert all(d[j - 1] == np.abs(rs[j]['shares'] - rs[j - 1]['shares']).max() for j in (1, 2, 3))
    assert d[0] < d[1] < d[2] and rs[-1]['remaining'] == np.inf          # a creep, not a geometric approach
    q = RN.inner_k(L, K, f, np.array(pts, float), x, tol=1.0, share_tol=d[0] / 2, max_doublings=3)
    assert q['stop'] == 'cap'                                          # share_tol below every drift
    q = RN.inner_k(L, K, f, np.array(pts, float), x, share_tol=1e-7)   # B4's settings: the cap after the drifts turn
    dq = q['drifts']
    assert q['stop'] == 'cap' and len(dq) == 12 and dq[-1] < dq[-2]
    assert abs(q['remaining'] - dq[-1] ** 2 / (dq[-2] - dq[-1])) <= 1e-12 * q['remaining'] and q['remaining'] < 1e-5


def test_c3_outcomes():
    """C3's classification by the final counts: a capture tie before a capture, a capture before the count, two kernels at feature 0
    bound whatever the residual at the end (the finite mass's relaxation); with one kernel per feature and a residual, the runner-up's
    recorded distance decides (round 53): not departed (within 2σ throughout: C4's d_held, C2's d_max), stranded (it left), unsettled
    (no distance recorded)."""
    row = lambda nf, cap, rms, **kw: dict(n_fin=nf, captured=cap, rms=rms, **kw)
    assert TN.c3_outcome(row([1, 1, 2], [(163, 0, 2)], 1e-11)) == 'capture tie'
    assert TN.c3_outcome(row([1, 1, 1], [(163, 0, 2)], 1e-12)) == 'captured'
    assert TN.c3_outcome(row([2, 1, 1], [], 4e-8)) == 'bound'
    assert TN.c3_outcome(row([1, 1, 1], [], 3e-6)) == 'unsettled'
    assert TN.c3_outcome(row([1, 1, 1], [], 3e-6, d_held=0.44)) == 'not departed'
    assert TN.c3_outcome(row([1, 1, 1], [], 3e-6, d_max=0.39)) == 'not departed'
    assert TN.c3_outcome(row([1, 1, 1], [], 3e-6, d_held=2.5)) == 'stranded'
    assert TN.c3_outcome(row([1, 1, 1], [], 1e-12, d_held=0.4)) == 'one kernel'


def test_inner_job_reports_the_stop():
    """`inner_job` passes `inner_k`'s stop reason and the remaining-change estimate through (B4u records them: referee round 47)."""
    r = TN.inner_job(('edge', TN.A_B, (0.0, 4.0), TN.H / 2))
    assert r['stop'] == 'converged' and r['stationary'] and 'remaining' in r


def test_v_law_is_the_static_rest_as_the_unfitted_feature_vanishes():
    """(v)'s law with the pull projected off the entrant's φ (referee round 49) is exact to leading order: as A′ → 0 the static rest
    point's displacement over the law tends to 1 (d = 6σ, w = 0.3: within 0.5%); the unprojected pull's (d/σ)² − 2 misses by 2.9%."""
    from scipy.optimize import least_squares
    dd, w, Ap = 6.0, 0.3, 1e-3
    res, e = TN.v_residual(1, dd, Ap=Ap)
    p0 = np.r_[TN.A_B - w, 0.0, 1e-3 * TN.SIG]
    q = least_squares(res, p0, args=(w,), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000).x
    along = float(q[2] * e[0]) / TN.SIG
    law = TN.v_law(dd, w, q[0], Ap=Ap)
    assert abs(along / law - 1) < 5e-3
    old = float(np.sqrt(Ap * (dd ** 2 - 2) * np.exp(-dd ** 2 / 4) * q[0] / (w * (q[0] + w))))
    assert abs(along / old - 1) > 0.02


def test_c3_outcome_at_its_residual_cut_and_formed_by_feature():
    row = lambda nf, cap, rms: dict(n_fin=nf, captured=cap, rms=rms)
    assert TN.c3_outcome(row([1, 1, 1], [], 1.01e-8)) == 'unsettled' and TN.c3_outcome(row([1, 1, 1], [], 0.99e-8)) == 'one kernel'
    c0 = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    cs = np.array([[0.0, 0.0], [1.0, 0.0]])
    As = np.array([1.0, 0.5])
    w = np.array([0.5, 0.004, 1e-4])
    bk = TN.formed_by_feature(c0, cs, As, w, c0)                       # the second kernel's 0.004 is below θ of feature 0's 1 …
    assert bk['n_start'] == [1, 0] and bk['captured'] == []
    w2 = np.array([0.5, 0.2, 0.1])                                     # … the third kernel, started nearest feature 1, ends at feature 0
    bk = TN.formed_by_feature(c0, cs, As, w2, np.array([[0.0, 0.0], [1.0, 0.0], [0.1, 0.0]]))
    assert bk['n_start'] == [1, 2] and bk['n_fin'] == [2, 1] and bk['captured'] == [(2, 1, 0)]


def test_b4u_note_with_the_continuation():
    """B4's caption on hand-built reruns: the stops' reasons, the counts, the shares' changes and the continuation's (B4x)."""
    sh = lambda a: np.array([1 - a, a, 0.0, 0.0])
    B4 = dict(rows=[dict(x=(5.0, 9.0), n=2, shares=sh(0.02), stationary=False), dict(x=(9.0, 9.0), n=1, shares=sh(0.009), stationary=False)])
    B4u = dict(rows=[dict(x=(5.0, 9.0), n=2, shares=sh(0.020001), stop='cap', drift=2e-7, remaining=np.inf),
                     dict(x=(9.0, 9.0), n=1, shares=sh(0.009), stop='converged', drift=1e-11, remaining=np.inf)])
    B4x = dict(rows=[dict(x=(5.0, 9.0), n=2, shares=sh(0.0200025), stop='converged', drifts=[1e-9, 4e-8, 1e-7, 3e-8, 1e-9])])
    s = TN.b4u_note(B4, B4u, B4x)
    assert "1 at the doubling cap (the last drift 2.0e-07)" in s and "1 settled on rerun" in s
    assert "their counts unchanged at 2 of 2, their shares moved by at most 1.0e-06" in s
    assert "the counts unchanged at 1 of 1, the shares moved after the first cap by at most 1.5e-06" in s
    assert "the drifts peaking at doublings 3–3, 1 settled, 0 at the new cap, 0 at the clock" in s
    assert 'doublings to reach' not in s and 'no count changes for' not in s


def test_the_fold_constant_and_its_level_in_the_report():
    """V3's c = w_f/ℓ(d), ℓ(d) = A′(d/σ)²e^{−d²/(4σ²)}, and `fold_level`'s c(d)ℓ(d) from V's folds (interpolated in d)."""
    V = [dict(n=1, d=5.0, Ap=TN.AP_V, rows=[dict(w=0.3, we=0.5, along=0.1, across=0.0)], fold=(0.1424, 0.1424), along_fold=0.94),
         dict(n=1, d=6.0, Ap=TN.AP_V, rows=[dict(w=0.3, we=0.5, along=0.08, across=0.0)], fold=(0.02066, 0.02066), along_fold=0.86)]
    lev = lambda a, d: a * d ** 2 * np.exp(-d ** 2 / 4)
    c5, c6 = 0.1424 / lev(TN.AP_V, 5.0), 0.02066 / lev(TN.AP_V, 6.0)
    assert abs(TN.fold_level(5.5, 0.7, V) - 0.5 * (c5 + c6) * lev(0.7, 5.5)) < 1e-12
    assert TN.fold_level(5.5, 0.7, None) is None
    V3 = [dict(n=1, d=6.0, Ap=1.4, rows=[], fold=(0.0403, 0.0403), along_fold=0.85)]
    md = TN.report(V=V, V3=V3)
    row = next(l for l in md.splitlines() if l.startswith('| 6 |') and 'A′' not in l and l.count('|') == len(TN.AP_V3) + 2)
    assert f"{0.0403 / lev(1.4, 6.0):.3f}" in row


def test_v_job_off_its_branch_and_the_rests_energy():
    """`v_job` (rounds 49, 53): at A′/A = 3.5 and 5σ the tail is unbound at the grid's top; continued from 0.99A the cold start's small
    entrant would leave for A′, a fresh start (the entrant moved along −e) finds the branch, which folds at c = 3.135; on a branch the
    rests' energy E*(w) increases with w at every step (the hold), as the report states, and on the far side it decreases at every step
    with no fold; C4's report counts the two tolerances' agreement and names the disagreement."""
    off = TN.v_job((1, 5.0, 2.8))
    assert abs(off['top'] - 0.99 * TN.A_B) < 1e-12 and off['start'].startswith('fresh') and not off['entrant_left']
    assert off['fold'] is not None and abs(np.sqrt(off['fold'][0] * off['fold'][1]) / (2.8 * 25 * np.exp(-25 / 4)) - 3.135) < 2e-3
    on = TN.v_job((1, 6.0, TN.AP_V))
    assert not on['entrant_left'] and on['top'] == TN.W_V[0] and on['fold'] is not None and on['start'] == 'cold'
    E, w = [q['E'] for q in on['rows']], [q['w'] for q in on['rows']]
    assert all((E[i + 1] - E[i]) / (w[i + 1] - w[i]) > 0 for i in range(len(E) - 1))
    md = TN.report(V=[on], V3=[off])
    assert '3.135 (fresh start)' in md and 'increases with w at every step' in md
    far = TN.v_job((1, 6.0, TN.AP_V, -1))
    Ef, wf = [q['E'] for q in far['rows']], [q['w'] for q in far['rows']]
    assert far['fold'] is None and len(far['rows']) == len(TN.W_V) and all(q['along'] < 0 for q in far['rows'])
    assert all((Ef[i + 1] - Ef[i]) / (wf[i + 1] - wf[i]) < 0 for i in range(len(Ef) - 1))
    md = TN.report(V=[on], VF=[far])
    assert f'| 1 | 6 | ' in md and f'none down to {TN.W_V[-1]:g} | {len(Ef) - 1} of {len(Ef) - 1},' in md
    base = dict(x=4.9, eps=1e-11, d2=6.08, A2=0.7, s_dep=None, s2=3.07e6, w_dep=None, w_max=0.0232, n_fin=[2, 1, 1], rms=5e-9, d_end=0.002,
                captured=[], w_end=0.0232)
    other = dict(base, s_dep=3.3e5, w_dep=0.0165, w_max=0.0215, n_fin=[1, 1, 1], rms=1e-16, d_end=6.08, captured=[(5, 0, 2)], w_end=0.7)
    md = TN.report(C4=[base], C4b=[other])
    assert 'the same outcome at 0 of 1 points; they differ at x = 4.9, ε = 1e-11: bound at 1e-12' in md
    md = TN.report(C4=[other], C4b=[other])
    assert 'the same outcome at 1 of 1 points.' in md


def test_t3_product_takes_the_complements():
    """T3's product rule: the corner 1_S gets Π_{i∈S}p_i Π_{i∉S}(1 − p_i), summing to 1 over the cube."""
    pv = [0.1, 0.3, 0.02]
    S = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    pr = TN.t3_product(pv, S)
    assert abs(pr.sum() - 1.0) < 1e-15
    assert abs(pr[S.index((0, 0, 0))] - 0.9 * 0.7 * 0.98) < 1e-15 and abs(pr[S.index((1, 1, 0))] - 0.1 * 0.3 * 0.98) < 1e-15


def test_the_t3_report_prints_the_frozen_tails():
    """T3's table (round 53): the unformed corners with shares ≥ 5·10⁻⁴ get the product's relative error at H/√2."""
    pv = {1.0: 0.3, 3.0: 0.05}
    S = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    x = (1.0, 1.0, 3.0)
    sh = TN.t3_product([pv[v] for v in x], S)
    sh[S.index((1, 1, 1))] *= 1.25                                           # 0.0045 → the product 20% below
    md = TN.report(T3=dict(edge=pv, cube=[dict(x=x, S=S, shares=sh.tolist(), stop='converged', drift=0.0)]))
    assert '| (1, 1, 3) | 7 (converged) | 7, ' in md and '| -20.0% (1) | — |' in md


def test_d2_report_shift_against_x_theta_at_the_offset():
    """D2's shift is the loop's crossing minus x_θ at its offset |t| (not at t = 0)."""
    xc = np.array([6.3, 6.2, 5.5])
    tt = np.array([0.0, TN.H / 8, TN.H / 2])
    x = np.linspace(0.0, 10.0, 11)
    T = dict(x=x, t=tt, pt=np.array([TN.TH_REL * (1 + 0.1 * (c - x)) for c in xc]))
    D2 = [dict(eps=1e-6, t=TN.H / 2, x=5.63), dict(eps=1e-6, t=0.0, x=6.43)]
    md = TN.report(T=T, D2=D2)
    assert '| 1e-06 | 0.500 | 5.630 | 5.500 | +0.130 |' in md and '| 1e-06 | 0.000 | 6.430 | 6.300 | +0.130 |' in md


def test_the_c3_fold_column_from_v():
    """C3's fold column for the 5σ variant: c(d)·A₂·(d/σ)²e^{−d²/(4σ²)} at feature 2's distance (V's c(d) interpolated)."""
    V = [dict(n=1, d=4.5, Ap=TN.AP_V, rows=[dict(w=0.3, we=0.5, along=0.1, across=0.0)], fold=(0.3, 0.3), along_fold=0.9),
         dict(n=1, d=5.5, Ap=TN.AP_V, rows=[dict(w=0.3, we=0.5, along=0.1, across=0.0)], fold=(0.08, 0.08), along_fold=0.9)]
    cs, As = TN.targets_C2()['near tie 5σ']
    d2 = float(np.linalg.norm(cs[2] - cs[0]) / TN.SIG)
    lev = lambda a, d: a * d ** 2 * np.exp(-d ** 2 / 4)
    c_ = np.interp(d2, [4.5, 5.5], [0.3 / lev(TN.AP_V, 4.5), 0.08 / lev(TN.AP_V, 5.5)])
    row = dict(name='near tie 5σ', eps=1e-9, mu=TN.MU_C2, w_tail=0.03, n_fin=[1, 1, 1], captured=[], rms=1e-12)
    md = TN.report(V=V, C3=[row])
    assert f"{c_ * As[2] * lev(1.0, d2):.3g} (d = {d2:.2f}σ)" in md


def test_the_secant_line_covers_v3():
    """The hold's sentence reads V's and V3's branches together: a V3 branch whose energy falls with w turns it."""
    up = dict(n=1, d=6.0, Ap=TN.AP_V, rows=[dict(w=0.3, we=0.5, along=0.1, across=0.0, E=2.0), dict(w=0.1, we=0.7, along=0.2, across=0.0, E=1.0)],
              fold=None, along_fold=None)
    down = dict(up, Ap=1.4, rows=[dict(w=0.3, we=0.5, along=0.1, across=0.0, E=1.0), dict(w=0.1, we=0.7, along=0.2, across=0.0, E=2.0)])
    assert 'increases with w at every step' in TN.report(V=[up])
    assert 'does NOT increase with w at every step' in TN.report(V=[up], V3=[down])
    assert '(1 steps;' in TN.report(V=[up], V3=[dict(up)]) and '(2 steps;' in TN.report(V=[up], V3=[dict(up, Ap=1.4)])   # each branch once


def test_the_race_job_both_sides():
    """R's ε-free race (round 53): at the gap 7 the runner-up, a tail (below θ|A| after the layer), goes toward feature 2 on the near
    tie and departs; on the mirrored target it goes to the far side, stays within 2σ and its coefficient grows."""
    r = TN.race_job(('near tie', 7.0, 1e-10, 3e5))
    assert r['side'] == 'toward' and r['s_d'] is not None and 5e4 < r['s_d'] < 2e5 and r['w_layer'] < TN.TH_REL * TN.A_B
    assert r['d_dep'] > 2.0 >= r['d_pre'] and abs(r['s_layer'] / 1e3 - 1) < 0.05            # the departure at 2σ; the layer's end at 10³
    assert r['s_d_grid'] / 1.04 < r['s_d'] <= r['s_d_grid'] and r['s_theta'] is None        # bisected within the grid step (round 58)
    assert r['e_closest'] is not None and r['e_closest'] > 0                                 # the collapse keeps the layer's side (a tail)
    m = TN.race_job(('near tie mirrored', 7.0, 1e-10, 3e5))
    assert m['side'] == 'far' and m['s_d'] is None and m['d_end'] < 2.0 and m['w_end'] > m['w_layer'] > 0
    assert m['closest'] < 1e-6 and r['closest'] < 1e-6                     # the inner layer's collapse
    assert r['e_end'] > TN.SIDE_R and m['e_end'] < -TN.SIDE_R               # the side, read again at the end


def test_the_race_side_of_a_formed_runner_up():
    """A formed runner-up (x = 2.5: its coefficient stays above θ|A|) rests about 0.13σ from feature 0, inside the old 0.2σ side
    threshold: the side is read at SIDE_R and agrees with its displacement along e at the end."""
    r = TN.race_job(('near tie', 2.5, 1e-9, 1e5))
    assert TN.SIDE_R < r['d_end'] < 0.2 and r['s_d'] is None and r['w_end'] > TN.TH_REL * TN.A_B
    assert r['side'] in ('toward', 'far') and (r['e_end'] > 0) == (r['side'] == 'toward') and abs(r['e_end']) > 0.95 * r['d_end']


def test_the_loop_records_the_side():
    """c4_job (round 53): at the gap 7 (a tail: the inner layer's collapse stops at ~5·10⁻⁹σ, above the rounding floor, so the side is
    the geometry's) the runner-up goes toward feature 2 on C2's target and away on the mirrored one; at ε = 10⁻⁶ neither departs by s₂."""
    for name, side in (('near tie', 'toward'), ('near tie mirrored', 'far')):
        r = TN.c4_job((7.0, 1e-6, 1e-10, name))
        assert r['side'] == side and r['closest'] > 1e-10 and TN.c3_outcome(r) == 'not departed' and r['s_dep'] is None
        assert abs(r['s2'] * np.sqrt(r['eps']) - (r['tau2'] - r['tau_e'])) < 1e-12 and TN.SIDE_R < r['d_held'] < 2.0


def test_race_rule_and_the_r_report():
    """The rule captured iff s_d < s₂, ties counted apart; the mirror and the scan tables; the FC line."""
    assert TN.race_rule(1e5, 3e5) == 'captured' and TN.race_rule(4e5, 3e5) == 'not departed' and TN.race_rule(None, 3e5) == 'not departed'
    mk = lambda x, side, sd, **kw: dict(dict(name='near tie', x=x, rtol=1e-13, s_end=3e8, status=0, side=side, closest=1e-9, s_d=sd,
                                             w_dep=0.01 if sd else None, w_layer=0.02, marks=[1e4, 1e5], w_marks=[0.02, 0.02],
                                             d_marks=[0.5, 0.6], s_theta=None, d_end=6.08, w_end=0.7), **kw)
    R = dict(main=[mk(4.9, 'toward', 3.3e5), mk(7.0, 'toward', 9.9e4)],
             mirror=[mk(7.0, 'far', None, name='near tie mirrored', s_theta=5.6e6, d_end=0.3, w_end=0.15)],
             scan=[mk(4.9, 'toward', None, rtol=1e-13), mk(4.9, 'far', None, rtol=1e-12)])
    c4 = lambda x, eps, out: dict(x=x, eps=eps, s2=9.7 / np.sqrt(eps), n_fin=[1, 1, 1] if out != 'bound' else [2, 1, 1],
                                  captured=[(5, 0, 2)] if out == 'captured' else [], rms=1e-16 if out == 'captured' else 5e-9, d_held=0.3)
    C4 = [c4(4.9, 1e-12, 'captured'), c4(7.0, 1e-8, 'bound'), c4(4.9, 1e-8, 'bound'),
          dict(c4(7.0, 1e-9, 'captured'), n_fin=[1, 1, 2])]                                # a capture tie: counted apart
    FC = [dict(name='near tie', side=-1, fold=(0.005, 0.005), at=(0.59, -0.58), w_end=0.003),
          dict(name='near tie', side=1, fold=(0.01668, 0.01668), at=(0.842, 0.840), w_end=0.003),
          dict(name='alone', side=1, fold=None, at=(0.59, -0.58), w_end=0.003)]
    md = '\n'.join(TN.r_lines(R, None, C4, None, FC))
    assert "against the loops' outcomes: 3 of 3 agree, 1 capture ties" in md
    assert '| C4 | 7 | 1e-08 | 9.7e+04 | 9.9e+04 | not departed | bound |' in md
    assert '| 7 | far | 1.0e-09 | 0.02 | 0.02, 0.02 | 0.300 | 5.6e+06 | none |' in md
    assert 'Runs not toward, per tolerance: 1e-12: 1 of 1; 1e-13: 0 of 1.' in md
    assert 'toward feature 2: the fold at w = 0.01668' in md and 'toward, feature 1 removed: no fold down to w = 0.003' in md
    assert abs(TN.fc_fold(FC) - 0.01668) < 1e-15 and TN.fc_fold(None) is None


def test_inner_k_clock_after_a_chunk_and_a_zero_drift(monkeypatch):
    """inner_k's clock read after a completed chunk (not only inside the integration), and the remaining change when the shares stop
    moving (one kernel: every drift zero → 0)."""
    L = X.law()
    K, _ = TN.local()
    pts, T = TN.TIES['edge']
    f = TN.A_B * K.phi(np.atleast_2d(T + np.array([0.0, TN.H / 2])))[:, 0]
    clock = [0.0]
    monkeypatch.setattr(RN.time, 'time', lambda: clock[0])
    real = RN.solve_ivp

    def chunk(*a, **k):
        out = real(*a, **k)
        clock[0] += 10.0
        return out
    monkeypatch.setattr(RN, 'solve_ivp', chunk)
    r = RN.inner_k(L, K, f, np.array(pts, float), np.array([0.0, 4.0]), tol=0.0, max_secs=5.0)
    assert r['stop'] == 'clock' and r['s'] > 0 and r['drifts'] == []
    monkeypatch.setattr(RN, 'solve_ivp', real)
    r1 = RN.inner_k(L, K, f, np.array(pts[:1], float), np.array([0.0]), tol=0.0, max_doublings=2)
    assert r1['stop'] == 'cap' and r1['drifts'] == [0.0, 0.0] and r1['remaining'] == 0.0


def test_the_third_moment_law_of_the_hold():
    """e3_law: odd in the side, ∝ w^{−1/2} at small w (E*₃·w^{1/2} → ∓(1/6)P₃(M₂*)^{3/2}), increasing in w toward A′ and decreasing away;
    V's measured secant tends to it as the displacement vanishes (d = 7σ: within 6% at w ≈ 0.28 on both sides)."""
    d = 6.0
    assert TN.e3_law(0.1, d, 1) == -TN.e3_law(0.1, d, -1)
    w = np.array([1e-6, 1e-7])
    l = TN.AP_V * d ** 2 * np.exp(-d ** 2 / 4)
    P3 = TN.AP_V * np.sqrt(np.pi) * d ** 3 * np.exp(-d ** 2 / 4) / (8 * TN.SIG ** 2)
    assert np.allclose(TN.e3_law(w, d, 1) * np.sqrt(w), -P3 * (TN.SIG ** 2 * l) ** 1.5 / 6, rtol=1e-5)
    assert TN.e3_law(0.2, d, 1) > TN.e3_law(0.1, d, 1) and TN.e3_law(0.2, d, -1) < TN.e3_law(0.1, d, -1)
    for side in (1, -1):
        r = TN.v_job((1, 7.0, TN.AP_V, side))
        q = TN.e3_ratios(r)
        assert q and abs(q[0][0] - 0.28) < 0.05 and abs(q[0][1] - 1) < 0.06 and np.sign(q[0][2]) == side


def test_the_report_reproduces_the_committed_md():
    """Golden (referee round 58): the report from the committed raw is the committed md, byte for byte."""
    out = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
    old = E.load(os.path.join(out, 'transit_nd_raw.npz'))
    assert TN.report(*[old.get(k) for k in TN.KEYS]) == open(os.path.join(out, 'transit_nd.md')).read()


def test_the_mirror_race_grows_past_theta_and_a_formed_runner_up_does_not():
    """R's far side (round 58): at the gap 7 the held tail's coefficient passes θ|A| near s = 5.7·10⁶; a formed runner-up (x = 3) is above
    θ|A| after the layer, so it has no s_θ."""
    m = TN.race_job(('near tie mirrored', 7.0, TN.RTOL_R, 1e7))
    assert m['side'] == 'far' and m['s_theta'] is not None and 5e6 < m['s_theta'] < 6.5e6 and m['s_d'] is None
    f = TN.race_job(('near tie mirrored', 3.0, TN.RTOL_R, 1e5))
    assert f['w_layer'] > TN.TH_REL * TN.A_B and f['s_theta'] is None


def test_the_outcome_labels_of_held_runner_ups():
    """c3_outcome at d_held = 1.5σ and 2.5σ; C4M's held unformed tails labelled by their side whatever the residual (round 58)."""
    base = dict(n_fin=[1, 1, 1], captured=[], rms=1e-6)
    assert TN.c3_outcome(dict(base, d_held=1.5)) == 'not departed' and TN.c3_outcome(dict(base, d_held=2.5)) == 'stranded'
    for rms in (1e-12, 1e-6):
        assert TN.c4m_outcome(dict(base, rms=rms, d_held=0.3, w_end=0.005, side='far')) == 'held, unformed (far)'
    assert TN.c4m_outcome(dict(base, n_fin=[2, 1, 1], d_held=0.3, w_end=0.02, side='far')) == 'bound'


def test_pair_sums_are_the_critical_right_hand_side():
    """RP's `PairSums` (float64 and the long double) against `critical_rhs_nd` for the race's pair at a generic state."""
    L = TN.law()
    K, c0 = TN.full()
    cs, As = TN.targets_C2(TN.MU_R)['near tie']
    f = K.synth(cs, As)
    z = np.r_[3.0, 2.0, cs[0] + np.array([0.01, -0.02]), cs[0] + np.array([0.03, 0.01])]
    want = RN.critical_rhs_nd(L, K, f)(0.0, z)
    for dt in (np.float64, np.longdouble):
        assert np.allclose(TN.PairSums(L, K, cs, As, dt)(0.0, z), want, rtol=1e-11, atol=1e-16)


def test_the_pin_law_at_the_window_edge():
    """RP (round 58): on Part C's window the formed runner-up (x = 4.9) pins 1.55·10⁻¹⁴σ along −y, the law (the block's rest moving with W
    through the window's top edge, ∂c̄*/∂W = 1.14·10⁻¹²σ) within 3%, toward feature 2; with the edge 8.3σ away the pin and the shift vanish."""
    r = TN.rp_job(('near tie', 4.9, 1.0, 64))
    assert -1.7e-14 < r['D'][1][1] < -1.4e-14 and abs(r['law'][1][1] / r['D'][1][1] - 1) < 0.03 and r['side'] == 'toward'
    assert abs(r['dcdW'][1] / 1.14e-12 - 1) < 0.01 and abs(r['dcdW'][0]) < 1e-14 and r['s_d'] is not None
    e = TN.rp_job(('near tie', 4.9, 1.5, 64))
    assert np.linalg.norm(e['D'][1]) < 1e-15 and np.linalg.norm(e['dcdW']) < 1e-15
