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


def test_the_bisection_and_the_unbinding_level():
    """`crossing` finds where a falling share crosses TH_REL (T2, x_θ); `margins_diag`'s default ε is B's and B4's 10⁻¹⁰; C2's crude
    unbinding level A′(d/σ)e^{−d²/(4σ²)}, the largest over the later features."""
    assert abs(TN.crossing(lambda x: 0.5 * np.exp(-x), 0.0, 14.0, 50) - np.log(50)) < 1e-12
    K, c0 = TN.local()
    assert TN.margins_diag(K, c0, 1.0, 3.0) == TN.margins_diag(K, c0, 1.0, 3.0, 1e-10)
    cs = np.array([[0.0, 0.0], [6 * TN.SIG, 0.0], [0.0, 5 * TN.SIG]])
    lev = TN.w_unbind(cs, np.array([0.8, 0.7, -0.5]))
    assert abs(lev - max(0.7 * 6 * np.exp(-9.0), 0.5 * 5 * np.exp(-25 / 4))) < 1e-15


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
