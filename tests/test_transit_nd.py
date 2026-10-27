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
        # the exponential branch: the softmax, to O(ca/(m|A|))
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
    assert ri['stationary'] and np.abs(r['shares'] - ri['shares']).max() < 5e-4
    # not the softmax, not the edges' product
    assert np.abs(TN.softmax(np.r_[0.0, r['x_act'][1:]], L) - ri['shares']).max() > 1e-2
