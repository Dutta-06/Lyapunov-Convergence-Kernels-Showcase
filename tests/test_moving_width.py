"""Proposition D‴.59 (planned): the transit with moving widths (`src/core/moving_width.py`)."""
import os
import sys

import numpy as np

from src.core import moving_width as MW
from src.core import rich_transit_nd as RN

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))


def test_bumps_cw_derivatives_and_grams():
    """∂_cφ, ∂_σφ and ⟨r∇²_θφ⟩ against central differences of φ; the Gram blocks against their sums."""
    rng = np.random.default_rng(0)
    K = MW.BumpsCW(np.linspace(-1.5, 1.5, 201), 0.12)
    th = np.c_[rng.uniform(-0.5, 0.5, 4), rng.uniform(0.05, 0.3, 4)]
    r = rng.standard_normal(K.N)
    h = 1e-6
    for k in range(2):
        e = np.zeros(2); e[k] = h
        fd = (K.phi(th + e) - K.phi(th - e)) / (2 * h)
        assert np.abs(fd - K.dphi(th)[:, :, k]).max() < 1e-6 * np.abs(fd).max()
    g, hh, Q = K.corr(th, r, second=True)
    assert np.allclose(g, K.phi(th).T @ r / K.N, rtol=1e-12, atol=0)
    assert np.allclose(hh, np.einsum('s,sjk->jk', r, K.dphi(th)) / K.N, rtol=1e-12, atol=1e-15)
    for k in range(2):
        e = np.zeros(2); e[k] = h
        fd = (K.corr(th + e, r)[1] - K.corr(th - e, r)[1]) / (2 * h)
        assert np.abs(fd - Q[:, :, k]).max() < 1e-5 * np.abs(Q).max()
    P, D, H = K.grams(th)
    dP = K.dphi(th)
    assert np.allclose(P, K.phi(th).T @ K.phi(th) / K.N) and np.allclose(D[1, 2], K.phi(th)[:, 1] @ dP[:, 2, :] / K.N)
    assert np.allclose(H[0, 3], dP[:, 0, :].T @ dP[:, 3, :] / K.N)


def test_the_one_kernel_reduced_energy():
    """`overlap` against quadrature; the reduced energy vanishes at (c*, σ*) only and its closed-form gradient matches central differences."""
    x = np.linspace(-12, 12, 240001)
    dx = x[1] - x[0]
    q = np.exp(-(x - 0.3) ** 2 / (2 * 0.2 ** 2)) @ np.exp(-(x + 0.1) ** 2 / (2 * 0.45 ** 2)) * dx
    assert abs(MW.overlap(0.3, 0.2, -0.1, 0.45) / q - 1) < 1e-10
    A, cs, ss = 0.8, 0.03, 0.2
    assert abs(MW.reduced_energy(cs, ss, A, cs, ss)) < 1e-15
    rng = np.random.default_rng(1)
    for _ in range(10):
        c, s = rng.normal(0, 0.3), ss * np.exp(rng.normal(0, 0.7))
        assert MW.reduced_energy(c, s, A, cs, ss) > 0
        h = 1e-6
        fd = np.array([(MW.reduced_energy(c + h, s, A, cs, ss) - MW.reduced_energy(c - h, s, A, cs, ss)) / (2 * h),
                       (MW.reduced_energy(c, s + h, A, cs, ss) - MW.reduced_energy(c, s - h, A, cs, ss)) / (2 * h)])
        g = MW.one_kernel_grad(c, s, A, cs, ss)
        assert np.abs(fd - g).max() < 1e-6 * max(np.abs(g).max(), 1e-12)


def test_a_narrow_feature_is_formed_by_one_kernel_with_moving_widths():
    """At a quarter of the kernels' width the sliding path with moving widths enters one kernel and relaxes it to the feature exactly
    (centre, width, amplitude), within a small budget; with the widths fixed the same budget runs out."""
    import exp_moving_width as M
    L = M.law()
    ss = 0.25 * M.S0
    f = M.feature_sum([M.C_A], [ss], [M.A_A])
    Rm = RN.sliding_iss_nd(L, MW.BumpsCW(M.XS, M.S0), f, M.theta0(), np.ones(M.C0.size), max_events=40, exit_kkt=True, max_rhs=50_000)
    m = M.match(Rm, [M.C_A], [ss], [M.A_A])
    assert m['n'] == 1 and m['exact'] and not m['oob'] and m['entries'] == 1
    Rf = RN.sliding_iss_nd(L, RN.BumpsND(M.XS[:, None], M.S0), f, M.C0[:, None], np.ones(M.C0.size), max_events=40,
                           max_entries=M.MAX_ENT_FIXED, exit_kkt=True, max_rhs=50_000)
    assert Rf['out_of_budget'] and not M.match(Rf, [M.C_A], [ss], [M.A_A])['exact']
