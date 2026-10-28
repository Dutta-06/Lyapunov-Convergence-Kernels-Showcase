"""Proposition D‴.58: a feature wider than the kernels as a Gaussian quadrature of its mixing law (`src/core/wide_feature.py`)."""
from math import factorial

import numpy as np

from src.core import wide_feature as W

SIG = 0.12


def ss(r_):
    return SIG * np.sqrt(1 + r_ ** 2)


def test_half_squeeze_error_is_the_mehler_series():
    """(i): the half-squeeze rule's L² error (trapezoid) is the exact Mehler series, at τ/σ = 1/2, 1, 2 and K ≤ 6."""
    for r_ in (0.5, 1.0, 2.0):
        t = W.params(SIG, ss(r_))['t']
        for K in range(1, 7):
            e = W.l2_error(*W.half_squeeze(K, SIG, ss(r_)), SIG, ss(r_))
            s2, f2 = W.mehler_series(K, t)
            assert abs(e / np.sqrt(s2) - 1) < 1e-9
            assert 0.97 < f2 / s2 <= 1.0 and abs(f2 - W.first_term(K, t)) < 1e-12 * f2


def test_gh_defect_vanishes_below_2k_and_its_tail_ratios_are_integer_polynomials():
    """D_{n,K} = 0 for odd n and n < 2K, D_{2K,K} = 2^K K!, D_{2K+2,K} = −(K − 1)D_{2K,K}, D_{2K+4,K} = −(K² − 25K − 6)/2·D_{2K,K}."""
    for K in range(1, 13):
        assert all(W.gh_defect(n, K) == 0 for n in range(0, 2 * K)) and W.gh_defect(2 * K + 1, K) == 0
        D = W.gh_defect(2 * K, K)
        assert D == 2 ** K * factorial(K)
        assert W.gh_defect(2 * K + 2, K) == -(K - 1) * D
        assert 2 * W.gh_defect(2 * K + 4, K) == -(K * K - 25 * K - 6) * D


def test_half_squeeze_rule_is_the_gauss_hermite_rule_of_the_reweighted_law():
    """Nodes τ(σ/σ*)^{1/2}z_k, symmetric; amplitudes positive, A(σ*/σ)^{1/2} for K = 1; moment matching's nodes τz_k."""
    y, a = W.half_squeeze(1, SIG, ss(1.0), A=0.8)
    assert abs(y[0]) < 1e-15 and abs(a[0] - 0.8 * np.sqrt(ss(1.0) / SIG)) < 1e-14
    z, w = W.gh(5)
    y, a = W.half_squeeze(5, SIG, ss(2.0))
    tau = W.params(SIG, ss(2.0))['tau']
    assert np.allclose(y, tau * np.sqrt(SIG / ss(2.0)) * z, rtol=0, atol=1e-15) and np.allclose(y, -y[::-1], atol=1e-15)
    assert np.all(a > 0) and np.allclose(a, a[::-1], rtol=1e-13)
    assert np.allclose(W.moment_matching(5, SIG, ss(2.0))[0], tau * z, rtol=0, atol=1e-15)


def test_rates_per_kernel():
    """At τ/σ = 1 the half-squeeze rule's error falls by 5.64 at K = 8 (→ 1/t = 5.83, the O(1/K) deficit 1/(4Kt)), moment matching's
    by 1/λ = 3 to 0.2%."""
    s_ = ss(1.0)
    P = W.params(SIG, s_)
    e = [W.l2_error(*W.half_squeeze(K, SIG, s_), SIG, s_) for K in (7, 8)]
    m = [W.l2_error(*W.moment_matching(K, SIG, s_), SIG, s_) for K in (7, 8)]
    assert abs(e[0] / e[1] - (1 - 1 / (4 * 8)) / P['t']) < 0.01 / P['t']
    assert abs(m[0] / m[1] * P['lam'] - 1) < 2e-3


def test_optimal_fit_improves_slightly_on_the_rule():
    """The optimal 4-kernel fit at τ/σ = 1: 0.44% below the half-squeeze rule's error, its nodes within 0.5% of the rule's."""
    s_ = ss(1.0)
    y, a = W.half_squeeze(4, SIG, s_)
    e = W.l2_error(y, a, SIG, s_)
    yo, ao, eo = W.optimal_fit(y, SIG, s_)
    assert 0.99 * e < eo < 0.997 * e
    assert np.abs(np.sort(yo) / np.sort(y) - 1).max() < 5e-3
