"""
Proposition D‴.35 (PARADIGM.md §5d): the split of a class under the cohesion leak is an S_n-symmetric transcritical bifurcation.
On Proposition D‴.32's resolution instance (the lift of the three-atom approximant ϑ^b(q) of f_ε, the compact kernel at R = 1):
(i) along the splitting direction v the full particle drift is μξ_i + (c₃/2)(ξ_i² − S₂/n) up to third order (the cohesion, even, enters
at third order); (ii) the S_1 × S_{n−1} split state near the threshold is a saddle with the predicted eigenvalues; (iii) the particle
with the extreme sign(c₃)ξ_i(0) leaves first, on the side sign(c₃)v.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
import exp_split_bifurcation as SB                                        # noqa: E402
from src.core import cohesion_leak as CO                                  # noqa: E402

SQG = np.sqrt(E.G)


def _instance(counts):
    net, base, bq = E.approximant()
    w = np.asarray(counts) / np.sum(counts)
    at = E.gauge(net, bq, w)
    return net, base, at, w, np.repeat(at, counts, axis=0)


def test_the_reduced_equation_is_quadratic_with_the_residuals_third_derivative():
    counts = np.array([10, 10, 10])
    net, base, at, w, Th0 = _instance(counts)
    lam, v, c3s, _ = SB.coefficients(base, at, Th0, 1)
    assert abs(c3s[0] - c3s[1]) < 1e-5 * abs(c3s[0])
    sigma = 0.9 * lam / w[1]                                              # past the threshold: μ > 0
    mu = lam - sigma * w[1]
    loop = CO.CohesionLoop(base, 1.0, 'compact', sigma0=sigma)
    idx = np.arange(10, 20)
    xi_dir = np.random.default_rng(0).standard_normal(10); xi_dir -= xi_dir.mean(); xi_dir /= np.abs(xi_dir).max()
    errs = []
    for eps in (1e-3, 5e-4):
        xi = eps * xi_dir
        Th = Th0.copy(); Th[idx] += xi[:, None] * (SQG * v)[None, :]
        d = loop.drift(Th)[idx] / SQG[None, :] @ v                         # the drift's whitened component along v
        d = d - d.mean()
        pred = mu * xi + 0.5 * c3s[0] * (xi ** 2 - np.sum(xi ** 2) / 10)
        lin = mu * xi
        errs.append(np.abs(d - pred).max())
        assert np.abs(d - pred).max() < 0.1 * np.abs(pred - lin).max()  # the quadratic term is the residual's
    assert 5.0 < errs[0] / errs[1] < 11.0                                  # the remainder is cubic (ratio 8 for a halving)


def test_the_single_particle_split_state_is_a_saddle_with_the_predicted_eigenvalues():
    out = SB.b_job(5e-7)
    r = out['rows'][0]                                                     # q = 1
    n, mu = 10, out['mu']
    assert r['drift'] < 1e-15
    assert abs(r['groups'][0] / r['pred_groups'][0] - 1) < 0.01
    pos = np.sort(r['ev'])[1:]
    assert np.all(np.abs(pos / (mu * n / (n - 2)) - 1) < 0.01)            # the n − 2 unstable directions at +μn/(n − 2)
    assert abs(np.sort(r['ev'])[0] / (-mu) - 1) < 0.05                     # the radial one at −μ
    assert r['n_unstable'] == n - 2


def test_the_extreme_particle_leaves_first_on_the_side_of_the_quadratic_term():
    counts = np.array([4, 4, 4])
    net, base, at, w, Th0 = _instance(counts)
    lam, v, c3s, _ = SB.coefficients(base, at, Th0, 1)
    sigma = (lam - 4e-4) / w[1]
    rng = np.random.default_rng(5)
    Th = Th0 + 1e-3 * rng.standard_normal(Th0.shape) * SQG[None, :]
    idx = np.arange(4, 8)
    xi_of = lambda T: ((T[idx] - T[idx].mean(axis=0)) / SQG[None, :]) @ v
    first = []

    def cb(s_, T):
        x = xi_of(T)
        dev = np.abs(x - np.median(x))
        if not first and dev.max() > 0.1:
            i = int(np.argmax(dev)); first.append((i, np.sign(x[i] - np.median(x))))
    xi0 = xi_of(Th)
    CO.CohesionLoop(base, 1.0, 'compact', sigma0=sigma).run(Th, 30000.0, 0.25, rec_every=20, callback=cb)
    assert first, 'no particle left'
    assert first[0][0] == int(np.argmax(np.sign(c3s[0]) * xi0))
    assert first[0][1] == np.sign(c3s[0])
