"""Theorem D‴.18 (PARADIGM.md §5d): the constants of Theorem J are uniform in m in the per-particle norm |||·||| at a replicated configuration.
Checked on the averaged slow loop with the matching prior at the matching masses (the replicated configuration is the exact representation):
the second-order remainder of the drift is bounded by C_N|||u|||² with the same C_N for m = 6, 30, 150, 600; a separation perturbation of one
particle moves no other particle at first order and is the same at every m; and the averaged loop contracts |||u||| at the rate
min{λ_coll, λ_sep} with the same constant for m = 6, 30, 150."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import prior_design as PD

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846


def _setup(n_quad=96):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    exact = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    loop = PD.MixturePriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, exact, n_quad=n_quad)
    return net, exact, loop


def _norm3(u, lab, w):
    """|||u||| = max{(Σ_k w_k‖ū_k‖²_{G⁻¹})^{1/2}, max_i‖u_i − ū_{lab(i)}‖_{G⁻¹}}."""
    K = len(w)
    ubar = np.array([u[lab == k].mean(axis=0) for k in range(K)])
    coll = np.sqrt(np.sum(w * np.sum(ubar ** 2 / G[None, :], axis=1)))
    sep = np.sqrt(np.sum((u - ubar[lab]) ** 2 / G[None, :], axis=1)).max()
    return max(coll, sep)


def test_second_order_remainder_of_the_drift_is_uniform_in_m_in_the_per_particle_norm():
    net, exact, loop = _setup()
    w = np.ones(3) / 3
    rng = np.random.default_rng(0)
    CN, CNcoll = {}, {}
    h = 1e-5
    rem = lambda Th0, d0, u: loop.drift(Th0 + u) - d0 - (loop.drift(Th0 + h * u) - loop.drift(Th0 - h * u)) / (2 * h)   # second-order remainder
    ucoll = 0.2 * np.array([[0.6, 0.5, -0.3], [-0.4, 0.7, 0.2], [0.5, -0.2, 0.6]]) * SQG[None, :]                      # one displacement per class
    for m in (6, 30, 150, 600):
        lab = np.repeat(np.arange(3), m // 3)
        Th0 = exact[lab]
        d0 = loop.drift(Th0)
        assert np.abs(d0).max() < 1e-12                                     # the replicated configuration is stationary
        ratios = []
        for trial in range(4):
            u = rng.standard_normal((m, 3)) * SQG[None, :]
            u *= 0.2 / np.sqrt(np.sum(u ** 2 / G[None, :], axis=1)).max()   # every particle displaced by ≤ 0.2 (G⁻¹)
            N = rem(Th0, d0, u)
            ratios.append(np.sqrt(np.sum(N ** 2 / G[None, :], axis=1)).max() / _norm3(u, lab, w) ** 2)
        CN[m] = max(ratios)
        u = ucoll[lab]                                                       # a purely collective displacement: the K-atom remainder, lifted
        N = rem(Th0, d0, u)
        CNcoll[m] = np.sqrt(np.sum(N ** 2 / G[None, :], axis=1)).max() / _norm3(u, lab, w) ** 2
    vals = np.array(list(CN.values()))
    assert vals.max() < 0.1 and CN[600] <= CN[6]                             # ‖N_i(u)‖ ≤ C_N|||u|||² with one C_N for every m (not growing with m)
    vc = np.array(list(CNcoll.values()))
    assert vc.max() - vc.min() < 1e-6 * vc.max() and vc.max() < 0.1          # the collective remainder is exactly m-independent (0.064 here)


def test_a_separation_perturbation_moves_only_its_own_particle_at_first_order_and_is_the_same_at_every_m():
    net, exact, loop = _setup()
    resp = {}
    for m in (6, 60, 300):
        lab = np.repeat(np.arange(3), m // 3)
        Th0 = exact[lab]
        u = np.zeros((m, 3)); u[0] = np.array([0.1, 0.05, -0.05]); u[1] = -u[0]          # zero mean within class 0: a separation
        h = 1e-5
        lin = (loop.drift(Th0 + h * u) - loop.drift(Th0 - h * u)) / (2 * h)
        assert np.abs(lin[2:]).max() < 1e-9                                 # no other particle responds at first order
        assert np.abs(lin[0] + lin[1]).max() < 1e-9                          # the class mean does not move
        resp[m] = lin[0]
    assert np.abs(resp[6] - resp[60]).max() < 1e-9 and np.abs(resp[6] - resp[300]).max() < 1e-9
    # at zero residual the separation block is −σ₀I: the response is −σ₀u₀ (the fit's curvature vanishes with the residual)
    assert np.abs(resp[6] + SIGMA0 * np.array([0.1, 0.05, -0.05])).max() < 1e-8


def test_the_averaged_loop_contracts_the_per_particle_norm_at_an_m_independent_rate():
    net, exact, loop = _setup(n_quad=64)
    w = np.ones(3) / 3
    rng = np.random.default_rng(3)
    mu = 0.9 * SIGMA0                                                        # min{λ_coll, λ_sep} = σ₀ = 0.00707 here (λ_coll = 0.0074)
    consts = {}
    for m in (6, 30, 150):
        lab = np.repeat(np.arange(3), m // 3)
        u = rng.standard_normal((m, 3)) * SQG[None, :]
        u *= 0.5 / np.sqrt(np.sum(u ** 2 / G[None, :], axis=1)).max()
        Th = exact[lab] + u
        n0 = _norm3(u, lab, w)
        ds, worst = 1.0, 0.0
        for k in range(400):                                                 # midpoint rule in slow time, 400 slow units
            k1 = loop.drift(Th); Th = Th + ds * loop.drift(Th + 0.5 * ds * k1)
            worst = max(worst, _norm3(Th - exact[lab], lab, w) * np.exp(mu * (k + 1) * ds) / n0)
        consts[m] = worst
    vals = np.array(list(consts.values()))
    assert vals.max() < 1.5 and vals.max() < 1.3 * vals.min()               # |||u(s)||| ≤ C e^{−μs}|||u(0)||| with one C for every m
