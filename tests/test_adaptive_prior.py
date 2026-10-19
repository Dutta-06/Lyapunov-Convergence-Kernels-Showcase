"""Theorem D‴.23 (PARADIGM.md §5d): adaptive priors. The joint system is the gradient flow of 𝒥 = L̄ + (σ₀/2)W²(μ_Θ, ν_Z) in the stated metric
(the particles' drift −mG∇_ϑ𝒥 and the prior's −(η/(σ₀q_j))G∇_{z_j}𝒥, against finite differences); the exact representation with the prior at its
atoms is stationary, while the fixed prior's stationary point is not stationary for the adaptive prior; the linearisation at the exact
representation has the separation rate σ₀ and the collective pairs' rates λ_j^± of pair_rates; the fixed prior's bias is σ₀(H_L + σ₀I)⁻¹x⁰ to
first order; a gated run keeps the prior frozen while mass is in transit, with 𝒥 non-increasing; the joint linearisation is self-adjoint in
the gradient-flow metric (M times the Jacobian is minus the Hessian of 𝒥); and the settled gate opens only when the fixed-prior drift has
fallen below its tolerance, after the transit."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import ot_leak as OT
from src.core import adaptive_prior as AP

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846
Q = np.ones(3) / 3


def _setup(shift=0.15, n_quad=96):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    rep = lambda sh: net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + sh] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    exact, prior = rep(0.0), rep(shift)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=n_quad)
    return loop, exact, prior


def test_the_joint_system_is_the_gradient_flow_of_the_joint_energy():
    loop, exact, prior = _setup()
    eta = 3 * SIGMA0
    apl = AP.AdaptivePriorLoop(loop, eta)
    rng = np.random.default_rng(2)
    m = 12
    Th = np.repeat(exact, m // 3, axis=0) + 0.05 * rng.standard_normal((m, 3)) * np.sqrt(G)[None, :]
    Z = prior + 0.02 * rng.standard_normal((3, 3))
    gamma, _ = AP.plan_to(Th, Z, Q, G)
    dTh, dZ = apl.drift(Th, Z, gamma)
    J = lambda T_, Z_: apl.energy(T_, Z_)
    h = 1e-6
    gT = np.zeros_like(Th); gZ = np.zeros_like(Z)
    for i in range(m):
        for c in range(3):
            E = np.zeros_like(Th); E[i, c] = h
            gT[i, c] = (J(Th + E, Z) - J(Th - E, Z)) / (2 * h)
    for j in range(3):
        for c in range(3):
            E = np.zeros_like(Z); E[j, c] = h
            gZ[j, c] = (J(Th, Z + E) - J(Th, Z - E)) / (2 * h)
    assert np.abs(dTh + m * G[None, :] * gT).max() < 1e-7                          # dϑ_i/ds = −mG∇_{ϑ_i}𝒥
    assert np.abs(dZ + (eta / (SIGMA0 * Q[:, None])) * G[None, :] * gZ).max() < 1e-7  # dz_j/ds = −(η/(σ₀q_j))G∇_{z_j}𝒥


def test_the_exact_representation_is_stationary_and_the_fixed_priors_stationary_point_is_not():
    loop, exact, prior = _setup()
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    m = 15
    Th = np.repeat(exact, m // 3, axis=0)
    gamma, _ = AP.plan_to(Th, exact, Q, G)
    dTh, dZ = apl.drift(Th, exact, gamma)
    assert np.abs(dTh).max() < 1e-14 and np.abs(dZ).max() < 1e-14 and apl.energy(Th, exact) < 1e-28
    star, res, _ = loop.stationary_atoms(prior)                                  # the fixed prior's stationary point (biased)
    Th = np.repeat(star, m // 3, axis=0)
    gamma, _ = AP.plan_to(Th, prior, Q, G)
    dTh, dZ = apl.drift(Th, prior, gamma)
    assert res < 1e-12 and np.abs(dTh).max() < 1e-12                             # the particles are at rest …
    assert np.abs(dZ).max() > 1e-4                                               # … but the prior moves toward the class atoms


def test_the_linearisation_has_the_separation_rate_and_the_collective_pair_rates():
    loop, exact, prior = _setup()
    for eta in (SIGMA0, 10 * SIGMA0):
        apl = AP.AdaptivePriorLoop(loop, eta)
        m = 6
        Th0 = np.repeat(exact, m // 3, axis=0)
        gamma, _ = AP.plan_to(Th0, exact, Q, G)
        def F(v):
            T_ = v[:3 * m].reshape(m, 3); Z_ = v[3 * m:].reshape(3, 3)
            a, b = apl.drift(T_, Z_, gamma)
            return np.concatenate([a.ravel(), b.ravel()])
        v0 = np.concatenate([Th0.ravel(), exact.ravel()]); d = v0.size; h = 1e-6
        Jac = np.zeros((d, d))
        for j in range(d):
            e = np.zeros(d); e[j] = h
            Jac[:, j] = (F(v0 + e) - F(v0 - e)) / (2 * h)
        rates = np.sort(-np.linalg.eigvals(Jac).real)
        HL, hs, V = AP.gauss_newton_block(loop, exact, Q)
        lm, lp = AP.pair_rates(hs, SIGMA0, eta)
        pred = np.sort(np.concatenate([lm, lp, np.full(3 * (m - 3), SIGMA0)]))  # 9 pairs and the separations of three classes of two
        assert np.abs(rates - pred).max() < 1e-6 * max(1.0, pred.max())


def test_the_fixed_priors_bias_is_the_first_order_prediction():
    loop, exact, prior = _setup(shift=0.05)
    HL, hs, V = AP.gauss_newton_block(loop, exact, Q)
    star, res, _ = loop.stationary_atoms(prior)
    x0 = AP.whitened(prior - exact, Q, G); y = AP.whitened(star - exact, Q, G)
    pred = SIGMA0 * np.linalg.solve(HL + SIGMA0 * np.eye(9), x0)
    assert res < 1e-12 and np.linalg.norm(y - pred) < 0.05 * np.linalg.norm(y)   # first order in the prior's error (a smaller shift than the experiment's)
    assert np.linalg.norm(y) < 0.2 * np.linalg.norm(x0)                          # the stiff modes are suppressed …
    assert abs((V.T @ y)[0]) > 0.5 * abs((V.T @ x0)[0])                          # … the weakest mode keeps most of the error


def test_a_gated_run_freezes_the_prior_during_the_transit_and_decreases_the_energy():
    loop, exact, prior = _setup()
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    m = 15
    Th = np.repeat(exact, m // 3, axis=0)
    Th[0] = 0.4 * exact[0] + 0.6 * exact[1]                                      # a particle of class 0 inside cell 1: mass in transit
    gamma, _ = AP.plan_to(Th, prior, Q, G)
    assert apl.in_transit(Th, prior, gamma)
    rec = []
    apl.run(Th, prior, 300.0, 0.25, plan_every=4, rec_every=40, callback=lambda s, T_, Z_, ad: rec.append((s, Z_.copy(), ad, apl.energy(T_, Z_))), gate=True)
    Zs = [r[1] for r in rec]; ad = [r[2] for r in rec]; Js = np.array([r[3] for r in rec])
    first = ad.index(True)
    assert first > 0 and all(np.array_equal(Zs[k], prior) for k in range(first))   # frozen while in transit
    assert np.abs(Zs[-1] - prior).max() > 1e-4                                       # adapting afterwards
    assert np.all(np.diff(Js) <= 1e-12)                                              # 𝒥 non-increasing in both phases


def test_the_linearisation_is_self_adjoint_in_the_gradient_flow_metric():
    # (i)–(ii): the joint drift is −M⁻¹∇𝒥 with M = (1/m)G⁻¹ per particle and (σ₀q_j/η)G⁻¹ per prior atom, so M times its Jacobian at the exact
    # representation is minus the Hessian of 𝒥: symmetric and negative semidefinite (the symmetrised pair matrices of (ii))
    loop, exact, prior = _setup()
    for eta in (0.3 * SIGMA0, 3 * SIGMA0):
        apl = AP.AdaptivePriorLoop(loop, eta)
        m = 6
        Th0 = np.repeat(exact, m // 3, axis=0)
        gamma, _ = AP.plan_to(Th0, exact, Q, G)

        def F(v):
            T_ = v[:3 * m].reshape(m, 3); Z_ = v[3 * m:].reshape(3, 3)
            a, b = apl.drift(T_, Z_, gamma)
            return np.concatenate([a.ravel(), b.ravel()])
        v0 = np.concatenate([Th0.ravel(), exact.ravel()]); d = v0.size; h = 1e-6
        Jac = np.zeros((d, d))
        for j in range(d):
            e = np.zeros(d); e[j] = h
            Jac[:, j] = (F(v0 + e) - F(v0 - e)) / (2 * h)
        M = np.concatenate([np.tile(1.0 / (m * G), m), np.concatenate([(SIGMA0 * Q[j] / eta) / G for j in range(3)])])
        S = M[:, None] * Jac
        assert np.abs(S - S.T).max() < 1e-7 * np.abs(S).max()
        assert np.linalg.eigvalsh(-0.5 * (S + S.T)).min() > -1e-9 * np.abs(S).max()


def test_the_settled_gate_opens_at_the_fixed_prior_stationary_point():
    loop, exact, prior = _setup()
    m, tol = 15, 1e-3
    Th = np.repeat(exact, m // 3, axis=0)
    Th[0] = 0.4 * exact[0] + 0.6 * exact[1]                                      # mass in transit at the start
    t_transit = AP.AdaptivePriorLoop(loop, SIGMA0).run(Th, prior, 400.0, 0.25, plan_every=4, rec_every=10 ** 9, gate=True)[2]
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    rec = []

    def cb(s, T_, Z_, ad):
        g_, _ = AP.plan_to(T_, Z_, Q, G)
        rec.append((s, Z_.copy(), ad, apl.fixed_prior_drift(T_, Z_, g_), apl.in_transit(T_, Z_, g_), apl.energy(T_, Z_)))
    t_open = apl.run(Th, prior, 1200.0, 0.25, plan_every=4, rec_every=4, callback=cb, gate='settled', settle_tol=tol)[2]
    first = [r[2] for r in rec].index(True)
    assert np.isfinite(t_transit) and t_open > t_transit + 50.0                  # it waits beyond the end of the transit …
    assert all(np.array_equal(r[1], prior) for r in rec[:first + 1])             # … with the prior frozen until it opens
    assert rec[first][3] < tol and not rec[first][4] and rec[first - 1][3] >= tol  # it opens when the fixed-prior drift first falls below tol
    assert np.all(np.diff([r[5] for r in rec]) <= 1e-12)                         # 𝒥 non-increasing in both phases
