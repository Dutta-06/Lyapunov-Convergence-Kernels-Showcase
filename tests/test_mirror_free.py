"""Proposition D‴.45: the free mass–amplitude law as a mirror flow (src/core/mirror_free.py)."""
import os
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import mirror_free as MF


def _instance(seed=3, N=5, d=12):
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((N, d)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(d); idx = rng.choice(d, 2, replace=False); ws[idx] = rng.choice([-1, 1], 2) * rng.uniform(0.15, 0.3, 2)
    return Phi, Phi @ ws


def test_eta_closed_form_three_cases():
    """η against the quadrature of dα/(Gβa′) for q > c (also q ≫ c), |q| < c, and on a branch (q < −c, from a reference point)."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    f = lambda x, q: 1.0 / (L.G * L.beta(x, q) * L.amp(x)[1])
    for q in (3.0, 0.5, 1e4 * L.c, 0.03, 0.0, -0.02):
        for al in (-2.5, -0.4, 0.3, 1.7, 4.0):
            ref = quad(lambda x: f(x, q), 0.0, al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
            assert abs(L.eta(al, q)[0] - ref) <= 1e-9 * max(1.0, abs(ref))
    q = -0.05                                                     # a branch: β = 0 at α₀
    a0 = float(L.branch_point(q))
    assert a0 > 0 and abs(L.beta(a0, q)) < 1e-12
    for al in (a0 + 0.2, a0 + 1.0, a0 + 3.0):
        ref = quad(lambda x: f(x, q), a0 + 0.1, al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
        assert abs((L.eta(al, q) - L.eta(a0 + 0.1, q))[0] - ref) <= 1e-9 * max(1.0, abs(ref))


def test_potential_conjugate_and_metric():
    """ψ′ = η, ψ*(η) = (β − β⁰)/(γm) with (ψ*)′ = w, and dw/dη = mD (finite differences)."""
    L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    for q in (0.97, 0.01, -0.02):
        for al in (-1.3, 0.2, 2.2):
            hh = 1e-5
            dw = (L.w(al + hh, q) - L.w(al - hh, q)) / (2 * hh)
            de = (L.eta(al + hh, q) - L.eta(al - hh, q))[0] / (2 * hh)
            dpsi = (L.psi(al + hh, q) - L.psi(al - hh, q))[0] / (2 * hh)
            dps = (L.psi_star(al + hh, q) - L.psi_star(al - hh, q)) / (2 * hh)
            assert abs(dpsi / dw - L.eta(al, q)[0]) < 1e-6 * max(1, abs(L.eta(al, q)[0]))
            assert abs(dps / de - L.w(al, q)) < 1e-6
            assert abs(dw / de - L.m * L.relevance(al, q)) < 1e-6 * L.m * L.relevance(al, q)


def test_conservation_for_any_drive_and_conserved_law_shift():
    """Q constant along the loop with a smooth node-correlated noise; for the conserved law dQ_i/dt = −γβ_iλ."""
    Phi, y = _instance()
    L = MF.FreeLaw()
    rng = np.random.default_rng(1)
    om = np.exp(rng.uniform(-6, 1, 12)); ph = rng.uniform(0, 6.28, 12); Mx = 0.3 * rng.standard_normal((5, 12))
    xi = lambda t: Mx @ np.sin(om * t + ph)
    d = Phi.shape[1]
    z0 = np.r_[0.05 * rng.standard_normal(d), np.log(np.full(d, 0.3))]
    sol = solve_ivp(L.rhs(Phi, y, xi), (0, 3000), z0, method='DOP853', rtol=1e-11, atol=1e-13)
    q0 = L.Q(z0[:d], np.exp(z0[d:]))
    assert np.abs(L.Q(sol.y[:d, -1], np.exp(sol.y[d:, -1])) - q0).max() < 1e-8
    F = L.rhs(Phi, y, xi, conserved=True)
    z = sol.y[:, -1]
    dz = F(1.0, z)
    al, b = z[:d], np.exp(z[d:])
    a, da = L.amp(al)
    dQ = b * dz[d:] - L.c * (2 / L.A) * np.sinh(2 * al / L.A) * dz[:d]
    r = y - Phi @ (b * a / L.m) + xi(1.0); g = Phi.T @ r / Phi.shape[0]
    lam = b @ (a * g) / b.sum()
    assert np.allclose(dQ, -L.gamma * b * lam, rtol=1e-10, atol=1e-14)


def test_noiseless_limit_is_the_bregman_projection():
    Phi, y = _instance()
    L = MF.FreeLaw()
    d = Phi.shape[1]
    al0 = 0.05 * np.random.default_rng(5).standard_normal(d); b0 = np.full(d, 1.0)
    q = L.Q(al0, b0)
    w, al, it, res = L.bregman_projection(Phi, y, q, L.eta(al0, q))
    assert res < 1e-12
    sol = solve_ivp(L.rhs(Phi, y), (0, 4e4), np.r_[al0, np.log(b0)], method='LSODA', rtol=1e-11, atol=1e-13)
    wT = np.exp(sol.y[d:, -1]) * L.amp(sol.y[:d, -1])[0] / L.m
    assert np.abs(wT - w).max() < 1e-9


def test_euler_change_of_Q():
    """The Euler step's change of Q: second order in h with the closed-form coefficient."""
    L = MF.FreeLaw()
    rng = np.random.default_rng(2)
    al = rng.uniform(-2, 2, 8); b = rng.uniform(0.1, 3, 8); g = rng.standard_normal(8)
    for h in (1e-2, 1e-3):
        ex = L.dQ_exact(al, b, g, h); pr = L.dQ_euler(al, b, g ** 2, h)
        assert np.all(np.abs(ex / pr - 1) < 50 * h)


def test_real_loop_lyapunov():
    """V̇ = −k_ee² from the vector field at random states (the free law, a realisable target)."""
    L = MF.FreeLaw(m=6.0, G=1.0, gamma=0.5, a_max=1.0)
    rng = np.random.default_rng(4)
    d = 6
    cen, wid = np.linspace(-1.2, 1.2, d), np.full(d, 0.4)
    al0, b0 = 0.1 * rng.standard_normal(d), rng.uniform(0.5, 1.5, d)
    q = L.Q(al0, b0)
    wstar = rng.uniform(-0.15, 0.15, d)
    al_star = L.alpha_of_w(wstar, q)
    xr, xrd = (lambda t: np.sin(t)), (lambda t: np.cos(t))
    ke = 5.0
    F = L.real_rhs(cen, wid, wstar, xr, xrd, ke)
    for _ in range(5):
        al = al0 + 0.3 * rng.standard_normal(d)
        z = np.r_[rng.uniform(-1, 1), al, np.log(L.beta(al, q))]
        t = rng.uniform(0, 6)
        dz = F(t, z)
        hh = 1e-6
        V = lambda s: L.real_lyapunov(z + s * dz, t + s, xr, al_star, q)
        dV = (V(hh) - V(-hh)) / (2 * hh)
        e = z[0] - xr(t)
        assert abs(dV + ke * e ** 2) < 1e-6 * max(1.0, ke * e ** 2)


def test_mirror_step_keeps_Q_and_moves_eta_by_hg():
    L = MF.FreeLaw()
    rng = np.random.default_rng(6)
    al = rng.uniform(-2, 2, (4, 7)); q = rng.uniform(0.05, 2.0, 7)
    g = rng.standard_normal((4, 7))
    al1 = L.mirror_step(al, q, g, 0.3)
    assert np.allclose(L.eta(al1, np.broadcast_to(q, al.shape)) - L.eta(al, np.broadcast_to(q, al.shape)), 0.3 * g, atol=1e-11)
    assert np.allclose(L.Q(al1, L.beta(al1, q)), np.broadcast_to(q, al.shape), atol=1e-13)
