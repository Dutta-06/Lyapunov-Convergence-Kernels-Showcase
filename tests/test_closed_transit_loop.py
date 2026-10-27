"""Proposition D‴.49: the closed loop in C (`closed_loop`) against a tight Python integration of the same equations, Proposition D‴.46's
leak in it, and the exact averaged drive along the frozen periodic orbit (`rich_transit.FrozenOrbit`)."""
import numpy as np
from scipy.integrate import solve_ivp

from src.core import closed_loop as CL
from src.core import mirror_free as MF
from src.core import rich_transit as RT

L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
SIG, KE, OM = 0.12, 10.0, 0.2
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
C0 = np.linspace(-1, 1, 21)


def fstar(x):
    return np.exp(-(np.asarray(x, float)[..., None] - CSTAR) ** 2 / (2 * SIG ** 2)) @ ASTAR


def python_loop(eps, T, ts, metric='paradigm', lam=0.0):
    """The same loop on (x, α, log β, c) by LSODA at tight tolerances (the leak through mirror_free's η)."""
    d = C0.size
    q = L.Q(np.zeros(d), np.full(d, eps))

    def f(t, z):
        x = z[0]
        al, lb, c = z[1:1 + d], z[1 + d:1 + 2 * d], z[1 + 2 * d:]
        b = np.exp(lb); a, da = L.amp(al)
        ph = np.exp(-(x - c) ** 2 / (2 * SIG ** 2))
        e = x - np.sin(OM * t)
        u = e * ph - (lam * L.eta(al, q) if lam else 0.0)
        v = L.G * b * a if metric == 'paradigm' else (L.G * a if metric == 'wfr' else 0.0 * a)
        return np.r_[fstar(x) - ph @ (b * a / L.m) + OM * np.cos(OM * t) - KE * e, L.G * b * da * u, L.gamma * a * u,
                     v * e * ph * (x - c) / SIG ** 2]
    sol = solve_ivp(f, (0, T), np.r_[0.0, np.zeros(d), np.log(np.full(d, eps)), C0], method='LSODA', rtol=1e-10, atol=1e-13,
                    t_eval=ts, max_step=1.0)
    W = np.exp(sol.y[1 + d:1 + 2 * d].T) * L.amp(sol.y[1:1 + d].T)[0] / L.m
    return sol.y[0], W, sol.y[1 + 2 * d:].T


def test_c_loop_is_the_python_loop():
    """At ε = 10⁻² through the kernels' entries (t ≤ 1.2·10⁴): the state to 10⁻⁷, under the paradigm's metric and with fixed centres."""
    T = 1.2e4
    ts = np.linspace(0, T, 41)
    for metric in ('paradigm', 'fixed'):
        x, W, cc = python_loop(1e-2, T, ts, metric)
        r = CL.simulate(L, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), ts, metric=metric)
        assert np.abs(x - r['x']).max() < 1e-7 and np.abs(W - r['w']).max() < 1e-7 and np.abs(cc - r['c']).max() < 1e-7
        assert np.abs(W[-1]).max() > 0.05                                 # kernels have formed


def test_c_loop_leak_is_the_mirror_coordinate_leak():
    """Proposition D‴.46's leak −λη(α) on the drive (η on the curve from α = 0 with the mass ε, closed form) against mirror_free's η."""
    T = 6000.0
    ts = np.linspace(0, T, 31)
    x, W, cc = python_loop(1e-2, T, ts, 'fixed', lam=0.05)
    r = CL.simulate(L, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), ts, metric='fixed', lam=0.05)
    assert np.abs(W - r['w']).max() < 1e-7 and np.abs(x - r['x']).max() < 1e-7
    r0 = CL.simulate(L, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), ts, metric='fixed')
    assert np.abs(r['w'][-1]).sum() < np.abs(r0['w'][-1]).sum()        # the leak shrinks the coefficients


def test_frozen_orbit_is_the_periodic_solution():
    """FrozenOrbit's e_p solves ė = −k_ee + r(x_m + e) periodically (integrated over a period from e_p(0) it returns to e_p and follows
    the samples), and its truncation e_p ≈ L[r(x_m)] (OrbitBumps) holds to O(|r|²): halving r halves e_p up to a quadratic remainder."""
    T = 2 * np.pi / OM
    M = 512
    c = np.array([-0.5, 0.1, 0.6]); w = np.array([0.5, -0.3, 0.4])
    O = RT.FrozenOrbit(lambda t: np.sin(OM * t), T, M, SIG, KE, fstar)
    e = O.orbit(c, w)
    ts = T * np.arange(M) / M

    def rhs(t, y):
        x = np.sin(OM * t) + y[0]
        return [-KE * y[0] + fstar(np.array([x]))[0] - np.exp(-(x - c) ** 2 / (2 * SIG ** 2)) @ w]
    sol = solve_ivp(rhs, (0, T), [e[0]], rtol=1e-11, atol=1e-13, t_eval=np.r_[ts, T], method='DOP853')
    assert np.abs(sol.y[0][:-1] - e).max() < 1e-5 and abs(sol.y[0][-1] - e[0]) < 1e-5
    K = RT.OrbitBumps(lambda t: np.sin(OM * t), T, M, SIG, KE)
    lin = K.filt(fstar(K.xs) - K.phi(c) @ w)
    gap = np.abs(e - lin).max()
    assert 0 < gap < 0.2 * np.abs(e).max()                               # the truncation's error, not small at these parameters
    g, h, _ = O.drive(c, w)
    xp = np.sin(OM * ts) + e
    ph = np.exp(-(xp[:, None] - c[None, :]) ** 2 / (2 * SIG ** 2))
    assert np.allclose(g, ph.T @ e / M, atol=1e-12) and np.allclose(h, (ph * (xp[:, None] - c[None, :]) / SIG ** 2).T @ e / M, atol=1e-12)


def test_exact_averaged_loop_tracks_the_real_loop():
    """At ε = 10⁻⁴: the averaged loop along the frozen orbit forms the real loop's five kernels, entering within 0.5% of its times (the
    leading-order averaged loop forms six, D‴.49 A)."""
    eps = 1e-4
    d = C0.size
    Lam = RT.transit_scale(L, eps)
    O = RT.FrozenOrbit(lambda t: np.sin(OM * t), 2 * np.pi / OM, 512, SIG, KE, fstar)
    taus = np.linspace(0, 350, 1401)
    sol = solve_ivp(RT.frozen_orbit_rhs(L, O, np.full(d, eps)), (0, Lam * 350), np.r_[np.zeros(d), C0], method='LSODA', rtol=1e-9,
                    atol=1e-12, t_eval=Lam * taus)
    al = sol.y[:d].T
    W = (eps + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m
    th = 1e-2 * 0.9
    ent_avg = {i: taus[np.flatnonzero(np.abs(W[:, i]) > th)[0]] for i in range(d) if (np.abs(W[:, i]) > th).any()}
    r = CL.simulate(L, eps, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), np.linspace(0, Lam * 350, 11), th=th)
    ent_real = {i: t / Lam for i, t in enumerate(r['entries']) if np.isfinite(t)}
    assert set(ent_avg) == set(ent_real)
    assert len(ent_avg) == 5 and all(abs(ent_avg[i] / ent_real[i] - 1) < 0.005 for i in ent_avg)
