"""Proposition D‴.48: the rich transit with moving centres — the sliding inverse scale space, the frozen inactive kernels, the
Wasserstein–Fisher–Rao metric's condensation."""
import numpy as np
from scipy.integrate import solve_ivp

from src.core import mirror_free as MF
from src.core import rich_transit as RT

L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
K = RT.Bumps(np.linspace(-1, 1, 61), 0.12)
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
FSTAR = K.phi(CSTAR) @ ASTAR
C0 = np.linspace(-1, 1, 21)
D = C0.size


def flow(eps, tau_end, metric, n_t=221):
    Lam = RT.transit_scale(L, eps)
    taus = np.linspace(0, tau_end, n_t)
    sol = solve_ivp(RT.moving_rhs(L, K, FSTAR, metric), (0, Lam * tau_end), np.r_[np.zeros(D), np.log(np.full(D, eps)), C0],
                    method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
    al, lb, cc = sol.y[:D].T, sol.y[D:2 * D].T, sol.y[2 * D:].T
    return taus, np.exp(lb) * L.amp(al)[0] / L.m, cc, al, Lam


def test_critical_curve_and_the_depth_three_unit():
    """On Q = −c the coefficient is FreeLaw's, and w ≈ 2cα³/(ma_max²) at small α."""
    al = np.array([-1.3, -0.2, 0.01, 0.4, 2.0])
    assert np.allclose(RT.critical_w(L, al), [L.w(a, -L.c) for a in al], rtol=1e-12, atol=0)
    a = 1e-3
    assert abs(RT.critical_w(L, a) / (2 * L.c * a ** 3 / (L.m * L.A ** 2)) - 1) < 1e-5


def test_sliding_iss_recovers_the_bumps():
    """The sliding inverse scale space on three separated bumps: three entries (the kernels started nearest to the bumps, the largest
    amplitude first), each sliding onto its bump; the fit exact."""
    R = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D))
    assert [e[1] for e in R['events']] == [5, 16, 11] and all(e[2] == 'enter' for e in R['events'])
    order = np.argsort(R['c'][R['act']])
    assert np.abs(R['c'][R['act']][order] - CSTAR).max() < 1e-6
    assert np.abs(R['w'][R['act']][order] - ASTAR).max() < 1e-6
    assert R['hist'][-1][3] < 1e-9


def test_paradigm_flow_follows_the_sliding_path():
    """Proposition D‴.48(iii): at ε = 10⁻⁵ the loop forms the same three kernels, entering at the sliding path's times, at the bumps;
    the others stay where they started (their mass negligible)."""
    R = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D))
    t_ev = {e[1]: e[0] for e in R['events']}
    taus, W, cc, al, Lam = flow(1e-5, 1.3 * R['events'][-1][0], 'paradigm', n_t=441)
    act = np.abs(W[-1]) > 1e-2 * np.abs(ASTAR).max()
    assert sorted(np.flatnonzero(act)) == [5, 11, 16]
    for i in (5, 11, 16):
        k = np.flatnonzero(np.abs(W[:, i]) > 1e-2 * np.abs(ASTAR).max())[0]
        assert abs(taus[k] - t_ev[i]) < 0.15
    order = np.argsort(cc[-1][act])
    assert np.abs(cc[-1][act][order] - CSTAR).max() < 1e-3
    assert np.abs(cc[-1][~act] - C0[~act]).max() < 0.05 and np.abs(W[-1][~act]).max() < 1e-3


def test_inactive_kernels_obey_the_frozen_bound():
    """Proposition D‴.48(iii)(a): before the first entry each kernel's coefficient is below a_max ε^{3/2}μ^{−3}/(m(2c)^{1/2}) (μ its dual
    coordinate's margin), so its centre hardly moves."""
    eps = 1e-5
    taus, W, cc, al, Lam = flow(eps, 9.0, 'paradigm', n_t=91)
    q = L.Q(np.zeros(D), np.full(D, eps))
    lam1 = np.pi * L.A * np.sqrt((L.c - q) / (q + L.c)) / (4 * L.G * L.c)
    for k in range(1, len(taus)):
        mu = 1 - np.abs(L.eta(al[k], q)) / lam1
        assert mu.min() > 0
        bound = L.A * eps ** 1.5 / (L.m * np.sqrt(2 * L.c)) / mu ** 3
        assert np.all(np.abs(W[k]) <= 1.05 * bound)
    assert np.abs(cc[-1] - C0).max() < 1e-2


def test_wfr_metric_condenses_onto_the_limit():
    """Proposition D‴.48(iv): with the Wasserstein–Fisher–Rao metric the inactive kernels' centres follow the limit
    dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ⟩, du/dτ = ⟨f*φ⟩ before the first entry, gathering at three positions."""
    tau = 10.0
    lim = solve_ivp(RT.wfr_inactive_rhs(L, K, FSTAR, np.ones(D)), (0, tau), np.r_[np.zeros(D), C0], method='LSODA', rtol=1e-10,
                    atol=1e-12)
    c_lim = np.sort(lim.y[D:, -1])
    assert np.sum(np.diff(c_lim) > 1e-3) + 1 == 3
    taus, W, cc, al, Lam = flow(1e-6, tau, 'wfr', n_t=11)
    assert np.abs(np.sort(cc[-1]) - c_lim).max() < 1e-5
    assert np.abs(W[-1]).max() < 1e-3
