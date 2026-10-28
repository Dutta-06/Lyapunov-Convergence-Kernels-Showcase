"""Theorem D‴.53: the sliding limit — the inner problem is a gradient flow of the fit (an isolated entrant converges to its feature), and
the loop reaches the sliding path with the predicted first-order correction (the entries' pulls and shifts, the creep)."""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import exp_transit_limit as TL
from src.core import rich_transit as RT


def test_isolated_entrant_inner_problem_descends_the_fit_to_its_feature():
    """One bump A·φ(c*) and one entrant from ξ = 200: dξ/ds = −⟨rφ⟩, dc/ds = Λ₀Gm w⟨r∂_cφ⟩, w = 𝒲(ξ) — the fit E = ½⟨r²⟩ never increases
    (Lemma D‴.53a: a gradient flow in the mirror coordinate's metric) and the entrant ends on the feature, w = A, c = c*."""
    L, K, c0 = TL.setup()
    A, cs, ce = 0.8, 0.051, 0.1
    f = A * K.phi(np.array([cs]))[:, 0]
    lam0 = RT.transit_scale(L, 1.0)
    nx = K.xs.size

    def wc(z):
        xi, c = z
        w = float(RT.inner_coef(L, xi))
        return w, f - w * K.phi(np.array([c]))[:, 0]

    def rhs(s, z):
        w, r = wc(z)
        c = z[1]
        return [-(K.phi(np.array([c]))[:, 0] @ r) / nx, lam0 * L.G * L.m * w * (K.dphi(np.array([c]))[:, 0] @ r) / nx]
    ss = np.linspace(0, 6000, 601)
    sol = solve_ivp(rhs, (0, ss[-1]), [200.0, ce], method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ss)
    E = np.array([0.5 * np.mean(wc(z)[1] ** 2) for z in sol.y.T])
    assert np.all(np.diff(E) <= 1e-12 * E[0])
    assert E[-1] < 1e-14 * E[0]
    w_end, _ = wc(sol.y[:, -1])
    assert abs(w_end - A) < 1e-6 and abs(sol.y[1, -1] - cs) < 1e-6


def run(n, eps):
    """D‴.51's target n: its constants (`b_job`: converged pulls, s_θ, (H3)) and the loop against the path at the mass ε (`a_job`)."""
    cs, As = TL.TM.random_targets(TL.SEED_A, TL.N_POOL)[n]
    C = TL.b_job((n, cs, As))
    saved = TL.EPS_A
    try:
        TL.EPS_A = (eps,)
        r = TL.a_job(C)
    finally:
        TL.EPS_A = saved
    return cs, As, C, r['rows'][0]


def test_loop_reaches_the_path_with_its_first_order_correction():
    """A three-bump target whose path leaves margins ≥ 0.15 (D‴.51's target 39), the loop at ε = 10⁻¹⁰: exactly the path's kernels; the
    never-entering kernels on the core's profile; their centres displaced by the creep's εC_j; their duals off the path by ε^{1/2}·(the
    entries' pulls and shifts) + εD_j to 2·10⁻³ε^{1/2}; the entries' crossings at τ_k + ε^{1/2}(s_θ + σ_k) to 0.02ε^{1/2}."""
    cs, As, C, row = run(39, 1e-10)
    assert C['ok'] and C['mu0'] > 0.15
    eps = row['eps']
    assert row['n_form'] == 3
    assert row['block'] < 1e-10
    disp, creep = np.array(row['disp']), np.array(row['creep'])
    big = np.abs(creep) > 1e-3 * np.abs(creep).max()
    assert np.abs(disp[big] / creep[big] - 1).max() < 5e-3
    res = np.array([d[2] - d[3] - d[4] for d in row['dev']])
    assert np.abs(res).max() < 2e-3
    alone = np.array([d[2] - d[5] for d in row['dev']])
    assert np.abs(alone).max() > 0.03                                    # the shifts' term is needed
    assert max(abs(a - b) for _, a, b in row['cross']) < 0.02
    # the never-entering kernels' coefficients: the core's profile at their path margins, a(εb)^{3/2}/(m(2c)^{1/2})·cos(πμ/2)/sin³(πμ/2)
    L, K, c0, f, R, P = TL.path_of(cs, As)
    never = [j for j in range(c0.size) if j not in P['ever']]
    t_end = 1.3 * C['entries'][-1][0]
    k = len(P['taus']) - 1
    U = P['U'][k] + (t_end - P['taus'][k]) * (K.phi(P['centres'][k]).T @ P['res'][k] / K.xs.size)
    m = 1 - np.abs(U[never])
    prof = L.A * eps ** 1.5 / (L.m * np.sqrt(2 * L.c)) * np.cos(np.pi * m / 2) / np.sin(np.pi * m / 2) ** 3
    assert abs(row['w_never'] / prof.max() - 1) < 1e-2


def test_sign_per_entry_and_the_shifts_on_a_sliding_target():
    """D‴.51's target 41 at ε = 10⁻⁹ (referee round 41: its second entrant slides, its third arrives 4.8ε^{1/2} late): the duals to
    1.7·10⁻³ε^{1/2} and the crossings to 0.045ε^{1/2} (about 3× the observed 5.7·10⁻⁴ and 0.015) with converged predictions — the
    extrapolation in xi0 spread below 10⁻⁴ —; (H3): every entry's inner problem ends at the path's S_k; the block within 2·10⁻¹¹."""
    cs, As, C, row = run(41, 1e-9)
    assert max(C['s_spread']) < 1e-4 and max(C['p_spread']) < 1e-5
    assert max(C['h3w']) < 1e-9 and max(C['h3c']) < 1e-9 and all(C['stationary'])
    assert max(abs(s) for s in C['shifts']) > 1.0                       # an entry's shift matters here
    assert row['n_form'] == 3 and row['block'] < 2e-11
    res = np.array([d[2] - d[3] - d[4] for d in row['dev']])
    assert np.abs(res).max() < 1.7e-3
    assert max(abs(a - b) for _, a, b in row['cross']) < 0.045


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print(name, 'ok')
