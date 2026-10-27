"""Proposition D‴.51: the transit at a finite mass — the inner profile, the entrant's overshoot, the near-tie inner problem, the path's
duals and the predicted mass."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
from src.core import mirror_free as MF
from src.core import rich_transit as RT

L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)


def grid():
    return RT.Bumps(np.linspace(-1, 1, 61), 0.12), np.linspace(-1, 1, 21)


def test_inner_profile_is_the_small_mass_limit_of_the_mirror_map():
    """A kernel of mass ε at α has its dual (λ₁ − η(α))/(ε^{1/2}Λ) below the wall; as ε → 0 this is inner_xi(w(α)) — the critical
    curve's η_c, not the core's power law (which is off by more than 20% below ξ ≈ 15)."""
    al = np.array([0.05, 0.2, 0.5, 1.0, 3.3])
    for eps, tol in ((1e-7, 5e-3), (1e-9, 5e-4)):
        q = L.Q(np.zeros(1), np.array([eps]))
        Lam = RT.transit_scale(L, eps)
        lam1 = Lam * (1 + 0.0)
        xi_lib = np.array([(lam1 - float(L.eta(np.array([a]), q)[0])) / (np.sqrt(eps) * Lam) for a in al])
        w = np.array([float(L.w(np.array([a]), q)[0]) for a in al])
        assert np.allclose(RT.inner_xi(L, w), xi_lib, rtol=0, atol=tol)
    xi = np.array([-9.0, -1.0, 0.0, 2.0, 40.0, 500.0])
    assert np.allclose(RT.inner_xi(L, RT.inner_coef(L, xi)), xi, rtol=0, atol=1e-9)
    Kp3 = (8 / np.pi ** 3) * L.A / (L.m * np.sqrt(2 * L.c))
    assert abs(RT.inner_coef(L, 500.0) * 500.0 ** 3 / Kp3 - 1) < 2e-3                      # the core's power law far below the wall
    assert abs(RT.inner_coef(L, 15.0) * 15.0 ** 3 / Kp3 - 1) > 0.05


def test_overshoot_on_the_exponential_branch():
    """The overshoot −inner_xi(A) = η_c(α_A)/(Λ₀ρ); for large A the branch is exponential, w ≈ (ca/(2m))e^{2α/a}, so the overshoot grows like
    (a²/(4GcΛ₀))·log A, and two kernels whose duals differ by Δξ carry coefficients in the ratio e^{−κΔξ}, κ = πρ(2c)^{1/2}/a."""
    A = np.array([0.3, 0.8, 2.0, 20.0, 200.0])
    ov = RT.overshoot(L, A)
    assert np.allclose(ov, -RT.inner_xi(L, A))
    lam0 = RT.transit_scale(L, 1.0)
    slope = (ov[-1] - ov[-2]) / np.log(10.0)
    assert abs(slope / (L.A ** 2 / (4 * L.G * L.c * lam0)) - 1) < 1e-3
    kap = np.pi * np.sqrt(2 * L.c) / L.A
    w1, w2 = RT.inner_coef(L, -6.0), RT.inner_coef(L, -6.0 + 1.5)
    assert abs(np.log(w1 / w2) / (kap * 1.5) - 1) < 2e-3


def test_near_tie_inner_problem():
    """At an exact tie (x = 0) the inner problem splits the bump equally; at small x the split is the Gibbs one 1/(1 + e^{κx}); far behind,
    the runner-up keeps the frozen coefficient K'³ρ⁻³(x − Δ)⁻³, Δ the pull a frozen runner-up feels (the overshoot less I∞ > 0)."""
    K, c0 = grid()
    kap = np.pi * np.sqrt(2 * L.c) / L.A
    Kp3 = (8 / np.pi ** 3) * L.A / (L.m * np.sqrt(2 * L.c))
    r0 = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 0.0)
    assert abs(r0['share'] - 0.5) < 1e-6 and abs(r0['c_e'] - 0.05) < 1e-6 and abs(r0['c_r'] - 0.05) < 1e-6
    r1 = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 1.0)
    assert abs(r1['share'] * (1 + np.exp(kap)) - 1) < 0.03
    fr = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 200.0, both=False)
    delta = 200.0 - fr['xi_r']
    ov = float(RT.overshoot(L, 0.8))
    assert abs(fr['xi_e'] + ov) < 1e-4 and 0.5 < ov - delta < 2.0                         # the entrant ends at its overshoot; I∞ ∈ (0.5, 2)
    r3 = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 300.0)
    assert abs(r3['w_r'] * (300.0 - delta) ** 3 / Kp3 - 1) < 0.01
    rn = RT.near_tie_inner(L, K, -0.8, 0.05, 0.1, 0.0, 3.0)
    rp = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, 3.0)
    assert abs(rn['share'] / rp['share'] - 1) < 1e-6


def test_near_tie_loop_against_the_inner_problem():
    """The loop at ε = 10⁻⁹ on a bump 10⁻³ off the midpoint (x ≈ 110): the runner-up's share is the inner problem's at the same x to 2%;
    at ε = 10⁻⁶ (x ≈ 3.5) both kernels form, the runner-up with the Gibbs share to 10%."""
    import exp_mirror_transit as X
    K, c0 = grid()
    f = X.bumps_target(K, (0.051,), (0.8,))
    g = K.phi(c0).T @ f / K.xs.size
    mu = 1 - abs(g[10]) / abs(g[11])
    kap = np.pi * np.sqrt(2 * L.c) / L.A
    for eps, tol in ((1e-9, 0.02), (1e-6, 0.10)):
        taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, 1.32 / abs(g[11]), 'paradigm', n_t=111)
        share = abs(W[-1][10]) / (abs(W[-1][10]) + abs(W[-1][11]))
        x = mu / np.sqrt(eps)
        inner = RT.near_tie_inner(L, K, 0.8, 0.05, 0.1, 0.0, x)
        assert abs(share / inner['share'] - 1) < tol
        if eps == 1e-6:
            assert abs(share * (1 + np.exp(kap * x)) - 1) < 0.1


def three_bumps(n, seed=21):
    rng = np.random.default_rng(seed)
    for k in range(n + 1):
        while True:
            cs = np.sort(rng.uniform(-0.8, 0.8, 3))
            if np.diff(cs).min() >= 0.3:
                break
        As = rng.choice([-1, 1], 3) * rng.uniform(0.4, 1.0, 3)
    return cs, As


def test_path_duals_reach_the_walls_at_the_entries():
    """The reconstructed duals: every entrant's dual reaches its wall exactly at its entry (it is pinned there by the path's own entry
    rule), and the kernels that never enter stay inside their boxes all along."""
    import exp_mirror_transit as X
    K, c0 = grid()
    cs, As = three_bumps(10)
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(21))
    P = RT.path_duals(K, f, R)
    taus = np.array(P['taus'])
    nx = K.xs.size
    for t, j, kind, c in R['events']:
        if kind != 'enter':
            continue
        k = int(np.argmin(np.abs(taus - t)))
        u_before = P['U'][k - 1][j] + (taus[k] - taus[k - 1]) * (K.phi(P['centres'][k - 1]).T @ P['res'][k - 1] / nx)[j]
        assert abs(abs(u_before) - 1) < 1e-8
    never = [j for j in range(21) if j not in P['ever']]
    assert np.abs(P['U'][:, never]).max() < 1


def test_recovery_threshold_on_the_lost_near_tie():
    """D‴.48's target 10: the path leaves kernel 15 0.0065 inside its wall beside an unfitted bump — the predicted mass is the creep's,
    ≈ 5·10⁻⁹ (the loop forms the path's kernels at 2·10⁻⁹ and not at 10⁻⁸, `exp_transit_mass` C); the creep is monotone in ε."""
    import exp_mirror_transit as X
    K, c0 = grid()
    cs, As = three_bumps(10)
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(21))
    p = RT.recovery_threshold(L, K, f, c0, R, 1e-2 * np.abs(As).max())
    assert 3e-9 < p['eps_creep'] < 8e-9 and min(p['creep'], key=p['creep'].get) == 15
    assert p['eps_c'] == p['eps_creep'] and p['eps_entry'] > 1e-7
    P = RT.path_duals(K, f, R)
    assert RT.creep_reaches(L, K, P, 15, 2e-8, c0[15]) and not RT.creep_reaches(L, K, P, 15, 1e-9, c0[15])
