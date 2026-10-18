"""Theorem D‴ (PARADIGM.md §5d): the averaged slow loop, the three terms of the coupling identity, the residual-weighted
curvature against Proposition K, the whitened constants, and a short warm-start coupled run against the theorem's bound."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import periodic as PO
from src.core import chaos_learning as CL

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])


F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
X_REF_DOT = lambda t: np.array([2.0 * np.cos(t)])


def _instance(sigma0=0.014, n_quad=400):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    loop = CL.AveragedSlowLoop(net, F_VEC, X_REF, PER, G, sigma0, n_quad=n_quad)
    cl = PO.ClosedLoop(net, lambda Xa: F_VEC(Xa[:, 0])[:, None], X_REF, X_REF_DOT, PER, 25.0, G, sigma0)
    return net, loop, cl


def test_averaged_slow_loop_is_the_G_gradient_flow_of_J_m():
    """(𝔄_m): the drift equals −G∇J_m with J_m of periodic.ClosedLoop (same quadrature grid), to machine precision."""
    net, loop, cl = _instance()
    orb = cl.orbit(N=400)
    rng = np.random.default_rng(3)
    Th = CL.cold_start(net, 12, rng)
    assert np.abs(loop.grad_J(Th) - cl.grad_J_m(Th, orb)).max() < 1e-13
    assert np.abs(loop.drift(Th) + G[None, :] * cl.grad_J_m(Th, orb)).max() < 1e-13
    J = 12 * loop.L_bar(Th) + 0.5 * loop.sigma0 * float(np.sum(Th ** 2 / G[None, :]))
    assert abs(J - cl.J_m(Th, orb)) < 1e-12


def test_coupling_identity_three_terms_and_the_sign_of_the_interaction():
    """½dD/ds = T_A + T_B + T_C exactly; the interaction term T_B equals −‖f̂_Θ − f̂_Θ̄‖² up to O(D³/²) (flat convexity: it is
    non-positive at linear order) and obeys T_B ≤ (L₂^G)²D²/16 at every scale of δ."""
    net, loop, _ = _instance()
    rng = np.random.default_rng(3)
    Th = CL.cold_start(net, 12, rng)
    Thb = Th + 0.05 * rng.standard_normal(Th.shape)
    fh_mf = loop.f_hat(CL.cold_start(net, 500, rng))                     # a stand-in for f̂[μ_s]
    d = loop.decomposition(Th, Thb, fh_mf)
    assert abs(d['total'] - d['half_dD']) < 1e-15 * max(1.0, abs(d['half_dD']))
    Lc = CL.whitened_lipschitz_constants(net, G, 2.5, n_samples=20000)
    ratios = []
    for eps in (1.0, 0.3, 0.1, 0.03, 0.01):
        dd = loop.decomposition(Thb + eps * (Th - Thb), Thb)
        assert dd['T_C'] == 0.0
        assert dd['T_B'] <= Lc['L_2_G'] ** 2 * dd['D'] ** 2 / 16
        assert dd['T_B'] < 0
        ratios.append(dd['T_B'] / (-dd['h_norm'] ** 2))
    assert abs(ratios[-1] - 1) < 1e-3 and abs(ratios[0] - 1) < 0.05 and all(abs(r - 1) <= abs(ratios[0] - 1) + 1e-12 for r in ratios)


def test_residual_curvature_at_a_collapsed_configuration_is_proposition_K():
    """κ_R at the atoms of a collapsed configuration with the collapsed residual equals λ_max(G^{1/2}R₂G^{1/2}) of Proposition K,
    which at σ₀ = σ_crit equals σ₀ (the window's upper end); the collapsed critical point has zero drift."""
    net, loop, cl = _instance()
    orb = cl.orbit(N=400)
    g = net.unconstrained(np.array([[1.0]]), np.array([[0.0]]), np.array([0.8]))[0]
    s_crit, tb, R2 = PO.sigma_crit(cl, orb, g)
    loop_c = CL.AveragedSlowLoop(net, F_VEC, X_REF, PER, G, s_crit, n_quad=400)
    Thc = np.tile(tb[None, :], (5, 1))
    kr = loop_c.kappa_R(loop_c.residual(Thc), Thc)
    lam = np.linalg.eigvalsh(np.sqrt(G)[:, None] * R2 * np.sqrt(G)[None, :])[-1]
    assert np.abs(kr - lam).max() < 1e-8 and abs(lam - s_crit) < 1e-9
    assert np.abs(loop_c.drift(Thc)).max() < 1e-12
    assert loop_c.kappa_R_neighbourhood(loop_c.residual(Thc), Thc, 0.05, n_dirs=3, rng=np.random.default_rng(0)) >= lam - 1e-12


def test_whitened_constants_below_their_closed_forms_and_the_theorem_constants():
    """The searched L_ϑ^G, L₂^G lie below the closed-form bounds; theorem_constants reproduces C₁ = 4e^{2Λ₁}(L_ϑ^G)²K_Φ²/κ²
    and the thresholds are positive; Λ₁ of a profile that satisfies (H) from the start is zero."""
    net, loop, _ = _instance()
    Lc = CL.whitened_lipschitz_constants(net, G, 2.5, n_samples=20000)
    assert 0 < Lc['L_th_G'] <= Lc['L_th_G_closed'] * (1 + 1e-9) and 0 < Lc['L_2_G'] <= Lc['L_2_G_closed'] * (1 + 1e-9)
    c = CL.theorem_constants(loop, Lc['L_th_G'], Lc['L_2_G'], R0=3.0, kappa=0.005, Lambda1=0.0, rho_prime=0.3)
    assert abs(c['C1'] - 4 * Lc['L_th_G'] ** 2 * 16.0 / 0.005 ** 2) < 1e-6 * c['C1']
    assert c['rho'] > 0 and c['theta'] > 0 and c['Delta'] > 0 and c['m0'] > 0 and c['R_Theta'] >= 3.0
    s = np.linspace(0, 10, 101)
    assert CL.Lambda1_from_profile(s, 0.01 + 0 * s, 0.005, 5.0) == 0.0
    assert abs(CL.Lambda1_from_profile(s, -0.005 + 0 * s, 0.005, 5.0) - 0.05) < 1e-9


def test_warm_start_coupled_run_contracts_and_stays_below_the_bound():
    """A short warm-start run (m = 16, 64 against M = 512): along the run κ_R at the reference stays below σ₀ (hypothesis (H)
    with Λ₁ = 0), the coupled mean-square distance is far below C₁/m, and the function-level distance ‖f̂_Θ − f̂_Θ̄‖ is below
    L_ϑ^G·mean‖δ_i‖ (the Lipschitz bound used in the proof)."""
    net, loop, _ = _instance(sigma0=0.014)
    rng = np.random.default_rng(11)
    Th_all = CL.warm_start(net, BUMPS, 512, rng)
    rec = loop.coupled_run(Th_all, [16, 64], S_end=40.0, dt=0.02, rec_every=50, rng=rng)
    Lc = dict(L_th_G=4.0436, L_2_G=11.111)                                # searched values (test above)
    assert rec['kR_ref'].max() < loop.sigma0                              # (H) along the reference, margin κ := σ₀ − max κ_R
    kappa = loop.sigma0 - rec['kR_ref'].max()
    c = CL.theorem_constants(loop, Lc['L_th_G'], Lc['L_2_G'], R0=4.0, kappa=kappa, Lambda1=0.0, rho_prime=0.3)
    for m in (16, 64):
        D = rec['D'][m]
        assert D.max() < c['C1'] / m
        assert np.all(rec['hnorm'][m] <= Lc['L_th_G'] * rec['dev1'][m] * (1 + 1e-9))
        assert rec['hnorm_mf'][m].max() ** 2 <= 2 * Lc['L_th_G'] ** 2 * rec['D'][m].max() + 2 * 16.0 * (1.0 / m + 1.0 / 512)   # (c) with the M-reference's own fluctuation
    assert rec['D'][64].max() < rec['D'][16].max()                        # the m-law in the right direction
