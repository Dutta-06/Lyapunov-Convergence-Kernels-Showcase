"""
Proposition D‴.40 (PARADIGM.md §5d): the gauge under general noise. Checked here: the noise mask of the dithered integrator dithers only
the chosen coordinates (the amplitudes then move by their drift alone); the averaged fitness of the trace formula is, for isotropic noise
2τ diag(G), minus Proposition D‴.38's correction at κ_F = 1 (−(τ/2)∂_βℓ along M), the same for noise on the amplitudes only, and zero to
leading order for noise on the centres and widths only.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import gauge_noise as GN                                   # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def test_the_noise_mask_dithers_only_the_chosen_coordinates():
    # Proposition D‴.40's mask: with (0, 1, 1) the amplitudes move by their drift alone, the centres and widths as without the mask
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    T0 = E.exact_rep(net, np.ones(3) / 3) + 0.01
    out = {}
    for mask in (None, (0.0, 1.0, 1.0)):
        rec = []
        loop.run_dithered(T0, np.ones(3), 0.25, 0.25, 1e-3, np.random.default_rng(3), rec_every=1, callback=lambda s_, T, b: rec.append(T.copy()),
                          mask=mask)
        out[mask] = rec[-1]
    drift = loop.drift(T0, np.ones(3))[0]
    assert np.allclose(out[(0.0, 1.0, 1.0)][:, 0], T0[:, 0] + 0.25 * drift[:, 0], rtol=0, atol=1e-15)
    assert np.array_equal(out[(0.0, 1.0, 1.0)][:, 1:], out[None][:, 1:]) and not np.allclose(out[None][:, 0], out[(0.0, 1.0, 1.0)][:, 0])


def _instance(W):
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    beta = 3 * np.asarray(W, float)
    G9 = np.tile(E.G, 3)
    gradJ = lambda t: -loop.drift(t.reshape(3, 3), beta)[0].ravel() / G9
    u = lambda t: loop.fitness(loop.residual(t.reshape(3, 3), beta), t.reshape(3, 3))
    return net, base, loop, beta, G9, gradJ, u


def test_the_trace_formula_reduces_to_the_log_determinant_for_isotropic_noise_and_vanishes_without_amplitude_noise():
    tau = 1e-4
    W = np.array([0.30, 0.35, 0.35])
    net, base, loop, beta, G9, gradJ, u = _instance(W)
    t = E.exact_rep(net, W).ravel()
    iso = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(G9), G9)[0]
    fix = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=1.0, fixman_tau=tau).fixman_correction(t.reshape(3, 3), beta)
    flow = lambda ub: GN.replicator_flow(W, ub)
    assert np.abs(flow(iso) + flow(fix)).max() < 1e-5 * np.abs(flow(fix)).max()            # ū = −(τ/2)∂_βℓ up to a constant (2·10⁻⁷)
    amp = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(np.tile([1.0, 0.0, 0.0], 3)), G9)[0]
    assert np.abs(flow(amp) - flow(iso)).max() < 1e-5 * np.abs(flow(iso)).max()
    cs = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(np.tile([0.0, 0.25, 0.25], 3)), G9)[0]
    assert np.abs(flow(cs)).max() < 1e-5 * np.abs(flow(iso)).max()


def test_only_each_amplitudes_own_noise_variance_drives_its_class():
    # the generator identity: ū_k = −τ_k E[g′_k(α_k)], τ_k = Σ_{α_kα_k}/(2G_α), g_k = a_k/(β_k a′_k), so to leading order
    # ū_k = −(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k) — for a correlated noise on all nine coordinates, the correlations and the centres' and
    # widths' noise dropping out
    W = np.array([0.30, 0.35, 0.35])
    net, base, loop, beta, G9, gradJ, u = _instance(W)
    t = E.exact_rep(net, W).ravel()
    a = base.features(E.exact_rep(net, W))[0]
    q = (a / net.a_max) ** 2
    rng = np.random.default_rng(11)
    for _ in range(2):
        L = rng.standard_normal((9, 9))
        S = 2e-4 * 0.3 * (L @ L.T)
        tau_k = np.diag(S)[0::3] / (2 * G9[0])
        ub = GN.averaged_fitness(gradJ, u, t, S, G9)[0]
        assert np.abs(ub / (-(tau_k / 3) * (1 + q) / (W * (1 - q))) - 1).max() < 1e-5


def test_the_weighted_analytic_centre():
    import exp_gauge_fixman as Q
    import exp_gauge_noise as N
    assert np.abs(N.w_sharp_weighted([1.0, 1.0, 1.0]) - Q.w_sharp()).max() < 1e-12
    tk = N.levels(N.noises()['correlated 1'])
    w = N.w_sharp_weighted(tk)
    Wm = np.abs(E.A_BUMP) / 4.0
    grad = tk * (1 + Wm ** 2 / w ** 2) / (w - Wm ** 2 / w)                   # Lagrange: the weighted partials equal
    assert abs(w.sum() - 1) < 1e-12 and np.ptp(grad) < 1e-9 * grad.mean()


def test_the_correlated_noise_path_of_the_dithered_integrator():
    """`run_dithered` with `noise_L`: one step is Θ + Δt·drift + √Δt·L ξ with ξ the generator's next standard normal draws (referee round 31:
    the path was untested)."""
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    T0 = E.exact_rep(net, np.ones(3) / 3) + 0.01
    L = np.linalg.cholesky(np.eye(9) + 0.3 * np.ones((9, 9)))
    rec = []
    loop.run_dithered(T0, np.ones(3), 0.25, 0.25, 1e-3, np.random.default_rng(5), rec_every=1, callback=lambda s_, T, b: rec.append(T.copy()),
                      noise_L=L)
    xi = np.random.default_rng(5).standard_normal(9)
    expect = T0 + 0.25 * loop.drift(T0, np.ones(3))[0] + (np.sqrt(0.25) * L @ xi).reshape(3, 3)
    assert np.allclose(rec[-1], expect, rtol=0, atol=1e-14)


def test_the_levels_and_the_per_class_correction():
    import exp_gauge_noise as N
    tau = 1e-4
    assert np.allclose(N.levels(2 * tau * np.diag(np.tile(E.G, 3))), tau)                     # τ_k = Σ_{α_kα_k}/(2G_α)
    assert np.allclose(N.levels(N.noises(tau)['correlated E']), np.array(N.LEVELS_E) * tau)
    cE = N.corr_E()
    assert np.allclose(np.diag(cE)[1::3], E.G[1]) and np.allclose(np.diag(cE)[2::3], E.G[2]) and np.linalg.eigvalsh(cE).min() > 0
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35]); Th = E.exact_rep(net, W); beta = 3 * W
    A0 = base.features(Th)[0] * W
    one = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau, fixman_A=A0).fixman_correction(Th, beta)
    per = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=np.array(N.LEVELS_E) * tau,
                                    fixman_A=A0).fixman_correction(Th, beta)
    assert np.allclose(per, np.array(N.LEVELS_E) * one, rtol=1e-14)
    a = base.features(Th)[0]; q = (a / net.a_max) ** 2
    assert np.allclose(N.closed(W, [1.0, 2.0, 3.0]), -(np.array([1.0, 2.0, 3.0]) / 3) * (1 + q) / (W * (1 - q)))
    ub = np.array([0.3, -0.1, 0.2])
    assert np.allclose(GN.replicator_flow(W, ub), W * (ub - W @ ub)) and abs(GN.replicator_flow(W, ub).sum()) < 1e-15


def test_the_gauge_identity_pointwise_and_the_generator_applied_to_psi():
    """u_k = −g_k∂_{α_k}J_fit at every state (g_k = a_k/(β_ka′_k)); Part D's Lψ_k = u_k − g_k∂_{α_k}(J_coh + V) + τ_kg′_k equals
    ψ′b + (Σ_αα/2)ψ″ for ψ_k = (a_max²/4β_kG_α)cosh(2α_k/a_max); the chain's conditional expectation of ψ(α + bΔt + ν) in closed form."""
    import exp_gauge_plateau as P
    net, base = E.setup()
    beta = np.array([0.9, 1.05, 1.05])
    fit = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    rng = np.random.default_rng(2)
    Th = E.exact_rep(net, beta / 3) + 0.05 * rng.standard_normal((3, 3)); Th[0, 0] = 5.6        # α_1 outside the box
    a = base.features(Th)[0]; am = net.a_max; ap = 1 - (a / am) ** 2
    g = a / (beta * ap); gp = (1 + (a / am) ** 2) / (beta * ap)
    r = fit.residual(Th, beta); u = fit.fitness(r, Th)
    Jfit = lambda T: 0.5 * 3 * np.mean(fit.residual(T, beta) ** 2)                                # J = m·½⟨r²⟩, m = 3
    h = 1e-6
    dJ = np.array([(Jfit(Th + h * np.eye(3)[k][:, None] * np.eye(3)[0][None, :]) - Jfit(Th - h * np.eye(3)[k][:, None] * np.eye(3)[0][None, :]))
                   / (2 * h) for k in range(3)])
    assert np.allclose(u, -g * dJ, rtol=1e-6, atol=1e-12)
    Sig = 2e-3 * np.diag(np.tile(E.G, 3)); tk = np.diag(Sig)[0::3] / (2 * E.G[0])
    dT = loop.drift(Th, beta)[0]
    nonfit = -(dT[:, 0] - beta * E.G[0] * base.grad_weighted(r, Th)[:, 0]) / E.G[0]
    lpsi = u - g * nonfit + tk * gp
    c = am ** 2 / (4 * beta * E.G[0]); kk = 2 / am
    psi1 = c * kk * np.sinh(kk * Th[:, 0]); psi2 = c * kk * kk * np.cosh(kk * Th[:, 0])
    assert np.allclose(lpsi, psi1 * dT[:, 0] + 0.5 * np.diag(Sig)[0::3] * psi2, rtol=1e-10)
    assert np.allclose(psi1, g / E.G[0], rtol=1e-12) and np.allclose(psi2, gp / E.G[0], rtol=1e-12)
    x, wq = np.polynomial.hermite_e.hermegauss(60)                                                # E over ν ~ N(0, s²)
    s2, dt = np.diag(Sig)[0::3] * 0.25, 0.25
    num = np.array([np.sum(wq * c[k] * np.cosh(kk * (Th[k, 0] + dt * dT[k, 0] + np.sqrt(s2[k]) * x))) / np.sqrt(2 * np.pi) for k in range(3)])
    assert np.allclose(num, c * np.cosh(kk * (Th[:, 0] + dt * dT[:, 0])) * np.exp(0.5 * kk * kk * s2), rtol=1e-12)
