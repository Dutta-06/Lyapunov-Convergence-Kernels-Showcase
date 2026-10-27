"""
Proposition D‴.38 (PARADIGM.md §5d): a Fixman correction for the replicator under dither. Checked here: the loop's correction equals
(τ/2)∂_βℓ along the manifold of exact fits (ℓ = log det ∇²_ΘJ, by finite differences of the Hessian's determinant); it vanishes without
dither and keeps Σβ = m; in the slow-replicator limit it cancels the dither's drift; with κ_F = 2 the corrected field's zero is the unique
maximiser W♯ of ℓ (the D-optimal exact fit); near a saturation the one-dimensional zero-flux density behaves like δ^{κ_F − 1}, and below
the floor's cap the exponent falls back as stated; the tracked average relaxes at its time constant; kappa() is the leak's projected stiffness;
the long runs' and the paired analysis' statistics are right on synthetic records; the correction's general form,
through the exact fit's implicit function and computed from the current point only, reduces to the closed form.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
import exp_gauge_dither as GDD                                           # noqa: E402
import exp_gauge_drift as X                                              # noqa: E402
import exp_gauge_fixman as Q                                             # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def test_the_correction_is_the_log_determinants_derivative_along_the_manifold():
    net, base = E.setup()
    tau = 1e-4
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=1.0, fixman_tau=tau)
    beta = 3 * np.array([0.30, 0.35, 0.35]); hb = 1e-3
    grad = np.array([(X.log_det(beta + hb * e) - X.log_det(beta - hb * e)) / (2 * hb) for e in np.eye(3)])
    c = loop.fixman_correction(E.exact_rep(net, beta / 3), beta)
    assert np.abs(c / (0.5 * tau * grad) - 1).max() < 1e-4


def test_the_correction_vanishes_without_dither_and_keeps_the_mass():
    net, base = E.setup()
    rng = np.random.default_rng(3)
    Th = E.exact_rep(net, np.ones(3) / 3) + 0.05 * rng.standard_normal((3, 3))
    beta = np.array([0.8, 1.1, 1.1])
    d0 = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2).drift(Th, beta)
    dz = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0).drift(Th, beta)
    dc = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4).drift(Th, beta)
    assert np.array_equal(d0[1], dz[1]) and np.array_equal(d0[0], dc[0])
    assert abs(dc[1].sum()) < 1e-15 and np.abs(dc[1] - d0[1]).max() > 0


def test_the_slow_limit_cancellation():
    w = np.array([1 / 3, 1 / 3])
    b = Q.drift(w, 1e-5, 0.0)[0]; bc = Q.drift(w, 1e-5, 1.0)[0]
    assert np.abs(bc).max() < 3e-4 * np.abs(b).max()                        # the true ratio at γ = 10⁻⁵ is ≈ 7·10⁻⁵ (referee round 28)


def test_the_gain_selects_the_d_optimal_fit():
    ws = Q.w_sharp()
    f = lambda W, k: (1 + Q.W_MIN[k] ** 2 / W ** 2) / (W - Q.W_MIN[k] ** 2 / W)
    assert abs(ws.sum() - 1) < 1e-12 and np.ptp([f(ws[k], k) for k in range(3)]) < 1e-9     # Lagrange: equal partials of ℓ
    r = Q.b_job((1e-4, 2.0))
    assert r['ier'] == 1 and np.abs(r['zero'] - ws).max() < 1e-4 and np.all(r['eig'].real < 0)


def test_the_density_exponent_near_the_saturation():
    for kF in (0.0, 1.0, 2.0):
        assert abs(Q.c_job((kF, 0.003))['expo'] - (kF - 1)) < 0.03


def test_the_general_form_reduces_to_the_closed_form():
    # the correction's general form at a point of M: (τ/2)[∂_βℓ̂ − ∂_Θℓ̂ (∇²_ΘJ)⁻¹ ∇²_{Θβ}J], ℓ̂(Θ, β) = log det ∇²_ΘJ — the derivative of ℓ
    # along M through the exact fit's implicit function Θ*(β) — computed from the current point only; the closed form is its evaluation
    net, base = E.setup()
    tau = 1e-4
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=1.0, fixman_tau=tau)
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3).ravel()
    J = lambda t, b: loop.J(t.reshape(3, 3), b)
    h = 1e-4

    def hess(t, b):
        I = np.eye(9) * h; H = np.zeros((9, 9))
        for i in range(9):
            for j in range(i, 9):
                H[i, j] = H[j, i] = (J(t + I[i] + I[j], b) - J(t + I[i] - I[j], b) - J(t - I[i] + I[j], b) + J(t - I[i] - I[j], b)) / (4 * h * h)
        return H
    lhat = lambda t, b: np.linalg.slogdet(hess(t, b))[1]
    H = hess(Th, beta)
    hm = 1e-4                                                                 # the mixed Hessian ∇²_{Θβ}J by central differences of ∇_ΘJ
    gradT = lambda t, b: np.array([(J(t + hm * e, b) - J(t - hm * e, b)) / (2 * hm) for e in np.eye(9)])
    Hmix = np.array([(gradT(Th, beta + hm * e) - gradT(Th, beta - hm * e)) / (2 * hm) for e in np.eye(3)]).T     # 9 × 3
    hl = 2e-3
    dl_dT = np.array([(lhat(Th + hl * e, beta) - lhat(Th - hl * e, beta)) / (2 * hl) for e in np.eye(9)])
    dl_dB = np.array([(lhat(Th, beta + hl * e) - lhat(Th, beta - hl * e)) / (2 * hl) for e in np.eye(3)])
    general = 0.5 * tau * (dl_dB - dl_dT @ np.linalg.solve(H, Hmix))
    closed = loop.fixman_correction(Th.reshape(3, 3), beta)
    assert np.abs(general / closed - 1).max() < 2e-3


def test_the_masses_form_equals_the_plug_in_on_the_manifold_and_ignores_the_amplitudes_off_it():
    # the masses form takes W_min,k = |Â_k|/a_max from the classes' contributions Â_k: on M (Â_k = W_k a_k) it is the plug-in form; at fixed
    # masses it does not see an amplitude's excursion, which the plug-in form reads as a mass at its saturation
    net, base = E.setup()
    tau = 1e-4
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3)
    A_hat = beta / 3 * base.features(Th)[0]
    assert np.allclose(A_hat, E.A_BUMP, rtol=1e-12)
    plug = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau)
    mass = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau, fixman_A=A_hat)
    assert np.abs(mass.fixman_correction(Th, beta) / plug.fixman_correction(Th, beta) - 1).max() < 1e-12
    Tx = Th.copy(); Tx[0, 0] += 8.0                                           # the first amplitude far up its plateau, the masses fixed
    assert np.array_equal(mass.fixman_correction(Tx, beta), mass.fixman_correction(Th, beta))
    assert plug.fixman_correction(Tx, beta)[0] > 20 * plug.fixman_correction(Th, beta)[0]
    W_min = np.abs(E.A_BUMP) / net.a_max                                    # the closed form's derivative, (τ/2)∂_βℓ(W)
    W = beta / 3
    assert np.allclose(mass.fixman_correction(Tx, beta), 2.0 * (tau / 3) * (1 + W_min ** 2 / W ** 2) / (W - W_min ** 2 / W), rtol=1e-12)


def test_the_tracked_contributions_stay_put_on_the_manifold():
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, W)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0, fixman_A=E.A_BUMP.copy(), fixman_TA=1.0)
    loop.run_dithered(Th, 3 * W, 20.0, GDD.DT, 0.0, np.random.default_rng(0))
    assert np.abs(loop.fixman_A - E.A_BUMP).max() < 1e-6


def test_the_floor_caps_the_correction_and_its_exponent():
    # referee round 28 (D3aa-2): the denominator is floored at 0.01W_k — within δ_cap = W_min(1/√0.99 − 1) of the saturation the correction is
    # κ_F(τ/m)(1 + q)/(0.01W), and the exponent there is the unfloored one minus κ_F(1 − (1 − W_min²/W²)/0.01); without the floor it is not capped
    net, base = E.setup()
    tau, kF = 1e-4, 2.0
    W = np.array([0.2505, 0.37475, 0.37475])                                  # δ_1 = 0.0005 < δ_cap = 0.00126
    q = (Q.W_MIN / W) ** 2
    A_hat = Q.W_MIN * net.a_max * np.sign(E.A_BUMP)
    for floor in (0.01, 0.0):
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=kF, fixman_tau=tau, fixman_A=A_hat, fixman_floor=floor)
        c = loop.fixman_correction(E.exact_rep(net, W), 3 * W)
        den = np.maximum(W * (1 - q), floor * W)
        assert np.allclose(c, kF * (tau / 3) * (1 + q) / den, rtol=1e-12)
    assert abs(Q.delta_cap(0) - 0.25 * (1 / np.sqrt(0.99) - 1)) < 1e-15 and 0.00125 < Q.delta_cap(0) < 0.00127
    d = 0.0006
    e1, e0 = Q.c_job((kF, d, 0.01))['expo'], Q.c_job((kF, d, 0.0))['expo']
    cut = (1 - (Q.W_MIN[0] / (Q.W_MIN[0] + d)) ** 2) / 0.01
    assert abs(e1 - (e0 - kF * (1 - cut))) < 0.01 and e1 < 0.1 < 0.9 < e0


def test_the_tracked_average_relaxes_at_its_time_constant():
    # the running average Â ← Â + (dt/T_A)(W a − Â) on M without dither: a perturbation decays as (1 − dt/T_A)^n
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35]); Th = E.exact_rep(net, W)
    dA = np.array([0.05, -0.03, 0.02]); TA = 5.0; dt = GDD.DT
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0, fixman_A=E.A_BUMP + dA, fixman_TA=TA)
    loop.run_dithered(Th, 3 * W, 10 * dt, dt, 0.0, np.random.default_rng(0))
    assert np.allclose(loop.fixman_A - E.A_BUMP, dA * (1 - dt / TA) ** 10, rtol=1e-6, atol=1e-12)


def test_kappa_is_the_leaks_projected_stiffness():
    # kappa() against the Jacobian of ∂W·P·F_leak (limit_drift's F₁ with the leak as the perturbation) at w = 1/3
    w = np.ones(2) / 3; h = 1e-4; gamma, eps = 1e-2, 3e-4
    PF = lambda v: Q.drift(v, gamma, 0.0, eps=eps)[0] - Q.drift(v, gamma, 0.0)[0]
    Jl = np.array([(PF(w + h * e) - PF(w - h * e)) / (2 * h) for e in np.eye(2)]).T
    K = Q.kappa(np.ones(3) / 3, eps, gamma)
    assert np.abs(-Jl - K).max() < 1e-5 * np.abs(K).max()


def test_the_long_runs_and_the_paired_analysis_on_synthetic_records(monkeypatch):
    # summarise_long and part_G's statistics on records with known answers (the predictions stubbed out)
    rng = np.random.default_rng(4)
    n, S = 9, 41
    s = np.arange(S, dtype=float) * 1000
    base = np.zeros((n, S, 8)); base[:, :, 0] = s
    base[:, :, 1:4] = 2.0 + 0.1 * rng.standard_normal((n, S, 3))
    W = 1 / 3 + 0.01 * rng.standard_normal((n, S, 3)); W -= (W.sum(axis=2, keepdims=True) - 1) / 3
    base[:, :, 4:7] = W; base[:, :, 7] = 0.025
    sm = Q.summarise_long(base)
    h = S // 2
    assert np.allclose(sm['half'], W[:, h:].mean(axis=1)) and abs(sm['fit2'] - 0.025 ** 2) < 1e-15 and np.all(sm['below'] == 0)
    shift = np.array([0.002, -0.0015, -0.0005])
    R1 = base.copy(); R1[:, :, 4:7] += shift
    monkeypatch.setattr(Q, 'SEEDS_G', tuple(range(n)))
    monkeypatch.setattr(Q, 'drift', lambda w, g, k, floor=0.01, eps=0.0: (np.zeros(2), np.eye(2), None))
    monkeypatch.setattr(Q, 'second_order', lambda g, e, k, start=None, hH=None: dict(W0=np.zeros(3), parts={}, mean=np.full(3, float(k)), std=np.zeros(3),
                                                                                  resid=0.0, scale=1.0))
    out = {r['kF']: r for r in Q.part_G({0.0: base, 1.0: R1})}
    m0 = W[:, h:].mean(axis=1)
    assert np.allclose(out[0.0]['shift'], m0.mean(axis=0) - 1 / 3) and np.allclose(out[0.0]['se'], m0.std(axis=0, ddof=1) / np.sqrt(n))
    assert np.allclose(out[1.0]['diff'], shift, atol=1e-15) and np.all(out[1.0]['dse'] < 1e-15) and np.allclose(out[1.0]['dso'], 1.0)
    assert out[0.0]['esc'] == [] and out[1.0]['esc'] == []
