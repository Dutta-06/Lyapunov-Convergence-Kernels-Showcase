"""
Proposition D‴.38 (PARADIGM.md §5d): a Fixman correction for the replicator under dither. Checked here: the loop's correction equals
(τ/2)∂_βℓ along the manifold of exact fits (ℓ = log det ∇²_ΘJ, by finite differences of the Hessian's determinant); it vanishes without
dither and keeps Σβ = m; in the slow-replicator limit it cancels the dither's drift; with κ_F = 2 the corrected field's zero is the unique
maximiser W♯ of ℓ; near a saturation the one-dimensional zero-flux density behaves like δ^{κ_F − 1}; the correction's general form,
through the exact fit's implicit function and computed from the current point only, reduces to the closed form.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
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
    assert np.abs(bc).max() < 2e-3 * np.abs(b).max()


def test_the_gain_selects_the_best_conditioned_fit():
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
