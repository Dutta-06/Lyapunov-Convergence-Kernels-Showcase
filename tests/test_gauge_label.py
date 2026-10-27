"""
Proposition D‴.41 (PARADIGM.md §5d): measurement noise on the gauge. Checked here: the measured error's noise enters the loop along the
fields U_x of `label_noise.fields` (the drift with r + ξ minus the drift with r is (1/N)Σ_xξ(x)U_x, on the positions and the log-masses);
on M the fields are orthogonal to the gauge's tangent in the flow's metric (no diffusion along M); the Itô drift and the Wong–Zakai term
cancel at a point of M (the Stratonovich drift vanishes); the closed form's stationary point W_ℓ minimises Ψ_ℓ and is interior; the AR(1)
noise at ρ = 0 is the Euler scheme bit for bit, and Heun's scheme without noise is the explicit trapezoidal rule.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import label_noise as LN                                   # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def _loop(gamma=1e-2):
    net, base = E.setup()
    return net, base, MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)


def test_the_noise_enters_along_the_fields():
    net, base, loop = _loop()
    W = np.array([0.40, 0.25, 0.35]); b = 3 * W
    Th = E.exact_rep(net, W) + 0.03 * np.random.default_rng(0).standard_normal((3, 3))
    xi = np.random.default_rng(1).standard_normal(base.N)
    dT0, db0 = loop.drift(Th, b)
    dT1, db1 = loop.drift(Th, b, xi)
    UT, Ul = LN.fields(loop, Th, b)
    assert np.allclose(dT1 - dT0, np.einsum('x,xij->ij', xi, UT) / base.N, rtol=1e-12, atol=1e-15)
    assert np.allclose((db1 - db0) / b, xi @ Ul / base.N, rtol=1e-12, atol=1e-15)


def test_the_fields_are_orthogonal_to_the_gauge_in_the_flow_metric():
    # t ∈ T: δα_k = −(a_k/a′_k)δℓ_k, Σβ_kδℓ_k = 0; the metric Σ‖δϑ‖²_{G⁻¹} + Σ(β_k/γ)δℓ_k²
    net, base, loop = _loop()
    W = np.array([0.30, 0.35, 0.35]); b = 3 * W
    Th = E.exact_rep(net, W)
    a = base.features(Th)[0]; ap = 1 - (a / net.a_max) ** 2
    UT, Ul = LN.fields(loop, Th, b)
    for dl in (np.array([1.0, -1.0, 0.0]), np.array([0.5, 0.5, -1.0])):
        dl = dl - (b @ dl) / b.sum()
        t_alpha = -(a / ap) * dl
        inner = UT[:, :, 0] @ (t_alpha / E.G[0]) + Ul @ (b / loop.gamma * dl)
        scale = np.abs(UT[:, :, 0]).max() * np.abs(t_alpha).max() / E.G[0]
        assert np.abs(inner).max() < 1e-12 * scale


def test_the_wong_zakai_term_cancels_the_ito_drift():
    import exp_gauge_label as GL
    r = GL.a_job.__wrapped__((1e-2, (1 / 3, 1 / 3)))                     # uncached
    assert np.abs(r['ito']).max() > 1e-6
    assert np.abs(r['ito'] + r['wz']).max() < 1e-6 * np.abs(r['ito']).max()
    assert np.abs(r['D']).max() < 1e-12 * np.abs(r['ito']).max()


def test_the_closed_forms_target_is_interior_and_minimises_psi():
    net, base, loop = _loop()
    Th = E.exact_rep(net, np.ones(3) / 3)
    P = E.G[0] * np.mean(base.features(Th)[4] ** 2, axis=0)
    Wm = np.abs(E.A_BUMP) / net.a_max
    Wl = LN.w_label(P, Wm)
    assert abs(Wl.sum() - 1) < 1e-10 and np.all(Wl > Wm + 0.05)
    u = LN.selection(loop, E.exact_rep(net, Wl), 3 * Wl)
    assert np.ptp(u) < 1e-9 * np.abs(u).max()
    rng = np.random.default_rng(2)
    for _ in range(200):
        W = Wl + 0.05 * rng.standard_normal(3); W = W - (W.sum() - 1) / 3
        if np.all(W > Wm):
            assert LN.psi(W, P, Wm) >= LN.psi(Wl, P, Wm)
    # the closed form is −∂_βΨ_ℓ along M (Ψ_ℓ in β through W = β/3)
    W = np.array([0.40, 0.25, 0.35]); h = 1e-6
    u = LN.selection(loop, E.exact_rep(net, W), 3 * W)
    grad = np.array([(LN.psi(W + h * e / 3, P, Wm) - LN.psi(W - h * e / 3, P, Wm)) / (2 * h) for e in np.eye(3)])
    assert np.allclose(u, -grad, rtol=1e-6)


def test_the_schemes_reduce_as_they_should():
    net, base, loop = _loop()
    W = np.array([0.40, 0.25, 0.35]); T0 = E.exact_rep(net, W)
    run = lambda **kw: loop.run_dithered(T0, 3 * W, 2.5, 0.25, 0.0, np.random.default_rng(5), **kw)
    T_ito, b_ito = run(label_v=1e-4, label_scheme='ito')
    T_ar, b_ar = run(label_v=1e-4, label_scheme='ar1', label_rho=0.0)
    assert np.array_equal(T_ito, T_ar) and np.array_equal(b_ito, b_ar)
    T_h, b_h = run(label_v=1e-30, label_scheme='heun')                     # no noise to speak of: the explicit trapezoidal rule
    Th, lb = T0.copy(), np.log(3 * W)
    for _ in range(10):
        d1, g1 = loop.drift(Th, np.exp(lb))
        Tp, lp = Th + 0.25 * d1, lb + 0.25 * g1 / np.exp(lb)
        d2, g2 = loop.drift(Tp, np.exp(lp))
        Th, lb = Th + 0.125 * (d1 + d2), lb + 0.125 * (g1 / np.exp(lb) + g2 / np.exp(lp))
        lb = lb - np.log(np.exp(lb).sum() / 3)
    assert np.allclose(T_h, Th, rtol=0, atol=1e-12) and np.allclose(b_h, np.exp(lb), rtol=0, atol=1e-12)
    assert LN.ar1_factor(0.0) == 1.0 and abs(LN.ar1_factor(1 / 3) - 0.5) < 1e-15 and abs(LN.ar1_factor(-1 / 3) - 2.0) < 1e-15
