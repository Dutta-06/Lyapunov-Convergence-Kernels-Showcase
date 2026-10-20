"""Corollary D‴.26 and Lemma D‴.26a (PARADIGM.md §5d): perfect tracking with an adaptive prior in the real closed loop. At the lift of an exact
representation with the prior at its atoms and the plant on the reference, one period returns the state; the Floquet exponents of the lift
are σ₀/k_e exactly on the (m − K₀)p separations and the weighted K₀-atom loop's on the collapsed subspace; and the K₀-atom loop's are the
averaged prediction (1/k_e)λ_j^± of Theorem D‴.23(ii) up to a correction of second order in 1/k_e (slow time), explicit to that order with
a remainder of fourth order. Checked at the matching weights (η = ∞, m = 6) and at the unequal weights (2, 1, 1)/4 with the prior in the
state (η = σ₀, m = 4): the monodromy by central differences of the period map, the real loop's linearisation at the lift, and the
variational equation on the collapsed subspace (src/core/adaptive_floquet.CollapsedVariational)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import ot_leak as OT
from src.core import adaptive_prior as AP
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN
from src.core.periodic import floquet_exponents

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SIGMA0 = 0.15 * 0.047162567584665846
F_VEC = lambda X: sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
XR = lambda t: np.array([2.0 * np.sin(t)])
XRD = lambda t: np.array([2.0 * np.cos(t)])
KE = 25.0


def test_perfect_tracking_has_the_separation_exponent_and_the_K0_atom_spectrum_at_m_6():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    exact = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    q = np.ones(3) / 3
    out = {}
    for m, labels, weights in ((3, np.arange(3), q), (6, np.repeat(np.arange(3), 2), None)):
        cl = AF.AdaptiveRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, labels, 3, np.inf, weights=weights)
        z = cl.lift(exact)
        assert np.abs(AF.period_map(cl, z, m) - z).max() < 1e-12                             # perfect tracking is T₀-periodic
        fast, slow = floquet_exponents(AF.monodromy(cl, z, m), PER, 1)
        out[m] = np.sort(slow)
    s6 = out[6]
    order = np.argsort(np.abs(s6 - SIGMA0 / KE))
    sep, coll = s6[order[:9]], np.sort(s6[order[9:]])
    # the monodromy by central differences (step 10⁻⁵, rtol 10⁻¹³) resolves the multipliers to ~10⁻¹⁰, the slowest rate to ~10⁻⁵ relative
    assert np.abs(sep * KE / SIGMA0 - 1).max() < 1e-6                                   # the (m − K₀)p separations at σ₀/k_e (exactly, in theory)
    assert np.abs(coll / out[3] - 1).max() < 1e-4                                       # the collective exponents are the K₀-atom loop's
    loop = OT.OTPriorLoop(net, F_VEC, XR, PER, G, SIGMA0, exact, n_quad=256)
    h = AP.gauss_newton_block(loop, exact, q)[1]
    pred = AF.predicted_slow_exponents(h, SIGMA0, np.inf, KE, 3, 3, 3)
    assert np.abs(out[3] / pred - 1).max() < 0.2                                        # the averaged prediction, up to the O(k_e⁻²) slow-time correction (0.137)


W211 = np.array([2, 1, 1]) / 4                                                         # representable: the amplitudes A_k/w_k = 2, −2.8, 3.2


def _exact(net, w):
    return net.unconstrained(np.array([[A / wk] for (A, C, S), wk in zip(BUMPS, w)]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def test_the_real_loops_linearisation_at_the_lift_splits_into_separations_and_the_collapsed_variational_equation():
    """At finite η with the prior in the state and unequal weights (m = 4, classes (0, 0, 1, 2)): the central-difference Jacobian of the real
    loop's field at the lift maps a separation v to −(σ₀/k_e)v and the collapsed subspace into itself, where it is the variational equation of
    CollapsedVariational (in the whitened weighted coordinates), at every t — (a) and (b) of the corollary at the level of the generator."""
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    e = _exact(net, W211)
    labels = np.array([0, 0, 1, 2]); m, p = 4, 3
    cl = AF.AdaptiveRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, labels, 3, SIGMA0)
    cv = AF.CollapsedVariational(net, XR, PER, G, SIGMA0, SIGMA0, e, W211)
    Dw = (np.sqrt(W211)[:, None] / np.sqrt(G)[None, :]).ravel()
    P = np.zeros((1 + m * p + 3 * p, 1 + 3 * p + 3 * p))                             # the collapsed coordinates (δe, δϑ_k, δz_k) → the full state
    P[0, 0] = 1.0
    for i, k in enumerate(labels):
        P[1 + i * p:1 + (i + 1) * p, 1 + k * p:1 + (k + 1) * p] = np.eye(p)
    P[1 + m * p:, 1 + 3 * p:] = np.eye(3 * p)
    Tw = np.concatenate([[1.0], Dw, Dw])
    for t in (0.0, 1.3, 2.9, 4.4):
        z = cl.lift(e); z[0] = XR(t)[0]
        Dz = len(z); Jr = np.zeros((Dz, Dz))
        for j in range(Dz):
            h = 1e-6 * max(1.0, abs(z[j])); dz = np.zeros(Dz); dz[j] = h
            Jr[:, j] = (cl.rhs(t, z + dz, m) - cl.rhs(t, z - dz, m)) / (2 * h)
        for c in range(p):                                                          # (a): class 0's separations decay at σ₀/k_e
            v = np.zeros(Dz); v[1 + c] = 1.0; v[1 + p + c] = -1.0
            assert np.abs(Jr @ v + (SIGMA0 / KE) * v).max() < 1e-9
        Mk = (Tw[:, None] ** -1) * cv.generator(t, KE) * Tw[None, :]                # the collapsed generator in the original coordinates
        assert np.abs(Jr @ P - P @ Mk).max() < 1e-8 * np.abs(Mk).max()              # (b): invariant, and the variational equation


def test_at_finite_eta_and_unequal_weights_the_lift_has_the_separation_exponent_and_the_K0_atom_spectrum_at_m_4():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    e = _exact(net, W211)
    cl = AF.AdaptiveRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, np.array([0, 0, 1, 2]), 3, SIGMA0)
    z = cl.lift(e)
    assert np.abs(AF.period_map(cl, z, 4) - z).max() < 1e-12                             # perfect tracking is T₀-periodic
    fast, slow = floquet_exponents(AF.monodromy(cl, z, 4), PER, 1)
    slow = np.sort(slow)
    order = np.argsort(np.abs(slow - SIGMA0 / KE))
    sep, coll = slow[order[:3]], np.sort(slow[order[3:]])
    assert len(coll) == 18                                                              # 9 collective pairs: the prior is in the state
    assert np.abs(sep * KE / SIGMA0 - 1).max() < 1e-6                                   # the (m − K₀)p = 3 separations at σ₀/k_e
    cv = AF.CollapsedVariational(net, XR, PER, G, SIGMA0, SIGMA0, e, W211)
    var, err = cv.exponents(KE)
    assert np.abs(coll / var - 1).max() < 1e-4                                          # the collective exponents: the collapsed variational ones
    gens = cv.averaged()
    p1, p2 = cv.predicted(KE, 1, gens=gens), cv.predicted(KE, 2, gens=gens)
    assert np.abs(var / p1 - 1).max() < 0.27                                            # λ_j^±(w)/k_e up to the O(k_e⁻²) correction (0.133 here)
    assert np.max(np.abs(var - p2) / p1) < 0.05                                         # the second-order prediction (0.026 here)
    assert np.all(var < p1)        # an instance property (every c_j < 0 here, the filter's part dominating); the sign of c_j is not fixed in general


def test_the_collective_exponents_are_even_in_one_over_ke_with_the_explicit_second_order_coefficient():
    """An asymmetric reference, where the first-order averaging correction Ā₁ does not vanish: Ā₁ is skew in the metric diag(I, (σ₀/η)I),
    the residual of the second-order prediction falls like k_e⁻⁴ (relative: the odd orders vanish), the measured second-order coefficient of
    the slowest mode is the explicit c₁, and the error channel's exponent (Liouville) is k_e − tr H/k_e up to O(k_e⁻³)."""
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    e = _exact(net, W211)
    xr = lambda t: np.array([2.0 * np.sin(t) + 0.8 * np.sin(2.0 * t)])
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, SIGMA0, e, W211)
    A0, A1, A2, H, filt = gens = cv.averaged()
    Wm = np.diag(np.r_[np.ones(9), np.ones(9)])                                         # η = σ₀: the metric diag(I, (σ₀/η)I) = I
    assert np.abs(A1).max() > 1e-3 and np.abs(Wm @ A1 + (Wm @ A1).T).max() < 1e-13      # non-zero and skew
    lam, c, cf = cv.second_order(gens=gens)
    res, meas = [], []
    for ke in (100.0, 200.0):
        var, err = cv.exponents(ke)
        p1, p2 = cv.predicted(ke, 1, gens=gens), cv.predicted(ke, 2, gens=gens)
        res.append(np.max(np.abs(var - p2) / p1))
        meas.append((var[0] * ke - lam[0]) * ke ** 2)
        assert abs(err - (ke - np.trace(H) / ke)) < 2e-4                               # the error channel by Liouville's formula
    assert res[0] / res[1] > 12 and res[1] < 1e-4                                       # O(k_e⁻⁴): the third order vanishes as well
    assert abs(meas[1] / c[0] - 1) < 0.01 and c[0] < 0 and cf[0] < 0                    # c₁ explicit (negative here), with the filter's part
    # the mechanism of the vanishing first order: ℓ_jĀ₁r_j = 0 for every mode, Ā₁ being skew in the metric of Ā₀'s self-adjointness; at
    # η = 10σ₀ that metric diag(I, (σ₀/η)I) is not the identity, and Ā₀ is self-adjoint in it only
    lv, Rv = np.linalg.eig(A0); Lv = np.linalg.inv(Rv)
    assert np.abs(np.diag(Lv @ A1 @ Rv)).max() < 1e-10 * np.abs(A1).max()
    cv10 = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, 10 * SIGMA0, e, W211)
    B0, B1 = cv10.averaged()[:2]
    W10 = np.diag(np.r_[np.ones(9), np.full(9, 0.1)])
    assert np.abs(W10 @ B0 - (W10 @ B0).T).max() < 1e-12 * np.abs(B0).max() and np.abs(B0 - B0.T).max() > 0.1 * np.abs(B0).max()
    assert np.abs(W10 @ B1 + (W10 @ B1).T).max() < 1e-12 and np.abs(B1).max() > 1e-3


FPRIME = lambda x: sum(-A * (x - C) / S ** 2 * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def _asym_cv(eta):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    e = _exact(net, W211)
    xr = lambda t: np.array([2.0 * np.sin(t) + 0.8 * np.sin(2.0 * t)])
    return AF.CollapsedVariational(net, xr, PER, G, SIGMA0, eta, e, W211)


def test_chang_slow_subspace_carries_exactly_the_slow_multipliers():
    """Lemma D‴.26a(i)–(ii): the contraction gives the periodic P with an invariance residual at rounding; the slow block's exponents are the
    full variational equation's (to the monodromies' accuracy), and P − P₂ scales like ε³ (the expansion of (i))."""
    cv = _asym_cv(SIGMA0)
    t, K, Lam, om = cv.chang_grid(512)
    der = lambda X, k=1: np.real(np.fft.ifft(((1j * om) ** k).reshape(-1, 1, 1) * np.fft.fft(X, axis=0), axis=0))
    tr = []
    for ke in (25.0, 50.0):
        eps = 1.0 / ke
        P, A, it, res = cv.chang(eps)
        assert res < 1e-12 and it < 60
        rates = np.sort(np.linalg.eigvals(-cv.slow_generator(eps)[0]).real)
        full = np.sort(cv.exponents(ke)[0]) * ke
        assert np.abs(rates / full - 1).max() < 1e-8
        P2 = -K + eps * der(K) - eps ** 2 * (der(K, 2) + np.einsum('tij,tkj,tkl->til', K, K, K) + K @ Lam)
        tr.append(np.abs(P - P2).max())
    assert 6 < tr[0] / tr[1] < 10                                                     # O(ε³)


def test_the_slow_rates_are_exactly_even_in_epsilon_and_a_mismatch_breaks_the_parity():
    """Lemma D‴.26a(iv): the slow block at −ε (the same equation for P) has the slow-time rates of the slow block at ε — exactly, while their
    deviation from the first-order prediction is 12%; with the plant's slope as a mismatch in the error channel the equality fails at first
    order. At η = ∞ as well."""
    for eta in (SIGMA0, np.inf):
        cv = _asym_cv(eta)
        A0 = cv.averaged()[0]
        nu0 = np.sort(np.linalg.eigvals(-A0).real)
        eps = 1.0 / 25.0
        rp = np.sort(np.linalg.eigvals(-cv.slow_generator(eps)[0]).real)
        rm = np.sort(np.linalg.eigvals(-cv.slow_generator(-eps)[0]).real)
        assert np.abs(rp - rm).max() < 1e-10 * np.abs(rp).max() and np.abs(rp - nu0).max() > 1e-2 * np.abs(nu0).max()
        tg = cv.chang_grid(512)[0]
        a = np.array([FPRIME(float(cv.xr(s)[0])) for s in tg])
        mp = np.sort(np.linalg.eigvals(-cv.slow_generator(eps, mismatch=a)[0]).real)
        mm = np.sort(np.linalg.eigvals(-cv.slow_generator(-eps, mismatch=a)[0]).real)
        assert np.abs(mp - mm).max() > 1e-2 * np.abs(mp).max()


def test_the_slow_noise_is_decoupled_exactly_and_high_passed():
    """Lemma D‴.26a(v): Y′ = Y + εSv is decoupled from v (the residual at rounding), its noise b → −(K̇ + aK)ᵀ at O(ε) (halving per doubling of
    k_e), and the exact periodic covariance at t = 0 has Cov(v) → N₀/(2ε³), E[vY′ᵀ] → −N₀bᵀ and Cov(Y′) → Σ̄ (the averaged slow covariance)."""
    cv = _asym_cv(SIGMA0)
    dev = []
    for ke in (50.0, 100.0):
        P, S, b, b0, coup = MN.noise_decoupling(cv, FPRIME, 1.0 / ke)
        assert coup < 1e-14
        dev.append(np.abs(b - b0).max() / np.abs(b0).max())
    assert 1.7 < dev[0] / dev[1] < 2.3 and dev[1] < 0.15
    o = MN.decoupled_covariance(cv, FPRIME, 100.0)
    assert abs(o['Vv'][0, 0] * 2e-6 - 1) < 1e-2
    assert np.abs(o['C'] + o['b'].T).max() < 3e-2 * np.abs(o['b']).max()
    assert np.abs(o['CYp'] - o['Sigma']).max() < 3e-2 * np.abs(o['Sigma']).max()

