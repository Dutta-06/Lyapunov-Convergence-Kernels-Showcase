"""Corollary D‴.26 (PARADIGM.md §5d): perfect tracking with an adaptive prior in the real closed loop. At the lift of an exact
representation (here the matching one, w = q = 1/3, m = 6) with the prior at its atoms and the plant on the reference, one period returns
the state; the Floquet exponents of the lift are σ₀/k_e exactly on the (m − K₀)p separations and the weighted K₀-atom loop's on the
collapsed subspace; and the K₀-atom loop's are the averaged prediction (1/k_e)λ_j^± of Theorem D‴.23(ii) up to the averaging correction."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import ot_leak as OT
from src.core import adaptive_prior as AP
from src.core import adaptive_floquet as AF
from src.core.periodic import floquet_exponents

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SIGMA0 = 0.15 * 0.047162567584665846
F_VEC = lambda X: sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
XR = lambda t: np.array([2.0 * np.sin(t)])
XRD = lambda t: np.array([2.0 * np.cos(t)])
KE = 25.0


def test_perfect_tracking_has_the_separation_exponent_and_the_K0_atom_spectrum_for_every_m():
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
    assert np.abs(out[3] / pred - 1).max() < 0.2                                        # the averaged prediction, up to the O(1/k_e) correction
