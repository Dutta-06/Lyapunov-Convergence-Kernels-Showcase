"""
Proposition D‴.39 (PARADIGM.md §5d): dissipativity of the dithered gauge. Checked here: the dead zone V = (σ_B/2)Σ_i‖ϑ_i − Π_B(ϑ_i)‖²_{G⁻¹}
leaves the drift and the energy unchanged inside the box and adds minus its G-gradient, −σ_B(ϑ_i − Π_B(ϑ_i)), outside; J is bounded on
the unconstrained coordinates (so e^{−J/τ} is not integrable) and J + V grows like (σ_B/2)dist²; the box's wall lies inside the saturation,
W♯ is representable in it, and the polytope's vertices have two classes at their walls.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402

BOX = (5.0, 3.0, 4.0)


def test_the_dead_zone_is_inactive_inside_the_box():
    net, base = E.setup()
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3) + 0.05 * np.random.default_rng(1).standard_normal((3, 3))
    assert np.all(np.abs(Th) < np.array(BOX)[None, :])
    free = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    boxed = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=BOX, box_sigma=0.3)
    assert np.array_equal(free.drift(Th, beta)[0], boxed.drift(Th, beta)[0]) and free.J(Th, beta) == boxed.J(Th, beta)


def test_the_dead_zone_force_is_minus_the_g_gradient_of_its_energy():
    net, base = E.setup()
    beta = 3 * np.array([0.30, 0.35, 0.35])
    Th = E.exact_rep(net, beta / 3)
    Th[0, 0] = 7.5; Th[2, 1] = -3.4; Th[1, 2] = 4.6                          # three coordinates beyond the box
    free = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    boxed = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=BOX, box_sigma=0.3)
    dV = boxed.drift(Th, beta)[0] - free.drift(Th, beta)[0]
    V = lambda T: boxed.J(T, beta) - free.J(T, beta)
    h = 1e-5
    grad = np.array([(V(Th + h * e.reshape(3, 3)) - V(Th - h * e.reshape(3, 3))) / (2 * h) for e in np.eye(9)]).reshape(3, 3)
    assert np.allclose(dV, -E.G[None, :] * grad, atol=1e-8)
    exc = np.zeros((3, 3)); exc[0, 0] = 2.5; exc[2, 1] = -0.4; exc[1, 2] = 0.6
    assert np.allclose(dV, -0.3 * exc, atol=1e-12)


def test_the_energy_is_bounded_on_the_unconstrained_coordinates_and_the_dead_zone_confines():
    # (i): the physical coordinates are bounded and the cohesion energy W ≤ R²/8, so J is bounded — e^{−J/τ} is not integrable; with the dead
    # zone, J + V ≥ (σ_B/2)dist_{G⁻¹}(Θ, B)² − sup|J|
    net, base = E.setup()
    beta = 3 * np.array([0.30, 0.35, 0.35])
    free = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    boxed = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=BOX, box_sigma=0.3)
    m = 3
    sup_r2 = float(np.mean((np.abs(base.F) + net.a_max) ** 2))                # |f̂| ≤ (1/m)Σβ_i a_max = a_max
    bound = 0.5 * m * sup_r2 + 0.5 * free.sigma0 * m * free.R ** 2 / 8
    rng = np.random.default_rng(5)
    for scale in (1.0, 10.0, 1e3, 1e6):
        for _ in range(50):
            Th = scale * rng.standard_normal((3, 3))
            J = free.J(Th, beta)
            assert 0.0 <= J <= bound
            d2 = float(np.sum(boxed.box_excess(Th) ** 2 / E.G[None, :]))
            assert boxed.J(Th, beta) >= 0.15 * d2 - 1e-12


def test_the_box_keeps_the_best_conditioned_fit_and_moves_the_wall_inside_the_saturation():
    import exp_gauge_fixman as Q
    import exp_gauge_plateau as P
    ws = Q.w_sharp()
    assert np.all(P.W_BOX > P.W_MIN) and np.all(ws > P.W_BOX)                # the wall inside the saturation, W♯ representable in the box
    a_max = E.setup()[0].a_max
    assert np.all(np.abs(E.A_BUMP) / ws < a_max * np.tanh(BOX[0] / a_max))   # its amplitudes inside the box
    V = P.vertices()
    assert np.allclose(V.sum(axis=1), 1) and np.all(V >= P.W_BOX[None, :] - 1e-15)
    assert np.allclose(np.sort(V, axis=None)[:6], np.sort(np.r_[P.W_BOX, P.W_BOX]))
