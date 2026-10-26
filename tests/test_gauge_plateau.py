"""
Proposition D‴.39 (PARADIGM.md §5d): dissipativity of the dithered gauge. Checked here: the dead zone V = (σ_B/2)Σ_i‖ϑ_i − Π_B(ϑ_i)‖²_{G⁻¹}
leaves the drift and the energy unchanged inside the box and adds minus its G-gradient, −σ_B(ϑ_i − Π_B(ϑ_i)), outside.
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
