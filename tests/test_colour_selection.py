"""
Proposition D‴.43 (PARADIGM.md §5d): the colour of the measurement noise selects the kernels. Checked here: det(I + τK) against its
subset (Cauchy–Binet) expansion; the potential's gradient in the relevances D is (4/N)ω; the potential is the trace for white draws and
(1/τ)log det(I + τK) for a slow Markov noise, the step-h filter tending to it; R is linear in D for white draws and concave for coloured
ones; the support count τR − k log(4τ/N) → log(det(Φ_SᵀΦ_S)Π|w_i|); the projection keeps Φw and matches the direct formula, also when
ΦDΦᵀ is singular; the reduced flow stays on M, its P decays at the rate ω exactly, the white flow ends at weighted basis pursuit and a
slow flow, from the same start, at a sparser exact fit; the gradient flow's interpolant against a small-step descent; one step of
label-noise gradient descent against a hand-coded step with AR(1) draws.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import colour_selection as CS                              # noqa: E402


def _design(N, d, seed):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((N, d)) / np.sqrt(N), rng


def test_the_subset_expansion_and_the_gradient():
    Phi, rng = _design(4, 7, 1)
    D = rng.uniform(0.05, 1.5, 7)
    for tau in (0.3, 10.0, 1e3):
        assert abs(CS.det_subsets(Phi, D, tau) / np.linalg.det(np.eye(4) + tau * CS.gram(Phi, D)) - 1) < 1e-12
        e = 1e-6
        g = np.array([(CS.potential(Phi, D + e * np.eye(7)[i], 0, ('ou', tau)) - CS.potential(Phi, D - e * np.eye(7)[i], 0, ('ou', tau)))
                      / (2 * e) for i in range(7)])
        assert np.abs(g - (4 / 4) * CS.weights(Phi, D, 0, ('ou', tau))).max() < 1e-7 * np.abs(g).max()
    h, rho = 0.25, 0.9                                                    # at a step: ∂R/∂D = (4/N)φ_iᵀf(K)φ_i with the scheme's filter
    g = np.array([(CS.potential(Phi, D + 1e-6 * np.eye(7)[i], h, ('ar1', rho)) - CS.potential(Phi, D - 1e-6 * np.eye(7)[i], h, ('ar1', rho)))
                  / 2e-6 for i in range(7)])
    assert np.abs(g - CS.weights(Phi, D, h, ('ar1', rho))).max() < 1e-6 * np.abs(g).max()


def test_the_potentials_white_and_slow():
    Phi, rng = _design(5, 9, 2)
    D = rng.uniform(0.05, 1.5, 9)
    K = CS.gram(Phi, D)
    assert abs(CS.potential(Phi, D, 0, ('fresh',)) - np.trace(K)) < 1e-12
    assert np.allclose(CS.weights(Phi, D, 0, ('fresh',)), (Phi ** 2).sum(axis=0), rtol=1e-12)
    tc = 10.0
    assert abs(CS.potential(Phi, D, 0, ('ou', tc)) - np.log(np.linalg.det(np.eye(5) + tc * K)) / tc) < 1e-12
    h = 1e-4                                                              # finely sampled AR(1) draws → the continuum filter
    assert abs(CS.potential(Phi, D, h, ('ar1', np.exp(-h / tc))) / CS.potential(Phi, D, 0, ('ou', tc)) - 1) < 1e-3
    for _ in range(10):                                                   # R along a segment in D: linear (white), concave (coloured)
        D1, D2 = rng.uniform(0.05, 1.5, 9), rng.uniform(0.05, 1.5, 9)
        sec = lambda corr: sum(c * CS.potential(Phi, D1 + t * (D2 - D1), 0, corr) for c, t in ((1, 0.0), (-2, 0.5), (1, 1.0)))
        assert abs(sec(('fresh',))) < 1e-12 and sec(('ou', tc)) < 0


def test_the_support_count():
    Phi, rng = _design(10, 40, 3)
    w = np.zeros(40); S = [3, 17, 25]; w[S] = (1.0, -0.6, 0.8)
    ref = np.log(np.linalg.det(Phi[:, S].T @ Phi[:, S]) * np.prod(np.abs(w[S])))
    vals = [tc * CS.potential(Phi, np.abs(w), 0, ('ou', tc)) - 3 * np.log(4 * tc / 10) for tc in (1e2, 1e4, 1e6)]
    assert abs(vals[2] - ref) < 1e-4 and abs(vals[2] - ref) < abs(vals[1] - ref) < abs(vals[0] - ref)


def test_the_projection():
    Phi, rng = _design(6, 15, 4)
    D = rng.uniform(0.1, 1.0, 15); q = rng.standard_normal(15)
    direct = q - D * (Phi.T @ np.linalg.solve((Phi * D) @ Phi.T, Phi @ q))
    p = CS.project(Phi, D, q)
    assert np.abs(p - direct).max() < 1e-12 and np.abs(Phi @ p).max() < 1e-12
    assert np.abs(CS.project(Phi, D, p) - p).max() < 1e-12                # a projection
    D[3:] = 1e-200                                                        # three atoms alive: ΦDΦᵀ singular to working precision;
    q = rng.uniform(0.5, 2, 15) * D * rng.uniform(-1, 1, 15)              # in the flow q = ω⊙w with |w_i| ≤ D_i
    p = CS.project(Phi, D, q)
    assert np.all(np.isfinite(p)) and np.abs(Phi @ p).max() < 1e-12 * np.abs(q).max()
    assert np.abs(p[3:]).max() < 1e-150                                   # the dead atoms do not move


def test_the_reduced_flow():
    Phi, rng = _design(10, 40, 0)
    ws = np.zeros(40); idx = rng.choice(40, 4, replace=False); ws[idx] = rng.choice([-1, 1], 4) * rng.uniform(0.5, 1.5, 4)
    y = Phi @ ws
    w0, P0 = CS.flow_interpolant(Phi, y, np.full(40, 0.5), np.full(40, 0.5))
    assert np.abs(Phi @ w0 - y).max() < 1e-12 and np.allclose(P0, 0.25)
    s, W, P = CS.reduced_flow(Phi, w0, P0, 5.0, 0, ('ou', 10.0), t_eval=np.linspace(0, 5.0, 6))
    assert np.abs(Phi @ W - y[:, None]).max() < 1e-7                      # the flow stays on M
    z = np.r_[w0, np.log(P0)]
    om = CS.weights(Phi, np.sqrt(w0 ** 2 + 4 * P0 ** 2), 0, ('ou', 10.0))
    assert np.allclose(CS.rhs(0.0, z, Phi, 0, ('ou', 10.0))[40:], -om)    # d log P/ds = −ω
    c = (Phi ** 2).sum(axis=0)
    wbp = CS.basis_pursuit(Phi, y, c)
    _, Ww, _ = CS.reduced_flow(Phi, w0, P0, 2000.0, 0, ('fresh',))
    assert np.abs(Ww[:, -1] - wbp).max() < 1e-5                          # white: weighted basis pursuit (here not the truth)
    assert np.abs(wbp - ws).max() > 0.5
    _, Wt, _ = CS.reduced_flow(Phi, w0, P0, 400.0 + 60 * np.sqrt(101.0), 0, ('ou', 100.0))
    assert np.abs(Wt[:, -1] - ws).max() < 1e-6                            # slow: the 4-sparse truth
    assert (np.abs(Wt[:, -1]) > 1e-3).sum() == 4 < (np.abs(wbp) > 1e-3).sum()


def test_the_flow_interpolant_and_one_noisy_step():
    Phi, rng = _design(8, 20, 5)
    y = Phi @ rng.standard_normal(20)
    u, v = rng.uniform(0.2, 1, 20), rng.uniform(0.2, 1, 20)
    wf, Pf = CS.flow_interpolant(Phi, y, u, v)
    wg, _ = CS.gd_interpolant(Phi, y, u=u, v=v, h=0.002)
    assert np.abs(Phi @ wf - y).max() < 1e-12 and np.allclose(Pf, u * v) and np.abs(wf - wg).max() < 2e-4

    class Rec:
        def __init__(self):
            self.r, self.out = np.random.default_rng(9), []

        def standard_normal(self, shape):
            self.out.append(self.r.standard_normal(shape)); return self.out[-1]
    h, sig, rho, rec = 0.1, 0.3, 0.8, Rec()
    u1, v1, _ = CS.label_noise_gd(Phi, y, u, v, h, sig, rho, 2, rec, 3)
    xi = sig * rec.out[0]; uu, vv = np.tile(u, (3, 1)), np.tile(v, (3, 1))
    for k in range(2):
        g = ((uu * uu - vv * vv) @ Phi.T - y[None, :] - xi) @ Phi / 8
        uu, vv = uu - 2 * h * g * uu, vv + 2 * h * g * vv
        xi = rho * xi + np.sqrt(1 - rho * rho) * sig * rec.out[k + 1]
    assert np.allclose(u1, uu, rtol=0, atol=1e-14) and np.allclose(v1, vv, rtol=0, atol=1e-14)


def test_the_kernel_dictionary():
    Phi = CS.rbf_dictionary(np.linspace(0, 1, 15), np.linspace(0, 1, 40), 0.05)
    assert Phi.shape == (15, 40) and np.allclose(np.linalg.norm(Phi, axis=0), 1.0)
