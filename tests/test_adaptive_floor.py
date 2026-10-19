"""Proposition D‴.24 (PARADIGM.md §5d): the finite-size floor is the leak's bias at the sampled weights, and an adaptive prior removes it
without a mass constraint. The nearest-atom leak toward an adaptive prior is the gradient flow of L̄ + (σ₀/2)(1/m)Σ_i min_j‖ϑ_i − z_j‖²_{G⁻¹}
in the metric with the cell masses, and that energy is the minimum of the OT energy over the prior's masses (semi-relaxed transport); at a
matching prior the function-level weight-sensitivity of Proposition D‴.9 is the weight-derivative of the leak's bias, Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰,
its floor law the sum of the mode terms h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j, below the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ); and the exact representation of
sampled weights, with the prior at its atoms, is stationary for the adaptive prior while the fixed prior's stationary point is biased; and
the fixed prior's weight sensitivity is H_L(H_L + σ₀I)⁻¹ times the gauge's, the rest being the bias."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import ot_leak as OT
from src.core import prior_design as PD
from src.core import adaptive_prior as AP
from src.core import floor_bias as FB

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
F_VEC = lambda X: sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
X_REF = lambda t: np.array([2.0 * np.sin(t)])
SIGMA0 = 0.15 * 0.047162567584665846
Q = np.ones(3) / 3


def _setup(shift=0.0, n_quad=96):
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    rep = lambda sh: net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + sh] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    exact, prior = rep(0.0), rep(shift)
    loop = OT.OTPriorLoop(net, F_VEC, X_REF, PER, G, SIGMA0, prior, n_quad=n_quad)
    return loop, exact, prior


def test_the_nearest_atom_adaptive_flow_is_the_gradient_flow_of_the_free_mass_energy():
    loop, exact, prior = _setup(0.15)
    eta = 2 * SIGMA0
    apl = AP.AdaptivePriorLoop(loop, eta, plan='nearest')
    rng = np.random.default_rng(3)
    counts = np.array([5, 4, 3]); m = int(counts.sum())
    Th = np.repeat(exact, counts, axis=0) + 0.05 * rng.standard_normal((m, 3)) * SQG[None, :]
    Z = prior + 0.02 * rng.standard_normal((3, 3))
    gamma, _ = apl.plan(Th, Z)
    w = gamma.sum(axis=0)
    assert np.allclose(w, counts / m)                                         # the cell masses are the sampled weights
    dTh, dZ = apl.drift(Th, Z, gamma)
    J = lambda T_, Z_: apl.energy(T_, Z_)
    h = 1e-6; gT = np.zeros_like(Th); gZ = np.zeros_like(Z)
    for i in range(m):
        for c in range(3):
            E = np.zeros_like(Th); E[i, c] = h
            gT[i, c] = (J(Th + E, Z) - J(Th - E, Z)) / (2 * h)
    for j in range(3):
        for c in range(3):
            E = np.zeros_like(Z); E[j, c] = h
            gZ[j, c] = (J(Th, Z + E) - J(Th, Z - E)) / (2 * h)
    assert np.abs(dTh + m * G[None, :] * gT).max() < 1e-7                      # dϑ_i/ds = −mG∇_{ϑ_i}𝒥_∘
    assert np.abs(dZ + (eta / (SIGMA0 * w[:, None])) * G[None, :] * gZ).max() < 1e-7   # dz_j/ds = −(η/(σ₀w_j))G∇_{z_j}𝒥_∘, w the cell masses


def test_the_free_mass_energy_is_the_minimum_of_the_OT_energy_over_the_masses():
    loop, exact, prior = _setup(0.0)
    rng = np.random.default_rng(5)
    m = 12
    for trial in range(4):
        Th = exact[rng.integers(0, 3, m)] + 1.2 * rng.standard_normal((m, 3)) * SQG[None, :]
        Z = exact + 0.3 * rng.standard_normal((3, 3))
        gN, dist = AP.nearest_plan(Th, Z, G)
        wN = gN.sum(axis=0)
        if np.any(wN == 0):
            continue
        assert abs(AP.plan_to(Th, Z, wN, G)[1] - dist) < 1e-12                 # attained at the cell masses …
        for _ in range(5):
            qq = rng.dirichlet(np.ones(3))
            assert AP.plan_to(Th, Z, qq, G)[1] >= dist - 1e-12               # … and a minimum over the masses


def test_the_floor_law_is_the_weight_derivative_of_the_leaks_bias():
    loop, exact, prior = _setup(0.0)
    out = FB.bias_floor_law(loop, Q, exact, SIGMA0)
    S = PD.weight_sensitivity_prior(loop.nearest, Q, exact)
    F2 = PD.function_sensitivity(loop.nearest, Q, exact, S) @ np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
    assert np.linalg.norm(out['F_bias'] - F2) < 1e-5 * np.linalg.norm(F2)       # F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰
    numeric = float(np.trace(F2 @ out['C2'] @ F2.T)) / F2.shape[0]
    assert abs(out['law'] - numeric) < 1e-5 * numeric                         # the mode terms sum to the law
    assert out['law'] < out['bound']                                           # each mode's weight ≤ σ₀/4
    share = out['terms'] / out['terms'].sum()
    near = np.abs(np.log(out['h'] / SIGMA0)) < np.log(5.0)                     # the modes with h_j within a factor 5 of σ₀ …
    assert share[near].sum() > 0.8                                            # … carry most of the floor


def test_the_exact_representation_of_sampled_weights_is_stationary_for_the_adaptive_prior():
    loop, exact, prior = _setup(0.0)
    counts = np.array([6, 5, 4]); m = int(counts.sum()); w = counts / m
    te, res, fit = FB.exact_weighted_atoms(loop, w, exact)
    assert res < 1e-12 and fit < 1e-12
    assert np.allclose(loop.features(te)[0], np.array([A for A, C, S in BUMPS]) / w)   # the amplitudes A_k/w_k (the gauge)
    Th = np.repeat(te, counts, axis=0)
    apl = AP.AdaptivePriorLoop(loop, SIGMA0, plan='nearest')
    gamma, _ = apl.plan(Th, te)
    dTh, dZ = apl.drift(Th, te, gamma)
    assert np.abs(dTh).max() < 1e-12 and np.abs(dZ).max() < 1e-12 and loop.fit_error(Th) < 1e-12
    Tw, resw = PD.weighted_stationary_atoms_prior(loop.nearest, w, exact, w0=Q, n_cont=4)   # the fixed (matching) prior at the same weights
    assert resw < 1e-12 and loop.fit_error(np.repeat(Tw, counts, axis=0)) > 1e-3       # is biased: the floor


def test_the_weight_fluctuation_splits_between_the_parameters_and_the_bias():
    # (iii): the fixed prior's weight sensitivity, whitened and weighted, is H_L(H_L + σ₀I)⁻¹ times the gauge's S^e = ∂_wϑ^e; the adaptive
    # prior carries all of S^e, so its parameter-level (coupled-distance) floor is the larger one
    loop, exact, prior = _setup(0.0)
    out = FB.bias_floor_law(loop, Q, exact, SIGMA0)
    S = PD.weight_sensitivity_prior(loop.nearest, Q, exact)
    V2 = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
    Sy = np.column_stack([(np.sqrt(Q)[:, None] * np.einsum('kcj,j->kc', S, V2[:, c]) / SQG[None, :]).ravel() for c in range(2)])
    Sey = -out['dx0']
    split = out['HL'] @ np.linalg.solve(out['HL'] + SIGMA0 * np.eye(9), Sey)
    assert np.linalg.norm(Sy - split) < 1e-5 * np.linalg.norm(Sy)
    assert np.linalg.norm(Sy - Sey) > 0.1 * np.linalg.norm(Sy)                 # the bias part is not negligible
    D_fixed, D_gauge = np.trace(Sy @ out['C2'] @ Sy.T), np.trace(Sey @ out['C2'] @ Sey.T)
    assert D_gauge > 1.5 * D_fixed
