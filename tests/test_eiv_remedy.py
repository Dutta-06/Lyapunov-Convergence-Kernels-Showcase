"""Proposition D‴.31 (PARADIGM.md §5d): the adaptive prior is a washout leak — its leak target lags the class barycentre by the barycentre's
exponentially forgotten increments, so the leak acts on the separations in full and on the barycentre through the high-pass iω/(iω + η) —
and along every affine time-invariant prior the stationary response R to a persistent force and the bias B toward a misplaced anchor satisfy
HR + B = I; the function-level response to a force is at most the norm of the smallest target perturbation exerting it, under every leak
S ≥ 0; and the desired regressor (the adaptive law's regressor at the reference) leaves the lift, the linearisation and the first-order
noise response unchanged while turning the second-order force into the least-squares force of g = k_eE[e]."""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SIGMA0 = 0.15 * 0.047162567584665846
Q3 = np.ones(3) / 3
F_VEC = lambda X: sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
FP = lambda x: sum(-A * (x - C) / S ** 2 * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
FPP = lambda x: sum(A * ((x - C) ** 2 / S ** 4 - 1 / S ** 2) * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
XR = lambda t: np.array([2.0 * np.sin(t)])
XRD = lambda t: np.array([2.0 * np.cos(t)])
KE = 25.0


def _setup():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    atoms = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    cv = AF.CollapsedVariational(net, XR, PER, G, SIGMA0, np.inf, atoms, Q3)
    return net, atoms, cv


def test_the_prior_lags_the_barycentre_by_its_forgotten_increments():
    """(i) on the real loop (finite η, a perturbed start, classes (0, 0, 1, 2)): T̂ − z = ∫₀ᵗe^{−(η/k_e)(t−τ)}dT̂(τ) + e^{−(η/k_e)t}(T̂(0) − z(0)),
    the convolution integrated along the solution from the barycentre's velocity (the class mean of Θ̇)."""
    net, atoms, cv = _setup()
    eta = 5 * SIGMA0
    labels = np.array([0, 0, 1, 2]); m, p = 4, 3
    cl = AF.AdaptiveRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, labels, 3, eta, weights=np.array([1, 1, 2, 2]) / 6)
    z0 = cl.lift(atoms)
    rng = np.random.default_rng(3)
    z0[1:1 + m * p] += 0.05 * rng.standard_normal(m * p)
    z0[1 + m * p:] += 0.05 * rng.standard_normal(3 * p)
    c = eta / KE
    T = 3 * PER

    def aug(t, y):
        z = y[:cl.state_dim(m)]
        dz = cl.rhs(t, z, m)
        dT = cl.class_means(dz[1:1 + m * p].reshape(m, p))
        return np.concatenate([dz, (dT - c * y[cl.state_dim(m):].reshape(3, p)).ravel()])
    y0 = np.concatenate([z0, np.zeros(3 * p)])
    sol = solve_ivp(aug, (0, T), y0, method='DOP853', rtol=1e-11, atol=1e-13)
    yT = sol.y[:, -1]
    Th0 = z0[1:1 + m * p].reshape(m, p); Z0 = z0[1 + m * p:].reshape(3, p)
    ThT = yT[1:1 + m * p].reshape(m, p); ZT = yT[1 + m * p:1 + m * p + 3 * p].reshape(3, p); W = yT[cl.state_dim(m):].reshape(3, p)
    lhs = cl.class_means(ThT) - ZT
    rhs = W + np.exp(-c * T) * (cl.class_means(Th0) - Z0)
    assert np.abs(lhs - rhs).max() < 1e-9 and np.abs(W).max() > 1e-4


def test_the_pair_block_is_the_washed_out_leak_and_the_stationary_identity_holds_along_the_prior_family():
    """(i) linearised: the transfer function of the averaged collective block from a force on y to y is, mode by mode,
    1/(iω + h_j + σ₀(iω + κ)/(iω + η + κ)); (ii) at zero frequency HR + B = I with R = (H + S₀)⁻¹, S₀ = σ₀κ/(η + κ), for the adaptive prior
    (κ = 0: S₀ = 0, B = 0, R = H⁻¹), the anchored priors and the fixed prior."""
    net, atoms, cv = _setup()
    H, Q = MN.noise_matrices(cv, FP)
    h, V = np.linalg.eigh(H)
    d = len(h); eta = SIGMA0
    for kap in (0.0, 0.1 * eta, eta, 10 * eta):
        A = MN.collective_block(H, SIGMA0, eta, kap)
        for om in (0.0, 0.3 * eta, 3 * eta, 30 * eta):
            Rw = np.linalg.solve(1j * om * np.eye(A.shape[0]) - A, np.vstack([np.eye(d), np.zeros((d, d))]))[:d]
            pred = V @ np.diag(MN.washout_transfer(h, SIGMA0, eta, om, kap)) @ V.T
            assert np.abs(Rw - pred).max() < 1e-9 * np.abs(pred).max()
        R, Bm = MN.stationary_operators(H, SIGMA0, eta, kap)
        S0 = MN.zero_frequency_leak(SIGMA0, eta, kap)
        assert np.abs(H @ R + Bm - np.eye(d)).max() < 1e-12
        assert np.abs(R - np.linalg.inv(H + S0 * np.eye(d))).max() < 1e-10 * np.abs(R).max()
    R, Bm = MN.stationary_operators(H, SIGMA0, 0.0)
    assert np.abs(H @ R + Bm - np.eye(d)).max() < 1e-12 and np.abs(Bm - SIGMA0 * np.linalg.inv(H + SIGMA0 * np.eye(d))).max() < 1e-12
    R, Bm = MN.stationary_operators(H, SIGMA0, eta, 0.0)
    assert np.abs(Bm).max() < 1e-12 and np.abs(R - np.linalg.inv(H)).max() < 1e-8 * np.abs(R).max()


def test_the_response_to_a_force_is_at_most_its_smallest_equivalent_target_perturbation_under_every_leak():
    """(iii): Fᵀ(H + S)⁻¹H(H + S)⁻¹F ≤ FᵀH⁻¹F for S ≥ 0, with equality at S = 0; FᵀH⁻¹F is the squared norm of g = ĴH⁻¹F, the smallest g with
    ⟨Ĵᵀg⟩ = F; a least-squares force ⟨Ĵᵀg⟩ has a response at most ‖g‖²_ρ, and at S = 0 the response Ĵy is the projection of g."""
    net, atoms, cv = _setup()
    N = 512
    t = np.arange(N) * PER / N
    J = np.array([cv.J(s)[0] for s in t])
    H = J.T @ J / N
    rng = np.random.default_rng(0)
    d = H.shape[0]
    for trial in range(20):
        F = rng.standard_normal(d)
        M = rng.standard_normal((d, d))
        S = [0.0, SIGMA0, 10 * SIGMA0, 0.01 * M @ M.T][trial % 4]
        assert MN.response_size(H, F, S) <= MN.equivalent_perturbation(H, F) * (1 + 1e-10)
        g_min = J @ np.linalg.solve(H, F)
        assert abs(np.mean(g_min ** 2) / MN.equivalent_perturbation(H, F) - 1) < 1e-9 and np.abs(J.T @ g_min / N - F).max() < 1e-9 * np.abs(F).max()
        g = rng.standard_normal(N)
        Fg = J.T @ g / N
        assert MN.response_size(H, Fg, S) <= np.mean(g ** 2)
        y = np.linalg.solve(H, Fg)
        Pg = J @ np.linalg.solve(J.T @ J, J.T @ g)
        assert np.abs(J @ y - Pg).max() < 1e-8 * np.abs(Pg).max()
    assert abs(MN.response_size(H, F, 0.0) / MN.equivalent_perturbation(H, F) - 1) < 1e-10


def test_the_desired_regressor_keeps_the_lift_and_the_first_order_response():
    """(iv)(a)–(b): without noise the lift stays exact under DR and DR*; on a common noise path DR's run differs from MR's by the second-order
    terms only — the difference falls like the noise's variance while each deviation falls like its standard deviation."""
    net, atoms, cv = _setup()
    grid = 2 * np.sin(np.arange(64) * PER / 64)
    runs = {}
    for reg, comp in (('measured', 'measured'), ('reference', 'measured'), ('reference', 'reference')):
        cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, KE, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3, regressor=reg, compensation=comp)
        z0 = cl.lift(atoms)
        rec = lambda t, z: np.array([cl.f_hat(np.array([x]), z[1:10].reshape(3, 3))[0] for x in grid]) - F_VEC(grid)
        clean = MN.simulate_noisy(cl, z0, 3, 2.0, 1e-3, 100.0, 0.0, np.random.default_rng(0), 2000, rec)
        assert np.abs(clean[-1]).max() < 1e-12                                       # the lift is a solution
        for std in (1e-3, 1e-2):
            runs[(reg, comp, std)] = MN.simulate_noisy(cl, z0, 3, 5.0, 1e-3, 100.0, std ** 2, np.random.default_rng(7), 5000, rec)[-1]
    dev = {s: np.sqrt(np.mean(runs[('measured', 'measured', s)] ** 2)) for s in (1e-3, 1e-2)}
    dif = {s: np.sqrt(np.mean((runs[('measured', 'measured', s)] - runs[('reference', 'measured', s)]) ** 2)) for s in (1e-3, 1e-2)}
    assert 7 < dev[1e-2] / dev[1e-3] < 14                                              # first order: ∝ std
    assert 70 < dif[1e-2] / dif[1e-3] < 140 and dif[1e-3] < 0.05 * dev[1e-3]            # MR − DR: ∝ std², small at small noise
    dif_star = np.sqrt(np.mean((runs[('measured', 'measured', 1e-3)] - runs[('reference', 'reference', 1e-3)]) ** 2))
    assert dif_star < 0.1 * dev[1e-3]                                                 # DR*: the same response up to O(a/k_e)


def test_the_desired_regressor_turns_the_errors_in_variables_force_into_a_least_squares_force():
    """(iv)(c): DR's and DR*'s forces are ⟨Ĵᵀg⟩/k_e (real time) with g = k_eE[e], |g| ≤ ½|f″|Eν²(1 + 2|a|/k_e) along the orbit (DR) and
    ≤ ½|f″|Eν²(k_e/(k_e − a))³ (DR*), their adaptive-prior responses at most ‖g‖²_ρ; MR's force exceeds its least-squares part by
    F₁ = k_e⟨Ĵ′ᵀE[ε_m²]⟩, whose projected mean tracking error F₁ᵀH⁻¹F₁ is three orders of magnitude above DR's bias here."""
    net, atoms, cv = _setup()
    H, Q = MN.noise_matrices(cv, FP)
    N = 1024
    t = np.arange(N) * PER / N
    J = np.array([cv.J(s)[0] for s in t])
    x = 2 * np.sin(t); a = FP(x); b2 = FPP(x)
    for ke, B in ((25.0, 100.0), (100.0, 1000.0)):
        Fr, g = MN.eiv_force(cv, FP, FPP, ke, B, 1.0, variant='reference', return_g=True)
        assert np.abs(Fr - J.T @ (g / ke) / N).max() < 1e-12
        assert np.all(np.abs(g) <= 0.5 * np.abs(b2) * (1 + 2 * np.abs(a) / ke) + 1e-15)
        Fb, gb = MN.eiv_force(cv, FP, FPP, ke, B, 1.0, variant='reference_both', return_g=True)
        assert np.all(np.abs(gb) <= 0.5 * np.abs(b2) * (ke / (ke - a)) ** 3 + 1e-15)
        for FF, gg in ((Fr, g), (Fb, gb)):
            assert MN.equivalent_perturbation(H, ke * FF) <= np.mean(gg ** 2) * (1 + 1e-9)
        Fm = MN.eiv_force(cv, FP, FPP, ke, B, 1.0, variant='measured')
        F1 = ke * (Fm - Fr)
        assert MN.equivalent_perturbation(H, F1) > 1e3 * MN.equivalent_perturbation(H, ke * Fr)
