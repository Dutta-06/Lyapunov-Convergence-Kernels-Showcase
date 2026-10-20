"""Proposition D‴.30 (PARADIGM.md §5d): measurement noise at perfect tracking. The averaged collective block's stationary covariance under
the diffusion (N₀/k_e)Q gives the function-level floor (N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) — at η = ∞ half the trace of Q whatever
H, increasing in η; the exactly linearised loop with white measurement noise has that floor plus the fast part, at η = 0, σ₀ and ∞; the
noise is coherent (a replicated start stays replicated, for every m, at every order) but the replicated state attracts only below a noise
level (the separations' second-order drift); the OU noise of the simulations has its stationary variance and correlation, and the law sees
the measurement; a constant force — the errors-in-variables term — is answered by (H + σ₀I)⁻¹ under the fixed prior and by H⁻¹ under the
adaptive prior for every η > 0, and the quasi-static force is the full second-order mean's noise part."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN

PER = 2 * np.pi
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
G = np.array([1.0, 0.25, 0.25])
SIGMA0 = 0.15 * 0.047162567584665846
F_VEC = lambda X: sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
FPRIME = lambda x: sum(-A * (x - C) / S ** 2 * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
XR = lambda t: np.array([2.0 * np.sin(t)])
XRD = lambda t: np.array([2.0 * np.cos(t)])
Q3 = np.ones(3) / 3


def _setup():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
    exact = net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))
    return net, exact


def test_the_closed_form_floor_is_the_lyapunov_solution_and_the_adaptive_floor_is_half_the_trace():
    rng = np.random.default_rng(0)
    U = np.linalg.qr(rng.standard_normal((9, 9)))[0]
    H = U @ np.diag(np.geomspace(1e-4, 0.2, 9)) @ U.T
    R = rng.standard_normal((9, 9)); Q = R @ R.T
    floors, pvars = [], []
    for eta in (0.0, SIGMA0, 10 * SIGMA0, np.inf):
        Sig = MN.averaged_covariance(H, Q, SIGMA0, eta, 2.0, 25.0)
        fl, sh, pv = MN.floor_closed_form(H, Q, SIGMA0, eta, 2.0, 25.0)
        assert abs(np.trace(H @ Sig) / fl - 1) < 1e-10                          # the closed form is the Lyapunov solution's
        floors.append(fl); pvars.append(pv)
    assert abs(floors[-1] / (2.0 / 50.0 * np.trace(Q)) - 1) < 1e-12              # η = ∞: (N₀/2k_e)tr Q
    assert np.all(np.diff(floors) > 0) and np.all(np.diff(np.array(pvars), axis=0) > 0)   # increasing in η, mode by mode


def test_the_linearised_loop_with_white_noise_has_the_predicted_floors():
    net, exact = _setup()
    ke = 100.0
    for eta in (0.0, SIGMA0, np.inf):
        cv = AF.CollapsedVariational(net, XR, PER, G, SIGMA0, np.inf if eta == 0 else eta, exact, Q3)
        gen, D = (MN.fixed_generator(cv, ke), 1 + cv.d) if eta == 0 else ((lambda t: cv.generator(t, ke)), cv.dim)
        b0 = MN.collapsed_noise_input(cv, FPRIME, ke)
        b = lambda t: b0(t)[:D]
        P0 = MN.periodic_covariance(gen, b, PER, D, 1.0)
        ts = np.linspace(0, PER, 33)
        tt, Ps = MN.transient_covariance(gen, b, D, 1.0, PER, ts, P0=P0, rtol=1e-10)
        H, Q = MN.noise_matrices(cv, FPRIME)
        trH2, J4 = MN.fast_parts(cv)
        orbit = np.trapezoid([np.trace(H @ P[1:10, 1:10]) for P in Ps], tt) / PER
        point = np.trapezoid([cv.J(t)[0] @ P[1:10, 1:10] @ cv.J(t)[0] for t, P in zip(tt, Ps)], tt) / PER
        po, pp = MN.floors(H, Q, trH2, J4, SIGMA0, eta, 1.0, ke)
        assert abs(orbit / po - 1) < 5e-3 and abs(point / pp - 1) < 5e-3          # the floors to 0.5% (0.1–0.2% here)
        assert np.max(np.abs(Ps[-1] - P0)) < 1e-6 * np.abs(P0).max()                 # the covariance is periodic


def test_measurement_noise_is_coherent_a_replicated_start_stays_replicated():
    net, exact = _setup()
    runs = []
    for labels, weights, m in ((np.repeat(np.arange(3), 2), None, 6), (np.arange(3), Q3, 3)):
        cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, 25.0, G, SIGMA0, labels, 3, np.inf, weights=weights)
        rec = lambda t, z, m=m, cl=cl: (z[1:1 + 3 * m].reshape(m, 3).copy())
        runs.append(MN.simulate_noisy(cl, cl.lift(exact), m, 10.0, 1e-3, 100.0, 1e-4, np.random.default_rng(7), 2000, rec))
    for Th6, Th3 in zip(*runs):
        assert np.abs(Th6[0::2] - Th6[1::2]).max() == 0.0                         # the classes stay exactly replicated …
        assert np.abs(Th6[0::2] - Th3).max() < 1e-12                               # … and follow the K₀-atom loop driven by the same noise
    assert np.abs(runs[1][-1] - exact).max() > 1e-6                                # the noise does move the parameters


def test_invariance_is_not_attraction_the_separations_grow_above_a_noise_level():
    """A separation of the third class (m = 6, along the top eigenvector of its second-order drift S₃) decays at σ₀/k_e − λ_max(S₃)Eν²:
    it decays at the noise std 0.01 and grows at 0.05, both rates within 10⁻⁴ of the prediction (100 time units, one noise path each)."""
    net, exact = _setup()
    cv = AF.CollapsedVariational(net, XR, PER, G, SIGMA0, np.inf, exact, Q3)
    fpp = lambda x: sum(A * ((x - C) ** 2 / S ** 4 - 1 / S ** 2) * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
    S3 = MN.separation_drift(cv, FPRIME, fpp, 25.0, 100.0, 1.0, N=256)[2]
    lam, U = np.linalg.eigh(S3)
    assert np.allclose(S3, S3.T) and lam[-1] > 0.1                                   # a positive drift on the weakest class
    labels = np.repeat(np.arange(3), 2)
    cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, 25.0, G, SIGMA0, labels, 3, np.inf)
    rates = {}
    for std in (0.01, 0.05):
        z0 = cl.lift(exact)
        z0[13:16] += 1e-6 * np.sqrt(G) * U[:, -1]; z0[16:19] -= 1e-6 * np.sqrt(G) * U[:, -1]
        rec = lambda t, z: np.sqrt(np.sum(((z[13:16] - z[16:19]) / np.sqrt(G)) ** 2))
        out = np.array(MN.simulate_noisy(cl, z0, 6, 100.0, 2e-3, 100.0, std ** 2, np.random.default_rng(11), 2500, rec))
        rates[std] = np.polyfit(np.arange(len(out))[2:] * 5.0, np.log(out[2:]), 1)[0]
        assert abs(rates[std] - (-SIGMA0 / 25.0 + lam[-1] * std ** 2)) < 1e-4
    assert rates[0.01] < 0 < rates[0.05]


def test_the_ou_noise_has_its_stationary_variance_and_the_law_sees_the_measurement():
    """`simulate_noisy`'s ν is the OU chain of rate B and variance var held within a step (its sample variance and lag-one correlation
    e^{−Bdt}); and at the lift the adaptive law's increment under ν is G∇_ϑ(νΦ̃(x_m + ν; ϑ)) — the measured error at the measured point."""
    class Stub:
        class net:
            n = 1

        def __init__(self):
            self.nus = []

        def rhs_noisy(self, t, z, m, nu):
            self.nus.append(float(nu[0]))
            return np.zeros_like(z)
    st = Stub()
    MN.simulate_noisy(st, np.zeros(1), 1, 100.0, 1e-3, 100.0, 4e-4, np.random.default_rng(5), 10 ** 9, lambda t, z: None)
    nu = np.array(st.nus[::4])
    assert abs(nu.var() / 4e-4 - 1) < 0.08 and abs(np.corrcoef(nu[:-1], nu[1:])[0, 1] - np.exp(-0.1)) < 0.01
    net, exact = _setup()
    cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, 25.0, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3)
    z = cl.lift(exact); t = 1.1; z[0] = XR(t)[0]
    for nu in (0.0, 0.03):
        law = cl.rhs_noisy(t, z, 3, np.array([nu]))[1:10].reshape(3, 3)
        assert np.abs(law - G[None, :] * net.grad_g(XR(t) + nu, np.array([nu]), exact)).max() < 1e-12
    assert np.abs(law).max() > 1e-4


def test_a_constant_force_is_answered_by_the_inverse_curvature_under_the_adaptive_prior():
    rng = np.random.default_rng(1)
    U = np.linalg.qr(rng.standard_normal((6, 6)))[0]
    h = np.geomspace(1e-4, 0.2, 6); H = U @ np.diag(h) @ U.T
    F = rng.standard_normal(6)
    s_end = 1e8
    y_fixed = MN.mean_response(H, F, SIGMA0, 0.0, s_end)
    assert np.allclose(y_fixed, np.linalg.solve(H + SIGMA0 * np.eye(6), F), rtol=1e-8, atol=1e-12)
    for eta in (SIGMA0, np.inf):
        assert np.allclose(MN.mean_response(H, F, SIGMA0, eta, s_end), np.linalg.solve(H, F), rtol=1e-6, atol=1e-10)
    net, exact = _setup()
    cv = AF.CollapsedVariational(net, XR, PER, G, SIGMA0, np.inf, exact, Q3)
    fpp = lambda x: sum(A * ((x - C) ** 2 / S ** 4 - 1 / S ** 2) * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
    # the quasi-static force (the closed forms of E[ε_m²] and E[e]) against the full second-order mean's noise part — the periodic covariance
    # contracted with the true field's second derivatives — and the parameters' fluctuation's part, negligible here for the measured regressor
    Fq = 25.0 * MN.eiv_force(cv, FPRIME, fpp, 25.0, 100.0, 1.0)
    cl = MN.NoisyRealLoop(net, F_VEC, XR, XRD, PER, 25.0, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3)
    out = MN.second_order_mean(cl, cv, FPRIME, 100.0, 1.0, N=128)
    assert np.linalg.norm(out['F_noise'] - Fq) < 2e-3 * np.linalg.norm(Fq)
    H = MN.noise_matrices(cv, FPRIME)[0]
    b = MN.mean_biases(H, out)
    assert abs(b['full'] / b['noise'] - 1) < 1e-4 and abs(b['noise'] / MN.equivalent_perturbation(H, Fq) - 1) < 5e-3
