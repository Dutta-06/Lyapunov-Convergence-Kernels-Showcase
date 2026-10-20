"""
Proposition D‴.37 (PARADIGM.md §5d): the dither's drift along a manifold of equilibria. Checked here: the second-order identities of the
limit map against flows whose limit map is known in closed form (a curved coordinate system on a linear flow; a curved manifold), and
against the limit map computed by integration for a nonlinear gradient flow in a non-orthogonal metric; on Proposition D‴.36's instance,
the diffusion against D‴.36's projection formula, and the drift at small γ against the replicator flow of the log-determinant.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
from src.core import gauge_drift as GDR                                    # noqa: E402
import exp_gauge_dither as GD                                             # noqa: E402
import exp_gauge_drift as X                                               # noqa: E402


def test_curved_coordinates_on_a_linear_flow():
    # x′ = 0, z′ = −λz written in u = x + κz², v = z: Φ(u, v) = (u − κv², 0), so ½∂²Φ_u[Σ] = −κΣ_vv
    lam, kap = 0.7, 0.4
    F = lambda y: np.array([-2 * kap * lam * y[1] ** 2, -lam * y[1]])
    chart = lambda w: np.array([w[0], 0.0])
    coords = lambda y: np.array([y[0]])
    S = np.diag([3e-4, 2e-4])
    d, D, ev = GDR.limit_drift(F, chart, coords, np.array([0.2]), S)
    assert abs(d[0] / (-kap * S[1, 1]) - 1) < 1e-6 and abs(D[0, 0] / S[0, 0] - 1) < 1e-8 and abs(ev[0] + lam) < 1e-8


def test_a_curved_manifold():
    # x′ = 0, z′ = −λz in u = x, v = z + qx²: M = {v = qu²}, Φ = (u, qu²); the drift of v is qΣ_uu (the chart's curvature), of u zero
    lam, q, u0 = 0.5, 0.8, 0.3
    F = lambda y: np.array([0.0, -lam * (y[1] - q * y[0] ** 2)])
    chart = lambda w: np.array([w[0], q * w[0] ** 2])
    S = np.diag([2e-4, 1e-4])
    dv, _, _ = GDR.limit_drift(F, chart, lambda y: np.array([y[1]]), np.array([u0]), S)
    du, _, _ = GDR.limit_drift(F, chart, lambda y: np.array([y[0]]), np.array([u0]), S)
    assert abs(dv[0] / (q * S[0, 0]) - 1) < 1e-6 and abs(du[0]) < 1e-10


def test_against_the_integrated_limit_map():
    # a gradient flow of J = ½(1 + u/2)(v − qu²)² in the metric g: the limit map by integration, its second differences along the noise
    q = 0.6
    g = np.array([[1.0, 0.3], [0.3, 2.0]]); gi = np.linalg.inv(g)

    def F(y):
        u, v = y
        r = v - q * u * u
        grad = np.array([0.25 * r * r - (1 + 0.5 * u) * r * 2 * q * u, (1 + 0.5 * u) * r])
        return -gi @ grad
    chart = lambda w: np.array([w[0], q * w[0] ** 2])
    coords = lambda y: np.array([y[0]])
    u0 = 0.4
    S = np.diag([2e-4, 1e-4])
    d, _, _ = GDR.limit_drift(F, chart, coords, np.array([u0]), S)
    x0 = chart(np.array([u0])); hh = 2e-2
    bf = 0.0
    for i in range(2):
        e = np.eye(2)[i] * hh
        Pp = GDR.limit_map(F, x0 + e, 400.0)[0][0]; Pm = GDR.limit_map(F, x0 - e, 400.0)[0][0]
        bf += 0.5 * S[i, i] * (Pp - 2 * u0 + Pm) / hh ** 2
    assert abs(d[0] / bf - 1) < 2e-3


def test_the_instances_diffusion_is_the_projection_formula():
    w = np.array([0.30, 0.35])
    d, D, ev = X.drift_at(w)
    D3 = GD.gauge_diffusion(np.r_[w, 1 - w.sum()], X.GAMMA, X.TAU)
    assert np.abs(D - D3[:2, :2]).max() < 1e-4 * np.abs(D3).max() and np.all(ev < 0)


def test_the_slow_replicator_limit_is_the_log_determinants_replicator_flow():
    for w in (np.array([1 / 3, 1 / 3]), np.array([0.40, 0.25])):
        d = X.drift_at(w, gamma=1e-5)[0] / (1e-5 * X.TAU)
        law = X.log_det_law(np.r_[w, 1 - w.sum()])[:2]
        assert np.abs(d - law).max() < 1e-2 * np.abs(law).max()
    d2 = X.drift_at(np.array([1 / 3, 1 / 3]), gamma=1e-2)[0] / (1e-2 * X.TAU)    # at γ = 10⁻² the finite-γ correction is O(1)
    assert np.abs(d2 - X.log_det_law(np.ones(3) / 3)[:2]).max() > 0.1


def test_the_log_determinant_does_not_see_the_split_within_a_class():
    # two coincident particles in the first class: the fit's block is rank p along the class's β-weighted mean, its Schur complement the
    # cohesion's on the separation, so det ∇²_ΘJ depends on the class's total mass only
    import exp_cohesion_leak as E
    from src.core import mass_transfer as MT
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    T3 = E.exact_rep(net, np.ones(3) / 3)
    x0 = np.vstack([T3[0], T3[0], T3[1], T3[2]]).ravel()
    out = []
    for t in (0.5, 0.2):
        beta = np.array([t, 1 - t, 1.0, 1.0]); beta *= 4 / beta.sum()
        J = lambda x: loop.J(x.reshape(4, 3), beta)
        h = 1e-4; I = np.eye(12) * h; H = np.zeros((12, 12))
        for i in range(12):
            for j in range(i, 12):
                H[i, j] = H[j, i] = (J(x0 + I[i] + I[j]) - J(x0 + I[i] - I[j]) - J(x0 - I[i] + I[j]) + J(x0 - I[i] - I[j])) / (4 * h * h)
        out.append(np.linalg.slogdet(H))
    assert out[0][0] > 0 and out[1][0] > 0 and abs(out[0][1] - out[1][1]) < 1e-5


def test_the_ito_term_of_a_nonlinear_observable():
    # round 27 (D3z-8): on the first toy, observe w = u²: dw = 2u·(−κΣ_vv) + ½·2·Σ_uu (the Itô term of the observable)
    lam, kap, u0 = 0.7, 0.4, 0.2
    F = lambda y: np.array([-2 * kap * lam * y[1] ** 2, -lam * y[1]])
    chart = lambda w: np.array([np.sqrt(w[0]), 0.0])
    coords = lambda y: np.array([y[0] ** 2])
    S = np.diag([3e-4, 2e-4])
    d, D, _ = GDR.limit_drift(F, chart, coords, np.array([u0 ** 2]), S)
    assert abs(d[0] / (2 * u0 * (-kap * S[1, 1]) + S[0, 0]) - 1) < 1e-5 and abs(D[0, 0] / (4 * u0 ** 2 * S[0, 0]) - 1) < 1e-6


def test_the_instance_at_gamma_1e2_against_the_integrated_limit_map():
    # round 27 (D3z-8): Part A's brute-force value at W = 1/3 (LSODA to 1.2·10⁵, second differences at ±0.02), drift of (W_1, W_2)
    d = X.drift_at(np.array([1 / 3, 1 / 3]))[0]
    assert np.abs(d - np.array([-2.05497e-7, 1.57975e-7])).max() < 2e-3 * 2.05e-7


def test_the_log_determinant_closed_form():
    # round 27 (D3z-6): with one particle per class ℓ = 2Σ log(W_k − W_min,k²/W_k) + const
    Wm = X.W_MIN
    vals = [X.log_det(3 * W) - 2 * np.sum(np.log(W - Wm ** 2 / W)) for W in (np.ones(3) / 3, np.array([0.27, 0.36, 0.37]))]
    assert abs(vals[0] - vals[1]) < 1e-6
