"""
Proposition D‴.41 (PARADIGM.md §5d): measurement noise on the gauge — the selection is set by the noise's per-update variance
==========================================================================================================================

The measured error r + ξ drives both laws of (𝔉), so the noise enters along the descent's own vector fields: the drift is
⟨rU⟩_ρ = (1/N)Σ_x r(x)U_x and the noise (1/N)Σ_x ξ(x)U_x, with U_x = mM∇f̂(x), M = diag(G, γ/β) on (Θ, log β) projected on Σβ = m:

    U_x = (β_iG∇_ϑΦ̃(x; ϑ_i),  γ(Φ̃(x; ϑ_i) − f̂(x)))_i.

On the manifold M of exact fits (r ≡ 0) the fields are orthogonal to M in the flow's metric, and the limit map's invariance makes the
Wong–Zakai term cancel the Itô drift: with Stratonovich (or any smooth, coloured) noise the gauge does not move at order v, with noise
drawn afresh at each update it moves by −(v/2)∂Φ⟨(U·∇)U⟩_ρ, which on the masses in the slow-replicator limit is the replicator flow of

    ū_k = −(v/2)β_kG_α⟨Φ_k²⟩_ρ(1 − (a_k/a_max)⁴) = −∂_{β_k}Ψ_ℓ,   Ψ_ℓ = (v/4)Σ_kG_α⟨Φ_k²⟩_ρ(β_ka′_k)²,

the amplitudes' part of (v/4)tr(G∇²J) (the rest of the trace is constant on M); its minimiser on the simplex is interior.

This module: the fields, the Wong–Zakai term, the closed-form selection, Ψ_ℓ and its minimiser, the per-step factor of AR(1) noise.
"""
import numpy as np


def _sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def fields(loop, Th, beta):
    """U_x for every quadrature node in the coordinates (Θ, log β): (UT (N, m, p), Ul (N, m)); `loop` a `ReplicatorCohesionLoop`."""
    base = loop.base
    net = base.net
    a, c, s, diff, phi = base.features(Th)
    da = 1.0 - (a / net.a_max) ** 2
    dc = 1.0 - (c / net.c_max) ** 2
    sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
    gT = np.empty((base.N,) + Th.shape)
    gT[:, :, 0] = phi * da
    gT[:, :, 1] = a * dc / s ** 2 * phi * diff
    gT[:, :, 2] = a * ds / s ** 3 * phi * diff ** 2
    UT = beta[None, :, None] * loop.G[None, None, :] * gT
    Phit = phi * a[None, :]
    fh = Phit @ beta / Th.shape[0]
    Ul = loop.gamma * (Phit - fh[:, None])
    return UT, Ul


def selection(loop, Th, beta, v=1.0):
    """The closed form on M in the slow-replicator limit: ū_k = −(v/2)β_kG_α⟨Φ_k²⟩_ρ(1 − q_k²), q_k = (a_k/a_max)², (m,)."""
    a, c, s, diff, phi = loop.base.features(Th)
    q = (a / loop.base.net.a_max) ** 2
    return -0.5 * v * beta * loop.G[0] * np.mean(phi ** 2, axis=0) * (1 - q ** 2)


def psi(W, P, W_min, v=1.0, m=3):
    """Ψ_ℓ at the masses W (Σ W = 1) on M: (v/4)Σ_kP_k(β_ka′_k)², β_ka′_k = m(W_k − W_min,k²/W_k), P_k = G_α⟨Φ_k²⟩_ρ."""
    W = np.asarray(W, float)
    return 0.25 * v * float(np.sum(np.asarray(P, float) * (m * (W - np.asarray(W_min, float) ** 2 / W)) ** 2))


def w_label(P, W_min):
    """The minimiser W_ℓ of Ψ_ℓ on the simplex: P_k(W_k − W_min,k²/W_k)(1 + W_min,k²/W_k²) equal (increasing in W_k from 0 at W_min,k),
    so W_ℓ,k > W_min,k for every k when Σ_kW_min,k < 1."""
    from scipy.optimize import brentq
    P = np.asarray(P, float); Wm = np.asarray(W_min, float)
    g = lambda W, k: P[k] * (W - Wm[k] ** 2 / W) * (1 + Wm[k] ** 2 / W ** 2)
    W_of = lambda mu, k: brentq(lambda W: g(W, k) - mu, Wm[k], 1.0 + 10 * Wm[k])
    hi = min(g(1.0, k) for k in range(len(P)))
    mu = brentq(lambda mu: sum(W_of(mu, k) for k in range(len(P))) - 1.0, 1e-12 * hi, hi)
    return np.array([W_of(mu, k) for k in range(len(P))])


def ar1_factor(rho):
    """The selection of AR(1) noise (successive draws correlated by ρ) over that of noise drawn afresh, at the same intensity: (1 − ρ)/(1 + ρ)
    (the selection is set by the per-step variance, which the correlation lowers at a fixed intensity)."""
    return (1.0 - rho) / (1.0 + rho)


def flat_fields(loop, x, beta_of, m=3):
    """U_x in the flat coordinates (Θ (m·p), log(β_j/β_m), j < m) of `exp_gauge_drift.instance`: (N, m·p + m − 1)."""
    Th = x[:m * 3].reshape(m, 3)
    UT, Ul = fields(loop, Th, beta_of(x[m * 3:]))
    return np.concatenate([UT.reshape(UT.shape[0], -1), Ul[:, :m - 1] - Ul[:, m - 1:]], axis=1)


def wong_zakai(loop, x, beta_of, h=1e-4, m=3, weights=None):
    """⟨w(U·∇)U⟩_ρ = (1/N)Σ_x w_x ∂_{U_x}U_x at the flat point x (central differences), (m·p + m − 1,); w ≡ 1 unless `weights` (the noise's
    variance profile over the nodes, heteroscedastic noise) is given."""
    U0 = flat_fields(loop, x, beta_of, m)
    wt = np.ones(U0.shape[0]) if weights is None else np.asarray(weights, float)
    out = np.zeros(U0.shape[1])
    for i in range(U0.shape[0]):
        out += wt[i] * (flat_fields(loop, x + h * U0[i], beta_of, m)[i] - flat_fields(loop, x - h * U0[i], beta_of, m)[i]) / (2 * h)
    return out / U0.shape[0]
