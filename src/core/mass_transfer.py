"""
Proposition D‴.34 (PARADIGM.md §5d): mass transfer — the replicator weight law with the cohesion leak
=====================================================================================================

The particles carry weights β_i > 0, the network f̂ = (1/m)Σ_iβ_iΦ̃(·; ϑ_i) and its measure μ_{Θ,β} = (1/m)Σ_iβ_iδ_{ϑ_i}, and the weights
adapt by the replicator law

    dβ_i/ds = γβ_i(u_i − λ),   u_i := ⟨r Φ̃(·; ϑ_i)⟩_ρ,   λ := Σ_jβ_ju_j/Σ_jβ_j,                                    (𝔉)

(real time: β̇_i = γβ_ieᵀ(Φ̃(x; ϑ_i) − f̂(x)) when Σβ = m), which keeps Σ_iβ_i = m — the network's output bounded by the saturation as for
equal weights —, beside the transport law with the cohesion leak (Proposition D‴.32),

    dϑ_i/ds = β_iG⟨∇_ϑ(r Φ̃(·; ϑ_i))⟩_ρ − (σ₀/m)Σ_jψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²)(ϑ_i − ϑ_j).

The pair is the gradient flow of J(Θ, β) = mL̄(μ_{Θ,β}) + (σ₀m/2)∫∫W dμ_Θdμ_Θ (the cohesion energy of the unweighted cloud) in the metric
Σ_i‖δϑ_i‖²_{G⁻¹} + Σ_iδβ_i²/(γβ_i) restricted to Σ_iδβ_i = 0 (Euclidean on the positions, Fisher–Rao on the weights): ∂_{β_i}(mL̄) = −u_i and
∂_{ϑ_i}(mL̄) = −β_i⟨r∇_ϑΦ̃(·; ϑ_i)⟩_ρ. With β ≡ 1 and γ = 0 it is the loop of Proposition D‴.32. In real time with the positions fixed at
the atoms of an exact representation (weights β*, Σβ* = m) the replicator law has the Lyapunov function ½|e|² + (1/(γm))Σ_iβ*_i log(β*_i/β_i),
whose derivative is −k_e|e|².

This module: the weighted averaged loop on top of an `AveragedSlowLoop` (its quadrature, features and gradients), the energy, the run (the
midpoint rule on (Θ, log β)), the class masses, and the real-loop law with the KL Lyapunov function.
"""
import numpy as np

from . import cohesion_leak as CO


class ReplicatorCohesionLoop:
    """The cohesion leak (range R, kind, strength σ₀) with the replicator weight law of rate γ on top of an `AveragedSlowLoop` `base`;
    `balanced` keeps Σβ = m (λ the β-mean of u), otherwise λ = 0."""

    def __init__(self, base, R, kind='gauss', sigma0=None, gamma=1.0, balanced=True, eps=0.0):
        self.base = base
        self.G = base.G
        self.R, self.kind = float(R), kind
        self.sigma0 = base.sigma0 if sigma0 is None else float(sigma0)
        self.gamma, self.balanced = float(gamma), bool(balanced)
        self.eps = float(eps)                     # Proposition D‴.36: the entropic leak on the masses, energy (ε/γ)Σβ_i log β_i

    def f_hat(self, Th, beta):
        a, c, s, diff, phi = self.base.features(Th)
        return phi @ (beta * a) / Th.shape[0]

    def residual(self, Th, beta):
        return self.base.F - self.f_hat(Th, beta)

    def fitness(self, r, Th):
        """u_i = ⟨r Φ̃(·; ϑ_i)⟩_ρ, (m,)."""
        a, c, s, diff, phi = self.base.features(Th)
        return a * (phi.T @ r) / self.base.N

    def drift(self, Th, beta):
        r = self.residual(Th, beta)
        dTh = beta[:, None] * self.G[None, :] * self.base.grad_weighted(r, Th) + CO.cohesion_force(Th, self.G, self.R, self.sigma0, self.kind)
        u = self.fitness(r, Th)
        if self.eps:                              # −∂_β of (ε/γ)Σβ log β, up to the multiplier
            u = u - (self.eps / self.gamma) * np.log(beta)
        lam = float(beta @ u / beta.sum()) if self.balanced else 0.0
        return dTh, self.gamma * beta * (u - lam)

    def J(self, Th, beta):
        m = Th.shape[0]
        r = self.residual(Th, beta)
        out = 0.5 * m * float(np.mean(r ** 2)) + 0.5 * self.sigma0 * m * CO.interaction_energy(Th, self.G, self.R, self.kind)
        if self.eps:
            out += (self.eps / self.gamma) * float(np.sum(beta * np.log(beta)))
        return out

    def fit_error(self, Th, beta):
        r = self.residual(Th, beta)
        return float(np.sqrt(np.mean(r ** 2) / np.mean(self.base.F ** 2)))

    def run(self, Th, beta, S_end, dt, rec_every=200, callback=None):
        """The midpoint rule in slow time on (Θ, log β) (the weights stay positive; with `balanced` they are renormalised to Σβ = m after
        each step, which the exact flow conserves); callback(s, Θ, β) at the records. Returns (Θ, β) at S_end."""
        Th, lb = np.array(Th, float), np.log(np.array(beta, float))
        m = Th.shape[0]
        n = int(round(S_end / dt))
        for k in range(n + 1):
            b = np.exp(lb)
            if callback is not None and k % rec_every == 0:
                callback(k * dt, Th, b)
            if k == n:
                break
            d1, g1 = self.drift(Th, b)
            Tm, lm = Th + 0.5 * dt * d1, lb + 0.5 * dt * g1 / b
            d2, g2 = self.drift(Tm, np.exp(lm))
            Th, lb = Th + dt * d2, lb + dt * g2 / np.exp(lm)
            if self.balanced:
                lb = lb - np.log(np.exp(lb).sum() / m)
        return Th, np.exp(lb)


    def run_dithered(self, Th, beta, S_end, dt, tau, rng, rec_every=4000, callback=None):
        """Proposition D‴.36: the positions dithered (Euler–Maruyama, the noise √(2τG) per particle), the weights by the midpoint-free Euler
        step of (𝔉) on log β, renormalised to Σβ = m; callback(s, Θ, β) at the records."""
        Th, lb = np.array(Th, float), np.log(np.array(beta, float))
        m = Th.shape[0]
        sq = np.sqrt(2.0 * tau * dt * self.G)[None, :]
        n = int(round(S_end / dt))
        for k in range(n + 1):
            b = np.exp(lb)
            if callback is not None and k % rec_every == 0:
                callback(k * dt, Th, b)
            if k == n:
                break
            dT, db = self.drift(Th, b)
            Th = Th + dt * dT + sq * rng.standard_normal(Th.shape)
            lb = lb + dt * db / b
            if self.balanced:
                lb = lb - np.log(np.exp(lb).sum() / m)
        return Th, np.exp(lb)


def class_masses(Th, beta, G, thr):
    """The clusters of Θ (single linkage at thr) with their masses Σ_{i∈k}β_i/m and their β-weighted atoms: (labels, masses (K,), atoms (K, p))."""
    lab = CO.clusters(Th, G, thr)
    K = int(lab.max() + 1)
    m = Th.shape[0]
    W = np.array([beta[lab == k].sum() / m for k in range(K)])
    A = np.array([(beta[lab == k, None] * Th[lab == k]).sum(axis=0) / beta[lab == k].sum() for k in range(K)])
    return lab, W, A


def replicator_real_rhs(t, z, net, f_vec, x_ref, x_ref_dot, ke, Th, gamma, m):
    """The real closed loop with the positions Θ frozen and the replicator law on the weights: z = (x, β); ẋ = f(x) + ẋ_m − k_ee − f̂(x),
    β̇_i = γβ_ie(Φ̃(x; ϑ_i) − f̂(x)) (n = 1, Σβ = m conserved)."""
    x, beta = z[:1], z[1:]
    e = x - x_ref(t)
    a, c, s = net.physical(Th)
    phi = np.exp(-np.sum((x[None, :] - c) ** 2, axis=1) / (2 * s ** 2))
    Phi = phi * a[:, 0]
    fh = float(beta @ Phi) / m
    xdot = f_vec(x[None, :])[0] + x_ref_dot(t) - ke * e - fh
    return np.concatenate([xdot, gamma * beta * float(e[0]) * (Phi - fh)])


def kl_lyapunov(z, x_ref, t, beta_star, gamma, m):
    """½e² + (1/(γm))Σ_iβ*_i log(β*_i/β_i) (Σβ = Σβ* = m)."""
    e = z[0] - x_ref(t)[0]
    beta = z[1:]
    return 0.5 * float(e ** 2) + float(np.sum(beta_star * np.log(beta_star / beta))) / (gamma * m)
