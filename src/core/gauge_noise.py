"""
Proposition D‴.40 (PARADIGM.md §5d): the gauge under general noise
================================================================

At an exact fit with the masses held fixed, white noise of covariance Σ (per unit slow time) on the positions makes them fluctuate about
the fit: to leading order a Gaussian of covariance C, the solution of the Lyapunov equation

    (−GH)C + C(−GH)ᵀ + Σ = 0,   H = ∇²_ΘJ at the fit, G the positions' metric (dϑ = −G∇J ds + noise),

displaced by the mean shift μ = −½H⁻¹∇_Θ tr(H(Θ)C)|_C (the cubic term of J). The masses' fitness u_k = −∂_{β_k}J (zero at the fit)
averaged over these fluctuations is

    ū_k = ∇u_k·μ + ½ tr(∇²u_k C),

and in the slow-replicator limit the masses follow the replicator flow of ū (dW_k/ds = γW_k(ū_k − Σ_lW_lū_l)) — for Σ = 2τ diag(G) the
descent of the entropic free energy (τ/2) log det H (Proposition D‴.37(ii)), for other Σ a different selection, in general not a gradient.

This module: the averaged fitness by central differences (`averaged_fitness`), from callables for ∇_ΘJ and the fitness vector.
"""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov


def averaged_fitness(gradJ, u, theta, Sigma, Gvec, h=1e-4, hH=1e-5):
    """ū (the fitness averaged over the positions' stationary fluctuations at fixed masses), C and μ at the exact fit `theta` (n,): `gradJ(t)`
    the gradient of J in the positions (n,), `u(t)` the fitness vector (K,), `Sigma` the positions' noise covariance (n × n), `Gvec` the
    metric's diagonal (n,). Central differences: steps `hH` for the Hessian of J, `h` for the third derivatives and those of u."""
    theta = np.asarray(theta, float)
    n = len(theta)
    I = np.eye(n)
    H = np.array([(gradJ(theta + hH * e) - gradJ(theta - hH * e)) / (2 * hH) for e in I]).T
    H = 0.5 * (H + H.T)
    A = -(np.asarray(Gvec, float)[:, None] * H)
    C = solve_continuous_lyapunov(A, -np.asarray(Sigma, float))
    C = 0.5 * (C + C.T)
    lam, V = np.linalg.eigh(C)

    def trHC(t):                                                              # tr(H(t)C) = Σ_m λ_m ∂²_{v_m}J(t)
        return sum(l * v @ (gradJ(t + h * v) - gradJ(t - h * v)) / (2 * h) for l, v in zip(lam, V.T) if l != 0.0)
    g = np.array([(trHC(theta + h * e) - trHC(theta - h * e)) / (2 * h) for e in I])
    mu = -0.5 * np.linalg.solve(H, g)
    gu = np.array([(u(theta + h * e) - u(theta - h * e)) / (2 * h) for e in I]).T          # K × n
    K = gu.shape[0]
    Hu = np.zeros((K, n, n))
    for i in range(n):
        for j in range(i, n):
            ei, ej = I[i] * h, I[j] * h
            d = (u(theta + ei + ej) - u(theta + ei - ej) - u(theta - ei + ej) + u(theta - ei - ej)) / (4 * h * h)
            Hu[:, i, j] = d
            Hu[:, j, i] = d
    return gu @ mu + 0.5 * np.einsum('kij,ij->k', Hu, C), C, mu


def replicator_flow(W, ubar):
    """dW_k/ds per unit γ for the fitness ū: W_k(ū_k − Σ_lW_lū_l)."""
    W = np.asarray(W, float)
    return W * (ubar - W @ ubar)
