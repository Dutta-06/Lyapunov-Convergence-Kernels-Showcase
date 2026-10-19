"""
Perfect tracking with an adaptive prior in the real closed loop: the Floquet spectrum (Corollary D‴.26, PARADIGM.md §5d)
=========================================================================================================================

The real closed loop (CL) of Theorem J with the nearest-atom leak toward an adaptive prior: the particles leak at the rate σ = σ₀/k_e toward
the atom z_{π(i)} of their class, and the atoms move toward the class means at the rate η/k_e,

    ẋ = f(x) + ẋ_m − k_e(x − x_m) − f̂_Θ(x),   ϑ̇_i = G∇_ϑ(eᵀΦ̃(x; ϑ_i)) − (σ₀/k_e)(ϑ_i − z_{π(i)}),   ż_j = −(η/k_e)(z_j − ϑ̄_j),

ϑ̄_j the mean of class j (η = ∞: z_j = ϑ̄_j, eliminated). The classes are frozen: exact near a replicated configuration whose particles lie
inside their cells. With Θ the lift of an exact representation ϑ^e(w) of the sampled weights and Z = ϑ^e(w), the plant on the reference is
a T₀-periodic solution with e ≡ 0 (Proposition D‴.24(vi)); this module gives the loop with the prior in the state, the lift, and the
predicted slow Floquet exponents (1/k_e)·{σ₀ on the separations, λ_j^±(w) on the collective pairs} of Theorem D‴.23(ii) at the weights w.
`lifted_orbit.monodromy_fd` and `periodic.floquet_exponents` compute the actual ones.
"""
import numpy as np
from scipy.integrate import solve_ivp

from .periodic import ClosedLoop
from .adaptive_prior import pair_rates


class AdaptiveRealLoop(ClosedLoop):
    """(CL) with the nearest-atom leak toward an adaptive prior; the state is (x, Θ, Z) for finite η and (x, Θ) for η = ∞; `labels` (m,)
    the frozen classes, K0 the number of atoms; `weights` (m,) the weights of the particles in f̂ (default 1/m; the weighted K₀-atom loop has
    one particle per class with the weights w)."""

    def __init__(self, net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0, labels, K0, eta, weights=None):
        super().__init__(net, f_vec, x_ref, x_ref_dot, period, ke, G_diag, sigma0)
        self.labels = np.asarray(labels, int)
        self.K0 = int(K0)
        self.eta = float(eta)
        self.counts = np.bincount(self.labels, minlength=self.K0)
        assert np.all(self.counts > 0)
        self.weights = None if weights is None else np.asarray(weights, float)

    def f_hat(self, x, Th):
        if self.weights is None:
            return self.net.f_hat(x, Th)
        a, c, s = self.net.physical(Th)
        phi = np.exp(-np.sum((x[None, :] - c) ** 2, axis=1) / (2 * s ** 2))
        return (self.weights[:, None] * phi[:, None] * a).sum(axis=0)

    def class_means(self, Th):
        out = np.zeros((self.K0, self.p))
        np.add.at(out, self.labels, Th)
        return out / self.counts[:, None]

    def state_dim(self, m):
        return self.net.n + m * self.p + (0 if np.isinf(self.eta) else self.K0 * self.p)

    def rhs(self, t, z, m):
        n, p = self.net.n, self.p
        x, Th = z[:n], z[n:n + m * p].reshape(m, p)
        Z = self.class_means(Th) if np.isinf(self.eta) else z[n + m * p:].reshape(self.K0, p)
        xr, xrd = self.xr(t), self.xrd(t)
        e = x - xr
        fx = self.f_vec(x[None, :])[0]
        xdot = fx + xrd - self.ke * e - self.f_hat(x, Th)
        Thdot = self.G[None, :] * self.net.grad_g(x, e, Th) - (self.sigma0 / self.ke) * (Th - Z[self.labels])
        if np.isinf(self.eta):
            return np.concatenate([xdot, Thdot.ravel()])
        Zdot = -(self.eta / self.ke) * (Z - self.class_means(Th))
        return np.concatenate([xdot, Thdot.ravel(), Zdot.ravel()])

    def lift(self, atoms):
        """The state of perfect tracking at t = 0: x = x_m(0), the lift of the atoms by the classes, Z = the atoms."""
        atoms = np.asarray(atoms, float)
        parts = [np.atleast_1d(self.xr(0.0)), atoms[self.labels].ravel()]
        if not np.isinf(self.eta):
            parts.append(atoms.ravel())
        return np.concatenate(parts)


def predicted_slow_exponents(h, sigma0, eta, ke, m, K0, p):
    """The predicted slow Floquet exponents (rates per unit time, ascending) at a lifted exact representation: σ₀/k_e with multiplicity
    (m − K₀)p (the separations), and λ_j^±/k_e for the K₀p collective pairs (η = ∞: h_j/k_e only, the prior eliminated)."""
    h = np.asarray(h, float)
    lm, lp = pair_rates(h, sigma0, eta)
    coll = list(lm) if np.isinf(eta) else list(lm) + list(lp)
    sep = [sigma0] * ((m - K0) * p)
    return np.sort(np.array(coll + sep)) / ke


def period_map(cl, z0, m, rtol=1e-13, atol=1e-15):
    """One period of the loop from t = 0 (DOP853 at tight tolerances)."""
    sol = solve_ivp(lambda t, z: cl.rhs(t, z, m), (0.0, cl.T0), np.asarray(z0, float), method='DOP853', rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1]


def monodromy(cl, z_p, m, fd=1e-5, rtol=1e-13, atol=1e-15):
    """The monodromy matrix at the periodic state z_p by central differences of the period map, with the relative step fd·max(1, |z_j|)
    (truncation O(fd²); the integrator's error enters divided by the step, so a larger step and tight tolerances resolve the multipliers
    near 1 far better than lifted_orbit.monodromy_fd's 10⁻⁷ step at rtol 10⁻¹¹)."""
    d = len(z_p)
    DP = np.zeros((d, d))
    for j in range(d):
        h = fd * max(1.0, abs(z_p[j]))
        e = np.zeros(d); e[j] = h
        DP[:, j] = (period_map(cl, z_p + e, m, rtol, atol) - period_map(cl, z_p - e, m, rtol, atol)) / (2 * h)
    return DP

