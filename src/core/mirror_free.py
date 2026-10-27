"""
Proposition D‴.45 (PARADIGM.md §5d): the free mass–amplitude law is a mirror flow — a conservation law per particle, the implicit bias as
a Bregman projection, the noise's invariance in continuous time, the sampled-data drift, the real loop's Lyapunov function
====================================================================================================================================

The replicator–amplitude law of Proposition D‴.34 without the mass constraint (free masses: λ ≡ 0), one particle per class (or the classes
beyond the cohesion's range), each particle driven by one scalar signal g_i — ⟨(r + ξ)φ_i⟩ on a dictionary, eφ_i(x) in the real loop:

    dα_i = Gβ_ia′_ig_i dt,    dβ_i = γβ_ia_ig_i dt,    a_i = a_max tanh(α_i/a_max),   w_i = β_ia_i/m.

(i)   Q_i := β_i − c cosh(2α_i/a_max), c := γa_max²/(4G), is constant for every drive (a noise path read in Stratonovich's sense).
(ii)  On a level set, η_i := ∫dα_i/(Gβ_ia′_i) obeys dη_i = g_i dt and dw_i/dη_i = mD_i, D_i = (β_i/m²)(Gβ_ia′_i² + γa_i²): w follows the
      mirror flow of the separable potential Ψ_q(w) = Σ_iψ_{q_i}(w_i), ψ′ = η, whose conjugate is the mass gained, ψ*(η) = (β − β⁰)/(γm)
      (β⁰ = q + c the mass at α = 0), and ψ(w) = ηw − (β − β⁰)/(γm).
(iii) On a dictionary the noiseless loop's limit, when it reaches an exact fit, is the Bregman projection argmin{Ψ_q(w) − ⟨η(0), w⟩ :
      Φw = y}; noise in continuous time (any colour, any node covariance) moves η only within η(0) + range(Φᵀ), so the projection of the
      current state never moves.
(iv)  An Euler step on (α, log β) changes Q_i by (h²g_i²γβ_i/2)(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²)) + O(h³) (`dQ_euler`; `dQ_exact` the
      whole change); the mirror update η⁺ = η + hg on each particle's curve (`mirror_step`) keeps Q and the projection at every step.
(v)   In the real loop ė = −k_ee + f(x) − f̂(x), f = Σ_iw*_iφ_i, V = ½e² + Σ_i[w*_i(η*_i − η_i) − (β*_i − β_i)/(γm)] (the Bregman divergence,
      (η*, β*) on the particle's own curve at w*_i) has V̇ = −k_ee².
(vi)  The conserved law (Σβ = m) shifts every level: dQ_i = −γβ_iλ dt, λ = Σ_jw_jg_j.

This module: the curves, the mirror coordinates and their inverse (closed forms, the three cases q > c, |q| < c, and the branch q < −c on
which the sign of w is frozen), the potential and its conjugate, the Bregman projection on a dictionary (Newton on the dual), the
loop's vector fields (free and conserved, with a drive), the Euler step's change of Q, and the real loop with its Lyapunov function.
"""
import numpy as np
from scipy.optimize import brentq


class FreeLaw:
    """The free mass–amplitude law's geometry for constants (m, G, γ, a_max); c = γa_max²/(4G)."""

    def __init__(self, m=12.0, G=1.0, gamma=0.1, a_max=1.0):
        self.m, self.G, self.gamma, self.A = float(m), float(G), float(gamma), float(a_max)
        self.c = self.gamma * self.A ** 2 / (4 * self.G)

    # --- the state ---------------------------------------------------------------------------------------------------------------
    def amp(self, al):
        """(a, a′) at α."""
        a = self.A * np.tanh(np.asarray(al, float) / self.A)
        return a, 1.0 - (a / self.A) ** 2

    def Q(self, al, beta):
        """The conserved quantity β − c cosh(2α/a_max)."""
        return np.asarray(beta, float) - self.c * np.cosh(2 * np.asarray(al, float) / self.A)

    def beta(self, al, q):
        """The mass on the level set Q = q at α."""
        return np.asarray(q, float) + self.c * np.cosh(2 * np.asarray(al, float) / self.A)

    def w(self, al, q):
        return self.beta(al, q) * self.amp(al)[0] / self.m

    def relevance(self, al, q):
        """D = (β/m²)(Gβa′² + γa²) on the curve: dw/dη = mD."""
        b = self.beta(al, q)
        a, da = self.amp(al)
        return (b / self.m ** 2) * (self.G * b * da ** 2 + self.gamma * a ** 2)

    def branch_point(self, q):
        """α₀ ≥ 0 with β(±α₀) = 0 when q ≤ −c (the sign of w frozen beyond it), else 0."""
        q = np.asarray(q, float)
        return np.where(q + self.c < 0, 0.5 * self.A * np.arccosh(np.maximum(-q / self.c, 1.0)), 0.0)

    # --- the mirror coordinates ----------------------------------------------------------------------------------------------------
    def eta(self, al, q):
        """η = ∫dα/(Gβa′) along the curve: from α = 0 when q > −c; on a branch (q < −c) a primitive that → −∞ at the branch point.
        Closed form: η = (α − a_max S)/(2Gc), S = κ artanh(κ t) (κ² = (q − c)/(q + c) ∈ (0, 1)), −κ′arctan(κ′t) (κ′² = −κ² > 0),
        κ·½log((κt + 1)/(κt − 1)) (κ > 1, the branch), t = tanh(α/a_max) = a/a_max; for κ near 1 the cancellation is avoided."""
        A, G, c = self.A, self.G, self.c
        al = np.atleast_1d(np.asarray(al, float)).copy()
        q = np.broadcast_to(np.asarray(q, float), al.shape)
        t = np.tanh(al / A)
        out = np.empty_like(al)
        k2 = (q - c) / (q + c)
        ok = q + c > 0
        mid = ok & (k2 > 0)                                                     # q > c: κ ∈ (0, 1)
        if np.any(mid):
            km, tm, am = np.sqrt(k2[mid]), t[mid], al[mid]
            dl = (2 * c / (q[mid] + c)) / (1 + km)                              # 1 − κ, without cancellation
            ktm = km * tm
            stable = np.abs(tm) < 0.9
            # α − Aκ artanh(κt) = (A/2)[log1p(δt/(1 + κt)) − log1p(−δt/(1 − κt)) + δ log((1 + κt)/(1 − κt))] (α = A artanh t)
            v = np.where(stable, 0.5 * A * (np.log1p(dl * tm / (1 + ktm)) - np.log1p(-dl * tm / (1 - ktm))
                                            + dl * np.log((1 + ktm) / (1 - ktm))), 0.0)
            v = np.where(stable, v, am - A * km * np.arctanh(ktm))
            out[mid] = v / (2 * G * c)
        low = ok & (k2 <= 0)                                                    # |q| ≤ c: κ′ = √(−κ²)
        if np.any(low):
            kp = np.sqrt(-k2[low])
            out[low] = (al[low] + A * kp * np.arctan(kp * t[low])) / (2 * G * c)
        br = ~ok                                                                # the branch: κ > 1, |κt| > 1
        if np.any(br):
            kb = np.sqrt(k2[br]); x = kb * t[br]
            with np.errstate(divide='ignore', invalid='ignore'):
                out[br] = (al[br] - A * kb * 0.5 * np.log((x + 1) / (x - 1))) / (2 * G * c)
        return out

    def alpha_of_eta(self, e, q, sgn=None):
        """The inverse of η on each particle's curve (the side of a branch from sgn, default +)."""
        e = np.atleast_1d(np.asarray(e, float)); q = np.broadcast_to(np.asarray(q, float), e.shape)
        sgn = np.ones_like(e) if sgn is None else np.broadcast_to(np.asarray(sgn, float), e.shape)
        a0 = self.branch_point(q)
        out = np.empty_like(e)
        for j in range(e.size):
            f = lambda x: float(self.eta(np.array([x]), q[j])[0]) - e[j]
            if q[j] + self.c > 0:
                lo, hi = -1.0, 1.0
            else:
                eps = 1e-9 * (1.0 + a0[j])
                lo, hi = (a0[j] + eps, a0[j] + 1.0) if sgn[j] >= 0 else (-a0[j] - 1.0, -a0[j] - eps)
            while f(lo) > 0 and lo > -340:
                lo = 2 * lo if lo < 0 else lo - 1.0
            while f(hi) < 0 and hi < 340:
                hi = 2 * hi if hi > 0 else hi + 1.0
            lo, hi = max(lo, -340.0), min(hi, 340.0)
            if f(lo) > 0:
                out[j] = lo
            elif f(hi) < 0:
                out[j] = hi
            else:
                out[j] = brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=500)
        return out

    def psi(self, al, q):
        """ψ_q(w(α)) = ηw − (β − β⁰)/(γm), β⁰ = q + c (q > −c; on a branch up to the constant of η)."""
        return self.eta(al, q) * self.w(al, q) - (self.beta(al, q) - (q + self.c)) / (self.gamma * self.m)

    def psi_star(self, al, q):
        """ψ*_q(η(α)) = (β − β⁰)/(γm): the conjugate potential is the mass gained."""
        return (self.beta(al, q) - (q + self.c)) / (self.gamma * self.m)

    def bregman_div(self, al_star, al, q):
        """B(w*, w) = w*(η* − η) − (β* − β)/(γm) on each curve (η* = η(α*), β* = β(α*)), summed."""
        return float(np.sum(self.w(al_star, q) * (self.eta(al_star, q) - self.eta(al, q))
                            - (self.beta(al_star, q) - self.beta(al, q)) / (self.gamma * self.m)))

    def alpha_of_w(self, w, q, sgn=None):
        """α on the curve with w(α) = w (w is increasing in α on each side)."""
        w = np.atleast_1d(np.asarray(w, float)); q = np.broadcast_to(np.asarray(q, float), w.shape)
        a0 = self.branch_point(q)
        out = np.empty_like(w)
        for j in range(w.size):
            f = lambda x: float(self.w(x, q[j])) - w[j]
            if q[j] + self.c > 0:
                lo, hi = -1.0, 1.0
            else:
                lo, hi = (a0[j], a0[j] + 1.0) if (w[j] >= 0) else (-a0[j] - 1.0, -a0[j])
            while f(lo) > 0 and lo > -340:
                lo = 2 * lo if lo < 0 else lo - 1.0
            while f(hi) < 0 and hi < 340:
                hi = 2 * hi if hi > 0 else hi + 1.0
            out[j] = brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=500)
        return out

    # --- the loop on a dictionary ------------------------------------------------------------------------------------------------
    def rhs(self, Phi, y, xi=None, conserved=False):
        """The loop on (α, log β) on the design Phi (N × d), measured error r + ξ(t); `conserved` adds the replicator's −λ."""
        N, d = Phi.shape
        G, gam, m = self.G, self.gamma, self.m

        def f(t, z):
            al, lb = z[:d], z[d:]
            b = np.exp(lb); a, da = self.amp(al)
            r = y - Phi @ (b * a / m)
            if xi is not None:
                r = r + xi(t)
            g = Phi.T @ r / N
            u = a * g
            lam = (b @ u) / b.sum() if conserved else 0.0
            return np.r_[G * b * da * g, gam * (u - lam)]
        return f

    def bregman_projection(self, Phi, y, q, eta0, sgn=None, tol=1e-13, it_max=300):
        """argmin{Σψ_{q_i}(w_i) − ⟨η0, w⟩ : Φw = y}: w_i = w(α(η0_i + (Φᵀz)_i)), Newton on the dual z (Jacobian Φ diag(mD)Φᵀ), damped.
        Returns (w, α, iterations, the residual's norm)."""
        z = np.zeros(Phi.shape[0])

        def at(z):
            al = self.alpha_of_eta(eta0 + Phi.T @ z, q, sgn)
            return self.w(al, q), al
        w, al = at(z)
        res = Phi @ w - y
        for it in range(it_max):
            if np.linalg.norm(res) < tol:
                break
            step = np.linalg.lstsq(Phi @ ((self.m * self.relevance(al, q))[:, None] * Phi.T), res, rcond=None)[0]
            lam = 1.0
            while lam > 1e-10:
                w2, al2 = at(z - lam * step)
                res2 = Phi @ w2 - y
                if np.linalg.norm(res2) < np.linalg.norm(res):
                    break
                lam /= 2
            z, w, al, res = z - lam * step, w2, al2, res2
        return w, al, it, float(np.linalg.norm(res))

    # --- sampled data -------------------------------------------------------------------------------------------------------------
    def dQ_exact(self, al, beta, g, h):
        """Q's change under one Euler step on (α, log β) with the drive g."""
        a, da = self.amp(al)
        al1 = al + h * self.G * beta * da * g
        b1 = beta * np.exp(h * self.gamma * a * g)
        return self.Q(al1, b1) - self.Q(al, beta)

    def dQ_euler(self, al, beta, g2, h):
        """Its second-order part (h²g²γβ/2)(γa² − Gβa′(1 + a²/a_max²)) (g2 = g², or its mean)."""
        a, da = self.amp(al)
        return 0.5 * h ** 2 * g2 * self.gamma * beta * (self.gamma * a ** 2 - self.G * beta * da * (1 + (a / self.A) ** 2))

    def mirror_step(self, al, q, g, h, iters=8):
        """The mirror update: η⁺ = η + hg, each particle placed on its own curve at η⁺ (Newton in α from the current α, the step
        capped at a_max); the mass is β(α⁺, q) — Q is kept exactly, and η moves by exactly hg."""
        q = np.broadcast_to(np.asarray(q, float), np.shape(al))
        e1 = self.eta(al, q) + h * np.asarray(g, float)
        a = np.array(al, float)
        for _ in range(iters):
            step = (self.eta(a, q) - e1) * self.G * self.beta(a, q) * self.amp(a)[1]
            a = a - np.clip(step, -self.A, self.A)
        return a

    # --- the real loop --------------------------------------------------------------------------------------------------------------
    def real_rhs(self, centres, widths, wstar, x_ref, x_ref_dot, ke, conserved=False):
        """The real loop with Gaussian atoms at fixed centres and widths: z = (x, α, log β), e = x − x_m(t),
        ẋ = f(x) + ẋ_m − k_ee − f̂(x), f = Σw*_iφ_i, f̂ = Σw_iφ_i; α̇_i = Gβ_ia′_ieφ_i(x), (log β_i)˙ = γ(a_ieφ_i(x) − λ)."""
        centres, widths, wstar = (np.asarray(v, float) for v in (centres, widths, wstar))
        d = centres.size
        G, gam, m = self.G, self.gamma, self.m

        def f(t, z):
            x, al, lb = z[0], z[1:1 + d], z[1 + d:]
            b = np.exp(lb); a, da = self.amp(al)
            phi = np.exp(-(x - centres) ** 2 / (2 * widths ** 2))
            e = x - x_ref(t)
            fh = phi @ (b * a / m)
            u = a * e * phi
            lam = (b @ u) / b.sum() if conserved else 0.0
            return np.r_[phi @ wstar + x_ref_dot(t) - ke * e - fh, G * b * da * e * phi, gam * (u - lam)]
        return f

    def real_lyapunov(self, z, t, x_ref, al_star, q):
        """V = ½e² + Σ_iB_{q_i}(w*_i, w_i) at the state z (α* the curves' points at w*)."""
        d = q.size
        return 0.5 * (z[0] - x_ref(t)) ** 2 + self.bregman_div(al_star, z[1:1 + d], q)
