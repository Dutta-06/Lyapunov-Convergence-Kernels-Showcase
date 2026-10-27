"""
Proposition D‴.48 (PARADIGM.md §5d): the rich transit — from small masses the kernels enter one at a time
====================================================================================================================================

The free mass–amplitude law of Proposition D‴.45 (dα_i = Gβ_ia′_ig_i, dβ_i = γβ_ia_ig_i, w_i = β_ia_i/m, Q_i conserved) started at α = 0
with masses β⁰_i = εb_i follows, in the time τ = t/Λ(ε), Λ = (πa_max/(4G))(2/(cε))^{1/2}, the inverse scale space flow of
Σ_i|w_i|/√b_i (`mirror_free.iss_path`): the kernels enter one at a time. With the kernels' centres moving by the paradigm's law
(Proposition D‴.34's positions, Euclidean: ċ_i = Gβ_ia_i⟨r∂_cφ_i⟩ = Gmw_i⟨r∂_cφ_i⟩) an inactive kernel — coefficient O(ε^{3/2}) — does
not move, and an entered one slides: the limit is the sliding inverse scale space (`sliding_iss`): the inactive kernels' dual coordinates
integrate the residual's correlation at their starting positions, an entrant joins the active block, the block relaxes on the fast time
(on the critical curve Q = −c, an ε-independent flow) to a stationary point of the fit over its coefficients and centres. With a
Wasserstein–Fisher–Rao metric on the positions instead (the mass dropped from the centre's speed, ċ_i = Ga_i⟨r∂_cφ_i⟩) an inactive kernel's
amplitude is O(ε^{1/2}) and its centre moves at O(1) in τ — dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨r∂_cφ⟩, du/dτ = ⟨rφ⟩ (`wfr_inactive_rhs`):
the kernels climb the correlation's peaks before entering — condensation.

This module: Gaussian kernels on a sample grid (`Bumps`; on a reference orbit with the closed loop's error filter, `OrbitBumps` —
Proposition D‴.49), the free law with moving centres under either metric (`moving_rhs`) and with moving widths too
(`moving_rhs_widths`), the active block's fast flow on the critical curve (`critical_rhs`), the sliding inverse scale space
(`sliding_iss`) and the inactive kernels' limit under the Wasserstein–Fisher–Rao metric (`wfr_inactive_rhs`).
"""
import numpy as np
from scipy.integrate import solve_ivp


class Bumps:
    """Gaussian kernels exp(−(x − c)²/(2σ²)) of one width on the samples xs; ⟨f g⟩ is the mean over the samples."""

    def __init__(self, xs, sig):
        self.xs, self.sig = np.asarray(xs, float), float(sig)

    def phi(self, c):
        return np.exp(-(self.xs[:, None] - np.atleast_1d(c)[None, :]) ** 2 / (2 * self.sig ** 2))

    def dphi(self, c):
        """∂_cφ."""
        c = np.atleast_1d(c)
        return self.phi(c) * (self.xs[:, None] - c[None, :]) / self.sig ** 2

    def mean(self, v):
        return np.asarray(v).sum(axis=0) / self.xs.size

    def filt(self, r):
        """The residual as the drive sees it: itself on a dictionary (`OrbitBumps`: the closed loop's error filter)."""
        return r


class OrbitBumps(Bumps):
    """The closed loop's averaged drive (Proposition D‴.49): the samples x_m(t_s) of a T-periodic reference over one period (M points) and
    the error filter L = 1/(s + k_e) — the periodic steady state of ė = −k_ee + r, by FFT —, so that ⟨φ_i filt(r)⟩ = T⁻¹∮φ_i(x_m)L[r(x_m)]dt
    is the averaged drive ⟨eφ_i(x)⟩ of the slow adaptation."""

    def __init__(self, xm, T, M, sig, ke):
        ts = T * np.arange(M) / M
        super().__init__(xm(ts), sig)
        om = 2 * np.pi * np.fft.fftfreq(M, d=T / M)
        self.H, self.ke, self.T = 1.0 / (1j * om + ke), float(ke), float(T)

    def filt(self, r):
        return np.real(np.fft.ifft(np.fft.fft(r) * self.H))


def transit_scale(L, eps):
    """Λ(ε) = (πa_max/(4G))(2/(cε))^{1/2}: λ₁ = Λ·b^{−1/2}(1 + O(ε)) for the mass εb (Proposition D‴.45(ii))."""
    return (np.pi * L.A / (4 * L.G)) * np.sqrt(2 / (L.c * eps))


def moving_rhs(L, K, fstar, metric='paradigm'):
    """The free law with moving centres, z = (α, log β, c): α̇ = Gβa′g, (log β)˙ = γag, ċ = Gβa⟨r∂_cφ⟩ (metric 'paradigm') or Ga⟨r∂_cφ⟩
    ('wfr'); g = ⟨rφ⟩, r = f* − Σ_i(β_ia_i/m)φ_i."""
    if metric not in ('paradigm', 'wfr'):
        raise ValueError(metric)

    def f(t, z):
        d = z.size // 3
        al, lb, c = z[:d], z[d:2 * d], z[2 * d:]
        b = np.exp(lb); a, da = L.amp(al)
        Ph = K.phi(c)
        r = K.filt(fstar - Ph @ (b * a / L.m))
        g = Ph.T @ r / K.xs.size
        h = K.dphi(c).T @ r / K.xs.size
        speed = L.G * b * a if metric == 'paradigm' else L.G * a
        return np.r_[L.G * b * da * g, L.gamma * a * g, speed * h]
    return f


def moving_rhs_widths(L, K, fstar):
    """The free law with moving centres and widths (the paradigm's metric), z = (α, log β, c, σ): ċ = Gβa⟨r∂_cφ⟩, σ̇ = Gβa⟨r∂_σφ⟩,
    φ_i = exp(−(x − c_i)²/(2σ_i²)) on K's samples (K's own width unused)."""
    def f(t, z):
        d = z.size // 4
        al, lb, c, sg = z[:d], z[d:2 * d], z[2 * d:3 * d], z[3 * d:]
        b = np.exp(lb); a, da = L.amp(al)
        dx = K.xs[:, None] - c[None, :]
        Ph = np.exp(-dx ** 2 / (2 * sg ** 2))
        r = K.filt(fstar - Ph @ (b * a / L.m))
        g = Ph.T @ r / K.xs.size
        h = (Ph * dx / sg ** 2).T @ r / K.xs.size
        hs = (Ph * dx ** 2 / sg ** 3).T @ r / K.xs.size
        sp = L.G * b * a
        return np.r_[L.G * b * da * g, L.gamma * a * g, sp * h, sp * hs]
    return f


def critical_w(L, al, b0=0.0):
    """The coefficient on the curve from α = 0 with the mass b0 (b0 = 0: the critical curve Q = −c): β = b0 + 2c sinh²(α/a_max) — Q kept
    exactly, without cancellation —, w = βa/m (≈ 2cα³/(ma_max²) at small α on the critical curve)."""
    return (b0 + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m


def curve_mobility(L, al, b0=0.0):
    """dα/dη = Gβa′ on the curve from α = 0 with the mass b0: G(b0 sech²(α/a_max) + 2c tanh²(α/a_max)), no cancellation at large α."""
    x = np.minimum(np.abs(al) / L.A, 350.0)
    return L.G * (b0 / np.cosh(x) ** 2 + 2 * L.c * np.tanh(x) ** 2)


def critical_rhs(L, K, fstar, metric='paradigm', b0=0.0):
    """The free law with moving centres on its curves, z = (α, c): α̇ = Gβa′g, ċ = Gβa⟨r∂_cφ⟩ (metric 'paradigm') or Ga⟨r∂_cφ⟩ ('wfr'),
    β = b0 + 2c sinh²(α/a_max) — each kernel kept on its curve exactly (round 37: integrating (α, log β) drifts off Q at small masses).
    b0 = 0: the active block's fast flow on the critical curve; b0 = εb: the whole loop from α = 0 with the masses εb."""
    if metric not in ('paradigm', 'wfr'):
        raise ValueError(metric)

    def f(t, z):
        d = z.size // 2
        al, c = z[:d], z[d:]
        b = b0 + 2 * L.c * np.sinh(al / L.A) ** 2
        a = L.A * np.tanh(al / L.A)
        Ph = K.phi(c)
        r = K.filt(fstar - Ph @ (b * a / L.m))
        g = Ph.T @ r / K.xs.size
        h = K.dphi(c).T @ r / K.xs.size
        speed = L.G * b * a if metric == 'paradigm' else L.G * a
        return np.r_[curve_mobility(L, al, b0) * g, speed * h]
    return f


def sliding_iss(L, K, fstar, c0, rho, seed=1e-7, t_fast=1e5, tol=1e-9, exit_tol=1e-6, max_events=60, tie_tol=1e-9, gtol=1e-9,
                max_doublings=14):
    """The sliding inverse scale space (the paradigm's metric). The inactive kernels stay where they are (their start, or where they
    left the block), their u_j growing at ⟨rφ(c_j)⟩ (r the active block's residual); the entrants — every kernel whose |u_j| reaches ρ_j
    within the relative tie_tol of the first, so that exact ties enter together — join the active block with the coefficient `seed` (the
    sign of u_j); the block relaxes by `critical_rhs` (the fast flow on the critical curve) over t_fast, then over the elapsed time again
    (doubling it) until it is stationary — on every kernel of the block with a coefficient above exit_tol·max, |⟨rφ_i⟩| and
    |w_i⟨r∂_cφ_i⟩| below gtol·max_j|⟨f*φ_j⟩|; on the others, no outward correlation (s_i⟨rφ_i⟩ ≤ gtol·max) — or max_doublings is reached.
    A kernel of the block whose coefficient has fallen below exit_tol·max leaves, frozen where it is, only if its correlation points into
    its box (s_i⟨rφ_i⟩ < 0); one with an outward correlation stays on its wall. Stops when the rms residual is below tol. Returns
    dict(events = [(τ, j, 'enter'/'exit', c_j)], hist = [(τ, c, w, rms)], c, w, act, relaxed = [whether each relaxation became
    stationary])."""
    d = len(c0)
    rho = np.broadcast_to(np.asarray(rho, float), (d,))
    u = np.zeros(d); c = np.array(c0, float); al = np.zeros(d)
    act = np.zeros(d, bool)
    tau, ev, hist, relaxed = 0.0, [], [], []
    nx = K.xs.size
    scale = float(np.abs(K.phi(c).T @ K.filt(fstar)).max()) / nx
    rhs = critical_rhs(L, K, fstar)
    for _ in range(max_events):
        w = np.where(act, critical_w(L, al), 0.0)
        r = fstar - K.phi(c) @ w
        rms = float(np.sqrt(np.mean(r ** 2)))
        hist.append((tau, c.copy(), w.copy(), rms))
        if rms < tol:
            break
        v = K.phi(c).T @ K.filt(r) / nx
        dt = np.full(d, np.inf)
        ina = ~act & (v != 0)
        dt[ina] = (rho[ina] - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
        dmin = dt.min()
        if not np.isfinite(dmin):
            break
        ent = np.flatnonzero(dt <= dmin * (1 + tie_tol) + tie_tol * np.finfo(float).tiny)
        u[~act] += dmin * v[~act]
        tau += dmin
        for j in ent:
            u[j] = np.sign(v[j]) * rho[j]
            act[j] = True
            al[j] = np.sign(v[j]) * (seed * L.m * L.A ** 2 / (2 * L.c)) ** (1 / 3)
            ev.append((tau, int(j), 'enter', float(c[j])))
        idx = np.flatnonzero(act)
        n = idx.size
        z = np.r_[al[idx], c[idx]]
        elapsed, chunk, ok = 0.0, float(t_fast), False
        for _k in range(max_doublings + 1):
            sol = solve_ivp(rhs, (0, chunk), z, method='LSODA', rtol=1e-10, atol=1e-13)
            z = sol.y[:, -1]
            elapsed += chunk
            wi = critical_w(L, z[:n])
            ri = K.filt(fstar - K.phi(z[n:]) @ wi)
            g = K.phi(z[n:]).T @ ri / nx
            h = K.dphi(z[n:]).T @ ri / nx
            big = np.abs(wi) > exit_tol * np.abs(wi).max()
            si = np.sign(u[idx])
            if (np.all(np.abs(g[big]) <= gtol * scale) and np.all(np.abs(wi[big] * h[big]) <= gtol * scale)
                    and np.all(si[~big] * g[~big] <= gtol * scale)):
                ok = True
                break
            chunk = elapsed
        relaxed.append(ok)
        al[idx] = z[:n]; c[idx] = z[n:]
        for k, i in enumerate(idx):
            if not big[k] and si[k] * g[k] < 0:
                act[i] = False; al[i] = 0.0
                ev.append((tau, int(i), 'exit', float(c[i])))
    return dict(events=ev, hist=hist, c=c, w=np.where(act, critical_w(L, al), 0.0), act=act, relaxed=relaxed)


def wfr_inactive_rhs(L, K, fstar, rho):
    """The inactive kernels' limit under the Wasserstein–Fisher–Rao metric before any entry (r = f*), z = (u, c), masses εb (ρ = b^{−1/2}):
    in the core a ≈ α = a_max(εb/(2c))^{1/2}tan(πη/(2λ₁)), so du/dτ = ⟨f*φ(c)⟩, dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ(c)⟩."""
    rho = np.asarray(rho, float)
    k = np.pi * L.A ** 2 / (4 * L.c)

    def f(t, z):
        d = z.size // 2
        u, c = z[:d], z[d:]
        g = K.phi(c).T @ fstar / K.xs.size
        h = K.dphi(c).T @ fstar / K.xs.size
        return np.r_[g, k * np.tan(np.pi * u / (2 * rho)) * h / rho]
    return f
