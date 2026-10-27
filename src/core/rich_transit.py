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
Proposition D‴.49's leading order —, and the closed loop's exact averaged drive along the frozen periodic orbit, `FrozenOrbit` with
`frozen_orbit_rhs`), the free law with moving centres under either metric (`moving_rhs`) and with moving widths too
(`moving_rhs_widths`), the active block's fast flow on the critical curve (`critical_rhs`), the sliding inverse scale space
(`sliding_iss`) and the inactive kernels' limit under the Wasserstein–Fisher–Rao metric (`wfr_inactive_rhs`). Proposition D‴.51 (the
transit at a finite mass): the critical curve's inner profile near the wall (`eta_crit`, `coef_crit`, `inner_coef`, `inner_xi`), the
entrant's overshoot (`overshoot`), the ε-free inner problem of a near tie (`near_tie_inner`), the sliding path's duals (`path_duals`), the
creep of a kernel left near its wall along the path (`creep_reaches`) and the mass below which the loop forms exactly the path's kernels,
predicted from the path alone (`recovery_threshold`).
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


class FrozenOrbit:
    """The closed loop's exact averaged drive (Proposition D‴.49): with the parameters frozen, the loop ė = −k_ee + r(x_m + e),
    r = f − f̂, has for k_e > Lip(r) a unique, exponentially stable T-periodic solution e_p = L[r(x_m + e_p)] (L = 1/(s + k_e); the map is
    a contraction of the periodic functions with factor Lip(r)/k_e) — computed on M samples of a period by the fixed point through the FFT,
    warm-started from the last call. The slow adaptation averages the regressor along the frozen orbit x_p = x_m + e_p: g_i = ⟨φ_i(x_p)e_p⟩,
    h_i = ⟨∂_cφ_i(x_p)e_p⟩ (the period's mean). `OrbitBumps` is its truncation at x_p = x_m, e_p = L[r(x_m)]."""

    def __init__(self, xm, T, M, sig, ke, f, tol=1e-12, it_max=2000):
        self.ts = T * np.arange(M) / M
        self.xm = xm(self.ts)
        om = 2 * np.pi * np.fft.fftfreq(M, d=T / M)
        self.H = 1.0 / (1j * om + ke)
        self.sig, self.ke, self.f, self.tol, self.it_max = float(sig), float(ke), f, tol, it_max
        self.e = np.zeros(M)
        self.iters = 0

    def orbit(self, c, w):
        """The frozen orbit's error e_p (raises if the fixed point does not converge — k_e too small against Lip(r))."""
        e = self.e
        for k in range(self.it_max):
            x = self.xm + e
            ph = np.exp(-(x[:, None] - c[None, :]) ** 2 / (2 * self.sig ** 2))
            en = np.real(np.fft.ifft(np.fft.fft(self.f(x) - ph @ w) * self.H))
            dif = np.abs(en - e).max()
            e = en
            if dif <= self.tol * max(1.0, np.abs(e).max()):
                break
        else:
            raise RuntimeError('FrozenOrbit: the periodic orbit did not converge (k_e ≤ Lip(f − f̂)?)')
        self.iters += k + 1
        self.e = e
        return e

    def drive(self, c, w):
        """(g, h, e_p): the averaged correlations ⟨φ_i(x_p)e_p⟩, ⟨∂_cφ_i(x_p)e_p⟩ and the orbit's error."""
        e = self.orbit(c, w)
        x = self.xm + e
        dx = x[:, None] - c[None, :]
        ph = np.exp(-dx ** 2 / (2 * self.sig ** 2))
        return ph.T @ e / e.size, (ph * dx / self.sig ** 2).T @ e / e.size, e


def frozen_orbit_rhs(L, O, b0, metric='paradigm'):
    """The averaged closed loop on the free law's curves from α = 0 with the masses b0 (β = b0 + 2c sinh²(α/a_max), Q kept), z = (α, c):
    α̇ = Gβa′g, ċ = Gβa h (metric 'paradigm'), Ga h ('wfr') or 0 ('fixed'); (g, h) the frozen orbit's drive (`FrozenOrbit`)."""
    if metric not in ('paradigm', 'wfr', 'fixed'):
        raise ValueError(metric)
    b0 = np.asarray(b0, float)

    def f(t, z):
        d = z.size // 2
        al, c = z[:d], z[d:]
        b = b0 + 2 * L.c * np.sinh(al / L.A) ** 2
        a = L.A * np.tanh(al / L.A)
        g, h, _ = O.drive(c, b * a / L.m)
        speed = L.G * b * a if metric == 'paradigm' else (L.G * a if metric == 'wfr' else 0.0 * a)
        return np.r_[curve_mobility(L, al, b0) * g, speed * h]
    return f


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


def polish_block(K, fstar, w, c, s):
    """The active block's fit from (w, c) by a sign-constrained least-squares solve (trust-region reflective): min ½⟨(f* − Φ(c)w)²⟩ over
    (w, c) with s_iw_i ≥ 0. On a dictionary (K.filt the identity) the fast flow is a gradient flow of this fit, so its end is a stationary
    point of it; the solve reaches the point where the flow crawls along the block's slow modes (two kernels Δ apart: the fit's Hessian has
    eigenvalues ∝ 1, 1, Δ², Δ⁶ — Proposition D‴.52). Returns (w, c)."""
    from scipy.optimize import least_squares
    n, nx = w.size, K.xs.size
    lo = np.r_[np.where(s > 0, 0.0, -np.inf), np.full(n, -np.inf)]
    hi = np.r_[np.where(s > 0, np.inf, 0.0), np.full(n, np.inf)]
    z0 = np.r_[np.where(s > 0, np.maximum(w, 1e-300), np.minimum(w, -1e-300)), c]
    sol = least_squares(lambda z: (fstar - K.phi(z[n:]) @ z[:n]) / np.sqrt(nx), z0,
                        jac=lambda z: -np.c_[K.phi(z[n:]), K.dphi(z[n:]) * z[:n][None, :]] / np.sqrt(nx),
                        bounds=(lo, hi), method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
    return sol.x[:n], sol.x[n:]


def sliding_iss(L, K, fstar, c0, rho, seed=1e-7, t_fast=1e5, tol=1e-9, exit_tol=1e-6, max_events=60, tie_tol=1e-9, gtol=1e-9,
                max_doublings=14, polish=False):
    """The sliding inverse scale space (the paradigm's metric). The inactive kernels stay where they are (their start, or where they
    left the block), their u_j growing at ⟨rφ(c_j)⟩ (r the active block's residual); the entrants — every kernel whose |u_j| reaches ρ_j
    within the relative tie_tol of the first, so that exact ties enter together — join the active block with the coefficient `seed` (the
    sign of u_j); the block relaxes by `critical_rhs` (the fast flow on the critical curve) over t_fast, then over the elapsed time again
    (doubling it) until it is stationary — on every kernel of the block with a coefficient above exit_tol·max, |⟨rφ_i⟩| and
    |w_i⟨r∂_cφ_i⟩| below gtol·max_j|⟨f*φ_j⟩|; on the others, no outward correlation (s_i⟨rφ_i⟩ ≤ gtol·max) — or max_doublings is reached.
    A kernel of the block whose coefficient has fallen below exit_tol·max leaves, frozen where it is, only if its correlation points into
    its box (s_i⟨rφ_i⟩ < 0); one with an outward correlation stays on its wall. Stops when the rms residual is below tol. Returns
    dict(events = [(τ, j, 'enter'/'exit', c_j)], hist = [(τ, c, w, rms)], c, w, act, relaxed = [whether each relaxation became
    stationary]). polish=True (a dictionary only): after the first t_fast the block is moved to its fit's stationary point by
    `polish_block` — the flow's end when no other stationary point lies on its way — and checked there, the flow resuming if needed."""
    d = len(c0)
    rho = np.broadcast_to(np.asarray(rho, float), (d,))
    u = np.zeros(d); c = np.array(c0, float); al = np.zeros(d)
    act = np.zeros(d, bool)
    tau, ev, hist, relaxed = 0.0, [], [], []
    nx = K.xs.size
    scale = float(np.abs(K.phi(c).T @ K.filt(fstar)).max()) / nx
    rhs = critical_rhs(L, K, fstar)
    if polish and type(K).filt.__qualname__ != 'Bumps.filt':
        raise ValueError('polish needs a dictionary (the identity filter)')
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
            if polish and _k == 0:
                wp, cp = polish_block(K, fstar, critical_w(L, z[:n]), z[n:], np.sign(u[idx]))
                z = np.r_[np.sign(wp) * _alpha_of(L, np.abs(wp), coef_crit), cp]
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


# ------------------------------------------------------------------------------------------------------------------------------------
# Proposition D‴.51 (PARADIGM.md §5d): the transit at a finite mass — near ties, the entrant's overshoot, the creep, the predicted ε
# ------------------------------------------------------------------------------------------------------------------------------------

def eta_crit(L, al):
    """The critical curve's mirror coordinate relative to the wall, η_c(α) = lim_{ε→0}(η_ε(α) − λ₁) = (α − a coth(α/a))/(2Gc):
    −a²/(2Gcα) near α = 0 (the core's tangent), (α − a)/(2Gc) for α ≫ a (the exponential branch)."""
    x = np.asarray(al, float) / L.A
    return (L.A / (2 * L.G * L.c)) * (x - 1 / np.tanh(x))


def coef_crit(L, al):
    """The coefficient w = βa/m on the critical curve β = 2c sinh²(α/a), a = a_max tanh(α/a) (|α|, overflow-safe)."""
    x = np.minimum(np.abs(np.asarray(al, float)) / L.A, 350.0)
    return 2 * L.c * np.sinh(x) ** 2 * L.A * np.tanh(x) / L.m


def _alpha_of(L, target, fn, lo=1e-12, hi=1e3, it=200):
    """α > 0 with fn(L, α) = target (fn increasing), by bisection in log α, vectorised."""
    t = np.asarray(target, float)
    a, b = np.full(t.shape, np.log(lo)), np.full(t.shape, np.log(hi))
    for _ in range(it):
        m_ = 0.5 * (a + b)
        up = fn(L, np.exp(m_)) < t
        a, b = np.where(up, m_, a), np.where(up, b, m_)
    return np.exp(0.5 * (a + b))


def inner_coef(L, xi, rho=1.0):
    """The inner profile: a kernel whose dual coordinate is u = ρ(1 − ε^{1/2}ξ) carries, as ε → 0 at fixed ξ, the coefficient
    𝒲(ξ) = w_c(α) with η_c(α) = −Λ₀ρξ, Λ₀ = Λ(ε)ε^{1/2} = (πa/(4G))(2/c)^{1/2}: ≈ K'³ρ⁻³ξ⁻³ for ξ ≫ 1 (the core's profile,
    K'³ = (8/π³)a/(m(2c)^{1/2})), ∝ exp(−πρ(2c)^{1/2}ξ) for ξ ≪ −1 (the exponential branch)."""
    lam0 = transit_scale(L, 1.0)
    al = _alpha_of(L, -lam0 * rho * np.asarray(xi, float), eta_crit)
    return coef_crit(L, al)


def inner_xi(L, w, rho=1.0):
    """The inverse of `inner_coef`: the inner dual coordinate ξ of a kernel carrying |w| near its wall."""
    al = _alpha_of(L, np.abs(np.asarray(w, float)), coef_crit)
    return -eta_crit(L, al) / (transit_scale(L, 1.0) * rho)


def overshoot(L, A, rho=1.0):
    """The entrant's overshoot: an entered kernel carrying the coefficient A has its dual coordinate beyond its wall by
    η_c(α_A)/Λ = ε^{1/2}·(−inner_xi(A)) — what a runner-up's dual is pulled along by while the entrant grows (Proposition D‴.51(i))."""
    return -inner_xi(L, A, rho)


def near_tie_inner(L, K, A, cstar, c_e, c_r, x, rho=1.0, xi0=400.0, s_end=5000.0, both=True, tol=1e-9, max_doublings=8):
    """The ε-free inner problem of a near tie (Proposition D‴.51(i)): an entrant at c_e and a runner-up at c_r, a bump A at cstar, in
    s = (τ − τ_e)/ε^{1/2} and ξ_i = (ρ − s_Au_i)/(ρε^{1/2}) (s_A the bump's sign): dξ_i/ds = −s_A⟨rφ_i⟩/ρ, dc_i/ds = Λ₀Gm w_i⟨r∂_cφ_i⟩,
    w_i = s_A𝒲(ξ_i), r = Aφ(c*) − Σ_iw_iφ(c_i); the runner-up starts x behind the entrant (ξ_r − ξ_e = x, from ξ_e = xi0). both=False
    freezes the runner-up (its coefficient dropped): the entrant alone, the runner-up's ξ then recording the pull. The start xi0 behind the
    wall truncates the entrant's approach, where its core coefficient 𝒲 ∝ ξ⁻³ already acts on the runner-up: an error ∝ xi0⁻² in the
    pull. Integrated over the entrant's outer arrival time s₀ = xi0ρ/|⟨fφ_e⟩| plus s_end, then over the elapsed time again (doubling) until
    stationary — every rate below tol·|A| (a small |A| relaxes slowly) — or max_doublings is reached. Returns the end state with the
    runner-up's share |w_r|/(|w_e| + |w_r|), whether it became stationary, and (both=False) the matched pull — the runner-up's outer margin
    at the entrant's arrival, x + s₀(d_e − d_r) (d the drifts before the entry), less ξ_r at the end, its drift after the entry (the
    end residual's) subtracted over s − s₀ (`entry_pull`'s definition; = x − ξ_r(∞) at an exact tie onto a fitted bump)."""
    lam0 = transit_scale(L, 1.0)
    sA = np.sign(A)
    f_ = A * K.phi(np.array([cstar]))[:, 0]
    nx = K.xs.size

    def W(xi):
        return float(inner_coef(L, xi, rho))

    def rhs(s, z):
        xe, xr, ce, cr = z
        we = sA * W(xe)
        wr = sA * W(xr) if both else 0.0
        pe, pr = K.phi(np.array([ce]))[:, 0], K.phi(np.array([cr]))[:, 0]
        r = f_ - we * pe - wr * pr
        de, dr = K.dphi(np.array([ce]))[:, 0], K.dphi(np.array([cr]))[:, 0]
        return [-sA * (pe @ r) / nx / rho, -sA * (pr @ r) / nx / rho, lam0 * L.G * L.m * we * (de @ r) / nx,
                lam0 * L.G * L.m * wr * (dr @ r) / nx]
    pe0, pr0 = K.phi(np.array([c_e]))[:, 0], K.phi(np.array([c_r]))[:, 0]
    d_e, d_r = sA * (pe0 @ f_) / nx / rho, sA * (pr0 @ f_) / nx / rho
    s_arr = xi0 / d_e
    z, s0, chunk, ok = np.array([xi0, xi0 + x, c_e, c_r], float), 0.0, s_arr + float(s_end), False
    for _ in range(max_doublings + 1):
        sol = solve_ivp(rhs, (s0, s0 + chunk), z, method='LSODA', rtol=1e-10, atol=1e-12)
        z, s0 = sol.y[:, -1], s0 + chunk
        if np.all(np.abs(rhs(s0, z)) <= tol * abs(A)):
            ok = True
            break
        chunk = s0
    xe, xr, ce, cr = z
    we, wr = W(xe), (W(xr) if both else 0.0)
    d_end = -rhs(s0, z)[1]
    pull = (x + s_arr * (d_e - d_r)) - xr - d_end * (s0 - s_arr) if not both else np.nan
    return dict(xi_e=float(xe), xi_r=float(xr), c_e=float(ce), c_r=float(cr), w_e=we, w_r=wr, share=wr / (we + wr),
                status=int(sol.status), stationary=ok, s=s0, pull=float(pull), s_arrival=float(s_arr))


def path_duals(K, fstar, R):
    """The sliding path's (`sliding_iss`) dual coordinates reconstructed at its event times: on each interval the block's residual r_k is
    fixed, the inactive kernels' u grows at ⟨r_kφ(c_j)⟩, an entrant is pinned at its wall (±1 in units of ρ — equal masses). Returns
    dict(taus, res (r_k after the k-th relaxation), centres, U (the duals after the events at taus[k]), ever (the kernels that enter))."""
    hist = R['hist']
    taus = [h[0] for h in hist]
    res = [fstar - K.phi(h[1]) @ h[2] for h in hist]
    cs = [h[1] for h in hist]
    d = len(cs[0])
    nx = K.xs.size
    ever = sorted({j for t, j, k, c in R['events'] if k == 'enter'})
    u = np.zeros(d); act = np.zeros(d, bool); U = [u.copy()]
    ev = R['events']; e = 0
    for k in range(len(taus) - 1):
        v = K.phi(cs[k]).T @ K.filt(res[k]) / nx
        u = u + np.where(act, 0.0, (taus[k + 1] - taus[k]) * v)
        while e < len(ev) and ev[e][0] <= taus[k + 1] + 1e-12:
            t, j, kind, c = ev[e]
            if kind == 'enter':
                act[j] = True; u[j] = np.sign(u[j])
            else:
                act[j] = False
            e += 1
        U.append(u.copy())
    return dict(taus=taus, res=res, centres=cs, U=np.array(U), ever=ever)


def creep_reaches(L, K, P, j, eps, c0j, b=1.0):
    """Whether kernel j (never entering on the path), started at c0j and integrated along the path's residuals with its creep law
    dc/dτ = ε(πa²/(4c))b^{3/2}cos(πμ/2)/sin³(πμ/2)·s⟨r∂_cφ⟩ (Proposition D‴.48(iii)(a)) and du/dτ = ⟨rφ(c)⟩, comes within the inner
    layer 3ε^{1/2} of its wall before the path's last event (P from `path_duals`)."""
    pref = np.pi * L.A ** 2 / (4 * L.c) * b ** 1.5
    nx = K.xs.size
    mu_in = 3 * np.sqrt(eps)
    z = np.array([0.0, c0j])
    for k in range(len(P['taus']) - 1):
        r = K.filt(P['res'][k])

        def f(t, z):
            u, c = z
            mu = max(1 - abs(u), mu_in)
            ph, dph = K.phi(np.array([c]))[:, 0], K.dphi(np.array([c]))[:, 0]
            return [ph @ r / nx, eps * pref * np.cos(np.pi * mu / 2) / np.sin(np.pi * mu / 2) ** 3 * np.sign(u) * (dph @ r) / nx]

        def hit(t, z):
            return 1 - abs(z[0]) - mu_in
        hit.terminal = True
        sol = solve_ivp(f, (P['taus'][k], P['taus'][k + 1]), z, method='LSODA', rtol=1e-8, atol=1e-11, events=hit)
        if sol.status == 1:
            return True
        z = sol.y[:, -1]
    return False


def entry_pull(L, K, fstar, R, P, k, i, tracers, rho=1.0, xi0=400.0, s_end=200.0, tol=1e-10, max_doublings=12):
    """The pull an entry exerts on kernels that do not enter, measured on the entry's ε-free inner problem (Proposition D‴.51(iv)): in
    s = (τ − τ_e)/ε^{1/2} the active block of the path before the event k (index into P's times) and the entrant i follow the loop on the
    critical curve — dα/ds = Λ₀·2Gc tanh²(α/a)⟨rφ⟩, dc/ds = Λ₀Gm w⟨r∂_cφ⟩, w = 2c sinh²(α/a)·a tanh(α/a)/m —, the block from its
    coefficients, the entrant from the inner coordinate xi0 behind its wall; the tracers j stay where they are (inactive kernels do not
    move), their inner coordinates ξ_j = μ_j/ε^{1/2} moving at dξ_j/ds = −s_j⟨rφ_j⟩/ρ. The pull Δ_j = lim_s (ξ_j^out(s) − ξ_j(s)), ξ_j^out
    the outer solution — the drift of the residual before the entry until the entrant's outer arrival at its wall, s₀ = xi0ρ/(s_i⟨r⁻φ_i⟩),
    the inner problem's final residual's after: kernel j forms at the mass ε iff μ_j/ε^{1/2} − Δ_j < ξ_θ (𝒲(ξ_θ) = θ). Integrated over
    s₀ + s_end, then doubling until the block's rates fall below tol (or max_doublings). Returns dict(pull = {j: Δ_j}, stationary, s,
    w_end (the block's and the entrant's coefficients at the end))."""
    lam0 = transit_scale(L, 1.0)
    nx = K.xs.size
    hist = R['hist']
    cb, wb = hist[k - 1][1], hist[k - 1][2]
    blk = [int(j) for j in np.flatnonzero(wb != 0) if j != i]
    idx = blk + [int(i)]
    n = len(idx)
    si = float(np.sign(P['U'][k][i])) or 1.0
    tr = [int(j) for j in tracers]
    st = np.array([np.sign(P['U'][k][j]) or 1.0 for j in tr])
    ctr = np.array([cb[j] for j in tr])
    Ptr = K.phi(ctr) if tr else np.zeros((nx, 0))
    al0 = np.r_[[np.sign(wb[j]) * _alpha_of(L, abs(wb[j]), coef_crit) for j in blk],
                si * _alpha_of(L, -lam0 * rho * xi0, eta_crit)]
    z0 = np.r_[al0, [cb[j] for j in idx], np.zeros(len(tr))]
    r_pre = K.filt(fstar - K.phi(cb) @ wb)
    d_pre = st * (Ptr.T @ r_pre) / nx / rho
    s0 = xi0 * rho / (si * (K.phi(np.array([cb[i]]))[:, 0] @ r_pre) / nx)

    def resid(z):
        al, c = z[:n], z[n:2 * n]
        w = coef_crit(L, al) * np.sign(al)
        return K.filt(fstar - K.phi(c) @ w), al, c, w

    def rhs(s, z):
        r, al, c, w = resid(z)
        g = K.phi(c).T @ r / nx
        h = K.dphi(c).T @ r / nx
        x = np.minimum(np.abs(al) / L.A, 350.0)
        return np.r_[lam0 * 2 * L.G * L.c * np.tanh(x) ** 2 * g, lam0 * L.G * L.m * w * h, -st * (Ptr.T @ r) / nx / rho]
    z, s_, chunk, ok = z0, 0.0, s0 + float(s_end), False
    for _ in range(max_doublings + 1):
        sol = solve_ivp(rhs, (s_, s_ + chunk), z, method='LSODA', rtol=1e-9, atol=1e-12)
        z, s_ = sol.y[:, -1], s_ + chunk
        if np.all(np.abs(rhs(s_, z)[:2 * n]) <= tol):
            ok = True
            break
        chunk = s_
    r_end, al, c, w = resid(z)
    d_post = st * (Ptr.T @ r_end) / nx / rho
    xi_out = -(d_pre * min(s_, s0) + d_post * max(s_ - s0, 0.0))
    pull = {j: float(xi_out[q] - z[2 * n + q]) for q, j in enumerate(tr)}
    return dict(pull=pull, stationary=ok, s=s_, w_end=dict(zip(idx, w.tolist())))


def recovery_threshold(L, K, fstar, c0, R, th_abs, n_creep=4, eps_lo=1e-14, eps_hi=1e-3, n_bisect=14, measured_pull=True):
    """The mass below which the loop forms exactly the sliding path's kernels, predicted from the path alone (Proposition D‴.51(iii)):
    (a) at each entry i, every kernel j that never enters is pulled toward its wall by the entrant's overshoot — it is formed when its
    margin μ_j(τ_i) < ε^{1/2}x_θ(A_i)κ_ij, x_θ(A) = ln(|A|/θ − 1)/κ (κ = πρ(2c)^{1/2} the inner layer's Gibbs rate, θ the formation
    threshold th_abs), κ_ij = s_j⟨Δr_iφ_j⟩/(s_i⟨Δr_iφ_i⟩) (Δr_i the residual the entry removes, φ_i the entrant after its relaxation) —
    the near-tie law's estimate, right for the entrant's near-tie partner, which slides onto the entrant's feature and shares it. The pull
    on a kernel that stays where it is is measured on the entry's ε-free inner problem (`entry_pull`, Δ_ij): x_ij = Δ_ij + ξ_θ
    (𝒲(ξ_θ) = θ), raised to x_θκ_ij when j is pulled (Δ_ij > 0) and receives at most the entrant's own share of the residual change
    (0 < κ_ij ≤ 1: a near-tie partner); j can form only if x_ij > 0. Where an entry rearranges the block — an active kernel sent on to
    another feature, the residual changing far from the entrant, ⟨Δr_iφ_i⟩ small — κ_ij is large and the near-tie estimate would
    overstate the pull many times. ε_entry = min (μ_j/x_ij)²; (b) the n_creep never-entering kernels closest to their walls, integrated with their creep along the
    path (`creep_reaches`): ε_creep the smallest ε (bisection in log ε) at which one reaches the inner layer. Returns dict(eps_entry,
    entry=(i, j, μ, κ_ij, x_θ, Δ_ij), eps_creep, creep={j: ε_j}, eps_c = min). measured_pull=False: the near-tie estimate alone."""
    P = path_duals(K, fstar, R)
    taus, res, cs, U = P['taus'], P['res'], P['centres'], P['U']
    d = len(c0)
    nx = K.xs.size
    never = [j for j in range(d) if j not in P['ever']]
    kap = np.pi * np.sqrt(2 * L.c)
    xi_th = float(inner_xi(L, th_abs))
    eps_entry, arg = np.inf, None
    for t, i, kind, c in R['events']:
        if kind != 'enter':
            continue
        k = int(np.argmin(np.abs(np.array(taus) - t)))
        dr = K.filt(res[k - 1] - res[k])
        gi = K.phi(np.array([cs[k][i]]))[:, 0] @ dr / nx
        Ai = abs(float(R['w'][i]))                                   # the entrant's coefficient at the path's end (0 if it left)
        if gi == 0 or Ai <= th_abs:
            continue
        xth = np.log(Ai / th_abs - 1) / kap if Ai > 2 * th_abs else 0.0
        pull = entry_pull(L, K, fstar, R, P, k, i, never)['pull'] if (measured_pull and never) else {}
        for j in never:
            sj = np.sign(U[k][j]) if U[k][j] != 0 else 1.0
            kij = sj * (K.phi(np.array([c0[j]]))[:, 0] @ dr / nx) / abs(gi)
            if j in pull:
                x_ij = pull[j] + xi_th
                if pull[j] > 0 and 0 < kij <= 1:
                    x_ij = max(x_ij, xth * kij)
            else:
                x_ij = xth * kij if (kij > 0 and xth > 0) else -np.inf
            if x_ij > 0:
                e_ = ((1 - abs(U[k][j])) / x_ij) ** 2
                if e_ < eps_entry:
                    eps_entry, arg = float(e_), (int(i), int(j), float(1 - abs(U[k][j])), float(kij), float(xth),
                                                 float(pull.get(j, np.nan)))
    closest = sorted(never, key=lambda j: min(1 - abs(U[k][j]) for k in range(len(taus))))[:n_creep]
    creep = {}
    for j in closest:
        if not creep_reaches(L, K, P, j, eps_hi, c0[j]):
            creep[j] = np.inf
            continue
        lo, hi = np.log(eps_lo), np.log(eps_hi)
        if creep_reaches(L, K, P, j, eps_lo, c0[j]):
            creep[j] = eps_lo
            continue
        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            if creep_reaches(L, K, P, j, np.exp(mid), c0[j]):
                hi = mid
            else:
                lo = mid
        creep[j] = float(np.exp(hi))
    eps_creep = min(creep.values()) if creep else np.inf
    return dict(eps_entry=eps_entry, entry=arg, eps_creep=eps_creep, creep=creep, eps_c=min(eps_entry, eps_creep))
