"""
Explicit strict Lyapunov function for the skew-symmetric mixed error model (PARADIGM.md Lemma E.0)
===================================================================================================

The linearisation of the closed loop at an ideal configuration (Theorem E, proof step 3) is

    ė = −k e − Ψ(t) ζ,      ζ̇ = Ψ(t)ᵀ e,        e ∈ ℝⁿ,  ζ ∈ ℝ^q,  Ψ(t) ∈ ℝ^{n×q},            (LTV)

with Ψ bounded and C¹ with bounded derivative, and persistently exciting:

    ‖Ψ(t)‖ ≤ Ψ_M,   ‖Ψ̇(t)‖ ≤ Ψ_D,   ∫_t^{t+T} Ψ(r)ᵀΨ(r) dr ⪰ α I_q   for all t ≥ 0.                 (PE)

Version 3 of PARADIGM.md proved uniform exponential stability of (LTV) by citing Morgan–Narendra (1977) for the
qualitative statement; the rate and the overshoot were existence constants.  Lemma E.0 replaces that citation by
a self-contained strict Lyapunov function with explicit constants.  With M(r) := Ψ(r)ᵀΨ(r),

    Q(t) := ∫_t^{t+T} (t + T − r) M(r) dr,         0 ⪯ Q(t) ⪯ Q_M I,   Q_M := T²Ψ_M²/2,
    Q̇(t)  = ∫_t^{t+T} M(r) dr − T M(t),

    V(t, e, ζ) := ½‖e‖² + ½‖ζ‖² − ε ζᵀQ(t)ζ + εT eᵀΨ(t)ζ.

Along (LTV),

    d/dt[−ζᵀQζ]  = −ζᵀ(∫_t^{t+T}M)ζ + T‖Ψζ‖² − 2ζᵀQΨᵀe  ≤  −α‖ζ‖² + T‖Ψζ‖² − 2ζᵀQΨᵀe,
    d/dt[eᵀΨζ]   = −k eᵀΨζ − ‖Ψζ‖² + eᵀΨ̇ζ + ‖Ψᵀe‖²,

so the ‖Ψζ‖² terms cancel in V̇ and

    V̇ ≤ −(k − εTΨ_M²)‖e‖² − εα‖ζ‖² + ε c₁ ‖e‖‖ζ‖,      c₁ := T²Ψ_M³ + T k Ψ_M + T Ψ_D.

Young's inequality and the closed-form choice

    ε := min{ 1 / (4(Q_M + TΨ_M/2)),   k / (4TΨ_M²),   kα / c₁² }

give  ¼‖w‖² ≤ V ≤ ¾‖w‖²  and  V̇ ≤ −γ‖w‖²  with  γ := min{k/4, εα/2},  w := (e, ζ).  Hence

    ‖w(t)‖ ≤ √3 · exp(−(2γ/3)(t − t₀)) · ‖w(t₀)‖.

Every ε with ε(Q_M + TΨ_M/2) < ½ gives a valid Lyapunov function (see strict_lyapunov_constants / optimise_eps).

Every constant is a closed-form function of (k, α, T, Ψ_M, Ψ_D).  The rate vanishes as Ψ_D → ∞, as it must
(Ψ = sin t² is bounded and PE but not uniformly exponentially stabilising; PARADIGM.md §6 step 3), and it scales
like ω⁻² under the time rescaling Ψ(ωt), which is the scaling of the exact Floquet rates in Remark E.1.

The construction is of the Mazenc–Malisoff type (strict Lyapunov functions with a backward–forward integral of the
PE Gramian); the contribution here is the explicit, self-contained version for the mixed class (LTV).
"""
import numpy as np


# ---------------------------------------------------------------- constants of Lemma E.0
def strict_lyapunov_constants(k, alpha, T, Psi_M, Psi_D, eps=None):
    """
    Constants of Lemma E.0 for the data (k, α, T, Ψ_M, Ψ_D).

    With Q_M := T²Ψ_M²/2 and c₁ := T²Ψ_M³ + TkΨ_M + TΨ_D, every ε with ε(Q_M + TΨ_M/2) < ½ is admissible, and

        V̇ ≤ −a‖e‖² − b‖ζ‖² + c‖e‖‖ζ‖,     a := k − εTΨ_M²,  b := εα,  c := εc₁,
        p_lo‖w‖² ≤ V ≤ p_hi‖w‖²,           p_lo/hi := ½ ∓ ε(Q_M + TΨ_M/2),

    so V̇ ≤ −γ‖w‖² with γ := λ_min[[a, −c/2], [−c/2, b]] (when positive), and

        ‖w(t)‖ ≤ √(p_hi/p_lo) · exp(−(γ/(2p_hi))(t − t₀)) · ‖w(t₀)‖.

    The closed-form choice of the lemma is ε := min{1/(4(Q_M + TΨ_M/2)), k/(4TΨ_M²), kα/c₁²}, for which
    p_lo ≥ ¼, p_hi ≤ ¾, γ ≥ γ_cf := min{k/4, εα/2}, the overshoot is at most √3 and the rate at least (2/3)γ_cf.
    Pass eps to evaluate another admissible ε (see optimise_eps).
    """
    k, alpha, T, Psi_M, Psi_D = map(float, (k, alpha, T, Psi_M, Psi_D))
    if min(k, alpha, T, Psi_M) <= 0 or Psi_D < 0:
        raise ValueError("need k, α, T, Ψ_M > 0 and Ψ_D ≥ 0")
    Q_M = T ** 2 * Psi_M ** 2 / 2.0
    c1 = T ** 2 * Psi_M ** 3 + T * k * Psi_M + T * Psi_D
    eps_V = 1.0 / (4.0 * (Q_M + T * Psi_M / 2.0))
    closed_form = eps is None
    if closed_form:
        eps = min(eps_V, k / (4.0 * T * Psi_M ** 2), k * alpha / c1 ** 2)
    eps = float(eps)
    if not (0 < eps * (Q_M + T * Psi_M / 2.0) < 0.5):
        raise ValueError("ε must satisfy 0 < ε(Q_M + TΨ_M/2) < 1/2 for V to be positive definite")
    a, b, c = k - eps * T * Psi_M ** 2, eps * alpha, eps * c1
    lam_max = 0.5 * (a + b) + np.sqrt(0.25 * (a - b) ** 2 + 0.25 * c ** 2)
    gamma = (a * b - 0.25 * c ** 2) / lam_max            # λ_min of the 2×2 form as det/λ_max: no cancellation when b, c ≪ a (round-5 review, E0-1)
    p_lo = 0.5 - eps * (Q_M + T * Psi_M / 2.0)
    p_hi = 0.5 + eps * (Q_M + T * Psi_M / 2.0)
    gamma_cf = min(k / 4.0, eps * alpha / 2.0) if closed_form else np.nan
    return dict(eps=eps, gamma=float(gamma), gamma_cf=float(gamma_cf), Q_M=Q_M, c1=c1, a=a, b=b, c=c,
                p_lo=p_lo, p_hi=p_hi,
                rate=float(gamma / (2.0 * p_hi)),                 # ‖w(t)‖ ≤ overshoot·e^{−rate (t−t₀)}‖w(t₀)‖
                overshoot=float(np.sqrt(p_hi / p_lo)),
                grad_V=1.0 + eps * (2.0 * Q_M + T * Psi_M))        # ‖∂V/∂w‖ ≤ grad_V · ‖w‖  (‖H − I‖ ≤ 2εQ_M + εTΨ_M)


def optimise_eps(k, alpha, T, Psi_M, Psi_D, n_grid=400):
    """The admissible ε maximising the guaranteed rate γ/(2p_hi) over a log grid that always contains the closed-form ε
    (every admissible ε gives a valid Lyapunov function, so the maximum over the grid is a valid bound for the given
    data; the data (α, Ψ_M, Ψ_D) themselves are what the caller supplies).  Never returns None: the closed-form constants
    are returned when no grid point beats them."""
    Q_M = T ** 2 * Psi_M ** 2 / 2.0
    eps_V = 0.5 / (Q_M + T * Psi_M / 2.0)
    best = strict_lyapunov_constants(k, alpha, T, Psi_M, Psi_D)                  # closed form: γ ≥ γ_cf > 0
    lo = min(eps_V * 1e-9, 1e-3 * best['eps'])
    for eps in np.geomspace(lo, eps_V * (1 - 1e-6), n_grid):
        c = strict_lyapunov_constants(k, alpha, T, Psi_M, Psi_D, eps=eps)
        if c['gamma'] > 0 and c['rate'] > best['rate']:
            best = c
    return best


# ---------------------------------------------------------------- Gramians by Gauss–Legendre quadrature
def _gl_nodes(a, b, n_nodes):
    g, w = np.polynomial.legendre.leggauss(n_nodes)
    return 0.5 * (b - a) * g + 0.5 * (b + a), 0.5 * (b - a) * w


def gramian(Psi_fun, t, T, n_nodes=200):
    """∫_t^{t+T} Ψ(r)ᵀΨ(r) dr."""
    r, w = _gl_nodes(t, t + T, n_nodes)
    out = None
    for ri, wi in zip(r, w):
        P = Psi_fun(ri)
        M = P.T @ P
        out = wi * M if out is None else out + wi * M
    return out


def Q_matrix(Psi_fun, t, T, n_nodes=200):
    """Q(t) = ∫_t^{t+T} (t + T − r) Ψ(r)ᵀΨ(r) dr  (the weight makes Q̇ = ∫_t^{t+T}M − T M(t))."""
    r, w = _gl_nodes(t, t + T, n_nodes)
    out = None
    for ri, wi in zip(r, w):
        P = Psi_fun(ri)
        M = (t + T - ri) * (P.T @ P)
        out = wi * M if out is None else out + wi * M
    return out


def pe_level(Psi_fun, T, t_grid, n_nodes=200, check=True, rtol=1e-6):
    """α(T) := min over t in t_grid of λ_min(∫_t^{t+T} ΨᵀΨ)  (the PE level on the grid).  With check=True the quadrature
    is repeated with 2·n_nodes at the minimising t and a ValueError is raised if the two values differ by more than
    rtol: the Gauss–Legendre rule must resolve every oscillation of Ψ in the window (round-5 review, E0-3)."""
    vals = [(np.linalg.eigvalsh(gramian(Psi_fun, t, T, n_nodes))[0], t) for t in t_grid]
    a, t_min = min(vals, key=lambda v: v[0])
    if check:
        a2 = np.linalg.eigvalsh(gramian(Psi_fun, t_min, T, 2 * n_nodes))[0]
        if abs(a2 - a) > rtol * max(abs(a), abs(a2), 1e-300):
            raise ValueError(f"pe_level: quadrature not converged ({a:.6e} with {n_nodes} nodes, {a2:.6e} with {2 * n_nodes}); increase n_nodes")
    return float(a)


def sup_norms(Psi_fun, Psidot_fun, t_grid):
    """(Ψ_M, Ψ_D) as maxima of the operator norms over the grid."""
    Psi_M = max(np.linalg.norm(Psi_fun(t), 2) for t in t_grid)
    Psi_D = max(np.linalg.norm(Psidot_fun(t), 2) for t in t_grid)
    return float(Psi_M), float(Psi_D)


# ---------------------------------------------------------------- the strict Lyapunov function
class StrictLyapunov:
    """V(t, e, ζ) of Lemma E.0 for given Ψ(·), Ψ̇(·) and constants; evaluates V, V̇ along (LTV), and the bounds."""

    def __init__(self, Psi_fun, Psidot_fun, k, T, consts, n_nodes=200):
        self.Psi, self.Psidot, self.k, self.T, self.c, self.n_nodes = Psi_fun, Psidot_fun, float(k), float(T), consts, n_nodes

    def Q(self, t):
        return Q_matrix(self.Psi, t, self.T, self.n_nodes)

    def Qdot(self, t):
        P = self.Psi(t)
        return gramian(self.Psi, t, self.T, self.n_nodes) - self.T * (P.T @ P)

    def V(self, t, e, zeta):
        eps, T = self.c['eps'], self.T
        P = self.Psi(t)
        return 0.5 * e @ e + 0.5 * zeta @ zeta - eps * zeta @ self.Q(t) @ zeta + eps * T * e @ P @ zeta

    def field(self, t, e, zeta):
        P = self.Psi(t)
        return -self.k * e - P @ zeta, P.T @ e

    def Vdot(self, t, e, zeta):
        """Exact time derivative of V along (LTV) by the chain rule (Q̇ analytic)."""
        eps, T = self.c['eps'], self.T
        P, Pd = self.Psi(t), self.Psidot(t)
        ed, zd = self.field(t, e, zeta)
        dV = e @ ed + zeta @ zd
        dV += -eps * (zeta @ self.Qdot(t) @ zeta + 2.0 * zeta @ self.Q(t) @ zd)
        dV += eps * T * (ed @ P @ zeta + e @ Pd @ zeta + e @ P @ zd)
        return float(dV)

    def check_point(self, t, e, zeta):
        """Returns (V, V̇, γ‖w‖², p_lo‖w‖², p_hi‖w‖²) at one state: the lemma asserts V̇ ≤ −γ‖w‖² and p_lo‖w‖² ≤ V ≤ p_hi‖w‖²."""
        w2 = float(e @ e + zeta @ zeta)
        return self.V(t, e, zeta), self.Vdot(t, e, zeta), self.c['gamma'] * w2, self.c['p_lo'] * w2, self.c['p_hi'] * w2


# ---------------------------------------------------------------- exact rates for periodic Ψ (reference values)
def ltv_matrix(Psi, k):
    n, q = Psi.shape
    return np.block([[-k * np.eye(n), -Psi], [Psi.T, np.zeros((q, q))]])


def floquet_rate(Psi_fun, k, period, n_steps=4000):
    """Exact exponential rate of (LTV) for a `period`-periodic Ψ: −max_i log|λ_i(Φ(period))| / period, Φ the monodromy
    matrix (RK4 on the matrix equation).  This is the sharp rate the lemma's bound must stay below."""
    P0 = Psi_fun(0.0)
    d = P0.shape[0] + P0.shape[1]
    Phi = np.eye(d)
    dt = period / n_steps
    for j in range(n_steps):
        t = j * dt
        A1 = ltv_matrix(Psi_fun(t), k); A2 = ltv_matrix(Psi_fun(t + 0.5 * dt), k); A4 = ltv_matrix(Psi_fun(t + dt), k)
        k1 = A1 @ Phi; k2 = A2 @ (Phi + 0.5 * dt * k1); k3 = A2 @ (Phi + 0.5 * dt * k2); k4 = A4 @ (Phi + dt * k3)
        Phi = Phi + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    lam = np.linalg.eigvals(Phi)
    return float(-np.log(np.max(np.abs(lam))) / period)


def simulate_ltv(Psi_fun, k, w0, T_end, dt):
    """RK4 on (LTV); returns (t, ‖w(t)‖)."""
    n = len(w0)
    w = np.array(w0, float)
    ts, norms = [0.0], [np.linalg.norm(w)]
    nst = int(round(T_end / dt))
    for j in range(nst):
        t = j * dt
        A1 = ltv_matrix(Psi_fun(t), k); A2 = ltv_matrix(Psi_fun(t + 0.5 * dt), k); A4 = ltv_matrix(Psi_fun(t + dt), k)
        k1 = A1 @ w; k2 = A2 @ (w + 0.5 * dt * k1); k3 = A2 @ (w + 0.5 * dt * k2); k4 = A4 @ (w + dt * k3)
        w = w + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        ts.append((j + 1) * dt); norms.append(np.linalg.norm(w))
    return np.array(ts), np.array(norms)


# ---------------------------------------------------------------- the regressor of Theorem E
class TheoremERegressor:
    """
    Ψ_G(t) = m^{−1/2} [∇_ϑΦ(x_m(t);ϑ*_1), …, ∇_ϑΦ(x_m(t);ϑ*_m)] (I_m ⊗ G^{1/2})  ∈ ℝ^{n×mp}, and its exact time
    derivative, for Gaussian features Φ(x;ϑ) = a φ(x;c,s), φ = exp(−‖x−c‖²/2s²).  These are the w-coordinates of
    Theorem E step 3: ė = −k_e e − Ψ_G ζ, ζ̇ = Ψ_Gᵀ e.
    """

    def __init__(self, P_star, G_diag, x_ref, x_ref_dot):
        self.P = np.atleast_2d(np.asarray(P_star, float))
        self.m, self.p = self.P.shape
        self.n = (self.p - 1) // 2
        self.sqG = np.sqrt(np.asarray(G_diag, float))
        self.xr, self.xrd = x_ref, x_ref_dot

    def _blocks(self, t):
        n = self.n
        x, xd = np.asarray(self.xr(t), float), np.asarray(self.xrd(t), float)
        Psi = np.zeros((n, self.m * self.p)); Psid = np.zeros_like(Psi)
        for j in range(self.m):
            a, c, s = self.P[j, :n], self.P[j, n:2 * n], self.P[j, 2 * n]
            u = x - c
            phi = np.exp(-u @ u / (2 * s * s))
            phid = phi * (-(u @ xd) / s ** 2)
            gc = phi * u / s ** 2;                          gcd = phid * u / s ** 2 + phi * xd / s ** 2
            gs = phi * (u @ u) / s ** 3;                    gsd = phid * (u @ u) / s ** 3 + phi * 2 * (u @ xd) / s ** 3
            B = np.zeros((n, self.p)); Bd = np.zeros((n, self.p))
            B[:, :n] = phi * np.eye(n);                     Bd[:, :n] = phid * np.eye(n)
            B[:, n:2 * n] = np.outer(a, gc);                Bd[:, n:2 * n] = np.outer(a, gcd)
            B[:, 2 * n] = a * gs;                           Bd[:, 2 * n] = a * gsd
            Psi[:, j * self.p:(j + 1) * self.p] = B * self.sqG[None, :] / np.sqrt(self.m)
            Psid[:, j * self.p:(j + 1) * self.p] = Bd * self.sqG[None, :] / np.sqrt(self.m)
        return Psi, Psid

    def Psi(self, t):
        return self._blocks(t)[0]

    def Psidot(self, t):
        return self._blocks(t)[1]


def explicit_rate(Psi_fun, Psidot_fun, k, T_candidates, t_grid, n_nodes=200):
    """
    Lemma E.0's rate, maximised over the PE window T (any window length is admissible with its own α(T)).
    Returns (best rate, dict of the data and constants at the best T, table of all candidates).
    """
    Psi_M, Psi_D = sup_norms(Psi_fun, Psidot_fun, t_grid)
    rows = []
    for T in T_candidates:
        alpha = pe_level(Psi_fun, T, t_grid, n_nodes)
        if alpha <= 0:
            rows.append(dict(T=T, alpha=alpha, rate=0.0)); continue
        c = strict_lyapunov_constants(k, alpha, T, Psi_M, Psi_D)
        rows.append(dict(T=T, alpha=alpha, rate=c['rate'], eps=c['eps'], gamma=c['gamma']))
    best = max(rows, key=lambda r: r['rate'])
    best = dict(best, Psi_M=Psi_M, Psi_D=Psi_D, k=k)
    return best['rate'], best, rows


# ---------------------------------------------------------------- the nonlinear radius of Theorem E (explicit form)
def theorem_E_constants(k, alpha, T, Psi_M, Psi_D, m, G_diag, r, L_J, eps=None):
    """
    Explicit constants of Theorem E (revised).  Inputs: the LTV data (k, α, T, Ψ_M, Ψ_D) of Ψ_G, the number of
    kernels m, the gain diagonal, r := min(ρ, δ) (radius of the ball D in z = (e, ϑ̃)), and L_J := a Lipschitz
    constant of the Jacobian of the closed-loop field on D, uniform in t, so that ‖R(t,z)‖ ≤ (L_J/2)‖z‖².

    With S := diag(I, m^{−1/2}(I⊗G^{−1/2})) the similarity to w-coordinates, the remainder in w satisfies
    ‖R̃(t,w)‖ ≤ C̃_R‖w‖² with C̃_R := ‖S‖(L_J/2)‖S⁻¹‖², and D contains {‖w‖ < r/‖S⁻¹‖}.  On the ball
    ‖w‖ ≤ r_w := min{r/‖S⁻¹‖, γ/(2 g_V C̃_R)}, where ‖∂V/∂w‖ ≤ g_V‖w‖, the strict Lyapunov function gives
    V̇ ≤ −(γ/2)‖w‖², and {V ≤ p_lo r_w²} ⊂ {‖w‖ ≤ r_w} is forward invariant.  Hence for ‖w(t₀)‖ ≤ r₀ := r_w√(p_lo/p_hi),

        ‖w(t)‖ ≤ √(p_hi/p_lo) · exp(−(γ/(4p_hi))(t − t₀)) · ‖w(t₀)‖.

    Transfer to (‖e‖, W).  In w-coordinates ‖ζ‖² = (1/m)Σ‖ϑ̃_i‖²_{G⁻¹} = W_π², so ‖w‖² = 2V_π exactly.  Hence
    V_π(t₀) ≤ c₀ := r₀²/2 gives ‖w(t₀)‖ ≤ r₀, and ‖e‖ + W ≤ ‖e‖ + W_π ≤ √2‖w‖ with ‖w(t₀)‖ ≤ ‖e(t₀)‖ + W(t₀) when π
    is optimal at t₀, so C := √(2p_hi/p_lo) ≤ √6, free of m and G.  r0_z := r0_w/‖S‖ is a sufficient radius in z.
    """
    c = strict_lyapunov_constants(k, alpha, T, Psi_M, Psi_D, eps=eps)
    G = np.asarray(G_diag, float)
    lam_G, lam_min = float(G.max()), float(G.min())
    S_norm = max(1.0, (m * lam_min) ** -0.5)
    Sinv_norm = max(1.0, (m * lam_G) ** 0.5)
    C_R = S_norm * (L_J / 2.0) * Sinv_norm ** 2
    r_w = r / Sinv_norm if C_R <= 0 else min(r / Sinv_norm, c['gamma'] / (2.0 * c['grad_V'] * C_R))
    r0_w = r_w * np.sqrt(c['p_lo'] / c['p_hi'])
    r0_z = r0_w / S_norm                                   # ‖w‖ ≤ ‖S‖‖z‖, so ‖z‖ ≤ r0_z ⇒ ‖w‖ ≤ r0_w
    c0 = r0_w ** 2 / 2.0
    C = np.sqrt(2.0) * c['overshoot']
    return dict(c, r_w=r_w, r0_w=r0_w, r0_z=r0_z, c0=c0, C=C, C_R=C_R, S_norm=S_norm, Sinv_norm=Sinv_norm,
                rate_nonlinear=float(c['gamma'] / (4.0 * c['p_hi'])))


def jacobian_lipschitz_estimate(field, z_dim, r, t_grid, rng, n_samples=200, h=1e-5):
    """
    Numerical estimate of L_J on the ball ‖z‖ < r: the largest ratio ‖DF(t,z) − DF(t,z')‖/‖z − z'‖ over random pairs
    (central differences for DF).  This is an estimate from samples, not a certified bound; Theorem E's statement
    takes L_J as the true Lipschitz constant, and a closed-form bound would follow from the third-derivative
    suprema of the Gaussian (not computed here).
    """
    def DF(t, z):
        J = np.zeros((z_dim, z_dim))
        for j in range(z_dim):
            ej = np.zeros(z_dim); ej[j] = h
            J[:, j] = (field(z + ej, t) - field(z - ej, t)) / (2 * h)
        return J
    worst = 0.0
    for _ in range(n_samples):
        t = rng.choice(t_grid)
        z = rng.normal(size=z_dim); z *= rng.uniform(0, r) / np.linalg.norm(z)
        zp = rng.normal(size=z_dim); zp *= rng.uniform(0, r) / np.linalg.norm(zp)
        worst = max(worst, np.linalg.norm(DF(t, z) - DF(t, zp), 2) / np.linalg.norm(z - zp))
    return float(worst)
