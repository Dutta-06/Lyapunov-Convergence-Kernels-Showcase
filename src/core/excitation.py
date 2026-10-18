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

Lemma E.0′ (SlowRegimeLyapunov, at the end of this module) is the sharp companion for the slow-adaptation regime of a
T-periodic Ψ: a second-order slaving of e to its slow manifold and a second-order averaging correction give a strict
Lyapunov function whose rate is the averaged rate λ_min⟨ΨᵀΨ⟩/k up to a factor 1 − δ, for k above a threshold set by the
peak factors of Ψ, Ψ̇, Ψ̈; the inequality V̇ ≤ −γV is certified pointwise by a generalised eigenvalue problem.
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
    matrix (RK4 on the matrix equation).  This is the sharp rate the lemma's bound must stay below.  The step count is
    raised if needed so that k·dt ≤ 1/4 (RK4 is unstable on the fast channel beyond k·dt ≈ 2.8)."""
    P0 = Psi_fun(0.0)
    d = P0.shape[0] + P0.shape[1]
    Phi = np.eye(d)
    n_steps = max(int(n_steps), int(np.ceil(4.0 * k * period)))
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


# ---------------------------------------------------------------- Lemma E.0′: the slow-adaptation regime
def zero_mean_antiderivative(F, T):
    """
    F: (N, …) samples of a smooth T-periodic function on the uniform grid t_j = jT/N.  Returns the samples of the
    zero-mean antiderivative of F − ⟨F⟩ on the same grid (division by iω in Fourier space: spectral accuracy for smooth
    periodic data).  A zero-mean antiderivative g of a zero-mean f obeys ‖g‖_∞ ≤ (T/4)‖f‖_∞.
    """
    F = np.asarray(F, float)
    Fh = np.fft.fft(F - F.mean(axis=0), axis=0)
    om = 2 * np.pi * np.fft.fftfreq(F.shape[0], d=T / F.shape[0])
    shape = (-1,) + (1,) * (F.ndim - 1)
    om = om.reshape(shape)
    with np.errstate(divide='ignore', invalid='ignore'):
        Gh = np.where(om != 0, Fh / (1j * om), 0.0)
    return np.real(np.fft.ifft(Gh, axis=0))


class SlowRegimeLyapunov:
    """
    Lemma E.0′ (PARADIGM.md §6): for a T-periodic regressor Ψ with M := ΨᵀΨ and M̄ := ⟨M⟩ ≻ 0, the strict Lyapunov
    function

        V := ½ ζᵀ(I + εP₁(t) + ε²P₂(t))ζ + (β/2)‖w‖²,   w := e − e₂(ζ,t),   e₂ := −εΨζ + ε²Ψ̇ζ − ε³ΨMζ,   ε := 1/k,

    with P₁ the zero-mean antiderivative of 2(M − M̄) and P₂ := −(M − M̄) + S, S the zero-mean antiderivative of
    F := P₁M + MP₁ (whose mean N₂ = ½⟨P₁M + MP₁⟩ vanishes identically), satisfies along (LTV)

        V̇ ≤ −(1/k)(1 − δ − r_ζ k) ζᵀM̄ζ − (βk/4)‖w‖² ≤ −γ V,    γ = min{2c_ζλ_min(M̄)/(k(1 + εp₁ + ε²p₂)), k/2},

    for k above a threshold fixed by the *peak factors* ρ, ρ_D, ρ_DD of Ψ, Ψ̇, Ψ̈ against M̄ (conditions (A)–(D) of the
    lemma).  The rate γ/2 ≈ (1 − δ)λ_min(M̄)/k is the averaged rate of Remark E.1: the bound pays for the ill-conditioning
    of M̄ once, where Lemma E.0 pays twice.  Everything is computed from samples of Ψ, Ψ̇ on a uniform grid over one
    period (Ψ̈ by central differences of Ψ̇).
    """

    def __init__(self, Psi_fun, Psidot_fun, T, N=2048, h=1e-4):
        self.T, self.N = float(T), int(N)
        self.ts = np.arange(self.N) * self.T / self.N
        self.Ps = np.array([Psi_fun(t) for t in self.ts])
        self.Pd = np.array([Psidot_fun(t) for t in self.ts])
        self.Pdd = np.array([(Psidot_fun(t + h) - Psidot_fun(t - h)) / (2 * h) for t in self.ts])
        self.n, self.q = self.Ps.shape[1:]
        self.Ms = np.einsum('ikj,ikl->ijl', self.Ps, self.Ps)
        self.Mdot = np.einsum('ikj,ikl->ijl', self.Pd, self.Ps) + np.einsum('ikj,ikl->ijl', self.Ps, self.Pd)
        self.Mbar = self.Ms.mean(0)
        lam, U = np.linalg.eigh(self.Mbar)
        if lam[0] <= 0:
            raise ValueError("M̄ is not positive definite: the regressor is not persistently exciting over the period")
        self.lam_min, self.lam_max = float(lam[0]), float(lam[-1])
        Mh = U @ np.diag(lam ** -0.5) @ U.T
        peak = lambda Xs: float(max(np.linalg.eigvalsh(Mh @ x.T @ x @ Mh)[-1] for x in Xs))
        self.rho, self.rho_D, self.rho_DD = peak(self.Ps), peak(self.Pd), peak(self.Pdd)
        self.Psi_M = float(max(np.linalg.norm(p, 2) for p in self.Ps))
        self.Psi_D = float(max(np.linalg.norm(p, 2) for p in self.Pd))
        self.Psi_DD = float(max(np.linalg.norm(p, 2) for p in self.Pdd))
        self.P1 = zero_mean_antiderivative(2 * (self.Ms - self.Mbar), self.T)
        PM = np.einsum('ijk,ikl->ijl', self.P1, self.Ms)
        self.F = PM + np.swapaxes(PM, 1, 2)
        self.N2 = 0.5 * self.F.mean(0)                     # = 0 exactly; kept as a numerical check
        self.S = zero_mean_antiderivative(self.F, self.T)
        self.P2 = -(self.Ms - self.Mbar) + self.S
        self.P2dot = self.F - self.Mdot
        self.p1_meas = float(max(np.linalg.norm(p, 2) for p in self.P1))
        self.p2_meas = float(max(np.linalg.norm(p, 2) for p in self.P2))

    # ---- the constants of the lemma
    def constants(self, k, delta=0.1):
        k = float(k); e = 1.0 / k; T = self.T
        PM, PD = self.Psi_M, self.Psi_D
        sr, srD, srDD = np.sqrt(self.rho), np.sqrt(self.rho_D), np.sqrt(self.rho_DD)
        pi1 = 0.5 * T * PM * (sr + 1)                                   # ‖P₁v‖ ≤ π₁‖v‖_M̄
        pi2 = PM * (sr + 1) + 0.25 * T ** 2 * PM ** 3 * (2 * sr + 1)   # ‖P₂v‖ ≤ π₂‖v‖_M̄
        p1 = 0.5 * T * PM ** 2                                          # ‖P₁‖ ≤ p₁
        p2 = PM ** 2 + 0.5 * T ** 2 * PM ** 4                           # ‖P₂‖ ≤ p₂
        r_z = (e ** 3 * (pi2 * PM * sr + pi1 * PM * srD + PM ** 2 * self.rho)
               + e ** 4 * (pi2 * PM * srD + pi1 * PM ** 3 * sr) + e ** 5 * pi2 * PM ** 3 * sr)
        c1 = sr + e * PM * pi1 + e ** 2 * PM * pi2
        c2 = (e ** 2 * srDD + e ** 3 * (3 * PD * PM * sr + 2 * PM ** 2 * srD) + e ** 4 * (2 * PM ** 4 * sr + PD * PM * srD)
              + e ** 5 * (PM ** 4 * srD + PD * PM ** 3 * sr) + e ** 6 * PM ** 5 * sr)
        d_w = e * PM ** 2 + e ** 2 * PD * PM + e ** 3 * PM ** 4
        beta = 2 * c1 ** 2 / delta
        c_z = (1 - delta) - r_z / e
        hi, lo = 1 + e * p1 + e ** 2 * p2, 1 - e * p1 - e ** 2 * p2
        conds = dict(A=bool(c1 * c2 <= delta / 2), B=bool(d_w <= k / 4), C=bool(c_z > 0), D=bool(lo > 0))
        ok = all(conds.values())
        gamma = min(2 * e * c_z * self.lam_min / hi, k / 2) if ok else float('nan')
        L_M = e * PM + e ** 2 * PD + e ** 3 * PM ** 3                     # ‖e₂‖ ≤ L_M ‖ζ‖
        over = (1 + L_M) ** 2 * np.sqrt(max(hi, beta) / min(lo, beta)) if ok else float('nan')
        return dict(k=k, delta=delta, beta=beta, c1=c1, c2=c2, c1c2=c1 * c2, d_w=d_w, r_z=r_z, c_z=c_z, p1=p1, p2=p2,
                    pi1=pi1, pi2=pi2, hi=hi, lo=lo, conds=conds, ok=ok, gamma=gamma, rate=gamma / 2 if ok else float('nan'),
                    overshoot=over, averaged_rate=self.lam_min / k)

    def threshold(self, delta=0.1, k_max=1e6):
        """The smallest k on a fine log grid at which conditions (A)–(D) hold (nan if none below k_max)."""
        for k in np.geomspace(1.0, k_max, 2000):
            if self.constants(k, delta)['ok']:
                return float(k)
        return float('nan')

    # ---- the certificate: sup over states of V̇/V on the time grid
    def quadratic_forms(self, i, k, beta):
        """(H, N_sym) at grid index i: V = ½xᵀHx and V̇ = xᵀN_sym x along (LTV) for x = (e, ζ)."""
        e = 1.0 / k; n, q = self.n, self.q
        Ps, Pd, Pdd, M, Md = self.Ps[i], self.Pd[i], self.Pdd[i], self.Ms[i], self.Mdot[i]
        L = -e * Ps + e ** 2 * Pd - e ** 3 * Ps @ M
        Ldot = -e * Pd + e ** 2 * Pdd - e ** 3 * (Pd @ M + Ps @ Md)
        Hz = np.eye(q) + e * self.P1[i] + e ** 2 * self.P2[i]
        Hzdot = 2 * e * (M - self.Mbar) + e ** 2 * self.P2dot[i]
        H = np.block([[beta * np.eye(n), -beta * L], [-beta * L.T, Hz + beta * L.T @ L]])
        Hdot = np.block([[np.zeros((n, n)), -beta * Ldot], [-beta * Ldot.T, Hzdot + beta * (Ldot.T @ L + L.T @ Ldot)]])
        A = ltv_matrix(Ps, k)
        return H, 0.5 * (H @ A + A.T @ H + Hdot)

    def certificate(self, k, delta=0.1, n_t=256):
        """sup over the grid of t and over all states of V̇/V (a generalised eigenvalue problem at each t); the lemma
        asserts that it is ≤ −γ.  Returns (sup V̇/V, −γ)."""
        c = self.constants(k, delta)
        worst = -np.inf
        for i in np.linspace(0, self.N - 1, n_t).astype(int):
            H, Nsym = self.quadratic_forms(i, k, c['beta'])
            lamH, UH = np.linalg.eigh(H)
            Hm = UH @ np.diag(lamH ** -0.5) @ UH.T
            worst = max(worst, 2 * float(np.linalg.eigvalsh(Hm @ Nsym @ Hm)[-1]))
        return worst, -c['gamma']


def slow_regime_theorem_E_constants(slow, k, delta, m, G_diag, r, L_J):
    """
    Corollary E.0′.2 (Theorem E in the slow regime): the constants of Theorem E's steps 4–5 with Lemma E.0′'s V in place
    of Lemma E.0's.  `slow` is a SlowRegimeLyapunov built from Ψ_G, Ψ̇_G; (k, δ) must satisfy its conditions.  With
    2p_lo ≤ inf_t λ_min(𝐇(t)) and 2p_hi ≥ sup_t λ_max(𝐇(t)) (explicit bounds below; the grid values are also returned),
    V := ½wᵀ𝐇w obeys V̇ ≤ −γ′V along (LTV), and along the nonlinear loop V̇ ≤ −(γ′/2)V on the forward-invariant set
    {V ≤ c₀}, c₀ := min{[γ′p_lo/(4√p_hi C̃_R)]², p_lo r²/‖S⁻¹‖²}; hence for ‖w(t₀)‖ ≤ r₀ := (c₀/p_hi)^{1/2}

        ‖w(t)‖ ≤ (p_hi/p_lo)^{1/2} · exp(−(γ′/4)(t − t₀)) · ‖w(t₀)‖,

    and Theorem E's conclusions hold with λ := γ′/4 and C := (2p_hi/p_lo)^{1/2}.
    """
    c = slow.constants(k, delta)
    if not c['ok']:
        raise ValueError(f"the conditions of Lemma E.0′ fail at k = {k}, δ = {delta}: {c['conds']}")
    e = 1.0 / k
    L_M = e * slow.Psi_M + e ** 2 * slow.Psi_D + e ** 3 * slow.Psi_M ** 3
    beta, lo, hi = c['beta'], c['lo'], c['hi']
    p_lo = 0.5 * min(beta / 2.0, lo - beta * L_M ** 2)
    p_hi = 0.5 * max(beta * (1 + L_M), hi + beta * L_M * (1 + L_M))
    if p_lo <= 0:
        raise ValueError("the explicit lower bound on the quadratic form is not positive (βL_M² ≥ 1 − εp₁ − ε²p₂)")
    lam_lo, lam_hi = np.inf, 0.0                                    # grid values of λ_min(𝐇), λ_max(𝐇) for information
    for i in np.linspace(0, slow.N - 1, 128).astype(int):
        H, _ = slow.quadratic_forms(i, k, beta)
        ev = np.linalg.eigvalsh(H); lam_lo, lam_hi = min(lam_lo, ev[0]), max(lam_hi, ev[-1])
    G = np.asarray(G_diag, float)
    S_norm = max(1.0, (m * float(G.min())) ** -0.5)
    Sinv_norm = max(1.0, (m * float(G.max())) ** 0.5)
    C_R = S_norm * (L_J / 2.0) * Sinv_norm ** 2
    gamma = c['gamma']
    v0_sq = (gamma * p_lo / (4.0 * np.sqrt(p_hi) * C_R)) ** 2 if C_R > 0 else np.inf
    c0 = min(v0_sq, p_lo * r ** 2 / Sinv_norm ** 2)
    r0_w = np.sqrt(c0 / p_hi)
    return dict(c, p_lo=p_lo, p_hi=p_hi, p_lo_grid=0.5 * lam_lo, p_hi_grid=0.5 * lam_hi, L_M=L_M, C_R=C_R, S_norm=S_norm,
                Sinv_norm=Sinv_norm, c0=c0, r0_w=r0_w, r0_z=r0_w / S_norm, c0_VOT=r0_w ** 2 / 2.0,
                rate_nonlinear=gamma / 4.0, C=np.sqrt(2.0 * p_hi / p_lo))

