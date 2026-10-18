"""
Theorem D‴ (PARADIGM.md §5d): uniform-in-time propagation of chaos for the averaged slow loop in the learning regime
=======================================================================================================================

The averaged slow dynamics of the finite-m closed loop (Theorem J′'s leading order L̄₀ with the dither switched off; slow
time s = t/k_e, leak σ₀ = k_eσ) is the G-preconditioned gradient flow of J_m(Θ) = m L̄(μ_Θ) + (σ₀/2)Σ‖ϑ_i − ϑ⁰‖²_{G⁻¹}:

    dϑ_i/ds = −G∇_{ϑ_i}J_m(Θ) = G⟨∇_ϑ(r_Θ Φ̃(·;ϑ_i))⟩_ρ − σ₀(ϑ_i − ϑ⁰),    r_Θ := f − f̂_Θ,                       (𝔄_m)

with ρ the occupation measure of the reference (here: the time grid of one period).  Its mean-field limit is the Wasserstein
gradient flow (𝔄) of §6c at τ = 0.  Under the synchronous coupling with i.i.d. reference particles ϑ̄_i(s) = X_s(ϑ_i(0))
(X_s the characteristic flow of the mean-field equation), the mean-square distance D := (1/m)Σ‖ϑ_i − ϑ̄_i‖²_{G⁻¹} obeys

    ½ dD/ds = T_A + T_B + T_C,
    T_A := −(1/m)Σ⟨δ_i, ∇𝒱[μ_Θ̄](ϑ_i) − ∇𝒱[μ_Θ̄](ϑ̄_i)⟩ ≤ −(σ₀ − κ_R)D         (self term; residual-weighted curvature),
    T_B := −(1/m)Σ⟨δ_i, ∇V_L[μ_Θ](ϑ_i) − ∇V_L[μ_Θ̄](ϑ_i)⟩ = −‖h‖² + O(D²)    (interaction; flat convexity makes it ≤ 0 at
                                                                            linear order; h := f̂_Θ − f̂_Θ̄),
    T_C := −(1/m)Σ⟨δ_i, ∇V_L[μ_Θ̄](ϑ̄_i) − ∇V_L[μ_s](ϑ̄_i)⟩ ≤ L_ϑ^G √D · Y    (the reference's sampling fluctuation Y).

κ_R[μ](A) := sup_{ϑ∈A} λ_max(G^{1/2}R[μ](ϑ)G^{1/2}) with R[μ](ϑ) := ⟨r_μ ∇²_ϑΦ̃(·;ϑ)⟩_ρ (Proposition K's R₂ at a collapsed
configuration).  This module implements (𝔄_m) for n = 1, the three terms, κ_R, the whitened Lipschitz constants, the
theorem's constant C₁ = 4e^{2Λ₁}(L_ϑ^G)²K_Φ²/κ², and the coupled runs of `experiments/exp_chaos_learning.py`.
"""
import numpy as np


def _sig(z):
    return 1.0 / (1.0 + np.exp(-z))


class AveragedSlowLoop:
    """(𝔄_m) for a saturated network with n = 1: plant f (vectorised on (N,) arrays), reference x_ref(t) → (1,), period,
    diagonal gains G (p,), leak σ₀ (slow-time normalisation σ₀ = k_eσ), prior ϑ⁰; ρ is the uniform time grid of one period
    (n_quad points)."""

    def __init__(self, net, f_vec, x_ref, period, G, sigma0, theta0=None, n_quad=256):
        if net.n != 1:
            raise NotImplementedError("the averaged slow loop is implemented for n = 1")
        self.net = net
        self.G = np.asarray(G, float)
        self.sqG = np.sqrt(self.G)
        self.sigma0 = float(sigma0)
        self.theta0 = np.zeros(3) if theta0 is None else np.asarray(theta0, float)
        self.period = float(period)
        self.N = int(n_quad)
        ts = np.arange(self.N) * self.period / self.N
        self.X = np.array([x_ref(t)[0] for t in ts])                     # (N,)
        self.F = np.asarray(f_vec(self.X), float)                        # (N,)
        self.F_rho = float(np.max(np.abs(self.F)))

    # ------------------------------------------------------------------ features and the objective
    def features(self, Th):
        a, c, s = self.net.physical(Th)
        diff = self.X[:, None] - c[None, :, 0]                            # (N, m)
        phi = np.exp(-diff ** 2 / (2 * s[None, :] ** 2))
        return a[:, 0], c[:, 0], s, diff, phi

    def f_hat(self, Th):
        a, c, s, diff, phi = self.features(Th)
        return phi @ a / Th.shape[0]

    def residual(self, Th):
        return self.F - self.f_hat(Th)

    def grad_weighted(self, r, Th):
        """(1/N)Σ_x r(x) ∇_ϑΦ̃(x; ϑ_i) in the unconstrained coordinates (chain rule through the saturations), (m, 3)."""
        a, c, s, diff, phi = self.features(Th)
        net = self.net
        da = 1.0 - (a / net.a_max) ** 2
        dc = 1.0 - (c / net.c_max) ** 2
        sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
        out = np.empty_like(Th)
        out[:, 0] = (phi.T @ r) * da
        out[:, 1] = a * dc / s ** 2 * ((phi * diff).T @ r)
        out[:, 2] = a * ds / s ** 3 * ((phi * diff ** 2).T @ r)
        return out / self.N

    def grad_V_L(self, r, Th):
        """∇_ϑV_L[μ](ϑ_i) = −⟨r_μ ∇_ϑΦ̃(·;ϑ_i)⟩_ρ for a frozen residual r, (m, 3)."""
        return -self.grad_weighted(r, Th)

    def grad_J(self, Th):
        """∇_{ϑ_i}J_m(Θ) = ∇V_L[μ_Θ](ϑ_i) + σ₀G⁻¹(ϑ_i − ϑ⁰)  (the same object as periodic.ClosedLoop.grad_J_m)."""
        return self.grad_V_L(self.residual(Th), Th) + self.sigma0 * (Th - self.theta0[None, :]) / self.G[None, :]

    def drift(self, Th):
        """dΘ/ds = −G∇J_m(Θ) = G⟨r_Θ∇_ϑΦ̃(·;ϑ_i)⟩_ρ − σ₀(ϑ_i − ϑ⁰), (m, 3); the features are evaluated once."""
        a, c, s, diff, phi = self.features(Th)
        m = Th.shape[0]
        r = self.F - phi @ a / m
        net = self.net
        da = 1.0 - (a / net.a_max) ** 2
        dc = 1.0 - (c / net.c_max) ** 2
        sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
        g = np.empty_like(Th)
        g[:, 0] = (phi.T @ r) * da
        g[:, 1] = a * dc / s ** 2 * ((phi * diff).T @ r)
        g[:, 2] = a * ds / s ** 3 * ((phi * diff ** 2).T @ r)
        return self.G[None, :] * g / self.N - self.sigma0 * (Th - self.theta0[None, :])

    def L_bar(self, Th):
        r = self.residual(Th)
        return 0.5 * float(np.mean(r ** 2))

    def fit_error(self, Th):
        """‖f − f̂_Θ‖_{L²(ρ)}/‖f‖_{L²(ρ)}."""
        r = self.residual(Th)
        return float(np.sqrt(np.mean(r ** 2) / np.mean(self.F ** 2)))

    # ------------------------------------------------------------------ the residual-weighted curvature
    def residual_curvature(self, r, Th, h=1e-5):
        """R[μ](ϑ_i) = ∇²_ϑ⟨r Φ̃(·;ϑ)⟩_ρ at ϑ = ϑ_i with r frozen, (m, 3, 3): central differences of grad_weighted."""
        m = Th.shape[0]
        R = np.zeros((m, 3, 3))
        for j in range(3):
            E = np.zeros_like(Th); E[:, j] = h
            R[:, :, j] = (self.grad_weighted(r, Th + E) - self.grad_weighted(r, Th - E)) / (2 * h)
        return 0.5 * (R + np.swapaxes(R, 1, 2))

    def kappa_R(self, r, Th, h=1e-5):
        """λ_max(G^{1/2}R[μ](ϑ_i)G^{1/2}) for every row, (m,)."""
        R = self.residual_curvature(r, Th, h)
        Rw = self.sqG[None, :, None] * R * self.sqG[None, None, :]
        return np.linalg.eigvalsh(Rw)[:, -1]

    def kappa_R_neighbourhood(self, r, Th, rho, n_dirs=8, rng=None, h=1e-5):
        """sup of λ_max(G^{1/2}RG^{1/2}) over the rows of Th and over n_dirs random points of the G⁻¹-ball of radius rho
        around each row (a Monte Carlo estimate of κ_R[μ](N_ρ(supp μ_Θ)) from below; the maximum over the rows is exact)."""
        rng = np.random.default_rng(0) if rng is None else rng
        best = float(self.kappa_R(r, Th, h).max())
        for _ in range(n_dirs):
            u = rng.standard_normal(Th.shape)
            u /= np.linalg.norm(u, axis=1, keepdims=True)
            rad = rho * rng.uniform(0, 1, (Th.shape[0], 1)) ** (1.0 / 3.0)
            best = max(best, float(self.kappa_R(r, Th + rad * u * self.sqG[None, :], h).max()))
        return best

    # ------------------------------------------------------------------ the three terms of the coupling identity
    def decomposition(self, Th, Thb, fhat_mf=None):
        """
        For the m-system Θ and the reference Θ̄ (same shape) with f̂[μ_s] = fhat_mf on the grid (default f̂_Θ̄, in which case
        T_C = 0): the three terms T_A, T_B, T_C, their sum, the exact ½dD/ds = −(1/m)Σ⟨δ_i, G[∇𝒱[μ_Θ](ϑ_i) − ∇𝒱[μ_s](ϑ̄_i)]⟩_{G⁻¹}
        (the mean-field drift at ϑ̄_i uses the residual f − fhat_mf), D, ‖h‖_{L²(ρ)} with h := f̂_Θ − f̂_Θ̄, and the bound
        (L₂^G)²D²/16 on T_B requires L₂^G (returned as bound_B_coeff·D² with bound_B_coeff = 1/16; multiply by (L₂^G)²).
        """
        d = Th - Thb
        rT, rB = self.residual(Th), self.residual(Thb)
        r_mf = rB if fhat_mf is None else self.F - np.asarray(fhat_mf, float)
        D = float(np.mean(np.sum(d ** 2 / self.G[None, :], axis=1)))
        # T_A: self term with the reference residual (leak included)
        gA = self.grad_V_L(rB, Th) - self.grad_V_L(rB, Thb) + self.sigma0 * d / self.G[None, :]
        T_A = -float(np.mean(np.sum(d * gA, axis=1)))
        # T_B: interaction, −⟨δ_i, ∇V_L[μ_Θ](ϑ_i) − ∇V_L[μ_Θ̄](ϑ_i)⟩ = −⟨δ_i, −⟨(rT − rB)∇Φ̃(ϑ_i)⟩⟩
        gB = self.grad_V_L(rT, Th) - self.grad_V_L(rB, Th)
        T_B = -float(np.mean(np.sum(d * gB, axis=1)))
        # T_C: sampling fluctuation of the reference
        gC = self.grad_V_L(rB, Thb) - self.grad_V_L(r_mf, Thb)
        T_C = -float(np.mean(np.sum(d * gC, axis=1)))
        # exact half-derivative: drift difference (G-preconditioned) paired with δ in the G⁻¹ metric
        drift_m = -self.G[None, :] * (self.grad_V_L(rT, Th) + self.sigma0 * (Th - self.theta0[None, :]) / self.G[None, :])
        drift_mf = -self.G[None, :] * (self.grad_V_L(r_mf, Thb) + self.sigma0 * (Thb - self.theta0[None, :]) / self.G[None, :])
        half_dD = float(np.mean(np.sum(d * (drift_m - drift_mf) / self.G[None, :], axis=1)))
        h = self.f_hat(Th) - self.f_hat(Thb)
        return dict(T_A=T_A, T_B=T_B, T_C=T_C, total=T_A + T_B + T_C, half_dD=half_dD, D=D,
                    h_norm=float(np.sqrt(np.mean(h ** 2))), bound_B_coeff=1.0 / 16.0)

    # ------------------------------------------------------------------ the coupled run
    def coupled_run(self, Th_all, m_list, S_end, dt, rec_every=25, rho_nbhd=None, n_dirs=4, rng=None, record_final=True):
        """
        The m-systems (first m rows of Th_all) coupled synchronously to the M-row reference (M = Th_all.shape[0], the
        mean-field flow's stand-in, as in Theorems D′/D″), integrated by the midpoint rule in slow time.  Records, every
        rec_every steps: s, D_m (mean square, G⁻¹ metric), the reference's fit error, physical rms spread and rms drift,
        κ_R at the reference particles (max) and, with rho_nbhd, over their G⁻¹-neighbourhood of that radius (Monte Carlo),
        κ_R at the m-particles (max), the function-level distances ‖f̂_Θ − f̂_Θ̄‖_{L²(ρ)} (the same m particles; bounded by
        L_ϑ^G·mean‖δ_i‖) and ‖f̂_Θ − f̂_ref‖_{L²(ρ)} (against the M-reference, the stand-in for f̂[μ_s]), and the mean of ‖δ_i‖_{G⁻¹}.
        """
        rng = np.random.default_rng(0) if rng is None else rng
        M = Th_all.shape[0]
        m_list = [int(m) for m in m_list]
        sys_ = {m: Th_all[:m].copy() for m in m_list + [M]}
        n_steps = int(round(S_end / dt))
        rec = dict(s=[], fit=[], spread=[], drift_rms=[], kR_ref=[], kR_ref_nbhd=[],
                   D={m: [] for m in m_list}, kR_m={m: [] for m in m_list}, hnorm={m: [] for m in m_list},
                   hnorm_mf={m: [] for m in m_list}, dev1={m: [] for m in m_list})
        for k in range(n_steps):
            new = {}
            for m, Th in sys_.items():
                k1 = self.drift(Th)
                k2 = self.drift(Th + 0.5 * dt * k1)
                new[m] = Th + dt * k2
            sys_ = new
            if k % rec_every == 0 or k == n_steps - 1:
                s = (k + 1) * dt
                ref = sys_[M]
                r_ref = self.residual(ref)
                rec['s'].append(s)
                rec['fit'].append(self.fit_error(ref))
                a, c, sw = self.net.physical(ref); P = np.column_stack([a[:, 0], c[:, 0], sw])
                rec['spread'].append(float(np.sqrt(np.mean(np.sum((P - P.mean(0)) ** 2, axis=1)))))
                rec['drift_rms'].append(float(np.sqrt(np.mean(np.sum(self.drift(ref) ** 2 / self.G[None, :], axis=1)))))
                kr = self.kappa_R(r_ref, ref)
                rec['kR_ref'].append(float(kr.max()))
                rec['kR_ref_nbhd'].append(self.kappa_R_neighbourhood(r_ref, ref, rho_nbhd, n_dirs, rng) if rho_nbhd else float('nan'))
                fh_ref = self.f_hat(ref)
                for m in m_list:
                    d = sys_[m] - ref[:m]
                    dn = np.sqrt(np.sum(d ** 2 / self.G[None, :], axis=1))
                    rec['D'][m].append(float(np.mean(dn ** 2)))
                    rec['dev1'][m].append(float(np.mean(dn)))
                    rec['kR_m'][m].append(float(self.kappa_R(r_ref, sys_[m]).max()))
                    fh_m = self.f_hat(sys_[m])
                    rec['hnorm'][m].append(float(np.sqrt(np.mean((fh_m - self.f_hat(ref[:m])) ** 2))))      # ‖f̂_Θ − f̂_Θ̄‖ (same m particles)
                    rec['hnorm_mf'][m].append(float(np.sqrt(np.mean((fh_m - fh_ref) ** 2))))              # ‖f̂_Θ − f̂_ref‖ (the M-reference as f̂[μ_s])
        for key in ('s', 'fit', 'spread', 'drift_rms', 'kR_ref', 'kR_ref_nbhd'):
            rec[key] = np.array(rec[key])
        for key in ('D', 'kR_m', 'hnorm', 'hnorm_mf', 'dev1'):
            rec[key] = {m: np.array(v) for m, v in rec[key].items()}
        if record_final:
            ref = sys_[M]
            rec['final'] = sys_
            rec['dev_final'] = {m: np.sqrt(np.sum((sys_[m] - ref[:m]) ** 2 / self.G[None, :], axis=1)) for m in m_list}
        return rec


# ---------------------------------------------------------------------- constants of the theorem
def whitened_lipschitz_constants(net, G, x_box, n_samples=60000, seed=0, h=1e-5):
    """
    L_ϑ^G := sup_{x,ϑ}‖∇_ϑΦ̃(x;ϑ)G^{1/2}‖ and L₂^G := sup_{x,ϑ}‖G^{1/2}∇²_ϑΦ̃(x;ϑ)G^{1/2}‖ (n = 1) for the saturated features,
    by a search over sampled ϑ (physical and unconstrained) and a grid in |x| ≤ x_box with a Nelder–Mead polish, together
    with closed-form upper bounds: L_ϑ^G ≤ (G_a + G_c(a_max e^{−1/2}/s_min)² + G_s(a_max e^{−1}(s_max − s_min)/(2s_min))²)^{1/2}
    (per-coordinate suprema |∂_αΦ̃| ≤ 1, |∂_γΦ̃| ≤ a_max e^{−1/2}/s_min, |∂_κΦ̃| ≤ a_max·2e^{−1}/s_min·(s_max − s_min)/4) and
    L₂^G ≤ λ_G·a_max/s_min² (the closed form of `exp_chaos_uniform.lipschitz_constants_exact`).
    """
    from scipy.optimize import minimize
    G = np.asarray(G, float); sqG = np.sqrt(G)
    rng = np.random.default_rng(seed)
    a = rng.uniform(-net.a_max, net.a_max, (n_samples // 2, 1)); c = rng.uniform(-x_box, x_box, (n_samples // 2, 1))
    s = np.exp(rng.uniform(np.log(net.s_min), np.log(net.s_max), n_samples // 2))
    Th = np.vstack([net.unconstrained(a, c, s), rng.normal(0, 3, (n_samples // 2, 3))])

    def grad(x, T):                                                        # (m,3) ∇_ϑΦ̃(x;ϑ)
        return net.grad_g(np.array([x]), np.array([1.0]), T)

    def L1(x, T):
        return np.linalg.norm(grad(x, T) * sqG[None, :], axis=1)

    def L2(x, T):
        Hs = []
        for j in range(3):
            E = np.zeros_like(T); E[:, j] = h
            Hs.append((grad(x, T + E) - grad(x, T - E)) / (2 * h))
        H = np.stack(Hs, axis=2); H = 0.5 * (H + np.swapaxes(H, 1, 2))
        H = sqG[None, :, None] * H * sqG[None, None, :]
        return np.linalg.norm(H, ord=2, axis=(1, 2))
    out = {}
    for name, fun in (('L_th_G', L1), ('L_2_G', L2)):
        best, arg = 0.0, None
        for x in np.linspace(-x_box, x_box, 121):
            v = fun(x, Th)
            i = int(np.argmax(v))
            if v[i] > best:
                best, arg = float(v[i]), (x, Th[i].copy())
        neg = lambda z: -float(fun(float(np.clip(z[0], -x_box, x_box)), z[1:][None, :])[0])
        r = minimize(neg, np.concatenate([[arg[0]], arg[1]]), method='Nelder-Mead', options=dict(xatol=1e-9, fatol=1e-12, maxiter=20000))
        out[name] = max(best, -float(r.fun))
    ds_max = (net.s_max - net.s_min) / 4.0
    out['L_th_G_closed'] = float(np.sqrt(G[0] * 1.0 + G[1] * (net.a_max * np.exp(-0.5) / net.s_min) ** 2
                                         + G[2] * (net.a_max * 2 * np.exp(-1.0) / net.s_min * ds_max) ** 2))
    out['L_2_G_closed'] = float(G.max() * net.a_max / net.s_min ** 2)
    return out


def theorem_constants(loop, L_th_G, L_2_G, R0, kappa, Lambda1, rho_prime):
    """
    The constants of Theorem D‴ for the averaged slow loop `loop`: K_Φ = a_max (n = 1), the a-priori radius
    R_Θ = max{R₀, L_ϑ^G(F_ρ + K_Φ)/σ₀} and D_max = 4R_Θ², the speed bound V_max = L_ϑ^G(F_ρ + K_Φ) + σ₀R_Θ, the thresholds
    ϱ and θ of the bootstrap, the grid step Δ = θ/(4L_ϑ^G V_max), the main constant C₁ = 4e^{2Λ₁}(L_ϑ^G)²K_Φ²/κ² (so that
    sup_s E D ≤ C₁/m + tail), the sample size m₀ = 16K_Φ²/θ² from which the tail bound applies, and the tail exponent
    θ²/(32K_Φ²) (P(𝒢^c) ≤ (S/Δ + 2)exp(−m θ²/(32K_Φ²)) over a horizon S).
    """
    K = loop.net.K_Phi
    F = loop.F_rho
    R_Theta = max(R0, L_th_G * (F + K) / loop.sigma0)
    D_max = 4 * R_Theta ** 2
    V_max = L_th_G * (F + K) + loop.sigma0 * R_Theta
    rho = min(np.sqrt(2 * kappa) / L_2_G, rho_prime * kappa * np.exp(-Lambda1) / (4 * L_th_G ** 2))
    theta = min(kappa / (8 * L_2_G), rho * kappa * np.exp(-Lambda1) / (4 * L_th_G), rho_prime * kappa * np.exp(-Lambda1) / (4 * L_th_G))
    Delta = theta / (4 * L_th_G * V_max)
    C1 = 4 * np.exp(2 * Lambda1) * L_th_G ** 2 * K ** 2 / kappa ** 2
    return dict(K_Phi=K, F_rho=F, R_Theta=R_Theta, D_max=D_max, V_max=V_max, rho=rho, theta=theta, Delta=Delta, C1=C1,
                m0=16 * K ** 2 / theta ** 2, tail_exponent=theta ** 2 / (32 * K ** 2))


def Lambda1_from_profile(s, a, kappa, s1):
    """Λ₁ := ∫₀^{s₁}(κ − a(s))₊ ds from a sampled profile a(s) (trapezoid)."""
    s = np.asarray(s, float); a = np.asarray(a, float)
    mask = s <= s1
    if mask.sum() < 2:
        return 0.0
    return float(np.trapezoid(np.maximum(kappa - a[mask], 0.0), s[mask]))


# ---------------------------------------------------------------------- initial clouds
def cold_start(net, M, rng, a_sd=0.8, c_range=2.0, s_range=(0.4, 1.0)):
    """A spread i.i.d. cloud: a ~ N(0, a_sd²), c ~ U(−c_range, c_range), s ~ U(s_range)."""
    return net.unconstrained(rng.normal(0, a_sd, (M, 1)), rng.uniform(-c_range, c_range, (M, 1)), rng.uniform(*s_range, M))


def warm_start(net, bumps, M, rng, a_sd=0.15, c_sd=0.1, s_logsd=0.1):
    """An i.i.d. cloud near an exact representation of a sum of K Gaussian bumps (A_k, C_k, S_k) by K equal clusters: each
    particle picks a bump uniformly and sits at (K·A_k, C_k, S_k) up to the given jitter (the weights carry the factor K
    because the network averages over its kernels)."""
    K = len(bumps)
    kk = rng.integers(0, K, M)
    A = np.array([b[0] for b in bumps]); C = np.array([b[1] for b in bumps]); S = np.array([b[2] for b in bumps])
    a0 = K * A[kk] + rng.normal(0, a_sd, M); c0 = C[kk] + rng.normal(0, c_sd, M); s0 = S[kk] * np.exp(rng.normal(0, s_logsd, M))
    return net.unconstrained(a0[:, None], c0[:, None], s0)
