"""
Measurement noise at perfect tracking: a coherent, high-passed noise, and the noise floors of the fixed and adaptive priors
============================================================================================================================

The real closed loop (CL) of Theorem J with the measurement y = x + ν in place of x in the controller and the adaptive law,
u = ẋ_m − k_e(y − x_m) − f̂_Θ(y), ϑ̇_i = G∇_ϑ((y − x_m)ᵀΦ̃(y; ϑ_i)) − (σ₀/k_e)(ϑ_i − T_i), ν a stationary noise of spectral density N₀ at the
frequencies of the regressor (white at those frequencies). Linearised in ν at the lift of an exact representation (Corollary D‴.26), in the
whitened weighted coordinates of `adaptive_floquet.CollapsedVariational`,

    δė = −k_eδe − Ĵ(t)y − (k_e + a(t))ν,   ẏ = Ĵ(t)ᵀ(δe + ν) − (σ₀/k_e)(y − ζ),   ζ̇ = −(η/k_e)(ζ − y),   a(t) := ∂_xf̂(x_m(t)):

the measured error's noise is ((d/dt − a)/(d/dt + k_e))ν, high-passed by the loop, and in the slow time s = t/k_e the collective
coordinates diffuse with the matrix (N₀/k_e)Q, Q := ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩ (Proposition D‴.30). Every particle of a class receives the same
noise increment, so the separations are not excited at first order: the noise is coherent and its floor does not average out with m; at
second order the separations' stability needs the noise below a level (`separation_drift`, `second_order_transient`). This module gives
Q and H = ⟨ĴᵀĴ⟩, the stationary covariance of the averaged collective block (the Lyapunov equation), the closed-form function-level floors
(N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) in the eigenbasis of H, and the exact periodic stationary covariance of the linearised loop.

Proposition D‴.31: the adaptive prior is a washout leak — on the class barycentre the leak acts through the high-pass iω/(iω + η) — and
along the prior family dz/ds = −η(z − T̂) − κ(z − z_a) the zero-frequency response R and bias B satisfy RH + B = I (`collective_block`
with κ, `zero_frequency_leak`, `stationary_operators`, `washout_transfer`; a prior fixed on some coordinates, an active face or a matrix
anchor gain, the operators from the block: `anchored_block`, `anchored_operators`, `anchored_leak`; a general zero-frequency leak by its
definition: `zero_frequency_operators`); a
persistent force's function-level response is at most FᵀH⁻¹F, the smallest equivalent target perturbation (`equivalent_perturbation`,
`response_size`); and with the regressor at the reference (`NoisyRealLoop(regressor='reference')`, and the compensation too) the second-
order force of the noise's own moments is the least-squares force of g = k_eE[e] (`eiv_force(variant=…)`, `mean_error`; the linearisation
with the compensation at the reference: `reference_compensation`).

The full second-order mean (round 21 of the review): `second_order_mean` contracts the periodic covariance of the linearised loop (the OU
noise a state) with the true field's second derivatives — the noise part is the quasi-static force above, the parameters' part the
compensation's curvature and the law's ∇²_ϑΦ̃E[ε_mδϑ] — and `mean_biases` gives the stationary biases of the mean parameters and of the
mean function; `separation_drift` the separations' second-order drift S_k at the lift, whose largest eigenvalue against σ₀/k_e decides
whether the replicated configuration attracts at the start; along the noisy trajectory the mean error also carries the mean tracking error
of the errors-in-variables deviation as it develops (round 22 of the review, ID-1): `second_order_transient` gives the second-order mean and
moments from the lift and `separation_drift_series` the drift along them.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov, solve_discrete_lyapunov

from . import adaptive_floquet as AF


def noise_matrices(cv, fprime, N=1024):
    """(H, Q) of the collapsed variational equation `cv` (a CollapsedVariational, n = 1): H = ⟨ĴᵀĴ⟩ and Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩ with
    a(t) = fprime(x_m(t)) = ∂_xf̂ along the orbit (at an exact representation ∂_xf), from N samples (spectral derivative)."""
    t = np.arange(N) * cv.T0 / N
    J = np.array([cv.J(s)[0] for s in t])
    k = np.fft.fftfreq(N, d=1.0 / N) * (2 * np.pi / cv.T0)
    Jd = np.real(np.fft.ifft((1j * k)[:, None] * np.fft.fft(J, axis=0), axis=0))
    a = np.array([float(np.atleast_1d(fprime(np.atleast_1d(cv.xr(s))[0]))[0]) for s in t])
    V = Jd + a[:, None] * J
    return J.T @ J / N, V.T @ V / N


def fast_parts(cv, N=1024):
    """(tr H², ⟨|Ĵ|⁴⟩): the fast part of y — the measured error's high frequencies pass the loop, y_fast ≈ Ĵ(t)ᵀξ(t) with ξ the integrated
    high-passed noise, of variance N₀/(2k_e) for white noise — adds (N₀/2k_e)tr H² to the function-level floor E‖f̂ − f‖²_ρ (the error over
    the orbit at a fixed time) and (N₀/2k_e)⟨|Ĵ|⁴⟩ to the compensation error at the plant's current point, whatever the prior."""
    t = np.arange(N) * cv.T0 / N
    J = np.array([cv.J(s)[0] for s in t])
    H = J.T @ J / N
    return float(np.trace(H @ H)), float(np.mean(np.sum(J ** 2, axis=1) ** 2))


def floors(H, Q, trH2, J4, sigma0, eta, N0, ke):
    """(the function-level floor E‖f̂ − f‖²_ρ, the compensation error's floor at the current point) for white measurement noise:
    (N₀/2k_e)[Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) + tr H²] and (N₀/2k_e)[Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) + ⟨|Ĵ|⁴⟩], up to O(k_e⁻²)."""
    slow = floor_closed_form(H, Q, sigma0, eta, N0, ke)[0]
    c = N0 / (2 * ke)
    return slow + c * trH2, slow + c * J4


def collective_block(H, sigma0, eta, kappa=0.0):
    """The averaged collective block in slow time (the drift matrix of (y, ζ); η = ∞: −H on y alone; η = 0: the fixed prior, −(H + σ₀I)).
    `kappa` > 0 (finite η): the prior leaked toward its fixed anchor as well, dz/ds = −η(z − T̂) − κ(z − z_a) (Proposition D‴.31(ii))."""
    d = H.shape[0]
    if np.isinf(eta):
        return -H
    if eta == 0:
        return -(H + sigma0 * np.eye(d))
    I = np.eye(d)
    return np.block([[-(H + sigma0 * I), sigma0 * I], [eta * I, -(eta + kappa) * I]])


def zero_frequency_leak(sigma0, eta, kappa=0.0):
    """The leak's zero-frequency gain S₀ of the prior family dz/ds = −η(z − T̂) − κ(z − z_a): the leak on the barycentre is the σ-modification
    passed through S(iω)/σ₀ = (iω + κ)/(iω + η + κ) — at κ = 0 the washout iω/(iω + η) — so S₀ = σ₀κ/(η + κ) (η = 0: σ₀, the fixed prior)."""
    if eta == 0:
        return sigma0
    if np.isinf(eta):
        return 0.0
    return sigma0 * kappa / (eta + kappa)


def stationary_operators(H, sigma0, eta, kappa=0.0):
    """(R, B) of the averaged collective block at zero frequency: the response y = RF to a persistent force F on y, and the bias y = Bδ toward
    a misplaced anchor (the prior's fixed point displaced by δ: the fixed prior at δ, or the anchor z_a = δ), from the block itself — so that
    RH + B = I can be checked rather than assumed (here S₀ is a scalar, so HR + B = I as well)."""
    d = H.shape[0]
    A = collective_block(H, sigma0, eta, kappa)
    D = A.shape[0]
    R = -np.linalg.solve(A, np.vstack([np.eye(d), np.zeros((D - d, d))]))[:d]
    if eta == 0:
        Bm = -np.linalg.solve(A, sigma0 * np.eye(d))
    elif np.isinf(eta):
        Bm = np.zeros((d, d))
    else:
        Bm = -np.linalg.solve(A, np.vstack([np.zeros((d, d)), kappa * np.eye(d)]))[:d]
    return R, Bm


def zero_frequency_operators(H, S0):
    """(R, B) = ((H + S₀)⁻¹, (H + S₀)⁻¹S₀) for a general zero-frequency leak S₀ (a matrix: a prior fixed on some coordinates, the linearised
    projected prior with an active face): RH + B = I always, HR + B = I only when S₀ commutes with H; per mode h_jr_jj + b_jj = 1 in the
    eigenbasis of H whatever S₀ (r_jj = v_jᵀRv_j, b_jj = v_jᵀBv_j)."""
    R = np.linalg.inv(H + S0)
    return R, R @ S0


def anchored_block(H, sigma0, eta, K=None, fixed=None):
    """The averaged collective block of a prior that is held at its anchor on the coordinates `fixed` (a boolean mask: a prior fixed there,
    or a projected prior's active face) and elsewhere moves by dz_f/ds = −η(z_f − y_f) − K(z_f − z_a,f), K a (f, f) anchor gain (None: 0),
    in whitened coordinates (round 22, ID-5). Returns (A, Bin): the drift matrix of (y, z_f) (η = ∞: z_f = y_f, the state y alone) and the
    input of the anchor displacement δ (the anchor z_a = δ: σ₀δ_x on y through the fixed coordinates, Kδ_f on z_f)."""
    d = H.shape[0]
    fixed = np.zeros(d, bool) if fixed is None else np.asarray(fixed, bool)
    Ex, Ef = np.eye(d)[:, fixed], np.eye(d)[:, ~fixed]
    f = Ef.shape[1]
    K = np.zeros((f, f)) if K is None else np.asarray(K, float)
    if np.isinf(eta):
        return -(H + sigma0 * Ex @ Ex.T), sigma0 * Ex @ Ex.T
    A = np.block([[-(H + sigma0 * np.eye(d)), sigma0 * Ef], [eta * Ef.T, -(eta * np.eye(f) + K)]])
    Bin = np.vstack([sigma0 * Ex @ Ex.T, K @ Ef.T])
    return A, Bin


def anchored_operators(H, sigma0, eta, K=None, fixed=None):
    """(R, B) of `anchored_block` at zero frequency, from the block itself: the response y = RF to a persistent force F on y and the bias
    y = Bδ toward the displaced anchor."""
    d = H.shape[0]
    A, Bin = anchored_block(H, sigma0, eta, K, fixed)
    D = A.shape[0]
    R = -np.linalg.solve(A, np.vstack([np.eye(d), np.zeros((D - d, d))]))[:d]
    return R, -np.linalg.solve(A, Bin)[:d]


def anchored_leak(H, sigma0, eta, K=None, fixed=None):
    """The zero-frequency leak of `anchored_block` by elimination of z_f: S₀ = σ₀[E_xE_xᵀ + E_f(ηI + K)⁻¹KE_fᵀ] (η = ∞: σ₀E_xE_xᵀ), a matrix
    that commutes with H only exceptionally; S₀ + S₀ᵀ ⪰ 0 when K + Kᵀ ⪰ 0."""
    d = H.shape[0]
    fixed = np.zeros(d, bool) if fixed is None else np.asarray(fixed, bool)
    Ex, Ef = np.eye(d)[:, fixed], np.eye(d)[:, ~fixed]
    S0 = sigma0 * Ex @ Ex.T
    if not np.isinf(eta) and K is not None:
        f = Ef.shape[1]
        S0 = S0 + sigma0 * Ef @ np.linalg.solve(eta * np.eye(f) + K, K) @ Ef.T
    return S0


def washout_transfer(h, sigma0, eta, omega, kappa=0.0):
    """The transfer function of the pair block of a mode (stiffness h) from a force on y at the slow frequency ω to y:
    1/(iω + h + σ₀(iω + κ)/(iω + η + κ)) — the stiffness h + S(iω), the leak washed out at zero frequency when κ = 0."""
    w = 1j * np.asarray(omega, float)
    return 1.0 / (w + h + sigma0 * (w + kappa) / (w + eta + kappa))


def averaged_covariance(H, Q, sigma0, eta, N0, ke):
    """The stationary covariance of y under the averaged collective block with the diffusion (N₀/k_e)Q on y (the Lyapunov equation)."""
    A = collective_block(H, sigma0, eta)
    d = H.shape[0]
    Qb = np.zeros_like(A); Qb[:d, :d] = (N0 / ke) * Q
    return solve_continuous_lyapunov(A, -Qb)[:d, :d]


def floor_closed_form(H, Q, sigma0, eta, N0, ke):
    """The function-level noise floor tr(HΣ_yy) = (N₀/2k_e)Σ_jQ̃_jj(h_j + η)/(h_j + η + σ₀) (Q̃ = VᵀQV in the eigenbasis of H; η = ∞: the
    factor 1, so (N₀/2k_e)tr Q whatever H; η = 0: the fixed prior's Tikhonov filter factor h_j/(h_j + σ₀)), and the per-mode parameter
    variances (N₀/2k_e)Q̃_jj(h_j + η)/(h_j(h_j + η + σ₀)). Returns (floor, per-mode function-level shares, per-mode parameter variances)."""
    h, V = np.linalg.eigh(H)
    q = np.einsum('ij,ik,kj->j', V, Q, V)
    c = N0 / (2 * ke)
    fac = np.ones_like(h) if np.isinf(eta) else (h + eta) / (h + eta + sigma0)
    shares = c * q * fac
    return float(shares.sum()), shares, c * q * fac / h


def periodic_covariance(gen, bvec, period, D, N0, rtol=1e-11):
    """The exact stationary T₀-periodic covariance at t = 0 of dX = M(t)X dt + b(t)√N₀ dW (gen(t) = M(t), bvec(t) = b(t), dimension D):
    P(0) = ΦP(0)Φᵀ + Q_T, Φ the monodromy and Q_T the covariance accumulated over a period from 0 (a discrete Lyapunov equation)."""
    f = lambda t, y: (gen(t) @ y.reshape(D, D)).ravel()
    Phi = solve_ivp(f, (0.0, period), np.eye(D).ravel(), method='DOP853', rtol=1e-12, atol=1e-15).y[:, -1].reshape(D, D)

    def rhs(t, p):
        P = p.reshape(D, D); M = gen(t); b = bvec(t)
        return (M @ P + P @ M.T + N0 * np.outer(b, b)).ravel()
    QT = solve_ivp(rhs, (0.0, period), np.zeros(D * D), method='DOP853', rtol=rtol, atol=1e-16).y[:, -1].reshape(D, D)
    return solve_discrete_lyapunov(Phi, QT)


def transient_covariance(gen, bvec, D, N0, t_end, t_eval, P0=None, rtol=1e-9):
    """The covariance P(t) of the linearised loop from P(0) = P0 (default 0) at the times t_eval (the covariance ODE)."""
    def rhs(t, p):
        P = p.reshape(D, D); M = gen(t); b = bvec(t)
        return (M @ P + P @ M.T + N0 * np.outer(b, b)).ravel()
    p0 = np.zeros(D * D) if P0 is None else np.asarray(P0, float).ravel()
    sol = solve_ivp(rhs, (0.0, t_end), p0, method='LSODA', t_eval=t_eval, rtol=rtol, atol=1e-18)
    return sol.t, sol.y.T.reshape(-1, D, D)


def collapsed_noise_input(cv, fprime, ke, compensation='measured'):
    """b(t) of the collapsed linearised loop (δe, y[, ζ]): (−(k_e + a(t)), Ĵ(t)ᵀ, 0); with the compensation at the reference, (−k_e, Ĵ(t)ᵀ, 0)
    (the compensation no longer sees ν). The regressor's point does not enter at first order."""
    def b(t):
        J = cv.J(t)[0]
        out = np.zeros(cv.dim)
        a = 0.0 if compensation == 'reference' else float(np.atleast_1d(fprime(np.atleast_1d(cv.xr(t))[0]))[0])
        out[0] = -(ke + a)
        out[1:1 + cv.d] = J
        return out
    return b


def reference_compensation(gen, cv, fprime, ke):
    """The linearised loop with the compensation at the reference (u ∋ −f̂_Θ(x_m)): the error channel −(k_e − a(t))δe in place of −k_eδe
    (f(x) − f̂(x_m) = aδe − Ĵy), for `gen` either collapsed generator; the noise input's error entry is then −k_e
    (`collapsed_noise_input(…, compensation='reference')`)."""
    def M(t):
        out = gen(t).copy()
        out[0, 0] += float(np.atleast_1d(fprime(np.atleast_1d(cv.xr(t))[0]))[0])
        return out
    return M


def fixed_generator(cv, ke):
    """The collapsed linearised loop with a fixed prior at the atoms (η = 0, no ζ): ẏ = Ĵᵀδe − (σ₀/k_e)y."""
    def M(t):
        J = cv.J(t)[0]; d = cv.d
        out = np.zeros((1 + d, 1 + d))
        out[0, 0] = -ke; out[0, 1:] = -J; out[1:, 0] = J; out[1:, 1:] = -(cv.sigma0 / ke) * np.eye(d)
        return out
    return M


# ---------------------------------------------------------------------- the real loop with measurement noise
class NoisyRealLoop(AF.AdaptiveRealLoop):
    """`AdaptiveRealLoop` with the measurement y = x + ν used by the controller and the adaptive law: ẋ = f(x) + ẋ_m − k_e(y − x_m) − f̂_Θ(x_c),
    ϑ̇_i = G∇_ϑ((y − x_m)ᵀΦ̃(x_r; ϑ_i)) − (σ₀/k_e)(ϑ_i − z_{π(i)}), the prior as in AdaptiveRealLoop (η = 0: a fixed prior at the initial Z).
    The regressor's point x_r and the compensation's point x_c are the measurement y (`regressor`, `compensation` = 'measured', the default:
    Proposition D‴.30) or the reference x_m ('reference': the desired regressor of Proposition D‴.31, and with both the desired compensation
    adaptation law)."""

    def __init__(self, *args, regressor='measured', compensation='measured', **kw):
        super().__init__(*args, **kw)
        assert regressor in ('measured', 'reference') and compensation in ('measured', 'reference')
        self.regressor, self.compensation = regressor, compensation

    def rhs_noisy(self, t, z, m, nu):
        n, p = self.net.n, self.p
        x, Th = z[:n], z[n:n + m * p].reshape(m, p)
        if np.isinf(self.eta):
            Z = self.class_means(Th)
        else:
            Z = z[n + m * p:].reshape(self.K0, p)
        y = x + nu
        xr, xrd = self.xr(t), self.xrd(t)
        em = y - xr
        fx = self.f_vec(x[None, :])[0]
        xc = y if getattr(self, 'compensation', 'measured') == 'measured' else xr
        xa = y if getattr(self, 'regressor', 'measured') == 'measured' else xr
        xdot = fx + xrd - self.ke * em - self.f_hat(xc, Th)
        Thdot = self.G[None, :] * self.net.grad_g(xa, em, Th) - (self.sigma0 / self.ke) * (Th - Z[self.labels])
        if np.isinf(self.eta):
            return np.concatenate([xdot, Thdot.ravel()])
        Zdot = -(self.eta / self.ke) * (Z - self.class_means(Th))
        return np.concatenate([xdot, Thdot.ravel(), Zdot.ravel()])


def simulate_noisy(cl, z0, m, T, dt, B, var, rng, rec_every, record, sign=1.0):
    """RK4 with the measurement noise an Ornstein–Uhlenbeck process ν (rate B, stationary variance var, so the low-frequency spectral
    density N₀ = 2var/B), sampled exactly at the steps and held within a step; ν(0) from the stationary law. `record(t, z)` every
    `rec_every` steps. sign = −1 runs the antithetic path −ν of the same draws (the OU recursion is linear). Returns the list of records."""
    z = np.asarray(z0, float).copy()
    n = cl.net.n
    rho = np.exp(-B * dt); sd = sign * np.sqrt(var * (1.0 - rho ** 2))
    nu = sign * np.sqrt(var) * rng.standard_normal(n)
    out = []
    N = int(round(T / dt))
    t = 0.0
    for k in range(N + 1):
        if k % rec_every == 0:
            out.append(record(t, z))
        if k == N:
            break
        f = lambda tt, zz: cl.rhs_noisy(tt, zz, m, nu)
        k1 = f(t, z); k2 = f(t + dt / 2, z + dt / 2 * k1); k3 = f(t + dt / 2, z + dt / 2 * k2); k4 = f(t + dt, z + dt * k3)
        z = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        nu = rho * nu + sd * rng.standard_normal(n)
    return out


def ou_augmented(gen, bvec, D, B, N0):
    """The linearised loop driven by an Ornstein–Uhlenbeck measurement noise (rate B, low-frequency density N₀): the augmented generator of
    (X, ν) and its noise input — for transient_covariance with the initial covariance diag(0, …, 0, N₀B/2)."""
    def Ma(t):
        out = np.zeros((D + 1, D + 1))
        out[:D, :D] = gen(t); out[:D, D] = bvec(t); out[D, D] = -B
        return out
    c = np.zeros(D + 1); c[D] = B
    return Ma, (lambda t: c)


def mean_error(a, b2, ke, B, var, variant='measured'):
    """E[e] at second order, quasi-statically at a point of the orbit (a = ∂_xf, b2 = ∂²_xf there; f̂ = f near the orbit), for the OU
    measurement noise (rate B, variance var): with the compensation f̂(y) ('measured', 'reference': the regressor does not enter E[e]),
    −(b2/(2k_e))var(B − k_e − 2a)/(k_e + B) — the compensation's rectification −½b2·var of the noise and −b2E[eν], its curvature times the
    state error's correlation with the noise (the plant's and the compensation's rectifications of the state error cancel, f̂ = f), rejected
    with the gain k_e; with the compensation at the reference ('reference_both'), b2·Var(e)/(2(k_e − a)), Var(e) = k_e²var/((k_e − a)(k_e − a + B))
    — the plant's rectification of its state jitter alone. It needs ∂²_xf̂ = ∂²_xf on the orbit (an exact representation supplies it). These
    are the noise's own second moments at fixed parameters; the parameters' fluctuation adds terms of the same order (`second_order_mean`)."""
    if variant == 'reference_both':
        c = ke - a
        return 0.5 * b2 * ke ** 2 * var / (c ** 2 * (c + B))
    return -(b2 / (2 * ke)) * var * (B - ke - 2 * a) / (ke + B)


def eiv_force(cv, fprime, fsecond, ke, B, var, N=1024, hx=1e-5, variant='measured', return_g=False):
    """The second-order mean force of the measurement noise on the collective coordinates (real time, whitened). `variant` = 'measured'
    (Proposition D‴.30: the regressor and the compensation at y): the errors in variables E[ε_m ∇Φ̃(y)] at second order,
    ⟨Ĵ′(x_m)ᵀE[ε_m²] + Ĵ(x_m)ᵀE[e]⟩ with the OU noise (rate B, variance var) and, quasi-statically along the orbit,
    E[ε_m²] = var(a² + k_eB)/(k_e(k_e + B)) and E[e] = mean_error(…) (Ĵ′ = ∂_xĴ by central differences); 'reference' (the desired regressor,
    Proposition D‴.31(iv)) and 'reference_both' (and the compensation at the reference): ⟨Ĵ(x_m)ᵀE[e]⟩ alone — a least-squares force, the
    learning force of the target perturbation g = k_eE[e]. In the slow time the force is k_e times this. With `return_g`, also the samples
    of g = k_eE[e] along the orbit (N, at t_k = kT₀/N)."""
    assert variant in ('measured', 'reference', 'reference_both')
    t = np.arange(N) * cv.T0 / N
    F = np.zeros(cv.d)
    g = np.zeros(N)
    for k, s in enumerate(t):
        x = float(np.atleast_1d(cv.xr(s))[0])
        Jx = AF.whitened_jacobian(cv.net, np.array([x]), cv.atoms, cv.w, cv.G)[0]
        a = float(np.atleast_1d(fprime(x))[0]); b2 = float(np.atleast_1d(fsecond(x))[0])
        Ee = mean_error(a, b2, ke, B, var, variant)
        g[k] = ke * Ee
        F += Jx * Ee
        if variant == 'measured':
            Jp = (AF.whitened_jacobian(cv.net, np.array([x + hx]), cv.atoms, cv.w, cv.G)[0]
                  - AF.whitened_jacobian(cv.net, np.array([x - hx]), cv.atoms, cv.w, cv.G)[0]) / (2 * hx)
            F += Jp * var * (a ** 2 + ke * B) / (ke * (ke + B))
    return (F / N, g) if return_g else F / N


def equivalent_perturbation(H, F):
    """FᵀH⁻¹F = min{‖g‖²_ρ : ⟨Ĵᵀg⟩ = F}: the squared norm of the smallest target perturbation whose learning force is F — the function-level
    stationary response to F under the adaptive prior, and a bound on it under every leak S with S + Sᵀ ⪰ 0 (Proposition D‴.31(iii))."""
    return float(F @ np.linalg.solve(H, F))


def response_size(H, F, S):
    """The function-level size ‖Ĵ(H + S)⁻¹F‖²_ρ = Fᵀ(H + S)⁻ᵀH(H + S)⁻¹F of the stationary response to a persistent force F under the
    zero-frequency leak S (a matrix, symmetric or not, or a scalar times I)."""
    S = S * np.eye(H.shape[0]) if np.isscalar(S) else np.asarray(S, float)
    y = np.linalg.solve(H + S, F)
    return float(y @ H @ y)


def mean_response(H, F_slow, sigma0, eta, s):
    """The mean collective deviation ȳ(s) under a constant slow-time force F_slow from ȳ(0) = 0 (and the prior at rest): the fixed prior
    (η = 0) relaxes at H + σ₀I to (H + σ₀I)⁻¹F, the adaptive prior (η = ∞) at H to H⁻¹F; a finite η reaches H⁻¹F as well (its stationary
    point: ζ = y, Hy = F), through the pair dynamics (integrated exactly by the matrix exponential)."""
    from scipy.linalg import expm
    d = H.shape[0]
    A = collective_block(H, sigma0, eta)
    rhs = np.zeros(A.shape[0]); rhs[:d] = F_slow
    Yinf = -np.linalg.solve(A, rhs)
    return (Yinf - expm(A * s) @ Yinf)[:d]


# ---------------------------------------------------------------------- the full second-order mean (the parameters' fluctuation included)
def _whitening(cv, cl):
    """(T, Tinv): the physical state (x, Θ[, Z]) of the K₀-atom loop `cl` from the whitened (δe, y[, ζ]) of `cv`, δϑ_k = √G y_k/√w_k (T is
    Dp × D); under a fixed prior (η = 0) Z is in the state but constant, and the whitened state has no ζ."""
    d = cv.d
    sc = (np.sqrt(np.asarray(cv.G, float))[None, :] / np.sqrt(np.asarray(cv.w, float))[:, None]).ravel()
    Dp = cl.state_dim(cl.K0)
    D = 1 + d if (cl.eta == 0 or np.isinf(cl.eta)) else 1 + 2 * d
    T = np.zeros((Dp, D)); Tinv = np.zeros((D, Dp))
    T[0, 0] = Tinv[0, 0] = 1.0
    T[1:1 + d, 1:1 + d] = np.diag(sc); Tinv[1:1 + d, 1:1 + d] = np.diag(1 / sc)
    if D == 1 + 2 * d:
        T[1 + d:, 1 + d:] = np.diag(sc); Tinv[1 + d:, 1 + d:] = np.diag(1 / sc)
    return T, Tinv


def second_order_mean(cl, cv, fprime, B, var, N=256, h=2e-3, rtol=1e-10, periodic=False):
    """The full second-order small-noise mean of the weighted K₀-atom loop at the lift (D3r-1 of the round-21 review): `cl` a NoisyRealLoop
    with one particle per class (labels 0, …, K₀ − 1, the weights cv.w; its regressor and compensation the variant, its η and k_e the case),
    `cv` the CollapsedVariational at the same atoms (η = ∞ for a fixed prior), an OU measurement noise (rate B, variance var).

    The periodic covariance P(t) of the linearised loop with the noise as a state, mapped to the physical coordinates and contracted with the
    second derivatives of the true vector field (a fourth-order stencil along P's eigenvectors), gives the second-order mean forcing
    f(t) = ½∂²F[P(t)]; it splits exactly into the noise part (P's (δe, ν) block: the noise's own second moments at fixed parameters,
    Proposition D‴.30(iii)'s force) and the parameter part (the blocks with δΘ: the compensation's curvature ½D²_Θf̂[Σ_ΘΘ] and its
    slope–parameter term in the error equation, the correlation ∇²_ϑΦ̃E[ε_mδϑ] in the learning law). Returns the slow-time whitened forces,
    each split into the error equation's part ⟨Ĵᵀg⟩, g = k_ef_e/c (c(t) the error channel's rate), and the learning law's k_e⟨f_y⟩
    (Fe_noise, Fl_noise, Fe_param, Fl_param; F_noise and F_param their sums), the samples of g (g_noise, g_param) and of Ĵ (J) at
    t_k = kT₀/N, the mean function's curvature c̄ = ½D²_Θf̂[Σ̄_ΘΘ] on the orbit's samples x_m(t_k) (Σ̄ the period average, `Sig`, physical
    coordinates), and the relative periodicity error of P; with `periodic`, also the exact periodic second-order mean m(t_k) of the
    linearised loop (whitened) under the full forcing."""
    ke, eta = cl.ke, cl.eta
    d = cv.d
    gen = fixed_generator(cv, ke) if eta == 0 else (lambda t: cv.generator(t, ke))
    if cl.compensation == 'reference':
        gen = reference_compensation(gen, cv, fprime, ke)
    T, Tinv = _whitening(cv, cl)
    Dp, D = T.shape
    b0 = collapsed_noise_input(cv, fprime, ke, compensation=cl.compensation)
    N0 = 2 * var / B
    Ma, cfun = ou_augmented(gen, lambda t: b0(t)[:D], D, B, N0)
    P0 = periodic_covariance(Ma, cfun, cv.T0, D + 1, N0)
    tk = np.arange(N) * cv.T0 / N
    _, Ps = transient_covariance(Ma, cfun, D + 1, N0, cv.T0, np.append(tk, cv.T0), P0=P0, rtol=rtol)
    per_err = float(np.abs(Ps[-1] - P0).max() / np.abs(P0).max())
    Ta = np.zeros((Dp + 1, D + 1)); Ta[:Dp, :D] = T; Ta[Dp, D] = 1.0
    z_lift = cl.lift(cv.atoms)
    noise = np.array([0, Dp])
    w4 = np.array([-1.0, 16.0, -30.0, 16.0, -1.0]) / 12.0

    def field(t, xi):
        return cl.rhs_noisy(t, xi[:Dp], cl.K0, xi[Dp:Dp + 1])

    def contract(t, xi0, Pm, F0):
        lam, U = np.linalg.eigh(0.5 * (Pm + Pm.T))
        out = np.zeros(Dp)
        for l, u in zip(lam, U.T):
            if l != 0.0:
                out += 0.5 * l * sum(c * (F0 if j == 0 else field(t, xi0 + j * h * u)) for c, j in zip(w4, (-2, -1, 0, 1, 2))) / h ** 2
        return out
    fn, fp = np.zeros((N, D)), np.zeros((N, D))
    Jk, ck = np.zeros((N, d)), np.zeros(N)
    for k, t in enumerate(tk):
        Pp = Ta @ Ps[k] @ Ta.T
        Pn = np.zeros_like(Pp); Pn[np.ix_(noise, noise)] = Pp[np.ix_(noise, noise)]
        xi0 = np.append(z_lift, 0.0); xi0[0] = float(np.atleast_1d(cv.xr(t))[0])
        F0 = field(t, xi0)
        fn[k] = Tinv @ contract(t, xi0, Pn, F0)
        fp[k] = Tinv @ contract(t, xi0, Pp - Pn, F0)
        Jk[k] = cv.J(t)[0]
        ck[k] = -gen(t)[0, 0]
    part_e = lambda f: np.mean(ke * Jk * (f[:, 0] / ck)[:, None], axis=0)
    part_l = lambda f: np.mean(ke * f[:, 1:1 + d], axis=0)
    Sig = np.array([(Ta @ P @ Ta.T)[1:1 + d, 1:1 + d] for P in Ps[:N]]).mean(axis=0)
    lam, U = np.linalg.eigh(0.5 * (Sig + Sig.T))
    Th0 = np.asarray(cv.atoms, float)
    xs = np.array([float(np.atleast_1d(cv.xr(t))[0]) for t in tk])
    fh = lambda Th: np.array([cl.f_hat(np.array([x]), Th)[0] for x in xs])
    f00 = fh(Th0)
    c_bar = np.zeros(N)
    for l, u in zip(lam, U.T):
        c_bar += 0.5 * l * sum(c * (f00 if j == 0 else fh(Th0 + j * h * u.reshape(Th0.shape))) for c, j in zip(w4, (-2, -1, 0, 1, 2))) / h ** 2
    out = dict(Fe_noise=part_e(fn), Fl_noise=part_l(fn), Fe_param=part_e(fp), Fl_param=part_l(fp), g_noise=ke * fn[:, 0] / ck,
               g_param=ke * fp[:, 0] / ck, c_bar=c_bar, J=Jk, t=tk, per_err=per_err, Sig=Sig)
    out['F_noise'] = out['Fe_noise'] + out['Fl_noise']
    out['F_param'] = out['Fe_param'] + out['Fl_param']
    if periodic:
        fhat = np.fft.fft(fn + fp, axis=0) / N
        om = np.fft.fftfreq(N, d=1.0 / N) * (2 * np.pi / cv.T0)
        force = lambda t: np.real(np.exp(1j * om * t) @ fhat)
        rhs = lambda t, y: np.concatenate([(gen(t) @ y[:D * D].reshape(D, D)).ravel(), gen(t) @ y[D * D:] + force(t)])
        sol = solve_ivp(rhs, (0.0, cv.T0), np.concatenate([np.eye(D).ravel(), np.zeros(D)]), method='DOP853', rtol=1e-12, atol=1e-16)
        Phi, pT = sol.y[:D * D, -1].reshape(D, D), sol.y[D * D:, -1]
        m0 = np.linalg.solve(np.eye(D) - Phi, pT)
        sol2 = solve_ivp(lambda t, y: gen(t) @ y + force(t), (0.0, cv.T0), m0, method='DOP853', rtol=1e-12, atol=1e-16,
                         t_eval=np.append(tk, cv.T0))
        out['m'] = sol2.y.T[:N]
        out['m_per_err'] = float(np.abs(sol2.y[:, -1] - m0).max() / np.abs(m0).max())
    return out


def mean_biases(H, out, S=None):
    """From `second_order_mean`'s output for a prior with the zero-frequency leak S (None: an adaptive prior, every η > 0; σ₀: the fixed prior
    — `out` computed with that prior, whose covariance sets the parameters' part): the stationary response ȳ = (H + S)⁻¹F to the noise part
    and to the full force, and the function-level biases — of the mean parameters, ‖Ĵȳ‖²_ρ = ȳᵀHȳ, and of the mean function, ‖Ĵȳ + c̄‖²_ρ
    (the curvature c̄ = ½D²_Θf̂[Σ̄_ΘΘ] of the network under the parameters' fluctuation shifts the mean function; the loop removes its tangent
    part) — the least-squares bounds ‖g‖²_ρ of the error equation's part (noise, full), and the responses to the parameters' part alone and
    to its learning-law part alone. Returns a dict."""
    Rm = np.linalg.inv(H if S is None else H + (S * np.eye(H.shape[0]) if np.isscalar(S) else S))
    yn, yf = Rm @ out['F_noise'], Rm @ (out['F_noise'] + out['F_param'])
    J = out['J']
    yp, yl = Rm @ out['F_param'], Rm @ out['Fl_param']
    return dict(noise=float(yn @ H @ yn), full=float(yf @ H @ yf), function=float(np.mean((J @ yf + out['c_bar']) ** 2)),
                g2_noise=float(np.mean(out['g_noise'] ** 2)), g2_full=float(np.mean((out['g_noise'] + out['g_param']) ** 2)),
                param=float(yp @ H @ yp), law=float(yl @ H @ yl), c2=float(np.mean(out['c_bar'] ** 2)))


def separation_drift(cv, fprime, fsecond, ke, B, var, N=512, hx=1e-5, hth=1e-5, variant='measured'):
    """S_k (K₀, p, p), real time, per class k: the averaged second-order drift of a separation of class k under the measured regressor at the
    lift, ⟨E[ε_m²]G^{1/2}∂_x∇²_ϑΦ̃(x_m; ϑ_k)G^{1/2} + E[e]G^{1/2}∇²_ϑΦ̃(x_m; ϑ_k)G^{1/2}⟩ with the quasi-static moments of `eiv_force` (the
    noise's own, at fixed parameters): a separation u = G^{−1/2}δϑ_i of class k obeys u̇ = (G^{1/2}ε_m∇²_ϑΦ̃(y; ϑ_k)G^{1/2} − σ₀/k_e)u along
    the noisy class trajectory, so the replicated configuration attracts at the start while σ₀/k_e exceeds the largest eigenvalue of S_k
    (D3q-3): invariance under every path, attraction only below a noise level. Along the trajectory E[e] also carries the mean tracking
    error of the developing errors-in-variables deviation (round 22, ID-1): `separation_drift_series` with `second_order_transient`.
    `variant` 'reference' (the desired regressor) or 'reference_both' (and the compensation at x_m): the law's Jacobian is then
    ε_m∇²_ϑΦ̃(x_m; ϑ_k), with no ∂_x term, so S_k = ⟨E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}⟩ alone, E[e] the variant's (`mean_error`)."""
    assert variant in ('measured', 'reference', 'reference_both')
    t = np.arange(N) * cv.T0 / N
    p = cv.atoms.shape[1]
    g = np.sqrt(np.asarray(cv.G, float))
    S = np.zeros((len(cv.atoms), p, p))

    def hess(x, th):
        Hm = np.zeros((p, p))
        for j in range(p):
            dp = np.zeros(p); dp[j] = hth
            Hm[:, j] = (cv.net.grad_g(np.array([x]), np.array([1.0]), (th + dp)[None, :])[0]
                        - cv.net.grad_g(np.array([x]), np.array([1.0]), (th - dp)[None, :])[0]) / (2 * hth)
        return 0.5 * (Hm + Hm.T)
    for s in t:
        x = float(np.atleast_1d(cv.xr(s))[0])
        a = float(np.atleast_1d(fprime(x))[0]); b2 = float(np.atleast_1d(fsecond(x))[0])
        Ee = mean_error(a, b2, ke, B, var, variant)
        Eem2 = var * (a ** 2 + ke * B) / (ke * (ke + B)) if variant == 'measured' else 0.0
        for k, th in enumerate(cv.atoms):
            Hx = (hess(x + hx, th) - hess(x - hx, th)) / (2 * hx) if Eem2 else 0.0
            S[k] += g[:, None] * (Eem2 * Hx + Ee * hess(x, th)) * g[None, :]
    return S / N


def second_order_transient(cl, cv, fprime, B, var, t_end, dt=0.05, h=2e-3, rtol=1e-9):
    """The second-order small-noise transient of the weighted K₀-atom loop `cl` (one particle per class, as in `second_order_mean`) from its
    lift, the OU noise (rate B, variance var) stationary at t = 0 and the state's mean and covariance otherwise zero (round 22, ID-1): the
    covariance P(t) of the augmented linearised loop, the second-order forcing ½∂²F[P(t)] of the true field (the stencil of
    `second_order_mean`) on the grid t_k = k·dt, and the mean m(t) of ṁ = M(t)m + force(t), m(0) = 0 (RK4 on the grid, the forcing linear
    in between). Returns dict(t, m (whitened (δe, y[, ζ])), Ee = E[e](t) = m_δe — the noise's own mean error plus the mean tracking error of
    the developing errors-in-variables deviation —, Eem2 = E[ε_m²](t) = E[(δe + ν)²] at second order)."""
    ke, eta = cl.ke, cl.eta
    gen = fixed_generator(cv, ke) if eta == 0 else (lambda t: cv.generator(t, ke))
    if cl.compensation == 'reference':
        gen = reference_compensation(gen, cv, fprime, ke)
    T, Tinv = _whitening(cv, cl)
    Dp, D = T.shape
    b0 = collapsed_noise_input(cv, fprime, ke, compensation=cl.compensation)
    N0 = 2 * var / B
    Ma, cfun = ou_augmented(gen, lambda t: b0(t)[:D], D, B, N0)
    n_t = int(round(t_end / dt))
    tk = np.arange(n_t + 1) * dt
    P0 = np.zeros((D + 1, D + 1)); P0[D, D] = var
    _, Ps = transient_covariance(Ma, cfun, D + 1, N0, t_end, tk, P0=P0, rtol=rtol)
    Ta = np.zeros((Dp + 1, D + 1)); Ta[:Dp, :D] = T; Ta[Dp, D] = 1.0
    z_lift = cl.lift(cv.atoms)
    w4 = np.array([-1.0, 16.0, -30.0, 16.0, -1.0]) / 12.0
    field = lambda t, xi: cl.rhs_noisy(t, xi[:Dp], cl.K0, xi[Dp:Dp + 1])
    force = np.zeros((len(tk), D))
    for k, t in enumerate(tk):
        Pp = Ta @ Ps[k] @ Ta.T
        xi0 = np.append(z_lift, 0.0); xi0[0] = float(np.atleast_1d(cv.xr(t))[0])
        F0 = field(t, xi0)
        lam, U = np.linalg.eigh(0.5 * (Pp + Pp.T))
        out = np.zeros(Dp)
        for l, u in zip(lam, U.T):
            if l > 0.0:
                out += 0.5 * l * sum(c * (F0 if j == 0 else field(t, xi0 + j * h * u)) for c, j in zip(w4, (-2, -1, 0, 1, 2))) / h ** 2
        force[k] = Tinv @ out
    m = np.zeros((len(tk), D))
    for k in range(n_t):
        t = tk[k]; fm = 0.5 * (force[k] + force[k + 1])
        g1 = gen(t) @ m[k] + force[k]
        g2 = gen(t + dt / 2) @ (m[k] + dt / 2 * g1) + fm
        g3 = gen(t + dt / 2) @ (m[k] + dt / 2 * g2) + fm
        g4 = gen(t + dt) @ (m[k] + dt * g3) + force[k + 1]
        m[k + 1] = m[k] + dt / 6 * (g1 + 2 * g2 + 2 * g3 + g4)
    Eem2 = np.array([P[0, 0] + 2 * P[0, D] + P[D, D] for P in Ps])
    return dict(t=tk, m=m, Ee=m[:, 0].copy(), Eem2=Eem2)


def separation_drift_series(cv, ts, Ee, Eem2, hx=1e-5, hth=1e-5):
    """S_k(t) (len(ts), K₀, p, p), real time: the separations' second-order drift E[ε_m²](t)G^{1/2}∂_x∇²_ϑΦ̃(x_m(t); ϑ_k)G^{1/2} +
    E[e](t)G^{1/2}∇²_ϑΦ̃(x_m(t); ϑ_k)G^{1/2} with the moments given along a trajectory (`second_order_transient`; `separation_drift` is its
    period average with the noise's own quasi-static moments at the lift). The parameters' cross term Σ_lE[ε_mδϑ_kl]∂_{ϑ_l}(·) is omitted
    (on the instance below 10⁻³ of the drift, referee round 22)."""
    p = cv.atoms.shape[1]
    g = np.sqrt(np.asarray(cv.G, float))

    def hess(x, th):
        Hm = np.zeros((p, p))
        for j in range(p):
            dp = np.zeros(p); dp[j] = hth
            Hm[:, j] = (cv.net.grad_g(np.array([x]), np.array([1.0]), (th + dp)[None, :])[0]
                        - cv.net.grad_g(np.array([x]), np.array([1.0]), (th - dp)[None, :])[0]) / (2 * hth)
        return 0.5 * (Hm + Hm.T)
    out = np.zeros((len(ts), len(cv.atoms), p, p))
    for i, s in enumerate(ts):
        x = float(np.atleast_1d(cv.xr(s))[0])
        for k, th in enumerate(cv.atoms):
            Hx = (hess(x + hx, th) - hess(x - hx, th)) / (2 * hx)
            out[i, k] = g[:, None] * (Eem2[i] * Hx + Ee[i] * hess(x, th)) * g[None, :]
    return out


# ---------------------------------------------------------------------- Lemma D‴.26a(v): the slow noise, decoupled exactly
def noise_decoupling(cv, fprime, eps, N=512, tol=1e-18, maxit=500, with_b1=False):
    """Chang's second transformation for the linearised loop with white measurement noise (Lemma D‴.26a(v)): S(t; ε), the T₀-periodic
    solution of the anti-causal S − εṠ = εKᵀ − ε²(SPKᵀ + AS) (P, A from `CollapsedVariational.chang`), with which Y′ = Y + εSv (v = w − PY,
    δe = εw) obeys dY′ = εAY′dt + εb√N₀dW exactly — the drift decoupled from v, the noise shared with it —, b = ε⁻¹[Kᵀ − ε⁻¹(1 + εa)S −
    εSPKᵀ]; and b₀ = −(K̇ + aK)ᵀ, its leading term (the next, εb₁ with b₁ = −(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ, with with_b1=True: referee round 22,
    D3s-5). Returns (P, S, b, b₀ (N, q, n) each, the coupling residual max|εKᵀ + εṠ − S − ε²(SPKᵀ + AS)|[, b₁])."""
    t, K, Lam, om = cv.chang_grid(N)
    P, A, _, _ = cv.chang(eps, N)
    Kt = np.transpose(K, (0, 2, 1))
    PKt = np.einsum('tij,tkj->tik', P, K)
    a = np.array([float(np.atleast_1d(fprime(np.atleast_1d(cv.xr(s))[0]))[0]) for s in t])
    S = eps * Kt.copy()
    for it in range(maxit):
        rhs = eps * Kt - eps ** 2 * (np.einsum('tij,tjk->tik', S, PKt) + np.einsum('tij,tjk->tik', A, S))
        Sn = np.real(np.fft.ifft(np.fft.fft(rhs, axis=0) / (1 - 1j * eps * om)[:, None, None], axis=0))
        dS = np.abs(Sn - S).max(); S = Sn
        if dS < tol:
            break
    d = lambda X: np.real(np.fft.ifft((1j * om).reshape(-1, 1, 1) * np.fft.fft(X, axis=0), axis=0))
    coupling = float(np.abs(eps * Kt + eps * d(S) - S - eps ** 2 * (np.einsum('tij,tjk->tik', S, PKt) + np.einsum('tij,tjk->tik', A, S))).max())
    b = (Kt - (S / eps) * (1 + eps * a)[:, None, None] - eps * np.einsum('tij,tjk->tik', S, PKt)) / eps
    b0 = -np.transpose(d(K) + a[:, None, None] * K, (0, 2, 1))
    if with_b1:
        KKtK = np.einsum('tij,tkj,tkl->til', K, K, K)
        b1 = -np.transpose(d(d(K)) + KKtK + np.einsum('tij,kj->tik', K, Lam) + a[:, None, None] * d(K), (0, 2, 1))
        return P, S, b, b0, coupling, b1
    return P, S, b, b0, coupling


def decoupled_covariance(cv, fprime, ke, N=512, period_average=False, n_avg=256):
    """The exact periodic stationary covariance of the linearised loop under white measurement noise (N₀ = 1) at t = 0, in the coordinates
    (v, Y′) of Lemma D‴.26a(v), against the lemma's leading terms: returns dict(Vv = Cov(v), C = E[vY′ᵀ], CY′ = Cov(Y′), CY = Cov(Y), the
    averaged slow covariance Σ̄ (Ā₀Σ̄ + Σ̄Ā₀ᵀ + εQ̄ = 0), the fast jitter (ε/2)K(0)ᵀK(0), the next term of Cov(Y) − Cov(Y′),
    −(ε²/2)(d(KᵀK)/dt + 2aKᵀK)(0) — the cross term −ε(SC + CᵀSᵀ) ≈ −ε²(d(KᵀK)/dt + 2aKᵀK) plus the jitter's own next term, half of it with
    the opposite sign (referee round 22, D3s-2) —, b(0)); with period_average=True also the
    period average, the minimum and the maximum over n_avg times of the floor tr(H Cov y(t)) and of its prediction tr(H(Σ̄ + (ε/2)Ĵ(t)ᵀĴ(t)))
    (the lemma's floor is a period average: at time t the fast part is (ε/2)tr(HĴ(t)ᵀĴ(t)), an oscillation of the floor's own order). At
    η = 0 (the fixed prior) K = Ĵ, Y = y, Λ = σ₀I."""
    eps = 1.0 / ke
    n = cv.n
    P, S, b, b0, _ = noise_decoupling(cv, fprime, eps, N)
    D = cv.dim
    Pc = periodic_covariance(lambda t: cv.generator(t, ke), collapsed_noise_input(cv, fprime, ke), cv.T0, D, 1.0)
    Tm = np.zeros((D, D))
    Tm[:n, :n] = np.eye(n) / eps; Tm[:n, n:] = -P[0]
    Tm[n:, :n] = S[0]; Tm[n:, n:] = np.eye(D - n) - eps * S[0] @ P[0]
    Cn = Tm @ Pc @ Tm.T
    H, Q = noise_matrices(cv, fprime)
    A0 = collective_block(H, cv.sigma0, cv.eta)
    Qb = np.zeros_like(A0); Qb[:cv.d, :cv.d] = Q
    tg, Kg, Lam, om = cv.chang_grid(N)
    K0 = Kg[0]
    Kd0 = np.real(np.fft.ifft((1j * om)[:, None, None] * np.fft.fft(Kg, axis=0), axis=0))[0]
    a0 = float(np.atleast_1d(fprime(np.atleast_1d(cv.xr(0.0))[0]))[0])
    Sigma = solve_continuous_lyapunov(A0, -eps * Qb)
    out = dict(Vv=Cn[:n, :n], C=Cn[:n, n:], CYp=Cn[n:, n:], CY=Pc[n:, n:], Sigma=Sigma, fast=(eps / 2) * K0.T @ K0,
               next=-0.5 * eps ** 2 * (Kd0.T @ K0 + K0.T @ Kd0 + 2 * a0 * K0.T @ K0), b=b[0], b0=b0[0])
    if period_average:
        d = cv.d
        gen = lambda t: cv.generator(t, ke)
        bv = collapsed_noise_input(cv, fprime, ke)

        def rhs(t, p):
            Pm = p.reshape(D, D); M = gen(t); bb = bv(t)
            return (M @ Pm + Pm @ M.T + np.outer(bb, bb)).ravel()
        ts = np.arange(n_avg) * cv.T0 / n_avg
        sol = solve_ivp(rhs, (0.0, cv.T0), Pc.ravel(), method='DOP853', t_eval=ts, rtol=1e-11, atol=1e-16)
        Hd = H
        fl = np.array([np.trace(Hd @ sol.y[:, i].reshape(D, D)[n:n + d, n:n + d]) for i in range(len(ts))])
        pr = np.array([np.trace(Hd @ (Sigma[:d, :d] + (eps / 2) * cv.J(t).T @ cv.J(t))) for t in ts])
        out.update(floor_avg=float(fl.mean()), floor_min=float(fl.min()), floor_max=float(fl.max()),
                   pred_avg=float(pr.mean()), pred_min=float(pr.min()), pred_max=float(pr.max()))
    return out

