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
noise increment, so the separations are not excited: the noise is coherent and its floor does not average out with m. This module gives
Q and H = ⟨ĴᵀĴ⟩, the stationary covariance of the averaged collective block (the Lyapunov equation), the closed-form function-level floors
(N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) in the eigenbasis of H, and the exact periodic stationary covariance of the linearised loop.
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
    HR + B = I can be checked rather than assumed."""
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


def simulate_noisy(cl, z0, m, T, dt, B, var, rng, rec_every, record):
    """RK4 with the measurement noise an Ornstein–Uhlenbeck process ν (rate B, stationary variance var, so the low-frequency spectral
    density N₀ = 2var/B), sampled exactly at the steps and held within a step; ν(0) from the stationary law. `record(t, z)` every
    `rec_every` steps. Returns the list of records."""
    z = np.asarray(z0, float).copy()
    n = cl.net.n
    rho = np.exp(-B * dt); sd = np.sqrt(var * (1.0 - rho ** 2))
    nu = np.sqrt(var) * rng.standard_normal(n)
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
    −(b2/(2k_e))var(B − k_e − 2a)/(k_e + B) — the controller's rectification ½b2·var of the noise, rejected with the gain k_e, and the plant's
    of its own jitter; with the compensation at the reference ('reference_both'), b2·Var(e)/(2(k_e − a)), Var(e) = k_e²var/((k_e − a)(k_e − a + B))
    — the plant's rectification of its state jitter alone."""
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
    stationary response to F under the adaptive prior, and a bound on it under every leak S ≥ 0 (Proposition D‴.31(iii))."""
    return float(F @ np.linalg.solve(H, F))


def response_size(H, F, S):
    """The function-level size ‖Ĵ(H + S)⁻¹F‖²_ρ = Fᵀ(H + S)⁻¹H(H + S)⁻¹F of the stationary response to a persistent force F under the
    zero-frequency leak S (a matrix, or a scalar times I)."""
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
