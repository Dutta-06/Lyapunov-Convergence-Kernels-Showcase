"""
Proposition D‴.41: measurement noise on the gauge — no selection in continuous time; per update, the increment's second moment selects
=====================================================================================================================================

Proposition D‴.36's instance (three bumps, one particle per class, m = 3). The measured error r + ξ drives both laws (the positions and
the masses), ξ independent over the quadrature nodes, of intensity v (⟨ξh⟩_ρ has the variance v⟨h²⟩_ρ per unit slow time): the noise
enters along the descent's own fields U_x = (β_iG∇_ϑΦ̃(x; ϑ_i), γ(Φ̃(x; ϑ_i) − f̂(x))) (`label_noise.fields`), the drift being ⟨rU⟩_ρ.

A. The identity at three points of M and γ = 10⁻², 10⁻⁵ (Katzenberger's limit with the covariance Σ = v⟨UUᵀ⟩_ρ at the point, in the flat
   coordinates of exp_gauge_drift.py): the Itô drift ½∂²Φ[Σ], the Wong–Zakai term (v/2)∂Φ⟨(U·∇)U⟩_ρ, their sum (the Stratonovich drift,
   zero by the proposition), the diffusion along M (zero), and the closed form on the masses in the slow-replicator limit,
   ū_k = −(v/2)β_kG_α⟨Φ_k²⟩_ρ(1 − (a_k/a_max)⁴), against the Itô drift; and the same sum for heteroscedastic noise (the variance along the
   orbit modulated by 1 + 0.8 sin t). A2: the noise's structure over the nodes at W₀ and 1/3, γ = 10⁻¹ … 10⁻⁵ — heteroscedastic profiles
   (1 + 0.8 sin t, 1 + 0.8 cos 2t), a measurement offset (the same ξ at every node), a noise correlated along the orbit (a real-time
   Ornstein–Uhlenbeck noise of correlation time τ_n, C_ij = Σ_pR(t_i − t_j + pT), its one-sided kernel's Lévy-area part included): the
   Stratonovich drift over its own Itô drift (O(γ): in the slow-replicator limit a second mechanism cancels it); and the coloured-noise
   identity at first order in the correlation time, ⟨∇ψ·AU⟩_ρ = −½P⟨∂²F[U, U]⟩_ρ + P⟨(AU·∇)U⟩_ρ = 0 (A = ∂F on M).
B. The selection's target: P_k = G_α⟨Φ_k²⟩_ρ, the minimiser W_ℓ of Ψ_ℓ = (v/4)Σ_kP_k(β_ka′_k)² on the simplex (against Proposition
   D‴.38's W♯), and the closed form's replicator flow linearised there; the zero of the Itô drift at γ = 10⁻² (W_ℓ is the slow limit's);
   the update's coordinates: an Euler step in the amplitude a instead of α (Katzenberger's drift of the Itô SDE in those coordinates at
   γ = 10⁻⁵, against the closed form with (1 − q)(1 + 3q) for 1 − q², and its target).
C. The full loop with measurement noise only (τ = 0, γ = 10⁻², v = V_C, from W₀ = (0.40, 0.25, 0.35), SEEDS_C × S_C slow units): the noise
   drawn afresh at each step (Euler, at the steps 0.25 and 0.125), AR(1) noise at ρ = −1/2, −1/3, 1/3, 0.6, 0.8, 0.9, 0.95 (the same
   intensity), and Heun's scheme (Stratonovich) with fresh and with AR(1) (ρ = 0.9) draws. The gauge coordinate y(s) = Φ(z(s)) (the limit
   map of the noiseless flow from the state, every REC_PHI slow units): its least-squares slope and its end displacement, against
   (1 − ρ)/(1 + ρ) times part A's Itô drift at W₀ and against the linear-response prediction — Katzenberger's drift of the filtered
   covariance Σ_f (`label_noise.filtered_covariance`: the increment's second moment, the cross term h²∂²Φ[F, ⟨Uξ⟩] included) times the
   per-step variance's factor; zero for Heun. (The raw masses' window means carry the O(v) shift of the normal fluctuations' mean, which
   builds up from the start on M.)
D. The example: the scalar factorisation f(u, v) = uv fitted to y = 1 with label noise, u̇ = (ξ − r)v, v̇ = (ξ − r)u, r = uv − y. The
   rescaling symmetry's Noether charge Q = u² − v² is conserved by the gradient flow and, pathwise, by every noise field (v, u) — so by the
   Stratonovich (Heun) dynamics —; an Euler step gives exactly Q_{n+1} = Q_n(1 − (ξ_n − r_n)²h²) whatever the noise sequence, so E[Q]
   decays at h·E[(ξ − r)²] per unit time: σ² for noise drawn afresh (the Itô SDE's exact rate), (1 − ρ)/(1 + ρ)σ² times the linear-
   response factor φ(hλ) (λ = u² + v² = √(Q² + 4) on uv = 1; `label_noise.lr_factor`) for AR(1) noise at the same intensity — the
   residual r follows a slow noise and absorbs it (φ → 1/(1 + λτ_c)) —; N_PATHS_D paths, the decay rate −log(E[Q(T)]/Q(0))/(σ²T) against
   the per-step factor and against the linear-response prediction (dQ/dt = −σ²(1 − ρ)/(1 + ρ)φ(hλ(Q))Q integrated).
E. A diagonal linear network, w = u⊙u − v⊙v ∈ ℝ^d fitted to n < d samples (X standard normal, w* sparse), from u = v = α to the
   gradient flow's interpolant, then with label noise of intensity V_E (the residual r_j + ξ_j, ξ independent over the samples). Each
   product P_i = u_iv_i — the effective initialisation scale, which sets how close the interpolant is to the ℓ¹-minimal one — is conserved
   by the gradient flow and by every noise field, so by the Stratonovich dynamics; an Euler step multiplies it by exactly
   1 − 4h²(g_i − ξ̃_i)² (g = ∇_wL, ξ̃ = Xᵀξ/n), so E[P_i] decays at 4h²E[(X_iᵀe)²]/n² per step, e = r − ξ: 4V_Ec_i/n (c_i = Σ_jX_ji²)
   times the per-step factor for white noise, and for correlated noise times Σ_jφ(hκ_j)(X_iᵀq_j)²/c_i over the residual's modes (κ_j, q_j)
   of K = (4/n)X diag(u² + v²)Xᵀ — coordinate by coordinate: the colour changes the selection's direction —; the measured rates against
   the linear-response prediction along the run (K from the paths' mean every T_E/25), and the interpolant the noiseless gradient flow
   reaches from the end state (its ℓ¹ norm, its distance from w*, against the start's).
F. A generic overparametrised network: f(x) = Σ_h a_h tanh(w_hx + b_h), H_F units, fitted to N_F points of sin(πx) on [−1, 1], started on
   the interpolating manifold Γ (a Gauss–Newton projection of a random draw); label noise of intensity V_F; the sharpness
   T = (1/n)Σ_j‖∇f(x_j)‖² (tr of the Gauss–Newton matrix) of the projection onto Γ along the run: its initial slope (a quadratic fit)
   against the leading-order Itô prediction −(V_F/4)|P_T∇T|² (P_T the orthogonal projection onto TΓ, ∇T by central differences) times the
   per-step factor, and against the linear-response prediction −(V_F/4)(1 − ρ)/(1 + ρ)∇T·P_T∇trΨ(K) (K = JᵀJ/n, Ψ′ = φ: label noise's
   potential tr Ψ(K), the trace for white noise, (1/τ_c)log det(I + τ_cK) for a slow one); no exact identity here.
G. The noise's colour moves the target. (G1) The zero of the per-update selection on the simplex — Katzenberger's drift of the filtered
   covariance Σ_f at the step DT for AR(1) draws of correlation time τ_c — continued from τ_c = 0 (W_ℓ's finite-γ version) to τ_c = 2·10⁴
   at γ = 10⁻² and 10⁻⁵, with its stability. (G2) The slow limit: for label noise Σ = −mAM on M (A = ∂F, M the mobility), so
   φ(A) → (τ_cA)⁻¹ turns τ_cΣ_f into −A_DΣQᵀ = m·QMQᵀ — the flow's own mobility on the normal space, a thermal noise in the flow's metric —:
   at three points, that limit's residual and symmetry, and the direction of the slow noise's selection against the dither's (Σ = 2τG,
   Proposition D‴.37) and against the fresh noise's (the cosines; the ratio of the dither's drift to its normal part's).

Output: results/chaos/gauge_label.md and gauge_label_raw.npz (GL_WORKERS; `--from-raw` regenerates the md; `--only A,B,C,D,E,F,G`, A
includes A2).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mass_transfer as MT
from src.core import gauge_drift as GDR
from src.core import label_noise as LN
from src.core.runcache import cached
import exp_cohesion_leak as E
import exp_gauge_dither as GD
import exp_gauge_drift as X
import exp_gauge_fixman as Q

OUT = os.environ.get('GL_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GL_WORKERS', 4))
A_MAX, DT = 4.0, GD.DT
W_A = ((1 / 3, 1 / 3), (0.30, 0.35), (0.40, 0.25))
GAMMAS_A = (1e-2, 1e-5)
GAMMA_C, V_C, W0_C, S_C, REC_S, REC_PHI = 1e-2, 5e-5, (0.40, 0.25, 0.35), 2e5, 1000.0, 2e4   # slow units between records, limit maps
SEEDS_C = tuple(range(800, 804))
SCHEMES_C = (('ito', 0.0, DT), ('ito', 0.0, DT / 2), ('ar1', -0.5, DT), ('ar1', -1 / 3, DT), ('ar1', 1 / 3, DT), ('ar1', 0.6, DT),
             ('ar1', 0.8, DT), ('ar1', 0.9, DT), ('ar1', 0.95, DT), ('heun', 0.0, DT), ('heun', 0.9, DT))
GAMMAS_A2, W_A2, TAUN_A2 = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5), ((0.40, 0.25), (1 / 3, 1 / 3)), (0.016, 0.16, 1.6)   # τ_n over the period
SIG2_D, DT_D, T_D, N_PATHS_D, UV0_D, SEED_D = 0.01, 1e-3, 100.0, 400, (2.0, 0.5), 900
SCHEMES_D = (('ito', 0.0), ('ar1', -1 / 3), ('ar1', 1 / 3), ('ar1', 0.6), ('block', 4), ('ar1', 0.9), ('ar1', 0.99), ('ar1', 0.999),
             ('heun', 0.0), ('heun', 0.99))                                   # 'block': redrawn every k steps
D_E, N_E, ALPHA_E, V_E, DT_E, T_E, P_E, SEED_E = 40, 10, 0.5, 0.001, 1e-3, 250.0, 8, 950
H_F, N_F, V_F, DT_F, T_F, P_F, SEED_F, SEED_NF = 20, 8, 1e-3, 2e-3, 200.0, 128, 15, 3
L_MAP = 5e4                                                                    # the limit map's horizon (the slowest normal rate ≈ 3·10⁻⁴: e⁻¹⁵)


def beta_of(z):
    e = np.exp(np.r_[z, 0.0]); return 3 * e / e.sum()


@cached
def a_job(args):
    """Part A at one point and γ: the Itô drift, the Wong–Zakai term, the diffusion, the closed form (v = 1)."""
    gamma, w = args
    F, chart, coords, S0 = X.instance(gamma)
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)
    w = np.asarray(w, float)
    x0 = chart(w)
    U0 = LN.flat_fields(loop, x0, beta_of)
    Sig = U0.T @ U0 / U0.shape[0]
    ito, D, lam = GDR.limit_drift(F, chart, coords, w, Sig, h=X.H_ID, h2=X.H2_ID)
    A = GDR._jac(F, x0, X.H_ID); Jg = GDR._jac(chart, w, X.H_ID)               # the projection ∂Φ = P as in gauge_drift.limit_drift
    ev, V = np.linalg.eig(A); o = np.argsort(np.abs(ev)); V = np.c_[Jg, V[:, o][:, 2:].real]
    P = V[:, :2] @ np.linalg.inv(V)[:2, :]
    Jw = GDR._jac(coords, x0, X.H_ID)
    wz = Jw @ P @ (0.5 * LN.wong_zakai(loop, x0, beta_of))
    W = np.r_[w, 1 - w.sum()]
    u = LN.selection(loop, E.exact_rep(net, W), 3 * W)
    closed = gamma * W * (u - W @ u)
    wt = 1 + 0.8 * np.sin(2 * np.pi * np.arange(U0.shape[0]) / U0.shape[0]); wt = wt / wt.mean()   # heteroscedastic noise
    ito_h = GDR.limit_drift(F, chart, coords, w, (U0 * wt[:, None]).T @ U0 / U0.shape[0], h=X.H_ID, h2=X.H2_ID)[0]
    wz_h = Jw @ P @ (0.5 * LN.wong_zakai(loop, x0, beta_of, weights=wt))
    return dict(gamma=gamma, w=w, ito=np.r_[ito, -ito.sum()], wz=np.r_[wz, -wz.sum()], D=D, lam=lam, closed=closed,
                fit=float(np.abs(loop.residual(E.exact_rep(net, W), 3 * W)).max()), ito_h=np.r_[ito_h, -ito_h.sum()], wz_h=np.r_[wz_h, -wz_h.sum()])


def node_structures(N):
    """Part A2's noise structures over the N nodes (t_i = 2πi/N, the orbit's period T = 2π): name → the one-sided kernel K (K_ij: the noise at
    node j before node i; the covariance C = K + Kᵀ, normalised to tr C = N, as white noise's)."""
    t = 2 * np.pi * np.arange(N) / N
    out = {'white (independent, homoscedastic)': np.eye(N) / 2}
    for name, wt in (('1 + 0.8 sin t', 1 + 0.8 * np.sin(t)), ('1 + 0.8 cos 2t', 1 + 0.8 * np.cos(2 * t))):
        out[name] = np.diag(wt / wt.mean()) / 2
    out['an offset (one ξ at every node)'] = np.ones((N, N)) / 2
    for tn in TAUN_A2:
        e1 = np.exp(-2 * np.pi / (tn * 2 * np.pi))
        K = np.exp(-((t[:, None] - t[None, :]) % (2 * np.pi)) / (tn * 2 * np.pi)) / (1 - e1)      # Σ_{p≥0}R(t_i − t_j + pT), R(u) = e^{−u/τ_n}
        np.fill_diagonal(K, 0.5 * (1 + e1) / (1 - e1))                                            # half the lag-0 term
        out[f'correlated along the orbit, τ_n/T = {tn:g}'] = K * N / np.trace(K + K.T)
    return out


@cached
def a2_job(args):
    """Part A2 at one point and γ: per structure of `node_structures`, the Itô drift ½∂²Φ[Σ_C], the Wong–Zakai term of C and the one-sided
    kernel's Lévy-area part (W_1, W_2, W_3, per unit v); the coloured-noise identity at first order, −½P⟨∂²F[U, U]⟩ + P⟨(AU·∇)U⟩, and the
    size of its second term."""
    gamma, w = args
    F, chart, coords, S0 = X.instance(gamma)
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)
    w = np.asarray(w, float); x0 = chart(w)
    U0 = LN.flat_fields(loop, x0, beta_of); N = U0.shape[0]
    A = GDR._jac(F, x0, X.H_ID)
    P = GDR.spectral_projection(A, GDR._jac(chart, w, X.H_ID))[0]
    Jw = GDR._jac(coords, x0, X.H_ID)
    Dd = LN.field_derivatives(loop, x0, beta_of)
    full = lambda d: np.r_[d, -d.sum()]
    rows = {}
    for name, K in node_structures(N).items():
        C = K + K.T
        ito = GDR.limit_drift(F, chart, coords, w, U0.T @ C @ U0 / N, h=X.H_ID, h2=X.H2_ID)[0]
        rows[name] = dict(ito=full(ito), sym=full(Jw @ P @ (0.5 * np.einsum('ij,jik->k', C, Dd) / N)),
                          lev=full(Jw @ P @ (0.5 * np.einsum('ij,jik->k', K - K.T, Dd) / N)))
    t1 = np.mean([GDR._d2(F, x0, U0[y], U0[y], X.H2_ID) for y in range(N)], axis=0)
    hh = 1e-4
    t2 = np.mean([(LN.flat_fields(loop, x0 + hh * (A @ U0[y]), beta_of)[y] - LN.flat_fields(loop, x0 - hh * (A @ U0[y]), beta_of)[y]) / (2 * hh)
                  for y in range(N)], axis=0)
    return dict(gamma=gamma, w=w, rows=rows, c1=full(Jw @ P @ (t2 - 0.5 * t1)), c1_scale=full(Jw @ P @ t2))


def ito_zero(gamma, start):
    """The zero of part A's Itô drift on the simplex at γ (from `start`), and the drift's size there per unit γ."""
    from scipy.optimize import fsolve
    F, chart, coords, S0 = X.instance(gamma)
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)

    def ito(w):
        U0 = LN.flat_fields(loop, chart(w), beta_of)
        return GDR.limit_drift(F, chart, coords, w, U0.T @ U0 / U0.shape[0], h=X.H_ID, h2=X.H2_ID)[0] / gamma
    w = fsolve(ito, np.asarray(start, float)[:2], xtol=1e-12)
    return np.r_[w, 1 - w.sum()], float(np.abs(ito(w)).max())


def a_coords(gamma, w):
    """Katzenberger's drift (W_1, W_2, W_3 per unit v) of the Itô SDE in the coordinates with the amplitudes a = a_max tanh(α/a_max) for α — an
    Euler step in a —, and the closed form with (1 − q)(1 + 3q)."""
    F, chart, coords, S0 = X.instance(gamma)
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)
    am, ia = net.a_max, np.arange(0, 9, 3)

    def to_a(x):
        y = np.array(x, float); y[ia] = am * np.tanh(y[ia] / am); return y

    def from_a(y):
        x = np.array(y, float); x[ia] = am * np.arctanh(x[ia] / am); return x

    def Fa(y):
        f = np.array(F(from_a(y)), float); f[ia] *= 1 - (y[ia] / am) ** 2; return f
    w = np.asarray(w, float)
    x0 = chart(w)
    Ua = LN.flat_fields(loop, x0, beta_of); Ua[:, ia] *= 1 - (to_a(x0)[ia] / am) ** 2
    d = GDR.limit_drift(Fa, lambda w_: to_a(chart(w_)), lambda y: coords(from_a(y)), w, Ua.T @ Ua / Ua.shape[0], h=X.H_ID, h2=X.H2_ID)[0]
    W = np.r_[w, 1 - w.sum()]
    u = LN.selection(loop, E.exact_rep(net, W), 3 * W, coord='a')
    return np.r_[d, -d.sum()], gamma * W * (u - W @ u)


def g_drift(gamma, w, corr, h=None, parts=False):
    """Katzenberger's drift (W_1, W_2) of the filtered covariance at γ, the chart point w, for the draws `corr` at the step h (DT); with
    `parts` also A, P, Σ and the instance."""
    F, chart, coords, S0 = X.instance(gamma)
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)
    w = np.asarray(w, float); x0 = chart(w)
    U0 = LN.flat_fields(loop, x0, beta_of); Sig = U0.T @ U0 / U0.shape[0]
    A = GDR._jac(F, x0, X.H_ID); P = GDR.spectral_projection(A, GDR._jac(chart, w, X.H_ID))[0]
    d = GDR.limit_drift(F, chart, coords, w, LN.filtered_covariance(A, P, Sig, DT if h is None else h, corr), h=X.H_ID, h2=X.H2_ID)[0]
    return (d, A, P, Sig, (F, chart, coords, S0)) if parts else d


TAUC_G = (0.0,) + tuple(np.geomspace(5.0, 2e4, 22))


@cached
def g1_job(gamma):
    """Part G1: the zero of the filtered selection continued in τ_c (rows τ_c, W_1..3, its two stability exponents per unit drift scale);
    stops where no zero is left with every W_k above 1.0005 W_min,k."""
    from scipy.optimize import fsolve
    Wm = np.abs(E.A_BUMP) / A_MAX

    def dr(w, corr):
        W = np.r_[w, 1 - w.sum()]
        return np.full(2, np.nan) if np.any(W <= 1.0005 * Wm) else g_drift(gamma, w, corr)
    net, base = E.setup()
    P = E.G[0] * np.mean(base.features(E.exact_rep(net, np.ones(3) / 3))[4] ** 2, axis=0)
    start = ito_zero(gamma, LN.w_label(P, Wm))[0][:2]
    rows, last = [], None
    for tc in TAUC_G:
        corr = ('fresh',) if tc == 0 else ('ar1', float(np.exp(-DT / tc)))
        sc = np.abs(dr(start, corr)).max()
        w = fsolve(lambda w_: np.nan_to_num(dr(w_, corr) / sc, nan=1e3), start, xtol=1e-12)
        d = dr(w, corr)
        if not np.all(np.isfinite(d)) or np.abs(d).max() > 1e-6 * sc:
            last = dict(tc=tc, drift=np.r_[dr(start, corr), -dr(start, corr).sum()] / sc)
            break
        lam = np.linalg.eigvals(np.array([(dr(w + 1e-6 * e, corr) - dr(w - 1e-6 * e, corr)) / 2e-6 for e in np.eye(2)]).T / sc).real
        rows.append(np.r_[tc, w, 1 - w.sum(), np.sort(lam)])
        start = w
    return dict(gamma=gamma, rows=np.array(rows), stop=last, W_min=Wm)


@cached
def g2_job(args):
    """Part G2 at one point: the slow limit τ_cΣ_f against −A_DΣQᵀ (= m·QMQᵀ when Σ = −mAM) at τ_c = 10⁶, the latter's asymmetry, and the
    drifts (W_1..3) of the slow noise, the dither (Σ = 2τG), the dither's normal part and the fresh noise."""
    gamma, w = args
    tc = 1e6
    d_slow, A, P, Sig, (F, chart, coords, S0) = g_drift(gamma, w, ('ar1', float(np.exp(-DT / tc))), parts=True)
    n = A.shape[0]; Q = np.eye(n) - P
    AD = np.linalg.solve(A - np.abs(np.linalg.eigvals(A)).max() * P, Q)
    Mn = -AD @ Sig @ Q.T
    Sf = LN.filtered_covariance(A, P, Sig, DT, ('ar1', float(np.exp(-DT / tc))))
    lim = lambda S_: GDR.limit_drift(F, chart, coords, np.asarray(w, float), S_, h=X.H_ID, h2=X.H2_ID)[0]
    full = lambda d: np.r_[d, -d.sum()]
    return dict(gamma=gamma, w=np.asarray(w, float), lim_res=float(np.abs(tc * Sf - Mn).max() / np.abs(Mn).max()),
                asym=float(np.abs(Mn - Mn.T).max() / np.abs(Mn).max()), slow=full(d_slow), dither=full(lim(S0)), dither_n=full(lim(Q @ S0 @ Q.T)),
                fresh=full(lim(Sig)))


def part_B():
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_C)
    Th = E.exact_rep(net, np.ones(3) / 3)
    P = E.G[0] * np.mean(base.features(Th)[4] ** 2, axis=0)
    Wm = np.abs(E.A_BUMP) / A_MAX
    Wl = LN.w_label(P, Wm)

    def flow(w):                                                               # dW/ds per γv in the chart (W_1, W_2)
        W = np.r_[w, 1 - w.sum()]; u = LN.selection(loop, E.exact_rep(net, W), 3 * W)
        return (W * (u - W @ u))[:2]
    h = 1e-6
    Jac = np.array([(flow(Wl[:2] + h * e) - flow(Wl[:2] - h * e)) / (2 * h) for e in np.eye(2)]).T
    u_l = LN.selection(loop, E.exact_rep(net, Wl), 3 * Wl)
    z2, z2r = ito_zero(GAMMA_C, Wl)
    ac = [(w, *a_coords(1e-5, w)) for w in ((0.40, 0.25), (1 / 3, 1 / 3))]
    return dict(P=P, W_min=Wm, W_l=Wl, amp=np.abs(E.A_BUMP) / Wl / A_MAX, u_l=u_l, eig=np.linalg.eigvals(Jac),
                psi_l=LN.psi(Wl, P, Wm), psi_u=LN.psi(np.ones(3) / 3, P, Wm), psi_0=LN.psi(np.asarray(W0_C), P, Wm), W_sharp=Q.w_sharp(),
                zero2=z2, zero2_res=z2r, W_la=LN.w_label(P, Wm, coord='a'), acoords=ac)


@cached
def c_series(args):
    """The full loop with measurement noise only: (scheme, ρ, step, seed) → records (s, α_1..3, W_1..3, the fit error) and the gauge
    coordinate y = Φ(z) every REC_PHI slow units (rows s, y_1..3)."""
    scheme, rho, dt, seed, S, v, gamma, W0 = args
    t0 = time.time()
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)
    W0 = np.asarray(W0, float)
    rec, full = [], []
    every = int(round(REC_PHI / REC_S))

    def cb(s_, T, bb):
        if len(rec) % every == 0:
            full.append((s_, T.copy(), bb.copy()))
        rec.append(np.r_[s_, T[:, 0], bb / 3, loop.fit_error(T, bb)])
    loop.run_dithered(E.exact_rep(net, W0), 3 * W0, S, dt, 0.0, np.random.default_rng(seed), rec_every=int(round(REC_S / dt)), callback=cb,
                      label_v=v, label_scheme=scheme, label_rho=rho)
    F, chart, coords, S0 = X.instance(gamma)
    y, res = [], 0.0
    for s_, T, bb in full:
        xe, rs = GDR.limit_map(F, np.r_[T.ravel(), np.log(bb[:2] / bb[2])], L_MAP)
        we = coords(xe); y.append(np.r_[s_, we, 1 - we.sum()]); res = max(res, rs)
    print(f"[{time.time() - t0:.0f}s] C {scheme} rho {rho:+.3f} dt {dt:g} seed {seed}", flush=True)
    return dict(scheme=scheme, rho=rho, dt=dt, seed=seed, rec=np.array(rec), y=np.array(y), res=res)


def d_job(args):
    """Part D: the scalar factorisation with label noise, one scheme: the paths' mean of Q = u² − v² every 10 time units, its std, and the
    largest |Q − Q(0)| of a path (for Heun)."""
    scheme, rho = args
    rng = np.random.default_rng(SEED_D)
    n = int(round(T_D / DT_D)); h = DT_D
    u = np.full(N_PATHS_D, UV0_D[0]); v = np.full(N_PATHS_D, UV0_D[1])
    rho_ = rho if scheme in ('ar1', 'heun') else 0.0
    kb = int(rho) if scheme == 'block' else 1
    sx = np.sqrt(SIG2_D / h * (1 - rho_) / (1 + rho_) / kb)               # the same intensity σ² for every scheme
    xi = sx * rng.standard_normal(N_PATHS_D)
    Q0 = UV0_D[0] ** 2 - UV0_D[1] ** 2
    every = int(round(10.0 / h)); rec = []
    for k in range(n + 1):
        if k % every == 0:
            Q = u * u - v * v
            rec.append((k * h, Q.mean(), Q.std(), np.abs(Q - Q0).max(), np.abs(u * v - 1).mean()))
        if k == n:
            break
        if k > 0 and k % kb == 0:
            xi = rho_ * xi + np.sqrt(1 - rho_ ** 2) * sx * rng.standard_normal(N_PATHS_D)
        c = xi - (u * v - 1.0)
        if scheme == 'heun':
            up, vp = u + h * c * v, v + h * c * u
            cp = xi - (up * vp - 1.0)
            u, v = u + 0.5 * h * (c * v + cp * vp), v + 0.5 * h * (c * u + cp * up)
        else:
            u, v = u + h * c * v, v + h * c * u
    rec = np.array(rec)
    return dict(scheme=scheme, rho=rho, rec=rec, rate=-np.log(rec[-1, 1] / Q0) / (SIG2_D * T_D), lr=d_prediction(scheme, rho))


def d_prediction(scheme, rho):
    """Part D's linear-response rate: dQ/dt = −σ²(the per-step factor)φ(hλ(Q))Q, λ = √(Q² + 4) on uv = 1, integrated over [0, T_D] (zero for
    Heun); −log(Q(T)/Q(0))/(σ²T)."""
    if scheme == 'heun':
        return 0.0
    from scipy.integrate import solve_ivp
    f, c = factor_C(scheme, rho), corr_of(scheme, rho)
    Q0 = UV0_D[0] ** 2 - UV0_D[1] ** 2
    sol = solve_ivp(lambda t, Q: -SIG2_D * f * LN.lr_factor(np.sqrt(Q ** 2 + 4), DT_D, c) * Q, (0.0, T_D), [Q0], rtol=1e-10, atol=1e-12)
    return float(-np.log(sol.y[0, -1] / Q0) / (SIG2_D * T_D))


def diag_problem():
    rng = np.random.default_rng(SEED_E)
    X = rng.standard_normal((N_E, D_E))
    ws = np.zeros(D_E); ws[:3] = (1.0, -0.5, 0.8)
    return X, X @ ws, ws


def e_job(args):
    """Part E, one scheme: the gradient flow to the interpolant, then label noise; the products P_i = u_iv_i's log-decay rates over [0, T_E]
    (the paths' mean of P_i), against 4V_Ec_i/n, and ‖w‖₁, ‖w − w*‖₁, the residual at the start and the end."""
    scheme, rho = args
    X, y, ws = diag_problem()
    n, h = N_E, DT_E
    u = np.full((P_E, D_E), ALPHA_E); v = np.full((P_E, D_E), ALPHA_E)

    def grad(u, v, xi):
        g = ((u * u - v * v) @ X.T - y[None, :] - xi) @ X / n                   # ∂L/∂w, L = (1/2n)Σ_j(⟨x_j, w⟩ − y_j)² with noisy labels
        return 2 * g * u, -2 * g * v

    def flow(u, v):                                                            # the noiseless gradient flow to an interpolant
        for k_ in range(400000):
            gu, gv = grad(u, v, 0.0); u, v = u - h * gu, v - h * gv
            if k_ % 1000 == 0 and np.abs((u * u - v * v) @ X.T - y).max() < 1e-9:
                break
        return u, v
    u, v = flow(u, v)
    rng = np.random.default_rng(SEED_E + 1)
    rho_ = rho if scheme in ('ar1', 'heun') else 0.0
    kb = int(rho) if scheme == 'block' else 1
    sx = np.sqrt(V_E * n / h * (1 - rho_) / (1 + rho_) / kb)
    xi = sx * rng.standard_normal((P_E, n))
    every = int(round(T_E / h / 25)); Dm = []                                  # the paths' mean of u² + v², for the prediction
    P0 = (u * v).mean(axis=0); w0 = (u * u - v * v)
    c = (X ** 2).sum(axis=0)
    trH = lambda u, v: float(((4 / n) * (c[None, :] * (u * u + v * v)).sum(axis=1)).mean())     # tr ∇²L at an interpolant (Gauss–Newton)
    st = dict(l1=float(np.abs(w0).sum(axis=1).mean()), err=float(np.abs(w0 - ws).sum(axis=1).mean()), res=float(np.abs(w0 @ X.T - y).max()),
              tr=trH(u, v))
    for k in range(int(round(T_E / h)) + 1):
        if k % every == 0:
            Dm.append((u * u + v * v).mean(axis=0))
        if k == int(round(T_E / h)):
            break
        if k > 0 and k % kb == 0:
            xi = rho_ * xi + np.sqrt(1 - rho_ ** 2) * sx * rng.standard_normal((P_E, n))
        gu, gv = grad(u, v, xi)
        if scheme == 'heun':
            up, vp = u - h * gu, v - h * gv
            gu2, gv2 = grad(up, vp, xi)
            u, v = u - 0.5 * h * (gu + gu2), v - 0.5 * h * (gv + gv2)
        else:
            u, v = u - h * gu, v - h * gv
    P1 = (u * v).mean(axis=0)
    rate = -np.log(P1 / P0) / T_E
    lr = np.zeros(D_E)
    if scheme != 'heun':                                                       # Σ_jφ(hκ_j)(X_iᵀq_j)²/c_i along the run, times the factor
        for i, D_ in enumerate(Dm):
            kap, qv = np.linalg.eigh((4 / n) * (X * D_[None, :]) @ X.T)
            s_ = (LN.lr_factor(kap, h, corr_of(scheme, rho))[:, None] * (qv.T @ X) ** 2).sum(axis=0) / c
            lr += s_ * (0.5 if i in (0, len(Dm) - 1) else 1.0) / (len(Dm) - 1)
        lr *= factor_C(scheme, rho)
    u, v = flow(u, v); w1 = u * u - v * v                                      # the interpolant reached from the end state
    return dict(scheme=scheme, rho=rho, ratio=rate / (4 * V_E * c / n), lr=lr, start=st,
                end=dict(l1=float(np.abs(w1).sum(axis=1).mean()), err=float(np.abs(w1 - ws).sum(axis=1).mean()), res=float(np.abs(w1 @ X.T - y).max()),
                         tr=trH(u, v)))


def mlp_parts():
    xs = np.linspace(-1, 1, N_F); ys = np.sin(np.pi * xs)
    rng = np.random.default_rng(SEED_F)
    t0 = np.r_[0.5 * rng.standard_normal(H_F), 4 * rng.standard_normal(H_F), 2 * rng.standard_normal(H_F)]

    def jac(t):                                                                # f (P, n) and ∇_θf (P, 3H, n)
        a, w, b = t[..., :H_F], t[..., H_F:2 * H_F], t[..., 2 * H_F:]
        th = np.tanh(w[..., :, None] * xs[None, None, :] + b[..., :, None]); s2 = 1 - th ** 2
        f = np.einsum('ph,phn->pn', a, th)
        return f, np.concatenate([th, a[..., :, None] * s2 * xs[None, None, :], a[..., :, None] * s2], axis=1)

    def project(t, iters=30):                                                  # Gauss–Newton (minimum-norm) steps onto Γ
        for _ in range(iters):
            f, J = jac(t)
            d = np.linalg.solve(np.einsum('pin,pim->pnm', J, J), (f - ys[None, :])[..., None])[..., 0]
            t = t - np.einsum('pin,pn->pi', J, d)
        return t
    sharp = lambda t: (jac(t)[1] ** 2).sum(axis=(1, 2)) / N_F
    return xs, ys, t0, jac, project, sharp


def f_prediction():
    """The leading-order Itô slope of T at the start, −(V_F/4)|P_T∇T|², and the normal relaxation rates there."""
    xs, ys, t0, jac, project, sharp = mlp_parts()
    t = project(t0[None, :])[0]
    hh = 1e-5
    gT = np.array([(sharp((t + hh * e)[None, :])[0] - sharp((t - hh * e)[None, :])[0]) / (2 * hh) for e in np.eye(3 * H_F)])
    J = jac(t[None, :])[1][0]                                                  # (3H, n)
    PT = np.eye(3 * H_F) - J @ np.linalg.solve(J.T @ J, J.T)
    lam = np.linalg.eigvalsh(J.T @ J / N_F)
    lr = []                                                                    # −(V/4)(the per-step factor)∇T·P_T∇trΨ(K), Ψ′ = φ
    for scheme, rho in SCHEMES_D:
        if scheme == 'heun':
            lr.append(0.0); continue
        trpsi = lambda t_: float(LN.lr_potential(np.linalg.eigvalsh((lambda J_: J_.T @ J_)(jac(t_[None, :])[1][0]) / N_F), DT_F,
                                                 corr_of(scheme, rho)).sum())
        gP = np.array([(trpsi(t + hh * e) - trpsi(t - hh * e)) / (2 * hh) for e in np.eye(3 * H_F)])
        lr.append(-(V_F / 4) * factor_C(scheme, rho) * float(gT @ PT @ gP))
    return dict(slope=-(V_F / 4) * float(gT @ PT @ gT), T0=float(sharp(t[None, :])[0]), lam=lam, lr=np.array(lr),
                resid=float(np.abs(jac(t[None, :])[0][0] - ys).max()))


def f_job(args):
    """Part F, one scheme: the sharpness of the projection onto Γ at 11 times, its initial slope (a quadratic fit)."""
    scheme, rho = args
    xs, ys, t0, jac, project, sharp = mlp_parts()
    t = project(np.repeat(t0[None, :], P_F, axis=0))
    rng = np.random.default_rng(SEED_NF)
    rho_ = rho if scheme in ('ar1', 'heun') else 0.0
    kb = int(rho) if scheme == 'block' else 1
    h = DT_F
    sx = np.sqrt(V_F * N_F / h * (1 - rho_) / (1 + rho_) / kb)
    xi = sx * rng.standard_normal((P_F, N_F))
    N = int(round(T_F / h)); every = N // 10; rec = []
    grad = lambda t, xi: (lambda f, J: np.einsum('pin,pn->pi', J, f - ys[None, :] - xi) / N_F)(*jac(t))
    for k in range(N + 1):
        if k % every == 0:
            rec.append((k * h, float(sharp(project(t.copy())).mean())))
        if k == N:
            break
        if k > 0 and k % kb == 0:
            xi = rho_ * xi + np.sqrt(1 - rho_ ** 2) * sx * rng.standard_normal((P_F, N_F))
        g = grad(t, xi)
        if scheme == 'heun':
            t = t - 0.5 * h * (g + grad(t - h * g, xi))
        else:
            t = t - h * g
    rec = np.array(rec)
    return dict(scheme=scheme, rho=rho, rec=rec, slope=float(np.polyfit(rec[:, 0], rec[:, 1], 2)[1]))


def name_C(scheme, rho, dt):
    out = {'ito': 'afresh each step', 'heun': 'Heun (Stratonovich)', 'block': f'redrawn every {int(rho)} steps'}.get(scheme, f'AR(1), ρ = {rho:+.3f}')
    if scheme == 'heun' and rho:
        out += f', AR(1) draws ρ = {rho:+.3f}'
    return out if dt == DT else f'{out}, step {dt:g}'


def factor_C(scheme, rho):
    """The per-step variance over the fresh noise's at the same intensity (the selection's ratio when λ_maxτ_c ≪ 1); zero for Heun."""
    return 0.0 if scheme == 'heun' else (1.0 / int(rho) if scheme == 'block' else LN.ar1_factor(rho))


def corr_of(scheme, rho):
    """The draws' step correlations for `label_noise.corr_sums`."""
    return {'ito': ('fresh',), 'ar1': ('ar1', rho), 'block': ('block', int(rho)), 'heun': ('ar1', rho)}[scheme]


def c_prediction():
    """Part C's linear-response predictions at W₀ per unit v, one per scheme of SCHEMES_C: the per-step factor times Katzenberger's drift of
    the filtered covariance Σ_f at the scheme's step (zero for Heun), (W_1, W_2, W_3)."""
    F, chart, coords, S0 = X.instance(GAMMA_C)
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_C)
    w = np.asarray(W0_C[:2], float); x0 = chart(w)
    U0 = LN.flat_fields(loop, x0, beta_of); Sig = U0.T @ U0 / U0.shape[0]
    A = GDR._jac(F, x0, X.H_ID); P = GDR.spectral_projection(A, GDR._jac(chart, w, X.H_ID))[0]
    out = []
    for scheme, rho, dt in SCHEMES_C:
        d = np.zeros(2) if scheme == 'heun' else factor_C(scheme, rho) * GDR.limit_drift(
            F, chart, coords, w, LN.filtered_covariance(A, P, Sig, dt, corr_of(scheme, rho)), h=X.H_ID, h2=X.H2_ID)[0]
        out.append(np.r_[d, -d.sum()])
    return np.array(out)


def summarise_C(runs, ito0, lr=None):
    """Per scheme: the gauge coordinate's least-squares slope (per unit slow time) and end displacement y(S_C) − W₀, mean ± s.e. over the
    seeds, against the prediction (1 − ρ)/(1 + ρ)·v·(Itô drift at W₀), and their ratios to v·(Itô drift at W₀) on its largest component;
    with `lr` (c_prediction) also against the linear-response prediction v·lr, the slope over it per component (± s.e.)."""
    out = []
    W0 = np.asarray(W0_C, float)
    se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M))
    k = int(np.argmax(np.abs(ito0)))
    for i_, (scheme, rho, dt) in enumerate(SCHEMES_C):
        rs = [r for r in runs if r['scheme'] == scheme and abs(r['rho'] - rho) < 1e-12 and r['dt'] == dt]
        R = np.array([r['rec'] for r in rs])
        slope = np.array([np.polyfit(r['y'][:, 0], r['y'][:, 1:], 1)[0] for r in rs])
        endd = np.array([r['y'][-1, 1:] for r in rs]) - W0[None, :]
        f = factor_C(scheme, rho)
        out.append(dict(scheme=scheme, rho=rho, dt=dt, factor=f, slope=slope.mean(axis=0), slope_se=se(slope), end=endd.mean(axis=0),
                        end_se=se(endd), pred=f * V_C * ito0, ratio=slope.mean(axis=0)[k] / (V_C * ito0[k]),
                        ratio_se=se(slope)[k] / abs(V_C * ito0[k]), k=k, almax=float(np.abs(R[:, :, 1:4]).max()), fit=float(R[:, :, 7].mean()),
                        res=max(r['res'] for r in rs), n=len(rs)))
        if lr is not None:
            out[-1].update(lr=V_C * lr[i_], over_lr=(slope.mean(axis=0) / (V_C * lr[i_]) if scheme != 'heun' else np.full(3, np.nan)),
                           over_lr_se=(se(slope) / np.abs(V_C * lr[i_]) if scheme != 'heun' else np.full(3, np.nan)))
    return out


def report(A=None, B=None, C=None, D=None, Ee=None, Ff=None, A2=None, G=None):
    fmt = E.fmt
    v3 = lambda v, p=4: ', '.join(f'{x:.{p}f}' for x in v)
    L_ = ["# Proposition D‴.41: measurement noise on the gauge (`exp_gauge_label.py`)", "",
          "Proposition D‴.36's instance (three bumps, one particle per class, m = 3). The measured error r + ξ drives both laws; ξ is "
          "independent over the quadrature nodes, of intensity v. Drifts are per unit slow time and per unit v.", ""]
    if A is not None:
        L_ += ["## A. The identity: Itô drift + Wong–Zakai term = the Stratonovich drift = 0", "",
               "| γ | W | the Itô drift ½∂²Φ[Σ] (W_1, W_2, W_3) | the Wong–Zakai term | their sum | its size over the Itô drift's | the diffusion along M (largest entry) | "
               "the closed form (slow-replicator limit) | its largest relative difference from the Itô drift |",
               "|---|---|---|---|---|---|---|---|---|"]
        for r in A:
            W = np.r_[r['w'], 1 - r['w'].sum()]; tot = r['ito'] + r['wz']
            rel = np.abs(r['closed'] - r['ito']).max() / np.abs(r['ito']).max()
            L_.append(f"| {r['gamma']:g} | {v3(W, 3)} | {', '.join(f'{v:+.4e}' for v in r['ito'])} | {', '.join(f'{v:+.4e}' for v in r['wz'])} | "
                      f"{', '.join(f'{v:+.1e}' for v in tot)} | {np.abs(tot).max() / np.abs(r['ito']).max():.1e} | {np.abs(r['D']).max():.1e} | "
                      f"{', '.join(f'{v:+.4e}' for v in r['closed'])} | {rel:.4f} |")
        L_ += ["", f"The exact fit's residual at these points: {max(r['fit'] for r in A):.1e}; the normal relaxation rates "
               f"{np.abs(A[0]['lam']).min():.3e} … {np.abs(A[0]['lam']).max():.3e}.", ""]
        if 'ito_h' in A[0]:
            L_ += ["Heteroscedastic noise (the variance along the orbit modulated by 1 + 0.8 sin t, mean 1): the Itô drift and the Stratonovich "
                   "drift (their sum with the Wong–Zakai term) —", "",
                   "| γ | W | the Itô drift | the Stratonovich drift | its size over the Itô drift's |", "|---|---|---|---|---|"]
            for r in A:
                W = np.r_[r['w'], 1 - r['w'].sum()]; tot = r['ito_h'] + r['wz_h']
                L_.append(f"| {r['gamma']:g} | {v3(W, 3)} | {', '.join(f'{x:+.4e}' for x in r['ito_h'])} | {', '.join(f'{x:+.3e}' for x in tot)} | "
                          f"{np.abs(tot).max() / np.abs(r['ito_h']).max():.2e} |")
            L_.append("")
    if A2 is not None:
        gs = sorted({r['gamma'] for r in A2}, reverse=True)
        L_ += ["## A2. The noise's structure over the nodes: the Stratonovich drift over its own Itô drift (largest component)", "",
               "The drift of the Stratonovich (smooth-noise) limit, ½∂²Φ[Σ_C] + the Wong–Zakai term of C + the one-sided kernel's Lévy-area "
               "part, over the Itô drift ½∂²Φ[Σ_C] of the same covariance C over the nodes. In brackets: the Lévy-area part alone over the Itô drift.", ""]
        for w in W_A2:
            rs = {r['gamma']: r for r in A2 if np.allclose(r['w'], w)}
            L_ += [f"W = {v3(np.r_[w, 1 - sum(w)], 3)}:", "", "| the noise over the nodes | " + ' | '.join(f'γ = {g:g}' for g in gs) + " |",
                   "|---|" + "---|" * len(gs)]
            for name in rs[gs[0]]['rows']:
                cells = []
                for g in gs:
                    q = rs[g]['rows'][name]; it = np.abs(q['ito']).max()
                    c = f"{np.abs(q['ito'] + q['sym'] + q['lev']).max() / it:.2e}"
                    if np.abs(q['lev']).max() > 0:
                        c += f" ({np.abs(q['lev']).max() / it:.1e})"
                    cells.append(c)
                L_.append(f"| {name} | " + ' | '.join(cells) + " |")
            L_.append("| the coloured-noise identity ⟨∇ψ·AU⟩ = −½P⟨∂²F[U, U]⟩ + P⟨(AU·∇)U⟩, over its second term | " +
                      ' | '.join(f"{np.abs(rs[g]['c1']).max() / np.abs(rs[g]['c1_scale']).max():.1e}" for g in gs) + " |")
            L_.append("")
    if B is not None:
        L_ += ["## B. The selection's target", "",
               f"P_k = G_α⟨Φ_k²⟩_ρ = {v3(B['P'])}; W_min = {v3(B['W_min'])}.", "",
               f"W_ℓ = argmin Ψ_ℓ = {v3(B['W_l'])} (|a_k|/a_max there {v3(B['amp'], 3)}; ū there {v3(B['u_l'], 5)}, equal); "
               f"Proposition D‴.38's W♯ = {v3(B['W_sharp'])}.", "",
               f"Ψ_ℓ/v: {B['psi_l']:.5f} at W_ℓ, {B['psi_u']:.5f} at 1/3, {B['psi_0']:.5f} at W₀ = {v3(W0_C, 2)}; the closed form's replicator "
               f"flow linearised at W_ℓ has the eigenvalues {', '.join(f'{v.real:+.5f}' for v in B['eig'])} (per γv, per unit slow time).", ""]
        if 'zero2' in B:
            L_ += [f"W_ℓ is the slow-replicator limit's target: at γ = {GAMMA_C:g} the Itô drift vanishes at {v3(B['zero2'])} (the drift there "
                   f"{B['zero2_res']:.1e} per unit γ).", "",
                   "The update's coordinates: an Euler step in the amplitude a = a_max tanh(α/a_max) instead of α (the Itô SDE in those coordinates; "
                   f"γ = 10⁻⁵) — the closed form's factor 1 − q² becomes (1 − q)(1 + 3q), its target W_ℓᵃ = {v3(B['W_la'])}:", "",
                   "| W | Katzenberger's drift, Euler in a | the closed form with (1 − q)(1 + 3q) | its largest relative difference |", "|---|---|---|---|"]
            for w, d, cl in B['acoords']:
                L_.append(f"| {v3(np.r_[w, 1 - sum(w)], 3)} | {', '.join(f'{x:+.4e}' for x in d)} | {', '.join(f'{x:+.4e}' for x in cl)} | "
                          f"{np.abs(d - cl).max() / np.abs(d).max():.1e} |")
            L_.append("")
    if C is not None:
        Ci, ito0 = C['sum'], C['ito0']
        L_ += [f"## C. The full loop with measurement noise only (τ = 0, γ = {GAMMA_C:g}, v = {V_C:g}, from W₀ = {v3(W0_C, 2)}, "
               f"{Ci[0]['n']} seeds, {S_C:g} slow units, steps of {DT:g} unless stated)", "",
               f"The gauge coordinate y(s) = Φ(z(s)): the noiseless flow ({L_MAP:g} slow units) from the state every {REC_PHI:g} slow units. Part A's "
               f"Itô drift at W₀ (γ = {GAMMA_C:g}, per unit v): {', '.join(f'{v:+.4e}' for v in ito0)}; the ratio is taken on W_{Ci[0]['k'] + 1}. "
               "The linear-response prediction: the per-step factor times Katzenberger's drift of the filtered covariance Σ_f at the scheme's step "
               "(the increment's second moment, `label_noise.filtered_covariance`), times v; zero for Heun's scheme.", ""]
        if 'lr' in Ci[0]:
            L_ += ["| noise | (1 − ρ)/(1 + ρ) | dy/ds, the least-squares slope (± s.e.) | predicted (1 − ρ)/(1 + ρ)·v·(Itô drift) | "
                   "the linear-response prediction | the ratio to v·(Itô drift) (± s.e.) | the slope over the linear-response prediction (± s.e.) | "
                   "y(S) − W₀ (± s.e.) | the largest |α| | the mean fit error |",
                   "|---|---|---|---|---|---|---|---|---|---|"]
            for r in Ci:
                ol = ('—' if r['scheme'] == 'heun' else ', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(r['over_lr'], r['over_lr_se'])))
                L_.append(f"| {name_C(r['scheme'], r['rho'], r['dt'])} | {r['factor']:.4f} | "
                          f"{', '.join(f'{a:+.3e} ± {b:.1e}' for a, b in zip(r['slope'], r['slope_se']))} | "
                          f"{', '.join(f'{v:+.3e}' for v in r['pred']) if r['scheme'] != 'heun' else '0'} | "
                          f"{', '.join(f'{v:+.3e}' for v in r['lr']) if r['scheme'] != 'heun' else '0'} | {r['ratio']:+.4f} ± {r['ratio_se']:.4f} | "
                          f"{ol} | {', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(r['end'], r['end_se']))} | {r['almax']:.2f} | {fmt(r['fit'])} |")
        else:
            L_ += ["| noise | (1 − ρ)/(1 + ρ) | dy/ds, the least-squares slope (± s.e.) | predicted (1 − ρ)/(1 + ρ)·v·(Itô drift) | the ratio to v·(Itô drift) (± s.e.) | "
                   "y(S) − W₀ (± s.e.) | the largest |α| | the mean fit error |",
                   "|---|---|---|---|---|---|---|---|"]
            for r in Ci:
                L_.append(f"| {name_C(r['scheme'], r['rho'], r['dt'])} | {r['factor']:.3f} | {', '.join(f'{a:+.3e} ± {b:.1e}' for a, b in zip(r['slope'], r['slope_se']))} | "
                          f"{', '.join(f'{v:+.3e}' for v in r['pred'])} | {r['ratio']:+.4f} ± {r['ratio_se']:.4f} | "
                          f"{', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(r['end'], r['end_se']))} | {r['almax']:.2f} | {fmt(r['fit'])} |")
        L_ += ["", f"The limit map's residual drift at its end: at most {max(r['res'] for r in Ci):.1e}.", ""]
    if D is not None:
        L_ += [f"## D. The example: the scalar factorisation f = uv with label noise (σ² = {SIG2_D:g}, step {DT_D:g}, T = {T_D:g}, {N_PATHS_D} paths, "
               f"(u, v)(0) = {UV0_D})", "",
               "Q = u² − v², Q(0) = " + f"{UV0_D[0] ** 2 - UV0_D[1] ** 2:g}" + ". The Euler step gives Q_{n+1} = Q_n(1 − (ξ_n − r_n)²h²) exactly: the decay is set by "
               "the per-step second moment of the error ξ − r, the per-step variance times the linear-response factor φ(hλ), λ = u² + v².", ""]
        has_lr = 'lr' in D[0]
        L_ += ["| noise | the per-step factor | E[Q] at t = 0, 20, 40, …, 100 | the decay rate −log(E[Q(T)]/Q(0))/(σ²T) | " +
               ("the linear-response rate | " if has_lr else "") + "the largest |Q − Q(0)| of a path | the mean |uv − 1| at T |",
               "|---|---|---|---|" + ("---|" if has_lr else "") + "---|---|"]
        for r in D:
            rec = r['rec']
            L_.append(f"| {name_C(r['scheme'], r['rho'], DT)} | {factor_C(r['scheme'], r['rho']):.4f} | "
                      f"{', '.join(f'{x:.4f}' for x in rec[::2, 1])} | {r['rate']:+.4f} | " + (f"{r['lr']:+.4f} | " if has_lr else "") +
                      f"{rec[-1, 3]:.2e} | {rec[-1, 4]:.4f} |")
        L_ += ["", "For noise drawn afresh at each step the Itô SDE's E[Q] decays exactly as e^{−σ²t} (dQ = −σ²Q dt + 0·dW); Heun's scheme keeps "
               "every path's Q to its discretisation error, with fresh or with correlated draws.", ""]
    if Ee is not None:
        has_lr = 'lr' in Ee[0]
        L_ += [f"## E. A diagonal linear network with label noise (d = {D_E}, n = {N_E}, α = {ALPHA_E:g}, V = {V_E:g}, step {DT_E:g}, T = {T_E:g}, "
               f"{P_E} paths)", "",
               "| noise | the per-step factor | the decay rates of P_i = u_iv_i over 4Vc_i/n: mean (range over i) | " +
               ("the linear-response prediction: mean (range) | measured over predicted: mean (range) | " if has_lr else "") +
               "the interpolant's tr ∇²L: start → end | its ‖w − w*‖₁: start → end | its ‖w‖₁: start → end | its residual |",
               "|---|---|---|" + ("---|---|" if has_lr else "") + "---|---|---|---|"]
        for r in Ee:
            lrc = ""
            if has_lr:
                lrc = (f"{r['lr'].mean():+.4f} ({r['lr'].min():+.4f} … {r['lr'].max():+.4f}) | " +
                       ("— | " if r['scheme'] == 'heun' else f"{(r['ratio'] / r['lr']).mean():.4f} ({(r['ratio'] / r['lr']).min():.4f} … {(r['ratio'] / r['lr']).max():.4f}) | "))
            L_.append(f"| {name_C(r['scheme'], r['rho'], DT)} | {factor_C(r['scheme'], r['rho']):.4f} | {r['ratio'].mean():+.4f} "
                      f"({r['ratio'].min():+.4f} … {r['ratio'].max():+.4f}) | " + lrc + f"{r['start']['tr']:.4f} → {r['end']['tr']:.4f} | "
                      f"{r['start']['err']:.4f} → {r['end']['err']:.4f} | {r['start']['l1']:.4f} → {r['end']['l1']:.4f} | {r['end']['res']:.1e} |")
        L_ += ["", f"‖w*‖₁ = {np.abs(diag_problem()[2]).sum():g}. The interpolant: the noiseless gradient flow from the state (at the start, after the "
               "gradient flow from u = v = α; at the end, after the noisy phase). The Euler step multiplies each P_i by exactly "
               "1 − 4h²(g_i − ξ̃_i)²: the label noise's selection — the effective initialisation scale shrinking, the interpolant descending "
               "tr ∇²L = (4/n)Σ_ic_i√(w_i² + 4P_i²), a weighted ℓ¹ norm Σ_ic_i|w_i| as P → 0 — runs at 4h²E[(X_iᵀe)²]/n² per step, e = r − ξ: the per-step "
               "variance times Σ_jφ(hκ_j)(X_iᵀq_j)²/c_i over the residual's modes (κ_j, q_j) of K = (4/n)X diag(u² + v²)Xᵀ, a factor that "
               "differs from coordinate to coordinate when the noise is slow; Heun's scheme keeps the P_i, and the interpolant does not move.", ""]
    if Ff is not None:
        pr = Ff['pred']
        has_lr = 'lr' in pr
        L_ += [f"## F. A generic network: {H_F} tanh units, {N_F} points of sin(πx), V = {V_F:g}, step {DT_F:g}, T = {T_F:g}, {P_F} paths", "",
               f"At the start (on Γ, residual {pr['resid']:.1e}): the sharpness T = {pr['T0']:.5f}; the normal relaxation rates (eigenvalues of "
               f"K = JᵀJ/n) {pr['lam'].min():.3e} … {pr['lam'].max():.3e}; the leading-order Itô slope −(V/4)|P_T∇T|² = {pr['slope']:.4e}." +
               (" The linear-response slope: −(V/4)(the per-step factor)∇T·P_T∇trΨ(K), Ψ′ = φ (`label_noise.lr_potential`)." if has_lr else ""), "",
               "| noise | the per-step factor | T of the projection onto Γ at t = 0, 40, …, 200 | dT/dt at 0 (a quadratic fit) | "
               "over the Itô prediction | over the measured Itô slope | " + ("the linear-response slope | measured over it | " if has_lr else ""),
               "|---|---|---|---|---|---|" + ("---|---|" if has_lr else "")]
        s0 = next(r['slope'] for r in Ff['runs'] if r['scheme'] == 'ito')
        for i_, r in enumerate(Ff['runs']):
            lrc = ""
            if has_lr:
                lr_ = pr['lr'][[k for k, sc in enumerate(SCHEMES_D) if sc == (r['scheme'], r['rho'])][0]]
                lrc = f"{lr_:+.3e} | " + ("— | " if r['scheme'] == 'heun' else f"{r['slope'] / lr_:+.3f} | ")
            L_.append(f"| {name_C(r['scheme'], r['rho'], DT)} | {factor_C(r['scheme'], r['rho']):.4f} | {', '.join(f'{x:.5f}' for x in r['rec'][::2, 1])} | "
                      f"{r['slope']:+.3e} | {r['slope'] / pr['slope']:+.3f} | {r['slope'] / s0:+.3f} | " + lrc)
        L_ += ["", "No exact identity holds here: the leading-order law (the per-step factor times the linear-response weights, Heun's zero) on a "
               "generic network.", ""]
    if G is not None:
        L_ += [f"## G. The noise's colour moves the target (AR(1) draws of correlation time τ_c at the step {DT:g})", "",
               "G1. The zero of the per-update selection on the simplex, continued in τ_c (τ_c = 0: fresh draws), with its two stability exponents "
               "(the eigenvalues of the drift's Jacobian over the drift's scale at the previous zero):", "",
               "| γ | τ_c | the zero W | W − W_min | the stability exponents |", "|---|---|---|---|---|"]
        for g in G['G1']:
            for r in g['rows']:
                L_.append(f"| {g['gamma']:g} | {r[0]:.4g} | {v3(r[1:4])} | {', '.join(f'{x:.2e}' for x in r[1:4] - g['W_min'])} | "
                          f"{', '.join(f'{x:+.3g}' for x in r[4:6])} |")
            if g['stop'] is not None:
                L_.append(f"| {g['gamma']:g} | {g['stop']['tc']:.4g} | no zero with W > 1.0005 W_min: the drift at the last zero "
                          f"{', '.join(f'{x:+.3f}' for x in g['stop']['drift'])} (over its scale) | | |")
        L_ += ["", "G2. The slow limit: τ_cΣ_f → −A_DΣQᵀ = m·QMQᵀ (A_D the group inverse of A, Q = I − P), the flow's own mobility on the normal "
               "space; the selection's direction against the dither's (Σ = 2τG) and the fresh noise's —", "",
               "| γ | W | τ_cΣ_f against −A_DΣQᵀ (τ_c = 10⁶) | −A_DΣQᵀ's asymmetry | the slow noise's drift (normalised) | the dither's | "
               "cos(slow, dither) | cos(slow, the dither's normal part) | cos(slow, fresh) | the dither's drift over its normal part's | "
               "the slow noise's drift over the dither's normal part's, × 2ττ_c/m (the effective temperature m·v_h/(2τ_c) over τ) |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        nrm = lambda d: d / np.abs(d).max()
        cos = lambda a, b: float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
        for r in G['G2']:
            L_.append(f"| {r['gamma']:g} | {v3(np.r_[r['w'], 1 - r['w'].sum()], 3)} | {r['lim_res']:.1e} | {r['asym']:.1e} | {v3(nrm(r['slow']))} | "
                      f"{v3(nrm(r['dither']))} | {cos(r['slow'], r['dither']):+.5f} | {cos(r['slow'], r['dither_n']):+.5f} | {cos(r['slow'], r['fresh']):+.4f} | "
                      f"{', '.join(f'{a / b:.4f}' for a, b in zip(r['dither'], r['dither_n']))} | "
                      f"{', '.join(f'{a / b * 2 * GD.TAU * 1e6 / 3:.4f}' for a, b in zip(r['slow'], r['dither_n']))} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_gauge_label.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_label_raw.npz')
    md = os.path.join(OUT, 'gauge_label.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEF'], A2=get('A2'), G=get('G')))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFG')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'C' in parts:
            h['C'] = pool.map_async(c_series, [(sc, rho, dt, sd, S_C, V_C, GAMMA_C, W0_C) for sc, rho, dt in SCHEMES_C for sd in SEEDS_C],
                                    chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, list(SCHEMES_D), chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(e_job, list(SCHEMES_D), chunksize=1)
        if 'F' in parts:
            h['F'] = pool.map_async(f_job, list(SCHEMES_D), chunksize=1)
        if 'A' in parts or 'C' in parts:
            h['A'] = pool.map_async(a_job, [(g, w) for g in GAMMAS_A for w in W_A], chunksize=1)
        if 'A' in parts:
            h['A2'] = pool.map_async(a2_job, [(g, w) for g in GAMMAS_A2 for w in W_A2], chunksize=1)
        if 'G' in parts:
            h['G1'] = pool.map_async(g1_job, [1e-2, 1e-5], chunksize=1)
            h['G2'] = pool.map_async(g2_job, [(g, w) for g in (1e-2, 1e-5) for w in ((0.40, 0.25), (1 / 3, 1 / 3), (0.30, 0.45))], chunksize=1)
        if 'A' in h:
            res['A'] = h['A'].get(); E.save(raw, **res); log('A')
        if 'A2' in h:
            res['A2'] = h['A2'].get(); E.save(raw, **res); log('A2')
        if 'B' in parts:
            res['B'] = part_B(); E.save(raw, **res); log('B')
        for k in 'DE':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
        if 'F' in h:
            res['F'] = dict(pred=f_prediction(), runs=h['F'].get()); E.save(raw, **res); log('F')
        if 'G1' in h:
            res['G'] = dict(G1=h['G1'].get(), G2=h['G2'].get()); E.save(raw, **res); log('G')
        if 'C' in h:
            runs = h['C'].get()
            ito0 = next(r['ito'] for r in res['A'] if r['gamma'] == GAMMA_C and np.allclose(np.r_[r['w'], 1 - r['w'].sum()], W0_C))
            res['C_runs'] = runs
            lr = c_prediction()
            res['C'] = dict(ito0=ito0, lr=lr, sum=summarise_C(runs, ito0, lr))
            E.save(raw, **res); log('C')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEF'], A2=res.get('A2'), G=res.get('G')))


if __name__ == '__main__':
    main()
