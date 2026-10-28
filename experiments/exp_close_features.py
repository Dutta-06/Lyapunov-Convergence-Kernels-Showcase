"""
Proposition D‴.52: close features — the transit counts them at the number-detection limit and its block places them at the support limit
====================================================================================================================================

The sliding inverse scale space (Proposition D‴.48(iii); `rich_transit.sliding_iss` with the block's least-squares polish) on two Gaussian
bumps of the kernels' width σ = 0.12 a distance Δ apart (61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, equal masses).
With A = A₁ + A₂ ≠ 0 the pair is a monopole at the centroid c̄ = (A₁c₁ + A₂c₂)/A (kept at 0.032) with the moments V = A₁A₂Δ²/A,
T₃ = A₁A₂(A₁ − A₂)Δ³/A²: f* = Aφ(c̄) + ½V∂²_cφ(c̄) + ⅙T₃∂³_cφ(c̄) + O(Δ⁴); the first entrant fits Aφ(c̄) and leaves ½Vψ,
ψ = ∂²_cφ − (⟨∂²_cφ φ⟩/⟨φφ⟩)φ, whose correlations ½V⟨ψφ_j⟩ drive the second entry: Δτ₂ = min_j (s − u_j)/v_j over the walls s = ±1
the duals u_j(τ₁) = ⟨f*φ_j⟩/|⟨f*φ_{j₁}⟩| reach at the rates v_j = ½V⟨ψφ_j⟩ — ∝ Δ⁻²; a dipole (A = 0) has its first entry and the second
∝ Δ⁻¹.

A. The second entry against Δ (no noise): A = (0.8, 0.6), (0.8, 0.8), (0.8, −0.4), (0.8, −0.8): τ₁, the entrants, Δτ₂ against the
   leading-order prediction (for the dipole the transit of its derivative shape −A₁∂_cφ, slowed by 1/Δ), the scaled delays Δτ₂Δ² (Δτ₂Δ
   for the dipole), and the path's end (two kernels at the features?). The expansion parameter is |A₁ − A₂|Δ/(|A|σ).
B. The two-kernel block at its exact fit: the Gauss–Newton Hessian's eigenvalues and the fast flow's relaxation rates (the critical
   curve's flow linearised at the fit) against Δ — ∝ 1, 1, Δ², Δ⁶.
C. Number detection under label noise: y = f* + ξ, ξ i.i.d. N(0, s²) on the samples, A = (0.8, 0.6); under H1 the pair, under H0 the
   one-bump fit closest to it (the single-kernel least-squares fit of f*); R_C draws each; the second entry's delay Δτ₂ as the statistic:
   the power at the false-alarm rate 5% — the threshold the H0 draws' quantile at the true offset: an oracle-calibrated simple null
   (referee round 55) — against Δ for each s, beside the Neyman–Pearson bound Φ(‖f₁ − f₀‖/s − z₀.₉₅) of the same simple hypotheses; the
   separation Δ₅₀ where the power crosses ½ (with bootstrap intervals), its slope in log s (½ predicted: Δ² ∝ s).
D. Placement under label noise: the two-kernel block after the second entry on the H1 draws — the rms errors of the positions and of the
   amplitude ratio against s/Δ² and s/Δ³ (the moments' singular values Δ, Δ³) and against their Cramér–Rao bounds at the truth; how many
   blocks hold both kernels within σ of the centroid, and how many second entrants are the first's neighbours.
E. The null law and its calibration across the grid offset (referee round 55): at the centroids OFFSETS_E within one grid cell and the
   cells CELLS_E near Δ₅₀, R0_E draws under H0 and R1_E under H1; after the first fit the duals move at ⟨r₁φ_j⟩, r₁ = Rξ to first order
   (R the projection off the fit's tangent span(φ, w∂_cφ)), so the null delay is min_j (1 − sign(n_j)u_j)/|n_j|, n ~ N(0, s²ΦᵀRΦ/N²),
   over the whole dictionary (`knot_law`). Its quantile at the truth against the H0 draws'; the test calibrated per draw by a parametric
   bootstrap of that law at the draw's own first fit (ĉ, ŵ): its size and power at every offset, against the oracle-calibrated test, the
   fixed threshold of the offset CBAR, and the Neyman–Pearson bound.

Output: results/chaos/close_features.md and close_features_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D,E`).
"""
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('CF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('CF_WORKERS', '3'))
CBAR = 0.032
PAIRS_A = ((0.8, 0.6), (0.8, 0.8), (0.8, -0.4), (0.8, -0.8))
DELTAS_A = (0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625)
DELTAS_B = (0.2, 0.1, 0.05, 0.025, 0.0125)
PAIR_C, SIGMAS_C, R_C, ALPHA_C = (0.8, 0.6), (0.001, 0.003, 0.01, 0.03), 200, 0.05
R_D, RHS_D = 40, 400_000          # D: the placement draws per (Δ, s) under H1 and their budget of right-hand-side evaluations each
DELTAS_C = (0.3, 0.2, 0.14, 0.1, 0.07, 0.05, 0.035, 0.025, 0.0175, 0.0125)
OFFSETS_E, CELLS_E = (0.0, 0.02, 0.032, 0.045, 0.049), ((0.0175, 0.003), (0.025, 0.01), (0.05, 0.03))
R0_E, R1_E, M_E, M_LAW = 1000, 400, 2000, 20000   # E: draws under H0 and H1, the bootstrap per draw, the law's samples at the truth
N_BOOT = 2000                                    # C: bootstrap resamples for Δ₅₀'s intervals


def setup():
    warnings.filterwarnings('ignore')
    L = X.law()
    K = RT.Bumps(np.linspace(-1, 1, X.NS_C), X.SIG_C)
    return L, K, np.linspace(-1, 1, X.D_C)


def pair_centres(A1, A2, D, cbar=None):
    """The pair's centres: the centroid at cbar (default CBAR; A ≠ 0), the midpoint at cbar (a dipole)."""
    cb = CBAR if cbar is None else cbar
    A = A1 + A2
    if abs(A) < 1e-12:
        return np.array([cb - D / 2, cb + D / 2])
    return np.array([cb - A2 * D / A, cb + A1 * D / A])


def d2phi(K, c):
    c = np.atleast_1d(c)
    dx = K.xs[:, None] - c[None, :]
    return K.phi(c) * (dx ** 2 / K.sig ** 4 - 1 / K.sig ** 2)


def taylor_basis(K, c):
    """e_n = ∂ⁿ_cφ(c)/n!, n = 0 … 3, on the samples."""
    dx = K.xs - c
    ph, s2 = np.exp(-dx ** 2 / (2 * K.sig ** 2)), K.sig ** 2
    return [ph, ph * dx / s2, ph * (dx ** 2 / s2 ** 2 - 1 / s2) / 2, ph * (dx ** 3 / s2 ** 3 - 3 * dx / s2 ** 2) / 6]


def perp_norm2(vs, k):
    """⟨(P⊥e_k)²⟩, P⊥ the projection orthogonal to e_0 … e_{k−1} (sample means)."""
    B = np.array(vs[:k]).T
    e = vs[k]
    r = e - B @ np.linalg.lstsq(B, e, rcond=None)[0]
    return float(np.mean(r ** 2))


def hessian_prediction(K, A1, A2, D):
    """The block's two small Hessian eigenvalues to leading order: the split v = (1, −1, Δ/(2w₁), Δ/(2w₂)) leaves Jv = (Δ³/2)e₃,
    λ₁ = Δ⁶⟨(P⊥e₃)²⟩/8, amplitude-free; the separation (0, 0, 1/w₁, −1/w₂) leaves Jv = −2Δe₂, λ₂ = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²)."""
    cs = pair_centres(A1, A2, D)
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    return D ** 6 * perp_norm2(vs, 3) / 8, 4 * D ** 2 * perp_norm2(vs, 2) * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2)


def predicted_delay(K, f, c0, A1, A2, D, order=2, cbar=None):
    """The leading-order second entry: the duals at the first entry u_j = ⟨fφ_j⟩/|⟨fφ_{j₁}⟩| moving at v_j = ½V⟨ψφ_j⟩ (ψ at the centroid)
    to the wall ±1 they meet first; returns (Δτ₂, j₂, the wall's sign). order = 3 (referee round 55): the residual
    P⊥(½V∂²_cφ + ⅙T₃∂³_cφ), P⊥ off the first fit's tangent span(φ, ∂_cφ) at the centroid, T₃ = A₁A₂(A₁ − A₂)Δ³/A² — the O(ϑ) term."""
    nx = K.xs.size
    A = A1 + A2
    cb = (A1 * pair_centres(A1, A2, D, cbar)[0] + A2 * pair_centres(A1, A2, D, cbar)[1]) / A
    V = A1 * A2 * D ** 2 / A
    ph, p2 = K.phi(np.array([cb]))[:, 0], d2phi(K, cb)[:, 0]
    if order == 3:
        e = taylor_basis(K, cb)                                       # e_n = ∂ⁿ_cφ/n!
        T = np.c_[e[0], e[1]]
        r = V * e[2] + A1 * A2 * (A1 - A2) * D ** 3 / A ** 2 * e[3]   # ½V∂²φ + ⅙T₃∂³φ
        psi = (r - T @ np.linalg.lstsq(T, r, rcond=None)[0]) / (0.5 * V)
    else:
        psi = p2 - (p2 @ ph) / (ph @ ph) * ph
    g = K.phi(c0).T @ f / nx
    j1 = int(np.argmax(np.abs(g)))
    u = g / abs(g[j1])
    v = 0.5 * V * (K.phi(c0).T @ psi) / nx                            # the duals' rates after the first fit (ρ = 1, per unit τ)
    best = (np.inf, None, 0)
    for j in range(c0.size):
        if j == j1 or v[j] == 0:
            continue
        wall = 1.0 if v[j] > 0 else -1.0
        t = (wall - u[j]) / v[j]
        if 0 < t < best[0]:
            best = (float(t), j, int(wall))
    return best


def a_job(args):
    A1, A2, D = args
    L, K, c0 = setup()
    cs = pair_centres(A1, A2, D)
    f = X.bumps_target(K, cs, (A1, A2))
    t0 = time.time()
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    ent = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
    act = np.flatnonzero(R['act'])
    out = dict(A=(A1, A2), D=D, cs=cs.tolist(), entries=ent, n_act=int(act.size), c_end=np.sort(R['c'][act]).tolist(),
               w_end=R['w'][act][np.argsort(R['c'][act])].tolist(), rms=float(R['hist'][-1][3]), relaxed=all(R['relaxed']),
               secs=time.time() - t0)
    if abs(A1 + A2) > 1e-12:
        out['pred'] = predicted_delay(K, f, c0, A1, A2, D)
    else:                                                              # a dipole: the transit of its derivative shape, slowed by 1/Δ
        shape = -A1 * K.dphi(np.array([CBAR]))[:, 0]
        Rs = RT.sliding_iss(L, K, shape, c0, np.ones(c0.size), max_events=2, polish=True)
        es = [(float(t), int(j)) for t, j, k, c in Rs['events'] if k == 'enter']
        out['shape'] = dict(tau1=es[0][0] / D, dtau2=(es[1][0] - es[0][0]) / D, kernels=[j for t, j in es])
    return out


def part_A():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, [(a1, a2, D) for a1, a2 in PAIRS_A for D in DELTAS_A])


def b_job(args):
    A1, A2, D = args
    L, K, c0 = setup()
    nx = K.xs.size
    cs = pair_centres(A1, A2, D)
    w = np.array([A1, A2])
    J = np.c_[K.phi(cs), K.dphi(cs) * w[None, :]] / np.sqrt(nx)
    hess = np.sort(np.linalg.eigvalsh(J.T @ J))
    f = X.bumps_target(K, cs, (A1, A2))
    rhs = RT.critical_rhs(L, K, f)
    al = np.sign(w) * RT._alpha_of(L, np.abs(w), RT.coef_crit)
    z = np.r_[al, cs]
    Jf = np.zeros((4, 4))
    for k in range(4):
        h = 1e-6 * max(1.0, abs(z[k]))
        e = np.zeros(4); e[k] = h
        Jf[:, k] = (rhs(0, z + e) - rhs(0, z - e)) / (2 * h)
    ev, V = np.linalg.eig(Jf)
    rates = np.sort(-np.real(ev))
    l1, l2 = hessian_prediction(K, A1, A2, D)
    # the flow's split and separation modes, identified by their eigenvectors in z = (α₁, α₂, c₁, c₂), against their leading order:
    # ẇ = −μ∂E/∂w (μ = w′(α)·curve_mobility, the amplitudes' mobility), ċ = −Gm∂E/∂c; the split's rate λ₁·2/(1/μ₁ + 1/μ₂) (its
    # Rayleigh quotient in the mobility's metric, the centres' share O(Δ²)); the separation's at frozen amplitudes (Gm Δ² ≫ μ)
    # Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) — the projection orthogonal to the translation only.
    h_ = 1e-7
    wp = np.array([(RT.coef_crit(L, abs(a) + h_) - RT.coef_crit(L, abs(a) - h_)) / (2 * h_) for a in al])
    mu = wp * RT.curve_mobility(L, al)
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    B1 = vs[1][:, None]
    e2p = vs[2] - B1 @ np.linalg.lstsq(B1, vs[2], rcond=None)[0]
    pred_split = l1 * 2 / (1 / mu[0] + 1 / mu[1])
    pred_sep = L.G * L.m * 4 * D ** 2 * float(np.mean(e2p ** 2)) * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2)

    def mode(s_):
        Vr = np.real(V)
        cos = np.abs(Vr.T @ s_) / (np.linalg.norm(Vr, axis=0) * np.linalg.norm(s_))
        k = int(np.argmax(cos))
        return float(-np.real(ev[k])), float(cos[k]), Vr[:, k]
    split = mode(np.r_[1 / wp[0], -1 / wp[1], D / (2 * A1), D / (2 * A2)])
    sep = mode(np.r_[0.0, 0.0, 1 / A1, -1 / A2])
    v = split[2]
    a = (v[0] * wp[0] - v[1] * wp[1]) / 2                             # the weight moved, w = ±a; the centres against Δ/(2w_k)·a
    comp = [float(v[2] / (D / (2 * A1) * a)), float(v[3] / (D / (2 * A2) * a))]
    return dict(A=(A1, A2), D=D, hess=hess.tolist(), rates=rates.tolist(), pred=(l1, l2), mu=mu.tolist(),
                split=split[:2], sep=sep[:2], pred_flow=(pred_split, pred_sep), split_c=comp)


DELTAS_SEP = (0.025, 0.0125, 0.00625, 0.003, 0.0015, 0.0005)


def sep_regimes(A1, A2, D):
    """(ii)'s separation rate in its two regimes (referee round 55): the exact linearisation −D·JᵀJ at the fit (D = diag(μ₁, μ₂, Gm, Gm)
    in (w, c), μ_k the amplitudes' mobility), its rate on the mode nearest the separation direction (0, 0, 1/w₁, −1/w₂), over the
    frozen-amplitude rate Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) and over the relaxed-amplitude one Gm·λ_sep (P⊥ off e₀ and e₁)."""
    L, K, c0 = setup()
    cs = pair_centres(A1, A2, D)
    w = np.array([A1, A2])
    J = np.c_[K.phi(cs), K.dphi(cs) * w[None, :]] / np.sqrt(K.xs.size)
    al = np.sign(w) * RT._alpha_of(L, np.abs(w), RT.coef_crit)
    h_ = 1e-7
    wp = np.array([(RT.coef_crit(L, abs(a) + h_) - RT.coef_crit(L, abs(a) - h_)) / (2 * h_) for a in al])
    d = np.r_[wp * RT.curve_mobility(L, al), L.G * L.m, L.G * L.m]
    lam, U = np.linalg.eigh(np.sqrt(d)[:, None] * (J.T @ J) * np.sqrt(d)[None, :])
    V = np.sqrt(d)[:, None] * U                                       # the eigenvectors of D·JᵀJ
    sv = np.r_[0.0, 0.0, 1 / A1, -1 / A2]
    k = int(np.argmax(np.abs(V.T @ sv) / np.linalg.norm(V, axis=0)))
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    B1 = vs[1][:, None]
    e2p = vs[2] - B1 @ np.linalg.lstsq(B1, vs[2], rcond=None)[0]
    wf = 4 * D ** 2 * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2) * L.G * L.m
    return float(lam[k] / (wf * np.mean(e2p ** 2))), float(lam[k] / (wf * perp_norm2(vs, 2))), float(np.min(d[:2] * np.diag(J.T @ J)[:2]))


def part_B():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(a1, a2, D) for a1, a2 in PAIRS_A for D in DELTAS_B])


def one_fit(K, f):
    """The single-kernel least-squares fit of f (the closest one-bump target): returns its centre and amplitude."""
    nx = K.xs.size
    g = K.phi(np.linspace(-1, 1, 201)).T @ f / nx
    c = np.linspace(-1, 1, 201)[int(np.argmax(np.abs(g)))]
    w = 0.0
    for _ in range(100):
        ph = K.phi(np.array([c]))[:, 0]
        w = (ph @ f) / (ph @ ph)
        dp = K.dphi(np.array([c]))[:, 0]
        r = f - w * ph
        gr = w * (dp @ r)
        hs = w * w * (dp @ dp)
        c = c + gr / hs
    return float(c), float(w)


def tangent_projector(K, c, w):
    """R = I − T(TᵀT)⁻¹Tᵀ, T = [φ(c), w∂_cφ(c)]: the projection off the one-kernel fit's tangent span."""
    T = np.c_[K.phi(np.array([c]))[:, 0], w * K.dphi(np.array([c]))[:, 0]]
    return np.eye(K.xs.size) - T @ np.linalg.solve(T.T @ T, T.T)


def knot_of(K, c0, c, w, Xi):
    """The second knot of the one-bump target wφ(c) under the noise rows Xi, first order in the noise, in the sliding limit (referee round
    55): the first entrant j₁ the largest |⟨yφ_j⟩| (equal masses: the duals at its entry u_j = ⟨yφ_j⟩/|⟨yφ_{j₁}⟩|), then the duals' rates
    n_j = ⟨Rξ, φ_j⟩ (r₁ = Rξ) to their walls: Δτ₂ = min_{j ≠ j₁} (1 − sign(n_j)u_j)/|n_j|."""
    N, M = K.xs.size, Xi.shape[0]
    Phi = K.phi(c0)
    v0 = (w * K.phi(np.array([c]))[:, 0][None, :] + Xi) @ Phi / N
    j1 = np.argmax(np.abs(v0), axis=1)
    u = v0 / np.abs(v0[np.arange(M), j1])[:, None]
    n = Xi @ (tangent_projector(K, c, w) @ Phi) / N
    with np.errstate(divide='ignore'):
        t = (1 - np.sign(n) * u) / np.abs(n)
    t[np.arange(M), j1] = np.inf
    return t.min(axis=1)


def knot_law(K, c0, c, w, s, M, rng):
    """The second knot's null law at the one-bump fit (c, w): `knot_of` on M draws of i.i.d. N(0, s²) noise."""
    return knot_of(K, c0, c, w, s * rng.standard_normal((M, K.xs.size)))


def np_power(K, f1, f0, s, alpha=0.05):
    """The Neyman–Pearson power of the simple hypotheses f₀ against f₁ under i.i.d. N(0, s²) labels: Φ(‖f₁ − f₀‖/s − z_{1−α})."""
    from scipy.stats import norm
    return float(norm.cdf(np.linalg.norm(f1 - f0) / s - norm.ppf(1 - alpha)))


def e_job(args):
    """E: at the centroid cbar, R0_E draws under H0 (the one-bump fit of the pair) and R1_E under H1: the second knot Δτ₂, the entrants,
    the first fit (ĉ, ŵ) and the knot's p-value under the law bootstrapped at that fit (M_E draws); the law at the truth (M_LAW draws),
    the Neyman–Pearson power and the projected noise's std over the unprojected one per kernel, (⟨(Rφ_j)²⟩/⟨φ_j²⟩)^{1/2}."""
    D, s, cbar, seed = args
    t0 = time.time()
    L, K, c0 = setup()
    N = K.xs.size
    f1 = X.bumps_target(K, pair_centres(*PAIR_C, D, cbar), PAIR_C)
    c1, w1 = one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    rows = {}
    for hyp, f, R_ in ((0, f0, R0_E), (1, f1, R1_E)):
        out = []
        for r in range(R_):
            y = f + s * rng.standard_normal(N)
            R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
            ent = [(t, j) for t, j, k, c in R['events'] if k == 'enter']
            dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.inf
            h = R['hist'][1]
            a = int(np.flatnonzero(h[2])[0])
            law = knot_law(K, c0, float(h[1][a]), float(h[2][a]), s, M_E, np.random.default_rng([seed, hyp, r]))
            out.append(dict(dt=float(dt), j=[j for t, j in ent], c1=float(h[1][a]), w1=float(h[2][a]), p=float(np.mean(law <= dt))))
        rows[hyp] = out
    law0 = knot_law(K, c0, c1, w1, s, M_LAW, np.random.default_rng([seed, 9]))
    Rm = tangent_projector(K, c1, w1)
    Phi = K.phi(c0)
    ratio = np.sqrt(np.mean((Rm @ Phi) ** 2, axis=0) / np.mean(Phi ** 2, axis=0))
    print(f"  E Δ {D:g} s {s:g} c̄ {cbar:g}: H0 q·s {np.quantile([q['dt'] for q in rows[0]], ALPHA_C) * s:.3f} (law "
          f"{np.quantile(law0, ALPHA_C) * s:.3f}) [{time.time() - t0:.0f}s]", flush=True)
    return dict(D=D, s=s, cbar=cbar, rows0=rows[0], rows1=rows[1], q_law=float(np.quantile(law0, ALPHA_C)),
                np=np_power(K, f1, f0, s), ratio=ratio.tolist(), c1=c1, w1=w1)


def part_E():
    from multiprocessing import Pool
    jobs = [(D, s, cb, 55000 + 100 * i + k) for i, (D, s) in enumerate(CELLS_E) for k, cb in enumerate(OFFSETS_E)]
    with Pool(N_WORKERS) as pool:
        return pool.map(e_job, jobs, chunksize=1)


def c_job(args):
    """C and D: R_C draws of y = f + ξ under H0 (the one-bump fit of the pair) or H1 (the pair): the path to its second entry, not
    relaxed after it (`max_entries` = 2: the delay Δτ₂ only); under H1 the first R_D draws once more to the two-kernel block's relaxation
    within RHS_D right-hand-side evaluations — its centres and coefficients if it relaxed, 'oob' if the budget ran out first (a block
    near coalescence relaxes at a rate ∝ Δ⁶, Part B; unbounded, such relaxations exhausted the memory)."""
    D, s, hyp, seed = args
    t0 = time.time()
    L, K, c0 = setup()
    cs = pair_centres(*PAIR_C, D)
    f1 = X.bumps_target(K, cs, PAIR_C)
    if hyp == 1:
        f = f1
    else:
        c1, w1 = one_fit(K, f1)
        f = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    out = []
    for r in range(R_C):
        y = f + s * rng.standard_normal(K.xs.size)
        R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [(t, j) for t, j, k, c in R['events'] if k == 'enter']
        dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.inf
        row = dict(dt=float(dt), j=[j for t, j in ent])
        if hyp == 1 and r < R_D and len(ent) > 1:
            R2 = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_rhs=RHS_D)
            act = np.flatnonzero(R2['act'])                          # the block after the second entry's relaxation
            row['oob'] = bool(R2['out_of_budget'])
            if not R2['out_of_budget'] and len(R2['relaxed']) > 1 and R2['relaxed'][1] and act.size == 2:
                o = np.argsort(R2['c'][act])
                row.update(c=R2['c'][act][o].tolist(), w=R2['w'][act][o].tolist())
        out.append(row)
    print(f"  C Δ {D:g} s {s:g} H{hyp}: second entries {sum(np.isfinite(q['dt']) for q in out)}/{len(out)}, placed "
          f"{sum('c' in q for q in out)}, out of budget {sum(q.get('oob', False) for q in out)} [{time.time() - t0:.0f}s]", flush=True)
    return dict(D=D, s=s, hyp=hyp, rows=out)


def part_C():
    from multiprocessing import Pool
    jobs = [(D, s, h, 1000 * k + 10 * h + i) for i, (D, s) in enumerate([(D, s) for D in DELTAS_C for s in SIGMAS_C])
            for h, k in ((0, 1), (1, 2))]
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, jobs, chunksize=1)


def crossing(Ds, pw, level=0.5):
    """The separation where the power, increasing in Δ, first crosses `level` (log-linear interpolation), NaN if it does not."""
    for k in range(len(Ds) - 1):
        if pw[k] < level <= pw[k + 1]:
            t = (level - pw[k]) / (pw[k + 1] - pw[k])
            return float(np.exp(np.log(Ds[k]) + t * (np.log(Ds[k + 1]) - np.log(Ds[k]))))
    return np.nan


def detection(C):
    """Per (Δ, s): the H0 quantile of Δτ₂ at ALPHA_C and the H1 power; Δ₅₀(s) by log-linear interpolation of the power."""
    tab = {}
    for r in C:
        tab.setdefault((r['D'], r['s']), {})[r['hyp']] = np.array([q['dt'] for q in r['rows']])
    out = {}
    for (D, s), d in tab.items():
        q = np.quantile(d[0], ALPHA_C)
        out[(D, s)] = dict(q0=float(q), power=float(np.mean(d[1] < q)), med1=float(np.median(d[1])), med0=float(np.median(d[0])))
    Ds = sorted(DELTAS_C)
    d50 = {s: crossing(Ds, [out[(D, s)]['power'] for D in Ds]) for s in SIGMAS_C}
    return out, d50


def d50_bootstrap(C, n=N_BOOT, seed=5500):
    """Δ₅₀'s and its slope's 95% intervals: both hypotheses' draws resampled per cell (the H0 quantile recomputed each time)."""
    tab = {}
    for r in C:
        tab.setdefault((r['D'], r['s']), {})[r['hyp']] = np.array([q['dt'] for q in r['rows']])
    rng = np.random.default_rng(seed)
    Ds = sorted(DELTAS_C)
    pw = {}
    for (D, s), d in sorted(tab.items()):
        q = np.quantile(d[0][rng.integers(0, d[0].size, (n, d[0].size))], ALPHA_C, axis=1)
        pw[(D, s)] = np.mean(d[1][rng.integers(0, d[1].size, (n, d[1].size))] < q[:, None], axis=1)
    d50 = {s: np.array([crossing(Ds, [pw[(D, s)][b] for D in Ds]) for b in range(n)]) for s in SIGMAS_C}
    sl = []
    for b in range(n):
        x = [(np.log(s), np.log(d50[s][b])) for s in SIGMAS_C if np.isfinite(d50[s][b])]
        sl.append(np.polyfit(*np.array(x).T, 1)[0] if len(x) >= 2 else np.nan)
    ci = lambda v: ((float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))) if np.isfinite(v).any() else (np.nan, np.nan)) + \
        (int(np.sum(~np.isfinite(v))),)
    return {s: ci(d50[s]) for s in SIGMAS_C}, ci(np.array(sl))


def crb(K, cs, w, s):
    """The Cramér–Rao bounds of the two-kernel fit at the truth, s²(JᵀJ)⁻¹ (J the samples' Jacobian in (w₁, w₂, c₁, c₂)): the positions'
    rms ((Var c₁ + Var c₂)/2)^{1/2} and the split w₁/(w₁ + w₂)'s std."""
    J = np.c_[K.phi(cs), K.dphi(cs) * w[None, :]]
    Cov = s ** 2 * np.linalg.inv(J.T @ J)
    W = w.sum()
    g = np.array([w[1] / W ** 2, -w[0] / W ** 2, 0.0, 0.0])
    return float(np.sqrt((Cov[2, 2] + Cov[3, 3]) / 2)), float(np.sqrt(g @ Cov @ g))


def report(A=None, B=None, C=None, E_=None):
    fmt = E.fmt
    L_ = ['# Proposition D‴.52: close features (`exp_close_features.py`)', '',
          f"The sliding inverse scale space (polished) on two bumps of the kernels' width {X.SIG_C} a distance Δ apart, the centroid (the "
          f"midpoint for a dipole) at {CBAR}; 61 samples of [−1, 1], 21 kernels on a grid of spacing 0.1, equal masses.", '']
    if A is not None:
        L_ += ['## A. The second entry', '',
               '| A₁, A₂ | Δ | entries (τ, kernel) | Δτ₂ | predicted (kernel, wall) | ratio | with T₃ | Δτ₂Δ² | Δτ₂Δ | the end: kernels, centres '
               '(truth), rms |', '|---|---|---|---|---|---|---|---|---|---|']
        La, Ka, c0a = setup()
        for r in A:
            ent = r['entries']
            dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.nan
            p, sh = r.get('pred'), r.get('shape')
            r3 = '—'
            if p:
                pr, ratio = f"{p[0]:.4g} ({p[1]}, {p[2]:+d})", f"{dt / p[0]:.4f}"
                p3 = predicted_delay(Ka, X.bumps_target(Ka, pair_centres(*r['A'], r['D']), r['A']), c0a, *r['A'], r['D'], order=3)
                r3 = f"{dt / p3[0]:.4f}" + (f" ({p3[1]}, {p3[2]:+d})" if p3[1] != p[1] else '')
            elif sh:
                pr = f"shape: τ₁ {sh['tau1']:.4g}, Δτ₂ {sh['dtau2']:.4g} ({sh['kernels']})"
                ratio = f"{ent[0][0] / sh['tau1']:.4f}, {dt / sh['dtau2']:.4f}"
            else:
                pr = ratio = '—'
            L_.append(f"| {r['A'][0]:g}, {r['A'][1]:g} | {r['D']:g} | {', '.join(f'({t:.4g}, {j})' for t, j in ent[:3])} | {dt:.4g} | {pr} | "
                      f"{ratio} | {r3} | {dt * r['D'] ** 2:.4g} | {dt * r['D']:.4g} | {r['n_act']}, {[round(v, 5) for v in r['c_end']]} "
                      f"({[round(v, 5) for v in r['cs']]}), {fmt(r['rms'])} |")
        L_.append('')
    if B is not None:
        L_ += ['## B. The two-kernel block at its exact fit', '',
               'The Gauss–Newton Hessian JᵀJ (J = ∂(Φ(c)w)/∂(w, c)/N^{1/2}) and the fast flow on the critical curve linearised at the fit '
               '(its relaxation rates): eigenvalues, the two smallest scaled by Δ² and Δ⁶ against their leading order (the split mode: '
               'λ₁ = Δ⁶⟨(P⊥e₃)²⟩/8, amplitude-free; the separation: λ₂ = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²); e_n = ∂ⁿ_cφ/n! at the pair, P⊥ '
               'orthogonal to the lower e_m), and the local slopes of the smallest two in log Δ. The flow\'s split and separation modes are '
               'identified by their eigenvectors\' overlap with (1/w′₁, −1/w′₂, Δ/(2w₁), Δ/(2w₂)) and (0, 0, 1/w₁, −1/w₂) in (α, c), and '
               'compared with λ₁·2/(1/μ₁ + 1/μ₂) (μ_k the amplitudes\' mobility, ẇ = −μ∂E/∂w) and, at frozen amplitudes, '
               'Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²).', '',
               '| A₁, A₂ | Δ | Hessian eigenvalues | λ₂/Δ², λ₁/Δ⁶ | predicted | slopes | flow rates | split: r/Δ⁶ (predicted; '
               'the mode\'s overlap; its centre components over Δ/(2w_k) per unit weight moved) | separation: r/Δ² (predicted; overlap) '
               '| slopes |',
               '|---|---|---|---|---|---|---|---|---|---|']
        prev = {}
        for r in B:
            key = tuple(r['A'])
            h, q = r['hess'], r['rates']
            sl = sq = '—'
            sp, se = r.get('split'), r.get('sep')
            fq = (sp[0], se[0]) if sp else (q[0], q[1])
            if key in prev:
                pD, ph, pq = prev[key]
                sl = ', '.join(f"{np.log(h[k] / ph[k]) / np.log(r['D'] / pD):.3f}" for k in (0, 1))
                sq = ', '.join(f"{np.log(fq[k] / pq[k]) / np.log(r['D'] / pD):.3f}" for k in (0, 1))
            prev[key] = (r['D'], h, fq)
            pl, pf = r.get('pred'), r.get('pred_flow')
            pt = f"{pl[1] / r['D'] ** 2:.4g}, {pl[0] / r['D'] ** 6:.4g}" if pl else '—'
            D6, D2 = r['D'] ** 6, r['D'] ** 2
            if sp:
                sc = r.get('split_c')
                fs = (f"{sp[0] / D6:.4g} ({pf[0] / D6:.4g}; {sp[1]:.3f}" + (f"; centres {sc[0]:.4f}, {sc[1]:.4f}" if sc is not None else '')
                      + f") | {se[0] / D2:.4g} ({pf[1] / D2:.4g}; {se[1]:.3f})")
            else:
                fs = f"{q[0] / D6:.4g} | {q[1] / D2:.4g}"
            L_.append(f"| {key[0]:g}, {key[1]:g} | {r['D']:g} | {', '.join(f'{v:.3e}' for v in h)} | {h[1] / D2:.4g}, "
                      f"{h[0] / D6:.4g} | {pt} | {sl} | {', '.join(f'{v:.3e}' for v in q)} | {fs} | {sq} |")
        L_ += ['', "The separation's rate in its two regimes (the exact linearisation −D·JᵀJ at the fit, the mode nearest the separation): "
               "over the frozen-amplitude rate Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) and over the relaxed-amplitude rate Gm·λ_sep "
               "(P⊥ off e₀ and e₁) — the first while the separation outpaces the amplitudes' mode μ_k⟨φ_k²⟩, the second once it is slower:", '',
               '| A₁, A₂ | ' + ' | '.join(f'Δ = {D:g}' for D in DELTAS_SEP) + ' | the amplitudes\' mode min μ_k⟨φ_k²⟩ |',
               '|---|' + '---|' * (len(DELTAS_SEP) + 1)]
        for A1, A2 in [p_ for p_ in PAIRS_A if abs(p_[0] + p_[1]) > 1e-12]:
            sr = [sep_regimes(A1, A2, D) for D in DELTAS_SEP]
            L_.append(f"| {A1:g}, {A2:g} | " + ' | '.join(f"{a:.3f}, {b:.3f}" for a, b, m in sr) + f" | {sr[0][2]:.3g} |")
        L_.append('')
    if C is not None:
        det, d50 = detection(C)
        L_ += ['## C. Number detection under label noise', '',
               f"A = {PAIR_C}; {R_C} draws per hypothesis; H0 the closest one-bump target; the power of Δτ₂ < its H0 quantile "
               f"{ALPHA_C:g} at the true offset (an oracle-calibrated simple null; E calibrates it), in brackets the Neyman–Pearson power "
               "Φ(‖f₁ − f₀‖/s − z₀.₉₅) of the same simple hypotheses (the medians of Δτ₂ under H1 and H0):", '',
               '| Δ | ' + ' | '.join(f's = {s:g}' for s in SIGMAS_C) + ' |', '|---|' + '---|' * len(SIGMAS_C)]
        Lc, Kc, c0c = setup()
        npw = {}
        for D in DELTAS_C:
            f1 = X.bumps_target(Kc, pair_centres(*PAIR_C, D), PAIR_C)
            c1_, w1_ = one_fit(Kc, f1)
            for s_ in SIGMAS_C:
                npw[(D, s_)] = np_power(Kc, f1, w1_ * Kc.phi(np.array([c1_]))[:, 0], s_, ALPHA_C)
        for D in DELTAS_C:
            L_.append(f"| {D:g} | " + ' | '.join(f"{det[(D, s)]['power']:.2f} [{npw[(D, s)]:.2f}] ({det[(D, s)]['med1']:.3g}, "
                                                f"{det[(D, s)]['med0']:.3g})" for s in SIGMAS_C) + ' |')
        ss = np.array([s for s in SIGMAS_C if np.isfinite(d50[s])])
        dd = np.array([d50[s] for s in ss])
        sl = np.polyfit(np.log(ss), np.log(dd), 1)[0] if ss.size >= 2 else np.nan
        # the prediction: under H1 the second entry comes after ≈ C₁/Δ² (C₁ the leading order's Δτ₂Δ² at the smallest Δ, (i)), so the
        # power is ½ where C₁/Δ² meets the H0 quantile q_α (∝ √N/s): Δ₅₀ = (C₁/q_α)^{1/2}, q_α interpolated in log Δ at the crossing
        Dm = min(DELTAS_C)
        C1 = predicted_delay(Kc, X.bumps_target(Kc, pair_centres(*PAIR_C, Dm), PAIR_C), c0c, *PAIR_C, Dm)[0] * Dm ** 2
        pred50 = {}
        for s_ in SIGMAS_C:
            Ds = np.array(sorted(DELTAS_C))
            gap = np.log(C1 / Ds ** 2) - np.log([det[(D, s_)]['q0'] for D in Ds])
            pred50[s_] = np.nan
            for k in range(len(Ds) - 1):
                if gap[k] > 0 >= gap[k + 1] or gap[k] <= 0 < gap[k + 1]:
                    t = gap[k] / (gap[k] - gap[k + 1])
                    pred50[s_] = float(np.exp(np.log(Ds[k]) + t * (np.log(Ds[k + 1]) - np.log(Ds[k]))))
                    break
        qs = {s_: float(np.median([det[(D, s_)]['q0'] for D in DELTAS_C]) * s_) for s_ in SIGMAS_C}
        ci50, cisl = d50_bootstrap(C)
        h0 = [(q['dt'], abs(q['j'][1] - q['j'][0]) <= 1) for r in C if r['hyp'] == 0 for q in r['rows'] if len(q['j']) > 1]
        early = [a for r in C if r['hyp'] == 0 for q, a in [(q, abs(q['j'][1] - q['j'][0]) <= 1) for q in r['rows'] if len(q['j']) > 1]
                 if q['dt'] <= np.quantile([q_['dt'] for q_ in r['rows']], ALPHA_C)]
        npd = {s_: crossing(sorted(DELTAS_C), [npw[(D, s_)] for D in sorted(DELTAS_C)]) for s_ in SIGMAS_C}
        L_ += ['', 'Δ₅₀ (the power ½): ' + ', '.join(f's = {s:g}: {fmt(d50[s])}' for s in SIGMAS_C) +
               f"; the slope of log Δ₅₀ in log s {sl:.3f} (½ predicted). Predicted from (i) and the H0 quantile, Δ₅₀ = (C₁/q_α)^{{1/2}} "
               f"(C₁ = {C1:.4f}, the leading order's Δτ₂Δ² at Δ = {Dm:g}): " + ', '.join(f's = {s:g}: {fmt(pred50[s])}' for s in SIGMAS_C) +
               "; the H0 quantile times s (its median over Δ; constant if q_α ∝ 1/s): " + ', '.join(f'{qs[s]:.4g}' for s in SIGMAS_C) + '. '
               f"The bootstrap's 95% intervals ({N_BOOT} resamples of both hypotheses' draws per cell): Δ₅₀ "
               + ', '.join(f"s = {s:g}: " + ('— (no crossing)' if ci50[s][2] == N_BOOT else f"[{ci50[s][0]:.4g}, {ci50[s][1]:.4g}]"
                                                 + (f" ({ci50[s][2]} resamples without a crossing)" if ci50[s][2] else ''))
                           for s in SIGMAS_C)
               + f"; the slope [{cisl[0]:.3f}, {cisl[1]:.3f}]. The Neyman–Pearson bound's Δ₅₀: "
               + ', '.join(f's = {s:g}: {fmt(npd[s])}' for s in SIGMAS_C) + '. '
               f"Under H0 the second entrant is the first's neighbour in {np.mean([a for t, a in h0]):.1%} of the {len(h0)} draws, in "
               f"{np.mean(early):.1%} of the earliest {ALPHA_C:.0%} of each cell (which set the quantile).", '']
        L_ += ['## D. Placement under label noise (the two-kernel block after the second entry, H1)', '',
               f"The first {R_D} H1 draws of each (Δ, s) relaxed after the second entry within {RHS_D} right-hand-side evaluations (a "
               "block near coalescence relaxes at a rate ∝ Δ⁶, B): the placed blocks (two kernels, relaxed) and those out of budget; "
               "the rms errors over the placed ones, against the Cramér–Rao bounds at the truth (s²(JᵀJ)⁻¹, `crb`); the placed blocks "
               "with both kernels within σ of their features, and the draws whose second entrant is the first's neighbour:", '',
               '| Δ | s | placed, out of budget (of the draws) | rms position error | × Δ²/s | / CRB | rms error of w₁/(w₁ + w₂) | × Δ³/s | '
               '/ CRB | both within σ | second entrant adjacent |',
               '|---|---|---|---|---|---|---|---|---|---|---|']
        cs_truth = {D: pair_centres(*PAIR_C, D) for D in DELTAS_C}
        for r in C:
            if r['hyp'] != 1:
                continue
            ok = [q for q in r['rows'] if 'c' in q]
            tried = [q for q in r['rows'] if 'oob' in q]
            adj = sum(len(q['j']) > 1 and abs(q['j'][1] - q['j'][0]) <= 1 for q in r['rows'][:R_D])
            if not ok:
                if tried:
                    L_.append(f"| {r['D']:g} | {r['s']:g} | 0, {sum(q['oob'] for q in tried)} ({len(tried)}) | — | — | — | — | — | — | — | "
                              f"{adj} |")
                continue
            pe = np.sqrt(np.mean([np.mean((np.array(q['c']) - cs_truth[r['D']]) ** 2) for q in ok]))
            fr = np.sqrt(np.mean([(q['w'][0] / (q['w'][0] + q['w'][1]) - PAIR_C[0] / sum(PAIR_C)) ** 2 for q in ok]))
            cp, cf_ = crb(Kc, cs_truth[r['D']], np.array(PAIR_C), r['s'])
            both = sum(bool(np.all(np.abs(np.array(q['c']) - cs_truth[r['D']]) < X.SIG_C)) for q in ok)
            L_.append(f"| {r['D']:g} | {r['s']:g} | {len(ok)}, {sum(q['oob'] for q in tried)} ({len(tried)}) | {pe:.3e} | "
                      f"{pe * r['D'] ** 2 / r['s']:.3g} | {pe / cp:.3g} | {fr:.3e} | {fr * r['D'] ** 3 / r['s']:.3g} | {fr / cf_:.3g} | "
                      f"{both} | {adj} |")
        L_.append('')
        D0 = min(DELTAS_C)
        cp, cf_ = crb(Kc, pair_centres(*PAIR_C, D0), np.array(PAIR_C), 1.0)
        good, big = [], None
        for r in C:
            ok = [q for q in r['rows'] if 'c' in q] if r['hyp'] == 1 else []
            if not ok:
                continue
            cs_ = cs_truth[r['D']]
            cp_, cf2 = crb(Kc, cs_, np.array(PAIR_C), r['s'])
            if cf2 <= 0.2:
                pe = np.sqrt(np.mean([np.mean((np.array(q['c']) - cs_) ** 2) for q in ok]))
                fr = np.sqrt(np.mean([(q['w'][0] / (q['w'][0] + q['w'][1]) - PAIR_C[0] / sum(PAIR_C)) ** 2 for q in ok]))
                good.append((pe / cp_, fr / cf2))
            for q in ok:
                if big is None or max(np.abs(q['w'])) > max(np.abs(big[2]['w'])):
                    big = (r['D'], r['s'], q)
        L_ += [f"The Cramér–Rao constants at Δ = {D0:g}: the positions' rms × Δ²/s {cp * D0 ** 2:.4f}, the split's × Δ³/s "
               f"{cf_ * D0 ** 3:.4f} — its std reaches 0.1 at s ≈ {0.1 / (cf_ * D0 ** 3):.1f}Δ³. In the {len(good)} cells where the split's "
               f"bound is at most 0.2 the measured errors are {min(g[0] for g in good):.2f}–{max(g[0] for g in good):.2f} (positions) and "
               f"{min(g[1] for g in good):.2f}–{max(g[1] for g in good):.2f} (split) of the bounds. The largest coefficient of a placed "
               f"block: Δ = {big[0]:g}, s = {big[1]:g}, w = ({big[2]['w'][0]:.4g}, {big[2]['w'][1]:.4g}) (their sum "
               f"{sum(big[2]['w']):.4g}), the centres {abs(big[2]['c'][1] - big[2]['c'][0]):.2g} apart, the entrants {big[2]['j']}.", '']
    if E_ is not None:
        L_ += ['## E. The null law and its calibration across the grid offset', '',
               f"At the centroids c̄ = {', '.join(f'{c:g}' for c in OFFSETS_E)} (the kernels on a grid of spacing 0.1: c̄ = 0.05 is a tie) and "
               f"the cells (Δ, s) = {', '.join(f'({D:g}, {s:g})' for D, s in CELLS_E)}, {R0_E} draws under H0 (the one-bump fit) and "
               f"{R1_E} under H1. The H0 quantile of Δτ₂ at {ALPHA_C:g} times s from the draws and from the null law at the truth "
               f"(`knot_law`, {M_LAW} samples: min_j (1 − sign(n_j)u_j)/|n_j|, n ~ N(0, s²ΦᵀRΦ/N²), R the projection off the fit's "
               "tangent span); the projected noise's std over the unprojected one, (⟨(Rφ_j)²⟩/⟨φ_j²⟩)^{1/2}, for the two kernels next to "
               "the fit; the share of the H0 second entrants next to the first; the false-alarm rate of the threshold set at "
               f"c̄ = {CBAR:g} and of the test calibrated per draw by the law bootstrapped at its own first fit ({M_E} samples); the power of "
               "the oracle-calibrated test (the offset's own H0 quantile), of the calibrated test, and the Neyman–Pearson bound:", '',
               '| Δ, s | c̄ | noiseless Δτ₂Δ² (leading order; kernel) | H0 quantile × s: draws, law | std ratio next to the fit | H0 second '
               f'entrants adjacent | false alarms: c̄ = {CBAR:g} threshold, calibrated | power: oracle, calibrated, Neyman–Pearson |',
               '|---|---|---|---|---|---|---|---|']
        Le, Ke, c0e = setup()
        for D, s_ in CELLS_E:
            es = sorted([e for e in E_ if e['D'] == D and e['s'] == s_], key=lambda e: e['cbar'])
            ref = next((e for e in es if abs(e['cbar'] - CBAR) < 1e-12), None)
            qref = np.quantile([q['dt'] for q in ref['rows0']], ALPHA_C) if ref is not None else np.nan
            for e in es:
                d0 = np.array([q['dt'] for q in e['rows0']]); d1 = np.array([q['dt'] for q in e['rows1']])
                p0 = np.array([q['p'] for q in e['rows0']]); p1 = np.array([q['p'] for q in e['rows1']])
                q0 = np.quantile(d0, ALPHA_C)
                cg = np.linspace(-1, 1, len(e['ratio']))
                nb = np.argsort(np.abs(cg - e['c1']))[:2]
                adj = np.mean([abs(q['j'][1] - q['j'][0]) <= 1 for q in e['rows0'] if len(q['j']) > 1])
                pd = predicted_delay(Ke, X.bumps_target(Ke, pair_centres(*PAIR_C, D, e['cbar']), PAIR_C), c0e, *PAIR_C, D,
                                     cbar=e['cbar'])
                L_.append(f"| {D:g}, {s_:g} | {e['cbar']:g} | {pd[0] * D ** 2:.3f} ({pd[1]}) | {q0 * s_:.3f}, {e['q_law'] * s_:.3f} | "
                          + ', '.join(f"{e['ratio'][k]:.3f}" for k in sorted(nb)) + f" | {adj:.2f} | {np.mean(d0 < qref):.3f}, "
                          f"{np.mean(p0 < ALPHA_C):.3f} | {np.mean(d1 < q0):.3f}, {np.mean(p1 < ALPHA_C):.3f}, {e['np']:.3f} |")
        L_.append('')
    L_.append('Numbers are produced by `experiments/exp_close_features.py`; nothing is edited by hand.')
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'close_features_raw.npz')
    md = os.path.join(OUT, 'close_features.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'E')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else {'A', 'B', 'C', 'D', 'E'}
    res = dict(old)
    t0 = time.time()
    for k, fn in (('B', part_B), ('A', part_A), ('C', part_C), ('E', part_E)):
        if k in parts or (k == 'C' and 'D' in parts):
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
