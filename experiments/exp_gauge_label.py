"""
Proposition D‴.41: measurement noise on the gauge — the selection is set by the noise's per-update variance
==========================================================================================================

Proposition D‴.36's instance (three bumps, one particle per class, m = 3). The measured error r + ξ drives both laws (the positions and
the masses), ξ independent over the quadrature nodes, of intensity v (⟨ξh⟩_ρ has the variance v⟨h²⟩_ρ per unit slow time): the noise
enters along the descent's own fields U_x = (β_iG∇_ϑΦ̃(x; ϑ_i), γ(Φ̃(x; ϑ_i) − f̂(x))) (`label_noise.fields`), the drift being ⟨rU⟩_ρ.

A. The identity at three points of M and γ = 10⁻², 10⁻⁵ (Katzenberger's limit with the covariance Σ = v⟨UUᵀ⟩_ρ at the point, in the flat
   coordinates of exp_gauge_drift.py): the Itô drift ½∂²Φ[Σ], the Wong–Zakai term (v/2)∂Φ⟨(U·∇)U⟩_ρ, their sum (the Stratonovich drift,
   zero by the proposition), the diffusion along M (zero), and the closed form on the masses in the slow-replicator limit,
   ū_k = −(v/2)β_kG_α⟨Φ_k²⟩_ρ(1 − (a_k/a_max)⁴), against the Itô drift.
B. The selection's target: P_k = G_α⟨Φ_k²⟩_ρ, the minimiser W_ℓ of Ψ_ℓ = (v/4)Σ_kP_k(β_ka′_k)² on the simplex (against Proposition
   D‴.38's W♯), and the closed form's replicator flow linearised there.
C. The full loop with measurement noise only (τ = 0, γ = 10⁻², v = V_C, from W₀ = (0.40, 0.25, 0.35), SEEDS_C × S_C slow units): the noise
   drawn afresh at each step (Euler, at the steps 0.25 and 0.125), AR(1) noise at ρ = −1/2, −1/3, 1/3, 0.6, 0.8 (the same intensity),
   and Heun's scheme (Stratonovich). The gauge coordinate y(s) = Φ(z(s)) (the limit map of the noiseless flow from the state, every
   REC_PHI slow units): its least-squares slope and its end displacement, against (1 − ρ)/(1 + ρ) times part A's Itô drift at W₀ (zero
   for Heun). (The raw masses' window means carry the O(v) shift of the normal fluctuations' mean, which builds up from the start on M.)
D. The example: the scalar factorisation f(u, v) = uv fitted to y = 1 with label noise, u̇ = (ξ − r)v, v̇ = (ξ − r)u, r = uv − y. The
   rescaling symmetry's Noether charge Q = u² − v² is conserved by the gradient flow and, pathwise, by every noise field (v, u) — so by the
   Stratonovich (Heun) dynamics —; an Euler step gives exactly Q_{n+1} = Q_n(1 − (ξ_n − r_n)²h²) whatever the noise sequence, so E[Q]
   decays at the rate set by the per-step variance: σ² for noise drawn afresh (the Itô SDE's exact rate), (1 − ρ)/(1 + ρ)σ² for AR(1)
   noise at the same intensity; N_PATHS_D paths, the decay rate −log(E[Q(T)]/Q(0))/(σ²T) against the factor.

Output: results/chaos/gauge_label.md and gauge_label_raw.npz (GL_WORKERS; `--from-raw` regenerates the md; `--only A,B,C,D`).
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
             ('ar1', 0.8, DT), ('heun', 0.0, DT))
SIG2_D, DT_D, T_D, N_PATHS_D, UV0_D, SEED_D = 0.01, 1e-3, 100.0, 400, (2.0, 0.5), 900
SCHEMES_D = (('ito', 0.0), ('ar1', -1 / 3), ('ar1', 1 / 3), ('ar1', 0.6), ('heun', 0.0))
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
    return dict(gamma=gamma, w=w, ito=np.r_[ito, -ito.sum()], wz=np.r_[wz, -wz.sum()], D=D, lam=lam, closed=closed,
                fit=float(np.abs(loop.residual(E.exact_rep(net, W), 3 * W)).max()))


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
    return dict(P=P, W_min=Wm, W_l=Wl, amp=np.abs(E.A_BUMP) / Wl / A_MAX, u_l=u_l, eig=np.linalg.eigvals(Jac),
                psi_l=LN.psi(Wl, P, Wm), psi_u=LN.psi(np.ones(3) / 3, P, Wm), psi_0=LN.psi(np.asarray(W0_C), P, Wm), W_sharp=Q.w_sharp())


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
    rho_ = rho if scheme == 'ar1' else 0.0
    sx = np.sqrt(SIG2_D / h * (1 - rho_) / (1 + rho_))
    xi = sx * rng.standard_normal(N_PATHS_D)
    Q0 = UV0_D[0] ** 2 - UV0_D[1] ** 2
    every = int(round(10.0 / h)); rec = []
    for k in range(n + 1):
        if k % every == 0:
            Q = u * u - v * v
            rec.append((k * h, Q.mean(), Q.std(), np.abs(Q - Q0).max(), np.abs(u * v - 1).mean()))
        if k == n:
            break
        if k > 0:
            xi = rho_ * xi + np.sqrt(1 - rho_ ** 2) * sx * rng.standard_normal(N_PATHS_D)
        c = xi - (u * v - 1.0)
        if scheme == 'heun':
            up, vp = u + h * c * v, v + h * c * u
            cp = xi - (up * vp - 1.0)
            u, v = u + 0.5 * h * (c * v + cp * vp), v + 0.5 * h * (c * u + cp * up)
        else:
            u, v = u + h * c * v, v + h * c * u
    rec = np.array(rec)
    return dict(scheme=scheme, rho=rho, rec=rec, rate=-np.log(rec[-1, 1] / Q0) / (SIG2_D * T_D))


def name_C(scheme, rho, dt):
    out = {'ito': 'afresh each step', 'heun': 'Heun (Stratonovich)'}.get(scheme, f'AR(1), ρ = {rho:+.3f}')
    return out if dt == DT else f'{out}, step {dt:g}'


def factor_C(scheme, rho):
    return 0.0 if scheme == 'heun' else LN.ar1_factor(rho)


def summarise_C(runs, ito0):
    """Per scheme: the gauge coordinate's least-squares slope (per unit slow time) and end displacement y(S_C) − W₀, mean ± s.e. over the
    seeds, against the prediction (1 − ρ)/(1 + ρ)·v·(Itô drift at W₀), and their ratios to v·(Itô drift at W₀) on its largest component."""
    out = []
    W0 = np.asarray(W0_C, float)
    se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M))
    k = int(np.argmax(np.abs(ito0)))
    for scheme, rho, dt in SCHEMES_C:
        rs = [r for r in runs if r['scheme'] == scheme and abs(r['rho'] - rho) < 1e-12 and r['dt'] == dt]
        R = np.array([r['rec'] for r in rs])
        slope = np.array([np.polyfit(r['y'][:, 0], r['y'][:, 1:], 1)[0] for r in rs])
        endd = np.array([r['y'][-1, 1:] for r in rs]) - W0[None, :]
        f = factor_C(scheme, rho)
        out.append(dict(scheme=scheme, rho=rho, dt=dt, factor=f, slope=slope.mean(axis=0), slope_se=se(slope), end=endd.mean(axis=0),
                        end_se=se(endd), pred=f * V_C * ito0, ratio=slope.mean(axis=0)[k] / (V_C * ito0[k]),
                        ratio_se=se(slope)[k] / abs(V_C * ito0[k]), k=k, almax=float(np.abs(R[:, :, 1:4]).max()), fit=float(R[:, :, 7].mean()),
                        res=max(r['res'] for r in rs), n=len(rs)))
    return out


def report(A=None, B=None, C=None, D=None):
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
    if B is not None:
        L_ += ["## B. The selection's target", "",
               f"P_k = G_α⟨Φ_k²⟩_ρ = {v3(B['P'])}; W_min = {v3(B['W_min'])}.", "",
               f"W_ℓ = argmin Ψ_ℓ = {v3(B['W_l'])} (|a_k|/a_max there {v3(B['amp'], 3)}; ū there {v3(B['u_l'], 5)}, equal); "
               f"Proposition D‴.38's W♯ = {v3(B['W_sharp'])}.", "",
               f"Ψ_ℓ/v: {B['psi_l']:.5f} at W_ℓ, {B['psi_u']:.5f} at 1/3, {B['psi_0']:.5f} at W₀ = {v3(W0_C, 2)}; the closed form's replicator "
               f"flow linearised at W_ℓ has the eigenvalues {', '.join(f'{v.real:+.5f}' for v in B['eig'])} (per γv, per unit slow time).", ""]
    if C is not None:
        Ci, ito0 = C['sum'], C['ito0']
        L_ += [f"## C. The full loop with measurement noise only (τ = 0, γ = {GAMMA_C:g}, v = {V_C:g}, from W₀ = {v3(W0_C, 2)}, "
               f"{Ci[0]['n']} seeds, {S_C:g} slow units, steps of {DT:g} unless stated)", "",
               f"The gauge coordinate y(s) = Φ(z(s)): the noiseless flow ({L_MAP:g} slow units) from the state every {REC_PHI:g} slow units. Part A's "
               f"Itô drift at W₀ (γ = {GAMMA_C:g}, per unit v): {', '.join(f'{v:+.4e}' for v in ito0)}; the ratio is taken on W_{Ci[0]['k'] + 1}.", "",
               "| noise | (1 − ρ)/(1 + ρ) | dy/ds, the least-squares slope (± s.e.) | predicted (1 − ρ)/(1 + ρ)·v·(Itô drift) | the ratio to v·(Itô drift) (± s.e.) | "
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
               "the per-step variance.", "",
               "| noise | (1 − ρ)/(1 + ρ) | E[Q] at t = 0, 20, 40, …, 100 | the decay rate −log(E[Q(T)]/Q(0))/(σ²T) | the largest |Q − Q(0)| of a path | "
               "the mean |uv − 1| at T |", "|---|---|---|---|---|---|"]
        for r in D:
            rec = r['rec']
            L_.append(f"| {name_C(r['scheme'], r['rho'], DT)} | {factor_C(r['scheme'], r['rho']):.3f} | "
                      f"{', '.join(f'{x:.4f}' for x in rec[::2, 1])} | {r['rate']:+.4f} | {rec[-1, 3]:.2e} | {rec[-1, 4]:.4f} |")
        L_ += ["", "For noise drawn afresh at each step the Itô SDE's E[Q] decays exactly as e^{−σ²t} (dQ = −σ²Q dt + 0·dW); Heun's scheme keeps "
               "every path's Q to its discretisation error.", ""]
    L_.append("Numbers are produced by `experiments/exp_gauge_label.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_label_raw.npz')
    md = os.path.join(OUT, 'gauge_label.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCD']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCD')
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
        if 'A' in parts or 'C' in parts:
            h['A'] = pool.map_async(a_job, [(g, w) for g in GAMMAS_A for w in W_A], chunksize=1)
        if 'A' in h:
            res['A'] = h['A'].get(); E.save(raw, **res); log('A')
        if 'B' in parts:
            res['B'] = part_B(); E.save(raw, **res); log('B')
        if 'D' in h:
            res['D'] = h['D'].get(); E.save(raw, **res); log('D')
        if 'C' in h:
            runs = h['C'].get()
            ito0 = next(r['ito'] for r in res['A'] if r['gamma'] == GAMMA_C and np.allclose(np.r_[r['w'], 1 - r['w'].sum()], W0_C))
            res['C_runs'] = runs
            res['C'] = dict(ito0=ito0, sum=summarise_C(runs, ito0))
            E.save(raw, **res); log('C')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCD']))


if __name__ == '__main__':
    main()
