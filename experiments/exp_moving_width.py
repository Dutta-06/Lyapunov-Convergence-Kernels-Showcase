"""
Proposition D‴.59: the transit with moving widths — one kernel per Gaussian feature, whatever its width (`src/core/moving_width.py`).

The paradigm's law moves each kernel's width beside its centre (σ̇ = Gβa⟨r∂_σφ⟩, ċ = Gβa⟨r∂_cφ⟩); on one-dimensional samples a kernel is
θ = (c, σ), and `rich_transit_nd`'s sliding path, loop and polish run on `moving_width.BumpsCW`. The kernels start on the grid −1, −0.9,
…, 1 with the width σ₀ = 0.12; the samples are 321 points of [−1.6, 1.6]. Against the same path with the widths fixed at σ₀
(`rich_transit_nd.BumpsND`), to MAX_ENT_FIXED entries and within MAX_RHS_FIXED right-hand-side evaluations: on a feature narrower than
the kernels the fixed widths' blocks are stiff and their relaxation runs for hours (a path stopped by the budget is reported so).

A. One feature (A = 0.8 at c* = 0.03) at the width ratios σ*/σ₀ = 1/4 … 4 (D‴.58's 1.031, 1.118, 1.414, 2.236 among them): moving widths
   against fixed widths — the kernels formed, the final (c, σ, w) against (c*, σ*, A), the residual.
B. The one-kernel reduced energy in the continuum: its local minima on a fine (c, σ) grid (one, at the exact fit), the closed-form
   gradient against finite differences.
C. Random three-feature targets (widths log-uniform in [0.05, 0.3], amplitudes ±U(0.4, 1), centres in [−0.7, 0.7] at least Δ apart):
   moving widths' exact recovery against fixed widths' kernel count and residual.
D. A narrow feature (σ = 0.06, A = 0.5) on a wide pedestal (σ = 0.3, A = 0.8) at the distances d: the kernels formed with moving widths.
E. Under noise (C's targets at Δ = 0.6, s = 10⁻³, 10⁻²): the discrepancy rule's stop (the residual's rms within (1 + δ)s), moving widths
   against fixed: the kernels and the error against the truth, against the oracle's s(3K/N)^{1/2}.
F. The real loop at small masses (ε = 10⁻⁶, the free law on its curves, `rich_transit_nd.loop_nd` with moving widths) on C's first
   targets: its formed kernels against the path's.
"""
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.core import rich_transit as RT
from src.core import rich_transit_nd as RN
from src.core import moving_width as MW
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('MW_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('MW_WORKERS', '2'))
XS = np.linspace(-1.6, 1.6, 321)
S0, C0 = 0.12, np.linspace(-1.0, 1.0, 21)
A_A, C_A = 0.8, 0.03
RATIOS_A = (0.25, 0.5, 0.75, float(np.sqrt(1 + 0.25 ** 2)), float(np.sqrt(1 + 0.5 ** 2)), float(np.sqrt(2.0)), float(np.sqrt(5.0)), 3.0, 4.0)
MAX_ENT_FIXED = 4                                                       # fixed widths: the path to 4 entries (a narrow feature's slow blocks)
MAX_RHS_FIXED, MAX_RHS_MOVING = 300_000, 1_000_000                      # the paths' budgets of right-hand-side evaluations (a moving-width
                                                                        # path to its exact fit spends ~2·10³)
SEED_C, N_C, DELTAS_C, K_C, WID_C = 590, 40, (0.4, 0.6), 3, (0.05, 0.3)
EXACT_TOL = 1e-6                                                        # exact: every feature's (c, σ, w) within it, the residual < 1e-8
D_D = (0.0, 0.1, 0.2, 0.3, 0.45, 0.6)
LEVELS_E, DELTA_E, N_E = (1e-3, 1e-2), 0.1, 20
EPS_F, N_F, TAU_F = 1e-6, 4, 1.3
EPS_G, TAU_G, S_G = (1e-6, 1e-8), 400.0, (1e3, 1e5, 1e7, 1e9, 1e11)  # G: the stalled blocks' flow to s = 1e11 and the loop to τ = 400
PAIR_TOL = 0.05                                                         # a coinciding pair: |Δc|/σ + |Δlog σ| below it (a same-sign
                                                                        # pair a little apart fits one Gaussian to O(Δc⁴))


def law():
    return X.law()


def theta0():
    return np.stack([C0, np.full(C0.size, S0)], -1)


def feature_sum(cs, ss, As, x=XS):
    return np.exp(-(x[:, None] - np.asarray(cs)[None, :]) ** 2 / (2 * np.asarray(ss)[None, :] ** 2)) @ np.asarray(As)


def run_paths(f, max_entries=None, max_events=40, fixed=True):
    """The sliding path with moving widths (to max_entries) and with fixed widths (to MAX_ENT_FIXED entries) on the samples' target f."""
    L = law()
    Rm = RN.sliding_iss_nd(L, MW.BumpsCW(XS, S0), f, theta0(), np.ones(C0.size), max_events=max_events, max_entries=max_entries,
                           exit_kkt=True, max_rhs=MAX_RHS_MOVING)
    Rf = (RN.sliding_iss_nd(L, RN.BumpsND(XS[:, None], S0), f, C0[:, None], np.ones(C0.size), max_events=max_events,
                            max_entries=MAX_ENT_FIXED, exit_kkt=True, max_rhs=MAX_RHS_FIXED) if fixed else None)
    return Rm, Rf


def match(R, cs, ss, As):
    """The active kernels against the features: the count, and for each feature the nearest active kernel — in (c, log σ) with moving
    widths, so that features at one centre are told apart by their widths — with its |Δc|, |Δσ| and |Δw|; exact if as many kernels as
    features, each within EXACT_TOL, the residual below 1e-8."""
    act = np.flatnonzero(R['act'])
    th, w = np.atleast_2d(R['c'][act]), R['w'][act]
    errs = []
    for c_, s_, a_ in zip(cs, ss, As):
        if not act.size:
            errs.append((np.inf, np.inf, np.inf))
            continue
        dist = np.abs(th[:, 0] - c_) / S0 + (np.abs(np.log(np.abs(th[:, 1]) / s_)) if th.shape[1] > 1 else 0.0)
        j = int(np.argmin(dist))
        ds = abs(th[j, 1] - s_) if th.shape[1] > 1 else abs(S0 - s_)
        errs.append((abs(th[j, 0] - c_), ds, abs(w[j] - a_)))
    errs = np.array(errs)
    rms = R['hist'][-1][3]
    return dict(n=int(act.size), err=errs.max(0).tolist(), rms=float(rms),
                exact=bool(act.size == len(cs) and errs.max() < EXACT_TOL and rms < 1e-8),
                entries=sum(e[2] == 'enter' for e in R['events']), exits=sum(e[2] == 'exit' for e in R['events']),
                theta=th.tolist(), w=w.tolist(), oob=bool(R.get('out_of_budget', False)), nfev=int(R.get('nfev', 0)))


def pairs_of(theta, w):
    """The coinciding pairs among active kernels: |Δc|/σ + |Δlog σ| < PAIR_TOL (σ their mean width), with their coefficients."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    out = []
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
            if abs(th[i, 0] - th[j, 0]) / sm + abs(np.log(abs(th[i, 1]) / abs(th[j, 1]))) < PAIR_TOL:
                out.append((float(w[i]), float(w[j])))
    return out


def outcome(m, As):
    """A moving-width path's outcome (referee-free classification of C's and D's rows): 'exact' (as many formed kernels — |w| above 1% of
    the smallest amplitude — as features, within EXACT_TOL, the residual below 1e-8; idle entrants below 1% allowed), 'coinciding pair'
    (the residual below 1e-8 with two kernels at one (c, σ): a split or a cancelling dipole), 'stopped by the budget' (a block that never
    settles: G), or 'not fitted'."""
    w = np.asarray(m['w'], float)
    formed = np.abs(w) > 1e-2 * np.abs(np.asarray(As)).min()
    if m.get('oob'):
        return 'stopped by the budget'
    if m['rms'] < 1e-8 and formed.sum() == len(As) and max(m['err']) < EXACT_TOL:
        return 'exact'
    if m['rms'] < 1e-8 and pairs_of(np.atleast_2d(m['theta'])[formed], w[formed]):
        return 'coinciding pair'
    return 'not fitted'


def a_job(ratio):
    t0 = time.time()
    ss = ratio * S0
    f = feature_sum([C_A], [ss], [A_A])
    Rm, Rf = run_paths(f)
    return dict(ratio=ratio, moving=match(Rm, [C_A], [ss], [A_A]), fixed=match(Rf, [C_A], [ss], [A_A]),
                tau=float(Rm['events'][0][0]) if Rm['events'] else np.nan, secs=time.time() - t0)


def part_A():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, RATIOS_A, chunksize=1)


def part_B():
    """The reduced energy on a (c, σ) grid around the exact fit for each ratio: its grid-local minima, and the gradient's check."""
    out = []
    cs_, A = 0.0, 1.0
    for ratio in RATIOS_A:
        ss = ratio * S0
        cg = np.linspace(-6 * ss, 6 * ss, 241)
        sg = np.geomspace(ss / 8, 8 * ss, 241)
        J = MW.reduced_energy(cg[:, None], sg[None, :], A, cs_, ss)
        inner = J[1:-1, 1:-1]
        mins = np.argwhere((inner < J[:-2, 1:-1]) & (inner < J[2:, 1:-1]) & (inner < J[1:-1, :-2]) & (inner < J[1:-1, 2:]))
        locs = [(float(cg[i + 1]), float(sg[j + 1])) for i, j in mins]
        rng = np.random.default_rng(SEED_C + int(100 * ratio))
        gerr = 0.0
        q = lambda c, s: MW.overlap(c, s, cs_, ss) ** 2 / (np.sqrt(np.pi) * s)        # J = ½A²(π^{1/2}σ* − q): differenced without
        for _ in range(20):                                                            # the constant, which would swamp a tail's gradient
            c, s = rng.normal(0, 2 * ss), ss * np.exp(rng.normal(0, 0.7))
            g = MW.one_kernel_grad(c, s, A, cs_, ss)
            h = 1e-6 * ss
            fd = -0.5 * A ** 2 * np.array([(q(c + h, s) - q(c - h, s)) / (2 * h), (q(c, s + h) - q(c, s - h)) / (2 * h)])
            gerr = max(gerr, float(np.abs(fd - g).max() / max(np.abs(g).max(), 1e-300)))
        out.append(dict(ratio=ratio, minima=locs, J_min=float(J.min()), grad_err=gerr))
    return out


def targets_C(rng, n, dmin):
    out = []
    while len(out) < n:
        cs = rng.uniform(-0.7, 0.7, K_C)
        if min(abs(cs[i] - cs[j]) for i in range(K_C) for j in range(i)) < dmin:
            continue
        ss = np.exp(rng.uniform(np.log(WID_C[0]), np.log(WID_C[1]), K_C))
        As = rng.choice([-1.0, 1.0], K_C) * rng.uniform(0.4, 1.0, K_C)
        out.append((cs, ss, As))
    return out


def c_job(args):
    q, dmin, (cs, ss, As) = args
    t0 = time.time()
    f = feature_sum(cs, ss, As)
    Rm, Rf = run_paths(f)
    return dict(q=q, dmin=dmin, cs=cs.tolist(), ss=ss.tolist(), As=As.tolist(), moving=match(Rm, cs, ss, As), fixed=match(Rf, cs, ss, As),
                secs=time.time() - t0)


def part_C():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    jobs = [(q, d, t) for d in DELTAS_C for q, t in enumerate(targets_C(rng, N_C, d))]
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, jobs, chunksize=2)


def d_job(d):
    cs, ss, As = [0.0, d], [0.3, 0.06], [0.8, 0.5]
    f = feature_sum(cs, ss, As)
    Rm, Rf = run_paths(f)
    return dict(d=d, moving=match(Rm, cs, ss, As), fixed=match(Rf, cs, ss, As))


def part_D():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(d_job, D_D, chunksize=1)


def blocks(R):
    """The relaxed blocks along a path: (kernels, the residual's rms) after each event, with the block."""
    out = []
    for tau, c, w, rms in R['hist'][1:]:
        act = np.flatnonzero(w != 0)
        out.append(dict(k=int(act.size), rms=float(rms), c=np.atleast_2d(c)[act], w=w[act]))
    return out


def e_job(args):
    q, s, (cs, ss, As) = args
    f = feature_sum(cs, ss, As)
    xi = s * np.random.default_rng(SEED_C + 1000 * q + int(1e4 * s)).standard_normal(XS.size)
    Rm, _ = run_paths(f + xi, max_entries=K_C + 2, fixed=False)
    out = dict(q=q, s=s)
    for name, R, moving in (('moving', Rm, True),):
        bl = [b for b in blocks(R) if b['k'] > 0]
        stop = next((b for b in bl if b['rms'] <= (1 + DELTA_E) * s), None)
        if stop is None:
            out[name] = dict(stopped=False)
            continue
        th = stop['c']
        P = (np.exp(-(XS[:, None] - th[:, 0][None, :]) ** 2 / (2 * th[:, 1][None, :] ** 2)) if moving else
             np.exp(-(XS[:, None] - th[:, 0][None, :]) ** 2 / (2 * S0 ** 2)))
        out[name] = dict(stopped=True, k=stop['k'], err=float(np.sqrt(np.mean((P @ stop['w'] - f) ** 2))))
    return out


def part_E():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    targets_C(rng, N_C, DELTAS_C[0])                                   # C's stream: the Δ = 0.4 targets first, then Δ = 0.6's
    tg = targets_C(rng, N_C, DELTAS_C[1])[:N_E]
    with Pool(N_WORKERS) as pool:
        return pool.map(e_job, [(q, s, t) for s in LEVELS_E for q, t in enumerate(tg)], chunksize=2)


def f_job(args):
    q, (cs, ss, As) = args
    t0 = time.time()
    L = law()
    K = MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    R = RN.sliding_iss_nd(L, K, f, theta0(), np.ones(C0.size), max_rhs=MAX_RHS_MOVING)
    taus = np.linspace(0, TAU_F * R['events'][-1][0], 41)
    _, W, C, Lam, al = RN.loop_nd(L, K, f, theta0(), EPS_F, taus, rtol=1e-9, atol=1e-12)
    big = np.abs(W[-1]) > 1e-2 * np.abs(As).min()
    path = sorted(int(j) for j in np.flatnonzero(R['act']))
    return dict(q=q, path=path, formed=sorted(int(j) for j in np.flatnonzero(big)),
                dtheta=float(np.abs(C[-1][path] - R['c'][path]).max()), dw=float(np.abs(W[-1][path] - R['w'][path]).max()),
                rest=float(np.abs(np.delete(W[-1], path)).max()), secs=time.time() - t0)


def part_F():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    targets_C(rng, N_C, DELTAS_C[0])
    tg = targets_C(rng, N_C, DELTAS_C[1])[:N_F]
    with Pool(N_WORKERS) as pool:
        return pool.map(f_job, list(enumerate(tg)), chunksize=1)


def g_job(args):
    """G: a path stopped by its budget (C) — its stalled block from the entry that started it (the path to that many entries, the
    entrant at its seed): the block's flow (the fast flow on the critical curve, LSODA with the exact Jacobian) to s = 1e11, the largest
    |w| and the closest opposite-sign pair (|Δc|/σ + |Δlog σ|) at S_G, the residual; and the real loop at the masses EPS_G to τ = TAU_G:
    its formed kernels, coinciding pairs with their coefficients, the largest |w| over the largest amplitude, the residual."""
    dmin, q, (cs, ss, As), n_ent = args
    from scipy.integrate import solve_ivp
    t0 = time.time()
    L = law()
    K = MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    R = RN.sliding_iss_nd(L, K, f, theta0(), np.ones(C0.size), max_entries=n_ent, exit_kkt=True)
    idx = np.flatnonzero(R['act'])
    al = np.sign(R['w'][idx]) * RT._alpha_of(L, np.abs(R['w'][idx]), RT.coef_crit)
    z0 = np.r_[al, R['c'][idx].ravel()]
    k = idx.size
    sol = solve_ivp(RN.critical_rhs_nd(L, K, f), (0.0, S_G[-1]), z0, method='LSODA', jac=RN.critical_jac_nd(L, K, f), rtol=1e-10,
                    atol=1e-14, t_eval=S_G)
    flow = []
    for i in range(sol.y.shape[1]):
        w = RT.critical_w(L, sol.y[:k, i])
        th = sol.y[k:, i].reshape(k, 2)
        sep = min((abs(th[a, 0] - th[b, 0]) / (0.5 * (abs(th[a, 1]) + abs(th[b, 1]))) + abs(np.log(abs(th[a, 1]) / abs(th[b, 1])))
                   for a in range(k) for b in range(a + 1, k) if w[a] * w[b] < 0), default=None)
        flow.append(dict(s=float(sol.t[i]), wmax=float(np.abs(w).max()), sep=None if sep is None else float(sep),
                         rms=float(np.sqrt(np.mean((f - K.synth(th, w)) ** 2)))))
    loops = []
    for eps in EPS_G:
        _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), eps, np.array([0.0, TAU_G]), rtol=1e-9, atol=1e-12)
        big = np.abs(W[-1]) > 1e-2 * np.abs(As).min()
        loops.append(dict(eps=eps, formed=int(big.sum()), pairs=pairs_of(C[-1][big], W[-1][big]),
                          wmax=float(np.abs(W[-1]).max() / np.abs(As).max()), rms=float(np.sqrt(np.mean((f - K.synth(C[-1], W[-1])) ** 2)))))
    return dict(dmin=dmin, q=q, entries=n_ent, flow=flow, status=int(sol.status), loops=loops, secs=time.time() - t0)


def dipole_exponents(flow):
    """The last two outputs' log-slopes of the largest |w| and of the opposite-sign pair's distance against s (the dipole: +1/3, −1/3 —
    on the free law's exponential branch the coefficient's mobility grows like |w| and the pair's energy excess falls like 1/w² at a
    fixed moment)."""
    a, b = flow[-2], flow[-1]
    if a['sep'] is None or b['sep'] is None:
        return None
    ls = np.log(b['s'] / a['s'])
    return float(np.log(b['wmax'] / a['wmax']) / ls), float(np.log(b['sep'] / a['sep']) / ls)


def part_G(C):
    """G: C's targets whose moving-width path was stopped by its budget."""
    from multiprocessing import Pool
    jobs = [(r['dmin'], r['q'], (np.array(r['cs']), np.array(r['ss']), np.array(r['As'])), int(r['moving']['entries'])) for r in C
            if r['moving'].get('oob')]
    with Pool(N_WORKERS) as pool:
        return pool.map(g_job, jobs, chunksize=1)


def report(A=None, B=None, C=None, D=None, Ee=None, F=None, G=None):
    N = XS.size
    L_ = ['# Proposition D‴.59: the transit with moving widths (`exp_moving_width.py`)', '',
          f"The sliding path with moving widths (θ = (c, σ), `moving_width.BumpsCW`) and with the widths fixed at σ₀ = {S0:g} "
          f"(`rich_transit_nd.BumpsND`); {C0.size} kernels started at −1, −0.9, …, 1; {N} samples of [{XS[0]:g}, {XS[-1]:g}]; exact: as many "
          f"kernels as features, each feature's centre, width and amplitude within {EXACT_TOL:g}, the residual's rms below 1e-08.", '']
    if A is not None:
        L_ += ['## A. One feature', '',
               f"A = {A_A:g} at c* = {C_A:g}, the width σ* = ratio·σ₀; moving widths: the kernels, the largest error of (c, σ, w), the "
               f"residual; fixed widths (up to {MAX_ENT_FIXED} entries, {MAX_RHS_FIXED:.0e} evaluations — † the budget spent first): the kernels and the residual:", '',
               '| σ*/σ₀ | moving: kernels | exact | max error (c, σ, w) | residual | fixed: kernels | residual |', '|---|---|---|---|---|---|---|']
        for r in A:
            m, f_ = r['moving'], r['fixed']
            L_.append(f"| {r['ratio']:.4g} | {m['n']} | {'yes' if m['exact'] else 'no'} | {', '.join(f'{v:.1e}' for v in m['err'])} | "
                      f"{m['rms']:.1e} | {f_['n']}{'†' if f_.get('oob') else ''} | {f_['rms']:.2e} |")
        L_.append('')
    if B is not None:
        L_ += ['## B. The one-kernel reduced energy (the continuum)', '',
               "J(c, σ) = ½‖f*‖² − ½⟨f*, φ_{c,σ}⟩²/‖φ_{c,σ}‖² on a 241 × 241 grid (c within ±6σ*, σ from σ*/8 to 8σ*, A = 1, c* = 0): its "
               "grid-local minima and the closed-form gradient's largest relative error against central differences (20 random points):", '',
               '| σ*/σ₀ | local minima (c, σ/σ*) | min J | gradient error |', '|---|---|---|---|']
        for r in B:
            ss = r['ratio'] * S0
            L_.append(f"| {r['ratio']:.4g} | " + ', '.join(f"({c:+.1e}, {s / ss:.3f})" for c, s in r['minima']) + f" | {r['J_min']:.1e} | "
                      f"{r['grad_err']:.1e} |")
        L_.append('')
    if C is not None:
        L_ += ['## C. Random three-feature targets', '',
               f"{N_C} targets per minimal separation (seed {SEED_C}; widths log-uniform in [{WID_C[0]:g}, {WID_C[1]:g}], amplitudes ±U(0.4, 1)):", '',
               '| Δ ≥ | moving widths: exact | kernels (median [range]) | residual (median [max]) | fixed widths: kernels (median [range]) | '
               'residual (median [max]) |', '|---|---|---|---|---|---|']
        for d in DELTAS_C:
            rs = [r for r in C if r['dmin'] == d]
            km, kf = [r['moving']['n'] for r in rs], [r['fixed']['n'] for r in rs]
            em, ef = [r['moving']['rms'] for r in rs], [r['fixed']['rms'] for r in rs]
            L_.append(f"| {d:g} | {sum(r['moving']['exact'] for r in rs)} of {len(rs)} | {np.median(km):g} [{min(km)}, {max(km)}] | "
                      f"{np.median(em):.1e} [{max(em):.1e}] | {np.median(kf):g} [{min(kf)}, {max(kf)}] ({sum(r['fixed'].get('oob', False) for r in rs)}† "
                      f"stopped by the budget) | {np.median(ef):.1e} [{max(ef):.1e}] |")
        cls = {}
        for r in C:
            cls.setdefault((r['dmin'], outcome(r['moving'], r['As'])), []).append(r)
        L_ += ['', "The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude): " + '; '.join(
            f"Δ ≥ {d:g}: " + ', '.join(f"{len(cls.get((d, o), []))} {o}" for o in ('exact', 'coinciding pair', 'stopped by the budget',
                                                                                    'not fitted')) for d in DELTAS_C) + '.']
        bad = [r for r in C if outcome(r['moving'], r['As']) != 'exact']
        L_ += ['', "Not exact with moving widths: " + ('; '.join(
            f"Δ ≥ {r['dmin']:g}, target {r['q']} (widths {', '.join(f'{v:.3f}' for v in r['ss'])}; centres "
            f"{', '.join(f'{v:+.2f}' for v in r['cs'])}; amplitudes {', '.join(f'{v:+.2f}' for v in r['As'])}): "
            f"{outcome(r['moving'], r['As'])}, {r['moving']['n']} kernels"
            + (f" (the pair's coefficients {', '.join(f'{a:+.3f}/{b:+.3f}' for a, b in pairs_of(r['moving']['theta'], r['moving']['w']))})"
               if pairs_of(r['moving']['theta'], r['moving']['w']) else '')
            + f", residual {r['moving']['rms']:.1e}" for r in bad) or 'none') + '.', '']
    if D is not None:
        L_ += ['## D. A narrow feature on a wide pedestal', '',
               "The pedestal 0.8φ(·; 0, 0.3) and the narrow feature 0.5φ(·; d, 0.06):", '',
               '| d | moving: kernels | exact | max error (c, σ, w) | residual | fixed: kernels | residual |', '|---|---|---|---|---|---|---|']
        for r in D:
            m, f_ = r['moving'], r['fixed']
            L_.append(f"| {r['d']:g} | {m['n']} | {'yes' if m['exact'] else 'no'} | {', '.join(f'{v:.1e}' for v in m['err'])} | {m['rms']:.1e} | "
                      f"{f_['n']}{'†' if f_.get('oob') else ''} | {f_['rms']:.2e} |")
        L_.append('')
    if Ee is not None:
        L_ += ['## E. Under noise: the discrepancy stop', '',
               f"C's first {N_E} targets at Δ ≥ {DELTAS_C[1]:g}, one draw each, the path with moving widths to {K_C + 2} entries (past "
               "the features a noise kernel can shrink onto a sample: the fit has no minimiser); the first block whose residual's rms (against "
               f"the data) is within (1 + {DELTA_E:g})s; its kernels and its error against the truth (rms on the samples), median [range], "
               "against the oracle's s(3K/N)^{1/2} (three parameters per kernel, K = 3):", '',
               '| s | stopped | kernels | error | the oracle |', '|---|---|---|---|---|']
        for s in LEVELS_E:
            rs = [r for r in Ee if r['s'] == s]
            row = f"| {s:g} |"
            for name in ('moving',):
                st = [r[name] for r in rs if r[name]['stopped']]
                ks, es = [x['k'] for x in st], [x['err'] for x in st]
                row += (f" {len(st)} of {len(rs)} | {np.median(ks):g} [{min(ks)}, {max(ks)}] | {np.median(es):.2e} [{min(es):.2e}, "
                        f"{max(es):.2e}] |" if st else f" 0 of {len(rs)} | — | — |")
                if name == 'moving':
                    row += f" {s * np.sqrt(3 * K_C / N):.2e} |"
            L_.append(row)
        L_.append('')
    if F is not None:
        L_ += ['## F. The real loop at a small mass', '',
               f"The free law with moving centres and widths from α = 0, masses ε = {EPS_F:g} (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to "
               f"{TAU_F:g} times the path's last entry, on C's first {N_F} targets at Δ ≥ {DELTAS_C[1]:g}: its formed kernels (|w| above 1% of "
               "the smallest amplitude) against the path's, and the path's kernels' largest deviation:", '',
               "| target | the path's kernels | formed by the loop | max \\|Δθ\\| | max \\|Δw\\| | the others' largest \\|w\\| |", '|---|---|---|---|---|---|']
        for r in F:
            L_.append(f"| {r['q']} | {r['path']} | {r['formed']} | {r['dtheta']:.1e} | {r['dw']:.1e} | {r['rest']:.1e} |")
        L_.append('')
    if G is not None:
        L_ += ['## G. The blocks that never settle: a dipole', '',
               f"C's paths stopped by their budget: the stalled block from its entry (the path to that many entries), its flow on the critical "
               f"curve to s = {S_G[-1]:g} — the largest |w| and the closest opposite-sign pair (|Δc|/σ + |Δlog σ|) at s = "
               f"{', '.join(f'{v:g}' for v in S_G)} —; and the real loop at ε = {', '.join(f'{e:g}' for e in EPS_G)} to τ = {TAU_G:g}: its "
               "formed kernels, the coinciding pairs' coefficients, the largest |w| over the largest amplitude and the residual:", '',
               '| Δ ≥ | target | entries | the flow: largest \\|w\\| | the opposite-sign pair\'s distance | the last factor 100 in s: the '
               'exponents of \\|w\\| and the distance | residual | the loop: formed, pairs, max\\|w\\|/max\\|A\\|, residual |',
               '|---|---|---|---|---|---|---|---|']
        for r in G:
            fl = r['flow']
            lp = '; '.join(f"ε = {x['eps']:g}: {x['formed']}, " + (', '.join(f"{a:+.3f}/{b:+.3f}" for a, b in x['pairs']) or 'none')
                           + f", {x['wmax']:.2f}, {x['rms']:.1e}" for x in r['loops'])
            wm = ', '.join(f"{x['wmax']:.3g}" for x in fl)
            ex = dipole_exponents(fl)
            L_.append(f"| {r['dmin']:g} | {r['q']} | {r['entries']} | {wm} | "
                      + ', '.join('—' if x['sep'] is None else f"{x['sep']:.2e}" for x in fl)
                      + (f" | {ex[0]:+.3f}, {ex[1]:+.3f}" if ex is not None else ' | —') + f" | {fl[-1]['rms']:.2e} | {lp} |")
        L_.append('')
    L_ += ['Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'moving_width_raw.npz')
    md = os.path.join(OUT, 'moving_width.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'E', 'F', 'G')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D), ('E', part_E), ('F', part_F),
                  ('G', lambda: part_G(res['C']))):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
