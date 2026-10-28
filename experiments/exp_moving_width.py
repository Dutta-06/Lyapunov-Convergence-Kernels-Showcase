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
   targets: its formed kernels against the path's, at TAUS_F times the path's last entry.
G. C's paths stopped by their budget: the stalled block's state when the budget ran out (the path as C's), its flow from the entry on
   the critical curve to s = 1e11 with the dipole law's constant (K = v²(E − E_inf), E_inf the optimal dipole fit's), and the real loop
   at EPS_G to the times TAUS_G (referee round 57).
H. Two features of one sign closer than the kernels' starting width (±0.08), symmetric and perturbed: the kernels of the other sign
   (referee round 57).
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
                                                                        # path to its exact fit spends 10³–10⁴)
SEED_C, N_C, DELTAS_C, K_C, WID_C = 590, 40, (0.4, 0.6), 3, (0.05, 0.3)
EXACT_TOL = 1e-6                                                        # exact: every feature's (c, σ, w) within it, the residual < 1e-8
D_D = (0.0, 0.1, 0.2, 0.3, 0.45, 0.6)
LEVELS_E, DELTA_E, N_E = (1e-3, 1e-2), 0.1, 20
EPS_F, N_F, TAUS_F = 1e-6, 4, (1.3, 6.0)
EPS_G, TAUS_G, S_G = (1e-6, 1e-8), (400.0, 1600.0, 6400.0), (1e3, 1e5, 1e7, 1e9, 3.2e10, 1e11)  # G: the flow to s = 1e11, the loop to τ
H_TARGETS = (('σ* 0.04, equal', (-0.08, 0.08), (0.04, 0.04), (1.0, 1.0)), ('σ* 0.06, equal', (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)),
             ('σ* 0.06, A₂ 0.99', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.99)), ('σ* 0.06, A₂ 0.95', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.95)),
             ('σ* 0.06, A₂ 0.9', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.9)), ('σ* 0.06, c₂ +0.01', (-0.08, 0.09), (0.06, 0.06), (1.0, 1.0)),
             ('σ* 0.06, σ₂ 0.065', (-0.08, 0.08), (0.06, 0.065), (1.0, 1.0)), ('σ* 0.04, A₂ 0.6', (-0.08, 0.08), (0.04, 0.04), (1.0, 0.6)),
             ('σ* 0.06, A₂ 0.6', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.6)))
PAIR_TOL = 0.01                                                         # a coinciding pair: |Δc|/σ + |Δlog σ| below it (a same-sign
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
    st = R.get('stalled')
    return dict(n=int(act.size), err=errs.max(0).tolist(), rms=float(rms),
                exact=bool(act.size == len(cs) and errs.max() < EXACT_TOL and rms < 1e-8),
                entries=sum(e[2] == 'enter' for e in R['events']), exits=sum(e[2] == 'exit' for e in R['events']),
                theta=th.tolist(), w=w.tolist(), oob=bool(R.get('out_of_budget', False)), nfev=int(R.get('nfev', 0)),
                stalled=None if st is None else dict(theta=np.atleast_2d(st['c']).tolist(), w=np.asarray(st['w']).tolist()))


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


def closest_pair(theta, w, opposite=False):
    """The closest pair of kernels (|Δc|/σ + |Δlog σ|, σ their mean width; of opposite signs if asked): (distance, i, j) or None."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    best = None
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            if opposite and w[i] * w[j] >= 0:
                continue
            sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
            d = abs(th[i, 0] - th[j, 0]) / sm + abs(np.log(abs(th[i, 1]) / abs(th[j, 1])))
            if best is None or d < best[0]:
                best = (float(d), i, j)
    return best


def width_split(theta, w, i, j):
    """A pair's width split against the third moment's (referee round 57): with Δc = c_j − c_i, W = w_i + w_j and σ their mean width the
    split σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ) cancels the pair's third moment, so the pair fits one Gaussian to O(Δc⁴). Returns (Δc,
    the observed split, the predicted)."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    if th[i, 0] > th[j, 0]:
        i, j = j, i
    dc = th[j, 0] - th[i, 0]
    sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
    return float(dc), float(abs(th[j, 1]) - abs(th[i, 1])), float((w[j] - w[i]) * dc ** 2 / (6 * (w[i] + w[j]) * sm))


def n_formed(m, As):
    """The formed kernels of a path's end: |w| above 1% of the smallest amplitude (as F)."""
    return int((np.abs(np.asarray(m['w'], float)) > 1e-2 * np.abs(np.asarray(As)).min()).sum())


def outcome(m, As):
    """A moving-width path's outcome: 'exact' (as many formed kernels — |w| above 1% of the smallest amplitude — as features, within
    EXACT_TOL, the residual below 1e-8; idle entrants below 1% allowed), 'coinciding pair' (the residual below 1e-8 with two kernels at one
    (c, σ), PAIR_TOL: a split or a cancelling pair), 'cancellation' (the residual below 1e-8 with a coefficient beyond 1.5 times the
    largest amplitude or a width below the samples' spacing: several kernels cancelling), 'dipole run-off' (stopped by the budget with
    an opposite-sign pair in the stalled block, |w| beyond 10 times the largest amplitude and 0.05 apart: G), 'stopped by the budget'
    (another stall), or 'not fitted' (referee round 57)."""
    w = np.asarray(m['w'], float)
    As = np.abs(np.asarray(As))
    formed = np.abs(w) > 1e-2 * As.min()
    if m.get('oob'):
        st = m.get('stalled')
        if st is not None:
            ws = np.asarray(st['w'], float)
            cp = closest_pair(st['theta'], ws, opposite=True)
            if cp is not None and np.abs(ws).max() > 10 * As.max() and cp[0] < 0.05:
                return 'dipole run-off'
        return 'stopped by the budget'
    th = np.atleast_2d(m['theta'])
    if m['rms'] < 1e-8 and formed.sum() == len(As) and max(m['err']) < EXACT_TOL:
        return 'exact'
    if m['rms'] < 1e-8 and pairs_of(th[formed], w[formed]):
        return 'coinciding pair'
    if m['rms'] < 1e-8 and formed.any() and (np.abs(w[formed]).max() > 1.5 * As.max() or np.abs(th[formed, 1]).min() < XS[1] - XS[0]):
        return 'cancellation'
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
    out = dict(q=q, s=s, oracle=oracle_error(f, xi, cs, ss, As))
    for name, R, moving in (('moving', Rm, True),):
        bl = [b for b in blocks(R) if b['k'] > 0]
        stop = next((b for b in bl if b['rms'] <= (1 + DELTA_E) * s), None)
        if stop is None:
            out[name] = dict(stopped=False, oob=bool(R.get('out_of_budget', False)), rms=bl[-1]['rms'] if bl else np.nan)
            continue
        th = stop['c']
        P = (np.exp(-(XS[:, None] - th[:, 0][None, :]) ** 2 / (2 * th[:, 1][None, :] ** 2)) if moving else
             np.exp(-(XS[:, None] - th[:, 0][None, :]) ** 2 / (2 * S0 ** 2)))
        out[name] = dict(stopped=True, k=stop['k'], err=float(np.sqrt(np.mean((P @ stop['w'] - f) ** 2))))
    return out


def oracle_error(f, xi, cs, ss, As):
    """The oracle fit: the K-kernel moving-width least-squares fit of f + ξ started at the truth (the right support); its error against f
    (rms on the samples) — per draw, against the stopped fit's (referee round 57)."""
    from scipy.optimize import least_squares
    k = len(As)
    K = MW.BumpsCW(XS, S0)
    y = f + xi
    res = lambda z: y - K.synth(z[k:].reshape(k, 2), z[:k])
    sol = least_squares(res, np.r_[As, np.c_[cs, ss].ravel()], method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    return float(np.sqrt(np.mean((K.synth(sol.x[k:].reshape(k, 2), sol.x[:k]) - f) ** 2)))


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
    t_last = R['events'][-1][0]
    taus = np.r_[0.0, np.array(TAUS_F) * t_last]
    _, W, C, Lam, al = RN.loop_nd(L, K, f, theta0(), EPS_F, taus, rtol=1e-9, atol=1e-12)
    path = sorted(int(j) for j in np.flatnonzero(R['act']))
    pp = closest_pair(R['c'][path], R['w'][path])
    snaps = []
    for t in range(1, taus.size):
        big = np.abs(W[t]) > 1e-2 * np.abs(As).min()
        fj = np.flatnonzero(big)
        cp = closest_pair(C[t][fj], W[t][fj])
        snaps.append(dict(times=float(TAUS_F[t - 1]), formed=sorted(int(j) for j in fj), dtheta=float(np.abs(C[t][path] - R['c'][path]).max()),
                          dw=float(np.abs(W[t][path] - R['w'][path]).max()), rest=float(np.abs(np.delete(W[t], path)).max()),
                          pair=None if cp is None else dict(d=cp[0], dc=float(abs(C[t][fj[cp[1]], 0] - C[t][fj[cp[2]], 0])),
                                                            ds=float(abs(C[t][fj[cp[1]], 1] - C[t][fj[cp[2]], 1]))),
                          rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2)))))
    return dict(q=q, path=path, path_pair=None if pp is None else pp[0], path_rms=float(R['hist'][-1][3]), snaps=snaps,
                formed=snaps[0]['formed'], dtheta=snaps[0]['dtheta'], dw=snaps[0]['dw'], rest=snaps[0]['rest'], secs=time.time() - t0)


def part_F():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    targets_C(rng, N_C, DELTAS_C[0])
    tg = targets_C(rng, N_C, DELTAS_C[1])[:N_F]
    with Pool(N_WORKERS) as pool:
        return pool.map(f_job, list(enumerate(tg)), chunksize=1)


def dipole_law(L, K, f, th, w):
    """The dipole's law at a block's state (θ, w) (referee round 57): its pair — the largest |w| and the closest kernel of the other sign —
    (a, b), v = |w_a − w_b|/2 (the
    half-difference, moved at the mobility μ₀|w|, μ₀ = 4Gc/A on the free law's exponential branch), the moment p = w_aθ_a + w_bθ_b −
    (w_a + w_b)θ̄ about the midpoint θ̄ (its σ-share |p_σ|/|p|), and the optimal dipole fit — the other kernels, Wφ(θ̄) and p·∇_θφ(θ̄)
    by least squares from the state —: E_inf = ½⟨r²⟩ there. Returns dict(v, E, E_inf, rms_inf, p, mu0) or None."""
    from scipy.optimize import least_squares
    th, w = np.atleast_2d(np.asarray(th, float)), np.asarray(w, float)
    a = int(np.argmax(np.abs(w)))
    opp = [j for j in range(w.size) if w[j] * w[a] < 0]
    if not opp:
        return None
    dist = lambda j: (abs(th[a, 0] - th[j, 0]) / (0.5 * (abs(th[a, 1]) + abs(th[j, 1]))) + abs(np.log(abs(th[a, 1]) / abs(th[j, 1]))))
    b = min(opp, key=dist)
    oth = [i for i in range(w.size) if i not in (a, b)]
    m = len(oth)
    tb = 0.5 * (th[a] + th[b])
    p0 = w[a] * (th[a] - tb) + w[b] * (th[b] - tb)

    def model(z):
        g = K.synth(z[m:3 * m].reshape(m, 2), z[:m]) if m else 0.0
        t_ = z[3 * m + 1:3 * m + 3]
        return g + z[3 * m] * K.phi(t_[None, :])[:, 0] + K.dphi(t_[None, :])[:, 0, :] @ z[3 * m + 3:]
    z0 = np.r_[w[oth], th[oth].ravel(), w[a] + w[b], tb, p0]
    sol = least_squares(lambda z: (f - model(z)) / np.sqrt(f.size), z0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    E_inf = float(np.mean((f - model(sol.x)) ** 2) / 2)
    return dict(v=float(abs(w[a] - w[b]) / 2), E=float(np.mean((f - K.synth(th, w)) ** 2) / 2), E_inf=E_inf,
                rms_inf=float(np.sqrt(2 * E_inf)), p=sol.x[3 * m + 3:].tolist(), mu0=float(4 * L.G * L.c / L.A),
                dc=float(abs(th[a, 0] - th[b, 0])), ds=float(abs(th[a, 1] - th[b, 1])))


def g_job(args):
    """G: a path stopped by its budget (C) — C's path again, the stalled block's state when the budget ran out (its largest |w|, the
    opposite-sign pair's distance, its residual) and the polishes of that block; the block from the entry that started it (the path to
    that many entries, the entrant at its seed): its flow (the fast flow on the critical curve, LSODA with the exact Jacobian) to
    s = 1e11, the largest |w|, the closest opposite-sign pair (|Δc|/σ + |Δlog σ|) and the dipole law (`dipole_law`) at S_G; and the
    real loop at the masses EPS_G at the times TAUS_G: its formed kernels, the closest pair's distance, coefficients and midpoint's offset
    from the nearest feature (in its width), the coinciding pairs, the largest |w| over the largest amplitude, the residual."""
    dmin, q, (cs, ss, As), n_ent = args
    from scipy.integrate import solve_ivp
    t0 = time.time()
    L = law()
    K = MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    Rc, _ = run_paths(f, fixed=False)
    st, stalled = Rc['stalled'], None
    if st is not None:
        cp = closest_pair(st['c'], st['w'], opposite=True)
        stalled = dict(wmax=float(np.abs(st['w']).max()), sep=None if cp is None else cp[0], n=int(np.size(st['w'])),
                       rms=float(np.sqrt(np.mean((f - K.synth(np.atleast_2d(st['c']), st['w'])) ** 2))),
                       polish=Rc['polish_log'][-1] if Rc['polish_log'] else [], nfev=int(Rc['nfev']))
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
                         rms=float(np.sqrt(np.mean((f - K.synth(th, w)) ** 2))), theta=th.tolist(), w=w.tolist()))
    law_ = [dipole_law(L, K, f, x['theta'], x['w']) for x in flow[-2:]] if sol.status == 0 else []
    loops = []
    for eps in EPS_G:
        _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, TAUS_G], rtol=1e-9, atol=1e-12)
        for t, tau in enumerate(TAUS_G, start=1):
            big = np.flatnonzero(np.abs(W[t]) > 1e-2 * np.abs(As).min())
            cp = closest_pair(C[t][big], W[t][big])
            pr = None
            if cp is not None:
                i, j = big[cp[1]], big[cp[2]]
                mid = 0.5 * (C[t][i] + C[t][j])
                k = int(np.argmin(np.abs(np.asarray(cs) - mid[0]) / np.asarray(ss)))
                pr = dict(d=cp[0], w=(float(W[t][i]), float(W[t][j])), off=float(abs(mid[0] - cs[k]) / ss[k]))
            loops.append(dict(eps=eps, tau=tau, formed=int(big.size), pair=pr, pairs=pairs_of(C[t][big], W[t][big]),
                              wmax=float(np.abs(W[t]).max() / np.abs(As).max()), rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2)))))
    return dict(dmin=dmin, q=q, entries=n_ent, stalled=stalled, flow=flow, law=law_, status=int(sol.status), loops=loops,
                secs=time.time() - t0)


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


def h_job(args):
    """H: two features of one sign closer than the kernels' starting width — the moving-width path's outcome, the most negative
    coefficient and the most kernels of the other sign in any block along the path, the end's formed kernels of the other sign, the
    smallest formed width, the largest |w| (referee round 57)."""
    name, cs, ss, As = args
    t0 = time.time()
    f = feature_sum(cs, ss, As)
    Rm, _ = run_paths(f, fixed=False)
    m = match(Rm, cs, ss, As)
    w, th = np.asarray(m['w'], float), np.atleast_2d(m['theta'])
    formed = np.abs(w) > 1e-2 * min(abs(a) for a in As)
    ws = [h[2][h[2] != 0] for h in Rm['hist'][1:]]
    return dict(name=name, cs=list(cs), ss=list(ss), As=list(As), moving=m, outcome=outcome(m, As),
                neg_path=int(max([int((x < -1e-2 * min(As)).sum()) for x in ws], default=0)),
                wmin_path=float(min([x.min() for x in ws if x.size], default=0.0)), neg_end=int((formed & (w < 0)).sum()),
                smin=float(np.abs(th[formed, 1]).min()) if formed.any() else float('nan'), wmax=float(np.abs(w).max()) if w.size else 0.0,
                secs=time.time() - t0)


def part_H():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, list(H_TARGETS), chunksize=1)


def report(A=None, B=None, C=None, D=None, Ee=None, F=None, G=None, H=None):
    N = XS.size
    L_ = ['# Proposition D‴.59: the transit with moving widths (`exp_moving_width.py`)', '',
          f"The sliding path with moving widths (θ = (c, σ), `moving_width.BumpsCW`) and with the widths fixed at σ₀ = {S0:g} "
          f"(`rich_transit_nd.BumpsND`); {C0.size} kernels started at −1, −0.9, …, 1; {N} samples of [{XS[0]:g}, {XS[-1]:g}]; exact: as many "
          f"formed kernels (|w| above 1% of the smallest amplitude) as features, each feature's centre, width and amplitude within "
          f"{EXACT_TOL:g}, the residual's rms below 1e-08.", '']
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
    gst = {(g['dmin'], g['q']): g.get('stalled') for g in (G if G is not None else [])}

    def oc(r):
        o = outcome(r['moving'], r['As'])
        st = gst.get((r['dmin'], r['q']))
        if o == 'stopped by the budget' and st is not None and st['sep'] is not None and st['wmax'] > 10 * np.abs(r['As']).max() \
                and st['sep'] < 0.05:
            return 'dipole run-off'
        return o

    def split_text(m, As):
        w = np.asarray(m['w'], float)
        fj = np.flatnonzero(np.abs(w) > 1e-2 * np.abs(np.asarray(As)).min())
        th = np.atleast_2d(m['theta'])[fj]
        out = []
        for i in range(fj.size):
            for j in range(i + 1, fj.size):
                if pairs_of(th[[i, j]], w[fj][[i, j]]):
                    dc, ds, pr = width_split(th, w[fj], i, j)
                    out.append(f"{w[fj][i]:+.3f}/{w[fj][j]:+.3f}, Δc {dc:.2e}, the width split {ds:+.4e} (predicted {pr:+.4e})")
        return '; '.join(out)
    if C is not None:
        L_ += ['## C. Random three-feature targets', '',
               f"{N_C} targets per minimal separation (seed {SEED_C}; widths log-uniform in [{WID_C[0]:g}, {WID_C[1]:g}], amplitudes ±U(0.4, 1)):", '',
               '| Δ ≥ | moving widths: exact | formed kernels (median [range]) | residual (median [max]) | fixed widths: kernels (median '
               '[range]) | residual (median [max]) |', '|---|---|---|---|---|---|']
        for d in DELTAS_C:
            rs = [r for r in C if r['dmin'] == d]
            km, kf = [n_formed(r['moving'], r['As']) for r in rs], [r['fixed']['n'] for r in rs]
            em, ef = [r['moving']['rms'] for r in rs], [r['fixed']['rms'] for r in rs]
            L_.append(f"| {d:g} | {sum(oc(r) == 'exact' for r in rs)} of {len(rs)} | {np.median(km):g} [{min(km)}, {max(km)}] | "
                      f"{np.median(em):.1e} [{max(em):.1e}] | {np.median(kf):g} [{min(kf)}, {max(kf)}] ({sum(r['fixed'].get('oob', False) for r in rs)}† "
                      f"stopped by the budget) | {np.median(ef):.1e} [{max(ef):.1e}] |")
        cls = {}
        for r in C:
            cls.setdefault((r['dmin'], oc(r)), []).append(r)
        L_ += ['', "The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude; a coinciding pair within "
               f"{PAIR_TOL:g} in |Δc|/σ + |Δlog σ|; a budget stop a dipole run-off if G's stalled block has an opposite-sign pair under 0.05 "
               "with |w| beyond 10 times the largest amplitude): " + '; '.join(
            f"Δ ≥ {d:g}: " + ', '.join(f"{len(cls.get((d, o), []))} {o}" for o in ('exact', 'coinciding pair', 'cancellation', 'dipole run-off', 'stopped by the budget', 'not fitted')) for d in DELTAS_C) + '.']
        bad = [r for r in C if oc(r) != 'exact']
        L_ += ['', "Not exact with moving widths (the state before a stalled block; G gives the block): " + ('; '.join(
            f"Δ ≥ {r['dmin']:g}, target {r['q']} (widths {', '.join(f'{v:.3f}' for v in r['ss'])}; centres "
            f"{', '.join(f'{v:+.2f}' for v in r['cs'])}; amplitudes {', '.join(f'{v:+.2f}' for v in r['As'])}): "
            f"{oc(r)}, {r['moving']['n']} kernels"
            + (f" (the pair: {split_text(r['moving'], r['As'])})" if split_text(r['moving'], r['As']) else '')
            + f", residual {r['moving']['rms']:.1e}" for r in bad) or 'none') + '. The width split (σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), '
            "the pair's third moment cancelled) against the observed.", '']
    if D is not None:
        L_ += ['## D. A narrow feature on a wide pedestal', '',
               "The pedestal 0.8φ(·; 0, 0.3) and the narrow feature 0.5φ(·; d, 0.06):", '',
               '| d | moving: formed kernels | outcome | max error (c, σ, w) | residual | fixed: kernels | residual |', '|---|---|---|---|---|---|---|']
        for r in D:
            m, f_ = r['moving'], r['fixed']
            As_ = [0.8, 0.5]
            sp = split_text(m, As_)
            oc_ = outcome(m, As_) + (f" ({sp})" if sp else '')
            L_.append(f"| {r['d']:g} | {n_formed(m, As_)} | {oc_} | {', '.join(f'{v:.1e}' for v in m['err'])} | {m['rms']:.1e} | "
                      f"{f_['n']}{'†' if f_.get('oob') else ''} | {f_['rms']:.2e} |")
        L_.append('')
    if Ee is not None:
        L_ += ['## E. Under noise: the discrepancy stop', '',
               f"C's first {N_E} targets at Δ ≥ {DELTAS_C[1]:g}, one draw each, the path with moving widths to {K_C + 2} entries (past "
               "the features a noise kernel can shrink onto a sample: the fit has no minimiser); the first block whose residual's rms (against "
               f"the data) is within (1 + {DELTA_E:g})s; its kernels and its error against the truth (rms on the samples), median [range] and "
               "the rms over the draws, against the oracle's s(3K/N)^{1/2} (three parameters per kernel, K = 3, the rms of the errors), and "
               "per draw against the oracle fit (the K-kernel least-squares fit of the data started at the truth):", '',
               '| s | stopped | kernels | error: median [range]; rms | the oracle | the stopped error over the oracle fit\'s: median [range] '
               '(within 10⁻⁴ of it) '
               '| not stopped: the last block\'s rms |', '|---|---|---|---|---|---|---|']
        for s in LEVELS_E:
            rs = [r for r in Ee if r['s'] == s]
            st = [(r['moving'], r.get('oracle')) for r in rs if r['moving']['stopped']]
            ks, es = [x['k'] for x, o in st], [x['err'] for x, o in st]
            qo = [x['err'] / o for x, o in st if o]
            ns = [r['moving'] for r in rs if not r['moving']['stopped']]
            L_.append(f"| {s:g} | {len(st)} of {len(rs)} | " + (f"{np.median(ks):g} [{min(ks)}, {max(ks)}] | {np.median(es):.2e} "
                      f"[{min(es):.2e}, {max(es):.2e}]; {np.sqrt(np.mean(np.square(es))):.3e}" if st else '— | —') + f" | "
                      f"{s * np.sqrt(3 * K_C / N):.3e} | " + (f"{np.median(qo):.4f} [{min(qo):.4f}, {max(qo):.4f}] "
                                                              f"({sum(abs(q - 1) < 1e-4 for q in qo)} of {len(qo)})" if qo else '—') + ' | '
                      + (', '.join(f"{x.get('rms', np.nan):.2e}" + (' (budget)' if x.get('oob') else '') for x in ns) or '—') + ' |')
        L_.append('')
    if F is not None:
        L_ += ['## F. The real loop at a small mass', '',
               f"The free law with moving centres and widths from α = 0, masses ε = {EPS_F:g} (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to "
               f"{', '.join(f'{t:g}' for t in TAUS_F)} times the path's last entry, on C's first {N_F} targets at Δ ≥ {DELTAS_C[1]:g}: its formed "
               "kernels (|w| above 1% of the smallest amplitude) against the path's, the path's kernels' largest deviation, the others' largest "
               "|w|, the closest pair among the formed kernels (|Δc|/σ + |Δlog σ|; its |Δc|, |Δσ|) and the residual; the path's closest pair "
               "and residual:", '',
               "| target | the path's kernels (closest pair; residual) | times | formed by the loop | max \\|Δθ\\| | max \\|Δw\\| | the others' "
               "largest \\|w\\| | the loop's closest pair: distance (\\|Δc\\|, \\|Δσ\\|) | residual |", '|---|---|---|---|---|---|---|---|---|']
        for r in F:
            for x in r.get('snaps', [dict(times=None, formed=r['formed'], dtheta=r['dtheta'], dw=r['dw'], rest=r['rest'], pair=None,
                                          rms=np.nan)]):
                pp = r.get('path_pair')
                pa = x.get('pair')
                L_.append(f"| {r['q']} | {r['path']} (" + ('—' if pp is None else f"{pp:.1e}") + f"; {r.get('path_rms', np.nan):.1e}) | "
                          + ('—' if x['times'] is None else f"{x['times']:g}") + f" | {x['formed']} | {x['dtheta']:.1e} | {x['dw']:.1e} | "
                          f"{x['rest']:.1e} | " + ('—' if pa is None else f"{pa['d']:.1e} ({pa['dc']:.1e}, {pa['ds']:.1e})")
                          + f" | {x['rms']:.1e} |")
        L_.append('')
    if G is not None:
        L_ += ['## G. The blocks that never settle: a dipole', '',
               "C's paths stopped by their budget. The stalled block (C's path again; its state when the budget ran out): its kernels, the "
               "largest |w|, the closest opposite-sign pair's distance (|Δc|/σ + |Δlog σ|), its residual against the optimal dipole fit's "
               "(below), its polishes (status/evaluations/accepted, the largest coefficient over the block's). The block from its entry, "
               f"flowed on the critical curve to s = {S_G[-1]:g}: the largest |w| and the pair's distance at s = {', '.join(f'{v:g}' for v in S_G)}; "
               "over the last factor in s the exponents of the largest |w|, the distance, the pair's |Δc| and |Δσ| and the residual's "
               "excess over the dipole infimum; the dipole law's constant K = v²(E − E_inf) at the last two outputs (v the pair's "
               "half-difference |w_a − w_b|/2, E_inf the optimal dipole fit's — the other kernels, Wφ(θ̄) and p·∇_θφ(θ̄) by least squares), "
               "the growth of v³ over the last interval against 3μ₀K (μ₀ = 4Gc/A: the coefficient's mobility μ₀|w| on the exponential "
               "branch), the moment's σ-share |p_σ|/|p|:", '',
               '| Δ ≥ | target | entries | the stalled block: kernels, max\\|w\\|, distance, residual (infimum), polishes | the flow: largest '
               '\\|w\\| | distance | exponents: \\|w\\|, distance, \\|Δc\\|, \\|Δσ\\|, the residual\'s excess | K | Δv³/Δs, 3μ₀K | \\|p_σ\\|/\\|p\\| |',
               '|---|---|---|---|---|---|---|---|---|---|']
        for r in G:
            fl, lw, st = r['flow'], r.get('law') or [], r.get('stalled')
            ls = np.log(fl[-1]['s'] / fl[-2]['s'])
            ex = dipole_exponents(fl)
            stc = '—'
            if st is not None:
                pol = ' '.join(f"{x['status']}/{x['nfev']}/{'a' if x['accepted'] else 'r'}/{x['wmax_ratio']:.3g}" for x in st['polish'])
                stc = (f"{st['n']}, {st['wmax']:.4g}, " + ('—' if st['sep'] is None else f"{st['sep']:.2e}") + f", {st['rms']:.6e}"
                       + (f" ({lw[-1]['rms_inf']:.6e})" if lw else '') + f", {pol or '—'}")
            exs, kk, vv, ps = '—', '—', '—', '—'
            if ex is not None and len(lw) == 2 and lw[0] is not None and lw[1] is not None:
                a_, b_ = lw
                exs = (f"{ex[0]:+.4f}, {ex[1]:+.4f}, {np.log(b_['dc'] / a_['dc']) / ls:+.4f}, {np.log(b_['ds'] / a_['ds']) / ls:+.4f}, "
                       f"{np.log((fl[-1]['rms'] - b_['rms_inf']) / (fl[-2]['rms'] - b_['rms_inf'])) / ls:+.4f}")
                K2 = [x['v'] ** 2 * (x['E'] - x['E_inf']) for x in lw]
                kk = ', '.join(f"{k_:.5e}" for k_ in K2)
                vv = f"{(b_['v'] ** 3 - a_['v'] ** 3) / (fl[-1]['s'] - fl[-2]['s']):.5e}, {3 * b_['mu0'] * K2[1]:.5e}"
                ps = f"{abs(b_['p'][1]) / np.hypot(*b_['p']):.3f}"
            L_.append(f"| {r['dmin']:g} | {r['q']} | {r['entries']} | {stc} | {', '.join(format(x['wmax'], '.3g') for x in fl)} | "
                      + ', '.join('—' if x['sep'] is None else f"{x['sep']:.2e}" for x in fl) + f" | {exs} | {kk} | {vv} | {ps} |")
        L_ += ['', f"The real loop at the masses ε from α = 0 at the times τ: its formed kernels, the closest pair among them (its distance, "
               f"coefficients, and its midpoint's distance from the nearest feature in that feature's width), the coinciding pairs (within "
               f"{PAIR_TOL:g}), the largest |w| over the largest amplitude, the residual:", '',
               '| Δ ≥ | target | ε | τ | formed | the closest pair: distance, coefficients, offset (σ*) | coinciding pairs | max\\|w\\|/max\\|A\\| | '
               'residual |', '|---|---|---|---|---|---|---|---|---|']
        for r in G:
            for x in r['loops']:
                pa = x.get('pair')
                L_.append(f"| {r['dmin']:g} | {r['q']} | {x['eps']:g} | {x.get('tau', 400.0):g} | {x['formed']} | "
                          + ('—' if pa is None else f"{pa['d']:.2e}, {pa['w'][0]:+.3f}/{pa['w'][1]:+.3f}, {pa['off']:.3f}") + ' | '
                          + (', '.join(f"{a:+.3f}/{b:+.3f}" for a, b in x['pairs']) or 'none') + f" | {x['wmax']:.2f} | {x['rms']:.1e} |")
        L_.append('')
    if H is not None:
        L_ += ['## H. Two features of one sign, close', '',
               "Two features φ(·; c, σ*) of the given amplitudes (both positive) at the centres ±0.08 unless stated — closer than the kernels' "
               f"starting width σ₀ = {S0:g} —, the moving-width path (as C's, its budget {MAX_RHS_MOVING:.0e}): the outcome, the most kernels of "
               "the other sign in a block along the path and the most negative coefficient, the end's formed kernels of the other sign, the "
               "smallest formed width, the largest |w|, the residual:", '',
               '| target | outcome | negative kernels along the path (most negative w) | at the end | smallest width | max\\|w\\| | residual |',
               '|---|---|---|---|---|---|---|']
        for r in H:
            L_.append(f"| {r['name']} | {r['outcome']} | {r['neg_path']} ({r['wmin_path']:+.3f}) | {r['neg_end']} | {r['smin']:.2e} | "
                      f"{r['wmax']:.3g} | {r['moving']['rms']:.1e} |")
        L_.append('')
    L_ += ['Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'moving_width_raw.npz')
    md = os.path.join(OUT, 'moving_width.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D), ('E', part_E), ('F', part_F),
                  ('G', lambda: part_G(res['C'])), ('H', part_H)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
