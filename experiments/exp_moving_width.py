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
F. The real loop at small masses (ε = 10⁻⁶, 10⁻⁸, the free law on its curves, `rich_transit_nd.loop_nd` with moving widths) on C's
   first targets: its formed kernels against the path's, at TAUS_F times the path's last entry (to 150: referee round 61).
G. C's paths stopped by their budget: the stalled block's state when the budget ran out (the path as C's), its flow from the entry on
   the critical curve to s = 1e11 with the dipole law's constant (K = v²(E − E_inf), E_inf the optimal dipole fit's), and the real loop
   at the masses EPS_G on the grid TAUS_G (referee round 57): the dipole pair's half-difference v, centre of mass and net coefficient
   against the features, its three phases — the dipole's growth v³ = 3μ₀KΛ(ε)(τ − τ_d) until the next kernel forms, the landing, the
   collapse — and the landed v against (3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3} ∝ ε^{−1/6} (referee round 61).
H. Two features of one sign at ±0.08 (1.33σ₀ apart), symmetric and perturbed: the kernels of the other sign (referee round 57), the
   end's kernels narrower than the samples' spacing or of negative width; H2: the σ* = 0.06 paths to the tolerance TOL_H2; H3: the
   separation scan SEP_H at σ* = 0.06, equal amplitudes; HL: the real loop on H's two equal targets (referee round 61).
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
EPS_F, N_F, TAUS_F = (1e-6, 1e-8), 4, (1.3, 6.0, 30.0, 150.0)            # F: the loop to 150 times the last entry (round 61)
EPS_G = (1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)                         # G: the loop at six masses (the landing law: round 61)
TAUS_G = tuple(float(x) for x in np.unique(np.r_[np.linspace(10.0, 50.0, 21), np.geomspace(50.0, 6400.0, 36)].round(6)))  # 50·2^{k/5}
TAUS_G_MD = (400.0, 1600.0, 6400.0)                                     # G: the loop's states the md tabulates
S_G = (1e3, 1e5, 1e7, 1e9, 3.2e10, 1e11)                                # G: the block's flow to s = 1e11
LAND_TOL = 0.01                                                         # G: landed — the pair's centre of mass within 1% of a feature's
                                                                        # width of it and its net coefficient within 1% of its amplitude
SEP_H, TOL_H2 = (0.05, 0.1, 0.12, 0.15, 0.2), 1e-13                     # H3: σ* = 0.06, equal, at ±SEP_H; H2: the paths to this tolerance
EPS_HL, TAUS_HL = (1e-6, 1e-8), (1600.0, 6400.0)                        # HL: the loop on H's two equal targets
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


def outcome(m, As, cs=None, ss=None):
    """A moving-width path's outcome: 'exact' (as many formed kernels — |w| above 1% of the smallest amplitude — as features, within
    EXACT_TOL, the residual below 1e-8; idle entrants below 1% allowed), 'exact, null cluster' (each feature's nearest kernel within
    EXACT_TOL and the residual below 1e-8, with more formed kernels: the others sum to nothing on the samples — referee round 61),
    'coinciding pair' (the residual below 1e-8 with two formed kernels at one (c, σ), PAIR_TOL, whose net coefficient and centre of mass
    fit a feature within LAND_TOL — of its amplitude, of its width — when the features (cs, ss) are given: a split or a cancelling pair),
    'null cluster' (such a pair fitting no feature), 'cancellation' (the residual below 1e-8 with a coefficient beyond 1.5 times the
    largest amplitude or a width below the samples' spacing: several kernels cancelling), 'dipole run-off' (stopped by the budget with an
    opposite-sign pair in the stalled block, |w| beyond 10 times the largest amplitude and 0.05 apart: G), 'stopped by the budget' (another
    stall), or 'not fitted' (referee round 57)."""
    w = np.asarray(m['w'], float)
    As_ = np.asarray(As, float)
    A_ = np.abs(As_)
    formed = np.abs(w) > 1e-2 * A_.min()
    if m.get('oob'):
        st = m.get('stalled')
        if st is not None:
            ws = np.asarray(st['w'], float)
            cp = closest_pair(st['theta'], ws, opposite=True)
            if cp is not None and np.abs(ws).max() > 10 * A_.max() and cp[0] < 0.05:
                return 'dipole run-off'
        return 'stopped by the budget'
    th = np.atleast_2d(m['theta'])
    if m['rms'] < 1e-8 and max(m['err']) < EXACT_TOL:
        return 'exact' if formed.sum() == len(As) else 'exact, null cluster'
    if m['rms'] < 1e-8:
        thf, wf = th[formed], w[formed]
        pr = [p for p in pair_index(thf, wf)]
        if pr:
            fit = cs is None or any(pair_fits(thf, wf, i, j, cs, ss, As_) for i, j in pr)
            return 'coinciding pair' if fit else 'null cluster'
    if m['rms'] < 1e-8 and formed.any() and (np.abs(w[formed]).max() > 1.5 * A_.max() or np.abs(th[formed, 1]).min() < XS[1] - XS[0]):
        return 'cancellation'
    return 'not fitted'


def pair_index(theta, w):
    """The coinciding pairs among the kernels (`pairs_of`'s rule) by their indices."""
    th = np.atleast_2d(np.asarray(theta, float))
    out = []
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
            if abs(th[i, 0] - th[j, 0]) / sm + abs(np.log(abs(th[i, 1]) / abs(th[j, 1]))) < PAIR_TOL:
                out.append((i, j))
    return out


def pair_fits(theta, w, i, j, cs, ss, As):
    """Whether the pair (i, j) fits a feature: its net coefficient W within LAND_TOL of some feature's amplitude and its centre of mass
    within LAND_TOL of that feature's width from its centre (any feature: features may share a centre, as D's at d = 0)."""
    th = np.atleast_2d(np.asarray(theta, float))
    W = w[i] + w[j]
    if W == 0:
        return False
    cm = (w[i] * th[i, 0] + w[j] * th[j, 0]) / W
    return any(abs(cm - c_) <= LAND_TOL * s_ and abs(W - a_) <= LAND_TOL * abs(a_) for c_, s_, a_ in zip(cs, ss, As))


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
    path = sorted(int(j) for j in np.flatnonzero(R['act']))
    pp = closest_pair(R['c'][path], R['w'][path])
    snaps = []
    for eps in EPS_F:
        _, W, C, Lam, al = RN.loop_nd(L, K, f, theta0(), eps, taus, rtol=1e-9, atol=1e-12)
        for t in range(1, taus.size):
            big = np.abs(W[t]) > 1e-2 * np.abs(As).min()
            fj = np.flatnonzero(big)
            cp = closest_pair(C[t][fj], W[t][fj])
            snaps.append(dict(eps=eps, times=float(TAUS_F[t - 1]), formed=sorted(int(j) for j in fj),
                              dtheta=float(np.abs(C[t][path] - R['c'][path]).max()), dw=float(np.abs(W[t][path] - R['w'][path]).max()),
                              rest=float(np.abs(np.delete(W[t], path)).max()),
                              pair=None if cp is None else dict(d=cp[0], dc=float(abs(C[t][fj[cp[1]], 0] - C[t][fj[cp[2]], 0])),
                                                                ds=float(abs(C[t][fj[cp[1]], 1] - C[t][fj[cp[2]], 1])),
                                                                w=(float(W[t][fj[cp[1]]]), float(W[t][fj[cp[2]]]))),
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
    stalled = stalled_summary(Rc, f, K)
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
            loops.append(dict(eps=eps, tau=float(tau), Lam=float(Lam), theta=C[t][big].tolist(), w=W[t][big].tolist(),
                              wmax=float(np.abs(W[t]).max()), rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2)))))
    return dict(dmin=dmin, q=q, cs=list(cs), ss=list(ss), As=list(As), entries=n_ent,
                tau_d=float(R['events'][-1][0]) if R['events'] else float('nan'), stalled=stalled, flow=flow, law=law_,
                status=int(sol.status), loops=loops, secs=time.time() - t0)


def stalled_summary(Rc, f, K):
    """G: a budget-stopped path's stalled block — its kernels, largest |w| (the accepted polish's: solver-dependent), the closest
    opposite-sign pair's distance, the residual, the block's polishes — or None."""
    st = Rc.get('stalled')
    if st is None:
        return None
    cp = closest_pair(st['c'], st['w'], opposite=True)
    return dict(wmax=float(np.abs(st['w']).max()), sep=None if cp is None else cp[0], n=int(np.size(st['w'])),
                rms=float(np.sqrt(np.mean((f - K.synth(np.atleast_2d(st['c']), st['w'])) ** 2))),
                polish=Rc['polish_log'][-1] if Rc['polish_log'] else [], nfev=int(Rc['nfev']))


def dipole_exponents(flow):
    """The last two outputs' log-slopes of the largest |w| and of the opposite-sign pair's distance against s (the dipole: +1/3, −1/3 —
    on the free law's exponential branch the coefficient's mobility grows like |w| and the pair's energy excess falls like 1/w² at a
    fixed moment)."""
    a, b = flow[-2], flow[-1]
    if a['sep'] is None or b['sep'] is None:
        return None
    ls = np.log(b['s'] / a['s'])
    return float(np.log(b['wmax'] / a['wmax']) / ls), float(np.log(b['sep'] / a['sep']) / ls)


def dipole_pair(theta, w):
    """The loop's dipole pair among formed kernels (`dipole_law`'s choice): the largest |w| and the closest kernel of the other sign —
    (a, b, v = |w_a − w_b|/2, its distance |Δc|/σ + |Δlog σ|) or None."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    if w.size < 2:
        return None
    a = int(np.argmax(np.abs(w)))
    opp = [j for j in range(w.size) if w[j] * w[a] < 0]
    if not opp:
        return None
    dist = lambda j: (abs(th[a, 0] - th[j, 0]) / (0.5 * (abs(th[a, 1]) + abs(th[j, 1]))) + abs(np.log(abs(th[a, 1]) / abs(th[j, 1]))))
    b = min(opp, key=dist)
    return a, b, float(abs(w[a] - w[b]) / 2), float(dist(b))


def landing(theta, w, cs, ss, As):
    """The dipole pair against the features: its centre of mass's offset from the nearest feature (in that feature's width) and its net
    coefficient's relative difference from the feature's amplitude, (v, distance, offset, (W − A)/A) — landed if both within LAND_TOL."""
    dp = dipole_pair(theta, w)
    if dp is None:
        return None
    a, b, v, d = dp
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    W = w[a] + w[b]
    cm = (w[a] * th[a, 0] + w[b] * th[b, 0]) / W if W != 0 else 0.5 * (th[a, 0] + th[b, 0])
    k = int(np.argmin(np.abs(np.asarray(cs) - cm) / np.asarray(ss)))
    return dict(v=v, d=d, off=float(abs(cm - cs[k]) / ss[k]), dW=float((W - As[k]) / abs(As[k])), w=(float(w[a]), float(w[b])))


def phases(r, eps):
    """G: the loop at the mass eps on the dipole's target r, in three phases (referee round 61). The dipole: from the block's entry τ_d
    (the path's) the pair's v grows as v³ = 3μ₀KΛ(ε)(τ − τ_d) (the block's law in its fast time Λ(ε)τ) until the next kernel forms, at
    τ₃ (the first grid time with more formed kernels than the block's); the slopes of v³ over the grid times strictly between τ_d and
    τ₃ against 3μ₀KΛ (median, range). The landing: v at τ = 400 against the law's (3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3}, and whether the pair is
    landed (`landing`) at 400, 1600, 6400. The collapse: over 1600 → 6400 the exponents of the pair's distance and of the residual in τ."""
    xs = sorted([x for x in r['loops'] if x['eps'] == eps], key=lambda x: x['tau'])
    lw = r.get('law') or []
    if not xs or not lw or lw[-1] is None:
        return None
    K_ = lw[-1]['v'] ** 2 * (lw[-1]['E'] - lw[-1]['E_inf'])
    rate = 3 * lw[-1]['mu0'] * K_ * xs[0]['Lam']
    n0, td = r['entries'], r['tau_d']
    t3 = next((x['tau'] for x in xs if x['tau'] > td and len(x['w']) > n0), None)
    at = {x['tau']: x for x in xs}
    win = [x for x in xs if td < x['tau'] < (t3 if t3 is not None else np.inf) and len(x['w']) == n0]
    vs = [(x['tau'], dipole_pair(x['theta'], x['w'])) for x in win]
    vs = [(tau, dp[2]) for tau, dp in vs if dp is not None]
    sl = [(b[1] ** 3 - a[1] ** 3) / (b[0] - a[0]) / rate for a, b in zip(vs[:-1], vs[1:])]
    land = {tau: landing(at[tau]['theta'], at[tau]['w'], r['cs'], r['ss'], r['As']) if tau in at else None for tau in TAUS_G_MD}
    v400 = land[400.0]['v'] if land.get(400.0) else None
    vstar = (rate * (t3 - td)) ** (1 / 3) if t3 is not None else None
    col = None
    a, b = land.get(1600.0), land.get(6400.0)
    if a and b and a['d'] > 0 and b['d'] > 0:
        lt = np.log(6400.0 / 1600.0)
        col = (float(np.log(b['d'] / a['d']) / lt), float(np.log(at[6400.0]['rms'] / at[1600.0]['rms']) / lt))
    return dict(rate=rate, K=K_, t3=t3, slopes=sl, v400=v400, vstar=vstar, land=land, collapse=col)


def part_G(C):
    """G: C's targets whose moving-width path was stopped by its budget."""
    from multiprocessing import Pool
    jobs = [(r['dmin'], r['q'], (np.array(r['cs']), np.array(r['ss']), np.array(r['As'])), int(r['moving']['entries'])) for r in C
            if r['moving'].get('oob')]
    with Pool(N_WORKERS) as pool:
        return pool.map(g_job, jobs, chunksize=1)


def h_stats(ws, m, As):
    """H's statistics of a path — ws: each block's active coefficients — and of its end m (referee rounds 57 and 61): the most kernels
    of the other sign (below −1% of the smallest amplitude) in a block and the most negative coefficient; the end's formed kernels of the
    other sign and smallest formed width; over all the end's active kernels, those narrower than the samples' spacing and those of
    negative width (the widths are unconstrained); the largest |w|."""
    w, th = np.asarray(m['w'], float), np.atleast_2d(np.asarray(m['theta'], float))
    amin = min(abs(a) for a in As)
    formed = np.abs(w) > 1e-2 * amin
    return dict(neg_path=int(max([int((x < -1e-2 * amin).sum()) for x in ws], default=0)),
                wmin_path=float(min([x.min() for x in ws if x.size], default=0.0)), neg_end=int((formed & (w < 0)).sum()),
                smin=float(np.abs(th[formed, 1]).min()) if formed.any() else float('nan'),
                sub=int((np.abs(th[:, 1]) < XS[1] - XS[0]).sum()) if w.size else 0, negw=int((th[:, 1] < 0).sum()) if w.size else 0,
                wmax=float(np.abs(w).max()) if w.size else 0.0)


def clusters(m, cs, ss, As):
    """H2: the end's formed kernels by the feature nearest each (in its width): per feature their count and their mass over the
    feature's, Σw_i|σ_i|/(A_kσ_k) (one kernel fitting it exactly: 1)."""
    w, th = np.asarray(m['w'], float), np.atleast_2d(np.asarray(m['theta'], float))
    formed = np.flatnonzero(np.abs(w) > 1e-2 * min(abs(a) for a in As))
    near = [int(np.argmin(np.abs(np.asarray(cs) - th[i, 0]) / np.asarray(ss))) for i in formed]
    return [(sum(k == j for k in near), float(sum(w[i] * abs(th[i, 1]) for i, k in zip(formed, near) if k == j) / (As[j] * ss[j])))
            for j in range(len(cs))]


def h_job(args):
    """H: two features of one sign — the moving-width path's outcome and `h_stats` (referee round 57); with a tolerance (H2, referee
    round 61) the path run to it instead of 1e-9: the last entry's τ, the last blocks' residuals and the formed kernels per feature."""
    name, cs, ss, As = args[:4]
    tol = args[4] if len(args) > 4 else None
    t0 = time.time()
    f = feature_sum(cs, ss, As)
    if tol is None:
        Rm, _ = run_paths(f, fixed=False)
    else:
        Rm = RN.sliding_iss_nd(law(), MW.BumpsCW(XS, S0), f, theta0(), np.ones(C0.size), max_events=40, exit_kkt=True,
                               max_rhs=MAX_RHS_MOVING, tol=tol)
    m = match(Rm, cs, ss, As)
    ws = [h[2][h[2] != 0] for h in Rm['hist'][1:]]
    out = dict(name=name, cs=list(cs), ss=list(ss), As=list(As), moving=m, outcome=outcome(m, As, cs, ss), tol=tol,
               secs=time.time() - t0, **h_stats(ws, m, As))
    if tol is not None:
        out.update(last_tau=float(Rm['events'][-1][0]) if Rm['events'] else float('nan'),
                   rms_blocks=[float(h[3]) for h in Rm['hist'][-6:]], clusters=clusters(m, cs, ss, As))
    return out


def part_H():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, list(H_TARGETS), chunksize=1)


def part_H2():
    """H2: the σ* = 0.06 targets' paths to TOL_H2 (referee round 61: below 1e-9 the 'exact' ends were the path's stopping tolerance)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, [tuple(x) + (TOL_H2,) for x in H_TARGETS if x[0].startswith('σ* 0.06')], chunksize=1)


def h3_targets():
    return [(f"σ* 0.06, ±{s:g}", (-s, s), (0.06, 0.06), (1.0, 1.0)) for s in SEP_H]


def part_H3():
    """H3: the separation scan at σ* = 0.06, equal amplitudes (H's ±0.08 besides; referee round 61)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, h3_targets(), chunksize=1)


def hl_job(args):
    """HL: the real loop on one of H's targets at the masses EPS_HL to the times TAUS_HL: the formed kernels (|w| above 1% of the
    smallest amplitude), their states, the largest |w|, the residual (referee round 61)."""
    name, cs, ss, As = args
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    out = []
    for eps in EPS_HL:
        _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, TAUS_HL], rtol=1e-9, atol=1e-12)
        for t, tau in enumerate(TAUS_HL, start=1):
            big = np.flatnonzero(np.abs(W[t]) > 1e-2 * min(abs(a) for a in As))
            out.append(dict(eps=eps, tau=float(tau), theta=C[t][big].tolist(), w=W[t][big].tolist(), wmax=float(np.abs(W[t]).max()),
                            rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2)))))
    return dict(name=name, cs=list(cs), ss=list(ss), As=list(As), loops=out, secs=time.time() - t0)


def part_HL():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(hl_job, [x for x in H_TARGETS if x[0] in ('σ* 0.04, equal', 'σ* 0.06, equal')], chunksize=1)


def report(A=None, B=None, C=None, D=None, Ee=None, F=None, G=None, H=None, H2=None, H3=None, HL=None):
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
    def oc(r):
        return outcome(r['moving'], r['As'], r['cs'], r['ss'])

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
        L_ += ['', "The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude; exact with a null "
               "cluster where each feature's nearest kernel is exact and the others sum to nothing; a coinciding pair within "
               f"{PAIR_TOL:g} in |Δc|/σ + |Δlog σ| whose net coefficient and centre of mass fit a feature within {LAND_TOL:g}, a null "
               "cluster where it fits none; a budget stop a dipole run-off if the stalled block has an opposite-sign pair under 0.05 "
               "with |w| beyond 10 times the largest amplitude): " + '; '.join(
            f"Δ ≥ {d:g}: " + ', '.join(f"{len(cls.get((d, o), []))} {o}" for o in ('exact', 'exact, null cluster', 'coinciding pair', 'null cluster', 'cancellation', 'dipole run-off', 'stopped by the budget', 'not fitted') if cls.get((d, o)) or o == 'exact')
            for d in DELTAS_C) + '.']
        bad = [r for r in C if oc(r) != 'exact']

        def stall(m):
            st = m.get('stalled')
            if not m.get('oob') or st is None:
                return ''
            cp = closest_pair(st['theta'], st['w'], opposite=True)
            return (f"; the stalled block: {len(st['w'])} kernels, max|w| {np.abs(st['w']).max():.4g} (the accepted polish's, "
                    "solver-dependent)" + ('' if cp is None else f", its closest opposite pair {cp[0]:.2e} apart"))
        L_ += ['', "Not exact with moving widths (the end's state; for a budget stop the state before the stalled block, and the "
               "stalled block): " + ('; '.join(
            f"Δ ≥ {r['dmin']:g}, target {r['q']} (widths {', '.join(f'{v:.3f}' for v in r['ss'])}; centres "
            f"{', '.join(f'{v:+.2f}' for v in r['cs'])}; amplitudes {', '.join(f'{v:+.2f}' for v in r['As'])}): "
            f"{oc(r)}, {r['moving']['n']} kernels"
            + (f" (the pair: {split_text(r['moving'], r['As'])})" if split_text(r['moving'], r['As']) else '')
            + f", residual {r['moving']['rms']:.1e}" + stall(r['moving']) for r in bad) or 'none') + '. The width split '
            "(σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), the pair's third moment cancelled) against the observed.", '']
    if D is not None:
        L_ += ['## D. A narrow feature on a wide pedestal', '',
               "The pedestal 0.8φ(·; 0, 0.3) and the narrow feature 0.5φ(·; d, 0.06):", '',
               '| d | moving: formed kernels | outcome | max error (c, σ, w) | residual | fixed: kernels | residual |', '|---|---|---|---|---|---|---|']
        for r in D:
            m, f_ = r['moving'], r['fixed']
            As_ = [0.8, 0.5]
            sp = split_text(m, As_)
            oc_ = outcome(m, As_, [0.0, r['d']], [0.3, 0.06]) + (f" ({sp})" if sp else '')
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
               '(within 10⁻⁴ of it: how many, the largest deviation among them) '
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
                                                              f"({sum(abs(q - 1) < 1e-4 for q in qo)} of {len(qo)}" + (f", {max(abs(q - 1) for q in qo if abs(q - 1) < 1e-4):.0e})" if any(abs(q - 1) < 1e-4 for q in qo) else ')') if qo else '—') + ' | '
                      + (', '.join(f"{x.get('rms', np.nan):.2e}" + (' (budget)' if x.get('oob') else '') for x in ns) or '—') + ' |')
        L_.append('')
    if F is not None:
        L_ += ['## F. The real loop at a small mass', '',
               f"The free law with moving centres and widths from α = 0, masses ε = {', '.join(f'{e:g}' for e in EPS_F)} (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to "
               f"{', '.join(f'{t:g}' for t in TAUS_F)} times the path's last entry, on C's first {N_F} targets at Δ ≥ {DELTAS_C[1]:g}: its formed "
               "kernels (|w| above 1% of the smallest amplitude) against the path's, the path's kernels' largest deviation, the others' largest "
               "|w|, the closest pair among the formed kernels (|Δc|/σ + |Δlog σ|; its |Δc|, |Δσ|) and the residual; the path's closest pair "
               "and residual:", '',
               "| target | the path's kernels (closest pair; residual) | ε | times | formed by the loop | max \\|Δθ\\| | max \\|Δw\\| | the others' "
               "largest \\|w\\| | the loop's closest pair: distance (\\|Δc\\|, \\|Δσ\\|; coefficients) | residual |", '|---|---|---|---|---|---|---|---|---|---|']
        for r in F:
            for x in r.get('snaps', [dict(times=None, formed=r['formed'], dtheta=r['dtheta'], dw=r['dw'], rest=r['rest'], pair=None,
                                          rms=np.nan)]):
                pp = r.get('path_pair')
                pa = x.get('pair')
                L_.append(f"| {r['q']} | {r['path']} (" + ('—' if pp is None else f"{pp:.1e}") + f"; {r.get('path_rms', np.nan):.1e}) | "
                          + f"{x.get('eps', EPS_F[0]):g} | " + ('—' if x['times'] is None else f"{x['times']:g}") + f" | {x['formed']} | "
                          f"{x['dtheta']:.1e} | {x['dw']:.1e} | {x['rest']:.1e} | "
                          + ('—' if pa is None else f"{pa['d']:.1e} ({pa['dc']:.1e}, {pa['ds']:.1e}"
                             + (f"; {pa['w'][0]:+.3f}/{pa['w'][1]:+.3f}" if 'w' in pa else '') + ")") + f" | {x['rms']:.1e} |")
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
        L_ += g_loop_lines(G)
    if H is not None:
        L_ += ['## H. Two features of one sign, close', '',
               "Two features φ(·; c, σ*) of the given amplitudes (both positive) at the centres ±0.08 — 0.16 apart, 1.33 times the kernels' "
               f"starting width σ₀ = {S0:g} — unless stated, the moving-width path (as C's: its tolerance 1e-09, its budget "
               f"{MAX_RHS_MOVING:.0e}): the outcome (as C's), the most kernels of the other sign in a block along the path and the most "
               "negative coefficient, the end's formed kernels of the other sign, the smallest formed width, over all the end's active kernels "
               "those narrower than the samples' spacing and those of negative width (the widths are not constrained), the largest |w|, the "
               "residual:", '']
        L_ += h_table(H)
    if H2 is not None:
        L_ += ['### H2. The σ* = 0.06 paths to a tolerance of ' + f"{TOL_H2:g}", '',
               "The same paths with the tolerance lowered (the 1e-09 ends were the tolerance's): the outcome, the entries and the last one's "
               "τ, the residuals of the last blocks, the formed kernels by the feature nearest each (their count and mass over the feature's, "
               "Σw_i|σ_i|/(Aσ*)), and "
               "H's columns:", '',
               '| target | outcome | entries (last τ) | the last blocks\' residuals | formed per feature (count: mass) | negative along the path '
               '(most negative w) | at the end | smallest width | sub-sample, negative widths | max\\|w\\| |', '|---|---|---|---|---|---|---|---|---|---|']
        for r in H2:
            L_.append(f"| {r['name']} | {r['outcome']} | {r['moving']['entries']} ({r['last_tau']:.2e}) | "
                      + ', '.join(f"{v:.2e}" for v in r['rms_blocks']) + ' | ' + ', '.join(f"{n}: {s_:.4f}" for n, s_ in r['clusters'])
                      + f" | {r['neg_path']} ({r['wmin_path']:+.3f}) | {r['neg_end']} | {r['smin']:.2e} | {r['sub']}, {r['negw']} | "
                      f"{r['wmax']:.3g} |")
        L_.append('')
    if H3 is not None:
        L_ += ['### H3. The separation scan', '',
               "σ* = 0.06, equal amplitudes 1, at ±c (2c in units of σ₀ in brackets; H's ±0.08 is 1.33σ₀), H's columns:", '']
        L_ += h_table(sorted(list(H3) + [r for r in (H if H is not None else []) if r['name'] == 'σ* 0.06, equal'],
                             key=lambda r: r['cs'][1]), sep=True)
    if HL is not None:
        L_ += ['### HL. The real loop on H\'s equal targets', '',
               "The free law from α = 0 at the masses ε (`loop_nd`, rtol 10⁻⁹) at the times τ: the formed kernels (|w| above 1% of the "
               "smallest amplitude), of the other sign, the smallest formed width, the formed kernels centred outside the samples' window, "
               "the largest |w|, the residual:", '',
               '| target | ε | τ | formed | of the other sign (most negative w) | smallest width | outside the window | max\\|w\\| | residual |',
               '|---|---|---|---|---|---|---|---|---|']
        for r in HL:
            for x in r['loops']:
                w, th = np.asarray(x['w'], float), np.atleast_2d(np.asarray(x['theta'], float)).reshape(-1, 2)
                L_.append(f"| {r['name']} | {x['eps']:g} | {x['tau']:g} | {w.size} | {int((w < 0).sum())}"
                          + (f" ({w.min():+.3f})" if (w < 0).any() else '') + f" | " + (f"{np.abs(th[:, 1]).min():.3f}" if w.size else '—')
                          + f" | {int((np.abs(th[:, 0]) > XS[-1]).sum())} | {x['wmax']:.3g} | {x['rms']:.1e} |")
        L_.append('')
    L_ += ['Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def h_table(H, sep=False):
    """H's table (H and H3)."""
    out = [('| separation | ' if sep else '| target | ') + 'outcome | negative kernels along the path (most negative w) | at the end | '
           'smallest width | sub-sample, negative widths (all active) | max\\|w\\| | residual |', '|---|---|---|---|---|---|---|---|']
    for r in H:
        lab = f"±{r['cs'][1]:g} ({(r['cs'][1] - r['cs'][0]) / S0:.2f}σ₀)" if sep else r['name']
        out.append(f"| {lab} | {r['outcome']} | {r['neg_path']} ({r['wmin_path']:+.3f}) | {r['neg_end']} | {r['smin']:.2e} | "
                   f"{r.get('sub', '—')}, {r.get('negw', '—')} | {r['wmax']:.3g} | {r['moving']['rms']:.1e} |")
    return out + ['']


def g_loop_lines(G):
    """G's loop tables (referee round 61): the states at TAUS_G_MD — the dipole pair (`dipole_pair`) with its distance, coefficients,
    centre of mass's offset from the nearest feature (in its width) and net coefficient against the feature's amplitude, the coinciding
    pairs, the largest |w| over the largest amplitude, the residual —; the three phases (`phases`); the landed v against ε."""
    if any('theta' not in x for r in G for x in r['loops']):
        return ['', "(G's loops of a run before referee round 61, without the loop's states: rerun G.)", '']
    L_ = ['', "The real loop at the masses ε from α = 0 (`loop_nd`, rtol 10⁻⁹) at τ = " + ', '.join(f"{t:g}" for t in TAUS_G_MD) + ": the "
          "formed kernels, the dipole pair (the largest |w| and the closest kernel of the other sign: its distance |Δc|/σ + |Δlog σ|, "
          "coefficients, its centre of mass's offset from the nearest feature in that feature's width and its net coefficient W against the "
          f"feature's amplitude, (W − A)/|A|; landed (L) if both within {LAND_TOL:g}), the coinciding pairs (within {PAIR_TOL:g}), the "
          "largest |w| over the largest amplitude, the residual:", '',
          '| Δ ≥ | target | ε | τ | formed | the dipole pair: distance, coefficients, offset (σ*), (W − A)/\\|A\\| | coinciding pairs | '
          'max\\|w\\|/max\\|A\\| | residual |', '|---|---|---|---|---|---|---|---|---|']
    for r in G:
        for x in sorted(r['loops'], key=lambda x: (x['eps'], x['tau'])):
            if x['tau'] not in TAUS_G_MD:
                continue
            lp = landing(x['theta'], x['w'], r['cs'], r['ss'], r['As'])
            th = np.atleast_2d(np.asarray(x['theta'], float)).reshape(-1, 2)
            L_.append(f"| {r['dmin']:g} | {r['q']} | {x['eps']:g} | {x['tau']:g} | {len(x['w'])} | "
                      + ('—' if lp is None else f"{lp['d']:.2e}, {lp['w'][0]:+.3f}/{lp['w'][1]:+.3f}, {lp['off']:.1e}, {lp['dW']:+.1e}"
                         + (' L' if lp['off'] <= LAND_TOL and abs(lp['dW']) <= LAND_TOL else '')) + ' | '
                      + (', '.join(f"{a:+.3f}/{b:+.3f}" for a, b in pairs_of(th, x['w'])) or 'none')
                      + f" | {x['wmax'] / np.abs(r['As']).max():.2f} | {x['rms']:.1e} |")
    spl = []                                                           # (ii)(a)'s width split on the landed pairs (referee round 61)
    for r in G:
        for x in r['loops']:
            lp = landing(x['theta'], x['w'], r['cs'], r['ss'], r['As']) if x['tau'] in TAUS_G_MD else None
            if lp is None or lp['off'] > LAND_TOL or abs(lp['dW']) > LAND_TOL:
                continue
            th, w = np.atleast_2d(np.asarray(x['theta'], float)).reshape(-1, 2), np.asarray(x['w'], float)
            a, b = dipole_pair(th, w)[:2]
            dc, ds, pr = width_split(th, w, a, b)
            spl.append((ds / pr, abs(w[a] - w[b]) * abs(dc) / (abs(w[a] + w[b]) * 0.5 * (abs(th[a, 1]) + abs(th[b, 1])))))
    if spl:
        rt, pm = np.array([s_[0] for s_ in spl]), np.array([s_[1] for s_ in spl])
        L_ += ['', f"(ii)(a)'s width split on the landed pairs (L, {len(spl)} rows): the observed σ_b − σ_a over (w_b − w_a)Δc²/(6Wσ) "
               f"{rt.min():.4f}–{rt.max():.4f} (median {np.median(rt):.4f}), at |w_a − w_b||Δc|/(|W|σ) = {pm.min():.2e}–{pm.max():.2e}; "
               f"its largest deviation {np.abs(rt - 1).max():.1e} at the parameter {pm[np.argmax(np.abs(rt - 1))]:.2e}."]
    L_ += ['', "The three phases (`phases`): the dipole block's entry τ_d (the path's) and the next kernel's formation τ₃ in the loop; over "
           "the grid times between them the slopes of v³ (v the dipole pair's half-difference) against 3μ₀KΛ(ε) — the block's law in its "
           "fast time Λ(ε)τ, K and μ₀ from the flow's last output — (median [range]); the landed v (at τ = 400) against the law's "
           "(3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3}; landed at τ = " + ', '.join(f"{t:g}" for t in TAUS_G_MD) + "; over 1600 → 6400 the exponents in τ of "
           "the pair's distance and of the residual:", '',
           '| Δ ≥ | target | ε | τ_d | τ₃ | d(v³)/dτ / 3μ₀KΛ | v at 400 | the law\'s v | landed | exponents: distance, residual |',
           '|---|---|---|---|---|---|---|---|---|---|']
    law_v, phs = {}, {}
    for r in G:
        for eps in EPS_G:
            ph = phases(r, eps)
            if ph is None:
                continue
            sl = ph['slopes']
            law_v[(r['dmin'], r['q'], eps)] = ph['v400']
            phs[(r['dmin'], r['q'], eps)] = ph
            L_.append(f"| {r['dmin']:g} | {r['q']} | {eps:g} | {r['tau_d']:.4g} | " + ('—' if ph['t3'] is None else f"{ph['t3']:g}") + ' | '
                      + (f"{np.median(sl):.3f} [{min(sl):.3f}, {max(sl):.3f}]" if sl else '—') + ' | '
                      + ('—' if ph['v400'] is None else f"{ph['v400']:.4g}") + ' | ' + ('—' if ph['vstar'] is None else f"{ph['vstar']:.4g}")
                      + ' | ' + ', '.join('—' if ph['land'][t] is None else ('yes' if ph['land'][t]['off'] <= LAND_TOL
                                                                                and abs(ph['land'][t]['dW']) <= LAND_TOL else 'no')
                                          for t in TAUS_G_MD) + ' | '
                      + ('—' if ph['collapse'] is None else f"{ph['collapse'][0]:+.3f}, {ph['collapse'][1]:+.3f}") + ' |')
    L_ += ['', "The landed v against the mass: v at τ = 400 and its ratio per decade of ε (10^{1/6} = " + f"{10 ** (1 / 6):.4f} if v ∝ "
           "ε^{−1/6}):", '', '| Δ ≥ | target | v at ε = ' + ', '.join(f"{e:g}" for e in EPS_G) + ' | ratios per decade |', '|---|---|---|---|']
    for r in G:
        vs = [law_v.get((r['dmin'], r['q'], e)) for e in EPS_G]
        rat = [f"{b / a:.4f}" if a and b else '—' for a, b in zip(vs[:-1], vs[1:])]
        L_.append(f"| {r['dmin']:g} | {r['q']} | " + ', '.join('—' if v is None else f"{v:.4g}" for v in vs) + f" | {', '.join(rat)} |")
    L_ += g_phase_summary(phs, law_v)
    return L_ + ['']


def g_phase_summary(phs, law_v):
    """G's three phases summed up (referee round 61): per ε the cells landed at τ = 6400; at the smallest ε the dipole phase's slope over
    the law (the cells' medians) and the landed v over the law's; the last decade's ratio of v; on the cells landed at both 1600 and 6400
    the collapse's exponents (distance, residual) and the frozen ones (|exponent| below 10⁻³)."""
    if not phs:
        return []
    ok = lambda lp: lp is not None and lp['off'] <= LAND_TOL and abs(lp['dW']) <= LAND_TOL
    e0 = min(EPS_G)
    n6 = {e: (sum(ok(ph['land'].get(6400.0)) for k, ph in phs.items() if k[2] == e), sum(k[2] == e for k in phs)) for e in EPS_G}
    sl = [np.median(ph['slopes']) for k, ph in phs.items() if k[2] == e0 and ph['slopes']]
    vr = [ph['v400'] / ph['vstar'] for k, ph in phs.items() if k[2] == e0 and ph['v400'] and ph['vstar']]
    dec = [law_v[k[:2] + (e0,)] / law_v[k[:2] + (10 * e0,)] for k in phs if k[2] == e0 and law_v.get(k[:2] + (10 * e0,)) and law_v.get(k)]
    both = [(k, ph['collapse']) for k, ph in phs.items() if ok(ph['land'].get(1600.0)) and ok(ph['land'].get(6400.0)) and ph['collapse']]
    fr = [k for k, c in both if abs(c[0]) < 1e-3]
    col = np.array([c for k, c in both if abs(c[0]) >= 1e-3])
    span = lambda v, f_: f"{min(v):{f_}}–{max(v):{f_}}" if len(v) else '—'
    return ['', "Over the six targets: landed at τ = 6400 in " + ', '.join(f"{a} of {b} at ε = {e:g}" for e, (a, b) in n6.items())
            + f"; at ε = {e0:g} the dipole phase's slope over the law {span(sl, '.3f')} (the cells' medians, {len(sl)} cells with a "
            f"phase), the landed v over the law's {span(vr, '.2f')}, v's ratio over the last decade {span(dec, '.4f')}; on the "
            f"{len(both)} cells landed at both 1600 and 6400, {len(fr)} frozen (" + ', '.join(f"{k[1]} at {k[2]:g}" for k in fr)
            + ")" + (f", the others collapsing with the exponents {col[:, 0].min():+.3f} … {col[:, 0].max():+.3f} (distance) and "
                     f"{col[:, 1].min():+.3f} … {col[:, 1].max():+.3f} "
                     "(residual)" if len(col) else '') + '.']


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'moving_width_raw.npz')
    md = os.path.join(OUT, 'moving_width.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'H2', 'H3', 'HL')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D), ('E', part_E), ('F', part_F),
                  ('G', lambda: part_G(res['C'])), ('H', part_H), ('H2', part_H2), ('H3', part_H3), ('HL', part_HL)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
