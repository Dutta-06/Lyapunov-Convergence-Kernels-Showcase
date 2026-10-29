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
EPS_HL, TAUS_HL = (1e-6, 1e-8), (1600.0, 6400.0, 2.5e4, 1e5)            # HL: the loop on H's targets (referee round 68: to 1e5)
HL_NAMES = ('σ* 0.04, equal', 'σ* 0.06, equal', 'σ* 0.06, A₂ 0.99', 'σ* 0.06, A₂ 0.95', 'σ* 0.06, A₂ 0.9', 'σ* 0.06, c₂ +0.01',
            'σ* 0.06, σ₂ 0.065')                                        # HL: the symmetric targets and the asymmetric ones
HP_NAMES = ('σ* 0.04, equal', 'σ* 0.06, equal', 'σ* 0.06, c₂ +0.01', 'σ* 0.06, A₂ 0.9')   # HP: H's path without the polish
TAUS_FINE = tuple(float(x) for x in np.round(np.arange(10.0, 140.05, 0.1), 6))  # G (referee round 68): the formations on a fine grid
EPS_GX = (1e-11, 1e-12)                                                 # G: two smaller masses (the fine phases only)
P_LAND = 0.28                                                           # G: (ii)(a)'s parameter p = |w_a − w_b||Δc|/(|W|σ) at the landing
SEQ_TOL = 0.5                                                           # G: the loop's intermediate slope over the path's law beyond 1 ± it:
                                                                        # another sequence (referee round 80: 27 at ε ≥ 10⁻⁹, 18 at 10⁻⁷)
CLOSURE_EXACT = 1e-12                                                   # GT: a closure fit whose residual's rms is below it is exact
                                                                        # (referee round 80: target 20's K′ and τ₄ were roundoff)
Q_GL, EPS_GL = (20, 9), 1e-10                                           # GL: two dipole targets' collapse (9: referee round 74), to TAUS_GL
TAUS_GL = tuple(float(x) for x in np.geomspace(400.0, 1e7, 29))
RATIOS_A2, WMUL_A2 = (0.25, 1.0, 4.0), (2.5, 3.0, 5.0, 10.0, 20.0)      # A2: one kernel from (c*, σ*) with w = WMUL_A2·A (rounds 74, 80)
S_A2 = tuple(float(x) for x in np.geomspace(1.0, 1e13, 27))             # A2: the block's flow's outputs (referee round 80: to 10¹³)
SPIKE_A2, ESC_A2 = 0.2, 0.5                                             # A2: on the spike σ/σ* below it; escaped once back above ESC_A2
HLL = (('σ* 0.06, equal', 1e-6), ('σ* 0.06, equal', 1e-8), ('σ* 0.06, A₂ 0.99', 1e-8), ('σ* 0.06, A₂ 0.95', 1e-6))
TAUS_HLL = (1e5, 1e6, 1e7, 1e8, 1e9)                                    # HLL: HL's cells to 10⁹ (referee rounds 74, 80)
VIS = 1e-10                                                             # visible on the samples: the largest φ there at least VIS
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


def moving_state(r):
    """C's moving-width path's formed kernels and residual — at its stalled block (G) where its budget stopped it (referee round 68: the
    table printed the block before it)."""
    m = r['moving']
    st = m.get('stalled')
    if st is None:
        return n_formed(m, r['As']), m['rms']
    w = np.asarray(st['w'], float)
    f = feature_sum(r['cs'], r['ss'], r['As'])
    rms = float(np.sqrt(np.mean((f - MW.BumpsCW(XS, S0).synth(np.atleast_2d(np.asarray(st['theta'], float)), w)) ** 2)))
    return int(np.sum(np.abs(w) > 1e-2 * min(abs(a) for a in r['As']))), rms


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


def a2_job(args):
    """A2 (referee round 74): one kernel's block on one feature (A's), started at (c*, σ*) with the coefficient w = wmul·A — E(start) =
    (1 − wmul)²·½⟨f*²⟩, below ½⟨f*²⟩ (i)'s condition) only for 0 < wmul < 2 —: its flow on the critical curve (LSODA, the exact
    Jacobian) to S_A2[-1], at the outputs S_A2 its w/A, σ/σ*, (c − c*)/σ* and E over ½⟨f*²⟩."""
    ratio, wmul = args
    from scipy.integrate import solve_ivp
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    ss = ratio * S0
    f = feature_sum([C_A], [ss], [A_A])
    w0 = np.array([wmul * A_A])
    z0 = np.r_[np.sign(w0) * RT._alpha_of(L, np.abs(w0), RT.coef_crit), C_A, ss]
    sol = solve_ivp(RN.critical_rhs_nd(L, K, f), (0.0, S_A2[-1]), z0, method='LSODA', jac=RN.critical_jac_nd(L, K, f), rtol=1e-10,
                    atol=1e-14, t_eval=S_A2)
    e0 = float(np.mean(f ** 2) / 2)
    rows = []
    for i in range(sol.y.shape[1]):
        w, th = RT.critical_w(L, sol.y[:1, i]), sol.y[1:, i].reshape(1, 2)
        rows.append(dict(s=float(sol.t[i]), w=float(w[0]) / A_A, sig=float(th[0, 1]) / ss, dc=float(th[0, 0] - C_A) / ss,
                         E=float(np.mean((f - K.synth(th, w)) ** 2) / 2) / e0))
    return dict(ratio=ratio, wmul=wmul, E0=(1 - wmul) ** 2, rows=rows, status=int(sol.status), secs=time.time() - t0)


def part_A2():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a2_job, [(r_, m_) for r_ in RATIOS_A2 for m_ in WMUL_A2], chunksize=1)


def a2_spike(rows):
    """A2's spike (referee round 80): the outputs on it (σ/σ* below SPIKE_A2), their E over ½⟨f*²⟩, and the first output after them back
    above ESC_A2 (the escape; None if none by the last output): (first s on it, last s on it, E's range, the escape's s) or None."""
    on = [i for i, x in enumerate(rows) if x['sig'] < SPIKE_A2]
    if not on:
        return None
    esc = next((rows[i]['s'] for i in range(on[-1] + 1, len(rows)) if rows[i]['sig'] > ESC_A2), None)
    es = [rows[i]['E'] for i in on]
    return rows[on[0]]['s'], rows[on[-1]]['s'], (min(es), max(es)), esc


def a2_lines(A2):
    """A2's table (referee rounds 74, 80): one kernel from (c*, σ*) with a large coefficient — its end, the smallest σ/σ*, and the spike:
    the outputs with σ/σ* below SPIKE_A2 (the width collapsed onto the nearest samples, J_N's σ → 0 limit), its E, and the escape."""
    if A2 is None or len(A2) == 0:                                      # raw-loaded: an object array
        return []
    out = ['### A2. One kernel from a large coefficient (referee rounds 74, 80)', '',
           f"A's feature at σ*/σ₀ = " + ', '.join(f"{r:g}" for r in RATIOS_A2) + ": one kernel at (c*, σ*) with w = " + ', '.join(
               f"{m:g}" for m in WMUL_A2) + f"·A — E(start) (1 − w/A)²·½⟨f*²⟩, (i)'s condition E(start) < ½‖f*‖² holding for w < 2A only —, "
           f"its block's flow on the critical curve to s = {S_A2[-1]:.0e}: at the end w/A, σ/σ*, (c − c*)/σ* and E over ½⟨f*²⟩; the "
           f"smallest σ/σ* on the way (at the outputs); the spike — the outputs with σ/σ* below {SPIKE_A2:g} —, its E over ½⟨f*²⟩, and the "
           f"first output after it with σ/σ* back above {ESC_A2:g} (the escape):", '',
           '| σ*/σ₀ | w(start)/A | E(start)/½⟨f*²⟩ | end: w/A, σ/σ*, (c − c*)/σ*, E/½⟨f*²⟩ | smallest σ/σ* | on the spike: s, E/½⟨f*²⟩ | '
           'escaped at s | status |', '|---|---|---|---|---|---|---|---|']
    for r in A2:
        e = r['rows'][-1] if r['rows'] else None
        sp = a2_spike(r['rows']) if r['rows'] else None
        out.append(f"| {r['ratio']:g} | {r['wmul']:g} | {r['E0']:g} | " + ('—' if e is None else
                   f"{e['w']:.4f}, {e['sig']:.4g}, {e['dc']:+.1e}, {e['E']:.2e} (s = {e['s']:.0e})") + ' | '
                   + (f"{min(x['sig'] for x in r['rows']):.3g}" if r['rows'] else '—') + ' | '
                   + ('—' if sp is None else f"{sp[0]:.0e} … {sp[1]:.0e}, {sp[2][0]:.3f}–{sp[2][1]:.3f}") + ' | '
                   + ('—' if sp is None else ('not by the end' if sp[3] is None else f"{sp[3]:.0e}")) + f" | {r['status']} |")
    return out + ['']


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
    by least squares from the state (`closure_fit`) —: E_inf = ½⟨r²⟩ there. Returns dict(v, E, E_inf, rms_inf, p, mu0) or None."""
    th, w = np.atleast_2d(np.asarray(th, float)), np.asarray(w, float)
    a = int(np.argmax(np.abs(w)))
    opp = [j for j in range(w.size) if w[j] * w[a] < 0]
    if not opp:
        return None
    dist = lambda j: (abs(th[a, 0] - th[j, 0]) / (0.5 * (abs(th[a, 1]) + abs(th[j, 1]))) + abs(np.log(abs(th[a, 1]) / abs(th[j, 1]))))
    b = min(opp, key=dist)
    oth = [i for i in range(w.size) if i not in (a, b)]
    tb = 0.5 * (th[a] + th[b])
    p0 = w[a] * (th[a] - tb) + w[b] * (th[b] - tb)
    cf = closure_fit(K, f, th[oth], w[oth], w[a] + w[b], tb, p0)
    E_inf = float(np.mean(cf['r'] ** 2) / 2)
    return dict(v=float(abs(w[a] - w[b]) / 2), E=float(np.mean((f - K.synth(th, w)) ** 2) / 2), E_inf=E_inf,
                rms_inf=float(np.sqrt(2 * E_inf)), p=cf['p'].tolist(), mu0=float(4 * L.G * L.c / L.A),
                dc=float(abs(th[a, 0] - th[b, 0])), ds=float(abs(th[a, 1] - th[b, 1])))


def closure_fit(K, f, th_o, w_o, W0, tb0, p0, max_nfev=20000):
    """The dipole closure's least-squares fit (referee rounds 57, 74): the other kernels (θ, w), Wφ(θ̄) and p·∇_θφ(θ̄) by
    Levenberg–Marquardt from the given start. Returns dict(r — the fit's residual —, tb, p, th, w — the other kernels, their widths as
    |σ| (referee round 80: φ is even in σ) —, W)."""
    from scipy.optimize import least_squares
    th_o, w_o = np.asarray(th_o, float).reshape(-1, 2), np.asarray(w_o, float).reshape(-1)
    m = w_o.size

    def model(z):
        g = K.synth(z[m:3 * m].reshape(m, 2), z[:m]) if m else 0.0
        t_ = z[3 * m + 1:3 * m + 3]
        return g + z[3 * m] * K.phi(t_[None, :])[:, 0] + K.dphi(t_[None, :])[:, 0, :] @ z[3 * m + 3:]
    z0 = np.r_[w_o, th_o.ravel(), W0, tb0, p0]
    sol = least_squares(lambda z: (f - model(z)) / np.sqrt(f.size), z0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15,
                        max_nfev=max_nfev)
    z = sol.x
    th = z[m:3 * m].reshape(m, 2).copy()
    th[:, 1] = np.abs(th[:, 1])
    return dict(r=f - model(z), tb=z[3 * m + 1:3 * m + 3], p=z[3 * m + 3:], th=th, w=z[:m], W=float(z[3 * m]))


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
    fine, loops = {}, []
    tg = np.unique(np.r_[TAUS_G, TAUS_FINE])                            # one run per mass: the extra outputs leave the steps alone
    for eps in EPS_G + EPS_GX:
        _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, tg], rtol=1e-9, atol=1e-12)
        sts = []
        for t, tau in enumerate(tg, start=1):
            big = np.flatnonzero(np.abs(W[t]) > 1e-2 * np.abs(As).min())
            sts.append((float(tau), C[t][big], W[t][big]))
            if eps in EPS_G and tau in TAUS_G:
                loops.append(dict(eps=eps, tau=float(tau), Lam=float(Lam), theta=C[t][big].tolist(), w=W[t][big].tolist(),
                                  wmax=float(np.abs(W[t]).max()), rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2)))))
        fine[eps] = fine_phases(sts, n_ent, float(R['events'][-1][0]) if R['events'] else float('nan'), cs, ss, As)
    return dict(dmin=dmin, q=q, cs=list(cs), ss=list(ss), As=list(As), entries=n_ent,
                tau_d=float(R['events'][-1][0]) if R['events'] else float('nan'), stalled=stalled, flow=flow, law=law_,
                status=int(sol.status), loops=loops, fine=fine, secs=time.time() - t0)


def path_duals_at(K, f, R):
    """The path's duals at the dipole block's entry τ_d (referee round 68): `path_duals_nd` gives them at the path's history times, and
    the path stops at τ_d, which has none, so they are advanced from the last history time with that interval's residual. Returns (u,
    act, the dictionary's θ, τ_d)."""
    P = RN.path_duals_nd(K, f, R)
    td = float(R['events'][-1][0])
    act = np.asarray(R['act'], bool)
    cen = np.asarray(P['centres'][-1])
    return P['U'][-1] + np.where(act, 0.0, (td - P['taus'][-1]) * K.corr(cen, P['res'][-1])[0]), act, cen, td


def first_wall(u, act, v):
    """The first inactive kernel whose dual u, moving at the correlations v, reaches its wall ±1: (its index, the time) or (None, inf)."""
    ina = ~np.asarray(act, bool) & (v != 0)
    if not ina.any():
        return None, float('inf')
    dt = np.full(v.size, np.inf)
    dt[ina] = (1 - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
    j = int(np.argmin(dt))
    return j, float(dt[j])


def t3_from_duals(K, f, R, theta, w):
    """τ₃ from the sliding path alone (referee round 68): the path's duals at the dipole block's entry τ_d (`path_duals_at`) moved at
    the correlations of the block's residual at the end of its flow on the critical curve (θ, w at s = 10¹¹: the dipole's asymptotic
    state, which the loop reaches in Λ(ε)τ); the first inactive kernel to reach its wall ±1 (`first_wall`)."""
    u, act, cen, td = path_duals_at(K, f, R)
    j, dt = first_wall(u, act, K.corr(cen, f - K.synth(np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)))[0])
    return None if j is None else float(td + dt)


def d3phi(tb, p, x=XS):
    """∇³_θφ(x; θ̄):p³, the third derivative of φ(x; θ̄ + tp) in t at 0 (referee round 74; analytic): with u = x − c and
    h = −u²/(2σ²) along the line, h′ = p_cu/σ² + p_σu²/σ³, h″ = −p_c²/σ² − 4p_cp_σu/σ³ − 3p_σ²u²/σ⁴, h‴ = 6p_c²p_σ/σ³ +
    18p_cp_σ²u/σ⁴ + 12p_σ³u²/σ⁵, and (e^h)‴ = (h‴ + 3h′h″ + h′³)e^h."""
    (c, sg), (pc, ps) = np.asarray(tb, float), np.asarray(p, float)
    u = np.asarray(x, float) - c
    h1 = pc * u / sg ** 2 + ps * u ** 2 / sg ** 3
    h2 = -pc ** 2 / sg ** 2 - 4 * pc * ps * u / sg ** 3 - 3 * ps ** 2 * u ** 2 / sg ** 4
    h3 = 6 * pc ** 2 * ps / sg ** 3 + 18 * pc * ps ** 2 * u / sg ** 4 + 12 * ps ** 3 * u ** 2 / sg ** 5
    return (h3 + 3 * h1 * h2 + h1 ** 3) * np.exp(-u ** 2 / (2 * sg ** 2))


def k_third(r, tb, p):
    """The dipole's energy coefficient from the third order (referee round 74): at the optimal dipole fit (residual r, midpoint θ̄,
    moment p) E(v) = E_inf + K/v² + O(v⁻⁴) with K = −⟨r, ∇³φ:p³⟩/24 (⟨·⟩ the samples' mean, as E = ½⟨r²⟩). Referee round 80: O(v⁻⁴),
    not O(v⁻³) — about the midpoint the pair's sum has only even powers of 1/v (W multiplies the even orders of Δθ = p/v, w_a − w_b = 2v
    the odd ones)."""
    return float(-np.mean(np.asarray(r, float) * d3phi(tb, p)) / 24)


def kappa_from_path(K, f, R, theta, w):
    """The dipole's intermediate phase from the sliding path alone (referee round 74). At the dipole block's asymptotic state (θ, w: its
    flow's end at s = 10¹¹) the optimal dipole fit (`closure_fit`) and K (`k_third`); τ₃ and its entrant from the duals at τ_d moved by
    that state's residual (`t3_from_duals`' rule); the entrant added at its seed with the coefficient 10⁻³ of its correlation's sign
    and the closure refitted from the first fit — the three-kernel block's optimal dipole fit and its K′ by the same formula —; the
    duals moved on from τ₃ at that fit's residual: τ₄, the first to reach a wall. While the three-kernel block holds, its pair follows
    the block's law with K′: v³ = 3μ₀Λ(ε)(K(τ₃ − τ_d) + K′(τ − τ₃)), so v_land = κv(τ₃) with κ³ = 1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)).
    Referee round 80: the three-kernel block's closure rms (rms2) — below CLOSURE_EXACT the block fits exactly: K′ = 0, no τ₄, κ = 1
    (exact) —; for K′ < 0 the pair retracts and v³ would reach 0 at t_diss = τ₃ + K(τ₃ − τ_d)/|K′| (the dipole dissolves unless τ₄ comes
    first); κ only where K > 0 and the bracket is positive.
    Returns dict(K_third, Kp, t3_path, t4_path, kappa_path, t_diss, rms2, exact), None where undefined."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    out = dict(K_third=None, Kp=None, t3_path=None, t4_path=None, kappa_path=None, t_diss=None, rms2=None, exact=False)
    dp = dipole_pair(th, w)
    if dp is None:
        return out
    a, b = dp[:2]
    oth = [i for i in range(w.size) if i not in (a, b)]
    tb = 0.5 * (th[a] + th[b])
    c1 = closure_fit(K, f, th[oth], w[oth], w[a] + w[b], tb, w[a] * (th[a] - tb) + w[b] * (th[b] - tb))
    out['K_third'] = k1 = k_third(c1['r'], c1['tb'], c1['p'])
    u, act, cen, td = path_duals_at(K, f, R)
    v1 = K.corr(cen, f - K.synth(th, w))[0]
    j3, d3 = first_wall(u, act, v1)
    if j3 is None:
        return out
    out['t3_path'] = td + d3
    u = u + np.where(act, 0.0, d3 * v1)
    act = act.copy()
    act[j3], u[j3] = True, np.sign(v1[j3])
    c2 = closure_fit(K, f, np.r_[c1['th'], cen[j3][None, :]], np.r_[c1['w'], np.sign(v1[j3]) * 1e-3], c1['W'], c1['tb'], c1['p'],
                     max_nfev=50000)
    out['Kp'] = k2 = k_third(c2['r'], c2['tb'], c2['p'])
    out['rms2'] = float(np.sqrt(np.mean(c2['r'] ** 2)))
    if out['rms2'] < CLOSURE_EXACT:                                    # referee round 80: an exact fit — K′ and τ₄ are roundoff
        out.update(Kp=0.0, exact=True, kappa_path=1.0)
        return out
    if k1 > 0 and k2 < 0:                                              # referee round 80: the pair retracts — v³ = 0 at t_diss
        out['t_diss'] = out['t3_path'] + k1 * d3 / abs(k2)
    j4, d4 = first_wall(u, act, K.corr(cen, c2['r'])[0])
    if j4 is None:
        return out
    out['t4_path'] = out['t3_path'] + d4
    br = 1 + k2 * d4 / (k1 * d3)
    out['kappa_path'] = float(np.cbrt(br)) if k1 > 0 and br > 0 else None   # referee round 80: no κ past the dissolution
    return out


def gt_job(args):
    """GT: `t3_from_duals` and `kappa_from_path` (referee round 74) for one of G's targets, its path rerun as `g_job`'s and its block's
    flow end from G."""
    q, cs, ss, As, n_ent, theta, w = args
    L = law()
    K = MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    R = RN.sliding_iss_nd(L, K, f, theta0(), np.ones(C0.size), max_entries=n_ent, exit_kkt=True)
    return dict(q=q, t3_pred=t3_from_duals(K, f, R, theta, w), **kappa_from_path(K, f, R, theta, w))


def part_GT(G):
    """GT (referee round 68's τ₃ from the path, corrected: G's stored t3_pred used the duals at the last history time, before τ_d)."""
    return [gt_job((r['q'], np.array(r['cs']), np.array(r['ss']), np.array(r['As']), int(r['entries']), r['flow'][-1]['theta'],
                    r['flow'][-1]['w'])) for r in G if r['status'] == 0]


def pair_p(th, w, a, b):
    """(ii)(a)'s expansion parameter of the pair (a, b): p = |w_a − w_b||Δθ|/(|W|σ̄), |Δθ| = (Δc² + Δσ²)^{1/2} (referee round 74: a
    σ-dipole's Δc alone understates it), σ̄ the pair's mean |σ|; inf at W = 0."""
    th, w = np.atleast_2d(np.asarray(th, float)), np.asarray(w, float)
    W = w[a] + w[b]
    if W == 0:
        return float('inf')
    return float(abs(w[a] - w[b]) * np.hypot(th[a, 0] - th[b, 0], th[a, 1] - th[b, 1]) / (abs(W) * 0.5 * (abs(th[a, 1]) + abs(th[b, 1]))))


def fine_phases(sts, n0, td, cs, ss, As):
    """G on the fine grid (referee rounds 68, 74): from the states (τ, formed θ, formed w), the block's dipole phase and what follows — τ₃
    (the first τ past τ_d with more than n0 formed kernels) and τ₄ (the next formation); the dipole pair's v over the times in
    (τ_d, τ₃) with n0 formed (every tenth, for the slopes of v³) and at the last of them (v3 at t_v3); over the intermediate phase,
    the times in [τ₃, τ₄) with n3 formed, v at every tenth and the last (dip2: the slope of v³ there); the intermediate phase at
    (τ₃ + τ₄)/2 (formed, the pair's offset from the nearest feature (σ*) and (W − A)/|A|); the first landing past τ₃ (`landing`:
    centre of mass and W within LAND_TOL) — its τ, v and (ii)(a)'s p (`pair_p`) —; and over every state with a dipole pair and p < 1
    whether p ≤ P_LAND and the landing agree (agree), and over those past τ₄ (agree4)."""
    ns = [(tau, len(w)) for tau, th, w in sts]
    t3 = next((tau for tau, n in ns if tau > td and n > n0), None)
    n3 = next((n for tau, n in ns if tau == t3), None)
    t4 = next((tau for tau, n in ns if t3 is not None and tau > t3 and n > n3), None)
    dip = [(tau, dipole_pair(th, w)) for tau, th, w in sts if td < tau < (t3 if t3 is not None else np.inf) and len(w) == n0]
    dip = [(tau, d[2]) for tau, d in dip if d is not None]
    mid = [(tau, dipole_pair(th, w)) for tau, th, w in sts if t4 is not None and t3 <= tau < t4 and len(w) == n3]
    mid = [(tau, d[2]) for tau, d in mid if d is not None]
    out = dict(t3=t3, t4=t4, v3=dip[-1][1] if dip else None, t_v3=dip[-1][0] if dip else None, dip=[(a, b) for a, b in dip[::10]],
               dip2=[(a, b) for a, b in mid[::10]] + ([mid[-1]] if mid and (len(mid) - 1) % 10 else []), inter=None, t_land=None,
               v_land=None, p_land=None, agree=None, agree4=None)
    if t3 is not None and t4 is not None:
        mid = next(((th, w) for tau, th, w in sts if tau >= 0.5 * (t3 + t4)), None)
        lp = landing(mid[0], mid[1], cs, ss, As) if mid else None
        out['inter'] = None if lp is None else dict(n=len(mid[1]), off=lp['off'], dW=lp['dW'])
    agree, tot, agree4, tot4 = 0, 0, 0, 0
    for tau, th, w in sts:
        lp = landing(th, w, cs, ss, As)
        if lp is None:
            continue
        a, b = dipole_pair(th, w)[:2]
        p = pair_p(th, w, a, b)
        ok = lp['off'] <= LAND_TOL and abs(lp['dW']) <= LAND_TOL
        if p < 1:
            tot += 1
            agree += (p <= P_LAND) == ok
            if t4 is not None and tau > t4:
                tot4 += 1
                agree4 += (p <= P_LAND) == ok
        if ok and out['t_land'] is None and t3 is not None and tau > t3:
            out.update(t_land=tau, v_land=lp['v'], p_land=float(p))
    out['agree'], out['agree4'] = (agree, tot), (agree4, tot4)
    return out


def gl_job(args):
    """GL (referee rounds 68, 74): the loop on one dipole target at the mass EPS_GL to TAUS_GL — per output the dipole pair among the
    formed kernels: its distance, centre of mass's offset and (W − A)/|A| (`landing`), |Δc|, w_a, w_b, its mean width and the nearest
    feature's amplitude, the residual, the energy ½⟨r²⟩ of all kernels and of the formed ones —, with Λ(ε), Gm and μ₀ for the collapse
    law's constant (Δc⁻⁶ ∝ τ − τ_L, Δc ∝ τ^{−1/6}, the residual ∝ τ^{−2/3}, W − A ∝ τ^{−1/3} as Δc², E ∝ Δc⁸ with the centre of mass
    pinned)."""
    q, cs, ss, As = args
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), EPS_GL, np.r_[0.0, TAUS_GL], rtol=1e-9, atol=1e-12)
    rows = []
    for t, tau in enumerate(TAUS_GL, start=1):
        big = np.abs(W[t]) > 1e-2 * np.abs(As).min()
        th, w = C[t][big], W[t][big]
        lp = landing(th, w, cs, ss, As)
        r = f - K.synth(C[t], W[t])
        row = dict(tau=float(tau), d=None if lp is None else lp['d'], dW=None if lp is None else lp['dW'], rms=float(np.sqrt(np.mean(r ** 2))),
                   E=float(np.mean(r ** 2) / 2), E_formed=float(np.mean((f - K.synth(th, w)) ** 2) / 2))
        if lp is not None:
            a, b = dipole_pair(th, w)[:2]
            cm = (w[a] * th[a, 0] + w[b] * th[b, 0]) / (w[a] + w[b])
            k = int(np.argmin(np.abs(np.asarray(cs) - cm) / np.asarray(ss)))
            row.update(dc=float(abs(th[a, 0] - th[b, 0])), wa=float(w[a]), wb=float(w[b]), sbar=float(0.5 * (abs(th[a, 1]) + abs(th[b, 1]))),
                       A=float(As[k]))
        rows.append(row)
    return dict(q=q, rows=rows, Lam=float(Lam), Gm=float(L.G * L.m), mu0=float(4 * L.G * L.c / L.A), secs=time.time() - t0)


def part_GL(G):
    """GL: G's targets Q_GL."""
    return [gl_job((q, np.array(r['cs']), np.array(r['ss']), np.array(r['As']))) for q in Q_GL for r in G if r['q'] == q]


def gl_constant(x0, x1, Gm, mu0, Lam):
    """The collapse law's constant over GL's interval (x0, x1) (referee round 74): the measured d(Δc⁻⁶)/dτ against the proof's
    48aGm(W²/(w_a² + w_b²))Λ, a = E/Δc⁸, and against it over 1 + C — C = Gm W²(W − A)²(1/|w_a| + 1/|w_b|)/(μ₀Δc²(w_a² + w_b²)) the
    slaved coefficients' share of the separation's metric (W − A ∝ Δc² moves at the coefficients' mobility μ₀|w|, far below the centres'
    Gm), (Gm/(16μ₀))W²p²/(vσ²) to leading order —, each averaged over the interval's ends. Returns (p at x1, ratio, ratio with C, C at
    x1)."""
    rate = lambda x: 48 * x['E'] / x['dc'] ** 8 * Gm * (x['wa'] + x['wb']) ** 2 / (x['wa'] ** 2 + x['wb'] ** 2) * Lam
    cc = lambda x: (Gm * (x['wa'] + x['wb']) ** 2 * (x['dW'] * abs(x['A'])) ** 2 * (1 / abs(x['wa']) + 1 / abs(x['wb']))
                    / (mu0 * x['dc'] ** 2 * (x['wa'] ** 2 + x['wb'] ** 2)))
    meas = (x1['dc'] ** -6 - x0['dc'] ** -6) / (x1['tau'] - x0['tau'])
    p = abs(x1['wa'] - x1['wb']) * x1['dc'] / (abs(x1['wa'] + x1['wb']) * x1['sbar'])
    return (float(p), float(meas / (0.5 * (rate(x0) + rate(x1)))),
            float(meas / (0.5 * (rate(x0) / (1 + cc(x0)) + rate(x1) / (1 + cc(x1))))), float(cc(x1)))


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


def visibility(th):
    """Each kernel's largest φ on the samples (0 at σ = 0): below VIS it is invisible there (referee rounds 68, 74)."""
    th = np.atleast_2d(np.asarray(th, float)).reshape(-1, 2)
    return np.array([float(np.exp(-((XS - t[0]) ** 2).min() / (2 * t[1] ** 2))) if t[1] != 0 else 0.0 for t in th])


def path_stop(Rm, tol, max_events=40):
    """Why a path ended (referee round 74): 'tolerance' (the last block's rms below tol), 'budget' (its evaluations spent), 'dictionary'
    (every kernel active), 'events' (max_events blocks) or 'no entrant'."""
    if Rm['hist'] and Rm['hist'][-1][3] < tol:
        return 'tolerance'
    if Rm.get('out_of_budget'):
        return 'budget'
    if np.all(Rm['act']):
        return 'dictionary'
    return 'events' if len(Rm['hist']) >= max_events else 'no entrant'


def h_stats(blocks, m, As):
    """H's statistics of a path — blocks: each block's active (θ, w) — and of its end m (referee rounds 57, 61, 74), counting only
    kernels visible on the samples (`visibility` ≥ VIS): the most kernels of the other sign (below −1% of the smallest amplitude) in a
    block and the most negative coefficient; the end's formed kernels of the other sign (and all of them, the invisible included) and
    the smallest formed width; over all the end's active kernels those narrower than the samples' spacing and those of negative width
    (the widths are unconstrained), and the formed ones invisible on the samples; the largest |w|."""
    w, th = np.asarray(m['w'], float), np.atleast_2d(np.asarray(m['theta'], float)).reshape(-1, 2)
    amin = min(abs(a) for a in As)
    formed = np.abs(w) > 1e-2 * amin
    vis = visibility(th) >= VIS
    seen = [(np.asarray(wb, float), visibility(tb) >= VIS) for tb, wb in blocks]
    return dict(neg_path=int(max([int(((wb < -1e-2 * amin) & vb).sum()) for wb, vb in seen], default=0)),
                wmin_path=float(min([wb[vb].min() for wb, vb in seen if vb.any()], default=0.0)),
                neg_end=int((formed & vis & (w < 0)).sum()), neg_end_all=int((formed & (w < 0)).sum()),
                smin=float(np.abs(th[formed & vis, 1]).min()) if (formed & vis).any() else float('nan'),
                sub=int((np.abs(th[:, 1]) < XS[1] - XS[0]).sum()) if w.size else 0, negw=int((th[:, 1] < 0).sum()) if w.size else 0,
                sub_formed=int((formed & (np.abs(th[:, 1]) < XS[1] - XS[0])).sum()) if w.size else 0,
                negw_formed=int((formed & (th[:, 1] < 0)).sum()) if w.size else 0,
                inv_formed=int((formed & ~vis).sum()), wmax=float(np.abs(w).max()) if w.size else 0.0)


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
    pol = bool(args[5]) if len(args) > 5 else True
    t0 = time.time()
    f = feature_sum(cs, ss, As)
    if tol is None and pol:
        Rm, _ = run_paths(f, fixed=False)
    else:
        Rm = RN.sliding_iss_nd(law(), MW.BumpsCW(XS, S0), f, theta0(), np.ones(C0.size), max_events=40, exit_kkt=True,
                               max_rhs=MAX_RHS_MOVING, tol=1e-9 if tol is None else tol, polish=pol)
    m = match(Rm, cs, ss, As)
    blocks = [(h[1][h[2] != 0], h[2][h[2] != 0]) for h in Rm['hist'][1:]]
    out = dict(name=name, cs=list(cs), ss=list(ss), As=list(As), moving=m, outcome=outcome(m, As, cs, ss), tol=tol, polish=pol,
               stop=path_stop(Rm, 1e-9 if tol is None else tol), secs=time.time() - t0, **h_stats(blocks, m, As), **stop_anatomy(Rm))
    if tol is not None:
        out.update(last_tau=float(Rm['events'][-1][0]) if Rm['events'] else float('nan'),
                   rms_blocks=[float(h[3]) for h in Rm['hist'][-6:]], clusters=clusters(m, cs, ss, As))
    return out


def stop_anatomy(Rm):
    """Referee round 68: why a path ended and what its polishes did — the active kernels (all of them: the dictionary exhausted), the
    accepted polishes and those of them stopped by their budget (status 0), the active kernels of negative width and those invisible on
    the samples (largest φ below 10⁻¹⁰)."""
    act = np.flatnonzero(Rm['act'])
    th = np.asarray(Rm['c'], float)[act]
    pl = [p_ for b in Rm.get('polish_log', []) for p_ in b if p_.get('accepted')]
    vis = visibility(th) if act.size else np.zeros(0)
    return dict(n_act=int(act.size), exhausted=bool(act.size == C0.size), pol_acc=len(pl),
                pol_capped=sum(p_.get('status') == 0 for p_ in pl), neg_sigma=int(np.sum(th[:, 1] < 0)) if act.size else 0,
                invisible=int(np.sum(vis < VIS)))


def part_HP():
    """HP (referee round 68): H's path without the least-squares polish on the symmetric targets and two asymmetric ones."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, [tuple(x) + (None, False) for x in H_TARGETS if x[0] in HP_NAMES], chunksize=1)


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
                            rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2))), neg=int(np.sum(W[t][big] < 0)),
                            smin=float(np.abs(C[t][big][:, 1]).min()) if big.size else float('nan')))
    return dict(name=name, cs=list(cs), ss=list(ss), As=list(As), loops=out, secs=time.time() - t0)


def part_HL():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(hl_job, [x for x in H_TARGETS if x[0] in HL_NAMES], chunksize=1)


def fit_split(th, w, cs, ss, As):
    """The fit's split (referee rounds 74, 80): per feature its main kernel — the formed kernel of its sign nearest it in c, each kernel
    used once — and that kernel's largest relative error in (c, σ, w) (c over σ*); the other formed kernels: their number, Σ|w|, their
    sum's largest value on the samples over max f* (a null cluster: small, beside exact main kernels; a cancellation: comparable to the
    features, the main kernels off), their widths' range and their centres' largest |c|. th, w: the formed kernels."""
    th, w = np.atleast_2d(np.asarray(th, float)).reshape(-1, 2), np.asarray(w, float).reshape(-1)
    main, err = [], []
    for c_, s_, A_ in zip(cs, ss, As):
        cand = [i for i in range(w.size) if np.sign(w[i]) == np.sign(A_) and i not in main]
        if not cand:
            err.append(None)
            continue
        i = min(cand, key=lambda i: abs(th[i, 0] - c_))
        main.append(i)
        err.append(float(max(abs(th[i, 0] - c_) / s_, abs(abs(th[i, 1]) / s_ - 1), abs(w[i] / A_ - 1))))
    rest = [i for i in range(w.size) if i not in main]
    foot = (float(np.max(np.abs(MW.BumpsCW(XS, S0).synth(th[rest], w[rest]))) / np.max(np.abs(feature_sum(cs, ss, As))))
            if rest else 0.0)
    return dict(main_err=err, n_rest=len(rest), rest_w=float(np.abs(w[rest]).sum()), rest_foot=foot,
                rest_sig=[float(np.abs(th[rest, 1]).min()), float(np.abs(th[rest, 1]).max())] if rest else None,
                rest_c=float(np.abs(th[rest, 0]).max()) if rest else None)


def split_cells(sp):
    """`fit_split`'s two table cells: the main kernels' errors; the others' n, Σ|w|, sum on the samples, widths, largest |c|."""
    return (', '.join('—' if v is None else f"{v:.1e}" for v in sp['main_err']) + ' | ' + f"{sp['n_rest']}, {sp['rest_w']:.3f}, "
            f"{sp['rest_foot']:.1e}, " + ('—' if sp['rest_sig'] is None else f"{sp['rest_sig'][0]:.3f}–{sp['rest_sig'][1]:.3f}") + ', '
            + ('—' if sp['rest_c'] is None else f"{sp['rest_c']:.2f}"))


def end_split(r):
    """`fit_split` of a path's end (H, HP: r['moving'], its formed kernels)."""
    m = r['moving']
    w, th = np.asarray(m['w'], float), np.atleast_2d(np.asarray(m['theta'], float)).reshape(-1, 2)
    fm = np.abs(w) > 1e-2 * min(abs(a) for a in r['As'])
    return fit_split(th[fm], w[fm], r['cs'], r['ss'], r['As'])


def hll_job(args):
    """HLL (referee rounds 74, 80): one of HL's cells continued to TAUS_HLL[-1]: at TAUS_HLL the formed kernels' states, the largest
    |w|, the residual, and the fit's split — per feature its main kernel (the formed kernel of its sign nearest it in c) and that kernel's
    largest relative error in (c, σ, w) (c over σ*); the other formed kernels' Σ|w|, their sum's largest value on the samples over
    max f* (a null cluster: small against the residual), their widths and their centres' largest |c| (`fit_split`)."""
    name, eps = args
    _, cs, ss, As = next(x for x in H_TARGETS if x[0] == name)
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    _, W, C, _, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, TAUS_HLL], rtol=1e-9, atol=1e-12)
    out = []
    for t, tau in enumerate(TAUS_HLL, start=1):
        big = np.flatnonzero(np.abs(W[t]) > 1e-2 * min(abs(a) for a in As))
        th, w = C[t][big], W[t][big]
        out.append(dict(tau=float(tau), theta=th.tolist(), w=w.tolist(), wmax=float(np.abs(W[t]).max()),
                        rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2))), **fit_split(th, w, cs, ss, As)))
    return dict(name=name, eps=eps, loops=out, secs=time.time() - t0)


def part_HLL():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(hll_job, list(HLL), chunksize=1)


def loop_row(name, eps, x):
    """HL's and HLL's row (referee round 74): the formed kernels, the visible ones (`visibility` ≥ VIS) of the other sign (the most
    negative w), the invisible formed ones, the smallest visible width, those centred outside the samples' window, the largest |w|, the
    residual."""
    w, th = np.asarray(x['w'], float), np.atleast_2d(np.asarray(x['theta'], float)).reshape(-1, 2)
    vis = visibility(th) >= VIS if w.size else np.zeros(0, bool)
    neg = (w < 0) & vis
    return (f"| {name} | {eps:g} | {x['tau']:g} | {w.size} | {int(neg.sum())}" + (f" ({w[neg].min():+.3f})" if neg.any() else '')
            + f" | {int((~vis).sum())} | " + (f"{np.abs(th[vis, 1]).min():.3f}" if vis.any() else '—')
            + f" | {int((np.abs(th[:, 0]) > XS[-1]).sum())} | {x['wmax']:.3g} | {x['rms']:.1e} |")


def report(A=None, B=None, C=None, D=None, Ee=None, F=None, G=None, H=None, H2=None, H3=None, HL=None, HP=None, GL=None, GT=None,
           A2=None, HLL=None):
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
    L_ += a2_lines(A2)
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
               f"{N_C} targets per minimal separation (seed {SEED_C}; widths log-uniform in [{WID_C[0]:g}, {WID_C[1]:g}], amplitudes ±U(0.4, 1); "
               "a moving-width path stopped by its budget counted at its stalled block — G — not at the block before it, referee round 68):", '',
               '| Δ ≥ | moving widths: exact | formed kernels (median [range]) | residual (median [max]) | fixed widths: kernels (median '
               '[range]) | residual (median [max]) |', '|---|---|---|---|---|---|']
        for d in DELTAS_C:
            rs = [r for r in C if r['dmin'] == d]
            km, kf = [moving_state(r)[0] for r in rs], [r['fixed']['n'] for r in rs]
            em, ef = [moving_state(r)[1] for r in rs], [r['fixed']['rms'] for r in rs]
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
        L_ += g_fine_lines(G, GT)
    if GL is not None:
        L_ += gl_lines(GL)
    if H is not None:
        L_ += ['## H. Two features of one sign, close', '',
               "Two features φ(·; c, σ*) of the given amplitudes (both positive) at the centres ±0.08 — 0.16 apart, 1.33 times the kernels' "
               f"starting width σ₀ = {S0:g} — unless stated, the moving-width path (as C's: its tolerance 1e-09, its budget "
               f"{MAX_RHS_MOVING:.0e}): the outcome (as C's), the most kernels of the other sign in a block along the path and the most "
               "negative coefficient, the end's formed kernels of the other sign and the smallest formed width — counting only kernels visible on "
               f"the samples (largest φ there at least {VIS:g}; referee round 74), the end's formed ones of the other sign including the "
               "invisible in brackets —, over all the end's active kernels those narrower than the samples' spacing and those of negative "
               "width (the widths are not constrained) and how many of them are formed, the formed kernels invisible on the samples, the "
               "largest |w|, the residual:", '']
        L_ += h_table(H)
    if H2 is not None:
        L_ += ['### H2. The σ* = 0.06 paths to a tolerance of ' + f"{TOL_H2:g}", '',
               "The same paths with the tolerance lowered (the 1e-09 ends were the tolerance's): the outcome, the entries and the last one's "
               "τ, the residuals of the last blocks, the formed kernels by the feature nearest each (their count and mass over the feature's, "
               "Σw_i|σ_i|/(Aσ*)), and "
               "H's columns:", '',
               '| target | outcome | entries (last τ) | the last blocks\' residuals | formed per feature (count: mass) | visible negative along '
               'the path (most negative w) | at the end (all formed) | smallest visible width | sub-sample, negative widths (all active; formed) | '
               'max\\|w\\| |', '|---|---|---|---|---|---|---|---|---|---|']
        for r in H2:
            L_.append(f"| {r['name']} | {r['outcome']} | {r['moving']['entries']} ({r['last_tau']:.2e}) | "
                      + ', '.join(f"{v:.2e}" for v in r['rms_blocks']) + ' | ' + ', '.join(f"{n}: {s_:.4f}" for n, s_ in r['clusters'])
                      + f" | {r['neg_path']} ({r['wmin_path']:+.3f}) | {r['neg_end']} ({r.get('neg_end_all', '—')}) | {r['smin']:.2e} | "
                      f"{r['sub']}, {r['negw']} ({r.get('sub_formed', '—')}, {r.get('negw_formed', '—')}) | {r['wmax']:.3g} |")
        if all('n_act' in r for r in H2):                               # referee round 68: why the paths ended
            L_ += ['', "Why they ended (referee rounds 68, 74): " + '; '.join(
                f"{r['name']}: the stop {r.get('stop', '—')}, {r['n_act']} of {C0.size} kernels active" + (' (the dictionary exhausted)' if r['exhausted'] else '')
                + f", {r['pol_capped']} of the {r['pol_acc']} accepted polishes stopped by their budget (status 0), {r['neg_sigma']} "
                f"active kernels of negative width, {r['invisible']} invisible on the samples (largest φ below 10⁻¹⁰)" for r in H2) + '.']
        L_.append('')
    if HP is not None:
        L_ += ['### HP. H\'s path without the polish (referee round 68)', '',
               "The moving-width path (tolerance 10⁻⁹) with the least-squares polish off, against H's (with it): the outcome, the "
               "negative kernels visible on the samples along the path and at the end (all formed at the end in brackets), the largest |w|, "
               "the residual, the active kernels and why the path stopped (referee round 74: tolerance, budget, dictionary exhausted), and the "
               "fit's split at the end (referee round 80; as HLL's):", '',
               '| target | outcome (H) | visible negative along the path, at the end (H) | max\\|w\\| (H) | residual (H) | active | stop | '
               'the main kernels\' errors | the others: n, Σ\\|w\\|, on the samples, widths, \\|c\\| max |',
               '|---|---|---|---|---|---|---|---|---|']
        Hn = {r['name']: r for r in (H if H is not None else [])}
        for r in HP:
            h = Hn.get(r['name'], {})
            L_.append(f"| {r['name']} | {r['outcome']} ({h.get('outcome', '—')}) | {r['neg_path']}, {r['neg_end']} [{r.get('neg_end_all', '—')}] "
                      f"({h.get('neg_path', '—')}, {h.get('neg_end', '—')} [{h.get('neg_end_all', '—')}]) | {r['wmax']:.3g} ("
                      + (f"{h['wmax']:.3g}" if 'wmax' in h else '—') + f") | {r['moving']['rms']:.1e} ("
                      + (f"{h['moving']['rms']:.1e}" if 'moving' in h else '—') + f") | {r['n_act']} | {r.get('stop', '—')} | "
                      f"{split_cells(end_split(r))} |")
        L_.append('')
    if H3 is not None:
        L_ += ['### H3. The separation scan', '',
               "σ* = 0.06, equal amplitudes 1, at ±c (2c in units of σ₀ in brackets; H's ±0.08 is 1.33σ₀), H's columns:", '']
        L_ += h_table(sorted(list(H3) + [r for r in (H if H is not None else []) if r['name'] == 'σ* 0.06, equal'],
                             key=lambda r: r['cs'][1]), sep=True)
    if HL is not None:
        L_ += ['### HL. The real loop on H\'s targets (referee round 68: the asymmetric ones, to τ = 10⁵)', '',
               "The free law from α = 0 at the masses ε (`loop_nd`, rtol 10⁻⁹) at the times τ: the formed kernels (|w| above 1% of the "
               f"smallest amplitude), those visible on the samples (largest φ there at least {VIS:g}; referee round 74) of the other sign, the "
               "invisible formed ones, the smallest visible width, the formed kernels centred outside the samples' window, the largest |w|, "
               "the residual, and the fit's split (referee round 80; as HLL's):", '',
               '| target | ε | τ | formed | visible of the other sign (most negative w) | invisible | smallest visible width | outside the window | '
               'max\\|w\\| | residual | the main kernels\' errors | the others: n, Σ\\|w\\|, on the samples, widths, \\|c\\| max |',
               '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for r in HL:
            L_ += [loop_row(r['name'], x['eps'], x) + f" {split_cells(fit_split(x['theta'], x['w'], r['cs'], r['ss'], r['As']))} |"
                   for x in r['loops']]
        L_.append('')
    if HLL is not None and len(HLL):
        L_ += ['### HLL. HL\'s cells to τ = 10⁹ (referee rounds 74, 80)', '', "HL's columns, and the fit's split: per feature its main kernel's "
               "largest relative error in (c, σ, w); the other formed kernels — their number, Σ\\|w\\|, their sum's largest value on the "
               "samples over max f* (a null cluster where it is small against the residual), their widths and their centres' largest \\|c\\|:", '',
               '| target | ε | τ | formed | visible of the other sign (most negative w) | invisible | smallest visible width | outside the window | '
               'max\\|w\\| | residual | the main kernels\' errors | the others: n, Σ\\|w\\|, on the samples, widths, \\|c\\| max |',
               '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for r in HLL:
            for x in r['loops']:
                L_.append(loop_row(r['name'], r['eps'], x) + f" {split_cells(x) if 'n_rest' in x else '— | —'} |")
        L_.append('')
    L_ += ['Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def v3_slopes(dp, rate):
    """The slopes of v³ between successive (τ, v) over rate."""
    return [(b[1] ** 3 - a[1] ** 3) / (b[0] - a[0]) / rate for a, b in zip(dp[:-1], dp[1:])]


def kappa_loop(ph, td, t3p=None):
    """The loop's κ (referee round 80): v_land over v(τ₃⁻), the loop's v at its last fine-grid time t_v3 before τ₃ carried by the law
    v³ ∝ τ − τ_d to the path's τ₃ t3p — v(t_v3) itself lags v(τ₃⁻) by the grid's offset ((τ₃ − τ_d)/(t_v3 − τ_d))^{1/3} —; without
    t3p, v_land/v(t_v3). Returns (κ, the offset) or (None, None)."""
    if not ph.get('v_land') or not ph.get('v3'):
        return None, None
    tv = ph.get('t_v3')
    off = ((t3p - td) / (tv - td)) ** (1 / 3) if t3p is not None and tv is not None and tv > td and t3p > td else 1.0
    return ph['v_land'] / (ph['v3'] * off), off


def g_fine_lines(G, GT=None):
    """G's phases on the fine grid (referee rounds 68, 74): per target and mass, τ_d, τ₃, τ₄, the dipole phase's slopes of v³ over the
    law (median), v before τ₃ over the law's at that time, the intermediate phase (formed, offset, (W − A)/|A|; its slopes of v³ over
    3μ₀K′Λ, K′ the path's), the first landing (τ, κ = v_land/v(τ₃⁻) by `kappa_loop`, p); per target the path's predictions
    (`kappa_from_path`, GT) against the loop's; κ against the mass; over every state whether p ≤ P_LAND and the landing agree."""
    if not any(r.get('fine') for r in G):
        return []
    gt = {x['q']: x for x in (GT if GT is not None else [])}
    out = ['', "The phases on a fine grid (τ = 10 … 140 by 0.1, with G's times; referee rounds 68, 74): τ₃ the first formation past the "
           "block's entry τ_d, τ₄ the next; over (τ_d, τ₃) the slopes of v³ over 3μ₀KΛ(ε) (median); v at the last time before τ₃ over the "
           "law's (3μ₀KΛ(ε)(t − τ_d))^{1/3} at that time t; at (τ₃ + τ₄)/2 the formed kernels, the dipole pair's offset from the nearest "
           "feature (σ*) and (W − A)/|A|, and over [τ₃, τ₄) the slopes of v³ over 3μ₀K′Λ(ε), K′ the three-kernel block's from the path "
           f"(below; median); the first landing past τ₃ (centre of mass and W within {LAND_TOL:g}): its τ, κ = v_land/v(τ₃⁻) with "
           "v(τ₃⁻) the loop's v at the last time before τ₃, κ′ with that v carried by the law to the path's τ₃ (`kappa_loop`; referee "
           "round 80: the grid's offset removed, below), and p = |w_a − w_b||Δθ|/(|W|σ̄); † another sequence — the intermediate slope "
           f"beyond 1 ± {SEQ_TOL:g} times the path's law —, its κ not the path's:", '',
           '| target | ε | τ_d | τ₃ | τ₄ | dipole: slope | v(τ₃⁻)/law | intermediate: formed, offset, (W − A)/\\|A\\|, slope | '
           "landing: τ, κ, κ′, p |", '|---|---|---|---|---|---|---|---|---|']
    agree, agree4, agree7, kap, offs = [0, 0], [0, 0], [0, 0], {}, {}
    for r in G:
        lw = r.get('law') or []
        if not r.get('fine') or not lw or lw[-1] is None:
            continue
        Kc, mu0 = lw[-1]['v'] ** 2 * (lw[-1]['E'] - lw[-1]['E_inf']), lw[-1]['mu0']
        Kp = gt.get(r['q'], {}).get('Kp')
        Lam = {x['eps']: x['Lam'] for x in r['loops']}
        for eps, ph in sorted(r['fine'].items(), key=lambda t: -t[0]):
            for acc, a_ in ((agree, ph.get('agree')), (agree4, ph.get('agree4')), (agree7, ph.get('agree') if eps <= 1e-7 else None)):
                if a_:
                    acc[0] += a_[0]
                    acc[1] += a_[1]
            lam = Lam.get(eps, RT.transit_scale(law(), eps))              # `loop_nd`'s Λ(ε) (the masses without G's tables)
            rate = 3 * mu0 * Kc * lam
            sl = v3_slopes(ph['dip'], rate)
            tv = ph.get('t_v3', ph['t3'])
            vlaw = (rate * (tv - r['tau_d'])) ** (1 / 3) if tv is not None and tv > r['tau_d'] else None
            sl2 = v3_slopes(ph.get('dip2') or [], 3 * mu0 * Kp * lam) if Kp else []
            k_, of_ = kappa_loop(ph, r['tau_d'], gt.get(r['q'], {}).get('t3_path'))
            k0 = kappa_loop(ph, r['tau_d'])[0]
            seq = bool(sl2) and abs(np.median(sl2) - 1) > SEQ_TOL         # referee round 80: the loop on another sequence
            kap.setdefault(r['q'], {})[eps] = (k_, seq)
            offs.setdefault(r['q'], {})[eps] = of_
            it = ph['inter']
            out.append(f"| {r['q']} | {eps:g} | {r['tau_d']:.4g} | " + ('—' if ph['t3'] is None else f"{ph['t3']:g}") + ' | '
                       + ('—' if ph['t4'] is None else f"{ph['t4']:g}") + ' | ' + (f"{np.median(sl):.4f}" if sl else '—') + ' | '
                       + (f"{ph['v3'] / vlaw:.4f}" if ph['v3'] and vlaw else '—') + ' | '
                       + ('—' if it is None else f"{it['n']}, {it['off']:.3g}, {it['dW']:+.3f}, " + (f"{np.median(sl2):.4f}" if sl2 else '—'))
                       + ' | ' + ('—' if ph['t_land'] is None else f"{ph['t_land']:g}, " + ('—' if k0 is None else f"{k0:.3f}") + ', '
                                  + ('—' if k_ is None else f"{k_:.3f}") + ('†' if seq else '') + f", {ph['p_land']:.3f}") + ' |')
    if gt:
        out += ['', "The intermediate phase from the path (`kappa_from_path`, referee rounds 74, 80): K from the third order, −⟨r_dip, "
                "∇³φ:p³⟩/24 at the block's optimal dipole fit, beside G's v²(E − E_inf); τ₃ from the duals; K′/K, K′ the three-kernel block's "
                f"(the τ₃ entrant added), with that block's closure rms (below {CLOSURE_EXACT:g} the block fits exactly: K′ = 0, no τ₄, "
                "κ = 1); τ₄ from the duals moved on by that block's residual; for K′ < 0 the dissolution time τ₃ + K(τ₃ − τ_d)/|K′|, at which "
                "v³ would reach 0 (the four stages need τ₄ before it); κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3}; the landing's v by the law "
                "v³ = 3μ₀Λ(ε)(K(τ₃ − τ_d) + K′(τ₄ − τ₃)) with the path's τ₃ and τ₄ — each against the loop's at the smallest mass (its κ as in "
                "the table above) —; and v at τ = 6400 over v_land at ε = 10⁻¹⁰, minus 1 (v after τ₄):", '',
                '| target | K: third order (G\'s) | τ₃: path (loop) | K′/K (the closure\'s rms) | τ₄: path (loop) | dissolution | '
                'κ: path (loop) | v_land: law (loop) | v(6400)/v_land − 1 |', '|---|---|---|---|---|---|---|---|---|']
        for r in G:
            g_, lw = gt.get(r['q']), r.get('law') or []
            if not g_ or not r.get('fine') or not lw or lw[-1] is None or g_.get('K_third') is None:
                continue
            Kc, mu0, e0 = lw[-1]['v'] ** 2 * (lw[-1]['E'] - lw[-1]['E_inf']), lw[-1]['mu0'], min(r['fine'])
            ph = r['fine'][e0]
            Lam = {x['eps']: x['Lam'] for x in r['loops']}
            lam0 = Lam.get(e0, RT.transit_scale(law(), e0))
            t3p, t4p, K3, Kp, ex = g_['t3_path'], g_.get('t4_path'), g_['K_third'], g_.get('Kp'), bool(g_.get('exact'))
            br = K3 * (t3p - r['tau_d']) + (0.0 if ex else Kp * (t4p - t3p) if t4p is not None and Kp is not None else np.nan)
            vl = (3 * mu0 * lam0 * br) ** (1 / 3) if np.isfinite(br) and br > 0 else None
            p10 = r['fine'].get(1e-10)
            x64 = next((x for x in r['loops'] if x['eps'] == 1e-10 and x['tau'] == 6400.0), None)
            dp64 = dipole_pair(x64['theta'], x64['w']) if x64 is not None else None
            fr = dp64[2] / p10['v_land'] - 1 if dp64 is not None and p10 and p10['v_land'] else None
            k_loop = kappa_loop(ph, r['tau_d'], t3p)[0]
            rm = '' if g_.get('rms2') is None else f" ({g_['rms2']:.1e})"
            kpc = f"0, exact{rm}" if ex else ('—' if Kp is None else f"{Kp / K3:+.4f}{rm}")
            t4c = '— (exact)' if ex else ('—' if t4p is None else f"{t4p:.4g}")
            out.append(f"| {r['q']} | {K3:.5e} ({Kc:.5e}) | {t3p:.4g} (" + ('—' if ph['t3'] is None else f"{ph['t3']:g}") + f") | {kpc} | "
                       + t4c + ' (' + ('—' if ph['t4'] is None else f"{ph['t4']:g}") + ') | '
                       + ('—' if g_.get('t_diss') is None else f"{g_['t_diss']:.4g}") + ' | '
                       + ('—' if g_.get('kappa_path') is None else f"{g_['kappa_path']:.4f}")
                       + ' (' + ('—' if k_loop is None else f"{k_loop:.4f}") + ') | ' + ('—' if vl is None else f"{vl:.4g}") + ' ('
                       + ('—' if not ph['v_land'] else f"{ph['v_land']:.4g}") + ') | ' + ('—' if fr is None else f"{fr:+.1e}") + ' |')
    out += ['', "κ′ against the mass (ε = " + ', '.join(f"{e:g}" for e in sorted({e for d in kap.values() for e in d}, reverse=True))
            + "; † another sequence): " + '; '.join(f"{q}: " + ', '.join('—' if d[e][0] is None else f"{d[e][0]:.3f}" + ('†' if d[e][1] else '')
                                                                         for e in sorted(d, reverse=True)) for q, d in kap.items())
            + ". The grid's offset ((τ₃ − τ_d)/(t_v3 − τ_d))^{1/3} − 1 that κ′ removes, at ε = " + f"{min(e for d in offs.values() for e in d):g}: "
            + ', '.join(f"{q}: " + ('—' if d[min(d)] is None else f"{d[min(d)] - 1:+.2%}") for q, d in offs.items())
            + f". Over every fine-grid state with a dipole pair and p < 1, p ≤ {P_LAND:g} and the landing agree in {agree[0]} of {agree[1]} "
            f"({agree7[0]} of {agree7[1]} at ε ≤ 10⁻⁷; past τ₄ {agree4[0]} of {agree4[1]}).", '']
    return out


def gl_lines(GL):
    """GL's collapse (referee rounds 68, 74): per target the local exponents in τ of the dipole pair's distance, (W − A) and the residual,
    and the collapse law's constant over successive rows (`gl_constant`)."""
    out = []
    for g in GL:
        rs = [x for x in g['rows'] if x['d'] is not None and x['d'] > 0]
        if len(rs) < 3:
            continue
        ex = lambda k, a, b: float(np.log(abs(b[k]) / abs(a[k])) / np.log(b['tau'] / a['tau'])) if a[k] and b[k] else float('nan')
        pts = rs[::4] + ([rs[-1]] if (len(rs) - 1) % 4 else [])
        out += ['', f"### GL. The collapse (referee rounds 68, 74): target {g['q']} at ε = {EPS_GL:g} to τ = {TAUS_GL[-1]:g}", '',
                "The dipole pair's distance, (W − A)/|A| and the residual, and their local exponents in τ against the law's −1/6, −1/3, "
                "−2/3 (E ∝ Δc⁸ with the centre of mass pinned):", '',
                '| τ | distance | (W − A)/\\|A\\| | residual | exponents since the previous row: distance, W − A, residual |',
                '|---|---|---|---|---|']
        prev = None
        for x in pts:
            out.append(f"| {x['tau']:.4g} | {x['d']:.3e} | {x['dW']:+.3e} | {x['rms']:.3e} | "
                       + ('—' if prev is None else f"{ex('d', prev, x):+.4f}, {ex('dW', prev, x):+.4f}, {ex('rms', prev, x):+.4f}") + ' |')
            prev = x
        if 'Gm' in g and all('dc' in x for x in rs):
            cons = [gl_constant(a, b, g['Gm'], g['mu0'], g['Lam']) for a, b in zip(rs[:-1], rs[1:])]
            ef = [x['E_formed'] / x['E'] for x in rs]
            sel = list(range(0, len(cons), 4)) + ([len(cons) - 1] if (len(cons) - 1) % 4 else [])
            out += ['', "The law's constant (referee round 74) over successive rows: the measured d(Δc⁻⁶)/dτ over 48aGm(W²/(w_a² + "
                    "w_b²))Λ(ε), a = E/Δc⁸ (the proof's metric: the centres alone), and over it divided by 1 + C, C the slaved coefficients' "
                    "share of the separation's metric (`gl_constant`); p = |w_a − w_b||Δc|/(|W|σ̄) and C at the interval's end:", '',
                    '| τ | p | ratio | with C | C |', '|---|---|---|---|---|']
            out += [f"| {rs[k + 1]['tau']:.4g} | {cons[k][0]:.4f} | {cons[k][1]:.4f} | {cons[k][2]:.4f} | {cons[k][3]:.4f} |" for k in sel]
            out += ['', f"Over the {len(cons)} intervals the ratio is {min(c[1] for c in cons):.4f}–{max(c[1] for c in cons):.4f}, with C "
                    f"{min(c[2] for c in cons):.4f}–{max(c[2] for c in cons):.4f}; the formed kernels' energy over all kernels' "
                    f"{min(ef):.6f}–{max(ef):.6f}."]
            pr = np.array([abs(x['wa'] - x['wb']) * x['dc'] / (abs(x['wa'] + x['wb']) * x['sbar']) for x in rs])
            av = np.array([x['E'] / x['dc'] ** 8 for x in rs])
            if pr.max() - pr.min() > 0.05:                              # referee round 74: what the constant with C leaves
                sl = [np.log(b['E_formed'] / a['E_formed']) / np.log(b['dc'] / a['dc']) for a, b in zip(rs[:-1], rs[1:])]
                ks = [(x - 8) / (0.5 * (pa + pb)) for x, pa, pb in zip(sl, pr[:-1], pr[1:])]   # referee round 80: the local slope
                rem = [c[2] / (x / 8) - 1 for c, x in zip(cons, sl)]
                out += [f"The remainder with C against p: (ratio with C − 1)/p {min((c[2] - 1) / c[0] for c in cons):.3f}–"
                        f"{max((c[2] - 1) / c[0] for c in cons):.3f}; a = E/Δc⁸ itself goes {av[0]:.4g} → {av[-1]:.4g} along the run. "
                        "E's next order in Δc from the local slope (referee round 80; the formed kernels' energy, the idle ones' floor "
                        f"left out): (d ln E/d ln Δc − 8)/p over the intervals (p their mean) {min(ks):.3f}–{max(ks):.3f}, median "
                        f"{np.median(ks):.3f}, which raises the force over 8aΔc⁷ by 1 + {np.median(ks) / 8:.3f}p; with each interval's own "
                        f"force factor (d ln E/d ln Δc)/8 the ratio with C is 1 {min(rem):+.1e} … {max(rem):+.1e} (the largest |·| at "
                        f"τ = {rs[1 + int(np.argmax(np.abs(rem)))]['tau']:.4g}); the mean width's share of the metric, p²/36, at most "
                        f"{pr.max() ** 2 / 36:.5f}."]
        out.append('')
    return out


def h_table(H, sep=False):
    """H's table (H and H3)."""
    out = [('| separation | ' if sep else '| target | ') + 'outcome | visible negative kernels along the path (most negative w) | at the end '
           '(all formed) | smallest visible width | sub-sample, negative widths (all active; formed) | formed invisible | max\\|w\\| | residual | '
           'the main kernels\' errors | the others: n, Σ\\|w\\|, on the samples, widths, \\|c\\| max |',
           '|---|---|---|---|---|---|---|---|---|---|---|']
    for r in H:
        lab = f"±{r['cs'][1]:g} ({(r['cs'][1] - r['cs'][0]) / S0:.2f}σ₀)" if sep else r['name']
        out.append(f"| {lab} | {r['outcome']} | {r['neg_path']} ({r['wmin_path']:+.3f}) | {r['neg_end']} ({r.get('neg_end_all', '—')}) | "
                   f"{r['smin']:.2e} | {r.get('sub', '—')}, {r.get('negw', '—')} ({r.get('sub_formed', '—')}, {r.get('negw_formed', '—')}) | "
                   f"{r.get('inv_formed', '—')} | {r['wmax']:.3g} | {r['moving']['rms']:.1e} | {split_cells(end_split(r))} |")
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
        L_ += ['', f"(ii)(a)'s width split on the landed pairs (L, {len(spl)} rows, {len({round(float(x), 6) for x in rt})} distinct to six "
               "decimals: frozen states repeat): the observed σ_b − σ_a over (w_b − w_a)Δc²/(6Wσ) "
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
    the collapse's exponents (distance, residual) and those not yet collapsing by τ = 6400 (|exponent| below 10⁻³: landed at a small p,
    they follow the same law later — referee round 68, GL)."""
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
            f"phase), the landed v over the law's at the grid's τ₃ {span(vr, '.2f')} (on the fine grid κ, below), v's ratio over the last "
            f"decade {span(dec, '.4f')}; on the {len(both)} cells landed at both 1600 and 6400, {len(fr)} not yet collapsing by τ = 6400 "
            "(|exponent| < 10⁻³: landed at a small p, the same law later — GL) (" + ', '.join(f"{k[1]} at {k[2]:g}" for k in fr)
            + ")" + (f", the others collapsing with the exponents {col[:, 0].min():+.3f} … {col[:, 0].max():+.3f} (distance) and "
                     f"{col[:, 1].min():+.3f} … {col[:, 1].max():+.3f} "
                     "(residual)" if len(col) else '') + '.']


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'moving_width_raw.npz')
    md = os.path.join(OUT, 'moving_width.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'H2', 'H3', 'HL', 'HP', 'GL', 'GT', 'A2', 'HLL')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D), ('E', part_E), ('F', part_F),
                  ('G', lambda: part_G(res['C'])), ('H', part_H), ('H2', part_H2), ('H3', part_H3), ('HL', part_HL), ('HP', part_HP),
                  ('GL', lambda: part_GL(res['G'])), ('GT', lambda: part_GT(res['G'])), ('A2', part_A2), ('HLL', part_HLL)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
