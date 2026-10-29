"""
Proposition D‴.58: the transit under width mismatch — a feature wider than the kernels is formed as a sequence of Gaussian quadratures,
the error falling by t = (σ* − σ)/(σ* + σ) per kernel
====================================================================================================================================

Kernels of width σ = 0.12 (the program's), a feature f* = Aφ_{σ*}(· − c*) with σ* = σ(1 + (τ/σ)²)^{1/2} (`src/core/wide_feature.py`):

A. The rules on ℝ: for τ/σ in TAUS_A and K ≤ K_A, the half-squeeze rule's L² error (trapezoid on a fine grid) against the exact Mehler
   series and its first term; the Gauss–Hermite rule of the mixing law (moment matching); the half-squeeze nodes with least-squares
   amplitudes; the optimal K-kernel fit (variable projection from the half-squeeze rule) — the ratios per kernel against 1/t and 1/λ;
   the optimum's nodes against the half-squeeze nodes. A2: the series' tail — the ratios D_{2K+2j,K}/D_{2K,K} for K ≤ K_A2 (integers;
   their j-th differences in K against j!·binom(1/2, j)(−2)^j) and the tail's share against ₂F₁(−½, −½; 1; t²) − 1. A3: the scaled
   Gauss–Hermite rules (`wide_feature.scaled_rule`) at the natural Fock and the Mercer-zero scalings against the half-squeeze one (their
   ratios per kernel against (2 − λ)/λ and 1/q), and the variance of least error (`best_scale`): one parameter's share of the optimum's
   gain over the rule (referee round 48). A4: (v), the frame continued to features narrower than the kernels (σ*/σ in NARROW_A): the
   Gauss–Hermite cosine sum's error (`narrow_rule`) against the series at |t|.
B. The transit: the sliding path (polished; the KKT exit rule) on f* (A_B, c* = C_B off the kernels' grid) sampled on X_B, in six
   numerical variants (VARIANTS_B: the polish step, the starting grid's offset, the polish's budget of evaluations, its repetition
   after the stationarity test and its gradient test; referee rounds 48, 54 and 60): after each entry the relaxed block's relative error against the optimal fit
   with as many kernels on the same samples (variable projection from the block and from the half-squeeze rule, the better), whether
   the relaxation passed the stationarity test, the first non-optimal block's polishes (status, evaluations, accepted), the block's
   nodes against the optimum's, the continuous half-squeeze error; the entries (the entrant's rate against the largest, its dual), the
   chunks redone after a stall or a failed solve, the polish's rejections. B2: variable projection from every non-optimal block. B3: each
   path's first non-optimal block, its first polish restarted from its recorded start one setting at a time (the budget, the gradient
   test). B4: the floor — the same Levenberg–Marquardt solve with the residual in double and in long double (referee round 70).
C. Under noise: y = f* + ξ (i.i.d. N(0, s²) on X_B), the path to K_C + 1 entries (its centres in the window ± BOX_C·σ, the thresholds
   on max|y|, the KKT exit rule): the discrepancy rule's stop (the first block with rms ≤ (1 + DELTA)s) against its prediction from (i)
   — the first K with e_K²‖f*‖²_rms ≤ s²((1 + δ)² − 1 + 2K/N), e_K the half-squeeze rule's L² error — and the probability of each stop
   from the rms's expectation and scatter (`stop_law`), the stopped fit's error against (e_K²‖f*‖²_rms + 2Ks²/N)^{1/2} (e_K B's optima
   on the samples); the kernels at the feature and elsewhere; the largest relative gap's K.

Usage: python experiments/exp_width_mismatch.py [--only A,B,C,D,B2,A5,B3,B4] [--variants i,j] [--from-raw] → results/chaos/width_mismatch.md
(+ _raw.npz). --variants: B for those variants only (indices into VARIANTS_B), merged into the raw's B. WM_WORKERS (default 2).
"""
import os
import sys
import time
import warnings
from fractions import Fraction
from math import comb, factorial

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
from src.core import wide_feature as W
import exp_cohesion_leak as E
import exp_transit_mass as TM

OUT = os.environ.get('WM_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('WM_WORKERS', '2'))
SIG = 0.12
TAUS_A, K_A, K_A2 = (0.25, 0.5, 1.0, 2.0), 10, 30
NARROW_A = (0.5, 0.8)                                                 # A4: σ*/σ < 1
E_FLOOR, G_FLOOR = 1e-10, 1e-13                                      # A: the optimum not computed below (its solve stops short), the grid's floor
TAUS_B, A_B, C_B, K_B = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0), 0.8, 0.03, 12
X_B = np.linspace(-1.6, 1.6, 321)
LEVELS_C, N_DRAW_C, TAU_C, K_C, DELTA = (1e-4, 1e-3, 1e-2), 10, 1.0, 10, 0.1
KS_D, J_MAX_D = (10, 20, 40, 80, 160, 320), 7                        # D: the tail ratio's K and the ratios' degree checked
MAX_DOUBLINGS = 24                                                    # chunks to 10⁵·2²⁴ ≈ 1.7·10¹² (beyond ~10¹¹ LSODA fails from relaxed
                                                                      # states and `relax_chunk` redoes the chunk by Radau: round 48)
VARIANTS_B = (('base', 1.0, 0.0, 5000, 0, 1e-15), ('polish 2σ', 2.0, 0.0, 5000, 0, 1e-15), ('grid +0.05', 1.0, 0.05, 5000, 0, 1e-15),
              ('polish 2σ, 5·10⁴ evaluations', 2.0, 0.0, 50_000, 0, 1e-15), ('polish 2σ, 5·10⁴, repeated', 2.0, 0.0, 50_000, 3, 1e-15),
              ('polish 2σ, 5·10⁴, repeated, no gradient test', 2.0, 0.0, 50_000, 3, None))
# B: (name, the polish step, the starting grid's shift, the polish's budget of evaluations, its repetitions after the stationarity test,
# its absolute gradient test — None: none; the fourth and fifth variants: referee round 54, the sixth: round 60)
OPT_TOL = 1e-3                                                        # B: a block is the optimal fit if its error is within 0.1% of it
UNCONV = 1e-9                                                         # A3: an optimum the one-parameter rule beats by more is unconverged
BOX_C = 2.0                                                           # C: the noisy path's centres in the samples' window ± BOX_C·σ


def sigma_star(r_):
    return SIG * np.sqrt(1 + r_ ** 2)


def part_A():
    rows = []
    for r_ in TAUS_A:
        ss = sigma_star(r_)
        P = W.params(SIG, ss)
        for K in range(1, K_A + 1):
            y, a = W.half_squeeze(K, SIG, ss)
            e = W.l2_error(y, a, SIG, ss)
            s2, f2 = W.mehler_series(K, P['t'])
            if np.sqrt(s2) < G_FLOOR:                                    # the trapezoid rule's double-precision floor
                break
            ym, am = W.moment_matching(K, SIG, ss)
            V = W.scalings(SIG, ss)
            vb, eb = W.best_scale(K, SIG, ss)
            row = dict(r=r_, K=K, t=P['t'], lam=P['lam'], e_hs=e, e_series=float(np.sqrt(s2)), e_first=float(np.sqrt(f2)),
                       e_mm=W.l2_error(ym, am, SIG, ss), e_ls=W.l2_error(y, W.ls_amplitudes(y, SIG, ss), SIG, ss), y_hs=y,
                       e_fock=W.l2_error(*W.scaled_rule(K, SIG, ss, V['fock']), SIG, ss),
                       e_mercer=W.l2_error(*W.scaled_rule(K, SIG, ss, V['mercer']), SIG, ss), e_best=eb, v_best=vb / V['hs'])
            if e > E_FLOOR:
                yo, ao, eo = W.optimal_fit(y, SIG, ss)
                row.update(e_opt=eo, y_opt=yo, a_opt=ao)
            rows.append(row)
    # A2: the tail's integer ratios and the tail's share
    ratios = {K: [W.gh_defect(2 * K + 2 * j, K) / W.gh_defect(2 * K, K) for j in range(1, 5)] for K in range(1, K_A2 + 1)}
    share = []
    for r_ in TAUS_A:
        t = W.params(SIG, sigma_star(r_))['t']
        for K in (5, 10, 20, 30):
            s2, f2 = W.mehler_series(K, t)
            share.append(dict(r=r_, K=K, share=s2 / f2 - 1))
    narrow = []
    for q in NARROW_A:
        ss, t = q * SIG, (1 - q) / (1 + q)
        for K in range(1, K_A + 1):
            s2, f2 = W.mehler_series(K, t)
            if np.sqrt(s2) < G_FLOOR:
                break
            narrow.append(dict(q=q, K=K, t=t, e=W.narrow_error(*W.narrow_rule(K, SIG, ss), SIG, ss), e_series=float(np.sqrt(s2))))
    return dict(rows=rows, ratios=ratios, share=share, narrow=narrow)


def blocks(R):
    """The relaxed blocks along a path: (k kernels, τ, centres, coefficients, rms, whether the relaxation became stationary) after each
    event (hist[i + 1] is iteration i's block, relaxed[i] its flag)."""
    out = []
    for i, h in enumerate(R['hist'][1:]):
        act = np.flatnonzero(h[2] != 0)
        o = np.argsort(h[1][act])
        out.append(dict(k=int(act.size), tau=float(h[0]), c=h[1][act][o], w=h[2][act][o], rms=float(h[3]),
                        relaxed=bool(R['relaxed'][i]) if i < len(R['relaxed']) else False))
    return out


def sample_fit(Phi_of, f, y0):
    """The least-squares fit with free centres on the samples by variable projection from y0: (centres, coefficients, rms)."""
    from scipy.optimize import least_squares

    def res(c):
        P = Phi_of(c)
        return f - P @ np.linalg.lstsq(P, f, rcond=None)[0]
    sol = least_squares(res, np.asarray(y0, float), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    c = np.sort(sol.x)
    P = Phi_of(c)
    a = np.linalg.lstsq(P, f, rcond=None)[0]
    return c, a, float(np.sqrt(np.mean((f - P @ a) ** 2)))


def vp_fit(Phi_of, f, y0, max_nfev=20000):
    """`sample_fit` with its solve's stop (referee round 65): (centres, coefficients, rms, status, evaluations) — MINPACK's status 1 the
    gradient test, 2 the decrease test, 3 the step test, 4 both, 0 the budget."""
    from scipy.optimize import least_squares

    def res(c):
        P = Phi_of(c)
        return f - P @ np.linalg.lstsq(P, f, rcond=None)[0]
    sol = least_squares(res, np.asarray(y0, float), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=max_nfev)
    c = np.sort(sol.x)
    P = Phi_of(c)
    a = np.linalg.lstsq(P, f, rcond=None)[0]
    return c, a, float(np.sqrt(np.mean((f - P @ a) ** 2))), int(sol.status), int(sol.nfev)


def vp_jacobian(K, f, c):
    """Variable projection's Jacobian in the centres, analytic (Golub–Pereyra): with P = Φ(c), the least-squares amplitudes a = P⁺f, the
    residual r = f − Pa and R = I − PP⁺, ∂r/∂c_i = −a_iR∂_cφ_i − (∂_cφ_i·r)(P⁺)ᵀe_i (referee round 65: central differences lose the
    small singular values)."""
    c = np.asarray(c, float)
    P, D = K.phi(c), K.dphi(c)
    Pp = np.linalg.pinv(P)
    a = Pp @ f
    r = f - P @ a
    return -(D - P @ (Pp @ D)) * a[None, :] - Pp.T * (D.T @ r)[None, :]


def b_job(args):
    r_, var = args
    name, step, shift, nfev, rep_, gt = next(v for v in VARIANTS_B if v[0] == var)
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L, _, c0 = TM.setup()
    c0 = c0 + shift
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, max_events=3 * K_B, tol=0.0, max_entries=K_B + 1,
                       max_doublings=MAX_DOUBLINGS, polish_step=step, exit_kkt=True, polish_nfev=nfev, repolish=rep_,
                       polish_gtol=gt)
    rf = float(np.sqrt(np.mean(f ** 2)))
    out = []
    for i, b in enumerate(blocks(R)):
        if b['k'] == 0:
            continue
        yh, _ = W.half_squeeze(b['k'], SIG, ss)
        cands = [sample_fit(K.phi, f, b['c']), sample_fit(K.phi, f, C_B + yh)]
        c_o, a_o, rms_o = min(cands, key=lambda z: z[2])
        out.append(dict(k=b['k'], tau=b['tau'], c=b['c'], w=b['w'], e=b['rms'] / rf, e_opt=rms_o / rf, c_opt=c_o, relaxed=b['relaxed'],
                        dc=float(np.abs(b['c'] - c_o).max()) if c_o.size == b['c'].size else np.nan,
                        dc_hs=float(np.abs(b['c'] - (C_B + np.sort(yh))).max()),
                        e_hs=W.l2_error(*W.half_squeeze(b['k'], SIG, ss), SIG, ss),
                        polish=R['polish_log'][i] if i < len(R['polish_log']) else []))
    print(f"  B tau/sigma {r_} ({var}): {len(out)} blocks, errors " + ', '.join(f"{o['e']:.2e}" for o in out) +
          f" [{time.time() - t0:.0f}s]", flush=True)
    return dict(r=r_, var=var, blocks=out, events=[(float(e[0]), int(e[1]), e[2], float(e[3])) for e in R['events']],
                relaxed=R['relaxed'], stiff=R['stiff_redone'], redone=R['solver_redone'], failed=R['solver_failed'],
                rejected=R['polish_rejected'], entries=R['entry_info'], secs=time.time() - t0)


def part_B(idx=None):
    """B's paths for the variants of VARIANTS_B with indices idx (all if None), in the order of VARIANTS_B, then TAUS_B."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(r_, v[0]) for i, v in enumerate(VARIANTS_B) if idx is None or i in idx for r_ in TAUS_B], chunksize=1)


def run_B(old, idx):
    """B's part: all paths (idx None), or the variants with indices idx rerun and merged into the stored ones (`merge_B`; referee round
    70: a partial rerun must not drop the others)."""
    return part_B() if idx is None else merge_B(old, part_B(idx))


def merge_B(old, new):
    """B's paths: the new ones replace the old of the same width and variant; in the order of VARIANTS_B, then TAUS_B."""
    names = [v[0] for v in VARIANTS_B]
    key = lambda b: (b['r'], b['var'])
    got = {key(b): b for b in (old if old is not None else [])}
    got.update({key(b): b for b in new})
    return sorted(got.values(), key=lambda b: (names.index(b['var']), TAUS_B.index(b['r'])))


def b2_job(args):
    """B2 (referee round 54): from one of B's non-optimal blocks, variable projection (Levenberg–Marquardt on the centres) started at the
    block's own centres and at the half-squeeze nodes — whether the block lies in the basin of the best optimum found (its relaxation
    unfinished) or at another stationary point; dc, the block's largest distance from the better optimum's nodes; and at the block's
    centres (least-squares amplitudes) the condition numbers of the two problems' Jacobians — the joint fit in (w, c), as the polish
    solves it, and variable projection's in c alone (`vp_jacobian`, analytic), both with the centres in units of σ, and the ratio of their
    smallest singular values (referee rounds 60, 65) — and the same in the polish's own units (the centres in x) and with every column
    scaled to unit norm (referee round 70: the comparison depends on the units); where variable projection from the block stops: its
    status, its closest nodes (σ) and the largest cosine between its residual and a column of its Jacobian (zero at a stationary point)."""
    r_, var, k, c = args
    warnings.filterwarnings('ignore')
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    rf = float(np.sqrt(np.mean(f ** 2)))
    yh, _ = W.half_squeeze(k, SIG, ss)
    cb, ab, eb, st, nf = vp_fit(K.phi, f, c)
    cr, _, er, _, _ = vp_fit(K.phi, f, C_B + yh)
    c = np.sort(np.asarray(c, float))
    P = K.phi(c)
    wl = np.linalg.lstsq(P, f, rcond=None)[0]
    sv = lambda J: np.linalg.svd(J, compute_uv=False)
    eq = lambda J: J / np.where(np.linalg.norm(J, axis=0) > 0, np.linalg.norm(J, axis=0), 1.0)[None, :]
    Jj, Jv = np.c_[P, K.dphi(c) * wl[None, :]], vp_jacobian(K, f, c)  # the polish's own units (the centres in x)
    sj, sp = sv(np.c_[P, K.dphi(c) * wl[None, :] * SIG]), sv(Jv * SIG)
    sjx, spx, sje, spe = sv(Jj), sv(Jv), sv(eq(Jj)), sv(eq(Jv))
    Je, re = vp_jacobian(K, f, cb), f - K.phi(cb) @ ab
    cos = float(np.max(np.abs(Je.T @ re) / (np.linalg.norm(Je, axis=0) * np.linalg.norm(re)))) if np.linalg.norm(re) > 0 else 0.0
    return dict(r=r_, var=var, k=k, e_block=eb / rf, e_rule=er / rf, dc=float(np.abs(c - (cb if eb <= er else cr)).max()),
                cond_joint=float(sj[0] / sj[-1]), cond_vp=float(sp[0] / sp[-1]), smin_ratio=float(sp[-1] / sj[-1]),
                cond_joint_x=float(sjx[0] / sjx[-1]), cond_vp_x=float(spx[0] / spx[-1]), smin_ratio_x=float(spx[-1] / sjx[-1]),
                smin_ratio_eq=float(spe[-1] / sje[-1]), vp_status=st, vp_nfev=nf,
                vp_gap=float(np.diff(cb).min() / SIG) if cb.size > 1 else np.nan, vp_cos=cos)


def seed_rule(K, f, w, c, s, cp, step, wfrac=1e-2):
    """`sliding_iss`'s acceptance of a polish from (w, c) to the centres cp: every kernel whose |w| is at most wfrac·max|w| (near its seed)
    lands within step kernel widths of its start or of its correlation's ascent target (`RT.ascent_target`). Returns whether it is
    accepted and the largest such kernel's move (σ; the nearer of the two)."""
    seedy = np.flatnonzero(np.abs(w) <= wfrac * np.abs(w).max())
    r0 = f - K.phi(c) @ w
    mv = [min(abs(cp[i] - c[i]), abs(cp[i] - RT.ascent_target(K, r0, c[i], s[i])[0])) / K.sig for i in seedy]
    return dict(accepted=bool(all(m <= step for m in mv)), move=float(max(mv)) if mv else 0.0)


B3_BIG = 500_000                                                      # B3 (referee round 70): the larger budget
B3_KEYS = ('path', 'budget', 'no test', 'budget, no test')            # B3's restarts: the path's settings, then one or both changed


class _StopPath(Exception):
    """`path_polishes`' recorder: the call asked for is recorded, the rest of the path is not needed."""


def path_polishes(r_, var, stop_after=None):
    """B's path (r_, var) rerun as `b_job`'s with every `polish_block` call recorded (referee round 70): its start (w, c, s) and settings,
    in the path's order (R['polish_log'] flattened: each call appends one record to its block's log). With stop_after = n the path stops
    once call n (from 0) is recorded — before that call's solve — and R is None. Returns (R, K, f, calls)."""
    name, step, shift, nfev, rep_, gt = next(v for v in VARIANTS_B if v[0] == var)
    L, _, c0 = TM.setup()
    K = RT.Bumps(X_B, SIG)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
    calls, real = [], RT.polish_block

    def rec(K_, fstar, w, c, s, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        calls.append(dict(w=np.array(w, float), c=np.array(c, float), s=np.array(s, float), max_nfev=max_nfev, gtol=gtol))
        if stop_after is not None and len(calls) > stop_after:
            raise _StopPath
        return real(K_, fstar, w, c, s, cbox=cbox, max_nfev=max_nfev, info=info, gtol=gtol)
    RT.polish_block = rec
    try:
        R = RT.sliding_iss(L, K, f, c0 + shift, np.ones(c0.size), polish=True, max_events=3 * K_B, tol=0.0, max_entries=K_B + 1,
                           max_doublings=MAX_DOUBLINGS, polish_step=step, exit_kkt=True, polish_nfev=nfev, repolish=rep_,
                           polish_gtol=gt)
    except _StopPath:
        R = None
    finally:
        RT.polish_block = real
    return R, K, f, calls


def first_nonoptimal(b):
    """B's path b (its stored blocks): the first block with polishes whose error is at least 1 + OPT_TOL times the optimum's, as
    dict(n0 — the number of polish calls before its first —, k, e, e_opt, n_polish, first — its first polish's stored stop —), or None
    (the stored blocks are those with kernels: a block without kernels has no polish, so n0 counts every call)."""
    n0 = 0
    for o in b['blocks']:
        n_i = len(o['polish']) if 'polish' in o else 0
        if o['k'] > 0 and n_i > 0 and abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL:
            return dict(n0=n0, k=int(o['k']), e=float(o['e']), e_opt=float(o['e_opt']), n_polish=n_i, first=dict(o['polish'][0]))
        n0 += n_i
    return None


def b3_job(args):
    """B3 (referee rounds 65, 70): the path's first non-optimal block (its error over the optimum's at least 1 + OPT_TOL), its first
    polish restarted from that polish's recorded start (`path_polishes`) by `polish_block`, one setting at a time — the path's budget and
    gradient test; the budget B3_BIG with the test; the path's budget without it; both changed —: each the solve's stop (TRF's status: 0
    the budget, 1 the gradient test, 2 the decrease test, 3 the step test), its evaluations, its error over the optimum's and the path's
    acceptance of its end (`seed_rule` with the path's polish step: the first polish has no improvement test)."""
    r_, var = args[:2]
    _, step, _, nfev, _, gt = next(v for v in VARIANTS_B if v[0] == var)
    warnings.filterwarnings('ignore')
    t0 = time.time()
    if len(args) > 2:                                                    # the block known from B's raw (`first_nonoptimal`): the path
        blk = args[2]                                                    # rerun only to its first polish
        _, K, f, calls = path_polishes(r_, var, stop_after=blk['n0'])
        if len(calls) <= blk['n0'] or calls[blk['n0']]['c'].size != blk['k']:
            raise RuntimeError(f"B3 {r_} {var}: the rerun's polish {blk['n0']} is not the stored block's ({len(calls)} calls)")
        out = b3_restart(r_, var, K, f, calls[blk['n0']], blk['k'], blk['e_opt'], blk['e'] / blk['e_opt'], blk['n_polish'], step,
                         nfev, gt, t0)
        fs = blk.get('first') or {}                                      # referee round 77: the 'path' restart is the stored polish
        if any(fs.get(k) is not None and out['path'][k] != fs[k] for k in ('status', 'nfev')):
            raise RuntimeError(f"B3 {r_} {var}: the 'path' restart {out['path']['status']}/{out['path']['nfev']} is not the stored "
                               f"first polish {fs.get('status')}/{fs.get('nfev')}")
        return out
    R, K, f, calls = path_polishes(r_, var)
    rf, ss, n0 = float(np.sqrt(np.mean(f ** 2))), sigma_star(r_), 0
    for i, b in enumerate(blocks(R)):
        n_i = len(R['polish_log'][i]) if i < len(R['polish_log']) else 0
        if b['k'] > 0 and n_i > 0:
            yh, _ = W.half_squeeze(b['k'], SIG, ss)
            e_opt = min(sample_fit(K.phi, f, b['c'])[2], sample_fit(K.phi, f, C_B + yh)[2]) / rf
            if abs(b['rms'] / rf / e_opt - 1) >= OPT_TOL:
                return b3_restart(r_, var, K, f, calls[n0], b['k'], e_opt, b['rms'] / rf / e_opt, n_i, step, nfev, gt, t0)
        n0 += n_i
    return dict(r=r_, var=var, k=None, secs=time.time() - t0)


def b3_restart(r_, var, K, f, st, k, e_opt, e_block, n_polish, step, nfev, gt, t0):
    """B3's restarts of one recorded polish start st, one setting at a time (B3_KEYS), each judged by `seed_rule` with the path's
    polish step."""
    rf = float(np.sqrt(np.mean(f ** 2)))
    out = dict(r=r_, var=var, k=k, e_opt=e_opt, e_block=e_block, n_polish=n_polish,
               e_start=float(np.sqrt(np.mean((f - K.phi(st['c']) @ st['w']) ** 2))) / rf / e_opt,
               settings=dict(nfev=st['max_nfev'], gtol=st['gtol']))
    for key, nf, g in zip(B3_KEYS, (nfev, B3_BIG, nfev, B3_BIG), (gt, gt, None, None)):
        wp, cp, info = RT.polish_block(K, f, st['w'], st['c'], st['s'], max_nfev=nf, info=True, gtol=g)
        out[key] = dict(info, ratio=float(np.sqrt(np.mean((f - K.phi(cp) @ wp) ** 2))) / rf / e_opt,
                        **seed_rule(K, f, st['w'], st['c'], st['s'], cp, step))
    out['secs'] = time.time() - t0
    print(f"  B3 tau/sigma {r_} ({var}) K {k}: start {out['e_start']:.4g}, " + ', '.join(
        f"{key} {out[key]['status']}/{out[key]['nfev']} {out[key]['ratio']:.5f}" for key in B3_KEYS)
          + f" [{out['secs']:.0f}s]", flush=True)
    return out


def b3_jobs(B):
    """B3's jobs: each of B's paths with a non-optimal block (its error at least 1 + OPT_TOL times the optimum's) — with the first such
    block that was polished (`first_nonoptimal`) where the stored blocks record their polishes, so that the rerun stops at its first
    polish."""
    out = []
    for b in B:
        if any(abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL for o in b['blocks']):
            blk = first_nonoptimal(b) if all('polish' in o and 'k' in o for o in b['blocks']) else None
            out.append((b['r'], b['var']) if blk is None else (b['r'], b['var'], blk))
    return out


def part_B3(B):
    """B3: `b3_job` for `b3_jobs(B)`."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b3_job, b3_jobs(B), chunksize=1)


def b3_binds(q, ref=1.0):
    """What bounds a B3 block's range end (referee rounds 70, 77), from the restarts' ends that reach the optimum — judged against ref, the
    optimum's error over the stored optimum's (below B4_MAX B4's long-double end) — and that the path's seed rule accepts: the fewest
    settings whose change gives such an end (the budget, the gradient test, either alone, or both); where every end that reaches it is
    rejected, the fewest settings that reach it "(and the radius)", and the radius alone where the path's own settings reach it; "neither"
    with the both-changed restart's stop (its budget spent, or the solve's own stop); an inconsistency where the path's own polish reaches
    it, accepted."""
    reach = lambda k: q[k]['ratio'] / ref < 1 + OPT_TOL
    ok = lambda k: reach(k) and q[k]['accepted']
    name = lambda ks: {('budget',): 'the budget', ('no test',): 'the gradient test'}.get(tuple(ks), 'the budget or the gradient test')
    if reach('path'):
        return 'the radius' if not q['path']['accepted'] else "nothing (inconsistent: the path's own polish reaches it, accepted)"
    for test in (ok, reach):
        one = [k for k in ('budget', 'no test') if test(k)]
        tail = '' if test is ok else ' (and the radius)'
        if one:
            return name(one) + tail
        if test('budget, no test'):
            return 'the budget and the gradient test' + tail
    st = q['budget, no test']['status']
    return 'neither (' + ('the budget spent' if st == 0 else f"the solve's stop, status {st}") + ')'


def b3_lines(B3, B4=None, B=None):
    """B3's table and summary (referee rounds 70, 77): per block the restarts' stops, errors and the seed rule, and what bounds its end —
    judged, below B4_MAX, against B4's long-double end (from the path's block), the better optimum there; with B, the 'path' restarts
    against the first polishes B stored (their stops and evaluations)."""
    if B3 is None or len(B3) == 0:                                    # raw-loaded: an object array
        return []
    rows = [q for q in B3 if q.get('k') is not None]
    refs = {(q['r'], q['k']): q['long double']['ratio'] for q in (B4 if B4 is not None else []) if q['start'] == 'block'}
    ref = lambda q: refs.get((q['r'], q['k']), 1.0)
    cell = lambda x: f"{x['status']}/{x['nfev']}, {x['ratio']:.5f}, {x['move']:.2f}σ {'a' if x['accepted'] else 'r'}"
    out = ['', "### B3. The range ends: the polish restarted from its recorded start, one setting at a time (referee rounds 65, 70, 77)", '',
           f"Each path's first non-optimal block (its error at least {1 + OPT_TOL:g} times the optimum's): the path rerun with every polish's "
           "start recorded, and from the block's first polish's start — its error over the optimum's — the polish's own solve "
           f"(`polish_block`, TRF) with the path's budget and gradient test, with {B3_BIG:g} evaluations and the test, with the path's budget "
           "without the test, and with both changed: each its stop (status 0: the budget; 1: the gradient test; 2: the decrease test; 3: the "
           "step test)/evaluations, its error over the optimum's and the path's seed rule on its end (the largest move of a kernel near its "
           f"seed, σ; a: accepted within the path's polish step, r: rejected); the block the path ended with; † the optimum below {E_FLOOR:g}:",
           '', "| τ/σ | variant | K | optimum | start | the path's block | " + ' | '.join(B3_KEYS) + ' | bounded by |',
           '|---|---|---|---|---|---|' + '---|' * len(B3_KEYS) + '---|']
    for q in rows:
        out.append(f"| {q['r']:g} | {q['var']} | {q['k']} | {q['e_opt']:.2e}{'†' if q['e_opt'] < E_FLOOR else ''} | {q['e_start']:.4g} | "
                   f"{q['e_block']:.4g} | " + ' | '.join(cell(q[k]) for k in B3_KEYS) + f" | {b3_binds(q, ref(q))} |")
    up = [q for q in rows if q['e_opt'] >= E_FLOOR]
    cnt = {}
    for q in up:
        cnt.setdefault(b3_binds(q, ref(q)), []).append(f"{q['r']:g} ({q['var']})")
    lo = [q for q in rows if q['e_opt'] < E_FLOOR]
    both = lambda q: q['budget, no test']
    if B is not None and len(B):
        first = {(b['r'], b['var']): (first_nonoptimal(b) or {}).get('first') for b in B}
        cmp = [(q, first.get((q['r'], q['var']))) for q in rows]
        cmp = [(q, fs) for q, fs in cmp if fs is not None]
        same = [q for q, fs in cmp if (q['path']['status'], q['path']['nfev']) == (fs.get('status'), fs.get('nfev'))]
        out += ['', f"The 'path' restart reproduces the first polish B stored — its stop and evaluations — at {len(same)} of {len(cmp)} blocks"
                + (" (not at " + '; '.join(f"{q['r']:g} ({q['var']})" for q, _ in cmp if q not in same) + ")" if len(same) < len(cmp) else '')
                + "."]
    out += ['', f"Above {E_FLOOR:g} ({len(up)} blocks: τ/σ (variant)) the end is bounded by", '']
    out += [f"- {k}: {len(v)} — {'; '.join(v)}" + (';' if i < len(cnt) - 1 else '.') for i, (k, v) in enumerate(sorted(cnt.items()))]
    if lo:
        out += ['', f"Below it ({len(lo)} blocks) the restart with both changed stops at " + ', '.join(
            f"{both(q)['ratio']:.4f} ({both(q)['status']}/{both(q)['nfev']}: {b3_binds(q, ref(q))})" for q in lo) + " times the stored "
            "optimum" + (", " + ', '.join(f"{both(q)['ratio'] / ref(q):.4f}" for q in lo) + " times B4's long-double end" if all(
                (q['r'], q['k']) in refs for q in lo) else '') + "."]
    by = {}
    for q in rows:
        by.setdefault(float(f"{q['e_opt']:.3g}"), []).append(both(q)['ratio'] / ref(q) - 1)
    out += ['', "The restart with both changed (the polish's TRF): its end over the optimum, minus 1 — below " + f"{B4_MAX:g} over B4's "
            "long-double end —, by the optimum (one solver from one kind of start; a negative value is the end below the stored optimum, "
            "converged only to that): " + '; '.join(
                f"{e:.3g}: {min(v):+.1e}" + (f" … {max(v):+.1e}" if max(v) > min(v) else '') for e, v in sorted(by.items(), reverse=True))
            + '.', '']
    return out



B4_CAP, B4_MAX = 400_000, 2e-9                                        # B4 (referee round 70): the solve's iterations; the blocks below
B4_PREC = (('double', np.float64), ('long double', np.longdouble))    # B4: the parameters' and the residual's arithmetic
B4_ANCHOR, B4_FLAT = 2000, 1e-10                                       # B4 (referee round 77): the anchored double's re-anchoring period
                                                                       # (accepted steps) and its stop: ½‖r‖² over a period falls by less
B4_MARKS = (100_000, 200_000, 300_000)                                 # B4 (round 77): the anchored double's end after these iterations


def lm_fit(f, w, c, dtype, max_iter=B4_CAP, mu0=1e-20, anchored=False, anchor_every=B4_ANCHOR, marks=None):
    """B4 (referee round 70): Levenberg–Marquardt on the joint fit in (w, c) on X_B — the parameters and the residual f − Φ(c)w in dtype
    (np.float64 or np.longdouble), the Jacobian and the damped step in double: the least-squares problem [J; μ^{1/2}D]d ≈ [−r; 0], D the
    columns' largest norms so far (its condition is the Jacobian's, not its square's). A step is accepted when the objective ½‖r‖²
    computed in dtype decreases; Nielsen's damping from μ = mu0, near Gauss–Newton — the slow modes' curvature is ~10⁻¹⁷ of the largest,
    and a larger damping shrinks their steps' predicted decrease below the objective's rounding. Stops after max_iter iterations ('cap')
    or once μ exceeds 10¹⁰⁰ without a decrease ('stall': no step the solve computes decreases the objective as its residual evaluates
    it). Referee round 77, anchored (double only): the residual about an anchor (w₀, c₀), r = r₀ − Σ_k φ₀_k[(w_k − w₀_k) +
    w_k·expm1(Δ_k)], Δ_k = (c_k − c₀_k)(2x − c_k − c₀_k)/(2σ²), r₀ = f − Φ(c₀)w₀ by direct subtraction — its rounding fixed, the change
    computed to its own relative accuracy —, re-anchored every anchor_every accepted steps and at a stall; it stops at a stall right
    after an anchor, or ('flat') once a period's accepted steps lower the objective, as that anchor computes it, by less than B4_FLAT of
    itself. marks: a dict whose keys are iteration counts, given the parameters after that many iterations. Returns (the parameters
    (w, c) in dtype, iterations, accepted steps, stop)."""
    x, fd = X_B.astype(dtype), np.asarray(f, float).astype(dtype)
    n = len(w)
    p = np.r_[np.asarray(w, float), np.asarray(c, float)].astype(dtype)
    s2 = 2 * dtype(SIG) ** 2
    anc = {}
    if anchored:
        if dtype is not np.float64:
            raise ValueError("lm_fit: the anchored residual is double's")

        def set_anchor(p_):
            anc['w'], anc['c'] = p_[:n].copy(), p_[n:].copy()
            anc['P'] = np.exp(-(X_B[:, None] - anc['c'][None, :]) ** 2 / (2 * SIG ** 2))
            anc['r'] = np.asarray(f, float) - anc['P'] @ anc['w']
            anc['acc'], anc['F'] = 0, float(anc['r'] @ anc['r'] / 2)

        def resid(p_):
            dl = (p_[None, n:] - anc['c'][None, :]) * (2 * X_B[:, None] - p_[None, n:] - anc['c'][None, :]) / (2 * SIG ** 2)
            return anc['r'] - (anc['P'] * ((p_[:n] - anc['w'])[None, :] + p_[None, :n] * np.expm1(dl))).sum(axis=1)
        set_anchor(p)
    else:
        resid = lambda p_: fd - np.exp(-(x[:, None] - p_[None, n:]) ** 2 / s2) @ p_[:n]

    def jac(p_):
        q = p_.astype(float)
        u = X_B[:, None] - q[None, n:]
        P = np.exp(-u ** 2 / (2 * SIG ** 2))
        return -np.c_[P, P * u / SIG ** 2 * q[None, :n]]
    r = resid(p)
    F, J = r @ r / 2, jac(p)
    D = np.linalg.norm(J, axis=0)
    mu, nu, it, acc = mu0, 2.0, 0, 0
    while it < max_iter:
        if marks is not None and it in marks:
            marks[it] = p.copy()
        it += 1
        r64 = r.astype(float)
        d = np.linalg.lstsq(np.r_[J, np.sqrt(mu) * np.diag(D)], np.r_[-r64, np.zeros(2 * n)], rcond=None)[0]
        Jd = J @ d
        pred = -(Jd @ r64) - 0.5 * (Jd @ Jd)
        pn = p + d.astype(dtype)
        rn = resid(pn)
        Fn = rn @ rn / 2
        if pred > 0 and Fn < F:
            rho = float(F - Fn) / pred
            p, r, F, acc = pn, rn, Fn, acc + 1
            mu, nu = max(mu * max(1 / 3, 1 - (2 * rho - 1) ** 3), 1e-40), 2.0
            J = jac(p)
            D = np.maximum(D, np.linalg.norm(J, axis=0))
            if anchored:
                anc['acc'] += 1
                if anc['acc'] >= anchor_every:                          # re-anchor: r₀ afresh at the current point
                    if anc['F'] - F < B4_FLAT * anc['F']:
                        return p, it, acc, 'flat'
                    set_anchor(p)
                    r = resid(p)
                    F = r @ r / 2
        else:
            mu, nu = mu * nu, 2 * nu
            if mu > 1e100:
                if anchored and anc['acc'] > 0:                         # a stall since the last anchor: re-anchor and go on
                    set_anchor(p)
                    r = resid(p)
                    F, mu, nu = r @ r / 2, mu0, 2.0
                    continue
                return p, it, acc, 'stall'
    return p, it, acc, 'cap'


def lm_rms(f, p, n):
    """The rms of f − Φ(c)w on X_B in long double, p = (w, c) with n kernels."""
    x, fd, p = X_B.astype(np.longdouble), np.asarray(f, float).astype(np.longdouble), np.asarray(p).astype(np.longdouble)
    r = fd - np.exp(-(x[:, None] - p[None, n:]) ** 2 / (2 * np.longdouble(SIG) ** 2)) @ p[:n]
    return float(np.sqrt(np.mean(r * r)))


def b4_job(args):
    """B4 (referee round 70): from one start — the path's block ('block') or the stored optimum's centres with least-squares amplitudes
    ('optimum') — `lm_fit` with the residual in double and in long double (B4_PREC): the start's error and each end's over the stored
    optimum's (the rms in long double), its iterations, accepted steps and stop, and long double's ε (an error if not extended); referee
    round 77: the same solve in double with the residual about an anchor ('anchored'), and the polish's own TRF with B3_BIG evaluations
    and no gradient test ('trf': its status and evaluations); the anchored end's error also after B4_MARKS iterations ('marks')."""
    r_, var, k, start, w0, c0, e_opt = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
    rf = float(np.sqrt(np.mean(f ** 2)))
    out = dict(r=r_, var=var, k=k, start=start, e_opt=e_opt, e_start=lm_rms(f, np.r_[w0, c0], k) / rf / e_opt,
               eps=float(np.finfo(np.longdouble).eps))
    if out['eps'] >= np.finfo(float).eps:
        raise RuntimeError("B4: long double is not extended on this platform")
    for key, dt in B4_PREC:
        p, it, acc, stop = lm_fit(f, w0, c0, dt)
        out[key] = dict(ratio=lm_rms(f, p, k) / rf / e_opt, it=it, acc=acc, stop=stop)
    marks = dict.fromkeys(B4_MARKS)
    p, it, acc, stop = lm_fit(f, w0, c0, np.float64, anchored=True, marks=marks)  # referee round 77: double, the residual about an anchor
    out['anchored'] = dict(ratio=lm_rms(f, p, k) / rf / e_opt, it=it, acc=acc, stop=stop,
                           marks=[[m, lm_rms(f, q, k) / rf / e_opt] for m, q in marks.items() if q is not None])
    K = RT.Bumps(X_B, SIG)                                              # the polish's own solve from the same start: TRF, B3's budget,
    wp, cp, info = RT.polish_block(K, f, np.asarray(w0, float), np.asarray(c0, float), np.ones(k), max_nfev=B3_BIG, info=True,
                                   gtol=None)                          # no gradient test
    out['trf'] = dict(ratio=lm_rms(f, np.r_[wp, cp], k) / rf / e_opt, status=int(info['status']), nfev=int(info['nfev']))
    out['secs'] = time.time() - t0
    print(f"  B4 tau/sigma {r_} ({var}) K {k} from the {start}: {out['e_start']:.6g} → " + ', '.join(
        f"{key} {out[key]['ratio']:.6f} ({out[key]['it']}, {out[key]['stop']})" for key in [k for k, _ in B4_PREC] + ['anchored'])
          + f", TRF {out['trf']['ratio']:.6f} ({out['trf']['status']}/{out['trf']['nfev']}) [{out['secs']:.0f}s]",
          flush=True)
    return out


def b4_jobs(B, B3):
    """B4's jobs: each distinct block (width, K) of B3 whose optimum is below B4_MAX — B3's first row for it, its block in B's path (the
    block with that K whose optimum is the row's) —, from the path's block and from the stored optimum."""
    out, seen = [], set()
    K = RT.Bumps(X_B, SIG)
    for q in (B3 if B3 is not None else []):
        if q.get('k') is None or q['e_opt'] >= B4_MAX or (q['r'], q['k']) in seen:
            continue
        seen.add((q['r'], q['k']))
        b = next(b for b in B if b['r'] == q['r'] and b['var'] == q['var'])
        o = min((o for o in b['blocks'] if o['k'] == q['k']), key=lambda o: abs(o['e_opt'] / q['e_opt'] - 1))
        f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(q['r']) ** 2))
        wo = np.linalg.lstsq(K.phi(o['c_opt']), f, rcond=None)[0]
        out += [(q['r'], q['var'], q['k'], 'block', np.asarray(o['w'], float), np.asarray(o['c'], float), q['e_opt']),
                (q['r'], q['var'], q['k'], 'optimum', wo, np.asarray(o['c_opt'], float), q['e_opt'])]
    return out


def part_B4(B, B3):
    """B4: `b4_job` for `b4_jobs(B, B3)`."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b4_job, b4_jobs(B, B3), chunksize=1)


def b4_lines(B4):
    """B4's table and summary (referee rounds 70, 77): the same solve in double and in long double from the blocks near the floor, in
    double about an anchor, and the polish's own TRF; every end against the long-double end from the path's block."""
    if B4 is None or len(B4) == 0:                                    # raw-loaded: an object array
        return []
    cell = lambda x: f"{x['ratio']:.6f}, {x['it']}/{x['acc']}, {x['stop']}"
    tcell = lambda x: f"{x['ratio']:.6f}, {x['status']}/{x['nfev']}"
    has = 'anchored' in B4[0]
    out = ['', "### B4. The floor: the same solve with the residual in double, in long double" + (" and about an anchor (referee rounds 70, 77)"
                                                                                                  if has else " (referee round 70)"), '',
           f"B3's blocks whose optimum is below {B4_MAX:g}, from the path's block and from the stored optimum (its centres with "
           "least-squares amplitudes): Levenberg–Marquardt on the joint fit in (w, c) (`lm_fit`: the Jacobian and the step in double; the "
           f"parameters and the residual f − Φ(c)w, by direct subtraction, in double or in long double (its ε {B4[0]['eps']:.3g})"
           + (f", or in double about an anchor (w₀, c₀) — r = r₀ − Σ_kφ₀_k[(w_k − w₀_k) + w_k·expm1(Δ_k)], Δ_k = (c_k − c₀_k)(2x − c_k − c₀_k)/"
              f"(2σ²), re-anchored every {B4_ANCHOR:g} accepted steps and at a stall —" if has else '')
           + "; a step accepted when the objective so computed decreases)" + (f" and the polish's own TRF (`polish_block`, {B3_BIG:g} "
           "evaluations, no gradient test)" if has else '') + ": the start's error and each end's over the optimum's (the rms in long double), "
           "its iterations/accepted steps and its stop (stall: no step the solve computes decreases the objective as its residual "
           f"evaluates it; cap: {B4_CAP:g} iterations" + (f"; flat: an anchor's {B4_ANCHOR:g} accepted steps lower it by less than "
           f"{B4_FLAT:g} of itself" if has else '') + ")" + ("; TRF's cell its status/evaluations" if has else '') + "; † the optimum "
           f"below {E_FLOOR:g}:", '',
           "| τ/σ | variant | K | optimum | start | its error | " + ' | '.join(k for k, _ in B4_PREC)
           + (' | double, anchored | TRF |' if has else ' |'), '|---|---|---|---|---|---|' + '---|' * (len(B4_PREC) + 2 * has)]
    for q in B4:
        out.append(f"| {q['r']:g} | {q['var']} | {q['k']} | {q['e_opt']:.2e}{'†' if q['e_opt'] < E_FLOOR else ''} | {q['start']} | "
                   f"{q['e_start']:.6g} | " + ' | '.join(cell(q[k]) for k, _ in B4_PREC)
                   + (f" | {cell(q['anchored'])} | {tcell(q['trf'])} |" if has else ' |'))
    key = lambda q: (q['r'], q['k'])
    blk = sorted((q for q in B4 if q['start'] == 'block'), key=lambda q: -q['e_opt'])
    opt = {key(q): q for q in B4 if q['start'] == 'optimum'}
    ld = {key(q): q['long double']['ratio'] for q in blk}
    lab = lambda q: f"{q['e_opt']:.3g} (τ/σ = {q['r']:g}, K = {q['k']})"
    both = [q for q in blk if key(q) in opt]
    over = lambda q, k: q[k]['ratio'] / ld[key(q)] - 1
    out += ['', "Each end over the long-double end from the path's block, minus 1, by the optimum — from the path's blocks, double: " + '; '.join(
        f"{lab(q)}: {over(q, 'double'):+.1e}" for q in blk) + (" — from the stored optima, double: " + '; '.join(
            f"{lab(q)}: {over(opt[key(q)], 'double'):+.1e}" for q in both) if both else '') + '.'
            + (" In long double the ends from the block and from the stored optimum agree within "
               f"{max(abs(q['long double']['ratio'] / opt[key(q)]['long double']['ratio'] - 1) for q in both):.1e}; the stored optimum over "
               "that end, minus 1: " + '; '.join(f"{lab(q)}: {1 / ld[key(q)] - 1:+.1e}" for q in both) + '.' if both else '')]
    if has:
        an = [over(q, 'anchored') for q in B4 if key(q) in ld]
        cap = [q for q in B4 if key(q) in ld and q['anchored']['stop'] == 'cap']
        late = lambda q: [x for m, x in q['anchored']['marks'] if m >= B4_MARKS[1]] + [q['anchored']['ratio']]
        spread = lambda q: (max(late(q)) - min(late(q))) / ld[key(q)]
        out += ['', f"The anchored double, from every start, over the long-double end from the path's block, minus 1: {min(an):+.1e} … "
                f"{max(an):+.1e} (its stops: " + ', '.join(f"{k} {sum(q['anchored']['stop'] == k for q in B4)}" for k in ('flat', 'stall', 'cap')
                                                            if any(q['anchored']['stop'] == k for q in B4)) + ")"
                + ("; at the cap, the same after " + ', '.join(f"{m:g}" for m, _ in cap[0]['anchored']['marks']) + f" and {B4_CAP:g} "
                   "iterations: " + '; '.join(f"{lab(q)} from the {q['start']}: " + ', '.join(
                       f"{x / ld[key(q)] - 1:+.2e}" for _, x in q['anchored']['marks']) + f", {over(q, 'anchored'):+.2e}" for q in cap)
                   + f" — over its last {B4_CAP - B4_MARKS[1]:g} iterations each end moves within " + f"{max(spread(q) for q in cap):.1e}"
                   if cap else '') + '.', '',
                "The polish's TRF from the path's blocks, over the long-double end, minus 1: " + '; '.join(
                    f"{lab(q)}: {over(q, 'trf'):+.1e} ({q['trf']['status']}/{q['trf']['nfev']})" for q in blk)
                + ("; from the stored optima: " + '; '.join(f"{lab(q)}: {over(opt[key(q)], 'trf'):+.1e}" for q in both) if both else '') + '.']
    out.append('')
    return out

def gh_defect_v(n, K, v):
    """D_{n,K}(v) = E He_n(√v Z) − Σ_k w_k He_n(√v z_k), exactly (v rational, n even; referee round 65): He_n(√v z) = Σ_m c_m v^{m/2}z^m
    reduced modulo He_K(z) — the rule integrates the remainder, of degree < K, exactly, and He_K vanishes at its nodes —, and
    E He_n(√v Z) = (v − 1)^{n/2}(n − 1)!!. At v = 2 it is `wide_feature.gh_defect`."""
    v = Fraction(v)
    c = W._he(n)
    p = [Fraction(c[m]) * v ** (m // 2) if m % 2 == 0 else Fraction(0) for m in range(len(c))]
    q = W._he(K)
    for d in range(len(p) - 1, K - 1, -1):
        f_ = p[d]
        if f_:
            for i in range(K + 1):
                p[d - K + i] -= f_ * q[i]
    QK = sum(p[m] * W._gauss_moment(m) for m in range(min(K, len(p))))
    return (v - 1) ** (n // 2) * W._gauss_moment(n) - QK


def frame(rho, sig, tau2):
    """The Mehler frame of ratio ρ of the kernel e^{−(y−y′)²/(4σ²)} (referee round 65): k(y, y′) = g(y)g(y′)Σ_n (ρⁿ/n!)He_n(y/ℓ)He_n(y′/ℓ)
    with ℓ² = 2σ²ρ/(1 − ρ²); the Gauss–Hermite rule of the frame's own Gaussian N(0, v_ρ) — g times N(0, τ²) —, v_ρ = 2σ²τ²/(2σ² +
    (1 − ρ)τ²), has no defect below 2K in it; the measure's variance in the frame's units v′ = v_ρ/ℓ². Returns (v_ρ, v′)."""
    v = 2 * sig ** 2 * tau2 / (2 * sig ** 2 + (1 - rho) * tau2)
    return v, v / (2 * sig ** 2 * rho / (1 - rho ** 2))


def part_A5(A):
    """A5 (referee round 65): per τ/σ the frames at ρ = t (the half-squeeze rule), q (the Mercer expansion) and λ (v′ = 1): v_ρ over the
    half-squeeze variance and v′; the Mercer frame's rule's (N(0, 2τ²/(β² + 1))) errors at A's K; and the first tail ratio in a frame of
    variance ratio v′, D_{2K+2}/D_{2K} against −v′K(K + 2) + (v′ − 1)(K + 1)(2K + 1), exactly."""
    rows = []
    for r_ in TAUS_A:
        ss = sigma_star(r_)
        P = W.params(SIG, ss)
        tau2 = P['tau'] ** 2
        b2 = np.sqrt(1 + 2 * tau2 / SIG ** 2)
        q = (b2 - 1) / (b2 + 1)
        vhs = tau2 * SIG / ss
        fr = {nm: (lambda vv: dict(rho=rho, v=vv[0] / vhs, vp=vv[1]))(frame(rho, SIG, tau2)) for nm, rho in (('t', P['t']), ('q', q),
                                                                                                          ('λ', P['lam']))}
        Ks = [x['K'] for x in A['rows'] if x['r'] == r_]
        vM = 2 * tau2 / (b2 + 1)
        rows.append(dict(r=r_, t=P['t'], q=q, frames=fr, vM=vM / vhs, Ks=Ks,
                         e_mframe=[W.l2_error(*W.scaled_rule(K, SIG, ss, vM), SIG, ss) for K in Ks]))
    chk = [dict(v=str(v), K=K, ok=bool(gh_defect_v(2 * K + 2, K, v) / gh_defect_v(2 * K, K, v)
                                       == -v * K * (K + 2) + (v - 1) * (K + 1) * (2 * K + 1)))
           for v in (Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3), Fraction(7, 4)) for K in range(1, 13)]
    return dict(rows=rows, check=chk)


def a5_lines(A5):
    """A5's table and the tail-ratio check."""
    if A5 is None:
        return []
    out = ['', '### A5. The frames: which Gauss–Hermite rule a Mehler expansion singles out', '',
           "The Mehler frames of the kernel e^{−(y−y′)²/(4σ²)} (ratio ρ, Hermite scale ℓ² = 2σ²ρ/(1 − ρ²)): the Gauss–Hermite rule of the "
           "frame's own Gaussian N(0, v_ρ), v_ρ = 2σ²τ²/(2σ² + (1 − ρ)τ²), has no defect below 2K in the frame; v′ = v_ρ/ℓ² the measure's "
           "variance in the frame's units. At ρ = t (the half-squeeze rule), q (the Mercer expansion) and λ: v_ρ over the half-squeeze "
           "variance, and v′; the Mercer frame's rule N(0, 2τ²/(β² + 1)) — not the Mercer-zero rule — and its ratio per kernel at the last "
           "three K against 1/t and 1/q:", '',
           '| τ/σ | ρ = t: v/v_hs, v′ | ρ = q: v/v_hs, v′ | ρ = λ: v/v_hs, v′ | the Mercer frame\'s rule: ratio per kernel (1/t, 1/q) |',
           '|---|---|---|---|---|']
    for x in A5['rows']:
        f_ = x['frames']
        e = x['e_mframe']
        out.append(f"| {x['r']:g} | " + ' | '.join(f"{f_[k]['v']:.6f}, {f_[k]['vp']:.6f}" for k in ('t', 'q', 'λ'))
                   + f" | {', '.join(f'{e[i - 1] / e[i]:.4f}' for i in range(len(e) - 3, len(e)))} at K = {x['Ks'][-3]}–{x['Ks'][-1]} "
                   f"({1 / x['t']:.4g}, {1 / x['q']:.4g}) |")
    on = [(x['r'], 2 / (x['q'] * (2 - x['frames']['q']['vp'])), x['Ks'][-1]) for x in A5['rows']]
    out += ['', "The Mercer frame's rule in its own frame: its first tail ratio grows like (v′ − 2)K², so its second term exceeds its first "
            "(the ratio of the two ≈ q²(2 − v′)²K²/4) from K ≈ 2/(q(2 − v′)): " + '; '.join(
                f"{k:.1f} at τ/σ = {r_:g} (the last K tested {kl})" for r_, k, kl in on) + " — the ratio per kernel shows its rate only "
            "beyond that K (referee round 77)."]
    ok = [c_ for c_ in A5['check'] if c_['ok']]
    out += ['', f"The first tail ratio D_{{2K+2}}/D_{{2K}} = −v′K(K + 2) + (v′ − 1)(K + 1)(2K + 1) (`gh_defect_v`, exact) at "
                f"v′ = {', '.join(sorted({c_['v'] for c_ in A5['check']}, key=lambda s: Fraction(s)))}, K = 1 … "
                f"{max(c_['K'] for c_ in A5['check'])}: " + (f"exact in all {len(ok)}" if len(ok) == len(A5['check'])
                                                              else f"FAILS in {len(A5['check']) - len(ok)} of {len(A5['check'])}")
            + " — −(K − 1) at v′ = 2, growing like (v′ − 2)K² at any other v′.", '']
    return out


def part_B2(B):
    """B2: `b2_job` from every non-optimal block of B."""
    from multiprocessing import Pool
    jobs = [(b['r'], b['var'], o['k'], np.asarray(o['c'])) for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL]
    with Pool(N_WORKERS) as pool:
        return pool.map(b2_job, jobs, chunksize=1)


def b2_line(B2):
    """B2's summary: the non-optimal blocks from which variable projection reaches the best optimum found (within OPT_TOL), counted where
    that optimum is at least E_FLOOR — below it the solve stops short of the optimum (A), and a miss says nothing of the block's basin
    (referee round 60); of the misses, those within 1% of the optimum (in its basin, unconverged) and the others with the block's largest
    distance from the optimum's nodes (σ); the same for each path's first non-optimal block."""
    if B2 is None or not len(B2):
        return []
    best = [min(q['e_block'], q['e_rule']) for q in B2]
    hit = [abs(q['e_block'] / b - 1) < OPT_TOL for q, b in zip(B2, best)]
    up = [b >= E_FLOOR for b in best]
    miss = [(q['e_block'] / b, q.get('dc', np.nan) / SIG, q) for q, b, h, u in zip(B2, best, hit, up) if u and not h]
    near = sorted(v for v, d, q in miss if v < 1.01)
    far = [(v, d, q) for v, d, q in miss if v >= 1.01]
    seen, first = set(), []                                            # each path's first non-optimal block (B2 keeps the paths' order)
    for i, q in enumerate(B2):
        if (q['r'], q['var']) not in seen:
            seen.add((q['r'], q['var']))
            first.append(i)
    fu = [i for i in first if up[i]]
    fw = [B2[i]['e_block'] / best[i] for i in fu if not hit[i]]
    ends = lambda n: f"{n} {'ends' if n == 1 else 'end'}"
    span = lambda v, f_: f"{min(v):{f_}}" if min(v) == max(v) else f"{min(v):{f_}}–{max(v):{f_}}"
    nb = len(first) - len(fu)
    sp_ = lambda k, f_: span([B2[i][k] for i in fu], f_)
    md_ = lambda k: f"{np.median([B2[i][k] for i in fu]):.3g}"
    return [f"Variable projection from the non-optimal block's own centres (B2) reaches the best optimum found in "
            f"{sum(h for h, u in zip(hit, up) if u)} of the {sum(up)} blocks whose optimum is at least {E_FLOOR:g}"
            + (f"; {ends(len(near))} within 1% of it ({', '.join(f'{v:.4g}' for v in near)} times: in its basin, unconverged)" if near else '')
            + (f"; {ends(len(far))} {span([v for v, d, q in far], '.4g')} times above it — the blocks' largest node distances from the "
               f"optimum's {span([d for v, d, q in far], '.2f')}σ —" + (
                   f" where variable projection stalls: its two closest nodes {span([q['vp_gap'] for v, d, q in far], '.3f')}σ apart, its "
                   f"end not stationary (the largest cosine of its residual with a Jacobian column "
                   f"{span([q['vp_cos'] for v, d, q in far], '.1e')}), its stop "
                   f"{', '.join(f'{q['vp_status']}/{q['vp_nfev']}' for v, d, q in far)} (status/evaluations: a miss says nothing of the "
                   f"block's basin)" if 'vp_gap' in far[0][2] else '') if far else '')
            + f". The {len(B2) - sum(up)} blocks with an optimum below {E_FLOOR:g}, where the solve stops short of the optimum (A), are not "
            f"counted ({sum(h for h, u in zip(hit, up) if not u)} reached). From the first non-optimal block of each path: "
            f"{sum(hit[i] for i in fu)} of {len(fu)}" + (f" (the others {', '.join(f'{v:.4g}' for v in sorted(fw))} times above it)" if fw else '')
            + (f"; on {nb} {'path' if nb == 1 else 'paths'} the first non-optimal block's optimum is below {E_FLOOR:g}" if nb else '') + "."
            + (f" There the Jacobian's condition number is {sp_('cond_joint_x', '.2g')} for the joint fit in (w, c) with the centres in x, "
               f"the polish's parameters ({sp_('cond_joint', '.2g')} with the centres in units of σ) and {sp_('cond_vp_x', '.2g')} for variable projection's "
               f"in c alone (analytic; the same in either unit of the centres); the ratio of their smallest singular values, variable "
               f"projection's over the joint fit's, is at least 1 in any common scaling of the centres for Kaufman's form of the former, "
               f"the centre columns projected off the amplitudes' (referee round 77: the joint fit's smallest singular value is at most "
               f"theirs — a Schur complement), from which the analytic Golub–Pereyra form used here differs by a term proportional to the "
               f"residual; its size depends on the units: {sp_('smin_ratio', '.3g')} (median {md_('smin_ratio')}) with the "
               f"centres in units of σ, {sp_('smin_ratio_x', '.3g')} (median {md_('smin_ratio_x')}) in x, "
               f"{sp_('smin_ratio_eq', '.3g')} (median {md_('smin_ratio_eq')}) with every column scaled to unit norm, which scales the two "
               f"problems' centre columns differently; TRF itself works in Coleman–Li's variables, the bounded amplitudes scaled by "
               f"v^{{1/2}}." if fu and 'smin_ratio_eq' in B2[fu[0]] else ''),
            '']


def cs_margin(r_, var, o, gtol=1e-9):
    """The stationarity test's threshold over its Cauchy–Schwarz bound — an effective cosine (referee round 65): the test is absolute,
    |⟨rφ_i⟩| and |w_i⟨r∂_cφ_i⟩| at most gtol·max_j|⟨f*φ_j⟩| (the starting grid's correlations, gtol the path's default), and by
    Cauchy–Schwarz |⟨rφ_i⟩| ≤ rms(r)·rms(φ_i), |⟨r∂_cφ_i⟩| ≤ rms(r)·rms(∂_cφ_i): the test passes whatever the state once the ratio is ≥ 1,
    and below 1 it asks the residual's cosine with the block's directions to be below the ratio."""
    shift = next(v for v in VARIANTS_B if v[0] == var)[2]
    K = RT.Bumps(X_B, SIG)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
    scale = float(np.abs(K.phi(TM.setup()[2] + shift).T @ f).max()) / X_B.size
    c, w = np.asarray(o['c'], float), np.asarray(o['w'], float)
    lim = np.maximum(np.sqrt(np.mean(K.phi(c) ** 2, axis=0)), np.abs(w) * np.sqrt(np.mean(K.dphi(c) ** 2, axis=0)))
    return float(gtol * scale / (o['e'] * float(np.sqrt(np.mean(f ** 2))) * lim.max()))


def cs_vacuous(r_, var, o, gtol=1e-9):
    """Whether a block passes `sliding_iss`'s stationarity test whatever its state (referee round 60): `cs_margin` ≥ 1."""
    return bool(cs_margin(r_, var, o, gtol) >= 1)


def gap_share(x):
    """A3: the share of the half-squeeze rule's gap to the optimal fit that the best one-parameter variance recovers; an optimum the
    one-parameter rule beats is unconverged (referee round 54)."""
    ref = min(x['e_opt'], x['e_best'])
    sh = (x['e_hs'] - x['e_best']) / (x['e_hs'] - ref) if x['e_hs'] > ref else np.nan
    return f"{sh:.1%}" + (', optimum unconverged' if x['e_best'] < x['e_opt'] * (1 - UNCONV) else '')


def gap_range(A):
    """The range of `gap_share` over A's rows with an optimum, and how many optima the variance beats."""
    xs = [x for x in A['rows'] if 'e_opt' in x and x['e_hs'] > min(x['e_opt'], x['e_best'])]
    if not xs:
        return '—'
    sh = [(x['e_hs'] - x['e_best']) / (x['e_hs'] - min(x['e_opt'], x['e_best'])) for x in xs]
    return f"{min(sh):.1%}–{max(sh):.1%} ({sum(x['e_best'] < x['e_opt'] * (1 - UNCONV) for x in xs)} optima unconverged)"


def polish_cell(o, rf=None):
    """A block's polishes: status/evaluations/accepted (status 0: the budget spent first; referee round 54), and for a repeated polish
    (its end's rms recorded) that end's error over the optimum's in brackets (referee round 65: a gradient-test stop at or above it;
    B3 restarts the solve from the block)."""
    return ' '.join(f"{p['status']}/{p['nfev']}/{'a' if p['accepted'] else 'r'}"
                    + (f" ({p['e1'] / (o['e_opt'] * rf):.4g})" if rf and 'e1' in p else '') for p in o.get('polish', [])) or '—'


def entry_lines(B):
    """B's entries (`sliding_iss`'s entry_info) by the block they join (referee round 60: the regimes apart) — none (a path's first entry:
    every dual zero, so the least time to the wall is the largest rate), an optimal block, a non-optimal one (block it − 1, the entries in
    the order of the events) — and over all: whether the entrant had the largest current rate, its rate over the largest, its dual before
    the entry (the entrant is the kernel of the least time to the wall (1 − |u_j|)/|⟨rφ_j⟩|: referee round 54), and where it starts
    against the block it joins: beyond its outermost node (how far) or between its nodes (their duals)."""
    rows = {'none (the first entry)': [], 'an optimal block': [], 'a non-optimal block': []}
    for b in B:
        info = b.get('entries') if b.get('entries') is not None else []
        shift = next(v for v in VARIANTS_B if v[0] == b['var'])[2]
        for e, ev in zip(info, [ev for ev in b['events'] if ev[2] == 'enter']):
            key = (b['r'], shift, e['it'], e['j'])
            if e['it'] == 0:
                rows['none (the first entry)'].append((e, None, key))
            elif e['it'] <= len(b['blocks']):
                o = b['blocks'][e['it'] - 1]
                nodes = np.asarray(o['c'])
                rows['an optimal block' if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL else 'a non-optimal block'].append(
                    (e, max(nodes.min() - ev[3], ev[3] - nodes.max()) / SIG, key))
    rows['all'] = [x for v in list(rows.values()) for x in v]
    if not rows['all']:
        return []
    groups = {}
    for e, d, key in rows['all']:
        groups.setdefault(key, []).append(abs(e['u0']))
    cid, spread = {}, 0.0                                              # an event: one key, its duals within 1e-4 (the variants part)
    for key, us in groups.items():
        us_, c_, lo_ = sorted(set(us)), 0, None
        for i, u in enumerate(us_):
            if i and u - us_[i - 1] > 1e-4:
                spread, c_, lo_ = max(spread, us_[i - 1] - lo_), c_ + 1, None
            lo_ = u if lo_ is None else lo_
            cid[(key, u)] = c_
        spread = max(spread, us_[-1] - lo_)
    out = ['', "The entries over the variants that record them, by the block they join — one event per width, starting grid, iteration and "
           "kernel whose duals agree within 10⁻⁴ — the same entrant at the same iteration, whether or not the variants' paths have parted "
           f"(the entries over the variants in brackets; the duals of one event agree to {spread:.1e}; referee rounds 65, 70) —: the entrant "
           "had the largest current rate |⟨rφ_j⟩|; its rate over "
           "the largest; its dual |u_j| at the block's start; whether it starts beyond the block's outermost node (how far) or between its "
           "nodes (their duals). The entrant is the kernel of the least time to the wall (1 − |u_j|)/|⟨rφ_j⟩| (`sliding_iss`'s rule), its "
           "dual charged during the earlier blocks.", '',
           '| joining | events (entries) | largest rate | rate / largest: median (range) | dual: median (least) | beyond: entrants, median '
           '(range) | between: entrants (duals) |', '|---|---|---|---|---|---|---|']
    counts = {}
    for name, xs_all in rows.items():
        if not xs_all:
            continue
        seen_, xs = set(), []
        for e, d, key in xs_all:
            k_ = (key, cid[(key, abs(e['u0']))])
            if k_ not in seen_:
                seen_.add(k_)
                xs.append((e, d))
        rel = np.array([abs(e['rate']) / e['rate_max'] for e, d in xs])
        u0 = np.array([abs(e['u0']) for e, d in xs])
        top = int(sum(abs(abs(e['rate']) - e['rate_max']) <= 1e-12 * e['rate_max'] for e, d in xs))
        out_ = [d for e, d in xs if d is not None and d > 0]
        ins = [abs(e['u0']) for e, d in xs if d is not None and d <= 0]
        counts[name] = len(xs)
        out.append(f"| {name} | {len(xs)} ({len(xs_all)}) | {top} | {np.median(rel):.2f} ({rel.min():.3f}–{rel.max():.2f}) | {np.median(u0):.3f} "
                   f"({u0.min():.3f}) | " + (f"{len(out_)}, {np.median(out_):.2f}σ ({min(out_):.2f}–{max(out_):.2f}σ)" if out_ else '—')
                   + ' | ' + (f"{len(ins)} ({min(ins):.3f}–{max(ins):.3f})" if ins else '—') + ' |')
    ov = sum(v for k, v in counts.items() if k != 'all') - counts['all']
    if ov:                                                             # one event, an optimal block in some variants, not in others
        out += ['', f"{ov} events join an optimal block in some variants and a non-optimal one in others: counted in both rows, once in "
                    "'all'."]
    return out + ['']


def c_job(args):
    s, q = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L, _, c0 = TM.setup()
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(TAU_C)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    xi = s * np.random.default_rng(5800 + 100 * LEVELS_C.index(s) + q).standard_normal(X_B.size)
    y = f + xi
    R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), polish=True, max_events=3 * K_C, tol=0.0, max_entries=K_C + 1,
                       max_doublings=MAX_DOUBLINGS, box=(X_B[0] - BOX_C * SIG, X_B[-1] + BOX_C * SIG),
                       w_scale=float(np.abs(y).max()), exit_kkt=True)
    rf = float(np.sqrt(np.mean(f ** 2)))
    bl = [b for b in blocks(R) if b['k'] > 0]
    ent = [e for e in R['events'] if e[2] == 'enter']
    taus = [e[0] for e in ent]
    out = dict(s=s, q=q, k=[b['k'] for b in bl], rms=[b['rms'] for b in bl],
               err=[float(np.sqrt(np.mean((K.phi(b['c']) @ b['w'] - f) ** 2)) / rf) for b in bl],
               at=[int(np.sum(np.abs(b['c'] - C_B) < 4 * ss)) for b in bl],
               gaps=[float((taus[i + 1] - taus[i]) / taus[i]) for i in range(len(taus) - 1)], exits=sum(e[2] == 'exit' for e in R['events']),
               stiff=R['stiff_redone'], runaway=R['runaway'], relaxed=[b['relaxed'] for b in bl], rf=rf, secs=time.time() - t0)
    return out


def part_C():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, [(s, q) for s in LEVELS_C for q in range(N_DRAW_C)], chunksize=1)


def d_job(args):
    """D: the tail ratio at (K, t), with enough terms for t^{2j} to fall below 10⁻¹⁷ (40/ln(1/t) of them, at least 40)."""
    K, t = args
    n = max(40, int(np.ceil(40 / np.log(1 / t))))
    return dict(K=K, t=t, ratio=W.tail_ratio(K, t, n_terms=n), n_terms=n)


def lead_coeffs(j_max):
    """The ratios P_j(K) = D_{2K+2j,K}/D_{2K,K} as polynomials in K: the j-th divided difference over K = 1 … j + 1 (the leading
    coefficient, exact) and the (j + 1)-th over K = 1 … j + 2 (zero for degree j)."""
    from fractions import Fraction

    def dd(ps):
        return ps[0][1] if len(ps) == 1 else (dd(ps[1:]) - dd(ps[:-1])) / (ps[-1][0] - ps[0][0])
    out = []
    for j in range(1, j_max + 1):
        pts = [(K, Fraction(W.gh_defects(K, 2 * K + 2 * j)[2 * K + 2 * j], W.gh_defects(K, 2 * K)[2 * K])) for K in range(1, j + 3)]
        lead = dd(pts[:j + 1])
        out.append(dict(j=j, lead=(lead.numerator, lead.denominator), next_dd=float(dd(pts))))
    return out


def part_D():
    from multiprocessing import Pool
    ts = tuple(float(W.params(SIG, sigma_star(r_))['t']) for r_ in (1.0, 2.0)) + (0.5, 0.8)
    with Pool(N_WORKERS) as pool:
        rows = pool.map(d_job, [(K, t) for t in ts for K in KS_D], chunksize=1)
    return dict(rows=rows, lead=lead_coeffs(J_MAX_D))


def flank(r_, c, w):
    """B: the residual's correlation with a kernel at c, C(c) = ⟨r, φ_σ(· − c)⟩ on the samples, after a block (c, w) on the feature of
    τ/σ = r_, on a fine grid of c: min C/max C (zero or positive where C ≥ 0) and its largest values' distances beyond the outermost
    nodes, left and right (in σ)."""
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    r = f - np.exp(-(X_B[:, None] - np.asarray(c)[None, :]) ** 2 / (2 * SIG ** 2)) @ np.asarray(w)
    cg = np.linspace(X_B[0], X_B[-1], 3201)
    C = np.exp(-(X_B[:, None] - cg[None, :]) ** 2 / (2 * SIG ** 2)).T @ r / X_B.size
    lo, hi = float(np.min(c)), float(np.max(c))
    L, R = cg < lo, cg > hi
    return float(C.min() / C.max()), float((lo - cg[L][np.argmax(C[L])]) / SIG), float((cg[R][np.argmax(C[R])] - hi) / SIG)


def predict_stop(s, rf, N):
    """(i)'s prediction of the discrepancy stop: the first K with e_K²rf² ≤ s²((1 + δ)² − 1 + 2K/N), and the error there."""
    ss = sigma_star(TAU_C)
    for K in range(1, 40):
        e = W.l2_error(*W.half_squeeze(K, SIG, ss), SIG, ss)
        if (e * rf) ** 2 <= s ** 2 * ((1 + DELTA) ** 2 - 1 + 2 * K / N):
            return K, float(np.sqrt(e ** 2 + 2 * K * s ** 2 / (N * rf ** 2)))
    return -1, np.nan


def stop_law(s, rf, N, k_max=12):
    """(iv)'s stop as a margin (referee round 54). The rms of the K-kernel fit to y = f* + ξ has the expectation
    m_K = s(b_K² + 1 − 2K/N)^{1/2}, b_K = e_K‖f*‖_N/s (2K parameters), and the scatter sd(rms²)/s² = (2(N − 2K)/N² + 4b_K²/N)^{1/2} (the
    noise's chi-square and its cross term with the bias); the rms falls along the path, so the rule stops at K with probability
    P(rms_K ≤ (1 + δ)s) − P(rms_{K−1} ≤ (1 + δ)s), in the normal approximation; e_K the half-squeeze rule's L² error, as in
    `predict_stop`. Also the first-term formula's K(s): the root of K = log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/(s(2δ + δ²)^{1/2}))/log(1/t)
    → (K(s), {K: P(stop at K)})."""
    from scipy.stats import norm
    ss = sigma_star(TAU_C)
    t = W.params(SIG, ss)['t']
    Kr = 1.0
    for _ in range(100):
        Kr = np.log((1 - t * t) ** 0.25 * (np.pi * Kr) ** 0.25 * rf / (s * np.sqrt(2 * DELTA + DELTA ** 2))) / np.log(1 / t)
    cdf = [0.0]
    for K in range(1, k_max + 1):
        b2 = (W.l2_error(*W.half_squeeze(K, SIG, ss), SIG, ss) * rf / s) ** 2
        m2 = b2 + 1 - 2 * K / N
        sd2 = np.sqrt(2 * (N - 2 * K) / N ** 2 + 4 * b2 / N)
        cdf.append(float(norm.cdf(((1 + DELTA) ** 2 - m2) / sd2)))
    return float(Kr), {K: cdf[K] - cdf[K - 1] for K in range(1, k_max + 1)}


def report(A=None, B=None, C=None, D=None, B2=None, A5=None, B3=None, B4=None):
    L_ = ['# Proposition D‴.58: the transit under width mismatch (`exp_width_mismatch.py`)', '',
          f"Kernels of width σ = {SIG}; a feature Aφ_σ*(· − c*), σ* = σ(1 + (τ/σ)²)^{{1/2}}; t = (σ* − σ)/(σ* + σ), λ = τ²/(σ² + σ*²) "
          "(`src/core/wide_feature.py`).", '']
    if A is not None:
        L_ += ['## A. The rules on ℝ', '',
               "Relative L² errors: the half-squeeze rule (nodes τ(σ/σ*)^{1/2}z_k, amplitudes A(σ*/σ)^{1/2}w_k exp(y_k²/(2σ(σ + σ*)))) by "
               "the trapezoid rule and by the Mehler series (their relative difference), the series' first term's share of the squared "
               "error; moment matching (nodes τz_k, amplitudes A(σ*/σ)w_k); the half-squeeze nodes with least-squares amplitudes; the "
               f"optimal K-kernel fit (variable projection; not computed below {E_FLOOR:g}, where its solve stops short of the optimum — A3 marks "
               "the optimum unconverged where the one-parameter rule beats it). Rows stop "
               f"where the series falls below {G_FLOOR:g} (the trapezoid rule's floor). Ratios per kernel e_{{K−1}}/e_K; the optimum's "
               "nodes over the half-squeeze nodes (range over the nodes).", '']
        for r_ in TAUS_A:
            rows = [x for x in A['rows'] if x['r'] == r_]
            t, lam = rows[0]['t'], rows[0]['lam']
            L_ += [f"### τ/σ = {r_:g} (1/t = {1 / t:.4f}, 1/λ = {1 / lam:.4f})", '',
                   '| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | '
                   'ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |', '|---|---|---|---|---|---|---|---|---|']
            for i, x in enumerate(rows):
                rat = (f"{rows[i - 1]['e_hs'] / x['e_hs']:.4f}, " + (f"{rows[i - 1]['e_opt'] / x['e_opt']:.4f}"
                       if 'e_opt' in x and 'e_opt' in rows[i - 1] else '—') + f", {rows[i - 1]['e_mm'] / x['e_mm']:.4f}") if i else '—'
                nod = '—'
                if 'y_opt' in x and x['K'] > 1:
                    m = np.abs(x['y_hs']) > 1e-9
                    q = np.sort(x['y_opt'])[m] / np.sort(x['y_hs'])[m]
                    nod = f"{q.min():.4f} … {q.max():.4f}"
                L_.append(f"| {x['K']} | {x['e_hs']:.4e} | {x['e_series']:.4e} ({abs(x['e_hs'] / x['e_series'] - 1):.0e}) | "
                          f"{(x['e_first'] / x['e_series']) ** 2:.5f} | {x['e_mm']:.3e} | {x['e_ls']:.4e} | "
                          + (f"{x['e_opt']:.4e}" if 'e_opt' in x else '—') + f" | {rat} | {nod} |")
            L_.append('')
        L_ += ['### A2. The series\' tail', '',
               "D_{2K+2j,K}/D_{2K,K} (integers) for j = 1, …, 4 and K = 1, …, " + f"{K_A2}" + "; their j-th differences in K, "
               "constant if they are polynomials of degree j, against j!·binom(1/2, j)(−2)^j = −1, −1, −3, −15.", '',
               '| j | K = 1, 2, 3, 4, 5, …, ' + f"{K_A2}" + ' | j-th differences (all K) | j!·binom(1/2, j)(−2)^j |', '|---|---|---|---|']
        ks = sorted(A['ratios'])
        for j in range(1, 5):
            v = np.array([A['ratios'][K][j - 1] for K in ks], float)
            dif = np.diff(v, n=j)
            lead = factorial(j) * comb(1, 1) * float(np.prod([0.5 - i for i in range(j)]) / factorial(j)) * (-2) ** j
            L_.append(f"| {j} | {', '.join(f'{int(x)}' for x in v[:5])}, …, {int(v[-1])} | "
                      f"{', '.join(sorted({f'{d:g}' for d in dif}))} | {lead:g} |")
        if 'e_fock' in A['rows'][0]:
            L_ += ['', '### A3. Scaled rules', '',
                   "The Gauss–Hermite rule of N(0, v) for the reweighted mixing law (`scaled_rule`) at the half-squeeze variance v_hs = τ²σ/σ*, "
                   "the natural Fock scaling v_hs/cosh r and the Mercer-zero scaling τ²(1 + 2τ²/σ²)^{−1/2}: their ratios per kernel at the last "
                   "K against 1/t, (2 − λ)/λ = (1 − t + t²)/t and 1/q, q = (β² − 1)/(β² + 1), β² = (1 + 2τ²/σ²)^{1/2}; the variance of least "
                   "error (`best_scale`): v/v_hs, its error and the optimal fit's over the rule's, at K = 1 and at the last K with an optimum "
                   "(in brackets the share of the gap to the optimum that the variance recovers, (e_hs − e_best)/(e_hs − min(e_opt, e_best)); "
                   f"'optimum unconverged' where e_best < e_opt by more than a relative {UNCONV:g}); over all K with an optimum the shares "
                   + gap_range(A) + ":", '',
                   '| τ/σ | half-squeeze (1/t) | Fock ((2 − λ)/λ) | Mercer-zero (1/q) | K = 1: v/v_hs, e_best/e_hs, e_opt/e_hs | the last K: K, '
                   'v/v_hs, e_best/e_hs, e_opt/e_hs |', '|---|---|---|---|---|---|']
            for r_ in TAUS_A:
                rows = [x for x in A['rows'] if x['r'] == r_]
                t, lam = rows[0]['t'], rows[0]['lam']
                b2 = np.sqrt(1 + 2 * r_ ** 2)
                q = (b2 - 1) / (b2 + 1)
                a, b = rows[-2], rows[-1]
                opt = [x for x in rows if 'e_opt' in x]
                cell = lambda x: (f"{x['v_best']:.5f}, {x['e_best'] / x['e_hs']:.5f}, " + (f"{x['e_opt'] / x['e_hs']:.5f}" if 'e_opt' in x else '—')
                                  + (f" ({gap_share(x)})" if 'e_opt' in x else ''))
                L_.append(f"| {r_:g} | {a['e_hs'] / b['e_hs']:.3f} ({1 / t:.3f}) | {a['e_fock'] / b['e_fock']:.3f} ({(2 - lam) / lam:.3f}) | "
                          f"{a['e_mercer'] / b['e_mercer']:.3f} ({1 / q:.3f}) | {cell(rows[0])} | "
                          + (f"{opt[-1]['K']}, {cell(opt[-1])}" if opt else '—') + ' |')
        if 'narrow' in A:
            L_ += ['', '### A4. Narrower features: the continued frame', '',
                   "σ* < σ: the nodes imaginary in conjugate pairs, the K-term sum e^{−x²/(2σ²)}Σ_k b_k cos(η_kx/σ²) at the Hermite-zero "
                   "frequencies η_k = (2σ²|sinh r|)^{1/2}z_k (`narrow_rule`); its relative L² error (trapezoid) against the Mehler series at "
                   "|t| = (σ − σ*)/(σ + σ*), and the ratios per kernel against 1/|t|:", '',
                   '| σ*/σ | 1/|t| | K: error (relative difference from the series) | ratio per kernel at the last K |', '|---|---|---|---|']
            for q in NARROW_A:
                rs = [x for x in A['narrow'] if x['q'] == q]
                L_.append(f"| {q:g} | {1 / rs[0]['t']:.4f} | " + ', '.join(f"{x['K']}: {x['e']:.3e} ({abs(x['e'] / x['e_series'] - 1):.0e})"
                                                                         for x in rs) + f" | {rs[-2]['e'] / rs[-1]['e']:.4f} |")
        L_ += a5_lines(A5)
        from scipy.special import hyp2f1
        L_ += ['', "The tail's share Σ_{n>2K}/(first term) against ₂F₁(−½, −½; 1; t²) − 1 (the share's limit if the ratios' leading "
               "coefficients govern the tail):", '', '| τ/σ | t | K = 5 | 10 | 20 | 30 | ₂F₁(−½, −½; 1; t²) − 1 |', '|---|---|---|---|---|---|---|']
        for r_ in TAUS_A:
            sh = [x['share'] for x in A['share'] if x['r'] == r_]
            t = W.params(SIG, sigma_star(r_))['t']
            L_.append(f"| {r_:g} | {t:.5f} | " + ' | '.join(f"{v:.5f}" for v in sh) + f" | {hyp2f1(-0.5, -0.5, 1, t ** 2) - 1:.5f} |")
    if B is not None:
        L_ += ['', '## B. The transit', '',
               f"The sliding path (polished; a kernel leaves only if its correlation stays inward with its own contribution removed) on "
               f"f* = {A_B}φ_σ*(· − {C_B}) sampled at {X_B.size} points of [{X_B[0]:g}, {X_B[-1]:g}], 21 kernels started at the grid "
               "−1, −0.9, …, 1, to " + f"{K_B + 1} entries, the fast flow's chunks up to 10⁵·2^{MAX_DOUBLINGS}: after each entry the relaxed "
               "block's relative rms error against the optimal fit with as many kernels on the same samples (the better of variable "
               f"projection from the block and from the half-squeeze nodes: the best local optimum found); a block is optimal within "
               f"{OPT_TOL:g}. The numerical variants: " + ', '.join(f"{n} (polish step {st:g}σ, grid shifted by {sh:g}, the polish's budget "
                                                                     f"{nf:g} evaluations" + (f", repeated up to {rp} times" if rp else '')
                                                                     + (f", its gradient test {gt:g}" if gt is not None else ', no gradient test')
                                                                     + ')' for n, st, sh, nf, rp, gt in VARIANTS_B)
               + ". 'Relaxed': the block passed the stationarity test. The first non-optimal block's polishes as status/evaluations/accepted "
               f"(status 0: the budget was spent first; 1: the gradient test), a repeated polish's end over the optimum in brackets. 'none': every block optimal, to the path's cap of {K_B} kernels "
               "(the range censored).", '',
               '| τ/σ | variant | optimal blocks (kernels) | the first non-optimal: kernels, its optimum\'s error, block / optimum, relaxed | its '
               'polishes | '
               'relaxed blocks | chunks redone (stall, failure), failed | polish rejected | exits |', '|---|---|---|---|---|---|---|---|---|']
        for r_ in TAUS_B:
            for v in VARIANTS_B:
                b = next((x for x in B if x['r'] == r_ and x['var'] == v[0]), None)
                if b is None:
                    continue
                opt = [o['k'] for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
                bad = next((o for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL), None)
                L_.append(f"| {r_:g} | {v[0]} | {', '.join(map(str, opt)) or '—'} | "
                          + (f"{bad['k']}, {bad['e_opt']:.2e}, {bad['e'] / bad['e_opt']:.4g}, {bad['relaxed']}" if bad else 'none') + f" | "
                          + (polish_cell(bad, float(np.sqrt(np.mean((A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))) ** 2))))
                             if bad else '—') + f" | "
                          f"{sum(o['relaxed'] for o in b['blocks'])}/{len(b['blocks'])} | {b['stiff']}, {b['redone']}, {b['failed']} | "
                          f"{b['rejected']} | {sum(e[2] == 'exit' for e in b['events'])} |")
        ok = [o['e_opt'] for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
        no = [o for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL]
        vac = [o['e'] for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL and o['relaxed']
               and cs_vacuous(b['r'], b['var'], o)]
        mg = [cs_margin(b['r'], b['var'], o) for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL and o['relaxed']]
        mg = [m for m in mg if m < 1]
        L_ += ['', f"Over every width and variant: the optimal blocks' optimum errors down to {min(ok):.2e}; the non-optimal blocks' "
               + (f"from {max(o['e_opt'] for o in no):.2e} down ({sum(not o['relaxed'] for o in no)} of {len(no)} failed the stationarity "
                  "test" + (f"; {len(vac)} of those that passed it pass it at any state by Cauchy–Schwarz, the test being absolute — their "
                            f"errors up to {max(vac):.1e}" if vac else '')
                  + (f"; for the other {len(mg)} its threshold over the Cauchy–Schwarz bound — the residual's cosine with the block's "
                     f"directions it allows — is {min(mg):.1e}–{max(mg):.2f} (median {np.median(mg):.1e}; above 0.1 in "
                     f"{sum(m > 0.1 for m in mg)}): a weak test" if mg else '') + "); "
                  f"non-optimal blocks with an optimum above 10⁻⁷: {sum(o['e_opt'] > 1e-7 for o in no)}; the non-optimal blocks' errors "
                  f"{min(o['e'] for o in no):.2e} … {max(o['e'] for o in no):.2e}." if no else "none."), '']
        L_ += b2_line(B2)
        L_ += b3_lines(B3, B4, B)
        L_ += b4_lines(B4)
        L_ += entry_lines(B)
        for b in [x for x in B if x['var'] == VARIANTS_B[0][0]]:
            L_ += [f"### τ/σ = {b['r']:g}, {b['var']} (events: " + ', '.join(f"{j}{'+' if k == 'enter' else '−'}" for t_, j, k, c in b['events'])
                   + ")", '',
                   '| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | '
                   'half-squeeze on ℝ | ratio per kernel |', '|---|---|---|---|---|---|---|---|---|---|']
            for i, o in enumerate(b['blocks']):
                rat = f"{b['blocks'][i - 1]['e'] / o['e']:.3f}" if i and b['blocks'][i - 1]['k'] == o['k'] - 1 else '—'
                L_.append(f"| {o['k']} | {o['tau']:.4g} | {o['e']:.4e} | {o['e_opt']:.4e} | {o['e'] / o['e_opt']:.6f} | {o['relaxed']} | "
                          f"{o['dc'] / SIG:.1e} | {o['dc_hs'] / SIG:.3f} | {o['e_hs']:.4e} | {rat} |")
            opt = [o for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
            fl = [flank(b['r'], o['c'], o['w']) for o in opt]
            if fl:
                L_ += ['', "The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = "
                       + ', '.join(str(o['k']) for o in opt) + f"): min C/max C within ±{max(abs(x[0]) for x in fl):.0e}, its peaks beyond "
                       "the outermost nodes at (left, right) " + ', '.join(f"({x[1]:.2f}, {x[2]:.2f})σ" for x in fl) + '.', '']
    if C is not None:
        N = X_B.size
        L_ += ['## C. Under noise', '',
               f"τ/σ = {TAU_C:g}, {N_DRAW_C} draws per level, the path to {K_C + 1} entries: the discrepancy rule's stop (kernels), its "
               "prediction from (i) and the predicted error there, the stopped fit's relative error (median [range]), the kernels within "
               "4σ* of the feature at the stop, the largest relative gap's kernel count (in how many draws it is the stop's).", '',
               '| s | stop: kernels (median [range]) | predicted | error at the stop: median [range] | predicted | kernels at the feature '
               '| largest gap at |', '|---|---|---|---|---|---|---|']
        L2 = ['', "At the stop (K kernels): the error over the oracle's s(2K/N)^{1/2} (the variance alone, the right support), over the "
              "bias–variance prediction (e_K²‖f*‖²_N + 2Ks²/N)^{1/2} — e_K the noiseless optimum's error with K kernels (B's base path "
              f"at τ/σ = {TAU_C:g}) — and over the discrepancy rule's bound s((1 + δ)² − 1 + 4K/N)^{{1/2}}; paths with a runaway kernel "
              "before the stop:", '', '| s | error / oracle: median [range] | / bias–variance | / bound: max | runaway |', '|---|---|---|---|---|']
        L3 = ['', "The stop as a margin (`stop_law`): the first-term formula's root K*(s) of K = log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/"
              "(s(2δ + δ²)^{1/2}))/log(1/t), whose ceiling is the predicted stop, and the probability of the stop at each K — the rms's "
              "expectation s(b_K² + 1 − 2K/N)^{1/2} and its scatter from the noise's chi-square and its cross term with the bias, "
              "b_K = e_K‖f*‖_N/s, e_K the half-squeeze rule's error, in the normal approximation (the rms falls along the path) —; the "
              "expected number of the draws stopping at the predicted K against the observed:", '',
              '| s | K*(s) | P(stop at K) | at the predicted K | expected | observed |', '|---|---|---|---|---|---|']
        for s in LEVELS_C:
            ds = [d for d in C if d['s'] == s]
            if not ds:
                continue
            ks, es, at, kg, same = [], [], [], [], 0
            for d in ds:
                i = next((i for i, r_ in enumerate(d['rms']) if r_ <= (1 + DELTA) * s), None)
                if i is None:
                    continue
                ks.append(d['k'][i]); es.append(d['err'][i]); at.append(d['at'][i])
                if d['gaps']:
                    kg.append(int(np.argmax(d['gaps'])) + 1)
                    same += kg[-1] == d['k'][i]
            Kp, ep = predict_stop(s, ds[0]['rf'], N)
            fmt = lambda v, f_: f"{np.median(v):{f_}} [{np.min(v):{f_}}, {np.max(v):{f_}}]" if len(v) else '—'
            L_.append(f"| {s:g} | {fmt(ks, 'g')} ({len(ks)}/{len(ds)} stopped) | {Kp} | {fmt(es, '.3e')} | {ep:.3e} | {fmt(at, 'g')} | "
                      f"{fmt(kg, 'g')}" + (f" (the stop in {same})" if kg else '') + " |")
            Kr, P = stop_law(s, ds[0]['rf'], N)
            L3.append(f"| {s:g} | {Kr:.2f} | " + ', '.join(f"{K}: {p:.3f}" for K, p in P.items() if p >= 5e-4) + f" | {P.get(Kp, 0.0):.3f} | "
                      f"{len(ds) * P.get(Kp, 0.0):.1f} of {len(ds)} | {sum(k == Kp for k in ks)} |")
            rf = ds[0]['rf']
            eopt = ({o['k']: o['e_opt'] for b in B for o in b['blocks'] if b['r'] == TAU_C and b['var'] == VARIANTS_B[0][0]}
                    if B is not None else {})
            r1 = [e / (s * np.sqrt(2 * k / N) / rf) for k, e in zip(ks, es)]
            r2 = [e / np.sqrt(eopt[k] ** 2 + 2 * k * s ** 2 / (N * rf ** 2)) for k, e in zip(ks, es) if k in eopt]
            r3 = [e / (s * np.sqrt((1 + DELTA) ** 2 - 1 + 4 * k / N) / rf) for k, e in zip(ks, es)]
            L2.append(f"| {s:g} | {fmt(r1, '.2f')} | {fmt(r2, '.3f')} | {max(r3):.3f} | {sum(d.get('runaway', False) for d in ds)} |"
                      if ks else f"| {s:g} | — | — | — | — |")
        L_ += L2 + L3
    if D is not None:
        from scipy.special import hyp2f1
        ts = sorted({r['t'] for r in D['rows']})
        Ks = sorted({r['K'] for r in D['rows']})
        by = {(r['t'], r['K']): r for r in D['rows']}
        L_ += ['', '## D. The tail of the exact series (the sharp constant, conjectural)', '',
               "The half-squeeze rule's squared error over the series' first term, Σ_j t^{2j}(2K)!/(2K + 2j)!·(D_{2K+2j,K}/D_{2K,K})² (the "
               "defects exact, by the remainders modulo He_K), against ₂F₁(−½, −½; 1; t²), the sum of the terms' limits binom(1/2, j)²t^{2j}; "
               "K times the gap in parentheses:", '',
               '| t | ₂F₁ | ' + ' | '.join(f'K = {K}' for K in Ks) + ' |', '|---|---|' + '---|' * len(Ks)]
        for t in ts:
            F = float(hyp2f1(-0.5, -0.5, 1.0, t * t))
            L_.append(f"| {t:.4g} | {F:.6f} | " + ' | '.join(f"{by[(t, K)]['ratio']:.6f} ({K * (by[(t, K)]['ratio'] - F):+.3f})"
                                                           for K in Ks) + ' |')
        cat = lambda j: comb(2 * (j - 1), j - 1) // j
        L_ += ['', 'The ratios D_{2K+2j,K}/D_{2K,K} are polynomials in K of degree j (the (j + 1)-th divided difference over K = 1 … j + 2: '
               + ', '.join(f"{r['next_dd']:g}" for r in D['lead']) + ') with the leading coefficients '
               + ', '.join(f"{r['lead'][0]}/{r['lead'][1]}" if r['lead'][1] != 1 else f"{r['lead'][0]}" for r in D['lead']) + ' at j = 1 … ' + str(len(D['lead'])) + ' — −C_{j−1}/2^{j−1} '
               + ('at every j' if all(r['lead'][0] * 2 ** (r['j'] - 1) == -cat(r['j']) * r['lead'][1] for r in D['lead']) else 'NOT at every j')
               + ' (C the Catalan numbers), binom(1/2, j)(−2)^j.']
    L_ += ['', 'Numbers are produced by `experiments/exp_width_mismatch.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'width_mismatch_raw.npz')
    md = os.path.join(OUT, 'width_mismatch.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'B2', 'A5', 'B3', 'B4')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    idx = ({int(i) for i in sys.argv[sys.argv.index('--variants') + 1].split(',')} if '--variants' in sys.argv else None)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', lambda: run_B(res.get('B'), idx)), ('C', part_C), ('D', part_D), ('B2', lambda: part_B2(res['B'])),
                  ('A5', lambda: part_A5(res['A'])), ('B3', lambda: part_B3(res['B'])), ('B4', lambda: part_B4(res['B'], res['B3']))):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
