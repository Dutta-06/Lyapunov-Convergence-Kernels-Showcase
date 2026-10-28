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
B. The transit: the sliding path (polished; the KKT exit rule) on f* (A_B, c* = C_B off the kernels' grid) sampled on X_B, in five
   numerical variants (VARIANTS_B: the polish step, the starting grid's offset, the polish's budget of evaluations and its repetition
   after the stationarity test; referee rounds 48 and 54): after each entry the relaxed block's relative error against the optimal fit
   with as many kernels on the same samples (variable projection from the block and from the half-squeeze rule, the better), whether
   the relaxation passed the stationarity test, the first non-optimal block's polishes (status, evaluations, accepted), the block's
   nodes against the optimum's, the continuous half-squeeze error; the entries (the entrant's rate against the largest, its dual), the
   chunks redone after a stall or a failed solve, the polish's rejections.
C. Under noise: y = f* + ξ (i.i.d. N(0, s²) on X_B), the path to K_C + 1 entries (its centres in the window ± BOX_C·σ, the thresholds
   on max|y|, the KKT exit rule): the discrepancy rule's stop (the first block with rms ≤ (1 + DELTA)s) against its prediction from (i)
   — the first K with e_K²‖f*‖²_rms ≤ s²((1 + δ)² − 1 + 2K/N), e_K the half-squeeze rule's L² error — and the probability of each stop
   from the rms's expectation and scatter (`stop_law`), the stopped fit's error against (e_K²‖f*‖²_rms + 2Ks²/N)^{1/2} (e_K B's optima
   on the samples); the kernels at the feature and elsewhere; the largest relative gap's K.

Usage: python experiments/exp_width_mismatch.py [--only A,B,C] [--from-raw] → results/chaos/width_mismatch.md (+ _raw.npz).
WM_WORKERS (default 2).
"""
import os
import sys
import time
import warnings
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
VARIANTS_B = (('base', 1.0, 0.0, 5000, 0), ('polish 2σ', 2.0, 0.0, 5000, 0), ('grid +0.05', 1.0, 0.05, 5000, 0),
              ('polish 2σ, 5·10⁴ evaluations', 2.0, 0.0, 50_000, 0), ('polish 2σ, 5·10⁴, repeated', 2.0, 0.0, 50_000, 3))
# B: (name, the polish step, the starting grid's shift, the polish's budget of evaluations, its repetitions after the stationarity test —
# the last two variants: referee round 54)
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


def b_job(args):
    r_, var = args
    name, step, shift, nfev, rep_ = next(v for v in VARIANTS_B if v[0] == var)
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L, _, c0 = TM.setup()
    c0 = c0 + shift
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, max_events=3 * K_B, tol=0.0, max_entries=K_B + 1,
                       max_doublings=MAX_DOUBLINGS, polish_step=step, exit_kkt=True, polish_nfev=nfev, repolish=rep_)
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


def part_B():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(r_, v[0]) for v in VARIANTS_B for r_ in TAUS_B], chunksize=1)


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


def polish_cell(o):
    """A block's polishes: status/evaluations/accepted (status 0: the budget spent first; referee round 54)."""
    return ' '.join(f"{p['status']}/{p['nfev']}/{'a' if p['accepted'] else 'r'}" for p in o.get('polish', [])) or '—'


def entry_lines(B):
    """B's entries (`sliding_iss`'s entry_info): the entrant's rate against the largest inactive rate, its dual before the entry, and
    whether the entrant had the largest rate — it is the kernel of the least time to the wall (1 − |u_j|)/|⟨rφ_j⟩| (referee round 54);
    where the entrant starts against the block it joins (block it − 1, the entries in the order of the events): beyond its outermost
    node (how far) or between its nodes (its dual)."""
    es, pos = [], []
    for b in B:
        info = b.get('entries') if b.get('entries') is not None else []
        for e, ev in zip(info, [ev for ev in b['events'] if ev[2] == 'enter']):
            es.append(e)
            if 1 <= e['it'] <= len(b['blocks']):
                nodes = np.asarray(b['blocks'][e['it'] - 1]['c'])
                pos.append((max(nodes.min() - ev[3], ev[3] - nodes.max()) / SIG, abs(e['u0'])))
    if not es:
        return []
    rel = np.array([abs(e['rate']) / e['rate_max'] for e in es])
    u0 = np.array([abs(e['u0']) for e in es])
    top = int(sum(abs(abs(e['rate']) - e['rate_max']) <= 1e-12 * e['rate_max'] for e in es))
    out = (f"The entries ({len(es)} over the variants that record them): the entrant had the largest current rate |⟨rφ_j⟩| in "
           f"{top}; its rate over the largest a median {np.median(rel):.2f} (range {rel.min():.3f}–{rel.max():.2f}), its dual |u_j| at the "
           f"block's start a median {np.median(u0):.3f} (least {u0.min():.3f}): the entrant is the kernel of the least time to the wall, "
           "(1 − |u_j|)/|⟨rφ_j⟩|, its dual charged during the earlier blocks.")
    out_ = [d for d, u in pos if d > 0]
    ins = [u for d, u in pos if d <= 0]
    if pos:
        out += (f" Of the {len(pos)} entrants joining a block, {len(out_)} start beyond its outermost node"
                + (f" (a median {np.median(out_):.2f}σ, range {min(out_):.2f}–{max(out_):.2f}σ, beyond it)" if out_ else '')
                + f" and {len(ins)} between its nodes" + (f" (their duals at the block's start {', '.join(f'{u:.3f}' for u in ins)})" if ins
                                                          else '') + '.')
    return ['', out, '']


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


def report(A=None, B=None, C=None, D=None):
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
                                                                     + ')' for n, st, sh, nf, rp in VARIANTS_B)
               + ". 'Relaxed': the block passed the stationarity test. The first non-optimal block's polishes as status/evaluations/accepted "
               "(status 0: the budget was spent first).", '',
               '| τ/σ | variant | optimal blocks (kernels) | the first non-optimal: kernels, its optimum\'s error, relaxed | its polishes | '
               'relaxed blocks | chunks redone (stall, failure), failed | polish rejected | exits |', '|---|---|---|---|---|---|---|---|---|']
        for r_ in TAUS_B:
            for v in VARIANTS_B:
                b = next((x for x in B if x['r'] == r_ and x['var'] == v[0]), None)
                if b is None:
                    continue
                opt = [o['k'] for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
                bad = next((o for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL), None)
                L_.append(f"| {r_:g} | {v[0]} | {', '.join(map(str, opt)) or '—'} | "
                          + (f"{bad['k']}, {bad['e_opt']:.2e}, {bad['relaxed']}" if bad else 'none') + f" | "
                          + (polish_cell(bad) if bad else '—') + f" | "
                          f"{sum(o['relaxed'] for o in b['blocks'])}/{len(b['blocks'])} | {b['stiff']}, {b['redone']}, {b['failed']} | "
                          f"{b['rejected']} | {sum(e[2] == 'exit' for e in b['events'])} |")
        ok = [o['e_opt'] for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
        no = [o for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL]
        L_ += ['', f"Over every width and variant: the optimal blocks' optimum errors down to {min(ok):.2e}; the non-optimal blocks' "
               + (f"from {max(o['e_opt'] for o in no):.2e} down ({sum(not o['relaxed'] for o in no)} of {len(no)} failed the stationarity "
                  "test); "
                  f"non-optimal blocks with an optimum above 10⁻⁷: {sum(o['e_opt'] > 1e-7 for o in no)}; the non-optimal blocks' errors "
                  f"{min(o['e'] for o in no):.2e} … {max(o['e'] for o in no):.2e}." if no else "none."), '']
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
        L3 = ['', "The stop as a margin (`stop_law`): the first-term formula's K(s) (the root of K = log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/"
              "(s(2δ + δ²)^{1/2}))/log(1/t), whose ceiling is the predicted stop), and the probability of the stop at each K — the rms's "
              "expectation s(b_K² + 1 − 2K/N)^{1/2} and its scatter from the noise's chi-square and its cross term with the bias, "
              "b_K = e_K‖f*‖_N/s, e_K the half-squeeze rule's error, in the normal approximation (the rms falls along the path) —; the "
              "expected number of the draws stopping at the predicted K against the observed:", '',
              '| s | K(s) | P(stop at K) | at the predicted K | expected | observed |', '|---|---|---|---|---|---|']
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
    keys = ('A', 'B', 'C', 'D')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
