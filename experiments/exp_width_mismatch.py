"""
Proposition D‴.58: the transit under width mismatch — a feature wider than the kernels is formed as a sequence of Gaussian quadratures,
the error falling by t = (σ* − σ)/(σ* + σ) per kernel
====================================================================================================================================

Kernels of width σ = 0.12 (the program's), a feature f* = Aφ_{σ*}(· − c*) with σ* = σ(1 + (τ/σ)²)^{1/2} (`src/core/wide_feature.py`):

A. The rules on ℝ: for τ/σ in TAUS_A and K ≤ K_A, the half-squeeze rule's L² error (trapezoid on a fine grid) against the exact Mehler
   series and its first term; the Gauss–Hermite rule of the mixing law (moment matching); the half-squeeze nodes with least-squares
   amplitudes; the optimal K-kernel fit (variable projection from the half-squeeze rule) — the ratios per kernel against 1/t and 1/λ;
   the optimum's nodes against the half-squeeze nodes. A2: the series' tail — the ratios D_{2K+2j,K}/D_{2K,K} for K ≤ K_A2 (integers;
   their j-th differences in K against j!·binom(1/2, j)(−2)^j) and the tail's share against ₂F₁(−½, −½; 1; t²) − 1.
B. The transit: the sliding path (polished) on f* (A_B, c* = C_B off the kernels' grid) sampled on X_B: after each entry the relaxed
   block's relative error against the optimal fit with as many kernels on the same samples (variable projection from the block and
   from the half-squeeze rule, the better), the block's nodes against the optimum's, the continuous half-squeeze error; the entries.
C. Under noise: y = f* + ξ (i.i.d. N(0, s²) on X_B), the path to K_C entries: the discrepancy rule's stop (the first relaxed block
   with rms ≤ (1 + DELTA)s) against its prediction from (i) — the first K with e_K²‖f*‖²_rms ≤ s²((1 + δ)² − 1 + 2K/N) —, the stopped
   fit's error against (e_K²‖f*‖²_rms + 2Ks²/N)^{1/2}; the kernels at the feature and elsewhere; the largest relative gap's K.

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
E_FLOOR, G_FLOOR = 1e-10, 1e-13                                      # below them double precision cannot place the optimum, the grid's error
TAUS_B, A_B, C_B, K_B = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0), 0.8, 0.03, 12
X_B = np.linspace(-1.6, 1.6, 321)
LEVELS_C, N_DRAW_C, TAU_C, K_C, DELTA = (1e-4, 1e-3, 1e-2), 10, 1.0, 10, 0.1
MAX_DOUBLINGS = 24                                                    # the fast flow to 10⁵·2²⁴ ≈ 1.7·10¹²: a block's slowest relaxations


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
            row = dict(r=r_, K=K, t=P['t'], lam=P['lam'], e_hs=e, e_series=float(np.sqrt(s2)), e_first=float(np.sqrt(f2)),
                       e_mm=W.l2_error(ym, am, SIG, ss), e_ls=W.l2_error(y, W.ls_amplitudes(y, SIG, ss), SIG, ss), y_hs=y)
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
    return dict(rows=rows, ratios=ratios, share=share)


def blocks(R):
    """The relaxed blocks along a path: (k kernels, τ, centres, coefficients, rms) after each event."""
    out = []
    for h in R['hist'][1:]:
        act = np.flatnonzero(h[2] != 0)
        o = np.argsort(h[1][act])
        out.append(dict(k=int(act.size), tau=float(h[0]), c=h[1][act][o], w=h[2][act][o], rms=float(h[3])))
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


def b_job(r_):
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L, _, c0 = TM.setup()
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, max_events=3 * K_B, tol=0.0, max_entries=K_B + 1,
                       max_doublings=MAX_DOUBLINGS)
    rf = float(np.sqrt(np.mean(f ** 2)))
    out = []
    for b in blocks(R):
        if b['k'] == 0:
            continue
        yh, _ = W.half_squeeze(b['k'], SIG, ss)
        cands = [sample_fit(K.phi, f, b['c']), sample_fit(K.phi, f, C_B + yh)]
        c_o, a_o, rms_o = min(cands, key=lambda z: z[2])
        out.append(dict(k=b['k'], tau=b['tau'], c=b['c'], w=b['w'], e=b['rms'] / rf, e_opt=rms_o / rf, c_opt=c_o,
                        dc=float(np.abs(b['c'] - c_o).max()) if c_o.size == b['c'].size else np.nan,
                        dc_hs=float(np.abs(b['c'] - (C_B + np.sort(yh))).max()),
                        e_hs=W.l2_error(*W.half_squeeze(b['k'], SIG, ss), SIG, ss)))
    print(f"  B tau/sigma {r_}: {len(out)} blocks, errors " + ', '.join(f"{o['e']:.2e}" for o in out) + f" [{time.time() - t0:.0f}s]",
          flush=True)
    return dict(r=r_, blocks=out, events=[(float(e[0]), int(e[1]), e[2], float(e[3])) for e in R['events']], relaxed=R['relaxed'],
                stiff=R['stiff_redone'], secs=time.time() - t0)


def part_B():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, TAUS_B, chunksize=1)


def c_job(args):
    s, q = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L, _, c0 = TM.setup()
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(TAU_C)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    xi = s * np.random.default_rng(5800 + 100 * LEVELS_C.index(s) + q).standard_normal(X_B.size)
    R = RT.sliding_iss(L, K, f + xi, c0, np.ones(c0.size), polish=True, max_events=3 * K_C, tol=0.0, max_entries=K_C + 1,
                       max_doublings=MAX_DOUBLINGS)
    rf = float(np.sqrt(np.mean(f ** 2)))
    bl = [b for b in blocks(R) if b['k'] > 0]
    ent = [e for e in R['events'] if e[2] == 'enter']
    taus = [e[0] for e in ent]
    out = dict(s=s, q=q, k=[b['k'] for b in bl], rms=[b['rms'] for b in bl],
               err=[float(np.sqrt(np.mean((K.phi(b['c']) @ b['w'] - f) ** 2)) / rf) for b in bl],
               at=[int(np.sum(np.abs(b['c'] - C_B) < 4 * ss)) for b in bl],
               gaps=[float((taus[i + 1] - taus[i]) / taus[i]) for i in range(len(taus) - 1)], exits=sum(e[2] == 'exit' for e in R['events']),
               stiff=R['stiff_redone'], rf=rf, secs=time.time() - t0)
    return out


def part_C():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, [(s, q) for s in LEVELS_C for q in range(N_DRAW_C)], chunksize=1)


def predict_stop(s, rf, N):
    """(i)'s prediction of the discrepancy stop: the first K with e_K²rf² ≤ s²((1 + δ)² − 1 + 2K/N), and the error there."""
    ss = sigma_star(TAU_C)
    for K in range(1, 40):
        e = W.l2_error(*W.half_squeeze(K, SIG, ss), SIG, ss)
        if (e * rf) ** 2 <= s ** 2 * ((1 + DELTA) ** 2 - 1 + 2 * K / N):
            return K, float(np.sqrt(e ** 2 + 2 * K * s ** 2 / (N * rf ** 2)))
    return -1, np.nan


def report(A=None, B=None, C=None):
    L_ = ['# Proposition D‴.58: the transit under width mismatch (`exp_width_mismatch.py`)', '',
          f"Kernels of width σ = {SIG}; a feature Aφ_σ*(· − c*), σ* = σ(1 + (τ/σ)²)^{{1/2}}; t = (σ* − σ)/(σ* + σ), λ = τ²/(σ² + σ*²) "
          "(`src/core/wide_feature.py`).", '']
    if A is not None:
        L_ += ['## A. The rules on ℝ', '',
               "Relative L² errors: the half-squeeze rule (nodes τ(σ/σ*)^{1/2}z_k, amplitudes A(σ*/σ)^{1/2}w_k exp(y_k²/(2σ(σ + σ*)))) by "
               "the trapezoid rule and by the Mehler series (their relative difference), the series' first term's share of the squared "
               "error; moment matching (nodes τz_k, amplitudes A(σ*/σ)w_k); the half-squeeze nodes with least-squares amplitudes; the "
               f"optimal K-kernel fit (variable projection; not computed below {E_FLOOR:g}, the double-precision floor). Rows stop "
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
        from scipy.special import hyp2f1
        L_ += ['', "The tail's share Σ_{n>2K}/(first term) against ₂F₁(−½, −½; 1; t²) − 1 (the share's limit if the ratios' leading "
               "coefficients govern the tail):", '', '| τ/σ | t | K = 5 | 10 | 20 | 30 | ₂F₁(−½, −½; 1; t²) − 1 |', '|---|---|---|---|---|---|---|']
        for r_ in TAUS_A:
            sh = [x['share'] for x in A['share'] if x['r'] == r_]
            t = W.params(SIG, sigma_star(r_))['t']
            L_.append(f"| {r_:g} | {t:.5f} | " + ' | '.join(f"{v:.5f}" for v in sh) + f" | {hyp2f1(-0.5, -0.5, 1, t ** 2) - 1:.5f} |")
    if B is not None:
        L_ += ['', '## B. The transit', '',
               f"The sliding path (polished) on f* = {A_B}φ_σ*(· − {C_B}) sampled at {X_B.size} points of [{X_B[0]:g}, {X_B[-1]:g}], 21 "
               "kernels started at the grid −1, −0.9, …, 1: after each entry the relaxed block's relative rms error, the optimal fit with "
               "as many kernels on the same samples, their ratio and the largest distance between their centres (in σ), the distance to "
               "the half-squeeze nodes (in σ), the continuous half-squeeze error; the block's entry time.", '']
        for b in B:
            L_ += [f"### τ/σ = {b['r']:g} (events: " + ', '.join(f"{j}{'+' if k == 'enter' else '−'}" for t_, j, k, c in b['events'])
                   + f"; relaxed {sum(b['relaxed'])}/{len(b['relaxed'])}; chunks redone by Radau {b['stiff']})", '',
                   '| kernels | τ | block | optimal fit | block / optimal | centres: block − optimal (σ) | block − half-squeeze (σ) | '
                   'half-squeeze on ℝ | ratio per kernel |', '|---|---|---|---|---|---|---|---|---|']
            for i, o in enumerate(b['blocks']):
                rat = f"{b['blocks'][i - 1]['e'] / o['e']:.3f}" if i and b['blocks'][i - 1]['k'] == o['k'] - 1 else '—'
                L_.append(f"| {o['k']} | {o['tau']:.4g} | {o['e']:.4e} | {o['e_opt']:.4e} | {o['e'] / o['e_opt']:.6f} | "
                          f"{o['dc'] / SIG:.1e} | {o['dc_hs'] / SIG:.3f} | {o['e_hs']:.4e} | {rat} |")
            L_.append('')
    if C is not None:
        N = X_B.size
        L_ += ['## C. Under noise', '',
               f"τ/σ = {TAU_C:g}, {N_DRAW_C} draws per level, the path to {K_C} entries: the discrepancy rule's stop (kernels), its "
               "prediction from (i) and the predicted error there, the stopped fit's relative error (median [range]), the kernels within "
               "4σ* of the feature at the stop, the largest relative gap's kernel count.", '',
               '| s | stop: kernels (median [range]) | predicted | error at the stop: median [range] | predicted | kernels at the feature '
               '| largest gap at |', '|---|---|---|---|---|---|---|']
        for s in LEVELS_C:
            ds = [d for d in C if d['s'] == s]
            ks, es, at, kg = [], [], [], []
            for d in ds:
                i = next((i for i, r_ in enumerate(d['rms']) if r_ <= (1 + DELTA) * s), None)
                if i is None:
                    continue
                ks.append(d['k'][i]); es.append(d['err'][i]); at.append(d['at'][i])
                if d['gaps']:
                    kg.append(int(np.argmax(d['gaps'])) + 1)
            Kp, ep = predict_stop(s, ds[0]['rf'], N)
            fmt = lambda v, f_: f"{np.median(v):{f_}} [{np.min(v):{f_}}, {np.max(v):{f_}}]" if len(v) else '—'
            L_.append(f"| {s:g} | {fmt(ks, 'g')} ({len(ks)}/{len(ds)} stopped) | {Kp} | {fmt(es, '.3e')} | {ep:.3e} | {fmt(at, 'g')} | "
                      f"{fmt(kg, 'g')} |")
    L_ += ['', 'Numbers are produced by `experiments/exp_width_mismatch.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'width_mismatch_raw.npz')
    md = os.path.join(OUT, 'width_mismatch.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
