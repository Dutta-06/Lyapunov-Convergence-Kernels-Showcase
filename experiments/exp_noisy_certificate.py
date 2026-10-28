"""
Proposition D‴.56: the transit under measurement noise — the path's linear response decides the support, the window and the stop
====================================================================================================================================

The sliding path (`sliding_iss`, polished) on noisy data y = f* + ξ against its first-order response along the noiseless path
(`rich_transit.path_response`): the duals at the entries u⁰ + Uξ, the entry times τ⁰ + Tξ and, after the K-th entry, the residual R_Kξ
— the projection of ξ off the K-kernel fit's tangent space —, which sets the first noise kernel's entry. Targets: D‴.55's recipe, three
features of the kernels' width at least DMIN apart, amplitudes ±U(0.4, 1), whose certificate has m > 0.8q (outside D‴.55's band, so
that without noise the path enters exactly the nearest kernels E).

A. The response is the first order: per target, one noise direction ξ at the levels h: the noisy path's entry times and its duals at the
   K-th entry against τ⁰ + hTξ and u⁰ + hUξ — their errors ∝ h².
B. The support: per target and noise level s, N_DRAW draws of i.i.d. Gaussian noise of standard deviation s on the samples: whether the
   path enters exactly E first, against the first-order criterion (per draw: max over the entries and the kernels j ∉ E of
   |u⁰_j + U_jξ| < 1) and against the linear model's probability of it (Monte Carlo over its Gaussian law).
C. The window: on the draws with the exact support, the first noise kernel's entry τ_{K+1} − τ_K against W = min_j (1 − sign(v_j)(u⁰_j +
   U_jξ))/|v_j|, v_j = ⟨R_Kξ, φ_j⟩, and the kernel that enters; W·s, whose law the first order makes independent of s.
D. Stopping and the oracle: the discrepancy rule (stop at the first entry after which the relaxed residual's rms is ≤ (1 + DELTA)s)
   against its first-order prediction (at the K-th iff rms(R_Kξ) ≤ (1 + DELTA)s; late with the probability P(χ²_{N−2K} > (1 + DELTA)²N))
   and the largest relative gap (τ_{k+1} − τ_k)/τ_k over the K + 1 entries; the fit stopped at the K-th entry — the K-kernel
   least-squares fit of y with free centres — against f*: its error against the projection P_Tξ of the noise on the tangent space and
   against s(2K/N)^{1/2}.
E. Past the first noise kernel (N_E targets, LEVELS_E, N_DRAW_E draws, the path to K + N_MORE entries, each under T_PATH s): the
   largest gap over all those entries, the discrepancy rule over the relaxed fits, and the fits the rules stop at against the oracle's.

Usage: python experiments/exp_noisy_certificate.py [--from-raw] [--only A,B,E] → results/chaos/noisy_certificate.md (+ _raw.npz).
NC_WORKERS (default 2).
"""
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_mirror_transit as X
import exp_transit_certificate as TC
import exp_transit_mass as TM

OUT = os.environ.get('NC_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('NC_WORKERS', '2'))
SEED, N_T, DMIN = 561, 40, 0.45
N_A, H_A = 10, (1e-3, 1e-4, 1e-5)
LEVELS, N_DRAW, N_MC, DELTA = (0.003, 0.01, 0.03, 0.1, 0.3), 10, 20000, 0.1
N_E, LEVELS_E, N_DRAW_E, N_MORE, T_PATH = 20, (0.01, 0.03, 0.1, 0.3), 5, 4, 180    # E: the path past the first noise kernel


def targets():
    rng = np.random.default_rng(SEED)
    L, K, c0 = TM.setup()
    out = []
    while len(out) < N_T:
        (cs, sg, a, u), = TC.targets(rng, 1, DMIN)
        As = sg * a
        c = TC.certificate(K, c0, cs, As)
        if c['unique'] and c['m'] > 0.8 * c['q']:
            out.append((cs, As))
    return out


def noiseless(cs, As):
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, cs, As)
    cert = TC.certificate(K, c0, cs, As)
    R0 = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    return L, K, c0, f, cert, R0, RT.path_response(K, f, R0, c0)


def noisy(L, K, c0, f, xi, nK, more=1):
    """The path on f* + ξ up to its (K + more)-th entry."""
    return RT.sliding_iss(L, K, f + xi, c0, np.ones(c0.size), polish=True, max_events=nK + more + 2, tol=0.0, max_entries=nK + more)


def a_job(args):
    q, cs, As = args
    L, K, c0, f, cert, R0, P = noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    xi = np.random.default_rng(1000 + q).standard_normal(nx)
    off = [j for j in range(c0.size) if j not in P['entrants']]
    rows = []
    for h in H_A:
        R = noisy(L, K, c0, f, h * xi, nK)
        ent = [e for e in R['events'] if e[2] == 'enter']
        t_obs = np.array([e[0] for e in ent[:nK]])
        u_obs = RT.path_duals(K, f + h * xi, R)['U'][nK]
        dt = t_obs - (P['taus'] + h * P['T'] @ xi)
        du = (u_obs - (P['u'] + h * P['U'] @ xi))[off]
        rows.append(dict(h=h, same=[e[1] for e in ent[:nK]] == P['entrants'], dt=float(np.abs(dt).max()),
                         dt1=float(h * np.abs(P['T'] @ xi).max()), du=float(np.abs(du).max()), du1=float(h * np.abs((P['U'] @ xi)[off]).max())))
    return dict(q=q, rows=rows)


def b_job(args):
    q, cs, As = args
    t0 = time.time()
    warnings.filterwarnings('ignore')
    L, K, c0, f, cert, R0, P = noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    E_ = sorted(cert['e'])
    off = [j for j in range(c0.size) if j not in E_]
    Phi0 = K.phi(c0)
    # the linear model's law of the support's criterion: the duals at every entry, the kernels off E
    Ust = np.vstack([P['U_all'][k][off] for k in range(nK)])
    ust = np.concatenate([P['u_all'][k][off] for k in range(nK)])
    rng = np.random.default_rng(2000 + q)
    out = dict(q=q, m=cert['m'], q_ct=cert['q'], path_exact=sorted(P['entrants']) == E_, levels=[])
    for s in LEVELS:
        Z = s * rng.standard_normal((nx, N_MC))
        p_mc = float(np.mean(np.abs(ust[:, None] + Ust @ Z).max(0) < 1))
        draws = []
        for k_ in range(N_DRAW):
            xi = s * rng.standard_normal(nx)
            R = noisy(L, K, c0, f, xi, nK)
            ent = [e for e in R['events'] if e[2] == 'enter']
            i1 = [i for i, e in enumerate(R['events']) if e[2] == 'enter'][nK] if len(ent) > nK else len(R['events'])
            exits = sum(e[2] == 'exit' for e in R['events'][:i1])          # up to the (K + 1)-th entry
            first = sorted(e[1] for e in ent[:nK])
            exact = first == E_ and exits == 0 and len(ent) >= nK
            crit = float(max(np.abs(P['u_all'][k][off] + P['U_all'][k][off] @ xi).max() for k in range(nK)))
            # the first-order residual after the K-th entry, R_Kξ, on every draw (the discrepancy rule's prediction)
            d = dict(exact=exact, crit=crit, exits=exits, rms_RK=float(np.sqrt(np.mean((P['Rlast'] @ xi) ** 2))))
            if exact:
                v = Phi0.T @ (P['Rlast'] @ xi) / nx
                u = P['u'] + P['U'] @ xi
                W = {j: (1 - np.sign(v[j]) * u[j]) / abs(v[j]) for j in off if v[j] != 0}
                jn = min(W, key=W.get)
                Ws = sorted(W.values())
                taus = [e[0] for e in ent]
                d.update(tauK=taus[nK - 1], W_pred=float(W[jn]), j_pred=int(jn), W_gap=float((Ws[1] - Ws[0]) / abs(Ws[0])),
                         W_obs=float(taus[nK] - taus[nK - 1]) if len(ent) > nK else np.nan, j_obs=int(ent[nK][1]) if len(ent) > nK else -1,
                         gaps=[float((taus[i + 1] - taus[i]) / taus[i]) for i in range(len(taus) - 1)])
                # the relaxed residuals and the fits' errors after entries 1, …, K (the (K + 1)-th entry is not relaxed)
                rms = [h[3] for h in R['hist'][1:]]
                d['k_disc'] = next((i + 1 for i, r_ in enumerate(rms) if r_ <= (1 + DELTA) * s), -1)
                d['rms_K'] = float(rms[nK - 1])
                d['errs'] = [float(np.sqrt(np.mean((K.phi(h[1]) @ h[2] - f) ** 2))) for h in R['hist'][1:]]
                d['exits_all'] = sum(e[2] == 'exit' for e in R['events'])
                hK = R['hist'][nK]
                fhat = K.phi(hK[1]) @ hK[2]
                d['err'] = float(np.sqrt(np.mean((fhat - f) ** 2)))
                d['err_T'] = float(np.sqrt(np.mean(((np.eye(nx) - P['Rlast']) @ xi) ** 2)))
            draws.append(d)
        out['levels'].append(dict(s=s, p_mc=p_mc, draws=draws))
    out['secs'] = time.time() - t0
    print(f"  target {q}: m {cert['m']:.3f}, q {cert['q']:.3f}, exact " +
          ', '.join(f"{sum(d['exact'] for d in lv['draws'])}" for lv in out['levels']) + f" [{out['secs']:.0f}s]", flush=True)
    return out


def e_job(args):
    """E: the path past the first noise kernel, to its (K + N_MORE)-th entry, each path under a wall-clock limit T_PATH: on the draws whose
    first K entries are E (no exits before the (K + 1)-th), the entry times, the relaxed residuals and the fits' errors after each entry."""
    import signal
    q, cs, As = args
    warnings.filterwarnings('ignore')
    L, K, c0, f, cert, R0, P = noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    E_ = sorted(cert['e'])
    rng = np.random.default_rng(3000 + q)

    def alarm(signum, frame):
        raise TimeoutError
    signal.signal(signal.SIGALRM, alarm)
    out = dict(q=q, levels=[])
    for s in LEVELS_E:
        draws = []
        for k_ in range(N_DRAW_E):
            xi = s * rng.standard_normal(nx)
            signal.alarm(T_PATH)
            try:
                R = noisy(L, K, c0, f, xi, nK, N_MORE)
            except TimeoutError:
                draws.append(dict(timeout=True))
                continue
            finally:
                signal.alarm(0)
            ent = [e for e in R['events'] if e[2] == 'enter']
            i1 = [i for i, e in enumerate(R['events']) if e[2] == 'enter'][nK] if len(ent) > nK else len(R['events'])
            exact = sorted(e[1] for e in ent[:nK]) == E_ and not any(e[2] == 'exit' for e in R['events'][:i1])
            d = dict(timeout=False, exact=exact, n_ent=len(ent), exits=sum(e[2] == 'exit' for e in R['events']), stiff=R['stiff_redone'])
            if exact:
                taus = [e[0] for e in ent]
                rms = [h[3] for h in R['hist'][1:]]
                d.update(gaps=[float((taus[i + 1] - taus[i]) / taus[i]) for i in range(len(taus) - 1)],
                         k_disc=next((i + 1 for i, r_ in enumerate(rms) if r_ <= (1 + DELTA) * s), -1),
                         errs=[float(np.sqrt(np.mean((K.phi(h[1]) @ h[2] - f) ** 2))) for h in R['hist'][1:]])
            draws.append(d)
        out['levels'].append(dict(s=s, draws=draws))
    print(f"  E target {q}: " + ', '.join(f"{sum(d.get('exact', False) for d in lv['draws'])}/{sum(not d['timeout'] for d in lv['draws'])}"
                                         for lv in out['levels']), flush=True)
    return out


def part_E():
    from multiprocessing import Pool
    tg = targets()[:N_E]
    with Pool(N_WORKERS) as pool:
        return pool.map(e_job, [(q, cs, As) for q, (cs, As) in enumerate(tg)], chunksize=1)


def part_A():
    from multiprocessing import Pool
    tg = targets()[:N_A]
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, [(q, cs, As) for q, (cs, As) in enumerate(tg)], chunksize=1)


def part_B():
    from multiprocessing import Pool
    tg = targets()
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(q, cs, As) for q, (cs, As) in enumerate(tg)], chunksize=1)


def report(A=None, B=None, E_=None):
    L_ = ['# Proposition D‴.56: the transit under measurement noise (`exp_noisy_certificate.py`)', '',
          "The sliding path on y = f* + ξ (21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]) against its "
          f"first-order response along the noiseless path (`path_response`); targets: three features at least {DMIN} apart (D‴.55's "
          "recipe), certificate margin m > 0.8q.", '']
    if A is not None:
        L_ += ['## A. The response is the first order', '',
               "Per target (one noise direction ξ, unit variance per sample) and level h: the noisy path's entry times and its duals at the "
               "K-th entry against the first-order prediction — the largest error, and the first-order term's own size.", '',
               '| target | h | same entries | entry times: error (first-order term) | duals: error (first-order term) |', '|---|---|---|---|---|']
        for a in A:
            for r in a['rows']:
                L_.append(f"| {a['q']} | {r['h']:g} | {r['same']} | {r['dt']:.2e} ({r['dt1']:.2e}) | {r['du']:.2e} ({r['du1']:.2e}) |")
        rat = [a['rows'][i]['du'] / a['rows'][i + 1]['du'] for a in A for i in range(len(a['rows']) - 1) if a['rows'][i + 1]['du'] > 0]
        L_ += ['', f"The duals' error falls by {min(rat):.1f}–{max(rat):.1f} per decade of h (second order: 100).", '']
    if B is not None:
        L_ += ['## B. The support', '',
               f"Per noise level s: {N_DRAW} draws on each of the {len(B)} targets ({sum(b['path_exact'] for b in B)} with the noiseless path "
               "exactly E). Observed: the path enters exactly E first; predicted per draw: the first-order criterion max|u⁰_j + U_jξ| < 1 "
               f"over the entries and the kernels off E; the expected count from the linear model's probability ({N_MC} Monte Carlo draws "
               "per target).", '',
               '| s | exact (observed) | criterion < 1 | agree per draw | expected (linear model) | observed exact, criterion ≥ 1 | '
               'not exact, criterion < 1 |', '|---|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            ds = [d for b in B for d in b['levels'][i]['draws']]
            ex = [d['exact'] for d in ds]
            cr = [d['crit'] < 1 for d in ds]
            L_.append(f"| {s:g} | {sum(ex)}/{len(ds)} | {sum(cr)} | {sum(a == b for a, b in zip(ex, cr))} | "
                      f"{sum(N_DRAW * b['levels'][i]['p_mc'] for b in B):.1f} | {sum(a and not b for a, b in zip(ex, cr))} | "
                      f"{sum(b and not a for a, b in zip(ex, cr))} |")
        L_ += ['', '## C. The window', '',
               "On the exact draws: the first noise kernel's entry after the K-th, τ_{K+1} − τ_K, against the first-order W — the ratio's median, "
               "range and 5–95% quantiles, and their half-width over s (constant if the relative error is O(s)) —; the kernel as predicted, and "
               "on the draws where it is not, the relative gap between the two smallest predicted windows (W₍₂₎ − W₍₁₎)/W₍₁₎ (median; all "
               "draws' median in brackets); W·s (median, quartiles).", '',
               '| s | exact draws | W_obs/W_pred: median [min, max] | 5–95% | half-width / s | next kernel as predicted | '
               'missed: (W₍₂₎ − W₍₁₎)/W₍₁₎ [all] | W·s: median [quartiles] |', '|---|---|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            ds = [d for b in B for d in b['levels'][i]['draws'] if d['exact'] and np.isfinite(d['W_obs'])]
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — | — | — | — |")
                continue
            r = np.array([d['W_obs'] / d['W_pred'] for d in ds])
            q5, q95 = np.percentile(r, [5, 95])
            ws = np.array([d['W_obs'] * s for d in ds])
            miss = [d['W_gap'] for d in ds if d['j_obs'] != d['j_pred']]
            mg = f"{np.median(miss):.3g}" if miss else '—'
            L_.append(f"| {s:g} | {len(ds)} | {np.median(r):.4f} [{r.min():.3f}, {r.max():.3f}] | [{q5:.4f}, {q95:.4f}] | "
                      f"{(q95 - q5) / 2 / s:.2f} | {sum(d['j_obs'] == d['j_pred'] for d in ds)}/{len(ds)} | {mg} "
                      f"[{np.median([d['W_gap'] for d in ds]):.3g}] | {np.median(ws):.3g} [{np.percentile(ws, 25):.3g}, "
                      f"{np.percentile(ws, 75):.3g}] |")
        nK, nx = 3, 61
        from scipy.stats import chi2
        p_late = float(chi2.sf((1 + DELTA) ** 2 * nx, nx - 2 * nK))
        L_ += ['', '## D. Stopping and the oracle', '',
               f"The discrepancy rule: the first entry after which the relaxed residual's rms is ≤ {1 + DELTA:g}s, on the exact draws: at the "
               f"K-th, not by the K-th (late), earlier; its first-order prediction — at the K-th iff rms(R_Kξ) ≤ {1 + DELTA:g}s, late otherwise, "
               f"R_Kξ the first-order residual (|R_Kξ|²/s² ~ χ²_(N−2K)) — and its agreement per draw; the first-order probability of a late "
               f"stop, P(χ²_(N−2K) > {1 + DELTA:g}²N) = {p_late:.4f}, times the draws (on all draws, exact or not, the count with rms(R_Kξ) > "
               "(1 + δ)s in brackets); the relaxed residual after the K-th entry against rms(R_Kξ).", '',
               '| s | exact draws | at K | late | early | first order: at K | agree per draw | late: expected [all draws] | '
               'rms_K / rms(R_Kξ): median [5–95%] |', '|---|---|---|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            al = [d for b in B for d in b['levels'][i]['draws']]
            ds = [d for d in al if d['exact']]
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — | — | — | — | — |")
                continue
            kd = np.array([d['k_disc'] for d in ds])
            fo = np.array([d['rms_RK'] <= (1 + DELTA) * s for d in ds])
            rr = np.array([d['rms_K'] / d['rms_RK'] for d in ds])
            L_.append(f"| {s:g} | {len(ds)} | {int(np.sum(kd == nK))} | {int(np.sum((kd > nK) | (kd < 0)))} | "
                      f"{int(np.sum((kd > 0) & (kd < nK)))} | {int(fo.sum())} | {int(np.sum(fo == (kd == nK)))} | "
                      f"{p_late * len(ds):.1f} [{sum(d['rms_RK'] > (1 + DELTA) * s for d in al)}/{len(al)}] | "
                      f"{np.median(rr):.4f} [{np.percentile(rr, 5):.4f}, {np.percentile(rr, 95):.4f}] |")
        L_ += ['', "The largest relative gap (τ_{k+1} − τ_k)/τ_k over the first K + 1 entries falling after the K-th; the fit stopped at the "
               "K-th entry: its rms error against the noise's projection on the tangent space, rms(P_Tξ), and against s(2K/N)^{1/2}.", '',
               '| s | exact draws | largest gap at K | error / rms(P_Tξ): median [range] | 5–95% | error / (s(2K/N)^{1/2}): median |',
               '|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            ds = [d for b in B for d in b['levels'][i]['draws'] if d['exact']]
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — | — |")
                continue
            e1 = np.array([d['err'] / d['err_T'] for d in ds])
            e2 = np.array([d['err'] / (s * np.sqrt(2 * nK / nx)) for d in ds])
            g1 = [int(np.argmax(d['gaps'][:nK])) == nK - 1 for d in ds if len(d['gaps']) >= nK]
            L_.append(f"| {s:g} | {len(ds)} | {sum(g1)}/{len(g1)} | {np.median(e1):.4f} [{e1.min():.3f}, {e1.max():.3f}] | "
                      f"[{np.percentile(e1, 5):.4f}, {np.percentile(e1, 95):.4f}] | {np.median(e2):.3f} |")
    if E_ is not None:
        nK = 3
        L_ += ['', '## E. Past the first noise kernel', '',
               f"The first {len(E_)} targets, {N_DRAW_E} draws per level, the path to its (K + {N_MORE})-th entry (each path under {T_PATH} s): "
               "on the draws whose first K entries are E, the largest relative gap over the K + "
               f"{N_MORE} entries falling after the K-th, the discrepancy rule over the relaxed fits after entries 1 … K + {N_MORE - 1} (at "
               "the K-th, later, never), and the fits the two rules stop at against the oracle's (the K-th) error. "
               + (f"Paths with a chunk of the fast flow redone by Radau (`relax_chunk`, LSODA's stall): "
                  f"{sum(d.get('stiff', 0) > 0 for e in E_ for lv in e['levels'] for d in lv['draws'])} of "
                  f"{sum(not d['timeout'] for e in E_ for lv in e['levels'] for d in lv['draws'])} finished."
                  if any('stiff' in d for e in E_ for lv in e['levels'] for d in lv['draws']) else ''), '',
               '| s | paths (over the limit) | exact | largest gap at K | discrepancy: at K, later, never | discrepancy-stopped error / oracle: '
               'median [max] | gap-stopped error / oracle: median [max] |', '|---|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS_E):
            al = [d for e in E_ for d in e['levels'][i]['draws']]
            ds = [d for d in al if not d['timeout'] and d['exact'] and len(d['gaps']) >= nK + N_MORE - 1]
            if not ds:
                L_.append(f"| {s:g} | {len(al)} ({sum(d['timeout'] for d in al)}) | 0 | — | — | — | — |")
                continue
            gN = [int(np.argmax(d['gaps'])) == nK - 1 for d in ds]
            kd = np.array([d['k_disc'] for d in ds])
            ed = np.array([d['errs'][d['k_disc'] - 1] / d['errs'][nK - 1] for d in ds if d['k_disc'] > 0])
            eg = np.array([d['errs'][int(np.argmax(d['gaps']))] / d['errs'][nK - 1] for d in ds])
            L_.append(f"| {s:g} | {len(al)} ({sum(d['timeout'] for d in al)}) | {len(ds)} | {sum(gN)}/{len(gN)} | "
                      f"{int(np.sum(kd == nK))}, {int(np.sum(kd > nK))}, {int(np.sum(kd < 0))} | "
                      + (f"{np.median(ed):.3f} [{ed.max():.3g}]" if ed.size else '—') + f" | {np.median(eg):.3f} [{eg.max():.3g}] |")
    L_ += ['', 'Numbers are produced by `experiments/exp_noisy_certificate.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'noisy_certificate_raw.npz')
    md = os.path.join(OUT, 'noisy_certificate.md')
    old = E.load(raw)
    keys = ('A', 'B', 'E')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('E', part_E)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
