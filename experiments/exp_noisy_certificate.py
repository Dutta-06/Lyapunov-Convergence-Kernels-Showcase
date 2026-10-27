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
D. Stopping and the oracle: the discrepancy rule (stop at the first entry after which the relaxed residual's rms is ≤ (1 + DELTA)s) and
   the largest relative gap (τ_{k+1} − τ_k)/τ_k over the K + 1 entries; the fit stopped at the K-th entry — the K-kernel least-squares
   fit of y with free centres — against f*: its error against the projection P_Tξ of the noise on the tangent space and against
   s(2K/N)^{1/2}.

Usage: python experiments/exp_noisy_certificate.py [--from-raw] [--only A,B] → results/chaos/noisy_certificate.md (+ _raw.npz).
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


def noisy(L, K, c0, f, xi, nK):
    return RT.sliding_iss(L, K, f + xi, c0, np.ones(c0.size), polish=True, max_events=nK + 3, tol=0.0, max_entries=nK + 1)


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
            exits = sum(e[2] == 'exit' for e in R['events'])
            first = sorted(e[1] for e in ent[:nK])
            exact = first == E_ and exits == 0 and len(ent) >= nK
            crit = float(max(np.abs(P['u_all'][k][off] + P['U_all'][k][off] @ xi).max() for k in range(nK)))
            d = dict(exact=exact, crit=crit, exits=exits)
            if exact:
                v = Phi0.T @ (P['Rlast'] @ xi) / nx
                u = P['u'] + P['U'] @ xi
                W = {j: (1 - np.sign(v[j]) * u[j]) / abs(v[j]) for j in off if v[j] != 0}
                jn = min(W, key=W.get)
                taus = [e[0] for e in ent]
                d.update(tauK=taus[nK - 1], W_pred=float(W[jn]), j_pred=int(jn),
                         W_obs=float(taus[nK] - taus[nK - 1]) if len(ent) > nK else np.nan, j_obs=int(ent[nK][1]) if len(ent) > nK else -1,
                         gaps=[float((taus[i + 1] - taus[i]) / taus[i]) for i in range(len(taus) - 1)])
                rms = [h[3] for h in R['hist'][1:nK + 1]]                    # the relaxed residuals after entries 1..K
                d['k_disc'] = next((i + 1 for i, r_ in enumerate(rms) if r_ <= (1 + DELTA) * s), -1)
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


def report(A=None, B=None):
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
               "On the exact draws: the first noise kernel's entry after the K-th, τ_{K+1} − τ_K, against the first-order W (the ratio's median "
               "and range), the kernel as predicted, and W·s (median, quartiles).", '',
               '| s | exact draws | W_obs/W_pred: median [min, max] | next kernel as predicted | W·s: median [quartiles] |', '|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            ds = [d for b in B for d in b['levels'][i]['draws'] if d['exact'] and np.isfinite(d['W_obs'])]
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — |")
                continue
            r = np.array([d['W_obs'] / d['W_pred'] for d in ds])
            ws = np.array([d['W_obs'] * s for d in ds])
            L_.append(f"| {s:g} | {len(ds)} | {np.median(r):.4f} [{r.min():.3f}, {r.max():.3f}] | "
                      f"{sum(d['j_obs'] == d['j_pred'] for d in ds)}/{len(ds)} | {np.median(ws):.3g} [{np.percentile(ws, 25):.3g}, "
                      f"{np.percentile(ws, 75):.3g}] |")
        L_ += ['', '## D. Stopping and the oracle', '',
               f"On the exact draws: the discrepancy rule (the first entry after which the relaxed residual's rms ≤ {1 + DELTA:g}s) stopping at "
               "the K-th; the largest relative gap (τ_{k+1} − τ_k)/τ_k over the K + 1 entries falling after the K-th; the fit stopped at the "
               "K-th entry: its rms error against the noise's projection on the tangent space, rms(P_Tξ), and against s(2K/N)^{1/2}.", '',
               '| s | exact draws | discrepancy rule at K | largest gap at K | error / rms(P_Tξ): median [range] | error / (s(2K/N)^{1/2}): median |',
               '|---|---|---|---|---|---|']
        for i, s in enumerate(LEVELS):
            ds = [d for b in B for d in b['levels'][i]['draws'] if d['exact']]
            if not ds:
                L_.append(f"| {s:g} | 0 | — | — | — | — |")
                continue
            nK = 3
            e1 = np.array([d['err'] / d['err_T'] for d in ds])
            e2 = np.array([d['err'] / (s * np.sqrt(2 * nK / 61)) for d in ds])
            gap = [int(np.argmax(d['gaps'])) == nK - 1 for d in ds if len(d['gaps']) >= nK]
            L_.append(f"| {s:g} | {len(ds)} | {sum(d['k_disc'] == nK for d in ds)}/{len(ds)} | {sum(gap)}/{len(gap)} | "
                      f"{np.median(e1):.4f} [{e1.min():.3f}, {e1.max():.3f}] | {np.median(e2):.3f} |")
    L_ += ['', 'Numbers are produced by `experiments/exp_noisy_certificate.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'noisy_certificate_raw.npz')
    md = os.path.join(OUT, 'noisy_certificate.md')
    old = E.load(raw)
    keys = ('A', 'B')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
