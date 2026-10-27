"""
Theorem D‴.53: the free law with moving centres converges to the sliding inverse scale space at the rate ε^{1/2}
====================================================================================================================================

The loop (the free law with moving centres from α = 0, masses ε; `rich_transit.loop_with_duals`, its dual coordinates integrated
alongside) against its sliding path (`sliding_iss`, polished) on targets whose path is exact, has simple entries, no exits and leaves
every kernel that never enters at a margin ≥ MU0:

A. The rates, per target and ε: the kernels that never enter — their largest coefficient (∝ ε^{3/2}) and their centres' displacement
   against the creep's εC_j (`creep_constants`); the block at the midpoints between entries against S_k (O(ε^{3/2})); the never-entering
   kernels' dual deviation s_j(u_j − u_j^path) against ε^{1/2}·(the summed pulls, `entry_pull`, extrapolated in xi0) + εD_j; the
   entries' crossing times (the entrant's |w| reaching half its block value) against τ_k + ε^{1/2}(s_θ − (its own summed pulls)/v_k).

Usage: python experiments/exp_transit_limit.py [--from-raw] → results/chaos/transit_limit.md (+ _raw.npz).
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
import exp_transit_mass as TM

OUT = os.environ.get('TL_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TL_WORKERS', '2'))
SEED_A, N_POOL, N_A, MU0 = 51, 48, 5, 0.05                             # D‴.51's out-of-sample targets: the N_A with the largest margins (≥ MU0)
EPS_A = (1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11)
XI0 = (400.0, 800.0)                                                   # the pulls' far-field tail ∝ 1/xi0: Richardson from these


def setup():
    return TM.setup()


def path_of(cs, As):
    L, K, c0 = setup()
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    P = RT.path_duals(K, f, R)
    return L, K, c0, f, R, P


def eligible(cs, As):
    """The path exact, entries simple (one kernel per event, distinct times), no exits; the never-entering kernels' smallest margin."""
    L, K, c0, f, R, P = path_of(cs, As)
    ev = R['events']
    if any(k != 'enter' for t, j, k, c in ev) or R['hist'][-1][3] > 1e-8:
        return False, 0.0
    ts = [t for t, j, k, c in ev]
    if len(set(np.round(ts, 12))) < len(ts):
        return False, 0.0
    never = [j for j in range(c0.size) if j not in P['ever']]
    mu = min(1 - abs(P['U'][k][j]) for k in range(len(P['taus'])) for j in never)
    return mu >= MU0, float(mu)


def select():
    """The N_A eligible targets with the largest margins."""
    out = []
    for n, (cs, As) in enumerate(TM.random_targets(SEED_A, N_POOL)):
        ok, mu = eligible(cs, As)
        if ok:
            out.append((n, cs, As, mu))
    return sorted(out, key=lambda t: -t[3])[:N_A]


def a_job(args):
    n, cs, As, mu0 = args
    t0 = time.time()
    L, K, c0, f, R, P = path_of(cs, As)
    d = c0.size
    T = P['taus']
    ents = [(k, float(T[k]), int(j)) for k in range(1, len(T)) for t, j, kind, c in R['events'] if kind == 'enter' and abs(t - T[k]) < 1e-12]
    never = [j for j in range(d) if j not in P['ever']]
    # the pulls on every kernel inactive after each entry, and the entrant's crossing time (half its block value), extrapolated in xi0;
    # the first-order deviation dev_j = (u_j − u_j^path)/ε^{1/2} by recursion over the entries (absolute units: a dual may change sign
    # between entries): entry k arrives shifted by σ_k = −s_k·dev_{j_k}/v_k (s_k its wall's sign, v_k = |a_{k−1,j_k}| its approach speed);
    # each tracer j gains s_j^{(k)}Δ_kj (its pull toward its wall, s_j^{(k)} its dual's sign at the entry) and σ_k(a_{k−1,j} − a_{k,j})
    # (a_{k,j} = ⟨r_kφ(c_j⁰)⟩ the path's rates: the time the shifted entry spends at the old rate)
    pulls, scross, theta, vk, sig, devs = [], [], {}, {}, [], []
    nx = K.xs.size
    dev = {jj: 0.0 for jj in range(d)}
    for q, (k, tk, j) in enumerate(ents):
        later = [jj for (_, _, jj) in ents[q + 1:]]
        tr = never + later
        th = 0.5 * abs(float(R['hist'][k][2][j]))
        theta[j] = th
        r_pre = f - K.phi(R['hist'][k - 1][1]) @ R['hist'][k - 1][2]
        r_post = f - K.phi(R['hist'][k][1]) @ R['hist'][k][2]
        vk[j] = abs(float(K.phi(np.array([c0[j]]))[:, 0] @ r_pre / nx))
        rr = [RT.entry_pull(L, K, f, R, P, k, j, tr, xi0=x0, cross=th) for x0 in XI0]
        pulls.append({jj: 2 * rr[1]['pull'][jj] - rr[0]['pull'][jj] for jj in tr})
        scross.append(2 * rr[1]['s_cross'] - rr[0]['s_cross'])
        sk = -(np.sign(P['U'][k][j]) or 1.0) * dev[j] / vk[j]
        sig.append(sk)
        for jj in tr:
            sj = np.sign(P['U'][k][jj]) or 1.0
            a0 = float(K.phi(np.array([c0[jj]]))[:, 0] @ r_pre / nx)
            a1 = float(K.phi(np.array([c0[jj]]))[:, 0] @ r_post / nx)
            dev[jj] += sj * pulls[-1][jj] + sk * (a0 - a1)
        devs.append(dict(dev))
    # midpoints between entries (the last: past the last entry by half the previous interval) and the end
    tk = [e[1] for e in ents]
    mids = [0.5 * (tk[q] + tk[q + 1]) for q in range(len(tk) - 1)] + [tk[-1] + 0.5 * (tk[-1] - (tk[-2] if len(tk) > 1 else 0.0))]
    t_end = 1.3 * tk[-1]
    grid = np.array(sorted(set(mids + [t_end])))
    Cn = {j: RT.creep_constants(L, K, P, j, c0[j], grid) for j in never}
    pred = []                                                          # the path at the midpoints
    for q, tm in enumerate(mids):
        k = ents[q][0]
        Uj = P['U'][k] + (tm - T[k]) * (K.phi(P['centres'][k]).T @ P['res'][k] / K.xs.size)
        pred.append(dict(tau=tm, w=np.array(R['hist'][k][2]), c=np.array(R['hist'][k][1]), U=Uj,
                         act=[e[2] for e in ents[:q + 1]], sum_pull={j: (np.sign(Uj[j]) or 1.0) * devs[q][j] for j in never},
                         pulls_only={j: (np.sign(Uj[j]) or 1.0) * sum((np.sign(P['U'][ents[l][0]][j]) or 1.0) * pulls[l][j]
                                                                      for l in range(q + 1)) for j in never}))
    rows = []
    for eps in EPS_A:
        t1 = time.time()
        lp = RT.loop_with_duals(L, K, f, c0, eps, grid, levels=theta, P=P)
        se = np.sqrt(eps)
        ix = {float(t): i for i, t in enumerate(grid)}
        blk, dev, dev_c = [], [], []
        for q, p in enumerate(pred):
            i = ix[float(p['tau'])]
            a = p['act']
            blk.append(max(np.abs(lp['w'][i][a] - p['w'][a]).max(), np.abs(lp['c'][i][a] - p['c'][a]).max()))
            for j in never:
                s = np.sign(p['U'][j]) or 1.0
                C, D = Cn[j][0][list(grid).index(p['tau'])], Cn[j][1][list(grid).index(p['tau'])]
                obs = s * lp['du'][i][j] / se
                dev.append((j, q, float(obs), float(p['sum_pull'][j]), float(s * eps * D / se), float(p['pulls_only'][j])))
        ie = ix[float(t_end)]
        wn = float(np.abs(lp['w'][ie][never]).max())
        disp = np.array([lp['c'][ie][j] - c0[j] for j in never])
        Cend = np.array([Cn[j][0][-1] for j in never])
        cr = []
        for q, (k, tk_, j) in enumerate(ents):
            cr.append((j, (lp['cross'][j] - tk_) / se if lp['cross'][j] is not None else np.nan, scross[q] + sig[q]))
        n_form = int(np.sum(np.abs(lp['w'][ie]) > X.TH_ACT * np.abs(As).max()))
        rows.append(dict(eps=eps, block=float(max(blk)), w_never=wn, disp=disp.tolist(), creep=(eps * Cend).tolist(), dev=dev, cross=cr,
                         n_form=n_form, secs=time.time() - t1))
        print(f"  target {n} eps {eps:g}: formed {n_form}, block {max(blk):.2e}, w_never {wn:.2e} [{time.time() - t1:.0f}s]", flush=True)
    return dict(n=n, cs=list(cs), As=list(As), mu0=mu0, entries=[(e[1], e[2]) for e in ents], scross=scross, shifts=sig, pulls=[
        {int(k): float(v) for k, v in p.items()} for p in pulls], rows=rows, secs=time.time() - t0)


def part_A():
    from multiprocessing import Pool
    sel = select()
    print(f"selected {[s[0] for s in sel]} (margins {[round(s[3], 3) for s in sel]})", flush=True)
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, sel)


def slope(x, y):
    x, y = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    return float(np.polyfit(x, y, 1)[0]) if x.size >= 2 else np.nan


def report(A=None):
    fmt = E.fmt
    L_ = ['# Theorem D‴.53: the sliding limit (`exp_transit_limit.py`)', '',
          "The free law with moving centres (m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025) from α = 0 with the masses ε, Gaussian kernels of "
          "width 0.12 on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, against its sliding path (polished); the "
          f"targets: of D‴.51's {N_POOL} out-of-sample three-bump targets (seed {SEED_A}) whose path is exact with simple entries and no exits, "
          f"the {N_A} with the largest margins (≥ {MU0}) for the kernels that never enter.", '']
    if A is not None:
        L_ += ['## A. The rates', '',
               "Per target and ε: the kernels formed; the block at the midpoints between entries against S_k (max over the coefficients and "
               "centres); the never-entering kernels' largest coefficient (× ε^{−3/2}); their centres' largest displacement against the "
               "creep's εC_j (the largest ratio's distance from 1); the duals: obs = s_j(u_j − u_j^path)/ε^{1/2} against the first-order "
               "deviation — the entries' pulls Δ_kj plus each entry's shift σ_k = −dev_{j_k}/v_k times the drop of j's rate at it, by "
               "recursion over the entries — plus the creep's ε^{1/2}D_j (the largest |obs − pred|, the median ratio obs/pred where "
               "|pred| > 0.1, and the largest |obs − ΣΔ| with the pulls alone); the entries' crossing times (τ_θ − τ_k)/ε^{1/2} against "
               "s_θ + σ_k (the largest difference).", '',
               '| target (margin) | ε | formed | block − S_k | max \\|w_never\\| (× ε^{−3/2}) | displacement / εC − 1 | duals: max \\|obs − pred\\| '
               '(median obs/pred) | the pulls alone | crossings: max \\|obs − pred\\| |',
               '|---|---|---|---|---|---|---|---|---|']
        for t in A:
            for r in t['rows']:
                e = r['eps']
                disp, cr_ = np.array(r['disp']), np.array(r['creep'])
                big = np.abs(cr_) > 1e-3 * np.abs(cr_).max()
                rel = float(np.abs(disp[big] / cr_[big] - 1).max()) if big.any() else np.nan
                obs = np.array([d[2] for d in r['dev']])
                pp = np.array([d[3] for d in r['dev']])
                cc = np.array([d[4] for d in r['dev']])
                po = np.array([d[5] for d in r['dev']])
                pr = pp + cc
                m = np.abs(pr) > 0.1
                med = float(np.median(obs[m] / pr[m])) if m.any() else np.nan
                xs = [abs(a - b) for _, a, b in r['cross'] if np.isfinite(a)]
                L_.append(f"| {t['n']} ({t['mu0']:.3f}) | {e:g} | {r['n_form']} | {fmt(r['block'])} | {fmt(r['w_never'])} "
                          f"({r['w_never'] / e ** 1.5:.3g}) | {rel:.3g} | {np.abs(obs - pr).max():.4f} ({med:.4f}) | "
                          f"{np.abs(obs - po).max():.4f} | {max(xs) if xs else np.nan:.4f} |")
        L_ += ['', 'Slopes in log ε (least squares over ε): ' + '; '.join(
            f"target {t['n']}: block {slope([r['eps'] for r in t['rows']], [r['block'] for r in t['rows']]):.3f}, w_never "
            f"{slope([r['eps'] for r in t['rows']], [r['w_never'] for r in t['rows']]):.3f}, displacement "
            f"{slope([r['eps'] for r in t['rows']], [np.abs(r['disp']).max() for r in t['rows']]):.3f}, duals' residual "
            f"{slope([r['eps'] for r in t['rows']], [max(abs(d[2] - d[3] - d[4]) for d in r['dev']) for r in t['rows']]):.3f}"
            for t in A) + ' (predicted 1.5, 1.5, 1 and, for the residual in units of ε^{1/2}, ≥ 0.25).', '']
        L_ += ['The entries (τ_k, kernel), s_θ and the pulls on the never-entering kernels (largest |Δ|): ' + '; '.join(
            f"target {t['n']}: " + ', '.join(f"({tk:.4g}, {j}) s_θ {s:.4g} max|Δ| {max(abs(v) for v in p.values()):.3f}"
                                             for (tk, j), s, p in zip(t['entries'], t['scross'], t['pulls'])) for t in A) + '.', '']
    L_.append('Numbers are produced by `experiments/exp_transit_limit.py`; nothing is edited by hand.')
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'transit_limit_raw.npz')
    md = os.path.join(OUT, 'transit_limit.md')
    old = E.load(raw)
    keys = ('A',)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A),):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
