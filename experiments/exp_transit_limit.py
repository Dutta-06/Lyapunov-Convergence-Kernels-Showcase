"""
Theorem D‴.53: the free law with moving centres converges to the sliding inverse scale space at the rate ε^{1/2}
====================================================================================================================================

The loop (the free law with moving centres from α = 0, masses ε; `rich_transit.loop_with_duals`, its dual coordinates integrated
alongside) against its sliding path (`sliding_iss`, polished) on D‴.51's out-of-sample three-bump targets (every one of them has an exact
path with simple entries and no exits):

B. The constants, for all the targets: per entry, the entrant's approach speed v_k = |⟨r_{k−1}φ(c⁰)⟩|, the sensitivity of that speed to
   its centre g′_k = |⟨r_{k−1}∂_cφ(c⁰)⟩| (the far field's creep tail ∝ g′²/v³ in the inner problem), its slide |c^k − c⁰|, the gap to the
   next entry, the inner problem's crossing time s_θ (the entrant's |w| at half its block value) with the entry's shift σ_k, the largest
   pull, the block's slowest relaxation rate λ_k in the inner time (`block_rate`), and (H3) — the inner problem's end against the path's
   S_k. `entry_pull` in its deviation form (the duals' departures from their outer solutions integrated) with Radau at the relative
   tolerance PULL_RTOL (converged to ~10⁻⁸ in s_θ; LSODA's steps collapse on some entries in this form); the pulls and s_θ extrapolated in xi0 — their far-field tail is a/xi0 + b/xi0² — by
   second-order Richardson from the last three of XI0 (xi0, 2xi0, 4xi0 ↦ (f₁ − 6f₂ + 8f₄)/3), its spread against the same from the
   first three reported.
A. The rates, on the N_A targets with the largest margins and the stress targets (the slowest entrant, the largest slide, the longest
   inner layer against its gap): per target and ε, the kernels that never enter — their largest coefficient (∝ ε^{3/2}) and their
   centres' displacement against the creep's εC_j (`creep_constants`); the block against S_k at the midpoints between entries and at the
   ends of (b)'s window [τ_k + C_Bε^{1/2}log(1/ε), τ_{k+1} − ε^{1/4}]; the never-entering kernels' dual deviation s_j(u_j − u_j^path)
   against ε^{1/2}·(the summed pulls and shifts) + εD_j; the entries' crossing times against τ_k + ε^{1/2}(s_θ + σ_k). The loop at the
   relative tolerance loop_rtol(ε) (10⁻¹³ for ε ≤ 10⁻¹⁰: referee round 41).

Usage: python experiments/exp_transit_limit.py [--from-raw] [--only B,A] → results/chaos/transit_limit.md (+ _raw.npz). TL_WORKERS.
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
SEED_A, N_POOL, N_A, MU0 = 51, 48, 5, 0.05                             # D‴.51's out-of-sample targets: the N_A with the largest margins
EPS_A = (1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11)
XI0 = (800.0, 1600.0, 3200.0, 6400.0)                                  # the pulls' tail a/xi0 + b/xi0²: second-order Richardson
PULL_RTOL, PULL_METHOD = 1e-10, 'Radau'
C_B = 30.0                                                             # (b)'s window starts τ_k + C_B·ε^{1/2}log(1/ε)


def loop_rtol(eps):
    return 1e-13 if eps <= 1e-10 else 1e-12 if eps <= 1e-9 else 1e-11


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


def entries(R, P):
    """(index into P's times, τ_k, entrant) for each entry."""
    T = P['taus']
    return [(k, float(T[k]), int(j)) for k in range(1, len(T)) for t, j, kind, c in R['events'] if kind == 'enter' and abs(t - T[k]) < 1e-12]


def b_job(args):
    """The constants of one target: per entry v, g′, the slide, the gap, s_θ (+ the spread of its extrapolation), the pulls on every kernel
    inactive after the entry (extrapolated), the shift σ_k and the recursion's first-order deviations dev (absolute units: a dual may
    change sign between entries; entry k arrives shifted by σ_k = −s_k·dev_{j_k}/v_k, each tracer j gains s_j^{(k)}Δ_kj and
    σ_k(a_{k−1,j} − a_{k,j}), a_{k,j} = ⟨r_kφ(c_j⁰)⟩ the path's rates), λ_k, and (H3)."""
    n, cs, As = args
    t0 = time.time()
    ok, mu0 = eligible(cs, As)
    L, K, c0, f, R, P = path_of(cs, As)
    d, nx, T = c0.size, K.xs.size, P['taus']
    ents = entries(R, P)
    never = [j for j in range(d) if j not in P['ever']]
    dev = {jj: 0.0 for jj in range(d)}
    out = dict(n=n, cs=list(cs), As=list(As), ok=ok, mu0=mu0, entries=[(e[1], e[2]) for e in ents], v=[], g=[], slide=[], gap=[], theta={},
               scross=[], s_spread=[], pulls=[], p_spread=[], shifts=[], devs=[], lam=[], h3w=[], h3c=[], stationary=[])
    for q, (k, tk, j) in enumerate(ents):
        tr = never + [jj for (_, _, jj) in ents[q + 1:]]
        wk, ck = R['hist'][k][2], R['hist'][k][1]
        th = 0.5 * abs(float(wk[j]))
        out['theta'][j] = th
        r_pre = f - K.phi(R['hist'][k - 1][1]) @ R['hist'][k - 1][2]
        r_post = f - K.phi(ck) @ wk
        v = abs(float(K.phi(np.array([c0[j]]))[:, 0] @ r_pre / nx))
        out['v'].append(v)
        out['g'].append(abs(float(K.dphi(np.array([c0[j]]))[:, 0] @ r_pre / nx)))
        out['slide'].append(abs(float(ck[j] - c0[j])))
        out['gap'].append(float(T[k + 1] - T[k]) if k + 1 < len(T) else np.inf)
        rr = [RT.entry_pull(L, K, f, R, P, k, j, tr, xi0=x0, cross=th, rtol=PULL_RTOL, method=PULL_METHOD, deviation=True) for x0 in XI0]
        ext = lambda a, b, c: (a - 6 * b + 8 * c) / 3
        p_hi = {jj: ext(*[r_['pull'][jj] for r_ in rr[1:]]) for jj in tr}
        p_lo = {jj: ext(*[r_['pull'][jj] for r_ in rr[:3]]) for jj in tr}
        out['pulls'].append(p_hi)
        out['p_spread'].append(max(abs(p_hi[jj] - p_lo[jj]) for jj in tr))
        out['scross'].append(ext(*[r_['s_cross'] for r_ in rr[1:]]))
        out['s_spread'].append(abs(out['scross'][-1] - ext(*[r_['s_cross'] for r_ in rr[:3]])))
        out['stationary'].append(all(r_['stationary'] for r_ in rr))
        blk = sorted(rr[-1]['w_end'])
        out['h3w'].append(max(abs(rr[-1]['w_end'][i] - wk[i]) for i in blk))
        out['h3c'].append(max(abs(rr[-1]['c_end'][i] - ck[i]) for i in blk))
        a = np.flatnonzero(wk)
        out['lam'].append(RT.block_rate(L, K, f, wk[a], ck[a]))
        sk = -(np.sign(P['U'][k][j]) or 1.0) * dev[j] / v
        out['shifts'].append(sk)
        for jj in tr:
            sj = np.sign(P['U'][k][jj]) or 1.0
            a0 = float(K.phi(np.array([c0[jj]]))[:, 0] @ r_pre / nx)
            a1 = float(K.phi(np.array([c0[jj]]))[:, 0] @ r_post / nx)
            dev[jj] += sj * p_hi[jj] + sk * (a0 - a1)
        out['devs'].append(dict(dev))
    out['secs'] = time.time() - t0
    print(f"  B target {n}: v {min(out['v']):.3g}, slide {max(out['slide']):.3g}, s_θ+σ "
          f"{max(abs(s + g) for s, g in zip(out['scross'], out['shifts'])):.4g} [{out['secs']:.0f}s]", flush=True)
    return out


def a_job(args):
    """The loop against the path and its first-order correction on one target, at every ε (the constants C from `b_job`)."""
    C = args
    n, cs, As, mu0 = C['n'], np.array(C['cs']), np.array(C['As']), C['mu0']
    t0 = time.time()
    L, K, c0, f, R, P = path_of(cs, As)
    d, T = c0.size, P['taus']
    ents = entries(R, P)
    never = [j for j in range(d) if j not in P['ever']]
    tk = [e[1] for e in ents]
    mids = [0.5 * (tk[q] + tk[q + 1]) for q in range(len(tk) - 1)] + [tk[-1] + 0.5 * (tk[-1] - (tk[-2] if len(tk) > 1 else 0.0))]
    t_end = 1.3 * tk[-1]
    base = sorted(set(mids + [t_end]))
    Cn = {j: RT.creep_constants(L, K, P, j, c0[j], base) for j in never}
    pred = []                                                          # the path at the midpoints
    for q, tm in enumerate(mids):
        k = ents[q][0]
        Uj = P['U'][k] + (tm - T[k]) * (K.phi(P['centres'][k]).T @ P['res'][k] / K.xs.size)
        pred.append(dict(tau=tm, k=k, act=[e[2] for e in ents[:q + 1]], U=Uj,
                         sum_pull={j: (np.sign(Uj[j]) or 1.0) * C['devs'][q][j] for j in never},
                         pulls_only={j: (np.sign(Uj[j]) or 1.0) * sum((np.sign(P['U'][ents[l][0]][j]) or 1.0) * C['pulls'][l][j]
                                                                      for l in range(q + 1)) for j in never}))
    rows = []
    for eps in EPS_A:
        t1 = time.time()
        se = np.sqrt(eps)
        win = []                                                       # (b)'s window: (q, τ, 'start'/'end')
        for q in range(len(tk)):
            a_ = tk[q] + C_B * se * np.log(1 / eps)
            b_ = tk[q + 1] - eps ** 0.25 if q + 1 < len(tk) else t_end
            if a_ < b_:
                win += [(q, a_, 'start')] + ([(q, b_, 'end')] if q + 1 < len(tk) else [])
        grid = np.array(sorted(set(base + [w[1] for w in win])))
        lp = RT.loop_with_duals(L, K, f, c0, eps, grid, levels=C['theta'], P=P, rtol=loop_rtol(eps))
        ix = {float(t): i for i, t in enumerate(grid)}

        def off(i, q):
            k, a = ents[q][0], [e[2] for e in ents[:q + 1]]
            return max(np.abs(lp['w'][i][a] - R['hist'][k][2][a]).max(), np.abs(lp['c'][i][a] - R['hist'][k][1][a]).max())
        blk = [off(ix[float(p['tau'])], q) for q, p in enumerate(pred)]
        wstart = [off(ix[float(t)], q) for q, t, kind in win if kind == 'start']
        wend = [off(ix[float(t)], q) for q, t, kind in win if kind == 'end']
        dev = []
        for q, p in enumerate(pred):
            i = ix[float(p['tau'])]
            for j in never:
                s = np.sign(p['U'][j]) or 1.0
                D = Cn[j][1][base.index(p['tau'])]
                dev.append((j, q, float(s * lp['du'][i][j] / se), float(p['sum_pull'][j]), float(s * eps * D / se), float(p['pulls_only'][j])))
        ie = ix[float(t_end)]
        wn = float(np.abs(lp['w'][ie][never]).max())
        disp = np.array([lp['c'][ie][j] - c0[j] for j in never])
        Cend = np.array([Cn[j][0][-1] for j in never])
        cr = [(j, (lp['cross'][j] - tk_) / se if lp['cross'][j] is not None else np.nan, C['scross'][q] + C['shifts'][q])
              for q, (k, tk_, j) in enumerate(ents)]
        n_form = int(np.sum(np.abs(lp['w'][ie]) > X.TH_ACT * np.abs(As).max()))
        rows.append(dict(eps=eps, block=float(max(blk)), wstart=float(max(wstart)) if wstart else np.nan,
                         wend=float(max(wend)) if wend else np.nan, w_never=wn, disp=disp.tolist(), creep=(eps * Cend).tolist(), dev=dev,
                         cross=cr, n_form=n_form, secs=time.time() - t1))
        print(f"  A target {n} eps {eps:g}: formed {n_form}, block {max(blk):.2e}, window {rows[-1]['wstart']:.2e}/{rows[-1]['wend']:.2e}, "
              f"w_never {wn:.2e} [{time.time() - t1:.0f}s]", flush=True)
    return dict(n=n, mu0=mu0, rows=rows, secs=time.time() - t0)


def part_B():
    from multiprocessing import Pool
    tg = TM.random_targets(SEED_A, N_POOL)
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(n, cs, As) for n, (cs, As) in enumerate(tg)], chunksize=1)


def chosen(B):
    """The N_A eligible targets with the largest margins, then the stress targets: the slowest entrant, the largest slide, the longest
    inner layer against its gap (|s_θ + σ|/gap), each the first not already chosen."""
    el = [b for b in B if b['ok']]
    top = [b['n'] for b in sorted(el, key=lambda b: -b['mu0'])[:N_A]]
    out = list(top)
    for key in (lambda b: min(b['v']), lambda b: -max(b['slide']),
                lambda b: -max(abs(s + g) / gp for s, g, gp in zip(b['scross'], b['shifts'], b['gap']))):
        for b in sorted(el, key=key):
            if b['n'] not in out:
                out.append(b['n'])
                break
    return out


def part_A(B):
    from multiprocessing import Pool
    sel = chosen(B)
    print(f"selected {sel}", flush=True)
    by = {b['n']: b for b in B}
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, [by[n] for n in sel], chunksize=1)


def slope(x, y):
    x, y = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    return float(np.polyfit(x, y, 1)[0]) if x.size >= 2 else np.nan


def report(B=None, A=None):
    fmt = E.fmt
    L_ = ['# Theorem D‴.53: the sliding limit (`exp_transit_limit.py`)', '',
          "The free law with moving centres (m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025) from α = 0 with equal masses ε, Gaussian kernels "
          "of width 0.12 on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, against its sliding path (polished); the "
          f"targets: D‴.51's {N_POOL} out-of-sample three-bump targets (seed {SEED_A}).", '']
    if B is not None:
        nel = sum(b['ok'] for b in B)
        L_ += ['## B. The constants of every target', '',
               f"{nel} of {len(B)} targets have an exact path with simple entries, no exits and margins ≥ {MU0}. Per entry (in the order of "
               "the entries): the approach speed v_k, its sensitivity to the entrant's centre g′_k, the slide |c^k − c⁰|, the gap to the "
               "next entry, s_θ + σ_k (the inner layer's crossing time with the entry's shift, in units of ε^{1/2}), the largest pull "
               "|Δ|, the block's slowest rate λ_k (inner time), and (H3): the inner problem's end against S_k (coefficients, centres). "
               f"s_θ and the pulls by second-order Richardson in xi0 from {XI0[1]:g}, {XI0[2]:g}, {XI0[3]:g} (`entry_pull`'s deviation form, "
               f"Radau, rtol {PULL_RTOL:g}); in brackets the spread against the same from {XI0[0]:g}, {XI0[1]:g}, {XI0[2]:g}.", '',
               '| target (margin) | v_k | g′_k | slide | gap | s_θ + σ_k [spread] | max \\|Δ\\| [spread] | λ_k | (H3) w, c |',
               '|---|---|---|---|---|---|---|---|---|']
        for b in B:
            L_.append(f"| {b['n']} ({b['mu0']:.3f}{'' if b['ok'] else ', excluded'}) | " + ', '.join(f'{x:.3g}' for x in b['v']) + ' | '
                      + ', '.join(f'{x:.3g}' for x in b['g']) + ' | ' + ', '.join(f'{x:.3g}' for x in b['slide']) + ' | '
                      + ', '.join('∞' if not np.isfinite(x) else f'{x:.3g}' for x in b['gap']) + ' | '
                      + ', '.join(f'{s + g:.5g} [{sp:.1e}]' for s, g, sp in zip(b['scross'], b['shifts'], b['s_spread'])) + ' | '
                      + ', '.join(f'{max(abs(x) for x in p.values()):.4g} [{sp:.1e}]' for p, sp in zip(b['pulls'], b['p_spread'])) + ' | '
                      + ', '.join(f'{x:.3g}' for x in b['lam']) + ' | '
                      + ', '.join(f'{fmt(w)}, {fmt(c)}' for w, c in zip(b['h3w'], b['h3c'])) + ' |')
        allv = [x for b in B for x in b['v']]
        alls = [x for b in B for x in b['slide']]
        sth = [abs(s + g) for b in B for s, g in zip(b['scross'], b['shifts'])]
        mp = [max(abs(x) for x in p.values()) for b in B for p in b['pulls']]
        L_ += ['', f"Over the {len(B)} targets: v_k {min(allv):.3g}–{max(allv):.3g}, slides {min(alls):.3g}–{max(alls):.3g}, |s_θ + σ_k| "
                   f"{min(sth):.4g}–{max(sth):.4g}, the largest pull per entry {min(mp):.3g}–{max(mp):.4g}; the extrapolation's spread at "
                   f"most {max(x for b in B for x in b['s_spread']):.1e} in s_θ and {max(x for b in B for x in b['p_spread']):.1e} in the "
                   f"pulls; λ_k {min(x for b in B for x in b['lam']):.3g}–{max(x for b in B for x in b['lam']):.3g}; (H3) to "
                   f"{max(x for b in B for x in b['h3w']):.1e} (coefficients), {max(x for b in B for x in b['h3c']):.1e} (centres); every "
                   f"inner problem stationary: {all(all(b['stationary']) for b in B)}. Part A's targets: {chosen(B)}.", '']
    if A is not None and all('wstart' in r for t in A for r in t['rows']):
        L_ += ['## A. The rates', '',
               "Per target and ε: the kernels formed; the block against S_k at the midpoints between entries and at (b)'s window's start "
               f"τ_k + {C_B:g}ε^{{1/2}}log(1/ε) and end τ_{{k+1}} − ε^{{1/4}} (max over the coefficients and centres); the never-entering "
               "kernels' largest coefficient (× ε^{−3/2}); their centres' largest displacement against the creep's εC_j (the largest ratio's "
               "distance from 1); the duals: obs = s_j(u_j − u_j^path)/ε^{1/2} against the first-order deviation — the entries' pulls Δ_kj "
               "plus each entry's shift σ_k = −dev_{j_k}/v_k times the drop of j's rate at it, by recursion over the entries — plus the "
               "creep's ε^{1/2}D_j (the largest |obs − pred|, the median ratio obs/pred where |pred| > 0.1, and the largest |obs − ΣΔ| with "
               "the pulls alone); the entries' crossing times (τ_θ − τ_k)/ε^{1/2} against s_θ + σ_k (the largest difference). The loop at "
               "rtol 10⁻¹¹ (ε ≥ 10⁻⁸), 10⁻¹² (10⁻⁹), 10⁻¹³ (≤ 10⁻¹⁰).", '',
               '| target (margin) | ε | formed | block − S_k: midpoints; window start, end | max \\|w_never\\| (× ε^{−3/2}) | displacement / εC − 1 '
               '| duals: max \\|obs − pred\\| (median obs/pred) | the pulls alone | crossings: max \\|obs − pred\\| |',
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
                L_.append(f"| {t['n']} ({t['mu0']:.3f}) | {e:g} | {r['n_form']} | {fmt(r['block'])}; {fmt(r['wstart'])}, {fmt(r['wend'])} | "
                          f"{fmt(r['w_never'])} ({r['w_never'] / e ** 1.5:.3g}) | {rel:.3g} | {np.abs(obs - pr).max():.4g} ({med:.4f}) | "
                          f"{np.abs(obs - po).max():.4g} | {max(xs) if xs else np.nan:.4g} |")
        L_ += ['', 'Slopes in log ε (least squares over ε ≤ 10⁻⁸): ' + '; '.join(
            f"target {t['n']}: block {slope(*zip(*[(r['eps'], r['block']) for r in t['rows'] if r['eps'] <= 1e-8])):.3f}, w_never "
            f"{slope(*zip(*[(r['eps'], r['w_never']) for r in t['rows'] if r['eps'] <= 1e-8])):.3f}, displacement "
            f"{slope(*zip(*[(r['eps'], np.abs(r['disp']).max()) for r in t['rows'] if r['eps'] <= 1e-8])):.3f}, duals' residual "
            f"{slope(*zip(*[(r['eps'], max(abs(d[2] - d[3] - d[4]) for d in r['dev'])) for r in t['rows'] if r['eps'] <= 1e-8])):.3f}"
            for t in A) + ' (predicted 1.5, 1.5, 1 and, for the residual in units of ε^{1/2}, ≥ 0.25).', '']
        L_ += ['The duals\' residual per decade (the ratio of successive largest |obs − pred|): ' + '; '.join(
            f"target {t['n']}: " + ', '.join(f"{a / b:.3g}" for a, b in zip(
                [max(abs(d[2] - d[3] - d[4]) for d in r['dev']) for r in t['rows'][:-1]],
                [max(abs(d[2] - d[3] - d[4]) for d in r['dev']) for r in t['rows'][1:]])) for t in A) + '.', '']
    L_.append('Numbers are produced by `experiments/exp_transit_limit.py`; nothing is edited by hand.')
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'transit_limit_raw.npz')
    md = os.path.join(OUT, 'transit_limit.md')
    old = E.load(raw)
    keys = ('B', 'A')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    if 'B' in parts:
        res['B'] = part_B(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] B", flush=True)
        open(md, 'w').write(report(*[res.get(k) for k in keys]))
    if 'A' in parts:
        res['A'] = part_A(res['B']); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] A", flush=True)
        open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
