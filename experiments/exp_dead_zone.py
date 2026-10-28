"""
Proposition D‴.57: the transit's dead zone — gated by the residual's energy the free law stops at the discrepancy point; gating every
particle freezes the last block before it has relaxed, gating only the unformed ones gives the oracle's fit
====================================================================================================================================

D‴.56's setting (21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]; three features at least 0.45 apart
with m > 0.8q; data y = f* + ξ, ξ i.i.d. N(0, s²)), the gate open while the residual's rms exceeds (1 + δ)s:

A. The sliding limit: the noisy path; the dead zone stops it at the first entry k whose relaxed residual's rms is ≤ (1 + δ)s (the
   discrepancy rule, D‴.56 D). With every particle gated (the energy dead zone) the k-th block's fast flow from its entrant's seed
   (`critical_rhs` on the block) stops where the rms falls to (1 + δ)s: the frozen fit's error against the oracle's (the k-th block
   relaxed, the k-kernel least-squares fit of y) and against s(4K/N + 2δ + δ²)^{1/2}; the frozen residual's tangent part u = r_c − r_K,
   |u|² against (1 + δ)²s²N − |r_K|², the cross term ⟨f̂_K − f*, u⟩ (mean zero predicted) and u's direction across a target's draws (the
   block's slowest mode). With only the unformed particles gated (the formation dead zone) the block relaxes: the oracle's fit.
B. The loop at the mass EPS_B from α = 0 on N_B targets × N_DRAW_B draws at LEVELS_B: ungated to τ_K + 3W (W the path's window to its
   first noise kernel), with the energy dead zone and with the formation dead zone (|w_i| ≤ THETA_B gated): the kernels formed (|w| >
   THETA_B) and the fits' errors against the sliding limit's.

C. The closed loop of Proposition D‴.49 under its measurement noise (Part E there: 16 tones in 10–30 rad/s of rms s, ε = EPS_C, to
   τ = TAU_C), the regressor at the measured state or at the reference (D‴.31's desired regressor): ungated, with D‴.49's mirror leak
   (λΛ = LEAK_C), and with the dead zone on the error's energy filtered over TF_C, gated between the rms GATE_C·s (every kernel, or
   only the unformed ones): the kernels formed, those pushed beyond the range, the coefficients and the tracking error at the end;
   and the noise-free loop with either regressor as the reference.

Usage: python experiments/exp_dead_zone.py [--only A,B,C] [--from-raw] → results/chaos/dead_zone.md (+ _raw.npz).
"""
import os
import sys
import time
import warnings

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_noisy_certificate as NC
import exp_transit_mass as TM

OUT = os.environ.get('DZ_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('DZ_WORKERS', '2'))
DELTA = 0.1
LEVELS_A, N_DRAW_A = (0.01, 0.03, 0.1), 10
EPS_B, N_B, N_DRAW_B, LEVELS_B, THETA_B = 1e-8, 16, 3, (0.01, 0.03, 0.1), 1e-3
T_RELAX = 1e8                                                          # the block's fast flow when the gate does not close
EPS_C, TAU_C, SIG_C, LEAK_C, TF_C, GATE_C = 1e-4, 600.0, (0.01, 0.02), 3e-3, 5.0, (1.2, 1.5)


def gated_relax(L, K, y, c_blk, w_blk, thr):
    """The block's fast flow on the critical curve from (c_blk, w_blk) until the residual's rms falls to thr (event), or T_RELAX.
    Returns (crossed, c, w)."""
    n = len(c_blk)
    al = np.sign(w_blk) * RT._alpha_of(L, np.abs(w_blk), RT.coef_crit)
    rhs = RT.critical_rhs(L, K, y)

    def ev(t, z):
        return float(np.sqrt(np.mean((y - K.phi(z[n:]) @ RT.critical_w(L, z[:n])) ** 2))) - thr
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, (0, T_RELAX), np.r_[al, c_blk], method='LSODA', rtol=1e-10, atol=1e-13, events=ev)
    z = sol.y[:, -1]
    return bool(sol.t_events[0].size), z[n:], RT.critical_w(L, z[:n])


def a_draw(L, K, c0, f, E_, xi, s, nK):
    """One draw: the stop, the energy dead zone's frozen fit and the oracle's, the tangent part u."""
    y = f + xi
    thr = (1 + DELTA) * s
    # to the (K + 1)-th entry, not relaxed (a noise kernel's relaxation is slow and not needed: a stop after the K-th counts as late)
    R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), polish=True, max_events=nK + 3, tol=0.0, max_entries=nK + 1)
    ent = [e for e in R['events'] if e[2] == 'enter']
    exits = sum(e[2] == 'exit' for e in R['events'])
    exact = sorted(e[1] for e in ent[:nK]) == E_ and exits == 0
    d = dict(exact=exact, exits=exits)
    if not exact:
        return d
    rms = [h[3] for h in R['hist'][1:nK + 1]]                          # relaxed after entries 1, …, K
    k = next((i + 1 for i, r_ in enumerate(rms) if r_ <= thr), -1)
    d['k'] = k
    if k < 1 or k > len(ent):
        return d
    # the k-th block before it relaxes: the (k − 1)-th block relaxed and the k-th entrant at its seed
    tau, c_prev, w_prev, _ = R['hist'][k - 1]
    j = ent[k - 1][1]
    r_prev = y - K.phi(c_prev) @ w_prev
    sg = np.sign(K.phi(c_prev)[:, j] @ r_prev)
    act = np.flatnonzero(w_prev != 0)
    blk = np.r_[act, j]
    w0 = np.r_[w_prev[act], sg * 1e-7]
    crossed, cc, wc = gated_relax(L, K, y, c_prev[blk], w0, thr)
    _, cK, wK, rmsK = R['hist'][k]                                     # the k-th block relaxed (polished): the oracle at k = K
    fK, fc = K.phi(cK) @ wK, K.phi(cc) @ wc
    u = fK - fc                                                        # r_c − r_K
    N = y.size
    d.update(crossed=crossed, err_dz=float(np.sqrt(np.mean((fc - f) ** 2))), err_or=float(np.sqrt(np.mean((fK - f) ** 2))),
             u2=float(u @ u), u2_pred=float(thr ** 2 * N - np.sum((y - fK) ** 2)),
             cross=float((fK - f) @ u / (np.linalg.norm(fK - f) * np.linalg.norm(u) + 1e-300)), u=u)
    return d


def a_job(args):
    q, cs, As = args
    t0 = time.time()
    warnings.filterwarnings('ignore')
    L, K, c0, f, cert, R0, P = NC.noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    E_ = sorted(cert['e'])
    rng = np.random.default_rng(5700 + q)
    out = dict(q=q, levels=[])
    for s in LEVELS_A:
        draws = [a_draw(L, K, c0, f, E_, s * rng.standard_normal(nx), s, nK) for _ in range(N_DRAW_A)]
        # u's direction across the draws (the block's slowest mode): |cos| against the first draw's
        us = [d['u'] for d in draws if d.get('k') == nK and 'u' in d]
        for d in draws:
            if 'u' in d:
                d['cos'] = float(abs(d['u'] @ us[0]) / (np.linalg.norm(d['u']) * np.linalg.norm(us[0]) + 1e-300)) if us else np.nan
                del d['u']
        out['levels'].append(dict(s=s, draws=draws))
    out['secs'] = time.time() - t0
    print(f"  A target {q}: " + ', '.join(f"{sum(d.get('k') == nK for d in lv['draws'])}/{len(lv['draws'])}" for lv in out['levels'])
          + f" [{out['secs']:.0f}s]", flush=True)
    return out


def part_A():
    from multiprocessing import Pool
    jobs = [(q, cs, As) for q, (cs, As) in enumerate(NC.targets())]
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, jobs, chunksize=1)


def loop_gated(L, K, y, c0, eps, tau_max, thr, mode):
    """The loop from α = 0 at the masses ε to τ_max (τ = t/Λ(ε)): mode 'none' ungated; 'energy' every particle's drive gated by the
    residual's rms (frozen once it falls to thr); 'formation' the unformed particles (|w| ≤ THETA_B) frozen at that moment, the formed
    block relaxing to τ_max. Returns (w, c, τ at which the gate closed or nan)."""
    d = len(c0)
    b0 = eps * np.ones(d)
    lam = RT.transit_scale(L, eps)
    rhs = RT.critical_rhs(L, K, y, b0=b0)
    z0 = np.r_[np.zeros(d), c0]
    wz = lambda z: RT.critical_w(L, z[:d], b0)
    if mode == 'none':
        sol = solve_ivp(rhs, (0, lam * tau_max), z0, method='LSODA', rtol=1e-9, atol=1e-12)
        z = sol.y[:, -1]
        return wz(z), z[d:], np.nan

    def ev(t, z):
        return float(np.sqrt(np.mean((y - K.phi(z[d:]) @ wz(z)) ** 2))) - thr
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(rhs, (0, lam * tau_max), z0, method='LSODA', rtol=1e-9, atol=1e-12, events=ev)
    z = sol.y[:, -1]
    if not sol.t_events[0].size:
        return wz(z), z[d:], np.nan
    t_c = float(sol.t_events[0][0])
    if mode == 'formation':
        frozen = np.abs(wz(z)) <= THETA_B
        mask = np.r_[~frozen, ~frozen].astype(float)
        sol2 = solve_ivp(lambda t, zz: mask * rhs(t, zz), (t_c, max(lam * tau_max, t_c * 1.5)), z, method='LSODA', rtol=1e-9, atol=1e-12)
        z = sol2.y[:, -1]
    return wz(z), z[d:], t_c / lam


def b_job(args):
    q, cs, As, s, k_ = args
    t0 = time.time()
    warnings.filterwarnings('ignore')
    L, K, c0, f, cert, R0, P = NC.noiseless(cs, As)
    nK, nx = len(cs), K.xs.size
    xi = s * np.random.default_rng(5800 + 97 * q + k_).standard_normal(nx)
    y = f + xi
    thr = (1 + DELTA) * s
    R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), polish=True, max_events=nK + 4, tol=0.0, max_entries=nK + 1)
    ent = [e for e in R['events'] if e[2] == 'enter']
    out = dict(q=q, s=s, k=k_, exact=sorted(e[1] for e in ent[:nK]) == sorted(cert['e']) and not any(e[2] == 'exit' for e in R['events']))
    if not out['exact'] or len(ent) <= nK:
        return out
    tau_K, W = ent[nK - 1][0], ent[nK][0] - ent[nK - 1][0]
    tau_max = tau_K + 3 * W
    lim = a_draw(L, K, c0, f, sorted(cert['e']), xi, s, nK)            # the sliding limit's stop and fits on the same draw
    out.update(tau_K=tau_K, W=W, k_lim=lim.get('k'), err_dz_lim=lim.get('err_dz', np.nan), err_or_lim=lim.get('err_or', np.nan))
    for mode in ('none', 'energy', 'formation'):
        w, c, tc = loop_gated(L, K, y, c0, EPS_B, tau_max, thr, mode)
        out[mode] = dict(formed=int((np.abs(w) > THETA_B).sum()), err=float(np.sqrt(np.mean((K.phi(c) @ w - f) ** 2))), tau_gate=tc,
                         rms=float(np.sqrt(np.mean((y - K.phi(c) @ w) ** 2))))
    out['secs'] = time.time() - t0
    print(f"  B target {q} s {s:g} draw {k_}: formed " + ', '.join(f"{m} {out[m]['formed']}" for m in ('none', 'energy', 'formation'))
          + f" [{out['secs']:.0f}s]", flush=True)
    return out


def part_B():
    from multiprocessing import Pool
    tg = NC.targets()[:N_B]
    jobs = [(q, cs, As, s, k_) for q, (cs, As) in enumerate(tg) for s in LEVELS_B for k_ in range(N_DRAW_B)]
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, jobs, chunksize=1)


def c_job(args):
    """C: one closed loop (D‴.49's setting, the C integrator) — mode 'none', 'leak', 'energy' or 'formation', the regressor at x̃ or x_m."""
    s, regx, mode = args
    import exp_closed_transit as CT
    from src.core import closed_loop as CL
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = CT.law()
    Lam = RT.transit_scale(L, EPS_C)
    T = Lam * TAU_C
    per = 2 * np.pi / CT.OM
    ts = np.unique(np.r_[np.linspace(0, T, 201), np.linspace(T - per, T, 60)])
    th = CT.th_abs()
    gate = (dict(mode=mode, Tf=TF_C, lo=GATE_C[0] * s, hi=GATE_C[1] * s, theta=th) if mode in ('energy', 'formation') and s > 0 else None)
    r = CL.simulate(L, EPS_C, CT.C0, CT.SIG, (CT.CSTAR, CT.ASTAR, CT.SIG), CT.KE, (0.0, 1.0, CT.OM, 0.0, np.sqrt(2) * CT.OM), ts,
                    noise=CL.tones(s, 10.0, 30.0, seed=7) if s > 0 else None, th=th, regx=regx,
                    lam=(LEAK_C / Lam if mode == 'leak' else 0.0), gate=gate)
    act = np.abs(r['w'][-1]) > th
    closed = float(ts[np.argmax(r['z'] < (GATE_C[0] * s) ** 2)] / Lam) if gate is not None and (r['z'] < (GATE_C[0] * s) ** 2).any() else np.nan
    xr = np.linspace(-1, 1, 201)                                       # the reference's range
    phi = lambda c: np.exp(-(xr[:, None] - np.asarray(c)[None, :]) ** 2 / (2 * CT.SIG ** 2))
    fit_rms = float(np.sqrt(np.mean((phi(r['c'][-1]) @ r['w'][-1] - phi(CT.CSTAR) @ CT.ASTAR) ** 2)))
    out = dict(s=s, regx=bool(regx), mode=mode, n_act=int(act.sum()), outside=int(np.sum(np.abs(r['c'][-1][act]) > CT.OUTSIDE)),
               w=r['w'][-1][act].tolist(), c=r['c'][-1][act].tolist(), wmax=float(np.abs(r['w'][-1]).max()), e_last=r['e_last'],
               e_true=r['e_true_rms_last'], e_true_max=r['e_true_last'], fit_rms=fit_rms, w_all=r['w'][-1].tolist(),
               c_all=r['c'][-1].tolist(), tau_closed=closed, secs=time.time() - t0)
    print(f"  C s {s:g} regx {regx} {mode}: {out['n_act']} kernels, {out['outside']} outside, max|w| {out['wmax']:.3g}, "
          f"max|e| {out['e_last']:.4f}, true error rms {out['e_true']:.4f}, fit rms {fit_rms:.4f} [{out['secs']:.0f}s]", flush=True)
    return out


def part_C():
    from multiprocessing import Pool
    jobs = ([(0.0, regx, 'none') for regx in (False, True)]                # the noise-free references
            + [(s, regx, mode) for regx in (False, True) for s in SIG_C for mode in ('none', 'leak', 'energy', 'formation')])
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, jobs, chunksize=1)


def report(A=None, B=None, C=None):
    nK, N = 3, 61
    L_ = ["# Proposition D‴.57: the transit's dead zone (`exp_dead_zone.py`)", "",
          "D‴.56's setting — 21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]; three features at least "
          "0.45 apart with m > 0.8q; y = f* + ξ, ξ i.i.d. N(0, s²) —, the gate open while the residual's rms exceeds (1 + δ)s, "
          f"δ = {DELTA:g}.", ""]
    if A is not None:
        L_ += ["## A. The sliding limit", "",
               f"Per level s: {N_DRAW_A} draws on each of the {len(A)} targets; on the draws whose noisy path enters exactly E first: the "
               "dead zone's stop (the first entry whose relaxed rms is ≤ (1 + δ)s) at the K-th entry, later or earlier; at the K-th, the "
               "energy dead zone's frozen fit against the oracle's (the K-th block relaxed): the errors' ratio, the mean square over its "
               "prediction given the stop at the K-th entry (± its standard error) and over s²(4K/N + 2δ + δ²) (the unconditional value), "
               "the oracle's over s²·2K/N, the unrelaxed tangent part |u|² over (1 + δ)²s²N − |r_K|², the mean |u|²/(s²N) over its "
               "prediction given the stop, the correlation of the oracle's error with u and |cos| of u against the target's first draw's:",
               "", "| s | exact draws | stop at K, later, earlier | err_dz/err_or: median [5–95%] | mean err_dz²: / prediction given the "
               "stop; / s²(4K/N + 2δ + δ²) | mean err_or² / (s²2K/N) | \\|u\\|² / (1 + δ)²s²N − \\|r_K\\|²: median [5–95%] | mean "
               "\\|u\\|²/(s²N) / prediction given the stop | corr(err_or, u): mean ± s.e. | \\|cos\\| of u: median |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        from scipy.stats import chi2
        c_ = (1 + DELTA) ** 2 * N
        p_K = chi2.cdf(c_, N - 2 * nK)
        m_c = (N - 2 * nK) * chi2.cdf(c_, N - 2 * nK + 2) / p_K             # E[χ²_{N−2K} | χ²_{N−2K} ≤ (1 + δ)²N]
        u_c = (c_ - m_c) / N                                                # E|u|²/(s²N) given the stop at K
        dz_c = (2 * nK + c_ - m_c) / N                                      # E err_dz²/s² given the stop at K
        for i, s in enumerate(LEVELS_A):
            ds = [d for t in A for d in t['levels'][i]['draws'] if d['exact']]
            atK = [d for d in ds if d.get('k') == nK and d.get('crossed')]
            later = sum(1 for d in ds if d.get('k', -1) > nK or d.get('k') == -1)
            earlier = sum(1 for d in ds if 0 < d.get('k', 0) < nK)
            if atK:
                rat = np.array([d['err_dz'] / d['err_or'] for d in atK])
                e2 = np.array([d['err_dz'] ** 2 for d in atK]) / s ** 2
                e_dz = e2.mean() / (4 * nK / N + 2 * DELTA + DELTA ** 2)
                e_or = np.mean([d['err_or'] ** 2 for d in atK]) / (s ** 2 * 2 * nK / N)
                u_r = np.array([d['u2'] / d['u2_pred'] for d in atK if d['u2_pred'] > 0])
                cr = np.array([d['cross'] for d in atK])
                cs_ = np.array([d['cos'] for d in atK if np.isfinite(d.get('cos', np.nan))])
                L_.append(f"| {s:g} | {len(ds)} | {len(atK)}, {later}, {earlier} | {np.median(rat):.3f} [{np.quantile(rat, 0.05):.3f}, "
                          f"{np.quantile(rat, 0.95):.3f}] | {e2.mean() / dz_c:.3f} ± {e2.std() / np.sqrt(e2.size) / dz_c:.3f}; {e_dz:.3f} | "
                          f"{e_or:.3f} | {np.median(u_r):.4f} [{np.quantile(u_r, 0.05):.4f}, {np.quantile(u_r, 0.95):.4f}] | "
                          f"{np.mean([d['u2'] / N for d in atK]) / s ** 2 / u_c:.3f} | "
                          f"{cr.mean():+.3f} ± {cr.std() / np.sqrt(cr.size):.3f} | {np.median(cs_):.4f} |")
            else:
                L_.append(f"| {s:g} | {len(ds)} | 0, {later}, {earlier} | — | — | — | — | — | — | — |")
        L_ += ["", f"At first order |u|²/N = (1 + δ)²s² − |r_K|²/N, of mean s²(2δ + δ² + 2K/N): |u|/(s√N) ≈ "
                   f"{np.sqrt(2 * DELTA + DELTA ** 2 + 2 * nK / N):.3f} against the oracle's s(2K/N)^{{1/2}} = {np.sqrt(2 * nK / N):.3f}s; without the "
                   f"cross term the errors' rms ratio would be ((4K/N + 2δ + δ²)/(2K/N))^{{1/2}} = "
                   f"{np.sqrt((4 * nK / N + 2 * DELTA + DELTA ** 2) / (2 * nK / N)):.3f}. The stop at the K-th entry is the event |r_K|² ≤ "
                   f"(1 + δ)²s²N, |r_K|²/s² ~ χ²_{{N−2K}} (probability {p_K:.3f}); given it E|r_K|² = s²(N − 2K)P(χ²_{{N−2K+2}} ≤ c)/"
                   f"P(χ²_{{N−2K}} ≤ c), c = (1 + δ)²N, = {m_c:.2f}s² against {N - 2 * nK}s², so E|u|²/N = {u_c:.4f}s² and E err_dz² = "
                   f"s²(2K + c − E[χ²_{{N−2K}} | ≤ c])/N = {dz_c:.4f}s² (×{dz_c / ((4 * nK / N + 2 * DELTA + DELTA ** 2)):.3f} the "
                   f"unconditional value; the rms ratio {np.sqrt(dz_c / (2 * nK / N)):.3f}).", ""]
    if B is not None:
        rs = [r for r in B if r.get('exact') and 'none' in r]
        L_ += ["## B. The loop", "",
               f"The loop at ε = {EPS_B:g} from α = 0 on {N_B} targets × {N_DRAW_B} draws at s = {', '.join(f'{s:g}' for s in LEVELS_B)} "
               "(the draws whose noisy path is exactly E), to τ_K + 3W: on the draws whose sliding limit stops at the K-th entry, the kernels "
               f"formed (|w| > {THETA_B:g}) and the fits' rms errors over the oracle's (the limit's K-th block relaxed), the energy dead zone's "
               "against the ungated loop's and the gated loops against the sliding limit's fits (the energy dead zone's frozen fit, the "
               "oracle's); on the late draws (the stop after the K-th entry) the kernels formed:", "",
               "| s | draws: stop at K, late | formed at K: ungated, energy, formation (median [range]) | err / err_or (medians): ungated, "
               "energy, formation | energy against the ungated (median) | the loop against the limit: max \\|err − err_dz\\| (energy), "
               "max \\|err − err_or\\| (formation) | the gate's τ / τ_K (median) | formed, late: ungated, energy, formation |",
               "|---|---|---|---|---|---|---|---|"]
        for s in LEVELS_B:
            ss = [r for r in rs if r['s'] == s]
            atK = [r for r in ss if r.get('k_lim') == nK]
            late = [r for r in ss if r.get('k_lim') != nK]
            if not atK:
                continue
            fm = {m: [r[m]['formed'] for r in atK] for m in ('none', 'energy', 'formation')}
            med = lambda v: f"{np.median(v):.3f}"
            L_.append(f"| {s:g} | {len(atK)}, {len(late)} | " + '; '.join(f"{int(np.median(fm[m]))} [{min(fm[m])}, {max(fm[m])}]" for m in fm)
                      + " | " + ', '.join(med([r[m]['err'] / r['err_or_lim'] for r in atK]) for m in ('none', 'energy', 'formation')) + " | "
                      + med([r['energy']['err'] / r['none']['err'] for r in atK]) + " | "
                      + f"{max(abs(r['energy']['err'] - r['err_dz_lim']) for r in atK):.1e}, "
                      + f"{max(abs(r['formation']['err'] - r['err_or_lim']) for r in atK):.1e} | "
                      + f"{np.median([r['energy']['tau_gate'] / r['tau_K'] for r in atK if np.isfinite(r['energy']['tau_gate'])]):.4f} | "
                      + ('; '.join(', '.join(str(r[m]['formed']) for m in ('none', 'energy', 'formation')) for r in late) or '—') + " |")
        L_.append("")
    if C is not None:
        lab = dict(none='ungated', leak=f'the mirror leak λΛ = {LEAK_C:g}', energy='the energy dead zone', formation='the formation dead zone')
        L_ += ["## C. The closed loop under measurement noise", "",
               f"D‴.49's Part E (k_e = 10, ω = 0.2, three bumps of width 0.12; 21 kernels from α = 0 at ε = {EPS_C:g} to τ = {TAU_C:g}; the "
               "noise 16 tones in 10–30 rad/s of rms s, the same phases); the dead zone on the error's energy filtered over "
               f"{TF_C:g} time units, the drive gated smoothly between the rms {GATE_C[0]:g}s and {GATE_C[1]:g}s; formed: |w| > 1% of max|A*|; "
               "outside: |c| > 1.2; the measured error's largest value over the last period, the true error x − x_m's rms there and the "
               "fit's rms error on the reference's range [−1, 1] (every kernel):", "",
               "| regressor | s | adaptation | formed | outside | the formed kernels' coefficients | their centres | max\\|ẽ\\|, last period | "
               "true error rms, last period | fit rms on [−1, 1] | the gate closed at τ |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(C, key=lambda r: (r['regx'], r['s'], ['none', 'leak', 'energy', 'formation'].index(r['mode']))):
            L_.append(f"| {'x_m' if r['regx'] else 'x̃'} | {r['s']:g} | {lab[r['mode']]} | {r['n_act']} | {r['outside']} | "
                      f"{', '.join(f'{v:.3f}' for v in r['w'])} | {', '.join(f'{v:.3f}' for v in r['c'])} | {r['e_last']:.4f} | "
                      + (f"{r['e_true']:.4f} | {r['fit_rms']:.4f}" if 'fit_rms' in r else '— | —') + " | "
                      + (f"{r['tau_closed']:.1f}" if np.isfinite(r['tau_closed']) else '—') + " |")
        import exp_closed_transit as CT
        from src.core import closed_loop as CL
        nz = CL.tones(1.0, 10.0, 30.0, seed=7)
        fl = float(np.sqrt(np.sum(nz['amp'] ** 2 / 2 * nz['w'] ** 2 / (nz['w'] ** 2 + CT.KE ** 2))))
        L_ += ["", f"Near the fit the loop high-passes the noise in the measured error (ẽ = e + n, e's response to n through −k_eẽ: "
                   f"|iω/(iω + k_e)|), a floor of rms {fl:.3f}s: the energy dead zone starts to close when the slow error's rms falls to "
                   f"{np.sqrt(GATE_C[1] ** 2 - fl ** 2):.2f}s and is shut below {np.sqrt(GATE_C[0] ** 2 - fl ** 2):.2f}s."]
        ref = {r['regx']: r for r in C if r['s'] == 0 and r['mode'] == 'none'}
        dev = []
        for s_ in SIG_C:
            for regx in (False, True):
                r = next((r for r in C if r['s'] == s_ and r['regx'] == regx and r['mode'] == 'none'), None)
                if r is None or regx not in ref:
                    continue
                r0 = ref[regx]
                if r['n_act'] != r0['n_act'] or r['outside'] or r0['outside']:
                    dev.append(f"{'x_m' if regx else 'x̃'}, s = {s_:g}: {r['n_act']} kernels against {r0['n_act']}")
                    continue
                o, o0 = np.argsort(r['c']), np.argsort(r0['c'])
                dev.append(f"{'x_m' if regx else 'x̃'}, s = {s_:g}: max|Δw| {np.abs(np.array(r['w'])[o] - np.array(r0['w'])[o0]).max():.2e}, "
                           f"max|Δc| {np.abs(np.array(r['c'])[o] - np.array(r0['c'])[o0]).max():.2e}")
        if dev:
            L_.append("The ungated loops against the noise-free loop with the same regressor (the formed kernels in the order of their "
                      "centres): " + '; '.join(dev) + ".")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_dead_zone.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'dead_zone_raw.npz')
    md = os.path.join(OUT, 'dead_zone.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(old.get('A'), old.get('B'), old.get('C')))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else {'A', 'B', 'C'}
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(res.get('A'), res.get('B'), res.get('C')))


if __name__ == '__main__':
    main()
