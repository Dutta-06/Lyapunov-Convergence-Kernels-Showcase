"""
Proposition D‴.49: the rich transit in the closed loop — from small masses the adaptive controller forms kernels one at a time
====================================================================================================================================

The scalar plant ẋ = f(x) + u, u = ẋ_m − k_ee − f̂(x), e = x − x_m, x_m = x₀ + a sin(ωt); f̂ = Σ_iw_iφ(·; c_i), Gaussian kernels of width 0.12
started at α = 0 with the masses ε on a grid of 21 centres in [−1, 1]; the paradigm's free law with moving centres driven by the error,
dα_i = Gβ_ia′_ieφ_i(x), dβ_i = γβ_ia_ieφ_i(x), ċ_i = Gβ_ia_ie∂_cφ_i(x) (Proposition D‴.48(iii); `rich_transit`). f = three bumps of the kernels'
width at −0.52, 0.07, 0.61 with amplitudes 0.9, −0.6, 0.75 (max|f′| = 4.55). The real loop is integrated in C (`closed_loop`, the kernels
kept on their curves); the averaged loop along the frozen periodic orbit (`rich_transit.FrozenOrbit`, exact for k_e > Lip(f − f̂)) and
its leading order in 1/k_e, the reference's own orbit seen through the error filter (`OrbitBumps`). Entry times (a kernel formed when
|w| > 10⁻²max|A*|) are interpolated in log|w| between output points; τ = t/Λ(ε).

A. The base point (k_e = 10, ω = 0.2, x_m = sin ωt): the exact averaged loop and its leading order for ε = 10⁻⁴ … 10⁻⁷, the real loop at
   ε = 10⁻⁴, 10⁻⁵, 10⁻⁶ — entries, kernels, coefficients, the inactive kernels' displacement, the tracking error.
B. The controller's knobs at ε = 10⁻⁵: k_e = 2 … 40, ω = 0.05 and 1, a = 0.6 and 0.8, a quasi-periodic reference — the real loop against
   the exact averaged loop (where the frozen orbit exists): the count, the kernels, the kernels pushed out of [−1.2, 1.2], the entered
   centres' swing within a period at the end.
C. A feature near the reference's range (x_m = 0.45 + 0.5 sin ωt, the range [−0.05, 0.95]): the first bump moved to 0.11, 0.19, 0.35 below
   it (0.9σ, 1.6σ, 2.9σ), the exact averaged loop at ε = 10⁻⁶ to τ = 5000: whether and when the bump is formed, against the averaged
   drive's size at its nearest kernel before any entry; and Part C's original bump at −0.52.
D. The positions' metric: the exact averaged loop under the Wasserstein–Fisher–Rao metric and with fixed centres, ε = 10⁻⁶ and 10⁻⁸.
E. Measurement noise (ε = 10⁻⁴; 16 tones in 10–30 rad/s of rms s in the measured x, seen by the error, the control and the regressor):
   s = 0.005, 0.01, 0.02 with the free law, and with Proposition D‴.46's leak −λη toward the start (λΛ(ε) = LEAK_E in the time τ — the
   leak's weight on the scaled dual u = η/Λ, a LASSO weight λΛρ on the averaged fit).

Output: results/chaos/closed_transit.md and closed_transit_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D,E`; progress to stdout).
"""
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mirror_free as MF
from src.core import rich_transit as RT
from src.core import closed_loop as CL
import exp_cohesion_leak as E

OUT = os.environ.get('CT_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('CT_WORKERS', '3'))
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)                      # c = 0.025
KE, OM, SIG = 10.0, 0.2, 0.12
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
C0 = np.linspace(-1, 1, 21)
M_ORBIT, TAU_MAX, TH = 512, 600.0, 1e-2                              # orbit samples; horizon in τ; formed: |w| > TH·max|A*|
EPS_A_AVG, EPS_A_REAL = (1e-4, 1e-5, 1e-6, 1e-7), (1e-4, 1e-5, 1e-6)
EPS_B = 1e-5
KE_B = (2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0, 40.0)
REF_B = (('ω = 0.05', dict(om=0.05)), ('ω = 1', dict(om=1.0)), ('a = 0.6', dict(a=0.6)), ('a = 0.8', dict(a=0.8)),
         ('quasi-periodic 0.6 sin ωt + 0.4 sin √2ωt', dict(a=0.6, a2=0.4)))
REF_HALF, D_C, EPS_C, TAU_C = (0.45, 0.5), (0.11, 0.19, 0.35), 1e-6, 5000.0
EPS_D, TAU_D = (1e-6, 1e-8), 3000.0
EPS_E, SIG_E, LEAK_E, TAU_E = 1e-4, (0.005, 0.01, 0.02), 3e-3, 600.0
LEAKS_E2 = (3e-4, 1e-3)                                               # E2: weaker leaks, without noise and at s = 0.01
OUTSIDE = 1.2                                                         # a kernel pushed out: |c| > OUTSIDE


def law():
    return MF.FreeLaw(**LAW)


def phi(x, c):
    return np.exp(-(x - c) ** 2 / (2 * SIG ** 2))


def make_f(cstar=CSTAR, astar=ASTAR):
    cstar, astar = np.asarray(cstar, float), np.asarray(astar, float)
    return lambda x: np.exp(-(np.asarray(x, float)[..., None] - cstar) ** 2 / (2 * SIG ** 2)) @ astar


def th_abs(astar=ASTAR):
    return TH * np.abs(astar).max()


def crossings(W, taus, th):
    """(τ, i) of each kernel's first |w| > th, interpolated in log|w| between output points."""
    out = []
    for i in range(W.shape[1]):
        k = np.flatnonzero(np.abs(W[:, i]) > th)
        if k.size:
            k = k[0]
            if k == 0:
                out.append((float(taus[0]), i))
                continue
            a, b = np.log(max(abs(W[k - 1, i]), 1e-300)), np.log(abs(W[k, i]))
            out.append((float(taus[k - 1] + (taus[k] - taus[k - 1]) * (np.log(th) - a) / (b - a)), i))
    return sorted(out)


def averaged(L, eps, ke=KE, om=OM, x0=0.0, a=1.0, cstar=CSTAR, astar=ASTAR, metric='paradigm', exact=True, tau_max=TAU_MAX, n_t=2401):
    """The averaged loop from α = 0, masses ε, to τ_max: exact (the frozen orbit) or its leading order (the reference's orbit through the
    filter). Returns None for the exact loop when the frozen orbit does not converge (k_e too small against Lip(f − f̂))."""
    d = C0.size
    Lam = RT.transit_scale(L, eps)
    taus = np.linspace(0, tau_max, n_t)
    f = make_f(cstar, astar)
    xm = lambda t: x0 + a * np.sin(om * t)
    if exact:
        O = RT.FrozenOrbit(xm, 2 * np.pi / om, M_ORBIT, SIG, ke, f)
        rhs = RT.frozen_orbit_rhs(L, O, np.full(d, eps), metric)
    else:
        K = RT.OrbitBumps(xm, 2 * np.pi / om, M_ORBIT, SIG, ke)
        rhs = RT.critical_rhs(L, K, f(K.xs), 'paradigm' if metric == 'fixed' else metric, b0=np.full(d, eps))
        if metric == 'fixed':
            base = rhs
            rhs = lambda t, z: np.r_[base(t, z)[:d], np.zeros(d)]
    try:
        sol = solve_ivp(rhs, (0, Lam * tau_max), np.r_[np.zeros(d), C0], method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
    except RuntimeError:
        return None
    al, cc = sol.y[:d].T, sol.y[d:].T
    W = (eps + 2 * L.c * np.sinh(al / L.A) ** 2) * L.A * np.tanh(al / L.A) / L.m
    th = th_abs(astar)
    act = np.abs(W[-1]) > th
    ent = crossings(W, taus, th)
    cl = np.sort(cc[-1][act])
    if exact:
        e = O.orbit(cc[-1], W[-1])
        rms = float(np.sqrt(np.mean(e ** 2)))
    else:
        rms = float(np.sqrt(np.mean((f(K.xs) - K.phi(cc[-1]) @ W[-1]) ** 2)))
    return dict(eps=eps, Lam=float(Lam), entries=ent, n_act=int(act.sum()), kernels=np.flatnonzero(act).tolist(),
                centres=cc[-1][act].tolist(), w=W[-1][act].tolist(), clusters=int(np.sum(np.diff(cl) > 1e-3) + 1) if cl.size else 0,
                rms=rms, disp_in=float(np.abs(cc[-1] - C0)[~act].max()) if (~act).any() else 0.0,
                t_last=float(Lam * max(t for t, i in ent)) if ent else np.nan, done=sol.status == 0)


def real(L, eps, T, ke=KE, om=OM, x0=0.0, a=1.0, a2=0.0, cstar=CSTAR, astar=ASTAR, noise=None, lam=0.0, n_end=6, n_out=401):
    """The real loop (C) from x(0) = x_m(0), α = 0, masses ε, to the time T; a dense record of the last n_end periods for the swing."""
    per = 2 * np.pi / om
    ts = np.unique(np.r_[np.linspace(0, T, n_out), np.linspace(max(T - n_end * per, 0), T, 60 * n_end)])
    ref = (x0, a, om, a2, np.sqrt(2) * om)
    th = th_abs(astar)
    t0 = time.time()
    r = CL.simulate(L, eps, C0, SIG, (cstar, astar, SIG), ke, ref, ts, noise=noise, lam=lam, th=th)
    Lam = RT.transit_scale(L, eps)
    act = np.abs(r['w'][-1]) > th
    ent = sorted((float(t / Lam), i) for i, t in enumerate(r['entries']) if np.isfinite(t))
    last = ts >= T - per
    swing = float(np.ptp(r['c'][last][:, act], axis=0).max()) if act.any() else 0.0
    return dict(eps=eps, T=float(T), entries=ent, n_act=int(act.sum()), kernels=np.flatnonzero(act).tolist(),
                centres=r['c'][-1][act].tolist(), w=r['w'][-1][act].tolist(), swing=swing,
                outside=int(np.sum(np.abs(r['c'][-1][act]) > OUTSIDE)), wmax=float(np.abs(r['w'][-1]).max()),
                disp_in=float(np.abs(r['c'][-1] - C0)[~act].max()) if (~act).any() else 0.0,
                e_first=r['e_first'], e_last=r['e_last'], secs=time.time() - t0, steps=r['steps'])


def _pool():
    from multiprocessing import Pool
    return Pool(N_WORKERS)


def a_avg(args):
    eps, exact = args
    warnings.filterwarnings('ignore')
    return averaged(law(), eps, exact=exact)


def a_real(eps):
    warnings.filterwarnings('ignore')
    L = law()
    av = averaged(L, eps)
    return real(L, eps, 1.3 * av['t_last'])


def part_A():
    with _pool() as pool:
        ex = pool.map(a_avg, [(e, True) for e in EPS_A_AVG])
        lo = pool.map(a_avg, [(e, False) for e in EPS_A_AVG])
        rl = pool.map(a_real, list(EPS_A_REAL))
    return dict(exact=ex, lead=lo, real=rl)


def b_job(args):
    label, kw = args
    warnings.filterwarnings('ignore')
    L = law()
    ke, om, a, a2 = kw.get('ke', KE), kw.get('om', OM), kw.get('a', 1.0), kw.get('a2', 0.0)
    av = averaged(L, EPS_B, ke=ke, om=om, a=a) if a2 == 0.0 else None
    tl = av['t_last'] if (av is not None and np.isfinite(av['t_last'])) else RT.transit_scale(L, EPS_B) * 300.0 * ke / KE
    rl = real(L, EPS_B, 1.3 * tl, ke=ke, om=om, a=a, a2=a2)
    print(f"  B {label}: averaged {av['n_act'] if av else '—'}, real {rl['n_act']} ({rl['secs']:.0f}s)", flush=True)
    return dict(label=label, kw=kw, avg=av, real=rl)


def part_B():
    jobs = [(f'k_e = {k:g}', dict(ke=k)) for k in KE_B] + list(REF_B)
    with _pool() as pool:
        return pool.map(b_job, jobs)


def c_job(dist):
    warnings.filterwarnings('ignore')
    L = law()
    lo = REF_HALF[0] - REF_HALF[1]
    cstar = CSTAR.copy()
    cstar[0] = lo - dist if dist is not None else CSTAR[0]
    av = averaged(L, EPS_C, x0=REF_HALF[0], a=REF_HALF[1], cstar=cstar, tau_max=TAU_C, n_t=20001)
    O = RT.FrozenOrbit(lambda t: REF_HALF[0] + REF_HALF[1] * np.sin(OM * t), 2 * np.pi / OM, M_ORBIT, SIG, KE, make_f(cstar))
    g, h, e = O.drive(C0, np.zeros(C0.size))
    j = int(np.argmin(np.abs(C0 - cstar[0])))
    first = [t for t, i in av['entries'] if abs(C0[i] - cstar[0]) <= 0.15]
    return dict(dist=dist, bump=float(cstar[0]), avg=av, drive_near=float(abs(g[j])), drive_max=float(np.abs(g).max()), kernel=j,
                t_first=float(min(first)) if first else np.nan,
                fitted=bool(any(abs(c - cstar[0]) < 0.01 for c in av['centres'])))


def part_C():
    with _pool() as pool:
        return pool.map(c_job, list(D_C) + [None])


def d_job(args):
    metric, eps = args
    warnings.filterwarnings('ignore')
    return averaged(law(), eps, metric=metric, tau_max=TAU_D, n_t=12001)


def part_D():
    with _pool() as pool:
        rows = pool.map(d_job, [(m, e) for m in ('wfr', 'fixed') for e in EPS_D])
    return [dict(metric=m, **r) for (m, e), r in zip([(m, e) for m in ('wfr', 'fixed') for e in EPS_D], rows)]


def e_job(args):
    s, lam = args
    warnings.filterwarnings('ignore')
    L = law()
    Lam = RT.transit_scale(L, EPS_E)
    r = real(L, EPS_E, Lam * TAU_E, noise=CL.tones(s, 10.0, 30.0, seed=7) if s > 0 else None, lam=lam / Lam)
    print(f"  E s {s} λ {lam}: {r['n_act']} kernels, {r['outside']} outside, max|w| {r['wmax']:.3g} ({r['secs']:.0f}s)", flush=True)
    return dict(s=s, lam=lam, **r)


def part_E():
    jobs = [(0.0, 0.0), (0.0, LEAK_E)] + [(s, lam) for lam in (0.0, LEAK_E) for s in SIG_E]
    with _pool() as pool:
        return pool.map(e_job, jobs)


def part_E2():
    with _pool() as pool:
        return pool.map(e_job, [(s, lam) for lam in LEAKS_E2 for s in (0.0, 0.01)])


def report(A=None, B=None, C=None, D=None, E_=None, E2=None):
    fmt = E.fmt
    L = law()
    ev = lambda ent: ', '.join(f"{i}: {t:.4g}" for t, i in ent)
    rnd = lambda v: [round(x, 4) for x in v]
    L_ = ["# Proposition D‴.49: the rich transit in the closed loop (`exp_closed_transit.py`)", "",
          f"ẋ = f(x) + u, u = ẋ_m − k_ee − f̂(x), x_m = x₀ + a sin(ωt) (base point k_e = {KE:g}, ω = {OM:g}, x₀ = 0, a = 1); f three bumps of "
          f"width {SIG} at {CSTAR.tolist()} with amplitudes {ASTAR.tolist()} (max|f′| = {np.abs(ASTAR).max() / (SIG * np.exp(0.5)):.3g}); "
          f"21 kernels of width {SIG} started on a grid of [−1, 1] at α = 0 with the masses ε; the free law (m = {LAW['m']:g}, G = "
          f"{LAW['G']:g}, γ = {LAW['gamma']:g}, a_max = {LAW['a_max']:g}) with moving centres; τ = t/Λ(ε), Λ = (πa_max/(4G))(2/(cε))^{{1/2}}. "
          f"A kernel has formed when |w| > {TH:g}·max|A*|; entries listed as kernel: τ (kernel i started at −1 + 0.1i), interpolated. The "
          "real loop in C (`closed_loop`); the exact averaged loop along the frozen periodic orbit (`FrozenOrbit`); its leading order in "
          "1/k_e along the reference's orbit (`OrbitBumps`).", ""]
    if A is not None:
        L_ += ["## A. The base point: the exact averaged loop, its leading order and the real loop", "",
               "| ε | loop | entries (kernel: τ) | formed | centres | coefficients | inactive: max \\|Δc\\| | max \\|e\\|, last period |",
               "|---|---|---|---|---|---|---|---|"]
        reals = {r['eps']: r for r in A['real']}
        for ex, lo in zip(A['exact'], A['lead']):
            L_.append(f"| {ex['eps']:g} | averaged, exact | {ev(ex['entries'])} | {ex['n_act']} | {rnd(ex['centres'])} | {rnd(ex['w'])} | "
                      f"{fmt(ex['disp_in'])} | — |")
            L_.append(f"| {lo['eps']:g} | averaged, leading order | {ev(lo['entries'])} | {lo['n_act']} | {rnd(lo['centres'])} | "
                      f"{rnd(lo['w'])} | {fmt(lo['disp_in'])} | — |")
            if ex['eps'] in reals:
                r = reals[ex['eps']]
                L_.append(f"| {r['eps']:g} | real (T = {r['T']:.3g}, {r['secs']:.0f} s) | {ev(r['entries'])} | {r['n_act']} | "
                          f"{rnd(r['centres'])} | {rnd(r['w'])} | {fmt(r['disp_in'])} | {fmt(r['e_last'])} |")
        L_.append("")
    if B is not None:
        L_ += [f"## B. The controller's knobs (ε = {EPS_B:g})", "",
               "The real loop against the exact averaged loop (— where the frozen orbit does not converge or the reference is not periodic): "
               f"the kernels formed, those pushed out of [−{OUTSIDE:g}, {OUTSIDE:g}], the largest |w|, the entered centres' largest swing "
               "within the last period, the last period's max|e|.", "",
               "| knob | averaged: kernels | real: entries (kernel: τ) | real: kernels | pushed out | max \\|w\\| | swing | max \\|e\\|, last period |",
               "|---|---|---|---|---|---|---|---|"]
        for b in B:
            av, r = b['avg'], b['real']
            L_.append(f"| {b['label']} | {av['kernels'] if av else '—'} | {ev(r['entries'])} | {r['kernels']} | {r['outside']} | "
                      f"{r['wmax']:.3g} | {fmt(r['swing'])} | {fmt(r['e_last'])} |")
        L_.append("")
    if C is not None:
        lo = REF_HALF[0] - REF_HALF[1]
        L_ += [f"## C. A feature near the range (x_m = {REF_HALF[0]} + {REF_HALF[1]} sin ωt, the range [{lo:g}, {REF_HALF[0] + REF_HALF[1]:g}]; "
               f"the exact averaged loop, ε = {EPS_C:g}, to τ = {TAU_C:g})", "",
               "| the first bump | below the range | the drive at its nearest kernel before any entry (of the largest) | entry near it (τ) | "
               "fitted | kernels | centres | coefficients |", "|---|---|---|---|---|---|---|---|"]
        for c in C:
            av = c['avg']
            below = f"{c['dist']:g} ({c['dist'] / SIG:.2g}σ)" if c['dist'] is not None else f"{lo - c['bump']:g} ({(lo - c['bump']) / SIG:.2g}σ)"
            L_.append(f"| {c['bump']:g} | {below} | {c['drive_near']:.3g} ({c['drive_near'] / c['drive_max']:.3g}) | "
                      f"{fmt(c['t_first']) if np.isfinite(c['t_first']) else 'none'} | {'yes' if c['fitted'] else 'no'} | {av['kernels']} | "
                      f"{rnd(av['centres'])} | {rnd(av['w'])} |")
        L_.append("")
    if D is not None:
        L_ += ["## D. The positions' metric (the exact averaged loop)", "",
               "| metric | ε | formed | clusters | centres | rms of the orbit's error |", "|---|---|---|---|---|---|"]
        for r in D:
            L_.append(f"| {r['metric']} | {r['eps']:g} | {r['n_act']} | {r['clusters']} | {sorted(set(round(v, 3) for v in r['centres']))} | "
                      f"{fmt(r['rms'])} |")
        L_.append("")
    if E_ is not None:
        L_ += [f"## E. Measurement noise (ε = {EPS_E:g}, τ ≤ {TAU_E:g}; 16 tones in 10–30 rad/s of rms s in the measured state)", "",
               "| s | leak λΛ | kernels | pushed out | max \\|w\\| | max \\|e\\|, last period | entries (kernel: τ) |", "|---|---|---|---|---|---|---|"]
        for r in list(E_) + list(E2 or []):
            L_.append(f"| {r['s']:g} | {r['lam']:g} | {r['kernels']} | {r['outside']} | {r['wmax']:.3g} | {fmt(r['e_last'])} | "
                      f"{ev(r['entries'])} |")
        if E2:
            L_ += ["", f"The rows with λΛ = {', '.join(f'{v:g}' for v in LEAKS_E2)} (E2): the leak's strength between none and {LEAK_E:g}, "
                   "without noise and at s = 0.01 — the mirror leak's bias (the amplitudes shrink toward the leak's centre) against its "
                   "hold on the kernels' positions."]
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_closed_transit.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'closed_transit_raw.npz')
    md = os.path.join(OUT, 'closed_transit.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'E', 'E2')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = {k: v for k, v in old.items() if k in keys and k not in parts}
    t0 = time.time()
    for k, fn in (('A', part_A), ('C', part_C), ('D', part_D), ('E', part_E), ('E2', part_E2), ('B', part_B)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
