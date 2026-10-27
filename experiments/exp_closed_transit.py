"""
Proposition D‴.49: the rich transit in the closed loop — from small masses the adaptive controller forms one kernel per feature
====================================================================================================================================

The scalar plant ẋ = f(x) + u, u = ẋ_m − k_ee − f̂(x), e = x − x_m, x_m = x₀ + a sin(ωt); f̂ = Σ_iw_iφ(·; c_i), Gaussian kernels of width 0.12
started at α = 0 with the masses ε on a grid of 21 centres in [−1, 1]; the paradigm's free law with moving centres driven by the error,
dα_i = Gβ_ia′_ieφ_i(x), dβ_i = γβ_ia_ieφ_i(x), ċ_i = Gβ_ia_ie∂_cφ_i(x) (Proposition D‴.48(iii); `rich_transit`). f = three bumps of the kernels'
width at −0.52, 0.07, 0.61 with amplitudes 0.9, −0.6, 0.75.

A. The averaged loop (the orbit x_m sampled over a period, the error filter L = 1/(s + k_e); `rich_transit.OrbitBumps`) for ε = 10⁻⁴ …
   10⁻⁷: the entries (τ = t/Λ(ε)), the formed kernels, their centres and coefficients, the last entry's time t = Λτ.
B. The real loop at ε = 10⁻⁴ and 10⁻⁵ against A: the entries, the formed kernels, the inactive kernels' displacement, the tracking error
   over the last period.
C. A reference that does not visit the bump at −0.52 (x_m = 0.45 + 0.5 sin ωt): the averaged loop at 10⁻⁶ and the real loop at 10⁻⁴.
D. The positions' metric in the averaged loop at ε = 10⁻⁶: the Wasserstein–Fisher–Rao metric (condensation) and fixed centres.

Output: results/chaos/closed_transit.md and closed_transit_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D`; progress to stdout).
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
import exp_cohesion_leak as E

OUT = os.environ.get('CT_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)                      # c = 0.025
KE, OM, SIG = 10.0, 0.2, 0.12
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
C0 = np.linspace(-1, 1, 21)
REF_FULL, REF_HALF = (0.0, 1.0), (0.45, 0.5)                          # x_m = x0 + a sin(ωt)
EPS_A, EPS_B, EPS_C, EPS_D = (1e-4, 1e-5, 1e-6, 1e-7), (1e-4, 1e-5), (1e-6, 1e-4), 1e-6
M_ORBIT, TAU_MAX = 512, 600.0
TH = 1e-2                                                             # formed: |w| > TH·max|A*|


def law():
    return MF.FreeLaw(**LAW)


def phi(x, c):
    return np.exp(-(x - c) ** 2 / (2 * SIG ** 2))


def fstar(x):
    return phi(np.asarray(x, float)[..., None], CSTAR) @ ASTAR


def orbit(ref):
    x0, a = ref
    return RT.OrbitBumps(lambda t: x0 + a * np.sin(OM * t), 2 * np.pi / OM, M_ORBIT, SIG, KE)


def entries(W, taus, th):
    out = []
    for i in range(W.shape[1]):
        k = np.flatnonzero(np.abs(W[:, i]) > th)
        if k.size:
            out.append((float(taus[k[0]]), i))
    return sorted(out)


def averaged(L, ref, eps, metric='paradigm', fixed=False, tau_max=TAU_MAX, n_t=1201):
    """The averaged loop from α = 0, masses ε, to τ_max (τ = t/Λ(ε))."""
    K = orbit(ref)
    fs = fstar(K.xs)
    d = C0.size
    Lam = RT.transit_scale(L, eps)
    taus = np.linspace(0, tau_max, n_t)
    z0 = np.r_[np.zeros(d), np.log(np.full(d, eps)), C0]
    f = RT.moving_rhs(L, K, fs, metric)
    rhs = (lambda t, z: np.r_[f(t, z)[:2 * d], np.zeros(d)]) if fixed else f
    sol = solve_ivp(rhs, (0, Lam * tau_max), z0, method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
    W = np.exp(sol.y[d:2 * d].T) * L.amp(sol.y[:d].T)[0] / L.m
    cc = sol.y[2 * d:].T
    th = TH * np.abs(ASTAR).max()
    act = np.abs(W[-1]) > th
    r = fs - K.phi(cc[-1]) @ W[-1]
    ent = entries(W, taus, th)
    cl = np.sort(cc[-1][act])
    return dict(eps=eps, Lam=float(Lam), entries=ent, n_act=int(act.sum()), centres=cc[-1][act].tolist(), w=W[-1][act].tolist(),
                clusters=int(np.sum(np.diff(cl) > 1e-3) + 1) if cl.size else 0, rms=float(np.sqrt(np.mean(r ** 2))),
                disp_in=float(np.abs(cc[-1] - C0)[~act].max()) if (~act).any() else 0.0,
                t_last=float(Lam * max(t for t, i in ent)) if ent else np.nan)


def real_loop(L, ref, eps, T, n_t=801):
    """The real loop (x, α, log β, c) from x(0) = x_m(0), α = 0, masses ε, to the time T."""
    x0, a = ref
    d = C0.size

    def f(t, z):
        x = z[0]
        al, lb, c = z[1:1 + d], z[1 + d:1 + 2 * d], z[1 + 2 * d:]
        b = np.exp(lb); amp, da = L.amp(al)
        ph = phi(x, c)
        e = x - (x0 + a * np.sin(OM * t))
        g = e * ph
        return np.r_[fstar(x) - ph @ (b * amp / L.m) + a * OM * np.cos(OM * t) - KE * e, L.G * b * da * g, L.gamma * amp * g,
                     L.G * b * amp * e * ph * (x - c) / SIG ** 2]
    z0 = np.r_[x0, np.zeros(d), np.log(np.full(d, eps)), C0]
    per = 2 * np.pi / OM
    ts = np.unique(np.r_[np.linspace(0, T, n_t), np.linspace(0, per, 200), np.linspace(T - per, T, 200)])   # dense first and last periods
    sol = solve_ivp(f, (0, T), z0, method='LSODA', rtol=1e-8, atol=1e-11, t_eval=ts, max_step=0.2 / OM)
    W = np.exp(sol.y[1 + d:1 + 2 * d].T) * L.amp(sol.y[1:1 + d].T)[0] / L.m
    cc = sol.y[1 + 2 * d:].T
    th = TH * np.abs(ASTAR).max()
    act = np.abs(W[-1]) > th
    Lam = RT.transit_scale(L, eps)
    ent = entries(W, sol.t / Lam, th)
    e_end = sol.y[0] - (x0 + a * np.sin(OM * sol.t))
    k_last = sol.t >= sol.t[-1] - per
    return dict(eps=eps, T=float(T), entries=ent, n_act=int(act.sum()), centres=cc[-1][act].tolist(), w=W[-1][act].tolist(),
                disp_in=float(np.abs(cc[-1] - C0)[~act].max()) if (~act).any() else 0.0,
                e_last=float(np.abs(e_end[k_last]).max()), e_first=float(np.abs(e_end[sol.t <= per]).max()))


def part_A():
    warnings.filterwarnings('ignore')
    L = law()
    rows = []
    for eps in EPS_A:
        t0 = time.time()
        rows.append(averaged(L, REF_FULL, eps))
        print(f"  A eps {eps:g}: {rows[-1]['n_act']} formed ({time.time() - t0:.0f}s)", flush=True)
    return rows


def part_B():
    warnings.filterwarnings('ignore')
    L = law()
    rows = []
    for eps in EPS_B:
        av = averaged(L, REF_FULL, eps)
        T = 1.3 * av['t_last']
        t0 = time.time()
        rl = real_loop(L, REF_FULL, eps, T)
        rl['secs'] = time.time() - t0
        rows.append(dict(avg=av, real=rl))
        print(f"  B eps {eps:g}: real loop {rl['n_act']} formed ({rl['secs']:.0f}s)", flush=True)
    return rows


def part_C():
    warnings.filterwarnings('ignore')
    L = law()
    av = averaged(L, REF_HALF, EPS_C[0])
    av4 = averaged(L, REF_HALF, EPS_C[1])
    t0 = time.time()
    rl = real_loop(L, REF_HALF, EPS_C[1], 1.3 * av4['t_last'])
    rl['secs'] = time.time() - t0
    print(f"  C: averaged {av['n_act']}, real loop {rl['n_act']} ({rl['secs']:.0f}s)", flush=True)
    return dict(avg=av, avg4=av4, real=rl)


def part_D():
    warnings.filterwarnings('ignore')
    L = law()
    return dict(wfr=averaged(L, REF_FULL, EPS_D, metric='wfr'), fixed=averaged(L, REF_FULL, EPS_D, fixed=True))


def report(A=None, B=None, C=None, D=None):
    fmt = E.fmt
    L = law()
    ev = lambda ent: ', '.join(f"{i}: {t:.4g}" for t, i in ent)
    L_ = ["# Proposition D‴.49: the rich transit in the closed loop (`exp_closed_transit.py`)", "",
          f"ẋ = f(x) + u, u = ẋ_m − k_ee − f̂(x), k_e = {KE:g}, x_m = x₀ + a sin(ωt), ω = {OM:g}; f three bumps of width {SIG} at "
          f"{CSTAR.tolist()} with amplitudes {ASTAR.tolist()}; 21 kernels of width {SIG} started on a grid of [−1, 1] at α = 0 with the masses ε; "
          f"the free law (m = {LAW['m']:g}, G = {LAW['G']:g}, γ = {LAW['gamma']:g}, a_max = {LAW['a_max']:g}) with moving centres; τ = t/Λ(ε), "
          f"Λ = (πa_max/(4G))(2/(cε))^{{1/2}}. A kernel has formed when |w| > {TH:g}·max|A*|; entries listed as kernel: τ (kernel i started at "
          f"−1 + 0.1i).", ""]
    if A is not None:
        L_ += ["## A. The averaged loop (the reference's orbit sampled over a period, the error filter 1/(s + k_e))", "",
               "| ε | entries (kernel: τ) | the last entry's time t | formed | centres | coefficients | rms residual | inactive: max |Δc| |",
               "|---|---|---|---|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | {ev(r['entries'])} | {r['t_last']:.3g} | {r['n_act']} | {[round(v, 4) for v in r['centres']]} | "
               f"{[round(v, 4) for v in r['w']]} | {fmt(r['rms'])} | {fmt(r['disp_in'])} |" for r in A]
        L_.append("")
    if B is not None:
        L_ += ["## B. The real loop against the averaged loop (x_m = sin ωt)", "",
               "| ε | | entries (kernel: τ) | formed | centres | coefficients | inactive: max |Δc| | max |e|: first, last period |",
               "|---|---|---|---|---|---|---|---|"]
        for r in B:
            a, q = r['avg'], r['real']
            L_.append(f"| {a['eps']:g} | averaged | {ev(a['entries'])} | {a['n_act']} | {[round(v, 4) for v in a['centres']]} | "
                      f"{[round(v, 4) for v in a['w']]} | {fmt(a['disp_in'])} | — |")
            L_.append(f"| {q['eps']:g} | real (T = {q['T']:.3g}, {q['secs']:.0f} s) | {ev(q['entries'])} | {q['n_act']} | "
                      f"{[round(v, 4) for v in q['centres']]} | {[round(v, 4) for v in q['w']]} | {fmt(q['disp_in'])} | "
                      f"{fmt(q['e_first'])}, {fmt(q['e_last'])} |")
        L_.append("")
    if C is not None:
        L_ += [f"## C. A reference that does not visit the bump at −0.52 (x_m = {REF_HALF[0]} + {REF_HALF[1]} sin ωt, the range "
               f"[{REF_HALF[0] - REF_HALF[1]:g}, {REF_HALF[0] + REF_HALF[1]:g}])", "",
               "| | ε | entries (kernel: τ) | formed | centres | coefficients |", "|---|---|---|---|---|---|"]
        for lab, r in (('averaged', C['avg']), ('averaged', C['avg4']), ('real', C['real'])):
            L_.append(f"| {lab} | {r['eps']:g} | {ev(r['entries'])} | {r['n_act']} | {[round(v, 4) for v in r['centres']]} | "
                      f"{[round(v, 4) for v in r['w']]} |")
        L_.append("")
    if D is not None:
        L_ += [f"## D. The positions' metric (the averaged loop, ε = {EPS_D:g})", ""]
        for lab, r in (('Wasserstein–Fisher–Rao (ċ = Ga e∂_cφ)', D['wfr']), ('fixed centres', D['fixed'])):
            L_.append(f"- {lab}: {r['n_act']} kernels formed in {r['clusters']} clusters at {sorted(set(round(v, 3) for v in r['centres']))}, "
                      f"rms residual {fmt(r['rms'])}.")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_closed_transit.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'closed_transit_raw.npz')
    md = os.path.join(OUT, 'closed_transit.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCD']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCD')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('D', part_D), ('C', part_C), ('B', part_B)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCD']))


if __name__ == '__main__':
    main()
