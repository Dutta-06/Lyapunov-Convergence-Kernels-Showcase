"""
Proposition D‴.39: dissipativity of the dithered gauge — the escape to saturated approximate fits, and a dead-zone confinement
==============================================================================================================================

Proposition D‴.36's instance (Proposition D‴.32's three bumps, one particle per class, m = 3, the cohesion leak), the positions dithered at
τ. Under the cohesion leak J is bounded on the unconstrained coordinates, so the positions' Gibbs law e^{−J/τ} is not normalisable; the
dither leaves the exact fit's well for approximate fits with saturated amplitudes (the plateau), on which the saturated coordinates diffuse
freely. The dead zone V = (σ_B/2)Σ_i‖ϑ_i − Π_B(ϑ_i)‖²_{G⁻¹} (`ReplicatorCohesionLoop(box=L, box_sigma=σ_B)`, B = Π_j[−L_j, L_j] in the
unconstrained coordinates (α, γ, κ)) makes it normalisable and leaves everything inside B unchanged.

A. The plateau: at the masses 1/3 and W♯ (Proposition D‴.38), for each class k the lowest configuration with α_k fixed at ±20 (the sign of
   A_k; L-BFGS over the other eight coordinates from N_A perturbed exact fits), its level ΔJ/τ above the exact fit (τ = 10⁻⁴) and its fit
   error; with α_k at ±ESC (the escape set's boundary, Part C's ΔJ); and with the box, α_k at its edge ±L_α.
B. The minimum-energy path from the exact fit to the plateau (the string method in the G⁻¹ metric, N_B images, the endpoint the lowest
   plateau of A with its saturated amplitude set to 10 and the rest relaxed): the profile ΔJ/τ along it.
C. The escapes at fixed masses (γ = 0, W = 1/3, τ = 1.5, 2, 3·10⁻⁴, 24 seeds, at most S_C slow units): the first time max_k|α_k| > ESC; the
   escape rate (escapes over the total time of the runs, censored runs included) against e^{−ΔJ/τ}, ΔJ the lowest level of J on the escape
   set's boundary {max_k|α_k| = ESC} (A; Freidlin–Wentzell: τ log of the mean exit time tends to it).
D. The dead zone at fixed masses (γ = 0, W = 1/3, τ = 2·10⁻⁴, 16 seeds, S_D slow units): without and with the box L = BOX (σ_B = 1): the
   fraction of the records with max|α| > ESC and the fit error by quarters, the largest |α|.
E. The full loop with the box (γ = 10⁻², τ = 10⁻⁴, the 12 seeds of exp_gauge_drift.py's Part H, 10⁶ slow units): κ_F = 0, and the masses
   form at κ_F = 2: the masses at 10⁵, 3·10⁵, 10⁶ against the vertices of the box's representable masses {W_k ≥ |A_k|/(a_max tanh(L_α/a_max))}
   (κ_F = 0) and W♯ (κ_F = 2); κ_F = 0 without the box is exp_gauge_drift.py's Part H.

Output: results/chaos/gauge_plateau.md and gauge_plateau_raw.npz (GP_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mass_transfer as MT
import exp_cohesion_leak as E
import exp_gauge_dither as GD
import exp_gauge_drift as X
import exp_gauge_fixman as Q

OUT = os.environ.get('GP_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GP_WORKERS', 4))
TAU, A_MAX, DT, REC = GD.TAU, 4.0, GD.DT, GD.REC
BOX, SIGMA_B = (5.0, 3.0, 4.0), 1.0                                         # L_α, L_γ, L_κ (|a| ≤ 3.39, |c| ≤ 2.28, s ∈ [0.33, 1.97])
AL_SAT, N_A = 20.0, 200
N_B, IT_B, DT_B, AL_B = 40, 30000, 0.5, 10.0
ESC = 6.0                                                                    # |a|/a_max > 0.905; in the well |α_k| < 4.5
S2_PLATEAU = 0.6                                                             # the second width on the plateaus' family (exact fit 0.4)
TAUS_C, SEEDS_C, S_C, REC_C = (1.5e-4, 2e-4, 3e-4), tuple(range(500, 524)), 4e5, 40
TAU_D, SEEDS_D, S_D = 2e-4, tuple(range(600, 616)), 4e5
GAMMA_E, SEEDS_E, S_E, T_E = 1e-2, X.SEEDS_H, 1e6, X.T_H
W_MIN = np.abs(E.A_BUMP) / A_MAX
W_BOX = np.abs(E.A_BUMP) / (A_MAX * np.tanh(BOX[0] / A_MAX))


def masses(which):
    return np.ones(3) / 3 if which == 'uniform' else Q.w_sharp()


def plateau_job(args):
    """The lowest J with α_k fixed at al (with the box: J + V), over N_A L-BFGS runs from the exact fit perturbed (seeded)."""
    from scipy.optimize import minimize
    which, k, at = args
    boxed = at == 'box'
    net, base = E.setup()
    W = masses(which); beta = 3 * W
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=(BOX if boxed else None), box_sigma=SIGMA_B)
    Th0 = E.exact_rep(net, W); J0 = loop.J(Th0, beta)
    al = np.sign(E.A_BUMP[k]) * dict(sat=AL_SAT, esc=ESC, box=BOX[0])[at]
    f = lambda x: loop.J(np.insert(x, 3 * k, al).reshape(3, 3), beta)
    rng = np.random.default_rng(1000 + 10 * k + ('sat', 'esc', 'box').index(at))
    best = None
    for _ in range(N_A):
        x0 = np.delete(Th0.ravel(), 3 * k) + rng.normal(0, 1.5, 8) * rng.uniform(0.2, 1.5)
        r = minimize(f, x0, method='L-BFGS-B', options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12))
        if best is None or r.fun < best.fun:
            best = r
    T = np.insert(best.x, 3 * k, al).reshape(3, 3)
    a, c, s = net.physical(T)
    return dict(which=which, k=k, at=at, dJ=(best.fun - J0) / TAU, fit=loop.fit_error(T, beta), a=a[:, 0], c=c[:, 0], s=s, Th=T)


def string_job(args):
    """The string method (the simplified string method of E–Ren–Vanden-Eijnden) in the G⁻¹ metric between the exact fit and the plateau's
    point `end`, the end relaxed with its saturated amplitude fixed at AL_B: the profile ΔJ/τ and the images."""
    from scipy.optimize import minimize
    from scipy.interpolate import interp1d
    which, k, end = args
    net, base = E.setup()
    W = masses(which); beta = 3 * W; G = E.G
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    A = E.exact_rep(net, W); J0 = loop.J(A, beta)
    f = lambda x: loop.J(np.insert(x, 3 * k, np.sign(E.A_BUMP[k]) * AL_B).reshape(3, 3), beta)
    r = minimize(f, np.delete(np.asarray(end, float).ravel(), 3 * k), method='L-BFGS-B', options=dict(maxiter=20000, ftol=1e-15, gtol=1e-12))
    B = np.insert(r.x, 3 * k, np.sign(E.A_BUMP[k]) * AL_B).reshape(3, 3)
    S = np.array([A + (B - A) * t for t in np.linspace(0, 1, N_B)])
    metric = np.sqrt(1 / G)[None, :]

    def repar(S):
        d = np.r_[0, np.cumsum(np.sqrt((((S[1:] - S[:-1]) * metric) ** 2).sum(axis=(1, 2))))]
        return interp1d(d / d[-1], S.reshape(N_B, -1), axis=0)(np.linspace(0, 1, N_B)).reshape(S.shape)
    for _ in range(IT_B):
        for i in range(1, N_B - 1):
            S[i] = S[i] + DT_B * loop.drift(S[i], beta)[0]
        S = repar(S)
    prof = np.array([(loop.J(x, beta) - J0) / TAU for x in S])
    fits = np.array([loop.fit_error(x, beta) for x in S])
    return dict(which=which, k=k, prof=prof, fits=fits, images=S)


class Escaped(Exception):
    pass


def escape_job(args):
    """The first record (every REC_C steps) with max|α_k| > ESC at fixed masses 1/3 (γ = 0), or None by S_C."""
    tau, seed = args
    t0 = time.time()
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    out = dict(tau=tau, seed=seed, t=None)

    def cb(s_, T, b):
        if np.abs(T[:, 0]).max() > ESC:
            out.update(t=s_, alpha=T[:, 0].copy(), fit=loop.fit_error(T, b)); raise Escaped
    try:
        loop.run_dithered(E.exact_rep(net, np.ones(3) / 3), np.ones(3), S_C, DT, tau, np.random.default_rng(seed), rec_every=REC_C, callback=cb)
    except Escaped:
        pass
    print(f"[{time.time() - t0:.0f}s] C tau {tau:g} seed {seed}: {out['t']}", flush=True)
    return out


def series(args):
    """Records (s, α_1..3, W_1..3, the fit error, c_1..3, s_1..3) every REC steps: (γ, τ, κ_F, form, boxed, seed, S)."""
    gamma, tau, kF, form, boxed, seed, S = args
    t0 = time.time()
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3)
    A0 = base.features(Th)[0] / 3 if form == 'masses' else None
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, fixman=kF, fixman_tau=tau, fixman_A=A0,
                                     box=(BOX if boxed else None), box_sigma=SIGMA_B)
    rec = []

    def cb(s_, T, b):
        a, c, s = net.physical(T)
        rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b), c[:, 0], s])
    loop.run_dithered(Th, np.ones(3), S, DT, tau, np.random.default_rng(seed), rec_every=REC, callback=cb)
    print(f"[{time.time() - t0:.0f}s] gamma {gamma:g} tau {tau:g} {form} kF {kF:g} box {boxed} seed {seed}", flush=True)
    return dict(gamma=gamma, tau=tau, kF=kF, form=form, boxed=boxed, seed=seed, rec=np.array(rec))


def part_C(runs, dJ):
    out = []
    for tau in TAUS_C:
        rr = [r for r in runs if r['tau'] == tau]
        t = np.array([r['t'] if r['t'] is not None else S_C for r in rr]); n = sum(r['t'] is not None for r in rr)
        rate = n / t.sum()
        fits = [r['fit'] for r in rr if r['t'] is not None]
        out.append(dict(tau=tau, n=n, N=len(rr), rate=rate, rate_se=rate / np.sqrt(max(n, 1)), med=np.median(t),
                        scaled=rate * np.exp(dJ * TAU / tau), fit=(np.mean(fits) if fits else np.nan)))
    ok = [r for r in out if r['n'] > 0]
    slope = np.polyfit([1 / r['tau'] for r in ok], np.log([r['rate'] for r in ok]), 1)[0] if len(ok) > 1 else np.nan
    return dict(rows=out, slope=slope, dJ=dJ * TAU)


def part_D(runs):
    out = []
    for boxed in (False, True):
        R = np.array([r['rec'] for r in runs if r['boxed'] == boxed]); q = np.array_split(np.arange(R.shape[1]), 4)
        al = np.abs(R[:, :, 1:4]).max(axis=2); pl = R[:, :, 12] > S2_PLATEAU
        out.append(dict(boxed=boxed, frac=[float(np.mean(al[:, ix] > ESC)) for ix in q], plat=[float(np.mean(pl[:, ix])) for ix in q],
                        fit=[float(R[:, ix, 7].mean()) for ix in q], almax=[float(al[:, ix].max()) for ix in q],
                        seeds=int(np.sum(pl.any(axis=1))), N=len(R)))
    return out


def summarise_E(R):
    s = R[0, :, 0]; W = R[:, :, 4:7]; al = R[:, :, 1:4]
    out = dict(times=[], half=W[:, R.shape[1] // 2:].mean(axis=1), fit=R[:, R.shape[1] // 2:, 7].mean(), almax=float(np.abs(al).max()),
               end=W[:, -1])
    for t in T_E:
        i = int(np.argmin(np.abs(s - t)))
        out['times'].append(dict(s=s[i], mean=W[:, i].mean(axis=0), se=W[:, i].std(axis=0, ddof=1) / np.sqrt(len(W)), std=W[:, i].std(axis=0, ddof=1)))
    return out


def vertices():
    return np.array([[1 - W_BOX[1] - W_BOX[2], W_BOX[1], W_BOX[2]], [W_BOX[0], 1 - W_BOX[0] - W_BOX[2], W_BOX[2]],
                     [W_BOX[0], W_BOX[1], 1 - W_BOX[0] - W_BOX[1]]])


def report(A=None, B=None, C=None, D=None, Ee=None):
    fmt = E.fmt
    v3 = lambda v, p=3: ', '.join(f'{x:.{p}f}' for x in v)
    L_ = ["# Proposition D‴.39: dissipativity of the dithered gauge (`exp_gauge_plateau.py`)", "",
          f"Proposition D‴.36's instance (three bumps, one particle per class, m = 3, the cohesion leak); energies over τ = {TAU:g}; the box "
          f"L = ({v3(BOX, 1)}) in (α, γ, κ), σ_B = {SIGMA_B:g}: |a| ≤ {A_MAX * np.tanh(BOX[0] / A_MAX):.3f}, the representable masses "
          f"W_k ≥ {v3(W_BOX, 4)} (the saturation's {v3(W_MIN, 4)}).", ""]
    if A is not None:
        L_ += ["## A. The plateau: the lowest configuration with one amplitude fixed far out", "",
               "| masses | class | α_k | ΔJ/τ above the exact fit | fit error | a | c | s |", "|---|---|---|---|---|---|---|---|"]
        for r in A:
            L_.append(f"| {r['which']} | {r['k'] + 1} | {dict(sat=f'±{AL_SAT:g}', esc=f'±{ESC:g}', box='±L_α (box)')[r['at']]} | {r['dJ']:.3f} | {r['fit']:.4f} | {v3(r['a'], 2)} | "
                      f"{v3(r['c'], 2)} | {v3(r['s'], 2)} |")
        L_.append("")
    if B is not None:
        L_ += ["## B. The minimum-energy path from the exact fit to the plateau (string method)", ""]
        for r in B:
            p = r['prof']
            s2 = [E.setup()[0].physical(T)[2][1] for T in r['images'][::4]]
            L_ += [f"Masses {r['which']}, toward class {r['k'] + 1}'s saturation: the largest ΔJ/τ along the path {p.max():.3f} (image {int(p.argmax())} of "
                   f"{len(p) - 1}); at every fourth image the profile {', '.join(f'{v:.2f}' for v in p[::4])}, the fit error "
                   f"{', '.join(f'{v:.4f}' for v in r['fits'][::4])}, the second width {', '.join(f'{v:.2f}' for v in s2)}.", ""]
    if C is not None:
        L_ += [f"## C. The escapes at fixed masses 1/3 (γ = 0; {len(SEEDS_C)} seeds, at most {S_C:g} slow units; an escape: max|α_k| > {ESC:g})", "",
               "| τ | escapes | the rate (per slow unit, ± s.e.) | the median time | the rate × e^{ΔJ/τ} | the fit error at the escape |", "|---|---|---|---|---|---|"]
        for r in C['rows']:
            L_.append(f"| {r['tau']:g} | {r['n']}/{r['N']} | {fmt(r['rate'])} ± {fmt(r['rate_se'])} | {r['med']:.3g} | {fmt(r['scaled'])} | {fmt(r['fit'])} |")
        L_ += ["", f"The slope of log(rate) against 1/τ: {fmt(-C['slope'])} (the lowest J on the escape set's boundary, ΔJ = {fmt(C['dJ'])}).", ""]
    if D is not None:
        L_ += [f"## D. The dead zone at fixed masses 1/3 (γ = 0, τ = {TAU_D:g}, {len(SEEDS_D)} seeds, {S_D:g} slow units)", "",
               f"A record is on the plateau when its second width exceeds {S2_PLATEAU:g} (A and B: the plateaus' configurations and the path's).", "",
               f"| box | seeds that reach the plateau | fraction of the records on it (s_2 > {S2_PLATEAU:g}), by quarters | fraction with max|α| > {ESC:g} | "
               "the largest |α|, by quarters | the fit error, by quarters |", "|---|---|---|---|---|---|"]
        for r in D:
            L_.append(f"| {'L = (' + v3(BOX, 1) + ')' if r['boxed'] else 'none'} | {r['seeds']}/{r['N']} | {', '.join(f'{v:.3f}' for v in r['plat'])} | "
                      f"{', '.join(f'{v:.3f}' for v in r['frac'])} | "
                      f"{', '.join(f'{v:.1f}' for v in r['almax'])} | {', '.join(f'{v:.4f}' for v in r['fit'])} |")
        L_.append("")
    if Ee is not None:
        ws = Q.w_sharp(); V = vertices()
        L_ += [f"## E. The full loop with the box (γ = {GAMMA_E:g}, τ = {TAU:g}, {len(SEEDS_E)} seeds, {S_E:g} slow units)", "",
               f"The box's vertices: {'; '.join(v3(v, 4) for v in V)}; W♯ = {v3(ws, 4)}.", "",
               "| κ_F | s | mean of W_1, W_2, W_3 (± s.e.) | std |", "|---|---|---|---|"]
        for r in Ee:
            for t in r['sum']['times']:
                L_.append(f"| {r['kF']:g} | {t['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(t['mean'], t['se']))} | {v3(t['std'], 4)} |")
        L_ += ["", "| κ_F | the masses at the end, by seed (the nearest vertex for κ_F = 0, its distance) | the largest |α| | the fit error (second half) |",
               "|---|---|---|---|"]
        for r in Ee:
            sm = r['sum']
            if r['kF'] == 0:
                d = np.linalg.norm(sm['end'][:, None, :] - V[None, :, :], axis=2)
                cell = '; '.join(f"{v3(w, 3)} (v{int(np.argmin(dd)) + 1}, {dd.min():.3f})" for w, dd in zip(sm['end'], d))
            else:
                cell = '; '.join(f"{v3(w, 3)} ({np.abs(w - ws).max():.3f})" for w in sm['end'])
            L_.append(f"| {r['kF']:g} | {cell} | {sm['almax']:.2f} | {fmt(sm['fit'])} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_gauge_plateau.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_plateau_raw.npz')
    md = os.path.join(OUT, 'gauge_plateau.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDE']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDE')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'A' in parts:
            res['A'] = pool.map(plateau_job, [(w, k, at) for at in ('sat', 'esc', 'box') for w in ('uniform', 'sharp') for k in range(3)], chunksize=1)
            E.save(raw, **res); log('A')
        if 'E' in parts:
            h['E'] = pool.map_async(series, [(GAMMA_E, TAU, k, f, True, s, S_E) for k, f in ((0.0, 'none'), (2.0, 'masses')) for s in SEEDS_E], chunksize=1)
        if 'B' in parts:
            ends = []
            for w in ('uniform', 'sharp'):
                rr = [r for r in res['A'] if r['which'] == w and r['at'] == 'sat']
                best = min(rr, key=lambda r: r['dJ']); ends.append((w, best['k'], best['Th']))
            h['B'] = pool.map_async(string_job, ends, chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(series, [(0.0, TAU_D, 0.0, 'none', b, s, S_D) for b in (False, True) for s in SEEDS_D], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(escape_job, [(t, s) for t in TAUS_C for s in SEEDS_C], chunksize=1)
        if 'B' in h:
            res['B'] = h['B'].get(); E.save(raw, **res); log('B')
        if 'C' in h:
            dJ = min(r['dJ'] for r in res['A'] if r['which'] == 'uniform' and r['at'] == 'esc')
            res['C_runs'] = h['C'].get(); res['C'] = part_C(res['C_runs'], dJ); E.save(raw, **res); log('C')
        if 'D' in h:
            res['D_rec'] = h['D'].get(); res['D'] = part_D(res['D_rec']); E.save(raw, **res); log('D')
        if 'E' in h:
            runs = h['E'].get()
            res['E_rec'] = {k: np.array([r['rec'] for r in runs if r['kF'] == k]) for k in (0.0, 2.0)}
            res['E'] = [dict(kF=k, sum=summarise_E(res['E_rec'][k])) for k in (0.0, 2.0)]
            E.save(raw, **res); log('E')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE']))


if __name__ == '__main__':
    main()
