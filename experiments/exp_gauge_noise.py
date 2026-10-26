"""
Proposition D‴.40: the gauge under general noise — only each amplitude's own noise intensity drives its class's mass
====================================================================================================================

Proposition D‴.36's instance (Proposition D‴.32's three bumps, one particle per class, m = 3). With white noise of covariance Σ (per unit slow
time) on the positions, the masses' fitness averaged over the positions' stationary fluctuations at fixed masses is, exactly,
ū_k = −τ_k E[g′_k(α_k)], τ_k := Σ_{α_kα_k}/(2G_α), g_k = a_k/(β_k a′_k) (the generator identity through the gauge identity
u_k = −g_k∂_{α_k}J), to leading order ū_k = −(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k): the replicator flow of the weighted log-determinant
Σ_kτ_k log(W_k a′_k).

A. The trace formula (`gauge_noise.averaged_fitness`: ū = ∇u·μ + ½tr(∇²u C), C the Lyapunov solution, μ the mean shift) against the
   closed form, for isotropic noise, noise on the amplitudes only, on the centres and widths only, per-class amplitude levels and three
   random correlated covariances, at three points of M.
B. Katzenberger's limit drift at γ = 10⁻⁵ (`gauge_drift.limit_drift` with the noise Σ) against the replicator flow of the closed form,
   for the same noises at two points; and its diffusion.
C. Noise on the centres and widths only at γ = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵: the masses' diffusion (zero: the gauge's tangent lies in the
   amplitude–mass coordinates and the metric is diagonal) and drift (O(γ)), against the isotropic noise's.
D. The identity in simulation at fixed masses 1/3 (γ = 0), inside the box of Proposition D‴.39 (so the stationary law exists): the time
   averages of u_k and of −τ_k g′_k(α_k) + g_k∂_{α_k}(J_coh + V) (the cohesion's and the dead zone's parts of the amplitude's force,
   which the gauge identity does not cover) over S_D slow units, for isotropic and correlated noise at
   τ = 10⁻³ and 10⁻², where the small-noise expansion is poor; their difference against its batch-means error; the amplitudes' noise
   intensity read from their quadratic variation, τ̂_k = Σ(Δα_k)²/(2G_αS).
E. The full loop at γ = 10⁻², τ = 10⁻⁴ (the 12 seeds of exp_gauge_drift.py's Part H, 10⁶ slow units): the centres and widths dithered only
   (the masses frozen); the same with the barrier of Proposition D‴.38's masses form at κ_F = 1 (no amplitude noise to cancel: the masses
   ascend ℓ to W♯); correlated noise on all coordinates ('correlated 1'), uncorrected, and with the masses form at the classes' levels τ_k
   and κ_F = 2 (the masses ascend Σ_kτ_k log(W_k a′_k) to its maximiser W♯_τ).

Output: results/chaos/gauge_noise.md and gauge_noise_raw.npz (GN_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`).
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
from src.core import gauge_drift as GDR
from src.core import gauge_noise as GN
import exp_cohesion_leak as E
import exp_gauge_dither as GD
import exp_gauge_drift as X
import exp_gauge_fixman as Q
import exp_gauge_plateau as P

OUT = os.environ.get('GN_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GN_WORKERS', 4))
TAU, A_MAX, DT, REC = GD.TAU, 4.0, GD.DT, GD.REC
G9 = np.tile(E.G, 3)
W_A = ((1 / 3, 1 / 3), (0.30, 0.35), (0.40, 0.25))
W_B = ((1 / 3, 1 / 3), (0.30, 0.35))
GAMMAS_C = (1e-2, 1e-3, 1e-4, 1e-5)
TAUS_D, S_D, SEED_D, N_BATCH = (1e-3, 1e-2), 2e5, 700, 20
GAMMA_E, SEEDS_E, S_E, T_E = 1e-2, X.SEEDS_H, 1e6, X.T_H
NOISE_SEEDS = (11, 12, 13)


def noises(tau=TAU):
    """name → Σ (9 × 9, per unit slow time) on the flattened positions (α, γ, κ per particle)."""
    out = {'isotropic': 2 * tau * np.diag(G9),
           'amplitudes only': 2 * tau * np.diag(np.tile([1.0, 0.0, 0.0], 3)),
           'centres and widths only': 2 * tau * np.diag(np.tile([0.0, E.G[1], E.G[2]], 3)),
           'amplitude levels (2, 1, 0.5)τ': 2 * tau * np.diag(np.repeat([2.0, 1.0, 0.5], 3) * np.tile([1.0, 0.0, 0.0], 3))}
    for i, sd in enumerate(NOISE_SEEDS):
        Lr = np.random.default_rng(sd).standard_normal((9, 9))
        out[f'correlated {i + 1}'] = 2 * tau * 0.3 * (Lr @ Lr.T)
    return out


def levels(Sig):
    """τ_k = Σ_{α_kα_k}/(2G_α)."""
    return np.diag(Sig)[0::3] / (2 * E.G[0])


def closed(W, tau_k):
    """−(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k) at the exact fit with the masses W."""
    W = np.asarray(W, float); q = (np.abs(E.A_BUMP) / W / A_MAX) ** 2
    return -(np.asarray(tau_k, float) / 3) * (1 + q) / (W * (1 - q))


def w_sharp_weighted(tk):
    """argmax of Σ_kτ_k log(W_k − W_min,k²/W_k) on the simplex (Lagrange: τ_k(1 + W_min,k²/W_k²)/(W_k − W_min,k²/W_k) equal) — the
    representable masses' analytic centre for the classes' noise levels; Q.w_sharp() for equal levels."""
    from scipy.optimize import brentq
    Wm = np.abs(E.A_BUMP) / A_MAX; tk = np.asarray(tk, float) / np.mean(tk)
    g = lambda W, k: tk[k] * (1 + Wm[k] ** 2 / W ** 2) / (W - Wm[k] ** 2 / W)        # decreasing in W
    W_of = lambda mu, k: brentq(lambda W: g(W, k) - mu, Wm[k] * (1 + 1e-12), 1.0)
    mu = brentq(lambda mu: sum(W_of(mu, k) for k in range(3)) - 1, max(g(1.0, k) for k in range(3)), 1e12)
    return np.array([W_of(mu, k) for k in range(3)])


def a_job(args):
    name, w = args
    W = np.r_[w, 1 - sum(w)]
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    beta = 3 * W
    gradJ = lambda t: -loop.drift(t.reshape(3, 3), beta)[0].ravel() / G9
    u = lambda t: loop.fitness(loop.residual(t.reshape(3, 3), beta), t.reshape(3, 3))
    Sig = noises()[name]
    ub = GN.averaged_fitness(gradJ, u, E.exact_rep(net, W).ravel(), Sig, G9)[0]
    cf = closed(W, levels(Sig))
    flow = lambda x: GN.replicator_flow(W, x)
    return dict(name=name, w=w, ubar=ub, closed=cf, rel=np.abs(flow(ub) - flow(cf)).max() / max(np.abs(flow(cf)).max(), 1e-30),
                flow=flow(ub) / TAU)


def b_job(args):
    name, w = args
    Sig = noises()[name]
    F, chart, coords, S = X.instance(1e-5)
    S11 = np.zeros((11, 11)); S11[:9, :9] = Sig
    b, D, _ = GDR.limit_drift(F, chart, coords, np.asarray(w, float), S11, h=X.H_ID, h2=X.H2_ID)
    W = np.r_[w, 1 - sum(w)]
    pred = GN.replicator_flow(W, closed(W, levels(Sig))) * 1e-5
    return dict(name=name, w=w, b=np.r_[b, -b.sum()], pred=pred, D=D)


def c_job(args):
    gamma, w = args
    F, chart, coords, S = X.instance(gamma)
    out = {}
    for name in ('isotropic', 'centres and widths only'):
        S11 = np.zeros((11, 11)); S11[:9, :9] = noises()[name]
        b, D, _ = GDR.limit_drift(F, chart, coords, np.asarray(w, float), S11, h=X.H_ID, h2=X.H2_ID)
        out[name] = dict(b=np.r_[b, -b.sum()] / (gamma * TAU), D=D / (gamma ** 2 * TAU))
    return dict(gamma=gamma, w=w, **out)


def d_job(args):
    """Fixed masses 1/3, the box, the noise Σ(τ): the time averages (by batches) of u_k and of −τ_k g′_k(α_k) + g_k ∂_{α_k}V, and the
    quadratic variation of the amplitudes."""
    name, tau = args
    t0 = time.time()
    net, base = E.setup()
    W = np.ones(3) / 3; beta = 3 * W
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    Sig = noises(tau)[name]; tk = levels(Sig)
    lam, V = np.linalg.eigh(Sig); Lc = V * np.sqrt(np.clip(lam, 0, None))            # Lc Lcᵀ = Σ
    n = int(round(S_D / DT)); nb = n // N_BATCH
    U = np.zeros((N_BATCH, 3)); Rhs = np.zeros((N_BATCH, 3)); qv = np.zeros(3)
    Th = E.exact_rep(net, W).copy(); rng = np.random.default_rng(SEED_D)
    sq = np.sqrt(DT) * Lc
    for k in range(nb * N_BATCH):
        a = base.features(Th)[0]
        ap = 1 - (a / A_MAX) ** 2
        g = a / (beta * ap); gp = (1 + (a / A_MAX) ** 2) / (beta * ap)
        r = loop.residual(Th, beta)
        u = loop.fitness(r, Th)
        dT = loop.drift(Th, beta)[0]
        nonfit = -(dT[:, 0] - beta * E.G[0] * base.grad_weighted(r, Th)[:, 0]) / E.G[0]   # ∂_{α_k}(J_coh + V): the cohesion and the dead zone
        i = k // nb
        U[i] += u; Rhs[i] += -tk * gp + g * nonfit
        dW = (sq @ rng.standard_normal(9)).reshape(3, 3)
        Th = Th + DT * dT + dW
        qv += (DT * dT[:, 0] + dW[:, 0]) ** 2
    U /= nb; Rhs /= nb
    print(f"[{time.time() - t0:.0f}s] D {name} tau {tau:g}", flush=True)
    return dict(name=name, tau=tau, tk=tk, U=U, Rhs=Rhs, tk_hat=qv / (2 * E.G[0] * nb * N_BATCH * DT), closed=closed(W, tk))


def e_series(args):
    """The full loop: (variant, seed) → records (s, α_1..3, W_1..3, the fit error)."""
    variant, seed = args
    t0 = time.time()
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3).copy()
    A0 = base.features(Th)[0] / 3
    Sig = noises()['correlated 1']
    kw = dict(mask=None, noise_L=None)
    if variant == 'nonlinear only':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E); kw['mask'] = (0.0, 1.0, 1.0)
    elif variant == 'nonlinear only + barrier':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, fixman=1.0, fixman_tau=TAU, fixman_A=A0); kw['mask'] = (0.0, 1.0, 1.0)
    elif variant == 'correlated':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E); kw['noise_L'] = np.linalg.cholesky(Sig)
    elif variant == 'correlated + class levels':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, fixman=2.0, fixman_tau=levels(Sig), fixman_A=A0)
        kw['noise_L'] = np.linalg.cholesky(Sig)
    else:
        raise ValueError(variant)
    rec = []
    loop.run_dithered(Th, np.ones(3), S_E, DT, TAU, np.random.default_rng(seed), rec_every=REC,
                      callback=lambda s_, T, b: rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b)]), **kw)
    print(f"[{time.time() - t0:.0f}s] E {variant} seed {seed}", flush=True)
    return dict(variant=variant, seed=seed, rec=np.array(rec))


VARIANTS_E = ('nonlinear only', 'nonlinear only + barrier', 'correlated', 'correlated + class levels')


def summarise_E(R):
    s = R[0, :, 0]; W = R[:, :, 4:7]
    out = dict(times=[], half=W[:, R.shape[1] // 2:].mean(axis=1), fit=R[:, R.shape[1] // 2:, 7].mean(), almax=float(np.abs(R[:, :, 1:4]).max()))
    for t in T_E:
        i = int(np.argmin(np.abs(s - t)))
        out['times'].append(dict(s=s[i], mean=W[:, i].mean(axis=0), se=W[:, i].std(axis=0, ddof=1) / np.sqrt(len(W)), std=W[:, i].std(axis=0, ddof=1)))
    return out


def report(A=None, B=None, C=None, D=None, Ee=None):
    fmt = E.fmt
    v3 = lambda v, p=4: ', '.join(f'{x:.{p}f}' for x in v)
    L_ = ["# Proposition D‴.40: the gauge under general noise (`exp_gauge_noise.py`)", "",
          f"Proposition D‴.36's instance (three bumps, one particle per class, m = 3); the noises at τ = {TAU:g}: isotropic 2τ diag(G); the "
          "amplitudes only; the centres and widths only; amplitude levels (2, 1, 0.5)τ; three correlated covariances 2τ·0.3·LLᵀ (L standard "
          f"normal, seeds {', '.join(str(s) for s in NOISE_SEEDS)}). τ_k = Σ_{{α_kα_k}}/(2G_α); the closed form −(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k).", ""]
    if A is not None:
        L_ += ["## A. The trace formula against the closed form", "",
               "| noise | W | τ_k / τ | ū (trace formula) | the closed form | their replicator flows' largest relative difference |", "|---|---|---|---|---|---|"]
        for r in A:
            zero = not np.any(levels(noises()[r['name']]))
            L_.append(f"| {r['name']} | {v3(np.r_[r['w'], 1 - sum(r['w'])], 3)} | {', '.join(f'{v:.3f}' for v in levels(noises()[r['name']]) / TAU)} | "
                      f"{', '.join(fmt(v) for v in r['ubar'])} | {', '.join(fmt(v) for v in r['closed'])} | "
                      f"{(fmt(np.abs(r['flow']).max() * TAU) + ' (the flow itself; the closed form is 0)') if zero else fmt(r['rel'])} |")
        L_.append("")
    if B is not None:
        L_ += ["## B. Katzenberger's limit drift at γ = 10⁻⁵ against the replicator flow of the closed form", "",
               "| noise | W | the drift per unit slow time | the closed form's flow (γ = 10⁻⁵) | largest relative difference | the diffusion's diagonal |",
               "|---|---|---|---|---|---|"]
        for r in B:
            den = max(np.abs(r['pred']).max(), 1e-30)
            L_.append(f"| {r['name']} | {v3(np.r_[r['w'], 1 - sum(r['w'])], 3)} | {', '.join(fmt(v) for v in r['b'])} | {', '.join(fmt(v) for v in r['pred'])} | "
                      f"{fmt(np.abs(r['b'] - r['pred']).max() / den) if np.abs(r['pred']).max() > 0 else fmt(np.abs(r['b']).max())} | "
                      f"{', '.join(fmt(v) for v in np.diag(r['D']))} |")
        L_.append("")
    if C is not None:
        L_ += ["## C. The centres and widths dithered only, at finite γ", "",
               "| γ | W | drift/γτ, isotropic | drift/γτ, centres and widths only | D/γ²τ, isotropic (diagonal) | D/γ²τ, centres and widths only |",
               "|---|---|---|---|---|---|"]
        for r in C:
            iso, cs = r['isotropic'], r['centres and widths only']
            L_.append(f"| {r['gamma']:g} | {v3(np.r_[r['w'], 1 - sum(r['w'])], 3)} | {', '.join(f'{v:+.5f}' for v in iso['b'])} | {', '.join(f'{v:+.2e}' for v in cs['b'])} | "
                      f"{', '.join(f'{v:.4f}' for v in np.diag(iso['D']))} | {np.abs(cs['D']).max():.1e} |")
        L_.append("")
    if D is not None:
        L_ += [f"## D. The identity in simulation (fixed masses 1/3, the box, {S_D:g} slow units, {N_BATCH} batches)", "",
               "| noise | τ | E[u_k] (± s.e.) | E[−τ_kg′_k(α_k) + g_k∂_{α_k}(J_coh + V)] (± s.e.) | the difference (± s.e.) | the closed form at the fit | τ_k / τ̂_k (quadratic variation) |",
               "|---|---|---|---|---|---|---|"]
        for r in D:
            se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M))
            U, Rh = r['U'], r['Rhs']
            L_.append(f"| {r['name']} | {r['tau']:g} | {', '.join(f'{a:+.3e} ± {b:.1e}' for a, b in zip(U.mean(axis=0), se(U)))} | "
                      f"{', '.join(f'{a:+.3e} ± {b:.1e}' for a, b in zip(Rh.mean(axis=0), se(Rh)))} | "
                      f"{', '.join(f'{a:+.1e} ± {b:.1e}' for a, b in zip((U - Rh).mean(axis=0), se(U - Rh)))} | {', '.join(f'{v:+.3e}' for v in r['closed'])} | "
                      f"{', '.join(f'{a:.4g} / {b:.4g}' for a, b in zip(r['tk'], r['tk_hat']))} |")
        L_.append("")
    if Ee is not None:
        ws = Q.w_sharp(); wt = w_sharp_weighted(levels(noises()['correlated 1']))
        target = {'nonlinear only': np.ones(3) / 3, 'nonlinear only + barrier': ws, 'correlated': None, 'correlated + class levels': wt}
        L_ += [f"## E. The full loop (γ = {GAMMA_E:g}, τ = {TAU:g}, {len(SEEDS_E)} seeds, {S_E:g} slow units)", "",
               f"Targets: the initial masses 1/3 (frozen), W♯ = {v3(ws)}, W♯_τ = {v3(wt)} (the maximiser of Σ_kτ_k log(W_k a′_k) for "
               f"'correlated 1', τ_k/τ = {', '.join(f'{v:.3f}' for v in levels(noises()['correlated 1']) / TAU)}).", "",
               "| variant | s | mean of W_1, W_2, W_3 (± s.e.) | std |", "|---|---|---|---|"]
        for r in Ee:
            for t in r['sum']['times']:
                L_.append(f"| {r['variant']} | {t['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(t['mean'], t['se']))} | {v3(t['std'])} |")
        L_ += ["", "| variant | the second half's mean (± s.e.) | its target | the distance from it (largest component) | the largest |α| | "
               "the fit error (second half) |", "|---|---|---|---|---|---|"]
        for r in Ee:
            hm = r['sum']['half']; se = hm.std(axis=0, ddof=1) / np.sqrt(len(hm)); tg = target[r['variant']]
            dist = '—' if tg is None else f"{np.abs(hm.mean(axis=0) - tg).max():.4f}"
            L_.append(f"| {r['variant']} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(hm.mean(axis=0), se))} | "
                      f"{'none' if tg is None else v3(tg)} | {dist} | {r['sum']['almax']:.2f} | {fmt(r['sum']['fit'])} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_gauge_noise.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_noise_raw.npz')
    md = os.path.join(OUT, 'gauge_noise.md')
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
    names = list(noises())
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'E' in parts:
            h['E'] = pool.map_async(e_series, [(v, s) for v in VARIANTS_E for s in SEEDS_E], chunksize=1)
        if 'A' in parts:
            h['A'] = pool.map_async(a_job, [(n, w) for n in names for w in W_A], chunksize=1)
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(n, w) for n in names for w in W_B], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(g, w) for g in GAMMAS_C for w in W_B], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, [(n, t) for n in ('isotropic', 'correlated 1') for t in TAUS_D], chunksize=1)
        for k in 'ABCD':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
        if 'E' in h:
            runs = h['E'].get()
            res['E_rec'] = {v: np.array([r['rec'] for r in runs if r['variant'] == v]) for v in VARIANTS_E}
            res['E'] = [dict(variant=v, sum=summarise_E(res['E_rec'][v])) for v in VARIANTS_E]
            E.save(raw, **res); log('E')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE']))


if __name__ == '__main__':
    main()
