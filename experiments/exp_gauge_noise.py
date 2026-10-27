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
   averages of u_k and of each term of the generator identity E[u_k] = −τ_kE[g′_k] + E[g_k∂_{α_k}J_coh] + E[g_k∂_{α_k}V] + E[Lψ_k] (ψ_k =
   ∫g_k dα_k/G_α; the last term the Euler chain's bias, measured exactly from the Gaussian increment's closed-form expectation), with
   (i)'s box-free residual E[u_k + τ_kg′_k] and its z, the share of steps outside the box, and the whole identity's residual against its
   martingale error; isotropic noise, the amplitudes only (the same τ_k, no centre or width noise), 'correlated 1' at τ = 10⁻³ and 10⁻²,
   and the centres and widths only at 10⁻²; the amplitudes' noise intensity from their quadratic variation, with and without the drift's
   part (referee round 31).
E. The full loop inside the box at γ = 10⁻², τ = 10⁻⁴ (the 12 seeds of exp_gauge_drift.py's Part H, 10⁶ slow units): the centres and
   widths dithered only (the masses frozen); the same with the barrier of Proposition D‴.38's masses form at κ_F = 1 (no amplitude noise
   to cancel: the masses ascend ℓ to W♯); a correlated noise on all coordinates whose amplitudes have the levels (2, 1, 0.5)τ and whose
   centres and widths stay at the dither's level ('correlated E'), uncorrected, and with the masses form at the classes' levels τ_k and
   κ_F = 2 (the masses ascend Σ_kτ_k log(W_k a′_k) to its maximiser W♯_τ, 0.084 from W♯).

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
from src.core.runcache import cached
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
TAUS_D, S_D, SEED_D, N_BATCH, N_BATCH_LONG = (1e-3, 1e-2), 2e5, 700, 40, 10
RUNS_D = tuple((n, t) for n in ('isotropic', 'amplitudes only', 'correlated 1') for t in TAUS_D) + (('centres and widths only', 1e-2),)
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
    out['correlated E'] = 2 * tau * corr_E()
    return out


LEVELS_E = (2.0, 1.0, 0.5)


def corr_E():
    """Part E's correlated noise over 2τ: the amplitudes at the levels LEVELS_E·G_α, the centres and widths at the dither's level G_γ, G_κ,
    correlated through the mixture ½I + ½R of a random correlation matrix R (seed 21) — so the amplitudes' levels differ and every
    coordinate's noise is correlated, while the centres and widths are not heated beyond the dither (referee round 31)."""
    Lr = np.random.default_rng(21).standard_normal((9, 9))
    C = Lr @ Lr.T; d = np.sqrt(np.diag(C))
    R = 0.5 * np.eye(9) + 0.5 * C / np.outer(d, d)
    sd = np.sqrt(np.tile(E.G, 3) * np.ravel([[lv, 1.0, 1.0] for lv in LEVELS_E]))
    return R * np.outer(sd, sd)


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


@cached
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


@cached
def b_job(args):
    name, w = args
    Sig = noises()[name]
    F, chart, coords, S = X.instance(1e-5)
    S11 = np.zeros((11, 11)); S11[:9, :9] = Sig
    b, D, _ = GDR.limit_drift(F, chart, coords, np.asarray(w, float), S11, h=X.H_ID, h2=X.H2_ID)
    W = np.r_[w, 1 - sum(w)]
    pred = GN.replicator_flow(W, closed(W, levels(Sig))) * 1e-5
    return dict(name=name, w=w, b=np.r_[b, -b.sum()], pred=pred, D=D)


@cached
def c_job(args):
    gamma, w = args
    F, chart, coords, S = X.instance(gamma)
    out = {}
    for name in ('isotropic', 'centres and widths only'):
        S11 = np.zeros((11, 11)); S11[:9, :9] = noises()[name]
        b, D, _ = GDR.limit_drift(F, chart, coords, np.asarray(w, float), S11, h=X.H_ID, h2=X.H2_ID)
        out[name] = dict(b=np.r_[b, -b.sum()] / (gamma * TAU), D=D / (gamma ** 2 * TAU))
    return dict(gamma=gamma, w=w, **out)


@cached
def d_job(args):
    """Fixed masses 1/3, the box, the noise Σ(τ): per batch the time averages of u_k and of the identity's terms −τ_kg′_k(α_k),
    g_k∂_{α_k}J_coh, g_k∂_{α_k}V (the box), Lψ_k (ψ_k = ∫g_k dα_k/G_α = (a_max²/4β_kG_α)cosh(2α_k/a_max), the generator applied to it: its
    stationary mean is the Euler chain's bias) and the chain's exact per-step bias E_n[ψ_k(α_{n+1}) − ψ_k(α_n)]/Δt − Lψ_k (the Gaussian
    increment's expectation in closed form); the martingale increments' sum of squares (the exact error of E[Lψ_k]); the share of steps with
    α_k outside the box; the amplitudes' quadratic variation with and without the drift's part."""
    name, tau = args
    t0 = time.time()
    net, base = E.setup()
    W = np.ones(3) / 3; beta = 3 * W
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    Sig = noises(tau)[name]; tk = levels(Sig)
    lam, V = np.linalg.eigh(Sig); Lc = V * np.sqrt(np.clip(lam, 0, None))            # Lc Lcᵀ = Σ
    n = int(round(S_D / DT)); nb = n // N_BATCH
    keys = ('u', 'gp', 'coh', 'box', 'lpsi', 'bias')
    acc = {k: np.zeros((N_BATCH, 3)) for k in keys}
    m2, out_box, qv, qv_noise = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)
    Th = E.exact_rep(net, W).copy(); rng = np.random.default_rng(SEED_D)
    sq = np.sqrt(DT) * Lc
    kk, cpsi, s2 = 2.0 / A_MAX, A_MAX ** 2 / (4 * beta * E.G[0]), np.diag(Sig)[0::3] * DT
    psi0 = cpsi * np.cosh(kk * Th[:, 0])
    for k in range(nb * N_BATCH):
        a = base.features(Th)[0]
        ap = 1 - (a / A_MAX) ** 2
        g = a / (beta * ap); gp = (1 + (a / A_MAX) ** 2) / (beta * ap)
        r = loop.residual(Th, beta)
        u = loop.fitness(r, Th)
        dT = loop.drift(Th, beta)[0]
        nonfit = -(dT[:, 0] - beta * E.G[0] * base.grad_weighted(r, Th)[:, 0]) / E.G[0]   # ∂_{α_k}(J_coh + V)
        boxt = P.SIGMA_B * loop.box_excess(Th)[:, 0] / E.G[0]                            # ∂_{α_k}V
        lpsi = u - g * nonfit + tk * gp                                                   # Lψ_k = (g_k/G_α)b_k + τ_kg′_k
        Epsi = cpsi * np.cosh(kk * (Th[:, 0] + DT * dT[:, 0])) * np.exp(0.5 * kk * kk * s2)
        i = k // nb
        acc['u'][i] += u; acc['gp'][i] += -tk * gp; acc['coh'][i] += g * (nonfit - boxt); acc['box'][i] += g * boxt
        acc['lpsi'][i] += lpsi; acc['bias'][i] += (Epsi - psi0) / DT - lpsi
        out_box += np.abs(Th[:, 0]) > P.BOX[0]
        dW = (sq @ rng.standard_normal(9)).reshape(3, 3)
        qv += (DT * dT[:, 0] + dW[:, 0]) ** 2; qv_noise += dW[:, 0] ** 2
        Th = Th + DT * dT + dW
        psi1 = cpsi * np.cosh(kk * Th[:, 0])
        m2 += (psi1 - Epsi) ** 2
        psi0 = psi1
    for k_ in keys:
        acc[k_] /= nb
    T = nb * N_BATCH * DT
    print(f"[{time.time() - t0:.0f}s] D {name} tau {tau:g}", flush=True)
    return dict(name=name, tau=tau, tk=tk, **acc, se_mart=np.sqrt(m2) / T, out_box=out_box / (nb * N_BATCH),
                tk_hat=qv / (2 * E.G[0] * T), tk_hat_noise=qv_noise / (2 * E.G[0] * T), closed=closed(W, tk))


@cached
def e_series(args):
    """The full loop inside the box of Proposition D‴.39 (so that a stationary law exists): (variant, seed) → records (s, α_1..3, W_1..3,
    the fit error)."""
    variant, seed = args
    t0 = time.time()
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3).copy()
    A0 = base.features(Th)[0] / 3
    Sig = noises()['correlated E']
    kw = dict(mask=None, noise_L=None)
    bx = dict(box=P.BOX, box_sigma=P.SIGMA_B)
    if variant == 'nonlinear only':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, **bx); kw['mask'] = (0.0, 1.0, 1.0)
    elif variant == 'nonlinear only + barrier':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, fixman=1.0, fixman_tau=TAU, fixman_A=A0, **bx)
        kw['mask'] = (0.0, 1.0, 1.0)
    elif variant == 'correlated E':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, **bx); kw['noise_L'] = np.linalg.cholesky(Sig)
    elif variant == 'correlated E + class levels':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, fixman=2.0, fixman_tau=levels(Sig), fixman_A=A0, **bx)
        kw['noise_L'] = np.linalg.cholesky(Sig)
    else:
        raise ValueError(variant)
    rec = []
    loop.run_dithered(Th, np.ones(3), S_E, DT, TAU, np.random.default_rng(seed), rec_every=REC,
                      callback=lambda s_, T, b: rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b)]), **kw)
    print(f"[{time.time() - t0:.0f}s] E {variant} seed {seed}", flush=True)
    return dict(variant=variant, seed=seed, rec=np.array(rec))


VARIANTS_E = ('nonlinear only', 'nonlinear only + barrier', 'correlated E', 'correlated E + class levels')


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
        se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M))
        lng = lambda M: M.reshape(N_BATCH_LONG, -1, 3).mean(axis=1)
        pm = lambda M: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(lng(M).mean(axis=0), se(lng(M))))
        L_ += [f"## D. The identity in simulation (fixed masses 1/3, the box of Proposition D‴.39, {S_D:g} slow units)", "",
               "The generator applied to ψ_k = ∫g_k dα_k/G_α = (a_max²/4β_kG_α)cosh(2α_k/a_max) gives, in the stationary law, "
               "E[u_k] = −τ_kE[g′_k] + E[g_k∂_{α_k}J_coh] + E[g_k∂_{α_k}V] + E[Lψ_k], the last the Euler chain's bias (zero for the diffusion); "
               "(i)'s box-free form leaves out the cohesion's and the box's terms. Errors: batch means over "
               f"{N_BATCH_LONG} batches of {S_D / N_BATCH_LONG:g} slow units (the plateau gives a slow tail: referee round 31).", "",
               "| noise | τ | τ_k/τ | E[u_k] | E[−τ_kg′_k] | E[g_k∂_{α_k}J_coh] | E[g_k∂_{α_k}V] (the box) | (i)'s residual E[u_k + τ_kg′_k] (z) | "
               "α_k outside the box (share of steps) |", "|---|---|---|---|---|---|---|---|---|"]
        zi = []
        for r in D:
            res = lng(r['u'] - r['gp']); z = res.mean(axis=0) / se(res); zi += list(z)
            L_.append(f"| {r['name']} | {r['tau']:g} | {', '.join(f'{v:.3f}' for v in r['tk'] / r['tau'])} | {pm(r['u'])} | {pm(r['gp'])} | "
                      f"{pm(r['coh'])} | {pm(r['box'])} | {', '.join(f'{a:+.2e} ({b:+.1f})' for a, b in zip(res.mean(axis=0), z))} | "
                      f"{', '.join(f'{v:.3f}' for v in r['out_box'])} |")
        L_ += ["", "The whole identity and the Euler chain's bias (the chain's exact per-step bias b_k = E_n[ψ_k(α_(n+1)) − ψ_k(α_n)]/Δt − Lψ_k, "
               "the Gaussian increment's expectation in closed form, so that E[Lψ_k] = −E[b_k] exactly; the error of E[Lψ_k] from the "
               "martingale increments ψ_k(α_(n+1)) − E_n[ψ_k(α_(n+1))]):", "",
               "| noise | τ | E[Lψ_k] (± its martingale error) | −E[b_k] | E[Lψ_k] + E[b_k] over that error | the bias E[b_k] over τ_kE[g′_k] | "
               "τ_k/τ̂_k: the amplitudes' quadratic variation, with the drift's part, without |", "|---|---|---|---|---|---|---|"]
        zf = []
        for r in D:
            lp, bi = r['lpsi'].mean(axis=0), r['bias'].mean(axis=0)
            z = (lp + bi) / r['se_mart']; zf += list(z)
            L_.append(f"| {r['name']} | {r['tau']:g} | {', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(lp, r['se_mart']))} | "
                      f"{', '.join(f'{-v:+.2e}' for v in bi)} | {', '.join(f'{v:+.2f}' for v in z)} | "
                      f"{', '.join(f'{v:+.4f}' for v in bi / -r['gp'].mean(axis=0))} | "
                      f"{', '.join(f'{a / b:.4f}, {a / c:.4f}' for a, b, c in zip(r['tk'], r['tk_hat'], r['tk_hat_noise']))} |")
        zi, zf = np.array(zi), np.array(zf)
        L_ += ["", f"The whole identity: χ² = {float(zf @ zf):.1f} for {len(zf)} values (largest |z| {np.abs(zf).max():.2f}). (i)'s box-free form: "
               f"the largest |z| {np.abs(zi).max():.1f}, its residual the box's term.", ""]
    if Ee is not None:
        ws = Q.w_sharp(); wt = w_sharp_weighted(levels(noises()['correlated E']))
        target = {'nonlinear only': np.ones(3) / 3, 'nonlinear only + barrier': ws, 'correlated E': None, 'correlated E + class levels': wt}
        L_ += [f"## E. The full loop inside the box (γ = {GAMMA_E:g}, τ = {TAU:g}, {len(SEEDS_E)} seeds, {S_E:g} slow units)", "",
               f"Targets: the initial masses 1/3 (frozen), W♯ = {v3(ws)}, W♯_τ = {v3(wt)} (the maximiser of Σ_kτ_k log(W_k a′_k) for "
               f"'correlated E', τ_k/τ = {', '.join(f'{v:.3f}' for v in levels(noises()['correlated E']) / TAU)}, its centres and widths at the "
               f"dither's level); the box's representable masses W_B = {v3(P.W_BOX)}.", "",
               "| variant | s | mean of W_1, W_2, W_3 (± s.e.) | std |", "|---|---|---|---|"]
        for r in Ee:
            for t in r['sum']['times']:
                L_.append(f"| {r['variant']} | {t['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(t['mean'], t['se']))} | {v3(t['std'])} |")
        L_ += ["", "| variant | the second half's mean (± s.e.) | its target | the distance from it (largest component) | the distance from W♯ | "
               "the largest |α| | the fit error (second half) |", "|---|---|---|---|---|---|---|"]
        for r in Ee:
            hm = r['sum']['half']; se = hm.std(axis=0, ddof=1) / np.sqrt(len(hm)); tg = target[r['variant']]
            dist = '—' if tg is None else f"{np.abs(hm.mean(axis=0) - tg).max():.4f}"
            L_.append(f"| {r['variant']} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(hm.mean(axis=0), se))} | "
                      f"{'none' if tg is None else v3(tg)} | {dist} | {np.abs(hm.mean(axis=0) - ws).max():.4f} | {r['sum']['almax']:.2f} | "
                      f"{fmt(r['sum']['fit'])} |")
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
        if 'A' in parts:
            h['A'] = pool.map_async(a_job, [(n, w) for n in names for w in W_A], chunksize=1)
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(n, w) for n in names for w in W_B], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(g, w) for g in GAMMAS_C for w in W_B], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, list(RUNS_D), chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(e_series, [(v, s) for v in VARIANTS_E for s in SEEDS_E], chunksize=1)
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
