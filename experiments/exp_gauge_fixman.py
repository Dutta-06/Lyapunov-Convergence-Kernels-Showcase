"""
Proposition D‴.38: a Fixman correction for the replicator under dither — the entropic drift cancelled, or reversed toward the D-optimal
(largest-determinant) exact fit; the saturation made non-attracting for the unfloored correction
=====================================================================================================================================

Proposition D‴.36's instance (Proposition D‴.32's three bumps, one particle per class at an exact representation, the replicator at γ, the
positions dithered at τ = 10⁻⁴). The correction adds κ_F(τ/2)∂_{β_k}ℓ to the replicator's fitness, ℓ = log det ∇²_ΘJ along the manifold M of
exact fits = 2Σ log(W_k − W_min,k²/W_k) + const (Proposition D‴.37), in one of two forms that agree on M: the plug-in form takes
W_min,k = W_k|a_k|/a_max at the current amplitudes, κ_F(τ/m)(1 + (a_k/a_max)²)/(W_k a′_k) (`ReplicatorCohesionLoop(fixman=κ_F,
fixman_tau=τ)`); the masses form takes W_min,k = |Â_k|/a_max from the classes' contributions Â_k, fixed at the initial exact fit (Â = A) or
tracked by a running average of W_k a_k (`fixman_A=Â`, `fixman_TA`). On the slow time scale an O(τ) perturbation F₁ of the drift adds
∂W·P·F₁ to Katzenberger's drift (`gauge_drift.limit_drift(F1=…)`), F₁ evaluated on M; here F₁ is the corrected flow minus the uncorrected one,
so Parts A–C hold for both forms.

A. The cancellation (κ_F = 1): the uncorrected drift over γτ against the log-determinant law, and the corrected drift over γτ, at three
   points for γ = 2·10⁻², …, 10⁻⁴.
B. The gain: the zero of the corrected field (Newton from W♯) and the eigenvalues of its Jacobian for κ_F = 1, 2 at γ = 10⁻², 10⁻³, 10⁻⁴,
   against W♯ = argmax ℓ on the simplex (ℓ = 2Σ log(W_k − W_min,k²/W_k), strictly concave).
C. The saturation: along W_2 = W_3 near the first class's saturation (W_1 − 0.25 = 0.002–0.004, γ = 10⁻², inside its layer), the
   exponent d log p/d log δ of the one-dimensional zero-flux density p ∝ D⁻¹exp∫2b/D of Katzenberger's limit (τ → 0) for κ_F = 0, 1, 2
   (predicted κ_F − 1).
D. The full loop over 10⁶ slow units at γ = 10⁻², the plug-in form (12 seeds, the seeds of exp_gauge_drift.py's Part H, whose records are the
   κ_F = 0 runs): κ_F = 1, 2 — the excursions into the saturation by quarters, the fraction of records beyond the saturation, the masses'
   means and spreads at 10⁵, 3·10⁵, 10⁶, the second half's mean against W♯ and against the deterministic mean path dW/ds = b_c(W) from 1/3
   (its second-half average), the rms fit error against equipartition's √(pτ/(m⟨F²⟩)), the effective number p_eff = m⟨F²⟩⟨e²⟩/τ.
E. As D with the masses form, Â fixed at the initial exact fit (the same seeds).
F. With the entropic leak at γ = 3·10⁻³, ε = 10⁻⁴ (16 seeds, 4·10⁵ slow units, the second half): no correction, and at κ_F = 2 the plug-in
   form, the masses form with Â fixed and with Â tracked (T_A = 10⁵) — the seeds whose amplitudes escape
   up the positions' plateau (some |α_k| > 10), the masses' range, the pinned state's shift (all seeds, and the seeds without an escape)
   against κ⁻¹b (b the corrected drift at w = 1/3, κ the leak's stiffness along M) and against the second order (below), the fit error.
G. With the entropic leak at γ = 10⁻², ε = 3·10⁻⁴ (inside Proposition D‴.37(iv)'s window; the 16 seeds of exp_gauge_drift.py's Part G, whose
   κ_F = 0 runs this repeats bit for bit), the masses form at κ_F = 1, 2: the shift against κ⁻¹b_c and against the second order, and the
   paired differences from κ_F = 0. The second order (referee round 28): the limit SDE dW = f(W)ds + D^{1/2}dB on M, f = b + P(F_leak + F_c),
   has the stationary mean E[W] = W₀ − J⁻¹·½∂²f[C] + …, W₀ the zero of f, J its Jacobian, JC + CJᵀ + D = 0 — the fluctuations acting
   through the curvature of f (mostly the leak's projected force, whose stiffness varies on the scale of the distance to the saturation).
H. The implemented floor: the denominator of the correction is floored at 0.01W_k, so the correction stops growing within
   δ_cap = W_min,1(1/√0.99 − 1) of the first class's saturation; the exponent of C at δ = 0.0015–0.0004 with the floor and without it,
   against −1 + κ_F·2δ/(0.01W_min,1) below δ_cap (κ_F − 1 without the floor).

Output: results/chaos/gauge_fixman.md and gauge_fixman_raw.npz (GF_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`;
`--reanalyse` with `--only` recomputes D–G's analysis from the stored records instead of running them).
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
from src.core.runcache import cached
import exp_cohesion_leak as E
import exp_gauge_dither as GD
import exp_gauge_drift as X

OUT = os.environ.get('GF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GF_WORKERS', 3))
TAU, A_MAX, DT = GD.TAU, 4.0, GD.DT
W_MIN = X.W_MIN
GAMMAS_A, W_A = (2e-2, 1e-2, 5e-3, 3e-3, 2e-3, 1e-3, 5e-4, 3e-4, 2e-4, 1e-4), ((1 / 3, 1 / 3), (0.40, 0.25), (0.36, 0.30))
GAMMAS_B, GAINS_B = (1e-2, 1e-3, 1e-4), (1.0, 2.0)
DELTAS_C, GAINS_C, GAMMA_C = (0.002, 0.0025, 0.003, 0.004), (0.0, 1.0, 2.0), 1e-2
DELTAS_H, GAINS_H, FLOORS_H = (0.0015, 0.0012, 0.0010, 0.0008, 0.0006, 0.0004), (1.0, 2.0), (0.01, 0.0)
H_SO = 2e-3                                                                # the second order's Hessian step on M
GAMMA_D, GAINS_D, S_D, SEEDS_D = 1e-2, (1.0, 2.0), 1e6, X.SEEDS_H
T_D = X.T_H
GAMMA_F, EPS_F, SEEDS_F, S_F, TA_F = 3e-3, 1e-4, tuple(range(400, 416)), 4e5, 1e5
MODES_F = (('none', 0.0), ('plug-in', 2.0), ('masses', 2.0), ('tracked', 2.0))
ESCAPE = 10.0                                                              # |α_k| beyond it: up the positions' plateau (|a_k|/a_max > 0.986)
GAMMA_G, EPS_G, GAINS_G, SEEDS_G, S_G = 1e-2, 3e-4, (0.0, 1.0, 2.0), tuple(range(200, 216)), 4e5   # exp_gauge_drift.py's Part G at ε = 3·10⁻⁴
REC = GD.REC


def instance(gamma, kF, floor=0.01, eps=0.0):
    """X.instance's flow with the correction of gain κ_F (and the entropic leak ε), and F₁ = the perturbed flow minus the plain one."""
    net, base = E.setup()
    loops = [MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, eps=eps, fixman=kF, fixman_tau=TAU, fixman_floor=floor),
             MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)]

    def beta_of(z):
        e = np.exp(np.r_[z, 0.0]); return 3 * e / e.sum()

    def flow(loop):
        def F(x):
            T = x[:9].reshape(3, 3); b = beta_of(x[9:]); dT, db = loop.drift(T, b); dl = db / b
            return np.concatenate([dT.ravel(), dl[:2] - dl[2]])
        return F
    Fk, F0 = flow(loops[0]), flow(loops[1])
    _, chart, coords, S = X.instance(gamma)
    return F0, (lambda x: Fk(x) - F0(x)), chart, coords, S


def drift(w, gamma, kF, floor=0.01, eps=0.0):
    F0, F1, chart, coords, S = instance(gamma, kF, floor, eps)
    return GDR.limit_drift(F0, chart, coords, np.asarray(w, float), S, h=X.H_ID, h2=X.H2_ID, F1=(F1 if (kF or eps) else None))


T3 = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]])


def second_order(gamma, eps, kF, start=None, hH=H_SO):
    """The stationary mean of the limit SDE dW = f(W)ds + D^{1/2}dB on M to second order in its fluctuations, f the drift with the leak and the
    correction (`drift(w, γ, κ_F, eps=ε)`): W₀ = the zero of f, J = ∂f(W₀), C from JC + CJᵀ + D(W₀) = 0, E[W] = W₀ − J⁻¹·½∂²f[C]. Returns the
    shifts from 1/3 of W₀ and of E[W], the fluctuation term split by f's parts (b, the leak, the correction), and the std √diag C."""
    import warnings
    from scipy.optimize import fsolve
    from scipy.linalg import solve_continuous_lyapunov
    w = np.ones(2) / 3
    fpart = {'b': lambda v: drift(v, gamma, 0.0)[0],
             'leak': lambda v: drift(v, gamma, 0.0, eps=eps)[0] - drift(v, gamma, 0.0)[0],
             'correction': lambda v: (drift(v, gamma, kF)[0] - drift(v, gamma, 0.0)[0]) if kF else np.zeros(2)}
    f = lambda v: sum(fp(v) for fp in fpart.values())
    with warnings.catch_warnings():                                         # xtol is below the finite differences' noise at κ_F = 1;
        warnings.simplefilter('ignore', RuntimeWarning)                     # the residual is returned instead
        W0 = fsolve(lambda v: f(v) * 1e7, w if start is None else w + np.asarray(start)[:2], xtol=1e-12)
    J = np.array([(f(W0 + 1e-4 * e) - f(W0 - 1e-4 * e)) / 2e-4 for e in np.eye(2)]).T
    C = solve_continuous_lyapunov(J, -drift(W0, gamma, kF, eps=eps)[1])
    e1, e2 = np.eye(2)

    def half_d2(fun):
        d2 = lambda u, v: (fun(W0 + hH * (u + v)) - fun(W0 + hH * (u - v)) - fun(W0 - hH * (u - v)) + fun(W0 - hH * (u + v))) / (4 * hH * hH)
        return 0.5 * (C[0, 0] * d2(e1, e1) + C[1, 1] * d2(e2, e2) + 2 * C[0, 1] * d2(e1, e2))
    parts = {k: -(T3 @ np.linalg.solve(J, half_d2(fp))) for k, fp in fpart.items()}
    W0_3 = T3 @ W0 + np.r_[0.0, 0.0, 1.0]
    return dict(W0=W0_3 - 1 / 3, parts=parts, mean=W0_3 - 1 / 3 + sum(parts.values()), std=np.sqrt(np.diag(T3 @ C @ T3.T)),
                resid=float(np.abs(f(W0)).max()), scale=float(np.abs(f(w)).max()))


def mean_path(gamma, kF, S, start=None, n_eval=51):
    """The deterministic mean path dW/ds = b_c(W) (Katzenberger's drift with the correction) from 1/3: its average over [S/2, S]."""
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, w: drift(w, gamma, kF)[0], (0.0, S), np.ones(2) / 3 if start is None else start, t_eval=np.linspace(S / 2, S, n_eval),
                    rtol=1e-6, atol=1e-10)
    W = np.c_[sol.y.T, 1 - sol.y.sum(axis=0)]
    return W.mean(axis=0)


def path_job(args):
    return mean_path(*args)


def w_sharp():
    """argmax of ℓ = 2Σ log(W_k − W_min,k²/W_k) on the simplex (Lagrange: ℓ's partial derivatives equal)."""
    from scipy.optimize import brentq
    g = lambda W, k: (1 + W_MIN[k] ** 2 / W ** 2) / (W - W_MIN[k] ** 2 / W)        # ½∂ℓ/∂W_k, decreasing in W

    def W_of(mu, k):
        return brentq(lambda W: g(W, k) - mu, W_MIN[k] * (1 + 1e-12), 50.0)
    mu = brentq(lambda mu: sum(W_of(mu, k) for k in range(3)) - 1, 0.05, 1e6)
    return np.array([W_of(mu, k) for k in range(3)])


def a_job(args):
    w, gamma = args
    law = X.log_det_law(np.r_[w, 1 - sum(w)])
    b = drift(w, gamma, 0.0)[0]; bc = drift(w, gamma, 1.0)[0]
    return dict(w=w, gamma=gamma, law=law, b=np.r_[b, -b.sum()] / (gamma * TAU), bc=np.r_[bc, -bc.sum()] / (gamma * TAU))


def conditioning():
    """W♯ maximises det ∇²_ΘJ along M (D-optimal); the condition number along M and its minimiser (Nelder–Mead from W♯), the smallest
    eigenvalue's maximiser, for comparison (referee round 28: W♯ is not the best-conditioned exact fit)."""
    from scipy.optimize import minimize
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    h = 1e-4

    def hess(w):
        W = np.r_[w, 1 - w.sum()]; beta = 3 * W
        Th = E.exact_rep(net, W).ravel()
        J = lambda t: loop.J(t.reshape(3, 3), beta)
        H = np.zeros((9, 9)); I = np.eye(9) * h
        for i in range(9):
            for j in range(i, 9):
                H[i, j] = H[j, i] = (J(Th + I[i] + I[j]) - J(Th + I[i] - I[j]) - J(Th - I[i] + I[j]) + J(Th - I[i] - I[j])) / (4 * h * h)
        return np.linalg.eigvalsh(H)
    ws = w_sharp()
    ev = hess(ws[:2])
    cond = lambda w: (lambda e: e[-1] / e[0])(hess(w)) if np.all(np.r_[w, 1 - w.sum()] > W_MIN + 1e-3) else 1e9
    lmin = lambda w: -hess(w)[0] if np.all(np.r_[w, 1 - w.sum()] > W_MIN + 1e-3) else 1e9
    rc = minimize(cond, ws[:2], method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-6))
    rl = minimize(lmin, ws[:2], method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-10))
    return dict(cond_sharp=ev[-1] / ev[0], lmin_sharp=ev[0], w_cond=np.r_[rc.x, 1 - rc.x.sum()], cond_min=rc.fun,
                w_lmin=np.r_[rl.x, 1 - rl.x.sum()], lmin_max=-rl.fun)


def b_job(args):
    gamma, kF = args
    from scipy.optimize import fsolve
    ws = w_sharp()

    def f(w):
        if not np.all(np.r_[w, 1 - np.sum(w)] > W_MIN + 1e-3):
            return np.full(2, 1e3)
        return drift(w, gamma, kF)[0] / (gamma * TAU)
    z, info, ier, msg = fsolve(f, ws[:2], full_output=True)
    h = 1e-4
    Jb = np.array([(f(z + h * e) - f(z - h * e)) / (2 * h) for e in np.eye(2)]).T * gamma * TAU
    std = None
    if np.all(np.linalg.eigvals(Jb).real < 0):                                # the linearised limit SDE's stationary spread there
        from scipy.linalg import solve_continuous_lyapunov
        Cw = solve_continuous_lyapunov(Jb, -drift(z, gamma, kF)[1])
        T3 = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]]); std = np.sqrt(np.diag(T3 @ Cw @ T3.T))
    return dict(gamma=gamma, kF=kF, zero=np.r_[z, 1 - z.sum()], ier=int(ier), eig=np.linalg.eigvals(Jb), resid=np.abs(f(z)).max(), std=std)


def c_job(args):
    kF, dl = args[:2]
    floor = args[2] if len(args) > 2 else 0.01
    W1 = W_MIN[0] + dl; w = np.array([W1, (1 - W1) / 2]); hh = 0.05 * dl
    b, D, _ = drift(w, GAMMA_C, kF, floor)
    Dp = drift(w + np.array([hh, -hh / 2]), GAMMA_C, kF, floor)[1]; Dm = drift(w - np.array([hh, -hh / 2]), GAMMA_C, kF, floor)[1]
    dlogD = (np.log(Dp[0, 0]) - np.log(Dm[0, 0])) / (2 * hh)
    return dict(kF=kF, delta=dl, floor=floor, expo=dl * 2 * b[0] / D[0, 0] - dl * dlogD)


def delta_cap(k=0, floor=0.01):
    """The distance from the k-th class's saturation within which the floor binds: W(1 − W_min²/W²) = floor·W."""
    return W_MIN[k] * (1 / np.sqrt(1 - floor) - 1)


@cached
def series(args):
    """GD.series with the correction in the given form ('none', 'plug-in', 'masses' with Â fixed at the initial exact fit, 'tracked' with Â
    the running average of W_k a_k over T_A): records (s, α_1..3, W_1..3, the fit error) every REC steps, and Â_1..3 for the masses forms."""
    gamma, eps, kF, seed, S = args[:5]
    form = args[5] if len(args) > 5 else 'plug-in'
    t0 = time.time()
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3).copy()
    A0 = base.features(Th)[0] / 3 if form in ('masses', 'tracked') else None          # Â = W_k a_k at the start: the bumps' amplitudes
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, eps=eps, fixman=(0.0 if form == 'none' else kF), fixman_tau=TAU,
                                     fixman_A=A0, fixman_TA=(TA_F if form == 'tracked' else None))
    rec = []
    extra = (lambda: loop.fixman_A) if A0 is not None else (lambda: np.empty(0))
    loop.run_dithered(Th, np.ones(3), S, DT, TAU, np.random.default_rng(seed), rec_every=REC,
                      callback=lambda s_, T, b: rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b), extra()]))
    print(f"[{time.time() - t0:.0f}s] gamma {gamma:g} eps {eps:g} {form} kF {kF:g} seed {seed}", flush=True)
    return dict(gamma=gamma, eps=eps, kF=kF, form=form, seed=seed, rec=np.array(rec))


def summarise_long(R):
    s = R[0, :, 0]; W = R[:, :, 4:7]; al = R[:, :, 1:4]
    q = np.array_split(np.arange(R.shape[1]), 4)
    out = dict(times=[], quarters=np.array([np.abs(al[:, ix]).max(axis=1).mean(axis=0) for ix in q]),
               big=np.array([np.mean(np.abs(al[:, ix]) > 10, axis=(0, 1)) for ix in q]), below=np.mean(W[:, q[-1]] < W_MIN, axis=(0, 1)),
               half=W[:, R.shape[1] // 2:].mean(axis=1), fit=R[:, R.shape[1] // 2:, 7].mean(), fit2=(R[:, R.shape[1] // 2:, 7] ** 2).mean(),
               almax=np.abs(al).max(axis=(0, 1)))
    for t in T_D:
        i = int(np.argmin(np.abs(s - t)))
        out['times'].append(dict(s=s[i], mean=W[:, i].mean(axis=0), se=W[:, i].std(axis=0, ddof=1) / np.sqrt(len(W)), std=W[:, i].std(axis=0, ddof=1),
                                 sat=np.abs(np.tanh(al[:, i] / A_MAX)).max(axis=0)))
    return out


def kappa(W, eps, gamma):
    A = E.A_BUMP; K = 3; m = 3; Ga = E.G[0]; W = np.asarray(W, float)
    a = A / W; ap = 1 - (a / A_MAX) ** 2
    T = np.zeros((2 * K, K - 1))
    for j in range(K - 1):
        dW = np.zeros(K); dW[j] = 1.0; dW[-1] = -1.0
        T[:K, j] = -(A / (W ** 2 * ap)) * dW; T[K:, j] = dW
    g = np.diag(np.r_[np.full(K, 1 / Ga), (m / gamma) / W])
    return np.linalg.solve(T.T @ g @ T, T.T @ np.diag(np.r_[np.zeros(K), eps * m / (gamma * W)]) @ T)


def part_F(rec):
    """rec: {(form, κ_F): records (seeds, n, 8 or 11)}. Per mode: the seeds with an escape up the positions' plateau (some |α_k| > ESCAPE)
    and the first such record, the masses' range, the second half's shift of the masses from 1/3 over all seeds and over the seeds without
    an escape (± s.e. over the seeds), the pooled within-seed std, the fit error, and the predicted shift κ⁻¹b."""
    w = np.ones(3) / 3
    pred, so = {}, {}
    out = []
    for form, kF in MODES_F:
        R = np.asarray(rec[(form, kF)], float); h = R.shape[1] // 2
        k_eff = 0.0 if form == 'none' else kF
        if k_eff not in pred:
            c = np.linalg.solve(kappa(w, EPS_F, GAMMA_F), drift(w[:2], GAMMA_F, k_eff)[0]); pred[k_eff] = np.r_[c, -c.sum()]
            so[k_eff] = second_order(GAMMA_F, EPS_F, k_eff, start=pred[k_eff])
        esc = np.abs(R[:, :, 1:4]).max(axis=(1, 2)) > ESCAPE
        first = [float(R[i, np.argmax(np.abs(R[i, :, 1:4]).max(axis=1) > ESCAPE), 0]) for i in np.flatnonzero(esc)]
        means = R[:, h:, 4:7].mean(axis=1); ok = ~esc
        se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M)) if len(M) > 1 else np.full(3, np.nan)
        out.append(dict(form=form, kF=kF, seeds_esc=[SEEDS_F[i] for i in np.flatnonzero(esc)], first=first,
                        Wmin=R[:, :, 4:7].min(axis=(0, 1)), Wmax=R[:, :, 4:7].max(axis=(0, 1)), almax=np.abs(R[:, :, 1:4]).max(),
                        shift=means.mean(axis=0) - 1 / 3, se=se(means), shift_ok=means[ok].mean(axis=0) - 1 / 3, se_ok=se(means[ok]),
                        std_ok=R[ok][:, h:, 4:7].std(axis=1).mean(axis=0), fit=R[:, h:, 7].mean(), fit_ok=R[ok][:, h:, 7].mean(),
                        fit_esc=(R[esc][:, h:, 7].mean() if esc.any() else np.nan), pred=pred[k_eff], so=so[k_eff],
                        A_end=(R[:, -1, 8:11] if R.shape[2] > 8 else None)))
    return out


def long_tables(D, title):
    fmt = E.fmt
    ws = w_sharp()
    L_ = [title, "", "| κ_F | s | mean of W_1, W_2, W_3 (± s.e.) | std | the largest |a|/a_max |", "|---|---|---|---|---|"]
    for r in D:
        for t in r['sum']['times']:
            L_.append(f"| {r['kF']:g} | {t['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(t['mean'], t['se']))} | "
                      f"{', '.join(f'{v:.4f}' for v in t['std'])} | {', '.join(f'{v:.4f}' for v in t['sat'])} |")
    L_ += ["", "| κ_F | mean largest |α_1|, |α_2|, |α_3| by quarters | fraction |α_k| > 10 in the last quarter | beyond the saturation (last quarter) | "
           "the second half's mean (± s.e. over the seeds) | W♯ | its distance from W♯ in s.e. | the mean path's second-half average | "
           "the distance from it in s.e. | the fit error (second half; rms) | the largest |α_k| over all records |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in D:
        sm = r['sum']; hm = sm['half']; se = hm.std(axis=0, ddof=1) / np.sqrt(len(hm)); pth = r.get('path')
        L_.append(f"| {r['kF']:g} | {'; '.join(', '.join(f'{v:.1f}' for v in q) for q in sm['quarters'])} | {', '.join(f'{v:.2f}' for v in sm['big'][-1])} | "
                  f"{', '.join(f'{v:.2f}' for v in sm['below'])} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(hm.mean(axis=0), se))} | "
                  f"{', '.join(f'{v:.4f}' for v in ws)} | {', '.join(f'{v:+.2f}' for v in (hm.mean(axis=0) - ws) / se)} | "
                  f"{'—' if pth is None else ', '.join(f'{v:.4f}' for v in pth)} | "
                  f"{'—' if pth is None else ', '.join(f'{v:+.2f}' for v in (hm.mean(axis=0) - pth) / se)} | "
                  f"{fmt(sm['fit'])}; {'—' if 'fit2' not in sm else fmt(float(np.sqrt(sm['fit2'])))} | "
                  f"{'—' if 'almax' not in sm else ', '.join(f'{v:.1f}' for v in sm['almax'])} |")
    return L_


def equipartition(D, F2):
    """At fixed masses the positions' law is Gibbs: each coordinate holding fit energy stores τ/2 of J = (m/2)⟨r²⟩, so ⟨e²⟩ = pτ/(m⟨F²⟩)
    (e² the squared relative fit error); read backwards, p_eff = m⟨F²⟩⟨e²⟩/τ and, with p = 9, τ̂ = m⟨F²⟩⟨e²⟩/9 (the rms, not the mean of e)."""
    return [(r['kF'], 3 * F2 * r['sum']['fit2'] / TAU, 3 * F2 * r['sum']['fit2'] / 9) for r in D if 'fit2' in r['sum']]


def part_G(rec):
    """rec: {κ_F: records (seeds, n, 8 or 11)} (κ_F = 0 uncorrected, else the masses form, the same seeds). The second half's shift of the
    masses from 1/3 (± s.e. over the seeds) against κ⁻¹b_c(w), and the paired differences from κ_F = 0 against κ⁻¹(b_c − b)."""
    from scipy.optimize import fsolve
    w = np.ones(3) / 3
    K = kappa(w, EPS_G, GAMMA_G)
    pred, nl = {}, {}
    for k in rec:
        c = np.linalg.solve(K, drift(w[:2], GAMMA_G, k)[0]); pred[k] = np.r_[c, -c.sum()]
        f = lambda c_: (drift(w[:2] + c_, GAMMA_G, k)[0] - kappa(np.r_[w[:2] + c_, 1 - (w[:2] + c_).sum()], EPS_G, GAMMA_G) @ c_) * 1e7
        cn = fsolve(f, c); nl[k] = np.r_[cn, -cn.sum()]                     # the nonlinear fixed point b_c(W) = κ(W)(W − w)
    so = {k: second_order(GAMMA_G, EPS_G, k, start=nl[k]) for k in rec}        # the second order (the fluctuation term)
    R0 = np.asarray(rec[0.0], float); h = R0.shape[1] // 2; m0 = R0[:, h:, 4:7].mean(axis=1)
    out = []
    for k in sorted(rec):
        R = np.asarray(rec[k], float); means = R[:, h:, 4:7].mean(axis=1); d = means - m0
        esc = np.abs(R[:, :, 1:4]).max(axis=(1, 2)) > ESCAPE
        out.append(dict(kF=k, shift=means.mean(axis=0) - 1 / 3, se=means.std(axis=0, ddof=1) / np.sqrt(len(means)), pred=pred[k], nl=nl[k],
                        diff=d.mean(axis=0), dse=d.std(axis=0, ddof=1) / np.sqrt(len(d)), dpred=pred[k] - pred[0.0], dnl=nl[k] - nl[0.0],
                        so=so[k], dso=so[k]['mean'] - so[0.0]['mean'],
                        esc=[SEEDS_G[i] for i in np.flatnonzero(esc)], std=R[:, h:, 4:7].std(axis=1).mean(axis=0), fit=R[:, h:, 7].mean()))
    return out


def part_G37(zg):
    """exp_gauge_drift.py's Part G (Proposition D‴.37 (G), no correction, ε = 10⁻⁴, 3·10⁻⁴, 10⁻³): the stored shifts against the second order."""
    return [dict(eps=g['eps'], shift=g['shift'], se=g['se'], lin=g['lin'], nl=g['nl'], so=second_order(GAMMA_G, g['eps'], 0.0, start=g['nl']))
            for g in zg]


def report(A=None, B=None, C=None, D=None, Ee=None, Ff=None, Gg=None, Hh=None, G37=None, Bc=None):
    fmt = E.fmt
    w3 = lambda w: ', '.join(f'{v:.3f}' for v in np.r_[w, 1 - sum(w)])
    L_ = ["# Proposition D‴.38: a Fixman correction for the replicator under dither (`exp_gauge_fixman.py`)", "",
          f"Proposition D‴.36's instance (three bumps, one particle per class at an exact representation, m = 3, the positions dithered at τ = {TAU:g}); "
          "the correction adds κ_F(τ/m)(1 + W_min,k²/W_k²)/(W_k − W_min,k²/W_k) to the replicator's fitness, with W_min,k = W_k|a_k|/a_max at the current "
          "amplitudes (the plug-in form, κ_F(τ/m)(1 + (a_k/a_max)²)/(W_k a′_k)) or W_min,k = |Â_k|/a_max from the classes' contributions (the masses form), "
          "equal on the manifold of exact fits (Parts A–C hold for both); drifts per unit slow time, W_k = β_k/3.", ""]
    if A is not None:
        L_ += ["## A. The cancellation (κ_F = 1)", "",
               "| W | γ | uncorrected drift / γτ | its distance from the log-determinant law | corrected drift / γτ | the corrected over the uncorrected (largest components) |",
               "|---|---|---|---|---|---|"]
        for r in A:
            L_.append(f"| {w3(r['w'])} | {r['gamma']:g} | {', '.join(f'{v:+.4f}' for v in r['b'])} | {np.abs(r['b'] - r['law']).max():.2e} | "
                      f"{', '.join(f'{v:+.5f}' for v in r['bc'])} | {np.abs(r['bc']).max() / np.abs(r['b']).max():.2e} |")
        L_.append("")
    if B is not None:
        ws = w_sharp()
        L_ += ["## B. The gain: the corrected field's zero", "", f"W♯ = argmax ℓ on the simplex: {', '.join(f'{v:.5f}' for v in ws)} "
               f"(|a|/a_max there {', '.join(f'{v:.3f}' for v in W_MIN / ws)}): the D-optimal exact fit (the largest det ∇²_ΘJ).", ""]
        if Bc is not None:
            L_ += [f"Not the best-conditioned: the condition number of ∇²_ΘJ along M is {Bc['cond_sharp']:.1f} at W♯ and {Bc['cond_min']:.1f} at its minimiser "
                   f"({', '.join(f'{v:.4f}' for v in Bc['w_cond'])}); the smallest eigenvalue is {Bc['lmin_sharp']:.6f} at W♯ and {Bc['lmin_max']:.6f} at its "
                   f"maximiser ({', '.join(f'{v:.4f}' for v in Bc['w_lmin'])}) — both nearly flat along M.", ""]
        L_ += [
               "| κ_F | γ | the zero (Newton from W♯) | its distance from W♯ | the Jacobian's eigenvalues (per unit slow time) | the residual | "
               "the stationary std of the limit SDE linearised there |", "|---|---|---|---|---|---|---|"]
        for r in B:
            sd = r.get('std')
            L_.append(f"| {r['kF']:g} | {r['gamma']:g} | {', '.join(f'{v:.5f}' for v in r['zero'])}{'' if r['ier'] == 1 else ' (no convergence)'} | "
                      f"{np.abs(r['zero'] - ws).max():.1e} | {', '.join(fmt(v.real) for v in r['eig'])} | {r['resid']:.1e} | "
                      f"{'—' if sd is None else ', '.join(f'{v:.4f}' for v in sd)} |")
        L_.append("")
    if C is not None:
        L_ += [f"## C. The saturation: the density's exponent near the first class's saturation (γ = {GAMMA_C:g}, W_2 = W_3)", "",
               "| κ_F | d log p / d log δ at δ = " + ', '.join(f'{d:g}' for d in DELTAS_C) + " | predicted |", "|---|---|---|"]
        for kF in GAINS_C:
            ex = [r['expo'] for r in C if r['kF'] == kF]
            L_.append(f"| {kF:g} | {', '.join(f'{v:+.3f}' for v in ex)} | {kF - 1:+g} |")
        L_.append("")
    if D is not None:
        L_ += long_tables(D, f"## D. The full loop over {S_D:g} slow units, the plug-in form (γ = {GAMMA_D:g}, {len(SEEDS_D)} seeds; κ_F = 0 from "
                             "exp_gauge_drift.py's Part H, the same seeds)")
        net, base = E.setup()
        F2 = float(np.mean(base.F ** 2))
        eq = equipartition(D, F2)
        L_ += ["", "Equipartition: at fixed masses the positions' law is Gibbs, each position coordinate that holds fit energy stores τ/2 of "
               "J = (m/2)⟨r²⟩, so the rms fit error is √(pτ/(m⟨F²⟩)) with p such coordinates — " + ', '.join(f"{p_}: {np.sqrt(p_ * TAU / (3 * F2)):.4f}" for p_ in (9, 8, 7))
               + f" (⟨F²⟩ = {F2:.4f}); a saturated amplitude (a′ → 0) holds none. Read backwards from the records' ⟨e²⟩ (the second half): "
               "p_eff = m⟨F²⟩⟨e²⟩/τ = " + ', '.join(f"{pe:.2f} (κ_F = {k:g})" for k, pe, _ in eq) + "; the dither's level with p = 9, "
               "τ̂ = m⟨F²⟩⟨e²⟩/9 = " + ', '.join(f"{th:.3g} (κ_F = {k:g})" for k, _, th in eq) + f" against τ = {TAU:g}. (The mean of e is lower than "
               "the rms by Jensen's inequality, E[χ_p]²/p = 0.946 for p = 9.)"]
        L_.append("")
    if Ee is not None:
        L_ += long_tables(Ee, f"## E. The full loop over {S_D:g} slow units, the masses form with Â fixed at the initial exact fit (γ = {GAMMA_D:g}, "
                              f"the seeds of D)")
        net, base = E.setup()
        eq = equipartition(Ee, float(np.mean(base.F ** 2)))
        if eq:
            L_ += ["", "Equipartition as in D: p_eff = " + ', '.join(f"{pe:.2f} (κ_F = {k:g})" for k, pe, _ in eq) + "; τ̂ = "
                   + ', '.join(f"{th:.3g} (κ_F = {k:g})" for k, _, th in eq) + "."]
        L_.append("")
    if Ff is not None:
        L_ += [f"## F. With the entropic leak (γ = {GAMMA_F:g}, ε = {EPS_F:g}, {len(SEEDS_F)} seeds, {S_F:g} slow units): the forms against the "
               f"positions' plateau", "",
               f"An escape: some |α_k| > {ESCAPE:g} (|a_k|/a_max > {np.tanh(ESCAPE / A_MAX):.3f}) in a record; the tracked form's T_A = {TA_F:g}.", "",
               "| correction | κ_F | seeds with an escape (first record beyond, slow units) | the largest |α| | the masses' range (min; max over the records) | "
               "the fit error: all seeds; without an escape; with one |", "|---|---|---|---|---|---|"]
        for r in Ff:
            esc = ', '.join(f"{s} ({t:g})" for s, t in zip(r['seeds_esc'], r['first'])) or 'none'
            L_.append(f"| {r['form']} | {r['kF']:g} | {esc} | {r['almax']:.1f} | {', '.join(f'{v:.3f}' for v in r['Wmin'])}; "
                      f"{', '.join(f'{v:.3f}' for v in r['Wmax'])} | {fmt(r['fit'])}; {fmt(r['fit_ok'])}; {fmt(r['fit_esc'])} |")
        L_ += ["", "| correction | κ_F | the masses' mean − 1/3, all seeds (± s.e.) | the seeds without an escape (± s.e.) | predicted κ⁻¹b | "
               "its distance from the prediction in s.e. (without an escape) | the second order | its distance in s.e. | the within-seed std "
               "(without an escape; the second order's) |", "|---|---|---|---|---|---|---|---|---|"]
        for r in Ff:
            so = r.get('so')
            L_.append(f"| {r['form']} | {r['kF']:g} | {', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift'], r['se']))} | "
                      f"{', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift_ok'], r['se_ok']))} | {', '.join(f'{v:+.5f}' for v in r['pred'])} | "
                      f"{', '.join(f'{v:+.2f}' for v in (r['shift_ok'] - r['pred']) / r['se_ok'])} | "
                      f"{'—' if so is None else ', '.join(f'{v:+.5f}' for v in so['mean'])} | "
                      f"{'—' if so is None else ', '.join(f'{v:+.2f}' for v in (r['shift_ok'] - so['mean']) / r['se_ok'])} | "
                      f"{', '.join(f'{v:.4f}' for v in r['std_ok'])}{'' if so is None else '; ' + ', '.join(f'{v:.4f}' for v in so['std'])} |")
        tr = [r for r in Ff if r['form'] == 'tracked']
        if tr:
            Ae = tr[0]['A_end']; esc = np.isin(SEEDS_F, tr[0]['seeds_esc'])
            L_ += ["", f"The tracked contributions Â at {S_F:g} (the bumps' amplitudes {', '.join(f'{v:+.3f}' for v in E.A_BUMP)}): the seeds without an "
                   f"escape {', '.join(f'{v:+.3f}' for v in Ae[~esc].mean(axis=0))} (mean)"
                   + (f"; with one {'; '.join(', '.join(f'{v:+.3f}' for v in a) for a in Ae[esc])}" if esc.any() else '') + "."]
        L_.append("")
    if Gg is not None:
        v3 = lambda v, p=5: ', '.join(f'{x:+.{p}f}' for x in v)
        L_ += [f"## G. With the entropic leak at γ = {GAMMA_G:g}, ε = {EPS_G:g} (the seeds of exp_gauge_drift.py's Part G, {S_G:g} slow units, the second "
               "half): the masses form, and its difference from κ_F = 0 on the same seeds", "",
               "The second order: E[W] = W₀ − J⁻¹·½∂²f[C] for the limit SDE on M with the leak and the correction (f its drift, W₀ its zero, J its "
               f"Jacobian, JC + CJᵀ + D = 0; Hessian step {H_SO:g}); the fluctuation term split by f's parts; |f(W₀)| ≤ "
               f"{max(r['so']['resid'] for r in Gg):.1e} against |f(1/3)| ≥ {min(r['so']['scale'] for r in Gg):.1e}.", "",
               "| κ_F | the masses' mean − 1/3 (± s.e.) | first order κ⁻¹b_c (the nonlinear fixed point W₀) | its distance in s.e. | "
               "the second order | its distance in s.e. | the fluctuation term: from b; the leak; the correction | the within-seed std; the second order's | "
               "seeds with an escape | the fit error |", "|---|---|---|---|---|---|---|---|---|---|"]
        for r in Gg:
            so = r['so']
            L_.append(f"| {r['kF']:g} | {', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift'], r['se']))} | {v3(r['pred'])} ({v3(r['nl'])}) | "
                      f"{', '.join(f'{v:+.2f}' for v in (r['shift'] - r['pred']) / r['se'])} | {v3(so['mean'])} | "
                      f"{', '.join(f'{v:+.2f}' for v in (r['shift'] - so['mean']) / r['se'])} | "
                      f"{v3(so['parts']['b'])}; {v3(so['parts']['leak'])}; {v3(so['parts']['correction'])} | "
                      f"{', '.join(f'{v:.4f}' for v in r['std'])}; {', '.join(f'{v:.4f}' for v in so['std'])} | {', '.join(str(x) for x in r['esc']) or 'none'} | "
                      f"{fmt(r['fit'])} |")
        L_ += ["", "| κ_F | the difference from κ_F = 0, the same seeds (± s.e.) | first order κ⁻¹(b_c − b) (nonlinear) | its distance in s.e. (nonlinear) | "
               "the second order | its distance in s.e. |", "|---|---|---|---|---|---|"]
        for r in Gg:
            if r['kF'] == 0:
                continue
            L_.append(f"| {r['kF']:g} | {', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['diff'], r['dse']))} | {v3(r['dpred'])} ({v3(r['dnl'])}) | "
                      f"{', '.join(f'{v:+.2f}' for v in (r['diff'] - r['dpred']) / r['dse'])} ({', '.join(f'{v:+.2f}' for v in (r['diff'] - r['dnl']) / r['dse'])}) | "
                      f"{v3(r['dso'])} | {', '.join(f'{v:+.2f}' for v in (r['diff'] - r['dso']) / r['dse'])} |")
        zg = E.load(os.path.join(OUT, 'gauge_drift_raw.npz')).get('G')
        old = [g for g in (zg if zg is not None else []) if g['eps'] == EPS_G]
        if old:
            d = np.abs(Gg[0]['shift'] - old[0]['shift']).max()
            L_ += ["", f"κ_F = 0 repeats exp_gauge_drift.py's Part G at ε = {EPS_G:g} (the same seeds and steps): its stored shift "
                   f"{', '.join(f'{v:+.5f}' for v in old[0]['shift'])}, the largest difference {d:.1e}."]
        if G37 is not None:
            L_ += ["", "exp_gauge_drift.py's Part G (no correction, 16 seeds, 4·10⁵ slow units) at its three ε against the second order:", "",
                   "| ε | the stored shift (± s.e.) | κ⁻¹b (nonlinear) | its distance in s.e. | the second order | its distance in s.e. | "
                   "the fluctuation term from the leak |", "|---|---|---|---|---|---|---|"]
            for r in G37:
                L_.append(f"| {r['eps']:g} | {', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift'], r['se']))} | {v3(r['lin'])} ({v3(r['nl'])}) | "
                          f"{', '.join(f'{v:+.2f}' for v in (r['shift'] - r['lin']) / r['se'])} | {v3(r['so']['mean'])} | "
                          f"{', '.join(f'{v:+.2f}' for v in (r['shift'] - r['so']['mean']) / r['se'])} | {v3(r['so']['parts']['leak'])} |")
        L_.append("")
    if Hh is not None:
        dc = delta_cap(0)
        L_ += [f"## H. The implemented floor: the exponent of C below δ_cap = W_min,1(1/√0.99 − 1) = {dc:.5f} (γ = {GAMMA_C:g}, W_2 = W_3)", "",
               "The correction's denominator is floored at 0.01W_k: within δ_cap of the saturation the correction stops growing, its contribution to the "
               "exponent falls from κ_F to κ_F(1 − W_min²/W²)/0.01 ≈ κ_F·200δ/W_min,1, and the τ → 0 density goes back to δ⁻¹ (log-divergent mass).", "",
               "| κ_F | floor | d log p / d log δ at δ = " + ', '.join(f'{d:g}' for d in DELTAS_H) + " | predicted −1 + κ_F·min(1, (1 − W_min²/W²)/floor) | "
               "the unfloored exponent minus κ_F(1 − min(1, (1 − W_min²/W²)/floor)) |", "|---|---|---|---|---|"]
        for kF in GAINS_H:
            ex0 = [r['expo'] for r in Hh if r['kF'] == kF and r['floor'] == 0.0]
            for fl in FLOORS_H:
                ex = [r['expo'] for r in Hh if r['kF'] == kF and r['floor'] == fl]
                cut = [min(1.0, (1 - (W_MIN[0] / (W_MIN[0] + d)) ** 2) / fl) if fl > 0 else 1.0 for d in DELTAS_H]
                pr = [-1 + kF * c for c in cut]
                pr0 = [e0 - kF * (1 - c) for e0, c in zip(ex0, cut)]
                L_.append(f"| {kF:g} | {fl:g} | {', '.join(f'{v:+.3f}' for v in ex)} | {', '.join(f'{v:+.3f}' for v in pr)} | "
                          f"{'—' if fl == 0 else ', '.join(f'{v:+.3f}' for v in pr0)} |")
        dev = [[r['expo'] for r in Hh if r['kF'] == kF and r['floor'] == 0.0][-1] - (kF - 1) for kF in GAINS_H]
        L_ += ["", f"Without the floor the exponent departs from κ_F − 1 at the smallest δ ({DELTAS_H[-1]:g}) by "
               + ', '.join(f'{v:+.3f} (κ_F = {k:g})' for v, k in zip(dev, GAINS_H)) + " — the layer's own departure from its asymptotic law at "
               f"γ = {GAMMA_C:g}, the same for both gains; the floor's effect is the stated one on top of it (the last column).", ""]
    L_.append("Numbers are produced by `experiments/exp_gauge_fixman.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_fixman_raw.npz')
    md = os.path.join(OUT, 'gauge_fixman.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEFGH'], G37=get('G37'), Bc=get('Bc')))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFGH')
    re_an = '--reanalyse' in sys.argv                                         # D–G from the stored records
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'A' in parts:
            h['A'] = pool.map_async(a_job, [(w, g) for w in W_A for g in GAMMAS_A], chunksize=1)
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(g, k) for k in GAINS_B for g in GAMMAS_B], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(k, d) for k in GAINS_C for d in DELTAS_C], chunksize=1)
        if 'H' in parts:
            h['H'] = pool.map_async(c_job, [(k, d, fl) for k in GAINS_H for fl in FLOORS_H for d in DELTAS_H], chunksize=1)
        if 'D' in parts and not re_an:
            h['D'] = pool.map_async(series, [(GAMMA_D, 0.0, k, s, S_D) for k in GAINS_D for s in SEEDS_D], chunksize=1)
        if 'E' in parts and not re_an:
            h['E'] = pool.map_async(series, [(GAMMA_D, 0.0, k, s, S_D, 'masses') for k in GAINS_D for s in SEEDS_D], chunksize=1)
        if 'F' in parts and not re_an:
            h['F'] = pool.map_async(series, [(GAMMA_F, EPS_F, k, s, S_F, f) for f, k in MODES_F for s in SEEDS_F], chunksize=1)
        if 'G' in parts and not re_an:
            h['G'] = pool.map_async(series, [(GAMMA_G, EPS_G, k, s, S_G, 'none' if k == 0 else 'masses') for k in GAINS_G for s in SEEDS_G],
                                    chunksize=1)
        paths = pool.map_async(path_job, [(GAMMA_D, g, S_D) for g in GAINS_D], chunksize=1) if ({'D', 'E'} & parts) else None
        for k in 'ABCH':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
        if 'B' in parts:
            res['Bc'] = conditioning(); E.save(raw, **res); log('B (conditioning)')
        z = E.load(os.path.join(OUT, 'gauge_drift_raw.npz')) if ({'D', 'E'} & parts) else None
        if paths is not None:
            paths = dict(zip(GAINS_D, paths.get()))
        for k, key in (('D', 'D_rec'), ('E', 'E_rec')):
            if k in parts:
                if k in h:
                    runs = h[k].get()
                    res[key] = {g: np.array([r['rec'] for r in runs if r['kF'] == g]) for g in GAINS_D}
                res[k] = [dict(kF=0.0, sum=summarise_long(np.asarray(z['H_rec'], float)))] + \
                         [dict(kF=g, sum=summarise_long(np.asarray(res[key][g], float)), path=paths[g]) for g in GAINS_D]
                E.save(raw, **res); log(k)
        if 'F' in parts:
            if 'F' in h:
                runs = h['F'].get()
                res['F_rec'] = {(f, k): np.array([r['rec'] for r in runs if r['form'] == f and r['kF'] == k]) for f, k in MODES_F}
            res['F'] = part_F(res['F_rec']); E.save(raw, **res); log('F')
        if 'G' in parts:
            if 'G' in h:
                runs = h['G'].get()
                res['G_rec'] = {k: np.array([r['rec'] for r in runs if r['kF'] == k]) for k in GAINS_G}
            res['G'] = part_G(res['G_rec'])
            zg = E.load(os.path.join(OUT, 'gauge_drift_raw.npz')).get('G')
            res['G37'] = part_G37(zg) if zg is not None else None
            E.save(raw, **res); log('G')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEFGH'], G37=res.get('G37'), Bc=res.get('Bc')))


if __name__ == '__main__':
    main()
