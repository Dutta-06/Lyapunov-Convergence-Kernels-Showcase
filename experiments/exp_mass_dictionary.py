"""
Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise — the conserved mass inflates the fit; a
null kernel restores the colour's selection
====================================================================================================================================

The loop of Proposition D‴.34 at fixed positions on a dictionary (Φ, N × d): masses β (Σβ = m, or free), amplitudes a = a_max tanh(α/a_max),
w = β⊙a/m, label noise on the measured error; the selection along the exact fits is the gradient flow of tr Ψ(K) in the flat coordinates
(α, √β), K = (m/N)Φ diag(D)Φᵀ (− (γ/N)yyᵀ with the conserved mass), D_i = (β_i/m²)(Gβ_ia′_i² + γa_i²), plus, for the paradigm's update on
log β, the replicator drift of the filtered second moment of the log-masses' increment (`mass_dictionary`).

A. Identities at random exact fits: the closed forms against Katzenberger's drift of the scheme (`gauge_drift.limit_drift` on a chart of
   M) for the flat and the paradigm's updates, fresh and AR(1) draws at a finite step, conserved and free masses; the cancellation of the
   masses' own selection by the paradigm's Itô term; K's spectrum against the flow's normal rates; the excess-mass bound; the gauge-minimal
   relevance D* against its two asymptotic laws.
B. The paradigm's noisy scheme (64 paths) against the reduced flow: the paths' mean ‖w‖₁ and its largest deviation from the flow's w at
   nine times, for the conserved, free and null-kernel laws under fresh and AR(1) draws.
C. The selection's ends on random instances (N = 5, d = 12 Gaussian designs, 2- and 3-sparse truths of ‖w*‖₁ < a_max): the conserved, free
   and null-kernel laws under white and slow (AR(1), τ_c = 10³, 10⁴, 10⁵) noise — ‖w‖₁/a_max, the excess, the atoms used and alive, the
   truth recovered, the distance to basis pursuit; for white noise also the flat update's end and the minimiser of Σ_i|φ_i|²D*(|w_i|).
D. The same on Gaussian kernel dictionaries (16 centres, 8 samples on [0, 1]).
E. The conserved law's white end for a_max = 0.5, 1, 2, 4: ‖w‖₁/a_max — the inflation goes to the representability boundary, whatever it is.

Output: results/chaos/mass_dictionary.md and mass_dictionary_raw.npz (MD_WORKERS; `--from-raw` regenerates the md; `--only A,B,C,D`).
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

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mass_dictionary as MD
from src.core import colour_selection as CS
from src.core import gauge_drift as GDR
from src.core.runcache import cached
import exp_cohesion_leak as E

OUT = os.environ.get('MD_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('MD_WORKERS', 2))
N_I, D_I, M_I, G_I, GAM_I, AMAX = 5, 12, 12.0, 1.0, 0.1, 1.0
MODES = ('conserved', 'free', 'null kernel')
TAUS = (0.0, 1e3, 1e4, 1e5)
S_END = 3e8
KS_C, N_INST_C = (2, 3), 40
KS_D, N_INST_D, N_D, D_D, WIDTH_D = (2, 3), 20, 8, 16, 0.08
AMAXES_E, N_INST_E = (0.5, 1.0, 2.0, 4.0), 12
H_B, SIG_B, SIG_B_SLOW, P_B, SEED_B = 1.0, 0.1, 0.03, 64, 4             # slow draws: the fit follows the noise, so a smaller σ keeps it near M


def instance(k, i, kind='gauss'):
    """A design with unit columns and a k-sparse truth (signs ±1, magnitudes in [0.15, 0.3], so ‖w*‖₁ < a_max = 1) — Gaussian (N_I × D_I)
    or Gaussian kernels on [0, 1] (N_D samples, D_D centres) — and the start's random amplitudes."""
    rng = np.random.default_rng((1000 if kind == 'gauss' else 5000) + 100 * k + i)
    if kind == 'gauss':
        Phi = rng.standard_normal((N_I, D_I)); Phi /= np.linalg.norm(Phi, axis=0)
    else:
        Phi = CS.rbf_dictionary(np.linspace(0.0, 1.0, N_D), np.linspace(0.0, 1.0, D_D), WIDTH_D)
    d = Phi.shape[1]
    ws = np.zeros(d); idx = rng.choice(d, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.15, 0.3, k)
    return Phi, ws, 0.05 * rng.standard_normal(d + 1)


def law(Phi, ws, mode, a_max=AMAX):
    """The MassDictionary for a mode (the null kernel: a zero column appended) and its design."""
    P2 = np.c_[Phi, np.zeros(Phi.shape[0])] if mode == 'null kernel' else Phi
    return MD.MassDictionary(P2, Phi @ ws, M_I, a_max=a_max, G=G_I, gamma=GAM_I, balanced=(mode != 'free')), P2


def start(X, al0):
    """The noiseless flow's end from the amplitudes al0 and equal masses m/(atoms) (the loop's own law, in (α, log β)), on M."""
    from scipy.integrate import solve_ivp
    d = X.d
    th = X.theta(al0[:d], np.full(d, X.m / d))

    def F(z):
        x = np.r_[z[:d], np.exp(z[d:] / 2)]
        _, _, b, a, da, _ = X.parts(x)
        g = X.Phi.T @ (X.y - X.Phi @ (b * a / X.m)) / X.N
        u = a * g
        lam = b @ u / b.sum() if X.balanced else 0.0
        return np.r_[X.G * b * da * g, X.gamma * (u - lam)]
    z = solve_ivp(lambda t, z: F(z), (0.0, 2e5), np.r_[th[:d], 2 * np.log(th[d:])], method='LSODA', rtol=1e-11, atol=1e-13).y[:, -1]
    return X.retract(np.r_[z[:d], np.exp(z[d:] / 2)])


def filt_of(tau):
    return (lambda k: np.ones_like(k)) if tau == 0 else (lambda k, t=tau: 1.0 / (1.0 + t * k))


def summary(X, th, ws, wbp, d0):
    """The end's ‖w‖₁/a_max, excess/m, atoms used (|w| > 10⁻³), alive (β > 10⁻³), the truth recovered (10⁻⁴), distances, the null mass."""
    w = X.w(th)[:d0]; b = X.parts(th)[2]
    return dict(l1=float(np.abs(w).sum() / X.a_max), exc=float(X.excess(th) / X.m) if X.balanced else np.nan,
                used=int((np.abs(w) > 1e-3).sum()), alive=int((b[:d0] > 1e-3).sum()), ok=bool(np.abs(w - ws).max() < 1e-4 * np.abs(ws).max()),
                d_truth=float(np.abs(w - ws).max()), d_bp=float(np.abs(w - wbp).max()), b0=float(b[d0]) if len(b) > d0 else np.nan,
                amax=float(np.abs(X.parts(th)[3][:d0]).max()))


def dstar_minimiser(X, d0):
    """The minimiser of Σ_i|φ_i|²D*(|w_i|) on {Φw = y} (free masses' white selection after the gauge), by SLSQP with D*'s envelope
    gradient; D* from the inner minimisation on a log grid, interpolated."""
    from scipy.optimize import minimize, minimize_scalar
    G, g, am, m = X.G, X.gamma, X.a_max, X.m
    xs = np.geomspace(1e-9, 2.0, 400)
    ts = np.array([minimize_scalar(lambda t: G * x * x * (1 - t * t) ** 2 / (am * am * t * t) + g * x * am * t / m, bounds=(1e-9, 1.0),
                                   method='bounded', options=dict(xatol=1e-13)).x for x in xs])
    tstar = lambda x: np.interp(np.log(np.maximum(x, 1e-9)), np.log(xs), ts)
    Phi, y = X.Phi[:, :d0], X.y
    c = (Phi ** 2).sum(axis=0)

    def f(w):
        x = np.abs(w); t = tstar(x)
        val = G * x * x * (1 - t * t) ** 2 / (am * am * t * t) + g * x * am * t / m
        dv = (2 * G * x * (1 - t * t) ** 2 / (am * am * t * t) + g * am * t / m) * np.sign(w)
        return float(c @ val), c * dv
    w0 = np.linalg.lstsq(Phi, y, rcond=None)[0]
    res = minimize(f, w0, jac=True, method='SLSQP', constraints=[dict(type='eq', fun=lambda w: Phi @ w - y, jac=lambda w: Phi)],
                   options=dict(ftol=1e-14, maxiter=2000))
    return res.x


def part_A():
    """The identities at random exact fits (a 5 × 12 Gaussian design)."""
    warnings.filterwarnings('ignore')
    rng = np.random.default_rng(3)
    Phi = rng.standard_normal((N_I, D_I)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(D_I); ws[[2, 7]] = [0.2, -0.17]
    rows, canc, spec, bound = [], [], [], []
    for mode in ('conserved', 'free'):
        X, _ = law(Phi, ws, mode)
        for rep in range(3):
            beta = rng.uniform(0.5, 1.5, D_I); beta *= M_I / beta.sum()
            th = X.retract(X.theta(0.3 * rng.standard_normal(D_I), beta))
            B = np.argsort(-np.abs(X.parts(th)[3]))[:N_I]
            for h, tc in ((0.0, 0.0), (0.5, 200.0)):
                corr = ('fresh',) if tc == 0 else ('ar1', float(np.exp(-h / tc)))
                filt = (lambda k: np.ones_like(k)) if tc == 0 else (lambda k, h_=h, c_=corr: MD.LN.lr_factor(k, h_, c_))
                for coords in ('flat', 'log'):
                    g, Sf = X.scheme_drift(th, B, h=h, corr=corr, coords=coords)
                    q = X.drift(th, filt, ito=(coords == 'log')) * X.m / 4
                    rows.append(dict(mode=mode, rep=rep, tc=tc, coords=coords, rel=float(np.abs(g - q).max() / np.abs(g).max())))
                if coords == 'log':
                    al, rho, b, a, da, dda = X.parts(th)
                    om = (X.m / X.N) * X.weights(th, filt); m2 = X.m ** 2
                    gr = np.r_[om * (b / m2) * (2 * X.G * b * da * dda + 2 * X.gamma * a * da), om * (2 * X.G * b * da ** 2 / m2) * 2 * rho]
                    full = X.drift(th, filt, ito=True)
                    canc.append(dict(mode=mode, rep=rep, tc=tc, rel=float(np.abs(full + X.project(th, X.Pm * gr)).max() / np.abs(full).max())))
            A = GDR._jac(lambda x: flat_field(X, x), th, 1e-6)
            ev = np.sort(-np.linalg.eigvals(A).real)[::-1][:N_I]
            kap = np.sort(np.linalg.eigvalsh(X.gram(th)))[::-1]
            spec.append(dict(mode=mode, rep=rep, rel=float(np.abs(ev - kap).max() / kap.max())))
            if X.balanced:
                E_ = X.excess(th)
                bound.append(dict(rep=rep, E=float(E_ / X.m), trK=float(X.potential(th, lambda k: k) + X.gamma / X.N * X.y @ X.y),
                                  lb=float(X.white_bound(E_)), sba=float(X.parts(th)[2] @ X.parts(th)[4] / X.m)))
    X, _ = law(Phi, ws, 'free')
    thr = X.gamma * X.a_max ** 3 / (X.G * X.m)
    xs = np.array([1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3]) * 1.0
    ds = X.relevance_min(xs)
    return dict(rows=rows, canc=canc, spec=spec, bound=bound, thr=thr, xs=xs, ds=ds,
                lin=X.gamma * X.a_max * xs / X.m, four3=1.5 * (2 * X.G) ** (1 / 3) * (X.gamma / X.m) ** (2 / 3) * xs ** (4 / 3))


def flat_field(X, x):
    """The noiseless flow in the flat coordinates (α, ϱ) (the replicator's λ for the conserved mass)."""
    d = X.d
    _, rho, b, a, da, _ = X.parts(x)
    g = X.Phi.T @ (X.y - X.Phi @ (b * a / X.m)) / X.N
    u = a * g
    lam = b @ u / b.sum() if X.balanced else 0.0
    return np.r_[X.G * b * da * g, 0.5 * rho * X.gamma * (u - lam)]


@cached
def b_job(args):
    """Part B, one (mode, τ_c): the paradigm's noisy scheme against the reduced flow at the step H_B (the filter of the step)."""
    warnings.filterwarnings('ignore')
    mode, tc = args
    t0 = time.time()
    Phi, ws, al0 = instance(2, SEED_B)
    X, P2 = law(Phi, ws, mode)
    th0 = start(X, al0)
    rho = 0.0 if tc == 0 else float(np.exp(-H_B / tc))
    corr = ('fresh',) if tc == 0 else ('ar1', rho)
    filt = lambda k: MD.LN.lr_factor(k, H_B, corr)
    sig = SIG_B if tc == 0 else SIG_B_SLOW
    v = H_B * sig ** 2 / X.N
    S = 240.0 if tc == 0 else 2400.0
    n = int(round(S / (X.m * v / 4) / H_B)); every = n // 8
    s_rec = np.arange(9) * every * H_B * X.m * v / 4
    _, Th = X.flow(th0, s_rec[-1], filt, ito=True, t_eval=s_rec, chunks=8)
    Wf = np.array([X.w(t) for t in Th])
    _, rho0, b0, a0, _, _ = X.parts(th0)
    rec, _, _ = X.noisy(th0[:X.d], b0, H_B, sig, rho, n, np.random.default_rng(91), P_B, every)
    print(f"[{time.time() - t0:.0f}s] B {mode} tau_c {tc:g}", flush=True)
    return dict(mode=mode, tc=tc, s=s_rec, n=n, l1_paths=np.abs(rec).sum(axis=2).mean(axis=1),
                l1_se=np.abs(rec).sum(axis=2).std(axis=1, ddof=1) / np.sqrt(P_B), l1_flow=np.abs(Wf).sum(axis=1),
                dev=np.abs(rec.mean(axis=1) - Wf).max(axis=1))


@cached
def c_job(args):
    """Parts C and D, one (kind, k, instance): the ends of the reduced flows (the paradigm's update) for every mode and colour; for white
    noise also the flat update's end; basis pursuit and the minimiser of Σ|φ_i|²D*(|w_i|)."""
    warnings.filterwarnings('ignore')
    kind, k, i = args
    t0 = time.time()
    Phi, ws, al0 = instance(k, i, kind)
    d0 = Phi.shape[1]
    wbp = CS.basis_pursuit(Phi, Phi @ ws, (Phi ** 2).sum(axis=0))
    out = dict(kind=kind, k=k, i=i, bp_ok=bool(np.abs(wbp - ws).max() < 1e-4 * np.abs(ws).max()), l1_truth=float(np.abs(ws).sum()))
    for mode in MODES:
        X, P2 = law(Phi, ws, mode)
        th0 = start(X, al0)
        out[(mode, 'start')] = summary(X, th0, ws, wbp, d0)
        for tau in TAUS:
            th = X.flow(th0, S_END, filt_of(tau), ito=True, chunks=60)[1][-1]
            out[(mode, tau)] = summary(X, th, ws, wbp, d0)
        if mode in ('conserved', 'free'):
            th = X.flow(th0, S_END, filt_of(0.0), ito=False, chunks=60)[1][-1]
            out[(mode, 'flat')] = summary(X, th, ws, wbp, d0)
            if mode == 'free':
                wd = dstar_minimiser(X, d0)
                out['dstar_vs_flat'] = float(np.abs(X.w(th)[:d0] - wd).max())
                out['dstar_vs_bp'] = float(np.abs(wd - wbp).max())
    print(f"[{time.time() - t0:.0f}s] C {kind} k {k} instance {i}", flush=True)
    return out


@cached
def e_job(args):
    """Part E, one (instance, a_max): the conserved law's white end against the representability boundary, for several a_max (the
    instance kept when its basis-pursuit solution is representable, ‖w_bp‖₁ < a_max)."""
    warnings.filterwarnings('ignore')
    i, am = args
    Phi, ws, al0 = instance(2, i)
    wbp = CS.basis_pursuit(Phi, Phi @ ws, (Phi ** 2).sum(axis=0))
    if np.abs(wbp).sum() >= am:
        return dict(i=i, am=am, skip=True)
    X, _ = law(Phi, ws, 'conserved', a_max=am)
    th0 = start(X, al0 * am)
    th = X.flow(th0, S_END, filt_of(0.0), ito=True, chunks=60)[1][-1]
    return dict(i=i, am=am, skip=False, l1_bp=float(np.abs(wbp).sum()), start=summary(X, th0, ws, wbp, Phi.shape[1]),
                end=summary(X, th, ws, wbp, Phi.shape[1]))


def report(A=None, B=None, C=None, D=None, Ee=None):
    L_ = ["# Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise (`exp_mass_dictionary.py`)", "",
          f"The loop at fixed positions: masses β (Σβ = m = {M_I:g}, or free), amplitudes a = a_max tanh(α/a_max) (a_max = {AMAX:g}), w = β⊙a/m, "
          f"G = {G_I:g}, γ = {GAM_I:g}; label noise on the measured error; the selection along the exact fits (`mass_dictionary`), in s = (m v/4)t. "
          "The paradigm's update: Euler on (α, log β), renormalised; the flat one: Euler on (α, √β), renormalised radially.", ""]
    if A is not None:
        L_ += ["## A. Identities at random exact fits (5 × 12)", "",
               "The closed-form drift against Katzenberger's drift of the scheme (`gauge_drift.limit_drift` on a chart of M), largest relative "
               "difference over three points:", "", "| masses | the draws | the flat update | the paradigm's update (with its Itô term) |",
               "|---|---|---|---|"]
        for mode in ('conserved', 'free'):
            for tc in (0.0, 200.0):
                r = lambda c: max(x['rel'] for x in A['rows'] if x['mode'] == mode and x['tc'] == tc and x['coords'] == c)
                L_.append(f"| {mode} | {'fresh, h → 0' if tc == 0 else 'AR(1), τ_c = 200, h = 0.5'} | {r('flat'):.1e} | {r('log'):.1e} |")
        L_ += ["", "The paradigm's drift against the flat drift without the masses' own selection (the γ-term of D_i differentiated in √β), "
               "largest relative difference: " + '; '.join(
                   f"{mode}, {'fresh' if tc == 0 else 'AR(1)'} {max(x['rel'] for x in A['canc'] if x['mode'] == mode and x['tc'] == tc):.1e}"
                   for mode in ('free', 'conserved') for tc in (0.0, 200.0)) + ".", "",
               "K's spectrum against the flow's normal rates (−eigenvalues of its Jacobian in (α, √β)): largest relative difference "
               + ', '.join(f"{mode} {max(x['rel'] for x in A['spec'] if x['mode'] == mode):.1e}" for mode in ('conserved', 'free')) + ".", "",
               "The excess-mass bound at the conserved law's points: | E/m | Σβa′/m | tr K | (G/(Nm))E²/Σ|φ_i|⁻² |", "|---|---|---|---|"]
        for r in A['bound']:
            L_.append(f"| {r['E']:.4f} | {r['sba']:.4f} | {r['trK']:.4e} | {r['lb']:.4e} |")
        L_ += ["", f"The gauge-minimal relevance D*(|w|) against (γa_max/m)|w| and the 4/3 law (the threshold γa_max³/(Gm) = {A['thr']:.4g}):", "",
               "| |w| | D* | over (γa_max/m)|w| | over (3/2)(2G)^{1/3}(γ/m)^{2/3}|w|^{4/3} |", "|---|---|---|---|"]
        for x, dsv, li, f3 in zip(A['xs'], A['ds'], A['lin'], A['four3']):
            L_.append(f"| {x:g} | {dsv:.4e} | {dsv / li:.4f} | {dsv / f3:.4f} |")
        L_.append("")
    if B is not None:
        L_ += [f"## B. The paradigm's noisy scheme against the reduced flow (5 × 12, the 2-sparse instance {SEED_B}, step {H_B:g}, σ = {SIG_B:g} "
               f"for fresh draws and {SIG_B_SLOW:g} for AR(1) draws, {P_B} paths)", "", "| masses | the draws | steps | ‖w‖₁: the paths' mean (± s.e.) at s = 0, S/2, S | the flow's | the paths' mean "
               "w against the flow's (largest deviation, at the nine times) |", "|---|---|---|---|---|---|"]
        for r in B:
            noise = 'fresh' if r['tc'] == 0 else f"AR(1), τ_c = {r['tc']:g}"
            lp, ls, lf = r['l1_paths'], r['l1_se'], r['l1_flow']
            L_.append(f"| {r['mode']} | {noise} | {r['n']} | " + ', '.join(f'{lp[j]:.4f} ± {ls[j]:.4f}' for j in (0, 4, 8)) + " | "
                      + ', '.join(f'{lf[j]:.4f}' for j in (0, 4, 8)) + " | " + ', '.join(f'{x:.3f}' for x in r['dev']) + " |")
        L_.append("")
    for name, R, ks in (('C', C, KS_C), ('D', D, KS_D)):
        if R is None:
            continue
        hdr = (f"## C. The selection's ends on Gaussian designs ({N_I} × {D_I}, {N_INST_C} instances per k)" if name == 'C' else
               f"## D. The selection's ends on Gaussian kernel dictionaries ({D_D} centres, {N_D} samples on [0, 1], width {WIDTH_D:g}, "
               f"{N_INST_D} instances per k)")
        L_ += [hdr, "", f"The paradigm's update; reduced flows to s = {S_END:.0e}; the start the noiseless flow's end from equal masses and "
               "amplitudes N(0, 0.05²). Means over the instances: ‖w‖₁/a_max (the truth's), the excess E/m (conserved law), the atoms used "
               "(|w| > 10⁻³), alive (β > 10⁻³); the fraction recovering the truth (within 10⁻⁴·max|w*|); the median distance to basis pursuit's "
               "solution; the null kernel's mass.", "",
               "| k | masses | the noise | ‖w‖₁/a_max (truth's) | E/m | used | alive | the truth | median distance to basis pursuit | the null kernel's mass |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for k in ks:
            rs = [r for r in R if r['k'] == k]
            L_.append(f"| {k} | — | basis pursuit | | | | | {np.mean([r['bp_ok'] for r in rs]):.2f} | | |")
            for mode in MODES:
                for key in ['start'] + list(TAUS) + (['flat'] if mode != 'null kernel' else []):
                    q = [r[(mode, key)] for r in rs]
                    lab = {'start': 'the start', 'flat': 'white, the flat update', 0.0: 'white'}.get(key, None) or f"slow, τ_c = {key:.0e}"
                    L_.append(f"| {k} | {mode} | {lab} | {np.mean([x['l1'] for x in q]):.3f} ({np.mean([r['l1_truth'] for r in rs]):.3f}) | "
                              f"{np.nanmean([x['exc'] for x in q]) if mode != 'free' else float('nan'):.3f} | {np.mean([x['used'] for x in q]):.1f} | "
                              f"{np.mean([x['alive'] for x in q]):.1f} | {np.mean([x['ok'] for x in q]):.2f} | {np.median([x['d_bp'] for x in q]):.3f} | "
                              + (f"{np.mean([x['b0'] for x in q]):.2f} |" if mode == 'null kernel' else "— |"))
            L_.append(f"| {k} | free | white: the flat update's end against the minimiser of Σ|φ_i|²D*(|w_i|) (median distance) | "
                      f"{np.median([r['dstar_vs_flat'] for r in rs]):.1e} | | | | | that minimiser against basis pursuit's: "
                      f"{np.median([r['dstar_vs_bp'] for r in rs]):.1e} | |")
        L_.append("")
    if Ee is not None:
        L_ += [f"## E. The conserved law's white end against the representability boundary ({N_I} × {D_I}, 2-sparse truths, {N_INST_E} instances; "
               "those whose basis-pursuit solution is not representable, ‖w_bp‖₁ ≥ a_max, left out)", "",
               "| a_max | instances | ‖w‖₁/a_max at the start | at the end | E/m at the end | the atoms used, alive | ‖w_bp‖₁/a_max |",
               "|---|---|---|---|---|---|---|"]
        for am in AMAXES_E:
            rs = [r for r in Ee if r['am'] == am and not r['skip']]
            if not rs:
                L_.append(f"| {am:g} | 0 | | | | | |"); continue
            L_.append(f"| {am:g} | {len(rs)} | {np.mean([r['start']['l1'] for r in rs]):.3f} | {np.mean([r['end']['l1'] for r in rs]):.3f} "
                      f"(min {np.min([r['end']['l1'] for r in rs]):.3f}) | {np.mean([r['end']['exc'] for r in rs]):.4f} | "
                      f"{np.mean([r['end']['used'] for r in rs]):.1f}, {np.mean([r['end']['alive'] for r in rs]):.1f} | "
                      f"{np.mean([r['l1_bp'] / am for r in rs]):.3f} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_mass_dictionary.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mass_dictionary_raw.npz')
    md = os.path.join(OUT, 'mass_dictionary.md')
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
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(mo, tc) for mo in MODES for tc in (0.0, 1e3)], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [('gauss', k, i) for k in KS_C for i in range(N_INST_C)], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(c_job, [('rbf', k, i) for k in KS_D for i in range(N_INST_D)], chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(e_job, [(i, am) for i in range(N_INST_E) for am in AMAXES_E], chunksize=1)
        if 'A' in parts:
            res['A'] = part_A(); E.save(raw, **res); log('A')
        for k in 'BCDE':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE']))


if __name__ == '__main__':
    main()
