"""
Proposition D‴.38: a Fixman correction for the replicator under dither — the entropic drift cancelled, or reversed toward the best-conditioned
fit; the saturation made non-attracting
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
   exponent d log p/d log δ of the one-dimensional zero-flux density p ∝ D⁻¹exp∫2b/D for κ_F = 0, 1, 2 (predicted κ_F − 1).
D. The full loop over 10⁶ slow units at γ = 10⁻², the plug-in form (12 seeds, the seeds of exp_gauge_drift.py's Part H, whose records are the
   κ_F = 0 runs): κ_F = 1, 2 — the excursions into the saturation by quarters, the fraction of records beyond the saturation, the masses'
   means and spreads at 10⁵, 3·10⁵, 10⁶, the second half's mean against W♯, the fit error against equipartition's √(pτ/(m⟨F²⟩)).
E. As D with the masses form, Â fixed at the initial exact fit (the same seeds).
F. With the entropic leak at γ = 3·10⁻³, ε = 10⁻⁴ (16 seeds, 4·10⁵ slow units, the second half): no correction, the plug-in form at
   κ_F = 1, 2, the masses form at κ_F = 1, 2 with Â fixed and at κ_F = 2 with Â tracked (T_A = 10⁵) — the seeds whose amplitudes escape
   up the positions' plateau (some |α_k| > 10), the masses' range, the pinned state's shift (all seeds, and the seeds without an escape)
   against κ⁻¹b (b the corrected drift at w = 1/3, κ the leak's stiffness along M), the fit error.
G. With the entropic leak at γ = 10⁻², ε = 3·10⁻⁴ (inside Proposition D‴.37(iv)'s window; the 16 seeds of exp_gauge_drift.py's Part G, whose
   κ_F = 0 runs this repeats bit for bit), the masses form at κ_F = 1, 2: the shift against κ⁻¹b_c, and the paired differences from κ_F = 0.

Output: results/chaos/gauge_fixman.md and gauge_fixman_raw.npz (GF_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`).
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
GAMMA_D, GAINS_D, S_D, SEEDS_D = 1e-2, (1.0, 2.0), 1e6, X.SEEDS_H
T_D = X.T_H
GAMMA_F, EPS_F, SEEDS_F, S_F, TA_F = 3e-3, 1e-4, tuple(range(400, 416)), 4e5, 1e5
MODES_F = (('none', 0.0), ('plug-in', 1.0), ('plug-in', 2.0), ('masses', 1.0), ('masses', 2.0), ('tracked', 2.0))
ESCAPE = 10.0                                                              # |α_k| beyond it: up the positions' plateau (|a_k|/a_max > 0.986)
GAMMA_G, EPS_G, GAINS_G, SEEDS_G, S_G = 1e-2, 3e-4, (0.0, 1.0, 2.0), tuple(range(200, 216)), 4e5   # exp_gauge_drift.py's Part G at ε = 3·10⁻⁴
REC = GD.REC


def instance(gamma, kF):
    """X.instance's flow with the correction of gain κ_F, and F₁ = the corrected flow minus the uncorrected one."""
    net, base = E.setup()
    loops = [MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, fixman=k, fixman_tau=TAU) for k in (kF, 0.0)]

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


def drift(w, gamma, kF):
    F0, F1, chart, coords, S = instance(gamma, kF)
    return GDR.limit_drift(F0, chart, coords, np.asarray(w, float), S, h=X.H_ID, h2=X.H2_ID, F1=(F1 if kF else None))


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
    kF, dl = args
    W1 = W_MIN[0] + dl; w = np.array([W1, (1 - W1) / 2]); hh = 0.05 * dl
    b, D, _ = drift(w, GAMMA_C, kF)
    Dp = drift(w + np.array([hh, -hh / 2]), GAMMA_C, kF)[1]; Dm = drift(w - np.array([hh, -hh / 2]), GAMMA_C, kF)[1]
    dlogD = (np.log(Dp[0, 0]) - np.log(Dm[0, 0])) / (2 * hh)
    return dict(kF=kF, delta=dl, expo=dl * 2 * b[0] / D[0, 0] - dl * dlogD)


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
               half=W[:, R.shape[1] // 2:].mean(axis=1), fit=R[:, R.shape[1] // 2:, 7].mean())
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
    pred = {}
    out = []
    for form, kF in MODES_F:
        R = np.asarray(rec[(form, kF)], float); h = R.shape[1] // 2
        k_eff = 0.0 if form == 'none' else kF
        if k_eff not in pred:
            c = np.linalg.solve(kappa(w, EPS_F, GAMMA_F), drift(w[:2], GAMMA_F, k_eff)[0]); pred[k_eff] = np.r_[c, -c.sum()]
        esc = np.abs(R[:, :, 1:4]).max(axis=(1, 2)) > ESCAPE
        first = [float(R[i, np.argmax(np.abs(R[i, :, 1:4]).max(axis=1) > ESCAPE), 0]) for i in np.flatnonzero(esc)]
        means = R[:, h:, 4:7].mean(axis=1); ok = ~esc
        se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M)) if len(M) > 1 else np.full(3, np.nan)
        out.append(dict(form=form, kF=kF, seeds_esc=[SEEDS_F[i] for i in np.flatnonzero(esc)], first=first,
                        Wmin=R[:, :, 4:7].min(axis=(0, 1)), Wmax=R[:, :, 4:7].max(axis=(0, 1)), almax=np.abs(R[:, :, 1:4]).max(),
                        shift=means.mean(axis=0) - 1 / 3, se=se(means), shift_ok=means[ok].mean(axis=0) - 1 / 3, se_ok=se(means[ok]),
                        std_ok=R[ok][:, h:, 4:7].std(axis=1).mean(axis=0), fit=R[:, h:, 7].mean(), fit_ok=R[ok][:, h:, 7].mean(),
                        fit_esc=(R[esc][:, h:, 7].mean() if esc.any() else np.nan), pred=pred[k_eff],
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
           "the second half's mean (± s.e. over the seeds) | W♯ | the fit error (second half) |", "|---|---|---|---|---|---|---|"]
    for r in D:
        sm = r['sum']; hm = sm['half']
        L_.append(f"| {r['kF']:g} | {'; '.join(', '.join(f'{v:.1f}' for v in q) for q in sm['quarters'])} | {', '.join(f'{v:.2f}' for v in sm['big'][-1])} | "
                  f"{', '.join(f'{v:.2f}' for v in sm['below'])} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(hm.mean(axis=0), hm.std(axis=0, ddof=1) / np.sqrt(len(hm))))} | "
                  f"{', '.join(f'{v:.4f}' for v in ws)} | {fmt(sm['fit'])} |")
    return L_


def part_G(rec):
    """rec: {κ_F: records (seeds, n, 8 or 11)} (κ_F = 0 uncorrected, else the masses form, the same seeds). The second half's shift of the
    masses from 1/3 (± s.e. over the seeds) against κ⁻¹b_c(w), and the paired differences from κ_F = 0 against κ⁻¹(b_c − b)."""
    w = np.ones(3) / 3
    K = kappa(w, EPS_G, GAMMA_G)
    pred = {}
    for k in rec:
        c = np.linalg.solve(K, drift(w[:2], GAMMA_G, k)[0]); pred[k] = np.r_[c, -c.sum()]
    R0 = np.asarray(rec[0.0], float); h = R0.shape[1] // 2; m0 = R0[:, h:, 4:7].mean(axis=1)
    out = []
    for k in sorted(rec):
        R = np.asarray(rec[k], float); means = R[:, h:, 4:7].mean(axis=1); d = means - m0
        esc = np.abs(R[:, :, 1:4]).max(axis=(1, 2)) > ESCAPE
        out.append(dict(kF=k, shift=means.mean(axis=0) - 1 / 3, se=means.std(axis=0, ddof=1) / np.sqrt(len(means)), pred=pred[k],
                        diff=d.mean(axis=0), dse=d.std(axis=0, ddof=1) / np.sqrt(len(d)), dpred=pred[k] - pred[0.0],
                        esc=[SEEDS_G[i] for i in np.flatnonzero(esc)], std=R[:, h:, 4:7].std(axis=1).mean(axis=0), fit=R[:, h:, 7].mean()))
    return out


def report(A=None, B=None, C=None, D=None, Ee=None, Ff=None, Gg=None):
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
               f"(|a|/a_max there {', '.join(f'{v:.3f}' for v in W_MIN / ws)}).", "",
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
        L_ += ["", "Equipartition: at fixed masses the positions' law is Gibbs, each position coordinate that holds fit energy stores τ/2 of "
               "J = (m/2)⟨r²⟩, so the fit error is √(pτ/(m⟨F²⟩)) with p such coordinates — " + ', '.join(f"{p_}: {np.sqrt(p_ * TAU / (3 * F2)):.4f}" for p_ in (9, 8, 7))
               + f" (⟨F²⟩ = {F2:.4f}); a saturated amplitude (a′ → 0) holds none. Read backwards, the dither's level from the fit error e: "
               "τ̂ = m⟨F²⟩e²/9 = " + ', '.join(f"{m_ * F2 * r['sum']['fit'] ** 2 / 9:.3g} (κ_F = {r['kF']:g})" for r in D for m_ in (3,))
               + f" against τ = {TAU:g}."]
        L_.append("")
    if Ee is not None:
        L_ += long_tables(Ee, f"## E. The full loop over {S_D:g} slow units, the masses form with Â fixed at the initial exact fit (γ = {GAMMA_D:g}, "
                              f"the seeds of D)")
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
               "the within-seed std (without an escape) |", "|---|---|---|---|---|---|"]
        for r in Ff:
            L_.append(f"| {r['form']} | {r['kF']:g} | {', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift'], r['se']))} | "
                      f"{', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift_ok'], r['se_ok']))} | {', '.join(f'{v:+.5f}' for v in r['pred'])} | "
                      f"{', '.join(f'{v:.4f}' for v in r['std_ok'])} |")
        tr = [r for r in Ff if r['form'] == 'tracked']
        if tr:
            Ae = tr[0]['A_end']; esc = np.isin(SEEDS_F, tr[0]['seeds_esc'])
            L_ += ["", f"The tracked contributions Â at {S_F:g} (the bumps' amplitudes {', '.join(f'{v:+.3f}' for v in E.A_BUMP)}): the seeds without an "
                   f"escape {', '.join(f'{v:+.3f}' for v in Ae[~esc].mean(axis=0))} (mean)"
                   + (f"; with one {'; '.join(', '.join(f'{v:+.3f}' for v in a) for a in Ae[esc])}" if esc.any() else '') + "."]
        L_.append("")
    if Gg is not None:
        L_ += [f"## G. With the entropic leak at γ = {GAMMA_G:g}, ε = {EPS_G:g} (the seeds of exp_gauge_drift.py's Part G, {S_G:g} slow units, the second "
               "half): the masses form, and its difference from κ_F = 0 on the same seeds", "",
               "| κ_F | the masses' mean − 1/3 (± s.e.) | predicted κ⁻¹b_c | the difference from κ_F = 0, the same seeds (± s.e.) | predicted κ⁻¹(b_c − b) | "
               "the within-seed std | seeds with an escape | the fit error |", "|---|---|---|---|---|---|---|---|"]
        for r in Gg:
            L_.append(f"| {r['kF']:g} | {', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift'], r['se']))} | {', '.join(f'{v:+.5f}' for v in r['pred'])} | "
                      f"{', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['diff'], r['dse']))} | {', '.join(f'{v:+.5f}' for v in r['dpred'])} | "
                      f"{', '.join(f'{v:.4f}' for v in r['std'])} | {', '.join(str(x) for x in r['esc']) or 'none'} | {fmt(r['fit'])} |")
        zg = E.load(os.path.join(OUT, 'gauge_drift_raw.npz')).get('G')
        old = [g for g in (zg if zg is not None else []) if g['eps'] == EPS_G]
        if old:
            d = np.abs(Gg[0]['shift'] - old[0]['shift']).max()
            L_ += ["", f"κ_F = 0 repeats exp_gauge_drift.py's Part G at ε = {EPS_G:g} (the same seeds and steps): its stored shift "
                   f"{', '.join(f'{v:+.5f}' for v in old[0]['shift'])}, the largest difference {d:.1e}."]
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_gauge_fixman.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_fixman_raw.npz')
    md = os.path.join(OUT, 'gauge_fixman.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEFG']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFG')
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
        if 'D' in parts:
            h['D'] = pool.map_async(series, [(GAMMA_D, 0.0, k, s, S_D) for k in GAINS_D for s in SEEDS_D], chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(series, [(GAMMA_D, 0.0, k, s, S_D, 'masses') for k in GAINS_D for s in SEEDS_D], chunksize=1)
        if 'F' in parts:
            h['F'] = pool.map_async(series, [(GAMMA_F, EPS_F, k, s, S_F, f) for f, k in MODES_F for s in SEEDS_F], chunksize=1)
        if 'G' in parts:
            h['G'] = pool.map_async(series, [(GAMMA_G, EPS_G, k, s, S_G, 'none' if k == 0 else 'masses') for k in GAINS_G for s in SEEDS_G],
                                    chunksize=1)
        for k in 'ABC':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
        z = E.load(os.path.join(OUT, 'gauge_drift_raw.npz')) if ('D' in h or 'E' in h) else None
        for k, key in (('D', 'D_rec'), ('E', 'E_rec')):
            if k in h:
                runs = h[k].get()
                res[key] = {g: np.array([r['rec'] for r in runs if r['kF'] == g]) for g in GAINS_D}
                res[k] = [dict(kF=0.0, sum=summarise_long(np.asarray(z['H_rec'], float)))] + \
                         [dict(kF=g, sum=summarise_long(res[key][g])) for g in GAINS_D]
                E.save(raw, **res); log(k)
        if 'F' in h:
            runs = h['F'].get()
            res['F_rec'] = {(f, k): np.array([r['rec'] for r in runs if r['form'] == f and r['kF'] == k]) for f, k in MODES_F}
            res['F'] = part_F(res['F_rec']); E.save(raw, **res); log('F')
        if 'G' in h:
            runs = h['G'].get()
            res['G_rec'] = {k: np.array([r['rec'] for r in runs if r['kF'] == k]) for k in GAINS_G}
            res['G'] = part_G(res['G_rec']); E.save(raw, **res); log('G')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEFG']))


if __name__ == '__main__':
    main()
