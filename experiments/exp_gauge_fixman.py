"""
Proposition D‴.38: a Fixman correction for the replicator under dither — the entropic drift cancelled, or reversed toward the best-conditioned
fit; the saturation made non-attracting
=====================================================================================================================================

Proposition D‴.36's instance (Proposition D‴.32's three bumps, one particle per class at an exact representation, the replicator at γ, the
positions dithered at τ = 10⁻⁴). The correction adds κ_F(τ/2)∂_{β_k}ℓ to the replicator's fitness, ℓ = log det ∇²_ΘJ along the manifold M of
exact fits, in Proposition D‴.37's closed form at the current amplitudes: κ_F(τ/m)(1 + (a_k/a_max)²)/(W_k a′_k)
(`ReplicatorCohesionLoop(fixman=κ_F, fixman_tau=τ)`). On the slow time scale an O(τ) perturbation F₁ of the drift adds ∂W·P·F₁ to
Katzenberger's drift (`gauge_drift.limit_drift(F1=…)`); here F₁ is the corrected flow minus the uncorrected one.

A. The cancellation (κ_F = 1): the uncorrected drift over γτ against the log-determinant law, and the corrected drift over γτ, at three
   points for γ = 2·10⁻², …, 10⁻⁴.
B. The gain: the zero of the corrected field (Newton from W♯) and the eigenvalues of its Jacobian for κ_F = 1, 2 at γ = 10⁻², 10⁻³, 10⁻⁴,
   against W♯ = argmax ℓ on the simplex (ℓ = 2Σ log(W_k − W_min,k²/W_k), strictly concave).
C. The saturation: along W_2 = W_3 near the first class's saturation (W_1 − 0.25 = 0.002–0.004, γ = 10⁻², inside its layer), the
   exponent d log p/d log δ of the one-dimensional zero-flux density p ∝ D⁻¹exp∫2b/D for κ_F = 0, 1, 2 (predicted κ_F − 1).
D. The full loop over 10⁶ slow units at γ = 10⁻² (12 seeds, the seeds of exp_gauge_drift.py's Part H, whose records are the κ_F = 0 runs):
   κ_F = 1, 2 — the excursions into the saturation by quarters, the fraction of records beyond the saturation, the masses' means and
   spreads at 10⁵, 3·10⁵, 10⁶, the second half's mean against W♯, the fit error against equipartition's √(pτ/(m⟨F²⟩)).
E. With the entropic leak at γ = 3·10⁻³, ε = 10⁻⁴ (16 seeds, 4·10⁵ slow units, the second half): the pinned state's shift for κ_F = 0, 1, 2
   against κ⁻¹b (b the corrected drift at w = 1/3, κ the leak's stiffness along M).

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
GAMMA_E, EPS_E, GAINS_E, SEEDS_E, S_E = 3e-3, 1e-4, (0.0, 1.0, 2.0), tuple(range(400, 416)), 4e5
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
    return dict(gamma=gamma, kF=kF, zero=np.r_[z, 1 - z.sum()], ier=int(ier), eig=np.linalg.eigvals(Jb), resid=np.abs(f(z)).max())


def c_job(args):
    kF, dl = args
    W1 = W_MIN[0] + dl; w = np.array([W1, (1 - W1) / 2]); hh = 0.05 * dl
    b, D, _ = drift(w, GAMMA_C, kF)
    Dp = drift(w + np.array([hh, -hh / 2]), GAMMA_C, kF)[1]; Dm = drift(w - np.array([hh, -hh / 2]), GAMMA_C, kF)[1]
    dlogD = (np.log(Dp[0, 0]) - np.log(Dm[0, 0])) / (2 * hh)
    return dict(kF=kF, delta=dl, expo=dl * 2 * b[0] / D[0, 0] - dl * dlogD)


def series(args):
    """GD.series with the correction: records (s, α_1..3, W_1..3, the fit error) every REC steps."""
    gamma, eps, kF, seed, S = args
    t0 = time.time()
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3).copy()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, eps=eps, fixman=kF, fixman_tau=TAU)
    rec = []
    loop.run_dithered(Th, np.ones(3), S, DT, TAU, np.random.default_rng(seed), rec_every=REC,
                      callback=lambda s_, T, b: rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b)]))
    print(f"[{time.time() - t0:.0f}s] gamma {gamma:g} eps {eps:g} kF {kF:g} seed {seed}", flush=True)
    return dict(gamma=gamma, eps=eps, kF=kF, seed=seed, rec=np.array(rec))


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


def part_E(runs):
    w = np.ones(3) / 3
    out = []
    for kF in GAINS_E:
        Rr = np.array([r['rec'] for r in runs if r['kF'] == kF]); h = Rr.shape[1] // 2
        means = Rr[:, h:, 4:7].mean(axis=1)
        b = drift(w[:2], GAMMA_E, kF)[0]
        c = np.linalg.solve(kappa(w, EPS_E, GAMMA_E), b)
        out.append(dict(kF=kF, shift=means.mean(axis=0) - 1 / 3, se=means.std(axis=0, ddof=1) / np.sqrt(len(means)), pred=np.r_[c, -c.sum()],
                        std=Rr[:, h:, 4:7].std(axis=(0, 1)), fit=Rr[:, h:, 7].mean()))
    return out


def report(A=None, B=None, C=None, D=None, Ee=None):
    fmt = E.fmt
    w3 = lambda w: ', '.join(f'{v:.3f}' for v in np.r_[w, 1 - sum(w)])
    L_ = ["# Proposition D‴.38: a Fixman correction for the replicator under dither (`exp_gauge_fixman.py`)", "",
          f"Proposition D‴.36's instance (three bumps, one particle per class at an exact representation, m = 3, the positions dithered at τ = {TAU:g}); "
          "the correction adds κ_F(τ/m)(1 + (a_k/a_max)²)/(W_k a′_k) to the replicator's fitness; drifts per unit slow time, W_k = β_k/3.", ""]
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
               "| κ_F | γ | the zero (Newton from W♯) | its distance from W♯ | the Jacobian's eigenvalues (per unit slow time) | the residual |", "|---|---|---|---|---|---|"]
        for r in B:
            L_.append(f"| {r['kF']:g} | {r['gamma']:g} | {', '.join(f'{v:.5f}' for v in r['zero'])}{'' if r['ier'] == 1 else ' (no convergence)'} | "
                      f"{np.abs(r['zero'] - ws).max():.1e} | {', '.join(fmt(v.real) for v in r['eig'])} | {r['resid']:.1e} |")
        L_.append("")
    if C is not None:
        L_ += [f"## C. The saturation: the density's exponent near the first class's saturation (γ = {GAMMA_C:g}, W_2 = W_3)", "",
               "| κ_F | d log p / d log δ at δ = " + ', '.join(f'{d:g}' for d in DELTAS_C) + " | predicted |", "|---|---|---|"]
        for kF in GAINS_C:
            ex = [r['expo'] for r in C if r['kF'] == kF]
            L_.append(f"| {kF:g} | {', '.join(f'{v:+.3f}' for v in ex)} | {kF - 1:+g} |")
        L_.append("")
    if D is not None:
        ws = w_sharp()
        L_ += [f"## D. The full loop over {S_D:g} slow units (γ = {GAMMA_D:g}, {len(SEEDS_D)} seeds; κ_F = 0 from exp_gauge_drift.py's Part H, the same seeds)", "",
               "| κ_F | s | mean of W_1, W_2, W_3 (± s.e.) | std | the largest |a|/a_max |", "|---|---|---|---|---|"]
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
        net, base = E.setup()
        F2 = float(np.mean(base.F ** 2))
        L_ += ["", "Equipartition: at fixed masses the positions' law is Gibbs, each position coordinate that holds fit energy stores τ/2 of "
               "J = (m/2)⟨r²⟩, so the fit error is √(pτ/(m⟨F²⟩)) with p such coordinates — " + ', '.join(f"{p_}: {np.sqrt(p_ * TAU / (3 * F2)):.4f}" for p_ in (9, 8, 7))
               + f" (⟨F²⟩ = {F2:.4f}); a saturated amplitude (a′ → 0) holds none."]
        L_.append("")
    if Ee is not None:
        L_ += [f"## E. With the entropic leak (γ = {GAMMA_E:g}, ε = {EPS_E:g}, {len(SEEDS_E)} seeds, {S_E:g} slow units, the second half)", "",
               "| κ_F | the masses' mean − 1/3 (± s.e.) | predicted κ⁻¹b | std | the fit error |", "|---|---|---|---|---|"]
        for r in Ee:
            L_.append(f"| {r['kF']:g} | {', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(r['shift'], r['se']))} | {', '.join(f'{v:+.5f}' for v in r['pred'])} | "
                      f"{', '.join(f'{v:.4f}' for v in r['std'])} | {fmt(r['fit'])} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_gauge_fixman.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_fixman_raw.npz')
    md = os.path.join(OUT, 'gauge_fixman.md')
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
            h['A'] = pool.map_async(a_job, [(w, g) for w in W_A for g in GAMMAS_A], chunksize=1)
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(g, k) for k in GAINS_B for g in GAMMAS_B], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(k, d) for k in GAINS_C for d in DELTAS_C], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(series, [(GAMMA_D, 0.0, k, s, S_D) for k in GAINS_D for s in SEEDS_D], chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(series, [(GAMMA_E, EPS_E, k, s, S_E) for k in GAINS_E for s in SEEDS_E], chunksize=1)
        for k in 'ABC':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
        if 'D' in h:
            runs = h['D'].get()
            z = E.load(os.path.join(OUT, 'gauge_drift_raw.npz'))
            res['D_rec'] = {k: np.array([r['rec'] for r in runs if r['kF'] == k]) for k in GAINS_D}
            res['D'] = [dict(kF=0.0, sum=summarise_long(np.asarray(z['H_rec'], float)))] + \
                       [dict(kF=k, sum=summarise_long(res['D_rec'][k])) for k in GAINS_D]
            E.save(raw, **res); log('D')
        if 'E' in h:
            res['E'] = part_E(h['E'].get()); E.save(raw, **res); log('E')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE']))


if __name__ == '__main__':
    main()
