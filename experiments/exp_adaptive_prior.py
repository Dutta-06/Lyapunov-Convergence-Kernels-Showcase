"""
Theorem D‴.23: adaptive priors — a learned prior removes the leak's bias, at the rate of the measure-space excitation
==================================================================================================================

The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit with the OT leak toward a prior whose centres are shifted by +0.15
(a misplaced prior: the fixed-prior stationary point Θ*(q) is biased), masses q = 1/3, and the adaptive prior of src/core/adaptive_prior.py.

A. The linear theory at the exact representation: the whitened weighted Gauss–Newton spectrum h_j of the loss, the fraction σ₀/(h_j + σ₀)
   of the prior's error that a fixed prior keeps in mode j, and the rates λ_j^− of the collective pairs for η = σ₀/10, σ₀, 10σ₀, ∞.
B. The fixed prior's bias: Θ*(q) by Newton against the first-order prediction y* = σ₀(H_L + σ₀I)⁻¹y⁰, mode by mode.
C. The K₀-atom adaptive loop (one particle per class; every replicated lift carries it) from Θ*(q) with Z = the shifted prior, for the four
   η, 6000 slow units: the mode amplitudes of the atoms' deviation from the exact representation against the exact solutions of the pairs,
   the fit error, the energy.
D. m-particle systems: m = 15, 60, 240 and M = 510 from the warm clouds of Proposition D‴.12 (four seeds), η = 0 (the fixed prior), σ₀ and ∞,
   3000 slow units, ungated and with two gates — 'transit' (the prior frozen while some particle is assigned to an atom that is not its
   nearest) and 'settled' (in addition, the first time, until the fixed-prior drift is below 10⁻⁵): the fit error, W²(μ_Θ, μ_exact), the
   class means' whitened weighted deviation from the exact representation against the K₀-atom run of Part C (per seed at s − t_open for
   the settled gate), the prior's distance, the final counts. `--only-gated transit|settled` recomputes one gate and merges it.
E. The ungated η = ∞ system with the largest ‖ȳ‖ at s = 3000, re-run from its seed and continued to 23000 slow units: the fit error, the
   energy, the class-centre spread and the drift; at the end the unregularised loss's Hessian at the collapsed atoms and a BFGS descent of
   that loss from them (`--only-E` recomputes Part E and merges it).

Output: results/chaos/adaptive_prior.md and adaptive_prior_raw.npz (AP_SEEDS, AP_WORKERS; `--from-raw` regenerates the md).
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
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import ot_leak as OT
from src.core import adaptive_prior as AP

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
SHIFT = 0.15
Q = np.ones(3) / 3
ETAS = (SIGMA0 / 10, SIGMA0, 10 * SIGMA0, np.inf)
ETAS_D = (0.0, SIGMA0, np.inf)
M_REF, MS = 510, (15, 60, 240)
JITTER = 0.02
S_C, S_D, DT, PLAN_EVERY = 6000.0, 3000.0, 0.25, 4
REC_C, REC_D = 100, 100                                                   # records every 25 slow units
FLOOR = 1e-6                                                              # mode amplitudes below it are at the second-order remainder's level
N_SEEDS = int(os.environ.get('AP_SEEDS', 4))
N_WORKERS = int(os.environ.get('AP_WORKERS', 3))


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def representation(net, shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def setup():
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    exact = representation(net, 0.0)
    prior = representation(net, SHIFT)
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, prior, n_quad=N_Q)
    return net, exact, prior, loop


def class_means(Th, gamma, q):
    return AP.barycentres(Th, gamma, q)


def linear_parts():
    net, exact, prior, loop = setup()
    HL, h, V = AP.gauss_newton_block(loop, exact, Q)
    star, res, marg = loop.stationary_atoms(prior)                          # the fixed prior's stationary point
    y0 = AP.whitened(prior - exact, Q, G)
    ya = AP.whitened(star - exact, Q, G)
    yp = SIGMA0 * np.linalg.solve(HL + SIGMA0 * np.eye(HL.shape[0]), y0)
    fit_star = loop.fit_error(np.repeat(star, 1, axis=0))
    rates = {eta: AP.pair_rates(h, SIGMA0, eta)[0] for eta in ETAS}
    return dict(h=h, V=V, HL=HL, star=star, newton_res=res, diag=marg, y0=y0, ya=ya, yp=yp, fit_star=loop.fit_error(star), rates=rates, exact=exact, prior=prior)


def k0_job(eta):
    """Part C: the K₀-atom adaptive loop from Θ*(q) with Z = the shifted prior."""
    net, exact, prior, loop = setup()
    lin = linear_parts()
    apl = AP.AdaptivePriorLoop(loop, eta)
    rec = dict(s=[], y=[], z=[], fit=[], J=[])

    def cb(s, Th, Z, adapting):
        rec['s'].append(s); rec['y'].append(AP.whitened(Th - exact, Q, G)); rec['z'].append(AP.whitened(Z - exact, Q, G))
        rec['fit'].append(loop.fit_error(Th)); rec['J'].append(apl.energy(Th, Z))
    apl.run(lin['star'].copy(), prior.copy(), S_C, DT, plan_every=1, rec_every=REC_C, callback=cb)
    return dict(eta=eta, **{k: np.array(v) for k, v in rec.items()})


def m_job(args):
    """Part D: the m-particle systems from the warm clouds."""
    seed, eta, gate = args
    t0 = time.time()
    net, exact, prior, loop = setup()
    rng = np.random.default_rng(seed)
    clouds = CL.mixture_start(exact, M_REF, rng, JITTER, G)
    out = dict(seed=seed, eta=eta, gate=gate, rec={}, t_open={})
    for m in list(MS) + [M_REF]:
        apl = AP.AdaptivePriorLoop(loop, eta if eta > 0 else 0.0)
        rec = dict(s=[], fit=[], W2ex=[], ybar=[], zdist=[], J=[], counts=[], adapting=[])

        def cb(s, Th, Z, adapting, rec=rec, apl=apl):
            gamma, W2 = AP.plan_to(Th, Z, Q, G)
            rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['W2ex'].append(OT.W2_empirical(Th, exact, G))
            rec['ybar'].append(AP.whitened(class_means(Th, gamma, Q) - exact, Q, G)); rec['zdist'].append(float(np.linalg.norm(AP.whitened(Z - exact, Q, G))))
            rec['J'].append(loop.L_bar(Th) + 0.5 * SIGMA0 * W2); rec['counts'].append(np.bincount(np.argmin(loop.cost_matrix(Th), axis=1), minlength=3))
            rec['adapting'].append(adapting)
        _, _, t_open = apl.run(clouds[:m].copy(), prior.copy(), S_D, DT, plan_every=PLAN_EVERY, rec_every=REC_D, callback=cb, gate=gate)
        out['rec'][m] = {k: np.array(v) for k, v in rec.items()}
        out['t_open'][m] = t_open
    print(f"[{time.time() - t0:.0f}s] seed {seed} eta {eta:g} gate {gate}", flush=True)
    return out


S_E, REC_E = 23000.0, 4000                                                 # Part E: continued to 23000 slow units, records every 1000


def trapped_job(D):
    """Part E: the ungated η = ∞ system with the largest ‖ȳ‖ at s = 3000 in Part D, re-run from its seed (the run is deterministic: the fit
    error at s = 3000 is compared with the stored one) and continued to S_E; every 2000 slow units the fit error, the energy, the largest
    spread of the class centres and the largest drift; at the end the collapsed class atoms, the eigenvalues of the unregularised loss's
    Hessian there (whitened, weighted) and a BFGS minimisation of the unregularised weighted loss from them."""
    from scipy.optimize import minimize
    cand = [(np.linalg.norm(d['rec'][m]['ybar'][int(np.argmin(np.abs(d['rec'][m]['s'] - 3000.0)))]), d['seed'], m, d)
            for d in D if np.isinf(d['eta']) and d.get('gate', False) is False for m in d['rec']]
    yb, seed, m, d = max(cand, key=lambda t: t[0])
    net, exact, prior, loop = setup()
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    apl = AP.AdaptivePriorLoop(loop, np.inf)
    rec = dict(s=[], fit=[], J=[], spread=[], drift=[])

    def cb(s, Th, Z, adapting):
        gamma, _ = AP.plan_to(Th, Z, Q, G)
        lab = np.argmax(gamma, axis=1); c = loop.features(Th)[1]
        rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['J'].append(apl.energy(Th, Z))
        rec['spread'].append(max(float(np.ptp(c[lab == j])) for j in range(3))); rec['drift'].append(float(np.abs(apl.drift(Th, Z, gamma)[0]).max()))
    Th, Z, _ = apl.run(clouds[:m].copy(), prior.copy(), S_E, DT, plan_every=PLAN_EVERY, rec_every=REC_E, callback=cb, gate=False)
    gamma, _ = AP.plan_to(Th, Z, Q, G); lab = np.argmax(gamma, axis=1)
    atoms = np.array([Th[lab == j].mean(axis=0) for j in range(3)])

    def L(v):
        at = v.reshape(3, 3); a_, c_, s_, d_, phi = loop.features(at)
        return 0.5 * float(np.mean((loop.F - phi @ (Q * a_)) ** 2))
    res = minimize(L, atoms.ravel(), method='BFGS', options=dict(gtol=1e-13, maxiter=20000))
    ev = AP.gauss_newton_block(loop, atoms, Q)[1]
    stored_fit = d['rec'][m]['fit'][int(np.argmin(np.abs(d['rec'][m]['s'] - 3000.0)))]
    feat = lambda at: tuple(np.asarray(x) for x in loop.features(at)[:3])
    return dict(seed=seed, m=m, ybar3000=yb, stored_fit3000=stored_fit, rec={k: np.array(v) for k, v in rec.items()}, atoms=feat(atoms),
                L_atoms=L(atoms.ravel()), hess_ev=ev, bfgs_L=float(res.fun), bfgs_atoms=feat(res.x.reshape(3, 3)), exact=feat(exact))


def report(lin, C, D, E=None, t0=None):
    h, V = lin['h'], lin['V']
    lines = ["# Theorem D‴.23: adaptive priors — a learned prior removes the leak's bias (`exp_adaptive_prior.py`)", "",
             f"The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f}; the OT leak toward the prior with centres shifted by +{SHIFT:g} (masses 1/3); "
             f"the adaptive prior moves its atoms toward the class barycentres at the rate η. Midpoint rule dt = {DT:g} slow units.", "",
             "## A. The linear theory at the exact representation", "",
             "The whitened weighted Gauss–Newton eigenvalues h_j of the loss (the collective block without the leak), the fraction σ₀/(h_j + σ₀) of the prior's error that a fixed prior "
             "keeps in mode j, and the slower rate λ_j^− of the collective pair for each η (η = ∞: λ_j^− = h_j). The separations decay at σ₀ in every case.", "",
             "| j | h_j | σ₀/(h_j + σ₀) | " + " | ".join(f"λ_j^− at η = {('∞' if np.isinf(e) else f'{e:.5f}')}" for e in ETAS) + " |",
             "|---|---|---|" + "---|" * len(ETAS)]
    for j in range(len(h)):
        lines.append(f"| {j + 1} | {h[j]:.4e} | {SIGMA0 / (h[j] + SIGMA0):.3f} | " + " | ".join(f"{lin['rates'][e][j]:.3e}" for e in ETAS) + " |")
    lines += ["", f"A fixed prior's collective rates are h_j + σ₀ = {h[0] + SIGMA0:.4e}–{h[-1] + SIGMA0:.4e}.", "",
              "## B. The fixed prior's bias (shifted prior)", "",
              f"Θ*(q) by Newton (residual {lin['newton_res']:.1e}; diagonal plan ε = {lin['diag'][0]:.3f} < Δ_ν/2 = {lin['diag'][1]:.3f}); fit error at Θ*(q) {lin['fit_star']:.4e} (0 at the exact representation). "
              f"Mode amplitudes (components on the eigenvectors of H_L, whitened weighted coordinates) of the prior's error y⁰, of the actual bias y = Θ*(q) − exact, and of the first-order prediction "
              f"σ₀(H_L + σ₀I)⁻¹y⁰; ‖y − prediction‖/‖y‖ = {np.linalg.norm(lin['ya'] - lin['yp']) / np.linalg.norm(lin['ya']):.3e}.", "",
              "| j | y⁰_j | y_j (actual) | predicted | actual/y⁰ |", "|---|---|---|---|---|"]
    a0, aa, ap = V.T @ lin['y0'], V.T @ lin['ya'], V.T @ lin['yp']
    for j in range(len(h)):
        lines.append(f"| {j + 1} | {a0[j]:+.4e} | {aa[j]:+.4e} | {ap[j]:+.4e} | {aa[j] / a0[j]:.3f} |")
    # C
    lines += ["", "## C. The K₀-atom adaptive loop from Θ*(q) with the prior at the shifted centres", "",
              "Mode amplitudes of the atoms' deviation from the exact representation at s = 0, 1000, 3000, 6000 against the exact solution of the linear pairs (started from the actual Θ*(q) "
              "and the actual prior); the fit error and the energy.", "",
              "| η | mode | amplitude at s = 0, 1000, 3000, 6000 | linear prediction at the same times | ratio at s = 1000, 3000, 6000 | log-slope on [3000, 6000] / λ_j^− | fit error at s = 0, 1000, 3000, 6000 | 𝒥 at s = 6000 |",
              "|---|---|---|---|---|---|---|---|"]
    low = []                                                               # amplitudes below FLOOR at s ≥ 3000 (the remainder's level)
    for c in C:
        eta = c['eta']
        A = c['y'] @ V; Zm = c['z'] @ V
        late = np.abs(A[c['s'] >= 3000.0]); low += list(late[late < FLOOR])
        idx = [int(np.argmin(np.abs(c['s'] - t))) for t in (0.0, 1000.0, 3000.0, 6000.0)]
        lam = AP.pair_rates(h, SIGMA0, eta)[0]
        for j in (0, 1, 2, 3, len(h) - 1):
            pred_y, _ = AP.pair_solution(h[j], SIGMA0, eta, A[0, j], Zm[0, j], c['s'][idx])
            ratio = ", ".join(f"{A[i, j] / p_:.3f}" if min(abs(A[i, j]), abs(p_)) >= FLOOR else "—" for i, p_ in zip(idx[1:], pred_y[1:]))
            slope = (-(np.log(abs(A[idx[3], j])) - np.log(abs(A[idx[2], j]))) / (c['s'][idx[3]] - c['s'][idx[2]]) / lam[j]
                     if min(abs(A[idx[2], j]), abs(A[idx[3], j])) >= FLOOR else float('nan'))
            lines.append(f"| {'∞' if np.isinf(eta) else f'{eta:.5f}'} | {j + 1} | " + ", ".join(f"{A[i, j]:+.3e}" for i in idx) + " | " + ", ".join(f"{x:+.3e}" for x in pred_y)
                         + f" | {ratio} | " + (f"{slope:.3f}" if np.isfinite(slope) else "—")
                         + " | " + (", ".join(f"{c['fit'][i]:.3e}" for i in idx) if j == 0 else "") + " | " + (f"{c['J'][-1]:.3e}" if j == 0 else "") + " |")
    lines += ["", "Ratio: the amplitude over the linear prediction (constant in time when the rate is the predicted one; its offset from 1 is set in the nonlinear first phase, "
              f"the shift being finite). Log-slope: the decay rate of the amplitude on [3000, 6000] over λ_j^−. Both only where the amplitudes are above {FLOOR:.0e}: below it the "
              f"faster modes stop at the level of the second-order remainder driven by the slow modes (the amplitudes below {FLOOR:.0e} at s ≥ 3000, every mode and η: "
              f"{min(low):.1e}–{max(low):.1e})."]
    # D
    lines += ["", "## D. m-particle systems from the warm clouds (shifted prior)", "",
              f"Sizes m = {', '.join(map(str, MS))} and M = {M_REF}, {N_SEEDS} seeds, {S_D:g} slow units, plan every {PLAN_EVERY} steps. Means over the seeds (min–max). ‖ȳ‖: the whitened weighted norm of "
              "the class means' deviation from the exact representation (classes by the plan); W²: to the exact representation's measure. Gated: the prior is frozen while mass is in transit "
              "(some particle assigned to an atom that is not its nearest), so that the class barycentres it moves to contain no particle of another cluster ('transit'), and in addition, the first time, "
              "until the fixed-prior drift is below 10⁻⁵ ('settled': the adaptation starts from the fixed-prior stationary point); t_open is the first time it adapts. "
              "In parentheses: the K₀-atom run of Part C at s (ungated) or at s − t_open (gated, mean t_open).", "",
              "| η | gated | m | t_open (seeds) | fit error at s = 0, 500, 3000 | W²(μ_Θ, μ_exact) at s = 0, 500, 3000 | ‖ȳ‖ at s = 500, 3000 | prior's distance at s = 3000 | final counts |",
              "|---|---|---|---|---|---|---|---|---|"]
    Cmap = {c['eta']: c for c in C}
    combos = [(0.0, False), (SIGMA0, False), (np.inf, False), (SIGMA0, True), (np.inf, True), (SIGMA0, 'settled'), (np.inf, 'settled')]
    for eta, gate in combos:
        for m in list(MS) + [M_REF]:
            rs = [d for d in D if d['eta'] == eta and d.get('gate', False) == gate]
            if not rs:
                continue
            recs = [d['rec'][m] for d in rs]
            at = lambda key, t: np.array([r[key][int(np.argmin(np.abs(r['s'] - t)))] for r in recs])
            fit = [at('fit', t) for t in (0.0, 500.0, 3000.0)]
            w2 = [at('W2ex', t) for t in (0.0, 500.0, 3000.0)]
            yb = [np.array([np.linalg.norm(v) for v in at('ybar', t)]) for t in (500.0, 3000.0)]
            tops = [d.get('t_open', {}).get(m, 0.0) for d in rs] if gate else [0.0 for d in rs]
            ref = ""
            if eta in Cmap:
                c = Cmap[eta]; sh = float(np.mean(tops))
                ref = " (" + ", ".join(f"{np.linalg.norm(c['y'][int(np.argmin(np.abs(c['s'] - max(t - sh, 0.0))))]):.3e}" for t in (500.0, 3000.0)) + ")"
            elif eta == 0.0:
                ref = f" (Θ*(q): {np.linalg.norm(lin['ya']):.3e})"
            zd = at('zdist', 3000.0)
            cnt = "; ".join(str(tuple(int(x) for x in r['counts'][-1])) for r in recs)
            to = ", ".join(f"{x:g}" for x in tops) if gate else "—"
            gname = {False: 'no', True: 'transit', 'settled': 'settled'}[gate]
            lines.append(f"| {('∞' if np.isinf(eta) else f'{eta:.5f}')} | {gname} | {m} | {to} | " + ", ".join(f"{x.mean():.3e}" for x in fit) + " | " + ", ".join(f"{x.mean():.3e}" for x in w2)
                         + " | " + ", ".join(f"{x.mean():.3e} ({x.min():.2e}–{x.max():.2e})" for x in yb) + ref + f" | {zd.mean():.3e} | {cnt} |")
    # the settled gate against the K₀-atom run of Part C, per seed (the K₀-atom records interpolated log-linearly at 3000 − t_open)
    sl = []
    for eta in (SIGMA0, np.inf):
        rs = [d for d in D if d['eta'] == eta and d.get('gate', False) == 'settled']
        if not rs or eta not in Cmap:
            continue
        c = Cmap[eta]
        ly = np.log([np.linalg.norm(v) for v in c['y']]); lf = np.log(c['fit'])
        for m in list(MS) + [M_REF]:
            ry, rf = [], []
            for d in rs:
                r = d['rec'][m]; t_o = d['t_open'][m]
                i3 = int(np.argmin(np.abs(r['s'] - 3000.0)))
                ry.append(np.linalg.norm(r['ybar'][i3]) / np.exp(np.interp(r['s'][i3] - t_o, c['s'], ly)))
                rf.append(r['fit'][i3] / np.exp(np.interp(r['s'][i3] - t_o, c['s'], lf)))
            sl.append(f"η = {'∞' if np.isinf(eta) else f'{eta:.5f}'}, m = {m}: ‖ȳ‖ {min(ry):.3f}–{max(ry):.3f}, fit error {min(rf):.3f}–{max(rf):.3f}")
    if sl:
        lines += ["", "The settled gate against the K₀-atom run of Part C: per seed, ‖ȳ‖ and the fit error at s = 3000 over the K₀-atom run's at s = 3000 − t_open (min–max over the seeds): "
                  + "; ".join(sl) + "."]
    if E:
        r = E['rec']
        i3 = int(np.flatnonzero(r['s'] == 3000.0)[0])
        fa = lambda t: "a = (" + ", ".join(f"{x:+.3f}" for x in t[0]) + "), c = (" + ", ".join(f"{x:+.3f}" for x in t[1]) + "), s = (" + ", ".join(f"{x:.3f}" for x in t[2]) + ")"
        lines += ["", "## E. The ungated η = ∞ system with the largest ‖ȳ‖ at s = 3000, continued", "",
                  f"Seed {E['seed']}, m = {E['m']} (‖ȳ‖ = {E['ybar3000']:.3f} at s = 3000), re-run from its seed (fit error at s = 3000: {r['fit'][i3]:.6e}, stored {E['stored_fit3000']:.6e}) "
                  f"and continued to s = {S_E:g}.", "",
                  "| s | fit error | 𝒥 | largest spread of the class centres | largest drift |", "|---|---|---|---|---|"]
        for k in range(len(r['s'])):
            if r['s'][k] <= 3000.0 or (r['s'][k] - 3000.0) % 4000.0 == 0.0 or k == len(r['s']) - 1:
                lines.append(f"| {r['s'][k]:.0f} | {r['fit'][k]:.4e} | {r['J'][k]:.6e} | {r['spread'][k]:.4f} | {r['drift'][k]:.2e} |")
        lines += ["", f"At s = {S_E:g} the classes are collapsed onto the atoms {fa(E['atoms'])} (exact representation: {fa(E['exact'])}); the unregularised weighted loss "
                  f"there is {E['L_atoms']:.4e}, the eigenvalues of its Hessian (whitened, weighted) are " + ", ".join(f"{x:.2e}" for x in E['hess_ev'][:3]) + " (smallest three) — a plateau, not "
                  f"a minimum — and a BFGS minimisation of the unregularised loss from those atoms ends at {E['bfgs_L']:.1e} with {fa(E['bfgs_atoms'])}."]
    lines += ["", "Numbers are produced by `experiments/exp_adaptive_prior.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'adaptive_prior.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'adaptive_prior_raw.npz')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        report(z['lin'].item(), list(z['C']), list(z['D']), z['E'].item() if 'E' in z.files else None, t0)
        return
    if '--only-E' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        lin, C, D = z['lin'].item(), list(z['C']), list(z['D'])
        E = trapped_job(D)
        np.savez(raw, lin=np.array(lin, dtype=object), C=np.array(C, dtype=object), D=np.array(D, dtype=object), E=np.array(E, dtype=object))
        report(lin, C, D, E, t0)
        return
    if '--only-gated' in sys.argv:
        which = sys.argv[sys.argv.index('--only-gated') + 1] if len(sys.argv) > sys.argv.index('--only-gated') + 1 else 'transit'
        gate = True if which == 'transit' else which
        z = np.load(raw, allow_pickle=True)
        lin, C, D = z['lin'].item(), list(z['C']), [d for d in z['D'] if d.get('gate', False) != gate]
        with mp.Pool(N_WORKERS) as pool:
            D += pool.map(m_job, [(seed, eta, gate) for seed in range(N_SEEDS) for eta in ETAS_D if eta > 0])
        E = z['E'].item() if 'E' in z.files else None
        save = dict(lin=np.array(lin, dtype=object), C=np.array(C, dtype=object), D=np.array(D, dtype=object))
        if E is not None:
            save['E'] = np.array(E, dtype=object)
        np.savez(raw, **save)
        report(lin, C, D, E, t0)
        return
    lin = linear_parts()
    with mp.Pool(N_WORKERS) as pool:
        C = pool.map(k0_job, list(ETAS))
        D = pool.map(m_job, [(seed, eta, False) for seed in range(N_SEEDS) for eta in ETAS_D] + [(seed, eta, g) for g in (True, 'settled') for seed in range(N_SEEDS) for eta in ETAS_D if eta > 0])
    E = trapped_job(D)
    np.savez(raw, lin=np.array(lin, dtype=object), C=np.array(C, dtype=object), D=np.array(D, dtype=object), E=np.array(E, dtype=object))
    report(lin, C, D, E, t0)


if __name__ == '__main__':
    main()
