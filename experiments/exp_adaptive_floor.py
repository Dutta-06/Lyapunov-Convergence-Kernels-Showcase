"""
Proposition D‴.24: the finite-size floor is the leak's bias at the sampled weights; an adaptive prior removes it without a mass constraint
====================================================================================================================================

The averaged slow loop on the three-bump instance (σ₀ = 0.15σ_crit unless stated), the warm clouds of Proposition D‴.12 (m = 15, 60, 240
against M = 510, the m-system the first m rows, four seeds), midpoint rule dt = 0.25, the plan (cell or OT) every 4 steps.

A. The floor as a bias (no dynamics), matching prior, σ₀ = 0.05, 0.15, 0.5 σ_crit: Proposition D‴.9's function-level sensitivity F (central
   differences of the fixed-prior nearest-atom stationary point) against the bias formula Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ (src/core/floor_bias.py), the
   floor law per unit of 1/m − 1/M, its terms per mode h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j and the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ); the predicted floors at
   m = 15, 60, 240. The admissible weights of the instance (an exact representation needs |A_k|/w_k < a_max) and, per seed and size, the
   sampled weights' margin min_k(w_k − |A_k|/a_max) and the exact representation's residual.
B. Dynamics at the matching prior, the nearest-atom leak toward the fixed prior (η = 0) and toward the adaptive prior (η = σ₀, ∞; the cell plan,
   the cell masses), 6000 slow units: ‖f̂_m − f̂_M‖² and the fit error every 50 slow units, the counts, the class atoms' whitened weighted
   distance from the exact representation of the sampled weights, and at the end W²(μ_m, μ_M) against the W² between the exact
   representations of the two systems' weights (the parameter-level floor of the amplitude–mass gauge).
C. The design table at the misplaced prior (centres +0.15), 3000 slow units: the nearest-atom and the OT leak, each toward the fixed and the
   adaptive prior (η = σ₀; the OT one with the settled gate of Theorem D‴.23): the fit error (the bias) and ‖f̂_m − f̂_M‖² (the floor) at the end.
D. The real closed loop (k_e = 25, σ = σ₀/k_e, the prior's rate η/k_e; RK4 plant with the parameters frozen over the step, Euler parameters,
   dt = 10⁻³, the cells every 100 steps), matching prior, nearest-atom leak, m = 15 and 60 at the sampled weights of seed 0: from the lift of
   the fixed prior's stationary point Θ*(w), toward the fixed prior and toward the adaptive prior (η = ∞), 3000 time units; and from the exact
   representation of w with the adaptive prior at its atoms. The rms tracking error over windows of five periods against (1/k_e)‖f − f̂‖_ρ of
   the averaged loop at s = t/k_e from the same start (`--only D` computes it and merges it).

Output: results/chaos/adaptive_floor.md and adaptive_floor_raw.npz (AF_SEEDS, AF_WORKERS, AF_SB, AF_SC, AF_TREAL, AF_OUT; `--from-raw` regenerates
the md).
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
from src.core import prior_design as PD
from src.core import adaptive_prior as AP
from src.core import floor_bias as FB
from src.core.ot_leak import _cost_matrix, transport_plan

OUT = os.environ.get('AF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
A_BUMP = np.array([b[0] for b in BUMPS])
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q = np.ones(3) / 3
MS, M_REF = (15, 60, 240), 510
JITTER = 0.02
DT, PLAN_EVERY, REC = 0.25, 4, 200                                        # records every 50 slow units
S_B, S_C = float(os.environ.get('AF_SB', 6000.0)), float(os.environ.get('AF_SC', 3000.0))
SHIFT_C = 0.15
N_SEEDS = int(os.environ.get('AF_SEEDS', 4))
N_WORKERS = int(os.environ.get('AF_WORKERS', 2))


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def representation(net, shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def setup(shift=0.0, sigma0=SIGMA0):
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    exact = representation(net, 0.0)
    prior = representation(net, shift)
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, sigma0, prior, n_quad=N_Q)
    return net, exact, prior, loop


# ---------------------------------------------------------------------------- A. the floor as a bias
def part_A():
    out = dict(sig=[])
    for frac in (0.05, 0.15, 0.5):
        net, exact, prior, loop = setup(0.0, frac * SIGMA_CRIT)
        law = FB.bias_floor_law(loop, Q, exact, frac * SIGMA_CRIT)
        S = PD.weight_sensitivity_prior(loop.nearest, Q, exact)
        F2 = PD.function_sensitivity(loop.nearest, Q, exact, S) @ np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
        numeric = float(np.trace(F2 @ law['C2'] @ F2.T)) / F2.shape[0]
        # the parameter level: the fixed prior's sensitivity S (whitened weighted) against H_L(H_L + σ₀I)⁻¹S^e, S^e = ∂_wϑ^e (the gauge's);
        # Proposition D‴.5's coupled-distance law per unit of 1/m − 1/M for both
        V2 = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
        Sy = np.column_stack([(np.sqrt(Q)[:, None] * np.einsum('kcj,j->kc', S, V2[:, c]) / np.sqrt(G)[None, :]).ravel() for c in range(2)])
        Sey = -law['dx0']
        split = law['HL'] @ np.linalg.solve(law['HL'] + frac * SIGMA_CRIT * np.eye(9), Sey)
        out['sig'].append(dict(frac=frac, sigma0=frac * SIGMA_CRIT, rel=float(np.linalg.norm(law['F_bias'] - F2) / np.linalg.norm(F2)),
                               law=law['law'], numeric=numeric, bound=law['bound'], h=law['h'], terms=law['terms'],
                               split_rel=float(np.linalg.norm(Sy - split) / np.linalg.norm(Sy)),
                               D_fixed=float(np.trace(Sy @ law['C2'] @ Sy.T)), D_gauge=float(np.trace(Sey @ law['C2'] @ Sey.T))))
    # admissibility of the sampled weights and the exact representations
    net, exact, prior, loop = setup(0.0)
    adm = {}
    for seed in range(N_SEEDS):
        clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
        for m in MS + (M_REF,):
            counts = np.bincount(np.argmin(_cost_matrix(clouds[:m], exact, G), axis=1), minlength=3)
            w = counts / m
            margin = float(np.min(w - np.abs(A_BUMP) / net.a_max))
            if margin > 0:
                te, res, fit = FB.exact_weighted_atoms(loop, w, exact)
            else:
                res, fit = float('nan'), float('nan')
            adm[(seed, m)] = dict(counts=counts, margin=margin, res=res, fit=fit)
    out['adm'] = adm
    return out


# ---------------------------------------------------------------------------- B, C. dynamics
def class_deviation(Th, Z, w_exact_atoms, w):
    """The whitened weighted norm of the cell means' deviation from the exact representation of the weights w (cells of Z)."""
    lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
    means = np.array([Th[lab == k].mean(axis=0) if np.any(lab == k) else w_exact_atoms[k] for k in range(3)])
    return float(np.linalg.norm(AP.whitened(means - w_exact_atoms, w, G)))


def run_job(args):
    part, seed, plan, eta, gate, shift, S_end = args
    t0 = time.time()
    net, exact, prior, loop = setup(shift)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    sysm = {m: clouds[:m].copy() for m in MS + (M_REF,)}
    out = dict(part=part, seed=seed, plan=plan, eta=eta, gate=gate, shift=shift, rec={}, final={}, t_open={})
    for m, Th0 in sysm.items():
        counts0 = np.bincount(np.argmin(_cost_matrix(Th0, exact, G), axis=1), minlength=3)
        w0 = counts0 / m
        te = FB.exact_weighted_atoms(loop, w0, exact)[0] if np.min(w0 - np.abs(A_BUMP) / net.a_max) > 0 else exact
        apl = AP.AdaptivePriorLoop(loop, eta, plan=plan)
        rec = dict(s=[], fit=[], dev=[], counts=[])

        def cb(s, Th, Z, adapting, rec=rec, te=te, w0=w0):
            rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['dev'].append(class_deviation(Th, Z, te, w0))
            rec['counts'].append(np.bincount(np.argmin(_cost_matrix(Th, Z, G), axis=1), minlength=3))
        Th, Z, t_open = apl.run(Th0, prior.copy(), S_end, DT, plan_every=PLAN_EVERY, rec_every=REC, callback=cb, gate=gate)
        out['rec'][m] = {k: np.array(v) for k, v in rec.items()}
        out['final'][m] = Th.astype(np.float32)
        out['t_open'][m] = t_open
    # the function-level distance to the reference, along the records: rerun is not needed — the f̂'s at the records are not stored, so the
    # distance is computed from the stored finals at the end and, along the way, from the fits (both systems converge to f when exact)
    fM = loop.f_hat(np.asarray(out['final'][M_REF], float))
    out['h2_end'] = {m: float(np.mean((loop.f_hat(np.asarray(out['final'][m], float)) - fM) ** 2)) for m in MS}
    print(f"[{time.time() - t0:.0f}s] {part} seed {seed} {plan} eta {eta:g} gate {gate}", flush=True)
    return out


def run_job_h2(args):
    """As run_job, but recording ‖f̂_m − f̂_M‖² along the way (all sizes integrated in lockstep)."""
    part, seed, plan, eta, gate, shift, S_end = args
    t0 = time.time()
    net, exact, prior, loop = setup(shift)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    sizes = MS + (M_REF,)
    Th = {m: clouds[:m].copy() for m in sizes}
    Z = {m: prior.copy() for m in sizes}
    apl = {m: AP.AdaptivePriorLoop(loop, eta, plan=plan) for m in sizes}
    te, w0 = {}, {}
    for m in sizes:
        c0 = np.bincount(np.argmin(_cost_matrix(Th[m], exact, G), axis=1), minlength=3); w0[m] = c0 / m
        te[m] = FB.exact_weighted_atoms(loop, w0[m], exact)[0] if np.min(w0[m] - np.abs(A_BUMP) / net.a_max) > 0 else exact
    rec = dict(s=[], h2={m: [] for m in MS}, fit={m: [] for m in sizes}, dev={m: [] for m in sizes}, counts={m: [] for m in sizes})
    gamma = {}
    n = int(round(S_end / DT))
    for k in range(n + 1):
        if k % PLAN_EVERY == 0:
            for m in sizes:
                gamma[m], _ = apl[m].plan(Th[m], Z[m])
                if np.isinf(eta):
                    Z[m] = apl[m].bary(Th[m], gamma[m], Z[m]); gamma[m], _ = apl[m].plan(Th[m], Z[m])
        if k % REC == 0:
            rec['s'].append(k * DT)
            fM = loop.f_hat(Th[M_REF])
            for m in sizes:
                rec['fit'][m].append(loop.fit_error(Th[m])); rec['dev'][m].append(class_deviation(Th[m], Z[m], te[m], w0[m]))
                rec['counts'][m].append(np.bincount(np.argmin(_cost_matrix(Th[m], Z[m], G), axis=1), minlength=3))
            for m in MS:
                rec['h2'][m].append(float(np.mean((loop.f_hat(Th[m]) - fM) ** 2)))
        if k == n:
            break
        for m in sizes:
            a = apl[m]
            k1T, k1Z = a.drift(Th[m], Z[m], gamma[m])
            k2T, k2Z = a.drift(Th[m] + 0.5 * DT * k1T, Z[m] + 0.5 * DT * k1Z, gamma[m])
            Th[m] = Th[m] + DT * k2T; Z[m] = Z[m] + DT * k2Z
            if np.isinf(eta):
                Z[m] = a.bary(Th[m], gamma[m], Z[m])
    out = dict(part=part, seed=seed, plan=plan, eta=eta, gate=gate, shift=shift, s=np.array(rec['s']),
               h2={m: np.array(v) for m, v in rec['h2'].items()}, fit={m: np.array(v) for m, v in rec['fit'].items()},
               dev={m: np.array(v) for m, v in rec['dev'].items()}, counts={m: np.array(v) for m, v in rec['counts'].items()},
               final={m: Th[m].astype(np.float32) for m in sizes}, w0=w0)
    # the parameter level at the end: W²(μ_m, μ_M) against W² between the exact representations of the two weight vectors
    out['W2_end'] = {m: OT.W2_empirical(Th[m], Th[M_REF], G) for m in MS}
    out['W2_exact'] = {m: float(np.sum(transport_plan(_cost_matrix(te[m], te[M_REF], G), w0[m], w0[M_REF]) * _cost_matrix(te[m], te[M_REF], G))) for m in MS}
    print(f"[{time.time() - t0:.0f}s] {part} seed {seed} {plan} eta {eta:g}", flush=True)
    return out


# ---------------------------------------------------------------------------- D. the real closed loop
KE, DT_REAL, T_REAL = 25.0, 1e-3, float(os.environ.get('AF_TREAL', 3000.0))


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t)])


def real_loop(net, Th0, Z0, eta, T, adaptive):
    """The real closed loop with the nearest-atom leak toward Z (fixed, or adaptive at the rate η/k_e; η = ∞: the cell means): e² every 0.1 time
    units, the final state."""
    x = x_ref(0.0).copy(); Th = np.array(Th0, float); Z = np.array(Z0, float)
    sig = SIGMA0 / KE
    n = int(round(T / DT_REAL)); t = 0.0
    ts, e2 = [], []
    for i in range(n):
        if i % 100 == 0:
            lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
            if adaptive and np.isinf(eta):
                for k in range(Z.shape[0]):
                    if np.any(lab == k):
                        Z[k] = Th[lab == k].mean(axis=0)
            tgt = Z[lab]
            ts.append(t); e2.append(float((x - x_ref(t))[0] ** 2))
        e = x - x_ref(t)
        xdot = lambda y, tau: f_vec(y) + x_ref_dot(tau) - KE * (y - x_ref(tau)) - net.f_hat(y, Th)
        k1 = xdot(x, t); k2 = xdot(x + 0.5 * DT_REAL * k1, t + 0.5 * DT_REAL); k3 = xdot(x + 0.5 * DT_REAL * k2, t + 0.5 * DT_REAL)
        k4 = xdot(x + DT_REAL * k3, t + DT_REAL)
        Thn = Th + DT_REAL * (G[None, :] * net.grad_g(x, e, Th) - sig * (Th - tgt))
        if adaptive and np.isfinite(eta) and eta > 0:
            means = np.array([Th[lab == k].mean(axis=0) if np.any(lab == k) else Z[k] for k in range(Z.shape[0])])
            Z = Z - DT_REAL * (eta / KE) * (Z - means)
        x = x + DT_REAL / 6 * (k1 + 2 * k2 + 2 * k3 + k4); Th = Thn
        t = (i + 1) * DT_REAL
    return np.array(ts), np.array(e2), Th, Z


def part_D_job(args):
    m, case = args
    t0 = time.time()
    net, exact, prior, loop = setup(0.0)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(0), JITTER, G)
    counts = np.bincount(np.argmin(_cost_matrix(clouds[:m], exact, G), axis=1), minlength=3); w = counts / m
    Tw, res = PD.weighted_stationary_atoms_prior(loop.nearest, w, exact, w0=Q, n_cont=4)    # the fixed prior's stationary point at w
    te, rese, fite = FB.exact_weighted_atoms(loop, w, exact)                                  # the exact representation at w
    if case == 'exact':
        Th0, Z0, adaptive, eta = np.repeat(te, counts, axis=0), te.copy(), True, np.inf
    else:
        Th0, Z0, adaptive, eta = np.repeat(Tw, counts, axis=0), exact.copy(), case == 'adaptive', (np.inf if case == 'adaptive' else 0.0)
    ts, e2, Th, Z = real_loop(net, Th0, Z0, eta, T_REAL, adaptive)
    # the averaged loop from the same start, at s = t/k_e
    apl = AP.AdaptivePriorLoop(loop, eta if adaptive else 0.0, plan='nearest')
    rec = dict(s=[], fit=[])
    apl.run(Th0, Z0, T_REAL / KE, 0.05, plan_every=1, rec_every=20, callback=lambda s_, T_, Z_, ad: (rec['s'].append(s_), rec['fit'].append(loop.fit_error(T_))))
    normF = float(np.sqrt(np.mean(loop.F ** 2)))
    print(f"[{time.time() - t0:.0f}s] D m {m} {case}", flush=True)
    return dict(m=m, case=case, counts=counts, ts=ts, e2=e2, s_avg=np.array(rec['s']), fit_avg=np.array(rec['fit']), normF=normF,
                res_fixed=res, res_exact=rese, fit_exact=fite, final=Th.astype(np.float32))


def post_B(B):
    """Part B's parameter-level columns from the stored finals: the exact representation of each system's sampled weights (none on the
    boundary of the admissible set), W² between the exact representations of the m-system's and the reference's weights, and the end
    distance of the class means (classes by the exact atoms' cells) from the exact representation."""
    net, exact, prior, loop = setup(0.0)
    cache = {}
    for r in B:
        te = {}
        for m in MS + (M_REF,):
            w = np.asarray(r['w0'][m], float)
            key = tuple(np.round(w * m).astype(int))
            if key not in cache:
                cache[key] = FB.exact_weighted_atoms(loop, w, exact) if np.min(w - np.abs(A_BUMP) / net.a_max) > 0 else None
            te[m] = cache[key]
        r['W2_exact'], r['dev_end'] = {}, {}
        for m in MS:
            ok = te[m] is not None and te[M_REF] is not None and te[m][2] < 1e-10 and te[M_REF][2] < 1e-10
            if ok:
                C = _cost_matrix(te[m][0], te[M_REF][0], G)
                r['W2_exact'][m] = float(np.sum(transport_plan(C, r['w0'][m], r['w0'][M_REF]) * C))
                r['dev_end'][m] = class_deviation(np.asarray(r['final'][m], float), exact, te[m][0], np.asarray(r['w0'][m], float))
            else:
                r['W2_exact'][m] = float('nan'); r['dev_end'][m] = float('nan')
    return B


# ---------------------------------------------------------------------------- report
def report(A, B, C, D=None, t0=None):
    lines = ["# Proposition D‴.24: the finite-size floor is the leak's bias; an adaptive prior removes it without a mass constraint (`exp_adaptive_floor.py`)", "",
             f"The three-bump instance (σ_crit = {SIGMA_CRIT:.5f}); the warm clouds of Proposition D‴.12 (width {JITTER:g}), m = {', '.join(map(str, MS))} against M = {M_REF}, {N_SEEDS} seeds; "
             f"midpoint rule dt = {DT:g}, the plan every {PLAN_EVERY} steps.", "",
             "## A. The floor as the weight-derivative of the leak's bias (matching prior, no dynamics)", "",
             "F: Proposition D‴.9's function-level sensitivity along e₁ − e₃, e₂ − e₃ (central differences of the fixed-prior nearest-atom stationary point); F_bias = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰. "
             "The law per unit of 1/m − 1/M, tr(FC̃Fᵀ)/N, from F and as the sum of the mode terms h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j; the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ).", "",
             "| σ₀/σ_crit | ‖F − F_bias‖/‖F‖ | law from F | sum of the mode terms | bound | law/bound | floor (1/m − 1/M)·law at m = 15, 60, 240 |",
             "|---|---|---|---|---|---|---|"]
    for r in A['sig']:
        fl = ", ".join(f"{(1.0 / m - 1.0 / M_REF) * r['numeric']:.2e}" for m in MS)
        lines.append(f"| {r['frac']:g} | {r['rel']:.1e} | {r['numeric']:.4e} | {r['law']:.4e} | {r['bound']:.4e} | {r['law'] / r['bound']:.3f} | {fl} |")
    lines += ["", "The parameter level: the fixed prior's weight sensitivity S (whitened weighted) against H_L(H_L + σ₀I)⁻¹S^e with S^e = ∂_wϑ^e (per mode the fraction "
              "h_j/(h_j + σ₀) of the gauge's displacement, the rest σ₀/(h_j + σ₀) being the bias), and Proposition D‴.5's coupled-distance law per unit of 1/m − 1/M, "
              "Σ tr(S̃C̃S̃ᵀ), for the fixed prior and for the adaptive one (S^e: the whole displacement in the parameters).", "",
              "| σ₀/σ_crit | ‖S − H_L(H_L + σ₀I)⁻¹S^e‖/‖S‖ | coupled-distance law, fixed prior | adaptive prior | ratio |", "|---|---|---|---|---|"]
    for r in A['sig']:
        lines.append(f"| {r['frac']:g} | {r['split_rel']:.1e} | {r['D_fixed']:.4e} | {r['D_gauge']:.4e} | {r['D_gauge'] / r['D_fixed']:.3f} |")
    lines += ["", "Mode by mode (whitened weighted Gauss–Newton eigenvalues h_j; the share of each term in the law):", "",
              "| σ₀/σ_crit | " + " | ".join(f"mode {j + 1}" for j in range(9)) + " |", "|---|" + "---|" * 9]
    for r in A['sig']:
        lines.append(f"| {r['frac']:g}: h_j/σ₀ | " + " | ".join(f"{x:.2f}" for x in r['h'] / r['sigma0']) + " |")
        lines.append(f"| {r['frac']:g}: share | " + " | ".join(f"{x:.3f}" for x in r['terms'] / r['terms'].sum()) + " |")
    lines += ["", f"Admissible weights: an exact representation at the weights w needs |A_k|/w_k < a_max = 4, i.e. w > ({', '.join(f'{x:.3f}' for x in np.abs(A_BUMP) / 4.0)}). "
              "Per seed and size: the sampled counts, the margin min_k(w_k − |A_k|/a_max), and the exact representation's gradient residual and fit error.", "",
              "| m | " + " | ".join(f"seed {s}" for s in range(N_SEEDS)) + " |", "|---|" + "---|" * N_SEEDS]
    for m in MS + (M_REF,):
        cells = []
        for s in range(N_SEEDS):
            a = A['adm'][(s, m)]
            cells.append(f"{tuple(int(x) for x in a['counts'])}: {a['margin']:+.4f}" + (f", {a['res']:.0e}, {a['fit']:.0e}" if np.isfinite(a['res']) else " (none)"))
        lines.append(f"| {m} | " + " | ".join(cells) + " |")
    # B
    lines += ["", f"## B. The nearest-atom leak at the matching prior: fixed (η = 0) against adaptive priors ({S_B:g} slow units)", "",
              "‖f̂_m − f̂_M‖² (the function-level floor) and the fit error at s = 600, 3000, 6000, per seed; W²(μ_m, μ_M) at the end against the W² between the exact "
              "representations of the two systems' sampled weights (the parameter-level floor); the class means' whitened weighted distance from the exact representation of "
              "the sampled weights at the end; the counts (all conserved). —: no exact representation at those weights (on the boundary of the admissible set).", "",
              "| η | m | ‖f̂_m − f̂_M‖² at s = 600, 3000, 6000 (per seed) | fit error at 6000 (per seed) | W²(μ_m, μ_M) at the end / exact-representation W² (per seed) | distance from ϑ^e(w) at the end | counts conserved |",
              "|---|---|---|---|---|---|---|"]
    for eta in (0.0, SIGMA0, np.inf):
        rs = sorted([r for r in B if r['eta'] == eta], key=lambda r: r['seed'])
        for m in MS:
            h2 = "; ".join(", ".join(f"{r['h2'][m][int(np.argmin(np.abs(r['s'] - t)))]:.1e}" for t in (600.0, 3000.0, S_B)) for r in rs)
            fit = "; ".join(f"{r['fit'][m][-1]:.1e}" for r in rs)
            w2 = "; ".join(f"{r['W2_end'][m]:.2e} / " + (f"{r['W2_exact'][m]:.2e}" if np.isfinite(r['W2_exact'][m]) else "—") for r in rs)
            dev = "; ".join((f"{r['dev_end'][m]:.1e}" if np.isfinite(r['dev_end'][m]) else "—") if 'dev_end' in r else f"{r['dev'][m][-1]:.1e}" for r in rs)
            cons = all(np.all(r['counts'][m] == r['counts'][m][0]) for r in rs)
            lines.append(f"| {'∞' if np.isinf(eta) else f'{eta:.5f}'} | {m} | {h2} | {fit} | {w2} | {dev} | {'yes' if cons else 'no'} |")
    # C
    lines += ["", f"## C. The design table at the misplaced prior (centres +{SHIFT_C:g}), {S_C:g} slow units", "",
              "The fit error (the bias; 0 at an exact representation) and ‖f̂_m − f̂_M‖² (the floor) at the end, means over the seeds (min–max); the OT leak toward the adaptive "
              "prior with the settled gate of Theorem D‴.23.", "",
              "| leak | prior | m | fit error at the end | ‖f̂_m − f̂_M‖² at the end |", "|---|---|---|---|---|"]
    for plan, eta, gate in (('nearest', 0.0, False), ('nearest', SIGMA0, False), ('ot', 0.0, False), ('ot', SIGMA0, 'settled')):
        rs = [r for r in C if r['plan'] == plan and r['eta'] == eta and r['gate'] == gate]
        for m in MS:
            fit = np.array([r['rec'][m]['fit'][-1] for r in rs]); h2 = np.array([r['h2_end'][m] for r in rs])
            lines.append(f"| {'OT' if plan == 'ot' else 'nearest atom'} | {'adaptive (η = σ₀' + (', settled gate)' if gate else ')') if eta > 0 else 'fixed'} | {m} | "
                         f"{fit.mean():.2e} ({fit.min():.1e}–{fit.max():.1e}) | {h2.mean():.2e} ({h2.min():.1e}–{h2.max():.1e}) |")
    if D:
        win = lambda ts, e2, t: float(np.sqrt(np.mean(e2[(ts > t - 5 * PER) & (ts <= t)])))
        times = (5 * PER, 1000.0, 2000.0, T_REAL)
        lines += ["", f"## D. The real closed loop (k_e = {KE:g}, dt = {DT_REAL:g}, {T_REAL:g} time units), matching prior, nearest-atom leak, seed 0's sampled weights", "",
                  "The rms tracking error over the five periods ending at t = " + ", ".join(f"{x:.0f}" for x in times) + ", against (1/k_e)‖f − f̂‖_ρ of the averaged loop "
                  "from the same start at s = t/k_e (the slaved error of Theorem J at leading order); from the exact representation of w with the adaptive prior at its atoms "
                  "the tracking error stays at the integrator's level.", "",
                  "| m | counts | start, prior | e_rms at those times | (1/k_e)‖f − f̂‖_ρ averaged at s = t/k_e | ratio at the end |", "|---|---|---|---|---|---|"]
        name = {'fixed': 'Θ*(w), fixed prior', 'adaptive': 'Θ*(w), adaptive prior (η = ∞)', 'exact': 'ϑ^e(w), adaptive prior at ϑ^e(w)'}
        for d in sorted(D, key=lambda d: (d['m'], ['fixed', 'adaptive', 'exact'].index(d['case']))):
            er = [win(d['ts'], d['e2'], t) for t in times]
            pr = [float(np.interp(t / KE, d['s_avg'], d['fit_avg'])) * d['normF'] / KE for t in times]
            lines.append(f"| {d['m']} | {tuple(int(x) for x in d['counts'])} | {name[d['case']]} | " + ", ".join(f"{x:.3e}" for x in er) + " | "
                         + ", ".join(f"{x:.3e}" for x in pr) + f" | {er[-1] / pr[-1]:.3f} |" if pr[-1] > 1e-12 else
                         f"| {d['m']} | {tuple(int(x) for x in d['counts'])} | {name[d['case']]} | " + ", ".join(f"{x:.1e}" for x in er) + " | " + ", ".join(f"{x:.1e}" for x in pr) + " | — |")
    lines += ["", "Numbers are produced by `experiments/exp_adaptive_floor.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'adaptive_floor.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'adaptive_floor_raw.npz')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        report(z['A'].item(), list(z['B']), list(z['C']), list(z['D']) if 'D' in z.files else None, t0)
        return
    if '--only' in sys.argv and sys.argv[sys.argv.index('--only') + 1] == 'A':          # Part A and Part B's post-processing
        with np.load(raw, allow_pickle=True) as z:                                   # read everything before the file is rewritten
            save = {k: z[k] for k in z.files}
        A = part_A(); B = post_B(list(save['B']))
        save['A'] = np.array(A, dtype=object); save['B'] = np.array(B, dtype=object)
        np.savez(raw, **save)
        report(A, B, list(save['C']), list(save['D']) if 'D' in save else None, t0)
        return
    if '--only' in sys.argv and sys.argv[sys.argv.index('--only') + 1] == 'D':
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(N_WORKERS) as pool:
            D = pool.map(part_D_job, [(m, case) for m in (15, 60) for case in ('fixed', 'adaptive', 'exact')], chunksize=1)
        save['D'] = np.array(D, dtype=object)
        np.savez(raw, **save)
        report(save['A'].item(), list(save['B']), list(save['C']), D, t0)
        return
    A = part_A()
    print("[A done]", flush=True)
    jobsB = [('B', seed, 'nearest', eta, False, 0.0, S_B) for eta in (0.0, SIGMA0, np.inf) for seed in range(N_SEEDS)]
    jobsC = [('C', seed, plan, eta, gate, SHIFT_C, S_C) for plan, eta, gate in (('nearest', 0.0, False), ('nearest', SIGMA0, False), ('ot', 0.0, False), ('ot', SIGMA0, 'settled'))
             for seed in range(N_SEEDS)]
    with mp.Pool(N_WORKERS) as pool:
        B = post_B(pool.map(run_job_h2, jobsB, chunksize=1))
        C = pool.map(run_job, jobsC, chunksize=1)
    np.savez(raw, A=np.array(A, dtype=object), B=np.array(B, dtype=object), C=np.array(C, dtype=object))
    report(A, B, C, None, t0)


if __name__ == '__main__':
    main()
