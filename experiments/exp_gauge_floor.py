"""
Proposition D‴.28: the adaptive prior's floor removal is a gauge property — it survives an approximation error
==========================================================================================================

The three-bump instance with a fourth bump that three kernels cannot represent, f_ε = f + ε·exp(−(x − 0.8)²/(2·0.3²)), ε ∈ EPS; K₀ = 3, the
masses q = 1/3, σ₀ = 0.15σ_crit. ϑ^b(q): the best weighted K₀-atom approximant (a nondegenerate critical point of the unregularised weighted
loss L_q, continued in ε from the exact representation of f_0). The amplitude–mass gauge G_w (α_k ↦ a⁻¹(q_ka(α_k)/w_k)) leaves f̂ unchanged,
so L_w∘G_w = L_q.

A. No dynamics, per ε: the approximant's fit error; the whitened unregularised Hessian Ã = H_L − R̃ (its smallest eigenvalue, ‖R̃‖/‖H_L‖);
   the gauge set W_g = {w_k > q_k|a^b_k|/a_max}; the gauge check — Newton on the unregularised weighted gradient at sampled weights from
   G_wϑ^b(q): its distance to G_wϑ^b(q) and its f̂ against f̂_b(q); the function-level sensitivity of the adaptive prior's stationary point
   (central differences in w: ≈ 0) and of the fixed matching prior's (ϑ⁰ = ϑ^b(q)) against σ₀Ĵ(Ã + σ₀I)⁻¹g̃ and against D‴.24's form with H_L
   in place of Ã; the fixed prior's floor law tr(FC̃Fᵀ)/N; the residual curvature at ϑ^b(q) — its amplitude row (zero at stationarity) and
   λ⁺_k := λ_max(R^G_k) — and the adaptive prior's splitting margins σ₀ − (q_k/w_k)λ⁺_k against the margins computed at G_wϑ^b(q), over all
   count vectors of m = 15, 60 with the probability of a negative margin.
B. Dynamics at ε = 0.1 and 0.3: the nearest-atom leak toward the fixed matching prior ϑ^b(q) (η = 0) and the adaptive prior (η = σ₀, ∞) from
   warm clouds around ϑ^b(q) (width 0.02), m = 15, 60, 240 against M = 510, four seeds, S_B slow units: ‖f̂_m − f̂_M‖² and the fit error at
   s = 600, 3000, S_B, the counts, the class means' distance from G_wϑ^b(q), and the end diameters of the classes.

Output: results/chaos/gauge_floor.md and gauge_floor_raw.npz (GF_WORKERS, GF_SB, GF_OUT; `--only-B` computes Part B, saving after
every job and skipping the jobs already saved; `--from-raw` regenerates the md).
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
from src.core.ot_leak import _cost_matrix

OUT = os.environ.get('GF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
EXTRA = (0.8, 0.3)                                                         # the fourth bump's centre and width
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q = np.ones(3) / 3
EPS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5)
EPS_B = (0.1, 0.3)
MS, M_REF = (15, 60, 240), 510
JITTER = 0.02
DT, PLAN_EVERY, REC = 0.25, 4, 200
S_B = float(os.environ.get('GF_SB', 6000.0))
N_SEEDS = 4
N_WORKERS = int(os.environ.get('GF_WORKERS', 2))
H_FD = 1e-5


def target(eps):
    def f(X):
        X = np.asarray(X, float)
        return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS) + eps * np.exp(-(X - EXTRA[0]) ** 2 / (2 * EXTRA[1] ** 2))
    return f


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def net_():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact0(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def best_approximant(eps, n_steps=None):
    """ϑ^b(q) at ε: the critical point of the unregularised weighted loss at q, continued in ε from f_0's exact representation."""
    net = net_()
    th = exact0(net)
    n_steps = n_steps or max(1, int(np.ceil(eps / 0.025)))
    for k in range(1, n_steps + 1):
        e = eps * k / n_steps
        loop = OT.OTPriorLoop(net, target(e), x_ref, PER, G, SIGMA0, th, n_quad=N_Q)
        th, res, fit = FB.exact_weighted_atoms(loop, Q, th, n_cont=1, w0=Q, fit_tol=1.0)
    return net, th, res, fit


def setup(eps):
    """The loop with the fixed matching prior ϑ⁰ = ϑ^b(q) and the approximant."""
    net, tb, res, fit = best_approximant(eps)
    loop = OT.OTPriorLoop(net, target(eps), x_ref, PER, G, SIGMA0, tb, n_quad=N_Q)
    return net, tb, loop, res, fit


def gauge(net, tb, w):
    """G_wϑ^b(q): the amplitudes q_ka^b_k/w_k (None outside W_g)."""
    a, c, s = net.physical(tb)
    an = Q[:, None] * a / np.asarray(w, float)[:, None]
    if np.any(np.abs(an) >= net.a_max):
        return None
    return net.unconstrained(an, c, s)


def fhat(loop, T, w):
    a, c, s, d, phi = loop.features(T)
    return phi @ (np.asarray(w, float) * a)


def hessians(loop, tb, h=1e-6):
    """At ϑ^b(q): the whitened weighted unregularised Hessian Ã, the Gauss–Newton block H_L, Ĵ and the whitened gauge displacement g̃ along
    e₁ − e₃, e₂ − e₃ (g̃ = D^{1/2}∂_w(G_wϑ^b) with the sign of D‴.24's ∂_wx⁰ reversed: f̂'s direct weight effect is Ĵg̃)."""
    K, p = tb.shape; d = K * p
    gradL = lambda v: (Q[:, None] * FB.unregularised_grad(loop, Q, v.reshape(K, p))).ravel()          # ∇_{ϑ_k}L_q = q_k∇V_L(ϑ_k)
    A = np.zeros((d, d))
    for j in range(d):
        E = np.zeros(d); E[j] = h
        A[:, j] = (gradL(tb.ravel() + E) - gradL(tb.ravel() - E)) / (2 * h)
    A = 0.5 * (A + A.T)
    Dh = np.kron(np.sqrt(Q), np.sqrt(1.0 / G))                              # D^{1/2}, D = diag(q_kG⁻¹)
    At = A / Dh[:, None] / Dh[None, :]                                      # whitened: D^{−1/2}AD^{−1/2}
    Jh = FB.whitened_jacobian(loop, Q, tb)
    HL = Jh.T @ Jh / loop.N
    net = loop.net
    Vd = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
    g = np.column_stack([((gauge(net, tb, Q + h * Vd[:, c]) - gauge(net, tb, Q - h * Vd[:, c])) / (2 * h)).ravel() for c in range(2)])
    gt = -Dh[:, None] * g                                                   # g̃ := −D^{1/2}∂_w(G_wϑ^b): f̂'s direct weight effect is Ĵg̃
    return dict(A=A, At=At, HL=HL, Jh=Jh, gt=gt, Vd=Vd)


def part_A_eps(eps):
    t0 = time.time()
    net, tb, loop, res, fit = setup(eps)
    Hs = hessians(loop, tb)
    At, HL, Jh, gt, Vd = Hs['At'], Hs['HL'], Hs['Jh'], Hs['gt'], Hs['Vd']
    Rt = HL - At
    ev = np.linalg.eigvalsh(At)
    a_b = net.physical(tb)[0][:, 0]
    lo = Q * np.abs(a_b) / net.a_max                                        # W_g: w_k > lo_k
    # the gauge check at sampled weights
    rng = np.random.default_rng(1)
    gchk = []
    for _ in range(6):
        while True:
            w = rng.dirichlet(30 * np.ones(3))
            if np.all(w > lo + 0.01):
                break
        tg = gauge(net, tb, w)
        tn, rn, fn = FB.exact_weighted_atoms(loop, w, tg, n_cont=1, w0=w, fit_tol=1.0)
        gchk.append(dict(w=w, dist=float(np.abs(tn - tg).max()), dfh=float(np.abs(fhat(loop, tn, w) - fhat(loop, tb, Q)).max()), res=rn))
    # sensitivities: adaptive (critical points of L_w) and fixed (Θ*(w) toward ϑ⁰ = ϑ^b(q)), central differences
    Fa, Ff = [], []
    for c in range(2):
        o = {}
        for sg in (+1, -1):
            w = Q + sg * H_FD * Vd[:, c]
            ta = FB.exact_weighted_atoms(loop, w, gauge(net, tb, w), n_cont=1, w0=w, fit_tol=1.0)[0]
            tf, rf = PD.weighted_stationary_atoms_prior(loop.nearest, w, tb)
            o[sg] = (fhat(loop, ta, w), fhat(loop, tf, w))
        Fa.append((o[1][0] - o[-1][0]) / (2 * H_FD)); Ff.append((o[1][1] - o[-1][1]) / (2 * H_FD))
    Fa, Ff = np.array(Fa).T, np.array(Ff).T
    F_gauge = SIGMA0 * Jh @ np.linalg.solve(At + SIGMA0 * np.eye(9), gt)    # σ₀Ĵ(Ã + σ₀I)⁻¹g̃
    F_gn = SIGMA0 * Jh @ np.linalg.solve(HL + SIGMA0 * np.eye(9), gt)       # D‴.24's form with H_L
    Ct = (np.diag(Q) - np.outer(Q, Q))[:2, :2]
    law = lambda F: float(np.trace(F @ Ct @ F.T)) / loop.N
    # the residual curvature at ϑ^b(q): its amplitude row and λ⁺_k; the adaptive prior's margins
    a, c_, s_, d_, phi = loop.features(tb)
    r = loop.F - phi @ (Q * a)
    Rk = loop.residual_curvature(r, tb)                                     # (K, 3, 3), per particle
    sqG = np.sqrt(G)
    RG = sqG[None, :, None] * Rk * sqG[None, None, :]
    amp_row = float(np.abs(Rk[:, 0, :]).max())
    lam = np.array([np.linalg.eigvalsh(0.5 * (RG[k] + RG[k].T))[-1] for k in range(3)])
    print(f"[{time.time() - t0:.0f}s] A eps {eps:g}", flush=True)
    return dict(eps=eps, fit=fit, res=res, ev=ev, R_rel=float(np.linalg.norm(Rt) / np.linalg.norm(HL)), lo=lo, a_b=a_b, gauge=gchk,
                Fa=float(np.linalg.norm(Fa) / np.linalg.norm(Ff)), law_fixed=law(Ff), law_gauge=law(F_gauge), law_gn=law(F_gn),
                rel_gauge=float(np.linalg.norm(F_gauge - Ff) / np.linalg.norm(Ff)), rel_gn=float(np.linalg.norm(F_gn - Ff) / np.linalg.norm(Ff)),
                amp_row=amp_row, lam=lam, R_norm=float(np.abs(Rk).max()))


def margins_job(args):
    """The adaptive prior's splitting margins over all count vectors of m: predicted σ₀ − (q_k/w_k)λ⁺_k against computed at G_wϑ^b(q)."""
    eps, m, lam = args
    from scipy.stats import multinomial
    net, tb, loop, res, fit = setup(eps)
    rows = []
    for n1 in range(1, m - 1):
        for n2 in range(1, m - n1):
            n = np.array([n1, n2, m - n1 - n2]); w = n / m
            tg = gauge(net, tb, w)
            if tg is None:
                continue
            a, c, s, d, phi = loop.features(tg)
            sep = SIGMA0 - loop.kappa_R(loop.F - phi @ (w * a), tg)
            pred = SIGMA0 - (Q / w) * lam
            rows.append(dict(n=n, p=float(multinomial.pmf(n, m, Q)), sep=sep, pred=pred))
    return (eps, m), rows


# ---------------------------------------------------------------------------- B. dynamics
def class_diameter(Th, lab, K=3):
    out = 0.0
    for k in range(K):
        idx = lab == k
        if idx.sum() >= 2:
            out = max(out, float(np.sqrt(_cost_matrix(Th[idx], Th[idx], G).max())))
    return out


def run_job(args):
    """All sizes in lockstep (midpoint rule, the cell plan every PLAN_EVERY steps): the nearest-atom leak toward ϑ^b(q), fixed or adaptive."""
    eps, seed, eta = args
    t0 = time.time()
    net, tb, loop, res, fit = setup(eps)
    clouds = CL.mixture_start(tb, M_REF, np.random.default_rng(seed), JITTER, G)
    lab0 = np.argmin(_cost_matrix(clouds, tb, G), axis=1)
    sizes = MS + (M_REF,)
    Th = {m: clouds[:m].copy() for m in sizes}
    Z = {m: tb.copy() for m in sizes}
    apl = {m: AP.AdaptivePriorLoop(loop, eta, plan='nearest') for m in sizes}
    w0 = {m: np.bincount(lab0[:m], minlength=3) / m for m in sizes}
    rec = dict(s=[], h2={m: [] for m in MS}, fit={m: [] for m in sizes}, counts={m: [] for m in sizes})
    gamma = {}
    n = int(round(S_B / DT))
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
                rec['fit'][m].append(loop.fit_error(Th[m]))
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
    out = dict(eps=eps, seed=seed, eta=eta, s=np.array(rec['s']), h2={m: np.array(v) for m, v in rec['h2'].items()},
               fit={m: np.array(v) for m, v in rec['fit'].items()}, counts={m: np.array(v) for m, v in rec['counts'].items()}, w0=w0,
               final={m: Th[m].astype(np.float32) for m in sizes}, approx_fit=fit)
    out['dist'], out['diam'] = {}, {}
    for m in sizes:
        tg = gauge(net, tb, w0[m])
        lab = lab0[:m]
        means = np.array([Th[m][lab == j].mean(axis=0) for j in range(3)])
        out['dist'][m] = float(np.linalg.norm(AP.whitened(means - tg, w0[m], G))) if tg is not None else float('nan')
        out['diam'][m] = class_diameter(Th[m], lab)
    print(f"[{time.time() - t0:.0f}s] B eps {eps:g} seed {seed} eta {eta:g}", flush=True)
    return out


# ---------------------------------------------------------------------------- report
def report(A, M, B, t0=None):
    es = lambda e: '∞' if np.isinf(e) else f'{e:.5f}'
    lines = ["# Proposition D‴.28: the adaptive prior's floor removal is a gauge property (`exp_gauge_floor.py`)", "",
             f"The three-bump instance with a fourth bump ε·exp(−(x − {EXTRA[0]:g})²/(2·{EXTRA[1]:g}²)) that K₀ = 3 kernels cannot represent; masses q = 1/3, "
             f"σ₀ = 0.15σ_crit = {SIGMA0:.5f}; ϑ^b(q) the best weighted approximant (a critical point of the unregularised loss, continued in ε).", "",
             "## A. The approximant, the gauge and the sensitivities (no dynamics)", "",
             "Per ε: the approximant's relative fit error and gradient residual; the smallest eigenvalue of the whitened unregularised Hessian Ã = H_L − R̃ and "
             "‖R̃‖/‖H_L‖; the gauge check at six weight vectors of W_g (the largest distance of Newton's critical point from G_wϑ^b(q), and of its f̂ from "
             "f̂_b(q)); the adaptive prior's function-level sensitivity relative to the fixed prior's (central differences); the fixed matching prior's "
             "sensitivity against σ₀Ĵ(Ã + σ₀I)⁻¹g̃ and against D‴.24's form with H_L; the fixed prior's floor law per unit of 1/m − 1/M.", "",
             "| ε | fit error | residual | min eig Ã | ‖R̃‖/‖H_L‖ | gauge: distance, Δf̂ | ‖F_ad‖/‖F_fixed‖ | F_fixed vs σ₀Ĵ(Ã + σ₀I)⁻¹g̃ | vs H_L form | law (fixed) |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in A:
        gd = max(x['dist'] for x in r['gauge']); gf = max(x['dfh'] for x in r['gauge'])
        lines.append(f"| {r['eps']:g} | {r['fit']:.3e} | {r['res']:.1e} | {r['ev'][0]:.3e} | {r['R_rel']:.2e} | {gd:.1e}, {gf:.1e} | {r['Fa']:.1e} | "
                     f"{r['rel_gauge']:.1e} | {r['rel_gn']:.1e} | {r['law_fixed']:.4e} |")
    lines += ["", "The residual curvature at ϑ^b(q): the largest entry of its amplitude row and column (zero at a critical point), and λ⁺_k = λ_max(G^{1/2}R_kG^{1/2}) per "
              "class; W_g's thresholds q_k|a^b_k|/a_max.", "", "| ε | amplitude row | largest entry of R | λ⁺_k | W_g: w > |", "|---|---|---|---|---|"]
    for r in A:
        lines.append(f"| {r['eps']:g} | {r['amp_row']:.1e} | {r['R_norm']:.2e} | " + ", ".join(f"{x:.2e}" for x in r['lam']) + " | (" + ", ".join(f"{x:.4f}" for x in r['lo']) + ") |")
    if M:
        lines += ["", "The adaptive prior's splitting margins at G_wϑ^b(q) over every count vector of m in W_g (with no empty class): the largest difference from the "
                  "prediction σ₀ − (q_k/w_k)λ⁺_k, the probability of a negative margin (some class can split), and the smallest margin.", "",
                  "| ε | m | count vectors | largest |computed − predicted| | P(negative margin) | smallest margin |", "|---|---|---|---|---|---|"]
        for (eps, m) in sorted(M):
            rows = M[(eps, m)]
            dmax = max(float(np.abs(x['sep'] - x['pred']).max()) for x in rows)
            pneg = sum(x['p'] for x in rows if np.min(x['sep']) <= 0)
            smin = min(float(np.min(x['sep'])) for x in rows)
            lines.append(f"| {eps:g} | {m} | {len(rows)} | {dmax:.1e} | {pneg:.4f} | {smin:+.2e} |")
    if B:
        lines += ["", f"## B. The nearest-atom leak from warm clouds around ϑ^b(q): fixed (η = 0) against adaptive priors ({S_B:g} slow units)", "",
                  "Per ε, η and m: ‖f̂_m − f̂_M‖² at s = 600, 3000, " + f"{S_B:g} (min–max over the seeds), the fit error at the end (min–max) against the approximant's, "
                  "the class means' whitened weighted distance from G_wϑ^b(q) at the end (max), the largest class diameter at the end (max), and the counts conserved.", "",
                  "| ε | η | m | ‖f̂_m − f̂_M‖² at those s | fit error at the end (approximant) | distance from G_wϑ^b | class diameter | counts conserved |",
                  "|---|---|---|---|---|---|---|---|"]
        for eps in sorted(set(r['eps'] for r in B)):
            for eta in (0.0, SIGMA0, np.inf):
                rs = [r for r in B if r['eps'] == eps and r['eta'] == eta]
                if not rs:
                    continue
                for m in MS:
                    h2 = [", ".join(f"{r['h2'][m][int(np.argmin(np.abs(r['s'] - t)))]:.1e}" for t in (600.0, 3000.0, S_B)) for r in rs]
                    h2v = np.array([[r['h2'][m][int(np.argmin(np.abs(r['s'] - t)))] for t in (600.0, 3000.0, S_B)] for r in rs])
                    fits = np.array([r['fit'][m][-1] for r in rs])
                    dist = np.array([r['dist'][m] for r in rs]); diam = np.array([r['diam'][m] for r in rs])
                    cons = all(np.all(r['counts'][m] == r['counts'][m][0]) for r in rs)
                    lines.append(f"| {eps:g} | {es(eta)} | {m} | " + ", ".join(f"{h2v[:, j].min():.1e}–{h2v[:, j].max():.1e}" for j in range(3)) + " | "
                                 f"{fits.min():.4e}–{fits.max():.4e} ({rs[0]['approx_fit']:.4e}) | {np.nanmax(dist):.1e} | {diam.max():.1e} | {'yes' if cons else 'no'} |")
        # the samples by their margin to the boundary of W_g
        Aeps = {r['eps']: r for r in A}
        near, excess = [], {m: [] for m in MS}
        for eps in sorted(set(r['eps'] for r in B)):
            lo = np.asarray(Aeps[eps]['lo'], float)
            for seed in sorted(set(r['seed'] for r in B)):
                byeta = {r['eta']: r for r in B if r['eps'] == eps and r['seed'] == seed}
                for m in MS:
                    w = np.asarray(byeta[0.0]['w0'][m], float)
                    marg = float(np.min(w - lo))
                    ad = [byeta[e] for e in (SIGMA0, np.inf) if e in byeta]
                    if marg < 0.02:
                        near.append((eps, m, seed, w, marg, [byeta[e]['fit'][m][-1] for e in (0.0, SIGMA0, np.inf) if e in byeta],
                                     [r['dist'][m] for r in ad], byeta[0.0]['approx_fit']))
                    else:
                        excess[m] += [r['fit'][m][-1] / r['approx_fit'] - 1 for r in ad]
        lines += ["", "The samples within 0.02 of the boundary of W_g or outside it (the margin min_k(w_k − q_k|a^b_k|/a_max); negative: outside, no gauge transform): "
                  "the fit error at the end under the fixed prior and the adaptive priors (η = σ₀, ∞) against the approximant's, and the adaptive class means' distance "
                  "from G_wϑ^b(q).", "",
                  "| ε | m | seed | w | margin | fit error: fixed, η = σ₀, η = ∞ (approximant) | distance from G_wϑ^b: η = σ₀, ∞ |", "|---|---|---|---|---|---|---|"]
        for eps, m, seed, w, marg, fits, dists, af in near:
            lines.append(f"| {eps:g} | {m} | {seed} | ({', '.join(f'{x:.3f}' for x in w)}) | {marg:+.4f} | " + ", ".join(f"{x:.4e}" for x in fits)
                         + f" ({af:.4e}) | " + ", ".join('—' if not np.isfinite(x) else f"{x:.2e}" for x in dists) + " |")
        if any(excess.values()):
            lines += ["", "At the other samples (margin ≥ 0.02) the adaptive priors' fit error at the end exceeds the approximant's by at most "
                      + ", ".join(f"{max(v):.1e} (m = {m}, median {np.median(v):.1e})" for m, v in excess.items() if v) + " (relative): still converging toward the gauge "
                      "transform, slowest at small m, where the samples lie nearest the boundary."]
    lines += ["", "Numbers are produced by `experiments/exp_gauge_floor.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'gauge_floor.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def save_raw(raw, **arrays):
    tmp = raw[:-4] + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, raw)


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'gauge_floor_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            report(list(z['A']), z['M'].item(), list(z['B']) if 'B' in z.files else None, t0)
        return
    jobsB = [(eps, seed, eta) for eps in EPS_B for eta in (0.0, SIGMA0, np.inf) for seed in range(N_SEEDS)]
    with mp.Pool(N_WORKERS) as pool:
        if '--only-B' in sys.argv:                                         # Parts A and the margins from the raw file
            with np.load(raw, allow_pickle=True) as z:
                A, M = list(z['A']), z['M'].item()
                done = list(z['B']) if 'B' in z.files else []
        else:
            A = pool.map(part_A_eps, list(EPS), chunksize=1)
            lam = {r['eps']: r['lam'] for r in A}
            M = dict(pool.map(margins_job, [(eps, m, lam[eps]) for eps in EPS if eps > 0 for m in (15, 60)], chunksize=1))
            save_raw(raw, A=np.array(A, dtype=object), M=np.array(M, dtype=object))
            done = []
        B = None
        if '--no-B' not in sys.argv:
            have = {(r['eps'], r['seed'], r['eta']) for r in done}
            B = list(done)
            for r in pool.imap_unordered(run_job, [j for j in jobsB if j not in have], chunksize=1):
                B.append(r)                                                # saved after every job: an interrupted run keeps its results
                save_raw(raw, A=np.array(A, dtype=object), M=np.array(M, dtype=object), B=np.array(B, dtype=object))
            B.sort(key=lambda r: (r['eps'], r['eta'], r['seed']))
            save_raw(raw, A=np.array(A, dtype=object), M=np.array(M, dtype=object), B=np.array(B, dtype=object))
    report(A, M, B, t0)


if __name__ == '__main__':
    main()
