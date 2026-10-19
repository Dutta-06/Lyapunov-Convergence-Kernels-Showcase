"""
An application: cogging-torque compensation in a position servo — kernel adaptive control with transport leaks
=================================================================================================================

The servo ẍ = −τ_cog(x) + u (unit inertia, viscous friction compensated), τ_cog(x) = 0.4 sin(3x + 0.3) + 0.1 sin 7x unknown, tracks
x_m = sin t under the sliding controller u = ẍ_m − λė − k_s s − f̂_Θ(x), s := ė + λe (src/core/sliding.py: the error channel is
ṡ = −k_s s + r_Θ(x), so the averaged slow loop is Theorem J's along the position orbit x_m(t) ∈ [−1, 1]). The network: K₀ = 4 saturated
Gaussian kernels on the position (a_max = 4, c_max = 1.6, widths in [0.05, 1]); τ_cog has no exact four-kernel representation, so the
analysis is Proposition D‴.28's: the best weighted approximant ϑ^b(q) (q = 1/4), a nondegenerate critical point with a 1.5% fit error.

A. No dynamics: ϑ^b(q), its fit error, the whitened unregularised Hessian Ã (its spectrum: the weakest mode), the residual curvature's
   λ⁺_k and the gauge set W_g; the fixed prior's floor law (D‴.28(ii)); the collective rates — the fixed prior's ã_j + σ₀ against the
   adaptive prior's λ_j^− (η = σ₀, ∞) — and the misplaced prior's bias (an engineer's cogging model with the phase 0.5 in place of 0.3):
   the fixed prior's stationary point and its fit error.
B. The averaged loop from warm clouds around ϑ^b(q) (width 0.02), m = 16, 64, 128 against M = 256, three seeds, S_B slow units, for the
   matching and the misplaced prior: the nearest-atom leak toward the fixed and the adaptive prior (η = ∞ and σ₀; at η = ∞ the prior's
   initial position enters only through the first cells), the OT leak toward the fixed prior, and the classical σ-modification (the leak
   toward ϑ⁰ = 0): the fit error (the bias) and ‖f̂_m − f̂_M‖² (the finite-size floor).
C. The real servo (k_s = 25, λ = 5, RK4 dt = 2·10⁻³, T_C time units, m = 16 and 64, seed 0): from the stationary points of the averaged
   designs, the position-tracking error's rms over five periods against the averaged prediction — the slaved s = r(x_m)/k_s filtered by
   (d/dt + λ)⁻¹ along the orbit.

Output: results/servo/servo_cogging.md and servo_cogging_raw.npz (SV_WORKERS, SV_SB, SV_TC; `--from-raw` regenerates the md).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import ot_leak as OT
from src.core import prior_design as PD
from src.core import adaptive_prior as AP
from src.core import floor_bias as FB
from src.core.ot_leak import _cost_matrix

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'servo')
NET_KW = dict(a_max=4.0, c_max=1.6, s_min=0.05, s_max=1.0)
G = np.array([1.0, 0.25, 0.25])
PER = 2 * np.pi
N_Q = 400
K0 = 4
Q = np.ones(K0) / K0
SIGMA0 = 0.01
PHASE, PHASE_PRIOR = 0.3, 0.5
MS, M_REF = (16, 64, 128), 256
JITTER = 0.02
DT, PLAN_EVERY, REC = 0.25, 4, 200
S_B = float(os.environ.get('SV_SB', 2000.0))
N_SEEDS = 3
N_WORKERS = int(os.environ.get('SV_WORKERS', 2))
KS, LAM, DT_REAL = 25.0, 5.0, 2e-3
T_C = float(os.environ.get('SV_TC', 1000.0))


def tau_cog(X, phase=PHASE):
    X = np.asarray(X, float)
    return 0.4 * np.sin(3 * X + phase) + 0.1 * np.sin(7 * X)


def f_vec(X, phase=PHASE):
    return -tau_cog(X, phase)


def x_ref(t):
    return np.array([np.sin(t)])


def net_():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


START = (np.array([1.0, 1.4, -1.9, -0.2]), np.array([-0.9, -0.33, 0.3, 0.78]), np.array([0.25, 0.25, 0.27, 0.1]))


def approximant(phase):
    """The best weighted K₀-kernel approximant at q of −τ_cog with the given phase: BFGS on the unregularised weighted loss from START,
    polished by Newton; returns (Θ, gradient residual, relative fit error)."""
    net = net_()
    th0 = net.unconstrained(START[0][:, None], START[1][:, None], START[2])
    loop = OT.OTPriorLoop(net, lambda X: f_vec(X, phase), x_ref, PER, G, SIGMA0, th0, n_quad=N_Q)

    def fun(v):
        T = v.reshape(K0, 3); a, c, s, d, ph = loop.features(T); r = loop.F - ph @ (Q * a)
        return 0.5 * float(np.mean(r ** 2)), (Q[:, None] * loop.grad_V_L(r, T)).ravel()
    res = minimize(fun, th0.ravel(), jac=True, method='BFGS', options=dict(gtol=1e-14, maxiter=40000))
    th, rr, fit = FB.exact_weighted_atoms(loop, Q, res.x.reshape(K0, 3), n_cont=1, w0=Q, fit_tol=1.0)
    return th, rr, fit


def setup(prior_phase=PHASE, sigma0=SIGMA0, zero_prior=False):
    """The loop for the true target with the prior at the approximant of the cogging model with prior_phase (the matching prior when it
    equals PHASE), or at ϑ⁰ = 0 (σ-modification)."""
    net = net_()
    tb = approximant(PHASE)[0]
    prior = np.zeros((1, 3)) if zero_prior else approximant(prior_phase)[0]
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, sigma0, prior, n_quad=N_Q)
    return net, tb, prior, loop


# ---------------------------------------------------------------------------- A
def part_A():
    net, tb, prior, loop = setup(PHASE)
    th, rr, fit = approximant(PHASE)
    gradL = lambda v: (Q[:, None] * FB.unregularised_grad(loop, Q, v.reshape(K0, 3))).ravel()
    A = np.zeros((3 * K0, 3 * K0)); e = 1e-6
    for j in range(3 * K0):
        E = np.zeros(3 * K0); E[j] = e
        A[:, j] = (gradL(tb.ravel() + E) - gradL(tb.ravel() - E)) / (2 * e)
    Dh = np.kron(np.sqrt(Q), np.sqrt(1.0 / G))
    At = 0.5 * (A + A.T) / Dh[:, None] / Dh[None, :]
    ev = np.linalg.eigvalsh(At)
    a, c, s, d, ph = loop.features(tb)
    r = loop.F - ph @ (Q * a)
    lam = loop.kappa_R(r, tb)
    lo = Q * np.abs(a) / net.a_max
    rates = {eta: np.sort(AP.pair_rates(ev, SIGMA0, eta)[0]) for eta in (SIGMA0, np.inf)}
    # the fixed matching prior's floor law: σ₀Ĵ(Ã + σ₀I)⁻¹g̃, the gauge displacement along e_j − e_K
    Jh = FB.whitened_jacobian(loop, Q, tb)
    h = 1e-6
    Vd = np.array([np.eye(K0)[j] - np.eye(K0)[K0 - 1] for j in range(K0 - 1)]).T
    def gauge(w):
        a_, c_, s_ = net.physical(tb)
        return net.unconstrained(Q[:, None] * a_ / np.asarray(w, float)[:, None], c_, s_)
    gt = np.column_stack([-Dh * ((gauge(Q + h * Vd[:, k]) - gauge(Q - h * Vd[:, k])) / (2 * h)).ravel() for k in range(K0 - 1)])
    F = SIGMA0 * Jh @ np.linalg.solve(At + SIGMA0 * np.eye(3 * K0), gt)
    Ct = (np.diag(Q) - np.outer(Q, Q))[:K0 - 1, :K0 - 1]
    law = float(np.trace(F @ Ct @ F.T)) / loop.N
    # the misplaced prior: the fixed prior's stationary point and its fit error
    netp, tbp, priorp, loopp = setup(PHASE_PRIOR)
    Tw, resw = PD.weighted_stationary_atoms_prior(loopp.nearest, Q, priorp)
    fitw = float(np.sqrt(np.mean((loopp.F - loopp.features(Tw)[4] @ (Q * loopp.features(Tw)[0])) ** 2) / np.mean(loopp.F ** 2)))
    fit_prior = float(np.sqrt(np.mean((loopp.F - loopp.features(priorp)[4] @ (Q * loopp.features(priorp)[0])) ** 2) / np.mean(loopp.F ** 2)))
    # the same replicated point from the approximant, and its D‴.5 splitting margin (the fixed prior's replicated point a saddle if negative)
    Tw2, resw2 = PD.weighted_stationary_atoms_prior(loopp.nearest, Q, tbp)
    marg_mis, _, kap_mis = PD.splitting_margin_prior(loopp.nearest, Q, priorp)
    marg_match = float(SIGMA0 - np.max(lam))
    return dict(fit=fit, res=rr, ev=ev, lam=lam, lo=lo, rates=rates, law=law, phys=tuple(x for x in net.physical(tb)), fit_misplaced=fitw,
                res_misplaced=resw, fit_prior=fit_prior, same_point=float(np.abs(Tw2 - Tw).max()), margin_misplaced=marg_mis,
                kappa_misplaced=np.asarray(kap_mis), margin_matching=marg_match)


# ---------------------------------------------------------------------------- B
DESIGNS = ('nearest-fixed', 'nearest-adaptive', 'nearest-adaptive-slow', 'ot-fixed', 'sigma-mod')


def design_loop(design, prior_phase):
    zero = design == 'sigma-mod'
    net, tb, prior, loop = setup(prior_phase, zero_prior=zero)
    eta = np.inf if design == 'nearest-adaptive' else (SIGMA0 if design == 'nearest-adaptive-slow' else 0.0)
    plan = 'ot' if design == 'ot-fixed' else 'nearest'
    return net, tb, prior, loop, eta, plan


def run_job(args):
    design, prior_phase, seed = args
    t0 = time.time()
    net, tb, prior, loop, eta, plan = design_loop(design, prior_phase)
    clouds = CL.mixture_start(tb, M_REF, np.random.default_rng(seed), JITTER, G)
    out = dict(design=design, prior_phase=prior_phase, seed=seed, rec={}, final={})
    for m in MS + (M_REF,):
        apl = AP.AdaptivePriorLoop(loop, eta, plan=plan)
        rec = dict(s=[], fit=[], fh=[])
        Z0 = prior.copy()

        def cb(s, Th, Z, adapting, rec=rec):
            rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['fh'].append(loop.f_hat(Th).astype(np.float32))
        Th, Z, _ = apl.run(clouds[:m].copy(), Z0, S_B, DT, plan_every=PLAN_EVERY, rec_every=REC, callback=cb, gate=False)
        out['rec'][m] = {k: np.array(v) for k, v in rec.items()}
        out['final'][m] = (Th.astype(np.float32), np.asarray(Z, np.float32))
    print(f"[{time.time() - t0:.0f}s] B {design} prior phase {prior_phase} seed {seed}", flush=True)
    return out


# ---------------------------------------------------------------------------- C
def real_job(args):
    """The real servo from the end state of a Part B design (seed 0), the design's leak, T_C time units: e² every 0.1 time units."""
    design, prior_phase, m, Th0, Z0 = args
    t0 = time.time()
    net, tb, prior, loop, eta, plan = design_loop(design, prior_phase)
    Th = np.asarray(Th0, float).copy(); Z = np.asarray(Z0, float).copy()
    x, v, t = 0.0, 1.0, 0.0                                                # on the reference: x_m(0) = 0, ẋ_m(0) = 1
    sig = SIGMA0 / KS

    def targets(Th, Z):
        if design == 'ot-fixed':
            return loop.leak_target(Th)
        lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
        if np.isinf(eta):
            Z = np.array([Th[lab == k].mean(axis=0) if np.any(lab == k) else Z[k] for k in range(Z.shape[0])])
        return Z[lab]

    def F(t, x, v, Th, T):
        e = x - np.sin(t); ed = v - np.cos(t); s = ed + LAM * e
        fh = float(net.f_hat(np.array([x]), Th)[0])
        u = -np.sin(t) - LAM * ed - KS * s - fh
        vd = float(f_vec(np.array([x]))[0]) + u
        return v, vd, G[None, :] * net.grad_g(np.array([x]), np.array([s]), Th) - sig * (Th - T)
    n = int(round(T_C / DT_REAL)); ts, e2 = [], []
    marks = {int(round(tt / DT_REAL)): tt for tt in (5 * PER, 500.0, T_C)}
    snaps = {}
    T = targets(Th, Z)
    for i in range(n + 1):
        if i in marks:
            snaps[marks[i]] = Th.copy()
        if i == n:
            break
        if i % 100 == 0:
            T = targets(Th, Z)
            ts.append(t); e2.append((x - np.sin(t)) ** 2)
        k1 = F(t, x, v, Th, T); k2 = F(t + DT_REAL / 2, x + DT_REAL / 2 * k1[0], v + DT_REAL / 2 * k1[1], Th + DT_REAL / 2 * k1[2], T)
        k3 = F(t + DT_REAL / 2, x + DT_REAL / 2 * k2[0], v + DT_REAL / 2 * k2[1], Th + DT_REAL / 2 * k2[2], T)
        k4 = F(t + DT_REAL, x + DT_REAL * k3[0], v + DT_REAL * k3[1], Th + DT_REAL * k3[2], T)
        x += DT_REAL / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); v += DT_REAL / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Th = Th + DT_REAL / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]); t += DT_REAL
    # the averaged prediction at the parameters reached at each window's end: s = r(x_m)/k_s, e = (d/dt + λ)⁻¹s along the orbit
    from src.core.sliding import filtered_error
    pred = {}
    for tt, Ts in snaps.items():
        s_pred, e_pred = filtered_error(loop.F - loop.f_hat(Ts), PER, KS, LAM)
        pred[tt] = float(np.sqrt(np.mean(e_pred ** 2)))
    print(f"[{time.time() - t0:.0f}s] C {design} prior phase {prior_phase} m {m}", flush=True)
    return dict(design=design, prior_phase=prior_phase, m=m, ts=np.array(ts), e2=np.array(e2), e_pred=pred, final=Th.astype(np.float32))


# ---------------------------------------------------------------------------- report
NAMES = {'nearest-fixed': 'nearest atom, fixed prior', 'nearest-adaptive': 'nearest atom, adaptive prior (η = ∞)',
         'nearest-adaptive-slow': 'nearest atom, adaptive prior (η = σ₀)', 'ot-fixed': 'OT, fixed prior',
         'sigma-mod': 'σ-modification (ϑ⁰ = 0)'}


def class_diameter(Th, Z):
    """The largest whitened diameter of the classes of Θ by the nearest atom of Z."""
    Th = np.asarray(Th, float); Z = np.asarray(Z, float)
    lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
    out = 0.0
    for k in range(Z.shape[0]):
        X = Th[lab == k]
        if len(X) > 1:
            out = max(out, float(np.sqrt(np.max(np.sum((X[:, None, :] - X[None, :, :]) ** 2 / G[None, None, :], axis=2)))))
    return out


def report(A, B, C=None, t0=None):
    from src.core.sliding import filtered_error
    a, c, s = A['phys']
    lines = ["# An application: cogging-torque compensation in a position servo (`exp_servo_cogging.py`)", "",
             "The servo ẍ = −τ_cog(x) + u, τ_cog(x) = 0.4 sin(3x + 0.3) + 0.1 sin 7x, tracking x_m = sin t under the sliding controller (s = ė + λe); "
             f"K₀ = {K0} saturated Gaussian kernels on the position, masses 1/{K0}, σ₀ = {SIGMA0:g}; the averaged slow loop along the position orbit.", "",
             "## A. The approximant and the rates (no dynamics)", "",
             f"The best weighted approximant ϑ^b(q): a = ({', '.join(f'{x:.3f}' for x in a[:, 0])}), c = ({', '.join(f'{x:.3f}' for x in c[:, 0])}), "
             f"s = ({', '.join(f'{x:.3f}' for x in s)}); gradient residual {A['res']:.1e}, relative fit error {A['fit']:.4e} (no exact representation). "
             "The whitened unregularised Hessian Ã's spectrum (units of σ₀): " + ", ".join(f"{x / SIGMA0:.2e}" for x in A['ev']) + ".", "",
             "The collective rates (units of σ₀; the slowest three): the fixed prior's ã_j + σ₀: " + ", ".join(f"{(x + SIGMA0) / SIGMA0:.4f}" for x in A['ev'][:3])
             + "; the adaptive prior's λ_j^− at η = σ₀: " + ", ".join(f"{x / SIGMA0:.2e}" for x in A['rates'][SIGMA0][:3]) + ", at η = ∞: "
             + ", ".join(f"{x / SIGMA0:.2e}" for x in A['rates'][np.inf][:3]) + ".", "",
             "The residual curvature: λ⁺_k/σ₀ = " + ", ".join(f"{x / SIGMA0:.4f}" for x in A['lam']) + " (Proposition D‴.28(iii): margins σ₀ − (q_k/w_k)λ⁺_k); the gauge set "
             "W_g: w > (" + ", ".join(f"{x:.4f}" for x in A['lo']) + f"). The fixed matching prior's floor law (D‴.28(ii)) per unit of 1/m − 1/M: {A['law']:.4e}.", "",
             f"The misplaced prior (the approximant of the cogging model with the phase {PHASE_PRIOR:g} in place of {PHASE:g}): its own fit error on the true "
             f"target {A['fit_prior']:.4e}; the fixed prior's replicated stationary point there (residual {A['res_misplaced']:.1e}; the same from the prior's atoms "
             f"and from ϑ^b(q), to {A.get('same_point', float('nan')):.1e}) has the fit error {A['fit_misplaced']:.4e} — the leak's bias, which the adaptive prior "
             "removes" + (f" — and Proposition D‴.5's splitting margin {A['margin_misplaced']:+.2e} (κ_R/σ₀ = "
                          + ", ".join(f"{x / SIGMA0:.3f}" for x in A['kappa_misplaced']) + "): a saddle, which the m-systems leave by splitting their classes "
                          f"(Part B). At the matching prior the replicated point ϑ^b(q) has the margin {A['margin_matching']:+.2e} at w = q."
                          if 'margin_misplaced' in A else "."), "",
             f"## B. The averaged loop from warm clouds around ϑ^b(q) ({S_B:g} slow units, {N_SEEDS} seeds)", "",
             "The fit error at s = 600, 1500 and the end (mean over the seeds), ‖f̂_m − f̂_M‖² at the end (mean, min–max), and the largest whitened class "
             "diameter at the end over the seeds (the classes by the nearest atom of the design's prior at the end; —: σ-modification, a single atom).", "",
             "| prior | design | m | fit error at 600, 1500, end | ‖f̂_m − f̂_M‖² at the end | class diameter |", "|---|---|---|---|---|---|"]
    for pp in (PHASE, PHASE_PRIOR):
        for dsg in DESIGNS:
            rs = [r for r in B if r['design'] == dsg and r['prior_phase'] == pp]
            if not rs:
                continue
            if dsg == 'sigma-mod' and pp != PHASE:
                continue
            for m in MS:
                fits = np.array([[r['rec'][m]['fit'][int(np.argmin(np.abs(r['rec'][m]['s'] - t)))] for t in (600.0, 1500.0, S_B)] for r in rs])
                h2 = np.array([float(np.mean((r['rec'][m]['fh'][-1].astype(float) - r['rec'][M_REF]['fh'][-1].astype(float)) ** 2)) for r in rs])
                pname = 'σ-mod' if dsg == 'sigma-mod' else ('matching' if pp == PHASE else f'misplaced (phase {pp:g})')
                diam = '—' if dsg == 'sigma-mod' else f"{max(class_diameter(*r['final'][m]) for r in rs):.2f}"
                lines.append(f"| {pname} | {NAMES[dsg]} | {m} | " + ", ".join(f"{x:.3e}" for x in fits.mean(axis=0)) + f" | {h2.mean():.2e} ({h2.min():.1e}–{h2.max():.1e}) | {diam} |")
    if C:
        lines += ["", f"## C. The real servo (k_s = {KS:g}, λ = {LAM:g}, dt = {DT_REAL:g}, {T_C:g} time units), seed 0", "",
                  "From the end state of each averaged design: the position-tracking error's rms over the five periods ending at the times below, against "
                  "the averaged prediction at the parameters reached then (the slaved s = r(x_m)/k_s filtered by (d/dt + λ)⁻¹ along the orbit), and at the end "
                  "against the prediction with the error filter kept exactly, s = (d/dt + k_s)⁻¹r(x_m).", "",
                  "| prior | design | m | e_rms at t = 31, 500, end | predicted | ratios | exact filter at the end: predicted, ratio |", "|---|---|---|---|---|---|---|"]
        win = 5 * PER
        for cc in sorted(C, key=lambda c: (c['prior_phase'], DESIGNS.index(c['design']), c['m'])):
            er = []
            for tt in (win, 500.0, T_C):
                msk = (cc['ts'] > tt - win) & (cc['ts'] <= tt + 1e-9)
                er.append(float(np.sqrt(np.mean(cc['e2'][msk]))))
            pname = 'σ-mod' if cc['design'] == 'sigma-mod' else ('matching' if cc['prior_phase'] == PHASE else 'misplaced')
            pr = [cc['e_pred'][tt] for tt in (win, 500.0, T_C)]
            net_, tb_, prior_, loop_, eta_, plan_ = design_loop(cc['design'], cc['prior_phase'])
            e_ex = filtered_error(loop_.F - loop_.f_hat(np.asarray(cc['final'], float)), PER, KS, LAM, exact=True)[1]
            pe = float(np.sqrt(np.mean(e_ex ** 2)))
            lines.append(f"| {pname} | {NAMES[cc['design']]} | {cc['m']} | " + ", ".join(f"{x:.3e}" for x in er) + " | " + ", ".join(f"{x:.3e}" for x in pr) + " | "
                         + ", ".join(f"{a / b:.3f}" for a, b in zip(er, pr)) + f" | {pe:.3e}, {er[-1] / pe:.4f} |")
    lines += ["", "Numbers are produced by `experiments/exp_servo_cogging.py`; nothing is edited by hand."]
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'servo_cogging.md'), 'w', encoding='utf-8') as f:
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
    os.makedirs(OUT, exist_ok=True)
    raw = os.path.join(OUT, 'servo_cogging_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            report(z['A'].item(), list(z['B']), list(z['C']) if 'C' in z.files else None, t0)
        return
    if '--only-A' in sys.argv:                                              # Part A recomputed, merged with the stored B and C
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        save['A'] = np.array(part_A(), dtype=object)
        save_raw(raw, **save)
        report(save['A'].item(), list(save['B']), list(save['C']) if 'C' in save else None, t0)
        return
    A = part_A()
    jobs = [(d, pp, seed) for pp in (PHASE, PHASE_PRIOR) for d in DESIGNS if not (d == 'sigma-mod' and pp != PHASE) for seed in range(N_SEEDS)]
    with mp.Pool(N_WORKERS) as pool:
        B = pool.map(run_job, jobs, chunksize=1)
        save_raw(raw, A=np.array(A, dtype=object), B=np.array(B, dtype=object))
        cj = [(r['design'], r['prior_phase'], m, r['final'][m][0], r['final'][m][1]) for r in B if r['seed'] == 0 for m in (16, 64)]
        C = pool.map(real_job, cj, chunksize=1)
    save_raw(raw, A=np.array(A, dtype=object), B=np.array(B, dtype=object), C=np.array(C, dtype=object))
    report(A, B, C, t0)


if __name__ == '__main__':
    main()
