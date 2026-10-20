"""
Proposition D‴.32: the cohesion leak — a leak on the unexcited subspace, without a prior, that selects the number of kernels
==========================================================================================================================

The three-bump instance (σ₀ = 0.15σ_crit), the averaged slow loop with the cohesion leak (src/core/cohesion_leak.py): the particles attract
each other with the kernel ψ(‖ϑ_i − ϑ_j‖²_{G⁻¹}/R²), compact ((1 − z)³₊, range R) or Gaussian (e^{−z/2}, scale R), instead of leaking toward a
prior. Midpoint rule in slow time, dt = 0.25.

A. The linear theory (no dynamics). (1) At every count vector of m = 15 in W the lift of the exact representation ϑ^e(w) (amplitudes A_k/w_k,
   the gauge in closed form): the largest drift (compact and Gaussian, R = 1), the smallest inter-atom distance, and the Jacobian's spectrum
   (central differences) against Proposition D‴.32(iii): −σ₀w_k + eig(R^G_k) on the separations of class k (R = 0 at an exact
   representation) and the unregularised weighted K-atom flow on the class means. (2) The Gaussian kernel's bias: the replicated stationary
   point near ϑ^e(w) by Newton on the class drift, R = 1, 1.5, 2, 3, its whitened distance from ϑ^e(w) against the first-order
   prediction Ã_w⁻¹b and its fit error. (3) The splitting margins with an approximation error: the fourth bump of Proposition D‴.28
   (ε = 0.3), the approximant ϑ^b(q) and its gauge transforms; σ₀w_k − (q_k/w_k)λ⁺_k against the Jacobian at three count vectors of m = 15,
   at σ₀ and at σ₀/8, where a margin turns negative; runs of 100000 slow units from the replicated configuration with a small jitter at
   σ₀/8 (counts (5, 5, 5): every margin positive; (6, 4, 5): one negative) and at σ₀ — whether the predicted class splits, and what the
   split does to the fit. (4) Degenerate representations: twins (bump 3 carried by two classes with the same centre and width, masses 3/15 and 2/15, its
   weighted amplitude shared 3:1)
   and a dead class (amplitude 0, centre 2.4, width 0.6, mass 2/15): runs with the Gaussian kernel (R = 1 for the twins, R = 2 for the dead
   class) and with the compact kernel at R = 1, the twins' distance, the dead class's distance from the nearest live class and the fit error.
B. The floor at warm starts (Proposition D‴.24's clouds: width 0.02 around the exact representation at q, m = 15, 60, 240 against M = 510, the
   m-system the first m rows, four seeds), 6000 slow units: the nearest-atom leak toward the fixed matching prior (η = 0, the floor) and
   toward the adaptive prior (η = ∞), the cohesion leak (compact and Gaussian, R = 1): ‖f̂_m − f̂_M‖² and the fit error along the run, the
   class diameters at the end, the class means' whitened distance from ϑ^e(w^m).
C. The unknown model order (broad starts: a ~ U(−0.2, 0.2), c ~ U(−2.5, 2.5), s ~ U(0.4, 1), m = 60, six seeds; m = 120, three seeds):
   the Gaussian cohesion leak at R = 2 for 20000 slow units and then at R = 1 for 20000 (the annealed scale); the compact one at R = 1;
   no leak; the σ-modification toward 0; the nearest-atom leak toward an adaptive prior (η = ∞) with K = 2, 3, 4 atoms started at the first
   K particles. At the end: the number of clusters (single linkage at 0.25), their masses, centres, widths and |a|/a_max, the fit error.
D. The real closed loop (k_e = 25, σ = σ₀/k_e): (1) the Floquet exponents of the lift of ϑ^e(w), w = (6, 5, 4)/15, under the compact
   cohesion leak (R = 1; the monodromy by central differences of the period map, DOP853 at rtol 10⁻¹³) against σ₀w_k/k_e on the separations
   and the collapsed variational equation without a leak (Corollary D‴.26 at η = ∞) on the collapsed subspace; (2) from the warm clouds of
   seed 0 (m = 15 and 60; RK4 plant, Euler parameters, dt = 10⁻³, 3000 time units) the rms tracking error under the cohesion leak and under
   the nearest-atom leak toward the fixed matching prior.

E. Dither (the averaged loop with an independent Brownian forcing √(2τG) per particle, Euler–Maruyama at dt = 0.25, τ = 10⁻⁵), from the lift
   of the exact representation at q, m = 60, 20000 slow units, three seeds: no leak, the cohesion leak (compact and Gaussian, R = 1), the
   nearest-atom leak toward the fixed matching prior, the σ-modification toward 0: the within-class rms spread (whitened, per coordinate)
   against τ/(σ₀w_k) for the leaks that act on the separations and 2τs for none, the class means' whitened deviation, and the fit error
   averaged over the second half.

F. An application: the cogging servo of Proposition D‴.29 (experiments/exp_servo_cogging.py's instance and averaged loop on the position
   orbit, σ₀ = 0.01), whose torque four kernels fit to 1.46% at best in the interior (ϑ^b(q)) and to 0.64% at the saturation. Broad starts
   (a ~ U(−0.2, 0.2), c ~ U(−1.2, 1.2), s ~ U(0.1, 0.6)), m = 64, three seeds: the Gaussian cohesion leak at R = 0.7 for 20000 slow units,
   then at R = 0.35 for 60000; at the end the clusters (masses, |a|/a_max) and the fit error, against the best equal-weight K-atom fits
   (K = 4, 6, 8; 40 random starts each; no model of the torque anywhere).

Output: results/chaos/cohesion_leak.md and cohesion_leak_raw.npz (CO_WORKERS, CO_OUT; `--only A,B,C,D,E,F` runs the listed parts and merges
them into the raw file; `--from-raw` regenerates the md).
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
from src.core import adaptive_prior as AP
from src.core import adaptive_floquet as AF
from src.core import cohesion_leak as CO
from src.core import periodic as P
from src.core.chaos_learning import AveragedSlowLoop
from src.core.ot_leak import _cost_matrix

OUT = os.environ.get('CO_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
A_BUMP = np.array([b[0] for b in BUMPS]); C_BUMP = np.array([b[1] for b in BUMPS]); S_BUMP = np.array([b[2] for b in BUMPS])
EXTRA = (0.8, 0.3)
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q = np.ones(3) / 3
THR = np.array([10, 7, 8])                                                 # |A_k|/a_max = THR_k/40 exactly
DT = 0.25
R_WARM = 1.0
MS, M_REF, JITTER, N_SEEDS, S_B, REC_B = (15, 60, 240), 510, 0.02, 4, 6000.0, 200
R_BROAD, R_FINE, S_C1, S_C2 = 2.0, 1.0, 20000.0, 20000.0
SEEDS_C = {60: tuple(range(6)), 120: (0, 1, 2)}
KE, DT_REAL, T_REAL = 25.0, 1e-3, 3000.0
N_WORKERS = int(os.environ.get('CO_WORKERS', 2))
CLUSTER_THR = 0.25
TAU_E, S_E, SEEDS_E, REC_E = 1e-5, 20000.0, (0, 1, 2), 400
S_SPLIT = 100000.0


def f_vec(X, eps=0.0):
    X = np.asarray(X, float)
    out = sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)
    return out + eps * np.exp(-(X - EXTRA[0]) ** 2 / (2 * EXTRA[1] ** 2)) if eps else out


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t)])


def setup(eps=0.0, sigma0=SIGMA0):
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    base = AveragedSlowLoop(net, lambda X: f_vec(X, eps), x_ref, PER, G, sigma0, n_quad=N_Q)
    return net, base


def exact_rep(net, w):
    """ϑ^e(w): the bumps with the amplitudes A_k/w_k (the amplitude–mass gauge in closed form)."""
    return net.unconstrained((A_BUMP / np.asarray(w, float))[:, None], C_BUMP[:, None], S_BUMP)


def representable(counts):
    n = np.asarray(counts, int)
    return bool(np.all(40 * n - THR * int(n.sum()) > 0))


def count_vectors(m):
    return [np.array([a, b, m - a - b]) for a in range(1, m) for b in range(1, m - a) if representable((a, b, m - a - b))]


def class_drift(base, w, atoms, R, kind, sigma0):
    """The drift of a replicated configuration's classes (one row per class): learning with the weights w plus the attraction of the other
    classes, −σ₀Σ_lw_lψ(d²_kl/R²)(ϑ_k − ϑ_l)."""
    a, c, s, diff, phi = base.features(atoms)
    learn = G[None, :] * base.grad_weighted(base.F - phi @ (w * a), atoms)
    D, d2 = CO.pair_d2(atoms, G)
    return learn - sigma0 * np.einsum('l,kl,klp->kp', w, CO.psi(d2 / R ** 2, kind), D)


def whiten(dev, w):
    return (np.sqrt(np.asarray(w, float))[:, None] * dev / np.sqrt(G)[None, :]).ravel()


def collective_hessian(base, w, atoms, h=1e-6):
    """Ã_w: minus the Jacobian of the unregularised weighted K-atom drift in the whitened weighted coordinates."""
    K = atoms.shape[0]; d = 3 * K
    J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (class_drift(base, w, (atoms.ravel() + e).reshape(K, 3), 1.0, 'compact', 0.0)
                   - class_drift(base, w, (atoms.ravel() - e).reshape(K, 3), 1.0, 'compact', 0.0)).ravel() / (2 * h)
    S = np.repeat(np.sqrt(w), 3) / np.tile(np.sqrt(G), K)
    return -(S[:, None] * J / S[None, :])


def newton_classes(base, w, atoms0, R, kind, sigma0, iters=40, h=1e-7):
    """The replicated stationary point near atoms0: Newton on the class drift (finite-difference Jacobian)."""
    x = np.array(atoms0, float).ravel(); K = len(w)
    F = lambda v: class_drift(base, w, v.reshape(K, 3), R, kind, sigma0).ravel()
    for _ in range(iters):
        f0 = F(x)
        if np.abs(f0).max() < 1e-15:
            break
        J = np.zeros((x.size, x.size))
        for j in range(x.size):
            e = np.zeros(x.size); e[j] = h
            J[:, j] = (F(x + e) - F(x - e)) / (2 * h)
        x = x - np.linalg.solve(J, f0)
    return x.reshape(K, 3), float(np.abs(F(x)).max())


# ---------------------------------------------------------------------------- A. the linear theory
def part_A():
    t0 = time.time()
    net, base = setup()
    out = {}
    # (1) every count vector of m = 15 in W
    rows = []
    for counts in count_vectors(15):
        w = counts / counts.sum(); at = exact_rep(net, w); Th = np.repeat(at, counts, axis=0)
        D, d2 = CO.pair_d2(at, G)
        row = dict(counts=counts, dmin=float(np.sqrt(d2[np.triu_indices(3, 1)]).min()))
        for kind in ('compact', 'gauss'):
            loop = CO.CohesionLoop(base, R_WARM, kind)
            ev = np.sort(np.linalg.eigvals(loop.jacobian(Th)).real)
            sep, coll = CO.predicted_spectrum(base, at, counts, SIGMA0)
            row[kind] = dict(drift=float(np.abs(loop.drift(Th)).max()), spec_err=float(np.abs(ev - np.sort(np.concatenate([sep, coll]))).max()),
                             slowest=float(-ev[-1]), sep_pred=-SIGMA0 * w, coll_min=float(-coll[-1]))
        rows.append(row)
    out['cv'] = rows
    # (2) the Gaussian kernel's bias
    bias = []
    for counts in ([5, 5, 5], [6, 5, 4]):
        counts = np.array(counts); w = counts / 15; at = exact_rep(net, w)
        At = collective_hessian(base, w, at)
        for R in (1.0, 1.5, 2.0, 3.0):
            st, res = newton_classes(base, w, at, R, 'gauss', SIGMA0)
            y = whiten(st - at, w)
            b = whiten(class_drift(base, w, at, R, 'gauss', SIGMA0), w)            # the attraction at ϑ^e(w), whitened (learning there is zero)
            S = np.repeat(np.sqrt(w), 3) / np.tile(np.sqrt(G), 3)
            y1 = np.linalg.solve(At, b)                                             # first order: Ã y = b in the whitened coordinates
            bias.append(dict(counts=counts, R=R, res=res, y=float(np.linalg.norm(y)), y1=float(np.linalg.norm(y1)),
                             rel=float(np.linalg.norm(y - y1) / max(np.linalg.norm(y), 1e-300)), fit=base.fit_error(np.repeat(st, counts, axis=0))))
    out['bias'] = bias
    # (3) the splitting margins with an approximation error (Proposition D‴.28's fourth bump, ε = 0.3)
    out['split'] = split_margins()
    # (4) degenerate representations
    out['degen'] = degenerate_runs()
    print(f"[{time.time() - t0:.0f}s] A", flush=True)
    return out


def approximant(eps_target=0.3, steps=6):
    """ϑ^b(q): the weighted three-atom approximant of f_ε, continued in ε from the exact representation (BFGS on the unregularised loss)."""
    net, _ = setup()
    x = exact_rep(net, Q).ravel()
    for eps in np.linspace(0, eps_target, steps + 1)[1:]:
        _, base = setup(eps)

        def fun(v):
            Th = v.reshape(3, 3); a, c, s, diff, phi = base.features(Th)
            return 0.5 * np.mean((base.F - phi @ (Q * a)) ** 2)

        def grad(v):
            Th = v.reshape(3, 3); a, c, s, diff, phi = base.features(Th)
            return -(Q[:, None] * base.grad_weighted(base.F - phi @ (Q * a), Th)).ravel()
        x = minimize(fun, x, jac=grad, method='BFGS', options=dict(gtol=1e-14, maxiter=50000)).x
    return net, base, x.reshape(3, 3)


def gauge(net, atoms_q, w):
    a, c, s = net.physical(atoms_q)
    return net.unconstrained(Q[:, None] * a / np.asarray(w, float)[:, None], c, s)


def split_margins():
    net, base, bq = approximant()
    a, c, s, diff, phi = base.features(bq)
    r = base.F - phi @ (Q * a)
    lam = base.kappa_R(r, bq)                                                   # λ⁺_k at ϑ^b(q)
    gnorm = float(np.abs(Q[:, None] * base.grad_weighted(r, bq)).max())
    rows = []
    for counts in ([5, 5, 5], [6, 4, 5], [4, 6, 5]):
        counts = np.array(counts); w = counts / 15; at = gauge(net, bq, w); Th = np.repeat(at, counts, axis=0)
        for frac in (1.0, 1.0 / 8):
            sig = SIGMA0 * frac
            loop = CO.CohesionLoop(base, R_WARM, 'compact', sigma0=sig)
            ev = np.linalg.eigvals(loop.jacobian(Th)).real
            sep, coll = CO.predicted_spectrum(base, at, counts, sig)
            pred_margin = sig * w - (Q / w) * lam
            rows.append(dict(counts=counts, frac=frac, margins=pred_margin, top_sep=float(sep[-1]), top_ev=float(np.sort(ev)[-1]),
                             spec_err=float(np.abs(np.sort(ev) - np.sort(np.concatenate([sep, coll]))).max()), drift=float(np.abs(loop.drift(Th)).max())))
    # runs from the replicated configuration with a small jitter: does the predicted class split, and what does the split do?
    rng = np.random.default_rng(0)
    runs = []
    for counts, frac in (((5, 5, 5), 1.0 / 8), ((6, 4, 5), 1.0 / 8), ((6, 4, 5), 1.0)):
        counts = np.array(counts); w = counts / 15; at = gauge(net, bq, w)
        Th0 = np.repeat(at, counts, axis=0) + 1e-3 * rng.standard_normal((15, 3)) * np.sqrt(G)[None, :]
        loop = CO.CohesionLoop(base, R_WARM, 'compact', sigma0=SIGMA0 * frac)
        rec = []

        def cb(s_, T, rec=rec, counts=counts):
            lab = np.repeat(np.arange(3), counts)
            rec.append([s_] + [float(np.sqrt(np.sum((T[lab == k] - T[lab == k].mean(axis=0)) ** 2 / G[None, :], axis=1)).max()) for k in range(3)]
                       + [base.fit_error(T), CO.clusters(T, G, CLUSTER_THR).max() + 1])
        Th = loop.run(Th0, S_SPLIT, DT, rec_every=int(S_SPLIT / DT / 10), callback=cb)
        runs.append(dict(counts=counts, frac=frac, rec=np.array(rec), end=describe(net, base, Th)))
    return dict(lam=lam, grad=gnorm, fit_b=base.fit_error(np.repeat(bq, 5, axis=0)), rows=rows, runs=runs)


def degenerate_runs():
    net, base = setup()
    out = {}
    # twins: bump 3 by two classes with the same centre and width (masses 3/15 and 2/15), its weighted amplitude shared 3:1
    counts = np.array([5, 5, 3, 2]); w = counts / 15
    amps = np.array([A_BUMP[0] / w[0], A_BUMP[1] / w[1], 0.75 * A_BUMP[2] / w[2], 0.25 * A_BUMP[2] / w[3]])
    at = net.unconstrained(amps[:, None], np.r_[C_BUMP, C_BUMP[2]][:, None], np.r_[S_BUMP, S_BUMP[2]])
    Th0 = np.repeat(at, counts, axis=0)
    lab = np.repeat(np.arange(4), counts)
    for kind, R in (('gauss', 1.0), ('compact', 1.0)):
        loop = CO.CohesionLoop(base, R, kind)
        rec = []

        def cb(s_, T, rec=rec):
            m2, m3 = T[lab == 2].mean(axis=0), T[lab == 3].mean(axis=0)
            rec.append((s_, float(np.sqrt(np.sum((m2 - m3) ** 2 / G))), base.fit_error(T)))
        loop.run(Th0, 40000.0, DT, rec_every=8000, callback=cb)
        out[('twins', kind)] = np.array(rec)
    out['twins_d0'] = float(np.sqrt(np.sum((at[2] - at[3]) ** 2 / G)))
    # a dead class: amplitude 0 at c = 2.4, s = 0.6, mass 2/15; the bumps by classes of masses (5, 4, 4)/15
    counts = np.array([5, 4, 4, 2]); w = counts / 15
    amps = np.array([A_BUMP[0] / w[0], A_BUMP[1] / w[1], A_BUMP[2] / w[2], 0.0])
    at = net.unconstrained(amps[:, None], np.r_[C_BUMP, 2.4][:, None], np.r_[S_BUMP, 0.6])
    Th0 = np.repeat(at, counts, axis=0)
    lab = np.repeat(np.arange(4), counts)
    for kind, R in (('gauss', 2.0), ('compact', 1.0)):
        loop = CO.CohesionLoop(base, R, kind)
        rec = []

        def cb(s_, T, rec=rec):
            md = T[lab == 3].mean(axis=0)
            live = np.array([T[lab == k].mean(axis=0) for k in range(3)])
            rec.append((s_, float(np.sqrt(np.sum((live - md) ** 2 / G, axis=1)).min()), base.fit_error(T), float(np.sqrt(np.sum(T[lab == 3].var(axis=0) / G)))))
        loop.run(Th0, 40000.0, DT, rec_every=8000, callback=cb)
        out[('dead', kind)] = np.array(rec)
    out['dead_d0'] = float(np.sqrt(np.sum((at[:3] - at[3]) ** 2 / G, axis=1)).min())
    return out


# ---------------------------------------------------------------------------- B. the floor at warm starts
def b_job(args):
    seed, leak = args
    t0 = time.time()
    net, base = setup()
    exact = exact_rep(net, Q)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    sizes = MS + (M_REF,)
    Th = {m: clouds[:m].copy() for m in sizes}
    w0 = {m: np.bincount(np.argmin(_cost_matrix(Th[m], exact, G), axis=1), minlength=3) / Th[m].shape[0] for m in sizes}
    if leak in ('fixed', 'adaptive'):
        oloop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact, n_quad=N_Q)
        apl = AP.AdaptivePriorLoop(oloop, 0.0 if leak == 'fixed' else np.inf, plan='nearest')
        Z = {m: exact.copy() for m in sizes}
    else:
        loop = CO.CohesionLoop(base, R_WARM, 'compact' if leak == 'coh_c' else 'gauss')
    rec = dict(s=[], h2={m: [] for m in MS}, fit={m: [] for m in sizes})
    n = int(round(S_B / DT))
    gamma = {}
    for k in range(n + 1):
        if leak in ('fixed', 'adaptive') and k % 4 == 0:
            for m in sizes:
                gamma[m], _ = apl.plan(Th[m], Z[m])
                if leak == 'adaptive':
                    Z[m] = apl.bary(Th[m], gamma[m], Z[m]); gamma[m], _ = apl.plan(Th[m], Z[m])
        if k % REC_B == 0:
            rec['s'].append(k * DT)
            fM = base.f_hat(Th[M_REF])
            for m in sizes:
                rec['fit'][m].append(base.fit_error(Th[m]))
            for m in MS:
                rec['h2'][m].append(float(np.mean((base.f_hat(Th[m]) - fM) ** 2)))
        if k == n:
            break
        for m in sizes:
            if leak in ('fixed', 'adaptive'):
                k1, _ = apl.drift(Th[m], Z[m], gamma[m]); k2, _ = apl.drift(Th[m] + 0.5 * DT * k1, Z[m], gamma[m])
                Th[m] = Th[m] + DT * k2
                if leak == 'adaptive':
                    Z[m] = apl.bary(Th[m], gamma[m], Z[m])
            else:
                k1 = loop.drift(Th[m]); Th[m] = Th[m] + DT * loop.drift(Th[m] + 0.5 * DT * k1)
    end = {}
    for m in sizes:
        lab = np.argmin(_cost_matrix(Th[m], exact, G), axis=1)
        diam = max(float(np.sqrt(np.sum((Th[m][lab == j] - Th[m][lab == j].mean(axis=0)) ** 2 / G[None, :], axis=1)).max()) for j in range(3) if np.any(lab == j))
        means = np.array([Th[m][lab == j].mean(axis=0) for j in range(3)])
        dev = float(np.linalg.norm(whiten(means - exact_rep(net, w0[m]), w0[m]))) if representable(np.round(w0[m] * Th[m].shape[0])) else float('nan')
        end[m] = dict(diam=diam, dev=dev, K=int(CO.clusters(Th[m], G, CLUSTER_THR).max() + 1))
    print(f"[{time.time() - t0:.0f}s] B seed {seed} {leak}", flush=True)
    return dict(seed=seed, leak=leak, w0=w0, s=np.array(rec['s']), h2={m: np.array(v) for m, v in rec['h2'].items()},
                fit={m: np.array(v) for m, v in rec['fit'].items()}, end=end)


# ---------------------------------------------------------------------------- C. the unknown model order
def broad(net, m, seed):
    rng = np.random.default_rng(seed)
    return net.unconstrained(rng.uniform(-0.2, 0.2, (m, 1)), rng.uniform(-2.5, 2.5, (m, 1)), rng.uniform(0.4, 1.0, m))


def describe(net, base, Th):
    lab = CO.clusters(Th, G, CLUSTER_THR)
    cl = []
    for k in np.argsort(-np.bincount(lab)):
        idx = lab == k
        a, c, s = net.physical(Th[idx].mean(axis=0, keepdims=True))
        cl.append((float(idx.mean()), float(a[0, 0]), float(c[0, 0]), float(s[0]), float(np.sqrt(np.sum(Th[idx].var(axis=0) / G)))))
    return dict(K=len(cl), clusters=np.array(cl), fit=base.fit_error(Th))


def c_job(args):
    m, seed, leak = args
    t0 = time.time()
    net, base = setup()
    Th0 = broad(net, m, seed)
    trace = []
    cb = lambda s_, T: trace.append((s_, base.fit_error(T), CO.clusters(T, G, CLUSTER_THR).max() + 1))
    if leak == 'anneal':
        Th = CO.CohesionLoop(base, R_BROAD, 'gauss').run(Th0, S_C1, DT, rec_every=4000, callback=cb)
        mid = describe(net, base, Th)
        Th = CO.CohesionLoop(base, R_FINE, 'gauss').run(Th, S_C2, DT, rec_every=4000, callback=lambda s_, T: cb(s_ + S_C1, T))
    elif leak == 'coh_c':
        mid = None; Th = CO.CohesionLoop(base, R_FINE, 'compact').run(Th0, S_C1 + S_C2, DT, rec_every=4000, callback=cb)
    elif leak in ('none', 'sigma'):
        mid = None
        sig = SIGMA0 if leak == 'sigma' else 0.0
        drift = lambda T: G[None, :] * base.grad_weighted(base.residual(T), T) - sig * T
        Th = np.array(Th0, float); n = int(round((S_C1 + S_C2) / DT))
        for k in range(n + 1):
            if k % 4000 == 0:
                cb(k * DT, Th)
            if k == n:
                break
            k1 = drift(Th); Th = Th + DT * drift(Th + 0.5 * DT * k1)
    else:                                                                       # 'adaptiveK': the nearest-atom adaptive prior, η = ∞, K atoms
        K = int(leak[-1]); mid = None
        oloop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, Th0[:K].copy(), q=np.ones(K) / K, n_quad=N_Q)
        apl = AP.AdaptivePriorLoop(oloop, np.inf, plan='nearest')
        Th, Z, _ = apl.run(Th0, Th0[:K].copy(), S_C1 + S_C2, DT, plan_every=4, rec_every=4000, callback=lambda s_, T, Z_, ad: cb(s_, T))
    print(f"[{time.time() - t0:.0f}s] C m {m} seed {seed} {leak}", flush=True)
    return dict(m=m, seed=seed, leak=leak, end=describe(net, base, Th), mid=mid, trace=np.array(trace, float))


# ---------------------------------------------------------------------------- D. the real closed loop
def part_D_floquet():
    t0 = time.time()
    net, base = setup()
    counts = np.array([6, 5, 4]); w = counts / 15; at = exact_rep(net, w); m = 15
    cl = CO.CohesionRealLoop(net, lambda X: f_vec(X), x_ref, x_ref_dot, PER, KE, G, SIGMA0, R_WARM, 'compact')
    z0 = cl.lift(at, counts)
    ret = float(np.abs(AF.period_map(cl, z0, m) - z0).max())
    M = AF.monodromy(cl, z0, m)
    mult = np.linalg.eigvals(M)
    rates = np.sort(-np.log(np.abs(mult)) / PER)
    cv = AF.CollapsedVariational(net, x_ref, PER, G, SIGMA0, np.inf, at, w)
    coll = np.sort(cv.exponents(KE)[0])                                       # the slow exponents (the error channel's by Liouville apart)
    sep_pred = np.sort(np.repeat(SIGMA0 * w / KE, 3 * (counts - 1)))
    print(f"[{time.time() - t0:.0f}s] D floquet", flush=True)
    return dict(ret=ret, rates=rates, coll=coll, sep_pred=sep_pred, counts=counts)


def real_run(net, Th0, Z0, leak, T):
    """(CL) with the cohesion leak (compact, R = 1) or the nearest-atom leak toward the fixed prior Z0: e² every 0.1 time units."""
    x = x_ref(0.0).copy(); Th = np.array(Th0, float); Z = np.array(Z0, float)
    n = int(round(T / DT_REAL)); t = 0.0; ts, e2 = [], []
    for i in range(n):
        if i % 100 == 0:
            if leak == 'fixed':
                tgt = Z[np.argmin(_cost_matrix(Th, Z, G), axis=1)]
            ts.append(t); e2.append(float((x - x_ref(t))[0] ** 2))
        e = x - x_ref(t)
        xdot = lambda y, tau: f_vec(y) + x_ref_dot(tau) - KE * (y - x_ref(tau)) - net.f_hat(y, Th)
        k1 = xdot(x, t); k2 = xdot(x + 0.5 * DT_REAL * k1, t + 0.5 * DT_REAL); k3 = xdot(x + 0.5 * DT_REAL * k2, t + 0.5 * DT_REAL)
        k4 = xdot(x + DT_REAL * k3, t + DT_REAL)
        if leak == 'fixed':
            lk = -(SIGMA0 / KE) * (Th - tgt)
        else:
            lk = CO.cohesion_force(Th, G, R_WARM, SIGMA0 / KE, 'compact')
        Th = Th + DT_REAL * (G[None, :] * net.grad_g(x, e, Th) + lk)
        x = x + DT_REAL / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t = (i + 1) * DT_REAL
    return np.array(ts), np.array(e2), Th


def d_job(args):
    m, leak = args
    t0 = time.time()
    net, base = setup()
    exact = exact_rep(net, Q)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(0), JITTER, G)
    Th0 = clouds[:m].copy()
    ts, e2, Th = real_run(net, Th0, exact, leak, T_REAL)
    counts = np.bincount(np.argmin(_cost_matrix(Th0, exact, G), axis=1), minlength=3)
    print(f"[{time.time() - t0:.0f}s] D m {m} {leak}", flush=True)
    return dict(m=m, leak=leak, ts=ts, e2=e2, counts=counts, fit_end=base.fit_error(Th), K_end=int(CO.clusters(Th, G, CLUSTER_THR).max() + 1))


# ---------------------------------------------------------------------------- E. dither
def e_job(args):
    seed, leak = args
    t0 = time.time()
    net, base = setup()
    exact = exact_rep(net, Q)
    m = 60
    counts = np.array([20, 20, 20]); lab = np.repeat(np.arange(3), counts)
    Th = np.repeat(exact, counts, axis=0)
    rng = np.random.default_rng(100 + seed)
    if leak in ('coh_c', 'coh_g'):
        loop = CO.CohesionLoop(base, R_WARM, 'compact' if leak == 'coh_c' else 'gauss')
        drift = loop.drift
    elif leak == 'fixed':
        drift = lambda T: G[None, :] * base.grad_weighted(base.residual(T), T) - SIGMA0 * (T - exact[np.argmin(_cost_matrix(T, exact, G), axis=1)])
    elif leak == 'sigma':
        drift = lambda T: G[None, :] * base.grad_weighted(base.residual(T), T) - SIGMA0 * T
    else:
        drift = lambda T: G[None, :] * base.grad_weighted(base.residual(T), T)
    n = int(round(S_E / DT)); amp = np.sqrt(2 * TAU_E * DT) * np.sqrt(G)[None, :]
    rec = []
    for k in range(n + 1):
        if k % REC_E == 0:
            means = np.array([Th[lab == j].mean(axis=0) for j in range(3)])
            spread = float(np.sqrt(np.mean([np.mean(np.sum((Th[lab == j] - means[j]) ** 2 / G[None, :], axis=1)) / 3 for j in range(3)])))
            dev = float(np.linalg.norm(whiten(means - exact, Q)))
            rec.append((k * DT, spread, dev, base.fit_error(Th)))
        if k == n:
            break
        Th = Th + DT * drift(Th) + amp * rng.standard_normal(Th.shape)
    print(f"[{time.time() - t0:.0f}s] E seed {seed} {leak}", flush=True)
    return dict(seed=seed, leak=leak, rec=np.array(rec))


# ---------------------------------------------------------------------------- F. the cogging servo
SV_NET = dict(a_max=4.0, c_max=1.6, s_min=0.05, s_max=1.0)
SV_SIGMA0, SV_NQ, SV_R0, SV_R1, SV_S1, SV_S2, SV_M = 0.01, 400, 0.7, 0.35, 20000.0, 60000.0, 64


def sv_setup():
    tau = lambda X: 0.4 * np.sin(3 * np.asarray(X, float) + 0.3) + 0.1 * np.sin(7 * np.asarray(X, float))
    net = L.SaturatedLiftedNetwork(1, **SV_NET)
    base = AveragedSlowLoop(net, lambda X: -tau(X), lambda t: np.array([np.sin(t)]), PER, G, SV_SIGMA0, n_quad=SV_NQ)
    return net, base


def f_job(args):
    kind, seed = args
    t0 = time.time()
    net, base = sv_setup()
    rng = np.random.default_rng(seed)
    if kind == 'best':                                                          # the best equal-weight K-atom fit, K = seed
        K = seed; w = np.ones(K) / K; best = None
        for trial in range(40):
            th0 = net.unconstrained(rng.uniform(-1, 1, (K, 1)), rng.uniform(-1.2, 1.2, (K, 1)), rng.uniform(0.1, 0.6, K))

            def fun(v):
                T = v.reshape(K, 3); a, c, s_, d, ph = base.features(T); r = base.F - ph @ (w * a)
                return 0.5 * float(np.mean(r ** 2)), -(w[:, None] * base.grad_weighted(r, T)).ravel()
            T = minimize(fun, th0.ravel(), jac=True, method='BFGS', options=dict(gtol=1e-12, maxiter=20000)).x.reshape(K, 3)
            a, c, s_, d, ph = base.features(T)
            fit = float(np.sqrt(np.mean((base.F - ph @ (w * a)) ** 2) / np.mean(base.F ** 2)))
            if best is None or fit < best[0]:
                best = (fit, float(np.abs(a).max() / net.a_max))
        print(f"[{time.time() - t0:.0f}s] F best K = {K}", flush=True)
        return dict(kind='best', K=K, fit=best[0], amax=best[1])
    Th0 = net.unconstrained(rng.uniform(-0.2, 0.2, (SV_M, 1)), rng.uniform(-1.2, 1.2, (SV_M, 1)), rng.uniform(0.1, 0.6, SV_M))
    trace = []
    cb = lambda s_, T: trace.append((s_, base.fit_error(T), CO.clusters(T, G, CLUSTER_THR).max() + 1))
    Th = CO.CohesionLoop(base, SV_R0, 'gauss').run(Th0, SV_S1, DT, rec_every=8000, callback=cb)
    Th = CO.CohesionLoop(base, SV_R1, 'gauss').run(Th, SV_S2, DT, rec_every=8000, callback=lambda s_, T: cb(s_ + SV_S1, T))
    print(f"[{time.time() - t0:.0f}s] F seed {seed}", flush=True)
    return dict(kind='coh', seed=seed, end=describe(net, base, Th), trace=np.array(trace, float))


# ---------------------------------------------------------------------------- the report
def rms_windows(ts, e2, T0=PER, n_per=5):
    W = n_per * T0
    edges = np.arange(0.0, ts[-1] + 1e-9, W)
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (ts >= a) & (ts < b)
        if np.any(sel):
            out.append((a, float(np.sqrt(np.mean(e2[sel])))))
    return np.array(out)


def fmt(x, p=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return '—'
    return f'{x:.{p}g}' if abs(x) >= 1e-3 or x == 0 else f'{x:.{p - 1}e}'.replace('e-0', 'e-').replace('e+0', 'e+')


def report(A, B, C, D, E=None, F=None):
    L_ = ['# Proposition D‴.32: the cohesion leak — a leak on the unexcited subspace, without a prior', '',
          f'The three-bump instance (σ₀ = 0.15σ_crit = {SIGMA0:.6g}); the averaged slow loop with the cohesion leak dϑ_i/ds = G⟨∇(r Φ̃)⟩ − (σ₀/m)Σ_jψ(‖ϑ_i − ϑ_j‖²/R²)(ϑ_i − ϑ_j), ψ compact (1 − z)³₊ or Gaussian e^(−z/2); midpoint rule dt = {DT}.', '']
    if A is not None:
        L_ += ['## A. The linear theory', '', '### A1. Every count vector of m = 15 in W (the lift of the exact representation ϑ^e(w))', '']
        cv = A['cv']
        dmin = min(r['dmin'] for r in cv)
        for kind in ('compact', 'gauss'):
            dr = max(r[kind]['drift'] for r in cv); se = max(r[kind]['spec_err'] for r in cv)
            L_.append(f'- {kind} (R = {R_WARM:g}): {len(cv)} count vectors; the largest drift {fmt(dr)}; the Jacobian\'s spectrum against (iii) to ≤ {fmt(se)} (absolute; σ₀ = {SIGMA0:.4g}).')
        sl = [r['compact']['slowest'] for r in cv]
        L_.append(f'- The smallest inter-atom distance {dmin:.3f} (> R); the slowest rate (the weakest collective mode, the leak absent there) from {fmt(min(sl))} to {fmt(max(sl))}; the separations at σ₀w_k from {fmt(SIGMA0 * min(r["counts"].min() for r in cv) / 15)} to {fmt(SIGMA0 * max(r["counts"].max() for r in cv) / 15)}.')
        L_ += ['', '### A2. The Gaussian kernel\'s bias (the replicated stationary point near ϑ^e(w))', '', '| counts | R | Newton residual | ‖y‖ (whitened) | first order ‖Ã⁻¹b‖ | relative difference | fit error |', '|---|---|---|---|---|---|---|']
        for r in A['bias']:
            L_.append(f"| {tuple(int(v) for v in r['counts'])} | {r['R']:g} | {fmt(r['res'])} | {fmt(r['y'])} | {fmt(r['y1'])} | {fmt(r['rel'])} | {fmt(r['fit'])} |")
        S = A['split']
        L_ += ['', '### A3. The splitting margins with an approximation error (the fourth bump, ε = 0.3)', '',
               f"ϑ^b(q): gradient {fmt(S['grad'])}, fit error {fmt(S['fit_b'])}; λ⁺_k = {', '.join(fmt(v) for v in S['lam'])}.", '',
               '| counts | σ₀ factor | predicted margins σ₀w_k − (q_k/w_k)λ⁺_k | top separation eigenvalue (predicted) | top eigenvalue of the Jacobian | spectrum error | largest drift |', '|---|---|---|---|---|---|---|']
        for r in S['rows']:
            L_.append(f"| {tuple(int(v) for v in r['counts'])} | {r['frac']:g} | {', '.join(fmt(v) for v in r['margins'])} | {fmt(r['top_sep'])} | {fmt(r['top_ev'])} | {fmt(r['spec_err'])} | {fmt(r['drift'])} |")
        for run in S.get('runs', []):
            d = run['rec']
            L_.append(f"- σ₀ × {run['frac']:g}, counts {tuple(int(v) for v in run['counts'])}, jitter 10⁻³, {S_SPLIT:g} slow units: the classes' diameters at s = "
                      + ', '.join(f'{int(x[0])}' for x in d[::2]) + ': ' + '; '.join('(' + ', '.join(fmt(v) for v in x[1:4]) + ')' for x in d[::2])
                      + f"; clusters {int(d[0][5])} → {int(d[-1][5])}; fit error {fmt(d[0][4])} → {fmt(d[-1][4])}; at the end (mass, a, c, s): "
                      + '; '.join(f'({v[0]:.3f}, {v[1]:+.3f}, {v[2]:+.3f}, {v[3]:.3f})' for v in run['end']['clusters']) + '.')
        Dg = A['degen']
        L_ += ['', '### A4. Degenerate representations', '']
        for kind in ('gauss', 'compact'):
            tw = Dg[('twins', kind)]
            L_.append(f"- Twins, {kind}: the twins' distance {fmt(Dg['twins_d0'])} → " + ', '.join(f'{fmt(x[1])} (s = {int(x[0])})' for x in tw[1:]) + f"; fit error {fmt(tw[0][2])} → {fmt(tw[-1][2])}.")
        for kind in ('gauss', 'compact'):
            dd = Dg[('dead', kind)]
            L_.append(f"- Dead class, {kind}: its distance from the nearest live class {fmt(Dg['dead_d0'])} → " + ', '.join(f'{fmt(x[1])} (s = {int(x[0])})' for x in dd[1:]) + f"; fit error {fmt(dd[0][2])} → {fmt(dd[-1][2])}.")
        L_.append('')
    if B is not None:
        L_ += ['## B. The floor at warm starts (‖f̂_m − f̂_M‖² and the fit error)', '', '| leak | m | ‖f̂_m − f̂_M‖² at s = 600 | 3000 | 6000 | fit error at 6000 (median over seeds) | class diameter at the end (max) | class means from ϑ^e(w) (max) | clusters at the end |', '|---|---|---|---|---|---|---|---|---|']
        for leak in ('fixed', 'adaptive', 'coh_c', 'coh_g'):
            rs = [r for r in B if r['leak'] == leak]
            if not rs:
                continue
            s = rs[0]['s']
            for m in MS:
                vals = []
                for st in (600.0, 3000.0, 6000.0):
                    i = int(np.argmin(np.abs(s - st)))
                    v = [r['h2'][m][i] for r in rs]
                    vals.append(f'{fmt(min(v))}–{fmt(max(v))}')
                fit = np.median([r['fit'][m][-1] for r in rs])
                L_.append(f"| {leak} | {m} | {' | '.join(vals)} | {fmt(fit)} | {fmt(max(r['end'][m]['diam'] for r in rs))} | {fmt(np.nanmax([r['end'][m]['dev'] for r in rs]))} | {', '.join(str(r['end'][m]['K']) for r in rs)} |")
        L_.append('')
    if C is not None:
        L_ += ['## C. The unknown model order (broad starts)', '', '| leak | m | seed | clusters | fit error | masses | |a|/a_max (max) | centres | widths |', '|---|---|---|---|---|---|---|---|---|']
        for r in sorted(C, key=lambda r: (r['leak'], r['m'], r['seed'])):
            e = r['end']; cl = e['clusters']
            show = cl[:6]
            L_.append(f"| {r['leak']} | {r['m']} | {r['seed']} | {e['K']} | {fmt(e['fit'])} | {', '.join(f'{v:.3f}' for v in show[:, 0])}{' …' if len(cl) > 6 else ''} | {np.abs(cl[:, 1]).max() / NET_KW['a_max']:.3f} | {', '.join(f'{v:+.3f}' for v in show[:, 2])} | {', '.join(f'{v:.3f}' for v in show[:, 3])} |")
        L_.append('')
    if D is not None:
        F = D['floquet']
        L_ += ['## D. The real closed loop (k_e = 25)', '', f"### D1. The Floquet exponents of the lift (counts (6, 5, 4), compact R = {R_WARM:g})", '',
               f"The lift returns after one period to {fmt(F['ret'])}. Slow exponents (rates per unit time):", '']
        rates = F['rates']
        sep = F['sep_pred']; nsep = len(sep)
        slow = rates[rates < 1.0]
        # the separations: the exponents nearest σ₀w_k/k_e
        coll = F['coll']; ncoll = len(coll[coll < 1.0])
        L_.append(f"- {len(slow)} slow exponents ({nsep} separations predicted at σ₀w_k/k_e = {', '.join(fmt(v) for v in sorted(set(np.round(sep, 12))))}, {ncoll} collective).")
        # match: sort the slow exponents and the predictions together
        pred = np.sort(np.concatenate([sep, coll[coll < 1.0]]))
        if len(pred) == len(slow):
            rel = np.abs(slow - pred) / pred
            L_.append(f"- Against the prediction (the separations σ₀w_k/k_e and the collapsed variational equation without a leak): relative errors ≤ {fmt(rel.max())} (median {fmt(np.median(rel))}).")
        L_ += ['', '### D2. From the warm clouds (seed 0): the rms tracking error over windows of five periods', '', '| leak | m | counts | rms at t ≈ 0 | 500 | 1500 | 2900 | fit error at the end | clusters |', '|---|---|---|---|---|---|---|---|---|']
        for r in D['runs']:
            w_ = rms_windows(r['ts'], r['e2'])
            pick = [w_[int(np.argmin(np.abs(w_[:, 0] - t)))][1] for t in (0.0, 500.0, 1500.0, 2900.0)]
            L_.append(f"| {r['leak']} | {r['m']} | {tuple(int(v) for v in r['counts'])} | {' | '.join(fmt(v) for v in pick)} | {fmt(r['fit_end'])} | {r['K_end']} |")
        L_.append('')
    if E is not None:
        L_ += ['## E. Dither (τ = ' + f'{TAU_E:g}' + ', m = 60 from the lift of ϑ^e(q), three seeds)', '',
               f'Predictions: the within-class spread per coordinate √(τ/(σ₀w_k)) = {np.sqrt(TAU_E / (SIGMA0 / 3)):.4g} (w_k = 1/3) for the cohesion leak, √(τ/σ₀) = {np.sqrt(TAU_E / SIGMA0):.4g} for the leaks toward a prior, √(2τs) without a leak ({np.sqrt(2 * TAU_E * S_E / 2):.4g} at s = {S_E / 2:g}, {np.sqrt(2 * TAU_E * S_E):.4g} at s = {S_E:g}).', '',
               '| leak | spread at s = 2000 | 10000 | 20000 | class means\' deviation (whitened, median over the second half) | fit error (mean over the second half) |', '|---|---|---|---|---|---|']
        for leak in ('none', 'coh_c', 'coh_g', 'fixed', 'sigma'):
            rs = [r for r in E if r['leak'] == leak]
            if not rs:
                continue
            cols = []
            for st in (2000.0, 10000.0, 20000.0):
                v = [r['rec'][int(np.argmin(np.abs(r['rec'][:, 0] - st))), 1] for r in rs]
                cols.append(f'{fmt(min(v))}–{fmt(max(v))}')
            half = [r['rec'][r['rec'][:, 0] >= S_E / 2] for r in rs]
            dev = np.median(np.concatenate([h[:, 2] for h in half])); fit = np.mean(np.concatenate([h[:, 3] for h in half]))
            L_.append(f"| {leak} | {' | '.join(cols)} | {fmt(dev)} | {fmt(fit)} |")
        L_.append('')
    if F is not None:
        L_ += ['## F. The cogging servo (no model of the torque)', '', 'The best equal-weight K-atom fits (40 random starts each): ' +
               '; '.join(f"K = {r['K']}: {fmt(r['fit'])} (largest |a|/a_max {r['amax']:.3f})" for r in sorted([r for r in F if r['kind'] == 'best'], key=lambda r: r['K'])) +
               '. The four-kernel designs of Proposition D‴.29: ϑ^b(q) 1.46%, the saturated infimum 0.64%.', '',
               f'| seed | clusters | fit error at s = {SV_S1:g} | at the end ({SV_S1 + SV_S2:g}) | masses | largest |a|/a_max |', '|---|---|---|---|---|---|']
        for r in sorted([r for r in F if r['kind'] == 'coh'], key=lambda r: r['seed']):
            tr = r['trace']; i1 = int(np.argmin(np.abs(tr[:, 0] - SV_S1)))
            cl = r['end']['clusters']
            L_.append(f"| {r['seed']} | {r['end']['K']} | {fmt(tr[i1, 1])} | {fmt(r['end']['fit'])} | {', '.join(f'{v:.3f}' for v in cl[:, 0])} | {np.abs(cl[:, 1]).max() / SV_NET['a_max']:.3f} |")
        L_.append('')
    return '\n'.join(L_) + '\n'


# ---------------------------------------------------------------------------- driver
def save(path, **parts):
    np.savez_compressed(path, **{k: np.array([v], dtype=object) for k, v in parts.items() if v is not None})


def load(path):
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    return {k: z[k][0] for k in z.files}


def main():
    raw = os.path.join(OUT, 'cohesion_leak_raw.npz')
    md = os.path.join(OUT, 'cohesion_leak.md')
    args = sys.argv[1:]
    old = load(raw)
    if '--from-raw' in args:
        open(md, 'w').write(report(old.get('A'), old.get('B'), old.get('C'), old.get('D'), old.get('E'), old.get('F')))
        return
    parts = set('ABCDEF')
    if '--only' in args:
        parts = set(args[args.index('--only') + 1].split(','))
    from multiprocessing import Pool
    res = dict(old)
    with Pool(N_WORKERS) as pool:
        asyncs = {}
        if 'B' in parts:
            asyncs['B'] = pool.map_async(b_job, [(s, leak) for leak in ('coh_g', 'coh_c', 'fixed', 'adaptive') for s in range(N_SEEDS)])
        if 'C' in parts:
            jobs = [(m, s, leak) for m, seeds in SEEDS_C.items() for s in seeds
                    for leak in (('anneal', 'coh_c', 'none', 'sigma', 'adaptive2', 'adaptive3', 'adaptive4') if m == 60 else ('anneal',))]
            asyncs['C'] = pool.map_async(c_job, jobs)
        if 'D' in parts:
            asyncs['Druns'] = pool.map_async(d_job, [(m, leak) for m in (15, 60) for leak in ('coh', 'fixed')])
        if 'E' in parts:
            asyncs['E'] = pool.map_async(e_job, [(s, leak) for leak in ('none', 'coh_c', 'coh_g', 'fixed', 'sigma') for s in SEEDS_E])
        if 'F' in parts:
            asyncs['F'] = pool.map_async(f_job, [('coh', s) for s in (0, 1, 2)] + [('best', K) for K in (4, 6, 8)])
        if 'A' in parts:
            res['A'] = part_A()
            save(raw, **res)
        if 'D' in parts:
            res['D'] = dict(floquet=part_D_floquet(), runs=asyncs['Druns'].get())
            save(raw, **res)
        if 'B' in parts:
            res['B'] = asyncs['B'].get(); save(raw, **res)
        if 'C' in parts:
            res['C'] = asyncs['C'].get(); save(raw, **res)
        if 'E' in parts:
            res['E'] = asyncs['E'].get(); save(raw, **res)
        if 'F' in parts:
            res['F'] = asyncs['F'].get(); save(raw, **res)
    open(md, 'w').write(report(res.get('A'), res.get('B'), res.get('C'), res.get('D'), res.get('E'), res.get('F')))


if __name__ == '__main__':
    main()
