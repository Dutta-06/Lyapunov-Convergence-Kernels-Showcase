"""
Proposition D‴.20: the entropic OT leak — a soft mass constraint with a critical excess
======================================================================================

The averaged slow loop of Theorem D‴ on the three-bump instance (σ₀ = 0.15σ_crit), with the entropic OT leak toward the matching prior
(three atoms at the exact representation, masses q = 1/3) at regularisations ε, against the hard OT leak (ε → 0) of Proposition D‴.12.

1. The critical excess (no dynamics). For the weights w_δ = q + δ(1.5e₁ − 0.5) the weighted stationary point under the entropic leak, the
   separation block of the over-weighted class σ₀(I − (2/ε)Cov^G₁) − R^G₁ (Cov^G₁ the whitened covariance of the atoms under the class's
   conditional plan) and its smallest eigenvalue, against ε; the critical ε_c(δ) by bisection, and the rank-one prediction
   ε_c ≈ 2σ₀p(1 − p)Δ̂ᵀ(σ₀I − R^G₁)⁻¹Δ̂ with p the excess fraction sent away and Δ̂ the whitened displacement to the receiving atom.
2. Dynamics. m = 15, 60, 240 against M = 510 from the warm clouds of `exp_ot_leak.py` (the same seeds), ε = 1, 4, 16, 64, 600 slow units:
   the cluster weights by nearest atom (equalised or not), the fit error, W²_{G⁻¹}(μ_m, ν) and the function-level distance to the reference.

Output: results/chaos/entropic_leak.md and entropic_leak_raw.npz (EL_S_END, EL_SEEDS, CLW_WORKERS).
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
from src.core import entropic_leak as EL
from src.core import ot_leak as OT

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q = np.ones(3) / 3
EPSS = (1.0, 4.0, 16.0, 64.0)
DELTAS = (0.02, 0.05, 0.1, 0.2)
MS, M_REF = (15, 60, 240), 510
S_END, DT_S, REC, PLAN_EVERY = float(os.environ.get('EL_S_END', 600.0)), 0.02, 50, 10
N_SEEDS = int(os.environ.get('EL_SEEDS', 4))
N_WORKERS = int(os.environ.get('CLW_WORKERS', 3))
JITTER = 0.02


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def make_loop(net, eps):
    return EL.EntropicPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact_representation(net), eps, n_quad=N_Q)


# ---------------------------------------------------------------------------- 1. the critical excess
def critical_row(delta):
    net = make_net(); exact = exact_representation(net)
    w = Q + delta * (1.5 * np.eye(3)[0] - 0.5)
    def margin(eps):
        lp = make_loop(net, eps)
        atoms, res = EL.weighted_stationary_atoms_entropic(lp, w, exact, eps, Q)
        blocks, mins, cov = EL.separation_blocks_entropic(lp, w, atoms, eps, Q)
        return mins[0], atoms, res, blocks, cov, lp
    grid = []
    for eps in (1.0, 4.0, 16.0, 64.0, 256.0):
        mg, atoms, res, blocks, cov, lp = margin(eps)
        grid.append((eps, mg, res))
    # bisection in log ε between the last negative and the first positive grid point
    lo = max([e for e, mg, r in grid if mg < 0], default=None); hi = min([e for e, mg, r in grid if mg > 0 and (lo is None or e > lo)], default=None)
    eps_c = float('nan'); pred = float('nan'); p = float('nan')
    if lo is not None and hi is not None:
        a, b = np.log(lo), np.log(hi)
        for _ in range(30):
            c = 0.5 * (a + b)
            if margin(np.exp(c))[0] < 0:
                a = c
            else:
                b = c
        eps_c = float(np.exp(0.5 * (a + b)))
        mg, atoms, res, blocks, cov, lp = margin(eps_c)
        # the rank-one prediction at ε_c: the excess fraction p sent away and the receiving atom from the class's conditional plan
        pi, T, covs = EL.weighted_targets(lp.prior, w, atoms, Q, eps_c, G)
        cond = pi[0] / pi[0].sum()
        p = float(1.0 - cond[0])
        j = int(np.argmax(cond[1:]) + 1)
        dh = (exact[j] - exact[0]) / SQG
        a_, c_, s_, diff_, phi_ = lp.features(atoms); r = lp.F - phi_ @ (w * a_)
        Rw = SQG[:, None] * lp.residual_curvature(r, atoms)[0] * SQG[None, :]
        A = SIGMA0 * np.eye(3) - Rw
        pred = float(2 * SIGMA0 * p * (1 - p) * dh @ np.linalg.solve(A, dh))
    return dict(delta=delta, w=w, grid=grid, eps_c=eps_c, p=p, pred=pred)


# ---------------------------------------------------------------------------- 2. dynamics
def job(args):
    seed, eps = args
    net = make_net(); exact = exact_representation(net)
    rng = np.random.default_rng(seed)
    Th_all = CL.mixture_start(exact, M_REF, rng, JITTER, G)                 # the clouds of exp_ot_leak.py (same seeds, same generator calls)
    sysm = {m: Th_all[:m].copy() for m in MS + (M_REF,)}
    loops = {m: make_loop(net, eps) for m in sysm}
    T = {}
    n = int(round(S_END / DT_S))
    rec = dict(s=[], w={m: [] for m in sysm}, fit={m: [] for m in sysm}, W2nu={m: [] for m in sysm}, h2={m: [] for m in MS})
    t0 = time.time()
    for k in range(n):
        if k % PLAN_EVERY == 0:
            for m, Th in sysm.items():
                T[m] = loops[m].leak_target(Th)
        new = {}
        for m, Th in sysm.items():
            lp = loops[m]
            k1 = lp.drift(Th, T[m]); k2 = lp.drift(Th + 0.5 * DT_S * k1, T[m]); new[m] = Th + DT_S * k2
        sysm = new
        if k % REC == 0 or k == n - 1:
            rec['s'].append((k + 1) * DT_S)
            fM = loops[M_REF].f_hat(sysm[M_REF])
            for m, Th in sysm.items():
                rec['w'][m].append(CL.cluster_weights(Th, exact, G)[1])
                rec['fit'][m].append(loops[m].fit_error(Th))
                rec['W2nu'][m].append(_w2nu(Th, exact))
            for m in MS:
                rec['h2'][m].append(float(np.mean((loops[m].f_hat(sysm[m]) - fM) ** 2)))
    out = dict(seed=seed, eps=eps, s=np.array(rec['s']))
    for key in ('w', 'fit', 'W2nu', 'h2'):
        out[key] = {m: np.array(v) for m, v in rec[key].items()}
    out['final'] = {m: Th.astype(np.float32) for m, Th in sysm.items()}
    print(f"seed {seed}, ε = {eps}: {time.time() - t0:.0f}s", flush=True)
    return out


def _w2nu(Th, prior):
    """W²_{G⁻¹}(μ_Θ, ν) for ν uniform on the three prior atoms (the hard transport cost, a diagnostic for every leak)."""
    m = Th.shape[0]
    C = OT._cost_matrix(Th, prior, G)
    if m % 3 == 0:
        rep = np.repeat(np.arange(3), m // 3)
        from scipy.optimize import linear_sum_assignment
        r, c = linear_sum_assignment(C[:, rep]); return float(C[r, rep[c]].mean())
    gamma = OT.transport_plan(C, np.full(m, 1.0 / m), Q); return float(np.sum(gamma * C))


def report(crit, runs, t0=None):
    lines = ["# Proposition D‴.20: the entropic OT leak — a soft mass constraint with a critical excess (`exp_entropic_leak.py`)", "",
             f"The three-bump instance with the matching prior (masses q = 1/3) at σ₀ = 0.15σ_crit = {SIGMA0:.5f}; the entropic plan by log-domain Sinkhorn (tolerance 10⁻¹³).", "",
             "## 1. The critical excess: the separation margin of the over-weighted class of Θ*(w_δ), w_δ = q + δ(1.5e₁ − 0.5), against ε", "",
             "| δ | w_δ | smallest eigenvalue of the separation block of class 1 at ε = 1, 4, 16, 64, 256 | ε_c (bisection) | excess fraction p sent away at ε_c | rank-one prediction 2σ₀p(1 − p)Δ̂ᵀ(σ₀I − R^G)⁻¹Δ̂ at ε_c |",
             "|---|---|---|---|---|---|"]
    for c in crit:
        lines.append(f"| {c['delta']:g} | ({', '.join(f'{x:.4f}' for x in c['w'])}) | " + ", ".join(f"{mg:+.5f}" for e, mg, r in c['grid']) + f" | {c['eps_c']:.3f} | {c['p']:.4f} | {c['pred']:.3f} |")
    lines += ["", "Below ε_c the over-weighted collapsed class is a saddle — it splits, and the split-off particles carry the excess to the under-weighted atom (the hard OT leak's transit); above ε_c it is stable and the "
              "excess is absorbed by a bias of the class atom toward its target T₁ = Σ_j(π_{1j}/w₁)ϑ⁰_j."]
    lines += ["", f"## 2. Dynamics: m = {', '.join(str(m) for m in MS)} against M = {M_REF} from the warm clouds of `exp_ot_leak.py`, {S_END:g} slow units, {N_SEEDS} seeds", "",
              "| ε | m | cluster weights at the end: max |w − q|_∞ over seeds (start) | fit error at the end (mean over seeds) | W²(μ_m, ν) at the end (mean) | ‖f̂_m − f̂_M‖² on the last 150 slow units (mean) |",
              "|---|---|---|---|---|---|"]
    for eps in EPSS:
        R = [r for r in runs if r['eps'] == eps]
        for m in MS + (M_REF,):
            w0 = max(float(np.abs(r['w'][m][0] - Q).max()) for r in R); w1 = max(float(np.abs(r['w'][m][-1] - Q).max()) for r in R)
            fit = np.mean([r['fit'][m][-1] for r in R]); W2 = np.mean([r['W2nu'][m][-1] for r in R])
            h2 = np.mean([float(np.mean(r['h2'][m][r['s'] > S_END - 150])) for r in R]) if m in MS else float('nan')
            lines.append(f"| {eps:g} | {m}{' (reference)' if m == M_REF else ''} | {w1:.4f} ({w0:.4f}) | {fit:.3e} | {W2:.3e} | " + (f"{h2:.3e}" if np.isfinite(h2) else "—") + " |")
    lines += ["", "Numbers are produced by `experiments/exp_entropic_leak.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'entropic_leak.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    if '--from-raw' in sys.argv:
        z = np.load(os.path.join(OUT, 'entropic_leak_raw.npz'), allow_pickle=True)
        report(list(z['crit']), list(z['runs']), t0)
        return
    crit = [critical_row(d) for d in DELTAS]
    print(f"[critical excess done, {time.time() - t0:.0f}s]", flush=True)
    jobs = [(seed, eps) for eps in EPSS for seed in range(N_SEEDS)]
    with mp.Pool(N_WORKERS) as pool:
        runs = pool.map(job, jobs, chunksize=1)
    np.savez(os.path.join(OUT, 'entropic_leak_raw.npz'), crit=np.array(crit, dtype=object), runs=np.array(runs, dtype=object), allow_pickle=True)
    report(crit, runs, t0)


if __name__ == '__main__':
    main()
