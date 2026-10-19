"""
Theorem D‴.18: Theorem J uniformly in m at a replicated configuration — the per-particle basin and the m-independent rate
=====================================================================================================================

The real closed loop (CL) of Theorem J on the three-bump instance (k_e = 25, σ = σ₀/k_e, σ₀ = 0.15σ_crit; RK4 plant with the control and the
reference evaluated at the stages and the parameters frozen over the step, Euler parameters, dt = 5·10⁻³), integrated for many systems at
once (all particles in one array, the plant states in another; each system's f̂ by a segment sum).

A. The matching prior at the matching masses q = (1/3, 1/3, 1/3): the lift is the exact representation, a constant orbit with zero tracking
   error. Sizes m = 3, 15, 60, 240; every particle displaced from its class atom by ρ in a random direction of the G⁻¹-unit sphere,
   ρ = 0.3, 1.0, 2.5 (inside the Voronoi cells, whose distance from the atoms is ≥ Δ_ν/2 = 3.13) and 4.5 (across a cell boundary for some
   particles); the same displacements under the nearest-atom leak and under the OT leak. Records: d_∞ = max_i min_k‖ϑ_i − ϑ⁰_k‖_{G⁻¹}, the
   class counts by nearest atom, W²_{G⁻¹}(μ_Θ, ν), the rms tracking error.
B. A moving orbit: the OT leak toward the prior with centres shifted by +0.15 (the residual is not zero, so the weighted 3-atom orbit Γ_q moves
   and the tracking error is not zero), masses q; the reference is the 3-atom loop started at Θ*(q), integrated alongside; sizes m = 15, 60,
   240 started at the lift of Θ*(q) displaced per particle by ρ = 0.3, 1.0; d_∞ = max_i‖ϑ_i − ϑ^{ref}_{π₀(i)}‖_{G⁻¹} against the lift of the
   reference trajectory.

The prediction of Theorem D‴.18(b): every d_∞ inside the per-particle ball decays at the rate min{λ_coll, λ_sep}/k_e — the same for every m —
with an m-independent constant.

Output: results/chaos/uniform_basin.md and uniform_basin_raw.npz (UB_T for the horizon).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import prior_design as PD
from src.core import ot_leak as OT

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
KE = 25.0
SIG = SIGMA0 / KE
DT = 5e-3
T_END = float(os.environ.get('UB_T', 8000.0))
PLAN_EVERY, REC_EVERY, E_EVERY = 200, 1000, 20
MS_A, RHOS_A = (3, 15, 60, 240), (0.3, 1.0, 2.5, 4.5)
MS_B, RHOS_B = (15, 60, 240), (0.3, 1.0)
Q = np.ones(3) / 3
A_MAX, C_MAX, S_MIN, S_MAX = 4.0, 3.0, 0.3, 2.0


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net, c_shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + c_shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def sphere(rng, n):
    v = rng.standard_normal((n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


class Batch:
    """Many real-loop systems integrated together: particles Th (P, 3) with owner (P,), plant states X (S,), per-system leak ('nearest' or
    'ot') and prior atoms (S, K, 3); the step is the scheme of exp_ot_real_loop.py."""

    def __init__(self, systems):
        self.S = len(systems)
        self.Th = np.vstack([s['Th0'] for s in systems])
        self.owner = np.concatenate([np.full(s['Th0'].shape[0], k) for k, s in enumerate(systems)])
        self.m = np.array([s['Th0'].shape[0] for s in systems])
        self.inv_m = 1.0 / self.m[self.owner]
        self.leak = [s['leak'] for s in systems]
        self.prior = np.array([s['prior'] for s in systems])
        self.X = np.full(self.S, 0.0)                               # x_m(0) = 0: e(0) = 0
        self.slices = np.concatenate([[0], np.cumsum(self.m)])
        self.T = self.targets()

    def physical(self, Th):
        a = A_MAX * np.tanh(Th[:, 0] / A_MAX)
        c = C_MAX * np.tanh(Th[:, 1] / C_MAX)
        s = S_MIN + (S_MAX - S_MIN) / (1.0 + np.exp(-np.clip(Th[:, 2], -500, 500)))
        return a, c, s

    def f_hat(self, Y, a, c, s):
        y = Y[self.owner]
        phi = np.exp(-(y - c) ** 2 / (2 * s ** 2))
        return np.bincount(self.owner, weights=self.inv_m * a * phi, minlength=self.S)

    def targets(self):
        T = np.empty_like(self.Th)
        for k in range(self.S):
            sl = slice(self.slices[k], self.slices[k + 1])
            Th = self.Th[sl]; P = self.prior[k]
            C = np.sum((Th[:, None, :] - P[None, :, :]) ** 2 / G[None, None, :], axis=2)
            if self.leak[k] == 'nearest':
                T[sl] = P[np.argmin(C, axis=1)]
            else:
                m = Th.shape[0]; rep = np.repeat(np.arange(P.shape[0]), m // P.shape[0])
                rows, cols = linear_sum_assignment(C[:, rep])
                lab = np.empty(m, int); lab[rows] = rep[cols]
                T[sl] = P[lab]
        return T

    def step(self, t, refresh):
        if refresh:
            self.T = self.targets()
        a, c, s = self.physical(self.Th)
        X = self.X
        xr = lambda tau: 2.0 * np.sin(tau); xrd = lambda tau: 2.0 * np.cos(tau)
        xdot = lambda Y, tau: f_vec(Y) + xrd(tau) - KE * (Y - xr(tau)) - self.f_hat(Y, a, c, s)
        k1 = xdot(X, t); k2 = xdot(X + 0.5 * DT * k1, t + 0.5 * DT); k3 = xdot(X + 0.5 * DT * k2, t + 0.5 * DT); k4 = xdot(X + DT * k3, t + DT)
        e = X - xr(t)
        # grad_g per particle with its system's (x, e)
        x = X[self.owner]; ee = e[self.owner]
        diff = x - c
        phi = np.exp(-diff ** 2 / (2 * s ** 2))
        da = 1.0 - (a / A_MAX) ** 2; dc = 1.0 - (c / C_MAX) ** 2
        sk = 1.0 / (1.0 + np.exp(-np.clip(self.Th[:, 2], -500, 500))); ds = (S_MAX - S_MIN) * sk * (1 - sk)
        g = np.empty_like(self.Th)
        g[:, 0] = phi * ee * da
        g[:, 1] = a * ee * phi / s ** 2 * diff * dc
        g[:, 2] = a * ee * phi * diff ** 2 / s ** 3 * ds
        self.Th = self.Th + DT * (G[None, :] * g - SIG * (self.Th - self.T))
        self.X = X + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return e

    def system(self, k):
        return self.Th[self.slices[k]:self.slices[k + 1]]


def w2_to_prior(Th, P):
    """W²_{G⁻¹}(μ_Θ, ν) for ν uniform on the rows of P (m divisible by K): the assignment problem on the replicated atoms."""
    m = Th.shape[0]; rep = np.repeat(np.arange(P.shape[0]), m // P.shape[0])
    C = np.sum((Th[:, None, :] - P[None, rep, :]) ** 2 / G[None, None, :], axis=2)
    r, cidx = linear_sum_assignment(C)
    return float(C[r, cidx].mean())


def run_batch(part, systems, meta):
    B = Batch(systems)
    n = int(round(T_END / DT))
    rec = {k: dict(t=[], dinf=[], counts=[], w2=[]) for k in range(B.S)}
    e2acc = np.zeros(B.S); e2n = 0
    t0 = time.time()
    for i in range(n):
        t = i * DT
        if i % REC_EVERY == 0:
            for k in range(B.S):
                Th = B.system(k); P = B.prior[k]
                d = np.sqrt(np.sum((Th[:, None, :] - P[None, :, :]) ** 2 / G[None, None, :], axis=2))
                lab = np.argmin(d, axis=1)
                rec[k]['t'].append(t); rec[k]['counts'].append(np.bincount(lab, minlength=P.shape[0]))
                if part == 'A':
                    rec[k]['dinf'].append(float(d.min(axis=1).max()))
                    rec[k]['w2'].append(w2_to_prior(Th, P))
                else:
                    ref = B.system(meta['ref_index'])
                    lab0 = systems[k]['labels0']
                    rec[k]['dinf'].append(float(np.sqrt(np.sum((Th - ref[lab0]) ** 2 / G[None, :], axis=1)).max()))
                    rec[k]['w2'].append(float('nan'))
        e = B.step(t, i % PLAN_EVERY == 0)
        if i % E_EVERY == 0 and t > T_END - 10 * PER:
            e2acc += e ** 2; e2n += 1
    out = []
    for k, s in enumerate(systems):
        r = {kk: np.array(v) for kk, v in rec[k].items()}
        r.update(dict(leak=s['leak'], m=int(B.m[k]), rho=s['rho'], seed=s['seed'], e_rms=float(np.sqrt(e2acc[k] / max(e2n, 1))), final=B.system(k).copy()))
        out.append(r)
    print(f"[part {part}] {B.S} systems, {B.Th.shape[0]} particles, {time.time() - t0:.0f}s", flush=True)
    return out


def predictions(net, prior, w=Q):
    """λ_sep = min_k λ_min(σ₀I − R^G_k) and λ_coll = λ_min(H_w) at the weighted stationary point Θ*(w) of the nearest-atom loop toward `prior`."""
    f_x = lambda X: f_vec(X)
    loop = PD.MixturePriorLoop(net, f_x, lambda t: np.array([2.0 * np.sin(t)]), PER, G, SIGMA0, prior, n_quad=256)
    atoms, res = PD.weighted_stationary_atoms_prior(loop, w, prior)
    a, c, s, diff, phi = loop.features(atoms); r = loop.F - phi @ (w * a)
    R = loop.residual_curvature(r, atoms); Rw = SQG[None, :, None] * R * SQG[None, None, :]
    lam_sep = min(float(np.linalg.eigvalsh(SIGMA0 * np.eye(3) - Rw[k]).min()) for k in range(3))
    K = 3; d = 3 * K; v = atoms.ravel(); Jm = np.zeros((d, d)); h = 1e-6
    gf = lambda V: PD.weighted_energy_grad_prior(loop, w, V.reshape(K, 3)).ravel()
    for j in range(d):
        ej = np.zeros(d); ej[j] = h
        Jm[:, j] = (gf(v + ej) - gf(v - ej)) / (2 * h)
    sw = np.repeat(np.sqrt(w), 3); sg = np.tile(SQG, K)
    Hs = (sw * sg)[:, None] * Jm * (sg / sw)[None, :]; Hs = 0.5 * (Hs + Hs.T)
    lam_coll = float(np.linalg.eigvalsh(Hs).min())
    return atoms, res, lam_sep, lam_coll


def job(part):
    net = make_net()
    if part == 'A':
        prior = exact_representation(net)
        systems = []
        for m in MS_A:
            for rho in RHOS_A:
                seed = 1000 * m + int(10 * rho)
                rng = np.random.default_rng(seed)
                Th0 = np.repeat(prior, m // 3, axis=0) + rho * sphere(rng, m) * SQG[None, :]
                for leak in ('nearest', 'ot'):
                    systems.append(dict(Th0=Th0.copy(), leak=leak, prior=prior, rho=rho, seed=seed))
        return dict(part='A', runs=run_batch('A', systems, {}), prior=prior)
    prior = exact_representation(net, c_shift=0.15)
    atoms, res, lam_sep, lam_coll = predictions(net, prior)
    systems = [dict(Th0=atoms.copy(), leak='ot', prior=prior, rho=0.0, seed=0, labels0=np.arange(3))]
    for m in MS_B:
        for rho in RHOS_B:
            seed = 2000 * m + int(10 * rho)
            rng = np.random.default_rng(seed)
            lab0 = np.repeat(np.arange(3), m // 3)
            systems.append(dict(Th0=atoms[lab0] + rho * sphere(rng, m) * SQG[None, :], leak='ot', prior=prior, rho=rho, seed=seed, labels0=lab0))
    return dict(part='B', runs=run_batch('B', systems, dict(ref_index=0)), prior=prior, atoms=atoms, res=res, lam_sep=lam_sep, lam_coll=lam_coll)


def decay_rate(t, d, frac=0.5):
    msk = (t > frac * t[-1]) & (d > 0)
    if msk.sum() < 5:
        return float('nan')
    return float(-np.polyfit(t[msk], np.log(d[msk]), 1)[0])


def report(A, B, t0=None):
    net = make_net()
    priorA = exact_representation(net)
    _, _, lsepA, lcollA = predictions(net, priorA)
    muA = min(lsepA, lcollA) / KE; muB = min(B['lam_sep'], B['lam_coll']) / KE
    lines = ["# Theorem D‴.18: Theorem J uniformly in m at a replicated configuration — the per-particle basin and the m-independent rate (`exp_uniform_basin.py`)", "",
             f"The real closed loop of Theorem J on the three-bump instance: k_e = {KE:g}, σ = σ₀/k_e with σ₀ = 0.15σ_crit = {SIGMA0:.5f}, RK4 plant with the control and the reference at the stages and the parameters "
             f"frozen over the step, Euler parameters, dt = {DT:g}, horizon {T_END:g} time units ({T_END / PER:.0f} periods, {T_END / KE:.0f} slow units); the leak target refreshed every {PLAN_EVERY} steps; records every "
             f"{REC_EVERY * DT:g} time units. Fitted rates: least-squares slope of log d_∞ over the second half of the run.", "",
             f"## A. The matching prior at the matching masses (the lift is the exact representation, a constant orbit)", "",
             f"Prediction: λ_sep = {lsepA:.5f} (= σ₀ at zero residual), λ_coll = {lcollA:.5f}, rate min{{λ_coll, λ_sep}}/k_e = {muA:.3e} per time unit. Every particle is displaced from its class atom by ρ "
             "in a random direction of the G⁻¹-unit sphere (the same displacements under both leaks); the Voronoi cells extend at least Δ_ν/2 = 3.13 from the atoms.", "",
             "| leak | ρ | m | final counts by nearest atom | d_∞(0) → d_∞(end) | fitted rate of d_∞ (second half) | C = max_t d_∞(t)e^{μt}/d_∞(0) (μ the predicted rate) | W²(μ_Θ, ν): start → end | rms tracking error (last 10 periods) |",
             "|---|---|---|---|---|---|---|---|---|"]
    for leak in ('nearest', 'ot'):
        for rho in RHOS_A:
            for m in MS_A:
                r = next(x for x in A['runs'] if x['leak'] == leak and x['rho'] == rho and x['m'] == m)
                rate = decay_rate(r['t'], r['dinf'])
                Cc = float(np.max(r['dinf'] * np.exp(muA * r['t'])) / r['dinf'][0])
                lines.append(f"| {leak} | {rho:g} | {m} | {tuple(int(x) for x in r['counts'][-1])} | {r['dinf'][0]:.3f} → {r['dinf'][-1]:.2e} | {rate:.3e} | {Cc:.2f} | {r['w2'][0]:.3e} → {r['w2'][-1]:.2e} | {r['e_rms']:.1e} |")
    inside = [x for x in A['runs'] if x['rho'] < 3.13]
    rates = [decay_rate(x['t'], x['dinf']) for x in inside]
    Cs = [float(np.max(x['dinf'] * np.exp(muA * x['t'])) / x['dinf'][0]) for x in inside]
    out45 = {(x['leak'], x['m']): x for x in A['runs'] if x['rho'] == 4.5}
    lines += ["", f"Inside the cells (ρ = 0.3, 1.0, 2.5; both leaks; all m): fitted rates {min(rates):.3e}–{max(rates):.3e} against the predicted {muA:.3e} (ratios {min(rates) / muA:.3f}–{max(rates) / muA:.3f}); "
              f"constants C = {min(Cs):.2f}–{max(Cs):.2f}. At ρ = 4.5: the nearest-atom leak ends with counts " + ", ".join(f"{tuple(int(v) for v in out45[('nearest', m)]['counts'][-1])} (m = {m})" for m in MS_A)
              + " and the OT leak with " + ", ".join(f"{tuple(int(v) for v in out45[('ot', m)]['counts'][-1])}" for m in MS_A) + "; W²(μ_Θ, ν) at the end " + ", ".join(f"{out45[('nearest', m)]['w2'][-1]:.2e} / {out45[('ot', m)]['w2'][-1]:.2e}" for m in MS_A) + " (nearest / OT)."]
    lines += ["", "## B. A moving orbit: the OT leak toward the prior with centres shifted by +0.15 (non-zero residual), masses q", "",
              f"Θ*(q) by Newton (residual {B['res']:.1e}); prediction λ_sep = {B['lam_sep']:.5f}, λ_coll = {B['lam_coll']:.5f}, rate {muB:.3e} per time unit. The reference is the 3-atom loop from Θ*(q); d_∞ is measured "
              "against the lift of the reference trajectory.", "",
              "| ρ | m | d_∞(0) → d_∞(end) | fitted rate | C | rms tracking error (m-system / reference) |", "|---|---|---|---|---|---|"]
    ref = B['runs'][0]
    for r in B['runs'][1:]:
        rate = decay_rate(r['t'], r['dinf'])
        Cc = float(np.max(r['dinf'] * np.exp(muB * r['t'])) / r['dinf'][0])
        lines.append(f"| {r['rho']:g} | {r['m']} | {r['dinf'][0]:.3f} → {r['dinf'][-1]:.2e} | {rate:.3e} | {Cc:.2f} | {r['e_rms']:.2e} / {ref['e_rms']:.2e} |")
    ratesB = [decay_rate(r['t'], r['dinf']) for r in B['runs'][1:]]
    lines += ["", f"Fitted rates {min(ratesB):.3e}–{max(ratesB):.3e} against the predicted {muB:.3e} (ratios {min(ratesB) / muB:.3f}–{max(ratesB) / muB:.3f}).", "",
              "Numbers are produced by `experiments/exp_uniform_basin.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'uniform_basin.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    if '--from-raw' in sys.argv:
        z = np.load(os.path.join(OUT, 'uniform_basin_raw.npz'), allow_pickle=True)
        report(z['A'].item(), z['B'].item(), t0)
        return
    with mp.Pool(2) as pool:
        A, B = pool.map(job, ['A', 'B'])
    for part in (A, B):
        for r in part['runs']:
            r['final'] = r['final'].astype(np.float32)
    np.savez(os.path.join(OUT, 'uniform_basin_raw.npz'), A=np.array(A, dtype=object), B=np.array(B, dtype=object), allow_pickle=True)
    report(A, B, t0)


if __name__ == '__main__':
    main()
