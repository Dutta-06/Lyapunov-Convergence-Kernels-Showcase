"""
Proposition D‴.34: mass transfer — the replicator weight law with the cohesion leak
===================================================================================

The weighted network f̂ = (1/m)Σβ_iΦ̃(·; ϑ_i), Σβ = m, the replicator law dβ_i/ds = γβ_i(u_i − λ) beside the transport law with the cohesion
leak (Proposition D‴.32), on the instances of `exp_cohesion_leak.py` (σ₀ = 0.15σ_crit, midpoint rule dt = 0.25 on (Θ, log β)):

A. Fixed positions (Proposition D‴.34(iv)): the five-bump target's classes at their atoms (K = 5, one particle per class, the masses
   (0.15, 0.30, 0.20, 0.25, 0.10) short of representability in two classes), the amplitudes adapted (the centres and widths frozen), γ = 0
   (transport alone) and γ > 0; and the target with its amplitudes scaled by 1.25 (Σ_k|A_k|/a_max = 1.09: no mass vector fits).
B. The five-bump target from broad starts (Part H of `exp_cohesion_leak.py`, where the transit's masses were not representable): the
   annealed Gaussian cohesion leak (R = 2, then 1) with γ = 0 (Proposition D‴.32's loop) and γ > 0, three seeds.
C. The resolution's finest stage (Part G: the fourth bump, compact R = 1, σ₀/64): from the jittered lift of ϑ^b(q) (counts 10, 10, 10)
   with γ = 0 and γ > 0, three seeds.
D. The instance from broad starts (Part C of `exp_cohesion_leak.py`): the annealed schedule with γ > 0, three seeds.
E. The real loop at the instance's atoms (m = 3, the replicator on the weights from β(0) = (1.6, 0.7, 0.7)): the relative-entropy Lyapunov
   function V = ½e² + (1/(γm))Σβ*log(β*/β) against its derivative −k_ee², Σβ, and the weights' limit; and the linearisation of the averaged
   loop (Θ, log β) at the lift of the exact representation (counts (6, 5, 4), compact R = 1): its zero modes and the rest of its spectrum
   against Proposition D‴.32's.

Output: results/chaos/mass_transfer.md and mass_transfer_raw.npz (MT_WORKERS; `--from-raw` regenerates the md; `--only A,B,…` recomputes parts).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import langevin as L
from src.core import cohesion_leak as CO
from src.core import mass_transfer as MT
from src.core.chaos_learning import AveragedSlowLoop
import exp_cohesion_leak as E                                           # the instances of Proposition D‴.32
import exp_measurement_noise as XM                                       # the real loop's instance (f, x_m, the exact representation of q = 1/3)

OUT = os.environ.get('MT_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('MT_WORKERS', 3))
GAMMAS = (0.0, 1e-3, 1e-2)
SEEDS = (0, 1, 2)
S_A, S_C = 40000.0, 100000.0
W0_A = np.array([0.15, 0.30, 0.20, 0.25, 0.10])
GAMMA_E, KE_E, T_E = 2.0, 25.0, 200.0


def target_loop(bumps, s_min=None):
    fv = lambda X: sum(A * np.exp(-(np.asarray(X, float) - C) ** 2 / (2 * S ** 2)) for A, C, S in bumps)
    net = L.SaturatedLiftedNetwork(1, **dict(E.NET_KW, s_min=E.NET_KW['s_min'] if s_min is None else s_min))
    return net, AveragedSlowLoop(net, fv, E.x_ref, E.PER, E.G, E.SIGMA0, n_quad=E.N_Q)


def summary(net, loop, Th, beta):
    """The clusters (single linkage at 0.25) by centre: masses Σβ/m, β-weighted atoms (a, c, s), particle counts; the fit error."""
    lab, W, A = MT.class_masses(Th, beta, E.G, E.CLUSTER_THR)
    a, c, s = net.physical(A)
    order = np.argsort(c[:, 0])
    return dict(K=int(len(W)), K1=int(np.sum(W > 0.01)), masses=W[order], a=a[order, 0], c=c[order, 0], s=s[order],
                n=np.array([np.sum(lab == k) for k in order]), fit=loop.fit_error(Th, beta), beta_min=float(beta.min()), beta_max=float(beta.max()))


class Frozen(MT.ReplicatorCohesionLoop):
    """The positions' centres and widths frozen (the amplitudes and the weights adapt)."""

    def drift(self, Th, beta):
        dT, db = super().drift(Th, beta)
        dT = dT.copy(); dT[:, 1:] = 0.0
        return dT, db


# ---------------------------------------------------------------------------- A. fixed positions
def a_job(args):
    gamma, scale = args
    t0 = time.time()
    bumps = [(scale * A, C, S) for A, C, S in E.TARGETS_H['five']]
    net, base = target_loop(bumps, s_min=0.2)
    K = len(bumps)
    beta = K * W0_A
    Th = net.unconstrained(np.full((K, 1), 0.5), np.array([[C] for A, C, S in bumps]), np.array([S for A, C, S in bumps]))
    fl = Frozen(base, 1.0, 'compact', sigma0=0.0, gamma=gamma)
    rec = []
    Th, beta = fl.run(Th, beta, S_A, E.DT, rec_every=4000,
                      callback=lambda s_, T, b: rec.append([s_, fl.fit_error(T, b)] + list(b / K) + list(np.abs(net.physical(T)[0][:, 0]) / net.a_max)))
    print(f"[{time.time() - t0:.0f}s] A gamma {gamma:g} scale {scale:g}", flush=True)
    return dict(gamma=gamma, scale=scale, need=np.array([abs(A) / net.a_max for A, C, S in bumps]), rec=np.array(rec),
                end=summary(net, fl, Th, beta))


# ---------------------------------------------------------------------------- B. five bumps, D. the instance (broad starts, annealed)
def broad_job(args):
    part, gamma, seed = args[:3]
    name = args[3] if len(args) > 3 else None
    t0 = time.time()
    if part == 'B':
        net, base = target_loop(E.TARGETS_H['five'], s_min=0.2)
    elif part == 'F':
        net, base = target_loop(E.TARGETS_K[name])
    else:
        net, base = target_loop(E.BUMPS)
    Th, beta = E.broad(net, E.M_H, seed), np.ones(E.M_H)
    rec = []
    fitf = MT.ReplicatorCohesionLoop(base, 1.0).fit_error
    cb = lambda s_, T, b: rec.append((s_, fitf(T, b), len(MT.class_masses(T, b, E.G, E.CLUSTER_THR)[1])))
    Th, beta = MT.ReplicatorCohesionLoop(base, E.R_BROAD, 'gauss', gamma=gamma).run(Th, beta, E.S_C1, E.DT, rec_every=4000, callback=cb)
    mid = summary(net, MT.ReplicatorCohesionLoop(base, E.R_BROAD, 'gauss', gamma=gamma), Th, beta)
    loop = MT.ReplicatorCohesionLoop(base, E.R_FINE, 'gauss', gamma=gamma)
    Th, beta = loop.run(Th, beta, E.S_C2, E.DT, rec_every=4000, callback=lambda s_, T, b: cb(s_ + E.S_C1, T, b))
    print(f"[{time.time() - t0:.0f}s] {part} {name or ''} gamma {gamma:g} seed {seed}", flush=True)
    return dict(part=part, name=name, gamma=gamma, seed=seed, mid=dict(K=mid['K'], fit=mid['fit']), end=summary(net, loop, Th, beta), rec=np.array(rec))


# ---------------------------------------------------------------------------- C. the finest resolution stage
def c_job(args):
    gamma, seed = args
    t0 = time.time()
    net, base, bq = E.approximant()
    counts = np.array([10, 10, 10]); w = counts / 30
    rng = np.random.default_rng(700 + seed)
    Th = np.repeat(E.gauge(net, bq, w), counts, axis=0) + 1e-3 * rng.standard_normal((30, 3)) * np.sqrt(E.G)[None, :]
    loop = MT.ReplicatorCohesionLoop(base, E.R_WARM, 'compact', sigma0=E.SIGMA0 / 64, gamma=gamma)
    rec = []
    Th, beta = loop.run(Th, np.ones(30), S_C, E.DT, rec_every=8000,
                        callback=lambda s_, T, b: rec.append((s_, loop.fit_error(T, b), len(MT.class_masses(T, b, E.G, E.CLUSTER_THR)[1]))))
    print(f"[{time.time() - t0:.0f}s] C gamma {gamma:g} seed {seed}", flush=True)
    return dict(gamma=gamma, seed=seed, end=summary(net, loop, Th, beta), rec=np.array(rec))


# ---------------------------------------------------------------------------- E. the real loop and the linearisation
def part_E():
    t0 = time.time()
    net, atoms = XM.setup()
    m = 3
    bstar = np.ones(m)
    z0 = np.concatenate([XM.xr(0.0) + 0.05, [1.6, 0.7, 0.7]])
    rhs = lambda t, z: MT.replicator_real_rhs(t, z, net, XM.f_vec, XM.xr, XM.xrd, KE_E, atoms, GAMMA_E, m)
    ts = np.linspace(0.0, T_E, 8001)
    sol = solve_ivp(rhs, (0.0, T_E), z0, method='DOP853', rtol=1e-12, atol=1e-14, t_eval=ts)
    V = np.array([MT.kl_lyapunov(sol.y[:, i], XM.xr, t, bstar, GAMMA_E, m) for i, t in enumerate(ts)])
    e = sol.y[0] - XM.xr(ts)[0]
    dV = np.empty_like(ts)                                                # V̇ along the solution from the vector field itself
    for i, t in enumerate(ts):
        zd = rhs(t, sol.y[:, i])
        dV[i] = e[i] * (zd[0] - XM.xrd(t)[0]) - float(np.sum(bstar * zd[1:] / sol.y[1:, i])) / (GAMMA_E * m)
    real = dict(V0=float(V[0]), Vend=float(V[-1]), monotone=bool(np.all(np.diff(V) <= 1e-13)),
                dVerr=float(np.abs(dV + KE_E * e ** 2).max() / np.abs(KE_E * e ** 2).max()),
                sum_err=float(np.abs(sol.y[1:].sum(axis=0) - m).max()), beta_end=sol.y[1:, -1].tolist())
    # the linearisation of the averaged loop at the lift of ϑ^e(w), counts (6, 5, 4), compact R = 1
    netc, base = E.setup()
    counts = np.array([6, 5, 4]); w = counts / 15
    Th0 = np.repeat(E.exact_rep(netc, w), counts, axis=0)
    out = {}
    for gamma in (0.0, 1e-2):
        loop = MT.ReplicatorCohesionLoop(base, E.R_WARM, 'compact', gamma=gamma)
        x0 = np.concatenate([Th0.ravel(), np.zeros(15)])

        def F(x):
            T = x[:45].reshape(15, 3); b = np.exp(x[45:])
            dT, db = loop.drift(T, b)
            return np.concatenate([dT.ravel(), db / b])
        Jm = np.zeros((60, 60)); h = 1e-6
        for j in range(60):
            dx = np.zeros(60); dx[j] = h
            Jm[:, j] = (F(x0 + dx) - F(x0 - dx)) / (2 * h)
        ev = np.sort(np.real(np.linalg.eigvals(Jm)))
        out[gamma] = dict(ev=ev, drift=float(np.abs(F(x0)).max()))
    tol = 1e-9
    n0 = {g: int(np.sum(np.abs(o['ev']) < tol)) for g, o in out.items()}
    pred = np.sort(np.concatenate(CO.predicted_spectrum(base, E.exact_rep(netc, w), counts, E.SIGMA0)))
    print(f"[{time.time() - t0:.0f}s] E", flush=True)
    return dict(real=real, lin={g: dict(ev=o['ev'], drift=o['drift'], zero=n0[g]) for g, o in out.items()}, pred=pred)


def report(A=None, B=None, C=None, D=None, Ep=None, F=None):
    fmt = E.fmt
    L_ = ["# Proposition D‴.34: mass transfer — the replicator weight law with the cohesion leak (`exp_mass_transfer.py`)", "",
          f"The instances of `exp_cohesion_leak.py` (σ₀ = {E.SIGMA0:.6g}); the weighted network f̂ = (1/m)Σβ_iΦ̃(·;ϑ_i), Σβ = m, the replicator law "
          "dβ_i/ds = γβ_i(u_i − λ), u_i = ⟨rΦ̃(·;ϑ_i)⟩, beside the transport law with the cohesion leak (γ = 0: Proposition D‴.32's loop); masses are Σβ/m per "
          "cluster (single linkage at 0.25), the classes listed by centre.", ""]
    if A is not None:
        L_ += ["## A. Fixed positions (the five-bump atoms, one particle per class, the amplitudes adapted; 40000 slow units)", "",
               "| γ | amplitude scale | Σ|A_k|/a_max | fit error at s = 4000, 20000, 40000 | masses at the end | needed |A_k|/a_max | |a|/a_max at the end |",
               "|---|---|---|---|---|---|---|"]
        for r in sorted(A, key=lambda r: (r['scale'], r['gamma'])):
            e = r['end']; rec = r['rec']
            pick = [rec[int(np.argmin(np.abs(rec[:, 0] - s_))), 1] for s_ in (4000.0, 20000.0, 40000.0)]
            L_.append(f"| {r['gamma']:g} | {r['scale']:g} | {r['need'].sum():.3f} | {', '.join(fmt(v) for v in pick)} | {', '.join(f'{v:.3f}' for v in e['masses'])} | "
                      f"{', '.join(f'{v:.3f}' for v in r['need'])} | {', '.join(f'{abs(v) / 4:.3f}' for v in e['a'])} |")
        L_.append("")
    for part, rows, title in (('B', B, 'The five-bump target from broad starts (m = 60, the annealed Gaussian cohesion leak)'),
                              ('D', D, 'The instance from broad starts (m = 60, the annealed Gaussian cohesion leak)')):
        if rows is None:
            continue
        L_ += [f"## {part}. {title}", "",
               "| γ | seed | clusters after R = 2 | fit error | clusters at the end (> 1% of the mass) | fit error | masses | centres | particles | |a|/a_max (max) |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(rows, key=lambda r: (r['gamma'], r['seed'])):
            e = r['end']
            L_.append(f"| {r['gamma']:g} | {r['seed']} | {r['mid']['K']} | {fmt(r['mid']['fit'])} | {e['K']} ({e['K1']}) | {fmt(e['fit'])} | "
                      f"{', '.join(f'{v:.3f}' for v in e['masses'][:10])}{' …' if e['K'] > 10 else ''} | {', '.join(f'{v:+.3f}' for v in e['c'][:10])} | "
                      f"{', '.join(str(int(v)) for v in e['n'][:10])} | {np.abs(e['a']).max() / 4:.3f} |")
        L_.append("")
    if C is not None:
        L_ += ["## C. The finest resolution stage (the fourth bump, compact R = 1, σ₀/64, from the jittered lift of ϑ^b(q), 100000 slow units)", "",
               "| γ | seed | clusters at the end (> 1% of the mass) | fit error | masses | centres | particles |", "|---|---|---|---|---|---|---|"]
        for r in sorted(C, key=lambda r: (r['gamma'], r['seed'])):
            e = r['end']
            L_.append(f"| {r['gamma']:g} | {r['seed']} | {e['K']} ({e['K1']}) | {fmt(e['fit'])} | {', '.join(f'{v:.3f}' for v in e['masses'])} | "
                      f"{', '.join(f'{v:+.3f}' for v in e['c'])} | {', '.join(str(int(v)) for v in e['n'])} |")
        L_.append("")
    if F is not None:
        L_ += ["## F. The targets where Proposition D‴.32's count failed (Part K of `exp_cohesion_leak.py`; the annealed schedule, broad starts, m = 60, γ = 0.01)", "",
               "| target | Σ|A_k|/a_max | seed | clusters at the end (> 1% of the mass) | fit error | masses | centres | needed |A_k|/a_max | |a|/a_max (max) |",
               "|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(F, key=lambda r: (r['name'], r['seed'])):
            e = r['end']
            need = [abs(b[0]) / 4.0 for b in E.TARGETS_K[r['name']]]
            L_.append(f"| {r['name']} | {sum(need):.3f} | {r['seed']} | {e['K']} ({e['K1']}) | {fmt(e['fit'])} | {', '.join(f'{v:.3f}' for v in e['masses'][:8])}{' …' if e['K'] > 8 else ''} | "
                      f"{', '.join(f'{v:+.3f}' for v in e['c'][:8])}{' …' if e['K'] > 8 else ''} | {', '.join(f'{v:.3f}' for v in need)} | {np.abs(e['a']).max() / 4:.3f} |")
        L_.append("")
    if Ep is not None:
        r = Ep['real']
        L_ += ["## E. The real loop and the linearisation", "",
               f"The real loop (k_e = {KE_E:g}, γ = {GAMMA_E:g}, the positions at the instance's atoms, β(0) = (1.6, 0.7, 0.7), {T_E:g} time units): V = ½e² + "
               f"(1/(γm))Σβ*log(β*/β) from {fmt(r['V0'])} to {fmt(r['Vend'])}, {'monotone' if r['monotone'] else 'NOT monotone'}; V̇ (from the vector field) + k_ee² "
               f"within {fmt(r['dVerr'])} of max k_ee²; Σβ − m ≤ {fmt(r['sum_err'])}; β at the end "
               + ", ".join(f'{v:.5f}' for v in r['beta_end']) + ".", ""]
        for g, o in Ep['lin'].items():
            ev = o['ev']
            L_.append(f"- The averaged loop's linearisation in (Θ, log β) at the lift of ϑ^e(w), counts (6, 5, 4), γ = {g:g}: the largest drift {fmt(o['drift'])}; "
                      f"{o['zero']} eigenvalues below 10⁻⁹ in modulus; the others from {ev[0]:.4e} to {ev[np.abs(ev) >= 1e-9].max():.4e}.")
        g0 = Ep['lin'][0.0]['ev']
        nz0 = np.sort(g0[np.abs(g0) >= 1e-9]); pr = np.sort(Ep['pred'])
        err = fmt(float(np.abs(nz0 - pr).max())) if len(nz0) == len(pr) else '— (counts differ)'
        L_.append(f"- At γ = 0 the {len(nz0)} nonzero eigenvalues are Proposition D‴.32's predicted spectrum (the separations at −σ₀w_k, the collective modes) to {err}.")
        L_.append("")
    L_ += ["Numbers are produced by `experiments/exp_mass_transfer.py`; nothing is edited by hand."]
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mass_transfer_raw.npz')
    md = os.path.join(OUT, 'mass_transfer.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(old.get('A'), old.get('B'), old.get('C'), old.get('D'), old.get('E'), old.get('F')))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEF')
    from multiprocessing import Pool
    res = dict(old)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'A' in parts:
            h['A'] = pool.map_async(a_job, [(g, s) for s in (1.0, 1.25) for g in GAMMAS])
        if 'B' in parts:
            h['B'] = pool.map_async(broad_job, [('B', g, s) for g in GAMMAS for s in SEEDS])
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(g, s) for g in (0.0, 1e-2) for s in SEEDS])
        if 'D' in parts:
            h['D'] = pool.map_async(broad_job, [('D', 1e-2, s) for s in SEEDS])
        if 'E' in parts:
            h['E'] = pool.apply_async(part_E)
        if 'F' in parts:
            h['F'] = pool.map_async(broad_job, [('F', 1e-2, s, name) for name in ('three heavy', 'two close') for s in SEEDS])
        for k in 'ABCDEF':
            if k in h:
                res[k] = h[k].get()
                E.save(raw, **res)
    open(md, 'w').write(report(res.get('A'), res.get('B'), res.get('C'), res.get('D'), res.get('E'), res.get('F')))


if __name__ == '__main__':
    main()
