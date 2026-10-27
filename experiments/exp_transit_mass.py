"""
Proposition D‴.51: the transit at a finite mass — near ties share a feature, the entrant's overshoot, the creep, the predicted mass
====================================================================================================================================

The free law with moving centres (Proposition D‴.48(iii); m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025; Gaussian kernels of width 0.12
on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, masses ε) at a finite mass, against the ε-free inner problem of a
near tie and the predictions made from the sliding path alone (`rich_transit`).

A. The near-tie law: one bump between the kernels at 0.0 and 0.1, δ from the midpoint (amplitudes 0.3, 0.8, 2.0; δ = 10⁻⁴, 10⁻³, 10⁻²;
   ε = 10⁻⁵ … 10⁻¹⁰): the runner-up's share of the bump against x = μ*/ε^{1/2} (μ* its margin when the entrant enters), against the
   inner problem's share at the same x (`near_tie_inner`, at the exact tie: the limit ε → 0 at fixed x), the Gibbs split 1/(1 + e^{κx}) (κ = πρ(2c)^{1/2}/a) and the
   frozen coefficient K'³ρ⁻³/(|A|(x − Δ)³) (Δ the inner problem's pull, B).
B. The entrant's overshoot η_c(α_A)/(Λ₀ρ) (`overshoot`) against the pull Δ(A) the inner problem exerts on a frozen runner-up (x = 200) and
   the loop's (the runner-up's final inner coordinate read from its coefficient, δ = 10⁻³, ε = 10⁻¹⁰); I∞ = overshoot − Δ.
C. The creep: target 10 of Proposition D‴.48's random targets (C4, C6), whose path leaves a kernel 0.0065 inside its wall beside an
   unfitted bump — its predicted creep threshold (`recovery_threshold`) bracketed by the loop.
D. Prediction from the path, out of sample: 48 random three-bump targets (D‴.48's recipe, a fresh seed); the predicted ε_c against the loop
   at ε = 10⁻⁵, 10⁻⁶, 10⁻⁷, 10⁻⁸ (formed: |w| > 10⁻²max|A*|; recovered: exactly the path's kernels formed).
E. The predicted ε_c over 400 random targets: the fraction predicted not to recover at ε, entry- and creep-limited, and its slope in log ε.

Output: results/chaos/transit_mass.md and transit_mass_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D,E`).
"""
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('TM_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TM_WORKERS', '3'))
AMPS_A, OFFS_A = (0.3, 0.8, 2.0), (1e-4, 1e-3, 1e-2)
EPS_A = (1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
TAU_A = 1.32                                                          # the run to 1.32 τ_entry (the relaxation after it)
AMPS_B = (0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 4.0)
X_B, OFF_B, EPS_B = 200.0, 1e-3, 1e-10
CREEP_TARGET, EPS_C = 10, (1e-9, 2e-9, 3.5e-9, 7e-9, 1.4e-8, 3e-8)
N_D, SEED_D, EPS_D = 48, 51, (1e-5, 1e-6, 1e-7, 1e-8)
N_E, SEED_E = 400, 52
EPS_GRID_E = tuple(10.0 ** -k for k in np.arange(4, 12.5, 0.5))


def setup():
    warnings.filterwarnings('ignore')
    L = X.law()
    K = RT.Bumps(np.linspace(-1, 1, X.NS_C), X.SIG_C)
    c0 = np.linspace(-1, 1, X.D_C)
    return L, K, c0


def random_targets(seed, n):
    """D‴.48's recipe: three bumps of the kernels' width, centres in [−0.8, 0.8] at least 0.3 apart, amplitudes ±U(0.4, 1)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        while True:
            cs = np.sort(rng.uniform(-0.8, 0.8, 3))
            if np.diff(cs).min() >= 0.3:
                break
        As = rng.choice([-1, 1], 3) * rng.uniform(0.4, 1.0, 3)
        out.append((cs, As))
    return out


def near_tie_case(A, off):
    L, K, c0 = setup()
    cs = 0.05 + off
    f = X.bumps_target(K, (cs,), (A,))
    g = K.phi(c0).T @ f / K.xs.size
    order = np.argsort(-np.abs(g))
    w_, n_ = int(order[0]), int(order[1])
    return L, K, c0, cs, f, g, w_, n_, 1 - abs(g[n_]) / abs(g[w_])


def a_run(args):
    A, off, eps = args
    L, K, c0, cs, f, g, w_, n_, mu = near_tie_case(A, off)
    t0 = time.time()
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, TAU_A / abs(g[w_]), 'paradigm', n_t=331)
    x = mu / np.sqrt(eps)
    inner = RT.near_tie_inner(L, K, A, 0.5 * (c0[w_] + c0[n_]), c0[w_], c0[n_], x)     # the limit ε → 0 at fixed x: the exact tie
    wn, ww = abs(W[-1][n_]), abs(W[-1][w_])
    return dict(A=A, off=off, eps=eps, mu=float(mu), x=float(x), entrant=w_, runner=n_, share=float(wn / (ww + wn)), w_r=float(wn),
                c_r=float(cc[-1][n_]), c_e=float(cc[-1][w_]), share_inner=float(inner['share']), secs=time.time() - t0)


def part_A():
    from multiprocessing import Pool
    jobs = [(A, off, eps) for A in AMPS_A for off in OFFS_A for eps in EPS_A]
    with Pool(N_WORKERS) as pool:
        rows = pool.map(a_run, jobs)
    return rows


def b_run(A):
    L, K, c0 = setup()
    ov = float(RT.overshoot(L, A))
    fr = RT.near_tie_inner(L, K, A, 0.05, 0.1, 0.0, X_B, both=False)
    L, K, c0, cs, f, g, w_, n_, mu = near_tie_case(A, OFF_B)
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, EPS_B, TAU_A / abs(g[w_]), 'paradigm', n_t=331)
    x = mu / np.sqrt(EPS_B)
    xi_r = float(RT.inner_xi(L, abs(W[-1][n_])))
    fr_off = RT.near_tie_inner(L, K, A, cs, c0[w_], c0[n_], X_B, both=False)
    return dict(A=A, overshoot=ov, delta_inner=float(X_B - fr['xi_r']), xi_e=fr['xi_e'], delta_inner_off=float(X_B - fr_off['xi_r']),
                x_loop=float(x), delta_loop=float(x - xi_r))


def part_B():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_run, list(AMPS_B))


def c_run(args):
    eps, cs, As = args
    L, K, c0 = setup()
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C))
    te = 1.3 * R['events'][-1][0]
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, te, 'paradigm')
    sm = X.summarize(K, f, c0, taus, W, cc, X.TH_ACT * np.abs(As).max())
    path = X.net_entries(R)
    return dict(eps=eps, n_act=sm['n_act'], order=sm['order'], ok=bool(sm['n_act'] == len(path) and set(sm['order']) == set(path)),
                disp=float(np.abs(cc[-1] - c0).max()))


def part_C():
    from multiprocessing import Pool
    L, K, c0 = setup()
    cs, As = random_targets(21, X.N_RAND)[CREEP_TARGET]
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C))
    pred = RT.recovery_threshold(L, K, f, c0, R, X.TH_ACT * np.abs(As).max())
    with Pool(N_WORKERS) as pool:
        rows = pool.map(c_run, [(e, cs, As) for e in EPS_C])
    return dict(cs=cs.tolist(), As=As.tolist(), path=X.net_entries(R), pred=pred, rows=rows)


def d_pred(args):
    n, cs, As = args
    L, K, c0 = setup()
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(X.D_C))
    pred = RT.recovery_threshold(L, K, f, c0, R, X.TH_ACT * np.abs(As).max())
    return dict(n=n, cs=cs.tolist(), As=As.tolist(), path=X.net_entries(R), slide_n=int(R['act'].sum()),
                slide_rms=float(R['hist'][-1][3]), te=float(1.3 * R['events'][-1][0]), pred=pred)


def d_loop(args):
    n, cs, As, eps, te, path = args
    L, K, c0 = setup()
    f = X.bumps_target(K, cs, As)
    taus, W, cc, Lam = X.run_moving(L, K, f, c0, eps, te, 'paradigm')
    sm = X.summarize(K, f, c0, taus, W, cc, X.TH_ACT * np.abs(As).max())
    return dict(n=n, eps=eps, n_act=sm['n_act'], order=sm['order'], ok=bool(sm['n_act'] == len(path) and set(sm['order']) == set(path)))


def part_D():
    from multiprocessing import Pool
    tg = random_targets(SEED_D, N_D)
    with Pool(N_WORKERS) as pool:
        preds = pool.map(d_pred, [(n, cs, As) for n, (cs, As) in enumerate(tg)])
        jobs = [(p['n'], np.array(p['cs']), np.array(p['As']), eps, p['te'], p['path']) for p in preds for eps in EPS_D]
        loops = pool.map(d_loop, jobs)
    return dict(preds=preds, loops=loops)


def part_E():
    from multiprocessing import Pool
    tg = random_targets(SEED_E, N_E)
    with Pool(N_WORKERS) as pool:
        preds = pool.map(d_pred, [(n, cs, As) for n, (cs, As) in enumerate(tg)])
    return dict(preds=[dict(n=p['n'], slide_n=p['slide_n'], eps_entry=p['pred']['eps_entry'], eps_creep=p['pred']['eps_creep'])
                       for p in preds])


def report(A=None, B=None, C=None, D=None, Ee=None):
    fmt = E.fmt
    L, K, c0 = setup()
    kap = np.pi * np.sqrt(2 * L.c) / L.A
    Kp3 = (8 / np.pi ** 3) * L.A / (L.m * np.sqrt(2 * L.c))
    L_ = ['# Proposition D‴.51: the transit at a finite mass (`exp_transit_mass.py`)', '',
          "The free law with moving centres (Proposition D‴.48(iii); m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025), Gaussian kernels of width "
          "0.12 on 61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, masses ε; the inner problem, the overshoot and the "
          "predictions from `rich_transit`.", '']
    dl = {r['A']: r['delta_inner'] for r in B} if B is not None else {}
    if A is not None:
        L_ += ['## A. The near-tie law', '',
               "One bump between the kernels at 0.0 and 0.1, δ from their midpoint; μ* the runner-up's margin when the entrant enters, "
               "x = μ*/ε^{1/2}; the runner-up's share |w_r|/(|w_e| + |w_r|) of the bump in the loop (at 1.32 τ_entry) against the ε-free inner "
               "problem's at the same x (its limit ε → 0 at fixed x: the bump at the midpoint, the runner-up x behind; the loop's ratio to it → 1 "
               f"as δ ∝ ε^{{1/2}}x → 0), the Gibbs split 1/(1 + e^{{κx}}) (κ = πρ(2c)^{{1/2}}/a = {kap:.4f}) and the frozen "
               f"coefficient K'³ρ⁻³/(|A|(x − Δ)³) (K'³ = {Kp3:.4f}; Δ from B):", '',
               "| A | δ | ε | x | share: loop | inner problem | loop/inner | Gibbs split | frozen | the runner-up's centre |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(A, key=lambda r: (r['A'], r['x'])):
            gib = 1 / (1 + np.exp(kap * r['x']))
            dd = dl.get(r['A'])
            fro = Kp3 / (r['A'] * (r['x'] - dd) ** 3) if dd is not None and r['x'] > dd + 1 else float('nan')
            L_.append(f"| {r['A']:g} | {r['off']:g} | {r['eps']:g} | {r['x']:.4g} | {fmt(r['share'])} | {fmt(r['share_inner'])} | "
                      f"{r['share'] / r['share_inner']:.4f} | {fmt(gib)} | {fmt(fro)} | {r['c_r']:.4f} |")
        rat = np.array([r['share'] / r['share_inner'] for r in A])
        L_ += ['', f"The loop's share over the inner problem's: {rat.min():.4f}–{rat.max():.4f} over all {len(A)} runs "
               f"(x from {min(r['x'] for r in A):.3g} to {max(r['x'] for r in A):.4g}).", '']
    if B is not None:
        L_ += ["## B. The entrant's overshoot and the runner-up's pull", '',
               "The overshoot η_c(α_A)/(Λ₀ρ) (an entered kernel's dual beyond its wall, in units of ε^{1/2}ρ) against the pull Δ the inner "
               f"problem exerts on a frozen runner-up (x = {X_B:g}; at the exact tie and at δ = {OFF_B:g}) and the loop's (δ = {OFF_B:g}, "
               f"ε = {EPS_B:g}: Δ = x − ξ_r(∞), ξ_r read from the runner-up's coefficient); I∞ = overshoot − Δ:", '',
               "| A | overshoot | Δ: inner, tie | I∞ | Δ: inner, δ | Δ: loop (x) |", "|---|---|---|---|---|---|"]
        L_ += [f"| {r['A']:g} | {r['overshoot']:.4f} | {r['delta_inner']:.4f} | {r['overshoot'] - r['delta_inner']:.4f} | "
               f"{r['delta_inner_off']:.4f} | {r['delta_loop']:.4f} ({r['x_loop']:.4g}) |" for r in B]
        L_.append('')
    if C is not None:
        p = C['pred']
        L_ += ['## C. The creep', '',
               f"Target {CREEP_TARGET} of D‴.48's random targets (centres {[round(v, 3) for v in C['cs']]}, amplitudes "
               f"{[round(v, 3) for v in C['As']]}; the path's kernels {C['path']}). Predicted from the path: the entry threshold "
               f"{fmt(p['eps_entry'])} (kernel {p['entry'][1] if p['entry'] else '—'} at the entry of {p['entry'][0] if p['entry'] else '—'}, "
               f"margin {fmt(p['entry'][2]) if p['entry'] else '—'}), the creep thresholds "
               f"{', '.join(f'{j}: {fmt(v)}' for j, v in p['creep'].items())}: ε_c = {fmt(p['eps_c'])}. The loop:", '',
               "| ε | ε/ε_c | formed | entrants | the path's kernels | the inactive kernels' largest displacement |", "|---|---|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | {r['eps'] / p['eps_c']:.3g} | {r['n_act']} | {r['order']} | {'yes' if r['ok'] else 'no'} | {fmt(r['disp'])} |"
               for r in C['rows']]
        L_.append('')
    if D is not None:
        preds = {p['n']: p for p in D['preds']}
        L_ += ['## D. Prediction from the path (out of sample)', '',
               f"{len(preds)} random three-bump targets (seed {SEED_D}); ε_c = min(ε_entry, ε_creep) predicted from the sliding path; the loop "
               "recovers when it forms exactly the path's kernels. Predicted (ε < ε_c) against observed:", '',
               "| ε | recovered | predicted to recover | agree | recovered though not predicted | not recovered though predicted |",
               "|---|---|---|---|---|---|"]
        miss = []
        for eps in EPS_D:
            rr = [r for r in D['loops'] if r['eps'] == eps]
            pr = [eps < preds[r['n']]['pred']['eps_c'] for r in rr]
            ob = [r['ok'] for r in rr]
            L_.append(f"| {eps:g} | {sum(ob)} of {len(ob)} | {sum(pr)} | {sum(p_ == o for p_, o in zip(pr, ob))} | "
                      f"{sum((not p_) and o for p_, o in zip(pr, ob))} | {sum(p_ and (not o) for p_, o in zip(pr, ob))} |")
            miss += [(r['n'], eps) for r, p_, o in zip(rr, pr, ob) if p_ != o]
        dis = ', '.join(f"({n}, {e:g}, {preds[n]['pred']['eps_c']:.3g})" for n, e in miss) or 'none'
        L_ += ['', f"The disagreements (target, ε, ε_c): {dis}; the "
               f"path itself exact (three kernels, rms residual < 10⁻⁸) on {sum(1 for p in preds.values() if p['slide_n'] == 3 and p['slide_rms'] < 1e-8)} "
               f"of {len(preds)}.", '']
    if Ee is not None:
        pe = Ee['preds']
        ee = np.array([p['eps_entry'] for p in pe]); ec = np.array([p['eps_creep'] for p in pe]); eps_c = np.minimum(ee, ec)
        L_ += ['## E. The predicted mass over many targets', '',
               f"{len(pe)} random three-bump targets (seed {SEED_E}): the fraction predicted not to recover at ε (ε ≥ ε_c), and the fractions "
               "limited by an entry's near tie (ε ≥ ε_entry) and by the creep (ε ≥ ε_creep); the local slope of each in log ε:", '',
               "| ε | not recovered | by a near tie | by the creep | slope (all) |", "|---|---|---|---|---|"]
        fr = [np.mean(eps_c <= e) for e in EPS_GRID_E]
        for k, e in enumerate(EPS_GRID_E):
            sl = (np.log(fr[k - 1] / fr[k]) / np.log(EPS_GRID_E[k - 1] / e)) if k > 0 and fr[k] > 0 and fr[k - 1] > 0 else float('nan')
            L_.append(f"| {e:.3g} | {fr[k]:.4f} | {np.mean(ee <= e):.4f} | {np.mean(ec <= e):.4f} | {fmt(sl)} |")
        L_.append('')
    L_.append('Numbers are produced by `experiments/exp_transit_mass.py`; nothing is edited by hand.')
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'transit_mass_raw.npz')
    md = os.path.join(OUT, 'transit_mass.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDE']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDE')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('B', part_B), ('A', part_A), ('C', part_C), ('D', part_D), ('E', part_E)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE']))


if __name__ == '__main__':
    main()
