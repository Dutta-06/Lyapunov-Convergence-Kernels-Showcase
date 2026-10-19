"""
Theorem D‴.21: a measure-level basin for the barycentric OT leak — the Polyak–Łojasiewicz inequality along the runs
================================================================================================================

The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit, as in `exp_ot_leak.py`: commensurate sizes m = 15, 60, 240 and
M = 510 from the warm clouds of Proposition D‴.12 (the jittered mixture of width 0.02 around the exact representation, one cloud per seed,
the m-system its first m rows), four seeds, 600 slow units, midpoint rule with dt = 0.02 and the leak target recomputed every 10 steps, under
the OT leak and — for the contrast — the nearest-atom leak. Two priors: the matching prior (the exact representation; r* = 0, μ* = ν) and the
prior with centres shifted by +0.15 (μ* = Θ*(q) by Newton, r* ≠ 0).

Records every 5 slow units, for every system, against μ* (src/core/ot_basin.py): the excess energy E = 𝒥(μ_Θ) − 𝒥* (𝒥 := L̄ + (σ₀/2)W²(·, ν),
the OT functional, evaluated along both leaks' trajectories), the class-respecting distance D_a, W²(μ_Θ, μ*), ‖Δf̂‖², the exact effective
margins κ_E and κ_Q of the identities (E) and (P) and their residuals, the observed PL ratio |∇J_a|²_G/(mE) (= −d log E/ds along the OT flow)
and its lower bound Q²/(DE), and the plan changes; every 25 slow units the plan-segment curvatures κ*, κ*_−, κ^Θ (9 points per segment),
the margin α = σ₀ − (κ* + κ^Θ)/2 and the theorem's rate λ = α·min{2, 2α/(σ₀ + κ*_−)}.

The report adds: the a-posteriori bound E(0)exp(−∫λ₊) against E(s); the a-priori basin E₀ = κ²/(2(L₂^G)²) with κ := σ₀ − sup λ_max R^G[μ*]
over a box containing every trajectory (by search); for the nearest-atom runs, the end energy against the energy of the replicated
configuration of the sampled weights (the spurious minima).

Output: results/chaos/ot_basin.md and ot_basin_raw.npz (OTB_SEEDS, OTB_S_END, OTB_WORKERS to change the pool; `--from-raw` regenerates the md).
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
from src.core import prior_design as PD
from src.core import ot_leak as OT
from src.core import ot_basin as OB

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
M_REF, MS = 510, (15, 60, 240)
S_END = float(os.environ.get('OTB_S_END', 600.0))
DT, PLAN_EVERY, REC, SEG = 0.02, 10, 250, 1250                          # records every 5 slow units, segment curvatures every 25
JITTER = 0.02
N_SEEDS = int(os.environ.get('OTB_SEEDS', 4))
N_WORKERS = int(os.environ.get('OTB_WORKERS', 2))
Q = np.ones(3) / 3
SHIFTS = {'matching': 0.0, 'shifted': 0.15}
TIMES = (0.0, 25.0, 100.0, 300.0, 600.0)


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def representation(net, shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def setup(prior_name):
    net = make_net()
    prior = representation(net, SHIFTS[prior_name])
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, prior, n_quad=N_Q)
    atoms, res, marg = loop.stationary_atoms(prior)
    return net, loop, atoms, res, marg


def record(B, Th, prev_lab):
    t = B.terms(Th)
    W2s = OT.W2_empirical(Th, B.atoms, G)
    sw = int(np.sum(t['lab'] != prev_lab)) if prev_lab is not None else 0
    rel = lambda a, b: abs(a - b) / max(abs(a), 1e-300)
    return dict(E=t['E'], D=t['D'], W2s=W2s, ndf2=t['ndf2'], kE=t['kappa_E'], kQ=t['kappa_Q'], ratio=t['ratio'], lam_eff=t['lam_eff'],
                idE=rel(t['E'], t['E_ident']), idQ=rel(t['Q'], t['Q_ident']), switches=sw), t['lab'], t['r']


def job(args):
    seed, prior_name, leak = args
    t0 = time.time()
    net, loop, atoms, res, marg = setup(prior_name)
    B = OB.BasinDiagnostics(loop, atoms)
    exact = representation(net, 0.0)
    rng = np.random.default_rng(seed)
    clouds = CL.mixture_start(exact, M_REF, rng, JITTER, G)
    sys_ = {m: clouds[:m].copy() for m in list(MS) + [M_REF]}
    n_steps = int(round(S_END / DT))
    keys = ('E', 'D', 'W2s', 'ndf2', 'kE', 'kQ', 'ratio', 'lam_eff', 'idE', 'idQ', 'switches')
    rec = {m: dict(s=[], **{k: [] for k in keys}) for m in sys_}
    seg = {m: dict(s=[], kstar=[], kstar_minus=[], ktheta=[], alpha=[], lam=[]) for m in sys_}
    prev = {m: None for m in sys_}
    T = {}
    for k in range(n_steps + 1):
        if k % REC == 0:
            for m, Th in sys_.items():
                d, lab, r = record(B, Th, prev[m]); prev[m] = lab
                rec[m]['s'].append(k * DT)
                for kk in keys:
                    rec[m][kk].append(d[kk])
                if k % SEG == 0:
                    sc = B.segment_curvatures(Th, lab, n_t=9, r=r)
                    al, lam = B.theorem_rate(**sc)
                    seg[m]['s'].append(k * DT); seg[m]['kstar'].append(sc['kappa_star']); seg[m]['kstar_minus'].append(sc['kappa_star_minus'])
                    seg[m]['ktheta'].append(sc['kappa_theta']); seg[m]['alpha'].append(al); seg[m]['lam'].append(lam)
        if k == n_steps:
            break
        if k % PLAN_EVERY == 0:
            for m, Th in sys_.items():
                T[m] = loop.leak_target(Th) if leak == 'ot' else loop.nearest.leak_target(Th)
        for m, Th in list(sys_.items()):
            k1 = loop.drift(Th, T[m])
            sys_[m] = Th + DT * loop.drift(Th + 0.5 * DT * k1, T[m])
    out = dict(seed=seed, prior=prior_name, leak=leak, atoms=atoms, newton_res=res, diag_margin=marg, Jstar=B.Jstar, stationarity=B.stationarity,
               rec={m: {k: np.array(v) for k, v in r.items()} for m, r in rec.items()},
               seg={m: {k: np.array(v) for k, v in r.items()} for m, r in seg.items()},
               final={m: Th.astype(np.float32) for m, Th in sys_.items()})
    # the nearest-atom contrast: the replicated configuration of the sampled weights (a critical point of the nearest-atom energy)
    if leak == 'nearest':
        spur = {}
        for m, Th in sys_.items():
            lab = np.argmin(loop.cost_matrix(Th), axis=1)
            counts = np.bincount(lab, minlength=3)
            if np.all(counts > 0):
                w = counts / m
                Tw, resw = PD.weighted_stationary_atoms_prior(loop.nearest, w, atoms, w0=Q, n_cont=4)
                lift = np.repeat(Tw, counts, axis=0)
                d0 = float(np.abs(loop.nearest.drift(lift)).max())
                dOT = float(np.abs(loop.drift(lift)).max())
                spur[m] = dict(counts=counts, newton_res=resw, E=B.terms(lift)['E'], drift_nearest=d0, drift_ot=dOT)
        out['spurious'] = spur
    print(f"[{time.time() - t0:.0f}s] seed {seed} {prior_name} {leak}", flush=True)
    return out


def a_priori_constants(runs):
    """L₂^G (search, as in exp_chaos_learning.py) and, for each prior, κ := σ₀ − sup λ_max R^G[μ*] over the bounding box of every
    configuration recorded at the end of the runs and of the prior and stationary atoms, enlarged by 10% (random search, 20000 points)."""
    net = make_net()
    Lc = CL.whitened_lipschitz_constants(net, G, 2.5, n_samples=60000)
    out = dict(L2G=Lc['L_2_G'])
    for name in SHIFTS:
        net, loop, atoms, res, marg = setup(name)
        B = OB.BasinDiagnostics(loop, atoms)
        pts = [atoms, loop.prior] + [np.asarray(r['final'][m], float) for r in runs if r['prior'] == name for m in r['final']]
        P = np.vstack(pts)
        lo, hi = P.min(axis=0), P.max(axis=0); pad = 0.1 * (hi - lo)
        rng = np.random.default_rng(11)
        X = rng.uniform(lo - pad, hi + pad, (20000, 3))
        R = loop.residual_curvature(B.r_s, X)
        Rw = SQG[None, :, None] * R * SQG[None, None, :]
        sup = float(np.linalg.eigvalsh(Rw)[:, -1].max())
        kappa = SIGMA0 - sup
        out[name] = dict(sup_kappa_R=sup, kappa=kappa, E0=(kappa ** 2 / (2 * Lc['L_2_G'] ** 2) if kappa > 0 else float('nan')), box=(lo - pad, hi + pad))
    return out


def moved_fraction(ms=(30, 300), ks=(1, 2, 3, 5), t_move=0.6, s_run=300.0, ds=0.5):
    """The measure-level basin at the matching prior: the lift of q with k particles of class 0 (k = 1, 2, 3, 5 at m = 30, ten times as many at
    m = 300: the same fractions) moved t_move of the way to atom 1, into cell 1. At s = 0: the plan's labels of the moved particles, κ^Θ, the
    margin α and λ. For the smallest fraction a run of s_run slow units under both leaks (midpoint rule, plan every step): E at the end, the
    observed ratio at the end, the a-posteriori bound with λ held at the smaller endpoint over every 5 slow units."""
    net, loop, atoms, res, marg = setup('matching')
    B = OB.BasinDiagnostics(loop, atoms)
    rows, runs = [], {}
    for m in ms:
        lab0 = np.repeat(np.arange(3), m // 3)
        for k0 in ks:
            k = k0 * m // ms[0]
            Th = atoms[lab0].copy()
            Th[:k] = (1.0 - t_move) * atoms[0] + t_move * atoms[1]
            t = B.terms(Th)
            sc = B.segment_curvatures(Th, t['lab'], n_t=9)
            al, lam = B.theorem_rate(**sc)
            rows.append(dict(m=m, k=k, p=k / m, kept=bool(np.all(t['lab'][:k] == 0)), E=t['E'], D=t['D'], ktheta=sc['kappa_theta'], alpha=al, lam=lam, ratio=t['ratio']))
        for leak in ('ot', 'nearest'):
            Th = atoms[lab0].copy(); Th[0] = (1.0 - t_move) * atoms[0] + t_move * atoms[1]
            t = B.terms(Th); al, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab'], n_t=5))
            E0, bound, lam_prev = t['E'], 0.0, lam
            n = int(round(s_run / ds)); every = int(round(5.0 / ds))
            for j in range(n):
                T = loop.leak_target(Th) if leak == 'ot' else loop.nearest.leak_target(Th)
                k1 = loop.drift(Th, T); Th = Th + ds * loop.drift(Th + 0.5 * ds * k1, T)
                if j % every == every - 1:
                    t = B.terms(Th); al, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab'], n_t=5))
                    bound += every * ds * min(lam_prev, lam); lam_prev = lam
            runs[(m, leak)] = dict(E0=E0, E=t['E'], ratio=t['ratio'], bound=E0 * np.exp(-bound), lab0=int(t['lab'][0]), cell=int(np.argmin(loop.cost_matrix(Th[:1])[0])))
    return rows, runs


# ---------------------------------------------------------------------------- report
def at(s, v, t):
    return float(v[int(np.argmin(np.abs(s - t)))])


def fmt_ms(vals, p=2):
    vals = np.asarray(vals, float)
    return f"{vals.mean():.{p}e}" + (f" ({vals.min():.1e}–{vals.max():.1e})" if len(vals) > 1 else "")


def report(runs, consts, t0=None):
    lines = ["# Theorem D‴.21: a measure-level basin for the OT leak — the Polyak–Łojasiewicz inequality along the runs (`exp_ot_basin.py`)", "",
             f"The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f} (2σ₀ = {2 * SIGMA0:.5f}); sizes m = {', '.join(map(str, MS))} and M = {M_REF} from the warm clouds "
             f"of Proposition D‴.12 (width {JITTER:g} around the exact representation; the m-system the first m rows); {N_SEEDS} seeds; {S_END:g} slow units, midpoint rule dt = {DT:g}, leak target every "
             f"{PLAN_EVERY} steps. E = 𝒥(μ_Θ) − 𝒥* with 𝒥 = L̄ + (σ₀/2)W²(·, ν) the OT functional, evaluated along both leaks' trajectories; D_a the class-respecting distance to μ*; the observed PL ratio "
             "|∇J_a|²_G/(mE), which along the OT flow is −d log E/ds; λ the theorem's rate from the plan-segment curvatures (9 points per segment, every 25 slow units). "
             "Entries: mean over the seeds (min–max).", ""]
    for name in SHIFTS:
        rr = [r for r in runs if r['prior'] == name]
        r0 = rr[0]
        lin = None
        lines += [f"## Prior: {name} (centres shifted by {SHIFTS[name]:g})", "",
                  f"μ*: Newton residual {r0['newton_res']:.1e}, stationarity of the atoms {r0['stationarity']:.1e}, 𝒥* = {r0['Jstar']:.4e}; diagonal plan: ε = {r0['diag_margin'][0]:.3f} < Δ_ν/2 = {r0['diag_margin'][1]:.3f}. "
                  f"Identities (E), (P) along every run: largest relative residual {max(float(np.max(r['rec'][m]['idE'])) for r in rr for m in r['rec']):.1e} (E) and {max(float(np.max(r['rec'][m]['idQ'])) for r in rr for m in r['rec']):.1e} (P).", "",
                  "| leak | m | E at s = " + ", ".join(f"{t:g}" for t in TIMES) + " | D_a at s = 0, 100, 600 | W²(μ_Θ, μ*) at s = 0, 100, 600 | plan changes (particle relabellings) |",
                  "|---|---|---|---|---|---|"]
        for leak in ('ot', 'nearest'):
            for m in list(MS) + [M_REF]:
                rs = [r for r in rr if r['leak'] == leak]
                Es = [" / ".join(f"{at(r['rec'][m]['s'], r['rec'][m]['E'], t):.2e}" for r in rs) for t in TIMES]
                Em = [np.mean([at(r['rec'][m]['s'], r['rec'][m]['E'], t) for r in rs]) for t in TIMES]
                Dm = [np.mean([at(r['rec'][m]['s'], r['rec'][m]['D'], t) for r in rs]) for t in (0.0, 100.0, 600.0)]
                Wm = [np.mean([at(r['rec'][m]['s'], r['rec'][m]['W2s'], t) for r in rs]) for t in (0.0, 100.0, 600.0)]
                sw = [int(np.sum(r['rec'][m]['switches'])) for r in rs]
                lines.append(f"| {leak} | {m} | " + ", ".join(f"{x:.2e}" for x in Em) + " | " + ", ".join(f"{x:.2e}" for x in Dm) + " | " + ", ".join(f"{x:.2e}" for x in Wm)
                             + f" | {min(sw)}–{max(sw)} |")
        # the PL inequality along the OT runs
        lines += ["", "The PL inequality along the OT runs. s₊: the first segment record after which the margin α stays positive; the observed ratio and the bounds over the records "
                  "with α > 0 (the chain ratio ≥ Q²/(DE) ≥ λ must hold); the late values: means over s ∈ [500, 600] (records with E > 10⁻²⁶).", "",
                  "| m | s₊ (seeds) | E(s₊)/E(0) | min over α > 0 of ratio/λ and Q²/(DE)/λ | late observed ratio | late λ | late α | late κ*, κ^Θ | a-posteriori bound E(0)e^{−∫λ₊} at s = 600 against E(600) |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for m in list(MS) + [M_REF]:
            rs = [r for r in rr if r['leak'] == 'ot']
            splus, Erat, mins, late_ratio, late_lam, late_al, late_k, bounds = [], [], [], [], [], [], [], []
            for r in rs:
                sg, rc = r['seg'][m], r['rec'][m]
                pos = sg['alpha'] > 0
                idx = None
                for j in range(len(pos)):
                    if np.all(pos[j:]):
                        idx = j; break
                sp = float(sg['s'][idx]) if idx is not None else float('nan')
                splus.append(sp)
                Erat.append(at(rc['s'], rc['E'], sp) / rc['E'][0] if idx is not None else float('nan'))
                # ratio and lam_eff at the segment records
                ok = []
                for j, s_ in enumerate(sg['s']):
                    if sg['alpha'][j] > 0 and sg['lam'][j] > 0:
                        i = int(np.argmin(np.abs(rc['s'] - s_)))
                        if rc['E'][i] > 1e-26:
                            ok.append((rc['ratio'][i] / sg['lam'][j], rc['lam_eff'][i] / sg['lam'][j]))
                mins.append((min(x[0] for x in ok), min(x[1] for x in ok)) if ok else (float('nan'), float('nan')))
                late = (rc['s'] >= 500.0) & (rc['E'] > 1e-26)
                late_ratio.append(float(np.mean(rc['ratio'][late])) if late.any() else float('nan'))
                lateg = sg['s'] >= 500.0
                late_lam.append(float(np.mean(sg['lam'][lateg]))); late_al.append(float(np.mean(sg['alpha'][lateg])))
                late_k.append((float(np.mean(sg['kstar'][lateg])), float(np.mean(sg['ktheta'][lateg]))))
                # a-posteriori bound with λ piecewise constant between the segment records
                lam_t = np.interp(rc['s'], sg['s'], np.maximum(sg['lam'], 0.0))
                integral = float(np.trapezoid(lam_t, rc['s'])) if hasattr(np, 'trapezoid') else float(np.trapz(lam_t, rc['s']))
                bounds.append((rc['E'][0] * np.exp(-integral), rc['E'][-1]))
            lines.append(f"| {m} | " + ", ".join(f"{x:g}" for x in splus) + " | " + ", ".join(f"{x:.1e}" for x in Erat) + " | "
                         + f"{min(x[0] for x in mins):.2f}, {min(x[1] for x in mins):.2f} | {np.nanmean(late_ratio):.5f} | {np.mean(late_lam):.5f} | {np.mean(late_al):.5f} | "
                         + f"{np.mean([k[0] for k in late_k]):.1e}, {np.mean([k[1] for k in late_k]):.1e} | " + "; ".join(f"{b:.1e} ≥ {e:.1e}" for b, e in bounds) + " |")
        # nearest-atom spurious minima
        lines += ["", "The nearest-atom contrast: E at the end of the nearest-atom runs against the E of the replicated configuration of their sampled weights (a critical point of the "
                  "nearest-atom energy; its drift under the nearest-atom leak and under the OT leak).", "",
                  "| m | sampled counts (seeds) | E(600) | E of the replicated configuration of the sampled weights | its drift: nearest-atom / OT leak |", "|---|---|---|---|---|"]
        for m in list(MS) + [M_REF]:
            rs = [r for r in rr if r['leak'] == 'nearest']
            cs = "; ".join(str(tuple(int(x) for x in r['spurious'][m]['counts'])) for r in rs if m in r['spurious'])
            Ee = "; ".join(f"{r['rec'][m]['E'][-1]:.2e}" for r in rs)
            Es = "; ".join(f"{r['spurious'][m]['E']:.2e}" for r in rs if m in r['spurious'])
            dr = "; ".join(f"{r['spurious'][m]['drift_nearest']:.0e} / {r['spurious'][m]['drift_ot']:.1e}" for r in rs if m in r['spurious'])
            lines.append(f"| {m} | {cs} | {Ee} | {Es} | {dr} |")
        c = consts[name]
        lines += ["", f"The a-priori basin of Theorem D‴.21(iv): sup λ_max R^G[μ*] over the box of the runs = {c['sup_kappa_R']:.4e}, so κ = σ₀ − sup = {c['kappa']:.5f}"
                  + (f" and E₀ = κ²/(2(L₂^G)²) = {c['E0']:.2e} (L₂^G = {consts['L2G']:.3f})." if np.isfinite(c['E0']) else " ≤ 0: hypothesis (H*) fails on that box, and only the a-posteriori form applies."), ""]
    rows, mruns = moved_fraction()
    lines += ["## The measure-level basin: particles moved into another cell (matching prior)", "",
              "The lift of q with k particles of class 0 moved 0.6 of the way to atom 1 (into cell 1); the same fractions p = k/m at m = 30 and 300. At s = 0: whether the plan keeps the "
              "moved particles in class 0, E, D_π, the current-residual curvature κ^Θ along the plan segments, the margin α and the rate λ (κ* = 0 at the matching prior).", "",
              "| m | k | p | kept in class 0 | E | D_π | κ^Θ | α | λ | observed ratio |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['m']} | {r['k']} | {r['p']:.4f} | {r['kept']} | {r['E']:.3e} | {r['D']:.3e} | {r['ktheta']:.3e} | {r['alpha']:+.5f} | {r['lam']:.5f} | {r['ratio']:.5f} |")
    lines += ["", "One particle moved, 300 slow units (midpoint rule, ds = 0.5, plan every step): ", ""]
    for (m, leak), r in mruns.items():
        bd = f" (a-posteriori bound {r['bound']:.3e})" if leak == 'ot' else " (E is the OT functional's excess, not this leak's Lyapunov function)"
        lines.append(f"- m = {m}, {leak}: E {r['E0']:.3e} → {r['E']:.3e}{bd}; observed ratio at the end {r['ratio']:.5f}; the moved particle labelled {r['lab0']} by the OT plan, "
                     f"in cell {r['cell']} at the end.")
    lines += ["", "Numbers are produced by `experiments/exp_ot_basin.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'ot_basin.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'ot_basin_raw.npz')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        report(list(z['runs']), z['consts'].item(), t0)
        return
    jobs = [(seed, name, leak) for seed in range(N_SEEDS) for name in SHIFTS for leak in ('ot', 'nearest')]
    with mp.Pool(N_WORKERS) as pool:
        runs = pool.map(job, jobs)
    consts = a_priori_constants(runs)
    np.savez(raw, runs=np.array(runs, dtype=object), consts=np.array(consts, dtype=object))
    report(runs, consts, t0)


if __name__ == '__main__':
    main()
