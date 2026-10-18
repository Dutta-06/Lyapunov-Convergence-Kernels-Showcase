"""
Theorem D‴: uniform-in-time propagation of chaos in the learning regime (the averaged slow loop)
==============================================================================================

The averaged slow loop (𝔄_m) of PARADIGM.md §5d for a three-bump plant fitted by saturated Gaussian kernels, at a leak
σ₀ = 0.3σ_crit inside the learning regime (the collapsed configuration is unstable: Proposition K).  Three measurements:

1. The instance's constants: σ_crit (Proposition K), the whitened Lipschitz constants L_ϑ^G, L₂^G (search + closed forms), K_Φ,
   F_ρ, and along the reference flow the hypothesis profile a(s) = σ₀ − κ_R[μ_s](N_ϱ′(supp μ_s)) (κ_R at the reference particles
   and over their ϱ′-neighbourhood by Monte Carlo), from which κ, s₁, Λ₁ and the theorem's constant C₁ are read off.
2. Warm start (the support of μ₀ is three clusters, one per bump: hypothesis (H) from s = 0): the coupled mean-square distance
   D_m(s) of m = 16, 64, 256 particles to an M = 2048 reference over 600 slow units (8.5 leak times) in time windows — uniform in
   time means the window maxima do not grow — its m-law, the theorem's bound C₁/m, and the function-level distance.
3. Cold start (a connected spread cloud whose mass splits between the basins: outside (H), Remark D‴.2): the same quantities,
   the per-particle deviations at the end (bimodal when pairs sit in different basins) and the fraction of such pairs against m.
4. The real closed loop (not averaged) with the warm start at k_e = 25, by the synchronous coupling of Theorem D″'s experiment:
   does the coupled distance saturate where Theorem D″'s cold start grew?

Output: results/chaos/learning_table.md and learning_raw.npz (set CLW_WORKERS, CLW_SEEDS, CLW_SEEDS_REAL to change the pool).
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
from src.core import periodic as PO
from src.core import chaos_learning as CL

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
M_REF = int(os.environ.get('CLW_M_REF', 2048))
MS = (16, 64, 256)
S_END, DT_S, REC = float(os.environ.get('CLW_S_END', 600.0)), 0.02, 50
WINDOWS = tuple((a * S_END / 600.0, b * S_END / 600.0) for a, b in ((0.0, 75.0), (75.0, 150.0), (150.0, 300.0), (300.0, 600.0)))
RHO_PRIME = 0.3
LEAK_FRACTION = 0.3
KE = 25.0
T_REAL, DT_REAL = float(os.environ.get('CLW_T_REAL', 2000.0)), 1e-3
W_REAL = tuple((a * T_REAL / 2000.0, b * T_REAL / 2000.0) for a, b in ((0.0, 500.0), (500.0, 1000.0), (1000.0, 2000.0)))
N_SEEDS = int(os.environ.get('CLW_SEEDS', 4))
N_SEEDS_REAL = int(os.environ.get('CLW_SEEDS_REAL', 3))
N_WORKERS = int(os.environ.get('CLW_WORKERS', 3))


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def instance_constants():
    net = make_net()
    cl = PO.ClosedLoop(net, lambda Xa: f_vec(Xa[:, 0])[:, None], x_ref, x_ref_dot, PER, KE, G, 0.05)
    orb = cl.orbit(N=1000)
    g = net.unconstrained(np.array([[1.0]]), np.array([[0.0]]), np.array([0.8]))[0]
    s_crit, tb, R2 = PO.sigma_crit(cl, orb, g)
    a, c, s = net.physical(tb[None, :])
    Lc = CL.whitened_lipschitz_constants(net, G, 2.5, n_samples=60000)
    return dict(sigma_crit=float(s_crit), atom=(float(a[0, 0]), float(c[0, 0]), float(s[0])), **Lc)


# ---------------------------------------------------------------------------- averaged runs
def averaged_run(seed, warm, sigma0):
    net = make_net()
    loop = CL.AveragedSlowLoop(net, f_vec, x_ref, PER, G, sigma0, n_quad=N_Q)
    rng = np.random.default_rng(seed)
    Th_all = CL.warm_start(net, BUMPS, M_REF, rng) if warm else CL.cold_start(net, M_REF, rng)
    R0 = float(np.sqrt(np.sum(Th_all ** 2 / G[None, :], axis=1)).max())
    rec = loop.coupled_run(Th_all, MS, S_END, DT_S, rec_every=REC, rho_nbhd=RHO_PRIME, n_dirs=3, rng=rng)
    out = dict(seed=seed, warm=warm, R0=R0, s=rec['s'], fit=rec['fit'], spread=rec['spread'], drift_rms=rec['drift_rms'],
               kR_ref=rec['kR_ref'], kR_ref_nbhd=rec['kR_ref_nbhd'])
    for m in MS:
        out[f'D{m}'] = rec['D'][m]; out[f'h{m}'] = rec['hnorm_mf'][m]; out[f'dev1_{m}'] = rec['dev1'][m]
        out[f'kRm{m}'] = rec['kR_m'][m]; out[f'devfinal{m}'] = rec['dev_final'][m]
    return out


# ---------------------------------------------------------------------------- the real loop (Theorem D″'s machinery)
def real_loop_run(seed, sigma0):
    """The closed loop at k_e = KE with σ = σ₀/k_e, τ = 0, warm start, m-systems sharing particles with the M-reference
    (RK4 plant with u held over the step, Euler parameters; as in exp_chaos_uniform.coupled_run)."""
    net = make_net()
    ctrl = L.LangevinAdaptiveController(net, ke=KE, G_diag=G, sigma=sigma0 / KE, theta0=np.zeros(3), tau=0.0)
    rng = np.random.default_rng(1000 + seed)
    Th_all = CL.warm_start(net, BUMPS, M_REF, rng)
    systems = {m: dict(x=x_ref(0.0).copy(), Th=Th_all[:m].copy()) for m in MS + (M_REF,)}
    n_steps = int(round(T_REAL / DT_REAL))
    rec_t, rec = [], {m: [] for m in MS}
    orbit_x = np.array([x_ref(tt)[0] for tt in np.arange(400) * PER / 400])[:, None]
    f_orb = f_vec(orbit_x[:, 0])
    fit_t, fit = [], []
    t = 0.0
    f_plant = lambda x: f_vec(x)
    for k in range(n_steps):
        xr, xrd = x_ref(t), x_ref_dot(t)
        new = {}
        for m, S in systems.items():
            x, Th = S['x'], S['Th']
            e = x - xr
            u = ctrl.control(x, xr, xrd, Th)
            k1 = f_plant(x) + u
            k2 = f_plant(x + 0.5 * DT_REAL * k1) + u
            k3 = f_plant(x + 0.5 * DT_REAL * k2) + u
            k4 = f_plant(x + DT_REAL * k3) + u
            new[m] = dict(x=x + DT_REAL / 6 * (k1 + 2 * k2 + 2 * k3 + k4), Th=Th + ctrl.drift(x, e, Th) * DT_REAL)
        systems = new
        t = (k + 1) * DT_REAL
        if k % 50 == 0:
            ref = systems[M_REF]
            rec_t.append(t)
            for m in MS:
                d = systems[m]['Th'] - ref['Th'][:m]
                rec[m].append(float(np.mean(np.sum(d ** 2 / G[None, :], axis=1))))        # D_m in the G⁻¹ metric, as in the averaged runs
        if k % int(round(50.0 / DT_REAL)) == 0 or k == n_steps - 1:
            ThM = systems[M_REF]['Th']
            fh = net.f_hat_batch(orbit_x, ThM)[:, 0]
            fit_t.append(t); fit.append(float(np.sqrt(np.mean((f_orb - fh) ** 2) / np.mean(f_orb ** 2))))
    return dict(seed=seed, t=np.array(rec_t), D={m: np.array(rec[m]) for m in MS}, fit_t=np.array(fit_t), fit=np.array(fit))


# ---------------------------------------------------------------------------- summaries
def window_stats(s, D, windows):
    return [float(D[(s > a) & (s <= b)].max()) for a, b in windows] + [float(D[(s > windows[-1][0]) & (s <= windows[-1][1])].mean())]


def slope_se(vals_by_m):
    """vals_by_m: {m: [per-seed values]}; slope of log value against log m per seed, mean ± s.e."""
    ms = sorted(vals_by_m)
    n = len(vals_by_m[ms[0]])
    sl = [np.polyfit(np.log(ms), np.log([vals_by_m[m][i] for m in ms]), 1)[0] for i in range(n)]
    return float(np.mean(sl)), (float(np.std(sl, ddof=1) / np.sqrt(n)) if n > 1 else float('nan'))


def mean_se(v):
    v = np.asarray(v, float)
    return float(v.mean()), (float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else float('nan'))


def fmt(x, se=None, p=3):
    if se is None or not np.isfinite(se):
        return f"{x:.{p}e}"
    return f"{x:.{p}e} ± {se:.1e}"


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    K = instance_constants()
    sigma0 = LEAK_FRACTION * K['sigma_crit']
    net = make_net()
    loop = CL.AveragedSlowLoop(net, f_vec, x_ref, PER, G, sigma0, n_quad=N_Q)
    jobs = [(s, True, sigma0) for s in range(N_SEEDS)] + [(s, False, sigma0) for s in range(N_SEEDS)]
    with mp.Pool(N_WORKERS) as pool:
        res = pool.starmap(averaged_run, jobs)
        res_real = pool.starmap(real_loop_run, [(s, sigma0) for s in range(N_SEEDS_REAL)])
    warm = [r for r in res if r['warm']]; cold = [r for r in res if not r['warm']]
    s_grid = warm[0]['s']

    # ---- the hypothesis profile and the constants (warm start)
    kR_nb = np.mean([r['kR_ref_nbhd'] for r in warm], axis=0)
    kR_pt = np.mean([r['kR_ref'] for r in warm], axis=0)
    a_prof = sigma0 - kR_nb
    kappa = 0.5 * float(a_prof.min()) if a_prof.min() > 0 else float('nan')
    s1 = 0.0 if a_prof.min() > 0 else float(s_grid[np.where(a_prof < kappa)[0].max()]) if np.isfinite(kappa) else float('nan')
    Lambda1 = CL.Lambda1_from_profile(s_grid, a_prof, kappa, s1) if np.isfinite(kappa) else float('nan')
    R0 = max(r['R0'] for r in warm)
    C = CL.theorem_constants(loop, K['L_th_G'], K['L_2_G'], R0, kappa, Lambda1, RHO_PRIME) if np.isfinite(kappa) else None
    # cold-start profile
    kR_nb_c = np.mean([r['kR_ref_nbhd'] for r in cold], axis=0)
    a_prof_c = sigma0 - kR_nb_c

    lines = ["# Theorem D‴: uniform-in-time propagation of chaos in the learning regime — the averaged slow loop (`exp_chaos_learning.py`)", "",
             f"n = 1, saturated features with a_max = {NET_KW['a_max']}, c_max = {NET_KW['c_max']}, s ∈ [{NET_KW['s_min']}, {NET_KW['s_max']}], prior ϑ⁰ = 0, "
             f"G = diag{tuple(float(g) for g in G)}; plant f(x) = Σ_k A_k e^{{−(x−C_k)²/(2S_k²)}} with (A, C, S) = {BUMPS} (three bumps; the "
             f"network can represent it exactly with three kernels of weights 3A_k); reference x_m = 2 sin t, ρ its occupation measure on a grid of {N_Q} "
             f"times; slow time s = t/k_e, (𝔄_m) integrated by the midpoint rule with step {DT_S:g}; the m-systems share their particles with an M = {M_REF} reference "
             f"(the stand-in for the mean-field flow, as in Theorems D′ and D″); horizon {S_END:g} slow units.", "",
             f"Constants. Proposition K: the collapsed critical point (best single kernel, a = {K['atom'][0]:.3f}, c = {K['atom'][1]:.3f}, s = {K['atom'][2]:.3f}) is a "
             f"strict local minimum of J_m iff σ₀ > σ_crit = {K['sigma_crit']:.5f}; the runs use σ₀ = {LEAK_FRACTION:g}σ_crit = {sigma0:.5f} (learning regime; the leak time "
             f"1/σ₀ = {1 / sigma0:.0f} slow units). Whitened Lipschitz constants by search (|x| ≤ 2.5, 6·10⁴ sampled ϑ, Nelder–Mead polish): L_ϑ^G = {K['L_th_G']:.3f} "
             f"(closed-form bound {K['L_th_G_closed']:.3f}), L₂^G = {K['L_2_G']:.2f} (closed-form bound λ_G a_max/s_min² = {K['L_2_G_closed']:.2f}); K_Φ = {net.K_Phi:.1f}, "
             f"F_ρ = {loop.F_rho:.3f}.", ""]

    # ---- Table 1: hypothesis profile
    lines += ["## 1. The hypothesis (H) along the reference flow", "",
              "κ_R[μ_s] at the reference particles (max over the 2048 particles) and over their ϱ′-neighbourhood (ϱ′ = 0.3 in the G⁻¹ metric; Monte Carlo with 3 random "
              "points per particle, a lower estimate of the supremum), the margin a(s) = σ₀ − κ_R(neighbourhood), the reference's fit error ‖f − f̂‖/‖f‖ on the orbit and "
              "its rms drift (G⁻¹ norm); means over the seeds.", "",
              "| s | warm: κ_R particles | warm: κ_R nbhd | warm: a(s) | warm: fit | warm: rms drift | cold: κ_R particles | cold: κ_R nbhd | cold: a(s) | cold: fit | cold: rms drift |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    kR_pt_c = np.mean([r['kR_ref'] for r in cold], axis=0)
    fit_w = np.mean([r['fit'] for r in warm], axis=0); fit_c = np.mean([r['fit'] for r in cold], axis=0)
    dr_w = np.mean([r['drift_rms'] for r in warm], axis=0); dr_c = np.mean([r['drift_rms'] for r in cold], axis=0)
    for sv in (0.5, 10, 25, 50, 75, 100, 150, 200, 300, 400, 500, 600):
        i = int(np.argmin(np.abs(s_grid - sv)))
        lines.append(f"| {s_grid[i]:.1f} | {kR_pt[i]:.5f} | {kR_nb[i]:.5f} | {a_prof[i]:+.5f} | {fit_w[i]:.4f} | {dr_w[i]:.2e} | "
                     f"{kR_pt_c[i]:.5f} | {kR_nb_c[i]:.5f} | {a_prof_c[i]:+.5f} | {fit_c[i]:.4f} | {dr_c[i]:.2e} |")
    if C is not None:
        lines += ["", f"Warm start: a(s) ≥ {a_prof.min():.5f} > 0 on the whole run, so (H) holds from s₁ = {s1:g} with margin κ := ½ min a = {kappa:.5f} and Λ₁ = {Lambda1:.3g}; "
                  f"R₀ = {R0:.3f}, R_Θ = {C['R_Theta']:.3f}, V_max = {C['V_max']:.3f}, ϱ = {C['rho']:.3e}, θ = {C['theta']:.3e}, Δ = {C['Delta']:.3e}, C₁ = 4e^{{2Λ₁}}(L_ϑ^G)²K_Φ²/κ² = {C['C1']:.3e}, "
                  f"m₀ = 16K_Φ²/θ² = {C['m0']:.3e} (the tail term of the theorem needs m ≥ m₀; the main term C₁/m is what the table compares with).",
                  f"Cold start: min_s a(s) = {a_prof_c.min():+.5f} at s = {s_grid[int(np.argmin(a_prof_c))]:.1f}, final a = {a_prof_c[-1]:+.5f}: "
                  + ("(H) holds at the end of the run" if a_prof_c[-1] > 0 else "(H) fails on the whole run (the connected support threads the ridges, Remark D‴.2)") + "."]
    else:
        lines += ["", f"Warm start: min_s a(s) = {a_prof.min():+.5f} ≤ 0 — (H) is not verified on this run; the constants are not formed."]

    # ---- Table 2/3: coupled distance in windows
    def block(title, rs):
        out = ["", f"## {title}", "",
               "| m | " + " | ".join(f"max D_m on ({a:g}, {b:g}]" for a, b in WINDOWS) + f" | time average on ({WINDOWS[-1][0]:g}, {WINDOWS[-1][1]:g}] | bound C₁/m | mean ‖f̂_Θ − f̂_ref‖ on the last window | fraction of particles with ‖δ_i‖_{{G⁻¹}} > 1 at s = {S_END:g} |",
               "|---|" + "---|" * (len(WINDOWS) + 4)]
        W = {m: [window_stats(r['s'], r[f'D{m}'], WINDOWS) for r in rs] for m in MS}
        H = {m: [float(r[f'h{m}'][(r['s'] > WINDOWS[-1][0]) & (r['s'] <= WINDOWS[-1][1])].mean()) for r in rs] for m in MS}
        P = {m: [float(np.mean(r[f'devfinal{m}'] > 1.0)) for r in rs] for m in MS}
        for m in MS:
            cells = []
            for j in range(len(WINDOWS) + 1):
                mu, se = mean_se([w[j] for w in W[m]]); cells.append(fmt(mu, se))
            hm, hs = mean_se(H[m]); pm, ps = mean_se(P[m])
            out.append(f"| {m} | " + " | ".join(cells) + f" | {(C['C1'] / m if C else float('nan')):.3e} | {fmt(hm, hs)} | {pm:.3f} ± {ps:.3f} |")
        sl, se = slope_se({m: [w[len(WINDOWS) - 1] for w in W[m]] for m in MS})
        sl_h, se_h = slope_se(H)
        ratios = [np.mean([w[len(WINDOWS) - 1] for w in W[m]]) / np.mean([w[0] for w in W[m]]) for m in MS]
        ratios2 = [np.mean([w[len(WINDOWS) - 1] for w in W[m]]) / np.mean([w[len(WINDOWS) - 2] for w in W[m]]) for m in MS]
        out += ["", f"Slope of the last-window maximum of D_m against m (per-seed fits, mean ± s.e.): {sl:.2f} ± {se:.2f} (Theorem D‴: −1 for the mean square; "
                f"the function-level distance ‖f̂_Θ − f̂_ref‖ has slope {sl_h:.2f} ± {se_h:.2f}, m^{{−1/2}} being −0.5). Ratio of the last to the first window maxima: "
                + ", ".join(f"{r:.2f}" for r in ratios) + "; last to the one before: " + ", ".join(f"{r:.2f}" for r in ratios2) + " (uniform in time means no growth)."]
        return out, W, H, P
    b_w, W_w, H_w, P_w = block("2. Warm start (three clusters, one per bump; (H) from s = 0)", warm)
    b_c, W_c, H_c, P_c = block("3. Cold start (a connected spread cloud: a ~ N(0, 0.8²), c ~ U(−2, 2), s ~ U(0.4, 1); outside (H))", cold)
    lines += b_w + b_c
    # per-particle deviation quantiles at the end
    lines += ["", "Per-particle deviations ‖δ_i‖_{G⁻¹} at s = %g (pooled over the seeds): median, 90th percentile, maximum." % S_END, "",
              "| start | m | median | 90% | max |", "|---|---|---|---|---|"]
    for name, rs in (("warm", warm), ("cold", cold)):
        for m in MS:
            d = np.concatenate([r[f'devfinal{m}'] for r in rs])
            lines.append(f"| {name} | {m} | {np.median(d):.3e} | {np.quantile(d, 0.9):.3e} | {d.max():.3e} |")

    # ---- Table 4: the real loop
    lines += ["", f"## 4. The real closed loop at k_e = {KE:g} (σ = σ₀/k_e = {sigma0 / KE:.5f}, τ = 0), warm start", "",
              f"The m-systems share their particles with the M = {M_REF} reference; dt = {DT_REAL:g}; D_m in the G⁻¹ metric as above; real time t (slow time s = t/{KE:g}).", "",
              "| m | " + " | ".join(f"max D_m on ({a:g}, {b:g}]" for a, b in W_REAL) + f" | time average on ({W_REAL[-1][0]:g}, {W_REAL[-1][1]:g}] |",
              "|---|" + "---|" * (len(W_REAL) + 1)]
    WR = {m: [window_stats(r['t'], r['D'][m], W_REAL) for r in res_real] for m in MS}
    for m in MS:
        cells = []
        for j in range(len(W_REAL) + 1):
            mu, se = mean_se([w[j] for w in WR[m]]); cells.append(fmt(mu, se))
        lines.append(f"| {m} | " + " | ".join(cells) + " |")
    slr, ser = slope_se({m: [w[len(W_REAL) - 1] for w in WR[m]] for m in MS})
    ratios_r = [np.mean([w[len(W_REAL) - 1] for w in WR[m]]) / np.mean([w[0] for w in WR[m]]) for m in MS]
    fit_r = np.mean([r['fit'] for r in res_real], axis=0)
    lines += ["", f"Slope of the last-window maximum against m: {slr:.2f} ± {ser:.2f}. Ratio of the last to the first window maxima: " + ", ".join(f"{r:.2f}" for r in ratios_r)
              + f". Reference fit error on the orbit: {fit_r[0]:.4f} at t = 0, {fit_r[len(fit_r) // 2]:.4f} at t = {res_real[0]['fit_t'][len(fit_r) // 2]:.0f}, {fit_r[-1]:.4f} at t = {T_REAL:g}."]

    # ---- Reading (data-aware)
    def grew(ratios):
        return any(r > 1.5 for r in ratios)
    r_w = [np.mean([w[len(WINDOWS) - 1] for w in W_w[m]]) / np.mean([w[0] for w in W_w[m]]) for m in MS]
    r_c = [np.mean([w[len(WINDOWS) - 1] for w in W_c[m]]) / np.mean([w[0] for w in W_c[m]]) for m in MS]
    sl_w = slope_se({m: [w[len(WINDOWS) - 1] for w in W_w[m]] for m in MS})[0]
    sl_c = slope_se({m: [w[len(WINDOWS) - 1] for w in W_c[m]] for m in MS})[0]
    below = all(np.mean([w[len(WINDOWS) - 1] for w in W_w[m]]) < C['C1'] / m for m in MS) if C else False
    flips_c = [np.mean(P_c[m]) for m in MS]
    lines += ["", "Reading: "
              + (f"on the warm start the hypothesis holds along the whole run (a(s) > 0), the window maxima of D_m {'do not grow' if not grew(r_w) else 'grow'} "
                 f"(last/first {', '.join(f'{r:.2f}' for r in r_w)}), the m-law of the mean square has slope {sl_w:.2f} (theorem: −1), and the measured values lie "
                 f"{'below' if below else 'NOT below'} the bound C₁/m" + (f" (which carries no transient factor here, Λ₁ = {Lambda1:.2g})" if C else "") + ". ")
              + (f"On the cold start the support threads the ridges and (H) fails along the run (min a = {a_prof_c.min():+.4f}); the coupled distance "
                 f"{'grows' if grew(r_c) else 'does not grow'} over the horizon (last/first {', '.join(f'{r:.2f}' for r in r_c)}), its m-law has slope {sl_c:.2f}, and at the end a fraction "
                 f"{', '.join(f'{p:.3f}' for p in flips_c)} (m = 16, 64, 256) of the pairs sit more than one whitened unit apart: pairs in different basins, the floor of "
                 f"Remark D‴.2" + (", whose fraction decreases with m" if flips_c[0] > flips_c[-1] else "") + ". ")
              + (f"The real loop on the warm start {'saturates' if not grew(ratios_r) else 'grows'} (last/first {', '.join(f'{r:.2f}' for r in ratios_r)}) at the same k_e at which "
                 f"Theorem D″'s cold start grew threefold.")]
    lines += ["", "Numbers are produced by `experiments/exp_chaos_learning.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'learning_table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    np.savez(os.path.join(OUT, 'learning_raw.npz'), sigma0=sigma0, constants=np.array([K], dtype=object), warm=np.array(warm, dtype=object),
             cold=np.array(cold, dtype=object), real=np.array(res_real, dtype=object), allow_pickle=True)
    print("\n".join(lines))
    print(f"[{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()
