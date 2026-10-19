"""
Proposition D‴.30: measurement noise at perfect tracking — a coherent, high-passed noise; the floors of the fixed and the adaptive prior
=======================================================================================================================================

The three-bump instance (σ₀ = 0.15σ_crit, x_m = 2 sin t), the real closed loop with the nearest-atom leak toward a fixed (η = 0) or an
adaptive prior (η > 0), at the lift of the exact representation of the matching weights q = 1/3 with the prior at its atoms (perfect
tracking), and the measurement y = x + ν used by the controller and the adaptive law (src/core/measurement_noise.py).

A. The linear theory (no dynamics): H = ⟨ĴᵀĴ⟩ and Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩ in the eigenbasis of H (a = ∂_xf̂ on the orbit); per mode the
   function-level floor's share (N₀/2k_e)Q_jj(h_j + η)/(h_j + η + σ₀) and the parameter variance at η = 0, σ₀, 10σ₀, ∞, with the
   bias-removal rate λ_j^−(η); the second-order (errors-in-variables) force per unit of the noise's variance and the stationary mean
   deviations it produces, (H + σ₀I)⁻¹ against H⁻¹ times k_e times the force.
B. The exact linearised loop driven by white measurement noise: its periodic stationary covariance against the averaged prediction, at
   k_e = 25, 50, 100, 200 and η = 0, σ₀, ∞ (the function-level floor ⟨Ĵ P_yy Ĵᵀ⟩ over a period).
C. The real nonlinear loop (the weighted K₀-atom loop; RK4, dt = 10⁻³, k_e = 25) with an Ornstein–Uhlenbeck measurement noise (rate
   B = 100, exact at the steps), N_MC independent runs per case, T_C time units, at two noise levels: E‖f̂ − f‖²_ρ(t) against the linear
   prediction — the transient covariance of the exactly linearised loop with the OU noise, plus the second-order mean deviation — for the
   fixed and the adaptive prior (η = ∞). And one run of the lift m = 15 driven by the same noise path as a K₀-atom run: the same f̂ to
   rounding (the noise is coherent: a replicated start stays replicated).

Output: results/chaos/measurement_noise.md and measurement_noise_raw.npz (MN_WORKERS, MN_RUNS, MN_TC; `--from-raw` regenerates the md,
`--only-A` recomputes Part A and merges it).
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
from src.core import adaptive_prior as AP
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q3 = np.ones(3) / 3
KE = 25.0
ETAS_A = (0.0, SIGMA0, 10 * SIGMA0, np.inf)
KES_B = (25.0, 50.0, 100.0, 200.0)
ETAS_B = (0.0, SIGMA0, np.inf)
B_OU, DT_C = 100.0, 1e-3
STDS_C = (0.003, 0.01)
N_MC = int(os.environ.get('MN_RUNS', 12))
T_C = float(os.environ.get('MN_TC', 300.0))
REC_C = 25.0
N_WORKERS = int(os.environ.get('MN_WORKERS', 3))
N_GRID = 256


def f_vec(X):
    X = np.asarray(X, float)
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def fprime(x):
    return sum(-A * (x - C) / S ** 2 * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def fsecond(x):
    return sum(A * ((x - C) ** 2 / S ** 4 - 1 / S ** 2) * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


xr = lambda t: np.array([2.0 * np.sin(t)])
xrd = lambda t: np.array([2.0 * np.cos(t)])


def setup():
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    A = np.array([b[0] for b in BUMPS]); C = np.array([b[1] for b in BUMPS]); S = np.array([b[2] for b in BUMPS])
    atoms = net.unconstrained((A / Q3)[:, None], C[:, None], S)
    return net, atoms


def save_raw(path, **arrays):
    tmp = path + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------- A
def part_A():
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, np.inf, atoms, Q3)
    H, Q = MN.noise_matrices(cv, fprime)
    h, V = np.linalg.eigh(H)
    qd = np.einsum('ij,ik,kj->j', V, Q, V)
    t = np.arange(1024) * PER / 1024
    J = np.array([cv.J(s)[0] for s in t]); k = np.fft.fftfreq(1024, d=1.0 / 1024)
    Jd = np.real(np.fft.ifft((1j * k)[:, None] * np.fft.fft(J, axis=0), axis=0))
    a = np.array([fprime(2 * np.sin(s)) for s in t])
    parts = dict(deriv=float(np.trace(Jd.T @ Jd) / 1024), cross=float(2 * np.trace(Jd.T @ (a[:, None] * J)) / 1024),
                 slope=float(np.trace((a[:, None] * J).T @ (a[:, None] * J)) / 1024))
    rows = {}
    for eta in ETAS_A:
        fl, sh, pv = MN.floor_closed_form(H, Q, SIGMA0, eta, 1.0, KE)
        lm = (h + SIGMA0) * 0 if eta == 0 else AP.pair_rates(h, SIGMA0, eta)[0]
        rows[eta] = dict(floor=fl, shares=sh, pvar=pv, rate=lm)
    F1 = MN.eiv_force(cv, fprime, fsecond, KE, B_OU, 1.0)          # per unit of the noise's variance (real time, whitened)
    Fs = KE * F1                                                    # slow time
    fv = V.T @ Fs
    bias = {0.0: float(np.sum(fv ** 2 * h / (h + SIGMA0) ** 2)), np.inf: float(np.sum(fv ** 2 / h))}
    trH2, J4 = MN.fast_parts(cv)
    print("[A]", flush=True)
    return dict(H=H, Q=Q, h=h, V=V, qd=qd, parts=parts, rows=rows, F1=F1, fv=fv, bias=bias, trH2=trH2, J4=J4)


# ---------------------------------------------------------------------------- B
def job_B(args):
    ke, eta = args
    t0 = time.time()
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, eta if eta > 0 else np.inf, atoms, Q3)
    if eta == 0:
        gen = MN.fixed_generator(cv, ke); D = 1 + cv.d
    else:
        gen = lambda t: cv.generator(t, ke); D = cv.dim
    b0 = MN.collapsed_noise_input(cv, fprime, ke)
    bv = (lambda t: b0(t)[:D])
    P0 = MN.periodic_covariance(gen, bv, PER, D, 1.0)
    # the function-level floor over a period: ⟨Ĵ(t)P_yy(t)Ĵ(t)ᵀ⟩ along the periodic solution
    ts = np.linspace(0, PER, 65)
    tt, Ps = MN.transient_covariance(gen, bv, D, 1.0, PER, ts, P0=P0, rtol=1e-10)
    H, Q = MN.noise_matrices(cv, fprime)
    trH2, J4 = MN.fast_parts(cv)
    point = [cv.J(t)[0] @ P[1:1 + cv.d, 1:1 + cv.d] @ cv.J(t)[0] for t, P in zip(tt, Ps)]
    orbit = [np.trace(H @ P[1:1 + cv.d, 1:1 + cv.d]) for P in Ps]
    exact_point = float(np.trapezoid(point, tt) / PER); exact_orbit = float(np.trapezoid(orbit, tt) / PER)
    slow = MN.floor_closed_form(H, Q, SIGMA0, eta, 1.0, ke)[0]
    pred_orbit, pred_point = MN.floors(H, Q, trH2, J4, SIGMA0, eta, 1.0, ke)
    print(f"[{time.time() - t0:.0f}s] B ke {ke:g} eta {eta:g}", flush=True)
    return dict(ke=ke, eta=eta, exact_orbit=exact_orbit, exact_point=exact_point, slow=slow, pred_orbit=pred_orbit, pred_point=pred_point)


# ---------------------------------------------------------------------------- C
def job_C(args):
    kind, eta, std, seed = args
    t0 = time.time()
    net, atoms = setup()
    grid = 2 * np.sin(np.arange(N_GRID) * PER / N_GRID); Fg = f_vec(grid)
    if kind == 'lift':
        labels, weights, m = np.repeat(np.arange(3), 5), None, 15
    else:
        labels, weights, m = np.arange(3), Q3, 3
    cl = MN.NoisyRealLoop(net, f_vec, xr, xrd, PER, KE, G, SIGMA0, labels, 3, eta, weights=weights)
    z0 = cl.lift(atoms)

    def rec(t, z):
        Th = z[1:1 + m * 3].reshape(m, 3)
        fh = cl.f_hat(grid[:, None], Th) if False else np.array([cl.f_hat(np.array([x]), Th)[0] for x in grid])
        return np.concatenate([[np.mean((fh - Fg) ** 2)], fh])
    out = MN.simulate_noisy(cl, z0, m, T_C, DT_C, B_OU, std ** 2, np.random.default_rng(1000 + seed), int(round(REC_C / DT_C)), rec)
    out = np.array(out)
    print(f"[{time.time() - t0:.0f}s] C {kind} eta {eta:g} std {std:g} seed {seed}", flush=True)
    return dict(kind=kind, eta=eta, std=std, seed=seed, err=out[:, 0], fh=out[:, 1:].astype(np.float32) if kind == 'lift' or seed == 0 else None)


def predictions_C(std, eta):
    """The linear prediction of E‖f̂ − f‖²_ρ at the records: the transient covariance of the exactly linearised loop with the OU noise (from
    the stationary noise and P = 0) plus the second-order mean deviation."""
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, np.inf, atoms, Q3)
    var = std ** 2; N0 = 2 * var / B_OU
    if eta == 0:
        gen = MN.fixed_generator(cv, KE); D = 1 + cv.d
    else:
        gen = lambda t: cv.generator(t, KE); D = cv.dim
    b0 = MN.collapsed_noise_input(cv, fprime, KE)
    Ma, cfun = MN.ou_augmented(gen, lambda t: b0(t)[:D], D, B_OU, N0)
    P0 = np.zeros((D + 1, D + 1)); P0[-1, -1] = var
    tev = np.arange(0.0, T_C + 1e-9, REC_C)
    tt, Ps = MN.transient_covariance(Ma, cfun, D + 1, N0, T_C, tev, P0=P0)
    H, Q = MN.noise_matrices(cv, fprime)
    cov = np.array([np.trace(H @ P[1:1 + cv.d, 1:1 + cv.d]) for P in Ps])
    F = KE * MN.eiv_force(cv, fprime, fsecond, KE, B_OU, var)
    mean = np.array([(lambda y: y @ H @ y)(MN.mean_response(H, F, SIGMA0, eta, t / KE)) for t in tev])
    return dict(t=tev, cov=cov, mean=mean)


# ---------------------------------------------------------------------------- report
def report(A, Bres, Cres, Cpred, t0=None):
    es = lambda e: '∞' if np.isinf(e) else ('0 (fixed)' if e == 0 else f'{e:.5f}')
    h, qd = A['h'], A['qd']
    lines = ["# Proposition D‴.30: measurement noise at perfect tracking (`exp_measurement_noise.py`)", "",
             f"The three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f}, x_m = 2 sin t, the nearest-atom leak toward a fixed (η = 0) or an adaptive prior, at the "
             "lift of the exact representation of q = 1/3 with the prior at its atoms; the measurement y = x + ν in the controller and the adaptive law.", "",
             "## A. The linear theory and the second-order force (no dynamics)", "",
             f"k_e = {KE:g}, N₀ = 1. Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩, a = ∂_xf̂ on the orbit: tr Q = {np.trace(A['Q']):.4f} = {A['parts']['deriv']:.4f} (⟨Ĵ̇ᵀĴ̇⟩) "
             f"{A['parts']['cross']:+.4f} (cross) {A['parts']['slope']:+.4f} (⟨a²ĴᵀĴ⟩). Per mode (the eigenbasis of H): h_j, Q_jj, Q_jj/h_j, and the "
             "function-level shares of the floor, (N₀/2k_e)Q_jj(h_j + η)/(h_j + η + σ₀), for each η.", "",
             "| j | h_j | Q_jj | Q_jj/h_j | " + " | ".join(f"share at η = {es(e)}" for e in ETAS_A) + " |", "|---|---|---|---|" + "---|" * len(ETAS_A)]
    for j in range(len(h)):
        lines.append(f"| {j + 1} | {h[j]:.4e} | {qd[j]:.4e} | {qd[j] / h[j]:.2f} | " + " | ".join(f"{A['rows'][e]['shares'][j]:.4e}" for e in ETAS_A) + " |")
    lines += ["", "| η | floor (N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) | relative to η = ∞ | the weakest mode's parameter variance | its bias-removal rate λ_1^− |",
              "|---|---|---|---|---|"]
    finf = A['rows'][np.inf]['floor']
    for e in ETAS_A:
        r = A['rows'][e]
        lines.append(f"| {es(e)} | {r['floor']:.5e} | {r['floor'] / finf:.4f} | {r['pvar'][0]:.4e} | {r['rate'][0]:.3e} |")
    lines += ["", f"At η = ∞ the floor is (N₀/2k_e)tr Q = {np.trace(A['Q']) / (2 * KE):.5e} whatever H. The second-order force (errors in variables; OU noise of rate "
              f"B = {B_OU:g}), slow time, per unit of the noise's variance, in the eigenbasis of H: " + ", ".join(f"{x:+.3e}" for x in A['fv']) + ". The stationary "
              f"mean deviation's function-level size per unit of the variance squared: fixed prior Σ_jF_j²h_j/(h_j + σ₀)² = {A['bias'][0.0]:.4e}, adaptive prior "
              f"(every η > 0) Σ_jF_j²/h_j = {A['bias'][np.inf]:.4e} (ratio {A['bias'][np.inf] / A['bias'][0.0]:.2f})."]
    if 'trH2' in A:
        xo = {e: A['rows'][e]['floor'] * 2.0 / B_OU / A['bias'][e] for e in (0.0, np.inf)}
        lines += ["", f"The fast parts: tr H² = {A['trH2']:.4f} (the function-level floor) and ⟨|Ĵ|⁴⟩ = {A['J4']:.4f} (the compensation error at the current point), "
                  f"each times N₀/(2k_e). With the OU noise (N₀ = 2Eν²/B) the errors-in-variables bias equals the floor of the parameters at Eν² = "
                  f"{xo[0.0]:.2e} (std {np.sqrt(xo[0.0]):.4f}) under the fixed prior and {xo[np.inf]:.2e} (std {np.sqrt(xo[np.inf]):.4f}) under the adaptive one, "
                  "and dominates above."]
    if Bres:
        lines += ["", "## B. The exactly linearised loop with white measurement noise against the averaged prediction", "",
                  "The periodic stationary covariance P(t) (the monodromy and a discrete Lyapunov equation; N₀ = 1), averaged over a period: the function-level floor "
                  "E‖f̂ − f‖²_ρ = ⟨tr(HP_yy)⟩ and the compensation error at the plant's current point ⟨ĴP_yyĴᵀ⟩, against the prediction — the slow (parameter) part "
                  "(N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) plus the fast part (N₀/2k_e)tr H² and (N₀/2k_e)⟨|Ĵ|⁴⟩ respectively.", "",
                  "| k_e | η | E‖f̂ − f‖²_ρ: exact | predicted | ratio | at the current point: exact | predicted | ratio | the slow part alone |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(Bres, key=lambda r: (r['eta'] if np.isfinite(r['eta']) else 1e9, r['ke'])):
            lines.append(f"| {r['ke']:g} | {es(r['eta'])} | {r['exact_orbit']:.5e} | {r['pred_orbit']:.5e} | {r['exact_orbit'] / r['pred_orbit']:.4f} | "
                         f"{r['exact_point']:.5e} | {r['pred_point']:.5e} | {r['exact_point'] / r['pred_point']:.4f} | {r['slow']:.5e} |")
    if Cres:
        lines += ["", f"## C. The real loop with an Ornstein–Uhlenbeck measurement noise (k_e = {KE:g}, B = {B_OU:g}, RK4 dt = {DT_C:g}, {T_C:g} time units)", "",
                  "E‖f̂ − f‖²_ρ at the records (mean ± standard error over the runs) against the linear prediction: the transient covariance of the exactly linearised "
                  "loop with the OU noise plus the second-order mean deviation (their sum, and the two parts).", ""]
        for std in STDS_C:
            for eta in (0.0, np.inf):
                rs = [r for r in Cres if r['kind'] == 'K0' and r['std'] == std and r['eta'] == eta]
                if not rs:
                    continue
                E = np.array([r['err'] for r in rs])
                pr = Cpred[(std, eta)]
                lines += [f"Noise std {std:g} (variance {std ** 2:.1e}), η = {es(eta)}, {len(rs)} runs:", "",
                          "| t | Monte Carlo | prediction | covariance part | mean part |", "|---|---|---|---|---|"]
                for k_, t in enumerate(pr['t']):
                    if k_ == 0:
                        continue
                    lines.append(f"| {t:g} | {E[:, k_].mean():.4e} ± {E[:, k_].std(ddof=1) / np.sqrt(len(rs)):.1e} | {pr['cov'][k_] + pr['mean'][k_]:.4e} | "
                                 f"{pr['cov'][k_]:.4e} | {pr['mean'][k_]:.4e} |")
                lines += [""]
        lift = [r for r in Cres if r['kind'] == 'lift']
        if lift:
            r = lift[0]
            twin = [x for x in Cres if x['kind'] == 'K0' and x['seed'] == r['seed'] and x['std'] == r['std'] and x['eta'] == r['eta']][0]
            lines.append(f"The lift m = 15 (counts 5, 5, 5) driven by the same noise path as the K₀-atom run (seed {r['seed']}, std {r['std']:g}, η = {es(r['eta'])}): "
                         f"the largest difference of f̂ on the grid over the records {np.abs(r['fh'].astype(float) - twin['fh'].astype(float)).max():.1e}, of E‖f̂ − f‖² "
                         f"{np.abs(r['err'] - twin['err']).max():.1e} — the noise is coherent: a replicated start stays replicated.")
    lines += ["", "Numbers are produced by `experiments/exp_measurement_noise.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'measurement_noise.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'measurement_noise_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            report(z['A'].item(), list(z['B']), list(z['C']), z['P'].item(), t0)
        return
    if '--only-A' in sys.argv:                                              # Part A recomputed, merged with the stored B, C and P
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        save['A'] = np.array(part_A(), dtype=object)
        save_raw(raw, **save)
        report(save['A'].item(), list(save['B']), list(save['C']), save['P'].item(), t0)
        return
    A = part_A()
    with mp.Pool(N_WORKERS) as pool:
        Bres = pool.map(job_B, [(ke, eta) for ke in KES_B for eta in ETAS_B], chunksize=1)
        save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bres, dtype=object))
        jobs = [('K0', eta, std, seed) for std in STDS_C for eta in (0.0, np.inf) for seed in range(N_MC)] + [('lift', np.inf, STDS_C[1], 0)]
        Cres = pool.map(job_C, jobs, chunksize=1)
    Cpred = {(std, eta): predictions_C(std, eta) for std in STDS_C for eta in (0.0, np.inf)}
    save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bres, dtype=object), C=np.array(Cres, dtype=object), P=np.array(Cpred, dtype=object))
    report(A, Bres, Cres, Cpred, t0)


if __name__ == '__main__':
    main()
