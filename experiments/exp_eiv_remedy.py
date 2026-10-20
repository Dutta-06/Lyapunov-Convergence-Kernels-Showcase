"""
Proposition D‴.31: the adaptive prior is a washout leak, and the desired regressor turns the noise's persistent force into a least-squares one
==============================================================================================================================================

The instance of Proposition D‴.30 (`exp_measurement_noise.py`: the three bumps, σ₀ = 0.15σ_crit, x_m = 2 sin t, the lift of the exact
representation of q = 1/3 with the prior at its atoms) and three versions of the controller and the adaptive law under the measurement
y = x + ν: the measured regressor (MR, Proposition D‴.30: ϑ̇ ∋ G∇_ϑ((y − x_m)Φ̃(y; ϑ)), the compensation f̂(y)), the desired regressor (DR: the
regressor at the reference, ∇_ϑΦ̃(x_m; ϑ), the compensation f̂(y)) and the desired compensation adaptation law (DR*: both at the reference).

A. The second-order force and its stationary responses (no dynamics), at k_e = 25, 50, 100, 200 and the Ornstein–Uhlenbeck rates B = 100,
   1000: the function-level stationary bias per unit (Eν²)² under the fixed prior, Fᵀ(H + σ₀I)⁻¹H(H + σ₀I)⁻¹F, and under the adaptive prior,
   FᵀH⁻¹F — the squared norm of the smallest target perturbation whose learning force is F — for the force of the noise's own second
   moments (the quasi-static `eiv_force`) and for the full second-order force (`second_order_mean`: the parameters' fluctuation
   included; the bias of the mean parameters and of the mean function); for DR and DR* the bound ‖g‖²_ρ, g = k_eE[e]; the rescaled
   projected mean tracking error k_e²‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁ under the adaptive prior and the mean tracking error itself; the noise level at
   which the bias equals the floor of the parameters.
B. The stationary identity RH + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a) (η = σ₀; κ/η = 0, 0.1, 1, 10, and the fixed
   prior): the leak's zero-frequency gain S₀, the weakest mode's response and bias factor, the errors-in-variables bias of MR and DR, the
   bias of an anchor misplaced by a unit along the weakest and along the strongest mode, and the identity's residual (from the block); and
   priors whose zero-frequency leak S₀ does not commute with H, R and B again from their dynamic blocks (`anchored_block`: a prior fixed on
   the amplitudes or on the centres and adaptive at η = σ₀ elsewhere, one active face of the projected adaptive prior at η = ∞, a matrix
   anchor gain K = ηP_amp, a random anchor gain with K + Kᵀ ⪰ 0), against (H + S₀)⁻¹ and (H + S₀)⁻¹S₀ with S₀ by elimination: RH + B = I
   holds and HR + B = I does not (round 22 of the review, ID-5: the identity derived from the block, not from the definitions).
C. The real loop (RK4, dt = 10⁻³, k_e = 25, B = 100, the noise's std 0.01, η = ∞, N_RUNS runs of T_C time units): MR, DR and DR* driven by
   the same noise paths (common random numbers). E‖f̂ − f‖²_ρ against the prediction — the transient covariance of the linearised loop
   (MR and DR share it: their linearisations coincide) plus the second-order mean deviation — and the run mean of f̂_MR − f̂_DR on the
   orbit's grid against its prediction Ĵ(ȳ_MR − ȳ_DR) (the first-order responses cancel path by path).
D. The separations under the desired regressor (Corollary D‴.33): the separations' second-order drift S_k at the lift for MR, DR and DR*
   (with the regressor at x_m only the mean error's term remains), the noise thresholds under the adaptive prior and the cohesion leak, the
   drift at DR's and DR*'s stationary points (k_eE[e] = (I − P)g) and along DR's second-order transient; the real loop — the lift m = 15
   with class 3 separated along MR's top eigenvector, N_D antithetic pairs per law and std (0.05, 0.1) on common paths — the separation's
   rate under MR and DR.

Output: results/chaos/eiv_remedy.md and eiv_remedy_raw.npz (ER_WORKERS, ER_RUNS, ER_TC; `--from-raw` regenerates the md; `--only-AB`
recomputes Parts A and B and merges them, `--only-B` Part B, `--only-D` Part D, `--only-Dpred` Part D's predictions; ER_ND the pairs of
Part D).
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
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN
import exp_measurement_noise as X                                   # the instance of Proposition D‴.30

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
VARIANTS = ('measured', 'reference', 'reference_both')
NAMES = {'measured': 'MR', 'reference': 'DR', 'reference_both': 'DR*'}
LAW = {'measured': ('measured', 'measured'), 'reference': ('reference', 'measured'), 'reference_both': ('reference', 'reference')}
KES_A = (25.0, 50.0, 100.0, 200.0)
BS_A = (100.0, 1000.0)
KAPPAS_B = (0.0, 0.1, 1.0, 10.0)                                   # κ/η at η = σ₀; then the fixed prior
KE, B_OU, DT_C, STD_C = X.KE, X.B_OU, X.DT_C, 0.01
N_RUNS = int(os.environ.get('ER_RUNS', 12))
T_C = float(os.environ.get('ER_TC', 300.0))
REC_C = 25.0
N_WORKERS = int(os.environ.get('ER_WORKERS', 3))
STDS_D = (0.05, 0.1)
N_D = int(os.environ.get('ER_ND', 4))                                 # antithetic pairs per (law, std)


def collapsed():
    net, atoms = X.setup()
    return net, atoms, AF.CollapsedVariational(net, X.xr, X.PER, X.G, X.SIGMA0, np.inf, atoms, X.Q3)


# ---------------------------------------------------------------------------- A
def job_A(args):
    """The full second-order mean (`second_order_mean`) of variant v at (k_e, B) under the adaptive (η = ∞) or the fixed prior (η = 0): its
    biases (`mean_biases`) per unit (Eν²)²."""
    ke, B, v, eta = args
    t0 = time.time()
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, ke, X.G, X.SIGMA0, np.arange(3), 3, eta, weights=X.Q3, regressor=reg, compensation=comp)
    out = MN.second_order_mean(cl, cv, X.fprime, B, 1.0)
    res = MN.mean_biases(H, out, X.SIGMA0 if eta == 0 else None)
    print(f"[{time.time() - t0:.0f}s] A {NAMES[v]} ke {ke:g} B {B:g} eta {eta:g}", flush=True)
    return dict(ke=ke, B=B, v=v, eta=eta, per=out['per_err'], **res)


def part_A(pool=None):
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    jobs = [(ke, B, v, eta) for ke in KES_A for B in BS_A for v in VARIANTS for eta in (np.inf, 0.0)]
    full = {j[:4]: r for j, r in zip(jobs, (pool.map(job_A, jobs, chunksize=1) if pool is not None else map(job_A, jobs)))}
    rows = []
    for ke in KES_A:
        for B in BS_A:
            for v in VARIANTS:
                F, g = MN.eiv_force(cv, X.fprime, X.fsecond, ke, B, 1.0, variant=v, return_g=True)
                Fs = ke * F                                           # slow time, per unit of the noise's variance
                Fl = ke * MN.eiv_force(cv, X.fprime, X.fsecond, ke, B, 1.0, variant='reference' if v == 'measured' else v)   # ⟨Ĵᵀg⟩
                F1 = Fs - Fl                                          # the regressor's errors in variables proper (MR only)
                adapt = MN.equivalent_perturbation(H, Fs)
                fixed = MN.response_size(H, Fs, X.SIGMA0)
                merr = MN.equivalent_perturbation(H, F1)              # ‖P(g − Ĵȳ_∞)‖²_ρ under the adaptive prior: F₁ᵀH⁻¹F₁
                floor = MN.floor_closed_form(H, Q, X.SIGMA0, np.inf, 2.0 / B, ke)[0]   # the parameters' floor per unit Eν² (N₀ = 2Eν²/B)
                fa, ff = full[(ke, B, v, np.inf)], full[(ke, B, v, 0.0)]
                rows.append(dict(ke=ke, B=B, v=v, adapt=adapt, fixed=fixed, g2=float(np.mean(g ** 2)), merr=merr, floor=floor,
                                 cross=floor / fa['full'], full_adapt=fa, full_fixed=ff))
    print("[A]", flush=True)
    return dict(rows=rows, h=np.linalg.eigvalsh(H))


# ---------------------------------------------------------------------------- B
def part_B():
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    h, V = np.linalg.eigh(H)
    F = {v: KE * MN.eiv_force(cv, X.fprime, X.fsecond, KE, B_OU, 1.0, variant=v) for v in ('measured', 'reference')}
    eta = X.SIGMA0
    cases = [(eta, k * eta) for k in KAPPAS_B] + [(0.0, 0.0)]
    rows = []
    for e_, kap in cases:
        R, Bm = MN.stationary_operators(H, X.SIGMA0, e_, kap)
        S0 = MN.zero_frequency_leak(X.SIGMA0, e_, kap)
        resid = float(np.abs(H @ R + Bm - np.eye(H.shape[0])).max())
        r1 = float(V[:, 0] @ R @ V[:, 0]); b1 = float(V[:, 0] @ Bm @ V[:, 0])
        eiv = {v: float((R @ F[v]) @ H @ (R @ F[v])) for v in F}
        anc = [float((Bm @ V[:, j]) @ H @ (Bm @ V[:, j])) for j in (0, len(h) - 1)]
        rows.append(dict(eta=e_, kappa=kap, S0=S0, resid=resid, resid_RH=float(np.abs(R @ H + Bm - np.eye(H.shape[0])).max()), r1=r1, b1=b1,
                         eiv=eiv, anchor=anc))
    # priors whose zero-frequency leak does not commute with H, R and B from their dynamic blocks (whitened coordinates (k, l) ↦ 3k + l:
    # l = 0 the amplitude, 1 the centre, 2 the width), against the definitions with S₀ by elimination (round 22, ID-5)
    d = H.shape[0]
    mask = lambda idx: np.isin(np.arange(d), idx)
    Mr = np.random.default_rng(7).standard_normal((d, d))
    Kr = eta * (Mr @ Mr.T / np.linalg.norm(Mr @ Mr.T, 2) + 0.5 * (Mr - Mr.T) / np.linalg.norm(Mr - Mr.T, 2))   # K + Kᵀ ⪰ 0, not symmetric
    cases = [('fixed on the amplitudes, η = σ₀ elsewhere', eta, None, mask([0, 3, 6])),
             ('fixed on the centres, η = σ₀ elsewhere', eta, None, mask([1, 4, 7])),
             ("one active face (atom 1's amplitude), η = ∞", np.inf, None, mask([0])),
             ('anchor gain K = ηP_amp', eta, eta * np.diag(mask([0, 3, 6]).astype(float)), None),
             ('a random anchor gain, K + Kᵀ ⪰ 0, ‖K‖ ≈ η', eta, Kr, None)]
    gen_rows = []
    I = np.eye(d)
    for name, e_, K, fixed in cases:
        R, Bm = MN.anchored_operators(H, X.SIGMA0, e_, K, fixed)
        S0m = MN.anchored_leak(H, X.SIGMA0, e_, K, fixed)
        R0, B0 = MN.zero_frequency_operators(H, S0m)
        per_mode = np.abs(h * np.einsum('ij,ik,kj->j', V, R, V) + np.einsum('ij,ik,kj->j', V, Bm, V) - 1).max()
        eiv = {v: float((R @ F[v]) @ H @ (R @ F[v])) for v in F}
        gen_rows.append(dict(name=name, RH=float(np.abs(R @ H + Bm - I).max()), HR=float(np.abs(H @ R + Bm - I).max()), per_mode=float(per_mode),
                             comm=float(np.abs(H @ S0m - S0m @ H).max()), eiv=eiv,
                             vs_def=float(max(np.abs(R - R0).max() / np.abs(R0).max(), np.abs(Bm - B0).max())),
                             dim=int(MN.anchored_block(H, X.SIGMA0, e_, K, fixed)[0].shape[0])))
    print("[B]", flush=True)
    return dict(rows=rows, h=h, general=gen_rows, n_weak=int(np.sum(h < X.SIGMA0)))


# ---------------------------------------------------------------------------- C
def job_C(args):
    v, seed = args
    t0 = time.time()
    net, atoms = X.setup()
    grid = 2 * np.sin(np.arange(X.N_GRID) * X.PER / X.N_GRID); Fg = X.f_vec(grid)
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, np.arange(3), 3, np.inf, weights=X.Q3,
                          regressor=reg, compensation=comp)
    z0 = cl.lift(atoms)

    def rec(t, z):
        Th = z[1:1 + 9].reshape(3, 3)
        fh = np.array([cl.f_hat(np.array([x]), Th)[0] for x in grid])
        return np.concatenate([[np.mean((fh - Fg) ** 2)], fh - Fg])
    out = np.array(MN.simulate_noisy(cl, z0, 3, T_C, DT_C, B_OU, STD_C ** 2, np.random.default_rng(5000 + seed), int(round(REC_C / DT_C)), rec))
    print(f"[{time.time() - t0:.0f}s] C {NAMES[v]} seed {seed}", flush=True)
    return dict(v=v, seed=seed, err=out[:, 0], dev=out[:, 1:])


def predictions_C():
    """The linear prediction at the records for each variant: the transient covariance part tr(HP_yy) (the linearised loop with the OU noise,
    MR and DR identical; DR* with the error channel −(k_e − a) and the noise input −k_e), the mean part ȳᵀHȳ, and the mean deviation Ĵȳ on
    the grid."""
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    var = STD_C ** 2; N0 = 2 * var / B_OU
    tev = np.arange(0.0, T_C + 1e-9, REC_C)
    tg = np.arange(X.N_GRID) * X.PER / X.N_GRID
    Jg = np.array([cv.J(s)[0] for s in tg])                              # the whitened Jacobian on the grid (x = 2 sin t_k)
    out = {}
    for v in VARIANTS:
        comp = LAW[v][1]
        gen = lambda t: cv.generator(t, KE)
        if comp == 'reference':
            gen = MN.reference_compensation(gen, cv, X.fprime, KE)
        D = cv.dim
        b0 = MN.collapsed_noise_input(cv, X.fprime, KE, compensation=comp)
        Ma, cfun = MN.ou_augmented(gen, b0, D, B_OU, N0)
        P0 = np.zeros((D + 1, D + 1)); P0[-1, -1] = var
        tt, Ps = MN.transient_covariance(Ma, cfun, D + 1, N0, T_C, tev, P0=P0)
        cov = np.array([np.trace(H @ P[1:1 + cv.d, 1:1 + cv.d]) for P in Ps])
        F = KE * MN.eiv_force(cv, X.fprime, X.fsecond, KE, B_OU, var, variant=v)
        ys = np.array([MN.mean_response(H, F, X.SIGMA0, np.inf, t / KE) for t in tev])
        out[v] = dict(t=tev, cov=cov, mean=np.einsum('ti,ij,tj->t', ys, H, ys), dev=ys @ Jg.T)
    return out


# ---------------------------------------------------------------------------- report
# ---------------------------------------------------------------------------- D (Corollary D‴.33)
def part_D_pred():
    """The separations' second-order drift S_k (real time, per unit Eν²; `separation_drift`) at the lift for MR, DR and DR*: its extreme
    eigenvalues per class and the noise thresholds — the std at which σ₀/k_e (the adaptive prior) or σ₀w_k/k_e (the cohesion leak of
    Proposition D‴.32, w_k = 1/3) equals the largest eigenvalue; for DR and DR* at the adaptive prior's stationary point (k_eE[e] = g − Ĵȳ,
    ȳ = H⁻¹F): the residual mean error and the drift; and along DR's second-order transient from the lift (`second_order_transient`, the
    regressor at x_m so no ∂_x term) the drift's largest eigenvalue per 50-unit window and its maximum over the samples."""
    t0 = time.time()
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    lift = {}
    u = None
    for v in VARIANTS:
        S = MN.separation_drift(cv, X.fprime, X.fsecond, KE, B_OU, 1.0, variant=v)
        ev = [np.linalg.eigvalsh(sk) for sk in S]
        u = np.linalg.eigh(S[2])[1][:, -1] if u is None else u                   # MR's top eigenvector of S₃: the real loop's start
        lift[v] = dict(lo=[float(e[0]) for e in ev], hi=[float(e[-1]) for e in ev], rq=float(u @ S[2] @ u))
    scan = {}                                                             # the drift's largest eigenvalue over the classes against k_e (B = 4k_e)
    for ke in KES_A:
        for v in VARIANTS:
            S = MN.separation_drift(cv, X.fprime, X.fsecond, ke, 4 * ke, 1.0, variant=v)
            scan[(ke, v)] = float(max(np.linalg.eigvalsh(sk)[-1] for sk in S))
    N = 512
    t = np.arange(N) * X.PER / N
    J = np.array([AF.whitened_jacobian(net, np.array([float(X.xr(s_)[0])]), cv.atoms, cv.w, cv.G)[0] for s_ in t])
    stat = {}
    for v in ('reference', 'reference_both'):
        F, g = MN.eiv_force(cv, X.fprime, X.fsecond, KE, B_OU, 1.0, N=N, variant=v, return_g=True)
        Ee = (g - J @ np.linalg.solve(H, KE * F)) / KE
        ev = [np.linalg.eigvalsh(sk) for sk in MN.separation_drift_series(cv, t, Ee, np.zeros(N)).mean(axis=0)]
        stat[v] = dict(rms_lift=float(np.sqrt(np.mean((g / KE) ** 2))), rms=float(np.sqrt(np.mean(Ee ** 2))),
                       lo=[float(e[0]) for e in ev], hi=[float(e[-1]) for e in ev])
    reg, comp = LAW['reference']
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, np.arange(3), 3, np.inf, weights=X.Q3, regressor=reg, compensation=comp)
    o = MN.second_order_transient(cl, cv, X.fprime, B_OU, 1.0, X.T_D, dt=0.05)
    Sd = MN.separation_drift_series(cv, o['t'], o['Ee'], np.zeros_like(o['Ee']))
    edges = np.arange(0.0, X.T_D + 1e-9, 50.0)
    win = [dict(a=a, b=b, lam=[float(np.linalg.eigvalsh(Sd[(o['t'] >= a) & (o['t'] < b)].mean(axis=0)[k])[-1]) for k in range(3)],
                rmsE=float(np.sqrt(np.mean(o['Ee'][(o['t'] >= a) & (o['t'] < b)] ** 2)))) for a, b in zip(edges[:-1], edges[1:])]
    pmax = [float(max(np.linalg.eigvalsh(np.mean(Sd[i:i + 126, k], axis=0))[-1] for i in range(0, len(o['t']) - 126, 21))) for k in range(3)]
    print(f"[{time.time() - t0:.0f}s] D pred", flush=True)
    return dict(lift=lift, stat=stat, win=win, period_max=pmax, w=X.Q3.tolist(), scan=scan)


def job_D(args):
    """The lift m = 15 with class 3's particles separated by c_i·10⁻⁶ along G^{1/2}u (u the top eigenvector of MR's S₃ at the lift, Part D of
    `exp_measurement_noise.py`), one noise path under the law v (std, seed, sign: ±1 the antithetic pair; the paths common to the laws): the
    separation's growth rate, least squares on log‖separation‖ over t ≥ 50."""
    v, std, seed, sign = args
    t0 = time.time()
    net, atoms, cv = collapsed()
    u = np.linalg.eigh(MN.separation_drift(cv, X.fprime, X.fsecond, KE, B_OU, 1.0)[2])[1][:, -1]
    labels = np.repeat(np.arange(3), 5); m = 15
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, labels, 3, np.inf, regressor=reg, compensation=comp)
    z0 = cl.lift(atoms)
    for i, c in zip(range(10, 15), (2.0, 1.0, 0.0, -1.0, -2.0)):
        z0[1 + 3 * i:4 + 3 * i] += 1e-6 * c * np.sqrt(X.G) * u

    def rec(t, z):
        Th = z[1:1 + m * 3].reshape(m, 3)[10:15]
        return np.array([t, np.sqrt(np.sum(((Th - Th.mean(axis=0)) / np.sqrt(X.G)) ** 2))])
    out = np.array(MN.simulate_noisy(cl, z0, m, X.T_D, DT_C, B_OU, std ** 2, np.random.default_rng(5000 + seed), int(round(5.0 / DT_C)), rec,
                                     sign=sign))
    msk = out[:, 0] >= 50.0
    rate = float(np.polyfit(out[msk, 0], np.log(out[msk, 1]), 1)[0])
    print(f"[{time.time() - t0:.0f}s] D {NAMES[v]} std {std:g} seed {seed} sign {sign:+g}", flush=True)
    return dict(v=v, std=std, seed=seed, sign=sign, rate=rate, t=out[:, 0], sep=out[:, 1])


def report(A, Bp, Cres, Cpred, t0=None, D=None):
    es = lambda e: '∞' if np.isinf(e) else ('0 (fixed)' if e == 0 else f'{e:.5f}')
    L = ["# Proposition D‴.31: the adaptive prior is a washout leak; the desired regressor (`exp_eiv_remedy.py`)", "",
         f"The instance of Proposition D‴.30 (σ₀ = {X.SIGMA0:.5f}, x_m = 2 sin t, the lift of the exact representation of q = 1/3 with the prior at its atoms) "
         "under the measurement y = x + ν, with an Ornstein–Uhlenbeck noise of rate B. MR: the measured regressor and compensation (Proposition D‴.30); DR: the "
         "regressor at the reference, the compensation f̂(y); DR*: both at the reference (the desired compensation adaptation law).", "",
         "## A. The second-order force and its stationary responses (no dynamics)", "",
         "Per unit (Eν²)², function level, the averaged stationary response: the bias of the mean parameters under the fixed prior and under the adaptive prior "
         "(η = ∞, the noise part's for every η > 0; = FᵀH⁻¹F, the squared norm of the smallest target perturbation exerting F) for the force of the noise's own second moments (quasi-static, "
         "`eiv_force`) and for the full second-order force (the parameters' fluctuation included, `second_order_mean`), the bias of the mean function under the "
         "adaptive prior with the full force (‖Ĵȳ + ½D²_Θf̂[Σ_ΘΘ]‖²_ρ), the least-squares bound ‖g‖²_ρ (g = k_eE[e] of the noise part; DR, DR*), the rescaled "
         "projected mean tracking error at the adaptive prior's stationary point, k_e²‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁ (F₁ the regressor's errors in variables; zero by "
         "construction for DR and DR*), and the mean tracking error itself ‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁/k_e². The last column is the noise variance at which the adaptive "
         "prior's full bias equals its floor of the parameters (N₀/2k_e)tr Q, N₀ = 2Eν²/B (std above 0.1 is outside the small-noise regime).", "",
         "| k_e | B | law | fixed: noise part | fixed: full | adaptive: noise part | adaptive: full | adaptive: mean function | ‖g‖²_ρ | k_e²‖P E[e]‖²_ρ | ‖P E[e]‖²_ρ | Eν² at bias = floor (std) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in A['rows']:
        g2 = f"{r['g2']:.4e}" if r['v'] != 'measured' else '—'
        fa, ff = r['full_adapt'], r['full_fixed']
        L.append(f"| {r['ke']:g} | {r['B']:g} | {NAMES[r['v']]} | {r['fixed']:.4e} | {ff['full']:.4e} | {r['adapt']:.4e} | {fa['full']:.4e} | "
                 f"{fa['function']:.4e} | {g2} | {r['merr']:.3e} | {r['merr'] / r['ke'] ** 2:.3e} | {r['cross']:.2e} ({np.sqrt(r['cross']):.4f}) |")
    rel = lambda r: r['full_adapt']['full'] / r['adapt'] - 1
    L += ["", "The parameters' fluctuation's share of the adaptive prior's bias (full over noise part − 1): " + "; ".join(
        f"{NAMES[v]} " + ", ".join(f"{rel(r):+.2%}" for r in A['rows'] if r['v'] == v) for v in VARIANTS) + " (the rows' order); of the mean function's: "
        + "; ".join(f"{NAMES[v]} " + ", ".join(f"{r['full_adapt']['function'] / r['adapt'] - 1:+.2%}" for r in A['rows'] if r['v'] == v) for v in VARIANTS) + "."]
    if 'param' in A['rows'][0]['full_adapt']:
        L += ["", "The parameters' part's own function-level response under the adaptive prior (and its learning-law part's), per unit (Eν²)²: " + "; ".join(
            f"{NAMES[v]} " + ", ".join(f"{r['full_adapt']['param']:.2e} ({r['full_adapt']['law']:.2e})" for r in A['rows'] if r['v'] == v) for v in VARIANTS)
            + " (the rows' order) — set by the parameters' covariance, falling with k_e."]
    mr = {(r['ke'], r['B']): r['full_adapt']['full'] for r in A['rows'] if r['v'] == 'measured'}
    dr = {(r['ke'], r['B']): r['full_adapt']['full'] for r in A['rows'] if r['v'] == 'reference'}
    mr25, mr200 = mr[(25.0, 1000.0)], mr[(200.0, 1000.0)]
    law = lambda ke, B: (ke * B / (ke + B)) ** 2
    L += ["", "The adaptive prior's full bias, MR over DR: " + ", ".join(f"{mr[k] / dr[k]:.3g} (k_e = {k[0]:g}, B = {k[1]:g})" for k in sorted(mr)) + ". "
          f"MR's grows from k_e = 25 to 200 at B = 1000 by {mr200 / mr25:.1f} — the law (k_eB/(k_e + B))² gives {law(200.0, 1000.0) / law(25.0, 1000.0):.1f}, k_e² 64; "
          "the noise parts of DR's and DR*'s stay below ¼‖∂²_xf‖²_ρ = "
          f"{0.25 * np.mean(X.fsecond(2 * np.sin(np.arange(4096) * X.PER / 4096)) ** 2):.4e} for every k_e and B, and their full biases at most "
          f"{max(r['full_adapt']['full'] for r in A['rows'] if r['v'] != 'measured'):.4f} here."]
    L += ["", "## B. The stationary identity RH + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a)", "",
          f"η = σ₀; R and B from the averaged collective block at zero frequency (k_e = {KE:g}, B = {B_OU:g}). h₁ = {Bp['h'][0]:.4e}, h_max = {Bp['h'][-1]:.4e}. "
          "The errors-in-variables biases are function-level per unit (Eν²)²; the anchor's biases are function-level for a unit misplacement along the weakest and "
          "the strongest mode.", "",
          "| prior | S₀ (zero-frequency leak) | weakest mode: response r₁ | bias factor b₁ | h₁r₁ + b₁ | EIV bias, MR | EIV bias, DR | anchor bias, weakest | anchor bias, strongest | max|RH + B − I| |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in Bp['rows']:
        name = 'fixed (η = 0)' if r['eta'] == 0 else (f"adaptive, κ = 0" if r['kappa'] == 0 else f"κ = {r['kappa'] / r['eta']:g}η")
        L.append(f"| {name} | {r['S0']:.4e} | {r['r1']:.4e} | {r['b1']:.4f} | {Bp['h'][0] * r['r1'] + r['b1']:.6f} | {r['eiv']['measured']:.4e} | "
                 f"{r['eiv']['reference']:.4e} | {r['anchor'][0]:.3e} | {r['anchor'][1]:.3e} | {r.get('resid_RH', r['resid']):.1e} |")
    ratios = [r['eiv']['measured'] / r['eiv']['reference'] for r in Bp['rows']]
    L += ["", "Along the family the zero-frequency leak S₀ = σ₀κ/(η + κ) trades the anchor's bias for the response to persistent forces, mode by mode, at "
          f"h_jr_j + b_j = 1; the trade binds on the modes with h_j < σ₀ ({Bp.get('n_weak', '—')} of {len(Bp['h'])} here), where 1/h_j exceeds the fixed "
          "prior's 1/σ₀. The desired regressor changes the force instead: MR's bias over DR's (noise parts) along the family "
          + ", ".join(f"{x:.3g}" for x in ratios) + "."]
    if 'general' in Bp:
        if 'vs_def' in Bp['general'][0]:
            L += ["", "Priors whose zero-frequency leak S₀ does not commute with H: R and B from the dynamic block of (y, z_f) at zero frequency "
                  "(`anchored_block`; the block's dimension), against (H + S₀)⁻¹ and (H + S₀)⁻¹S₀ with S₀ = σ₀[E_xE_xᵀ + E_f(ηI + K)⁻¹KE_fᵀ] by "
                  "elimination of z_f (the largest relative difference of R, the largest difference of B): RH + B = I holds, HR + B = I does not; "
                  "per mode h_jr_jj + b_jj = 1 holds whatever S₀.", "",
                  "| prior | block dimension | block against the definitions | max|HS₀ − S₀H| | max|RH + B − I| | max|HR + B − I| | max_j|h_jr_jj + b_jj − 1| | EIV bias, MR | EIV bias, DR |",
                  "|---|---|---|---|---|---|---|---|---|"]
            for r in Bp['general']:
                L.append(f"| {r['name']} | {r['dim']} | {r['vs_def']:.1e} | {r['comm']:.2e} | {r['RH']:.1e} | {r['HR']:.3f} | {r['per_mode']:.1e} | "
                         f"{r['eiv']['measured']:.4e} | {r['eiv']['reference']:.4e} |")
        else:
            L += ["", "Zero-frequency leaks S₀ that do not commute with H (R = (H + S₀)⁻¹, B = (H + S₀)⁻¹S₀): RH + B = I holds, HR + B = I does not; per mode "
                  "h_jr_jj + b_jj = 1 holds whatever S₀.", "",
                  "| S₀ | max|HS₀ − S₀H| | max|RH + B − I| | max|HR + B − I| | max_j|h_jr_jj + b_jj − 1| | EIV bias, MR | EIV bias, DR |", "|---|---|---|---|---|---|---|"]
            for r in Bp['general']:
                L.append(f"| {r['name']} | {r['comm']:.2e} | {r['RH']:.1e} | {r['HR']:.3f} | {r['per_mode']:.1e} | {r['eiv']['measured']:.4e} | {r['eiv']['reference']:.4e} |")
    if Cres:
        L += ["", f"## C. The real loop: MR, DR and DR* driven by the same noise paths (k_e = {KE:g}, B = {B_OU:g}, std {STD_C:g}, η = ∞, RK4 dt = {DT_C:g})", "",
              f"E‖f̂ − f‖²_ρ at the records (mean ± standard error over {N_RUNS} runs) against the prediction (covariance part + mean part); MR and DR share the "
              "covariance part (their linearisations coincide).", "",
              "| t | MR: Monte Carlo | MR: predicted | DR: Monte Carlo | DR: predicted | DR*: Monte Carlo | DR*: predicted | covariance part (MR, DR) | mean part: MR | DR | DR* |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        E = {v: np.array([r['err'] for r in sorted([r for r in Cres if r['v'] == v], key=lambda r: r['seed'])]) for v in VARIANTS}
        Dv = {v: np.array([r['dev'] for r in sorted([r for r in Cres if r['v'] == v], key=lambda r: r['seed'])]) for v in VARIANTS}
        for k_, t in enumerate(Cpred['measured']['t']):
            if k_ == 0:
                continue
            cells = []
            for v in VARIANTS:
                p = Cpred[v]
                cells += [f"{E[v][:, k_].mean():.4e} ± {E[v][:, k_].std(ddof=1) / np.sqrt(len(E[v])):.1e}", f"{p['cov'][k_] + p['mean'][k_]:.4e}"]
            L.append(f"| {t:g} | " + " | ".join(cells) + f" | {Cpred['measured']['cov'][k_]:.4e} | {Cpred['measured']['mean'][k_]:.4e} | "
                     f"{Cpred['reference']['mean'][k_]:.4e} | {Cpred['reference_both']['mean'][k_]:.4e} |")
        diff = Dv['measured'] - Dv['reference']                           # path by path: the first-order responses cancel
        pd = Cpred['measured']['dev'] - Cpred['reference']['dev']
        L += ["", "The run mean of f̂_MR − f̂_DR on the orbit's grid (common random numbers) against Ĵ(ȳ_MR − ȳ_DR): its mean square (with the run-to-run "
              "standard error of the mean function, in the same units), the prediction's, and the relative L²_ρ error of the prediction.", "",
              "| t | ⟨(run mean)²⟩ | predicted | relative error ‖run mean − Ĵ(ȳ_MR − ȳ_DR)‖/‖Ĵ(ȳ_MR − ȳ_DR)‖ | sampling error of the run mean (rms) | relative error over the relative sampling error |",
              "|---|---|---|---|---|---|"]
        for k_, t in enumerate(Cpred['measured']['t']):
            if k_ == 0:
                continue
            m_ = diff[:, k_].mean(axis=0); se = np.sqrt(np.mean(diff[:, k_].var(axis=0, ddof=1) / len(diff)))
            re_ = np.sqrt(np.mean((m_ - pd[k_]) ** 2) / np.mean(pd[k_] ** 2))
            L.append(f"| {t:g} | {np.mean(m_ ** 2):.4e} | {np.mean(pd[k_] ** 2):.4e} | {re_:.4f} | {se:.2e} | {re_ / (se / np.sqrt(np.mean(pd[k_] ** 2))):.2f} |")
    if D:
        pr = D['pred']
        thr = lambda lam, stiff: np.sqrt(stiff / KE / lam) if lam > 0 else np.inf
        fs = lambda x: '—' if np.isinf(x) else f'{x:.3f}'
        L += ["", "## D. The separations under the desired regressor (Corollary D‴.33)", "",
              f"The separations' second-order drift S_k at the lift (real time, per unit Eν²; k_e = {KE:g}, B = {B_OU:g}): its extreme eigenvalues per class and "
              "the noise std at which its largest eigenvalue reaches the separations' stiffness — σ₀/k_e under the adaptive prior, σ₀w_k/k_e under the cohesion "
              "leak (Proposition D‴.32, w_k = 1/3); '—' where every eigenvalue is negative.", "",
              "| law | class 1: [min, max] | class 2 | class 3 | threshold std, adaptive prior | threshold std, cohesion leak |", "|---|---|---|---|---|---|"]
        for v in VARIANTS:
            r = pr['lift'][v]; lam = max(r['hi'])
            L.append(f"| {NAMES[v]} | " + " | ".join(f"[{lo:+.5f}, {hi:+.5f}]" for lo, hi in zip(r['lo'], r['hi'])) + f" | {fs(thr(lam, X.SIGMA0))} | "
                     f"{fs(thr(lam, X.SIGMA0 * pr['w'][0]))} |")
        if 'scan' in pr:
            L += ["", "The drift's largest eigenvalue over the classes against k_e at B = 4k_e (the noise white over the loop's bandwidth), and k_e times it: "
                  "the measured regressor's drift is O(1), so its threshold std falls like k_e^{−1/2}; with the regressor at the reference the drift is "
                  "O(1/k_e), as the stiffness σ₀/k_e is, so the threshold does not fall with k_e.", "",
                  "| k_e | MR: λ_max | DR: λ_max; k_e·λ_max | DR*: λ_max; k_e·λ_max | threshold std (adaptive prior): MR; DR; DR* |", "|---|---|---|---|---|"]
            for ke in KES_A:
                sc = {v: pr['scan'][(ke, v)] for v in VARIANTS}
                L.append(f"| {ke:g} | {sc['measured']:+.4e} | {sc['reference']:+.3e}; {ke * sc['reference']:+.4f} | {sc['reference_both']:+.3e}; "
                         f"{ke * sc['reference_both']:+.4f} | " + "; ".join(fs(np.sqrt(X.SIGMA0 / ke / sc[v]) if sc[v] > 0 else np.inf) for v in VARIANTS) + " |")
        L += ["", "With the regressor at the reference the drift is the mean error's term ⟨E[e]G^{1/2}∇²_ϑΦ̃G^{1/2}⟩ alone, and k_eE[e] = g − Ĵȳ moves from g at "
              "the lift to (I − P)g at the adaptive prior's stationary point (ȳ = H⁻¹F): the rms of E[e] at the lift and at the stationary point, and the "
              "drift there.", "",
              "| law | rms E[e]: lift | stationary point | drift at the stationary point: class 1, 2, 3 [min, max] |", "|---|---|---|---|"]
        for v in ('reference', 'reference_both'):
            r = pr['stat'][v]
            L.append(f"| {NAMES[v]} | {r['rms_lift']:.4e} | {r['rms']:.4e} | " + "; ".join(f"[{lo:+.1e}, {hi:+.1e}]" for lo, hi in zip(r['lo'], r['hi'])) + " |")
        L += ["", f"Along DR's second-order transient from the lift (Eν² = 1): per 50-unit window the mean error's rms and the largest eigenvalue of the "
              "window-averaged S_k (classes 1, 2, 3); the largest over every one-period average: " + ", ".join(f"{x:+.2e}" for x in pr['period_max']) + ".", "",
              "| window | rms E[e] | λ_max(S₁), λ_max(S₂), λ_max(S₃) |", "|---|---|---|"]
        for w_ in pr['win']:
            L.append(f"| [{w_['a']:g}, {w_['b']:g}) | {w_['rmsE']:.4f} | " + ", ".join(f"{x:+.2e}" for x in w_['lam']) + " |")
        L += ["", f"The real loop: the lift m = 15 with class 3 separated by 10⁻⁶·(2, 1, 0, −1, −2) along MR's top eigenvector of S₃, {N_D} antithetic pairs per "
              f"(law, std) on paths common to the laws, {X.T_D:g} time units, the separation's rate (least squares over t ≥ 50) against −σ₀/k_e = "
              f"{-X.SIGMA0 / KE:.3e} plus Eν² times the law's S₃ at the lift — along the start's direction u (uᵀS₃u) and at its largest eigenvalue "
              "(the rate the separation tends to once it has turned onto the top eigenvector).", "",
              "| law | std | measured rate (pair means ± s.e.) | −σ₀/k_e + Eν²uᵀS₃u | −σ₀/k_e + Eν²λ_max(S₃) | runs growing |", "|---|---|---|---|---|---|"]
        for v in ('measured', 'reference'):
            for std in STDS_D:
                rs = [r for r in D['runs'] if r['v'] == v and r['std'] == std]
                if not rs:
                    continue
                seeds = sorted(set(r['seed'] for r in rs))
                pm = np.array([np.mean([r['rate'] for r in rs if r['seed'] == sd]) for sd in seeds])
                rq = pr['lift'][v].get('rq', np.nan)
                L.append(f"| {NAMES[v]} | {std:g} | {pm.mean():+.4e} ± {pm.std(ddof=1) / np.sqrt(len(pm)):.1e} | {-X.SIGMA0 / KE + std ** 2 * rq:+.4e} | "
                         f"{-X.SIGMA0 / KE + std ** 2 * pr['lift'][v]['hi'][2]:+.4e} | {sum(r['rate'] > 0 for r in rs)} of {len(rs)} |")
    L += ["", "Numbers are produced by `experiments/exp_eiv_remedy.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'eiv_remedy.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'eiv_remedy_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            report(z['A'].item(), z['B'].item(), list(z['C']), z['P'].item(), t0, D=z['D'].item() if 'D' in z.files else None)
        return
    if '--only-Dpred' in sys.argv:                                           # Part D's predictions recomputed, its runs kept
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        Dd = save['D'].item()
        Dd['pred'] = part_D_pred()
        save['D'] = np.array(Dd, dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=Dd)
        return
    if '--only-D' in sys.argv:                                               # Part D (Corollary D‴.33) computed, merged with the stored parts
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(N_WORKERS) as pool:
            pr = pool.apply_async(part_D_pred)
            runs = pool.map(job_D, [(v, std, seed, sign) for v in ('measured', 'reference') for std in STDS_D for seed in range(N_D)
                                    for sign in (1.0, -1.0)], chunksize=1)
            save['D'] = np.array(dict(pred=pr.get(), runs=runs), dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=save['D'].item())
        return
    if '--only-B' in sys.argv:                                               # Part B recomputed, merged with the stored A, C and P
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        save['B'] = np.array(part_B(), dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=save['D'].item() if 'D' in save else None)
        return
    if '--only-AB' in sys.argv:                                              # Parts A and B recomputed, merged with the stored C and P
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(N_WORKERS) as pool:
            save['A'] = np.array(part_A(pool), dtype=object)
        save['B'] = np.array(part_B(), dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=save['D'].item() if 'D' in save else None)
        return
    with mp.Pool(N_WORKERS) as pool:
        A, Bp = part_A(pool), part_B()
        X.save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bp, dtype=object))
        Cres = pool.map(job_C, [(v, seed) for seed in range(N_RUNS) for v in VARIANTS], chunksize=1)
        pr = pool.apply_async(part_D_pred)
        Druns = pool.map(job_D, [(v, std, seed, sign) for v in ('measured', 'reference') for std in STDS_D for seed in range(N_D)
                                 for sign in (1.0, -1.0)], chunksize=1)
        Dd = dict(pred=pr.get(), runs=Druns)
    Cpred = predictions_C()
    X.save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bp, dtype=object), C=np.array(Cres, dtype=object), P=np.array(Cpred, dtype=object),
               D=np.array(Dd, dtype=object))
    report(A, Bp, Cres, Cpred, t0, D=Dd)


if __name__ == '__main__':
    main()
