"""
Proposition D‴.31: the adaptive prior is a washout leak, and the desired regressor turns the noise's persistent force into a least-squares one
==============================================================================================================================================

The instance of Proposition D‴.30 (`exp_measurement_noise.py`: the three bumps, σ₀ = 0.15σ_crit, x_m = 2 sin t, the lift of the exact
representation of q = 1/3 with the prior at its atoms) and three versions of the controller and the adaptive law under the measurement
y = x + ν: the measured regressor (MR, Proposition D‴.30: ϑ̇ ∋ G∇_ϑ((y − x_m)Φ̃(y; ϑ)), the compensation f̂(y)), the desired regressor (DR: the
regressor at the reference, ∇_ϑΦ̃(x_m; ϑ), the compensation f̂(y)) and the desired compensation adaptation law (DR*: both at the reference).

A. The second-order force and its stationary responses (no dynamics), at k_e = 25, 50, 100, 200 and the Ornstein–Uhlenbeck rates B = 100,
   1000: the function-level stationary bias per unit (Eν²)² under the fixed prior, Fᵀ(H + σ₀I)⁻¹H(H + σ₀I)⁻¹F, and under the adaptive prior,
   FᵀH⁻¹F — the squared norm of the smallest target perturbation whose learning force is F; for DR and DR* the bound ‖g‖²_ρ, g = k_eE[e];
   the projected stationary mean tracking error under the adaptive prior; the noise level at which the bias equals the floor of the
   parameters.
B. The stationary identity HR + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a) (η = σ₀; κ/η = 0, 0.1, 1, 10, and the fixed
   prior): the leak's zero-frequency gain S₀, the weakest mode's response and bias factor, the errors-in-variables bias of MR and DR, the
   bias of an anchor misplaced by a unit along the weakest and along the strongest mode, and the identity's residual (from the block).
C. The real loop (RK4, dt = 10⁻³, k_e = 25, B = 100, the noise's std 0.01, η = ∞, N_RUNS runs of T_C time units): MR, DR and DR* driven by
   the same noise paths (common random numbers). E‖f̂ − f‖²_ρ against the prediction — the transient covariance of the linearised loop
   (MR and DR share it: their linearisations coincide) plus the second-order mean deviation — and the run mean of f̂_MR − f̂_DR on the
   orbit's grid against its prediction Ĵ(ȳ_MR − ȳ_DR) (the first-order responses cancel path by path).

Output: results/chaos/eiv_remedy.md and eiv_remedy_raw.npz (ER_WORKERS, ER_RUNS, ER_TC; `--from-raw` regenerates the md).
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


def collapsed():
    net, atoms = X.setup()
    return net, atoms, AF.CollapsedVariational(net, X.xr, X.PER, X.G, X.SIGMA0, np.inf, atoms, X.Q3)


# ---------------------------------------------------------------------------- A
def part_A():
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
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
                rows.append(dict(ke=ke, B=B, v=v, adapt=adapt, fixed=fixed, g2=float(np.mean(g ** 2)), merr=merr, floor=floor,
                                 cross=floor / adapt))
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
        rows.append(dict(eta=e_, kappa=kap, S0=S0, resid=resid, r1=r1, b1=b1, eiv=eiv, anchor=anc))
    print("[B]", flush=True)
    return dict(rows=rows, h=h)


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
def report(A, Bp, Cres, Cpred, t0=None):
    es = lambda e: '∞' if np.isinf(e) else ('0 (fixed)' if e == 0 else f'{e:.5f}')
    L = ["# Proposition D‴.31: the adaptive prior is a washout leak; the desired regressor (`exp_eiv_remedy.py`)", "",
         f"The instance of Proposition D‴.30 (σ₀ = {X.SIGMA0:.5f}, x_m = 2 sin t, the lift of the exact representation of q = 1/3 with the prior at its atoms) "
         "under the measurement y = x + ν, with an Ornstein–Uhlenbeck noise of rate B. MR: the measured regressor and compensation (Proposition D‴.30); DR: the "
         "regressor at the reference, the compensation f̂(y); DR*: both at the reference (the desired compensation adaptation law).", "",
         "## A. The second-order force and its stationary responses (no dynamics)", "",
         "Per unit (Eν²)², function level: the stationary bias under the fixed prior and under the adaptive prior (every η > 0; = FᵀH⁻¹F, the squared norm of the "
         "smallest target perturbation exerting F), the least-squares bound ‖g‖²_ρ (g = k_eE[e]; DR, DR*), and the projected mean tracking error at the adaptive "
         "prior's stationary point, ‖P(g − Ĵȳ_∞)‖²_ρ = F₁ᵀH⁻¹F₁ (F₁ the regressor's errors in variables, zero for DR and DR*). The last column is the noise "
         "variance at which the adaptive prior's bias equals its floor of the parameters (N₀/2k_e)tr Q, N₀ = 2Eν²/B; the bias dominates above it.", "",
         "| k_e | B | law | fixed prior | adaptive prior | ‖g‖²_ρ | projected mean error (adaptive) | Eν² at bias = floor (std) |", "|---|---|---|---|---|---|---|---|"]
    for r in A['rows']:
        g2 = f"{r['g2']:.4e}" if r['v'] != 'measured' else '—'
        L.append(f"| {r['ke']:g} | {r['B']:g} | {NAMES[r['v']]} | {r['fixed']:.4e} | {r['adapt']:.4e} | {g2} | {r['merr']:.3e} | "
                 f"{r['cross']:.2e} ({np.sqrt(r['cross']):.4f}) |")
    mr = {(r['ke'], r['B']): r['adapt'] for r in A['rows'] if r['v'] == 'measured'}
    dr = {(r['ke'], r['B']): r['adapt'] for r in A['rows'] if r['v'] == 'reference'}
    L += ["", "The adaptive prior's bias, MR over DR: " + ", ".join(f"{mr[k] / dr[k]:.3g} (k_e = {k[0]:g}, B = {k[1]:g})" for k in sorted(mr)) + ". "
          "MR's grows like k_e² where B ≫ k_e; DR's and DR*'s stay below ¼‖∂²_xf‖²_ρ(1 + O(1/k_e)) = "
          f"{0.25 * np.mean(X.fsecond(2 * np.sin(np.arange(4096) * X.PER / 4096)) ** 2):.4e} for every k_e and B."]
    L += ["", "## B. The stationary identity HR + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a)", "",
          f"η = σ₀; R and B from the averaged collective block at zero frequency (k_e = {KE:g}, B = {B_OU:g}). h₁ = {Bp['h'][0]:.4e}, h_max = {Bp['h'][-1]:.4e}. "
          "The errors-in-variables biases are function-level per unit (Eν²)²; the anchor's biases are function-level for a unit misplacement along the weakest and "
          "the strongest mode.", "",
          "| prior | S₀ (zero-frequency leak) | weakest mode: response r₁ | bias factor b₁ | h₁r₁ + b₁ | EIV bias, MR | EIV bias, DR | anchor bias, weakest | anchor bias, strongest | max|HR + B − I| |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in Bp['rows']:
        name = 'fixed (η = 0)' if r['eta'] == 0 else (f"adaptive, κ = 0" if r['kappa'] == 0 else f"κ = {r['kappa'] / r['eta']:g}η")
        L.append(f"| {name} | {r['S0']:.4e} | {r['r1']:.4e} | {r['b1']:.4f} | {Bp['h'][0] * r['r1'] + r['b1']:.6f} | {r['eiv']['measured']:.4e} | "
                 f"{r['eiv']['reference']:.4e} | {r['anchor'][0]:.3e} | {r['anchor'][1]:.3e} | {r['resid']:.1e} |")
    L += ["", "Along the family the zero-frequency leak S₀ = σ₀κ/(η + κ) trades the anchor's bias for the response to persistent forces, mode by mode, at "
          "h_jr_j + b_j = 1: no prior removes the bias and bounds the response at once. The desired regressor changes the force instead: its bias is three to "
          "four orders of magnitude below MR's under every prior of the family."]
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
              "| t | ⟨(run mean)²⟩ | predicted | relative error ‖run mean − Ĵ(ȳ_MR − ȳ_DR)‖/‖Ĵ(ȳ_MR − ȳ_DR)‖ | sampling error of the run mean (rms) |", "|---|---|---|---|---|"]
        for k_, t in enumerate(Cpred['measured']['t']):
            if k_ == 0:
                continue
            m_ = diff[:, k_].mean(axis=0); se = np.sqrt(np.mean(diff[:, k_].var(axis=0, ddof=1) / len(diff)))
            L.append(f"| {t:g} | {np.mean(m_ ** 2):.4e} | {np.mean(pd[k_] ** 2):.4e} | {np.sqrt(np.mean((m_ - pd[k_]) ** 2) / np.mean(pd[k_] ** 2)):.4f} | {se:.2e} |")
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
            report(z['A'].item(), z['B'].item(), list(z['C']), z['P'].item(), t0)
        return
    A, Bp = part_A(), part_B()
    X.save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bp, dtype=object))
    with mp.Pool(N_WORKERS) as pool:
        Cres = pool.map(job_C, [(v, seed) for seed in range(N_RUNS) for v in VARIANTS], chunksize=1)
    Cpred = predictions_C()
    X.save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bp, dtype=object), C=np.array(Cres, dtype=object), P=np.array(Cpred, dtype=object))
    report(A, Bp, Cres, Cpred, t0)


if __name__ == '__main__':
    main()
