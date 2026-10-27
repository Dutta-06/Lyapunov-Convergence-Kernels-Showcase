"""
Proposition D‴.50: the rich transit under label noise — the dual is a Langevin diffusion
====================================================================================================================================

The free mass–amplitude law on the 6 × 14 Gaussian designs of Proposition D‴.48 (`exp_mirror_transit`: 3-sparse truths, masses εb_i, b_i
log-uniform in [1/4, 4]) with an isotropic white label noise of intensity σ: dζ = (y − Φw)dt/N + (σ/N)dW (`dual_langevin`).

A. The Gibbs law on a small well-conditioned dictionary (4 × 8, orthonormal rows scaled by √8, masses 2b_i with b_i ∈ U(1, 2) — the
   lazy side, where every particle moves and the dual mixes in a few time units — and σ = 1, far from Gaussian): the mirror update (the
   dual's Euler–Maruyama, at two steps) against an independent random-walk Metropolis sampler of π ∝ exp(−F/T) — the coefficients' means
   and standard deviations, E[Φw] − y; the paradigm's law on its curves read in Stratonovich's sense by the stochastic Heun scheme in α
   follows the mirror update pathwise (the same noise; the strong error ∝ h) but drifts along the dual's conserved directions
   P_{ker Φ}η over long times (a discretisation that does not keep them samples another law); read in Itô's sense (Euler on
   (α, log β), the paradigm's update) the conserved quantities drift at the predicted rate.
B. The fluctuation–dissipation relation on the support (ε = 10⁻⁶, three designs, σ = 0.03, 0.1, 0.3): from the noiseless transit's end,
   the coefficients' covariance over a window against T·dw/dη (T = σ²/(2N)), their correlations against the design's, their mean against
   weighted basis pursuit, the inactive coefficients.
C. The noisy transit (the first design, ε = 10⁻⁴, 10⁻⁶, 10⁻⁸, σ = 0.1, 0.3, 1, four noise paths): the noise's share sup_τ|u − u_det|
   (u_det the noiseless run on the same grid) against κ = σ/(NΛ^{1/2}) and against its first-order response along the inverse scale space
   path (δν ← P_T(δν + κdW), the same noise path: for entries only κP_{T(τ)}W(τ)) — its error over all τ, outside windows of ±0.05 and
   ±0.1 around the path's events (where an entry's O(κ) shift makes the response O(1)-nonlinear), and in L²; the entries' order and times.
D. The stationary mean in the rich limit at larger noise (the first design, ε = 10⁻⁴, 10⁻⁶, 10⁻⁸, σ = 0.3, 1): the support's mean against
   weighted basis pursuit with its standard error, the inactive coefficients.

Output: results/chaos/transit_noise.md and transit_noise_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D`).
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
from src.core import mirror_free as MF
from src.core import rich_transit as RT
from src.core import dual_langevin as DL
from src.core.runcache import cached
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('TN_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TN_WORKERS', '2'))
NA_A, DA_A, EPS_A, SIG_A = 4, 8, 2.0, 1.0                                   # A: the small dictionary, masses 2b_i, the noise
H_A, T_A, N_RWM, BURN_A = (0.01, 0.005), 2.0e4, 1_000_000, 0.1               # A: mirror-update steps, horizon; Metropolis steps; burn-in
H_PATH_A, T_PATH_A, T_DRIFT_A = (0.01, 0.0025), 50.0, 2000.0                  # A: Heun against the mirror update, pathwise and long
SEEDS_B, EPS_B, SIGS_B, H_B, TW_B = (3, 4, 5), 1e-6, (0.03, 0.1, 0.3), 0.1, 4.0e4
SEED_CD = 3                                                              # C, D: the first design of exp_mirror_transit
EPS_C, SIGS_C, PATHS_C, H_C, REC_C = (1e-4, 1e-6, 1e-8), (0.1, 0.3, 1.0), 4, 0.25, 20
WIN_C = (0.05, 0.1)                                                       # C: windows around the path's events (in τ)
EPS_D, SIGS_D, H_D, TW_D = (1e-4, 1e-6, 1e-8), (0.3, 1.0), 0.1, 4.0e4
NB = 20                                                                  # batch means


def setup(seed, eps):
    L = X.law()
    Phi, ws, b = X.design(seed)
    y = Phi @ ws
    rho = b ** -0.5
    path = MF.iss_path(Phi, y, rho)
    q = L.Q(np.zeros(b.size), eps * b)
    return L, Phi, y, b, rho, path, q


def batch_se(x, nb=NB):
    """Standard error of the mean of a time series by batch means."""
    x = np.asarray(x)
    n = (x.shape[0] // nb) * nb
    bm = x[:n].reshape(nb, n // nb, *x.shape[1:]).mean(1)
    return bm.std(0, ddof=1) / np.sqrt(nb)


def plateau(L, Phi, y, q, path, eps):
    """The noiseless transit's end (1.5τ*), as the start of a noisy window."""
    F = DL.FastFree(L, q)
    Lam = RT.transit_scale(L, eps)
    ts, Z, W, al = DL.mirror_update(F, Phi, y, 0.0, 0.25, 1.5 * path['taus'][-1] * Lam, rec_every=10 ** 9)
    return F, Z[-1], al, W[-1]


# --- A ---------------------------------------------------------------------------------------------------------------------------------
def small_dictionary():
    """4 × 8, orthonormal rows scaled by √8 (ΦΦᵀ = 8I), y = Φ(0.6, −0.4, 0, …), masses 2b_i, b_i ∈ U(1, 2)."""
    L = X.law()
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((DA_A, DA_A)))
    Phi = np.sqrt(DA_A) * Q[:NA_A]
    y = Phi @ np.r_[0.6, -0.4, np.zeros(DA_A - 2)]
    b = rng.uniform(1.0, 2.0, DA_A)
    return L, Phi, y, b, L.Q(np.zeros(DA_A), EPS_A * b)


@cached
def a_job(args):
    kind, h, seed = args
    warnings.filterwarnings('ignore')
    L, Phi, y, b, q = small_dictionary()
    N, d = Phi.shape
    T = SIG_A ** 2 / (2 * N)
    F = DL.FastFree(L, q)
    t0 = time.time()
    if kind == 'rwm':
        ts, Z, W, al = DL.mirror_update(F, Phi, y, 0.0, 0.05, 2000.0, rec_every=10 ** 9)   # the mode: the Bregman projection's dual
        D = F.dw(al)
        cov = (2.38 ** 2 / N) * T * np.linalg.inv(Phi @ (D[:, None] * Phi.T))
        Zs, Ws, acc = DL.gibbs_rwm(F, Phi, y, T, Z[-1], cov, N_RWM, seed=seed, thin=10)
        Ws = Ws[int(BURN_A * len(Ws)):]
        return dict(kind=kind, h=h, mean=Ws.mean(0), sd=Ws.std(0), se=batch_se(Ws), fit=(Ws @ Phi.T).mean(0) - y, acc=acc,
                    secs=time.time() - t0, w_proj=W[-1], rates=np.linalg.eigvalsh(Phi @ (D[:, None] * Phi.T) / N))
    if kind == 'mirror':
        ts, Z, Ws, al = DL.mirror_update(F, Phi, y, SIG_A, h, T_A, seed=seed, rec_every=10)
        Ws = Ws[int(BURN_A * len(Ws)):]
        return dict(kind=kind, h=h, mean=Ws.mean(0), sd=Ws.std(0), se=batch_se(Ws), fit=(Ws @ Phi.T).mean(0) - y, secs=time.time() - t0)
    if kind == 'path':                                                     # Heun on α against the mirror update, the same noise
        Pk = np.eye(d) - np.linalg.pinv(Phi) @ Phi                          # the projector onto ker Φ: P_{ker Φ}η is conserved
        out = dict(kind=kind, h=h)
        for name, T_end in (('short', T_PATH_A), ('long', T_DRIFT_A)):
            rec = max(1, int(round(1.0 / h)))
            ts, Z, Wm, al = DL.mirror_update(F, Phi, y, SIG_A, h, T_end, seed=seed, rec_every=rec)
            ts2, A, B = DL.primal_sde(L, Phi, y, np.zeros(d), np.log(EPS_A * b), SIG_A, h, T_end, seed=seed, rec_every=rec,
                                      scheme='stratonovich')
            eta_h = F.eta(np.abs(A)) * np.sign(A)
            out[name] = dict(dw=float(np.abs(np.exp(B) * L.amp(A)[0] / L.m - Wm).max()), deta=float(np.abs(eta_h - Z @ Phi).max()),
                             inv=np.abs(eta_h @ Pk).max(1), ts=ts)
        out['secs'] = time.time() - t0
        return out
    # 'ito': the paradigm's Euler update on (α, log β)
    ts, A, B = DL.primal_sde(L, Phi, y, np.zeros(d), np.log(EPS_A * b), SIG_A, h, T_A / 10, seed=seed, rec_every=10, scheme='ito')
    Ws = np.exp(B) * L.amp(A)[0] / L.m
    Qs = L.Q(A, np.exp(B))
    k1 = max(1, int(round(1.0 / (10 * h))))                                # t = 1
    return dict(kind=kind, h=h, mean=Ws[len(Ws) // 2:].mean(0), sd=Ws[len(Ws) // 2:].std(0), dq_end=Qs[-1] - q,
                rate1=(Qs[k1] - q) / ts[k1], drift0=DL.ito_q_drift(L, np.zeros(d), EPS_A * b, Phi, SIG_A), t_end=float(ts[-1]),
                secs=time.time() - t0)


def part_A(pool):
    jobs = [('rwm', 0.0, 0)] + [('mirror', h, 2) for h in H_A] + [('path', h, 5) for h in H_PATH_A] + [('ito', H_A[0], 3)]
    out = {}
    for r in pool.map(a_job, jobs):
        out[f"{r['kind']} {r['h']}"] = r
        print(f"  A {r['kind']} h {r['h']}: {r['secs']:.0f}s", flush=True)
    return out


# --- B ---------------------------------------------------------------------------------------------------------------------------------
@cached
def b_job(args):
    seed, sig, h = args
    warnings.filterwarnings('ignore')
    L, Phi, y, b, rho, path, q = setup(seed, EPS_B)
    N, d = Phi.shape
    S = np.flatnonzero(path['w_end'])
    F, z0, al0, w0 = plateau(L, Phi, y, q, path, EPS_B)
    T = sig ** 2 / (2 * N)
    D = F.dw(al0)
    t0 = time.time()
    ts, Z, Ws, al = DL.mirror_update(F, Phi, y, sig, h, TW_B, seed=100 + seed, rec_every=1, zeta0=z0, al0=al0)
    Ws = Ws[len(Ws) // 20:]
    C = np.cov(Ws[:, S].T)
    G = Phi[:, S].T @ Phi[:, S]
    sd = np.sqrt(np.diag(C))
    off = np.setdiff1d(np.arange(d), S)
    rates = np.sort(np.linalg.eigvals((G / N) @ np.diag(D[S])).real)
    return dict(seed=seed, sig=sig, h=h, S=S, ratio=np.diag(C) / (T * D[S]), corr=(C / np.outer(sd, sd))[np.triu_indices(S.size, 1)],
                dcorr=(G / np.sqrt(np.outer(np.diag(G), np.diag(G))))[np.triu_indices(S.size, 1)],
                bias=Ws[:, S].mean(0) - path['w_end'][S], se=batch_se(Ws[:, S]), sd=sd, w_bp=path['w_end'][S],
                inactive=float(np.abs(Ws[:, off]).max()) if off.size else 0.0, start=float(np.abs(w0 - path['w_end']).max()),
                rates=rates, n_relax=float(TW_B * rates[0]), secs=time.time() - t0)


def part_B(pool):
    jobs = [(s, g, H_B) for s in SEEDS_B for g in SIGS_B] + [(SEEDS_B[0], 0.1, H_B / 2)]
    rows = pool.map(b_job, jobs)
    for r in rows:
        print(f"  B design {r['seed']} σ {r['sig']} h {r['h']}: Var/(T·dw/dη) {np.round(r['ratio'], 3)} ({r['secs']:.0f}s)", flush=True)
    return rows


# --- C ---------------------------------------------------------------------------------------------------------------------------------
@cached
def c_job(args):
    eps, sig, k = args                                                     # sig = 0: the noiseless run on the same grid
    warnings.filterwarnings('ignore')
    L, Phi, y, b, rho, path, q = setup(SEED_CD, eps)
    N, d = Phi.shape
    F = DL.FastFree(L, q)
    Lam = RT.transit_scale(L, eps)
    t_end = 1.5 * path['taus'][-1] * Lam
    nst = int(round(t_end / H_C))
    xi = np.random.default_rng(1000 + k).standard_normal((nst + 1, N))
    t0 = time.time()
    ts, Z, W, al = DL.mirror_update(F, Phi, y, sig, H_C, t_end, rec_every=REC_C, noise=lambda j: xi[j])
    taus = ts / Lam
    U = Z @ Phi / Lam
    th = 1e-3 * np.abs(path['w_end']).max()
    ent = X.entries_of(W, taus, th)
    out = dict(eps=eps, sig=sig, k=k, taus=taus, U=U, entries=ent, secs=time.time() - t0, Lam=Lam,
               w_end=W[-1], iss=float(np.abs(U - np.array([MF.iss_eval(path, t, Phi, y)[0] for t in taus])).max()))
    if sig > 0:
        # the first-order response on the record grid: increments of W in τ, summed over each record block
        blocks = np.add.reduceat(xi[:nst], np.arange(0, nst, REC_C), axis=0)[:taus.size - 1] * np.sqrt(H_C / Lam)
        out['U_lin'] = DL.linear_response(path, Phi, taus, blocks, sig / (N * np.sqrt(Lam)))
    return out


def part_C(pool):
    jobs = [(e, 0.0, 0) for e in EPS_C] + [(e, s, k) for e in EPS_C for s in SIGS_C for k in range(PATHS_C)]
    res = pool.map(c_job, jobs)
    det = {r['eps']: r for r in res if r['sig'] == 0}
    rows = []
    for r in res:
        if r['sig'] == 0:
            continue
        dd = det[r['eps']]
        dU = r['U'] - dd['U']
        err = np.abs(dU - r['U_lin'])
        ev = [t for t, i, kind in setup(SEED_CD, r['eps'])[5]['events']]
        out = {}
        for w in WIN_C:
            keep = np.all(np.abs(r['taus'][:, None] - np.array(ev)[None, :]) >= w, axis=1)
            out[w] = float(err[keep].max() / np.abs(dU).max())
        kap = r['sig'] / (X.NA * np.sqrt(r['Lam']))
        row = dict(eps=r['eps'], sig=r['sig'], k=r['k'], share=float(np.abs(dU).max()), kappa=float(kap),
                   lin_err=float(err.max() / np.abs(dU).max()), lin_out=out, lin_l2=float(np.sqrt(np.mean(err ** 2) / np.mean(dU ** 2))),
                   iss=r['iss'], iss_det=dd['iss'],
                   order=[i for i in np.argsort(r['entries']) if np.isfinite(r['entries'][i])] ==
                         [i for i in np.argsort(dd['entries']) if np.isfinite(dd['entries'][i])],
                   dt_entry=float(np.nanmax(np.abs(np.where(np.isfinite(r['entries']), r['entries'], np.nan) -
                                                    np.where(np.isfinite(dd['entries']), dd['entries'], np.nan)))),
                   end=float(np.abs(r['w_end'] - dd['w_end']).max()))
        rows.append(row)
        print(f"  C ε {r['eps']:g} σ {r['sig']} path {r['k']}: share {row['share']:.2e} = {row['share'] / kap:.2f}κ, linear response "
              f"error {row['lin_err']:.3f} (outside the windows {', '.join(f'{v:.3f}' for v in out.values())}; L² {row['lin_l2']:.3f}), "
              f"order kept {row['order']} ({r['secs']:.0f}s)", flush=True)
    return dict(rows=rows, det={e: dict(iss=det[e]['iss'], entries=det[e]['entries']) for e in det})


# --- D ---------------------------------------------------------------------------------------------------------------------------------
@cached
def d_job(args):
    eps, sig = args
    warnings.filterwarnings('ignore')
    L, Phi, y, b, rho, path, q = setup(SEED_CD, eps)
    N, d = Phi.shape
    S = np.flatnonzero(path['w_end'])
    off = np.setdiff1d(np.arange(d), S)
    F, z0, al0, w0 = plateau(L, Phi, y, q, path, eps)
    t0 = time.time()
    ts, Z, Ws, al = DL.mirror_update(F, Phi, y, sig, H_D, TW_D, seed=7, rec_every=1, zeta0=z0, al0=al0)
    Ws = Ws[len(Ws) // 20:]
    th = 1e-3 * np.abs(path['w_end']).max()
    return dict(eps=eps, sig=sig, bias=Ws[:, S].mean(0) - path['w_end'][S], se=batch_se(Ws[:, S]), sd=Ws[:, S].std(0),
                inactive=float(np.abs(Ws[:, off]).max()), frac=float(np.mean(np.abs(Ws[:, off]).max(1) > th)),
                fit=float(np.abs(Ws.mean(0) @ Phi.T - y).max()), start=float(np.abs(w0 - path['w_end']).max()), secs=time.time() - t0)


def part_D(pool):
    rows = pool.map(d_job, [(e, s) for e in EPS_D for s in SIGS_D])
    for r in rows:
        print(f"  D ε {r['eps']:g} σ {r['sig']}: bias/SE {np.round(r['bias'] / r['se'], 2)}, inactive max {r['inactive']:.1e} "
              f"({r['secs']:.0f}s)", flush=True)
    return rows


# --- report ----------------------------------------------------------------------------------------------------------------------------
def fmt(v, p=3):
    return '(' + ', '.join(f'{x:.{p}g}' for x in np.atleast_1d(v)) + ')'


def report(A=None, B=None, C=None, D=None):
    s = ['# Proposition D‴.50: the rich transit under label noise — the dual is a Langevin diffusion', '',
         'Generated by `experiments/exp_transit_noise.py` (the designs of `exp_mirror_transit`; m = 12, G = 1, γ = 0.1, a_max = 1, '
         'c = 0.025; the label noise dζ = (y − Φw)dt/N + (σ/N)dW, T = σ²/(2N)).', '']
    if A is not None:
        L, Phi, y, b, q = small_dictionary()
        r = A['rwm 0.0']
        s += ['## A. The Gibbs law (a 4 × 8 dictionary with orthonormal rows, masses 2b_i, σ = 1, T = 0.125)', '',
              f"Random-walk Metropolis on π ∝ exp(−F/T) ({N_RWM} steps, acceptance {r['acc']:.2f}); the mirror update over t = {T_A:g} after "
              f"a burn-in of {BURN_A:.0%}; the linearised dual's rates at the projection {fmt(r['rates'])}. The coefficients' means and "
              'standard deviations, the largest deviation from the sampler\'s in units of the combined standard error (batch means) and '
              'relatively, the mean fit.', '',
              '| dynamics | E[w] | sd(w) | max_i |ΔE[w_i]|/se | max_i |Δsd_i|/sd_i | |E[Φw] − y|_∞ |', '|---|---|---|---|---|---|',
              f"| Metropolis (π) | {fmt(r['mean'])} | {fmt(r['sd'])} | — | — | {np.abs(r['fit']).max():.1e} |"]
        for h in H_A:
            x = A[f'mirror {h}']
            z = np.abs(x['mean'] - r['mean']) / np.sqrt(x['se'] ** 2 + r['se'] ** 2)
            s.append(f"| mirror update, h = {h} | {fmt(x['mean'])} | {fmt(x['sd'])} | {z.max():.1f} | {np.abs(x['sd'] / r['sd'] - 1).max():.3f} | "
                     f"{np.abs(x['fit']).max():.1e} |")
        s += ['', f"The noiseless projection: {fmt(r['w_proj'])}.", '',
              'The paradigm\'s law read in Stratonovich\'s sense (stochastic Heun on α, β on its curve) against the mirror update with the '
              'same noise: the largest differences in w and η, and the conserved directions\' drift |P_{ker Φ}η|_∞ at the end —', '',
              f'| h | t ≤ {T_PATH_A:g}: max|Δw|, max|Δη| | t ≤ {T_DRIFT_A:g}: max|Δw|, max|Δη| | |P_{{ker Φ}}η| at t = {T_PATH_A:g}, '
              f'{T_DRIFT_A:g} |', '|---|---|---|---|']
        for h in H_PATH_A:
            x = A[f'path {h}']
            s.append(f"| {h} | {x['short']['dw']:.1e}, {x['short']['deta']:.1e} | {x['long']['dw']:.1e}, {x['long']['deta']:.1e} | "
                     f"{x['short']['inv'][-1]:.1e}, {x['long']['inv'][-1]:.1e} |")
        x = A[f'ito {H_A[0]}']
        s += ['', f"Read in Itô's sense (Euler on (α, log β), h = {x['h']}): Q(t) − Q(0) at t = {x['t_end']:g} from {x['dq_end'].min():.3g} to "
              f"{x['dq_end'].max():.3g}; the rate over the first time unit against (γβ/2)(γa² − Gβa′(1 + a²/a_max²))σ²|φ|²/N² at the start: "
              f"ratios {fmt(x['rate1'] / x['drift0'], 3)}; the coefficients' means {fmt(x['mean'])}, standard deviations {fmt(x['sd'])}.", '']
    if B is not None:
        s += ['## B. Fluctuation–dissipation on the support (ε = 10⁻⁶)', '',
              f'From the noiseless transit\'s end, a window t = {TW_B:g} (the first 5 % dropped). Var(w_i)/(T·dw_i/dη), the correlations '
              'of the fluctuations against the design\'s (the cosines of the support\'s columns), the mean\'s offset from weighted basis '
              'pursuit in standard errors, the inactive coefficients\' largest |w|, the window in relaxation times of the slowest mode.', '',
              '| design | σ | h | support | Var/(T·dw/dη) | corr | design corr | (E[w] − w_BP)/se | max inactive |w| | relaxation times |',
              '|---|---|---|---|---|---|---|---|---|---|']
        for r in B:
            s.append(f"| {r['seed']} | {r['sig']} | {r['h']} | {r['S'].tolist()} | {fmt(r['ratio'])} | {fmt(r['corr'], 2)} | "
                     f"{fmt(r['dcorr'], 2)} | {fmt(r['bias'] / r['se'], 2)} | {r['inactive']:.1e} | {r['n_relax']:.0f} |")
        s.append('')
    if C is not None:
        rows = C['rows']
        s += ['## C. The noisy transit (the first design)', '',
              'Per ε and σ over four noise paths: the noise\'s share sup_τ|u − u_det| (mean, and divided by κ = σ/(NΛ^{1/2})), the '
              'first-order response\'s (δν ← P_T(δν + κdW) along the path) relative error |δu − δu_lin|/sup|δu| — its sup over all τ, '
              f'its sup outside ±{WIN_C[0]:g} and ±{WIN_C[1]:g} of the path\'s events, its L² over τ (the largest over the paths) —, the '
              'entries\' order kept, the largest entry-time shift; the noiseless run\'s sup|u_det − u_ISS| for comparison.', '',
              f'| ε | σ | κ | share (mean) | share/κ | linear response: sup | outside ±{WIN_C[0]:g} | outside ±{WIN_C[1]:g} | L² | order kept | '
              'max entry shift | sup|u_det − u_ISS| |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for e in EPS_C:
            for g in SIGS_C:
                rr = [r for r in rows if r['eps'] == e and r['sig'] == g]
                sh = np.array([r['share'] for r in rr])
                s.append(f"| {e:g} | {g} | {rr[0]['kappa']:.2e} | {sh.mean():.2e} | {np.mean(sh / rr[0]['kappa']):.2f} | "
                         f"{max(r['lin_err'] for r in rr):.3f} | {max(r['lin_out'][WIN_C[0]] for r in rr):.3f} | "
                         f"{max(r['lin_out'][WIN_C[1]] for r in rr):.3f} | {max(r['lin_l2'] for r in rr):.3f} | "
                         f"{sum(r['order'] for r in rr)}/{len(rr)} | {max(r['dt_entry'] for r in rr):.2e} | {rr[0]['iss_det']:.2e} |")
        sh = {g: np.array([np.mean([r['share'] for r in rows if r['eps'] == e and r['sig'] == g]) for e in EPS_C]) for g in SIGS_C}
        sl = {g: np.polyfit(np.log(EPS_C), np.log(sh[g]), 1)[0] for g in SIGS_C}
        s += ['', 'Log-log slope of the share in ε: ' + ', '.join(f'σ = {g}: {sl[g]:.3f}' for g in SIGS_C) + ' (κ ∝ ε^{1/4}).', '']
    if D is not None:
        s += ['## D. The stationary mean in the rich limit (the first design)', '',
              f'From the noiseless transit\'s end, a window t = {TW_D:g}: the support\'s mean against weighted basis pursuit (in standard '
              'errors and absolutely), the standard deviations, the inactive coefficients (largest |w| and the fraction of the time any '
              'exceeds 10⁻³max|w_BP|), the time-averaged fit.', '',
              '| ε | σ | (E[w_S] − w_BP)/se | |E[w_S] − w_BP|_∞ | sd(w_S) | max inactive |w| | fraction | |Φw̄ − y|_∞ |',
              '|---|---|---|---|---|---|---|---|']
        for r in D:
            s.append(f"| {r['eps']:g} | {r['sig']} | {fmt(r['bias'] / r['se'], 2)} | {np.abs(r['bias']).max():.1e} | {fmt(r['sd'])} | "
                     f"{r['inactive']:.1e} | {r['frac']:.3f} | {r['fit']:.1e} |")
        s.append('')
    return '\n'.join(s)


def main():
    raw = os.path.join(OUT, 'transit_noise_raw.npz')
    md = os.path.join(OUT, 'transit_noise.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCD']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCD')
    res = dict(old)
    t0 = time.time()
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D)):
            if k in parts:
                res[k] = fn(pool); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCD']))


if __name__ == '__main__':
    main()
