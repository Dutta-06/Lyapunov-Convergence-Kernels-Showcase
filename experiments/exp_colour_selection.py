"""
Proposition D‴.43: the colour of the measurement noise selects the kernels — from the trace (a weighted ℓ¹ norm) to the log-determinant
(a support count)
=====================================================================================================================================

A kernel dictionary sampled at N < d points (the design Φ), the coefficients w = u⊙u − v⊙v (the unbalanced Fisher–Rao law of a signed
kernel measure's two parts on a fixed dictionary), fitted by gradient descent with label noise; along the manifold M of exact fits the
parameters follow the reduced flow of `colour_selection` (R = tr Ψ(K): the trace for white draws, (1/τ_c)log det(I + τ_cK) for a slow
AR(1) noise).

A. The potential's identities at random points: det(I + τK) against its subset (Cauchy–Binet) expansion; ∂R/∂D_i = (4/N)ω_i; R's second
   differences along random segments in D (zero for white draws, negative for coloured ones: concave); on a k-sparse exact fit
   τ_cR − k log(4τ_c/N) → log(det(Φ_SᵀΦ_S)Π_{i∈S}|w_i|) as τ_c grows (the support's size times log τ_c).
B. Label-noise gradient descent (P_B paths, step H_B, per-step variance σ_B²) against the reduced flow with the scheme's filter at the step,
   on N_B × D_B Gaussian designs with a K_B-sparse truth: fresh draws and AR(1) draws of correlation time τ_c = 10, 100 (in time units);
   the paths' mean against the flow at nine times of the rescaled time s = 4v_ht/N, and the interpolant the noiseless gradient flow
   reaches from each path's end state (its error and support), against weighted basis pursuit (white's prediction). The start: the
   noiseless flow's interpolant from u = v = α (exactly: the hyperbolic entropy's minimiser on M, P = α²).
C. Recovery against sparsity (the reduced flows, continuum filters): the fraction of instances whose flow ends within 10⁻⁴ (relative)
   of the truth, for weighted basis pursuit, the white flow, slow flows (τ_c = 10, 100, 1000), and the white flow followed by the slow one
   (τ_c = 100) — N × d = 10 × 40 and 20 × 80.
D. Gaussian kernel dictionaries on [0, 1] (D_D centres, N_D equispaced samples, two widths — neighbouring atoms nearly collinear),
   random k-atom truths (signs ±1, magnitudes in [0.5, 1.5]): the recovery rates of basis pursuit and of the white and slow flows, and,
   on the first instances where basis pursuit fails, label-noise descent (fresh and AR(1) draws, τ_c = 100) from the same start.
E. Other colours at the step H_B: AR(1), a smooth noise (Gaussian step correlations) and a draw kept k steps, at comparable
   correlation times, against fresh draws (the reduced flow with the scheme's filter), on part C's 10 × 40 instances at k = 3, 4.
F. The start: the gradient flow's interpolant from u = v = α for α = 0.1, 0.5, 2 (the white flow's end does not depend on it, a slow
   flow's does: the potential is concave and has many local minima), on part C's 10 × 40 instances at k = 3, 4.

G. References on part C's 10 × 40 instances: reweighted ℓ¹ (Candès–Wakin–Boyd: weights 1/(|w_i| + ε), from basis pursuit's solution)
   and the majorisation–minimisation of the slow noise's own potential (R concave in D, so its linearisation majorises it: weighted ℓ¹
   with the weights ω(|w|), τ_c = 100 — Fazel–Hindi–Boyd's log-det heuristic), from basis pursuit's solution and from the flows' start:
   the same penalty minimised by linear programmes instead of descended by the noise.

H. Which exact fits the selection holds: at a k-sparse exact fit (P = 0, support S) the inactive atoms decay when the weighted dual
   certificate |c_i| < ω_i holds, c = Φ_Iᵀ(Φ_S⁺)ᵀ(ω_S⊙sign w_S) (sufficient); the largest |c_i|/ω_i at the truth and at basis pursuit's
   solution, for white draws and slow ones (τ_c = 10, 100, 1000), on part C's 10 × 40 instances — for white draws the weighted ℓ¹
   recovery condition, for slow ones ∝ 1/τ_c at a support smaller than N.

Output: results/chaos/colour_selection.md and colour_selection_raw.npz (CS_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`).
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
from src.core import colour_selection as CS
from src.core.runcache import cached
import exp_cohesion_leak as E

OUT = os.environ.get('CS_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('CS_WORKERS', 4))
ALPHA, TOL = 0.5, 1e-4                                                        # the start u = v = α; success within TOL·max|w*|
N_B, D_B, K_B, H_B, SIG_B, P_B, SEEDS_B = 10, 40, 4, 0.5, 0.1, 8, (0, 1, 2, 3)
TAUS_B, S_B = (0.0, 10.0, 100.0), {0.0: 300.0, 10.0: 150.0, 100.0: 400.0}  # the rescaled time each colour runs
SIZES_C, KS_C, N_INST_C = ((10, 40), (20, 80)), {(10, 40): (1, 2, 3, 4, 5, 6), (20, 80): (2, 4, 6, 8, 10)}, {(10, 40): 40, (20, 80): 20}
TAUS_C, S_WHITE, S_SLOW = (10.0, 100.0, 1000.0), 2000.0, lambda tc: 400.0 + 60.0 * np.sqrt(1 + tc)
D_D, N_D, WIDTHS_D, KS_D, N_INST_D, N_NOISY_D = 40, 20, (0.03, 0.05), (3, 4, 5, 6), 40, 3
KS_E, N_INST_E, TC_E = (3, 4), 20, 10.0
ALPHAS_F = (0.1, 0.5, 2.0)
KS_G, N_INST_G, EPS_G, IT_G, TC_G = (2, 3, 4, 5, 6), 40, 0.1, 12, 100.0
TAUS_H = (0.0, 10.0, 100.0, 1000.0)


def instance(N, d, k, seed):
    """A Gaussian design (columns of norm ≈ 1) and a k-sparse truth with signs ±1 and magnitudes in [0.5, 1.5]."""
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((N, d)) / np.sqrt(N)
    ws = np.zeros(d); idx = rng.choice(d, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.5, 1.5, k)
    return Phi, ws


def start(Phi, y, alpha):
    """The noiseless gradient flow's interpolant from u = v = α (exactly, `colour_selection.flow_interpolant`): P = α²."""
    d = Phi.shape[1]
    return CS.flow_interpolant(Phi, y, np.full(d, alpha), np.full(d, alpha))


def seed_of(N, k, i):
    return 100000 * N + 1000 * k + i


def ok(w, ws):
    return bool(np.abs(w - ws).max() < TOL * np.abs(ws).max())


def corr_ct(tc):
    return ('fresh',) if tc == 0 else ('ou', tc)


def part_A():
    rng = np.random.default_rng(1)
    rows = []
    for N, d in ((4, 7), (5, 9)):
        Phi = rng.standard_normal((N, d)) / np.sqrt(N); D = rng.uniform(0.05, 1.5, d)
        for tau in (0.3, 10.0, 1e3):
            K = CS.gram(Phi, D)
            exact = np.linalg.det(np.eye(N) + tau * K)
            e = 1e-6
            g = np.array([(CS.potential(Phi, D + e * np.eye(d)[i], 0, ('ou', tau)) - CS.potential(Phi, D - e * np.eye(d)[i], 0, ('ou', tau)))
                          / (2 * e) for i in range(d)])
            sec = []
            for _ in range(20):                                               # second differences along random segments
                D1, D2 = rng.uniform(0.05, 1.5, d), rng.uniform(0.05, 1.5, d)
                f = [CS.potential(Phi, D1 + t * (D2 - D1), 0, ('ou', tau)) for t in (0.0, 0.5, 1.0)]
                sec.append(f[0] + f[2] - 2 * f[1])
            fw = [CS.potential(Phi, D1 + t * (D2 - D1), 0, ('fresh',)) for t in (0.0, 0.5, 1.0)]
            rows.append(dict(N=N, d=d, tau=tau, cb=abs(CS.det_subsets(Phi, D, tau) / exact - 1),
                             grad=float(np.abs(g - (4 / N) * CS.weights(Phi, D, 0, ('ou', tau))).max() / np.abs(g).max()),
                             sec_max=float(max(sec)), sec_min=float(min(sec)), sec_white=float(abs(fw[0] + fw[2] - 2 * fw[1]))))
    lim = []
    for N, d, k in ((10, 40, 3), (10, 40, 6), (20, 80, 5)):
        Phi, ws = instance(N, d, k, 7 + k)
        S = np.flatnonzero(ws)
        ref = np.log(np.linalg.det(Phi[:, S].T @ Phi[:, S]) * np.prod(np.abs(ws[S])))
        lim.append(dict(N=N, d=d, k=k, ref=ref, vals=[(tc, tc * CS.potential(Phi, np.abs(ws), 0, ('ou', tc)) - k * np.log(4 * tc / N))
                                                      for tc in (1e2, 1e4, 1e6)]))
    return dict(rows=rows, lim=lim)


@cached
def b_job(args):
    """Part B, one (seed, τ_c): the paths' mean against the reduced flow at nine times, the end interpolants, basis pursuit."""
    seed, tc = args
    t0 = time.time()
    Phi, ws = instance(N_B, D_B, K_B, seed)
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    u0, v0 = CS.uv_of(w0, P0)
    corr = ('fresh',) if tc == 0 else ('ar1', float(np.exp(-H_B / tc)))
    rho = 0.0 if tc == 0 else float(np.exp(-H_B / tc))
    vh = H_B * SIG_B ** 2 / N_B
    S = S_B[tc]
    n = int(round(S * N_B / (4 * vh) / H_B)); every = n // 8
    s_rec = np.arange(9) * every * H_B * 4 * vh / N_B
    _, Wf, _ = CS.reduced_flow(Phi, w0, P0, s_rec[-1], H_B, corr, t_eval=s_rec)
    u, v, rec = CS.label_noise_gd(Phi, y, u0, v0, H_B, SIG_B, rho, n, np.random.default_rng(seed + 77), P_B, rec_every=every)
    ends, conv = CS.end_interpolants(Phi, y, u, v)
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    print(f"[{time.time() - t0:.0f}s] B seed {seed} tau_c {tc:g}", flush=True)
    return dict(seed=seed, tc=tc, s=s_rec, dev=np.abs(rec - Wf.T).max(axis=1), flow_err=np.abs(Wf.T - ws[None, :]).max(axis=1),
                path_err=np.abs(rec - ws[None, :]).max(axis=1), end_err=np.abs(ends - ws[None, :]).max(axis=1),
                end_nnz=np.where(conv, (np.abs(ends) > 1e-3).sum(axis=1), -1), bp_err=float(np.abs(wbp - ws).max()), n=n, vh=vh,
                wmax=float(np.abs(ws).max()))


@cached
def c_job(args):
    """Part C, one instance: success of basis pursuit, the white flow, the slow flows and the annealed one; the supports' sizes."""
    N, d, k, i = args
    Phi, ws = instance(N, d, k, seed_of(N, k, i))
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    nnz = lambda w: int((np.abs(w) > 1e-3 * np.abs(ws).max()).sum())
    out = {}
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0)); out['bp'] = (ok(wbp, ws), nnz(wbp))
    _, Ww, Pw = CS.reduced_flow(Phi, w0, P0, S_WHITE, 0, ('fresh',)); out['white'] = (ok(Ww[:, -1], ws), nnz(Ww[:, -1]))
    for tc in TAUS_C:
        _, Wt, _ = CS.reduced_flow(Phi, w0, P0, S_SLOW(tc), 0, ('ou', tc)); out[f'slow {tc:g}'] = (ok(Wt[:, -1], ws), nnz(Wt[:, -1]))
    _, Wa, _ = CS.reduced_flow(Phi, Ww[:, -1], np.maximum(Pw[:, -1], 1e-300), S_SLOW(100.0), 0, ('ou', 100.0))
    out['white, then slow 100'] = (ok(Wa[:, -1], ws), nnz(Wa[:, -1]))
    return dict(N=N, d=d, k=k, i=i, out=out)


def dict_instance(width, k, i):
    """A Gaussian kernel dictionary (D_D centres on [0, 1], N_D equispaced samples, unit columns) and a random k-atom truth."""
    Phi = CS.rbf_dictionary(np.linspace(0.0, 1.0, N_D), np.linspace(0.0, 1.0, D_D), width)
    rng = np.random.default_rng(9000 + 100 * k + i)
    ws = np.zeros(D_D); idx = rng.choice(D_D, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.5, 1.5, k)
    return Phi, ws


@cached
def d_job(args):
    """Part D, one (width, k, instance): success of basis pursuit and of the white and slow flows; the supports' sizes."""
    width, k, i = args
    Phi, ws = dict_instance(width, k, i)
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    nnz = lambda w: int((np.abs(w) > 1e-3 * np.abs(ws).max()).sum())
    out = {}
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0)); out['bp'] = (ok(wbp, ws), nnz(wbp))
    _, Ww, _ = CS.reduced_flow(Phi, w0, P0, S_WHITE, 0, ('fresh',)); out['white'] = (ok(Ww[:, -1], ws), nnz(Ww[:, -1]))
    for tc in (10.0, 100.0):
        _, Wt, _ = CS.reduced_flow(Phi, w0, P0, S_SLOW(tc), 0, ('ou', tc)); out[f'slow {tc:g}'] = (ok(Wt[:, -1], ws), nnz(Wt[:, -1]))
    return dict(width=width, k=k, i=i, out=out)


@cached
def d_noisy_job(args):
    """Part D, label-noise descent on one instance (fresh draws and AR(1) τ_c = 100, P_B paths, step H_B, σ_B): the end interpolants'
    errors and supports, against the reduced flows'."""
    width, k, i = args
    t0 = time.time()
    Phi, ws = dict_instance(width, k, i)
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    u0, v0 = CS.uv_of(w0, P0)
    lam = np.linalg.eigvalsh(CS.gram(Phi, np.sqrt(w0 ** 2 + 4 * P0 ** 2))).max()
    h = min(H_B, 0.5 / lam)                                                   # hλ_max ≤ 1/2 at the start
    vh = h * SIG_B ** 2 / N_D
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    out = dict(width=width, k=k, i=i, h=h)
    for tc in (0.0, 100.0):
        rho = 0.0 if tc == 0 else float(np.exp(-h / tc))
        S = S_B[tc]
        n = int(round(S * N_D / (4 * vh) / h))
        u, v, _ = CS.label_noise_gd(Phi, y, u0, v0, h, SIG_B, rho, n, np.random.default_rng(1000 * k + i), P_B)
        ends, conv = CS.end_interpolants(Phi, y, u, v)
        raw = u * u - v * v
        flow = CS.reduced_flow(Phi, w0, P0, S, h, ('fresh',) if tc == 0 else ('ar1', rho))[1][:, -1]
        out[tc] = dict(n=n, err=np.abs(ends - ws[None, :]).max(axis=1),
                       nnz=np.where(conv, (np.abs(ends) > 1e-3 * np.abs(ws).max()).sum(axis=1), -1),
                       raw_err=np.abs(raw - ws[None, :]).max(axis=1), raw_bp=np.abs(raw - wbp[None, :]).max(axis=1),
                       alive=((u * u + v * v) > 1e-8 * (u * u + v * v).max(axis=1, keepdims=True)).sum(axis=1),
                       flow_err=float(np.abs(flow - ws).max()), flow_nnz=int((np.abs(flow) > 1e-3 * np.abs(ws).max()).sum()),
                       flow_bp=float(np.abs(flow - wbp).max()))
    print(f"[{time.time() - t0:.0f}s] D noisy width {width:g} k {k} instance {i}", flush=True)
    return out


def d_noisy_cases():
    """The first N_NOISY_D instances of width WIDTHS_D[1], k = 4 where basis pursuit fails."""
    cases, i = [], 0
    while len(cases) < N_NOISY_D and i < N_INST_D:
        Phi, ws = dict_instance(WIDTHS_D[1], 4, i)
        if not ok(CS.basis_pursuit(Phi, Phi @ ws, (Phi ** 2).sum(axis=0)), ws):
            cases.append((WIDTHS_D[1], 4, i))
        i += 1
    return cases


def e_corrs():
    k = int(round(TC_E / H_B))
    return {'fresh': ('fresh',), f'AR(1), τ_c = {TC_E:g}': ('ar1', float(np.exp(-H_B / TC_E))),
            f'smooth (Gaussian), τ_c = {TC_E:g}': ('gauss', TC_E / H_B), f'kept {2 * k} steps': ('block', 2 * k)}


@cached
def e_job(args):
    """Part E, one instance: success of the reduced flow with each colour's filter at the step H_B."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    out = {}
    for name, corr in e_corrs().items():
        _, W, _ = CS.reduced_flow(Phi, w0, P0, S_WHITE if corr[0] == 'fresh' else S_SLOW(TC_E), H_B, corr)
        out[name] = ok(W[:, -1], ws)
    return dict(k=k, i=i, out=out)


@cached
def f_job(args):
    """Part F, one instance: success of the white and slow (τ_c = 100) flows from the interpolants of u = v = α."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    out = {}
    for a in ALPHAS_F:
        w0, P0 = start(Phi, y, a)
        _, Ww, _ = CS.reduced_flow(Phi, w0, P0, S_WHITE, 0, ('fresh',))
        _, Wt, _ = CS.reduced_flow(Phi, w0, P0, S_SLOW(100.0), 0, ('ou', 100.0))
        out[a] = (ok(w0, ws), ok(Ww[:, -1], ws), ok(Wt[:, -1], ws))
    return dict(k=k, i=i, out=out)


def mm_logdet(Phi, y, w, tc, iters=IT_G):
    """The majorisation–minimisation of R = (1/τ_c)log det(I + τ_cK) on M at P = 0: w ← argmin Σ_iω_i(|w|)|w_i| subject to Φw = y."""
    for _ in range(iters):
        w = CS.basis_pursuit(Phi, y, CS.weights(Phi, np.abs(w) + 1e-12, 0, ('ou', tc)))
    return w


def reweighted_l1(Phi, y, w, eps=EPS_G, iters=IT_G):
    """Candès–Wakin–Boyd: w ← argmin Σ_i|w_i|/(|w_i^{old}| + ε) subject to Φw = y."""
    for _ in range(iters):
        w = CS.basis_pursuit(Phi, y, 1.0 / (np.abs(w) + eps))
    return w


@cached
def g_job(args):
    """Part G, one instance: success of reweighted ℓ¹ and of the log-det's majorisation–minimisation (from basis pursuit and from the
    flows' start), with their support sizes."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    w0, _ = start(Phi, y, ALPHA)
    nnz = lambda w: int((np.abs(w) > 1e-3 * np.abs(ws).max()).sum())
    eps = EPS_G * np.abs(ws).max()
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    out = {'bp': (ok(wbp, ws), nnz(wbp))}
    w = reweighted_l1(Phi, y, wbp, eps); out['reweighted ℓ¹'] = (ok(w, ws), nnz(w))
    w = mm_logdet(Phi, y, wbp, TC_G); out['log-det MM from basis pursuit'] = (ok(w, ws), nnz(w))
    w = mm_logdet(Phi, y, w0, TC_G); out['log-det MM from the start'] = (ok(w, ws), nnz(w))
    return dict(k=k, i=i, out=out)


def certificate(Phi, w, tc):
    """The largest |c_i|/ω_i over the atoms off w's support, c = Φ_Iᵀ(Φ_S⁺)ᵀ(ω_S⊙sign w_S), ω at D = |w| (white: tc = 0)."""
    S = np.flatnonzero(np.abs(w) > 1e-9 * np.abs(w).max()); I = np.setdiff1d(np.arange(Phi.shape[1]), S)
    om = CS.weights(Phi, np.abs(w) + 1e-300, 0, corr_ct(tc))
    c = Phi[:, I].T @ np.linalg.pinv(Phi[:, S]).T @ (om[S] * np.sign(w[S]))
    return float(np.max(np.abs(c) / om[I])), len(S)


@cached
def h_job(args):
    """Part H, one instance: the certificate ratios at the truth and at basis pursuit's solution, for each colour; basis pursuit's success."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    wbp = np.where(np.abs(wbp) > 1e-9 * np.abs(wbp).max(), wbp, 0.0)
    return dict(k=k, i=i, bp=ok(wbp, ws), truth={tc: certificate(Phi, ws, tc)[0] for tc in TAUS_H},
                at_bp={tc: certificate(Phi, wbp, tc)[0] for tc in TAUS_H}, bp_nnz=certificate(Phi, wbp, 0.0)[1])


def fmt_end(x):
    return 'n.c.' if not np.isfinite(x) else f'{x:.0e}'


def report(A=None, B=None, C=None, D=None, Ee=None, Ff=None, G=None, H=None):
    L_ = ["# Proposition D‴.43: the colour of the measurement noise selects the kernels (`exp_colour_selection.py`)", "",
          "w = u⊙u − v⊙v on a fixed dictionary Φ (N × d), fitted by gradient descent with label noise; along the exact fits the reduced flow "
          "dP/ds = −ω⊙P, dw/ds = −Π_D(ω⊙w), ω_i = φ_iᵀf(K)φ_i, the gradient flow of R = tr Ψ(K) (`colour_selection`). Success: within "
          f"{TOL:g}·max|w*| of the truth. The start: the noiseless gradient flow's interpolant from u = v = {ALPHA:g} (the hyperbolic "
          "entropy's minimiser on M, P = α²); the end interpolants likewise from each path's end state.", ""]
    if A is not None:
        L_ += ["## A. The potential's identities", "",
               "| N × d | τ | det(I + τK) against its subset expansion (relative) | ∂R/∂D against (4/N)ω (relative) | R's second differences "
               "along 20 random segments (max, min) | white: |second difference| |", "|---|---|---|---|---|---|"]
        for r in A['rows']:
            L_.append(f"| {r['N']} × {r['d']} | {r['tau']:g} | {r['cb']:.1e} | {r['grad']:.1e} | {r['sec_max']:+.2e}, {r['sec_min']:+.2e} | "
                      f"{r['sec_white']:.1e} |")
        L_ += ["", "On a k-sparse exact fit (P = 0, D = |w|): τ_cR − k log(4τ_c/N) against log(det(Φ_SᵀΦ_S)Π_{i∈S}|w_i|) —", "",
               "| N × d | k | the limit | τ_c = 10², 10⁴, 10⁶ |", "|---|---|---|---|"]
        for r in A['lim']:
            L_.append(f"| {r['N']} × {r['d']} | {r['k']} | {r['ref']:+.6f} | {', '.join(f'{v:+.6f}' for _, v in r['vals'])} |")
        L_.append("")
    if B is not None:
        L_ += [f"## B. Label-noise descent against the reduced flow ({N_B} × {D_B}, a {K_B}-sparse truth, step {H_B:g}, σ = {SIG_B:g}, {P_B} paths)", "",
               "| seed | the noise | steps | s at the end | the paths' mean against the flow (max over w, at s = 0, S/8, …, S) | the flow's error "
               "at S | the paths' mean's error at S | the end interpolants' errors (per path) | their supports' sizes | basis pursuit's error |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for r in B:
            noise = 'fresh' if r['tc'] == 0 else f"AR(1), τ_c = {r['tc']:g}"
            L_.append(f"| {r['seed']} | {noise} | {r['n']} | {r['s'][-1]:.0f} | "
                      f"{', '.join(f'{x:.3f}' for x in r['dev'])} | {r['flow_err'][-1]:.1e} | {r['path_err'][-1]:.1e} | "
                      f"{', '.join(fmt_end(x) for x in r['end_err'])} | {', '.join(str(int(x)) if x >= 0 else 'n.c.' for x in r['end_nnz'])} | "
                      f"{r['bp_err']:.2f} |")
        L_.append("")
    if C is not None:
        meth = list(C[0]['out'])
        L_ += ["## C. Recovery against sparsity (the reduced flows, continuum filters)", "",
               "| N × d | k | instances | " + ' | '.join(meth) + " | the slow flow's (τ_c = 100) mean support size | basis pursuit's |",
               "|---|---|---|" + "---|" * len(meth) + "---|---|"]
        keys = sorted({(r['N'], r['d'], r['k']) for r in C})
        for N, d, k in keys:
            rs = [r for r in C if (r['N'], r['d'], r['k']) == (N, d, k)]
            L_.append(f"| {N} × {d} | {k} | {len(rs)} | " + ' | '.join(f"{np.mean([r['out'][m][0] for r in rs]):.2f}" for m in meth) +
                      f" | {np.mean([r['out']['slow 100'][1] for r in rs]):.1f} | {np.mean([r['out']['bp'][1] for r in rs]):.1f} |")
        L_.append("")
    if D is not None:
        Dr, Dn = D['rates'], D['noisy']
        meth = list(Dr[0]['out'])
        L_ += [f"## D. Gaussian kernel dictionaries on [0, 1] ({D_D} centres, {N_D} equispaced samples), random k-atom truths ({N_INST_D} per row)", "",
               "| width | the largest off-diagonal |φ_iᵀφ_j| | k | " + ' | '.join(meth) + " | mean support size: basis pursuit, slow 100 |",
               "|---|---|---|" + "---|" * len(meth) + "---|"]
        for width in WIDTHS_D:
            Phi = dict_instance(width, 3, 0)[0]
            coh = float(np.abs(Phi.T @ Phi - np.eye(D_D)).max())
            for k in KS_D:
                rs = [r for r in Dr if r['width'] == width and r['k'] == k]
                L_.append(f"| {width:g} | {coh:.2f} | {k} | " + ' | '.join(f"{np.mean([r['out'][m][0] for r in rs]):.2f}" for m in meth) +
                          f" | {np.mean([r['out']['bp'][1] for r in rs]):.1f}, {np.mean([r['out']['slow 100'][1] for r in rs]):.1f} |")
        L_ += ["", f"Label-noise descent on the first instances where basis pursuit fails (width {WIDTHS_D[1]:g}, k = 4; {P_B} paths, σ = {SIG_B:g}, "
               "the step min(0.5, 1/(2λ_max)) at the start):", "",
               "| instance | step | the noise | steps | the end interpolants' errors (per path) | their supports' sizes | the ends' atoms "
               "alive (D > 10⁻⁸ max D) | the raw ends' distance to the truth, to basis pursuit's solution (max over paths) | the reduced "
               "flow's error, support, distance to basis pursuit's solution |", "|---|---|---|---|---|---|---|---|---|"]
        for r in Dn:
            for tc in (0.0, 100.0):
                q = r[tc]
                L_.append(f"| {r['i']} | {r['h']:.3g} | {'fresh' if tc == 0 else 'AR(1), τ_c = 100'} | {q['n']} | {', '.join(fmt_end(x) for x in q['err'])} | "
                          f"{', '.join(str(int(x)) if x >= 0 else 'n.c.' for x in q['nnz'])} | {', '.join(str(int(x)) for x in q['alive'])} | "
                          f"{q['raw_err'].max():.2f}, {q['raw_bp'].max():.2f} | {q['flow_err']:.1e}, {q['flow_nnz']}, {q['flow_bp']:.1e} |")
        L_ += ["", "n.c.: the end's flow interpolant does not converge — the exact fit needs atoms whose masses have underflowed (fewer than N "
               "alive, a residual outside their span), its maximiser is at infinity; the raw end (residual O(σ)) is reported instead.", ""]
    if Ee is not None:
        names = list(Ee[0]['out'])
        L_ += [f"## E. Other colours at the step {H_B:g} (the reduced flow with each scheme's filter; 10 × 40, {N_INST_E} instances)", "",
               "| k | " + ' | '.join(names) + " |", "|---|" + "---|" * len(names)]
        for k in KS_E:
            rs = [r for r in Ee if r['k'] == k]
            L_.append(f"| {k} | " + ' | '.join(f"{np.mean([r['out'][n_] for r in rs]):.2f}" for n_ in names) + " |")
        L_.append("")
    if Ff is not None:
        L_ += ["## F. The start: the interpolant of u = v = α (10 × 40)", "",
               "| k | α | the interpolant itself | the white flow | the slow flow (τ_c = 100) |", "|---|---|---|---|---|"]
        for k in KS_E:
            rs = [r for r in Ff if r['k'] == k]
            for a in ALPHAS_F:
                L_.append(f"| {k} | {a:g} | " + ' | '.join(f"{np.mean([r['out'][a][j] for r in rs]):.2f}" for j in range(3)) + " |")
        L_.append("")
    if G is not None:
        meth = list(G[0]['out'])
        L_ += [f"## G. References: the same penalties minimised by linear programmes (10 × 40, {N_INST_G} instances; {IT_G} iterations; "
               f"ε = {EPS_G:g}·max|w*|; τ_c = {TC_G:g})", "",
               "| k | " + ' | '.join(meth) + " | their mean support sizes |", "|---|" + "---|" * len(meth) + "---|"]
        for k in KS_G:
            rs = [r for r in G if r['k'] == k]
            L_.append(f"| {k} | " + ' | '.join(f"{np.mean([r['out'][m][0] for r in rs]):.2f}" for m in meth) + " | " +
                      ', '.join(f"{np.mean([r['out'][m][1] for r in rs]):.1f}" for m in meth) + " |")
        L_.append("")
    if H is not None:
        L_ += ["## H. Which exact fits the selection holds: the dual certificate max_{i∉S}|c_i|/ω_i (below 1: the fit is stable)", "",
               "| k | instances | at the truth: the fraction certified, white; τ_c = 10, 100, 1000 | the median ratio at the truth, white; τ_c = 10, 100, 1000 | "
               "where basis pursuit fails: its solution's support size, and the fraction of its solutions certified for τ_c = 10, 100, 1000 |",
               "|---|---|---|---|---|"]
        for k in KS_G:
            rs = [r for r in H if r['k'] == k]
            fails = [r for r in rs if not r['bp']]
            frac_truth = '; '.join(f"{np.mean([r['truth'][tc] < 1 for r in rs]):.2f}" for tc in TAUS_H)
            med_truth = '; '.join(f"{np.median([r['truth'][tc] for r in rs]):.3g}" for tc in TAUS_H)
            at_bp = (f"{np.mean([r['bp_nnz'] for r in fails]):.1f}; " + '; '.join(f"{np.mean([r['at_bp'][tc] < 1 for r in fails]):.2f}" for tc in TAUS_H[1:])
                     if fails else '—')
            L_.append(f"| {k} | {len(rs)} | {frac_truth} | {med_truth} | {at_bp} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_colour_selection.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'colour_selection_raw.npz')
    md = os.path.join(OUT, 'colour_selection.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEFGH']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFGH')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(s, tc) for s in SEEDS_B for tc in TAUS_B], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, [(w_, k, i) for w_ in WIDTHS_D for k in KS_D for i in range(N_INST_D)], chunksize=1)
            h['Dn'] = pool.map_async(d_noisy_job, d_noisy_cases(), chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(N, d, k, i) for N, d in SIZES_C for k in KS_C[(N, d)] for i in range(N_INST_C[(N, d)])], chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(e_job, [(k, i) for k in KS_E for i in range(N_INST_E)], chunksize=1)
        if 'F' in parts:
            h['F'] = pool.map_async(f_job, [(k, i) for k in KS_E for i in range(N_INST_E)], chunksize=1)
        if 'G' in parts:
            h['G'] = pool.map_async(g_job, [(k, i) for k in KS_G for i in range(N_INST_G)], chunksize=1)
        if 'H' in parts:
            h['H'] = pool.map_async(h_job, [(k, i) for k in KS_G for i in range(N_INST_G)], chunksize=1)
        if 'A' in parts:
            res['A'] = part_A(); E.save(raw, **res); log('A')
        for k in 'BDCEFGH':
            if k in h:
                res[k] = h[k].get() if k != 'D' else dict(rates=h['D'].get(), noisy=h['Dn'].get())
                E.save(raw, **res); log(k)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEFGH']))


if __name__ == '__main__':
    main()
