"""
Proposition D‴.52: close features — the transit counts them at the number-detection limit and its block places them at the support limit
====================================================================================================================================

The sliding inverse scale space (Proposition D‴.48(iii); `rich_transit.sliding_iss` with the block's least-squares polish) on two Gaussian
bumps of the kernels' width σ = 0.12 a distance Δ apart (61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, equal masses).
With A = A₁ + A₂ ≠ 0 the pair is a monopole at the centroid c̄ = (A₁c₁ + A₂c₂)/A (kept at 0.032) with the moments V = A₁A₂Δ²/A,
T₃ = A₁A₂(A₁ − A₂)Δ³/A²: f* = Aφ(c̄) + ½V∂²_cφ(c̄) + ⅙T₃∂³_cφ(c̄) + O(Δ⁴); the first entrant fits Aφ(c̄) and leaves ½Vψ,
ψ = ∂²_cφ − (⟨∂²_cφ φ⟩/⟨φφ⟩)φ, whose correlations ½V⟨ψφ_j⟩ drive the second entry: Δτ₂ = min_j (s − u_j)/v_j over the walls s = ±1
the duals u_j(τ₁) = ⟨f*φ_j⟩/|⟨f*φ_{j₁}⟩| reach at the rates v_j = ½V⟨ψφ_j⟩ — ∝ Δ⁻²; a dipole (A = 0) has its first entry and the second
∝ Δ⁻¹.

A. The second entry against Δ (no noise): A = (0.8, 0.6), (0.8, 0.8), (0.8, −0.4), (0.8, −0.8): τ₁, the entrants, Δτ₂ against the
   leading-order prediction (for the dipole the transit of its derivative shape −A₁∂_cφ, slowed by 1/Δ), the scaled delays Δτ₂Δ² (Δτ₂Δ
   for the dipole), and the path's end (two kernels at the features?). The expansion parameter is |A₁ − A₂|Δ/(|A|σ).
B. The two-kernel block at its exact fit: the Gauss–Newton Hessian's eigenvalues and the fast flow's relaxation rates (the critical
   curve's flow linearised at the fit) against Δ — ∝ 1, 1, Δ², Δ⁶.
C. Number detection under label noise: y = f* + ξ, ξ i.i.d. N(0, s²) on the samples, A = (0.8, 0.6); under H1 the pair, under H0 the
   one-bump fit closest to it (the single-kernel least-squares fit of f*); R_C draws each; the second entry's delay Δτ₂ as the statistic:
   the power at the false-alarm rate 5% against Δ for each s, the separation Δ₅₀ where the power crosses ½, its slope in log s (½
   predicted: Δ² ∝ s).
D. Placement under label noise: the two-kernel block after the second entry on the H1 draws — the rms errors of the positions and of the
   amplitude ratio against s/Δ² and s/Δ³ (the moments' singular values Δ, Δ³).

Output: results/chaos/close_features.md and close_features_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D`).
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

OUT = os.environ.get('CF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('CF_WORKERS', '3'))
CBAR = 0.032
PAIRS_A = ((0.8, 0.6), (0.8, 0.8), (0.8, -0.4), (0.8, -0.8))
DELTAS_A = (0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625)
DELTAS_B = (0.2, 0.1, 0.05, 0.025, 0.0125)
PAIR_C, SIGMAS_C, R_C, ALPHA_C = (0.8, 0.6), (0.001, 0.003, 0.01, 0.03), 200, 0.05
DELTAS_C = (0.3, 0.2, 0.14, 0.1, 0.07, 0.05, 0.035, 0.025, 0.0175, 0.0125)


def setup():
    warnings.filterwarnings('ignore')
    L = X.law()
    K = RT.Bumps(np.linspace(-1, 1, X.NS_C), X.SIG_C)
    return L, K, np.linspace(-1, 1, X.D_C)


def pair_centres(A1, A2, D):
    """The pair's centres: the centroid at CBAR (A ≠ 0), the midpoint at CBAR (a dipole)."""
    A = A1 + A2
    if abs(A) < 1e-12:
        return np.array([CBAR - D / 2, CBAR + D / 2])
    return np.array([CBAR - A2 * D / A, CBAR + A1 * D / A])


def d2phi(K, c):
    c = np.atleast_1d(c)
    dx = K.xs[:, None] - c[None, :]
    return K.phi(c) * (dx ** 2 / K.sig ** 4 - 1 / K.sig ** 2)


def taylor_basis(K, c):
    """e_n = ∂ⁿ_cφ(c)/n!, n = 0 … 3, on the samples."""
    dx = K.xs - c
    ph, s2 = np.exp(-dx ** 2 / (2 * K.sig ** 2)), K.sig ** 2
    return [ph, ph * dx / s2, ph * (dx ** 2 / s2 ** 2 - 1 / s2) / 2, ph * (dx ** 3 / s2 ** 3 - 3 * dx / s2 ** 2) / 6]


def perp_norm2(vs, k):
    """⟨(P⊥e_k)²⟩, P⊥ the projection orthogonal to e_0 … e_{k−1} (sample means)."""
    B = np.array(vs[:k]).T
    e = vs[k]
    r = e - B @ np.linalg.lstsq(B, e, rcond=None)[0]
    return float(np.mean(r ** 2))


def hessian_prediction(K, A1, A2, D):
    """The block's two small Hessian eigenvalues to leading order: the split v = (1, −1, Δ/(2w₁), Δ/(2w₂)) leaves Jv = (Δ³/2)e₃,
    λ₁ = Δ⁶⟨(P⊥e₃)²⟩/8, amplitude-free; the separation (0, 0, 1/w₁, −1/w₂) leaves Jv = −2Δe₂, λ₂ = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²)."""
    cs = pair_centres(A1, A2, D)
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    return D ** 6 * perp_norm2(vs, 3) / 8, 4 * D ** 2 * perp_norm2(vs, 2) * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2)


def predicted_delay(K, f, c0, A1, A2, D):
    """The leading-order second entry: the duals at the first entry u_j = ⟨fφ_j⟩/|⟨fφ_{j₁}⟩| moving at v_j = ½V⟨ψφ_j⟩ (ψ at the centroid)
    to the wall ±1 they meet first; returns (Δτ₂, j₂, the wall's sign)."""
    nx = K.xs.size
    A = A1 + A2
    cb = (A1 * pair_centres(A1, A2, D)[0] + A2 * pair_centres(A1, A2, D)[1]) / A
    V = A1 * A2 * D ** 2 / A
    ph, p2 = K.phi(np.array([cb]))[:, 0], d2phi(K, cb)[:, 0]
    psi = p2 - (p2 @ ph) / (ph @ ph) * ph
    g = K.phi(c0).T @ f / nx
    j1 = int(np.argmax(np.abs(g)))
    u = g / abs(g[j1])
    v = 0.5 * V * (K.phi(c0).T @ psi) / nx                            # the duals' rates after the first fit (ρ = 1, per unit τ)
    best = (np.inf, None, 0)
    for j in range(c0.size):
        if j == j1 or v[j] == 0:
            continue
        wall = 1.0 if v[j] > 0 else -1.0
        t = (wall - u[j]) / v[j]
        if 0 < t < best[0]:
            best = (float(t), j, int(wall))
    return best


def a_job(args):
    A1, A2, D = args
    L, K, c0 = setup()
    cs = pair_centres(A1, A2, D)
    f = X.bumps_target(K, cs, (A1, A2))
    t0 = time.time()
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    ent = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
    act = np.flatnonzero(R['act'])
    out = dict(A=(A1, A2), D=D, cs=cs.tolist(), entries=ent, n_act=int(act.size), c_end=np.sort(R['c'][act]).tolist(),
               w_end=R['w'][act][np.argsort(R['c'][act])].tolist(), rms=float(R['hist'][-1][3]), relaxed=all(R['relaxed']),
               secs=time.time() - t0)
    if abs(A1 + A2) > 1e-12:
        out['pred'] = predicted_delay(K, f, c0, A1, A2, D)
    else:                                                              # a dipole: the transit of its derivative shape, slowed by 1/Δ
        shape = -A1 * K.dphi(np.array([CBAR]))[:, 0]
        Rs = RT.sliding_iss(L, K, shape, c0, np.ones(c0.size), max_events=2, polish=True)
        es = [(float(t), int(j)) for t, j, k, c in Rs['events'] if k == 'enter']
        out['shape'] = dict(tau1=es[0][0] / D, dtau2=(es[1][0] - es[0][0]) / D, kernels=[j for t, j in es])
    return out


def part_A():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, [(a1, a2, D) for a1, a2 in PAIRS_A for D in DELTAS_A])


def b_job(args):
    A1, A2, D = args
    L, K, c0 = setup()
    nx = K.xs.size
    cs = pair_centres(A1, A2, D)
    w = np.array([A1, A2])
    J = np.c_[K.phi(cs), K.dphi(cs) * w[None, :]] / np.sqrt(nx)
    hess = np.sort(np.linalg.eigvalsh(J.T @ J))
    f = X.bumps_target(K, cs, (A1, A2))
    rhs = RT.critical_rhs(L, K, f)
    al = np.sign(w) * RT._alpha_of(L, np.abs(w), RT.coef_crit)
    z = np.r_[al, cs]
    Jf = np.zeros((4, 4))
    for k in range(4):
        h = 1e-6 * max(1.0, abs(z[k]))
        e = np.zeros(4); e[k] = h
        Jf[:, k] = (rhs(0, z + e) - rhs(0, z - e)) / (2 * h)
    ev, V = np.linalg.eig(Jf)
    rates = np.sort(-np.real(ev))
    l1, l2 = hessian_prediction(K, A1, A2, D)
    # the flow's split and separation modes, identified by their eigenvectors in z = (α₁, α₂, c₁, c₂), against their leading order:
    # ẇ = −μ∂E/∂w (μ = w′(α)·curve_mobility, the amplitudes' mobility), ċ = −Gm∂E/∂c; the split's rate λ₁·2/(1/μ₁ + 1/μ₂) (its
    # Rayleigh quotient in the mobility's metric, the centres' share O(Δ²)); the separation's at frozen amplitudes (Gm Δ² ≫ μ)
    # Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) — the projection orthogonal to the translation only.
    h_ = 1e-7
    wp = np.array([(RT.coef_crit(L, abs(a) + h_) - RT.coef_crit(L, abs(a) - h_)) / (2 * h_) for a in al])
    mu = wp * RT.curve_mobility(L, al)
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    B1 = vs[1][:, None]
    e2p = vs[2] - B1 @ np.linalg.lstsq(B1, vs[2], rcond=None)[0]
    pred_split = l1 * 2 / (1 / mu[0] + 1 / mu[1])
    pred_sep = L.G * L.m * 4 * D ** 2 * float(np.mean(e2p ** 2)) * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2)

    def mode(s_):
        Vr = np.real(V)
        cos = np.abs(Vr.T @ s_) / (np.linalg.norm(Vr, axis=0) * np.linalg.norm(s_))
        k = int(np.argmax(cos))
        return float(-np.real(ev[k])), float(cos[k])
    split = mode(np.r_[1 / wp[0], -1 / wp[1], D / (2 * A1), D / (2 * A2)])
    sep = mode(np.r_[0.0, 0.0, 1 / A1, -1 / A2])
    return dict(A=(A1, A2), D=D, hess=hess.tolist(), rates=rates.tolist(), pred=(l1, l2), mu=mu.tolist(),
                split=split, sep=sep, pred_flow=(pred_split, pred_sep))


def part_B():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(a1, a2, D) for a1, a2 in PAIRS_A for D in DELTAS_B])


def one_fit(K, f):
    """The single-kernel least-squares fit of f (the closest one-bump target): returns its centre and amplitude."""
    nx = K.xs.size
    g = K.phi(np.linspace(-1, 1, 201)).T @ f / nx
    c = np.linspace(-1, 1, 201)[int(np.argmax(np.abs(g)))]
    w = 0.0
    for _ in range(100):
        ph = K.phi(np.array([c]))[:, 0]
        w = (ph @ f) / (ph @ ph)
        dp = K.dphi(np.array([c]))[:, 0]
        r = f - w * ph
        gr = w * (dp @ r)
        hs = w * w * (dp @ dp)
        c = c + gr / hs
    return float(c), float(w)


def c_job(args):
    D, s, hyp, seed = args
    L, K, c0 = setup()
    cs = pair_centres(*PAIR_C, D)
    f1 = X.bumps_target(K, cs, PAIR_C)
    if hyp == 1:
        f = f1
    else:
        c1, w1 = one_fit(K, f1)
        f = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    out = []
    for r in range(R_C):
        y = f + s * rng.standard_normal(K.xs.size)
        R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True)
        ent = [(t, j) for t, j, k, c in R['events'] if k == 'enter']
        dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.inf
        row = dict(dt=float(dt), j=[j for t, j in ent])
        if hyp == 1 and len(ent) > 1:
            act = np.flatnonzero(R['act'])                          # the block after the second entry's relaxation
            if act.size == 2:
                o = np.argsort(R['c'][act])
                row.update(c=R['c'][act][o].tolist(), w=R['w'][act][o].tolist())
        out.append(row)
    return dict(D=D, s=s, hyp=hyp, rows=out)


def part_C():
    from multiprocessing import Pool
    jobs = [(D, s, h, 1000 * k + 10 * h + i) for i, (D, s) in enumerate([(D, s) for D in DELTAS_C for s in SIGMAS_C])
            for h, k in ((0, 1), (1, 2))]
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, jobs)


def detection(C):
    """Per (Δ, s): the H0 quantile of Δτ₂ at ALPHA_C and the H1 power; Δ₅₀(s) by log-linear interpolation of the power."""
    tab = {}
    for r in C:
        tab.setdefault((r['D'], r['s']), {})[r['hyp']] = np.array([q['dt'] for q in r['rows']])
    out = {}
    for (D, s), d in tab.items():
        q = np.quantile(d[0], ALPHA_C)
        out[(D, s)] = dict(q0=float(q), power=float(np.mean(d[1] < q)), med1=float(np.median(d[1])), med0=float(np.median(d[0])))
    d50 = {}
    for s in SIGMAS_C:
        Ds = sorted(DELTAS_C)
        pw = [out[(D, s)]['power'] for D in Ds]
        d50[s] = np.nan
        for k in range(len(Ds) - 1):
            if pw[k] < 0.5 <= pw[k + 1]:
                t = (0.5 - pw[k]) / (pw[k + 1] - pw[k])
                d50[s] = float(np.exp(np.log(Ds[k]) + t * (np.log(Ds[k + 1]) - np.log(Ds[k]))))
                break
    return out, d50


def report(A=None, B=None, C=None):
    fmt = E.fmt
    L_ = ['# Proposition D‴.52: close features (`exp_close_features.py`)', '',
          f"The sliding inverse scale space (polished) on two bumps of the kernels' width {X.SIG_C} a distance Δ apart, the centroid (the "
          f"midpoint for a dipole) at {CBAR}; 61 samples of [−1, 1], 21 kernels on a grid of spacing 0.1, equal masses.", '']
    if A is not None:
        L_ += ['## A. The second entry', '',
               '| A₁, A₂ | Δ | entries (τ, kernel) | Δτ₂ | predicted (kernel, wall) | ratio | Δτ₂Δ² | Δτ₂Δ | the end: kernels, centres (truth), rms |',
               '|---|---|---|---|---|---|---|---|---|']
        for r in A:
            ent = r['entries']
            dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.nan
            p, sh = r.get('pred'), r.get('shape')
            if p:
                pr, ratio = f"{p[0]:.4g} ({p[1]}, {p[2]:+d})", f"{dt / p[0]:.4f}"
            elif sh:
                pr = f"shape: τ₁ {sh['tau1']:.4g}, Δτ₂ {sh['dtau2']:.4g} ({sh['kernels']})"
                ratio = f"{ent[0][0] / sh['tau1']:.4f}, {dt / sh['dtau2']:.4f}"
            else:
                pr = ratio = '—'
            L_.append(f"| {r['A'][0]:g}, {r['A'][1]:g} | {r['D']:g} | {', '.join(f'({t:.4g}, {j})' for t, j in ent[:3])} | {dt:.4g} | {pr} | "
                      f"{ratio} | {dt * r['D'] ** 2:.4g} | {dt * r['D']:.4g} | {r['n_act']}, {[round(v, 5) for v in r['c_end']]} "
                      f"({[round(v, 5) for v in r['cs']]}), {fmt(r['rms'])} |")
        L_.append('')
    if B is not None:
        L_ += ['## B. The two-kernel block at its exact fit', '',
               'The Gauss–Newton Hessian JᵀJ (J = ∂(Φ(c)w)/∂(w, c)/N^{1/2}) and the fast flow on the critical curve linearised at the fit '
               '(its relaxation rates): eigenvalues, the two smallest scaled by Δ² and Δ⁶ against their leading order (the split mode: '
               'λ₁ = Δ⁶⟨(P⊥e₃)²⟩/8, amplitude-free; the separation: λ₂ = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²); e_n = ∂ⁿ_cφ/n! at the pair, P⊥ '
               'orthogonal to the lower e_m), and the local slopes of the smallest two in log Δ. The flow\'s split and separation modes are '
               'identified by their eigenvectors\' overlap with (1/w′₁, −1/w′₂, Δ/(2w₁), Δ/(2w₂)) and (0, 0, 1/w₁, −1/w₂) in (α, c), and '
               'compared with λ₁·2/(1/μ₁ + 1/μ₂) (μ_k the amplitudes\' mobility, ẇ = −μ∂E/∂w) and, at frozen amplitudes, '
               'Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²).', '',
               '| A₁, A₂ | Δ | Hessian eigenvalues | λ₂/Δ², λ₁/Δ⁶ | predicted | slopes | flow rates | split: r/Δ⁶ (predicted; '
               'the mode\'s overlap) | separation: r/Δ² (predicted; overlap) | slopes |',
               '|---|---|---|---|---|---|---|---|---|---|']
        prev = {}
        for r in B:
            key = tuple(r['A'])
            h, q = r['hess'], r['rates']
            sl = sq = '—'
            sp, se = r.get('split'), r.get('sep')
            fq = (sp[0], se[0]) if sp else (q[0], q[1])
            if key in prev:
                pD, ph, pq = prev[key]
                sl = ', '.join(f"{np.log(h[k] / ph[k]) / np.log(r['D'] / pD):.3f}" for k in (0, 1))
                sq = ', '.join(f"{np.log(fq[k] / pq[k]) / np.log(r['D'] / pD):.3f}" for k in (0, 1))
            prev[key] = (r['D'], h, fq)
            pl, pf = r.get('pred'), r.get('pred_flow')
            pt = f"{pl[1] / r['D'] ** 2:.4g}, {pl[0] / r['D'] ** 6:.4g}" if pl else '—'
            D6, D2 = r['D'] ** 6, r['D'] ** 2
            if sp:
                fs = (f"{sp[0] / D6:.4g} ({pf[0] / D6:.4g}; {sp[1]:.3f}) | {se[0] / D2:.4g} ({pf[1] / D2:.4g}; {se[1]:.3f})")
            else:
                fs = f"{q[0] / D6:.4g} | {q[1] / D2:.4g}"
            L_.append(f"| {key[0]:g}, {key[1]:g} | {r['D']:g} | {', '.join(f'{v:.3e}' for v in h)} | {h[1] / D2:.4g}, "
                      f"{h[0] / D6:.4g} | {pt} | {sl} | {', '.join(f'{v:.3e}' for v in q)} | {fs} | {sq} |")
        L_.append('')
    if C is not None:
        det, d50 = detection(C)
        L_ += ['## C. Number detection under label noise', '',
               f"A = {PAIR_C}; {R_C} draws per hypothesis; H0 the closest one-bump target; the power of Δτ₂ < its H0 quantile "
               f"{ALPHA_C:g} (the medians of Δτ₂ under H1 and H0):", '',
               '| Δ | ' + ' | '.join(f's = {s:g}' for s in SIGMAS_C) + ' |', '|---|' + '---|' * len(SIGMAS_C)]
        for D in DELTAS_C:
            L_.append(f"| {D:g} | " + ' | '.join(f"{det[(D, s)]['power']:.2f} ({det[(D, s)]['med1']:.3g}, {det[(D, s)]['med0']:.3g})"
                                                for s in SIGMAS_C) + ' |')
        ss = np.array([s for s in SIGMAS_C if np.isfinite(d50[s])])
        dd = np.array([d50[s] for s in ss])
        sl = np.polyfit(np.log(ss), np.log(dd), 1)[0] if ss.size >= 2 else np.nan
        # the prediction: under H1 the second entry comes after ≈ C₁/Δ² (C₁ the leading order's Δτ₂Δ² at the smallest Δ, (i)), so the
        # power is ½ where C₁/Δ² meets the H0 quantile q_α (∝ √N/s): Δ₅₀ = (C₁/q_α)^{1/2}, q_α interpolated in log Δ at the crossing
        Lc, Kc, c0c = setup()
        Dm = min(DELTAS_C)
        C1 = predicted_delay(Kc, X.bumps_target(Kc, pair_centres(*PAIR_C, Dm), PAIR_C), c0c, *PAIR_C, Dm)[0] * Dm ** 2
        pred50 = {}
        for s_ in SIGMAS_C:
            Ds = np.array(sorted(DELTAS_C))
            gap = np.log(C1 / Ds ** 2) - np.log([det[(D, s_)]['q0'] for D in Ds])
            pred50[s_] = np.nan
            for k in range(len(Ds) - 1):
                if gap[k] > 0 >= gap[k + 1] or gap[k] <= 0 < gap[k + 1]:
                    t = gap[k] / (gap[k] - gap[k + 1])
                    pred50[s_] = float(np.exp(np.log(Ds[k]) + t * (np.log(Ds[k + 1]) - np.log(Ds[k]))))
                    break
        qs = {s_: float(np.median([det[(D, s_)]['q0'] for D in DELTAS_C]) * s_) for s_ in SIGMAS_C}
        L_ += ['', 'Δ₅₀ (the power ½): ' + ', '.join(f's = {s:g}: {fmt(d50[s])}' for s in SIGMAS_C) +
               f"; the slope of log Δ₅₀ in log s {sl:.3f} (½ predicted). Predicted from (i) and the H0 quantile, Δ₅₀ = (C₁/q_α)^{{1/2}} "
               f"(C₁ = {C1:.4f}, the leading order's Δτ₂Δ² at Δ = {Dm:g}): " + ', '.join(f's = {s:g}: {fmt(pred50[s])}' for s in SIGMAS_C) +
               "; the H0 quantile times s (its median over Δ; constant if q_α ∝ 1/s): " + ', '.join(f'{qs[s]:.4g}' for s in SIGMAS_C) + '.', '']
        L_ += ['## D. Placement under label noise (the two-kernel block after the second entry, H1)', '',
               '| Δ | s | runs | rms position error | × Δ²/s | rms error of w₁/(w₁ + w₂) | × Δ³/s |', '|---|---|---|---|---|---|---|']
        cs_truth = {D: pair_centres(*PAIR_C, D) for D in DELTAS_C}
        for r in C:
            if r['hyp'] != 1:
                continue
            ok = [q for q in r['rows'] if 'c' in q]
            if not ok:
                continue
            pe = np.sqrt(np.mean([np.mean((np.array(q['c']) - cs_truth[r['D']]) ** 2) for q in ok]))
            fr = np.sqrt(np.mean([(q['w'][0] / (q['w'][0] + q['w'][1]) - PAIR_C[0] / sum(PAIR_C)) ** 2 for q in ok]))
            L_.append(f"| {r['D']:g} | {r['s']:g} | {len(ok)} | {pe:.3e} | {pe * r['D'] ** 2 / r['s']:.3g} | {fr:.3e} | "
                      f"{fr * r['D'] ** 3 / r['s']:.3g} |")
        L_.append('')
    L_.append('Numbers are produced by `experiments/exp_close_features.py`; nothing is edited by hand.')
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'close_features_raw.npz')
    md = os.path.join(OUT, 'close_features.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else {'A', 'B', 'C', 'D'}
    res = dict(old)
    t0 = time.time()
    for k, fn in (('B', part_B), ('A', part_A), ('C', part_C)):
        if k in parts or (k == 'C' and 'D' in parts):
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
