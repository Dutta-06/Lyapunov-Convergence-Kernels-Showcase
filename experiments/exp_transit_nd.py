"""
Proposition D‴.54: the transit in ℝⁿ — one feature on a lattice of kernels (the k-fold near tie and its shares, the formation map at a
vertex) and several features in the plane
====================================================================================================================================

The free law with moving centres (Proposition D‴.48(iii); m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025), Gaussian kernels of width SIG in
the plane (`rich_transit_nd`, the samples on a tensor grid), masses ε, from α = 0.

A. The k-fold near tie. One feature of amplitude A near a tie of the kernels' square lattice of spacing H: an edge (the kernels at (0, 0)
   and (H, 0), the feature near (H/2, 0): a 2-fold tie) and a vertex ((0, 0), (H, 0), (0, H), (H, H), the feature near (H/2, H/2): a
   4-fold tie), placed so that the runner-ups' margins when the entrant enters are μ_j = ε^{1/2}x_j for the given gaps (the diagonal's
   gap then follows: x₁₁ ≈ x₁₀ + x₀₁). The loop to TAU_A τ_entry at ε = EPS_A; the shares |w_j|/Σ|w| against the ε-free inner problem
   with the k kernels at the same gaps (the feature at the tie; `inner_k`), the softmax e^{−κx_j}/Σe^{−κx_i} (κ = πρ(2c)^{1/2}) and, at
   the vertex, the product of the two edges' inner problems (the factorisation of the 4-fold split).
B. The formation map at a vertex: the number of kernels formed (|w| > TH_REL·|A|) over gaps (x₁₀, x₀₁) on a grid, the loop at EPS_B
   against the 4-kernel inner problem and the independent rule #{j: x_j < x_θ} (x_θ: the 2-fold inner problem's crossing of TH_REL);
   the probabilities for a feature placed uniformly in a Voronoi cell, P(N ≥ 2) ∝ ε^{1/2}, P(N ≥ 3) and P(N = 4) ∝ ε, from the map's areas.
C. Several features in the plane: random three-feature targets on the lattice over [−1, 1]², the sliding path against the loop at ε =
   EPS_C: the kernels formed, the never-entering kernels' largest coefficient (∝ ε^{3/2}) and displacement (∝ ε), the block against the
   path's.
T. The 2-kernel (edge) inner problem's runner-up share p(x) on a grid of gaps — the product rule: at a vertex the kernel e + h1_S has
   the share Π_{i∈S}p(x_i)Π_{i∉S}(1 − p(x_i)), formed when it exceeds TH_REL —, used by B's map and D.
D. The probabilities: features placed uniformly in a Voronoi cell (N_D positions) at ε = EPS_D, the loop's count against the product
   rule at the positions' gaps and against the asymptotic P(N ≥ 2) = 4w/h − 4w²/h², w = 2σ²x_θε^{1/2}/h the layer's width.

Usage: python experiments/exp_transit_nd.py [--only A,B,C,T,D] [--from-raw] → results/chaos/transit_nd.md (+ _raw.npz).
"""
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
from src.core import rich_transit_nd as RN
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('TN_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TND_WORKERS', '2'))
SIG, H, DS = 0.15, 0.125, 0.05                                        # the kernels' width, the lattice's spacing, the samples' spacing
TH_REL = 1e-2                                                          # formed: |w| > TH_REL·|A|
TAU_A = 1.32                                                           # the loop to 1.32 τ_entry (the relaxation after the entry)
AMPS_A = (0.3, 1.0)
EPS_A = (1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
EDGE_X = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0)
VERTEX_X = ((0.0, 0.0), (0.5, 0.5), (1.0, 1.0), (1.0, 3.0), (2.0, 2.0), (2.0, 5.0), (4.0, 4.0))
A_B, EPS_B, GRID_B = 0.8, (1e-6, 1e-8, 1e-10), np.arange(0.0, 10.5, 1.0)
SEED_C, N_C, EPS_C = 54, 6, (1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
X_T = np.round(np.arange(0.0, 14.05, 0.1), 10)                         # the 2-kernel (edge) inner problem's share, tabulated
SEED_D, N_D, EPS_D = 541, 600, (1e-5, 1e-6)


def law():
    warnings.filterwarnings('ignore')
    return X.law()


def local():
    """One feature: the samples symmetric about the vertex (H/2, H/2) (so that the ties are exact), the kernels on the lattice H·(i, j),
    i, j = −3 … 4, around it."""
    xs = H / 2 + DS * np.arange(-18, 19)
    K = RN.BumpsGrid2(xs, xs, SIG)
    g = H * np.arange(-3, 5)
    return K, np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)


def kidx(c0, p):
    return int(np.argmin(((c0 - np.asarray(p, float)) ** 2).sum(1)))


TIES = {'edge': ([(0, 0), (H, 0)], np.array([H / 2, 0.0])),
        'vertex': ([(0, 0), (H, 0), (0, H), (H, H)], np.array([H / 2, H / 2]))}


def margins(K, c0, A, cstar, idx):
    """The margins μ_j = 1 − v_j/v_e (v = |⟨fφ_j⟩| at the start, e the largest) of the kernels idx (the first the entrant)."""
    v = np.abs(K.corr(c0[idx], A * K.phi(np.atleast_2d(cstar))[:, 0])[0])
    return 1 - v / v[0], v[0]


def place(K, c0, A, kind, x, eps):
    """The feature's centre near the tie with the runner-ups' margins ε^{1/2}x (edge: x₁₀; vertex: x₁₀, x₀₁ — the diagonal's follows), by
    Newton's method on the offset from the tie toward the entrant."""
    pts, T = TIES[kind]
    idx = [kidx(c0, p) for p in pts]
    m = 1 if kind == 'edge' else 2
    tgt = np.sqrt(eps) * np.asarray(x, float)[:m]
    d = np.full(m, 1e-12) + tgt * 2 * SIG ** 2 / H                    # the continuum's linear map as the start
    for _ in range(30):
        def mu(dd):
            off = np.zeros(2); off[:m] = dd
            return margins(K, c0, A, T - off, idx)[0][1:1 + m]
        f0 = mu(d) - tgt
        J = np.empty((m, m))
        for i in range(m):
            e = np.zeros(m); e[i] = 1e-8
            J[:, i] = (mu(d + e) - mu(d - e)) / (2 * e[i])
        step = np.linalg.solve(J, f0)
        d = d - step
        if np.all(np.abs(step) <= 1e-15 + 1e-12 * np.abs(d)):
            break
    off = np.zeros(2); off[:m] = d
    return T - off, idx


def softmax(x, L):
    kap = np.pi * (2 * L.c) ** 0.5
    s = np.exp(-kap * (np.asarray(x, float) - np.min(x)))
    return s / s.sum()


def a_job(args):
    """One loop: (kind, A, x, ε) → the shares of the tie's kernels, their margins in units of ε^{1/2}, the others' largest |w|."""
    kind, A, x, eps = args
    t0 = time.time()
    L = law(); K, c0 = local()
    cstar, idx = place(K, c0, A, kind, x, eps)
    mu, v_e = margins(K, c0, A, cstar, idx)
    f = A * K.phi(np.atleast_2d(cstar))[:, 0]
    taus = np.linspace(0, TAU_A / v_e, 23)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus)
    w = np.abs(W[-1][idx])
    rest = np.delete(np.abs(W[-1]), idx)
    return dict(kind=kind, A=A, x=list(x), eps=eps, cstar=cstar, x_act=mu / np.sqrt(eps), shares=w / w.sum(), w=W[-1][idx],
                c=C[-1][idx], rest=float(rest.max()), n_formed=int((np.abs(W[-1]) > TH_REL * abs(A)).sum()), secs=time.time() - t0)


def inner_job(args):
    """The ε-free inner problem at the tie: (kind, A, gaps of the tie's kernels) → shares."""
    kind, A, xs = args
    L = law(); K, c0 = local()
    pts, T = TIES[kind]
    f = A * K.phi(np.atleast_2d(T))[:, 0]
    out = RN.inner_k(L, K, f, np.array(pts, float), np.asarray(xs, float))
    return dict(kind=kind, A=A, xs=list(xs), shares=out['shares'], stationary=out['stationary'], s=out['s'], c=out['c'])


def part_A():
    from multiprocessing import Pool
    loops = [('edge', A, (x,), e) for A in AMPS_A for x in EDGE_X for e in EPS_A]
    loops += [('vertex', A, x, e) for A in AMPS_A for x in VERTEX_X for e in EPS_A]
    with Pool(N_WORKERS) as pool:
        rl = pool.map(a_job, loops, chunksize=1)
        # the inner problem at the loop's gaps at the smallest ε (the limit's), and the edges' problems for the factorisation
        jobs = []
        for r in rl:
            if r['eps'] == min(EPS_A):
                jobs.append((r['kind'], r['A'], tuple(np.r_[0.0, r['x_act'][1:]])))
                if r['kind'] == 'vertex':
                    jobs += [('edge', r['A'], (0.0, r['x_act'][1])), ('edge', r['A'], (0.0, r['x_act'][2]))]
        ri = pool.map(inner_job, jobs, chunksize=1)
    return dict(loops=rl, inner=ri)


def threshold_x(L, K, A, lo=0.0, hi=14.0, it=30):
    """x_θ: the gap at which the 2-fold (edge) inner problem's runner-up share crosses TH_REL."""
    for _ in range(it):
        m_ = 0.5 * (lo + hi)
        if inner_job(('edge', A, (0.0, m_)))['shares'][1] > TH_REL:
            lo = m_
        else:
            hi = m_
    return 0.5 * (lo + hi)


def b_job(args):
    x1, x2, eps = args
    r = a_job(('vertex', A_B, (x1, x2), eps))
    return dict(x=(x1, x2), eps=eps, n=r['n_formed'], x_act=r['x_act'], shares=r['shares'])


def b_inner(args):
    x1, x2 = args
    L = law(); K, c0 = local()
    xd = margins_diag(K, c0, x1, x2)
    r = inner_job(('vertex', A_B, (0.0, x1, x2, xd)))
    return dict(x=(x1, x2), xd=xd, n=int((r['shares'] > TH_REL).sum()), shares=r['shares'], stationary=r['stationary'])


def margins_diag(K, c0, x1, x2, eps=1e-10):
    """The diagonal's gap at the vertex for the edge gaps (x₁₀, x₀₁) (at a small ε, where the map is linear)."""
    cstar, idx = place(K, c0, A_B, 'vertex', (x1, x2), eps)
    return float(margins(K, c0, A_B, cstar, idx)[0][3] / np.sqrt(eps))


def part_B():
    from multiprocessing import Pool
    pts = [(a, b) for a in GRID_B for b in GRID_B if a <= b]
    with Pool(N_WORKERS) as pool:
        rl = pool.map(b_job, [(a, b, e) for e in EPS_B for a, b in pts], chunksize=1)
        ri = pool.map(b_inner, pts, chunksize=1)
    L = law(); K, c0 = local()
    return dict(loops=rl, inner=ri, x_theta=threshold_x(L, K, A_B))


def targets_C(seed, n):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        cs = rng.uniform(-0.7, 0.7, (3, 2))
        if min(np.linalg.norm(cs[i] - cs[j]) for i in range(3) for j in range(i)) < 3 * SIG:
            continue
        out.append((cs, rng.choice([-1.0, 1.0], 3) * rng.uniform(0.4, 1.0, 3)))
    return out


def full():
    xs = np.arange(-1, 1 + 1e-9, DS)
    K = RN.BumpsGrid2(xs, xs, SIG)
    g = np.arange(-1, 1 + 1e-9, H)
    return K, np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)


def c_job(args):
    q, eps = args
    t0 = time.time()
    L = law(); K, c0 = full()
    cs, As = targets_C(SEED_C, N_C)[q]
    f = K.synth(cs, As)
    R = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)))
    path = sorted(int(j) for j in np.flatnonzero(R['act']))
    t_last = R['events'][-1][0]
    taus = np.linspace(0, 1.3 * t_last, 27)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus)
    big = np.abs(W[-1]) > TH_REL * np.abs(As).max()
    never = np.array([j for j in range(len(c0)) if j not in path])
    rms = float(np.sqrt(np.mean((f - K.synth(C[-1], W[-1])) ** 2)))
    return dict(q=q, eps=eps, path=path, events=[(float(t), int(j), k) for t, j, k, c in R['events']], path_rms=R['hist'][-1][3],
                formed=sorted(int(j) for j in np.flatnonzero(big)), w_never=float(np.abs(W[-1][never]).max()),
                disp_never=float(np.linalg.norm(C[-1][never] - c0[never], axis=1).max()),
                dw=float(np.abs(W[-1][path] - R['w'][path]).max()), dc=float(np.abs(C[-1][path] - R['c'][path]).max()), rms=rms,
                secs=time.time() - t0)


def part_C():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, [(q, e) for q in range(N_C) for e in EPS_C], chunksize=1)


def t_job(x):
    return float(inner_job(('edge', A_B, (0.0, float(x))))['shares'][1])


def part_T():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return dict(x=X_T, p=np.array(pool.map(t_job, list(X_T), chunksize=4)))


def product_count(T, x1, x2):
    """The product rule's count at a vertex with the edge gaps x1, x2 (p interpolated in the table T)."""
    p1, p2 = np.interp([x1, x2], T['x'], T['p'], right=0.0)
    shares = [(1 - p1) * (1 - p2), p1 * (1 - p2), (1 - p1) * p2, p1 * p2]
    return int(sum(s > TH_REL for s in shares))


def d_job(args):
    """A feature at a uniform position in the Voronoi cell of (0, 0): the loop's count and the gaps to the two nearest faces."""
    q, eps = args
    cstar = np.random.default_rng(SEED_D + q).uniform(-H / 2, H / 2, 2)
    L = law(); K, c0 = local()
    f = A_B * K.phi(np.atleast_2d(cstar))[:, 0]
    v = np.abs(K.corr(c0, f)[0])
    e = int(np.argmax(v))
    taus = np.linspace(0, TAU_A / v[e], 23)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus)
    gaps = []
    for ax in range(2):                                                # the nearer face along each axis
        nb = [kidx(c0, c0[e] + s * H * np.eye(2)[ax]) for s in (-1, 1)]
        gaps.append(min((1 - v[j] / v[e]) / np.sqrt(eps) for j in nb))
    return dict(q=q, eps=eps, cstar=cstar, e=e, n=int((np.abs(W[-1]) > TH_REL * A_B).sum()), gaps=gaps)


def part_D():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(d_job, [(q, e) for e in EPS_D for q in range(N_D)], chunksize=4)


def slope(x, y):
    return float(np.polyfit(np.log10(x), np.log10(y), 1)[0])


def report(A=None, B=None, C=None, T=None, D=None):
    L = law()
    kap = np.pi * (2 * L.c) ** 0.5
    L_ = ["# Proposition D‴.54: the transit in ℝⁿ (`exp_transit_nd.py`)", "",
          f"The free law with moving centres (m = {L.m:g}, G = {L.G:g}, γ = {L.gamma:g}, a_max = {L.A:g}, c = {L.c:g}), Gaussian kernels of width "
          f"{SIG:g} in the plane on a square lattice of spacing {H:g}, the samples on a grid of spacing {DS:g}; masses ε; formed: |w| > "
          f"{TH_REL:g}|A|; κ = πρ(2c)^{{1/2}} = {kap:.4f}.", ""]
    if A is not None:
        L_ += ["## A. The k-fold near tie", "",
               f"One feature near an edge (the kernels at (0, 0), ({H:g}, 0)) or a vertex (and ({H:g}, 0), (0, {H:g}), ({H:g}, {H:g})) of the "
               f"lattice — 8 × 8 kernels around it, the samples symmetric about the vertex —, its margins ε^{{1/2}}x at the entry; the loop to "
               f"{TAU_A:g} τ_entry. The gaps x (the runner-ups', as placed; at the vertex the diagonal's as it follows), the shares of the "
               "tie's kernels in the loop against the inner problem's at the gaps of the smallest ε (the feature at the tie), the softmax at "
               "those gaps and, at the vertex, the product of the two edges' inner problems:", "",
               "| tie | A | x | shares: loop (ε = " + ', '.join(f'{e:g}' for e in EPS_A) + ") | inner problem | softmax | the edges' product | "
               "the others' largest \\|w\\| (ε = " + ', '.join(f'{e:g}' for e in EPS_A) + ") |", "|---|---|---|---|---|---|---|---|"]
        inner = {(r['kind'], r['A'], tuple(np.round(r['xs'], 9))): r for r in A['inner']}
        groups = {}
        for r in A['loops']:
            groups.setdefault((r['kind'], r['A'], tuple(r['x'])), []).append(r)
        worst = {}
        for key, rs in groups.items():
            rs = sorted(rs, key=lambda r: -r['eps'])
            kind, Aa, x = key
            xa = np.r_[0.0, rs[-1]['x_act'][1:]]
            ri = inner.get((kind, Aa, tuple(np.round(xa, 9))))
            sm = softmax(xa, L)
            prod = '—'
            if kind == 'vertex':
                e1 = inner.get(('edge', Aa, tuple(np.round([0.0, xa[1]], 9))))
                e2 = inner.get(('edge', Aa, tuple(np.round([0.0, xa[2]], 9))))
                if e1 is not None and e2 is not None:
                    p1, p2 = e1['shares'][1], e2['shares'][1]
                    pv = np.array([(1 - p1) * (1 - p2), p1 * (1 - p2), (1 - p1) * p2, p1 * p2])
                    prod = ', '.join(f'{v:.4g}' for v in pv)
                    if ri is not None:
                        worst.setdefault(('prod', kind), []).append(float(np.abs(pv - ri['shares']).max()))
            if ri is not None:
                worst.setdefault(('loop', kind), []).append(float(np.abs(rs[-1]['shares'] - ri['shares']).max()))
                worst.setdefault(('soft', kind), []).append(float(np.abs(sm - ri['shares']).max()))
            L_.append(f"| {kind} | {Aa:g} | {', '.join(f'{v:.3f}' for v in xa[1:])} | "
                      + '; '.join(', '.join(f'{v:.4g}' for v in r['shares']) for r in rs) + " | "
                      + (', '.join(f'{v:.4g}' for v in ri['shares']) + ('' if ri['stationary'] else ' (not stationary)') if ri is not None else '—')
                      + f" | {', '.join(f'{v:.4g}' for v in sm)} | {prod} | " + ', '.join(f"{r['rest']:.2e}" for r in rs) + " |")
        L_ += ["", "The largest difference of the shares from the inner problem's (at the smallest ε): " + '; '.join(
            f"{k[1]} — {dict(loop='the loop', soft='the softmax', prod='the edges’ product')[k[0]]} {max(v):.2e}" for k, v in sorted(worst.items())) + ".", ""]
    if B is not None:
        xt = B['x_theta']
        L_ += ["## B. The formation map at a vertex", "",
               f"A = {A_B:g}; the gaps (x₁₀, x₀₁) of the edge runner-ups on a grid (the diagonal's ≈ their sum), x₁₀ ≤ x₀₁; the number of kernels "
               f"formed by the loop at ε = {', '.join(f'{e:g}' for e in EPS_B)}, by the 4-kernel inner problem, and by the independent rule "
               f"#{{j: x_j < x_θ}} + 1 with x_θ = {xt:.3f} (the 2-fold inner problem's share crossing {TH_REL:g}; the Gibbs split's "
               f"ln(1/{TH_REL:g} − 1)/κ = {np.log(1 / TH_REL - 1) / kap:.3f}):", "",
               "| x₁₀ | x₀₁ | x₁₁ | loop | inner problem | the independent rule | the product rule (T) |", "|---|---|---|---|---|---|---|"]
        byx = {}
        for r in B['loops']:
            byx.setdefault(tuple(r['x']), {})[r['eps']] = r['n']
        inn = {tuple(r['x']): r for r in B['inner']}
        agree = {e: 0 for e in EPS_B}; agree_rule = 0; agree_prod = 0; tot = 0
        for x in sorted(byx):
            ri = inn[x]
            rule = 1 + int(x[0] < xt) + int(x[1] < xt) + int(ri['xd'] < xt)
            tot += 1; agree_rule += int(rule == ri['n'])
            for e in EPS_B:
                agree[e] += int(byx[x].get(e) == ri['n'])
            pr = product_count(T, *x) if T is not None else None
            agree_prod += int(pr == ri['n'])
            L_.append(f"| {x[0]:g} | {x[1]:g} | {ri['xd']:.2f} | {', '.join(str(byx[x].get(e, '—')) for e in EPS_B)} | {ri['n']} | {rule} | "
                      f"{pr if pr is not None else '—'} |")
        L_ += ["", f"The loop's count equals the inner problem's at {', '.join(f'{agree[e]} of {tot} (ε = {e:g})' for e in EPS_B)}; the "
               f"independent rule at {agree_rule} of {tot}" + (f", the product rule at {agree_prod} of {tot}" if T is not None else '') + ".", ""]
        if T is not None:
            # the areas in the gaps' quarter-plane (x₁₀, x₀₁ ≥ 0) by the product rule on a fine grid; per vertex the offsets from it are
            # (2σ²ε^{1/2}/H)·(x₁₀, x₀₁), and a cell has four vertices: P = 4(2σ²/H²)²·area·ε; the faces' layers: 4·(2σ²/H²)·x_θ·ε^{1/2}
            g = np.arange(0.0, 14.0, 0.02) + 0.01
            P1, P2 = np.meshgrid(np.interp(g, T['x'], T['p'], right=0.0), np.interp(g, T['x'], T['p'], right=0.0), indexing='ij')
            third = np.minimum(P1 * (1 - P2), (1 - P1) * P2)
            a3, a4 = float((third > TH_REL).sum() * 0.02 ** 2), float((P1 * P2 > TH_REL).sum() * 0.02 ** 2)
            xt_T = float(np.interp(-TH_REL, -T['p'], T['x']))
            k = 2 * SIG ** 2 / H ** 2
            L_ += [f"The product rule's regions in the gaps' quarter-plane at a vertex: N ≥ 3 of area {a3:.1f}, N = 4 of {a4:.1f} (x_θ = "
                   f"{xt_T:.3f} from the table; a square x_θ² = {xt_T ** 2:.1f} and a triangle x_θ²/2 = {xt_T ** 2 / 2:.1f} at the crudest). For a "
                   f"feature uniform in a Voronoi cell the offsets from a vertex are (2σ²ε^{{1/2}}/H)·(x₁₀, x₀₁) at leading order, so P(N ≥ 2) = "
                   f"4·(2σ²/H²)x_θε^{{1/2}} = {4 * k * xt_T:.1f}·ε^{{1/2}} (the layers along the four faces, less O(ε) where they overlap), "
                   f"P(N ≥ 3) = 4(2σ²/H²)²·{a3:.1f}·ε = {4 * k * k * a3:.0f}·ε and P(N = 4) = {4 * k * k * a4:.0f}·ε.", ""]
    if C is not None:
        L_ += ["## C. Several features in the plane", "",
               f"Random three-feature targets (seed {SEED_C}; centres in [−0.7, 0.7]² at least 3σ apart, amplitudes ±U(0.4, 1)) on the "
               f"lattice over [−1, 1]² ({int(round(2 / H)) + 1}² kernels; the samples' grid {int(round(2 / DS)) + 1}²); the sliding path against "
               "the loop to 1.3 times the path's last entry:", "",
               "| target | the path (entries) | ε | formed = the path's | the never-entering kernels' largest \\|w\\|·ε^{−3/2} | their largest "
               "displacement/ε | the block: max \\|w − w_path\\|, max \\|c − c_path\\| | rms residual |", "|---|---|---|---|---|---|---|---|"]
        for q in sorted(set(r['q'] for r in C)):
            rs = sorted([r for r in C if r['q'] == q], key=lambda r: -r['eps'])
            for r in rs:
                L_.append(f"| {q} | {r['path']} ({', '.join(f'{t:.2f}' for t, j, k in r['events'])}) | {r['eps']:g} | "
                          f"{'yes' if r['formed'] == r['path'] else 'no: ' + str(r['formed'])} | {r['w_never'] / r['eps'] ** 1.5:.4g} | "
                          f"{r['disp_never'] / r['eps']:.4g} | {r['dw']:.2e}, {r['dc']:.2e} | {r['rms']:.2e} |")
            es = [r['eps'] for r in rs]
            L_[-1] += ''
            L_.append(f"| {q} | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | {slope(es[-3:], [r['w_never'] for r in rs][-3:]):.3f} | "
                      f"{slope(es[-3:], [r['disp_never'] for r in rs][-3:]):.3f} | {slope(es[-3:], [r['dw'] for r in rs][-3:]):.3f}, "
                      f"{slope(es[-3:], [r['dc'] for r in rs][-3:]):.3f} | |")
        L_.append("")
    if D is not None and T is not None:
        xt_T = float(np.interp(-TH_REL, -T['p'], T['x']))
        k = 2 * SIG ** 2 / H ** 2
        L_ += ["## D. The probabilities over uniformly placed features", "",
               f"A = {A_B:g} at {N_D} uniform positions in the Voronoi cell of the kernel at (0, 0) (seed {SEED_D}); the loop's count N against "
               "the product rule at each position's gaps to its two nearest faces, and the fractions against the asymptotics (w = "
               f"2σ²x_θε^{{1/2}}/H, P(N ≥ 2) = 4w/H − 4w²/H² from the faces' layers; P(N ≥ 3) from B's areas):", "",
               "| ε | N = 1, 2, 3, 4 (loop) | the product rule agrees | P(N ≥ 2): loop (± s.e.) | asymptotic | P(N ≥ 3): loop | asymptotic |",
               "|---|---|---|---|---|---|---|"]
        g = np.arange(0.0, 14.0, 0.02) + 0.01
        P1, P2 = np.meshgrid(np.interp(g, T['x'], T['p'], right=0.0), np.interp(g, T['x'], T['p'], right=0.0), indexing='ij')
        a3 = float((np.minimum(P1 * (1 - P2), (1 - P1) * P2) > TH_REL).sum() * 0.02 ** 2)
        for e in EPS_D:
            rs = [r for r in D if r['eps'] == e]
            n = np.array([r['n'] for r in rs])
            pr = np.array([product_count(T, *r['gaps']) for r in rs])
            w = k * xt_T * np.sqrt(e)
            p2, p3 = float(np.mean(n >= 2)), float(np.mean(n >= 3))
            L_.append(f"| {e:g} | {', '.join(str(int((n == j).sum())) for j in (1, 2, 3, 4))} | {int((n == pr).sum())} of {len(rs)} | "
                      f"{p2:.4f} ± {np.sqrt(p2 * (1 - p2) / len(rs)):.4f} | {4 * w - 4 * w * w:.4f} | {p3:.4f} | {4 * k * k * a3 * e:.4f} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_transit_nd.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'transit_nd_raw.npz')
    md = os.path.join(OUT, 'transit_nd.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in 'ABCTD']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCTD')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('T', part_T), ('D', part_D)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in 'ABCTD']))


if __name__ == '__main__':
    main()
