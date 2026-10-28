"""
Theorem D‴.55: the transit's certificate — exact recovery of separated features by their starting kernels of largest correlation
(the nearest in the continuum), the path's skeleton and the fitted skeleton as predictors of the outcome at a finite separation
====================================================================================================================================

The sliding inverse scale space (`rich_transit_nd.sliding_iss_nd`, in one dimension on a 1-D `BumpsND` — the same equations as
`rich_transit.sliding_iss` with the exact Jacobian, which keeps a block of a large kernel and a seed from stalling the integrator
(referee round 46): the limit of the free law with moving centres from small masses) on random targets of K separated Gaussian
features, against the certificate on the starting set (the samples' inner product; e_k the starting point of largest correlation with
the k-th feature, the nearest in the continuum) and the path's skeletons; the unpolished flow (polish off) as the audit of the
least-squares polish on A and D(b):

A. One dimension (D‴.48's dictionary: 21 kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]): three features at
   least DMIN_A apart, amplitudes ±U(0.4, 1): the certificate's margin m = 1 − max_{j∉E}|η(c_j)| against the path's outcome — exact with
   the nearest kernels (the active set E, the residual zero, no exits), exact with others (three kernels, the residual zero), or not.
B. The amplitudes' range: the same positions and signs with |A_k| = 10^{−r·u_k} (u_k uniform in [0, 1]) for r = R_B — the certificate
   does not change; does the outcome?
C. The plane (D‴.54's lattice: 17 × 17 kernels of width 0.15, spacing 0.125, 41² samples of [−1, 1]²): three features at least
   DMIN_C apart, amplitudes ±U(0.4, 1).
D. Designed families (referee round 43): pairs on starting points with amplitude ratios, pairs at s·σ over a grid cell, every order of
   the amplitudes.
E. The two-feature map: the first feature on a starting point, the second at the distance D with the amplitude ±10^ℓ, over a grid of
   (D, ℓ): the fitted skeleton's outcome everywhere (no flow), the threshold separation D*(ℓ), the other margins against it, the linear
   skeleton's error against the fitted one's in the ratio (its order), and the path at sampled points against it.
F. A wide range: three features at 5σ and 3σ, the amplitudes 1, 10⁻³, 10⁻⁶ in all six orders, same and alternating signs (exits only
   below 10⁻¹² of the largest coefficient).

Usage: python experiments/exp_transit_certificate.py [--only A,B,C,D,E,F] [--from-raw] → results/chaos/transit_certificate.md (+ _raw.npz).
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
import exp_transit_mass as TM
import exp_transit_nd as TN

OUT = os.environ.get('TC_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TC_WORKERS', '2'))
K_F = 3
SEED_A, N_A, DMIN_A = 551, 80, (0.25, 0.3, 0.35, 0.45, 0.6)
SEED_B, N_B, DMIN_B, R_B = 552, 40, 0.35, (0, 1, 2, 3)
SEED_C, N_C, DMIN_C = 553, 40, (0.3, 0.45, 0.6)
RMS_TOL = 1e-8
BANDS = (0.5, 0.8, 1.0)                                                # the band |m| ≤ Cq around the certificate's threshold
T_MAX = 600                                                            # D: a path's wall-clock limit (s)
LAT_OFF, LAT_D, RATIOS_D = 1e-6, (0.3, 0.4, 0.5), (1e-3, 0.02, 0.1, 1 / 3, 1.0, 3.0, 10.0, 50.0, 1e3)
SEED_D, SEPS_D, N_OFF_D = 554, (1.5, 2.0, 2.5, 3.0, 3.5, 4.0), 20       # D(b): pairs at s·σ, the midpoint uniform over a grid cell
N_ORD_D, RANGE_D = 30, 10.0                                            # D(c): A's targets outside 0.8q, the amplitudes 1, 10^{−1/2}, 10^{−1} in every order
D_E, L_E = np.round(np.arange(0.12, 0.6001, 0.01), 10), np.round(np.arange(-3.0, 3.0001, 0.125), 10)  # E: the map (D, log₁₀ ratio)
SEED_E, N_EV, LMAX_EV, T_EV = 555, 120, 2.5, 300                      # E: paths per sign checking the map (|log₁₀ ratio| ≤ LMAX_EV)
M0_E = 0.05                                                            # E: the second feature's own near ties (m₀ ≤ M0_E) set apart
F_SEPS, F_MAGS, F_OFF = (0.6, 0.36), (1.0, 1e-3, 1e-6), 0.013          # F: 5σ and 3σ, the range 10⁶, the features off the grid


def entrants(K, c0, cs):
    """E: for each feature the starting point of largest correlation with it on the samples (the nearest in the continuum; on samples
    they differ near the window's edge — referee round 46), and the correlations ⟨φ(c_j)φ(c*_k)⟩."""
    c0 = np.asarray(c0, float).reshape(len(c0), -1)
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    P0 = K.phi(c0 if c0.shape[1] > 1 else c0[:, 0])
    Pf = K.phi(cs if cs.shape[1] > 1 else cs[:, 0])
    gram = P0.T @ Pf / P0.shape[0]
    return [int(np.argmax(gram[:, k])) for k in range(len(cs))], gram


def certificate(K, c0, cs, As):
    """η at every starting point (the samples' inner product), E (`entrants`), the margin m off E, and the cross-talk q (continuum,
    q = max_k Σ_{l≠k} G((Δ_kl − 2R)₊)/G(R), G(d) = e^{−d²/(4σ²)})."""
    c0 = np.asarray(c0, float).reshape(len(c0), -1)
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    d0 = np.linalg.norm(c0[:, None, :] - cs[None, :, :], axis=2)
    e, gram = entrants(K, c0, cs)
    own = np.array([gram[e[k], k] for k in range(len(cs))])
    eta = (np.sign(As)[None, :] * gram / own[None, :]).sum(1)
    off = np.array([j for j in range(len(c0)) if j not in e])
    sig = K.sig
    G = lambda d: np.exp(-np.asarray(d, float) ** 2 / (4 * sig ** 2))
    R = float(max(d0[e[k], k] for k in range(len(cs))))
    D = np.linalg.norm(cs[:, None, :] - cs[None, :, :], axis=2)
    q = float(max(sum(G(max(D[k, l] - 2 * R, 0.0)) for l in range(len(cs)) if l != k) for k in range(len(cs)))) / float(G(R))
    # the relative margin (referee round 43): the entrants' own values s_kη(c_{e_k}) — 1 up to the cross-talk — against the largest off E;
    # the single-feature (Voronoi) margin: each feature alone, no signs or sums
    m_rel = float(min(np.sign(As[k]) * eta[e[k]] for k in range(len(cs))) - np.abs(eta[off]).max())
    m_vor = float(min(1 - max(gram[j, k] for j in range(len(c0)) if j != e[k]) / gram[e[k], k] for k in range(len(cs))))
    return dict(eta=eta, e=e, m=float(1 - np.abs(eta[off]).max()), argm=int(off[np.argmax(np.abs(eta[off]))]), q=q, R=R,
                unique=len(set(e)) == len(e), m_rel=m_rel, m_vor=m_vor)


def skeleton(K, c0, cs, As, fit=False, atoms=False):
    """The sliding path's skeleton (referee round 43's screening): each block replaced by its exact fit — the fitted features' kernels at
    their centres, the residual between entries the target's projection off their tangent space span{φ(c*_m), ∇_cφ(c*_m)} — the duals
    piecewise linear in τ at ⟨rφ(c_j)⟩ from the starting points. It runs the entries in order (each the first dual to reach ±1) until all K
    features are fitted or a kernel off E enters. Returns the entrants, whether they are E (ok), and the margin: 1 − max |u_j| over the
    entries and the kernels off E, or, if a kernel off E enters, minus the rightful entrants' least distance from their walls then.
    fit: the fitted skeleton — each block at its own sign-constrained least-squares fit started at its features (`polish_block`: the
    block's stationary point near them, the projection's nonlinear completion) — the skeleton without its linearisation. atoms: the
    baseline without the tangent directions — the fitted features' atoms A_mφ(c*_m) subtracted (referee round 46). 'phase': the number
    of features fitted when the outcome is decided (the intruder's entry, or the phase of the least margin)."""
    c0 = np.asarray(c0, float).reshape(len(c0), -1)
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    sq = lambda c: c if c.shape[1] > 1 else c[:, 0]
    nx = K.xs.size if hasattr(K, 'xs') and np.ndim(K.xs) == 1 else K.N
    Pf = K.phi(sq(cs))
    f = Pf @ np.asarray(As, float)
    P0 = K.phi(sq(c0))
    e = entrants(K, c0, cs)[0]
    off = [j for j in range(len(c0)) if j not in e]
    u, act, fitted, ent, worst, wphase = np.zeros(len(c0)), np.zeros(len(c0), bool), [], [], 0.0, 0
    while len(fitted) < len(cs):
        if fitted and atoms:
            r = f - Pf[:, fitted] @ np.asarray(As, float)[fitted]
        elif fitted and fit:
            a_, s_ = np.asarray(As, float)[fitted], np.sign(np.asarray(As, float)[fitted])
            if cs.shape[1] > 1:
                w_, c_ = RN.polish_block_nd(K, f, a_, cs[fitted], s_)
                r = f - K.synth(c_, w_)
            else:
                w_, c_ = RT.polish_block(K, f, a_, cs[fitted, 0], s_)
                r = f - K.phi(c_) @ w_
        elif fitted:
            Dp = K.dphi(sq(cs[fitted])).reshape(nx, -1)
            B = np.c_[Pf[:, fitted], Dp]
            r = f - B @ np.linalg.lstsq(B, f, rcond=None)[0]
        else:
            r = f
        v = P0.T @ r / nx
        ina = ~act & (v != 0)
        if not ina.any():
            break
        dt = np.full(len(c0), np.inf)
        dt[ina] = (1 - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
        j = int(np.argmin(dt))
        u[ina] += dt[j] * v[ina]
        u[j] = np.sign(u[j])
        act[j] = True
        ent.append(j)
        if j not in e or e.index(j) in fitted:
            # an intruder: the margin is minus the rightful entrants' distance from their walls at its entry (a near tie when small)
            gap = min(1 - abs(u[e[m]]) for m in range(len(cs)) if m not in fitted)
            return dict(entrants=ent, ok=False, margin=float(-gap), phase=len(fitted))
        cur = float(np.abs(u[[i for i in off if not act[i]]]).max()) if any(not act[i] for i in off) else 0.0
        if cur > worst:
            worst, wphase = cur, len(fitted)
        fitted.append(e.index(j))
    return dict(entrants=ent, ok=True, margin=float(1 - worst), phase=wphase)


def outcome(R, e, K):
    act = sorted(int(j) for j in np.flatnonzero(R['act']))
    clean = R['hist'][-1][3] < RMS_TOL and all(ev[2] == 'enter' for ev in R['events'])
    if clean and act == sorted(e) and len(R['events']) == K:
        return 'nearest'
    if R['hist'][-1][3] < RMS_TOL and len(act) == K:
        return 'other'
    return 'not'


def targets(rng, n, dmin, lo=-0.7, hi=0.7, dim=1):
    out = []
    while len(out) < n:
        cs = rng.uniform(lo, hi, (K_F, dim))
        if min(np.linalg.norm(cs[i] - cs[j]) for i in range(K_F) for j in range(i)) < dmin:
            continue
        out.append((cs[:, 0] if dim == 1 else cs, rng.choice([-1.0, 1.0], K_F), rng.uniform(0.4, 1.0, K_F), rng.uniform(0, 1, K_F)))
    return out


def path_1d(L, K, f, c0, polish=True, exit_tol=1e-6):
    """The sliding path in one dimension by `sliding_iss_nd` on the 1-D samples (the exact Jacobian)."""
    return RN.sliding_iss_nd(L, RN.BumpsND(K.xs[:, None], K.sig), f, np.asarray(c0, float)[:, None], np.ones(len(c0)), polish=polish,
                             exit_tol=exit_tol)


def job_1d(args):
    """One target: the certificate and the path (polished); with flow=True also the unpolished flow's outcome (the polish's audit)."""
    cs, As = args[:2]
    flow = len(args) > 2 and args[2]
    exit_tol = args[3] if len(args) > 3 else 1e-6
    warnings.filterwarnings('ignore')
    L, K, c0 = TM.setup()
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(cs=cs, As=As, skip=True)
    f = X.bumps_target(K, cs, As)
    R = path_1d(L, K, f, c0, exit_tol=exit_tol)
    out = dict(cs=cs, As=As, skip=False, m=cert['m'], q=cert['q'], argm=cert['argm'], e=cert['e'], out=outcome(R, cert['e'], len(cs)),
               order=[int(j) for t, j, k, c in R['events']], rms=R['hist'][-1][3], m_rel=cert['m_rel'], m_vor=cert['m_vor'],
               rej=int(R['polish_rejected']))
    if flow:
        t0 = time.time()
        Rf = path_1d(L, K, f, c0, polish=False, exit_tol=exit_tol)
        out.update(out_flow=outcome(Rf, cert['e'], len(cs)), order_flow=[int(j) for t, j, k, c in Rf['events']],
                   flow_relaxed=bool(all(Rf['relaxed'])), flow_secs=time.time() - t0)
    return out


def job_d(args):
    """D: `job_1d` under a wall-clock limit T_MAX (out = 'timeout', the certificate still recorded); the family's tag kept."""
    import signal
    tag, cs, As = args[:3]
    tmax = int(args[3]) if len(args) > 3 else T_MAX
    flow = len(args) > 4 and args[4]
    exit_tol = args[5] if len(args) > 5 else 1e-6

    def alarm(signum, frame):
        raise TimeoutError
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(tmax)
    t0 = time.time()
    try:
        r = job_1d((cs, As, flow, exit_tol))
    except TimeoutError:
        L, K, c0 = TM.setup()
        cert = certificate(K, c0, cs, As)
        r = dict(cs=cs, As=As, skip=not cert['unique'], m=cert['m'], q=cert['q'], argm=cert['argm'], e=cert['e'], out='timeout', order=[],
                 rms=np.nan, m_rel=cert['m_rel'], m_vor=cert['m_vor'], rej=-1)
    finally:
        signal.alarm(0)
    r.update(tag=tag, secs=time.time() - t0)
    return r


def part_D(A):
    """Designed families (referee round 43): (a) pairs on starting points (LAT_OFF off) at the distances LAT_D, same and opposite signs,
    the second amplitude the first's times RATIOS_D; (b) pairs at s·σ (SEPS_D), equal amplitudes, same and opposite signs, the midpoint
    uniform over a grid cell; (c) N_ORD_D of A's targets with |m| > 0.8q, the amplitudes 1, RANGE_D^{−1/2}, RANGE_D^{−1} in every order."""
    from itertools import permutations
    from multiprocessing import Pool
    jobs = []
    for D in LAT_D:
        for sg in (1.0, -1.0):
            for rt in RATIOS_D:
                jobs.append((('a', D, sg, rt), np.array([-0.1, -0.1 + D]) + LAT_OFF, np.array([1.0, sg * rt])))
    rng = np.random.default_rng(SEED_D)
    for s_ in SEPS_D:
        for sg in (1.0, -1.0):
            for i in range(N_OFF_D):
                mid = rng.uniform(-0.05, 0.05)
                jobs.append((('b', s_, sg, i), np.array([mid - s_ * 0.06, mid + s_ * 0.06]), np.array([1.0, sg]), T_MAX, True))
    out_ = [r for r in A if not r['skip'] and abs(r['m']) > 0.8 * r['q']]
    pick = sorted(rng.choice(len(out_), size=min(N_ORD_D, len(out_)), replace=False))
    mags = RANGE_D ** -np.array([0.0, 0.5, 1.0])
    for t in pick:
        r = out_[t]
        for perm in permutations(range(3)):
            jobs.append((('c', int(t), perm), np.array(r['cs']), np.sign(r['As']) * mags[list(perm)]))
    with Pool(N_WORKERS) as pool:
        return pool.map(job_d, jobs, chunksize=1)


def e_point(args):
    """E: one point of the two-feature map — the margins (the certificate's, the relative, the single-feature) and both skeletons'."""
    D, l, sg = args
    warnings.filterwarnings('ignore')
    L, K, c0 = TM.setup()
    cs = np.array([-0.1, -0.1 + D]) + LAT_OFF
    As = np.array([1.0, sg * 10.0 ** l])
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(D=D, l=l, sg=sg, skip=True)
    sk, fs = skeleton(K, c0, cs, As), skeleton(K, c0, cs, As, fit=True)
    return dict(D=D, l=l, sg=sg, skip=False, m=cert['m'], q=cert['q'], m_rel=cert['m_rel'], m_vor=cert['m_vor'], m_skel=sk['margin'],
                m_fit=fs['margin'], fit_ok=bool(fs['ok']), m_atoms=skeleton(K, c0, cs, As, atoms=True)['margin'], phase=fs['phase'])


def part_E():
    """The two-feature map on the grid D_E × L_E for both signs; the path (under T_EV s) at N_EV sampled points per sign, |ℓ| ≤ LMAX_EV."""
    from multiprocessing import Pool
    pts = [(float(D), float(l), sg) for sg in (1.0, -1.0) for D in D_E for l in L_E]
    rng = np.random.default_rng(SEED_E)
    pick = []
    for sg in (1.0, -1.0):
        ch = [p for p in pts if p[2] == sg and abs(p[1]) <= LMAX_EV]
        pick += [ch[i] for i in sorted(rng.choice(len(ch), size=N_EV, replace=False))]
    with Pool(N_WORKERS) as pool:
        grid = pool.map(e_point, pts, chunksize=40)
        paths = pool.map(job_d, [(('e', D, l, sg), np.array([-0.1, -0.1 + D]) + LAT_OFF, np.array([1.0, sg * 10.0 ** l]), T_EV)
                                 for D, l, sg in pick], chunksize=1)
    return dict(grid=grid, paths=paths)


def job_2d(args):
    cs, As = args
    warnings.filterwarnings('ignore')
    L = X.law()
    K, c0 = TN.full()
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(cs=cs, As=As, skip=True)
    R = RN.sliding_iss_nd(L, K, K.synth(cs, As), c0, np.ones(len(c0)))
    return dict(cs=cs, As=As, skip=False, m=cert['m'], q=cert['q'], argm=cert['argm'], e=cert['e'], out=outcome(R, cert['e'], len(cs)),
                order=[int(j) for t, j, k, c in R['events']], rms=R['hist'][-1][3], rej=int(R['polish_rejected']))


def part_A():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_A)
    jobs, tags = [], []
    for dmin in DMIN_A:
        for cs, s, a, u in targets(rng, N_A, dmin):
            jobs.append((cs, s * a, True)); tags.append(dmin)
    with Pool(N_WORKERS) as pool:
        rs = pool.map(job_1d, jobs, chunksize=4)
    return [dict(dmin=t, **r) for t, r in zip(tags, rs)]


def part_F():
    """Three features at 5σ and 3σ (off the grid by F_OFF), the amplitudes F_MAGS in all six orders, signs + + + and + − +; exits only
    below 10⁻¹² of the largest coefficient."""
    from itertools import permutations
    from multiprocessing import Pool
    jobs = []
    for s_ in F_SEPS:
        for sg in ((1.0, 1.0, 1.0), (1.0, -1.0, 1.0)):
            for perm in permutations(range(3)):
                jobs.append((('f', s_, sg, perm), np.array([-s_, 0.0, s_]) + F_OFF, np.array(sg) * np.array(F_MAGS)[list(perm)], T_MAX,
                             False, 1e-12))
    with Pool(N_WORKERS) as pool:
        return pool.map(job_d, jobs, chunksize=1)


def part_B():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_B)
    jobs, tags = [], []
    for i, (cs, s, a, u) in enumerate(targets(rng, N_B, DMIN_B)):
        for r in R_B:
            jobs.append((cs, s * 10.0 ** (-r * u))); tags.append((i, r))
    with Pool(N_WORKERS) as pool:
        rs = pool.map(job_1d, jobs, chunksize=4)
    return [dict(cfg=t[0], r=t[1], **x) for t, x in zip(tags, rs)]


def part_C():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    jobs, tags = [], []
    for dmin in DMIN_C:
        for cs, s, a, u in targets(rng, N_C, dmin, dim=2):
            jobs.append((cs, s * a)); tags.append(dmin)
    with Pool(N_WORKERS) as pool:
        rs = pool.map(job_2d, jobs, chunksize=2)
    return [dict(dmin=t, **r) for t, r in zip(tags, rs)]


def confusion(rows):
    rows = [r for r in rows if not r['skip']]
    pos = [r for r in rows if r['m'] > 0]
    neg = [r for r in rows if r['m'] <= 0]
    cnt = lambda rs, o: sum(r['out'] == o for r in rs)
    return rows, pos, neg, cnt


def fill_margins(rows, dim=1):
    """The skeletons' margins recomputed from the targets at report time (a stored margin would hide a change of `skeleton`: referee
    round 46) — the linear skeleton, the fitted one (and the phase at which it decides) and the atoms-only baseline —, and the relative and
    single-feature margins where a row lacks them."""
    if dim == 1:
        L, K, c0 = TM.setup()
    else:
        K, c0 = TN.full()
    for r in rows:
        if r['skip']:
            continue
        cs, As = np.asarray(r['cs']), np.asarray(r['As'])
        if 'm_rel' not in r:
            c = certificate(K, c0, cs, As)
            r['m_rel'], r['m_vor'] = c['m_rel'], c['m_vor']
        r['m_skel'] = skeleton(K, c0, cs, As)['margin']
        fs = skeleton(K, c0, cs, As, fit=True)
        r['m_fit'], r['phase'] = fs['margin'], fs['phase']
        r['m_atoms'] = skeleton(K, c0, cs, As, atoms=True)['margin']
    return rows


PREDICTORS = (('m', 'the certificate, 1 − max off E'), ('m_rel', 'the relative margin'), ('m_vor', 'the single feature (no signs or sums)'),
              ('m_atoms', 'the atoms subtracted (no tangent directions)'), ('m_skel', "the path's skeleton"), ('m_fit', 'the fitted skeleton'))
SHORT = dict(m='the certificate', m_rel='relative', m_vor='single-feature', m_atoms='atoms subtracted', m_skel='the skeleton',
             m_fit='the fitted skeleton')                              # D(b)'s column labels


def predictors(rows):
    """How often each predictor's sign (and the baselines 'always E' and 'E iff the path's first entrant is in E') matches exactness with
    E; where computed, how often the unpolished flow's outcome is the polished path's."""
    rs = [r for r in rows if not r['skip'] and r['out'] != 'timeout']
    ok = [r['out'] == 'nearest' for r in rs]
    parts = [f"{lab} {sum((r[k] > 0) == o for r, o in zip(rs, ok))}" for k, lab in PREDICTORS if all(k in r for r in rs)]
    first = [bool(r['order']) and r['order'][0] in r['e'] for r in rs]
    txt = (f"of {len(rs)}: " + ', '.join(parts) + f", always E {sum(ok)}, E iff the first entrant is in E {sum(f == o for f, o in zip(first, ok))} "
           f"(the first entrant off E — the outcome decided at the first entry — on {sum(not f for f in first)}")
    if all('phase' in r for r in rs):
        txt += f"; the fitted skeleton decided at its first phase on {sum(r['phase'] == 0 for r in rs)}"
    txt += ")"
    fl = [r for r in rs if 'out_flow' in r]
    if fl:
        txt += (f"; the unpolished flow's outcome the polished path's on {sum(r['out_flow'] == r['out'] for r in fl)} of {len(fl)} (its "
                f"blocks relaxed on {sum(r['flow_relaxed'] for r in fl)})")
    return txt


def band(rs):
    """Outside the band |m| ≤ Cq: how many targets, and on how many the certificate's sign decides; the paths' ends."""
    parts = []
    for C in BANDS:
        ob = [r for r in rs if abs(r['m']) > C * r['q']]
        parts.append(f"C = {C:g}: {sum((r['m'] > 0) == (r['out'] == 'nearest') for r in ob)} of {len(ob)}")
    return ("Outside the band |m| ≤ Cq the sign decides on " + '; '.join(parts) + f". Every path ends at an exact fit with "
            f"{K_F} kernels on {sum(r['out'] != 'not' for r in rs)} of {len(rs)} targets.")


def report(A=None, B=None, C=None, D=None, E_=None, F=None):
    L_ = ["# Theorem D‴.55: the transit's certificate (`exp_transit_certificate.py`)", "",
          "The sliding inverse scale space (the free law's limit from small masses, moving centres) on random targets of three Gaussian "
          "features of the kernels' width; the certificate η(c_j) = Σ_k sign(A_k)⟨φ(c*_k)φ(c_j)⟩/⟨φ(c*_k)φ(c_{e_k})⟩ on the starting points "
          "(the samples' inner product, e_k the starting point of largest correlation with feature k — the nearest in the continuum —, E "
          "their set), its margin m = 1 − max_{j∉E}|η(c_j)| and "
          "the cross-talk q = max_k Σ_{l≠k}G((Δ_kl − 2R)₊)/G(R). The path's outcome: exact with E (its active set E, the residual below "
          f"{RMS_TOL:g}, no exits), exact with other kernels (three kernels, the residual below {RMS_TOL:g}), or not exact. The paths by "
          "`sliding_iss_nd` (in one dimension on the 1-D samples, with the exact Jacobian), the block polished by its least-squares fit "
          "unless a kernel below 10⁻² of the block's largest coefficient would land beyond one kernel width of both its start and its "
          "correlation's ascent's end (a rejected polish; the bound does not check single-kernel blocks or the block's largest kernel); the "
          "skeletons recomputed from the targets.", ""]

    def table(rows, key, label):
        out = [f"| {label} | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | "
               "the disagreements' m (q) |", "|---|---|---|---|---|---|"]
        for v in sorted(set(r[key] for r in rows)):
            rs, pos, neg, cnt = confusion([r for r in rows if r[key] == v])
            bad = [r for r in rs if (r['m'] > 0) != (r['out'] == 'nearest')]
            out.append(f"| {v:g} | {len(rs)} | {cnt(pos, 'nearest')}, {cnt(pos, 'other')}, {cnt(pos, 'not')} | {cnt(neg, 'nearest')}, "
                       f"{cnt(neg, 'other')}, {cnt(neg, 'not')} | {np.median([r['q'] for r in rs]):.3f} | "
                       + (', '.join(f"{r['m']:+.4f} ({r['q']:.2f})" for r in bad) or '—') + " |")
        return out
    if A is not None:
        L_ += ["## A. One dimension", "", f"{N_A} targets per minimal separation (seed {SEED_A}; σ = 0.12, the grid's spacing 0.1); "
               "amplitudes ±U(0.4, 1):", ""] + table(A, 'dmin', 'Δ ≥')
        rs, pos, neg, cnt = confusion(A)
        bad = [r for r in rs if (r['m'] > 0) != (r['out'] == 'nearest')]
        L_ += ["", f"Over all {len(rs)} targets the certificate's sign predicts exactness with E on {len(rs) - len(bad)}; the "
               f"disagreements' margins |m| ≤ {max([abs(r['m']) for r in bad], default=0):.4f}, at most "
               f"{max([abs(r['m']) / r['q'] for r in bad], default=0):.4f} times their cross-talk q. " + band(rs),
               "", "Exactness with E predicted by the sign of each margin (the relative margin min_k s_kη(c_{e_k}) − max_{j∉E}|η(c_j)|; the "
               "single-feature margin min_k(1 − max_{j≠e_k}⟨φ_kφ_j⟩/⟨φ_kφ_{e_k}⟩)), " + predictors(fill_margins(A)) + ".", ""]
    if B is not None:
        cfgs = sorted(set(r['cfg'] for r in B))
        same, n_ok, flip = 0, 0, []
        for c in cfgs:
            rs = [r for r in B if r['cfg'] == c and not r['skip']]
            outs = [r['out'] for r in sorted(rs, key=lambda r: r['r'])]
            if len(set(outs)) == 1:
                same += 1
            else:
                flip.append((c, outs, rs[0]['m']))
            n_ok += int(rs[0]['m'] > 0)
        orders = sum(len(set(tuple(r['order']) for r in B if r['cfg'] == c and not r['skip'])) > 1 for c in cfgs)
        ratio = {c: min(abs(r['m']) / r['q'] for r in B if r['cfg'] == c and not r['skip']) for c in cfgs}
        rng_r = {rr: [float(np.abs(r['As']).max() / np.abs(r['As']).min()) for r in B if r['r'] == rr and not r['skip']] for rr in R_B}
        keep = [c for c in cfgs if ratio[c] > 0.8 and c not in [f[0] for f in flip]]
        L_ += ["## B. The amplitudes' range", "",
               f"{N_B} configurations (seed {SEED_B}; positions at least {DMIN_B:g} apart, signs) with |A_k| = 10^{{−r·u_k}}, r = "
               f"{', '.join(str(r) for r in R_B)} (the amplitudes' range up to 10^{max(R_B)}); the certificate depends on the positions and "
               f"signs only ({n_ok} of {len(cfgs)} with m > 0). The outcome is the same at every range on {same} of {len(cfgs)} "
               f"configurations, while the order of the entries changes with the range on {orders}"
               + (f"; it changes on {len(flip)}: " + '; '.join(f"configuration {c} (m = {m:+.4f}): {', '.join(o)}" for c, o, m in flip)
                  if flip else '') + (f". The configurations whose outcome changes have |m|/q ≤ {max(ratio[c] for c, o, m in flip):.2f}; of "
                                      f"those with |m| > 0.8q, {len(keep)} of {sum(ratio[c] > 0.8 for c in cfgs)} keep their outcome at every "
                                      "range" if flip else '') + ". The ranges realised, max|A|/min|A|, median [max]: " +
               ', '.join(f"r = {rr}: {np.median(v):.3g} [{max(v):.3g}]" for rr, v in rng_r.items()) + "; the amplitudes' order is the same "
               "for every r > 0 (u_k fixed), so this Part varies the range, not the order (D(c) varies the order).",
               "", "Exactness with E predicted by the margins' signs, " + predictors(fill_margins(B)) + ".", ""]
    if C is not None:
        L_ += ["## C. The plane", "", f"{N_C} targets per minimal separation (seed {SEED_C}; σ = 0.15, the lattice's spacing 0.125; "
               "centres in [−0.7, 0.7]²), amplitudes ±U(0.4, 1):", ""] + table(C, 'dmin', 'Δ ≥')
        rs, pos, neg, cnt = confusion(C)
        bad = [r for r in rs if (r['m'] > 0) != (r['out'] == 'nearest')]
        L_ += ["", f"Over all {len(rs)} targets the certificate's sign predicts exactness with E on {len(rs) - len(bad)}. " + band(rs),
               "", "Exactness with E predicted by the margins' signs, " + predictors(fill_margins(C, dim=2)) + ".", ""]
    if D is not None:
        code = {'nearest': 'E', 'other': 'O', 'not': '×', 'timeout': 't'}
        D = fill_margins(D)                                            # the skeleton's margin (and older rows' relative margins)
        da = [r for r in D if r['tag'][0] == 'a']
        db = [r for r in D if r['tag'][0] == 'b']
        dc = [r for r in D if r['tag'][0] == 'c']
        L_ += ["## D. Designed families (referee round 43)", "",
               f"(a) Two features on starting points ({LAT_OFF:g} off) at −0.1 and −0.1 + D, the first's amplitude 1, the second's ±ρ_A: the "
               "margins (amplitude-free) and the outcome at each ρ_A (E exact with E, O exact with others, × not exact, t over "
               f"{T_MAX} s):", "",
               "| D | signs | m | m/q | relative | single-feature | " + ' | '.join(f"{x:.3g}" for x in RATIOS_D) + " |",
               '|---|---|---|---|---|---|' + '---|' * len(RATIOS_D)]
        for Dd in LAT_D:
            for sg in (1.0, -1.0):
                rs = sorted([r for r in da if r['tag'][1] == Dd and r['tag'][2] == sg], key=lambda r: r['tag'][3])
                r0 = rs[0]
                L_.append(f"| {Dd:g} | {'+ +' if sg > 0 else '+ −'} | {r0['m']:+.4f} | {r0['m'] / r0['q']:+.3f} | {r0['m_rel']:+.4f} | "
                          f"{r0['m_vor']:+.4f} | " + ' | '.join(code[r['out']] for r in rs) + " |")
        L_ += ["", f"(b) Two features at s·σ, equal amplitudes, the midpoint uniform over a grid cell ({N_OFF_D} per row): the outcomes and "
               "how often each margin's sign is right:", "",
               "| s | signs | E, O, ×, t | m > 0 | " + ' | '.join(f"{SHORT[k]} right" for k, _ in PREDICTORS) + " | always E right |",
               '|---' * (5 + len(PREDICTORS)) + '|']
        for s_ in SEPS_D:
            for sg in (1.0, -1.0):
                rs = [r for r in db if r['tag'][1] == s_ and r['tag'][2] == sg and not r['skip']]
                rr = [r for r in rs if r['out'] != 'timeout']
                ok = [r['out'] == 'nearest' for r in rr]
                L_.append(f"| {s_:g} | {'+ +' if sg > 0 else '+ −'} | " + ', '.join(str(sum(r['out'] == o for r in rs)) for o in code) +
                          f" | {sum(r['m'] > 0 for r in rs)} | " + ' | '.join(str(sum((r[k] > 0) == o for r, o in zip(rr, ok)))
                                                                         for k, _ in PREDICTORS) + f" | {sum(ok)} |")
        tg = sorted(set(r['tag'][1] for r in dc))
        vary = [t for t in tg if len(set(r['out'] for r in dc if r['tag'][1] == t)) > 1]
        L_ += ["", f"(c) {len(tg)} of A's targets with |m| > 0.8q, the amplitude magnitudes 1, {RANGE_D:g}^{{−1/2}}, {RANGE_D:g}^{{−1}} in all six "
               f"orders (the signs kept): the outcome changes with the order on {len(vary)} of them" +
               (': ' + '; '.join(f"m/q = {[r for r in dc if r['tag'][1] == t][0]['m'] / [r for r in dc if r['tag'][1] == t][0]['q']:+.2f}: "
                                 + ''.join(code[r['out']] for r in sorted([r for r in dc if r['tag'][1] == t], key=lambda r: r['tag'][2]))
                                 for t in vary) if vary else '') + ".", "",
               "Over D, exactness with E predicted by the margins' signs: (a) " + predictors(da) + "; (b) " + predictors(db) + "; (c) " +
               predictors(dc) + f". Paths over {T_MAX} s: {sum(r['out'] == 'timeout' for r in D)}; polishes rejected: "
               f"{sum(max(r.get('rej', 0), 0) for r in D)} on {sum(r.get('rej', 0) > 0 for r in D)} paths.", ""]
        bad_b = [r for r in db if not r['skip'] and r['out'] != 'nearest' and r['m_vor'] > M0_E]
        by_s = [(s_, sg, [r for r in bad_b if r['tag'][1] == s_ and r['tag'][2] == sg]) for s_ in SEPS_D for sg in (1.0, -1.0)]
        L_ += [f"D(b)'s failures that are not the second feature's own near ties (m₀ > {M0_E:g}), per row: " + ', '.join(
            f"{s_:g}σ {'+ +' if sg > 0 else '+ −'} {len(rs)}" + (" (m₀ = " + ', '.join(f"{r['m_vor']:.3f}" for r in rs) + ")"
                                                              if 0 < len(rs) <= 4 else '') for s_, sg, rs in by_s) + ".", ""]
        dd = [r for r in D if not r['skip'] and r['out'] != 'timeout']
        sent = []
        for C_ in BANDS:                                               # the band of A and C, tested on the designed families
            ob = [r for r in dd if abs(r['m']) > C_ * r['q']]
            wr = [r for r in ob if (r['m'] > 0) != (r['out'] == 'nearest')]
            sent.append(f"|m| > {C_:g}q: wrong on {len(wr)} of {len(ob)}" + (f" (|m|/q up to {max(abs(r['m']) / r['q'] for r in wr):.3f}, "
                                                                          f"{sum(r['m'] > 0 for r in wr)} of them with m > 0)" if wr else ""))
        L_ += ["Outside the band, the certificate's sign on the designed targets: " + '; '.join(sent) + ".", ""]
    if E_ is not None:
        g = [r for r in E_['grid'] if not r['skip']]
        sig = 0.12
        L_ += ["## E. The two-feature map", "",
               f"The first feature (amplitude 1) on the starting point −0.1 ({LAT_OFF:g} off), the second at the distance D = {D_E[0]:g} … "
               f"{D_E[-1]:g} (step {D_E[1] - D_E[0]:g}; σ = {sig:g}) with the amplitude ±10^ℓ, ℓ = {L_E[0]:g} … {L_E[-1]:g} (step "
               f"{L_E[1] - L_E[0]:g}): the fitted skeleton's outcome at every point (E: exact with the nearest kernels, ·: not), shown at "
               "every 0.06 in D and every ½ in ℓ:", ""]
        lv = [l for l in L_E if abs(l * 2 - round(l * 2)) < 1e-9]
        dv = [D for D in D_E if abs(round((D - D_E[0]) / 0.06) * 0.06 - (D - D_E[0])) < 1e-9]
        at = {(round(r['D'], 9), round(r['l'], 9), r['sg']): r for r in g}
        for sg in (1.0, -1.0):
            L_ += [f"{'Same' if sg > 0 else 'Opposite'} signs:", "", "| D/σ \\ ℓ | " + ' | '.join(f"{l:g}" for l in lv) + " |",
                   '|---|' + '---|' * len(lv)]
            for D in dv:
                L_.append(f"| {D / sig:.2f} | " + ' | '.join(('E' if at[(round(D, 9), round(l, 9), sg)]['fit_ok'] else '·')
                                                        if (round(D, 9), round(l, 9), sg) in at else ' ' for l in lv) + " |")
            L_.append("")
        thr, thr0 = {}, {}
        for sg in (1.0, -1.0):
            for l in L_E:
                col = sorted([r for r in g if r['sg'] == sg and abs(r['l'] - l) < 1e-9], key=lambda r: r['D'])
                bad = [r['D'] for r in col if not r['fit_ok']]
                thr[(sg, round(l, 9))] = (max(bad) + (D_E[1] - D_E[0])) if bad else col[0]['D']
                # without the second feature's own near ties (its single-feature margin m₀ ≤ M0_E: near a midpoint of the grid)
                bad0 = [r['D'] for r in col if not r['fit_ok'] and r['m_vor'] > M0_E]
                thr0[(sg, round(l, 9))] = (max(bad0) + (D_E[1] - D_E[0])) if bad0 else col[0]['D']
        show = [-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0]
        far = max(r['D'] for r in g if not r['fit_ok'] and r['m_vor'] > M0_E)  # the farthest failure that is not a near tie of its own
        tie_D = sorted(set(round(r['D'], 3) for r in g if not r['fit_ok'] and r['D'] > far))
        L_ += ["The fitted skeleton's threshold D*(ℓ) on this grid (the first feature on a starting point, |ℓ| ≤ 3: not (iii)'s Δ₀) — every "
               "D of the grid at or beyond it exact with E —, in units of σ: " + '; '.join(
            f"{'same' if sg > 0 else 'opposite'} signs " + ', '.join(f"ℓ = {l:g}: {thr[(sg, l)] / sig:.2f}" for l in show) for sg in (1.0, -1.0))
               + f"; over all ℓ at most {max(thr.values()) / sig:.2f}σ. Beyond {far / sig:.2f}σ every failure is the second feature's own near "
               f"tie (its single-feature margin m₀ ≤ {M0_E:g}, near a midpoint of the grid: D = {', '.join(f'{d:g}' for d in tie_D)}); without "
               "them D*(ℓ) is " + '; '.join(f"{'same' if sg > 0 else 'opposite'} signs " + ', '.join(
                   f"ℓ = {l:g}: {thr0[(sg, l)] / sig:.2f}" for l in show) for sg in (1.0, -1.0)) + f", at most {max(thr0.values()) / sig:.2f}σ."]
        ends = []                                                      # (v): the outcome at both ends of the ratio
        for sg in (1.0, -1.0):
            same_lo = same_hi = differ = 0
            for D in D_E:
                col = {round(r['l'], 9): r['fit_ok'] for r in g if r['sg'] == sg and abs(r['D'] - D) < 1e-9}
                lo = {col[k] for k in col if k <= -2 + 1e-9}
                hi = {col[k] for k in col if k >= 2 - 1e-9}
                same_lo += len(lo) == 1
                same_hi += len(hi) == 1
                differ += len(lo) == 1 and len(hi) == 1 and lo != hi
            ends.append(f"{'same' if sg > 0 else 'opposite'} signs: constant for ℓ ≤ −2 at {same_lo} and for ℓ ≥ 2 at {same_hi} of {len(D_E)} "
                        f"separations, the two ends' outcomes different at {differ}")
        L_ += ["", "The two ends of the ratio ((v)): " + '; '.join(ends) + "."]
        agree = {k: sum((r[k] > 0) == r['fit_ok'] for r in g) for k, _ in PREDICTORS if k != 'm_fit' and all(k in r for r in g)}
        L_ += ["", f"The margins' signs against the fitted skeleton on the map's {len(g)} points: " + ', '.join(
            f"{lab} {agree[k]}" for k, lab in PREDICTORS if k != 'm_fit' and all(k in r for r in g)) +
               f"; exact with E on {sum(r['fit_ok'] for r in g)}; the linear skeleton's |m_sk| > 0.01 where its sign differs from the fitted "
               f"skeleton's at {sum((r['m_skel'] > 0) != r['fit_ok'] and abs(r['m_skel']) > 0.01 for r in g)} points (|m_sk| up to "
               f"{max([abs(r['m_skel']) for r in g if (r['m_skel'] > 0) != r['fit_ok']], default=0):.3f})."]
        # (ii)'s order: the linear skeleton's margin against the fitted one's along the ratio, the tail being the smaller feature
        slopes = []
        for sg in (1.0, -1.0):
            for D in D_E:
                if not any(abs(D / sig - x) < 1e-6 for x in (2.5, 3.0, 3.5, 4.0)):
                    continue
                for side in (-1, 1):
                    col = [r for r in g if r['sg'] == sg and abs(r['D'] - D) < 1e-9 and 1 - 1e-9 <= side * r['l'] <= 3 + 1e-9
                           and r['fit_ok'] and r['m_skel'] > 0 and abs(r['m_skel'] - r['m_fit']) > 1e-13]
                    if len(col) >= 6:
                        x = np.array([abs(r['l']) for r in col]) * np.log(10)
                        y = np.log([abs(r['m_skel'] - r['m_fit']) for r in col])
                        slopes.append((D / sig, sg, side, float(-np.polyfit(x, y, 1)[0])))
        if slopes:
            L_ += ["", "(ii)'s order: |m_sk − m_fit| against the ratio of the smaller feature's amplitude to the larger's, its log–log "
                   "slope over 10⁻³ … 10⁻¹ (1: first order in the tail, 2: second): " + ', '.join(
                       f"{d:.2f}σ {'+ +' if sg > 0 else '+ −'} {'second small' if side < 0 else 'first small'} {sl:.2f}"
                       for d, sg, side, sl in slopes) + "."]
        P = [r for r in E_['paths'] if not r['skip']]
        done = [r for r in P if r['out'] != 'timeout']
        key = lambda r: (round(r['tag'][1], 9), round(r['tag'][2], 9), r['tag'][3])
        ok = sum((r['out'] == 'nearest') == at[key(r)]['fit_ok'] for r in done if key(r) in at)
        L_ += ["", f"The path at {len(P)} sampled points (|ℓ| ≤ {LMAX_EV:g}; under {T_EV} s each): it agrees with the fitted skeleton at {ok} of "
               f"{len(done)} ({len(P) - len(done)} over the limit); exact with E at {sum(r['out'] == 'nearest' for r in done)}; polishes "
               f"rejected {sum(max(r.get('rej', 0), 0) for r in done)} on {sum(r.get('rej', 0) > 0 for r in done)} paths.", ""]
    if F is not None:
        code = {'nearest': 'E', 'other': 'O', 'not': '×', 'timeout': 't'}
        L_ += ["## F. A wide range", "", f"Three features at −s, 0, s (+{F_OFF:g}; s = {', '.join(f'{x:g}' for x in F_SEPS)}: "
               f"{', '.join(f'{x / 0.12:g}σ' for x in F_SEPS)}), the amplitudes {', '.join(f'{x:g}' for x in F_MAGS)} in all six orders, "
               "the signs + + + and + − +; exits only below 10⁻¹² of the largest coefficient: the outcome per order (E exact with E, O with "
               "others, × not, t over the limit), the polishes rejected, the skeletons' signs:", "",
               "| s | signs | outcomes (orders 012, 021, 102, 120, 201, 210) | polishes rejected | the skeleton right | the fitted skeleton right |",
               "|---|---|---|---|---|---|"]
        F = fill_margins(list(F))
        for s_ in F_SEPS:
            for sg in ((1.0, 1.0, 1.0), (1.0, -1.0, 1.0)):
                rs = sorted([r for r in F if r['tag'][1] == s_ and tuple(r['tag'][2]) == sg], key=lambda r: tuple(r['tag'][3]))
                rr = [r for r in rs if not r['skip'] and r['out'] != 'timeout']
                L_.append(f"| {s_ / 0.12:g}σ | {' '.join('+' if v > 0 else '−' for v in sg)} | " + ''.join(code[r['out']] for r in rs) +
                          f" | {sum(max(r.get('rej', 0), 0) for r in rs)} | {sum((r['m_skel'] > 0) == (r['out'] == 'nearest') for r in rr)} of "
                          f"{len(rr)} | {sum((r['m_fit'] > 0) == (r['out'] == 'nearest') for r in rr)} of {len(rr)} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_transit_certificate.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'transit_certificate_raw.npz')
    md = os.path.join(OUT, 'transit_certificate.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in 'ABCDEF']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEF')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', lambda: part_D(res['A'])), ('E', part_E), ('F', part_F)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEF']))


if __name__ == '__main__':
    main()
