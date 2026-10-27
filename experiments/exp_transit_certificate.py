"""
Theorem D‴.55: the transit's certificate — exact recovery of separated features by their nearest starting kernels iff
η(c_j) = Σ_k sign(A_k)⟨φ(c*_k)φ(c_j)⟩/⟨φ(c*_k)φ(c_{e_k})⟩ stays below one in modulus on the other starting points
====================================================================================================================================

The sliding inverse scale space (`rich_transit.sliding_iss`, `rich_transit_nd.sliding_iss_nd`: the limit of the free law with moving
centres from small masses) on random targets of K separated Gaussian features, against the certificate on the starting set (the
samples' inner product; e_k the starting point nearest to the k-th feature):

A. One dimension (D‴.48's dictionary: 21 kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]): three features at
   least DMIN_A apart, amplitudes ±U(0.4, 1): the certificate's margin m = 1 − max_{j∉E}|η(c_j)| against the path's outcome — exact with
   the nearest kernels (the active set E, the residual zero, no exits), exact with others (three kernels, the residual zero), or not.
B. The amplitudes' range: the same positions and signs with |A_k| = 10^{−r·u_k} (u_k uniform in [0, 1]) for r = R_B — the certificate
   does not change; does the outcome?
C. The plane (D‴.54's lattice: 17 × 17 kernels of width 0.15, spacing 0.125, 41² samples of [−1, 1]²): three features at least
   DMIN_C apart, amplitudes ±U(0.4, 1).

Usage: python experiments/exp_transit_certificate.py [--only A,B,C] [--from-raw] → results/chaos/transit_certificate.md (+ _raw.npz).
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


def certificate(K, c0, cs, As):
    """η at every starting point (the samples' inner product), the nearest points E, the margin m off E, and the cross-talk q (continuum,
    q = max_k Σ_{l≠k} G((Δ_kl − 2R)₊)/G(R), G(d) = e^{−d²/(4σ²)})."""
    c0 = np.asarray(c0, float).reshape(len(c0), -1)
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    d0 = np.linalg.norm(c0[:, None, :] - cs[None, :, :], axis=2)
    e = [int(np.argmin(d0[:, k])) for k in range(len(cs))]
    P0 = K.phi(c0 if c0.shape[1] > 1 else c0[:, 0])
    Pf = K.phi(cs if cs.shape[1] > 1 else cs[:, 0])
    gram = P0.T @ Pf / P0.shape[0]
    own = np.array([gram[e[k], k] for k in range(len(cs))])
    eta = (np.sign(As)[None, :] * gram / own[None, :]).sum(1)
    off = np.array([j for j in range(len(c0)) if j not in e])
    sig = K.sig
    G = lambda d: np.exp(-np.asarray(d, float) ** 2 / (4 * sig ** 2))
    R = float(max(d0[e[k], k] for k in range(len(cs))))
    D = np.linalg.norm(cs[:, None, :] - cs[None, :, :], axis=2)
    q = float(max(sum(G(max(D[k, l] - 2 * R, 0.0)) for l in range(len(cs)) if l != k) for k in range(len(cs)))) / float(G(R))
    return dict(eta=eta, e=e, m=float(1 - np.abs(eta[off]).max()), argm=int(off[np.argmax(np.abs(eta[off]))]), q=q, R=R,
                unique=len(set(e)) == len(e))


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


def job_1d(args):
    cs, As = args
    warnings.filterwarnings('ignore')
    L, K, c0 = TM.setup()
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(cs=cs, As=As, skip=True)
    R = RT.sliding_iss(L, K, X.bumps_target(K, cs, As), c0, np.ones(c0.size), polish=True)
    return dict(cs=cs, As=As, skip=False, m=cert['m'], q=cert['q'], argm=cert['argm'], e=cert['e'], out=outcome(R, cert['e'], len(cs)),
                order=[int(j) for t, j, k, c in R['events']], rms=R['hist'][-1][3])


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
                order=[int(j) for t, j, k, c in R['events']], rms=R['hist'][-1][3])


def part_A():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_A)
    jobs, tags = [], []
    for dmin in DMIN_A:
        for cs, s, a, u in targets(rng, N_A, dmin):
            jobs.append((cs, s * a)); tags.append(dmin)
    with Pool(N_WORKERS) as pool:
        rs = pool.map(job_1d, jobs, chunksize=4)
    return [dict(dmin=t, **r) for t, r in zip(tags, rs)]


def part_B():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_B)
    jobs, tags = [], []
    for q, (cs, s, a, u) in enumerate(targets(rng, N_B, DMIN_B)):
        for r in R_B:
            jobs.append((cs, s * 10.0 ** (-r * u))); tags.append((q, r))
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


def band(rs):
    """Outside the band |m| ≤ Cq: how many targets, and on how many the certificate's sign decides; the paths' ends."""
    parts = []
    for C in BANDS:
        ob = [r for r in rs if abs(r['m']) > C * r['q']]
        parts.append(f"C = {C:g}: {sum((r['m'] > 0) == (r['out'] == 'nearest') for r in ob)} of {len(ob)}")
    return ("Outside the band |m| ≤ Cq the sign decides on " + '; '.join(parts) + f". Every path ends at an exact fit with "
            f"{K_F} kernels on {sum(r['out'] != 'not' for r in rs)} of {len(rs)} targets.")


def report(A=None, B=None, C=None):
    L_ = ["# Theorem D‴.55: the transit's certificate (`exp_transit_certificate.py`)", "",
          "The sliding inverse scale space (the free law's limit from small masses, moving centres) on random targets of three Gaussian "
          "features of the kernels' width; the certificate η(c_j) = Σ_k sign(A_k)⟨φ(c*_k)φ(c_j)⟩/⟨φ(c*_k)φ(c_{e_k})⟩ on the starting points "
          "(the samples' inner product, e_k the nearest starting point to feature k, E their set), its margin m = 1 − max_{j∉E}|η(c_j)| and "
          "the cross-talk q = max_k Σ_{l≠k}G((Δ_kl − 2R)₊)/G(R). The path's outcome: exact with E (its active set E, the residual below "
          f"{RMS_TOL:g}, no exits), exact with other kernels (three kernels, the residual below {RMS_TOL:g}), or not exact.", ""]

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
               f"{max([abs(r['m']) / r['q'] for r in bad], default=0):.3f} times their cross-talk q. " + band(rs), ""]
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
        L_ += ["## B. The amplitudes' range", "",
               f"{N_B} configurations (seed {SEED_B}; positions at least {DMIN_B:g} apart, signs) with |A_k| = 10^{{−r·u_k}}, r = "
               f"{', '.join(str(r) for r in R_B)} (the amplitudes' range up to 10^{max(R_B)}); the certificate depends on the positions and "
               f"signs only ({n_ok} of {len(cfgs)} with m > 0). The outcome is the same at every range on {same} of {len(cfgs)} "
               f"configurations, while the order of the entries changes with the range on {orders}"
               + (f"; it changes on {len(flip)}: " + '; '.join(f"configuration {c} (m = {m:+.4f}): {', '.join(o)}" for c, o, m in flip)
                  if flip else '') + (f". The configurations whose outcome changes have |m|/q ≤ {max(ratio[c] for c, o, m in flip):.2f}; those "
                                      f"with |m| > 0.8q ({sum(ratio[c] > 0.8 for c in cfgs)}) keep their outcome at every range" if flip else '') + ".", ""]
    if C is not None:
        L_ += ["## C. The plane", "", f"{N_C} targets per minimal separation (seed {SEED_C}; σ = 0.15, the lattice's spacing 0.125; "
               "centres in [−0.7, 0.7]²), amplitudes ±U(0.4, 1):", ""] + table(C, 'dmin', 'Δ ≥')
        rs, pos, neg, cnt = confusion(C)
        bad = [r for r in rs if (r['m'] > 0) != (r['out'] == 'nearest')]
        L_ += ["", f"Over all {len(rs)} targets the certificate's sign predicts exactness with E on {len(rs) - len(bad)}. " + band(rs), ""]
    L_.append("Numbers are produced by `experiments/exp_transit_certificate.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'transit_certificate_raw.npz')
    md = os.path.join(OUT, 'transit_certificate.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in 'ABC']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABC')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in 'ABC']))


if __name__ == '__main__':
    main()
