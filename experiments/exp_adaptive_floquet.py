"""
Corollary D‴.26: perfect tracking with an adaptive prior in the real closed loop — the Floquet spectrum, uniformly in m
========================================================================================================================

The three-bump instance (σ₀ = 0.15σ_crit, the reference x_m = 2 sin t, T₀ = 2π), the nearest-atom leak toward an adaptive prior in the real
closed loop of Theorem J (src/core/adaptive_floquet.py), and the sampled weights w = (6, 5, 4)/15, admissible: the lift of their exact
representation ϑ^e(w) (the amplitudes A_k/w_k) with the prior at its atoms and the plant on the reference is a T₀-periodic solution with
e ≡ 0 (Proposition D‴.24(vi)). For the weighted K₀-atom loop (one particle per class, the weights w) and the lifts m = 15 and 30 (the
counts 6, 5, 4 times 1 and 2), η = σ₀ and ∞, and k_e = 25 and 50: the return residual of one period, the monodromy matrix by central
differences of the period map (step 10⁻⁵ relative, DOP853 at rtol 10⁻¹³) and its Floquet exponents (rates −log|λ|/T₀). Checked: the
separations' exponents are σ₀/k_e exactly ((m − K₀)p of them), the lifts' other exponents are the K₀-atom loop's, and the K₀-atom loop's are
the averaged prediction
(1/k_e)λ_j^±(w) of Theorem D‴.23(ii) at the weights w (η = ∞: h_j(w)/k_e) up to the averaging correction, h_j(w) the whitened weighted
Gauss–Newton spectrum at ϑ^e(w); the correction's scaling from k_e = 25 to 50.

Output: results/chaos/adaptive_floquet.md and adaptive_floquet_raw.npz (AFL_WORKERS; `--from-raw` regenerates the md).
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
from src.core import ot_leak as OT
from src.core import adaptive_prior as AP
from src.core import adaptive_floquet as AF
from src.core.periodic import floquet_exponents

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
COUNTS = np.array([6, 5, 4])
MULTS = (0, 1, 2)                                                          # 0: the weighted K₀-atom loop
ETAS = (SIGMA0, np.inf)
KES = (25.0, 50.0)
N_WORKERS = int(os.environ.get('AFL_WORKERS', 3))


def f_vec(X):
    X = np.asarray(X, float)
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def setup():
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    w = COUNTS / COUNTS.sum()
    A = np.array([b[0] for b in BUMPS]); C = np.array([b[1] for b in BUMPS]); S = np.array([b[2] for b in BUMPS])
    exact_w = net.unconstrained((A / w)[:, None], C[:, None], S)
    q_exact = net.unconstrained((3 * A)[:, None], C[:, None], S)
    loop = OT.OTPriorLoop(net, lambda X: f_vec(X), lambda t: np.array([2.0 * np.sin(t)]), PER, G, SIGMA0, q_exact, n_quad=256)
    return net, w, exact_w, loop


def job(args):
    mult, eta, ke = args
    t0 = time.time()
    net, w, exact_w, loop = setup()
    if mult == 0:
        labels, weights = np.arange(3), w
    else:
        labels, weights = np.repeat(np.arange(3), COUNTS * mult), None
    m = len(labels)
    cl = AF.AdaptiveRealLoop(net, lambda X: f_vec(X), lambda t: np.array([2.0 * np.sin(t)]), lambda t: np.array([2.0 * np.cos(t)]), PER, ke, G,
                             SIGMA0, labels, 3, eta, weights=weights)
    z_p = cl.lift(exact_w)
    residual = float(np.abs(AF.period_map(cl, z_p, m) - z_p).max())
    DP = AF.monodromy(cl, z_p, m)
    fast, slow = floquet_exponents(DP, PER, 1)
    h = AP.gauss_newton_block(loop, exact_w, w)[1]
    pred = AF.predicted_slow_exponents(h, SIGMA0, eta, ke, m, 3, 3)
    print(f"[{time.time() - t0:.0f}s] m {m} eta {eta:g} ke {ke:g}", flush=True)
    return dict(mult=mult, m=m, eta=eta, ke=ke, residual=residual, fast=fast, slow=slow, pred=pred, h=h)


def split(r):
    """The (m − K₀)p exponents nearest σ₀/k_e (the separations) and the others (the collective ones), both sorted."""
    s = np.sort(r['slow']); n_sep = (r['m'] - 3) * 3
    order = np.argsort(np.abs(s - SIGMA0 / r['ke']))
    sep = np.sort(s[order[:n_sep]]); coll = np.sort(s[order[n_sep:]])
    return sep, coll


def report(R, t0=None):
    eta_s = lambda e: '∞' if np.isinf(e) else f'{e:.5f}'
    net, w, exact_w, loop = setup()
    lines = ["# Corollary D‴.26: perfect tracking with an adaptive prior in the real loop — the Floquet spectrum (`exp_adaptive_floquet.py`)", "",
             f"The three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f}, x_m = 2 sin t, the nearest-atom leak toward an adaptive prior in the real closed loop; the sampled "
             f"weights w = {tuple(int(c) for c in COUNTS)}/{int(COUNTS.sum())} (the amplitudes A_k/w_k = " + ", ".join(f"{b[0] / wk:.3f}" for b, wk in zip(BUMPS, w)) + "), the "
             "lift of their exact representation with the prior at its atoms and the plant on the reference. Floquet exponents: rates −log|λ|/T₀ of the monodromy matrix (central "
             "differences of the period map, step 10⁻⁵ relative, DOP853 at rtol 10⁻¹³); the error channel's multiplier (≈ e^{−k_eT₀}) is below the resolution of the differences and is left out.", "",
             "The whitened weighted Gauss–Newton spectrum at ϑ^e(w): h_j(w) = " + ", ".join(f"{x:.4e}" for x in R[0]['h']) + ".", "",
             "## The weighted K₀-atom loop against the averaged prediction", "",
             "| η | k_e | return residual | collective exponents × k_e (actual / predicted (1/k_e)λ_j^±(w) × k_e) | largest relative deviation |", "|---|---|---|---|---|"]
    K = {(r['eta'], r['ke']): r for r in R if r['mult'] == 0}
    for (eta, ke), r in sorted(K.items(), key=lambda t: (np.isinf(t[0][0]), t[0][1])):
        a = np.sort(r['slow']); p = np.sort(r['pred'])
        lines.append(f"| {eta_s(eta)} | {ke:g} | {r['residual']:.1e} | " + ", ".join(f"{x * ke:.4e}/{y * ke:.4e}" for x, y in zip(a, p)) + f" | {np.max(np.abs(a - p) / p):.2e} |")
    lines += [""]
    for eta in ETAS:
        if (eta, KES[0]) in K and (eta, KES[1]) in K:
            d = [np.max(np.abs(np.sort(K[(eta, ke)]['slow']) - np.sort(K[(eta, ke)]['pred'])) * ke) for ke in KES]
            rl = [np.max(np.abs(np.sort(K[(eta, ke)]['slow']) - np.sort(K[(eta, ke)]['pred'])) / np.sort(K[(eta, ke)]['pred'])) for ke in KES]
            lines.append(f"η = {eta_s(eta)}: the largest deviation of the slow-time rates (exponent × k_e) from λ_j^±(w) is {d[0]:.2e} at k_e = {KES[0]:g} and {d[1]:.2e} at "
                         f"k_e = {KES[1]:g} (ratio {d[0] / d[1]:.2f}, i.e. k_e^{{−{np.log(d[0] / d[1]) / np.log(KES[1] / KES[0]):.2f}}}; relative {rl[0]:.2e} and {rl[1]:.2e}) — the "
                         "averaging correction is of second order in 1/k_e in slow time, so the exponents are (1/k_e)λ_j^±(w) + O(k_e⁻³) (Corollary D‴.26(c)).")
    lines += ["", "## The lifts: separations at σ₀/k_e, the collective exponents the K₀-atom loop's", "",
              "| m | η | k_e | return residual | separations: count, largest relative distance from σ₀/k_e | collective exponents: count, largest relative difference from the K₀-atom loop's |",
              "|---|---|---|---|---|---|"]
    for r in sorted([r for r in R if r['mult'] > 0], key=lambda r: (np.isinf(r['eta']), r['ke'], r['m'])):
        sep, coll = split(r)
        ref = np.sort(K[(r['eta'], r['ke'])]['slow'])
        lines.append(f"| {r['m']} | {eta_s(r['eta'])} | {r['ke']:g} | {r['residual']:.1e} | {len(sep)}, {np.max(np.abs(sep * r['ke'] / SIGMA0 - 1)):.1e} | "
                     f"{len(coll)}, {np.max(np.abs(coll - ref) / ref):.1e} |")
    lines += ["", "Numbers are produced by `experiments/exp_adaptive_floquet.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'adaptive_floquet.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'adaptive_floquet_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            R = list(z['R'])
        report(R, t0)
        return
    jobs = [(mult, eta, ke) for mult in MULTS for eta in ETAS for ke in KES]
    jobs.sort(key=lambda j: -j[0])                                          # the large systems first
    with mp.Pool(N_WORKERS) as pool:
        R = pool.map(job, jobs, chunksize=1)
    np.savez(raw, R=np.array(R, dtype=object))
    report(R, t0)


if __name__ == '__main__':
    main()
