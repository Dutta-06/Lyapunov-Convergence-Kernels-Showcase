"""
Corollary D‴.26: perfect tracking with an adaptive prior in the real closed loop — the Floquet spectrum, uniformly in m
========================================================================================================================

The three-bump instance (σ₀ = 0.15σ_crit, the reference x_m = 2 sin t, T₀ = 2π), the nearest-atom leak toward an adaptive prior in the real
closed loop of Theorem J (src/core/adaptive_floquet.py), and the sampled weights w = (6, 5, 4)/15, representable (and a saddle of the fixed
prior, Proposition D‴.24's numerics): the lift of their exact
representation ϑ^e(w) (the amplitudes A_k/w_k) with the prior at its atoms and the plant on the reference is a T₀-periodic solution with
e ≡ 0 (Proposition D‴.24(vi)). For the weighted K₀-atom loop (one particle per class, the weights w) and the lifts m = 15 and 30 (the
counts 6, 5, 4 times 1 and 2), η = σ₀ and ∞, and k_e = 25 and 50: the return residual of one period, the monodromy matrix by central
differences of the period map (step 10⁻⁵ relative, DOP853 at rtol 10⁻¹³) and its Floquet exponents (rates −log|λ|/T₀). Checked: the
separations' exponents are σ₀/k_e exactly ((m − K₀)p of them), the lifts' other exponents are the K₀-atom loop's, and the K₀-atom loop's are
the averaged prediction
(1/k_e)λ_j^±(w) of Theorem D‴.23(ii) at the weights w (η = ∞: h_j(w)/k_e) up to the averaging correction, h_j(w) the whitened weighted
Gauss–Newton spectrum at ϑ^e(w); the correction's scaling from k_e = 25 to 50.

Part C (`--only-C`): the variational equation along the lift on the collapsed subspace (src/core/adaptive_floquet.CollapsedVariational: the
whitened weighted coordinates, the monodromy by integrating it, no differences) — checked against the real loop's linearisation and against
Part A's differences — at k_e = 12.5 … 400, for x_m = 2 sin t (where the first-order averaging correction Ā₁ vanishes identically: in one
dimension and by the symmetry t → π − t) and for the asymmetric x_m = 2 sin t + 0.8 sin 2t (Ā₁ ≠ 0, skew): the deviations of the collective
exponents from the first-order prediction λ_j^±(w)/k_e (O(k_e⁻²) relative) and from the second-order one (O(k_e⁻⁴): the odd orders vanish),
the second-order coefficients c_j and their error-filter part, and the error channel's exponent by Liouville's formula.

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


def save_raw(path, **arrays):
    """np.savez to a temporary file, then an atomic rename: an interrupted save never truncates the raw file."""
    tmp = path + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, path)
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
KES_C = (12.5, 25.0, 50.0, 100.0, 200.0, 400.0)                           # Part C
REFS_C = ('sym', 'asym')
N_WORKERS = int(os.environ.get('AFL_WORKERS', 3))


def x_ref_C(key):
    if key == 'sym':
        return lambda t: np.array([2.0 * np.sin(t)])
    return lambda t: np.array([2.0 * np.sin(t) + 0.8 * np.sin(2.0 * t)])


def x_ref_dot_C(key):
    if key == 'sym':
        return lambda t: np.array([2.0 * np.cos(t)])
    return lambda t: np.array([2.0 * np.cos(t) + 1.6 * np.cos(2.0 * t)])


REF_NAME = {'sym': '2 sin t', 'asym': '2 sin t + 0.8 sin 2t'}


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


def job_C(args):
    """Part C: one reference and one η — the variational exponents at KES_C against the first- and second-order predictions, the coefficients
    c_j, the norms of Ā₁ (and of its self-adjoint part in the metric diag(I, (σ₀/η)I)), and the variational generator against the central
    differences of the real loop's field (AdaptiveRealLoop.rhs, the weighted K₀-atom loop with the prior in the state) at six times."""
    key, eta = args
    t0 = time.time()
    net, w, exact_w, loop = setup()
    cv = AF.CollapsedVariational(net, x_ref_C(key), PER, G, SIGMA0, eta, exact_w, w)
    gens = cv.averaged()
    A0, A1, A2, H, filt = gens
    d = cv.d
    Wm = np.eye(A0.shape[0]) if np.isinf(eta) else np.diag(np.r_[np.ones(d), np.full(d, SIGMA0 / eta)])
    lam, c, cf = cv.second_order(gens=gens)
    rows = []
    for ke in KES_C:
        slow, err = cv.exponents(ke)
        rows.append(dict(ke=ke, slow=slow, err=err, p1=cv.predicted(ke, 1, gens=gens), p2=cv.predicted(ke, 2, gens=gens)))
    # the generator against the real loop's linearisation (the K₀-atom loop with the weights w; whitened coordinates)
    cl = AF.AdaptiveRealLoop(net, lambda X: f_vec(X), x_ref_C(key), x_ref_dot_C(key), PER, KES[0], G, SIGMA0, np.arange(3), 3, eta, weights=w)
    Dw = (np.sqrt(w)[:, None] / np.sqrt(G)[None, :]).ravel()
    T = np.concatenate([[1.0], Dw] + ([] if np.isinf(eta) else [Dw]))
    lin = 0.0
    for t in np.linspace(0.0, PER, 6, endpoint=False):
        z = cl.lift(exact_w); z[0] = x_ref_C(key)(t)[0]
        Dz = len(z); Jr = np.zeros((Dz, Dz))
        for j in range(Dz):
            hh = 1e-6 * max(1.0, abs(z[j])); e = np.zeros(Dz); e[j] = hh
            Jr[:, j] = (cl.rhs(t, z + e, 3) - cl.rhs(t, z - e, 3)) / (2 * hh)
        Mw = T[:, None] * Jr / T[None, :]
        M = cv.generator(t, KES[0])
        lin = max(lin, float(np.abs(Mw - M).max() / np.abs(M).max()))
    print(f"[{time.time() - t0:.0f}s] C {key} eta {eta:g}", flush=True)
    return dict(key=key, eta=eta, h=np.linalg.eigvalsh(H), lam=lam, c=c, cf=cf, rows=rows, A1=float(np.abs(A1).max()),
                A1_selfadj=float(np.abs(Wm @ A1 + (Wm @ A1).T).max()), lin=lin, trH=float(np.trace(H)))


def split(r):
    """The (m − K₀)p exponents nearest σ₀/k_e (the separations) and the others (the collective ones), both sorted."""
    s = np.sort(r['slow']); n_sep = (r['m'] - 3) * 3
    order = np.argsort(np.abs(s - SIGMA0 / r['ke']))
    sep = np.sort(s[order[:n_sep]]); coll = np.sort(s[order[n_sep:]])
    return sep, coll


def report(R, t0=None, RC=None):
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
    if RC:
        lines += ["", "## C. The variational equation: the remainder of (c), to k_e = 400", "",
                  "The variational equation along the lift on the collapsed subspace, in the whitened weighted coordinates (δe, y, ζ): δė = −k_eδe − Ĵ(t)y, ẏ = Ĵ(t)ᵀδe − "
                  "(σ₀/k_e)(y − ζ), ζ̇ = −(η/k_e)(ζ − y) (η = ∞: ζ = y); its monodromy by integrating it (DOP853, rtol 10⁻¹³; no differences), the error channel's "
                  "exponent by Liouville's formula (the exponents sum to k_e + (σ₀ + η)·9/k_e; η = ∞: k_e). Two references: x_m = 2 sin t, where the first-order "
                  "averaging correction Ā₁ vanishes identically (in one dimension ⟨ĴᵀĴ̇⟩ is the integral of an exact differential, and the commutator term vanishes "
                  "by the symmetry t → π − t), and the asymmetric x_m = 2 sin t + 0.8 sin 2t, where it does not.", ""]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            lines.append(f"- {REF_NAME[r['key']]}, η = {eta_s(r['eta'])}: ‖Ā₁‖_max = {r['A1']:.2e}, its self-adjoint part {r['A1_selfadj']:.1e} (Ā₁ is skew); "
                         f"the generator against the central differences of the real loop's field at six times: {r['lin']:.1e} (relative); h_j = "
                         + ", ".join(f"{x:.4e}" for x in r['h']) + ".")
        fdd = []
        for r in RC:
            if r['key'] != 'sym':
                continue
            for row in r['rows']:
                if (r['eta'], row['ke']) in K:
                    fdd.append(np.max(np.abs(np.sort(K[(r['eta'], row['ke'])]['slow']) - row['slow']) / row['slow']))
        if fdd:
            lines += ["", f"Against the finite differences of the first table (x_m = 2 sin t, k_e = 25 and 50, both η) the variational exponents differ by at most "
                      f"{max(fdd):.1e} (relative): the resolution of the differences.", ""]
        lines += ["| reference | η | k_e | slowest mode: (actual − first order)/first order | largest over the modes | largest \|actual − second order\|/first order | "
                  "error channel (Liouville) | k_e − tr H/k_e |", "|---|---|---|---|---|---|---|---|"]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            for row in r['rows']:
                d1 = (row['slow'] - row['p1']) / row['p1']; d2 = np.abs(row['slow'] - row['p2']) / row['p1']
                lines.append(f"| {REF_NAME[r['key']]} | {eta_s(r['eta'])} | {row['ke']:g} | {d1[0]:+.4e} | {np.max(np.abs(d1)):.3e} | {np.max(d2):.2e} | "
                             f"{row['err']:.6f} | {row['ke'] - r['trH'] / row['ke']:.6f} |")
        lines += [""]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            d1 = [np.max(np.abs(row['slow'] - row['p1']) / row['p1']) for row in r['rows']]
            d2 = [np.max(np.abs(row['slow'] - row['p2']) / row['p1']) for row in r['rows']]
            neg = all(np.all(row['slow'] < row['p1']) for row in r['rows'])
            lines.append(f"- {REF_NAME[r['key']]}, η = {eta_s(r['eta'])}: per doubling of k_e the largest first-order deviation falls by "
                         + ", ".join(f"{a / b:.2f}" for a, b in zip(d1[:-1], d1[1:])) + " (→ 4: O(k_e⁻²) relative) and the second-order residual by "
                         + ", ".join(f"{a / b:.2f}" for a, b in zip(d2[:-1], d2[1:])) + " (→ 16: O(k_e⁻⁴), the odd orders vanishing); every exponent below its first-order "
                         f"prediction at every k_e: {neg}.")
        lines += ["", "The second-order coefficients (the slow-time rates are λ_j + c_j/k_e² + O(k_e⁻⁴)), c_j/λ_j in the order of λ_j, the error filter's part of "
                  "the slowest mode's (at η = ∞ minus the Rayleigh quotient v₁ᵀ⟨Ĵ̇ᵀĴ̇⟩v₁/h₁ = Σ_n n²|Ĵ_nv₁|²/Σ_n|Ĵ_nv₁|², the mode's mean squared frequency in the "
                  "regressor), and the measured (actual·k_e − λ₁)k_e²/λ₁ of the slowest mode at k_e = 100, 200, 400:", ""]
        for r in sorted(RC, key=lambda r: (r['key'] != 'sym', np.isinf(r['eta']))):
            meas = [(row['slow'][0] * row['ke'] - r['lam'][0]) * row['ke'] ** 2 / r['lam'][0] for row in r['rows'] if row['ke'] >= 100]
            lines.append(f"- {REF_NAME[r['key']]}, η = {eta_s(r['eta'])}: c_j/λ_j = " + ", ".join(f"{x:.2f}" for x in r['c'] / r['lam'])
                         + f"; the slowest c₁/λ₁ = {r['c'][0] / r['lam'][0]:.2f}, of which the filter's {r['cf'][0] / r['lam'][0]:.2f} = −ω₁², "
                         f"ω₁ = {np.sqrt(-r['cf'][0] / r['lam'][0]):.2f} (the slowest mode's root-mean-square frequency in the regressor); measured "
                         + ", ".join(f"{x:.2f}" for x in meas) + ".")
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
            RC = list(z['RC']) if 'RC' in z.files else None
        report(R, t0, RC)
        return
    if '--only-C' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            R = list(z['R'])
        with mp.Pool(min(N_WORKERS, 4)) as pool:
            RC = pool.map(job_C, [(key, eta) for key in REFS_C for eta in ETAS], chunksize=1)
        save_raw(raw, R=np.array(R, dtype=object), RC=np.array(RC, dtype=object))
        report(R, t0, RC)
        return
    jobs = [(mult, eta, ke) for mult in MULTS for eta in ETAS for ke in KES]
    jobs.sort(key=lambda j: -j[0])                                          # the large systems first
    with mp.Pool(N_WORKERS) as pool:
        R = pool.map(job, jobs, chunksize=1)
        RC = pool.map(job_C, [(key, eta) for key in REFS_C for eta in ETAS], chunksize=1)
    save_raw(raw, R=np.array(R, dtype=object), RC=np.array(RC, dtype=object))
    report(R, t0, RC)


if __name__ == '__main__':
    main()
