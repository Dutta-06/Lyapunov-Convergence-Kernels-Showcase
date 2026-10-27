"""
Proposition D‴.46: a leak in the mirror coordinates — the free mass–amplitude law made robust, its bias a Bregman-regularised fit (a ridge
from large starting masses, a LASSO from small ones)
====================================================================================================================================

The free law of Proposition D‴.45 with the leak −σ(η_i − η_p,i) added to each particle's scalar drive (Q_i kept, η explicit):
dα_i = Gβ_ia′_i[g_i − σ(η_i − η_p,i)], dβ_i = γβ_ia_i[g_i − σ(η_i − η_p,i)] (`mirror_free.FreeLaw.real_rhs_leak`, `leak_equilibrium`).

A. Identities: V̇ = −k_ee² + ed − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] in the real loop with a disturbance, at random states; the averaged
   equilibrium's residual; its ridge limit (large starting masses) and LASSO limit (small ones) on random positive definite H.
B. The real loop under an approximation error and a disturbance (a scalar plant, six Gaussian atoms, x_m = sin t + 0.4 sin 2.7t): V, the
   largest |w| and the error over time with σ = 0, 0.01, 0.1 — against the ultimate bound (d̄²/(2k_e) + σB(w*, w_p))/min(k_e, σ).
C. Slow adaptation: the real loop's late mean of w against the averaged loop's equilibrium w_σ (H, d̄_φ the orbit's averages), from a
   lazy and from an intermediate start.
D. The regularisation paths on the orbit's Gram (twelve atoms, a 3-sparse truth): the equilibrium against the ridge (lazy start) and the
   LASSO (rich start, penalty σλ₁(β⁰)) for a range of σ, the supports, the truth's atoms kept.

Output: results/chaos/mirror_leak.md and mirror_leak_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D`).
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
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mirror_free as MF
import exp_cohesion_leak as E

OUT = os.environ.get('ML_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
KE, M6, G6, GAM6, AMAX = 5.0, 6.0, 1.0, 0.5, 1.0
XR = (lambda t: np.sin(t) + 0.4 * np.sin(2.7 * t))
XRD = (lambda t: np.cos(t) + 1.08 * np.cos(2.7 * t))
CEN6, WID6 = np.linspace(-1.2, 1.2, 6), np.full(6, 0.4)
SIGMAS_B, T_B = (0.0, 0.01, 0.1), 3000.0
DIST_B = (lambda t: 0.2 * np.sin(3.1 * t) + 0.1)                             # a disturbance with a bias
APPROX_B = 0.1                                                                  # f = Σw*φ + 0.1 sin(4x): not representable
KE_C, G_C, GAM_C, T_C = 20.0, 0.5, 0.05, 40000.0
CASES_C = (('lazy', 3.0, 5e-4), ('intermediate', 0.3, 1e-4))     # (start, β⁰, σ): equilibria 0.037 and 0.12 from w*
D12, WID12, SIGMAS_D = 12, 0.25, (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2)


def phi_of(x, cen, wid):
    return np.exp(-(np.asarray(x)[..., None] - cen) ** 2 / (2 * wid ** 2))


def orbit_gram(cen, wid, T=2 * np.pi * 10, n=200001):
    """H = ⟨φφᵀ⟩ over x_m(t), t ∈ [0, T] (the quasi-periodic orbit's average), and the time grid's φ."""
    ts = np.linspace(0, T, n)
    P = phi_of(XR(ts), cen, wid)
    return P.T @ P / n, ts, P


def part_A():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M6, G6, GAM6, AMAX)
    rng = np.random.default_rng(4)
    al0, b0 = 0.1 * rng.standard_normal(6), rng.uniform(0.5, 1.5, 6); q = L.Q(al0, b0)
    wstar = rng.uniform(-0.15, 0.15, 6); al_star = L.alpha_of_w(wstar, q)
    eta_p = L.eta(al0, q)
    ftrue = lambda x: phi_of(x, CEN6, WID6) @ wstar
    errs = []
    for sig in (0.0, 0.3):
        F = L.real_rhs_leak(CEN6, WID6, ftrue, XR, XRD, KE, q, eta_p, sig, dist=DIST_B)
        for _ in range(6):
            al = al0 + 0.3 * rng.standard_normal(6); z = np.r_[rng.uniform(-1, 1), al, np.log(L.beta(al, q))]; t = rng.uniform(0, 6)
            dz = F(t, z); hh = 1e-6
            V = lambda s: L.real_lyapunov(z + s * dz, t + s, XR, al_star, q)
            dV = (V(hh) - V(-hh)) / (2 * hh)
            e = z[0] - XR(t)
            B = lambda a1, a2: L.bregman_div(a1, a2, q)
            pred = -KE * e ** 2 + e * DIST_B(t) - sig * (B(al_star, al) + B(al, al0) - B(al_star, al0))
            errs.append(abs(dV - pred) / max(1.0, abs(pred)))
    Lz = MF.FreeLaw(12.0, 1.0, 0.1, 1.0)
    lim = []
    for rep in range(4):
        A = rng.standard_normal((30, 10)); H = A.T @ A / 30
        ws = np.zeros(10); ws[rng.choice(10, 2, replace=False)] = rng.choice([-1, 1], 2) * rng.uniform(0.15, 0.3, 2)
        for b0 in (1e2, 1e4):
            q = Lz.Q(np.zeros(10), np.full(10, b0))
            w, _ = Lz.leak_equilibrium(H, ws, q, np.zeros(10), 1e-2, 1.0)
            ridge = np.linalg.solve(H + 1e-2 * Lz.m / (Lz.G * b0 ** 2) * np.eye(10), H @ ws)
            lim.append(dict(kind='ridge', b0=b0, d=float(np.abs(w - ridge).max() / np.abs(ridge).max())))
        for b0 in (1e-4, 1e-6, 1e-8):
            q = Lz.Q(np.zeros(10), np.full(10, b0))
            kap = np.sqrt((2 * Lz.c - b0) / b0); lam1 = np.pi * Lz.A * kap / (4 * Lz.G * Lz.c)
            sig = 0.02 / lam1                                                     # the LASSO penalty σλ₁ = 0.02 at every β⁰
            w, _ = Lz.leak_equilibrium(H, ws, q, np.zeros(10), sig, 1.0)
            ls = MF.lasso(H, ws, np.full(10, sig * lam1), 1.0)
            lim.append(dict(kind='lasso', b0=b0, d=float(np.abs(w - ls).max() / np.abs(ls).max()),
                            same_support=bool(np.array_equal(np.abs(w) > 1e-3 * np.abs(w).max(), np.abs(ls) > 1e-12))))
    return dict(dV=float(max(errs)), lim=lim)


def part_B():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M6, G6, GAM6, AMAX)
    rng = np.random.default_rng(11)
    al0, b0 = 0.1 * rng.standard_normal(6), rng.uniform(0.5, 1.5, 6); q = L.Q(al0, b0)
    wstar = rng.uniform(-0.15, 0.15, 6)
    ftrue = lambda x: phi_of(x, CEN6, WID6) @ wstar + APPROX_B * np.sin(4 * x)
    eta_p = L.eta(al0, q)
    al_star = L.alpha_of_w(wstar, q)
    out = []
    for sig in SIGMAS_B:
        F = L.real_rhs_leak(CEN6, WID6, ftrue, XR, XRD, KE, q, eta_p, sig, dist=DIST_B)
        z0 = np.r_[XR(0.0) + 0.3, al0, np.log(b0)]
        ts = np.linspace(0, T_B, 3001)
        sol = solve_ivp(F, (0, T_B), z0, method='DOP853', rtol=1e-10, atol=1e-12, t_eval=ts)
        Z = sol.y
        W = np.exp(Z[7:]) * L.amp(Z[1:7])[0] / L.m
        V = np.array([L.real_lyapunov(Z[:, j], t, XR, al_star, q) for j, t in enumerate(ts)])
        e = Z[0] - XR(ts)
        qs = np.array([L.Q(Z[1:7, j], np.exp(Z[7:, j])) for j in range(len(ts))])
        quart = [slice(k * 750, (k + 1) * 750 + 1) for k in range(4)]
        out.append(dict(sig=sig, ok=bool(sol.success), wmax=[float(np.abs(W[:, s]).max()) for s in quart],
                        Vmax=[float(V[s].max()) for s in quart], erms=[float(np.sqrt(np.mean(e[s] ** 2))) for s in quart],
                        dQ=float(np.abs(qs - q).max()), V0=float(V[0]), Bp=float(L.bregman_div(al_star, al0, q))))
    return dict(rows=out, dbar=float(0.2 + 0.1 + APPROX_B), ke=KE)


def part_C():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M6, G_C, GAM_C, AMAX)
    rng = np.random.default_rng(13)
    wstar = rng.uniform(-0.15, 0.15, 6)
    ftrue = lambda x: phi_of(x, CEN6, WID6) @ wstar
    H, ts_o, P = orbit_gram(CEN6, WID6)
    rows = []
    for lab, b0v, sig in CASES_C:
        al0 = np.zeros(6); b0 = np.full(6, b0v); q = L.Q(al0, b0); eta_p = L.eta(al0, q)
        F = L.real_rhs_leak(CEN6, WID6, ftrue, XR, XRD, KE_C, q, eta_p, sig)
        z0 = np.r_[XR(0.0), al0, np.log(b0)]
        ts = np.linspace(0, T_C, 20001)
        sol = solve_ivp(F, (0, T_C), z0, method='LSODA', rtol=1e-9, atol=1e-12, t_eval=ts)
        W = np.exp(sol.y[7:]) * L.amp(sol.y[1:7])[0] / L.m
        late = W[:, ts > 0.75 * T_C].mean(axis=1)
        drift = float(np.abs(W[:, ts > 0.75 * T_C][:, -1] - W[:, ts > 0.75 * T_C][:, 0]).max())
        w_eq, _ = L.leak_equilibrium(H, wstar, q, eta_p, sig, KE_C)
        rows.append(dict(lab=lab, b0=b0v, sig=sig, ok=bool(sol.success), late=late.tolist(), weq=w_eq.tolist(), wstar=wstar.tolist(),
                         d=float(np.abs(late - w_eq).max()), d_star=float(np.abs(late - wstar).max()), drift=drift))
    return rows


def part_D():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(12.0, 1.0, 0.1, AMAX)
    cen = np.linspace(-1.4, 1.4, D12); wid = np.full(D12, WID12)
    H, _, _ = orbit_gram(cen, wid)
    rng = np.random.default_rng(17)
    ws = np.zeros(D12); idx = np.array([2, 6, 9]); ws[idx] = [0.25, -0.2, 0.3]
    out = dict(eig=np.linalg.eigvalsh(H).tolist(), truth=idx.tolist(), rows=[])
    for lab, b0 in (('lazy', 30.0), ('rich', 1e-6)):
        q = L.Q(np.zeros(D12), np.full(D12, b0))
        kap = np.sqrt(abs(2 * L.c - b0) / b0); lam1 = np.pi * L.A * kap / (4 * L.G * L.c)
        for mu in SIGMAS_D:
            sig = mu / lam1 if lab == 'rich' else mu
            w, _ = L.leak_equilibrium(H, ws, q, np.zeros(D12), sig, 1.0)
            if lab == 'rich':
                ref = MF.lasso(H, ws, np.full(D12, mu), 1.0)
            else:
                ref = np.linalg.solve(H + sig * L.m / (L.G * b0 ** 2) * np.eye(D12), H @ ws)
            sup = np.flatnonzero(np.abs(w) > 1e-3 * max(np.abs(w).max(), 1e-12))
            out['rows'].append(dict(lab=lab, mu=mu, d=float(np.abs(w - ref).max()), sup=sup.tolist(),
                                    truth_kept=bool(set(idx) <= set(sup)), extra=int(len(set(sup) - set(idx))),
                                    err=float(np.abs(w - ws).max())))
    return out


def report(A=None, B=None, C=None, D=None):
    fmt = E.fmt
    L_ = ["# Proposition D‴.46: a leak in the mirror coordinates (`exp_mirror_leak.py`)", "",
          "The free law with the leak −σ(η_i − η_p,i) on each particle's scalar drive (Q kept); η, B the mirror coordinate and Bregman "
          "divergence of Proposition D‴.45 (`mirror_free`).", ""]
    if A is not None:
        L_ += ["## A. Identities", "", f"V̇ against −k_ee² + ed − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] in the real loop (k_e = {KE:g}, "
               f"a disturbance, σ = 0 and 0.3, 12 random states): largest relative difference {fmt(A['dV'])}.", "",
               "The averaged equilibrium against its limits on random positive definite H (10 × 10, 4 draws): the ridge "
               "(1/k_e)(w − w*)ᵀH(w − w*) + σ(m/(Gβ⁰²))|w|² at large starting masses (σ = 10⁻²), the LASSO with penalty σλ₁(β⁰) = 0.02 at "
               "small ones (relative distance, largest over the draws; the supports equal):", "", "| the limit | β⁰ | distance | supports equal |",
               "|---|---|---|---|"]
        for kind in ('ridge', 'lasso'):
            for b0 in sorted({r['b0'] for r in A['lim'] if r['kind'] == kind}, reverse=True):
                rs = [r for r in A['lim'] if r['kind'] == kind and r['b0'] == b0]
                L_.append(f"| {kind} | {b0:g} | {fmt(max(r['d'] for r in rs))} | "
                          + (f"{sum(r['same_support'] for r in rs)} of {len(rs)} |" if kind == 'lasso' else "— |"))
        L_.append("")
    if B is not None:
        L_ += [f"## B. The real loop under an approximation error ({APPROX_B:g} sin 4x) and a disturbance (0.2 sin 3.1t + 0.1), {T_B:g} time "
               f"units (k_e = {KE:g}, m = {M6:g}, G = {G6:g}, γ = {GAM6:g}; six atoms of width 0.4)", "",
               "By quarters: the largest |w|, the largest V (V = ½e² + ΣB(w*, w), w* the representable part), the rms error; the largest "
               f"|Q − Q(0)|; and the ultimate bound max(V(0), (d̄²/(2k_e) + σB(w*, w_p))/min(k_e, σ)) with d̄ = {B['dbar']:g} (the disturbance "
               "and the approximation error):", "",
               "| σ | max |w| by quarters | max V by quarters | rms e by quarters | max |Q − Q(0)| | the bound |", "|---|---|---|---|---|---|"]
        for r in B['rows']:
            bound = (max(r['V0'], (B['dbar'] ** 2 / (2 * B['ke']) + r['sig'] * r['Bp']) / min(B['ke'], r['sig'])) if r['sig'] > 0 else float('nan'))
            L_.append(f"| {r['sig']:g} | " + ', '.join(f"{x:.3f}" for x in r['wmax']) + " | " + ', '.join(fmt(x) for x in r['Vmax'])
                      + " | " + ', '.join(f"{x:.4f}" for x in r['erms']) + f" | {fmt(r['dQ'])} | " + (fmt(bound) if r['sig'] > 0 else "—") + " |")
        L_.append("")
    if C is not None:
        L_ += [f"## C. Slow adaptation (k_e = {KE_C:g}, G = {G_C:g}, γ = {GAM_C:g}, {T_C:g} time units, from α = 0 and equal masses β⁰, the leak "
               "toward the start): the late mean of w (the last quarter) against the averaged loop's equilibrium w_σ", "",
               "| start | β⁰ | σ | max |w̄ − w_σ| | max |w_σ − w*| | max |w̄ − w*| | the last quarter's drift |", "|---|---|---|---|---|---|---|"]
        for r in C:
            L_.append(f"| {r['lab']} | {r['b0']:g} | {r['sig']:g} | {fmt(r['d'])} | {fmt(float(np.abs(np.array(r['weq']) - np.array(r['wstar'])).max()))} | "
                      f"{fmt(r['d_star'])} | {fmt(r['drift'])} |")
        L_.append("")
    if D is not None:
        L_ += [f"## D. The regularisation paths on the orbit's Gram ({D12} atoms of width {WID12:g} on [−1.4, 1.4]; its eigenvalues "
               f"{fmt(min(D['eig']))} to {fmt(max(D['eig']))}; a 3-sparse truth on atoms {D['truth']})", "",
               "The equilibrium against the ridge (lazy start, β⁰ = 30) and the LASSO of penalty μ = σλ₁(β⁰) (rich start, β⁰ = 10⁻⁶); its "
               "support (|w| > 10⁻³ max|w|), whether it keeps the truth's atoms, the atoms beyond them, max |w − w*|:", "",
               "| start | μ (rich) or σ (lazy) | distance to the reference | support | the truth kept | extra atoms | max |w − w*| |",
               "|---|---|---|---|---|---|---|"]
        for r in D['rows']:
            L_.append(f"| {r['lab']} | {r['mu']:g} | {fmt(r['d'])} | {r['sup']} | {'yes' if r['truth_kept'] else 'no'} | {r['extra']} | {fmt(r['err'])} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_mirror_leak.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_leak_raw.npz')
    md = os.path.join(OUT, 'mirror_leak.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCD']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCD')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('D', part_D), ('B', part_B), ('C', part_C)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCD']))


if __name__ == '__main__':
    main()
