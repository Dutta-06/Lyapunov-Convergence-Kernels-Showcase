"""
Proposition D‴.45: the free mass–amplitude law is a mirror flow — a conservation law per particle, the implicit bias as a Bregman
projection, the noise's invariance in continuous time, the sampled-data drift, the real loop's Lyapunov function
====================================================================================================================================

The replicator–amplitude law of Proposition D‴.34 with free masses (λ ≡ 0), each particle driven by one scalar signal g_i:
dα_i = Gβ_ia′_ig_i dt, dβ_i = γβ_ia_ig_i dt; Q_i = β_i − c cosh(2α_i/a_max) (c = γa_max²/(4G)) is conserved for every drive, and on a
dictionary w = β⊙a/m is the mirror flow dη = Φᵀ(r + ξ)dt/N of the potential Ψ_q, ψ′ = η = ∫dα/(Gβa′), ψ* = (β − β⁰)/(γm)
(`mirror_free`).

A. Identities: η's closed form against quadrature (the three cases and the branch); the conjugate and the metric; Q along the loop,
   noiseless and driven by a smooth node-correlated coloured noise; the conserved law's shift dQ_i = −γβ_iλ dt; Q's drift along
   Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's update) against the Euler step's closed form.
B. The implicit bias: from α = 0 and equal masses s·m/d, the noiseless loop's end against the Bregman projection, and the projection
   against the minimum-norm fit (the lazy regime) and basis pursuit (the rich regime), for s = 10 … 10⁻⁴ (N × d = 5 × 12, 2- and 3-sparse
   truths); and from small random amplitudes at tiny masses (the branches q < −c: the signs frozen).
C. Noise: (1) in continuous time a smooth noise of large amplitude, correlated over the nodes and coloured, moves neither Q nor the state's
   projection; (2) the Euler step's mean change of Q against the closed form; (3) along the paradigm's reduced flow (D‴.44's selection,
   free masses, white noise) from the lazy start: Q decreases and the state moves toward basis pursuit; (4) the sampled-data loop with
   label noise: the paradigm's Euler update moves Q and the projection, the mirror update on the curves moves neither.
D. The real loop (a scalar plant, Gaussian atoms at fixed centres, two particles on one atom, x_m = sin t + 0.4 sin 2.7t): V = ½e² + Σ Bregman divergences monotone
   with V̇ = −k_ee²; the pair's mirror coordinates move together (η_i − η_i⁰ = η_j − η_j⁰), with and without a smooth sensor noise; the
   pair's split at the end against the Bregman split.
E. The conserved law (Σβ = m) from the same start: not a mirror flow — Q moves, and its end differs from the free law's projection.

Output: results/chaos/mirror_free.md and mirror_free_raw.npz (MF_WORKERS; `--from-raw` regenerates the md; `--only A,B,C,D,E`).
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
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mirror_free as MF
from src.core import mass_dictionary as MD
from src.core import colour_selection as CS
from src.core.runcache import cached
import exp_cohesion_leak as E

OUT = os.environ.get('MF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('MF_WORKERS', 2))
N_I, D_I, M_I, G_I, GAM_I, AMAX = 5, 12, 12.0, 1.0, 0.1, 1.0
SCALES_B, KS_B, N_INST_B, T_B = (10.0, 1.0, 0.1, 1e-2, 1e-3, 1e-4), (2, 3), 20, 1e7
EPS_RICH, SCALES_LAZY = (1e-3, 1e-5, 1e-7), (1e2, 1e4)                     # unequal masses: ε·r and scale·r
FLOW_MIN_B = 1e-2                                                           # the flow is integrated for s ≥ this (slower below)
S_C3, N_INST_C3 = (0.0, 1e5, 1e6, 1e7, 1e8), 6
KE_D, GAM_D, M_D, T_D = 5.0, 0.5, 6.0, 2000.0
H_C4, SIG_C4, N_C4, P_C4 = 0.5, 0.3, 20000, 16


def law():
    return MF.FreeLaw(M_I, G_I, GAM_I, AMAX)


def instance(k, i):
    """A 5 × 12 Gaussian design with unit columns and a k-sparse truth (signs ±1, magnitudes in [0.15, 0.3])."""
    rng = np.random.default_rng(7000 + 100 * k + i)
    Phi = rng.standard_normal((N_I, D_I)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(D_I); idx = rng.choice(D_I, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.15, 0.3, k)
    return Phi, ws


def smooth_noise(rng, N, K=30, amp=0.3, lo=1e-3, hi=3.0):
    """A smooth noise on the nodes: K sinusoids of log-uniform frequencies in [lo, hi], mixed over the nodes by a random matrix."""
    om = np.exp(rng.uniform(np.log(lo), np.log(hi), K)); ph = rng.uniform(0, 2 * np.pi, K); Mx = amp * rng.standard_normal((N, K))
    return lambda t: Mx @ np.sin(om * t + ph)


def flow_end(L, Phi, y, z0, T, conserved=False, xi=None):
    sol = solve_ivp(L.rhs(Phi, y, xi, conserved), (0.0, T), z0, method='LSODA', rtol=1e-10, atol=1e-14)
    return sol.y[:, -1]


def part_A():
    warnings.filterwarnings('ignore')
    L = law()
    f = lambda x, q: 1.0 / (L.G * L.beta(x, q) * L.amp(x)[1])
    qs = (100 * L.c, 3.0, 0.5, L.c, 0.5 * L.c, 0.0, -0.5 * L.c)
    eta_err = 0.0
    for q in qs:
        for al in np.linspace(-4, 4, 17):
            ref = quad(lambda x: f(x, q), 0.0, al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
            eta_err = max(eta_err, abs(L.eta(al, q)[0] - ref) / max(1.0, abs(ref)))
    qb = -2 * L.c; a0 = float(L.branch_point(qb))                              # the branch, from a reference point
    for al in a0 + np.array([0.2, 0.5, 1.0, 2.0, 3.5]):
        ref = quad(lambda x: f(x, qb), a0 + 0.1, al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
        eta_err = max(eta_err, abs((L.eta(al, qb) - L.eta(a0 + 0.1, qb))[0] - ref) / max(1.0, abs(ref)))
    conj_err = metric_err = 0.0
    for q in (0.97, 0.01, -0.02):
        for al in (-1.3, 0.2, 2.2):
            hh = 1e-5
            dw = (L.w(al + hh, q) - L.w(al - hh, q)) / (2 * hh)
            de = (L.eta(al + hh, q) - L.eta(al - hh, q))[0] / (2 * hh)
            dps = (L.psi_star(al + hh, q) - L.psi_star(al - hh, q)) / (2 * hh)
            conj_err = max(conj_err, abs(dps / de - L.w(al, q)) / max(abs(L.w(al, q)), 1e-3))
            metric_err = max(metric_err, abs(dw / de - L.m * L.relevance(al, q)) / (L.m * L.relevance(al, q)))
    cons = []
    for i in range(3):
        Phi, ws = instance(2, i); y = Phi @ ws; d = D_I
        rng = np.random.default_rng(300 + i)
        z0 = np.r_[0.05 * rng.standard_normal(d), np.log(np.full(d, 0.3))]
        q0 = L.Q(z0[:d], np.exp(z0[d:]))
        out = {}
        for lab, xi, T in (('noiseless', None, 1e4), ('smooth noise', smooth_noise(rng, N_I), 3e4)):
            sol = solve_ivp(L.rhs(Phi, y, xi), (0.0, T), z0, method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
            ts = np.linspace(0, T, 60)
            out[lab] = float(max(np.abs(L.Q(sol.sol(t)[:d], np.exp(sol.sol(t)[d:])) - q0).max() for t in ts))
        Fc = L.rhs(Phi, y, smooth_noise(np.random.default_rng(9), N_I), conserved=True)
        z = sol.y[:, -1]; dz = Fc(1.0, z); al, b = z[:d], np.exp(z[d:])
        dQ = b * dz[d:] - L.c * (2 / L.A) * np.sinh(2 * al / L.A) * dz[:d]
        a = L.amp(al)[0]
        xi9 = smooth_noise(np.random.default_rng(9), N_I)
        g = Phi.T @ (y - Phi @ (b * a / L.m) + xi9(1.0)) / N_I
        lam = b @ (a * g) / b.sum()
        out['shift'] = float(np.abs(dQ + L.gamma * b * lam).max() / np.abs(L.gamma * b * lam).max())
        cons.append(out)
    # Q's drift along D‴.44's reduced flow (free masses, white noise, the paradigm's update) against the Euler step's closed form
    red = []
    for i in range(3):
        Phi, ws = instance(2, i); y = Phi @ ws
        X = MD.MassDictionary(Phi, y, M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
        al0 = np.zeros(D_I); q = L.Q(al0, np.full(D_I, M_I / D_I))
        _, al, _, _ = L.bregman_projection(Phi, y, q, L.eta(al0, q))
        th = X.theta(al, L.beta(al, q))
        dr = X.drift(th, lambda k: np.ones_like(k), ito=True)
        b = L.beta(al, q); a, da = L.amp(al)
        dQ = -L.c * (2 / L.A) * np.sinh(2 * al / L.A) * dr[:D_I] + 2 * np.sqrt(b) * dr[D_I:]
        pred = (2 * (Phi ** 2).sum(axis=0) / (L.m * N_I)) * L.gamma * b * (L.gamma * a ** 2 - L.G * b * da * (1 + (a / L.A) ** 2))
        red.append(float(np.abs(dQ - pred).max() / np.abs(pred).max()))
    return dict(eta_err=eta_err, conj_err=conj_err, metric_err=metric_err, cons=cons, red=red)


@cached
def b_job(args):
    """Part B, one (k, instance): the projection for every start scale, the flow's end for s ≥ FLOW_MIN_B; the branches' variant."""
    warnings.filterwarnings('ignore')
    k, i = args
    t0 = time.time()
    L = law()
    Phi, ws = instance(k, i); y = Phi @ ws
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    wl2 = np.linalg.pinv(Phi) @ y
    out = dict(k=k, i=i, bp_ok=bool(np.abs(wbp - ws).max() < 1e-4 * np.abs(ws).max()), l1_bp=float(np.abs(wbp).sum()))
    for s in SCALES_B:
        al0 = np.zeros(D_I); b0 = np.full(D_I, s * M_I / D_I)
        q = L.Q(al0, b0)
        w, al, it, res = L.bregman_projection(Phi, y, q, L.eta(al0, q))
        r = dict(l1=float(np.abs(w).sum()), used=int((np.abs(w) > 1e-3).sum()), res=res, it=it,
                 d_bp=float(np.abs(w - wbp).max()), d_l2=float(np.abs(w - wl2).max() / np.abs(wl2).max()),
                 ok=bool(np.abs(w - ws).max() < 1e-3 * np.abs(ws).max()), amax=float(np.abs(L.amp(al)[0]).max() / AMAX),
                 bmax=float(L.beta(al, q).max()))
        if s >= FLOW_MIN_B:
            z = flow_end(L, Phi, y, np.r_[al0, np.log(b0)], T_B)
            wT = np.exp(z[D_I:]) * L.amp(z[:D_I])[0] / L.m
            r.update(flow_res=float(np.linalg.norm(Phi @ wT - y)), flow_d=float(np.abs(wT - w).max()))
        out[s] = r
    # the branches: small random amplitudes at tiny masses
    for s in (1e-2, 1e-4):
        al0 = 0.05 * np.random.default_rng(900 + 100 * k + i).standard_normal(D_I); b0 = np.full(D_I, s * M_I / D_I)
        q = L.Q(al0, b0); br = q + L.c < 0
        w, al, it, res = L.bregman_projection(Phi, y, q, L.eta(al0, q), sgn=np.sign(al0))
        out[('branch', s)] = dict(n_branch=int(br.sum()), used=int((np.abs(w) > 1e-3).sum()),
                                  used_br=int(((np.abs(w) > 1e-3) & br).sum()), res=res,
                                  signs_kept=bool(np.all((np.sign(w[br]) == np.sign(al0[br])) | (np.abs(w[br]) < 1e-12))),
                                  d_bp=float(np.abs(w - wbp).max()))
    # unequal masses (ratios r in [1/4, 4]): the rich limit against weighted basis pursuit (weights 1/√β⁰), the lazy one against the
    # weighted minimum norm (weights 1/β⁰²)
    from scipy.optimize import linprog
    r = np.exp(np.random.default_rng(800 + 100 * k + i).uniform(np.log(0.25), np.log(4.0), D_I))
    cw = 1.0 / np.sqrt(r)
    lp = linprog(np.r_[cw, cw], A_eq=np.c_[Phi, -Phi], b_eq=y, bounds=(0, None), method='highs')
    wl1 = lp.x[:D_I] - lp.x[D_I:]
    for eps in EPS_RICH:
        q = L.Q(np.zeros(D_I), eps * r)
        w = L.bregman_projection(Phi, y, q, L.eta(np.zeros(D_I), q))[0]
        out[('rich', eps)] = dict(d=float(np.abs(w - wl1).max()), d_unw=float(np.abs(w - wbp).max()))
    for big in SCALES_LAZY:
        b0 = big * r; q = L.Q(np.zeros(D_I), b0)
        w = L.bregman_projection(Phi, y, q, L.eta(np.zeros(D_I), q))[0]
        wl2 = b0 ** 2 * (Phi.T @ np.linalg.solve(Phi @ ((b0 ** 2)[:, None] * Phi.T), y))
        out[('lazy', big)] = dict(d=float(np.abs(w - wl2).max() / np.abs(wl2).max()))
    print(f"[{time.time() - t0:.0f}s] B k {k} instance {i}", flush=True)
    return out


def part_C12():
    """C1: the smooth noise in continuous time; C2: the Euler step's mean change of Q."""
    warnings.filterwarnings('ignore')
    L = law()
    c1 = []
    for i in range(4):
        Phi, ws = instance(2, i); y = Phi @ ws; d = D_I
        rng = np.random.default_rng(400 + i)
        al0 = np.zeros(d); b0 = np.full(d, 0.1 * M_I / d); z0 = np.r_[al0, np.log(b0)]
        q = L.Q(al0, b0)
        w_ref = L.bregman_projection(Phi, y, q, L.eta(al0, q))[0]
        xi = smooth_noise(rng, N_I)
        sol = solve_ivp(L.rhs(Phi, y, xi), (0.0, 3e4), z0, method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
        rows = []
        for t in (1e3, 1e4, 3e4):
            z = sol.sol(t); al = z[:d]
            wi = L.bregman_projection(Phi, y, q, L.eta(al, q), sgn=np.sign(al))[0]
            rows.append(dict(t=t, dQ=float(np.abs(L.Q(al, np.exp(z[d:])) - q).max()),
                             res=float(np.linalg.norm(Phi @ (np.exp(z[d:]) * L.amp(al)[0] / L.m) - y)), d_int=float(np.abs(wi - w_ref).max())))
        c1.append(rows)
    c2 = []
    for h in (0.5, 0.05):
        for i in range(3):
            Phi, ws = instance(2, i); y = Phi @ ws; d = D_I
            rng = np.random.default_rng(500 + i)
            al = 0.4 * rng.standard_normal(d); b = rng.uniform(0.2, 2.0, d)
            sig, P = 0.3, 400000
            a, da = L.amp(al)
            r0 = y - Phi @ (b * a / L.m)
            g = (Phi.T @ (r0[:, None] + sig * rng.standard_normal((N_I, P)))) / N_I
            dQ = L.dQ_exact(al[:, None], b[:, None], g, h)
            Eg2 = (Phi.T @ r0 / N_I) ** 2 + sig ** 2 * (Phi ** 2).sum(axis=0) / N_I ** 2
            pred = L.dQ_euler(al, b, Eg2, h)
            se = dQ.std(axis=1) / np.sqrt(P)
            c2.append(dict(h=h, rel=float(np.abs(dQ.mean(axis=1) - pred).max() / np.abs(pred).max()),
                           z=float(np.abs((dQ.mean(axis=1) - pred) / se).max())))
    return dict(c1=c1, c2=c2)


def part_C4():
    """C4: the sampled-data loop at the step H_C4 with label noise drawn afresh (σ = SIG_C4, P_C4 paths, N_C4 steps) from the lazy
    projection (α = 0, β = m/d): the paradigm's Euler update on (α, log β) against the mirror update η⁺ = η + hg on each particle's curve."""
    warnings.filterwarnings('ignore')
    L = law()
    rows = []
    for i in range(2):
        Phi, ws = instance(2, i); y = Phi @ ws
        wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
        al00 = np.zeros(D_I); q0 = L.Q(al00, np.full(D_I, M_I / D_I))
        w0, al0, _, _ = L.bregman_projection(Phi, y, q0, L.eta(al00, q0))
        rng = np.random.default_rng(600 + i)
        out = dict(i=i, d_bp0=float(np.abs(w0 - wbp).max()), l1_0=float(np.abs(w0).sum()))
        for upd in ('Euler', 'mirror'):
            al = np.tile(al0, (P_C4, 1)); b = np.tile(L.beta(al0, q0), (P_C4, 1))
            for n in range(N_C4):
                a, da = L.amp(al)
                r = y[None, :] - (b * a / L.m) @ Phi.T + SIG_C4 * rng.standard_normal((P_C4, N_I))
                g = r @ Phi / N_I
                if upd == 'Euler':
                    al, b = al + H_C4 * L.G * b * da * g, b * np.exp(H_C4 * L.gamma * a * g)
                else:
                    al = L.mirror_step(al, q0, g, H_C4); b = L.beta(al, q0)
            Qp = L.Q(al, b)
            proj = []
            for p_ in range(P_C4):
                wp = L.bregman_projection(Phi, y, Qp[p_], L.eta(al[p_], Qp[p_]), sgn=np.sign(al[p_]))[0]
                proj.append(wp)
            proj = np.array(proj)
            out[upd] = dict(dQ_mean=float((Qp - q0).mean()), dQ_max=float(np.abs(Qp - q0).max()),
                            d_proj=float(np.abs(proj - w0).max(axis=1).mean()), d_bp=float(np.abs(proj - wbp).max(axis=1).mean()),
                            l1=float(np.abs(proj).sum(axis=1).mean()), res=float(np.linalg.norm((b * L.amp(al)[0] / L.m) @ Phi.T - y, axis=1).mean()))
        rows.append(out)
    return rows


@cached
def c3_job(i):
    """Part C3, one instance: D‴.44's reduced flow (free masses, white noise, the paradigm's update) from the lazy projection (s = 1)."""
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law()
    Phi, ws = instance(2, i); y = Phi @ ws
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    X = MD.MassDictionary(Phi, y, M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
    al0 = np.zeros(D_I); q0 = L.Q(al0, np.full(D_I, M_I / D_I))
    _, al, _, _ = L.bregman_projection(Phi, y, q0, L.eta(al0, q0))
    th0 = X.theta(al, L.beta(al, q0))
    _, Th = X.flow(th0, S_C3[-1], lambda k: np.ones_like(k), ito=True, t_eval=np.array(S_C3), chunks=60)
    rows = []
    for th in Th:
        alx, _, b, a, _, _ = X.parts(th); w = X.w(th); Qx = L.Q(alx, b)
        rows.append(dict(l1=float(np.abs(w).sum()), d_bp=float(np.abs(w - wbp).max()), used=int((np.abs(w) > 1e-3).sum()),
                         qmean=float(Qx.mean()), qmin=float(Qx.min()), qmax=float(Qx.max())))
    print(f"[{time.time() - t0:.0f}s] C3 instance {i}", flush=True)
    return dict(i=i, rows=rows, bp_ok=bool(np.abs(wbp - ws).max() < 1e-4 * np.abs(ws).max()))


def part_D():
    """The real loop: Gaussian atoms at fixed centres (six, width 0.4) plus a second particle on the third atom; x_m = sin t."""
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M_D, G_I, GAM_D, AMAX)
    rng = np.random.default_rng(11)
    cen6 = np.linspace(-1.2, 1.2, 6)
    cen = np.r_[cen6, cen6[2]]; wid = np.full(7, 0.4)
    al0 = 0.1 * rng.standard_normal(7); b0 = rng.uniform(0.5, 1.5, 7)
    q = L.Q(al0, b0)
    wtrue6 = rng.uniform(-0.15, 0.15, 6)
    # the pair's split: the Bregman split of wtrue6[2] given the start (η_i − η_i⁰ = η_j − η_j⁰ and w_i + w_j = wtrue6[2])
    from scipy.optimize import brentq
    e0 = L.eta(al0, q)
    pair_w = lambda dl: L.w(L.alpha_of_eta(np.array([e0[2] + dl]), q[2:3])[0], q[2]) + \
        L.w(L.alpha_of_eta(np.array([e0[6] + dl]), q[6:7])[0], q[6])
    dl = brentq(lambda x: pair_w(x) - wtrue6[2], -200, 200, xtol=1e-14)
    wi = L.w(L.alpha_of_eta(np.array([e0[2] + dl]), q[2:3])[0], q[2]); wj = wtrue6[2] - wi
    wstar = np.r_[wtrue6[:2], wi, wtrue6[3:], wj]                     # the target's coefficients, the pair at its Bregman split
    al_star = L.alpha_of_w(wstar, q)
    xr, xrd = (lambda t: np.sin(t) + 0.4 * np.sin(2.7 * t)), (lambda t: np.cos(t) + 1.08 * np.cos(2.7 * t))
    out = {}
    for lab, noise in (('noiseless', None), ('smooth sensor noise', smooth_noise(np.random.default_rng(12), 1, K=20, amp=0.1, lo=0.05, hi=20.0))):
        F0 = L.real_rhs(cen, wid, wstar, xr, xrd, KE_D)
        if noise is None:
            F = F0
        else:
            def F(t, z, F0=F0, noise=noise):                                 # the adaptation sees e + n(t); the plant does not
                x, al, lb = z[0], z[1:8], z[8:]
                b = np.exp(lb); a, da = L.amp(al)
                phi = np.exp(-(x - cen) ** 2 / (2 * wid ** 2))
                em = x - xr(t) + noise(t)[0]
                return np.r_[F0(t, z)[0], L.G * b * da * em * phi, L.gamma * a * em * phi]
        z0 = np.r_[xr(0.0) + 0.3, al0, np.log(b0)]
        ts = np.linspace(0, T_D, 8001)
        sol = solve_ivp(F, (0, T_D), z0, method='DOP853', rtol=1e-12, atol=1e-14, t_eval=ts)
        Z = sol.y
        V = np.array([L.real_lyapunov(Z[:, j], t, xr, al_star, q) for j, t in enumerate(ts)])
        e = Z[0] - xr(ts)
        dV = np.empty_like(ts)
        for j, t in enumerate(ts):
            zd = F0(t, Z[:, j]); hh = 1e-6
            dV[j] = (L.real_lyapunov(Z[:, j] + hh * zd, t + hh, xr, al_star, q) - L.real_lyapunov(Z[:, j] - hh * zd, t - hh, xr, al_star, q)) / (2 * hh)
        etas = np.array([L.eta(Z[1:8, j], q) for j in range(len(ts))])
        pair = np.abs((etas[:, 2] - e0[2]) - (etas[:, 6] - e0[6])).max() / max(1.0, np.abs(etas[:, 2] - e0[2]).max())
        Qs = np.array([L.Q(Z[1:8, j], np.exp(Z[8:, j])) for j in range(len(ts))])
        wend = np.exp(Z[8:, -1]) * L.amp(Z[1:8, -1])[0] / L.m
        out[lab] = dict(V0=float(V[0]), Vend=float(V[-1]), monotone=bool(np.all(np.diff(V) <= 1e-12 * max(1.0, V[0]))),
                        dVerr=float(np.abs(dV + KE_D * e ** 2).max() / np.abs(KE_D * e ** 2).max()) if noise is None else np.nan,
                        pair=float(pair), dQ=float(np.abs(Qs - q).max()), e_end=float(np.abs(e[-400:]).max()),
                        d6=float(np.abs(np.r_[wend[:2], wend[3:6]] - np.r_[wtrue6[:2], wtrue6[3:]]).max()),
                        dsum=float(abs(wend[2] + wend[6] - wtrue6[2])), dsplit=float(abs(wend[2] - wi)))
    return dict(out=out, wi=float(wi), wj=float(wj), wsum=float(wtrue6[2]))


def part_E():
    """The conserved law (Σβ = m) against the free law from α = 0, β = m/d."""
    warnings.filterwarnings('ignore')
    L = law()
    rows = []
    for i in range(6):
        Phi, ws = instance(2, i); y = Phi @ ws
        al0 = np.zeros(D_I); b0 = np.full(D_I, M_I / D_I); z0 = np.r_[al0, np.log(b0)]
        q0 = L.Q(al0, b0)
        wf = L.bregman_projection(Phi, y, q0, L.eta(al0, q0))[0]
        z = flow_end(L, Phi, y, z0, 1e6, conserved=True)
        wc = np.exp(z[D_I:]) * L.amp(z[:D_I])[0] / L.m
        rows.append(dict(res=float(np.linalg.norm(Phi @ wc - y)), d=float(np.abs(wc - wf).max()), dQ=float(np.abs(L.Q(z[:D_I], np.exp(z[D_I:])) - q0).max()),
                         l1c=float(np.abs(wc).sum()), l1f=float(np.abs(wf).sum()), mass=float(np.exp(z[D_I:]).sum())))
    return rows


def report(A=None, B=None, C=None, D=None, Ee=None):
    fmt = E.fmt
    L = law()
    L_ = ["# Proposition D‴.45: the free mass–amplitude law is a mirror flow (`exp_mirror_free.py`)", "",
          f"Free masses (λ ≡ 0): dα_i = Gβ_ia′_ig_i, dβ_i = γβ_ia_ig_i, a = a_max tanh(α/a_max), w = β⊙a/m; m = {M_I:g}, G = {G_I:g}, "
          f"γ = {GAM_I:g}, a_max = {AMAX:g} (c = γa_max²/(4G) = {L.c:g}) unless said otherwise; Q_i = β_i − c cosh(2α_i/a_max); on a "
          "dictionary the mirror flow dη = Φᵀ(r + ξ)dt/N of Ψ_q, ψ′ = η = ∫dα/(Gβa′) (`mirror_free`).", ""]
    if A is not None:
        L_ += ["## A. Identities", "",
               f"η's closed form against quadrature (q = 100c, 3, 0.5, c, c/2, 0, −c/2 and the branch q = −2c, 17 and 5 points): largest "
               f"relative difference {fmt(A['eta_err'])}; (ψ*)′ = w against finite differences {fmt(A['conj_err'])}; dw/dη = mD "
               f"{fmt(A['metric_err'])}.", "",
               "Q along the loop (5 × 12, 2-sparse truths, β(0) = 0.3, α(0) ~ N(0, 0.05²)), largest |Q(t) − Q(0)| over 60 times; the "
               "conserved law's shift dQ_i/dt = −γβ_iλ at the end state (relative):", "",
               "| instance | noiseless, t ≤ 10⁴ | a smooth noise (30 sinusoids, frequencies 10⁻³–3, mixed over the nodes, amplitude 0.3), t ≤ 3·10⁴ "
               "| the conserved law's shift |", "|---|---|---|---|"]
        for j, r in enumerate(A['cons']):
            L_.append(f"| {j} | {fmt(r['noiseless'])} | {fmt(r['smooth noise'])} | {fmt(r['shift'])} |")
        L_ += ["", "Q's drift along Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's update; at the lazy projection) "
               "against the Euler step's closed form (2|φ_i|²/(mN))γβ_i(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²)) per unit s: largest relative "
               "difference " + ', '.join(fmt(x) for x in A['red']) + ".", ""]
    if B is not None:
        L_ += [f"## B. The implicit bias: the start's masses set the regime ({N_I} × {D_I}, {N_INST_B} instances per k)", "",
               "From α = 0 and equal masses β = s·m/d: the Bregman projection (Newton on the dual) — ‖w‖₁ over basis pursuit's, the atoms "
               "used (|w| > 10⁻³), the fraction recovering the truth (10⁻³·max|w*|), the median distances to basis pursuit and (relative) to "
               f"the minimum-norm fit Φ⁺y, the largest |a|/a_max and β; for s ≥ {FLOW_MIN_B:g} the noiseless loop's end at t = {T_B:.0e} "
               "against the projection (its largest residual and distance).", "",
               "| k | s | ‖w‖₁/‖w_bp‖₁ | used | the truth | median ‖w − w_bp‖∞ | median ‖w − Φ⁺y‖∞/‖Φ⁺y‖∞ | max |a|/a_max, max β | the loop's end: "
               "residual, distance to the projection |", "|---|---|---|---|---|---|---|---|---|"]
        for k in KS_B:
            rs = [r for r in B if r['k'] == k]
            L_.append(f"| {k} | basis pursuit | 1 | | {np.mean([r['bp_ok'] for r in rs]):.2f} | | | | |")
            for s in SCALES_B:
                q = [r[s] for r in rs]
                fl = (f"{fmt(max(x['flow_res'] for x in q))}, {fmt(max(x['flow_d'] for x in q))}" if 'flow_d' in q[0] else "—")
                L_.append(f"| {k} | {s:g} | {np.mean([x['l1'] / r['l1_bp'] for x, r in zip(q, rs)]):.4f} | {np.mean([x['used'] for x in q]):.1f} | "
                          f"{np.mean([x['ok'] for x in q]):.2f} | {fmt(np.median([x['d_bp'] for x in q]))} | {fmt(np.median([x['d_l2'] for x in q]))} | "
                          f"{np.max([x['amax'] for x in q]):.3f}, {np.max([x['bmax'] for x in q]):.3g} | {fl} |")
        L_ += ["", "Small random amplitudes (α(0) ~ N(0, 0.05²)) at tiny masses: the particles with q < −c sit on a branch (the sign of w "
               "frozen); the projection's atoms used, those on a branch, the fraction of instances whose branches kept their start's signs, "
               "the distance to basis pursuit:", "",
               "| k | s | branches (mean) | used | used on a branch | signs kept | median ‖w − w_bp‖∞ |",
               "|---|---|---|---|---|---|---|"]
        for k in KS_B:
            rs = [r for r in B if r['k'] == k]
            for s in (1e-2, 1e-4):
                q = [r[('branch', s)] for r in rs]
                L_.append(f"| {k} | {s:g} | {np.mean([x['n_branch'] for x in q]):.1f} | {np.mean([x['used'] for x in q]):.1f} | "
                          f"{np.mean([x['used_br'] for x in q]):.2f} | {np.mean([x['signs_kept'] for x in q]):.2f} | {fmt(np.median([x['d_bp'] for x in q]))} |")
        L_ += ["", "Unequal masses β⁰_i = s·r_i (r_i log-uniform in [1/4, 4]), α = 0: the projection against weighted basis pursuit "
               "(min Σ|w_i|/√r_i) as s → 0 — and against the unweighted one — and against the weighted minimum norm (min Σw_i²/r_i²) as "
               "s → ∞ (relative), medians over the instances:", "", "| k | " + " | ".join(f"s = {e:g}: weighted, unweighted ℓ¹" for e in EPS_RICH)
               + " | " + " | ".join(f"s = {b:g}: weighted ℓ²" for b in SCALES_LAZY) + " |", "|---|" + "---|" * (len(EPS_RICH) + len(SCALES_LAZY))]
        for k in KS_B:
            rs = [r for r in B if r['k'] == k]
            L_.append(f"| {k} | " + " | ".join(f"{fmt(np.median([r[('rich', e)]['d'] for r in rs]))}, {fmt(np.median([r[('rich', e)]['d_unw'] for r in rs]))}"
                                               for e in EPS_RICH) + " | "
                      + " | ".join(fmt(np.median([r[('lazy', b)]['d'] for r in rs])) for b in SCALES_LAZY) + " |")
        L_.append("")
    if C is not None:
        L_ += ["## C. Noise", "", "(1) Continuous time: the loop from α = 0, β = 0.1·m/d driven by a smooth noise (30 sinusoids, frequencies "
               "10⁻³–3, mixed over the nodes, amplitude 0.3) — |Q − Q(0)|, the residual, and the distance of the state's projection to the "
               "noiseless loop's limit:", "", "| instance | t | max |Q − Q(0)| | |Φw − y| | the projection's distance |", "|---|---|---|---|---|"]
        for j, rows in enumerate(C['c1']):
            for r in rows:
                L_.append(f"| {j} | {r['t']:g} | {fmt(r['dQ'])} | {r['res']:.3f} | {fmt(r['d_int'])} |")
        L_ += ["", "(2) The Euler step on (α, log β) with label noise (σ = 0.3, 4·10⁵ draws at random states): the mean change of Q against "
               "the closed form (h²E[g²]γβ/2)(γa² − Gβa′(1 + a²/a_max²)) — largest relative difference (largest deviation in standard errors): "
               + '; '.join(f"h = {r['h']:g}: {fmt(r['rel'])} ({r['z']:.1f})" for r in C['c2']) + ".", ""]
        if 'c3' in C:
            L_ += ["(3) The selection in the mirror coordinates: Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's "
                   "update) from the lazy projection (s = 1) — ‖w‖₁, the distance to basis pursuit, the atoms used, Q's mean and range:", "",
                   "| instance (basis pursuit = truth) | s | ‖w‖₁ | ‖w − w_bp‖∞ | used | Q: mean, min, max |", "|---|---|---|---|---|---|"]
            for r in C['c3']:
                for s, x in zip(S_C3, r['rows']):
                    L_.append(f"| {r['i']} ({'yes' if r['bp_ok'] else 'no'}) | {s:.0e} | {x['l1']:.4f} | {fmt(x['d_bp'])} | {x['used']} | "
                              f"{x['qmean']:.4g}, {x['qmin']:.4g}, {x['qmax']:.4g} |")
            L_.append("")
        if 'c4' in C:
            L_ += [f"(4) The sampled-data loop at the step h = {H_C4:g} with label noise drawn afresh (σ = {SIG_C4:g}, {P_C4} paths, "
                   f"{N_C4} steps) from the lazy projection (α = 0, β = m/d): the paradigm's Euler update on (α, log β) against the mirror "
                   "update η⁺ = η + hg on each particle's curve — Q's change (mean over the particles and paths, largest), the projection of "
                   "the end state: its distance to the start's (mean over the paths), to basis pursuit's, its ‖w‖₁; the end's residual:", "",
                   "| instance | the start: ‖w‖₁, ‖w − w_bp‖∞ | the update | Q − Q(0): mean, max | ‖P − P(0)‖∞ | ‖P − w_bp‖∞ | ‖P‖₁ | |Φw − y| |",
                   "|---|---|---|---|---|---|---|---|"]
            for r in C['c4']:
                for upd in ('Euler', 'mirror'):
                    x = r[upd]
                    L_.append(f"| {r['i']} | {r['l1_0']:.4f}, {fmt(r['d_bp0'])} | {upd} | {x['dQ_mean']:+.3e}, {fmt(x['dQ_max'])} | "
                              f"{fmt(x['d_proj'])} | {fmt(x['d_bp'])} | {x['l1']:.4f} | {x['res']:.3f} |")
            L_.append("")
    if D is not None:
        L_ += [f"## D. The real loop (k_e = {KE_D:g}, γ = {GAM_D:g}, m = {M_D:g}, six Gaussian atoms of width 0.4 on [−1.2, 1.2] and a second "
               f"particle on the third; x_m = sin t + 0.4 sin 2.7t; {T_D:g} time units)", "",
               f"The pair's Bregman split of w*₃ = {D['wsum']:.5f}: {D['wi']:.5f} + {D['wj']:.5f}.", "",
               "| the adaptation's error | V(0) | V(T) | V monotone | max |V̇ + k_ee²|/max k_ee² | the pair's mirror coordinates "
               "(max |Δη_i − Δη_j|, relative) | max |Q − Q(0)| | max |e| (last 10%) | the distinct atoms' |w − w*| | the pair's sum, split |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for lab, r in D['out'].items():
            L_.append(f"| {lab} | {fmt(r['V0'])} | {fmt(r['Vend'])} | {'yes' if r['monotone'] else 'no'} | {fmt(r['dVerr'])} | {fmt(r['pair'])} | "
                      f"{fmt(r['dQ'])} | {fmt(r['e_end'])} | {fmt(r['d6'])} | {fmt(r['dsum'])}, {fmt(r['dsplit'])} |")
        L_.append("")
    if Ee is not None:
        L_ += ["## E. The conserved law (Σβ = m) from the same start (α = 0, β = m/d; t = 10⁶)", "",
               "| instance | its residual | max |Q − Q(0)| | ‖w‖₁: conserved, free (the projection) | the ends' distance | Σβ |", "|---|---|---|---|---|---|"]
        for j, r in enumerate(Ee):
            L_.append(f"| {j} | {fmt(r['res'])} | {fmt(r['dQ'])} | {r['l1c']:.4f}, {r['l1f']:.4f} | {fmt(r['d'])} | {r['mass']:.6f} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_mirror_free.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_free_raw.npz')
    md = os.path.join(OUT, 'mirror_free.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDE']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDE')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(k, i) for k in KS_B for i in range(N_INST_B)], chunksize=1)
        if 'C' in parts:
            h['C3'] = pool.map_async(c3_job, list(range(N_INST_C3)), chunksize=1)
        if 'A' in parts:
            res['A'] = part_A(); E.save(raw, **res); log('A')
        if 'D' in parts:
            res['D'] = part_D(); E.save(raw, **res); log('D')
        if 'E' in parts:
            res['E'] = part_E(); E.save(raw, **res); log('E')
        if 'C' in parts:
            C = part_C12(); C['c4'] = part_C4(); C['c3'] = h['C3'].get(); res['C'] = C; E.save(raw, **res); log('C')
        if 'B' in parts:
            res['B'] = h['B'].get(); E.save(raw, **res); log('B')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE']))


if __name__ == '__main__':
    main()
