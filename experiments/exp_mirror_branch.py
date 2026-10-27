"""
Proposition D‴.47: the free law's branches — an exponentiated gradient with frozen signs; the limit for every start is the Bregman
projection onto the sign-constrained least-squares fits, the wrong-signed particles going extinct; positive kernel dictionaries
====================================================================================================================================

The free mass–amplitude law of Proposition D‴.45 with particles on branches (Q_i < −c: a start with β_i(0) < 2c sinh²(α_i(0)/a_max)),
represented without cancellation near their branch points (`mirror_free.ParticleSet`).

A. Identities: the slope dη/dlog|w| = a/(γa² + Gβa′²) at random states (finite differences); the temperature 1/(γa₀) at branch points;
   B(w̄, w) → β/(γm) as w̄ → 0 on a branch; the fraction of small random starts on branches; a ± pair's coefficient against
   2u sinh(γa₀Δ) and its potential against the hyperbolic entropy; the critical curve's power law.
B. The limit on random Gaussian designs (6 × 12) from small random starts (all on branches, random or positive signs; and a mixed start):
   the flow to t = 3·10⁷ against the Bregman projection onto S (computed independently, `sign_ls_set` and `projection`), J − min_K J,
   the extinct masses and their decay rates against γa₀|g_i(∞)|, B(w̄, w(t)) nonincreasing.
C. The rich limits on one design (8 × 20): branch starts at masses ε → the sign-constrained basis pursuit weighted by 1/|a(0)|; starts at
   α = 0 → basis pursuit; a mixed start (half at α = 0).
D. A pair per atom: the rich limit of ± pairs → basis pursuit weighted by 1/a₀.
E. Positive kernel dictionaries (Gaussian atoms on a grid, few samples, a positive 3-sparse target with two adjacent atoms): the
   nonnegative certificate's margin; basis pursuit with equal and with column-norm weights; the loop from positive branches at four scales
   and from α = 0 at two, to t = 10⁹; the other atoms' mass against the certificate's bound.
F. The real loop on positive branches: V̇ = −k_ee² with w* vanishing on unused atoms (random states); a run (twelve atoms, a positive
   target on three): V, the error, the unused atoms' masses by quarters.

Output: results/chaos/mirror_branch.md and mirror_branch_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D,E,F`).
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
from scipy.integrate import solve_ivp, quad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mirror_free as MF
import exp_cohesion_leak as E

OUT = os.environ.get('MB_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
LAW_B = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)                 # c = 0.025
LAW_E = dict(m=12.0, G=0.25, gamma=1.0, a_max=1.0)                # c = 1
SEEDS_B, NB, DB, T_B = range(8), 6, 12, 3e7
NC, DC, EPS_BR, EPS_FREE = 8, 20, (1e-4, 1e-8, 1e-16, 1e-32, 1e-64, 1e-128), (1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
U_PAIR = (1e-4, 1e-8, 1e-16, 1e-32, 1e-64, 1e-128)
CASES_E = ((21, 0.15, 7), (41, 0.15, 7), (41, 0.1, 11))
STARTS_E = (('branch', 1.0, 1e-6), ('branch', 1.0, 1e-2), ('branch', 1.0, 1.0), ('branch', 2.0, 10.0), ('alpha=0', 0.0, 1e-6),
            ('alpha=0', 0.0, 1e2))
T_E = 1e9
LAW_F = dict(m=12.0, G=1.0, gamma=1.0, a_max=1.0)                # c = 0.25
KE_F, T_F = 5.0, 4000.0
XR = (lambda t: np.sin(t) + 0.4 * np.sin(2.7 * t))
XRD = (lambda t: np.cos(t) + 1.08 * np.cos(2.7 * t))


def part_A():
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(0)
    # A1: the slope at random states, free and on branches (finite differences along each particle's state)
    al_i = rng.uniform(-1.5, 1.5, 40); b_i = 10 ** rng.uniform(-6, 1, 40)
    P = MF.ParticleSet(L, al_i, b_i)
    x = P.x0.copy(); h = 1e-6
    err = 0.0
    for j in range(P.d):
        xp, xm = x.copy(), x.copy(); xp[j] += h; xm[j] -= h
        num = (P.eta(xp)[j] - P.eta(xm)[j]) / (np.log(abs(P.w(xp)[j])) - np.log(abs(P.w(xm)[j])))
        err = max(err, abs(num / P.temperature(x)[j] - 1))
    # A2: the temperature at branch points
    tmp = []
    for a0 in (0.05, 0.3, 1.0, 2.0):
        Pb = MF.ParticleSet(L, np.array([a0 + 1e-30]), np.array([2 * L.c * np.sinh(a0) * np.sinh(1e-30) * 2]))
        xb = np.array([np.log(1e-40)])
        Pb.al0[:] = a0
        tmp.append(dict(a0=a0, T=float(Pb.temperature(xb)[0]), ref=1 / (L.gamma * L.A * np.tanh(a0 / L.A))))
    # A3: B(w̄, w) → β/(γm) as w̄ → 0 on a branch
    Pc = MF.ParticleSet(L, np.array([0.4]), np.array([1e-3]))
    xw = np.array([np.log(0.05)])
    b0 = []
    for dbar in (1e-6, 1e-12, 1e-30, 1e-100):
        xbar = np.array([np.log(dbar)])
        Bv = Pc.bdiv(Pc.w(xbar), xbar, xw)
        b0.append(dict(dbar=dbar, wbar=float(Pc.w(xbar)[0]), B=Bv, ref=float(Pc.beta(xw)[0] / (L.gamma * L.m))))
    # A4: the fraction of small random starts on branches
    frac = []
    for b in (1e-3, 1e-5, 1e-7, 1e-9):
        al = rng.uniform(-0.3, 0.3, 20000)
        frac.append(dict(b=b, frac=float(np.mean(b < 2 * L.c * np.sinh(al / L.A) ** 2))))
    # A5: a ± pair: W(Δ) against 2u sinh(γa₀Δ), the potential against the hyperbolic entropy (relative errors by |W|/u)
    al0, bb = 0.3, 1e-7
    Pp = MF.ParticleSet(L, np.array([al0, -al0]), np.array([bb, bb]))
    e0 = Pp.eta(Pp.x0); u = float(Pp.w(Pp.x0)[0]); a0 = L.A * np.tanh(Pp.al0[0] / L.A)
    pair = []
    for D in (1.0, 10.0, 30.0, 60.0, 90.0, 120.0):
        xs = Pp.x_of_eta(e0 + D)
        W = float(Pp.w(xs).sum())
        ref = 2 * u * np.sinh(L.gamma * a0 * D)
        pair.append(dict(D=D, ratio=W / u, rel=abs(W / ref - 1)))
    # A6: the critical curve Q = −c: η·w^{1/3} → −(a_max²/(2Gc))(2c/(ma_max²))^{1/3} as α → 0⁺ (η from α to 1 by quadrature, shifted)
    A, G, c, m = L.A, L.G, L.c, L.m
    beta_c = lambda a_: 2 * c * np.sinh(a_ / A) ** 2
    integrand = lambda a_: 1 / (G * beta_c(a_) * (1 - np.tanh(a_ / A) ** 2))
    crit = []
    for a_ in (1e-2, 1e-3, 1e-4):
        eta_ = -quad(integrand, a_, 1.0, limit=200, epsabs=0, epsrel=1e-12)[0]
        w_ = beta_c(a_) * A * np.tanh(a_ / A) / m
        crit.append(dict(alpha=a_, val=eta_ * w_ ** (1 / 3), ref=-(A ** 2 / (2 * G * c)) * (2 * c / (m * A ** 2)) ** (1 / 3)))
    return dict(slope=err, temp=tmp, b0=b0, frac=frac, pair=pair, u=u, a0=a0, crit=crit)


def gauss_design(seed, N, d):
    return np.random.default_rng(seed).standard_normal((N, d))


def flow_B(P, Phi, y, T, n=60):
    ts = np.concatenate([[0], np.geomspace(1e-1, T, n)])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sol = solve_ivp(P.rhs(Phi, y), (0, T), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ts)
    return sol


def b_case(seed, kind):
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(100 + seed)
    Phi = gauss_design(seed, NB, DB)
    ws = np.zeros(DB); ws[[0, 3, 5, 8]] = [0.9, -0.7, 0.6, -0.5]; y = Phi @ ws
    amp = rng.uniform(0.2, 0.4, DB)
    if kind == 'positive':
        al_i, b_i = amp, np.full(DB, 1e-3)
    elif kind == 'random signs':
        al_i, b_i = amp * rng.choice([-1, 1], DB), np.full(DB, 1e-3)
    else:                                                                    # mixed: half the particles at α = 0 (free signs)
        al_i = amp * rng.choice([-1, 1], DB); al_i[::2] = 0.0; b_i = np.full(DB, 1e-3)
    P = MF.ParticleSet(L, al_i, b_i)
    eta0 = P.eta(P.x0)
    yh, Z, _ = MF.sign_ls_set(Phi, y, P.s)
    keep = np.setdiff1d(np.arange(DB), Z)
    wbar, it, res = MF.projection(P, Phi, yh, eta0, keep)
    sol = flow_B(P, Phi, y, T_B)
    W = np.array([P.w(sol.y[:, k]) for k in range(sol.y.shape[1])])
    Bt = np.array([P.beta(sol.y[:, k]) for k in range(sol.y.shape[1])])
    # B(w̄, w(t)) along the flow (w̄ = the projection; its zeros on Z read as β/(γm)); w̄'s own state off Z by root finding
    xb = P.state_of_w(wbar, keep)
    Bv = np.array([P.bdiv(wbar, xb, sol.y[:, k]) for k in range(sol.y.shape[1])])
    incr = float(max(0.0, np.max(np.diff(Bv)) / max(Bv[0], 1e-300)))
    J = 0.5 * NB ** -1 * L.m * np.sum((W @ Phi.T - y) ** 2, axis=1)
    Jstar = 0.5 * L.m / NB * np.sum((yh - y) ** 2)
    g_inf = -(Phi.T @ (Phi @ wbar - y)) / NB
    a0 = L.A * np.tanh(P.al0 / L.A)
    rates = []
    for i in Z:
        sel = (Bt[:, i] > 1e-250) & (Bt[:, i] < 1e-12)
        if sel.sum() >= 3:
            sl = np.polyfit(sol.t[sel], np.log(Bt[sel, i]), 1)[0]
            rates.append(dict(i=int(i), fit=float(-sl), ref=float(L.gamma * a0[i] * abs(g_inf[i])), sg=float(P.s[i] * g_inf[i])))
    return dict(seed=seed, kind=kind, nbr=int(P.br.sum()), yres=float(np.linalg.norm(y - yh)), Z=[int(i) for i in Z], it=it, res=res,
                status=int(sol.status), dist=[float(np.abs(W[k] - wbar).max()) for k in range(W.shape[0])], t=sol.t.tolist(),
                dJ=[float(J[k] - Jstar) for k in range(W.shape[0])], bZ=[float(Bt[k, Z].max()) if len(Z) else 0.0 for k in range(W.shape[0])],
                Bincr=incr, rates=rates)


def part_B():
    return [b_case(s, kind) for kind in ('positive', 'random signs', 'mixed') for s in SEEDS_B]


def part_C():
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(7)
    Phi = rng.standard_normal((NC, DC))
    ws = np.zeros(DC); ws[[1, 5, 9]] = [0.8, -0.6, 0.5]; y = Phi @ ws
    sg = rng.choice([-1.0, 1.0], DC); sg[[1, 5, 9]] = np.sign(ws[[1, 5, 9]])
    amp0 = rng.uniform(0.2, 0.6, DC)
    a_init = L.A * np.tanh(amp0 / L.A)
    w_scbp = MF.weighted_bp(Phi, y, 1 / a_init, signs=sg)
    w_bp = MF.weighted_bp(Phi, y, np.ones(DC))
    br = []
    for eps in EPS_BR:
        P = MF.ParticleSet(L, sg * amp0, np.full(DC, eps))
        wb, it, res = MF.projection(P, Phi, y, P.eta(P.x0))
        br.append(dict(eps=eps, all_branch=bool(P.br.all()), it=it, res=res, d=float(np.abs(wb - w_scbp).max())))
    fr = []
    for eps in EPS_FREE:
        P = MF.ParticleSet(L, np.zeros(DC), np.full(DC, eps))
        wb, it, res = MF.projection(P, Phi, y, P.eta(P.x0))
        fr.append(dict(eps=eps, it=it, res=res, d=float(np.abs(wb - w_bp).max())))
    # a mixed start: the even atoms at α = 0 (free), the odd ones on branches with the truth's signs where it has them
    mx = []
    odd = np.arange(DC) % 2 == 1
    w_scbp_odd = np.zeros(DC)
    sub = MF.weighted_bp(Phi[:, odd], y, 1 / a_init[odd], signs=sg[odd])
    w_scbp_odd[odd] = sub
    for eps in (1e-4, 1e-8, 1e-12):
        al_i = np.where(odd, sg * amp0, 0.0)
        P = MF.ParticleSet(L, al_i, np.full(DC, eps))
        wb, it, res = MF.projection(P, Phi, y, P.eta(P.x0))
        mx.append(dict(eps=eps, it=it, res=res, free=float(np.abs(wb[~odd]).sum()), d=float(np.abs(wb - w_scbp_odd).max())))
    return dict(scbp_supp=[int(i) for i in np.flatnonzero(np.abs(w_scbp) > 1e-9)], scbp_is_truth=bool(np.abs(w_scbp - ws).max() < 1e-8),
                bp_supp=[int(i) for i in np.flatnonzero(np.abs(w_bp) > 1e-9)], br=br, fr=fr, mx=mx,
                odd_feasible=bool(np.all(np.isfinite(sub))))


def part_D():
    L = MF.FreeLaw(**LAW_B)
    rng = np.random.default_rng(7)
    Phi = rng.standard_normal((NC, DC))
    ws = np.zeros(DC); ws[[1, 5, 9]] = [0.8, -0.6, 0.5]; y = Phi @ ws
    al = 0.4
    Phi2 = np.hstack([Phi, Phi])
    rows = []
    for u in U_PAIR:
        b = u * L.m / (L.A * np.tanh(al / L.A))                               # w(0) = u for each member
        P = MF.ParticleSet(L, np.r_[np.full(DC, al), np.full(DC, -al)], np.full(2 * DC, b))
        assert P.br.all(), u
        wb, it, res = MF.projection(P, Phi2, y, P.eta(P.x0))
        W = wb[:DC] + wb[DC:]
        a0 = L.A * np.tanh(P.al0[0] / L.A)
        wref = MF.weighted_bp(Phi, y, np.full(DC, 1 / a0))
        rows.append(dict(u=u, it=it, res=res, d=float(np.abs(W - wref).max()), both=int(np.sum((np.abs(wb[:DC]) > 1e-3 * np.abs(W).max())
                                                                                             & (np.abs(wb[DC:]) > 1e-3 * np.abs(W).max())))))
    return dict(rows=rows, bp_supp=[int(i) for i in np.flatnonzero(np.abs(MF.weighted_bp(Phi, y, np.ones(DC))) > 1e-9)])


def e_case(d, wid, N):
    L = MF.FreeLaw(**LAW_E)
    cen = np.linspace(-1, 1, d); xs = np.linspace(-1.1, 1.1, N)
    Phi = np.exp(-(xs[:, None] - cen[None, :]) ** 2 / (2 * wid ** 2))
    i0 = d // 4; supp = np.array([i0, i0 + 1, 3 * d // 4])
    ws = np.zeros(d); ws[supp] = [1.0, 0.7, 0.5]; y = Phi @ ws
    mu, v = MF.nonneg_certificate(Phi, supp)
    w_bp = MF.weighted_bp(Phi, y, np.ones(d)); w_bpn = MF.weighted_bp(Phi, y, np.linalg.norm(Phi, axis=0))
    rows = []
    for lab, al_v, b_v in STARTS_E:
        P = MF.ParticleSet(L, np.full(d, al_v), np.full(d, b_v))
        ts = np.concatenate([[0], np.geomspace(1.0, T_E, 36)])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sol = solve_ivp(P.rhs(Phi, y), (0, T_E), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ts)
        Wt = np.array([P.w(sol.y[:, k]) for k in range(sol.y.shape[1])])
        off = np.setdiff1d(np.arange(d), supp)
        row = dict(lab=lab, al=al_v, b=b_v, branch=bool(P.br.all()), t=sol.t.tolist(), err=np.abs(Wt - ws).max(axis=1).tolist(),
                   offmass=np.abs(Wt[:, off]).sum(axis=1).tolist(), res=float(np.linalg.norm(Phi @ Wt[-1] - y)),
                   supp=[int(i) for i in np.flatnonzero(np.abs(Wt[-1]) > 1e-4 * np.abs(Wt[-1]).max())])
        if P.br.all() and v is not None:
            # the certificate's bound (|v|/μ)(2N·B(w*, w(0))/t)^{1/2}: B(w*, w(0)) with w*'s zeros read as β/(γm)
            xs_ = P.state_of_w(ws, supp)
            B0 = P.bdiv(ws, xs_, P.x0)
            row['bound'] = [float(np.linalg.norm(v) / mu * np.sqrt(2 * N * B0 / t)) if t > 0 else float('inf') for t in sol.t]
        rows.append(row)
    return dict(d=d, wid=wid, N=N, supp=supp.tolist(), mu=float(mu), bp=float(np.abs(w_bp - ws).max()), bpn=float(np.abs(w_bpn - ws).max()),
                rows=rows)


def part_E():
    return [e_case(*c) for c in CASES_E]


def part_F():
    L = MF.FreeLaw(**LAW_F)
    rng = np.random.default_rng(5)
    d = 12
    cen, wid = np.linspace(-1.4, 1.4, d), np.full(d, 0.25)
    wstar = np.zeros(d); wstar[[2, 6, 9]] = [0.12, 0.08, 0.1]
    P = MF.ParticleSet(L, np.full(d, 0.5), np.full(d, 1e-3))
    assert P.br.all()
    f = P.real_rhs(cen, wid, wstar, XR, XRD, KE_F)
    # the identity V̇ = −k_ee² at random states (w* on the curves of the used atoms, β/(γm) on the others)
    used = np.flatnonzero(wstar)
    xs_ = P.state_of_w(wstar, used)
    V = lambda t, z: 0.5 * (z[0] - XR(t)) ** 2 + P.bdiv(wstar, xs_, z[1:])
    worst = 0.0
    for k in range(12):
        t = rng.uniform(0, 50); z = np.r_[XR(t) + rng.normal(0, 0.3), P.x0 + rng.normal(0, 0.5, d)]
        dz = f(t, z); h = 1e-6
        num = (V(t + h, z + h * dz) - V(t - h, z - h * dz)) / (2 * h)
        e = z[0] - XR(t)
        worst = max(worst, abs(num + KE_F * e ** 2) / max(abs(KE_F * e ** 2), 1e-12))
    z0 = np.r_[XR(0.0), P.x0]
    ts = np.linspace(0, T_F, 8001)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sol = solve_ivp(f, (0, T_F), z0, method='LSODA', rtol=1e-9, atol=1e-12, t_eval=ts)
    Vt = np.array([V(t, sol.y[:, k]) for k, t in enumerate(sol.t)])
    e = sol.y[0] - XR(sol.t)
    unused = np.setdiff1d(np.arange(d), used)
    q = np.array_split(np.arange(sol.t.size), 4)
    Bu = np.array([P.beta(sol.y[1:, k])[unused].max() for k in range(sol.t.size)])
    Wend = P.w(sol.y[1:, -1])
    return dict(dV=worst, V=[float(Vt[i].max()) for i in q], Vinc=float(max(0.0, np.max(np.diff(Vt)))), erms=[float(np.sqrt(np.mean(e[i] ** 2))) for i in q],
                bunused=[float(Bu[i].max()) for i in q], err_end=float(np.abs(Wend - wstar).max()), V0=float(Vt[0]))


def report(A=None, B=None, C=None, D=None, Ee=None, F=None):
    fmt = E.fmt
    L_ = ["# Proposition D‴.47: the free law's branches (`exp_mirror_branch.py`)", "",
          "The free mass–amplitude law (Proposition D‴.45) with particles on branches (Q < −c), represented without cancellation near their "
          "branch points (`mirror_free.ParticleSet`).", ""]
    if A is not None:
        L_ += ["## A. Identities", "",
               f"dη/dlog|w| against a/(γa² + Gβa′²) at 40 random states (free and on branches; finite differences): largest relative "
               f"difference {fmt(A['slope'])}.", "", "The temperature at the branch point against 1/(γa₀):", "", "| α₀ | a/(γa² + Gβa′²) | 1/(γa₀) |",
               "|---|---|---|"]
        L_ += [f"| {r['a0']:g} | {r['T']:.10g} | {r['ref']:.10g} |" for r in A['temp']]
        L_ += ["", "B(w̄, w) on a branch as w̄ → 0 (δ̄ the distance to the branch point) against β/(γm):", "",
               "| δ̄ | w̄ | B(w̄, w) | β/(γm) |", "|---|---|---|---|"]
        L_ += [f"| {r['dbar']:g} | {fmt(r['wbar'])} | {r['B']:.12g} | {r['ref']:.12g} |" for r in A['b0']]
        L_ += ["", "The fraction of starts on branches, α(0) uniform on [−0.3, 0.3]: "
               + ", ".join(f"{r['frac']:.4f} at β(0) = {r['b']:g}" for r in A['frac']) + ".", "",
               f"A ± pair (α(0) = ±0.3, u = w₊(0) = {fmt(A['u'])}, a₀ = {A['a0']:.6g}): the atom's coefficient W(Δ) against 2u sinh(γa₀Δ):", "",
               "| Δ | W/u | relative difference |", "|---|---|---|"]
        L_ += [f"| {r['D']:g} | {fmt(r['ratio'])} | {fmt(r['rel'])} |" for r in A['pair']]
        L_ += ["", "The critical curve Q = −c: η·w^{1/3} (η by quadrature, up to its constant) against −(a_max²/(2Gc))(2c/(ma_max²))^{1/3} = "
               f"{A['crit'][0]['ref']:.6g}: " + ", ".join(f"{r['val']:.6g} at α = {r['alpha']:g}" for r in A['crit']) + ".", ""]
    if B is not None:
        L_ += [f"## B. The limit on random Gaussian designs ({NB} × {DB}, the target on four atoms with signs +, −, +, −; masses 10⁻³, amplitudes "
               "0.2–0.4)", "",
               "The flow (LSODA) to t = 3·10⁷ against the Bregman projection of its start onto S (computed independently); |y − ŷ| the "
               "sign-constrained residual, |Z| the coordinates vanishing on S; J − min_K J at the end; the largest mass on Z at the end; "
               "B(w̄, w(t))'s largest relative increase; the extinction rates (fitted where 10⁻²⁵⁰ < β < 10⁻¹²) against γa₀|g(∞)|:", "",
               "| start | seed | branches | |y − ŷ| | |Z| | projection (it, residual) | max |w(T) − w̄| | J − J* | max β on Z | B's increase | rates: largest relative error (n) |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in B:
            rr = max((abs(x['fit'] / x['ref'] - 1) for x in r['rates']), default=float('nan'))
            L_.append(f"| {r['kind']} | {r['seed']} | {r['nbr']} | {fmt(r['yres'])} | {len(r['Z'])} | {r['it']}, {fmt(r['res'])} | {fmt(r['dist'][-1])} | "
                      f"{fmt(r['dJ'][-1])} | {fmt(r['bZ'][-1])} | {fmt(r['Bincr'])} | {fmt(rr)} ({len(r['rates'])}) |")
        L_.append("")
    if C is not None:
        L_ += [f"## C. The rich limits on one design ({NC} × {DC}; the target on atoms [1, 5, 9])", "",
               f"Sign-constrained basis pursuit weighted by 1/|a(0)| (the frozen signs: the target's on its atoms, random elsewhere): support "
               f"{C['scbp_supp']} (the target: {'yes' if C['scbp_is_truth'] else 'no'}); basis pursuit: support {C['bp_supp']}.", "",
               "| start | ε | all on branches | projection (it, residual) | distance to its limit |", "|---|---|---|---|---|"]
        L_ += [f"| branches (amplitudes 0.2–0.6) | {r['eps']:g} | {'yes' if r['all_branch'] else 'no'} | {r['it']}, {fmt(r['res'])} | {fmt(r['d'])} |" for r in C['br']]
        L_ += [f"| α(0) = 0 | {r['eps']:g} | no | {r['it']}, {fmt(r['res'])} | {fmt(r['d'])} |" for r in C['fr']]
        L_ += ["", "A mixed start (the even atoms at α = 0, the odd ones on branches): the free particles' total |w| and the distance to the "
               "branch atoms' own sign-constrained weighted basis pursuit:", "", "| ε | projection (it, residual) | Σ|w| on free particles | distance |",
               "|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | {r['it']}, {fmt(r['res'])} | {fmt(r['free'])} | {fmt(r['d'])} |" for r in C['mx']]
        L_.append("")
    if D is not None:
        L_ += ["## D. A pair per atom (± 0.4, equal masses; the design of C)", "",
               "The atom's coefficient W = w₊ + w₋ of the projection against basis pursuit weighted by 1/a₀; the atoms where both members "
               f"carry more than 10⁻³ of the largest |W| (basis pursuit's support {D['bp_supp']}):", "",
               "| u = w(0) | projection (it, residual) | distance | atoms with both members active |", "|---|---|---|---|"]
        L_ += [f"| {r['u']:g} | {r['it']}, {fmt(r['res'])} | {fmt(r['d'])} | {r['both']} |" for r in D['rows']]
        L_.append("")
    if Ee is not None:
        L_ += ["## E. Positive kernel dictionaries (Gaussian atoms of one width on a grid of [−1, 1], samples evenly on [−1.1, 1.1], the target "
               "1.0, 0.7, 0.5 on two adjacent atoms and one more)", ""]
        for c in Ee:
            L_ += [f"### {c['d']} atoms, width {c['wid']:g}, {c['N']} samples (the target on {c['supp']})", "",
                   f"The nonnegative certificate's margin μ = {fmt(c['mu'])}; basis pursuit's max |w − w*|: {fmt(c['bp'])} with equal weights, "
                   f"{fmt(c['bpn'])} with the columns' norms.", "",
                   "| start | α(0) | β(0) | on branches | max |w − w*| at t = 10³, 10⁵, 10⁷, 10⁹ | the other atoms' mass at 10⁹ (the certificate's bound) | "
                   "|Φw − y| | support at 10⁹ |", "|---|---|---|---|---|---|---|---|"]
            for r in c['rows']:
                ts = np.array(r['t'])
                pick = [int(np.argmin(np.abs(np.log10(np.maximum(ts, 1e-300)) - k))) for k in (3, 5, 7, 9)]
                bd = f" ({fmt(r['bound'][-1])})" if 'bound' in r else ""
                L_.append(f"| {r['lab']} | {r['al']:g} | {r['b']:g} | {'yes' if r['branch'] else 'no'} | "
                          + ", ".join(fmt(r['err'][k]) for k in pick) + f" | {fmt(r['offmass'][-1])}{bd} | {fmt(r['res'])} | {r['supp']} |")
            L_.append("")
    if F is not None:
        L_ += [f"## F. The real loop on positive branches (k_e = {KE_F:g}, G = γ = 1, twelve atoms of width 0.25, the target on atoms [2, 6, 9], "
               f"α(0) = 0.5, β(0) = 10⁻³, {T_F:g} time units)", "",
               f"V̇ against −k_ee² at 12 random states (w* vanishing on nine atoms, their terms β/(γm)): largest relative difference {fmt(F['dV'])}.", "",
               "By quarters: max V, the rms error, the largest mass on the unused atoms; V's largest increase between samples "
               f"{fmt(F['Vinc'])} (V(0) = {fmt(F['V0'])}); max |w − w*| at the end {fmt(F['err_end'])}:", "",
               "| | quarter 1 | quarter 2 | quarter 3 | quarter 4 |", "|---|---|---|---|---|",
               "| max V | " + " | ".join(fmt(x) for x in F['V']) + " |",
               "| rms e | " + " | ".join(fmt(x) for x in F['erms']) + " |",
               "| max β, unused atoms | " + " | ".join(fmt(x) for x in F['bunused']) + " |", ""]
    L_.append("Numbers are produced by `experiments/exp_mirror_branch.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_branch_raw.npz')
    md = os.path.join(OUT, 'mirror_branch.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEF']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEF')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('C', part_C), ('D', part_D), ('F', part_F), ('B', part_B), ('E', part_E)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEF']))


if __name__ == '__main__':
    main()
