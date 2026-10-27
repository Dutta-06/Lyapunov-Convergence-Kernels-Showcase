"""
Proposition D‴.48: the rich transit — from small masses the kernels enter one at a time
====================================================================================================================================

The free mass–amplitude law (Proposition D‴.45; `mirror_free`, `rich_transit`) from small masses, in the time rescaled by the transit's
scale.

A. From α = 0 (masses εb_i, three 6 × 14 Gaussian designs, unequal b_i): the inverse scale space path of Σ|w_i|/√b_i (`iss_path`) against
   the loop in τ = t/Λ(ε), Λ = (πa_max/(4G))(2/(cε))^{1/2}, u = η/Λ, for ε = 10⁻² … 10⁻⁸: sup_τ|u − u_ISS| (with its log-log slope),
   the residual's rms difference over τ, the entry times, the end's distance to weighted basis pursuit.
B. On branches (amplitudes s_i·U(0.2, 0.6), masses ε): the one-sided path (frozen signs, weights 1/(γ|a_i(0)|)) against the loop in
   τ = t/log(1/ε), u = (η − η(0))/log(1/ε), for ε = 10⁻⁴ … 10⁻⁶⁴; a ± pair per atom (amplitudes ±0.4, masses with w(0) = u) against the
   two-sided path (weights 1/(γa₀)).
C. Moving centres (Gaussian kernels of width 0.12 on 61 samples of [−1, 1]; 21 kernels on a grid; the target three bumps of the same
   width off the grid): (1) the paradigm's law (ċ = Gβa⟨r∂_cφ⟩) against the sliding inverse scale space (`sliding_iss`) for ε = 10⁻³ …
   10⁻⁷ — entries, the formed kernels, the inactive kernels' displacement and mass; (2) the Wasserstein–Fisher–Rao metric (ċ = Ga⟨r∂_cφ⟩):
   the inactive kernels' centres against their limit (`wfr_inactive_rhs`) before the first entry — condensation —, the formed kernels;
   (3) the centres held fixed (the grid dictionary's transit); (4) random three-bump targets: the formed kernels under each metric and with
   fixed centres; (5) two bumps closer than the kernels' width, and a target of wider bumps than the kernels.

Output: results/chaos/mirror_transit.md and mirror_transit_raw.npz (`--from-raw` regenerates the md; `--only A,B,C`).
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
from src.core import rich_transit as RT
import exp_cohesion_leak as E

OUT = os.environ.get('MT_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)                       # c = 0.025
NA, DA, SEEDS_A = 6, 14, (3, 4, 5)
EPS_A = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
NT_A = 2001                                                           # τ-grid points on [0, 1.5τ*]
EPS_B, U_B = (1e-4, 1e-8, 1e-16, 1e-32, 1e-64), (1e-4, 1e-8, 1e-16, 1e-32, 1e-64)
NS_C, SIG_C, D_C = 61, 0.12, 21
CSTAR, ASTAR = (-0.52, 0.07, 0.61), (0.9, -0.6, 0.75)
EPS_C1, EPS_C2 = (1e-3, 1e-4, 1e-5, 1e-6, 1e-7), (1e-4, 1e-5, 1e-6)
TAU_PRE = 10.0                                                        # before the first entry (C2)
N_RAND, EPS_RAND, EPS_RAND2, EPS_RAND3 = 16, 1e-5, 1e-6, 1e-8                # C4: the paradigm at all three, WFR and fixed at the second
CREEP_TARGET, EPS_CREEP = 10, (1e-8, 1e-10)                                  # C6: a lost near-tie's creep
TH_ACT = 1e-2                                                         # a formed kernel: |w| > TH_ACT·max|A*|


def law():
    return MF.FreeLaw(**LAW)


def design(seed):
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((NA, DA))
    ws = np.zeros(DA); idx = rng.choice(DA, 3, replace=False); ws[idx] = rng.choice([-1, 1], 3) * rng.uniform(0.3, 0.9, 3)
    b = np.exp(rng.uniform(np.log(0.25), np.log(4.0), DA))
    return Phi, ws, b


def entries_of(W, taus, th):
    """The first τ at which |w_i| exceeds th, per coordinate (inf if never)."""
    out = np.full(W.shape[1], np.inf)
    for i in range(W.shape[1]):
        k = np.flatnonzero(np.abs(W[:, i]) > th)
        if k.size:
            out[i] = taus[k[0]]
    return out


def part_A():
    """The inverse scale space limit from α = 0."""
    warnings.filterwarnings('ignore')
    L = law()
    res = []
    for seed in SEEDS_A:
        Phi, ws, b = design(seed); y = Phi @ ws
        rho = b ** -0.5
        path = MF.iss_path(Phi, y, rho)
        wbp = MF.weighted_bp(Phi, y, rho)
        tau_end = 1.5 * path['taus'][-1]
        taus = np.linspace(0, tau_end, NT_A)
        Uref = np.array([MF.iss_eval(path, t, Phi, y)[0] for t in taus])
        Wref = np.array([MF.iss_eval(path, t, Phi, y)[1] for t in taus])
        th = 1e-3 * np.abs(wbp).max()
        ent_ref = entries_of(Wref, taus, th)                                    # the path's own crossings of the threshold
        rows = []
        for eps in EPS_A:
            t0 = time.time()
            Lam = RT.transit_scale(L, eps)
            b0 = eps * b; q = L.Q(np.zeros(DA), b0)
            sol = solve_ivp(L.rhs(Phi, y), (0, Lam * tau_end), np.r_[np.zeros(DA), np.log(b0)], method='LSODA', rtol=1e-10, atol=1e-13,
                            t_eval=Lam * taus)
            al = sol.y[:DA].T
            U = np.array([L.eta(al[k], q) for k in range(len(taus))]) / Lam
            W = np.array([L.w(al[k], q) for k in range(len(taus))])
            ent = entries_of(W, taus, th)
            both = np.isfinite(ent) & np.isfinite(ent_ref)
            rows.append(dict(eps=eps, Lam=float(Lam), du=float(np.abs(U - Uref).max()),
                             dres=float(np.sqrt(np.mean(np.sum((W @ Phi.T - Wref @ Phi.T) ** 2, axis=1)))),
                             dent=float(np.abs(ent - ent_ref)[both].max()) if both.any() else np.nan, n_only=int(np.sum(np.isfinite(ent) != np.isfinite(ent_ref))),
                             dend=float(np.abs(W[-1] - wbp).max()), secs=time.time() - t0))
        res.append(dict(seed=seed, events=[(float(t), int(i), k) for t, i, k in path['events']], tau_end=float(tau_end),
                        end_is_bp=float(np.abs(path['w_end'] - wbp).max()), rows=rows, dtau=float(taus[1] - taus[0])))
    return res


def part_B():
    """Branches (one-sided path, the scale log(1/ε)) and ± pairs (two-sided)."""
    warnings.filterwarnings('ignore')
    L = law()
    Phi, ws, b = design(SEEDS_A[0]); y = Phi @ ws
    rng = np.random.default_rng(11)
    s = rng.choice([-1.0, 1.0], DA); s[ws != 0] = np.sign(ws[ws != 0])
    amp = rng.uniform(0.2, 0.6, DA)
    a_init = L.A * np.tanh(amp / L.A)
    rho = 1 / (L.gamma * a_init)
    path = MF.iss_path(Phi, y, rho, signs=s)
    wlim = MF.weighted_bp(Phi, y, rho, signs=s)
    out = dict(events=[(float(t), int(i), k) for t, i, k in path['events']], end_is_lim=float(np.abs(path['w_end'] - wlim).max()))
    tau_end = 1.5 * path['taus'][-1]
    taus = np.linspace(0, tau_end, NT_A)
    Uref = np.array([MF.iss_eval(path, t, Phi, y)[0] for t in taus])
    Wref = np.array([MF.iss_eval(path, t, Phi, y)[1] for t in taus])
    th = 1e-3 * np.abs(wlim).max()
    ent_ref = entries_of(Wref, taus, th)
    rows = []
    for eps in EPS_B:
        t0 = time.time()
        Lg = np.log(1 / eps)
        P = MF.ParticleSet(L, s * amp, np.full(DA, eps))
        e0 = P.eta(P.x0)
        sol = solve_ivp(P.rhs(Phi, y), (0, Lg * tau_end), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=Lg * taus,
                        max_step=Lg * tau_end / NT_A)                    # a capped step: no trial state beyond the curves' range
        U = np.array([(P.eta(sol.y[:, k]) - e0) / Lg for k in range(len(taus))])
        W = np.array([P.w(sol.y[:, k]) for k in range(len(taus))])
        ent = entries_of(W, taus, th)
        both = np.isfinite(ent) & np.isfinite(ent_ref)
        rows.append(dict(eps=eps, all_branch=bool(P.br.all()), du=float(np.abs(U - Uref).max()),
                         dres=float(np.sqrt(np.mean(np.sum((W @ Phi.T - Wref @ Phi.T) ** 2, axis=1)))),
                         dent=float(np.abs(ent - ent_ref)[both].max()) if both.any() else np.nan, n_only=int(np.sum(np.isfinite(ent) != np.isfinite(ent_ref))),
                         dend=float(np.abs(W[-1] - wlim).max()), secs=time.time() - t0))
    out['rows'] = rows
    # a ± pair per atom: the two-sided path with the weights 1/(γa₀)
    al = 0.4
    Phi2 = np.hstack([Phi, Phi])
    rows2 = []
    for u0 in U_B:
        t0 = time.time()
        bb = u0 * L.m / (L.A * np.tanh(al / L.A))
        P = MF.ParticleSet(L, np.r_[np.full(DA, al), np.full(DA, -al)], np.full(2 * DA, bb))
        a0 = L.A * np.tanh(P.al0[0] / L.A)
        rho2 = np.full(DA, 1 / (L.gamma * a0))
        path2 = MF.iss_path(Phi, y, rho2)
        te = 1.5 * path2['taus'][-1]
        tg = np.linspace(0, te, NT_A)
        Lg = np.log(1 / u0)
        e0 = P.eta(P.x0)
        sol = solve_ivp(P.rhs(Phi2, y), (0, Lg * te), P.x0, method='LSODA', rtol=1e-10, atol=1e-12, t_eval=Lg * tg,
                        max_step=Lg * te / NT_A)
        Wp = np.array([P.w(sol.y[:, k]) for k in range(len(tg))])
        W = Wp[:, :DA] + Wp[:, DA:]
        Uplus = np.array([(P.eta(sol.y[:, k]) - e0)[:DA] / Lg for k in range(len(tg))])   # the atom's common increment / log(1/u)
        Uref2 = np.array([MF.iss_eval(path2, t, Phi, y)[0] for t in tg])
        Wref2 = np.array([MF.iss_eval(path2, t, Phi, y)[1] for t in tg])
        rows2.append(dict(u=u0, du=float(np.abs(Uplus - Uref2).max()),
                          dres=float(np.sqrt(np.mean(np.sum((W @ Phi.T - Wref2 @ Phi.T) ** 2, axis=1)))),
                          dend=float(np.abs(W[-1] - path2['w_end']).max()), n_ent=sum(1 for e in path2['events'] if e[2] == 'enter'),
                          secs=time.time() - t0))
    out['pairs'] = rows2
    return out


def bumps_target(K, cstar, astar):
    return K.phi(np.asarray(cstar, float)) @ np.asarray(astar, float)


def run_moving(L, K, fstar, c0, eps, tau_end, metric, n_t=441, fixed=False):
    """The free law with moving (or fixed) centres from α = 0, masses ε, to τ_end in τ = t/Λ(ε)."""
    d = len(c0)
    Lam = RT.transit_scale(L, eps)
    taus = np.linspace(0, tau_end, n_t)
    z0 = np.r_[np.zeros(d), np.log(np.full(d, eps)), c0]
    if fixed:                                                          # L.rhs: g = Φᵀr/N, N the samples — the same mean
        f0 = L.rhs(K.phi(np.asarray(c0, float)), fstar)
        rhs = lambda t, z: np.r_[f0(t, z[:2 * d]), np.zeros(d)]
    else:
        rhs = RT.moving_rhs(L, K, fstar, metric)
    sol = solve_ivp(rhs, (0, Lam * tau_end), z0, method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
    al, lb, cc = sol.y[:d].T, sol.y[d:2 * d].T, sol.y[2 * d:].T
    W = np.exp(lb) * L.amp(al)[0] / L.m
    return taus, W, cc, Lam


def summarize(K, fstar, c0, taus, W, cc, th):
    act = np.abs(W[-1]) > th
    r = fstar - K.phi(cc[-1]) @ W[-1]
    ent = entries_of(W, taus, th)
    order = [int(i) for i in np.argsort(ent) if np.isfinite(ent[i])]
    cl = np.sort(cc[-1][act])
    clusters = int(np.sum(np.diff(cl) > 1e-3) + 1) if cl.size else 0
    return dict(n_act=int(act.sum()), centres=cc[-1][act].tolist(), w=W[-1][act].tolist(), order=order, ent=[float(ent[i]) for i in order],
                rms=float(np.sqrt(np.mean(r ** 2))), disp_in=float(np.abs(cc[-1] - c0)[~act].max()) if (~act).any() else 0.0,
                w_in=float(np.abs(W[-1][~act]).max()) if (~act).any() else 0.0, clusters=clusters)


def part_C():
    warnings.filterwarnings('ignore')
    L = law()
    K = RT.Bumps(np.linspace(-1, 1, NS_C), SIG_C)
    c0 = np.linspace(-1, 1, D_C)
    fstar = bumps_target(K, CSTAR, ASTAR)
    th = TH_ACT * np.abs(ASTAR).max()
    out = {}
    # (1) the paradigm's metric against the sliding inverse scale space
    R = RT.sliding_iss(L, K, fstar, c0, np.ones(D_C))
    out['slide'] = dict(events=[(float(t), int(j), k, float(c)) for t, j, k, c in R['events']], centres=R['c'][R['act']].tolist(),
                        w=R['w'][R['act']].tolist(), rms=R['hist'][-1][3])
    tau_end = 1.3 * R['events'][-1][0]
    rows = []
    for eps in EPS_C1:
        t0 = time.time()
        taus, W, cc, Lam = run_moving(L, K, fstar, c0, eps, tau_end, 'paradigm')
        rows.append(dict(eps=eps, Lam=float(Lam), secs=time.time() - t0, **summarize(K, fstar, c0, taus, W, cc, th)))
    out['c1'] = rows
    out['dtau1'] = float(tau_end / 440)
    # (2) the Wasserstein–Fisher–Rao metric: condensation before the first entry, against the limit
    lim = solve_ivp(RT.wfr_inactive_rhs(L, K, fstar, np.ones(D_C)), (0, TAU_PRE), np.r_[np.zeros(D_C), c0], method='LSODA', rtol=1e-10,
                    atol=1e-12)
    c_lim = lim.y[D_C:, -1]
    g_peaks = K.phi(c_lim).T @ fstar / K.xs.size
    rows = []
    for eps in EPS_C2:
        t0 = time.time()
        taus, W, cc, Lam = run_moving(L, K, fstar, c0, eps, tau_end, 'wfr')
        k = int(np.searchsorted(taus, TAU_PRE))
        rows.append(dict(eps=eps, d_pre=float(np.abs(cc[k] - c_lim).max()), w_pre=float(np.abs(W[k]).max()), secs=time.time() - t0,
                         **summarize(K, fstar, c0, taus, W, cc, th)))
    cl = np.sort(c_lim)
    out['c2'] = dict(rows=rows, c_lim=c_lim.tolist(), lim_clusters=int(np.sum(np.diff(cl) > 1e-3) + 1))
    # (3) fixed centres (the grid dictionary)
    taus, W, cc, Lam = run_moving(L, K, fstar, c0, EPS_RAND2, tau_end * 3, 'paradigm', fixed=True)
    Phi = K.phi(c0)
    wls = np.linalg.lstsq(Phi, fstar, rcond=None)[0]
    out['c3'] = dict(eps=EPS_RAND2, tau=float(taus[-1]), **summarize(K, fstar, c0, taus, W, cc, th),
                     rms_ls=float(np.sqrt(np.mean((fstar - Phi @ wls) ** 2))))
    # (4) random three-bump targets
    rng = np.random.default_rng(21)
    rnd = []
    for n in range(N_RAND):
        while True:
            cs = np.sort(rng.uniform(-0.8, 0.8, 3))
            if np.diff(cs).min() >= 0.3:
                break
        As = rng.choice([-1, 1], 3) * rng.uniform(0.4, 1.0, 3)
        f_n = bumps_target(K, cs, As)
        th_n = TH_ACT * np.abs(As).max()
        Rn = RT.sliding_iss(L, K, f_n, c0, np.ones(D_C))
        te = 1.3 * Rn['events'][-1][0]
        rec = dict(cs=cs.tolist(), As=As.tolist(), slide_n=int(Rn['act'].sum()), slide_rms=Rn['hist'][-1][3],
                   slide_err=float(np.abs(np.sort(Rn['c'][Rn['act']]) - cs).max()) if Rn['act'].sum() == 3 else np.inf)
        rec['slide_order'] = [int(e[1]) for e in Rn['events'] if e[2] == 'enter']
        for key, metric, eps_n in (('paradigm', 'paradigm', EPS_RAND), ('paradigm2', 'paradigm', EPS_RAND2),
                                   ('paradigm3', 'paradigm', EPS_RAND3), ('wfr', 'wfr', EPS_RAND2)):
            taus, W, cc, Lam = run_moving(L, K, f_n, c0, eps_n, te, metric)
            sm = summarize(K, f_n, c0, taus, W, cc, th_n)
            rec[key] = dict(n_act=sm['n_act'], clusters=sm['clusters'], rms=sm['rms'], disp_in=sm['disp_in'], order=sm['order'],
                            err=float(np.abs(np.sort(sm['centres']) - cs).max()) if sm['n_act'] == 3 else np.inf)
        taus, W, cc, Lam = run_moving(L, K, f_n, c0, EPS_RAND2, 3 * te, 'paradigm', fixed=True)
        sm = summarize(K, f_n, c0, taus, W, cc, th_n)
        rec['fixed'] = dict(n_act=sm['n_act'], rms=sm['rms'])
        rnd.append(rec)
        print(f"  C4 target {n}: slide {rec['slide_n']} paradigm {rec['paradigm']['n_act']}, {rec['paradigm2']['n_act']}, "
              f"{rec['paradigm3']['n_act']} wfr {rec['wfr']['n_act']} fixed {rec['fixed']['n_act']}", flush=True)
    out['c4'] = rnd
    # (5) two close bumps; a target wider than the kernels
    hard = []
    for lab, cs, As, wid in (('two bumps 0.2 apart', (-0.1, 0.1), (0.8, 0.6), SIG_C), ('two bumps 0.3 apart', (-0.15, 0.15), (0.8, 0.6), SIG_C),
                             ('three wider bumps (width 0.2)', CSTAR, ASTAR, 0.2)):
        Kt = RT.Bumps(K.xs, wid)
        f_h = bumps_target(Kt, cs, As)
        Rh = RT.sliding_iss(L, K, f_h, c0, np.ones(D_C), max_events=30)
        te = 1.3 * Rh['events'][-1][0]
        taus, W, cc, Lam = run_moving(L, K, f_h, c0, EPS_RAND, te, 'paradigm')
        sm = summarize(K, f_h, c0, taus, W, cc, TH_ACT * np.abs(As).max())
        hard.append(dict(lab=lab, slide_n=int(Rh['act'].sum()), slide_c=Rh['c'][Rh['act']].tolist(), slide_rms=Rh['hist'][-1][3],
                         slide_events=len(Rh['events']), slide_done=bool(Rh['hist'][-1][3] < 1e-9),
                         **{k: sm[k] for k in ('n_act', 'centres', 'w', 'rms', 'clusters')}))
    out['c5'] = hard
    # (6) a lost near-tie: the target C4's CREEP_TARGET, the kernel left just inside its wall creeps at (πa_max²/(4c))εμ^{−3}⟨r∂_cφ⟩
    rc = rnd[CREEP_TARGET]
    cs, As = np.array(rc['cs']), np.array(rc['As'])
    f_c = bumps_target(K, cs, As)
    Rc = RT.sliding_iss(L, K, f_c, c0, np.ones(D_C))
    te = 1.3 * Rc['events'][-1][0]
    creep = dict(slide_order=[int(e[1]) for e in Rc['events'] if e[2] == 'enter'], slide_tau=[float(e[0]) for e in Rc['events']], rows=[])
    t2 = Rc['events'][1][0]; t3 = Rc['events'][2][0]
    for eps in EPS_CREEP:
        Lam = RT.transit_scale(L, eps)
        taus = np.linspace(0, te, 881)
        sol = solve_ivp(RT.moving_rhs(L, K, f_c, 'paradigm'), (0, Lam * te), np.r_[np.zeros(D_C), np.log(np.full(D_C, eps)), c0],
                        method='LSODA', rtol=1e-10, atol=1e-13, t_eval=Lam * taus)
        al, cc = sol.y[:D_C].T, sol.y[2 * D_C:].T
        W = np.exp(sol.y[D_C:2 * D_C].T) * L.amp(al)[0] / L.m
        q = L.Q(np.zeros(D_C), np.full(D_C, eps))
        lam1 = np.pi * L.A * np.sqrt((L.c - q) / (q + L.c)) / (4 * L.G * L.c)
        ent = entries_of(W, taus, TH_ACT * np.abs(As).max())
        order = [int(i) for i in np.argsort(ent) if np.isfinite(ent[i])]
        k2 = int(np.searchsorted(taus, t2 + 0.2)); kb = int(np.searchsorted(taus, min(t3, ent[order[2]]) - 0.1))
        mu = 1 - np.abs(L.eta(al[k2], q)) / lam1
        r2 = f_c - K.phi(cc[k2]) @ W[k2]
        h2 = K.dphi(cc[k2]).T @ r2 / K.xs.size
        near = [int(j) for j in np.flatnonzero((np.abs(W[k2]) < 1e-3) & (mu < 0.01))]   # the inactive kernels lingering at their wall
        creep['rows'].append(dict(eps=eps, order=order, centres=np.sort(cc[-1][np.abs(W[-1]) > TH_ACT * np.abs(As).max()]).tolist(),
                                  near=[dict(j=j, mu=float(mu[j]), h=float(h2[j]), disp=float(abs(cc[kb, j] - c0[j])),
                                             speed=float(np.pi * L.A ** 2 / (4 * L.c) * eps * mu[j] ** -3 * abs(h2[j]))) for j in near]))
    out['c6'] = creep
    return out


def report(A=None, B=None, C=None):
    fmt = E.fmt
    L = law()
    L_ = ["# Proposition D‴.48: the rich transit (`exp_mirror_transit.py`)", "",
          f"The free mass–amplitude law (Proposition D‴.45; m = {LAW['m']:g}, G = {LAW['G']:g}, γ = {LAW['gamma']:g}, a_max = {LAW['a_max']:g}, "
          f"c = {L.c:g}) from small masses, in the time rescaled by the transit's scale (`mirror_free.iss_path`, `rich_transit`).", ""]
    if A is not None:
        L_ += ["## A. From α = 0: the inverse scale space flow on the scale Λ(ε) = (πa_max/(4G))(2/(cε))^{1/2}", "",
               f"Three {NA} × {DA} Gaussian designs, 3-sparse truths, masses εb_i (b_i log-uniform in [1/4, 4]); the path of Σ|w_i|/√b_i against "
               "the loop in τ = t/Λ, u = η/Λ, on 2001 times in [0, 1.5τ*] (τ* the last event): sup_τ|u − u_ISS|, the residual's rms difference "
               "over τ, the largest difference of the entry times (the first sample where |w_i| > 10⁻³ max|w_bp|, on the loop and on the path; "
               "in brackets the atoms crossing on one of them only; the grid's step given), the end's distance to weighted basis pursuit.", ""]
        for a in A:
            ev = '; '.join(f"{k} {i} at {t:.4g}" for t, i, k in a['events'])
            L_ += [f"Design {a['seed']}: the path's events — {ev}; its end against weighted basis pursuit {fmt(a['end_is_bp'])}; the grid's step "
                   f"{fmt(a['dtau'])}.", "",
                   "| ε | Λ | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_bp|∞ |", "|---|---|---|---|---|---|"]
            L_ += [f"| {r['eps']:g} | {r['Lam']:.4g} | {fmt(r['du'])} | {fmt(r['dres'])} | {fmt(r['dent'])} ({r['n_only']}) | {fmt(r['dend'])} |"
                   for r in a['rows']]
            du = {r['eps']: r['du'] for r in a['rows']}
            L_ += ["", f"sup|u − u_ISS| between ε = 10⁻⁶ and 10⁻⁸: the slope {np.log(du[1e-6] / du[1e-8]) / np.log(100):.3f} in log ε (1/Λ ∝ ε^{{1/2}}: "
                   "0.5).", ""]
    if B is not None:
        L_ += ["## B. On branches: the logarithmic scale", "",
               f"Design {SEEDS_A[0]} of A; amplitudes s_i·U(0.2, 0.6) (the truth's signs on its atoms, random elsewhere), masses ε; the one-sided "
               "path (frozen signs, weights 1/(γ|a_i(0)|)) against the loop in τ = t/log(1/ε), u = (η − η(0))/log(1/ε). The path's events: "
               + '; '.join(f"{k} {i} at {t:.4g}" for t, i, k in B['events']) + f"; its end against the sign-constrained weighted basis pursuit "
               f"{fmt(B['end_is_lim'])}.", "",
               "| ε | all on branches | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_lim|∞ |",
               "|---|---|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | {'yes' if r['all_branch'] else 'no'} | {fmt(r['du'])} | {fmt(r['dres'])} | {fmt(r['dent'])} ({r['n_only']}) | "
               f"{fmt(r['dend'])} |" for r in B['rows']]
        L_ += ["", "A ± pair per atom (amplitudes ±0.4, equal masses, w(0) = u): the atom's common increment over log(1/u) against the two-sided "
               "path with the weights 1/(γa₀):", "", "| u | sup|u − u_ISS| | residual: rms difference | |W(1.5τ*) − w_end|∞ | entries on the path |",
               "|---|---|---|---|---|"]
        L_ += [f"| {r['u']:g} | {fmt(r['du'])} | {fmt(r['dres'])} | {fmt(r['dend'])} | {r['n_ent']} |" for r in B['pairs']]
        L_.append("")
    if C is not None:
        sl = C['slide']
        L_ += ["## C. Moving centres", "",
               f"Gaussian kernels of width {SIG_C} on {NS_C} samples of [−1, 1]; {D_C} kernels started on a grid of [−1, 1] (spacing 0.1), α = 0, "
               f"masses ε; the target {len(CSTAR)} bumps of the same width at {list(CSTAR)} with amplitudes {list(ASTAR)}. A kernel has formed when "
               f"|w| > {TH_ACT:g}·max|A*|.", "",
               "(1) The paradigm's metric (ċ = Gβa⟨r∂_cφ⟩). The sliding inverse scale space: " +
               '; '.join(f"{k} kernel {j} (start {c:+.1f}) at τ = {t:.4g}" for t, j, k, c in sl['events']) +
               f"; it ends with {len(sl['centres'])} kernels at {[round(v, 5) for v in sl['centres']]} with coefficients {[round(v, 5) for v in sl['w']]} "
               f"(rms residual {fmt(sl['rms'])}). The loop (the entries at the first sample — step {C['dtau1']:.3g} in τ — where |w| > {TH_ACT:g}·max|A*|):", "",
               "| ε | Λ | entries (kernel: τ) | formed | their centres | their coefficients | rms residual | inactive: max |Δc|, max |w| |",
               "|---|---|---|---|---|---|---|---|"]
        for r in C['c1']:
            ent = ', '.join(f"{i}: {t:.4g}" for i, t in zip(r['order'], r['ent']))
            L_.append(f"| {r['eps']:g} | {r['Lam']:.4g} | {ent} | {r['n_act']} | {[round(v, 4) for v in r['centres']]} | {[round(v, 4) for v in r['w']]} | "
                      f"{fmt(r['rms'])} | {fmt(r['disp_in'])}, {fmt(r['w_in'])} |")
        c2 = C['c2']
        L_ += ["", f"(2) The Wasserstein–Fisher–Rao metric (ċ = Ga⟨r∂_cφ⟩): before the first entry (τ = {TAU_PRE:g}) the inactive kernels' limit "
               f"(du/dτ = ⟨f*φ⟩, dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ⟩) gathers the {D_C} kernels at {c2['lim_clusters']} positions "
               f"{sorted(set(round(v, 3) for v in c2['c_lim']))}; the loop's centres at τ = {TAU_PRE:g} against that limit, and its end:", "",
               "| ε | max |c − c_lim| at τ = 10 | max |w| at τ = 10 | formed | their positions (clusters) | rms residual |", "|---|---|---|---|---|---|"]
        for r in c2['rows']:
            L_.append(f"| {r['eps']:g} | {fmt(r['d_pre'])} | {fmt(r['w_pre'])} | {r['n_act']} | {sorted(set(round(v, 3) for v in r['centres']))} "
                      f"({r['clusters']}) | {fmt(r['rms'])} |")
        c3 = C['c3']
        L_ += ["", f"(3) The centres held on the grid (ε = {c3['eps']:g}, to τ = {c3['tau']:.3g}): {c3['n_act']} kernels formed, rms residual "
               f"{fmt(c3['rms'])} (the grid's least-squares fit: {fmt(c3['rms_ls'])}).", ""]
        c4 = C['c4']
        n = len(c4)
        ok_s = sum(1 for r in c4 if r['slide_n'] == 3 and r['slide_err'] < 1e-3)
        ok = {k: sum(1 for r in c4 if r[k]['n_act'] == 3 and r[k]['err'] < 1e-3) for k in ('paradigm', 'paradigm2')}
        L_ += [f"(4) {n} random targets (three bumps of the kernels' width, centres in [−0.8, 0.8] at least 0.3 apart, amplitudes ±U(0.4, 1)): "
               f"the sliding inverse scale space ends with exactly three kernels within 10⁻³ of the bumps' centres on {ok_s} of {n}; the loop with "
               f"the paradigm's metric on {ok['paradigm']} of {n} at ε = {EPS_RAND:g} (formed kernels {sorted(r['paradigm']['n_act'] for r in c4)}, "
               f"the inactive kernels' largest displacement {fmt(max(r['paradigm']['disp_in'] for r in c4))}) and on {ok['paradigm2']} of {n} at "
               f"ε = {EPS_RAND2:g} (formed kernels {sorted(r['paradigm2']['n_act'] for r in c4)}, the inactive kernels' largest displacement "
               f"{fmt(max(r['paradigm2']['disp_in'] for r in c4))}; largest rms residual {fmt(max(r['paradigm2']['rms'] for r in c4))}) and on "
               f"{sum(1 for r in c4 if r['paradigm3']['n_act'] == 3 and r['paradigm3']['err'] < 1e-3)} of {n} at ε = {EPS_RAND3:g} (the entering "
               f"kernels those of the sliding path on {sum(1 for r in c4 if r['paradigm3']['order'] == r['slide_order'])} of {n}); with the "
               f"Wasserstein–Fisher–Rao metric (ε = {EPS_RAND2:g}) it forms {sorted(r['wfr']['n_act'] for r in c4)} kernels in "
               f"{sorted(r['wfr']['clusters'] for r in c4)} clusters (largest rms residual {fmt(max(r['wfr']['rms'] for r in c4))}); with the centres "
               f"held on the grid {sorted(r['fixed']['n_act'] for r in c4)} kernels (rms residual {fmt(min(r['fixed']['rms'] for r in c4))}–"
               f"{fmt(max(r['fixed']['rms'] for r in c4))}).", ""]
        L_ += ["(5) Harder targets (ε = 10⁻⁵, the paradigm's metric):", "",
               "| target | sliding ISS: kernels, their centres, rms residual (events) | loop: formed, centres, coefficients (clusters) | rms residual |",
               "|---|---|---|---|"]
        for r in C['c5']:
            L_.append(f"| {r['lab']} | {r['slide_n']}, {[round(v, 3) for v in r['slide_c']]}, {fmt(r['slide_rms'])} ({r['slide_events']}"
                      f"{'' if r['slide_done'] else ', not terminated'}) | "
                      f"{r['n_act']}, {[round(v, 3) for v in r['centres']]}, {[round(v, 3) for v in r['w']]} ({r['clusters']}) | {fmt(r['rms'])} |")
        L_.append("")
    if C is not None and 'c6' in C:
        c6 = C['c6']
        L_ += ["", f"(6) A lost near-tie (target {CREEP_TARGET} of (4)): the sliding path enters kernels {c6['slide_order']} at τ = "
               f"{[round(t, 3) for t in c6['slide_tau']]}. After the second entry the inactive kernels left near their walls linger at their margins μ "
               "and creep at ≈ (πa_max²/(4c))εμ^{−3}|⟨r∂_cφ⟩| in τ ((iii)(a)'s bound):", "",
               "| ε | the lingering kernels (margin < 0.01 after the second entry): kernel, μ, ⟨r∂_cφ⟩, the speed (πa_max²/(4c))εμ^{−3}|⟨r∂_cφ⟩|, the "
               "displacement before the third entry | the loop's entrants | the formed kernels' centres |", "|---|---|---|---|"]
        L_ += [f"| {r['eps']:g} | " + '; '.join(f"{n['j']}: {fmt(n['mu'])}, {fmt(n['h'])}, {fmt(n['speed'])}, {fmt(n['disp'])}" for n in r['near'])
               + f" | {r['order']} | {[round(v, 4) for v in r['centres']]} |" for r in c6['rows']]
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_mirror_transit.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_transit_raw.npz')
    md = os.path.join(OUT, 'mirror_transit.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABC']))
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
