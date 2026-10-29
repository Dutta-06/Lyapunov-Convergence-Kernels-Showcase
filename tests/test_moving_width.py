"""Proposition D‴.59: the transit with moving widths (`src/core/moving_width.py`)."""
import os
import sys

import numpy as np

from src.core import moving_width as MW
from src.core import rich_transit_nd as RN

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_moving_width as MWE  # noqa: E402


def test_bumps_cw_derivatives_and_grams():
    """∂_cφ, ∂_σφ and ⟨r∇²_θφ⟩ against central differences of φ; the Gram blocks against their sums."""
    rng = np.random.default_rng(0)
    K = MW.BumpsCW(np.linspace(-1.5, 1.5, 201), 0.12)
    th = np.c_[rng.uniform(-0.5, 0.5, 4), rng.uniform(0.05, 0.3, 4)]
    r = rng.standard_normal(K.N)
    h = 1e-6
    for k in range(2):
        e = np.zeros(2); e[k] = h
        fd = (K.phi(th + e) - K.phi(th - e)) / (2 * h)
        assert np.abs(fd - K.dphi(th)[:, :, k]).max() < 1e-6 * np.abs(fd).max()
    g, hh, Q = K.corr(th, r, second=True)
    assert np.allclose(g, K.phi(th).T @ r / K.N, rtol=1e-12, atol=0)
    assert np.allclose(hh, np.einsum('s,sjk->jk', r, K.dphi(th)) / K.N, rtol=1e-12, atol=1e-15)
    for k in range(2):
        e = np.zeros(2); e[k] = h
        fd = (K.corr(th + e, r)[1] - K.corr(th - e, r)[1]) / (2 * h)
        assert np.abs(fd - Q[:, :, k]).max() < 1e-5 * np.abs(Q).max()
    P, D, H = K.grams(th)
    dP = K.dphi(th)
    assert np.allclose(P, K.phi(th).T @ K.phi(th) / K.N) and np.allclose(D[1, 2], K.phi(th)[:, 1] @ dP[:, 2, :] / K.N)
    assert np.allclose(H[0, 3], dP[:, 0, :].T @ dP[:, 3, :] / K.N)


def test_the_one_kernel_reduced_energy():
    """`overlap` against quadrature; the reduced energy vanishes at (c*, σ*) only and its closed-form gradient matches central differences."""
    x = np.linspace(-12, 12, 240001)
    dx = x[1] - x[0]
    q = np.exp(-(x - 0.3) ** 2 / (2 * 0.2 ** 2)) @ np.exp(-(x + 0.1) ** 2 / (2 * 0.45 ** 2)) * dx
    assert abs(MW.overlap(0.3, 0.2, -0.1, 0.45) / q - 1) < 1e-10
    A, cs, ss = 0.8, 0.03, 0.2
    assert abs(MW.reduced_energy(cs, ss, A, cs, ss)) < 1e-15
    rng = np.random.default_rng(1)
    for _ in range(10):
        c, s = rng.normal(0, 0.3), ss * np.exp(rng.normal(0, 0.7))
        assert MW.reduced_energy(c, s, A, cs, ss) > 0
        h = 1e-6
        fd = np.array([(MW.reduced_energy(c + h, s, A, cs, ss) - MW.reduced_energy(c - h, s, A, cs, ss)) / (2 * h),
                       (MW.reduced_energy(c, s + h, A, cs, ss) - MW.reduced_energy(c, s - h, A, cs, ss)) / (2 * h)])
        g = MW.one_kernel_grad(c, s, A, cs, ss)
        assert np.abs(fd - g).max() < 1e-6 * max(np.abs(g).max(), 1e-12)


def test_a_narrow_feature_is_formed_by_one_kernel_with_moving_widths():
    """At a quarter of the kernels' width the sliding path with moving widths enters one kernel and relaxes it to the feature exactly
    (centre, width, amplitude), within a small budget; with the widths fixed the same budget runs out."""
    import exp_moving_width as M
    L = M.law()
    ss = 0.25 * M.S0
    f = M.feature_sum([M.C_A], [ss], [M.A_A])
    Rm = RN.sliding_iss_nd(L, MW.BumpsCW(M.XS, M.S0), f, M.theta0(), np.ones(M.C0.size), max_events=40, exit_kkt=True, max_rhs=50_000)
    m = M.match(Rm, [M.C_A], [ss], [M.A_A])
    assert m['n'] == 1 and m['exact'] and not m['oob'] and m['entries'] == 1
    Rf = RN.sliding_iss_nd(L, RN.BumpsND(M.XS[:, None], M.S0), f, M.C0[:, None], np.ones(M.C0.size), max_events=40,
                           max_entries=M.MAX_ENT_FIXED, exit_kkt=True, max_rhs=50_000)
    assert Rf['out_of_budget'] and not M.match(Rf, [M.C_A], [ss], [M.A_A])['exact']


def test_part_b_gradient_check_is_well_conditioned():
    """Part B checks the closed-form gradient against differences of the overlap term alone: in the far tail (c = 6.25σ*, where the
    gradient is ~10⁻¹⁵ and differences of J are swamped by its constant) the relative error stays at the differences' own level."""
    import exp_moving_width as M
    B = M.part_B()
    assert max(r['grad_err'] for r in B) < 1e-5
    assert all(len(r['minima']) == 1 and abs(r['minima'][0][0]) < 1e-12 for r in B)


def test_outcome_classes_and_pairs():
    """`pairs_of` finds kernels at one (c, σ) within PAIR_TOL of their width; `outcome` sorts a path's end into exact (idle entrants
    allowed), a coinciding pair, stopped by the budget, or not fitted."""
    import exp_moving_width as M
    th = [[0.35001, 0.2731707], [0.35097, 0.2731707], [-0.6, 0.08], [0.2, 0.1]]
    assert M.pairs_of(th, [0.33, 0.44, 0.6, -0.8]) == [(0.33, 0.44)] and M.pairs_of(th[2:], [0.6, -0.8]) == []
    As = [0.77, 0.6, -0.8]
    base = dict(rms=1e-13, err=[1e-9, 1e-9, 1e-9], oob=False)
    assert M.outcome(dict(base, theta=th[1:], w=[0.77, 0.6, -0.8]), As) == 'exact'
    assert M.outcome(dict(base, theta=th + [[0.9, 0.1]], w=[0.33, 0.44, 0.6, -0.8, 1e-6], err=[0.33, 1e-9, 1e-9]), As) == 'coinciding pair'
    assert M.outcome(dict(base, theta=th[1:], w=[0.77, 0.6, -0.8], oob=True), As) == 'stopped by the budget'
    assert M.outcome(dict(base, theta=th[1:], w=[0.77, 0.6, -0.8], rms=0.2), As) == 'not fitted'


def test_the_g_report():
    """G's tables (referee round 57): the stalled block (kernels, max|w|, distance, residual against the dipole infimum, its polishes);
    the flow's largest |w| and distance; the last interval's exponents; K = v²(E − E_inf) at the last two outputs; Δv³/Δs against
    3μ₀K; the moment's σ-share; the loop's snapshots with the closest pair's distance, coefficients and offset."""
    import exp_moving_width as M
    ws = (0.9, 4.6, 20.0, 92.0, 300.0, 425.0)
    seps = (5.9, 0.29, 0.062, 0.013, 0.0042, 0.0029)
    flow = [dict(s=s_, wmax=w, sep=d, rms=0.2 + 1 / w) for s_, w, d in zip(M.S_G, ws, seps)]
    law = [dict(v=v, E=E, E_inf=0.02, rms_inf=0.2, p=[0.1, -0.35], mu0=0.1, dc=dc, ds=ds)
           for v, E, dc, ds in ((300.0, 0.02 + 1e-5, 2e-4, 6e-4), (425.0, 0.02 + 5e-6, 1.4e-4, 4.2e-4))]
    stalled = dict(n=4, wmax=3910.0, sep=1.9e-4, rms=0.2000001, polish=[dict(status=2, nfev=2900, accepted=True, wmax_ratio=4300.0)], nfev=10)
    cs, ss, As = [-0.3, 0.2, 0.5], [0.1, 0.15, 0.08], [0.9, 0.54, -0.6]
    pair = lambda off, d: [[0.2 + off * 0.15 + 3.96 * d * 0.15 / 0.54, 0.15], [0.2 + off * 0.15 + 4.5 * d * 0.15 / 0.54, 0.15]]
    loops = [dict(eps=1e-8, tau=t, Lam=100.0, theta=pair(o, d) + [[-0.3, 0.1], [0.5, 0.08]], w=[4.5, -3.96, 0.9, -0.6], wmax=4.5, rms=r)
             for t, d, o, r in ((400.0, 4.6e-4, 0.008, 1.1e-11), (1600.0, 3e-4, 0.005, 5e-12), (6400.0, 2e-4, 0.002, 1e-12))]
    G = [dict(dmin=0.4, q=9, cs=cs, ss=ss, As=As, entries=2, tau_d=17.0, stalled=stalled, flow=flow, law=law, status=0, loops=loops)]
    md = M.report(G=G)
    ls = np.log(M.S_G[-1] / M.S_G[-2])
    ex = M.dipole_exponents(flow)
    exs = (f"{ex[0]:+.4f}, {ex[1]:+.4f}, {np.log(1.4e-4 / 2e-4) / ls:+.4f}, {np.log(4.2e-4 / 6e-4) / ls:+.4f}, "
           f"{np.log((flow[-1]['rms'] - 0.2) / (flow[-2]['rms'] - 0.2)) / ls:+.4f}")
    K2 = [300.0 ** 2 * 1e-5, 425.0 ** 2 * 5e-6]
    vv = f"{(425.0 ** 3 - 300.0 ** 3) / (M.S_G[-1] - M.S_G[-2]):.5e}, {3 * 0.1 * K2[1]:.5e}"
    assert (f"| 0.4 | 9 | 2 | 4, 3910, 1.90e-04, 2.000001e-01 (2.000000e-01), 2/2900/a/4.3e+03 | 0.9, 4.6, 20, 92, 300, 425 | "
            f"5.90e+00, 2.90e-01, 6.20e-02, 1.30e-02, 4.20e-03, 2.90e-03 | {exs} | {K2[0]:.5e}, {K2[1]:.5e} | {vv} | "
            f"{0.35 / np.hypot(0.1, 0.35):.3f} |") in md
    lp = M.landing(loops[1]['theta'], loops[1]['w'], cs, ss, As)                         # the dipole pair at τ = 1600 (round 61)
    assert abs(lp['off'] - 0.005) < 1e-12 and abs(lp['dW']) < 1e-12 and lp['w'] == (4.5, -3.96)
    assert (f"| 0.4 | 9 | 1e-08 | 1600 | 4 | {lp['d']:.2e}, +4.500/-3.960, 5.0e-03, {lp['dW']:+.1e} L | +4.500/-3.960 | 5.00 | "
            "5.0e-12 |") in md and abs(lp['d'] - 3e-4) < 1e-12
    ex = M.dipole_exponents([dict(s=1e9, wmax=1.0, sep=1.0), dict(s=1e11, wmax=100 ** (1 / 3), sep=100 ** (-1 / 3))])
    assert abs(ex[0] - 1 / 3) < 1e-12 and abs(ex[1] + 1 / 3) < 1e-12 and M.dipole_exponents(flow[:1] + [dict(flow[1], sep=None)]) is None


def test_the_formed_kernels_pairs_and_outcomes():
    """Referee round 57: `n_formed` counts |w| above 1% of the smallest |amplitude|; `pairs_of` at PAIR_TOL in |Δc|/σ + |Δlog σ|;
    `outcome`'s classes at their thresholds — exact (idle entrants allowed), a coinciding pair, a cancellation (a coefficient beyond 1.5
    times the largest amplitude, or a width below the samples' spacing), not fitted (an error 2·10⁻⁶ or a residual 2·10⁻⁸), a dipole
    run-off or another budget stop (the stalled block)."""
    import exp_moving_width as M
    assert M.n_formed(dict(w=[0.009, 0.0095, -0.5]), [1.0, -0.9]) == 2
    T = M.PAIR_TOL
    assert M.pairs_of([[0.0, 0.1], [0.1 * T * 0.99, 0.1]], [1.0, 2.0]) == [(1.0, 2.0)]
    assert M.pairs_of([[0.0, 0.1], [0.1 * T * 1.01, 0.1]], [1.0, 2.0]) == []
    assert M.pairs_of([[0.0, 0.1], [0.0, 0.1 * np.exp(T * 0.99)]], [1.0, 2.0]) and not M.pairs_of([[0.0, 0.1], [0.0, 0.1 * np.exp(T * 1.01)]], [1, 2])
    th = [[-0.5, 0.1], [0.0, 0.1], [0.5, 0.1]]
    m = lambda **kw: dict(dict(w=[1.0, 1.0, 1.0], theta=th, rms=1e-9, err=[1e-7, 1e-7, 1e-7], oob=False), **kw)
    As = [1.0, 1.0, 1.0]
    assert M.outcome(m(), As) == 'exact'
    assert M.outcome(m(w=[1.0, 1.0, 1.0, 1e-3, 1e-3], theta=th + [[0.9, 0.1], [0.9, 0.1]]), As) == 'exact'
    assert M.outcome(m(err=[2e-6, 1e-7, 1e-7]), As) == 'not fitted' and M.outcome(m(rms=2e-8), As) == 'not fitted'
    assert M.outcome(m(w=[1.0, 0.5, 0.5, 1.0], theta=[th[0], th[1], th[1], th[2]], err=[1e-7, 0.5, 1e-7]), As) == 'coinciding pair'
    assert M.outcome(m(w=[1.6, 1.0, 1.0], err=[0.6, 1e-7, 1e-7]), As) == 'cancellation'
    assert M.outcome(m(theta=[[-0.5, 0.005], th[1], th[2]], err=[0.1, 1e-7, 1e-7]), As) == 'cancellation'
    st = dict(theta=[[0.0, 0.1], [0.0002, 0.1], [0.5, 0.1]], w=[500.0, -499.0, 1.0])
    assert M.outcome(m(oob=True, stalled=st), As) == 'dipole run-off'
    assert M.outcome(m(oob=True, stalled=dict(st, w=[5.0, -4.0, 1.0])), As) == 'stopped by the budget'
    assert M.outcome(m(oob=True), As) == 'stopped by the budget'


def test_the_width_split_cancels_the_third_moment():
    """Referee round 57: a pair of unequal shares Δc apart fits one Gaussian to O(Δc⁴) when its widths split by (w_j − w_i)Δc²/(6Wσ);
    with equal widths the best single Gaussian leaves the third moment, O(Δc³)."""
    import exp_moving_width as M
    from scipy.optimize import least_squares
    x = np.linspace(-1, 1, 801)
    g = lambda c, s: np.exp(-(x - c) ** 2 / (2 * s ** 2))
    wi, wj, sb = 0.3, 0.5, 0.1
    W = wi + wj
    res = []
    for dc in (0.02, 0.01):
        ci, cj = -wj * dc / W, wi * dc / W
        _, _, pred = M.width_split([[ci, sb], [cj, sb]], [wi, wj], 0, 1)
        assert abs(pred - (wj - wi) * dc ** 2 / (6 * W * sb)) < 1e-15
        for d in (0.0, pred):
            f = wi * g(ci, sb - wj * d / W) + wj * g(cj, sb + wi * d / W)
            sol = least_squares(lambda z: f - z[0] * g(z[1], z[2]), [W, 0.0, sb], xtol=1e-15, ftol=1e-15, gtol=1e-15)
            res.append(np.sqrt(np.mean(sol.fun ** 2)))
    (e0, e1), (h0, h1) = res[:2], res[2:]
    assert e1 < 0.2 * e0 and 7 < e0 / h0 < 9 and 14 < e1 / h1 < 18               # O(Δc³) and O(Δc⁴): 8 and 16 per halving
    dc_, ds_, pr_ = M.width_split([[0.01, 0.11], [-0.01, 0.1]], [0.5, 0.3], 0, 1)        # sorted by centre: i the left kernel
    assert dc_ == 0.02 and abs(ds_ - 0.01) < 1e-15 and abs(pr_ - (0.5 - 0.3) * 0.02 ** 2 / (6 * 0.8 * 0.105)) < 1e-15


def test_the_dipole_law_finds_the_infimum():
    """Referee round 57: `dipole_law` on a pair ±v coalescing around a dipole p·∇_θφ beside two kernels: the pair (the largest |w| and its
    nearest kernel of the other sign), v, the moment, and E_inf (the dipole fit's) ≈ 0 for a target that is that dipole exactly (the pair's
    own error the third-order term p³/(24v²)·∇³φ)."""
    import exp_moving_width as M
    L = M.law()
    K = M.MW.BumpsCW(M.XS, M.S0)
    tb, p = np.array([0.1, 0.12]), np.array([0.002, -0.004])
    f = K.synth(tb[None, :], np.array([0.5])) + K.dphi(tb[None, :])[:, 0, :] @ p
    f = f + K.synth(np.array([[0.6, 0.1]]), np.array([0.4]))
    v = 10.0
    th = np.array([tb + p / (2 * v), tb - p / (2 * v), tb, [0.6, 0.1]])
    d = M.dipole_law(L, K, f, th, np.array([v, -v, 0.5, 0.4]))
    assert d['v'] == v and np.allclose(d['p'], p, rtol=1e-6) and d['E_inf'] < 1e-25 < d['E'] and d['mu0'] == 4 * L.G * L.c / L.A


def test_the_outcomes_at_their_thresholds():
    """Referee round 61: `outcome` with the features — exact with an idle null cluster (each feature's nearest kernel exact, the others
    summing to nothing); a coinciding pair only where its net coefficient and centre of mass fit a feature (else a null cluster), among
    formed kernels (|w| above 1% of the smallest amplitude) and only below the residual 1e-8; a budget stop a dipole run-off only for an
    opposite-sign pair closer than 0.05 with |w| beyond 10 times the largest amplitude."""
    import exp_moving_width as M
    cs, ss, As = [-0.5, 0.0, 0.5], [0.1, 0.1, 0.1], [1.0, -0.8, 0.6]
    th = [[c, s] for c, s in zip(cs, ss)]
    m = lambda **kw: dict(dict(w=list(As), theta=th, rms=1e-10, err=[1e-8] * 3, oob=False), **kw)
    oc = lambda mm: M.outcome(mm, As, cs, ss)
    assert oc(m()) == 'exact'
    nul = [[0.9, 0.05], [0.9, 0.05], [0.9, 0.05]]
    assert oc(m(w=list(As) + [-0.05, 0.1, -0.05], theta=th + nul)) == 'exact, null cluster'
    assert oc(m(w=list(As) + [-0.005, 0.005], theta=th + nul[:2])) == 'exact'                     # below 1%: idle
    pr = dict(w=[1.0, -0.3, -0.5, 0.6], theta=[th[0], th[1], [0.0, 0.1 * (1 + 1e-3)], th[2]], err=[1e-8, 0.5, 1e-8])
    assert oc(m(**pr)) == 'coinciding pair'
    assert oc(m(**dict(pr, w=[1.0, -0.3, -0.6, 0.6]))) == 'null cluster'                          # the pair's W −0.9 against −0.8
    assert oc(m(**dict(pr, theta=[th[0], [0.002, 0.1], [0.002, 0.1], th[2]]))) == 'null cluster'   # its centre 0.02σ* off
    assert oc(m(**dict(pr, rms=2e-8))) == 'not fitted'
    assert oc(m(**dict(pr, w=[1.0, -0.795, -0.005, 0.6]))) == 'not fitted'                       # the second kernel unformed
    st = dict(theta=[[0.0, 0.1], [0.004, 0.1], [0.5, 0.1]], w=[50.0, -49.0, 0.6])                # 0.04 apart, 50 = 50·max|A|
    assert oc(m(oob=True, stalled=st)) == 'dipole run-off'
    assert oc(m(oob=True, stalled=dict(st, theta=[[0.0, 0.1], [0.006, 0.1], [0.5, 0.1]]))) == 'stopped by the budget'
    assert oc(m(oob=True, stalled=dict(st, w=[9.9, -8.9, 0.6]))) == 'stopped by the budget'
    st2 = dict(theta=[[0.0, 0.1], [0.001, 0.1], [0.3, 0.1]], w=[50.0, 49.0, -98.0])                # the close pair of one sign
    assert oc(m(oob=True, stalled=st2)) == 'stopped by the budget'


def test_the_closest_pair_counts_the_widths():
    """`closest_pair`'s distance |Δc|/σ + |Δlog σ| (σ the mean width): two kernels at one centre are |log(σ_j/σ_i)| apart."""
    import exp_moving_width as M
    th, w = [[0.0, 0.1], [0.0, 0.2], [0.05, 0.1]], [1.0, -1.0, 1.0]
    d, i, j = M.closest_pair(th, w)
    assert (i, j) == (0, 2) and abs(d - 0.5) < 1e-12
    d, i, j = M.closest_pair(th, w, opposite=True)
    assert (i, j) == (0, 1) and abs(d - np.log(2)) < 1e-12


def test_the_dipole_law_picks_the_nearest_opposite_kernel_and_halves_the_energy():
    """Referee round 61: `dipole_law`'s pair is the largest |w| with the nearest kernel of the other sign (not a farther one); its E and
    E_inf are ½⟨r²⟩ — with a target the dipole model cannot fit (a narrow bump far away), E_inf = ½⟨g²⟩ of that bump."""
    import exp_moving_width as M
    L = M.law()
    K = M.MW.BumpsCW(M.XS, M.S0)
    tb, p, v = np.array([0.1, 0.12]), np.array([0.002, -0.004]), 10.0
    dip = K.synth(tb[None, :], np.array([0.5])) + K.dphi(tb[None, :])[:, 0, :] @ p
    th = np.array([tb + p / (2 * v), tb - p / (2 * v), [-0.9, 0.1], tb])
    d = M.dipole_law(L, K, dip + K.synth(np.array([[-0.9, 0.1]]), np.array([-0.3])), th, np.array([v, -v, -0.3, 0.5]))
    assert abs(d['dc'] - abs(p[0]) / v) < 1e-15
    g = K.synth(np.array([[1.3, 0.03]]), np.array([0.2]))
    d = M.dipole_law(L, K, dip + g, th[[0, 1, 3]], np.array([v, -v, 0.5]))
    assert abs(d['E_inf'] / (0.5 * np.mean(g ** 2)) - 1) < 1e-6 and abs(d['rms_inf'] ** 2 / 2 - d['E_inf']) < 1e-18
    assert abs(d['E'] - 0.5 * np.mean((dip + g - K.synth(th[[0, 1, 3]], np.array([v, -v, 0.5]))) ** 2)) < 1e-18


def test_h_statistics():
    """Referee round 61: H's statistics — the most kernels below −1% of the smallest amplitude in a block (a coefficient −0.005 is not
    one), the most negative coefficient, the end's formed kernels of the other sign, the smallest formed width (not an idle kernel's),
    and over all the end's kernels those narrower than the samples' spacing and those of negative width."""
    import exp_moving_width as M
    As = (1.0, 0.9)
    vis = lambda n: np.c_[np.linspace(-0.2, 0.2, n), np.full(n, 0.06)]
    blocks = [(vis(2), np.array([0.5, -0.005])), (vis(4), np.array([1.0, -0.02, -0.3, 0.9])), (vis(2), np.array([1.0, -0.4]))]
    m = dict(w=[1.0, -0.2, 0.9, 1e-4, -2e-4], theta=[[-0.08, 0.06], [0.0, 0.05], [0.08, 0.07], [0.3, 0.004], [0.4, -0.003]])
    h = M.h_stats(blocks, m, As)
    assert h == dict(neg_path=2, wmin_path=-0.4, neg_end=1, neg_end_all=1, smin=0.05, sub=2, negw=1, sub_formed=0, negw_formed=0,
                     inv_formed=0, wmax=1.0)
    cl = M.clusters(m, [-0.08, 0.08], [0.06, 0.06], As)
    assert [n for n, x in cl] == [2, 1] and abs(cl[0][1] - (0.06 - 0.2 * 0.05) / 0.06) < 1e-12 and abs(cl[1][1] - 0.9 * 0.07 / (0.9 * 0.06)) < 1e-12


def test_the_stalled_block_summary():
    """Referee round 61: G's stalled block — its closest pair of opposite signs (not a closer one of one sign), its largest |w|."""
    import exp_moving_width as M
    K = M.MW.BumpsCW(M.XS, M.S0)
    st = dict(c=[[0.0, 0.1], [0.001, 0.1], [0.3, 0.1]], w=np.array([50.0, 49.0, -98.0]))
    Rc = dict(stalled=st, polish_log=[[dict(status=2, nfev=10, accepted=True)]], nfev=7)
    f = np.zeros(M.XS.size)
    s = M.stalled_summary(Rc, f, K)
    assert s['wmax'] == 98.0 and abs(s['sep'] - 2.99) < 1e-12 and s['n'] == 3 and s['nfev'] == 7 and s['polish'][0]['nfev'] == 10
    assert M.stalled_summary(dict(Rc, stalled=None), f, K) is None


def test_the_loop_phases():
    """Referee round 61: `phases` on a synthetic loop whose dipole pair grows exactly as v³ = 3μ₀KΛ(τ − τ_d) until a third kernel forms
    at τ₃, then lands on the feature with v frozen: the slopes are 1, the landed v the law's (3μ₀KΛ(τ₃ − τ_d))^{1/3}; the collapse's
    exponents from the pair's distance and the residual over 1600 → 6400."""
    import exp_moving_width as M
    assert all(t in M.TAUS_G for t in M.TAUS_G_MD) and 30.0 in M.TAUS_G
    cs, ss, As = [0.0, 0.5], [0.1, 0.1], [0.6, -0.5]
    law = [None, dict(v=100.0, E=0.01 + 1e-6, E_inf=0.01, mu0=0.1)]
    K_, Lam, td, t3 = 1e4 * 1e-6, 50.0, 17.0, 30.0
    rate = 3 * 0.1 * K_ * Lam
    loops = []
    for tau in M.TAUS_G:
        v = (rate * (min(tau, t3) - td)) ** (1 / 3) if tau > td else 0.5
        d = 1e-3 * (tau / 400.0) ** -0.2 if tau >= t3 else 0.01
        W = 0.6 if tau >= t3 else 0.3
        wa, wb = W / 2 + v, W / 2 - v
        th = [[-wb * d * 0.1 / W, 0.1], [wa * d * 0.1 / W, 0.1]]
        ex = [[0.5, 0.1]] if tau >= t3 else []
        loops.append(dict(eps=1e-8, tau=tau, Lam=Lam, theta=th + ex, w=[wa, wb] + ([-0.5] if tau >= t3 else []), wmax=abs(wa),
                          rms=1e-6 * (tau / 400.0) ** -0.8))
    r = dict(cs=cs, ss=ss, As=As, entries=2, tau_d=td, law=law, loops=loops)
    ph = M.phases(r, 1e-8)
    assert ph['t3'] == 30.0 and len(ph['slopes']) >= 4 and max(abs(s - 1) for s in ph['slopes']) < 1e-9
    assert abs(ph['v400'] / ph['vstar'] - 1) < 1e-12 and abs(ph['rate'] - rate) < 1e-12
    assert all(ph['land'][t]['off'] < 1e-9 and abs(ph['land'][t]['dW']) < 1e-12 for t in M.TAUS_G_MD)
    assert abs(ph['collapse'][0] + 0.2) < 1e-9 and abs(ph['collapse'][1] + 0.8) < 1e-9
    md = '\n'.join(M.g_loop_lines([dict(r, dmin=0.4, q=9, loops=loops)]))
    assert (f"| 0.4 | 9 | 1e-08 | 17 | 30 | 1.000 [1.000, 1.000] | {ph['v400']:.4g} | {ph['vstar']:.4g} | yes, yes, yes | -0.200, -0.800 |"
            in md)


def test_a_pair_fits_any_feature_sharing_its_centre():
    """Referee round 61: a pair fits a feature if its net coefficient and centre of mass match any feature's — two features may share a
    centre (D at d = 0: the narrow feature split in two halves on the pedestal), and the nearest by centre alone is then ambiguous."""
    import exp_moving_width as M
    th = [[0.0, 0.06], [2e-4, 0.06]]
    assert M.pair_fits(th, [0.25, 0.25], 0, 1, [0.0, 0.0], [0.3, 0.06], [0.8, 0.5])
    assert not M.pair_fits(th, [0.25, 0.20], 0, 1, [0.0, 0.0], [0.3, 0.06], [0.8, 0.5])
    assert not M.pair_fits(th, [0.25, 0.25], 0, 1, [0.1, 0.1], [0.3, 0.06], [0.8, 0.5])


def test_the_oracle_fit_starts_at_the_truth():
    """Referee round 61 (E4): the oracle fit is the K-kernel least-squares fit of the data started at the truth — without noise it is
    the truth itself (error 0), with noise its error is below the noise's rms."""
    import exp_moving_width as M
    cs, ss, As = [-0.4, 0.3], [0.08, 0.15], [0.7, -0.5]
    K = MW.BumpsCW(M.XS, M.S0)
    f = K.synth(np.c_[cs, ss], np.array(As))
    assert M.oracle_error(f, np.zeros_like(f), cs, ss, As) < 1e-12
    xi = 1e-3 * np.random.default_rng(1).standard_normal(f.size)
    e = M.oracle_error(f, xi, cs, ss, As)
    assert 1e-5 < e < 1e-3


def test_pair_fits_and_pair_index_at_their_boundaries():
    """Referee round 68 (F2, F4, F5, I1): a pair fits a feature only with W of the feature's sign and size (not −A, not 0) and its centre
    of mass — not its midpoint — at the feature; `pair_index` pairs kernels only when both centre and width agree."""
    th = np.array([[0.0, 0.06], [0.02, 0.06]])
    w = np.array([1.5, -0.5])                                          # W = 1, centre of mass −0.01, midpoint +0.01
    assert MWE.pair_fits(th, w, 0, 1, [-0.01], [0.06], [1.0])
    assert not MWE.pair_fits(th, w, 0, 1, [0.01], [0.06], [1.0])      # the midpoint is not the centre of mass
    assert not MWE.pair_fits(th, -w, 0, 1, [-0.01], [0.06], [1.0])    # W = −1: the wrong sign
    assert not MWE.pair_fits(th, np.array([0.5, -0.5]), 0, 1, [0.0], [0.06], [0.0])   # W = 0 fits nothing
    assert MWE.pair_index(np.array([[0.1, 0.06], [0.1, 0.0601]]), [1.0, 1.0]) == [(0, 1)]
    assert MWE.pair_index(np.array([[0.1, 0.3], [0.1, 0.06]]), [1.0, 1.0]) == []      # one centre, two widths


def test_h_stats_counts_negatives_below_one_percent_of_the_smallest_amplitude():
    """Referee round 68 (H1, H7): a block's only negative at −0.5% of the smallest amplitude is not counted, one at −2% is; the
    threshold is the smallest amplitude's, not the largest's."""
    m = dict(w=[1.0, -0.005], theta=[[0.0, 0.06], [0.1, 0.06]])
    th = np.array(m['theta'])
    r = MWE.h_stats([(th, np.array([1.0, -0.005]))], m, (1.0, 0.6))
    assert r['neg_path'] == 0 and r['neg_end'] == 0 and r['wmin_path'] == -0.005
    r = MWE.h_stats([(th, np.array([1.0, -0.008]))], dict(m, w=[1.0, -0.008]), (1.0, 0.6))
    assert r['neg_path'] == 1 and r['neg_end'] == 1                    # −0.008 is beyond 1% of 0.6, not of 1.0


def test_stop_anatomy():
    """Referee round 68: `stop_anatomy` counts the active kernels (all: the dictionary exhausted), the accepted polishes and those stopped
    by their budget, the kernels of negative width and those invisible on the samples."""
    d = MWE.C0.size
    act = np.ones(d, bool)
    c = np.c_[np.linspace(-1, 1, d), np.full(d, 0.06)]
    c[3] = [0.775, -3e-4]                                               # negative width, invisible on the samples
    log = [[dict(status=0, accepted=True), dict(status=1, accepted=False)], [dict(status=2, accepted=True)]]
    r = MWE.stop_anatomy(dict(act=act, c=c, polish_log=log))
    assert r == dict(n_act=d, exhausted=True, pol_acc=2, pol_capped=1, neg_sigma=1, invisible=1)
    act[0] = False
    assert MWE.stop_anatomy(dict(act=act, c=c, polish_log=[]))['exhausted'] is False


def test_fine_phases_with_an_intermediate_phase():
    """Referee round 68: `fine_phases` finds τ₃ and τ₄ as formations, v before τ₃, the intermediate phase between them, the first
    landing after τ₃ with its p, and counts where p ≤ P_LAND agrees with the landing."""
    cs, ss, As = [0.0, 0.5], [0.06, 0.06], [1.0, 0.8]
    dip = lambda v, c: (np.array([[c - 0.01, 0.06], [c + 0.01, 0.06]]), np.array([1.0 + v, -v]))
    sts = [(10.0, *dip(0.5, 0.2)), (12.0, *dip(0.6, 0.2))]               # the dipole phase (n0 = 2), off the features
    th, w = dip(0.7, 0.2)
    sts.append((14.0, np.r_[th, [[0.5, 0.06]]], np.r_[w, 0.8]))          # τ₃: a third kernel forms
    sts.append((20.0, np.r_[th, [[0.5, 0.06]]], np.r_[w, 0.8]))          # the intermediate phase
    th2, w2 = np.array([[-0.0001, 0.06], [0.0003, 0.06]]), np.array([1.05, -0.05])
    sts.append((26.0, np.r_[th2, [[0.5, 0.06], [0.9, 0.06]]], np.r_[w2, 0.8, 0.01]))   # τ₄, landed
    r = MWE.fine_phases(sts, 2, 5.0, cs, ss, As)
    assert (r['t3'], r['t4'], r['v3']) == (14.0, 26.0, 1.1)                # v = |w_a − w_b|/2 at τ = 12
    assert r['inter']['n'] == 3 and r['inter']['off'] > 1 and r['t_land'] == 26.0 and abs(r['v_land'] - 0.55) < 1e-12
    p = 1.1 * 0.0004 / (1.0 * 0.06)
    assert abs(r['p_land'] - p) < 1e-12 and r['agree'][1] >= 1


def test_fine_phases_ignores_an_early_landing_and_an_exit():
    """Referee round 80 (FPf, FPg, FPh): a landing before τ₃ is not the landing; a state with fewer formed kernels after τ₃ (an exit) is
    not τ₄; the intermediate phase is read at the first state from (τ₃ + τ₄)/2, not at τ₃; the intermediate phase's v only at n₃
    formed."""
    cs, ss, As = [0.0, 0.5], [0.06, 0.06], [1.0, 0.8]
    off = lambda c: (np.array([[c - 0.01, 0.06], [c + 0.01, 0.06]]), np.array([1.5, -0.5]))
    landed = (np.array([[-0.0001, 0.06], [0.0003, 0.06]]), np.array([1.05, -0.05]))
    third = (np.array([[0.5, 0.06]]), np.array([0.8]))
    cat = lambda a, *b: (np.vstack([a[0]] + [x[0] for x in b]), np.concatenate([a[1]] + [x[1] for x in b]))
    sts = [(10.0, *off(0.2)), (12.0, *landed),                            # the dipole phase, landed early at τ = 12
           (14.0, *cat(off(0.2), third)),                                 # τ₃ (n₃ = 3), the pair at 0.2: 3σ* from feature 1
           (16.0, *off(0.2)),                                             # an exit: two formed
           (18.0, *cat(off(0.2), third)), (20.0, *cat(off(0.3), third)),   # the midpoint (14 + 26)/2: the pair at 0.3
           (26.0, *cat(landed, third, (np.array([[0.9, 0.06]]), np.array([0.01]))))]   # τ₄, landed
    r = MWE.fine_phases(sts, 2, 5.0, cs, ss, As)
    assert (r['t3'], r['t4'], r['t_land']) == (14.0, 26.0, 26.0)
    assert abs(r['inter']['off'] - abs(0.28 - 0.5) / 0.06) < 1e-9 and r['inter']['n'] == 3
    assert r['dip2'] == [(14.0, 1.0), (20.0, 1.0)]                    # every tenth of 14, 18, 20 (16 has two formed) and the last


RAW_MW = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'moving_width_raw.npz')


def test_t3_from_duals_reproduces_gt_and_the_loop():
    """Referee round 68: τ₃ from the sliding path alone — the duals advanced to the dipole block's entry τ_d, then moved by the block's
    residual at the end of its flow — reproduces GT's raw for target 9 and lies within the fine grid's step (0.1) of the loop's τ₃ at the
    smallest mass, for every G target (raw)."""
    raw = MWE.E.load(RAW_MW)
    gt = {x['q']: x['t3_pred'] for x in raw['GT']}
    for q in (9, 20):                                                   # referee round 80 (GT1): gt_job's every field against GT's raw
        g = next(r for r in raw['G'] if r['q'] == q)
        o = MWE.gt_job((q, np.array(g['cs']), np.array(g['ss']), np.array(g['As']), int(g['entries']), g['flow'][-1]['theta'],
                        g['flow'][-1]['w']))
        x = next(r for r in raw['GT'] if r['q'] == q)
        assert set(o) == set(x)
        for k_, v in o.items():
            assert (v is None and x[k_] is None) or (isinstance(v, bool) and v == x[k_]) or np.isclose(v, x[k_], rtol=1e-9, atol=0), k_
    assert next(r for r in raw['GT'] if r['q'] == 20)['exact'] is True
    n = 0
    for r in raw['G']:
        if r['q'] in gt and r.get('fine'):
            ph = r['fine'][min(r['fine'])]
            assert abs(gt[r['q']] - ph['t3']) < 0.1
            n += 1
    assert n == 6


def test_hl_job_on_a_cheap_case(monkeypatch):
    """Referee round 68 (R8): `hl_job` on H's separated target at one mass and time — the formed set is the loop's kernels above 1% of the
    smallest amplitude, with its count of negatives, smallest width, largest |w| over all kernels and the residual."""
    monkeypatch.setattr(MWE, 'EPS_HL', (1e-6,))
    monkeypatch.setattr(MWE, 'TAUS_HL', (60.0,))
    name, cs, ss, As = next(t for t in MWE.H_TARGETS if t[0] == 'σ* 0.06, A₂ 0.6')
    r = MWE.hl_job((name, cs, ss, As))
    L, K = MWE.law(), MWE.MW.BumpsCW(MWE.XS, MWE.S0)
    f = MWE.feature_sum(cs, ss, As)
    _, W, C, Lam, _ = MWE.RN.loop_nd(L, K, f, MWE.theta0(), 1e-6, np.r_[0.0, 60.0], rtol=1e-9, atol=1e-12)
    big = np.flatnonzero(np.abs(W[1]) > 1e-2 * 0.6)
    (q,) = r['loops']
    assert np.allclose(q['w'], W[1][big], rtol=1e-12, atol=0) and np.allclose(q['theta'], C[1][big], rtol=1e-12, atol=0)
    assert q['neg'] == int(np.sum(W[1][big] < 0)) and q['wmax'] == float(np.abs(W[1]).max())
    assert q['smin'] == float(np.abs(C[1][big][:, 1]).min())
    assert np.isclose(q['rms'], np.sqrt(np.mean((f - K.synth(C[1], W[1])) ** 2)), rtol=1e-12, atol=0)


def test_h_job_with_a_tolerance_and_h3_targets():
    """Referee round 68 (R7, R9, R10): `h_job` with a tolerance records the last entry's τ (the path's last event, not its first), the last
    blocks' residuals (the history's last six) and the clusters; H3's targets are equal amplitudes 1 at ±c for every separation."""
    name, cs, ss, As = next(t for t in MWE.H_TARGETS if t[0] == 'σ* 0.06, A₂ 0.6')
    r = MWE.h_job((name, cs, ss, As, 1e-13))
    Rm = MWE.RN.sliding_iss_nd(MWE.law(), MWE.MW.BumpsCW(MWE.XS, MWE.S0), MWE.feature_sum(cs, ss, As), MWE.theta0(),
                               np.ones(MWE.C0.size), max_events=40, exit_kkt=True, max_rhs=MWE.MAX_RHS_MOVING, tol=1e-13)
    ent = [e for e in Rm['events'] if e[1] is not None]
    assert len(ent) >= 2 and r['last_tau'] == float(Rm['events'][-1][0]) and r['last_tau'] > float(Rm['events'][0][0])
    assert r['rms_blocks'] == [float(h[3]) for h in Rm['hist'][-6:]]
    assert r['outcome'] == 'exact'
    for nm, c, s_, A in MWE.h3_targets():
        assert tuple(A) == (1.0, 1.0) and tuple(s_) == (0.06, 0.06) and c[0] == -c[1] and nm.endswith(f"±{c[1]:g}")


def test_oracle_error_against_a_direct_fit():
    """Referee round 68 (E2, E3): the oracle's error is the rms against f of the K-kernel least-squares fit of the data f + ξ (not of
    f + ½ξ, and not measured against the data), started at the truth — checked against a direct fit."""
    from scipy.optimize import least_squares
    cs, ss, As = (-0.3, 0.25), (0.05, 0.08), (1.0, -0.7)
    f = MWE.feature_sum(cs, ss, As)
    xi = 0.01 * np.random.default_rng(5).standard_normal(MWE.XS.size)
    K = MWE.MW.BumpsCW(MWE.XS, MWE.S0)
    sol = least_squares(lambda z: f + xi - K.synth(z[2:].reshape(2, 2), z[:2]), np.r_[As, np.c_[cs, ss].ravel()], method='lm',
                        xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    e = np.sqrt(np.mean((K.synth(sol.x[2:].reshape(2, 2), sol.x[:2]) - f) ** 2))
    assert np.isclose(MWE.oracle_error(f, xi, cs, ss, As), e, rtol=1e-10, atol=0)
    assert abs(e - np.sqrt(np.mean(xi ** 2))) > 1e-3


def test_c_table_counts_a_budget_stop_at_its_stalled_block():
    """Referee round 68 (D3bo-8): C's table counts a budget-stopped path at its stalled block — the residual G's `stalled_summary` records
    for the same target, the formed kernels the stalled block's (raw)."""
    raw = MWE.E.load(RAW_MW)
    g = {(r['dmin'], r['q']): r for r in raw['G']}
    n = 0
    for r in raw['C']:
        if r['moving'].get('stalled') is not None:
            k, rms = MWE.moving_state(r)
            st = g[(r['dmin'], r['q'])]['stalled']
            assert np.isclose(rms, st['rms'], rtol=1e-9, atol=0) and k == int(np.sum(np.abs(r['moving']['stalled']['w']) > 1e-2 * min(abs(a) for a in r['As'])))
            n += 1
    assert n == 6


def test_visibility_and_h_stats_count_only_visible_kernels():
    """Referee round 74 (D3bu-2): a kernel is visible on the samples if its largest φ there is at least VIS — a spike between two
    samples (σ = 3·10⁻⁴, 0.005 from the nearest) is not, one of largest φ ≈ 3.7·10⁻⁶ (σ = 10⁻³) is —; H's negative kernels along the
    path and at the end count the visible only (the end's formed ones of the other sign including the invisible in neg_end_all), the
    formed invisible ones are counted, and the end's sub-sample and negative-width kernels with how many are formed."""
    v = MWE.visibility([[0.005, 3e-4], [0.005, 1e-3], [0.0, 0.06], [0.2, 0.0]])
    assert v[0] < MWE.VIS < v[1] and abs(v[1] - np.exp(-0.005 ** 2 / (2e-6))) < 1e-18 and v[2] == 1.0 and v[3] == 0.0
    th = np.array([[-0.08, 0.06], [0.005, 3e-4], [0.08, 0.06], [0.3, -0.06]])
    w = np.array([1.2, -0.7, 1.1, -0.3])
    r = MWE.h_stats([(th, w)], dict(theta=th.tolist(), w=w.tolist()), (1.0, 1.0))
    assert (r['neg_path'], r['wmin_path'], r['neg_end'], r['neg_end_all'], r['inv_formed']) == (1, -0.3, 1, 2, 1)
    assert (r['sub'], r['sub_formed'], r['negw'], r['negw_formed']) == (1, 1, 1, 1) and r['smin'] == 0.06 and r['wmax'] == 1.2


def test_h_stats_does_not_count_an_unformed_invisible_kernel():
    """Referee round 80 (HS3): `h_stats`' formed invisible kernels are formed (|w| above 1% of the smallest amplitude) and invisible on
    the samples; an unformed invisible one (|w| 10⁻³, a spike between samples) is not counted."""
    m = dict(w=[1.0, 1e-3, 0.5], theta=[[0.0, 0.06], [0.0049, 1e-5], [0.2049, 1e-5]])
    st = MWE.h_stats([], m, [1.0])
    assert st['inv_formed'] == 1 and st['sub'] == 2 and st['sub_formed'] == 1


def test_path_stop_reasons():
    """Referee round 74: why a path ended — its tolerance, its budget, the dictionary exhausted, the events' cap, no entrant."""
    h = lambda rms: [(0.0, None, None, 1.0), (1.0, None, None, rms)]
    act = np.ones(21, bool)
    assert MWE.path_stop(dict(hist=h(1e-10), out_of_budget=True, act=act), 1e-9) == 'tolerance'
    assert MWE.path_stop(dict(hist=h(1e-6), out_of_budget=True, act=act), 1e-9) == 'budget'
    assert MWE.path_stop(dict(hist=h(1e-6), out_of_budget=False, act=act), 1e-9) == 'dictionary'
    act[3] = False
    assert MWE.path_stop(dict(hist=h(1e-6) * 20, out_of_budget=False, act=act), 1e-9) == 'events'
    assert MWE.path_stop(dict(hist=h(1e-6), out_of_budget=False, act=act), 1e-9) == 'no entrant'


def test_h_job_passes_its_polish_and_tolerance_and_clusters_unequal_amplitudes(monkeypatch):
    """Referee round 74 (HJ1, HJ2, HJ5, HJ6): `h_job` with (…, None, False) runs the path without the polish at the tolerance 1e-9, with
    a tolerance at that one with the polish; its stop and its clusters' masses Σw_i|σ_i|/(A_kσ_k) with unequal amplitudes."""
    seen = []
    th = np.array([[-0.08, 0.06], [0.081, 0.061], [0.5, 0.06]])
    w = np.array([1.0, 0.6, 0.0])

    def fake(L, K, f, c0, rho, **kw):
        seen.append(kw)
        act = np.array([True, True] + [False] * (c0.shape[0] - 2))
        c = np.array(c0, float); c[:2] = th[:2]
        ww = np.zeros(c0.shape[0]); ww[:2] = w[:2]
        return dict(act=act, c=c, w=ww, hist=[(0.0, c, np.zeros_like(ww), 1.0), (1.0, c, ww, 2e-10)], events=[(1.0, 0, 'enter', c[0])],
                    out_of_budget=False, nfev=10, polish_log=[], stalled=None)
    monkeypatch.setattr(MWE.RN, 'sliding_iss_nd', fake)
    name, cs, ss, As = 'x', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.6)
    r = MWE.h_job((name, cs, ss, As, None, False))
    assert seen[-1]['polish'] is False and seen[-1]['tol'] == 1e-9 and r['polish'] is False and r['stop'] == 'tolerance'
    r = MWE.h_job((name, cs, ss, As, 1e-13))
    assert seen[-1]['polish'] is True and seen[-1]['tol'] == 1e-13 and r['stop'] == 'no entrant'
    assert r['clusters'] == [(1, 1.0 * 0.06 / (1.0 * 0.06)), (1, 0.6 * 0.061 / (0.6 * 0.06))]


def test_stop_anatomy_counts_a_faint_kernel_as_visible():
    """Referee round 74 (SA4): a kernel of largest φ ≈ 10⁻⁶ on the samples is visible (VIS = 10⁻¹⁰), one between samples with
    σ = 3·10⁻⁴ is not."""
    c = np.array([[0.005, 1.2e-3], [0.005, 3e-4], [0.0, 0.06]])
    r = MWE.stop_anatomy(dict(act=np.array([True, True, True]), c=c, polish_log=[]))
    assert r['invisible'] == 1 and MWE.visibility(c[:1])[0] > 1e-7


def test_d3phi_and_the_third_order_k():
    """Referee rounds 74, 80: `d3phi` is the third derivative of φ(θ̄ + tp) in t (against differences), and `k_third` is the dipole's K: a
    pair at v, w = W/2 ± v, Δθ = p/v about the optimal dipole fit's θ̄ (`closure_fit`) has v²(E − E_inf) → K, the relative error
    falling like 1/v² — E − E_inf = K/v² + O(v⁻⁴), the pair's sum even in 1/v about the midpoint (D3ca-7)."""
    K = MW.BumpsCW(MWE.XS, MWE.S0)
    tb, p = np.array([0.05, 0.11]), np.array([0.013, -0.007])
    ph = lambda t: K.phi((tb + t * p)[None, :])[:, 0]
    h = 1e-3
    fd = (ph(2 * h) - 2 * ph(h) + 2 * ph(-h) - ph(-2 * h)) / (2 * h ** 3)
    assert np.abs(MWE.d3phi(tb, p) - fd).max() < 1e-4 * np.abs(fd).max()
    f = MWE.feature_sum([0.0, 0.03], [0.06, 0.05], [1.0, -0.6])
    cf = MWE.closure_fit(K, f, np.zeros((0, 2)), np.zeros(0), 0.4, np.array([0.015, 0.06]), np.array([0.015, 0.0]))
    k = MWE.k_third(cf['r'], cf['tb'], cf['p'])
    e_inf = np.mean(cf['r'] ** 2) / 2

    def est(v):
        d = cf['p'] / v
        g = K.synth(np.array([cf['tb'] + d / 2, cf['tb'] - d / 2]), np.array([cf['W'] / 2 + v, cf['W'] / 2 - v]))
        return v ** 2 * (np.mean((f - g) ** 2) / 2 - e_inf)
    e1, e2 = abs(est(30.0) / k - 1), abs(est(300.0) / k - 1)       # beyond v ~ 10³ the pair's cancellation rounds E − E_inf
    assert e2 < 0.05 * e1 and e2 < 1e-3 and k != 0 and np.sqrt(2 * e_inf) < 1e-3
    c30, c100 = (est(30.0) / k - 1) * 30.0 ** 2, (est(100.0) / k - 1) * 100.0 ** 2
    assert abs(c30 / c100 - 1) < 0.01 and c30 < 0                     # O(v⁻²) relative, not O(v⁻¹)


def test_closure_fit_returns_positive_widths():
    """Referee round 80 (D3ca-10): φ is even in σ, so a fit may end at a negative width; `closure_fit` returns |σ| — the other kernel
    started at σ = −0.06 is returned as the fit started at +0.06 is, at a width within 10⁻⁶ of 0.06, with the same residual."""
    K = MW.BumpsCW(MWE.XS, MWE.S0)
    f = MWE.feature_sum([0.0, 0.03, 0.5], [0.06, 0.05, 0.06], [1.0, -0.6, 0.7])
    fit = lambda s_: MWE.closure_fit(K, f, np.array([[0.5, s_]]), np.array([0.7]), 0.4, np.array([0.015, 0.06]), np.array([0.015, 0.0]))
    cn, cp = fit(-0.06), fit(0.06)
    assert np.array_equal(cn['th'], cp['th']) and np.array_equal(cn['r'], cp['r']) and abs(cn['th'][0, 1] - 0.06) < 1e-6


def test_first_wall_and_kappa_from_the_path():
    """Referee rounds 74, 80: `first_wall` — the first inactive dual to reach ±1 at the correlations v —; `kappa_from_path` on target 9
    (G's raw block): K from the third order is G's v²(E − E_inf), τ₃ and τ₄ within the fine grid's step of the loop's at the smallest
    mass, κ within 0.2% of the loop's κ′ (its v carried by the law to the path's τ₃: the grid's offset removed); K′ > 0, so no
    dissolution time; the three-kernel closure not exact."""
    j, t = MWE.first_wall(np.array([0.5, -0.2, 0.9]), np.array([False, True, False]), np.array([0.1, 1.0, -0.5]))
    assert j == 2 and abs(t - 3.8) < 1e-12
    assert MWE.first_wall(np.zeros(2), np.array([True, False]), np.array([1.0, 0.0])) == (None, float('inf'))
    raw = MWE.E.load(RAW_MW)
    g = next(r for r in raw['G'] if r['q'] == 9)
    K, L = MW.BumpsCW(MWE.XS, MWE.S0), MWE.law()
    f = MWE.feature_sum(g['cs'], g['ss'], g['As'])
    R = MWE.RN.sliding_iss_nd(L, K, f, MWE.theta0(), np.ones(MWE.C0.size), max_entries=int(g['entries']), exit_kkt=True)
    o = MWE.kappa_from_path(K, f, R, g['flow'][-1]['theta'], g['flow'][-1]['w'])
    lw = g['law'][-1]
    ph = g['fine'][min(g['fine'])]
    assert abs(o['K_third'] / (lw['v'] ** 2 * (lw['E'] - lw['E_inf'])) - 1) < 1e-4
    assert abs(o['t3_path'] - ph['t3']) < 0.1 and abs(o['t4_path'] - ph['t4']) < 0.2
    assert abs(o['kappa_path'] / MWE.kappa_loop(ph, g['tau_d'], o['t3_path'])[0] - 1) < 2e-3
    assert o['exact'] is False and o['t_diss'] is None and o['Kp'] > 0 and o['rms2'] > MWE.CLOSURE_EXACT
    assert abs(o['kappa_path'] ** 3 - 1 - o['Kp'] * (o['t4_path'] - o['t3_path']) / (o['K_third'] * (o['t3_path'] - g['tau_d']))) < 1e-9
    assert o['t3_path'] == MWE.t3_from_duals(K, f, R, g['flow'][-1]['theta'], g['flow'][-1]['w'])


def test_kappa_from_path_dissolution_and_an_exact_closure(monkeypatch):
    """Referee round 80 (D3ca-1, D3ca-6): `kappa_from_path` with its fits and walls stubbed — for K′ < 0 the dissolution time
    τ₃ + K(τ₃ − τ_d)/|K′| and κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3} while τ₄ comes first; no κ when τ₄ comes after it (the bracket
    negative); for K′ > 0 no dissolution; a closure fitting to below CLOSURE_EXACT is exact: K′ = 0, no τ₄, κ = 1."""
    class K_:
        def synth(self, th, w):
            return np.zeros(3)

        def corr(self, cen, r):
            return (np.zeros(4),)
    th, w = np.array([[0.0, 0.06], [0.01, 0.06], [0.5, 0.06]]), np.array([2.0, -1.5, 0.8])

    def run(k2, d4, r2=np.full(3, 1e-3)):
        ks = iter((2e-3, k2))
        walls = iter(((1, 10.0), (2, d4)))
        fits = iter((dict(r=np.full(3, 1e-2), th=th[2:], w=w[2:], W=0.5, tb=np.zeros(2), p=np.ones(2)),
                     dict(r=r2, th=th[1:], w=w[1:], W=0.5, tb=np.zeros(2), p=np.ones(2))))
        monkeypatch.setattr(MWE, 'k_third', lambda *a: next(ks))
        monkeypatch.setattr(MWE, 'first_wall', lambda *a: next(walls))
        monkeypatch.setattr(MWE, 'closure_fit', lambda *a, **k: next(fits))
        monkeypatch.setattr(MWE, 'path_duals_at', lambda *a: (np.zeros(4), np.array([True, False, False, False]), np.zeros((4, 2)), 10.0))
        return MWE.kappa_from_path(K_(), np.zeros(3), None, th, w)
    o = run(-4e-4, 30.0)                                                # τ₃ = 20, τ₄ = 50, dissolution at 20 + 2e-3·10/4e-4 = 70
    assert (o['t3_path'], o['t4_path']) == (20.0, 50.0) and abs(o['t_diss'] - 70.0) < 1e-12 and o['exact'] is False
    assert abs(o['kappa_path'] - 0.4 ** (1 / 3)) < 1e-12 and abs(o['rms2'] - 1e-3) < 1e-15
    o = run(-4e-4, 60.0)                                                # τ₄ = 80 after the dissolution: no κ
    assert o['t4_path'] == 80.0 and abs(o['t_diss'] - 70.0) < 1e-12 and o['kappa_path'] is None
    o = run(1e-4, 30.0)
    assert o['t_diss'] is None and abs(o['kappa_path'] - 1.15 ** (1 / 3)) < 1e-12
    o = run(1e-30, 30.0, r2=np.full(3, 1e-17))                          # an exact closure: its K′ and τ₄ are roundoff
    assert o['exact'] is True and o['Kp'] == 0.0 and o['t4_path'] is None and o['kappa_path'] == 1.0 and o['t_diss'] is None


def test_part_gt_takes_the_status_zero_entries(monkeypatch):
    """Referee round 74 (PG1): GT runs on G's entries whose block flow ended (status 0), not the others."""
    got = []
    monkeypatch.setattr(MWE, 'gt_job', lambda a: got.append(a[0]) or dict(q=a[0]))
    G = [dict(q=q, status=st, cs=[0.0], ss=[0.1], As=[1.0], entries=2, flow=[dict(theta=[[0.0, 0.1]], w=[1.0])]) for q, st in ((1, 0), (2, -1), (3, 0))]
    assert [x['q'] for x in MWE.part_GT(G)] == [1, 3] and got == [1, 3]


def test_fine_phases_counts_a_sigma_dipole_with_its_width_difference():
    """Referee round 74 (FP6, FP8, FP11, FP13): p takes |Δθ| — a landed σ-dipole (Δc = 0, Δσ ≠ 0) with 0.28 < p < 1 disagrees with the
    landing —, a state with 0.28 < p < 1 off the features agrees, states with p ≥ 1 are not counted; the dipole phase's v is kept at
    every tenth time with the last one's time; the intermediate phase's v and the sign of its (W − A)/|A|."""
    cs, ss, As = [0.0, 0.5], [0.06, 0.06], [1.0, 0.8]
    sts = []
    for k in range(25):                                                 # the dipole phase: p ≥ 1, off the features
        v = (0.5 * (k + 1)) ** (1 / 3)
        sts.append((6.0 + k, np.array([[0.2 - 0.05, 0.06], [0.2 + 0.05, 0.06]]), np.array([1.0 + v, -v])))
    third = [[0.5, 0.06]]
    for k in range(30):                                                 # τ₃ = 31 … 60: three formed, the pair at W = 0.5 off the features
        v = (sts[-1][2][0] - sts[-1][2][1]) / 2 + 0.01 * k
        sts.append((31.0 + k, np.r_[np.array([[0.3, 0.06], [0.31, 0.06]]), third], np.array([0.25 + v, 0.25 - v, 0.8])))
    land = np.array([[0.0, 0.06], [0.0001, 0.06], [0.5, 0.06], [0.9, 0.06]])
    sts.append((61.0, land, np.array([1.05, -0.05, 0.8, 0.02])))          # τ₄, landed (p small)
    sts.append((62.0, np.array([[0.0, 0.065], [0.0, 0.055], [0.5, 0.06], [0.9, 0.06]]), np.array([2.0, -1.0, 0.8, 0.02])))  # σ-dipole
    sts.append((63.0, np.array([[0.1, 0.06], [0.106, 0.06], [0.5, 0.06], [0.9, 0.06]]), np.array([2.0, -1.0, 0.8, 0.02])))  # off, p 0.3
    r = MWE.fine_phases(sts, 2, 5.0, cs, ss, As)
    assert (r['t3'], r['t4'], r['t_v3']) == (31.0, 61.0, 30.0) and [t for t, v in r['dip']] == [6.0, 16.0, 26.0]
    assert [t for t, v in r['dip2']][:3] == [31.0, 41.0, 51.0] and r['dip2'][-1][0] == 60.0
    assert r['inter']['n'] == 3 and r['inter']['dW'] < 0 and r['t_land'] == 61.0
    p62 = MWE.pair_p(sts[-2][1], sts[-2][2], 0, 1)
    assert abs(p62 - 3 * 0.01 / (1.0 * 0.06)) < 1e-12 and 0.28 < p62 < 1
    assert abs(MWE.pair_p(sts[-1][1], sts[-1][2], 0, 1) - 0.3) < 1e-9
    # counted: τ = 61 (landed, p small: agrees), 62 (landed, p 0.5: disagrees), 63 (not landed, p 0.3: agrees); p ≥ 1 elsewhere
    assert r['agree'] == (2, 3) and r['agree4'] == (1, 2)


def test_gl_job_and_its_constant(monkeypatch):
    """Referee rounds 74, 80 (GJ1, GJ3, GJ4, GJd): `gl_job` reads the loop's states (`loop_nd` monkeypatched) — the dipole pair's Δc,
    w_a, w_b, mean width (of both, at unequal widths) and nearest amplitude, E of all kernels and of the formed ones, Λ, Gm, μ₀ —;
    `gl_constant` on rows on the law
    Δc⁻⁶ = Δc_L⁻⁶ + 48aGm(W²/(w_a² + w_b²))Λτ is 1 (C = 0 at W = A) and its C is Gm W²(W − A)²(1/|w_a| + 1/|w_b|)/(μ₀Δc²(w_a² + w_b²))."""
    cs, ss, As = [0.0], [0.06], [1.0]
    taus = (400.0, 800.0)
    monkeypatch.setattr(MWE, 'TAUS_GL', taus)
    dcs = (2e-3, 1e-3)
    Wt = [np.zeros(3)] + [np.array([3.0, -2.0, 1e-4]) for _ in taus]
    Ct = [np.zeros((3, 2))] + [np.array([[-0.4 * d, 0.05], [0.6 * d, 0.07], [0.5, 0.06]]) for d in dcs]
    monkeypatch.setattr(MWE.RN, 'loop_nd', lambda *a, **k: (None, Wt, Ct, 77.0, None))
    g = MWE.gl_job((5, cs, ss, As))
    K, L = MW.BumpsCW(MWE.XS, MWE.S0), MWE.law()
    f = MWE.feature_sum(cs, ss, As)
    x = g['rows'][1]
    assert g['Lam'] == 77.0 and g['Gm'] == L.G * L.m and g['mu0'] == 4 * L.G * L.c / L.A and g['q'] == 5
    assert abs(x['dc'] - 1e-3) < 1e-15 and (x['wa'], x['wb'], x['A']) == (3.0, -2.0, 1.0) and abs(x['sbar'] - 0.06) < 1e-15
    assert np.isclose(x['E'], np.mean((f - K.synth(Ct[2], Wt[2])) ** 2) / 2) and np.isclose(x['E_formed'], np.mean((f - K.synth(Ct[2][:2], Wt[2][:2])) ** 2) / 2)
    Gm, mu0, Lam, a = 12.0, 0.1, 50.0, 3e6
    wa, wb = 3.0, -2.0
    R = 48 * a * Gm * 1.0 / (wa ** 2 + wb ** 2) * Lam
    row = lambda t, dW=0.0: dict(tau=t, dc=(1e14 + R * t) ** (-1 / 6), wa=wa, wb=wb, sbar=0.06, A=1.0, dW=dW,
                                 E=a * (1e14 + R * t) ** (-8 / 6))
    p, r0, r1, c = MWE.gl_constant(row(100.0), row(200.0), Gm, mu0, Lam)
    assert abs(r0 - 1) < 1e-9 and abs(r1 - 1) < 1e-9 and c == 0.0
    x1 = row(200.0, dW=-0.01)
    p, r0, r1, c = MWE.gl_constant(row(100.0, dW=-0.01), x1, Gm, mu0, Lam)
    x0 = row(100.0, dW=-0.01)
    cc = lambda x: Gm * 1.0 * 0.01 ** 2 * (1 / 3 + 1 / 2) / (mu0 * x['dc'] ** 2 * 13.0)
    rate = lambda x: 48 * x['E'] / x['dc'] ** 8 * Gm * 1.0 / 13.0 * Lam
    meas = (x1['dc'] ** -6 - x0['dc'] ** -6) / 100.0
    assert abs(c - cc(x1)) < 1e-9 * c and abs(r0 - 1) < 1e-9
    assert abs(r1 - meas / (0.5 * (rate(x0) / (1 + cc(x0)) + rate(x1) / (1 + cc(x1))))) < 1e-9 * r1 and r1 > r0
    assert abs(p - 5.0 * x1['dc'] / 0.06) < 1e-12


def test_a2_one_kernel_from_large_coefficients(monkeypatch):
    """Referee rounds 74, 80 (D3bu-5, D3ca-3): from (c*, σ*) with w = 3A (E(start) four times ½⟨f*²⟩, outside (i)'s condition) or 10A
    the one-kernel block collapses its width onto the samples at once — J_N's σ → 0 limit, E/½⟨f*²⟩ = 0.906 at σ*/σ₀ = 1 whatever w —,
    stays there (3A to s ≈ 3·10⁵, 10A to 10¹⁰) and escapes to the exact fit (10A by 3·10¹⁰); at σ*/σ₀ = 0.25 the end is σ/σ* = 1 (A2b);
    the flow's tolerances (A2d); A2's table and the spike's summary."""
    r3, r10, rq = MWE.a2_job((1.0, 3.0)), MWE.a2_job((1.0, 10.0)), MWE.a2_job((0.25, 2.5))
    for r in (r3, r10, rq):
        e = r['rows'][-1]
        assert abs(e['w'] - 1) < 1e-8 and abs(e['sig'] - 1) < 1e-8 and e['E'] < 1e-20 and r['status'] == 0
    assert r3['E0'] == 4.0 and r10['E0'] == 81.0
    sp3, sp10 = MWE.a2_spike(r3['rows']), MWE.a2_spike(r10['rows'])
    assert sp3[0] == 1.0 and np.isclose(sp3[1], 10 ** 5.5, rtol=1e-12) and np.isclose(sp3[3], 1e6, rtol=1e-12)
    assert sp10[0] == 1.0 and np.isclose(sp10[1], 1e10, rtol=1e-12) and np.isclose(sp10[3], 10 ** 10.5, rtol=1e-12)
    assert MWE.a2_spike(rq['rows']) is None
    assert abs(sp10[2][0] - 0.906) < 1e-3 and abs(sp10[2][1] - sp3[2][1]) < 1e-4    # the plateau's E is J_N's limit, for every w
    md = '\n'.join(MWE.a2_lines([r3, r10]))
    e3 = r3['rows'][-1]
    assert (f"| 1 | 3 | 4 | {e3['w']:.4f}, {e3['sig']:.4g}, {e3['dc']:+.1e}, {e3['E']:.2e} (s = 1e+13) | "
            f"{min(x['sig'] for x in r3['rows']):.3g} | 1e+00 … 3e+05, {sp3[2][0]:.3f}–{sp3[2][1]:.3f} | 1e+06 | 0 |") in md
    assert MWE.a2_lines(None) == [] and MWE.a2_lines([]) == []
    import scipy.integrate
    seen = {}
    real = scipy.integrate.solve_ivp
    monkeypatch.setattr(scipy.integrate, 'solve_ivp', lambda *a, **k: seen.update(k) or real(*a, **k))
    MWE.a2_job((1.0, 2.5))
    assert (seen['rtol'], seen['atol'], seen['method']) == (1e-10, 1e-14, 'LSODA')


def test_loop_row_counts_visible_kernels_and_hll_job(monkeypatch):
    """Referee round 74 (R8, D3bu-2): HL's row counts the visible formed kernels of the other sign and the invisible formed ones;
    `hll_job` keeps the loop's kernels above 1% of the smallest amplitude (one at 0.5% is not formed)."""
    x = dict(tau=1e5, theta=[[-0.08, 0.06], [0.005, 3e-4], [0.08, 0.06], [2.0, 0.06]], w=[1.2, -0.7, -0.3, 0.5], wmax=1.2, rms=2e-6)
    assert MWE.loop_row('t', 1e-8, x) == "| t | 1e-08 | 100000 | 4 | 1 (-0.300) | 1 | 0.060 | 1 | 1.2 | 2.0e-06 |"
    monkeypatch.setattr(MWE, 'TAUS_HLL', (10.0,))
    Wt = [np.zeros(3), np.array([1.0, -0.005, 0.02])]
    Ct = [np.zeros((3, 2)), np.array([[-0.08, 0.06], [0.0, 0.06], [0.08, 0.06]])]
    got = []
    monkeypatch.setattr(MWE.RN, 'loop_nd', lambda *a, **k: got.append(a[4]) or (None, Wt, Ct, 1.0, None))
    r = MWE.hll_job(('σ* 0.06, equal', 1e-6))
    assert r['loops'][0]['w'] == [1.0, 0.02] and r['loops'][0]['wmax'] == 1.0 and r['eps'] == 1e-6 and got == [1e-6]   # HLLb
    sp = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06]], [1.0, 0.02], (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0))
    assert {k_: r['loops'][0][k_] for k_ in sp} == sp                   # referee round 80: hll_job's split is `fit_split`'s
    monkeypatch.setattr(MWE, 'EPS_HL', (1e-7,))                         # HLa: hl_job's formed set at 1%, not 0.1%
    monkeypatch.setattr(MWE, 'TAUS_HL', (10.0,))
    got.clear()
    q = MWE.hl_job(('t', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.6)))
    assert q['loops'][0]['w'] == [1.0, 0.02] and got == [1e-7] and q['loops'][0]['neg'] == 0


def test_fit_split_tells_a_null_cluster_from_a_cancellation():
    """Referee round 80 (D3ca-2): `fit_split` — per feature the nearest formed kernel of its sign (each used once) and its largest
    relative error; the others' count, Σ|w|, sum on the samples over max f*, widths and largest |c| —; `split_cells` and `end_split`."""
    cs, ss, As = (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)
    null = MWE.fit_split([[-0.08, 0.06], [0.08, 0.06], [0.0, 1.5], [1.9, 0.3]], [1.0, 1.0, 0.4, -0.4], cs, ss, As)
    K = MW.BumpsCW(MWE.XS, MWE.S0)
    foot = np.abs(K.synth(np.array([[0.0, 1.5], [1.9, 0.3]]), np.array([0.4, -0.4]))).max() / np.abs(MWE.feature_sum(cs, ss, As)).max()
    assert null['main_err'] == [0.0, 0.0] and null['n_rest'] == 2 and null['rest_w'] == 0.8 and abs(null['rest_foot'] - foot) < 1e-15
    assert null['rest_sig'] == [0.3, 1.5] and null['rest_c'] == 1.9
    canc = MWE.fit_split([[-0.08, 0.06], [0.0, 0.05], [0.0, 0.04]], [1.0, 3.0, -2.4], cs, ss, As)
    assert canc['main_err'] == [0.0, 2.0] and canc['n_rest'] == 1 and canc['rest_foot'] > 1   # feature 2's: the +3 kernel at 0
    assert MWE.fit_split([[0.0, 0.06]], [-1.0], (0.0,), (0.06,), (1.0,))['main_err'] == [None]
    once = MWE.fit_split([[0.0, 0.06], [0.5, 0.06]], [1.0, 1.0], cs, ss, As)     # each kernel main to one feature only (FS1)
    assert once['main_err'][0] == abs(0.0 + 0.08) / 0.06 and abs(once['main_err'][1] - abs(0.5 - 0.08) / 0.06) < 1e-12
    assert MWE.split_cells(null) == f"0.0e+00, 0.0e+00 | 2, 0.800, {foot:.1e}, 0.300–1.500, 1.90"
    r = dict(moving=dict(theta=[[-0.08, 0.06], [0.08, 0.06], [0.5, 0.06]], w=[1.0, 1.0, 0.005]), cs=cs, ss=ss, As=As)
    assert MWE.end_split(r)['n_rest'] == 0                              # the unformed kernel (0.5%) left out


def test_g_fine_lines_and_gl_lines_on_synthetic_results():
    """Referee rounds 74, 80: the fine-grid table evaluates the dipole law at the time of v(τ₃⁻) and the intermediate phase's slopes over
    3μ₀K′Λ, κ and κ′ (v carried by the law to the path's τ₃), † for another sequence; the path's rows (GT: K, K′/K with the closure's
    rms, τ₃, τ₄, the dissolution time, κ, the landing's v by the law, v after τ₄; an exact closure); κ′ against the mass, the grid's
    offsets and the agreement's counts (all, ε ≤ 10⁻⁷, past τ₄); GL's constant table, its ranges and the local slope."""
    Lam, mu0, K_, Kp = 1e3, 0.1, 2e-3, -4e-4
    rate = 3 * mu0 * K_ * Lam
    t_d, t3, t4 = 10.0, 20.0, 30.0
    dip = [(t, (rate * (t - t_d)) ** (1 / 3)) for t in (11.0, 15.0, 19.9)]
    v3 = dip[-1][1]
    dip2 = [(t, (v3 ** 3 + 3 * mu0 * Kp * Lam * (t - t3)) ** (1 / 3)) for t in (20.0, 25.0, 29.9)]
    ph = dict(t3=t3, t4=t4, v3=v3, t_v3=19.9, dip=dip, dip2=dip2, inter=dict(n=3, off=2.5, dW=-0.4), t_land=31.0, v_land=0.97 * v3,
              p_land=0.1, agree=(5, 6), agree4=(2, 3))
    th = [[0.0, 0.06], [0.001, 0.06], [0.5, 0.06]]
    loops = [dict(eps=1e-10, tau=6400.0, Lam=Lam, theta=th, w=[0.97 * v3 + 0.5, -0.97 * v3 + 0.5, 0.8], wmax=1.0, rms=1e-9)]
    lam6 = MWE.RT.transit_scale(MWE.law(), 1e-6)
    dip2x = [(t, (v3 ** 3 - 3 * 3 * mu0 * Kp * lam6 * (t - t3)) ** (1 / 3)) for t in (20.0, 25.0, 29.9)]   # −3× the law: another sequence
    G = [dict(q=4, dmin=0.4, cs=[0.0, 0.5], ss=[0.06, 0.06], As=[1.0, 0.8], entries=2, tau_d=t_d, status=0, loops=loops,
              law=[dict(v=10.0, E=K_ / 100.0 + 0.5, E_inf=0.5, mu0=mu0)], fine={1e-6: dict(ph, agree=(1, 4), dip2=dip2x), 1e-10: ph})]
    G.append(dict(G[0], q=5, fine={1e-10: ph}))
    K3 = K_ * (1 + 1e-5)
    td_ = 19.95 + K3 * (19.95 - t_d) / abs(Kp)
    GT = [dict(q=4, t3_pred=19.95, K_third=K3, Kp=Kp, t3_path=19.95, t4_path=29.95, kappa_path=0.96, t_diss=td_, rms2=1e-3, exact=False),
          dict(q=5, t3_pred=19.95, K_third=K3, Kp=0.0, t3_path=19.95, t4_path=None, kappa_path=1.0, t_diss=None, rms2=4e-17, exact=True)]
    md = '\n'.join(MWE.g_fine_lines(G, GT))
    kp_ = 0.97 / ((19.95 - t_d) / (19.9 - t_d)) ** (1 / 3)
    assert f"| 4 | 1e-10 | 10 | 20 | 30 | 1.0000 | 1.0000 | 3, 2.5, -0.400, 1.0000 | 31, 0.970, {kp_:.3f}, 0.100 |" in md
    assert f"| 4 | 1e-06 | 10 | 20 | 30 | {(1.0 * Lam / lam6):.4f} |" in md
    assert f"| 3, 2.5, -0.400, -3.0000 | 31, 0.970, {kp_:.3f}†, 0.100 |" in md
    vl = (3 * mu0 * Lam * (K3 * (19.95 - t_d) + Kp * 10.0)) ** (1 / 3)
    vv = (0.97 * v3 + 0.5 - (-0.97 * v3 + 0.5)) / 2
    assert (f"| 4 | {K3:.5e} ({K_:.5e}) | 19.95 (20) | {Kp / K3:+.4f} (1.0e-03) | 29.95 (30) | {td_:.4g} | 0.9600 ({kp_:.4f}) | "
            f"{vl:.4g} ({0.97 * v3:.4g}) | {vv / (0.97 * v3) - 1:+.1e} |") in md
    vl5 = (3 * mu0 * Lam * K3 * (19.95 - t_d)) ** (1 / 3)
    assert f"| 5 | {K3:.5e} ({K_:.5e}) | 19.95 (20) | 0, exact (4.0e-17) | — (exact) (30) | — | 1.0000 ({kp_:.4f}) | {vl5:.4g} (" in md
    assert f"κ′ against the mass (ε = 1e-06, 1e-10; † another sequence): 4: {kp_:.3f}†, {kp_:.3f}; 5: {kp_:.3f}." in md
    off = ((19.95 - t_d) / (19.9 - t_d)) ** (1 / 3) - 1
    assert f"that κ′ removes, at ε = 1e-10: 4: {off:+.2%}, 5: {off:+.2%}." in md
    assert "agree in 11 of 16 (10 of 12 at ε ≤ 10⁻⁷; past τ₄ 6 of 9)." in md
    rows = [dict(tau=t, d=1e-3 * t ** (-1 / 6), dW=-0.01 * t ** (-1 / 3), rms=1e-6 * t ** (-2 / 3), E=1e-12, E_formed=0.9e-12,
                 dc=1e-3 * t ** (-1 / 6), wa=3.0, wb=-2.0, sbar=0.06, A=1.0) for t in (400.0, 800.0, 1600.0)]
    md = '\n'.join(MWE.gl_lines([dict(q=9, rows=rows, Lam=10.0, Gm=12.0, mu0=0.1)]))
    cons = [MWE.gl_constant(a, b, 12.0, 0.1, 10.0) for a, b in zip(rows[:-1], rows[1:])]
    assert "### GL. The collapse (referee rounds 68, 74): target 9" in md
    assert f"| 800 | {cons[0][0]:.4f} | {cons[0][1]:.4f} | {cons[0][2]:.4f} | {cons[0][3]:.4f} |" in md
    assert (f"Over the 2 intervals the ratio is {min(c[1] for c in cons):.4f}–{max(c[1] for c in cons):.4f}, with C "
            f"{min(c[2] for c in cons):.4f}–{max(c[2] for c in cons):.4f}; the formed kernels' energy over all kernels' 0.900000–0.900000.") in md
    assert "The remainder with C against p" not in md                   # p's range below 0.05: no fit
    # referee rounds 74, 80: E = a₈Δc⁸(1 + kp) along the run — the local slope (d ln E/d ln Δc − 8)/p of the formed kernels' E, the
    # force factor's remainder, the mean width's share p²/36
    k8, sb, wa, wb = 0.8, 0.06, 3.0, -2.0
    pr = lambda dc: abs(wa - wb) * dc / (abs(wa + wb) * sb)
    rows = [dict(tau=t, d=dc, dW=-0.01, rms=1e-6, E=1e3 * dc ** 8 * (1 + k8 * pr(dc)) + 1e-23, E_formed=1e3 * dc ** 8 * (1 + k8 * pr(dc)),
                 dc=dc, wa=wa, wb=wb, sbar=sb, A=1.0) for t, dc in ((400.0, 3e-3), (800.0, 2e-3), (1600.0, 1e-3))]   # E: an idle floor
    md = '\n'.join(MWE.gl_lines([dict(q=9, rows=rows, Lam=10.0, Gm=12.0, mu0=0.1)]))
    sl = [np.log(b['E_formed'] / a['E_formed']) / np.log(b['dc'] / a['dc']) for a, b in zip(rows[:-1], rows[1:])]
    ks = [(x - 8) / (0.5 * (pr(a['dc']) + pr(b['dc']))) for x, a, b in zip(sl, rows[:-1], rows[1:])]
    assert all(abs(k_ / k8 - 1) < 0.25 for k_ in ks)                   # the local slope sees E's next order, kp/(1 + kp)
    cons = [MWE.gl_constant(a, b, 12.0, 0.1, 10.0) for a, b in zip(rows[:-1], rows[1:])]
    rem = [c[2] / (x / 8) - 1 for c, x in zip(cons, sl)]
    assert f"(p their mean) {min(ks):.3f}–{max(ks):.3f}, median {np.median(ks):.3f}" in md and f"p²/36, at most {pr(3e-3) ** 2 / 36:.5f}." in md
    assert f"the ratio with C is 1 {min(rem):+.1e} … {max(rem):+.1e}" in md
