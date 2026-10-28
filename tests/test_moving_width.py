"""Proposition D‴.59: the transit with moving widths (`src/core/moving_width.py`)."""
import os
import sys

import numpy as np

from src.core import moving_width as MW
from src.core import rich_transit_nd as RN

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))


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
    ws = [np.array([0.5, -0.005]), np.array([1.0, -0.02, -0.3, 0.9]), np.array([1.0, -0.4])]
    m = dict(w=[1.0, -0.2, 0.9, 1e-4, -2e-4], theta=[[-0.08, 0.06], [0.0, 0.05], [0.08, 0.07], [0.3, 0.004], [0.4, -0.003]])
    h = M.h_stats(ws, m, As)
    assert h == dict(neg_path=2, wmin_path=-0.4, neg_end=1, smin=0.05, sub=2, negw=1, wmax=1.0)
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
