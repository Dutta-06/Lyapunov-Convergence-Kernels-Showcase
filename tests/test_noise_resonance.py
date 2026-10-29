"""Tests for experiments/exp_noise_resonance.py (Proposition D‴.60): the oscillator law's transfer function, the curve's mobilities, the
averaged loop's first order against the linearised loop, its tone averages against direct time averages, the loop's chunked period
means, the thresholds' bookkeeping and the report."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_noise_resonance as NR

CT, DZ = NR.CT, NR.DZ


def test_transfer_reduces_to_the_high_pass_without_adaptation():
    """At K = 0 and f̂′ = 0 the law's Y = 1/(iω + k_e) and the measured error's iωY = iω/(iω + k_e)."""
    om = np.array([3.0, 17.0])
    o = dict(f1=np.zeros(2), K=np.zeros(2))
    Y, Ee = NR.transfer(o, om, 10.0)
    assert np.allclose(Y, 1 / (1j * om[None, :] + 10.0))
    assert np.allclose(Ee, 1j * om[None, :] / (1j * om[None, :] + 10.0))
    o = dict(f1=np.array([2.0]), K=np.array([200.0]))
    Y, _ = NR.transfer(o, np.array([5.0]), 4.0)
    assert np.isclose(Y[0, 0], (5j - 2.0) / (200.0 - 25.0 + 20j))


def test_curve_inverts_and_its_derivatives():
    """alpha_of inverts w(α); curve's μ agrees with mobilities' and dμ/dw, dv/dw with finite differences along the curve."""
    L = NR.law()
    w = np.array([0.9, -0.6, 0.05])
    al = NR.alpha_of(L, NR.EPS, w, np.sign(w) * 3.0)
    for metric in NR.METRICS:
        ww, mu, dmu, v, dv = NR.curve(L, NR.EPS, al, metric)
        assert np.allclose(ww, w, rtol=1e-12, atol=1e-14)
        w0, mu0, v0 = NR.mobilities(L, NR.EPS, al, metric)
        assert np.allclose(mu, mu0) and np.allclose(v, v0)
        h = 1e-6
        wp, mup, _, vp, _ = NR.curve(L, NR.EPS, al + h, metric)
        wm, mum, _, vm, _ = NR.curve(L, NR.EPS, al - h, metric)
        assert np.allclose(dmu, (mup - mum) / (wp - wm), rtol=1e-6)
        assert np.allclose(dv, (vp - vm) / (wp - wm), rtol=1e-6)
    _, _, _, v, dv = NR.curve(L, NR.EPS, al, 'paradigm')
    assert np.allclose(v, L.G * L.m * w) and np.allclose(dv, L.G * L.m)              # v = Gβa = Gmw under the paradigm's metric


def test_jitter_filter():
    """The moving average over one period removes a linear trend and keeps a sinusoid of that period."""
    nwin = 101                                                                          # the sinusoid's period in samples (odd: centred)
    t = np.arange(4000) / nwin
    z = 0.3 + 0.02 * t + 0.1 * np.sin(2 * np.pi * t)
    assert np.allclose(NR.jitter(z, nwin), (0.1 * np.sin(2 * np.pi * t))[nwin:-nwin], atol=1e-12)
    assert np.abs(NR.jitter(0.3 + 0.02 * t, nwin)).max() < 1e-12


def _fit():
    L = NR.law()
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    return L, al0, c0, idx


def test_the_fitted_state_is_a_rest_of_the_averaged_loop_without_noise():
    L, al0, c0, idx = _fit()
    nz = dict(w=np.array([20.0]), ph=np.array([0.0]), amp=np.array([0.0]))
    for quasi in (False, True):
        rhs, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, quasi=quasi)
        y = np.zeros(n)
        y[1:4], y[4:7] = CT.ASTAR, CT.CSTAR
        for t in (0.0, 3.0, 11.0):
            assert np.abs(rhs(t, y)).max() < 1e-12


def test_the_averaged_loops_first_order_is_the_linearised_loop():
    """At a tiny amplitude the envelopes' Re(Θc e^{iωt}) is the linearised loop's δc (`variational`) from the same start."""
    L, al0, c0, idx = _fit()
    from scipy.integrate import solve_ivp
    om, amp, ph0 = 19.0, 1e-7, 0.4
    nz = dict(w=np.array([om]), ph=np.array([ph0]), amp=np.array([amp]))
    ts = np.linspace(9.0, 10.5, 301)                                                    # x_m passes the 0.61 kernel near t = 3.3, 12.4
    for metric in NR.METRICS:
        rhs, n = NR.envelope_rhs(L, NR.EPS, metric, nz, CT.KE, CT.CSTAR, CT.ASTAR)
        y0 = np.zeros(n)
        y0[1:4], y0[4:7] = CT.ASTAR, CT.CSTAR
        sol = solve_ivp(rhs, (9.0, 10.5), y0, t_eval=ts, method='DOP853', rtol=1e-10, atol=1e-16)
        Z = sol.y[7:7 + 7] + 1j * sol.y[14:21]                                          # one tone: (E, Θw×3, Θc×3)
        dc_env = np.real(Z[4:7] * np.exp(1j * om * ts)[None, :]).T
        _, dc_var = NR.variational(L, NR.EPS, al0[idx], c0[idx], metric, ts, om, amp, ph0, CT.KE)
        assert np.abs(dc_env - dc_var).max() < 1e-7 * np.abs(dc_var).max()


def test_the_tone_averages_are_time_averages():
    """For hand-built envelopes, the blur and the correlation terms equal their direct averages over one tone period."""
    L, al0, c0, idx = _fit()
    om = 15.0
    nz = dict(w=np.array([om]), ph=np.array([0.0]), amp=np.array([1e-3]))
    rhs_all, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR)
    rhs_no_blur, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl', 'cross', 'corr'))
    rhs_no_corr, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl', 'cross', 'blur'))
    rng = np.random.default_rng(3)
    Z = (rng.standard_normal(7) + 1j * rng.standard_normal(7)) * 1e-3
    y = np.zeros(n)
    y[1:4], y[4:7] = CT.ASTAR, CT.CSTAR
    y[7:14], y[14:21] = Z.real, Z.imag
    t = 2.0                                                                              # x_m = sin(0.4) ≈ 0.39
    blur = rhs_all(t, y)[0] - rhs_no_blur(t, y)[0]
    xm = np.sin(CT.OM * t)
    u = (xm - CT.CSTAR) / CT.SIG
    ph = np.exp(-u ** 2 / 2)
    d1, d2 = ph * u / CT.SIG, ph * (u ** 2 - 1) / CT.SIG ** 2
    s = np.linspace(0, 2 * np.pi, 4001)[:-1]
    rot = np.exp(1j * s)
    tw = np.real(Z[1:4, None] * rot[None, :])
    tc = np.real(Z[4:7, None] * rot[None, :])
    direct = -0.5 * np.mean((2 * d1[:, None] * tw * tc + CT.ASTAR[:, None] * d2[:, None] * tc ** 2).sum(0))
    assert np.isclose(blur, direct, rtol=1e-10, atol=1e-18)
    corr = rhs_all(t, y)[1:7] - rhs_no_corr(t, y)[1:7]
    al = NR.alpha_of(L, NR.EPS, CT.ASTAR, np.sign(CT.ASTAR) * 3.0)
    _, mu, dmu, v, dv = NR.curve(L, NR.EPS, al, 'paradigm')
    D = np.real((Z[0] + nz['amp'][0]) * rot)
    cw = np.mean((dmu[:, None] * ph[:, None] * tw + mu[:, None] * d1[:, None] * tc) * D[None, :], axis=1)
    cc = np.mean((dv[:, None] * d1[:, None] * tw + v[:, None] * d2[:, None] * tc) * D[None, :], axis=1)
    assert np.allclose(corr, np.r_[cw, cc], rtol=1e-10, atol=1e-18)


def test_quasi_static_envelopes_are_the_oscillator_law():
    """quasi=True: the measured error's envelope is A's law iωY·N and the correlation averages vanish."""
    L, al0, c0, idx = _fit()
    om = np.array([21.37])
    nz = dict(w=om, ph=np.array([0.0]), amp=np.array([1e-3]))
    rhs_q, nq = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, quasi=True)
    rhs_nc, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl', 'cross', 'blur'), quasi=True)
    assert nq == 7
    y = np.r_[0.0, CT.ASTAR, CT.CSTAR]
    for t in (1.0, 2.6, 5.0):
        assert np.allclose(rhs_q(t, y)[1:], rhs_nc(t, y)[1:], atol=1e-15)


def test_sim_means_do_not_depend_on_the_chunks():
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(0.005, dict(freqs=(21.37,)))
    W1, C1 = NR.sim_means(L, 'paradigm', nz, CT.KE, 2, chunk=1, k=400)
    W2, C2 = NR.sim_means(L, 'paradigm', nz, CT.KE, 2, chunk=2, k=400)
    assert W1.shape == (2, 3) and np.allclose(W1, W2, atol=1e-9) and np.allclose(C1, C2, atol=1e-9)


def test_left_periods_and_thresholds():
    dev = np.array([[0.1, 0.0, 0.2], [0.31, 0.0, 0.29], [0.5, 0.1, 0.4]])
    assert NR.left_periods(dev, 0.3) == [2, -1, 3]
    rs = [dict(s=0.02, left=[5, -1, -1]), dict(s=0.01, left=[-1, -1, -1]), dict(s=0.03, left=[4, -1, 7])]
    assert NR.thresholds(rs) == [0.02, None, 0.03]


def test_report_on_synthetic_results():
    A = [dict(metric=mt, om=om, sqrtK=15.0 if mt == 'paradigm' else 4.8, v=[1.0, 2.0, 3.0], mu=[0.1, 0.2, 0.3],
              sim_c=[1.0, 2.0, 3.0], var_c=[1.0, 2.0, 2.0], qs_c=[1.1, 2.0, 2.0] if om >= 19 else [1.5, 2.0, 2.0],
              sim_w=[1.0, 1.0, 1.0], var_w=[1.0, 1.0, 1.0], qs_w=[1.0, 1.0, 1.0]) for mt in NR.METRICS for om in (13.9, 19.1)]
    mk = len(NR.MARKS_B)
    b = lambda mt, fr, s, terms, quasi, f: dict(metric=mt, freqs=fr, s=s, terms=list(terms), quasi=quasi, status=0,
                                                 dw_sim=[[1e-3, -2e-3, 4e-3]] * mk, dc_sim=[[1e-4, 0.0, 2e-4]] * mk,
                                                 dw_env=[[f * 1e-3, -2e-3, 4e-3]] * mk, dc_env=[[1e-4, 0.0, 2e-4]] * mk)
    mt, fr, s = NR.B2_CASE
    B = [b(mt, list(fr), s, NR.TERMS, False, 1.0), b('wfr', None, s, NR.TERMS, False, 1.0)]
    B2 = [b(mt, list(fr), s, ('ctl', 'cross', 'blur'), False, 0.5), b(mt, list(fr), s, NR.TERMS, True, 0.25)]
    C = [dict(m=12.0, ke=10.0, s=x, left=[3 if x >= 0.02 else -1, -1, 5 if x >= 0.03 else -1], sqrtK=15.76) for x in (0.01, 0.02, 0.03)]
    md = NR.report(A, B, B2, C)
    assert '[1.0000, 1.5000]' in md                                                     # the loop over the linearisation, centres
    assert '[1.000, 1.500] at ω ≥ 10, [1.000, 1.100] at ω ≥ 19' in md
    assert "D‴.49's 16 tones" in md and 'one tone 21.37' in md
    assert f"| without {NR.TERM_NAMES['corr']} | +5.00e-04, -2.00e-03, +4.00e-03 | 0.500, 1.000, 1.000 |" in md
    assert '| quasi-static envelopes | +2.50e-04, -2.00e-03, +4.00e-03 | 0.250, 1.000, 1.000 |' in md
    assert '| 12 | 10 | 15.76 | 0.02 (1/1) / 0.02 | — / — | 0.03 (1/1) / 0.03 | 3 (0) | 3 [3, 5] |' in md


def test_the_leaving_levels_over_seeds():
    """Referee round 56: per kernel the least level at which a seed left (how many there), the level from which every seed left at
    every level above, and '*' where a seed left at one level and stayed at a higher one."""
    row = lambda x, sd, left: dict(m=12.0, ke=10.0, s=x, opts=dict(seed=sd), left=left, sqrtK=15.76)
    rs = [row(0.01, 7, [-1, -1, -1]), row(0.01, 1, [40, -1, -1]), row(0.02, 7, [30, 50, -1]), row(0.02, 1, [20, -1, -1]),
          row(0.03, 7, [25, -1, -1]), row(0.03, 1, [10, -1, 60])]
    assert NR.c_cells(rs) == ['0.01 (1/2) / 0.02', '0.02 (1/2) / — *', '0.03 (1/2) / —']


def test_the_control_averages_are_time_averages():
    """'ctl' and 'cross' equal the direct averages over one tone period of ½f″e₁² − ½f̂″(e₁ + n)² and −∇_θf̂′·θ₁(e₁ + n), at a mean
    error ē ≠ 0 (f and f̂ at x_m + ē)."""
    L, al0, c0, idx = _fit()
    om = 15.0
    nz = dict(w=np.array([om]), ph=np.array([0.7]), amp=np.array([1e-3]))
    mk = lambda terms: NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=terms)[0]
    rhs_all, rhs_no_ctl, rhs_no_cross = mk(NR.TERMS), mk(('cross', 'blur', 'corr')), mk(('ctl', 'blur', 'corr'))
    rng = np.random.default_rng(5)
    Z = (rng.standard_normal(7) + 1j * rng.standard_normal(7)) * 1e-3
    y = np.zeros(7 + 14)
    y[0], y[1:4], y[4:7] = 0.01, CT.ASTAR * 1.01, CT.CSTAR + 0.003
    y[7:14], y[14:21] = Z.real, Z.imag
    t = 2.0
    xs = np.sin(CT.OM * t) + y[0]
    rot = np.exp(1j * np.linspace(0, 2 * np.pi, 4001)[:-1])
    e1 = np.real(Z[0] * rot)
    n = np.real(nz['amp'][0] * np.exp(1j * nz['ph'][0]) * rot)
    u = (xs - CT.CSTAR) / CT.SIG
    g = np.exp(-u ** 2 / 2)
    f2 = (CT.ASTAR * g * (u ** 2 - 1) / CT.SIG ** 2).sum()
    w, c = y[1:4], y[4:7]
    v_ = (xs - c) / CT.SIG
    ph = np.exp(-v_ ** 2 / 2)
    d1, d2 = ph * v_ / CT.SIG, ph * (v_ ** 2 - 1) / CT.SIG ** 2
    fh2 = (w * d2).sum()
    ctl = np.mean(0.5 * f2 * e1 ** 2 - 0.5 * fh2 * (e1 + n) ** 2)
    assert np.isclose(rhs_all(t, y)[0] - rhs_no_ctl(t, y)[0], ctl, rtol=1e-9, atol=1e-18)
    tw, tc = np.real(Z[1:4, None] * rot[None, :]), np.real(Z[4:7, None] * rot[None, :])
    cross = np.mean((d1[:, None] * tw + w[:, None] * d2[:, None] * tc).sum(0) * (e1 + n))      # −∇_θf̂′·θ₁ = Σ(∂_cφ θw + w∂²_cφ θc)
    assert np.isclose(rhs_all(t, y)[0] - rhs_no_cross(t, y)[0], cross, rtol=1e-9, atol=1e-18)


def test_the_envelopes_are_driven_at_the_regressors_point():
    """With ē ≠ 0 the envelopes' drive uses the regressor at x_m, the fit's gradient at x_m + ē."""
    L, al0, c0, idx = _fit()
    om = 15.0
    N = 1e-3 * np.exp(0.2j)
    nz = dict(w=np.array([om]), ph=np.array([0.2]), amp=np.array([1e-3]))
    rhs, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR)
    y = np.zeros(n)
    y[0], y[1:4], y[4:7] = 0.03, CT.ASTAR, CT.CSTAR
    t = 3.0
    dZ = rhs(t, y)[7:14] + 1j * rhs(t, y)[14:21]
    xm = np.sin(CT.OM * t)
    u = (xm - CT.CSTAR) / CT.SIG
    ph = np.exp(-u ** 2 / 2)
    w, mu, v = NR.mobilities(L, NR.EPS, NR.alpha_of(L, NR.EPS, CT.ASTAR, np.sign(CT.ASTAR) * 3.0), 'paradigm')
    assert np.allclose(dZ[1:4], mu * ph * N, rtol=1e-10)
    assert np.allclose(dZ[4:7], v * ph * u / CT.SIG * N, rtol=1e-10)
    us = (xm + 0.03 - CT.CSTAR) / CT.SIG
    fh1 = -(CT.ASTAR * np.exp(-us ** 2 / 2) * us / CT.SIG).sum()
    assert np.isclose(dZ[0], -fh1 * N - CT.KE * N, rtol=1e-10)


def test_quasi_static_blur_is_the_oscillator_laws():
    """quasi=True at the fit: the blur average is −½|Y|²p·Σ_i[2μ_iv_iφ_i(∂_cφ_i)² + v_i²w_i(∂_cφ_i)²∂²_cφ_i] with A's Y."""
    L, al0, c0, idx = _fit()
    om = np.array([21.37])
    nz = dict(w=om, ph=np.array([0.0]), amp=np.array([1e-3]))
    rq, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('blur',), quasi=True)
    r0, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=(), quasi=True)
    y = np.r_[0.0, CT.ASTAR, CT.CSTAR]
    for t in (2.2, 2.6, 3.4):
        o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', np.array([np.sin(CT.OM * t)]))
        Y = NR.transfer(o, om, CT.KE)[0][0, 0]
        blur = (2 * o['mu'] * o['v'] * o['ph'][0] * o['d1'][0] ** 2 + o['v'] ** 2 * o['w'] * o['d1'][0] ** 2 * o['d2'][0]).sum()
        assert np.isclose(rq(t, y)[0] - r0(t, y)[0], -0.5 * abs(Y) ** 2 * (1e-3) ** 2 / 2 * blur, rtol=1e-9, atol=1e-18)


def test_sim_means_at_the_noiseless_fit_are_the_fit():
    L, al0, c0, idx = _fit()
    W, C = NR.sim_means(L, 'paradigm', None, CT.KE, 2, chunk=1, k=200)
    assert np.allclose(W, CT.ASTAR[None, :], atol=1e-10) and np.allclose(C, CT.CSTAR[None, :], atol=1e-10)


def test_quasi_static_control_average_has_the_laws_phase():
    """quasi=True at the fit: 'ctl' is ¼f″|D − N|² − ¼f̂″|D|² with A's D = iωY·N — the phase of (iω − f̂′) matters here, not in the blur."""
    L, al0, c0, idx = _fit()
    om = np.array([21.37])
    N = 1e-3
    nz = dict(w=om, ph=np.array([0.0]), amp=np.array([N]))
    rq, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=('ctl',), quasi=True)
    r0, _ = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=(), quasi=True)
    y = np.r_[0.0, CT.ASTAR, CT.CSTAR]
    for t in (2.2, 2.6, 3.4):
        o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', np.array([np.sin(CT.OM * t)]))
        D = NR.transfer(o, om, CT.KE)[1][0, 0] * N
        f2 = o['f2'][0]
        assert abs(o['f1'][0]) > 0.5                                                      # f̂′ ≠ 0 there: the sign is tested
        assert np.isclose(rq(t, y)[0] - r0(t, y)[0], 0.25 * f2 * abs(D - N) ** 2 - 0.25 * f2 * abs(D) ** 2, rtol=1e-9, atol=1e-18)


def test_the_second_orders_first_order_is_the_linearised_loop():
    """`second_order`'s first order is the linearised loop (`variational`, the measured error n at the start): their period means agree."""
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(19.0,)))
    W1, C1, W2, C2, stt = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', nz, CT.KE, 1, k=400)
    per = 2 * np.pi / CT.OM
    ts = per * np.arange(401) / 400
    dw, dc = NR.variational(L, NR.EPS, al0[idx], c0[idx], 'paradigm', ts, 19.0, float(nz['amp'][0]), float(nz['ph'][0]), CT.KE)
    mw, mc = np.trapezoid(dw, axis=0) / 400, np.trapezoid(dc, axis=0) / 400
    assert stt == 0 and W1.shape == (1, 3)
    assert np.abs(W1[0] - mw).max() < 1e-7 * np.abs(dw).max() and np.abs(C1[0] - mc).max() < 1e-7 * np.abs(dc).max()


def test_the_loops_even_part_is_the_second_order():
    """Under s·n and −s·n the loop's period means' half-sum less the fit is s²⟨w₂⟩ + O(s⁴) (`second_order`), its half-difference
    s⟨w₁⟩ + O(s³): at s = 0.002, one tone at 21.37, two periods, the even part within 1%; the odd part's remainder grows by ≈ 2³ from s
    to 2s (the first kernel's cubic remainder is 13% of s⟨w₁⟩ at s = 0.002)."""
    L, al0, c0, idx = _fit()
    W1, C1, W2, C2, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', DZ.noise_of(1.0, dict(freqs=(21.37,))), CT.KE, 2)
    rem = []
    for s in (0.002, 0.004):
        nz = DZ.noise_of(s, dict(freqs=(21.37,)))
        Wp, Cp = NR.sim_means(L, 'paradigm', nz, CT.KE, 2)
        Wm, Cm = NR.sim_means(L, 'paradigm', dict(nz, ph=nz['ph'] + np.pi), CT.KE, 2)
        if s == 0.002:
            even_w, even_c = (Wp + Wm) / 2 - CT.ASTAR, (Cp + Cm) / 2 - CT.CSTAR
            assert np.abs(even_w - s ** 2 * W2).max() < 0.01 * np.abs(s ** 2 * W2).max()
            assert np.abs(even_c - s ** 2 * C2).max() < 0.01 * np.abs(s ** 2 * C2).max()
        rem.append(np.abs((Wp - Wm) / 2 - s * W1).max())
    assert 6.0 < rem[1] / rem[0] < 10.0


def test_the_b3_report_on_synthetic_results():
    """B3's table: the even part over s²⟨w₂⟩ per s, the averaged loop (B) over it, the odd part and the ranges' line."""
    mk = len(NR.MARKS_B)
    w2, c2, w1 = [2.0, -1.0, 4.0], [1.0, 1.0, -1.0], [1e-3, 0.0, 2e-3]
    rows = [dict(s=s, even_w=[[s * s * 2.0 * 1.01, -s * s, s * s * 4.0]] * mk, odd_w=[[1e-6, 0.0, 2e-6]] * mk,
                 even_c=[[s * s, s * s, -s * s]] * mk, odd_c=[[0.0, 0.0, 0.0]] * mk) for s in NR.S_B3]
    B3 = [dict(metric='paradigm', freqs=[21.37], status=0, w1=[w1] * mk, c1=[[0.0] * 3] * mk, w2=[w2] * mk, c2=[c2] * mk, rows=rows)]
    s = NR.S_B3[-1]
    B = [dict(metric='paradigm', freqs=[21.37], s=s, terms=list(NR.TERMS), quasi=False, status=0, dw_sim=[[0.0] * 3] * mk,
              dc_sim=[[0.0] * 3] * mk, dw_env=[[s * s * 2.0 * 1.05, -s * s, s * s * 4.0]] * mk, dc_env=[[0.0] * 3] * mk)]
    md = '\n'.join(NR.b3_lines(B3, B))
    assert '## B3. The exact second order' in md and '| paradigm | one tone 21.37 | +2, -1, +4 |' in md
    assert '1.0100, 1.0000, 1.0000; 1.0100, 1.0000, 1.0000; 1.0100, 1.0000, 1.0000 | 1.050, 1.000, 1.000 |' in md
    assert f"[1.0000, 1.0100] at s = {s:g}" in md and 'paradigm, one tone 21.37: [1.0000, 1.0500]' in md
    assert NR.b3_lines(None) == [] and '## B3' in NR.report(B3=B3)


def test_the_envelopes_couple_through_the_mean_error():
    """Referee round 56: at ē ≠ 0 the envelopes' own terms (μ′φΘw + μ∂_cφΘc)ē and (v′∂_cφΘw + v∂²_cφΘc)ē are the adaptation law
    linearised at the mean state — central differences of μ(w)φ(x_m; c) and v(w)∂_cφ(x_m; c) along (Θw, Θc) —, the regressor at x_m."""
    L, al0, c0, idx = _fit()
    om, e = 15.0, 0.05
    nz = dict(w=np.array([om]), ph=np.array([0.0]), amp=np.array([0.0]))             # N = 0 and E = 0 below: D = 0
    rhs, n = NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR)
    rng = np.random.default_rng(11)
    Z = np.r_[0.0, rng.standard_normal(6)] + 1j * np.r_[0.0, rng.standard_normal(6)]
    w0, c0_ = CT.ASTAR * 1.02, CT.CSTAR + 0.004
    t = 2.0
    out = rhs(t, np.r_[e, w0, c0_, Z.real, Z.imag])
    dZ = out[7:14] + 1j * out[14:21]
    xm = np.sin(CT.OM * t)

    def law_terms(w, c):
        _, mu, v = NR.mobilities(L, NR.EPS, NR.alpha_of(L, NR.EPS, w, np.sign(w) * 3.0), 'paradigm')
        u = (xm - c) / CT.SIG
        ph = np.exp(-u ** 2 / 2)
        return mu * ph, v * ph * u / CT.SIG
    h = 1e-6
    for part in (np.real, np.imag):
        tw, tc = part(Z[1:4]), part(Z[4:7])
        (gp, kp), (gm, km) = law_terms(w0 + h * tw, c0_ + h * tc), law_terms(w0 - h * tw, c0_ - h * tc)
        assert np.allclose(part(dZ[1:4] + 1j * om * Z[1:4]), e * (gp - gm) / (2 * h), rtol=1e-6, atol=1e-12)
        assert np.allclose(part(dZ[4:7] + 1j * om * Z[4:7]), e * (kp - km) / (2 * h), rtol=1e-6, atol=1e-12)


def test_the_blur_and_the_correlation_at_a_mean_error():
    """Referee round 56: at ē ≠ 0 the blur is the direct tone average of −½θ₁ᵀ∇²_θf̂θ₁ at x_m + ē (where the control evaluates f̂) and
    the regressors' correlation that of ∂_θ(Γg)·θ₁(e₁ + n) with the regressor at x_m."""
    L, al0, c0, idx = _fit()
    om, e = 15.0, 0.04
    nz = dict(w=np.array([om]), ph=np.array([0.3]), amp=np.array([1e-3]))
    mk = lambda terms: NR.envelope_rhs(L, NR.EPS, 'paradigm', nz, CT.KE, CT.CSTAR, CT.ASTAR, terms=terms)[0]
    rhs_all, rhs_no_blur, rhs_no_corr = mk(NR.TERMS), mk(('ctl', 'cross', 'corr')), mk(('ctl', 'cross', 'blur'))
    rng = np.random.default_rng(4)
    Z = (rng.standard_normal(7) + 1j * rng.standard_normal(7)) * 1e-3
    w, c = CT.ASTAR * 0.99, CT.CSTAR - 0.002
    y = np.r_[e, w, c, Z.real, Z.imag]
    t = 2.0
    rot = np.exp(1j * np.linspace(0, 2 * np.pi, 4001)[:-1])
    tw, tc = np.real(Z[1:4, None] * rot[None, :]), np.real(Z[4:7, None] * rot[None, :])

    def kern(x):
        u = (x - c) / CT.SIG
        ph = np.exp(-u ** 2 / 2)
        return ph, ph * u / CT.SIG, ph * (u ** 2 - 1) / CT.SIG ** 2
    xm = np.sin(CT.OM * t)
    _, d1s, d2s = kern(xm + e)
    blur = -0.5 * np.mean((2 * d1s[:, None] * tw * tc + w[:, None] * d2s[:, None] * tc ** 2).sum(0))
    assert np.isclose(rhs_all(t, y)[0] - rhs_no_blur(t, y)[0], blur, rtol=1e-10, atol=1e-18)
    phr, d1r, d2r = kern(xm)
    _, mu, dmu, v, dv = NR.curve(L, NR.EPS, NR.alpha_of(L, NR.EPS, w, np.sign(w) * 3.0), 'paradigm')
    D = np.real((Z[0] + nz['amp'][0] * np.exp(1j * nz['ph'][0])) * rot)
    cw = np.mean((dmu[:, None] * phr[:, None] * tw + mu[:, None] * d1r[:, None] * tc) * D[None, :], axis=1)
    cc = np.mean((dv[:, None] * d1r[:, None] * tw + v[:, None] * d2r[:, None] * tc) * D[None, :], axis=1)
    assert np.allclose(rhs_all(t, y)[1:7] - rhs_no_corr(t, y)[1:7], np.r_[cw, cc], rtol=1e-10, atol=1e-18)
    _, d1x, _ = kern(xm + e)                                                            # the regressor point matters here
    assert np.abs(d1x - d1r).max() > 0.1 * np.abs(d1r).max()


def test_b2_drops_each_average_once_and_adds_the_quasi_static_envelopes():
    jobs = NR.b2_jobs()
    assert [NR.variant_label(dict(terms=list(j[3]), quasi=j[4])) for j in jobs] == \
        [f"without {NR.TERM_NAMES[t]}" for t in NR.TERMS] + ['quasi-static envelopes']
    assert all(j[:3] == NR.B2_CASE for j in jobs) and jobs[-1][3] == NR.TERMS


def test_b3_job_on_a_short_horizon(monkeypatch):
    """Referee round 56: `b3_job`'s bookkeeping on two periods — the half-sum of the loop under ±s·n less the fit over s²⟨w₂⟩ within 1%
    at s = 0.002, the half-difference near s⟨w₁⟩, the second order's means at the marks."""
    monkeypatch.setattr(NR, 'NP_B', 2)
    monkeypatch.setattr(NR, 'MARKS_B', (1, 2))
    monkeypatch.setattr(NR, 'S_B3', (0.002, 0.004))
    r = NR.b3_job(('paradigm', (21.37,)))
    w2, w1 = np.array(r['w2']), np.array(r['w1'])
    assert w2.shape == (2, 3) and r['status'] == 0 and [q['s'] for q in r['rows']] == [0.002, 0.004]
    q0 = r['rows'][0]
    assert np.abs(np.array(q0['even_w']) - 0.002 ** 2 * w2).max() < 0.01 * np.abs(0.002 ** 2 * w2).max()
    assert np.abs(np.array(q0['odd_w']) - 0.002 * w1).max() < 0.2 * np.abs(0.002 * w1).max()
    L, al0, c0, idx = _fit()
    W1, C1, W2, C2, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', DZ.noise_of(1.0, dict(freqs=(21.37,))), CT.KE, 2)
    assert np.allclose(w2, W2) and np.allclose(r['c2'], C2)


def test_the_c_report_uses_each_laws_own_fitted_state():
    """Referee round 56: the per-kernel √K in C's caption — each law's own fitted state (√K ∝ m^{1/2} when the centres' term
    dominates: m = 3 at half m = 12's), the first kernel's equal to the orbit's √K max."""
    C = [dict(m=m, ke=10.0, s=0.01, left=[-1, -1, -1], sqrtK=0.0) for m in NR.M_C]
    md = NR.report(C=C)
    cap = md.split('at its largest: ')[1].split(':\n')[0]
    vals = {float(p.split(':')[0].split('=')[1]): [float(x) for x in p.split(':')[1].split(',')] for p in cap.split('; ')}
    assert np.allclose(np.array(vals[3.0]) * 2, vals[12.0], rtol=1e-2)
    for m in NR.M_C:
        L = NR.law(m)
        al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
        o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', np.sin(np.linspace(0, 2 * np.pi, 2001)))
        assert abs(vals[m][0] - np.sqrt(o['K'].max())) < 0.02


def test_the_monodromy_and_the_floquet_means_are_the_linearised_loops():
    """Referee round 63: the Floquet recursion's period means over three periods equal the period means of `second_order`'s first order
    from the same start (in the first period to the trapezoid's error on the start's fast transient, later to 10⁻⁷), and the unforced
    multipliers lie inside the unit circle with the three slow coefficient modes first."""
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.4,)))
    a, ph = float(nz['amp'][0]), float(nz['ph'][0])
    M, P, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, 21.4, ph)
    Z = NR.floquet_means(M, P, a * G, a * H, 21.4, 3)
    W1, C1, _, _, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', nz, CT.KE, 3)
    for X, Y in ((Z[:, 1:4], W1), (Z[:, 4:7], C1)):
        assert np.abs(X[1:] - Y[1:]).max() < 1e-7 * np.abs(Y).max() and np.abs(X[0] - Y[0]).max() < 1e-3 * np.abs(Y).max()
    rho = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
    assert 0.99 < rho[2] <= rho[0] < 0.995 and rho[3] < 0.3


def test_f_job_and_its_report(monkeypatch):
    """F's job at m = 12 (two B4 tones, one commensurate tone, 20 periods): the multipliers and time constants, the Floquet sums against
    a synthetic B4 whose odd part is the sums plus a known error, the saturation at the commensurate tone equal to the recursion's
    limit."""
    monkeypatch.setattr(NR, 'FREQS_B4', ((21.4,), None, (21.3,)))
    monkeypatch.setattr(NR, 'OM_COMM_F', (21.4,))
    monkeypatch.setattr(NR, 'NP_B4', 20)
    r = NR.f_job(('paradigm', 12.0))
    assert len(r['rho']) == 7 and 120 < min(r['tau'][:3]) and max(r['tau'][:3]) < 180 and set(r['mode'][:3]) == {'w1', 'w2', 'w3'}
    fz = np.array(r['sums']['21.4']['w'])
    assert fz.shape == (20, 3)
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.4,)))
    a = float(nz['amp'][0])
    M, P, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, 21.4, float(nz['ph'][0]))
    far = NR.floquet_means(M, P, a * G, a * H, 21.4, 4000)[-1]
    assert np.allclose(r['comm']['21.4']['sat'], far, rtol=1e-3, atol=1e-6)
    B4 = [dict(freqs=[21.4], s=s, odd_w_all=(fz * (1 + 0.01 * (s / 0.00125) ** 2)).tolist()) for s in (0.00125, 0.0025)]
    md = '\n'.join(NR.f_lines([r], B4))
    assert '21.4: 1.00e-02 (s = 0.00125), 4.00e-02 (s = 0.0025), ×4.00' in md and '| 21.4 | 107 |' in md


def test_b4_job_on_two_periods(monkeypatch):
    """Referee round 63: `b4_job`'s bookkeeping on two periods — the even part per s² near ⟨w₂⟩ (`second_order`), the odd part per s
    near ⟨w₁⟩, the full arrays' marks equal to the marked rows, the centre's per-period odd part."""
    monkeypatch.setattr(NR, 'NP_B4', 2)
    monkeypatch.setattr(NR, 'MARKS_B4', (1, 2))
    r = NR.b4_job(((21.37,), 0.002))
    L, al0, c0, idx = _fit()
    W1, C1, W2, C2, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', DZ.noise_of(1.0, dict(freqs=(21.37,))), CT.KE, 2)
    assert np.abs(np.array(r['even_w']) - W2).max() < 0.02 * np.abs(W2).max()
    assert np.abs(np.array(r['odd_w']) - W1).max() < 0.25 * np.abs(W1).max()                 # 21.37's odd part: its O(s³) is large
    assert np.abs(np.array(r['even_c']) - C2).max() < 0.1 * np.abs(C2).max()                  # the centres' even part (round 63)
    assert np.allclose(r['even_w_all'], r['even_w']) and np.allclose(r['odd_w_all'], r['odd_w']) and len(r['odd_c_per']) == 2


def test_c_noise_and_c_job_pass_their_options(monkeypatch):
    """Referee round 63: `c_noise` and `c_job` hand the seed, the band, the tones, the phases, the metric and k_e to the loop."""
    o = NR.c_noise(0.02, dict(freqs=(18.9, 22.5), ph=(0.1, 0.2)))
    assert np.allclose(o['w'], [18.9, 22.5]) and np.allclose(o['ph'], [0.1, 0.2]) and np.allclose(o['amp'], 0.02)
    a, b = NR.c_noise(0.02, dict(seed=1)), NR.c_noise(0.02, dict(seed=2))
    assert len(a['w']) == 16 and not np.allclose(a['w'], b['w']) and 10 <= a['w'].min() and a['w'].max() <= 30
    c = NR.c_noise(0.02, dict(band=(25.0, 75.0), seed=3))
    assert 25 <= c['w'].min() and c['w'].max() <= 75
    d1, d2 = NR.c_noise(0.04, dict(freqs=(13.1,), seed=1)), NR.c_noise(0.04, dict(freqs=(13.1,), seed=2))
    assert d1['ph'][0] != d2['ph'][0] and np.isclose(d1['amp'][0], 0.04 * np.sqrt(2))
    seen = {}

    def fake(L, al0, c0, astar, ts, noise, metric, ke=None, x_init=None):
        seen.update(m=L.m, noise=noise, metric=metric, ke=ke, n=len(ts))
        n = len(ts)
        return dict(c=np.tile(c0, (n, 1)), w=np.zeros((n, c0.size)), al=np.zeros((n, c0.size)), x=np.zeros(n))
    monkeypatch.setattr(NR, 'simulate', fake)
    r = NR.c_job((6.0, 5.0, 0.03, dict(metric='wfr', freqs=(18.9, 22.4), seed=3)))
    assert seen['metric'] == 'wfr' and seen['ke'] == 5.0 and np.isclose(seen['m'], 6.0) and np.allclose(seen['noise']['w'], [18.9, 22.4])
    assert np.allclose(seen['noise']['ph'], NR.c_noise(0.03, dict(freqs=(18.9, 22.4), seed=3))['ph'])
    assert r['left'] == [-1, -1, -1] and r['w_left'] == [None] * 3 and r['opts']['freqs'] == [18.9, 22.4]


def test_c2_cases_and_their_report():
    """Referee round 63: C2's cases (four seeds per case, the draw's own pair with its phases), their grouping over the seeds and the
    jumps counted apart from the leaves."""
    cases = NR.c2_cases()
    labs = {c[0] for c in cases}
    assert labs == {'one tone', 'tone pair', "the draw's pair", 'band 25–75', 'm', 'WFR'}
    fr, ph = NR.draw_pair()
    assert abs(fr[0] - 18.9015) < 1e-3 and abs(fr[1] - 22.5019) < 1e-3 and abs((fr[1] - fr[0]) / CT.OM - 18.0019) < 1e-3
    assert sum(c[0] == "the draw's pair" for c in cases) == 4 and all('ph' in c[4] for c in cases if c[0] == "the draw's pair")
    pair = [c for c in cases if c[0] == 'tone pair' and c[4]['freqs'] == (18.9, 22.5) and c[3] == 0.02]
    assert sorted(c[4]['seed'] for c in pair) == sorted(NR.SEEDS_C)
    rs = [dict(label='tone pair', m=12.0, ke=10.0, s=0.08, opts=dict(freqs=[18.9, 22.5], seed=sd), left=lf, maxdev=md)
          for sd, lf, md in ((7, [124, -1, -1], [0.95, 0.1, 0.1]), (1, [69, -1, -1], [0.36, 0.1, 0.1]), (2, [-1, -1, -1], [0.2, 0.1, 0.1]))]
    md = '\n'.join(NR.c2_lines(rs))
    assert '| tone pair | paradigm | 12 | 10 | tones 18.9/22.5 (difference 18ω) | 0.08 | 3 | 2 (1), 0 (0), 0 (0) | 69–124 | 0.950, 0.100, 0.100 |' in md


def test_leave_stats_and_the_c3_report():
    """Referee round 63: at a leave the coefficient there and the mean offset over the periods leave − 30 … leave − 6; C3's rows."""
    Wp = np.tile(np.array([0.9, -0.6, 0.75]), (60, 1)) + np.linspace(0, 0.1, 60)[:, None] * np.array([1, 0, 0])[None, :]
    Cp = np.tile(CT.CSTAR, (60, 1))
    Cp[:, 0] += np.r_[np.full(40, -0.01), np.full(20, -0.9)]
    left = NR.left_periods(np.abs(Cp - CT.CSTAR), NR.LEFT_C)
    w_left, pre = NR.leave_stats(Wp, Cp, left)
    assert left == [41, -1, -1] and np.isclose(w_left[0], Wp[40, 0]) and np.isclose(pre[0], -0.01) and pre[1] is None
    md = '\n'.join(NR.c3_lines([dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=7), left=left, maxdev=[0.9, 0.0, 0.0], w_left=w_left,
                                      pre=pre, W=Wp.tolist(), C=Cp.tolist())]))
    assert "| D‴.49's 16 tones | 0.02 | 7 | 1 | 41 | 0.900 → 0.968 | -0.010 | 0.90 | -1.42 |" in md


def test_b5_splits_the_second_order_by_the_phase(monkeypatch):
    """Referee round 63: the second order at φ and φ + π/2 — the constant part (their half-sum) does not depend on φ, the 2ω part (their
    half-difference) turns by 2Δφ: at φ + π/4 the half-difference of the shifted pair is the first pair's turned by π/2."""
    monkeypatch.setattr(NR, 'NP_B', 2)
    monkeypatch.setattr(NR, 'MARKS_B', (1, 2))
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.37,)))
    run = lambda dph: NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', dict(nz, ph=np.asarray(nz['ph']) + dph), CT.KE, 2)[2]
    a, b, c, d = run(0.0), run(np.pi / 2), run(np.pi / 4), run(3 * np.pi / 4)
    assert np.allclose((a + b) / 2, (c + d) / 2, rtol=1e-6, atol=1e-9 * np.abs(a).max())
    r = NR.b5_job(('paradigm', 21.37))
    assert np.allclose(r['w2'], a) and np.allclose(r['w2_shift'], b) and len(r['w_env']) == 2
    md = '\n'.join(NR.b5_lines([r]))
    assert '| paradigm | 21.37 | 213.7 |' in md


def test_a2_rows_split_the_roll_off_from_the_resonance(monkeypatch):
    """Referee round 63: at K = 0 the jitter is the roll-off alone and the resonance factor at the smallest law nearer one; the peaks
    within 20% of √K_i."""
    monkeypatch.setattr(NR, 'M_A2', (3.0, 12.0))
    monkeypatch.setattr(NR, 'OM_A2', np.linspace(4.0, 40.0, 145))
    rows = NR.a2_rows()
    assert len(rows) == 4
    f = {(r['m'], tuple(r['band'])): np.array(r['jit']) / np.array(r['jit0']) for r in rows}
    assert np.all(f[(12.0, (10.0, 30.0))] > f[(3.0, (10.0, 30.0))]) and np.all(f[(3.0, (10.0, 30.0))] > 1)
    for r in rows:
        if r['m'] == 12.0:
            assert np.all(np.abs(np.array(r['peak']) / np.array(r['sqrtK']) - 1) < 0.2)
    md = '\n'.join(NR.a2_lines(rows))
    assert 'of which the roll-off' in md


def test_the_floquet_means_at_a_rotating_tone():
    """Referee round 69 (MO3, MO4, FM1–FM3): at 21.37 rad/s, where e^{iωT} is not ±1, the Floquet recursion's period means over three
    periods equal `second_order`'s first order (the rotation, Im G and Im H all enter) — to 2·10⁻³ of the largest: `second_order`'s
    trapezoid (2000 outputs a period) is exact for tones commensurate with the period (10⁻⁸ at 21.4) and errs by ~10⁻⁴ for the others,
    where a conjugated or dropped rotation errs by O(1)."""
    L, al0, c0, idx = _fit()
    nz = DZ.noise_of(1.0, dict(freqs=(21.37,)))
    a, ph = float(nz['amp'][0]), float(nz['ph'][0])
    M, P, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, 21.37, ph)
    Z = NR.floquet_means(M, P, a * G, a * H, 21.37, 3)
    W1, C1, _, _, _ = NR.second_order(L, NR.EPS, al0[idx], c0[idx], 'paradigm', nz, CT.KE, 3)
    for X, Y in ((Z[:, 1:4], W1), (Z[:, 4:7], C1)):
        assert np.abs(X - Y).max() < 2e-3 * np.abs(Y).max()
    assert abs(np.exp(1j * 21.37 * 2 * np.pi / CT.OM).real) < 0.99


def test_leave_stats_window_with_a_varying_offset():
    """Referee round 69 (LS1): the mean offset before a leave is over the periods leave − 30 … leave − 6 exactly."""
    n = 80
    Wp = np.tile(CT.ASTAR, (n, 1))
    Cp = np.tile(CT.CSTAR, (n, 1))
    Cp[:, 0] += np.r_[-0.001 * np.arange(50), np.full(n - 50, -0.9)]
    left = NR.left_periods(np.abs(Cp - CT.CSTAR), NR.LEFT_C)
    _, pre = NR.leave_stats(Wp, Cp, left)
    p = left[0]
    assert p == 51 and np.isclose(pre[0], np.mean(-0.001 * np.arange(p - 1 - 30, p - 1 - 5)))


def test_c_noise_negated_and_the_raised_start():
    """Referee rounds 69, 76: `c_noise`'s neg moves every phase by π (bands, seeds and tones alike); `raised` moves one formed kernel's
    coefficient by Δ in magnitude on its curve, keeping its sign, and nothing else."""
    for o in (dict(seed=3), dict(band=(25.0, 75.0), seed=2), dict(freqs=(13.1,), seed=7)):
        a, b = NR.c_noise(0.02, o), NR.c_noise(0.02, dict(o, neg=True))
        assert np.allclose(b['ph'], np.asarray(a['ph']) + np.pi) and np.allclose(a['w'], b['w']) and np.allclose(a['amp'], b['amp'])
    L, al0, c0, idx = _fit()
    for k in (0, 1, 2):
        al = NR.raised(L, al0, idx, k, 0.05)
        w0, w1 = NR.mobilities(L, NR.EPS, al0, 'paradigm')[0], NR.mobilities(L, NR.EPS, al, 'paradigm')[0]
        assert np.isclose(abs(w1[idx[k]]) - abs(w0[idx[k]]), 0.05, rtol=0, atol=1e-12) and np.sign(w1[idx[k]]) == np.sign(CT.ASTAR[k])
        assert np.array_equal(np.delete(al, idx[k]), np.delete(al0, idx[k]))
    al = NR.raised(L, al0, idx, 1, -0.05)
    assert np.isclose(abs(NR.mobilities(L, NR.EPS, al, 'paradigm')[0][idx[1]]), 0.55, rtol=0, atol=1e-9)


def test_the_basin_edge_brackets(monkeypatch):
    """Referee rounds 69, 76: without noise kernel 1 returns from a raise of 0.07 and is ejected from 0.075 (`edge_run`); with trace the
    anatomy on TRACE_K outputs a period — at 0.075 the centre stays 0.05 off its feature from t/T ≈ 0.594 on (non-returning), x_m below
    0, and ends below −1; at 0.07 it passes 0.05 off on both sides and comes back —; kernel 3 raised detaches at t/T ≈ 0.40, x_m above
    0; the reference's frequency reaches the loop (at ω_ref = 0.4 a raise of 0.075 returns, at 0.1 one of 0.05 is ejected); `edge_job`
    brackets the edge by the grid and bisections, every run on the bracket's side its outcome says, and records ω_ref and the anatomy at
    hi and at lo; EDGE_CASES hold the rows outside the fitted window."""
    p0, x0 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.07)
    p1, x1 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.075)
    assert p0 == -1 and x0 < 0.01 and 1 <= p1 <= 5 and x1 > NR.JUMP_C
    t1 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.075, NP=3, trace=True)
    t0 = NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.07, NP=3, trace=True)
    assert 0.55 < t1[2] < 0.65 and t1[3] < 0 and t1[5] < -1.0 and t1[6] < -0.5 and len(t1) == 8
    assert t0[2] is None and t0[3] is None and t0[6] < -0.05 and t0[7] > 0.04 and abs(t0[5] - CT.CSTAR[0]) < 0.01
    t3 = NR.edge_run(12.0, 10.0, 'paradigm', 2, 0.075, NP=3, trace=True)
    assert 0.35 < t3[2] < 0.45 and t3[3] > 0 and t3[0] >= 0
    assert NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.075, om=0.4)[0] == -1 and NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.05, om=0.1)[0] >= 0
    assert NR.edge_run(12.0, 10.0, 'paradigm', 0, 0.05)[0] == -1
    monkeypatch.setattr(NR, 'EDGE_GRID', (0.05, 0.1))
    monkeypatch.setattr(NR, 'BISECT_EDGE', 2)
    e = NR.edge_job((12.0, 10.0, 'paradigm', 0))
    assert np.isclose(e['lo'], 0.0625) and np.isclose(e['hi'], 0.075) and e['om'] == CT.OM
    for x in e['runs']:
        assert (x['left'] >= 0) == (x['d'] >= e['hi'])
    assert e['hi_trace']['t_detach'] is not None and e['lo_trace']['t_detach'] is None and e['t_detach'] == e['hi_trace']['t_detach']
    assert np.isclose(e['hi_trace']['t_detach'], t1[2]) and np.isclose(e['hi_trace']['c_end'], t1[5])
    monkeypatch.setattr(NR, 'EDGE_GRID', (0.1, 0.15))
    monkeypatch.setattr(NR, 'BISECT_EDGE', 0)
    e = NR.edge_job((12.0, 10.0, 'paradigm', 0, 0.4))
    assert e['om'] == 0.4 and e['lo'] == 0.1 and e['hi'] == 0.15
    for c in [(12.0, 2.5, 'paradigm', 0), (12.0, 80.0, 'paradigm', 0), (1.5, 10.0, 'paradigm', 0), (24.0, 10.0, 'paradigm', 0),
              (12.0, 20.0, 'wfr', 0), (12.0, 40.0, 'wfr', 0), (12.0, 10.0, 'paradigm', 0, 0.1), (3.0, 10.0, 'paradigm', 0, 0.4)]:
        assert c in NR.EDGE_CASES


def test_c3_job_with_a_fake_loop_and_a_raised_start(monkeypatch):
    """Referee rounds 69, 76 (C3's mutants): `c3_job` takes the period means of 8 outputs a period of the formed kernels — the coefficients
    as W, the centres as C —, the leave, the largest distance, the coefficient at the leave and the mean offset before it; opts rise starts
    the loop from `raised`."""
    seen = {}
    NPc, k = NR.NP_C, 8

    def fake(L, al0, c0, astar, ts, noise, metric, ke=None, x_init=None, om=None):
        seen.update(al0=np.array(al0), n=len(ts))
        n = len(ts)
        c = np.tile(c0, (n, 1)).astype(float)
        w = np.zeros((n, c0.size))
        _, _, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
        w[:, idx] = CT.ASTAR
        per = np.minimum(np.arange(n) // k, NPc - 1)
        c[:, idx[0]] += np.where(per >= 100, 0.5, 0.0)
        w[:, idx[0]] += 0.001 * per
        return dict(c=c, w=w, al=np.zeros((n, c0.size)), x=np.zeros(n))
    monkeypatch.setattr(NR, 'simulate', fake)
    r = NR.c3_job((12.0, 10.0, 0.02, dict(seed=7, rise=[0, 0.04])))
    L, al0, c0, idx = _fit()
    assert np.allclose(seen['al0'], NR.raised(NR.law(12.0), al0, idx, 0, 0.04)) and seen['n'] == NPc * k + 1
    assert r['left'] == [101, -1, -1] and np.isclose(r['maxdev'][0], 0.5) and np.isclose(r['w_left'][0], CT.ASTAR[0] + 0.1)
    assert np.isclose(r['pre'][0], 0.0) and np.allclose(np.asarray(r['W'])[:, 0], CT.ASTAR[0] + 0.001 * np.arange(NPc))
    C = np.asarray(r['C'])
    assert np.allclose(C[:100, 0], CT.CSTAR[0]) and np.allclose(C[100:, 0], CT.CSTAR[0] + 0.5) and np.allclose(C[:, 1:], CT.CSTAR[1:])


def test_w2_job_is_the_even_part(monkeypatch):
    """Referee round 69: `w2_job` is ((W₊ + W₋)/2 − A)/s² of the loop under ±S_W2·n, the −n run's phases moved by π."""
    calls = []

    def fake(L, metric, noise, ke, NP, chunk=10, k=2000):
        calls.append((metric, ke, NP, np.asarray(noise['ph'], float).copy()))
        sgn = 1.0 if len(calls) == 1 else -1.0
        W = np.tile(CT.ASTAR, (NP, 1)) + sgn * 0.01 + 2e-6 * np.arange(1, NP + 1)[:, None]
        return W, np.tile(CT.CSTAR, (NP, 1))
    monkeypatch.setattr(NR, 'sim_means', fake)
    monkeypatch.setattr(NR, 'NP_W2', 3)
    r = NR.w2_job((3.0, 5.0, 'wfr', dict(seed=2)))
    assert [c[:3] for c in calls] == [('wfr', 5.0, 3)] * 2 and np.allclose(calls[1][3], calls[0][3] + np.pi)
    assert np.allclose(np.asarray(r['even_w_all']), 2e-6 * np.arange(1, 4)[:, None] / NR.S_W2 ** 2 * np.ones(3))


def test_the_like_for_like_the_edge_and_the_raised_reports():
    """Referee rounds 69, 76: C4's row — the rise at the leave, the twin's half-sum (not where the twin left by the leave) and the draw's
    own s²⟨w₂⟩ — and its correlation (signed); EDGE's bracket mapped to the level s* = (edge/⟨w₂⟩)^{1/2} against the observed first jump,
    over the testable cases, with the permutation null; and C6's rows, on synthetic results."""
    n = 300
    W = np.tile(CT.ASTAR, (n, 1)); W[:, 0] += np.linspace(0, 0.09, n)
    Wm = np.tile(CT.ASTAR, (n, 1)); Wm[:, 0] += np.linspace(0, 0.07, n)
    C3 = [dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=5), left=[150, -1, -1], maxdev=[0.9, 0, 0], W=W.tolist())]
    C4 = [dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=5, neg=True), left=[-1, -1, -1], maxdev=[0.1, 0, 0], W=Wm.tolist())]
    ew = np.zeros((n, 3)); ew[:, 0] = np.linspace(100.0, 200.0, n)
    W2 = [dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=5), even_w_all=ew.tolist())]
    md = '\n'.join(NR.c4_lines(C3, C4, W2))
    rise, half, own = W[149, 0] - CT.ASTAR[0], (W[149, 0] + Wm[149, 0]) / 2 - CT.ASTAR[0], 0.02 ** 2 * ew[149, 0]
    assert (f"| 0.02 | 5 | 150 | {rise:+.4f} | -1 | {half:+.4f} ({half / rise:.3f}) | {own:+.4f} ({own / (0.02 ** 2 * ew[-1, 0]):.2f}) | "
            f"{rise / own:.3f} |") in md
    C4p = [dict(C4[0], left=[150, -1, -1])]
    assert f"| 0.02 | 5 | 150 | {rise:+.4f} | 150 | — |" in '\n'.join(NR.c4_lines(C3, C4p, W2))
    C3c, W2c = [], []
    for sd, p, e2 in ((1, 100, 300.0), (2, 140, 200.0), (3, 180, 100.0)):
        C3c.append(dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=sd), left=[p, -1, -1], maxdev=[0.9, 0, 0], W=W.tolist()))
        W2c.append(dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=sd), even_w_all=np.full((n, 3), e2).tolist()))
    rs = np.array([W[p - 1, 0] - CT.ASTAR[0] for p in (100, 140, 180)])
    cc = np.corrcoef(rs, 0.02 ** 2 * np.array([300.0, 200.0, 100.0]))[0, 1]
    assert cc < -0.99 and f"across the 3 draws {cc:.2f}" in '\n'.join(NR.c4_lines(C3c, C4, W2c))
    E_ = [dict(m=12.0, ke=10.0, metric='paradigm', k=0, lo=0.07, hi=0.075, runs=[dict(d=0.075, left=2, maxdev=0.9)])]
    W2b = [dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=sd), even_w_all=np.full((n, 3), 200.0).tolist()) for sd in NR.W2_SEEDS]
    C = [dict(m=12.0, ke=10.0, s=s_, opts=dict(seed=sd), left=[5 if s_ >= 0.02 else -1, -1, -1], maxdev=[0.9 if s_ >= 0.02 else 0, 0, 0])
         for sd in NR.W2_SEEDS for s_ in (0.015, 0.02)]
    md = '\n'.join(NR.edge_lines(E_, W2b, C, None))
    a, b = np.sqrt(0.07 / 200.0), np.sqrt(0.075 / 200.0)
    assert f"| 12 | 10 | paradigm | 0.2 | 1 | D‴.49's band | (0.07, 0.075] | 2 | — | ({a:.4f}, {b:.4f}], 0.02 |" in md
    assert ("consistent in 4 of the 4 testable draw-configurations (those that jumped, and those that did not where a level above s* "
            "was run), later in 0, earlier in 0; of those that never jumped, s* lay below the largest level run in 0. With the draws' ⟨w₂⟩ "
            "permuted within each configuration (200 permutations) the consistent count is 4.00 on average (4–4; at least 4 in 200)") in md
    W2c = W2b + [dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(band=(25.0, 75.0), seed=7), even_w_all=np.full((n, 3), 400.0).tolist())]
    C2 = [dict(m=12.0, ke=10.0, s=s_, opts=dict(band=[25.0, 75.0], seed=7), left=[5 if s_ >= 0.06 else -1, -1, -1],
               maxdev=[0.9 if s_ >= 0.06 else 0, 0, 0]) for s_ in (0.04, 0.06)]
    md = '\n'.join(NR.edge_lines(E_, W2c, C, C2))
    a2, b2 = np.sqrt(0.07 / 400.0), np.sqrt(0.075 / 400.0)
    assert f"| 12 | 10 | paradigm | 0.2 | 1 | 25–75 | (0.07, 0.075] | 2 | — | ({a2:.4f}, {b2:.4f}], 0.06 | —, — |" in md
    assert ("consistent in 4 of the 5 testable draw-configurations (those that jumped, and those that did not where a level above s* was "
            f"run), later in 1 (m = 12, k_e = 10, 25–75 (1; {0.06 / (0.5 * (a2 + b2)):.2f}–{0.06 / (0.5 * (a2 + b2)):.2f} times s*)), "
            "earlier in 0") in md
    C6 = [dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=5, rise=[0, 0.05]), left=[40, -1, -1], maxdev=[0.95, 0, 0])]
    assert "| 16 tones, draw 5 | 0.02 | 0.05 | 40 (0.950) | 150 |" in '\n'.join(NR.c6_lines(C6, C3))


def test_f_job_sums_and_saturations_at_a_rotating_tone(monkeypatch):
    """Referee round 69 (FJ4–FJ6): `f_job`'s Floquet sums for a tone are `floquet_means` with the tone's phase and amplitude (rms 1) over
    NP_B4 periods; at a commensurate tone its first three periods and the saturation P(I − M)⁻¹Re(aG) + Re(aH)."""
    monkeypatch.setattr(NR, 'FREQS_B4', ((21.37,),))
    monkeypatch.setattr(NR, 'NP_B4', 5)
    monkeypatch.setattr(NR, 'OM_COMM_F', (21.4,))
    r = NR.f_job(('paradigm', 12.0))
    L, al0, c0, idx = _fit()
    M, P = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE)
    for om, NP in ((21.37, 5), (21.4, 3)):
        nz = DZ.noise_of(1.0, dict(freqs=(om,)))
        a = float(nz['amp'][0])
        _, _, G, H = NR.monodromy(L, NR.EPS, al0[idx], c0[idx], 'paradigm', CT.KE, om, float(nz['ph'][0]))
        Z = NR.floquet_means(M, P, a * G, a * H, om, NP)
        if om == 21.37:
            assert np.allclose(r['sums']['21.37']['w'], Z[:, 1:4], rtol=1e-10, atol=0) and np.allclose(r['sums']['21.37']['c'], Z[:, 4:7])
        else:
            assert np.allclose(r['comm']['21.4']['first'], Z, rtol=1e-10, atol=0) and len(r['comm']['21.4']['first']) == 3
            sat = P @ np.linalg.solve(np.eye(M.shape[0]) - M, np.real(a * G)) + np.real(a * H)
            assert np.allclose(r['comm']['21.4']['sat'], sat, rtol=1e-10, atol=0)


def test_b5_job_phases_and_the_averaged_loop_per_s2(monkeypatch):
    """Referee round 69 (B5b–B5d): `b5_job` runs the second order at the tone's phase and at φ + π/2 and stores both runs' w and c, and the
    averaged loop at s = S_B[−1] and at S_B5 per s²; the report prints the small level's ratios."""
    calls, env = [], []
    nP = len(NR.MARKS_B)

    def fake_so(L, eps, al, c, metric, noise, ke, NP, k=2000):
        calls.append(float(np.asarray(noise['ph'])[0]))
        q = len(calls)
        return None, None, np.full((NP, 3), float(q)), np.full((NP, 3), 10.0 * q), 0

    def fake_env(L, eps, metric, noise, ke, NP, k=400, terms=NR.TERMS, quasi=False):
        env.append(float(np.asarray(noise['amp'])[0]))
        return np.tile(CT.ASTAR, (NP, 1)) + 0.5, np.tile(CT.CSTAR, (NP, 1)) - 0.25, 0
    monkeypatch.setattr(NR, 'second_order', fake_so)
    monkeypatch.setattr(NR, 'envelope_means', fake_env)
    r = NR.b5_job(('paradigm', 13.1))
    ph = float(DZ.noise_of(1.0, dict(freqs=(13.1,)))['ph'][0])
    assert np.isclose(calls[0], ph) and np.isclose(calls[1], ph + np.pi / 2)
    s = NR.S_B[-1]
    assert np.isclose(env[0], s * np.sqrt(2.0))
    assert np.allclose(r['w2'], 1.0) and np.allclose(r['w2_shift'], 2.0) and np.allclose(r['c2'], 10.0) and np.allclose(r['c2_shift'], 20.0)
    assert np.allclose(r['w_env'], 0.5 / s ** 2) and np.allclose(r['c_env'], -0.25 / s ** 2) and len(r['w2']) == nP
    assert len(env) == 2 and np.isclose(env[1], NR.S_B5 * np.sqrt(2.0))
    assert np.allclose(r['w_env_small'], 0.5 / NR.S_B5 ** 2) and np.allclose(r['c_env_small'], -0.25 / NR.S_B5 ** 2)
    md = '\n'.join(NR.b5_lines([r]))
    q = 0.5 / NR.S_B5 ** 2 / 1.5
    assert f"13.1: {q:.4f}–{q:.4f} (over the exact second order {0.5 / NR.S_B5 ** 2:.4f}–{0.5 / NR.S_B5 ** 2:.4f})" in md


def test_the_sensitivity_and_bodes_integral():
    """Referee round 69 (D3bp-4): A2's factor per tone and orbit point is the adaptation loop's sensitivity — the quasi-static transfer
    with K over that at K = 0, the plant's slope cancelling —; |S|² − 1 = K(2ω² − K)/|K − ω² + ik_eω|², so |S| > 1 exactly above
    √(K/2); |S(√K)| = (1 + K/k_e²)^{1/2}; Bode's integral ∫₀^∞ ln|S(iω)|dω = 0; the report prints them."""
    K, ke = 100.0, 10.0
    om = np.array([5.0, 7.0, 7.08, 10.0, 15.0, 30.0])
    S = NR.sens(K, ke, om)
    assert np.allclose(S ** 2 - 1, K * (2 * om ** 2 - K) / np.abs(K - om ** 2 + 1j * ke * om) ** 2, rtol=1e-12, atol=1e-14)
    assert np.array_equal(S > 1, om > np.sqrt(K / 2))
    assert np.isclose(NR.sens(K, ke, np.sqrt(K)), np.sqrt(1 + K / ke ** 2), rtol=1e-14)
    assert np.allclose(NR.sens(K, ke, [10.0, 15.0, 20.0, 30.0]), [1.41421, 1.38493, 1.24035, 1.11035], atol=1e-5)
    o = dict(f1=np.array([-3.0, 0.0, 2.5]), K=np.array([K, 40.0, 250.0]))
    o0 = dict(f1=o['f1'], K=np.zeros(3))
    for w_ in (4.0, 13.1, 37.0):
        r = np.abs(NR.transfer(o, np.array([w_]), ke)[0][:, 0] / NR.transfer(o0, np.array([w_]), ke)[0][:, 0])
        assert np.allclose(r, NR.sens(o['K'], ke, w_), rtol=1e-12)
    for K_, ke_ in ((248.3, 10.0), (100.0, 5.0), (10.0, 40.0), (1.0, 1.0)):
        assert abs(NR.bode_integral(K_, ke_)) < 1e-8
    rows = [dict(m=m, band=[10.0, 30.0], sqrtK=[sk, 1.0, 1.0], peak=[pk, 4.0, 5.0], jit=[1.0, 1.0, 1.0], jit0=[1.0, 1.0, 1.0])
            for m, sk, pk in ((4.83, 10.0, 8.35), (12.0, 15.76, 14.15))]
    md = '\n'.join(NR.a2_lines(rows))
    assert '|S(√K)| = (1 + K/k_e²)^{1/2} = 1.414, 1.866 at m = 4.83, 12' in md
    assert 'at m = 4.83 (√K₁ = 10.00) |S| = 1.41, 1.38, 1.24, 1.11 at ω = 10, 15, 20, 30' in md
    assert '4.00* (4.000)' in md and '8.35 (0.835)' in md and '5.00 (5.000)' in md and '5.00*' not in md


def test_c_noise_shifted_is_the_draw_from_a_later_start():
    """Referee round 69 (C7): opts shift starts the noise shift periods of the reference into the draw — each phase moved by ωt₀ —,
    so the shifted noise at t is the draw's at t + t₀; with neg and freqs as before."""
    T = 2 * np.pi / CT.OM
    a, b = NR.c_noise(0.02, dict(seed=7)), NR.c_noise(0.02, dict(seed=7, shift=25))
    n = lambda nz, t: float(np.sum(nz['amp'] * np.cos(nz['w'] * t + nz['ph'])))
    for t in (0.0, 1.3, 17.0):
        assert np.isclose(n(b, t), n(a, t + 25 * T), atol=1e-12)
    assert np.allclose(b['w'], a['w']) and np.allclose(b['amp'], a['amp'])
    t1 = NR.c_noise(0.02, dict(freqs=(13.1,), seed=7, shift=3))
    t0 = NR.c_noise(0.02, dict(freqs=(13.1,), seed=7))
    assert np.isclose(t1['ph'][0] - t0['ph'][0], 13.1 * 3 * T)
    assert NR.c_noise(0.0, dict(seed=7, shift=50)) is None and NR.c_noise(0.02, dict(seed=7, shift=0))['ph'] is not None


def test_bracket_edge_on_a_step_and_its_rechecks():
    """Referee round 69: `bracket_edge` goes up the grid to the first ejecting raise, rechecks the next `above` raises, bisects nb times,
    and keeps every run; a step at 0.0731 is bracketed to the grid step over 2^nb; nothing ejecting gives hi None."""
    calls = []

    def run(d):
        calls.append(d)
        return (2 if d > 0.0731 else -1), (0.9 if d > 0.0731 else 0.001)
    lo, hi, runs = NR.bracket_edge(run, grid=(0.02, 0.04, 0.06, 0.08, 0.1, 0.13), nb=6, above=2)
    assert lo < 0.0731 < hi and np.isclose(hi - lo, 0.02 / 64)
    assert calls[:6] == [0.02, 0.04, 0.06, 0.08, 0.1, 0.13] and len(calls) == 12
    assert [x[0] for x in runs] == sorted(calls)
    lo, hi, runs = NR.bracket_edge(lambda d: (-1, 0.0), grid=(0.1, 0.2), nb=3)
    assert hi is None and lo == 0.2 and len(runs) == 2
    lo, hi, runs = NR.bracket_edge(lambda d: (1, 1.0), grid=(0.1, 0.2), nb=2, above=1)
    assert lo == 0.0 and np.isclose(hi, 0.025) and [x[0] for x in runs] == [0.025, 0.05, 0.1, 0.2]


def test_the_edges_scaling_with_the_centre_mobility():
    """Referee rounds 69, 76 (EDGE): `v_of` is a formed kernel's centre mobility at the fit — Gmw under the paradigm's metric (∝ m), G·a
    under WFR's —; on edges following Δ_e = 0.08(k_e/v₁)^{1/2} exactly the window's power law has α = β = 1/2 (kernel 3, WFR and ω_ref
    rows kept out of the fit; the bracket's midpoint), and each row outside the window prints against the law 0.08(k_e/v)^{1/2} and,
    for the paradigm's kernel 1, the power law; the ω_ref slopes; `star_verdict` sorts a first jump against s*'s bracket."""
    v12, v3 = NR.v1_of(12.0, 'paradigm'), NR.v1_of(3.0, 'paradigm')
    assert np.isclose(v12, 12.0 * CT.ASTAR[0], rtol=1e-6) and np.isclose(v12 / v3, 4.0, rtol=1e-9)
    assert 0.9 < NR.v1_of(12.0, 'wfr') < 1.1 and np.isclose(NR.v_of(12.0, 'paradigm', 2), 12.0 * CT.ASTAR[2], rtol=1e-6)
    row = lambda m, ke, d, w=1e-7, **kw: dict(dict(m=m, ke=ke, metric='paradigm', k=0, lo=d - w, hi=d + w, runs=[]), **kw)
    E_ = [row(m, ke, 0.08 * np.sqrt(ke / NR.v1_of(m, 'paradigm')), w=0.004) for m in (3.0, 12.0) for ke in (5.0, 10.0, 20.0)]
    d3 = 0.5 * np.sqrt(10.0 / NR.v_of(12.0, 'paradigm', 2))
    E_ += [row(12.0, 10.0, d3, k=2), row(12.0, 80.0, 1.5 * 0.08 * np.sqrt(80.0 / v12), w=0.01),
           dict(m=12.0, ke=40.0, metric='wfr', k=0, lo=1.0, hi=None, runs=[]),
           row(12.0, 10.0, 0.08 * np.sqrt(10.0 / v12) * 0.5 ** 0.75, om=0.1), row(12.0, 10.0, 0.08 * np.sqrt(10.0 / v12) * 2 ** 0.75, om=0.4),
           row(24.0, 10.0, 0.9 * 0.08 * np.sqrt(10.0 / NR.v1_of(24.0, 'paradigm')), w=0.002)]
    md = '\n'.join(NR.edge_lines(E_, None, None, None))
    assert f"m = 24, k_e = 10: {0.9 * 0.08 * np.sqrt(10.0 / NR.v1_of(24.0, 'paradigm')):.4f}, 0.90 times the law, the power law's 0.90" in md
    assert 'over the paradigm\'s 6 kernel-1 rows in m = 3–12, k_e = 5–40' in md and 'has α = 0.500 ±' in md and 'β = 0.500 ±' in md
    assert 'Δ_e·(v₁/k_e)^{1/2} = 0.0800 … 0.0800 over the window (median 0.0800)' in md
    assert f"m = 12, k_e = 80: {1.5 * 0.08 * np.sqrt(80.0 / v12):.4f}, 1.50 times the law, the power law's 1.50" in md
    assert f"m = 12, k_e = 10, kernel 3: {d3:.4f}, 6.25 times the law" in md
    assert f"m = 12, k_e = 40, wfr: > 1 (the law {0.08 * np.sqrt(40.0 / NR.v1_of(12.0, 'wfr')):.3f})" in md
    assert "the log-slope by least squares (pairwise): m = 12, k_e = 10: 0.75 (0.75, 0.75)" in md
    fit = NR.edge_fit(E_)
    assert len(fit['rows']) == 6 and np.isclose(fit['alpha'], 0.5) and np.isclose(fit['beta'], 0.5)
    assert NR.edge_fit(E_[6:]) is None
    assert NR.star_verdict(0.018, 0.019, 0.02, [0.015, 0.02]) == 'consistent'
    assert NR.star_verdict(0.035, 0.036, 0.08, [0.04, 0.06, 0.08]) == 'later'
    assert NR.star_verdict(0.008, 0.009, 0.0075, [0.005, 0.0075]) == 'earlier'
    assert NR.star_verdict(0.07, 0.072, None, [0.04, 0.06, 0.08]) == 'none, predicted'
    assert NR.star_verdict(0.09, 0.1, None, [0.04, 0.06, 0.08]) == 'untested'


def test_c2_single_tones_rise_on_the_raw():
    """Referee round 69: C2 prints, per level, the first coefficient's rise over the fit by the last period in the single-tone runs (the
    paradigm's metric) whose first kernel did not jump — the drift that never reaches the tones' moving edge."""
    raw = NR.E.load(os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'noise_resonance_raw.npz'))
    md = '\n'.join(NR.c2_lines(raw['C2']))
    rise = {}
    for r in raw['C2']:
        o = r['opts']
        if 'freqs' in o and len(o['freqs']) == 1 and o.get('metric', 'paradigm') == 'paradigm' and r['m'] == 12.0 and r['ke'] == 10.0:
            if not (r['left'][0] >= 0 and r['maxdev'][0] > NR.JUMP_C):
                rise.setdefault(r['s'], []).append(float(r['w_end'][0] - CT.ASTAR[0]))
    assert rise and all(len(v) == 24 for v in rise.values())
    for s_, v in rise.items():
        assert f"s = {s_:g}: {min(v):+.3f} … {max(v):+.3f} ({len(v)} runs)" in md


def test_star_cases_count_and_null():
    """Referee round 76 (s*, ED01–ED05): kernel 3's s* from its own ⟨w₂⟩ and its own jumps; `star_count`'s verdicts with replaced ⟨w₂⟩;
    `star_null` permutes the draws' ⟨w₂⟩ within each configuration only."""
    n = 300
    E_ = [dict(m=12.0, ke=10.0, metric='paradigm', k=0, lo=0.07, hi=0.075, runs=[]),
          dict(m=12.0, ke=10.0, metric='paradigm', k=2, lo=0.074, hi=0.075, runs=[])]
    W2 = []
    for i, sd in enumerate(NR.W2_SEEDS):
        ew = np.zeros((n, 3)); ew[:, 0] = 200.0; ew[:, 2] = 100.0 * (i + 1)
        W2.append(dict(m=12.0, ke=10.0, metric='paradigm', opts=dict(seed=sd), even_w_all=ew.tolist()))
    C = [dict(m=12.0, ke=10.0, s=s_, opts=dict(seed=sd), left=[5 if s_ >= 0.02 else -1, -1, 5 if s_ >= 0.03 else -1],
              maxdev=[0.9 if s_ >= 0.02 else 0, 0, 0.9 if s_ >= 0.03 else 0]) for sd in NR.W2_SEEDS for s_ in (0.015, 0.02, 0.03)]
    cases = NR.star_cases(E_, W2, C, None)
    k3 = [q for q in cases if q['key'][3] == 2]
    assert [q['w2'] for q in k3] == [100.0, 200.0, 300.0, 400.0] and all(q['obs'] == 0.03 for q in k3)
    assert all(q['obs'] == 0.02 for q in cases if q['key'][3] == 0 and q['key'][4] is None)
    assert all(q['obs'] is None and not q['run'] for q in cases if q['key'][4] is not None)   # W2's band rows: no runs here
    rec = NR.star_count(cases)
    v3 = [v for q, v, _ in rec if q['key'][3] == 2]
    sstar = [np.sqrt(0.0745 / x) for x in (100.0, 200.0, 300.0, 400.0)]       # 0.0273, 0.0193, 0.0158, 0.0136
    assert v3 == ['consistent', 'later', 'later', 'later'] and sstar[0] > 0.02 > sstar[1]
    alt = NR.star_count(cases, [q['w2'] if q['key'][3] == 0 else 100.0 for q in cases])
    assert [v for q, v, _ in alt if q['key'][3] == 2] == ['consistent'] * 4
    null = NR.star_null(cases, n=50, seed=1)
    assert len(null) == 50 and all(x == 5 for x in null)
    md = '\n'.join(NR.edge_lines(E_, W2, C, None))
    w1, w3 = 0.005 / (0.5 * (np.sqrt(0.07 / 200) + np.sqrt(0.075 / 200))), 0.01 / (0.5 * (np.sqrt(0.074 / 100) + np.sqrt(0.075 / 100)))
    assert f"are {w1:.2f}–{w3:.2f} of s* wide" in md and (f"{w1:.2f}", f"{w3:.2f}") == ('0.26', '0.37')
    assert "the consistent count is 5.00 on average (5–5; at least 5 in 200)" in md


def test_edge_fit_standard_errors():
    """Referee round 76: `edge_fit`'s coefficients and standard errors are ordinary least squares' with n − 3 degrees of freedom, over
    the window's rows only."""
    rng = np.random.default_rng(3)
    rows = []
    for m in (3.0, 6.0, 12.0):
        for ke in (5.0, 10.0, 20.0, 40.0):
            d = 0.05 * ke ** 0.5 * m ** -0.5 * np.exp(0.05 * rng.standard_normal())
            rows.append(dict(m=m, ke=ke, metric='paradigm', k=0, lo=d - 1e-9, hi=d + 1e-9, runs=[]))
    fit = NR.edge_fit(rows + [dict(rows[0], m=24.0), dict(rows[1], ke=80.0), dict(rows[2], metric='wfr'), dict(rows[3], k=2)])
    X = np.array([[1.0, np.log(r['ke']), np.log(r['m'])] for r in rows])
    y = np.log([0.5 * (r['lo'] + r['hi']) for r in rows])
    cf = np.linalg.solve(X.T @ X, X.T @ y)
    res = y - X @ cf
    se = np.sqrt(np.diag(res @ res / (len(y) - 3) * np.linalg.inv(X.T @ X)))
    assert len(fit['rows']) == 12 and np.allclose([fit['c'], fit['alpha'], -fit['beta']], cf, rtol=1e-12)
    assert np.allclose(fit['se'], se, rtol=1e-9) and np.allclose(fit['res'], res, atol=1e-12)


def test_bands_of_and_scan_sets():
    """Referee round 76: `bands_of` gives the maximal runs of flagged raises; `scan_sets` the leaving and the jumping bands of a scan."""
    ds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    assert NR.bands_of(ds, [False, True, False, True, True, True]) == [(0.1, 0.1), (0.3, 0.5)]
    assert NR.bands_of(ds, [False] * 6) == [] and NR.bands_of(ds, [True] * 6) == [(0.0, 0.5)]
    r = dict(runs=[dict(d=d, left=(3 if d in (0.1, 0.3, 0.4, 0.5) else -1), maxdev=(0.9 if d >= 0.4 else 0.4)) for d in reversed(ds)])
    assert NR.scan_sets(r) == ([(0.1, 0.1), (0.3, 0.5)], [(0.4, 0.5)])
    assert NR.fmt_bands([(0.1, 0.1), (0.3, 0.5)]) == '0.1, 0.3–0.5' and NR.fmt_bands([]) == 'none'


def test_c7_job_scans_and_lines_with_a_fake_run(monkeypatch):
    """Referee rounds 69, 76 (C7): `c7_job` runs `edge_run` under the case's noise (the window's shift, the tone's sign) over NP_C7 periods
    at every raise of SCAN_C7 (kernel 1; m = 12, k_e = 10, the paradigm's metric unless given); `c7_lines` prints the noise-free sets, the
    tones' least raises under ±n with the even part, E and the odd part ('>1' kept), and the windows' least raises with their band counts,
    their least and median ('>1' above every raise) and the windows that are not an interval."""
    seen = []
    monkeypatch.setattr(NR, 'SCAN_C7', (0.0, 0.065, 0.0675, 0.07, 0.0725, 0.1, 0.2, 0.5, 1.0))
    monkeypatch.setattr(NR, 'c_noise', lambda s, o: None if s == 0 else dict(
        tag=(tuple(o.get('freqs', ())), bool(o.get('neg')), o.get('shift', 0), s)))

    def fake(m, ke, metric, k, d, NP=NR.NP_EDGE, trace=False, noise=None, om=None):
        seen.append((m, ke, metric, k, NP, None if noise is None else noise['tag']))
        if noise is None:
            ej, md = d >= 0.0725, 0.9
        else:
            fr, neg, sh, s = noise['tag']
            if fr:
                th = {(False, 0.01): 0.1, (True, 0.01): 0.2, (False, 0.02): 2.0, (True, 0.02): 0.5}[(neg, s)]
                ej, md = d >= th, 0.9
            else:
                ej = sh == 50 and (d == 0.0675 or d >= 0.1)
                md = 0.4 if d == 0.0675 else 0.9
        return (3 if ej else -1), (md if ej else 0.01)
    monkeypatch.setattr(NR, 'edge_run', fake)
    cases = [(dict(seed=7), 0.0), (dict(freqs=(13.1,), seed=7), 0.01), (dict(freqs=(13.1,), seed=7, neg=True), 0.01),
             (dict(freqs=(13.1,), seed=7), 0.02), (dict(freqs=(13.1,), seed=7, neg=True), 0.02),
             (dict(seed=7, shift=0), 0.02), (dict(seed=7, shift=50), 0.02)]
    rs = [NR.c7_job(c) for c in cases]
    assert all(x[:5] == (12.0, 10.0, 'paradigm', 0, NR.NP_C7) for x in seen) and len(seen) == 9 * len(cases)
    assert {x[5] for x in seen if x[5] is not None and x[5][0] == ()} == {((), False, 0, 0.02), ((), False, 50, 0.02)}
    assert rs[6]['opts'] == dict(seed=7, shift=50) and [x['d'] for x in rs[6]['runs']] == list(NR.SCAN_C7)
    md = '\n'.join(NR.c7_lines(rs))
    assert "| 12 | 10 | paradigm | 0.0725–1 | 0.0725–1 |" in md
    assert "| one tone 13.1 | 0.01 | 0.1 (0.1–1) | 0.2 (0.2–1) | 0.15 | +775 | -5.000 |" in md
    assert "| one tone 13.1 | 0.02 | >1 (none) | 0.5 (0.5–1) | >1 | — | — |" in md
    assert "E per tone and level (the even part's rise per s²; '—' where a least raise is beyond the grid): one tone 13.1: +775 at 0.01." in md
    assert "| 12 | 10 | paradigm | 0.02 | >1, 0.0675 [2] | 0.0675, >1 | 1 of 2 | 0.1, >1 |" in md
    assert ("Without noise the leaving set is an interval from the edge to the grid's end at every configuration; under the band's noise it is not an interval in 1 of the 2 windows: m = 12, k_e = 10, s = 0.02, from period 50: 0.0675, 0.1–1.") in md
    rr = NR.c7_job((dict(seed=7, shift=25), 0.04, (12.0, 20.0, 'paradigm')))
    assert seen[-1][:4] == (12.0, 20.0, 'paradigm', 0) and rr['ke'] == 20.0 and rr['s'] == 0.04


def test_c8_verdict_confusion_and_flags():
    """Referee round 76 (C8): a run's first flag against its first leave — hit in the leave's window (p < P ≤ p + STEP_C8), early, miss;
    for a run that never left correct or a false alarm —; the windows' confusion counts up to the leave's window; the rules' flags."""
    st = list(range(0, 300, 10))
    assert NR.c8_verdict(133, [130]) == ('hit', 130) and NR.c8_verdict(130, [120]) == ('hit', 120) and NR.c8_verdict(131, [130]) == ('hit', 130)
    assert NR.c8_verdict(133, [120, 130]) == ('early', 120) and NR.c8_verdict(133, [140, 150]) == ('miss', None)
    assert NR.c8_verdict(-1, []) == ('correct', None) and NR.c8_verdict(-1, [50, 60]) == ('false alarm', 50)
    assert list(NR.c8_confusion(133, [120, 130, 200], st)) == [1, 0, 1, 12]
    assert list(NR.c8_confusion(133, [], st)) == [0, 1, 0, 13] and list(NR.c8_confusion(-1, [50], st)) == [0, 0, 1, 29]
    r = dict(wins=[dict(p=0, d=0.0, left=-1, maxdev=0.01), dict(p=10, d=0.08, left=4, maxdev=0.4), dict(p=20, d=0.06, left=2, maxdev=0.9)])
    assert NR.c8_flags(r, 'leave') == [10, 20] and NR.c8_flags(r, 'jump') == [20]
    assert NR.c8_flags(r, 'the noise-free edge', 0.07) == [10] and NR.c8_flags(r, 'the noise-free edge', None) == []


def test_c8_job_and_lines_with_a_fake_run(monkeypatch):
    """Referee round 76 (C8): `c8_job` runs C's loop (`c_run`, the metric passed) and every STEP_C8 periods p `edge_run` with kernel 1 raised
    by the run's own δ(p) = |W₁| over period p minus its fitted value (0 at p = 0) under the run's noise shifted by p, over NP_C7 periods;
    `c8_lines` prints each run's row and the rules' tallies in and out of sample."""
    L, al0, c0, idx = _fit()
    w0 = abs(float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]]))
    Wp = np.tile(CT.ASTAR, (NR.NP_C, 1)) + np.linspace(0.0, 0.15, NR.NP_C)[:, None] * np.array([1.0, 0, 0])[None, :]
    Cp = np.tile(CT.CSTAR, (NR.NP_C, 1)); Cp[132:, 0] -= 0.9
    got = {}

    def fake_run(m, ke, s, opts):
        got.update(m=m, ke=ke, s=s, opts=dict(opts))
        return None, L, al0, c0, idx, Wp, Cp
    calls = []
    monkeypatch.setattr(NR, 'c_run', fake_run)
    monkeypatch.setattr(NR, 'c_noise', lambda s, o: dict(shift=o.get('shift', 0), s=s))
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=NR.NP_EDGE, trace=False, noise=None, om=None: (
        calls.append((m, ke, mt, k, NP, noise['shift'], d)), ((3, 0.9) if d >= 0.055 else (-1, 0.01)))[1])
    r = NR.c8_job((12.0, 10.0, 'wfr', 0.02, dict(seed=3)))
    assert got == dict(m=12.0, ke=10.0, s=0.02, opts=dict(seed=3, metric='wfr')) and r['left'] == [133, -1, -1] and np.isclose(r['w0'], w0)
    assert [c[5] for c in calls] == list(range(0, NR.NP_C, NR.STEP_C8)) and all(c[:5] == (12.0, 10.0, 'wfr', 0, NR.NP_C7) for c in calls)
    assert calls[0][6] == 0.0 and np.isclose(calls[13][6], abs(Wp[129, 0]) - w0) and len(r['wins']) == len(calls)
    first = next(x['p'] for x in r['wins'] if x['left'] >= 0)
    rr = [r, dict(r, opts=dict(seed=4)), dict(r, opts=dict(seed=7), left=[-1, -1, -1], maxdev=[0.1, 0, 0])]
    E_ = [dict(m=12.0, ke=10.0, metric='wfr', k=0, lo=0.09, hi=0.11, runs=[])]
    md = '\n'.join(NR.c8_lines(rr, E_))
    assert first == 120 and "| 12 | 10 | wfr | 0.02 | 3 | 133 j |" in md
    assert "120, 130, 140, 150 … (18): early" in md and "| leave | out | 0 | 2 | 0 | 0 | 0 |" in md and "| 12 | 10 | wfr | 0.02 | 4 | 133 j |" in md
    assert "| leave | in | 0 | 0 | 0 | 0 | 1 |" in md and "| the noise-free edge | out |" in md
    assert ("The totals: in sample 1 runs, 0 of them leaving, and 30 earlier windows; out of sample 2 runs, 2 of them leaving, and 26 "
            "earlier windows.") in md


def test_b5_lines_summary_over_every_mark():
    """Referee round 76 (B5L1): B5's summary at the small level spans every mark and kernel, not the last mark's."""
    nP = len(NR.MARKS_B)
    a = np.ones((nP, 3)); es = np.ones((nP, 3)); es[0, 0] = 0.9; es[-1, 2] = 1.1
    r = dict(metric='paradigm', om=13.1, w2=a.tolist(), w2_shift=a.tolist(), w_env=a.tolist(), w_env_small=es.tolist())
    md = '\n'.join(NR.b5_lines([r]))
    assert "13.1: 0.9000–1.1000 (over the exact second order 0.9000–1.1000)" in md


def test_a2_rows_reference_includes_the_plant_slope(monkeypatch):
    """Referee round 76 (AR1): A2's K = 0 reference is the roll-off (iω − f̂′)/(−ω² + ik_eω), the plant's slope f̂′ included."""
    from src.core import closed_loop as CL
    monkeypatch.setattr(NR, 'M_A2', (12.0,))
    monkeypatch.setattr(NR, 'OM_A2', np.linspace(4.0, 40.0, 13))
    rows = NR.a2_rows()
    L = NR.law(12.0)
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    xm = np.sin(np.linspace(0, 2 * np.pi, 2001)[:-1])
    o = NR.orbit_terms(L, NR.EPS, al0[idx], c0[idx], 'paradigm', xm)
    nz = CL.tones(1.0, *NR.BANDS_A2[0], seed=7)
    om, amp = np.asarray(nz['w'], float), np.asarray(nz['amp'], float)
    Y0 = (1j * om[None, :] - o['f1'][:, None]) / (-om[None, :] ** 2 + 1j * CT.KE * om[None, :])
    g = o['v'][None, :] * o['d1']
    j0 = np.sqrt(np.mean((np.abs(g[:, :, None] * Y0[:, None, :] * amp[None, None, :]) ** 2 / 2).sum(2), axis=0))
    assert np.allclose(rows[0]['jit0'], j0, rtol=1e-12) and np.abs(o['f1']).max() > 0.1


def test_w2_cases_cover_the_band_draws_and_every_configuration():
    """Referee round 76 (W2C): W2's draws are C3's eleven band draws at m = 12, k_e = 10 and W2_SEEDS for every other configuration."""
    cs = NR.w2_cases()
    band = [c for c in cs if c[:3] == (12.0, 10.0, 'paradigm') and 'band' not in c[3]]
    assert sorted(c[3]['seed'] for c in band) == list(range(1, 12))
    others = {(m, ke, mt, b) for m, ke, mt, b in NR.W2_CONFIGS} - {(12.0, 10.0, 'paradigm', None)}
    assert len(cs) == 11 + len(NR.W2_SEEDS) * len(others)


def _c9_runs():
    """Synthetic C3/C6 runs of draw 3 for C9: kernel 1's |w| = w₀ + 0.001·period (+ add), leaving at P; and runs C9 must leave out."""
    L = NR.law(NR.C7_CFG[0])
    al0, c0, idx = DZ.fitted_state(L, NR.EPS, CT.CSTAR, CT.ASTAR)
    w0 = float(NR.mobilities(L, NR.EPS, al0, 'paradigm')[0][idx[0]])

    def run(opts, P, s=0.02, add=0.0):
        W = np.tile(CT.ASTAR, (NR.NP_C, 1)).astype(float)
        W[:, 0] = np.sign(w0) * (abs(w0) + 0.001 * np.arange(NR.NP_C) + add)
        return dict(opts=opts, s=s, left=[P, -1, -1], W=W.tolist())
    C3 = [run(dict(seed=3), 157), run(dict(seed=3), 155, s=0.03), run(dict(freqs=(13.1,), seed=3), 155),
          run(dict(seed=3, neg=True), 155), run(dict(seed=7), 155)]
    C6 = [run(dict(seed=3, rise=[0, 0.05]), 166, add=0.02), run(dict(seed=3, rise=[0, 0.03]), 120, add=0.01),
          run(dict(seed=3, rise=[0, 0.02]), 150, add=0.005), run(dict(seed=3, rise=[0, 0.04]), 160, add=0.015)]
    return C3, C6, abs(w0)


def test_c9_pairs_own_job_and_part(monkeypatch):
    """Referee round 76 (C9D): `c9_pairs` takes, per window (seed, p), the draw's s = 0.02 runs without tones or −n that have not left by
    p (a leave at p itself excluded), with δ(p) = |W₁| over the period before p minus the fitted value; `c9_own_job` runs `edge_run` at
    C7_CFG with kernel 1 raised by δ under the draw's noise from p for NP_C7 periods; `part_C9D` maps it over every pair of the default
    windows."""
    C3, C6, w0 = _c9_runs()
    pr = NR.c9_pairs(C3, C6, [(3, 150), (3, 160), (3, 165), (3, 100)])
    assert [(x[0], x[1], x[2], x[4]) for x in pr] == [(3, 150, 0.0, 157), (3, 150, 0.04, 160), (3, 150, 0.05, 166), (3, 160, 0.05, 166),
                                                      (3, 165, 0.05, 166), (3, 100, 0.0, 157), (3, 100, 0.02, 150), (3, 100, 0.03, 120),
                                                      (3, 100, 0.04, 160), (3, 100, 0.05, 166)]
    assert np.allclose([x[3] for x in pr], [0.149, 0.164, 0.169, 0.179, 0.184, 0.099, 0.104, 0.109, 0.114, 0.119], atol=1e-12)
    calls = []
    monkeypatch.setattr(NR, 'c_noise', lambda s, o: dict(s=s, **o))
    monkeypatch.setattr(NR, 'edge_run', lambda m, ke, mt, k, d, NP=NR.NP_EDGE, trace=False, noise=None, om=None: (
        calls.append((m, ke, mt, k, d, NP, noise)), (4, 0.7))[1])
    r = NR.c9_own_job((3, 150, 0.05, 0.169))
    assert calls == [(12.0, 10.0, 'paradigm', 0, 0.169, NR.NP_C7, dict(s=0.02, seed=3, shift=150))]
    assert r == dict(seed=3, p=150, rise=0.05, d=0.169, left=4, maxdev=0.7)
    got = []
    monkeypatch.setattr(NR, '_pool_map', lambda fn, jobs: got.append((fn, list(jobs))) or [])
    NR.part_C9D(C3, C6)
    assert got[0][0] is NR.c9_own_job and got[0][1] == [x[:4] for x in NR.c9_pairs(C3, C6)]
    assert {(x[0], x[1]) for x in got[0][1]} <= {(o['seed'], int(o['shift'])) for o, _ in NR.C9_CASES} and len(got[0][1]) > 0


def test_c9_lines_counts_and_the_runs_flagged():
    """Referee round 76 (C9, C9D): the table's cells (δ, whether the nearest raise leaves — a leave without a jump counts —, the loop at
    δ, the run — its leave at the window's end counts —), the 2×2 counts against the run's leave by the nearest raise and by the loop at
    δ, the misses, and the runs flagged in some window containing their leave (an earlier window's flag counts)."""
    C3, C6, _ = _c9_runs()
    sc = lambda sh, lv, md=None: dict(opts=dict(seed=3, shift=sh), s=0.02, runs=[
        dict(d=d, left=l_, maxdev=(md or {}).get(d, 0.9 if l_ >= 0 else 0.01)) for d, l_ in zip((0.14, 0.15, 0.16, 0.17, 0.18, 0.19), lv)])
    C9 = [sc(165, (-1, -1, -1, -1, -1, 1)), sc(160, (-1, -1, -1, -1, 2, 2)), sc(150, (3, 3, -1, 2, -1, -1), {0.17: 0.4})]
    md = '\n'.join(NR.c9_lines(C9, C3, C6))
    assert "| 3 | 150 | 0.14–0.15, 0.17 | 0: +0.1490, yes, 157; 0.04: +0.1640, no, 160; 0.05: +0.1690, yes, — |" in md
    assert "| 3 | 160 | 0.18–0.19 | 0.05: +0.1790, yes, 166 |" in md and "| 3 | 165 | 0.19 | 0.05: +0.1840, no, 166 |" in md
    assert md.index("| 3 | 150 |") < md.index("| 3 | 160 |") < md.index("| 3 | 165 |")
    assert ("Of the 5 (run, window) pairs, by the grid's nearest raise: 2 in the leaving set with the run leaving there, 1 in the set but "
            "staying, 2 leaving though not in the set, 0 neither.") in md and "By the window's loop" not in md
    D = [dict(seed=3, p=150, rise=0.0, d=0.149, left=6, maxdev=0.9), dict(seed=3, p=150, rise=0.04, d=0.164, left=9, maxdev=0.9),
         dict(seed=3, p=150, rise=0.05, d=0.169, left=-1, maxdev=0.01), dict(seed=3, p=160, rise=0.05, d=0.179, left=6, maxdev=0.9),
         dict(seed=3, p=165, rise=0.05, d=0.184, left=-1, maxdev=0.01)]
    md = '\n'.join(NR.c9_lines(C9, C3, C6, D))
    assert "| 3 | 150 | 0.14–0.15, 0.17 | 0: +0.1490, yes, 156, 157; 0.04: +0.1640, no, 159, 160; 0.05: +0.1690, yes, —, — |" in md
    assert "| 3 | 160 | 0.18–0.19 | 0.05: +0.1790, yes, 166, 166 |" in md and "| 3 | 165 | 0.19 | 0.05: +0.1840, no, —, 166 |" in md
    assert ("By the window's loop at δ itself: 3 leaving with the run, 0 leaving while the run stays, 1 staying while the run leaves, "
            "1 neither; the run leaving though the loop at δ stays: draw 3 raised by 0.05 in the window from period 165 (it leaves at "
            "166). Of the 3 runs whose leave a window contains, 3 are flagged — the loop at δ leaving — in at least one of those windows.") in md
    D[3] = dict(D[3], left=-1)
    md = '\n'.join(NR.c9_lines(C9, C3, C6, D))
    assert ("2 leaving with the run, 0 leaving while the run stays, 2 staying while the run leaves, 1 neither; the run leaving though the "
            "loop at δ stays: draw 3 raised by 0.05 in the window from period 160 (it leaves at 166); draw 3 raised by 0.05 in the window "
            "from period 165 (it leaves at 166). Of the 3 runs whose leave a window contains, 2 are flagged — the loop at δ leaving — in at "
            "least one of those windows: not draw 3 raised by 0.05.") in md
    D[2] = dict(D[2], left=2)
    assert "2 leaving with the run, 1 leaving while the run stays, 2 staying while the run leaves, 0 neither" in '\n'.join(
        NR.c9_lines(C9, C3, C6, D))
    assert NR.c9_lines([], C3, C6, D) == [] and NR.c9_lines(None) == []

def test_main_passes_c3_and_c6_to_c9d(monkeypatch, tmp_path):
    """Referee round 76: `main --only C9D` computes C9D from the raw's C3 and C6 and keeps every other part."""
    old = dict(C3=['c3'], C6=['c6'], C9=['c9'])
    saved, got = [], []
    monkeypatch.setattr(NR.E, 'load', lambda p: dict(old))
    monkeypatch.setattr(NR.E, 'save', lambda p, **kw: saved.append(kw))
    monkeypatch.setattr(NR, 'OUT', str(tmp_path))
    monkeypatch.setattr(NR, 'report', lambda **kw: '')
    monkeypatch.setattr(NR, 'part_C9D', lambda C3, C6: got.append((C3, C6)) or ['d'])
    monkeypatch.setattr(sys, 'argv', ['exp_noise_resonance.py', '--only', 'C9D'])
    NR.main()
    assert got == [(['c3'], ['c6'])] and saved[-1] == dict(C3=['c3'], C6=['c6'], C9=['c9'], C9D=['d'])
