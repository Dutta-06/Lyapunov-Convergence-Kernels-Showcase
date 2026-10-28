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
    assert '| 12 | 10 | 15.76 | 0.02 (1/1) / 0.02 | — / — | 0.03 (1/1) / 0.03 | 3 (3 leaves) |' in md


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
