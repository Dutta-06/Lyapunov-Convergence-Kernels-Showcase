"""Proposition D‴.58: the experiment's analysis (`experiments/exp_width_mismatch.py`)."""
import os
import sys
import warnings
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

import exp_width_mismatch as WM
from src.core import wide_feature as W


def test_flank_peaks_of_the_one_kernel_optimum():
    """B's `flank`: after the optimal one-kernel fit on the samples (centre at the feature, least-squares amplitude) the residual's
    correlation C is nonnegative to rounding and peaks symmetrically ≈ 2σ beyond the node, on both flanks."""
    ss = WM.sigma_star(1.0)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    phi = np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.SIG ** 2))
    w = float(phi @ f / (phi @ phi))
    mn, left, right = WM.flank(1.0, np.array([WM.C_B]), np.array([w]))
    assert abs(mn) < 1e-12 and abs(left - right) < 0.02 and 2.0 < left < 2.5


def test_lead_coefficients_and_the_tail_job():
    """D: the ratios' leading coefficients −C_{j−1}/2^{j−1} (j ≤ 3) with a vanishing next divided difference; `d_job` is the library's
    tail ratio."""
    lc = WM.lead_coeffs(3)
    assert [Fraction(*r['lead']) for r in lc] == [Fraction(-1), Fraction(-1, 2), Fraction(-1, 2)] and all(r['next_dd'] == 0 for r in lc)
    r = WM.d_job((6, 0.5))
    assert r['n_terms'] == 58 and abs(r['ratio'] - W.tail_ratio(6, 0.5, n_terms=58)) < 1e-15


def test_the_discrepancy_stop_prediction():
    """C's prediction: the first K whose rule error e_K satisfies e_K²rf² ≤ s²((1 + δ)² − 1 + 2K/N); larger noise stops earlier."""
    rf, N = 0.4, WM.X_B.size
    ks = [WM.predict_stop(s, rf, N)[0] for s in (1e-4, 1e-3, 1e-2)]
    assert ks[0] > ks[1] > ks[2] >= 1
    K, err = WM.predict_stop(1e-3, rf, N)
    ss = WM.sigma_star(WM.TAU_C)
    e = W.l2_error(*W.half_squeeze(K, WM.SIG, ss), WM.SIG, ss)
    e_prev = W.l2_error(*W.half_squeeze(K - 1, WM.SIG, ss), WM.SIG, ss)
    bound = lambda k: 1e-6 * ((1 + WM.DELTA) ** 2 - 1 + 2 * k / N)
    assert (e * rf) ** 2 <= bound(K) and (e_prev * rf) ** 2 > bound(K - 1)
    assert abs(err - np.sqrt(e ** 2 + 2 * K * 1e-6 / (N * rf ** 2))) < 1e-15


def test_sigma_star_and_the_blocks():
    """σ* = σ(1 + (τ/σ)²)^{1/2}; `blocks`: hist[i + 1] is iteration i's block (its active kernels sorted by centre), relaxed[i] its flag
    (False past the flags' end)."""
    assert abs(WM.sigma_star(1.0) - WM.SIG * np.sqrt(2)) < 1e-15 and abs(WM.sigma_star(0.0) - WM.SIG) < 1e-15
    c = np.array([0.5, -0.2, 0.1])
    R = dict(hist=[(0.0, c, np.zeros(3), 1.0), (1.0, c, np.array([0.0, 0.3, 0.0]), 0.5), (2.0, c, np.array([0.2, 0.3, 0.0]), 0.1)],
             relaxed=[True])
    b = WM.blocks(R)
    assert [x['k'] for x in b] == [1, 2] and [x['relaxed'] for x in b] == [True, False]
    assert np.allclose(b[1]['c'], [-0.2, 0.5]) and np.allclose(b[1]['w'], [0.3, 0.2]) and b[1]['rms'] == 0.1 and b[0]['tau'] == 1.0


def test_sample_fit_sorts_and_projects():
    """`sample_fit` from an unsorted start returns sorted centres, the least-squares coefficients and their rms residual."""
    x = WM.X_B
    phi = lambda c: np.exp(-(x[:, None] - np.asarray(c)[None, :]) ** 2 / (2 * WM.SIG ** 2))
    f = phi([-0.3, 0.2]) @ np.array([0.5, -0.4])
    c, a, rms = WM.sample_fit(phi, f, [0.25, -0.35])
    assert np.all(np.diff(c) > 0) and np.allclose(c, [-0.3, 0.2], atol=1e-8) and np.allclose(a, [0.5, -0.4], atol=1e-8) and rms < 1e-10


def test_the_stop_prediction_uses_2k_over_n():
    """A case where the bound with 2K/N stops at K = 3 and one with K/N would not (N = 4: the variance term dominant)."""
    rf, N = 0.25, 4
    ss_ = WM.sigma_star(WM.TAU_C)
    e3 = W.l2_error(*W.half_squeeze(3, WM.SIG, ss_), WM.SIG, ss_)
    s = e3 * rf / np.sqrt((1 + WM.DELTA) ** 2 - 1 + 1.5 * 3 / N)
    K, err = WM.predict_stop(s, rf, N)
    assert K == 3 and abs(err - np.sqrt(e3 ** 2 + 2 * 3 * s ** 2 / (N * rf ** 2))) < 1e-15


def test_the_tail_job_has_enough_terms_and_the_degree_check_can_fail(monkeypatch):
    """`d_job` at t = 0.8 sums 40/ln(1/t) terms (180), its ratio converged; `lead_coeffs`' (j + 1)-th divided difference is nonzero for
    ratios of degree j + 1 (a fake defect table with P_1(K) = K²)."""
    r = WM.d_job((10, 0.8))
    assert r['n_terms'] == 180 and abs(r['ratio'] - W.tail_ratio(10, 0.8, n_terms=400)) < 1e-12
    fake = lambda K, n: [0] * n + [K * K if n == 2 * K + 2 else 1]
    monkeypatch.setattr(WM.W, 'gh_defects', fake)
    lc = WM.lead_coeffs(1)
    assert lc[0]['next_dd'] != 0 and Fraction(*lc[0]['lead']) == 3


def test_the_report_on_synthetic_results():
    """A2's leading coefficients −1, −1, −3, −15; B's optimal blocks within OPT_TOL and the first non-optimal one; C's stop at the first
    rms within (1 + δ)s; D's Catalan verdict, false for a wrong coefficient."""
    rows = []
    for r_ in WM.TAUS_A:
        P = W.params(WM.SIG, WM.sigma_star(r_))
        rows.append(dict(r=r_, K=1, t=P['t'], lam=P['lam'], e_hs=0.1, e_series=0.1, e_first=0.099, e_mm=0.2, e_ls=0.09, y_hs=np.array([0.0])))
    ratios = {K: [W.gh_defect(2 * K + 2 * j, K) // W.gh_defect(2 * K, K) for j in range(1, 5)] for K in range(1, 7)}
    share = [dict(r=r_, K=K, share=0.01) for r_ in WM.TAUS_A for K in (5, 10, 20, 30)]
    md = WM.report(A=dict(rows=rows, ratios=ratios, share=share))
    lines = md.splitlines()
    i = next(i for i, l in enumerate(lines) if l.startswith('| j | K = 1, 2, 3, 4, 5'))
    assert [l.rsplit('|', 2)[1].strip() for l in lines[i + 2:i + 6]] == ['-1', '-1', '-3', '-15']
    blk = lambda k, e, e_opt, rel: dict(k=k, tau=float(k), c=WM.C_B + np.linspace(-0.1, 0.1, k), w=np.full(k, 0.8 / k), e=e, e_opt=e_opt,
                                        c_opt=None, relaxed=rel, dc=0.0, dc_hs=0.0, e_hs=e)
    B = [dict(r=0.25, var='base', blocks=[blk(1, 1e-2 * (1 + 5e-4), 1e-2, True), blk(2, 1e-3 * (1 + 2e-3), 1e-3, False)],
              events=[(1.0, 10, 'enter', 0.0), (2.0, 11, 'enter', 0.0)], relaxed=[True, False], stiff=0, redone=0, failed=0, rejected=0)]
    md = WM.report(B=B)
    assert '| 0.25 | base | 1 | 2, 1.00e-03, 1.002, False | — | 1/2 | 0, 0, 0 | 0 | 0 |' in md
    B[0]['blocks'][1]['polish'] = [dict(status=0, nfev=5000, accepted=True), dict(status=1, nfev=12, accepted=False, e0=1.0, e1=1.0)]
    B[0]['entries'] = [dict(it=0, j=10, u0=0.0, rate=0.5, rate_max=0.5, t_wall=2.0, t_min=2.0),
                       dict(it=1, j=11, u0=0.9, rate=0.1, rate_max=0.4, t_wall=1.0, t_min=1.0)]
    md = WM.report(B=B)
    assert '| 2, 1.00e-03, 1.002, False | 0/5000/a 1/12/r | 1/2 |' in md                 # the first non-optimal block's polishes (round 54)
    assert ('the entrant had the largest current rate |⟨rφ_j⟩| in 1; its rate over the largest a median 0.62 (range 0.250–1.00), its '
            'dual |u_j| at the block\'s start a median 0.450 (least 0.000)') in md
    assert ('Of the 1 entrants joining a block, 1 start beyond its outermost node (a median 0.58σ, range 0.58–0.58σ, beyond it) and 0 '
            'between its nodes.') in md
    x = dict(e_hs=1.0, e_best=0.5, e_opt=0.4)
    assert WM.gap_share(x) == '83.3%' and WM.gap_share(dict(x, e_opt=0.6)) == '100.0%, optimum unconverged'
    assert WM.gap_share(dict(x, e_opt=0.5 * (1 + 1e-12))) == '100.0%'                 # within UNCONV: not flagged
    C = [dict(s=1e-3, q=0, k=[1, 2, 3], rms=[5e-3, 1.05e-3, 0.9e-3], err=[0.02, 0.005, 0.004], at=[1, 2, 3], gaps=[0.5, 0.2], exits=0,
              stiff=0, runaway=False, relaxed=[True] * 3, rf=0.25)]
    md = WM.report(C=C)
    assert '| 0.001 | 2 [2, 2] (1/1 stopped) |' in md and '| 1 [1, 1] (the stop in 0) |' in md
    D = dict(rows=[dict(K=10, t=0.5, ratio=1.0, n_terms=58)], lead=[dict(j=1, lead=(-1, 1), next_dd=0.0), dict(j=2, lead=(-1, 2), next_dd=0.0)])
    assert 'at every j' in WM.report(D=D) and 'NOT' not in WM.report(D=D)
    D['lead'].append(dict(j=3, lead=(1, 2), next_dd=0.0))
    assert 'NOT at every j' in WM.report(D=D)


def test_the_quoted_cells_on_synthetic_results():
    """Referee round 54: the cells the text quotes, on synthetic results. A3: the ratios per kernel at the last K against their rates, the
    best variance at K = 1 and at the last K with an optimum, the gap shares and their range (one optimum unconverged); A4: the errors and
    the last ratio. B: the stall, failure and failed-chunk counts (distinct), the first non-optimal block with its polishes, the summary
    over every width and variant, the base path's table (an exit among the events, ratios only after one more kernel) and the flank of
    its optimal blocks. C: the stop at the first rms within (1 + δ)s, the error there over the oracle, over the bias–variance prediction
    from B's base path (not the shifted grid's) and over the bound, the largest gap's kernel count, the runaway count, a level where no
    draw stops."""
    rows, narrow = [], []
    for r_ in WM.TAUS_A:
        P = W.params(WM.SIG, WM.sigma_star(r_))
        x0 = dict(r=r_, t=P['t'], lam=P['lam'], y_hs=np.array([0.0]))
        rows += [dict(x0, K=1, e_hs=0.1, e_series=0.1, e_first=0.099, e_mm=0.5, e_ls=0.09, e_fock=0.2, e_mercer=0.3, e_best=0.09,
                      v_best=0.98765, e_opt=0.08),
                 dict(x0, K=2, e_hs=0.04, e_series=0.04, e_first=0.0399, e_mm=0.25, e_ls=0.035, e_fock=0.05, e_mercer=0.1, e_best=0.03,
                      v_best=0.97531, e_opt=0.035 if r_ == 2.0 else 0.025),
                 dict(x0, K=3, e_hs=0.016, e_series=0.016, e_first=0.0159, e_mm=0.125, e_ls=0.015, e_fock=0.0125, e_mercer=0.0333,
                      e_best=0.012, v_best=0.96)]
    for q in WM.NARROW_A:
        narrow += [dict(q=q, K=K, t=(1 - q) / (1 + q), e=e, e_series=e) for K, e in zip((1, 2, 3), (0.1 * q, 0.03 * q, 0.009 * q * q))]
    ratios = {K: [W.gh_defect(2 * K + 2 * j, K) // W.gh_defect(2 * K, K) for j in range(1, 5)] for K in range(1, 7)}
    share = [dict(r=r_, K=K, share=0.01) for r_ in WM.TAUS_A for K in (5, 10, 20, 30)]
    md = WM.report(A=dict(rows=rows, ratios=ratios, share=share, narrow=narrow))
    assert 'over all K with an optimum the shares 50.0%–100.0% (1 optima unconverged):' in md
    for r_ in WM.TAUS_A:
        P = W.params(WM.SIG, WM.sigma_star(r_))
        b2 = np.sqrt(1 + 2 * r_ ** 2)
        last = '0.87500 (100.0%, optimum unconverged)' if r_ == 2.0 else '0.62500 (66.7%)'
        assert (f"| {r_:g} | 2.500 ({1 / P['t']:.3f}) | 4.000 ({(2 - P['lam']) / P['lam']:.3f}) | 3.003 ({(b2 + 1) / (b2 - 1):.3f}) | "
                f"0.98765, 0.90000, 0.80000 (50.0%) | 2, 0.97531, 0.75000, {last} |") in md
    for q in WM.NARROW_A:
        assert (f"| {q:g} | {(1 + q) / (1 - q):.4f} | 1: {0.1 * q:.3e} (0e+00), 2: {0.03 * q:.3e} (0e+00), 3: {0.009 * q * q:.3e} (0e+00) | "
                f"{0.03 / (0.009 * q):.4f} |") in md

    def blk(k, tau, e, e_opt, rel, dc=0.0, dc_hs=0.0, e_hs=None, shift=0.0, **kw):
        return dict(k=k, tau=tau, c=WM.C_B + shift + 0.2 * WM.SIG * (np.arange(k) - (k - 1) / 2), w=np.full(k, WM.A_B / k), e=e,
                    e_opt=e_opt, c_opt=None, relaxed=rel, dc=dc, dc_hs=dc_hs, e_hs=e if e_hs is None else e_hs, **kw)
    pol = [dict(status=0, nfev=5000, accepted=True), dict(status=1, nfev=12, accepted=False, e0=1.0, e1=1.0)]
    path = lambda r_, var, es: dict(r=r_, var=var, blocks=[blk(k, 10.0 + k, e * (1 + 1e-4), e, True) for k, e in es],
                                    events=[(10.0 + k, 9 + k, 'enter', 0.0) for k, e in es], relaxed=[True] * len(es), stiff=0, redone=0,
                                    failed=0, rejected=0)
    B = [dict(r=0.25, var='base', blocks=[blk(1, 10.72, 1.0005e-2, 1e-2, True, dc=1.2e-8, e_hs=1.1e-2, shift=0.03),
                                          blk(2, 11.87, 2.0002e-3, 2e-3, True, dc=2.4e-8, dc_hs=0.012, e_hs=2.5e-3, shift=-0.02),
                                          blk(3, 12.85, 3.006e-8, 3e-8, False, dc=3.6e-6, dc_hs=0.024, e_hs=5e-4, polish=pol),
                                          blk(2, 13.5, 2.1e-3, 2e-3, True, e_hs=2.5e-3)],
              events=[(10.72, 10, 'enter', 0.0), (11.87, 11, 'enter', 0.1), (12.85, 12, 'enter', 0.2), (13.5, 11, 'exit', 0.1)],
              relaxed=[True, True, False, True], stiff=1, redone=2, failed=3, rejected=4),
         dict(r=0.25, var='grid +0.05', blocks=[blk(1, 10.7, 1.0003e-2, 1e-2, True), blk(2, 11.9, 6e-7, 5e-7, True)],
              events=[(10.7, 10, 'enter', 0.05), (11.9, 11, 'enter', 0.15)], relaxed=[True, True], stiff=0, redone=0, failed=0, rejected=0),
         path(1.0, 'base', ((1, 0.02), (2, 0.005), (3, 0.004))), path(1.0, 'grid +0.05', ((1, 0.03), (2, 0.007), (3, 0.006)))]
    rf, N, s = 0.25, WM.X_B.size, 1e-3
    C = [dict(s=s, k=[1, 2, 3, 4], rms=[5e-3, 1.15e-3, 1.08e-3, 0.95e-3], err=[0.02, 0.006, 0.0042, 0.004], at=[1, 2, 3, 3],
              gaps=[0.5, 0.9, 0.2], runaway=False, rf=rf),
         dict(s=s, k=[1, 2, 3, 4], rms=[4e-3, 1.09e-3, 1.0e-3, 0.9e-3], err=[0.03, 0.0055, 0.005, 0.0049], at=[1, 2, 2, 2],
              gaps=[0.2, 0.9, 0.5], runaway=True, rf=rf),
         dict(s=s, k=[1, 2, 3, 4], rms=[5e-3, 3e-3, 2e-3, 1.2e-3], err=[0.02, 0.01, 0.008, 0.007], at=[1, 1, 1, 1], gaps=[0.4, 0.3, 0.2],
              runaway=False, rf=rf),
         dict(s=1e-2, k=[1, 2], rms=[0.05, 0.02], err=[0.2, 0.1], at=[1, 1], gaps=[0.5], runaway=False, rf=rf)]
    md = WM.report(B=B, C=C)
    assert '| 0.25 | base | 1, 2 | 3, 3.00e-08, 1.002, False | 0/5000/a 1/12/r | 3/4 | 1, 2, 3 | 4 | 1 |' in md
    assert '| 0.25 | grid +0.05 | 1 | 2, 5.00e-07, 1.2, True | — | 2/2 | 0, 0, 0 | 0 | 0 |' in md
    assert '| 1 | base | 1, 2, 3 | none | — | 3/3 | 0, 0, 0 | 0 | 0 |' in md
    assert ("Over every width and variant: the optimal blocks' optimum errors down to 2.00e-03; the non-optimal blocks' from 2.00e-03 down "
            "(1 of 3 failed the stationarity test); non-optimal blocks with an optimum above 10⁻⁷: 2; the non-optimal blocks' errors "
            "3.01e-08 … 2.10e-03.") in md
    assert '### τ/σ = 0.25, base (events: 10+, 11+, 12+, 11−)' in md
    for row in ('| 1 | 10.72 | 1.0005e-02 | 1.0000e-02 | 1.000500 | True | 1.0e-07 | 0.000 | 1.1000e-02 | — |',
                f'| 2 | 11.87 | 2.0002e-03 | 2.0000e-03 | 1.000100 | True | 2.0e-07 | 0.100 | 2.5000e-03 | {1.0005e-2 / 2.0002e-3:.3f} |',
                f'| 3 | 12.85 | 3.0060e-08 | 3.0000e-08 | 1.002000 | False | 3.0e-05 | 0.200 | 5.0000e-04 | {2.0002e-3 / 3.006e-8:.3f} |',
                '| 2 | 13.5 | 2.1000e-03 | 2.0000e-03 | 1.050000 | True | 0.0e+00 | 0.000 | 2.5000e-03 | — |'):
        assert row in md

    def flank(r_, c, w):                                   # C(c) on the fine grid, by a loop
        x = WM.X_B
        r = WM.A_B * np.exp(-(x - WM.C_B) ** 2 / (2 * WM.sigma_star(r_) ** 2))
        for ck, wk in zip(c, w):
            r = r - wk * np.exp(-(x - ck) ** 2 / (2 * WM.SIG ** 2))
        cg = np.linspace(x[0], x[-1], 3201)
        Cg = np.array([np.mean(r * np.exp(-(x - g) ** 2 / (2 * WM.SIG ** 2))) for g in cg])
        lft, rgt = cg < min(c), cg > max(c)
        return (Cg.min() / Cg.max(), (min(c) - cg[lft][np.argmax(Cg[lft])]) / WM.SIG,
                (cg[rgt][np.argmax(Cg[rgt])] - max(c)) / WM.SIG)
    fl = [flank(0.25, o['c'], o['w']) for o in B[0]['blocks'][:2]]
    assert all(f"{x[1]:.2f}" != f"{x[2]:.2f}" for x in fl)                      # off-centre blocks: the two flanks differ
    assert ("after the optimal blocks (K = 1, 2): min C/max C within ±" + f"{max(abs(x[0]) for x in fl):.0e}, its peaks beyond the "
            "outermost nodes at (left, right) " + ', '.join(f"({x[1]:.2f}, {x[2]:.2f})σ" for x in fl) + '.') in md
    fmt = lambda v, f_: f"{np.median(v):{f_}} [{np.min(v):{f_}}, {np.max(v):{f_}}]"
    ks, es = [3, 2], [0.0042, 0.0055]
    eopt = {2: 0.005, 3: 0.004}
    Kp, ep = WM.predict_stop(s, rf, N)
    assert (f"| 0.001 | 2.5 [2, 3] (2/3 stopped) | {Kp} | 4.850e-03 [4.200e-03, 5.500e-03] | {ep:.3e} | 2.5 [2, 3] | 2 [2, 2] (the stop "
            "in 1) |") in md
    r1 = [e / (s * np.sqrt(2 * k / N) / rf) for k, e in zip(ks, es)]
    r2 = [e / np.sqrt(eopt[k] ** 2 + 2 * k * s ** 2 / (N * rf ** 2)) for k, e in zip(ks, es)]
    r3 = [e / (s * np.sqrt((1 + WM.DELTA) ** 2 - 1 + 4 * k / N) / rf) for k, e in zip(ks, es)]
    assert f"| 0.001 | {fmt(r1, '.2f')} | {fmt(r2, '.3f')} | {max(r3):.3f} | 1 |" in md
    Kp, ep = WM.predict_stop(1e-2, rf, N)
    assert f"| 0.01 | — (0/1 stopped) | {Kp} | — | {ep:.3e} | — | — |" in md and '| 0.01 | — | — | — | — |' in md


def test_the_stop_law():
    """Referee round 54: `stop_law`'s K(s) is the root of the first-term formula, its ceiling the predicted stop at the three levels; the
    stop's probabilities sum to one but the chance of no stop by K = 12 and are largest at the predicted K; the report prints them with the expected and observed counts."""
    ss = WM.sigma_star(WM.TAU_C)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    rf, N, t = float(np.sqrt(np.mean(f ** 2))), WM.X_B.size, W.params(WM.SIG, ss)['t']
    for s in WM.LEVELS_C:
        Kr, P = WM.stop_law(s, rf, N)
        rhs = np.log((1 - t * t) ** 0.25 * (np.pi * Kr) ** 0.25 * rf / (s * np.sqrt(2 * WM.DELTA + WM.DELTA ** 2))) / np.log(1 / t)
        Kp = WM.predict_stop(s, rf, N)[0]
        assert abs(Kr - rhs) < 1e-12 and int(np.ceil(Kr)) == Kp and max(P, key=P.get) == Kp
        assert 1 - 1e-3 < sum(P.values()) <= 1 + 1e-12 and min(P.values()) >= 0           # the rest: no stop by K = 12
    C = [dict(s=1e-3, k=[3, 4], rms=[2e-3, 1e-3], err=[0.01, 0.005], at=[3, 4], gaps=[0.5], runaway=False, rf=rf) for _ in range(2)]
    Kr, P = WM.stop_law(1e-3, rf, N)
    md = WM.report(C=C)
    assert (f"| 0.001 | {Kr:.2f} | " + ', '.join(f"{K}: {p:.3f}" for K, p in P.items() if p >= 5e-4) + f" | {P[4]:.3f} | "
            f"{2 * P[4]:.1f} of 2 | 2 |") in md


def test_b2_variable_projection_from_the_blocks():
    """B2 (round 54): from the half-squeeze nodes themselves both starts reach one optimum; the summary counts the blocks in the optimum's
    basin overall and each path's first non-optimal block (B2 in the paths' order)."""
    yh, _ = WM.W.half_squeeze(2, WM.SIG, WM.sigma_star(1.0))
    q = WM.b2_job((1.0, 'base', 2, WM.C_B + np.sort(yh)))
    assert abs(q['e_block'] / q['e_rule'] - 1) < 1e-9 and 0 < q['e_rule'] < 0.1
    B2 = [dict(r=1.0, var='base', k=10, e_block=1.0, e_rule=1.0), dict(r=1.0, var='base', k=11, e_block=5.0, e_rule=1.0),
          dict(r=0.5, var='base', k=6, e_block=3.0, e_rule=1.0)]
    line = WM.b2_line(B2)[0]
    assert line == ("Variable projection from the non-optimal block's own centres (B2) reaches the best optimum found in 1 of 3; in the others it "
                    "ends 3–5 times above it. From the first non-optimal block of each path: 1 of 2 (the others 3 times above it).")
    assert WM.b2_line(None) == []
