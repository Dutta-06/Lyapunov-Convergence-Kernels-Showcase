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
    rfB = float(np.sqrt(np.mean((WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))) ** 2)))
    B[0]['blocks'][1]['polish'] = [dict(status=0, nfev=5000, accepted=True),
                                   dict(status=1, nfev=12, accepted=False, e0=1.0, e1=1.002e-3 * rfB)]  # its end over the optimum (round 65)
    B[0]['entries'] = [dict(it=0, j=10, u0=0.0, rate=0.5, rate_max=0.5, t_wall=2.0, t_min=2.0),
                       dict(it=1, j=11, u0=0.9, rate=0.1, rate_max=0.4, t_wall=1.0, t_min=1.0)]
    md = WM.report(B=B)
    assert '| 2, 1.00e-03, 1.002, False | 0/5000/a 1/12/r (1.002) | 1/2 |' in md        # the first non-optimal block's polishes (round 54)
    assert '| none (the first entry) | 1 (1) | 1 | 1.00 (1.000–1.00) | 0.000 (0.000) | — | — |' in md   # the entries by regime (round 60)
    assert '| an optimal block | 1 (1) | 0 | 0.25 (0.250–0.25) | 0.900 (0.900) | 1, 0.58σ (0.58–0.58σ) | — |' in md
    assert '| all | 2 (2) | 1 | 0.62 (0.250–1.00) | 0.450 (0.000) | 1, 0.58σ (0.58–0.58σ) | — |' in md and '| a non-optimal block |' not in md
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
    rfB = float(np.sqrt(np.mean((WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))) ** 2)))
    pol = [dict(status=0, nfev=5000, accepted=True), dict(status=1, nfev=12, accepted=False, e0=1.0, e1=3.006e-8 * rfB)]
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
    assert '| 0.25 | base | 1, 2 | 3, 3.00e-08, 1.002, False | 0/5000/a 1/12/r (1.002) | 3/4 | 1, 2, 3 | 4 | 1 |' in md
    assert '| 0.25 | grid +0.05 | 1 | 2, 5.00e-07, 1.2, True | — | 2/2 | 0, 0, 0 | 0 | 0 |' in md
    assert '| 1 | base | 1, 2, 3 | none | — | 3/3 | 0, 0, 0 | 0 | 0 |' in md
    mg = sorted(WM.cs_margin(b['r'], b['var'], o) for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= WM.OPT_TOL and o['relaxed'])
    assert ("Over every width and variant: the optimal blocks' optimum errors down to 2.00e-03; the non-optimal blocks' from 2.00e-03 down "
            f"(1 of 3 failed the stationarity test; for the other 2 its threshold over the Cauchy–Schwarz bound — the residual's cosine with "
            f"the block's directions it allows — is {mg[0]:.1e}–{mg[1]:.2f} (median {np.median(mg):.1e}; above 0.1 in "
            f"{sum(m > 0.1 for m in mg)}): a weak test); non-optimal blocks with an optimum above 10⁻⁷: 2; the non-optimal blocks' errors "
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
    """B2 (rounds 54 and 60): from the half-squeeze nodes themselves both starts reach one optimum, at the block's own nodes; the summary
    counts the blocks in the optimum's basin only where the optimum is at least E_FLOOR (below it the solve stops short), the misses
    within 1% (in its basin, unconverged) apart from the far ones with their nodes' distance, and each path's first non-optimal block —
    a path being a width and a variant (B2 in the paths' order)."""
    from src.core import rich_transit as RT
    yh, _ = WM.W.half_squeeze(2, WM.SIG, WM.sigma_star(1.0))
    c = WM.C_B + np.sort(yh)
    q = WM.b2_job((1.0, 'base', 2, c))
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    assert abs(q['e_block'] / q['e_rule'] - 1) < 1e-9 and 0 < q['e_rule'] < 0.1
    assert q['dc'] == float(np.abs(c - WM.sample_fit(RT.Bumps(WM.X_B, WM.SIG).phi, f, c)[0]).max()) and 0 < q['dc'] < 0.01 * WM.SIG
    assert 1 < q['cond_vp'] < q['cond_joint'] < 1e12
    B2 = [dict(r=1.0, var='base', k=10, e_block=1e-6, e_rule=1.5e-6, dc=0.0),
          dict(r=1.0, var='base', k=11, e_block=1.005e-6, e_rule=1e-6, dc=0.001),
          dict(r=1.0, var='grid +0.05', k=10, e_block=3e-6, e_rule=1e-6, dc=0.3),
          dict(r=0.5, var='base', k=8, e_block=5e-12, e_rule=1e-12, dc=0.2),
          dict(r=0.5, var='base', k=9, e_block=2e-11, e_rule=2e-11, dc=0.0),
          dict(r=0.5, var='polish 2σ', k=8, e_block=4e-6 * (1 + 1e-4), e_rule=4e-6, dc=0.0)]
    assert WM.b2_line(B2)[0] == (
        "Variable projection from the non-optimal block's own centres (B2) reaches the best optimum found in 2 of the 4 blocks whose "
        "optimum is at least 1e-10; 1 ends within 1% of it (1.005 times: in its basin, unconverged); 1 ends 3 times above it — the "
        "blocks' largest node distances from the optimum's 2.50σ —. The 2 blocks with an optimum below 1e-10, where the solve stops short of the "
        "optimum (A), are not counted (1 reached). From the first non-optimal block of each path: 2 of 3 (the others 3 times above it); on "
        "1 path the first non-optimal block's optimum is below 1e-10.")
    assert WM.b2_line(None) == []
    for q, (cj, cv, sm) in zip(B2, ((1e7, 1e4, 1.5), (2e7, 2e4, 9.0), (4e8, 3e5, 2.0), (9e9, 9e9, 7.0), (9e9, 9e9, 7.0), (5e7, 2e5, 1.0))):
        q.update(cond_joint=cj, cond_vp=cv, smin_ratio=sm, vp_status=3, vp_nfev=77, vp_gap=0.002, vp_cos=3e-4)
    md = WM.b2_line(B2)[0]
    assert md.endswith("on 1 path the first non-optimal block's optimum is below 1e-10.")
    for q, (cx, vx, rx, re) in zip(B2, ((1e6, 2e3, 2.5, 80.0), (2e6, 3e3, 9.0, 9.0), (5e8, 4e4, 1.2, 2000.0), (9e9, 9e9, 7.0, 7.0),
                                        (9e9, 9e9, 7.0, 7.0), (3e6, 1e4, 1.6, 49.0))):
        q.update(cond_joint_x=cx, cond_vp_x=vx, smin_ratio_x=rx, smin_ratio_eq=re)
    md = WM.b2_line(B2)[0]
    assert md.endswith(" There the Jacobian's condition number is 1e+06–5e+08 for the joint fit in (w, c) with the centres in x, the "
                       "polish's parameters (1e+07–4e+08 with the centres in units of σ) and 2e+03–4e+04 for variable projection's in c "
                       "alone (analytic; the same in either unit of the centres); the ratio of their smallest singular values, variable "
                       "projection's over the joint fit's, is at least 1 in any common scaling of the centres for Kaufman's form of the "
                       "former, the centre columns projected off the amplitudes' (referee round 77: the joint fit's smallest singular value "
                       "is at most theirs — a Schur complement), from which the analytic Golub–Pereyra form used here differs by a term "
                       "proportional to the residual; its size depends on the units: 1–2 (median 1.5) with the centres in units of σ, "
                       "1.2–2.5 (median 1.6) in x, 49–2e+03 (median 80) with every column scaled to unit norm, which scales the two "
                       "problems' centre columns differently; TRF itself works in Coleman–Li's variables, the bounded amplitudes scaled by "
                       "v^{1/2}.")
    assert ("where variable projection stalls: its two closest nodes 0.002σ apart, its end not stationary (the largest cosine of its "
            "residual with a Jacobian column 3.0e-04), its stop 3/77") in md          # the far miss (grid +0.05)


def test_the_entries_by_the_block_they_join():
    """Referee round 60: `entry_lines` apart by regime — the first entry, entries joining an optimal block and a non-optimal one — on
    asymmetric data: the largest rate exactly (a rate 0.95 of it is not), medians (not means), beyond the outermost node at any distance
    (0.3σ is beyond) and between the nodes with the duals."""
    blk = lambda k, opt, c: dict(k=k, tau=float(k), c=np.array(c), w=np.full(k, 0.1), e=1e-3 * (1 + (5e-4 if opt else 5e-3)), e_opt=1e-3,
                                 c_opt=None, relaxed=True, dc=0.0, dc_hs=0.0, e_hs=1e-3)
    B = [dict(r=1.0, var='base', blocks=[blk(1, True, [0.0]), blk(2, False, [-0.1, 0.1]), blk(3, True, [-0.2, 0.0, 0.2]),
                                         blk(4, True, [-0.3, -0.1, 0.1, 0.3])],
              events=[(1.0, 10, 'enter', 0.0), (2.0, 11, 'enter', 0.036), (3.0, 12, 'enter', 0.05), (4.0, 13, 'enter', 0.5),
                      (5.0, 14, 'enter', 0.42)],
              entries=[dict(it=0, j=10, u0=0.0, rate=0.5, rate_max=0.5), dict(it=1, j=11, u0=0.7, rate=0.38, rate_max=0.4),
                       dict(it=2, j=12, u0=0.97, rate=0.1, rate_max=0.4), dict(it=3, j=13, u0=-0.1, rate=-0.2, rate_max=0.2),
                       dict(it=4, j=14, u0=0.9, rate=0.08, rate_max=0.4)])]
    md = '\n'.join(WM.entry_lines(B))
    assert '| none (the first entry) | 1 (1) | 1 | 1.00 (1.000–1.00) | 0.000 (0.000) | — | — |' in md
    assert '| an optimal block | 3 (3) | 1 | 0.95 (0.200–1.00) | 0.700 (0.100) | 3, 1.00σ (0.30–2.50σ) | — |' in md
    assert '| a non-optimal block | 1 (1) | 0 | 0.25 (0.250–0.25) | 0.970 (0.970) | — | 1 (0.970–0.970) |' in md
    assert '| all | 5 (5) | 2 | 0.95 (0.200–1.00) | 0.700 (0.000) | 3, 1.00σ (0.30–2.50σ) | 1 (0.970–0.970) |' in md
    B2_ = B + [dict(B[0], var='polish 2σ', entries=[dict(e, u0=e['u0'] + 5e-6) for e in B[0]['entries']]),
               dict(B[0], var='grid +0.05')]                       # a repeat (same events) and another grid (other events)
    md2 = '\n'.join(WM.entry_lines(B2_))
    assert '| an optimal block | 6 (9) |' in md2 and 'the duals of one event agree to 5.0e-06' in md2
    B3_ = B + [dict(B[0], var='polish 2σ', blocks=[B[0]['blocks'][0], blk(2, True, [-0.1, 0.1])] + B[0]['blocks'][2:])]
    md3 = '\n'.join(WM.entry_lines(B3_))                              # one event joins an optimal block in one variant only (round 65)
    assert "1 events join an optimal block in some variants and a non-optimal one in others" in md3 and '| all | 5 (10) |' in md3
    assert '| an optimal block | 4 (7) |' in md3 and 'events join an optimal block' not in md + md2
    assert WM.entry_lines([dict(B[0], entries=None)]) == []


def test_the_stop_laws_moments_are_the_projections():
    """Referee round 60: `stop_law`'s probabilities are the normal law of rms²/s² with the mean b_K² + 1 − 2K/N and the sd
    (2(N − 2K)/N² + 4b_K²/N)^{1/2}, b_K = e_K‖f*‖_N/s (recomputed here from the half-squeeze errors); and these are the exact first two
    moments of |R(f + ξ)|²/(Ns²), R the projection off a 2K-dimensional subspace and Rf of rms b_Ks (Monte Carlo, 4·10⁴ draws)."""
    from scipy.stats import norm
    ss = WM.sigma_star(WM.TAU_C)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    rf, N = float(np.sqrt(np.mean(f ** 2))), WM.X_B.size
    for s in WM.LEVELS_C:
        cdf = [0.0]
        for K in range(1, 13):
            b2 = (W.l2_error(*W.half_squeeze(K, WM.SIG, ss), WM.SIG, ss) * rf / s) ** 2
            cdf.append(norm.cdf(((1 + WM.DELTA) ** 2 - b2 - 1 + 2 * K / N) / np.sqrt(2 * (N - 2 * K) / N ** 2 + 4 * b2 / N)))
        P = WM.stop_law(s, rf, N)[1]
        assert max(abs(P[K] - (cdf[K] - cdf[K - 1])) for K in range(1, 13)) < 1e-12
    rng = np.random.default_rng(0)
    K, b, s = 4, 0.7, 1.0
    Q = np.linalg.qr(rng.standard_normal((N, 2 * K)))[0]
    bias = rng.standard_normal(N)
    bias -= Q @ (Q.T @ bias)
    bias *= b * s * np.sqrt(N) / np.linalg.norm(bias)
    xi = rng.standard_normal((40000, N)) * s
    res = bias[None, :] + xi - (xi @ Q) @ Q.T
    m = np.sum(res ** 2, axis=1) / (N * s ** 2)
    assert abs(m.mean() - (b * b + 1 - 2 * K / N)) < 4 * m.std() / np.sqrt(m.size)
    assert abs(m.std() / np.sqrt(2 * (N - 2 * K) / N ** 2 + 4 * b * b / N) - 1) < 0.02


def test_the_stationarity_test_holds_by_cauchy_schwarz_below_an_error():
    """Referee round 60: `cs_vacuous` — a block passes the path's absolute stationarity test at any state once its rms error times the
    largest of rms(φ_i) and |w_i|·rms(∂_cφ_i) is at most 10⁻⁹·max_j|⟨f*φ_j⟩| over the starting grid (shifted in 'grid +0.05')."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    o = dict(c=np.array([WM.C_B - 0.05, WM.C_B + 0.05]), w=np.array([0.5, 0.4]))
    lim = max(np.sqrt(np.mean(K.phi(o['c']) ** 2, axis=0)).max(), (o['w'] * np.sqrt(np.mean(K.dphi(o['c']) ** 2, axis=0))).max())
    for r_ in (0.5, 2.0):
        f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(r_) ** 2))
        rf = float(np.sqrt(np.mean(f ** 2)))
        for var, sh in (('base', 0.0), ('grid +0.05', 0.05)):
            e_ = 1e-9 * float(np.abs(K.phi(np.linspace(-1, 1, 21) + sh).T @ f).max()) / WM.X_B.size / (rf * lim)
            assert WM.cs_vacuous(r_, var, dict(o, e=0.999 * e_)) and not WM.cs_vacuous(r_, var, dict(o, e=1.001 * e_))


def test_the_summary_counts_the_vacuous_stationarity_tests():
    """Referee round 60: B's summary counts the relaxed non-optimal blocks whose stationarity test holds at any state."""
    blk = lambda k, e, q: dict(k=k, tau=float(k), c=WM.C_B + 0.1 * np.arange(k), w=np.full(k, 0.8 / k), e=e * q, e_opt=e, c_opt=None,
                               relaxed=True, dc=0.0, dc_hs=0.0, e_hs=e)
    B = [dict(r=1.0, var='base', blocks=[blk(1, 1e-2, 1.0001), blk(2, 1e-14, 1.5), blk(3, 1e-6, 1.5), blk(4, 1e-7, 1.5)], events=[],
              relaxed=[True] * 4, stiff=0, redone=0, failed=0, rejected=0)]
    md = WM.report(B=B)
    assert ("(0 of 3 failed the stationarity test; 1 of those that passed it pass it at any state by Cauchy–Schwarz, the test being "
            "absolute — their errors up to 1.5e-14; for the other 2 its threshold over the Cauchy–Schwarz bound") in md
    mg = sorted(WM.cs_margin(1.0, 'base', o) for o in B[0]['blocks'][2:])
    assert f"is {mg[0]:.1e}–{mg[1]:.2f} (median {np.median(mg):.1e}; above 0.1 in {sum(m > 0.1 for m in mg)}): a weak test" in md


def test_merge_b_replaces_and_orders_the_paths():
    """Round 60's `--variants`: the new paths replace the raw's of the same width and variant, the others are kept, in the order of
    VARIANTS_B, then TAUS_B."""
    v5 = WM.VARIANTS_B[5][0]
    old = np.array([dict(r=1.0, var=v5, x=0), dict(r=0.5, var='base', x=0), dict(r=0.25, var=v5, x=0)], dtype=object)
    m = WM.merge_B(old, [dict(r=2.0, var=v5, x=1), dict(r=1.0, var=v5, x=1)])
    assert [(b['r'], b['var'], b['x']) for b in m] == [(0.5, 'base', 0), (0.25, v5, 0), (1.0, v5, 1), (2.0, v5, 1)]
    assert [(b['r'], b['x']) for b in WM.merge_B(None, [dict(r=2.0, var=v5, x=1)])] == [(2.0, 1)]
    assert v5 == 'polish 2σ, 5·10⁴, repeated, no gradient test' and WM.VARIANTS_B[5][5] is None
    assert all(v[5] == 1e-15 for v in WM.VARIANTS_B[:5])


def test_b_paths_reproduce_the_raw():
    """Referee round 60: `b_job` reproduces the committed raw on B's two fastest paths at τ/σ = 2 — the base variant (rejected first
    polishes) and the repeated one (rejected repeats): the events, the blocks' kernels, errors and optima, the polish logs, the counts
    and the entries."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'width_mismatch_raw.npz'))
    for var in ('base', 'polish 2σ, 5·10⁴, repeated'):
        ref = next(b for b in raw['B'] if b['r'] == 2.0 and b['var'] == var)
        b = WM.b_job((2.0, var))
        assert [(j, k) for t, j, k, c in b['events']] == [(j, k) for t, j, k, c in ref['events']]
        assert np.allclose([e[0] for e in b['events']], [e[0] for e in ref['events']], rtol=1e-9, atol=0)
        assert [o['k'] for o in b['blocks']] == [o['k'] for o in ref['blocks']]
        for key in ('e', 'e_opt'):
            assert np.allclose([o[key] for o in b['blocks']], [o[key] for o in ref['blocks']], rtol=1e-9, atol=0)
        log = lambda x: [[(p['status'], p['nfev'], p['accepted']) for p in o['polish']] for o in x['blocks']]
        assert log(b) == log(ref) and any(p for o in b['blocks'] for p in o['polish'])
        assert [b[k] for k in ('rejected', 'stiff', 'redone', 'failed')] == [ref[k] for k in ('rejected', 'stiff', 'redone', 'failed')]
        assert [e['j'] for e in b['entries']] == [e['j'] for e in ref['entries']]
        assert max(abs(e['u0'] - f_['u0']) + abs(e['rate'] / f_['rate'] - 1) for e, f_ in zip(b['entries'], ref['entries'])) < 1e-9
    assert ref['rejected'] == 0 and any('e0' in p for o in ref['blocks'] for p in o['polish'])
    assert next(b for b in raw['B'] if b['r'] == 2.0 and b['var'] == 'base')['rejected'] == 13


def test_b_job_passes_the_variants_polish_settings(monkeypatch):
    """Referee round 65 (X01): `b_job` hands `sliding_iss` each variant's polish step, budget, repetitions and gradient test (None: the
    sixth variant), and the starting grid's shift."""
    from src.core import rich_transit as RT
    seen = []

    def stub(L, K, f, c0, rho, **kw):
        seen.append((np.asarray(c0).copy(), kw))
        return dict(hist=[(0.0, c0, np.zeros(len(c0)), 1.0)], relaxed=[], polish_log=[], events=[], stiff_redone=0, solver_redone=0,
                    solver_failed=0, polish_rejected=0, entry_info=[])
    monkeypatch.setattr(RT, 'sliding_iss', stub)
    for name, step, shift, nfev, rep_, gt in WM.VARIANTS_B:
        WM.b_job((0.25, name))
        c0, kw = seen[-1]
        assert (kw['polish_step'], kw['polish_nfev'], kw['repolish'], kw['polish_gtol']) == (step, nfev, rep_, gt)
        assert abs(c0[0] - (-1.0 + shift)) < 1e-12
    assert seen[-1][1]['polish_gtol'] is None


def test_the_polish_defaults_are_pinned():
    """Referee round 65 (R01): `polish_block`'s and `sliding_iss`'s gradient test default to 10⁻¹⁵ — other experiments call the polish
    without it."""
    import inspect
    from src.core import rich_transit as RT
    assert inspect.signature(RT.polish_block).parameters['gtol'].default == 1e-15
    assert inspect.signature(RT.sliding_iss).parameters['polish_gtol'].default == 1e-15


def test_b2_takes_dc_from_the_better_start(monkeypatch):
    """Referee round 65 (X02): `b2_job`'s dc is the block's distance from the better of the two variable-projection ends."""
    ends = iter([(np.array([0.0, 0.1]), np.array([1.0, 1.0]), 2e-6, 3, 50), (np.array([0.02, 0.13]), np.array([1.0, 1.0]), 1e-6, 2, 40)])
    monkeypatch.setattr(WM, 'vp_fit', lambda Phi_of, f, y0, max_nfev=20000: next(ends))
    q = WM.b2_job((1.0, 'base', 2, np.array([0.0, 0.1])))
    assert abs(q['dc'] - 0.03) < 1e-12 and q['vp_status'] == 3 and q['vp_nfev'] == 50 and abs(q['vp_gap'] - 0.1 / WM.SIG) < 1e-12


def test_the_vp_jacobian_and_the_condition_numbers():
    """Referee round 65 (X04, X05, X07): `vp_jacobian` is variable projection's Jacobian (against central differences); `b2_job`'s
    condition numbers are those of the joint Jacobian with the least-squares amplitudes and of the analytic one, the centres in units of
    σ, and the ratio of their smallest singular values."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    c = np.array([-0.15, 0.03, 0.2])

    def vp(cc):
        P = K.phi(cc)
        return f - P @ np.linalg.lstsq(P, f, rcond=None)[0]
    h = 1e-6
    Jd = np.stack([(vp(c + h * e) - vp(c - h * e)) / (2 * h) for e in np.eye(3)], axis=1)
    J = WM.vp_jacobian(K, f, c)
    assert np.linalg.norm(J - Jd) < 1e-6 * np.linalg.norm(J)
    q = WM.b2_job((1.0, 'base', 3, c))
    wl = np.linalg.lstsq(K.phi(c), f, rcond=None)[0]
    sj = np.linalg.svd(np.c_[K.phi(c), K.dphi(c) * wl[None, :] * WM.SIG], compute_uv=False)
    sp = np.linalg.svd(J * WM.SIG, compute_uv=False)
    assert abs(q['cond_joint'] / (sj[0] / sj[-1]) - 1) < 1e-9 and abs(q['cond_vp'] / (sp[0] / sp[-1]) - 1) < 1e-9
    assert abs(q['smin_ratio'] / (sp[-1] / sj[-1]) - 1) < 1e-9
    Jj = np.c_[K.phi(c), K.dphi(c) * wl[None, :]]
    sx, px = np.linalg.svd(Jj, compute_uv=False), np.linalg.svd(J, compute_uv=False)
    assert abs(q['cond_joint_x'] / (sx[0] / sx[-1]) - 1) < 1e-9 and abs(q['cond_vp_x'] / (px[0] / px[-1]) - 1) < 1e-9
    assert abs(q['smin_ratio_x'] / (px[-1] / sx[-1]) - 1) < 1e-9 and q['cond_joint_x'] != q['cond_joint']
    e1 = np.linalg.svd(Jj / np.linalg.norm(Jj, axis=0), compute_uv=False)
    e2 = np.linalg.svd(J / np.linalg.norm(J, axis=0), compute_uv=False)
    assert abs(q['smin_ratio_eq'] / (e2[-1] / e1[-1]) - 1) < 1e-9


def test_main_variants_merges_the_new_paths(monkeypatch, tmp_path):
    """Referee round 65 (X28–X30): `--variants i,j` runs B for those variants only (in VARIANTS_B's order) and merges them into the raw's B."""
    import exp_cohesion_leak as E
    v = [x[0] for x in WM.VARIANTS_B]
    E.save(str(tmp_path / 'width_mismatch_raw.npz'), B=[dict(r=1.0, var=v[0], x=0), dict(r=1.0, var=v[5], x=0)])
    calls, saved = [], []
    monkeypatch.setattr(WM, 'OUT', str(tmp_path))
    monkeypatch.setattr(WM, 'b_job', lambda a: (calls.append(a), dict(r=a[0], var=a[1], x=1))[1])
    monkeypatch.setattr(WM, 'report', lambda *a, **k: '')
    monkeypatch.setattr(E, 'save', lambda path, **kw: saved.append(kw))
    import multiprocessing

    class FakePool:
        def __init__(self, n):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, jobs, chunksize=1):
            return [fn(j) for j in jobs]
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'B', '--variants', '0,5'])
    WM.main()
    assert calls == [(r_, v[0]) for r_ in WM.TAUS_B] + [(r_, v[5]) for r_ in WM.TAUS_B]
    B = saved[-1]['B']
    assert len(B) == 2 * len(WM.TAUS_B) and all(b['x'] == 1 for b in B) and [b['var'] for b in B][:len(WM.TAUS_B)] == [v[0]] * len(WM.TAUS_B)


def test_b2_line_at_its_boundaries():
    """Referee round 65 (X09, X14): a miss exactly 1.01 times the optimum is far, not near."""
    B2 = [dict(r=1.0, var='base', k=10, e_block=0.505, e_rule=0.5, dc=0.0)]
    md = WM.b2_line(B2)[0]
    assert '1 ends 1.01 times above it' in md and 'within 1% of it' not in md


def test_the_entries_at_their_boundaries():
    """Referee round 65 (X23, X24): an entrant at the block's outermost node is between its nodes (d = 0), and a rate within a relative
    10⁻¹² of the largest counts as the largest."""
    blk = lambda k, c: dict(k=k, tau=float(k), c=np.array(c), w=np.full(k, 0.1), e=1e-3, e_opt=1e-3, c_opt=None, relaxed=True, dc=0.0,
                            dc_hs=0.0, e_hs=1e-3)
    B = [dict(r=1.0, var='base', blocks=[blk(1, [0.0]), blk(2, [-0.1, 0.1])],
              events=[(1.0, 10, 'enter', 0.0), (2.0, 11, 'enter', 0.0), (3.0, 12, 'enter', 0.1)],
              entries=[dict(it=0, j=10, u0=0.0, rate=0.5, rate_max=0.5), dict(it=1, j=11, u0=0.9, rate=0.4 * (1 - 5e-13), rate_max=0.4),
                       dict(it=2, j=12, u0=0.95, rate=0.1, rate_max=0.4)])]
    md = '\n'.join(WM.entry_lines(B))
    assert '| an optimal block | 2 (2) | 1 |' in md and '| — | 2 (0.900–0.950) |' in md


def test_frame_series_against_the_exact_defects():
    """Referee round 81 (D3cb-5): `gh_nodes` is the Gauss–Hermite rule of N(0, 1) (He_K at its nodes zero, the weights summing to 1, to
    the arithmetic's 110 digits); `frame_series` at rational ρ and v is Σ_{n ≥ 2K} (ρⁿ/n!)D_n(v)² with `gh_defect_v`'s exact defects."""
    from decimal import Decimal, localcontext
    from fractions import Fraction
    from math import factorial
    for K in (1, 3, 6):
        zs, ws = WM.gh_nodes(K)
        with localcontext() as ctx:
            ctx.prec = WM.A5_PREC
            assert abs(sum(ws) - 1) < Decimal(10) ** -100 and all(abs(WM.he_all(z, K)[K]) < Decimal(10) ** -90 for z in zs)
    for K, rho, v in ((1, Fraction(1, 3), Fraction(3, 2)), (2, Fraction(1, 4), Fraction(2)), (3, Fraction(2, 5), Fraction(5, 4))):
        with localcontext() as ctx:
            ctx.prec = WM.A5_PREC
            rd, vd = Decimal(rho.numerator) / rho.denominator, Decimal(v.numerator) / v.denominator
        S, terms, low = WM.frame_series(K, rd, vd)
        ref = [rho ** n / factorial(n) * WM.gh_defect_v(n, K, v) ** 2 for n in range(2 * K, 2 * K + 2 * len(terms), 2)]
        assert all(abs(Fraction(t) / r_ - 1) < Fraction(1, 10 ** 60) for t, r_ in zip(terms, ref) if r_ != 0)
        assert low < Decimal(10) ** -80


def test_the_frames_and_the_first_tail_ratio():
    """Referee round 65: `gh_defect_v` at v = 2 is `gh_defect`, and the first tail ratio is −vK(K + 2) + (v − 1)(K + 1)(2K + 1) exactly;
    the frame at ρ = t is the half-squeeze rule's (v_ρ = τ²σ/σ*, v′ = 2), at ρ = λ v′ = 1, at ρ = q the variance 2τ²/(β² + 1)."""
    from fractions import Fraction
    for K in range(1, 8):
        for n in range(2 * K, 2 * K + 6, 2):
            assert WM.gh_defect_v(n, K, 2) == W.gh_defect(n, K)
        for v in (Fraction(1, 2), Fraction(5, 3)):
            assert WM.gh_defect_v(2 * K + 2, K, v) / WM.gh_defect_v(2 * K, K, v) == -v * K * (K + 2) + (v - 1) * (K + 1) * (2 * K + 1)
    for r_ in (0.25, 1.0, 2.0):
        ss = WM.sigma_star(r_)
        P = W.params(WM.SIG, ss)
        tau2 = P['tau'] ** 2
        v, vp = WM.frame(P['t'], WM.SIG, tau2)
        assert abs(v / (tau2 * WM.SIG / ss) - 1) < 1e-12 and abs(vp - 2) < 1e-12
        assert abs(WM.frame(P['lam'], WM.SIG, tau2)[1] - 1) < 1e-12
        b2 = np.sqrt(1 + 2 * tau2 / WM.SIG ** 2)
        assert abs(WM.frame((b2 - 1) / (b2 + 1), WM.SIG, tau2)[0] / (2 * tau2 / (b2 + 1)) - 1) < 1e-12


def test_b3_restarts_the_polish_from_its_recorded_start(monkeypatch):
    """Referee rounds 65, 70: `b3_job` finds the path's first non-optimal block and restarts its first polish from that polish's recorded
    start, one setting at a time — (the path's budget, its gradient test), (B3_BIG, the test), (the path's budget, none), (B3_BIG, none) —,
    with a variant that has a gradient test; `seed_rule` judges each end from the start; `path_polishes` records the calls in the path's
    order and restores `polish_block`."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    rf = np.sqrt(np.mean(f ** 2))
    n = 4
    cs = np.linspace(-0.3, 0.3, n)
    hist = [(0.0, cs.copy(), np.zeros(n), 1.0)]
    c1 = cs.copy(); c1[1] = 0.03; w1 = np.zeros(n); w1[1] = 0.8
    r1 = np.sqrt(np.mean((f - K.phi(np.array([0.03])) @ np.array([0.8])) ** 2))
    hist.append((1.0, c1, w1, r1))                                     # block 0: one kernel, not the optimum
    R = dict(hist=hist, relaxed=[True], polish_log=[[dict(status=0, nfev=1, accepted=True)]])
    start = dict(w=np.array([0.7]), c=np.array([0.0]), s=np.array([1.0]), max_nfev=5000, gtol=1e-15)
    monkeypatch.setattr(WM, 'path_polishes', lambda r_, var: (R, K, f, [start]))
    calls = []

    def stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        calls.append((max_nfev, gtol, tuple(w_), tuple(c_)))
        return w_ * 1.1, c_ + 0.01, dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', stub)
    v = WM.VARIANTS_B[3]
    q = WM.b3_job((1.0, v[0]))
    assert [x[:2] for x in calls] == [(v[3], v[5]), (WM.B3_BIG, v[5]), (v[3], None), (WM.B3_BIG, None)]
    assert all(x[2:] == ((0.7,), (0.0,)) for x in calls)
    e_opt = min(WM.sample_fit(K.phi, f, c1[[1]])[2], WM.sample_fit(K.phi, f, WM.C_B + W.half_squeeze(1, WM.SIG, WM.sigma_star(1.0))[0])[2]) / rf
    assert np.isclose(q['e_opt'], e_opt) and np.isclose(q['e_block'], r1 / rf / e_opt)
    assert np.isclose(q['e_start'], np.sqrt(np.mean((f - K.phi(np.array([0.0])) @ np.array([0.7])) ** 2)) / rf / e_opt)
    sr = WM.seed_rule(K, f, np.array([0.7]), np.array([0.0]), np.array([1.0]), np.array([0.01]), v[1])
    assert q['path']['move'] == sr['move'] and q['path']['accepted'] == sr['accepted']


def test_path_polishes_records_and_restores(monkeypatch):
    """Referee round 70: `path_polishes` runs the variant's path with `polish_block` wrapped — every call's start and settings recorded
    in order, the real solve's result passed through — and restores `polish_block` even when the path fails."""
    from src.core import rich_transit as RT
    real, seen = RT.polish_block, {}

    def fake_path(L, K, f, c0, rho, **kw):
        seen.update(kw)
        out = [RT.polish_block(K, f, np.array([0.5]), np.array([0.01 * i]), np.array([1.0]), max_nfev=7 + i, info=True, gtol=None)
               for i in range(2)]
        return dict(out=out)
    monkeypatch.setattr(RT, 'sliding_iss', fake_path)
    v = WM.VARIANTS_B[4]
    R, K, f, calls = WM.path_polishes(1.0, v[0])
    assert RT.polish_block is real and len(calls) == 2
    assert [(c['max_nfev'], c['gtol'], float(c['c'][0])) for c in calls] == [(7, None, 0.0), (8, None, 0.01)]
    assert seen['polish_step'] == v[1] and seen['polish_nfev'] == v[3] and seen['repolish'] == v[4] and seen['polish_gtol'] == v[5]
    assert len(R['out'][0]) == 3 and R['out'][0][2]['nfev'] <= 7

    def boom(*a, **k):
        raise RuntimeError('path failed')
    monkeypatch.setattr(RT, 'sliding_iss', boom)
    try:
        WM.path_polishes(1.0, v[0])
    except RuntimeError:
        pass
    assert RT.polish_block is real


def test_b3_binds_and_lines_on_synthetic_results():
    """Referee rounds 70, 77: `b3_binds` labels from the ends that reach the optimum and are accepted — the budget, the gradient test,
    either, both —, "(and the radius)" where every reaching end is rejected, the radius alone, "neither" with the stop, the
    inconsistency —, judged against its ref; B3's table and summary, below B4_MAX against B4's long-double end."""
    P = lambda ra, acc=True, st=0, nf=5000, mv=0.5: dict(status=st, nfev=nf, ratio=ra, move=mv, accepted=acc)
    mk = lambda r_, var, e_opt, a, b, c, d: dict(r=r_, var=var, k=6, e_opt=e_opt, e_start=17.0, e_block=1.46, n_polish=1,
                                                 **dict(zip(WM.B3_KEYS, (a, b, c, d))))
    B3 = [mk(0.5, 'polish 2σ', 6e-8, P(1.46), P(1.0, st=1, nf=6635), P(1.46), P(1.0, st=3, nf=6654)),
          mk(0.25, 'base', 6e-8, P(1.002, st=1), P(1.0017, st=1), P(1.0000007, st=3), P(1.0, st=3)),
          mk(1.0, 'grid +0.05', 3e-7, P(1.0, acc=False, mv=1.98), P(1.0, acc=False, mv=1.98), P(1.0, acc=False), P(1.0, acc=False)),
          mk(0.75, 'x', 2e-9, P(1.3), P(1.2), P(1.2), P(1.1)),
          mk(1.5, 'y', 4e-9, P(1.0), P(1.0), P(1.0), P(1.0)),
          mk(0.25, 'z', 2.5e-11, P(3.7), P(3.7), P(3.7), P(1.0147))]
    labs = [WM.b3_binds(q) for q in B3]
    assert labs == ['the budget', 'the gradient test', 'the radius', 'neither (the budget spent)',
                    "nothing (inconsistent: the path's own polish reaches it, accepted)", 'neither (the budget spent)']
    alt = lambda i, **kw: dict(B3[i], **kw)
    assert WM.b3_binds(alt(0, budget=P(1.0, acc=False))) == 'the budget and the gradient test'
    assert WM.b3_binds(alt(0, **{'budget, no test': P(1.0, acc=False)})) == 'the budget'
    assert WM.b3_binds(alt(1, budget=P(1.0, acc=False))) == 'the gradient test'
    assert WM.b3_binds(alt(1, budget=P(1.0))) == 'the budget or the gradient test'
    assert WM.b3_binds(alt(1, budget=P(1.0, acc=False), **{'no test': P(1.0, acc=False), 'budget, no test': P(1.0, acc=False)})) == \
        'the budget or the gradient test (and the radius)'
    assert WM.b3_binds(alt(0, budget=P(1.2), **{'budget, no test': P(1.0, acc=False)})) == 'the budget and the gradient test (and the radius)'
    assert WM.b3_binds(alt(3, **{'budget, no test': P(1.1, st=3)})) == "neither (the solve's stop, status 3)"
    assert WM.b3_binds(alt(3, **{'budget, no test': P(1.1, st=1)})) == "neither (the solve's stop, status 1)"
    near = alt(5, **{'budget, no test': P(1.0015)})
    assert WM.b3_binds(near) == 'neither (the budget spent)' and WM.b3_binds(near, ref=1.001) == 'the budget and the gradient test'
    s = '\n'.join(WM.b3_lines(B3 + [dict(r=2.0, var='base', k=None)]))
    assert ("| 0.5 | polish 2σ | 6 | 6.00e-08 | 17 | 1.46 | 0/5000, 1.46000, 0.50σ a | 1/6635, 1.00000, 0.50σ a | 0/5000, 1.46000, 0.50σ a | "
            "3/6654, 1.00000, 0.50σ a | the budget |") in s
    assert "| 0.25 | z | 6 | 2.50e-11† |" in s and s.endswith('\n')
    assert ("\n\nBelow it (1 blocks) the restart with both changed stops at 1.0147 (0/5000: neither (the budget spent)) times the stored "
            "optimum.") in s
    assert "Above 1e-10 (5 blocks: τ/σ (variant)) the end is bounded by" in s
    assert ("- neither (the budget spent): 1 — 0.75 (x);\n- nothing (inconsistent: the path's own polish reaches it, accepted): 1 — 1.5 (y);"
            "\n- the budget: 1 — 0.5 (polish 2σ);\n- the gradient test: 1 — 0.25 (base);\n- the radius: 1 — 1 (grid +0.05).") in s
    assert ("by the optimum (one solver from one kind of start; above B4_MAX a negative value is the end below the stored optimum, "
            "converged only to that): 3e-07: +0.0e+00; 6e-08: +0.0e+00; 4e-09: +0.0e+00; 2e-09: +1.0e-01; 2.5e-11: +1.5e-02.") in s
    B4 = [dict(r=0.25, k=6, start='block', e_opt=2.5e-11, **{'long double': dict(ratio=0.99)}),
          dict(r=0.25, k=6, start='optimum', e_opt=2.5e-11, **{'long double': dict(ratio=0.995)})]
    s4 = '\n'.join(WM.b3_lines([B3[5], alt(5, r=0.75, var='w', **{'budget, no test': P(0.99)})], B4))
    assert "stops at 1.0147 (0/5000: neither (the budget spent)), 0.9900 (0/5000: the budget and the gradient test) times the stored optimum." in s4
    s4 = '\n'.join(WM.b3_lines([B3[5]], B4))
    assert "stops at 1.0147 (0/5000: neither (the budget spent)) times the stored optimum, 1.0249 times B4's best end." in s4
    assert "2.5e-11: +2.5e-02." in s4
    z2 = alt(5, var='z2', e_opt=2.6e-11, **{'budget, no test': P(1.0)})  # referee round 81 (D3cb-6): the reference over each row's own optimum
    s4 = '\n'.join(WM.b3_lines([B3[5], z2], B4))
    assert f"1.0249, {1.0 / (0.99 * 2.5e-11 / 2.6e-11):.4f} times B4's best end." in s4
    Bs = [dict(r=0.25, var='z', blocks=[dict(k=6, e=1.1, e_opt=1.0, polish=[dict(status=0, nfev=5000)])]),
          dict(r=0.25, var='z2', blocks=[dict(k=6, e=1.1, e_opt=1.0, polish=[dict(status=1, nfev=77)])])]
    s5 = '\n'.join(WM.b3_lines([B3[5], z2], B4, Bs))
    assert "The 'path' restart reproduces the first polish B stored — its stop and evaluations — at 1 of 2 blocks (not at 0.25 (z2))." in s5
    q2 = [dict(B3[1]), dict(B3[0])]
    q2[0]['budget, no test'] = P(0.99998, st=3)
    s2 = '\n'.join(WM.b3_lines(q2))
    assert "6e-08: -2.0e-05 … +0.0e+00." in s2
    assert WM.b3_lines([]) == [] and WM.b3_lines(None) == []


def test_seed_rule_with_two_seed_kernels_and_its_default():
    """Referee round 70 (S02, S05–S08): with two kernels near their seeds the rule takes the larger of their moves, each the nearer of its
    start and its ascent target — one decided by the target —; the seed set is |w| ≤ 10⁻²·max|w|, `sliding_iss`'s default."""
    import inspect
    from src.core import rich_transit as RT
    assert (inspect.signature(WM.seed_rule).parameters['wfrac'].default
            == inspect.signature(RT.sliding_iss).parameters['polish_wfrac'].default == 1e-2)
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    c, w = np.array([-0.4, 0.03, 0.5]), np.array([5e-3, 0.8, 7e-3])
    r0 = f - K.phi(c) @ w
    t0, t2 = (RT.ascent_target(K, r0, c[i], 1.0)[0] for i in (0, 2))
    cp = c.copy(); cp[0] = t0 + 0.2 * WM.SIG; cp[2] = c[2] + 0.7 * WM.SIG
    r = WM.seed_rule(K, f, w, c, np.ones(3), cp, 1.0)
    m0 = min(abs(cp[0] - c[0]), abs(cp[0] - t0)) / WM.SIG
    m2 = min(abs(cp[2] - c[2]), abs(cp[2] - t2)) / WM.SIG
    assert abs(cp[0] - c[0]) / WM.SIG > 1.0 and np.isclose(m0, 0.2)
    assert np.isclose(r['move'], max(m0, m2)) and r['accepted']
    assert not WM.seed_rule(K, f, w, c, np.ones(3), cp, 0.1)['accepted']
    assert WM.seed_rule(K, f, np.array([5e-3, 0.8, 0.5]), c, np.ones(3), cp, 1.0)['move'] == m0


def test_part_a5_on_its_own_numbers():
    """Referee rounds 70, 81 (A01–A05, A09; D3cb-5): A5's frames at ρ = t, q, λ and the Mercer frame's rule N(0, 2τ²/(β² + 1)) — not the
    Mercer-zero rule — with its errors at the given K, by hand; the tail-ratio checks all exact; the exact frame series (`a5_series`) is
    the rule's error (its ratio per kernel the actual errors' to 10⁻⁴), its tail above the first term from K = 6 at τ/σ = 2, and its
    lines."""
    A = dict(rows=[dict(r=1.0, K=3), dict(r=1.0, K=4), dict(r=2.0, K=5)])
    A5 = WM.part_A5(A)
    x = next(r for r in A5['rows'] if r['r'] == 1.0)
    ss = WM.sigma_star(1.0)
    P = W.params(WM.SIG, ss)
    tau2 = P['tau'] ** 2
    b2 = np.sqrt(1 + 2 * tau2 / WM.SIG ** 2)
    vhs, vM = tau2 * WM.SIG / ss, 2 * tau2 / (b2 + 1)
    assert x['Ks'] == [3, 4] and np.isclose(x['vM'], vM / vhs) and np.isclose(x['q'], (b2 - 1) / (b2 + 1))
    assert np.isclose(x['frames']['t']['v'], 1.0) and np.isclose(x['frames']['λ']['vp'], 1.0)
    assert np.allclose(x['e_mframe'], [W.l2_error(*W.scaled_rule(K, WM.SIG, ss, vM), WM.SIG, ss) for K in (3, 4)])
    assert all(c_['ok'] for c_ in A5['check'])
    sr = {y['r']: y for y in A5['series']}
    by = {y['K']: y for y in sr[1.0]['frames']['Mercer']}
    rat = lambda K: float(np.sqrt(by[K - 1]['S'] / by[K]['S']))
    e = [W.l2_error(*W.scaled_rule(K, WM.SIG, ss, vM), WM.SIG, ss) for K in (3, 4)]
    assert abs(rat(4) / (e[0] / e[1]) - 1) < 1e-4                       # referee round 81: the exact series is the rule's error
    assert sr[2.0]['frames']['Mercer'][4]['tail'] < 1 < sr[2.0]['frames']['Mercer'][5]['tail']   # the tail above the first from K = 6
    A5['rows'] = [dict(x, Ks=[3, 4, 5, 6], e_mframe=[1.0, 0.25, 0.05, 0.01]) for x in A5['rows']]
    md = '\n'.join(WM.a5_lines(A5))
    assert f"| 4.0000, 5.0000, 5.0000 at K = 4–6 ({1 / x['t']:.4g}, {1 / x['q']:.4g}) |" in md and "ratio per kernel at the last three K" in md
    assert "2/(q(2 − v′))" not in md and "| 2 | Mercer | " in md and "| 1 | half-squeeze | " in md
    assert (f"| 1 | Mercer | " + ', '.join(f"{by[K]['tail']:.3g}" for K in (5, 10, 20, 40)) + " | 8 | "
            + ', '.join(f"{rat(K):.4f}" for K in (10, 20, 30, 40))) in md


def test_b2_job_where_variable_projection_leaves_the_block():
    """Referee round 70 (B04–B06, B09, V02, V03): `b2_job` from a two-kernel block away from the optimum — variable projection from the
    block and from the half-squeeze nodes, the better end's distance, the stop, the closest nodes and the stationarity cosine — against a
    direct `vp_fit`; `vp_fit` returns sorted centres, their least-squares coefficients, the rms, the status and the evaluations."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    r_ = 1.0
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(r_) ** 2))
    rf = np.sqrt(np.mean(f ** 2))
    c = np.array([0.2, -0.15])
    cv, av, ev, st, nf = WM.vp_fit(K.phi, f, c)
    assert np.all(np.diff(cv) > 0) and np.allclose(av, np.linalg.lstsq(K.phi(cv), f, rcond=None)[0])
    assert np.isclose(ev, np.sqrt(np.mean((f - K.phi(cv) @ av) ** 2))) and isinstance(st, int) and isinstance(nf, int) and nf > 0
    q = WM.b2_job((r_, 'base', 2, c))
    yh, _ = W.half_squeeze(2, WM.SIG, WM.sigma_star(r_))
    cr, _, er, _, _ = WM.vp_fit(K.phi, f, WM.C_B + yh)
    assert np.isclose(q['e_block'], ev / rf) and np.isclose(q['e_rule'], er / rf) and q['vp_status'] == st and q['vp_nfev'] == nf
    assert np.isclose(q['dc'], np.abs(np.sort(c) - (cv if ev <= er else cr)).max()) and q['dc'] > 0.01
    assert np.isclose(q['vp_gap'], np.diff(cv).min() / WM.SIG)
    Je = WM.vp_jacobian(K, f, cv)
    re = f - K.phi(cv) @ av
    assert np.isclose(q['vp_cos'], np.max(np.abs(Je.T @ re) / (np.linalg.norm(Je, axis=0) * np.linalg.norm(re))))


def test_b3_jobs_and_run_b_with_variants(monkeypatch):
    """Referee round 70 (P01, X02): B3's jobs are the paths with a non-optimal block; `run_B` with variants reruns those and merges them
    into the stored paths (none dropped), without variants reruns all."""
    o = lambda e: dict(e=e, e_opt=1.0)
    B = [dict(r=0.5, var='base', blocks=[o(1.0), o(1.0005)]), dict(r=1.0, var='base', blocks=[o(1.0), o(1.2)]),
         dict(r=0.5, var='polish 2σ', blocks=[o(2.0)])]
    assert WM.b3_jobs(B) == [(1.0, 'base'), (0.5, 'polish 2σ')]
    new = [dict(r=0.5, var='polish 2σ', blocks=[o(1.0)])]
    monkeypatch.setattr(WM, 'part_B', lambda idx=None: new if idx is not None else ['all'])
    got = WM.run_B(B, {1})
    assert len(got) == 3 and next(b for b in got if b['var'] == 'polish 2σ')['blocks'] == [o(1.0)]
    assert [(b['var'], b['r']) for b in got] == [('base', 0.5), ('base', 1.0), ('polish 2σ', 0.5)]
    assert WM.run_B(B, None) == ['all']


def test_b3_stops_the_path_at_the_stored_blocks_first_polish(monkeypatch):
    """Referee rounds 70, 77 (B3's cost; FN01, BR08): `first_nonoptimal` finds in B's stored blocks the first polished non-optimal block,
    the number of polish calls before it (unequal counts, a block without polishes) and its first polish's stored stop; `path_polishes`
    with stop_after = n records call n and stops the path before its solve (R None, the solve not called, `polish_block` restored);
    `b3_job` with the block restarts that recorded start, with the stored block's errors over its optimum's, and refuses a rerun whose call
    n is not the block's size or whose 'path' restart is not the stored first polish."""
    from src.core import rich_transit as RT
    o = lambda k, e, n: dict(k=k, e=e, e_opt=0.5, polish=[dict(status=1, nfev=1, accepted=True)] * n)
    b = dict(r=0.5, var='base', blocks=[o(1, 0.5, 3), o(2, 1.5, 0), o(3, 0.75, 3), o(4, 1.0, 1)])
    assert WM.first_nonoptimal(b) == dict(n0=3, k=3, e=0.75, e_opt=0.5, n_polish=3, first=dict(status=1, nfev=1, accepted=True))
    assert WM.first_nonoptimal(dict(blocks=[o(1, 0.5, 2), o(2, 0.50025, 1)])) is None
    assert WM.b3_jobs([b, dict(r=1.0, var='base', blocks=[o(1, 0.5, 1)])]) == [(0.5, 'base', WM.first_nonoptimal(b))]
    solved = []

    def real_stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        solved.append(c_.size)
        return w_, c_, dict(status=1, nfev=1)
    monkeypatch.setattr(RT, 'polish_block', real_stub)

    def fake_path(L, K, f, c0, rho, **kw):
        for k in (1, 1, 1, 3, 3, 4):
            RT.polish_block(K, f, np.full(k, 0.5), np.linspace(-0.1, 0.1, k), np.ones(k), info=True)
        return 'full'
    monkeypatch.setattr(RT, 'sliding_iss', fake_path)
    R, K, f, calls = WM.path_polishes(0.5, 'base', stop_after=3)
    assert R is None and len(calls) == 4 and calls[3]['c'].size == 3 and solved == [1, 1, 1] and RT.polish_block is real_stub
    R, K, f, calls = WM.path_polishes(0.5, 'base')
    assert R == 'full' and len(calls) == 6 and solved == [1, 1, 1, 1, 1, 1, 3, 3, 4]
    q = WM.b3_job((0.5, 'base', WM.first_nonoptimal(b)))
    assert q['k'] == 3 and q['e_opt'] == 0.5 and q['e_block'] == 1.5 and q['n_polish'] == 3
    assert q['settings'] == dict(nfev=5000, gtol=1e-15) and all(q[k]['nfev'] == 1 for k in WM.B3_KEYS)
    rf = np.sqrt(np.mean(f ** 2))
    st = calls[3]
    assert np.isclose(q['e_start'], np.sqrt(np.mean((f - K.phi(st['c']) @ st['w']) ** 2)) / rf / 0.5)
    assert np.isclose(q['path']['ratio'], np.sqrt(np.mean((f - K.phi(st['c']) @ st['w']) ** 2)) / rf / 0.5)
    for bad, msg in ((dict(WM.first_nonoptimal(b), k=4), 'not the stored block'),
                     (dict(WM.first_nonoptimal(b), first=dict(status=0, nfev=1)), 'is not the stored first polish'),
                     (dict(WM.first_nonoptimal(b), first=dict(status=1, nfev=2)), 'is not the stored first polish')):
        try:
            WM.b3_job((0.5, 'base', bad))
            raise AssertionError('a mismatched block was accepted')
        except RuntimeError as e:
            assert msg in str(e)
    assert WM.b3_job((0.5, 'base', dict(WM.first_nonoptimal(b), first=dict(status=1))))['k'] == 3
    solved.clear()                                                      # referee round 81: acceptance compared, before the larger restarts
    try:
        WM.b3_job((0.5, 'base', dict(WM.first_nonoptimal(b), first=dict(status=1, nfev=1, accepted=False))))
        raise AssertionError('a mismatched acceptance was accepted')
    except RuntimeError as e:
        assert 'is not the stored first polish' in str(e) and solved == [1, 1, 1, 3]


def test_lm_fit_in_double_and_long_double():
    """Referee rounds 70, 77 (B4): `lm_fit` from variable projection's end at τ/σ = 1/2, K = 8 (its optimum ≈ 2.1·10⁻¹⁰): in double the
    solve stalls with nothing gained (its steps lost to the residual's rounding), in long double the same solve goes on 10⁻⁵ below, and so
    does the solve in double with the residual about an anchor; the iteration cap stops it; the anchored solve's marks, and (referee
    round 81) the anchored solve in long double from the anchored double's end, ending no higher; on an exact two-kernel sum from a
    perturbed start both reach the sum; `lm_rms` is the long-double rms; the defaults are B4's."""
    import inspect
    import pytest
    from src.core import rich_transit as RT
    if np.finfo(np.longdouble).eps >= 1e-18:
        pytest.skip('long double is not extended on this platform')
    K = RT.Bumps(WM.X_B, WM.SIG)
    ss = WM.sigma_star(0.5)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    co, ao, eo = WM.sample_fit(K.phi, f, WM.C_B + W.half_squeeze(8, WM.SIG, ss)[0])
    pd, itd, accd, std = WM.lm_fit(f, ao, co, np.float64)
    pl, itl, accl, stl = WM.lm_fit(f, ao, co, np.longdouble)
    ed, el = WM.lm_rms(f, pd, 8) / eo, WM.lm_rms(f, pl, 8) / eo
    assert pd.dtype == np.float64 and pl.dtype == np.longdouble and std == stl == 'stall' and accd <= 5 and accl > 100
    assert ed > 1 - 1e-6 and el < ed - 5e-6
    pa, ita, acca, sta = WM.lm_fit(f, ao, co, np.float64, anchored=True, max_iter=1000)
    ea = WM.lm_rms(f, pa, 8) / eo
    assert pa.dtype == np.float64 and ita == 1000 and sta == 'cap' and acca > 100 and abs(ea / el - 1) < 1e-8
    marks = {10: None, 10 ** 9: None}
    WM.lm_fit(f, ao, co, np.float64, anchored=True, max_iter=50, marks=marks)
    assert np.array_equal(marks[10], WM.lm_fit(f, ao, co, np.float64, anchored=True, max_iter=10)[0]) and marks[10 ** 9] is None
    pq, itq, accq, stq = WM.lm_fit(f, pa[:8], pa[8:], np.longdouble, anchored=True, max_iter=300)
    assert pq.dtype == np.longdouble and WM.exact_rms(f, pq, 8) <= WM.exact_rms(f, pa, 8) * (1 + 1e-12)
    p, it, acc, st = WM.lm_fit(f, ao, co, np.longdouble, max_iter=50)
    assert it == 50 and st == 'cap' and WM.lm_rms(f, p, 8) / eo > el
    w, c = np.array([0.5, 0.3]), np.array([-0.05, 0.2])
    f2 = K.phi(c) @ w
    for kw in (dict(), dict(anchored=True)):
        p2, it2, acc2, st2 = WM.lm_fit(f2, w * np.array([1.02, 0.97]), c + np.array([0.01, -0.008]), np.float64, **kw)
        assert st2 == 'stall' and acc2 > 3 and np.allclose(p2, np.r_[w, c], rtol=0, atol=1e-12) and WM.lm_rms(f2, p2, 2) < 1e-15
    xl = WM.X_B.astype(np.longdouble)
    pr = np.r_[w, c].astype(np.longdouble)
    rl = (f2 + 1e-3).astype(np.longdouble) - np.exp(-(xl[:, None] - pr[None, 2:]) ** 2 / (2 * np.longdouble(WM.SIG) ** 2)) @ pr[:2]
    assert WM.lm_rms(f2 + 1e-3, np.r_[w, c], 2) == float(np.sqrt(np.mean(rl * rl)))
    sig = inspect.signature(WM.lm_fit).parameters
    assert sig['max_iter'].default == WM.B4_CAP == 400_000 and sig['mu0'].default == 1e-20 and WM.B4_MAX == 2e-9
    assert sig['anchored'].default is False and sig['anchor_every'].default == WM.B4_ANCHOR == 2000 and WM.B4_FLAT == 1e-10
    assert sig['stall_reanchor'].default is True and WM.EXACT_PREC == 60
    assert WM.B4_MARKS == (100_000, 200_000, 300_000) and sig['marks'].default is None
    assert [k for k, _ in WM.B4_PREC] == ['double', 'long double'] and WM.B4_PREC[1][1] is np.longdouble


def test_b4_jobs_job_and_lines(monkeypatch):
    """Referee rounds 70, 77, 81 (B4; BS06, BS07; D3cb-2): the jobs are B3's distinct blocks below B4_MAX (its first row for each — not the
    first path of that width —, the block of that K nearest its optimum), from the block and from the stored optimum with least-squares
    amplitudes; `b4_job` (its solves faked, each call's values distinct) reports the start's and each end's exact and long-double error
    over the optimum's, for both arithmetics, the anchored double with its marks and its predicted floor, without stall re-anchoring, the
    long double's anchored solve from the anchored double's end, and the polish's TRF without the gradient test; `b4_lines` prints the
    table against each block's best end and the summaries."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    o = lambda k, e_opt, c: dict(k=k, e=2 * e_opt, e_opt=e_opt, c=np.array(c), w=np.full(len(c), 0.4), c_opt=np.array(c) + 0.01)
    B = [dict(r=0.5, var='v1', blocks=[o(1, 1e-3, [0.0]), o(2, 3e-10, [-0.1, 0.1]), o(2, 5e-10, [-0.2, 0.2])]),
         dict(r=0.5, var='v2', blocks=[o(2, 3e-10, [-0.15, 0.15])]),
         dict(r=1.0, var='v1', blocks=[o(3, 1.5e-9, [-0.1, 0.0, 0.1])]),
         dict(r=1.0, var='v2', blocks=[o(3, 1e-8, [-0.12, 0.0, 0.12])])]
    q = lambda r_, var, k, e_opt: dict(r=r_, var=var, k=k, e_opt=e_opt)
    B3 = [q(0.5, 'v2', 2, 3e-10), q(0.5, 'v1', 2, 5e-10), q(1.0, 'v2', 3, 1e-8), q(1.0, 'v1', 3, 1.5e-9), dict(r=2.0, var='v1', k=None)]
    jobs = WM.b4_jobs(B, B3)
    assert [(j[0], j[1], j[2], j[3]) for j in jobs] == [(0.5, 'v2', 2, 'block'), (0.5, 'v2', 2, 'optimum'), (1.0, 'v1', 3, 'block'),
                                                         (1.0, 'v1', 3, 'optimum')]
    assert np.array_equal(jobs[0][5], [-0.15, 0.15]) and np.allclose(jobs[1][5], [-0.14, 0.16], rtol=0, atol=1e-15) and jobs[0][6] == 3e-10
    assert np.array_equal(jobs[2][5], [-0.1, 0.0, 0.1]) and jobs[2][6] == 1.5e-9
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.5) ** 2))
    assert np.allclose(jobs[1][4], np.linalg.lstsq(K.phi(np.array([-0.14, 0.16])), f, rcond=None)[0])
    assert WM.b4_jobs(B, None) == [] and WM.b4_jobs(B, [q(0.5, 'v1', 2, 2e-9)]) == []
    seen, pol = [], []
    calls = iter(range(10))

    def fake(f_, w_, c_, dtype, max_iter=WM.B4_CAP, mu0=1e-20, anchored=False, anchor_every=WM.B4_ANCHOR, marks=None, stall_reanchor=True):
        i = next(calls)                                                 # referee round 81: distinct values per call, moved amplitudes
        seen.append((dtype, anchored, marks is not None, stall_reanchor, np.r_[w_, c_].astype(float)))
        if marks is not None:
            marks[WM.B4_MARKS[0]] = np.r_[w_, c_ + 0.002]
        return np.r_[w_ * (1 + 0.001 * i), c_ + 0.001 * (i + 1)].astype(dtype), 7 + i, 5 + i, ('stall', 'cap', 'flat')[i % 3]
    monkeypatch.setattr(WM, 'lm_fit', fake)

    def pb(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        pol.append((max_nfev, gtol, info))
        return w_ * 1.01, c_ + 0.004, dict(status=3, nfev=1234)
    monkeypatch.setattr(RT, 'polish_block', pb)
    rf = float(np.sqrt(np.mean(f ** 2)))
    w0, c0 = jobs[0][4], jobs[0][5]
    ex = lambda w_, c_: WM.exact_rms(f, np.r_[w_, c_], 2) / rf / 3e-10
    out = WM.b4_job(jobs[0])
    assert [x[:4] for x in seen] == [(np.float64, False, False, True), (np.longdouble, False, False, True), (np.float64, True, True, True),
                                     (np.float64, True, False, False), (np.longdouble, True, False, True)] and pol == [(WM.B3_BIG, None, True)]
    pa = np.r_[w0 * 1.002, c0 + 0.003]                                  # the anchored double's end (call 2) starts the long double's
    assert np.allclose(seen[4][4], pa, rtol=0, atol=1e-15) and all(np.array_equal(x[4], np.r_[w0, c0]) for x in seen[:4])
    assert out['start'] == 'block' and out['k'] == 2 and out['e_start'] == ex(w0, c0)
    for i, k_ in enumerate(('double', 'long double', 'anchored', 'nostall')):
        q = out[k_]
        assert (q['it'], q['acc'], q['stop']) == (7 + i, 5 + i, ('stall', 'cap', 'flat')[i % 3]), k_
        assert np.isclose(q['ratio'], ex(w0 * (1 + 0.001 * i), c0 + 0.001 * (i + 1)), rtol=1e-12), k_
        assert np.isclose(q['ratio_ld'], WM.lm_rms(f, np.r_[w0 * (1 + 0.001 * i), c0 + 0.001 * (i + 1)], 2) / rf / 3e-10, rtol=1e-12)
    assert np.isclose(out['anchored ld']['ratio'], ex(pa[:2] * 1.004, pa[2:] + 0.005), rtol=1e-12)
    assert out['double']['p_lo'] == [0.0] * 4 and np.allclose(out['double']['p'], np.r_[w0, c0 + 0.001], rtol=0, atol=1e-15)
    assert len(out['anchored']['marks']) == 1 and out['anchored']['marks'][0][0] == WM.B4_MARKS[0]
    assert np.isclose(out['anchored']['marks'][0][1], ex(w0, c0 + 0.002), rtol=1e-12)
    assert out['anchored']['floor'] == WM.anchor_floor(f, pa, 2)
    assert out['trf']['status'] == 3 and out['trf']['nfev'] == 1234 and np.isclose(out['trf']['ratio'], ex(w0 * 1.01, c0 + 0.004), rtol=1e-12)
    assert out['eps'] == float(np.finfo(np.longdouble).eps)
    R = lambda ra, rl, it=100, acc=90, st='stall': dict(ratio=ra, ratio_ld=rl, it=it, acc=acc, stop=st)
    T = lambda ra, rl: dict(ratio=ra, ratio_ld=rl, status=3, nfev=2500)
    cell = lambda r_, k, start, e_opt, e_start, d, ld, an, ns, al, tr, fl=7e-13, mk=(): dict(
        r=r_, var='x', k=k, start=start, e_opt=e_opt, e_start=e_start, eps=1.0842e-19, double=d, nostall=ns, trf=tr,
        anchored=dict(an, floor=fl, marks=[list(x) for x in mk]), **{'long double': ld, 'anchored ld': al})
    rows = [cell(0.25, 6, 'block', 2.5e-11, 4.95, R(1.0019, 1.0019), R(0.99997, 0.99996), R(0.99996, 0.99996, 5000, 4000, 'flat'),
                 R(0.999962, 0.999962, 900, 800), R(0.9999599, 0.9999599, 60, 50), T(1.01, 1.01), mk=((100000, 0.999961),)),
            cell(0.25, 6, 'optimum', 2.5e-11, 1.00003, R(0.99997, 0.99997), R(0.99996, 0.999959), R(0.99996, 0.99996, 3000, 2000, 'flat'),
                 R(0.99996, 0.99996, 40, 30), R(0.99996, 0.99996, 70, 60), T(0.99999, 0.99999), fl=2e-13),
            cell(1.0, 12, 'block', 1.6e-9, 1.47, R(0.999997, 0.999997), R(0.999996, 0.999996), R(0.999996, 0.999996, 8000, 7000, 'stall'),
                 R(0.999996, 0.999996), R(0.999996, 0.999996), T(1.0, 1.0)),
            cell(1.0, 12, 'optimum', 1.6e-9, 1.0, R(0.999996, 0.999996), R(0.999996, 0.999996), R(0.999996, 0.999996, 9000, 8000, 'flat'),
                 R(0.999996, 0.999996), R(0.999996, 0.999996), T(0.999996, 0.999996))]
    s = '\n'.join(WM.b4_lines(rows))
    assert "in double or in long double (its ε 1.08e-19), or about an anchor" in s and "(`polish_block`, 500000 evaluations, no gradient test)" in s
    ex4 = lambda ra: f"{ra / 0.9999599 - 1:.1e}"
    assert (f"| 0.25 | x | 6 | 2.50e-11† | block | 4.95 | {0.9999599:.12f} | {ex4(1.0019)}, 100/90, stall | {ex4(0.99997)}, 100/90, stall | "
            f"{ex4(0.99996)}, 5000/4000, flat | {ex4(0.999962)}, 900/800, stall | {ex4(0.9999599)}, 60/50, stall | {ex4(1.01)}, 3/2500 |") in s
    ev = [r_[k_]['ratio_ld'] / r_[k_]['ratio'] - 1 for r_ in rows for k_, _ in WM.B4_ENDS]
    evl = [r_['long double']['ratio_ld'] / r_['long double']['ratio'] - 1 for r_ in rows]
    assert (f"against the exact one, minus 1: {min(ev):+.1e} … {max(ev):+.1e} over every end, {min(evl):+.1e} … {max(evl):+.1e} at the "
            "long-double solve's ends (2 of 4 below it") in s
    assert "The anchored double's ends from the block and from the stored optimum agree within 0.0e+00; its stops: flat 3, stall 1, after 3000 … 9000 iterations" in s
    assert "the floor the anchored double's frozen anchor rounding predicts, ‖P_Je₀‖²/(2‖r‖²) at its end (`anchor_floor`): 2.0e-13 … 7.0e-13" in s
    assert f"The stored optimum over the best end, minus 1: 1.6e-09 (τ/σ = 1, K = 12): {1.0 / 0.999996 - 1:+.1e}; 2.5e-11 (τ/σ = 0.25, K = 6): {1.00003 / 0.9999599 - 1:+.1e}." in s
    assert f"after 100000, 200000, 300000 iterations, over the best end, minus 1: 2.5e-11 (τ/σ = 0.25, K = 6) from the block: {0.999961 / 0.9999599 - 1:+.1e}." in s
    best = WM.b4_best(rows)
    assert best == {(0.25, 6): 0.9999599 * 2.5e-11, (1.0, 12): 0.999996 * 1.6e-9}
    assert WM.b4_lines(None) == [] and WM.b4_lines([]) == []


def test_b3_restart_judges_the_seed_from_the_start_with_the_paths_step(monkeypatch):
    """Referee round 77 (BR02, BR03): `b3_restart` judges every end by `seed_rule` from the start's amplitudes — a seed kernel that grows
    by the end is still judged — with the path's polish step: a 2σ path accepts a seed's move between 1σ and 2σ, a 1σ path rejects it."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(1.0) ** 2))
    w0, c0 = np.array([0.7, 0.005]), np.array([0.0, 0.3])
    end = c0 + np.array([0.0, 1.5 * WM.SIG])
    t = RT.ascent_target(K, f - K.phi(c0) @ w0, c0[1], 1.0)[0]
    assert abs(end[1] - t) / WM.SIG > 2

    def stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        return np.array([0.7, 0.3]), end.copy(), dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', stub)
    st = dict(w=w0, c=c0, s=np.ones(2), max_nfev=5000, gtol=1e-15)
    for step, ok in ((2.0, True), (1.0, False)):
        q = WM.b3_restart(1.0, 'x', K, f, st, 2, 1e-3, 2.0, 1, step, 5000, 1e-15, 0.0)
        assert all(np.isclose(q[k]['move'], 1.5) and q[k]['accepted'] is ok for k in WM.B3_KEYS)


def test_b3_job_on_a_variant_without_the_gradient_test(monkeypatch):
    """Referee round 77 (BR07): on the variant without a gradient test B3's restarts have none — the path's budget and B3_BIG, then the
    same two again."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.5) ** 2))
    v = WM.VARIANTS_B[5]
    assert v[5] is None
    start = dict(w=np.array([0.7]), c=np.array([0.0]), s=np.array([1.0]), max_nfev=v[3], gtol=None)
    monkeypatch.setattr(WM, 'path_polishes', lambda r_, var, stop_after=None: (None, K, f, [start]))
    calls = []

    def stub(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        calls.append((max_nfev, gtol))
        return w_, c_, dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', stub)
    q = WM.b3_job((0.5, v[0], dict(n0=0, k=1, e=2e-3, e_opt=1e-3, n_polish=1, first=dict(status=0, nfev=v[3]))))
    assert calls == [(v[3], None), (WM.B3_BIG, None), (v[3], None), (WM.B3_BIG, None)] and q['settings'] == dict(nfev=v[3], gtol=None)


def test_path_polishes_passes_the_settings_to_the_real_solve(monkeypatch):
    """Referee round 77 (PP04): `path_polishes`' recorder passes each call's own budget and gradient test to the solve it wraps."""
    from src.core import rich_transit as RT
    got = []

    def solve(K_, f_, w_, c_, s_, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        got.append((max_nfev, gtol))
        return w_, c_, dict(status=0, nfev=max_nfev)
    monkeypatch.setattr(RT, 'polish_block', solve)

    def fake_path(L, K, f, c0, rho, **kw):
        for i, g in enumerate((None, 1e-12)):
            RT.polish_block(K, f, np.array([0.5]), np.array([0.0]), np.array([1.0]), max_nfev=7 + i, info=True, gtol=g)
        return 'full'
    monkeypatch.setattr(RT, 'sliding_iss', fake_path)
    R, K, f, calls = WM.path_polishes(1.0, 'base')
    assert R == 'full' and got == [(7, None), (8, 1e-12)] and [(c['max_nfev'], c['gtol']) for c in calls] == got
    assert RT.polish_block is solve


def test_lm_fit_flat_stop(monkeypatch):
    """Referee rounds 77, 81: the anchored solve stops ('flat') once the exact rms over a period of anchor_every accepted steps falls by
    less than B4_FLAT of itself — here a period of 20 and B4_FLAT = ½ —; the period's accepted steps are counted across the re-anchors
    at stalls (the exact rms evaluated once at the start and once per period); without stall re-anchoring the solve stops at its first
    stall, earlier."""
    from src.core import rich_transit as RT
    K = RT.Bumps(WM.X_B, WM.SIG)
    ss = WM.sigma_star(0.5)
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * ss ** 2))
    co, ao, eo = WM.sample_fit(K.phi, f, WM.C_B + W.half_squeeze(8, WM.SIG, ss)[0])
    monkeypatch.setattr(WM, 'B4_FLAT', 0.5)
    p, it, acc, st = WM.lm_fit(f, ao, co, np.float64, anchored=True, anchor_every=20)
    assert st == 'flat' and acc == 20 and it < 100
    monkeypatch.setattr(WM, 'B4_FLAT', -1.0)                            # never flat
    real, n_ex, stalls = WM.exact_rms, [], []
    monkeypatch.setattr(WM, 'exact_rms', lambda *a: n_ex.append(1) or real(*a))
    p, it, acc, st = WM.lm_fit(f, ao, co, np.float64, anchored=True, anchor_every=30, max_iter=1500)
    assert st in ('cap', 'stall') and acc >= 60 and len(n_ex) == 1 + acc // 30
    pn, itn, accn, stn = WM.lm_fit(f, ao, co, np.float64, anchored=True, anchor_every=30, max_iter=1500, stall_reanchor=False)
    assert stn == 'stall' and itn < it and accn < acc
    pd_, itd, accd, std_ = WM.lm_fit(f, ao, co, np.float64, max_iter=1500)
    assert std_ == 'stall' and accd <= accn                              # direct subtraction stalls first: its rounding is not frozen


def test_the_exact_residual_and_the_anchor_floor(monkeypatch):
    """Referee round 81 (D3cb-1): `exact_resid` — every double and long double taken exactly, exp in 60-digit decimal — agrees with 90
    digits to the double's rounding, the long double entering as its exact hi + lo; at B4's 2.55·10⁻¹¹ optimum (τ/σ = 1/4, K = 6) the
    long-double rms by direct subtraction (`lm_rms`) is off by 10⁻¹³ … 10⁻⁹ relative, the residual by direct subtraction in double within
    10⁻¹⁵ of max|f| (a few of its units in the last place); `anchor_floor` is ‖P_Je₀‖²/(2‖r‖²) for a rounding e₀ (here planted: its part in the Jacobian's range and a part
    orthogonal to it)."""
    from decimal import Decimal
    from fractions import Fraction
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    j = next(j for j in WM.b4_jobs(raw['B'], raw['B3']) if j[0] == 0.25 and j[2] == 6 and j[3] == 'optimum')
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))
    p = np.r_[j[4], j[5]]
    r60, e60 = WM.exact_resid(f, p, 6)
    r90, e90 = WM.exact_resid(f, p, 6, prec=90)
    assert np.array_equal(r60, r90) and e60 == e90
    rd = f - np.exp(-(WM.X_B[:, None] - p[None, 6:]) ** 2 / (2 * WM.SIG ** 2)) @ p[:6]
    assert np.abs(rd - r60).max() < 1e-15 * np.abs(f).max() and np.abs(rd - r60).max() > 0
    rel = abs(WM.lm_rms(f, p, 6) / e60 - 1)
    assert 1e-13 < rel < 1e-9
    v = np.longdouble(1) / np.longdouble(3)
    hi = float(v)
    assert Fraction(WM._dec(v)) == Fraction(hi) + Fraction(float(v - np.longdouble(hi))) and WM._dec(v) != Decimal(hi)   # exact, in any context
    u = WM.X_B[:, None] - p[None, 6:]
    P = np.exp(-u ** 2 / (2 * WM.SIG ** 2))
    J = -np.c_[P, P * u / WM.SIG ** 2 * p[None, :6]]
    a = J[:, 0] * 1e-15
    b = np.random.default_rng(3).standard_normal(len(WM.X_B))
    b -= J @ np.linalg.lstsq(J, b, rcond=None)[0]                       # orthogonal to the range
    b *= 1e-16 / np.abs(b).max()
    monkeypatch.setattr(WM, 'exact_resid', lambda f_, p_, k_, prec=60: (rd - a - b, 0.0))
    fl = WM.anchor_floor(f, p, 6)
    assert np.isclose(fl, (a @ a) / (2 * ((rd - a - b) @ (rd - a - b))), rtol=1e-6)


def test_b4_job_on_two_cells_against_the_raw():
    """Referee round 81 (D3cb-2): the real `b4_job` — no fakes — on B4's two fastest cells (τ/σ = 1/4, K = 5, from the stored optimum and
    from the path's block) reproduces the raw's row: every solve's exact and long-double errors, iterations, accepted steps, stop and
    parameters, the anchored double's floor, TRF's status and evaluations."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    for start in ('optimum', 'block'):
        j = next(j for j in WM.b4_jobs(raw['B'], raw['B3']) if j[0] == 0.25 and j[2] == 5 and j[3] == start)
        out = WM.b4_job(j)
        x = next(q for q in raw['B4'] if (q['r'], q['k'], q['start']) == (0.25, 5, start))
        assert out['e_start'] == x['e_start']
        for k_, _ in WM.B4_ENDS:
            for f_ in ('ratio', 'ratio_ld', 'p', 'p_lo') + (('status', 'nfev') if k_ == 'trf' else ('it', 'acc', 'stop')):
                assert out[k_][f_] == x[k_][f_], (start, k_, f_)
        assert out['anchored']['floor'] == x['anchored']['floor']


def test_lm_fit_long_descent_from_the_paths_block():
    """Referee round 77 (LM15): from B4's farthest start — the path's block at τ/σ = 1/4, K = 5, ten times the optimum's error — the
    double solve reaches the optimum within 10⁻⁶ in about 8000 iterations, its damping relaxing after its rejections."""
    import exp_cohesion_leak as E
    raw = E.load(os.path.join(WM.OUT, 'width_mismatch_raw.npz'))
    j = next(j for j in WM.b4_jobs(raw['B'], raw['B3']) if j[0] == 0.25 and j[2] == 5 and j[3] == 'block')
    f = WM.A_B * np.exp(-(WM.X_B - WM.C_B) ** 2 / (2 * WM.sigma_star(0.25) ** 2))
    rf = float(np.sqrt(np.mean(f ** 2)))
    assert WM.lm_rms(f, np.r_[j[4], j[5]], 5) / rf / j[6] > 10
    p, it, acc, st = WM.lm_fit(f, j[4], j[5], np.float64, max_iter=20_000)
    assert st == 'stall' and 5000 < it < 20_000 and acc > 5000 and abs(WM.lm_rms(f, p, 5) / rf / j[6] - 1) < 1e-6


def test_part_b4_and_main_only_b4(monkeypatch, tmp_path):
    """Referee round 77 (PB01, MN01): `part_B4` maps `b4_job` over every job; `main --only B4` passes the raw's B and B3 to it and stores
    its result."""
    import multiprocessing
    import exp_cohesion_leak as E

    class FakePool:
        def __init__(self, n):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, jobs, chunksize=1):
            return [fn(j) for j in jobs]
    monkeypatch.setattr(multiprocessing, 'Pool', FakePool)
    monkeypatch.setattr(WM, 'b4_jobs', lambda B, B3: [('j', i) for i in range(5)])
    monkeypatch.setattr(WM, 'b4_job', lambda j: dict(job=j))
    assert WM.part_B4('B', 'B3') == [dict(job=('j', i)) for i in range(5)]
    E.save(str(tmp_path / 'width_mismatch_raw.npz'), B=[dict(r=1.0, var='base')], B3=[dict(r=1.0, var='base', k=None)])
    got, saved = [], []
    monkeypatch.setattr(WM, 'OUT', str(tmp_path))
    monkeypatch.setattr(WM, 'part_B4', lambda B, B3: (got.append((B, B3)), ['b4'])[1])
    monkeypatch.setattr(WM, 'report', lambda *a, **k: '')
    monkeypatch.setattr(E, 'save', lambda path, **kw: saved.append(kw))
    monkeypatch.setattr(sys, 'argv', ['x', '--only', 'B4'])
    WM.main()
    assert len(got) == 1 and list(got[0][0]) == [dict(r=1.0, var='base')] and list(got[0][1]) == [dict(r=1.0, var='base', k=None)]
    assert saved[-1]['B4'] == ['b4'] and list(saved[-1]['B']) == [dict(r=1.0, var='base')]
