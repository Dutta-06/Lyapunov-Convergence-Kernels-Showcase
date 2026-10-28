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
    assert '| 0.25 | base | 1 | 2, 1.00e-03, False | 1/2 | 0, 0, 0 | 0 | 0 |' in md
    C = [dict(s=1e-3, q=0, k=[1, 2, 3], rms=[5e-3, 1.05e-3, 0.9e-3], err=[0.02, 0.005, 0.004], at=[1, 2, 3], gaps=[0.5, 0.2], exits=0,
              stiff=0, runaway=False, relaxed=[True] * 3, rf=0.25)]
    md = WM.report(C=C)
    assert '| 0.001 | 2 [2, 2] (1/1 stopped) |' in md
    D = dict(rows=[dict(K=10, t=0.5, ratio=1.0, n_terms=58)], lead=[dict(j=1, lead=(-1, 1), next_dd=0.0), dict(j=2, lead=(-1, 2), next_dd=0.0)])
    assert 'at every j' in WM.report(D=D) and 'NOT' not in WM.report(D=D)
    D['lead'].append(dict(j=3, lead=(1, 2), next_dd=0.0))
    assert 'NOT at every j' in WM.report(D=D)
