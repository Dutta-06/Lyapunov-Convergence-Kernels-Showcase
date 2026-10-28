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
