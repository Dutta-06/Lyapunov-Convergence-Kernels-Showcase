"""Theorem D‴.55: the transit's certificate (`experiments/exp_transit_certificate.py`)."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

import exp_transit_certificate as TC
import exp_transit_mass as TM


def test_one_feature_certificate_is_the_normalised_profile():
    L, K, c0 = TM.setup()
    cs, As = np.array([0.137]), np.array([-0.7])
    cert = TC.certificate(K, c0, cs, As)
    e = cert['e'][0]
    assert e == int(np.argmin(np.abs(c0 - 0.137)))
    # η(c_e) = 1, and η elsewhere is the sample-Gram profile normalised by the nearest kernel's, signed
    g = K.phi(c0).T @ K.phi(cs)[:, 0] / K.xs.size
    assert abs(cert['eta'][e] + 1.0) < 1e-12                           # sign(A) = −1
    assert np.allclose(cert['eta'], -g / g[e])
    assert 0 < cert['m'] < 1 and cert['q'] == 0.0


def test_certificate_depends_on_signs_and_positions_only():
    L, K, c0 = TM.setup()
    cs = np.array([-0.5, 0.05, 0.62])
    a = TC.certificate(K, c0, cs, np.array([0.9, -0.6, 0.75]))
    b = TC.certificate(K, c0, cs, np.array([1e-3, -50.0, 2.0]))
    assert np.allclose(a['eta'], b['eta']) and a['e'] == b['e'] and a['m'] == b['m']


def test_separated_features_are_recovered_by_their_nearest_kernels():
    r = TC.job_1d((np.array([-0.52, 0.02, 0.61]), np.array([0.8, -0.05, 0.9])))
    assert r['m'] > 0.8 * r['q'] and r['out'] == 'nearest' and r['rms'] < TC.RMS_TOL


def test_a_point_between_two_same_sign_features_intrudes():
    # two same-sign features 0.27 (2.25σ) apart, a starting point midway: the certificate there exceeds one and that kernel enters
    r = TC.job_1d((np.array([-0.135, 0.135, 0.7]), np.array([0.8, 0.8, 0.5])))
    assert r['m'] < -0.8 * r['q'] and r['out'] != 'nearest'
    L, K, c0 = TM.setup()
    assert r['argm'] == int(np.argmin(np.abs(c0)))                    # the kernel at 0, midway
    assert r['argm'] in r['order']


def test_skeleton_for_one_feature_is_the_single_feature_margin():
    """K = 1: the skeleton enters the feature's kernel and its margin is 1 − max_{j≠e}⟨φφ_j⟩/⟨φφ_e⟩, the certificate's and the
    single-feature margin alike."""
    L, K, c0 = TM.setup()
    cs, As = np.array([0.137]), np.array([-0.7])
    cert = TC.certificate(K, c0, cs, As)
    sk = TC.skeleton(K, c0, cs, As)
    assert sk['ok'] and sk['entrants'] == cert['e']
    assert abs(sk['margin'] - cert['m']) < 1e-12 and abs(cert['m_vor'] - cert['m']) < 1e-12 and abs(cert['m_rel'] - cert['m']) < 1e-12


def test_skeleton_follows_the_amplitude_ratio_where_the_certificate_cannot():
    """Round 43's counterexample: two opposite-sign features on starting points (10⁻⁶ off), 0.3 apart. The amplitude-free certificate is
    the same at every ratio (m = 1.057q); the path is exact with E at equal amplitudes and not when the second is halved (the kernel at 0.3
    slides onto the feature at 0.2) — and the skeleton, which follows the order of the entries and the screening, says the same."""
    L, K, c0 = TM.setup()
    cs = np.array([-0.1, 0.2]) + 1e-6
    out = {}
    for a2 in (-1.0, -0.5):
        As = np.array([1.0, a2])
        r = TC.job_1d((cs, As))
        sk = TC.skeleton(K, c0, cs, As)
        assert abs(r['m'] / r['q'] - 1.057) < 2e-3
        assert (sk['margin'] > 0) == (r['out'] == 'nearest') == sk['ok']
        out[a2] = r['out']
    assert out == {-1.0: 'nearest', -0.5: 'other'}


def test_fitted_skeleton_resolves_the_linear_skeletons_near_tie():
    """Part D(a)'s lattice pair 0.3 apart, same signs, the amplitude ratio 0.02 (and 50): the path is exact with E, the linear skeleton
    lets the intruder in by a near tie (m_sk = −0.0030), the fitted skeleton — each block at its own least-squares fit — the rightful
    kernel by +10⁻⁴: the linear skeleton's misses are its linearisation."""
    L, K, c0 = TM.setup()
    cs = np.array([-0.1, 0.2]) + TC.LAT_OFF
    for As in (np.array([1.0, 0.02]), np.array([1.0, 50.0])):
        lin, fit = TC.skeleton(K, c0, cs, As), TC.skeleton(K, c0, cs, As, fit=True)
        assert not lin['ok'] and -4e-3 < lin['margin'] < -2e-3
        assert fit['ok'] and 0 < fit['margin'] < 1e-3 and sorted(fit['entrants']) == [9, 12]
    assert TC.job_1d((cs, np.array([1.0, 0.02])))['out'] == 'nearest'


def test_relative_margin_compares_with_the_entrants_own_values():
    """m_rel = min_k s_kη(c_{e_k}) − max_{j∉E}|η(c_j)|: for an opposite-sign pair η(c_{e_k}) < 1, so m_rel < m."""
    L, K, c0 = TM.setup()
    cs, As = np.array([-0.1, 0.2]) + 1e-6, np.array([1.0, -1.0])
    c = TC.certificate(K, c0, cs, As)
    off = [j for j in range(c0.size) if j not in c['e']]
    own = min(np.sign(As[k]) * c['eta'][c['e'][k]] for k in range(2))
    assert own < 1 and abs(c['m_rel'] - (own - np.abs(c['eta'][off]).max())) < 1e-14 and c['m_rel'] < c['m']


def test_outcome_needs_the_nearest_kernels_and_no_exits():
    ev = lambda *xs: [(float(i), j, k, 0.0) for i, (j, k) in enumerate(xs)]
    act = np.zeros(21, bool); act[[3, 9]] = True
    hist = [(0.0, None, None, 1.0), (1.0, None, None, 1e-12)]
    assert TC.outcome(dict(events=ev((3, 'enter'), (9, 'enter')), act=act, hist=hist), [3, 9], 2) == 'nearest'
    assert TC.outcome(dict(events=ev((3, 'enter'), (5, 'enter'), (5, 'exit'), (9, 'enter')), act=act, hist=hist), [3, 9], 2) == 'other'
    assert TC.outcome(dict(events=ev((3, 'enter'), (9, 'enter')), act=act, hist=hist[:1]), [3, 9], 2) == 'not'


def test_bounded_polish_keeps_a_seed_level_entrant_from_jumping():
    """Round 43's target 267 of Part A with a short fast flow (t_fast = 10³): unbounded, the polish carries the second entrant, still at
    its seed coefficient, 6.7σ onto another feature and the path ends with other kernels; bounded (a near-seed kernel may move by one
    kernel width at most), it ends exactly with E, as with the default fast flow."""
    rng = np.random.default_rng(TC.SEED_A)
    tg = [(cs, s * a) for dmin in TC.DMIN_A for cs, s, a, u in TC.targets(rng, TC.N_A, dmin)]
    cs, As = tg[267]
    L, K, c0 = TM.setup()
    cert = TC.certificate(K, c0, cs, As)
    f = TC.X.bumps_target(K, cs, As)
    R_free = TC.RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, t_fast=1e3, polish_step=np.inf)
    R_bnd = TC.RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, t_fast=1e3)
    assert TC.outcome(R_free, cert['e'], 3) == 'other' and TC.outcome(R_bnd, cert['e'], 3) == 'nearest'
    assert R_bnd['polish_rejected'] > 0


def test_band_counts_the_targets_outside_it():
    rows = [dict(m=0.5, q=0.1, out='nearest'), dict(m=-0.5, q=0.1, out='other'), dict(m=0.01, q=0.1, out='other'),
            dict(m=0.2, q=0.1, out='other')]
    txt = TC.band(rows)
    assert 'C = 0.5: 2 of 3' in txt and 'C = 1: 2 of 3' in txt
