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




def test_fitted_skeleton_fits_negative_features():
    """D(a)'s opposite-sign pair 0.3 apart at the ratio 1 (exact with E): the fitted skeleton's boxes follow the features' signs — with every
    box positive (referee round 46's mutant K7) the negative feature could not be fitted."""
    L, K, c0 = TM.setup()
    cs, As = np.array([-0.1, 0.2]) + TC.LAT_OFF, np.array([1.0, -1.0])
    fit = TC.skeleton(K, c0, cs, As, fit=True)
    assert fit['ok'] and fit['margin'] > 0 and sorted(fit['entrants']) == sorted(TC.entrants(K, c0, cs)[0])
    assert TC.job_1d((cs, As))['out'] == 'nearest'


def test_atoms_baseline_drops_the_tangent_directions():
    """The atoms-only baseline subtracts the fitted features' atoms and nothing else: for one feature its residual after the entry is zero,
    and on D(a)'s lattice pair at the ratio 0.02 its margin differs from the linear skeleton's (which also removes the tangent directions)."""
    L, K, c0 = TM.setup()
    a = TC.skeleton(K, c0, np.array([0.137]), np.array([-0.7]), atoms=True)
    assert a['ok'] and a['margin'] > 0
    cs, As = np.array([-0.1, 0.2]) + TC.LAT_OFF, np.array([1.0, 0.02])
    assert abs(TC.skeleton(K, c0, cs, As, atoms=True)['margin'] - TC.skeleton(K, c0, cs, As)['margin']) > 1e-4


def test_entrants_are_the_largest_correlations_near_the_window_edge():
    """Near the window's edge the samples' largest correlation is not the nearest kernel (referee round 46): a feature at 0.97 is entered by
    the kernel at 0.9 (index 19), which E must name for the path to count as exact with E."""
    L, K, c0 = TM.setup()
    cs, As = np.array([0.97]), np.array([0.8])
    e = TC.entrants(K, c0, cs)[0]
    assert e == [19] and int(np.argmin(np.abs(c0 - 0.97))) == 20
    assert TC.job_1d((cs, As))['out'] == 'nearest'


def test_ascent_target_ends_at_the_feature():
    """`ascent_target`: a seed of the feature's sign started 0.08 away climbs its correlation with the feature to the feature's centre
    (within the step σ/40); of the opposite sign it climbs away."""
    L, K, c0 = TM.setup()
    r = TC.X.bumps_target(K, np.array([0.137]), np.array([0.7]))
    c = TC.RT.ascent_target(K, r, np.array([0.137 + 0.08]), 1.0)
    assert abs(float(c[0]) - 0.137) < K.sig / 40 + 1e-12
    assert abs(float(TC.RT.ascent_target(K, r, np.array([0.137 + 0.08]), -1.0)[0]) - 0.137) > 0.08


def test_outcome_tolerances():
    """`outcome`: exact with E needs the residual below RMS_TOL, the active set E, no exits and exactly K events."""
    act = np.zeros(21, bool); act[[3, 9]] = True
    ev = [(0.0, 3, 'enter', 0.0), (1.0, 9, 'enter', 0.0)]
    R = lambda rms, events=ev: dict(events=events, act=act, hist=[(0.0, None, None, 1.0), (1.0, None, None, rms)])
    assert TC.outcome(R(0.5 * TC.RMS_TOL), [3, 9], 2) == 'nearest'
    assert TC.outcome(R(2 * TC.RMS_TOL), [3, 9], 2) == 'not'
    assert TC.outcome(R(0.5 * TC.RMS_TOL, ev + [(2.0, 5, 'enter', 0.0)]), [3, 9], 2) == 'other'
    assert TC.outcome(R(0.5 * TC.RMS_TOL), [3, 8], 2) == 'other'


def test_skeletons_on_the_targets_that_break_their_mutants():
    """Round 46's surviving mutants, each on a stored target it gets wrong: the linear skeleton's tangent directions at the features (not
    the starting points; B's configuration 31 at r = 3, exact with E); the fitted skeleton's boxes with the features' signs (D(a)'s
    opposite-sign pair 0.5 apart at the ratio 50) and its fit started at the features (B's configuration 11 at r = 3)."""
    L, K, c0 = TM.setup()
    cs = np.array([0.628354315580911, -0.04456980956766876, -0.6300330902222451])
    As = np.array([-0.5978787114809792, -0.31619518893970555, -0.0021571116965809])
    sk = TC.skeleton(K, c0, cs, As)
    assert sk['ok'] and sk['margin'] > 0.01                            # tangents at the starting points: −0.97
    assert TC.job_1d((cs, As))['out'] == 'nearest'
    fit = TC.skeleton(K, c0, np.array([-0.1, 0.4]) + TC.LAT_OFF, np.array([1.0, -50.0]), fit=True)
    assert fit['ok'] and fit['margin'] > 0.1                           # every box positive: −0.99
    cs = np.array([-0.26146688489923686, 0.6985887704301363, -0.6172569386174308])
    As = np.array([0.9608083993684658, 0.001283451133584242, 0.0018132196297868347])
    fit = TC.skeleton(K, c0, cs, As, fit=True)
    assert fit['ok'] and fit['margin'] > 0.01                          # started at the starting points: −0.34
    assert TC.job_1d((cs, As))['out'] == 'nearest'


def test_outcome_counts_kernels_and_its_tolerance_is_absolute():
    """`outcome`: the residual's tolerance is 10⁻⁸ (absolute), and an exact fit with more kernels than features is neither exact with E
    nor with other kernels (Part F's four-kernel fit)."""
    act = np.zeros(21, bool); act[[3, 9]] = True
    ev = [(0.0, 3, 'enter', 0.0), (1.0, 9, 'enter', 0.0)]
    R = lambda rms, a=act, events=ev: dict(events=events, act=a, hist=[(0.0, None, None, 1.0), (1.0, None, None, rms)])
    assert TC.outcome(R(5e-9), [3, 9], 2) == 'nearest' and TC.outcome(R(2e-8), [3, 9], 2) == 'not'
    a3 = act.copy(); a3[5] = True
    assert TC.outcome(R(1e-12, a3, ev + [(2.0, 5, 'enter', 0.0)]), [3, 9], 2) == 'not'


def test_the_polish_safeguard_on_stored_paths():
    """Round 46's surviving safeguard mutants, on stored paths (both `sliding_iss` and `sliding_iss_nd`, the experiment's): A's target 114
    accepts a near-seed kernel's polish 3.8σ from its start because it lands at its correlation's ascent (bounded by the start only, by
    the ascent of the wrong sign or by a quarter width it is rejected three times); A's target 12 and B's configuration 28 at r = 3 reject
    2 and 4 polishes, each retried (B's configuration 14 at r = 2 ends exact with E only with the near-seed scale relative to the block)."""
    L, K, c0 = TM.setup()
    T = {'A114': ([0.6541004699596475, 0.14717182832824738, -0.1741542573559407],
                  [-0.9993392739038046, -0.5218652653923312, -0.8701696976682551], 'other', 0),
         'A12': ([-0.3482759573852491, -0.0002953631766601239, -0.6856589204161185],
                 [0.8834332947114262, 0.9508709095889194, 0.40545651538789224], 'other', 2),
         'B115': ([0.23946960475117063, -0.19182582990096742, -0.6226746645823192],
                  [-0.017550238780768175, 0.0041097263697728054, 0.0032936667939692126], 'nearest', 4),
         'B58': ([-0.05502623411688401, -0.4293485291166097, 0.5138963933303842],
                 [0.044839474500752985, 0.08864282863781996, 0.03386534213286307], 'nearest', 1)}
    for name, (cs, As, out, rej) in T.items():
        cs, As = np.array(cs), np.array(As)
        e = TC.entrants(K, c0, cs)[0]
        R1 = TC.RT.sliding_iss(L, K, TC.X.bumps_target(K, cs, As), c0, np.ones(c0.size), polish=True)
        Rn = TC.job_1d((cs, As))
        assert (TC.outcome(R1, e, 3), int(R1['polish_rejected'])) == (out, rej), name
        assert (Rn['out'], Rn['rej']) == (out, rej), name
