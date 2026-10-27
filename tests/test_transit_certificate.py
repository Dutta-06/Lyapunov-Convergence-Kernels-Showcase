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
