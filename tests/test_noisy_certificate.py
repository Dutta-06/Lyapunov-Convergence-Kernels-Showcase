"""Proposition D‴.56: the transit under measurement noise — the sliding path's linear response (`rich_transit.path_response`), the window
and the stopped fit (`experiments/exp_noisy_certificate.py`)."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

from src.core import rich_transit as RT
import exp_noisy_certificate as NC
import exp_transit_mass as TM
import exp_mirror_transit as X

CS, AS = np.array([-0.52, 0.02, 0.61]), np.array([0.8, -0.35, 0.9])


def test_block_response_at_an_exact_fit_is_the_projection_off_the_tangent_space():
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, CS, AS)
    Rm = RT.block_response(K, f, AS, CS)
    J = np.c_[K.phi(CS), K.dphi(CS) * AS[None, :]]
    assert np.abs(Rm @ Rm - Rm).max() < 1e-10 and np.abs(Rm - Rm.T).max() < 1e-10
    assert np.abs(Rm @ J).max() < 1e-10
    assert abs(np.trace(Rm) - (K.xs.size - 2 * CS.size)) < 1e-8                 # rank N − 2K


def test_sliding_iss_stops_at_max_entries_before_relaxing():
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, CS, AS)
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, max_entries=2)
    assert [e[2] for e in R['events']] == ['enter', 'enter'] and len(R['relaxed']) == 1 and len(R['hist']) == 2


def test_path_response_is_the_derivative_of_the_noisy_path():
    """Central differences of the noisy path's entry times and of its duals at the K-th entry against T·ξ and U·ξ."""
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, CS, AS)
    R0 = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    P = RT.path_response(K, f, R0, c0)
    xi = np.random.default_rng(3).standard_normal(K.xs.size)
    h, nK = 1e-5, CS.size
    kw = dict(polish=True, max_events=nK + 3, tol=0.0, max_entries=nK + 1)
    Rp = RT.sliding_iss(L, K, f + h * xi, c0, np.ones(c0.size), **kw)
    Rm = RT.sliding_iss(L, K, f - h * xi, c0, np.ones(c0.size), **kw)
    tp, tm = (np.array([e[0] for e in R_['events'][:nK]]) for R_ in (Rp, Rm))
    assert np.abs((tp - tm) / (2 * h) - P['T'] @ xi).max() < 1e-5 * np.abs(P['T'] @ xi).max()
    off = [j for j in range(c0.size) if j not in P['entrants']]
    up, um = (RT.path_duals(K, f + s_ * h * xi, R_)['U'][nK] for s_, R_ in ((1, Rp), (-1, Rm)))
    fd = ((up - um) / (2 * h))[off]
    assert np.abs(fd - (P['U'] @ xi)[off]).max() < 1e-5 * np.abs((P['U'] @ xi)[off]).max()
    assert np.abs(P['u_all'][-1] - P['u']).max() == 0 and np.abs(P['U_all'][-1] - P['U']).max() == 0


def test_window_next_kernel_and_oracle_at_small_noise():
    """s = 0.001 on a separated target: the path enters exactly E, the first noise kernel enters after the first-order window W (to 1%:
    the relative error is O(s), the residual R_Kξ + O(|ξ|²)) and is the predicted one, and the fit stopped at the K-th entry is the oracle
    — its error is the noise's projection on the tangent space (to 1%)."""
    L, K, c0, f, cert, R0, P = NC.noiseless(CS, AS)
    nK, nx = CS.size, K.xs.size
    off = [j for j in range(c0.size) if j not in cert['e']]
    for seed in (1, 2):
        xi = 0.001 * np.random.default_rng(seed).standard_normal(nx)
        R = NC.noisy(L, K, c0, f, xi, nK)
        ent = [e for e in R['events'] if e[2] == 'enter']
        assert sorted(e[1] for e in ent[:nK]) == sorted(cert['e']) and len(ent) == nK + 1
        v = K.phi(c0).T @ (P['Rlast'] @ xi) / nx
        u = P['u'] + P['U'] @ xi
        W = {j: (1 - np.sign(v[j]) * u[j]) / abs(v[j]) for j in off}
        jn = min(W, key=W.get)
        assert ent[nK][1] == jn and abs((ent[nK][0] - ent[nK - 1][0]) / W[jn] - 1) < 1e-2
        hK = R['hist'][nK]
        err = np.sqrt(np.mean((K.phi(hK[1]) @ hK[2] - f) ** 2))
        err_T = np.sqrt(np.mean(((np.eye(nx) - P['Rlast']) @ xi) ** 2))
        assert abs(err / err_T - 1) < 1e-2
