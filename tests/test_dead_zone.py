"""Proposition D‴.57: the transit's dead zone — the energy dead zone's frozen block in the sliding limit, the gated loop
(`experiments/exp_dead_zone.py`), and the dead zone in the closed loop's C integrator (`closed_loop.simulate(gate=...)`)."""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
warnings.filterwarnings('ignore')

from src.core import closed_loop as CL
from src.core import mirror_free as MF
from src.core import rich_transit as RT
import exp_dead_zone as DZ
import exp_mirror_transit as X
import exp_noisy_certificate as NC
import exp_transit_mass as TM


def test_gated_relax_stops_where_the_residuals_rms_reaches_the_threshold():
    """The noiseless path on target 0: the (K − 1)-th block relaxed and the K-th entrant at its seed, gated at half the (K − 1)-th
    relaxed rms: the block's fast flow crosses the threshold and stops on it."""
    cs, As = NC.targets()[0]
    L, K, c0 = TM.setup()
    f = X.bumps_target(K, cs, As)
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    nK = cs.size
    _, c_prev, w_prev, rms_prev = R['hist'][nK - 1]
    j = [e for e in R['events'] if e[2] == 'enter'][nK - 1][1]
    act = np.flatnonzero(w_prev != 0)
    sg = np.sign(K.phi(c_prev)[:, j] @ (f - K.phi(c_prev) @ w_prev))
    thr = 0.5 * rms_prev
    crossed, cc, wc = DZ.gated_relax(L, K, f, c_prev[np.r_[act, j]], np.r_[w_prev[act], sg * 1e-7], thr)
    assert crossed
    assert abs(np.sqrt(np.mean((f - K.phi(cc) @ wc) ** 2)) / thr - 1) < 1e-8
    assert np.abs(wc[-1]) > 1e-3                                         # the entrant has grown


def test_energy_dead_zone_freezes_the_last_block_on_the_energy_identity():
    """Target 0, s = 0.01 (Part A's first draw): the dead zone stops at the K-th entry; the frozen block's unrelaxed part u obeys |u|² =
    (1 + δ)²s²N − |r_K|² to 10⁻³ (first order), and the frozen fit is the oracle's displaced by u (an identity of the bookkeeping)."""
    cs, As = NC.targets()[0]
    L, K, c0, f, cert, R0, P = NC.noiseless(cs, As)
    nK, nx = cs.size, K.xs.size
    s = 0.01
    xi = s * np.random.default_rng(5700).standard_normal(nx)
    d = DZ.a_draw(L, K, c0, f, sorted(cert['e']), xi, s, nK)
    assert d['exact'] and d['k'] == nK and d['crossed']
    assert abs(d['u2'] / d['u2_pred'] - 1) < 1e-3
    assert d['err_dz'] > d['err_or']
    u = np.sqrt(d['u2'] / nx)
    assert abs(d['err_dz'] ** 2 - (d['err_or'] ** 2 + u ** 2 - 2 * d['cross'] * d['err_or'] * u)) < 1e-12 * d['err_dz'] ** 2 + 1e-18


def test_gated_loops_form_the_features_kernels_and_follow_the_sliding_limit():
    """Part B's draw (target 0, s = 0.03, draw 0) at ε = 10⁻⁸ to τ_K + 3W: ungated the loop has formed noise kernels; either gate forms
    the three features' kernels only, the energy dead zone's fit is the limit's frozen fit and the formation dead zone's the oracle's."""
    cs, As = NC.targets()[0]
    out = DZ.b_job((0, cs, As, 0.03, 0))
    assert out['exact'] and out['k_lim'] == cs.size
    assert out['none']['formed'] > cs.size and out['energy']['formed'] == out['formation']['formed'] == cs.size
    assert abs(out['energy']['err'] - out['err_dz_lim']) < 1e-8 and abs(out['formation']['err'] - out['err_or_lim']) < 1e-9
    assert out['energy']['err'] > 1.5 * out['formation']['err']


LAW = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
SIG, KE, OM = 0.12, 10.0, 0.2
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
C0 = np.linspace(-1, 1, 21)
TS = np.linspace(0, 1.2e4, 41)
OPEN_TOL = 1e-8                                                          # measured: 2.5·10⁻¹¹ on the centres


def test_closed_dead_zone_freezes_the_adaptation_exactly():
    """A gate that never opens (the filtered energy starts at 0 below the thresholds and the filter is slow): every kernel's drive is
    zero — the coefficients stay at 0 and the centres at their start, exactly, under either gate."""
    for mode in ('energy', 'formation'):
        r = CL.simulate(LAW, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), TS,
                        gate=dict(mode=mode, Tf=1e12, lo=10.0, hi=20.0, z0=0.0, theta=0.01))
        assert np.abs(r['w']).max() == 0.0 and np.abs(r['c'] - C0).max() == 0.0
        assert 0.0 < r['z'][-1] < 1e-6


def test_open_dead_zone_is_the_ungated_loop():
    """A gate that never closes (thresholds far below the error's energy): the ungated loop, up to the integrator's error control on the
    extra state."""
    r0 = CL.simulate(LAW, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), TS)
    assert np.abs(r0['w'][-1]).max() > 0.05                               # kernels have formed
    for mode in ('energy', 'formation'):
        r = CL.simulate(LAW, 1e-2, C0, SIG, (CSTAR, ASTAR, SIG), KE, (0.0, 1.0, OM), TS,
                        gate=dict(mode=mode, Tf=1e12, lo=0.5e-9, hi=1e-9, theta=0.01))
        assert np.abs(r['w'] - r0['w']).max() < OPEN_TOL and np.abs(r['c'] - r0['c']).max() < OPEN_TOL
