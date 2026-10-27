"""
Proposition D‴.42 (PARADIGM.md §5d): several particles per class. Checked here: at an exact fit with each class's particles together,
log det ∇²_ΘJ minus its one-particle-per-class value at the same class masses is a constant, whatever the split within the classes; the
class masses' drift and diffusion in Katzenberger's limit are 3/m times the one-particle-per-class ones in the slow-replicator limit, and
the fractions within a class do not diffuse; the correction per class (`fixman_classes`) is the one-particle-per-class correction at the
class masses times 3/m, given to every particle of the class, and its tracked contributions are summed per class.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
import exp_gauge_classes as C                                            # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def test_the_determinant_does_not_see_the_split():
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    Wc = np.array([0.30, 0.35, 0.35])
    l1 = C.log_det(loop, E.exact_rep(net, Wc), 3 * Wc)
    d = [C.log_det(loop, *C.exact_m(net, Wc, s, 2)) - l1 for s in ((0.5, 0.5, 0.5), (0.1, 0.8, 0.35))]
    Wc2 = np.array([0.40, 0.25, 0.35])
    d.append(C.log_det(loop, *C.exact_m(net, Wc2, (0.3, 0.6, 0.5), 2)) - C.log_det(loop, E.exact_rep(net, Wc2), 3 * Wc2))
    assert np.ptp(d) < 1e-5


def test_the_class_masses_follow_the_law_at_m_and_the_fractions_do_not_diffuse():
    r = C.b_job((2, 1e-5, (0.30, 0.35), (0.3, 0.6, 0.8)))
    assert np.allclose(r['cls'] / r['cls3'], 0.5, rtol=2e-4) and np.allclose(r['Dcls'] / r['Dcls3'], 0.5, rtol=1e-6)
    assert np.abs(r['frac']).max() < 1e-3 * np.abs(r['cls']).max() and r['Dfrac'] < 1e-12 * r['Dcls'].max()


def test_the_correction_per_class():
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35])
    loop3 = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4, fixman_A=E.A_BUMP.copy())
    c3 = loop3.fixman_correction(E.exact_rep(net, W), 3 * W)
    lab = np.repeat(np.arange(3), 2)
    loop6 = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4, fixman_A=E.A_BUMP.copy(), fixman_classes=lab)
    Th, beta = C.exact_m(net, W, (0.2, 0.6, 0.5), 2)
    c6 = loop6.fixman_correction(Th, beta)
    assert np.allclose(c6, 0.5 * c3[lab], rtol=1e-12)
    # the plug-in form per class reads the class's current contribution: on M the same
    plug = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=1e-4, fixman_classes=lab)
    assert np.allclose(plug.fixman_correction(Th, beta), c6, rtol=1e-12)
    # the tracked contributions are summed per class: on M without dither they stay at the bumps' amplitudes
    tr = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=0.0, fixman_A=E.A_BUMP + 0.01, fixman_TA=5.0,
                                   fixman_classes=lab)
    tr.run_dithered(Th, beta, 2.5, 0.25, 0.0, np.random.default_rng(0))
    assert np.allclose(tr.fixman_A - E.A_BUMP, 0.01 * (1 - 0.25 / 5.0) ** 10, rtol=1e-6)
