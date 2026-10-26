"""
Proposition D‴.40 (PARADIGM.md §5d): the gauge under general noise. Checked here: the noise mask of the dithered integrator dithers only
the chosen coordinates (the amplitudes then move by their drift alone); the averaged fitness of the trace formula is, for isotropic noise
2τ diag(G), minus Proposition D‴.38's correction at κ_F = 1 (−(τ/2)∂_βℓ along M), the same for noise on the amplitudes only, and zero to
leading order for noise on the centres and widths only.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import gauge_noise as GN                                   # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def test_the_noise_mask_dithers_only_the_chosen_coordinates():
    # Proposition D‴.40's mask: with (0, 1, 1) the amplitudes move by their drift alone, the centres and widths as without the mask
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    T0 = E.exact_rep(net, np.ones(3) / 3) + 0.01
    out = {}
    for mask in (None, (0.0, 1.0, 1.0)):
        rec = []
        loop.run_dithered(T0, np.ones(3), 0.25, 0.25, 1e-3, np.random.default_rng(3), rec_every=1, callback=lambda s_, T, b: rec.append(T.copy()),
                          mask=mask)
        out[mask] = rec[-1]
    drift = loop.drift(T0, np.ones(3))[0]
    assert np.allclose(out[(0.0, 1.0, 1.0)][:, 0], T0[:, 0] + 0.25 * drift[:, 0], rtol=0, atol=1e-15)
    assert np.array_equal(out[(0.0, 1.0, 1.0)][:, 1:], out[None][:, 1:]) and not np.allclose(out[None][:, 0], out[(0.0, 1.0, 1.0)][:, 0])


def _instance(W):
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    beta = 3 * np.asarray(W, float)
    G9 = np.tile(E.G, 3)
    gradJ = lambda t: -loop.drift(t.reshape(3, 3), beta)[0].ravel() / G9
    u = lambda t: loop.fitness(loop.residual(t.reshape(3, 3), beta), t.reshape(3, 3))
    return net, base, loop, beta, G9, gradJ, u


def test_the_trace_formula_reduces_to_the_log_determinant_for_isotropic_noise_and_vanishes_without_amplitude_noise():
    tau = 1e-4
    W = np.array([0.30, 0.35, 0.35])
    net, base, loop, beta, G9, gradJ, u = _instance(W)
    t = E.exact_rep(net, W).ravel()
    iso = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(G9), G9)[0]
    fix = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=1.0, fixman_tau=tau).fixman_correction(t.reshape(3, 3), beta)
    flow = lambda ub: GN.replicator_flow(W, ub)
    assert np.abs(flow(iso) + flow(fix)).max() < 1e-5 * np.abs(flow(fix)).max()            # ū = −(τ/2)∂_βℓ up to a constant (2·10⁻⁷)
    amp = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(np.tile([1.0, 0.0, 0.0], 3)), G9)[0]
    assert np.abs(flow(amp) - flow(iso)).max() < 1e-5 * np.abs(flow(iso)).max()
    cs = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(np.tile([0.0, 0.25, 0.25], 3)), G9)[0]
    assert np.abs(flow(cs)).max() < 1e-5 * np.abs(flow(iso)).max()
