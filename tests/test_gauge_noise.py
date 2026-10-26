"""
Proposition D‴.40 (PARADIGM.md §5d): the gauge under general noise. Checked here: the noise mask of the dithered integrator dithers only
the chosen coordinates (the amplitudes then move by their drift alone).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
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
