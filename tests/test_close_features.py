"""Proposition D‴.52: close features — the second entry's delay ∝ Δ⁻² (a dipole's ∝ Δ⁻¹, the transit of its derivative shape), the
two-kernel block's stiffnesses ∝ Δ² and Δ⁶ in closed form, and the fast flow's rates along those modes."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import exp_close_features as CF
import exp_mirror_transit as X
from src.core import rich_transit as RT


def test_block_stiffnesses_are_the_closed_forms():
    """At the exact fit the Gauss–Newton Hessian's two small eigenvalues are λ_split = Δ⁶⟨(P⊥e₃)²⟩/8 — the same for every pair of
    amplitudes — and λ_sep = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²), to O(Δ)."""
    L, K, c0 = CF.setup()
    D = 0.0125
    split = []
    for A1, A2 in ((0.8, 0.6), (0.8, -0.4), (0.8, 0.8)):
        r = CF.b_job((A1, A2, D))
        l1, l2 = r['pred']
        assert abs(r['hess'][0] / l1 - 1) < 5e-3
        assert abs(r['hess'][1] / l2 - 1) < 5e-2
        assert r['hess'][2] > 10 * r['hess'][1]                           # the other two are of order one
        split.append(r['hess'][0] / D ** 6)
    assert np.ptp(split) / np.mean(split) < 1e-3                         # amplitude-free
    vs = CF.taylor_basis(K, CF.CBAR)                                     # the continuum's √π/(48σ⁵ℓ) and √π/(8σ³ℓ), ℓ = 2, to the grid
    s = K.sig
    assert abs(CF.perp_norm2(vs, 3) / (np.sqrt(np.pi) / (48 * s ** 5 * 2)) - 1) < 2e-2
    assert abs(CF.perp_norm2(vs, 2) / (np.sqrt(np.pi) / (8 * s ** 3 * 2)) - 1) < 2e-2


def test_flow_relaxes_the_split_at_the_mobility_weighted_stiffness():
    """The fast flow's slowest mode is the split, at λ_split·2/(1/μ₁ + 1/μ₂) (μ_k the amplitudes' mobility); the separation relaxes at
    Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) while the amplitudes are frozen on its scale; both modes identified by their eigenvectors."""
    r = CF.b_job((0.8, 0.6, 0.0125))
    (rs, os_), (rp, op) = r['split'], r['sep']
    ps, pp = r['pred_flow']
    assert os_ > 0.999 and op > 0.999
    assert abs(rs / ps - 1) < 1e-2
    assert abs(rp / pp - 1) < 5e-2
    assert rs == min(r['rates'])
    L = X.law()
    mu_branch = 4 * L.G * L.c * np.array([0.8, 0.6]) / L.A                 # μ ≈ 4Gc|w|/a on the exponential branch
    assert np.allclose(r['mu'], mu_branch, rtol=1e-3)


def test_second_entry_delay_is_the_leading_order():
    """A same-sign pair Δ = 0.025 apart: the second entry comes after Δτ₂ = min_j (s_j − u_j)/v_j, v_j = ½V⟨ψφ_j⟩, to 0.5%; the path ends at
    the exact two-kernel fit."""
    r = CF.a_job((0.8, 0.6, 0.025))
    (t1, j1), (t2, j2) = r['entries'][:2]
    dt, j, wall = r['pred']
    assert j2 == j and wall == 1
    assert abs((t2 - t1) / dt - 1) < 5e-3
    assert r['n_act'] == 2 and np.allclose(r['c_end'], r['cs'], atol=1e-6) and r['rms'] < 1e-10


def test_dipole_is_its_derivative_shape_slowed_by_the_separation():
    """A dipole (A₂ = −A₁) Δ = 0.025 apart: τ₁ and Δτ₂ are the derivative shape −A₁∂_cφ's, divided by Δ, to 0.5% and 1.5%, the same
    kernels entering."""
    r = CF.a_job((0.8, -0.8, 0.025))
    (t1, j1), (t2, j2) = r['entries'][:2]
    sh = r['shape']
    assert [j1, j2] == sh['kernels']
    assert abs(t1 / sh['tau1'] - 1) < 5e-3
    assert abs((t2 - t1) / sh['dtau2'] - 1) < 1.5e-2


def test_predicted_delay_scales_as_the_inverse_square_separation():
    """The leading-order delay is ∝ Δ⁻²: halving Δ multiplies it by 4 (to O(Δ))."""
    L, K, c0 = CF.setup()
    d = []
    for D in (0.02, 0.01):
        f = X.bumps_target(K, CF.pair_centres(0.8, 0.6, D), (0.8, 0.6))
        d.append(CF.predicted_delay(K, f, c0, 0.8, 0.6, D)[0])
    assert abs(d[1] / d[0] / 4 - 1) < 2e-2


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn(); print(name, 'ok')
