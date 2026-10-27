"""Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise (src/core/mass_dictionary.py)."""
import numpy as np

from src.core import mass_dictionary as MD
from src.core import gauge_drift as GDR

ONE = lambda k: np.ones_like(k)


def _instance(N=3, d=6, seed=0, balanced=True, m=None):
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((N, d)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(d); ws[[1, 4]] = [0.2, -0.15]
    X = MD.MassDictionary(Phi, Phi @ ws, float(d if m is None else m), balanced=balanced)
    beta = rng.uniform(0.5, 1.5, d); beta *= X.m / beta.sum()
    th = X.retract(X.theta(0.3 * rng.standard_normal(d), beta))
    return X, th, rng


def test_on_the_manifold_and_the_relevance():
    for bal in (True, False):
        X, th, _ = _instance(balanced=bal)
        assert np.abs(X.constraints(th)).max() < 1e-12
        al, rho, b, a, da, dda = X.parts(th)
        assert np.allclose(X.w(th), b * a / X.m)
        J = X.cjac(th)[:X.N]                                                  # the fit's Jacobian in θ, mobility-weighted
        Kd = (X.m / X.N) * (X.Phi * X.relevance(th)[None, :]) @ X.Phi.T
        assert np.allclose(Kd, (X.m / X.N) * (J * X.Pm[None, :]) @ J.T)      # K = (m/N)·J P Jᵀ, D_i its diagonal weights
        if bal:                                                               # the conserved mass removes the direction y
            assert np.allclose(X.gram(th), Kd - (X.gamma / X.N) * np.outer(X.y, X.y))
            assert np.linalg.eigvalsh(X.gram(th)).min() > -1e-12


def test_the_normal_rates_are_the_spectrum_of_K():
    """K's nonzero eigenvalues are the flow's normal rates at a point of M (the flat coordinates, the 0-homogeneous extension)."""
    for bal in (True, False):
        X, th, _ = _instance(balanced=bal)
        d, N = X.d, X.N

        def F(x):
            _, rho, b, a, da, _ = X.parts(x)
            g = X.Phi.T @ (X.y - X.Phi @ (b * a / X.m)) / N
            u = a * g
            lam = b @ u / b.sum() if bal else 0.0
            return np.r_[X.G * b * da * g, 0.5 * rho * X.gamma * (u - lam)]
        A = GDR._jac(F, th, 1e-6)
        ev = np.sort(-np.linalg.eigvals(A).real)[::-1][:N]
        assert np.allclose(ev, np.sort(np.linalg.eigvalsh(X.gram(th)))[::-1], rtol=1e-5, atol=1e-9)


def test_the_closed_forms_are_katzenbergers_drift():
    """The selection against `gauge_drift.limit_drift` on a chart of M: the flat scheme (radially renormalised), the paradigm's (Euler on
    log β, renormalised) with its Itô term; fresh draws and a slow AR(1) at a finite step; conserved and free masses."""
    for bal in (True, False):
        X, th, _ = _instance(balanced=bal)
        B = np.argsort(-np.abs(X.parts(th)[3]))[:X.N]
        for coords in ('flat', 'log'):
            g, _ = X.scheme_drift(th, B, coords=coords)
            assert np.abs(g - X.drift(th, ONE, ito=coords == 'log') * X.m / 4).max() < 1e-4 * np.abs(g).max()
        h, corr = 0.5, ('ar1', float(np.exp(-0.5 / 200.0)))
        filt = lambda k: MD.LN.lr_factor(k, h, corr)
        g, Sf = X.scheme_drift(th, B, h=h, corr=corr, coords='log')
        assert np.abs(g - X.drift(th, filt, ito=True) * X.m / 4).max() < 1e-4 * np.abs(g).max()        # the filtered Itô term's
        assert np.allclose(X.ito_term(th, filt=filt), X.ito_term(th, np.diag(Sf)[X.d:]), rtol=1e-7, atol=0)   # closed form


def test_the_excess_mass_bound():
    X, th, rng = _instance(balanced=True)
    for _ in range(20):
        beta = rng.uniform(0.2, 2.0, X.d); beta *= X.m / beta.sum()
        t = X.retract(X.theta(0.5 * rng.standard_normal(X.d), beta))
        if np.abs(X.constraints(t)).max() > 1e-10:
            continue
        al, rho, b, a, da, dda = X.parts(t)
        E = X.excess(t)
        assert E >= 0 and np.sum(b * da) >= E - 1e-12                        # Σβa′ ≥ E
        assert X.potential(t, lambda k: k) + (X.gamma / X.N) * X.y @ X.y >= X.white_bound(E) - 1e-12


def test_the_gauge_minimal_relevance():
    X, _, _ = _instance(balanced=False, m=12.0)
    big = np.array([0.1, 0.3]); thr = X.gamma * X.a_max ** 3 / (X.G * X.m)
    assert np.allclose(X.relevance_min(big), X.gamma * X.a_max * big / X.m, rtol=2e-2)         # ℓ¹ above the threshold
    small = np.array([1e-6, 1e-5]) * thr
    law = 1.5 * (2 * X.G) ** (1 / 3) * (X.gamma / X.m) ** (2 / 3) * small ** (4 / 3)
    assert np.allclose(X.relevance_min(small), law, rtol=2e-2)                                   # the 4/3 law below it


def test_the_null_kernel_holds_mass_without_cost():
    X, th, rng = _instance(balanced=True)
    Phi0 = np.c_[X.Phi, np.zeros(X.N)]
    Y = MD.MassDictionary(Phi0, X.y, X.m, balanced=True)
    b = X.parts(th)[2]
    for b0 in (0.5, 2.0):
        beta = np.r_[b * (X.m - b0) / X.m, b0]
        t = Y.retract(Y.theta(np.r_[th[:X.d], 0.0], beta))
        g = Y.grad(t, ONE)
        assert abs(g[X.d]) < 1e-14 and abs(g[-1]) < 1e-14                   # R does not depend on the null kernel's state


def test_the_noisy_scheme_keeps_the_mass_and_fits():
    X, th, rng = _instance(balanced=True)
    al, rho, b, a, da, dda = X.parts(th)
    rec, al1, b1 = X.noisy(al, b, 0.5, 0.0, 0.0, 200, rng, 3, 100)
    assert np.allclose(b1.sum(axis=1), X.m)
    assert np.abs(rec[-1] - X.w(th)[None, :]).max() < 1e-10                   # without noise a point of M stays put
