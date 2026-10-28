"""Proposition D‴.44: the replicator's own masses on a kernel dictionary under measurement noise (src/core/mass_dictionary.py)."""
import os
import sys

import numpy as np

from src.core import mass_dictionary as MD
from src.core import gauge_drift as GDR

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))

ONE = lambda k: np.ones_like(k)
SLOW = lambda k: 1.0 / (1.0 + 300.0 * k)


def _instance(N=3, d=6, seed=0, balanced=True, m=None, a_max=1.0, unit=True):
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((N, d)); Phi /= np.linalg.norm(Phi, axis=0)
    if not unit:
        Phi *= rng.uniform(0.5, 1.6, d)[None, :]                                # columns of unequal norms
    ws = np.zeros(d); ws[[1, 4]] = [0.2, -0.15]
    X = MD.MassDictionary(Phi, Phi @ ws, float(d if m is None else m), a_max=a_max, balanced=balanced)
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
    """tr K (the conserved law's K, with its −(γ/N)yyᵀ) ≥ (G/(Nm))E²/Σ|φ_i|⁻² at random exact fits — also near saturation, with columns of
    unequal norms and a_max ≠ 1 —, with equality at y = 0, a = 0, β_i ∝ |φ_i|⁻² (both Cauchy–Schwarz steps tight)."""
    for am in (1.0, 2.0):
        X, th, rng = _instance(balanced=True, a_max=am, unit=False, d=7)
        n = 0
        for scale in (0.5, 3.0):                                              # moderate and saturated amplitudes
            for _ in range(20):
                beta = rng.uniform(0.2, 2.0, X.d); beta *= X.m / beta.sum()
                t = X.retract(X.theta(scale * am * rng.standard_normal(X.d), beta))
                if np.abs(X.constraints(t)).max() > 1e-10:
                    continue
                n += 1
                al, rho, b, a, da, dda = X.parts(t)
                E = X.excess(t)
                assert np.isclose(E, X.m * (1 - np.abs(X.w(t)).sum() / am), rtol=1e-14, atol=0)
                assert E >= 0 and np.sum(b * da) >= E - 1e-12                    # Σβa′ ≥ E
                assert X.potential(t, lambda k: k) >= X.white_bound(E) - 1e-12
        assert n >= 20
    X, _, _ = _instance(balanced=True, unit=False, d=7, a_max=2.0)
    Y = MD.MassDictionary(X.Phi, np.zeros(X.N), X.m, a_max=2.0, balanced=True)
    c2 = (Y.Phi ** 2).sum(axis=0)
    beta = Y.m * (1 / c2) / np.sum(1 / c2)
    t0 = Y.theta(np.zeros(Y.d), beta)
    assert np.abs(Y.constraints(t0)).max() < 1e-13 and np.isclose(Y.excess(t0), Y.m, rtol=0, atol=1e-13)
    assert np.isclose(Y.potential(t0, lambda k: k), Y.white_bound(Y.m), rtol=1e-12, atol=0)


def test_the_gauge_minimal_relevance():
    """D* against its two laws (atol = 0: D* is 10⁻¹² there) and its minimiser: D at the root of 2Gm|w|(1 − t⁴) = γa_max³t³ is D*."""
    for G in (1.0, 2.0):
        X = MD.MassDictionary(np.eye(2, 3) + 0.1, np.ones(2), 12.0, G=G, balanced=False)
        big = np.array([0.1, 0.3]); thr = X.gamma * X.a_max ** 3 / (X.G * X.m)
        assert np.allclose(X.relevance_min(big), X.gamma * X.a_max * big / X.m, rtol=2e-2, atol=0)         # ℓ¹ above the threshold
        small = np.array([1e-6, 1e-5]) * thr
        law = 1.5 * (2 * X.G) ** (1 / 3) * (X.gamma / X.m) ** (2 / 3) * small ** (4 / 3)
        assert np.allclose(X.relevance_min(small), law, rtol=2e-2, atol=0)                                    # the 4/3 law below it
        xs = np.array([1e-5, 1e-3, 0.03, 0.3])
        t = X.gauge_t(xs, paradigm=False)
        assert np.allclose(2 * X.G * X.m * xs * (1 - t ** 4), X.gamma * X.a_max ** 3 * t ** 3, rtol=1e-12, atol=0)
        D_at = X.G * xs ** 2 * (1 - t ** 2) ** 2 / (X.a_max ** 2 * t ** 2) + X.gamma * xs * X.a_max * t / X.m
        assert np.allclose(D_at, X.relevance_min(xs), rtol=1e-9, atol=0)
        h = 1e-4 * xs                                                             # D*′ = (2γa_max/m)t/(1 + t²): increasing, D* convex
        fd = (X.relevance_min(xs + h) - X.relevance_min(xs - h)) / (2 * h)
        assert np.allclose(fd, 2 * X.gamma * X.a_max / X.m * t / (1 + t * t), rtol=1e-7, atol=0)


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


def test_the_null_kernels_potential_is_the_free_ones_only_for_white_noise():
    """K_null = K_free − (γ/N)yyᵀ: R_null − R_free = −(γ/N)|y|² for white noise at every point, (1/τ)log(1 − τ(γ/N)yᵀ(I + τK_free)⁻¹y) —
    not a constant — for Ψ = log(1 + τκ)/τ (referee round 34)."""
    X, _, rng = _instance(balanced=False, d=6)
    Y = MD.MassDictionary(np.c_[X.Phi, np.zeros(X.N)], X.y, X.m, balanced=True)
    tau = 1e4
    slow = lambda k: np.log1p(tau * k) / tau
    vals = []
    for _ in range(5):
        beta = rng.uniform(0.5, 1.5, Y.d); beta *= Y.m / beta.sum()
        t = Y.retract(Y.theta(0.3 * rng.standard_normal(Y.d), beta))
        tr = np.r_[t[:X.d], t[Y.d:Y.d + X.d]]
        assert np.abs(X.constraints(tr)).max() < 1e-12
        assert np.isclose(Y.potential(t, lambda k: k) - X.potential(tr, lambda k: k), -X.gamma / X.N * X.y @ X.y, rtol=1e-10, atol=0)
        Kf = X.gram(tr)
        lemma = np.log(1 - tau * X.gamma / X.N * X.y @ np.linalg.solve(np.eye(X.N) + tau * Kf, X.y)) / tau
        dv = Y.potential(t, slow) - X.potential(tr, slow)
        assert np.isclose(dv, lemma, rtol=1e-8, atol=0)
        vals.append(dv)
    assert np.ptp(vals) > 1e-3 * np.abs(vals).max()


def test_the_noisy_scheme_keeps_the_mass_and_fits():
    X, th, rng = _instance(balanced=True)
    al, rho, b, a, da, dda = X.parts(th)
    rec, al1, b1 = X.noisy(al, b, 0.5, 0.0, 0.0, 200, rng, 3, 100)
    assert np.allclose(b1.sum(axis=1), X.m)
    assert np.abs(rec[-1] - X.w(th)[None, :]).max() < 1e-10                   # without noise a point of M stays put


def test_one_noiseless_step_is_the_law():
    """Off M, without noise: α₁ = α + hGβa′g, log β₁ = log β + hγ(u − λ) (λ = Σβu/Σβ, then renormalised to Σβ = m), g = Φᵀ(y − Φw)/N."""
    for bal in (True, False):
        X, th, rng = _instance(balanced=bal)
        al = th[:X.d] + 0.2 * rng.standard_normal(X.d)
        beta = X.parts(th)[2] * np.exp(0.1 * rng.standard_normal(X.d))
        if bal:
            beta *= X.m / beta.sum()
        h = 0.3
        _, al1, b1, kept = X.noisy(al, beta, h, 0.0, 0.0, 1, rng, 2, 1, keep=(1,))
        a = X.a_max * np.tanh(al / X.a_max); da = 1 - (a / X.a_max) ** 2
        g = X.Phi.T @ (X.y - X.Phi @ (beta * a / X.m)) / X.N
        u = a * g
        lb = np.log(beta) + h * X.gamma * (u - (beta @ u / beta.sum() if bal else 0.0))
        if bal:
            lb -= np.log(np.exp(lb).sum() / X.m)
        assert np.allclose(al1[0], al + h * X.G * beta * da * g, rtol=1e-13, atol=1e-15)
        assert np.allclose(np.log(b1[0]), lb, rtol=1e-13, atol=1e-14)
        assert np.array_equal(kept[0][0], al1) and np.array_equal(kept[0][1], b1)


def test_one_noisy_step_moments_and_the_ar1_noise():
    """From a point of M (free masses) with a short step the increments are the noise's: covariance (hσ/N)²·VΦᵀΦVᵀ per atom pair (V the
    fields Gβa′ in α, γa in log β), mean zero; the AR(1) draws stay stationary (the increments' variance at steps 1 and 40 alike) with
    lag-one correlation ρ."""
    X, th, rng = _instance(balanced=False)
    al, rho_, b, a, da, _ = X.parts(th)
    h, sig, P = 1e-3, 0.5, 40000
    for corr in (0.0, 0.6):
        _, _, _, kept = X.noisy(al, b, h, sig, corr, 41, np.random.default_rng(7), P, 41, keep=(0, 1, 2, 40, 41))
        z = [np.c_[k_[0], np.log(k_[1])] for k_ in kept]                       # (P, 2d) states (α, log β)
        d1, d2, d40 = z[1] - z[0], z[2] - z[1], z[4] - z[3]
        V = np.r_[X.G * b * da, X.gamma * a]
        Cth = (h * sig / X.N) ** 2 * np.outer(V, V) * np.tile(X.Phi.T @ X.Phi, (2, 2))
        C1 = np.cov(d1.T)
        assert np.abs(C1 - Cth).max() < 0.05 * np.abs(Cth).max()
        assert np.abs(d1.mean(axis=0)).max() < 5 * np.sqrt(np.diag(Cth).max() / P)
        assert np.abs(np.var(d40, axis=0) / np.var(d1, axis=0) - 1).max() < 0.05
        c12 = np.array([np.corrcoef(d1[:, j], d2[:, j])[0, 1] for j in range(2 * X.d)])
        assert np.abs(c12 - corr).max() < 0.03


def test_the_gauge_invariant_drift_closed_form():
    """Free masses: ∇Q·(the selection's drift) is (ω_iγβ_i/(Nm))(κγa_i² − 2Gβ_i(1 − t_i⁴)) — κ = 1 flat, 2 the paradigm's — for white and
    slow noise, also against Katzenberger's drift of the scheme; Q has no quadratic variation (the noise moves an atom along its curve)."""
    X, th, _ = _instance(balanced=False, unit=False, a_max=1.3)
    al, rho, b, a, da, _ = X.parts(th)
    c = X.gamma * X.a_max ** 2 / (4 * X.G)
    gQ = lambda q: -c * (2 / X.a_max) * np.sinh(2 * al / X.a_max) * q[:X.d] + 2 * rho * q[X.d:]
    for filt in (ONE, SLOW):
        for ito in (False, True):
            qd = X.q_drift(th, filt, ito)
            assert np.abs(gQ(X.drift(th, filt, ito)) - qd).max() < 1e-12 * np.abs(qd).max()
    B = np.argsort(-np.abs(a))[:X.N]
    g, Sf = X.scheme_drift(th, B, coords='log')
    assert np.abs(gQ(g) - X.q_drift(th, ONE, True) * X.m / 4).max() < 1e-5 * np.abs(X.q_drift(th, ONE, True) * X.m / 4).max()
    gz = np.r_[-c * (2 / X.a_max) * np.sinh(2 * al / X.a_max), b]
    for i in range(X.d):
        j = [i, X.d + i]
        assert abs(gz[j] @ Sf[np.ix_(j, j)] @ gz[j]) < 1e-12 * (gz[j] @ gz[j]) * np.abs(Sf[np.ix_(j, j)]).max()
    assert np.allclose(X.gauge_invariant(th), b - c * np.cosh(2 * al / X.a_max), rtol=1e-15, atol=0)


def test_the_paradigms_free_selection_is_the_euler_lag():
    X, th, _ = _instance(balanced=False, unit=False)
    for filt in (ONE, SLOW):
        full = X.drift(th, filt, ito=True)
        assert np.abs(X.project(th, X.euler_lag(th, filt)) - full).max() < 1e-12 * np.abs(full).max()


def test_the_gauge_equilibria_and_the_mirror_rate():
    """At the paradigm's gauge equilibrium (Gm|w|(1 − t⁴) = γa_max³t³) Q's paradigm drift vanishes, at D*'s split the flat one; there the
    mirror coordinate η = (α + a_maxκ′arctan(κ′t))/(2Gc), κ′² = (c − Q)/(c + Q), drifts by ω·`mirror_rate_eq`(|w|)."""
    X = MD.MassDictionary(np.array([[1.0, 0.3], [0.2, 1.0]]), np.array([0.1, 0.2]), 12.0, a_max=1.0, balanced=False)
    G, g, am, m, N = X.G, X.gamma, X.a_max, X.m, X.N
    c = g * am ** 2 / (4 * G)
    for x in (1e-3, 0.05, 0.3):
        for par in (True, False):
            t = X.gauge_t([x, x], paradigm=par)
            th = X.theta(am * np.arctanh(t), m * x / (am * t))
            assert np.allclose(X.w(th), x, rtol=1e-12)
            assert np.allclose(X.gauge_ratio(th), 1.0 if par else 0.5, rtol=1e-10)
            assert np.abs(X.q_drift(th, ONE, ito=par)).max() < 1e-12 * np.abs(X.q_drift(th, ONE, ito=not par)).max()
        t = X.gauge_t([x], paradigm=True)[0]

        def eta(al, lb):
            b = np.exp(lb); Q = b - c * np.cosh(2 * al / am); kp = np.sqrt((c - Q) / (c + Q))
            return (al + am * kp * np.arctan(kp * np.tanh(al / am))) / (2 * G * c)
        al0, lb0 = am * np.arctanh(t), np.log(m * x / (am * t))
        th = X.theta([al0, al0], np.exp([lb0, lb0]))
        lag = X.euler_lag(th, ONE) / X.weights(th, ONE)[0]                       # per unit ω, in θ
        u = np.array([lag[0], 2 * lag[2] / np.sqrt(np.exp(lb0))])                # (α, log β) components of atom 0
        hh = 1e-6
        grad = np.array([(eta(al0 + hh, lb0) - eta(al0 - hh, lb0)) / (2 * hh), (eta(al0, lb0 + hh) - eta(al0, lb0 - hh)) / (2 * hh)])
        assert np.isclose(grad @ u, X.mirror_rate_eq([x])[0], rtol=1e-5)
    big = np.array([0.3, 1.0])
    assert np.allclose(X.mirror_rate_eq(big), -g ** 2 * am ** 4 / (2 * N * G * m ** 2 * big), rtol=0.05)   # the logarithmic law


def test_the_paradigms_drift_is_not_a_gradient():
    """Loop integrals of P⁻¹·drift on M: the flat (white) drift's vanish (a gradient in the flow's metric), the paradigm's scale as ε²."""
    for bal in (False, True):
        X, th, rng = _instance(balanced=bal, d=6)
        q1, q2 = X.project(th, rng.standard_normal(2 * X.d)), X.project(th, rng.standard_normal(2 * X.d))
        ip = lambda x, y: x @ (y / X.Pm)
        u = q1 / np.sqrt(ip(q1, q1)); q2 = q2 - ip(q2, u) * u; v = q2 / np.sqrt(ip(q2, q2))

        def loop(eps, ito, n=48):
            ph = 2 * np.pi * np.arange(n) / n
            P_ = np.array([X.retract(th + eps * (np.cos(p) * u + np.sin(p) * v)) for p in ph])
            k = np.fft.fftfreq(n, 1.0 / n)
            dP = np.real(np.fft.ifft(1j * k[:, None] * np.fft.fft(P_, axis=0), axis=0))
            return sum((X.drift(p, ONE, ito=ito) / X.Pm) @ dp for p, dp in zip(P_, dP)) * 2 * np.pi / n
        I1, I2 = loop(1e-2, True), loop(5e-3, True)
        assert abs(I1 / I2 - 4) < 0.05                                            # ε²: a nonzero 2-form
        assert abs(loop(1e-2, False)) < 1e-6 * abs(I1)


def test_the_conserved_leftover_is_the_cross_term():
    X, th, _ = _instance(balanced=True, unit=False)
    import exp_mass_dictionary as M
    for filt in (ONE, SLOW):
        left = X.drift(th, filt, ito=True) - M.own_removed(X, th, filt)
        cr = X.project(th, X.cross_term(th, filt) * 4 / X.m)
        assert np.abs(left - cr).max() < 1e-10 * np.abs(left).max()


def test_summary_on_hand_built_states():
    """Used: |w| > 10⁻³ (an atom at 5·10⁻³ counts), alive: β > 10⁻³ (an atom at 5·10⁻³ counts), the excess over the real atoms."""
    import exp_mass_dictionary as M
    Phi = np.array([[1.0, 0.0, 0.6, 0.3], [0.0, 1.0, 0.8, -0.5]])
    ws = np.array([0.3, -0.2, 0.0, 0.005])
    Y = MD.MassDictionary(np.c_[Phi, np.zeros(2)], Phi @ ws, 6.0, a_max=2.0, balanced=True)
    t1 = np.array([0.5, -0.4, 0.0005, 0.3, 0.7])                                  # amplitudes t = a/a_max (the null atom's 0.7)
    beta = np.r_[6.0 * 0.3 / (2.0 * 0.5), 6.0 * 0.2 / (2.0 * 0.4), 5e-3, 0.05, 0.0]
    beta[4] = 6.0 - beta[:4].sum()
    th = Y.theta(2.0 * np.arctanh(t1), beta)
    wbp = np.array([0.3, -0.2, 0.0, 0.0])
    s = M.summary(Y, th, ws, wbp, 4)
    w = beta[:4] * 2.0 * t1[:4] / 6.0
    assert np.allclose(s['w'], w) and np.isclose(s['l1'], np.abs(w).sum() / 2.0) and np.isclose(w[3], 0.005)
    assert np.isclose(s['exc'], 1 - s['l1'])                                      # over the real atoms: the null atom's |w_0| left out
    assert s['used'] == 3 and s['alive'] == 4 and s['ok'] and np.isclose(s['b0'], beta[4]) and np.isclose(s['amax'], 0.5)
    assert np.isclose(s['d_bp'], 0.005) and s['gr'] is None
    X = MD.MassDictionary(Phi, Phi @ ws, 6.0, a_max=2.0, balanced=False)
    th2 = X.theta(2.0 * np.arctanh(t1[:4]), beta[:4])
    s2 = M.summary(X, th2, ws + np.r_[0, 0, 1e-3, 0], wbp, 4)
    assert np.isnan(s2['exc']) and not s2['ok'] and s2['gr'].shape == (4,)
    assert np.allclose(s2['gr'], X.gauge_ratio(th2))
    assert np.allclose(M.filt_of(0.0)(np.array([0.0, 3.0])), 1.0)
    assert np.allclose(M.filt_of(100.0)(np.array([0.0, 0.03])), [1.0, 0.25])


def test_report_rows_on_a_synthetic_result():
    """Part C's table on two hand-made instances: one row per (mode, colour) over the finished flows, the fractions split by bp_ok."""
    import exp_mass_dictionary as M

    def sm(l1, ok, d_bp, used=2, alive=3, w=(0.1, 0.0), b0=np.nan, exc=np.nan, gr=None):
        return dict(l1=l1, exc=exc, used=used, alive=alive, ok=ok, d_truth=0.0 if ok else 0.1, d_bp=d_bp, b0=b0, amax=0.9,
                    w=np.array(w), gr=gr)
    R = []
    for i, bp_ok, d_bp in ((0, True, 0.0), (1, False, 0.2), (2, True, 0.05)):
        r = dict(kind='gauss', k=2, i=i, bp_ok=bp_ok, l1_truth=0.4, l1_bp=0.4 + 0.1 * i, dstar_vs_flat=1e-4, dstar_vs_bp=1e-3)
        for mode in M.MODES:
            for key in ['start'] + list(M.TAUS) + (['flat'] if mode != 'null kernel' else []):
                r[(mode, key)] = sm(0.5 + 0.1 * i, i == 0, d_bp, exc=0.5 if mode != 'free' else np.nan,
                                    b0=4.0 if mode == 'null kernel' else np.nan, gr=np.array([1.0, 1.0]) if mode == 'free' else None)
        R.append(r)
    R[1][('conserved', 0.0)] = dict(stopped='over the limit', s=2.5e6, last=sm(0.9, False, 0.3))
    md = M.report(C=R)
    assert '| 2 | — | basis pursuit | 3 | 0.500 (0.400) | | | | 2/3 | | | |' in md
    assert '| 2 | conserved | white | 2 | 0.600 (0.400) | 0.500 | 2.0 (2) | 3.0 | 1/2, 0/0 | 1/2, 0/0 | 0.050 | — |' in md
    assert '| 2 | free | white | 3 | 0.600 (0.400) | — | 2.0 (2) | 3.0 | 1/2, 0/1 | 1/2, 0/1 | 0.200 | — |' in md
    assert 'k = 2, instance 1, conserved, 0.0: over the limit at s = 2.5e+06' in md
    assert '| 2 | white | 0 of 3 (—) | 0 of 1 (—) |' in md


def test_the_flow_keeps_a_dropped_atom_at_the_floor():
    """The retraction's correction is additive in β, so on a dropped atom it can cross zero; the flow carries such an atom at the floor
    e^{−700} into the next piece instead of a log of a negative mass."""
    X, th, _ = _instance(balanced=False)
    real, calls = X.retract, [0]

    def crossing(t, **kw):
        t = real(t, **kw).copy()
        calls[0] += 1
        if calls[0] > 1:                                                          # after the start's retraction
            t[X.d] = -1e-300                                                      # atom 0's √β pushed just below zero
        return t
    X.retract = crossing
    s, Th = X.flow(th, 5.0, ONE, chunks=3)
    assert np.all(np.isfinite(Th)) and len(s) == 4 and calls[0] > 4


def test_settle_is_the_loops_noiseless_law():
    """`settle`: the loop's own law in (α, log β) — the replicator's centring for the conserved mass —, retracted onto M."""
    from scipy.integrate import solve_ivp
    for bal in (True, False):
        X, th, rng = _instance(balanced=bal)
        al = th[:X.d] + 0.3 * rng.standard_normal(X.d)
        beta = X.parts(th)[2] * np.exp(0.2 * rng.standard_normal(X.d))
        if bal:
            beta *= X.m / beta.sum()

        def F(t, z):
            a = X.a_max * np.tanh(z[:X.d] / X.a_max); b = np.exp(z[X.d:]); da = 1 - (a / X.a_max) ** 2
            g = X.Phi.T @ (X.y - X.Phi @ (b * a / X.m)) / X.N
            u = a * g
            return np.r_[X.G * b * da * g, X.gamma * (u - (b @ u / b.sum() if bal else 0.0))]
        z = solve_ivp(F, (0.0, 50.0), np.r_[al, np.log(beta)], method='DOP853', rtol=1e-12, atol=1e-14).y[:, -1]
        ref = X.retract(np.r_[z[:X.d], np.exp(z[X.d:] / 2)])
        assert np.abs(X.settle(al, beta, T=50.0) - ref).max() < 1e-7


def test_projection_is_accurate_with_nearly_dead_atoms_and_the_flow_finishes():
    # masses down to e^{−60}: the normal equations J P Jᵀ are singular to rounding there (their condition is A's squared; the old
    # projection's idempotency error 1.8·10⁻⁷ at this point, the new one's 1.3·10⁻¹¹), and the
    # free masses' flat flow on this instance crawled for hours; on A = P^{1/2}Jᵀ the projection is a projection and the flow finishes
    import time
    import exp_mass_dictionary as M
    Phi, ws, al0 = M.instance(2, 1, 'gauss')
    X, _ = M.law(Phi, ws, 'free')
    rng = np.random.default_rng(3)
    th = np.r_[rng.normal(size=X.d), np.exp(-np.linspace(0.0, 30.0, X.d))]
    q = rng.normal(size=2 * X.d)
    p = X.project(th, q)
    assert np.abs(X.project(th, p) - p).max() < 1e-9 * np.abs(p).max()
    t0 = time.time()
    th_end = X.flow(M.start(X, al0), M.S_END, M.filt_of(0.0), ito=False, chunks=60)[1][-1]
    assert time.time() - t0 < 120
    w = X.w(th_end)
    assert set(np.flatnonzero(np.abs(w) > 1e-6)) == set(np.flatnonzero(ws))


def test_the_flow_stops_at_its_wall_clock_limit_and_on_failure(monkeypatch):
    """`flow(max_secs=…)` raises FlowTimeout once over its limit (checked at every evaluation, so one long piece cannot overrun it), with
    the last state reached; within the limit the flow is the unlimited one, bit for bit; a piece that LSODA and Radau both fail raises
    FlowFailed (solve_ivp's status, referee round 34)."""
    import pytest
    import scipy.integrate
    X, th, _ = _instance()
    with pytest.raises(MD.FlowTimeout) as e:
        X.flow(th, 10.0, ONE, max_secs=0.0)
    assert e.value.s == 0.0 and np.allclose(e.value.th, X.retract(th))
    assert np.array_equal(X.flow(th, 10.0, ONE, max_secs=60.0)[1], X.flow(th, 10.0, ONE)[1])
    real = scipy.integrate.solve_ivp
    calls = []

    def failing(fun, t_span, y0, method='RK45', **kw):
        sol = real(fun, t_span, y0, method=method, **kw)
        calls.append(method)
        if t_span[0] > 0:
            sol.success, sol.message = False, 'mock failure'
        return sol
    monkeypatch.setattr(scipy.integrate, 'solve_ivp', failing)
    with pytest.raises(MD.FlowFailed) as e:
        X.flow(th, 10.0, ONE, chunks=4)
    assert calls == ['LSODA', 'LSODA', 'Radau'] and e.value.s == 2.5 and np.all(np.isfinite(e.value.th))
