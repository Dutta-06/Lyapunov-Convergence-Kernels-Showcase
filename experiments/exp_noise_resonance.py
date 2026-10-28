"""
Proposition D‴.60 (planned): measurement noise and the formed kernels — the adaptation's resonance and the noise-driven drift
===========================================================================================================================

D‴.49's closed loop — ẋ = f(x) + u, u = ẋ_m − k_eẽ − f̂(x̃), ẽ = x̃ − x_m, x̃ = x + n; the regressor at x_m; 21 kernels of width 0.12, the
free law at the masses ε = EPS; three features — started from the fitted state (each feature's kernel at its centre with its amplitude,
D‴.57's `fitted_state`) under measurement noise n. Near the fit, y = Σ_i(φ_iδw_i + w_i∂_cφ_iδc_i) (the fit's perturbation at x_m)
obeys y″ + k_eẏ + K(t)y = K(ṅ − f̂′n), K(t) = Σ_i[μ_iφ_i² + v_iw_i(∂_cφ_i)²] at x_m(t), and the parameters jitter by
δw_i = (μ_iφ_i/K)y, δc_i = (v_i∂_cφ_i/K)y — μ_i the coefficient's mobility on its curve, v_i the centre's: Gβ_ia_i = Gmw_i under the
paradigm's metric, G·a_i under the Wasserstein–Fisher–Rao metric.

A. The linear response: one tone of amplitude AP_A at the frequencies OM_A, both metrics, NP_A periods of the reference: per formed kernel
   the rms jitter (the deviation from the moving average over one tone period) of the loop, of the loop linearised at the fit along the
   orbit (`variational`) and of the quasi-static oscillator law (`transfer`), through the same filter; √K's largest value per metric.
B. The second order: one tone at OM_B or D‴.49's 16 tones, of rms s ∈ S_B, both metrics, NP_B periods: the formed kernels' period means
   of w and c (dense outputs) at MARKS_B against the loop averaged over the fast tones (`envelope_means`: the mean state and the
   first-order complex envelopes, exact along the orbit). B2: at one tone the averaged loop without each of its tone averages (the
   control's curvature, its slope moved by the jitter, the blur, the regressors' correlation with the drive) and with quasi-static
   envelopes.
B3. The exact second order: the loop expanded about the fit to second order in the noise along the orbit, without averaging
   (`second_order`), against the loop's even part in the noise (its period means under s·n and −s·n, half-summed: the drift) and its odd
   part (half-differenced), at s ∈ S_B3 and B's tones, both metrics; the averaged loop (B) against it.
B4. The long horizon and commensurability (referee round 56): the loop's even and odd parts per s² and s at s ∈ S_B4 over NP_B4 periods,
   at one tone (21.37; 21.4 = 107ω and 21.3 = 106.5ω, ω the reference's; 13.1 = 65.5ω) and D‴.49's 16 tones: whether the even part
   saturates and whether the odd part stays bounded.
C. The leaving levels: the levels of D‴.49's noise at which a formed kernel leaves its feature within NP_C periods (D‴.57's Part D),
   the paradigm's metric, against the error gain k_e (KE_C) and the law's m (M_C: the centres' mobility Gmw), under the noise's seeds
   SEEDS_C (the tones' frequencies and phases) — with √K_max's place against the noise band 10–30 rad/s. C2: controls at k_e = 10 —
   single tones, the band moved to 25–75 rad/s (above √K_max), m between 3 and 12, the Wasserstein–Fisher–Rao metric, tone pairs whose
   difference is commensurate with the reference or not (referee round 56).

Usage: python experiments/exp_noise_resonance.py [--only A,B,B2,B3,B4,C,C2] [--from-raw] → results/chaos/noise_resonance.md (+ _raw.npz).
"""
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import exp_closed_transit as CT
import exp_cohesion_leak as E
import exp_dead_zone as DZ
from src.core import mirror_free as MF

OUT = os.environ.get('NR_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('NR_WORKERS', '2'))
EPS = DZ.EPS_C
METRICS = ('paradigm', 'wfr')
AP_A, NP_A, PH_A = 2e-3, 4, 0.3
OM_A = (4.13, 7.3, 10.7, 13.9, 16.3, 19.1, 23.7, 29.3, 37.1)
S_B, NP_B, MARKS_B, OM_B = (0.0025, 0.005), 30, (3, 10, 30), (13.1, 21.37, 29.3)
B2_CASE = ('paradigm', (21.37,), 0.005)
KE_C, M_C, NP_C = (5.0, 10.0, 20.0, 40.0), (3.0, 12.0), 300
SEEDS_C = (7, 1, 2, 3)                                                 # C: the noise's seeds (7 is D‴.49's), referee round 56
C2_CASES = ([('one tone', 12.0, 10.0, s_, dict(freqs=(f,))) for f in (13.1, 15.7, 21.3) for s_ in (0.02, 0.04)]
            + [('band 25–75', m_, 10.0, s_, dict(band=(25.0, 75.0))) for m_ in (3.0, 12.0) for s_ in (0.02, 0.04, 0.06, 0.08)]
            + [('m', m_, 10.0, s_, {}) for m_ in (4.8, 6.0, 8.0) for s_ in (0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08)]
            + [('WFR', 12.0, ke_, s_, dict(metric='wfr')) for ke_ in (5.0, 10.0) for s_ in (0.02, 0.04, 0.08)]
            + [('tone pair', 12.0, 10.0, 0.02, dict(freqs=fr, seed=sd)) for fr in ((18.9, 22.5), (18.9, 22.4)) for sd in SEEDS_C])
NP_B4, MARKS_B4, S_B4 = 300, (1, 3, 10, 30, 60, 100, 150, 200, 250, 300), (0.00125, 0.0025)
FREQS_B4 = ((21.37,), None, (21.4,), (21.3,), (13.1,))
LEV_C = (0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08)
LEFT_C = 0.3
TERMS = ('ctl', 'cross', 'blur', 'corr')
S_B3 = (0.00125, 0.0025, 0.005)
TERM_NAMES = dict(ctl="the control's curvature", cross="the control's slope moved by the jitter", blur='the blur',
                  corr="the regressors' correlation with the drive")


def law(m=None):
    return MF.FreeLaw(**dict(CT.LAW, **({} if m is None else dict(m=m))))


def mobilities(L, eps, al, metric):
    """On the curves from the masses eps: w, the coefficient's mobility μ (ẇ = μẽφ) and the centre's v (ċ = vẽ∂_cφ)."""
    th, sh = np.tanh(al / L.A), np.sinh(al / L.A)
    beta = eps + 2 * L.c * sh ** 2
    a, da = L.A * th, 1 - th ** 2
    w = beta * a / L.m
    dwda = ((2 * L.c / L.A) * np.sinh(2 * al / L.A) * a + beta * da) / L.m
    return w, L.G * beta * da * dwda, (L.G * beta * a if metric == 'paradigm' else L.G * a)


def orbit_terms(L, eps, al, c, metric, xm):
    """At the orbit's points xm: the formed kernels' φ_i, ∂_cφ_i, ∂²_cφ_i, the gain K, the mobilities and the fit's f̂′, f̂″ there."""
    w, mu, v = mobilities(L, eps, al, metric)
    u = (xm[:, None] - c[None, :]) / CT.SIG
    ph = np.exp(-u ** 2 / 2)
    d1 = ph * u / CT.SIG
    d2 = ph * (u ** 2 - 1) / CT.SIG ** 2
    return dict(w=w, mu=mu, v=v, ph=ph, d1=d1, d2=d2, K=(mu * ph ** 2 + v * w * d1 ** 2).sum(1),
                f1=-(w * d1).sum(1), f2=(w * d2).sum(1))


def transfer(o, om, ke):
    """Per orbit point and tone: Y = ∫ẽ per unit noise amplitude, (iω − f̂′)/(K − ω² + ik_eω), and ẽ's, iωY."""
    Y = (1j * om[None, :] - o['f1'][:, None]) / (o['K'][:, None] - om[None, :] ** 2 + 1j * ke * om[None, :])
    return Y, 1j * om[None, :] * Y


def simulate(L, al0, c0, astar, ts, noise, metric, ke=None, x_init=None):
    from src.core import closed_loop as CL
    return CL.simulate(L, EPS, c0, CT.SIG, (CT.CSTAR, astar, CT.SIG), CT.KE if ke is None else ke,
                       (0.0, 1.0, CT.OM, 0.0, np.sqrt(2) * CT.OM), ts, noise=noise, th=CT.th_abs(astar), regx=True, al0=al0, metric=metric,
                       x_init=x_init)


def variational(L, eps, al, c, metric, ts, om, amp, ph0, ke):
    """The loop linearised at the exact fit along x_m(t) (the formed kernels; the others' mobilities vanish at α = 0 to O(ε²)), driven by
    the tone amp·cos(ωt + ph0): e′ = ṅ − f′(x_m)n − Σ_i(φ_iδw_i + w_i∂_cφ_iδc_i) − k_e e, δw_i′ = μ_iφ_i e, δc_i′ = v_i∂_cφ_i e, e the
    measured error, starting with the true error zero (e = n, as the loop starts at x = x_m). Returns δw, δc at ts."""
    from scipy.integrate import solve_ivp
    w, mu, v = mobilities(L, eps, al, metric)

    def rhs(t, y):
        o = orbit_terms(L, eps, al, c, metric, np.array([np.sin(CT.OM * t)]))
        n, dn = amp * np.cos(om * t + ph0), -amp * om * np.sin(om * t + ph0)
        e = y[0]
        return np.r_[dn - o['f1'][0] * n - (o['ph'][0] * y[1:1 + c.size] + w * o['d1'][0] * y[1 + c.size:]).sum() - ke * e,
                     mu * o['ph'][0] * e, v * o['d1'][0] * e]
    y0 = np.zeros(1 + 2 * c.size)
    y0[0] = amp * np.cos(om * ts[0] + ph0)
    sol = solve_ivp(rhs, (ts[0], ts[-1]), y0, method='DOP853', t_eval=ts, rtol=1e-10, atol=1e-14, max_step=0.2 / om)
    return sol.y[1:1 + c.size].T, sol.y[1 + c.size:].T


def jitter(z, nwin):
    """The deviation from the moving average over nwin samples (one tone period), the ends dropped."""
    return (z - np.convolve(z, np.ones(nwin) / nwin, mode='same'))[nwin:-nwin]


def a_job(args):
    """A: one tone of amplitude AP_A at ω from the fitted state: per formed kernel the rms jitter of the loop, of its linearisation along
    the orbit (`variational`) and of the quasi-static oscillator law δθ = Re(g·Y·N e^{iωt}) (`transfer`), through the same filter."""
    metric, om = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law()
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    per = 2 * np.pi / CT.OM
    ts = np.linspace(0, NP_A * per, NP_A * 4000 + 1)
    r = simulate(L, al0, c0, CT.ASTAR, ts, dict(w=np.array([om]), ph=np.array([PH_A]), amp=np.array([AP_A])), metric)
    dwv, dcv = variational(L, EPS, al0[idx], c0[idx], metric, ts, om, AP_A, PH_A, CT.KE)
    o = orbit_terms(L, EPS, al0[idx], c0[idx], metric, np.sin(CT.OM * ts))
    YN = transfer(o, np.array([om]), CT.KE)[0][:, 0] * AP_A * np.exp(1j * (om * ts + PH_A))
    dcq, dwq = np.real(o['v'] * o['d1'] * YN[:, None]), np.real(o['mu'] * o['ph'] * YN[:, None])
    nwin = max(3, int(round(2 * np.pi / om / (ts[1] - ts[0]))))
    rms = lambda z: float(np.sqrt(np.mean(jitter(z, nwin) ** 2)))
    out = dict(metric=metric, om=om, sqrtK=float(np.sqrt(o['K'].max())), v=o['v'].tolist(), mu=o['mu'].tolist(),
               sim_c=[rms(r['c'][:, i]) for i in idx], var_c=[rms(dcv[:, q]) for q in range(3)], qs_c=[rms(dcq[:, q]) for q in range(3)],
               sim_w=[rms(r['w'][:, i]) for i in idx], var_w=[rms(dwv[:, q]) for q in range(3)], qs_w=[rms(dwq[:, q]) for q in range(3)],
               secs=time.time() - t0)
    print(f"  A {metric} ω {om}: centres loop/linearised {', '.join(f'{a / b:.4f}' for a, b in zip(out['sim_c'], out['var_c']))}, "
          f"quasi-static/linearised {', '.join(f'{a / b:.3f}' for a, b in zip(out['qs_c'], out['var_c']))} [{out['secs']:.0f}s]", flush=True)
    return out


def alpha_of(L, eps, w, al):
    """α on the curves from the masses eps with w(α) = w, by Newton's method from al (w(α) is increasing)."""
    for _ in range(60):
        th, sh = np.tanh(al / L.A), np.sinh(al / L.A)
        beta = eps + 2 * L.c * sh ** 2
        a = L.A * th
        f = beta * a / L.m - w
        df = ((2 * L.c / L.A) * np.sinh(2 * al / L.A) * a + beta * (1 - th ** 2)) / L.m
        step = f / df
        al = al - step
        if np.all(np.abs(step) <= 1e-14 * (1 + np.abs(al))):
            break
    return al


def curve(L, eps, al, metric):
    """On the curves from the masses eps at α: w, the coefficient's mobility μ (ẇ = μẽφ) and dμ/dw, the centre's v (ċ = vẽ∂_cφ) and
    dv/dw."""
    A, c, m, G = L.A, L.c, L.m, L.G
    th, sh = np.tanh(al / A), np.sinh(al / A)
    beta = eps + 2 * c * sh ** 2
    b1, b2 = (2 * c / A) * np.sinh(2 * al / A), (4 * c / A ** 2) * np.cosh(2 * al / A)
    a, a1, a2 = A * th, 1 - th ** 2, -2 * th * (1 - th ** 2) / A
    w = beta * a / m
    w1, w2 = (b1 * a + beta * a1) / m, (b2 * a + 2 * b1 * a1 + beta * a2) / m
    mu = G * beta * a1 * w1
    mu_a = G * (b1 * a1 * w1 + beta * a2 * w1 + beta * a1 * w2)
    v, v_a = (G * beta * a, G * (b1 * a + beta * a1)) if metric == 'paradigm' else (G * a, G * a1)
    return w, mu, mu_a / w1, v, v_a / w1


def envelope_rhs(L, eps, metric, noise, ke, cstar, astar, terms=TERMS, quasi=False):
    """The loop averaged over the fast tones. With e = x − x_m the true error and n = Σ_k Re(N_k e^{iω_kt}):
        ė = f(x_m + e) − f̂(x_m + e + n) − k_e(e + n),   ẇ_i = μ_iφ_i(x_m)(e + n),   ċ_i = v_i∂_cφ_i(x_m)(e + n);
    the state is X̄ + Σ_k Re(Z_k e^{iω_kt}) + O(s²) with the envelopes Z_k = (E_k, Θw_k, Θc_k) of the first order — the loop linearised
    at X̄ in the frame rotating with ω_k, exact along the orbit —
        Ė = −iωE + (f′ − f̂′)E − f̂′N − Σ_i(φ_iΘw_i + w_i∂_cφ_iΘc_i) − k_eD,   D = E + N the measured error's envelope,
        Θ̇w_i = −iωΘw_i + μ_iφ_iD + (μ′_iφ_iΘw_i + μ_i∂_cφ_iΘc_i)ē,   Θ̇c_i = −iωΘc_i + v_i∂_cφ_iD + (v′_i∂_cφ_iΘw_i + v_i∂²_cφ_iΘc_i)ē,
    and X̄ = (ē, w̄, c̄) the loop plus the tone averages ⟨Re(Ae^{iωt})Re(Be^{iωt})⟩ = ½Re(AB̄) of its quadratic terms: in ė the
    plant's and the control's curvature ½f″⟨e₁²⟩ − ½f̂″⟨(e₁ + n)²⟩ ('ctl'), the control's slope moved by the jitter −⟨∇_θf̂′·θ₁(e₁ + n)⟩
    ('cross') and the blur −½⟨θ₁ᵀ∇²_θf̂θ₁⟩ ('blur'); in ẇ, ċ the regressors' and mobilities' correlation with the drive
    ⟨∂_θ(Γg)·θ₁(e₁ + n)⟩ ('corr'); f̂ and f at x_m + ē, the regressor at x_m; the tones' cross averages dropped (distinct frequencies).
    `terms` selects the averages; quasi=True replaces the envelopes by their quasi-static values (Ż = 0 at ē = 0:
    D = iω(iω − f′)N/(K − ω² + ik_eω), Θ = Γg·D/(iω): A's oscillator law). Returns (rhs, the state's size)."""
    om = np.asarray(noise['w'], float)
    Nc = np.asarray(noise['amp'], float) * np.exp(1j * np.asarray(noise['ph'], float))
    T, nk = om.size, len(cstar)
    cstar, astar = np.asarray(cstar, float), np.asarray(astar, float)
    st = dict(al=None)

    def kern(x, c):
        u = (x - c) / CT.SIG
        ph = np.exp(-u ** 2 / 2)
        return ph, ph * u / CT.SIG, ph * (u ** 2 - 1) / CT.SIG ** 2

    def h(A_, B_):
        return 0.5 * np.real(A_ * np.conj(B_))

    def rhs(t, y):
        e, w, c = y[0], y[1:1 + nk], y[1 + nk:1 + 2 * nk]
        al = alpha_of(L, eps, w, st['al'] if st['al'] is not None else np.sign(w) * 3.0)
        st['al'] = al
        _, mu, dmu, v, dv = curve(L, eps, al, metric)
        xm = np.sin(CT.OM * t)
        xs = xm + e
        g_, g1, g2 = kern(xs, cstar)
        fs, f1, f2 = (astar * g_).sum(), -(astar * g1).sum(), (astar * g2).sum()
        phs, d1s, d2s = kern(xs, c)
        phr, d1r, d2r = kern(xm, c)
        fh, fh1, fh2 = (w * phs).sum(), -(w * d1s).sum(), (w * d2s).sum()
        if quasi:
            Keff = (phs * mu * phr + w * d1s * v * d1r).sum()
            D = 1j * om * (1j * om - f1) * Nc / (Keff - om ** 2 + 1j * ke * om)
            E1 = D - Nc
            Tw, Tc = (mu * phr)[None, :] * (D / (1j * om))[:, None], (v * d1r)[None, :] * (D / (1j * om))[:, None]
            dZ = None
        else:
            z = y[1 + 2 * nk:].reshape(2, T, 1 + 2 * nk)
            Z = z[0] + 1j * z[1]
            E1, Tw, Tc = Z[:, 0], Z[:, 1:1 + nk], Z[:, 1 + nk:]
            D = E1 + Nc
            dE = -1j * om * E1 + (f1 - fh1) * E1 - fh1 * Nc - (phs * Tw + w * d1s * Tc).sum(1) - ke * D
            dTw = -1j * om[:, None] * Tw + (mu * phr)[None, :] * D[:, None] + (dmu * phr * Tw + mu * d1r * Tc) * e
            dTc = -1j * om[:, None] * Tc + (v * d1r)[None, :] * D[:, None] + (dv * d1r * Tw + v * d2r * Tc) * e
            dZ = np.concatenate([dE[:, None], dTw, dTc], 1)
        Qe = 0.0
        if 'ctl' in terms:
            Qe += (0.5 * f2 * h(E1, E1) - 0.5 * fh2 * h(D, D)).sum()
        if 'cross' in terms:
            Qe += h(d1s * Tw + w * d2s * Tc, D[:, None]).sum()
        if 'blur' in terms:
            Qe -= 0.5 * (2 * d1s * h(Tw, Tc) + w * d2s * h(Tc, Tc)).sum()
        corr = 'corr' in terms
        Qw = h(dmu * phr * Tw + mu * d1r * Tc, D[:, None]).sum(0) if corr else 0.0
        Qc = h(dv * d1r * Tw + v * d2r * Tc, D[:, None]).sum(0) if corr else 0.0
        head = np.r_[fs - fh - ke * e + Qe, mu * phr * e + Qw, v * d1r * e + Qc]
        return head if dZ is None else np.r_[head, dZ.real.ravel(), dZ.imag.ravel()]
    return rhs, 1 + 2 * nk + (0 if quasi else 2 * T * (1 + 2 * nk))


def envelope_means(L, eps, metric, noise, ke, NP, k=400, terms=TERMS, quasi=False):
    """The averaged loop from the fitted state (X̄ = (0, A*, c*), Z = 0: the true error zero as in the loop) for NP periods of the
    reference: the formed kernels' period means of w̄ and c̄ (trapezoid on k points a period), the solver's status."""
    from scipy.integrate import solve_ivp
    rhs, n = envelope_rhs(L, eps, metric, noise, ke, CT.CSTAR, CT.ASTAR, terms, quasi)
    y0 = np.zeros(n)
    y0[1:4], y0[4:7] = CT.ASTAR, CT.CSTAR
    per = 2 * np.pi / CT.OM
    ts = per * np.arange(NP * k + 1) / k
    sol = solve_ivp(rhs, (0.0, ts[-1]), y0, method='LSODA', t_eval=ts, rtol=1e-8, atol=1e-11)
    Y = sol.y.T
    n_p = min(NP, (len(Y) - 1) // k)
    pm = lambda cols: np.array([np.trapezoid(Y[j * k:(j + 1) * k + 1, cols], axis=0) / k for j in range(n_p)])
    return pm(slice(1, 4)), pm(slice(4, 7)), int(sol.status)


def sim_means(L, metric, noise, ke, NP, chunk=10, k=2000):
    """The loop from the fitted state for NP periods: the formed kernels' period means of w and c (trapezoid on k dense outputs a
    period — the jitter averaged, not aliased), in chunks of `chunk` periods (bounded memory)."""
    al, c, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    per = 2 * np.pi / CT.OM
    x, Wp, Cp = None, [], []
    for p0 in range(0, NP, chunk):
        n = min(chunk, NP - p0)
        ts = per * (p0 + np.arange(n * k + 1) / k)
        r = simulate(L, al, c, CT.ASTAR, ts, noise, metric, ke, x_init=x)
        for j in range(n):
            sl = slice(j * k, (j + 1) * k + 1)
            Wp.append(np.trapezoid(r['w'][sl][:, idx], axis=0) / k)
            Cp.append(np.trapezoid(r['c'][sl][:, idx], axis=0) / k)
        x, al, c = r['x'][-1], r['al'][-1].copy(), r['c'][-1].copy()
    return np.array(Wp), np.array(Cp)


def second_order(L, eps, al, c, metric, noise, ke, NP, k=2000):
    """The loop expanded about the exact fit along x_m(t) to second order in the noise (the formed kernels; no averaging): with
    n = Σ_k amp_k cos(ω_kt + ph_k), d₁ = e₁ + n the measured error's first order, f′ = f̂′ and f″ = f̂″ at the fit,
        e₁′ = −f′n − Σ_i(φ_iw₁_i + w_i∂_cφ_ic₁_i) − k_ed₁,   w₁_i′ = μ_iφ_id₁,   c₁_i′ = v_i∂_cφ_id₁,
        e₂′ = −Σ_i(φ_iw₂_i + w_i∂_cφ_ic₂_i) − k_ee₂ − ½f″n(2e₁ + n) + Σ_i(∂_cφ_iw₁_i + w_i∂²_cφ_ic₁_i)d₁
              − Σ_i(∂_cφ_iw₁_ic₁_i + ½w_i∂²_cφ_ic₁_i²),
        w₂_i′ = μ_iφ_ie₂ + (μ′_iφ_iw₁_i + μ_i∂_cφ_ic₁_i)d₁,   c₂_i′ = v_i∂_cφ_ie₂ + (v′_i∂_cφ_iw₁_i + v_i∂²_cφ_ic₁_i)d₁,
    from zero (the loop starts at x = x_m) — the terms of `envelope_rhs`'s tone averages ('ctl', 'cross', 'blur'; 'corr' in w₂′, c₂′)
    before averaging. The loop's period means are then those of s·(w₁, c₁) + s²·(w₂, c₂) + O(s³) for the noise s·n. Returns the
    period means (trapezoid on k points a period, as `sim_means`) of w₁, c₁, w₂, c₂ over NP periods of the reference."""
    from scipy.integrate import solve_ivp
    w, mu, dmu, v, dv = curve(L, eps, al, metric)
    om, ph0, amp = (np.asarray(noise[q], float) for q in ('w', 'ph', 'amp'))
    nk = c.size

    def rhs(t, y):
        xm = np.sin(CT.OM * t)
        u = (xm - c) / CT.SIG
        ph = np.exp(-u ** 2 / 2)
        d1, d2 = ph * u / CT.SIG, ph * (u ** 2 - 1) / CT.SIG ** 2
        f1, f2 = -(w * d1).sum(), (w * d2).sum()
        n = (amp * np.cos(om * t + ph0)).sum()
        e1, w1, c1 = y[0], y[1:1 + nk], y[1 + nk:1 + 2 * nk]
        e2, w2, c2 = y[1 + 2 * nk], y[2 + 2 * nk:2 + 3 * nk], y[2 + 3 * nk:]
        dm = e1 + n
        q = -0.5 * f2 * n * (2 * e1 + n) + ((d1 * w1 + w * d2 * c1) * dm).sum() - (d1 * w1 * c1 + 0.5 * w * d2 * c1 ** 2).sum()
        return np.r_[-f1 * n - (ph * w1 + w * d1 * c1).sum() - ke * dm, mu * ph * dm, v * d1 * dm,
                     -(ph * w2 + w * d1 * c2).sum() - ke * e2 + q, mu * ph * e2 + (dmu * ph * w1 + mu * d1 * c1) * dm,
                     v * d1 * e2 + (dv * d1 * w1 + v * d2 * c1) * dm]
    per = 2 * np.pi / CT.OM
    ts = per * np.arange(NP * k + 1) / k
    sol = solve_ivp(rhs, (0.0, ts[-1]), np.zeros(2 + 4 * nk), method='DOP853', t_eval=ts, rtol=1e-10, atol=1e-14,
                    max_step=0.2 / om.max())
    Y = sol.y.T
    pm = lambda cols: np.array([np.trapezoid(Y[j * k:(j + 1) * k + 1, cols], axis=0) / k for j in range(NP)])
    return (pm(slice(1, 1 + nk)), pm(slice(1 + nk, 1 + 2 * nk)), pm(slice(2 + 2 * nk, 2 + 3 * nk)), pm(slice(2 + 3 * nk, 2 + 4 * nk)),
            int(sol.status))


def b3_job(args):
    """B3: the exact second order against the loop — the loop's period means under the noise s·n and −s·n (the tones' phases moved by
    π) at the marks MARKS_B: their half-difference (the odd orders: s·⟨w₁⟩ + O(s³)) and half-sum less the fit (the even orders:
    s²⟨w₂⟩ + O(s⁴)), at s ∈ S_B3, against `second_order`'s period means for the unit noise (tones `freqs`, or D‴.49's 16 if None)."""
    metric, freqs = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law()
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    opts = {} if freqs is None else dict(freqs=freqs)
    W1, C1, W2, C2, stt = second_order(L, EPS, al0[idx], c0[idx], metric, DZ.noise_of(1.0, opts), CT.KE, NP_B)
    t1 = time.time()
    mk = [m - 1 for m in MARKS_B]
    rows = []
    for s in S_B3:
        nz = DZ.noise_of(s, opts)
        Wp, Cp = sim_means(L, metric, nz, CT.KE, NP_B)
        Wm, Cm = sim_means(L, metric, dict(nz, ph=np.asarray(nz['ph']) + np.pi), CT.KE, NP_B)
        rows.append(dict(s=s, even_w=((Wp + Wm)[mk] / 2 - CT.ASTAR).tolist(), odd_w=((Wp - Wm)[mk] / 2).tolist(),
                         even_c=((Cp + Cm)[mk] / 2 - CT.CSTAR).tolist(), odd_c=((Cp - Cm)[mk] / 2).tolist()))
    out = dict(metric=metric, freqs=None if freqs is None else list(freqs), status=stt, w1=W1[mk].tolist(), c1=C1[mk].tolist(),
               w2=W2[mk].tolist(), c2=C2[mk].tolist(), rows=rows, secs_theory=t1 - t0, secs=time.time() - t0)
    print(f"  B3 {metric} {noise_label(out['freqs'])}: even/s²⟨w₂⟩ at {MARKS_B[-1]}: "
          + '; '.join(f"s {r['s']:g}: " + ', '.join(f"{a / (r['s'] ** 2 * b):.4f}" for a, b in zip(r['even_w'][-1], out['w2'][-1]))
                      for r in rows) + f" [{out['secs_theory']:.0f}s + {out['secs'] - out['secs_theory']:.0f}s]", flush=True)
    return out


def b_job(args):
    """B: the loop's period means against the averaged loop's (`envelope_means`; `terms`, `quasi` its variants) from the fitted state, at
    the marks MARKS_B, under tones `freqs` of total rms s (equal amplitudes; D‴.49's 16 tones if None)."""
    metric, freqs, s, terms, quasi = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law()
    nz = DZ.noise_of(s, {} if freqs is None else dict(freqs=freqs))
    We, Ce, stt = envelope_means(L, EPS, metric, nz, CT.KE, NP_B, terms=terms, quasi=quasi)
    t1 = time.time()
    Ws, Cs = sim_means(L, metric, nz, CT.KE, NP_B)
    mk = [m - 1 for m in MARKS_B]
    out = dict(metric=metric, freqs=None if freqs is None else list(freqs), s=s, terms=list(terms), quasi=quasi, status=stt,
               dw_sim=(Ws[mk] - CT.ASTAR).tolist(), dc_sim=(Cs[mk] - CT.CSTAR).tolist(),
               dw_env=(We[mk] - CT.ASTAR).tolist(), dc_env=(Ce[mk] - CT.CSTAR).tolist(), secs_env=t1 - t0, secs=time.time() - t0)
    print(f"  B {metric} {noise_label(out['freqs'])} s {s:g} {variant_label(out)}: Δw at {MARKS_B[-1]}: loop "
          f"{', '.join(f'{x:+.2e}' for x in out['dw_sim'][-1])}, averaged {', '.join(f'{x:+.2e}' for x in out['dw_env'][-1])} "
          f"[{out['secs']:.0f}s]", flush=True)
    return out


def noise_label(freqs):
    return "D‴.49's 16 tones" if freqs is None else 'one tone ' + ', '.join(f'{f:g}' for f in freqs)


def variant_label(r):
    if r['quasi']:
        return 'quasi-static envelopes'
    miss = [t for t in TERMS if t not in r['terms']]
    return 'all averages' if not miss else 'without ' + ', '.join(TERM_NAMES[t] for t in miss)


def c_noise(s, opts):
    """C's noise of rms s: D‴.49's 16 tones in 10–30 rad/s (opts seed), in another band (opts band), or the tones opts freqs (equal
    amplitudes, the seed's phases)."""
    from src.core import closed_loop as CL
    if 'band' in opts:
        return CL.tones(s, *opts['band'], seed=opts.get('seed', 7))
    return DZ.noise_of(s, {k: v for k, v in opts.items() if k in ('freqs', 'seed')})


def c_job(args):
    """C: the metric (opts metric, default the paradigm's) at the law's m and the error gain k_e, the noise of rms s (`c_noise`) from the
    fitted state for NP_C periods: per formed kernel the first period whose period-averaged centre is LEFT_C from its feature (−1: never),
    the largest period-averaged distance, the coefficients at the end."""
    m, ke, s = args[:3]
    opts = dict(args[3]) if len(args) > 3 else {}
    metric = opts.get('metric', 'paradigm')
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law(m)
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    k = 8
    ts = np.linspace(0, NP_C * 2 * np.pi / CT.OM, NP_C * k + 1)
    r = simulate(L, al0, c0, CT.ASTAR, ts, c_noise(s, opts), metric, ke)
    dev = np.abs(r['c'][:-1].reshape(NP_C, k, c0.size).mean(1)[:, idx] - CT.CSTAR)
    o = orbit_terms(L, EPS, al0[idx], c0[idx], metric, np.sin(np.linspace(0, 2 * np.pi, 2001)))
    out = dict(m=m, ke=ke, s=s, opts={k_: (list(v) if isinstance(v, tuple) else v) for k_, v in opts.items()}, left=left_periods(dev, LEFT_C),
               maxdev=dev.max(0).tolist(), w_end=[float(r['w'][-1, i]) for i in idx], sqrtK=float(np.sqrt(o['K'].max())),
               secs=time.time() - t0)
    print(f"  C m {m:g} k_e {ke:g} s {s:g} {out['opts']}: left {out['left']} [{out['secs']:.0f}s]", flush=True)
    return out


def b4_job(args):
    """B4: the loop under the noise s·n and −s·n for NP_B4 periods (`sim_means`): the even and odd parts of the period means at
    MARKS_B4, per s² and per s."""
    freqs, s = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law()
    nz = DZ.noise_of(s, {} if freqs is None else dict(freqs=freqs))
    Wp, Cp = sim_means(L, 'paradigm', nz, CT.KE, NP_B4)
    Wm, Cm = sim_means(L, 'paradigm', dict(nz, ph=np.asarray(nz['ph']) + np.pi), CT.KE, NP_B4)
    mk = [m - 1 for m in MARKS_B4]
    out = dict(freqs=None if freqs is None else list(freqs), s=s, even_w=(((Wp + Wm) / 2 - CT.ASTAR)[mk] / s ** 2).tolist(),
               odd_w=(((Wp - Wm) / 2)[mk] / s).tolist(), even_c=(((Cp + Cm) / 2 - CT.CSTAR)[mk] / s ** 2).tolist(),
               odd_c=(((Cp - Cm) / 2)[mk] / s).tolist(), odd_c_per=((Cp - Cm)[:12, 0] / (2 * s)).tolist(), secs=time.time() - t0)
    print(f"  B4 {noise_label(out['freqs'])} s {s:g}: even/s² at {MARKS_B4[-1]}: {', '.join(f'{x:+.4g}' for x in out['even_w'][-1])} "
          f"[{out['secs']:.0f}s]", flush=True)
    return out


def left_periods(dev, cut):
    """Per kernel (column) the first period (from 1) whose period-averaged distance exceeds cut, −1 if none."""
    return [int(np.argmax(dev[:, q] > cut)) + 1 if (dev[:, q] > cut).any() else -1 for q in range(dev.shape[1])]


def thresholds(rs, nq=3):
    """Per kernel the smallest level among the runs rs at which it left (None: at none)."""
    rs = sorted(rs, key=lambda r: r['s'])
    return [next((r['s'] for r in rs if r['left'][q] >= 0), None) for q in range(nq)]


def b_jobs():
    return [(mt, (om,), s, TERMS, False) for mt in METRICS for om in OM_B for s in S_B] + [(mt, None, S_B[-1], TERMS, False)
                                                                                          for mt in METRICS]


def b3_jobs():
    return [(mt, (om,)) for mt in METRICS for om in OM_B] + [(mt, None) for mt in METRICS]


def b2_jobs():
    mt, fr, s = B2_CASE
    return [(mt, fr, s, tuple(t for t in TERMS if t != drop), False) for drop in TERMS] + [(mt, fr, s, TERMS, True)]


def _pool_map(fn, jobs):
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(fn, jobs, chunksize=1)


def part_A():
    return _pool_map(a_job, [(mt, om) for mt in METRICS for om in OM_A])


def part_B():
    return _pool_map(b_job, b_jobs())


def part_B2():
    return _pool_map(b_job, b2_jobs())


def part_B3():
    return _pool_map(b3_job, b3_jobs())


def part_C():
    return _pool_map(c_job, [(m, ke, s, dict(seed=sd)) for m in M_C for ke in KE_C for s in LEV_C for sd in SEEDS_C])


def part_C2():
    return [dict(r, label=lab) for (lab, *_), r in zip(C2_CASES, _pool_map(c_job, [cs[1:] for cs in C2_CASES]))]


def part_B4():
    return _pool_map(b4_job, [(fr, s) for fr in FREQS_B4 for s in S_B4])


def ratio_range(rows, num, den, cond=lambda r: True):
    """The least and largest ratio num/den over the rows passing cond and their kernels."""
    q = [a / b for r in rows if cond(r) for a, b in zip(r[num], r[den])]
    return (min(q), max(q)) if q else (np.nan, np.nan)


def b3_lines(B3, B=None):
    """B3's report: the loop's even and odd parts against the exact second order, and the averaged loop (B) against it."""
    if B3 is None:
        return []
    out = ["## B3. The exact second order", "",
           f"The loop expanded about the fit to second order in the noise along the orbit, without averaging (`second_order`: the first "
           "order and the second driven by its products — the terms B averages), against the loop under the noise s·n and −s·n (the tones' "
           "phases moved by π): its period means' half-sum less the fit (the even orders, s²⟨w₂⟩ + O(s⁴)) and half-difference (the odd "
           f"orders, s⟨w₁⟩ + O(s³)), at period {MARKS_B[-1]}; the averaged loop's Δw (B, at the largest s) over s²⟨w₂⟩:", "",
           "| metric | noise | ⟨w₂⟩ | even part / s²⟨w₂⟩ at s = " + ', '.join(f"{q:g}" for q in S_B3) + " | the averaged loop (B) / s²⟨w₂⟩ "
           f"| odd part at s = {S_B3[-1]:g}: loop (s⟨w₁⟩) | ⟨c₂⟩ | even part / s²⟨c₂⟩ at s = " + ', '.join(f"{q:g}" for q in S_B3) + " |",
           "|---|---|---|---|---|---|---|---|"]
    key = lambda r: (r['metric'] != 'paradigm', r['freqs'] is None, r['freqs'] or [0])
    q_even, q_avg, q_even_c = {q: [] for q in S_B3}, {}, {q: [] for q in S_B3}
    for r in sorted(B3, key=key):
        w2, c2 = r['w2'][-1], r['c2'][-1]
        cells = []
        for row in r['rows']:
            qs = [a / (row['s'] ** 2 * b) for a, b in zip(row['even_w'][-1], w2)]
            q_even[row['s']] += qs
            cells.append(', '.join(f"{x:.4f}" for x in qs))
        b = next((q for q in (B if B is not None else []) if q['metric'] == r['metric'] and q['freqs'] == r['freqs'] and q['s'] == S_B3[-1]
                  and len(q['terms']) == len(TERMS) and not q['quasi']), None)
        avg = '—'
        if b is not None:
            qa = [a / (S_B3[-1] ** 2 * c) for a, c in zip(b['dw_env'][-1], w2)]
            q_avg[(r['metric'], noise_label(r['freqs']))] = qa
            avg = ', '.join(f"{x:.3f}" for x in qa)
        top = r['rows'][-1]
        odd = ', '.join(f"{a:+.2e} ({S_B3[-1] * c:+.2e})" for a, c in zip(top['odd_w'][-1], r['w1'][-1]))
        ecs = []
        for row in r['rows']:
            qc = [a / (row['s'] ** 2 * c) for a, c in zip(row['even_c'][-1], c2)]
            q_even_c[row['s']] += qc
            ecs.append(', '.join(f"{x:.4f}" for x in qc))
        ec = '; '.join(ecs)
        out.append(f"| {r['metric']} | {noise_label(r['freqs'])} | {', '.join(f'{x:+.4g}' for x in w2)} | {'; '.join(cells)} | {avg} | {odd} | "
                   f"{', '.join(f'{x:+.4g}' for x in c2)} | {ec} |")
    rng = lambda xs: f"[{min(xs):.4f}, {max(xs):.4f}]" if xs else '—'
    out += ["", "The even part over the exact second order: " + '; '.join(f"{rng(q_even[q])} at s = {q:g}" for q in S_B3)
            + " (the relative error falls as s²); the centres' even part over s²⟨c₂⟩: " + '; '.join(f"{rng(q_even_c[q])} at s = {q:g}"
                                                                                                 for q in S_B3)
            + "; the averaged loop over it at s = " + f"{S_B3[-1]:g}: "
            + '; '.join(f"{k[0]}, {k[1]}: {rng(v)}" for k, v in q_avg.items()) + ".", ""]
    return out


def c_cells(rs, nq=3):
    """C's cells over the seeds: per kernel the least level at which a seed left (the seeds left there), the least level from which every
    seed left at every level above ('—' if none), '*' where a seed that left at one level stays at a higher one."""
    lv = sorted({r['s'] for r in rs})
    out = []
    for q in range(nq):
        frac = [(s_, sum(r['left'][q] >= 0 for r in rs if r['s'] == s_), sum(r['s'] == s_ for r in rs)) for s_ in lv]
        first = next(((s_, k, n) for s_, k, n in frac if k), None)
        allf = next((s_ for i, (s_, k, n) in enumerate(frac) if all(k2 == n2 for _, k2, n2 in frac[i:])), None)
        seeds = sorted({r.get('opts', {}).get('seed', 7) for r in rs})
        mono = all(not (r1['left'][q] >= 0 and r2['left'][q] < 0) for sd in seeds for r1 in rs for r2 in rs
                   if r1.get('opts', {}).get('seed', 7) == sd == r2.get('opts', {}).get('seed', 7) and r1['s'] < r2['s'])
        out.append(('—' if first is None else f"{first[0]:g} ({first[1]}/{first[2]})") + ' / ' + ('—' if allf is None else f"{allf:g}")
                   + ('' if mono else ' *'))
    return out


def report(A=None, B=None, B2=None, C=None, B3=None, B4=None, C2=None):
    L_ = ["# Proposition D‴.60: measurement noise and the formed kernels (`exp_noise_resonance.py`)", "",
          f"D‴.49's closed loop (k_e = {CT.KE:g}, ω = {CT.OM:g}, three bumps of width {CT.SIG:g}, 21 kernels at the masses ε = {EPS:g}, the "
          "regressor at x_m) from the fitted state — each feature's kernel at its centre with its amplitude, the others at α = 0 —; the "
          "paradigm's metric (the centres' mobility Gβa = Gmw) and the Wasserstein–Fisher–Rao metric (G·a).", ""]
    if A is not None:
        L_ += ["## A. The linear response", "",
               f"One tone of amplitude {AP_A:g} at ω for {NP_A} periods of the reference: per formed kernel the rms jitter (the deviation "
               "from the moving average over one tone period) of the loop and of the quasi-static oscillator law — y″ + k_eẏ + K(t)y = "
               "K(ṅ − f̂′n), δc_i = (v_i∂_cφ_i/K)y, δw_i = (μ_iφ_i/K)y — over the loop linearised at the fit along the orbit, through the "
               "same filter; the linearised centre jitter per unit noise amplitude:", "",
               "| metric | √K max | ω | centres: loop / linearised | centres: quasi-static / linearised | coefficients: loop / linearised "
               "| coefficients: quasi-static / linearised | linearised centre jitter / amplitude |", "|---|---|---|---|---|---|---|---|"]
        f3 = lambda r, a, b, fm: ', '.join(f"{x / y:{fm}}" for x, y in zip(r[a], r[b]))
        for r in sorted(A, key=lambda r: (r['metric'] != 'paradigm', r['om'])):
            L_.append(f"| {r['metric']} | {r['sqrtK']:.2f} | {r['om']:g} | {f3(r, 'sim_c', 'var_c', '.4f')} | {f3(r, 'qs_c', 'var_c', '.3f')} | "
                      f"{f3(r, 'sim_w', 'var_w', '.4f')} | {f3(r, 'qs_w', 'var_w', '.3f')} | {', '.join(f'{x / AP_A:.3f}' for x in r['var_c'])} |")
        lo, hi = ratio_range(A, 'sim_c', 'var_c')
        lw, hw = ratio_range(A, 'sim_w', 'var_w')
        q10 = ratio_range(A, 'qs_c', 'var_c', lambda r: r['om'] >= 10)
        q19 = ratio_range(A, 'qs_c', 'var_c', lambda r: r['om'] >= 19)
        mob = '; '.join(f"{r['metric']}: v = {', '.join(f'{x:.3g}' for x in r['v'])}, μ = {', '.join(f'{x:.3g}' for x in r['mu'])}"
                        for r in [next(q for q in A if q['metric'] == mt) for mt in METRICS if any(q['metric'] == mt for q in A)])
        per_m = '; '.join(f"{mt}: centres [{ratio_range(A, 'sim_c', 'var_c', lambda r, m_=mt: r['metric'] == m_)[0]:.4f}, "
                          f"{ratio_range(A, 'sim_c', 'var_c', lambda r, m_=mt: r['metric'] == m_)[1]:.4f}], coefficients "
                          f"[{ratio_range(A, 'sim_w', 'var_w', lambda r, m_=mt: r['metric'] == m_)[0]:.4f}, "
                          f"{ratio_range(A, 'sim_w', 'var_w', lambda r, m_=mt: r['metric'] == m_)[1]:.4f}]"
                          for mt in METRICS if any(q['metric'] == mt for q in A))
        L_ += ["", f"The loop over its linearisation: centres [{lo:.4f}, {hi:.4f}], coefficients [{lw:.4f}, {hw:.4f}] (per metric: {per_m}); "
                   f"the quasi-static law over the linearisation, centres: [{q10[0]:.3f}, {q10[1]:.3f}] at ω ≥ 10, [{q19[0]:.3f}, "
                   f"{q19[1]:.3f}] at ω ≥ 19. The mobilities at the fit: {mob}.", ""]
    if B is not None:
        evn = {}
        for r3 in (B3 if B3 is not None else []):
            for row in r3['rows']:
                evn[(r3['metric'], tuple(r3['freqs'] or ()), row['s'])] = row['even_w'][-1]
        L_ += ["## B. The second order: the averaged loop", "",
               f"Tones of total rms s (one tone, or D‴.49's 16 in 10–30 rad/s), {NP_B} periods from the fitted state: the formed kernels' "
               "period means of w and c (the loop's from dense outputs) against the loop averaged over the fast tones (the mean state and "
               "the complex envelopes of the first order, exact along the orbit; the quadratic terms' tone averages): Δw = w̄ − A, Δc = "
               f"c̄ − c* at the periods {', '.join(str(m) for m in MARKS_B)}:", "",
               "| metric | noise | s | Δw: loop (averaged), at " + ', '.join(str(m) for m in MARKS_B) + " | averaged / the loop's even "
               f"part (B3; the loop's own Δw includes its odd part), at {MARKS_B[-1]} | Δc at {MARKS_B[-1]}: loop (averaged) |",
               "|---|---|---|---|---|---|"]
        key = lambda r: (r['metric'] != 'paradigm', r['freqs'] is None, r['freqs'] or [0], r['s'])
        for r in sorted(B, key=key):
            cells = '; '.join(', '.join(f"{a:+.2e} ({b:+.2e})" for a, b in zip(r['dw_sim'][j], r['dw_env'][j])) for j in range(len(MARKS_B)))
            ev = evn.get((r['metric'], tuple(r['freqs'] or ()), r['s']))
            L_.append(f"| {r['metric']} | {noise_label(r['freqs'])} | {r['s']:g} | {cells} | "
                      + ('—' if ev is None else ', '.join(f'{b / a:.3f}' for a, b in zip(ev, r['dw_env'][-1]))) + ' | '
                      + ', '.join(f"{a:+.2e} ({b:+.2e})" for a, b in zip(r['dc_sim'][-1], r['dc_env'][-1])) + " |")
        L_.append("")
    if B2 is not None:
        mt, fr, s = B2_CASE
        full = next((r for r in (B if B is not None else []) if r['metric'] == mt and r['freqs'] == list(fr) and r['s'] == s), None)
        L_ += ["## B2. Which averages drive the drift", "",
               f"{mt}, {noise_label(list(fr))}, s = {s:g}: the averaged loop's Δw and Δc at period {MARKS_B[-1]} without each tone average, "
               "and with the envelopes replaced by their quasi-static values (A's oscillator law), against the full average"
               + ("" if full is None else f" (Δw {', '.join(f'{x:+.2e}' for x in full['dw_env'][-1])}; the loop's "
                  f"{', '.join(f'{x:+.2e}' for x in full['dw_sim'][-1])})") + ":", "",
               "| variant | Δw | Δw / the full average's | Δc |", "|---|---|---|---|"]
        for r in B2:
            rel = ('—' if full is None else ', '.join(f"{a / b:.3f}" for a, b in zip(r['dw_env'][-1], full['dw_env'][-1])))
            L_.append(f"| {variant_label(r)} | {', '.join(f'{x:+.2e}' for x in r['dw_env'][-1])} | {rel} | "
                      f"{', '.join(f'{x:+.2e}' for x in r['dc_env'][-1])} |")
        L_.append("")
    L_ += b3_lines(B3, B)
    if B4 is not None:
        L_ += ["## B4. The long horizon and commensurability", "",
               f"The paradigm's metric, the loop under s·n and −s·n for {NP_B4} periods: the even part of the coefficients' period means per s² "
               f"and the odd part per s at the periods {', '.join(str(m) for m in MARKS_B4)} (kernel by kernel), at s = "
               f"{', '.join(f'{q:g}' for q in S_B4)}; the reference's ω = {CT.OM:g} (21.4 = 107ω, 21.3 = 106.5ω, 13.1 = 65.5ω):", "",
               "| noise | s | even part / s² (w), at the periods | odd part / s (w), at the periods | the centres' even part / s² at the "
               "last period | the first kernel's centre, odd part / s, periods 1–12 |", "|---|---|---|---|---|---|"]
        for r in sorted(B4, key=lambda r: (r['freqs'] is None, r['freqs'] or [0], r['s'])):
            ev = '; '.join(', '.join(f"{x:+.3g}" for x in row) for row in r['even_w'])
            od = '; '.join(', '.join(f"{x:+.3g}" for x in row) for row in r['odd_w'])
            L_.append(f"| {noise_label(r['freqs'])} | {r['s']:g} | {ev} | {od} | {', '.join(f'{x:+.3g}' for x in r['even_c'][-1])} | "
                      f"{', '.join(f'{x:+.3f}' for x in r['odd_c_per'])} |")
        L_.append("")
    if C is not None:
        xm = np.sin(np.linspace(0, 2 * np.pi, 2001))
        kk = []
        for m_ in M_C:                                                 # each law's own fitted state (as `c_job`'s)
            Lk = law(m_)
            al0, c0, idx = DZ.fitted_state(Lk, EPS, CT.CSTAR, CT.ASTAR)
            o = orbit_terms(Lk, EPS, al0[idx], c0[idx], 'paradigm', xm)
            kk.append(f"m = {m_:g}: " + ', '.join(f"{v:.2f}" for v in np.sqrt((o['mu'] * o['ph'] ** 2 + o['v'] * o['w'] * o['d1'] ** 2).max(0))))
        L_ += ["## C. Leaving the features: the error gain, the centres' mobility, the seeds", "",
               f"The paradigm's metric, D‴.49's noise of rms s (16 tones in 10–30 rad/s; the seeds {', '.join(map(str, SEEDS_C))} draw their "
               f"frequencies and phases, 7 is D‴.49's), {NP_C} periods from the fitted state: per formed kernel the least level (of "
               f"{', '.join(f'{s:g}' for s in LEV_C)}) at which a seed left its feature ({LEFT_C:g} away, period-averaged) with the seeds "
               "that left there, and the least level from which every seed left at every level above ('—': none; '*': a seed left at a "
               "level and stayed at a higher one), against k_e and the law's m (the centres' mobility Gmw); √K's largest value along the "
               "orbit, and per kernel (μ_iφ_i² + v_iw_i(∂_cφ_i)²)^{1/2} at its largest: " + '; '.join(kk) + ":", "",
               "| m | k_e | √K max | " + ' | '.join(f"A = {a:g}" for a in CT.ASTAR) + " | the median period of leaving (all levels) |",
               "|---|---|---|---|---|---|---|"]
        for m in M_C:
            for ke in KE_C:
                rs = [r for r in C if r['m'] == m and r['ke'] == ke]
                if not rs:
                    continue
                per_ = [r['left'][q] for r in rs for q in range(3) if r['left'][q] >= 0]
                L_.append(f"| {m:g} | {ke:g} | {rs[0]['sqrtK']:.2f} | " + ' | '.join(c_cells(rs)) + f" | "
                          + (f"{np.median(per_):.0f} ({len(per_)} leaves)" if per_ else '—') + " |")
        L_.append("")
    if C2 is not None:
        L12 = law(12.0)
        al_, c_, ix_ = DZ.fitted_state(L12, EPS, CT.CSTAR, CT.ASTAR)
        o_ = orbit_terms(L12, EPS, al_[ix_], c_[ix_], 'paradigm', np.sin(np.linspace(0, 2 * np.pi, 2001)))
        S_ = float((o_['w'] ** 2 * o_['d1'] ** 2).sum(1).max())
        L_ += ["## C2. Controls (k_e and m as stated)", "",
               f"The same loop and horizon: one tone; the 16 tones moved to 25–75 rad/s (above √K_max); m between 3 and 12 (the design "
               f"bound m ≲ ω_lo²/(G max_t Σ_iw_i²(∂_cφ_i)²) = {100.0 / (L12.G * S_):.2f} at ω_lo = 10, the band's lower edge: √K_max at "
               "ω_lo); the Wasserstein–Fisher–Rao metric; tone pairs 18.9/22.5 (the difference 3.6 = 18ω) and 18.9/22.4 (17.5ω) under four "
               "seeds' phases. Per run the period each kernel left (−1: never) and its largest period-averaged distance from the feature:", "",
               "| control | m | k_e | noise | s | left at the period | largest distance |", "|---|---|---|---|---|---|---|"]
        for r in C2:
            o = r['opts']
            nz = ('16 tones in ' + '–'.join(f"{x:g}" for x in o['band']) if 'band' in o else
                  (noise_label(o['freqs']) if 'freqs' in o else "D‴.49's 16 tones")) + (f", seed {o['seed']}" if 'seed' in o else '')
            L_.append(f"| {r['label']} ({o.get('metric', 'paradigm')}) | {r['m']:g} | {r['ke']:g} | {nz} | {r['s']:g} | "
                      f"{', '.join(map(str, r['left']))} | {', '.join(f'{x:.3f}' for x in r['maxdev'])} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_noise_resonance.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'noise_resonance_raw.npz')
    md = os.path.join(OUT, 'noise_resonance.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in ('A', 'B', 'B2', 'C', 'B3', 'B4', 'C2')]))
        return
    parts = (set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else
             {'A', 'B', 'B2', 'B3', 'B4', 'C', 'C2'})
    res = {k: v for k, v in old.items() if k not in parts}
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('B2', part_B2), ('B3', part_B3), ('B4', part_B4), ('C', part_C), ('C2', part_C2)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k_) for k_ in ('A', 'B', 'B2', 'C', 'B3', 'B4', 'C2')]))


if __name__ == '__main__':
    main()
