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
C. The thresholds: the smallest level of D‴.49's noise at which a formed kernel leaves its feature within NP_C periods (D‴.57's Part D),
   the paradigm's metric, against the error gain k_e (KE_C) and the law's m (M_C: the centres' mobility Gmw) — with √K_max's place
   against the noise band 10–30 rad/s.

Usage: python experiments/exp_noise_resonance.py [--only A,B,B2,C] [--from-raw] → results/chaos/noise_resonance.md (+ _raw.npz).
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
LEV_C = (0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.025, 0.03, 0.04, 0.06, 0.08)
LEFT_C = 0.3
TERMS = ('ctl', 'cross', 'blur', 'corr')
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


def c_job(args):
    """C: the paradigm's metric at the law's m and the error gain k_e, D‴.49's noise of rms s from the fitted state for NP_C periods:
    per formed kernel the first period whose period-averaged centre is LEFT_C from its feature (−1: never), the coefficients at the end."""
    m, ke, s = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law(m)
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    k = 8
    ts = np.linspace(0, NP_C * 2 * np.pi / CT.OM, NP_C * k + 1)
    r = simulate(L, al0, c0, CT.ASTAR, ts, DZ.noise_of(s, {}), 'paradigm', ke)
    dev = np.abs(r['c'][:-1].reshape(NP_C, k, c0.size).mean(1)[:, idx] - CT.CSTAR)
    o = orbit_terms(L, EPS, al0[idx], c0[idx], 'paradigm', np.sin(np.linspace(0, 2 * np.pi, 2001)))
    out = dict(m=m, ke=ke, s=s, left=left_periods(dev, LEFT_C), w_end=[float(r['w'][-1, i]) for i in idx],
               sqrtK=float(np.sqrt(o['K'].max())), secs=time.time() - t0)
    print(f"  C m {m:g} k_e {ke:g} s {s:g}: left {out['left']} [{out['secs']:.0f}s]", flush=True)
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


def part_C():
    return _pool_map(c_job, [(m, ke, s) for m in M_C for ke in KE_C for s in LEV_C])


def ratio_range(rows, num, den, cond=lambda r: True):
    """The least and largest ratio num/den over the rows passing cond and their kernels."""
    q = [a / b for r in rows if cond(r) for a, b in zip(r[num], r[den])]
    return (min(q), max(q)) if q else (np.nan, np.nan)


def report(A=None, B=None, B2=None, C=None):
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
        L_ += ["", f"The loop over its linearisation: centres [{lo:.4f}, {hi:.4f}], coefficients [{lw:.4f}, {hw:.4f}]; the quasi-static law "
                   f"over the linearisation, centres: [{q10[0]:.3f}, {q10[1]:.3f}] at ω ≥ 10, [{q19[0]:.3f}, {q19[1]:.3f}] at ω ≥ 19. The "
                   f"mobilities at the fit: {mob}.", ""]
    if B is not None:
        L_ += ["## B. The second order: the averaged loop", "",
               f"Tones of total rms s (one tone, or D‴.49's 16 in 10–30 rad/s), {NP_B} periods from the fitted state: the formed kernels' "
               "period means of w and c (the loop's from dense outputs) against the loop averaged over the fast tones (the mean state and "
               "the complex envelopes of the first order, exact along the orbit; the quadratic terms' tone averages): Δw = w̄ − A, Δc = "
               f"c̄ − c* at the periods {', '.join(str(m) for m in MARKS_B)}:", "",
               "| metric | noise | s | Δw: loop (averaged), at " + ', '.join(str(m) for m in MARKS_B) + " | averaged / loop: Δw at "
               f"{MARKS_B[-1]} | Δc at {MARKS_B[-1]}: loop (averaged) |", "|---|---|---|---|---|---|"]
        key = lambda r: (r['metric'] != 'paradigm', r['freqs'] is None, r['freqs'] or [0], r['s'])
        for r in sorted(B, key=key):
            cells = '; '.join(', '.join(f"{a:+.2e} ({b:+.2e})" for a, b in zip(r['dw_sim'][j], r['dw_env'][j])) for j in range(len(MARKS_B)))
            L_.append(f"| {r['metric']} | {noise_label(r['freqs'])} | {r['s']:g} | {cells} | "
                      f"{', '.join(f'{b / a:.3f}' for a, b in zip(r['dw_sim'][-1], r['dw_env'][-1]))} | "
                      + ', '.join(f"{a:+.2e} ({b:+.2e})" for a, b in zip(r['dc_sim'][-1], r['dc_env'][-1])) + " |")
        L_.append("")
    if B2 is not None:
        mt, fr, s = B2_CASE
        full = next((r for r in (B or []) if r['metric'] == mt and r['freqs'] == list(fr) and r['s'] == s), None)
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
    if C is not None:
        L_ += ["## C. The thresholds against the error gain and the centres' mobility", "",
               f"The paradigm's metric, D‴.49's noise of rms s, {NP_C} periods from the fitted state: the smallest level (of "
               f"{', '.join(f'{s:g}' for s in LEV_C)}) at which each formed kernel left its feature ({LEFT_C:g} away, period-averaged), "
               "against k_e and the law's m (the centres' mobility Gmw), with √K's largest value along the orbit:", "",
               "| m | k_e | √K max | " + ' | '.join(f"A = {a:g}" for a in CT.ASTAR) + " |", "|---|---|---|---|---|---|"]
        for m in M_C:
            for ke in KE_C:
                rs = [r for r in C if r['m'] == m and r['ke'] == ke]
                if not rs:
                    continue
                top = max(r['s'] for r in rs)
                cells = [f"{x:g}" if x is not None else f"> {top:g}" for x in thresholds(rs)]
                L_.append(f"| {m:g} | {ke:g} | {rs[0]['sqrtK']:.2f} | " + ' | '.join(cells) + " |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_noise_resonance.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'noise_resonance_raw.npz')
    md = os.path.join(OUT, 'noise_resonance.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(old.get('A'), old.get('B'), old.get('B2'), old.get('C')))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else {'A', 'B', 'B2', 'C'}
    res = {k: v for k, v in old.items() if k not in parts}
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('B2', part_B2), ('C', part_C)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(res.get('A'), res.get('B'), res.get('B2'), res.get('C')))


if __name__ == '__main__':
    main()
