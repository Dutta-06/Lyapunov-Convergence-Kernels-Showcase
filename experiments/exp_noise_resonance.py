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
B4. The long horizon and commensurability (referee rounds 56 and 63): the loop's even and odd parts per s² and s at s ∈ S_B4 over NP_B4
   periods (≥ 5 time constants of the slow modes), at one tone (21.37; 21.4 = 107ω and 21.3 = 106.5ω, ω the reference's; 13.1 = 65.5ω
   and its neighbours 13.05, 13.15) and D‴.49's 16 tones: the even part's saturation, the odd part against the Floquet sums (F).
B5. The averaged loop against the exact second order's constant part (referee round 63): `second_order` at the tone's phase φ and
   φ + π/2 — their half-sum the constant part's response, their half-difference the 2ω part's — at 13.1 (2ω = 131ω) and its neighbours
   13.05, 13.15 (2ω = 130.5ω, 131.5ω), against `envelope_means`.
F. The Floquet multipliers of the loop linearised at the fit (referee round 63), per metric and law; the first order's period means by
   the Floquet recursion against B4's odd parts; the commensurate tones' saturation.
A2. Resonance against roll-off (referee round 63; quasi-static): the 16 tones' centre jitter with K and with K = 0 per law and band; the
   quasi-static peaks against √K_i.
C. The leaving levels: the levels of D‴.49's noise at which a formed kernel leaves its feature within NP_C periods (D‴.57's Part D),
   the paradigm's metric, against the error gain k_e (KE_C) and the law's m (M_C: the centres' mobility Gmw), under the noise's seeds
   SEEDS_C (the tones' frequencies and phases) — with √K_max's place against the noise band 10–30 rad/s. C2: controls at k_e = 10 —
   single tones, the band moved to 25–75 rad/s (above √K_max), m between 3 and 12, the Wasserstein–Fisher–Rao metric, tone pairs whose
   difference is commensurate with the reference or not (referee round 56) — under four seeds, the pairs at 0.02 … 0.08 and D‴.49's own
   pair with the draw's phases (referee round 63). C3: the anatomy of a leave — the formed kernels' period means through it.
C4, W2, EDGE, C6, C7 (referee round 69). C4: C3's band runs under the negated noise (the −n twins: the even part at the leave). W2: each
   draw's own second order over C's horizon (±S_W2·n), per configuration. EDGE: the noise-free loop's basin edge in a formed kernel's
   coefficient against m, k_e and the metric, the ejection's anatomy, and the drift's reach s* = (edge/⟨w₂⟩)^{1/2} against the first
   jumps. C6: C3's runs started with kernel 1's coefficient raised. C7: the edge under noise — single tones and the band's windows
   against the level.

Usage: python experiments/exp_noise_resonance.py [--only A,A2,B,B2,B3,B4,B5,C,C2,C3,F,C4,W2,EDGE,C6,C7,C8,C9,C9D] [--from-raw]
       → results/chaos/noise_resonance.md (+ _raw.npz).
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
TONES_C2 = (10.0, 13.1, 15.0, 15.7, 16.0, 21.3)                       # C2's single tones (50ω, 65.5ω, 75ω, 78.5ω, 80ω, 106.5ω)
PAIRS_C2 = ((18.9, 22.5), (18.9, 22.4), (18.9, 22.45), (20.0, 20.2))  # C2's pairs: differences 18ω, 17.5ω, 17.75ω, 1ω
M_SCAN_C2 = (3.5, 4.0, 4.5, 4.83, 5.5, 6.0, 8.0)                       # C2's laws between C's m = 3 and 12 (the bound 4.83)
NP_B4, MARKS_B4, S_B4 = 1200, (1, 3, 10, 30, 100, 300, 600, 900, 1200), (0.00125, 0.0025)
FREQS_B4 = ((21.37,), None, (21.4,), (21.3,), (13.1,), (13.05,), (13.15,))
BAND_B4 = (601, 1200)                                                  # B4: the 16 tones' bands over these periods
M_F, OM_COMM_F = (3.0, 3.5, 4.0, 4.5, 4.83, 5.5, 6.0, 8.0, 12.0), (10.0, 15.0, 21.4, 25.0, 29.4)
OM_B5, S_B5 = (13.05, 13.1, 13.15, 29.3), 0.00125                   # B5; the averaged loop also at S_B5 (referee round 69: its O(s⁴) out)
M_A2, BANDS_A2, OM_A2 = (2.0, 3.0, 4.0, 4.83, 6.0, 8.0, 12.0), ((10.0, 30.0), (25.0, 75.0)), np.linspace(4.0, 40.0, 721)
JUMP_C = 0.6                                                           # C, C2: a leave whose largest distance exceeds this is a jump
C3_CASES = ([(12.0, 10.0, s_, dict(seed=sd)) for s_ in (0.02, 0.03) for sd in range(1, 12)]
            + [(12.0, 10.0, 0.08, dict(freqs=(f,), seed=sd)) for f in (13.1, 15.0) for sd in (7, 1)])
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


def simulate(L, al0, c0, astar, ts, noise, metric, ke=None, x_init=None, om=None):
    from src.core import closed_loop as CL
    om = CT.OM if om is None else om                                     # referee round 76: the reference's frequency (EDGE)
    return CL.simulate(L, EPS, c0, CT.SIG, (CT.CSTAR, astar, CT.SIG), CT.KE if ke is None else ke,
                       (0.0, 1.0, om, 0.0, np.sqrt(2) * om), ts, noise=noise, th=CT.th_abs(astar), regx=True, al0=al0, metric=metric,
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
    amplitudes, the seed's phases, or the phases opts ph); with opts neg the noise negated (its phases moved by π: the −n twin, referee
    round 69)."""
    from src.core import closed_loop as CL
    if 'band' in opts:
        nz = CL.tones(s, *opts['band'], seed=opts.get('seed', 7))
    elif 'ph' in opts:
        fr = np.asarray(opts['freqs'], float)
        nz = dict(w=fr, ph=np.asarray(opts['ph'], float), amp=np.full(fr.size, s * np.sqrt(2.0 / fr.size)))
    else:
        nz = DZ.noise_of(s, {k: v for k, v in opts.items() if k in ('freqs', 'seed')})
    if opts.get('neg') and nz is not None:
        nz = dict(nz, ph=np.asarray(nz['ph'], float) + np.pi)
    if opts.get('shift') and nz is not None:                            # referee round 69 (C7): the noise from shift periods of the
        nz = dict(nz, ph=np.asarray(nz['ph'], float) + np.asarray(nz['w'], float) * opts['shift'] * 2 * np.pi / CT.OM)   # reference on
    return nz


def c_job(args):
    """C: the metric (opts metric, default the paradigm's) at the law's m and the error gain k_e, the noise of rms s (`c_noise`) from the
    fitted state for NP_C periods: per formed kernel the first period whose period-averaged centre is LEFT_C from its feature (−1: never),
    the largest period-averaged distance, the coefficients at the end."""
    m, ke, s = args[:3]
    opts = dict(args[3]) if len(args) > 3 else {}
    metric = opts.get('metric', 'paradigm')
    warnings.filterwarnings('ignore')
    t0 = time.time()
    r, L, al0, c0, idx, Wp, Cp = c_run(m, ke, s, opts)
    dev = np.abs(Cp - CT.CSTAR)
    o = orbit_terms(L, EPS, al0[idx], c0[idx], metric, np.sin(np.linspace(0, 2 * np.pi, 2001)))
    left = left_periods(dev, LEFT_C)
    w_left, pre = leave_stats(Wp, Cp, left)
    out = dict(m=m, ke=ke, s=s, opts={k_: (list(v) if isinstance(v, tuple) else v) for k_, v in opts.items()}, left=left,
               maxdev=dev.max(0).tolist(), w_end=[float(r['w'][-1, i]) for i in idx], sqrtK=float(np.sqrt(o['K'].max())),
               c_end=[float(r['c'][-1, i]) for i in idx], w_left=w_left, pre=pre, secs=time.time() - t0)
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
    ew, ow = ((Wp + Wm) / 2 - CT.ASTAR) / s ** 2, ((Wp - Wm) / 2) / s
    ec, oc = ((Cp + Cm) / 2 - CT.CSTAR) / s ** 2, ((Cp - Cm) / 2) / s
    out = dict(freqs=None if freqs is None else list(freqs), s=s, even_w=ew[mk].tolist(), odd_w=ow[mk].tolist(), even_c=ec[mk].tolist(),
               odd_c=oc[mk].tolist(), odd_c_per=oc[:12, 0].tolist(), even_w_all=ew.tolist(), odd_w_all=ow.tolist(),
               even_c_all=ec.tolist(), secs=time.time() - t0)
    print(f"  B4 {noise_label(out['freqs'])} s {s:g}: even/s² at {MARKS_B4[-1]}: {', '.join(f'{x:+.4g}' for x in out['even_w'][-1])} "
          f"[{out['secs']:.0f}s]", flush=True)
    return out


def monodromy(L, eps, al, c, metric, ke, om=None, ph=0.0):
    """F (referee round 63): the loop linearised at the exact fit along x_m(t) — `second_order`'s first order, the state (e₁, w₁, c₁) of
    the formed kernels, e₁ the true error and d₁ = e₁ + n the measured one — over one period T = 2π/ω_ref from t = 0: the monodromy
    M = Φ(T, 0), P = (1/T)∫₀ᵀΦ(t, 0)dt and, for the forcing n = Re(e^{i(ωt + ph)}), G, the state at T from zero, and H, its period mean
    (complex: the responses to cos and sin); the first order's period means follow by the recursion in `floquet_means`. Returns (M, P)
    or (M, P, G, H)."""
    from scipy.integrate import solve_ivp
    w, mu, _, v, _ = curve(L, eps, al, metric)
    nk = c.size
    d = 1 + 2 * nk
    T = 2 * np.pi / CT.OM

    def A_b(t):
        u = (np.sin(CT.OM * t) - c) / CT.SIG
        ph_ = np.exp(-u ** 2 / 2)
        d1 = ph_ * u / CT.SIG
        A = np.zeros((d, d))
        A[0, 0], A[0, 1:1 + nk], A[0, 1 + nk:] = -ke, -ph_, -w * d1
        A[1:1 + nk, 0], A[1 + nk:, 0] = mu * ph_, v * d1
        return A, np.r_[(w * d1).sum() - ke, mu * ph_, v * d1]          # the forcing's column: −f̂′ − k_e, μφ, v∂_cφ
    nf = om is not None

    def rhs(t, y):
        A, b = A_b(t)
        Phi = y[:d * d].reshape(d, d)
        out = [(A @ Phi).ravel(), Phi.ravel()]
        if nf:
            zc, zs = y[2 * d * d:2 * d * d + d], y[2 * d * d + 2 * d:2 * d * d + 3 * d]
            out += [A @ zc + b * np.cos(om * t + ph), zc, A @ zs + b * np.sin(om * t + ph), zs]
        return np.concatenate(out)
    y0 = np.zeros(2 * d * d + (4 * d if nf else 0))
    y0[:d * d] = np.eye(d).ravel()
    sol = solve_ivp(rhs, (0.0, T), y0, method='DOP853', rtol=1e-11, atol=1e-14, max_step=min(0.05, 0.2 / om) if nf else 0.05)
    y = sol.y[:, -1]
    M, P = y[:d * d].reshape(d, d), y[d * d:2 * d * d].reshape(d, d) / T
    if not nf:
        return M, P
    o = 2 * d * d
    return M, P, y[o:o + d] + 1j * y[o + 2 * d:o + 3 * d], (y[o + d:o + 2 * d] + 1j * y[o + 3 * d:o + 4 * d]) / T


def floquet_means(M, P, G, H, om, NP):
    """The first order's period means over NP periods from z₀ = 0 (the loop starts at the fit, e = 0) by the Floquet recursion: over the
    period j the tone is Re(e^{iωjT}·e^{i(ωt′ + ph)}), so z_{j+1} = Mz_j + Re(e^{iωjT}G) and ⟨z⟩_j = Pz_j + Re(e^{iωjT}H). On each slow
    mode of M (multiplier ρ) the increments rotate by ωT and add up to |δ||ρ^j − e^{iωjT}|/|ρ − e^{iωT}|: bounded, at commensurability
    the offset |δ|/(1 − ρ) reached with the mode's time constant."""
    T = 2 * np.pi / CT.OM
    z, out = np.zeros(M.shape[0]), []
    for j in range(NP):
        rot = np.exp(1j * om * j * T)
        out.append(P @ z + np.real(rot * H))
        z = M @ z + np.real(rot * G)
    return np.array(out)


def f_job(args):
    """F: per metric and law (m) the linearised loop's Floquet multipliers (by modulus), the time constants −1/ln|ρ| (periods) and each
    mode's largest component (e, w_i, c_i); for the paradigm's metric at m = 12 also the first order's period means by the Floquet
    recursion for B4's single tones (per unit rms, NP_B4 periods: the loop's odd part per s) and, at the commensurate tones OM_COMM_F, the
    period means of the first three periods and the saturation P(I − M)⁻¹Re(G) + Re(H) (per unit rms)."""
    metric, m = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law(m)
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    al, c = al0[idx], c0[idx]
    M, P = monodromy(L, EPS, al, c, metric, CT.KE)
    ev, V = np.linalg.eig(M)
    o = np.argsort(-np.abs(ev))
    names = ['e'] + [f'w{q + 1}' for q in range(c.size)] + [f'c{q + 1}' for q in range(c.size)]
    out = dict(metric=metric, m=m, rho=np.abs(ev[o]).tolist(), imag=float(np.abs(ev.imag).max()),
               tau=[float(-1 / np.log(abs(x))) for x in ev[o] if 1e-300 < abs(x) < 1],
               mode=[names[int(np.argmax(np.abs(V[:, k])))] for k in o])
    if metric == 'paradigm' and m == 12.0:
        sums, comm = {}, {}
        for fr in [f for f in FREQS_B4 if f is not None]:
            nz = DZ.noise_of(1.0, dict(freqs=fr))
            a = float(nz['amp'][0])
            _, _, G, H = monodromy(L, EPS, al, c, metric, CT.KE, fr[0], float(nz['ph'][0]))
            Z = floquet_means(M, P, a * G, a * H, fr[0], NP_B4)
            sums[f"{fr[0]:g}"] = dict(w=Z[:, 1:1 + c.size].tolist(), c=Z[:, 1 + c.size:].tolist())
        for om in OM_COMM_F:
            nz = DZ.noise_of(1.0, dict(freqs=(om,)))
            a = float(nz['amp'][0])
            _, _, G, H = monodromy(L, EPS, al, c, metric, CT.KE, om, float(nz['ph'][0]))
            Z = floquet_means(M, P, a * G, a * H, om, 3)
            sat = P @ np.linalg.solve(np.eye(M.shape[0]) - M, np.real(a * G)) + np.real(a * H)
            comm[f"{om:g}"] = dict(first=Z.tolist(), sat=sat.tolist())
        out.update(sums=sums, comm=comm)
    out['secs'] = time.time() - t0
    print(f"  F {metric} m {m:g}: |ρ| {', '.join(f'{x:.5f}' for x in out['rho'][:4])}, τ {', '.join(f'{x:.1f}' for x in out['tau'][:3])} "
          f"[{out['secs']:.0f}s]", flush=True)
    return out


def b5_job(args):
    """B5 (referee round 63): at one tone the exact second order (`second_order`) under the tone's phase φ and φ + π/2: the second order
    is the response to the products of the first order's tones, ½Re(AB̄) (φ-free) + ½Re(ABe^{2iφ}), so the half-sum of the two runs is the
    constant part's response and the half-difference the 2ω part's; against them the averaged loop (`envelope_means` at s = S_B[−1] and
    at S_B5, per s²), at MARKS_B."""
    metric, om = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law()
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    nz = DZ.noise_of(1.0, dict(freqs=(om,)))
    _, _, W2a, C2a, sa = second_order(L, EPS, al0[idx], c0[idx], metric, nz, CT.KE, NP_B)
    _, _, W2b, C2b, sb = second_order(L, EPS, al0[idx], c0[idx], metric, dict(nz, ph=np.asarray(nz['ph']) + np.pi / 2), CT.KE, NP_B)
    s = S_B[-1]
    We, Ce, se = envelope_means(L, EPS, metric, DZ.noise_of(s, dict(freqs=(om,))), CT.KE, NP_B)
    Wf, Cf, sf = envelope_means(L, EPS, metric, DZ.noise_of(S_B5, dict(freqs=(om,))), CT.KE, NP_B)
    mk = [q - 1 for q in MARKS_B]
    out = dict(metric=metric, om=om, status=[sa, sb, se, sf], w2=W2a[mk].tolist(), w2_shift=W2b[mk].tolist(), c2=C2a[mk].tolist(),
               c2_shift=C2b[mk].tolist(), w_env=((We[mk] - CT.ASTAR) / s ** 2).tolist(), c_env=((Ce[mk] - CT.CSTAR) / s ** 2).tolist(),
               w_env_small=((Wf[mk] - CT.ASTAR) / S_B5 ** 2).tolist(), c_env_small=((Cf[mk] - CT.CSTAR) / S_B5 ** 2).tolist(),
               secs=time.time() - t0)
    print(f"  B5 {metric} {om:g}: ⟨w₂⟩ at {MARKS_B[-1]} {', '.join(f'{x:+.4g}' for x in out['w2'][-1])}, shifted "
          f"{', '.join(f'{x:+.4g}' for x in out['w2_shift'][-1])} [{out['secs']:.0f}s]", flush=True)
    return out


def a2_rows():
    """A2 (referee round 63; quasi-static — A's oscillator law —, the paradigm's metric): per law (M_A2) each kernel's √K_i (its largest
    (μ_iφ_i² + v_iw_i(∂_cφ_i)²)^{1/2} on the orbit) and the frequency on OM_A2 where its centre jitter per unit amplitude — the rms over
    the orbit of |v_i∂_cφ_iY|/√2, Y = (iω − f̂′)/(K − ω² + ik_eω) — peaks; per band (BANDS_A2, seed 7's draw) the 16 tones' centre jitter
    per unit rms, and with K = 0 — the K-independent roll-off —, their ratio the resonance factor."""
    from src.core import closed_loop as CL
    xm = np.sin(np.linspace(0, 2 * np.pi, 2001)[:-1])
    rows = []
    for m in M_A2:
        L = law(m)
        al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
        o = orbit_terms(L, EPS, al0[idx], c0[idx], 'paradigm', xm)
        sqk = np.sqrt((o['mu'] * o['ph'] ** 2 + o['v'] * o['w'] * o['d1'] ** 2).max(0))
        g = o['v'][None, :] * o['d1']                                    # (orbit, kernel): v_i∂_cφ_i
        jit = lambda Y, amp: np.sqrt(np.mean((np.abs(g[:, :, None] * Y[:, None, :] * amp[None, None, :]) ** 2 / 2).sum(2), axis=0))
        per = np.array([jit(transfer(o, np.array([w_]), CT.KE)[0], np.ones(1)) for w_ in OM_A2])     # (ω, kernel)
        for band in (BANDS_A2 if m in (3.0, 12.0) else BANDS_A2[:1]):
            nz = CL.tones(1.0, *band, seed=7)
            om, amp = np.asarray(nz['w'], float), np.asarray(nz['amp'], float)
            Y0 = (1j * om[None, :] - o['f1'][:, None]) / (-om[None, :] ** 2 + 1j * CT.KE * om[None, :])
            rows.append(dict(m=m, band=list(band), sqrtK=sqk.tolist(), peak=OM_A2[np.argmax(per, axis=0)].tolist(),
                             jit=jit(transfer(o, om, CT.KE)[0], amp).tolist(), jit0=jit(Y0, amp).tolist()))
    return rows


def draw_pair():
    """D‴.49's draw (seed 7): its two tones nearest 18.9 and 22.5 rad/s and their phases — the difference 18.0019ω_ref, a beat of 521
    periods against the orbit's 18th harmonic (referee rounds 56 and 63)."""
    from src.core import closed_loop as CL
    nz = CL.tones(1.0, 10.0, 30.0, seed=7)
    j = [int(np.argmin(np.abs(nz['w'] - f))) for f in (18.9, 22.5)]
    return tuple(float(nz['w'][q]) for q in j), tuple(float(nz['ph'][q]) for q in j)


def c_run(m, ke, s, opts):
    """C's loop: the metric (opts metric, default the paradigm's) at the law's m and the error gain k_e, the noise of rms s (`c_noise`),
    from the fitted state for NP_C periods (8 outputs a period). Returns the run, the formed kernels' indices and their period means of
    w and c."""
    L = law(m)
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    if 'rise' in opts:                                                   # referee round 69: kernel k's coefficient raised by Δ at the start
        al0 = raised(L, al0, idx, int(opts['rise'][0]), float(opts['rise'][1]))
    k = 8
    ts = np.linspace(0, NP_C * 2 * np.pi / CT.OM, NP_C * k + 1)
    r = simulate(L, al0, c0, CT.ASTAR, ts, c_noise(s, opts), opts.get('metric', 'paradigm'), ke)
    per = lambda z: z[:-1].reshape(NP_C, k, c0.size).mean(1)[:, idx]
    return r, L, al0, c0, idx, per(r['w']), per(r['c'])


def leave_stats(Wp, Cp, left):
    """Per formed kernel at its leave (the first period LEFT_C from its feature): its period-mean coefficient there, and the mean signed
    offset of its period-mean centre over the periods leave − 30 … leave − 6 (None where it never left)."""
    out_w, out_pre = [], []
    for q, p in enumerate(left):
        if p < 0:
            out_w.append(None), out_pre.append(None)
            continue
        i = p - 1
        lo, hi = max(i - 30, 0), max(i - 5, 1)
        out_w.append(float(Wp[i, q])), out_pre.append(float(np.mean(Cp[lo:hi, q] - CT.CSTAR[q])))
    return out_w, out_pre


def c3_job(args):
    """C3 (referee round 63): the anatomy of a leave — C's loop at (m, k_e, s, opts) with the formed kernels' period means of w and c
    over NP_C periods, the leave periods, the largest distances, the coefficient at the leave and the mean offset before it."""
    m, ke, s, opts = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    r, L, al0, c0, idx, Wp, Cp = c_run(m, ke, s, dict(opts))
    dev = np.abs(Cp - CT.CSTAR)
    left = left_periods(dev, LEFT_C)
    w_left, pre = leave_stats(Wp, Cp, left)
    out = dict(m=m, ke=ke, s=s, opts={k_: (list(v) if isinstance(v, tuple) else v) for k_, v in dict(opts).items()}, left=left,
               maxdev=dev.max(0).tolist(), w_left=w_left, pre=pre, W=Wp.tolist(), C=Cp.tolist(), secs=time.time() - t0)
    print(f"  C3 m {m:g} k_e {ke:g} s {s:g} {out['opts']}: left {left} [{out['secs']:.0f}s]", flush=True)
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
    return ([(mt, (om,)) for mt in METRICS for om in OM_B] + [(mt, None) for mt in METRICS]
            + [('paradigm', (om,)) for om in OM_B5 if om not in OM_B])


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


def part_B3(old=None):
    """B3's jobs; those already in `old` (the raw's B3, keyed by metric and tones) are kept, not rerun (referee round 63 adds 13.1's
    neighbours)."""
    have = {(r['metric'], tuple(r['freqs']) if r['freqs'] is not None else None): r for r in (old if old is not None else [])}
    todo = [j for j in b3_jobs() if j not in have]
    new = dict(zip(todo, _pool_map(b3_job, todo))) if todo else {}
    return [have.get(j, new.get(j)) for j in b3_jobs()]


def part_C():
    return _pool_map(c_job, [(m, ke, s, dict(seed=sd)) for m in M_C for ke in KE_C for s in LEV_C for sd in SEEDS_C])


def c2_cases():
    """C2's cases (referee rounds 56 and 63): (label, m, k_e, s, opts), each under SEEDS_C but D‴.49's own pair (the draw's phases)."""
    fr, ph = draw_pair()
    return ([('one tone', 12.0, 10.0, s_, dict(freqs=(f,), seed=sd)) for f in TONES_C2 for s_ in (0.04, 0.06, 0.08) for sd in SEEDS_C]
            + [('tone pair', 12.0, 10.0, s_, dict(freqs=p, seed=sd)) for p in PAIRS_C2 for s_ in (0.02, 0.04, 0.06, 0.08) for sd in SEEDS_C]
            + [("the draw's pair", 12.0, 10.0, s_, dict(freqs=fr, ph=ph)) for s_ in (0.02, 0.04, 0.06, 0.08)]
            + [('band 25–75', m_, 10.0, s_, dict(band=(25.0, 75.0), seed=sd)) for m_ in (3.0, 12.0) for s_ in (0.02, 0.04, 0.06, 0.08)
               for sd in SEEDS_C]
            + [('m', m_, 10.0, s_, dict(seed=sd)) for m_ in M_SCAN_C2 for s_ in (0.02, 0.03, 0.04, 0.06, 0.08) for sd in SEEDS_C]
            + [('WFR', 12.0, ke_, s_, dict(metric='wfr', seed=sd)) for ke_ in (5.0, 10.0) for s_ in (0.02, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32)
               for sd in SEEDS_C])


def part_C2():
    cases = c2_cases()
    return [dict(r, label=lab) for (lab, *_), r in zip(cases, _pool_map(c_job, [cs[1:] for cs in cases]))]


def part_C3():
    fr, ph = draw_pair()
    return _pool_map(c3_job, list(C3_CASES) + [(12.0, 10.0, s_, dict(freqs=fr, ph=ph)) for s_ in (0.02, 0.04, 0.06, 0.08)])


def part_F():
    return _pool_map(f_job, [('paradigm', m) for m in M_F] + [('wfr', 12.0)])


def part_B5():
    return _pool_map(b5_job, [('paradigm', om) for om in OM_B5])


def part_A2():
    return a2_rows()


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
            + " (the relative error mostly ∝ s²: the remainders per doubling below); the centres' even part over s²⟨c₂⟩: " + '; '.join(f"{rng(q_even_c[q])} at s = {q:g}"
                                                                                                 for q in S_B3)
            + "; the averaged loop over it at s = " + f"{S_B3[-1]:g}: "
            + '; '.join(f"{k[0]}, {k[1]}: {rng(v)}" for k, v in q_avg.items()) + ".", ""]
    sc = []
    for r in sorted(B3, key=key):                                       # referee round 63: the remainders' scaling per doubling of s
        rows = sorted(r['rows'], key=lambda q: q['s'])
        ro = [np.array(q['odd_w'][-1]) - q['s'] * np.array(r['w1'][-1]) for q in rows]
        re_ = [np.array(q['even_w'][-1]) - q['s'] ** 2 * np.array(r['w2'][-1]) for q in rows]
        f = lambda R_: '; '.join(', '.join(f"×{a / b:.3g}" for a, b in zip(R_[i + 1], R_[i])) for i in range(len(R_) - 1))
        sc.append(f"{r['metric']}, {noise_label(r['freqs'])}: odd {f(ro)}; even {f(re_)}")
    out += [f"The remainders at period {MARKS_B[-1]} per doubling of s (s = " + ' → '.join(f"{q:g}" for q in S_B3) + "), per kernel — the odd part's "
            "less s⟨w₁⟩ (×8 if cubic) and the even part's less s²⟨w₂⟩ (×16 if quartic): " + ' | '.join(sc) + ".", ""]
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


def f_lines(F, B4=None):
    """F's report: the multipliers per metric and law; the Floquet sums against B4's odd parts; the commensurate tones."""
    if F is None:
        return []
    out = ["## F. The linearised loop's Floquet multipliers (referee round 63)", "",
           "The loop linearised at the fit along the orbit (e, the formed kernels' w and c) over one period of the reference: the monodromy's "
           "multipliers |ρ| by modulus, the slow modes' time constants −1/ln|ρ| in periods and each mode's largest component:", "",
           "| metric | m | \\|ρ\\| | τ of the three slowest (periods) | their largest components | the largest imaginary part |",
           "|---|---|---|---|---|---|"]
    for r in sorted(F, key=lambda r: (r['metric'] != 'paradigm', r['m'])):
        out.append(f"| {r['metric']} | {r['m']:g} | {', '.join(f'{x:.5f}' for x in r['rho'][:5])} | {', '.join(f'{x:.1f}' for x in r['tau'][:3])} | "
                   f"{', '.join(r['mode'][:3])} | {r['imag']:.1e} |")
    base = next((r for r in F if r['metric'] == 'paradigm' and r['m'] == 12.0 and 'sums' in r), None)
    if base is not None and B4 is not None:
        out += ["", "The first order's period means by the Floquet recursion (z_{j+1} = Mz_j + Re(e^{iωjT}G), ⟨z⟩_j = Pz_j + Re(e^{iωjT}H)) "
                f"against the loop's odd part per s (B4, all {NP_B4} periods, the coefficients): the largest difference over the largest "
                "Floquet mean, per s, and its growth per doubling of s (×4: the odd part's O(s³) remainder):", ""]
        cells = []
        for fr in [f for f in FREQS_B4 if f is not None]:
            sm = base['sums'].get(f"{fr[0]:g}")
            rows = sorted([r for r in B4 if r['freqs'] == list(fr) and 'odd_w_all' in r], key=lambda r: r['s'])
            if sm is None or not rows:
                continue
            fz = np.array(sm['w'])
            dif = [float(np.abs(np.array(r['odd_w_all']) - fz).max() / np.abs(fz).max()) for r in rows]
            cells.append(f"{fr[0]:g}: " + ', '.join(f"{d:.2e} (s = {r['s']:g})" for d, r in zip(dif, rows))
                         + (f", ×{dif[1] / dif[0]:.2f}" if len(dif) > 1 and dif[0] > 0 else ''))
        out.append('; '.join(cells) + ".")
    if base is not None:
        T = 2 * np.pi / CT.OM
        out += ["", "At the commensurate tones (ωT ∈ 2πℤ: the increments add on every slow mode and saturate at P(I − M)⁻¹Re(G) + Re(H)), "
                "per unit rms: the first kernel's centre period means at the periods 1–3, the coefficients' saturation and the centres':", "",
                "| ω | ω/ω_ref | centre 1, periods 1–3 | the coefficients' saturation | the centres' saturation |", "|---|---|---|---|---|"]
        nk = len(CT.CSTAR)
        for om, cm in sorted(base['comm'].items(), key=lambda kv: float(kv[0])):
            Z, sat = np.array(cm['first']), np.array(cm['sat'])
            out.append(f"| {om} | {float(om) / CT.OM:g} | {', '.join(f'{x:+.4f}' for x in Z[:, 1 + nk])} | "
                       f"{', '.join(f'{x:+.4f}' for x in sat[1:1 + nk])} | {', '.join(f'{x:+.4f}' for x in sat[1 + nk:])} |")
    out.append("")
    return out


def b4_lines(B4):
    """B4's report (referee round 63: NP_B4 periods, ≥ 5 time constants): the even part per s² and the odd part per s at the marks,
    the 16 tones' bands over BAND_B4, the first kernel's centre per period."""
    if B4 is None:
        return []
    out = ["## B4. The long horizon and commensurability", "",
           f"The paradigm's metric, the loop under s·n and −s·n for {NP_B4} periods: the even part of the coefficients' period means per s² "
           f"and the odd part per s at the periods {', '.join(str(m) for m in MARKS_B4)} (kernel by kernel), at s = "
           f"{', '.join(f'{q:g}' for q in S_B4)}; the reference's ω = {CT.OM:g} (21.4 = 107ω, 21.3 = 106.5ω, 13.1 = 65.5ω: its 2ω = 131ω; "
           "13.05, 13.15: 2ω = 130.5ω, 131.5ω):", "",
           "| noise | s | even part / s² (w), at the periods | odd part / s (w), at the periods | the centres' even part / s² at the "
           "last period | the first kernel's centre, odd part / s, periods 1–12 |", "|---|---|---|---|---|---|"]
    for r in sorted(B4, key=lambda r: (r['freqs'] is None, r['freqs'] or [0], r['s'])):
        ev = '; '.join(', '.join(f"{x:+.3g}" for x in row) for row in r['even_w'])
        od = '; '.join(', '.join(f"{x:+.3g}" for x in row) for row in r['odd_w'])
        out.append(f"| {noise_label(r['freqs'])} | {r['s']:g} | {ev} | {od} | {', '.join(f'{x:+.3g}' for x in r['even_c'][-1])} | "
                   f"{', '.join(f'{x:+.3f}' for x in r['odd_c_per'])} |")
    band = [r for r in B4 if r['freqs'] is None and 'even_w_all' in r]
    if band:
        a, b = BAND_B4
        out += ["", f"D‴.49's 16 tones over the periods {a}–{b}: per kernel the least and largest even part per s² of the coefficients "
                "and of the centres, and the odd part per s of the coefficients:", ""]
        for r in sorted(band, key=lambda r: r['s']):
            rng = lambda key: '; '.join(f"[{lo:+.4g}, {hi:+.4g}]" for lo, hi in zip(np.array(r[key])[a - 1:b].min(0),
                                                                               np.array(r[key])[a - 1:b].max(0)))
            out.append(f"s = {r['s']:g}: w {rng('even_w_all')}; c {rng('even_c_all')}; odd w {rng('odd_w_all')}.")
    out.append("")
    return out


def b5_lines(B5):
    """B5's report: the exact second order's constant and 2ω parts (the phase shifted by π/2) against the averaged loop."""
    if B5 is None:
        return []
    out = ["## B5. The averaged loop against the second order's constant part (referee round 63)", "",
           f"One tone, {NP_B} periods: the exact second order's ⟨w₂⟩ (`second_order`) at the tone's phase φ and at φ + π/2 — their half-sum "
           "the response to the constant part of the forcing (φ-free), their half-difference the 2ω part's (∝ e^{2iφ}) —, and the averaged "
           f"loop's Δw/s² at s = {S_B5:g} (its own O(s⁴) negligible; referee round 69) and at {S_B[-1]:g} over the exact second order and "
           f"over its constant part, at the periods {', '.join(str(m) for m in MARKS_B)}:", "",
           f"| metric | ω | 2ω/ω_ref | ⟨w₂⟩ (exact) | the constant part | the 2ω part | averaged / exact, s = {S_B5:g} | averaged / the "
           f"constant part, s = {S_B5:g} | averaged / the constant part, s = {S_B[-1]:g} |",
           "|---|---|---|---|---|---|---|---|---|"]
    lo_ = []
    for r in sorted(B5, key=lambda r: (r['metric'] != 'paradigm', r['om'])):
        a, b, e = np.array(r['w2']), np.array(r['w2_shift']), np.array(r['w_env'])
        es = np.array(r['w_env_small']) if 'w_env_small' in r else None
        const, two = (a + b) / 2, (a - b) / 2
        f = lambda X: '; '.join(', '.join(f"{x:+.4g}" for x in row) for row in X)
        g = lambda X: '; '.join(', '.join(f"{x:.4f}" for x in row) for row in X)
        out.append(f"| {r['metric']} | {r['om']:g} | {2 * r['om'] / CT.OM:g} | {f(a)} | {f(const)} | {f(two)} | "
                   + ('—' if es is None else g(es / a)) + ' | ' + ('—' if es is None else g(es / const)) + f" | {g(e / const)} |")
        if es is not None:
            lo_.append((r['om'], float((es / const).min()), float((es / const).max()), float((es / a).min()), float((es / a).max())))
    if lo_:
        out += ["", f"At s = {S_B5:g} the averaged loop over the constant part: " + '; '.join(
            f"{om:g}: {a:.4f}–{b:.4f} (over the exact second order {c:.4f}–{d:.4f})" for om, a, b, c, d in lo_) + "."]
    out.append("")
    return out


def sens(K, ke, om):
    """The adaptation loop's sensitivity |S(iω)| = |iω(iω + k_e)|/|K − ω² + ik_eω| — per tone, at a frozen K, A2's factor (the plant's
    slope cancels in the ratio of the quasi-static centre jitter with K to that at K = 0)."""
    om = np.asarray(om, float)
    return np.abs(1j * om * (1j * om + ke)) / np.abs(K - om ** 2 + 1j * ke * om)


def bode_integral(K, ke):
    """∫₀^∞ ln|S(iω)| dω by quadrature (Bode's sensitivity integral: 0 for the open loop K/(s(s + k_e)), relative degree 2, no pole
    in the open right half-plane)."""
    from scipy.integrate import quad
    f = lambda w: float(np.log(sens(K, ke, w)))
    return quad(f, 0.0, 1.0, limit=200)[0] + quad(f, 1.0, np.inf, limit=400)[0]


def a2_lines(A2):
    """A2's report: the quasi-static peaks against √K_i (a maximum at the grid's lower end marked); the 16 tones' jitter with and without
    K — their ratio the adaptation loop's sensitivity |S| over the band — per law and band, and the band's effect split into its roll-off
    and its resonance; |S| at the bound and at √K."""
    if A2 is None:
        return []
    out = ["## A2. Resonance against roll-off (quasi-static; referee round 63)", "",
           "A's oscillator law, the paradigm's metric: per law each kernel's √K_i and the frequency at which its centre jitter per unit "
           f"amplitude peaks (on ω = {OM_A2[0]:g} … {OM_A2[-1]:g}; * a maximum at the grid's lower end, the peak below it); seed 7's 16 tones "
           "in the band, per unit rms: the centre jitter (rms over the orbit and the tones), the same with K = 0 (the roll-off alone), and their "
           "ratio, the resonance factor — per tone and orbit point the adaptation loop's sensitivity |S| = |iω(iω + k_e)|/|K − ω² + ik_eω| "
           "(the plant's slope cancels), here averaged over the band and the orbit:", "",
           "| m | band | √K_i | the peaks (/ √K_i) | jitter | jitter at K = 0 | resonance factor |", "|---|---|---|---|---|---|---|"]
    for r in sorted(A2, key=lambda r: (r['band'], r['m'])):
        out.append(f"| {r['m']:g} | {'–'.join(f'{x:g}' for x in r['band'])} | {', '.join(f'{x:.2f}' for x in r['sqrtK'])} | "
                   + ', '.join(f"{p:.2f}{'*' if p <= OM_A2[0] else ''} ({p / k:.3f})" for p, k in zip(r['peak'], r['sqrtK']))
                   + f" | {', '.join(f'{x:.4f}' for x in r['jit'])} | "
                   f"{', '.join(f'{x:.4f}' for x in r['jit0'])} | {', '.join(f'{a / b:.3f}' for a, b in zip(r['jit'], r['jit0']))} |")
    lines = []
    for m in sorted({r['m'] for r in A2}):
        lo = next((r for r in A2 if r['m'] == m and r['band'] == list(BANDS_A2[0])), None)
        hi = next((r for r in A2 if r['m'] == m and r['band'] == list(BANDS_A2[1])), None)
        if lo is None or hi is None:
            continue
        tot = [a / b for a, b in zip(lo['jit'], hi['jit'])]
        roll = [a / b for a, b in zip(lo['jit0'], hi['jit0'])]
        lines.append(f"m = {m:g}: the jitter falls by {', '.join(f'{x:.2f}' for x in tot)} from {'–'.join(f'{x:g}' for x in BANDS_A2[0])} to "
                     f"{'–'.join(f'{x:g}' for x in BANDS_A2[1])}, of which the roll-off {', '.join(f'{x:.2f}' for x in roll)} and the resonance "
                     f"{', '.join(f'{a / b:.2f}' for a, b in zip(tot, roll))}")
    if lines:
        out += ["", '; '.join(lines) + "."]
    k1 = sorted({(r['m'], r['sqrtK'][0]) for r in A2})
    if k1:
        bnd = min(k1, key=lambda x: abs(x[1] - 10.0))
        oms = (10.0, 15.0, 20.0, 30.0)
        out += ["", f"The sensitivity (k_e = {CT.KE:g}, the first kernel's K): |S|² − 1 = K(2ω² − K)/|K − ω² + ik_eω|², so |S| > 1 at every ω > "
                    f"√(K/2) whatever k_e — the amplification above is the price of the attenuation below (∫₀^∞ ln|S(iω)|dω = 0, Bode's "
                    f"integral; by quadrature {max(abs(bode_integral(k ** 2, CT.KE)) for _, k in k1):.1e} at most over the laws) —; |S(√K)| = "
                    f"(1 + K/k_e²)^{{1/2}} = {', '.join(f'{float(sens(k ** 2, CT.KE, k)):.3f}' for _, k in k1)} at m = "
                    f"{', '.join(f'{m:g}' for m, _ in k1)}; at m = {bnd[0]:g} (√K₁ = {bnd[1]:.2f}) |S| = "
                    f"{', '.join(f'{x:.2f}' for x in sens(bnd[1] ** 2, CT.KE, oms))} at ω = {', '.join(f'{w:g}' for w in oms)}."]
    out.append("")
    return out


def leaves_by_level(rs, jump=JUMP_C):
    """Per level: the leaves (kernel-runs that left), those that jumped (largest distance > jump) and the leave periods."""
    out = {}
    for r in rs:
        d = out.setdefault(r['s'], dict(n=0, leaves=0, jumps=0, periods=[]))
        d['n'] += 1
        for q, p in enumerate(r['left']):
            if p >= 0:
                d['leaves'] += 1
                d['jumps'] += int(r['maxdev'][q] > jump)
                d['periods'].append(p)
    return out


def c2_groups(C2):
    """C2's runs grouped by control, law, gain, noise (the seeds pooled) and level."""
    gs = {}
    for r in C2:
        o = dict(r['opts'])
        o.pop('seed', None)
        key = (r['label'], o.get('metric', 'paradigm'), r['m'], r['ke'], c2_noise_label(o), r['s'])
        gs.setdefault(key, []).append(r)
    return gs


def c2_noise_label(o):
    if 'band' in o:
        return '16 tones in ' + '–'.join(f"{x:g}" for x in o['band'])
    if 'freqs' in o:
        fr = list(o['freqs'])
        lab = noise_label(fr) if len(fr) == 1 else 'tones ' + '/'.join(f"{f:g}" for f in fr)
        if len(fr) == 2:
            lab += f" (difference {abs(fr[1] - fr[0]) / CT.OM:.4g}ω)"
        return lab + (", the draw's phases" if 'ph' in o else '')
    return "D‴.49's 16 tones"


def c2_lines(C2):
    """C2's report (referee round 63): per control, law, gain, noise and level, over the seeds: the runs, the kernels that left (0.3 away)
    and those that jumped (> JUMP_C), the leave periods and the largest distances."""
    if C2 is None:
        return []
    L12 = law(12.0)
    al_, c_, ix_ = DZ.fitted_state(L12, EPS, CT.CSTAR, CT.ASTAR)
    o_ = orbit_terms(L12, EPS, al_[ix_], c_[ix_], 'paradigm', np.sin(np.linspace(0, 2 * np.pi, 2001)))
    S_ = float((o_['w'] ** 2 * o_['d1'] ** 2).sum(1).max())
    out = ["## C2. Controls (k_e and m as stated)", "",
           f"The same loop and horizon, each case under the seeds {', '.join(map(str, SEEDS_C))} (the tones' phases; for a band also its "
           "frequencies) unless the draw's phases are given: one tone; tone pairs of total rms s — their difference a multiple of the "
           "reference's ω (18ω, 1ω), half of one (17.5ω) or between (17.75ω), and D‴.49's own pair 18.9015/22.5019 (18.0019ω) with the "
           "draw's phases —; the 16 tones moved to 25–75 rad/s (above √K_max); m between 3 and 12 (the design bound m ≲ "
           f"ω_lo²/(G max_t Σ_iw_i²(∂_cφ_i)²) = {100.0 / (L12.G * S_):.2f} at ω_lo = 10); the Wasserstein–Fisher–Rao metric to s = 0.32. Per "
           f"group the kernels that left their feature ({LEFT_C:g} away, period-averaged) — 'jumped' if their largest distance exceeds "
           f"{JUMP_C:g} (out of its feature) —, the leave periods and the largest distances over the runs:", "",
           "| control | metric | m | k_e | noise | s | runs | left: per kernel (jumped) | leave periods | largest distances (per kernel, max over runs) |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    order = {'one tone': 0, 'tone pair': 1, "the draw's pair": 2, 'band 25–75': 3, 'm': 4, 'WFR': 5}
    for key, rs in sorted(c2_groups(C2).items(), key=lambda kv: (order.get(kv[0][0], 9), kv[0][1:])):
        lab, mt, m, ke, nz, s = key
        nq = len(rs[0]['left'])
        lefts = [sum(r['left'][q] >= 0 for r in rs) for q in range(nq)]
        jumps = [sum(r['left'][q] >= 0 and r['maxdev'][q] > JUMP_C for r in rs) for q in range(nq)]
        per = sorted(p for r in rs for p in r['left'] if p >= 0)
        out.append(f"| {lab} | {mt} | {m:g} | {ke:g} | {nz} | {s:g} | {len(rs)} | "
                   + ', '.join(f"{a} ({b})" for a, b in zip(lefts, jumps)) + " | "
                   + ('—' if not per else (f"{per[0]}" if len(per) == 1 else f"{per[0]}–{per[-1]}")) + " | "
                   + ', '.join(f"{max(r['maxdev'][q] for r in rs):.3f}" for q in range(nq)) + " |")
    out.append("")
    rise = {}                                                          # referee round 69: how far single tones drift without jumping
    for key, rs in c2_groups(C2).items():
        if key[0] == 'one tone' and key[1] == 'paradigm':
            for r in rs:
                if 'w_end' in r and not (r['left'][0] >= 0 and r['maxdev'][0] > JUMP_C):
                    rise.setdefault(key[5], []).append(float(r['w_end'][0] - CT.ASTAR[0]))
    if rise:
        out += ["Single tones (the paradigm's metric), the runs whose first kernel did not jump: its coefficient's rise over the fit by the "
                f"last period, per level — " + '; '.join(f"s = {s_:g}: {min(v):+.3f} … {max(v):+.3f} ({len(v)} runs)"
                                                          for s_, v in sorted(rise.items())) + '.', '']
    return out


def c3_lines(C3):
    """C3's report: per run the leaves, the coefficient's drift to the leave, the mean offset before it and where the kernel ends."""
    if C3 is None:
        return []
    out = ["## C3. The anatomy of a leave (referee round 63)", "",
           f"C's loop (m = 12, k_e = 10 unless stated), the formed kernels' period means over {NP_C} periods: per kernel that left, the "
           "leave period, the coefficient there against its start, the mean signed offset of the centre over the periods leave − 30 … "
           f"leave − 6, its largest distance (> {JUMP_C:g}: a jump) and the centre at the end:", "",
           "| noise | s | seed | kernel | leave period | w: start → leave | mean offset before | largest distance | centre at the end |",
           "|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(C3, key=lambda r: (c2_noise_label({k: v for k, v in r['opts'].items() if k != 'seed'}), r['s'], r['opts'].get('seed', 7))):
        nz = c2_noise_label({k: v for k, v in r['opts'].items() if k != 'seed'})
        W, C = np.array(r['W']), np.array(r['C'])
        for q, p in enumerate(r['left']):
            if p < 0:
                continue
            out.append(f"| {nz} | {r['s']:g} | {r['opts'].get('seed', '—')} | {q + 1} | {p} | {W[0, q]:.3f} → {r['w_left'][q]:.3f} | "
                       f"{r['pre'][q]:+.3f} | {r['maxdev'][q]:.2f} | {C[-1, q]:+.2f} |")
    out.append("")
    grp = {}                                                           # the leaves summed up: the band's runs per level, the single tones
    for r in C3:
        tone = 'freqs' in r['opts'] and len(r['opts']['freqs']) == 1
        if 'freqs' in r['opts'] and not tone:
            continue
        W, dev = np.array(r['W']), np.abs(np.array(r['C']) - CT.CSTAR)
        for q, p in enumerate(r['left']):
            if p >= 0:
                key = 'single tones' if tone else f"the 16 tones at s = {r['s']:g}"
                pj = int(np.argmax(dev[:, q] > JUMP_C)) + 1 if (dev[:, q] > JUMP_C).any() else None
                grp.setdefault(key, []).append((p, r['w_left'][q] / W[0, q], r['pre'][q], r['maxdev'][q], q,
                                                r['opts'].get('seed'), None if pj is None else pj - p, float(np.array(r['C'])[-1, q])))
    span = lambda v, f_: f"{min(v):{f_}}" if min(v) == max(v) else f"{min(v):{f_}} … {max(v):{f_}}"
    parts = []
    for key in sorted(grp, key=lambda k: (k == 'single tones', k)):
        xs = grp[key]
        k1 = [x for x in xs if x[4] == 0]
        runs = {x[5] for x in xs}
        parts.append(f"{key}: {len(xs)} leaves ({len(k1)} of the first kernel" + (f", in {len(runs)} runs" if key != 'single tones' else '')
                     + f"), periods {span([x[0] for x in xs], 'd')}, the coefficient at the leave over its start "
                     f"{span([x[1] for x in xs], '.3f')}, the mean offset before {span([x[2] for x in xs], '+.3f')}, the largest distance "
                     f"{span([x[3] for x in xs], '.2f')} ({sum(x[3] > JUMP_C for x in xs)} jumps"
                     + (f", beyond {JUMP_C:g} {span([x[6] for x in xs if x[6] is not None], 'd')} periods after the leave"
                        if any(x[6] is not None for x in xs) else '') + ")")
    if parts:
        out += ["Summed up — " + '; '.join(parts) + ".", ""]
    per = []                                                           # referee round 69: the same per kernel, and where the jumps end
    for key in sorted(grp, key=lambda k: (k == 'single tones', k)):
        for q in sorted({x[4] for x in grp[key]}):
            xs = [x for x in grp[key] if x[4] == q]
            jm = [x for x in xs if x[3] > JUMP_C]
            per.append(f"{key}, kernel {q + 1}: {len(xs)} leave{'s' if len(xs) != 1 else ''}, periods {span([x[0] for x in xs], 'd')}, the coefficient over its start "
                       f"{span([x[1] for x in xs], '.3f')}, the mean offset before {span([x[2] for x in xs], '+.3f')}, the largest distance "
                       f"{span([x[3] for x in xs], '.2f')} ({len(jm)} jumps" + (f", {sum(abs(x[7]) > 1 for x in jm)} of them ending beyond ±1"
                                                                               if jm else '') + ")")
    if per:
        out += ["Per kernel — " + '; '.join(per) + ".", ""]
    return out


S_W2, NP_W2, W2_SEEDS = 0.00125, 300, (7, 1, 2, 3)                  # W2 (referee round 69): each draw's own second order, C's horizon
W2_CONFIGS = ([(12.0, ke, 'paradigm', None) for ke in KE_C] + [(3.0, ke, 'paradigm', None) for ke in KE_C]
              + [(m, 10.0, 'paradigm', None) for m in M_SCAN_C2] + [(12.0, ke, 'wfr', None) for ke in (5.0, 10.0)]
              + [(m, 10.0, 'paradigm', (25.0, 75.0)) for m in (12.0, 3.0)])
EDGE_GRID, NP_EDGE, BISECT_EDGE = (0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.16, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0), 100, 4
EDGE_CASES = [(m, ke, mt, 0) for m, ke, mt, band in W2_CONFIGS if band is None] + [(12.0, 10.0, 'paradigm', 2)]
EDGE_OUT = [(12.0, 2.5, 'paradigm', 0), (12.0, 80.0, 'paradigm', 0), (1.5, 10.0, 'paradigm', 0), (24.0, 10.0, 'paradigm', 0),
            (12.0, 20.0, 'wfr', 0), (12.0, 40.0, 'wfr', 0), (12.0, 5.0, 'paradigm', 2), (12.0, 20.0, 'paradigm', 2), (3.0, 10.0, 'paradigm', 2)]
OM_EDGE = (0.1, 0.4)                                                    # referee round 76: the edge outside its fitted window (EDGE_OUT)
EDGE_CASES += EDGE_OUT + [(m, ke, 'paradigm', 0, om) for m, ke in ((12.0, 10.0), (12.0, 20.0), (3.0, 10.0)) for om in OM_EDGE]
TRACE_K = 1024                                                          # the ejection's anatomy: outputs a period (referee round 76)
C6_CASES = ([(12.0, 10.0, 0.02, dict(seed=sd, rise=[0, d])) for sd in (7, 3, 10) for d in (0.03, 0.05, 0.065)]
            + [(12.0, 10.0, 0.02, dict(freqs=(f,), seed=7, rise=[0, d])) for f in (13.1, 15.0) for d in (0.1, 0.2, 0.3)])
S_C7, NP_C7, BISECT_C7 = (0.0025, 0.005, 0.01, 0.02, 0.03), 10, 6         # C7 (referee round 69): the edge under noise
S_C7_HI = (0.04, 0.06, 0.08)                                            # referee round 76: the tones also at C2's levels, and under −n
SCAN_C7 = tuple(float(x) for x in np.round(np.r_[np.arange(0, 0.2, 0.0025), np.arange(0.2, 0.4, 0.005), np.arange(0.4, 1.0 + 1e-9, 0.02)],
                                           6))                          # referee round 76: the raises scanned (the ejecting set)
SHIFT_C7, TONES_C7 = (0, 25, 50, 75, 100, 125, 150, 175), (10.0, 13.1, 15.0, 21.3)
C7_CFG = (12.0, 10.0, 'paradigm')                                       # C7's configuration (m, k_e, metric); then those where the drift's
C7_MORE = (((12.0, 20.0, 'paradigm'), (0.02, 0.04, 0.06, 0.08)),         # reach of the noise-free edge misses the first jump (EDGE)
           ((4.0, 10.0, 'paradigm'), (0.04, 0.06, 0.08)), ((12.0, 10.0, 'wfr'), (0.08, 0.16, 0.24)))
C7_CASES = ([(dict(seed=7), 0.0)] + [(dict(freqs=(f,), seed=7, **ng), s_) for f in TONES_C7 for s_ in S_C7 + S_C7_HI
                                        for ng in ({}, dict(neg=True))]
            + [(dict(seed=7, shift=sh), s_) for sh in SHIFT_C7 for s_ in S_C7]
            + [(dict(seed=7), 0.0, cfg) for cfg, _ in C7_MORE] + [(dict(seed=7, shift=sh), s_, cfg) for cfg, lv in C7_MORE
                                                               for sh in SHIFT_C7 for s_ in lv])
C9_CASES = ([(dict(seed=3, shift=sh), 0.02) for sh in (140, 145, 150, 155, 160, 165, 185, 190)]
            + [(dict(seed=10, shift=sh), 0.02) for sh in (100, 105, 110)]
            + [(dict(seed=7, shift=sh), 0.02) for sh in (55, 110, 115, 125, 130)])   # C9 (referee round 76): C6's draws' windows
STEP_C8 = 10                                                            # C8 (referee round 76): a window every STEP_C8 periods
C8_MORE = (((12.0, 20.0, 'paradigm'), (0.04, 0.06, 0.08)), ((4.0, 10.0, 'paradigm'), (0.04, 0.06, 0.08)),
           ((12.0, 10.0, 'wfr'), (0.16, 0.24, 0.32)))
C8_CASES = ([(12.0, 10.0, 'paradigm', s_, dict(seed=sd)) for s_ in (0.02, 0.03) for sd in range(1, 12)]
            + [(12.0, 10.0, 'paradigm', s_, dict(seed=sd)) for s_ in (0.01, 0.015) for sd in W2_SEEDS]
            + [(cf[0], cf[1], cf[2], s_, dict(seed=sd)) for cf, lv in C8_MORE for s_ in lv for sd in W2_SEEDS])


def raised(L, al0, idx, k, delta):
    """al0 with formed kernel k's coefficient raised by delta in magnitude on its curve from the masses EPS (referee round 69)."""
    i = idx[k]
    w0 = mobilities(L, EPS, al0, 'paradigm')[0][i]
    al = np.array(al0, float)
    al[i] = alpha_of(L, EPS, w0 + np.sign(w0) * delta, al0[i])
    return al


def edge_run(m, ke, metric, k, delta, NP=NP_EDGE, trace=False, noise=None, om=None):
    """The noise-free loop (or under the noise `noise`, C7) from the fitted state with formed kernel k's coefficient raised by delta (the
    plant unchanged; referee round 69: the fit's basin in the coefficient), NP periods (64 outputs a period) of the reference at om
    (default ω_ref; referee round 76): the period kernel k leaves its feature (LEFT_C, the period means; −1: never) and its largest
    period-averaged distance; with trace (TRACE_K outputs a period; referee round 76) also the ejection's anatomy — the first time (in
    periods) from which its centre stays at least 0.05 from its feature to the end (non-returning; None if it is back by the end), the
    reference x_m there, the centre at the end of that period and at the end, and the centre's least and largest offset from its feature."""
    L = law(m)
    om = CT.OM if om is None else om
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    kk = TRACE_K if trace else 64
    T = 2 * np.pi / om
    ts = np.linspace(0, NP * T, NP * kk + 1)
    r = simulate(L, raised(L, al0, idx, k, delta), c0, CT.ASTAR, ts, noise, metric, ke, om=om)
    dev = np.abs(r['c'][:-1].reshape(NP, kk, c0.size).mean(1)[:, idx] - CT.CSTAR)
    out = (left_periods(dev, LEFT_C)[k], float(dev[:, k].max()))
    if not trace:
        return out
    off = r['c'][:, idx[k]] - CT.CSTAR[k]
    near = np.flatnonzero(np.abs(off) <= 0.05)
    ext = (float(r['c'][-1, idx[k]]), float(off.min()), float(off.max()))
    if near.size and near[-1] == len(off) - 1:
        return out + (None, None, None) + ext
    j = int(near[-1]) + 1 if near.size else 0
    p = int(ts[j] // T)
    return out + (float(ts[j] / T), float(np.sin(om * ts[j])), float(r['c'][min((p + 1) * kk, len(ts) - 1), idx[k]])) + ext


def bracket_edge(run, grid=None, nb=None, above=0):
    """The edge of run(Δ) — (the leave period, −1: none; the largest distance) — up the grid (default EDGE_GRID) to the first Δ that
    ejects, then nb bisections (default BISECT_EDGE), and the next `above` grid raises beyond the first ejecting one rechecked: (lo, hi,
    the runs sorted by Δ); hi None if no raise of the grid ejects."""
    grid, nb = (EDGE_GRID if grid is None else grid), (BISECT_EDGE if nb is None else nb)
    runs = []
    for d in grid:
        runs.append((d,) + tuple(run(d)))
        if runs[-1][1] >= 0:
            break
    hi = next((d for d, p, _ in runs if p >= 0), None)
    lo = max([d for d, p, _ in runs if p < 0 and (hi is None or d < hi)], default=0.0)
    if hi is not None:
        for d in [g for g in grid if g > hi][:above]:
            runs.append((d,) + tuple(run(d)))
        for _ in range(nb):
            mid = 0.5 * (lo + hi)
            runs.append((mid,) + tuple(run(mid)))
            lo, hi = (lo, mid) if runs[-1][1] >= 0 else (mid, hi)
    return lo, hi, sorted(runs)


def edge_job(args):
    """EDGE: the basin's edge in kernel k's coefficient for the law's m, the gain k_e, the metric and (referee round 76) the reference's
    frequency om (default ω_ref) — `edge_run` up EDGE_GRID to the first raise that ejects the kernel, then BISECT_EDGE bisections
    (`bracket_edge`): the edge's bracket (lo, hi], every run, and the anatomy (`edge_run` with trace, 3 periods) at hi and at lo."""
    m, ke, metric, k = args[:4]
    om = args[4] if len(args) > 4 else CT.OM
    warnings.filterwarnings('ignore')
    t0 = time.time()
    lo, hi, runs = bracket_edge(lambda d: edge_run(m, ke, metric, k, d, om=om))
    out = dict(m=m, ke=ke, metric=metric, k=k, om=om, lo=lo, hi=hi, runs=[dict(d=d, left=p, maxdev=x) for d, p, x in runs],
               secs=time.time() - t0)
    keys = ('t_detach', 'x_detach', 'c_after', 'c_end', 'off_min', 'off_max')
    for tag, d in (('hi', hi), ('lo', lo)):
        if d is not None:
            out[tag + '_trace'] = dict(zip(keys, edge_run(m, ke, metric, k, d, NP=3, trace=True, om=om)[2:]))
    if hi is not None:                                                  # the round-69 fields, now from the fine grid
        out.update({x: out['hi_trace'][x] for x in ('t_detach', 'x_detach', 'c_after')})
    print(f"  EDGE m {m:g} k_e {ke:g} {metric} kernel {k + 1} om {om:g}: ({lo:.4g}, {hi}] [{out['secs']:.0f}s]", flush=True)
    return out


def w2_job(args):
    """W2 (referee round 69): one draw's second order at C's horizon — the loop under ±S_W2·n (the law's m, the gain k_e, the metric, the
    noise `c_noise(S_W2, opts)`) for NP_W2 periods (`sim_means`): per period the even part of the formed kernels' period-mean
    coefficients per s², ((W₊ + W₋)/2 − A)/s² (B4's construction)."""
    m, ke, metric, opts = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law(m)
    nz = c_noise(S_W2, dict(opts))
    Wp, _ = sim_means(L, metric, nz, ke, NP_W2)
    Wm, _ = sim_means(L, metric, dict(nz, ph=np.asarray(nz['ph'], float) + np.pi), ke, NP_W2)
    ew = ((Wp + Wm) / 2 - CT.ASTAR) / S_W2 ** 2
    out = dict(m=m, ke=ke, metric=metric, opts={k_: (list(v) if isinstance(v, tuple) else v) for k_, v in dict(opts).items()},
               even_w_all=ew.tolist(), secs=time.time() - t0)
    print(f"  W2 m {m:g} k_e {ke:g} {metric} {out['opts']}: even/s² at {NP_W2} {', '.join(f'{x:+.4g}' for x in ew[-1])} "
          f"[{out['secs']:.0f}s]", flush=True)
    return out


def w2_cases():
    """W2's draws: C3's eleven band draws at m = 12, k_e = 10 and, for every other configuration of C and C2 with the band's tones,
    the draws W2_SEEDS."""
    cases = [(12.0, 10.0, 'paradigm', dict(seed=sd)) for sd in range(1, 12)]
    for m, ke, mt, band in W2_CONFIGS:
        if (m, ke, mt, band) == (12.0, 10.0, 'paradigm', None):
            continue
        cases += [(m, ke, mt, dict(seed=sd) if band is None else dict(band=band, seed=sd)) for sd in W2_SEEDS]
    return cases


def part_W2():
    return _pool_map(w2_job, w2_cases())


def part_EDGE():
    return _pool_map(edge_job, EDGE_CASES)


def part_C4():
    """C4 (referee round 69): the −n twins of C3's band runs — the noise negated (its phases moved by π)."""
    return _pool_map(c3_job, [(12.0, 10.0, s_, dict(seed=sd, neg=True)) for s_ in (0.02, 0.03) for sd in range(1, 12)])


def part_C6():
    """C6 (referee round 69): C3's runs started with kernel 1's coefficient raised (opts rise) — bands below the noise-free edge, single
    tones beyond it."""
    return _pool_map(c3_job, C6_CASES)


def c7_job(args):
    """C7 (referee rounds 69, 76): the basin in kernel 1's coefficient under the noise `c_noise(s, opts)` at the configuration (m, k_e,
    metric) — args (opts, s) or (opts, s, (m, k_e, metric)), default C7_CFG —: `edge_run` over NP_C7 periods at every raise of SCAN_C7
    (referee round 76: under noise the ejecting set need not be an interval), each raise's leave period and largest distance."""
    opts, s = args[:2]
    m, ke, mt = args[2] if len(args) > 2 else C7_CFG
    warnings.filterwarnings('ignore')
    t0 = time.time()
    nz = c_noise(s, dict(opts))
    runs = [(d,) + tuple(edge_run(m, ke, mt, 0, d, NP=NP_C7, noise=nz)) for d in SCAN_C7]
    out = dict(m=m, ke=ke, metric=mt, s=s, opts={k_: (list(v) if isinstance(v, tuple) else v) for k_, v in dict(opts).items()},
               runs=[dict(d=d, left=p, maxdev=x) for d, p, x in runs], secs=time.time() - t0)
    lv = [d for d, p, _ in runs if p >= 0]
    print(f"  C7 m {m:g} k_e {ke:g} {mt} {out['opts']} s {s:g}: least leaving raise {min(lv) if lv else None}, "
          f"{len(lv)} of {len(runs)} [{out['secs']:.0f}s]", flush=True)
    return out


def part_C7():
    """C7 (referee round 69): `c7_job` for C7_CASES."""
    return _pool_map(c7_job, C7_CASES)


def c8_job(args):
    """C8 (referee round 76): the edge under noise tested in time, run by run — C's run at (m, k_e, metric, s, opts) for NP_C periods
    (`c_run`), and every STEP_C8 periods p a window: the noise-free loop's fitted state with kernel 1's coefficient raised by the run's own
    offset there, δ(p) = |w₁| over the period before p minus its fitted value (0 at p = 0), under the run's own noise from period p
    (`c_noise` with shift p) for NP_C7 periods (`edge_run`). Returns the run's leave periods and largest distances, kernel 1's period-mean
    coefficients, and per window (p, δ(p), its leave period, its largest distance)."""
    m, ke, mt, s, opts = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    r, L, al0, c0, idx, Wp, Cp = c_run(m, ke, s, dict(opts, metric=mt))
    dev = np.abs(Cp - CT.CSTAR)
    w0 = abs(float(mobilities(L, EPS, al0, 'paradigm')[0][idx[0]]))
    wins = []
    for p in range(0, NP_C, STEP_C8):
        dl = 0.0 if p == 0 else abs(float(Wp[p - 1, 0])) - w0
        lp, md = edge_run(m, ke, mt, 0, dl, NP=NP_C7, noise=c_noise(s, dict(opts, shift=p)))
        wins.append(dict(p=p, d=dl, left=lp, maxdev=md))
    out = dict(m=m, ke=ke, metric=mt, s=s, opts={k_: (list(v) if isinstance(v, tuple) else v) for k_, v in dict(opts).items()},
               left=left_periods(dev, LEFT_C), maxdev=dev.max(0).tolist(), w=Wp[:, 0].tolist(), w0=w0, wins=wins, secs=time.time() - t0)
    fl = [x['p'] for x in wins if x['left'] >= 0]
    print(f"  C8 m {m:g} k_e {ke:g} {mt} s {s:g} {out['opts']}: left {out['left'][0]}, flagged {fl[:4]} [{out['secs']:.0f}s]", flush=True)
    return out


def part_C9():
    """C9 (referee round 76): `c7_job` (the scan of SCAN_C7) for the windows of C6's draws around their leaves, C9_CASES."""
    return _pool_map(c7_job, C9_CASES)


def c9_pairs(C3, C6, windows=None):
    """C9's (run, window) pairs (referee round 76): for each window (seed, from period p) — default C9_CASES' — and each run of its draw
    at s = 0.02, un-raised (C3) or raised by Δ at the start (C6), that has not left by p: (seed, p, Δ, the run's offset δ(p) — kernel 1's
    |w| over the period before p minus its fitted value —, the run's leave period)."""
    L = law(C7_CFG[0])
    al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
    w0 = abs(float(mobilities(L, EPS, al0, 'paradigm')[0][idx[0]]))
    runs = [r for r in list(C3 if C3 is not None else []) + list(C6 if C6 is not None else [])
            if r['s'] == 0.02 and 'freqs' not in r['opts'] and 'seed' in r['opts'] and not r['opts'].get('neg')]
    wins = windows if windows is not None else [(o['seed'], int(o['shift'])) for o, _ in C9_CASES]
    out = []
    for sd, p in wins:
        for q in sorted([q for q in runs if q['opts'].get('seed') == sd], key=lambda q: q['opts'].get('rise', [0, 0.0])[1]):
            P = int(q['left'][0])
            if p < 1 or 0 <= P <= p:
                continue
            out.append((sd, p, float(q['opts'].get('rise', [0, 0.0])[1]), abs(float(np.asarray(q['W'], float)[p - 1, 0])) - w0, P))
    return out


def c9_own_job(args):
    """C9D (referee round 76): the window's loop at a run's own offset — `edge_run` at C7_CFG with kernel 1 raised by δ, under the draw's
    noise at s = 0.02 from period p, NP_C7 periods; args (seed, p, Δ, δ). Its leave period counts from the window's start (−1: none)."""
    sd, p, rise, d = args
    m, ke, mt = C7_CFG
    warnings.filterwarnings('ignore')
    lp, md = edge_run(m, ke, mt, 0, d, NP=NP_C7, noise=c_noise(0.02, dict(seed=sd, shift=p)))
    return dict(seed=sd, p=p, rise=rise, d=d, left=lp, maxdev=md)


def part_C9D(C3=None, C6=None):
    """C9D (referee round 76): `c9_own_job` for every pair of `c9_pairs` (C3's and C6's runs)."""
    return _pool_map(c9_own_job, [x[:4] for x in c9_pairs(C3, C6)])


def part_C8():
    """C8 (referee round 76): `c8_job` for C8_CASES."""
    return _pool_map(c8_job, C8_CASES)


def w2_of(W2, m, ke, metric, band, seed):
    """W2's even part per s² (per period, per kernel) for a configuration and draw, or None."""
    for r in (W2 if W2 is not None else []):
        o = r['opts']
        if (r['m'], r['ke'], r['metric'], tuple(o['band']) if 'band' in o else None, o.get('seed')) == (m, ke, metric, band, seed):
            return np.asarray(r['even_w_all'], float)
    return None


def c4_lines(C3, C4, W2):
    """C4 (referee round 69): like for like — each C3 band run's rise of kernel 1 at its leave against its own draw's even part at the same
    period and level: from its −n twin, (W₊ + W₋)/2 − A (while the twin has not left), and from W2's second order, s²⟨w₂⟩."""
    if C3 is None or C4 is None or len(C3) == 0 or len(C4) == 0:
        return []
    tw = {(r['s'], r['opts'].get('seed')): r for r in C4}
    out = ['', "### C4. Like for like (referee round 69): the rise at the leave against the same draw's even part", '',
           f"C3's band runs (m = 12, k_e = 10) that lost kernel 1: its leave period p, the rise W₊ − A of its period-mean coefficient at "
           f"p, the −n twin's leave, the half-sum (W₊ + W₋)/2 − A at p (where the twin has not left by p) over the rise, the draw's own second "
           f"order s²⟨w₂⟩ at p (W2: ±{S_W2:g}·n over {NP_W2} periods) and its share of the value at period {NP_W2}, the rise over s²⟨w₂⟩(p):", '',
           f"| s | draw | p | rise | the twin's leave | half-sum / rise | s²⟨w₂⟩(p) (share of period {NP_W2}'s) | rise / s²⟨w₂⟩(p) |",
           '|---|---|---|---|---|---|---|---|']
    stats = {}
    for r in sorted([r for r in C3 if 'seed' in r['opts'] and 'freqs' not in r['opts']], key=lambda r: (r['s'], r['opts']['seed'])):
        p, sd, s_ = r['left'][0], r['opts']['seed'], r['s']
        if p < 0:
            continue
        W = np.asarray(r['W'], float)
        rise = float(W[p - 1, 0] - CT.ASTAR[0])
        t = tw.get((s_, sd))
        half = None
        if t is not None and (t['left'][0] < 0 or t['left'][0] > p):
            half = float((W[p - 1, 0] + np.asarray(t['W'], float)[p - 1, 0]) / 2 - CT.ASTAR[0])
        ew = w2_of(W2, 12.0, 10.0, 'paradigm', None, sd)
        own = None if ew is None else (s_ ** 2 * float(ew[p - 1, 0]), s_ ** 2 * float(ew[-1, 0]))
        st = stats.setdefault(s_, dict(hr=[], ro=[], sh=[], pairs=[]))
        if half is not None:
            st['hr'].append(half / rise)
        if own is not None:
            st['ro'].append(rise / own[0]); st['sh'].append(own[0] / own[1]); st['pairs'].append((rise, own[0]))
        out.append(f"| {s_:g} | {sd} | {p} | {rise:+.4f} | {t['left'][0] if t is not None else '—'} | "
                   + ('—' if half is None else f"{half:+.4f} ({half / rise:.3f})") + ' | '
                   + ('—' if own is None else f"{own[0]:+.4f} ({own[0] / own[1]:.2f})") + ' | '
                   + ('—' if own is None else f"{rise / own[0]:.3f}") + ' |')
    sm = []
    for s_, st in sorted(stats.items()):
        cc = f"{float(np.corrcoef(np.array(st['pairs']).T)[0, 1]):.2f}" if len(st['pairs']) > 2 else '—'
        sm.append(f"at s = {s_:g} the half-sum over the rise {span_(st['hr'], '.3f')} ({len(st['hr'])} draws whose twin had not left), "
                  f"the rise over s²⟨w₂⟩(p) {span_(st['ro'], '.3f')}" + (f" (median {np.median(st['ro']):.3f})" if st['ro'] else '')
                  + f", s²⟨w₂⟩(p) over its period-300 value "
                  f"{span_(st['sh'], '.2f')}, the rise's correlation with s²⟨w₂⟩(p) across the {len(st['pairs'])} draws {cc}")
    return out + ['', "Over the draws: " + '; '.join(sm) + '.', '']


def span_(v, f_):
    return f"{min(v):{f_}}–{max(v):{f_}}" if len(v) else '—'


def first_jump(rs, k):
    """The smallest level among the runs rs at which kernel k jumped (left and its largest distance above JUMP_C), None if none."""
    js = [r['s'] for r in rs if r['left'][k] >= 0 and r['maxdev'][k] > JUMP_C]
    return min(js) if js else None


_V1 = {}


def v_of(m, metric, k=0):
    """Formed kernel k's centre mobility at the fit (`mobilities`: Gmw_k under the paradigm's metric, G·a_k under WFR's)."""
    if (m, metric, k) not in _V1:
        L = law(m)
        al0, c0, idx = DZ.fitted_state(L, EPS, CT.CSTAR, CT.ASTAR)
        _V1[(m, metric, k)] = float(mobilities(L, EPS, al0, metric)[2][idx[k]])
    return _V1[(m, metric, k)]


def v1_of(m, metric):
    """The first formed kernel's centre mobility at the fit (`v_of` at k = 0)."""
    return v_of(m, metric, 0)


def star_verdict(a, b, obs, levs):
    """The first jump obs (None: none at the levels levs run) against s*'s bracket (a, b]: 'consistent' where (the level run just below
    obs, obs] meets it, 'later' where that lower level is already at least b, 'earlier' where obs ≤ a; without a jump 'none, predicted'
    where a level above b was run, else 'untested'."""
    if obs is not None:
        lo = max([x for x in levs if x < obs], default=0.0)
        return 'consistent' if (lo < b and obs > a) else ('later' if lo >= b else 'earlier')
    return 'none, predicted' if (levs and max(levs) > b) else 'untested'


def star_cases(EDGE, W2, C, C2):
    """The s* cases (referee rounds 69, 76): per noise-free edge at ω_ref (each EDGE row, and kernel 1's for W2's bands) and per draw of
    W2_SEEDS with a second order and runs: the configuration key, the edge's bracket, the draw's ⟨w₂⟩ at period NP_W2 for the kernel, the
    levels run and the first jump."""
    allruns = (list(C) if C is not None else []) + (list(C2) if C2 is not None else [])
    bkey = lambda o: tuple(float(x) for x in o['band']) if 'band' in o else None
    rows_ = [(e, None) for e in EDGE if e.get('om', CT.OM) == CT.OM] + [
        (e, b) for e in EDGE for (m_, ke_, mt_, b) in W2_CONFIGS
        if b is not None and e['k'] == 0 and e.get('om', CT.OM) == CT.OM and (m_, ke_, mt_) == (e['m'], e['ke'], e['metric'])]
    out = []
    for e, band in rows_:
        m, ke, mt, k = e['m'], e['ke'], e['metric'], e['k']
        for sd in W2_SEEDS:
            ew = w2_of(W2, m, ke, mt, band, sd)
            rs = [r for r in allruns if r['m'] == m and r['ke'] == ke and r['opts'].get('metric', 'paradigm') == mt
                  and 'freqs' not in r['opts'] and bkey(r['opts']) == (None if band is None else tuple(map(float, band)))
                  and r['opts'].get('seed', 7) == sd]
            out.append(dict(key=(m, ke, mt, k, None if band is None else tuple(band)), e=e, seed=sd,
                            w2=None if ew is None else float(ew[-1, k]), levs=sorted({r['s'] for r in rs}),
                            obs=first_jump(rs, k) if rs else None, run=bool(rs)))
    return out


def star_count(cases, w2=None):
    """The s* verdicts (`star_verdict`) of `star_cases`, with the draws' ⟨w₂⟩ replaced by w2 (a list in the cases' order) if given:
    a list of (case, verdict, s*'s bracket) for the cases with an edge, a positive ⟨w₂⟩ and runs."""
    out = []
    for i, q in enumerate(cases):
        x = q['w2'] if w2 is None else w2[i]
        if q['e']['hi'] is None or x is None or x <= 0 or not q['levs']:
            continue
        a, b = np.sqrt(q['e']['lo'] / x), np.sqrt(q['e']['hi'] / x)
        out.append((q, star_verdict(a, b, q['obs'], q['levs']), (a, b)))
    return out


def star_null(cases, n=200, seed=0):
    """The permutation null of the s* count (referee round 76): the draws' ⟨w₂⟩ permuted within each configuration n times
    (np.random.default_rng(seed)); the 'consistent' count per permutation."""
    rng = np.random.default_rng(seed)
    groups = {}
    for i, q in enumerate(cases):
        groups.setdefault(q['key'], []).append(i)
    base = [q['w2'] for q in cases]
    res = []
    for _ in range(n):
        w2 = list(base)
        for ix in groups.values():
            vals = [base[i] for i in ix]
            for i, j in zip(ix, rng.permutation(len(ix))):
                w2[i] = vals[j]
        res.append(sum(v == 'consistent' for _, v, _ in star_count(cases, w2)))
    return res


def edge_fit(EDGE):
    """The edge's power law (referee rounds 69, 76): over the paradigm's kernel-1 rows at ω_ref with m in [3, 12] and k_e in [5, 40]
    (the window it is fitted in), log Δ_e = c + α log k_e − β log m by least squares (Δ_e the bracket's midpoint): (c, α, β), their
    standard errors, the rows and their residuals; None with fewer than four rows."""
    mid = lambda e: 0.5 * (e['lo'] + e['hi'])
    win = [e for e in EDGE if e['metric'] == 'paradigm' and e['k'] == 0 and e.get('om', CT.OM) == CT.OM and e['hi'] is not None
           and 3.0 <= e['m'] <= 12.0 and 5.0 <= e['ke'] <= 40.0]
    if len(win) < 4 or len({e['ke'] for e in win}) < 2 or len({e['m'] for e in win}) < 2:
        return None
    X = np.array([[1.0, np.log(e['ke']), np.log(e['m'])] for e in win])
    y = np.array([np.log(mid(e)) for e in win])
    cf, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ cf
    s2 = float(res @ res) / max(len(y) - 3, 1)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return dict(c=float(cf[0]), alpha=float(cf[1]), beta=float(-cf[2]), se=(float(se[0]), float(se[1]), float(se[2])), rows=win,
                res=res.tolist(), law=lambda m, ke: float(np.exp(cf[0] + cf[1] * np.log(ke) + cf[2] * np.log(m))))


def edge_lines(EDGE, W2, C, C2):
    """EDGE (referee rounds 69, 76): the basin's edge per configuration (and, outside the fitted window, at more gains, m, the metric and
    the reference's frequency), the ejection's anatomy, the power law with its errors in its window and the deviations outside it, and the
    drift's prediction of the first jump — s* = (edge/⟨w₂⟩(300))^{1/2} — against the observed, with its permutation null."""
    if EDGE is None or len(EDGE) == 0:
        return []
    mid = lambda e: 0.5 * (e['lo'] + e['hi'])
    out = ['', "### EDGE. The fit's basin in the coefficient, and the drift's reach (referee rounds 69, 76)", '',
           f"The noise-free loop from the fitted state with one formed kernel's coefficient raised in magnitude by Δ (the plant unchanged), "
           f"{NP_EDGE} periods of the reference at ω_ref (ω_ref = {CT.OM:g} unless stated): the smallest Δ of the grid "
           f"{', '.join(f'{d:g}' for d in EDGE_GRID)} that ejects the kernel (its period-averaged centre {LEFT_C:g} from its feature) and "
           f"{BISECT_EDGE} bisections — the edge's bracket (lo, hi] and the ejection's period at hi —; the anatomy over 3 periods on "
           f"{TRACE_K} outputs a period — at hi the first time from which the centre stays 0.05 off its feature (t/T), x_m then, the centre at "
           "the end of that period and at the end; at lo the centre's least and largest offset from its feature —; then per draw the level "
           f"s* = (edge/⟨w₂⟩)^{{1/2}} at which the draw's even drift at period {NP_W2} (W2, per s²) reaches the edge (the bracket mapped) "
           "against the smallest level at which C's or C2's runs of that draw saw the kernel jump (— none up to the largest level run):", '',
           '| m | k_e | metric | ω_ref | kernel | noise | edge (lo, hi] | ejected at period (at hi: t/T, x_m; the centre after, at the end) '
           '| at lo: the offsets | ' + ' | '.join(f'draw {sd}: s* (lo, hi], observed' for sd in W2_SEEDS) + ' |',
           '|---|---|---|---|---|---|---|---|---|' + '---|' * len(W2_SEEDS)]
    cases = star_cases(EDGE, W2, C, C2)
    byrow = {}
    for q in cases:
        byrow.setdefault((id(q['e']), q['key'][4]), []).append(q)
    rows_ = [(e, None) for e in EDGE] + [(e, b) for e in EDGE for (m_, ke_, mt_, b) in W2_CONFIGS
                                          if b is not None and e['k'] == 0 and e.get('om', CT.OM) == CT.OM
                                          and (m_, ke_, mt_) == (e['m'], e['ke'], e['metric'])]
    for e, band in rows_:
        m, ke, mt, k, om = e['m'], e['ke'], e['metric'], e['k'], e.get('om', CT.OM)
        p_hi = next((x['left'] for x in e['runs'] if x['d'] == e['hi']), None)
        cells = []
        for q in byrow.get((id(e), None if band is None else tuple(band)), []):
            x = q['w2']
            if x is None or e['hi'] is None or x <= 0:
                cells.append(f"—, {q['obs'] if q['obs'] is not None else '—'}")
                continue
            a, b = np.sqrt(e['lo'] / x), np.sqrt(e['hi'] / x)
            cells.append(f"({a:.4f}, {b:.4f}], {q['obs'] if q['obs'] is not None else '—'}")
        cells += ['—'] * (len(W2_SEEDS) - len(cells))
        h, l_ = e.get('hi_trace'), e.get('lo_trace')
        an = ('' if not h or h.get('t_detach') is None else
              f" ({h['t_detach']:.3f}, {h['x_detach']:+.3f}; {h['c_after']:+.3f}, {h['c_end']:+.3f})")
        lo_ = '—' if not l_ else f"{l_['off_min']:+.3f}, {l_['off_max']:+.3f}"
        lab = "D‴.49's band" if band is None else f"{band[0]:g}–{band[1]:g}"
        edge = f"({e['lo']:.4g}, {e['hi']:.4g}]" if e['hi'] is not None else f"> {EDGE_GRID[-1]:g}"
        out.append(f"| {m:g} | {ke:g} | {mt} | {om:g} | {k + 1} | {lab} | {edge} | {p_hi if e['hi'] is not None else '—'}{an} | {lo_} | "
                   + ' | '.join(cells) + ' |')
    recs = star_count(cases)
    cnt = lambda v: sum(1 for _, x, _ in recs if x == v)
    testable = sum(1 for _, x, _ in recs if x != 'untested')

    def where(v):
        g = {}
        for q, x, (a, b) in recs:
            if x == v:
                g.setdefault(q['key'], []).append(None if q['obs'] is None else q['obs'] / (0.5 * (a + b)))
        return ', '.join(f"m = {m:g}, k_e = {ke:g}{'' if mt == 'paradigm' else ', ' + mt}{'' if k == 0 else f', kernel {k + 1}'}"
                         f"{'' if band is None else f', {band[0]:g}–{band[1]:g}'} ({len(v_)}"
                         + (f"; {min(v_):.2f}–{max(v_):.2f} times s*" if v_[0] is not None else '') + ')'
                         for (m, ke, mt, k, band), v_ in g.items())
    if recs:
        null = star_null(cases)
        wid = [(q['obs'] - max([x for x in q['levs'] if x < q['obs']], default=0.0)) / (0.5 * (a + b))
               for q, x, (a, b) in recs if x == 'consistent']
        out += ['', f"The first jump against the drift's reach — (the level run just below it, the first jump] meeting s*'s bracket —: "
                    f"consistent in {cnt('consistent')} of the {testable} testable draw-configurations (those that jumped, and those that did "
                    f"not where a level above s* was run), later in {cnt('later')}" + (f" ({where('later')})" if cnt('later') else '')
                    + f", earlier in {cnt('earlier')}" + (f" ({where('earlier')})" if cnt('earlier') else '')
                    + f"; of those that never jumped, s* lay below the largest level run in {cnt('none, predicted')}"
                    + (f" ({where('none, predicted')})" if cnt('none, predicted') else '')
                    + f". With the draws' ⟨w₂⟩ permuted within each configuration ({len(null)} permutations) the consistent count is "
                    f"{np.mean(null):.2f} on average ({min(null)}–{max(null)}; at least {cnt('consistent')} in "
                    f"{sum(x >= cnt('consistent') for x in null)}): the prediction is the configuration's, not the draw's. The consistent "
                    f"cases' windows — the level below the jump to the jump — are {span_(wid, '.2f')} of s* wide.", '']
    fit = edge_fit(EDGE)
    if fit is not None:
        cc = lambda e: mid(e) * np.sqrt(v_of(e['m'], e['metric'], e['k']) / e['ke'])
        k1w = [cc(e) for e in fit['rows']]
        c0 = float(np.median(k1w))
        oth = [e for e in EDGE if e not in fit['rows']]
        fmt_o = lambda e: (f"m = {e['m']:g}, k_e = {e['ke']:g}" + ('' if e['metric'] == 'paradigm' else f", {e['metric']}")
                           + ('' if e.get('om', CT.OM) == CT.OM else f", ω_ref = {e['om']:g}") + ('' if e['k'] == 0 else f", kernel {e['k'] + 1}"))
        dv = []
        for e in oth:
            pred = c0 * np.sqrt(e['ke'] / v_of(e['m'], e['metric'], e['k']))
            if e['hi'] is None:
                dv.append(f"{fmt_o(e)}: > {EDGE_GRID[-1]:g} (the law {pred:.3f})")
            else:
                pw = (f", the power law's {mid(e) / fit['law'](e['m'], e['ke']):.2f}" if (e['metric'] == 'paradigm' and e['k'] == 0) else '')
                dv.append(f"{fmt_o(e)}: {mid(e):.4f}, {mid(e) / pred:.2f} times the law{pw}")
        oms = {}
        for e in EDGE:
            if e['metric'] == 'paradigm' and e['k'] == 0 and e['hi'] is not None:
                oms.setdefault((e['m'], e['ke']), {})[e.get('om', CT.OM)] = mid(e)
        sl = []
        for (m, ke), d in sorted(oms.items()):
            if len(d) >= 3:
                x, y = np.log(sorted(d)), np.log([d[o] for o in sorted(d)])
                sl.append(f"m = {m:g}, k_e = {ke:g}: {np.polyfit(x, y, 1)[0]:.2f} (" + ', '.join(
                    f"{np.diff(y)[i] / np.diff(x)[i]:.2f}" for i in range(len(x) - 1)) + ')')
        se = fit['se']
        out += [f"The edge's scaling (the bracket's midpoint Δ_e): over the paradigm's {len(fit['rows'])} kernel-1 rows in m = 3–12, k_e = 5–40 at "
                f"ω_ref = {CT.OM:g} (one target), the least-squares Δ_e ≈ c·k_e^α·m^(−β) has α = {fit['alpha']:.3f} ± {se[1]:.3f} and β = "
                f"{fit['beta']:.3f} ± {se[2]:.3f} (standard errors), its residuals {min(fit['res']):+.3f} … {max(fit['res']):+.3f} in the log "
                f"(the corners {', '.join(f'{r:+.2f}' for e, r in zip(fit['rows'], fit['res']) if e['m'] in (3.0, 12.0) and e['ke'] in (5.0, 40.0))}); "
                f"in its (k_e/v₁)^{{1/2}} form, v₁ the kernel's centre mobility at the fit (Gmw₁ under the paradigm's metric — "
                f"{v1_of(12.0, 'paradigm'):.4g} at m = 12 —, G·a₁ under WFR's — {v1_of(12.0, 'wfr'):.4g}), Δ_e·(v₁/k_e)^{{1/2}} = "
                f"{min(k1w):.4f} … {max(k1w):.4f} over the window (median {c0:.4f}). Outside the window, each edge against that law "
                f"{c0:.4f}·(k_e/v)^{{1/2}} (v the kernel's own centre mobility): " + '; '.join(dv) + '.'
                + (" The edge against ω_ref (0.1, 0.2, 0.4), the log-slope by least squares (pairwise): " + '; '.join(sl) + '.' if sl else ''),
                '']
    return out


def c6_lines(C6, C3):
    """C6 (referee round 69): runs started with kernel 1's coefficient raised by Δ — its leave and largest distance against the un-raised
    run's (C3) where there is one."""
    if C6 is None or len(C6) == 0:
        return []
    base = {(r['s'], r['opts'].get('seed')): r for r in (C3 if C3 is not None else []) if 'freqs' not in r['opts']}
    out = ['', "### C6. Started above the fit (referee round 69)", '',
           "C3's runs with kernel 1's coefficient raised by Δ at the start (below the noise-free edge for the band's draws, beyond it for "
           "single tones): kernel 1's leave and largest distance, and the un-raised run's leave:", '',
           '| noise | s | Δ | leave (largest distance) | un-raised: leave |', '|---|---|---|---|---|']
    for r in C6:
        o = r['opts']
        lab = noise_label(o.get('freqs')) if 'freqs' in o else f"16 tones, draw {o.get('seed')}"
        b = base.get((r['s'], o.get('seed'))) if 'freqs' not in o else None
        lv = r['left'][0]
        out.append(f"| {lab} | {r['s']:g} | {o['rise'][1]:g} | {lv if lv >= 0 else 'none'} ({r['maxdev'][0]:.3f}) | "
                   + ('—' if b is None else f"{b['left'][0]}") + ' |')
    return out + ['']


def c7_label(o):
    if 'freqs' in o:
        return noise_label(o['freqs']) + (', −n' if o.get('neg') else '')
    return "D‴.49's 16 tones" + (f", from period {o['shift']:g}" if o.get('shift') else '')


def bands_of(ds, flags):
    """The maximal runs of consecutive raises of the sorted grid ds whose flags are true, as (first, last) pairs (referee round 76)."""
    out, cur = [], None
    for d, f in zip(ds, flags):
        if f:
            cur = (d, d) if cur is None else (cur[0], d)
        elif cur is not None:
            out.append(cur)
            cur = None
    return out + ([cur] if cur is not None else [])


def scan_sets(r):
    """A C7 scan's ejecting sets (referee round 76): the bands of raises whose run left (LEFT_C) and of those that jumped (JUMP_C)."""
    rs = sorted(r['runs'], key=lambda x: x['d'])
    ds = [x['d'] for x in rs]
    return bands_of(ds, [x['left'] >= 0 for x in rs]), bands_of(ds, [x['maxdev'] > JUMP_C for x in rs])


def fmt_bands(b):
    return ', '.join(f"{a:g}" if a == c else f"{a:g}–{c:g}" for a, c in b) if b else 'none'


def c7_lines(C7, W2=None, C=None, C2=None):
    """C7 (referee rounds 69, 76): the basin in kernel 1's coefficient under noise, from the scans of SCAN_C7 — per case the ejecting
    set's bands (left, jumped); for single tones under +n and −n the least ejecting raise, the even part (the half-sum) and its rise over
    the noise-free one per s², and the odd part per s; for the band's windows per configuration and level each window's least ejecting
    raise and band count, their least and median ('>1': none up to the grid's end, counted as above every raise), and the windows whose
    ejecting set is not an interval."""
    if C7 is None or len(C7) == 0:
        return []
    cfg = lambda r: (r.get('m', C7_CFG[0]), r.get('ke', C7_CFG[1]), r.get('metric', C7_CFG[2]))
    least = lambda b: b[0][0] if b else np.inf
    top = max(SCAN_C7)
    f_ = lambda x: f"{x:g}" if np.isfinite(x) else f">{top:g}"
    sets = [(r,) + scan_sets(r) for r in C7]
    e0s = {cfg(r): least(lv) for r, lv, _ in sets if r['s'] == 0}
    out = ['', "### C7. The basin under noise: the ejecting set (referee rounds 69, 76)", '',
           f"As EDGE (kernel 1) over {NP_C7} periods under the noise of rms s, at every raise Δ of the grid 0–0.2 by 0.0025, 0.2–0.4 by "
           f"0.005 and 0.4–{top:g} by 0.02: the raises whose run leaves its feature ({LEFT_C:g}) — the ejecting set, as bands of the grid — "
           f"and those whose run jumps (> {JUMP_C:g}). Without noise (s = 0, per configuration):", '',
           '| m | k_e | metric | leaving raises | jumping raises |', '|---|---|---|---|---|']
    for r, lv, jp in sorted([x for x in sets if x[0]['s'] == 0], key=lambda x: cfg(x[0])):
        m, ke, mt = cfg(r)
        out.append(f"| {m:g} | {ke:g} | {mt} | {fmt_bands(lv)} | {fmt_bands(jp)} |")
    e0 = e0s.get(C7_CFG)
    tones = {}
    for r, lv, jp in sets:
        if cfg(r) == C7_CFG and 'freqs' in r['opts'] and r['s'] > 0:
            tones.setdefault((tuple(r['opts']['freqs']), r['s']), {})['−' if r['opts'].get('neg') else '+'] = (lv, jp)
    if tones:
        out += ['', f"Single tones at m = {C7_CFG[0]:g}, k_e = {C7_CFG[1]:g} (the seed's phase; −n: the tone negated): the least leaving raise "
                    "and the leaving bands under +n and under −n, the even part — the half-sum of the two least raises — and its rise over "
                    f"the noise-free {f_(e0) if e0 is not None else '—'} per s² (E), and the odd part — the half-difference — per s:", '',
                '| tone | s | +n: least (bands) | −n: least (bands) | even | E | odd / s |', '|---|---|---|---|---|---|---|']
        Es = {}
        for (fr, s_), d in sorted(tones.items()):
            if '+' not in d or '−' not in d:
                continue
            a, b = least(d['+'][0]), least(d['−'][0])
            ev = 0.5 * (a + b)
            E_ = (ev - e0) / s_ ** 2 if (e0 is not None and np.isfinite(ev)) else None
            if E_ is not None:
                Es.setdefault(fr, []).append((s_, E_))
            od = (a - b) / (2 * s_) if np.isfinite(a) and np.isfinite(b) else None
            out.append(f"| {noise_label(fr)} | {s_:g} | {f_(a)} ({fmt_bands(d['+'][0])}) | {f_(b)} ({fmt_bands(d['−'][0])}) | {f_(ev)} | "
                       + ('—' if E_ is None else f"{E_:+.4g}") + ' | ' + ('—' if od is None else f"{od:+.3f}") + ' |')
        if Es:
            out += ['', "E per tone and level (the even part's rise per s²; '—' where a least raise is beyond the grid): " + '; '.join(
                f"{noise_label(fr)}: " + ', '.join(f"{e:+.4g} at {s_:g}" for s_, e in sorted(v)) for fr, v in sorted(Es.items())) + '.']
    band = [x for x in sets if 'freqs' not in x[0]['opts'] and x[0]['s'] > 0]
    out += ['', f"D‴.49's 16 tones (draw 7) started at the periods {', '.join(f'{x:g}' for x in SHIFT_C7)} of the draw (its windows), per "
                "configuration and level: each window's least leaving raise (in brackets the number of leaving bands where more than one), "
                "their least and median, the windows whose leaving set is not an interval, and the least jumping raise's least and median:", '',
                '| m | k_e | metric | s | the windows\' least leaving raises | least, median | not an interval | jumping: least, median |',
                '|---|---|---|---|---|---|---|---|']
    gaps = []
    for key in sorted({(cfg(x[0]), x[0]['s']) for x in band}, key=lambda q: (q[0] != C7_CFG, q[0], q[1])):
        (m, ke, mt), s_ = key
        ws = sorted([x for x in band if (cfg(x[0]), x[0]['s']) == key], key=lambda x: x[0]['opts'].get('shift', 0))
        el = [least(lv) for _, lv, _ in ws]
        ej = [least(jp) for _, _, jp in ws]
        ni = [x for x in ws if len(x[1]) > 1]
        gaps += [(key, x[0]['opts'].get('shift', 0), x[1]) for x in ni]
        out.append(f"| {m:g} | {ke:g} | {mt} | {s_:g} | " + ', '.join(f_(e) + (f" [{len(x[1])}]" if len(x[1]) > 1 else '')
                                                                for e, x in zip(el, ws))
                   + f" | {f_(min(el))}, {f_(float(np.median(el)))} | {len(ni)} of {len(ws)} | {f_(min(ej))}, {f_(float(np.median(ej)))} |")
    nw = len(band)
    nf = [x for x in sets if x[0]['s'] == 0]
    bad = [cfg(r) for r, lv, _ in nf if not (len(lv) == 1 and lv[0][1] == top)]
    out += ['', ("Without noise the leaving set is an interval from the edge to the grid's end at every configuration" if not bad else
                 "Without noise the leaving set is not an interval from the edge to the grid's end at " + '; '.join(
                     f"m = {m:g}, k_e = {ke:g}, {mt}" for m, ke, mt in bad)) + f"; under the band's noise it is not an interval in "
                f"{len(gaps)} of the {nw} windows" + (": " + '; '.join(
                    f"m = {k_[0][0]:g}, k_e = {k_[0][1]:g}{'' if k_[0][2] == 'paradigm' else ', ' + k_[0][2]}, s = {k_[1]:g}, from period "
                    f"{sh:g}: {fmt_bands(lv)}" for k_, sh, lv in gaps[:12]) + (f"; and {len(gaps) - 12} more" if len(gaps) > 12 else '')
                    if gaps else '') + '.', '']
    return out


def c9_lines(C9, C3=None, C6=None, C9D=None):
    """C9 (referee round 76): the windows of C6's draws around their leaves (s = 0.02; m = 12, k_e = 10) — each window's leaving set
    (SCAN_C7's bands) — and each run of the draw (`c9_pairs`): its offset δ(p) at the window's start, whether the grid's raise nearest to
    it leaves, the leave of the window's loop at δ itself (C9D, `c9_own_job`), and the run's own leave in that window ((p, p + NP_C7]);
    the 2×2 counts against the run's leave for both, the pairs where the run left though the window's loop at δ stays, and per run whose
    leave some window contains whether one of those windows' loops at δ leaves."""
    if C9 is None or len(C9) == 0:
        return []
    own = {(int(x['seed']), int(x['p']), float(x['rise'])): x for x in (C9D if C9D is not None else [])}
    out = ['', "### C9. The windows of C6's draws around their leaves (referee round 76)", '',
           f"D‴.49's 16 tones at s = 0.02 (m = {C7_CFG[0]:g}, k_e = {C7_CFG[1]:g}): each window's leaving set (as C7, bands of the grid), and each "
           "run of the draw — un-raised (C3) or raised by Δ at the start (C6) — that has not left by the window's start: its offset δ there, "
           "whether the grid's raise nearest to δ leaves (yes/no), the period the window's loop at δ itself leaves (C9D; '—': not within "
           f"the window's {NP_C7} periods), and the period the run leaves ('—': not within the window):", '',
           "| draw | from period | leaving set | the runs (Δ: δ, nearest raise, the window's loop at δ, the run) |", '|---|---|---|---|']
    agree, miss, held = [], [], {}
    for r in sorted(C9, key=lambda r: (r['opts'].get('seed'), r['opts'].get('shift', 0))):
        sd, p = int(r['opts'].get('seed')), int(r['opts'].get('shift', 0))
        lv, _ = scan_sets(r)
        rs = sorted(r['runs'], key=lambda x: x['d'])
        cells = []
        for _, _, rise, d, P in c9_pairs(C3, C6, [(sd, p)]):
            inset = min(rs, key=lambda x: abs(x['d'] - d))['left'] >= 0
            there = p < P <= p + NP_C7
            o = own.get((sd, p, rise))
            ol = None if o is None else o['left'] >= 0
            agree.append((inset, there, d < 0, ol))
            if ol is False and there:
                miss.append((sd, p, rise, P))
            if there and ol is not None:
                held[(sd, rise)] = held.get((sd, rise), False) or ol
            cells.append(f"{rise:g}: {d:+.4f}, {'yes' if inset else 'no'}, "
                         + ('' if o is None else (f"{p + int(o['left'])}" if ol else '—') + ', ') + (f"{P}" if there else '—'))
        out.append(f"| {sd} | {p} | {fmt_bands(lv)} | " + ('; '.join(cells) if cells else '—') + ' |')
    cnt = lambda i, a, b: sum(1 for x in agree if (x[i], x[1]) == (a, b))
    out += ['', f"Of the {len(agree)} (run, window) pairs, by the grid's nearest raise: {cnt(0, True, True)} in the leaving set with the run "
                f"leaving there, {cnt(0, True, False)} in the set but staying, {cnt(0, False, True)} leaving though not in the set, "
                f"{cnt(0, False, False)} neither"
                + (f" ({sum(1 for x in agree if x[2])} with δ < 0, compared with the raise 0)" if any(x[2] for x in agree) else '') + '.']
    if all(x[3] is not None for x in agree):
        out += ['', f"By the window's loop at δ itself: {cnt(3, True, True)} leaving with the run, {cnt(3, True, False)} leaving while the run "
                    f"stays, {cnt(3, False, True)} staying while the run leaves, {cnt(3, False, False)} neither"
                    + ('; the run leaving though the loop at δ stays: ' + '; '.join(
                        f"draw {sd} raised by {rise:g} in the window from period {p} (it leaves at {P})" for sd, p, rise, P in miss)
                       if miss else '') + f". Of the {len(held)} runs whose leave a window contains, {sum(held.values())} are flagged — the "
                    "loop at δ leaving — in at least one of those windows" + (": not " + '; '.join(
                        f"draw {sd} raised by {rise:g}" for (sd, rise), v in sorted(held.items()) if not v) if not all(held.values()) else '')
                    + '.']
    return out + ['']


def c8_verdict(P, flags, step=None):
    """A run's first flag against its first leave P (−1: none; referee round 76): ('hit', p) in the leave's window, ('early', p) before
    it, ('miss', None) none up to it; for a run that never left ('false alarm', p) or ('correct', None)."""
    step = STEP_C8 if step is None else step
    if P < 0:
        return ('false alarm', flags[0]) if flags else ('correct', None)
    pl = step * ((P - 1) // step)
    f = next((x for x in flags if x <= pl), None)
    return ('miss', None) if f is None else (('hit', f) if f == pl else ('early', f))


def c8_confusion(P, flags, starts, step=None):
    """The windows' confusion counts (tp, fn, fp, tn) of one run up to its leave's window (all windows if it never left): positive the
    leave's window, negative every earlier one (referee round 76)."""
    step = STEP_C8 if step is None else step
    pl = None if P < 0 else step * ((P - 1) // step)
    c = np.zeros(4, int)
    for p in starts:
        if pl is not None and p > pl:
            break
        pos, fl = (pl is not None and p == pl), p in flags
        c += np.array([pos and fl, pos and not fl, (not pos) and fl, (not pos) and not fl], int)
    return c


C8_RULES = ('leave', 'jump', 'the noise-free edge')


def c8_flags(r, rule, e0=None):
    """The windows of a C8 run flagged by a rule: 'leave' — the window's run left; 'jump' — it jumped; 'the noise-free edge' — δ(p) at
    or beyond e0 (referee round 76)."""
    w = r['wins']
    if rule == 'leave':
        return [x['p'] for x in w if x['left'] >= 0]
    if rule == 'jump':
        return [x['p'] for x in w if x['maxdev'] > JUMP_C]
    return [x['p'] for x in w if e0 is not None and x['d'] >= e0]


def c8_lines(C8, EDGE=None):
    """C8 (referee round 76): the edge under noise tested in time, run by run — per run kernel 1's first leave (and whether it jumped) and
    the windows each rule flags (the window from period p ejecting the kernel raised by the run's own δ(p): a leave, a jump; the null rule:
    δ(p) at or beyond the noise-free edge, EDGE's bracket midpoint at ω_ref), each rule's first flag against the leave, and per rule over
    the runs of draw 7 (C7's, in sample) and of the other draws (out of sample) the verdicts and the windows' confusion counts."""
    if C8 is None or len(C8) == 0:
        return []
    e0 = {(e['m'], e['ke'], e['metric']): 0.5 * (e['lo'] + e['hi']) for e in (EDGE if EDGE is not None else [])
          if e['k'] == 0 and e.get('om', CT.OM) == CT.OM and e['hi'] is not None}
    out = ['', "### C8. The edge under noise in time, run by run (referee round 76)", '',
           f"C's run (kernel 1) at the configuration and level under a draw of D‴.49's 16 tones, {NP_C} periods, and every {STEP_C8} periods p a "
           f"window: the noise-free loop's fitted state with kernel 1's coefficient raised by the run's own offset δ(p) — its period mean "
           f"over the period before p minus its fitted value —, under the run's own noise from period p, for {NP_C7} periods. A window is "
           f"flagged by the rule 'leave' if it leaves its feature ({LEFT_C:g}), by 'jump' if it jumps (> {JUMP_C:g}); the null rule flags it "
           "if δ(p) is at or beyond the noise-free edge (EDGE). Per run its first leave (j: it jumped), δ in the leave's window, and per rule "
           "the flagged windows (up to four) with the first flag's verdict (hit: in the leave's window; early; miss; for a run that never "
           "left: correct, or a false alarm):", '',
           '| m | k_e | metric | s | draw | leave | δ at the leave\'s window | ' + ' | '.join(C8_RULES) + ' |',
           '|---|---|---|---|---|---|---|' + '---|' * len(C8_RULES)]
    tally = {}
    for r in sorted(C8, key=lambda r: ((r['m'], r['ke'], r['metric']) != C7_CFG, r['m'], r['ke'], r['metric'], r['s'],
                                       r['opts'].get('seed', 7))):
        key = (r['m'], r['ke'], r['metric'])
        P, sd = r['left'][0], r['opts'].get('seed', 7)
        starts = [x['p'] for x in r['wins']]
        pl = None if P < 0 else STEP_C8 * ((P - 1) // STEP_C8)
        dl = next((x['d'] for x in r['wins'] if x['p'] == pl), None)
        cells = []
        for rule in C8_RULES:
            fl = c8_flags(r, rule, e0.get(key))
            v, f = c8_verdict(P, fl)
            t = tally.setdefault((rule, 'in' if sd == 7 else 'out'), dict(v={}, c=np.zeros(4, int)))
            t['v'][v] = t['v'].get(v, 0) + 1
            t['c'] += c8_confusion(P, fl, starts)
            cells.append((', '.join(f"{x:g}" for x in fl[:4]) + (f" … ({len(fl)})" if len(fl) > 4 else '') if fl else 'none') + f": {v}")
        jl = P >= 0 and r['maxdev'][0] > JUMP_C
        out.append(f"| {r['m']:g} | {r['ke']:g} | {r['metric']} | {r['s']:g} | {sd} | " + (f"{P}{' j' if jl else ''}" if P >= 0 else 'none')
                   + ' | ' + ('—' if dl is None else f"{dl:+.4f}") + ' | ' + ' | '.join(cells) + ' |')
    out += ['', "Per rule — in sample (draw 7) and out of sample (the other draws): the runs' verdicts, and the windows' confusion counts "
                "(the leave's window flagged / not flagged, an earlier window flagged / not flagged; a run that never left: all its windows "
                "earlier):", '', '| rule | sample | hit | early | miss | correct | false alarm | leave\'s window flagged | not flagged | '
                "earlier flagged | earlier not flagged |", '|---|---|---|---|---|---|---|---|---|---|---|']
    for rule in C8_RULES:
        for smp in ('in', 'out'):
            t = tally.get((rule, smp))
            if t is None:
                continue
            out.append(f"| {rule} | {smp} | " + ' | '.join(str(t['v'].get(v, 0)) for v in ('hit', 'early', 'miss', 'correct', 'false alarm'))
                       + ' | ' + ' | '.join(str(int(x)) for x in t['c']) + ' |')
    tot = [(smp, tally[(C8_RULES[0], smp)]) for smp in ('in', 'out') if (C8_RULES[0], smp) in tally]
    out += ['', 'The totals: ' + '; '.join(f"{'in sample' if smp == 'in' else 'out of sample'} {sum(t['v'].values())} runs, "
                                          f"{int(t['c'][0] + t['c'][1])} of them leaving, and {int(t['c'][2] + t['c'][3])} earlier windows"
                                          for smp, t in tot) + '.']
    return out + ['']


def report(A=None, B=None, B2=None, C=None, B3=None, B4=None, C2=None, F=None, B5=None, A2=None, C3=None, C4=None, W2=None, EDGE=None,
           C6=None, C7=None, C8=None, C9=None, C9D=None):
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
    L_ += a2_lines(A2)
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
    L_ += b5_lines(B5)
    L_ += b4_lines(B4)
    L_ += f_lines(F, B4)
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
               "| m | k_e | √K max | " + ' | '.join(f"A = {a:g}" for a in CT.ASTAR) + f" | leaves (jumps: the largest distance > {JUMP_C:g}) "
               "| the leave periods: median [min, max] (all levels) |", "|---|---|---|---|---|---|---|---|"]
        for m in M_C:
            for ke in KE_C:
                rs = [r for r in C if r['m'] == m and r['ke'] == ke]
                if not rs:
                    continue
                per_ = [r['left'][q] for r in rs for q in range(3) if r['left'][q] >= 0]
                jmp = sum(r['left'][q] >= 0 and r['maxdev'][q] > JUMP_C for r in rs if 'maxdev' in r for q in range(3))
                L_.append(f"| {m:g} | {ke:g} | {rs[0]['sqrtK']:.2f} | " + ' | '.join(c_cells(rs)) + f" | {len(per_)} ({jmp}) | "
                          + (f"{np.median(per_):.0f} [{min(per_)}, {max(per_)}]" if per_ else '—') + " |")
        lv = []
        for m, ke in ((12.0, 10.0), (3.0, 5.0)):
            byl = leaves_by_level([r for r in C if r['m'] == m and r['ke'] == ke and 'maxdev' in r])
            lv.append(f"m = {m:g}, k_e = {ke:g}: " + '; '.join(
                f"s = {s_:g}: {d['leaves']} ({d['jumps']}) of {d['n']} runs, periods {min(d['periods'])}–{max(d['periods'])}"
                for s_, d in sorted(byl.items()) if d['leaves']))
        L_ += ["", "Per level the leaves over the seeds and kernels (the jumps among them) and their periods — against the slow modes' "
               "time constants (F) —: " + ' | '.join(lv) + "."]
        L_.append("")
    L_ += c2_lines(C2)
    L_ += c3_lines(C3)
    L_ += c4_lines(C3, C4, W2)
    L_ += edge_lines(EDGE, W2, C, C2)
    L_ += c6_lines(C6, C3)
    L_ += c7_lines(C7, W2, C, C2)
    L_ += c8_lines(C8, EDGE)
    L_ += c9_lines(C9, C3, C6, C9D)
    L_.append("Numbers are produced by `experiments/exp_noise_resonance.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


KEYS = ('A', 'A2', 'F', 'B', 'B2', 'B3', 'B5', 'B4', 'C', 'C2', 'C3', 'C4', 'W2', 'EDGE', 'C6', 'C7', 'C8', 'C9', 'C9D')


def main():
    raw = os.path.join(OUT, 'noise_resonance_raw.npz')
    md = os.path.join(OUT, 'noise_resonance.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(**{k: old.get(k) for k in KEYS}))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(KEYS)
    res = {k: v for k, v in old.items() if k not in parts}
    t0 = time.time()
    fns = dict(A=part_A, A2=part_A2, B=part_B, B2=part_B2, B3=part_B3, B4=part_B4, B5=part_B5, C=part_C, C2=part_C2, C3=part_C3,
               F=part_F, C4=part_C4, W2=part_W2, EDGE=part_EDGE, C6=part_C6, C7=part_C7, C8=part_C8, C9=part_C9, C9D=part_C9D)
    for k in [k for k in KEYS if k in parts]:
        res[k] = (fns[k](old.get('B3')) if k == 'B3' else fns[k](res.get('C3'), res.get('C6')) if k == 'C9D' else fns[k]())
        E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
        open(md, 'w').write(report(**{k_: res.get(k_) for k_ in KEYS}))


if __name__ == '__main__':
    main()
