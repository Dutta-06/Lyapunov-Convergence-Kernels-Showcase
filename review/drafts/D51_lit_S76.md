### S76. The transit at a finite mass (Proposition D‴.51): near ties, the overshoot, the creep

*Verdict: the vanishing-initialisation limit and its jump process are known, and smoothing families of the inverse scale space exist;
the free law's finite-mass structure at the walls — its power-law tail below the wall, the near-tie inner problem, the creep and the
migrating losers it causes — is not, as far as we found.* Pesme–Flammarion (NeurIPS 2023, "Saddle-to-saddle dynamics in diagonal linear
networks") and Berthier (JMLR 2023, incremental learning) prove that as the initialisation vanishes the time-rescaled trajectory converges
to a jump process visiting saddles (a LARS-like recursion); the alternating gradient flows of Kunin et al. (NeurIPS 2025) describe the
same small-initialisation regime for two-layer networks as dormant units maximising a utility and active ones minimising the loss, and
prove convergence to gradient flow for diagonal linear networks as the initialisation vanishes. These work in the limit and in general
position; we found no treatment of what a small but finite initialisation does at a near tie (two coordinates reaching the threshold
within the smoothing width), of the coefficient an entering coordinate's dual carries beyond its threshold, or of a margin condition for
the finite-initialisation trajectory to keep the limit's support (web search, October 2026; flag: the appendices of these papers not
re-read). Gissin–Shalev-Shwartz–Daniely (ICLR 2020, "The implicit bias of depth") show a depth separation for incremental learning in
deep linear models — an exponentially small initialisation needed at depth 2, a polynomially small one sufficing at depth ≥ 3 —; the
margin ε^{1/2} of D‴.51(ii) is its counterpart for the free law, whose coefficient below the wall, ∝ ε^{3/2} (D‴.48(o)), is a depth-3
unit's. The linearised Bregman inverse scale space (Osher–Ruan–Xiong–Yao–Yin, ACHA 2016) is a smoothing family of the inverse scale space
with a damping κ (the flow as κ → ∞), with sign-consistency and early-stopping results at finite κ; its coefficients are κ·shrink(z) —
exactly zero below the threshold and linear beyond it —, so it has no tail below the wall, no creep and no migrating losers, and its near
ties are shared linearly, not by a Gibbs split. Exact regularisation (Friedlander–Tseng, SIAM J. Optim. 2007; Yin's analysis of
linearised Bregman, SIAM J. Imaging Sci. 2010): a finite smoothing parameter beyond which the sparse solution is recovered exactly, set
by a Lagrange multiplier (a dual certificate) — the static counterpart of a mass for exact recovery read off the path. The Gibbs split is
the exponentiated-gradient / Hedge weight ratio (Littlestone–Warmuth; Freund–Schapire; Kivinen–Warmuth): two experts whose cumulative gains
differ by δ carry weights in the ratio e^{−ηδ} — here it appears inside the inner layer of the paradigm's transit, the free law's critical
curve being an exponential (entropic) mirror map beyond the wall at the temperature a/(4Gc) in its mirror coordinate (D‴.47 found the
same on branches at 1/(γa₀)). Matched asymptotic expansions are classical (Kevorkian–Cole; Holmes). The creep is the paradigm's own (the
mass in the position law's drift, D‴.48(iii)(a)); its scales, the partners' differential approach creep that shifts the near-tie
threshold, the losers' migration, and the estimate of the mass for exact recovery from the sliding path — with its error measured against
the loop's own bisected thresholds — are new to this program, a finite-mass counterpart of super-resolution's grid/off-the-grid
questions (Duval–Peyré's grid-based BLASSO recovering pairs of spikes around an off-grid one: there by the regularisation, here by the
transit's near ties).
