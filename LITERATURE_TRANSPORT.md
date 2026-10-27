# Literature audit for the transport paradigm (PARADIGM.md)

**Method.** Six independent primary-source sweeps, each a different search angle: optimal transport in
control, mean-field networks, neuro-adaptive control, leakage and permutation symmetry, propagation of
chaos, and Wasserstein geometry. A synthesis agent deduplicated their sources into
`review/LITERATURE_refs.md` (89 entries with URLs; 79 unsourced entries dropped). The per-sweep verdicts
are in `review/LITERATURE_sweeps.md`. The final merge agent hit a session limit, so this document is the
merge, written from the sweep verdicts. It adds sources I verified directly for Theorems G and H, marked
[S1]–[S9], and, in round 3, sources for Propositions B.6, H.7 and H.8, marked [S10]–[S14]; round 5 adds [S15]–[S17] for
Lemma E.0, and round 6 adds [S23]–[S34] for Theorem J′ (bibliographic details from memory, not re-verified offline). Numbers [1]–[89]
refer to `review/LITERATURE_refs.md`.

The audit was adversarial toward novelty. A verdict of "novel" means only that none of the searches
found the construction.

---

## Supplementary sources verified for Theorems G and H

- [S1] P. Marion, R. Berthier. *Leveraging the two-timescale regime to demonstrate convergence of neural networks*. NeurIPS 2023; arXiv:2304.09576.
- [S2] R. Barboni, G. Peyré, F.-X. Vialard. *Ultra-fast feature learning for the training of two-layer neural networks in the two-timescale regime*. arXiv:2504.18208 (2025).
- [S3] L. Chizat. *Mean-field Langevin dynamics: exponential convergence and annealing*. TMLR 2022; arXiv:2202.01009. Assumptions 1–3, Theorem 3.2, Corollary 3.3, Theorem 4.1 and Proposition 5.1 were read in the PDF.
- [S4] A. Nitanda, D. Wu, T. Suzuki. *Convex analysis of the mean field Langevin dynamics*. AISTATS 2022, PMLR 151:9741–9757; arXiv:2201.10469.
- [S5] F. Chen, Z. Ren, S. Wang. *Uniform-in-time propagation of chaos for mean field Langevin dynamics*. arXiv:2212.03050 (v4, 2025).
- [S6] L. Chizat, F. Bach. NeurIPS 2018 (= [29]). Assumptions 3.4(iii), the boundary conditions for the partially 1-homogeneous case, were read in the PDF. They require Neumann conditions on a bounded convex domain, or C¹-uniform behaviour at infinity. Gaussian kernels with free centres and widths satisfy neither (see the PARADIGM.md §6c preamble). Its Theorem 3.5 also needs the support of μ_0 to separate {−r₀}×Θ from {r₀}×Θ; a Dirac start violates this.
- [S7] G. H. Golub, V. Pereyra. *Separable nonlinear least squares: the variable projection method and its applications*. Inverse Problems 19 (2003).
- [S8] A. P. Morgan, K. S. Narendra. *On the stability of nonautonomous differential equations ẋ = [A + B(t)]x with skew symmetric matrix B(t)*. SIAM J. Control Optim. 15(1):163–176 (1977). The round-1 Theorem E judge identified this as the correct citation for the mixed error model. B. D. O. Anderson, IEEE TAC 22 (1977), covers the identifier equations.
- [S9] Standard references: Ambrosio–Gigli–Savaré, *Gradient Flows* (2008); Khalil, *Nonlinear Systems* (3rd ed.); Kokotović–Khalil–O'Reilly, *Singular Perturbation Methods in Control* (1986); Sastry–Bodson, *Adaptive Control* (1989); Graf–Luschgy, *Foundations of Quantization for Probability Distributions* (2000); Holley–Stroock (1987); Bakry–Émery (1985); Otto–Villani, J. Funct. Anal. 173 (2000) (LSI implies Talagrand); Ambrosio–Savaré (2007), Handbook of Differential Equations: Evolutionary Equations 3; Sanders–Verhulst–Murdock, *Averaging Methods in Nonlinear Dynamical Systems* (2nd ed., 2007).
- [S10] Y. Zhang, Q. Cai, Z. Yang, Y. Chen, Z. Wang. *Can temporal-difference and Q-learning learn representation? A mean-field theory*. NeurIPS 2020; arXiv:2006.04761. A mean-field analysis of TD learning with two-layer networks, whose update is not the gradient of any objective.
- [S11] J. N. Tsitsiklis, B. Van Roy. *An analysis of temporal-difference learning with function approximation*. IEEE TAC 42(5):674–690 (1997). TD is a semi-gradient method: it holds the bootstrapped target fixed. The data distribution is that of a fixed policy and does not depend on the parameters.
- [S12] M. Fruchart, R. Hanai, P. B. Littlewood, V. Vitelli. *Non-reciprocal phase transitions*. Nature 592:363–369 (2021). Non-reciprocal interactions in many-body systems.
- [S14] J. C. Perdomo, T. Zrnic, C. Mendler-Dünner, M. Hardt. *Performative prediction*. ICML 2020, PMLR 119; arXiv:2002.06673. Here predictions shift the data distribution, and retraining with the distribution held at its current value converges to performatively stable points rather than optimal ones. This is the closest structural analogue of PARADIGM.md H.7(c).
- [S15] F. Mazenc, M. Malisoff. *Constructions of Strict Lyapunov Functions*. Springer, Communications and Control Engineering, 2009. Strict Lyapunov functions for persistently excited systems built from a backward–forward integral of the PE Gramian. (Bibliographic details from memory; not re-verified against the publisher's record in this session.)
- [S16] F. Mazenc, M. de Queiroz, M. Malisoff. *Uniform global asymptotic stability of a class of adaptively controlled nonlinear systems*. IEEE Trans. Automat. Control 54(5), 2009. (Details from memory; not re-verified here.)
- [S17] A. Loría, E. Panteley. *Uniform exponential stability of linear time-varying systems: revisited*. Systems & Control Letters 47(1):13–24, 2002; and B. D. O. Anderson, *Exponential stability of linear equations arising in adaptive identification*, IEEE Trans. Automat. Control 22(1):83–88, 1977. Explicit UES bounds for persistently excited identifiers. (Details from memory; not re-verified here.)
- [S18] J. Cortés, S. Martínez, T. Karatas, F. Bullo. *Coverage control for mobile sensing networks*. IEEE Trans. Robotics and Automation 20(2):243–255, 2004. Continuous-time Lloyd descent ṗ_i = −k(p_i − C_{V_i}) toward the centroids of Voronoi cells, with its LaSalle convergence proof. (Details from memory; not re-verified here.)
- [S19] M. Balzer, T. Schlömer, O. Deussen. *Capacity-constrained point distributions*. ACM SIGGRAPH 2009; Q. Mérigot, *A multiscale approach to optimal transport*, Computer Graphics Forum 30(5), 2011; F. de Goes, K. Breeden, V. Ostromoukhov, M. Desbrun, *Blue noise through optimal transport*, ACM Trans. Graphics 31(6), 2012. Equal-mass (capacity-constrained) quantisation through optimal transport plans and their barycentres. (Details from memory; not re-verified here.)
- [S20] S. P. Lloyd. *Least squares quantization in PCM*. IEEE Trans. Inform. Theory 28(2):129–137, 1982; Q. Du, V. Faber, M. Gunzburger, *Centroidal Voronoi tessellations: applications and algorithms*, SIAM Review 41(4):637–676, 1999. Lloyd's algorithm and centroidal Voronoi tessellations. (Details from memory; not re-verified here.)
- [S21] A. Saberi, H. Khalil. *Quadratic-type Lyapunov functions for singularly perturbed systems*. IEEE Trans. Automat. Control 29(6):542–550, 1984; and H. K. Khalil, *Nonlinear Systems*, 3rd ed., Prentice Hall, 2002, Ch. 11 (composite Lyapunov functions for two-time-scale systems). (Details from memory; not re-verified here.)
- [S22] D. Aeyels, J. Peuteman. *A new asymptotic stability criterion for nonlinear time-variant differential equations*. IEEE Trans. Automat. Control 43(7):968–971, 1998; J. Peuteman, D. Aeyels, *Exponential stability of slowly time-varying nonlinear systems*, Math. Control Signals Systems 15:202–228, 2002. Stability of slowly varying and averaged systems through the averaged vector field. (Details from memory; not re-verified here.)
- [S13] J. Chevallier. *Uniform decomposition of probability measures: quantization, clustering and rate of convergence*. J. Appl. Probab. 55(4):1037–1045 (2018). Equal-weight approximation rates in W₂ (used in PARADIGM.md Corollary B.6), with Fournier–Guillin, Probab. Theory Relat. Fields 162 (2015), for i.i.d. atoms.

---
- [S23] E. Pardoux, A. Yu. Veretennikov. *On the Poisson equation and diffusion approximation. I*. Ann. Probab. 29(3):1061–1085, 2001 (and *II*, Ann. Probab. 31(3):1166–1192, 2003). Existence, uniqueness and polynomial growth bounds for solutions of the Poisson equation Lu = f − π(f) of a uniformly elliptic diffusion with a recurrence (confining) condition; the tool behind the invariant-measure expansion of Theorem J′. (Details from memory; not re-verified here.)
- [S24] D. Gilbarg, N. S. Trudinger. *Elliptic Partial Differential Equations of Second Order*. Springer, 2nd ed., 1983 (reprinted 2001), Theorems 6.2 and 6.17: interior Schauder estimates and higher regularity. (Details from memory; not re-verified here.)
- [S25] L. Hörmander. *Hypoelliptic second order differential equations*. Acta Math. 119:147–171, 1967; and M. Hairer, *On Malliavin's proof of Hörmander's theorem*, Bull. Sci. Math. 135(6–7):650–666, 2011 (a self-contained probabilistic proof for diffusions with smooth coefficients of bounded derivatives). (Details from memory; not re-verified here.)
- [S26] P. Cattiaux, L. Mesnager. *Hypoelliptic non-homogeneous diffusions*. Probab. Theory Relat. Fields 123:453–483, 2002. Smooth transition densities for time-inhomogeneous diffusions under the parabolic Hörmander condition. (Details from memory; not re-verified here.)
- [S27] D. W. Stroock, S. R. S. Varadhan. *On the support of diffusion processes with applications to the strong maximum principle*. Proc. Sixth Berkeley Symp. Math. Statist. Prob. III:333–359, 1972. The support of the law of a diffusion is the closure of the skeleton (controlled-ODE) trajectories. (Details from memory; not re-verified here.)
- [S28] G. Da Prato, J. Zabczyk. *Ergodicity for Infinite Dimensional Systems*. Cambridge University Press, LMS Lecture Notes 229, 1996, Theorem 4.2.1 (Doob–Khasminskii: strong Feller + irreducible ⇒ at most one invariant measure, and convergence of the laws) and §4.3. (Details from memory; not re-verified here.)
- [S29] U. Seifert. *Stochastic thermodynamics, fluctuation theorems and molecular machines*. Rep. Prog. Phys. 75:126001, 2012; and C. Maes, K. Netočný, *Time-reversal and entropy*, J. Stat. Phys. 110:269–310, 2003. Entropy production rate of a stationary diffusion, ∫jᵀD⁻¹j/ρ, and its vanishing exactly in detailed balance. (Details from memory; not re-verified here.)
- [S30] G. A. Pavliotis, A. M. Stuart. *Multiscale Methods: Averaging and Homogenization*. Springer, Texts in Applied Mathematics 53, 2008, Chapters 10–11 (averaging for SDEs with fast deterministic or stochastic variables via the perturbation expansion of the backward Kolmogorov equation; the solvability conditions used in Theorem J′'s Step 3). (Details from memory; not re-verified here.)
- [S31] J. C. Mattingly, A. M. Stuart, M. V. Tretyakov. *Convergence of numerical time-averaging and stationary measures via Poisson equations*. SIAM J. Numer. Anal. 48(2):552–577, 2010. Perturbation of invariant measures through the Poisson equation of the limiting generator; the structure of Theorem J′'s Step 5 (there for numerical discretisations, here for the averaging of a hypoelliptic closed loop). (Details from memory; not re-verified here.)
- [S32] M. Hairer. *Ergodic properties of Markov processes*. Lecture notes, University of Warwick, 2006 (Lemma 1.13, Corollary 1.14, Theorem 1.16: strong Feller + a point in the support of every transition probability ⇒ uniqueness of the invariant measure; minorisation on compacts). (Details from memory; not re-verified here.)
- [S33] M. Hairer, J. C. Mattingly. *Yet another look at Harris' ergodic theorem for Markov chains*. In: Seminar on Stochastic Analysis, Random Fields and Applications VI, Progr. Probab. 63:109–117, Birkhäuser, 2011. Geometric ergodicity in a weighted total variation norm from a Lyapunov function and a minorisation on its sublevel sets. (Details from memory; not re-verified here.)
- [S34] A. Dembo, O. Zeitouni. *Large Deviations Techniques and Applications*. Springer, 2nd ed., 1998, Theorem 6.2.10 (Sanov's theorem) and Theorem 4.3.1 (Varadhan's lemma). The large-deviation principle for the empirical measure under a mean-field Gibbs measure with a bounded continuous energy, used in Corollary J′.3. (Details from memory; not re-verified here.)
- [S35] R. L. Dobrushin. *Vlasov equations*. Functional Analysis and Its Applications 13:115–123, 1979. Existence, uniqueness and the representation by characteristics of the Vlasov equation with a bounded Lipschitz vector field depending Lipschitz-continuously (bounded-Lipschitz / W₁ metric) on the measure; stability in the measure. Used for the mean-field equation (𝔄₀) of Theorem D‴. (Details from memory; not re-verified here.)
- [S36] C. McDiarmid. *On the method of bounded differences*. In: Surveys in Combinatorics 1989, London Math. Soc. Lecture Note Ser. 141:148–188, Cambridge Univ. Press, 1989. The bounded-differences (McDiarmid) inequality: a function of m independent variables that changes by at most c_j when the j-th is changed satisfies P(f − Ef ≥ t) ≤ exp(−2t²/Σc_j²). Used for the reference's sampling fluctuation in Theorem D‴, Step 8. (Details from memory; not re-verified here.)
- [S37] F. Malrieu. *Logarithmic Sobolev inequalities for some nonlinear PDE's*. Stochastic Processes and their Applications 95(1):109–132, 2001. Uniform-in-time propagation of chaos for McKean–Vlasov diffusions with a uniformly convex confinement and a convex interaction, by synchronous coupling and functional inequalities. The convexity there is global; Theorem D‴ needs it on a neighbourhood of the support of the mean-field flow and gets the interaction's sign from flat convexity. (Details from memory; not re-verified here.)
- [S38] A. Durmus, A. Eberle, A. Guillin, R. Zimmer. *An elementary approach to uniform in time propagation of chaos*. Proceedings of the AMS 148:5387–5398, 2020. Uniform-in-time propagation of chaos by reflection coupling for McKean–Vlasov diffusions whose confinement is convex outside a ball and whose interaction is Lipschitz with a small constant. The interaction of the fitted network is not small in the learning regime; its sign (Theorem D‴, Step 3) replaces the smallness. (Details from memory; not re-verified here.)
- [S39] C. Liu, L. Zhu, M. Belkin. *Loss landscapes and optimization in over-parameterized non-linear systems and neural networks*. Applied and Computational Harmonic Analysis 59:85–116, 2022. The decomposition of the Hessian of a least-squares loss into the Gram (tangent-kernel) term and the residual-weighted second-derivative term, the latter vanishing at interpolation; the PL condition near interpolation. The same decomposition, per particle, is the self term of Theorem D‴'s coupling identity. (Details from memory; not re-verified here.)
- [S40] A. J. Hoffman, J. B. Kruskal. *Integral boundary points of convex polyhedra*. In: Linear Inequalities and Related Systems (H. W. Kuhn, A. W. Tucker, eds.), Annals of Mathematics Studies 38, Princeton, 1956, 223–246 (total unimodularity: with integral right-hand sides every vertex is integral; the transportation constraints are the incidence matrix of a bipartite graph). V. Klee, C. Witzgall. *Facets and vertices of transportation polytopes*. In: Mathematics of the Decision Sciences, Part 1, Lectures in Applied Mathematics 11, AMS, 1968, 257–282 (the vertices of the transportation polytope are the plans supported on spanning forests). A. Schrijver, *Theory of Linear and Integer Programming*, Wiley, 1986, Ch. 19. Used in Proposition D‴.12(ii), (iv). (Details from memory; not re-verified here.)
- [S41] C. Villani. *Optimal Transport: Old and New*. Grundlehren der mathematischen Wissenschaften 338, Springer, 2009, Theorem 5.10 (c-cyclical monotonicity characterises the optimal plans); F. Santambrogio. *Optimal Transport for Applied Mathematicians*. Birkhäuser, 2015, Ch. 1 (the discrete case: the transportation linear program and its vertices). Used in Proposition D‴.12(iii).
- [S42] F. Aurenhammer, F. Hoffmann, B. Aronov. *Minkowski-type theorems and least-squares clustering*. Algorithmica 20:61–76, 1998 (semi-discrete optimal transport: the optimal plan from an absolutely continuous measure to an atomic one is induced by a power diagram — Laguerre cells — with prescribed masses). J. Kitagawa, Q. Mérigot, B. Thibert. *Convergence of a Newton algorithm for semi-discrete optimal transport*. J. Eur. Math. Soc. 21(9):2603–2651, 2019. Used for the mean-field form of the barycentric leak toward an atomic prior (Proposition D‴.12). (Details from memory; not re-verified here.)
- [S43] J. Bolte, A. Daniilidis, A. Lewis. *The Łojasiewicz inequality for nonsmooth subanalytic functions with applications to subgradient dynamical systems*. SIAM J. Optim. 17(4):1205–1223, 2007 (the Kurdyka–Łojasiewicz inequality for Clarke subgradients of subanalytic functions; bounded subgradient trajectories have finite length and converge). J. Bolte, A. Daniilidis, A. Lewis, M. Shiota. *Clarke subgradients of stratifiable functions*. SIAM J. Optim. 18(2):556–572, 2007 (the extension to functions definable in an o-minimal structure). D. Davis, D. Drusvyatskiy, S. Kakade, J. D. Lee. *Stochastic subgradient method converges on tame functions*. Found. Comput. Math. 20:119–154, 2020 (the chain rule for Clarke subgradients of Whitney-stratifiable functions along absolutely continuous curves, hence the descent identity along subgradient trajectories). Used in Proposition D‴.17. (Details from memory; not re-verified here.)
- [S44] K. Kurdyka. *On gradients of functions definable in o-minimal structures*. Ann. Inst. Fourier 48(3):769–783, 1998 (the gradient inequality for C¹ definable functions and the convergence of their gradient trajectories). L. van den Dries, A. Macintyre, D. Marker. *The elementary theory of restricted analytic fields with exponentiation*. Ann. of Math. 140:183–205, 1994 (ℝ_{an,exp} is o-minimal: restricted analytic functions and the exponential are definable). Used in Proposition D‴.17. (Details from memory; not re-verified here.)
- [S45] M. Cuturi. *Sinkhorn distances: lightspeed computation of optimal transport*. NeurIPS 2013; G. Peyré, M. Cuturi. *Computational Optimal Transport*. Foundations and Trends in Machine Learning 11(5–6), 2019 (the entropic plan, its duality and Sinkhorn's algorithm); R. Cominetti, J. San Martín. *Asymptotic analysis of the exponential penalty trajectory in linear programming*. Math. Program. 67:169–187, 1994 (as ε → 0 the entropic solution converges to an optimal plan); N. Courty, R. Flamary, D. Tuia, A. Rakotomamonjy. *Optimal transport for domain adaptation*. IEEE TPAMI 39(9):1853–1865, 2017 (barycentric mappings of entropic plans); A. Ramdas, N. García Trillos, M. Cuturi. *On Wasserstein two-sample testing and related families of nonparametric tests*. Entropy 19(2):47, 2017; A. Genevay, G. Peyré, M. Cuturi. *Learning generative models with Sinkhorn divergences*. AISTATS 2018; J. Feydy, T. Séjourné, F.-X. Vialard, S. Amari, A. Trouvé, G. Peyré. *Interpolating between optimal transport and MMD using Sinkhorn divergences*. AISTATS 2019 (the entropic cost between OT (ε → 0) and an MMD-type discrepancy (ε → ∞), and the debiasing that removes the entropic bias); Y. M. Asano, C. Rupprecht, A. Vedaldi. *Self-labelling via simultaneous clustering and representation learning*. ICLR 2020 (equipartitioned soft assignment by Sinkhorn). Used in Proposition D‴.20. (Details from memory; not re-verified here.)
- [S46] K. Rose, E. Gurewitz, G. C. Fox. *Statistical mechanics and phase transitions in clustering*. Phys. Rev. Lett. 65(8):945–948, 1990; K. Rose. *Deterministic annealing for clustering, compression, classification, regression, and related optimization problems*. Proc. IEEE 86(11):2210–2239, 1998 (a cluster's codevector splits when the temperature falls below twice the largest eigenvalue of the conditional covariance of its data). The precedent of the critical regularisation of Proposition D‴.20. (Details from memory; not re-verified here.)
- [S47] B. T. Polyak. *Gradient methods for minimizing functionals*. Zh. Vychisl. Mat. Mat. Fiz. 3(4):643–653, 1963 (the gradient-dominance inequality ‖∇f‖² ≥ 2μ(f − f*) and the linear convergence of gradient descent under it); H. Karimi, J. Nutini, M. Schmidt. *Linear convergence of gradient and proximal-gradient methods under the Polyak–Łojasiewicz condition*. ECML PKDD 2016 (its place among the conditions weaker than strong convexity). The form of the inequality of Theorem D‴.21. (Details from memory; not re-verified here.)
- [S48] A. Blanchet, J. Bolte. *A family of functional inequalities: Łojasiewicz inequalities and displacement convex functions*. J. Funct. Anal. 275(7):1650–1673, 2018 (Łojasiewicz inequalities in Wasserstein space; for displacement-convex functionals they are equivalent to Talagrand-type inequalities and give the convergence of the Wasserstein gradient flow). Also S. Chewi, T. Maunu, P. Rigollet, A. J. Stromme. *Gradient descent algorithms for Bures–Wasserstein barycenters*. COLT 2020 (a PL inequality for the barycentre functional through optimal maps); Q. Mérigot, F. Santambrogio, C. Sarrazin. *Non-asymptotic convergence bounds for Wasserstein approximation using point clouds*. NeurIPS 2021 (a modified PL inequality for ½W² between a uniform point cloud and a density — the capacity-constrained Lloyd case); and the generalised geodesics of L. Ambrosio, N. Gigli, G. Savaré, *Gradient Flows in Metric Spaces and in the Space of Probability Measures*, Birkhäuser 2008, §9.2 (curves built from two plans with a common base, along which ½W²(·, ν) is 1-convex). The Wasserstein PL theory around Theorem D‴.21: its λ-convex form does not apply (the energy is semiconvex along generalised geodesics only with a negative constant), and the theorem's segments are the generalised geodesic based at the prior. (Details from memory; not re-verified here.)
- [S50] F. Chen, Z. Ren, S. Wang. *Uniform-in-time propagation of chaos for mean field Langevin dynamics*. arXiv:2212.03050, 2022; T. Suzuki, A. Nitanda, D. Wu. *Uniform-in-time propagation of chaos for the mean-field gradient Langevin dynamics*. ICLR 2023 (an O(1/N) error uniform in time under a uniform log-Sobolev inequality, for entropy-regularised mean-field training). The comparison for Remark D‴.22, which gives the coincidence of the long-time limits, not a uniform-in-time estimate. (Details from memory; not re-verified here.)
- [S51] J. Sirignano, K. Spiliopoulos. *Mean field analysis of neural networks: a central limit theorem*. Stochastic Processes and their Applications 130(3):1820–1852, 2020; Z. Chen, G. M. Rotskoff, J. Bruna, E. Vanden-Eijnden. *A dynamical central limit theorem for shallow neural networks*. NeurIPS 2020 (the fluctuations of trained mean-field networks around their limit, of order m^{−1/2} in the parameters and 1/m in the mean-square function, by central limit theorems for the training dynamics; Chen et al. also show that if the mean-field dynamics converges to a measure that interpolates the training data, the asymptotic deviation eventually vanishes in the CLT scaling — the mechanism by which the adaptive prior of Proposition D‴.24 removes the function-level floor, its replicated stationary points being unregularised interpolants). The context of Proposition D‴.24, which identifies the stationary function-level floor of a regularised adaptive law as the weight-slope of the regulariser's bias, and restores the interpolation property under regularisation with a leak toward a learned prior. (Details from memory and, for the interpolation result, from the abstract as found by the round-18 referee; not re-verified here.)
- [S52] K. S. Narendra, A. M. Annaswamy. *A new adaptive law for robust adaptation without persistent excitation*. IEEE TAC 32(2):134–145, 1987 (the e-modification: the leak scaled by |e|, vanishing at perfect tracking). P. A. Ioannou. *General robust adaptive law for tuning controllers*. Proc. 25th IEEE CDC, 474–476, 1986 (the switching σ-modification "guarantees robustness without destroying the ideal properties of the adaptive control schemes": quoted from the abstract, read through a web search in referee round 19), and P. A. Ioannou, J. Sun, *Theory and design of robust direct and indirect adaptive-control schemes*, Int. J. Control 47(3):775–813, 1988 (the comparison of the dead-zone, fixed-σ, ε₁- and switching-σ modifications). W. Magnus. *On the exponential solution of differential equations for a linear operator*. Comm. Pure Appl. Math. 7(4):649–673, 1954; S. Blanes, F. Casas, J. A. Oteo, J. Ros. *The Magnus expansion and some of its applications*. Phys. Rep. 470(5–6):151–238, 2009 (the second term of the expansion is a commutator). R. L. Kosut, B. D. O. Anderson, I. M. Y. Mareels. *Stability theory for adaptive systems: method of averaging and persistency of excitation*. IEEE TAC 32(1):26–34, 1987 (the averaging of adaptive loops; the rate set by the averaged regressor's excitation). The precedents of Corollary D‴.26 and Proposition D‴.27 (bibliographic details checked by a web search; the content from memory except where quoted).
- [S53] R. Z. Khas'minskii. *A limit theorem for the solutions of differential equations with random right-hand sides*. Theory Probab. Appl. 11(3):390–406, 1966 (stochastic averaging). N. Sadegh, R. Horowitz. *Stability and robustness analysis of a class of adaptive controllers for robotic manipulators*. Int. J. Robotics Research 9(3):74–92, 1990 (the desired compensation adaptation law: the adaptation law and the compensation use the desired trajectory, not the measured one). T. Söderström. *Errors-in-Variables Methods in System Identification*. Springer, 2018 (bias compensation, instrumental variables). The precedents of Proposition D‴.30 (bibliographic details checked by a web search; the content from memory).
- [S54] M. A. Hassouneh, H.-C. Lee, E. H. Abed. *Washout filters in feedback control: benefits, limitations and extensions*. Proc. 2004 American Control Conference, paper FrA02.2 (feedback through washout filters preserves the open loop's equilibria; their limitations, and generalised washout filters with coupled auxiliary states). The precedent of Proposition D‴.31's washout reading (bibliographic details checked by a web search; the content from the abstract and from memory).
- [S55] L. Dinh, R. Pascanu, S. Bengio, Y. Bengio. *Sharp minima can generalize for deep nets*. ICML 2017 (rescalings of rectifier networks give equivalent models with arbitrarily sharper minima: the curvature at a minimum is not reparametrisation-invariant). B. Neyshabur, R. Salakhutdinov, N. Srebro. *Path-SGD: path-normalized optimization in deep neural networks*. NeurIPS 2015 (the rescaling invariance of ReLU networks). B. Şimşek, F. Ged, A. Jacot, F. Spadaro, C. Hongler, W. Gerstner, J. Brea. *Geometry of the loss landscape in overparameterized neural networks: symmetries and invariances*. ICML 2021 (symmetry-induced critical points). The precedents of Proposition D‴.28's gauge (the first and the third checked by the round-20 referee through a web search; the content from memory).
- [S56] E. Hairer, C. Lubich, G. Wanner. *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations*. 2nd ed., Springer, 2006 (symmetric methods and reversibility: an expansion in a small parameter is even when the problem is invariant under a symmetry reversing the parameter's sign, as for the error expansion of symmetric integrators). With [S52]'s Blanes–Casas–Oteo–Ros (the time symmetry of the Magnus exponent), the precedent of Corollary D‴.26(c)'s parity argument (from memory).
- [S57] B. Widrow, J. M. McCool, M. G. Larimore, C. R. Johnson Jr. *Stationary and nonstationary learning characteristics of the LMS adaptive filter*. Proc. IEEE 64(8):1151–1162, 1976 (the misadjustment of LMS, μσ²tr R/2, depends on the input's spectrum only through a trace: the mechanism of the η = ∞ floor (N₀/2k_e)tr Q). K. Mayyas, T. Aboulnasr. *Leaky LMS algorithm: MSE analysis for Gaussian data*. IEEE Trans. Signal Process. 45(4):927–934, 1997 (the steady state of a leaky adaptive filter: the fixed prior's factor h_j/(h_j + σ₀)). A. Benveniste, M. Métivier, P. Priouret. *Adaptive Algorithms and Stochastic Approximations*. Springer, 1990 (the Gaussian diffusion approximation of a constant-gain adaptive algorithm and its Lyapunov covariance). R. Carmona, F. Delarue, D. Lacker. *Mean field games with common noise*. Ann. Probab. 44(6):3740–3803, 2016 (a noise common to all particles does not average out in the mean-field limit: the conditional McKean–Vlasov limit). B. D. O. Anderson. *Adaptive systems, lack of persistency of excitation and bursting phenomena*. Automatica 21(3):247–258, 1985 (the drift of the weakly excited directions under disturbances without persistent excitation). The precedents of Proposition D‴.30 named by the round-21 referee (bibliographic details checked by a web search; the content from memory and from the abstracts). A. Dieuleveut, A. Durmus, F. Bach. *Bridging the gap between constant step size stochastic gradient descent and Markov chains*. Ann. Statist. 48(3):1348–1382, 2020, doi:10.1214/19-AOS1850 (constant-step SGD as a Markov chain: on a non-quadratic objective its stationary mean is O(γ) off the optimum — the second-order bias of the iterates' own fluctuation — and Richardson–Romberg extrapolation removes it). M. J. Box. *Bias in nonlinear estimation*. J. R. Stat. Soc. B 33(2):171–190, with discussion 190–201, 1971 (the O(σ²) bias of nonlinear least squares from the model's curvature in its parameters: the analogue of ½D²_Θf̂[Σ_ΘΘ] and its least-squares force). The precedents of Proposition D‴.31's parameters' errors in variables named by the round-22 referee (bibliographic details checked by a web search; the content from the abstracts and from memory).
- [S58] T. Yucelen, W. M. Haddad. *Low-frequency learning and fast adaptation in model reference adaptive control*. IEEE Trans. Automat. Control 58(4):1080–1085, 2013 (a modification term keeping the weights near a low-pass filtered copy of themselves — a leak toward a filtered estimate, which removes the adaptation's high-frequency content and leaves the low-frequency learning: the washout reading). With [S49]'s elastic averaging (the centre a moving average of the workers), the precedents of Proposition D‴.31(i) named by the round-21 referee (bibliographic details checked by a web search; the content from the abstract).
- [S59] K. W. Chang. *Singular perturbations of a general boundary value problem*. SIAM J. Math. Anal. 3(3), 1972, doi:10.1137/0503050 (the decoupling transformation of a linear singularly perturbed system into its slow and fast blocks; the Chang transformation of Kokotović–Khalil–O'Reilly, [S9]). P. C. Moan, J. Niesen. *Convergence of the Magnus series*. Found. Comput. Math. 8(3):291–301, 2008 (the Magnus series of Ẏ = A(t)Y converges when ∫‖A‖₂ < π, the bound sharp). T. Kato. *Perturbation Theory for Linear Operators*. Springer, 1966 (the analytic continuation of a simple eigenvalue). The closest precedents, named by the round-22 referee: P. V. Kokotović, H. K. Khalil, J. O'Reilly, [S9], Ch. 4 (stochastic linear filtering and control: Chang's transformation with white-noise inputs, the slow subsystem's transformed noise and the expansion of its covariance), Ch. 5 (linear time-varying systems) and App. A; A. H. Haddad, P. V. Kokotović. *Stochastic control of linear singularly perturbed systems*. IEEE TAC 22(5):815–821, 1977; X. Yang, J. J. Zhu. *A generalization of Chang transformation for linear time-varying systems*. 49th IEEE CDC, 2010, and *Chang transformation for decoupling of singularly-perturbed linear slowly time-varying systems*. 51st IEEE CDC, 2012 (the transformation's differential Riccati and Sylvester equations solved by a contraction — the structure of (i) and (v); the periodic case is the specialisation in which a bounded solution unique in a ball is periodic); N. Berglund, B. Gentz. *Geometric singular perturbation theory for stochastic differential equations*. J. Differential Equations 191(1):1–54, 2003, doi:10.1016/S0022-0396(03)00020-2, and *Noise-Induced Phenomena in Slow-Fast Dynamical Systems*. Springer, 2006 (stochastic slow manifolds and their covariance). The tools of Lemma D‴.26a (Moan–Niesen, Chang's DOI, Berglund–Gentz and the Yang–Zhu titles checked by a web search, Chang's pages not; Haddad–Kokotović from reference lists; the content from memory and the abstracts).
- [S60] K. Fukunaga, L. D. Hostetler. *The estimation of the gradient of a density function, with applications in pattern recognition*. IEEE Trans. Inf. Theory 21(1):32–40, 1975; Y. Cheng. *Mean shift, mode seeking, and clustering*. IEEE TPAMI 17(8):790–799, 1995, doi:10.1109/34.400568 (the mean-shift step is the gradient of the kernel density with the shadow kernel; k-means-like clustering as a special case). M. Á. Carreira-Perpiñán, C. K. I. Williams. *On the number of modes of a Gaussian mixture*. Scale-Space Methods in Computer Vision, LNCS 2695, 625–640, 2003, and *An isotropic Gaussian mixture can have more modes than components*, EDI-INF-RR-0185, Edinburgh, 2003 (three equal isotropic Gaussians at the vertices of an equilateral triangle have a mode at the centre for a range of variances: extra modes need overlapping components). J. A. Carrillo, R. J. McCann, C. Villani. *Kinetic equilibration rates for granular media and related equations: entropy dissipation and mass transportation estimates*. Rev. Mat. Iberoam. 19(3):971–1018, 2003; A. L. Bertozzi, T. Laurent, J. Rosado. *L^p theory for the multidimensional aggregation equation*. Comm. Pure Appl. Math. 64(1):45–83, 2011, doi:10.1002/cpa.20334 (aggregation equations ∂_tu + ∇·(uv) = 0, v = −∇K*u, as Wasserstein gradient flows of interaction energies). R. Hegselmann, U. Krause. *Opinion dynamics and bounded confidence: models, analysis and simulation*. JASSS 5(3), 2002 (agents averaging with those within a confidence radius form clusters at least the radius apart). L. Wu, D. Wang, Q. Liu. *Splitting steepest descent for growing neural architectures*. NeurIPS 2019, arXiv:1910.02366 (a neuron is split along the minimum eigenvector of its splitting matrix, the residual-weighted Hessian of the neuron, when that eigenvalue is negative — the residual curvature R_k of Proposition D‴.5 with the sign convention of a loss). J. Lian, Y. Lee, S. D. Sudhoff, S. H. Żak. *Self-organizing radial basis function network for real-time approximation of continuous-time dynamical systems*. IEEE Trans. Neural Netw. 19(3):460–474, 2008 (RBF units added and removed on line by accuracy criteria, with a Lyapunov bound on the approximation error). R. Kulhavý. *Restricted exponential forgetting in real-time identification*. Automatica 23(5):589–600, 1987 (directional forgetting: forgetting only in the excited directions). L. Bungert, T. Roith, P. Wacker. *Polarized consensus-based dynamics for optimization and sampling*. Math. Program. 211:125–155, 2025, doi:10.1007/s10107-024-02095-y, arXiv:2211.05238 (each particle drifts toward its own kernel-weighted local mean: a bounded-confidence dynamics with a common objective); B. Geshkovski, C. Letrouit, Y. Polyanskiy, P. Rigollet. *The emergence of clusters in self-attention dynamics*. NeurIPS 2023, arXiv:2305.05465 (particles under softmax-weighted mutual attraction cluster). K. Fukumizu, S. Amari. *Local minima and plateaus in hierarchical structures of multilayer perceptrons*. Neural Netw. 13(3):317–327, 2000 (duplicating a hidden unit turns a critical point of the smaller network into lines of critical points of the larger one, with segments of minima and saddles: the twins' curves of (v)). T. D. Hocking, A. Joulin, F. Bach, J.-P. Vert. *Clusterpath: an algorithm for clustering using convex fusion penalties*. ICML 2011, 745–752 (convex clustering, with Gaussian-kernel weights); H. Ohlsson, L. Ljung, S. Boyd. *Segmentation of ARX-models using sum-of-norms regularization*. Automatica 46(6):1107–1111, 2010 (one constant trades the fit against the number of distinct parameter values). K. Rose. *Deterministic annealing for clustering, compression, classification, regression, and related optimization problems*. Proc. IEEE 86(11):2210–2239, 1998 (a cluster splits when the annealing parameter crosses a critical value set by the cluster's covariance; the annealed scale of Numerics (C) is such a schedule). R. Kulhavý, M. B. Zarrop. *On a general concept of forgetting*. Int. J. Control 58(4):905–924, 1993 (stabilized forgetting: forgetting toward a prior on the unexcited directions, the analogue of Theorem D‴.23's adaptive prior rather than of this leak). The precedents of Proposition D‴.32 (Cheng, Carreira-Perpiñán–Williams, Bertozzi–Laurent–Rosado, Wu–Wang–Liu, Lian et al., Bungert–Roith–Wacker, Geshkovski et al., Fukumizu–Amari, Hocking et al., Ohlsson–Ljung–Boyd and Rose checked by a web search, the contents of the last five from abstracts and citing texts; Kulhavý–Zarrop's details checked, its content from memory; the others from memory).
- [S61] J. Teramae, D. Tanaka. *Robustness of the noise-induced phase synchronization in a general class of limit cycle oscillators*. Phys. Rev. Lett. 93:204103, 2004, doi:10.1103/PhysRevLett.93.204103 (uncoupled identical oscillators driven by a common weak noise phase-lock: the Lyapunov exponent of the synchronised state, derived to second order in the noise, is negative). D. S. Goldobin, A. Pikovsky. *Synchronization and desynchronization of self-sustained oscillators by common noise*. Phys. Rev. E 71:045201, 2005, doi:10.1103/PhysRevE.71.045201 (weak common noise makes the largest Lyapunov exponent negative; a strong one can desynchronise nonisochronous oscillators). H. Nakao, K. Arai, Y. Kawamura. *Noise-induced synchronization and clustering in ensembles of uncoupled limit-cycle oscillators*. Phys. Rev. Lett. 98:184101, 2007, doi:10.1103/PhysRevLett.98.184101 (common and independent noise: the stationary distribution of the phase differences, clustering). The precedents of Corollary D‴.33's frame (bibliographic details checked by a web search; the content from the abstracts).
- [S62] M. Liero, A. Mielke, G. Savaré. *Optimal entropy-transport problems and a new Hellinger–Kantorovich distance between positive measures*. Invent. Math. 211(3):969–1117, 2018, doi:10.1007/s00222-017-0759-8 (the Hellinger–Kantorovich, or Wasserstein–Fisher–Rao, distance: transport and creation of mass); L. Chizat. *Sparse optimization on measures with over-parameterized gradient descent*. Math. Program. 194:487–532, 2022, doi:10.1007/s10107-021-01636-z (conic particle gradient descent on the positions and the weights of particles; local exponential and global convergence for sparse optimisation on measures under non-degeneracy); G. Rotskoff, S. Jelassi, J. Bruna, E. Vanden-Eijnden. *Neuron birth-death dynamics accelerates gradient descent and converges asymptotically*. ICML 2019, PMLR 97:5508–5517 (a non-local mass-transport term added to the mean-field dynamics of a one-hidden-layer network, the same minimiser, faster convergence); J. Kivinen, M. K. Warmuth. *Exponentiated gradient versus gradient descent for linear predictors*. Inform. and Comput. 132(1):1–63, 1997 (the multiplicative EG update); N. M. Boffi, J.-J. E. Slotine. *Implicit regularization and momentum algorithms in nonlinearly parameterized adaptive control and prediction*. Neural Comput. 33(3):590–673, 2021, doi:10.1162/neco_a_01360 (adaptation laws of mirror-descent type with a Bregman divergence in place of the quadratic parameter error in the Lyapunov function); J. Hofbauer, K. Sigmund. *Evolutionary Games and Population Dynamics*. Cambridge University Press, 1998 (the replicator equation; the relative entropy as its Lyapunov function at an evolutionarily stable state). S. Kondratyev, L. Monsaingeon, D. Vorotnikov. *A new optimal transport distance on the space of finite Radon measures*. Adv. Differential Equations 21(11/12):1117–1164, 2016, arXiv:1505.07746, and L. Chizat, G. Peyré, B. Schmitzer, F.-X. Vialard. *An interpolating distance between optimal transport and Fisher–Rao metrics*. Found. Comput. Math. 18(1):1–44, 2018, doi:10.1007/s10208-016-9331-y (the other two origins of the Wasserstein–Fisher–Rao distance: transport with a Fisher–Rao source term); L. Chizat, F. Bach. *On the global convergence of gradient descent for over-parameterized models using optimal transport*. NeurIPS 2018 (gradient descent on the positions and the square-root weights of a two-layer network as a Wasserstein–Fisher–Rao-type flow: the lift (√β, Θ) of (i)); S. Shahshahani. *A new mathematical framework for the study of linkage and selection*. Mem. Amer. Math. Soc. 17(211), 1979 (the Riemannian metric on the simplex in which the replicator equation of selection is a gradient flow); A. P. Dempster, N. M. Laird, D. B. Rubin. *Maximum likelihood from incomplete data via the EM algorithm*. J. R. Stat. Soc. B 39(1):1–38, 1977 (EM's multiplicative update of mixture weights). (ii) is the negative-entropy instance, on the simplex, of the mirror-descent adaptation laws of Boffi–Slotine. The precedents of Proposition D‴.34 (Liero–Mielke–Savaré, Chizat, Rotskoff et al., Kivinen–Warmuth, Boffi–Slotine, Kondratyev–Monsaingeon–Vorotnikov, Chizat–Peyré–Schmitzer–Vialard (its volume from memory) and Shahshahani checked by a web search; Hofbauer–Sigmund, Chizat–Bach and Dempster–Laird–Rubin from memory).
- [S63] M. Golubitsky, I. Stewart. *The Symmetry Perspective: From Equilibrium to Chaos in Phase Space and Physical Space*. Progress in Mathematics 200, Birkhäuser, 2002 (S_N-equivariant bifurcation, the invariant theory of S_N, and the speciation model: identical coupled phenotypes split by a symmetry-breaking bifurcation whose generic form, through the standard representation's quadratic equivariant, is a jump); T. Elmhirst. *S_N-equivariant symmetry-breaking bifurcations*. Int. J. Bifurcation and Chaos 14(3):1017–1036, 2004 (the normal form of the Liapunov–Schmidt reduced S_N-equivariant problem, the branches with their stability conditions); Y. Arjevani, M. Field. *Equivariant bifurcation, quadratic equivariants, and symmetry breaking for the standard representation of S_k*. Nonlinearity 35(6):2809–2857, 2022, arXiv:2107.02422 (generic S_k-equivariant bifurcation with quadratic equivariants, extending Ihrig–Golubitsky and Chossat–Lauterbach–Melbourne; motivated by the spurious minima of student–teacher networks). I. Stewart, T. Elmhirst, J. Cohen. *Symmetry-breaking as an origin of species*. In *Bifurcation, Symmetry and Patterns* (J. Buescu, S. B. S. D. Castro, A. P. S. Dias, I. S. Labouriau, eds.), Birkhäuser, 2003, 3–54 (speciation as an all-to-all coupled S_N-symmetric jump bifurcation); J. Brea, B. Şimşek, B. Illing, W. Gerstner. *Weight-space symmetry in deep networks gives rise to permutation saddles, connected by equal-loss valleys across the loss landscape*. arXiv:1907.02911, 2019; B. Şimşek, F. Ged, A. Jacot, F. Spadaro, C. Hongler, W. Gerstner, J. Brea. *Geometry of the loss landscape in overparameterized neural networks: symmetries and invariances*. ICML 2021 (the symmetry-induced critical points of networks with duplicated units, of which the replicated class is one). The precedents of Proposition D‴.35 (Golubitsky–Stewart, Elmhirst, Arjevani–Field and, by round 25's referee, Stewart–Elmhirst–Cohen checked by a web search for their bibliographic details and abstracts; the Golubitsky–Stewart book's speciation chapter from its table of contents; Brea et al. and Şimşek et al. from memory; Ihrig–Golubitsky and Chossat–Lauterbach–Melbourne not seen).
- [S64] G. S. Katzenberger. *Solutions of a stochastic differential equation forced onto a manifold by a large drift*. Ann. Probab. 19(4):1587–1628, 1991, doi:10.1214/aop/1176990225 (semimartingales driven by a drift that pushes them onto a manifold of fixed points of a deterministic flow: relative compactness in the Skorohod topology and the limit's stochastic integral equation, through the flow's limit map Φ — its first derivative projecting the noise onto the manifold, its second giving the drift). G. Blanc, N. Gupta, G. Valiant, P. Valiant. *Implicit regularization for deep neural networks driven by an Ornstein–Uhlenbeck like process*. COLT 2020, PMLR 125:483–513 (SGD with label noise near a zero-loss parameter: a local implicit regulariser, the sum over the data of the model's squared parameter gradients). Z. Li, T. Wang, S. Arora. *What happens after SGD reaches zero loss? A mathematical framework*. ICLR 2022, arXiv:2110.06914 (after reaching the manifold of minimisers, small-step SGD follows a Katzenberger-limit SDE on that manifold, for arbitrary noise covariance, on the time scale η⁻²). T. L. Parsons, T. Rogers. *Dimension reduction for stochastic dynamical systems forced onto a manifold by large drift: a constructive approach with examples from theoretical biology*. J. Phys. A 50(41):415601, 2017, arXiv:1510.07031 (the projection computed constructively for population-genetics models). G. W. A. Constable, A. J. McKane, T. Rogers. *Stochastic dynamics on slow manifolds*. J. Phys. A 46(29):295002, 2013, arXiv:1301.7697 (the projection of noise onto a slow manifold of a fast–slow system, with the noise-induced drift). The precedents of Proposition D‴.36 (Katzenberger, Blanc et al. and Li–Wang–Arora checked by a web search for their bibliographic details and abstracts; Parsons–Rogers and Constable–McKane–Rogers by round 26's referee; Constable–McKane–Rogers' journal details from memory; Funaki, Kurtz, Berglund–Gentz and the entropy-regularised replicator not seen).
- [S65] R. Z. Khasminskii. *On the principle of averaging the Itô's stochastic differential equations*. Kybernetika 4(3):260–279, 1968 (the averaging principle: a slow variable driven by a fast ergodic diffusion follows the drift averaged over the fast invariant law); G. A. Pavliotis, A. M. Stuart. *Multiscale Methods: Averaging and Homogenization*. Texts in Applied Mathematics 53, Springer, 2008 (averaging for SDEs, the effective dynamics of slow variables). J. G. Kirkwood. *Statistical mechanics of fluid mixtures*. J. Chem. Phys. 3:300–313, 1935 (the potential of mean force: a slow coordinate in equilibrium with fast ones feels the free energy −τ log Z, entropy included). S. Hochreiter, J. Schmidhuber. *Flat minima*. Neural Computation 9(1):1–42, 1997 (flat minima — large regions of nearly constant error — as simple networks); P. Chaudhari, A. Choromanska, S. Soatto, Y. LeCun, C. Baldassi, C. Borgs, J. Chayes, L. Sagun, R. Zecchina. *Entropy-SGD: biasing gradient descent into wide valleys*. ICLR 2017, arXiv:1611.01838 (the local entropy, a log-partition function of a Langevin inner loop, as an objective favouring flat regions). The precedents of Proposition D‴.37 (Khasminskii, Hochreiter–Schmidhuber and Chaudhari et al. checked by a web search for their bibliographic details and abstracts; Pavliotis–Stuart and Kirkwood from memory).
- [S66] M. Fixman. *Classical statistical mechanics of constraints: a theorem and application to polymers*. Proc. Natl. Acad. Sci. USA 71(8):3050–3053, 1974, doi:10.1073/pnas.71.8.3050 (holonomic constraints make the kinetic energy depend on the coordinates; the metric determinant that results, and its evaluation for distance constraints — the origin of the "Fixman potential", the log-determinant correction that makes a rigidly constrained system sample like the stiff flexible one). E. J. Hinch. *Brownian motion with stiff bonds and rigid constraints*. J. Fluid Mech. 271:219–234, 1994, doi:10.1017/S0022112094001746, and P. S. Grassia, E. J. Hinch, L. C. Nitsche. *Computer simulations of Brownian motion of complex systems*. J. Fluid Mech. 282:373–403, 1995 (the Brownian dynamics of stiff bonds against rigid constraints and the drift corrections that reconcile them). G. Ciccotti, T. Lelièvre, E. Vanden-Eijnden. *Projection of diffusions on submanifolds: application to mean force computation*. Comm. Pure Appl. Math. 61(3):371–408, 2008, doi:10.1002/cpa.20210 (constrained diffusions ergodic for the Boltzmann–Gibbs law restricted to a submanifold, and the free energy's gradient — the mean force — from them). The precedents of Proposition D‴.38 (the three checked by a web search for their bibliographic details and abstracts — Hinch's by its listing only; Grassia–Hinch–Nitsche from the same search; a 2024 preprint, arXiv:2403.13185, disputes the completeness of Fixman's limiting distribution, not seen beyond its abstract). Further precedents (referee round 28, by search at abstract or citation level unless marked): adaptive biasing force — E. Darve, A. Pohorille, J. Chem. Phys. 115:9169, 2001 (from memory); T. Lelièvre, M. Rousset, G. Stoltz, *Long-time convergence of an adaptive biasing force method*, Nonlinearity 21:1155–1181, 2008 — the online removal of a free-energy drift along a slow coordinate; well-tempered metadynamics' bias factor, A. Barducci, G. Bussi, M. Parrinello, Phys. Rev. Lett. 100:020603, 2008 (from memory) — a gain on such a bias; the log-barrier, its maximiser and the max-det selection — G. Sonnevend (1986) and Y. Nesterov, A. Nemirovski (1994) (from memory), L. Vandenberghe, S. Boyd, S.-P. Wu, *Determinant maximization with linear matrix inequality constraints*, SIAM J. Matrix Anal. Appl. 19(2):499–533, 1998; the stiff limit's zero-point (log-det) term — E. Helfand, J. Chem. Phys. 71:5000, 1979 (via citations) — and, in the overdamped setting, R. Waszkiewicz, M. Lisicki, *The trimer paradox: the effect of stiff constraints on equilibrium distributions in overdamped dynamics*, J. Chem. Phys. 2025 (arXiv:2403.13185; abstract only) — for J quadratic in the residual, ℓ is the log pseudo-determinant of the constraint metric up to a constant, which reconciles the two names; noise-induced (Itô) drifts of the replicator — D. Fudenberg, C. Harris, J. Econ. Theory 57:420–441, 1992 (from memory), P. Mertikopoulos, Y. Viossat, *Imitation dynamics with payoff shocks*, Int. J. Game Theory 45:291–320, 2016.
- [S67] J. C. Mattingly, A. M. Stuart, D. J. Higham. *Ergodicity for SDEs and approximations: locally Lipschitz vector fields and degenerate noise*. Stochastic Process. Appl. 101(2):185–232, 2002, doi:10.1016/S0304-4149(02)00150-3 (ergodicity from a minorisation and a Lyapunov condition — the dissipativity a confining potential provides —, gradient systems and the Langevin equation among the examples). M. Raginsky, A. Rakhlin, M. Telgarsky. *Non-convex learning via Stochastic Gradient Langevin Dynamics: a nonasymptotic analysis*. COLT 2017, PMLR 65:1674–1703 (the Langevin dynamics of learning under a dissipativity assumption on the loss, converging to its Gibbs law). M. I. Freidlin, A. D. Wentzell. *Random Perturbations of Dynamical Systems*. Grundlehren 260, Springer, 3rd ed. 2012, doi:10.1007/978-3-642-25847-3 (the exit problem: the exit time's exponential order is the quasi-potential's minimum on the boundary, for a gradient diffusion the energy's); A. Bovier, M. Eckhoff, V. Gayrard, M. Klein. *Metastability in reversible diffusion processes I: sharp asymptotics for capacities and exit times*. J. Eur. Math. Soc. 6(4):399–424, 2004, doi:10.4171/JEMS/14 (metastable exit times through capacities; the Eyring–Kramers formula). P. A. Ioannou, J. Sun. *Robust Adaptive Control*. Prentice Hall, 1996 (parameter projection, σ-modification and the dead zone as robustness modifications of adaptive laws); J.-B. Pomet, L. Praly. *Adaptive nonlinear regulation: estimation from the Lyapunov equation*. IEEE Trans. Automat. Control 37(6):729–740, 1992 (usually credited with the Lipschitz projection, a boundary layer around the convex set); E. Lavretsky, T. E. Gibson, A. M. Annaswamy. *Projection operator in adaptive systems*. arXiv:1112.4232 (a survey of the projection operator). A. Durmus, É. Moulines, M. Pereyra. *Efficient Bayesian computation by proximal Markov chain Monte Carlo: when Langevin meets Moreau*. SIAM J. Imaging Sci. 11(1):473–506, 2018, arXiv:1612.07471 (MYULA: the unadjusted Langevin algorithm on the Moreau–Yosida envelope of a nonsmooth potential); N. Brosse, A. Durmus, É. Moulines, M. Pereyra. *Sampling from a log-concave distribution with compact support with proximal Langevin Monte Carlo*. COLT 2017, PMLR 65:319–342, arXiv:1705.08964 (a convex body's indicator replaced by its Moreau–Yosida envelope — the squared distance to the set, a box's dead zone —, with explicit convergence bounds; both search-verified, added after referee round 31). The precedents of Proposition D‴.39 (Mattingly–Stuart–Higham, Raginsky–Rakhlin–Telgarsky, Bovier et al. and the Freidlin–Wentzell edition checked by a web search for their bibliographic details and abstracts; Pomet–Praly's details from the search, its credit for the Lipschitz projection the usual one, not checked against its text; Ioannou–Sun and the survey from memory).
- [S69] E. Wong, M. Zakai. *On the convergence of ordinary integrals to stochastic integrals*. Ann. Math. Statist. 36(5):1560–1564, 1965 (integrals driven by smooth approximations of Brownian motion converge to Stratonovich integrals: the white limit of physical noise is Stratonovich's). N. G. van Kampen. *Itô versus Stratonovich*. J. Stat. Phys. 24:175–187, 1981, doi:10.1007/BF01007642 (the nonlinear Langevin equation needs an interpretation rule; for external noise Stratonovich's applies). W. Horsthemke, R. Lefever. *Noise-Induced Transitions*. Springer Series in Synergetics 15, 1984 (noise-induced drifts under multiplicative noise). G. Blanc, N. Gupta, G. Valiant, P. Valiant. *Implicit regularization for deep neural networks driven by an Ornstein–Uhlenbeck like process*. COLT 2020, PMLR 125:483–513 (labels perturbed by fresh independent noise at each step: near zero training error the dynamics follow an implicit regulariser, the sum over the training points of the model's squared parameter gradient). A. Damian, T. Ma, J. D. Lee. *Label noise SGD provably prefers flat global minimizers*. NeurIPS 2021 (arXiv:2106.06530). J. Z. HaoChen, C. Wei, J. D. Lee, T. Ma. *Shape matters: understanding the implicit bias of the noise covariance*. COLT 2021, PMLR 134:2315–2357 (parameter-dependent noise biases toward minima with smaller noise variance, spherical noise does not). A. Shalova, A. Schlichting, M. Peletier. *Singular-limit analysis of gradient descent with noise injection*. J. Mach. Learn. Res. 27, no. 146, 2026 (arXiv:2404.12293; abstract and summaries only: label noise, dropout and minibatch noise evolve along the zero-loss set on different time scales, minibatch noise trivially). Q. Li, C. Tai, W. E. *Stochastic modified equations and dynamics of stochastic gradient algorithms I: mathematical foundations*. J. Mach. Learn. Res. 20(40):1–47, 2019 (SGD's weak approximation by Itô SDEs). Z. Li, S. Malladi, S. Arora. *On the validity of modeling SGD with stochastic differential equations (SDEs)*. NeurIPS 2021. A. Orvieto, H. Kersting, F. Proske, F. Bach, A. Lucchi. *Anticorrelated noise injection for improved generalization*. ICML 2022, PMLR 162:17094–17116 (anti-correlated perturbations of the parameters favour wide minima — a perturbed-gradient mechanism, not the error's noise). The precedents of Proposition D‴.41 (the bibliographic details checked by a web search; the contents at abstract level, Shalova–Schlichting–Peletier's label-noise theorem not seen).
- [S70] M. E. Tipping. *Sparse Bayesian learning and the relevance vector machine*. J. Mach. Learn. Res. 1:211–244, 2001 (automatic relevance determination: a Gaussian prior with one variance per basis function, fitted by the marginal likelihood, whose −log is log|σ²I + ΦΓΦᵀ| + yᵀ(σ²I + ΦΓΦᵀ)⁻¹y; from memory). D. P. Wipf, B. D. Rao. *Sparse Bayesian learning for basis selection*. IEEE Trans. Signal Process. 52(8):2153–2164, 2004, and *ℓ₀-norm minimization for basis selection*. NeurIPS 17, 2004 (the SBL cost keeps its global minimum at the maximally sparse solution with fewer troublesome local minima than FOCUSS-type penalties; search-verified at abstract level, the journal paper's volume and pages from memory). M. Fazel, H. Hindi, S. Boyd. *Log-det heuristic for matrix rank minimization with applications to Hankel and Euclidean distance matrices*. Proc. American Control Conference, 2003 (log det(X + δI) as a smooth concave surrogate of rank, minimised by iterated trace minimisations whose first step is the trace heuristic; search-verified). E. J. Candès, M. B. Wakin, S. P. Boyd. *Enhancing sparsity by reweighted ℓ₁ minimization*. J. Fourier Anal. Appl. 14:877–905, 2008 (the log penalty Σ log(|w_i| + ε) by reweighted ℓ¹; from memory). D. Wipf, S. Nagarajan. *Iterative reweighted ℓ₁ and ℓ₂ methods for finding sparse solutions*. IEEE J. Sel. Topics Signal Process. 4(2):317–329, 2010, doi:10.1109/JSTSP.2010.2042413 (SBL's non-separable penalty as reweighted ℓ¹ with weights [φ_iᵀ(αI + ΦW̃X̃Φᵀ)⁻¹φ_i]^{1/2}, their eq. (32); Theorems 1–2: such non-separable reweighting never does worse than ℓ¹ and succeeds where ℓ¹ fails on a set of positive measure — the phenomenon the slow selection shows; read from the authors' PDF by referee round 32). J. Z. HaoChen, C. Wei, J. D. Lee, T. Ma. *Shape matters: understanding the implicit bias of the noise covariance*. COLT 2021, PMLR 134 (label-noise SGD on a quadratically parametrised linear model recovers the sparse ground truth, Gaussian parameter noise does not; search-verified). S. Pesme, L. Pillaud-Vivien, N. Flammarion. *Implicit bias of SGD for diagonal linear networks: a provable benefit of stochasticity*. NeurIPS 2021, arXiv:2106.09524 (SGD's noise as a smaller effective initialisation; search-verified). L. Pillaud-Vivien, J. Reygner, N. Flammarion's label-noise Lasso (COLT 2022) and Damian–Ma–Lee's finite-step regulariser are [S69]'s; Woodworth et al.'s kernel and rich regimes [S72]'s. F. Chen, D. Kunin, A. Yamamura, S. Ganguli. *Stochastic collapse: how gradient noise attracts SGD dynamics towards simpler subnetworks*. NeurIPS 2023, arXiv:2306.04251 (invariant sets of vanishing or identical neurons attract under gradient noise; search-verified). L. Chizat. *Sparse optimization on measures with over-parameterized gradient descent*. Math. Program. 194:487–532, 2022 (the conic, u²-parametrised Wasserstein–Fisher–Rao particle flows for sparse deconvolution on measures; from memory). L. Chizat, F. Bach. *On the global convergence of gradient descent for over-parameterized models using optimal transport*. NeurIPS 2018, arXiv:1805.09545 (search-verified). J.-J. Fuchs. *On sparse representations in arbitrary redundant bases*. IEEE Trans. Inf. Theory 50(6):1341–1344, 2004, and J. A. Tropp, *Greed is good*, IEEE Trans. Inf. Theory 50(10):2231–2242, 2004 (the dual certificate / exact recovery condition for ℓ¹; from memory). The precedents of Proposition D‴.43 (a web search found no work on the temporal correlation of label noise as a sparsity prior, nor on a log-determinant regulariser from coloured label noise; closest: the ℓ¹/Lasso selection of white label noise, sparse Bayesian learning's log-determinant, the log-det heuristic).
- [S72] Z. Li, T. Wang, J. D. Lee, S. Arora. *Implicit bias of gradient descent on reparametrized models: on equivalence to mirror descent*. NeurIPS 2022, arXiv:2207.04036 (gradient flow under a commuting parametrisation is a mirror flow whose mirror map depends only on the start and the parametrisation, not on the loss; conversely every mirror flow is such a gradient flow). S. Marcotte, R. Gribonval, G. Peyré. *Abide by the law and follow the flow: conservation laws for gradient flows*. NeurIPS 2023, arXiv:2307.00144 (the conserved functions of gradient flows of reparametrised models — one per unit for a 2 → 1 map; finite steps break exact conservation). S. Azulay, E. Moroshko, M. S. Nacson, B. Woodworth, N. Srebro, A. Globerson, D. Soudry. *On the implicit bias of initialization shape: beyond infinitesimal mirror descent*. ICML 2021, PMLR 139:468–477, arXiv:2102.09769 (u ⊙ v with δ = u² − v² conserved; the initialisation's shape sets the separable potential). B. Woodworth, S. Gunasekar, J. D. Lee, E. Moroshko, P. Savarese, I. Golan, D. Soudry, N. Srebro. *Kernel and rich regimes in overparametrized models*. COLT 2020, PMLR 125:3635–3673, arXiv:1906.05827 (the initialisation's scale moves the implicit bias from the kernel (ℓ²) to the rich (ℓ¹) regime; depth D gives polynomial rich limits). U. Ghai, E. Hazan, Y. Singer. *Exponentiated gradient meets gradient descent*. ALT 2020, PMLR 117 (the hyperbolic entropy interpolating ℓ² and the entropy). E. Amid, M. K. Warmuth. *Reparameterizing mirror descent as gradient descent*. NeurIPS 2020, arXiv:2002.10487 (continuous unnormalised exponentiated gradient as gradient descent on a reparametrisation). J. Kivinen, M. K. Warmuth. *Exponentiated gradient versus gradient descent for linear predictors*. Information and Computation 132(1):1–63, 1997 (the EG and unnormalised EGU updates; from memory). S. Gunasekar, J. Lee, D. Soudry, N. Srebro. *Characterizing implicit bias in terms of optimization geometry*. ICML 2018 (mirror descent on underdetermined least squares converges to the start's Bregman projection; the stochastic version keeps it). N. Azizan, B. Hassibi. *Stochastic gradient/mirror descent: minimax optimality and implicit regularization*. ICLR 2019 (stochastic mirror descent's limit is the start's Bregman projection; from memory). N. M. Boffi, J.-J. E. Slotine. *Implicit regularization and momentum algorithms in nonlinearly parameterized adaptive control and prediction*. Neural Computation 33(3), 2021, arXiv:1912.13154 (mirror-descent adaptation laws and their Bregman Lyapunov functions; implicit regularisation among the parameters consistent with the data). S. Pesme, L. Pillaud-Vivien, N. Flammarion. *Implicit bias of SGD for diagonal linear networks: a provable benefit of stochasticity*. NeurIPS 2021, arXiv:2106.09524 (SGD's noise, through the Itô stochastic gradient flow, shrinks the effective initialisation and selects a sparser zero-loss solution). M. Even, S. Pesme, S. Gunasekar, N. Flammarion. *(S)GD over diagonal linear networks: implicit bias, large stepsizes and edge of stability*. NeurIPS 2023 (the per-step change of the effective initialisation in closed form). J. Z. HaoChen, C. Wei, J. D. Lee, T. Ma. *Shape matters: understanding the implicit bias of the noise covariance*. COLT 2021, PMLR 134 (parameter-dependent label noise selects, isotropic Gaussian noise does not). R. Cominetti, J. San Martín (Math. Program. 67, 1994) and J. Weed (COLT 2018): the exponential (entropic) penalty in linear programming, approaching a vertex exponentially in the inverse temperature. The precedents of Proposition D‴.45 (search-verified for their bibliographic details and abstracts by this round's referee or by us, their proofs not seen; Kivinen–Warmuth and Azizan–Hassibi from memory).
- [S73] P. A. Ioannou, P. V. Kokotović. *Adaptive Systems with Reduced Models*. Springer, 1983; P. A. Ioannou, J. Sun. *Robust Adaptive Control*. Prentice Hall, 1996 (the σ-modification — a leakage term in the adaptation law — and the ultimate boundedness of the parameters under bounded disturbances; the classical drift of the unmodified law under a decaying regressor and a bounded disturbance; from memory). M. Bodson, S. Sastry, B. D. O. Anderson, I. Mareels, R. R. Bitmead. *Nonlinear averaging theorems, and the determination of parameter convergence rates in adaptive control*. Systems & Control Letters 7(3), 1986 (two-time-scale averaging of adaptive systems); B. D. O. Anderson, R. R. Bitmead, C. R. Johnson, P. V. Kokotović, R. L. Kosut, I. Mareels, L. Praly, B. D. Riedle. *Stability of Adaptive Systems: Passivity and Averaging Analysis*. MIT Press, 1986; S. Sastry, M. Bodson. *Adaptive Control: Stability, Convergence and Robustness*. Prentice Hall, 1989 (the averaged system with the filtered regressor's cross-correlation; from memory). B. Gao, L. Pavel. *Discounted mirror descent dynamics in concave games*. CDC 2019, 5942–5947; A. R. Romano, L. Pavel. *Continuous-time discounted mirror-descent dynamics in monotone concave games*. IEEE TAC 66(11):5451–5458, 2021, arXiv:1912.03460 (a leak in the dual coordinates of a continuous-time mirror flow: its rest points are the equilibria perturbed by the mirror map's regulariser, unique for monotone games, converging as the leak vanishes, with Fenchel-coupling Lyapunov functions). P. Coucheney, B. Gaujal, P. Mertikopoulos. *Penalty-regulated dynamics and robust learning procedures in games*. Math. Oper. Res. 40(3):611–633, 2015 (the replicator's entropic case: exponentially discounted payoffs with a smooth best response). J. Pérolat et al. *From Poincaré recurrence to convergence in imperfect information games: finding equilibrium via regularization*. ICML 2021; K. Abe et al. *Adaptively perturbed mirror descent for learning in games*. ICML 2024 (regularisation toward a reference with a decreasing Lyapunov function; from this round's referee's search). N. Azizan, S. Lale, B. Hassibi. *Explicit regularization via regularizer mirror descent*. arXiv:2202.10788 (explicit regularisation through mirror descent). F. Alvarez, J. Bolte, O. Brahic (SICON 2004, [S74]); H. Attouch (1996, from memory): Tikhonov (viscosity) selection as the regularisation vanishes. J. Duchi, S. Shalev-Shwartz, Y. Singer, A. Tewari. *Composite objective mirror descent*. COLT 2010 (a proximal method — not a dual leak); L. Xiao. *Dual averaging methods for regularized stochastic learning and online optimization*. JMLR 11:2543–2596, 2010. R. Tibshirani. *Regression shrinkage and selection via the lasso*. J. R. Stat. Soc. B 58(1):267–288, 1996; P. Zhao, B. Yu. *On model selection consistency of Lasso*. JMLR 7:2541–2563, 2006; M. J. Wainwright. *Sharp thresholds for high-dimensional and noisy sparsity recovery using ℓ₁-constrained quadratic programming (Lasso)*. IEEE Trans. Inf. Theory 55(5):2183–2202, 2009 (the irrepresentable condition; from memory). Mirror-descent adaptation laws with Bregman Lyapunov functions (Boffi–Slotine, [S72]). The precedents of Proposition D‴.46 (search-verified for their bibliographic details and abstracts where stated): a leak in a mirror flow's dual coordinates, its perturbed rest point, its Bregman Lyapunov function and the vanishing-leak selection are established; what is specific here is their realisation on the paradigm's (α, β) law with Q kept, its bound in the adaptive loop, and the starting masses' LASSO weights.
- [S74] J. Kivinen, M. K. Warmuth. *Exponentiated gradient versus gradient descent for linear predictors*. Inform. Comput. 132(1):1–63, 1997 (multiplicative updates as mirror descent with the relative entropy; their sign-constrained and ± pair (EG±) forms; for the free law's branches, Proposition D‴.47, also Hessian–Riemannian flows on constrained domains — Alvarez–Bolte–Brahic —, nonnegative least squares — Lawson–Hanson — and the hyperbolic entropy of u² − v² parametrisations — Ghai–Hazan–Singer, Woodworth et al.).
- [S75] S. Pesme, N. Flammarion. *Saddle-to-saddle dynamics in diagonal linear networks*. NeurIPS 2023, arXiv:2304.00488 (at vanishing initialisation the time-rescaled gradient flow of a depth-2 diagonal network follows a piecewise-constant jump process ending at the minimum-ℓ¹ interpolator, its recursion the adaptive inverse scale space method of M. Burger, M. Möller, M. Benning, S. Osher, Math. Comp. 82:269–299, 2013; for the rich transit, Proposition D‴.48, also Berthier's incremental learning, JMLR 2023, and Berthier–Pillaud-Vivien's inverse scale space limits of mirror flows, 2026).
- [S49] R. T. Rockafellar. *Monotone operators and the proximal point algorithm*. SIAM J. Control Optim. 14(5):877–898, 1976 (the proximal point method: a quadratic penalty toward a centre moved to the last iterate — iterated Tikhonov regularisation; the leak toward a moving centre is its continuous-time form; for a maximal monotone operator the iterates are bounded iff a zero exists, the precedent of the adaptive prior's lost a-priori bound). M. Hanke, C. W. Groetsch. *Nonstationary iterated Tikhonov regularization*. J. Optim. Theory Appl. 98(1):37–53, 1998 (iterated Tikhonov regularisation: on a quadratic, each step contracts the error in the eigendirection of curvature h by σ/(h + σ), the complement 1 − h/(h + σ) of the Tikhonov filter factor h/(h + σ), i.e. the residual factor; the bound h/(h + σ)² ≤ 1/(4σ) is ridge regression's AM–GM bound on the filter factors). S. Zhang, A. Choromanska, Y. LeCun. *Deep learning with elastic averaging SGD*. NeurIPS 2015 (workers pulled elastically toward a centre variable that moves toward their average: the coupled pair dynamics of the class means and the prior's atoms; its stability analysis of the worker–centre map on a quadratic, in the learning and the moving rate, is the discrete-time form of the pair matrix). S. J. Nowlan, G. E. Hinton. *Simplifying neural networks by soft weight-sharing*. Neural Computation 4(4):473–493, 1992; K. Ullrich, E. Meeds, M. Welling. *Soft weight-sharing for neural network compression*. ICLR 2017 (a mixture prior on the weights whose centres are learned with them; its hard limit is a K-means quantiser of the weights). The precedents of Theorem D‴.23. (Details from memory; not re-verified here.)

## Claim-by-claim verdicts

### A. The adaptive law is a Wasserstein gradient flow of a linear potential

- **Closest prior.**
  - Particle gradient flows coincide with Wasserstein gradient flows from atomic data: [29] Prop. B.1 and [28]. The potential-energy flow transports particles along −∇V: [71] §4.5.
  - RBF units as Wasserstein particles in confined training: [41]. Projection as a support constraint: [46] Prop. 3.1, in a control venue.
  - Time-dependent linear potential with confinement and exogenous data: [38], [48].
  - Adaptive laws read as gradient or Onsager flows of the law of one parameter vector: [25], [17], [18]. An ensemble of estimate particles in Wasserstein space: [26], whose full text could not be read.
- **Distinction.** No source states that the classical centre/width laws with 1/m gains are the characteristics of the continuity equation for a μ-linear potential whose coefficient e(t) comes from plant feedback.
- **Strongest objection.** For a μ-linear potential the "gradient flow" is the continuity equation with velocity −∇V. Its characteristics are the particle ODEs by definition. The identification is therefore essentially tautological, and its ingredients are published.
- **Verdict: partially anticipated.** The contribution is interpretive. It is a change of viewpoint that licenses Theorems B–H, not a theorem with new mathematics. PARADIGM.md says so.

### B. Optimal-assignment Lyapunov function V_OT = ½‖e‖² + ½W²

- **Closest prior.**
  - The minimum over permutations of matched parameter error equals W₂ between neuron measures: [62], [63], [64], [66].
  - Labelled Taylor remainders with explicit curvature constants in neuro-adaptive control: [2], [19], [20], [21], [23].
  - Label-free local analyses with width-independent constants in machine learning: [49], [50], [54].
  - OT costs inside Lyapunov arguments for swarms, where the measure is the plant: [88], [89].
- **Distinction.** No closed-loop Lyapunov function couples the tracking error to W²(μ_t, μ*). None handles assignment switching through Dini derivatives. None has the two-regime remainder with a drift term V_*W.
- **Strongest objection.** The labelled inequality of [20] holds for every relabelling, so minimising over permutations gives the matched bound almost for free. The derivative of W² along transport flows is standard (AGS Thm 8.4.7). m-uniformity comes from the normalisation, not from OT; the round-1 review confirmed this.
- **Verdict: partially anticipated.** The packaging is new; every ingredient is known. The genuine content is uniformity over labellings and across assignment switches, which fixed labellings cannot provide.

### C. OT σ-modification

- **Closest prior.**
  - Labelled leakage: [3], [4], [5].
  - Offline permutation-aware regularisation toward a reference network: [61], [64], [65], [67], [68].
  - Proximal OT to the previous iterate: [45].
  - Curvature obstruction and assignment switching: [71] §4.4 two-atom example, [72], [73] (AGS Example 9.1.5, Lemma 9.2.1, Thm 7.3.2).
  - Local convexification near a target measure: [76].
- **Distinction.** No source found has online-recomputed assignment leakage inside a Lyapunov-based adaptive law, a strict σ-modification inequality with bias W(μ*, μ⁰), an assignment-gap condition, or a reshuffling result for the labelled leak (PARADIGM.md C(v)).
- **Strongest objection.** C(i) and the obstruction are textbook. The gap condition is LP sensitivity of the optimal assignment. Three probes in this direction were blocked.
- **Verdict: novel as a control construction, medium confidence.** C(i) and Remark C.2 are cited as known.

### D. Mean-field closed loop and Dobrushin estimate

- **Closest prior.**
  - ODE coupled to a continuity equation, with a Dobrushin-type estimate and mean-field limit: [77] Prop. 3.2 and Thm 3.3, [78] Thm 4.1, [79], [80].
  - Dobrushin's estimate and exact empirical solutions: [82], [83].
  - Particle systems that act as a controller: [85]. Closed-loop mean-field learners: [43], [44].
  - Uniform-in-time propagation of chaos: [36], [37], [S5].
  - Limit-interchange caveat: [81], [86].
- **Distinction.** None of these couples a tracking plant to the continuity equation of the controller's own parameter measure.
- **Strongest objection.** D is an instance of [77]–[79], with W₂ in place of W₁. Its constant e^{ΛT} is bounded-horizon only.
- **Verdict: anticipated in technique.** It is stated as an application.

### E. Label-free local parameter convergence under tangential excitation

- **Closest prior.**
  - Static analogues: [49] (non-degeneracy, local exponential convergence, label-free), [50], [52], [54] (lifted Gaussian weight–centre–width space), [57].
  - Control-side Jacobian-excitation convergence for a labelled vector: [12], [20], [22].
  - Identifiability only up to permutation: [58], [59], [60].
- **Distinction.** Here PE is a condition on the time-varying tangent regressor along the closed-loop reference. The result is stated for the matched error jointly with e, and the labelled error provably does not converge.
- **Strongest objection.** This is [49] Assumption 5 transplanted into a closed loop, with PE replacing non-degeneracy.
- **Verdict: partially anticipated.**

### F. Transport cost of adaptation (Wasserstein length ≤ λ_G^{1/2} L₁ ∫‖e‖)

- **Prior.** The length bound for absolutely continuous curves is AGS Thm 8.3.1. No sweep targeted F.
- **Verdict: elementary.** It is presented as a corollary with a control interpretation, not as a contribution.

### G. Variable-projection reduction of the slow geometry

- **Closest prior.** Two-timescale (variable-projection) training of two-layer networks: [S1] (global convergence in a univariate setting) and [S2] (mean-field variable projection with rates). Variable projection itself: [S7]. Averaging and singular perturbation in adaptive control: [S9].
- **Distinction (corrected in round 2).**
  - What the closed loop supplies: the error-channel and averaging scales, and the data distribution (the occupation measure of the reference). The loss is the occupation-measure loss, and the controller learns on its own excitation.
  - What the designer supplies: the separation between weights and geometry is a gain choice, δ ≪ g_a·min(m₀/‖Γ^{1/2}∇²L*Γ^{1/2}‖, 1/‖Da*Γ^{1/2}‖²) with (g_c, g_s) = δ(γ_c, γ_s), exactly as in [S1], [S2]. Version 2's claim that the loop produces this separation was wrong.
  - What is new: the leading-order slow flow, valid for every slope S, is the gradient flow of the profile loss in the metric induced on the slow manifold (PARADIGM.md G(b′)). It reduces to plain variable projection only when the slope S = (δ/g_a)‖Da*Γ^{1/2}‖² is small. `exp_slow_manifold_slope.py` gives a reproducible counterexample to plain VarPro under the rate condition alone.
  - L* is the *unregularised* profile loss: the λ → 0 limit of λℒ^λ of [S2] (time rescaled). It is not [S2]'s λ → 0 functional ℒ⁰, which is min-norm interpolation with an ultra-fast-diffusion flow. Version 2's statement that "the geometry marginal follows the Wasserstein gradient flow of the variable-projection functional of [S2]" was false: as a functional of the geometry marginal, L* depends only on the support, so it is not W-continuous and has no mean-field form.
- **Strongest objection.** Each reduction step is a standard theorem on its own time scale, and their composition on the O(k_e/(mδ)) horizon is formal. The induced-metric correction is a standard geometric singular-perturbation computation.
- **Verdict: a fixed-m, closed-loop analogue of the two-timescale/VarPro training of [S1], [S2].** The induced-metric slow flow and the slope condition (H2c), with its counterexample, are the content added here. My searches were standard, not an exhaustive sweep, so confidence is medium-low.

### H. Langevin adaptive control: the averaged mean-field loop is mean-field Langevin dynamics

- **Closest prior.**
  - [25] (Patil, Oct 2026): projected adaptation as an Onsager flow, σ-modification as confinement, and reflected Langevin as a constrained free-energy flow, for linearly parametrised models and the law of one parameter vector. The functional is linear in the law, and no global result is claimed.
  - [24] (LyLA-Therm 2025): Langevin parameter updates in a labelled deep-network controller, with ultimate boundedness in probability.
  - [38]: online mean-field learning with a time-varying linear potential and uniform-in-time propagation of chaos; the data are exogenous and there is no plant.
  - Mean-field Langevin theory: [S3], [S4], [S5].
- **Distinction.** The measure here is the empirical measure of kernels with nonlinear mean-field coupling, and the functional is convex but not linear in μ. The loss is the occupation-measure loss of the closed loop. Three control devices map onto the hypotheses of [S3]: saturated lifted features give bounded features, σ-modification in the G⁻¹ metric gives Bakry–Émery confinement, and parameter dither gives entropy. The log-Sobolev constant comes out in control quantities (σ, τ, k_e, a_max, F_ρ).
- **Strongest objection.** The global convergence is [S3]/[S4] applied verbatim. The step from the real closed loop to the averaged system is formal (PARADIGM.md Remark H.0), and the statement for the real loop is Conjecture H.6, restated in round 2 in the scaling σ = σ₀/k_e, τ = τ₀/k_e.
- **Verdict: novel as a bridge.** None of the searches found it in adaptive control.

### H.7 Structure of the exactly averaged closed loop: semi-gradient and orbit-area circulation

- **Closest prior.**
  - Second-order averaging of periodically forced systems [S9] (Sanders–Verhulst–Murdock).
  - Performative prediction [S14]: gradient steps that hold the model-dependent data distribution fixed. This is the structural analogue of H.7(c).
  - Semi-gradient TD [S11], which holds the bootstrapped target fixed, and its mean-field analysis with two-layer networks [S10].
  - Non-reciprocal interactions and the failure of free-energy descriptions in many-body systems [S12].
- **Distinction.** No source found states three things together.
  - The first correction to the averaged adaptive closed loop is a semi-gradient of the closed-loop loss on the closed-loop state distribution, and its performative defect is non-reciprocal.
  - It includes an antisymmetric kernel–kernel interaction for n ≥ 2. Its kernel is an integral of Σ_k det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ')] over the region that the reference orbit encloses (Stokes).
  - So at finite feedback gain the averaged adaptive dynamics is, in general, not the W_{G⁻¹}-gradient flow of any free energy. This is checked for n = 1 and n = 2; the state-shift (semi-gradient) part needs no enclosed area.
- **Search record (round 3).** Three targeted searches found no combined precedent:
  - "non-reciprocal mean-field interaction adaptive control averaging closed-loop neural network parameters not gradient flow";
  - "adaptive control as semi-gradient learning closed-loop state distribution depends on parameters two-timescale averaging non-gradient correction";
  - "averaging adaptive systems periodic reference second-order correction antisymmetric term area enclosed orbit geometric phase parameter drift".

  A fourth search confirmed [S10].
- **Strongest objection.** The expansion is routine second-order averaging at a frozen measure. The term from the intra-period parameter oscillation (PARADIGM.md Remark H.7a) is the classical second-order averaging term. The semi-gradient reading is an interpretation, and the circulation kernel is Stokes' theorem.
- **What is specific here.**
  - The closed-loop consequence: the learned kernels depend on the direction in which the reference orbit is traversed, at order 1/k_e, although the occupation measure, and with it the averaged mean-field system, is unchanged.
  - The first-order theory predicts this effect quantitatively in the real closed loop: to 3% in size, with cosine 0.994, at k_e = 40.
- **Verdict: no precedent found; the ingredients are known.** Confidence is low to medium, because these were three targeted searches, not a sweep.

### E.0 Explicit strict Lyapunov function for the skew-symmetric mixed error model

- **Closest prior.** Strict Lyapunov functions for persistently excited systems built from a backward–forward integral
  of the PE Gramian: [S15], and [S16] for adaptively controlled systems. Explicit UES bounds for the identifier
  ζ̇ = −ΨᵀΨζ: [S17]. The qualitative statement for the mixed class is Morgan–Narendra [S8], cited by version 3.
- **Distinction.** The lemma gives, for the mixed class ė = −ke − Ψζ, ζ̇ = Ψᵀe, a closed-form rate and overshoot in
  (k, α, T, Ψ_M, Ψ_D), with the cancellation of the ‖Ψζ‖² terms between the Gramian cross term and the eᵀΨζ cross term
  as its only idea, and it quantifies its own conservativeness (Remark E.0a: a factor 3κ(1 + Ψ_D/(kΨ_M))² in the
  slow-adaptation regime, κ the condition number of the averaged Gramian).
- **Verdict: a known construction made explicit for this class.** Its role is rigour, not novelty: Theorem E no longer
  depends on an unverified citation. The ω⁻² law of Remark E.1 now has a proof.

### E.0′ Sharp explicit rate in the slow-adaptation regime

- **Closest prior.** Composite (quadratic-type) Lyapunov functions for singularly perturbed systems [S21], where a fast
  variable is slaved to a slow manifold and the interconnection terms are absorbed with sup-norm constants; averaging
  Lyapunov functions for slowly time-varying systems [S22], where a zero-mean antiderivative of the fluctuation of the
  vector field is added as a correction term. Lemma E.0′ combines the two for the mixed class (LTV): a second-order
  slaving of the error (e₂ = −εΨζ + ε²Ψ̇ζ − ε³ΨMζ) and two averaging corrections (P₁, P₂) of the slow quadratic form.
- **Distinction.** The interconnection and remainder terms are bounded not in sup norm but per direction against the
  average excitation ‖v‖_{M̄} — the *peak factors* ρ, ρ_D, ρ_DD — which is what makes the rate λ_min⟨ΨᵀΨ⟩/k (the exact
  averaged rate) rather than quadratic in λ_min; the identity N₂ = ½⟨P₁M + MP₁⟩ = 0 shows that the second-order averaged
  correction to the symmetric part vanishes for this class, which explains the measured k⁻² relative deviation of the
  exact rate from λ_min/k. No source stating an explicit rate for the mixed PE class that is sharp in the slow limit, or
  using peak factors of the regressor's derivatives as the threshold data, was found (the averaging and PE literatures
  were not swept exhaustively).
- **Verdict: a combination of two classical constructions with a new choice of norms; the sharp explicit rate for the
  mixed class appears new.** Its scope is the slow regime above an explicit threshold; the fast regime stays with
  Lemma E.0.

### E.0″ The slow regime for non-periodic PE regressors

- **Claim.** Lemma E.0′'s construction with the period mean replaced by a forward exponential average over a window T,
  whose derivative identity TṀ̄ = M̄ − M gives the antiderivatives P₁ = −2T(M̄ − M̂) and S = −T(F̄ − F̂) in closed form; the
  periodic cancellation N₂ = 0 is replaced by a residual N̄₂ whose price δ₂ = εn₂/λ is O(1/k); the rate is (1 − δ − δ₂)λ_T/k with
  λ_T the exponentially weighted PE level of the window.
- **Prior.** Strict Lyapunov functions for PE systems from windowed or exponentially weighted integrals of the regressor are
  the Mazenc–Malisoff construction [S15], [S16]; averaged (auxiliary) systems with sliding or weighted windows appear in the
  stability theory of slowly varying systems [S22]; the exponentially weighted PE condition is equivalent to classical PE
  (elementary). Verdict: the ingredients are standard; what is new is the same as for Lemma E.0′ — the rate is the averaged
  rate λ_T/k up to 1 − δ − δ₂, with all constants computable from the regressor — now without periodicity, which is what the
  Lorenz-tracking experiments of Part 1 need.

### I. Reciprocity criterion for the exactly averaged loop

- **Closest prior.** The coboundary (cyclic-sum) criterion for reciprocal pairwise interactions is standard in the
  non-reciprocal many-body literature [S12]; the injectivity of the Gaussian convolution (Weierstrass transform) and
  Stokes' theorem with winding numbers are classical analysis.
- **Distinction.** No source found states that for a feature family containing the zero feature the criterion collapses
  to the vanishing of the antisymmetric kernel, nor the characterisations: for n = 1 the averaged adaptive loop is
  reciprocal at order k_e⁻² iff the residual vanishes on the orbit; for n ≥ 2 at perfect fit iff the orbit's 1-current
  vanishes and the residual's Jacobian is symmetric along the orbit; and no open set of measures is reciprocal. The
  decomposition of the performative defect into the curl of the residual and a Lie bracket of the feature fields was
  not found either.
- **Verdict: new, with classical tools; the content is specific to the closed-loop learner.** The result settles
  README open problem 3 of version 3 (a general non-reciprocity theorem beyond checked instances) for n = 1 and for the
  perfect-fit case, and reduces the general case to an explicit bilinear form.

### C′. The barycentric OT leak toward a prior with a different number of atoms

- **Closest prior.** The pure law ϑ̇_i = −σ_m(ϑ_i − T̄_i) is the continuous-time Lloyd descent of coverage control
  [S18] with the Voronoi cell replaced by the cell of an equal-mass optimal plan, that is, capacity-constrained
  quantisation through optimal transport [S19]; OT-based coverage control is [89]; Lloyd's algorithm and centroidal
  Voronoi tessellations are [S20]. The repository's earlier pointer [14] is an application of CVT to adaptive
  estimation, not the quantisation theory itself.
- **Distinction.** The law is used as a σ-modification leak inside an adaptive controller, and the analysis is the
  closed-loop one: the dissipation inequality up to the within-cell variance, the induced-coupling bound √C_bar ≤ b + 2η,
  the Lyapunov inequality with the bias b_bar (the ideal's centroid defect), and the transfer reading (anchoring m kernels
  to a learned measure of m₀ ≠ m kernels). The ω-limit statement for the pure leak is proved by LP stability for
  generalised solutions (round 5, C1), not by the LaSalle argument of [S18], whose Voronoi cells vary continuously.
- **Verdict: the construction is known in coverage control and quantisation; new as a leak in adaptive control, with
  the induced-coupling and Lyapunov analysis.** The §0 and README rows say so.

### D′. Propagation of chaos for the dithered closed loop

- **Prior.** Sznitman's synchronous coupling for McKean–Vlasov systems with Lipschitz drifts, whose standard output is
  the m^{−1/2} rate for the particles and for averages of bounded functions — a dimension-free exponent — with the
  dimension entering only through the empirical-measure term (Fournier–Guillin, [S13]).
- **Distinction.** The plant is in the loop: the mean field enters through a deterministic state x̄ driven by the
  empirical mean of the features, and the constants are control quantities (K_Φ, the Lipschitz bounds of the saturated
  features, k_e enters only through the tube). The per-particle sup bound is the pathwise Grönwall estimate (round 5, D1).
- **Verdict: the technique and the dimension-free exponent are standard; the closed-loop statement with its constants
  is new in this setting.** Uniformity in time is open (PARADIGM.md §8, branch 3).
- **Theorem D″ (uniform in time in the leak-dominated regime).** The contraction argument with a confining potential is the
  standard route to uniform-in-time propagation of chaos for McKean–Vlasov systems (e.g. the small-interaction regime of
  Malrieu 2001/2003, and [S5] for mean-field Langevin dynamics); the closed-loop version adds the plant as a second
  contracting channel and a small-gain condition. Its certified regime is the one in which the leak stops the learning
  (round 6, DPP-2, DPP-3): the result is a statement about where the method stops. Verdict: standard technique, new
    object; the learning regime, where the long-horizon numerics show the coupled distance growing, is open.

### D‴. Uniform-in-time propagation of chaos in the learning regime (the averaged slow loop)

- **Claim.** For the G-preconditioned gradient flow of J_m = mL̄(μ_Θ) + leak (the averaged slow loop at τ = 0), synchronously
  coupled to the characteristic flow of its mean-field equation, the mean-square coupling distance D obeys ½Ḋ ≤ −(σ₀ − κ_R)D +
  (L₂^G)²D²/16 + L_ϑ^G√D·Y (+ L₂^GYD), the interaction term being −‖f̂_Θ − f̂_Θ̄‖² at linear order; under the strict second-order
  optimality condition along the mean-field flow (a neighbourhood of the support, margin κ, from time s₁), sup_s E D ≤
  4e^{2Λ₁}(L_ϑ^G)²K_Φ²/(κ²m) + an exponentially small tail (PARADIGM.md §5d).
- **Prior.** Uniform-in-time propagation of chaos by synchronous coupling under *global* convexity of the confinement and
  convexity of the interaction is Malrieu's [S37]; non-convex confinements with *small* Lipschitz interactions are handled by
  reflection coupling [S38]; flat-convex energies with dither by the entropy method [S5]. The Hessian decomposition
  "Gram + residual-weighted feature curvature" is the standard tool of the interpolation / PL-condition literature [S39] and
  of lazy-training analyses, where the second term vanishes at interpolation.
- **Distinction.** (i) The interaction of the fitted network is neither small nor convex in the particle coordinates; its
  sign comes from the flat convexity of L̄ (f̂ linear in the measure), which makes −‖f̂_Θ − f̂_Θ̄‖² the leading interaction
  term under the coupling — the finite-m, two-measure version of [S39]'s decomposition. (ii) The convexity needed is local:
  σ₀I − G^{1/2}R[μ_s]G^{1/2} ⪰ κI on a neighbourhood of the support of the mean-field flow, i.e. the strict second-order
  optimality condition of the fitted configuration, evaluated at the actual residual — not a global convexity (which would
  force the collapsed configuration to be stable and the network not to learn, Theorem D″). (iii) The window κ_R[μ_∞] <
  σ₀ < σ_crit between the collapsed instability (Proposition K) and the fitted stability is a design statement about
  the leak that has no counterpart in the cited work. (iv) The limits are stated with the result: a connected initial
  support that splits between basins violates the hypothesis for all times (the characteristic flow is a homeomorphism),
  and the coupling then has a floor of pairs in different basins; the real loop at finite k_e is not covered.
- **Proposition D‴.5 (the weight-sensitivity law).** At finite m the stationary coupled distance is the multinomial fluctuation
  of the cluster weights propagated through the implicit function w ↦ Θ*(w): a delta-method statement (the Gaussian limit of a
  smooth function of a multinomial), standard in statistics; what is specific is its object (the shared-particle coupling of a
  particle approximation to a stationary K-atom measure, with the covariance (1/m − 1/M)(diag w̄ − w̄w̄ᵀ)) and its role as the
  sharp form of the theorem's constant (stiffness⁻¹ × feature correlations against κ⁻¹L_ϑ^GK_Φ).
- **Verdict: the technique (synchronous coupling + bootstrap + bounded differences) is standard; the sign structure of the
  interaction and the local, residual-weighted hypothesis — and with them a uniform-in-time statement in a regime where the
  network learns — are new in propagation of chaos for adaptive control; the theorem's regime is asymptotic in m and the
  finite-m law is a delta-method statement, new in its object.** Numerics: `exp_chaos_learning.py`.
- **Proposition D‴.9 (measure-valued priors).** Leaks toward a single prior (σ-modification, Ioannou–Kokotović) are the standard of
  robust adaptive control; regularisation toward the nearest of several centres is standard in machine learning (mixture priors,
  prototype regularisers) and the paradigm's own Theorem C′ leaks toward a prior *measure* through an optimal plan. The nearest-atom
  leak is the simplest measure-valued prior; what is new is its role in the window of Theorem D‴ — the prior's placement sets the
  stationary residual, hence the lower end of the window and the relaxation rate — and the finding that removing the bias raises
  the coupling's floors (on the instance the weight sensitivity rises with the amplitudes): a bias–variance trade-off stated for the
  particle system.
  Verdict: standard object, new design statement; the function-level sensitivity law is the delta method once more.
- **Proposition D‴.12 (the mass constraint).** Capacity-constrained and equal-mass quantisation through optimal plans, and semi-discrete
  optimal transport with its Laguerre cells, are classical ([S19], [S42]); the integrality of the transportation polytope and the forest
  structure of its vertices are textbook ([S40]); c-cyclical monotonicity is the standard optimality criterion ([S41]). What is new is the
  observation that the finite-m floor of the mean-field limit of an adaptive law is carried by a conservation law of its leak (the cluster
  masses under the nearest-atom leak), and that a design of the leak — the mass constraint of the OT plan — breaks it: on the commensurate
  sizes the stationary measure of the particle system coincides with the mean-field one, and the finite-m fluctuation of the plan is counted
  (at most K₀ − 1 split particles) rather than estimated. Propagation of chaos with an exactly m-independent stationary law is not in the
  literature we know of; the mean-field form of the barycentric leak for an atomic prior (the Laguerre-cell leak) was §8's open item.
  Numerics: `exp_ot_leak.py`.
  Verdict: textbook ingredients, new mechanism and new object (a conservation law of the leak as the source of the finite-size floor, broken by the plan's mass constraint).
- **Proposition D‴.17 (convergence of the OT dynamics).** The Kurdyka–Łojasiewicz route to the convergence of (sub)gradient trajectories of
  definable functions is standard ([S43], [S44]); its application to an adaptive law whose energy is a semiconcave transport cost plus a smooth
  loss is a routine verification (the saturated Gaussian features are exponential-analytic, the plan is a minimum over finitely many vertices).
  Verdict: imported tool, routine verification; the content is that the mass-constrained law converges without a rate, which separates the
  qualitative question (settled) from Theorem D‴'s quantitative one (open for the OT dynamics).
- **Proposition D‴.20 (the entropic OT leak).** Entropic OT, Sinkhorn's algorithm and barycentric mappings of entropic plans are standard
  ([S45]), and so is the interpolation of the entropic cost between optimal transport and a kernel (MMD-type) discrepancy as ε runs from 0
  to ∞ (Sinkhorn divergences, [S45]); size-constrained soft assignment by Sinkhorn (equipartition) is related prior art for a soft mass
  constraint. A critical temperature equal to twice the largest eigenvalue of a conditional covariance is deterministic annealing's ([S46]),
  with the same roles — the particles as the codevectors, the prior's atoms as the data; here the curvature of the fit enters, and the
  threshold is the largest generalised eigenvalue of the conditional covariance relative to I − R^G/σ₀. What is new is the soft mass
  constraint as an adaptive law's leak — a smooth gradient flow between the OT leak and the σ-modification toward the prior's barycentre —
  in which the conditional covariance of a collapsed class is an effect of its mass excess, and its consequence for propagation of chaos:
  the masses equalise only down to a critical excess, so the finite-size floor returns at large m for every positive regularisation.
  Verdict: standard ingredients, new object and new consequence (on the commensurate sizes only the hard constraint is floor-free for every m).
- **Theorem D‴.21 (a measure-level basin for the OT leak).** The Polyak–Łojasiewicz inequality and its linear convergence are classical ([S47]);
  in Wasserstein space Łojasiewicz inequalities are known for λ-convex and displacement-convex functionals, and PL inequalities through
  optimal maps for some that are not geodesically convex — barycentre functionals, the distance of uniform point clouds to a density ([S48]);
  the leak alone, ½W²(·, ν), satisfies one with constant 2 (its metric slope is W(·, ν); AGS §10.4). The theorem's plan segments are the
  generalised geodesic from μ* to μ_Θ based at the prior (AGS §9.2), along which ½W²(·, ν) is 1-convex, and turning convexity along such a
  curve into a slope bound is AGS's device. The adaptive law's energy is semiconvex along generalised geodesics only with a negative constant
  (the loss's curvature at the actual residual), so a global argument is unavailable. What is new is the PL inequality along that geodesic
  for this energy — the loss linearised exactly in the fitted function with the curvature taken at the stationary residual, where the fixed
  counts of the mass constraint make the geodesic a particle interpolation and the fitted function's increment the average of the
  per-particle increments: its constant is set by the residual-weighted curvature along the segments and does not depend on m (sharp, 2σ₀,
  at a matching prior); an m-independent energy basin, and uniqueness of the global minimum among empirical measures under a curvature
  hypothesis at the stationary residual; and the converse for the nearest-atom leak, whose energy has a spurious minimum at every
  admissible count vector — the energy-landscape form of the finite-size floor. Numerics: `exp_ot_basin.py`.
  Verdict: classical inequality and geometry (the generalised geodesic), new mechanism for a new energy (the counts as the source of
  uniformity in m).
- **Theorem D‴.23 (adaptive priors).** A quadratic penalty toward a centre that is moved toward the iterate is the proximal point method
  (iterated Tikhonov regularisation) and, in continuous time with a centre drawn to the mean of several coupled variables, elastic
  averaging ([S49]); a learned mixture prior on the weights is soft weight-sharing, whose hard limit is a K-means regulariser ([S49]),
  and the capacity-constrained form of that limit is [S19]'s. The per-mode formulas are the classical quadratic analyses: σ₀/(h_j + σ₀),
  the complement of the Tikhonov filter factor h_j/(h_j + σ₀) (the residual factor), iterated Tikhonov's contraction by it per step (Hanke–Groetsch, [S49]; in continuous time the gradient
  flow of the Moreau envelope, which is the slow-prior limit ηh_j/(h_j + σ₀) of the rate), and elastic averaging's worker–centre analysis on
  a quadratic, whose continuous-time form is the pair matrix [[−(h + σ₀), σ₀], [η, −η]]; empirical Bayes for the mean of a Gaussian prior
  is the same proximal iteration. What is new is their transfer to the mass-constrained particle system of the OT leak — the prior's atoms
  learned, its masses kept: H_L identified as the measure-space excitation spectrum, the separations decoupled at σ₀, and hence the
  trade-off between the fixed prior's bias σ₀/(h_j + σ₀) and the adaptive prior's rate ηh_j/(h_j + σ₀ + η) < λ_j^− < min{h_j, η} uniformly
  in m — the rate half of the classical trade-off between a leak's robustness and persistent excitation, whose robustness half (the
  leak's a-priori bound) the adaptive prior gives up (Proposition D‴.27 restores it by projecting the prior); and the transit's drag, with the gate that starts the adaptation from the fixed-prior
  stationary point. The limits of the bias removal are largely the precedents' own properties — the lost a-priori bound is the proximal
  point method's behaviour without a minimiser ([S49]), the K₀-cluster bias the quantisation error of soft weight-sharing's hard limit —
  and specific here only in their particle, mass-constrained form and in the run-off from a fit that the network attains exactly. Numerics:
  `exp_adaptive_prior.py`. Verdict: standard ingredients and per-mode formulas (proximal point and iterated Tikhonov, elastic averaging,
  soft weight-sharing), a new mechanism (their transfer to the mass-constrained particle system, uniform in m) and a new failure mode with
  its remedy (the transit's drag and the gate).
- **Proposition D‴.24 (the floor is the leak's bias).** The m^{−1}-size fluctuations of trained mean-field networks are classical ([30],
  [S51]) and Proposition D‴.5/D‴.9's floor law is the delta method; the amplitude–mass degeneracy is the homogeneity of the conic view
  ([49], [55]); semi-relaxed transport (free target masses) reducing to the nearest-centre assignment is elementary. That the fluctuations' asymptotic deviation vanishes at an interpolating
  limit is Chen et al.'s ([S51]); the derivation of F is the chain rule applied to Theorem D‴.23(iii)'s first-order bias at x⁰(q) = 0, and the
  weight h_jσ₀²/(h_j + σ₀)² is the squared complement of the Tikhonov filter factor times the mode's visibility, bounded by ridge regression's
  AM–GM bound. What is new is the identification of the function-level floor of a regularised adaptive law with the weight-derivative of
  the regulariser's bias at the sampled weights — at a matching prior F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰, decomposed by modes with the band-pass
  weight, whose bound is σ₀/4 times the adaptive prior's coupled-distance law — the split of the sample's weight fluctuation between the
  parameters (h_j/(h_j + σ₀)) and the bias (σ₀/(h_j + σ₀)), and a leak toward a learned prior (the OT leak with the masses minimised out) as
  the way to restore the interpolation property under regularisation, which also removes the fixed prior's splitting instability; the
  trade-off is a larger coupled distance of the gauge, the leak's a-priori bound (restored by projecting the prior, Proposition D‴.27)
  and the collective rate. On the instance the law (B) is
  asymptotic (the exact floor 2.8 times it at m = 15), and the samples without an exact representation are not rare at small m (0.36 at
  m = 15). Numerics: `exp_adaptive_floor.py`.
  Verdict: classical ingredients (the delta method, the Gauss–Newton spectrum, semi-relaxed transport, the filter-factor bound, the
  interpolation mechanism), new mechanism (the floor as the bias's slope, mode by mode, and its removal by a learned prior under
  regularisation) and a resolution of Remark D‴.10's bias–variance trade-off at the function level, at the price of the collective rate.
- **Corollary D‴.26 (perfect tracking with an adaptive prior: the Floquet spectrum).** Floquet theory and the averaging of a periodic
  variational equation are classical (first- and second-order averaging; the Magnus expansion, whose second term is a commutator, so that
  for a self-adjoint generator it is skew [S52]); the averaging analysis of adaptive loops, with the rate set by the averaged regressor's
  persistency of excitation, is Kosut–Anderson–Mareels' and Sastry–Bodson's ([S52], [S9]); asymptotic tracking under exact
  parameterisation is the textbook goal of adaptive control, which σ-modification gives up for robustness ([3], [4]), and leaks
  compatible with it are known — the e-modification scales the leak by |e| and vanishes at perfect tracking ([S52]), a centre that follows
  the parameters has unbiased fixed points ([S49]). (a) and (b) are Corollary D‴.14(ii) and Proposition D‴.15 at zero residual, with the
  prior's coordinates added; that an expansion is even when the problem is invariant under a symmetry reversing the small parameter's sign
  is the classical mechanism of the even error expansions of symmetric integrators and of the time symmetry of the Magnus exponent
  ([S56], [S52]). What is specific here is the combination: with the adaptive prior the zero-residual lift exists at every representable
  sample (not only at w = q under the OT leak), its Floquet spectrum is exactly m-independent with collective rates λ_j^±(w)/k_e, and its
  expansion in 1/k_e has an explicit second-order coefficient — the error filter's attenuation against the slaving's speed-up, its sign
  not fixed (positive coefficients appear when the adaptation gain is raised) — which a derivative mismatch ∂_x(f − f̂) ≠ 0 would spoil at
  first order. Numerics: `exp_adaptive_floquet.py`. Verdict: classical tools, a new structural statement (an exact, m-independent spectrum
  at every representable sample, with the explicit second-order coefficient of this loop).
- **Lemma D‴.26a (the exact slow subspace).** Chang's decoupling of linear singularly perturbed systems — the transformation pair that
  block-diagonalises the system into its slow and fast parts — is classical ([S59], [S9]), for time-varying systems by a contraction
  (Yang–Zhu, [S59]), with white-noise inputs and the slow subsystem's transformed noise and covariance (Kokotović–Khalil–O'Reilly Ch. 4,
  Haddad–Kokotović, [S59]); stochastic slow manifolds and their covariance are Berglund and Gentz's ([S59]); the convergence of the Magnus
  series is Moan and Niesen's ([S59]); the evenness of an expansion under a symmetry reversing the small parameter is the mechanism of
  symmetric integrators, and Φ_K(−ε) = Φ_K̃(ε)⁻¹ is their adjoint relation ([S56]). What is specific is the application to the
  perfect-tracking linearisation of the transport-leak loop: the slaving of Corollary D‴.26 becomes an exact block-triangularisation, so the
  slow multipliers are exactly the slow block's; the parity becomes the exact identity eig Ā(−ε) = eig Ā(ε), the slow block at negative ε
  defined by the same (anti-causal) contraction; and the cancellation Kᵀ − (1 + εa)(Kᵀ + εK̇ᵀ) = −ε(K̇ + aK)ᵀ makes the slow noise of
  Proposition D‴.30 εbν, at leading order the high-passed −ε(K̇ + aK)ᵀν (the drift decoupled exactly, the noise shared with the fast part).
  Numerics: `exp_adaptive_floquet.py` Part D. Verdict: classical tools; a rigorous foundation for D‴.26(c) and D‴.30(ii) (the referee
  round 22's audit: fair once the stochastic and time-varying precedents are named).
- **Proposition D‴.27 (a projected adaptive prior).** Projection of the parameter estimates onto a convex set (Pomet–Praly), the exterior
  quadratic penalty and the switching σ-modification are the classical robustness modifications of adaptive laws ([4]), and the switching
  σ-modification keeps the ideal (unbiased) properties of the law when the true parameter lies in its set ([4], [S52]); the a-priori bound
  by a projection inequality and a comparison argument is standard. On collapsed classes at η = ∞ the construction is exactly that law,
  toward the set, with the leak off on it. What is specific here is its transfer to the learned centre of a mass-constrained K₀-atom prior at
  finite η — a projected adaptive law for the prior's atoms (whose velocity jumps at the faces, while the particles' law keeps its smooth
  gradient form) — and its compatibility with Theorem D‴.23's bias removal wherever the target is representable inside the set; at η = ∞
  the plan is Lloyd-consistent rather than jointly optimal (as in Lloyd's algorithm for K-means, whose fixed points are local), so the
  monotone quantity is the selected energy, and Proposition D‴.17's convergence argument does not apply as it stands. Numerics:
  `exp_adaptive_prior.py` Part H. Verdict: classical ingredients and, at η = ∞, a classical law; a new transfer (the learned centre of the
  mass-constrained prior at finite η) that restores the robustness half of the trade-off which the adaptive prior gave up.
- **Proposition D‴.28 (the floor removal is a gauge property).** The amplitude–mass degeneracy is the positive homogeneity that the conic
  (Wasserstein–Fisher–Rao) view of mean-field networks exploits ([49], [55]), and rescaling symmetries of network parameters — which map
  critical points to critical points with the same function while changing their curvature — are classical ([S55]);
  the implicit-function derivation of the fixed prior's sensitivity is the one of Proposition D‴.24(i). What is specific here is that the
  gauge, not exact representability, is the mechanism of the adaptive prior's floor removal — so it survives an approximation error, the
  regime of real plants, at a strict local minimum in the interior — together with the exact scaling R_k(G_wϑ) = (q_k/w_k)R_k(ϑ) of the
  residual curvature at a critical point for a scalar amplitude per kernel, which turns the adaptive prior's splitting margins into a
  closed formula: the scalar-amplitude instance of the reparametrisation dependence of curvature [S55]. Numerics: `exp_gauge_floor.py`.
  Verdict: classical symmetry, a new identification of the mechanism and an exact margin formula for a scalar amplitude.
- **Proposition D‴.29 (second-order plants).** The sliding-variable (filtered-error) construction is the classical route of adaptive control
  to higher-order plants, and Gaussian networks with it are Sanner and Slotine's ([1]); singular-perturbation averaging with a stable fast
  linear subsystem is standard ([S9]); the exact filter (d/dt + k_s)⁻¹ in the position-error prediction is the filtered error itself, and
  the weak excitation of Gaussian networks along an orbit is the classical partial persistence of excitation ([6]–[8]). What is specific
  here is the transfer of the transport leaks of §5d — their finite-size floors, the mass constraint, the adaptive and projected priors and
  the gauge — to such plants through the occupation measure of the regressor's orbit, and the exact fibre of isotropic phase-plane kernels
  on a circular orbit (von Mises restrictions: in general position one zero Gauss–Newton eigenvalue per kernel). The numerics are
  noise-free: the leaks' robustness is not exercised there, and no leak at all fits as well. Numerics: `exp_servo_cogging.py`. Verdict: a
  classical reduction; new applications of §5d's mechanisms to a standard servo problem, whose case rests on the noisy setting still to
  be run.
- **Proposition D‴.30 (measurement noise at perfect tracking).** Measurement noise in adaptive control, the robust modifications against
  bounded disturbances ([4]) and the errors-in-variables bias of estimates from noisy regressors (bias compensation and instrumental
  variables in identification, [S53]) are classical, and so is the remedy of putting the reference in the regressor (the desired
  compensation adaptation law, [S53]); stochastic averaging is Khas'minskii's ([S53]). So are the mechanisms behind the floor and the
  amplification ([S57]): the Lyapunov covariance as the diffusion approximation of a constant-gain algorithm (Benveniste–Métivier–
  Priouret), a misadjustment depending on the input's spectrum only through a trace (Widrow et al.: the η = ∞ floor's tr Q), the leaky
  filter's steady state (Mayyas–Aboulnasr: the fixed prior's factor h_j/(h_j + σ₀)), a common noise not averaging out in the mean-field
  limit (Carmona–Delarue–Lacker: the coherence) and the drift of weakly excited directions under a persistent disturbance (Anderson: the
  H⁻¹ response). What is specific here is their combination in the transport-leak setting: the floor's interpolation
  (h_j + η)/(h_j + η + σ₀) in the prior's mobility (an elementary 2×2 Lyapunov computation), the driving matrix Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩
  of the sensitivity-filtered measurement the law sees (built on the derivative Gram of Corollary D‴.26's filter correction), and the
  coherence of the K₀-atom identification with its noise threshold at the lift (invariance, not attraction; along the noisy trajectory
  the separations' drift grows with the developing errors-in-variables deviation, and the long-time threshold is open). Numerics:
  `exp_measurement_noise.py`.
  Verdict: classical ingredients; a quantified structural statement for the transport leaks under noise.
- **Proposition D‴.31 (the washout prior and the desired regressor).** Washout filters — feedback through a high-pass filter, which leaves
  the open loop's equilibria where they are and therefore cannot oppose a persistent input — are classical in flight control and in
  bifurcation control ([S54]); an adaptive law kept near a low-pass filtered copy of itself is Yucelen and Haddad's low-frequency learning
  ([S58]), which has the washout structure, and a centre moving with its workers' average is elastic averaging's ([S49]); the desired
  compensation adaptation law, whose regressor and compensation use the reference trajectory, is Sadegh and Horowitz's ([S53]; its
  nonlinear damping, which DR* omits, is for global stability), and its regressor's freedom from measurement noise is part of its
  classical motivation. That a quadratic regularisation's stationary response and its bias toward the anchor sum to the identity —
  RH + B = I, per mode the Tikhonov filter factor plus the residual factor, the Kalman-gain split of a regularised estimate — is
  elementary, as are the minimum-norm bound (iii) and the split of a force into the regressor's errors in variables, which bias, and the
  output's, of which only the mean biases. What is specific here is their transfer to the transport leaks: the reading of Theorem D‴.23's
  adaptive prior as a washout σ-modification of the class barycentres (for integral plans; split rows leave a non-washed-out part), uniform
  in m, and the quantification on the instance. The parameters' own errors in variables is the classical second-order bias of a
  constant-gain stochastic approximation from its own fluctuation (Benveniste–Métivier–Priouret; Dieuleveut–Durmus–Bach for constant-step
  SGD; Box's curvature bias of nonlinear least squares: [S57]); what is added is its quantification here — small in absolute terms
  (≤ 3.9·10⁻⁴ per (Eν²)²), −8.7% to +3.3% of DR's and DR*'s biases and +56% of DR's at k_e = B = 100, where DR's noise part nearly
  cancels; its own response falling like k_e⁻² on the one-dimensional reference, where the law's leading term averages out — and the
  observation that the desired regressor does not remove it. Numerics: `exp_eiv_remedy.py`. Verdict: classical ingredients; a structural
  reading of the adaptive prior and a quantified, partial remedy for the transport leaks.
- **Proposition D‴.32 (the cohesion leak).** The leak is the mean-shift step of the parameters' own kernel density (Fukunaga–Hostetler,
  Cheng: mode seeking with a shadow kernel), so its mean-field form is an aggregation equation in the fit's external potential with the fit's
  and the cohesion's interactions — a Wasserstein gradient flow of an interaction energy (Carrillo–McCann–Villani, Bertozzi–Laurent–Rosado);
  particles drifting toward their own kernel-weighted local means are polarized consensus-based dynamics (Bungert–Roith–Wacker), and
  clustering under softmax-weighted attraction is the self-attention dynamics of Geshkovski–Letrouit–Polyanskiy–Rigollet; agents averaging
  within a confidence radius form clusters at least the radius apart (Hegselmann–Krause), which is what the finite range does to the
  classes; a duplicated hidden unit turns a critical point into a line of critical points (Fukumizu–Amari), the twins' curves of (v); a
  neuron split along the top eigenvector of its residual-weighted Hessian is splitting steepest descent (Wu–Wang–Liu); a prior-free
  pairwise attraction whose strength trades the fit against the number of distinct parameters is convex clustering and sum-of-norms
  identification (Hocking–Joulin–Bach–Vert, Ohlsson–Ljung–Boyd), and a cluster splitting when an annealed parameter crosses a value set by
  its covariance is deterministic annealing (Rose), the structure of (iv)'s resolution and of the annealed scale; RBF units added and
  removed by accuracy criteria are self-organising RBF networks (Lian et al.); forgetting only in the excited directions is directional
  forgetting (Kulhavý), and forgetting toward a prior on the unexcited ones stabilized forgetting (Kulhavý–Zarrop); the extra modes of an
  isotropic Gaussian mixture need overlapping components (Carreira-Perpiñán–Williams), which (v) uses ([S60]). What is specific here is the
  combination as a leak in mean-field adaptive control: the exact blocks at clustered configurations with the mass factor σ₀w_k, the
  gauge-scaled splitting margin σ₀w_k − (q_k/w_k)λ⁺_k, the twins' contraction rate, and the real loop's periodic lift with the separation
  exponents σ₀w_k/k_e. The number of kernels is not selected by the leak (Numerics (C), (J), (K)). Numerics: `exp_cohesion_leak.py`.
  Verdict: classical ingredients; a new prior-free leak on the unexcited subspace.
- **Corollary D‴.33 (the desired regressor and the coherence under noise).** That a regressor evaluated on the reference trajectory
  removes the measurement noise from the adaptation's regressor is part of the classical motivation of the desired compensation
  adaptation law (Sadegh and Horowitz, [S53]). The corollary's question — whether identical subsystems driven by a common noise stay
  together — is that of noise-induced synchronisation: identical uncoupled oscillators driven by a common weak noise synchronise when the
  transverse Lyapunov exponent, to second order in the noise, is negative, and a strong noise can make it positive (Teramae–Tanaka,
  Goldobin–Pikovsky; with independent noise, clustering: Nakao–Arai–Kawamura; [S61]). Proposition D‴.30(i)'s threshold and this corollary
  are that exponent for the replicated kernels, the leak their transverse stiffness, the variance term the exponent's own second-order
  correction. That the regressor's term leaves the separations' drift is the one-line consequence of a regressor that does not see the
  measurement. What is specific here: that consequence for the K₀-atom identification — the drift O(1/k_e) like the stiffness, a
  threshold bounded below uniformly in k_e and constant when the noise's bandwidth scales with the gain — and its quantification on the
  instance: the sign at the lift (negative in the centre directions of separated kernels; the whole matrix only for separated
  alternating-sign kernels), its loss along the trajectory as the mean error relaxes with the collective modes, and the stationary point,
  where the parameters' part of the force dominates the mean error. Numerics: `exp_eiv_remedy.py` Part D. Verdict: a classical remedy in a
  classical frame; a new structural consequence for the transport leaks (and for the cohesion leak, whose separations are the softest).
- **Proposition D‴.34 (mass transfer: the replicator weight law).** Transport and creation of mass together are the Wasserstein–Fisher–Rao
  (Hellinger–Kantorovich) geometry of positive measures (Liero–Mielke–Savaré); particles that move and change weight are conic particle
  gradient descent, with global convergence for sparse optimisation on measures under non-degeneracy (Chizat), and a non-local
  birth–death term added to the mean-field dynamics of a one-hidden-layer network keeps the minimiser and speeds the convergence
  (Rotskoff–Jelassi–Bruna–Vanden-Eijnden); the multiplicative weight update is the exponentiated-gradient law (Kivinen–Warmuth) and EM's update of mixture weights (Dempster–Laird–Rubin), adaptation laws of mirror-descent type replace the quadratic parameter error of the Lyapunov function by a Bregman divergence (Boffi–Slotine) — (ii) is their negative-entropy instance —, and the replicator equation is the gradient flow of Shahshahani's metric with the relative entropy as its Lyapunov function at a stable state (Shahshahani, Hofbauer–Sigmund); the Wasserstein–Fisher–Rao distance has two more origins (Kondratyev–Monsaingeon–Vorotnikov, Chizat–Peyré–Schmitzer–Vialard), and gradient descent on the positions and square-root weights of a two-layer network is its particle form (Chizat–Bach) ([S62]). The law here is gradient descent on (√β, Θ) — Euclidean on the positions, not the WFR metric's mass-weighted
  one —, with the cohesion energy of the unweighted cloud. What is specific: the replicator as the mass transfer of the cohesion leak in mean-field adaptive control — the stationary set unchanged at clustered configurations under the finite range (under the Gaussian kernel the replicator selects the masses), the saturation's bound kept without a projection, the saturation trap of the
  transit's masses removed, the representability criterion Σ_k|A_k|/a_max < 1 at fixed positions for bounded trajectories, the zero modes at the lift — and the
  numerics (the five-bump and three-heavy-bump targets found with their kernels; bumps closer than the broad scale and the shedding not
  repaired). Numerics: `exp_mass_transfer.py`. Verdict: classical ingredients; a new remedy for the transit's masses.
- **Proposition D‴.35 (the split as an S_n-symmetric transcritical bifurcation).** Bifurcation from a fully symmetric state of N identical
  coupled units is S_N-equivariant; the standard representation admits a quadratic equivariant, so the generic bifurcation is transcritical,
  with branches of isotropy S_q × S_{N−q} and no stable small branch — the jump of the speciation models of Stewart–Elmhirst–Cohen and
  Golubitsky–Stewart, the normal forms and stability of Elmhirst, the generic theory with quadratic equivariants of Arjevani–Field, who
  motivate it by the permutation-symmetric critical points of two-layer networks; duplicated units give such critical points and saddles
  (Fukumizu–Amari [S60], Brea et al., Şimşek et al.) ([S63]); a cluster splits at a critical parameter set by its covariance in deterministic
  annealing (Rose, [S60]). Order preservation is a consequence of the equivariance, and the escape time is the Riccati solution. What is
  specific here: the identification for the cohesion leak — the quadratic coefficient the residual's third derivative along the splitting
  direction, the cohesion and the mean term entering only at third order or in the common direction —, the margin after a departure (the
  stiffness's drop and the curvature's gauge rise), and the instance's quantification (the extreme particle leaving first, alone or with a
  near-tied neighbour; the mass-free leak delaying the cascade without stopping it). Numerics: `exp_split_bifurcation.py`. Verdict:
  classical bifurcation theory; a new identification for the cohesion leak's over-splitting.
- **Proposition D‴.36 (the gauge under dither).** A diffusion pushed onto a manifold of equilibria by a large drift converges, on the slow time
  scale, to a diffusion on the manifold whose noise is the original projected by the derivative of the flow's limit map (Katzenberger [S64]);
  the noise-driven motion of SGD along the manifold of minimisers of an overparameterised model is the same limit (Blanc et al.; Li–Wang–Arora
  [S64]), and the constructive projection is Parsons–Rogers'. The entropic leak on the weights is the negative-entropy mirror map's
  regulariser ([S62]). What is specific here: the gauge of the amplitude–mass representation (Proposition D‴.28) as such a manifold for the
  replicator of Proposition D‴.34; the projection computed in the flow's own metric (Euclidean–G⁻¹ on the positions, Fisher–Rao on the masses), with
  the two-class closed form; the excursion into the saturation, where the limit map's derivative degenerates; and the entropic leak on the
  masses as a leak on the unexcited subspace, which pins the gauge without bias where the particle fractions are representable. Numerics:
  `exp_gauge_dither.py`. Verdict: a classical limit theorem applied; the replicator's unpinned masses under dither, and their remedy.
- **Proposition D‴.37 (the dither's drift along the gauge).** The second derivative of the limit map and the drift it gives a diffusion
  forced onto a manifold of equilibria are Katzenberger's; for isotropic noise the drift along a manifold of minima is −½ of the gradient
  of the log-pseudo-determinant of the Hessian (Li–Wang–Arora's Corollary 5.1, [S64]; checked by round 27's referee on its summary only);
  the averaging of a slow variable over a fast ergodic one is Khasminskii's, and a slow coordinate in equilibrium with fast ones feels their
  free energy, the potential of mean force ([S65]); that noisy gradient dynamics prefers flat minima is the flat-minima picture of
  Hochreiter–Schmidhuber, and Entropy-SGD makes the local entropy the objective ([S65]). Not seen, named by the referee: Fixman's and
  Bornemann–Schütte's corrections for stiff constraints, the label-noise analyses of Damian–Ma–Lee and HaoChen et al., noisy replicator
  dynamics (Fudenberg–Harris, Foster–Young), Rose's deterministic annealing. What is specific here: the replicator's masses as the slow
  variable on the amplitude–mass gauge, with noise on the positions only, and its free energy in closed form (2Σ log(W_ka′_k)); the
  consequence that dither drives the masses into the saturation — the trap the replicator removed (Proposition D‴.34) returning through the
  noise —; the boundary layer of width O(√γ), where the replicator outruns the saturating amplitude, which slows but does not stop it; and
  the entropic leak's bias κ⁻¹b(w). Numerics: `exp_gauge_drift.py`. Verdict: classical limits (Katzenberger, Li–Wang–Arora, averaging)
  applied; a new mechanism for the replicator's gauge under dither.
- **Proposition D‴.38 (a Fixman correction for the replicator under dither).** That a stiff system's free energy along its soft directions
  carries a log-determinant, and that adding the opposite term makes the constrained dynamics sample as intended, is Fixman's (the
  log-det term the stiff limit's zero-point term, Helfand; disputed in the overdamped setting by Waszkiewicz–Lisicki); the drift
  corrections between stiff bonds and rigid constraints are Hinch's, and the projected diffusions and the mean force on submanifolds are
  Ciccotti–Lelièvre–Vanden-Eijnden's ([S66]); removing a free-energy drift online is the adaptive biasing force method, a gain on it
  well-tempered metadynamics' bias factor, and the max-det maximiser of a log-barrier is classical in convex optimisation ([S66]); the
  log-determinant drift itself is Li–Wang–Arora's for isotropic noise ([S64]) and Proposition D‴.37's for the replicator. What is
  specific here: the correction for the replicator's gauge, applied online to the masses only, in closed form on the masses — a
  log-barrier on the representable masses, the corrected loop a gradient flow — and the distinction between that masses form and the
  form at the current amplitudes, which reads the positions' escapes off the manifold as masses at their saturation; the selection of
  the D-optimal exact fit W♯ with the gain and the confinement exponent δ^{κ_F − 1} are applications; the leak's pinned state predicted
  to second order (the fluctuations through the leak's curvature, which the paired seeds resolve). Numerics: `exp_gauge_fixman.py`.
  Verdict: a classical correction carried to a new setting, with a new design knob.
- **Proposition D‴.39 (dissipativity of the dithered gauge).** That a Langevin dynamics needs a confining potential to have an
  invariant probability law is classical (Mattingly–Stuart–Higham; in learning, Raginsky–Rakhlin–Telgarsky's dissipativity), as are the
  exit times' exponential order (Freidlin–Wentzell; sharp, Bovier–Eckhoff–Gayrard–Klein) and parameter projection with its Lipschitz
  boundary layer (Ioannou–Sun, Pomet–Praly), and the dead zone itself is the Moreau–Yosida envelope of the box's indicator, with which
  MYULA samples on a convex set (Durmus–Moulines–Pereyra; Brosse et al.) — a switching σ-modification in adaptive control ([S67]). What
  is specific here: the cohesion leak's removal of the prior (Proposition D‴.32)
  removes with it the dissipativity the prior leak gave the dithered loop; on the instance the dither's way out of the exact fit's well is
  a monotone climb, without a saddle, to a plateau of approximate fits with saturated amplitudes; the dead zone on the unconstrained
  coordinates replaces the gauge's saturation by a non-degenerate wall, so the entropic drift of Proposition D‴.37 ends at the vertices of
  the representable polytope (at leading order: the envelope is C^{1,1}, so the drift's description has a boundary layer at the faces),
  and keeps Proposition D‴.38's selection of W♯; it leaves the deterministic objects inside the box unchanged, not the stationary
  averages, which weigh the plateau where it meets the box. Numerics: `exp_gauge_plateau.py`. Verdict: classical
  ingredients; the loss of dissipativity under the cohesion leak and its consequences for the gauge are new.
- **Proposition D‴.41 (measurement noise on the gauge).** That isotropic or label noise near a manifold of minimisers produces a drift
  along it, computed by Katzenberger's theorem, is Li–Wang–Arora's ([S64]), the isotropic noise's log-determinant drift among it; that
  label noise drawn afresh at each step selects flat minima (the trace of the Gauss–Newton matrix), with a strength set by step × variance
  and, at a finite step, the regulariser −(1/2η)tr log(I − η∇²L/2), is Blanc–Gupta–Valiant–Valiant's and Damian–Ma–Lee's, the noise
  covariance's shape mattering HaoChen–Wei–Lee–Ma's, the singular limits of noise injection Shalova–Schlichting–Peletier's; that the white
  limit of physical noise is Stratonovich's, that Itô's and Stratonovich's drifts differ by the noise-induced term, and that a coloured
  noise's limit is Itô's, Stratonovich's or neither according to its correlation time against the relaxation time, is classical
  (Wong–Zakai, van Kampen, Horsthemke–Lefever, Kupferman–Pavliotis–Stuart, Pavliotis–Stuart); the continuum filter 1/(1 + kτ) is the
  active Ornstein–Uhlenbeck particle's effective temperature and det(I + τ∇²U) the unified coloured-noise approximation's factor
  (Jung–Hänggi, Maggi et al.); SGD's SDE models are Itô's (Li–Tai–E, Li–Malladi–Arora); correlated and anti-correlated noise act on the
  Hessian's trace, a weighted trace or the weight variances (Orvieto et al., Lucchi et al., Beneventano, Kühn–Rosenow, and an abstract by
  Kim et al. with a correlation-time trade-off); the examples' Itô laws are Ziyin–Li–Ueda's law of balance and Ziyin–Wang–Li–Wu's noise
  equilibria, the diagonal network's label-noise Lasso Pillaud-Vivien–Reygner–Flammarion's ([S69]). What is specific here: for noise
  entering through the error, the identity that the Wong–Zakai term cancels Katzenberger's Itô drift exactly on the manifold — from the
  antisymmetry the limit map's invariance imposes on ⟨∇r ⊗ ∇(∂Φ·U)⟩ — and its extension to smooth noise of any temporal spectrum; the
  coloured filter of label noise's drift along M — the increment's second moment, the cross term ⟨∇ψ·nξ⟩ cancelling by the same
  antisymmetry —, which contains Damian–Ma–Lee's finite-step regulariser as its white case and passes to (1/τ_c)log det(I + τ_c∇²J) for a
  slow Markov noise, a concave potential for every colour, a slow Markov noise acting as a dither of the flow's own mobility; the closed
  forms on the replicator's masses, with their coordinate dependence, and the coloured selection's multistability there. We did not find
  the cancellation or the coloured filter of the drift along M stated. Numerics: `exp_gauge_label.py`. Verdict: classical ingredients, a
  new exact identity with a practical consequence (continuous-time adaptation does not select through its measurement noise) and a new
  law for correlated label noise (its colour sets the implicit regulariser, concave and multistable when coloured).
- **Proposition D‴.42 (several particles per class).** The change to a cluster's barycentre and internal (Jacobi) coordinates, with unit
  Jacobian, the block structure of a translation-invariant interaction, and the non-Gibbs stationary law of a linear system driven by a
  noise not proportional to its mobility are classical; that the positions' Hessian at an exact fit factors into the barycentres' fit part
  and the internal cohesion part follows. What is specific here: the consequence for the replicator's gauge — under the isotropic dither
  the log-determinant selection and its Fixman correction act on the class masses exactly as with one particle per class, with m the
  number of particles, and the within-class fractions neither drift at Katzenberger's order nor diffuse at any γ; an anisotropic dither
  makes them drift at O(γτ) and moves the class masses off the per-class law at unequal splits; under measurement noise the equal split
  is an invariant subspace carrying the one-per-class loop, and the cohesion breaks the cancellation elsewhere. Numerics:
  `exp_gauge_classes.py`. Verdict: a robustness result closing an open item, with its limits.
- **Proposition D‴.43 (the colour of the measurement noise selects the kernels).** That label noise drawn afresh selects a weighted ℓ¹ norm
  on quadratic parametrisations — and so recovers sparse ground truths under ℓ¹'s conditions —, and at a finite σ a Lasso, is
  HaoChen–Wei–Lee–Ma's and Pillaud-Vivien–Reygner–Flammarion's; that a finite step turns the selected trace into a strictly convex
  regulariser, Damian–Ma–Lee's; that the hyperbolic entropy's interpolant interpolates between the kernel and the rich regimes, Woodworth
  et al.'s; that sparse Bayesian learning's cost log det(λI + ΦΓΦᵀ) + yᵀ(λI + ΦΓΦᵀ)⁻¹y keeps its global minimum at the sparsest solution
  with few local minima, Tipping's and Wipf–Rao's, and that its non-separable reweighting succeeds where ℓ¹ fails on a set of positive
  measure, Wipf–Nagarajan's; that log det is a concave surrogate of rank minimised by reweighted trace steps, Fazel–Hindi–Boyd's, and
  reweighted ℓ¹ Candès–Wakin–Boyd's; that gradient noise attracts sparse invariant sets, Chen et al.'s; that Fisher–Rao and
  Wasserstein–Fisher–Rao particle flows solve sparse problems on measures, Chizat's; the dual certificate for ℓ¹, Fuchs' and Tropp's ([S69],
  [S70], [S72]). What is specific here: that the temporal correlation of the label noise sets the sparsity prior of the explicit (Euler)
  update's selection along the exact fits — the trace for white noise, Fazel's log-det heuristic on Φ diag|w|Φᵀ for a slow Markov noise
  (sparse Bayesian learning's log-det term at γ ∝ |w| and λ = 1/τ_c, not its penalty: Wipf–Rao's results on its local minima do not
  transfer), concave for every colour in the continuum — through Proposition D‴.41's filter; the stability law of sparse fits under it (a
  weighted certificate, sufficient, ∝ 1/τ_c beyond a correlation time that grows with the dictionary's coherence); and that this selection
  recovers sparse kernel expansions where basis pursuit fails — the reduced flow, and label-noise descent at small σ —, as Wipf–Nagarajan's
  theorem leads one to expect, at rates sparse Bayesian learning itself matches or beats. We did not find it stated. Numerics:
  `exp_colour_selection.py`. Verdict: a new mechanism — the noise's spectrum as the sparsity prior of a sampled-data adaptation, at a rate
  ∝ hσ² (none in continuous time; Heun's scheme removes most of it) —, not a better sparse solver; for self-adaptive kernel networks with an
  explicit update a slow measurement noise prunes the dictionary, slowly (millions of updates here).
- **Proposition D‴.45 (the free mass–amplitude law as a mirror flow).** In Proposition D‴.44's flat coordinates the free law is the
  Euclidean gradient flow of a separable 2 → 1 reparametrisation, so its being a mirror flow with a start-dependent potential is
  Li–Wang–Lee–Arora's theorem, with one conserved function per unit (Marcotte–Gribonval–Peyré; Azulay et al.'s u ⊙ v); the kernel and
  rich regimes and the ℓ²-to-entropy shape (Woodworth et al.; Ghai–Hazan–Singer), exponentiated gradient as mirror descent (Kivinen–
  Warmuth; Amid–Warmuth), the invariance of stochastic mirror descent's limit under noise in the dual span (Gunasekar et al.; Azizan–
  Hassibi), SGD's Itô noise shrinking the effective initialisation (Pesme–Pillaud-Vivien–Flammarion; Even et al.) and isotropic label noise
  not selecting (HaoChen–Wei–Lee–Ma), and Bregman Lyapunov functions for mirror-descent adaptation laws (Boffi–Slotine) are classical
  ([S72]). What is specific here: the paradigm's saturating unit made explicit — the conserved Q, the mirror coordinate in every regime,
  the conjugate as the mass gained, the saturation's soft wall and entropy tail, the rich limit's weights 1/√β⁰ (a power law in the
  starting mass) with their concave correction, the branches —; the Euler update's per-step change of Q identified with the conserved
  share of Proposition D‴.44's selection, and the mirror update that removes it. Numerics: `exp_mirror_free.py`. Verdict: a classical
  mechanism made explicit for the paradigm's law, with design consequences — the starting masses as a sparsity prior (paid for in time),
  the integrator as a switch for the noise's selection.
- **Proposition D‴.46 (a leak in the mirror coordinates).** A leak in the dual coordinates of a continuous-time mirror flow, its
  regulariser-perturbed rest point, its Bregman (Fenchel-coupling) Lyapunov function and the vanishing-leak selection are established
  (Gao–Pavel; Romano–Pavel; Coucheney–Gaujal–Mertikopoulos for the replicator; Pérolat et al.; Abe et al.; Alvarez–Bolte–Brahic), as are
  the σ-modification's ultimate bound (Ioannou–Kokotović; Ioannou–Sun), the averaging of adaptive systems with the filtered regressor
  (Anderson et al.; Bodson et al.; Sastry–Bodson), the LASSO and its irrepresentable condition (Tibshirani; Zhao–Yu; Wainwright) and
  mirror-descent adaptation laws (Boffi–Slotine) ([S73]). What is specific here: the leak realised on the paradigm's (α, β) law with
  the conservation law kept, its ultimate bound in the adaptive loop with this law, the averaged equilibrium with the error filter's
  Gram as a monotone equation for this law, and the starting masses setting the regime — a ridge, or a LASSO with penalty σλ₁ ∝ σ/√β⁰
  that keeps the target's kernels under the irrepresentable condition. Numerics: `exp_mirror_leak.py`. Verdict: classical ingredients
  assembled into a robust version of the paradigm's mass law, whose small-start limit selects the target's kernels under the (weighted) irrepresentable condition.
- **Proposition D‴.47 (the free law's branches).** Exponentiated gradient and its Bregman-projection bias (Kivinen–Warmuth; Gunasekar
  et al.), Hessian–Riemannian gradient flows on constrained domains (Alvarez–Bolte–Brahic), nonnegative least squares by
  overparametrised gradient flow (Chou–Maly–Verdun et al.), the hyperbolic entropy of u² − v² (Woodworth et al.; Ghai–Hazan–Singer), the
  uniqueness of nonnegative sparse solutions and nonnegative super-resolution without separation (Donoho–Tanner; Bruckstein–Elad–
  Zibulevsky; Schiebinger–Robeva–Recht; Eftekhari et al.) and the exponential penalty in linear programming (Cominetti–San Martín; Weed)
  are classical ([S74]). What is specific here: the paradigm's own law becomes an exponentiated gradient on its branches with a temperature
  read off the conservation law, so that small starts freeze signs; the limit for every start — mixed free and branch particles, inexact
  fits — as the Bregman projection onto the sign-constrained least-squares fits with the extinction of the vanishing coordinates
  (completing D‴.45(iii)); B(0, w) = β/(γm); the two rich limits' temperatures; and the scale-independent selection of a positive target's
  kernels. Numerics: `exp_mirror_branch.py`. Verdict: classical mechanisms identified in the paradigm's law, closing its implicit-bias
  picture for every start.
- **Proposition D‴.48 (the rich transit).** The inverse scale space flow and its adaptive algorithm (Burger–Gilboa–Osher–Xu; Burger–
  Möller–Benning–Osher), saddle-to-saddle dynamics of diagonal linear networks at vanishing initialisation and their identification with
  the adaptive inverse scale space method (Pesme–Flammarion), incremental learning (Berthier; Gissin–Shalev-Shwartz–Daniely), inverse scale
  space limits of mirror flows by Mosco convergence (Berthier–Pillaud-Vivien, 2026), depth-D potentials and polynomial scales (Woodworth et
  al.; Jacot et al.; Chou–Maly–Rauhut), sliding Frank–Wolfe and ADCG, conic particle gradient descent, condensation and early alignment at
  small initialisation, and alternating gradient flows are known ([S75]). What is specific here: the paradigm's free mass–amplitude unit
  as an effective depth-3 unit with an explicit potential, so that its transit has the power scale (mass)^{−1/2} from zero amplitudes and
  the logarithmic one on its branches — a continuum in between, set by α(0)/(mass)^{1/2} —; the sliding transit with the paradigm's moving
  centres (the inactive kernels frozen by the mass in the position law's drift, the creep's profile and the margin condition ε ≪ μ³, ties
  shared); and the positions' metric deciding where the massless kernels go (frozen, or condensed at the correlation's peaks) and how small
  the masses must be for one kernel per feature — not the limit's count, which is one per feature under both (round 37 refuted the
  dichotomy of counts). Unbalanced initialisations (Kunin et al. 2024), deep diagonal networks as mirror flows (Labarrière et al. 2024) and
  feature quantisation at small initialisation (Maennel–Bousquet–Gelly 2018) are related. Numerics: `exp_mirror_transit.py`. Verdict: the
  limit theorem Berthier–Pillaud-Vivien's route for an ε-dependent family, (ii) their entropic case; the sliding transit, as a gradient-flow
  limit, apparently new (formal).

### J. Periodic orbits of the finite-m closed loop near the critical points of the averaged loss

- **Prior.** Two-time-scale averaging under attraction (Khalil Ch. 10–11, Sanders–Verhulst–Murdock Ch. 5,
  Sastry–Bodson Ch. 4): hyperbolic equilibria of the averaged system correspond to periodic orbits of the original
  system with the same stability type. The periodic Riccati block-triangularisation of the monodromy matrix is the
  standard slow-manifold reduction for singularly perturbed linear periodic systems.
- **Distinction.** The theorem identifies the averaged system with the gradient flow of the m-particle regularised
  occupation loss J_m (the finite-m form of Theorem H's functional at τ = 0, in the scaling σ = σ₀/k_e), gives the
  O(1/k_e) location of the orbit including the second-order oscillation term of Remark H.7a, and the k_e-uniform basin.
  No source applying this to an adaptive controller read as a measure was found.
- **Verdict: classical tools, applied to a new object; the first rigorous real-loop ↔ averaged-system statement of the
  program.** The remainder regularity was supplied in round 5 (J-1).

### K. Collapsed configurations and the orientation of the orbit

- **Prior.** The isotypic (symmetry-adapted) splitting of a Hessian at a permutation-symmetric point is standard
  equivariant linear algebra; the stability condition σ₀G⁻¹ ≻ R₂ is the local form of the Chizat–Bach first-variation
  optimality condition at a Dirac mass [S3]-type analyses; higher-order averaging corrections to Floquet exponents are
  classical (Sanders–Verhulst–Murdock).
- **Distinction.** The reading of the orientation of the reference orbit as a design input — an orientation-odd O(k_e⁻²)
  shift of the separation exponents that decides the stability of the collapsed configuration in a window around the
  critical leak, computable from the averaged data — was not found. (The orientation effect itself is Remark H.7a's; the
  round-3 search found no precedent.)
- **Verdict: new as a design statement; the mathematics is classical.** Checked in one instance.

### J′. The stationary law of the dithered finite-m closed loop

- **Claim.** With dither τ = τ₀/k_e and leak σ = σ₀/k_e, the finite-m closed loop has an invariant law whose
  parameter marginal is the Gibbs measure exp(−J_m/τ₀) up to O(1/k_e), with an explicit first-order correction driven by
  the first-order averaged drift of Theorem J (Proposition H.7 + Remark H.7a) as a vector field on parameter space; the
  stationary state carries a probability current of order 1/k_e that vanishes exactly when that drift is reciprocal, and the
  stationary parameter process is not time-reversible.
- **Ingredients, all known:** averaging of SDEs with fast variables through the backward-Kolmogorov expansion [S30];
  perturbation of invariant measures through the Poisson equation of the limiting generator [S31], [S23]; hypoellipticity
  and uniqueness [S25]–[S28]; the thermodynamics of non-equilibrium steady states (currents, entropy production,
  detailed balance) [S29]; second-order averaging [Sanders–Verhulst–Murdock].
- **Verdict:** no source found that applies the invariant-measure expansion to an adaptive control loop (a hypoelliptic
  system with a deterministic fast error channel, a periodic reference and a slow Langevin parameter process), nor one that
  identifies the first-order non-equilibrium current of a learning loop with a reciprocity defect of the averaged dynamics.
  The closest objects are the stationary measures of performative/strategic learning dynamics [S14] and the non-reciprocal
  many-body steady states of [S12]; neither has a plant in the loop or a 1/k_e expansion. The proof technique is standard
  (Pavliotis–Stuart; Mattingly–Stuart–Tretyakov); the result is new in its object and in the identification of the
  correction with Theorem I's defect. As with Theorem J, the non-commutativity term ½⟨[∇²𝒥, ∇²𝒰]_G⟩ in the antisymmetric part
  of the second-order averaged drift of a periodic gradient flow is presumably known in the averaging literature; no source
  was located.
- **Corollary J′.3** (the Gibbs measures π_m are chaotic with limit the minimiser of Theorem H's free energy) is Sanov's
  theorem with Varadhan's lemma [S34] applied to a bounded continuous energy; the statement is classical in form
  (Messer–Spohn 1982 for mean-field Gibbs measures); new only as the bridge between Theorems J′ and H.
- **Not checked numerically:** the uniqueness of the invariant measure (Step 1); the circulation of the real loop (below the
  Monte Carlo resolution at the affordable cost; checked on the first-order averaged SDE).

### H.8 Input-to-state stability of the averaged system

- **Prior.** The proof is the energy-dissipation argument of [S3] Thm 3.2, with the perturbation absorbed by Young's inequality and the limit transferred by Talagrand (Otto–Villani).
- **Verdict: elementary given [S3].** Its role is to reduce Conjecture H.6 to an averaging estimate.

---

## What this means for the claims in PARADIGM.md

The paradigm's contribution is a viewpoint with a set of correct, machine-checked theorems that connect
adaptive control of nonlinearly parametrised approximators to optimal transport and mean-field learning
theory. Within it, the constructions that appear new are:

- the closed-loop OT Lyapunov function (B, its packaging);
- OT σ-modification with its reshuffling result (C);
- the induced-metric form of the occupation-measure variable-projection reduction (G);
- the mean-field Langevin bridge with its global convergence statement (H);
- the frozen-measure structure of the exactly averaged loop (H.7), and the characterisation of its non-reciprocity (I);
- the finite-m averaging theorem for the real closed loop (J) and the orientation as a design input for the collapsed
  configuration (K);
- the explicit strict Lyapunov functions for the mixed PE class (E.0, and E.0′ with the sharp slow-regime rate through
  peak factors);
- the closed-loop propagation-of-chaos statement (D′) and the barycentric leak's Lyapunov analysis (C′), both built on
  standard constructions (Sznitman's coupling; the Lloyd descent of coverage control).

The tools are known. The open problems in PARADIGM.md §8 are where new mathematics is required.

## Open questions this audit could not settle

1. [26] (ICML 2026) could not be read in full. Its relation to A and B is unknown.
2. Assignment-free formation control is a likely source of min-over-permutation Lyapunov functions with Dini derivatives. Candidates are Zavlanos–Pappas and Turpin–Michael–Kumar (CAPT). The probes were blocked. This could reduce B's novelty further.
3. The averaging and singular-perturbation literature on neuro-adaptive control was not swept exhaustively for G.
4. Whether [S5]'s hypotheses hold verbatim for the averaged closed loop of H (the finite-m statement, now the conditional Remark H.5).
5. Whether the averaging literature on periodically forced gradient systems already contains the orbit-area (Stokes) term of Proposition H.7 in some form, for example as geometric-phase or Hannay-type corrections. This literature was not swept.
