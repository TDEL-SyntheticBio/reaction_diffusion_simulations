# Supplementary Note 3: Dimensionless Parameterization and Fixed Point Structure for Numerical Simulations

> Converted from `SupplementaryNote3_May21.pdf`. Text is faithful to the source;
> equations have been transcribed to LaTeX and keep their original numbering.

This note serves as a prerequisite for the numerical simulations. JAPI (juxtacrine-activator paracrine-inhibitor) and PAPI (paracrine-activator paracrine-inhibitor) are the two reaction-diffusion architectures compared throughout this paper; full definitions are given in Supplementary Note 1. Here we reduce the dimensional JAPI system to a minimal set of dimensionless groups, establish the fixed point structure in that parameter space, and define the parameter regimes explored in the numerical results. The PAPI nondimensionalization is given alongside for parallel reference. The linear stability analysis of Supplementary Note 1 is formulated in dimensional parameters and does not require this reduction; however, all numerical simulations are reported in the dimensionless units defined here.

## 1. Nondimensionalization

**Dimensional JAPI system.** In JAPI the activator is membrane-tethered, while the inhibitor is a diffusible field. Cells respond only to ligand presented on adjacent cells, not their own, so the autoactivation is non-cell-autonomous. The activator on cell $j$ is induced only by activator-ligand presented on its neighboring cells through synNotch-ligand contact. Both species are produced by the same cell-autonomous response function $f$, evaluated on the signal received from neighbors and on the local inhibitor. With cells indexed by $j$ on a lattice, the dimensional dynamics are:

$$\frac{da_j}{dt} = \beta_a f\!\left((Ka)_j, i_j\right) - \mu_a a_j, \tag{1}$$

$$\frac{\partial i}{\partial t} = D_i \frac{\partial^2 i}{\partial x^2} + \beta_i f\!\left((Ka)_j, i_j\right) - \gamma_i i, \tag{2}$$

where $i_j \equiv i(x_j, t)$, the Hill response function is (Alon, 2007)

$$f(a,i) = \frac{(a/k_a)^{n_a}}{1 + (a/k_a)^{n_a} + (i/k_i)^{n_i}}, \tag{3}$$

and the juxtacrine kernel sums activator over neighbors with no self-sensing:

$$(Ka)_j = \sum_l K_{jl} a_l, \qquad K_{jj} = 0. \tag{4}$$

Writing the kernel as a function of the integer cell separation $r = l - j$ (the coupling is translation-invariant on the lattice, so $K_{jl} = K(r)$), symmetric nearest-neighbor coupling on a 1D lattice gives $K(\pm 1) = 1/2$ and $K(r) = 0$ otherwise; each cell receives the average of its two neighbors. The support of $K$ sets the spatial range over which the activator can act: the nearest-neighbor choice fixes this activator interaction range to a single cell spacing by design, so it does not introduce a length scale competing with the inhibitor diffusion length. The total coupling weight is $\kappa \equiv \sum_l K_{jl} = 1$, so that homogeneous activator levels propagate unchanged through the kernel: $(Ka)_j = \kappa a_0 = a_0$ at any homogeneous state. Spatial position on the lattice is measured in units of cell-cell spacing, so the dimensionless lattice spacing is one by convention; the Fourier kernel eigenvalue $\hat{K}(q) = \cos(q)$ (Supplementary Note 1) follows from this convention, with $q$ measured in radians per cell. Equations (1)-(2) write the activator on the discrete lattice and the inhibitor as a continuous diffusive field for analytical convenience and to match the formulation of the linear stability analysis; in numerical simulations both species are implemented on the same lattice with a discrete Laplacian for the inhibitor.

The dimensional system has seven parameters with units $(\beta_a, \beta_i, \mu_a, \gamma_i, D_i, k_a, k_i)$, plus two already-dimensionless Hill exponents ($n_a$ and $n_i$).

**Choice of scales.** We rescale using the natural scales of the system (Murray, 2003):

$$\tilde{a} = \frac{a}{k_a}, \quad \tilde{i} = \frac{i}{k_i}, \quad \tilde{t} = \mu_a t, \quad \tilde{x} = \frac{x}{\ell}, \quad \ell = \sqrt{\frac{D_i}{\mu_a}}. \tag{5}$$

The time scale $1/\mu_a$ is set by the activator's degradation rate. The two concentration scales $k_a$ and $k_i$ are the Hill thresholds for activation and inhibition. The characteristic length $\ell = \sqrt{D_i/\mu_a}$ is the distance an inhibitor molecule diffuses on the activator time scale; it is not the same as the dimensional inhibitor diffusion length $\lambda_i = \sqrt{D_i/\gamma_i}$, which represents the distance an inhibitor diffuses in its own lifetime. The latter appears in dimensionless units as $1/\sqrt{\gamma}$ (Section 3).

**Substitution: activator equation.** Substituting (5) into (1) and dividing through by $k_a\mu_a$:

$$\frac{da_j}{dt} = \underbrace{\frac{\beta_a}{k_a \mu_a}}_{= \; r_a} f\!\left((Ka)_j, i_j\right) - a_j, \tag{6}$$

where tildes have been dropped: $a$ and $i$ now denote the dimensionless variables of (5), with the Hill thresholds $k_a$ and $k_i$ absorbed into their definition (equivalently, $k_a = k_i = 1$ in these units). The response function (3) accordingly reduces to

$$f(a,i) = \frac{a^{n_a}}{1 + a^{n_a} + i^{n_i}}, \tag{7}$$

which we use throughout the remainder of the note. The kernel weights $K_{jl}$ are dimensionless lattice quantities set by geometry and so are unchanged by the rescaling.

**Substitution: inhibitor equation.** Substituting (5) into (2), the diffusion coefficient picks up a factor $D_i/(\mu_a \ell^2)$:

$$\frac{\partial i}{\partial t} = \underbrace{\frac{D_i}{\mu_a \ell^2}}_{= \; 1 \text{ by choice of } \ell} \frac{\partial^2 i}{\partial x^2} + \underbrace{\frac{\beta_i}{k_i \mu_a}}_{= \; r_i} f\!\left((Ka)_j, i_j\right) - \underbrace{\frac{\gamma_i}{\mu_a}}_{= \; \gamma} i. \tag{8}$$

The choice $\ell = \sqrt{D_i/\mu_a}$ makes $D_i/(\mu_a\ell^2) = 1$ exactly. This is what is meant by "$D_i$ has been absorbed into the length scale": $D_i$ does not appear as a free parameter in the dimensionless equations because it has been used up defining $\ell$; the entire inhibitor diffusion strength now lives in the choice of length unit.

**Dimensionless JAPI system.** Collecting:

$$\frac{da_j}{dt} = r_a f\!\left((Ka)_j, i_j\right) - a_j, \tag{9}$$

$$\frac{\partial i}{\partial t} = \frac{\partial^2 i}{\partial x^2} + r_i f\!\left((Ka)_j, i_j\right) - \gamma i, \tag{10}$$

with five dimensionless groups:

$$r_a = \frac{\beta_a}{k_a \mu_a}, \quad r_i = \frac{\beta_i}{k_i \mu_a}, \quad \gamma = \frac{\gamma_i}{\mu_a}, \quad n_a, \quad n_i. \tag{11}$$

$r_a$ is the dimensionless activator drive (the production-to-degradation ratio); $r_i$ is the analogous quantity for the inhibitor; $\gamma$ is the ratio of inhibitor to activator degradation rates, governing how quickly the inhibitor decays relative to the activator; $n_a$ and $n_i$ govern the steepness of the activation and inhibition responses. The activator diffusion coefficient $D_a$ does not appear: it has been eliminated by the JAPI architecture itself, where the activator is membrane-bound rather than diffusible (Supplementary Note 1).

**Comparison with PAPI.** In PAPI, the activator is also a diffusible field, with diffusion coefficient $D_a$. The dimensional inventory adds $D_a$, and the activator equation gains a Laplacian term. Applying the same rescaling (which uses $\ell = \sqrt{D_i/\mu_a}$ as the unit of length in both architectures) yields the dimensionless PAPI system:

$$\frac{\partial a}{\partial t} = \frac{1}{D}\frac{\partial^2 a}{\partial x^2} + r_a f(a,i) - a, \tag{12}$$

$$\frac{\partial i}{\partial t} = \frac{\partial^2 i}{\partial x^2} + r_i f(a,i) - \gamma i, \tag{13}$$

with one additional dimensionless group:

$$D = \frac{D_i}{D_a}, \tag{14}$$

the inhibitor-to-activator diffusion ratio. We adopt the convention $D = D_i/D_a$ (rather than its reciprocal) so that the Turing-favorable regime corresponds to $D > 1$, consistent with classical reaction-diffusion literature. PAPI therefore has six dimensionless groups: $r_a, r_i, \gamma, D, n_a, n_i$. The same rescaling absorbs $D_i$ into the length scale in both architectures; the parameter $D_a$ that distinguishes PAPI cannot be absorbed into the same length scale and survives as the residual ratio $D = D_i/D_a$.

## 2. Steady State Analysis in Dimensionless Form

Setting spatial gradients to zero, the dimensionless steady states satisfy:

$$r_a f(a_0, i_0) = a_0 \tag{15}$$

$$r_i f(a_0, i_0) = \gamma i_0 \tag{16}$$

**Trivial fixed point.** $(a_0, i_0) = (0,0)$ always satisfies (15), (16). For $n_a > 1$, the activator production $f(a,i) \sim a^{n_a}$ vanishes faster than linearly near the origin, so the linearized dynamics carry no production term and the trivial fixed point is linearly stable; this state therefore does not support Turing instability in the regimes considered here (Supplementary Note 1). The $n_a = 1$ case is parameter-dependent and is not considered.

**Reduction to a scalar root-finding problem.** Dividing (15) by (16) gives

$$i_0 = \frac{r_i}{\gamma r_a} a_0, \tag{17}$$

so any nontrivial fixed point lies on this line in the $(a_0, i_0)$ plane. The two steady-state equations therefore reduce to a scalar root-finding problem along (17); the dynamics themselves remain two-dimensional.

**Nontrivial fixed points.** Substituting (17) into (15) and dividing through by the overall factor of $a$ (which separates the trivial root $a_0 = 0$ already accounted for above), the nontrivial fixed points are the positive roots of

$$g(a) \equiv 1 + a^{n_a} + \left(\frac{r_i}{\gamma r_a}\right)^{n_i} a^{n_i} - r_a a^{n_a - 1} = 0. \tag{18}$$

Note that $g(0) = 1 \neq 0$: the trivial root lives in the full equation $a \cdot g(a) = 0$, not in $g$ itself. For integer Hill exponents $n_a, n_i$, Descartes' rule of signs admits 0, 1 (a double root), or 2 positive real roots, corresponding to the monostable, saddle-node, and bistable cases respectively; for non-integer Hill exponents the same count holds by a generalized argument on the sign changes of $g(a)$.

**Saddle-node boundary.** In the $(r_a, r_i/\gamma)$ parameter plane (at fixed $n_a, n_i$), the boundary between monostable and bistable regimes is the set on which (18) admits a positive double root, i.e. $g(a) = 0$ and $g'(a) = 0$ hold simultaneously. We solve this pair numerically for each $(n_a, n_i)$ to map the bistable region used in the simulations.

**Activator-only limit.** In the limit $r_i \to 0$, (18) reduces to $1 + a^{n_a} - r_a a^{n_a - 1} = 0$, and the saddle-node tangency between the production curve $r_a f(a,0) = r_a a^{n_a}/(1 + a^{n_a})$ and the degradation line $a$ (Strogatz, 1994) gives the closed-form thresholds:

$$a_c = (n_a - 1)^{1/n_a} \tag{19}$$

$$r_a^c = \frac{n_a (n_a - 1)^{1/n_a}}{n_a - 1} \tag{20}$$

$r_a^c$ is therefore the saddle-node boundary in the activator-only slice $r_i = 0$, and serves as a lower-bound calibration for the full inhibited boundary: for $r_i > 0$, the saddle-node boundary lies above $r_a^c$ in the $(r_a, r_i/\gamma)$ plane. The critical value $r_a^c$ decreases monotonically toward 1 as $n_a$ increases, so higher cooperativity makes bistability easier to achieve in this limit.

**Fixed-point structure in the bistable region.** In the bistable region of the inhibited $(r_a, r_i/\gamma)$ plane, the system has three fixed points along the line (17): the stable off state at $(0,0)$, an unstable intermediate fixed point at the smaller positive root of $g$, and a stable activated state at the larger positive root, with $a_0 \to r_a$ and $i_0 \to r_i/\gamma$ for $r_a$ well above the boundary. Stability classifications follow from the standard $2 \times 2$ Jacobian analysis of (15)-(16) (Strogatz, 1994) and are not reproduced here.

**Effect of inhibition.** Increasing $r_i/\gamma$ shifts the saddle-node boundary toward larger $r_a$, shrinking the bistable region. Mechanically, the inhibitor suppresses the production curve at any nontrivial fixed point relative to the activator-only case (Gierer & Meinhardt, 1972), raising the value of the unstable intermediate fixed point in $a$ and shifting the separatrix between the two basins of attraction toward the activated state.

## 3. Parameter Space for Numerical Simulations

The five dimensionless groups (11) define the parameter space explored numerically; for PAPI, the additional group $D = D_i/D_a$ is included, fixed at $D = 10$ across the scan (a standard Turing-favorable regime in classical reaction-diffusion literature). Simulations span both the monostable and bistable regimes, separated by the saddle-node boundary in the $(r_a, r_i/\gamma)$ plane (Section 2), which reduces to $r_a = r_a^c$ in the activator-only limit $r_i \to 0$. The scan is systematic across the full parameter space (cf. Marcon et al., 2016); visualizations of patterning outcomes are necessarily lower-dimensional projections, with two or three of the groups held fixed to render a 2D plane or 3D volume. The choice of which dimensions to display is dictated by visualization rather than by physics, and the held-fixed values are reported alongside the corresponding figure. Full numerical setup (lattice dimensions, boundary conditions, initial-condition protocol, noise amplitude, time horizon, solver and tolerances, and pattern-classification criteria) is specified in Methods.

JAPI and PAPI share the same homogeneous fixed-point structure (Section 2) but differ in their spatial dispersion: in JAPI the activator coupling enters multiplicatively as $r_a f_a \hat{K}(q)$, while in PAPI it enters additively as $-(1/D)\Lambda(q)$. At small $q$ the JAPI cosine kernel can be locally matched to an effective activator diffusion coefficient ($\hat{K}(q) = 1 - q^2/2 + O(q^4)$, equivalent to a continuum diffusion $D_a^{\mathrm{eff}} = r_a f_a / 2$ at long wavelengths; Supplementary Note 1), but the full dispersion relations and nonlinear couplings are not equivalent at any single value of $D$, and the two systems are scanned independently.

Two distinct lengths appear in this formulation: the rescaling length $\ell = \sqrt{D_i/\mu_a}$ used to define dimensionless space (Section 1), and the dimensional inhibitor diffusion length $\lambda_i = \sqrt{D_i/\gamma_i}$, the distance an inhibitor molecule diffuses in its own lifetime. In dimensionless units, $\lambda_i/\ell = 1/\sqrt{\gamma}$, which sets the natural spatial scale of domain size in the irregular patterning regime (Supplementary Note 2). All spatial outputs (domain widths, inter-feature spacing, pattern wavelengths) are reported in dimensionless units, i.e. in multiples of $\ell$.

## References

1. Gierer, A. & Meinhardt, H. (1972). A theory of biological pattern formation. *Kybernetik*, 12(1), 30-39.
2. Strogatz, S. H. (1994). *Nonlinear Dynamics and Chaos: With Applications to Physics, Biology, Chemistry, and Engineering*. Addison-Wesley, Reading, MA.
3. Murray, J. D. (2003). *Mathematical Biology II: Spatial Models and Biomedical Applications* (3rd ed.). Springer-Verlag, New York.
4. Alon, U. (2007). *An Introduction to Systems Biology: Design Principles of Biological Circuits*. Chapman & Hall/CRC, Boca Raton, FL.
5. Marcon, L., Diego, X., Sharpe, J. & Müller, P. (2016). High-throughput mathematical analysis identifies Turing networks for patterning with equally diffusing signals. *eLife*, 5, e14022.
