# Supplementary Note 1: Linear Stability Analysis of JAPI

> Converted from `SupplementaryNote1_May21.pdf`. Text is faithful to the source;
> equations have been transcribed to LaTeX and keep their original numbering.

## Model

This note analyzes the linear stability of a juxtacrine-activator paracrine-inhibitor (JAPI) reaction-diffusion system, in which the activator is membrane-tethered and propagates by cell-cell contact, while the inhibitor is paracrine and diffuses. JAPI is contrasted with classical paracrine-activator paracrine-inhibitor (PAPI) systems, the prototypical implementation of local-activation/lateral-inhibition (LALI) reaction-diffusion patterning, in which both species diffuse. We work in units where the lattice spacing is unity, and assume periodic boundary conditions on a one-dimensional lattice (taken either infinite or sufficiently large that finite-size effects are negligible; small finite-lattice corrections are discussed where relevant).

We consider a JAPI system in which the activator is membrane-tethered and propagates through neighbor-mediated relay, while the inhibitor is diffusible. The system is formulated on a one-dimensional lattice of cells indexed by integer $j$ along the lattice (so that $j$ serves simultaneously as cell label and, in unit-spacing units, as spatial coordinate), with the understanding that all cells are geometrically identical and that the lattice is translationally invariant. The juxtacrine interaction kernel is written as $K_{jl} = K(j-l)$, where $K(r)$ is a symmetric, nonnegative, row-normalized function of the cell-to-cell displacement $r$; its specific form (nearest-neighbor cosine kernel) is given below. In the main text we use the shorthand $Ka$ for the kernel-weighted activator input; here we work with the explicit lattice form $(Ka)_j \equiv \sum_l K_{jl} a_l$. Under these assumptions the dynamics take the form

$$\frac{da_j}{dt} = \beta_a f\!\left(\sum_l K_{jl} a_l,\; i_j\right) - \mu_a a_j \tag{1}$$

$$\frac{di_j}{dt} = \beta_i f\!\left(\sum_l K_{jl} a_l,\; i_j\right) - \gamma_i i_j + D_i \Delta i_j \tag{2}$$

where $a_j$ and $i_j$ denote activator and inhibitor levels in cell $j$; $\beta_a, \beta_i$ are production rates; $\mu_a, \gamma_i$ are degradation rates; and $D_i$ is the inhibitor diffusion coefficient. The symbol $\Delta$ denotes the nearest-neighbor discrete Laplacian on the cell index,

$$\Delta i_j \equiv i_{j+1} - 2 i_j + i_{j-1}, \tag{3}$$

a finite-difference operator between cells (not a continuum derivative within a cell). The Fourier representation of $\Delta$ on the lattice is given in the Fourier Decomposition section below. Both species share the same nonlinear production function $f$, evaluated at the juxtacrine-weighted activator input $(Ka)_j$ and the local inhibitor level $i_j$. This coupling structure ensures that the Fourier transform of the kernel, $\hat{K}(q) = \sum_r K(r) e^{-iqr}$ (defined precisely in the Fourier Decomposition section), appears in both rows of the linear operator derived below (Turing, 1952; Murray, 2003).

**Production function.** For the steady state analysis below we take $f$ to be a two-input Hill function representing competitive inhibition,

$$f(a,i) = \frac{(a/k_a)^{n_a}}{1 + (a/k_a)^{n_a} + (i/k_i)^{n_i}} \tag{4}$$

where $k_a$ and $k_i$ are activation and inhibition thresholds and $n_a, n_i$ are Hill coefficients. We take $f$ to be shared between activator and inhibitor because this reflects the experimental implementation (see Fig. 2 in main text), in which both species are driven by a single transcriptional channel: the synNotch receptor activation directly drives expression of both the membrane-tethered activator and the secreted inhibitor through a common promoter. The form of the mode-dependent linear operator extends to the case of distinct production functions $f_a \neq f_i$ with notational changes; the compact determinant identity in the existence proof below is written for the shared-promoter case used here.

**Scope of the linear stability analysis.** The specific functional form of the production functions determines the location and number of homogeneous steady states, which we analyze for the competitive inhibition case in Supplementary Note 3. The linear stability analysis that follows is, however, agnostic to the specific form of the production functions: it requires only that the steady state $(a_0, i_0)$ exists and that the partial derivatives $f_a$ and $f_i$ evaluated there have the appropriate signs, positive for self-activation and negative for inhibition. The results therefore apply to any reaction-diffusion system with this architecture, not only the competitive inhibition case (Gierer & Meinhardt, 1972; Kondo & Miura, 2010).

## Steady State Analysis

Before linearizing, we establish the fixed point structure of the system. Seeking spatially homogeneous solutions $a_j = a$, $i_j = i$, the discrete Laplacian vanishes and the juxtacrine sum reduces to $\sum_l K_{jl} a_l = \kappa a$, where

$$\kappa \equiv \sum_l K_{jl} \tag{5}$$

is the total juxtacrine coupling weight, independent of $j$ by translational invariance. With the kernel $K_{jl} = K(j-l)$ assumed symmetric ($K(r) = K(-r)$), nonnegative ($K(r) \ge 0$), and row-normalized ($\sum_r K(r) = 1$), the Fourier transform $\hat{K}(q)$ is real, even, and bounded above by $\hat{K}(0) = \kappa = 1$. For unnormalized kernels, the local activator amplification at $q=0$ becomes $\kappa \beta_a f_{a0}$, where $f_{a0} \equiv \partial f/\partial a$ evaluated at the activated homogeneous steady state (its explicit form for the competitive Hill function is given in Eq. (15) below), and the conditions derived below should be replaced accordingly with the substitution $\alpha \to \kappa\alpha$. The system reduces to the single-cell ODE:

$$\frac{da}{dt} = \beta_a f(\kappa a, i) - \mu_a a \tag{6}$$

$$\frac{di}{dt} = \beta_i f(\kappa a, i) - \gamma_i i \tag{7}$$

**Trivial fixed point.** The origin $(a,i) = (0,0)$ is always a fixed point. For $n_a > 1$, $f(0,0) = 0$ and $\partial f/\partial a|_{(0,0)} = 0$, so the Jacobian at the origin is diagonal with eigenvalues $-\mu_a$ and $-\gamma_i$, both strictly negative. The off state is therefore linearly stable and does not support Turing instability whenever $n_a > 1$. (For $n_a = 1$ the Jacobian acquires a finite off-diagonal contribution and origin stability becomes parameter-dependent; we do not consider this case here.)

**Nontrivial fixed points.** Existence of a nontrivial activated fixed point $(a_0, i_0)$ with $a_0 > 0$ depends on both activator and inhibitor parameters in the coupled two-species system. A useful lower bound on the activator production rate required for such a fixed point follows from the reduction of the system to the activation branch alone (called JA, see also main text). A nontrivial fixed point of the JA system satisfies:

$$\beta_a \frac{(\kappa a_0/k_a)^{n_a}}{1 + (\kappa a_0/k_a)^{n_a}} = \mu_a a_0 \tag{8}$$

The onset of bistability in the JA reduction, the parameter value at which a pair of nontrivial fixed points appears, is found by requiring this equation to hold simultaneously with the tangency condition that the production and degradation curves have equal slopes. Solving these two conditions jointly gives the critical activator concentration:

$$a_c = \frac{k_a}{\kappa}(n_a - 1)^{1/n_a} \tag{9}$$

and the critical production rate:

$$\beta_a^c = \frac{\mu_a k_a}{\kappa} \cdot \frac{n_a (n_a - 1)^{1/n_a}}{n_a - 1} \tag{10}$$

Bistability of the JA system requires $n_a > 1$ and $\beta_a > \beta_a^c$. We note that the full bifurcation condition for the coupled $(a,i)$ system is more complex and depends on the inhibitor parameters as well; the JA bistability condition above provides a useful lower bound. A complete analysis of the two-species fixed point structure in dimensionless parameter space, including the dependence on Hill coefficients and the full two-species phase portrait, is provided in Supplementary Note 3.

When the JA bistability condition holds, the one-dimensional JA reduction has three fixed points along the activator axis: the stable off state $(0,0)$, an unstable threshold state, and a stable activated state. In the full two-species system, the inhibitor tracks the activator at steady state through the ratio $i_0/a_0 = \beta_i \mu_a/(\beta_a \gamma_i)$, which follows from dividing the two steady-state equations.

All subsequent linear stability analysis is performed around the nontrivial activated fixed point $(a_0, i_0)$, whose existence we assume in what follows.

## Homogeneous Steady State

We now identify the activated steady state $(a_0, i_0)$ to linearize around. The steady-state values satisfy:

$$\beta_a f(\kappa a_0, i_0) = \mu_a a_0 \tag{11}$$

$$\beta_i f(\kappa a_0, i_0) = \gamma_i i_0 \tag{12}$$

## Partial Derivatives of the Production Function

Linearization requires the partial derivatives of $f$ with respect to each argument. Writing $A = (a/k_a)^{n_a}$ and $I = (i/k_i)^{n_i}$, so that $f = A/(1 + A + I)$, direct differentiation gives

$$f_a \equiv \frac{\partial f}{\partial a} = \frac{n_a A (1+I)}{a (1 + A + I)^2} \tag{13}$$

$$f_i \equiv \frac{\partial f}{\partial i} = -\frac{n_i A I}{i (1 + A + I)^2} \tag{14}$$

Evaluated at the homogeneous state, with $A_0 = (\kappa a_0/k_a)^{n_a}$ and $I_0 = (i_0/k_i)^{n_i}$,

$$f_{a0} = \frac{n_a A_0 (1 + I_0)}{(\kappa a_0)(1 + A_0 + I_0)^2} > 0 \tag{15}$$

$$f_{i0} = -\frac{n_i A_0 I_0}{i_0 (1 + A_0 + I_0)^2} < 0 \tag{16}$$

Note that $f_{a0}$ is evaluated at the juxtacrine input $\kappa a_0$ rather than $a_0$ alone, reflecting the fact that the effective activation signal experienced by each cell in the homogeneous state is amplified by the total coupling weight $\kappa$.

## Linearization

We perturb around the homogeneous state, writing $a_j = a_0 + \delta a_j$ and $i_j = i_0 + \delta i_j$. The perturbed juxtacrine input is

$$\sum_l K_{jl} a_l = \kappa a_0 + \sum_l K_{jl}\, \delta a_l .$$

Expanding $f$ to first order and canceling steady-state terms using the steady-state equations yields the linearized system

$$\frac{d\,\delta a_j}{dt} = \beta_a f_{a0} \sum_l K_{jl}\, \delta a_l + \beta_a f_{i0}\, \delta i_j - \mu_a\, \delta a_j \tag{17}$$

$$\frac{d\,\delta i_j}{dt} = \beta_i f_{a0} \sum_l K_{jl}\, \delta a_l + \beta_i f_{i0}\, \delta i_j - \gamma_i\, \delta i_j + D_i \Delta\, \delta i_j \tag{18}$$

The architecture-specific structure is already apparent: the activator couples spatially through the juxtacrine kernel $K$, while the inhibitor couples through diffusion.

## Fourier Decomposition

Cells are indexed by integers $j$ along the 1D lattice, so $j$ serves simultaneously as cell label and (in unit-spacing units) spatial coordinate; both the discrete Laplacian defined above and the Fourier ansatz below rely on this 1D ordered labeling. Extensions to higher-dimensional lattices replace $j$ with a multi-index without changing the structure of the analysis.

We seek normal-mode solutions of the form

$$\delta a_j(t) = \hat{a}\, e^{\lambda t + iqj}, \qquad \delta i_j(t) = \hat{\imath}\, e^{\lambda t + iqj} \tag{19}$$

where $q \in (-\pi, \pi]$ is the lattice wavenumber and $\lambda$ is the complex growth rate. Under this ansatz the convolution with $K$ diagonalizes,

$$\sum_l K_{jl}\, \delta a_l = \hat{K}(q)\, \delta a_j, \qquad \hat{K}(q) = \sum_r K(r)\, e^{-iqr} \tag{20}$$

At $q = 0$ this reduces to $\hat{K}(0) = \sum_r K(r) = \kappa$, consistent with the homogeneous-state analysis above.

**Discrete Laplacian.** For nearest-neighbor diffusion on a 1D lattice the discrete Laplacian acting on a Fourier mode gives

$$\Delta e^{iqj} = -\Lambda(q) e^{iqj}, \qquad \Lambda(q) = 2(1 - \cos q) \ge 0. \tag{21}$$

In the long-wavelength limit $q \to 0$ we have $\Lambda(q) = q^2 + O(q^4)$, recovering the continuum Laplacian eigenvalue and ensuring consistency with reaction-diffusion theory at large spatial scales.

**Nearest-neighbor juxtacrine kernel.** For a symmetric nearest-neighbor kernel $K(\pm 1) = \tfrac{1}{2}$, $K(r) = 0$ otherwise, the Fourier transform evaluates to

$$\hat{K}(q) = \cos q. \tag{22}$$

This is not a delta function (whose transform would be flat); rather it is a finite-range kernel whose Fourier transform is positive at low $q$ (coherent amplification of long-wavelength modes) and becomes negative near $q = \pi$ (suppression of checkerboard modes). The sign change at high wavenumber is a qualitative feature of discrete juxtacrine coupling that has no analog in diffusive activator spread, though at low wavenumbers the cosine kernel admits an effective-diffusion rewriting (analyzed in the comparison section below). The nearest-neighbor kernel is the natural minimal model for juxtacrine signaling in which cells interact only with immediate neighbors; more extended or graded kernels (e.g. a Gaussian $K(x) \propto e^{-x^2/2s^2}$) can be introduced if a continuously tunable interaction range is desired, but introduce an additional phenomenological length scale not analyzed here.

**Eigenvalue problem.** Substituting the Laplacian and kernel transforms into the linearized dynamics, the system reduces for each mode $q$ to the $2 \times 2$ eigenvalue problem

$$\lambda \begin{pmatrix} \hat{a} \\ \hat{\imath} \end{pmatrix} = A(q) \begin{pmatrix} \hat{a} \\ \hat{\imath} \end{pmatrix} \tag{23}$$

with linear operator

$$A(q) = \begin{pmatrix} -\mu_a + \beta_a f_{a0} \hat{K}(q) & \beta_a f_{i0} \\[2pt] \beta_i f_{a0} \hat{K}(q) & -\gamma_i + \beta_i f_{i0} - D_i \Lambda(q) \end{pmatrix} \tag{24}$$

This is the central linear operator for JAPI.

## Dispersion Relation

The eigenvalues of $A(q)$ are found from the characteristic equation

$$\lambda^2 - \tau(q)\lambda + \Delta(q) = 0 \tag{25}$$

where the trace and determinant are

$$\tau(q) = \left[-\mu_a + \beta_a f_{a0}\hat{K}(q)\right] + \left[-\gamma_i + \beta_i f_{i0} - D_i \Lambda(q)\right] \tag{26}$$

$$\Delta(q) = \left[-\mu_a + \beta_a f_{a0}\hat{K}(q)\right]\left[-\gamma_i + \beta_i f_{i0} - D_i \Lambda(q)\right] - \beta_a f_{i0} \cdot \beta_i f_{a0} \hat{K}(q) \tag{27}$$

The two eigenvalues are

$$\lambda_\pm(q) = \tfrac{1}{2}\left[\tau(q) \pm \sqrt{\tau(q)^2 - 4\Delta(q)}\right] \tag{28}$$

and the relevant dispersion relation is $\lambda_{\max}(q) = \max\{\mathrm{Re}(\lambda_+), \mathrm{Re}(\lambda_-)\}$. A patterned state arises when $\lambda_{\max}(q) > 0$ for at least one mode $q \neq 0$ while the $q = 0$ mode remains stable, and the dominant pattern wavelength at onset is set by the mode $q^*$ that maximizes $\lambda_{\max}(q)$.

## Stability Conditions

**Homogeneous stability.** The $q = 0$ mode corresponds to spatially uniform perturbations. Since $\Lambda(0) = 0$ and $\hat{K}(0) = \kappa$, stability requires

$$\tau(0) < 0 \quad \text{and} \quad \Delta(0) > 0. \tag{29}$$

These are the standard Routh-Hurwitz conditions for the local reaction system, ensuring that the homogeneous steady state is stable to uniform perturbations.

**Spatial instability.** A patterned state emerges when the homogeneous conditions hold yet $\lambda_{\max}(q) > 0$ for some $q \neq 0$. Under the trace-negative conditions verified in the next subsection, this is equivalent to $\Delta(q)$ becoming negative at some nonzero wavenumber, here controlled by the interplay between $\hat{K}(q)$ and $\Lambda(q)$ rather than by two diffusion coefficients. The existence of parameter regimes where this occurs is proven in the next subsection.

## Existence of Turing Instabilities in JAPI

We show that the JAPI operator $A(q)$ derived above admits finite-wavenumber instabilities of the homogeneous on state under the standard local activator-inhibitor prerequisites, establishing that the Turing regime is accessible to the JAPI architecture. The closed-form threshold, mode selection, and unstable band are the object of future theoretical work; here we establish existence.

**Setup and prerequisites.** Introduce the abbreviations

$$\alpha \equiv \beta_a f_{a0} > 0, \qquad \nu \equiv -\beta_i f_{i0} > 0, \tag{30}$$

where positivity follows from $f_{a0} > 0$ and $f_{i0} < 0$ at the activated steady state.

*Sign convention.* The minus sign in the definition of $\nu$ is chosen so that $\nu$ is positive; with this convention $\beta_i f_{i0} = -\nu$, so the (2,2) entry of $A(q)$ reads $-\gamma_i - \nu - D_i\Lambda(q)$. All expressions below (in particular the determinant identity in Eq. (32)) follow this convention. An alternative convention in which $\nu$ is defined as $\beta_i f_{i0}$ (negative-valued) would give equivalent results with sign flips wherever $\nu$ appears; e.g., the leading term of Eq. (32) would read $\mu_a(\gamma_i - \nu)$ rather than $\mu_a(\gamma_i + \nu)$, with the same numerical value because of the sign flip in $\nu$.

The local prerequisites for a Turing-type instability of the on state are:

- **(i)** $\alpha > \mu_a$ — local activator self-amplification at $q = 0$
- **(ii)** $\alpha - \mu_a - \gamma_i - \nu < 0$ — Routh-Hurwitz trace condition at $q = 0$
- **(iii)** $\mu_a \nu > \gamma_i(\alpha - \mu_a)$ — Routh-Hurwitz determinant condition at $q = 0$

Conditions (ii) and (iii) together ensure homogeneous stability of the on state.

We additionally assume that the lattice kernel and Laplacian eigenvalue satisfy the normalization properties

$$\hat{K}(0) = 1, \quad \Lambda(0) = 0, \quad \Lambda(q) > 0 \text{ for } q \neq 0, \quad \hat{K}(q) \le 1, \tag{31}$$

and that the lattice spectrum contains at least one allowed nonzero mode in a sufficiently small neighborhood of $q = 0$ where the continuity arguments below apply. This is automatic for an infinite lattice or for a sufficiently large finite lattice, but should be checked mode-by-mode on small finite lattices. For the cosine kernel $\hat{K}(q) = \cos q$, this holds in the continuous-$q$ or sufficiently large finite-lattice limit whenever $\alpha > \mu_a$; on small finite lattices the allowed modes must be checked explicitly.

**Determinant structure.** Direct expansion of $\det A(q)$ using Eq. (24), after using $\beta_a f_{i0} \cdot \beta_i f_{a0} = -\alpha\nu$, yields

$$\Delta(q) = \mu_a(\gamma_i + \nu) - \alpha\gamma_i \hat{K}(q) + D_i \Lambda(q)\left[\mu_a - \alpha\hat{K}(q)\right]. \tag{32}$$

The right-hand side separates into a $q$-dependent local part, $\mu_a(\gamma_i + \nu) - \alpha\gamma_i\hat{K}(q)$, and a spatial-coupling part $D_i\Lambda(q)(\mu_a - \alpha\hat{K}(q))$ that is linear in $D_i$. The coefficient of $D_i$ is determined by the sign of $\mu_a - \alpha\hat{K}(q)$: for modes where $\alpha\hat{K}(q) > \mu_a$, the coefficient is negative, and $\Delta(q)$ decreases linearly as $D_i$ increases.

**All modes are stable at $D_i = 0$.** Before showing that finite $D_i$ can destabilize a nonzero mode, we verify that no spatial instability exists in the absence of inhibitor diffusion. At $D_i = 0$,

$$\Delta(q, 0) = \mu_a(\gamma_i + \nu) - \alpha\gamma_i\hat{K}(q). \tag{33}$$

Using $\hat{K}(q) \le 1$,

$$\Delta(q,0) \ge \mu_a(\gamma_i + \nu) - \alpha\gamma_i = \mu_a\nu - \gamma_i(\alpha - \mu_a) > 0 \tag{34}$$

by condition (iii). The trace is

$$\mathrm{tr}\,A(q,0) = -\mu_a + \alpha\hat{K}(q) - \gamma_i - \nu \le -\mu_a + \alpha - \gamma_i - \nu < 0 \tag{35}$$

by condition (ii). All modes are therefore linearly stable at $D_i = 0$, regardless of wavenumber. Spatial instability in JAPI is consequently diffusion-driven: it requires $D_i > 0$ and emerges as $D_i$ increases beyond a finite threshold, in the standard Turing sense.

**Existence of a destabilizing finite-wavenumber mode.** By continuity of $\hat{K}(q)$ and condition (i), there exists $q^* \neq 0$ sufficiently close to $q = 0$ such that both

$$\alpha\hat{K}(q^*) > \mu_a \tag{36}$$

$$N(q^*) \equiv \mu_a(\gamma_i + \nu) - \alpha\gamma_i\hat{K}(q^*) > 0 \tag{37}$$

hold. The first follows from $\alpha\hat{K}(0) = \alpha > \mu_a$, and the second follows from $N(0) = \mu_a\nu - \gamma_i(\alpha - \mu_a) > 0$ (condition iii), both by continuity.

For such a $q^*$, write

$$\Delta(q^*) = N(q^*) - D_i M(q^*), \tag{38}$$

where $M(q^*) \equiv \Lambda(q^*)\left[\alpha\hat{K}(q^*) - \mu_a\right] > 0$. The determinant $\Delta(q^*)$ is therefore a strictly decreasing linear function of $D_i$ with positive intercept $N(q^*)$ and positive slope magnitude $M(q^*)$. It vanishes at the finite positive threshold

$$D_{\mathrm{crit}}(q^*) = \frac{N(q^*)}{M(q^*)} > 0, \tag{39}$$

and becomes negative for all $D_i > D_{\mathrm{crit}}(q^*)$. At this mode, the JAPI operator therefore acquires a positive-real-part eigenvalue, while the homogeneous mode remains stable because $q = 0$ is unaffected by $D_i$.

**Existence statement.** Under conditions (i)-(iii) and the lattice spectral assumptions stated above, the JAPI operator $A(q)$ admits at least one wavenumber $q^* \neq 0$ at which the homogeneous on state is destabilized by sufficiently strong inhibitor diffusion, while the homogeneous mode remains stable. This is the defining signature of a Turing-type instability.

Because the off state $(0,0)$ remains linearly stable for $n_a > 1$ (see Trivial fixed point analysis above), this argument establishes linear instability of the activated homogeneous branch only, not global convergence of the nonlinear system to a patterned state. Nonlinear pattern selection, saturation, and the possibility of basins of attraction containing the off state require numerical simulation or weakly nonlinear analysis, and are not the object of this note.

The JAPI architecture therefore admits Turing instabilities under the same type of local activator-inhibitor prerequisites familiar from classical PAPI, a comparison developed in detail in the next subsection, with the destabilizing spatial threshold provided by $D_i$ alone rather than by a ratio of two diffusion coefficients. The precise threshold

$$D_T^{\mathrm{JAPI}} = \inf_{\substack{q \neq 0 \\ \alpha\hat{K}(q) > \mu_a}} D_{\mathrm{crit}}(q) \tag{40}$$

(replaced by a minimum over allowed lattice modes $q_m = 2\pi m/N$ on a finite periodic lattice of $N$ cells) depends on the lattice geometry through $\hat{K}(q)$ and $\Lambda(q)$, and on the lattice spectrum on finite lattices. For the 1D nearest-neighbor cosine kernel in the continuous-$q$ limit, the infimum is attained in the interior of the allowed band under the stated conditions and admits a direct closed-form expression by minimization of $D_{\mathrm{crit}}(q)$ over $q$. We do not pursue this in detail here; closed-form analysis of $D_T^{\mathrm{JAPI}}$ for general kernels, finite-lattice corrections, higher-dimensional geometries, and the relationship to the classical PAPI threshold, are the object of future theoretical work.

**Worked example.** The following parameter set is presented as a local-Jacobian demonstration of existence: it satisfies conditions (i)-(iii) and the kernel assumptions, and so admits a Turing instability of the operator $A(q)$. We do not claim it corresponds to a particular set of Hill-function parameters $(\beta_a, \beta_i, k_a, k_i, n_a, n_i)$ satisfying the homogeneous steady-state equations; for the competitive Hill production function, the steady-state constraint imposes $\alpha/\mu_a < n_a$ and $\nu/\gamma_i < n_i$, so realizing the example below from the Hill model would require $n_a > 5$ and $n_i > 10$, which are higher than the experimentally measured Hill coefficients in the implementation analyzed in the main text. A Hill-realizable parameter set in the irregular regime relevant to the experimental implementation is analyzed in future theoretical work; the example here serves only to exhibit a Jacobian-level parameter set in which the existence proof applies.

Take

$$\alpha = 5, \quad \mu_a = 1, \quad \gamma_i = 0.5, \quad \nu = 5, \tag{41}$$

with the 1D nearest-neighbor kernel $\hat{K}(q) = \cos q$. Verification of (i)-(iii):

- **(i)** $\alpha - \mu_a = 4 > 0$ ✓
- **(ii)** $\alpha - \mu_a - \gamma_i - \nu = -1.5 < 0$ ✓
- **(iii)** $\mu_a\nu - \gamma_i(\alpha - \mu_a) = 3 > 0$ ✓

Pick $q^*$ with $\cos q^* = 0.8$. Then $\alpha\hat{K}(q^*) = 4 > \mu_a$, $\Lambda(q^*) = 0.4$, and

$$N(q^*) = 1 \cdot 5.5 - 5 \cdot 0.5 \cdot 0.8 = 3.5 > 0$$
$$M(q^*) = 0.4 \cdot (4 - 1) = 1.2$$
$$D_{\mathrm{crit}}(q^*) = 3.5/1.2 \approx 2.92.$$

For any $D_i > 2.92$, the mode $q^*$ satisfies $\Delta(q^*) < 0$. At $q = 0$, both Routh-Hurwitz conditions hold regardless of $D_i$: $\tau(0) = -1.5 < 0$ and $\Delta(0) = 3 > 0$, so the homogeneous mode remains stable. The system is therefore in a Turing regime for $D_i > 2.92$ at this parameter set. The chosen $q^*$ (with $\cos q^* = 0.8$) is illustrative; direct numerical minimization of $D_{\mathrm{crit}}(q)$ over the allowed band gives the true onset $D_T \approx 2.46$ for this parameter set, at the optimal mode $\cos q_T \approx 0.65$, $q_T \approx 0.86$ rad (wavelength $\approx 7.3$ cells).

**Comment on the geometric threshold.** For the 1D nearest-neighbor cosine kernel, the expression

$$D_{\mathrm{geom}} = \frac{\alpha\gamma_i}{2(\alpha - \mu_a)} = 0.3125 \tag{42}$$

for the parameters above, obtained by requiring that the determinant minimum lie within the physical interval of accessible wavenumbers, has appeared in informal analyses as a candidate threshold. $D_{\mathrm{geom}}$ is strictly necessary for finite-$q$ instability (it gives the lower bound on $D_i$ for the minimum of $\Delta(c)$ over $c = \cos q$ to fall in the physical band) but is not the threshold itself: $D_T^{\mathrm{JAPI}} > D_{\mathrm{geom}}$ in general. A direct numerical minimization of $D_{\mathrm{crit}}(q)$ gives $D_T \approx 2.46$ for the parameter set above. A closed-form expression for $D_T^{\mathrm{JAPI}}$, including for general kernels, is developed in future theoretical work.

## Structure of the Linear Operator: Comparison with Classical LALI

We now compare the JAPI linear operator derived above with the corresponding operator for a classical paracrine activator-paracrine inhibitor (PAPI) system. The comparison clarifies which features are shared between the two architectures and which are specific to JAPI.

For PAPI, an identical linearization procedure yields

$$A_{\mathrm{PAPI}}(q) = \begin{pmatrix} -\mu_a + \beta_a f_{a0} - D_a \Lambda(q) & \beta_a f_{i0} \\[2pt] \beta_i f_{a0} & -\gamma_i + \beta_i f_{i0} - D_i \Lambda(q) \end{pmatrix} \tag{43}$$

where $D_a$ is the activator diffusion coefficient and $\Lambda(q)$ is the same discrete Laplacian eigenvalue as in the JAPI case. A note on notation: $f_{a0}$ in the PAPI operator is evaluated at the homogeneous state $a_0$ (no juxtacrine input), and so is formally distinct from $f_{a0}$ in JAPI which is evaluated at the juxtacrine input $\kappa a_0$. For the symmetric nearest-neighbor kernel $K(\pm 1) = 1/2$ used here, $\kappa = 1$ and the two evaluations coincide; the structural comparison below holds in the general case as well.

Placing the two operators side by side:

$$A_{\mathrm{JAPI}}(q) = \begin{pmatrix} -\mu_a + \beta_a f_{a0}\hat{K}(q) & \beta_a f_{i0} \\[2pt] \beta_i f_{a0}\hat{K}(q) & -\gamma_i + \beta_i f_{i0} - D_i\Lambda(q) \end{pmatrix} \tag{44}$$

The two operators are identical except in the entries involving the activator's spatial coupling. Specifically:

- **Entry (1,1).** In PAPI, the activator self-term carries $-D_a\Lambda(q)$, a diffusive penalty whose $q$-dependent shape $\Lambda(q)$ is fixed by the lattice geometry and whose amplitude $D_a$ is an independent biophysical parameter, tunable without affecting the local reaction kinetics. In JAPI, this is replaced by $\beta_a f_{a0}\hat{K}(q)$, where the $q$-dependent shape $\hat{K}(q)$ is again fixed by the lattice geometry but the amplitude is set by $\beta_a f_{a0}$, the same coupling that determines the local self-amplification at $q = 0$. The kernel shapes are fixed by lattice geometry in both architectures; what differs is whether the amplitude is independent of the reaction kinetics ($D_a$ in PAPI, yes) or entangled with them ($\beta_a f_{a0}$ in JAPI, no).

- **Entry (2,1).** In PAPI, the inhibitor cross-term $\beta_i f_{a0}$ carries no $\hat{K}(q)$ factor because inhibitor production responds to the local activator perturbation rather than to a juxtacrine-weighted activator input. In JAPI, this entry carries $\beta_i f_{a0}\hat{K}(q)$ because the inhibitor is produced in response to the juxtacrine activator input, which is itself mode-dependent. The kernel therefore modulates not only the activator self-term but also the inhibitor cross-coupling.

- **Entries (2,2) and (1,2).** Identical in both systems. The inhibitor diffusion term $-D_i\Lambda(q)$ and the inhibitory feedback $\beta_a f_{i0}$ are unchanged.

The key consequence is that $D_a$ drops out of the JAPI patterning problem entirely (Marcon et al., 2016). The instability mechanism is qualitatively preserved: both architectures combine local activator self-amplification with inhibitor-mediated long-range suppression. However, JAPI implements the activator-side spatial structure through the juxtacrine kernel rather than through an independently tunable activator diffusion coefficient.

For convenience in the discussion that follows, we denote the diagonal entries of the JAPI operator as $P(q) \equiv -\mu_a + \beta_a f_{a0}\hat{K}(q)$ and $Q(q) \equiv -\gamma_i + \beta_i f_{i0} - D_i\Lambda(q)$ (activator and inhibitor self-couplings, respectively). The inhibitor self-coupling $Q(q)$ is identical in both architectures and always negative: $f_{i0} < 0$ ensures the inhibitor self-loop is stabilizing, and $-D_i\Lambda(q)$ makes $Q(q)$ more negative at higher $q$. The instability therefore requires $P(q) > 0$ at some nonzero wavenumber, meaning $\beta_a f_{a0}\hat{K}(q) > \mu_a$ at those modes.

What differs between the two systems is the parametric structure of $P(q)$. In PAPI, $P_{\mathrm{PAPI}}(q) = -\mu_a + \beta_a f_{a0} - D_a\Lambda(q)$, and the amplitude on the spatial term is $D_a$, an independent molecular parameter that can be tuned without affecting the local reaction kinetics: decreasing $D_a$ keeps $P_{\mathrm{PAPI}}(q)$ positive over a wider range of wavenumbers and widens the unstable band. In JAPI, $P_{\mathrm{JAPI}}(q) = -\mu_a + \beta_a f_{a0}\hat{K}(q)$, and the amplitude on the spatial term is $\beta_a f_{a0}$, which is not independent of the reaction kinetics: the same coupling sets both the $q = 0$ self-amplification (via $\beta_a f_{a0}$ directly) and the $q$-dependent activator coupling (via $\beta_a f_{a0}\hat{K}(q)$). This entanglement is the precise sense in which JAPI is a more constrained system: the qualitative instability mechanism is the same, but the independent tuning of activator spatial spread through a parameter decoupled from the reaction kinetics, which $D_a$ provides in PAPI, has been removed by construction.

For the nearest-neighbor kernel, $\hat{K}(q) = \cos q = 1 - \Lambda(q)/2$. Substituting into the JAPI activator self-coupling gives

$$P_{\mathrm{JAPI}}(q) = -\mu_a + \alpha\hat{K}(q) = (-\mu_a + \alpha) - \frac{\alpha}{2}\Lambda(q), \tag{45}$$

which has exactly the same $q$-dependence as a PAPI activator self-coupling with an effective activator diffusion coefficient $D_a^{\mathrm{eff}} = \alpha/2$. In this sense, the JAPI (1,1) entry is diffusion-like in the nearest-neighbor case, but with the effective diffusion amplitude slaved to the local activator gain $\alpha = \beta_a f_{a0}$ rather than tunable independently. The genuinely architecture-specific structural distinction is therefore not the form of the (1,1) entry, but the kernel filtering of the (2,1) off-diagonal entry: in PAPI the activator-to-inhibitor coupling carries no $q$-dependence, while in JAPI it carries a $\hat{K}(q)$ factor inherited from the juxtacrine activator input. For more extended kernels, the (1,1) entry can also deviate from the standard diffusive $q^2$ form at intermediate and high wavenumbers; in the 1D nearest-neighbor case analyzed here, the deviation is small at low $q$ but the sign change of $\cos q$ near $q = \pi$ remains a qualitative difference from Fickian diffusion. Closed-form analysis of selected wavenumber, instability bandwidth, and parameter space geometry that follow from these distinctions, including the consequences of the (2,1) entry's kernel filtering, are the object of future theoretical work.

## Summary

Linear stability analysis of JAPI yields a $2 \times 2$ mode-dependent growth operator $A(q)$ in which juxtacrine relay enters through the Fourier transform of the neighbor kernel $\hat{K}(q)$, while inhibitor diffusion contributes the stabilizing penalty $-D_i\Lambda(q)$. Under the standard local activator-inhibitor prerequisites for Turing patterning, we prove that the JAPI operator admits finite-wavenumber instabilities, establishing that the Turing regime is accessible to the JAPI architecture. The result is a sufficient condition stated at the Jacobian level: it establishes that the JAPI architecture admits Turing-type instability whenever the local activator-inhibitor prerequisites are met at some activated steady state, independent of whether a given Hill parameterization realizes those prerequisites. Direct comparison with the classical PAPI operator shows that the two systems differ in exactly two entries, both involving the activator spatial coupling, while sharing identical inhibitor diffusion, inhibitory feedback, and local reaction structure.

The instability mechanism is qualitatively preserved in JAPI: Turing-type patterning requires local activator self-amplification ($P(q) > 0$ at some nonzero wavenumber) combined with long-range inhibitory suppression ($Q(q) < 0$ driven by $D_i\Lambda(q)$). What JAPI loses relative to PAPI is $D_a$ as an independent tuning parameter for the width of the unstable band. The activator spatial term is instead fixed by kernel geometry and circuit gain $\beta_a f_{a0}$, which are entangled with the rest of the reaction network. The closed-form threshold expression for the JAPI Turing onset, its precise relationship to the classical PAPI threshold, and the structure of the unstable band are the object of future theoretical work.

## References

1. Turing, A. M. (1952). The chemical basis of morphogenesis. *Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences*, 237(641), 37-72.
2. Gierer, A. & Meinhardt, H. (1972). A theory of biological pattern formation. *Kybernetik*, 12(1), 30-39.
3. Murray, J. D. (2003). *Mathematical Biology II: Spatial Models and Biomedical Applications* (3rd ed.). Springer-Verlag, New York.
4. Kondo, S. & Miura, T. (2010). Reaction-diffusion model as a framework for understanding biological pattern formation. *Science*, 329(5999), 1616-1620.
5. Marcon, L., Diego, X., Sharpe, J. & Müller, P. (2016). High-throughput mathematical analysis identifies Turing networks for patterning with equally diffusing signals. *eLife*, 5, e14022.
