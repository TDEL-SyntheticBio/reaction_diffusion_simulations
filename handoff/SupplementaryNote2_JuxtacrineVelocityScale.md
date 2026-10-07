# Supplementary Note 2 — Physical Interpretation of the Juxtacrine Activator Velocity Scale

> Markdown transcription of `SupplementaryNote2_May21.pdf`. Equation numbers match the
> PDF. Equations are in LaTeX; `$...$` inline, `$$...$$` display.

Supplementary Note 1 establishes that the activator diffusion coefficient $D_a$ does not
appear in the JAPI patterning problem, and proves the existence of Turing-type
finite-wavenumber instabilities under the standard local prerequisites. The Note 1 analysis
identifies the conditions under which spatially patterned states can arise. The present note
addresses a complementary question: **how rapidly an activated domain expands through the
tissue before inhibitor-mediated arrest.** To this end, we develop the physical
interpretation of what governs the activator's spatial propagation in JAPI: an effective
propagation velocity $v_a$ characterizing the cell-to-cell relay of the membrane-tethered
activator. We derive a minimal threshold-crossing estimate for $v_a$ in an
inhibitor-suppressed limit; a full traveling-wave analysis incorporating the coupled
neighbor dynamics, lattice-discreteness effects, and inhibitor-mediated corrections is the
subject of future theoretical work.

---

## Setup

For reference, the full PAPI and JAPI dynamics take the forms:

**PAPI:**

$$\frac{\partial a}{\partial t} = D_a \nabla^2 a + \beta_a f(a, i) - \mu_a a \tag{1}$$

$$\frac{\partial i}{\partial t} = D_i \nabla^2 i + \beta_i f(a, i) - \gamma_i i \tag{2}$$

**JAPI:**

$$\frac{da_j}{dt} = \beta_a f\!\left(\sum_l K_{jl} a_l,\; i_j\right) - \mu_a a_j \tag{3}$$

$$\frac{di_j}{dt} = \beta_i f\!\left(\sum_l K_{jl} a_l,\; i_j\right) - \gamma_i i_j + D_i \Delta i_j \tag{4}$$

The PAPI activator equation contains the diffusion term $D_a \nabla^2 a$; the JAPI activator
equation replaces it with the discrete kernel sum over neighboring cells. The PAPI equations
are written in continuum, appropriate for the diffusive regime; the JAPI equations are
written on a lattice, reflecting the cell-by-cell architecture of the relay (Murray, 2003).
The juxtacrine kernel $K_{jl}$ is taken to be nonnegative, symmetric, row-normalized, and
self-input free ($K_{jj} = 0$), as stated in Supplementary Note 1.

The full coupled four-equation system is analyzed in Supplementary Note 1 (linearization and
stability) and Supplementary Note 4 (arrest dynamics). For the velocity-scale analysis that
follows, we work in an **inhibitor-suppressed limit** appropriate to the early expansion
phase of an activated domain: before significant inhibitor accumulation, the inhibitor
argument of $f$ contributes negligibly to activator production. Operationally, we set
$f(\cdot, i_j) \approx f(\cdot, 0)$, so the activator equation reduces to a one-species
front-propagation problem driven by activator self-amplification and degradation alone. This
limit also describes exactly the **transceiver** configuration where the inhibitor branch is
genetically absent. The circuit parameters governing $v_a$ in this limit are derived and
discussed below; the inhibitor parameters re-enter the dynamics when the system transitions
from expansion to arrest, which is the subject of Supplementary Note 4.

---

## Physical picture of juxtacrine relay

In PAPI, the activator spreads through a tissue by diffusion. Its reach at any position is
characterized by the diffusion length $\lambda_a = \sqrt{D_a/\mu_a}$, which captures the
steady-state balance between diffusion and degradation in the linearized reaction-diffusion
equation.

In JAPI, the membrane-tethered activator does not diffuse. Instead, an activated cell
presents activator ligand to its immediate neighbors through cell-cell contact. When a
neighbor receives a juxtacrine signal, it begins producing its own activator ligand, which
it presents to its own neighbors. The sensitivity of a cell to its kernel-weighted input
$s_j = \sum_l K_{jl} a_l$ is set by $k_a$, the input concentration at which production is
half-maximal. The activator "spreads" through this relay, advancing one cell at a time as
activator accumulates and drives production in successive neighbors.

The rate of this relay defines an effective propagation velocity $v_a$ for the activation
front. **Rather than a length scale set by diffusion, JAPI has a velocity set by relay
kinetics.** Note that, unlike paracrine diffusion, the relay in JAPI is intrinsically
discrete: the cell is the unit of advance, and the front advances in integer cell-diameter
steps rather than as a continuous spatial process. This discrete relay mode is a lattice
analog of front propagation in reaction-diffusion systems (Kolmogorov et al., 1937). The
precise classification (Fisher-KPP, Zeldovich, or other type), and how it differs from the
corresponding classical PAPI fronts given JAPI's lattice discreteness and the absence of an
independent activator diffusion coefficient, is the subject of future theoretical work.

---

## Threshold-crossing estimate

To make $v_a$ concrete, we consider a single relay step: a previously inactive cell $j$
exposed to juxtacrine signal from already-activated neighbors. Writing the kernel-weighted
incoming signal to cell $j$ as $s_j(t) = \sum_l K_{jl} a_l(t)$, the internal activator
$a_j(t)$ evolves (in the inhibitor-suppressed limit) as

$$\frac{da_j}{dt} = \beta_a F(s_j) - \mu_a a_j \tag{5}$$

where $F(s) = s^{n_a}/(k_a^{n_a} + s^{n_a})$ is the activation-side Hill function and $k_a$
is the half-maximal input concentration of the receiver.

### Two thresholds

Two distinct thresholds appear naturally in a relay step, because the input to a downstream
cell is kernel-weighted rather than equal to the activator level of the upstream cell. Cell
$j+1$ activates when its input $s_{j+1} = \sum_l K_{j+1,l} a_l$ crosses $k_a$. If the only
contribution to $s_{j+1}$ comes from cell $j$ (worst case, no contribution from the other
side), then cell $j$ must accumulate an activator level $a_j = k_a / K_{j+1,j}$ before
$s_{j+1}$ reaches $k_a$; this elevated level is the **output threshold** $\theta_a$. More
generally, $\theta_a \sim k_a / \sum_l K_{lj}$, with the precise factor depending on lattice
geometry and front orientation. For the symmetric nearest-neighbor coupling used here
($K(\pm 1) = 1/2$), the two thresholds differ by an $O(1)$ geometric prefactor, and we
identify $a_j$'s threshold-crossing condition with $a_j = k_a$ for the one-step estimate that
follows. (For identical cells with a delta-function kernel $K_{j+1,j} = 1$ the two thresholds
coincide; the distinction arises whenever the kernel distributes the activator signal across
multiple neighbors.) A self-consistent treatment distinguishing input and output thresholds
enters the full traveling-wave analysis and is deferred to future theoretical work.

### Quasi-static approximation

We assume that during the threshold-crossing event of cell $j$, the input $s_j(t)$ is
approximately constant, with the activated neighbors at their quasi-steady-state activator
levels on the timescale of $j$'s activation. Concretely, $s_j^* = \sum_{l \in A} K_{jl}
a_l^*$, where $A$ is the set of already-active neighbors and $a_l^* \approx \beta_a/\mu_a$
in the well-saturated limit ($F(s_l^*) \to 1$). This is a leading-order approximation; a
self-consistent treatment of the leading-edge cells, where upstream activator levels and
their inputs evolve together, is part of the full traveling-wave analysis deferred to future
theoretical work.

Under this approximation, define $\beta_{\text{eff}} \equiv \beta_a F(s_j^*)$, where
$s_j^*$ is the (constant) input to cell $j$. Starting from $a_j(0) = 0$, the cell ODE has
the closed-form solution

$$a_j(t) = \frac{\beta_{\text{eff}}}{\mu_a}\left(1 - e^{-\mu_a t}\right) \tag{6}$$

The threshold-crossing time $t_{\text{th}}$ at which $a_j(t_{\text{th}}) = k_a$ is

$$t_{\text{th}} = -\frac{1}{\mu_a} \ln\!\left(1 - \frac{\mu_a k_a}{\beta_{\text{eff}}}\right) \tag{7}$$

which is finite and positive provided $\beta_{\text{eff}} > \mu_a k_a$ (else the activated
steady state of cell $j$, $\beta_{\text{eff}}/\mu_a$, lies below $k_a$ and the relay fails
altogether). Including the additional intracellular latency $\tau_{\text{cell}}$ from
receptor activation to surface ligand presentation, the total time for one relay step is
$T_{\text{step}} = \tau_{\text{cell}} + t_{\text{th}}$, and the resulting relay velocity
scale is

$$v_a = \frac{\ell_{\text{cell}}}{\tau_{\text{cell}} + t_{\text{th}}} \tag{8}$$

in physical length units (with $\ell_{\text{cell}}$ the cell-to-cell distance), or
$v_a = 1/(\tau_{\text{cell}} + t_{\text{th}})$ in cell-per-time units.

### Limit behaviors

Three regimes follow directly from the form of $t_{\text{th}}$:

- **Well above threshold** ($\beta_{\text{eff}} \gg \mu_a k_a$): Taylor expansion of the
  logarithm gives $t_{\text{th}} \approx k_a/\beta_{\text{eff}}$, so $v_a \approx
  \ell_{\text{cell}}/(\tau_{\text{cell}} + k_a/\beta_{\text{eff}})$. The relay rate is
  dominated by production rate $\beta_{\text{eff}}$ and sensitivity $k_a$; degradation enters
  only as a higher-order correction.
- **Near threshold** ($\beta_{\text{eff}} \to \mu_a k_a$ from above): the argument of the
  logarithm approaches zero, $t_{\text{th}} \to \infty$, and $v_a \to 0$. The relay slows and
  eventually fails.
- **Subthreshold** ($\beta_{\text{eff}} < \mu_a k_a$): the activated steady state
  $\beta_{\text{eff}}/\mu_a$ lies below $k_a$, and the relay cannot proceed at all.

This estimate captures the regime-dependent role of each circuit parameter, discussed in the
next subsection.

---

## Circuit parameters governing $v_a$

For fixed lattice geometry, contact topology, kernel weights, and front orientation, the
threshold-crossing time $t_{\text{th}}$, and hence the relay velocity scale $v_a$, depend on
**five circuit parameters**, all properties of the activation circuit rather than of the
medium in which signaling occurs. The sensitivities below are computed at fixed $s_j^*$; in
the full coupled dynamics, the same parameters also enter through the upstream activated-cell
levels $a_l^*$ and hence through $s_j^*$ itself.

### Activator production rate ($\beta_a$)

Sets the effective production amplitude $\beta_{\text{eff}} = \beta_a F(s_j^*)$. Higher
$\beta_a$ raises $\beta_{\text{eff}}$, reduces $t_{\text{th}}$, and accelerates the relay.
Well above threshold, $t_{\text{th}} \approx k_a/\beta_{\text{eff}} \propto 1/\beta_a$, so
$v_a$ grows monotonically with $\beta_a$ and saturates at $v_a \to
\ell_{\text{cell}}/\tau_{\text{cell}}$ as $\beta_a \to \infty$.

### Cooperativity of activation ($n_a$)

Sets the sharpness of the Hill response $F(s)$ around $k_a$, and shapes cell behavior on
both sides of $k_a$, not only above it. At fixed input $s_j^*$ (the quasi-static regime
adopted here), only the **value** $F(s_j^*)$ enters $\beta_{\text{eff}} = \beta_a F(s_j^*)$,
not the slope $\partial F/\partial s$, and the effect of $n_a$ on the relay velocity depends
on whether $s_j^*$ sits above or below $k_a$:

- $s_j^* > k_a$: $F(s_j^*)$ rises toward 1 faster as $n_a$ increases, raising
  $\beta_{\text{eff}}$ and **accelerating** the relay;
- $s_j^* = k_a$: $F(s_j^*) = 1/2$ for all $n_a$, so cooperativity has **no effect** at the
  half-maximal input;
- $s_j^* < k_a$: $F(s_j^*)$ falls toward 0 faster as $n_a$ increases, lowering
  $\beta_{\text{eff}}$ and **slowing or preventing** the relay.

The local slope $\partial F/\partial s$ enters in a full dynamic front problem where
$s_j(t)$ evolves during the threshold-crossing event, but does not contribute directly in
the constant-input estimate.

### Activation sensitivity ($k_a$)

The half-maximal input level of the Hill function (with kernel-geometry prefactors absorbed,
as discussed in the Threshold-crossing estimate above). Higher $k_a$ means lower sensitivity
to a given juxtacrine input level, so cell $j$ requires longer to accumulate enough activator
to drive production in its next neighbor. This appears in the relay-time formula as
$t_{\text{th}} = -(1/\mu_a)\ln(1 - \mu_a k_a/\beta_{\text{eff}})$, which increases
monotonically with $k_a$. $k_a$ also sets the **propagation-failure boundary**: when $\mu_a
k_a$ exceeds $\beta_{\text{eff}}$, the activated steady state $\beta_{\text{eff}}/\mu_a$
falls below $k_a$ and the relay cannot proceed.

### Activator degradation rate ($\mu_a$)

Enters in two distinct ways. In the well-above-threshold regime, $t_{\text{th}} \approx
k_a/\beta_{\text{eff}}$ and $\mu_a$ contributes only as a higher-order correction; the relay
rate is essentially independent of degradation. Near threshold, higher $\mu_a$ lowers the
activated steady state $\beta_{\text{eff}}/\mu_a$ toward $k_a$, increases $t_{\text{th}}$
(through the logarithm), and slows the relay. Sufficiently high $\mu_a$ drives the activated
steady state below $k_a$ (equivalently, $\beta_{\text{eff}} < \mu_a k_a$), and the relay can
no longer proceed. **The role of $\mu_a$ is therefore best characterized as setting the
regime in which the relay operates rather than as a leading-order multiplier of velocity.**

### Cell response time ($\tau_{\text{cell}}$)

The latency from receiving juxtacrine signal to presenting sufficient surface ligand to a
neighbor. This summarizes the intracellular kinetics of the relay step — receptor
activation, downstream signaling, transcription, translation, and membrane trafficking of the
new ligand. $\tau_{\text{cell}}$ enters the relay step time additively: $T_{\text{step}} =
\tau_{\text{cell}} + t_{\text{th}}$. The effect of $\tau_{\text{cell}}$ on $v_a$ therefore
depends on the regime: when $\tau_{\text{cell}} \gg t_{\text{th}}$, the relay rate is set
primarily by intracellular latency; when $t_{\text{th}} \gg \tau_{\text{cell}}$, by the
threshold-crossing dynamics.

### Genetic tunability

These five circuit parameters are genetically tunable, in principle: production rate through
promoter strength and translation efficiency, cooperativity through receptor architecture,
sensitivity through receptor affinity, degradation rate through degron tags or fusion
partners, and cell response time through signaling pathway design. The realized tunability in
any particular experimental system depends on cell type, trafficking, and other
context-dependent factors. Stochastic integration and multiplicity of infection of lentiviral
vectors primarily shape variation in production rate; the other parameters would need to be
modified through changes in the protein sequences themselves.

---

## Relationship to $D_a$ and to circuit control

In PAPI, the activator's spatial behavior depends on both the molecular diffusion coefficient
$D_a$ (set by the molecular identity of the activator) and on circuit-encoded kinetic
parameters (such as $\mu_a$, which contributes to the diffusion length $\sqrt{D_a/\mu_a}$).
$D_a$ is a biophysical property of the activator molecule and is not accessible to genetic
manipulation without changing the molecule itself — the central protein-engineering challenge
for synthetic PAPI in mammalian cells.

In JAPI, $D_a$ is absent by construction: the membrane-tethered activator does not diffuse.
For fixed tissue geometry and fixed juxtacrine contact kernel $K_{jl}$, the kinetic
parameters governing the activator's spatial behavior are all circuit-encoded ($\beta_a$,
$n_a$, $\mu_a$, $k_a$, $\tau_{\text{cell}}$). This is a reduction in the total number of
parameters available to control activator dynamics — **JAPI offers fewer handles than PAPI —
but the parameter that has been eliminated is the one that has historically been hardest to
engineer.** The remaining handles, in JAPI, are all genetically tunable rather than tied to
molecular biophysics.

This architectural feature is what makes JAPI experimentally accessible in mammalian cells
without requiring molecular engineering of differential diffusion coefficients.

---

## Dynamic regime and scope

The propagation velocity $v_a$ as derived above characterizes the relay rate in the
inhibitor-suppressed limit $f(\cdot, i_j) \approx f(\cdot, 0)$, corresponding either to the
early expansion phase of an activated domain before significant inhibitor accumulation, or to
the transceiver configuration where the inhibitor branch is absent. As inhibitor accumulates
around an expanding domain, $f(\cdot, i_j)$ is suppressed at the front by the local inhibitor
level, the effective production $\beta_{\text{eff}}$ drops, $t_{\text{th}}$ grows, and the
front eventually arrests. The arrest condition therefore depends critically on the inhibitor
branch: on $\beta_i$, $k_i$, $n_i$, $\gamma_i$, and $D_i$, which together determine the
spatial profile of inhibitor around an activated domain and the threshold at which
$\beta_{\text{eff}}$ falls below $\mu_a k_a$ and relay can no longer advance. The full arrest
dynamics, including the determination of arrested domain size from the inhibitor profile, are
addressed in Supplementary Note 4.

A full traveling-wave analysis yielding $v_a$ in closed form would relax the quasi-static
approximation for $s_j(t)$, treat the coupled dynamics of activator levels across the lattice
self-consistently, incorporate lattice-discreteness effects (such as front pinning,
anisotropy, and propagation failure under thresholding), and include corrections from
inhibitor accumulation and explicit dependence on the inhibitor-side parameters. This is the
subject of future theoretical work.

---

## References

1. Kolmogorov, A. N., Petrovsky, I. G. & Piskunov, N. S. (1937). A study of the diffusion
   equation with increase in the amount of substance, and its application to a biological
   problem. *Bulletin of Moscow State University, Series A: Mathematics and Mechanics*,
   1(6), 1–25.
2. Murray, J. D. (2003). *Mathematical Biology II: Spatial Models and Biomedical
   Applications* (3rd ed.). Springer-Verlag, New York.
