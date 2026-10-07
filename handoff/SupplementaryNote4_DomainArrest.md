# Supplementary Note 4 — Domain Arrest in the Irregular Patterning Regime of Reaction-Diffusion Systems

> Markdown transcription of `SupplementaryNote4_May21.pdf`. Equation numbers match the
> PDF. Equations are in LaTeX; `$...$` inline, `$$...$$` display.
>
> **This is the note that defines "irregular" mechanistically.** Eq. 5 is the existence
> condition for arrest; Eq. 6/7 give the arrested domain half-width.

This note addresses a question raised by the experimental observation that JAPI circuits
implemented in mammalian cells produce robust patterns **despite operating in a parameter
regime where linear stability analysis classifies the homogeneous activated state as
stable.** The same regime is accessible numerically in both JAPI and PAPI architectures and
gives rise to non-periodic, irregular patterns formed of finite-sized activated domains
distributed in space (see main text Fig. 2D for example). This note provides a heuristic
backbone of an analytical mechanism by which this occurs.

The analysis is **general to reaction-diffusion systems with a diffusing inhibitor** and
applies to both JAPI and PAPI architectures: the inhibitor-field calculation requires only
$\beta_i$, $D_i$, $\gamma_i$, and the activation threshold $I_c$, not the form of the
activator spatial coupling. Architecture-specific content enters the result in two places:

1. through the activator-side parameters that set $I_c$ (the inhibitor concentration that
   prevents activation at the front, which depends on the receiver cell's production rate,
   degradation rate, and Hill response), and
2. through the activator front velocity $v$, which determines whether the quasi-static
   approximation underlying the closed-form result is valid.

Within the quasi-static regime (defined below), the closed-form arrested domain size depends
on $(\lambda_i, \beta_i, I_c)$ only; $v$ does not appear at leading order. Throughout this
note, we treat the activator front velocity as an abstract quantity $v$ representing the rate
at which the activated domain expands. In the JAPI architecture, $v$ reaches its maximum
value, $v_a$ (Supplementary Note 2), in the inhibitor-free limit; during pattern formation
$v$ is bounded above by $v_a$ and decreases as inhibitor accumulates around the expanding
domain.

The analysis is **restricted to the irregular regime**, where individual domains form and
stabilize in isolation. The irregular regime arises from both monostable and bistable
parameter conditions in JAPI/PAPI architectures; we comment on the implications for stability
at the end of the existence-condition section. **It does not apply to the Turing regime**,
where domain boundaries are set by the dispersion relation rather than by self-consistent
inhibitor accumulation.

---

## Setup: an isolated activated domain in steady state

In the irregular regime, patterns form by **nucleation**: local fluctuations cross a
threshold, an activated domain appears, expands by relay or diffusion of the activator, and
produces inhibitor that accumulates around it. We analyze a single such domain treated in
isolation, assuming it has reached a steady state in which:

- The domain has half-width $R$ and is bounded by an activation front
- Inside the domain, the activator is at its activated steady state, so all cells within the
  domain produce inhibitor at the maximum rate $\beta_i$
- Outside the domain, no cells are activated and no inhibitor is produced locally
- The inhibitor field is treated as **quasi-static**, instantly tracking the slowly evolving
  domain geometry

The quasi-static treatment corresponds to the **slow-front limit** $v \ll \sqrt{D_i
\gamma_i}$, where the inhibitor relaxation time $1/\gamma_i$ is much shorter than the
timescale on which the activation front advances by one diffusion length. This is a
simplifying limit adopted to enable the closed-form solution presented below; the validity of
this assumption for the experimental system, the corrections at finite $v$, and the breakdown
of the quasi-static framework at high $v$ are addressed in the final section of this note.

We additionally treat the activation front as a **locally planar one-dimensional boundary**.
This approximation is exact for stripe-like domains and accurate for compact two- or
three-dimensional domains when the boundary curvature radius is large compared to the
inhibitor diffusion length $\lambda_i$. For compact two-dimensional domains with $R \sim
\lambda_i$, curvature corrections modify the boundary inhibitor concentration; the
leading-order result of the planar calculation captures the qualitative dependence on
$\beta_i$, $\lambda_i$, and $I_c$, while quantitative agreement requires the full
Green's-function solution. We adopt the planar approximation for transparency and analytical
tractability; numerical simulations in the main text (Fig. 2L-M) do not rely on it.

In one dimension, treating the domain interior as a region of constant inhibitor production
and the exterior as a region with no production, the inhibitor concentration $I(x)$ satisfies
the linear equation:

$$D_i \frac{d^2 I}{dx^2} - \gamma_i I + \beta_i \mathbb{1}_{|x| < R} = 0 \tag{1}$$

where $\mathbb{1}_{|x| < R}$ is the indicator function of the domain. The natural length
scale is the inhibitor diffusion length:

$$\lambda_i = \sqrt{D_i/\gamma_i} \tag{2}$$

the characteristic distance over which the inhibitory halo around the domain decays.

### Table 1: Parameters used in the arrest analysis

| symbol | meaning |
|---|---|
| $R$ | domain half-width |
| $I(x)$ | inhibitor concentration |
| $\beta_i$ | inhibitor production rate inside activated cells |
| $D_i$ | inhibitor diffusion coefficient |
| $\gamma_i$ | inhibitor degradation rate |
| $\lambda_i = \sqrt{D_i/\gamma_i}$ | inhibitor diffusion length |
| $I_c$ | inhibitor concentration that blocks activation at the domain boundary |
| $v_a$ | inhibitor-free activator front velocity |
| $\rho = 2\gamma_i I_c/\beta_i$ | dimensionless arrest threshold (introduced below) |

---

## Solving for the boundary inhibitor concentration

Equation (1) is solvable exactly. By symmetry the solution is even in $x$. Inside the domain,
the general solution combines a particular solution $I_p = \beta_i/\gamma_i$ with the even
homogeneous mode $\cosh(x/\lambda_i)$. Outside the domain, the solution decays exponentially
as $e^{-(|x|-R)/\lambda_i}$ to satisfy $I(x) \to 0$ at infinity. Matching $I$ and $I'$ at $x
= R$ and solving the resulting linear system gives:

$$I(R) = \frac{\beta_i}{2\gamma_i}\left(1 - e^{-2R/\lambda_i}\right) \tag{3}$$

The boundary value has two clean limits:

- **Small domains** ($R \ll \lambda_i$): Taylor expansion gives $I(R) \approx
  (\beta_i/\gamma_i)\cdot(R/\lambda_i)$, so the boundary inhibitor grows **linearly** with
  domain size, since the entire domain contributes to the inhibitor field at the edge.
- **Large domains** ($R \gg \lambda_i$): the exponential vanishes and $I(R) \to
  \beta_i/(2\gamma_i)$, so the boundary value **saturates**, because cells more than a few
  $\lambda_i$ away from the boundary do not contribute to it. The factor of one-half reflects
  the symmetric leakage of inhibitor outward from the boundary.

---

## The arrest condition

The expanding front stalls when the inhibitor concentration at the boundary reaches a
critical threshold $I_c$ sufficient to prevent further activation: $I(R) = I_c$. $I_c$ is
itself a derived quantity that depends on the activator-side parameters governing the
receiver cell's response (production rate, degradation rate, Hill parameters of activation):
the inhibitor level at which the receiver fails to activate is set by the balance between
inhibition and the activator signal arriving from the neighbor. Through this dependence,
activator-side parameters enter the arrest condition (Eq. 6) by setting $I_c$, separately
from the velocity channel discussed in the regime structure section below. **This note takes
$I_c$ as an empirical input to the arrest analysis; deriving $I_c$ from circuit parameters is
part of the future theoretical work.** Substituting (3) gives the arrest condition:

$$\beta_i\left(1 - e^{-2R/\lambda_i}\right) = 2\gamma_i I_c \tag{4}$$

which can be written compactly as $\beta_i h(R/\lambda_i) = \gamma_i I_c$ with the
dimensionless geometric factor $h(x) = (1 - e^{-2x})/2$.

We note that **nucleation**, the formation of the initial activated domain that this analysis
assumes as its starting point, requires perturbations above a finite threshold, and is
documented numerically in main text Fig. S7C-D rather than analytically here.

---

## Existence condition for arrest

Equation (4) has an immediate consequence. Since $1 - e^{-2R/\lambda_i} < 1$ for any finite
$R$, arrest is only possible when:

$$I_c < \frac{\beta_i}{2\gamma_i} \tag{5}$$

If the activation threshold $I_c$ exceeds half the maximum interior inhibitor concentration
$\beta_i/\gamma_i$, the boundary inhibitor can **never** reach $I_c$ regardless of how large
the domain grows. The front never stalls and the system proceeds to uniform activation. This
existence condition therefore defines a quantitative boundary in parameter space between
conditions where isolated domains arrest at finite size and conditions where activation
propagates indefinitely until the system reaches uniform activation.

We note that equation (6), derived in the next section, identifies the half-width at which
the arrest condition is satisfied; stability of the arrested state is addressed there.

---

## Domain size when arrest occurs

When (5) is satisfied, solving (4) explicitly for $R$:

$$R = \frac{\lambda_i}{2} \ln\!\left(\frac{\beta_i}{\beta_i - 2\gamma_i I_c}\right) \tag{6}$$

In dimensionless form, defining $\rho = 2\gamma_i I_c/\beta_i$ as the **dimensionless arrest
threshold** (satisfying $0 < \rho < 1$ when the existence condition holds), Eq. (6) becomes

$$\frac{R}{\lambda_i} = \frac{1}{2}\ln\!\left(\frac{1}{1-\rho}\right) \tag{7}$$

This separates the length scale $\lambda_i$ from the dimensionless arrest threshold $\rho$,
and makes transparent that the existence-condition divergence ($R \to \infty$ as $\rho \to
1^-$) corresponds to approaching the boundary $I_c \to \beta_i/(2\gamma_i)$ from below.

The qualitative dependencies on independent parameters follow directly:

- increasing $D_i$ at fixed $\beta_i, \gamma_i, I_c$ increases $\lambda_i$ while leaving
  $\rho$ unchanged, so **$R$ increases proportionally**;
- increasing $\beta_i$ at fixed $D_i, \gamma_i, I_c$ decreases $\rho$, so **$R$ decreases**;
- increasing $I_c$ at fixed $\beta_i, \gamma_i, D_i$ increases $\rho$, so **$R$ increases**.

The effect of changing $\gamma_i$ in isolation is **ambiguous** because $\gamma_i$ enters
both $\lambda_i$ and $\rho$; the dimensionless form clarifies that the net effect depends on
which other quantities are held fixed.

### Stability of the arrested state

Equation (6) identifies the half-width at which the arrest condition $I(R) = I_c$ is
satisfied; it does not by itself establish dynamical stability of the arrested state. Local
stability requires a boundary-motion law (e.g., $dR/dt = V(I_c - I(R))$ with $V(0) = 0$ and
$V'(0) > 0$) and analysis of small perturbations to $R$. Because $dI/dR > 0$ in the planar
quasi-static calculation, the **geometric condition for dynamical stability is satisfied**: a
small increase in $R$ raises the boundary inhibitor above $I_c$, suppressing further
expansion; a small decrease lowers it below $I_c$, releasing the front. A formal stability
analysis with an explicit front-motion law is the subject of future theoretical work.

Empirically, arrested domains persist in numerical simulations across the parameter ranges
examined here (main text Fig. 2 and Fig. S7), consistent with dynamical stability in **both
the bistable and monostable regimes** accessible to the system. In the bistable case, the
arrested configuration sits between two homogeneous steady states; in the monostable case,
persistence of arrested domains arises despite the absence of an alternative stable
homogeneous state, **distinguishing the irregular regime from a simple bistable picture.**

---

## Activator front velocity and the regime structure of arrest

The closed-form arrest condition derived above (Eq. 6) makes **no reference to the activator
front velocity $v$**: $R$ is determined entirely by inhibitor parameters and the activation
threshold. In the quasi-static limit, the arrested domain size follows from the geometric
self-consistency condition $I(R) = I_c$, which involves only the inhibitor field; the role of
$v$ is to set how quickly arrest is reached, not what size it stalls at.

This $v$-independence holds only in the quasi-static regime, $v \ll \sqrt{D_i \gamma_i}$.
Because $v$ varies in time during expansion (reaching its maximum $v_a$ in the inhibitor-free
limit at the start of the front's advance, and approaching zero near arrest), we use $v_a$ as
the relevant $v$ in regime comparisons: it is the maximum value $v$ takes during expansion,
and serves as a conservative upper bound. If $v_a/\sqrt{D_i\gamma_i}$ places the system in
regime (a), the actual $v(t)/\sqrt{D_i\gamma_i}$ does so at all times during expansion.
Three regimes can be distinguished:

**(a) Quasi-static, $v_a \ll \sqrt{D_i\gamma_i}$.** Equation (6) applies at leading order;
$R$ does not depend on $v$. A heuristic order-of-magnitude estimate for the subleading
correction is $\Delta R \sim v_a/\gamma_i$, valid in the bulk of regime (a) where $dI/dR$
remains of order $\beta_i/(\gamma_i\lambda_i)$. This estimate breaks down near the existence
boundary ($\rho \to 1^-$), where $dI/dR$ vanishes and small finite-velocity reductions in
boundary inhibitor can produce disproportionately large corrections to $R$. A rigorous
moving-front derivation including the velocity-dependent boundary inhibitor concentration is
the subject of future theoretical work.

**(b) Lagging, $v_a \sim \sqrt{D_i\gamma_i}$.** The boundary inhibitor field is substantially
below its quasi-static value at any given $R$. $R$ remains finite but is larger than equation
(6) predicts and depends on $v$ in a nontrivial nonlinear way.

**(c) Runaway, $v_a \gtrsim v_{\text{crit}} \sim \sqrt{D_i\gamma_i}$.** The boundary
inhibitor saturates at a value below $I_c$, the arrest condition is never satisfied, and the
domain expands without bound until external constraints (tissue size, depletion of resources)
intervene.

The transition from (b) to (c) is set by a critical velocity $v_{\text{crit}}$ at which the
dynamically saturated boundary inhibitor falls below $I_c$. Determining $v_{\text{crit}}$
explicitly, and the full functional form $R(v)$ in regime (b), requires a traveling-wave
analysis and is left to future theoretical work.

In regime (a) the dependence of $R$ on inhibitor parameters takes a simple form: $R$
increases monotonically with $\lambda_i$ and decreases monotonically with $\beta_i$. The
dependence on $v$ is regime-specific: $R$ is approximately $v$-independent at leading order
in regime (a), depends nontrivially on $v$ in regime (b), and is undefined in regime (c).

### Empirical placement of JAPI

For the JAPI experimental system, the propagation velocity in L929 fibroblast transceivers (a
configuration that lacks the inhibitor branch and therefore reports $v_a$ directly) is $v_a
\sim 0.13$ mm/day (Santorelli et al., 2024), equivalent to $\sim 1.5 \times 10^{-3}$ µm/s.
Using inhibitor parameters consistent with engineered secreted morphogens (diffusion $D_i
\sim 10$ µm²/s, half-life $t_{1/2} \sim 2$–6 hours, comparable to measured values for
canonical Lefty in zebrafish; Müller et al., 2012), we estimate the inhibitor degradation
rate as $\gamma_i = \ln(2)/t_{1/2} \approx 3.2 \times 10^{-5}$ to $9.6 \times 10^{-5}$
s⁻¹, giving $\sqrt{D_i\gamma_i} \approx 1.8 \times 10^{-2}$ to $3.1 \times 10^{-2}$ µm/s.

The dimensionless ratio $v_a/\sqrt{D_i\gamma_i}$ therefore falls in the range **$\sim
0.05$–0.08, placing the system in regime (a)**: equation (6) applies as the leading-order
description and arrest occurs robustly.

Finite-velocity corrections are bounded in the bulk of regime (a): the heuristic lag scale
$v_a/\gamma_i$ evaluates to $\sim 16$–47 µm, a fraction of the domain half-width $R \sim
50$–100 µm (the half-width corresponding to arrested domain diameters of $\sim 100$–200 µm),
suggesting bulk corrections to $R$ of order $\sim 30$–50% when evaluated with the
conservative upper-bound velocity $v_a$. The actual correction is likely smaller still, since
the front velocity near arrest is well below $v_a$; corrections near the existence boundary
$\rho \to 1$ require the moving-front analysis deferred to future theoretical work. These are
order-of-magnitude estimates; precise placement would require direct measurement of the
diffusion coefficient and degradation rate of the specific inhibitor used in the JAPI
implementation.

The architectural distinction between JAPI and PAPI enters here through both $I_c$ (set by
activator-side parameters of each architecture) and $v_a$: in PAPI, $v_a$ scales as the
square root of the activator diffusion coefficient times the activator's linear growth rate
(the Fisher-KPP front scaling) and therefore depends on $D_a$, while in JAPI, $D_a$ is absent
and $v_a$ is set by the parameters of Supplementary Note 2.

---

## Biological summary

The arrest mechanism analyzed here gives a simple picture of what controls domain size in the
irregular regime of reaction-diffusion systems:

- **Stronger inhibitor production** (higher $\beta_i$) → boundary inhibitor reaches $I_c$ at
  smaller $R$ → **smaller domains**.
- **Wider inhibitor reach** (higher $\lambda_i$ at fixed dimensionless threshold $\rho =
  2\gamma_i I_c/\beta_i$) → inhibitor halo extends further → **larger domains**. Increasing
  $D_i$ at fixed $\beta_i, \gamma_i, I_c$ increases $R$ proportionally to $\lambda_i$; the
  effect of changing $\gamma_i$ in isolation depends on which other quantities are held
  fixed, since $\gamma_i$ enters both $\lambda_i$ and the dimensionless threshold $\rho$.
- **Front velocity** affects domain size in a regime-specific way: $R$ is approximately
  independent of $v$ in the quasi-static regime ($v_a \ll \sqrt{D_i\gamma_i}$), depends
  nontrivially on $v$ when $v_a$ is comparable to $\sqrt{D_i\gamma_i}$, and arrest fails
  entirely when $v_a$ exceeds a critical value. The relevant comparison uses $v_a$, the
  inhibitor-free upper bound on the front velocity; the actual front velocity during
  expansion is bounded above by $v_a$ and approaches zero near arrest.
- **Arrest is only possible when $I_c < \beta_i/(2\gamma_i)$.** When this condition fails,
  the front cannot be stalled and the system proceeds to uniform activation.

This analysis applies to both JAPI and PAPI architectures, with architecture entering through
two channels: $I_c$ (set by the activator-side parameters governing the receiver cell's
response) and $v_a$ (set by how the activator's spatial coupling generates a front velocity).
The inhibitor-field calculation itself is architecture-independent. Within regime (a), the
closed-form result for $R$ (Eq. 6) depends on $(\lambda_i, \beta_i, I_c)$; $v_a$ does not
appear at leading order, but determines regime placement. A formal derivation of $v_a$ in
terms of activator-side circuit parameters, together with the activator-side derivation of
$I_c$ and corrections from inhibitor coupling, is the subject of future theoretical work. In
JAPI, $v_a$ is set by a smaller number of parameters since $D_a$ is absent, making domain
size control more parsimonious through this channel.

---

## References

1. Müller, P., Rogers, K. W., Jordan, B. M., Lee, J. S., Robson, D., Ramanathan, S. &
   Schier, A. F. (2012). Differential diffusivity of Nodal and Lefty underlies a
   reaction-diffusion patterning system. *Science*, 336(6082), 721–724.
2. Santorelli, M., Bhamidipati, P. S., Courte, J., Swedlund, B., Jain, N., Poon, K.,
   Schildknecht, D., Kavanagh, A., MacKrell, V. A., Sondkar, T., Malaguti, M., Quadrato, G.,
   Lowell, S., Thomson, M. & Morsut, L. (2024). Control of spatio-temporal patterning via
   cell growth in a multicellular synthetic gene circuit. *Nature Communications*, 15, 9867.
