# JAPI: what "irregular" means, and how to measure it

Handoff from a chat session, 6-7 October 2026. Everything below was computed
against Ben's own engine (`github.com/BenSwedlund/reaction_diffusion_simulations`,
commit 07a328b), not a reimplementation.

---

## 1. The conceptual answer

"Irregular" is doing at least four jobs in the manuscript and they are
independent of each other:

1. a classification rule in the screen (a central spike expands to finite width and stops);
2. a linear-stability statement (no finite-q instability of the activated homogeneous state);
3. a spatial-statistics phenotype (broad Fourier, monotonic autocorrelation);
4. a mechanism hypothesis (stochastic nucleation, front propagation, arrest).

The key structural fact: **whether a finite-q instability exists, and whether a
stable localized domain exists, are two independent yes/no axes.** The current
classifier assigns one label by precedence and therefore cannot name the cell
where both are true.

A localized domain survives where a homogeneous activated state cannot, because
inhibitor leaks out sideways from a finite patch and the interior therefore sees
less inhibitor than a uniform field would. The domain grows until the interior
can no longer shed inhibitor fast enough, and the front stalls. That gives the
equilibrium width, and it is why the entire region with no activated homogeneous
state (about 30% of the published grid, over 60% at the implementation's
inhibitor range) is still available to arrested domains.

**Operational definition.** An irregular pattern is a field of arrested domains.
Each domain is a deterministic object with a parameter-set width. The
irregularity lives in where and when they nucleated, not in the object.

---

## 2. The regime map, in two gains

Exact closed forms at any activated steady state, verified against finite
differences to eight decimals:

- **alpha = n_a (1 - P)** , the activator gain
- **nu = gamma n_i Q** , the inhibitor gain, with **G = n_i Q**

where P, Q, S are the fractions of the promoter occupied by activator, blocked by
inhibitor, and free (P + Q + S = 1).

The dispersion relation reduces to

**Delta(q)/gamma = (1 + lambda_i^2 Lambda(q))(1 - alpha Khat(q)) + G**

with **lambda_i = sqrt(D_i/gamma)**, the inhibitor range in cell diameters. Gamma
appears nowhere else, which is why it is not an independent knob except through
the trace condition at small gamma.

Three walls:
- **alpha = 1**: no self-amplification below it.
- **G = alpha - 1**: this is the **existence boundary** (saddle-node), not an
  all-on boundary. On the activated branch G - (alpha - 1) is never negative and
  goes to zero exactly at the saddle-node (verified on 2,295 states). Uniform-ON
  comes from the trace condition tau(0) = alpha - 1 - gamma(1 + G) >= 0 instead.
- **G < W(lambda_i, alpha)**: the finite-q instability. The window opens only when
  lambda_i^2 > alpha/(2(alpha - 1)).

**The parts ceiling:** G <= (n_i/n_a) alpha exactly, since Q <= Q + S. Combined
with G > alpha - 1 this gives, when n_a > n_i,

**alpha < n_a / (n_a - n_i)**

At n_a = 10, n_i = 4 the cap is 1.667 and the maximum reachable alpha on the
activated branch is 1.660. That is why the published screen finds zero Turing
instabilities at the measured Hill pair: it is a structural impossibility, not a
sampling gap. No promoter, copy number, degron or inhibitor range fixes it,
because none of those move n_a or n_i.

### Partition of the published 273,600-point grid

Reconstructed exactly from the Methods (10 b_a x 10 b_i x 19 gamma x 12 n_a x 12
n_i = 273,600). Closed-form classifier, JAPI, 1D, D_i = 10:

| class | sets | % |
|---|---|---|
| n_a = 1, excluded by Supp Note 1 | 22,800 | 8.33 |
| no activated homogeneous state | 83,478 | 30.51 |
| uniform ON (trace unstable) | 2,590 | 0.95 |
| **periodic (Turing)** | **10,697** | **3.91** |
| stable, no finite-q instability | 154,035 | 56.30 |

Split by inhibitor range, periodic accessibility is **non-monotonic in
lambda_i**: it peaks at 11.5% around lambda_i = 3.5 cells, falls to about 3.5% at
the implementation's 6.3 cells, and to 0.2% at 10 cells. The implementation sits
on the far side of the peak, so shortening the inhibitor range moves toward the
optimum. The Nodal Finger 1 fusion is therefore a regime knob as well as a
domain-size knob.

Restricted to lambda_i >= 4 cells, 67% of the nucleate-and-arrest sets have no
activated state at all; below lambda_i = 2 cells only 12% do. The aggregate
number is dominated by a short-inhibitor-range corner of the grid that no
implementation in the paper occupies.

---

## 3. Measuring order: what to use and what not to

**Do not compare nearest-neighbour spacing CV against 0.52.** That is the
Rayleigh value for dimensionless Poisson points. Segmented domains have finite
size and touching ones merge into one component, so a minimum separation is built
into the measurement and CV comes out low with no mechanism at all.

**Use g(r) on segmented centroids**, normalized against a null generated in the
same window with the identical segmentation: g(r) = H_data(r) / mean H_null(r).
This absorbs edge effects and irregular masks, so no analytic edge correction is
needed. Pool across replicates by summing raw counts and then dividing, never by
averaging per-replicate curves.

### The null ladder

Push every null through the identical rasterization and segmentation.

1. **CSR** (uniform points). g = 1 by construction.
2. **Germ-grain**: Poisson nuclei grown to discs of the measured mean domain
   radius, merging allowed, no inhibition. Hole at twice the grain radius, flat
   beyond. **This is a model of the JA control.**
3. **Hard-core at r_ex**: minimum centre separation 2 r_ex. Sweep r_ex to find
   what reproduces the observed hole.
4. **Equilibrium hard discs**: Monte Carlo randomization of the observed
   configuration keeping only the exclusion constraint. Stricter than random
   sequential adsorption, and feasible at any density (RSA cannot pack above
   about 55% and will hang).
5. **Selective thinning (Matern II)**: M nuclei, delete any with an older
   neighbour within d, tune d so the survivor count matches.

### Critical confound

**The g(r) hole is not a readout of the inhibitor range.** It is confounded by
how many domains died. A field thinning from 81 nuclei to 19 survivors produces a
hole at 14.5 cells even though lambda_i is 4.5, and the required thinning
distance is set by the survival ratio alone. One parameter fitted to the survivor
*count* correctly predicted the *hole*, which is a genuine test and it passed.

To separate interaction range from thinning history you need the initial
nucleation density, which means timelapse rather than endpoint images. Without
it, a fitted exclusion length is an upper bound on interaction range, not a
measurement of it.

---

## 4. Results

### Which protocol each result used

Two initial-condition protocols appear in the paper and both appear below. Each
result was run under the protocol its own published panel uses.

| protocol | published panels | results below |
|---|---|---|
| `random_uniform_over0`, spike 1.0, no nucleation (synchronous start) | Fig 1C, 1F, 1G | Fig 1C spot panels, b_i sweep, Turing-vs-arrest, inhibitor field |
| `all_off` + `nucleation_rate` 0.01 (continuous asynchronous nucleation) | all Fig 3 panels, matches the experiment | S8D/3H reproduction, JA vs JAPI |

That leaves the g(r) analysis, the null ladder and the thinning test measured
only under the synchronous start, not under the protocol that matches the
experiment. See open item 1: it is the gap most likely to change conclusions, and
the one time the experiment-matched protocol was used it immediately produced an
artifact the synchronous runs never showed (section 5, size filtering).

### Fig 1C spot panels, synchronous start, 3 seeds

| b_i | label | coverage | n | area | hole | CV at t=500 |
|---|---|---|---|---|---|---|
| 5 | "regular spots" | 0.204 | 87 | 23.6 | 8.8 | 0.070 |
| 12 | "irregular spots" | 0.031 | 19 | 16.8 | 14.8 | 0.082 |

Both are far more ordered than Poisson, and their CVs are indistinguishable. The
published visual distinction between them is not a difference in positional
order; it is coverage, count and spacing.

**Neither has converged at the end of the published run.** CV falls monotonically
through every timepoint with no plateau. At t = 50, the run length used for every
experiment-matched panel, the two differ twofold (0.112 vs 0.223). At t = 500,
the run length used for Fig 1C itself, they are the same. **The
regular-versus-irregular distinction is a transient that the figure runs long
enough to erase.**

### b_i sweep, 3 seeds, equilibrium hard-disc null

| b_i | cov | n(t10) | n(t500) | area | morph | hole | trough obs/null | 2nd peak obs/null |
|---|---|---|---|---|---|---|---|---|
| 3 | 0.431 | 13 | 21 | 206 | labyrinth | 8.8 | 0.00 / 0.68 | 1.81 / 1.32 |
| 4 | 0.284 | 70 | 66 | 42.9 | spots | 9.2 | 0.50 / 0.78 | 1.42 / 1.09 |
| 5 | 0.204 | 93 | 87 | 23.6 | spots | 8.8 | 0.35 / 0.68 | 1.51 / 1.18 |
| 6 | 0.160 | 97 | 87 | 18.3 | spots | 9.2 | 0.42 / 0.60 | 1.53 / 1.11 |
| 8 | 0.106 | 88 | 79 | 13.5 | spots | 9.2 | 0.36 / 0.61 | 1.44 / 1.16 |
| 12 | 0.031 | 19 | 19 | 16.8 | spots | 14.8 | 0.31 / 0.56 | 1.67 / 1.34 |

Hex-kernel classifier: b_i = 3 and 3.5 are Turing; b_i >= 4 has no activated
state. So rows 4 through 12 are pure arrest.

**The hole sits at 2 lambda_i = 8.9 cells throughout, except at b_i = 12** where
heavy early thinning (81 -> 19 before t = 10) inflates it.

**Second-shell order exceeds the hard-disc null at every b_i and does not decay**
with distance from the Turing boundary. Leading explanation, untested: the real
interaction is a soft exponential repulsion at range lambda_i, not a hard core,
and soft repulsion orders more than hard discs at the same density. The right
null is a Gibbs process with an exponential pair potential.

### Turing spots vs arrest spots at matched coverage

| | cov | n | hole | 1st peak | trough obs/null | **2nd peak obs/null** |
|---|---|---|---|---|---|---|
| Turing, b_a = 20, b_i = 26 | 0.154 | 99 | 8.2 | 3.10 @ 9.8 | 0.31 / 0.65 | **1.52 / 1.11** |
| Arrest, b_a = 5, b_i = 6 | 0.160 | 87 | 9.2 | 2.87 @ 10.8 | 0.42 / 0.64 | **1.52 / 1.09** |

**Second-shell order does not distinguish mechanism.** This is the third
independent route to the October finding that mechanism does not predict
appearance.

Note that this series has no Turing *spot* regime at b_a = 5: b_i = 3 is a
labyrinth, so point statistics cannot be compared between regimes within the
published series.

### Inhibitor field

Inhibitor in the gaps versus inside domains is **0.76 (arrest) and 0.79
(Turing)**. It only drops about a fifth, so the halos overlap heavily and the
field is a modest modulation on a high plateau rather than isolated bumps in
empty space. Peak wavelength from the activator and inhibitor agree (8.5 vs 9.6
Turing, 10.7 vs 10.9 arrest), so the inhibitor gives a better-sampled spectral
estimate exactly where the activator is too sparse for point statistics.

At steady state the inhibitor is always the Green's function convolution of the
production field in both regimes, so the activator-inhibitor relationship carries
**no mechanism information**. Do not design an experiment around it.

### JA vs JAPI at matched mean spot size, 5 seeds, continuous nucleation

| | mean area | area CV | n | g(r) hole | 2 x grain |
|---|---|---|---|---|---|
| JAPI (b_i = 20) | 71 | **0.42** | 17 | 11.8 | 8.8 |
| JA (b_i = 0) | 80 | **0.76** | 12 | 9.8 | 9.4 |

Two independent signatures of the inhibitor, both measurable on the existing Fig
3L data:

1. **Size uniformity.** Area CV 0.42 against 0.76 at the same mean. Deterministic
   arrest versus whatever size the front happened to reach. This needs no null
   model and no pair correlation.
2. **Exclusion beyond geometry.** JAPI's hole exceeds twice the grain radius by a
   factor of 1.34 and sits at 2 lambda_i = 12.6. JA sits at its geometric floor
   (ratio 1.04).

Worth noting: the circuit **without** the inhibitor is the one that looks
irregular (lumpy, merged, broad size spread). The circuit with it produces round,
uniform, evenly separated spots.

---

## 5. Methodological findings that must carry forward

**Size filtering is mandatory under continuous nucleation.** The raw
connected-component set is dominated by transient single-cell nucleation events.
Unfiltered, area CV came out at 1.6 to 2.2 and the g(r) hole collapsed to 1.8
cells for both JA and JAPI, destroying the comparison. All JA/JAPI numbers above
use components of at least 12 cells. The hole moves from 11.8 at a 12-cell cut to
13.8 at a 16-cell cut, so the cut must be reported with the number. **This
applies equally to image analysis**, where fresh activation events will do the
same thing.

**Thresholding on uniform fields inverts.** A uniformly-off field has near-zero
range, so a relative-amplitude guard labels it ON. Gate on variance first, then
classify uniform fields by level. Same failure mode Nicole documented for Otsu.

**RSA cannot pack above about 55% area fraction** and will hang. Use Monte Carlo
randomization of the observed configuration instead.

**Use the fastest-growing mode, not the determinant minimum**, to predict
wavelength. They differ by about a factor of two away from threshold. Validated:
predicted 10.6, 9.2, 8.5 against measured 9.8-10.8, 8.9-9.8, 8.2-8.9.

---

## 6. Discrepancies found

### In the deposited code

- **No regime classifier exists in the repo.** The procedure is described step by
  step in Methods ("Defining Patterning Regimes"), but nothing implements it. The
  only automated labeling is `analyze_patterns.py`, which is 1D and returns a
  binary `pattern_flag` from a spectral score, with no Turing-vs-irregular
  distinction. The four-way labels behind Fig 1H are not reproducible from the
  deposit.
- **Hex parity bug confirmed.** `simulation_2D.py` shifts even rows
  (`even_offsets` for `r % 2 == 0`). `visualize_2D.py` line 48 shifts odd rows,
  and its own docstring on line 37 claims even-r. A second, flat-top convention
  appears at line 350. Published hex figures misrepresent adjacency.
- **Hill function differs between files.** `simulation_2D.py` computes
  `(aa + basal)/(1 + basal + aa + ii)`; `finding_steady_states.py` computes
  `basal + aa/(1 + basal + aa + ii)`. Identical at `basal = 0.0`, which is what
  both configs use, so nothing published is affected. They diverge the moment
  anyone sets a nonzero basal, and the steady-state solver would then be solving
  a different model than the simulator.
- **Nucleation is not reproducible.** Spike indices are seeded
  (`init_seed = 2`), but nucleation uses the unseeded global legacy RNG.

### Between Methods, Supp Table 1 and code

- **Nucleation increment.** Methods say the activator is increased by the steady
  state or `spike_value`. The code does `+= 2 * a_ss`. Factor of two.
- **Fig 1C initial amplitude.** Methods say random uniform 0 to 2; Supp Table 1
  says `spike_value` 1.00.
- **Screen spike amplitude unreported.** Methods give only "a user-defined
  height". Recovered by matching published labels: below 2.0 nothing nucleates at
  all, and from 2.0 to 20 the Fig 1F sequence reproduces identically. So the
  classification is insensitive to it on that line, but a single global amplitude
  will under-count irregular at parameter sets whose nucleation threshold exceeds
  it.
- **The "multiple regular peaks" clause.** Methods state that sets producing
  multiple regular peaks from a local perturbation, not already classified as
  Turing, were added to the irregular category. So the published irregular class
  is a union of two different outcomes: single arrested domains, and
  multiple-peak outcomes that failed the Turing test.

### Against the October 2026 session notes

- Fig 1C counts per seed are 24, 12, 21 at b_i = 12, not 19 each. The 19.0 is a
  mean. The count is frozen per seed from t = 10 onward at a seed-dependent value.
- b_a = 20, b_i = 26 gives coverage 0.152 here, not 0.334. Probably an
  initial-condition difference; worth resolving since the matched-pair argument
  rests on it.
- b_a = 20, b_i = 28 comes out **trace-unstable**, not linearly stable, on the
  simulator's own hex kernel. The margin is small (tau(0) = +0.023 at b_i = 27),
  so this is near-marginal rather than a gross disagreement, but it means that
  matched pair may straddle a trace boundary rather than the Turing boundary.

---

## 7. Open items, in priority order

1. **Rerun the g(r) analysis under continuous nucleation.** Everything in section
   4 except the JA/JAPI comparison used the synchronous start. The experiment uses
   continuous nucleation. This is the gap most likely to change conclusions.
2. **Soft-repulsion null.** A Gibbs process with an exponential pair potential at
   range lambda_i. If it reproduces the trough and second peak, the entire series
   is explained with no wavelength selection anywhere and the question closes.
3. **Reclassify the published screen** with the closed-form classifier and report
   the two independent axes rather than one precedence-ordered label. Quantify how
   much of the irregular class is the multiple-regular-peaks case.
4. **Statistical hardening**: more seeds, larger fields, threshold sensitivity
   sweep, and synthetic ground-truth tests for the measurement pipeline (perfect
   hex lattice should give CV near zero, Poisson 0.52, hard-core whatever it
   gives). The pipeline has never been run against a known answer.
5. **Pin mu_a from timelapse** so the timescale columns can be compared against
   experiment duration. Without it, growth rates and arrest times are internally
   consistent but not cashable.
6. **Image-side implementation** of g(r) with mask-aware nulls, size filtering and
   survival-fraction covariates. A separate brief exists for this.

---

## 8. Reproducing

Ben's engine at `github.com/BenSwedlund/reaction_diffusion_simulations`, commit
07a328b, used unmodified. Conventions read from source rather than assumed:

- hex neighbours in **even-r** offset, zero-flux boundaries, `nbr_count` varies at
  edges
- Laplacian is the raw sum form, `sum(neighbours) - count * i`, so the hex
  prefactor is **6**, not the continuum-matched 4; `dx = 1.0`
- juxtacrine signal is the row-normalized mean of the 6 neighbours
- `act_decay_rate = act_half_sat = inh_half_sat = 1.0` in all configs, so
  r_a = `act_prod_rate`, r_i = `inh_prod_rate`, gamma = `inh_decay_rate`
- `basal_prod = 0.0` in both committed configs

Supplementary Table 1 is the complete per-panel reproduction key: 59 rows with
panel, dimensions, domain size, initial condition, spike value, nucleation rate,
max steps and all five circuit parameters.

Analysis scripts from this session are throwaway and should be rewritten properly
in the repo rather than ported.
