# JAPI patterning: computational findings, October 2026 session

Starting position for the manuscript pass. Everything below was computed in this
session against the engine in `simulation_2D.py`. Code in `analysis_code/`, raw
results in `analysis_data/`, figures alongside this file.

Context: aimed at Reviewer 1's "elephant in the room", namely why periodic
patterning was never observed in vitro, and at Reviewer 2's request for
principled criteria separating irregular patterning from unpatterned noise.

---

## 1. The main result: irregularity comes from asynchronous initiation, not from the parameter regime

Same circuit, same parameters throughout (the parametrized Fig. 2L set:
b_a = 5, b_i = 25, gamma = 0.5, D_i = 20, n_a = 10, n_i = 4). Only the
initiation protocol differs.

| | t=10 | t=50 | t=100 | t=200 | t=350 | t=500 |
|---|---|---|---|---|---|---|
| ongoing nucleation, CV | 0.483 | 0.481 | 0.455 | 0.466 | 0.458 | 0.466 |
| ongoing nucleation, spacing | 4.08 | 4.32 | 4.18 | 4.26 | 4.56 | 4.55 |
| synchronized start, CV | 0.196 | 0.223 | 0.204 | 0.193 | 0.190 | 0.190 |
| synchronized start, spacing | 9.13 | 11.05 | 11.59 | 11.97 | 11.98 | 11.98 |

Under continuous stochastic activation (rate 0.01, the experimental protocol),
nearest-neighbour spacing CV sits at 0.46 to 0.48 for the whole run. Random
Poisson placement is 0.52. The pattern is statistically indistinguishable from
randomly positioned spots, permanently, and spacing is jammed at ~4.5 cells,
about a third of what the inhibitor sets.

Start the identical circuit synchronized and it orders within ten time units:
CV 0.19, spacing relaxing to 12 cells.

Interpretation: each new nucleation event fires wherever a cell happens to
switch on, inside the inhibitor halos of existing spots. This both compresses
the spacing and injects disorder faster than the field anneals it out. The
circuit is not in a regime that produces irregular patterns; it is in a regime
that produces orderable patterns, continuously disrupted by how they start.

**Dose response is not monotonic.** Lowering the firing rate does not recover
order, it only thins the field: rate 0.003 gives CV 0.48 with 46 spots, rate
0.001 gives CV 0.58 with 10 spots, rate 0.0003 nucleates essentially nothing.
Asynchrony is the problem, not its amount.

**Prediction.** Regular patterns should require synchronizing activation, a
pulsed or inducible initiation compressing onset into a window short compared
to pattern formation, with no change to receptor architecture, Hill
coefficients, promoter strengths or protein stability.

**Caveats.** One seed per condition (effect is large, 0.47 vs 0.19, but should
be replicated). The synchronized run needed noise amplitude 2 to nucleate at
these parameters, so the comparison is not amplitude-matched.

Figure: `elephant_nucleation.png`

---

## 2. The protocol discrepancy in Supplementary Table 1

| figures | initialization | nucleation | steps |
|---|---|---|---|
| 1C regime panels | random_uniform_over0, amplitude 1 | none | 50,000 |
| 2L, 3B, 3H, 3L, S8D (experiment-matched) | all_off | rate 0.01, kick 2 | 5,000 |

Every panel compared against cells runs a tenth as long, from a different
initial condition, with continuous nucleation throughout. At 5,000 steps
(t = 50) every pattern in this family is still disordered and still ordering.
So the regime characterization in Fig. 1C and the experiment comparisons are
not the same computational experiment.

---

## 3. Mechanism does not predict appearance

The Turing versus nucleation distinction is real, computable and was validated
cleanly. It controls exactly one thing: whether infinitesimal noise suffices to
initiate. It does not predict what the pattern looks like or how it evolves.

**Validation that stands.** 2D linear stability against 2D simulation at 20
points: 20/20 agreement on periodic-or-not. Predicted wavelength matches the
measured onset wavelength to within FFT bin resolution (predicted 12.1, 11.4,
11.1, 11.1; measured 12.5, 11.1, 11.1, 11.1). In 1D, the solver reproduces all
four Fig. 1F labels exactly, and a single burst gives 4 to 6 domains at Turing
points versus a median of 1 arrested domain at non-Turing points, with no
overlap.

**The null result.** Matched spot pair at b_a = 20, D_i = 10, n_a = n_i = 3,
four seeds each, b_i = 26 (Turing, coverage 0.334) against b_i = 28 (linearly
stable, coverage 0.317):

| | spacing t20 | t100 | t500 | growth | CV t500 |
|---|---|---|---|---|---|
| Turing, b_i=26 | 8.25±0.19 | 9.22±0.16 | 9.70±0.20 | +17.6% | 0.077±0.003 |
| non-Turing, b_i=28 | 8.21±0.16 | 9.27±0.23 | 9.72±0.18 | +18.4% | 0.074±0.006 |

Indistinguishable at every time point, with the difference an order of
magnitude below seed-to-seed spread. Crossing the linear-stability boundary at
matched morphology and density changes nothing observable.

This comparison was only possible once a Turing **spot** regime was found. At
b_a = 5 the Turing window sits entirely at 47-57% coverage, so it is labyrinth
by construction and both Fig. 1C spot panels are non-Turing. Turing spots exist
at b_a = 10 (b_i = 9.6) and b_a = 20 (b_i = 26).

Figures: `lsa_validation.png`, `phase_2d_confirmation.png`,
`matched_pair_dynamics.png`, `full_phase_map.png`

---

## 4. Where the Fig. 1C boundaries actually sit

At b_a = 5, n_a = n_i = 3, D_i = 10, gamma = 0.5:

- Turing instability ends between b_i = 3.6 and 3.7
- activated homogeneous state vanishes (saddle-node) between 3.7 and 3.8
- stripe-to-spot morphological crossover around b_i = 4.2 to 4.6

So the mechanistic and morphological transitions are about 0.5 apart in b_i and
are not the same transition. Scanning across the stability boundary, coverage
goes 0.468 to 0.455, domain count 45 to 51, eccentricity 0.862 to 0.844,
spacing 11.1 to 10.9. No kink in anything.

Only the stripes panel (b_i = 3) is Turing-unstable. Both spot panels,
"reg spots" at b_i = 5 and "irreg spots" at b_i = 12, have (0,0) as their only
homogeneous state and are nucleation-driven.

Figure: `boundary_morphology.png`

---

## 5. Hill coefficients: periodic is inaccessible at the measured pair

Scanning 28,125 parameter sets at n_a = 10, n_i = 4 spanning b_a 1 to 30,
gamma 0.1 to 10, inhibitor range 0.5 to 30 and amplitude 0.1 to 1000: zero
Turing instabilities. Controls: n_a = 4, n_i = 4 gives 3.2% of the same space,
n_a = 10, n_i = 10 gives 1.6%.

**Inhibitor halflife is not an independent knob.** Writing the inhibitor
equation in terms of range sqrt(D_i/gamma) and amplitude b_i/gamma, the
determinant scales as gamma times a gamma-independent bracket, so the sign of
the Turing condition cannot depend on gamma. Confirmed numerically: the Turing
region at gamma = 0.5 and gamma = 2.0 is identical to four decimal places
(0.0864 both). The only gamma effect is at small gamma (below ~0.35 on the
tested line), where the uniform mode destabilizes and the system goes all-on.

Consequence: fixing D_i = 10 while sweeping gamma from 0.1 to 10 traces a
diagonal that spans the relevant reduced space anyway, so the published screen
is more robust to that choice than it looks.

**Why the mCherry receptor did not help.** At the fitted operating point
(b_i/gamma = 50) there is no activated homogeneous state for any n_a from 2 to
10, so swapping the receptor cannot produce periodic patterns there. The
periodic window sits at b_i/gamma of roughly 3 to 7 across every Hill pair that
has one, while the implementation sits at 50. The operating point is
inhibitor-dominated about tenfold.

**Which genetic edits move it.** A PEST on the TF scales activation and
inhibition together, leaving b_i/(gamma·b_a) unchanged at 10, and never enters
the Turing region. A PEST on the secreted inhibitor does nothing (clearance is
set by the medium). What works: weakening inhibitor expression 8 to 12 fold, or
raising activator expression ~8 fold. Both move the ratio toward 1.

Figures: `periodic_accessibility_analysis.png`, `mcherry_slice_analysis.png`,
`genetic_knobs_analysis.png`

---

## 6. Things I claimed and then had to withdraw

Recorded so they are not re-derived.

- **"Coverage controls the coarsening rate."** Correlation −0.92, but driven
  entirely by labyrinth points. Within the spot class it is +0.13. Dead.
- **"Turing patterns lock their wavelength, nucleated ones coarsen."** The
  Turing wavelength is selected at onset at every amplitude tested, and
  everything coarsens eventually; the near-infinitesimal run just coarsens later
  (t ~ 150 to 200). What looked like locking was a too-short run plus integer
  FFT bins. Also, in stripes the coarsening is discrete (100/9 to 100/7) and
  stochastic, so single-seed trajectories are not interpretable.
- **"The published protocol starts at the noise scale."** Measurement artifact
  from evaluating at t = 0 before any pattern existed.

---

## 7. Open, worth doing

1. Replicate the nucleation result with several seeds and amplitude-matched
   initiation.
2. Deeply-Turing against deeply-nucleation-driven at matched coverage, moving
   b_a and b_i together. The matched pair tested only the immediate
   neighbourhood of the boundary.
3. Measure n_i for the anti-mCherry circuit. Only n_a (~4) was measured; n_i was
   measured only for anti-GFP (~4). One soluble-inhibitor titration, same assay
   as S6C.
4. A separate, measured order axis for the manuscript, grounded in spacing CV
   against the 0.52 random baseline and in whether ordering arrests. Spectral
   prominence cleanly separates pattern from noise (1.8 for a noise field
   against 400 to 900 for real patterns), which answers Reviewer 2 point 2.

---

## 8. Earlier in this session: the summer students' work

Deliverable is `ReactionDiffusionAnalysesSummary_v2.pdf`, Nicole's write-up with
sections added by me (Otsu failure on uniform fields and the variance-gate
approach, prominence versus sharpness, spacing CV against the random baseline,
an open-questions section, b_i and regime columns on the threshold table).

What that work established, all of which the above builds on:

- Fig. 1C reproduced for JAPI and PAPI. The published recipe is
  random_uniform_over0 with spike_value 1, no nucleation. The Methods text
  saying "random uniform 0 to 2" contradicts Supplementary Table 1 and lands
  outside the productive window; the irregular regime only nucleates in a narrow
  amplitude window around 1.
- Termination: mean per-cell change below ~1e-3, which stops runs at 10,000 to
  20,000 steps rather than 50,000. Known issue: the mean is taken over all
  cells, so sparse fields are diluted by their inactive background and stop
  earlier in their own evolution. Fix is to compute over active cells or
  normalize by activated fraction. Not yet done.
- Thresholding: Otsu fails on uniform fields (all-on read as 0% active, all-off
  as 100%). Adopted approach is a variance gate first, then Otsu on anything
  patterned. Component count is robust to both threshold value and method.
- Patterns form by a coalesce, split, relax cascade, and spot spacing continues
  to widen and regularize over timescales roughly 100x longer than formation.
