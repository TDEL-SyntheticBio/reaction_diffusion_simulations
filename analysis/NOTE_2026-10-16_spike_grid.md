# NOTE 2026-10-16: 1D single-spike arrest map over the (n_a, n_i) plane

Record of what was run and what came out. No framing; the arrest map is meant to be put beside the Turing map.
Everything here is 1D (`1D_simuations/simulation.py` conventions: juxtacrine signal = mean of the two neighbours,
Neumann Laplacian with prefactor 1, inhibitor range sqrt(D_i / gamma) = 10, 7.07, 5.77, 5, 4.47, 4.08, 3.78, 3.54, 3.33,
3.16 cells at gamma 0.1 to 1.0 and 2.24 to 1.0 cells at gamma 2 to 10). No 2D results or conventions are used.

## Files

| file | content |
|---|---|
| `analysis/spike_grid.py` | driver: audit, bit-exact integrator replica, labelling, validation, spike-2 subsample, aggregation |
| `analysis/summaries/state_audit_turing.csv` | 273,600 sets: `truth_has_state`, `truth_a_star`, `truth_i_star`, `truth_H`, `n_roots`, `trace0`, `det0`, `turing`, `q_fastest`, `sigma_max`, `wavelength_fastest`, `q_unstable_min/max`, `N_modes` (= 200) |
| `analysis/summaries/spike_grid_full.csv` | one row per set: parameters, `step_returned` (the engine's loop index at the stop), `steps_used` (= step + 1 Euler steps), `hit_cap`, `label`, `n_domains`, `largest_domain`, `cells_above`, `fwhm`, `domain_start`, `a_min/a_max/a_mean`, `i_mean/i_max`, `dom_lambda`, `score`, `amp_sd`, `pattern_flag`, the audit columns, `lambda_i` and the sizes divided by it |
| `analysis/summaries/spike_grid_by_hill.csv` | per (n_a, n_i): counts of the four outcomes, capped runs, audit counts (`n_truth_has_state`, `n_turing`), `pattern_flag` count, one-domain size median / quartiles / max in cells and median / quartiles in lambda_i, the > 150-cell and <= 5-cell counts, and the cross-counts between labels and the audit |
| `analysis/summaries/spike_grid_validation.csv`, `spike_grid_validation_by_hill.csv` | the random_tight validation sample and its per-cell disagreement rates |
| `analysis/summaries/spike_height_2.csv` | the spike = 2.0 subsample |
| `results/spike_grid/full/cell_na??_ni??.{npz,csv}` | per-cell checkpoints (gitignored): final activator and inhibitor profiles of every set as float32 (`a_final`, `i_final`, 1,900 x 200 each), the stop step, and the per-set measurement rows. 115 MB of compressed npz plus 57 MB of per-cell CSV for 144 cells (165 MiB on disk). Regenerate with `python analysis/spike_grid.py run 4` (1.6 h on 4 cores; finished cells are skipped) and `python analysis/spike_grid.py aggregate`. |

## Conventions

- Grid exactly as the committed `1D_batch/config.yaml` sweep: b_a, b_i = linspace(1, 10, 10); gamma = linspace(0.1, 1, 10) +
  linspace(2, 10, 9), 19 values; n_a, n_i = 1..12; 273,600 sets. **The gamma axis was not reduced.**
- Fixed: N = 200, steps = 50,000, dt = 0.01, dx = 1, save_every = 200, stopping_threshold = 1e-4, min_steps = 500, juxtacrine,
  act_decay = act_half_sat = inh_half_sat = 1, basal = 0, D_i = 10, act_diffusion = 1 (ignored), spike 10 at cell N // 2 = 100,
  inhibitor zero. Hill coefficients are passed as floats, as `make_param_grid` passes swept values.
- Integration: a replica of `run_coupled_neumann`'s loop vectorised over the 1,900 sets of one (n_a, n_i) cell, same expression
  order, same stopping rule (every 200 steps, once step > 500, stop when the L1 change per cell-field since the previous save is
  below 1e-4). It returns the state at the true final step. `run_simulation` is not called; the deposited
  `fast_stable_steady_state` is never called (activator_spike needs no steady state). Bit-exact check: 30 runs (15 at N = 100,
  15 at N = 200, four of them at n_a = n_i = 3) against `run_simulation`, comparing the engine's frame and step: 30 of 30
  identical. `steps_used` in the CSV is the engine's `step` + 1; `step_returned` is the engine's value.
- Labels from contiguous runs of cells with a > 0.3 b_a on the final activator: 0 runs "all off"; 1 run of >= 190 cells
  (0.95 N) "all on"; 1 run "one domain"; >= 2 runs "multiple domains". No peak finding. `fwhm` is the length of the contiguous
  run containing the maximum with a >= a_max / 2.
- `analyze_pattern` is the deposited function, executed from the source of `1D_batch/analyze_patterns.py` (that module runs a
  CSV script at import time, so the function is extracted by `ast` and executed verbatim), called as the deposit calls it:
  `analyze_pattern(a, dx=1.0, plot=False)`.
- Audit: H* = largest root above 1e-5 of H = A / (A + I + 1) with A = (b_a H)^n_a, I = (b_i H / gamma)^n_i, from a scan on
  20,001 uniform H points plus a geometric approach to H = 1 (1 - 1e-6 to 1 - 1e-16, for saturated states at high n_a) and 60
  bisections; a* = b_a H*, i* = (b_i / gamma) H*. Stability from the 2x2 reaction Jacobian (trace < 0 and det > 0).
  Turing from the dispersion relation over the accessible modes q = pi m / N, m = 0..199, K = cos q, Lambda = 2 (1 - cos q):
  unstable iff det < 0 or trace > 0; Turing = stable at q = 0 and unstable at some q > 0. `sigma_max` is the largest real part
  of the unstable eigenvalue over those modes. The closed-form inequality is not used anywhere.
- Cross-check at n_a 10, n_i 4, N = 100 (all-on rule >= 95 cells): all off 390, all on 1293, one domain 217, multiple 0,
  exactly as required.

## Runtime

- 30 single runs at N = 200 (8 of them at n_a = 3): mean 0.37 s, max 2.26 s, one hit the cap. Naive projection
  273,600 x 0.37 s / 4 cores = 7.1 h (<= 8 h, so the full grid was run). The batched integrator is faster: the Turing-heavy
  cell n_a = n_i = 3 (1,900 sets, 158 capped) took 190 s, projecting 1.9 h; the full grid took **5,679 s = 1.58 h** on 4 workers
  (144 cells, one checkpoint per cell, 1,900 sets each).
- Validation sample and spike-2 subsample: a few minutes each.

## Results, full grid (273,600 sets, N = 200)

| label | sets | hit the 50,000-step cap | truth_has_state | turing |
|---|---|---|---|---|
| all on | 148,785 | 205 | 148,785 | 21 |
| all off | 90,224 | 2,811 | 21,808 | 1,999 |
| one domain | 24,298 | 1,222 | 7,214 | 1,377 |
| multiple domains | 10,293 | 8,775 | 7,447 | 7,300 |

- Every "all on" set has a stable activated state in the audit. 2,846 "multiple domains" sets have neither a stable state nor a
  Turing call; 2,993 are not Turing.
- The deposited `pattern_flag` is True for 24,184 of 24,298 one-domain sets, 10,282 of 10,293 multiple-domain sets, 8,679
  all-off sets and 2,797 all-on sets.
- Per-cell one-domain counts (rows n_a, columns n_i; the full table with quartiles is `spike_grid_by_hill.csv`):

| n_a \ n_i | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2 | 105 | 65 | 20 | 12 | 10 | 6 | 3 | 2 | 1 | 1 | 0 | 0 |
| 3 | 112 | 153 | 191 | 164 | 89 | 55 | 39 | 39 | 34 | 31 | 27 | 26 |
| 4 | 75 | 179 | 186 | 224 | 231 | 204 | 147 | 106 | 83 | 65 | 52 | 47 |
| 5 | 67 | 148 | 183 | 202 | 242 | 258 | 252 | 228 | 188 | 145 | 118 | 98 |
| 6 | 56 | 113 | 202 | 206 | 222 | 239 | 267 | 273 | 271 | 248 | 219 | 179 |
| 7 | 37 | 90 | 162 | 221 | 215 | 229 | 247 | 267 | 278 | 284 | 276 | 263 |
| 8 | 24 | 76 | 135 | 197 | 226 | 224 | 237 | 258 | 268 | 278 | 291 | 292 |
| 9 | 14 | 81 | 136 | 208 | 266 | 293 | 300 | 303 | 319 | 332 | 345 | 358 |
| 10 | 13 | 88 | 141 | 203 | 258 | 294 | 320 | 316 | 322 | 336 | 343 | 355 |
| 11 | 17 | 87 | 140 | 192 | 245 | 291 | 310 | 337 | 322 | 330 | 341 | 347 |
| 12 | 21 | 87 | 141 | 175 | 231 | 275 | 305 | 325 | 341 | 332 | 335 | 343 |

- Turing counts per cell from the audit (same layout): zero at n_a = 1; n_a = 2 row 18, 127, 171, 164, 152, 134, 120, 113, 99,
  95, 87, 82; n_a = 3 row 2, 47, 134, 190, 214, 221, 212, 199, 188, 172, 163, 142; falling with n_a to the n_a = 12 row 0, 0, 0,
  3, 1, 7, 3, 7, 22, 23, 33, 51. Multiple-domain counts per cell are in the CSV.
- The measured cell n_a 10, n_i 4 at N = 200: all on 1,308, all off 389, one domain 203, multiple 0 (N = 100: 1,293 / 390 / 217 / 0).

## Caveats to carry

- In cell n_a = n_i = 1 ten sets have truth_a_star equal to 0.3 b_a to 1e-15, so their uniform activated state sits exactly on
  the segmentation level; the 1 one-domain and 7 multiple-domain labels in that cell (and two of its all-on / all-off labels)
  are decided by residuals below 1e-4 of a field that is uniform to that precision. They reproduce from the float64 engine but
  are not domains in any physical sense.

- **Boundary truncation is not cleared at N = 200.** Of the 24,298 one-domain sets, 590 have a domain above 150 cells
  (maximum 189; 989 are >= 100 cells; 399 are 101-150). At the measured cell the largest one-domain size is 188 at N = 200
  versus 94 at N = 100: the same sets reach the box again. By gamma: 106 of the 1,019 one-domain sets at gamma 0.1 exceed
  150 cells, 93 of 2,074 at 0.2, 72 at 0.3, falling to 2 at gamma 5 and 0 at gamma >= 6. Of the 1,222 one-domain runs that hit
  the step cap, the median domain is 129 cells (10th-90th percentile 55-180). These sizes are box-limited and should not be
  quoted as intrinsic.
- **The 5-cell floor is a lattice limit and is present at every gamma, not only gamma >= 2.** The one-domain size median is 5
  cells at every gamma except 0.1 (median 9); 17,101 of 24,298 one-domain sets are <= 5 cells, 14,210 of the 21,074 at gamma <= 1
  and 2,891 of the 3,224 at gamma >= 2. At gamma >= 6 the maximum one-domain size is 5 cells. N = 200 does not change this.
- Size in units of the 1D inhibitor range sqrt(D_i / gamma): one-domain median 1.22 lambda_i (quartiles 0.73-2.2); by gamma the
  median rises from 0.7-0.9 lambda_i at gamma 0.1-0.3 to 1.58 at gamma 1 and 5.0 at gamma 10, where lambda_i is 1 cell and the
  floor dominates.
- 2,811 "all off" sets hit the step cap with a sub-threshold activator still evolving (a_max between 0.12 and 1.0 of the 0.3 b_a
  level, median 0.65; 267 of them below 0.3); 1,637 of them are Turing sets and they sit at n_a 2-5 (1,409 / 1,015 / 343 / 44). "All off" there means
  below the 0.3 b_a level, not zero.

## Validation sample (random_tight around truth_a_star)

Sample: per (n_a, n_i) cell, up to 9 Turing sets and up to 9 stable non-Turing sets drawn at random (seed 0) from the audit:
2,225 sets (929 Turing, 1,296 stable non-Turing; every cell holds 9-18: the 29 cells with no Turing set in the audit, all of n_a = 1 and 17 cells at low n_i, hold only their 9 stable picks). Protocol: the engine's random_tight draw (seed 2,
the same uniform stream for every set) scaled to +-5 % of (truth_a_star, truth_i_star) in both fields, N = 200, the run
settings above (stopping rule on). Outcome from the final activator relative to a*: "patterned" if max - min > 0.05 a*,
"off" if the range is smaller and the mean is below a* / 2, "uniform" otherwise. Agreement: Turing sets should be patterned,
stable non-Turing sets uniform. Every first-pass disagreement was rerun for the full 50,000 steps (t = 500, stopping rule off)
at +-5 % and at +-0.5 %.

- Stable non-Turing sets: 1,296 of 1,296 stayed uniform at +-5 % (range below 0.03 a*; mean within 3.1e-5 of a*, 8 beyond 1e-5, all at n_a = 1). **No finite-amplitude collapse
  occurred in this sample** (the measured cell's 2 of 1,494 did not recur here).
- Turing sets: 929; first-pass disagreements 34. Eleven were slow modes (sigma_max 0.002-0.010) at which the stopping rule
  fired at step 1,600-25,400 with the amplitude at 0.025-0.050 a* (one at 4e-4 a*); on the full-length run they are patterned,
  so they agree. Twenty-two are Turing-unstable states that fell to the all-off state (activator mean / a* at t = 500 between 1e-205 and 3e-160 from +-5 %
  and between 7e-201 and 3e-93 from +-0.5 %): all at b_a = 2 (one at b_a = 3), a* = 1.3-2.6, spread over the cells n_a 6-12 listed below. The linear
  call "unstable" holds for them; the outcome is not a pattern. One set (n_a 3, n_i 10, b_a 8, b_i 6, gamma 0.7, sigma_max =
  0.0020) stayed uniform through t = 500 at both amplitudes; its expected amplification over the run is exp(sigma t) = 2.7,
  so the run is too short to decide it.
- Disagreement rate after the retests, counting the fell-to-off sets as disagreements: 23 of 2,225 (1.0 %), in 14 of the
  144 cells; zero in the other 130. Per cell (n_a, n_i: disagreements / sample): 3,10: 1/18 (the slow set); 6,5: 1/18;
  6,6: 1/18; 7,3: 1/15; 8,1: 2/11; 8,5: 3/18; 9,2: 1/10; 9,4: 3/14; 9,6: 3/18; 10,7: 3/18; 11,6: 1/14; 11,8: 1/18;
  12,4: 1/12; 12,9: 1/18. The full per-cell table is `spike_grid_validation_by_hill.csv`; per-set rows with both retests and a
  `verdict` column are in `spike_grid_validation.csv`.
- Counting the fell-to-off sets as agreeing with "unstable" (which is what the audit asserts), the analytic calls disagree
  with simulation on 1 of 2,225 sets, the slow one above.

## Spike height 2.0 (side question)

Subsample b_a, b_i in {1, 3, 5, 7, 10}, gamma in {0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 5, 10}: 225 sets per (n_a, n_i) cell, 32,400
runs at N = 200, spike 2.0 at the centre cell, everything else as above.

| n_a | all off | all on | multiple | one domain | smallest b_a with any response |
|---|---|---|---|---|---|
| 1 | 1,728 | 971 | 1 | 0 | 3 |
| 2 | 1,560 | 1,007 | 100 | 33 | 3 |
| 3 | 1,602 | 832 | 191 | 75 | 5 |
| 4 | 1,839 | 599 | 173 | 89 | 7 |
| 5 | 1,887 | 554 | 153 | 106 | 7 |
| 6 | 2,230 | 351 | 58 | 61 | 10 |
| 7 | 2,264 | 333 | 39 | 64 | 10 |
| 8-12 | 2,700 each | 0 | 0 | 0 | none |

Spike 2.0 is dead on the whole subsample for n_a >= 8 (13,500 sets, all "all off") and responds at n_a <= 7 only above a
b_a floor that rises with n_a (b_a >= 3 at n_a 1-2, >= 5 at n_a 3, >= 7 at n_a 4-5, >= 10 at n_a 6-7). Per-cell counts are in
`spike_height_2.csv`.
