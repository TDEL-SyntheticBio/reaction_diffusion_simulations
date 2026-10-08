# CLAUDE.md — JAPI / PAPI reaction-diffusion engine: conventions

Fork of the simulation code for the JAPI patterning paper (in revision). The engine is unchanged
from upstream commit `07a328b`. Everything below is how the code behaves, verified against the
source on 2026-10-07 and re-checked by independent passes; each item carries a `file:line`. Do not
take the handoff documents' descriptions of the engine on trust. Results, tables and findings
live in `analysis/REPORT_2026-10-07.md`; add new findings there, not here.

`handoff/` holds the analysis handoff (its engine-convention list is **section 8**, not 2), the
October session notes, Supplementary Notes 1–4 and the Supp Table 1/2 per-panel keys, whose
column names are the code parameter names. A blank `nucleation_rate` there means 0.0 (see §9).

## 1. Model as implemented (2D, `2D_simulations/simulation_2D.py`)

Forward Euler in `_vectorized_step` (`:145-195`); the two-pair engine repeats the formulas in
`_pair_hill` (`:309-345`) and `_vectorized_step_two_pairs_cross` (`:348-411`), so any change to
the Hill form must be made at both `:181` and `:344`.

```
hill  = (A + basal) / (A + I + 1 + basal),   A = (act_signal/ka)^n_a,  I = (inh/ki)^n_i
a    += dt * (act_prod_rate * hill - act_decay_rate * a) [+ act_diffusion * dt/dx^2 * lap(a)  if paracrine]
i    += dt * (inh_prod_rate * hill - inh_decay_rate * i) + inh_diffusion * dt/dx^2 * lap(i)
```

- juxtacrine: `act_signal` = mean of `a` over the *valid* neighbours, self excluded (`:173`);
  edge cells average over fewer neighbours. `act_diffusion` is ignored (`:192-193`).
- paracrine: `act_signal = a` and the activator diffuses with the same raw-sum Laplacian (`:175,189`).
- `activator_type` is matched by exact string: 2D tests `== "juxtacrine"` for the signal and
  `== "paracrine"` for diffusion, so any other string gives a self-sensing, non-diffusing hybrid;
  1D treats anything but `"paracrine"` as juxtacrine. Use only the two exact strings.
- `basal_prod` must be in the parameter dict (`:164` indexes it; the solver defaults it).
- `dt = 0.01`, `dx = 1.0` everywhere: **t = steps × 0.01**. `dt` also scales the nucleation
  probability (`rate*dt`) and the noise (`sqrt(dt)`).
- Order per iteration: Euler update → nucleation kick → noise + clip (activator only, `:287`) →
  every `save_every` steps: save + stopping check. Hill terms are zeroed for non-positive
  signals (`:179-180`); nothing else is clipped.
- Two-pair engine: `p1.cross_inhibition_rate * hill2` is added to pair 1's inhibitor production
  and vice versa (`:375-386`); per-pair steady states ignore it, so §10 does not apply with it.

## 2. Hex lattice: even-r offset, zero-flux

`_build_hex_neighbor_arrays` (`:108-142`): **even rows are shifted half a cell to the right**.
Even-row cell (r,c) has upper neighbours (r−1,c),(r−1,c+1); odd-row cell (r−1,c−1),(r−1,c).
Physical centres: `x = c + 0.5*(r % 2 == 0)`, `y = r*sqrt(3)/2`, unit centre-to-centre spacing
(row pitch √3/2, cell area √3/2). Out-of-bounds neighbours are masked; `nbr_count` (6 interior,
2–5 at edges) divides the juxtacrine signal and multiplies the Laplacian self term, so at every
cell `lap = nbr_count * (mean_of_valid_neighbours − u)`: zero flux.

**Parity is wrong elsewhere in the deposit.** `visualize_2D._hex_polygons` (`visualize_2D.py:48`)
shifts odd rows, so every PNG/MP4 draws the adjacency of the column-reversed field;
`_add_overlay_hex_field_rgb` (`:350-355`) uses a third, non-tiling flat-top geometry.
`res_analysis/npz_feature_distribution.py get_neighbors` (`:69-109`) returns the engine's mirror
under its default `--layout even-r` (its `odd-r` branch matches the engine), and its
`approx_diameter_hexagons` is in hex side lengths. `radial_autocor_npz.py`, `crosscorr_npz.py`,
`fft_analysis_2D.py`, `npz_overlap_analysis.py` treat the index array as a square lattice.
Build adjacency from `_build_hex_neighbor_arrays` and distances from the coordinates above;
render with `analysis/plot_hex.py`.

## 3. Laplacian prefactor and inhibitor range

`lap_i = sum(valid neighbours) − nbr_count * i` (`:186`) times `inh_diffusion*dt/dx**2`, with no
hex geometric factor. The raw six-neighbour sum is **1.5× the continuum Laplacian** (the
continuum-matched form would be (2/3)× this, self-coefficient 4). Symbols: Λ_hex(q) =
6 − Σ_j cos(q·n_j) ∈ [0, 9], K̂_hex = (1/6) Σ_j cos(q·n_j) ∈ [−1/2, 1], and **Λ_hex = 6(1 − K̂_hex)**
exactly in the bulk (1D: Λ = 2(1 − cos q), prefactor 1). The dispersion relation written with the
lattice symbol and λ² = D_i/γ is consistent, but **the inhibitor range in cells is √(1.5·D_i/γ)
in 2D** (5.48 at D_i = 10, γ = 0.5; 7.75 at D_i = 20) and √(D_i/γ) in 1D; both confirmed on
arrested-domain halos to 0.1%. The handoff's 2D λ_i values are √1.5 too small.

## 4. Forward-Euler stability and positivity

Most oscillatory inhibitor mode multiplier 1 − dt(γ + 9·D_i): stable iff **D_i < (2/dt − γ)/9 =
22.2** at dt = 0.01, γ = 0.5 (D_i = 22 bounded, 23 blows up; same bound for paracrine
`act_diffusion`; 1D ≈ 50). Reaction stiffness at high b_i or n_i shrinks the margin. Blow-up
shows as a sign-alternating inhibitor of growing magnitude while `np.isfinite` stays True: test
max|i| and negativity, not finiteness. The stencil is non-monotone for dt(γ + 6·D_i) > 1 (true at
D_i = 20), so the inhibitor can go transiently negative and is never clipped.

## 5. Parameters and nondimensionalisation

Every committed parameter/config file has `act_decay_rate = act_half_sat = inh_half_sat = 1.0`,
`basal_prod = 0.0`, so `r_a = act_prod_rate`, `r_i = inh_prod_rate`, `γ = inh_decay_rate`,
`D_i = inh_diffusion` (lattice units) in Supplementary Note 3's notation. `1D_simuations/tests.py`
is stale and unrunnable; do not use it as a reference. The Hill function differs between the
engines (`(A+basal)/(A+I+1+basal)`) and both copies of `finding_steady_states.hill_with_grads`
(`basal + A/(basal+1+A+I)`) whenever basal ≠ 0; never set basal ≠ 0 without reconciling them.
The two solver copies (1D and 2D) share a module name: run 1D and 2D code in separate processes.

## 6. Steady-state solver and the nucleation kick

`run_coupled_hex` (`:228-234`) calls `fast_stable_steady_state`; if it raises or returns a
non-positive or non-finite `a_ss` **or `i_ss`**, the engine sets `a_ss = i_ss = spike_value`,
which also feeds `spike_steady_state` / `random_tight` inits and the a_ss/i_ss metadata. The
nucleation kick is `2*a_ss`: **2·spike_value for every committed b_a = 5 panel** (no activated
state; the JA control hits the fallback through i_ss = 0), 2·a_ss where a state exists
(Fig 3B labyrinth). The solver returns "none" in three situations: no activated state, a
trace-unstable activated state (uniform-ON regime), or a saturated state with H* in (1−1e-9, 1);
its no-root detection is fragile at large b_i/γ. It is not a regime classifier; refine roots
exactly for closed-form work (it rounds to 3 dp); always log a_ss/i_ss.

## 7. Nucleation and randomness

- Each step every cell fires with probability `nucleation_rate*dt` (`:263,281`); a firing cell
  gets `a += 2*a_ss` only if `a < a_ss` (`:282-283`). Nucleation and noise act on the activator.
- Nucleation, noise and the 2D init modes `random_uniform_over0`, `random_tight`, `near_zero`
  use the **unseeded** numpy global RNG; only the `activator_random_spikes` indices are seeded
  (`init_seed = 2`, hard-coded). `batch_runner_2D.py:42` seeds Python's `random`, which nothing
  uses. No 2D run in the deposit is reproducible. `np.random.seed(s)` immediately before
  `run_coupled_hex` makes a single-process run reproducible (`analysis/engine_wrap.py` does this).
- 1D seeds `default_rng(2)` with the seed hard-coded (`simulation.py:179`): every 1D random run
  is the same realisation, and 1D has no nucleation or noise.
- 2D init modes: `all_off`, `random_uniform_over0` (a ~ U(0, spike_value)), `spike_steady_state`
  (centre cell at (a_ss, i_ss), both = spike_value under the fallback), `activator_random_spikes`,
  `both_on`, `random_tight`, `near_zero`, `inhibitor_spike`. **2D has no `activator_spike`**; 1D has it.

## 8. Stopping criterion and frame bookkeeping

Every `save_every` steps (loop indices 0, 200, …, each **after** the update) a frame is appended
and the run stops when `step > min_steps` and the L1 change per cell-field since the previous
frame is below `stopping_threshold` (`:289-298`; 1e-4 committed). `history[0]` is the initial
state; `history[k≥1]` is the state after (k−1)·save_every + 1 Euler steps. **On a full-length run
the last frame is not the final state** (step 4800 for a 5,000-step run) although `main_2D.py:174`
labels it `final_step_004999` and both batch runners write it as final; frame labels in
`main_2D.py --step-size` and `animate_histories` overstate by save_every − 1. Use
steps = n·save_every + 1 (as `engine_wrap.run` does) or save the final state yourself. Under
nucleation the criterion never fires; without it, all-on and all-off panels stop early. Always
record `steps_used`; elapsed time is (returned step + 1)·dt.

## 9. Entry-point traps

- `batch_runner_2D.py:92` defaults `nucleation_rate` to 0.02 and the committed
  `2D_batch/config.yaml` has no such key, so the Fig-1C-like template nucleates. Engine defaults
  differ (`run_coupled_hex` 0.0, `run_two_pair_hex` 0.01). **Always pass `nucleation_rate`, including 0.0.**
- `batch_coupled_2D.py:64-80` never forwards `nucleation_rate`, `noise_amplitude` or
  `set_peak_height`; `config_coupled.yaml:17` is dead. Both batch runners default `init_mode` to
  the nonexistent `"activator_spike"`, `min_steps` to 10,000, `spike_value` to 5.0.
- npz keys differ: batch `A_final`/`R_final`; `main_2D.py` `a`/`i`; coupled `a1,i1,a2,i2`.
- Methods vs code: the kick is `2*a_ss` (Methods: a_ss or spike_value); Fig 1C amplitude is
  U(0, 1) per Supp Table 1 (Methods: U(0, 2)); 1D `activator_spike` rows use 2.0.

## 10. Closed forms (basal = 0, μ_a = 1, κ = 1, single pair)

With P = A/(1+A+I), Q = I/(1+A+I) at an activated state (A = a0^n_a, I = i0^n_i):
α = n_a(1 − P), ν = γ·n_i·Q, G = ν/γ (exact; checked against `hill_with_grads`).
Δ(q)/γ = (1 + λ²Λ(q))(1 − αK̂(q)) + G with λ² = D_i/γ and the lattice symbols of §3.
q = 0: G > α − 1 (saddle-node at equality) and τ(0) = α − 1 − γ(1 + G) < 0. Finite q with
Λ = p(1 − K̂), p = 2 (1D) or 6 (hex), L = pλ²: the minimum is interior iff λ² > α/(p(α−1))
(necessary only); instability iff **G < [α(1+L) + L]²/(4αL) − (1 + L)**. 1D and hex windows
differ, so never reuse 1D classifier output for 2D. The determinant-minimum mode is not the
fastest-growing one; use the latter for wavelengths. Paracrine: K̂ ≡ 1 and Supp Note 1 Eq. 43 applies.

## 11. Deposit analysis code and environment

No regime classifier exists in the repo; `1D_batch/analyze_patterns.py` is 1D-only and emits
a True/False/NaN `pattern_flag`. The `res_analysis/` scripts carry the §2 geometry caveats.
Python 3.13; `scipy`, `matplotlib`, `joblib`, `tqdm` must be pip-installed each session.
100×100 juxtacrine: ≈ 8 s per 5,000 steps on one core; 4 cores. `results/`, `runs/`, `*.npz`,
`*.mp4`, `*.png` are gitignored except `analysis/summaries/*.png`.

## 12. Working and measurement rules

- Wrap the engine; never edit `simulation_2D.py` / `simulation.py` to reproduce panels.
  Fixes (seeding, final frame, parity) go in `analysis/`, opt-in, recorded per run.
- Report seed count, grid, `steps_used`, threshold and size cut with every number. Run
  `python analysis/validate_pipeline.py` (asserted ground-truth checks) before trusting a new number.
- Segmentation: variance gate, then level. **Under continuous nucleation never threshold at half
  the field maximum** (it follows the newest kick and cut through every domain at b_i = 12): use a
  fixed level 0.3·act_prod_rate with persistence across two frames t = 10 apart. Components are
  bimodal (single-cell transients and domains); counting transients as spots is what produced the
  October notes' near-Poisson CV. Set the size cut from the measured domain-size distribution.
- Spacing CV: compare with a same-n, same-window Poisson reference (`pointstats.csr_reference_cv`),
  not 0.52. Pool g(r) by summing counts; report hole/peak with bootstrap SD and pair counts; quote
  peaks only when significant; use the lattice-site null for components under ~10 cells.
- Quote inhibitor ranges as √(1.5·D_i/γ) in 2D and √(D_i/γ) in 1D, and say which.
- "Arrested" needs the isolated-seed test run to t ≥ 600 (`analysis/isolated_seed_fate.py`);
  a domain constant in area for 40 time units can still divide later.
- Use the fastest-growing mode for wavelength predictions.

## 13. Analysis package (`analysis/`)

| module | purpose |
|---|---|
| `hexgeom.py` | coordinates, `label_components` on the engine's adjacency, `component_table`, `segment`, `rasterize_discs`, `equivalent_radius` |
| `pointstats.py` | `nn_cv`, `border_cv`, `csr_reference_cv`, `PooledG` (null-normalised g(r), bootstrap SDs, gated peaks) |
| `synthetic.py`, `validate_pipeline.py` | known-answer patterns and the asserted validation suite (includes the engine fate test) |
| `engine_wrap.py` | seeded, protocol-explicit `run()`; last frame is the final state |
| `run_panels.py`, `measure_panels.py`, `domain_tracking.py` | published sets under both protocols; measurements; lifetimes |
| `single_domain_halo.py`, `isolated_seed_fate.py`, `arrest_onset_sweep.py`, `nucleus_threshold.py` | isolated-domain experiments |
| `plot_hex.py` | the only hex-correct renderer in the repo |

Raw runs go to `results/` (gitignored); summary CSVs and figures to `analysis/summaries/`.
