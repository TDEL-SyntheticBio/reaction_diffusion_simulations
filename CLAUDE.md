# CLAUDE.md — JAPI / PAPI reaction-diffusion engine

Fork of the simulation code for the JAPI patterning paper (in revision). The engine
is unchanged from upstream commit `07a328b`. Everything below was verified against
the source on 2026-10-07 by reading the code and running numeric checks; each item
carries a `file:line` so it can be re-verified. Do not take the handoff documents'
descriptions of the engine on trust; take these.

Supporting material lives in `handoff/`:

- `JAPI_irregular_handoff_2026-10-07.md` — regime map, g(r) methodology, results, open items
- `JAPI_session_findings_Oct2026.md` — earlier session (asynchronous nucleation result)
- `SupplementaryNote1..4_*.md` — the paper's supplementary notes (LSA, velocity, nondimensionalisation, arrest)
- `supp_table1_panel_parameters.csv`, `supp_table2_coupled_panel_parameters.csv` — per-panel reproduction keys

## 1. Model as implemented (2D, `2D_simulations/simulation_2D.py`)

Forward Euler, one step per call of `_vectorized_step` (`simulation_2D.py:145-195`):

```
hill  = (A + basal) / (A + I + 1 + basal),   A = (act_signal/ka)^n_a,  I = (inh/ki)^n_i
a    += dt * (act_prod_rate * hill - act_decay_rate * a) [+ act_diffusion * dt/dx^2 * lap(a)  if paracrine]
i    += dt * (inh_prod_rate * hill - inh_decay_rate * i) + inh_diffusion * dt/dx^2 * lap(i)
```

- **juxtacrine** `act_signal` = mean of `a` over the *valid* neighbours, self excluded
  (`simulation_2D.py:173`). Interior weight 1/6 per neighbour; edge cells average over
  fewer neighbours (`nbr_count`), never over self.
- **paracrine** `act_signal = a` (self-sensing) and the activator diffuses with the same
  Laplacian (`simulation_2D.py:175,189-191`).
- `dt = 0.01`, `dx = 1.0` everywhere. **Time unit: t = steps × 0.01**, so 5,000 steps = t 50,
  50,000 steps = t 500.
- Negative values are not clipped except after additive noise (`simulation_2D.py:287`).
  Hill terms are zeroed for non-positive signals (`simulation_2D.py:179-180`).

The two-pair engine (`run_two_pair_hex`, `simulation_2D.py:414-580`) uses identical
conventions; cross-inhibition adds `p1.cross_inhibition_rate * hill2` to pair 1's inhibitor
production and vice versa (`simulation_2D.py:375-386`).

## 2. Hex lattice: even-r offset, zero-flux

`_build_hex_neighbor_arrays` (`simulation_2D.py:108-142`):

- **even-r**: *even* rows are shifted half a cell to the **right**. Even-row cell (r,c) has
  upper neighbours (r-1,c),(r-1,c+1); odd-row cell has (r-1,c-1),(r-1,c). Verified: 0/168
  cells disagree with an explicit even-r coordinate construction; 168/168 disagree with odd-r.
- Physical centre coordinates for any new plotting or distance measurement:
  `x = c + 0.5*(r % 2 == 0)`, `y = r*sqrt(3)/2` (unit centre-to-centre spacing).
- **Zero-flux boundaries**: out-of-bounds neighbours are masked; `nbr_count` = number of
  valid neighbours (6 interior; 2,3,4,5 at edges and corners) and divides the juxtacrine
  signal and multiplies the self term of the Laplacian. Neighbour relation is symmetric.

### Visualiser parity bug (do not trust rendered adjacency)

`visualize_2D._hex_polygons` (`visualize_2D.py:46-50`) shifts **odd** rows right (odd-r)
although its docstring (`:37`) says even-r. Published hex renderings therefore draw each
odd row displaced by one cell relative to the simulated adjacency. `_add_overlay_hex_field_rgb`
(`visualize_2D.py:350-355`) uses a third, flat-top geometry (orientation 30°, dx = 1.5R,
dy = √3R) that is not a consistent tiling of the same lattice. Fix the plotting, never the
engine, and use the coordinates in §2 when you do.

## 3. Laplacian: raw sum, prefactor 6, effective D is 1.5× nominal

`lap_i = sum(valid neighbours) - nbr_count * i` (`simulation_2D.py:186`), scaled by
`inh_diffusion * dt / dx**2`. This is the raw neighbour sum, **not** the continuum-matched
hex Laplacian (which would be (2/3)× this, self-coefficient 4). Consequences, all verified
numerically:

| quantity | value |
|---|---|
| raw-sum Laplacian of x²+y² (continuum value 4) | 6.000 → prefactor **1.5** |
| Fourier symbol Λ_hex(q) = 6 − Σ_j cos(q·n_j) | range [0, 9]; Λ_hex ≈ 1.5\|q\|² at small q |
| juxtacrine kernel symbol K̂_hex(q) = (1/6) Σ_j cos(q·n_j) | range [−1/2, 1] |
| identity | **Λ_hex = 6 (1 − K̂_hex)** exactly (same neighbour set); in 1D Λ = 2 (1 − cos q) |
| inhibitor decay length on the lattice, D_i=10, γ=0.5 | fitted **5.5 cells** = √(1.5·D_i/γ); nominal √(D_i/γ) = 4.5 |
| same, D_i=20 | fitted **7.9 cells**; nominal 6.3 |
| same, D_i=5 | fitted 3.9 cells; nominal 3.2 |

So: the dispersion relation written with the lattice symbol, Δ/γ = (1 + λ²Λ(q))(1 − αK̂(q)) + G
with λ² = D_i/γ, is internally consistent (the 1.5 lives inside Λ_hex). But **any statement
of the inhibitor range in cells must use √(1.5·D_i/γ)**, not √(D_i/γ). The handoff's
"λ_i = 4.5 cells (D_i=10) / 6.3 cells (D_i=20)" understates the engine's range by 22%.

The 1D engine (`1D_simuations/simulation.py:101-112`) has prefactor 1 (two neighbours), so
√(D_i/γ) is correct in 1D.

## 4. Forward-Euler stability limit

Most oscillatory inhibitor mode multiplier = 1 − dt(γ + 9 D_i). At dt = 0.01, γ = 0.5 the
bound is **D_i < 22.2**. Verified empirically: D_i = 22 stays bounded, **D_i = 23 blows up**
(inhibitor max 2×10⁹ in 400 steps). The committed D_i = 20 (Fig 2L/3/S8D) has ~10% margin.
Any sweep of `inh_diffusion` above 22 at dt = 0.01 is numerically meaningless; reduce dt.
In 1D the bound is D_i < (2/dt − γ)/4 ≈ 50.

## 5. Parameters and nondimensionalisation

In **every** committed parameter/config file (`2D_simulations/parameters_2D.py`,
`parameters_coupled_2D.py` p1 and p2, `2D_batch/config.yaml`, `2D_batch/config_coupled.yaml`
p1 and p2, `1D_simuations/parameters.py`, `1D_batch/config.yaml`):

```
act_decay_rate = act_half_sat = inh_half_sat = 1.0,   basal_prod = 0.0
```

so, in the notation of Supplementary Note 3: `r_a = act_prod_rate`, `r_i = inh_prod_rate`,
`γ = inh_decay_rate`, `n_a = act_hill_coeff`, `n_i = inh_hill_coeff`, `D_i = inh_diffusion`
(lattice units, see §3). Supp Table 1 column names are the code parameter names.

**Hill function differs between files when basal ≠ 0.** `simulation_2D.py:181` computes
`(A+basal)/(A+I+1+basal)`; `finding_steady_states.hill_with_grads` (`2D_simulations/finding_steady_states.py:17-18`)
computes `basal + A/(basal+1+A+I)` (can exceed 1). Identical at basal = 0 (max diff 2e-16
over 2000 random points), max diff = basal otherwise. Never set basal ≠ 0 without
reconciling the two.

## 6. Steady-state solver and its fallback — this sets the nucleation kick

`run_coupled_hex` (`simulation_2D.py:228-234`) calls `fast_stable_steady_state`; if it
returns a non-positive or non-finite `a_ss` **or `i_ss`**, the engine sets
`a_ss = i_ss = spike_value`. Solver output for the published panels (JAPI, κ = 1):

| panel | b_a | b_i | γ | n_a | n_i | solver a_ss, i_ss | a_ss the engine uses | kick = 2·a_ss |
|---|---|---|---|---|---|---|---|---|
| 1C all on | 5 | 1 | 0.5 | 3 | 3 | 4.658, 1.863 | 4.658 | 9.3 |
| 1C stripes | 5 | 3 | 0.5 | 3 | 3 | 1.706, 2.048 | 1.706 | 3.4 |
| 1C reg / irreg / off spots | 5 | 5,12,14 | 0.5 | 3 | 3 | none | 1.0 (spike_value) | 2.0 |
| 2L JA | 5 | 0 | 0.5 | 10 | 4 | 5.0, **0.0** → fails i_ss>0 | 2.0 (spike_value) | **4.0** |
| 2L JAPI, 3B spots, 3H, 3L, S8D | 5 | 15–50 | 0.5 | 10 | 4 | none | 2.0 (spike_value) | **4.0** |
| 3B labyrinth | 6 | 15 | 0.5 | 10 | 4 | 5.911, 29.55 | 5.911 | **11.8** |
| S5F irregular bistable | 5 | 5 | 0.4 | 6 | 3 | 4.043, 10.11 | 4.043 | 8.1 |
| S5F irregular monostable | 5 | 10 | 0.5 | 6 | 3 | none | 1.0 | 2.0 |

So every experiment-matched panel at b_a = 5 nucleates with a kick of **4.0** (= 2 × spike_value),
not 2 as the October session notes say and not a_ss as Methods say; the Fig 3B labyrinth
panel kicks with 11.8. The JA control reaches the same 4.0 only via the fallback.

## 7. Nucleation and randomness

- Each step every cell fires with probability `nucleation_rate * dt` (`simulation_2D.py:263,281`).
  A firing cell gets `a += 2 * a_ss` **only if its activator is below a_ss** (`:282-283`).
- Nucleation, `noise_amplitude` noise, and the 2D init modes `random_uniform_over0`,
  `random_tight`, `near_zero` all draw from the **unseeded** numpy global RNG
  (`simulation_2D.py:60-99,281,286`). Two calls of `random_uniform_over0` differ. The only
  seeded draw is the spike-index choice for `activator_random_spikes` (`init_seed = 2`, `:236-244`).
  `batch_runner_2D.py:42` seeds Python's `random`, which nothing uses: **no 2D run in the
  deposit is reproducible**, including the Fig 1C initial conditions. New code must seed
  explicitly (pass a `numpy.random.Generator`), never rely on the engine.
- 1D (`1D_simuations/simulation.py:33-41`) seeds its initial condition (`seed=2`) and has
  no nucleation and no noise, so 1D runs are reproducible.
- Init modes (2D): `all_off` both fields 0; `random_uniform_over0` a ~ U(0, spike_value),
  i = 0 (Fig 1C uses spike_value 1 → U(0,1)); `spike_steady_state` centre cell at (a_ss, i_ss);
  `activator_random_spikes` n_points cells at `set_peak_height`; `both_on`; `random_tight`
  (±5% around a_ss, i_ss); `near_zero`; `inhibitor_spike`.

## 8. Stopping criterion

Every `save_every` steps: `diff` = L1 change of a and i since the last saved frame; stop when
`step > min_steps` and `diff/(2·Ny·Nx) < stopping_threshold` (`simulation_2D.py:289-298`).
Committed values: `stopping_threshold = 1e-4`, `save_every = 200`, `min_steps` 500–1000.
Under nucleation at 0.01 on 100×100 the criterion never fires (measured diff ≈ 0.3 per cell per
window), so those runs go the full `max_steps`. Without nucleation a Fig 1C-like run was still at
1.3e-3 after 20,000 steps, so at 1e-4 those runs also tend to run long; the October notes'
"~1e-3" threshold is not what is committed. Always record `steps_used`.

## 9. Analysis code in the deposit

- No regime classifier (Turing / irregular / all-on / all-off) exists anywhere in the repo.
  `1D_simuations/1D_batch/analyze_patterns.py` is 1D-only and emits a binary `pattern_flag`
  (spectral peak/background score > 8 and amplitude check). `2D_simulations/res_analysis/`
  holds FFT, cross-correlation, radial autocorrelation, hex connected components and overlap
  metrics, none of which assign a regime label.
- **`res_analysis/npz_feature_distribution.py` has the hex parity mirrored.** Its `get_neighbors`
  (`:69-100`) labels "even rows shifted right" as `even-r` but the offsets it returns are the
  engine's *odd-r* mirror: with the default `--layout even-r`, 168/168 cells get a neighbour
  set different from the engine; `--layout odd-r` matches the engine on 168/168. Any component
  count, size or diameter produced by the deposit with the default layout used the wrong
  adjacency. To match the engine pass `--layout odd-r`, or build components from
  `simulation_2D._build_hex_neighbor_arrays` directly (preferred).
- `radial_autocor_npz.py`, `crosscorr_npz.py` and `fft_analysis_2D.py` treat the (Ny, Nx)
  index array as a square lattice: no √3/2 row spacing, no half-cell row shift. Radial
  distances along rows are overstated by 15% and the FFT is anisotropically distorted.
  Use the §2 coordinates for any distance-based statistic.

## 10. Closed forms that are safe to use (verified against the code's Hill derivatives)

At any activated homogeneous state with k_a = k_i = μ_a = 1, κ = 1, with
P = A/(1+A+I), Q = I/(1+A+I):

- α ≡ β_a f_a0 = **n_a (1 − P)**, ν ≡ −β_i f_i0 = **γ n_i Q**, G ≡ ν/γ = n_i Q
  (match `hill_with_grads` to 6 decimals on four parameter sets).
- Δ(q)/γ = (1 + λ²Λ(q))(1 − αK̂(q)) + G, λ² = D_i/γ, Λ and K̂ the *lattice* symbols of §3.
- q = 0: determinant condition G > α − 1 (saddle-node at equality); trace condition
  τ(0) = α − 1 − γ(1 + G) < 0.
- Cross-check: these closed forms reproduce the handoff's matched-pair numbers exactly
  (b_a = 20, n = 3, γ = 0.5: τ(0) = −0.006 at b_i = 26, **+0.023 at b_i = 27**, +0.056 at 28), so
  the b_i = 28 member of the matched pair is trace-unstable, not linearly stable, on either
  lattice; the pair straddles the trace boundary. Fig 1C series (b_a = 5): τ(0) turns positive
  between b_i = 3.6 and 3.7 while the finite-q determinant minimum is still deeply negative, and
  the activated state vanishes between 3.7 and 3.8.
- 1D finite-q window opens iff λ² > α/(2(α−1)) (minimum of Δ over c = cos q interior).
  On the hex lattice the same algebra with Λ = 6(1 − K̂) gives λ² > α/(6(α−1)); the 1D
  and hex windows differ, so do not reuse 1D classifier output for 2D claims.

## 11. Known discrepancies between Methods, Supp Table 1 and code

- Nucleation increment: Methods say a_ss or spike_value; code does `2*a_ss` with the §6 fallback.
- Fig 1C amplitude: Methods say U(0,2); Supp Table 1 and `2D_batch/config.yaml` say spike_value 1 → U(0,1).
- Screen spike amplitude for 1D `activator_spike` is 2.0 in Supp Table 1 (Fig 1F); 1D batch
  config template has 10.
- Handoff §8 lists the engine conventions (the handoff's §2 is the regime map, not conventions).

## 12. Environment and performance

- Python 3.13; `numpy`, `pandas`, `pyyaml` preinstalled; `scipy`, `matplotlib`, `joblib`, `tqdm`
  must be pip-installed each session (cloud container is ephemeral).
- 100×100 juxtacrine run: ≈ 8 s per 5,000 steps (t 50), ≈ 1.4 min per 50,000 steps (t 500),
  single core; 4 cores available.
- `results/`, `runs/`, `*.png`, `*.npz`, `*.mp4` are gitignored. Put new analysis code under a
  tracked directory (e.g. `analysis/`), raw outputs under `results/`, and commit summary CSVs
  with the cut sizes and seeds used.

## 13. Working rules for this repo

- Do not modify `simulation_2D.py` / `simulation.py` when reproducing published panels;
  wrap them. If a bug fix is needed (RNG seeding, visualiser parity), make it opt-in and
  record which runs used it.
- Report the size-filter cut, seed count, grid size and `steps_used` with every number.
- Threshold uniform fields by variance gate first, then by level (Otsu inverts on uniform fields).
- Use the fastest-growing mode, not the determinant minimum, for wavelength predictions.
- Pool g(r) by summing counts across replicates, never by averaging curves.
