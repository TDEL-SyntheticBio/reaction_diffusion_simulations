# CLAUDE.md — JAPI / PAPI reaction-diffusion engine

Fork of the simulation code for the JAPI patterning paper (in revision). The engine is
unchanged from upstream commit `07a328b`. Everything below was verified against the source
on 2026-10-07: read by hand, checked numerically against the repo's own functions, and then
re-verified by twelve independent adversarial passes (two per claim bundle). Each item carries
a `file:line` so it can be re-verified. Do not take the handoff documents' descriptions of the
engine on trust; take these.

Supporting material lives in `handoff/`:

- `JAPI_irregular_handoff_2026-10-07.md` — regime map, g(r) methodology, results, open items.
  Its engine-convention list is **section 8** ("Reproducing"); section 2 is the regime map.
- `JAPI_session_findings_Oct2026.md` — earlier session (asynchronous nucleation result)
- `SupplementaryNote1..4_*.md` — the paper's supplementary notes (LSA, velocity, nondimensionalisation, arrest)
- `supp_table1_panel_parameters.csv`, `supp_table2_coupled_panel_parameters.csv` — per-panel
  reproduction keys. Column names are the code parameter names. A blank `nucleation_rate`
  means none and must be passed as `0.0` explicitly (see §9).

## 1. Model as implemented (2D, `2D_simulations/simulation_2D.py`)

Forward Euler, one step per call of `_vectorized_step` (`simulation_2D.py:145-195`); the
two-pair engine repeats the same formulas in `_pair_hill` (`:309-345`) and
`_vectorized_step_two_pairs_cross` (`:348-411`). Any change to the Hill form must be made in
**both** places (`:181` and `:344`).

```
hill  = (A + basal) / (A + I + 1 + basal),   A = (act_signal/ka)^n_a,  I = (inh/ki)^n_i
a    += dt * (act_prod_rate * hill - act_decay_rate * a) [+ act_diffusion * dt/dx^2 * lap(a)  if paracrine]
i    += dt * (inh_prod_rate * hill - inh_decay_rate * i) + inh_diffusion * dt/dx^2 * lap(i)
```

- **juxtacrine**: `act_signal` = mean of `a` over the *valid* neighbours, self excluded
  (`:173`). Interior weight 1/6 per neighbour; edge cells average over fewer neighbours
  (`nbr_count`), never over self. `act_diffusion` is ignored (`:192-193`) even though several
  configs set it to 1.0.
- **paracrine**: `act_signal = a` (self-sensing) and the activator diffuses with the same
  raw-sum Laplacian (`:175,189-191`), so it also carries the 1.5 lattice factor of §3.
- `activator_type` is matched by exact string. In 2D the signal branch tests `== "juxtacrine"`
  and the diffusion branch tests `== "paracrine"`; any other string ("soluble", "Juxtacrine")
  silently gives a self-sensing, non-diffusing hybrid. In 1D anything not `"paracrine"` is
  juxtacrine. Only the two exact strings are safe.
- `basal_prod` must be present in the parameter dict: the engine indexes it (`:164`) while
  the solver defaults it (`finding_steady_states.py:46`), so a config without the key gets an
  a_ss and then crashes on the first step.
- `dt = 0.01`, `dx = 1.0` everywhere. **Time unit: t = steps × 0.01**: 5,000 steps = t 50,
  50,000 steps = t 500. `dt` also sets the nucleation probability (`rate*dt`) and the noise
  scale (`sqrt(dt)`), so lowering dt keeps both consistent per unit time.
- Negative values are never clipped except the activator after additive noise (`:287`).
  Hill terms are zeroed for non-positive signals (`:179-180`), so negatives do not make NaN.
- Order within one loop iteration: Euler update → nucleation kick (mask evaluated on the
  post-Euler field) → noise + clip → every `save_every` steps: save + stopping check.
- Cross-inhibition (two-pair engine): `p1.cross_inhibition_rate * hill2` is added to pair 1's
  **inhibitor** production and vice versa (`:375-386`). The per-pair steady states are
  computed without it (`:457,465`); with nonzero cross rates the 2×2 closed forms of §10 do
  not apply.

## 2. Hex lattice: even-r offset, zero-flux

`_build_hex_neighbor_arrays` (`simulation_2D.py:108-142`):

- **even-r**: *even* rows are shifted half a cell to the **right**. Even-row cell (r,c) has
  upper neighbours (r−1,c),(r−1,c+1); odd-row cell has (r−1,c−1),(r−1,c). Verified: 0/168 cells
  disagree with an explicit even-r coordinate construction; 168/168 disagree with odd-r.
- Physical centre coordinates for any new plotting or distance measurement:
  `x = c + 0.5*(r % 2 == 0)`, `y = r*sqrt(3)/2` (unit centre-to-centre spacing). "Lattice
  units" means centre-to-centre spacing; the row pitch is √3/2, not 1.
- **Zero-flux boundaries**: out-of-bounds neighbours are masked; `nbr_count` = number of
  valid neighbours (6 interior; 5 on edges; 2, 3 or 4 at corners depending on row parity).
  It divides the juxtacrine signal and multiplies the self term of the Laplacian, so at every
  cell `lap = nbr_count * (mean_of_valid_neighbours − u)`. The neighbour relation is symmetric.

### Hex parity is wrong in the visualiser and in the component finder

- `visualize_2D._hex_polygons` (`visualize_2D.py:46-50`) shifts **odd** rows right (odd-r)
  although its docstring (`:37`) says even-r. The drawn adjacency equals the engine adjacency
  of the column-reversed field, so every published PNG/MP4 from `plot_one_frame`,
  `animate_histories`, `animate_four_histories`, `plot_four_frames` misdraws adjacency at the
  single-cell level. Fixing it needs `:48` (`r % 2 == 0`), `_hex_bounds :66` (constant +0.5)
  and the docstring. `_add_overlay_hex_field_rgb` (`:350-355`) uses a third, flat-top
  geometry (orientation 30°, dx = 1.5R, dy = √3R) that is not a consistent tiling and is
  stretched 4/3 relative to the four flanking panels.
- `res_analysis/npz_feature_distribution.py get_neighbors` (`:69-109`): its `even-r` branch
  (the **default** `--layout`) returns the engine's *odd-r* mirror: 168/168 cells get a
  different neighbour set from the engine; `--layout odd-r` matches the engine on 168/168.
  Every component count, size or diameter produced by the deposit with the default layout used
  the wrong adjacency. Pass `--layout odd-r`, or build components from
  `simulation_2D._build_hex_neighbor_arrays` directly (preferred). Its
  `approx_diameter_hexagons` returns a diameter in hex **side** lengths (a single hex → 1.82),
  not cell diameters (flat-to-flat = √3 sides).
- `radial_autocor_npz.py`, `crosscorr_npz.py`, `fft_analysis_2D.py`, `npz_overlap_analysis.py`
  treat the (Ny, Nx) index array as a square lattice (no √3/2 row pitch, no row stagger):
  distances along rows are overstated by 15% and the FFT is anisotropically distorted. Use
  the coordinates above for any distance-based statistic.

## 3. Laplacian: raw sum, prefactor 6, effective D is 1.5× nominal

`lap_i = sum(valid neighbours) − nbr_count * i` (`simulation_2D.py:186`), scaled by
`inh_diffusion * dt / dx**2`, with no hex geometric factor. This is the raw neighbour sum,
**not** the continuum-matched hex Laplacian (which would be (2/3)× this, self-coefficient 4).
Verified numerically:

| quantity | value |
|---|---|
| raw-sum Laplacian of x²+y² (continuum value 4) | 6.000 → prefactor **1.5** |
| Fourier symbol Λ_hex(q) = 6 − Σ_j cos(q·n_j) | range [0, 9]; Λ_hex ≈ 1.5\|q\|² at small q |
| juxtacrine kernel symbol K̂_hex(q) = (1/6) Σ_j cos(q·n_j) | range [−1/2, 1] |
| identity (bulk) | **Λ_hex = 6 (1 − K̂_hex)** exactly (same neighbour set); 1D: Λ = 2 (1 − cos q) |
| inhibitor decay length on the lattice, D_i=10, γ=0.5 | fitted **5.5–5.7 cells** ≈ √(1.5·D_i/γ) = 5.48; nominal √(D_i/γ) = 4.47 |
| same, D_i=20 | fitted 7.9–8.3 cells ≈ 7.75; nominal 6.32 |
| same, D_i=5 | fitted 3.9 cells ≈ 3.87; nominal 3.16 |

So: the dispersion relation written with the lattice symbol, Δ/γ = (1 + λ²Λ(q))(1 − αK̂(q)) + G
with λ² = D_i/γ, is internally consistent (the 1.5 lives inside Λ_hex). But **any statement of
the inhibitor range in cells must use √(1.5·D_i/γ)**. The handoff's "λ_i = 4.5 cells (D_i=10),
6.3 cells (D_i=20)" understates the engine's range by 22%, and its "g(r) hole sits at 2λ_i"
statements were made with the nominal value. Fig 1F/1G (1D, D_i=10) and Fig 1C (2D, D_i=10)
are therefore not at the same continuum inhibitor range.

The 1D engine (`1D_simuations/simulation.py:101-112`) has prefactor 1 (two neighbours; Λ ∈ [0, 4]),
so √(D_i/γ) is correct in 1D (exact lattice decay length 1/acosh(1 + γ/2D) = 4.48 at D=10, γ=0.5).

## 4. Forward-Euler stability and positivity

- Most oscillatory inhibitor mode multiplier = 1 − dt(γ + 9 D_i). At dt = 0.01, γ = 0.5 the
  bound is **D_i < 22.2** (the finite zero-flux grid's top eigenvalue is 8.998 at 100×100,
  8.956 at 20×20, so marginally higher there). Verified empirically: D_i = 22 stays bounded,
  **D_i = 23 blows up**. The committed D_i = 20 (Fig 2L/3/S8D) has ~10% margin. Reaction
  stiffness (γ − b_i·∂H/∂i ≥ γ at an activated state) shrinks the margin further at high b_i
  or n_i. Any sweep of `inh_diffusion` above ~22 at dt = 0.01 is numerically meaningless:
  reduce dt. Same bound for a paracrine `act_diffusion` (D_a < 22.1). 1D bound ≈ 50.
- Blow-up signature: the inhibitor alternates sign with magnitude growing to 1e171 while
  `np.isfinite` is still True and the activator stays bounded (Hill terms are zeroed for
  non-positive signals). Test for negative inhibitor or growth of max|i|, not finiteness.
- Positivity: the monotone-stencil condition dt(γ + 6 D_i) ≤ 1 fails at D_i = 20 (self
  coefficient −0.205), so the inhibitor can go transiently negative near sharp features and
  is never clipped. At D_i = 10 the stencil is monotone.

## 5. Parameters and nondimensionalisation

In **every** committed parameter/config file (`2D_simulations/parameters_2D.py`,
`parameters_coupled_2D.py` p1 and p2, `2D_batch/config.yaml`, `2D_batch/config_coupled.yaml`
p1 and p2, `1D_simuations/parameters.py`, `1D_batch/config.yaml`):

```
act_decay_rate = act_half_sat = inh_half_sat = 1.0,   basal_prod = 0.0
```

so, in the notation of Supplementary Note 3: `r_a = act_prod_rate`, `r_i = inh_prod_rate`,
`γ = inh_decay_rate`, `n_a = act_hill_coeff`, `n_i = inh_hill_coeff`, `D_i = inh_diffusion`
(lattice units, see §3). The only exception anywhere is `1D_simuations/tests.py`, which sets
`act_decay_rate = 0.0` locally and is stale and unrunnable (wrong arity, nonexistent init mode,
nonstandard activator_type strings): do not use it as a reference. `2D_batch/grid.py:78` casts
swept values to float (Hill coefficients become 10.0; harmless).

**Hill function differs between files when basal ≠ 0.** `simulation_2D.py:181,344` and
`1D simulation.py:30` compute `(A+basal)/(A+I+1+basal)`; both copies of
`finding_steady_states.hill_with_grads` compute `basal + A/(basal+1+A+I)` (can exceed 1).
Identical at basal = 0 (max diff 2e-16 over 2000 random points); the pointwise gap is
basal·(basal+A+I)/(basal+1+A+I), approaching basal as A → ∞. Never set basal ≠ 0 without
reconciling the two; none of the §10 closed forms hold with basal ≠ 0.

## 6. Steady-state solver and its fallback — this sets the nucleation kick

`run_coupled_hex` (`simulation_2D.py:228-234`; two-pair `:456-470`; 1D `simulation.py:166-173`)
calls `fast_stable_steady_state`; if it raises or returns a non-positive or non-finite `a_ss`
**or `i_ss`**, the engine sets `a_ss = i_ss = spike_value`. The fallback also feeds the
`spike_steady_state` and `random_tight` initial conditions and the a_ss/i_ss metadata written
to npz/CSV. Solver output for the published panels (JAPI, κ = 1, values rounded to 3 dp by the
solver):

| panel | b_a | b_i | γ | n_a | n_i | solver a_ss, i_ss | a_ss the engine uses | kick = 2·a_ss |
|---|---|---|---|---|---|---|---|---|
| 1C all on | 5 | 1 | 0.5 | 3 | 3 | 4.658, 1.863 | 4.658 | 9.3 |
| 1C stripes | 5 | 3 | 0.5 | 3 | 3 | 1.706, 2.048 | 1.706 | 3.4 |
| 1C reg / irreg / off spots | 5 | 5,12,14 | 0.5 | 3 | 3 | none | 1.0 (spike_value) | 2.0 |
| 2L JA | 5 | 0 | 0.5 | 10 | 4 | 5.0, **0.0** → fails i_ss>0 | 2.0 (spike_value) | **4.0** |
| 2L JAPI, 3B spots, 3H, 3L, S8D | 5 | 15–50 | 0.5 | 10 | 4 | none (exact scan confirms) | 2.0 (spike_value) | **4.0** |
| 3B labyrinth | 6 | 15 | 0.5 | 10 | 4 | 5.911, 29.554 (exact 5.914, 29.568) | 5.911 | **11.8** |
| S5F irregular bistable | 5 | 5 | 0.4 | 6 | 3 | 4.043, 10.11 | 4.043 | 8.1 |
| S5F irregular monostable | 5 | 10 | 0.5 | 6 | 3 | none | 1.0 | 2.0 |
| coupled p1 (4B/4L/5B) | 5 | 8 | 0.5 | 4 | 4 | none | 2.0 | 4.0 |

So every experiment-matched panel at b_a = 5 nucleates with a kick of **4.0** (= 2 × spike_value),
not 2 as the October session notes say and not a_ss as Methods say; the Fig 3B labyrinth panel
kicks with 11.8. The JA control reaches 4.0 only via the fallback. At b_a = 5, n_a = 10, n_i = 4
an activated state exists only below b_i ≈ 14.2, so the b_i = 15 panels sit just past the
saddle-node; b_i = 10 would have an activated state and a different kick.

**The solver is not a regime classifier.** It returns (0, 0, 0) in three distinct situations:
no activated state; an activated state that is **trace-unstable** (uniform-ON regime; e.g.
b_a = 20, b_i = 27–28, n = 3, γ = 0.5), because `_is_reaction_stable` requires both eigenvalues
negative; and a **saturated** activated state whose H* lies in (1−1e-9, 1), outside its bracket
(230 of 1,579 random sets with roots, including sets inside the published grid). Its "no root"
detection is also fragile at large b_i/γ (margin 3e-5 at I = 100). For any closed-form or
classifier work refine roots exactly (the solver uses tol 5e-4 in H and rounds to 3 dp; the
error in a scales as A·5e-4, up to 1.5e-2 in α). Always log a_ss/i_ss for every run (the batch
CSVs already do).

Two copies of the solver exist with the same module name (`1D_simuations/` and
`2D_simulations/finding_steady_states.py`); only the `_mini_brent` safeguard differs, giving
3rd-decimal differences in a_ss (0 qualitative disagreements over 1,584 sets). Importing both
engines in one interpreter binds both to whichever copy was imported first: run 1D and 2D
checks in separate processes.

## 7. Nucleation and randomness

- Each step every cell fires with probability `nucleation_rate * dt` (`simulation_2D.py:263,281`).
  A firing cell gets `a += 2 * a_ss` **only if its activator is below a_ss** (`:282-283`), so a
  cell just under a_ss lands just under 3·a_ss. Nucleation and noise act on the activator only.
  On 100×100 at rate 0.01 that is ~1 kick per step while most cells are off (0.94 measured);
  ~0.04 per step on the 20×20 Fig 3L grid; falling to ~0.56 in the labyrinth run where 44% of
  cells end above a_ss.
- Nucleation, `noise_amplitude` noise, and the 2D init modes `random_uniform_over0`,
  `random_tight`, `near_zero` all draw from the **unseeded** numpy global RNG
  (`simulation_2D.py:60-99,281,286` and `:535-549`). The only seeded draw is the spike-index
  choice for `activator_random_spikes` (`init_seed = 2` hard-coded, `:236-244`, `:472-480`; both
  pairs of the two-pair engine receive the same indices). `batch_runner_2D.py:42` seeds Python's
  `random`, which nothing uses, and joblib workers each carry their own unseeded state: **no 2D
  run in the deposit is reproducible**, including the Fig 1C initial conditions. Recipe that
  works in a single process: call `np.random.seed(s)` immediately before `run_coupled_hex`
  (the per-step draw count is fixed); better, give new code an explicit `numpy.random.Generator`.
- 1D (`1D_simuations/simulation.py:33-41,179`) seeds `default_rng(2)` with the seed
  **hard-coded** and not exposed, has no nucleation and no noise: every 1D random run uses the
  one same noise realisation, and the S5F noise-level series are that field rescaled by
  spike_value, not independent samples. Replicate seeds in 1D need a code change.
- Init modes (2D): `all_off` both fields 0 (with basal = 0 and no noise, nucleation then
  supplies everything); `random_uniform_over0` a ~ U(0, spike_value), i = 0 (Fig 1C spike_value 1
  → U(0,1)); `spike_steady_state` centre cell at (a_ss, i_ss); `activator_random_spikes`
  n_points cells at `set_peak_height`; `both_on`; `random_tight` (±5% around a_ss, i_ss);
  `near_zero`; `inhibitor_spike`. 1D has `activator_spike` (Supp Table 1's 1D rows), `random_both`,
  `activator_on`, `inhibitor_on`, `two_activator_spikes` instead; **2D has no `activator_spike`**.

## 8. Stopping criterion and frame bookkeeping

- Every `save_every` steps (loop indices 0, 200, 400, …, each **after** the update) the engine
  appends a frame and sets `diff` = L1 change of a and i since the previous frame; it stops when
  `step > min_steps` and `diff/(2·Ny·Nx) < stopping_threshold` (`simulation_2D.py:289-298`;
  two-pair normalises by 4·Ny·Nx, 1D by 2·N). Committed values: threshold 1e-4, save_every 200,
  min_steps 500–1000.
- **The "final" frame is not the final state.** `history[0]` is the initial condition;
  `history[k≥1]` is the state after (k−1)·save_every + 1 Euler steps; on a full-length run the
  last frame is at loop step floor((steps−1)/save_every)·save_every, i.e. **step 4800 (t 48.01)
  for a 5,000-step run** and 49,800 (t 498.01) for 50,000, while `main_2D.py:174` labels it
  `final_step_004999` and both batch runners write `hist[-1]` as the final field. Frame labels
  in `main_2D.py --step-size` (`:188`) and `animate_histories` (`visualize_2D.py:179`) overstate
  by save_every − 1. With an early stop the break follows the append, so the last frame is
  current. The returned `step` is steps−1 on a full run or the loop index at the break; elapsed
  time is (step+1)·dt. Set `save_every` to divide (steps−1) or save the final state yourself.
- Under nucleation at 0.01 the criterion never fires (window diff ≥ 0.05 per cell, median 0.13
  on 100×100, 0.10 on 20×20), so those runs go the full `max_steps`. Without nucleation (Fig 1C
  protocol, 50,000 max): all-on (b_i = 1) and all-off (b_i = 14) stopped at steps 2,000–2,200
  while stripes, regular and irregular spots (b_i = 3, 5, 12) ran the full 50,000 with final
  window diffs 1.7e-3, 1.2e-3, 4.5e-4. The October notes' "~1e-3" threshold is not what is
  committed. Always record `steps_used`.

## 9. Entry-point traps: defaults the batch runners substitute silently

- `2D_batch/batch_runner_2D.py:92` passes `nucleation_rate = params.get("nucleation_rate", 0.02)`
  and the committed `2D_batch/config.yaml` has **no `nucleation_rate` key**, so the Fig-1C-like
  template nucleates at 0.02 although Supp Table 1 lists no nucleation for 1C. Engine defaults
  differ again: `run_coupled_hex` 0.0 (`:212`), `run_two_pair_hex` 0.01 (`:430`), `parameters_2D.py` 0.01.
  **Always pass `nucleation_rate` explicitly, including 0.0.**
- `batch_coupled_2D.py:64-80` never forwards `nucleation_rate`, `noise_amplitude` or
  `set_peak_height`: `config_coupled.yaml:17` is dead (the engine default 0.01 coincides), a YAML
  sweep over it does nothing, and co-initiation spikes use 25.0 instead of the 20 in
  `parameters_coupled_2D.py`. `batch_runner_2D.py` never forwards `n_points`/`set_peak_height`.
- Both batch runners default `init_mode` to `"activator_spike"`, which does not exist in 2D
  (ValueError), `min_steps` to 10,000, `spike_value` to 5.0, `save_every` to 100 when a YAML omits
  them. The committed YAMLs set all of these; a new YAML that drops one changes behaviour silently.
- npz keys differ by writer: `batch_runner_2D.py` writes `A_final`/`R_final`; `main_2D.py`
  writes `a`/`i` (+ `step`, `a_ss`, `i_ss`); the coupled scripts write `a1,i1,a2,i2`.
  `npz_feature_distribution.py` defaults to `A_final`.
- Reproducing Fig 1C via `main_2D.py` requires editing `parameters_2D.py` (init_mode → `random_uniform_over0`,
  spike_value → 1.0, nucleation_rate → 0.0, steps → 50,000).

## 10. Closed forms that are safe to use (verified against the code's Hill derivatives)

Preconditions: basal = 0, μ_a = act_decay_rate = 1, κ = 1 (true for every committed set; k_a = k_i = 1
is not needed). At any activated homogeneous state with P = A/(1+A+I), Q = I/(1+A+I),
A = a0^n_a, I = i0^n_i (Hill terms, not the solver's production ratios):

- α ≡ β_a f_a0 = **n_a (1 − P)**, ν ≡ −β_i f_i0 = **γ n_i Q**, G ≡ ν/γ = n_i Q
  (match `hill_with_grads` to 6 decimals on four parameter sets; general μ_a: α = μ_a n_a (1−P)).
- Δ(q)/γ = (1 + λ²Λ(q))(1 − αK̂(q)) + G, λ² = D_i/γ, Λ and K̂ the *lattice* symbols of §3
  (bulk/periodic result; edges break translation invariance).
- q = 0: determinant condition G > α − 1 (saddle-node at equality; automatic on the upper
  branch); trace condition τ(0) = α − 1 − γ(1 + G) < 0.
- Finite-q: with Λ = p(1 − K̂), p = 2 (1D) or 6 (hex), and L = pλ², Δ/γ is quadratic in
  c = K̂ with minimum at c* = 1/2 + 1/(2α) + 1/(2L). The minimum is interior (c* < 1) iff
  **λ² > α/(p(α−1))** — this is *necessary only* (Supp Note 1's D_geom). The *sufficient*
  finite-q instability condition is **G < W = [α(1+L) + L]² / (4αL) − (1+L)**; check: W = 10.00 = G
  at Supp Note 1's worked-example threshold D_T ≈ 2.46. 1D and hex windows differ (in
  physical length λ_phys² = 1.5 D_i/γ the hex necessary condition reads λ_phys² > α/(4(α−1))),
  so do not reuse 1D classifier output for 2D claims. c* is the determinant-minimum mode, not
  the fastest-growing one; use the latter for wavelength prediction.
- Cross-check: these closed forms reproduce the handoff's matched-pair numbers exactly
  (b_a = 20, n = 3, γ = 0.5: τ(0) = −0.006 at b_i = 26, **+0.023 at b_i = 27**, +0.056 at 28), so
  the b_i = 28 member of the matched pair is trace-unstable, not linearly stable, on either
  lattice; the pair straddles the trace boundary. Fig 1C series (b_a = 5): τ(0) turns positive
  between b_i = 3.6 and 3.7 while the finite-q determinant minimum is still deeply negative, and
  the activated state vanishes between 3.7 and 3.8.
- Paracrine: K̂ ≡ 1 in the (1,1) and (2,1) entries and Supp Note 1 Eq. 43 applies with Λ_hex.

## 11. Analysis code in the deposit

- No regime classifier (Turing / irregular / all-on / all-off) exists anywhere in the repo.
  `1D_simuations/1D_batch/analyze_patterns.py` is 1D-only and emits `pattern_flag` ∈
  {True, False, NaN} (spectral peak/background score > 8 and amplitude check; NaN for uniform
  fields); it runs at import, needs scipy/pandas/matplotlib, and raises a length-mismatch
  ValueError if any `activator_final` fails to parse. The 2D `res_analysis/` scripts compute
  FFTs, cross-correlations, radial autocorrelation, hex connected components and overlap
  metrics, none of which assign a regime label, and all with the geometry caveats of §2.

## 12. Known discrepancies between Methods, Supp Table 1 and code

- Nucleation increment: Methods say a_ss or spike_value; code does `2*a_ss` with the §6 fallback.
- Fig 1C amplitude: Methods say U(0,2); Supp Table 1 and `2D_batch/config.yaml` say spike_value 1 → U(0,1).
- Supp Table 1's 1D `activator_spike` rows use 2.0 (Fig 1F); the 1D batch template has 10.
  The 1D batch template's `inh_decay_rate = 1.0` matches no panel.

## 13. Environment and performance

- Python 3.13; `numpy`, `pandas`, `pyyaml` preinstalled; `scipy`, `matplotlib`, `joblib`, `tqdm`
  must be pip-installed each session (cloud container is ephemeral).
- 100×100 juxtacrine run: ≈ 8 s per 5,000 steps (t 50), ≈ 1.4 min per 50,000 steps (t 500),
  single core; 4 cores available.
- `results/`, `runs/`, `*.png`, `*.npz`, `*.mp4` are gitignored. Put new analysis code under a
  tracked directory (e.g. `analysis/`), raw outputs under `results/`, and commit summary CSVs
  with the cut sizes and seeds used.

## 14. Working rules for this repo

- Do not modify `simulation_2D.py` / `simulation.py` when reproducing published panels;
  wrap them. If a bug fix is needed (RNG seeding, final-frame save, visualiser parity), make it
  opt-in and record which runs used it.
- Report the size-filter cut, seed count, grid size and `steps_used` with every number.
- Threshold uniform fields by variance gate first, then by level (Otsu inverts on uniform fields).
- Use the fastest-growing mode, not the determinant minimum, for wavelength predictions.
- Pool g(r) by summing counts across replicates, never by averaging curves.
- Quote inhibitor ranges in cells as √(1.5·D_i/γ) in 2D and √(D_i/γ) in 1D, and say which.
