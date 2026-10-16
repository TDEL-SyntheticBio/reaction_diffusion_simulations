# Nucleation-rate sweep and day-4 spatial order, 15 October 2026

Added here so the two runs are reproducible against a named commit. The full write-up, with
the comparisons against experiment and the list of what it corrects, is the project doc
`claude/model_vs_invitro_calibration_Oct12.md` (amended 15 Oct).

## What is new in this commit

| file | what it is |
|---|---|
| `analysis/nucleation_sweep.py` | the sweep driver. 8 rates x 3 replicates at the published Fig 2L set, to t = 30, frames every t = 0.5 |
| `analysis/measure_nucleation_sweep.py` | measures it at day 4 and prints the comparison table |
| `analysis/spatial_order.py` | the g(r) and S(k) estimator used on both simulations and images |
| `analysis/summaries/nucleation_sweep_per_run.csv` | per-run output of the sweep |
| `analysis/summaries/invitro_well_summary_2026-10-14.csv` | 37 segmented whole wells: counts, density, diameter, size CV |
| `analysis/summaries/gr_pooled_day4.csv`, `sk_pooled_day4.csv` | pooled day-4 g(r) and S(k), in vitro and simulation, with nulls |
| `analysis/summaries/per_well_order_day4.csv` | per-well order statistics |
| `analysis/summaries/log2sd_by_sizecut.csv` | log2 area SD at the 2-cell and 6-cell cuts |

To re-run:

```
python3 analysis/nucleation_sweep.py          # ~7 min on 2 cores, writes results/nucleation_sweep/
python3 analysis/measure_nucleation_sweep.py  # prints the table, rewrites the summary CSV
```

Set `NUCSWEEP_OUT` to put the fields elsewhere. `results/` is gitignored, so the raw `.npz`
are not in the repo; they are archived separately (see below) and are reproducible from this
commit with the seeds in the driver.

## Measurement conventions, which matter more than they look

These are the conventions used throughout, and two of them are corrections to earlier work.

- **Threshold** `0.3 * act_prod_rate` on the activator, fixed, not half-max. Same as
  `analysis/fig1c_batch.py`.
- **Adjacency** the engine's own even-r hex neighbours, via `analysis/hexgeom.py`. Not
  `res_analysis/npz_feature_distribution.py`, whose default layout is the mirror image.
- **Minimum object size: 2 cells, not 6.** The experimental masks cut at 40 um equivalent
  diameter. On this lattice at 28 um spacing one cell is 29 um and two cells are 41 um, so
  two cells is the matched cut. Using 6 cells instead shifts the mean diameter +23 um and
  drops the log2 area SD from 1.66 to 0.69, which is what produced the retracted claim that
  the model's size distribution was about twice too narrow.
- **lambda_i on the hex lattice is sqrt(1.5 * D_i / gamma)**, not sqrt(D_i / gamma).
- **Time** 1 simulation unit = 5.4 h, from the JA front speed against Santorelli 2024 Fig 3D.
  Day 4 = t = 17.8. The overnight batches read t = 600, which is day 135.

## g(r) and S(k)

`spatial_order.py` uses a translation-corrected (Ohser) pair correlation so no point is
discarded to a border guard, with separate numerator and denominator so windows of different
size and density pool exactly. Separations are in units of each window's own mean spacing
a = 1/sqrt(density), which is what lets wells of different density be pooled without smearing
the first shell. Three nulls are provided: uniform points, hard-core discs at the observed
size distribution, and a germ-grain model (random placement with overlap, re-segmented) at
matched coverage. The germ-grain one is the null that matters for the exclusion zone.

Validated on CSR: translation-corrected g(r) returns 0.994 over 12 realisations of 800 points.

## Headline numbers

Sweep, published parameter set, read at day 4:

| nucleation rate | domains/mm2 | diameter | log2 area SD |
|---|---|---|---|
| 0.005 | 0.64 | 160 um | 1.56 |
| **0.01 (published)** | **2.45** | 166 um | 1.64 |
| 0.02 | 5.65 | 149 um | 1.68 |
| 0.04 | 11.6 | 115 um | 1.78 |

The eight measured day-4 JAPI wells average 2.37 domains/mm2, so the published rate is right
to 3%. But across the measured density range the model's diameter moves only 162 to 170 um
while the wells span 101 to 235 um, so nucleation rate does not account for the well-to-well
size variation.

Order, 8 wells pooled, 4075 domains: one correlation shell at 0.81 of the mean spacing
(g = 1.236 +- 0.032), no second shell (0.018 +- 0.021), no peak in S(k) at any wavenumber.
Testing each well against its own coverage-matched null, JAPI and JA are indistinguishable
(p = 0.80).
