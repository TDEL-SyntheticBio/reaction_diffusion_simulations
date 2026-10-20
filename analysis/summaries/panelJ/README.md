# Panel J read-off and the 10/4 dip check (1D)

Two working tables from the Fig 1 panel J reconstruction (16 October 2026). 1D only.

## panelJ_readoff.csv

Values read off the published panel J image, one row per (n_a, n_i) circle.

| column | meaning |
|---|---|
| n_a, n_i | Hill pair |
| area | circle area in image pixels |
| n_sets | circle size converted to a count of parameter sets, via the radius calibration count = 34.94 r - 12.5 |
| frac | colour read off the colour bar |
| err | read-off uncertainty estimate for that circle |

Findings from comparing it with analysis/summaries/spike_grid_full.csv: circle size tracks the count of the deposit's
pattern_flag from a spike (Spearman 0.984); colour tracks the share of patterned sets with several domains versus one
(Spearman 0.976, against 0.905 for the Turing count). n_a = 1 is absent from panel J. At n 10/4 the published colour is
about 0.05; our contiguous-domain labelling gives 0.

## dip_check_10_4.csv

Every set at n 10/4 that ends as one domain, from a single-spike run at N = 100 with the deposit's settings.

| column | meaning |
|---|---|
| t | set index within the 10/4 run |
| b_a, b_i, gamma | parameters |
| size | domain size in cells (contiguous, above 0.3 b_a) |
| dip_rel | depth of the dip in the middle of the domain's top, relative to its maximum |
| peaks_f64, peaks_f32 | peaks found by scipy find_peaks on the final profile in float64 and in float32 |

Summary: 113 of 217 single domains have a dip, 69 deeper than 0.1%, none deeper than 1%. find_peaks reports 2 or more
peaks in 1 set at float64 (b_a 3, b_i 10, gamma 2, an 11-cell domain with a 0.25% dip) and in 42 sets at float32.
Hypothesis, not yet checked against Ben's own files: the 0.05 in panel J at 10/4 comes from peak counting of such
domains.
