# NOTE 2026-10-24: architecture screen extension, run record

Run record only. The extension (seven more Hill pairs, eight key circuits) is the "extension" section of
`analysis/archscreen.py` at commit 5c18b1e; nothing was changed and no bug had to be fixed. No interpretation here.

## Commands and timing (cloud container, 4 cores, 3 workers throughout)

| step | command | started (UTC) | wall time | result |
|---|---|---|---|---|
| 1 | `python analysis/archscreen.py check` | 2026-10-24 04:49 | seconds | 0.0 differences for JAPI and for PAPI |
| 2 | `python analysis/archscreen.py ext_count` | 04:49 | seconds | 56 configs, 1026 jobs |
| 3 | `python analysis/archscreen.py ext_thresholds 3` | 04:52 | 192 s | 56 thresholds; copied to `analysis/summaries/archscreen_ext_thresholds.json` (commit 41ee940) |
| 4 | `setsid nohup python analysis/archscreen.py ext_run 3` | 04:57 (first job 05:00) | 0.56 h (ended 05:34) | 1026 of 1026 runs, N = 100 |
| 5 | `python analysis/archscreen.py ext_aggregate`, then `ext_sheets` | 05:34 | about 1 min | CSVs and 42 sheets |

Restarts: none; the run completed in one go. Save point: 7e799d8 (523 runs aggregated, CSVs only). Final commit ae6aef5
(1,026 runs, with the sheets). Total time from the first check to the final push: about 47 minutes.

## Failures and numerical health

- `results/archscreen/ext/failed/`: 0 files.
- Non-finite values (`finite` False): 0 of 1,026 runs.
- Inhibitor minimum below -1 (`b_min < -1`): 0 of 1,026; the minimum `b_min` is 0.0, no run has a negative inhibitor.
- Largest values: `a_max` 93.4, `b_max` 161.3. Per-run time (`seconds`): mean 5.9 s, max 10.1 s.

## Thresholds (single-cell ignition, 56 configurations)

22 configurations ignite from one cell; 34 never ignite up to a kick of 100 and the driver uses kick = 100 for them.
Never igniting: LI (w 0) and A intra / B membrane (w 0) at every Hill pair except n 4/10 and n 10/10 (where both ignite
at 1); JAJI (w 0.25) at all seven pairs; JAJI (w 1) at n 2/10, 4/2 and 10/2; PAPI at n 10/2; PAJI at every pair except
n 10/10 (threshold 80.6); single-morphogen RD at all seven pairs. JAPI ignites at all seven pairs (threshold 6.26 at
n_a 2, 8.98 at 4/2, 8.66 at 4/10, 10.0 at 10/2, 7.50 at 10/10). Full table: `archscreen_ext_thresholds.csv` and the json copy.

## Outputs

- `analysis/summaries/archscreen_ext_runs.csv` (1,026 rows, 44 columns; starts noise / nuc001 / seeds, 342 each, replicate
  0; 144 runs per Hill pair, 162 at n 10/2 where JAJI (w 0.25) has an extra edge position), `archscreen_ext_thresholds.csv`,
  `archscreen_ext_thresholds.json`, `archscreen_ext_edges.csv`.
- Sheets: `analysis/summaries/archscreen_ext_sheets/`, 42 PNGs (t = 50 and day-4 sets, one per Hill pair and start; 20 MB),
  force-added because PNGs in subfolders are gitignored.
- Raw output on the cloud container only (`results/` is gitignored): `results/archscreen/ext/` = 25,299,302 bytes (30 MB
  on disk): `rows/`, `frames/`, `thresholds.json`, `run.log`, `thresholds.log`. Not committed.
