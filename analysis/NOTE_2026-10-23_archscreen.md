# NOTE 2026-10-23: architecture screen, run record

Run record only. Design, circuits, starts, metrics and labels are those of `analysis/archscreen.py` (commit 7b27d51);
nothing was changed and no bug had to be fixed. No interpretation here.

## Commands and timing (cloud container, 4 cores, 3 workers throughout)

| step | command | started (UTC) | wall time | result |
|---|---|---|---|---|
| 1 | `python analysis/archscreen.py check` | 2026-10-22 23:48 | seconds | 0.0 differences for JAPI and for PAPI |
| 2 | `python analysis/archscreen.py count` | 23:48 | seconds | 108 configs, 2592 main jobs, 72 check200 jobs |
| 3 | `python analysis/archscreen.py thresholds 3` | 23:49 | 353 s | 108 thresholds; copied to `analysis/summaries/archscreen_thresholds.json` (commit 061cfb2) |
| 4 | `setsid nohup python analysis/archscreen.py run 3` | 23:55 | 1.48 h (ended 01:25) | 2592 of 2592 runs, N = 100 |
| 5 | `python analysis/archscreen.py check200 3` | 01:25 | 0.20 h (ended 01:38) | 72 of 72 runs, N = 200 |
| 6 | `python analysis/archscreen.py aggregate`, then `sheets` | 01:38 | about 2 min | CSVs and 72 sheets |

Restarts: none. The main process and check200 ran once each to completion. Save points (aggregate + commit of the
CSVs, no sheets): 352a270 (511 runs), 0c36da1 (1,053), d5ab958 (1,509), 3fa33d7 (2,051), 191f9f3 (2,512), then the
final commit 349cdab (2,592 main + 72 check200, with the sheets).

## Failures and numerical health

- `results/archscreen/failed/`: 0 files. `results/archscreen/check200/failed/`: 0 files.
- Non-finite values (`finite` column False): 0 of 2,592 main runs, 0 of 72 check200 runs.
- Inhibitor minimum below -1 (`b_min < -1`): 0 of 2,592 and 0 of 72. The minimum `b_min` is 0.0 in both tables; no run has
  a negative inhibitor anywhere.
- Per-run time (`seconds`): main mean 6.1 s, max 13.0 s; check200 mean 29.4 s, max 42.5 s.

## Thresholds (single-cell ignition, 108 configurations)

36 configurations ignite from one cell; 72 never ignite up to a kick of 100 and the driver uses kick = 100 for them.
Never igniting: every configuration with w = 0 (36 of 36); with w = 0.25: all nine circuits of the separate promoter at
n 10/4, and all but A intra / B membrane at n 4/4 for both promoters; with w = 1: A intra with any B, A membrane / B intra
(shared), A diffusible / B intra and B membrane (shared, both Hill pairs), and at n 10/4 everything except JAPI and PAPI.
Full table with thresholds and kicks: `analysis/summaries/archscreen_thresholds.csv` (and the json copy).

## Outputs

- `analysis/summaries/archscreen_runs.csv` (2,592 rows, 44 columns; starts noise / nuc001 / seeds, 864 each, replicate 0),
  `archscreen_check200_runs.csv` (72 rows), `archscreen_thresholds.csv`, `archscreen_thresholds.json`,
  `archscreen_edges.csv`.
- Sheets: `analysis/summaries/archscreen_sheets/` holds 72 PNGs (t = 50 and day-4 sets, one per Hill pair, promoter, w and
  start; 32 MB), force-added because PNGs in subfolders are gitignored.
- Raw output on the cloud container only (`results/` is gitignored): `results/archscreen/` = 46,542,878 bytes (58 MB on
  disk): `rows/`, `frames/`, `check200/`, `thresholds.json`, `run.log`, `check200.log`, `thresholds.log`. Not committed.
