# NOTE 2026-10-17: 2D initiation sweep, run record

Run record only. Design, deviations from the deposited engine, measurements and labels are those documented in the
docstring of `analysis/sweep2d_initiation.py` (commit 24e0d45); nothing was changed. No interpretation here.

## Pre-run checks (`python analysis/sweep2d_checks.py`, 2026-10-17 02:06 UTC)

Checks 1-5 PASS. Check 4: 17 of 18 ignition tests ignited; the one exception is the documented n 10/2 at 1.0x the edge
(no ignition). Check 6 on this machine: 27.2 s per full-size run, 3,456 runs = 26.2 core-hours, projected 8.7 h on the 3
workers used.

## Execution

| item | value |
|---|---|
| command | `nohup python analysis/sweep2d_initiation.py run 3` (detached with `setsid`), cloud container, 4 cores |
| workers | 3 |
| launched | 2026-10-17 02:12 UTC |
| container restart | 08:12 UTC: the container hosting the session restarted and killed the sweep at 2,252 completed runs; the result tree survived on disk; the same command was relaunched at 08:12 UTC and resumed at the 2,253rd run (earlier log kept as `results/sweep2d/run_before_restart_0812.log`) |
| finished | 14:22 UTC |
| wall time | 12.2 h in two segments: 5.98 h for 2,252 runs, 6.16 h for 1,204 runs (the restarted container ran the first hours of the second segment at about half speed, then recovered) |
| per-run time | mean 37.9 s, median 29.5 s, max 79.2 s (the `seconds` column of `sweep2d_runs.csv`) |
| checkpoints | commits 3dc0cf4 (1,303 runs; replicate 0 complete), 42dc969 (1,808), 9d42117 (2,298), 0608a1e (2,801), cd12e25 (3,300), and the final commit (3,456) |

## Completion

| replicate | runs completed | of |
|---|---|---|
| 0 | 1,296 | 1,296 (every condition: 216 seed1 + 1,080 with replicates) |
| 1 | 1,080 | 1,080 |
| 2 | 1,080 | 1,080 |
| total | 3,456 | 3,456 |

Per start: seed1 216, seeds 648, rf1 648, rf2 648, nuc01 648, nuc001 648. Failures: `results/sweep2d/failed/` holds 0 files.
Labels at t = 50 over all runs: spots 1,492, sheet 941, dead 680, network 343 (per condition in `sweep2d_by_condition.csv`).

## Numerical health

- Runs with `i_min` below -1: 0. Minimum `i_min` over all runs is 0.0; no run has a negative inhibitor minimum.
- Non-finite values: none (no inf anywhere; `a_max`, `i_max`, `i_min`, `coverage` finite in all 3,456 rows). Maximum
  `a_max` 12.0 (= H_ON), maximum `i_max` 38.3.
- NaN appears only where the metric is undefined by construction: `t_first` 650, `t90` 687, `size_*` and
  `diam_over_lambda` 680 (the dead runs), `uniformity_t50` 2,113 and `uniformity_day4` 1,609 (fewer than 10 spots).

## Outputs

- `analysis/summaries/sweep2d_runs.csv` (3,456 rows, 28 columns) and `analysis/summaries/sweep2d_by_condition.csv`
  (1,296 rows), written by `python analysis/sweep2d_initiation.py aggregate`.
- `analysis/summaries/sweep2d_rows/` : copies of the 3,456 per-run JSON row files (14 MB on disk, about 2.5 MB of content).
- Raw output on the cloud container only (`results/` is gitignored): `results/sweep2d/` = 746 MB (761,744,645 bytes):
  `frames/` 3,456 npz files (activator at the day-1 to day-4 equivalents, t = 30 and t = 50 in float16, inhibitor at
  t = 50; about 212 KB each), `rows/`, `edges.json`, the two run logs. Not committed. Regenerate with the same command
  (resumable; finished runs are skipped).
