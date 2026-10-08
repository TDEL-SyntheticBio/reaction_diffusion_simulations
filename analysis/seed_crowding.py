"""
Does the fate of one isolated seed predict what a whole field does, and how long does each take to settle?

Parameter sets: Fig 1C line (b_a 5, n 3/3, D_i 10, gamma 0.5) at b_i = 5, 9, 12.
Six initial conditions per set, kick level 5.0 wherever the mode allows it:
  seed1    spike_steady_state, spike_value 5.0 (deterministic; one run)
  seed9 / seed36 / seed144   activator_random_spikes, n_points 9/36/144, set_peak_height 5.0
  random   random_uniform_over0, spike_value 1.0, no nucleation      (Fig 1C protocol)
  nucl     all_off, spike_value 2.0, nucleation_rate 0.01            (Fig 2L protocol, kick 4.0)
Three replicate RNG seeds for everything but seed1. 141x141, dt 0.01, to t = 600, convergence check
disabled (stopping_threshold 0, min_steps > steps). The engine saves every 200 steps so its own
convergence diagnostic (L1 change per cell-field over a 200-step window, CLAUDE.md s8) is recorded for
every window; measurements are taken on the frames at t = 10, 20, ..., 600.

Engine caveat (CLAUDE.md s7): activator_random_spikes draws its positions from the hard-coded
init_seed = 2, so replicates would be identical through the public API. The worker overrides
numpy.random.default_rng during the call so that replicate r draws positions from seed 1000*r + 7;
the seed used is recorded in the CSV. Everything else is seeded with np.random.seed(r).

Usage: python analysis/seed_crowding.py run [n_workers]     (~1 h on 4 cores; writes results/seed_crowding/*.npz
                                                             and analysis/summaries/seed_crowding_frames.csv)
       python analysis/seed_crowding.py analyze               (settling table, criterion table, figures)
"""
from __future__ import annotations

import contextlib, csv, io, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "2D_simulations"))
from analysis.hexgeom import label_components  # noqa: E402

N, DT, SAVE_EVERY, T_END = 141, 0.01, 200, 600.0
STEPS = int(round(T_END / DT)) + 1                      # 60001: last frame is the final state
THR = 1.5                                               # 0.3 * r_a
BIS = (5, 9, 12)
CONDITIONS = {
    "seed1":   dict(init_mode="spike_steady_state", spike_value=5.0, n_points=0, set_peak_height=5.0, nucleation_rate=0.0),
    "seed9":   dict(init_mode="activator_random_spikes", spike_value=5.0, n_points=9, set_peak_height=5.0, nucleation_rate=0.0),
    "seed36":  dict(init_mode="activator_random_spikes", spike_value=5.0, n_points=36, set_peak_height=5.0, nucleation_rate=0.0),
    "seed144": dict(init_mode="activator_random_spikes", spike_value=5.0, n_points=144, set_peak_height=5.0, nucleation_rate=0.0),
    "random":  dict(init_mode="random_uniform_over0", spike_value=1.0, n_points=0, set_peak_height=5.0, nucleation_rate=0.0),
    "nucl":    dict(init_mode="all_off", spike_value=2.0, n_points=0, set_peak_height=5.0, nucleation_rate=0.01),
}
REPS = {"seed1": [0], "seed9": [1, 2, 3], "seed36": [1, 2, 3], "seed144": [1, 2, 3], "random": [1, 2, 3], "nucl": [1, 2, 3]}
OUT = ROOT / "results" / "seed_crowding"
SUMM = ROOT / "analysis" / "summaries"
FRAMES_CSV = SUMM / "seed_crowding_frames.csv"
SNAP_T = (60, 150, 300, 600)


def params(bi):
    return dict(act_half_sat=1.0, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=1.0, act_prod_rate=5.0,
                inh_prod_rate=float(bi), inh_decay_rate=0.5, act_hill_coeff=3, inh_hill_coeff=3, inh_diffusion=10.0)


def frame_time(k):
    """Frame k (k >= 1) is the state after (k-1)*SAVE_EVERY + 1 Euler steps (CLAUDE.md s8); frame 0 is t = 0."""
    return 0.0 if k == 0 else ((k - 1) * SAVE_EVERY + 1) * DT


def measure(a):
    mask = a > THR
    labels, n = label_components(mask) if mask.any() else (None, 0)
    largest = int(np.bincount(labels[labels >= 0]).max()) if n else 0
    return dict(cells_on=int(mask.sum()), blobs=n, largest_blob=largest, coverage=float(mask.mean()), total_activator=float(a.sum()))


def run_one(job):
    bi, cond, rep = job
    from simulation_2D import run_coupled_hex
    c = CONDITIONS[cond]
    spike_seed = ""
    np.random.seed(rep)
    orig_rng = np.random.default_rng
    if c["init_mode"] == "activator_random_spikes":
        spike_seed = 1000 * rep + 7
        np.random.default_rng = lambda seed=None, _s=spike_seed: orig_rng(_s)   # override the engine's hard-coded init_seed
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            A, I, step, a_ss, i_ss = run_coupled_hex(N, N, STEPS, DT, 1.0, params(bi), 0.0, 10**9, init_mode=c["init_mode"],
                                                     activator_type="juxtacrine", n_points=c["n_points"], spike_value=c["spike_value"],
                                                     save_every=SAVE_EVERY, nucleation_rate=c["nucleation_rate"], noise_amplitude=0.0,
                                                     set_peak_height=c["set_peak_height"])
    finally:
        np.random.default_rng = orig_rng
    assert len(A) == STEPS // SAVE_EVERY + 2 and abs(frame_time(len(A) - 1) - (T_END + DT)) < 1e-9
    # the engine's own convergence diagnostic for every 200-step window (what stopping_threshold is compared with)
    crit = np.array([(np.abs(A[k] - A[k - 1]).sum() + np.abs(I[k] - I[k - 1]).sum()) / (2 * N * N) for k in range(1, len(A))])
    rows, snaps = [], {}
    for k in range(len(A)):
        t = frame_time(k)
        if k == 0 or (k - 1) % 5 == 0:                 # t = 0, 0.01, 10.01, 20.01, ... 600.01
            m = measure(A[k])
            rows.append(dict(bi=bi, condition=cond, replicate=rep, spike_seed=spike_seed, frame=k, t=round(t, 2), **m))
        if k >= 1 and round(t) in SNAP_T and abs(t - round(t)) < 0.05:
            snaps[int(round(t))] = A[k].astype(np.float32)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"bi{bi}_{cond}_rep{rep}.npz", crit=crit, crit_t=np.array([frame_time(k) for k in range(1, len(A))]),
                        a_ss=a_ss, i_ss=i_ss, **{f"a_t{t}": s for t, s in snaps.items()})
    return rows, dict(bi=bi, condition=cond, replicate=rep, a_ss=a_ss, i_ss=i_ss, wall_s=round(time.perf_counter() - t0, 1))


def run(n_workers=4):
    jobs = [(bi, cond, rep) for bi in BIS for cond in CONDITIONS for rep in REPS[cond]]
    print(f"{len(jobs)} runs on {n_workers} workers", flush=True)
    rows, meta = [], []
    with Pool(n_workers) as pool:
        for r, m in pool.imap_unordered(run_one, jobs, chunksize=1):
            rows.extend(r); meta.append(m); print("done", m, flush=True)
    rows.sort(key=lambda r: (r["bi"], list(CONDITIONS).index(r["condition"]), r["replicate"], r["frame"]))
    SUMM.mkdir(parents=True, exist_ok=True)
    with open(FRAMES_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    with open(SUMM / "seed_crowding_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(meta[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(sorted(meta, key=lambda m: (m["bi"], m["condition"], m["replicate"])))
    print(f"wrote {FRAMES_CSV}", flush=True)


def settle_time(series_t, series_a, series_b):
    """First t after which both integer series stay unchanged through the last frame (at least two frames); None if never."""
    last = len(series_t) - 1
    for j in range(last, -1, -1):
        if series_a[j] != series_a[last] or series_b[j] != series_b[last]:
            return series_t[j + 1] if j + 1 < last else None
    return series_t[0]


def analyze():
    rows = list(csv.DictReader(open(FRAMES_CSV)))
    for r in rows:
        for k in ("bi", "replicate", "frame", "cells_on", "blobs", "largest_blob"):
            r[k] = int(r[k])
        r["t"] = float(r["t"]); r["coverage"] = float(r["coverage"]); r["total_activator"] = float(r["total_activator"])
    runs = {}
    for r in rows:
        runs.setdefault((r["bi"], r["condition"], r["replicate"]), []).append(r)
    out = []
    for key, rs in sorted(runs.items(), key=lambda kv: (kv[0][0], list(CONDITIONS).index(kv[0][1]), kv[0][2])):
        rs = [r for r in rs if r["t"] >= 10]                     # measured frames t = 10 ... 600
        ts = [r["t"] for r in rs]
        st = settle_time(ts, [r["blobs"] for r in rs], [r["cells_on"] for r in rs])
        d = np.load(OUT / f"bi{key[0]}_{key[1]}_rep{key[2]}.npz")
        crit, ct = d["crit"], d["crit_t"]
        # published criterion: diff < 1e-4 at a save point with loop index > min_steps (1000 in 2D_batch/config.yaml)
        loop_index = np.arange(len(crit)) * SAVE_EVERY            # window k (crit[k-1]) is checked at loop index k*SAVE_EVERY
        hit = np.nonzero((crit < 1e-4) & (loop_index > 1000))[0]
        t_stop = float(ct[hit[0]]) if hit.size else None
        final = rs[-1]
        b = [r["blobs"] for r in rs]; c = [r["cells_on"] for r in rs]
        j = len(b) - 1
        while j > 0 and b[j - 1] == b[-1]:
            j -= 1
        settle_blobs = ts[j] if j < len(b) - 1 else None
        out.append(dict(bi=key[0], condition=key[1], replicate=key[2], settle_t=st if st is not None else "not settled",
                        settle_blobs_t=settle_blobs if settle_blobs is not None else "not settled", cells_drift_500_600=c[-1] - c[-11],
                        criterion_stop_t=t_stop if t_stop is not None else "never", crit_min=f"{crit[5:].min():.1e}", crit_final=f"{crit[-1]:.1e}",
                        blobs_final=final["blobs"], cells_on_final=final["cells_on"], largest_final=final["largest_blob"],
                        coverage_final=round(final["coverage"], 4), total_activator_final=round(final["total_activator"], 1),
                        blobs_t60=[r for r in rs if r["t"] < 61][-1]["blobs"], blobs_t300=[r for r in rs if r["t"] < 301][-1]["blobs"]))
    with open(SUMM / "seed_crowding_settling.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(out)
    print(f"{'bi':>3} {'condition':<8}{'rep':>4}{'settle':>12}{'blobs settle':>13}{'drift':>7}{'1e-4 stop':>11}{'crit min':>9}{'crit end':>9}{'blobs@60':>9}{'@300':>6}{'@600':>6}{'cells@600':>10}{'largest':>8}{'cov':>7}")
    for o in out:
        print(f"{o['bi']:>3} {o['condition']:<8}{o['replicate']:>4}{str(o['settle_t']):>12}{str(o['settle_blobs_t']):>13}{o['cells_drift_500_600']:>+7}{str(o['criterion_stop_t']):>11}{o['crit_min']:>9}{o['crit_final']:>9}"
              f"{o['blobs_t60']:>9}{o['blobs_t300']:>6}{o['blobs_final']:>6}{o['cells_on_final']:>10}{o['largest_final']:>8}{o['coverage_final']:>7.3f}")
    # figures: one grid per parameter set, rows = conditions (replicate 1 or 0), columns = snapshot times
    from analysis.plot_hex import save_panels
    for bi in BIS:
        fields, titles = [], []
        for cond in CONDITIONS:
            rep = REPS[cond][0]
            d = np.load(OUT / f"bi{bi}_{cond}_rep{rep}.npz")
            for t in SNAP_T:
                fields.append(d[f"a_t{t}"]); titles.append(f"{cond} rep{rep}, t={t}")
        save_panels(fields, titles, SUMM / f"seed_crowding_bi{bi}.png", ncols=len(SNAP_T),
                    suptitle=f"b_a=5, b_i={bi}, n=3/3, D_i=10, gamma=0.5; 141x141, threshold {THR}; activator, hex-correct")
    print("figures written")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif sys.argv[1] == "analyze":
        analyze()
    else:
        raise SystemExit(__doc__)
