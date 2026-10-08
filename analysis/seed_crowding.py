"""
Does the fate of one isolated seed predict what a whole field does, and how long does each take to settle?

Two parameter groups (CLI: `run <group> [workers]`, `analyze <group>`):
  fig1c       Fig 1C line (b_a 5, n 3/3, D_i 10, gamma 0.5) at b_i = 5, 9, 12. Seed level 5.0; synchronous start U(0,1).
  expmatched  the experiment-matched Supp Table 1 sets (n 10/4, gamma 0.5): 2L JAPI and JA, 3B labyrinth/big/small,
              3H low D_i, 3L, S8D low r_i. A single cell at level 5 dies on every one of these sets and level 10
              nucleates on all of them (domain identical at 10, 20, 40), so seeded modalities use level 10 and one
              extra lone-seed run at level 5 documents the failure. Synchronous start U(0,2) (U(0,1) does not
              nucleate at n_a = 10; REPORT s5). Lone seeds run to t = 1000, fields to t = 600.
Conditions per set: lone seed, 9 / 36 / 144 seeds (activator_random_spikes, peak = seed level), random
(random_uniform_over0, no nucleation) and nucleation (all_off, spike 2.0, rate 0.01, kick 4.0: the published
Fig 2L/3 protocol). Three replicate seeds for everything but the lone seeds. 141x141, dt 0.01, convergence check
disabled (stopping_threshold 0, min_steps > steps). The engine saves every 200 steps so its own convergence
diagnostic (CLAUDE.md s8) is recorded for every window; measurements are taken at t = 10, 20, ..., end:
cells above 0.3*r_a, blob count on the engine's own neighbour offsets, largest blob, coverage, total activator.

Engine caveats (CLAUDE.md s7): activator_random_spikes places its seeds from a hard-coded init_seed, so the worker
overrides numpy.random.default_rng during the call (replicate r -> seed 1000*r + 7, recorded as spike_seed); the
lone seed uses the same mode with n_points = 1 and a seed chosen so the cell lies within 3 cells of the centre
(activator only, no inhibitor preload). In the fig1c group the lone seed is spike_steady_state at 5.0 (a = i = 5
at the centre), as run on 2026-10-08. Everything else is seeded with np.random.seed(replicate).
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

N, DT, SAVE_EVERY = 141, 0.01, 200
SNAP_T = (60, 150, 300, 600)
SUMM = ROOT / "analysis" / "summaries"


def circuit(ba, bi, g, na, ni, D):
    return dict(act_half_sat=1.0, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=1.0, act_prod_rate=float(ba),
                inh_prod_rate=float(bi), inh_decay_rate=float(g), act_hill_coeff=na, inh_hill_coeff=ni, inh_diffusion=float(D))


GROUPS = {
    "fig1c": dict(
        sets={"bi5": circuit(5, 5, 0.5, 3, 3, 10), "bi9": circuit(5, 9, 0.5, 3, 3, 10), "bi12": circuit(5, 12, 0.5, 3, 3, 10)},
        seed_level=5.0, sync_spike=1.0, t_field=600.0, t_lone=600.0, lone_mode="spike_steady_state", extra_lone_levels=(),
        results="seed_crowding", prefix="seed_crowding"),
    "expmatched": dict(
        sets={"2L_JAPI": circuit(5, 25, 0.5, 10, 4, 20), "2L_JA": circuit(5, 0, 0.5, 10, 4, 20), "3B_labyrinth": circuit(6, 15, 0.5, 10, 4, 20),
              "3B_big": circuit(5, 15, 0.5, 10, 4, 20), "3B_small": circuit(5, 30, 0.5, 10, 4, 20), "3H_lowDi": circuit(5, 25, 0.5, 10, 4, 5),
              "3L": circuit(5, 20, 0.5, 10, 4, 20), "S8D_lowri": circuit(5, 50, 0.5, 10, 4, 20)},
        seed_level=10.0, sync_spike=2.0, t_field=600.0, t_lone=1000.0, lone_mode="activator_random_spikes", extra_lone_levels=(5.0,),
        results="seed_crowding_exp", prefix="seed_crowding_exp"),
}
REPS = {"seed1": [0], "seed1_L5": [0], "seed9": [1, 2, 3], "seed36": [1, 2, 3], "seed144": [1, 2, 3], "random": [1, 2, 3], "nucl": [1, 2, 3]}


def conditions(group):
    G = GROUPS[group]; L = G["seed_level"]
    c = {"seed1": dict(init_mode=G["lone_mode"], spike_value=L, n_points=1, set_peak_height=L, nucleation_rate=0.0, t_end=G["t_lone"])}
    for extra in G["extra_lone_levels"]:
        c[f"seed1_L{extra:g}"] = dict(init_mode=G["lone_mode"], spike_value=extra, n_points=1, set_peak_height=extra, nucleation_rate=0.0, t_end=G["t_lone"])
    for n in (9, 36, 144):
        c[f"seed{n}"] = dict(init_mode="activator_random_spikes", spike_value=L, n_points=n, set_peak_height=L, nucleation_rate=0.0, t_end=G["t_field"])
    c["random"] = dict(init_mode="random_uniform_over0", spike_value=G["sync_spike"], n_points=0, set_peak_height=L, nucleation_rate=0.0, t_end=G["t_field"])
    c["nucl"] = dict(init_mode="all_off", spike_value=2.0, n_points=0, set_peak_height=L, nucleation_rate=0.01, t_end=G["t_field"])
    return c


def central_seed(n, tol=3):
    """Replicate seed for which the engine's activator_random_spikes places its single cell within tol of the centre."""
    for s in range(1, 20000):
        idx = int(np.random.default_rng(s).choice(n * n, size=1, replace=False)[0])
        r, c = divmod(idx, n)
        if abs(r - n // 2) <= tol and abs(c - n // 2) <= tol:
            return s
    raise RuntimeError("no central seed found")


def frame_time(k):
    """Frame k (k >= 1) is the state after (k-1)*SAVE_EVERY + 1 Euler steps (CLAUDE.md s8); frame 0 is t = 0."""
    return 0.0 if k == 0 else ((k - 1) * SAVE_EVERY + 1) * DT


def measure(a, thr):
    mask = a > thr
    labels, n = label_components(mask) if mask.any() else (None, 0)
    largest = int(np.bincount(labels[labels >= 0]).max()) if n else 0
    return dict(cells_on=int(mask.sum()), blobs=n, largest_blob=largest, coverage=float(mask.mean()), total_activator=float(a.sum()))


def run_one(job):
    group, setname, cond, rep = job
    from simulation_2D import run_coupled_hex
    G = GROUPS[group]; p = G["sets"][setname]; c = conditions(group)[cond]
    thr = 0.3 * p["act_prod_rate"]
    steps = int(round(c["t_end"] / DT)) + 1
    out = ROOT / "results" / G["results"]
    spike_seed = ""
    np.random.seed(rep)
    orig_rng = np.random.default_rng
    if c["init_mode"] == "activator_random_spikes":
        spike_seed = central_seed(N) if c["n_points"] == 1 else 1000 * rep + 7
        np.random.default_rng = lambda seed=None, _s=spike_seed: orig_rng(_s)   # override the engine's hard-coded init_seed
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            A, I, step, a_ss, i_ss = run_coupled_hex(N, N, steps, DT, 1.0, p, 0.0, 10**9, init_mode=c["init_mode"], activator_type="juxtacrine",
                                                     n_points=c["n_points"], spike_value=c["spike_value"], save_every=SAVE_EVERY,
                                                     nucleation_rate=c["nucleation_rate"], noise_amplitude=0.0, set_peak_height=c["set_peak_height"])
    finally:
        np.random.default_rng = orig_rng
    assert len(A) == steps // SAVE_EVERY + 2 and abs(frame_time(len(A) - 1) - (c["t_end"] + DT)) < 1e-9
    crit = np.array([(np.abs(A[k] - A[k - 1]).sum() + np.abs(I[k] - I[k - 1]).sum()) / (2 * N * N) for k in range(1, len(A))])
    rows, snaps = [], {}
    for k in range(len(A)):
        t = frame_time(k)
        if k == 0 or (k - 1) % 5 == 0:
            rows.append(dict(set=setname, bi=p["inh_prod_rate"], condition=cond, replicate=rep, spike_seed=spike_seed, frame=k, t=round(t, 2), **measure(A[k], thr)))
        if k >= 1 and round(t) in SNAP_T + (1000,) and abs(t - round(t)) < 0.05:
            snaps[int(round(t))] = A[k].astype(np.float32)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / f"{setname}_{cond}_rep{rep}.npz", crit=crit, crit_t=np.array([frame_time(k) for k in range(1, len(A))]),
                        a_ss=a_ss, i_ss=i_ss, threshold=thr, **{f"a_t{t}": s for t, s in snaps.items()})
    return rows, dict(set=setname, condition=cond, replicate=rep, a_ss=a_ss, i_ss=i_ss, t_end=c["t_end"], seed_level=c["set_peak_height"],
                      wall_s=round(time.perf_counter() - t0, 1))


def run(group, n_workers=4):
    G = GROUPS[group]; conds = conditions(group)
    jobs = [(group, s, cond, rep) for s in G["sets"] for cond in conds for rep in REPS[cond]]
    print(f"{len(jobs)} runs on {n_workers} workers", flush=True)
    rows, meta = [], []
    with Pool(n_workers) as pool:
        for r, m in pool.imap_unordered(run_one, jobs, chunksize=1):
            rows.extend(r); meta.append(m); print("done", m, flush=True)
    order = {s: i for i, s in enumerate(G["sets"])}; corder = {c: i for i, c in enumerate(conds)}
    rows.sort(key=lambda r: (order[r["set"]], corder[r["condition"]], r["replicate"], r["frame"]))
    SUMM.mkdir(parents=True, exist_ok=True)
    with open(SUMM / f"{G['prefix']}_frames.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    with open(SUMM / f"{G['prefix']}_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(meta[0].keys()), lineterminator="\n"); w.writeheader()
        w.writerows(sorted(meta, key=lambda m: (order[m["set"]], corder[m["condition"]], m["replicate"])))
    print(f"wrote {SUMM / (G['prefix'] + '_frames.csv')}", flush=True)


def settle_time(series_t, series_a, series_b):
    """First t after which both integer series stay unchanged through the last frame (at least two frames); None if never."""
    last = len(series_t) - 1
    for j in range(last, -1, -1):
        if series_a[j] != series_a[last] or series_b[j] != series_b[last]:
            return series_t[j + 1] if j + 1 < last else None
    return series_t[0]


def analyze(group):
    G = GROUPS[group]; conds = list(conditions(group)); out_dir = ROOT / "results" / G["results"]
    rows = list(csv.DictReader(open(SUMM / f"{G['prefix']}_frames.csv")))
    for r in rows:
        r.setdefault("set", f"bi{int(float(r['bi']))}")
        for k in ("replicate", "frame", "cells_on", "blobs", "largest_blob"):
            r[k] = int(r[k])
        r["bi"] = float(r["bi"]); r["t"] = float(r["t"]); r["coverage"] = float(r["coverage"]); r["total_activator"] = float(r["total_activator"])
    runs = {}
    for r in rows:
        runs.setdefault((r["set"], r["condition"], r["replicate"]), []).append(r)
    order = {s: i for i, s in enumerate(G["sets"])}
    out = []
    for key, rs in sorted(runs.items(), key=lambda kv: (order.get(kv[0][0], 99), conds.index(kv[0][1]), kv[0][2])):
        rs = [r for r in rs if r["t"] >= 10]
        ts = [r["t"] for r in rs]; b = [r["blobs"] for r in rs]; c = [r["cells_on"] for r in rs]
        st = settle_time(ts, b, c)
        j = len(b) - 1
        while j > 0 and b[j - 1] == b[-1]:
            j -= 1
        settle_blobs = ts[j] if j < len(b) - 1 else None
        fname = out_dir / f"{key[0]}_{key[1]}_rep{key[2]}.npz"
        if not fname.exists():                                   # fig1c group was run with bi-named files
            fname = out_dir / f"bi{int(rs[0]['bi'])}_{key[1]}_rep{key[2]}.npz"
        d = np.load(fname); crit, ct = d["crit"], d["crit_t"]
        loop_index = np.arange(len(crit)) * SAVE_EVERY
        hit = np.nonzero((crit < 1e-4) & (loop_index > 1000))[0]
        final = rs[-1]
        at = lambda tt: [r for r in rs if r["t"] < tt + 1][-1]
        out.append(dict(set=key[0], condition=key[1], replicate=key[2], t_end=ts[-1], settle_t=st if st is not None else "not settled",
                        settle_blobs_t=settle_blobs if settle_blobs is not None else "not settled", cells_drift_last100=c[-1] - c[-11],
                        criterion_stop_t=float(ct[hit[0]]) if hit.size else "never", crit_min=f"{crit[5:].min():.1e}", crit_final=f"{crit[-1]:.1e}",
                        blobs_t60=at(60)["blobs"], blobs_t300=at(300)["blobs"], blobs_t600=at(600)["blobs"], blobs_final=final["blobs"],
                        cells_on_final=final["cells_on"], largest_final=final["largest_blob"], coverage_final=round(final["coverage"], 4),
                        total_activator_final=round(final["total_activator"], 1)))
    with open(SUMM / f"{G['prefix']}_settling.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(out)
    print(f"{'set':<13}{'condition':<9}{'rep':>4}{'t_end':>6}{'settle':>12}{'blobs settle':>13}{'drift':>7}{'1e-4 stop':>11}{'crit min':>9}{'crit end':>9}{'b@60':>6}{'b@300':>6}{'b@600':>6}{'b@end':>6}{'cells':>7}{'largest':>8}{'cov':>7}")
    for o in out:
        print(f"{o['set']:<13}{o['condition']:<9}{o['replicate']:>4}{o['t_end']:>6.0f}{str(o['settle_t']):>12}{str(o['settle_blobs_t']):>13}{o['cells_drift_last100']:>+7}{str(o['criterion_stop_t']):>11}"
              f"{o['crit_min']:>9}{o['crit_final']:>9}{o['blobs_t60']:>6}{o['blobs_t300']:>6}{o['blobs_t600']:>6}{o['blobs_final']:>6}{o['cells_on_final']:>7}{o['largest_final']:>8}{o['coverage_final']:>7.3f}")
    from analysis.plot_hex import save_panels
    for setname in G["sets"]:
        fields, titles = [], []
        for cond in conds:
            rep = REPS[cond][0]
            fname = out_dir / f"{setname}_{cond}_rep{rep}.npz"
            if not fname.exists():
                fname = out_dir / f"bi{int(G['sets'][setname]['inh_prod_rate'])}_{cond}_rep{rep}.npz"
            d = np.load(fname)
            for t in SNAP_T:
                fields.append(d[f"a_t{t}"]); titles.append(f"{cond} rep{rep}, t={t}")
        p = G["sets"][setname]
        save_panels(fields, titles, SUMM / f"{G['prefix']}_{setname}.png", ncols=len(SNAP_T),
                    suptitle=f"{setname}: b_a={p['act_prod_rate']:g}, b_i={p['inh_prod_rate']:g}, n={p['act_hill_coeff']}/{p['inh_hill_coeff']}, D_i={p['inh_diffusion']:g}, gamma={p['inh_decay_rate']:g}; "
                             f"141x141, threshold {0.3*p['act_prod_rate']:.2g}; activator, hex-correct")
    print("figures written")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4)
    elif sys.argv[1] == "analyze":
        analyze(sys.argv[2])
    else:
        raise SystemExit(__doc__)
