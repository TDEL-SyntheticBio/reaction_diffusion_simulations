"""
2D initiation sweep (JAPI): how the architecture turns different starts into discrete spots.

Architecture study for Fig 1. Nothing here is anchored to measured values.

AXES
  Hill pair   n_a, n_i each in {2, 4, 10}                                   9 pairs
  edge pos    b_i = f x edge(n_a, n_i), f in {0.1, 0.2, 0.4, 0.7, 1, 1.4, 2, 3}   8 values
              edge = the b_i at which the uniform activated state stops existing
              (a single-cell property: same in 1D and 2D, same for JAPI and PAPI)
  D_i         5, 10, 20 (the values used in the paper's 2D panels)        3 values
  start       seed1     one cell at the centre set to H_ON                 1 run
              seeds     SEED_DENSITY single cells at H_ON, random positions 3 replicates
              rf1, rf2  activator uniform in [0, M x act_half_sat], M = 1, 2 3 replicates each
              nuc01, nuc001  nucleation: each step every cell below b_a is set to H_ON with
                        probability rate x dt, rate = 0.01, 0.001          3 replicates each
FIXED
  b_a 3, gamma 0.5, act_half_sat = inh_half_sat = 1, act_decay 1, basal 0, juxtacrine,
  H_ON = 4 b_a = 12 (above the 2D single-cell ignition threshold at every Hill pair: 6.0 at
  n_a 2, 8.8 at n_a 4, 7.4 at n_a 10, measured on 2026-10-16; the one exception is n 10/2
  at and beyond the edge: no ignition from a single cell of 12 at 1.0x, none even at 80 at
  1.4x, only a 61-cell domain at 0.4x), 200 x 200 hex, dt 0.01, t = 50.
  Inhibitor starts at zero everywhere.

  Random streams are common across parameters: replicate r of a start uses the same seed
  positions / noise pattern / kick stream at every Hill pair, edge position and D_i, so
  differences between parameter sets are not differences between random draws.

DEVIATIONS FROM THE DEPOSIT (deliberate, all recorded here)
  - Nucleation sets a fired cell to H_ON if it is below b_a. The engine adds 2 a_ss if below
    a_ss, where a_ss comes from the steady-state solver and falls back to spike_value when the
    solver finds no state, i.e. it jumps exactly at the edge this sweep straddles.
  - Seeds are built here; run_coupled_hex draws seed positions from a hard-coded generator.
  - The integration loop is ours, calling the engine's _vectorized_step and
    _build_hex_neighbor_arrays unchanged, so the last frame is the true final state.

MEASUREMENTS (every t = 1, on the activator, level 0.3 b_a, engine hex adjacency)
  spot     connected component of >= 7 cells
  stable   a spot at t that overlaps a spot at t - 5 (failed kicks of 12 fall below 0.3 b_a
           within about 3 time units, so a 5-unit persistence window excludes them)
  label at t = 50: dead (no spot) | sheet (coverage >= 0.9) | network (largest component
           spans the field between opposite borders) | spots (otherwise)
  speed    t_first: first t with a stable spot; t90: first t at which the stable-spot count
           reaches 90% of its value at t = 50
  size     median spot size in cells and median equivalent diameter / lambda_i,
           lambda_i = sqrt(1.5 D_i / gamma) (2D, raw six-neighbour Laplacian)
  uniform  nearest-neighbour spacing CV of spot centroids divided by the same statistic for the
           same number of uniform random points in the same window (1 = random, < 1 more even);
           only with >= 10 spots. Recorded at t = 50 and at the day-4 equivalent.
Saved per run: the time series above, activator at t ~ days 1-4 equivalents (1 unit ~ 5.4 h),
t = 30 and t = 50 (float16), inhibitor at t = 50.

USAGE
  python analysis/sweep2d_initiation.py edges
  python analysis/sweep2d_checks.py                       # trap checks and timing, run first
  python analysis/sweep2d_initiation.py run [workers]     # full sweep, resumable
  python analysis/sweep2d_initiation.py aggregate         # -> analysis/summaries/sweep2d_*.csv
Env: SWEEP2D_OUT (default results/sweep2d), SWEEP2D_N (default 200).
"""
from __future__ import annotations

import itertools, json, os, sys, time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "2D_simulations"))
sys.path.insert(0, str(ROOT / "analysis"))
from simulation_2D import _vectorized_step, _build_hex_neighbor_arrays  # noqa: E402
from hexgeom import label_components, cell_centres, window, CELL_AREA     # noqa: E402
from pointstats import nn_cv, csr_reference_cv                             # noqa: E402

# ---------------------------------------------------------------- configuration
BA, GA, KA = 3.0, 0.5, 1.0
H_ON = 4.0 * BA
THR = 0.3 * BA
MIN_SPOT = 7
PERSIST = 5.0
DT, T_END = 0.01, 50.0
N_DEFAULT = int(os.environ.get("SWEEP2D_N", 200))
OUT = Path(os.environ.get("SWEEP2D_OUT", ROOT / "results" / "sweep2d"))

HILL = [(na, ni) for na in (2, 4, 10) for ni in (2, 4, 10)]
FRACS = [0.1, 0.2, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0]
DIS = [5.0, 10.0, 20.0]
SEED_DENSITY = 70 / 40000.0                    # seeds per cell: 70 on a 200 x 200 field
STARTS = {                                     # name: (id, replicates)
    "seed1": (0, 1), "seeds": (1, 3), "rf1": (2, 3), "rf2": (3, 3), "nuc01": (4, 3), "nuc001": (5, 3),
}
RATES = {"nuc01": 0.01, "nuc001": 0.001}
RF_MAX = {"rf1": 1.0 * KA, "rf2": 2.0 * KA}
DAY = 24.0 / 5.4                               # model units per day (1 unit ~ 5.4 h)
SAVE_T = [DAY, 2 * DAY, 3 * DAY, 4 * DAY, 30.0, 50.0]


# ---------------------------------------------------------------- the edge
def stable_state(ba, bi, ga, na, ni):
    """Largest stable activated root of the uniform equations, or None."""
    c = bi / (ga * ba)

    def g(u):
        a = u * ba
        A = a ** na
        I = (c * a) ** ni
        return A / (A + I + 1.0) - u
    U = np.concatenate([np.logspace(-6, -2, 200), np.linspace(0.0101, 1.0, 3000)])
    G = np.array([g(u) for u in U])
    roots = [brentq(g, U[k], U[k + 1]) * ba for k in range(len(U) - 1) if np.sign(G[k]) != np.sign(G[k + 1])]
    for a in sorted(roots, reverse=True):
        if a <= 1e-5:
            continue
        i = c * a
        A = a ** na
        I = i ** ni
        S = A + I + 1.0
        dHa = na * A / a * (I + 1.0) / S ** 2
        dHi = -ni * I / i * A / S ** 2
        j00, j01, j10, j11 = ba * dHa - 1.0, ba * dHi, bi * dHa, bi * dHi - ga
        if j00 + j11 < 0 and j00 * j11 - j01 * j10 > 0:
            return a
    return None


def edge_bi(na, ni, ba=BA, ga=GA):
    grid = np.logspace(-4, 5, 400)
    has = [stable_state(ba, b, ga, na, ni) is not None for b in grid]
    k = max(i for i, h in enumerate(has) if h)
    lo, hi = grid[k], grid[k + 1]
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        if stable_state(ba, mid, ga, na, ni) is not None:
            lo = mid
        else:
            hi = mid
    return float(lo)


_EDGES: dict = {}


def edges():
    if not _EDGES:
        f = OUT / "edges.json"
        if f.exists():
            _EDGES.update({tuple(map(int, k.split("_"))): v for k, v in json.loads(f.read_text()).items()})
        else:
            _EDGES.update({pr: edge_bi(*pr) for pr in HILL})
            OUT.mkdir(parents=True, exist_ok=True)
            f.write_text(json.dumps({f"{k[0]}_{k[1]}": v for k, v in _EDGES.items()}, indent=1))
    return _EDGES


# ---------------------------------------------------------------- jobs
def all_jobs():
    """Replicate is the outermost loop: every condition once, then the second replicates, then the
    third. A sweep stopped early still gives complete maps, with fewer replicates."""
    jobs = []
    for rep in range(max(n for _, n in STARTS.values())):
        for (na, ni), fi, di in itertools.product(HILL, range(len(FRACS)), range(len(DIS))):
            for start, (_, nrep) in STARTS.items():
                if rep < nrep:
                    jobs.append((na, ni, fi, di, start, rep))
    return jobs


def job_name(job):
    na, ni, fi, di, start, rep = job
    return f"na{na:02d}_ni{ni:02d}_f{fi}_D{int(DIS[di]):02d}_{start}_r{rep}"


def params(na, ni, bi, D):
    return dict(act_half_sat=KA, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=0.0,
                inh_diffusion=float(D), act_prod_rate=BA, inh_prod_rate=float(bi), inh_decay_rate=GA,
                act_hill_coeff=na, inh_hill_coeff=ni)


def stream(start, rep):
    """Common random stream for (start, replicate), independent of the parameters."""
    return np.random.default_rng(np.random.SeedSequence([20261016, STARTS[start][0], rep]))


def initial(start, rng, ny, nx):
    a = np.zeros((ny, nx))
    if start == "seed1":
        a[ny // 2, nx // 2] = H_ON
    elif start == "seeds":
        n = max(1, int(round(SEED_DENSITY * ny * nx)))
        idx = rng.choice(ny * nx, size=n, replace=False)
        a.ravel()[idx] = H_ON
    elif start in RF_MAX:
        a = rng.random((ny, nx)) * RF_MAX[start]
    return a, np.zeros((ny, nx))


def kick(a, rng, p_fire):
    """Nucleation, in place: every cell below b_a fires with probability p_fire and is set to H_ON."""
    fired = (rng.random(a.shape) < p_fire) & (a < BA)
    if fired.any():
        a[fired] = H_ON
    return int(fired.sum())


# ---------------------------------------------------------------- measurement helpers
_GEOM: dict = {}


def geom(ny, nx):
    if (ny, nx) not in _GEOM:
        x, y = cell_centres(ny, nx)
        _GEOM[(ny, nx)] = (_build_hex_neighbor_arrays(ny, nx), x, y, window(ny, nx))
    return _GEOM[(ny, nx)]


def spots(a, ny, nx):
    """Returns (labels, list of component index arrays >= MIN_SPOT, coverage, largest-component mask)."""
    m = a > THR
    lab, n = label_components(m)
    if n == 0:
        return lab, [], float(m.mean()), None
    sizes = np.bincount(lab[lab >= 0], minlength=n)
    big = np.nonzero(sizes >= MIN_SPOT)[0]
    largest = lab == int(np.argmax(sizes))
    return lab, [lab == k for k in big], float(m.mean()), largest


def spans(mask):
    r = np.nonzero(mask.any(axis=1))[0]
    c = np.nonzero(mask.any(axis=0))[0]
    ny, nx = mask.shape
    return bool((r.min() == 0 and r.max() == ny - 1) or (c.min() == 0 and c.max() == nx - 1))


def uniformity(comps, x, y, win, rng):
    if len(comps) < 10:
        return np.nan, len(comps)
    xy = np.array([[x[m].mean(), y[m].mean()] for m in comps])
    ref, _ = csr_reference_cv(len(xy), win, rng, n_real=100)
    return float(nn_cv(xy) / ref), len(comps)


# ---------------------------------------------------------------- one run
def run_one(job, ny=None, nx=None, save=True, out=None):
    ny = ny or N_DEFAULT
    nx = nx or N_DEFAULT
    out = Path(out or OUT)
    na, ni, fi, di, start, rep = job
    D = DIS[di]
    bi = FRACS[fi] * edges()[(na, ni)]
    p = params(na, ni, bi, D)
    (nr, nc, nm, ncnt), x, y, win = geom(ny, nx)
    rng = stream(start, rep)
    a, i = initial(start, rng, ny, nx)
    rate = RATES.get(start, 0.0)
    p_fire = rate * DT
    steps = int(round(T_END / DT))
    per = int(round(1.0 / DT))
    save_steps = {int(round(t / DT)): t for t in SAVE_T}
    lag = int(round(PERSIST))
    ts, cov, nspot, nstable, largest_frac = [], [], [], [], []
    hist_masks = {}
    frames, n_kicks = {}, 0
    imin, imax = 0.0, 0.0
    u_day4, n_day4 = np.nan, 0
    t0 = time.time()
    for s in range(1, steps + 1):
        a, i = _vectorized_step(a, i, DT, 1.0, p, "juxtacrine", nr, nc, nm, ncnt)
        if p_fire > 0.0:
            n_kicks += kick(a, rng, p_fire)
        if s in save_steps:
            frames[save_steps[s]] = a.astype(np.float16)
        if s % per == 0:
            t = s // per
            imin = min(imin, float(i.min())); imax = max(imax, float(i.max()))
            _, comps, c, largest = spots(a, ny, nx)
            allmask = np.zeros((ny, nx), bool)
            for m in comps:
                allmask |= m
            prev = hist_masks.get(t - lag)
            ns = int(sum(1 for m in comps if prev is not None and (m & prev).any())) if prev is not None else 0
            hist_masks[t] = allmask
            hist_masks.pop(t - lag - 1, None)
            ts.append(t); cov.append(c); nspot.append(len(comps)); nstable.append(ns)
            largest_frac.append(float(largest.sum() / max((a > THR).sum(), 1)) if largest is not None else 0.0)
            if t == int(round(4 * DAY)):
                u_day4, n_day4 = uniformity(comps, x, y, win, np.random.default_rng(1))
    # final measurements
    _, comps, c50, largest = spots(a, ny, nx)
    if not comps:
        label = "dead"
    elif c50 >= 0.9:
        label = "sheet"
    elif largest is not None and largest.sum() >= MIN_SPOT and spans(largest):
        label = "network"
    else:
        label = "spots"
    sizes = np.array([m.sum() for m in comps]) if comps else np.array([])
    lam = np.sqrt(1.5 * D / GA)
    eqd = 2.0 * np.sqrt(sizes * CELL_AREA / np.pi) if len(sizes) else np.array([])
    u50, n50 = uniformity(comps, x, y, win, np.random.default_rng(1))
    nst = np.array(nstable)
    final_st = int(nst[-1])
    t_first = float(np.array(ts)[nst > 0][0]) if (nst > 0).any() else np.nan
    t90 = float(np.array(ts)[nst >= 0.9 * final_st][0]) if final_st > 0 else np.nan
    row = dict(name=job_name(job), n_a=na, n_i=ni, f_edge=FRACS[fi], b_i=bi, D_i=D, start=start, rep=rep,
               label=label, coverage=c50, n_spots=len(comps), n_stable=final_st,
               largest_frac=largest_frac[-1], t_first=t_first, t90=t90,
               size_median=float(np.median(sizes)) if len(sizes) else np.nan,
               size_q25=float(np.percentile(sizes, 25)) if len(sizes) else np.nan,
               size_q75=float(np.percentile(sizes, 75)) if len(sizes) else np.nan,
               diam_over_lambda=float(np.median(eqd) / lam) if len(eqd) else np.nan, lambda_i=lam,
               uniformity_t50=u50, uniformity_day4=u_day4, n_spots_day4=n_day4,
               n_kicks=n_kicks, i_min=imin, i_max=imax, a_max=float(a.max()),
               seconds=time.time() - t0)
    if save:
        (out / "rows").mkdir(parents=True, exist_ok=True)
        (out / "frames").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out / "frames" / f"{row['name']}.npz",
                            t=np.array(ts), coverage=np.array(cov), n_spots=np.array(nspot),
                            n_stable=np.array(nstable), largest_frac=np.array(largest_frac),
                            saved_t=np.array(sorted(frames)), a=np.array([frames[k] for k in sorted(frames)]),
                            i_final=i.astype(np.float16))
        (out / "rows" / f"{row['name']}.json").write_text(json.dumps(row))
    return row


def _run_safe(job):
    try:
        return run_one(job)
    except Exception as e:                    # keep the sweep going; record the failure
        (OUT / "failed").mkdir(parents=True, exist_ok=True)
        (OUT / "failed" / f"{job_name(job)}.txt").write_text(repr(e))
        return None


# ---------------------------------------------------------------- commands
def cmd_run(workers):
    from multiprocessing import Pool
    edges()
    jobs = [j for j in all_jobs() if not (OUT / "rows" / f"{job_name(j)}.json").exists()]
    print(f"{len(all_jobs())} jobs, {len(jobs)} still to run, {workers} workers -> {OUT}", flush=True)
    t0 = time.time()
    done = 0
    with Pool(workers) as pool:
        for r in pool.imap_unordered(_run_safe, jobs, chunksize=1):
            done += 1
            if done % 50 == 0 or done == len(jobs):
                el = time.time() - t0
                print(f"  {done}/{len(jobs)}  {el/3600:.2f} h elapsed, ~{el/done*(len(jobs)-done)/3600:.2f} h left",
                      flush=True)


def cmd_aggregate():
    import pandas as pd
    rows = [json.loads(f.read_text()) for f in sorted((OUT / "rows").glob("*.json"))]
    df = pd.DataFrame(rows)
    summ = ROOT / "analysis" / "summaries"
    summ.mkdir(parents=True, exist_ok=True)
    df.to_csv(summ / "sweep2d_runs.csv", index=False)
    g = (df.groupby(["n_a", "n_i", "f_edge", "D_i", "start"])
           .agg(runs=("label", "size"), dead=("label", lambda s: int((s == "dead").sum())),
                sheet=("label", lambda s: int((s == "sheet").sum())),
                network=("label", lambda s: int((s == "network").sum())),
                spots=("label", lambda s: int((s == "spots").sum())),
                n_stable=("n_stable", "median"), t_first=("t_first", "median"), t90=("t90", "median"),
                size_median=("size_median", "median"), diam_over_lambda=("diam_over_lambda", "median"),
                uniformity_t50=("uniformity_t50", "median"), uniformity_day4=("uniformity_day4", "median"))
           .reset_index())
    g.to_csv(summ / "sweep2d_by_condition.csv", index=False)
    nf = len(list((OUT / "failed").glob("*.txt"))) if (OUT / "failed").exists() else 0
    print(f"{len(df)} runs aggregated, {nf} failed -> {summ}/sweep2d_runs.csv, sweep2d_by_condition.csv")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "edges":
        for k, v in edges().items():
            print(f"n {k[0]}/{k[1]}: edge b_i = {v:.4f}")
    elif cmd == "run":
        cmd_run(int(sys.argv[2]) if len(sys.argv) > 2 else (os.cpu_count() or 2))
    elif cmd == "aggregate":
        cmd_aggregate()
    else:
        raise SystemExit(f"unknown command {cmd}; checks and timing live in analysis/sweep2d_checks.py")
