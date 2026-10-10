"""
Overnight batch (10-11 October 2026): the ladder design moved to the Fig 1C line, b_a 5, n_a = n_i = 3, D_i 10, gamma 0.5
(lambda_i = sqrt(1.5 D_i / gamma) = 5.477 cells throughout), seven b_i values and five initial conditions, 200x200, dt 0.01, to
t = 600, convergence check disabled, frames every t = 10, both channels. Everything is recorded; nothing is interpreted here.

Sets (b_i): 0 (JA, no inhibitor), 2.5 (inside the Turing band), 3.2 (Turing band, near its edge), 5 (replicating; the published
"regular spots" panel), 9 (divides late), 12 (arrests; the published "irregular spots" panel), 20 (arrests, small domains).
Conditions: seed1 (one activator-only cell at level 10 near the centre, replicate 0), seed36 (activator_random_spikes at level
10, placement seed overridden per replicate), random (random_uniform_over0, U(0, 1): the Supp Table 1 amplitude), nucl (all_off,
nucleation 0.01, kick 2 a_ss gated at a < a_ss, a_ss from the engine's solver or its fallback = spike_value 2), smallpert
(random_tight: both fields uniform within +-5 % of the activated steady state; the state is the upper-branch root of the closed
form (analysis.regime_map.activated_states), handed to the engine by replacing its solver for this condition only; at b_i 0 the
inhibitor state is 0 and is passed as 1e-12 so the engine does not fall back; where no activated state exists the engine's own
fallback a_ss = i_ss = spike_value = 2 is used and recorded as a_ss_source = fallback).  Replicates 1-3 on the stochastic
conditions; np.random.seed(replicate) before every engine call.

Measurements: the ladder driver's metric set (analysis/ladder_batch.measure_run: per-frame counts, coverage, areas, NN spacing
with the same-n CSR reference, eccentricity, convergence diagnostic; g(r) and hex-correct power spectra at t = 50, 100, 200, 400,
600; domain tracking; both channels, activator at the fixed level 0.3 r_a with persistence over the frame and the t - 10 frame,
inhibitor at half its own frame maximum, real domains >= 6 cells, transients smaller) plus, new here, a block of measures at
THREE time points per run: the first recorded frame (t >= 10) whose persistent activator coverage reaches 10 %, the first reaching 20 %,
and t = 600 (the activator's crossing frames are used for both channels; each channel's own crossing frames are recorded too):
  block_domains   count (real, transient), coverage, equivalent radius mean/median in cells and divided by lambda_i,
                  area distribution mean/SD/CV/median/p10/p90 (real and all), areas list
  block_spacing   NN spacing mean, SD, CV, Hanisch border-corrected CV, with three matched nulls at the measured count in the
                  same window, 200 realisations each: (a) CSR points, (b) hard-core discs with the measured per-domain
                  equivalent radii (random sequential placement, centre distance >= r_i + r_j), (c) hard-core points with
                  minimum separation = the field's own g(r) hole radius; each null's mean and SD of the CV and the z-score
                  (measured - null mean) / null SD, plain and border-corrected; the number of realisations that could be placed
  block_order     psi6 bond-orientational order (Delaunay = Voronoi neighbours, interior points: farther than the mean NN
                  spacing from the window edge) and the fraction of interior domains with exactly six Voronoi neighbours, with
                  the same three nulls (mean, SD, z)
  block_gr        g(r) of the field (PooledG, CSR null, lattice-site null when the median real domain is under 10 cells):
                  hole radius, first and second peak position and height with Poisson-bootstrap SDs, trough, pair counts
  block_fft       hex-correct radially averaged power spectrum: peak wavelength, prominence, width, relative FWHM
  settling        per channel, from the recorded frames: settle_t_crit = first t from which the L1 change per cell between
                  consecutive recorded frames (t = 10 apart) stays below 1e-3 to the end; settle_t_count = first t from which the
                  real-domain count stays constant; settle_t_cells = first t from which the persistent cells-on changes by less
                  than 1 % per frame; NaN where the condition never holds.
  exponent        per b_i, point and channel: OLS slope of ln(spectral peak wavelength) on ln(real-domain count) across every run
                  of the five conditions (runs with >= 2 real domains and a spectral peak), with its standard error, n and r^2.
Under nucleation real domains are components of at least 6 cells and transients are reported separately, on the engine's
even-r adjacency, fixed threshold 0.3 r_a with persistence over two frames.

Usage:
  python analysis/fig1c_batch.py time                        # one run (bi12, nucl) with the full pipeline into results/fig1c/timing
  python analysis/fig1c_batch.py run [workers] [set ...] [--conds c1,c2]
  python analysis/fig1c_batch.py collect                     # parts -> analysis/summaries/fig1c_*.csv (+ the exponent table)
  python analysis/fig1c_batch.py remeasure [workers]         # recompute all measurements from the stored frames, then collect
  python analysis/fig1c_batch.py figures [set ...]           # hex-correct grids per set and channel
  python analysis/fig1c_batch.py exponent                    # the exponent table and figure from the collected CSVs
Outputs: results/fig1c/<set>_<cond>_rep<r>.npz (frames of both channels, snapshots, final state in float64 for continuation);
  results/fig1c/parts/*; analysis/summaries/fig1c_{runs,frames,gr,gr_curves,fft,fft_curves,tracks,tracking,block_domains,
  block_spacing,block_order,block_gr,block_fft,settling,exponent}.csv; analysis/summaries/fig1c_<set>_{a,i}.png, fig1c_exponent.png.
"""
from __future__ import annotations

import contextlib, csv, io, json, resource, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay, QhullError, cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "2D_simulations"))
from analysis.hexgeom import cell_centres, window  # noqa: E402
from analysis.ladder_batch import (FRAME_EVERY, GR_T, MIN_REAL, SNAP_T, FIG_T, components, area_stats, gr_features,  # noqa: E402
                                   hex_power_spectrum, measure_frame, measure_run, write_parts)
from analysis.pointstats import border_cv, nn_spacing  # noqa: E402
from analysis.regime_map import activated_states  # noqa: E402
from analysis.seed_crowding import DT, SAVE_EVERY, central_seed, circuit, frame_time  # noqa: E402

N, T_END = 200, 600.0
BA, GAMMA, NA, NI, DI, LEVEL = 5.0, 0.5, 3, 3, 10.0, 10.0
LAM = float(np.sqrt(1.5 * DI / GAMMA))
SETS = {"bi0": (0.0, "JA, no inhibitor"), "bi2.5": (2.5, "inside the Turing band"), "bi3.2": (3.2, "Turing band, near its edge"),
        "bi5": (5.0, "replicating; Fig 1C regular spots"), "bi9": (9.0, "divides late"), "bi12": (12.0, "arrests; Fig 1C irregular spots"),
        "bi20": (20.0, "arrests, small domains")}
ORDER = ["bi5", "bi12", "bi0", "bi2.5", "bi9", "bi20", "bi3.2"]        # b_i 3.2 last, then 20: the first to be dropped if time runs out
CONDS = {
    "seed1": dict(init_mode="activator_random_spikes", spike_value=LEVEL, n_points=1, set_peak_height=LEVEL, nucleation_rate=0.0),
    "seed36": dict(init_mode="activator_random_spikes", spike_value=LEVEL, n_points=36, set_peak_height=LEVEL, nucleation_rate=0.0),
    "random": dict(init_mode="random_uniform_over0", spike_value=1.0, n_points=0, set_peak_height=LEVEL, nucleation_rate=0.0),
    "nucl": dict(init_mode="all_off", spike_value=2.0, n_points=0, set_peak_height=LEVEL, nucleation_rate=0.01),
    "smallpert": dict(init_mode="random_tight", spike_value=2.0, n_points=0, set_peak_height=LEVEL, nucleation_rate=0.0),
}
REPS = {"seed1": [0], "seed36": [1, 2, 3], "random": [1, 2, 3], "nucl": [1, 2, 3], "smallpert": [1, 2, 3]}
COV_POINTS = ((0.10, "cov10"), (0.20, "cov20"))
NULL_REAL, NULL_SEED = 200, 20261010
SETTLE_CRIT = 1e-3                                  # L1 change per cell between recorded frames (t = 10 apart)
RES = ROOT / "results" / "fig1c"; PARTS = RES / "parts"; SUMM = ROOT / "analysis" / "summaries"
FAMILIES = ("runs", "frames", "gr", "gr_curves", "fft", "fft_curves", "tracks", "tracking",
            "block_domains", "block_spacing", "block_order", "block_gr", "block_fft", "settling")


def set_params(name):
    bi, anchor = SETS[name]
    return bi, anchor, DI


def closed_form_state(bi):
    """Upper-branch activated state (a0, i0, alpha, G, tau0) of the closed form at this b_i, or None."""
    roots = [r for r in activated_states(BA, bi, GAMMA, NA, NI) if r[3] > r[2] - 1]
    if not roots:
        return None
    return max(roots, key=lambda r: r[0])


def req(area):
    return np.sqrt(np.asarray(area, dtype=float) * np.sqrt(3) / 2 / np.pi)


# ---------------------------------------------------------------------------------------------- matched nulls
def place_hard_core(n, radii, win, rng, max_tries=500):
    """Random sequential placement of n centres in the window with centre distance >= radii[i] + radii[j]; None when jammed."""
    xmin, xmax, ymin, ymax = win
    pts = np.empty((n, 2)); rad = np.broadcast_to(np.asarray(radii, dtype=float), (n,))
    for i in range(n):
        for _ in range(max_tries):
            c = np.array([rng.uniform(xmin, xmax), rng.uniform(ymin, ymax)])
            if i == 0 or np.all(np.hypot(pts[:i, 0] - c[0], pts[:i, 1] - c[1]) >= rad[:i] + rad[i]):
                pts[i] = c; break
        else:
            return None
    return pts


def order_params(xy, win, margin):
    """psi6 (mean over interior points of |mean exp(6 i theta)| over Delaunay neighbours) and the fraction of interior points with
    exactly six neighbours; interior = farther than margin from the window edge. (nan, nan, 0) when fewer than 4 points."""
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 4:
        return np.nan, np.nan, 0
    try:
        tri = Delaunay(xy)
    except (QhullError, ValueError):
        return np.nan, np.nan, 0
    indptr, indices = tri.vertex_neighbor_vertices
    xmin, xmax, ymin, ymax = win
    interior = (xy[:, 0] > xmin + margin) & (xy[:, 0] < xmax - margin) & (xy[:, 1] > ymin + margin) & (xy[:, 1] < ymax - margin)
    psi, six = [], []
    for i in np.nonzero(interior)[0]:
        nb = indices[indptr[i]:indptr[i + 1]]
        if nb.size == 0:
            continue
        th = np.arctan2(xy[nb, 1] - xy[i, 1], xy[nb, 0] - xy[i, 0])
        psi.append(abs(np.exp(6j * th).mean())); six.append(nb.size == 6)
    if not psi:
        return np.nan, np.nan, 0
    return float(np.mean(psi)), float(np.mean(six)), len(psi)


def point_stats(xy, win, margin):
    d = nn_spacing(xy)
    cv = float(d.std() / d.mean()) if d.size > 1 and d.mean() > 0 else np.nan
    cvb, nb = border_cv(xy, win) if len(xy) >= 3 else (np.nan, 0)
    psi6, frac6, n_int = order_params(xy, win, margin)
    return dict(nn_mean=float(d.mean()) if d.size else np.nan, nn_sd=float(d.std()) if d.size > 1 else np.nan, nn_cv=cv, nn_cv_border=cvb,
                nn_n_border=nb, psi6=psi6, frac6=frac6, n_interior=n_int)


def matched_nulls(xy, radii, hole, win, margin, seed):
    """The three nulls at the measured count: (a) CSR, (b) hard-core discs with the measured per-domain radii, (c) hard-core at the
    hole radius. Returns per null: n placed, mean/SD of nn_cv, nn_cv_border, psi6, frac6 over the realisations."""
    n = len(xy); out = {}
    xmin, xmax, ymin, ymax = win
    for tag, kind in (("csr", "a"), ("disc", "b"), ("hole", "c")):
        rng = np.random.default_rng([NULL_SEED, seed, {"csr": 1, "disc": 2, "hole": 3}[tag]])
        vals = {k: [] for k in ("nn_cv", "nn_cv_border", "psi6", "frac6")}; placed = 0
        if n >= 3 and (kind != "c" or np.isfinite(hole)):
            for _ in range(NULL_REAL):
                if kind == "a":
                    pts = np.column_stack([rng.uniform(xmin, xmax, n), rng.uniform(ymin, ymax, n)])
                elif kind == "b":
                    pts = place_hard_core(n, rng.permutation(radii), win, rng)
                else:
                    pts = place_hard_core(n, 0.5 * hole, win, rng)
                if pts is None:
                    continue
                placed += 1; s = point_stats(pts, win, margin)
                for k in vals:
                    vals[k].append(s[k])
        out[f"null_{tag}_n"] = placed
        for k, v in vals.items():
            v = np.asarray(v, dtype=float); fin = v[np.isfinite(v)]
            out[f"null_{tag}_{k}_mean"] = float(fin.mean()) if fin.size >= 20 else np.nan
            out[f"null_{tag}_{k}_sd"] = float(fin.std(ddof=1)) if fin.size >= 20 else np.nan
    return out


def z_cols(meas, nulls):
    out = {}
    for tag in ("csr", "disc", "hole"):
        for k in ("nn_cv", "nn_cv_border", "psi6", "frac6"):
            m, s = nulls[f"null_{tag}_{k}_mean"], nulls[f"null_{tag}_{k}_sd"]
            out[f"z_{tag}_{k}"] = float((meas[k] - m) / s) if np.isfinite(m) and np.isfinite(s) and s > 0 and np.isfinite(meas[k]) else np.nan
    return out


# ---------------------------------------------------------------------------------------------- the block at three time points
def crossing_frames(frame_rows, ch):
    rows = sorted([r for r in frame_rows if r["channel"] == ch and r["frame"] >= FRAME_EVERY + 1], key=lambda r: r["frame"])   # recorded frames t >= 10
    out = {}
    for level, name in COV_POINTS:
        hit = [r for r in rows if r["coverage"] >= level]
        out[name] = (hit[0]["frame"], hit[0]["t"]) if hit else (None, None)
    return out


def block_measures(fields, rec, thr_fixed, win, lattice_xy, meta, rep, frame_rows):
    dom_rows, sp_rows, ord_rows, gr_rows, fft_rows = [], [], [], [], []
    cross_a = crossing_frames(frame_rows, "a")
    points = [(name, *cross_a[name]) for _, name in COV_POINTS] + [("t600", rec[-1], round(frame_time(rec[-1]), 2))]
    for ch, H in fields.items():
        cross_ch = crossing_frames(frame_rows, ch)
        i_floor = thr_fixed["i"] / 30.0
        def level(k):
            return thr_fixed["a"] if ch == "a" else (0.5 * float(H[k].max()) if H[k].max() >= i_floor else np.inf)
        for name, k, t in points:
            base = dict(**meta, channel=ch, point=name, frame=k if k is not None else np.nan, t=t if t is not None else np.nan,
                        reached=int(k is not None), own_channel_frame=cross_ch.get(name, (k, t))[0] if name != "t600" else k,
                        own_channel_t=cross_ch.get(name, (k, t))[1] if name != "t600" else t, lambda_i=LAM)
            if k is None:
                dom_rows.append(dict(**base, note=f"activator coverage never reached {name[3:]} %"))
                continue
            kp = k - FRAME_EVERY if k - FRAME_EVERY >= FRAME_EVERY + 1 else None
            m, comp, xy = measure_frame(H[k], H[kp] > level(kp) if kp is not None else None, level(k), win, None if False else _NoCsr())
            real = comp["area"] >= MIN_REAL
            areas_real = comp["area"][real]; radii = req(areas_real)
            dom_rows.append(dict(**base, threshold=level(k), persistence=int(kp is not None), n_real=m["blobs_real"], n_transient=m["n_transient"],
                                 n_all=m["blobs_all"], coverage=m["coverage"], coverage_raw=m["coverage_raw"], cells_on=m["cells_on"], largest=m["largest"],
                                 req_mean=float(radii.mean()) if radii.size else np.nan, req_median=float(np.median(radii)) if radii.size else np.nan,
                                 req_mean_over_lambda=float(radii.mean() / LAM) if radii.size else np.nan,
                                 req_median_over_lambda=float(np.median(radii) / LAM) if radii.size else np.nan,
                                 **{kk: m[kk] for kk in m if kk.startswith("area_")}, ecc_mean=m["ecc_mean"], ecc_sd=m["ecc_sd"], areas=m["areas"]))
            # g(r) of the field
            med = m["area_real_median"]; null = "lattice" if (np.isfinite(med) and med < 10) else "uniform"
            feat, _curve = gr_features(xy, win, lattice_xy if null == "lattice" else None, seed=rep * 1000 + int(round(t)) + 7)
            gr_rows.append(dict(**base, n_points=len(xy), null=null, **feat))
            # spacing and order with the matched nulls
            margin = float(nn_spacing(xy).mean()) if len(xy) >= 2 else 0.0
            meas = point_stats(xy, win, margin) if len(xy) else dict(nn_mean=np.nan, nn_sd=np.nan, nn_cv=np.nan, nn_cv_border=np.nan, nn_n_border=0,
                                                                   psi6=np.nan, frac6=np.nan, n_interior=0)
            nulls = matched_nulls(xy, radii, feat["hole"], win, margin, seed=rep * 100000 + int(round(t)) * 10 + (1 if ch == "a" else 2))
            z = z_cols(meas, nulls)
            sp_rows.append(dict(**base, n_points=len(xy), margin=margin, hole_used=feat["hole"], **{kk: meas[kk] for kk in ("nn_mean", "nn_sd", "nn_cv", "nn_cv_border", "nn_n_border")},
                                **{kk: v for kk, v in nulls.items() if "psi6" not in kk and "frac6" not in kk},
                                **{kk: v for kk, v in z.items() if "psi6" not in kk and "frac6" not in kk}))
            ord_rows.append(dict(**base, n_points=len(xy), margin=margin, psi6=meas["psi6"], frac6=meas["frac6"], n_interior=meas["n_interior"],
                                 **{kk: v for kk, v in nulls.items() if "psi6" in kk or "frac6" in kk or kk.endswith("_n")},
                                 **{kk: v for kk, v in z.items() if "psi6" in kk or "frac6" in kk}))
            ff, _fc = hex_power_spectrum(H[k])
            fft_rows.append(dict(**base, n_real=m["blobs_real"], field_mean=float(H[k].mean()), field_var=float(H[k].var()), **ff))
    return dom_rows, sp_rows, ord_rows, gr_rows, fft_rows


class _NoCsr:
    """measure_frame wants a CSR-reference callable; the block carries its own nulls, so this returns blanks."""
    def __call__(self, n):
        return dict(csr_cv=np.nan, csr_cv_border=np.nan, csr_cv_sd=np.nan, csr_cv_border_sd=np.nan, csr_n_real=0)


def settling_rows(fields, rec, frame_rows, meta):
    rows = []
    n2 = fields["a"][rec[0]].size
    for ch, H in fields.items():
        ks = [k for k in rec if k >= FRAME_EVERY + 1]           # t = 10, 20, ..., 600
        ts = [round(frame_time(k), 2) for k in ks]
        l1 = np.array([np.abs(H[k] - H[kp]).sum() / n2 for kp, k in zip(ks[:-1], ks[1:])])       # change over each t = 10 interval, at ts[1:]
        fr = {r["frame"]: r for r in frame_rows if r["channel"] == ch}
        counts = np.array([fr[k]["blobs_real"] for k in ks]); cells = np.array([fr[k]["cells_on"] for k in ks], dtype=float)
        def first_from(ok, times):          # first time from which ok holds to the end
            if not ok.size or not ok[-1]:
                return np.nan
            j = len(ok) - 1
            while j > 0 and ok[j - 1]:
                j -= 1
            return float(times[j])
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.abs(np.diff(cells)) / np.maximum(cells[1:], 1.0)
        rows.append(dict(**meta, channel=ch, settle_t_crit=first_from(l1 < SETTLE_CRIT, ts[1:]), settle_crit_level=SETTLE_CRIT,
                         settle_t_count=first_from(np.diff(counts) == 0, ts[1:]), settle_t_cells=first_from(rel < 0.01, ts[1:]),
                         l1_final=float(l1[-1]) if l1.size else np.nan, l1_max_after_t300=float(l1[np.array(ts[1:]) >= 300].max()) if (np.array(ts[1:]) >= 300).any() else np.nan,
                         count_final=int(counts[-1]), count_max=int(counts.max()), t_count_max=float(ts[int(np.argmax(counts))]),
                         cells_final=float(cells[-1])))
    return rows


# ---------------------------------------------------------------------------------------------- one run
def run_one(job):
    setname, cond, rep, out_dir = job
    try:
        return _run_one(job)
    except Exception as exc:
        import traceback
        (out_dir / "parts").mkdir(parents=True, exist_ok=True)
        (out_dir / "parts" / f"{setname}_{cond}_rep{rep}_error.txt").write_text(traceback.format_exc())
        print(f"RUN FAILED {setname}_{cond}_rep{rep}: {exc!r}", flush=True)
        return dict(set=setname, condition=cond, replicate=rep, error=repr(exc))


def _run_one(job):
    setname, cond, rep, out_dir = job
    import simulation_2D
    from simulation_2D import run_coupled_hex
    bi, anchor, D = set_params(setname)
    p = circuit(BA, bi, GAMMA, NA, NI, D); c = CONDS[cond]
    thr_fixed = {"a": 0.3 * p["act_prod_rate"], "i": 0.3 * p["inh_prod_rate"] / p["inh_decay_rate"]}
    steps = int(round(T_END / DT)) + 1
    out_dir.mkdir(parents=True, exist_ok=True); (out_dir / "parts").mkdir(exist_ok=True)
    tag = f"{setname}_{cond}_rep{rep}"
    meta = dict(set=setname, b_i=bi, D_i=D, n_a=NA, n_i=NI, anchor=anchor, condition=cond, replicate=rep)
    spike_seed = ""
    np.random.seed(rep)
    orig_rng = np.random.default_rng; orig_solver = simulation_2D.fast_stable_steady_state
    cf = closed_form_state(bi); cf_cols = dict(cf_a_ss=cf[0] if cf else np.nan, cf_i_ss=cf[1] if cf else np.nan, cf_alpha=cf[2] if cf else np.nan,
                                                cf_G=cf[3] if cf else np.nan, cf_tau0=cf[4] if cf else np.nan, cf_state=int(cf is not None))
    if c["init_mode"] == "activator_random_spikes":
        spike_seed = central_seed(N) if c["n_points"] == 1 else 1000 * rep + 7
        np.random.default_rng = lambda seed=None, _s=spike_seed: orig_rng(_s)
    if cond == "smallpert":          # the engine's solver is replaced by the closed form for this condition only
        state = (float(cf[0]), max(float(cf[1]), 1e-12), 0.0) if cf else (0.0, 0.0, 0.0)
        simulation_2D.fast_stable_steady_state = lambda *a, **k: state
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            A, I, step, a_ss, i_ss = run_coupled_hex(N, N, steps, DT, 1.0, p, 0.0, 10 ** 9, init_mode=c["init_mode"], activator_type="juxtacrine",
                                                     n_points=c["n_points"], spike_value=c["spike_value"], save_every=SAVE_EVERY,
                                                     nucleation_rate=c["nucleation_rate"], noise_amplitude=0.0, set_peak_height=c["set_peak_height"])
    finally:
        np.random.default_rng = orig_rng; simulation_2D.fast_stable_steady_state = orig_solver
    t_engine = time.perf_counter() - t0
    assert len(A) == steps // SAVE_EVERY + 2 and abs(frame_time(len(A) - 1) - (T_END + DT)) < 1e-9
    A = [np.asarray(f, dtype=np.float64) for f in A]; I = [np.asarray(f, dtype=np.float64) for f in I]
    n2 = N * N
    crit_a = np.array([np.abs(A[k] - A[k - 1]).sum() / n2 for k in range(1, len(A))])
    crit_i = np.array([np.abs(I[k] - I[k - 1]).sum() / n2 for k in range(1, len(I))])
    crit = 0.5 * (crit_a + crit_i); crit_t = np.array([frame_time(k) for k in range(1, len(A))])
    win = window(N, N); lattice_xy = np.column_stack([c_.ravel() for c_ in cell_centres(N, N)])
    fields = {"a": A, "i": I}
    rec = [0] + [k for k in range(1, len(A)) if (k - 1) % FRAME_EVERY == 0]
    snaps = {f"{ch}_t{round(frame_time(k))}": H[k].astype(np.float32) for ch, H in fields.items() for k in rec if k >= 1 and round(frame_time(k)) in SNAP_T}
    np.savez_compressed(out_dir / f"{tag}.npz", a_final=A[-1], i_final=I[-1], t_final=frame_time(len(A) - 1), steps=steps, step_returned=step,
                        crit=crit, crit_a=crit_a, crit_i=crit_i, crit_t=crit_t, a_ss=a_ss, i_ss=i_ss, thr_a=thr_fixed["a"], thr_i_fixed=thr_fixed["i"],
                        params=json.dumps(p), condition=json.dumps(c), spike_seed=str(spike_seed), np_seed=rep, grid=N, closed_form=json.dumps(cf_cols),
                        frame_t=np.array([frame_time(k) for k in rec]), frame_k=np.array(rec),
                        a_frames=np.stack([A[k] for k in rec]).astype(np.float32), i_frames=np.stack([I[k] for k in rec]).astype(np.float32), **snaps)
    t_saved = time.perf_counter() - t0
    try:
        rows_out = measure_all(fields, rec, crit, crit_a, crit_i, thr_fixed, win, lattice_xy, meta, rep)
    except Exception as exc:
        import traceback
        (out_dir / "parts" / f"{tag}_error.txt").write_text(traceback.format_exc())
        print(f"MEASUREMENT FAILED {tag}: {exc!r} (engine output saved)", flush=True)
        return dict(**meta, error=repr(exc), wall_engine_s=round(t_engine, 1))
    t_total = time.perf_counter() - t0
    run_row = dict(**meta, init_mode=c["init_mode"], n_points=c["n_points"], spike_value=c["spike_value"], set_peak_height=c["set_peak_height"],
                   nucleation_rate=c["nucleation_rate"], spike_seed=spike_seed, np_seed=rep, a_ss=float(a_ss), i_ss=float(i_ss), **protocol_cols(a_ss, i_ss, c, cf),
                   **cf_cols, thr_a=thr_fixed["a"], thr_i_fixed=thr_fixed["i"], grid=N, steps=steps, t_final=round(frame_time(len(A) - 1), 2),
                   wall_engine_s=round(t_engine, 1), wall_save_s=round(t_saved - t_engine, 1), wall_total_s=round(t_total, 1),
                   crit_final=float(crit[-1]), **last_cols(rows_out[0], rec))
    write_all(out_dir, tag, rows_out, run_row)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(f"done {tag}: engine {t_engine:.0f} s, total {t_total:.0f} s, peak RSS {rss:.2f} GB, a: {run_row['a_blobs_real']} real / {run_row['a_n_transient']} transient, "
          f"cov {run_row['a_coverage']:.3f}; i: {run_row['i_blobs_real']} real, cov {run_row['i_coverage']:.3f}", flush=True)
    return run_row


def protocol_cols(a_ss, i_ss, c, cf):
    fallback = abs(float(a_ss) - float(c["spike_value"])) < 1e-12 and abs(float(i_ss) - float(c["spike_value"])) < 1e-12
    src = "fallback (= spike_value)" if fallback else ("closed form" if c["init_mode"] == "random_tight" else "solver")
    return dict(a_ss_source=src, nucl_kick=2 * float(a_ss) if c["nucleation_rate"] > 0 else 0.0,
                nucl_gate=float(a_ss) if c["nucleation_rate"] > 0 else np.nan,
                smallpert_note=("" if c["init_mode"] != "random_tight" else
                                ("both fields uniform within +-5 % of the closed-form activated state" if cf and not fallback else
                                 "no activated state at this b_i: engine fallback, both fields uniform within +-5 % of spike_value")))


def last_cols(frame_rows, rec):
    last = {r["channel"]: r for r in frame_rows if r["frame"] == rec[-1]}
    return dict(thr_i_halfmax_final=last["i"]["threshold"],
                **{f"{ch}_{k}": last[ch][k] for ch in ("a", "i") for k in ("cells_on", "coverage", "blobs_all", "blobs_real", "n_transient", "largest",
                                                                          "area_real_median", "area_real_cv", "nn_cv_border", "csr_cv_border", "ecc_mean", "total")})


def measure_all(fields, rec, crit, crit_a, crit_i, thr_fixed, win, lattice_xy, meta, rep):
    base = measure_run(fields, rec, crit, crit_a, crit_i, thr_fixed, win, lattice_xy, meta, rep)
    frame_rows = base[0]
    block = block_measures(fields, rec, thr_fixed, win, lattice_xy, meta, rep, frame_rows)
    settle = settling_rows(fields, rec, frame_rows, meta)
    return (*base, *block, settle)


def write_all(out_dir, tag, rows_out, run_row):
    frame_rows, gr_rows, gr_curves, fft_rows, fft_curves, track_rows_all, track_summ, dom, sp, od, bgr, bfft, settle = rows_out
    write_parts(out_dir, tag, frame_rows, gr_rows, gr_curves, fft_rows, fft_curves, track_rows_all, track_summ, run_row)
    for name, rows in (("block_domains", dom), ("block_spacing", sp), ("block_order", od), ("block_gr", bgr), ("block_fft", bfft), ("settling", settle)):
        if rows:
            cols = list(dict.fromkeys(k for r in rows for k in r))
            with open(out_dir / "parts" / f"{tag}_{name}.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n"); w.writeheader(); w.writerows(rows)


def remeasure(path):
    try:
        _remeasure(path)
    except Exception as exc:
        import traceback
        (path.parent / "parts" / f"{path.stem}_remeasure_error.txt").write_text(traceback.format_exc())
        print(f"REMEASURE FAILED {path.name}: {exc!r}", flush=True)


def _remeasure(path):
    d = np.load(path); tag = path.stem; out_dir = path.parent
    setname, cond, rep = tag.rsplit("_rep", 1)[0].rsplit("_", 1)[0], tag.rsplit("_rep", 1)[0].rsplit("_", 1)[1], int(tag.rsplit("_rep", 1)[1])
    bi, anchor, D = set_params(setname); meta = dict(set=setname, b_i=bi, D_i=D, n_a=NA, n_i=NI, anchor=anchor, condition=cond, replicate=rep)
    rec = [int(k) for k in d["frame_k"]]
    fields = {"a": dict(zip(rec, d["a_frames"].astype(np.float64))), "i": dict(zip(rec, d["i_frames"].astype(np.float64)))}
    n = int(d["grid"]); win = window(n, n); lattice_xy = np.column_stack([c_.ravel() for c_ in cell_centres(n, n)])
    thr_fixed = {"a": float(d["thr_a"]), "i": float(d["thr_i_fixed"])}
    out = measure_all(fields, rec, d["crit"], d["crit_a"], d["crit_i"], thr_fixed, win, lattice_xy, meta, rep)
    runs_part = out_dir / "parts" / f"{tag}_runs.csv"
    run_row = dict(next(csv.DictReader(open(runs_part)))) if runs_part.exists() else dict(**meta)
    c = json.loads(str(d["condition"])); cf = closed_form_state(bi)
    run_row.update(**protocol_cols(float(d["a_ss"]), float(d["i_ss"]), c, cf), **last_cols(out[0], rec))
    write_all(out_dir, tag, out, run_row)
    print(f"remeasured {tag}", flush=True)


# ---------------------------------------------------------------------------------------------- batch, collect, exponent, figures
def parse_sets(args):
    bad = [a for a in args if a not in SETS]
    if bad:
        raise SystemExit(f"unknown set(s) {bad}; known: {list(SETS)}")
    return list(args) or ORDER


def jobs_for(sets, conds=None):
    return [(s, cond, rep, RES) for s in sets for cond in (conds or CONDS) for rep in REPS[cond]]


def collect(out_dir=RES, prefix="fig1c"):
    parts = out_dir / "parts"
    for name in FAMILIES:
        rows, cols = [], []
        for f in sorted(parts.glob(f"*_{name}.csv")):
            with open(f) as fh:
                rd = csv.DictReader(fh)
                rows += list(rd)
                for ccol in rd.fieldnames or []:
                    if ccol not in cols:
                        cols.append(ccol)
        if rows:
            with open(SUMM / f"{prefix}_{name}.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n"); w.writeheader(); w.writerows(rows)
            print(f"{prefix}_{name}.csv: {len(rows)} rows")
    try:
        exponent()
    except Exception as exc:
        print(f"exponent table not written: {exc!r}")


def exponent(prefix="fig1c"):
    """OLS of ln(spectral peak wavelength) on ln(real-domain count) across the runs of each b_i, per point and channel."""
    import pandas as pd
    fft = pd.read_csv(SUMM / f"{prefix}_block_fft.csv")
    rows = []
    for (s, bi, point, ch), g in fft.groupby(["set", "b_i", "point", "channel"]):
        g = g[(g.n_real >= 2) & np.isfinite(g.fft_peak_wavelength)]
        x, y = np.log(g.n_real.values.astype(float)), np.log(g.fft_peak_wavelength.values.astype(float))
        row = dict(set=s, b_i=bi, point=point, channel=ch, n_runs=int(len(g)), conditions=" ".join(sorted(set(g.condition))),
                   count_min=float(g.n_real.min()) if len(g) else np.nan, count_max=float(g.n_real.max()) if len(g) else np.nan,
                   wavelength_min=float(g.fft_peak_wavelength.min()) if len(g) else np.nan, wavelength_max=float(g.fft_peak_wavelength.max()) if len(g) else np.nan,
                   slope=np.nan, slope_se=np.nan, intercept=np.nan, r2=np.nan)
        if len(g) >= 3 and np.ptp(x) > 0:
            X = np.column_stack([np.ones_like(x), x]); beta, res, *_ = np.linalg.lstsq(X, y, rcond=None)
            yhat = X @ beta; rss = float(((y - yhat) ** 2).sum()); tss = float(((y - y.mean()) ** 2).sum())
            dof = len(g) - 2; s2 = rss / dof if dof > 0 else np.nan
            cov = s2 * np.linalg.inv(X.T @ X) if dof > 0 else np.full((2, 2), np.nan)
            row.update(slope=float(beta[1]), slope_se=float(np.sqrt(cov[1, 1])), intercept=float(beta[0]), r2=float(1 - rss / tss) if tss > 0 else np.nan)
        rows.append(row)
    out = pd.DataFrame(rows).sort_values(["point", "channel", "b_i"])
    out.to_csv(SUMM / f"{prefix}_exponent.csv", index=False, lineterminator="\n")
    exponent_figure(fft, out, prefix)
    print(f"{prefix}_exponent.csv: {len(out)} rows")
    return out


def exponent_figure(fft, out, prefix):
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {"seed1": "#2a78d6", "seed36": "#eb6834", "random": "#eda100", "nucl": "#e87ba4", "smallpert": "#1baf7a"}
    bis = sorted(SETS[s][0] for s in SETS if s in set(fft.set))
    for ch in ("a", "i"):
        fig, axes = plt.subplots(3, len(bis), figsize=(2.6 * len(bis), 8.2), squeeze=False); fig.patch.set_facecolor("#fcfcfb")
        for r, point in enumerate(("cov10", "cov20", "t600")):
            for cidx, bi in enumerate(bis):
                ax = axes[r, cidx]; ax.set_facecolor("#fcfcfb"); ax.grid(True, color="#e6e5e1", linewidth=0.6); ax.set_axisbelow(True)
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
                g = fft[(fft.b_i == bi) & (fft.point == point) & (fft.channel == ch) & (fft.n_real >= 2) & np.isfinite(fft.fft_peak_wavelength)]
                for cond, gg in g.groupby("condition"):
                    ax.plot(gg.n_real, gg.fft_peak_wavelength, "o", color=col.get(cond, "#52514e"), markersize=4.5, markerfacecolor="#fcfcfb", markeredgewidth=1.4, label=cond)
                e = out[(out.b_i == bi) & (out.point == point) & (out.channel == ch)]
                if len(e) and np.isfinite(e.slope.iloc[0]):
                    xs = np.linspace(np.log(e.count_min.iloc[0]), np.log(e.count_max.iloc[0]), 20)
                    ax.plot(np.exp(xs), np.exp(e.intercept.iloc[0] + e.slope.iloc[0] * xs), "-", color="#52514e", linewidth=1.2)
                    ax.text(0.03, 0.04, f"slope {e.slope.iloc[0]:.2f} ± {e.slope_se.iloc[0]:.2f}, n {int(e.n_runs.iloc[0])}", transform=ax.transAxes, fontsize=7.5, color="#52514e")
                ax.set_xscale("log"); ax.set_yscale("log")
                ax.set_title(f"b_i {bi:g}, {point}", fontsize=9, loc="left")
                if r == 2:
                    ax.set_xlabel("real domains", fontsize=8, color="#52514e")
                if cidx == 0:
                    ax.set_ylabel("spectral peak wavelength (cells)", fontsize=8, color="#52514e")
                if r == 0 and cidx == 0:
                    ax.legend(fontsize=7, frameon=False)
        fig.suptitle(f"Fig 1C line (b_a 5, n 3/3, D_i 10, gamma 0.5): spectral peak wavelength against domain count, {'activator' if ch == 'a' else 'inhibitor'}; "
                     f"line = OLS fit of ln(wavelength) on ln(count)", fontsize=9.5)
        fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(SUMM / f"{prefix}_exponent_{ch}.png", dpi=150, facecolor=fig.get_facecolor()); plt.close(fig)


def figures(sets):
    from analysis.plot_hex import save_panels
    rows_cond = [("seed1", 0), ("seed36", 1), ("random", 1), ("nucl", 1), ("smallpert", 1)]
    for s in sets:
        bi, anchor, D = set_params(s)
        for ch, label in (("a", "activator"), ("i", "inhibitor")):
            fields, titles = [], []
            for cond, rep in rows_cond:
                f = RES / f"{s}_{cond}_rep{rep}.npz"
                if not f.exists():
                    continue
                d = np.load(f)
                for t in FIG_T:
                    fields.append(d[f"{ch}_t{t}"].astype(float)); titles.append(f"{cond} rep{rep}, t={t}")
            if fields:
                save_panels(fields, titles, SUMM / f"fig1c_{s}_{ch}.png", ncols=len(FIG_T),
                            suptitle=f"{s} ({anchor}): b_a=5, b_i={bi:g}, n=3/3, D_i={D:g}, gamma=0.5; 200x200; {label}, hex-correct, colour 0..max")
                print(f"fig1c_{s}_{ch}.png")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "time":
        out = RES / "timing"
        t0 = time.perf_counter()
        r = run_one(("bi12", "nucl", 1, out))
        print(f"wall {time.perf_counter() - t0:.0f} s for one 200x200 nucleation run with the full pipeline (engine {r.get('wall_engine_s')} s)")
    elif cmd == "run":
        args = sys.argv[2:]
        conds = None
        if "--conds" in args:
            k = args.index("--conds"); conds = args[k + 1].split(","); args = args[:k] + args[k + 2:]
            bad = [c for c in conds if c not in CONDS]
            if bad:
                raise SystemExit(f"unknown condition(s) {bad}; known: {list(CONDS)}")
        n_workers = int(args[0]) if args and args[0].isdigit() else 4
        sets = parse_sets(args[1:] if args and args[0].isdigit() else args)
        jobs = [j for j in jobs_for(sets, conds) if not ((RES / f"{j[0]}_{j[1]}_rep{j[2]}.npz").exists() and (PARTS / f"{j[0]}_{j[1]}_rep{j[2]}_runs.csv").exists())]
        print(f"{len(jobs)} runs on {n_workers} workers: sets {sets} (finished runs are skipped)", flush=True)
        t0 = time.perf_counter()
        try:
            with Pool(n_workers) as pool:
                for k, r in enumerate(pool.imap_unordered(run_one, jobs, chunksize=1), 1):
                    print(f"[{k}/{len(jobs)}] {time.perf_counter() - t0:.0f} s elapsed" + (f"  ERROR in {r['set']}_{r['condition']}_rep{r['replicate']}" if "error" in r else ""), flush=True)
        finally:
            collect()
    elif cmd == "collect":
        collect()
    elif cmd == "remeasure":
        files = sorted(RES.glob("*.npz"))
        try:
            with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as pool:
                pool.map(remeasure, files, chunksize=1)
        finally:
            collect()
    elif cmd == "figures":
        figures(parse_sets(sys.argv[2:]))
    elif cmd == "exponent":
        exponent()
