"""
Overnight ladder batch (9-10 October 2026): a ladder in b_i at the measured Hill pair (b_a 5, n 10/4, gamma 0.5, D_i 20),
six initial conditions per set as in REPORT section 11, 200x200, to t = 600, frames every t = 10, convergence check
disabled. Everything is recorded; nothing is interpreted here.

Sets (b_i): 13.5 (below the Turing sliver), 14.1 (inside the engine-confirmed sliver), 15 (3B_big), 20 (3L), 25 (2L_JAPI),
30 (3B_small), 50 (S8D_lowri). Note: at b_i = 13.5 and 14.1 the engine's solver finds an activated state (a_ss 4.58 / 4.35), so there
the nucleation kick is 2 a_ss = 9.2 / 8.7 gated at a < a_ss, while the five sets without a state use the fallback a_ss = spike_value = 2
(kick 4, gate a < 2); the run table records a_ss_source, nucl_kick and nucl_gate. Nothing is forced to a common kick. Conditions: seed1 (one activator-only cell at level 10 near the centre), seed9/36/144
(activator_random_spikes at level 10, placement seed overridden per replicate), random (random_uniform_over0, U(0, 2)),
nucl (all_off, nucleation 0.01, kick 2 a_ss); replicates 1-3 for the field conditions, 0 for seed1.

Measurements (both channels: a = activator at the fixed threshold 0.3 r_a; i = inhibitor at half its own frame maximum, because the
fixed level 0.3 r_i / gamma lies above the diffusion-flattened inhibitor field in most runs; the inhibitor is also counted at that fixed
level in the i_fixed_* columns), components on the engine's hex adjacency, persistence over the frame and the previous frame (t - 10),
real domains = components of at least 6 cells, transients = smaller ones:
  per frame (t = 0, 0.01, 10, 20, ..., 600): cells on and coverage (raw and persistent), blob counts >= 1 and >= 6 cells, largest,
    total field, per-domain area list, area mean/SD/CV/median/p10/p90 (real and all), nearest-neighbour spacing of real-domain
    centroids (mean, SD, CV, Hanisch border-corrected CV) with the same-window same-n CSR reference for both, eccentricity mean/SD,
    the engine's convergence diagnostic for the 200-step window ending at the frame;
  at t = 50, 100, 200, 400, 600: g(r) features and curves (PooledG; lattice-site null when the median real domain is under 10 cells)
    with Poisson-bootstrap SDs, hex-correct radially averaged power spectrum with peak wavelength, prominence and width;
  across frames: domain tracks (nearest-centroid linking within 3 cells): births, deaths, lifetimes, centroid drift;
  snapshots of both fields at t = 50, 100, 150, 200, 300, 400, 600 (float32) and the final state (float64) for continuation.

Usage:
  python analysis/ladder_batch.py time                      # one representative run (bi25, nucl) with the full pipeline; prints timings
  python analysis/ladder_batch.py run [workers] [set ...] [--conds c1,c2]   # the batch (default: all sets in ORDER, all conditions)
  python analysis/ladder_batch.py collect                   # concatenate results/ladder/parts/* into analysis/summaries/ladder_*.csv
  python analysis/ladder_batch.py remeasure [workers]       # recompute all measurements from the stored frames (after a code change), then collect
  python analysis/ladder_batch.py figures [set ...]         # hex-correct grids per set and channel from the saved snapshots
Outputs: results/ladder/<set>_<cond>_rep<r>.npz; results/ladder/parts/*; analysis/summaries/ladder_{runs,frames,gr,gr_curves,fft,
  fft_curves,tracks,tracking}.csv; analysis/summaries/ladder_<set>_{a,i}.png.
"""
from __future__ import annotations

import contextlib, csv, io, json, os, resource, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "2D_simulations"))
from analysis.hexgeom import cell_centres, label_components, window  # noqa: E402
from analysis.pointstats import PooledG, border_cv, nn_spacing  # noqa: E402
from analysis.seed_crowding import DT, SAVE_EVERY, central_seed, circuit, frame_time  # noqa: E402

N, T_END = 200, 600.0
BA, GAMMA, NA, NI, DI, LEVEL = 5.0, 0.5, 10, 4, 20.0, 10.0
SETS = {"bi13.5": (13.5, "below the Turing sliver"), "bi14.1": (14.1, "inside the Turing sliver"), "bi15": (15.0, "3B_big"),
        "bi20": (20.0, "3L"), "bi25": (25.0, "2L_JAPI"), "bi30": (30.0, "3B_small"), "bi50": (50.0, "S8D_lowri"),
        # D_i sweep at the 2L set (b_i 25): the inhibitor range sqrt(1.5 D_i / gamma) is 3.87, 5.48 and 7.75 cells at D_i 5, 10, 20 (bi25)
        "bi25_D5": (25.0, "2L_JAPI at D_i 5", 5.0), "bi25_D10": (25.0, "2L_JAPI at D_i 10", 10.0)}


def set_params(name):
    """(b_i, anchor, D_i) of a set; D_i defaults to the ladder's 20."""
    v = SETS[name]
    return v[0], v[1], (v[2] if len(v) > 2 else DI)
ORDER = ["bi13.5", "bi14.1", "bi15", "bi25", "bi50", "bi20", "bi30"]      # the middle of the ladder runs last
CONDS = {
    "seed1": dict(init_mode="activator_random_spikes", spike_value=LEVEL, n_points=1, set_peak_height=LEVEL, nucleation_rate=0.0),
    "seed9": dict(init_mode="activator_random_spikes", spike_value=LEVEL, n_points=9, set_peak_height=LEVEL, nucleation_rate=0.0),
    "seed36": dict(init_mode="activator_random_spikes", spike_value=LEVEL, n_points=36, set_peak_height=LEVEL, nucleation_rate=0.0),
    "seed144": dict(init_mode="activator_random_spikes", spike_value=LEVEL, n_points=144, set_peak_height=LEVEL, nucleation_rate=0.0),
    "random": dict(init_mode="random_uniform_over0", spike_value=2.0, n_points=0, set_peak_height=LEVEL, nucleation_rate=0.0),
    "nucl": dict(init_mode="all_off", spike_value=2.0, n_points=0, set_peak_height=LEVEL, nucleation_rate=0.01),
}
REPS = {"seed1": [0], "seed9": [1, 2, 3], "seed36": [1, 2, 3], "seed144": [1, 2, 3], "random": [1, 2, 3], "nucl": [1, 2, 3]}
SNAP_T = (50, 100, 150, 200, 300, 400, 600)      # saved fields
GR_T = (50, 100, 200, 400, 600)                   # g(r) and power spectrum
FIG_T = (50, 150, 300, 600)
MIN_REAL, LINK_R, FRAME_EVERY = 6, 3.0, 5         # FRAME_EVERY saved frames (200 steps each) per t = 10 record
CSR_REAL, GR_NULL, GR_BOOT = 100, 100, 200
RES = ROOT / "results" / "ladder"; PARTS = RES / "parts"; SUMM = ROOT / "analysis" / "summaries"


# ---------------------------------------------------------------------------------------------- geometry helpers
def components(mask):
    """Per-component area, centroid and eccentricity of a boolean mask on the engine's adjacency."""
    if not mask.any():
        return dict(labels=None, area=np.zeros(0), cx=np.zeros(0), cy=np.zeros(0), ecc=np.zeros(0))
    labels, n = label_components(mask)
    x, y = cell_centres(*mask.shape)
    lab, xs, ys = labels[mask], x[mask], y[mask]
    area = np.bincount(lab, minlength=n).astype(float)
    cx = np.bincount(lab, weights=xs, minlength=n) / area; cy = np.bincount(lab, weights=ys, minlength=n) / area
    sxx = np.bincount(lab, weights=xs * xs, minlength=n) / area - cx ** 2
    syy = np.bincount(lab, weights=ys * ys, minlength=n) / area - cy ** 2
    sxy = np.bincount(lab, weights=xs * ys, minlength=n) / area - cx * cy
    half, det = 0.5 * (sxx + syy), sxx * syy - sxy ** 2
    disc = np.sqrt(np.maximum(half ** 2 - det, 0.0)); l1, l2 = half + disc, half - disc
    with np.errstate(divide="ignore", invalid="ignore"):
        ecc = np.where((area >= 3) & (l1 > 0), np.sqrt(np.clip(1.0 - l2 / l1, 0.0, 1.0)), np.nan)   # 0 round, 1 line
    return dict(labels=labels, area=area, cx=cx, cy=cy, ecc=ecc)


def area_stats(areas, prefix):
    a = np.asarray(areas, dtype=float)
    if a.size == 0:
        return {f"{prefix}_n": 0, **{f"{prefix}_{k}": np.nan for k in ("mean", "sd", "cv", "median", "p10", "p90")}}
    sd = float(a.std(ddof=1)) if a.size > 1 else 0.0
    return {f"{prefix}_n": int(a.size), f"{prefix}_mean": float(a.mean()), f"{prefix}_sd": sd, f"{prefix}_cv": sd / float(a.mean()),
            f"{prefix}_median": float(np.median(a)), f"{prefix}_p10": float(np.percentile(a, 10)), f"{prefix}_p90": float(np.percentile(a, 90))}


CSR_SEED = 20261009
_CSR_CACHE = {}


class CsrCache:
    """Same-window, same-n CSR reference (plain and border-corrected NN CV): a pure function of (n, window) seeded by n, so every
    channel, replicate, set, frame and remeasure gets the identical value; n_real = max(CSR_REAL, 10000 // n) realisations (<= 1000);
    the Monte-Carlo SD of each reference is returned with it."""
    def __init__(self, win, seed=None):
        self.win = win

    def __call__(self, n):
        if n < 3:
            return dict(csr_cv=np.nan, csr_cv_border=np.nan, csr_cv_sd=np.nan, csr_cv_border_sd=np.nan, csr_n_real=0)
        key = (n, self.win)
        if key not in _CSR_CACHE:
            rng = np.random.default_rng([CSR_SEED, n]); n_real = min(1000, max(CSR_REAL, 10000 // n))
            xmin, xmax, ymin, ymax = self.win
            plain, border = [], []
            for _ in range(n_real):
                pts = np.column_stack([rng.uniform(xmin, xmax, n), rng.uniform(ymin, ymax, n)])
                d = nn_spacing(pts); plain.append(d.std() / d.mean()); border.append(border_cv(pts, self.win)[0])
            plain, border = np.array(plain), np.array(border, dtype=float)
            _CSR_CACHE[key] = dict(csr_cv=float(plain.mean()), csr_cv_border=float(np.nanmean(border)), csr_cv_sd=float(plain.std()),
                                   csr_cv_border_sd=float(np.nanstd(border)), csr_n_real=n_real)
        return _CSR_CACHE[key]


def measure_frame(F, prev_mask, thr, win, csr):
    raw = F > thr
    mask = raw & prev_mask if prev_mask is not None else raw
    comp = components(mask)
    real = comp["area"] >= MIN_REAL
    areas_all = np.sort(comp["area"])[::-1].astype(int); areas_real = areas_all[areas_all >= MIN_REAL]
    xy = np.column_stack([comp["cx"][real], comp["cy"][real]]) if real.any() else np.empty((0, 2))
    d = nn_spacing(xy)
    cvb, nb = border_cv(xy, win) if len(xy) >= 3 else (np.nan, 0)
    csr_ref = csr(len(xy))
    ecc = comp["ecc"][real]; ecc = ecc[np.isfinite(ecc)]
    return dict(cells_on_raw=int(raw.sum()), coverage_raw=float(raw.mean()), cells_on=int(mask.sum()), coverage=float(mask.mean()),
                blobs_all=int(len(areas_all)), blobs_real=int(len(areas_real)), n_transient=int(len(areas_all) - len(areas_real)),
                largest=int(areas_all[0]) if len(areas_all) else 0, total=float(F.sum()), field_max=float(F.max()),
                **area_stats(areas_real, "area_real"), **area_stats(areas_all, "area_all"),
                nn_mean=float(d.mean()) if d.size else np.nan, nn_sd=float(d.std()) if d.size > 1 else np.nan,     # population SD: nn_cv = nn_sd / nn_mean, as in csr_reference_cv
                nn_cv=float(d.std() / d.mean()) if d.size > 1 and d.mean() > 0 else np.nan, nn_cv_border=cvb, nn_n_border=nb,
                **csr_ref,
                ecc_mean=float(ecc.mean()) if ecc.size else np.nan, ecc_sd=float(ecc.std(ddof=1)) if ecc.size > 1 else np.nan,
                areas=" ".join(map(str, areas_all))), comp, xy


# ---------------------------------------------------------------------------------------------- g(r) and power spectrum
def gr_features(xy, win, lattice_xy, seed):
    pg = PooledG(r_max=40.0, dr=0.5, n_null=GR_NULL, seed=seed, lattice_xy=lattice_xy)
    pg.add(xy, win)
    feat = pg.summary(n_boot=GR_BOOT)
    curve = [dict(r=float(r), g=float(g), h_data=float(hd), h_null=float(hn)) for r, g, hd, hn in zip(pg.r, pg.g(), pg.h_data, pg.h_null)]
    return feat, curve


def hex_power_spectrum(F, k_min_wavelength=80.0):
    """
    Radially averaged power spectrum of a field on the engine's even-r hex lattice, computed at the physical cell positions
    (x = c + 0.5 on even rows, y = r sqrt(3)/2): a row-wise FFT in x combined with the exact phase factors of the row offsets and
    row positions in y, after mean removal and a separable Hann window. Bins of width 2 pi / nx in |q| up to pi.
    The q = 0 mode (the DC residual of the windowed, mean-removed field) is excluded from the bins and from the totals.
    Returns the curve (k, P) and the strongest interior peak with wavelength 2 pi / k, its prominence, its full width at half
    PROMINENCE (scipy peak_widths, rel_height 0.5: the contour at P_peak - prominence/2, equal to the FWHM only when
    fft_peak_prominence_rel = 1; fft_power_in_band is the fraction of the power inside that band) and an explicit FWHM
    (fft_peak_fwhm_k: the contour at P_peak/2, walked outwards from the peak; fft_peak_fwhm_border = 1 when a side reaches the band edge).
    Checked against a direct DFT at the hex positions (agreement to 1e-15) and on a triangular arrangement of spots with
    nearest-neighbour spacing L, whose first reciprocal shell gives peak wavelength L sqrt(3)/2 (17.4 for L = 20), not L.
    """
    ny, nx = F.shape
    wy, wx = np.hanning(ny), np.hanning(nx)
    f = (F - F.mean()) * wy[:, None] * wx[None, :]
    Fx = np.fft.fft(f, axis=1)
    qx = 2 * np.pi * np.fft.fftfreq(nx)
    qy = 2 * np.pi * np.fft.fftfreq(ny, d=np.sqrt(3) / 2)
    yr = np.arange(ny) * np.sqrt(3) / 2; off = 0.5 * (np.arange(ny) % 2 == 0)
    G = np.exp(-1j * np.outer(qy, yr)) @ (Fx * np.exp(-1j * np.outer(off, qx)))
    S = np.abs(G) ** 2 / ((wy ** 2).sum() * (wx ** 2).sum())
    K = np.sqrt(qx[None, :] ** 2 + qy[:, None] ** 2)
    dk = 2 * np.pi / nx; edges = np.arange(0.0, np.pi + dk, dk)
    idx = np.digitize(K.ravel(), edges) - 1; ok = (idx >= 0) & (idx < len(edges) - 1) & (K.ravel() > 0)
    cnt = np.bincount(idx[ok], minlength=len(edges) - 1)
    P = np.bincount(idx[ok], weights=S.ravel()[ok], minlength=len(edges) - 1) / np.maximum(cnt, 1)
    k = 0.5 * (edges[:-1] + edges[1:])
    total = float(S.ravel()[ok].sum())
    feat = dict(fft_peak_k=np.nan, fft_peak_wavelength=np.nan, fft_peak_power=np.nan, fft_peak_prominence=np.nan, fft_peak_prominence_rel=np.nan,
                fft_peak_width_k=np.nan, fft_peak_width_rel=np.nan, fft_peak_fwhm_k=np.nan, fft_peak_fwhm_rel=np.nan, fft_peak_fwhm_border=np.nan,
                fft_total_power=total, fft_power_in_band=np.nan)
    sel = k >= 2 * np.pi / k_min_wavelength
    Ps = P[sel]; ks = k[sel]
    peaks, _ = find_peaks(Ps)
    if peaks.size:
        prom = peak_prominences(Ps, peaks)[0]
        p = peaks[int(np.argmax(Ps[peaks]))]
        pr = float(prom[list(peaks).index(p)])
        w = peak_widths(Ps, [p], rel_height=0.5)
        width_k = float(w[0][0] * dk); lo, hi = ks[int(round(w[2][0]))] if 0 <= int(round(w[2][0])) < ks.size else ks[0], ks[min(int(round(w[3][0])), ks.size - 1)]
        band = (k >= lo) & (k <= hi)
        half = 0.5 * Ps[p]; L = p; R = p                     # explicit FWHM: walk out from the peak to the half-height crossings
        while L > 0 and Ps[L] > half:
            L -= 1
        while R < Ps.size - 1 and Ps[R] > half:
            R += 1
        border = int(L == 0 and Ps[L] > half or R == Ps.size - 1 and Ps[R] > half)
        kl = ks[L] + (half - Ps[L]) / (Ps[L + 1] - Ps[L]) * dk if Ps[L] <= half < Ps[L + 1] else ks[L]
        kr = ks[R] - (half - Ps[R]) / (Ps[R - 1] - Ps[R]) * dk if Ps[R] <= half < Ps[R - 1] else ks[R]
        feat.update(fft_peak_k=float(ks[p]), fft_peak_wavelength=float(2 * np.pi / ks[p]), fft_peak_power=float(Ps[p]), fft_peak_prominence=pr,
                    fft_peak_prominence_rel=pr / float(Ps[p]) if Ps[p] > 0 else np.nan, fft_peak_width_k=width_k,
                    fft_peak_width_rel=width_k / float(ks[p]), fft_peak_fwhm_k=float(kr - kl), fft_peak_fwhm_rel=float((kr - kl) / ks[p]),
                    fft_peak_fwhm_border=border, fft_power_in_band=float((P[band] * cnt[band]).sum() / max(total, 1e-300)))
    curve = [dict(k=float(kk), wavelength=(np.inf if j == 0 else float(2 * np.pi / kk)), P=float(pp), n_modes=int(cc)) for j, (kk, pp, cc) in enumerate(zip(k, P, cnt))]
    return feat, curve


# ---------------------------------------------------------------------------------------------- domain tracking
def track(frames_xy, frames_area, times):
    """Greedy nearest-centroid linking (one-to-one, within LINK_R) of real domains across consecutive recorded frames."""
    tracks, prev_ids, prev_xy, next_id = {}, [], np.empty((0, 2)), 0
    births, deaths = [], []
    for k, (xy, area, t) in enumerate(zip(frames_xy, frames_area, times)):
        ids = [-1] * len(xy)
        if len(prev_xy) and len(xy):
            dist, j = cKDTree(prev_xy).query(xy, k=1)
            used = set()
            for n in np.argsort(dist):
                if dist[n] <= LINK_R and j[n] not in used:
                    ids[n] = prev_ids[j[n]]; used.add(j[n])
        nb = 0
        for n in range(len(xy)):
            if ids[n] == -1:
                ids[n] = next_id; tracks[next_id] = dict(birth=t, last=t, areas=[], xy=[]); next_id += 1; nb += 1
            tr = tracks[ids[n]]; tr["last"] = t; tr["areas"].append(float(area[n])); tr["xy"].append((float(xy[n, 0]), float(xy[n, 1])))
        if k > 0:
            births.append((t, nb)); deaths.append((t, len(set(prev_ids) - set(ids))))
        prev_ids, prev_xy = ids, xy
    return tracks, births, deaths


def track_rows(tracks, births, deaths, t_last, meta):
    rows, alive_drift, completed = [], [], []
    for tid, tr in tracks.items():
        p = np.array(tr["xy"]); steps = np.linalg.norm(np.diff(p, axis=0), axis=1) if len(p) > 1 else np.zeros(0)
        alive = tr["last"] >= t_last - 1e-6
        life = tr["last"] - tr["birth"] + 10.0
        drift = float(np.linalg.norm(p[-1] - p[0]))
        rows.append(dict(**meta, track=tid, birth_t=tr["birth"], last_t=tr["last"], alive_at_end=int(alive), lifetime=life, n_frames=len(p),
                         area_first=tr["areas"][0], area_max=max(tr["areas"]), area_last=tr["areas"][-1], x0=p[0, 0], y0=p[0, 1],
                         x_last=p[-1, 0], y_last=p[-1, 1], drift=drift, path_length=float(steps.sum()), max_step=float(steps.max()) if steps.size else 0.0))
        (alive_drift if alive else completed).append(drift if alive else life)
    bt = np.array([b for t, b in births if t >= 50 - 1e-6]); dt = np.array([d for t, d in deaths if t >= 50 - 1e-6])
    long = [r for r in rows if r["alive_at_end"] and r["birth_t"] <= 60.0]
    summary = dict(**meta, n_tracks=len(tracks), n_alive_end=int(sum(r["alive_at_end"] for r in rows)),
                   births_per_10_t50plus=float(bt.mean()) if bt.size else np.nan, deaths_per_10_t50plus=float(dt.mean()) if dt.size else np.nan,
                   births_total=int(sum(b for _, b in births)), deaths_total=int(sum(d for _, d in deaths)),
                   n_completed=len(completed), median_lifetime_completed=float(np.median(completed)) if completed else np.nan,
                   frac_completed_le20=float(np.mean(np.array(completed) <= 20)) if completed else np.nan,
                   median_drift_alive=float(np.median(alive_drift)) if alive_drift else np.nan, max_drift_alive=float(max(alive_drift)) if alive_drift else np.nan,
                   n_alive_since_t60=len(long), median_drift_alive_since_t60=float(np.median([r["drift"] for r in long])) if long else np.nan,
                   max_drift_alive_since_t60=float(max(r["drift"] for r in long)) if long else np.nan)
    return rows, summary


# ---------------------------------------------------------------------------------------------- one run
def run_one(job):
    """Wrapper: a failure in one run is written to parts/<tag>_error.txt and returned as a row, so the pool keeps going."""
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
    from simulation_2D import run_coupled_hex
    bi, anchor, D = set_params(setname)
    p = circuit(BA, bi, GAMMA, NA, NI, D); c = CONDS[cond]
    thr_fixed = {"a": 0.3 * p["act_prod_rate"], "i": 0.3 * p["inh_prod_rate"] / p["inh_decay_rate"]}
    steps = int(round(T_END / DT)) + 1
    out_dir.mkdir(parents=True, exist_ok=True); (out_dir / "parts").mkdir(exist_ok=True)
    tag = f"{setname}_{cond}_rep{rep}"
    meta = dict(set=setname, b_i=bi, D_i=D, anchor=anchor, condition=cond, replicate=rep)
    spike_seed = ""
    np.random.seed(rep)
    orig_rng = np.random.default_rng
    if c["init_mode"] == "activator_random_spikes":
        spike_seed = central_seed(N) if c["n_points"] == 1 else 1000 * rep + 7
        np.random.default_rng = lambda seed=None, _s=spike_seed: orig_rng(_s)
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            A, I, step, a_ss, i_ss = run_coupled_hex(N, N, steps, DT, 1.0, p, 0.0, 10 ** 9, init_mode=c["init_mode"], activator_type="juxtacrine",
                                                     n_points=c["n_points"], spike_value=c["spike_value"], save_every=SAVE_EVERY,
                                                     nucleation_rate=c["nucleation_rate"], noise_amplitude=0.0, set_peak_height=c["set_peak_height"])
    finally:
        np.random.default_rng = orig_rng
    t_engine = time.perf_counter() - t0
    assert len(A) == steps // SAVE_EVERY + 2 and abs(frame_time(len(A) - 1) - (T_END + DT)) < 1e-9
    A = [np.asarray(f, dtype=np.float64) for f in A]; I = [np.asarray(f, dtype=np.float64) for f in I]
    n2 = N * N
    crit_a = np.array([np.abs(A[k] - A[k - 1]).sum() / n2 for k in range(1, len(A))])
    crit_i = np.array([np.abs(I[k] - I[k - 1]).sum() / n2 for k in range(1, len(I))])
    crit = 0.5 * (crit_a + crit_i)                       # the engine's stopping quantity: L1 change per cell-field over the window
    crit_t = np.array([frame_time(k) for k in range(1, len(A))])
    win = window(N, N); lattice_xy = np.column_stack([c_.ravel() for c_ in cell_centres(N, N)])
    fields = {"a": A, "i": I}
    rec = [0] + [k for k in range(1, len(A)) if (k - 1) % FRAME_EVERY == 0]
    # the engine output is saved before any measurement: every recorded frame (t = 0, 0.01, 10, ..., 600) of both channels in float32,
    # the snapshots at SNAP_T, the final state in float64 and the convergence series, so that everything below can be recomputed
    snaps = {f"{ch}_t{round(frame_time(k))}": H[k].astype(np.float32) for ch, H in fields.items() for k in rec if k >= 1 and round(frame_time(k)) in SNAP_T}
    np.savez_compressed(out_dir / f"{tag}.npz", a_final=A[-1], i_final=I[-1], t_final=frame_time(len(A) - 1), steps=steps, step_returned=step,
                        crit=crit, crit_a=crit_a, crit_i=crit_i, crit_t=crit_t, a_ss=a_ss, i_ss=i_ss, thr_a=thr_fixed["a"], thr_i_fixed=thr_fixed["i"],
                        params=json.dumps(p), condition=json.dumps(c), spike_seed=str(spike_seed), np_seed=rep, grid=N,
                        frame_t=np.array([frame_time(k) for k in rec]), frame_k=np.array(rec),
                        a_frames=np.stack([A[k] for k in rec]).astype(np.float32), i_frames=np.stack([I[k] for k in rec]).astype(np.float32), **snaps)
    t_saved = time.perf_counter() - t0
    try:
        rows_out = measure_run(fields, rec, crit, crit_a, crit_i, thr_fixed, win, lattice_xy, meta, rep)
    except Exception as exc:           # the run is on disk; record the failure and keep the batch going
        import traceback
        (out_dir / "parts" / f"{tag}_error.txt").write_text(traceback.format_exc())
        print(f"MEASUREMENT FAILED {tag}: {exc!r} (engine output saved)", flush=True)
        return dict(**meta, error=repr(exc), wall_engine_s=round(t_engine, 1))
    frame_rows, gr_rows, gr_curves, fft_rows, fft_curves, track_rows_all, track_summ = rows_out
    t_total = time.perf_counter() - t0
    last = {r["channel"]: r for r in frame_rows if r["frame"] == rec[-1]}
    last_i_level = last["i"]["threshold"]
    run_row = dict(**meta, init_mode=c["init_mode"], n_points=c["n_points"], spike_value=c["spike_value"], set_peak_height=c["set_peak_height"],
                   nucleation_rate=c["nucleation_rate"], spike_seed=spike_seed, np_seed=rep, a_ss=float(a_ss), i_ss=float(i_ss), **protocol_cols(a_ss, c),
                   thr_a=thr_fixed["a"], thr_i_fixed=thr_fixed["i"], thr_i_halfmax_final=last_i_level,
                   grid=N, steps=steps, t_final=round(frame_time(len(A) - 1), 2), wall_engine_s=round(t_engine, 1), wall_save_s=round(t_saved - t_engine, 1),
                   wall_total_s=round(t_total, 1), crit_final=float(crit[-1]), crit_min_after_t10=float(crit[5:].min()),
                   **{f"{ch}_{k}": last[ch][k] for ch in ("a", "i") for k in ("cells_on", "coverage", "blobs_all", "blobs_real", "n_transient", "largest",
                                                                             "area_real_median", "area_real_cv", "nn_cv_border", "csr_cv_border", "ecc_mean", "total")})
    write_parts(out_dir, tag, frame_rows, gr_rows, gr_curves, fft_rows, fft_curves, track_rows_all, track_summ, run_row)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(f"done {tag}: engine {t_engine:.0f} s, total {t_total:.0f} s, peak RSS {rss:.2f} GB, a: {last['a']['blobs_real']} real / {last['a']['n_transient']} transient, "
          f"cov {last['a']['coverage']:.3f}; i: {last['i']['blobs_real']} real, cov {last['i']['coverage']:.3f}", flush=True)
    return run_row


def protocol_cols(a_ss, c):
    """What the engine actually does with the nucleation and spike settings: kick 2 a_ss gated at a < a_ss; a_ss from the solver or the fallback."""
    fallback = abs(float(a_ss) - float(c["spike_value"])) < 1e-12
    return dict(a_ss_source="fallback (= spike_value)" if fallback else "solver", nucl_kick=2 * float(a_ss) if c["nucleation_rate"] > 0 else 0.0,
                nucl_gate=float(a_ss) if c["nucleation_rate"] > 0 else np.nan)


def write_parts(out_dir, tag, frame_rows, gr_rows, gr_curves, fft_rows, fft_curves, track_rows_all, track_summ, run_row):
    for name, rows in (("frames", frame_rows), ("gr", gr_rows), ("gr_curves", gr_curves), ("fft", fft_rows), ("fft_curves", fft_curves),
                       ("tracks", track_rows_all), ("tracking", track_summ), ("runs", [run_row])):
        if rows:
            cols = list(dict.fromkeys(k for r in rows for k in r))        # ordered union: inhibitor rows carry extra i_fixed_* columns
            with open(out_dir / "parts" / f"{tag}_{name}.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n"); w.writeheader(); w.writerows(rows)


def remeasure(path):
    """Recompute every measurement of a saved run from its stored frames (after a change to the measurement code); the runs part keeps
    its engine columns and gets the measurement columns refreshed. A failure is written to parts/<tag>_remeasure_error.txt."""
    try:
        _remeasure(path)
    except Exception as exc:
        import traceback
        (path.parent / "parts" / f"{path.stem}_remeasure_error.txt").write_text(traceback.format_exc())
        print(f"REMEASURE FAILED {path.name}: {exc!r}", flush=True)


def _remeasure(path):
    d = np.load(path); tag = path.stem; out_dir = path.parent
    setname, cond, rep = tag.rsplit("_rep", 1)[0].rsplit("_", 1)[0], tag.rsplit("_rep", 1)[0].rsplit("_", 1)[1], int(tag.rsplit("_rep", 1)[1])
    bi, anchor, D = set_params(setname); meta = dict(set=setname, b_i=bi, D_i=D, anchor=anchor, condition=cond, replicate=rep)
    rec = [int(k) for k in d["frame_k"]]
    fields = {"a": dict(zip(rec, d["a_frames"].astype(np.float64))), "i": dict(zip(rec, d["i_frames"].astype(np.float64)))}
    n = int(d["grid"]); win = window(n, n); lattice_xy = np.column_stack([c_.ravel() for c_ in cell_centres(n, n)])
    thr_fixed = {"a": float(d["thr_a"]), "i": float(d["thr_i_fixed"])}
    out = measure_run(fields, rec, d["crit"], d["crit_a"], d["crit_i"], thr_fixed, win, lattice_xy, meta, rep)
    frame_rows = out[0]
    last = {r["channel"]: r for r in frame_rows if r["frame"] == rec[-1]}
    runs_part = out_dir / "parts" / f"{tag}_runs.csv"
    run_row = dict(next(csv.DictReader(open(runs_part)))) if runs_part.exists() else dict(**meta)
    run_row.update(**protocol_cols(float(d["a_ss"]), json.loads(str(d["condition"]))), thr_i_halfmax_final=last["i"]["threshold"],
                   **{f"{ch}_{k}": last[ch][k] for ch in ("a", "i") for k in ("cells_on", "coverage", "blobs_all", "blobs_real", "n_transient", "largest",
                                                                             "area_real_median", "area_real_cv", "nn_cv_border", "csr_cv_border", "ecc_mean", "total")})
    write_parts(out_dir, tag, *out, run_row)
    print(f"remeasured {tag}", flush=True)


def measure_run(fields, rec, crit, crit_a, crit_i, thr_fixed, win, lattice_xy, meta, rep):
    """All measurements of one run from its saved frames; returns the row lists for the parts files."""
    frame_rows, gr_rows, gr_curves, fft_rows, fft_curves, track_rows_all, track_summ = [], [], [], [], [], [], []
    for ch, H in fields.items():
        csr = CsrCache(win)
        xy_seq, area_seq, t_seq = [], [], []
        i_floor = thr_fixed["i"] / 30.0          # 0.01 r_i / gamma: an inhibitor field below this is off (no half-max segmentation of a decayed field)
        def level(k):      # activator: fixed 0.3 r_a; inhibitor: half of the frame's own maximum (persistence uses the previous frame's own level)
            return thr_fixed["a"] if ch == "a" else (0.5 * float(H[k].max()) if H[k].max() >= i_floor else np.inf)
        for k in rec:
            t = frame_time(k); tr = round(t)
            kp = k - FRAME_EVERY if k - FRAME_EVERY >= FRAME_EVERY + 1 else None      # partner = the t - 10 frame; the t = 10 row has none (t = 0.01 is pre-pattern)
            m, comp, xy = measure_frame(H[k], H[kp] > level(kp) if kp is not None else None, level(k), win, csr)
            row = dict(**meta, channel=ch, threshold=level(k), threshold_kind="fixed 0.3 r_a" if ch == "a" else "half of frame max", frame=k, t=round(t, 2),
                       persistence=int(kp is not None), crit=float(crit[k - 1]) if k >= 1 else np.nan,
                       crit_channel=float((crit_a if ch == "a" else crit_i)[k - 1]) if k >= 1 else np.nan,
                       crit_window_steps=(1 if k == 1 else SAVE_EVERY) if k >= 1 else 0, **m)
            if ch == "i":
                fm = H[k] > thr_fixed["i"]
                fc = components(fm)
                row.update(i_fixed_level=thr_fixed["i"], i_fixed_cells=int(fm.sum()), i_fixed_coverage=float(fm.mean()), i_fixed_blobs=int(len(fc["area"])),
                           i_fixed_blobs_real=int((fc["area"] >= MIN_REAL).sum()), i_fixed_largest=int(fc["area"].max()) if fc["area"].size else 0)
            frame_rows.append(row)
            if k >= FRAME_EVERY + 1:                 # tracking from t = 10 on
                xy_seq.append(xy); area_seq.append(comp["area"][comp["area"] >= MIN_REAL] if comp["area"].size else np.zeros(0)); t_seq.append(tr)
            if k >= 1 and tr in GR_T:
                med = m["area_real_median"]
                null = "lattice" if (np.isfinite(med) and med < 10) else "uniform"
                feat, curve = gr_features(xy, win, lattice_xy if null == "lattice" else None, seed=rep * 1000 + tr)
                gr_rows.append(dict(**meta, channel=ch, t=tr, n_points=len(xy), null=null, **feat))
                gr_curves += [dict(**meta, channel=ch, t=tr, **cv) for cv in curve]
                ff, fc = hex_power_spectrum(H[k])
                fft_rows.append(dict(**meta, channel=ch, t=tr, field_mean=float(H[k].mean()), field_var=float(H[k].var()), **ff))
                fft_curves += [dict(**meta, channel=ch, t=tr, **cv) for cv in fc]
        tracks, births, deaths = track(xy_seq, area_seq, t_seq)
        rows, summ = track_rows(tracks, births, deaths, t_seq[-1], dict(**meta, channel=ch))
        track_rows_all += rows; track_summ.append(summ)
    return frame_rows, gr_rows, gr_curves, fft_rows, fft_curves, track_rows_all, track_summ


def parse_sets(args):
    bad = [a for a in args if a not in SETS]
    if bad:
        raise SystemExit(f"unknown set(s) {bad}; known: {list(SETS)}")
    return list(args) or ORDER


def jobs_for(sets, conds=None):
    return [(s, cond, rep, RES) for s in sets for cond in (conds or CONDS) for rep in REPS[cond]]


def collect(out_dir=RES, prefix="ladder"):
    parts = out_dir / "parts"
    for name in ("runs", "frames", "gr", "gr_curves", "fft", "fft_curves", "tracks", "tracking"):
        rows, cols = [], []
        for f in sorted(parts.glob(f"*_{name}.csv")):
            with open(f) as fh:
                rd = csv.DictReader(fh)
                for r in rd:
                    rows.append(r)
                for ccol in rd.fieldnames or []:
                    if ccol not in cols:
                        cols.append(ccol)
        if rows:
            with open(SUMM / f"{prefix}_{name}.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n"); w.writeheader(); w.writerows(rows)
            print(f"{prefix}_{name}.csv: {len(rows)} rows")


def figures(sets):
    from analysis.plot_hex import save_panels
    rows_cond = [("seed1", 0), ("seed9", 1), ("seed36", 1), ("seed144", 1), ("random", 1), ("nucl", 1)]
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
                save_panels(fields, titles, SUMM / f"ladder_{s}_{ch}.png", ncols=len(FIG_T),
                            suptitle=f"{s} ({anchor}): b_a=5, b_i={bi:g}, n=10/4, D_i={D:g}, gamma=0.5; 200x200; {label}, hex-correct, colour 0..max")
                print(f"ladder_{s}_{ch}.png")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "time":
        out = RES / "timing"
        t0 = time.perf_counter()
        r = run_one(("bi25", "nucl", 1, out))
        print(f"wall {time.perf_counter() - t0:.0f} s for one 200x200 nucleation run with the full pipeline (engine {r['wall_engine_s']} s)")
    elif cmd == "run":
        args = sys.argv[2:]
        conds = None
        if "--conds" in args:                       # e.g. --conds seed1,seed36,random,nucl
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
