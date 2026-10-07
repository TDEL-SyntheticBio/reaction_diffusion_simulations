"""
Measure the Fig 1C panel runs produced by run_fig1c_panels.py with the validated hex-aware pipeline.

Per run and time point: coverage, component count (two size cuts), area mean/CV, NN-spacing CV
(plain and border-corrected) with a CSR reference computed for the SAME n in the SAME window.
Per (protocol, b_i, t): g(r) pooled over seeds by summing counts; hole, first peak, trough, second
peak; compared with 2*lambda_i nominal (sqrt(D/g)) and lattice-corrected (sqrt(1.5 D/g)).
Usage: python analysis/measure_panels.py <results subdir>   (default fig1c_panels)
Outputs: analysis/summaries/<subdir>_measurements.csv, <subdir>_gr.csv, <subdir>_gr_<protocol>_bi<bi>_t<t>.csv
"""
from __future__ import annotations

import csv, sys, zlib
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from analysis.hexgeom import component_table, segment, window  # noqa: E402
from analysis.pointstats import PooledG, nn_cv, nn_spacing  # noqa: E402

PANEL = sys.argv[1] if len(sys.argv) > 1 else "fig1c_panels"
RUNS = ROOT / "results" / PANEL
T_POINTS = (10, 50, 100, 200, 350, 500)
MIN_SIZE = {"synchronous": 1, "nucleation": 12}      # handoff §5: size filter mandatory under continuous nucleation


def border_cv(xy, win):
    xmin, xmax, ymin, ymax = win
    if len(xy) < 3:
        return np.nan, 0
    d = nn_spacing(xy)
    edge = np.minimum.reduce([xy[:, 0] - xmin, xmax - xy[:, 0], xy[:, 1] - ymin, ymax - xy[:, 1]])
    keep = d < edge
    return (float(d[keep].std() / d[keep].mean()) if keep.sum() > 2 else np.nan), int(keep.sum())


def csr_reference_cv(n, win, rng, n_real=300):
    """Plain and border-corrected NN CV of n uniform points in this window (finite-n, finite-window reference)."""
    xmin, xmax, ymin, ymax = win
    plain, border = [], []
    for _ in range(n_real):
        pts = np.column_stack([rng.uniform(xmin, xmax, n), rng.uniform(ymin, ymax, n)])
        plain.append(nn_cv(pts)); border.append(border_cv(pts, win)[0])
    return float(np.nanmean(plain)), float(np.nanmean(border))


def measure_frame(a, protocol):
    ny, nx = a.shape
    mask, thr = segment(a, method="half_max")
    cov = float(mask.mean())
    tab_all = component_table(mask, min_size=1)
    tab = component_table(mask, min_size=MIN_SIZE[protocol])
    xy = np.column_stack([tab["cx"], tab["cy"]])
    win = window(ny, nx)
    bcv, nkeep = border_cv(xy, win)
    return dict(coverage=cov, thr=thr, n_all=tab_all["n_kept"], n=tab["n_kept"], area_mean=float(tab["area"].mean()) if tab["n_kept"] else np.nan,
                area_cv=float(tab["area"].std() / tab["area"].mean()) if tab["n_kept"] > 1 else np.nan,
                cv_plain=nn_cv(xy) if len(xy) > 2 else np.nan, cv_border=bcv, n_border=nkeep, a_max=float(a.max())), xy, win


def main():
    files = sorted(RUNS.glob("*.npz"))
    if not files:
        sys.exit("no runs found")
    rng = np.random.default_rng(0)
    rows, pooled = [], {}
    csr_cache = {}
    for f in files:
        d = np.load(f, allow_pickle=False)
        protocol = str(d["protocol"]); bi = int(d["inh_prod_rate"]); seed = int(d["seed"])
        D_I, GAMMA = float(d["inh_diffusion"]), float(d["inh_decay_rate"])
        t = d["t"]
        for tp in T_POINTS:
            k = int(np.argmin(np.abs(t - tp)))
            if abs(t[k] - tp) > 0.02:
                continue
            m, xy, win = measure_frame(d["a"][k].astype(float), protocol)
            n = len(xy)
            if n > 2:
                if n not in csr_cache:
                    csr_cache[n] = csr_reference_cv(n, win, rng)
                m["csr_cv_plain"], m["csr_cv_border"] = csr_cache[n]
            else:
                m["csr_cv_plain"] = m["csr_cv_border"] = np.nan
            rows.append(dict(protocol=protocol, bi=bi, seed=seed, t=tp, **{k_: (round(v, 4) if isinstance(v, float) else v) for k_, v in m.items()}))
            if tp in (50, 500):
                key = (protocol, bi, tp)
                pooled.setdefault(key, PooledG(r_max=35, dr=0.5, n_null=150, seed=zlib.crc32(repr(key).encode())))   # deterministic across processes
                if n > 2:
                    pooled[key].add(xy, win)
    out = ROOT / "analysis" / "summaries" / f"{PANEL}_measurements.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)

    g_rows = []
    lam_nom, lam_lat = np.sqrt(D_I / GAMMA), np.sqrt(1.5 * D_I / GAMMA)
    print(f"2*lambda_i nominal = {2*lam_nom:.2f} cells, lattice-corrected = {2*lam_lat:.2f} cells\n")
    print(f"{'protocol':<12}{'b_i':>4}{'t':>5}{'seeds':>6}{'n/seed':>7}{'hole':>7}{'peak1_r':>8}{'peak1_g':>8}{'trough_r':>9}{'trough_g':>9}{'peak2_r':>8}{'peak2_g':>8}")
    for key in sorted(pooled):
        G = pooled[key]; s = G.summary(); protocol, bi, tp = key
        g_rows.append(dict(protocol=protocol, bi=bi, t=tp, seeds=len(G.n_points), n_mean=float(np.mean(G.n_points)) if G.n_points else 0, **{k_: round(v, 3) for k_, v in s.items()}))
        print(f"{protocol:<12}{bi:>4}{tp:>5}{len(G.n_points):>6}{np.mean(G.n_points) if G.n_points else 0:>7.1f}{s['hole']:>7.2f}{s['peak1_r']:>8.2f}{s['peak1_g']:>8.2f}{s['trough_r']:>9.2f}{s['trough_g']:>9.2f}{s['peak2_r']:>8.2f}{s['peak2_g']:>8.2f}")
        with open(ROOT / "analysis" / "summaries" / f"{PANEL}_gr_{protocol}_bi{bi}_t{tp}.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n"); w.writerow(["r", "g", "h_data", "h_null_mean"])
            for r_, g_, hd, hn in zip(G.r, G.g(), G.h_data, G.h_null):
                w.writerow([r_, g_, hd, hn])
    with open(ROOT / "analysis" / "summaries" / f"{PANEL}_gr.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(g_rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(g_rows)

    print("\nper-run time series (mean over seeds):")
    import collections
    agg = collections.defaultdict(list)
    for r in rows:
        agg[(r["protocol"], r["bi"], r["t"])].append(r)
    print(f"{'protocol':<12}{'b_i':>4}{'t':>5}{'cov':>7}{'n_all':>7}{'n_cut':>7}{'area':>7}{'areaCV':>8}{'CV':>7}{'CVborder':>9}{'CSR CV':>8}{'CSR bord':>9}")
    for key in sorted(agg):
        rs = agg[key]; mean = lambda k_: np.nanmean([r[k_] for r in rs])
        print(f"{key[0]:<12}{key[1]:>4}{key[2]:>5}{mean('coverage'):>7.3f}{mean('n_all'):>7.1f}{mean('n'):>7.1f}{mean('area_mean'):>7.1f}{mean('area_cv'):>8.3f}{mean('cv_plain'):>7.3f}{mean('cv_border'):>9.3f}{mean('csr_cv_plain'):>8.3f}{mean('csr_cv_border'):>9.3f}")


if __name__ == "__main__":
    main()
