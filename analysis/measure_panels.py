"""
Measure the Fig 1C panel runs produced by run_fig1c_panels.py with the validated hex-aware pipeline.

Per run and time point: coverage, component count (two size cuts), area mean/CV, NN-spacing CV
(plain and border-corrected) with a CSR reference computed for the SAME n in the SAME window.
Per (protocol, b_i, t): g(r) pooled over seeds by summing counts; hole, first peak, trough, second
peak; compared with 2*lambda_i nominal (sqrt(D/g)) and lattice-corrected (sqrt(1.5 D/g)).
Usage: python analysis/measure_panels.py <results subdir> [halfmax|fixed:<frac of r_a>] [persist]
  halfmax      threshold 0.5*max(a) per frame (default; under continuous nucleation fresh kicks set max(a))
  fixed:0.3    threshold 0.3*act_prod_rate, independent of the kick amplitude
  persist      a cell counts only if above threshold in this frame AND the previous saved frame (t-10),
               which removes transient nucleation events (they decay within ~1 time unit)
Outputs carry a suffix naming the mode, e.g. <subdir>_fixed0.3_persist_measurements.csv
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
from analysis.pointstats import PooledG, border_cv, csr_reference_cv, nn_cv  # noqa: E402

PANEL = sys.argv[1] if len(sys.argv) > 1 else "fig1c_panels"
RUNS = ROOT / "results" / PANEL
SEG = sys.argv[2] if len(sys.argv) > 2 else "halfmax"
PERSIST = "persist" in sys.argv[3:]
SUFFIX = "" if (SEG == "halfmax" and not PERSIST) else "_" + SEG.replace(":", "") + ("_persist" if PERSIST else "")
T_POINTS = (10, 50, 100, 200, 350, 500)
MIN_SIZE = {"synchronous": 1, "nucleation": 12}      # handoff §5: size filter mandatory under continuous nucleation


def segment_mode(a, r_a):
    if SEG == "halfmax":
        return segment(a, method="half_max")
    if SEG.startswith("fixed:"):
        return segment(a, method="fixed", level=float(SEG.split(":")[1]) * r_a)
    raise ValueError(SEG)


def measure_frame(a, protocol, a_prev=None, r_a=5.0):
    ny, nx = a.shape
    mask, thr = segment_mode(a, r_a)
    if PERSIST and a_prev is not None:
        mask &= segment_mode(a_prev, r_a)[0]
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
            a_prev = d["a"][k - 1].astype(float) if k > 0 else None
            m, xy, win = measure_frame(d["a"][k].astype(float), protocol, a_prev, r_a=float(d["act_prod_rate"]))
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
    out = ROOT / "analysis" / "summaries" / f"{PANEL}{SUFFIX}_measurements.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)

    g_rows = []
    lam_nom, lam_lat = np.sqrt(D_I / GAMMA), np.sqrt(1.5 * D_I / GAMMA)
    print(f"segmentation: {SEG}{' + persistence' if PERSIST else ''}; size cut {MIN_SIZE}")
    print(f"2*lambda_i nominal = {2*lam_nom:.2f} cells, lattice-corrected = {2*lam_lat:.2f} cells\n")
    print(f"{'protocol':<12}{'b_i':>4}{'t':>5}{'seeds':>6}{'n/seed':>7}{'pairs':>7}{'hole':>7}{'+-':>5}{'peak1_r':>8}{'+-':>5}{'peak1_g':>8}{'pairs':>6}{'trough_r':>9}{'trough_g':>9}{'peak2_r':>8}{'peak2_g':>8}")
    for key in sorted(pooled):
        G = pooled[key]; s = G.summary(); protocol, bi, tp = key
        g_rows.append(dict(protocol=protocol, bi=bi, t=tp, **{k_: round(v, 3) for k_, v in s.items()}))
        print(f"{protocol:<12}{bi:>4}{tp:>5}{len(G.n_points):>6}{np.mean(G.n_points) if G.n_points else 0:>7.1f}{s['total_pairs']:>7.0f}"
              f"{s['hole']:>7.2f}{s['hole_sd']:>5.2f}{s['peak1_r']:>8.2f}{s['peak1_r_sd']:>5.2f}{s['peak1_g']:>8.2f}{s['peak1_pairs']:>6.0f}{s['trough_r']:>9.2f}{s['trough_g']:>9.2f}{s['peak2_r']:>8.2f}{s['peak2_g']:>8.2f}")
        with open(ROOT / "analysis" / "summaries" / f"{PANEL}{SUFFIX}_gr_{protocol}_bi{bi}_t{tp}.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n"); w.writerow(["r", "g", "h_data", "h_null_mean"])
            for r_, g_, hd, hn in zip(G.r, G.g(), G.h_data, G.h_null):
                w.writerow([r_, g_, hd, hn])
    with open(ROOT / "analysis" / "summaries" / f"{PANEL}{SUFFIX}_gr.csv", "w", newline="") as fh:
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
