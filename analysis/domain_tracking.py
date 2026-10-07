"""
Track domains across saved frames (t spacing 10) to measure lifetimes and turnover under continuous
nucleation, where a snapshot count cannot tell a frozen pattern from a dynamic steady state.

Segmentation: fixed threshold 0.3*act_prod_rate, components >= MIN_CELLS cells. Linking: nearest
centroid within LINK_R cells between consecutive frames (domains do not move; LINK_R allows boundary
flicker). Reports, per run: births and deaths per 10 time units in the steady window, the lifetime
distribution of completed tracks, and the age of domains alive at the last frame.
Usage: python analysis/domain_tracking.py <results subdir> [<subdir> ...]
Output: analysis/summaries/domain_tracking.csv (appended rows per run)
"""
from __future__ import annotations

import csv, sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from analysis.hexgeom import component_table  # noqa: E402

MIN_CELLS, LINK_R, T_START = 6, 3.0, 50.0


def track_run(path: Path):
    d = np.load(path)
    t = d["t"]; r_a = float(d["act_prod_rate"])
    frames = [k for k in range(len(t)) if t[k] >= T_START - 1e-6]
    tracks = {}            # id -> dict(birth, last, areas)
    prev_ids, prev_xy = [], np.empty((0, 2))
    next_id = 0
    births, deaths = [], []
    for k in frames:
        tab = component_table(d["a"][k].astype(float) > 0.3 * r_a, min_size=MIN_CELLS)
        xy = np.column_stack([tab["cx"], tab["cy"]])
        ids = [-1] * len(xy)
        if len(prev_xy) and len(xy):
            dist, j = cKDTree(prev_xy).query(xy, k=1)
            used = set()
            for n in np.argsort(dist):               # greedy nearest linking, one-to-one
                if dist[n] <= LINK_R and j[n] not in used:
                    ids[n] = prev_ids[j[n]]; used.add(j[n])
        nb = 0
        for n in range(len(xy)):
            if ids[n] == -1:
                ids[n] = next_id; tracks[next_id] = dict(birth=t[k], last=t[k], areas=[]); next_id += 1; nb += 1
            tracks[ids[n]]["last"] = t[k]; tracks[ids[n]]["areas"].append(float(tab["area"][n]))
        nd = len(set(prev_ids) - set(ids))
        if k != frames[0]:
            births.append(nb); deaths.append(nd)
        prev_ids, prev_xy = ids, xy
    t_last = t[frames[-1]]
    completed = [tr["last"] - tr["birth"] + 10 for tr in tracks.values() if tr["last"] < t_last - 1e-6]
    alive = [tr for tr in tracks.values() if tr["last"] >= t_last - 1e-6]
    ages = [t_last - tr["birth"] + 10 for tr in alive]
    return dict(file=path.name, n_frames=len(frames), n_tracks=len(tracks), mean_alive=float(np.mean([len(set(ids_)) for ids_ in [prev_ids]])) if prev_ids else 0,
                n_alive_end=len(alive), births_per_10=float(np.mean(births)) if births else 0.0, deaths_per_10=float(np.mean(deaths)) if deaths else 0.0,
                n_completed=len(completed), median_lifetime_completed=float(np.median(completed)) if completed else np.nan,
                frac_completed_le20=float(np.mean(np.array(completed) <= 20)) if completed else np.nan,
                frac_alive_since_t50=float(np.mean([tr["birth"] <= T_START + 0.1 for tr in alive])) if alive else np.nan,
                median_age_alive=float(np.median(ages)) if ages else np.nan,
                mean_area_alive=float(np.mean([np.mean(tr["areas"]) for tr in alive])) if alive else np.nan)


if __name__ == "__main__":
    rows = []
    for sub in sys.argv[1:]:
        for p in sorted((ROOT / "results" / sub).glob("*.npz")):
            r = track_run(p); r["subdir"] = sub; rows.append(r)
    cols = ["subdir", "file", "n_alive_end", "births_per_10", "deaths_per_10", "n_tracks", "n_completed", "median_lifetime_completed",
            "frac_completed_le20", "frac_alive_since_t50", "median_age_alive", "mean_area_alive"]
    out = ROOT / "analysis" / "summaries" / "domain_tracking.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    print(f"{'run':<34}{'alive@end':>10}{'births/10':>10}{'deaths/10':>10}{'tracks':>7}{'done':>6}{'med life':>9}{'<=20':>6}{'old':>6}{'med age':>8}{'area':>7}")
    for r in rows:
        print(f"{r['subdir'][:12]+'/'+r['file'][:21]:<34}{r['n_alive_end']:>10}{r['births_per_10']:>10.2f}{r['deaths_per_10']:>10.2f}{r['n_tracks']:>7}{r['n_completed']:>6}"
              f"{r['median_lifetime_completed']:>9.0f}{r['frac_completed_le20']:>6.2f}{r['frac_alive_since_t50']:>6.2f}{r['median_age_alive']:>8.0f}{r['mean_area_alive']:>7.1f}")
    print("\ncolumns: alive@end = domains at t=500; births/deaths per 10 time units over t=50..500; tracks = distinct domains seen;")
    print("done = tracks that ended before t=500 with their median lifetime and the fraction lasting <= 20; old = fraction of domains alive at")
    print("t=500 that already existed at t=50; med age = median age of domains alive at t=500 (max 460 = present since t=50).")
