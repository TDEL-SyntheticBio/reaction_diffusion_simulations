"""
Domain-size distribution of the saved t = 600 activator snapshots of the seed-crowding experiments (REPORT §10, §11),
which recorded only count, total area, largest blob, coverage and total activator per frame.
Components of (a > 0.3*r_a) on the engine's hex adjacency (no persistence: only the t = 600 frame was saved);
"real" = components of at least 6 cells. Nothing is re-run.
Usage: python analysis/snapshot_areas.py <results subdir> [<subdir> ...]   (e.g. seed_crowding_exp seed_crowding)
Output: analysis/summaries/<subdir>_areas_t600.csv, one row per run, with the full per-component area list.
"""
from __future__ import annotations
import csv, re, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from analysis.hexgeom import component_table  # noqa: E402


def stats(areas: np.ndarray, prefix: str) -> dict:
    if len(areas) == 0:
        return {f"{prefix}_{k}": np.nan for k in ("n", "mean", "sd", "cv", "median", "p10", "p90", "max")} | {f"{prefix}_n": 0}
    return {f"{prefix}_n": int(len(areas)), f"{prefix}_mean": float(areas.mean()), f"{prefix}_sd": float(areas.std(ddof=1)) if len(areas) > 1 else 0.0,
            f"{prefix}_cv": float(areas.std(ddof=1) / areas.mean()) if len(areas) > 1 else 0.0, f"{prefix}_median": float(np.median(areas)),
            f"{prefix}_p10": float(np.percentile(areas, 10)), f"{prefix}_p90": float(np.percentile(areas, 90)), f"{prefix}_max": int(areas.max())}


def main(subdirs):
    for sub in subdirs:
        rows = []
        for path in sorted((ROOT / "results" / sub).glob("*.npz")):
            m = re.match(r"(.+)_(seed1_L5|seed1|seed9|seed36|seed144|random|nucl)_rep(\d+)\.npz", path.name)
            if not m:
                continue
            d = np.load(path)
            if "a_t600" not in d.files:
                continue
            a = d["a_t600"]; thr = float(d["threshold"]) if "threshold" in d.files else 1.5   # older files: all b_a = 5
            areas = np.sort(component_table(a > thr, min_size=1)["area"].astype(int))[::-1]
            rows.append(dict(set=m.group(1), condition=m.group(2), replicate=int(m.group(3)), threshold=thr, grid=a.shape[0], t=600,
                             cells_on=int((a > thr).sum()), coverage=float((a > thr).mean()),
                             **stats(areas[areas >= 6], "real"), **stats(areas, "all"), n_transient=int((areas < 6).sum()),
                             areas=" ".join(map(str, areas))))
        out = ROOT / "analysis" / "summaries" / f"{sub}_areas_t600.csv"
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
        print(f"{out}: {len(rows)} runs")


if __name__ == "__main__":
    main(sys.argv[1:])
