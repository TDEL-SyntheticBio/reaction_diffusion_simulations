"""
Fate of an ISOLATED 19-cell seed along the Fig 1C line (b_a=5, gamma=0.5, n=3/3, D_i=10), run to t=600 on
141x141 with no nucleation. Records area and component count every t=25 and snapshots at t=100, 300, 600,
to classify each b_i as: stable disc (1 component, constant area), ring (1 component, hollow), fragmented
(>1 component, growth stopped) or self-replicating colony (component count still rising).
Outputs: analysis/summaries/isolated_seed_fate.csv (time series), results/isolated_seed/bi<bi>.npz,
         analysis/summaries/isolated_seed_fate.png
"""
from __future__ import annotations

import csv, sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "2D_simulations"))
from simulation_2D import _vectorized_step  # noqa: E402
from analysis.hexgeom import cell_centres, neighbor_arrays, rasterize_discs, label_components  # noqa: E402
from analysis.single_domain_halo import params  # noqa: E402

OUT = ROOT / "results" / "isolated_seed"
N = 141
T_END, DT, EVERY = 600.0, 0.01, 25.0
SNAPS = (100.0, 300.0, 600.0)


def one(bi):
    p = params(5, bi, 0.5, 3, 3, 10)
    nbr = neighbor_arrays(N, N)
    x, y = cell_centres(N, N)
    a = np.where(rasterize_discs(N, N, np.array([[x[N // 2, N // 2], y[N // 2, N // 2]]]), 2.5), 5.0, 0.0)
    i = np.zeros_like(a)
    series, snaps = [], {}
    for step in range(1, int(round(T_END / DT)) + 1):
        a, i = _vectorized_step(a, i, DT, 1.0, p, "juxtacrine", *nbr)
        t = round(step * DT, 2)
        if step % int(EVERY / DT) == 0:
            mask = a > 0.5 * a.max() if a.max() > 0.05 else np.zeros_like(a, bool)
            ncomp = label_components(mask)[1] if mask.any() else 0
            # hollow test: fraction of the mask's bounding-disc interior that is active
            if mask.any():
                xs, ys = x[mask], y[mask]; cx, cy = xs.mean(), ys.mean(); rmax = np.hypot(xs - cx, ys - cy).max()
                inside = np.hypot(x - cx, y - cy) <= rmax; fill = mask[inside].mean()
            else:
                fill = 0.0
            touches = bool(mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any())
            series.append(dict(bi=bi, t=t, area=int(mask.sum()), ncomp=ncomp, fill_fraction=round(float(fill), 3), touches_edge=touches, a_max=round(float(a.max()), 3)))
        if t in SNAPS:
            snaps[t] = a.copy()
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"bi{bi}.npz", **{f"a_t{int(t)}": s for t, s in snaps.items()})
    return series, snaps


if __name__ == "__main__":
    bis = [6, 7, 8, 9, 10, 11]
    with Pool(2) as pool:
        res = pool.map(one, bis, chunksize=1)
    rows = [r for series, _ in res for r in series]
    with open(ROOT / "analysis" / "summaries" / "isolated_seed_fate.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"{'b_i':>4} " + " ".join(f"{'t='+str(int(t)):>12}" for t in (100, 200, 300, 400, 500, 600)) + "   (area/ncomp/fill)")
    for (series, _), bi in zip(res, bis):
        cells = {r["t"]: r for r in series}
        print(f"{bi:>4} " + " ".join(f"{cells[t]['area']:>5}/{cells[t]['ncomp']:>2}/{cells[t]['fill_fraction']:.2f}" for t in (100.0, 200.0, 300.0, 400.0, 500.0, 600.0)))
    from analysis.plot_hex import save_panels
    fields, titles = [], []
    for (series, snaps), bi in zip(res, bis):
        for t in SNAPS:
            fields.append(snaps[t]); titles.append(f"b_i={bi}, t={int(t)}")
    save_panels(fields, titles, ROOT / "analysis" / "summaries" / "isolated_seed_fate.png", ncols=3,
                suptitle="Isolated 19-cell seed at level 5, b_a=5, gamma=0.5, n=3/3, D_i=10, no nucleation (141x141, hex-correct)")
    print("done")
