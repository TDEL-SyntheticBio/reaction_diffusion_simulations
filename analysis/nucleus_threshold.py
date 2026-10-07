"""
(1) Critical nucleus size at the sparse end of the Fig 1C line: plant discs of increasing radius at level 5
    (b_i = 12) and see which survive; also check persistence of the arrested 19-cell domain at b_i = 12, 14,
    16, 20 to t = 300.
(2) Render hex-correct pictures of the isolated-seed outcome at b_i = 5 (spreading colony) and b_i = 8
    (arrested disc) at t = 150 and t = 300.
Outputs: analysis/summaries/nucleus_threshold.csv, results/single_domain/*.npz, analysis/summaries/single_domain_bi5_bi8.png
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

OUT = ROOT / "results" / "single_domain"


def evolve(p, ny, nx, seed_radius, seed_level, t_end, snap_t=()):
    nbr = neighbor_arrays(ny, nx)
    x, y = cell_centres(ny, nx)
    cx, cy = x[ny // 2, nx // 2], y[ny // 2, nx // 2]
    a = np.where(rasterize_discs(ny, nx, np.array([[cx, cy]]), seed_radius), float(seed_level), 0.0)
    i = np.zeros_like(a)
    snaps, areas = {}, []
    steps = int(round(t_end / 0.01))
    for step in range(1, steps + 1):
        a, i = _vectorized_step(a, i, 0.01, 1.0, p, "juxtacrine", *nbr)
        if step % 1000 == 0:
            areas.append(int((a > 0.5 * a.max()).sum()) if a.max() > 0.05 else 0)
        if round(step * 0.01, 2) in snap_t:
            snaps[round(step * 0.01, 2)] = (a.copy(), i.copy())
    return a, i, areas, snaps


def threshold_job(job):
    bi, radius, t_end = job
    a, i, areas, _ = evolve(params(5, bi, 0.5, 3, 3, 10), 101, 101, radius, 5.0, t_end)
    mask = a > 0.5 * a.max() if a.max() > 0.05 else np.zeros_like(a, bool)
    _, ncomp = label_components(mask) if mask.any() else (None, 0)
    seed_cells = int(rasterize_discs(101, 101, np.array([[cell_centres(101, 101)[0][50, 50], cell_centres(101, 101)[1][50, 50]]]), radius).sum())
    outcome = "died" if not mask.any() else ("arrested" if len(set(areas[-5:])) == 1 else "changing")
    return dict(bi=bi, seed_radius=radius, seed_cells=seed_cells, t_end=t_end, outcome=outcome, final_area=int(mask.sum()), ncomp=ncomp,
                a_max=round(float(a.max()), 3), area_hist=str(areas[-5:]))


def picture_job(bi):
    p = params(5, bi, 0.5, 3, 3, 10)
    a, i, areas, snaps = evolve(p, 141, 141, 2.5, 5.0, 300.0, snap_t=(50.0, 150.0, 300.0))
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"single_bi{bi}.npz", **{f"a_t{int(t)}": s[0] for t, s in snaps.items()}, **{f"i_t{int(t)}": s[1] for t, s in snaps.items()},
                        area_hist=np.array(areas))
    return bi, snaps


if __name__ == "__main__":
    jobs = [(12, r, 300.0) for r in (0.6, 1.0, 1.5, 2.0, 2.5)] + [(bi, 2.5, 300.0) for bi in (14, 16, 20)]
    with Pool(2) as pool:
        rows = pool.map(threshold_job, jobs, chunksize=1)
        pics = pool.map(picture_job, [5, 8], chunksize=1)
    with open(ROOT / "analysis" / "summaries" / "nucleus_threshold.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
    from analysis.plot_hex import save_panels
    fields, titles = [], []
    for bi, snaps in pics:
        for t in (50.0, 150.0, 300.0):
            fields.append(snaps[t][0]); titles.append(f"b_i={bi}, t={int(t)}: activator")
    save_panels(fields, titles, ROOT / "analysis" / "summaries" / "single_domain_bi5_bi8.png", ncols=3,
                suptitle="Isolated 19-cell seed, b_a=5, gamma=0.5, n=3/3, D_i=10, no nucleation (141x141, hex-correct rendering)")
    # inhibitor too, for the colony case
    bi5 = dict(pics)[5]
    save_panels([bi5[300.0][0], bi5[300.0][1]], ["b_i=5, t=300: activator", "b_i=5, t=300: inhibitor"],
                ROOT / "analysis" / "summaries" / "single_domain_bi5_t300_a_i.png", cmaps=["Greens", "Blues"])
    print("pictures written")
