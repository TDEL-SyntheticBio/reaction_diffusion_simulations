"""
Does an ISOLATED domain arrest at the Fig 1C spot parameters? Sweep b_i at b_a=5, gamma=0.5, n=3/3, D_i=10
(the handoff's b_i sweep line), single disc seed, no nucleation, up to t=600; report arrested size or growth
and a shape descriptor (ratio of principal second moments; 1 = disc, large = finger/labyrinth).
Output: analysis/summaries/arrest_onset.csv
"""
from __future__ import annotations

import csv, sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from analysis.single_domain_halo import params, run_single_domain_2d, fit_exterior_2d  # noqa: E402
from analysis.hexgeom import cell_centres, label_components  # noqa: E402


def shape(mask):
    ny, nx = mask.shape
    x, y = cell_centres(ny, nx)
    xs, ys = x[mask], y[mask]
    cov = np.cov(np.vstack([xs, ys]))
    ev = np.sort(np.linalg.eigvalsh(cov))
    _, ncomp = label_components(mask)
    return float(np.sqrt(ev[1] / max(ev[0], 1e-12))), ncomp


def one(bi):
    p = params(5, bi, 0.5, 3, 3, 10)
    a, i, step, hist, c = run_single_domain_2d(p, ny=141, nx=141, max_steps=60000, check_every=1000)
    if a.max() < 0.1:
        return dict(bi=bi, outcome="died", t_end=(step + 1) / 100, area=0, R=np.nan, aspect=np.nan, ncomp=0, lam=np.nan, area_hist=str(hist[-5:]))
    mask = a > 0.5 * a.max()
    arrested = len(set(hist[-4:])) == 1
    asp, ncomp = shape(mask)
    touches_edge = mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any()
    try:
        lam = fit_exterior_2d(a, i, c, r_max=60)["lam_k0"] if arrested and not touches_edge else np.nan
    except Exception:
        lam = np.nan
    return dict(bi=bi, outcome="arrested" if arrested else ("filled/edge" if touches_edge else "growing"), t_end=(step + 1) / 100,
                area=int(mask.sum()), R=round(float(np.sqrt(mask.sum() / np.pi)), 2), aspect=round(asp, 2), ncomp=ncomp,
                lam=round(float(lam), 2) if np.isfinite(lam) else np.nan, area_hist=str(hist[-5:]))


if __name__ == "__main__":
    bis = [4, 5, 6, 7, 8, 10, 12, 14]
    with Pool(1) as pool:
        rows = pool.map(one, bis, chunksize=1)
    out = ROOT / "analysis" / "summaries" / "arrest_onset.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
