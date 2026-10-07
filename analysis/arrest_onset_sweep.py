"""
Does an ISOLATED domain arrest at the Fig 1C spot parameters? Sweep b_i at b_a=5, gamma=0.5, n=3/3, D_i=10
(the handoff's b_i sweep line), single disc seed, no nucleation, up to t=600; report arrested size or growth
and the aspect ratio of the largest component. NOTE: "arrested" here means 4 equal half-max-area samples
1000 steps apart; isolated_seed_fate.py showed that b_i = 8-10 pass this test and divide later, so use that
script's time series for the fate and this one only for the halo decay length of truly stable domains.
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
from analysis.hexgeom import cell_centres, label_components, equivalent_radius  # noqa: E402


def shape(mask):
    """Principal-moment aspect ratio of the LARGEST component (1 = disc; nan below 3 cells) and the component count."""
    ny, nx = mask.shape
    x, y = cell_centres(ny, nx)
    labels, ncomp = label_components(mask)
    if ncomp == 0:
        return np.nan, 0
    areas = np.bincount(labels[labels >= 0], minlength=ncomp)
    big = labels == int(np.argmax(areas))
    if big.sum() < 3:
        return np.nan, ncomp
    cov = np.cov(np.vstack([x[big], y[big]]))
    ev = np.sort(np.linalg.eigvalsh(cov))
    return float(np.sqrt(ev[1] / max(ev[0], 1e-12))), ncomp


def one(bi):
    p = params(5, bi, 0.5, 3, 3, 10)
    a, i, step, hist, c = run_single_domain_2d(p, ny=141, nx=141, max_steps=60000, check_every=1000)
    if a.max() < 0.1:
        return dict(bi=bi, outcome="died", touches_edge=False, t_end=(step + 1) / 100, area=0, R=np.nan, aspect=np.nan, ncomp=0, lam=np.nan, area_hist=str(hist[-5:]))
    mask = a > 0.5 * a.max()
    arrested = len(set(hist[-4:])) == 1
    asp, ncomp = shape(mask)
    touches_edge = mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any()
    try:
        lam = fit_exterior_2d(a, i, c, r_max=60)["lam_k0"] if arrested and not touches_edge else np.nan
    except Exception:
        lam = np.nan
    outcome = "filled/edge" if touches_edge else ("arrested" if arrested else "growing")
    return dict(bi=bi, outcome=outcome, touches_edge=touches_edge, t_end=(step + 1) / 100,
                area=int(mask.sum()), R=round(equivalent_radius(mask.sum()), 2), aspect=round(asp, 2) if np.isfinite(asp) else np.nan, ncomp=ncomp,
                lam=round(float(lam), 2) if np.isfinite(lam) else np.nan, area_hist=str(hist[-5:]))


if __name__ == "__main__":
    bis = [4, 5, 6, 7, 8, 10, 12, 14]
    with Pool(1) as pool:
        rows = pool.map(one, bis, chunksize=1)
    out = ROOT / "analysis" / "summaries" / "arrest_onset.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
