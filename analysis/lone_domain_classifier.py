"""
Lone-domain classifier: what does one isolated activated domain do at a given parameter set?

Procedure (one engine step loop around simulation_2D._vectorized_step; nothing in the engine is modified):
  1. Plant a 19-cell disc (radius 2.5) at level r_a = act_prod_rate at the centre of an NxN grid, no nucleation.
  2. Relax for t_relax, recording the active area (cells above 0.3 r_a) every 10 t.
       - area reaches 0 -> "dies"
       - the domain touches the border, or coverage exceeds 0.5 -> "spreading" ("uniform" if coverage > 0.9)
       - area still changing over the last 50 t -> "spreading" (growing or fragmenting colony)
       - otherwise the domain is quasi-static: go to 3
  3. Multiply the activator by (1 + eps * N(0,1)), eps = 1e-4, seeded, and run t_probe more, recording the
     L1 deviation per cell from the pre-perturbation field every 2 t and the blob count every 10 t.
       - blob count > 1 at any time -> "divides"
       - fit log(deviation) over the second half of the probe window: sigma > sigma_min -> "unstable" (will divide),
         else "stable" (arrested domain).
The growth rate sigma is the discriminator at the arrest/replication boundary, where the actual division of a
round-off-seeded domain can take hundreds of time units (REPORT s10, b_i = 9: sigma 0.07, division at t ~ 350).
Usage: python analysis/lone_domain_classifier.py validate     (Fig 1C line, known answers)
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
from analysis.hexgeom import cell_centres, label_components, neighbor_arrays, rasterize_discs  # noqa: E402

DT = 0.01


def circuit(ba, bi, g, na, ni, D):
    return dict(act_half_sat=1.0, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=0.0, act_prod_rate=float(ba),
                inh_prod_rate=float(bi), inh_decay_rate=float(g), act_hill_coeff=na, inh_hill_coeff=ni, inh_diffusion=float(D))


def classify(p, n=81, t_relax=150.0, t_probe=200.0, eps=1e-4, sigma_min=0.01, seed=0, seed_radius=2.5):
    nbr = neighbor_arrays(n, n)
    x, y = cell_centres(n, n)
    cx, cy = x[n // 2, n // 2], y[n // 2, n // 2]
    thr = 0.3 * p["act_prod_rate"]
    a = np.where(rasterize_discs(n, n, np.array([[cx, cy]]), seed_radius), p["act_prod_rate"], 0.0)
    i = np.zeros_like(a)
    areas = []
    res = dict(outcome=None, area=0, blobs=0, sigma=np.nan, t_div=np.nan, coverage=0.0, t_relax=t_relax, t_probe=t_probe)
    for step in range(1, int(round(t_relax / DT)) + 1):
        a, i = _vectorized_step(a, i, DT, 1.0, p, "juxtacrine", *nbr)
        if step % 1000 == 0:
            mask = a > thr
            areas.append(int(mask.sum()))
            if not mask.any():
                res.update(outcome="dies"); return res
            if mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any() or mask.mean() > 0.5:
                res.update(outcome="uniform" if mask.mean() > 0.9 else "spreading", area=int(mask.sum()), coverage=float(mask.mean()),
                           blobs=label_components(mask)[1]); return res
    mask = a > thr
    res.update(area=int(mask.sum()), coverage=float(mask.mean()), blobs=label_components(mask)[1])
    if len(set(areas[-5:])) > 1 or res["blobs"] != 1:
        res.update(outcome="spreading"); return res
    # probe: perturb and watch
    rng = np.random.default_rng(seed)
    a_ref, i_ref = a.copy(), i.copy()
    a = a * (1 + eps * rng.standard_normal(a.shape))
    dev_t, dev = [], []
    for step in range(1, int(round(t_probe / DT)) + 1):
        a, i = _vectorized_step(a, i, DT, 1.0, p, "juxtacrine", *nbr)
        if step % 200 == 0:
            dev_t.append(step * DT); dev.append((np.abs(a - a_ref).sum() + np.abs(i - i_ref).sum()) / (2 * n * n))
        if step % 1000 == 0:
            mask = a > thr
            nb = label_components(mask)[1] if mask.any() else 0
            if nb > 1:
                res.update(outcome="divides", t_div=step * DT, blobs=nb, area=int(mask.sum())); return res
    dev_t, dev = np.array(dev_t), np.array(dev)
    half = dev_t > t_probe / 2
    sigma = float(np.polyfit(dev_t[half], np.log(np.maximum(dev[half], 1e-300)), 1)[0]) if half.sum() > 3 else np.nan
    res.update(sigma=sigma, outcome="unstable" if sigma > sigma_min else "stable", dev_start=float(dev[0]), dev_end=float(dev[-1]))
    return res


def _job(args):
    name, p = args
    r = classify(p); r["name"] = name; r["expected"] = EXPECTED.get(name, "")
    return r


EXPECTED = {"bi5": "spreading", "bi7": "spreading", "bi8": "spreading", "bi9": "unstable", "bi10": "unstable", "bi11": "stable?", "bi12": "stable", "bi14": "stable"}

if __name__ == "__main__":
    if sys.argv[1] == "validate":
        jobs = [(f"bi{bi}", circuit(5, bi, 0.5, 3, 3, 10)) for bi in (5, 7, 8, 9, 10, 11, 12, 14)]
        with Pool(2) as pool:
            out = pool.map(_job, jobs, chunksize=1)
        with open(ROOT / "analysis" / "summaries" / "lone_domain_classifier_validation.csv", "w", newline="") as f:
            cols = ["name", "expected", "outcome", "area", "blobs", "coverage", "sigma", "t_div", "dev_start", "dev_end", "t_relax", "t_probe"]
            w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(out)
        print(f"{'set':<6}{'expected':<10}{'outcome':<11}{'area':>6}{'blobs':>6}{'sigma':>9}{'t_div':>7}{'dev start->end':>22}")
        for r in out:
            print(f"{r['name']:<6}{r['expected']:<10}{r['outcome']:<11}{r['area']:>6}{r['blobs']:>6}{r['sigma']:>9.4f}{r['t_div']:>7.0f}  {r.get('dev_start', float('nan')):.1e} -> {r.get('dev_end', float('nan')):.1e}")
