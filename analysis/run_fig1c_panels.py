"""
Run the Fig 1C spot panels (b_a=5, gamma=0.5, n=3/3, D_i=10; b_i=5 'regular', b_i=12 'irregular') under
both published initiation protocols, several seeds each, to t=500, saving a frame every t=10.
Output: results/fig1c_panels/<protocol>_bi<bi>_seed<s>.npz (gitignored) + analysis/summaries/fig1c_runs.csv
"""
from __future__ import annotations

import csv, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from analysis.engine_wrap import circuit, run  # noqa: E402

OUT = ROOT / "results" / "fig1c_panels"
SUMMARY = ROOT / "analysis" / "summaries" / "fig1c_runs.csv"
T_END, SAVE_EVERY = 500.0, 1000     # frame every t = 10


def one(job):
    protocol, bi, seed = job
    t0 = time.perf_counter()
    res = run(circuit(5, bi, 0.5, 3, 3, 10), protocol=protocol, seed=seed, t_end=T_END, save_every=SAVE_EVERY)
    steps, A, I = zip(*res["frames"])
    path = OUT / f"{protocol}_bi{bi}_seed{seed}.npz"
    np.savez_compressed(path, euler_step=np.array(steps), t=np.array(steps) * res["dt"], a=np.array(A, dtype=np.float32),
                        i=np.array(I, dtype=np.float32), **{k: np.array(v) for k, v in res["settings"].items() if not isinstance(v, str)},
                        protocol=protocol, init_mode=res["settings"]["init_mode"], a_ss=res["a_ss"], i_ss=res["i_ss"])
    return dict(protocol=protocol, bi=bi, seed=seed, steps_used=res["steps_used"], a_ss=res["a_ss"], i_ss=res["i_ss"],
                final_a_max=float(A[-1].max()), wall_s=round(time.perf_counter() - t0, 1), file=path.name)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(prot, bi, seed) for prot in ("synchronous", "nucleation") for bi in (5, 12) for seed in (1, 2, 3)]
    with Pool(3) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    with open(SUMMARY, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
    print(f"wrote {SUMMARY}")
