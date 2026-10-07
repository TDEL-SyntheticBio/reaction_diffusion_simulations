"""
Run a published parameter set under both initiation protocols with several seeds (generalises
run_fig1c_panels.py). Usage: python analysis/run_panels.py <panel> [n_seeds] [t_end] [spike_value_override] [only_bi]
With an override the results go to results/<panel>_spike<value>/ (kick under nucleation = 2*spike_value, CLAUDE.md §6).
Output: results/<panel>/<protocol>_bi<bi>_seed<s>.npz + analysis/summaries/<panel>_runs.csv
"""
from __future__ import annotations

import csv, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from analysis.engine_wrap import circuit, run  # noqa: E402

# (b_a, b_i, gamma, n_a, n_i, D_i), protocols with their spike_value, grid
PANELS = {
    "fig1c": dict(circuits={5: circuit(5, 5, 0.5, 3, 3, 10), 12: circuit(5, 12, 0.5, 3, 3, 10)},
                  protocols={"synchronous": 1.0, "nucleation": 2.0}, grid=(100, 100)),
    # Fig 2L JAPI (b_i=25) and JA (b_i=0); synchronous start needs amplitude 2 to nucleate (session findings §1)
    "fig2L": dict(circuits={25: circuit(5, 25, 0.5, 10, 4, 20), 0: circuit(5, 0, 0.5, 10, 4, 20)},
                  protocols={"nucleation": 2.0, "synchronous": 2.0}, grid=(100, 100)),
}
SAVE_EVERY = 1000   # frame every t = 10


def one(job):
    panel, protocol, bi, seed, t_end, spike, tag = job
    spec = PANELS[panel]
    t0 = time.perf_counter()
    res = run(spec["circuits"][bi], protocol=protocol, seed=seed, t_end=t_end, save_every=SAVE_EVERY,
              spike_value=spike if spike is not None else spec["protocols"][protocol], ny=spec["grid"][0], nx=spec["grid"][1])
    steps, A, I = zip(*res["frames"])
    out = ROOT / "results" / tag
    path = out / f"{protocol}_bi{bi}_seed{seed}.npz"
    np.savez_compressed(path, euler_step=np.array(steps), t=np.array(steps) * res["dt"], a=np.array(A, dtype=np.float32),
                        i=np.array(I, dtype=np.float32), **{k: np.array(v) for k, v in res["settings"].items() if not isinstance(v, str)},
                        protocol=protocol, init_mode=res["settings"]["init_mode"], a_ss=res["a_ss"], i_ss=res["i_ss"])
    return dict(panel=panel, protocol=protocol, bi=bi, seed=seed, steps_used=res["steps_used"], a_ss=res["a_ss"], i_ss=res["i_ss"],
                final_a_max=float(A[-1].max()), wall_s=round(time.perf_counter() - t0, 1), file=path.name)


if __name__ == "__main__":
    panel = sys.argv[1]
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    t_end = float(sys.argv[3]) if len(sys.argv) > 3 else 500.0
    spike = float(sys.argv[4]) if len(sys.argv) > 4 else None
    only_bi = int(sys.argv[5]) if len(sys.argv) > 5 else None
    spec = PANELS[panel]
    tag = panel if spike is None else f"{panel}_spike{spike:g}"
    (ROOT / "results" / tag).mkdir(parents=True, exist_ok=True)
    protocols = list(spec["protocols"]) if spike is None else ["nucleation"]
    jobs = [(panel, prot, bi, seed, t_end, spike, tag) for prot in protocols for bi in spec["circuits"] for seed in range(1, n_seeds + 1)
            if not (bi == 0 and prot == "synchronous") and (only_bi is None or bi == only_bi)]   # JA control only under the published nucleation protocol
    with Pool(3) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    summary = ROOT / "analysis" / "summaries" / f"{tag}_runs.csv"
    with open(summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
