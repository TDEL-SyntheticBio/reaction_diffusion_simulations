"""
Is the Turing regime reachable at the measured Hill pair (n_a 10, n_i 4)? The closed forms (CLAUDE.md s10) put a
narrow finite-q instability of the activated state just below the saddle-node; here the engine is started at the
solver's activated homogeneous state with +-5% noise (init_mode random_tight) and the growth of spatial structure is
recorded. Output: analysis/summaries/turing_sliver_check.csv
"""
from __future__ import annotations

import contextlib, csv, io, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "2D_simulations"))
from simulation_2D import run_coupled_hex  # noqa: E402
from finding_steady_states import fast_stable_steady_state  # noqa: E402
from analysis.regime_map import closed_form_class  # noqa: E402
from analysis.lone_domain_classifier import circuit  # noqa: E402

N, T_END = 60, 80.0
CASES = [(5, bi, D) for D in (20, 10) for bi in (13.0, 13.5, 13.8, 13.9, 14.0, 14.1, 14.15)] + [(10, bi, 20) for bi in (40, 44, 46, 47, 48)]

if __name__ == "__main__":
    rows = []
    for ba, bi, D in CASES:
        p = circuit(ba, bi, 0.5, 10, 4, D)
        cf, alpha, G = closed_form_class(ba, bi, 0.5, 10, 4, D)
        a_ss, i_ss, _ = fast_stable_steady_state(p)
        if not (a_ss > 0 and i_ss > 0):
            rows.append(dict(b_a=ba, b_i=bi, D_i=D, closed_form=cf, alpha=alpha, G=G, a_ss=a_ss, engine="no activated state (solver)")); continue
        np.random.seed(1)
        with contextlib.redirect_stdout(io.StringIO()):
            A, I, *_ = run_coupled_hex(N, N, int(T_END / 0.01) + 1, 0.01, 1.0, p, 0.0, 10**9, init_mode="random_tight", spike_value=1.0, save_every=2000, nucleation_rate=0.0)
        rel = [float(np.std(a) / max(np.mean(a), 1e-12)) for a in A]
        grows = rel[-1] > 10 * rel[1]
        rows.append(dict(b_a=ba, b_i=bi, D_i=D, closed_form=cf, alpha=round(alpha, 3), G=round(G, 3), a_ss=a_ss,
                         spread_t0=round(rel[1], 3), spread_t40=round(rel[len(rel) // 2], 3), spread_end=round(rel[-1], 3),
                         final_min=round(float(A[-1].min()), 2), final_max=round(float(A[-1].max()), 2), engine="pattern grows" if grows else "perturbation decays"))
    with open(ROOT / "analysis" / "summaries" / "turing_sliver_check.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[-1].keys()), lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
