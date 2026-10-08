"""
Regime map in the (b_a, b_i) plane from the lone-domain classifier, with the closed-form linear-stability
classes overlaid (CLAUDE.md s10). Panels: Hill pairs 3/3, 4/4 and 10/4 at gamma = 0.5, D_i = 10 (the Fig 1C
value), plus 10/4 at D_i = 20 (the experimental value). b_a in B_A, b_i log-spaced 1..100 (13 values).

Closed-form classes per point (JAPI, hex lattice symbols): "no activated state" (no root with G > alpha-1
on the upper branch), "uniform ON" (tau(0) > 0), "Turing" (Delta(q) < 0 for some q with tau(0) < 0),
"linearly stable" otherwise. Lone-domain classes: dies / spreading / uniform / divides / unstable / stable,
with the growth rate sigma where measured.
Usage: python analysis/regime_map.py run [workers]   ->  analysis/summaries/regime_map.csv
       python analysis/regime_map.py plot            ->  analysis/summaries/regime_map.png
"""
from __future__ import annotations

import csv, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "2D_simulations"))
from analysis.lone_domain_classifier import circuit, classify  # noqa: E402

PANELS = {"n3_3_D10": (3, 3, 10.0), "n4_4_D10": (4, 4, 10.0), "n10_4_D10": (10, 4, 10.0), "n10_4_D20": (10, 4, 20.0)}
B_A = [2.5, 3, 4, 5, 7, 10, 15, 20]
B_I = [float(f"{v:.3g}") for v in np.logspace(0, 2, 13)]
GAMMA = 0.5
OUT = ROOT / "analysis" / "summaries" / "regime_map.csv"


# ---- closed forms (CLAUDE.md s10) --------------------------------------------------------------------
def hill_terms(a, i, na, ni):
    A = a ** na if a > 0 else 0.0; I = i ** ni if i > 0 else 0.0
    return A, I


def activated_states(ba, bi, g, na, ni):
    """Positive roots of H = f(A H, I H) with A = b_a, I = b_i/gamma; returns list of (a0, i0, alpha, G, tau0)."""
    A, I = ba, bi / g
    Hs = np.linspace(1e-6, 1 - 1e-6, 40001)
    def gfun(h):
        Ah, Ih = hill_terms(A * h, I * h, na, ni)
        return Ah / (1 + Ah + Ih) - h
    gv = np.array([gfun(h) for h in Hs])
    roots = []
    for k in np.nonzero(np.sign(gv[:-1]) != np.sign(gv[1:]))[0]:
        lo, hi = Hs[k], Hs[k + 1]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if gfun(lo) * gfun(mid) <= 0:
                hi = mid
            else:
                lo = mid
        h = 0.5 * (lo + hi); a0, i0 = A * h, I * h
        Ah, Ih = hill_terms(a0, i0, na, ni); den = 1 + Ah + Ih
        P, Q = Ah / den, Ih / den
        alpha, G = na * (1 - P), ni * Q
        roots.append((a0, i0, alpha, G, alpha - 1 - g * (1 + G)))
    return roots


def closed_form_class(ba, bi, g, na, ni, D):
    roots = [r for r in activated_states(ba, bi, g, na, ni) if r[3] > r[2] - 1]    # upper branch: G > alpha - 1
    if not roots:
        return "no activated state", np.nan, np.nan
    a0, i0, alpha, G, tau0 = max(roots, key=lambda r: r[0])
    if tau0 > 0:
        return "uniform ON", alpha, G
    lam2 = D / g
    K = np.linspace(-0.5, 1, 3001)                       # hex kernel symbol; Lambda_hex = 6 (1 - K)
    delta = (1 + lam2 * 6 * (1 - K)) * (1 - alpha * K) + G
    return ("Turing" if delta.min() < 0 else "linearly stable"), alpha, G


def _job(args):
    panel, ba, bi = args
    na, ni, D = PANELS[panel]
    t0 = time.perf_counter()
    r = classify(circuit(ba, bi, GAMMA, na, ni, D))
    cls, alpha, G = closed_form_class(ba, bi, GAMMA, na, ni, D)
    r.update(panel=panel, n_a=na, n_i=ni, D_i=D, b_a=ba, b_i=bi, closed_form=cls, alpha=alpha, G=G, wall_s=round(time.perf_counter() - t0, 1))
    return r


def run(n_workers=4):
    jobs = [(panel, ba, bi) for panel in PANELS for ba in B_A for bi in B_I]
    print(f"{len(jobs)} classifications on {n_workers} workers", flush=True)
    rows = []
    with Pool(n_workers) as pool:
        for r in pool.imap_unordered(_job, jobs, chunksize=1):
            rows.append(r); print("done", {k: r[k] for k in ("panel", "b_a", "b_i", "outcome", "sigma", "closed_form", "wall_s")}, flush=True)
    rows.sort(key=lambda r: (list(PANELS).index(r["panel"]), r["b_a"], r["b_i"]))
    cols = ["panel", "n_a", "n_i", "D_i", "b_a", "b_i", "closed_form", "alpha", "G", "outcome", "sigma", "rising", "area", "blobs", "coverage", "t_div",
            "dev_start", "dev_mid", "dev_end", "t_relax", "t_probe", "wall_s"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT}", flush=True)


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = list(csv.DictReader(open(OUT)))
    colors = {"dies": "#bbbbbb", "spreading": "#1b7837", "uniform": "#00441b", "divides": "#7fbc41", "unstable": "#fdb863", "stable": "#b2182b"}
    cf_marker = {"no activated state": "", "uniform ON": "U", "Turing": "T", "linearly stable": "S"}
    fig, axes = plt.subplots(1, len(PANELS), figsize=(4.2 * len(PANELS), 4.4), squeeze=False)
    for ax, panel in zip(axes[0], PANELS):
        na, ni, D = PANELS[panel]
        for r in rows:
            if r["panel"] != panel:
                continue
            ba, bi = float(r["b_a"]), float(r["b_i"])
            ax.scatter(bi, ba, s=110, c=colors.get(r["outcome"], "k"), marker="s", edgecolors="k", linewidths=0.3)
            ax.text(bi, ba, cf_marker.get(r["closed_form"], "?"), ha="center", va="center", fontsize=6)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("b_i"); ax.set_ylabel("b_a")
        ax.set_title(f"n_a={na}, n_i={ni}, D_i={D:g}, gamma={GAMMA}", fontsize=9)
    handles = [plt.Line2D([], [], marker="s", ls="", color=c, markeredgecolor="k", label=k) for k, c in colors.items()]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Lone-domain fate (colour) and closed-form linear class (letter: T Turing, U uniform ON, S linearly stable, blank no activated state)", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95)); fig.savefig(ROOT / "analysis" / "summaries" / "regime_map.png", dpi=160); plt.close(fig)
    print("regime_map.png written")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif sys.argv[1] == "plot":
        plot()
