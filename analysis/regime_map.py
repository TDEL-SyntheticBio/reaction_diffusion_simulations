"""
Regime map in the (b_a, b_i) plane from the lone-domain classifier, with the closed-form linear-stability
classes overlaid (CLAUDE.md s10). Panels: Hill pairs 3/3, 4/4 and 10/4 at gamma = 0.5, D_i = 10 (the Fig 1C
value), plus 10/4 at D_i = 20 (the experimental value). b_a in B_A, b_i log-spaced 1..100 (13 values).

Closed-form classes per point (JAPI, hex lattice symbols): "no activated state" (no root with G > alpha-1
on the upper branch), "uniform ON" (tau(0) > 0), "Turing" (Delta(q) < 0 for some q with tau(0) < 0),
"linearly stable" otherwise. Lone-domain classes: dies / spreading / uniform / divides / unstable / stable,
with the growth rate sigma where measured.
Usage: python analysis/regime_map.py run [workers]   ->  analysis/summaries/regime_map.csv
       python analysis/regime_map.py background      ->  analysis/summaries/regime_map_closed_form.csv (fine closed-form grid)
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
    a, i = np.asarray(a, dtype=float), np.asarray(i, dtype=float)      # scalars or arrays
    A = np.where(a > 0, np.abs(a) ** na, 0.0); I = np.where(i > 0, np.abs(i) ** ni, 0.0)
    return A, I


def activated_states(ba, bi, g, na, ni):
    """Positive roots of H = f(A H, I H) with A = b_a, I = b_i/gamma; returns list of (a0, i0, alpha, G, tau0)."""
    A, I = ba, bi / g
    # H grid: uniform on (0, 1) plus a geometric approach to H = 1, because saturated activated states have
    # 1 - H* ~ (1 + I^n_i) / A^n_a, i.e. below 1e-6 at n_a = 10 for b_a >= ~5.3 at low b_i (and 1e-7 at b_i = 0, b_a = 5).
    Hs = np.unique(np.concatenate([np.linspace(1e-6, 1 - 1e-6, 40001), 1 - np.logspace(-6, -16, 2001)]))
    def gfun(h):
        Ah, Ih = hill_terms(A * h, I * h, na, ni)
        return Ah / (1 + Ah + Ih) - h
    gv = gfun(Hs)
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


BG_B_A = [float(f"{v:.3g}") for v in np.logspace(np.log10(2.5), np.log10(20), 41)]
BG_B_I = [float(f"{v:.3g}") for v in np.logspace(0, 2, 81)]
BG_OUT = ROOT / "analysis" / "summaries" / "regime_map_closed_form.csv"


def background():
    rows = []
    for panel, (na, ni, D) in PANELS.items():
        for ba in BG_B_A:
            for bi in BG_B_I:
                cls, alpha, G = closed_form_class(ba, bi, GAMMA, na, ni, D)
                rows.append(dict(panel=panel, b_a=ba, b_i=bi, closed_form=cls, alpha=alpha, G=G))
        print(panel, {c: sum(r["closed_form"] == c for r in rows if r["panel"] == panel) for c in ("no activated state", "uniform ON", "Turing", "linearly stable")}, flush=True)
    with open(BG_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["panel", "b_a", "b_i", "closed_form", "alpha", "G"], lineterminator="\n"); w.writeheader(); w.writerows(rows)


def refresh_closed_form():
    """Rewrite the closed_form / alpha / G columns of regime_map.csv from closed_form_class (the classifier outcomes are untouched)."""
    if not OUT.exists():
        return
    rows = list(csv.DictReader(open(OUT)))
    for r in rows:
        cls, alpha, G = closed_form_class(float(r["b_a"]), float(r["b_i"]), GAMMA, int(r["n_a"]), int(r["n_i"]), float(r["D_i"]))
        r["closed_form"], r["alpha"], r["G"] = cls, alpha, G
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"refreshed closed-form columns of {OUT}", flush=True)


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


# Published sets drawn on the map (open circles with a label): REPORT §10 (Fig 1C line) and §11 (experiment-matched).
PUBLISHED = {
    "n3_3_D10": [(5, 5, "1", "Fig 1C b_i=5: replicating"), (5, 9, "2", "Fig 1C b_i=9: divides"), (5, 12, "3", "Fig 1C b_i=12: arrested")],
    "n10_4_D20": [(5, 15, "4", "3B_big: lone disc arrested, field worms"), (6, 15, "5", "3B_labyrinth: labyrinth"), (5, 20, "6", "3L: arrested"),
                  (5, 25, "7", "2L_JAPI: arrested"), (5, 30, "8", "3B_small: arrested"), (5, 50, "9", "S8D_lowri: arrested")],
}


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = list(csv.DictReader(open(OUT))) if OUT.exists() else []
    bg = list(csv.DictReader(open(BG_OUT))) if BG_OUT.exists() else []
    colors = {"dies": "#bbbbbb", "spreading": "#1b7837", "uniform": "#00441b", "divides": "#7fbc41", "unstable": "#fdb863", "stable": "#b2182b"}
    bg_colors = {"no activated state": "#ffffff", "linearly stable": "#e0e0e0", "uniform ON": "#c6dbef", "Turing": "#fdd0a2"}
    cf_marker = {"no activated state": "", "uniform ON": "U", "Turing": "T", "linearly stable": "S"}
    fig, axes = plt.subplots(1, len(PANELS), figsize=(4.2 * len(PANELS), 4.4), squeeze=False)
    for ax, panel in zip(axes[0], PANELS):
        na, ni, D = PANELS[panel]
        pb = [r for r in bg if r["panel"] == panel]
        if pb:
            bas = sorted({float(r["b_a"]) for r in pb}); bis = sorted({float(r["b_i"]) for r in pb})
            idx = {c: k for k, c in enumerate(bg_colors)}
            Z = np.full((len(bas), len(bis)), np.nan)
            for r in pb:
                Z[bas.index(float(r["b_a"])), bis.index(float(r["b_i"]))] = idx[r["closed_form"]]
            from matplotlib.colors import ListedColormap
            ax.pcolormesh(bis, bas, Z, cmap=ListedColormap(list(bg_colors.values())), vmin=-0.5, vmax=len(bg_colors) - 0.5, shading="nearest")
        for r in rows:
            if r["panel"] != panel:
                continue
            ba, bi = float(r["b_a"]), float(r["b_i"])
            ax.scatter(bi, ba, s=110, c=colors.get(r["outcome"], "k"), marker="s", edgecolors="k", linewidths=0.3)
            ax.text(bi, ba, cf_marker.get(r["closed_form"], "?"), ha="center", va="center", fontsize=6)
        for ba, bi, num, _ in PUBLISHED.get(panel, []):
            ax.scatter(bi, ba, s=260, facecolors="none", edgecolors="k", linewidths=1.2, zorder=5)
            ax.annotate(num, (bi, ba), xytext=(0, 9), textcoords="offset points", fontsize=7, ha="center", zorder=6)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("b_i"); ax.set_ylabel("b_a")
        ax.set_title(f"n_a={na}, n_i={ni}, D_i={D:g}, gamma={GAMMA}", fontsize=9)
    handles = [plt.Line2D([], [], marker="s", ls="", color=c, markeredgecolor="k", label=k) for k, c in colors.items()]
    handles += [plt.Rectangle((0, 0), 1, 1, color=c, label="background: " + k) for k, c in bg_colors.items() if k != "no activated state"]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.045))
    key = "; ".join(f"{num} {label}" for pts in PUBLISHED.values() for _, _, num, label in pts)
    fig.text(0.5, 0.012, "circles, published sets (REPORT §10, §11): " + key, ha="center", fontsize=6.5)
    fig.suptitle("Lone-domain fate (squares; letter = closed-form class of the activated state) over the closed-form background (white = no activated state)", fontsize=9)
    fig.tight_layout(rect=(0, 0.09, 1, 0.95)); fig.savefig(ROOT / "analysis" / "summaries" / "regime_map.png", dpi=160); plt.close(fig)
    print("regime_map.png written")


# Long unperturbed isolated-seed runs for the marginal "unstable" squares (sigma 0.005-0.03) of the map, to tell a slow
# division from a slow hop between neighbouring arrested shapes: classifier seeding (19-cell disc at level r_a, 81x81), to t = 1500,
# area, blob count and per-cell change every 50 t.  Output: analysis/summaries/regime_map_marginal.csv
MARGINAL = [("n4_4_D10", 2.5, 3.16), ("n4_4_D10", 5, 10), ("n4_4_D10", 7, 14.7), ("n4_4_D10", 10, 21.5), ("n4_4_D10", 15, 31.6),      # 31-cell "unstable"
            ("n4_4_D10", 4, 14.7), ("n4_4_D10", 5, 21.5), ("n4_4_D10", 7, 31.6), ("n4_4_D10", 10, 46.4), ("n4_4_D10", 15, 68.1), ("n4_4_D10", 20, 100),   # 7-cell "unstable"
            ("n3_3_D10", 3, 3.16),                                                                      # the n 3/3 "unstable" square
            ("n10_4_D20", 5, 46.4), ("n10_4_D10", 5, 31.6), ("n10_4_D10", 7, 68.1)]                     # 10/4: the "unstable" square and the two largest |sigma| stable ones
MARGINAL_OUT = ROOT / "analysis" / "summaries" / "regime_map_marginal.csv"


def _long_run(pt, n=81, t_end=1500.0, every=50.0):
    from analysis.hexgeom import cell_centres, label_components, neighbor_arrays, rasterize_discs
    from simulation_2D import _vectorized_step
    from analysis.lone_domain_classifier import DT
    panel, ba, bi = pt
    na, ni, D = PANELS[panel]
    p = circuit(ba, bi, GAMMA, na, ni, D)
    nbr = neighbor_arrays(n, n); x, y = cell_centres(n, n); cx, cy = x[n // 2, n // 2], y[n // 2, n // 2]
    thr = 0.3 * p["act_prod_rate"]
    a = np.where(rasterize_discs(n, n, np.array([[cx, cy]]), 2.5), p["act_prod_rate"], 0.0); i = np.zeros_like(a)
    traj, prev = [], None
    for step in range(1, int(round(t_end / DT)) + 1):
        a, i = _vectorized_step(a, i, DT, 1.0, p, "juxtacrine", *nbr)
        if step % int(round(every / DT)) == 0:
            m = a > thr; area = int(m.sum()); nb = label_components(m)[1] if area else 0
            ch = float((np.abs(a - prev[0]).sum() + np.abs(i - prev[1]).sum()) / (2 * n * n)) if prev is not None else np.nan
            traj.append((step * DT, area, nb, ch)); prev = (a.copy(), i.copy())
            if area == 0 or m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any():
                break
    t_last, area_end, blobs_end, ch_end = traj[-1]
    fate = "dies" if area_end == 0 else ("spreads to border" if t_last < t_end else ("divides" if blobs_end > 1 else "one domain"))
    area0 = traj[0][1]
    t_change = next((t for t, A, b, c in traj if A != area0 or b != 1), np.nan)   # first frame whose area or blob count differs from t = 50
    return dict(panel=panel, n_a=na, n_i=ni, D_i=D, b_a=ba, b_i=bi, fate=fate, area_t50=area0, t_first_change=t_change, area_end=area_end,
                blobs_end=blobs_end, change_end=ch_end, t_end=t_last,
                trajectory=" ".join(f"{int(t)}:{A}/{b}" for t, A, b, c in traj))


def marginal(n_workers=4):
    rows = []
    with Pool(n_workers) as pool:
        for r in pool.imap_unordered(_long_run, MARGINAL):
            r["sigma"] = next((float(x["sigma"]) for x in csv.DictReader(open(OUT)) if x["panel"] == r["panel"] and float(x["b_a"]) == r["b_a"] and float(x["b_i"]) == r["b_i"]), np.nan) if OUT.exists() else np.nan
            rows.append(r); print("done", {k: r[k] for k in ("panel", "b_a", "b_i", "sigma", "fate", "area_t50", "t_first_change", "area_end", "blobs_end")}, flush=True)
    rows.sort(key=lambda r: (list(PANELS).index(r["panel"]), r["b_a"], r["b_i"]))
    cols = ["panel", "n_a", "n_i", "D_i", "b_a", "b_i", "sigma", "fate", "area_t50", "t_first_change", "area_end", "blobs_end", "change_end", "t_end", "trajectory"]
    with open(MARGINAL_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"wrote {MARGINAL_OUT}", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "marginal":
        marginal(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif sys.argv[1] == "run":
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif sys.argv[1] == "background":
        background(); refresh_closed_form()
    elif sys.argv[1] == "plot":
        plot()
