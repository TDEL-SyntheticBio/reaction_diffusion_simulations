"""
D_i sweep at the 2L set (b_a 5, b_i 25, gamma 0.5, n 10/4): do the pattern's lengths scale with the inhibitor range
lambda_i = sqrt(1.5 D_i / gamma) = 3.87, 5.48, 7.75 cells at D_i = 5, 10, 20?  Sets bi25_D5, bi25_D10 (this sweep) and
bi25 (the ladder batch, D_i 20), all 200x200 to t = 600 under the ladder driver; conditions seed1, seed36, random, nucl.

Lengths at t = 600, in cells and divided by lambda_i, per set and condition (range over the three replicates):
  R_lone      equivalent radius of the lone seed's arrested disc, sqrt(area * sqrt(3)/2 / pi)
  nn_mean     mean nearest-neighbour spacing of real-domain centroids (ladder_frames.csv)
  hole        g(r) hole radius, per run (ladder_gr.csv) and pooled over the three replicates (PooledG on the t = 600 real-domain
              centroids under the same persistence mask as the per-run rows: a > 0.3 r_a at t = 600 and at t = 590)
  peak1_r     first g(r) peak position, per run and pooled (NaN where the 3-sigma / 20-pair gate is not passed)
  R_eq_mean   mean over real domains of the per-domain equivalent radius (from the areas lists); R_eq_median likewise
  fft_wavelength  hex-correct power-spectrum peak wavelength of the activator field (ladder_fft.csv); its relative FWHM is kept as a non-length extra
Grid-size control: the step-1 set 3H_lowDi (same parameters at D_i 5 on 141x141, seed_crowding_exp_frames.csv) is written as set
3H_lowDi_141 for the lone-seed radius only; it is not part of the sweep.
Outputs: analysis/summaries/di_sweep_lengths.csv (long: one row per set, condition, quantity, with cells, lambda and ratio),
         analysis/summaries/di_sweep_collapse.png (ratios versus D_i), analysis/summaries/di_sweep_gr.png (pooled g against r / lambda_i).
Usage: python analysis/di_sweep_collapse.py
"""
from __future__ import annotations
import csv, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from analysis.hexgeom import cell_centres, component_table, window  # noqa: E402
from analysis.pointstats import PooledG  # noqa: E402

SUMM = ROOT / "analysis" / "summaries"; RES = ROOT / "results" / "ladder"
GAMMA = 0.5
SETS = {"bi25_D5": 5.0, "bi25_D10": 10.0, "bi25": 20.0}
CONDS = ["seed1", "seed36", "random", "nucl"]
LAM = {s: float(np.sqrt(1.5 * D / GAMMA)) for s, D in SETS.items()}
CELL_AREA = np.sqrt(3) / 2
NOT_LENGTH = ("n_real", "pooled_pairs", "peak1_g_pooled", "fft_fwhm_rel")


def req(area):
    return np.sqrt(np.asarray(area, dtype=float) * CELL_AREA / np.pi)


def rng_str(v, fmt="{:.2f}"):
    v = np.asarray([x for x in v if np.isfinite(x)], dtype=float)
    if v.size == 0:
        return "n/a"
    return fmt.format(v.min()) if v.min() == v.max() else f"{fmt.format(v.min())}–{fmt.format(v.max())}"


def main():
    frames = pd.read_csv(SUMM / "ladder_frames.csv"); frames = frames[(frames.channel == "a") & (np.abs(frames.t - 600.01) < 0.05)]
    gr = pd.read_csv(SUMM / "ladder_gr.csv"); gr = gr[(gr.channel == "a") & (gr.t == 600)]
    fft = pd.read_csv(SUMM / "ladder_fft.csv"); fft = fft[(fft.channel == "a") & (fft.t == 600)]
    rows, curves = [], {}
    win = window(200, 200); lattice_xy = np.column_stack([c.ravel() for c in cell_centres(200, 200)])
    for s, D in SETS.items():
        lam = LAM[s]
        for cond in CONDS:
            f = frames[(frames.set == s) & (frames.condition == cond)]
            if f.empty:
                continue
            g = gr[(gr.set == s) & (gr.condition == cond)]
            ff = fft[(fft.set == s) & (fft.condition == cond)]
            vals = {}
            if cond == "seed1":
                vals["R_lone"] = [float(req(a)) for a in f.largest]
            else:
                vals["nn_mean"] = list(f.nn_mean)
                vals["hole_run"] = list(g.hole); vals["peak1_r_run"] = list(g.peak1_r)
                re_mean, re_median = [], []
                for areas in f.areas:
                    ar = np.array([int(v) for v in str(areas).split()]); ar = ar[ar >= 6]
                    re_mean.append(float(req(ar).mean())); re_median.append(float(np.median(req(ar))))
                vals["R_eq_mean"] = re_mean; vals["R_eq_median"] = re_median
                vals["n_real"] = list(f.blobs_real.astype(float))
                vals["fft_wavelength"] = list(ff.fft_peak_wavelength); vals["fft_fwhm_rel"] = list(ff.fft_peak_fwhm_rel)
                # pooled g(r) over the three replicates on the stored t = 600 activator snapshots
                pg = None
                for r in sorted(f.replicate):
                    d = np.load(RES / f"{s}_{cond}_rep{int(r)}.npz"); thr = float(d["thr_a"]) if "thr_a" in d else 1.5
                    ft = d["frame_t"]; k600 = int(np.argmin(np.abs(ft - 600.01))); k590 = int(np.argmin(np.abs(ft - 590.01)))
                    assert abs(ft[k600] - 600.01) < 0.05 and abs(ft[k590] - 590.01) < 0.05
                    mask = (d["a_frames"][k600] > thr) & (d["a_frames"][k590] > thr)  # persistence mask, as ladder_frames.csv
                    tab = component_table(mask, min_size=6); xy = np.column_stack([tab["cx"], tab["cy"]])
                    if pg is None:
                        pg = PooledG(r_max=40.0, dr=0.5, n_null=200, seed=1, lattice_xy=lattice_xy if np.median(tab["area"]) < 10 else None)
                    pg.add(xy, win)
                summ = pg.summary(n_boot=300)
                vals["hole_pooled"] = [summ["hole"]]; vals["peak1_r_pooled"] = [summ["peak1_r"]]
                vals["peak1_g_pooled"] = [summ["peak1_g"]]; vals["pooled_pairs"] = [summ["total_pairs"]]
                vals["hole_pooled_sd"] = [summ["hole_sd"]]; vals["peak1_r_pooled_sd"] = [summ["peak1_r_sd"]]
                curves[(s, cond)] = (pg.r, pg.g(), summ)
            for q, v in vals.items():
                v = np.asarray(v, dtype=float); fin = v[np.isfinite(v)]
                is_len = q not in NOT_LENGTH
                rows.append(dict(set=s, D_i=D, lambda_i=round(lam, 3), condition=cond, quantity=q, n=int(fin.size),
                                 cells_min=float(fin.min()) if fin.size else np.nan, cells_max=float(fin.max()) if fin.size else np.nan,
                                 cells_mean=float(fin.mean()) if fin.size else np.nan,
                                 ratio_min=float(fin.min() / lam) if (fin.size and is_len) else np.nan,
                                 ratio_max=float(fin.max() / lam) if (fin.size and is_len) else np.nan,
                                 ratio_mean=float(fin.mean() / lam) if (fin.size and is_len) else np.nan,
                                 values=" ".join(f"{x:.4g}" for x in v)))
    # grid-size control: 3H_lowDi lone seed (D_i 5, 141x141, step 1) at t = 600
    ctrl = pd.read_csv(SUMM / "seed_crowding_exp_frames.csv")
    ctrl = ctrl[(ctrl.set == "3H_lowDi") & (ctrl.condition == "seed1") & (np.abs(ctrl.t - 600.01) < 0.05)]
    if not ctrl.empty:
        lam = LAM["bi25_D5"]; v = np.asarray([float(req(a)) for a in ctrl.largest_blob])
        rows.append(dict(set="3H_lowDi_141", D_i=5.0, lambda_i=round(lam, 3), condition="seed1", quantity="R_lone", n=int(v.size),
                         cells_min=float(v.min()), cells_max=float(v.max()), cells_mean=float(v.mean()),
                         ratio_min=float(v.min() / lam), ratio_max=float(v.max() / lam), ratio_mean=float(v.mean() / lam),
                         values=" ".join(f"{x:.4g}" for x in v)))
    out = pd.DataFrame(rows); out.to_csv(SUMM / "di_sweep_lengths.csv", index=False, lineterminator="\n")
    # ---- printed table: cells and ratio, range over replicates ----
    print(f"lambda_i: " + ", ".join(f"D_i {D:g} -> {LAM[s]:.2f}" for s, D in SETS.items()))
    for cond in CONDS:
        print(f"\n== {cond} ==")
        qs = ["R_lone"] if cond == "seed1" else ["nn_mean", "hole_run", "hole_pooled", "peak1_r_run", "peak1_r_pooled", "R_eq_mean", "R_eq_median",
                                                  "fft_wavelength", "fft_fwhm_rel", "n_real"]
        print(f"{'quantity':<16}" + "".join(f"{'D_i ' + str(int(D)) + ' cells':>22}{'/lambda':>14}" for D in SETS.values()))
        for q in qs:
            line = f"{q:<16}"
            for s in SETS:
                r = out[(out.set == s) & (out.condition == cond) & (out.quantity == q)]
                if r.empty:
                    line += f"{'-':>22}{'-':>14}"; continue
                v = np.array([float(x) for x in r["values"].iloc[0].split()])
                line += f"{rng_str(v):>22}{(rng_str(v / LAM[s]) if q not in NOT_LENGTH else ''):>14}"
            print(line)
        if cond == "seed1":
            c = out[out.set == "3H_lowDi_141"]
            if not c.empty:
                print(f"grid-size control 3H_lowDi (D_i 5, 141x141): R_lone {c.cells_mean.iloc[0]:.2f} cells, /lambda {c.ratio_mean.iloc[0]:.2f}")
    plots(out, curves)


def plots(out, curves):
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {"seed36": "#eb6834", "random": "#eda100", "nucl": "#e87ba4", "seed1": "#2a78d6"}
    label = {"seed36": "36 seeds", "random": "random U(0,2)", "nucl": "nucleation", "seed1": "lone seed"}
    Ds = [5.0, 10.0, 20.0]
    quantities = [("R_lone", "lone-seed disc radius R_eq"), ("nn_mean", "mean NN spacing"), ("hole_pooled", "g(r) hole radius (pooled)"),
                  ("peak1_r_pooled", "first g(r) peak (pooled)"), ("R_eq_mean", "mean domain radius R_eq"), ("R_eq_median", "median domain radius R_eq")]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2)); fig.patch.set_facecolor("#fcfcfb")
    for ax, (q, title) in zip(axes.ravel(), quantities):
        ax.set_facecolor("#fcfcfb"); ax.grid(True, color="#e6e5e1", linewidth=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for cond in (["seed1"] if q == "R_lone" else ["seed36", "random", "nucl"]):
            r = out[(out.condition == cond) & (out.quantity == q) & (out.set != "3H_lowDi_141")].sort_values("D_i")
            if r.empty or not np.isfinite(r.ratio_mean).any():
                continue
            x = list(r.D_i); y = list(r.ratio_mean); lo = list(r.ratio_min); hi = list(r.ratio_max)
            ax.plot(x, y, "-o", color=col[cond], linewidth=2, markersize=5, markerfacecolor="#fcfcfb", markeredgewidth=1.6, label=label[cond])
            ax.fill_between(x, lo, hi, color=col[cond], alpha=0.15, linewidth=0)
        if q == "R_lone":
            c = out[(out.set == "3H_lowDi_141") & (out.quantity == q)]
            if not c.empty:
                ax.plot(c.D_i, c.ratio_mean, marker="x", linestyle="none", color="#52514e", markersize=8, markeredgewidth=1.6,
                        label="3H_lowDi, 141x141 (grid-size control)")
        ax.set_xscale("log"); ax.set_xticks(Ds); ax.set_xticklabels([f"{d:g}" for d in Ds]); ax.minorticks_off()
        ax.set_xlabel("D_i", fontsize=9, color="#52514e"); ax.set_ylabel("length / lambda_i", fontsize=9, color="#52514e")
        ax.set_title(title, fontsize=10, loc="left"); ax.set_ylim(bottom=0); ax.legend(fontsize=8, frameon=False)
    fig.suptitle("D_i sweep at b_a 5, b_i 25, n 10/4, gamma 0.5 (200x200, t = 600): lengths in units of lambda_i = sqrt(1.5 D_i/gamma); bands = replicate range", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(SUMM / "di_sweep_collapse.png", dpi=160, facecolor=fig.get_facecolor()); plt.close(fig)
    # pooled g(r) against r / lambda
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2)); fig.patch.set_facecolor("#fcfcfb")
    dcol = {5.0: "#2a78d6", 10.0: "#1baf7a", 20.0: "#e87ba4"}
    for row, scaled in enumerate((False, True)):
        for ax, cond in zip(axes[row], ["seed36", "random", "nucl"]):
            ax.set_facecolor("#fcfcfb"); ax.grid(True, color="#e6e5e1", linewidth=0.6); ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for s, D in SETS.items():
                if (s, cond) not in curves:
                    continue
                r, g, summ = curves[(s, cond)]
                ax.plot(r / LAM[s] if scaled else r, g, color=dcol[D], linewidth=1.8,
                        label=f"D_i {D:g} (lambda {LAM[s]:.2f}; {int(summ['total_pairs'])} pairs)")
            gmax = max([float(np.nanmax(curves[(s, cond)][1])) for s in SETS if (s, cond) in curves] + [3.0])
            ax.axhline(1.0, color="#c3c2b7", linewidth=1); ax.set_ylim(0, 1.06 * gmax)
            ax.set_xlim(0, 40 / LAM["bi25_D5"] if scaled else 40)
            ax.set_xlabel("r / lambda_i" if scaled else "r (cells)", fontsize=9, color="#52514e")
            ax.set_ylabel("g(r), pooled over 3 replicates", fontsize=9, color="#52514e")
            ax.set_title(label[cond] + (" — r / lambda_i" if scaled else " — r in cells"), fontsize=10, loc="left"); ax.legend(fontsize=7.5, frameon=False)
    fig.suptitle("Pooled pair correlation of activator real domains at t = 600 (CSR null): top in cells, bottom in units of lambda_i", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(SUMM / "di_sweep_gr.png", dpi=160, facecolor=fig.get_facecolor()); plt.close(fig)
    print("di_sweep_collapse.png, di_sweep_gr.png written")


if __name__ == "__main__":
    main()
