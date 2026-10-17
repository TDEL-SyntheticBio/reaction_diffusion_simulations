"""
Picture sheets of the 2D initiation sweep from the saved frames (results/sweep2d/frames/*.npz; nothing is rerun).

One PNG per (start, D_i) for replicate 0: rows = the 9 Hill pairs (2/2 ... 10/10), columns = the 8 edge positions
(0.1 ... 3 x edge). Two sets: the activator at t = 50 and at the day-4 equivalent (t = 17.8, the fourth saved frame).
Colour scale fixed at 0 .. b_a = 3 in every panel. Panel title: label (at t = 50, from the row file), n_spots, n_stable
(at t = 50 from the row file for the t = 50 set; at t = 18, the nearest recorded time, from the run's time series for the
day-4 set).

Hex-correct, fast raster (CLAUDE.md section 2, even-r): each cell is 2 pixels wide, even rows are shifted right by one
pixel, and rows are drawn with aspect sqrt(3)/2 (row pitch sqrt(3) pixels for a 2-pixel cell). visualize_2D.py is not used.

Usage: python analysis/sweep2d_sheets.py            -> analysis/summaries/sweep2d_sheets/{t50,day4}_<start>_D<D_i>.png
"""
from __future__ import annotations

import json, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FRAMES = ROOT / "results" / "sweep2d" / "frames"; ROWS = ROOT / "results" / "sweep2d" / "rows"
OUT = ROOT / "analysis" / "summaries" / "sweep2d_sheets"
BA = 3.0
HILL = [(2, 2), (2, 4), (2, 10), (4, 2), (4, 4), (4, 10), (10, 2), (10, 4), (10, 10)]
F_EDGE = [0.1, 0.2, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0]
D_I = [5, 10, 20]
STARTS = ["seed1", "seeds", "rf1", "rf2", "nuc01", "nuc001"]
SETS = {"t50": (5, 50.0, "t = 50"), "day4": (3, 17.8, "day-4 equivalent, t = 17.8")}


def hex_raster(field):
    """(ny, 2 nx + 1) image: cell (r, c) covers pixels [2c + s, 2c + s + 2) with s = 1 on even rows (even-r offset)."""
    ny, nx = field.shape
    img = np.full((ny, 2 * nx + 1), np.nan)
    for s, rows in ((1, slice(0, ny, 2)), (0, slice(1, ny, 2))):
        f = field[rows]
        img[rows, s:s + 2 * nx:2] = f
        img[rows, s + 1:s + 1 + 2 * nx:2] = f
    return img


def sheet(start, D, key):
    idx, t_lab, t_text = SETS[key]
    fig, axes = plt.subplots(len(HILL), len(F_EDGE), figsize=(2.05 * len(F_EDGE), 1.95 * len(HILL) + 0.6), squeeze=False)
    fig.patch.set_facecolor("white")
    missing = 0
    for r, (na, ni) in enumerate(HILL):
        for c, f in enumerate(F_EDGE):
            ax = axes[r, c]; ax.set_xticks([]); ax.set_yticks([])
            name = f"na{na:02d}_ni{ni:02d}_f{c}_D{D:02d}_{start}_r0"
            fp = FRAMES / f"{name}.npz"; rp = ROWS / f"{name}.json"
            if not fp.exists() or not rp.exists():
                ax.set_title("missing", fontsize=6); ax.set_facecolor("#eeeeee"); missing += 1; continue
            d = np.load(fp); row = json.load(open(rp))
            assert abs(row["f_edge"] - f) < 1e-9 and abs(row["D_i"] - D) < 1e-9 and row["start"] == start and row["rep"] == 0, name
            a = d["a"][idx].astype(np.float32)
            assert abs(float(d["saved_t"][idx]) - t_lab) < 0.05, (name, d["saved_t"][idx])
            img = hex_raster(a)
            ny, nx = a.shape
            ax.imshow(img, cmap="Greens", vmin=0.0, vmax=BA, origin="lower", interpolation="nearest",
                      extent=(0, img.shape[1], 0, ny * np.sqrt(3)), aspect="equal")
            if key == "t50":
                ns, nst = row["n_spots"], row["n_stable"]
            else:
                k = int(np.argmin(np.abs(d["t"] - t_lab))); ns, nst = int(d["n_spots"][k]), int(d["n_stable"][k])
            ax.set_title(f"{row['label']}, {ns}, {nst}", fontsize=6.5, pad=2)
            if c == 0:
                ax.set_ylabel(f"n {na}/{ni}", fontsize=7)
            if r == len(HILL) - 1:
                ax.set_xlabel(f"{f:g}x edge", fontsize=7)
            for sp in ax.spines.values():
                sp.set_linewidth(0.4)
    fig.suptitle(f"2D initiation sweep, start {start}, D_i {D}, replicate 0: activator at {t_text}, colour 0..b_a = 3 fixed; "
                 f"panel title = label at t = 50, n_spots, n_stable" + ("" if key == "t50" else " (at t = 18)") + "; hex-correct even-r raster", fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{key}_{start}_D{D}.png"
    fig.savefig(out, dpi=100, facecolor="white"); plt.close(fig)
    return out, missing


if __name__ == "__main__":
    keys = sys.argv[1:] or list(SETS)
    total_missing = 0
    for key in keys:
        for start in STARTS:
            for D in D_I:
                out, m = sheet(start, D, key); total_missing += m
                print(f"{out.name}: {out.stat().st_size / 1e6:.2f} MB" + (f", {m} missing panels" if m else ""), flush=True)
    print(f"missing panels in total: {total_missing}")
