"""Hex-correct rendering of (Ny, Nx) fields on the engine's even-r lattice (CLAUDE.md §2; visualize_2D.py draws odd-r)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from analysis.hexgeom import cell_centres  # noqa: E402

_ANG = np.deg2rad([30, 90, 150, 210, 270, 330])
_R = 1 / np.sqrt(3)          # circumradius for unit centre-to-centre spacing (pointy-top)


def hex_field(ax, field, cmap="Greens", vmin=0.0, vmax=None, title=None):
    field = np.asarray(field, dtype=float)
    ny, nx = field.shape
    x, y = cell_centres(ny, nx)
    verts = np.stack([x.ravel()[:, None] + _R * np.cos(_ANG), y.ravel()[:, None] + _R * np.sin(_ANG)], axis=-1)
    pc = PolyCollection(verts, array=field.ravel(), cmap=cmap, edgecolors="none", clim=(vmin, vmax if vmax is not None else max(field.max(), 1e-12)))
    ax.add_collection(pc)
    ax.set_xlim(-1, nx + 0.5); ax.set_ylim(-1, (ny - 1) * np.sqrt(3) / 2 + 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=9)
    return pc


def save_panels(fields, titles, path, cmaps=None, suptitle=None, ncols=None):
    n = len(fields); ncols = ncols or n; nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.0 * nrows), squeeze=False)
    for k, ax in enumerate(axes.ravel()):
        if k < n:
            hex_field(ax, fields[k], cmap=(cmaps[k] if cmaps else "Greens"), title=titles[k])
        else:
            ax.axis("off")
    if suptitle:
        fig.suptitle(suptitle, fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)
