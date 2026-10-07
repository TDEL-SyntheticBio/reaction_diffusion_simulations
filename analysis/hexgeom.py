"""
Geometry of the engine's even-r hex lattice (CLAUDE.md §2) and hex-aware segmentation.

Everything here is built on simulation_2D._build_hex_neighbor_arrays so that adjacency can
never drift from the engine's. Physical coordinates use unit centre-to-centre spacing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

_ENGINE_DIR = Path(__file__).resolve().parent.parent / "2D_simulations"
if str(_ENGINE_DIR) not in sys.path:
    sys.path.insert(0, str(_ENGINE_DIR))
from simulation_2D import _build_hex_neighbor_arrays  # noqa: E402

SQRT3_2 = np.sqrt(3.0) / 2.0


def cell_centres(ny: int, nx: int) -> tuple[np.ndarray, np.ndarray]:
    """Physical (x, y) of every cell, even rows shifted right by half a cell."""
    r, c = np.mgrid[0:ny, 0:nx]
    x = c + 0.5 * (r % 2 == 0)
    y = r * SQRT3_2
    return x.astype(float), y.astype(float)


def window(ny: int, nx: int) -> tuple[float, float, float, float]:
    """Bounding box (xmin, xmax, ymin, ymax) of the cell centres."""
    return 0.0, nx - 1 + 0.5, 0.0, (ny - 1) * SQRT3_2


def window_area(ny: int, nx: int) -> float:
    xmin, xmax, ymin, ymax = window(ny, nx)
    return (xmax - xmin) * (ymax - ymin)


_NBR_CACHE: dict[tuple[int, int], tuple] = {}


def neighbor_arrays(ny: int, nx: int):
    key = (ny, nx)
    if key not in _NBR_CACHE:
        _NBR_CACHE[key] = _build_hex_neighbor_arrays(ny, nx)
    return _NBR_CACHE[key]


def label_components(mask: np.ndarray) -> tuple[np.ndarray, int]:
    """
    Connected components of a boolean mask under the ENGINE's hex adjacency.

    Returns (labels, n) with labels == -1 outside the mask and 0..n-1 inside.
    """
    mask = np.asarray(mask, dtype=bool)
    ny, nx = mask.shape
    nbr_r, nbr_c, nbr_mask, _ = neighbor_arrays(ny, nx)
    idx = np.arange(ny * nx).reshape(ny, nx)
    src = np.broadcast_to(idx[:, :, None], nbr_r.shape)
    dst = idx[nbr_r, nbr_c]
    keep = nbr_mask & mask[:, :, None] & mask[nbr_r, nbr_c]
    s, d = src[keep], dst[keep]
    n_cells = ny * nx
    adj = coo_matrix((np.ones(s.size, dtype=np.int8), (s, d)), shape=(n_cells, n_cells))
    n_all, lab_all = connected_components(adj, directed=False)
    labels = np.full(ny * nx, -1, dtype=np.int64)
    flat_mask = mask.ravel()
    # relabel only the components that contain mask cells, consecutively
    comp_ids = np.unique(lab_all[flat_mask])
    remap = {cid: k for k, cid in enumerate(comp_ids)}
    labels[flat_mask] = [remap[v] for v in lab_all[flat_mask]]
    return labels.reshape(ny, nx), len(comp_ids)


def component_table(mask: np.ndarray, min_size: int = 1) -> dict[str, np.ndarray]:
    """
    Per-component area (cells) and physical centroid, dropping components smaller than min_size.
    """
    labels, n = label_components(mask)
    ny, nx = mask.shape
    x, y = cell_centres(ny, nx)
    areas = np.bincount(labels[labels >= 0], minlength=n).astype(float)
    cx = np.bincount(labels[labels >= 0], weights=x[labels >= 0], minlength=n) / np.maximum(areas, 1)
    cy = np.bincount(labels[labels >= 0], weights=y[labels >= 0], minlength=n) / np.maximum(areas, 1)
    keep = areas >= min_size
    return {"area": areas[keep], "cx": cx[keep], "cy": cy[keep], "n_total": n, "n_kept": int(keep.sum())}


def otsu_threshold(values: np.ndarray, nbins: int = 256) -> float:
    v = np.asarray(values, dtype=float).ravel()
    hist, edges = np.histogram(v, bins=nbins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    w = hist.astype(float)
    p = w / w.sum()
    omega = np.cumsum(p)
    mu = np.cumsum(p * centres)
    mu_t = mu[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b = (mu_t * omega - mu) ** 2 / (omega * (1 - omega))
    sigma_b[~np.isfinite(sigma_b)] = -1
    return float(centres[int(np.argmax(sigma_b))])


def segment(field: np.ndarray, method: str = "half_max", rel_var_gate: float = 1e-3, level: float | None = None) -> tuple[np.ndarray, float]:
    """
    Variance-gated segmentation (CLAUDE.md §14). Uniform fields return an all-False or all-True
    mask by LEVEL (relative to 1.0 in dimensionless units), never by Otsu.
    Returns (mask, threshold_used).
    """
    f = np.asarray(field, dtype=float)
    fmax = float(f.max())
    if fmax <= 0:
        return np.zeros_like(f, dtype=bool), np.nan
    if f.std() / max(fmax, 1e-12) < rel_var_gate:          # uniform field: classify by level
        on = f.mean() > 0.5
        return np.full_like(f, on, dtype=bool), np.nan
    if method == "half_max":
        thr = 0.5 * fmax                     # fragile under continuous nucleation: fresh kicks set fmax (see measure_panels.py)
    elif method == "otsu":
        thr = otsu_threshold(f)
    elif method == "fixed":
        if level is None:
            raise ValueError("fixed segmentation needs level")
        thr = float(level)
    else:
        raise ValueError(method)
    return f > thr, thr


def rasterize_discs(ny: int, nx: int, centres_xy: np.ndarray, radius: float | np.ndarray) -> np.ndarray:
    """Boolean field: cell active if its centre lies within `radius` of any disc centre."""
    x, y = cell_centres(ny, nx)
    mask = np.zeros((ny, nx), dtype=bool)
    radius = np.broadcast_to(np.asarray(radius, dtype=float), (len(centres_xy),))
    for (px, py), rad in zip(centres_xy, radius):
        mask |= (x - px) ** 2 + (y - py) ** 2 <= rad ** 2
    return mask
