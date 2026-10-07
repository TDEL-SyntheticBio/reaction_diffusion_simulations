"""
Point-pattern statistics on physical coordinates: nearest-neighbour spacing CV and a pair
correlation g(r) normalised against a null generated in the same window (CLAUDE.md §14,
handoff §3). Pooling across replicates is by summing counts, never by averaging curves.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

RAYLEIGH_CV = float(np.sqrt(4 / np.pi - 1))   # NN-spacing CV of a 2D Poisson process, 0.5227


def nn_spacing(xy: np.ndarray) -> np.ndarray:
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 2:
        return np.array([])
    d, _ = cKDTree(xy).query(xy, k=2)
    return d[:, 1]


def nn_cv(xy: np.ndarray) -> float:
    d = nn_spacing(xy)
    return float(d.std() / d.mean()) if d.size > 1 and d.mean() > 0 else np.nan


def pair_counts(xy: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Histogram of all unordered pairwise distances."""
    xy = np.asarray(xy, dtype=float)
    n = len(xy)
    if n < 2:
        return np.zeros(len(edges) - 1)
    iu = np.triu_indices(n, k=1)
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))[iu]
    return np.histogram(d, bins=edges)[0].astype(float)


def csr_null_counts(n: int, win: tuple[float, float, float, float], edges: np.ndarray,
                    n_real: int, rng: np.random.Generator) -> np.ndarray:
    """Mean pair-count histogram of n uniform points in the window, over n_real realisations."""
    xmin, xmax, ymin, ymax = win
    acc = np.zeros(len(edges) - 1)
    for _ in range(n_real):
        pts = np.column_stack([rng.uniform(xmin, xmax, n), rng.uniform(ymin, ymax, n)])
        acc += pair_counts(pts, edges)
    return acc / n_real


class PooledG:
    """Accumulate data and null pair counts across replicates; g = sum(H_data)/sum(mean H_null)."""

    def __init__(self, r_max: float = 30.0, dr: float = 0.5, n_null: int = 200, seed: int = 0):
        self.edges = np.arange(0.0, r_max + dr, dr)
        self.r = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.h_data = np.zeros(len(self.r))
        self.h_null = np.zeros(len(self.r))
        self.n_null = n_null
        self.rng = np.random.default_rng(seed)
        self.n_points = []

    def add(self, xy: np.ndarray, win: tuple[float, float, float, float]) -> None:
        xy = np.asarray(xy, dtype=float)
        self.n_points.append(len(xy))
        self.h_data += pair_counts(xy, self.edges)
        self.h_null += csr_null_counts(len(xy), win, self.edges, self.n_null, self.rng)

    def g(self) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.h_null > 0, self.h_data / self.h_null, np.nan)

    def summary(self, rise_level: float = 0.5) -> dict[str, float]:
        """
        hole: r at which g first rises through rise_level (linear interpolation); nan if never.
        peak1_r / peak1_g: position and height of the first local maximum after the hole.
        trough_g: minimum of g between the first and second peaks (nan if no second peak).
        """
        g = self.g(); r = self.r
        out = {"hole": np.nan, "peak1_r": np.nan, "peak1_g": np.nan, "trough_r": np.nan,
               "trough_g": np.nan, "peak2_r": np.nan, "peak2_g": np.nan}
        valid = np.isfinite(g)
        idx = np.nonzero(valid & (g >= rise_level))[0]
        if idx.size == 0:
            return out
        k = idx[0]
        if k > 0 and np.isfinite(g[k - 1]):
            out["hole"] = float(r[k - 1] + (rise_level - g[k - 1]) / (g[k] - g[k - 1]) * (r[k] - r[k - 1]))
        else:
            out["hole"] = float(r[k])
        gs = np.where(valid, g, -np.inf)
        peaks = [i for i in range(k, len(g) - 1) if gs[i] >= gs[i - 1] and gs[i] > gs[i + 1]]
        if peaks:
            p1 = peaks[0]; out["peak1_r"], out["peak1_g"] = float(r[p1]), float(g[p1])
            later = [p for p in peaks[1:] if g[p] > 1.0 or p == peaks[-1]]
            if len(peaks) > 1:
                p2 = peaks[1]
                seg = slice(p1, p2 + 1)
                t = p1 + int(np.nanargmin(np.where(valid[seg], g[seg], np.inf)))
                out.update(trough_r=float(r[t]), trough_g=float(g[t]), peak2_r=float(r[p2]), peak2_g=float(g[p2]))
        return out
