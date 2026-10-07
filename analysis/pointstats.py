"""
Point-pattern statistics on physical coordinates: nearest-neighbour spacing CV (plain, border-corrected,
and a same-n same-window CSR reference) and a pair correlation g(r) normalised against a null generated
in the same window (CLAUDE.md §14, handoff §3). Pooling across replicates is by summing counts, never by
averaging curves. Feature positions carry Poisson-bootstrap uncertainties and pair counts.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

RAYLEIGH_CV = float(np.sqrt(4 / np.pi - 1))   # NN-spacing CV of a 2D Poisson process as n -> inf, 0.5227.
# In a finite window the plain CV of a Poisson pattern is higher (0.55 at n = 20-90 in a 100x100 window):
# compare against csr_reference_cv(n, window), never against this constant.


def nn_spacing(xy: np.ndarray) -> np.ndarray:
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 2:
        return np.array([])
    d, _ = cKDTree(xy).query(xy, k=2)
    return d[:, 1]


def nn_cv(xy: np.ndarray) -> float:
    d = nn_spacing(xy)
    return float(d.std() / d.mean()) if d.size > 1 and d.mean() > 0 else np.nan


def border_cv(xy: np.ndarray, win: tuple[float, float, float, float]) -> tuple[float, int]:
    """Hanisch border correction: keep points whose NN distance is below their distance to the window edge."""
    xy = np.asarray(xy, dtype=float)
    xmin, xmax, ymin, ymax = win
    if len(xy) < 3:
        return np.nan, 0
    d = nn_spacing(xy)
    edge = np.minimum.reduce([xy[:, 0] - xmin, xmax - xy[:, 0], xy[:, 1] - ymin, ymax - xy[:, 1]])
    keep = d < edge
    return (float(d[keep].std() / d[keep].mean()) if keep.sum() > 2 else np.nan), int(keep.sum())


def csr_reference_cv(n: int, win: tuple[float, float, float, float], rng: np.random.Generator, n_real: int = 300) -> tuple[float, float]:
    """Plain and border-corrected NN CV of n uniform points in this window: the finite-n, finite-window Poisson reference."""
    xmin, xmax, ymin, ymax = win
    plain, border = [], []
    for _ in range(n_real):
        pts = np.column_stack([rng.uniform(xmin, xmax, n), rng.uniform(ymin, ymax, n)])
        plain.append(nn_cv(pts)); border.append(border_cv(pts, win)[0])
    return float(np.nanmean(plain)), float(np.nanmean(border))


def pair_counts(xy: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Histogram of all unordered pairwise distances (bins [a, b), last bin closed)."""
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 2:
        return np.zeros(len(edges) - 1)
    return np.histogram(pdist(xy), bins=edges)[0].astype(float)


def csr_null_counts(n: int, win: tuple[float, float, float, float], edges: np.ndarray, n_real: int,
                    rng: np.random.Generator, lattice_xy: np.ndarray | None = None) -> np.ndarray:
    """
    Mean pair-count histogram of n random points in the window over n_real realisations.
    With lattice_xy (the (N,2) cell centres) the points are random distinct lattice cells instead of
    uniform continuous points: use this when components are small (< ~10 cells), because centroids of
    small components sit on or near lattice sites and a continuous null does not cancel the lattice shells.
    """
    xmin, xmax, ymin, ymax = win
    acc = np.zeros(len(edges) - 1)
    for _ in range(n_real):
        if lattice_xy is None:
            pts = np.column_stack([rng.uniform(xmin, xmax, n), rng.uniform(ymin, ymax, n)])
        else:
            pts = lattice_xy[rng.choice(len(lattice_xy), size=min(n, len(lattice_xy)), replace=False)]
        acc += pair_counts(pts, edges)
    return acc / n_real


class PooledG:
    """Accumulate data and null pair counts across replicates; g = sum(H_data) / sum(mean H_null)."""

    def __init__(self, r_max: float = 30.0, dr: float = 0.5, n_null: int = 200, seed: int = 0, lattice_xy: np.ndarray | None = None):
        self.edges = np.arange(0.0, r_max + dr, dr)
        self.r = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.dr = dr
        self.h_data = np.zeros(len(self.r))
        self.h_null = np.zeros(len(self.r))
        self.n_null = n_null
        self.rng = np.random.default_rng(seed)
        self.lattice_xy = lattice_xy
        self.n_points: list[int] = []

    def add(self, xy: np.ndarray, win: tuple[float, float, float, float]) -> None:
        xy = np.asarray(xy, dtype=float)
        self.n_points.append(len(xy))
        self.h_data += pair_counts(xy, self.edges)
        self.h_null += csr_null_counts(len(xy), win, self.edges, self.n_null, self.rng, self.lattice_xy)

    def g(self, h_data: np.ndarray | None = None) -> np.ndarray:
        h = self.h_data if h_data is None else h_data
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.h_null > 0, h / self.h_null, np.nan)

    # ---- feature extraction -------------------------------------------------------------------------
    def _features(self, h_data: np.ndarray, rise_level: float, min_pairs_bin: int, z_min: float) -> dict[str, float]:
        g = self.g(h_data); r = self.r; hn = self.h_null
        out = {"hole": np.nan, "peak1_r": np.nan, "peak1_g": np.nan, "peak1_pairs": 0.0, "trough_r": np.nan, "trough_g": np.nan,
               "peak2_r": np.nan, "peak2_g": np.nan, "peak2_pairs": 0.0}
        valid = np.isfinite(g)
        idx = np.nonzero(valid & (g >= rise_level))[0]
        if idx.size == 0:
            return out
        k = idx[0]
        if k == 0:
            out["hole"] = 0.0
        elif np.isfinite(g[k - 1]):
            # linear interpolation between bin centres; biased low by up to ~dr for a spike-like rise
            out["hole"] = float(r[k - 1] + (rise_level - g[k - 1]) / (g[k] - g[k - 1]) * (r[k] - r[k - 1]))
        else:
            out["hole"] = float(r[k])

        def accept_peak(lo: int, hi: int, sign: float):
            """Strict interior extremum of sign*g in bins [lo, hi] with enough pairs and a >z_min excess over the null."""
            hi = min(hi, len(g) - 2); lo = max(lo, 1)
            if hi <= lo:
                return None
            gs = np.where(valid, sign * g, -np.inf)
            p = lo + int(np.argmax(gs[lo:hi + 1]))
            if not (gs[p] > gs[p - 1] and gs[p] > gs[p + 1]):
                return None
            if sign > 0:
                if h_data[p] < min_pairs_bin or (h_data[p] - hn[p]) < z_min * np.sqrt(max(hn[p], 1e-12)):
                    return None
            return p

        # first coordination shell: capped at 1.4*hole + 1 so the sqrt(3) shell of a triangular arrangement is excluded
        hi1 = int(np.searchsorted(r, 1.4 * out["hole"] + 1.0))
        p1 = accept_peak(k, hi1, +1.0)
        if p1 is None:
            return out
        out.update(peak1_r=float(r[p1]), peak1_g=float(g[p1]), peak1_pairs=float(h_data[p1]))
        lo2 = int(np.searchsorted(r, 1.5 * r[p1])); hi2 = int(np.searchsorted(r, 2.5 * r[p1]))
        p2 = accept_peak(lo2, hi2, +1.0)
        if p2 is not None:
            out.update(peak2_r=float(r[p2]), peak2_g=float(g[p2]), peak2_pairs=float(h_data[p2]))
            t = accept_peak(p1 + 1, p2 - 1, -1.0)
            if t is not None:
                out.update(trough_r=float(r[t]), trough_g=float(g[t]))
        return out

    def summary(self, rise_level: float = 0.5, min_pairs_bin: int = 20, z_min: float = 3.0, n_boot: int = 400,
                min_total_pairs: int = 100) -> dict[str, float]:
        """
        hole: r where g first rises through rise_level (interpolated between bin centres; 0 if the first bin already does).
        peak1: strict interior maximum of g in [hole, 1.4*hole + 1] with >= min_pairs_bin pairs and a > z_min sigma excess
               over the null; peak2 likewise in [1.5, 2.5] x peak1_r; trough: interior minimum between them.
        *_sd: Poisson-bootstrap standard deviation (h_data resampled n_boot times). Everything is nan when the pooled
        data hold fewer than min_total_pairs pairs.
        """
        out = self._features(self.h_data, rise_level, min_pairs_bin, z_min)
        out.update(total_pairs=float(self.h_data.sum()), n_replicates=float(len(self.n_points)),
                   n_points_mean=float(np.mean(self.n_points)) if self.n_points else 0.0, hole_sd=np.nan, peak1_r_sd=np.nan, peak2_r_sd=np.nan)
        if self.h_data.sum() < min_total_pairs:
            for key in ("hole", "peak1_r", "peak1_g", "trough_r", "trough_g", "peak2_r", "peak2_g"):
                out[key] = np.nan
            return out
        if n_boot > 0:
            rng = np.random.default_rng(12345)
            boots = [self._features(rng.poisson(self.h_data).astype(float), rise_level, min_pairs_bin, z_min) for _ in range(n_boot)]
            for key in ("hole", "peak1_r", "peak2_r"):
                vals = np.array([b[key] for b in boots], dtype=float)
                out[key + "_sd"] = float(np.nanstd(vals)) if np.isfinite(vals).sum() > 10 else np.nan
        return out
