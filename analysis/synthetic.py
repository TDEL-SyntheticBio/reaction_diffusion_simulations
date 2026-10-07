"""Synthetic point patterns with known answers, for validating the pipeline on the even-r lattice."""
from __future__ import annotations

import numpy as np

from .hexgeom import SQRT3_2, window


def hex_lattice_points(ny: int, nx: int, spacing: float, margin: float = 2.0) -> np.ndarray:
    """A perfect triangular lattice of points with the given spacing inside the window."""
    xmin, xmax, ymin, ymax = window(ny, nx)
    pts = []
    row = 0
    y = ymin + margin
    while y <= ymax - margin:
        x0 = xmin + margin + (0.5 * spacing if row % 2 else 0.0)
        x = x0
        while x <= xmax - margin:
            pts.append((x, y)); x += spacing
        y += spacing * SQRT3_2; row += 1
    return np.array(pts)


def csr_points(ny: int, nx: int, n: int, rng: np.random.Generator) -> np.ndarray:
    xmin, xmax, ymin, ymax = window(ny, nx)
    return np.column_stack([rng.uniform(xmin, xmax, n), rng.uniform(ymin, ymax, n)])


def hardcore_points(ny: int, nx: int, n: int, r_ex: float, rng: np.random.Generator,
                    max_tries: int = 200000) -> np.ndarray:
    """Random sequential adsorption with minimum centre separation 2*r_ex (keep density well below jamming)."""
    xmin, xmax, ymin, ymax = window(ny, nx)
    pts: list[tuple[float, float]] = []
    tries = 0
    while len(pts) < n and tries < max_tries:
        tries += 1
        p = (rng.uniform(xmin, xmax), rng.uniform(ymin, ymax))
        if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 >= (2 * r_ex) ** 2 for q in pts):
            pts.append(p)
    if len(pts) < n:
        raise RuntimeError(f"RSA placed only {len(pts)}/{n} points; lower n or r_ex")
    return np.array(pts)
