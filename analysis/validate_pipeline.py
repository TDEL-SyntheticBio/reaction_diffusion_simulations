"""
Ground-truth tests of the hex-aware pipeline (handoff open item 4: "the pipeline has never been
run against a known answer"). Each synthetic pattern is rasterised onto the engine's even-r
lattice as discs, segmented, re-centroided and measured exactly as a simulation field would be.

Expected answers:
  perfect hex lattice  -> NN CV ~ 0, g(r) sharp peaks at s, sqrt(3) s, 2 s, no mass below s
  Poisson (CSR)        -> NN CV ~ 0.52, g(r) ~ 1 everywhere (hole ~ 0)
  hard core r_ex       -> g(r) = 0 below 2 r_ex (hole ~ 2 r_ex), NN CV below 0.52
Also quantifies the error of treating the index array as a square lattice (what the deposit's
res_analysis scripts do) and of the mirrored component adjacency.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from analysis.hexgeom import cell_centres, component_table, label_components, rasterize_discs, window  # noqa: E402
from analysis.pointstats import RAYLEIGH_CV, PooledG, nn_cv  # noqa: E402
from analysis.synthetic import csr_points, hardcore_points, hex_lattice_points  # noqa: E402


def measure(mask: np.ndarray, min_size: int = 1):
    tab = component_table(mask, min_size=min_size)
    xy = np.column_stack([tab["cx"], tab["cy"]])
    return tab, xy


def index_coords_cv(mask: np.ndarray) -> float:
    """What you get if you centroid in (row, col) index units, as the deposit's scripts do."""
    labels, n = label_components(mask)
    r, c = np.mgrid[0:mask.shape[0], 0:mask.shape[1]]
    areas = np.bincount(labels[labels >= 0], minlength=n)
    cr = np.bincount(labels[labels >= 0], weights=r[labels >= 0], minlength=n) / areas
    cc = np.bincount(labels[labels >= 0], weights=c[labels >= 0], minlength=n) / areas
    return nn_cv(np.column_stack([cc, cr]))


def main() -> None:
    ny = nx = 100
    rng = np.random.default_rng(2026)
    disc_r = 2.2            # ~15-cell discs, like Fig 1C spots
    rows = []

    # 1. perfect hex lattice, spacing 9 cells
    s = 9.0
    pts = hex_lattice_points(ny, nx, s)
    mask = rasterize_discs(ny, nx, pts, disc_r)
    tab, xy = measure(mask)
    G = PooledG(r_max=30, dr=0.5, n_null=100, seed=1); G.add(xy, window(ny, nx)); sm = G.summary()
    cent_err = np.abs(np.sort(xy[:, 0] + 1j * xy[:, 1]) - np.sort(pts[:, 0] + 1j * pts[:, 1])).max()
    g = G.g(); below = np.nansum(g[G.r < s - 1.0])
    rows.append(("hex lattice s=9", len(pts), tab["n_kept"], nn_cv(xy), index_coords_cv(mask), sm["hole"], sm["peak1_r"],
                 f"centroid err {cent_err:.3f}; g mass below s-1: {below:.3f}; peaks expected at 9, 15.6, 18"))
    peaks_r = G.r[(g > 3) & np.isfinite(g)]
    print("   hex lattice g(r) > 3 at r =", np.round(peaks_r, 2))

    # 2. CSR, 5 replicates pooled (points far enough apart to stay separate discs most of the time)
    Gc = PooledG(r_max=30, dr=0.5, n_null=100, seed=2); cvs = []; cvs_idx = []; n_true = n_seg = 0
    for _ in range(5):
        pts = csr_points(ny, nx, 60, rng)
        mask = rasterize_discs(ny, nx, pts, disc_r)
        tab, xy = measure(mask)
        n_true += len(pts); n_seg += tab["n_kept"]
        cvs.append(nn_cv(xy)); cvs_idx.append(index_coords_cv(mask)); Gc.add(xy, window(ny, nx))
    sm = Gc.summary(); g = Gc.g()
    rows.append(("CSR n=60 x5 (discs)", n_true, n_seg, float(np.mean(cvs)), float(np.mean(cvs_idx)), sm["hole"], sm["peak1_r"],
                 f"g mean over r in [6,30]: {np.nanmean(g[(Gc.r > 6)]):.3f} (expect 1); Rayleigh CV {RAYLEIGH_CV:.3f}; merged discs lower n"))
    # 2b. CSR as pure points (no rasterisation) to isolate the statistic from segmentation
    Gp = PooledG(r_max=30, dr=0.5, n_null=100, seed=3); cvp = []
    for _ in range(20):
        pts = csr_points(ny, nx, 60, rng); cvp.append(nn_cv(pts)); Gp.add(pts, window(ny, nx))
    g = Gp.g()
    rows.append(("CSR n=60 x20 (points)", 1200, 1200, float(np.mean(cvp)), np.nan, Gp.summary()["hole"], np.nan,
                 f"g mean all r: {np.nanmean(g):.3f}, sd {np.nanstd(g):.3f} (expect 1 +- noise)"))

    # 3. hard core, exclusion radius 4 -> hole at 8 cells
    Gh = PooledG(r_max=30, dr=0.5, n_null=100, seed=4); cvh = []; n_true = n_seg = 0
    for _ in range(5):
        pts = hardcore_points(ny, nx, 60, 4.0, rng)
        mask = rasterize_discs(ny, nx, pts, disc_r)
        tab, xy = measure(mask)
        n_true += len(pts); n_seg += tab["n_kept"]; cvh.append(nn_cv(xy)); Gh.add(xy, window(ny, nx))
    sm = Gh.summary(); g = Gh.g()
    rows.append(("hard core r_ex=4 x5", n_true, n_seg, float(np.mean(cvh)), np.nan, sm["hole"], sm["peak1_r"],
                 f"expect hole ~8 (2 r_ex); g below 7.5: {np.nansum(g[Gh.r < 7.5]):.3f}; g mean r>12: {np.nanmean(g[Gh.r > 12]):.3f}"))

    # 4. adjacency check: a thin diagonal line is one component under engine adjacency; count under the mirror
    from analysis.hexgeom import neighbor_arrays
    line = np.zeros((ny, nx), dtype=bool)
    r0, c0 = 10, 10
    for k in range(30):                      # walk along the engine's (+1 row, +1 col on even rows) diagonal
        line[r0, c0] = True
        c0 += 1 if r0 % 2 == 0 else 0
        r0 += 1
    _, n_eng = label_components(line)
    mirror = line[:, ::-1]
    _, n_mir = label_components(mirror)       # same shape drawn in mirrored adjacency -> fragments
    rows.append(("diagonal line adjacency", 1, n_eng, np.nan, np.nan, np.nan, np.nan,
                 f"engine adjacency: {n_eng} component; the same line under mirrored (odd-r) adjacency: {n_mir} components"))

    print(f"\n{'pattern':<26}{'n_true':>7}{'n_seg':>6}{'CV phys':>9}{'CV index':>9}{'hole':>7}{'peak1':>7}  notes")
    for name, nt, ns, cv, cvi, hole, pk, note in rows:
        f = lambda v: "   nan" if not np.isfinite(v) else f"{v:6.3f}"
        print(f"{name:<26}{nt:>7}{ns:>6}{f(cv):>9}{f(cvi):>9}{f(hole):>7}{f(pk):>7}  {note}")


if __name__ == "__main__":
    main()
