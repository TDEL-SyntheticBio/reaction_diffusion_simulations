"""
Ground-truth tests of the hex-aware pipeline (handoff open item 4: "the pipeline has never been
run against a known answer"). Each synthetic pattern is rasterised onto the engine's even-r
lattice as discs, segmented, re-centroided and measured exactly as a simulation field would be.

Expected answers (asserted with tolerances; the script exits non-zero on any FAIL):
  perfect hex lattice      -> NN CV ~ 0, centroids within 0.1 cell, g(r) peaks at s, sqrt(3) s, 2 s, no mass below s
  Poisson points           -> NN CV equal to the same-n same-window CSR reference (0.55 here, not the
                              n->inf Rayleigh value 0.52), g(r) ~ 1 everywhere
  Poisson single cells     -> same, with the lattice-site null (centroids on lattice sites)
  Poisson discs r=2.2      -> discs closer than ~5 cells merge, so segmentation imposes a hard core: hole ~ 5,
                              NN CV well below the CSR reference (handoff §3: 'CV comes out low with no mechanism')
  hard core r_ex=4         -> g(r) = 0 below 2 r_ex, hole 8 +- 0.5
  mirrored adjacency       -> a 30-cell engine diagonal is 1 component; mirrored it fragments
Peak positions are bin centres (dr = 0.5), so a spacing on a bin edge reads +-dr/2.

The suite ends with an ENGINE fate check through the public API alone (run_coupled_hex with
init_mode='spike_steady_state', no custom seeding): at b_a=5, gamma=0.5, n=3/3, D_i=10 a single cell
at level 5 must give one 19-cell domain at a true fixed point for b_i=12, and for b_i=5 one growing
component that splits into twelve components in two six-fold orbits (six inner spots, six outer arcs:
the lattice's six-fold symmetry) and keeps growing. Frame k is the state after (k-1)*2000+1 Euler
steps, i.e. t = 20*(k-1) + 0.01; the deposit's own labels would call that frame t = 20*k.
This is the arrest-versus-replication split; it takes ~4 minutes. Pass --skip-engine to omit it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from analysis.hexgeom import cell_centres, component_table, label_components, rasterize_discs, window  # noqa: E402
from analysis.pointstats import PooledG, csr_reference_cv, nn_cv  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
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


FAILS: list[str] = []


def check(name: str, ok: bool, detail: str) -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILS.append(name)


def engine_fate_check() -> None:
    """Arrest vs self-replication of a single activated cell, through run_coupled_hex only."""
    import contextlib, io
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "2D_simulations"))
    from simulation_2D import run_coupled_hex
    N, steps, every = 141, 30001, 2000                 # frames after 1, 2001, ... Euler steps: t = 0.01 + 20*(k-1)
    frame = {60: 4, 120: 7, 180: 10, 240: 13, 300: 16}
    print("\nengine fate check (public API, spike_steady_state at level 5, threshold 1.5 = 0.3 r_a):")
    for bi in (12, 5):
        p = dict(act_half_sat=1.0, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=1.0, act_prod_rate=5.0,
                 inh_prod_rate=float(bi), inh_decay_rate=0.5, act_hill_coeff=3, inh_hill_coeff=3, inh_diffusion=10.0)
        with contextlib.redirect_stdout(io.StringIO()):
            A, I, step, a_ss, i_ss = run_coupled_hex(N, N, steps, 0.01, 1.0, p, 0.0, 10**9, init_mode="spike_steady_state",
                                                     activator_type="juxtacrine", spike_value=5.0, save_every=every,
                                                     nucleation_rate=0.0, noise_amplitude=0.0)
        check(f"b_i={bi} engine used the spike_value fallback", a_ss == 5.0 and i_ss == 5.0 and len(A) == 17, f"a_ss {a_ss}, i_ss {i_ss}, {len(A)} frames")
        comps = {}
        for t, k in frame.items():
            labels, n = label_components(A[k] > 1.5)
            areas = sorted(np.bincount(labels[labels >= 0], minlength=n).tolist()) if n else []
            comps[t] = (n, areas)
        per_tile = (np.abs(A[-1] - A[-2]).sum() + np.abs(I[-1] - I[-2]).sum()) / (2 * N * N)
        print("   " + "; ".join(f"t={t}: {n} comp, areas {ar if len(ar) <= 3 else [ar[0], '...', ar[-1]]}" for t, (n, ar) in comps.items()) + f"; per-tile change {per_tile:.1e}")
        if bi == 12:
            check("b_i=12 one 19-cell domain, identical t=60..300", all(n == 1 and ar == [19] for n, ar in comps.values()), f"{comps}")
            check("b_i=12 is at a fixed point", per_tile < 1e-10, f"per-tile change over the last 2000 steps {per_tile:.1e}")
        else:
            grow = [comps[t][1][0] for t in (60, 120, 180)]
            check("b_i=5 one growing component at t=60,120,180", all(comps[t][0] == 1 for t in (60, 120, 180)) and grow[0] < grow[1] < grow[2], f"areas {grow}")
            n240, ar240 = comps[240]
            orbits = {a: ar240.count(a) for a in set(ar240)}
            six_fold = n240 == 12 and len(orbits) == 2 and all(c == 6 for c in orbits.values())
            check("b_i=5 splits into twelve components in two six-fold orbits by t=240", six_fold, f"{n240} components, areas {ar240}")
            check("b_i=5 still growing at t=300", sum(comps[300][1]) > sum(ar240) and per_tile > 1e-3, f"total {sum(ar240)} -> {sum(comps[300][1])}, per-tile change {per_tile:.1e}")


def main() -> None:
    ny = nx = 100
    rng = np.random.default_rng(2026)
    disc_r = 2.2            # ~18-cell discs (16-19), like Fig 1C spots
    rows = []
    dr = 0.5

    # 1. perfect hex lattice, spacing 9 cells
    s = 9.0
    pts = hex_lattice_points(ny, nx, s)
    mask = rasterize_discs(ny, nx, pts, disc_r)
    tab, xy = measure(mask)
    G = PooledG(r_max=30, dr=dr, n_null=100, seed=1); G.add(xy, window(ny, nx)); sm = G.summary()
    cent_err = cKDTree(pts).query(xy)[0].max() if tab["n_kept"] == len(pts) else np.inf
    g = G.g(); below = np.nansum(g[G.r < s - 1.0])
    rows.append(("hex lattice s=9", len(pts), tab["n_kept"], nn_cv(xy), index_coords_cv(mask), sm["hole"], sm["peak1_r"],
                 f"centroid err {cent_err:.3f}; g mass below s-1: {below:.3f}; peaks expected at 9, 15.6, 18"))
    peaks_r = G.r[(g > 3) & np.isfinite(g)]
    print("   hex lattice g(r) > 3 at r =", np.round(peaks_r, 2))
    check("lattice count", tab["n_kept"] == len(pts), f"{tab['n_kept']} of {len(pts)} discs recovered")
    check("lattice centroids", cent_err < 0.1, f"max centroid error {cent_err:.3f} cells")
    check("lattice NN CV", nn_cv(xy) < 0.01, f"CV {nn_cv(xy):.4f}")
    check("lattice g below s", below == 0, f"g mass below s-1 = {below:.3f}")
    check("lattice first peak", abs(sm["peak1_r"] - s) <= dr, f"peak1 {sm['peak1_r']:.2f} vs {s} (+-{dr})")
    check("lattice hole", s - 1.0 <= sm["hole"] <= s, f"hole {sm['hole']:.2f} in [{s-1}, {s}]")
    for target in (np.sqrt(3) * s, 2 * s):
        check(f"lattice shell at {target:.1f}", np.any(np.abs(peaks_r - target) <= dr), f"peaks found at {np.round(peaks_r, 2)}")

    # 2. CSR, 5 replicates pooled (points far enough apart to stay separate discs most of the time)
    Gc = PooledG(r_max=30, dr=0.5, n_null=100, seed=2); cvs = []; cvs_idx = []; n_true = n_seg = 0
    for _ in range(5):
        pts = csr_points(ny, nx, 60, rng)
        mask = rasterize_discs(ny, nx, pts, disc_r)
        tab, xy = measure(mask)
        n_true += len(pts); n_seg += tab["n_kept"]
        cvs.append(nn_cv(xy)); cvs_idx.append(index_coords_cv(mask)); Gc.add(xy, window(ny, nx))
    sm = Gc.summary(); g = Gc.g()
    ref_plain, _ = csr_reference_cv(60, window(ny, nx), rng, n_real=200)
    rows.append(("CSR n=60 x5 (discs)", n_true, n_seg, float(np.mean(cvs)), float(np.mean(cvs_idx)), sm["hole"], sm["peak1_r"],
                 f"g mean over r in [6,30]: {np.nanmean(g[(Gc.r > 6)]):.3f} (expect 1); CSR reference CV {ref_plain:.3f}; merged discs impose a hard core"))
    check("Poisson discs merge -> hard core", 4.0 <= sm["hole"] <= 6.0 and np.mean(cvs) < ref_plain - 0.1, f"hole {sm['hole']:.2f}, CV {np.mean(cvs):.3f} vs CSR {ref_plain:.3f}")
    check("Poisson discs g far", abs(np.nanmean(g[Gc.r > 6]) - 1) < 0.1, f"mean g(r>6) = {np.nanmean(g[Gc.r > 6]):.3f}")
    # 2b. CSR as pure points (no rasterisation) to isolate the statistic from segmentation
    Gp = PooledG(r_max=30, dr=0.5, n_null=100, seed=3); cvp = []
    for _ in range(20):
        pts = csr_points(ny, nx, 60, rng); cvp.append(nn_cv(pts)); Gp.add(pts, window(ny, nx))
    g = Gp.g()
    rows.append(("CSR n=60 x20 (points)", 1200, 1200, float(np.mean(cvp)), np.nan, Gp.summary()["hole"], np.nan,
                 f"g mean all r: {np.nanmean(g):.3f}, sd {np.nanstd(g):.3f} (expect 1 +- noise); CSR reference CV {ref_plain:.3f}"))
    se = 0.071 / np.sqrt(20)      # sd of the n=60 CSR CV across realisations is ~0.071
    check("Poisson points CV = window reference", abs(np.mean(cvp) - ref_plain) < 3 * se, f"CV {np.mean(cvp):.3f} vs reference {ref_plain:.3f} (3 SE = {3*se:.3f})")
    check("Poisson points g = 1", abs(np.nanmean(g) - 1) < 0.03, f"mean g {np.nanmean(g):.3f}")
    # 2c. Poisson SINGLE CELLS on the lattice with the lattice-site null: centroids on lattice sites need this null
    from analysis.hexgeom import cell_centres
    X, Y = cell_centres(ny, nx); lattice_xy = np.column_stack([X.ravel(), Y.ravel()])
    Gl = PooledG(r_max=30, dr=0.5, n_null=100, seed=5, lattice_xy=lattice_xy); Gu = PooledG(r_max=30, dr=0.5, n_null=100, seed=6)
    for _ in range(20):
        cells = lattice_xy[rng.choice(len(lattice_xy), 120, replace=False)]
        Gl.add(cells, window(ny, nx)); Gu.add(cells, window(ny, nx))
    rms_l = np.sqrt(np.nanmean((Gl.g()[Gl.r < 10] - 1) ** 2)); rms_u = np.sqrt(np.nanmean((Gu.g()[Gu.r < 10] - 1) ** 2))
    rows.append(("CSR single cells x20", 2400, 2400, np.nan, np.nan, Gl.summary()["hole"], np.nan,
                 f"rms(g-1) for r<10: lattice null {rms_l:.3f}, continuous null {rms_u:.3f} (lattice shells leak through the continuous null)"))
    check("single cells need lattice null", rms_l < 0.1 and rms_u > 0.2, f"rms lattice {rms_l:.3f}, continuous {rms_u:.3f}")

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
    check("hard-core hole", abs(sm["hole"] - 8.0) <= 0.5, f"hole {sm['hole']:.2f} vs 8.0")
    check("hard-core g below 2 r_ex", np.nansum(g[Gh.r < 7.5]) == 0, f"g mass below 7.5 = {np.nansum(g[Gh.r < 7.5]):.3f}")
    check("hard-core g far", abs(np.nanmean(g[Gh.r > 12]) - 1) < 0.1, f"mean g(r>12) = {np.nanmean(g[Gh.r > 12]):.3f}")

    # 4. adjacency check: a thin diagonal line is one component under engine adjacency; count under the mirror
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
    check("engine diagonal is connected", n_eng == 1 and n_mir > 1, f"engine {n_eng}, mirrored {n_mir}")

    print(f"\n{'pattern':<26}{'n_true':>7}{'n_seg':>6}{'CV phys':>9}{'CV index':>9}{'hole':>7}{'peak1':>7}  notes")
    for name, nt, ns, cv, cvi, hole, pk, note in rows:
        f = lambda v: "   nan" if not np.isfinite(v) else f"{v:6.3f}"
        print(f"{name:<26}{nt:>7}{ns:>6}{f(cv):>9}{f(cvi):>9}{f(hole):>7}{f(pk):>7}  {note}")
    if "--skip-engine" not in sys.argv:
        engine_fate_check()
    print(f"\n{len(FAILS)} FAIL" + (": " + ", ".join(FAILS) if FAILS else ""))
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
