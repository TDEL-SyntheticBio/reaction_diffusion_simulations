"""
Measure the inhibitor decay length around ONE arrested domain in the engine, to settle whether
"lambda_i in cells" is sqrt(D_i/gamma) (handoff) or sqrt(1.5 D_i/gamma) (CLAUDE.md §3).

2D: a disc of activator is planted at the centre, no nucleation, run until the active area stops
changing; the exterior inhibitor is fitted to A*K0(r/lambda) (exact for a disc source in the
continuum) and, independently, to the asymptotic slope of log(i*sqrt(r)).
1D: the 1D engine with activator_spike at the centre; exterior fitted to exp(-x/lambda).
"""
from __future__ import annotations

import io, contextlib, sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import k0

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "2D_simulations"))
from simulation_2D import _vectorized_step  # noqa: E402
from analysis.hexgeom import cell_centres, neighbor_arrays, rasterize_discs, component_table, equivalent_radius  # noqa: E402


def params(ba, bi, g, na, ni, D):
    return dict(act_half_sat=1.0, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=0.0,
                inh_diffusion=D, act_prod_rate=ba, inh_prod_rate=bi, inh_decay_rate=g, act_hill_coeff=na, inh_hill_coeff=ni)


def run_single_domain_2d(p, ny=121, nx=121, seed_radius=2.5, seed_level=5.0, dt=0.01, max_steps=30000, check_every=500):
    nbr = neighbor_arrays(ny, nx)
    x, y = cell_centres(ny, nx)
    cx, cy = x[ny // 2, nx // 2], y[ny // 2, nx // 2]
    a = np.where(rasterize_discs(ny, nx, np.array([[cx, cy]]), seed_radius), seed_level, 0.0)
    i = np.zeros_like(a)
    area_hist = []
    for step in range(max_steps):
        a, i = _vectorized_step(a, i, dt, 1.0, p, "juxtacrine", *nbr)
        if step % check_every == 0:
            area_hist.append(int((a > 0.5 * a.max()).sum()) if a.max() > 0 else 0)
            if len(area_hist) >= 4 and len(set(area_hist[-4:])) == 1 and step > 4000:
                break
    return a, i, step, area_hist, (cx, cy)


def fit_exterior_2d(a, i, centre, r_min_offset=1.5, r_max=45.0):
    ny, nx = a.shape
    x, y = cell_centres(ny, nx)
    r = np.hypot(x - centre[0], y - centre[1])
    mask = a > 0.5 * a.max()
    R = equivalent_radius(mask.sum())                    # area-equivalent domain radius (cell area sqrt(3)/2)
    sel = (r > R + r_min_offset) & (r < r_max)
    rr, ii = r[sel], i[sel]
    try:
        (A, lam), _ = curve_fit(lambda r_, A_, l_: A_ * k0(r_ / l_), rr, ii, p0=(ii.max(), 5.0), maxfev=20000)
    except Exception:
        A, lam = np.nan, np.nan
    wall = min(centre[0], centre[1], (nx - 1 + 0.5) - centre[0], (ny - 1) * np.sqrt(3) / 2 - centre[1])
    lam_ref = lam if np.isfinite(lam) and 0 < lam < 20 else 6.0
    far = (r > R + max(6.0, 2 * lam_ref)) & (r < min(r_max, wall - 2 * lam_ref)) & (i > 0)   # cross-check only: biased ~1-2% high
    if far.sum() >= 10:
        slope, _ = np.polyfit(r[far], np.log(i[far] * np.sqrt(r[far])), 1)
    else:
        slope = np.nan
    i_in = i[r < max(R - 1.0, 0.5)].mean()
    nbr_r, nbr_c, nbr_mask, _ = neighbor_arrays(ny, nx)
    has_off_nbr = ((~mask[nbr_r, nbr_c]) & nbr_mask).any(axis=2)
    shell = mask & has_off_nbr                                   # outermost mask shell (engine adjacency)
    rim = (~mask) & (mask[nbr_r, nbr_c] & nbr_mask).any(axis=2)  # exterior cells touching the domain
    return dict(R=R, area=int(mask.sum()), lam_k0=lam, lam_slope=(-1 / slope) if np.isfinite(slope) and slope != 0 else np.nan, i_interior=i_in, i_shell=i[shell].mean(), i_rim=i[rim].mean(),
                i_boundary_over_interior=i[shell].mean() / i_in, a_interior=a[r < max(R - 1, 0.5)].mean())


def run_single_domain_1d(p, N=201, spike=2.0, dt=0.01, max_steps=30000):
    """
    1D engine (1D_simuations/simulation.py) with activator_spike at the centre, run in a SUBPROCESS so the
    1D copy of finding_steady_states never shadows the 2D copy in this interpreter (CLAUDE.md §6).
    N=201 is deliberate isolation, not the published N=100 (Supp Table 1 Fig 1F), so the exterior fit has room.
    The exponential fit is only reported when exactly one contiguous domain exists and it contains the seed.
    """
    import json, subprocess
    code = f"""
import io, contextlib, json, sys, numpy as np
sys.path.insert(0, {str(ROOT / "1D_simuations")!r})
from simulation import run_coupled_neumann
p = {dict(p)!r}
with contextlib.redirect_stdout(io.StringIO()):
    A, I, step, a_ss, i_ss = run_coupled_neumann({N}, {max_steps}, {dt}, 1.0, p, 1e-7, 2000, init_mode="activator_spike",
                                                 activator_type="juxtacrine", spike_value={spike}, save_every=500)
json.dump(dict(a=A[-1].tolist(), i=I[-1].tolist(), step=int(step)), sys.stdout)
"""
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout
    res = json.loads(out)
    a, i, step = np.array(res["a"]), np.array(res["i"]), res["step"]
    xg = np.arange(N) - N // 2
    mask = a > 0.5 * a.max() if a.max() > 0.05 else np.zeros(N, bool)
    runs = np.flatnonzero(np.diff(np.concatenate([[0], mask.astype(int), [0]])))
    n_domains = len(runs) // 2
    single = n_domains == 1 and mask[N // 2]
    R = mask.sum() / 2.0
    if single:
        sel = (np.abs(xg) > R + 3) & (np.abs(xg) < N // 2 - 5)
        slope, _ = np.polyfit(np.abs(xg[sel]), np.log(i[sel]), 1)
        lam = -1 / slope
    else:
        lam = np.nan
    return dict(step=step, n_domains=n_domains, single=single, half_width=R, lam_slope=lam,
                i_interior=i[np.abs(xg) < max(R - 1, 1)].mean(), i_boundary=i[np.abs(np.abs(xg) - R) <= 0.5].mean())


if __name__ == "__main__":
    print("2D single arrested domain (121x121, disc seed radius 2.5 at level 5, no nucleation)")
    print("R = area-equivalent radius in centre-spacing units (cell area sqrt(3)/2); I_shell = inhibitor on the outermost domain shell, I_rim just outside")
    print(f"{'params (b_a,b_i,g,n_a,n_i,D)':<30}{'steps':>7}{'R':>6}{'area':>6}{'lam K0-fit':>11}{'lam slope':>10}{'sqrt(D/g)':>10}{'sqrt(1.5D/g)':>13}{'Ishell/Iin':>11}{'Irim/Iin':>9}{'I_in':>7}")
    for label, pr in [("1C reg spots", params(5, 5, 0.5, 3, 3, 10)),
                      ("1C irreg spots", params(5, 12, 0.5, 3, 3, 10)),
                      ("2L JAPI", params(5, 25, 0.5, 10, 4, 20)),
                      ("3L JAPI", params(5, 20, 0.5, 10, 4, 20)),
                      ("3H low_di", params(5, 25, 0.5, 10, 4, 5))]:
        a, i, step, hist, c = run_single_domain_2d(pr)
        if a.max() < 0.1 or (a > 0.5 * a.max()).sum() < 3:
            print(f"{label:<30}{step:>7}  domain died (max a {a.max():.3f}); area history {hist}")
            continue
        f = fit_exterior_2d(a, i, c)
        D, g = pr["inh_diffusion"], pr["inh_decay_rate"]
        print(f"{label:<30}{step:>7}{f['R']:>6.2f}{f['area']:>6}{f['lam_k0']:>11.2f}{f['lam_slope']:>10.2f}{np.sqrt(D/g):>10.2f}{np.sqrt(1.5*D/g):>13.2f}{f['i_boundary_over_interior']:>11.3f}{f['i_rim']/f['i_interior']:>9.3f}{f['i_interior']:>7.2f}"
              + ("" if len(set(hist[-4:])) == 1 else f"   NOT ARRESTED, area history tail {hist[-6:]}"))
    print("\n1D single arrested domain (N=201, activator_spike 2.0 at centre)")
    for label, pr in [("1F irreg", params(5, 6, 0.5, 3, 3, 10)), ("1F periodic", params(5, 3, 0.5, 3, 3, 10))]:
        f = run_single_domain_1d(pr)
        lam_exact = 1 / np.arccosh(1 + pr["inh_decay_rate"] / (2 * pr["inh_diffusion"]))     # exact 1D lattice decay length
        if f["single"]:
            print(f"{label:<14} steps {f['step']:>6}  half-width {f['half_width']:5.1f}  lam slope {f['lam_slope']:5.2f}  "
                  f"sqrt(D/g) {np.sqrt(pr['inh_diffusion']/pr['inh_decay_rate']):.2f}  exact lattice {lam_exact:.2f}  I(R)/I_in {f['i_boundary']/f['i_interior']:.3f}")
        else:
            print(f"{label:<14} steps {f['step']:>6}  NOT A SINGLE DOMAIN ({f['n_domains']} domains): no decay-length fit")
