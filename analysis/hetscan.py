"""
Heterogeneity scan (2D hex): does fixed cell-to-cell expression heterogeneity, as modeled by Wang, Du, Linton and Elowitz
(aem7145, Methods 11.4), widen the range of b_i where patterns form, and on which side of the edge?

Heterogeneity: each cell gets two independent, fixed, mean-corrected lognormal multipliers,
    eta = exp(sigma z - sigma^2 / 2),  z ~ N(0, 1) iid,  E[eta] = 1,
one on b_a and one on b_i (same promoter, two separate viral integrations). sigma 0.3 is Elowitz's reference
value (linear-scale coefficient of variation about 0.31). For a given replicate the z fields are the same at every
b_i position and every sigma (common random numbers), so differences along those axes are not sampling noise.

Starts:
  elowitz  Elowitz's protocol: every cell at the highest stable uniform state of the homogeneous system, then an
           independent Gaussian perturbation with SD 2% of that value on each species. Above the edge that state is
           OFF (zero), and with no basal production the field stays exactly zero: this is expected, not a bug.
  nuc001   our in vitro-like start, as in bandscan.py: field at zero, each cell fires with probability 0.001 per unit
           time, kick = 1.5 x the homogeneous single-cell ignition threshold (bandscan_thresholds.csv, D_i 10).

Grid: JAPI and PAPI (shared promoter, w = 1), D_i 10, b_a 3, gamma 0.5, 100 x 100, t = 100, the band scan's 50 b_i
positions relative to the homogeneous edge, 3 replicates.
  main:    all 9 Hill pairs, sigma 0 and 0.3
  series:  Hill pairs 2/2, 4/4, 10/4, sigma 0.15 and 0.6
Recorded exactly as in bandscan.py (per unit time: area, domains >= 7 cells, persistent domains, field change rate,
mean activator; frames at day 4, t 50, t 100), plus the heterogeneity fields' realised spread.
"""
from __future__ import annotations

import sys, os, json, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import archscreen as X                                   # noqa: E402
import bandscan as BS                                    # noqa: E402
from hexgeom import label_components                     # noqa: E402

ROOT = X.ROOT
OUT = Path(os.environ.get("HET_OUT", ROOT / "results" / "hetscan"))
N = int(os.environ.get("HET_N", 100))
T_END = float(os.environ.get("HET_T", 100.0))
D_I = 10.0
CIRCUITS = BS.CIRCUITS
HILL = BS.HILL
HILL_SERIES = [(2, 2), (4, 4), (10, 4)]
SIG_MAIN = [0.0, 0.3]
SIG_SERIES = [0.15, 0.6]
STARTS = ["elowitz", "nuc001"]
REPS = [0, 1, 2]
FRACS = BS.FRACS


def job_name(job):
    A, B, na, ni, f, sig, start, rep = job
    return f"{'JAPI' if A == 'membrane' else 'PAPI'}_na{na:02d}_ni{ni:02d}_s{sig:.2f}_f{f:.4f}_{start}_r{rep}"


def all_jobs():
    jobs = []
    for rep in REPS:
        for (A, B, _) in CIRCUITS:
            for (na, ni) in HILL:
                sigs = SIG_MAIN + (SIG_SERIES if (na, ni) in HILL_SERIES else [])
                for sig in sigs:
                    for f in FRACS:
                        for st in STARTS:
                            jobs.append((A, B, na, ni, f, sig, st, rep))
    return jobs


def multipliers(sig, rep, n):
    """Two independent mean-corrected lognormal fields (on b_a, on b_i); same z for every sigma and position."""
    za = np.random.default_rng(np.random.SeedSequence([20261120, rep, 0])).standard_normal((n, n))
    zi = np.random.default_rng(np.random.SeedSequence([20261120, rep, 1])).standard_normal((n, n))
    return np.exp(sig * za - sig ** 2 / 2), np.exp(sig * zi - sig ** 2 / 2)


def initial(m, start, rng, n, na, ni, bi):
    a = np.zeros((n, n)); b = np.zeros((n, n))
    if start == "elowitz":
        st = [s for s in BS.states_fast(na, ni, bi) if s[2]]
        a0, b0, _ = max(st, key=lambda s: s[0])
        a = np.maximum(a0 * (1 + 0.02 * rng.standard_normal((n, n))), 0.0)
        b = np.maximum(b0 * (1 + 0.02 * rng.standard_normal((n, n))), 0.0)
    return a, b


def run(job, kick, out=OUT, n=N, t_end=T_END):
    A, B, na, ni, f, sig, start, rep = job
    E = BS.edge(na, ni); bi = f * E
    m = BS.model(A, B, na, ni, bi, D_I, n=n)
    eta_a, eta_i = multipliers(sig, rep, n)
    m.ba = X.BA * eta_a; m.bi = bi * eta_i                 # per-cell production rates (arrays broadcast in step)
    rng = np.random.default_rng(np.random.SeedSequence([20261120, 10 + STARTS.index(start), rep]))
    a, b = initial(m, start, rng, n, na, ni, bi)
    p_fire = 0.001 * X.DT
    steps = int(round(t_end / X.DT)); per = int(round(1 / X.DT))
    save = {int(round(4 * X.DAY / X.DT)): "day4", int(round(50 / X.DT)): "t50", steps: "t100"}
    thr = 0.3 * X.BA
    hist = []
    ts = {k: [] for k in ("t", "area", "n7", "area7", "np7", "areap7", "field_rate", "mean_a")}
    frames = {}; prev_frames = {}; a_prev = a.copy(); t0 = time.time()
    for s in range(1, steps + 1):
        a, b = m.step(a, b)
        if start == "nuc001":
            fire = (rng.random((n, n)) < p_fire) & (a < X.BA)
            a[fire] = kick
        if s % per == 0:
            on = a > thr; lab, k = label_components(on)
            big = np.zeros_like(on)
            if k:
                sz = np.bincount(lab[lab >= 0]); ids = np.nonzero(sz >= 7)[0]
                big = np.isin(lab, ids) if len(ids) else big
            else:
                ids = []
            ref = hist[0] if len(hist) == 5 else None
            pers = [j for j in ids if ref is not None and ((lab == j) & ref).any()]
            ts["t"].append(s // per); ts["area"].append(float(on.mean())); ts["n7"].append(int(len(ids)))
            ts["area7"].append(float(big.mean())); ts["np7"].append(int(len(pers)))
            ts["areap7"].append(float(np.isin(lab, pers).mean()) if pers else 0.0)
            ts["field_rate"].append(float(np.abs(a - a_prev).mean() / X.BA)); ts["mean_a"].append(float(a.mean()))
            a_prev = a.copy(); hist.append(big); hist = hist[-5:]
        for key, tag in save.items():
            if s == key - 500:
                prev_frames[tag] = a.copy()
            if s == key:
                frames[tag] = a.copy()
    meas = {k: X.measure(m, frames[k], prev_frames.get(k)) for k in frames}
    name = job_name(job)
    row = dict(name=name, circuit="JAPI" if A == "membrane" else "PAPI", n_a=na, n_i=ni, D_i=D_I, f=f, b_i=bi, edge=E,
               sigma=sig, start=start, rep=rep, kick=kick if start == "nuc001" else np.nan,
               turing_homogeneous=BS.turing_unstable(A, B, na, ni, bi, D_I),
               eta_a_cv=float(eta_a.std() / eta_a.mean()), eta_i_cv=float(eta_i.std() / eta_i.mean()),
               **{f"{k}_{t}": v for t, d in meas.items() for k, v in d.items()},
               b_min=float(b.min()), a_max=float(a.max()), finite=bool(np.isfinite(a).all() and np.isfinite(b).all()),
               seconds=time.time() - t0)
    out = Path(out); (out / "rows").mkdir(parents=True, exist_ok=True); (out / "frames").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "frames" / f"{name}.npz", **{k: v.astype(np.float16) for k, v in frames.items()},
                        **{k: np.array(v) for k, v in ts.items()})
    (out / "rows" / f"{name}.json").write_text(json.dumps(row))
    return row


def kicks():
    import pandas as pd
    t = pd.read_csv(ROOT / "analysis" / "summaries" / "bandscan_thresholds.csv")
    t = t[t.D_i == D_I]
    return {(r.A, r.B, int(r.n_a), int(r.n_i)): float(r.kick) for r in t.itertuples()}


def _run_job(args):
    job, kick = args
    try:
        return run(job, kick)
    except Exception as e:
        (OUT / "failed").mkdir(parents=True, exist_ok=True)
        (OUT / "failed" / f"{job_name(job)}.txt").write_text(repr(e))
        return None


def cmd_run(workers):
    from multiprocessing import Pool
    K = kicks(); jobs = all_jobs()
    todo = [(j, K[j[:4]]) for j in jobs if not (OUT / "rows" / f"{job_name(j)}.json").exists()]
    print(f"{len(jobs)} jobs, {len(todo)} to run, N = {N}, t = {T_END}, {workers} workers", flush=True)
    t0 = time.time()
    with Pool(workers) as p:
        for k, _ in enumerate(p.imap_unordered(_run_job, todo, chunksize=1), 1):
            if k % 100 == 0 or k == len(todo):
                el = time.time() - t0
                print(f"  {k}/{len(todo)}  {el/3600:.2f} h elapsed, about {el/k*(len(todo)-k)/3600:.2f} h left", flush=True)


def cmd_aggregate():
    import pandas as pd
    S = ROOT / "analysis" / "summaries"
    rows = [json.loads(f.read_text()) for f in sorted((OUT / "rows").glob("*.json"))]
    pd.DataFrame(rows).to_csv(S / "hetscan_runs.csv", index=False); print(f"{len(rows)} rows -> hetscan_runs.csv")
    recs = []
    for r in rows:
        z = np.load(OUT / "frames" / f"{r['name']}.npz")
        recs.append(pd.DataFrame({k: z[k] for k in ("t", "area", "n7", "area7", "np7", "areap7", "field_rate", "mean_a")}).assign(name=r["name"]))
    pd.concat(recs).to_csv(S / "hetscan_timeseries.csv.gz", index=False, compression="gzip"); print("time series -> hetscan_timeseries.csv.gz")


def growth_label(ts, t):
    """Domain-area growth over the previous 5 time units, %/day (as on the band scan sheets)."""
    i = list(ts["t"]).index(t); j = list(ts["t"]).index(t - 5)
    area, ap, bp = ts["area"][i], ts["areap7"][i], ts["areap7"][j]
    if area > 0.9: return "full field", "#787878"
    if ap == 0: return "no domains", "#787878"
    if bp == 0: return "forming", "#c81e1e"
    g = (ap - bp) / bp / 5 * X.DAY * 100
    if g > 50: return "forming", "#c81e1e"
    return (f"{g:+.0f}%/day" if abs(g) >= 1 else f"{g:+.1f}%/day"), ("#c81e1e" if g >= 5 else "#006e3c")


def cmd_sheets():
    """Replicate 0: one sheet per (circuit, start, sigma, time); rows = Hill pairs, columns = 13 positions;
    blue dashed line = homogeneous edge; orange frame = homogeneous Turing band; growth label under each panel."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    dest = ROOT / "analysis" / "summaries" / "hetscan_sheets"; dest.mkdir(parents=True, exist_ok=True)
    cols = [FRACS[i] for i in np.unique(np.searchsorted(FRACS, [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.4, 2.0, 2.8, 4.0]).clip(0, len(FRACS) - 1))]
    tval = {"day4": int(round(4 * X.DAY)), "t50": 50, "t100": 100}
    for (A, B, cname) in CIRCUITS:
        for st in STARTS:
            for sig in sorted(set(SIG_MAIN + SIG_SERIES)):
                hill = HILL if sig in SIG_MAIN else HILL_SERIES
                for tkey in ("day4", "t50", "t100"):
                    fig, axs = plt.subplots(len(hill), len(cols), figsize=(1.6 * len(cols), 1.95 * len(hill) + 1.0), squeeze=False)
                    for r, (na, ni) in enumerate(hill):
                        for c, f in enumerate(cols):
                            ax = axs[r, c]; ax.set_xticks([]); ax.set_yticks([])
                            if r == 0: ax.set_title(f"{f:.3g}x", fontsize=7)
                            if c == 0: ax.set_ylabel(f"n {na}/{ni}", fontsize=7)
                            nm = job_name((A, B, na, ni, f, sig, st, 0)); fp = OUT / "frames" / f"{nm}.npz"
                            if not fp.exists():
                                ax.set_facecolor("#eeeeee"); continue
                            z = np.load(fp); row = json.loads((OUT / "rows" / f"{nm}.json").read_text())
                            ax.imshow(X.hex_raster(z[tkey].astype(float)), cmap="Greens", vmin=0, vmax=X.BA, aspect=np.sqrt(3), interpolation="nearest")
                            if row["turing_homogeneous"]:
                                for sp in ax.spines.values(): sp.set_color("#d55e00"); sp.set_linewidth(1.5)
                            txt, col = growth_label(z, tval[tkey])
                            ax.set_xlabel(txt, fontsize=7, color=col, fontweight="bold", labelpad=2)
                    fig.suptitle(f"Heterogeneity scan, {cname}, start {st}, sigma {sig:g}, {tkey}, D_i 10, replicate 0. "
                                 f"Blue dashed = edge; orange frame = Turing band (homogeneous); label = domain-area growth, red >= 5%/day", fontsize=9)
                    fig.tight_layout(rect=[0, 0, 1, 0.96])
                    k = cols.index(next(x for x in cols if x > 1.0))
                    x = 0.5 * (axs[0, k - 1].get_position().x1 + axs[0, k].get_position().x0)
                    fig.add_artist(Line2D([x, x], [0.01, 0.95], transform=fig.transFigure, color="#005ac8", lw=2.5, ls="--"))
                    fig.savefig(dest / f"{cname}_{st}_s{sig:.2f}_{tkey}.png", dpi=80, facecolor="white"); plt.close(fig)
    print("sheets ->", dest)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "help"
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else (os.cpu_count() or 2)
    if cmd == "count":
        print("jobs", len(all_jobs()))
    elif cmd == "run":
        cmd_run(workers)
    elif cmd == "aggregate":
        cmd_aggregate()
    elif cmd == "sheets":
        cmd_sheets()
    else:
        print("usage: hetscan.py count | run [w] | aggregate | sheets")
