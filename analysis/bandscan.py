"""
Band scan (2D hex): how wide is the range of b_i where several persistent, slowing domains form, compared with the
exact Turing band at the same setting?

Circuits: JAPI (A membrane, B diffusible) and PAPI (A diffusible, B diffusible), shared promoter, full self-activation
(w = 1), i.e. the deposited engine's model (archscreen.Model reproduces _vectorized_step exactly; see archscreen.check()).
Axes: all nine Hill pairs {2, 4, 10}^2; b_i at 50 positions relative to the edge (40 log-spaced from 0.3x to 4x plus
10 linear between 0.9x and 1x); D_i 5, 10, 20; starts noise around the uniform state / scattered seeds / nucleation
0.001 (archscreen definitions); 2 replicates. b_a 3, gamma 0.5, D_a 1 for PAPI, 100 x 100, t = 100 (about day 22).
Kicks (seeds and nucleation): 1.5x the single-cell ignition threshold of each (circuit, Hill pair, D_i), capped at 100.

Recorded every time unit: activated area (cells above 0.3 b_a), number and area of domains >= 7 cells, number and area of
persistent domains (>= 7 cells overlapping a >= 7-cell domain 5 time units earlier), field change rate (mean |a(t) -
a(t-1)| / b_a per time unit), mean activator. Frames at day 4 (t = 17.8), t = 50 and t = 100. No label decides anything
here: the slowing-down score is computed afterwards from the time series, with thresholds chosen at analysis time.

Turing band: for every run, whether the highest stable uniform state is unstable to some finite wavelength (general
linear check on the hex lattice, contact signal symbol K in [-1/2, 1], diffusion symbol 6(1 - K)); and per setting the
band boundaries on a fine axis (turing_bands.csv).
"""
from __future__ import annotations

import sys, os, json, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import archscreen as X                                   # noqa: E402
from hexgeom import label_components                     # noqa: E402

ROOT = X.ROOT
OUT = Path(os.environ.get("BAND_OUT", ROOT / "results" / "bandscan"))
N = int(os.environ.get("BAND_N", 100))
T_END = float(os.environ.get("BAND_T", 100.0))
CIRCUITS = [("membrane", "diffusible", "JAPI"), ("diffusible", "diffusible", "PAPI")]
HILL = [(2, 2), (2, 4), (2, 10), (4, 2), (4, 4), (4, 10), (10, 2), (10, 4), (10, 10)]
DIS = [5.0, 10.0, 20.0]
STARTS = ["noise", "seeds", "nuc001"]
REPS = [0, 1]
FRACS = sorted(set(np.round(np.concatenate([np.geomspace(0.3, 4.0, 40), np.linspace(0.9, 1.0, 11)[:-1]]), 4)))
W, PROM = 1.0, "shared"


from functools import lru_cache


@lru_cache(maxsize=None)
def edge(na, ni):
    return X.ref_edge(na, ni)


def states_fast(na, ni, bi, ba=None, ga=None):
    """Uniform steady states for shared promoter, w = 1 (vectorised; checked against archscreen.uniform_states)."""
    ba = X.BA if ba is None else ba; ga = X.GA if ga is None else ga
    k = bi / (ga * ba)                                   # b = k a

    def g(a):
        a = np.asarray(a, float); pa = a ** na; inh = (k * a) ** ni
        return ba * pa / (pa + inh + 1) - a
    A = np.concatenate([[0.0], np.logspace(-6, np.log10(ba), 6000)]); G = g(A)
    idx = np.nonzero((np.sign(G[:-1]) != np.sign(G[1:])) & (G[:-1] != 0) & (G[1:] != 0))[0]
    lo, hi = A[idx].copy(), A[idx + 1].copy(); glo = G[idx]
    for _ in range(60):
        mid = 0.5 * (lo + hi); gm = g(mid); same = np.sign(gm) == np.sign(glo)
        lo = np.where(same, mid, lo); hi = np.where(same, hi, mid)
    roots = list(A[np.nonzero(G == 0)[0]]) + list(0.5 * (lo + hi))
    out = []
    for a in roots:
        b = k * a; e = 1e-7 * max(1.0, a, b)
        def f(x, y):
            px = x ** na; iy = y ** ni; h = px / (px + iy + 1)
            return np.array([ba * h - x, bi * h - ga * y])
        f0 = f(a, b); J = np.column_stack([(f(a + e, b) - f0) / e, (f(a, b + e) - f0) / e])
        out.append((float(a), float(b), bool(np.trace(J) < 0 and np.linalg.det(J) > 0)))
    return out


def model(A, B, na, ni, bi, D, n=N):
    return X.Model(A, B, W, PROM, na, ni, bi, di=D, n=n)


# ------------------------------------------------------------------ linear check (Turing) on the hex lattice
_K = np.linspace(-0.5, 1.0, 3001)


def turing_unstable(A, B, na, ni, bi, D):
    st = [s for s in states_fast(na, ni, bi) if s[2]]
    if not st:
        return False
    a, b, _ = max(st, key=lambda s: s[0])

    def prod(sA, sB):
        pa = sA ** na if sA > 0 else 0.0
        inh = sB ** ni if sB > 0 else 0.0
        h = pa / (pa + inh + 1)
        return np.array([X.BA * h, bi * h])
    e = 1e-7 * max(1.0, a, b); f0 = prod(a, b)
    dA = (prod(a + e, b) - f0) / e; dB = (prod(a, b + e) - f0) / e
    kA = _K if A == "membrane" else np.ones_like(_K); kB = _K if B == "membrane" else np.ones_like(_K)
    L = 6 * (1 - _K)
    DA = X.DA if A == "diffusible" else 0.0; DB = D if B == "diffusible" else 0.0
    j11 = dA[0] * kA - 1 - DA * L; j12 = dB[0] * kB; j21 = dA[1] * kA; j22 = dB[1] * kB - X.GA - DB * L
    tr = j11 + j22; det = j11 * j22 - j12 * j21
    re = (0.5 * (tr + np.sqrt((tr ** 2 - 4 * det).astype(complex)))).real
    fin = _K < 1 - 1e-9
    # exclude the marginal fold itself: require the uniform mode to be clearly stable
    return bool(re[fin].max() > 1e-6 and re[~fin].max() < -1e-6)


def turing_bands():
    rows = []
    fine = np.linspace(0.3, 1.0, 1401)
    for (A, B, name) in CIRCUITS:
        for (na, ni) in HILL:
            E = edge(na, ni)
            for D in DIS:
                t = np.array([turing_unstable(A, B, na, ni, f * E, D) for f in fine])
                rows.append(dict(circuit=name, n_a=na, n_i=ni, D_i=D, edge=E,
                                 turing_lo=float(fine[t].min()) if t.any() else np.nan,
                                 turing_hi=float(fine[t].max()) if t.any() else np.nan,
                                 turing_width_edge_units=float(t.mean() * 0.7)))
    return rows


# ------------------------------------------------------------------ ignition thresholds per (circuit, Hill pair, D_i)
def ignites(A, B, na, ni, bi, D, h, n=60, t=20.0):
    m = model(A, B, na, ni, bi, D, n=n)
    a = np.zeros((n, n)); b = np.zeros((n, n)); a0 = a.copy(); b0 = b.copy()
    a[n // 2, n // 2] = h
    for _ in range(int(round(t / X.DT))):
        a, b = m.step(a, b); a0, b0 = m.step(a0, b0)
    return int((np.abs(a - a0) > 0.3 * X.BA).sum()) >= 7


def _thr(cfg):
    A, B, na, ni, D = cfg
    bi = 0.7 * edge(na, ni)
    if ignites(A, B, na, ni, bi, D, X.KICK_LO):
        th = X.KICK_LO
    elif not ignites(A, B, na, ni, bi, D, X.KICK_CAP):
        th = None
    else:
        lo, hi = X.KICK_LO, X.KICK_CAP
        for _ in range(7):
            mid = np.sqrt(lo * hi)
            if ignites(A, B, na, ni, bi, D, mid): hi = mid
            else: lo = mid
        th = float(hi)
    return dict(A=A, B=B, n_a=na, n_i=ni, D_i=D, threshold=th, kick=X.kick_height(th))


# ------------------------------------------------------------------ one run
def job_name(job):
    A, B, na, ni, D, f, start, rep = job
    return f"{'JAPI' if A == 'membrane' else 'PAPI'}_na{na:02d}_ni{ni:02d}_D{int(D):02d}_f{f:.4f}_{start}_r{rep}"


def run(job, kick, out=OUT, n=N, t_end=T_END):
    A, B, na, ni, D, f, start, rep = job
    E = edge(na, ni); bi = f * E
    m = model(A, B, na, ni, bi, D, n=n)
    rng = np.random.default_rng(np.random.SeedSequence([20261030, STARTS.index(start), rep]))
    a, b = X.initial(m, start, rng, n, kick)
    p_fire = 0.001 * X.DT
    steps = int(round(t_end / X.DT)); per = int(round(1 / X.DT))
    save = {int(round(4 * X.DAY / X.DT)): "day4", int(round(50 / X.DT)): "t50", steps: "t100"}
    thr = 0.3 * X.BA
    hist = []                      # masks of >= 7-cell domains at the last 5 whole time units
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
        for key, step_at in save.items():
            if s == key - 500:
                prev_frames[step_at] = a.copy()
            if s == key:
                frames[step_at] = a.copy()
    meas = {k: X.measure(m, frames[k], prev_frames.get(k)) for k in frames}
    name = job_name(job)
    row = dict(name=name, circuit="JAPI" if A == "membrane" else "PAPI", n_a=na, n_i=ni, D_i=D, f=f, b_i=bi, edge=E,
               start=start, rep=rep, kick=kick, turing=turing_unstable(A, B, na, ni, bi, D),
               **{f"{k}_{t}": v for t, d in meas.items() for k, v in d.items()},
               b_min=float(b.min()), a_max=float(a.max()), finite=bool(np.isfinite(a).all() and np.isfinite(b).all()),
               seconds=time.time() - t0)
    out = Path(out); (out / "rows").mkdir(parents=True, exist_ok=True); (out / "frames").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "frames" / f"{name}.npz", **{k: v.astype(np.float16) for k, v in frames.items()},
                        **{k: np.array(v) for k, v in ts.items()})
    (out / "rows" / f"{name}.json").write_text(json.dumps(row))
    return row


def _run_job(args):
    job, kick = args
    try:
        return run(job, kick)
    except Exception as e:
        (OUT / "failed").mkdir(parents=True, exist_ok=True)
        (OUT / "failed" / f"{job_name(job)}.txt").write_text(repr(e))
        return None


def all_jobs():
    return [(A, B, na, ni, D, f, st, rep) for rep in REPS for (A, B, _) in CIRCUITS for (na, ni) in HILL for D in DIS
            for f in FRACS for st in STARTS]


# ------------------------------------------------------------------ commands
def cmd_thresholds(workers):
    from multiprocessing import Pool
    f = OUT / "thresholds.json"
    if f.exists():
        print("thresholds already computed:", f); return
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    cfgs = [(A, B, na, ni, D) for (A, B, _) in CIRCUITS for (na, ni) in HILL for D in DIS]
    with Pool(workers) as p:
        rows = p.map(_thr, cfgs, chunksize=1)
    f.write_text(json.dumps(rows, indent=1)); print(f"{len(rows)} thresholds in {time.time() - t0:.0f} s", flush=True)
    (OUT / "turing_bands.json").write_text(json.dumps(turing_bands(), indent=1)); print("Turing bands written", flush=True)


def cmd_run(workers):
    from multiprocessing import Pool
    K = {(r["A"], r["B"], r["n_a"], r["n_i"], r["D_i"]): r["kick"] for r in json.loads((OUT / "thresholds.json").read_text())}
    jobs = all_jobs()
    todo = [(j, K[j[:5]]) for j in jobs if not (OUT / "rows" / f"{job_name(j)}.json").exists()]
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
    pd.DataFrame(rows).to_csv(S / "bandscan_runs.csv", index=False); print(f"{len(rows)} rows -> bandscan_runs.csv")
    recs = []
    for r in rows:
        z = np.load(OUT / "frames" / f"{r['name']}.npz")
        recs.append(pd.DataFrame({k: z[k] for k in ("t", "area", "n7", "area7", "np7", "areap7", "field_rate", "mean_a")}).assign(name=r["name"]))
    pd.concat(recs).to_csv(S / "bandscan_timeseries.csv.gz", index=False, compression="gzip"); print("time series -> bandscan_timeseries.csv.gz")
    pd.DataFrame(json.loads((OUT / "thresholds.json").read_text())).to_csv(S / "bandscan_thresholds.csv", index=False)
    pd.DataFrame(json.loads((OUT / "turing_bands.json").read_text())).to_csv(S / "bandscan_turing_bands.csv", index=False)


def cmd_sheets():
    """Overview sheets, replicate 0, D_i 10: one per (circuit, start, time); rows = Hill pairs, columns = 13 positions."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    dest = ROOT / "analysis" / "summaries" / "bandscan_sheets"; dest.mkdir(parents=True, exist_ok=True)
    cols = [FRACS[i] for i in np.unique(np.searchsorted(FRACS, [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.4, 2.0, 2.8, 4.0]).clip(0, len(FRACS) - 1))]
    for (A, B, cname) in CIRCUITS:
        for st in STARTS:
            for tkey in ("day4", "t50", "t100"):
                fig, axs = plt.subplots(len(HILL), len(cols), figsize=(1.6 * len(cols), 1.75 * len(HILL) + 0.8), squeeze=False)
                for r, (na, ni) in enumerate(HILL):
                    for c, f in enumerate(cols):
                        ax = axs[r, c]; ax.set_xticks([]); ax.set_yticks([])
                        nm = job_name((A, B, na, ni, 10.0, f, st, 0)); fp = OUT / "frames" / f"{nm}.npz"
                        if not fp.exists():
                            ax.set_facecolor("#eeeeee"); continue
                        z = np.load(fp); row = json.loads((OUT / "rows" / f"{nm}.json").read_text())
                        ax.imshow(X.hex_raster(z[tkey].astype(float)), cmap="Greens", vmin=0, vmax=X.BA, aspect=np.sqrt(3), interpolation="nearest")
                        if row["turing"]:
                            for s in ax.spines.values(): s.set_color("#d55e00"); s.set_linewidth(1.5)
                        if r == 0: ax.set_title(f"{f:.3g}x", fontsize=7)
                        if c == 0: ax.set_ylabel(f"n {na}/{ni}", fontsize=7)
                fig.suptitle(f"Band scan, {cname}, start {X.start_name(st)}, {tkey}, D_i 10, replicate 0; orange frame = Turing-unstable", fontsize=9)
                fig.tight_layout(rect=[0, 0, 1, 0.97])
                fig.savefig(dest / f"{cname}_{st}_{tkey}.png", dpi=80, facecolor="white"); plt.close(fig)
    print("sheets ->", dest)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "help"
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else (os.cpu_count() or 2)
    if cmd == "count":
        print("positions", len(FRACS), "jobs", len(all_jobs()))
    elif cmd == "thresholds":
        cmd_thresholds(workers)
    elif cmd == "run":
        cmd_run(workers)
    elif cmd == "aggregate":
        cmd_aggregate()
    elif cmd == "sheets":
        cmd_sheets()
    else:
        print("usage: bandscan.py count | thresholds [w] | run [w] | aggregate | sheets")
