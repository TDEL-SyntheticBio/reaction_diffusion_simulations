"""
Architecture screen (2D hex): one kinetic framework, nine ways for the two signals to travel.

Each species is a node of one type (Supplementary Fig S15 legend):
  intra       acts only in its own cell (signal = own value)
  membrane    acts only on neighbouring cells (signal = mean over the valid neighbours, self excluded)
  diffusible  acts through its field, own cell included (signal = own value), and diffuses (raw six-neighbour Laplacian,
              the engine's convention; D_a for A, D_i for B)
A is the activator node, B the inhibitor node. A x B node types = 9 circuits, e.g.
  A membrane / B diffusible = JAPI,  A diffusible / B diffusible = PAPI,  A membrane / B intra (w = 0) = lateral inhibition,
  A membrane / B membrane = JAJI,    A diffusible / B membrane = PAJI,    A intra / B diffusible = single-morphogen RD.

Production (K_a = K_i = 1, act decay mu_a = 1, inhibitor decay gamma):
  act = (1 - w) c + w (s_A)^n_a        self-activation weight w; w = 1 is the current JAPI/PAPI model, w = 0 a constant drive c
  inh = (s_B)^n_i
  promoter "shared":    H = act / (act + inh + 1);  A made at b_a H, B made at b_i H         (JAPI-like: one promoter)
  promoter "separate":  A made at b_a act / (act + inh + 1);  B made at b_i s_A^n_a / (s_A^n_a + 1)   (LI-like)
With w = 1, promoter shared, the step is written in the same arithmetic order as the engine's _vectorized_step, and
check() verifies it reproduces the engine exactly for JAPI (juxtacrine) and PAPI (paracrine).

Uniform states do not depend on node types (neighbour mean = own value = field value when uniform), so for a given
(w, promoter, kinetics) all nine circuits share the same uniform states. edge(): the largest b_i at which the uniform
system has two stable states (bistable); None if it never is (e.g. w = 0).
"""
from __future__ import annotations

import sys, json, time, os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "2D_simulations"))
sys.path.insert(0, str(ROOT / "analysis"))
from simulation_2D import _vectorized_step, _build_hex_neighbor_arrays  # noqa: E402
from hexgeom import label_components                                    # noqa: E402

DT = 0.01
BA, GA, DI, DA, C0 = 3.0, 0.5, 10.0, 1.0, 1.0
TYPES = ["intra", "membrane", "diffusible"]
_NB = {}


def nbrs(n):
    if n not in _NB:
        _NB[n] = _build_hex_neighbor_arrays(n, n)
    return _NB[n]


class Model:
    def __init__(self, A_type, B_type, w, promoter, na, ni, bi, ba=BA, ga=GA, di=DI, da=DA, c=C0, n=100):
        self.A, self.B, self.w, self.prom = A_type, B_type, float(w), promoter
        self.na, self.ni, self.bi, self.ba, self.ga, self.di, self.da, self.c = na, ni, bi, ba, ga, di, da, c
        self.NR, self.NC, self.NM, self.NCNT = nbrs(n)

    def _nsum(self, x):
        return (x[self.NR, self.NC] * self.NM).sum(axis=2)

    def signal(self, x, typ):
        return self._nsum(x) / self.NCNT if typ == "membrane" else x

    def step(self, a, b):
        sA = self.signal(a, self.A); sB = self.signal(b, self.B)
        pa = np.where(sA > 0, sA ** self.na, 0.0)
        act = pa if self.w == 1.0 else (1.0 - self.w) * self.c + self.w * pa
        inh = np.where(sB > 0, sB ** self.ni, 0.0)
        h = act / (act + inh + 1.0)
        hb = h if self.prom == "shared" else pa / (pa + 1.0)
        ra = DT * (self.ba * h - 1.0 * a)
        rb = DT * (self.bi * hb - self.ga * b)
        da_ = self.da * DT / 1.0 ** 2 * (self._nsum(a) - self.NCNT * a) if self.A == "diffusible" else 0.0
        db_ = self.di * DT / 1.0 ** 2 * (self._nsum(b) - self.NCNT * b) if self.B == "diffusible" else 0.0
        return a + ra + da_, b + rb + db_


# ------------------------------------------------------------------ uniform states and the edge
def _uniform_rhs(a, b, w, prom, na, ni, bi, ba=BA, ga=GA, c=C0):
    pa = a ** na if a > 0 else 0.0
    act = pa if w == 1.0 else (1 - w) * c + w * pa
    inh = b ** ni if b > 0 else 0.0
    h = act / (act + inh + 1)
    hb = h if prom == "shared" else pa / (pa + 1)
    return ba * h - a, bi * hb - ga * b


def uniform_states(w, prom, na, ni, bi, ba=BA, ga=GA, c=C0):
    """All uniform steady states (a, b, stable) found by scanning a; b follows from the B equation."""
    def b_of(a):
        if prom == "shared":
            return (bi / (ga * ba)) * a            # b = (b_i/gamma) H and a = b_a H
        pa = a ** na
        return (bi / ga) * pa / (pa + 1)
    A = np.concatenate([[0.0], np.logspace(-6, np.log10(ba), 6000)])
    g = np.array([_uniform_rhs(x, b_of(x), w, prom, na, ni, bi, ba, ga, c)[0] for x in A])
    roots = []
    for k in range(len(A) - 1):
        if g[k] == 0:
            roots.append(A[k])
        elif np.sign(g[k]) != np.sign(g[k + 1]):
            lo, hi = A[k], A[k + 1]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                gm = _uniform_rhs(mid, b_of(mid), w, prom, na, ni, bi, ba, ga, c)[0]
                if np.sign(gm) == np.sign(g[k]): lo = mid
                else: hi = mid
            roots.append(0.5 * (lo + hi))
    out = []
    for a in roots:
        b = b_of(a); e = 1e-7 * max(1.0, a, b)
        f0 = np.array(_uniform_rhs(a, b, w, prom, na, ni, bi, ba, ga, c))
        J = np.column_stack([(np.array(_uniform_rhs(a + e, b, w, prom, na, ni, bi, ba, ga, c)) - f0) / e,
                             (np.array(_uniform_rhs(a, b + e, w, prom, na, ni, bi, ba, ga, c)) - f0) / e])
        tr, det = np.trace(J), np.linalg.det(J)
        out.append((float(a), float(b), bool(tr < 0 and det > 0)))
    return out


def n_stable(w, prom, na, ni, bi):
    return sum(s for _, _, s in uniform_states(w, prom, na, ni, bi))


def edge(w, prom, na, ni):
    """Largest b_i with two stable uniform states (bistable), by log scan then bisection; None if never bistable."""
    grid = np.logspace(-3, 3, 241)
    bis = [n_stable(w, prom, na, ni, x) >= 2 for x in grid]
    if not any(bis):
        return None
    k = max(i for i, v in enumerate(bis) if v)
    if k == len(grid) - 1:
        return float(grid[-1])
    lo, hi = grid[k], grid[k + 1]
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        if n_stable(w, prom, na, ni, mid) >= 2: lo = mid
        else: hi = mid
    return float(lo)


# ------------------------------------------------------------------ check against the deposited engine
def check(n=60, steps=400):
    rng = np.random.default_rng(1)
    a0 = rng.random((n, n)) * 3; b0 = rng.random((n, n))
    NR, NC, NM, NCNT = nbrs(n)
    res = {}
    for name, A_type, atype, dact in (("JAPI", "membrane", "juxtacrine", 0.0), ("PAPI", "diffusible", "paracrine", DA)):
        p = dict(act_half_sat=1.0, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=dact,
                 inh_diffusion=DI, act_prod_rate=BA, inh_prod_rate=1.7, inh_decay_rate=GA, act_hill_coeff=4, inh_hill_coeff=4)
        m = Model(A_type, "diffusible", 1.0, "shared", 4, 4, 1.7, n=n)
        a1, b1, a2, b2 = a0.copy(), b0.copy(), a0.copy(), b0.copy()
        for _ in range(steps):
            a1, b1 = _vectorized_step(a1, b1, DT, 1.0, p, atype, NR, NC, NM, NCNT)
            a2, b2 = m.step(a2, b2)
        res[name] = (float(np.abs(a1 - a2).max()), float(np.abs(b1 - b2).max()))
    return res



# ------------------------------------------------------------------ design
HILL = [(4, 4), (10, 4)]
W = [1.0, 0.25, 0.0]
PROMOTERS = ["shared", "separate"]
FRACS = [0.4, 0.7, 1.0, 1.4, 2.0, 3.0]
STARTS = ["noise", "seeds", "nuc001"]
CIRCUITS = [(A, B) for A in TYPES for B in TYPES]
SEED_DENSITY = 70 / 40000
DAY = 50 / 11.25
KICK_MULT, KICK_CAP, KICK_LO = 1.5, 100.0, 1.0
OUT = Path(os.environ.get("ARCH_OUT", ROOT / "results" / "archscreen"))
N_DEFAULT = int(os.environ.get("ARCH_N", 100))


def ref_edge(na, ni):
    return edge(1.0, "shared", na, ni)


def positions(w, prom, na, ni):
    """Union of the reference list (multiples of the shared w = 1 edge) and the version's own edge list, if it has one."""
    E_ref = ref_edge(na, ni); E_own = edge(w, prom, na, ni)
    pts = [("ref", f, f * E_ref) for f in FRACS]
    if E_own is not None and abs(E_own - E_ref) > 1e-6 * E_ref:
        pts += [("own", f, f * E_own) for f in FRACS]
    return pts, E_ref, E_own


# ------------------------------------------------------------------ ignition threshold per configuration
def ignites(A_type, B_type, w, prom, na, ni, bi, h, n=60, t=20.0):
    """A single cell at the centre set to h (everything else starting at zero) ignites if, at time t, at least 7 cells
    differ from the same run without the kick by more than 0.3 b_a."""
    m = Model(A_type, B_type, w, prom, na, ni, bi, n=n)
    a = np.zeros((n, n)); b = np.zeros((n, n)); a0 = a.copy(); b0 = b.copy()
    a[n // 2, n // 2] = h
    for _ in range(int(round(t / DT))):
        a, b = m.step(a, b); a0, b0 = m.step(a0, b0)
    return int((np.abs(a - a0) > 0.3 * BA).sum()) >= 7


def threshold(cfg):
    """Lowest single-cell kick that ignites, by log bisection on [1, 100] (7 steps); measured at 0.7x the version's own
    edge, or 0.7x the reference edge when the version has no edge. Returns (threshold or None, b_i used, note)."""
    A_type, B_type, w, prom, na, ni = cfg
    E = edge(w, prom, na, ni) or ref_edge(na, ni)
    bi = 0.7 * E
    if ignites(*cfg, bi, KICK_LO):
        return KICK_LO, bi, "ignites at 1"
    if not ignites(*cfg, bi, KICK_CAP):
        return None, bi, "no ignition up to 100"
    lo, hi = KICK_LO, KICK_CAP
    for _ in range(7):
        mid = np.sqrt(lo * hi)
        if ignites(*cfg, bi, mid): hi = mid
        else: lo = mid
    return float(hi), bi, ""


def kick_height(th):
    return KICK_CAP if th is None else float(min(KICK_MULT * th, KICK_CAP))


# ------------------------------------------------------------------ starts and measurements
def stream(start, rep):
    return np.random.default_rng(np.random.SeedSequence([20261022, STARTS.index(start), rep]))


def initial(m, start, rng, n, kick):
    a = np.zeros((n, n)); b = np.zeros((n, n))
    if start == "noise":
        st = uniform_states(m.w, m.prom, m.na, m.ni, m.bi)
        stab = [s for s in st if s[2]] or st
        a0, b0, _ = max(stab, key=lambda s: s[0])
        xi = rng.uniform(-1, 1, (2, n, n))
        a = a0 * (1 + 0.05 * xi[0]) + 0.01 * rng.random((n, n))
        b = b0 * (1 + 0.05 * xi[1])
    elif start == "seeds":
        k = max(1, int(round(SEED_DENSITY * n * n)))
        idx = rng.choice(n * n, size=k, replace=False)
        a.flat[idx] = kick
    return a, b


def square_raster(a):
    """Hex field (even-r) resampled on a square grid of 0.5-cell pixels in physical coordinates (nearest cell)."""
    n = a.shape[0]
    ys = np.arange(0, (n - 1) * np.sqrt(3) / 2 + 1e-9, 0.5); xs = np.arange(0, n, 0.5)
    r = np.clip(np.rint(ys / (np.sqrt(3) / 2)).astype(int), 0, n - 1)
    shift = (r % 2 == 0) * 0.5
    c = np.clip(np.rint(xs[None, :] - shift[:, None]).astype(int), 0, n - 1)
    return a[r[:, None], c]


def wavelength(a):
    img = square_raster(a); img = img - img.mean()
    if img.std() < 1e-9:
        return np.nan
    P = np.abs(np.fft.fft2(img)) ** 2
    ky = np.fft.fftfreq(img.shape[0], d=0.5); kx = np.fft.fftfreq(img.shape[1], d=0.5)
    K = np.sqrt(ky[:, None] ** 2 + kx[None, :] ** 2)
    kmin = 2.0 / (0.5 * min(img.shape))
    bins = np.linspace(kmin, K.max(), 120)
    idx = np.digitize(K.ravel(), bins)
    prof = np.bincount(idx, weights=P.ravel(), minlength=len(bins) + 1)[1:len(bins)]
    cnt = np.bincount(idx, minlength=len(bins) + 1)[1:len(bins)]
    prof = np.where(cnt > 0, prof / np.maximum(cnt, 1), 0)
    j = int(np.argmax(prof)); kc = 0.5 * (bins[j] + bins[j + 1])
    return float(1.0 / kc)


def moran(m, a):
    """Moran's I with the engine's hex adjacency: about 0 for independent cells, negative for salt-and-pepper,
    positive for domains."""
    z = a - a.mean(); den = (z * z).sum()
    if den <= 1e-18:
        return np.nan
    return float((a.size / m.NM.sum()) * (z * m._nsum(z)).sum() / den)


def measure(m, a, a_prev=None):
    """Raw numbers kept for every frame. The label is a working aid only and can be recomputed from the rows."""
    thr = 0.3 * BA
    on = a > thr
    lab, k = label_components(on)
    sizes = np.bincount(lab[lab >= 0]) if k else np.array([], int)
    big7 = [j for j in range(len(sizes)) if sizes[j] >= 7]
    persist = None
    if a_prev is not None:
        lp, kp = label_components(a_prev > thr)
        prev7 = np.zeros(a.shape, bool)
        if kp:
            sp = np.bincount(lp[lp >= 0])
            for j in np.nonzero(sp >= 7)[0]:
                prev7 |= lp == j
        persist = int(sum(1 for j in big7 if ((lab == j) & prev7).any()))
    span = False
    if big7:
        jb = int(np.argmax(sizes)); mk = lab == jb
        span = bool((mk[0].any() and mk[-1].any()) or (mk[:, 0].any() and mk[:, -1].any()))
    spread = float(a.std() / BA); wl = wavelength(a); mi = moran(m, a); cov = float(on.mean())
    n_dom = persist if persist is not None else len(big7)
    if spread < 0.02:
        label = "uniform"
    elif mi is not None and mi < -0.2:
        label = "fine (salt-and-pepper)"
    elif abs(mi) < 0.2:
        label = "speckle (cells independent)"
    elif cov >= 0.9:
        label = "sheet"
    elif n_dom == 0:
        label = "no persistent domains"
    elif span:
        label = "network"
    else:
        label = "spots"
    return dict(label=label, coverage=cov, mean_a=float(a.mean()), spread=spread, moran=mi, wavelength=wl,
                n_comp=int(k), n_comp7=len(big7), n_persist7=persist,
                size_median7=float(np.median(sizes[big7])) if big7 else np.nan,
                largest_frac=float(sizes.max() / max(on.sum(), 1)) if k else 0.0, spans=span)


def job_name(job):
    A_type, B_type, w, prom, na, ni, kind, f, bi, start, rep = job
    return f"{A_type[:4]}_{B_type[:4]}_w{w:g}_{prom[:3]}_na{na:02d}_ni{ni:02d}_{kind}{f:g}_{start}_r{rep}"


def run(job, kick, n=None, out=None, t_end=50.0, save=True):
    n = n or N_DEFAULT
    A_type, B_type, w, prom, na, ni, kind, f, bi, start, rep = job
    m = Model(A_type, B_type, w, prom, na, ni, bi, n=n)
    rng = stream(start, rep)
    a, b = initial(m, start, rng, n, kick)
    p_fire = 0.001 * DT
    steps = int(round(t_end / DT))
    marks = {int(round(4 * DAY / DT)) - 500: "day4_prev", int(round(4 * DAY / DT)): "day4", steps - 500: "t50_prev", steps: "t50"}
    ts = dict(t=[], coverage=[], n_comp7=[], mean_a=[]); keep = {}; t0 = time.time()
    for s in range(1, steps + 1):
        a, b = m.step(a, b)
        if start == "nuc001":
            fire = (rng.random((n, n)) < p_fire) & (a < BA)
            a[fire] = kick
        if s % 100 == 0:
            on = a > 0.3 * BA; lab, k = label_components(on)
            ts["t"].append(round(s * DT, 2)); ts["coverage"].append(float(on.mean())); ts["mean_a"].append(float(a.mean()))
            ts["n_comp7"].append(int((np.bincount(lab[lab >= 0]) >= 7).sum()) if k else 0)
        if s in marks:
            keep[marks[s]] = a.copy()
    md4 = measure(m, keep["day4"], keep["day4_prev"]); m50 = measure(m, keep["t50"], keep["t50_prev"])
    name = job_name(job)
    row = dict(name=name, A=A_type, B=B_type, circuit=f"{A_type}/{B_type}", w=w, promoter=prom, n_a=na, n_i=ni,
               position_kind=kind, f=f, b_i=bi, start=start, rep=rep, N=n, kick=kick,
               **{f"{k}_t50": v for k, v in m50.items()}, **{f"{k}_day4": v for k, v in md4.items()},
               b_min=float(b.min()), b_max=float(b.max()), a_max=float(a.max()),
               finite=bool(np.isfinite(a).all() and np.isfinite(b).all()), seconds=time.time() - t0)
    if out is not None and save:
        out = Path(out); (out / "rows").mkdir(parents=True, exist_ok=True); (out / "frames").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out / "frames" / f"{name}.npz", day4=keep["day4"].astype(np.float16), t50=keep["t50"].astype(np.float16),
                            b_t50=b.astype(np.float16), **{k: np.array(v) for k, v in ts.items()})
        (out / "rows" / f"{name}.json").write_text(json.dumps(row))
    return row


# ------------------------------------------------------------------ driver
CONFIGS = [(A, B, w, prom, na, ni) for (na, ni) in HILL for w in W for prom in PROMOTERS for (A, B) in CIRCUITS]
CHECK_CIRCUITS = [("membrane", "diffusible", 1.0, "shared"), ("diffusible", "diffusible", 1.0, "shared"),
                  ("intra", "diffusible", 1.0, "shared"), ("membrane", "intra", 0.0, "separate"),
                  ("membrane", "membrane", 1.0, "shared"), ("diffusible", "membrane", 1.0, "shared")]


def _thr_job(cfg):
    th, bi, note = threshold(cfg)
    return dict(A=cfg[0], B=cfg[1], w=cfg[2], promoter=cfg[3], n_a=cfg[4], n_i=cfg[5], threshold=th, b_i_used=bi,
                note=note, kick=kick_height(th))


def cmd_thresholds(workers):
    from multiprocessing import Pool
    f = OUT / "thresholds.json"
    if f.exists():
        print("thresholds already computed:", f); return
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with Pool(workers) as p:
        rows = p.map(_thr_job, CONFIGS, chunksize=1)
    f.write_text(json.dumps(rows, indent=1))
    print(f"{len(rows)} thresholds in {time.time() - t0:.0f} s -> {f}", flush=True)


def kicks():
    rows = json.loads((OUT / "thresholds.json").read_text())
    return {(r["A"], r["B"], r["w"], r["promoter"], r["n_a"], r["n_i"]): r["kick"] for r in rows}


def all_jobs(rep=0):
    jobs = []
    for (na, ni) in HILL:
        for w in W:
            for prom in PROMOTERS:
                pts, _, _ = positions(w, prom, na, ni)
                for (A, B) in CIRCUITS:
                    for kind, f, bi in pts:
                        for st in STARTS:
                            jobs.append((A, B, w, prom, na, ni, kind, f, bi, st, rep))
    return jobs


def check_jobs():
    jobs = []
    for (na, ni) in HILL:
        E = ref_edge(na, ni)
        for (A, B, w, prom) in CHECK_CIRCUITS:
            for f in (1.0, 1.4):
                for st in STARTS:
                    jobs.append((A, B, w, prom, na, ni, "ref", f, f * E, st, 0))
    return jobs


def _run_job(args):
    job, kick, n, out = args
    try:
        return run(job, kick, n=n, out=out)
    except Exception as e:  # keep going, record the failure
        (Path(out) / "failed").mkdir(parents=True, exist_ok=True)
        (Path(out) / "failed" / f"{job_name(job)}.txt").write_text(repr(e))
        return None


def cmd_run(workers, which="main"):
    from multiprocessing import Pool
    K = kicks()
    if which == "main":
        jobs, n, out = all_jobs(), N_DEFAULT, OUT
    else:
        jobs, n, out = check_jobs(), 200, OUT / "check200"
    todo = [(j, K[j[:6]], n, out) for j in jobs if not (Path(out) / "rows" / f"{job_name(j)}.json").exists()]
    print(f"{which}: {len(jobs)} jobs, {len(todo)} to run, N = {n}, {workers} workers", flush=True)
    t0 = time.time()
    with Pool(workers) as p:
        for k, r in enumerate(p.imap_unordered(_run_job, todo, chunksize=1), 1):
            if k % 50 == 0 or k == len(todo):
                el = time.time() - t0
                print(f"  {k}/{len(todo)}  {el/3600:.2f} h elapsed, about {el/k*(len(todo)-k)/3600:.2f} h left", flush=True)


def cmd_aggregate():
    import pandas as pd
    for sub, name in (("", "archscreen_runs.csv"), ("check200", "archscreen_check200_runs.csv")):
        rows = [json.loads(f.read_text()) for f in sorted((OUT / sub / "rows").glob("*.json"))] if (OUT / sub / "rows").exists() else []
        if rows:
            dest = ROOT / "analysis" / "summaries" / name
            pd.DataFrame(rows).to_csv(dest, index=False); print(f"{len(rows)} rows -> {dest}")
    thr = OUT / "thresholds.json"
    if thr.exists():
        dest = ROOT / "analysis" / "summaries" / "archscreen_thresholds.csv"
        pd.DataFrame(json.loads(thr.read_text())).to_csv(dest, index=False); print("thresholds ->", dest)
    edges = [dict(n_a=na, n_i=ni, w=w, promoter=p, edge=edge(w, p, na, ni), ref_edge=ref_edge(na, ni))
             for (na, ni) in HILL for w in W for p in PROMOTERS]
    pd.DataFrame(edges).to_csv(ROOT / "analysis" / "summaries" / "archscreen_edges.csv", index=False)


CIRCUIT_NAME = {("intra", "intra"): "no coupling", ("intra", "membrane"): "A intra / B membrane",
                ("intra", "diffusible"): "single-morphogen RD", ("membrane", "intra"): "LI (w 0) / JA + autonomous B",
                ("membrane", "membrane"): "JAJI", ("membrane", "diffusible"): "JAPI",
                ("diffusible", "intra"): "PA + autonomous B", ("diffusible", "membrane"): "PAJI", ("diffusible", "diffusible"): "PAPI"}


def hex_raster(field):
    ny, nx = field.shape; img = np.full((ny, 2 * nx + 1), np.nan)
    for r in range(ny):
        o = 1 if r % 2 == 0 else 0; img[r, o:o + 2 * nx] = np.repeat(field[r], 2)
    return img


def cmd_sheets():
    """One PNG per (Hill pair, promoter, w, start, time): rows = positions, columns = the 9 circuits; colour 0..b_a fixed."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    dest = ROOT / "analysis" / "summaries" / "archscreen_sheets"; dest.mkdir(parents=True, exist_ok=True)
    for (na, ni) in HILL:
        for w in W:
            for prom in PROMOTERS:
                pts, E_ref, E_own = positions(w, prom, na, ni)
                for st in STARTS:
                    for tkey in ("t50", "day4"):
                        fig, axs = plt.subplots(len(pts), 9, figsize=(18, 2.1 * len(pts) + 0.8), squeeze=False)
                        for r, (kind, f, bi) in enumerate(pts):
                            for c, (A, B) in enumerate(CIRCUITS):
                                ax = axs[r, c]; ax.set_xticks([]); ax.set_yticks([])
                                job = (A, B, w, prom, na, ni, kind, f, bi, st, 0); nm = job_name(job)
                                fp = OUT / "frames" / f"{nm}.npz"; rp = OUT / "rows" / f"{nm}.json"
                                if not fp.exists():
                                    ax.set_facecolor("#eeeeee"); ax.set_title("missing", fontsize=6); continue
                                z = np.load(fp); row = json.loads(rp.read_text())
                                ax.imshow(hex_raster(z[tkey].astype(float)), cmap="Greens", vmin=0, vmax=BA, aspect=np.sqrt(3), interpolation="nearest")
                                ax.set_title(f"{row['label_' + tkey]}\nI {row['moran_' + tkey]:.2f}, wl {row['wavelength_' + tkey]:.1f}", fontsize=6)
                                if r == 0: ax.text(0.5, 1.45, CIRCUIT_NAME[(A, B)], transform=ax.transAxes, ha="center", fontsize=7.5, weight="bold")
                                if c == 0: ax.set_ylabel(f"{f:g}x {'ref' if kind == 'ref' else 'own'} edge\nb_i {bi:.3g}", fontsize=7)
                        own = "none" if E_own is None else f"{E_own:.3g}"
                        fig.suptitle(f"Architecture screen, n {na}/{ni}, promoter {prom}, w {w:g}, start {start_name(st)}, "
                                     f"{'t = 50' if tkey == 't50' else 'day 4'}; ref edge {E_ref:.3g}, own edge {own}; colour 0..b_a; "
                                     f"title: label, Moran's I, wavelength (cells)", fontsize=9)
                        fig.tight_layout(rect=[0, 0, 1, 0.97])
                        fig.savefig(dest / f"{tkey}_na{na}_ni{ni}_{prom}_w{w:g}_{st}.png", dpi=80, facecolor="white"); plt.close(fig)
    print("sheets ->", dest)


def start_name(st):
    return {"noise": "noise around uniform state", "seeds": "scattered seeds", "nuc001": "nucleation 0.001"}[st]


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "help"
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else (os.cpu_count() or 2)
    if cmd == "check":
        print(check())
    elif cmd == "thresholds":
        cmd_thresholds(workers)
    elif cmd == "run":
        cmd_run(workers, "main")
    elif cmd == "check200":
        cmd_run(workers, "check200")
    elif cmd == "aggregate":
        cmd_aggregate()
    elif cmd == "sheets":
        cmd_sheets()
    elif cmd == "count":
        print("configs", len(CONFIGS), "main jobs", len(all_jobs()), "check200 jobs", len(check_jobs()))
    else:
        print("usage: archscreen.py check | thresholds [w] | run [w] | check200 [w] | aggregate | sheets | count")
