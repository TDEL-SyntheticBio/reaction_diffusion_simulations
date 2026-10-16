"""
1D single-spike arrest map over the full (n_a, n_i) plane, beside the analytic state / Turing audit.  1D only: the
1D_simuations/simulation.py engine (juxtacrine signal = mean of the two neighbours, Neumann Laplacian with prefactor 1,
inhibitor range sqrt(D_i / gamma)); nothing here uses 2D results or conventions.

Grid (the committed 1D_batch/config.yaml sweep): b_a, b_i = linspace(1, 10, 10); gamma = linspace(0.1, 1, 10) + linspace(2, 10, 9)
(19 values); n_a, n_i = 1..12; 273,600 sets.  Fixed: N = 200, steps = 50,000, dt = 0.01, dx = 1, save_every = 200,
stopping_threshold = 1e-4, min_steps = 500, juxtacrine, act_decay = act_half_sat = inh_half_sat = 1, basal = 0, D_i = 10,
act_diffusion = 1 (ignored under juxtacrine), spike 10 at cell N // 2 (the deposit's activator_spike mode), inhibitor zero.

Integration: a replica of run_coupled_neumann's loop (same expression order, same stopping rule: every save_every steps the
L1 change per cell-field since the previous save is compared with the threshold once step > min_steps), vectorised over
parameter sets within one (n_a, n_i) cell, returning the state at the true final step (run_simulation returns a frame
save_every - 1 steps stale on a capped run).  `check` verifies the replica bit-exact against run_simulation.  The deposited
fast_stable_steady_state is never called.

Audit: the uniform fixed point a = b_a H, i = (b_i / gamma) H, H = A / (A + I + 1), scanned and bisected in H for every set
(largest root above 1e-5), stability from the 2x2 reaction Jacobian (trace < 0 and det > 0), Turing from the dispersion
relation over the accessible lattice modes q = pi m / N, m = 0..N-1 (K = cos q, Lambda = 2 (1 - cos q)): a mode is unstable
iff det < 0 or trace > 0; Turing = stable at q = 0 and unstable at some q > 0.  No closed-form inequality is used.

Labels from contiguous domains of the final activator above 0.3 b_a: 0 domains "all off"; 1 domain spanning >= 0.95 N cells
"all on"; 1 domain "one domain"; >= 2 "multiple domains".  No peak finding.

Usage:
  python analysis/spike_grid.py audit [N]                 # analysis/summaries/state_audit_turing.csv (default N = 200)
  python analysis/spike_grid.py check                     # bit-exact replica check (30 runs) and the n_a 10 / n_i 4 cross-check at N = 100
  python analysis/spike_grid.py time                      # 30 single runs at N = 200 across the plane + one batched cell; projection
  python analysis/spike_grid.py run [workers] [--gamma10] # the full grid, one checkpoint per (n_a, n_i) cell in results/spike_grid/
  python analysis/spike_grid.py validate                  # ~2,000-set random_tight validation of the audit (results in spike_grid_validation.csv)
  python analysis/spike_grid.py spike2                    # spike = 2.0 subsample over the plane (spike_height_2.csv)
  python analysis/spike_grid.py aggregate                 # spike_grid_full.csv and spike_grid_by_hill.csv from the checkpoints
"""
from __future__ import annotations

import csv, io, contextlib, json, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SIM1D = ROOT / "1D_simuations"
sys.path.insert(0, str(SIM1D)); sys.path.insert(0, str(SIM1D / "1D_batch"))
SUMM = ROOT / "analysis" / "summaries"
RES = ROOT / "results" / "spike_grid"

N_DEFAULT, STEPS, DT, DX, SAVE_EVERY, STOP_THR, MIN_STEPS = 200, 50000, 0.01, 1.0, 200, 1e-4, 500
D_I, SPIKE = 10.0, 10.0
BA_VALUES = np.linspace(1, 10, 10); BI_VALUES = np.linspace(1, 10, 10)
GAMMA_19 = np.concatenate([np.linspace(0.1, 1.0, 10), np.linspace(2.0, 10.0, 9)])
GAMMA_10 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0, 5.0, 10.0])
HILL = np.arange(1, 13)
LABELS = ("all off", "all on", "one domain", "multiple domains")


def base_params(ba, bi, g, na, ni):
    """The committed 1D_batch base dict with the swept keys as floats, as make_param_grid builds them."""
    return dict(N=N_DEFAULT, steps=STEPS, dt=DT, dx=DX, save_every=SAVE_EVERY, stopping_threshold=STOP_THR, min_steps=MIN_STEPS,
                init_mode="activator_spike", activator_type="juxtacrine", spike_value=SPIKE, act_decay_rate=1.0, act_half_sat=1.0,
                inh_half_sat=1.0, basal_prod=0.0, act_hill_coeff=float(na), inh_hill_coeff=float(ni), act_prod_rate=float(ba),
                inh_prod_rate=float(bi), inh_decay_rate=float(g), act_diffusion=1.0, inh_diffusion=D_I)


# ------------------------------------------------------------------------------------------------ the audit
def fixed_point_scan(ba, bi, g, na, ni, n_uniform=20001):
    """Largest root H* > 1e-5 of H = A / (A + I + 1), A = (b_a H)^n_a, I = (b_i H / gamma)^n_i, for arrays of sets sharing (n_a, n_i).
    Scan on a uniform H grid plus a geometric approach to 1 (saturated states), then 60 bisections. Returns H* (nan if none) and the root count."""
    ba = np.asarray(ba, float)[:, None]; bi = np.asarray(bi, float)[:, None]; g = np.asarray(g, float)[:, None]
    Hs = np.unique(np.concatenate([np.linspace(1e-5, 1 - 1e-6, n_uniform), 1 - np.logspace(-6, -16, 1001)]))[None, :]

    def f(H):
        A = (ba * H) ** na; I = (bi * H / g) ** ni
        return A / (A + I + 1.0) - H
    fv = f(Hs)
    sgn = np.sign(fv)
    change = (sgn[:, :-1] * sgn[:, 1:] < 0) | (fv[:, 1:] == 0)
    n_roots = change.sum(1)
    Hstar = np.full(ba.shape[0], np.nan)
    has = n_roots > 0
    if has.any():
        last = np.array([np.nonzero(c)[0][-1] if c.any() else 0 for c in change])
        lo = Hs[0, last].copy(); hi = Hs[0, last + 1].copy()
        rows = np.arange(ba.shape[0])
        bal, bil, gl = ba[:, 0], bi[:, 0], g[:, 0]

        def f1(H):
            A = (bal * H) ** na; I = (bil * H / gl) ** ni
            return A / (A + I + 1.0) - H
        flo = f1(lo)
        for _ in range(60):
            mid = 0.5 * (lo + hi); fm = f1(mid)
            left = flo * fm <= 0
            hi = np.where(left, mid, hi); lo = np.where(left, lo, mid); flo = np.where(left, flo, fm)
        Hstar[has] = (0.5 * (lo + hi))[has]
    return Hstar, n_roots


def jacobian_terms(ba, bi, g, na, ni, H):
    a = ba * H; i = bi * H / g
    A = a ** na; I = i ** ni; S = A + I + 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        hA = na * A / a * (I + 1.0) / S ** 2
        hI = -ni * I / i * A / S ** 2
    return a, i, hA, hI


def audit_cell(na, ni, N, gammas=GAMMA_19):
    BA, BI, G = np.meshgrid(BA_VALUES, BI_VALUES, gammas, indexing="ij")
    ba, bi, g = BA.ravel(), BI.ravel(), G.ravel()
    H, n_roots = fixed_point_scan(ba, bi, g, na, ni)
    a, i, hA, hI = jacobian_terms(ba, bi, g, na, ni, H)
    j00 = ba * hA - 1.0; j01 = ba * hI; j10 = bi * hA; j11 = bi * hI - g
    tr = j00 + j11; det = j00 * j11 - j01 * j10
    stable = np.isfinite(H) & (tr < 0) & (det > 0)
    q = np.pi * np.arange(N) / N
    K = np.cos(q)[None, :]; Lam = 2.0 * (1.0 - np.cos(q))[None, :]
    J00 = ba[:, None] * hA[:, None] * K - 1.0; J01 = (ba * hI)[:, None]; J10 = bi[:, None] * hA[:, None] * K
    J11 = (bi * hI - g)[:, None] - D_I * Lam
    TR = J00 + J11; DET = J00 * J11 - J01 * J10
    unstable = (DET < 0) | (TR > 0)
    growth = np.where(unstable, 0.5 * (TR + np.sqrt(np.maximum(TR ** 2 - 4 * DET, 0.0))), -np.inf)      # largest real part where unstable
    turing = stable & unstable[:, 1:].any(1)
    q_fast = np.where(turing, q[np.argmax(np.where(unstable, growth, -np.inf), axis=1)], np.nan)
    sigma_max = np.where(turing, np.max(np.where(unstable, growth, -np.inf), axis=1), np.nan)      # largest linear growth rate over the unstable modes
    q_lo = np.full(len(ba), np.nan); q_hi = np.full(len(ba), np.nan)
    for r in np.nonzero(turing)[0]:
        m = np.nonzero(unstable[r, 1:])[0] + 1
        q_lo[r], q_hi[r] = q[m[0]], q[m[-1]]
    return dict(b_a=ba, b_i=bi, gamma=g, n_a=np.full(len(ba), na), n_i=np.full(len(ba), ni), truth_H=H, truth_a_star=a, truth_i_star=i,
                n_roots=n_roots, trace0=tr, det0=det, truth_has_state=stable, turing=turing, q_fastest=q_fast, sigma_max=sigma_max,
                wavelength_fastest=np.where(turing, 2 * np.pi / np.where(q_fast > 0, q_fast, np.nan), np.nan), q_unstable_min=q_lo, q_unstable_max=q_hi,
                N_modes=np.full(len(ba), N))


def audit(N=N_DEFAULT, out=None):
    rows = []
    for na in HILL:
        for ni in HILL:
            rows.append(audit_cell(int(na), int(ni), N))
    import pandas as pd
    df = pd.concat([pd.DataFrame(r) for r in rows], ignore_index=True)
    out = out or SUMM / "state_audit_turing.csv"
    df.to_csv(out, index=False, lineterminator="\n")
    print(f"{out}: {len(df)} sets; has_state {int(df.truth_has_state.sum())}, turing {int(df.turing.sum())}")
    return df


def load_audit(N=N_DEFAULT):
    import pandas as pd
    f = SUMM / "state_audit_turing.csv"
    if f.exists():
        df = pd.read_csv(f)
        if "N_modes" not in df or int(df.N_modes.iloc[0]) != N:
            print(f"audit on disk uses N = {df.N_modes.iloc[0] if 'N_modes' in df else '?'}; recomputing at N = {N}")
            df = audit(N)
        return df
    return audit(N)


# ------------------------------------------------------------------------------------------------ the integrator
def integrate_batch(ba, bi, g, na, ni, N, init_a, init_i, steps=STEPS, dt=DT, dx=DX, save_every=SAVE_EVERY, thr=STOP_THR, min_steps=MIN_STEPS):
    """Replica of run_coupled_neumann for M parameter sets at once (rows), same arithmetic and stopping rule per row; returns the final
    activator and inhibitor (true final step) and the engine's `step` value per row (the loop index at which it stopped, or steps - 1)."""
    M = len(ba)
    ba = np.asarray(ba, float)[:, None]; bi = np.asarray(bi, float)[:, None]; g = np.asarray(g, float)[:, None]
    na = float(na); ni = float(ni)
    a = np.array(init_a, float, copy=True); i = np.array(init_i, float, copy=True)
    assert a.shape == (M, N)
    out_a = np.empty_like(a); out_i = np.empty_like(i); out_step = np.full(M, steps - 1)
    active = np.arange(M)
    last_a = a.copy(); last_i = i.copy()
    coef = D_I * dt / dx ** 2
    sig = np.empty_like(a); lap = np.empty_like(i)
    for step in range(steps):
        sig[:, 1:-1] = 0.5 * (a[:, :-2] + a[:, 2:]); sig[:, 0] = a[:, 1]; sig[:, -1] = a[:, -2]
        lap[:, 1:-1] = i[:, 2:] - 2.0 * i[:, 1:-1] + i[:, :-2]; lap[:, 0] = i[:, 1] - i[:, 0]; lap[:, -1] = i[:, -2] - i[:, -1]
        act_term = np.power(np.maximum(sig, 0.0) / 1.0, na); inh_term = np.power(np.maximum(i, 0.0) / 1.0, ni)
        hill = (act_term + 0.0) / (act_term + inh_term + 1.0 + 0.0)
        a_new = (a + dt * (ba * hill - 1.0 * a)) + 0.0
        i_new = (i + dt * (bi * hill - g * i)) + coef * lap
        a, i = a_new, i_new
        if step % save_every == 0:
            diff = np.abs(a - last_a).sum(1) + np.abs(i - last_i).sum(1)
            last_a = a.copy(); last_i = i.copy()
            if step > min_steps:
                done = diff / (2 * N) < thr
                if done.any():
                    out_a[active[done]] = a[done]; out_i[active[done]] = i[done]; out_step[active[done]] = step
                    keep = ~done
                    if not keep.any():
                        return out_a, out_i, out_step
                    a, i, last_a, last_i = a[keep], i[keep], last_a[keep], last_i[keep]
                    ba, bi, g = ba[keep], bi[keep], g[keep]; active = active[keep]
                    sig = np.empty_like(a); lap = np.empty_like(i)
    out_a[active] = a; out_i[active] = i
    return out_a, out_i, out_step


def spike_init(M, N, spike=SPIKE):
    a = np.zeros((M, N)); a[:, N // 2] = spike
    return a, np.zeros((M, N))


def tight_init(a_star, i_star, N, frac=0.05, seed=2):
    """The engine's random_tight draw (seed 2, the same stream for every set), scaled per set: uniform within +-frac of the state."""
    rng = np.random.default_rng(seed)
    ra = rng.random(N); ri = rng.random(N)
    a_star = np.asarray(a_star, float)[:, None]; i_star = np.asarray(i_star, float)[:, None]
    lo, hi = (1 - frac) * a_star, (1 + frac) * a_star
    a0 = lo + (hi - lo) * ra[None, :]
    lo, hi = (1 - frac) * i_star, (1 + frac) * i_star
    i0 = lo + (hi - lo) * ri[None, :]
    return a0, i0


def run_from_centre(p):
    """Single-run replica (same loop as the batch, one row); for the bit-exact check and timing."""
    N = int(p["N"])
    a0, i0 = spike_init(1, N, p["spike_value"])
    a, i, step = integrate_batch([p["act_prod_rate"]], [p["inh_prod_rate"]], [p["inh_decay_rate"]], p["act_hill_coeff"], p["inh_hill_coeff"], N, a0, i0,
                                 steps=int(p["steps"]), dt=p["dt"], dx=p["dx"], save_every=int(p["save_every"]), thr=p["stopping_threshold"], min_steps=int(p["min_steps"]))
    return a[0], i[0], int(step[0])


def load_analyze_pattern():
    """The deposited analyze_pattern, pulled out of 1D_batch/analyze_patterns.py by source (that module runs its CSV script at
    import time, so it cannot be imported); the function body is executed verbatim with its own imports."""
    import ast
    src = (SIM1D / "1D_batch" / "analyze_patterns.py").read_text()
    tree = ast.parse(src)
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "analyze_pattern")
    ns = {}
    exec("import numpy as np\nfrom scipy.signal import detrend, find_peaks, windows\n", ns)
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(SIM1D / "1D_batch" / "analyze_patterns.py"), "exec"), ns)
    return ns["analyze_pattern"]


# ------------------------------------------------------------------------------------------------ measurement
def domains(mask):
    """Contiguous runs of True: (count, largest, total, [(start, length), ...])."""
    m = np.asarray(mask, bool)
    if not m.any():
        return 0, 0, 0, []
    d = np.diff(np.concatenate([[0], m.astype(int), [0]]))
    starts = np.nonzero(d == 1)[0]; ends = np.nonzero(d == -1)[0]
    lengths = ends - starts
    return len(starts), int(lengths.max()), int(m.sum()), list(zip(starts.tolist(), lengths.tolist()))


def label(n_dom, largest, N):
    if n_dom == 0:
        return "all off"
    if n_dom == 1 and largest >= 0.95 * N:
        return "all on"
    return "one domain" if n_dom == 1 else "multiple domains"


def measure(a, i, ba, N, step, analyze_pattern):
    thr = 0.3 * ba
    n_dom, largest, total, runs = domains(a > thr)
    amax = float(a.max())
    half = a >= 0.5 * amax if amax > 0 else np.zeros_like(a, bool)
    k = int(np.argmax(a)); fw = 0
    if amax > 0:
        l = k
        while l > 0 and half[l - 1]:
            l -= 1
        r = k
        while r < N - 1 and half[r + 1]:
            r += 1
        fw = r - l + 1
    with contextlib.redirect_stdout(io.StringIO()):
        dom_freq, dom_lambda, score, amp_sd, flag = analyze_pattern(a, dx=1.0, plot=False)
    return dict(step_returned=int(step), steps_used=int(step) + 1, hit_cap=int(step >= STEPS - 1), label=label(n_dom, largest, N), n_domains=n_dom,
                largest_domain=largest, cells_above=total, fwhm=fw, domain_start=runs[0][0] if runs else -1,
                a_min=float(a.min()), a_max=amax, a_mean=float(a.mean()), i_mean=float(i.mean()), i_max=float(i.max()),
                dom_lambda=float(dom_lambda) if np.isfinite(dom_lambda) else np.inf, score=float(score), amp_sd=float(amp_sd), pattern_flag=bool(flag))


# ------------------------------------------------------------------------------------------------ jobs
def _cell_job(args):
    na, ni, N, gammas, spike, tag = args
    analyze_pattern = load_analyze_pattern()
    out = RES / tag; out.mkdir(parents=True, exist_ok=True)
    f = out / f"cell_na{na:02d}_ni{ni:02d}.npz"
    if f.exists():
        return f
    BA, BI, G = np.meshgrid(BA_VALUES, BI_VALUES, gammas, indexing="ij")
    ba, bi, g = BA.ravel(), BI.ravel(), G.ravel()
    t0 = time.perf_counter()
    a0, i0 = spike_init(len(ba), N, spike)
    A, I, S = integrate_batch(ba, bi, g, na, ni, N, a0, i0)
    t_int = time.perf_counter() - t0
    rows = [dict(b_a=float(ba[k]), b_i=float(bi[k]), gamma=float(g[k]), n_a=na, n_i=ni, N=N, spike=spike, **measure(A[k], I[k], ba[k], N, S[k], analyze_pattern))
            for k in range(len(ba))]
    with open(out / f"cell_na{na:02d}_ni{ni:02d}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    np.savez_compressed(f, b_a=ba, b_i=bi, gamma=g, n_a=na, n_i=ni, N=N, spike=spike, step=S, a_final=A.astype(np.float32), i_final=I.astype(np.float32))
    counts = {lab: sum(r["label"] == lab for r in rows) for lab in LABELS}
    print(f"cell n_a {na:2d} n_i {ni:2d}: {len(ba)} sets, integrate {t_int:.0f} s, total {time.perf_counter() - t0:.0f} s, capped {sum(r['hit_cap'] for r in rows)}, {counts}", flush=True)
    return f


def run_grid(workers=4, gammas=GAMMA_19, tag="full", N=N_DEFAULT, spike=SPIKE):
    jobs = [(int(na), int(ni), N, gammas, spike, tag) for na in HILL for ni in HILL]
    t0 = time.perf_counter()
    with Pool(workers) as pool:
        for k, f in enumerate(pool.imap_unordered(_cell_job, jobs, chunksize=1), 1):
            print(f"[{k}/{len(jobs)}] {time.perf_counter() - t0:.0f} s elapsed", flush=True)


def aggregate(tag="full", N=N_DEFAULT):
    import pandas as pd
    parts = sorted((RES / tag).glob("cell_na*_ni*.csv"))
    df = pd.concat([pd.read_csv(f) for f in parts], ignore_index=True)
    aud = load_audit(N)
    key = ["b_a", "b_i", "gamma", "n_a", "n_i"]
    for c in key:
        df[c] = df[c].round(10); aud[c] = aud[c].round(10)
    df = df.merge(aud[key + ["truth_has_state", "truth_a_star", "truth_i_star", "turing", "wavelength_fastest"]], on=key, how="left")
    df["lambda_i"] = np.sqrt(D_I / df.gamma)
    df["largest_domain_over_lambda_i"] = df.largest_domain / df.lambda_i
    df["fwhm_over_lambda_i"] = df.fwhm / df.lambda_i
    df = df.sort_values(["n_a", "n_i", "b_a", "b_i", "gamma"]).reset_index(drop=True)
    df.to_csv(SUMM / "spike_grid_full.csv", index=False, lineterminator="\n")
    rows = []
    for (na, ni), g in df.groupby(["n_a", "n_i"]):
        one = g[g.label == "one domain"]
        row = dict(n_a=na, n_i=ni, n_sets=len(g), **{f"n_{lab.replace(' ', '_')}": int((g.label == lab).sum()) for lab in LABELS},
                   n_hit_cap=int(g.hit_cap.sum()), n_truth_has_state=int(g.truth_has_state.sum()), n_turing=int(g.turing.sum()),
                   n_pattern_flag=int(g.pattern_flag.sum()),
                   one_domain_size_median=float(one.largest_domain.median()) if len(one) else np.nan,
                   one_domain_size_q25=float(one.largest_domain.quantile(0.25)) if len(one) else np.nan,
                   one_domain_size_q75=float(one.largest_domain.quantile(0.75)) if len(one) else np.nan,
                   one_domain_size_max=int(one.largest_domain.max()) if len(one) else np.nan,
                   one_domain_over_150=int((one.largest_domain > 150).sum()), one_domain_le5=int((one.largest_domain <= 5).sum()),
                   one_domain_size_over_lambda_median=float(one.largest_domain_over_lambda_i.median()) if len(one) else np.nan,
                   one_domain_size_over_lambda_q25=float(one.largest_domain_over_lambda_i.quantile(0.25)) if len(one) else np.nan,
                   one_domain_size_over_lambda_q75=float(one.largest_domain_over_lambda_i.quantile(0.75)) if len(one) else np.nan,
                   n_one_domain_and_turing=int((one.turing).sum()), n_one_domain_and_has_state=int(one.truth_has_state.sum()),
                   n_multiple_and_turing=int(((g.label == "multiple domains") & g.turing).sum()),
                   n_all_on_and_has_state=int(((g.label == "all on") & g.truth_has_state).sum()),
                   n_all_on_and_no_state=int(((g.label == "all on") & ~g.truth_has_state.astype(bool)).sum()))
        rows.append(row)
    by = pd.DataFrame(rows)
    by.to_csv(SUMM / "spike_grid_by_hill.csv", index=False, lineterminator="\n")
    print(f"spike_grid_full.csv: {len(df)} rows; spike_grid_by_hill.csv: {len(by)} rows")
    print(df.label.value_counts().to_string())
    return df, by


# ------------------------------------------------------------------------------------------------ checks and timing
def check():
    from simulation import run_simulation
    analyze_pattern = load_analyze_pattern()
    rng = np.random.default_rng(0)
    print("bit-exact check of the replica against run_simulation (final activator/inhibitor of the stale frame and step)")
    bad = 0
    for N in (100, 200):
        for k in range(15):
            na, ni = int(rng.integers(1, 13)), int(rng.integers(1, 13))
            ba, bi, g = float(rng.choice(BA_VALUES)), float(rng.choice(BI_VALUES)), float(rng.choice(GAMMA_19))
            if k < 4:
                na, ni = 3, 3           # Turing-prone cells that may run to the cap
            p = base_params(ba, bi, g, na, ni); p["N"] = N
            with contextlib.redirect_stdout(io.StringIO()):
                r = run_simulation(p)
            a, i, step = run_from_centre(p)
            # run_simulation's final frame is the save at the last step % save_every == 0; replicate it for the comparison
            if step == STEPS - 1:
                # capped: the engine's frame is from step 49800; recompute our state at that step
                a_frame, i_frame, _ = integrate_batch([ba], [bi], [g], na, ni, N, *spike_init(1, N), steps=STEPS - (STEPS - 1) % SAVE_EVERY, min_steps=10 ** 9)
                a_cmp, i_cmp = a_frame[0], i_frame[0]
            else:
                a_cmp, i_cmp = a, i
            same = np.array_equal(a_cmp, r["activator_final"]) and np.array_equal(i_cmp, r["inhibitor_final"]) and step == r["steps_used"]
            bad += not same
            print(f"  N {N} n {na:2d}/{ni:2d} b_a {ba:4.1f} b_i {bi:4.1f} gamma {g:4.1f}: step {step} vs {r['steps_used']}, identical {same}")
    print(f"mismatches: {bad} of 30")
    # cross-check: n_a 10, n_i 4 at N = 100
    N = 100
    BA, BI, G = np.meshgrid(BA_VALUES, BI_VALUES, GAMMA_19, indexing="ij"); ba, bi, g = BA.ravel(), BI.ravel(), G.ravel()
    t0 = time.perf_counter()
    A, I, S = integrate_batch(ba, bi, g, 10, 4, N, *spike_init(len(ba), N))
    labs = [measure(A[k], I[k], ba[k], N, S[k], analyze_pattern)["label"] for k in range(len(ba))]
    counts = {lab: labs.count(lab) for lab in LABELS}
    print(f"cross-check n_a 10 n_i 4 N 100 ({time.perf_counter() - t0:.0f} s): {counts}  (expected all off 390, all on 1293, one domain 217, multiple 0)")
    return bad, counts


def timing():
    analyze_pattern = load_analyze_pattern()
    rng = np.random.default_rng(1)
    N = N_DEFAULT; times = []; capped = 0
    for k in range(30):
        na, ni = (3, int(rng.integers(1, 13))) if k < 8 else (int(rng.integers(1, 13)), int(rng.integers(1, 13)))
        ba, bi, g = float(rng.choice(BA_VALUES)), float(rng.choice(BI_VALUES)), float(rng.choice(GAMMA_19))
        p = base_params(ba, bi, g, na, ni)
        t0 = time.perf_counter(); a, i, step = run_from_centre(p); dtm = time.perf_counter() - t0
        times.append(dtm); capped += step >= STEPS - 1
        print(f"  n {na:2d}/{ni:2d} b_a {ba:4.1f} b_i {bi:4.1f} gamma {g:4.1f}: {dtm:6.2f} s, step {step}", flush=True)
    times = np.array(times)
    print(f"30 single runs at N = 200: mean {times.mean():.2f} s, max {times.max():.2f} s, {capped} capped; naive projection 273,600 x mean / 4 cores = {273600 * times.mean() / 4 / 3600:.1f} h")
    t0 = time.perf_counter()
    f = _cell_job((3, 3, N, GAMMA_19, SPIKE, "timing"))
    t_cell = time.perf_counter() - t0
    print(f"batched cell n_a 3 n_i 3 (1900 sets): {t_cell:.0f} s -> projection 144 cells / 4 workers = {144 * t_cell / 4 / 3600:.2f} h (upper bound: this cell is Turing-heavy)")
    return times, t_cell


def validate(n_target=2000, seed=0):
    """Stratified random_tight validation of the audit: per (n_a, n_i) cell, up to 7 turing sets and 7 has_state & ~turing sets."""
    import pandas as pd
    aud = load_audit(N_DEFAULT)
    rng = np.random.default_rng(seed)
    per = max(1, int(round(n_target / (144 * 2))) + 2)      # per-cell cap per stratum; many cells hold fewer Turing sets than the cap
    picks = []
    for (na, ni), g in aud.groupby(["n_a", "n_i"]):
        t = g[g.turing]; s = g[g.truth_has_state & ~g.turing]
        for sub in (t, s):
            if len(sub):
                picks.append(sub.iloc[rng.choice(len(sub), size=min(per, len(sub)), replace=False)])
    sample = pd.concat(picks, ignore_index=True)
    print(f"validation sample: {len(sample)} sets ({int(sample.turing.sum())} turing, {int((~sample.turing).sum())} stable non-Turing)")
    rows = []
    for (na, ni), g in sample.groupby(["n_a", "n_i"]):
        for frac in (0.05,):
            a0, i0 = tight_init(g.truth_a_star.values, g.truth_i_star.values, N_DEFAULT, frac=frac)
            A, I, S = integrate_batch(g.b_a.values, g.b_i.values, g.gamma.values, int(na), int(ni), N_DEFAULT, a0, i0)
            for k, (_, r) in enumerate(g.iterrows()):
                rows.append(dict(b_a=r.b_a, b_i=r.b_i, gamma=r.gamma, n_a=int(na), n_i=int(ni), truth_a_star=r.truth_a_star, truth_i_star=r.truth_i_star, truth_has_state=bool(r.truth_has_state),
                                 turing=bool(r.turing), sigma_max=float(r.sigma_max) if "sigma_max" in r else np.nan, frac=frac, step=int(S[k]), hit_cap=int(S[k] >= STEPS - 1),
                                 a_min=float(A[k].min()), a_max=float(A[k].max()),
                                 a_mean=float(A[k].mean()), rel_amp=float((A[k].max() - A[k].min()) / r.truth_a_star), rel_mean=float(A[k].mean() / r.truth_a_star)))
    df = pd.DataFrame(rows)
    df["outcome"] = outcome(df.rel_amp, df.rel_mean)
    df["agrees"] = np.where(df.turing, df.outcome == "patterned", df.outcome == "uniform")
    # retest every disagreement with the early-stopping rule disabled (the full 50,000 steps, t = 500), at +-5 % and at +-0.5 %:
    # the stopping rule fires while a slowly growing mode is still small, which the first pass cannot distinguish from stability
    dis = df[~df.agrees]
    rows2 = []
    for (na, ni), g in dis.groupby(["n_a", "n_i"]):
        extra = {}
        for frac, tag in ((0.05, "5"), (0.005, "05")):
            a0, i0 = tight_init(g.truth_a_star.values, g.truth_i_star.values, N_DEFAULT, frac=frac)
            A, I, S = integrate_batch(g.b_a.values, g.b_i.values, g.gamma.values, int(na), int(ni), N_DEFAULT, a0, i0, min_steps=STEPS)
            extra[tag] = (A, S)
        for k, (idx, r) in enumerate(g.iterrows()):
            A5, S5 = extra["5"]; A05, S05 = extra["05"]
            rows2.append(dict(index=idx, rel_amp_full5=float((A5[k].max() - A5[k].min()) / r.truth_a_star), rel_mean_full5=float(A5[k].mean() / r.truth_a_star),
                              rel_amp_full05=float((A05[k].max() - A05[k].min()) / r.truth_a_star), rel_mean_full05=float(A05[k].mean() / r.truth_a_star)))
    cols = ["rel_amp_full5", "rel_mean_full5", "rel_amp_full05", "rel_mean_full05"]
    if rows2:
        df = df.join(pd.DataFrame(rows2).set_index("index"))
    else:
        for c in cols:
            df[c] = np.nan
    df["outcome_full5"] = np.where(df[cols[0]].notna(), outcome(df.rel_amp_full5, df.rel_mean_full5), "")
    df["outcome_full05"] = np.where(df[cols[2]].notna(), outcome(df.rel_amp_full05, df.rel_mean_full05), "")
    # classification of the disagreements after the retests
    def verdict(r):
        if r.agrees:
            return "agrees"
        if r.turing:
            if r.outcome_full05 == "patterned" or r.outcome_full5 == "patterned":
                return "agrees after full-length run (slow mode; stopping rule fired early)"
            if r.outcome_full05 == "off" or r.outcome_full5 == "off":
                return "Turing-unstable, fell to all-off (not patterned)"
            return "Turing predicted, stayed uniform through t = 500 (expected amplification exp(sigma t) = %.2g)" % np.exp(r.sigma_max * 500) if np.isfinite(r.sigma_max) else "Turing predicted, stayed uniform"
        if r.outcome_full05 == "uniform":
            return "finite-amplitude collapse (+-5 % falls, +-0.5 % stays uniform)"
        return "stable predicted, did not stay uniform at +-0.5 %"
    df["verdict"] = df.apply(verdict, axis=1)
    df["final_agrees"] = df.verdict.str.startswith("agrees")
    df.to_csv(SUMM / "spike_grid_validation.csv", index=False, lineterminator="\n")
    by = df.groupby(["n_a", "n_i"]).agg(n=("agrees", "size"), n_turing=("turing", "sum"), n_stable_nonturing=("turing", lambda x: int((~x).sum())),
                                        disagree_first_pass=("agrees", lambda x: int((~x).sum())),
                                        disagree_after_retest=("final_agrees", lambda x: int((~x).sum())),
                                        turing_fell_to_off=("verdict", lambda x: int(x.str.startswith("Turing-unstable, fell").sum())),
                                        finite_collapses=("verdict", lambda x: int(x.str.startswith("finite-amplitude").sum()))).reset_index()
    by["disagreement_rate"] = by.disagree_after_retest / by.n
    by.to_csv(SUMM / "spike_grid_validation_by_hill.csv", index=False, lineterminator="\n")
    print(f"validation: {len(df)} sets ({int(df.turing.sum())} Turing, {int((~df.turing).sum())} stable non-Turing); disagreements at the first pass: {int((~df.agrees).sum())}; "
          f"after the full-length retests: {int((~df.final_agrees).sum())}")
    print(df.verdict.value_counts().to_string())
    return df


def outcome(rel_amp, rel_mean):
    """patterned: activator range > 5 % of a*; off: range small and mean below a* / 2; uniform otherwise."""
    rel_amp = np.asarray(rel_amp, float); rel_mean = np.asarray(rel_mean, float)
    return np.where(rel_amp > 0.05, "patterned", np.where(rel_mean < 0.5, "off", "uniform"))


def spike2():
    """Spike 2.0 over the plane at a reduced b_a / b_i / gamma subsample (5 x 5 x 9 = 225 sets per cell)."""
    import pandas as pd
    ba_v = np.array([1.0, 3.0, 5.0, 7.0, 10.0]); bi_v = ba_v.copy(); g_v = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 5.0, 10.0])
    analyze_pattern = load_analyze_pattern()
    rows = []
    for na in HILL:
        for ni in HILL:
            BA, BI, G = np.meshgrid(ba_v, bi_v, g_v, indexing="ij"); ba, bi, g = BA.ravel(), BI.ravel(), G.ravel()
            A, I, S = integrate_batch(ba, bi, g, int(na), int(ni), N_DEFAULT, *spike_init(len(ba), N_DEFAULT, 2.0))
            for k in range(len(ba)):
                rows.append(dict(b_a=float(ba[k]), b_i=float(bi[k]), gamma=float(g[k]), n_a=int(na), n_i=int(ni), N=N_DEFAULT, spike=2.0,
                                 **measure(A[k], I[k], ba[k], N_DEFAULT, S[k], analyze_pattern)))
        print(f"spike 2: n_a {na} done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(SUMM / "spike_height_2.csv", index=False, lineterminator="\n")
    by = df.groupby(["n_a", "n_i"]).label.value_counts().unstack(fill_value=0)
    print(by.to_string())
    return df


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "audit":
        audit(int(sys.argv[2]) if len(sys.argv) > 2 else N_DEFAULT)
    elif cmd == "check":
        check()
    elif cmd == "time":
        timing()
    elif cmd == "run":
        args = sys.argv[2:]; g10 = "--gamma10" in args; args = [a for a in args if a != "--gamma10"]
        run_grid(int(args[0]) if args else 4, gammas=GAMMA_10 if g10 else GAMMA_19)
    elif cmd == "validate":
        validate()
    elif cmd == "spike2":
        spike2()
    elif cmd == "aggregate":
        aggregate()
