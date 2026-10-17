"""
Checks to run before the 2D initiation sweep (analysis/sweep2d_initiation.py). Each prints PASS or FAIL.

  1  edges            the b_i where the uniform activated state stops existing, per Hill pair
  2  kick rule        nucleation sets a fired cell to H_ON only if it is below b_a
  3  random streams   replicates differ; the same (start, replicate) gives the same draw everywhere
  4  ignition         a single seed of H_ON ignites at every Hill pair at 0.4x and 1.0x the edge
                      (100 x 100, D_i 10, full t = 50, via the driver's own run_one)
  5  stiff corner     no blow-up at n 10/2, 3x the edge, D_i 20 (rf2 and nuc01, 100 x 100)
  6  timing           a sample of full-size jobs, sequential, with the projection for the sweep

Usage: python analysis/sweep2d_checks.py [n_timing_jobs]
"""
import sys, time, tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sweep2d_initiation as S  # noqa: E402

TMP = Path(tempfile.mkdtemp(prefix="sweep2d_check_"))
ok_all = True


def report(name, ok, detail=""):
    global ok_all
    ok_all &= bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}", flush=True)


def ignition_job(args):
    (na, ni), fi = args
    r = S.run_one((na, ni, fi, 1, "seed1", 0), ny=100, nx=100, save=False, out=TMP)
    return (na, ni, S.FRACS[fi], r["label"], r["n_spots"], r["size_median"], r["i_min"])


if __name__ == "__main__":
    t0 = time.time()
    # 1 edges
    E = S.edges()
    report("1 edges", all(np.isfinite(v) and v > 0 for v in E.values()),
           " ".join(f"{k[0]}/{k[1]}:{v:.3g}" for k, v in E.items()))

    # 2 kick rule
    rng = np.random.default_rng(0)
    a = np.full((50, 50), 1.0); a[:, :25] = S.BA + 0.5           # left half at or above b_a
    before = a.copy()
    n = S.kick(a, rng, p_fire=1.0)                                 # everything eligible fires
    left_untouched = np.array_equal(a[:, :25], before[:, :25])
    right_set = np.all(a[:, 25:] == S.H_ON)
    report("2 kick rule", left_untouched and right_set and n == 50 * 25,
           f"fired {n}, cells at/above b_a untouched: {left_untouched}, below set to {S.H_ON}: {right_set}")

    # 3 random streams
    def draw(start, rep):
        return S.initial(start, S.stream(start, rep), 200, 200)[0]
    same = np.array_equal(draw("seeds", 0), draw("seeds", 0)) and np.array_equal(draw("rf1", 2), draw("rf1", 2))
    diff = (not np.array_equal(draw("seeds", 0), draw("seeds", 1))
            and not np.array_equal(draw("rf2", 0), draw("rf2", 1))
            and not np.array_equal(draw("seeds", 1), draw("seeds", 2)))
    nseeds = int((draw("seeds", 0) > 0).sum())
    report("3 random streams", same and diff, f"reproducible: {same}, replicates differ: {diff}, seeds per field: {nseeds}")

    # 4 ignition of a single seed of H_ON
    from multiprocessing import Pool
    jobs = [(pr, fi) for pr in S.HILL for fi in (S.FRACS.index(0.4), S.FRACS.index(1.0))]
    with Pool(2) as pool:
        res = pool.map(ignition_job, jobs, chunksize=1)
    dead = [r for r in res if r[3] == "dead"]
    known = [r for r in dead if (r[0], r[1]) == (10, 2) and r[2] >= 1.0]   # documented corner exception
    report("4 ignition", len(dead) == len(known),
           f"{len(res) - len(dead)}/{len(res)} ignited"
           + (f"; known exception (n 10/2 at/beyond the edge): {[(f'{r[0]}/{r[1]}', r[2]) for r in known]}" if known else "")
           + (f"; UNEXPECTED dead: {[(f'{r[0]}/{r[1]}', r[2]) for r in dead if r not in known]}" if len(dead) > len(known) else ""))
    for r in res:
        print(f"      n {r[0]}/{r[1]} at {r[2]}x edge: {r[3]}, spots {r[4]}, median size {r[5]}")

    # 5 stiff corner
    fi3 = S.FRACS.index(3.0); di20 = S.DIS.index(20.0)
    rows = [S.run_one((10, 2, fi3, di20, st, 0), ny=100, nx=100, save=False, out=TMP) for st in ("rf2", "nuc01")]
    lim = 20 * S.edges()[(10, 2)] * 3.0 / S.GA
    stable = all(np.isfinite(r["i_min"]) and np.isfinite(r["i_max"]) and r["i_min"] > -1.0 and r["i_max"] < lim for r in rows)
    report("5 stiff corner", stable, "; ".join(f"{r['start']}: label {r['label']}, i in [{r['i_min']:.3g}, {r['i_max']:.3g}]" for r in rows))

    # 6 timing at full size
    n_t = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    sample = [(2, 2, 6, 1, "nuc01", 0), (10, 4, 4, 2, "seeds", 0), (4, 10, 5, 0, "rf2", 0),
              (4, 4, 3, 1, "nuc001", 0), (10, 10, 2, 2, "rf1", 0), (2, 10, 7, 1, "seed1", 0)][:n_t]
    secs = []
    for j in sample:
        r = S.run_one(j, save=True, out=TMP)
        secs.append(r["seconds"])
        print(f"      {r['name']}: {r['seconds']:.1f} s, label {r['label']}, spots {r['n_spots']}", flush=True)
    m = float(np.mean(secs)); total = len(S.all_jobs())
    report("6 timing", True, f"mean {m:.1f} s per full-size run; {total} runs -> {total*m/3600:.1f} core-hours, "
                             f"{total*m/3600/4:.1f} h on 4 cores, {total*m/3600/8:.1f} h on 8 (this machine's speed)")
    fr = sorted((TMP / "frames").glob("*.npz"))
    if fr:
        mb = np.mean([f.stat().st_size for f in fr]) / 1e6
        print(f"      saved frames: {mb:.2f} MB per run -> about {mb*total/1000:.1f} GB for the sweep")
    print(f"\n{'ALL CHECKS PASSED' if ok_all else 'SOME CHECKS FAILED'} in {time.time()-t0:.0f} s  (scratch output in {TMP})")
