"""
Nucleation-rate sweep at the published Fig 2L parameter set.

Everything is copied verbatim from analysis/ladder_batch CONDS['nucl'] and the
day-4 driver; the ONLY thing that varies is nucleation_rate and the replicate
seed. Published value is 0.01.

Fixed:  b_a 5, b_i 25, gamma 0.5, n_a 10, n_i 4, D_i 20   (set 'bi25', 2L_JAPI)
        protocol all_off, spike_value 2.0, set_peak_height 10
        200 x 200 cells = 27.2 mm2 at 28 um spacing
        t_end 30 (day 6.8), frame every t = 0.5, so days 1-4 are all readable
"""
import contextlib, io, os, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "2D_simulations"))
import analysis.ladder_batch as L

OUT = Path(os.environ.get("NUCSWEEP_OUT", ROOT / "results" / "nucleation_sweep"))
OUT.mkdir(parents=True, exist_ok=True)
N, DT, T_END, SAVE_EVERY = 200, 0.01, 30.0, 50
SETNAME = "bi25"
RATES = [0.00125, 0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16]
REPS = [1, 2, 3]


def one(rate, rep):
    import simulation_2D
    from simulation_2D import run_coupled_hex
    bi, anchor, D = L.set_params(SETNAME)
    p = L.circuit(L.BA, bi, L.GAMMA, L.NA, L.NI, D)
    c = dict(L.CONDS["nucl"])
    steps = int(round(T_END / DT)) + 1
    np.random.seed(rep)
    with contextlib.redirect_stdout(io.StringIO()):
        A, I, step, a_ss, i_ss = run_coupled_hex(
            N, N, steps, DT, 1.0, p, 0.0, 10**9, init_mode=c["init_mode"],
            activator_type="juxtacrine", n_points=c["n_points"],
            spike_value=c["spike_value"], save_every=SAVE_EVERY,
            nucleation_rate=rate, noise_amplitude=0.0,
            set_peak_height=c["set_peak_height"])
    ts = np.array([0.0] + [((k - 1) * SAVE_EVERY + 1) * DT for k in range(1, len(A))])
    np.savez_compressed(OUT / f"nucl_r{rate:g}_rep{rep}.npz",
                        a=np.asarray(A, np.float16), i=np.asarray(I, np.float16),
                        t=ts, b_i=bi, D_i=D, n_a=L.NA, n_i=L.NI,
                        rep=rep, nucleation_rate=rate, anchor=anchor)
    return len(A)


if __name__ == "__main__":
    t_all = time.time()
    for rate in RATES:
        for rep in REPS:
            f = OUT / f"nucl_r{rate:g}_rep{rep}.npz"
            if f.exists():
                continue
            t0 = time.time()
            try:
                n = one(rate, rep)
                print(f"[{time.strftime('%H:%M:%S')}] rate={rate:<8g} rep={rep}  "
                      f"{n} frames  {time.time()-t0:5.1f}s", flush=True)
            except Exception as e:
                print(f"FAILED rate={rate} rep={rep}: {e!r}", flush=True)
    print(f"[sweep done in {(time.time()-t_all)/60:.1f} min]", flush=True)
