"""
Seeded, protocol-explicit wrapper around the unmodified engine (CLAUDE.md §14: wrap, don't edit).

- Seeds the numpy global RNG immediately before run_coupled_hex, which makes init, nucleation and
  noise reproducible in a single process (CLAUDE.md §7). Run replicates in separate processes.
- Chooses steps = n_frames*save_every + 1 so the LAST saved frame is the true final state
  (CLAUDE.md §8), and records the Euler step of every frame.
- Passes nucleation_rate explicitly, including 0.0 (CLAUDE.md §9).
"""
from __future__ import annotations

import contextlib, io, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "2D_simulations"))
from simulation_2D import run_coupled_hex  # noqa: E402

FIXED = dict(act_half_sat=1.0, inh_half_sat=1.0, act_decay_rate=1.0, basal_prod=0.0, act_diffusion=0.0)


def circuit(ba, bi, gamma, na, ni, D):
    return dict(FIXED, inh_diffusion=float(D), act_prod_rate=float(ba), inh_prod_rate=float(bi),
                inh_decay_rate=float(gamma), act_hill_coeff=na, inh_hill_coeff=ni)


def run(p: dict, *, protocol: str, seed: int, t_end: float, ny: int = 100, nx: int = 100, dt: float = 0.01,
        save_every: int = 200, spike_value: float | None = None, nucleation_rate: float | None = None,
        noise_amplitude: float = 0.0):
    """
    protocol 'synchronous': random_uniform_over0, spike_value 1.0, nucleation 0.0   (Supp Table 1, Fig 1C; NOTE the
                            deposit's committed batch template would have run this with nucleation 0.02, CLAUDE.md §9)
    protocol 'nucleation' : all_off, spike_value 2.0, nucleation_rate 0.01           (Supp Table 1, Fig 2L/3/S8D)
    Frames land at t = 0, 0.01, 10.01, 20.01, ... and the final frame at t_end + 0.01; the deposit's own "final" frames
    sit 199 steps earlier (CLAUDE.md §8). The nucleation kick is 2*a_ss, which equals 2*spike_value only when the
    steady-state solver finds no activated state (true for every committed panel at b_a = 5, CLAUDE.md §6).
    Returns dict with frames (list of (euler_step, a, i)), steps_used, a_ss, i_ss and the settings.
    """
    if protocol == "synchronous":
        init, sv, nr = "random_uniform_over0", 1.0, 0.0
    elif protocol == "nucleation":
        init, sv, nr = "all_off", 2.0, 0.01
    else:
        raise ValueError(protocol)
    if spike_value is not None:
        sv = spike_value
    if nucleation_rate is not None:
        nr = nucleation_rate
    total = int(round(t_end / dt))
    if total % save_every:
        raise ValueError(f"t_end/dt = {total} must be a multiple of save_every = {save_every} so the last frame is the final state")
    n_frames = total // save_every
    steps = n_frames * save_every + 1
    np.random.seed(seed)
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        A, I, step, a_ss, i_ss = run_coupled_hex(ny, nx, steps, dt, 1.0, p, stopping_threshold=0.0,
                                                 min_steps=steps + 1, init_mode=init, activator_type="juxtacrine",
                                                 spike_value=sv, save_every=save_every, nucleation_rate=nr,
                                                 noise_amplitude=noise_amplitude)
    # frame k>=1 is the state after (k-1)*save_every + 1 Euler steps; frame 0 is the initial condition
    euler_steps = [0] + [(k - 1) * save_every + 1 for k in range(1, len(A))]
    assert euler_steps[-1] == steps, (euler_steps[-1], steps)
    return dict(frames=list(zip(euler_steps, A, I)), steps_used=step + 1, a_ss=a_ss, i_ss=i_ss, dt=dt,
                settings=dict(protocol=protocol, seed=seed, init_mode=init, spike_value=sv, nucleation_rate=nr,
                              noise_amplitude=noise_amplitude, ny=ny, nx=nx, steps=steps, save_every=save_every, **p))
