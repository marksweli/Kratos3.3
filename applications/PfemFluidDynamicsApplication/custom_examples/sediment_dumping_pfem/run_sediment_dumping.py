"""
run_sediment_dumping.py
=======================
Main script for the sediment-dumping (水平抛泥) PFEM two-fluid simulation.

Physical scenario
-----------------
A rectangular sediment block (3 cm × 2 cm, concentration α_s = 0.606)
is placed at the free surface of a water tank (0.7 m × 0.7 m) with its
upper face aligned with the water surface.  The block is released from
rest (sediment has an initial settling velocity of 15.4 cm/s downwards)
and falls through the water under gravity and inter-phase drag.

Method: PFEM (Particle Finite Element Method)
  • Lagrangian particle tracking
  • Delaunay triangulation + Alpha-Shape mesh reconstruction each step
  • Fractional-step FEM pressure–velocity coupling  (P1 triangular elements)
  • Schiller–Naumann (1935) drag model
  • Two-fluid continuity + momentum equations (Eqs. 4.50 – 4.52)

Run
---
    python run_sediment_dumping.py

Dependencies
------------
    pip install -r requirements.txt

Outputs
-------
    output/snapshot_t*.png   – concentration, velocity, pressure fields
    output/diagnostics.png   – time-evolution of max α_s and |u_s|
"""

import os
import sys
import time
import numpy as np

# Ensure the current directory is on the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mesh_generator import create_particle_grid
from pfem_two_fluid_solver import PFEMTwoFluidSolver
from visualizer import Visualizer

# ============================================================
# Simulation parameters
# ============================================================

PARAMS = {
    # --- Physical properties ---
    'rho_f':  1000.0,        # Water density             [kg/m³]
    'rho_s':  2650.0,        # Sediment density          [kg/m³]
    'nu_f':   1.0e-6,        # Water kin. viscosity      [m²/s]
    'nu_s':   1.0e-5,        # Sediment kin. viscosity   [m²/s]  (estimate)
    'd_p':    0.8e-3,        # Sediment particle diam.   [m]
    'g':      9.81,          # Gravitational accel.      [m/s²]

    # --- Domain ---
    'Lx':          0.70,     # Domain width              [m]
    'Ly':          0.70,     # Domain height             [m]
    'water_level': 0.50,     # Initial free-surface height [m]

    # --- Sediment block (initial condition) ---
    'sed_x_left':   0.335,   # Block left edge           [m]
    'sed_x_right':  0.365,   # Block right edge  (3 cm)  [m]
    'sed_y_top':    0.500,   # Block top = free surface  [m]
    'sed_y_bottom': 0.480,   # Block bottom  (2 cm)      [m]
    'alpha_s_0':    0.606,   # Sediment concentration in block
    'alpha_s_bg':   0.001,   # Background concentration (pure water)
    # The 15.4 cm/s is the *terminal* settling velocity (reference value).
    # The simulation starts from rest; sediment accelerates under gravity.
    'w_terminal':   0.154,   # Terminal settling speed [m/s]  (reference only)

    # --- Numerical ---
    'dx':          0.010,    # Particle spacing          [m]  (1 cm)
    'dt':          0.001,    # Time step                 [s]
    't_end':       0.200,    # End time                  [s]
    'alpha_shape': 80.0,     # Alpha-shape parameter     [m⁻¹]

    # --- Output ---
    'output_dir':   './output',
    'output_times': [0.02, 0.05, 0.08, 0.11, 0.14, 0.17, 0.20],
}


# ============================================================
# Initialisation
# ============================================================

def initialise(params):
    """
    Build the initial particle distribution and set field values.

    Returns
    -------
    solver : PFEMTwoFluidSolver
    """
    dx = params['dx']

    # ── Regular grid of particles in the water region ──────────────
    positions = create_particle_grid(
        x_range=(0.0, params['Lx']),
        y_range=(0.0, params['water_level']),
        dx=dx,
    )
    N = len(positions)
    print(f"[Init] {N} particles  (spacing = {dx*100:.1f} cm)")

    # ── Volume fractions ────────────────────────────────────────────
    alpha_s = np.full(N, params['alpha_s_bg'])

    in_block = (
        (positions[:, 0] >= params['sed_x_left'])  &
        (positions[:, 0] <= params['sed_x_right']) &
        (positions[:, 1] >= params['sed_y_bottom']) &
        (positions[:, 1] <= params['sed_y_top'])
    )
    alpha_s[in_block] = params['alpha_s_0']

    n_sed = int(np.sum(in_block))
    print(f"[Init] Sediment block: {n_sed} particles  "
          f"(x=[{params['sed_x_left']},{params['sed_x_right']}] m  "
          f"y=[{params['sed_y_bottom']},{params['sed_y_top']}] m)")

    # ── Initial velocities ──────────────────────────────────────────
    # Simulation starts from rest; sediment accelerates naturally under
    # gravity and inter-phase drag (the 15.4 cm/s terminal velocity is
    # a reference, not an initial condition).
    uf = np.zeros((N, 2))   # Water at rest
    us = np.zeros((N, 2))   # Sediment at rest

    # ── Create and configure solver ─────────────────────────────────
    solver = PFEMTwoFluidSolver(params)
    solver.setup_particles(positions, alpha_s, uf, us)
    return solver


# ============================================================
# Main simulation loop
# ============================================================

def run():
    """Entry point – runs the full simulation and saves outputs."""

    print("=" * 65)
    print("  Sediment Dumping PFEM Simulation  /  水平抛泥双流体数值模拟")
    print("=" * 65)

    os.makedirs(PARAMS['output_dir'], exist_ok=True)

    # ── Initialise ──────────────────────────────────────────────────
    solver = initialise(PARAMS)
    viz    = Visualizer(PARAMS)

    dt           = PARAMS['dt']
    t_end        = PARAMS['t_end']
    output_times = sorted(PARAMS['output_times'])
    out_index    = 0          # pointer into output_times

    # Diagnostic lists (for time-series plot)
    diag_times     = []
    diag_max_as    = []
    diag_max_vel_s = []

    # ── Save initial state ──────────────────────────────────────────
    print("\n[t=0.000] Saving initial snapshot …")
    viz.save_snapshot(solver, 0.0, PARAMS['output_dir'])

    # ── Time loop ───────────────────────────────────────────────────
    step        = 0
    wall_t0     = time.time()
    n_steps_tot = int(round(t_end / dt))

    print(f"\n[Loop]  {n_steps_tot} steps × dt={dt} s  →  t_end={t_end} s")
    print("-" * 65)

    while solver.time < t_end - 0.5 * dt:

        solver.step(dt)
        step += 1

        # ── Diagnostics ─────────────────────────────────────────────
        diag_times.append(solver.time)
        diag_max_as.append(float(np.max(solver.alpha_s)))
        diag_max_vel_s.append(float(np.max(np.linalg.norm(solver.us, axis=1))))

        # ── Progress printout (every 20 steps) ──────────────────────
        if step % 20 == 0:
            elapsed = time.time() - wall_t0
            pct     = solver.time / t_end
            eta     = elapsed / max(pct, 1e-9) * (1.0 - pct)
            n_tri   = len(solver._idx) if solver._idx is not None else 0
            print(f"  t={solver.time:.3f}s  step={step:4d}  "
                  f"max_αs={diag_max_as[-1]:.3f}  "
                  f"max|us|={diag_max_vel_s[-1]:.3f} m/s  "
                  f"n_tri={n_tri:5d}  "
                  f"elapsed={elapsed:.1f}s  ETA={eta:.0f}s")

        # ── Save snapshot at requested output times ──────────────────
        while (out_index < len(output_times)
               and solver.time >= output_times[out_index] - 0.5 * dt):
            t_out = output_times[out_index]
            print(f"\n  ── Snapshot at t ≈ {t_out:.2f} s ──")
            viz.save_snapshot(solver, solver.time, PARAMS['output_dir'])
            out_index += 1

    # ── Summary ─────────────────────────────────────────────────────
    total_wall = time.time() - wall_t0
    print("\n" + "=" * 65)
    print(f"  Simulation complete in {total_wall:.1f} s")
    print(f"  Output directory: {os.path.abspath(PARAMS['output_dir'])}")
    print("=" * 65)

    # ── Diagnostics plot ────────────────────────────────────────────
    viz.plot_statistics(diag_times, diag_max_as, diag_max_vel_s,
                        PARAMS['output_dir'])


# ============================================================
# Entry point
# ============================================================

if __name__ == '__main__':
    run()
