"""
visualizer.py
=============
Visualisation module for the PFEM sediment-dumping simulation.

Produces:
  - Scatter plot of particle positions coloured by sediment volume fraction α_s
  - Quiver plot of the water/sediment velocity fields
  - Combined figure saved as a PNG image

Usage
-----
    from visualizer import Visualizer
    viz = Visualizer(params)
    viz.save_snapshot(solver, t, output_dir)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')            # non-interactive backend for script use
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


class Visualizer:
    """
    Handles all visualisation tasks for the two-fluid PFEM simulation.

    Parameters
    ----------
    params : dict
        Simulation parameters (used for domain dimensions and labels).
    """

    def __init__(self, params):
        self.params = params
        self.Lx = params['Lx']
        self.Ly = params.get('water_level', params['Ly'])
        self.cmap_conc = 'YlOrRd'   # concentration colour map
        self.cmap_pres = 'RdBu_r'   # pressure colour map

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def save_snapshot(self, solver, t, output_dir):
        """
        Create and save a composite figure with three sub-panels:
          (a) Sediment concentration α_s
          (b) Sediment velocity magnitude
          (c) Pressure field

        Parameters
        ----------
        solver     : PFEMTwoFluidSolver  – current solver state
        t          : float               – current simulation time [s]
        output_dir : str                 – directory to write PNG
        """
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle(f'Sediment Dumping  –  PFEM Two-Fluid  |  t = {t:.3f} s',
                     fontsize=13, fontweight='bold')

        self._plot_concentration(axes[0], solver, t)
        self._plot_velocity(axes[1], solver, t)
        self._plot_pressure(axes[2], solver, t)

        plt.tight_layout()
        fname = os.path.join(output_dir, f'snapshot_t{t:.4f}.png')
        fig.savefig(fname, dpi=120, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved → {fname}")

    # ------------------------------------------------------------------
    # Sub-panel plotting helpers
    # ------------------------------------------------------------------

    def _scatter_base(self, ax, solver, field, cmap, vmin, vmax, label):
        """
        Common scatter-plot routine.

        Returns the PathCollection for colour-bar attachment.
        """
        sc = ax.scatter(
            solver.pos[:, 0], solver.pos[:, 1],
            c=field, cmap=cmap,
            vmin=vmin, vmax=vmax,
            s=6, linewidths=0, rasterized=True,
        )
        ax.set_xlim(0, self.Lx)
        ax.set_ylim(0, self.Ly * 1.05)
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]', fontsize=9)
        ax.set_ylabel('y [m]', fontsize=9)
        plt.colorbar(sc, ax=ax, label=label, fraction=0.046, pad=0.04)
        return sc

    def _plot_concentration(self, ax, solver, t):
        """
        Panel (a): sediment volume fraction α_s.
        """
        ax.set_title('Sediment concentration  α_s', fontsize=10)
        self._scatter_base(
            ax, solver,
            field=solver.alpha_s,
            cmap=self.cmap_conc,
            vmin=0.0, vmax=0.65,
            label='α_s',
        )
        # Annotate the initial water level
        wl = self.params.get('water_level', self.Ly)
        ax.axhline(wl, color='steelblue', lw=0.8, ls='--', alpha=0.5,
                   label='initial free surface')
        ax.legend(fontsize=7, loc='upper right')

    def _plot_velocity(self, ax, solver, t):
        """
        Panel (b): sediment velocity magnitude |u_s| with quiver arrows.
        """
        us_mag = np.linalg.norm(solver.us, axis=1)
        ax.set_title('Sediment velocity  |u_s|  [m/s]', fontsize=10)
        self._scatter_base(
            ax, solver,
            field=us_mag,
            cmap='viridis',
            vmin=0.0, vmax=max(0.3, np.max(us_mag) * 0.8),
            label='|u_s|  [m/s]',
        )

        # Quiver on a coarser sub-grid to avoid clutter
        skip = max(1, len(solver.pos) // 400)
        idx  = np.arange(0, len(solver.pos), skip)
        ax.quiver(
            solver.pos[idx, 0], solver.pos[idx, 1],
            solver.us[idx, 0],  solver.us[idx, 1],
            angles='xy', scale_units='xy',
            scale=3.0,   width=0.002,
            color='white', alpha=0.7,
        )

    def _plot_pressure(self, ax, solver, t):
        """
        Panel (c): dynamic pressure field.
        """
        p = solver.pressure
        p_range = max(np.max(np.abs(p)), 1.0)
        ax.set_title('Pressure  p  [Pa]', fontsize=10)
        self._scatter_base(
            ax, solver,
            field=p,
            cmap=self.cmap_pres,
            vmin=-p_range, vmax=p_range,
            label='p  [Pa]',
        )

    # ------------------------------------------------------------------
    # Convenience: plot summary statistics over time
    # ------------------------------------------------------------------

    def plot_statistics(self, times, max_alpha_s, max_vel_s, output_dir):
        """
        Save a figure with time-evolution of key diagnostics.

        Parameters
        ----------
        times       : list of float   – simulation times
        max_alpha_s : list of float   – max sediment concentration
        max_vel_s   : list of float   – max sediment speed
        output_dir  : str
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        fig.suptitle('Simulation Diagnostics', fontsize=12)

        ax1.plot(times, max_alpha_s, 'o-', color='darkorange', ms=3)
        ax1.set_ylabel('max α_s  [–]', fontsize=9)
        ax1.grid(True, alpha=0.3)

        ax2.plot(times, max_vel_s, 's-', color='steelblue', ms=3)
        ax2.set_ylabel('max |u_s|  [m/s]', fontsize=9)
        ax2.set_xlabel('time  [s]', fontsize=9)
        ax2.grid(True, alpha=0.3)

        for ax in (ax1, ax2):
            for t_ref in (0.08, 0.14):
                ax.axvline(t_ref, color='gray', lw=0.8, ls='--')

        plt.tight_layout()
        fname = os.path.join(output_dir, 'diagnostics.png')
        fig.savefig(fname, dpi=120, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved → {fname}")
