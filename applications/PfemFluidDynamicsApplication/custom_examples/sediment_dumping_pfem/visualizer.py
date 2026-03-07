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
        Create and save a composite figure with four sub-panels:
          (a) Full domain – sediment concentration α_s
          (b) Full domain – sediment velocity |u_s|
          (c) Zoomed view – α_s around the sediment cloud
          (d) Zoomed view – pressure field

        Parameters
        ----------
        solver     : PFEMTwoFluidSolver  – current solver state
        t          : float               – current simulation time [s]
        output_dir : str                 – directory to write PNG
        """
        fig = plt.figure(figsize=(16, 7))
        fig.suptitle(
            f'Sediment Dumping  –  PFEM Two-Fluid  |  t = {t:.3f} s',
            fontsize=13, fontweight='bold')

        ax1 = fig.add_subplot(1, 4, 1)
        ax2 = fig.add_subplot(1, 4, 2)
        ax3 = fig.add_subplot(1, 4, 3)
        ax4 = fig.add_subplot(1, 4, 4)

        self._plot_concentration(ax1, solver, t)
        self._plot_velocity(ax2, solver, t)
        self._plot_concentration_zoomed(ax3, solver, t)
        self._plot_pressure(ax4, solver, t)

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
        Sediment-rich particles are shown with larger markers.
        """
        ax.set_title('Sediment concentration  α_s', fontsize=10)
        thresh = 0.01
        bg     = solver.alpha_s < thresh
        fg     = ~bg

        # Background (water)
        ax.scatter(solver.pos[bg, 0], solver.pos[bg, 1],
                   c='#e8f4e8', s=3, linewidths=0, rasterized=True)
        # Sediment particles
        sc = ax.scatter(
            solver.pos[:, 0], solver.pos[:, 1],
            c=solver.alpha_s, cmap=self.cmap_conc,
            vmin=0.0, vmax=0.65,
            s=np.where(fg, 60, 4), linewidths=0, rasterized=True,
        )
        ax.set_xlim(0, self.Lx)
        ax.set_ylim(0, self.Ly * 1.05)
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]', fontsize=9)
        ax.set_ylabel('y [m]', fontsize=9)
        plt.colorbar(sc, ax=ax, label='α_s', fraction=0.046, pad=0.04)

        wl = self.params.get('water_level', self.Ly)
        ax.axhline(wl, color='steelblue', lw=0.8, ls='--', alpha=0.5,
                   label='initial free surface')
        ax.legend(fontsize=7, loc='upper right')

    def _plot_velocity(self, ax, solver, t):
        """
        Panel (b): sediment velocity magnitude |u_s| with quiver arrows.
        Only particles with significant sediment (α_s > threshold) are shown
        with large markers; the rest appear as small background dots.
        """
        us_mag  = np.linalg.norm(solver.us, axis=1)
        thresh  = 0.01   # α_s threshold to highlight sediment particles

        ax.set_title('Sediment velocity  |u_s|  [m/s]', fontsize=10)

        # Background particles (pure water)
        bg = solver.alpha_s < thresh
        ax.scatter(
            solver.pos[bg, 0], solver.pos[bg, 1],
            c='#d0e8f0', s=3, linewidths=0, rasterized=True,
        )

        # Sediment-rich particles
        fg = ~bg
        if np.any(fg):
            vmax_vel = max(0.3, np.max(us_mag[fg]) * 0.9)
            sc = ax.scatter(
                solver.pos[fg, 0], solver.pos[fg, 1],
                c=us_mag[fg], cmap='viridis',
                vmin=0.0, vmax=vmax_vel,
                s=60, linewidths=0, rasterized=True, zorder=5,
            )
            plt.colorbar(sc, ax=ax, label='|u_s|  [m/s]',
                         fraction=0.046, pad=0.04)
            # Quiver for sediment particles
            ax.quiver(
                solver.pos[fg, 0], solver.pos[fg, 1],
                solver.us[fg, 0],  solver.us[fg, 1],
                angles='xy', scale_units='xy',
                scale=4.0, width=0.004,
                color='yellow', alpha=0.85, zorder=6,
            )
        else:
            ax.scatter([], [], c=[], cmap='viridis', vmin=0, vmax=0.3, s=1)
            plt.colorbar(ax.collections[-1], ax=ax,
                         label='|u_s|  [m/s]', fraction=0.046, pad=0.04)

        ax.set_xlim(0, self.Lx)
        ax.set_ylim(0, self.Ly * 1.05)
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]', fontsize=9)
        ax.set_ylabel('y [m]', fontsize=9)

    def _plot_concentration_zoomed(self, ax, solver, t):
        """
        Panel (c): zoomed view of α_s around the sediment cloud centre.
        The zoom window tracks the centre-of-mass of the sediment cloud.
        """
        # Find centre of mass of the sediment cloud (weighted by α_s)
        as_ = solver.alpha_s
        total = np.sum(as_)
        if total > 1e-8:
            cx = np.dot(as_, solver.pos[:, 0]) / total
            cy = np.dot(as_, solver.pos[:, 1]) / total
        else:
            cx, cy = self.Lx / 2.0, self.Ly * 0.8

        # Zoom window  (±0.15 m around CoM, clamped to domain)
        hw = 0.18
        x0 = max(0.0,    cx - hw)
        x1 = min(self.Lx, cx + hw)
        y0 = max(0.0,    cy - hw)
        y1 = min(self.Ly * 1.1, cy + hw)

        ax.set_title('Zoomed: α_s (cloud CoM)', fontsize=10)
        sc = ax.scatter(
            solver.pos[:, 0], solver.pos[:, 1],
            c=solver.alpha_s, cmap=self.cmap_conc,
            vmin=0.0, vmax=0.65,
            s=30, linewidths=0, rasterized=True,
        )
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]', fontsize=9)
        ax.set_ylabel('y [m]', fontsize=9)
        plt.colorbar(sc, ax=ax, label='α_s', fraction=0.046, pad=0.04)

        # Centre-of-mass marker
        ax.plot(cx, cy, 'k+', ms=8, mew=1.5, label=f'CoM y={cy:.3f}m')
        ax.legend(fontsize=7, loc='upper right')

        # Initial free surface
        wl = self.params.get('water_level', self.Ly)
        ax.axhline(wl, color='steelblue', lw=0.8, ls='--', alpha=0.6)
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
