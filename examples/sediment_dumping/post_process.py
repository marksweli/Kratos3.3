"""
Post-processing script for the sediment dumping simulation
(water-sediment two-phase MPM).

Reads the binary VTK output files produced by Kratos MPMApplication and
generates:
  1. A 4-panel snapshot figure showing water + sediment MPs at key times.
  2. A trajectory & spreading plot (centroid descent + lateral spread).

Usage (run from the examples/sediment_dumping directory):
    python3 post_process.py
"""

import os
import glob
import struct

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ---------------------------------------------------------------------------
# Binary VTK parser (big-endian, format 4.0)
# ---------------------------------------------------------------------------

def _parse_vtk_binary(filepath):
    """Parse a Kratos-generated binary VTK unstructured grid file."""
    with open(filepath, "rb") as fh:
        raw = fh.read()

    pos = 0

    def _readline():
        nonlocal pos
        end = raw.find(b'\n', pos)
        if end == -1:
            end = len(raw)
        line = raw[pos:end].decode('ascii', errors='replace').strip()
        pos = end + 1
        return line

    def _skip_newline():
        nonlocal pos
        if pos < len(raw) and raw[pos:pos+1] == b'\n':
            pos += 1

    def _read_floats(n):
        nonlocal pos
        nb = n * 4
        if pos + nb > len(raw):
            raise ValueError(
                f"Unexpected end of file reading {n} floats at offset {pos} "
                f"in '{filepath}'.")
        vals = struct.unpack(f'>{n}f', raw[pos:pos+nb])
        pos += nb
        return vals

    def _read_ints(n):
        nonlocal pos
        nb = n * 4
        if pos + nb > len(raw):
            raise ValueError(
                f"Unexpected end of file reading {n} ints at offset {pos} "
                f"in '{filepath}'.")
        vals = struct.unpack(f'>{n}i', raw[pos:pos+nb])
        pos += nb
        return vals

    _readline()  # # vtk DataFile Version 4.0
    _readline()  # vtk output
    _readline()  # BINARY
    _readline()  # DATASET UNSTRUCTURED_GRID

    points_line = _readline()   # POINTS N float
    n_pts = int(points_line.split()[1])

    coords = _read_floats(n_pts * 3)
    _skip_newline()
    pts = np.array(coords, dtype=np.float32).reshape(n_pts, 3)[:, :2]

    cells_line = _readline()
    cells_parts = cells_line.split()
    n_cells, n_ints = int(cells_parts[1]), int(cells_parts[2])
    _read_ints(n_ints)
    _skip_newline()

    ct_line = _readline()
    _read_ints(int(ct_line.split()[1]))
    _skip_newline()

    _readline()   # CELL_DATA N

    field_line = _readline()
    n_fields = int(field_line.split()[2])

    result = {'points': pts}
    for _ in range(n_fields):
        fl = _readline()
        parts = fl.split()
        fname, n_comp, n_tuples = parts[0], int(parts[1]), int(parts[2])
        vals = _read_floats(n_comp * n_tuples)
        _skip_newline()
        if n_comp >= 2:
            result[fname] = np.array(vals, dtype=np.float32).reshape(n_tuples, n_comp)[:, :2]
        else:
            result[fname] = np.array(vals, dtype=np.float32)

    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _step_number(basename):
    """Return the integer step number from a filename like MPM_Material_0_42.vtk."""
    stem = os.path.splitext(basename)[0]
    parts = stem.rsplit('_', 1)
    if len(parts) != 2 or not parts[1].isdigit():
        raise ValueError(
            f"Cannot extract step number from VTK filename: '{basename}'. "
            "Expected format: 'MPM_Material_0_<step>.vtk'.")
    return int(parts[1])


def _collect_files(vtk_dir, dt=0.001):
    """Return a sorted list of (time, filepath) tuples for MPM_Material VTK files."""
    pattern = os.path.join(vtk_dir, "MPM_Material_0_*.vtk")
    files = sorted(glob.glob(pattern),
                   key=lambda p: _step_number(os.path.basename(p)))
    return [(_step_number(os.path.basename(p)) * dt, p) for p in files]


# Density threshold separating water (≤1200 kg/m³) from sediment (>1200 kg/m³).
# Chosen as the midpoint between pure water (1000) and sediment mixture (1800).
_SEDIMENT_DENSITY_THRESHOLD = 1200.0

# Height below which the sediment centroid is considered to have reached the bottom.
_IMPACT_HEIGHT_THRESHOLD = 0.05  # metres


def _split_phases(data):
    """Split MPs into water (density ≤ _SEDIMENT_DENSITY_THRESHOLD) and sediment."""
    pts = data['points']
    densities = data.get('MP_DENSITY', np.ones(len(pts), dtype=np.float32) * 1000.0)
    sed_mask = densities > _SEDIMENT_DENSITY_THRESHOLD
    return pts[sed_mask], pts[~sed_mask], sed_mask


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_snapshot(ax, data, time, domain=(0.7, 0.7)):
    """Render one time snapshot showing water + sediment MPs."""
    sed_pts, wat_pts, sed_mask = _split_phases(data)

    ax.set_facecolor('#d6eaf8')
    ax.set_xlim(0, domain[0])
    ax.set_ylim(0, domain[1])

    if len(wat_pts) > 0:
        ax.scatter(wat_pts[:, 0], wat_pts[:, 1],
                   c='#5b9bd5', s=18, alpha=0.5, edgecolors='none', zorder=2)

    if len(sed_pts) > 0:
        vx_abs = np.zeros(len(sed_pts))
        if 'MP_VELOCITY' in data:
            vels = data['MP_VELOCITY']
            vx_abs = np.abs(vels[sed_mask, 0])
        sc = ax.scatter(sed_pts[:, 0], sed_pts[:, 1],
                        c=vx_abs, cmap='YlOrRd', vmin=0, vmax=3,
                        s=40, edgecolors='k', linewidths=0.2, zorder=3)
        x_spread = float(np.max(sed_pts[:, 0]) - np.min(sed_pts[:, 0]))
        ax.text(0.02, 0.97, f'Δx = {x_spread:.3f} m',
                transform=ax.transAxes, fontsize=7, va='top',
                bbox=dict(facecolor='white', alpha=0.75, edgecolor='none', pad=2))

    ax.set_xlabel('x (m)', fontsize=8)
    ax.set_ylabel('y (m)', fontsize=8)
    ax.set_title(f't = {time:.3f} s', fontsize=9, fontweight='bold')
    ax.set_aspect('equal')
    ax.tick_params(labelsize=7)
    ax.grid(True, linestyle='--', linewidth=0.3, alpha=0.5)


def main():
    vtk_dir    = 'vtk_output'
    output_dir = 'results'
    os.makedirs(output_dir, exist_ok=True)

    time_files = _collect_files(vtk_dir, dt=0.001)
    if not time_files:
        print(f"No VTK files found in '{vtk_dir}'. Run the simulation first.")
        return

    all_times = [t for t, _ in time_files]
    print(f"Found {len(time_files)} snapshots: "
          f"t = {all_times[0]:.3f} … {all_times[-1]:.3f} s")

    # ------------------------------------------------------------------
    # 1) Four-panel snapshot: t = 0, 0.15, 0.35, 0.50 s
    # ------------------------------------------------------------------
    target_times = [0.0, 0.15, 0.35, 0.50]
    selected = []
    for tgt in target_times:
        closest_t, closest_fp = min(time_files, key=lambda tf: abs(tf[0] - tgt))
        selected.append((closest_t, closest_fp))

    fig, axes = plt.subplots(1, 4, figsize=(18, 5.5))
    for ax, (t, fp) in zip(axes, selected):
        _plot_snapshot(ax, _parse_vtk_binary(fp), t)

    sm = plt.cm.ScalarMappable(cmap='YlOrRd',
                                norm=plt.Normalize(vmin=0, vmax=3))
    sm.set_array([])
    fig.colorbar(sm, ax=axes.ravel().tolist(),
                 label='Sediment lateral velocity |vx| (m/s)',
                 shrink=0.55, pad=0.02)

    wat_p = mpatches.Patch(color='#5b9bd5', alpha=0.7, label='Water MPs')
    sed_p = mpatches.Patch(facecolor='#d35400', edgecolor='k',
                           linewidth=0.5, label='Sediment MPs')
    fig.legend(handles=[wat_p, sed_p], loc='lower left',
               bbox_to_anchor=(0.01, 0.01), fontsize=8)

    fig.suptitle(
        'Water–Sediment Two-Phase MPM Simulation: Sediment Cloud Settling and Spreading\n'
        r'$\rho_{\rm sed}=1800\ \mathrm{kg/m^3}$, '
        r'$\mu_{\rm sed}=0.05\ \mathrm{Pa\cdot s}$  |  '
        r'$\rho_{\rm wat}=1000\ \mathrm{kg/m^3}$, '
        r'$\mu_{\rm wat}=0.001\ \mathrm{Pa\cdot s}$  |  '
        'DispNewtonianFluidPlaneStrain2DLaw',
        fontsize=9, fontweight='bold', y=1.03)
    fig.tight_layout()

    snap_path = os.path.join(output_dir, 'sediment_snapshots.png')
    fig.savefig(snap_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Snapshot figure saved → '{snap_path}'")

    # ------------------------------------------------------------------
    # 2) Settling + spreading trajectories
    # ------------------------------------------------------------------
    t_list, cy_list, xspread_list = [], [], []
    for t, fp in time_files:
        d = _parse_vtk_binary(fp)
        sed_pts, _, _ = _split_phases(d)
        if len(sed_pts) == 0:
            continue
        cy_list.append(float(np.mean(sed_pts[:, 1])))
        xspread_list.append(float(np.max(sed_pts[:, 0]) - np.min(sed_pts[:, 0])))
        t_list.append(t)

    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(12, 4.5))

    ax2a.plot(t_list, cy_list, 'b-', linewidth=2, label='Sediment centroid (y)')
    ax2a.axhline(0, color='k', linestyle='--', linewidth=0.8, label='Bottom wall')
    ax2a.fill_between(t_list, 0, cy_list, alpha=0.15, color='blue')
    ax2a.set_xlabel('Time (s)', fontsize=11)
    ax2a.set_ylabel('Sediment centroid y-position (m)', fontsize=11)
    ax2a.set_title('Vertical settling of sediment cloud', fontsize=11)
    ax2a.legend(fontsize=9)
    ax2a.set_ylim(-0.03, 0.72)
    ax2a.grid(True, linestyle='--', linewidth=0.4)

    ax2b.plot(t_list, xspread_list, 'r-', linewidth=2,
              label='x-spread of sediment')
    impact_t = next((t for t, y in zip(t_list, cy_list)
                     if y < _IMPACT_HEIGHT_THRESHOLD), None)
    if impact_t:
        ax2b.axvline(x=impact_t, color='gray', linestyle=':', linewidth=1.2,
                     label=f'Impact (t≈{impact_t:.2f} s)')
    ax2b.fill_between(t_list, xspread_list[0], xspread_list, alpha=0.2, color='red')
    ax2b.set_xlabel('Time (s)', fontsize=11)
    ax2b.set_ylabel('Lateral spread Δx (m)', fontsize=11)
    ax2b.set_title('Lateral spreading of sediment cloud', fontsize=11)
    ax2b.legend(fontsize=9)
    ax2b.set_ylim(0, 0.72)
    ax2b.grid(True, linestyle='--', linewidth=0.4)

    fig2.suptitle(
        'Sediment Dumping: Gravitational Settling and Lateral Spreading\n'
        '(water-sediment two-phase MPM — DispNewtonianFluidPlaneStrain2DLaw)',
        fontsize=11, fontweight='bold')
    fig2.tight_layout()
    traj_path = os.path.join(output_dir, 'centroid_trajectory.png')
    fig2.savefig(traj_path, dpi=150)
    plt.close(fig2)
    print(f"Trajectory figure saved → '{traj_path}'")

    # ------------------------------------------------------------------
    # 3) Summary
    # ------------------------------------------------------------------
    print("\n--- Simulation Summary ---")
    print(f"Initial sediment centroid y  : {cy_list[0]:.4f} m")
    print(f"Final   sediment centroid y  : {cy_list[-1]:.4f} m")
    print(f"Initial x-spread             : {xspread_list[0]:.4f} m")
    print(f"Final   x-spread             : {xspread_list[-1]:.4f} m")
    print(f"Spreading ratio              : {xspread_list[-1]/xspread_list[0]:.2f}×")


if __name__ == '__main__':
    main()
