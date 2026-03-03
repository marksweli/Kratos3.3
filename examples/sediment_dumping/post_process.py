"""
Post-processing script for the sediment dumping simulation.

Reads the binary VTK output files produced by Kratos MPMApplication and
generates snapshot images of the sediment material-point positions at
specified times, together with a centroid-descent trajectory plot.

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
import matplotlib.cm as cm


# ---------------------------------------------------------------------------
# Binary VTK parser (big-endian, format 4.0)
# ---------------------------------------------------------------------------

def _parse_vtk_binary(filepath):
    """Parse a Kratos-generated binary VTK unstructured grid file.

    Returns a dict with keys:
      'points'          – (N, 2) float array of current MP positions
      'MP_VELOCITY'     – (N, 2) float array  [optional]
      'MP_DISPLACEMENT' – (N, 2) float array  [optional]
      'MP_DENSITY'      – (N,)   float array  [optional]
      'MP_VOLUME'       – (N,)   float array  [optional]
    """
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

    # --- header (ASCII) ---
    _readline()              # # vtk DataFile Version 4.0
    _readline()              # vtk output
    _readline()              # BINARY
    _readline()              # DATASET UNSTRUCTURED_GRID

    points_line = _readline()   # POINTS N float
    n_pts = int(points_line.split()[1])

    coords = _read_floats(n_pts * 3)
    _skip_newline()
    pts = np.array([(coords[i*3], coords[i*3+1]) for i in range(n_pts)],
                   dtype=np.float32)

    # CELLS N size
    cells_line = _readline()
    parts = cells_line.split()
    n_cells, n_ints = int(parts[1]), int(parts[2])
    _read_ints(n_ints)
    _skip_newline()

    # CELL_TYPES N
    ct_line = _readline()
    n_ct = int(ct_line.split()[1])
    _read_ints(n_ct)
    _skip_newline()

    # CELL_DATA N or POINT_DATA N
    _readline()

    # FIELD FieldData N
    field_line = _readline()
    n_fields = int(field_line.split()[2])

    result = {'points': pts}
    for _ in range(n_fields):
        fl = _readline()   # name n_comp n_tuples dtype
        parts = fl.split()
        fname, n_comp, n_tuples = parts[0], int(parts[1]), int(parts[2])
        vals = _read_floats(n_comp * n_tuples)
        _skip_newline()
        if n_comp >= 2:
            result[fname] = np.array(
                [(vals[i*n_comp], vals[i*n_comp+1]) for i in range(n_tuples)],
                dtype=np.float32)
        else:
            result[fname] = np.array(vals, dtype=np.float32)

    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _step_number(basename):
    """Return the integer step number from a filename like MPM_Material_0_42.vtk.

    Raises ValueError if the filename doesn't match the expected pattern.
    """
    stem = os.path.splitext(basename)[0]   # MPM_Material_0_42
    parts = stem.rsplit('_', 1)
    if len(parts) != 2 or not parts[1].isdigit():
        raise ValueError(
            f"Cannot extract step number from VTK filename: '{basename}'. "
            "Expected format: 'MPM_Material_0_<step>.vtk'.")
    return int(parts[1])


def _collect_files(vtk_dir, dt=0.001):
    """Return a sorted list of (time, filepath) tuples for MPM_Material VTK files."""
    pattern = os.path.join(vtk_dir, "MPM_Material_0_*.vtk")
    files = sorted(glob.glob(pattern), key=lambda p: _step_number(os.path.basename(p)))
    return [((_step_number(os.path.basename(p))) * dt, p) for p in files]


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _plot_snapshot(ax, data, time, vel_data=None, domain=(0.7, 0.7)):
    """Render one time snapshot on *ax*."""
    pts = data['points']          # current positions (Kratos writes updated coords)
    if pts.shape[0] == 0:
        ax.set_title(f"t = {time:.3f} s  (no data)")
        return

    # colour by vertical velocity magnitude if available
    if vel_data is not None and 'MP_VELOCITY' in vel_data:
        vy = np.abs(vel_data['MP_VELOCITY'][:, 1])
        c_vals = vy
        cmap = cm.YlOrRd
    else:
        c_vals = 'saddlebrown'
        cmap = None

    ax.set_facecolor('#cce8f4')   # water blue
    sc = ax.scatter(pts[:, 0], pts[:, 1],
                    c=c_vals, cmap=cmap,
                    s=55, edgecolors='k', linewidths=0.25,
                    zorder=3)
    if cmap is not None:
        plt.colorbar(sc, ax=ax, label='|vy| (m/s)', shrink=0.7, pad=0.02)

    ax.set_xlim(0, domain[0])
    ax.set_ylim(0, domain[1])
    ax.set_xlabel('x (m)', fontsize=9)
    ax.set_ylabel('y (m)', fontsize=9)
    ax.set_title(f't = {time:.3f} s', fontsize=10, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', linewidth=0.35, alpha=0.55)

    water_p = mpatches.Patch(color='#cce8f4', label='Water')
    sed_p   = mpatches.Patch(facecolor='saddlebrown', edgecolor='k',
                             linewidth=0.5, label='Sediment MPs')
    ax.legend(handles=[water_p, sed_p], fontsize=7, loc='lower right')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    vtk_dir    = 'vtk_output'
    output_dir = 'results'
    os.makedirs(output_dir, exist_ok=True)

    dt = 0.001          # time step used in the simulation
    output_interval = 0.01  # VTK output every 0.01 s  (every 10 steps)

    time_files = _collect_files(vtk_dir, dt=dt)
    if not time_files:
        print(f"No VTK files found in '{vtk_dir}'. Run the simulation first.")
        return

    all_times = [t for t, _ in time_files]
    print(f"Found {len(time_files)} snapshots: "
          f"t = {all_times[0]:.3f} … {all_times[-1]:.3f} s")

    # ------------------------------------------------------------------
    # 1) Four-panel snapshot figure: t = 0, 0.08, 0.14, 0.20 s
    # ------------------------------------------------------------------
    target_times = [0.0, 0.08, 0.14, 0.20]
    selected = []
    for tgt in target_times:
        closest_t, closest_fp = min(time_files, key=lambda tf: abs(tf[0] - tgt))
        selected.append((closest_t, closest_fp))

    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    for ax, (t, fp) in zip(axes, selected):
        d = _parse_vtk_binary(fp)
        _plot_snapshot(ax, d, t, vel_data=d)

    fig.suptitle('Sediment Dumping Simulation — Material Point Positions',
                 fontsize=13, fontweight='bold', y=1.01)
    fig.tight_layout()
    snap_path = os.path.join(output_dir, 'sediment_snapshots.png')
    fig.savefig(snap_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Snapshot figure saved → '{snap_path}'")

    # ------------------------------------------------------------------
    # 2) Centroid trajectory  (y vs t)
    # ------------------------------------------------------------------
    cy_list, vy_list, t_list = [], [], []
    for t, fp in time_files:
        d = _parse_vtk_binary(fp)
        pts = d['points']
        if pts.shape[0] == 0:
            continue
        cy_list.append(float(np.mean(pts[:, 1])))
        t_list.append(t)
        if 'MP_VELOCITY' in d:
            vy_list.append(float(np.mean(d['MP_VELOCITY'][:, 1])))

    # Analytical free-fall from y0 = 0.689 m with v0 = -0.154 m/s
    y0, v0, g = 0.689, -0.154, 9.81
    t_arr = np.linspace(0, max(t_list), 300)
    y_analytical = y0 + v0 * t_arr - 0.5 * g * t_arr**2

    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(12, 4))

    ax2a.plot(t_list, cy_list, 'b-o', markersize=3, linewidth=1.5,
              label='MPM centroid')
    ax2a.plot(t_arr, y_analytical, 'r--', linewidth=1.2,
              label='Free-fall (analytical)')
    ax2a.set_xlabel('Time (s)')
    ax2a.set_ylabel('Centroid y-position (m)')
    ax2a.set_title('Vertical descent of sediment centroid')
    ax2a.legend(fontsize=8)
    ax2a.grid(True, linestyle='--', linewidth=0.4)

    if vy_list:
        vy_analytical = v0 - g * np.array(t_list)
        ax2b.plot(t_list, vy_list, 'b-o', markersize=3, linewidth=1.5,
                  label='MPM mean vy')
        ax2b.plot(t_list, vy_analytical, 'r--', linewidth=1.2,
                  label='Free-fall vy (analytical)')
        ax2b.set_xlabel('Time (s)')
        ax2b.set_ylabel('Vertical velocity (m/s)')
        ax2b.set_title('Downward velocity of sediment centroid')
        ax2b.legend(fontsize=8)
        ax2b.grid(True, linestyle='--', linewidth=0.4)

    fig2.tight_layout()
    traj_path = os.path.join(output_dir, 'centroid_trajectory.png')
    fig2.savefig(traj_path, dpi=150)
    plt.close(fig2)
    print(f"Trajectory figure saved → '{traj_path}'")

    # ------------------------------------------------------------------
    # 3) Summary statistics
    # ------------------------------------------------------------------
    print("\n--- Summary ---")
    print(f"Initial centroid y : {cy_list[0]:.4f} m")
    print(f"Final   centroid y : {cy_list[-1]:.4f} m")
    print(f"Total descent      : {cy_list[0]-cy_list[-1]:.4f} m")
    if vy_list:
        print(f"Final velocity     : {vy_list[-1]:.3f} m/s")


if __name__ == '__main__':
    main()
