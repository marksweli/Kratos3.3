"""
Post-processing script for the sediment dumping simulation.

Reads VTK output files produced by Kratos MPMApplication and generates
snapshot images of the sediment blob position at specified time steps.

Usage:
    python3 post_process.py
"""

import os
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def parse_vtk_unstructured_grid(filepath):
    """Parse a legacy VTK unstructured grid file and return points and field data."""
    points = []
    mp_velocity = []
    mp_displacement = []
    mp_density = []

    with open(filepath, "r") as f:
        lines = f.readlines()

    i = 0
    n_points = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("POINTS"):
            parts = line.split()
            n_points = int(parts[1])
            i += 1
            while len(points) < n_points:
                vals = lines[i].strip().split()
                for k in range(0, len(vals), 3):
                    if len(points) < n_points:
                        points.append([float(vals[k]), float(vals[k+1])])
                i += 1
            continue

        if line.startswith("FIELD") or line.startswith("POINT_DATA"):
            i += 1
            continue

        if line.startswith("MP_VELOCITY"):
            parts = line.split()
            n = int(parts[1]) if len(parts) > 1 else n_points
            i += 1
            while len(mp_velocity) < n:
                vals = lines[i].strip().split()
                for k in range(0, len(vals), 3):
                    if len(mp_velocity) < n:
                        mp_velocity.append([float(vals[k]), float(vals[k+1])])
                i += 1
            continue

        if line.startswith("MP_DISPLACEMENT"):
            parts = line.split()
            n = int(parts[1]) if len(parts) > 1 else n_points
            i += 1
            while len(mp_displacement) < n:
                vals = lines[i].strip().split()
                for k in range(0, len(vals), 3):
                    if len(mp_displacement) < n:
                        mp_displacement.append([float(vals[k]), float(vals[k+1])])
                i += 1
            continue

        if line.startswith("MP_DENSITY"):
            parts = line.split()
            n = int(parts[1]) if len(parts) > 1 else n_points
            i += 1
            while len(mp_density) < n:
                vals = lines[i].strip().split()
                for v in vals:
                    if len(mp_density) < n:
                        mp_density.append(float(v))
                i += 1
            continue

        i += 1

    return {
        "points": np.array(points) if points else np.empty((0, 2)),
        "MP_VELOCITY": np.array(mp_velocity) if mp_velocity else np.empty((0, 2)),
        "MP_DISPLACEMENT": np.array(mp_displacement) if mp_displacement else np.empty((0, 2)),
        "MP_DENSITY": np.array(mp_density) if mp_density else np.empty(0),
    }


def get_time_from_filename(filename):
    """Extract simulation time from VTK filename (e.g. MPM_Material_0.08.vtk)."""
    match = re.search(r"_(\d+\.\d+)\.vtk$", filename)
    if match:
        return float(match.group(1))
    match = re.search(r"_(\d+)\.vtk$", filename)
    if match:
        return float(match.group(1))
    return None


def plot_snapshot(ax, data, time, domain=(0.7, 0.7)):
    """Plot a single time snapshot of material point positions."""
    pts = data["points"]
    disp = data["MP_DISPLACEMENT"]

    if pts.size == 0:
        ax.set_title(f"t = {time:.3f} s (no data)")
        return

    # If displacement data is available, compute current positions
    if disp.size > 0 and disp.shape[0] == pts.shape[0]:
        current_x = pts[:, 0] + disp[:, 0]
        current_y = pts[:, 1] + disp[:, 1]
    else:
        current_x = pts[:, 0]
        current_y = pts[:, 1]

    ax.set_facecolor("#d0e8f5")  # water color
    scatter = ax.scatter(
        current_x, current_y,
        c="#c8a05a", edgecolors="k", linewidths=0.3,
        s=60, zorder=3, label="Sediment MPs"
    )
    ax.set_xlim(0, domain[0])
    ax.set_ylim(0, domain[1])
    ax.set_xlabel("x (m)", fontsize=9)
    ax.set_ylabel("y (m)", fontsize=9)
    ax.set_title(f"t = {time:.3f} s", fontsize=10, fontweight="bold")
    ax.set_aspect("equal")
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)

    water_patch = mpatches.Patch(color="#d0e8f5", label="Water")
    sed_patch = mpatches.Patch(color="#c8a05a", label="Sediment")
    ax.legend(handles=[water_patch, sed_patch], fontsize=7, loc="lower right")


def main():
    vtk_dir = "vtk_output"
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    # Find all material-point VTK files (exclude background grid file)
    vtk_pattern = os.path.join(vtk_dir, "MPM_Material_*.vtk")
    vtk_files = sorted(glob.glob(vtk_pattern))

    if not vtk_files:
        print(f"No VTK files found in '{vtk_dir}'. Run the simulation first.")
        return

    # Map time -> filepath
    time_files = {}
    for fp in vtk_files:
        t = get_time_from_filename(os.path.basename(fp))
        if t is not None:
            time_files[t] = fp

    all_times = sorted(time_files.keys())
    print(f"Found {len(all_times)} VTK snapshots: t = {all_times[0]:.3f} .. {all_times[-1]:.3f} s")

    # --- Snapshot plot at t = 0.08 s and t = 0.14 s (or closest available) ---
    target_times = [0.08, 0.14]
    selected = []
    for tgt in target_times:
        closest = min(all_times, key=lambda t: abs(t - tgt))
        selected.append(closest)

    fig, axes = plt.subplots(1, len(selected), figsize=(5 * len(selected), 5))
    if len(selected) == 1:
        axes = [axes]

    for ax, t in zip(axes, selected):
        data = parse_vtk_unstructured_grid(time_files[t])
        plot_snapshot(ax, data, t)

    fig.suptitle("Sediment Dumping Simulation — MPM Snapshots", fontsize=12, fontweight="bold")
    fig.tight_layout()
    snapshot_path = os.path.join(output_dir, "sediment_snapshots.png")
    fig.savefig(snapshot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Snapshot image saved to '{snapshot_path}'")

    # --- Time series: centroid trajectory ---
    cx_list, cy_list = [], []
    t_list = []
    for t in all_times:
        data = parse_vtk_unstructured_grid(time_files[t])
        pts = data["points"]
        disp = data["MP_DISPLACEMENT"]
        if pts.size == 0:
            continue
        if disp.size > 0 and disp.shape[0] == pts.shape[0]:
            cx = np.mean(pts[:, 0] + disp[:, 0])
            cy = np.mean(pts[:, 1] + disp[:, 1])
        else:
            cx = np.mean(pts[:, 0])
            cy = np.mean(pts[:, 1])
        cx_list.append(cx)
        cy_list.append(cy)
        t_list.append(t)

    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(t_list, cy_list, "b-o", markersize=3, linewidth=1.5)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Centroid y-position (m)")
    ax2.set_title("Vertical descent of sediment centroid")
    ax2.grid(True, linestyle="--", linewidth=0.4)
    traj_path = os.path.join(output_dir, "centroid_trajectory.png")
    fig2.tight_layout()
    fig2.savefig(traj_path, dpi=150)
    plt.close(fig2)
    print(f"Centroid trajectory saved to '{traj_path}'")


if __name__ == "__main__":
    main()
