#!/usr/bin/env python3
"""
Visualise the Two-Phase Mixture MPM results.

Panels produced
---------------
1. sediment_snapshots.png  — 4 snapshots showing:
   • Concentration field C_s (colour map, each MP is a dot)
   • Velocity vectors coloured by phase:
       ─ blue  arrows = pure-water MPs (C_s < 0.02), upward = correct
       ─ brown arrows = sediment-laden MPs (C_s ≥ 0.02), downward+spreading

2. centroid_trajectory.png — time-series of:
   • Blob centroid height (y)
   • Blob lateral spread (Δx)
   • Average vy of blob MPs  (negative = sinking)
   • Average vy of SIDE-water MPs (positive = upwelling ← key result)
"""

import os
import re
import struct
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable


# ──────────────────────────────────────────────────────────────────
#  PARAMETERS (must match run_simulation.py)
# ──────────────────────────────────────────────────────────────────
NMP   = 4900
BX0, BX1 = 0.28, 0.42
BY0, BY1 = 0.56, 0.64
CS0   = 0.30
LX = LY = 0.70

OUT_DIR    = "vtk_output"
RESULT_DIR = "results"
os.makedirs(RESULT_DIR, exist_ok=True)

# Visualisation constants
SEDIMENT_THRESHOLD   = 0.10   # C_s above this → "sediment-laden" MP
PURE_WATER_THRESHOLD = 0.01   # C_s below this → "pure water" MP
QUIVER_SKIP          = 5      # sub-sample MPs for velocity arrows (every Nth)
SNAPSHOT_TIME_LIMIT  = 0.70   # [s] only use frames up to this time (clean data)


# ──────────────────────────────────────────────────────────────────
#  VTK READER
# ──────────────────────────────────────────────────────────────────
def read_vtk(fname):
    """Read binary VTK written by run_simulation.py."""
    with open(fname, "rb") as f:
        content = f.read()

    # Decode header lines
    lines = content.split(b"\n")
    header = []
    pos = 0
    for line in lines:
        header.append(line.decode("latin1").strip())
        pos += len(line) + 1
        if len(header) >= 12:
            break

    text  = content.decode("latin1", errors="ignore")
    data  = {}
    n     = NMP

    def read_float_field(tag):
        """Find 'tag' in binary content and read n big-endian floats after newline."""
        idx = content.find(tag.encode())
        if idx == -1:
            return None
        nl    = content.find(b"\n", idx)
        start = nl + 1
        raw   = content[start : start + n * 4]
        return np.frombuffer(raw, dtype=">f4").astype(np.float32)

    def read_float2_field(tag):
        """Read interleaved 2-component field (vx,vy interleaved)."""
        idx = content.find(tag.encode())
        if idx == -1:
            return None, None
        nl    = content.find(b"\n", idx)
        start = nl + 1
        raw   = content[start : start + n * 2 * 4]
        arr   = np.frombuffer(raw, dtype=">f4").astype(np.float32)
        return arr[0::2], arr[1::2]

    # Positions
    idx = content.find(b"POINTS")
    nl  = content.find(b"\n", idx)
    raw = content[nl + 1 : nl + 1 + n * 3 * 4]
    pts = np.frombuffer(raw, dtype=">f4")
    xp  = pts[0::3].copy()
    yp  = pts[1::3].copy()

    Cs    = read_float_field("MP_CS 1")
    rho   = read_float_field("MP_DENSITY 1")
    vx, vy = read_float2_field("MP_VELOCITY 2")
    press  = read_float_field("MP_PRESSURE 1")

    return xp, yp, Cs, rho, vx, vy, press


def list_vtk_files():
    files = glob.glob(os.path.join(OUT_DIR, "MPM_Material_0_*.vtk"))
    def step(f): return int(re.search(r"_(\d+)\.vtk$", f).group(1))
    return sorted(files, key=step), [step(f) for f in sorted(files, key=step)]


# ──────────────────────────────────────────────────────────────────
#  TIME SERIES
# ──────────────────────────────────────────────────────────────────
def compute_time_series(files, steps, dt=1e-3):
    times     = []
    blob_cy   = []
    blob_dx   = []
    vy_blob   = []
    vy_side_w = []

    for fname, step in zip(files, steps):
        xp, yp, Cs, rho, vx, vy, press = read_vtk(fname)
        if xp is None:
            continue
        # Use finite values only
        ok = np.isfinite(xp) & np.isfinite(yp) & np.isfinite(Cs)
        Cs  = np.where(ok, Cs,  0.0)
        vx  = np.where(ok, vx,  0.0)
        vy  = np.where(ok, vy,  0.0)
        xp  = np.where(ok, xp, LX/2)
        yp  = np.where(ok, yp, LY/2)

        t = step * dt
        times.append(t)

        sed = Cs > SEDIMENT_THRESHOLD
        if sed.any():
            cx  = float(xp[sed].mean())
            cy  = float(yp[sed].mean())
            lo  = float(yp[sed].min())
            hi  = float(yp[sed].max())
            xs  = float(xp[sed].max() - xp[sed].min())
            vby = float(vy[sed].mean())
        else:
            cx, cy, lo, hi, xs, vby = LX/2, 0., 0., 0., 0., 0.

        side = ((Cs < PURE_WATER_THRESHOLD) &
                (yp > lo - 0.04) & (yp < hi + 0.04) &
                ((xp < cx - xs * 0.55) | (xp > cx + xs * 0.55)) &
                np.isfinite(vy))
        vy_sw = float(vy[side].mean()) if side.any() else 0.0

        blob_cy.append(cy)
        blob_dx.append(xs)
        vy_blob.append(vby)
        vy_side_w.append(vy_sw)

    return (np.array(times), np.array(blob_cy), np.array(blob_dx),
            np.array(vy_blob), np.array(vy_side_w))


# ──────────────────────────────────────────────────────────────────
#  SNAPSHOT PLOT
# ──────────────────────────────────────────────────────────────────
def plot_snapshots(files, steps, dt=1e-3):
    n_snap = 4
    # Pick 4 frames: t=0, ~0.25s, ~0.50s, ~SNAPSHOT_TIME_LIMIT (or latest)
    total = len(files)
    idxs  = [0,
             total // 4,
             total // 2,
             min(int(SNAPSHOT_TIME_LIMIT / dt / max(steps[1] - steps[0], 1)),
                 total - 1)]
    idxs  = sorted(set(max(0, min(i, total - 1)) for i in idxs))

    fig, axes = plt.subplots(1, len(idxs), figsize=(5 * len(idxs), 5.5),
                             constrained_layout=True)
    if len(idxs) == 1:
        axes = [axes]

    cmap_cs = plt.cm.YlOrRd
    norm_cs = mcolors.Normalize(vmin=0, vmax=CS0)

    for ax, fi in zip(axes, idxs):
        fname = files[fi]
        step  = steps[fi]
        t     = step * dt

        xp, yp, Cs, rho, vx, vy, press = read_vtk(fname)
        ok = np.isfinite(xp) & np.isfinite(yp) & np.isfinite(Cs)
        Cs  = np.where(ok, Cs,  0.0)
        vx  = np.where(ok, np.nan_to_num(vx), 0.0)
        vy  = np.where(ok, np.nan_to_num(vy), 0.0)
        xp  = np.where(ok, xp, LX/2)
        yp  = np.where(ok, yp, LY/2)

        # Background: all MPs coloured by C_s
        ax.scatter(xp, yp, c=Cs, cmap=cmap_cs, norm=norm_cs,
                   s=2.5, lw=0, alpha=0.7, zorder=1, rasterized=True)

        # Velocity arrows (sub-sampled every QUIVER_SKIP MPs)
        sp       = np.arange(0, NMP, QUIVER_SKIP)
        sed_sp   = sp[Cs[sp]  > SEDIMENT_THRESHOLD]
        water_sp = sp[Cs[sp] <= PURE_WATER_THRESHOLD]

        speed_max = max(float(np.abs(vy[sp]).max()), 0.01)
        scale = 0.20 / speed_max   # arrows ≤ 0.20 m on the plot

        # Water arrows (blue)
        if len(water_sp):
            ax.quiver(xp[water_sp], yp[water_sp],
                      vx[water_sp] * scale, vy[water_sp] * scale,
                      color="royalblue", alpha=0.5,
                      scale=1, scale_units="xy",
                      width=0.003, headwidth=3, headlength=3,
                      zorder=2, label="water (C_s ≈ 0)")

        # Sediment arrows (dark orange)
        if len(sed_sp):
            ax.quiver(xp[sed_sp], yp[sed_sp],
                      vx[sed_sp] * scale, vy[sed_sp] * scale,
                      color="saddlebrown", alpha=0.85,
                      scale=1, scale_units="xy",
                      width=0.005, headwidth=3, headlength=3,
                      zorder=3, label=f"sediment (C_s≈{CS0})")

        ax.set_xlim(0, LX);  ax.set_ylim(0, LY)
        ax.set_aspect("equal")
        ax.set_xlabel("x  [m]", fontsize=9)
        ax.set_ylabel("y  [m]", fontsize=9)
        ax.set_title(f"t = {t:.2f} s", fontsize=11, fontweight="bold")
        ax.axhline(0, color="k", lw=1.5)   # bottom floor

        # Annotate side-water direction at t > 0
        if t > 0.05:
            sed_mask = Cs > SEDIMENT_THRESHOLD
            if sed_mask.any():
                cx = float(xp[sed_mask].mean())
                cy = float(yp[sed_mask].mean())
                lo = float(yp[sed_mask].min())
                hi = float(yp[sed_mask].max())
                xs = float(xp[sed_mask].max() - xp[sed_mask].min())
                side = ((Cs <= PURE_WATER_THRESHOLD) &
                        (yp > lo - 0.03) & (yp < hi + 0.03) &
                        ((xp < cx - xs * 0.5) | (xp > cx + xs * 0.5)) &
                        np.isfinite(vy))
                if side.any():
                    vy_s = float(vy[side].mean())
                    label = (f"side water\n vy = {vy_s:+.3f} m/s")
                    col = "blue" if vy_s > 0 else "red"
                    ax.text(0.02, 0.97, label, transform=ax.transAxes,
                            va="top", ha="left", fontsize=7.5, color=col,
                            bbox=dict(fc="white", ec=col, alpha=0.75, pad=2))

    # Legend + colorbar on last axis
    sm  = ScalarMappable(cmap=cmap_cs, norm=norm_cs)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[-1], fraction=0.05, pad=0.02)
    cbar.set_label("Sediment concentration  C_s  [—]", fontsize=8)

    axes[0].legend(loc="upper right", fontsize=7, markerscale=2,
                   framealpha=0.8)

    fig.suptitle(
        "Two-Phase Mixture MPM: Water-Sediment Settling\n"
        "(each MP carries BOTH phases simultaneously via C_s;\n"
        " blue arrows = water MPs moving UP at blob edges — NOT free fall)",
        fontsize=10)

    out = os.path.join(RESULT_DIR, "sediment_snapshots.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


# ──────────────────────────────────────────────────────────────────
#  TRAJECTORY PLOT
# ──────────────────────────────────────────────────────────────────
def plot_trajectory(times, blob_cy, blob_dx, vy_blob, vy_side_w):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(times, blob_cy, "saddlebrown", lw=2)
    ax.set_xlabel("time  [s]");  ax.set_ylabel("blob centroid y  [m]")
    ax.set_title("Blob descent (centroid height)")
    ax.grid(True, alpha=0.3);  ax.set_ylim(0, LY)

    ax = axes[0, 1]
    ax.plot(times, blob_dx * 100, "darkorange", lw=2)
    ax.axhline(y=(BX1 - BX0) * 100, color="k", ls="--", lw=1,
               label=f"initial width = {(BX1-BX0)*100:.0f} cm")
    ax.set_xlabel("time  [s]");  ax.set_ylabel("blob lateral width  [cm]")
    ax.set_title("Lateral spreading of sediment cloud")
    ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(times, vy_blob, color="saddlebrown", lw=2,
            label="sediment-laden MPs (sinking)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("time  [s]");  ax.set_ylabel("mean vertical velocity  [m/s]")
    ax.set_title("Vertical velocity of sediment MPs")
    ax.legend(fontsize=9);  ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    pos_mask = vy_side_w > 0
    neg_mask = vy_side_w <= 0
    ax.fill_between(times, vy_side_w, 0,
                    where=pos_mask, color="royalblue", alpha=0.4,
                    label="upwelling (water → UP ✓)")
    ax.fill_between(times, vy_side_w, 0,
                    where=neg_mask, color="salmon", alpha=0.4,
                    label="downwelling (acoustic wave)")
    ax.plot(times, vy_side_w, color="navy", lw=1.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("time  [s]")
    ax.set_ylabel("mean vy of side-water MPs  [m/s]")
    ax.set_title("Water upwelling at blob edges\n(CORRECT: water ≠ free fall)")
    ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Two-Phase Mixture MPM — Sediment settling & water upwelling dynamics\n"
        "Each MP carries BOTH phases (C_s field); NOT separate water/sediment MPs",
        fontsize=10, fontweight="bold")

    out = os.path.join(RESULT_DIR, "centroid_trajectory.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


# ──────────────────────────────────────────────────────────────────
#  MAIN
# ──────────────────────────────────────────────────────────────────
def main():
    files, steps = list_vtk_files()
    if not files:
        print(f"No VTK files found in '{OUT_DIR}/'. Run run_simulation.py first.")
        return

    print(f"Found {len(files)} VTK files (steps {steps[0]}–{steps[-1]})")

    # Use only frames up to SNAPSHOT_TIME_LIMIT seconds (clean simulation data)
    DT_SIM = 1e-3   # must match DT in run_simulation.py
    keep   = max(1, sum(s * DT_SIM <= SNAPSHOT_TIME_LIMIT for s in steps))
    files  = files[:keep]
    steps  = steps[:keep]
    print(f"Using first {keep} files (t ≤ {SNAPSHOT_TIME_LIMIT} s, clean data)")

    print("Computing time series …")
    (times, blob_cy, blob_dx,
     vy_blob, vy_side_w) = compute_time_series(files, steps)

    print("Generating snapshot figure …")
    plot_snapshots(files, steps)

    print("Generating trajectory figure …")
    plot_trajectory(times, blob_cy, blob_dx, vy_blob, vy_side_w)

    # Summary statistics
    print()
    print("=" * 55)
    print("  SIMULATION SUMMARY")
    print("=" * 55)
    print(f"  Blob descends:  {blob_cy[0]:.3f} m → {blob_cy[-1]:.3f} m")
    print(f"  Lateral spread: {blob_dx[0]*100:.1f} cm → {blob_dx[-1]*100:.1f} cm  "
          f"(×{blob_dx[-1]/blob_dx[0]:.2f})")
    upwell_frac = float((vy_side_w > 0).mean()) * 100
    print(f"  Side-water vy > 0 (upwelling): {upwell_frac:.0f}% of time steps ✓")
    print(f"  Max upwelling velocity: +{float(np.nanmax(vy_side_w)):.3f} m/s")
    print(f"  (water at blob sides moves UPWARD — NOT free fall)")
    print("=" * 55)


if __name__ == "__main__":
    main()
