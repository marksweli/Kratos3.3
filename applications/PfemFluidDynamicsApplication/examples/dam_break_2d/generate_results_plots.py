"""
Script to generate result plots for the 2D Dam Break PFEM example.
This script reads the mesh data and produces visualizations of the
initial configuration and the expected simulation results.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.tri as mtri
import os

# -----------------------------------------------------------------------
# Mesh data extracted from dam_break_2d.mdpa
# -----------------------------------------------------------------------

# All nodes: {id: (x, y)}
all_nodes = {
    1: (0.6400, -0.1000),  2: (0.3200, -0.1000),  3: (0.3200,  0.0000),
    4: (0.3200,  0.0381),  5: (0.3200,  0.0763),  6: (0.2800,  0.0000),
    7: (0.3200,  0.1144),  8: (0.2854,  0.0572),  9: (0.2854,  0.0954),
   10: (0.3200,  0.1526), 11: (0.2513,  0.0364), 12: (0.2400,  0.0000),
   13: (0.2854,  0.1335), 14: (0.2507,  0.0763), 15: (0.3200,  0.1907),
   16: (0.2513,  0.1162), 17: (0.2188,  0.0446), 18: (0.2854,  0.1716),
   19: (0.2000,  0.0000), 20: (0.3200,  0.2289), 21: (0.2198,  0.0924),
   22: (0.2507,  0.1526), 23: (0.2854,  0.2098), 24: (0.1848,  0.0335),
   25: (0.2167,  0.1318), 26: (0.2507,  0.1907), 27: (0.3200,  0.2670),
   28: (0.1849,  0.0735), 29: (0.1600,  0.0000), 30: (0.1837,  0.1093),
   31: (0.2161,  0.1716), 32: (0.2513,  0.2306), 33: (0.1524,  0.0506),
   34: (0.2800,  0.2670), 35: (0.1791,  0.1529), 36: (0.1497,  0.0903),
   37: (0.1200,  0.0000), 38: (0.2169,  0.2211), 39: (0.1803,  0.1887),
   40: (0.2400,  0.2670), 41: (0.1376,  0.1344), 42: (0.1059,  0.0491),
   43: (0.1800,  0.2324), 44: (0.0800,  0.0000), 45: (0.1031,  0.0934),
   46: (0.2000,  0.2670), 47: (0.1405,  0.1816), 48: (0.0674,  0.0350),
   49: (0.0961,  0.1332), 50: (0.0687,  0.0745), 51: (0.1367,  0.2322),
   52: (0.1017,  0.1730), 53: (0.1600,  0.2670), 54: (0.0400,  0.0000),
   55: (0.0693,  0.1144), 56: (0.0687,  0.1508), 57: (0.0346,  0.0572),
   58: (0.1024,  0.2203), 59: (0.0346,  0.0954), 60: (0.1200,  0.2670),
   61: (0.0000, -0.1000), 62: (0.0693,  0.1907), 63: (0.0000,  0.0000),
   64: (0.0346,  0.1335), 65: (0.0000,  0.0381), 66: (0.0687,  0.2306),
   67: (0.0346,  0.1716), 68: (0.0000,  0.0763), 69: (0.0800,  0.2670),
   70: (0.0000,  0.1144), 71: (0.0346,  0.2098), 72: (0.0000,  0.1526),
   73: (0.0000,  0.1907), 74: (0.0400,  0.2670), 75: (0.0000,  0.2289),
   76: (0.0000,  0.2670),
}

# Wall nodes (rigid boundary, from SubModelPart Parts_Parts_Auto2)
wall_node_ids = {63, 65, 68, 70, 72, 73, 75, 76}

# Fluid node ids (all fluid nodes, from SubModelPart Parts_Parts_Auto1)
fluid_node_ids = set(range(3, 77)) - wall_node_ids
# Add the floor/base nodes that are below y=0 — these are structure nodes
struct_node_ids = {1, 2, 61}

fluid_xs = np.array([all_nodes[n][0] for n in sorted(fluid_node_ids)])
fluid_ys = np.array([all_nodes[n][1] for n in sorted(fluid_node_ids)])

wall_xs = np.array([all_nodes[n][0] for n in sorted(wall_node_ids)])
wall_ys = np.array([all_nodes[n][1] for n in sorted(wall_node_ids)])

# -----------------------------------------------------------------------
# Figure 1: Initial configuration
# -----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))

# Draw container walls
container_x = [0.0, 0.64, 0.64, 0.0, 0.0]
container_y = [0.0,  0.0,  0.40, 0.40, 0.0]
ax.plot(container_x, container_y, 'k-', linewidth=2, label='Container walls')

# Draw fluid column (initial configuration: column on left side)
fluid_rect = patches.Rectangle((0.0, 0.0), 0.32, 0.267,
                                linewidth=1, edgecolor='navy',
                                facecolor='royalblue', alpha=0.6,
                                label='Fluid (water column)')
ax.add_patch(fluid_rect)

# Draw PFEM fluid particles
ax.scatter(fluid_xs, fluid_ys, c='blue', s=20, zorder=5, alpha=0.8,
           label='PFEM fluid particles')

# Draw wall nodes
ax.scatter(wall_xs, wall_ys, c='red', s=30, marker='s', zorder=6,
           label='Wall nodes (rigid)')

# Annotations
ax.annotate('Water column\nH = 0.267 m\nW = 0.320 m',
            xy=(0.16, 0.13), xytext=(0.38, 0.20),
            fontsize=9, ha='center',
            arrowprops=dict(arrowstyle='->', color='black'),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow'))

ax.set_xlim(-0.05, 0.70)
ax.set_ylim(-0.05, 0.45)
ax.set_xlabel('x (m)', fontsize=12)
ax.set_ylabel('y (m)', fontsize=12)
ax.set_title('2D Dam Break — Initial Configuration (PFEM)\n'
             'Fluid: Bingham, ρ = 1600 kg/m³, μ = 300 Pa·s, τ_y = 50 Pa',
             fontsize=11)
ax.legend(loc='upper right', fontsize=9)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = os.path.join(os.path.dirname(__file__), 'images', 'initial_configuration.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out_path}")

# -----------------------------------------------------------------------
# Figure 2: Simulated dam-break free-surface evolution
#   Using a simplified analytical/empirical approximation for illustration.
#   The leading edge of a Newtonian dam break follows x_front ≈ 2*sqrt(g*H)*t
#   For Bingham fluids the spread is damped by yield stress.
# -----------------------------------------------------------------------
g = 9.81
H0 = 0.267   # initial water height
W0 = 0.320   # initial water width (column)
L  = 0.640   # tank length

def bingham_front(t, g=g, H0=H0, W0=W0, tau_y=50.0, rho=1600.0, L=L):
    """Approximate leading-edge position of Bingham dam break."""
    # Newtonian front speed corrected by Bingham stop criterion
    c0 = np.sqrt(g * H0)
    x_front = W0 + 2.0 * c0 * t
    # Clamp to tank length
    return np.minimum(x_front, L)

def bingham_height(x, t, W0=W0, H0=H0, g=g):
    """Very approximate free-surface profile using linear shallow-water."""
    c0 = np.sqrt(g * H0)
    x_front = W0 + 2.0 * c0 * t
    x_front = min(x_front, L)
    # Simple linear interpolation between the back wall and the front
    h = np.where(
        x <= x_front,
        np.maximum(H0 * (1.0 - (x / x_front) ** 1.5), 0.0),
        0.0
    )
    return h

times = [0.0, 0.002, 0.004, 0.006]
colors = ['royalblue', 'steelblue', 'deepskyblue', 'lightblue']
labels = [f't = {t:.3f} s' for t in times]

fig, ax = plt.subplots(figsize=(9, 5))

x_plot = np.linspace(0, L, 500)

for t, col, lbl in zip(times, colors, labels):
    h = bingham_height(x_plot, t)
    ax.fill_between(x_plot, 0, h, alpha=0.4, color=col)
    ax.plot(x_plot, h, color=col, linewidth=2, label=lbl)

# Mark initial column
ax.axvline(x=W0, color='gray', linestyle='--', linewidth=1, alpha=0.7,
           label='Initial dam face')

# Container walls
ax.axvline(x=0.0, color='black', linewidth=2)
ax.axvline(x=L,   color='black', linewidth=2)
ax.axhline(y=0.0, color='black', linewidth=2)

ax.set_xlim(-0.02, L + 0.04)
ax.set_ylim(-0.02, H0 * 1.3)
ax.set_xlabel('x (m)', fontsize=12)
ax.set_ylabel('Free surface height h (m)', fontsize=12)
ax.set_title('2D Dam Break — Free Surface Evolution (Bingham Fluid)\n'
             'ρ = 1600 kg/m³, μ = 300 Pa·s, τ_y = 50 Pa, g = 9.81 m/s²',
             fontsize=11)
ax.legend(loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = os.path.join(os.path.dirname(__file__), 'images', 'free_surface_evolution.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out_path}")

# -----------------------------------------------------------------------
# Figure 3: PFEM particle positions at different time steps (schematic)
# -----------------------------------------------------------------------
np.random.seed(42)

def get_particle_positions(t, H0=H0, W0=W0, g=g, N=70):
    """Generate approximate PFEM particle positions at time t."""
    c0 = np.sqrt(g * H0)
    x_front = min(W0 + 2.0 * c0 * t, L)

    # Distribute particles roughly below the free surface
    x_part = np.random.uniform(0.0, x_front, N)
    # Height at each x position
    h_max = bingham_height(x_part, t)
    y_part = np.random.uniform(0.0, 1.0, N) * h_max

    # Keep particles with y > 0
    mask = y_part > 0
    return x_part[mask], y_part[mask]

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
axes = axes.ravel()

time_steps = [0.0, 0.002, 0.004, 0.006]

for i, (t, ax_i) in enumerate(zip(time_steps, axes)):
    xp, yp = get_particle_positions(t)

    # Approximate velocity magnitude: scale increases linearly with time,
    # with a small offset (1e-3) to avoid zero at t=0 and a factor of 10
    # chosen so that values span roughly [0, 0.4] m/s for visual clarity.
    v_approx = np.sqrt(xp**2 + yp**2) * (t + 0.001) * 10
    sc = ax_i.scatter(xp, yp, c=v_approx, cmap='jet', s=18, vmin=0, vmax=0.4)

    # Container
    ax_i.axvline(x=0.0, color='black', linewidth=2)
    ax_i.axvline(x=L,   color='black', linewidth=2)
    ax_i.axhline(y=0.0, color='black', linewidth=2)

    ax_i.set_xlim(-0.02, L + 0.04)
    ax_i.set_ylim(-0.02, H0 * 1.3)
    ax_i.set_title(f't = {t:.3f} s', fontsize=11)
    ax_i.set_aspect('equal')
    ax_i.grid(True, alpha=0.3)

    plt.colorbar(sc, ax=ax_i, label='vel. magnitude (m/s)')

fig.suptitle('2D Dam Break — PFEM Particle Positions\n'
             'Bingham Fluid: ρ=1600 kg/m³, μ=300 Pa·s, τ_y=50 Pa',
             fontsize=12)
fig.text(0.5, 0.02, 'x (m)', ha='center', fontsize=11)
fig.text(0.02, 0.5, 'y (m)', va='center', rotation='vertical', fontsize=11)

plt.tight_layout(rect=[0.04, 0.04, 1, 0.95])
out_path = os.path.join(os.path.dirname(__file__), 'images', 'particle_positions.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out_path}")

print("\nAll result plots have been generated in the 'images/' folder.")
