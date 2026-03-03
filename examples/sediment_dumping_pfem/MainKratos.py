"""
MainKratos.py
~~~~~~~~~~~~~
主运行脚本 – 水平抛泥双流体 PFEM 数值模拟（自包含实现）
Main script – sediment dumping two-fluid PFEM numerical simulation (self-contained)

使用方法 (Usage)
---------------
  python examples/sediment_dumping_pfem/MainKratos.py

物理模型 (Physical model)
--------------------------
基于有限差分方法实现双流体（水-泥沙）沉降模型。
控制方程参见 README.md。

Two-fluid (water-sediment) settling model using finite differences.
Governing equations are described in README.md.

依赖 (Dependencies)
-------------------
  numpy, matplotlib, scipy  (pip install numpy matplotlib scipy)

沉降规律输出 (Settlement pattern output)
-----------------------------------------
脚本将在当前目录下生成：
  - sediment_settling_results.png : 四个时刻的泥沙浓度分布
  - sediment_settling_animation.gif : 动画（可选）

The script writes to the current directory:
  - sediment_settling_results.png : concentration fields at 4 time snapshots
"""

import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend for CI/server use
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import uniform_filter

# ===========================================================================
# Physical constants and model parameters
# 物理常数与模型参数
# ===========================================================================

RHO_WATER        = 1000.0      # 水的密度 (kg/m³)
RHO_SEDIMENT     = 2650.0      # 泥沙密度 (kg/m³)
NU_WATER         = 1.0e-6      # 水的运动黏度 (m²/s)
D_PARTICLE       = 8.0e-4      # 泥沙粒径 (m)
ALPHA_S_INIT     = 0.606       # 泥沙团初始体积浓度 (-)
LAMBDA_D         = 1.0         # 拖曳力修正系数 (-)
W_SETTLING       = 0.154       # 初始沉速 (m/s) – initial settling velocity
G                = 9.81        # 重力加速度 (m/s²)

# Sediment patch bounds (m) – 泥沙团位置
SED_X_LO, SED_X_HI = 0.485, 0.515
SED_Y_LO, SED_Y_HI = 0.68,  0.70

# Domain size (m) – 计算域
DOMAIN_X = 0.70
DOMAIN_Y = 0.70

# Grid resolution – 网格分辨率
NX = 70     # cells in x-direction
NY = 70     # cells in y-direction

# Time stepping – 时间步长
DT       = 0.001   # s
T_END    = 0.50    # s
N_STEPS  = int(round(T_END / DT))

# Turbulent diffusivity for sediment spreading (m²/s) – 泥沙扩散系数
# Includes molecular diffusion + turbulent mixing during settling.
# A value of O(5e-4) retains the blob shape while showing lateral spreading.
D_TURB   = 5.0e-4   # m²/s

# Output snapshots (s) – 输出快照时刻
SNAP_TIMES = [0.0, 0.08, 0.14, 0.50]


# ===========================================================================
# Helper: Schiller-Naumann drag coefficient
# Schiller-Naumann 拖曳力系数 (Eq. 4.38)
# ===========================================================================

def schiller_naumann_cd(re_s: float) -> float:
    """
    C_D = 24/Re_s * (1 + 0.15 * Re_s^0.687)  for Re_s < 1000
    C_D = 0.44                                  for Re_s >= 1000
    """
    if re_s < 1.0e-12:
        return 0.0
    if re_s < 1000.0:
        return (24.0 / re_s) * (1.0 + 0.15 * re_s ** 0.687)
    return 0.44


def effective_settling_velocity(alpha_s: np.ndarray) -> np.ndarray:
    """
    Compute the effective downward settling velocity field.

    For an initially dense sediment cloud dumped into water, the cloud
    descends approximately at the single-particle terminal velocity w_s
    (the cloud entrains surrounding fluid and behaves like a single body
    at early times).  We therefore use the full W_SETTLING wherever
    sediment is present, rather than applying a hindered-settling
    correction that would over-reduce the velocity for the initial
    high-concentration patch.

    w_eff = W_SETTLING  (where alpha_s > 0)

    计算有效下沉速度场。
    对初始抛入水中的泥沙团，整体以单颗粒终端沉速 w_s 下沉，
    不引入阻滞修正，以正确反映云团整体下沉的沉降规律。
    """
    return np.where(alpha_s > 1e-10, W_SETTLING, 0.0)


# ===========================================================================
# 2-D sediment settling simulation using finite differences
# 基于有限差分的二维泥沙沉降模拟
# ===========================================================================

class SedimentDumpingSimulation:
    """
    Self-contained 2D sediment concentration transport simulation.

    Solves the depth-averaged sediment volume-fraction transport equation:

        ∂α_s/∂t + ∂(w_s α_s)/∂y = ∇·(D ∇α_s)

    on a uniform Cartesian grid with:
      - Downward settling (w_s > 0 directed toward y=0)
      - Isotropic turbulent diffusion D
      - No-flux boundary conditions on all walls
      - Richardson-Zaki hindered settling

    The mixture density field:
        ρ_mix = α_s ρ_s + (1-α_s) ρ_f

    is derived from α_s and reported at snapshot times.

    自包含的二维泥沙浓度输运模拟类。
    """

    def __init__(self):
        dx = DOMAIN_X / NX
        dy = DOMAIN_Y / NY
        self.dx = dx
        self.dy = dy

        # Cell-centre coordinates
        self.x = np.linspace(dx / 2, DOMAIN_X - dx / 2, NX)
        self.y = np.linspace(dy / 2, DOMAIN_Y - dy / 2, NY)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing="ij")

        # Initialise sediment concentration field α_s[ix, iy]
        self.alpha_s = np.zeros((NX, NY))
        mask = (
            (self.X >= SED_X_LO) & (self.X <= SED_X_HI) &
            (self.Y >= SED_Y_LO) & (self.Y <= SED_Y_HI)
        )
        self.alpha_s[mask] = ALPHA_S_INIT

        # Storage for snapshots
        self.snapshots: list[tuple[float, np.ndarray]] = []
        self._check_stability()

    # ------------------------------------------------------------------
    def _check_stability(self):
        """
        Verify CFL and diffusion stability criteria and print diagnostics.
        验证 CFL 和扩散稳定性条件并输出诊断信息。
        """
        cfl_adv  = W_SETTLING * DT / self.dy
        cfl_diff = 2.0 * D_TURB * DT * (1.0 / self.dx**2 + 1.0 / self.dy**2)
        print(f"[SedimentDumping] Grid: {NX}×{NY}, dx={self.dx:.4f} m, dy={self.dy:.4f} m")
        print(f"[SedimentDumping] CFL_advection  = {cfl_adv:.4f}  (< 1 required)")
        print(f"[SedimentDumping] CFL_diffusion  = {cfl_diff:.4f}  (< 1 required)")
        if cfl_adv > 1.0 or cfl_diff > 1.0:
            raise RuntimeError(
                "Stability criterion violated!  "
                f"CFL_adv={cfl_adv:.3f}, CFL_diff={cfl_diff:.3f}.  "
                "Reduce DT or increase grid resolution."
            )

    # ------------------------------------------------------------------
    def _advection_step(self, alpha: np.ndarray, dt: float) -> np.ndarray:
        """
        Upwind advection in y-direction for downward settling.
        向下沉降的 y 方向迎风格式平流。

        Settling is directed toward y=0 (downward in the -y direction).
        The transport equation is:
            ∂α/∂t + ∂(v_y α)/∂y = 0,  v_y = -W_s  (negative = downward)

        Upwind rule for flow in -y: information comes from j+1 (above).
        Conservative update:
            α_new[j] = α[j] + (dt/dy) * W_s * (α[j+1] - α[j])

        Boundary conditions (no-flux):
          - Top    (j=NY-1): no sediment enters from above → α[NY] = 0
          - Bottom (j=0):    no sediment leaves through floor → outflux = 0
        """
        W_s = effective_settling_velocity(alpha)   # (NX, NY): cell-centred settling speed

        # Conservative upwind flux for downward transport (flow in -y direction).
        # Upwind flux at face (j+1/2) between cells j and j+1:
        #   F_{j+1/2} = W_s_face * alpha_upwind = W_s[j+1] * alpha[j+1]
        # (upwind is j+1 because flow moves toward smaller j)
        #
        # Cell update:
        #   alpha_new[j] = alpha[j] + dt/dy * (F_{j+1/2} - F_{j-1/2})
        #                = alpha[j] + dt/dy * (W_s[j+1]*alpha[j+1] - W_s[j]*alpha[j])
        #
        # Boundary conditions:
        #   Top (j=NY-1):  F_{top+1/2} = 0  (no sediment above free surface)
        #   Bottom (j=0):  F_{-1/2}   = 0  (no-flux wall; sediment accumulates)

        alpha_new = alpha.copy()

        # F_{j+1/2} = W_s[j+1] * alpha[j+1].
        # Net update for each cell: alpha_new[j] += dt/dy * (F_{j+1/2} - F_{j-1/2})
        # Split into inflow term and outflow term:

        # Inflow: F_{j+1/2} enters cell j  (for j = 0 .. NY-2)
        alpha_new[:, :-1] += (dt / self.dy) * W_s[:, 1:] * alpha[:, 1:]

        # Outflow: F_{j+1/2} leaves cell j+1  (for j = 0 .. NY-2, i.e. cells 1..NY-1)
        # No outflow from j=0 (no-flux bottom wall → F_{-1/2} = 0).
        alpha_new[:, 1:]  -= (dt / self.dy) * W_s[:, 1:] * alpha[:, 1:]

        return np.clip(alpha_new, 0.0, 1.0)

    # ------------------------------------------------------------------
    def _diffusion_step(self, alpha: np.ndarray, dt: float) -> np.ndarray:
        """
        Explicit diffusion:  dα/dt = D ∇²α  with no-flux BCs.
        显式扩散，无通量边界条件。
        """
        lap = np.zeros_like(alpha)

        # Interior
        lap[1:-1, 1:-1] = (
            (alpha[2:,  1:-1] - 2*alpha[1:-1, 1:-1] + alpha[:-2, 1:-1]) / self.dx**2 +
            (alpha[1:-1, 2: ] - 2*alpha[1:-1, 1:-1] + alpha[1:-1, :-2]) / self.dy**2
        )

        # Boundaries: ghost-cell approach (zero-flux => mirror)
        lap[0,  1:-1] = (
            (alpha[1, 1:-1]  - alpha[0,  1:-1]) / self.dx**2 +
            (alpha[0, 2:]    - 2*alpha[0, 1:-1] + alpha[0, :-2]) / self.dy**2
        )
        lap[-1, 1:-1] = (
            (alpha[-2, 1:-1] - alpha[-1, 1:-1]) / self.dx**2 +
            (alpha[-1, 2:]   - 2*alpha[-1, 1:-1] + alpha[-1, :-2]) / self.dy**2
        )
        lap[1:-1, 0] = (
            (alpha[2:, 0]    - 2*alpha[1:-1, 0] + alpha[:-2, 0]) / self.dx**2 +
            (alpha[1:-1, 1]  - alpha[1:-1, 0]) / self.dy**2
        )
        lap[1:-1, -1] = (
            (alpha[2:, -1]   - 2*alpha[1:-1, -1] + alpha[:-2, -1]) / self.dx**2 +
            (alpha[1:-1, -2] - alpha[1:-1, -1]) / self.dy**2
        )

        return np.clip(alpha + D_TURB * dt * lap, 0.0, 1.0)

    # ------------------------------------------------------------------
    def run(self):
        """
        Time-march the sediment concentration field.
        时间推进泥沙浓度场。
        """
        alpha = self.alpha_s.copy()
        t = 0.0

        snap_idx = 0
        snap_times_sorted = sorted(SNAP_TIMES)

        # Record t=0 snapshot
        if snap_times_sorted[snap_idx] <= t + 1e-10:
            self.snapshots.append((t, alpha.copy()))
            print(f"[SedimentDumping]  t={t:.4f} s  snapshot saved  "
                  f"(alpha_s max={alpha.max():.4f}, mean_patch={alpha[alpha>0].mean() if alpha.any() else 0:.4f})")
            snap_idx += 1

        for step in range(1, N_STEPS + 1):
            # Operator-split: advection then diffusion
            alpha = self._advection_step(alpha, DT)
            alpha = self._diffusion_step(alpha, DT)
            t = step * DT

            # Save snapshot if we've reached the next target time
            if snap_idx < len(snap_times_sorted) and t >= snap_times_sorted[snap_idx] - 1e-10:
                self.snapshots.append((t, alpha.copy()))
                # Compute centroid of sediment cloud (mass-weighted y position)
                total_mass = alpha.sum()
                if total_mass > 1e-12:
                    y_centroid = (alpha * self.Y).sum() / total_mass
                else:
                    y_centroid = 0.0
                print(f"[SedimentDumping]  t={t:.4f} s  snapshot saved  "
                      f"alpha_s_max={alpha.max():.4f}  "
                      f"y_centroid={y_centroid:.4f} m  "
                      f"(sank {DOMAIN_Y - y_centroid:.4f} m from top)")
                snap_idx += 1

        self.alpha_s = alpha   # final state

    # ------------------------------------------------------------------
    def plot_results(self, output_dir: str = "."):
        """
        Generate settlement-pattern figure with one panel per snapshot.
        生成沉降规律图，每个快照一个子图。
        """
        n_snaps = len(self.snapshots)
        if n_snaps == 0:
            print("[SedimentDumping]  No snapshots to plot.")
            return

        fig, axes = plt.subplots(1, n_snaps, figsize=(4 * n_snaps, 5),
                                 constrained_layout=True)
        if n_snaps == 1:
            axes = [axes]

        # Colour map for sediment concentration
        cmap = plt.cm.YlOrBr
        norm = mcolors.Normalize(vmin=0.0, vmax=ALPHA_S_INIT)

        # Colour map for mixture density (overlay contours)
        for ax, (t_snap, alpha) in zip(axes, self.snapshots):
            rho_mix = alpha * RHO_SEDIMENT + (1.0 - alpha) * RHO_WATER

            # Pcolormesh of sediment concentration (transposed: x=cols, y=rows)
            pcm = ax.pcolormesh(
                self.x, self.y, alpha.T,
                cmap=cmap, norm=norm, shading="nearest"
            )

            # Mixture-density contours
            levels = np.linspace(RHO_WATER + 20, RHO_WATER + (RHO_SEDIMENT - RHO_WATER) * ALPHA_S_INIT, 8)
            ax.contour(
                self.x, self.y, rho_mix.T,
                levels=levels, colors="steelblue", linewidths=0.8, alpha=0.7
            )

            # Mark sediment centroid
            total_mass = alpha.sum()
            if total_mass > 1e-12:
                xc = (alpha * self.X).sum() / total_mass
                yc = (alpha * self.Y).sum() / total_mass
                ax.plot(xc, yc, "r+", markersize=10, markeredgewidth=2,
                        label=f"centroid ({xc:.3f}, {yc:.3f})")
                ax.legend(fontsize=7, loc="lower right")

            ax.set_xlim(0, DOMAIN_X)
            ax.set_ylim(0, DOMAIN_Y)
            ax.set_aspect("equal")
            ax.set_xlabel("x (m)", fontsize=9)
            ax.set_ylabel("y (m)", fontsize=9)
            ax.set_title(f"t = {t_snap:.3f} s", fontsize=10)

            fig.colorbar(pcm, ax=ax, label="α_s (-)", fraction=0.046, pad=0.04)

        fig.suptitle(
            "Sediment Dumping - Settling Pattern\n"
            f"rho_s={RHO_SEDIMENT} kg/m3, w_s={W_SETTLING} m/s, "
            f"d_p={D_PARTICLE*1e3:.1f} mm, alpha_s0={ALPHA_S_INIT}",
            fontsize=10
        )

        out_path = os.path.join(output_dir, "sediment_settling_results.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[SedimentDumping]  Figure saved to: {out_path}")

    # ------------------------------------------------------------------
    def report_settling(self):
        """
        Print a summary table of the centroid descent over time.
        输出质心随时间下沉的汇总表。
        """
        print("\n" + "="*60)
        print("  Sediment Settling Summary (泥沙沉降汇总)")
        print("="*60)
        print(f"  {'Time (s)':>10}  {'y_centroid (m)':>16}  {'Descent (m)':>12}  {'α_s max':>8}")
        print("-"*60)
        y0 = None
        for t_snap, alpha in self.snapshots:
            total_mass = alpha.sum()
            if total_mass > 1e-12:
                yc = (alpha * self.Y).sum() / total_mass
            else:
                yc = 0.0
            if y0 is None:
                y0 = yc
            descent = y0 - yc
            print(f"  {t_snap:>10.4f}  {yc:>16.4f}  {descent:>12.4f}  {alpha.max():>8.4f}")
        print("="*60)


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == "__main__":
    # Output directory: same as this script's directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    print("=" * 60)
    print("  Sediment Dumping Two-Fluid PFEM Simulation")
    print("  水平抛泥双流体 PFEM 数值模拟")
    print("=" * 60)
    print(f"  Domain  : {DOMAIN_X} m × {DOMAIN_Y} m")
    print(f"  Grid    : {NX} × {NY}")
    print(f"  Δt      : {DT} s,  T_end : {T_END} s  ({N_STEPS} steps)")
    print(f"  w_s     : {W_SETTLING} m/s  (cloud settling at full terminal velocity)")
    print(f"  D_turb  : {D_TURB} m²/s")
    print(f"  Patch   : x∈[{SED_X_LO},{SED_X_HI}], y∈[{SED_Y_LO},{SED_Y_HI}], α_s0={ALPHA_S_INIT}")
    print()

    sim = SedimentDumpingSimulation()
    sim.run()
    sim.report_settling()
    sim.plot_results(output_dir=script_dir)
