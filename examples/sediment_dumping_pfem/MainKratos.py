"""
MainKratos.py
~~~~~~~~~~~~~
主运行脚本 – 水平抛泥双流体 PFEM 数值模拟
Main script – sediment dumping two-fluid PFEM numerical simulation

使用方法 (Usage)
---------------
  cd <Kratos build/install dir>
  python examples/sediment_dumping_pfem/MainKratos.py

物理模型 (Physical model)
--------------------------
基于 KratosMultiphysics PfemFluidDynamicsApplication，实现双流体
（水-泥沙）模型，控制方程参见 README.md。
Based on KratosMultiphysics PfemFluidDynamicsApplication; implements a
two-fluid (water-sediment) model.  Governing equations are described
in README.md.
"""

import os
import sys
import math

import KratosMultiphysics
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis import (
    PfemFluidDynamicsAnalysis,
)

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
K_SPREADING      = 0.5         # 浓度扩散速率系数 (1/s) – concentration spreading rate coefficient

# Sediment patch bounds (m) – 泥沙团位置
SED_X_LO, SED_X_HI = 0.485, 0.515
SED_Y_LO, SED_Y_HI = 0.68,  0.70


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


# ===========================================================================
# Custom analysis stage that adds two-fluid physics each step
# 自定义分析阶段，在每个时间步中添加双流体物理
# ===========================================================================

class SedimentDumpingAnalysis(PfemFluidDynamicsAnalysis):
    """
    Extends PfemFluidDynamicsAnalysis with:
      1. Initialisation of sediment concentration (POROSITY) on all nodes.
      2. Per-step update of:
         - inter-phase drag force added to body force (BODY_FORCE)
         - mixture density / viscosity (DENSITY, VISCOSITY)
         - sediment-concentration transport (explicit advection via DISPLACEMENT)

    此类继承自 PfemFluidDynamicsAnalysis，在每个时间步：
      1. 初始化泥沙浓度字段（POROSITY）
      2. 更新相间拖曳力（加入体力 BODY_FORCE）
      3. 更新混合密度/黏度（DENSITY, VISCOSITY）
      4. 泥沙浓度的显式平流输运
    """

    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        self._dt = parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()

    # ------------------------------------------------------------------
    def Initialize(self):
        super().Initialize()
        self._initialize_sediment_concentration()

    # ------------------------------------------------------------------
    def _initialize_sediment_concentration(self):
        """
        Set initial sediment volume fraction (POROSITY) on every node.
        Nodes inside the sediment patch get ALPHA_S_INIT; all others get 0.
        为每个节点设置初始泥沙体积分数（POROSITY）。
        位于泥沙团内的节点赋值 ALPHA_S_INIT，其余节点赋 0。
        """
        model_part = self.main_model_part
        for node in model_part.Nodes:
            x, y = node.X, node.Y
            if (SED_X_LO <= x <= SED_X_HI) and (SED_Y_LO <= y <= SED_Y_HI):
                alpha_s = ALPHA_S_INIT
            else:
                alpha_s = 0.0
            node.SetSolutionStepValue(KratosMultiphysics.POROSITY, 0, alpha_s)
        KratosMultiphysics.Logger.Print(
            "::[SedimentDumping]:: Sediment concentration initialised.",
            label=""
        )

    # ------------------------------------------------------------------
    def InitializeSolutionStep(self):
        """
        Before each solve: update mixture properties and drag body force.
        每次求解前：更新混合物性质和拖曳体力。
        """
        super().InitializeSolutionStep()
        self._update_mixture_properties()

    # ------------------------------------------------------------------
    def FinalizeSolutionStep(self):
        """
        After each solve: advect sediment concentration.
        每次求解后：平流输运泥沙浓度。
        """
        super().FinalizeSolutionStep()
        self._advect_sediment_concentration()

    # ------------------------------------------------------------------
    def _update_mixture_properties(self):
        """
        Update effective density, viscosity and add inter-phase drag body force.
        更新有效密度、黏度，并添加相间拖曳体力。

        The two-fluid momentum equations (Eq. 4.52) introduce an
        inter-phase drag acceleration γ * Δu that acts as a body force
        on each phase.  In single-mesh PFEM we represent this as an
        effective body force on the mixture.

        双流体动量方程 (Eq. 4.52) 引入了相间拖曳加速度 γ*Δu，
        在单网格 PFEM 中将其表示为混合物的等效体力。
        """
        model_part = self.main_model_part
        gravity_y  = -9.81   # m/s²

        for node in model_part.Nodes:
            alpha_s = node.GetSolutionStepValue(KratosMultiphysics.POROSITY, 0)
            alpha_s = max(0.0, min(1.0, alpha_s))
            alpha_f = 1.0 - alpha_s

            # ---- Mixture density (线性混合密度) ----
            rho_mix = alpha_s * RHO_SEDIMENT + alpha_f * RHO_WATER

            # ---- Mixture dynamic viscosity (简单混合黏度) ----
            # Einstein relation: μ_mix = μ_water * (1 + 2.5*α_s)
            mu_water = RHO_WATER * NU_WATER          # 0.001 Pa·s
            mu_mix   = mu_water * (1.0 + 2.5 * alpha_s)

            node.SetSolutionStepValue(KratosMultiphysics.DENSITY,   0, rho_mix)
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY,  0, mu_mix)

            # ---- Inter-phase drag body force ----
            # Velocity difference (water velocity = mixture velocity in this approach)
            vx = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0)
            vy = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0)

            # Sediment settling velocity (w_s) – terminal velocity model
            # w_s = W_SETTLING m/s (given); applied as relative velocity in y-direction
            ws = W_SETTLING

            # Relative velocity: fluid relative to sediment
            dv_x = 0.0
            dv_y = ws * alpha_s      # sediment falls at ws relative to fluid

            dv_mag = math.sqrt(dv_x ** 2 + dv_y ** 2)

            # Sediment Reynolds number (Eq. 4.38 context)
            re_s = dv_mag * D_PARTICLE / NU_WATER if dv_mag > 0 else 0.0
            cd   = schiller_naumann_cd(re_s)

            # Drag coefficient γ (Eq. 4.37)
            gamma = LAMBDA_D * (3.0 * cd * RHO_WATER) / (4.0 * D_PARTICLE) * dv_mag if dv_mag > 0 else 0.0

            # Drag body force on mixture (per unit mass)
            # f_drag = γ * Δu / ρ_mix
            if rho_mix > 0:
                bfx = gamma * dv_x / rho_mix
                bfy = gamma * dv_y / rho_mix + gravity_y
            else:
                bfx = 0.0
                bfy = gravity_y

            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, 0, bfx)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, 0, bfy)

    # ------------------------------------------------------------------
    def _advect_sediment_concentration(self):
        """
        Explicit first-order Lagrangian advection of sediment concentration.
        泥沙浓度的显式一阶拉格朗日平流输运。

        In PFEM the nodes move with the velocity field (Lagrangian frame),
        so the concentration scalar is simply advected by the node
        displacement that was computed this step.  We also apply a simple
        diffusion limiter to keep 0 ≤ α_s ≤ 1.

        在 PFEM 拉格朗日框架中，节点随速度场运动，
        因此浓度标量随节点位移自动平流。
        同时施加简单扩散限制保持 0 ≤ α_s ≤ 1。
        """
        model_part = self.main_model_part
        dt         = self._dt

        for node in model_part.Nodes:
            alpha_s_old = node.GetSolutionStepValue(KratosMultiphysics.POROSITY, 0)

            # Simple decay/spreading: d(alpha_s)/dt = -k * alpha_s * (1 - alpha_s)
            # where k is related to the settling-induced spreading
            d_alpha = -K_SPREADING * alpha_s_old * (1.0 - alpha_s_old) * dt

            alpha_s_new = max(0.0, min(1.0, alpha_s_old + d_alpha))
            node.SetSolutionStepValue(KratosMultiphysics.POROSITY, 0, alpha_s_new)


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == "__main__":
    # Path to this file's directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Load project parameters
    param_file = os.path.join(script_dir, "ProjectParameters.json")
    with open(param_file, "r") as f:
        parameters = KratosMultiphysics.Parameters(f.read())

    # Update paths to be relative to repo root so the solver can find them
    repo_root = os.path.dirname(os.path.dirname(script_dir))
    os.chdir(repo_root)

    model      = KratosMultiphysics.Model()
    simulation = SedimentDumpingAnalysis(model, parameters)
    simulation.Run()
