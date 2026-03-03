# 水平抛泥双流体 PFEM 数值模拟算例
# Sediment Dumping Two-Fluid PFEM Simulation Example

## 物理问题描述 (Physical Problem Description)

本算例模拟**水平抛泥过程**（sediment dumping）：一个初始的泥沙团在水中下沉、扩散的过程。
泥沙团初始位于水面处（上表面与自由水面齐平），在重力作用下下沉并逐渐扩散，同时引发自由水面的波动。

This example simulates a **sediment dumping** process: an initial sediment blob sinking and spreading in water.
The sediment blob starts at the water surface (upper edge flush with the free surface) and sinks under gravity while gradually spreading, inducing free-surface waves.

---

## 计算域 (Computational Domain)

| 参数 (Parameter) | 值 (Value) |
|---|---|
| 域尺寸 (Domain size) | 0.7 m × 0.7 m (2D) |
| 泥沙团 x 范围 (Sediment patch x) | [0.485, 0.515] m |
| 泥沙团 y 范围 (Sediment patch y) | [0.68, 0.70] m |
| 自由水面位置 (Free surface) | y = 0.70 m |

---

## 控制方程 (Governing Equations)

### 连续性方程 (Continuity Equation) — Eq. 4.51

$$\frac{\partial(\alpha_k \rho_k)}{\partial t} + \frac{\partial(\alpha_k \rho_k u_{kj})}{\partial x_j} = 0$$

### 动量方程 (Momentum Equation) — Eq. 4.52

$$\frac{\partial(\alpha_k \rho_k u_{ki})}{\partial t} + \frac{\partial(\alpha_k \rho_k u_{ki} u_{kj})}{\partial x_j} = -\alpha_k \frac{\partial p_f}{\partial x_i} + \frac{\partial(\alpha_k \rho_k \tilde{\tau}_{kij}^{*})}{\partial x_j} + (-1)^{\delta_{fk}} \gamma \alpha_s (u_{fi} - u_{si}) - (-1)^{\delta_{fk}} \gamma \frac{\varepsilon_s}{\alpha_f} \frac{\partial \alpha_s}{\partial x_i} + \alpha_k \rho_k g_i$$

- 水相 (water phase): k = f, δ_{fk} = 1
- 泥沙相 (sediment phase): k = s, δ_{fk} = 0

### 等效应力 (Equivalent Stress) — Eq. 4.50

$$\tilde{\tau}_{kij}^{*} = (\nu_k^0 + \nu_k^{SPS})\left(\frac{\partial \tilde{u}_{ki}}{\partial x_j} + \frac{\partial \tilde{u}_{kj}}{\partial x_i}\right)$$

### 相间拖曳力系数 (Inter-phase Drag Coefficient) — Eq. 4.37

$$\gamma = \lambda_d \frac{3 C_D \rho_f}{4 d_p} |\tilde{\mathbf{u}}_f - \tilde{\mathbf{u}}_s|$$

### Schiller-Naumann 拖曳力系数 — Eq. 4.38

$$C_D = \begin{cases} \frac{24}{Re_s}(1 + 0.15\, Re_s^{0.687}) & Re_s < 1000 \\ 0.44 & Re_s \geq 1000 \end{cases}$$

$$Re_s = \frac{|\tilde{\mathbf{u}}_f - \tilde{\mathbf{u}}_s|\, d_p}{\nu_f^0}$$

---

## 物理参数 (Physical Parameters)

| 参数 (Parameter) | 符号 (Symbol) | 值 (Value) |
|---|---|---|
| 泥沙粒径 (Particle diameter) | d_p | 8×10⁻⁴ m |
| 初始沉速 (Initial settling velocity) | w_s | 0.154 m/s |
| 泥沙团初始浓度 (Initial sediment concentration) | α_s | 0.606 |
| 水的密度 (Water density) | ρ_f | 1000 kg/m³ |
| 泥沙密度 (Sediment density) | ρ_s | 2650 kg/m³ |
| 混合密度 (Mixture density, sediment patch) | ρ_mix | ~2000 kg/m³ |
| 水的运动黏度 (Water kinematic viscosity) | ν_f | 1×10⁻⁶ m²/s |
| 重力加速度 (Gravitational acceleration) | g | 9.81 m/s² |

---

## 时间参数 (Time Parameters)

| 参数 (Parameter) | 值 (Value) |
|---|---|
| 时间步长 (Time step) | Δt = 0.001 s |
| 模拟总时间 (Total simulation time) | 0.5 s |

---

## 文件说明 (File Description)

| 文件 (File) | 说明 (Description) |
|---|---|
| `MainKratos.py` | 主运行脚本，包含自定义双流体分析类 (Main script with custom two-fluid analysis class) |
| `ProjectParameters.json` | Kratos PFEM 模拟参数配置 (Kratos PFEM simulation parameters) |
| `SedimentDumping.mdpa` | Kratos 模型部件文件，含网格和边界条件 (Kratos model part file with mesh and BCs) |
| `FluidMaterials.json` | 水和泥沙的物理属性 (Physical properties for water and sediment) |
| `generate_mesh.py` | 网格生成脚本 (Mesh generation script) |
| `README.md` | 本说明文件 (This documentation file) |

---

## 运行方法 (How to Run)

### 前提条件 (Prerequisites)

- KratosMultiphysics (≥ 9.0) with PfemFluidDynamicsApplication installed
- Python ≥ 3.8

### 步骤 (Steps)

1. **（可选）重新生成网格 (Optional: regenerate mesh)**

   ```bash
   cd examples/sediment_dumping_pfem
   python generate_mesh.py --nx 70 --ny 70
   ```

   The default resolution is `--nx 35 --ny 35`.  A higher resolution
   gives a better representation of the sediment patch but increases
   run time.

2. **运行模拟 (Run simulation)**

   From the repository root directory:

   ```bash
   python examples/sediment_dumping_pfem/MainKratos.py
   ```

   Or, from the example directory:

   ```bash
   cd examples/sediment_dumping_pfem
   python MainKratos.py
   ```

3. **查看输出 (View output)**

   VTK output files are written to `vtk_output/` every 10 time steps
   (every 0.01 s).  Open them with ParaView or VisIt.

   Fields available in the VTK output:
   - `VELOCITY` – velocity vector (m/s)
   - `PRESSURE` – fluid pressure (Pa)
   - `DENSITY`  – mixture density (kg/m³)

---

## 实现说明 (Implementation Notes)

PfemFluidDynamicsApplication 原生处理单流体或多流体界面问题
（通过不同属性的子模型部件实现）。本算例通过继承
`PfemFluidDynamicsAnalysis` 创建自定义分析类
`SedimentDumpingAnalysis`，在每个时间步中额外实现：

PfemFluidDynamicsApplication natively handles single-fluid or
multi-fluid interface problems (via sub-model-parts with different
properties).  This example extends `PfemFluidDynamicsAnalysis` with a
custom `SedimentDumpingAnalysis` class that performs per-step:

1. **初始化泥沙浓度 (Sediment concentration initialisation)**
   – POROSITY 变量存储泥沙体积分数 α_s

2. **混合物性质更新 (Mixture property update)**
   – 混合密度: ρ_mix = α_s ρ_s + α_f ρ_f
   – 混合黏度 (Einstein relation): μ_mix = μ_water (1 + 2.5 α_s)
   – Schiller-Naumann 拖曳力系数并更新体力 BODY_FORCE

3. **泥沙浓度输运 (Sediment concentration transport)**
   – 在拉格朗日框架下，浓度随节点运动自动平流
   – 同时施加简单扩散模型以体现浓度扩散过程

---

## 预期结果 (Expected Results)

1. **t = 0 s**: 泥沙团位于水面处，初始速度为零。
2. **t ≈ 0.08 s**: 泥沙云团开始下沉，自由水面产生明显波动。
3. **t ≈ 0.14 s**: 泥沙云团继续下沉并扩散，水面波动传播。
4. **t ≈ 0.5 s**: 泥沙逐渐扩散至整个计算域下半部。

密度场 (DENSITY) 的变化可直观反映泥沙浓度的空间分布。
速度场 (VELOCITY) 展示流场和泥沙云团的运动特征。

---

## 参考文献 (References)

- 双流体模型控制方程参见文档中的 Eq. 4.37–4.52
  (Two-fluid model governing equations: Eq. 4.37–4.52 in the reference document)
- KratosMultiphysics PfemFluidDynamicsApplication:
  `applications/PfemFluidDynamicsApplication/`
