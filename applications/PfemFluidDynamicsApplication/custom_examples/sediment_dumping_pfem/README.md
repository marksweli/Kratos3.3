# Sediment Dumping PFEM Two-Fluid Simulation
## 水平抛泥双流体数值模拟（PFEM方法）

This directory contains a self-contained Python implementation of the
**sediment-dumping** (水平抛泥) benchmark using the **Particle Finite Element
Method (PFEM)** with a two-fluid (water–sediment) continuum model.

The code runs **without** a compiled Kratos installation.  It uses only
`numpy`, `scipy`, and `matplotlib`.

---

## Physical Model

### Two-fluid governing equations

**Stress tensor (Eq. 4.50)**

$$\widetilde{\tau}^*_{kij} = \left(\nu_k^0 + \nu_k^\text{SPS}\right)
  \left(\frac{\partial \widetilde{u}_{ki}}{\partial x_j}
       + \frac{\partial \widetilde{u}_{kj}}{\partial x_i}\right)$$

**Continuity (Eq. 4.51)**

$$\frac{\partial(\alpha_k \rho_k)}{\partial t}
  + \frac{\partial(\alpha_k \rho_k u_{kj})}{\partial x_j} = 0$$

**Momentum (Eq. 4.52)**

$$\frac{\partial(\alpha_k \rho_k u_{ki})}{\partial t}
  + \frac{\partial(\alpha_k \rho_k u_{ki} u_{kj})}{\partial x_j} =
  -\alpha_k \frac{\partial p_f}{\partial x_i}
  + \frac{\partial(\alpha_k \rho_k \tau^*_{kij})}{\partial x_j}
  + (-1)^{\delta_{fk}}\,\gamma\,\alpha_s\,(u_{fi}-u_{si})
  - (-1)^{\delta_{fk}}\,\gamma\,\frac{\varepsilon_s}{\alpha_f}
    \frac{\partial \alpha_s}{\partial x_i}
  + \alpha_k \rho_k g_i$$

where $\delta_{fk}=1$ for the water phase ($k=f$) and $\delta_{fk}=0$ for
the sediment phase ($k=s$).

**Drag coefficient (Eqs. 4.37–4.38, Schiller–Naumann 1935)**

$$\gamma = \lambda_d\,\frac{3 C_D \rho_f}{4 d_p}\,|\widetilde{\boldsymbol{u}}_f - \widetilde{\boldsymbol{u}}_s|$$

$$C_D = \begin{cases}
  \dfrac{24}{Re_s}(1 + 0.15\,Re_s^{0.687}) & Re_s < 1000 \\[6pt]
  0.44 & Re_s \ge 1000
\end{cases}$$

---

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Domain | 0.7 m × 0.7 m |
| Sediment block | 3 cm × 2 cm (top at free surface) |
| Sediment diameter $d_p$ | 0.8 mm |
| Initial settling velocity | 15.4 cm/s (downward) |
| Initial concentration $\alpha_s$ | 0.606 |
| Water density $\rho_f$ | 1000 kg/m³ |
| Sediment density $\rho_s$ | 2650 kg/m³ |
| Water viscosity $\nu_f$ | 1×10⁻⁶ m²/s |
| Gravity | 9.81 m/s² |
| Particle spacing | 1 cm |
| Time step | 0.001 s |

---

## PFEM Algorithm

Each time step applies the **fractional-step projection** method:

1. **Mesh reconstruction** – Delaunay triangulation of particle positions
   filtered by the Alpha-Shape criterion (circumradius ≤ 1/α).
2. **Drag coupling** – Schiller–Naumann formula gives inter-phase drag γ.
3. **Velocity prediction** – explicit integration of viscosity + drag +
   gravity (no pressure gradient).
4. **Pressure Poisson** – assembled and solved on the Delaunay FEM mesh
   (P1 triangular elements); Dirichlet $p=0$ at free-surface nodes.
5. **Velocity correction** – subtract pressure-gradient contribution.
6. **Volume-fraction update** – from the sediment continuity equation.
7. **Lagrangian advection** – move particles by the mixture velocity.
8. **Boundary conditions** – no-penetration at solid walls.

---

## File Structure

```
sediment_dumping_pfem/
├── run_sediment_dumping.py    # Main script – entry point
├── pfem_two_fluid_solver.py   # PFEM two-fluid solver class
├── mesh_generator.py          # Delaunay triangulation & Alpha-Shape
├── visualizer.py              # Matplotlib visualisation
├── requirements.txt           # Python dependencies
└── README.md                  # This file
```

---

## Quick Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run the simulation
python run_sediment_dumping.py
```

Output PNG files are written to the `output/` sub-directory.  Key snapshots
are produced at $t = 0.08\,\text{s}$ and $t = 0.14\,\text{s}$ to compare
with the reference results in Chapter 4 (Fig. 5.19).

---

## Expected Results

| Time | Description |
|------|-------------|
| $t = 0.08$ s | Free surface dips at the sediment location; cloud detaches and begins to spread |
| $t = 0.14$ s | Sediment cloud fully detached from free surface; lateral spread visible |

---

## References

* Two-fluid model equations: Chapter 4, Eqs. 4.50–4.52
* Schiller, L. & Naumann, A. (1935). *Über die grundlegenden Berechnungen bei
  der Schwerkraftaufbereitung.* Z. Ver. Dtsch. Ing., 77, 318–320.
* Idelsohn, S.R., Oñate, E. & Del Pin, F. (2004). *The particle finite element
  method: a powerful tool to solve incompressible flows with free surfaces and
  breaking waves.* Int. J. Numer. Meth. Engng., 61, 964–989.
