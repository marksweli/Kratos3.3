# Sediment Dumping Simulation — 2D MPM Example

This example uses the **KratosMultiphysics MPMApplication** to simulate a
two-dimensional **sediment dumping** (horizontal sediment release) process.
A sediment blob is released at the water surface and falls under gravity,
illustrating the dynamics of a settling sediment cloud.

---

## 1. Problem Description

### Physical Setup

| Parameter | Value |
|-----------|-------|
| Domain | 0 – 0.7 m × 0 – 0.7 m (2-D) |
| Sediment blob shape | Rectangle, 3 cm (x) × 2 cm (y) |
| Blob initial position | Top surface at water surface (y = 0.70 m); y ∈ [0.68, 0.70] m |
| Blob centre x | 0.35 m |
| Particle diameter | $d_p = 0.8$ mm |
| Initial settling velocity | $w_s = 15.4$ cm/s (downward) |
| Initial sediment concentration | $\alpha_s = 0.606$ |
| Water density | $\rho_f = 1000$ kg/m³ |
| Sediment grain density | $\rho_s = 2650$ kg/m³ |
| Mixture density (bulk) | $\rho_{mix} = \alpha_s\rho_s + (1-\alpha_s)\rho_f \approx 2000$ kg/m³ |
| Gravity | $g = 9.81$ m/s² (downward) |

### Governing Equations

The dual-fluid system is governed by the following conservation laws for
each phase $k$ (fluid $f$ and sediment $s$).

**Mass conservation (4.51):**

$$\frac{\partial(\alpha_k \rho_k)}{\partial t} + \frac{\partial(\alpha_k \rho_k u_{kj})}{\partial x_j} = 0$$

**Momentum conservation (4.52):**

$$\frac{\partial(\alpha_k \rho_k u_{ki})}{\partial t}
+ \frac{\partial(\alpha_k \rho_k u_{ki} u_{kj})}{\partial x_j}
= -\alpha_k \frac{\partial p_f}{\partial x_i}
+ \frac{\partial(\alpha_k \rho_k \tilde{\tau}^*_{kij})}{\partial x_j}
+ (-1)^{\delta_{fk}} \gamma \alpha_s (u_{fi} - u_{si})
- (-1)^{\delta_{fk}} \gamma \frac{\varepsilon_s}{\alpha_f} \frac{\partial \alpha_s}{\partial x_i}
+ \alpha_k \rho_k g_i$$

where $\delta_{fk} = 1$ for the fluid phase and $\delta_{fk} = 0$ for the
sediment phase.

**Effective stress tensor (4.50):**

$$\tilde{\tau}^*_{kij} = (\nu^0_k + \nu^{SPS}_k)
\left(\frac{\partial \tilde{u}_{ki}}{\partial x_j}
+ \frac{\partial \tilde{u}_{kj}}{\partial x_i}\right)$$

**Inter-phase drag coefficient (4.37):**

$$\gamma = \lambda_d \frac{3 C_D \rho_f}{4 d_p} |\tilde{u}_f - \tilde{u}_s|$$

**Drag coefficient — Schiller–Naumann (4.38):**

$$C_D = \begin{cases}
\dfrac{24}{Re_s}\!\left(1 + 0.15\,Re_s^{0.687}\right) & Re_s < 1000 \\[6pt]
0.44 & Re_s \geq 1000
\end{cases}, \qquad
Re_s = \frac{|\tilde{u}_f - \tilde{u}_s|\,d_p}{\nu^0_f}$$

In this MPM simulation the sediment blob is represented as an equivalent
single-phase continuum with the bulk mixture density
$\rho_{mix} \approx 2000$ kg/m³ and an isotropic linear-elastic
constitutive law; full two-fluid coupling is beyond the scope of the
MPMApplication but the kinematic behaviour (free fall, deformation,
interaction with boundaries) is captured.

---

## 2. Directory Structure

```
examples/sediment_dumping/
├── ProjectParameters.json          # Kratos project parameters
├── SedimentDumping_Body.mdpa       # Material-point body (sediment blob)
├── SedimentDumping_Grid.mdpa       # Background Eulerian grid (70×70)
├── SedimentDumping_materials.json  # Constitutive law and material data
├── run_simulation.py               # Main simulation runner
├── post_process.py                 # Post-processing & visualisation
├── configure_and_build.sh          # Dependency installation & CMake build
└── README.md                       # This document
```

---

## 3. Background Grid

| Property | Value |
|----------|-------|
| Type | Structured quad mesh, `Element2D4N` |
| Domain | 0 – 0.7 m × 0 – 0.7 m |
| Element size | Δx = Δy = 0.01 m |
| Grid size | 70 × 70 elements, 71 × 71 nodes |
| Bottom BC | Fixed (zero displacement in x and y) |
| Left / Right BC | Free-slip (zero x-displacement) |
| Top BC | Free surface (no constraint) |

---

## 4. Compilation

### Prerequisites

- Ubuntu 20.04 or newer (or equivalent Debian-based distribution)
- CMake ≥ 3.16
- GCC ≥ 9 / G++ ≥ 9
- Python 3.8+
- Boost libraries

### Build

```bash
# From the repository root:
bash examples/sediment_dumping/configure_and_build.sh
```

The script will:
1. Install system dependencies with `apt-get`.
2. Set compiler and path environment variables.
3. Configure CMake with `MPMApplication` and `LinearSolversApplication`.
4. Build and install Kratos in `Release` mode using all available cores.

---

## 5. Running the Simulation

```bash
cd examples/sediment_dumping
python3 run_simulation.py
```

VTK output files are written to the `vtk_output/` sub-folder at every
0.01 s of simulation time.  The simulation covers 0 – 0.2 s with a time
step of Δt = 0.001 s.

---

## 6. Post-Processing

```bash
cd examples/sediment_dumping
python3 post_process.py
```

Two PNG images are generated in the `results/` directory:

| File | Content |
|------|---------|
| `results/sediment_snapshots.png` | Sediment material-point positions at t = 0.08 s and t = 0.14 s |
| `results/centroid_trajectory.png` | Vertical descent of the sediment cloud centroid over time |

Example snapshot images (generated after running the simulation):

![Sediment snapshots](results/sediment_snapshots.png)
![Centroid trajectory](results/centroid_trajectory.png)

---

## 7. Expected Physical Phenomena

| Time | Observed behaviour |
|------|--------------------|
| t = 0 s | Sediment blob sits at the water surface; initial downward velocity 0.154 m/s |
| t ≈ 0.08 s | Blob has descended ~1–2 cm; free surface is dragged downward by the falling sediment cloud |
| t ≈ 0.14 s | Blob has detached from the surface and continues to fall freely |
| t = 0.20 s | Blob approaches mid-domain depth; lateral spreading is visible |

The sediment cloud:
- Falls under gravity while decelerating due to fluid drag.
- Spreads laterally as it descends.
- Initially drags the free surface downward before detaching.

---

## 8. Key Files Reference

| File | Purpose |
|------|---------|
| `applications/MPMApplication/python_scripts/mpm_analysis.py` | Main MPM analysis class |
| `applications/MPMApplication/python_scripts/assign_gravity_to_material_point_process.py` | Gravity assignment |
| `applications/MPMApplication/python_scripts/assign_initial_velocity_to_material_point_process.py` | Initial velocity |
| `applications/MPMApplication/python_scripts/mpm_vtk_output_process.py` | VTK output |
| `applications/MPMApplication/tests/gravity_tests/` | Reference test case |
| `INSTALL.md` | Full Kratos build instructions |
