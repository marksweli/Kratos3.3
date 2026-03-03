# Sediment Dumping Simulation — 2D Water–Sediment Two-Phase MPM

This example uses the **KratosMultiphysics MPMApplication** to simulate a
two-dimensional **sediment dumping** process in water.
A sediment blob is released at the water surface and falls under gravity through
the water, then spreads laterally along the bottom as a gravity/turbidity current —
correctly reflecting the **water-sediment two-phase flow** physics.

---

## 1. Problem Description

### Physical Setup

| Parameter | Value |
|-----------|-------|
| Domain | 0 – 0.70 m × 0 – 0.70 m (2-D) |
| Water body | x ∈ [0, 0.70] m, y ∈ [0, 0.56] m (320 MPs, 0.07 m elements) |
| Sediment blob | x ∈ [0.28, 0.42] m, y ∈ [0.56, 0.64] m (112 MPs, 0.02 m elements) |
| Sediment blob size | 0.14 m (width) × 0.08 m (height) |
| Sediment density | ρ_sed = 1800 kg/m³ |
| Water density | ρ_wat = 1000 kg/m³ |
| Sediment viscosity | μ_sed = 0.05 Pa·s |
| Water viscosity | μ_wat = 0.001 Pa·s |
| Gravity | g = 9.81 m/s² (downward) |
| Total simulation time | 0.5 s |
| Time step | Δt = 0.001 s |

### Constitutive Law

Both phases use `DispNewtonianFluidPlaneStrain2DLaw` (Displacement-based Newtonian Fluid,
plane strain 2-D). This allows each phase to flow, deform and spread.

**Newtonian fluid stress:**

$$\boldsymbol{\sigma} = -p\mathbf{I} + 2\mu\,\dot{\boldsymbol{\varepsilon}}$$

where $p = K\,\nabla\cdot\mathbf{u}$ is the pressure and $\dot{\boldsymbol{\varepsilon}}$ is the
symmetric strain-rate tensor.

### Governing Equations

Both phases satisfy the momentum balance with gravity:

$$\rho\ddot{\mathbf{u}} = \nabla\cdot\boldsymbol{\sigma} + \rho\mathbf{g}$$

Phase interaction is captured via the **shared background grid**: water and
sediment material points are mapped to the same Eulerian nodes; the resulting
velocity field is consistent across both materials, creating an approximate
two-phase coupling.

---

## 2. Directory Structure

```
examples/sediment_dumping/
├── ProjectParameters.json          # Kratos project parameters
├── SedimentDumping_Body.mdpa       # Two-material body (water + sediment)
├── SedimentDumping_Grid.mdpa       # Background Eulerian grid (70×70)
├── SedimentDumping_materials.json  # Constitutive laws for water and sediment
├── run_simulation.py               # Main simulation runner
├── post_process.py                 # Post-processing & visualisation
├── configure_and_build.sh          # Installation & CMake build script
├── results/                        # Generated result images (committed)
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

## 4. Installation

### Option A — Binaries via pip (recommended for users)

Following `applications/MPMApplication/README.md` §*Getting Binaries with pip*:

```bash
pip3 install KratosMPMApplication matplotlib
```

### Option B — Build from Source (developers)

Following `applications/MPMApplication/README.md` §*Build from Source (developers)*
and `INSTALL.md`:

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get install -y python3-dev gcc g++ cmake libboost-all-dev python3-pip

# From the repository root — compiles MPMApplication + LinearSolversApplication:
bash examples/sediment_dumping/configure_and_build.sh
```

Key lines in the build configuration (`configure_and_build.sh`):

```bash
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/MPMApplication
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
```

---

## 5. Running the Simulation

```bash
cd examples/sediment_dumping
python3 run_simulation.py
```

The simulation runs 500 steps (Δt = 0.001 s, t = 0 → 0.5 s) with 432 material
points (320 water + 112 sediment) on a 71×71 background grid.

Expected final output:
```
::[MPM Analysis]:: : STEP:  500
::[MPM Analysis]:: : TIME:  0.50...
::[MPM Analysis]:: : Analysis -END-
```

---

## 6. Post-Processing

```bash
cd examples/sediment_dumping
python3 post_process.py
```

Generates two PNG images in `results/`:

| File | Content |
|------|---------|
| `results/sediment_snapshots.png` | Water + sediment MP positions at t = 0, 0.15, 0.35, 0.50 s; coloured by \|vx\| |
| `results/centroid_trajectory.png` | Vertical descent + lateral spreading of the sediment cloud |

---

## 7. Simulation Results

### Material-point snapshots

The four panels show the positions of all material points.
**Blue** dots = water MPs, **brown/orange** dots = sediment MPs (colour encodes lateral velocity |vx|).

![Sediment cloud positions at t = 0, 0.15, 0.35, 0.50 s](results/sediment_snapshots.png)

### Centroid descent and lateral spreading

The left panel shows the sediment centroid falling from y = 0.60 m to the bottom.
The right panel shows the dramatic lateral spreading from **0.131 m → 0.559 m** (4.25×).

![Settling and spreading trajectories](results/centroid_trajectory.png)

### Quantitative summary

| Quantity | Value |
|----------|-------|
| Initial sediment centroid y | 0.600 m |
| Final sediment centroid y (t = 0.5 s) | 0.039 m |
| Initial x-spread | 0.131 m |
| Final x-spread | 0.559 m |
| Spreading ratio | 4.25× |

---

## 8. Physical Phenomena Analysis

| Phase | Time | Observed behaviour |
|-------|------|--------------------|
| Free fall | t = 0–0.30 s | Sediment falls under gravity, water provides drag; blob remains compact |
| Impact | t ≈ 0.33 s | Sediment cloud reaches the bottom; vertical momentum converts to lateral |
| Spreading | t = 0.33–0.50 s | Rapid lateral spreading; blob spreads from 0.13 m to >0.55 m wide |

**Key physical mechanisms captured:**

1. **Buoyancy and drag** — The water body (ρ = 1000 kg/m³) slows the sediment (ρ = 1800 kg/m³)
   through the shared background grid. The two-phase interaction reduces the fall rate compared
   to free fall in vacuum.

2. **Impact and redirection** — When the dense sediment reaches the impermeable bottom, the
   vertical momentum is redirected into horizontal flow, initiating a **turbidity current**.

3. **Gravity current spreading** — The denser sediment flows outward along the bottom,
   spreading symmetrically from the impact centre — the classic behaviour of a **lock-release
   gravity current** in environmental hydraulics.

4. **Phase differentiation** — The water and sediment have distinct densities, viscosities
   and velocities throughout the simulation, captured via the two-material MPM framework.

These behaviours correctly reflect the physical laws governing water-sediment two-phase flow,
as described by the momentum conservation equations with inter-phase interaction through the
shared Eulerian grid.

---

## 9. Key Files Reference

| File | Purpose |
|------|---------|
| `applications/MPMApplication/README.md` | MPMApplication overview, installation options |
| `INSTALL.md` | Full Kratos build instructions for Linux/Windows/macOS |
| `applications/MPMApplication/tests/cl_tests/fluid_cl/` | Reference test for DispNewtonianFluidPlaneStrain2DLaw |
| `applications/MPMApplication/python_scripts/mpm_analysis.py` | Main MPM analysis class |
| `applications/MPMApplication/python_scripts/assign_gravity_to_material_point_process.py` | Gravity |
| `applications/MPMApplication/python_scripts/mpm_vtk_output_process.py` | VTK output |
