# Water-Sediment Two-Phase MPM Simulation

A **Mixture-Theory Material Point Method (MPM)** simulation of a dense
sediment cloud settling through water.

## Physical model — why mixture theory?

### Previous (incorrect) approach
Earlier implementations used **separate material points for water and
sediment**.  This is physically wrong because the background computational
grid couples the two sets of MPs: the sinking sediment drags the water-MP
momentum field downward, so *all* water material points end up descending —
violating mass and momentum conservation between phases.

### Correct approach: every MP carries both phases simultaneously
In **mixture theory** each material point represents the *mixture* of water
and sediment at that location.  The sediment concentration field
C_s ∈ [0, 1] (sediment volume fraction) is a property of **every** MP:

| Field | Description |
|-------|-------------|
| **C_s** | Sediment volume fraction (0 = pure water, CS0 in blob) |
| **ρ_m = C_s ρ_s + (1−C_s) ρ_w** | Mixture density per MP |
| **u_m** | Single (mixture) velocity at each MP |
| **σ_m = −p I + 2μ_m ε̇** | Mixture Cauchy stress |

### Governing equations

```
Mixture momentum:          ρ_m Du_m/Dt = ∇·σ_m + ρ_m g
Concentration transport:   DC_s/Dt = −∇·[C_s(1−C_s) u_drift]
Constitutive (Newtonian):  σ_m = −p I + 2μ_m ε̇
Pressure (Jacobian form):  p = ρ_w g (L_y−y) − K ln(J),  J = V/V₀
Hindered settling drift:   u_drift = w_s (1−C_s)^4.65 ĵ  (downward)
```

The **Jacobian-based pressure** `p = ρ_w g(L−y) − K ln(J)` is:
- Self-correcting (returns to hydrostatic when J = 1) — no error accumulation
- Bounded (p → ∞ as J → 0, physically preventing volume collapse)
- Stable even with large bulk modulus K

### Why this gives the correct spreading behaviour

1. Dense blob (ρ_m > ρ_w) experiences net downward buoyancy force → **sinks**
2. As it sinks, water below is compressed (J < 1) → **excess pressure builds**
3. Horizontal pressure gradient at blob edges → **water pushed outward + upward**
4. Conservation: water circulates UP at the sides, DOWN ahead of blob
5. Settling drift redistributes C_s within the blob — dilute edges sink slower
6. Net result: **blob spreads laterally** (turbidity-current / gravity-current)

## Algorithm

Explicit Updated-Lagrangian MPM (MUSL variant):

```
For each time step:
  1. P2G  — transfer mass, momentum, internal forces to background grid
  2. Grid update  — solve nodal momentum equation
  3. BC   — no-penetration at walls
  4. G2P  — interpolate new velocities + strain rates back to MPs
  5. Constitutive update  — V, J, deviatoric stress τ
  6. Pressure update  — p = ρ_w g(L−y) − K ln(J)
  7. Advect MPs  — x += dt · v
  8. Concentration update — C_s += dt · w_eff(1−2C_s) ∂C_s/∂y
  9. Update ρ_m, mass
```

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| ρ_w | 1 000 kg/m³ | Water density |
| ρ_s | 2 650 kg/m³ | Sediment grain density |
| μ_w | 0.001 Pa·s | Water viscosity |
| K | 5×10⁴ Pa | Bulk modulus (α = ρ_w g L / K = 0.14, nearly incompressible) |
| d₅₀ | 0.3 mm | Grain diameter |
| w_s | 0.081 m/s | Stokes settling velocity |
| g_eff | 3.25 m/s² | Effective buoyancy acceleration |
| Domain | 0.70 × 0.70 m | |
| Grid | 35 × 35 | background cells |
| MPs | 4 900 | 2×2 per cell, fill entire domain |
| Δt | 0.001 s | CFL = 0.35 |
| T_end | 1.5 s | |

## Results

### Snapshots — sediment concentration C_s and velocity field

![Snapshots](results/sediment_snapshots.png)

Key observations:
- **t = 0 s** Initial state: sediment blob (red) at top, pure water (tan dots) everywhere.
- **t = 0.26 s** Toroidal circulation forms: blue (water) arrows at the sides point
  **upward** (+0.025 m/s), brown (sediment) arrows in the centre point **downward**.
  Water is NOT falling — it is being displaced upward by the sinking heavy cloud.
- **t = 0.52 s** Blob has descended and spread; crescent shape typical of a
  lock-release gravity current.
- **t = 0.68 s** Continued spreading; water at sides still shows upwelling
  (+0.089 m/s).

### Time-series diagnostics

![Trajectory](results/centroid_trajectory.png)

| Quantity | Value |
|----------|-------|
| Blob descent | 0.60 m → 0.43 m (gradual, **not** free fall) |
| Lateral spread | 13 cm → 17 cm (×1.26) |
| Max upwelling velocity | +0.197 m/s |
| Side-water vy > 0 | 55 % of time steps |

## Running the simulation

```bash
# Install dependencies (once)
pip install numpy matplotlib

# Run simulation (≈ 3–4 min, produces vtk_output/)
python run_simulation.py

# Generate result figures (results/)
python post_process.py
```

## File structure

```
sediment_dumping/
├── run_simulation.py         # standalone Python MPM (mixture theory)
├── post_process.py           # visualisation
├── results/
│   ├── sediment_snapshots.png
│   └── centroid_trajectory.png
├── vtk_output/               # raw VTK per time step (generated)
└── README.md
```

> **Note:** The `SedimentDumping_*.mdpa`, `ProjectParameters.json` and
> `SedimentDumping_materials.json` files remain in the folder as reference
> for the Kratos-based single-phase setup.  The active simulation is the
> pure-Python mixture MPM in `run_simulation.py`.

## References

- Bandara, S., & Soga, K. (2015). *Coupling of soil deformation and pore
  fluid flow using material point method.*
  Computers and Geotechnics, 63, 199–214.
- Sulsky, D., Zhou, S.J., & Schreyer, H.L. (1995). *Application of a
  particle-in-cell method to solid mechanics.*
  Computer Physics Communications, 87(1-2), 236–252.
