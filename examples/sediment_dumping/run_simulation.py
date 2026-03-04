#!/usr/bin/env python3
"""
Two-Phase Water-Sediment MPM Simulation using Mixture Theory
=============================================================

CORRECT physical model:
  Each material point simultaneously carries information for BOTH water
  and sediment phases via the sediment volume fraction C_s.

  Observed behaviours:
  • Water MPs at the LATERAL SIDES of the sinking blob move UPWARD (pushed
    by the horizontal pressure gradient from the descending heavy cloud).
    NOT all water descends — pure-water MPs above the blob move upward.
  • The sediment cloud spreads as it sinks (turbidity-current behaviour).
  • After impacting the bottom the blob spreads outward along the floor.

Governing equations (Mixture Theory)
--------------------------------------
  Mixture momentum:          ρ_m Du_m/Dt = ∇·σ_m + ρ_m g
  Concentration transport:   DC_s/Dt = −∇·[C_s(1−C_s) u_drift]
  Constitutive (Newtonian):  σ_m = −p I + 2μ_m ε̇
  Jacobian pressure:         p = ρ_w g (L_y−y) − K ln(J)   J = V/V₀
      (K >> ρ_w g L_y → nearly incompressible; ln-form is bounded and
       self-correcting: no pressure accumulation error)
  Hindered settling:         u_drift = w_s(1−C_s)^4.65 ĵ  (downward)

Algorithm: Explicit Updated-Lagrangian MPM (MUSL)
Reference: Bandara & Soga (2015), Computers and Geotechnics 63, 199–214.
"""

import os
import math
import struct
import warnings

import numpy as np

warnings.filterwarnings("ignore")   # suppress numpy overflow warnings in division

# ══════════════════════════════════════════════════════════════════
#  PHYSICAL PARAMETERS
# ══════════════════════════════════════════════════════════════════
RHO_W  = 1000.0          # water density               [kg/m³]
RHO_S  = 2650.0          # sediment grain density      [kg/m³]
MU_W   = 0.001           # water dynamic viscosity     [Pa·s]
# K_BULK >> ρ_w g L_y ≈ 6 900 Pa  →  nearly incompressible water column
# Using Jacobian (ln-form) pressure: p = p_ref − K ln(J), bounded & stable
K_BULK = 5.0e4           # artificial bulk modulus     [Pa]  (α = 0.14)
MIN_J  = 0.10            # minimum Jacobian (volume ratio) — prevents log singularity
GRAV   = 9.81            # gravitational acceleration  [m/s²]
D50    = 3e-4            # grain diameter 0.30 mm      [m]

# Stokes terminal settling velocity (single grain, clear water)
W_S = D50 ** 2 * (RHO_S - RHO_W) * GRAV / (18.0 * MU_W)

# ══════════════════════════════════════════════════════════════════
#  DOMAIN & GRID
# ══════════════════════════════════════════════════════════════════
LX = LY = 0.70           # domain size                 [m]

NCX = NCY = 35           # background-grid cells per direction
DX  = LX / NCX           # cell width   = 0.02 m
DY  = LY / NCY

NNX = NCX + 1
NNY = NCY + 1
NN  = NNX * NNY

# Material points: 2 × 2 per cell
MPPD = 2
NMP  = NCX * NCY * MPPD * MPPD    # = 4 900
V0   = DX * DY / MPPD ** 2        # initial MP volume  [m²]

# ══════════════════════════════════════════════════════════════════
#  SEDIMENT BLOB
# ══════════════════════════════════════════════════════════════════
BX0, BX1 = 0.28, 0.42   # blob x-extent [m]
BY0, BY1 = 0.56, 0.64   # blob y-extent [m]
CS0      = 0.30          # initial sediment volume fraction in blob

# ══════════════════════════════════════════════════════════════════
#  TIME INTEGRATION
# ══════════════════════════════════════════════════════════════════
DT    = 1.0e-3           # time step  [s]   CFL = dt·c/DX = 0.32 for K=5e4
T_END = 1.50
NSTEP = int(T_END / DT)
NOUT  = 75               # output snapshots

OUT_DIR = "vtk_output"


# ──────────────────────────────────────────────────────────────────
#  CONSTITUTIVE HELPERS
# ──────────────────────────────────────────────────────────────────
def mixture_density(Cs):
    return Cs * RHO_S + (1.0 - Cs) * RHO_W


def mixture_viscosity(Cs):
    """Einstein: μ_m = μ_w (1 + 2.5 C_s)."""
    return MU_W * (1.0 + 2.5 * Cs)


def hindered_ws(Cs):
    """Richardson-Zaki: w_eff = w_s (1 − C_s)^4.65."""
    return W_S * (1.0 - Cs) ** 4.65


def jacobian_pressure(yp, J):
    """
    Bounded, self-correcting pressure from Jacobian (volume ratio).

    p = ρ_w g (L_y − y)   [hydrostatic reference at current depth]
      − K ln(J)            [compressibility correction; K ln(J)→∞ as J→0]

    Initially J = 1 → p = hydrostatic → zero net force on water MPs ✓
    J < 1 (compression) → extra pressure that pushes material outward ✓
    J > 1 (expansion)   → reduced pressure ✓
    No incremental accumulation → numerically stable ✓
    """
    J_clamped = np.maximum(J, MIN_J)  # prevent log singularity (J ≥ MIN_J)
    return RHO_W * GRAV * (LY - yp) - K_BULK * np.log(J_clamped)


# ──────────────────────────────────────────────────────────────────
#  NODE / CELL HELPERS
# ──────────────────────────────────────────────────────────────────
def cell_coords(xp, yp):
    ci  = np.clip((xp / DX).astype(np.int32), 0, NCX - 1)
    cj  = np.clip((yp / DY).astype(np.int32), 0, NCY - 1)
    xi  = xp / DX - ci
    eta = yp / DY - cj
    return ci, cj, xi, eta


def shape_functions(xi, eta):
    """Bilinear Q4; order (0,0),(1,0),(1,1),(0,1); returns N, dN/dx, dN/dy each (NMP,4)."""
    N = np.stack([
        (1 - xi) * (1 - eta),
         xi       * (1 - eta),
         xi       *  eta,
        (1 - xi)  *  eta,
    ], axis=1)
    dN_dx = np.stack([
        -(1 - eta) / DX,  (1 - eta) / DX,  eta / DX, -eta / DX,
    ], axis=1)
    dN_dy = np.stack([
        -(1 - xi) / DY, -xi / DY,  xi / DY,  (1 - xi) / DY,
    ], axis=1)
    return N, dN_dx, dN_dy


def corner_node_ids(ci, cj):
    return np.stack([
        ci       + cj       * NNX,
        (ci + 1) + cj       * NNX,
        (ci + 1) + (cj + 1) * NNX,
        ci       + (cj + 1) * NNX,
    ], axis=1).astype(np.int32)


# ──────────────────────────────────────────────────────────────────
#  INITIALISE MATERIAL POINTS
# ──────────────────────────────────────────────────────────────────
def init_mps():
    sub  = np.arange(NMP)
    cell = sub // (MPPD * MPPD)
    ci   = cell % NCX
    cj   = cell // NCX
    si   = (sub % (MPPD * MPPD)) % MPPD
    sj   = (sub % (MPPD * MPPD)) // MPPD

    xp = (ci + (si + 0.5) / MPPD) * DX
    yp = (cj + (sj + 0.5) / MPPD) * DY

    Cs    = np.where(
        (xp >= BX0) & (xp <= BX1) & (yp >= BY0) & (yp <= BY1),
        CS0, 0.0)

    rho_m = mixture_density(Cs)
    Vp    = np.full(NMP, V0)
    J     = np.ones(NMP)          # Jacobian = Vp / V0 = 1 initially
    mp    = rho_m * Vp

    vx = np.zeros(NMP)
    vy = np.zeros(NMP)

    # Deviatoric stress (zero initially)
    sxx = np.zeros(NMP)
    syy = np.zeros(NMP)
    sxy = np.zeros(NMP)

    # Pressure from Jacobian (hydrostatic initially)
    press = jacobian_pressure(yp, J)

    return xp, yp, Cs, rho_m, Vp, J, mp, vx, vy, press, sxx, syy, sxy


# ──────────────────────────────────────────────────────────────────
#  ONE MPM STEP
# ──────────────────────────────────────────────────────────────────
def mpm_step(xp, yp, Cs, rho_m, Vp, J, mp, vx, vy, press, sxx, syy, sxy):

    ci, cj, xi, eta = cell_coords(xp, yp)
    N, dNdx, dNdy   = shape_functions(xi, eta)
    nids             = corner_node_ids(ci, cj)   # (NMP, 4)

    # total Cauchy stress:  σ = −p I + τ_dev
    s_xx = -press + sxx
    s_yy = -press + syy
    s_xy =           sxy

    # ── P2G ───────────────────────────────────────────────────────
    gm  = np.zeros(NN);  gmx = np.zeros(NN);  gmy = np.zeros(NN)
    gfx = np.zeros(NN);  gfy = np.zeros(NN)

    for k in range(4):
        nk = nids[:, k]
        wm = N[:, k] * mp
        gm  += np.bincount(nk, weights=wm,       minlength=NN)
        gmx += np.bincount(nk, weights=wm * vx,  minlength=NN)
        gmy += np.bincount(nk, weights=wm * vy,  minlength=NN)
        # internal force: f_int = −(σ ∇N) V  (from weak form of div σ)
        gfx += np.bincount(nk,
            weights=-(dNdx[:, k] * s_xx + dNdy[:, k] * s_xy) * Vp,
            minlength=NN)
        gfy += np.bincount(nk,
            weights=(-(dNdx[:, k] * s_xy + dNdy[:, k] * s_yy) * Vp
                     + N[:, k] * mp * (-GRAV)),
            minlength=NN)

    # ── grid solve ────────────────────────────────────────────────
    valid = gm > 1e-14
    gvx_new = np.where(valid, gmx / gm, 0.0) + DT * np.where(valid, gfx / gm, 0.0)
    gvy_new = np.where(valid, gmy / gm, 0.0) + DT * np.where(valid, gfy / gm, 0.0)

    # ── boundary conditions ───────────────────────────────────────
    gvy_new[np.arange(NNX)]               = np.maximum(0.0, gvy_new[np.arange(NNX)])
    gvx_new[np.arange(NNY) * NNX]         = np.maximum(0.0, gvx_new[np.arange(NNY) * NNX])
    gvx_new[NCX + np.arange(NNY) * NNX]  = np.minimum(0.0, gvx_new[NCX + np.arange(NNY) * NNX])

    # ── G2P ───────────────────────────────────────────────────────
    vx_new = np.zeros(NMP);  vy_new = np.zeros(NMP)
    eps_xx = np.zeros(NMP);  eps_yy = np.zeros(NMP);  eps_xy = np.zeros(NMP)

    for k in range(4):
        nk = nids[:, k]
        gvx_k = gvx_new[nk];  gvy_k = gvy_new[nk]
        vx_new += N[:, k] * gvx_k
        vy_new += N[:, k] * gvy_k
        eps_xx += dNdx[:, k] * gvx_k
        eps_yy += dNdy[:, k] * gvy_k
        eps_xy += 0.5 * (dNdy[:, k] * gvx_k + dNdx[:, k] * gvy_k)

    div_v = eps_xx + eps_yy

    # ── constitutive update ───────────────────────────────────────
    # Volume: V^{n+1} = V^n (1 + dt div v)  [clamp to prevent collapse]
    Vp_new = np.maximum(Vp * (1.0 + DT * div_v), MIN_J * V0)
    J_new  = Vp_new / V0   # Jacobian  (volume ratio w.r.t. initial)

    # Deviatoric Newtonian stress
    mu_m    = mixture_viscosity(Cs)
    dev_xx  = eps_xx - div_v / 3.0
    dev_yy  = eps_yy - div_v / 3.0
    sxx_new = 2.0 * mu_m * dev_xx
    syy_new = 2.0 * mu_m * dev_yy
    sxy_new = 2.0 * mu_m * eps_xy

    # ── advect MPs ────────────────────────────────────────────────
    xp_new = np.clip(xp + DT * vx_new, 0.0, LX)
    yp_new = np.clip(yp + DT * vy_new, 0.0, LY)

    # Pressure computed DIRECTLY from current J (no incremental accumulation)
    press_new = jacobian_pressure(yp_new, J_new)

    # ── two-phase: concentration update ──────────────────────────
    # Map C_s to nodes (mass-weighted)
    g_cs_n = np.zeros(NN);  g_cs_d = np.zeros(NN)
    for k in range(4):
        nk = nids[:, k];  wm = N[:, k] * mp
        g_cs_n += np.bincount(nk, weights=wm * Cs, minlength=NN)
        g_cs_d += np.bincount(nk, weights=wm,       minlength=NN)
    g_cs = np.where(g_cs_d > 1e-14, g_cs_n / g_cs_d, 0.0)

    # ∂C_s/∂y at nodes (central differences)
    g2  = g_cs.reshape(NNY, NNX)
    ddy = np.empty_like(g2)
    ddy[1:-1, :] = (g2[2:, :] - g2[:-2, :]) / (2.0 * DY)
    ddy[0,    :] = (g2[1, :]  - g2[0, :])   / DY
    ddy[-1,   :] = (g2[-1, :] - g2[-2, :])  / DY
    ddy = ddy.ravel()

    # ∂C_s/∂y interpolated to MPs
    dcs_dy_p = np.zeros(NMP)
    for k in range(4):
        dcs_dy_p += N[:, k] * ddy[nids[:, k]]

    # DC_s/Dt = w_eff (1−2C_s) ∂C_s/∂y  (from settling flux divergence)
    w_eff  = hindered_ws(Cs)
    Cs_new = np.clip(Cs + DT * w_eff * (1.0 - 2.0 * Cs) * dcs_dy_p, 0.0, 1.0)

    # ── update mixture density and mass ──────────────────────────
    rho_new = mixture_density(Cs_new)
    mp_new  = rho_new * Vp_new

    return (xp_new, yp_new, Cs_new, rho_new,
            Vp_new, J_new, mp_new, vx_new, vy_new,
            press_new, sxx_new, syy_new, sxy_new)


# ──────────────────────────────────────────────────────────────────
#  VTK OUTPUT
# ──────────────────────────────────────────────────────────────────
def write_vtk(step, xp, yp, Cs, rho_m, vx, vy, press):
    fname = os.path.join(OUT_DIR, f"MPM_Material_0_{step}.vtk")
    n     = NMP
    with open(fname, "wb") as f:
        def ws(s): f.write(s.encode())
        def wf(a): f.write(struct.pack(f">{len(a)}f", *np.nan_to_num(a, nan=0.0, posinf=1e6, neginf=-1e6).astype(np.float32)))
        def wi(a): f.write(struct.pack(f">{len(a)}i", *a.astype(np.int32)))

        ws("# vtk DataFile Version 4.0\nvtk output\nBINARY\n")
        ws("DATASET UNSTRUCTURED_GRID\n")
        ws(f"POINTS {n} float\n")
        pts = np.zeros(3 * n, dtype=np.float32)
        pts[0::3] = np.nan_to_num(xp); pts[1::3] = np.nan_to_num(yp)
        wf(pts);  ws("\n")
        ws(f"CELLS {n} {2 * n}\n")
        cells = np.empty(2 * n, dtype=np.int32)
        cells[0::2] = 1;  cells[1::2] = np.arange(n)
        wi(cells);  ws("\n")
        ws(f"CELL_TYPES {n}\n")
        wi(np.ones(n, dtype=np.int32));  ws("\n")
        ws(f"CELL_DATA {n}\n")
        ws("FIELD FieldData 4\n")
        ws(f"MP_CS 1 {n} float\n");       wf(Cs);     ws("\n")
        ws(f"MP_DENSITY 1 {n} float\n");  wf(rho_m);  ws("\n")
        ws(f"MP_VELOCITY 2 {n} float\n")
        vel = np.empty(2 * n, dtype=np.float32)
        vel[0::2] = vx;  vel[1::2] = vy
        wf(vel);  ws("\n")
        ws(f"MP_PRESSURE 1 {n} float\n"); wf(press);  ws("\n")


# ──────────────────────────────────────────────────────────────────
#  MAIN
# ──────────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    (xp, yp, Cs, rho_m, Vp, J, mp,
     vx, vy, press, sxx, syy, sxy) = init_mps()

    c_sound = math.sqrt(K_BULK / RHO_W)
    cfl     = DT * c_sound / DX
    alpha   = RHO_W * GRAV * LY / K_BULK
    rho_b   = float(mixture_density(CS0))
    g_eff   = GRAV * (rho_b - RHO_W) / rho_b

    print("=" * 65)
    print("  Two-Phase Mixture MPM  —  Water-Sediment Two-Phase Flow")
    print("  Each MP carries BOTH phases (water + sediment) via C_s")
    print("  Pressure: p = ρ_w g (L-y) − K ln(J)  [bounded, stable]")
    print("=" * 65)
    print(f"  MPs      : {NMP} ({(Cs>0).sum()} sediment-laden, {(Cs==0.0).sum()} pure-water)")
    print(f"  K_bulk   : {K_BULK:.0e} Pa  →  c = {c_sound:.1f} m/s, CFL = {cfl:.3f}, α = {alpha:.3f}")
    print(f"  w_s      : {W_S:.4f} m/s,  g_eff = {g_eff:.3f} m/s²")
    print(f"  Time     : {NSTEP} steps × {DT} s = {T_END} s  ({NOUT} snapshots)")
    print()
    print(f"  {'t':>7}  {'blob_cy':>8}  {'Δx_blob':>8}  {'vy_blob':>8}"
          f"  {'vy_side_water':>14}  {'Cs_max':>7}")
    print("  " + "─" * 62)

    out_every = max(1, NSTEP // NOUT)
    write_vtk(0, xp, yp, Cs, rho_m, vx, vy, press)

    for step in range(1, NSTEP + 1):
        (xp, yp, Cs, rho_m,
         Vp, J, mp, vx, vy,
         press, sxx, syy, sxy) = mpm_step(
            xp, yp, Cs, rho_m, Vp, J, mp, vx, vy, press, sxx, syy, sxy)

        if step % out_every == 0:
            write_vtk(step, xp, yp, Cs, rho_m, vx, vy, press)
            t = step * DT

            sed = Cs > 0.10
            if sed.any():
                cx = float(xp[sed].mean())
                cy = float(yp[sed].mean())
                lo = float(yp[sed].min())
                hi = float(yp[sed].max())
                xs = float(xp[sed].max() - xp[sed].min())
                vy_blob = float(vy[sed].mean())
            else:
                cx, cy, lo, hi, xs, vy_blob = LX/2, 0., 0., 0., 0., 0.

            # Water MPs directly BESIDE the blob (same height band, outside x)
            side = ((Cs < 0.01) &
                    (yp > lo - 0.03) & (yp < hi + 0.03) &
                    ((xp < cx - xs * 0.55) | (xp > cx + xs * 0.55)))
            vy_side = float(vy[side].mean()) if side.any() else 0.0

            Cs_max = float(Cs.max()) if np.isfinite(Cs).all() else float('nan')
            print(f"  {t:7.3f}  {cy:8.4f}  {xs:8.4f}  {vy_blob:+8.4f}"
                  f"  {vy_side:+14.4f}  {Cs_max:7.4f}")

    print()
    print("Simulation complete.")


if __name__ == "__main__":
    main()
