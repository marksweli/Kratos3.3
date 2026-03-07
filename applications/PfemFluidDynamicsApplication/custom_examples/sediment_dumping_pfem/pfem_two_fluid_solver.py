"""
pfem_two_fluid_solver.py
========================
PFEM Two-Fluid Solver for the sediment-dumping (水平抛泥) simulation.

Physical model
--------------
Two-fluid continuum equations for a water–sediment mixture
(water phase k=f, sediment phase k=s):

Stress tensor (Eq. 4.50)
    τ*_kij = (ν_k^0 + ν_k^SPS) * (∂ũ_ki/∂x_j + ∂ũ_kj/∂x_i)

Continuity (Eq. 4.51)
    ∂(α_k ρ_k)/∂t + ∂(α_k ρ_k u_kj)/∂x_j = 0

Momentum (Eq. 4.52)
    ∂(α_k ρ_k u_ki)/∂t + ∂(α_k ρ_k u_ki u_kj)/∂x_j =
        −α_k ∂p_f/∂x_i
        + ∂(α_k ρ_k τ*_kij)/∂x_j
        + (−1)^{δ_fk} γ α_s (u_fi − u_si)
        − (−1)^{δ_fk} γ ε_s/α_f ∂α_s/∂x_i
        + α_k ρ_k g_i

    where δ_fk = 1 for k=f (water), δ_fk = 0 for k=s (sediment).

Drag coefficient (Eqs. 4.37–4.38, Schiller–Naumann 1935)
    γ = λ_d * (3 C_D ρ_f)/(4 d_p) * |ũ_f − ũ_s|

    C_D = 24/Re_s * (1 + 0.15 Re_s^{0.687})   if Re_s < 1000
    C_D = 0.44                                  if Re_s ≥ 1000
    Re_s = |ũ_f − ũ_s| d_p / ν_f

Numerical method
----------------
Fractional-step (projection) PFEM on P1 triangular elements:
  1. Delaunay triangulation + Alpha-Shape → valid FEM mesh
  2. Predict velocities (explicit advection + viscosity + drag + gravity)
  3. Solve pressure Poisson (FEM, P1 elements)
  4. Correct velocities
  5. Update volume fractions (continuity)
  6. Advect particles (Lagrangian step)
  7. Apply wall boundary conditions
"""

import warnings
import numpy as np
from scipy.sparse import coo_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

from mesh_generator import (
    delaunay_triangulate,
    alpha_shape_filter,
    build_element_data,
    find_free_surface_nodes,
)


class PFEMTwoFluidSolver:
    """
    PFEM solver for two-fluid (water–sediment) flow.

    Parameters
    ----------
    params : dict
        Physical and numerical parameters (see run_sediment_dumping.py).
    """

    def __init__(self, params):
        self.params = params

        # Physical constants
        self.rho_f = params['rho_f']
        self.rho_s = params['rho_s']
        self.nu_f  = params['nu_f']
        self.nu_s  = params['nu_s']
        self.d_p   = params['d_p']
        self.g     = params['g']
        self.Lx    = params['Lx']
        self.Ly    = params['Ly']
        self.alpha_param = params['alpha_shape']

        # Particle field arrays (set in setup_particles)
        self.pos      = None   # (N, 2)  positions
        self.uf       = None   # (N, 2)  water velocity
        self.us       = None   # (N, 2)  sediment velocity
        self.alpha_s  = None   # (N,)    sediment volume fraction
        self.pressure = None   # (N,)    dynamic pressure

        self.time = 0.0
        self.n    = 0

        # Cached element data (rebuilt each step)
        self._idx  = None
        self._area = None
        self._b    = None
        self._c    = None

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def setup_particles(self, positions, alpha_s_init,
                        uf_init=None, us_init=None):
        """
        Assign initial particle fields.

        Parameters
        ----------
        positions    : ndarray (N, 2)
        alpha_s_init : ndarray (N,)   – initial sediment volume fraction
        uf_init      : ndarray (N, 2) – initial water velocity (default 0)
        us_init      : ndarray (N, 2) – initial sediment velocity (default 0)
        """
        self.pos     = np.asarray(positions, dtype=float)
        self.n       = len(self.pos)
        self.alpha_s = np.clip(np.asarray(alpha_s_init, dtype=float),
                               0.0, 0.85)

        self.uf = (np.zeros((self.n, 2))
                   if uf_init is None else np.asarray(uf_init, dtype=float))
        self.us = (np.zeros((self.n, 2))
                   if us_init is None else np.asarray(us_init, dtype=float))

        self.pressure = np.zeros(self.n)

        print(f"[Solver] {self.n} particles  "
              f"max αs = {np.max(self.alpha_s):.3f}  "
              f"mean αs = {np.mean(self.alpha_s):.4f}")

    # ------------------------------------------------------------------
    # Drag force  (Schiller–Naumann, Eqs. 4.37–4.38)
    # ------------------------------------------------------------------

    def _drag_coefficient_cd(self, u_rel_mag):
        """
        Schiller–Naumann (1935) drag coefficient C_D (Eq. 4.38).

        C_D = 24/Re_s * (1 + 0.15 Re_s^{0.687})   Re_s < 1000
        C_D = 0.44                                  Re_s ≥ 1000
        """
        Re_s = np.maximum(u_rel_mag * self.d_p / self.nu_f, 1.0e-10)
        CD = np.where(
            Re_s < 1000.0,
            24.0 / Re_s * (1.0 + 0.15 * Re_s ** 0.687),
            0.44,
        )
        return CD

    def _gamma(self, u_rel, alpha_s):
        """
        Inter-phase drag coefficient γ (Eq. 4.37).

            γ = λ_d * (3 C_D ρ_f)/(4 d_p) * |u_f − u_s|

        λ_d = α_f^{−n} is the Richardson–Zaki hindered-settling correction
        with n = 3.7.
        """
        u_rel_mag = np.linalg.norm(u_rel, axis=1)   # (N,)
        CD = self._drag_coefficient_cd(u_rel_mag)

        alpha_f   = np.maximum(1.0 - alpha_s, 0.01)
        n_rz      = 3.7
        lambda_d  = alpha_f ** (-n_rz)               # hindered-settling factor

        gamma = (lambda_d * 3.0 * CD * self.rho_f
                 / (4.0 * self.d_p) * u_rel_mag)
        return gamma

    # ------------------------------------------------------------------
    # FEM mesh reconstruction
    # ------------------------------------------------------------------

    def _rebuild_mesh(self):
        """
        Perform Delaunay triangulation → Alpha-Shape filter →
        build vectorised element data.
        """
        if self.n < 3:
            self._idx  = np.empty((0, 3), dtype=int)
            self._area = np.empty(0)
            self._b    = np.empty((0, 3))
            self._c    = np.empty((0, 3))
            return False

        _, simplices = delaunay_triangulate(self.pos)
        simplices    = alpha_shape_filter(self.pos, simplices,
                                          self.alpha_param)

        if len(simplices) == 0:
            self._idx  = np.empty((0, 3), dtype=int)
            self._area = np.empty(0)
            self._b    = np.empty((0, 3))
            self._c    = np.empty((0, 3))
            return False

        self._idx, self._area, self._b, self._c = build_element_data(
            self.pos, simplices
        )
        return len(self._idx) > 0

    # ------------------------------------------------------------------
    # Vectorised FEM operators
    # ------------------------------------------------------------------

    def _assemble_laplacian(self):
        """
        Global pressure-Laplacian stiffness matrix (P1, COO → CSR).

            K_{mn} = Σ_e A_e (b_m b_n + c_m c_n)
        """
        ne = len(self._idx)
        if ne == 0:
            from scipy.sparse import csr_matrix
            return csr_matrix((self.n, self.n))

        A = self._area            # (ne,)
        b = self._b               # (ne, 3)
        c = self._c               # (ne, 3)

        # Local 3×3 stiffness tensors: (ne, 3, 3)
        K_e = A[:, None, None] * (
            b[:, :, None] * b[:, None, :]
            + c[:, :, None] * c[:, None, :]
        )

        # COO assembly
        rows = np.repeat(self._idx, 3, axis=1).ravel()   # ne*9
        cols = np.tile(self._idx, (1, 3)).ravel()         # ne*9
        data = K_e.ravel()

        K = coo_matrix((data, (rows, cols)),
                       shape=(self.n, self.n))
        return K.tocsr()

    def _lumped_mass(self):
        """
        Lumped mass vector  M_i = Σ_{e∋i} A_e / 3.
        """
        M = np.zeros(self.n)
        for k in range(3):
            np.add.at(M, self._idx[:, k], self._area / 3.0)
        return M

    def _divergence(self, u):
        """
        Node-centred divergence via area-weighted average of
        element-constant values  ∇·u|_e = Σ_i (b_i u_xi + c_i u_yi).

        Parameters
        ----------
        u : ndarray, shape (N, 2)

        Returns
        -------
        div : ndarray, shape (N,)
        """
        ux = u[self._idx, 0]   # (ne, 3)
        uy = u[self._idx, 1]

        div_e = (np.einsum('ei,ei->e', self._b, ux)
                 + np.einsum('ei,ei->e', self._c, uy))   # (ne,)

        div     = np.zeros(self.n)
        area_w  = np.zeros(self.n)
        for k in range(3):
            np.add.at(div,    self._idx[:, k], self._area * div_e)
            np.add.at(area_w, self._idx[:, k], self._area)

        mask = area_w > 1.0e-20
        div[mask] /= area_w[mask]
        return div

    def _gradient(self, f):
        """
        Node-centred gradient of a scalar field via area-weighted average.

        Parameters
        ----------
        f : ndarray, shape (N,)

        Returns
        -------
        grad : ndarray, shape (N, 2)
        """
        f_e = f[self._idx]   # (ne, 3)

        gx_e = np.einsum('ei,ei->e', self._b, f_e)
        gy_e = np.einsum('ei,ei->e', self._c, f_e)

        grad   = np.zeros((self.n, 2))
        area_w = np.zeros(self.n)
        for k in range(3):
            np.add.at(grad[:, 0], self._idx[:, k], self._area * gx_e)
            np.add.at(grad[:, 1], self._idx[:, k], self._area * gy_e)
            np.add.at(area_w,     self._idx[:, k], self._area)

        mask = area_w > 1.0e-20
        grad[mask] /= area_w[mask, None]
        return grad

    def _laplacian_u(self, u, nu):
        """
        Viscous acceleration  ν K u / M  (P1 FEM Laplacian).

        Parameters
        ----------
        u  : ndarray, shape (N, 2)
        nu : float – kinematic viscosity

        Returns
        -------
        a_visc : ndarray, shape (N, 2)
        """
        K = self._assemble_laplacian()
        M = self._lumped_mass()

        a_visc = np.zeros((self.n, 2))
        mask = M > 1.0e-20
        for dim in range(2):
            Ku = K @ u[:, dim]
            a_visc[mask, dim] = nu * Ku[mask] / M[mask]
        return a_visc

    # ------------------------------------------------------------------
    # Pressure solver  (fractional-step projection)
    # ------------------------------------------------------------------

    def _solve_pressure(self, uf_pred, us_pred, alpha_s, dt):
        """
        Solve the pressure Poisson equation for the DYNAMIC pressure.

        The correct incompressibility constraint for a two-fluid mixture is
        the volume conservation:

            ∇ · u_vol = 0,   u_vol = α_f u_f + α_s u_s

        Leading to the pressure Poisson equation:

            ∇²p_dyn = ρ_harm / Δt  ·  ∇ · u*_vol

        where  u*_vol = α_f u*_f + α_s u*_s   (volume-averaged)
               ρ_harm = 1 / (α_f/ρ_f + α_s/ρ_s)  (harmonic mean density)

        The reference hydrostatic pressure ρ_f g (y_ref − y) has already been
        subtracted in the prediction step, so p here is purely the dynamic part.

        Boundary conditions:
          - p_dyn = 0 at free-surface nodes   (Dirichlet)
          - ∂p_dyn/∂n = 0 at walls            (natural, FEM)

        Returns
        -------
        p : ndarray, shape (N,)
        """
        alpha_f   = np.maximum(1.0 - alpha_s, 0.01)

        # Volume-averaged predicted mixture velocity
        u_vol     = alpha_f[:, None] * uf_pred + alpha_s[:, None] * us_pred

        # Harmonic mean density  ρ_harm = 1 / (α_f/ρ_f + α_s/ρ_s)
        inv_rho   = alpha_f / self.rho_f + alpha_s / self.rho_s
        rho_harm  = 1.0 / inv_rho

        # Assemble Laplacian
        K = self._assemble_laplacian()

        # RHS:  f_i = Σ_e  (ρ_harm_e / Δt) · A_e/3 · (∇·u*_vol)_e
        ux    = u_vol[self._idx, 0]          # (ne, 3)
        uy    = u_vol[self._idx, 1]
        div_e = (np.einsum('ei,ei->e', self._b, ux)
                 + np.einsum('ei,ei->e', self._c, uy))

        rho_e = np.mean(rho_harm[self._idx], axis=1)    # (ne,)
        rhs_e = -rho_e / dt * div_e                      # (ne,)  NOTE: minus sign from weak form

        f = np.zeros(self.n)
        for k in range(3):
            np.add.at(f, self._idx[:, k], self._area / 3.0 * rhs_e)

        # Identify free-surface nodes (Dirichlet: p = 0)
        wall_tol = 2.0 * self.params['dx']
        fs_nodes = find_free_surface_nodes(
            self.pos, self._idx, self.Lx, self.Ly, wall_tol)

        # Apply Dirichlet BCs and fix isolated nodes
        M_lump = self._lumped_mass()
        K_mod  = K.tolil()
        for node in fs_nodes:
            K_mod[node, :] = 0.0
            K_mod[:, node] = 0.0
            K_mod[node, node] = 1.0
            f[node] = 0.0

        for node in range(self.n):
            if M_lump[node] < 1.0e-20:
                K_mod[node, :] = 0.0
                K_mod[node, node] = 1.0
                f[node] = 0.0

        K_mod = K_mod.tocsr()

        # Solve sparse linear system
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p = spsolve(K_mod, f)
        except Exception:
            p = np.zeros(self.n)

        if not np.all(np.isfinite(p)):
            p = np.zeros(self.n)

        return p

    # ------------------------------------------------------------------
    # Wall boundary conditions
    # ------------------------------------------------------------------

    def _apply_wall_bc(self):
        """
        Enforce no-penetration at the domain walls (y=0, x=0, x=Lx).
        Particles are reflected back only when they actually cross a wall.
        A tiny tolerance avoids false triggers at the first time step.
        """
        lo = 1.0e-6   # near-zero tolerance: only trigger at the actual wall

        # Left wall (x = 0)
        mask = self.pos[:, 0] < lo
        self.pos[mask, 0] = lo
        self.uf[mask, 0]  = np.maximum(self.uf[mask, 0], 0.0)
        self.us[mask, 0]  = np.maximum(self.us[mask, 0], 0.0)

        # Right wall (x = Lx)
        mask = self.pos[:, 0] > self.Lx - lo
        self.pos[mask, 0] = self.Lx - lo
        self.uf[mask, 0]  = np.minimum(self.uf[mask, 0], 0.0)
        self.us[mask, 0]  = np.minimum(self.us[mask, 0], 0.0)

        # Bottom wall (y = 0)
        mask = self.pos[:, 1] < lo
        self.pos[mask, 1] = lo
        self.uf[mask, 1]  = np.maximum(self.uf[mask, 1], 0.0)
        self.us[mask, 1]  = np.maximum(self.us[mask, 1], 0.0)

        # Top boundary: free surface – no hard wall constraint

    # ------------------------------------------------------------------
    # Main time-step
    # ------------------------------------------------------------------

    def step(self, dt):
        """
        Advance by one time step using the fractional-step PFEM algorithm:

        1.  Mesh reconstruction   (Delaunay + Alpha-Shape)
        2.  Viscosity forces      (FEM Laplacian, explicit)
        3.  Implicit drag coupling (2×2 per-particle system, Eq. 4.52)
            Solves for predicted velocities (uf*, us*) simultaneously,
            avoiding the explicit-drag time-step restriction:
                ρ_f/Δt (uf* − uf) = ρ_f a_visc_f + γ αs/αf (us* − uf*) + ρ_f g
                ρ_s/Δt (us* − us) = ρ_s a_visc_s − γ (us* − uf*)         + ρ_s g
        4.  Pressure Poisson     (FEM on Delaunay mesh)
        5.  Velocity correction  (pressure-gradient projection)
        6.  Volume-fraction update (continuity, Eq. 4.51)
        7.  Lagrangian particle advection
        8.  Wall boundary conditions

        Parameters
        ----------
        dt : float – time-step size [s]
        """
        # ── 1. Mesh reconstruction ──────────────────────────────────
        has_mesh = self._rebuild_mesh()

        if not has_mesh:
            # Degenerate: fall freely under gravity
            self.uf[:, 1] -= self.g * dt
            self.us[:, 1] -= self.g * dt
            alpha_f = np.maximum(1.0 - self.alpha_s, 0.01)
            rho_mix = alpha_f * self.rho_f + self.alpha_s * self.rho_s
            u_mix   = (alpha_f[:, None] * self.rho_f * self.uf
                       + self.alpha_s[:, None] * self.rho_s * self.us
                       ) / rho_mix[:, None]
            self.pos += dt * u_mix
            self._apply_wall_bc()
            self.time += dt
            return

        # ── 2. Viscosity (explicit FEM Laplacian) ───────────────────
        a_visc_f = self._laplacian_u(self.uf, self.nu_f)
        a_visc_s = self._laplacian_u(self.us, self.nu_s)

        # ── 3. Implicit drag coupling ───────────────────────────────
        # Compute drag coefficient γ from current relative velocity
        alpha_f = np.maximum(1.0 - self.alpha_s, 0.01)
        u_rel   = self.uf - self.us
        gamma   = self._gamma(u_rel, self.alpha_s)          # (N,)

        # ── Effective body forces  (fixed-reference hydrostatic) ────
        # We subtract the reference hydrostatic gradient  ∇p_ref = ρ_f g ê_y
        # from the body force so that pure water in equilibrium stays still:
        #   g_eff_f = g − ∇p_ref/ρ_f = (0,−g) − (0,−g) = 0
        #   g_eff_s = g − ∇p_ref/ρ_s = (0,−g) − (0,−ρ_f g/ρ_s)
        #           = (0, −g (1 − ρ_f/ρ_s))   ← buoyancy-adjusted gravity
        # These are CONSTANTS (independent of local concentration).
        g_eff_fy  = 0.0
        g_eff_sy  = -self.g * (1.0 - self.rho_f / self.rho_s)   # ≈ −6.11 m/s²
        g_eff_f   = np.array([0.0, g_eff_fy])
        g_eff_s   = np.array([0.0, g_eff_sy])

        # Per-particle 2×2 system  (solved analytically):
        #   A uf* + B us* = rf
        #   C uf* + D us* = rs
        # where the RHS already includes viscosity + buoyancy-reduced gravity.
        beta_f = gamma * self.alpha_s / (alpha_f + 1.0e-8)  # (N,)
        beta_s = gamma                                        # (N,)

        A = self.rho_f / dt + beta_f       # (N,)
        B = -beta_f                         # (N,)
        C = -beta_s                         # (N,)
        D = self.rho_s / dt + beta_s       # (N,)

        det = A * D - B * C                 # (N,) always positive

        # RHS: ρ/Δt * u^n + ρ * (a_visc + g_eff)  [buoyancy-reduced gravity]
        rf = (self.rho_f / dt * self.uf
              + self.rho_f * (a_visc_f + g_eff_f))           # (N, 2)
        rs = (self.rho_s / dt * self.us
              + self.rho_s * (a_visc_s + g_eff_s))           # (N, 2)

        # Cramer's rule (component-wise, same det for both x and y)
        uf_pred = (D[:, None] * rf - B[:, None] * rs) / det[:, None]
        us_pred = (A[:, None] * rs - C[:, None] * rf) / det[:, None]

        # ── 4. Pressure Poisson solve ───────────────────────────────
        p_new = self._solve_pressure(uf_pred, us_pred,
                                     self.alpha_s, dt)
        self.pressure = p_new

        # ── 5. Velocity correction ──────────────────────────────────
        grad_p = self._gradient(p_new)

        self.uf = uf_pred - (dt / self.rho_f) * grad_p
        self.us = us_pred - (dt / self.rho_s) * grad_p

        # ── 6. Volume-fraction update ───────────────────────────────
        # Continuity (Eq. 4.51):  Dα_s/Dt = −α_s ∇·u_s
        div_us       = self._divergence(self.us)
        self.alpha_s = np.clip(
            self.alpha_s - dt * self.alpha_s * div_us,
            0.0, 0.85,
        )

        # ── 7. Lagrangian particle advection ────────────────────────
        alpha_f = np.maximum(1.0 - self.alpha_s, 0.01)
        rho_mix = alpha_f * self.rho_f + self.alpha_s * self.rho_s

        u_mix = (alpha_f[:, None] * self.rho_f * self.uf
                 + self.alpha_s[:, None] * self.rho_s * self.us
                 ) / rho_mix[:, None]

        self.pos += dt * u_mix

        # ── 8. Boundary conditions ──────────────────────────────────
        self._apply_wall_bc()

        self.time += dt
