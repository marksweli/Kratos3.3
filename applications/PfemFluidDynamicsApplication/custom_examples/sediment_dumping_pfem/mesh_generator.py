"""
mesh_generator.py
=================
Mesh generation module for the PFEM sediment-dumping simulation.

Provides:
  - Regular-grid particle initialisation
  - Delaunay triangulation (via scipy)
  - Alpha-Shape filtering (circumradius criterion)
  - Boundary / free-surface node detection
"""

import numpy as np
from scipy.spatial import Delaunay


# ---------------------------------------------------------------------------
# Particle initialisation
# ---------------------------------------------------------------------------

def create_particle_grid(x_range, y_range, dx, dy=None):
    """
    Create a regular grid of particles.

    Parameters
    ----------
    x_range : (float, float)  –  (x_min, x_max)
    y_range : (float, float)  –  (y_min, y_max)
    dx      : float           –  horizontal particle spacing
    dy      : float, optional –  vertical spacing (default = dx)

    Returns
    -------
    positions : ndarray, shape (N, 2)
    """
    if dy is None:
        dy = dx
    x = np.arange(x_range[0] + dx / 2.0, x_range[1], dx)
    y = np.arange(y_range[0] + dy / 2.0, y_range[1], dy)
    X, Y = np.meshgrid(x, y)
    return np.column_stack([X.ravel(), Y.ravel()])


# ---------------------------------------------------------------------------
# Delaunay triangulation
# ---------------------------------------------------------------------------

def delaunay_triangulate(points):
    """
    Perform Delaunay triangulation of a 2-D point cloud.

    Parameters
    ----------
    points : ndarray, shape (N, 2)

    Returns
    -------
    tri      : scipy Delaunay object
    simplices: ndarray, shape (M, 3)  – triangle vertex indices
    """
    tri = Delaunay(points)
    return tri, tri.simplices.copy()


# ---------------------------------------------------------------------------
# Alpha-Shape filtering
# ---------------------------------------------------------------------------

def _circumradius(pts):
    """
    Circumradius of a 2-D triangle.

    Parameters
    ----------
    pts : ndarray, shape (3, 2)

    Returns
    -------
    R : float
    """
    v1 = pts[1] - pts[0]
    v2 = pts[2] - pts[0]
    area2 = abs(v1[0] * v2[1] - v1[1] * v2[0])   # 2 * area
    if area2 < 1.0e-20:
        return np.inf
    a = np.linalg.norm(pts[1] - pts[0])
    b = np.linalg.norm(pts[2] - pts[1])
    c = np.linalg.norm(pts[0] - pts[2])
    return (a * b * c) / (2.0 * area2)


def alpha_shape_filter(points, simplices, alpha):
    """
    Keep only triangles whose circumradius R satisfies R <= 1/alpha.

    Parameters
    ----------
    points    : ndarray, shape (N, 2)
    simplices : ndarray, shape (M, 3)
    alpha     : float   – alpha parameter [m⁻¹]; larger → more concave boundary

    Returns
    -------
    filtered : ndarray, shape (K, 3)
    """
    if len(simplices) == 0:
        return np.empty((0, 3), dtype=int)

    max_R = 1.0 / alpha
    keep = []
    for s in simplices:
        R = _circumradius(points[s])
        if R <= max_R:
            keep.append(s)

    return np.array(keep, dtype=int) if keep else np.empty((0, 3), dtype=int)


# ---------------------------------------------------------------------------
# Boundary / free-surface detection
# ---------------------------------------------------------------------------

def find_boundary_nodes(simplices, n_nodes):
    """
    Return indices of nodes that lie on the boundary of the triangulation
    (edges shared by exactly one triangle).

    Parameters
    ----------
    simplices : ndarray, shape (M, 3)
    n_nodes   : int

    Returns
    -------
    boundary : set of int
    """
    edge_count = {}
    for s in simplices:
        for a, b in ((s[0], s[1]), (s[1], s[2]), (s[0], s[2])):
            key = (min(a, b), max(a, b))
            edge_count[key] = edge_count.get(key, 0) + 1

    boundary = set()
    for (a, b), cnt in edge_count.items():
        if cnt == 1:
            boundary.add(a)
            boundary.add(b)
    return boundary


def find_free_surface_nodes(positions, simplices, Lx, Ly, wall_tol):
    """
    Identify free-surface nodes (boundary nodes that are not on a solid wall).

    Solid walls:  x ≈ 0,  x ≈ Lx,  y ≈ 0  (three sides of the tank).
    Free surface: top boundary of the particle cloud (y direction).

    Parameters
    ----------
    positions : ndarray, shape (N, 2)
    simplices : ndarray, shape (M, 3)
    Lx        : float  – domain width
    Ly        : float  – domain height
    wall_tol  : float  – tolerance for wall detection

    Returns
    -------
    free_surface : set of int
    """
    if len(simplices) == 0:
        return set()

    boundary = find_boundary_nodes(simplices, len(positions))

    free_surface = set()
    for node in boundary:
        x, y = positions[node]
        on_wall = (x < wall_tol or x > Lx - wall_tol or y < wall_tol)
        if not on_wall:
            free_surface.add(node)

    return free_surface


# ---------------------------------------------------------------------------
# Vectorised element-data extraction  (used by solver)
# ---------------------------------------------------------------------------

def build_element_data(points, simplices):
    """
    Pre-compute shape-function gradients, areas, and node indices for every
    valid triangular element.

    For a linear (P1) triangle with vertices at (x1,y1),(x2,y2),(x3,y3)
    and area A the shape-function gradients are constant:

        b_i = (y_{i+1} - y_{i+2}) / (2A)
        c_i = (x_{i+2} - x_{i+1}) / (2A)

    (indices taken mod 3).

    Parameters
    ----------
    points    : ndarray, shape (N, 2)
    simplices : ndarray, shape (M, 3)

    Returns
    -------
    idx  : ndarray, shape (K, 3)   – vertex indices of valid elements
    area : ndarray, shape (K,)
    b    : ndarray, shape (K, 3)   – x-derivatives of shape functions
    c    : ndarray, shape (K, 3)   – y-derivatives of shape functions
    """
    if len(simplices) == 0:
        empty = np.empty((0, 3))
        return (np.empty((0, 3), dtype=int),
                np.empty(0),
                empty.copy(), empty.copy())

    x = points[simplices, 0]   # (M, 3)
    y = points[simplices, 1]   # (M, 3)

    # Twice the signed area (positive for counter-clockwise ordering)
    two_A = ((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
             - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
    area = np.abs(two_A) / 2.0

    # Discard degenerate elements
    valid = area > 1.0e-20
    simplices = simplices[valid]
    area = area[valid]
    x = x[valid]
    y = y[valid]

    # Shape-function gradients
    b = np.column_stack([
        y[:, 1] - y[:, 2],
        y[:, 2] - y[:, 0],
        y[:, 0] - y[:, 1],
    ]) / (2.0 * area[:, np.newaxis])

    c = np.column_stack([
        x[:, 2] - x[:, 1],
        x[:, 0] - x[:, 2],
        x[:, 1] - x[:, 0],
    ]) / (2.0 * area[:, np.newaxis])

    return simplices, area, b, c
