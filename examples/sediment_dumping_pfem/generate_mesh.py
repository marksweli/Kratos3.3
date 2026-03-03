"""
generate_mesh.py
~~~~~~~~~~~~~~~~
Programmatic mesh generator for the 2D sediment dumping (水平抛泥) PFEM example.

Generates a structured triangular mesh for a 0.7m x 0.7m rectangular domain with:
  - Water (fluid) region: the entire domain filled with water up to y = 0.7m
  - Sediment patch:  x in [0.485, 0.515], y in [0.68, 0.70]  (upper-left of water body)
  - Rigid walls:     bottom (y=0), left (x=0), right (x=0.7) boundaries
  - Free surface:    top boundary (y=0.7)

Sub-model parts produced
------------------------
  Parts_Water        – fluid (water) nodes/elements
  Parts_Sediment     – sediment-concentration patch nodes/elements
  Parts_Walls        – no-slip wall condition nodes
  Parts_FreeSurface  – free-surface condition nodes

Usage
-----
  python generate_mesh.py               # uses default resolution
  python generate_mesh.py --nx 35 --ny 35

The script writes  SedimentDumping.mdpa  in the same directory.
"""

import argparse
import os
import sys

# ---------------------------------------------------------------------------
# Domain & patch parameters
# ---------------------------------------------------------------------------
DOMAIN_X   = (0.0, 0.7)   # m
DOMAIN_Y   = (0.0, 0.7)   # m

# Sediment patch (upper surface flush with free surface at y = 0.7)
SED_X      = (0.485, 0.515)  # m
SED_Y      = (0.68,  0.70)   # m


# ---------------------------------------------------------------------------
def build_mesh(nx: int, ny: int):
    """Return nodes, triangles and sub-model part memberships.

    Parameters
    ----------
    nx, ny : int
        Number of grid *intervals* in x- and y-directions.

    Returns
    -------
    nodes : list of (id, x, y)
    elements : list of (id, prop_id, n1, n2, n3)
    walls_node_ids : list of int   – bottom + left + right wall nodes
    free_surface_node_ids : list of int
    water_element_ids : list of int
    sediment_element_ids : list of int
    water_node_ids : list of int
    sediment_node_ids : list of int
    """
    x0, x1 = DOMAIN_X
    y0, y1 = DOMAIN_Y

    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny

    # ------------------------------------------------------------------
    # 1.  Build node grid
    # ------------------------------------------------------------------
    node_index = {}   # (ix, iy) -> node_id (1-based)
    nodes = []
    node_id = 1
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            x = x0 + ix * dx
            y = y0 + iy * dy
            node_index[(ix, iy)] = node_id
            nodes.append((node_id, x, y))
            node_id += 1

    # ------------------------------------------------------------------
    # 2.  Build triangular elements  (2 triangles per quad cell)
    #     Lower triangle: (i,j)-(i+1,j)-(i+1,j+1)
    #     Upper triangle: (i,j)-(i+1,j+1)-(i,j+1)
    # ------------------------------------------------------------------
    elements = []
    elem_id = 1

    # property 1 = water,  property 2 = sediment
    def _is_in_sediment_patch(cx, cy):
        return (SED_X[0] <= cx <= SED_X[1]) and (SED_Y[0] <= cy <= SED_Y[1])

    water_element_ids    = []
    sediment_element_ids = []

    for iy in range(ny):
        for ix in range(nx):
            n00 = node_index[(ix,   iy  )]
            n10 = node_index[(ix+1, iy  )]
            n11 = node_index[(ix+1, iy+1)]
            n01 = node_index[(ix,   iy+1)]

            # centroid of lower triangle
            cx_l = (nodes[n00-1][1] + nodes[n10-1][1] + nodes[n11-1][1]) / 3.0
            cy_l = (nodes[n00-1][2] + nodes[n10-1][2] + nodes[n11-1][2]) / 3.0
            prop_l = 2 if _is_in_sediment_patch(cx_l, cy_l) else 1

            # centroid of upper triangle
            cx_u = (nodes[n00-1][1] + nodes[n11-1][1] + nodes[n01-1][1]) / 3.0
            cy_u = (nodes[n00-1][2] + nodes[n11-1][2] + nodes[n01-1][2]) / 3.0
            prop_u = 2 if _is_in_sediment_patch(cx_u, cy_u) else 1

            elements.append((elem_id, prop_l, n00, n10, n11))
            if prop_l == 2:
                sediment_element_ids.append(elem_id)
            else:
                water_element_ids.append(elem_id)
            elem_id += 1

            elements.append((elem_id, prop_u, n00, n11, n01))
            if prop_u == 2:
                sediment_element_ids.append(elem_id)
            else:
                water_element_ids.append(elem_id)
            elem_id += 1

    # ------------------------------------------------------------------
    # 3.  Identify boundary node sets
    # ------------------------------------------------------------------
    walls_node_ids        = []
    free_surface_node_ids = []
    water_node_ids_set    = set()
    sediment_node_ids_set = set()

    for node_id_n, x, y in nodes:
        ix_pos = round((x - x0) / dx)
        iy_pos = round((y - y0) / dy)

        # Rigid walls: bottom, left, right
        if iy_pos == 0 or ix_pos == 0 or ix_pos == nx:
            walls_node_ids.append(node_id_n)

        # Free surface: top boundary
        if iy_pos == ny:
            free_surface_node_ids.append(node_id_n)

    # Assign nodes to water / sediment based on element connectivity
    for elem_id_e, prop_id, n1, n2, n3 in elements:
        if prop_id == 2:
            sediment_node_ids_set.update([n1, n2, n3])
        else:
            water_node_ids_set.update([n1, n2, n3])

    # Nodes that belong to sediment patch only (not shared with water)
    # Actually keep all nodes; water nodes = all non-sediment-exclusive nodes
    water_node_ids    = sorted(water_node_ids_set)
    sediment_node_ids = sorted(sediment_node_ids_set)

    return (nodes, elements,
            sorted(walls_node_ids), sorted(free_surface_node_ids),
            water_element_ids, sediment_element_ids,
            water_node_ids, sediment_node_ids)


# ---------------------------------------------------------------------------
def write_mdpa(output_path, nx, ny):
    """Write the Kratos .mdpa file."""
    (nodes, elements,
     walls_node_ids, free_surface_node_ids,
     water_element_ids, sediment_element_ids,
     water_node_ids, sediment_node_ids) = build_mesh(nx, ny)

    lines = []

    # --- Header ---
    lines.append("Begin ModelPartData")
    lines.append("//  VARIABLE_NAME value")
    lines.append("End ModelPartData")
    lines.append("")
    lines.append("Begin Properties 0")
    lines.append("End Properties")
    lines.append("")

    # --- Nodes ---
    lines.append("Begin Nodes")
    for nid, x, y in nodes:
        lines.append(f"    {nid:6d}   {x:.10f}   {y:.10f}   0.0000000000")
    lines.append("End Nodes")
    lines.append("")

    # --- Elements (fluid part only – wall is rigid, no elements) ---
    lines.append("Begin Elements TwoStepUpdatedLagrangianVPFluidElement2D")
    for eid, prop_id, n1, n2, n3 in elements:
        lines.append(f"    {eid:6d}   {prop_id}   {n1}   {n2}   {n3}")
    lines.append("End Elements")
    lines.append("")

    # ------------------------------------------------------------------
    # Sub-model parts
    # ------------------------------------------------------------------

    # ---- Parts_Water (fluid body – water) ----------------------------
    lines.append("Begin SubModelPart Parts_Water")
    lines.append("  Begin SubModelPartNodes")
    for nid in water_node_ids:
        lines.append(f"    {nid}")
    lines.append("  End SubModelPartNodes")
    lines.append("  Begin SubModelPartElements")
    for eid in water_element_ids:
        lines.append(f"    {eid}")
    lines.append("  End SubModelPartElements")
    lines.append("End SubModelPart")
    lines.append("")

    # ---- Parts_Sediment (sediment concentration patch) ---------------
    lines.append("Begin SubModelPart Parts_Sediment")
    lines.append("  Begin SubModelPartNodes")
    for nid in sediment_node_ids:
        lines.append(f"    {nid}")
    lines.append("  End SubModelPartNodes")
    lines.append("  Begin SubModelPartElements")
    for eid in sediment_element_ids:
        lines.append(f"    {eid}")
    lines.append("  End SubModelPartElements")
    lines.append("End SubModelPart")
    lines.append("")

    # ---- Parts_Walls (rigid boundary – no-slip) ----------------------
    # Rigid body: use a single dummy element (element 0 concept) or
    # represent as a node-only sub-model part for boundary conditions.
    lines.append("Begin SubModelPart Parts_Walls")
    lines.append("  Begin SubModelPartNodes")
    for nid in walls_node_ids:
        lines.append(f"    {nid}")
    lines.append("  End SubModelPartNodes")
    lines.append("End SubModelPart")
    lines.append("")

    # ---- VELOCITY BC – no-slip on walls ------------------------------
    lines.append("Begin SubModelPart VELOCITY_Velocity_Walls")
    lines.append("  Begin SubModelPartNodes")
    for nid in walls_node_ids:
        lines.append(f"    {nid}")
    lines.append("  End SubModelPartNodes")
    lines.append("End SubModelPart")
    lines.append("")

    with open(output_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"Wrote {output_path}")
    print(f"  Nodes    : {len(nodes)}")
    print(f"  Elements : {len(elements)}")
    print(f"    Water  : {len(water_element_ids)}")
    print(f"    Sediment: {len(sediment_element_ids)}")
    print(f"  Wall nodes      : {len(walls_node_ids)}")
    print(f"  Free-surface nds: {len(free_surface_node_ids)}")


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate SedimentDumping.mdpa for the PFEM sediment-dumping example"
    )
    parser.add_argument("--nx", type=int, default=35,
                        help="Number of grid intervals in x-direction (default: 35)")
    parser.add_argument("--ny", type=int, default=35,
                        help="Number of grid intervals in y-direction (default: 35)")
    parser.add_argument("--output", type=str,
                        default=os.path.join(os.path.dirname(__file__), "SedimentDumping.mdpa"),
                        help="Output .mdpa file path")
    args = parser.parse_args()

    write_mdpa(args.output, args.nx, args.ny)


if __name__ == "__main__":
    main()
