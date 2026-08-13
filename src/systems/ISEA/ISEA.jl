"""
    ISEA

Shared icosahedral machinery for the ISEA family of discrete global grids: the
spherical icosahedron (vertex/neighbor tables, per-base counterclockwise rings,
the [`Orientation`](@ref) placement rotation) and the Snyder equal-area
projection charts (per-face plane, forward/inverse, dev-frame slot maps).

Nothing here knows about any particular aperture or indexing scheme — this is
the common basis the ISEA systems build on: ISEA7H / IGEO7 today
([`IGeo7`](@ref)), ISEA3H / ISEA4H / ISEA9H later. Everything lives in the
*grid frame*, the standard ISEA placement; the identity orientation therefore
returns coordinates in the standard ISEA frame directly.

| file             | contents                                                |
|:-----------------|:--------------------------------------------------------|
| `icosahedron.jl` | sphere helpers, standard ISEA vertex/neighbor tables, `Orientation` |
| `snyder.jl`      | per-face Snyder ISEA plane + dev-frame slot maps        |

`snyder.jl` uses `icosahedron.jl`'s vertex tables and vector helpers and nothing
else; both are stdlib-only floating-point geometry. Provenance of every constant
and convention is cited inline and recorded in `docs/IGeo7/`.
"""
module ISEA

include("icosahedron.jl")
include("snyder.jl")

# Public surface: the names the ISEA grid systems and the test suites consume.
# Build-time helpers (`_make_faces`, `_corner_pos`, ...) and the intermediate
# Snyder constants stay internal.
export ADJ_DOT,
    COS_LG,
    DEV_CONE_DEG,
    DEV_SLOTS,
    DEV_SLOT_OF_FACE,
    DevSlot,
    FACES,
    FACE_TRIPLES,
    Face,
    ISEA_LAT_HI,
    ISEA_LON0,
    L_PLANE,
    NBASE,
    NBRS_CCW,
    NEIGHBORS,
    ORIENT_IDENTITY,
    Orientation,
    REFERENCE_EDGE,
    R_AUTHALIC,
    R_EA,
    R_EA2,
    SNY_COTT,
    SNY_SECTOR,
    SQRT3,
    TAN_LG,
    VERTICES,
    angdist,
    dev_angle_deg,
    dev_slot_index,
    dev_to_xyz,
    face_to_dev,
    from_grid,
    lonlat_to_xyz,
    nearest_vertex,
    snyder_fwd,
    snyder_inv_xyz,
    to_grid,
    vadd,
    vcross,
    vdot,
    vertex,
    vnorm,
    vnormalize,
    vscale,
    vsub,
    xyz_to_lonlat

end # module ISEA
