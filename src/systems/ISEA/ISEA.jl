"""
    ISEA

Shared icosahedron geometry and Snyder equal-area charts for ISEA grid systems.
All coordinates use the standard ISEA grid frame; [`Orientation`](@ref) maps
between that frame and world coordinates.
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
