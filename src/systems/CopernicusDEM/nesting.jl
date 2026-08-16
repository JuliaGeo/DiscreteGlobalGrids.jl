# Cross-resolution nesting: how one Copernicus DEM lattice sits inside another. A
# cross-SYSTEM relation, not a hierarchy edge — `parent`/`children` say nothing
# about GLO-30 inside GLO-90.

"""
    nesting_factor(coarse::CopernicusDEMSystem, fine::CopernicusDEMSystem) -> Int

`k`, the integer ratio `lat_intervals(fine) ÷ lat_intervals(coarse)`; `3` for the shipped
pair (GLO-30 inside GLO-90).

Throws an `ArgumentError` when the ratio is not an exact integer. The same check
covers the columns: `ncols_fine(r) == k * ncols_coarse(r)` in every tile row.
"""
@inline function nesting_factor(::CopernicusDEMSystem{NC},
        ::CopernicusDEMSystem{NF}) where {NC,NF}
    NF % NC == 0 || throw(ArgumentError(
        "a Copernicus DEM lattice with N = $NF is not an integer refinement of one " *
        "with N = $NC: $NF / $NC is not a whole number, so no cell of the coarse " *
        "lattice is a whole number of fine cells"))
    return Int(NF ÷ NC)
end

"""
    refine(coarse::CopernicusDEMSystem, fine::CopernicusDEMSystem, c) -> Vector{LevelIndex}

**Module-local**: call it as `DiscreteGlobalGrids.CopernicusDEM.refine`.

The cells of `fine` that `c`, a cell of `coarse`, refines into: the tile with the same
lower-left corner at level 0, and `k²` pixels at level 1, with
`k = `[`nesting_factor`](@ref) — `3` for the shipped GLO-30-inside-GLO-90 pair. The
pixels come in raster order, ascending in id but not contiguous: `k` runs of `k`
consecutive ids, `ncols(fine, r)` apart. The block never leaves the tile.

What nests is the POST lattice: the coarse post at `(j, i)` and the fine post at
`(k*j, k*i)` are the same point on the sphere.

!!! warning "The cell BOXES do not tile the coarse box"
    Both products are pixel-is-point, so the half-pixel outset is `k` times smaller
    on the fine side: the `k²` fine boxes tile a box of the right size translated
    south-east by `(k-1)/2` fine pixels. Exact cell-box nesting does not exist in
    this lattice; exact post nesting does, and that is what this function returns.

!!! warning "A grid hierarchy, not a value hierarchy"
    GLO-90 is resampled from WorldDEM, not a 3x3 mean of public GLO-30, so co-located
    posts carry different elevations. Use this pair to relate CELLS; aggregate VALUES
    yourself if you need them consistent.
"""
function refine(coarse::CopernicusDEMSystem, fine::CopernicusDEMSystem,
        c::DGG.LevelIndex)
    k = nesting_factor(coarse, fine)
    r, q, j, i = decode(coarse, c)
    # Level 0 is the same tile lattice on both sides, so this is the identity on the
    # id — decoded and re-encoded because the relation is between tile corners.
    DGG.level(c) == 0 && return [DGG.LevelIndex(0, tileordinal(r, q))]
    nc = ncols(fine, r)
    base = tilebase(fine, r, q) + Int64(k) * Int64(j) * nc + Int64(k) * Int64(i)
    out = Vector{DGG.LevelIndex}(undef, k * k)
    n = 0
    for a in 0:(k - 1)
        rowstart = base + Int64(a) * nc
        for b in 0:(k - 1)
            out[n += 1] = DGG.LevelIndex(1, rowstart + Int64(b))
        end
    end
    return out
end

"""
    coarsen(fine::CopernicusDEMSystem, coarse::CopernicusDEMSystem, c) -> LevelIndex

**Module-local**: call it as `DiscreteGlobalGrids.CopernicusDEM.coarsen`.

The single cell of `coarse` that `c`, a cell of `fine`, refines from: the tile with the
same lower-left corner at level 0, and the pixel `(j ÷ k, i ÷ k)` of that tile at level
1. The exact inverse of [`refine`](@ref) in both directions; [`refine`](@ref)'s two
warnings apply here too.
"""
function coarsen(fine::CopernicusDEMSystem, coarse::CopernicusDEMSystem,
        c::DGG.LevelIndex)
    k = nesting_factor(coarse, fine)
    r, q, j, i = decode(fine, c)
    DGG.level(c) == 0 && return DGG.LevelIndex(0, tileordinal(r, q))
    return DGG.LevelIndex(1, tilebase(coarse, r, q) +
                             Int64(j ÷ k) * ncols(coarse, r) + Int64(i ÷ k))
end
