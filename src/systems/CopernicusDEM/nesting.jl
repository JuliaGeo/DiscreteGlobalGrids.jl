# Cross-resolution nesting, separate from `parent` and `children`.

"""
    nesting_factor(coarse::CopernicusDEMSystem, fine::CopernicusDEMSystem) -> Int

The integer ratio `lat_intervals(fine) ÷ lat_intervals(coarse)`.

Throws `ArgumentError` unless the ratio is exact; column counts scale by the same factor.
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

Returns the same tile at level 0 or `k²` fine pixels at level 1, where `k` is
[`nesting_factor`](@ref). Pixels are in raster order as `k` non-contiguous runs.

The post at coarse `(j, i)` coincides with fine `(k*j, k*i)`.

!!! warning "The cell BOXES do not tile the coarse box"
    Pixel-is-point registration shifts the fine block southeast by `(k-1)/2`
    fine pixels. Posts nest exactly; cell boxes do not.

!!! warning "A grid hierarchy, not a value hierarchy"
    GLO-90 is not a 3x3 mean of public GLO-30. This relates cells, not values.
"""
function refine(coarse::CopernicusDEMSystem, fine::CopernicusDEMSystem,
        c::DGG.LevelIndex)
    k = nesting_factor(coarse, fine)
    r, q, j, i = decode(coarse, c)
    # Level 0 is the same tile lattice at both resolutions.
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

Returns the same tile at level 0 or coarse pixel `(j ÷ k, i ÷ k)` at level 1.
This is the exact inverse of [`refine`](@ref); its warnings also apply.
"""
function coarsen(fine::CopernicusDEMSystem, coarse::CopernicusDEMSystem,
        c::DGG.LevelIndex)
    k = nesting_factor(coarse, fine)
    r, q, j, i = decode(fine, c)
    DGG.level(c) == 0 && return DGG.LevelIndex(0, tileordinal(r, q))
    return DGG.LevelIndex(1, tilebase(coarse, r, q) +
                             Int64(j ÷ k) * ncols(coarse, r) + Int64(i ÷ k))
end
