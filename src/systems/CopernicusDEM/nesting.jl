# Cross-resolution nesting: how one Copernicus DEM lattice sits inside another.
#
# This is a cross-SYSTEM relation, not a hierarchy edge. `levels(sys) == 0:1` on both
# sides and neither system is an ancestor of the other; `parent`/`children` say nothing
# about GLO-30 inside GLO-90. What relates the two is that they are the same construction
# at two values of `N`: the same six-band table, the same tile lattice, and — the part
# that is measured rather than assumed — the same anchor for their pixel CENTRES.

"""
    nesting_factor(coarse::CopernicusDEMSystem, fine::CopernicusDEMSystem) -> Int

`k`, the integer ratio `lat_intervals(fine) ÷ lat_intervals(coarse)`; `3` for the shipped
pair (GLO-30 inside GLO-90).

Throws an `ArgumentError` naming both `N`s when the ratio is not an exact integer. One
check covers the columns too: `ncols = 2N / factor2` with the same six `factor2`s on both
sides, so `ncols_fine(r) == k * ncols_coarse(r)` in every tile row as soon as `k` is an
integer — the band table cannot go out of step.
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

**Module-local**, so call it as `DiscreteGlobalGrids.CopernicusDEM.refine`: it is not
exported and not a method of any `DiscreteGlobalGrids` generic. The module docstring
has why, and what it would collide with.

The cells of `fine` that `c`, a cell of `coarse`, refines into: one tile at level 0, and
`k²` pixels at level 1, where `k = lat_intervals(fine) ÷ lat_intervals(coarse)`
([`nesting_factor`](@ref)) — `k = 3` for the shipped GLO-30-inside-GLO-90 pair.

A level-0 tile refines to the tile of `fine` with the **same lower-left corner**; the two
lattices carry the same 64 800 tiles in the same order. A level-1 pixel at raster
`(j, i)` refines to the `k²` pixels `(k*j .. k*j+k-1, k*i .. k*i+k-1)` of that tile,
returned in this system's raster order — north row first, west to east — which is
**ascending** in id and (for `k > 1`) not contiguous: each of the `k` rows is a run of `k`
consecutive ids, and successive runs are `ncols(fine, r)` apart.

The block never leaves the tile (`k*j + k - 1 < k*N_coarse = N_fine`), so
`parent(fine, f) == only(refine(coarse, fine, parent(coarse, c)))` for every `f` returned.

# What nests is the POST lattice

Both products anchor their first pixel centre on the same integer-degree point, use the
same six-band reduction table, and have column counts that divide without remainder, so
the coarse post at `(j, i)` and the fine post at `(k*j, k*i)` are the same point on the
sphere. That co-location is the one thing here that is measured rather than argued: the
test suite sweeps it as `worst_post`, and asserts `< 1e-12` degrees against the 4.2e-4
degrees (1.5″) a half-pixel registration slip would move a post. The `k²` fine posts of
a block are the coarse post and the `k² - 1` fine posts immediately south and east of it.

!!! warning "The cell BOXES do not tile the coarse box"
    Both products are pixel-is-point, so a cell's box is its post ± half of **that
    product's own** pixel — and the half-pixel outset is `k` times smaller on the fine
    side. The `k²` fine boxes therefore tile a box of the right SIZE that is the coarse
    box translated south-east, on each axis, by `Δ_coarse·(1 - 1/k)/2` — half a coarse
    pixel less half a fine one, equivalently `(k-1)/2` fine pixels, which at `k = 3`
    comes out as exactly one whole GLO-30 pixel.

    That shift is a FRACTION of a coarse pixel, so its arcsecond value is per-band and
    no single figure states it. In longitude it is 1.0″ in the `[0, 50)` band, 1.5″ in
    `[50, 60)` and 10″ in `[85, 90)`; in latitude it is 1.0″ in every band, because
    `Δlat` does not vary by band. The test suite pins the fraction rather than the
    arcseconds — `worst_shift`, which measures the west and north edges against
    `(1 - 1/k)/2` and reads 2.3e-11 for the shipped pair.

    Their areas sum to that translated box exactly (to 3.6e-16 relative, logged as
    `worst_union`), and so miss the coarse cell's own area by about `tan(φ)·1″` — logged
    as `gap_lat50` and `gap_equator`: 6.0e-6 in the latitude-50 tile row and 8.5e-8 in
    the equator one.

    No uniform, tile-local index scheme fixes this. An exact box tiling would need the
    fine cells centred on the coarse box, which is `k*i - (k-1)/2 …` — not an integer
    block for even `k`, and ragged in any case: it would reach into the neighbouring
    tile at every tile edge, and across a band boundary the tile above has a different
    `ncols`, whose column edges do not fall on the coarse ones at all. The GLO-90 pixel
    at column 100 of an `N49` tile straddles latitude 50, and its west edge lands on
    GLO-30 column **199.5** of the `[50, 60)` band above it. Exact cell-box nesting does
    not exist in this lattice; exact post nesting does, and that is what this function
    returns.

    The `lat_s = 89` and `lat_s = -90` tile rows are further apart still, because
    [`cell_box`](@ref)'s clamp to `+90` and extension to `-90` are each half of the
    system's own pixel: there the block's area differs from the coarse cell's by a
    relative gap `|block - coarse| / coarse` — the convention every gap on this page is
    quoted in — of **1.78** in the `+90` row and **0.40** in the `-90` row, logged as
    `worst_pole_n` and `worst_pole_s`. Their posts still coincide; nothing else about
    them does.

!!! warning "A grid hierarchy, not a value hierarchy"
    GLO-90 is an independently produced mission product resampled from WorldDEM, not a
    3x3 mean of public GLO-30: co-located posts of the two products carry different
    elevations, by about a metre in the one pair that was read off AWS by hand. Nothing
    in `test/` reads a GLO-30 raster, so take that figure as an order of magnitude and
    not a measurement this repo reproduces. Use this pair to relate CELLS; aggregate
    VALUES yourself if you need them consistent.
"""
function refine(coarse::CopernicusDEMSystem, fine::CopernicusDEMSystem,
        c::DGG.LevelIndex)
    k = nesting_factor(coarse, fine)
    r, q, j, i = decode(coarse, c)
    # Level 0 is the same tile lattice on both sides — 180 rows of 360 tiles, north to
    # south then west to east — so this is the identity on the id. Decoded and re-encoded
    # rather than returned as-is, because what is being asserted is the CORNER.
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

**Module-local**, so call it as `DiscreteGlobalGrids.CopernicusDEM.coarsen`; same
reasons as [`refine`](@ref).

The single cell of `coarse` that `c`, a cell of `fine`, refines from: the tile with the
same lower-left corner at level 0, and the pixel `(j ÷ k, i ÷ k)` of that tile at level 1,
with `k = nesting_factor(coarse, fine)`.

The exact inverse of [`refine`](@ref) in both directions:
`c in refine(coarse, fine, coarsen(fine, coarse, c))` for every cell of `fine`, and
`coarsen(fine, coarse, f) == c` for every `f` in `refine(coarse, fine, c)`.

Read [`refine`](@ref)'s two warnings before using either: this maps posts, and neither
cell boxes nor elevation values.
"""
function coarsen(fine::CopernicusDEMSystem, coarse::CopernicusDEMSystem,
        c::DGG.LevelIndex)
    k = nesting_factor(coarse, fine)
    r, q, j, i = decode(fine, c)
    DGG.level(c) == 0 && return DGG.LevelIndex(0, tileordinal(r, q))
    return DGG.LevelIndex(1, tilebase(coarse, r, q) +
                             Int64(j ÷ k) * ncols(coarse, r) + Int64(i ÷ k))
end
