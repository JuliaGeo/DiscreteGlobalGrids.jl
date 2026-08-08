# test/core/test_globe_ids.jl — the lazy globe-complete id vector
# (src/core/globe_ids.jl; docs/design/full_globe_lookups.md §1.1, §1.4).
#
# Four claims, none of which any single system owns:
#
#   * `DGGSGlobeIds(system, level)` *is* the explicit ascending id vector. At a
#     level small enough to enumerate it is checked element for element against
#     `ordinal_to_cell`, in both directions, for every wired system.
#   * it is O(1) to build and to measure at any level, including levels whose
#     explicit vector cannot exist — H3 res 15 is 5.7e14 ids, some 4.6 PB. A
#     materializing implementation does not fail these tests slowly, it fails
#     them with an `OutOfMemoryError`, which is the point of running them here.
#   * it degrades: any non-scalar index falls through to `AbstractArray`'s
#     default and materializes an ordinary `Vector`. That, and nothing else, is
#     what turns a sliced globe lookup back into a partial one (§1.4).
#   * neither `strictly_increasing` nor `show` reads the elements — the two
#     places a globe would otherwise be walked or printed id by id.
#
# Suite lives in its own module: `num_cells` is deliberately unexported from the
# package (the system submodules own that name), so it is qualified here and
# never enters the shared test namespace.
module GlobeIdsTests

using Test

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids

# One entry per kernel-wired system, at a level whose every cell can be
# enumerated. The pair spans both id models: `Int64` ordinal ids (HEALPix) and
# `UInt64` structural ids (H3, A5, IGEO7), where the ordinal numbering is a
# native computation rather than `id + 1`.
const SMALL_GLOBES = (
    (H3DGGS(), 0),
    (A5DGGS(), 1),
    (IGEO7DGGS(), 1),
    (HEALPixDGGS(), 1),
)

# The same systems at a level whose explicit id vector is unrepresentable.
const HUGE_GLOBES = (
    (H3DGGS(), 15),
    (A5DGGS(), 25),
    (IGEO7DGGS(), 19),
    (HEALPixDGGS(), 29),
)

@testset "DGGSGlobeIds is the explicit id vector" begin
    for (system, level) in SMALL_GLOBES
        g = DGGSGlobeIds(system, level)
        total = Int(DGG.num_cells(system, level))
        @test g isa AbstractVector{DGG.cell_id_type(system)}
        @test size(g) == (total,)
        @test length(g) == total
        # The whole claim of the type, checked by enumeration: element `i` is
        # the `i`-th cell of the level in ascending canonical-id order.
        @test collect(g) == [DGG.ordinal_to_cell(system, level, i) for i in 1:total]
        @test issorted(g; lt=(<=))
        @test all(i -> DGG.cell_to_ordinal(system, level, g[i]) == i, eachindex(g))
        @test first(g) == DGG.ordinal_to_cell(system, level, 1)
        @test last(g) == DGG.ordinal_to_cell(system, level, total)
        # Out of range is the kernel's own `OrdinalRangeError`, not a
        # `BoundsError`: `getindex(::DGGSGlobeIds, ::Int)` is `ordinal_to_cell`,
        # whose guard names the level and its ordinal range.
        @test_throws DGG.OrdinalRangeError g[total + 1]
        @test_throws DGG.OrdinalRangeError g[0]
    end
    # H3's res-0 globe against the native enumeration, which is where the
    # ascending order the ordinal contract promises actually comes from.
    @test collect(DGGSGlobeIds(H3DGGS(), 0)) ==
          sort(collect(DGG.H3.H3Native.res0_cells()))
end

@testset "DGGSGlobeIds costs O(1) at any level" begin
    for (system, level) in HUGE_GLOBES
        g = DGGSGlobeIds(system, level)
        @test length(g) == Int(DGG.num_cells(system, level))
        @test g[1] == DGG.ordinal_to_cell(system, level, 1)
        @test g[end] == DGG.ordinal_to_cell(system, level, length(g))
        @test g[2] > g[1]
    end
    @test length(DGGSGlobeIds(H3DGGS(), 15)) == 569707381193162

    # Construction validates what it can without a wired kernel — the
    # `max_level` bound — and nothing that needs the cell count.
    @test_throws ArgumentError DGGSGlobeIds(H3DGGS(), 16)
    @test_throws ArgumentError DGGSGlobeIds(H3DGGS(), -1)
    # A5's res-30 encoding gap: a level the id format can hold but no uniform
    # grid occupies constructs here and throws at the first operation that asks
    # for the count, exactly as `DGGSGrid` does.
    @test DGGSGlobeIds(A5DGGS(), 30) isa DGGSGlobeIds
    @test_throws ArgumentError length(DGGSGlobeIds(A5DGGS(), 30))
    # ...and a registered-but-unwired system likewise constructs, then reports
    # the missing wiring rather than a silent wrong default.
    @test DGGSGlobeIds(RHEALPixDGGS(), 40) isa DGGSGlobeIds
    @test_throws DGG.NotPortedError length(DGGSGlobeIds(RHEALPixDGGS(), 40))
end

@testset "DGGSGlobeIds degrades to a materialized Vector" begin
    for (system, level) in SMALL_GLOBES
        g = DGGSGlobeIds(system, level)
        ID = DGG.cell_id_type(system)
        head = g[1:10]
        @test head isa Vector{ID}
        @test head == [DGG.ordinal_to_cell(system, level, i) for i in 1:10]
        @test g[[2, 5]] == [g[2], g[5]]
        @test collect(g) isa Vector{ID}
    end
    # The same on a globe no `collect` could survive: slicing is where the
    # laziness is meant to stop, so it must stop at the slice's own size.
    @test DGGSGlobeIds(H3DGGS(), 15)[1:10] isa Vector{UInt64}
end

@testset "DGGSGlobeIds is never walked element by element" begin
    # The O(N) check every `<X>Lookup` inner constructor runs, answered from
    # the type: `cell_to_ordinal` is contractually strictly monotone in the id,
    # so its inverse enumerated over `1:num_cells` is strictly ascending by
    # definition. Walking it would be 5.7e14 kernel calls.
    for (system, level) in HUGE_GLOBES
        @test DGG.Helpers.strictly_increasing(DGGSGlobeIds(system, level))
    end

    # `show` prints the description, never the ids. Base's `print_matrix`
    # allocates a row range spanning the whole array unless the caller set
    # `:limit`, so a leak here is an `OutOfMemoryError`, not a slow print.
    g = DGGSGlobeIds(H3DGGS(), 15)
    for text in (sprint(show, g), sprint(show, MIME"text/plain"(), g))
        @test occursin("DGGSGlobeIds", text)
        @test occursin("H3DGGS", text)          # the system, however qualified
        @test endswith(text, ", 15)")           # ...and the level
        @test length(text) < 60                 # ...and no ids at all
    end
end

end # module GlobeIdsTests
