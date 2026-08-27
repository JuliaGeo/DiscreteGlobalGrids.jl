# A stored cell axis IS a region, and answers the region verbs with the same
# code a computed one does.
#
# What these check is the bridge `idrank(grid, id) + 1 == globalindex(grid, c)`
# and nothing else: the store is written from a cube whose axis is a
# `CellLookup`, read back as a `ChunkedCellLookup`, and the two are required to
# answer identically. A conversion that dropped a cell, shifted an index by
# one, or left the runs unmerged fails here rather than somewhere downstream.

module DGGSIORegionTests

using Test
import DiscreteGlobalGrids as DGG

const HAS_ZARR = try
    @eval using Zarr
    true
catch err
    @warn "Zarr.jl is not loadable: the stored-region suite is skipped." exception = err
    false
end

if HAS_ZARR

import DimensionalData as DD
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, levelgrid, cellindex,
    ncells, Cells, CellVector, CellLookup, ChunkedCellLookup, dggread, dggwrite,
    mapneighbors, region, halo, border, interior, adjacency, PartialGrid
using DiscreteGlobalGrids.Engine: nwindows, windows

# The store, and the in-memory cube it was written from.
function roundtrip(dir, name, sys, level, cells; kw...)
    A = DD.DimArray(Float64.(eachindex(cells)),
        Cells(CellLookup(CellVector(sys, level, cells))); name=:e)
    path = joinpath(dir, name)
    dggwrite(path, DD.DimStack((; e=A)); kw...)
    return A, dggread(path)[:e]
end

# Every law that has to hold between a cube and the store it round-tripped
# through, whatever encoding the store used.
function same_region(label, A, B)
    @testset "$label" begin
        want, got = DD.lookup(A, Cells), DD.lookup(B, Cells)
        @test got isa ChunkedCellLookup
        cv = region(got)

        # The conversion names the same cells in the same order. Index order
        # is what lets a result computed through the twin be written back
        # against the store's own axis without a permutation.
        @test cv isa CellVector
        @test collect(cv) == collect(parent(want))
        @test collect(got) == collect(want)

        # As compressed as building it in memory. A conversion that emitted one
        # window per stored interval rather than merging the ones that abut
        # answers every question correctly and fails this.
        @test nwindows(windows(cv)) == nwindows(windows(parent(want)))

        # `region` is a memo, not a rebuild: the same object comes back.
        @test region(got) === cv

        # The four region verbs, and the sweep, agree cell for cell.
        @test collect(halo(got)) == collect(halo(want))
        @test collect(border(got)) == collect(border(want))
        @test collect(interior(got)) == collect(interior(want))
        @test adjacency(got) == adjacency(want)
        pg, pgw = PartialGrid(got), PartialGrid(want)
        @test DGG.ncells(pg) == DGG.ncells(pgw)
        @test [cellindex(pg, i) for i in 1:DGG.ncells(pg)] ==
              [cellindex(pgw, i) for i in 1:DGG.ncells(pgw)]

        f = (c, nbrs) -> length(nbrs) * DGG.localindex(c)
        @test mapneighbors(f, got) == mapneighbors(f, want)
        # `==` on cubes compares lookups too, and those differ by design.
        @test parent(mapneighbors(f, B)) == parent(mapneighbors(f, A))
    end
end

@testset "a stored axis answers the region verbs as the cube it came from" begin
    mktempdir() do dir
        sys, L = IGeo7System(), 4
        grid = levelgrid(sys, L)
        # Holes and separated runs: enough stored intervals that the merge is
        # exercised, and enough boundary that the halo is not trivial.
        cells = [cellindex(grid, p) for p in [50:149; 399:399; 699:848; 2000:2400]]

        A, B = roundtrip(dir, "ranges.zarr", sys, L, cells; encoding=:ranges)
        same_region("ranges", A, B)

        # The one encoding that has to read its ids to convert.
        A, B = roundtrip(dir, "dense.zarr", sys, L, cells; encoding=:dense)
        same_region("dense", A, B)

        # Implicit needs the whole level, which is also the single-window case.
        whole = [cellindex(grid, p) for p in 1:ncells(grid)]
        A, B = roundtrip(dir, "implicit.zarr", sys, L, whole; encoding=:implicit)
        same_region("implicit, the whole level", A, B)
        @test nwindows(windows(region(DD.lookup(B, Cells)))) == 1

        # A second id scheme, because the bridge is rank arithmetic and rank is
        # the system's.
        hsys, hL = HEALPixSystem(), 3
        hgrid = levelgrid(hsys, hL)
        hcells = [cellindex(hgrid, p) for p in [10:200; 400:600]]
        A, B = roundtrip(dir, "hpx.zarr", hsys, hL, hcells; encoding=:ranges)
        same_region("HEALPix, ranges", A, B)
    end
end

end # if HAS_ZARR

end # module DGGSIORegionTests
