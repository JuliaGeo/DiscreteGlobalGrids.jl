# The cross-language check: a store written with `target = :xdggs` opens in the
# real xdggs, through `xdggs.decode`, and its cell centres agree with this
# package's. Runs when `DGG_XDGGS_PYTHON` names a Python interpreter that has
# xarray, zarr and xdggs installed, and skips otherwise: the Julia suite has no
# Python of its own, and everything the store must contain for xdggs is also
# asserted byte-for-byte in `write.jl`.

module DGGSIOXdggsPythonTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: HEALPixSystem, CellVector, CellLookup, Cells, LevelIndex,
    levelgrid, ncells, cellindex, cell_centroid
import DimensionalData as DD

const PYTHON = get(ENV, "DGG_XDGGS_PYTHON", "")
const ZARR_LOADED = try
    @eval using Zarr
    true
catch
    false
end

if isempty(PYTHON) || !ZARR_LOADED
    @info "DGG_XDGGS_PYTHON is unset or Zarr.jl is unavailable: the xdggs decode check is skipped."
else

    const SCRIPT = joinpath(@__DIR__, "xdggs_decode.py")

    # `key=value` lines into a dictionary; the script's own `error=` line is
    # what a failing run reports through.
    function decode_with_xdggs(path)
        out = IOBuffer()
        err = IOBuffer()
        ok = success(pipeline(`$PYTHON $SCRIPT $path`; stdout=out, stderr=err))
        facts = Dict{String,String}()
        for line in eachline(IOBuffer(take!(out)))
            k, v = split(line, "="; limit=2)
            facts[k] = v
        end
        ok || error("xdggs could not open $path:\n" * get(facts, "error", "") *
                    "\n" * String(take!(err)))
        return facts
    end

    lonlat(p) = (atand(p[2], p[1]), asind(p[3]))

    @testset "xdggs opens a target = :xdggs store" begin
        sys = HEALPixSystem()
        grid = levelgrid(sys, 3)
        cells = [cellindex(grid, i) for i in 1:ncells(grid)]
        tas = DD.DimArray(Float32.(1:length(cells)),
            (Cells(CellLookup(CellVector(sys, 3, cells))),); name=:tas)
        path = DGG.dggwrite(joinpath(mktempdir(), "healpix.zarr"), tas;
            target=:xdggs, chunks=256)

        facts = decode_with_xdggs(path)
        @test facts["grid_name"] == "healpix"
        @test facts["level"] == "3"
        @test facts["indexing_scheme"] == "nested"
        @test facts["ncells"] == string(ncells(grid))
        @test facts["first_id"] == "0"
        @test occursin("tas", facts["data_vars"])

        lon, lat = lonlat(Tuple(cell_centroid(grid, LevelIndex(3, 0))))
        @test parse(Float64, facts["first_lon"]) ≈ lon atol = 1e-9
        @test parse(Float64, facts["first_lat"]) ≈ lat atol = 1e-9
    end

    @testset "the default encoding is what xdggs cannot open" begin
        # The guard the target exists for: `:auto` chooses ranges for this
        # axis, and xdggs's own convention finds no `cell_ids` there.
        sys = HEALPixSystem()
        grid = levelgrid(sys, 2)
        cells = [cellindex(grid, i) for i in 1:ncells(grid)]
        tas = DD.DimArray(Float32.(1:length(cells)),
            (Cells(CellLookup(CellVector(sys, 2, cells))),); name=:tas)
        path = DGG.dggwrite(joinpath(mktempdir(), "ranges.zarr"), tas; chunks=64)
        @test_throws ErrorException decode_with_xdggs(path)
    end

end # PYTHON && ZARR_LOADED

end # module
