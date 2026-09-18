# The CopDEM done ledger: what counts as written, and what the totals count.

module CopDEMStoreTests

using Test
import Printf: @sprintf
import DiscreteGlobalGrids as DGG

say(parts...) = nothing
stamp() = "t"
include(joinpath(@__DIR__, "..", "..", "scripts", "copdem_store.jl"))

@testset "done ledger" begin
    mktempdir() do dir
        path = joinpath(dir, "store.zarr.done.ndjson")
        log = DoneLog(path)
        record!(log, 7, 100, 10, 1.0, 1)
        record!(log, 7, 100, 40, 1.0, 2)     # recomputed after a resume
        record!(log, 9, 100, 0, 1.0, 1)
        close(log)

        # A chunk logged twice is one chunk, at its last line's values.
        led = doneledger(path)
        @test length(led) == 2 && led[7] == (cells = 100, nan = 40)

        # A written chunk with no ledger line is counted but unaccounted.
        t = storetotals(path, Set([7, 9, 11]))
        @test t == (chunks = 3, cells = 200, nan = 40, unaccounted = 1)

        # With no chunk listing on disk, the ledger alone says what is done.
        @test donechunks(path, joinpath(dir, "store.zarr"), "elevation") == Set([7, 9])
    end
end

@testset "chunk files name their chunk" begin
    mktempdir() do dir
        layer = joinpath(dir, "elevation")
        mkpath(layer)
        touch(joinpath(layer, ".zarray"))
        touch(joinpath(layer, "4.0"))
        @test storechunks(dir, "elevation") == Set([5])
        touch(joinpath(layer, "unexpected"))
        @test storechunks(dir, "elevation") === nothing
    end
end

@testset "chunk list round trip" begin
    mktempdir() do dir
        path = joinpath(dir, "store.zarr.columns.txt")
        @test load_chunklist(path) === nothing
        save_chunklist(path, 5, [3, 1759, 12235])
        @test load_chunklist(path) == [3, 1759, 12235]
    end
end

end # module
