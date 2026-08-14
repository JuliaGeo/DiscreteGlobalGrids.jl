# runtests.jl — DGGS-Zarr reader entry point.
#
# Part of the DiscreteGlobalGrids.jl test suite; run it from the package's test
# environment, either through the top-level test/runtests.jl or directly:
#   julia --project=<DiscreteGlobalGrids.jl> test/DGGSZarr/runtests.jl
#
# The suite lives in its own module, like every other, so the generic vocabulary
# (`cell_centers`, `cell_boundaries`, ...) cannot collide with a system's.
#
# Two tiers of test:
#
#   * `test_monotonic.jl` is self-contained — it builds its ids from the Z7
#     codecs and needs no data on disk, so it runs everywhere and is where the
#     number-line contract is actually pinned;
#   * `test_archives.jl` reads the real archives written by the Python IGEO7/Z7
#     tooling, which live outside this repository. It is skipped, loudly, when
#     they are absent, so a fresh checkout still runs green.

module DGGSZarrTestSuite

using Test

using DiscreteGlobalGrids
using DiscreteGlobalGrids: Helpers, ISEA, IGeo7
using DiscreteGlobalGrids.DGGSZarr

"""
    ARCHIVE_ROOT

Directory holding the reference archives, overridable through the
`DGGS_ZARR_TEST_DATA` environment variable so the suite can be pointed at a
copy. The default is the sibling paper repository the archives were written in.
"""
const ARCHIVE_ROOT = get(ENV, "DGGS_ZARR_TEST_DATA",
    joinpath(homedir(), "dev", "build", "igeo7_z7_xarray_paper", "data", "working"))

"The four reference archives, as (ranges, dense) pairs at one level each."
const ARCHIVE_PAIRS = (
    (level=10, ranges="pori_z7_r10_ranges.zarr", dense="pori_z7_r10.zarr"),
    (level=12, ranges="pori_z7_r12_ranges.zarr", dense="pori_z7_r12.zarr"),
)

archive_path(name) = joinpath(ARCHIVE_ROOT, name)

"`true` when every reference archive is present, so the on-disk tier can run."
have_archives() = all(ARCHIVE_PAIRS) do pair
    isdir(archive_path(pair.ranges)) && isdir(archive_path(pair.dense))
end

const TEST_FILES = (
    "test_monotonic.jl",
    "test_archives.jl",
)

@testset "DGGSZarr" begin
    for f in TEST_FILES
        path = joinpath(@__DIR__, f)
        if isfile(path)
            include(path)
        else
            @info "skipping $f (not present)"
        end
    end
end

end # module DGGSZarrTestSuite
