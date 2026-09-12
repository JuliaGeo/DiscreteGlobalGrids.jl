module RasterizeTests
using Test, Statistics
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeoInterface as GI

@testset "DGGS rasterization and masks" begin
    grid = DGG.levelgrid(DGG.HEALPixSystem(), 1)
    p = GI.Point((10.0, 25.0))
    q = GI.Point((-130.0, -35.0))
    points = [p, p, q]
    dest = DGG.rasterize(sum, points; to=grid, fill=[2, 4, 8])
    ip = only(DGG._raster_indices(grid, p))
    iq = only(DGG._raster_indices(grid, q))
    @test dest[ip] == 6
    @test dest[iq] == 8
    @test count(ismissing, dest) == length(dest) - 2
    @test DGG.rasterize(mean, points; to=grid, fill=[2,4,8])[ip] == 3
    @test DGG.rasterize(count, points; to=grid)[ip] == 2
    @test DGG.rasterize(first, points; to=grid, fill=[2,4,8])[ip] == 2
    @test DGG.rasterize(last, points; to=grid, fill=[2,4,8])[ip] == 4
    @test DGG.rasterize(extrema, points; to=grid, fill=[2,4,8])[ip] == (2,4)
    @test DGG.rasterize(x -> join(x, ","), points; to=grid, fill=["a","b","c"])[ip] == "a,b"
    @test_throws ArgumentError DGG.rasterize(points; to=grid, fill=1)
    @test_throws ArgumentError DGG.rasterize(sum, points; to=grid)
    @test_throws DimensionMismatch DGG.rasterize(sum, points; to=grid, fill=[1,2])
    @test_throws ArgumentError DGG.rasterize(sum, points; to=grid, fill=1, filename="no.tif")

    ids = DGG.rasterize(points; to=grid, op=vcat, fill=[[1],[2],[3]],
        eltype=Vector{Int}, init=Int[], missingval=Int[])
    @test ids[ip] == [1,2]
    @test ids[iq] == [3]
    emptyids = findall(isempty, parent(ids))
    push!(ids[first(emptyids)], 99)
    @test isempty(ids[last(emptyids)])
    mutated = DGG.rasterize(points; to=grid, op=append!, fill=[[1],[2],[3]],
        eltype=Vector{Int}, init=Int[], missingval=Int[])
    @test mutated[ip] == [1,2]
    @test mutated[iq] == [3]
    @test all(isempty, mutated[emptyids])
    DGG.rasterize!(mutated, [p]; op=append!, fill=[[4]], init=Int[], missingval=Int[])
    @test mutated[ip] == [1,2,4]
    @test mutated[iq] == [3]

    template = DD.DimArray(zeros(Int, 2, DGG.ncells(grid), 3),
        (DD.Ti(1:2), DGG.Cells(DGG.CellLookup(grid)), DD.Dim{:band}(1:3)))
    cube = DGG.rasterize(sum, points; to=template, fill=[2,4,8], missingval=0)
    @test size(cube) == size(template)
    @test all(==(6), parent(cube)[:,ip,:])
    @test all(==(8), parent(cube)[:,iq,:])
    DGG.rasterize!(cube, [p,p]; op=+, fill=[1,2], missingval=0)
    @test all(==(9), parent(cube)[:,ip,:])
    @test all(==(8), parent(cube)[:,iq,:])
    functional = DGG.rasterize(points; to=grid, fill=x -> x + 1, init=0, missingval=0)
    @test functional[ip] == 2
    @test functional[iq] == 1
    DGG.rasterize!(functional, [p]; fill=x -> x + 3)
    @test functional[ip] == 5

    stack = DGG.rasterize(sum, points; to=grid, fill=(population=[2,4,8], area=1))
    @test stack[:population][ip] == 6
    @test stack[:area][ip] == 2
    @test DGG.boolmask(points; to=grid)[ip]
    @test !DGG.boolmask([p]; to=grid)[iq]
    @test DGG.boolmask([p]; to=grid, invert=true)[iq]
    @test ismissing(DGG.missingmask([p]; to=grid)[iq])
    masked = DGG.mask(dest; with=[p])
    @test masked[ip] == 6
    @test ismissing(masked[iq])
    @test dest[iq] == 8
    @test all(iszero, DGG.rasterize(count, GI.Point[]; to=grid))
    numeric = DGG.rasterize(count, points; to=grid)
    DGG.rasterize!(count, numeric, [p,p])
    @test numeric[ip] == 4
    DGG.rasterize!(sum, numeric, [p,p]; fill=[3,4])
    @test numeric[ip] == 11
    @test DGG.rasterize(sum, [p,p]; to=grid, fill=Int8[100,100])[ip] === 200
    @test DGG.rasterize(sum, [p,p]; to=grid, op=*, fill=[3,4])[ip] == 12
    table = (geometry=points, population=[2,4,8], area=[1,1,1])
    @test DGG.rasterize(sum, table; to=grid, fill=:population)[ip] == 6
    tablelayers = DGG.rasterize(sum, table; to=grid, fill=(:population,:area))
    @test tablelayers[:population][ip] == 6
    @test tablelayers[:area][ip] == 2
    booldest = DD.DimArray(falses(length(dest)), DD.dims(dest))
    @test DGG.boolmask!(booldest, [p]) === booldest
    @test booldest[ip] && !booldest[iq]
    missdest = DD.DimArray(Union{Missing,Bool}[false for _ in dest], DD.dims(dest))
    DGG.missingmask!(missdest, [p])
    @test missdest[ip] === true
    @test missdest[iq] === missing
    @test DD.name(DGG.rasterize(sum, points; to=grid, fill=1, name=:burn)) == :burn
    @test isequal(DGG.rasterize(sum, points; to=grid, fill=[2,4,8], threaded=true), dest)
end
end
