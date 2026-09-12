module RasterAPIRegressionTests

using Test, Statistics
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeoInterface as GI
import Rasters as RA

struct EmptyPolygon end
GI.trait(::EmptyPolygon) = GI.PolygonTrait()
GI.isempty(::EmptyPolygon) = true

@testset "Raster API regression coverage" begin
    grid = DGG.levelgrid(DGG.HEALPixSystem(), 1)
    p = GI.Point((10.0, 25.0))
    q = GI.Point((-130.0, -35.0))
    ip = only(DGG._raster_indices(grid, p))
    iq = only(DGG._raster_indices(grid, q))

    @testset "iterable and mutable cell fills" begin
        tuplecells = DGG.rasterize(p; to=grid, fill=(1, 2))
        @test tuplecells[ip] == (1, 2)
        @test ismissing(tuplecells[iq])

        vectorcells = DGG.rasterize(p; to=grid, fill=[1, 2],
            eltype=Vector{Int}, missingval=Int[])
        @test vectorcells[ip] == [1, 2]
        emptyindices = findall(isempty, parent(vectorcells))
        @test length(emptyindices) > 1
        @test vectorcells[emptyindices[1]] !== vectorcells[emptyindices[2]]
        push!(vectorcells[emptyindices[1]], 99)
        @test isempty(vectorcells[emptyindices[2]])

        setcells = DGG.rasterize(p; to=grid, fill=Set((1, 2)),
            eltype=Set{Int}, missingval=Set{Int}())
        @test setcells[ip] == Set((1, 2))

        # Match Rasters: a length-one iterable repeats over a geometry collection.
        repeated = DGG.rasterize(sum, [p, p, q]; to=grid, fill=[7])
        @test repeated[ip] == 14
        @test repeated[iq] == 7

        update = x -> (push!(x, length(x) + 1); x)
        mutable = DGG.rasterize([p, p, q]; to=grid, fill=update,
            init=Int[], eltype=Vector{Int}, missingval=Int[], threaded=true)
        @test mutable[ip] == [1, 2]
        @test mutable[iq] == [1]
        untouched = findall(isempty, parent(mutable))
        @test mutable[untouched[1]] !== mutable[untouched[2]]

        objecttemplate = DD.DimArray(
            reshape([Int[] for _ in 1:(DGG.ncells(grid) * 2)],
                DGG.ncells(grid), 2),
            (DGG.Cells(DGG.CellLookup(grid)), DD.Ti(1:2)),
        )
        broadcastcells = DGG.rasterize(p; to=objecttemplate, fill=[3, 4],
            eltype=Vector{Int}, missingval=Int[], threaded=true)
        @test broadcastcells[ip, 1] == [3, 4]
        @test broadcastcells[ip, 2] == [3, 4]
        @test broadcastcells[ip, 1] !== broadcastcells[ip, 2]
        DGG.rasterize!(broadcastcells, p; fill=update, init=Int[], missingval=Int[])
        @test broadcastcells[ip, 1] == [3, 4, 3]
        @test broadcastcells[ip, 2] == [3, 4, 3]
        @test broadcastcells[ip, 1] !== broadcastcells[ip, 2]
    end

    @testset "Tables normalization" begin
        rows = [(geometry=p, value=2), (geometry=p, value=4),
            (geometry=q, value=8)]
        rowburn = DGG.rasterize(sum, rows; to=grid, fill=:value)
        @test rowburn[ip] == 6
        @test rowburn[iq] == 8

        columns = (longitude=[10.0, 10.0, -130.0],
            latitude=[25.0, 25.0, -35.0], value=[2, 4, 8])
        columnburn = DGG.rasterize(sum, columns; to=grid,
            geometrycolumn=(:longitude, :latitude), fill=:value)
        @test isequal(columnburn, rowburn)

        # A NamedTuple column table also has a GeoInterface FeatureTrait. It is
        # nevertheless three table rows, including for zonal and extraction.
        pointcolumns = (geometry=[p, missing, q], value=[2, 4, 8])
        @test DGG.rasterize(sum, pointcolumns; to=grid, fill=:value)[ip] == 2
        tableextract = DGG.extract(rowburn, pointcolumns; id=true, index=true)
        @test getproperty.(tableextract, :id) == [1, 2, 3]
        @test ismissing(tableextract[2].index)
        @test length(DGG.zonal(sum, rowburn; of=pointcolumns)) == 3
    end

    @testset "Rasters stack options and metadata" begin
        cellaxis = DGG.Cells(DGG.CellLookup(grid))
        a = RA.Raster(fill(-99.0, DGG.ncells(grid)), (cellaxis,);
            missingval=-99.0, name=:a, metadata=Dict("layer" => "a"))
        b = RA.Raster(fill(-88, DGG.ncells(grid)), (cellaxis,);
            missingval=-88, name=:b, metadata=Dict("layer" => "b"))
        refs = (DD.Ti(2000),)
        stack = RA.RasterStack((a=a, b=b);
            metadata=Dict("stack" => "kept"), refdims=refs)

        burned = DGG.rasterize(sum, [p]; to=stack, fill=1,
            init=(a=Float32(10), b=Int16(20)),
            eltype=(a=Float32, b=Int16),
            missingval=(a=Float32(-9), b=Int16(-8)))
        @test burned[:a][ip] === Float32(11)
        @test burned[:b][ip] === Int16(21)
        @test burned[:a][iq] === Float32(-9)
        @test burned[:b][iq] === Int16(-8)
        @test RA.missingval(burned) == (a=Float32(-9), b=Int16(-8))
        @test DD.metadata(burned) == DD.metadata(stack)
        @test DD.metadata(burned[:a]) == DD.metadata(a)
        @test DD.metadata(burned[:b]) == DD.metadata(b)
        @test DD.refdims(burned) == refs

        DGG.rasterize!(sum, burned, q; fill=1,
            init=(a=Float32(10), b=Int16(20)),
            missingval=(a=Float32(-9), b=Int16(-8)))
        @test burned[:a][iq] === Float32(11)
        @test burned[:b][iq] === Int16(21)

        a[ip] = 3
        b[ip] = 4
        zonedstack = DGG.zonal(sum, stack; of=[p])
        @test zonedstack isa RA.AbstractRasterStack
        @test zonedstack[:a][1] == 3
        @test zonedstack[:b][1] == 4
        @test DD.metadata(zonedstack) == DD.metadata(stack)
        @test DD.metadata(zonedstack[:a]) == DD.metadata(a)
        @test DD.metadata(zonedstack[:b]) == DD.metadata(b)
        @test DD.refdims(zonedstack) == refs

        nostoredmissing = RA.Raster(zeros(Int, DGG.ncells(grid)), (cellaxis,);
            missingval=nothing, name=:plain)
        inferred = DGG.rasterize(sum, [p]; to=nostoredmissing, fill=2)
        @test RA.missingval(inferred) === missing
        @test inferred[ip] == 2
        @test ismissing(inferred[iq])
        @test eltype(inferred) == Union{Missing, Int}

        @test_throws ArgumentError DGG.rasterize(sum, [p]; to=a, fill=2,
            crs=RA.EPSG(4326))
        @test_throws ArgumentError DGG.rasterize(sum, [p]; to=a, fill=2,
            mappedcrs=RA.EPSG(4326))
    end

    @testset "multidimensional extract and zonal empties" begin
        stored = sort([DGG.cellindex(grid, ip), DGG.cellindex(grid, iq)])
        partial = DGG.PartialGrid(DGG.system(grid), 1, stored)
        pip = DGG.localindex(partial, DGG.cellindex(grid, ip))
        piq = DGG.localindex(partial, DGG.cellindex(grid, iq))
        data = Array{Union{Missing,Int}}(missing, 2, 2, 2)
        data[:, pip, :] = [1 2; missing 4]
        data[:, piq, :] = [5 6; 7 8]
        cube = DD.DimArray(data,
            (DD.Ti(1:2), DGG.Cells(DGG.CellLookup(partial)), DD.Dim{:Band}(1:2));
            name=:signal)

        extracted = DGG.extract(cube, [p, q, missing];
            id=true, index=true, skipmissing=false)
        @test length(extracted) == 3
        @test DD.dims(extracted[1].signal) ==
            DD.dims(cube, (DD.Ti, DD.Dim{:Band}))
        @test isequal(parent(extracted[1].signal), data[:, pip, :])
        @test isequal(parent(extracted[2].signal), data[:, piq, :])
        @test ismissing(extracted[3].signal)
        @test getproperty.(DGG.extract(cube, [p, q, missing];
            id=true, skipmissing=true), :id) == [2]

        onecell = DGG.PartialGrid(DGG.system(grid), 1,
            [DGG.cellindex(grid, ip)])
        scalar = DD.DimArray([42], DGG.Cells(DGG.CellLookup(onecell)); name=:value)
        # GeoInterface's convenience wrapper constructor currently indexes the
        # first ring, so use its public concrete representation for zero rings.
        rings = Vector{GI.LinearRing}()
        emptywrapper = GI.Polygon{false,false,typeof(rings),Nothing,Nothing}(
            rings, nothing, nothing)
        for emptygeom in (EmptyPolygon(), emptywrapper)
            @test all(ismissing, DGG.rasterize(emptygeom; to=grid, fill=1))
            @test isempty(DGG.extract(scalar, emptygeom))
            @test DGG.zonal(sum, scalar; of=emptygeom) == 0
            @test DGG.zonal(mean, scalar; of=emptygeom, emptyval=-7) == -7
        end
        @test ismissing(DGG.zonal(sum, scalar; of=q))
        @test ismissing(DGG.zonal(sum, scalar; of=q, emptyval=-1))
        mixed = DGG.zonal(sum, scalar; of=[q, p], threaded=true)
        @test isequal(parent(mixed), Union{Missing,Int}[missing, 42])

        sliced = DD.DimArray(reshape([1, 2, 3, 4], 2, 1, 2),
            (DD.Ti(1:2), DGG.Cells(DGG.CellLookup(onecell)), DD.Dim{:Band}(1:2));
            name=:sliced)
        alloutside = DGG.zonal(sum, sliced; of=[q, q], threaded=true)
        @test size(alloutside) == (2, 2, 2)
        @test all(ismissing, alloutside)
        @test DD.dims(alloutside)[1:2] == DD.dims(sliced, (DD.Ti, DD.Dim{:Band}))

        delta = 1e-4
        tiny = GI.Polygon([[
            (10.0 - delta, 25.0 - delta), (10.0 + delta, 25.0 - delta),
            (10.0 + delta, 25.0 + delta), (10.0 - delta, 25.0 + delta),
            (10.0 - delta, 25.0 - delta),
        ]])
        @test isempty(DGG._raster_indices(onecell, tiny; boundary=:center))
        @test !isempty(DGG._raster_indices(onecell, tiny; boundary=:intersects))
        emptyfirst = DGG.zonal(sum, sliced; of=[tiny, p])
        @test all(iszero, parent(emptyfirst)[:, :, 1])
        @test parent(emptyfirst)[:, :, 2] == dropdims(parent(sliced); dims=2)
        @test all(==(-7), parent(DGG.zonal(sum, sliced;
            of=[tiny], emptyval=-7))[:, :, 1])

        mutableempty = DGG.zonal(collect, sliced; of=[tiny, tiny], emptyval=Int[])
        @test mutableempty[1, 1, 1] !== mutableempty[2, 1, 1]
        push!(mutableempty[1, 1, 1], 99)
        @test isempty(mutableempty[2, 1, 1])

        rr = RA.Raster(parent(sliced), DD.dims(sliced);
            missingval=-99, name=:sliced, metadata=Dict("units" => "x"))
        outside_raster = DGG.zonal(sum, rr; of=q)
        @test outside_raster isa RA.AbstractRaster
        @test size(outside_raster) == (2, 2)
        @test all(ismissing, outside_raster)
        @test RA.missingval(outside_raster) === missing
        mixed_raster = DGG.zonal(sum, rr; of=[q, p])
        @test mixed_raster isa RA.AbstractRaster
        @test RA.missingval(mixed_raster) === missing
        @test all(ismissing, parent(mixed_raster)[:, :, 1])
        @test parent(mixed_raster)[:, :, 2] == dropdims(parent(rr); dims=2)
        @test DD.metadata(mixed_raster) == DD.metadata(rr)
    end
end

end
