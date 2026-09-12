module RasterExtractZonalTests
using Test, Statistics
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeoInterface as GI

@testset "Extraction, normalization and multidimensional zonal" begin
    grid = DGG.levelgrid(DGG.HEALPixSystem(),2)
    p = GI.Point((10.0,25.0)); q = GI.Point((-130.0,-35.0))
    ip = only(DGG._raster_indices(grid,p)); iq = only(DGG._raster_indices(grid,q))
    cube = DD.DimArray(reshape(collect(1.0:2DGG.ncells(grid)*3),2,DGG.ncells(grid),3),
        (DD.Ti(1:2),DGG.Cells(DGG.CellLookup(grid)),DD.Dim{:Band}(1:3)); name=:temperature)
    rows = DGG.extract(cube,[p,q]; id=true,index=true)
    @test length(rows)==2
    @test rows[1].index==ip
    @test rows[1].geometry==(10.0,25.0)
    @test parent(rows[1].temperature)==parent(cube)[:,ip,:]
    @test DD.dims(rows[1].temperature)==DD.dims(cube,(DD.Ti,DD.Dim{:Band}))
    poly = GI.Polygon([[(0.0,10.0),(35.0,10.0),(35.0,40.0),(0.0,40.0),(0.0,10.0)]])
    extent = DGG.Extents.Extent(X=(0.0,35.0),Y=(10.0,40.0))
    @test DGG.zonal(sum,cube;of=extent,spatialslices=false)==sum(cube[:,DGG._raster_indices(grid,extent),:])
    pole = DGG.Extents.Extent(X=(-180.0,180.0),Y=(70.0,90.0))
    poleids = DGG._raster_indices(grid,pole)
    @test poleids==[i for i in 1:DGG.ncells(grid) if DGG._raster_center(grid,i)[2]>=70]
    ids = DGG._raster_indices(grid,poly)
    zones = DGG.zonal(sum,cube;of=[poly,p],threaded=false)
    @test size(zones)==(2,3,2)
    @test parent(zones)[:,:,1]==dropdims(sum(parent(cube)[:,ids,:];dims=2);dims=2)
    @test parent(zones)[:,:,2]==parent(cube)[:,ip,:]
    @test DGG.zonal(sum,cube;of=p,spatialslices=false)==sum(parent(cube)[:,ip,:])
    @test parent(DGG.zonal(sum,cube;of=p))==parent(cube)[:,ip,:]
    @test size(DGG.zonal(sum,cube;of=GI.Point[]))==(2,3,0)
    @test isequal(DGG.zonal(sum,cube;of=[poly,p],threaded=true),zones)
    @test size(DGG.zonal(sum,cube;of=[p],spatialslices=(DGG.Cells,DD.Ti)))==(3,1)
    @test_throws ArgumentError DGG.zonal(sum,cube;of=[p],spatialslices=(DD.Ti,))
    incomplete = DD.DimArray(Union{Missing,Float64}[missing,2.0],DGG.Cells(DGG.CellLookup(DGG.PartialGrid(DGG.system(grid),2,sort([DGG.cellindex(grid,ip),DGG.cellindex(grid,iq)])))))
    emptyzone=GI.Polygon([[(70.0,0.0),(71.0,0.0),(71.0,1.0),(70.0,1.0),(70.0,0.0)]])
    @test all(==(-1),DGG.zonal(mean,cube;of=[emptyzone],emptyval=-1))
    @test length(DGG.extract(cube,[p,q];flatten=false))==2
    @test length(DGG.extract(cube,[poly,poly];flatten=false))==2
    @test length(DGG.extract(cube,[missing];skipmissing=true))==0
    @test ismissing(only(DGG.extract(cube,[missing])).temperature)
    table=(geometry=[p,p,q],value=[1,2,3])
    @test size(DGG.zonal(sum,cube;of=table)) == (2,3,3)
    burned=DGG.rasterize(sum,table;to=grid,fill=:value)
    @test burned[ip]==3 && burned[iq]==3
    coordtable=(lon=[10.0,10.0],lat=[25.0,25.0],v=[4,5])
    @test DGG.rasterize(sum,coordtable;to=grid,geometrycolumn=(:lon,:lat),fill=:v)[ip]==9
    feats=GI.FeatureCollection([GI.Feature(p;properties=(v=7,)),GI.Feature(q;properties=(v=8,))])
    @test DGG.rasterize(sum,feats;to=grid,fill=:v)[ip]==7
    st=DD.DimStack((a=cube,b=DD.DimArray(collect(1:DGG.ncells(grid)),DGG.Cells(DGG.CellLookup(grid)))))
    zst=DGG.zonal(sum,st;of=[p,q])
    @test size(zst[:a])==(2,3,2)
    @test size(zst[:b])==(2,)
    @test only(DGG.extract(st,[p];name=:b)).b==ip
end
end
