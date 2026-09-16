module RasterRastersExtTests
using Test, Statistics
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import Rasters as RA
import GeoInterface as GI

@testset "Rasters wrappers and sentinels" begin
    grid = DGG.levelgrid(DGG.HEALPixSystem(),1)
    p = GI.Point((10.0,25.0))
    i = only(DGG._raster_indices(grid,p))
    A = RA.Raster(fill(-999.0,DGG.ncells(grid)),(DGG.Cells(DGG.CellLookup(grid)),);
        missingval=-999.0,name=:temperature,metadata=Dict("units"=>"C"))
    A[i]=2.0
    B=DGG.rasterize(sum,[p,p];to=A,fill=[3.0,4.0])
    @test B isa RA.AbstractRaster
    @test RA.missingval(B)==-999.0
    @test B[i]==7.0
    @test RA.metadata(B)==RA.metadata(A)
    @test DD.name(B)==DD.name(A)
    @test count(==(-999.0),B)==length(B)-1
    @test DGG.zonal(mean,A;of=[p],emptyval=NaN)[1]==2.0
    @test DGG.mask(A;with=[p]) isa RA.AbstractRaster
    @test RA.missingval(DGG.boolmask([p];to=A))===nothing
    st=RA.RasterStack((a=A,b=RA.rebuild(A;name=:b));metadata=Dict("title"=>"stack"))
    @test DGG.rasterize(sum,[p];to=st,fill=(a=1,b=2)) isa RA.AbstractRasterStack
    @test DGG.zonal(mean,st;of=[p],emptyval=NaN) isa RA.AbstractRasterStack
    @test DD.metadata(DGG.zonal(mean,st;of=[p],emptyval=NaN))==DD.metadata(st)
    @test length(DGG.extract(st,[p];skipmissing=true))==1
    A[i]=-999.0
    @test isempty(DGG.extract(A,[p];skipmissing=true))
    @test DGG.zonal(sum,A;of=p)==0
    @test DGG.zonal(mean,A;of=p,emptyval=NaN) |> isnan
end
end
