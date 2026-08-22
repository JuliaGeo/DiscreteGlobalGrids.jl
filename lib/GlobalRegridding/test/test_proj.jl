# Task-local Proj charts for projected RasterGrid coordinates.

import Proj

const GRProjExt = Base.get_extension(GlobalRegridding, :GlobalRegriddingProjExt)

tuple_isapprox(a, b) = all(isapprox(x, y) for (x, y) in zip(a, b))

function projected_snapshot(space, querycap)
    centre = cellcentroid(space, 5)
    ring = Tuple.(GI.getpoint(GI.getexterior(getcell(space, 5))))
    native = GR.chartcoords(space, centre)
    located = cellat(space, centre)
    candidates = GR.candidatechunks!(Int[], GR.chunkindex(space), querycap)
    return (; centre = Tuple(centre), ring, native, located, candidates)
end

@testset "task-local Proj raster transformations" begin
    @test GRProjExt !== nothing

    # Raster templates must use X/Y (easting/northing and longitude/latitude)
    # order. The selected, axis-normalized pipeline is cloned for each task.
    template = Proj.Transformation("EPSG:3857", "EPSG:4326"; always_xy = true)
    data = CountingChunked(zeros(3, 3), (1, 1))
    raster = DD.DimArray(data,
        (DD.X(-1_000_000.0:1_000_000.0:1_000_000.0),
         DD.Y(4_000_000.0:1_000_000.0:6_000_000.0)))
    space = RasterGrid(raster; native_to_unit_sphere = template)

    @test data.reads == 0
    @test hascellchart(space)
    @test space.tables === nothing
    @test space.xperiod === nothing
    @test space.native_to_unit_sphere.state === space.unit_sphere_to_native.state
    state = space.native_to_unit_sphere.state
    @test state.template === template
    @test state.direction == template.direction == Proj.PJ_FWD

    # Construction-time validation and winding already created this task's
    # pair. Later geometry reuses it, and neither clone is the shared template.
    pair = GRProjExt._task_pair(state)
    @test GRProjExt._task_pair(state) === pair
    @test task_local_storage()[state.task_key] === pair
    @test pair.ctx != C_NULL
    @test pair.forward.pj != C_NULL
    @test pair.reverse.pj != C_NULL
    @test pair.forward.pj != template.pj
    @test pair.reverse.pj != template.pj
    @test pair.forward.pj != pair.reverse.pj
    @test pair.forward.direction == state.direction
    @test pair.reverse.direction == inv(state.direction)
    prepared_forward = GR._task_prepared_raster_transform(space.native_to_unit_sphere)
    @test prepared_forward.transformation === pair.forward

    to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
    xaxis, yaxis = GR.chartaxes(space)
    x, y = first(xaxis), first(yaxis)
    expected = to_sphere(template((x, y)))
    @test tuple_isapprox(Tuple(cellcentroid(space, 1)), Tuple(expected))
    @test tuple_isapprox(GR.chartcoords(space, expected), (x, y))
    @test all(cellat(space, cellcentroid(space, i)) == i for i in 1:ncells(space))
    @test data.reads == 0

    querycap = SphericalCap(cellcentroid(space, 5), 0.5)
    reference = projected_snapshot(space, querycap)
    jobs = [Threads.@spawn begin
        owner = GRProjExt._task_pair(state)
        (; owner, result = projected_snapshot(space, querycap))
    end for _ in 1:4]
    task_results = fetch.(jobs)

    # Isolation follows Julia task identity, independently of which OS thread
    # happens to run a task.
    owners = getproperty.(task_results, :owner)
    @test length(unique(objectid.(owners))) == length(owners)
    @test length(unique(UInt(owner.ctx) for owner in owners)) == length(owners)
    @test length(unique(UInt(owner.forward.pj) for owner in owners)) == length(owners)
    @test length(unique(UInt(owner.reverse.pj) for owner in owners)) == length(owners)
    for task_result in task_results
        result = task_result.result
        @test tuple_isapprox(result.centre, reference.centre)
        @test all(tuple_isapprox(a, b) for (a, b) in zip(result.ring, reference.ring))
        @test tuple_isapprox(result.native, reference.native)
        @test result.located == reference.located
        @test result.candidates == reference.candidates
    end
    @test data.reads == 0

    # Explicit close uses the same idempotent path as the owner finalizer and
    # makes the task cache recreate a complete pair on its next use.
    close(pair)
    @test pair.closed
    @test pair.ctx == C_NULL
    @test pair.forward === nothing
    @test pair.reverse === nothing
    close(pair)
    replacement = GRProjExt._task_pair(state)
    @test replacement !== pair
    @test replacement.ctx != C_NULL
    @test cellat(space, cellcentroid(space, 5)) == 5

    # Exercise the registered finalizer itself, not only the shared close
    # implementation. Repeated finalization is harmless and TLS recovers.
    finalize(replacement)
    @test replacement.closed
    @test replacement.ctx == C_NULL
    @test replacement.forward === nothing
    @test replacement.reverse === nothing
    finalize(replacement)
    after_finalize = GRProjExt._task_pair(state)
    @test after_finalize !== replacement
    @test after_finalize.ctx != C_NULL

    # If the second clone fails, the already-owned first transformation is
    # finalized before the context. A failed construction never enters TLS.
    failure_template = Proj.Transformation(
        "EPSG:3857", "EPSG:4326"; always_xy = true)
    failure_state = GRProjExt._ProjChartState(failure_template)
    events = Symbol[]
    clone_calls = Ref(0)
    clone_object = function (pj, ctx)
        clone_calls[] += 1
        clone_calls[] == 1 && return Proj.proj_clone(pj, ctx)
        return Ptr{Proj.PJ}(C_NULL)
    end
    release_transformation = function (transformation)
        push!(events, :transformation)
        finalize(transformation)
        @test transformation.pj == C_NULL
    end
    destroy_context = function (ctx)
        push!(events, :context)
        destroyed = Proj.proj_context_destroy(ctx)
        @test destroyed == C_NULL
        return destroyed
    end
    @test_throws ErrorException GRProjExt._task_pair(failure_state;
        clone_object, release_transformation, destroy_context)
    @test events == [:transformation, :context]
    @test !haskey(task_local_storage(), failure_state.task_key)

    # A non-null clone whose Julia wrapper fails is still released through
    # Proj.jl before its context is destroyed.
    raw_template = Proj.Transformation("EPSG:3857", "EPSG:4326"; always_xy = true)
    raw_state = GRProjExt._ProjChartState(raw_template)
    raw_events = Symbol[]
    wrap_transformation = (ptr, direction) -> error("injected wrapper failure")
    destroy_raw = function (ptr)
        push!(raw_events, :raw)
        return Proj.proj_destroy(ptr)
    end
    raw_destroy_context = function (ctx)
        push!(raw_events, :context)
        return Proj.proj_context_destroy(ctx)
    end
    @test_throws ErrorException GRProjExt._new_task_pair(raw_state;
        wrap_transformation, destroy_raw, destroy_context = raw_destroy_context)
    @test raw_events == [:raw, :context]

    finalized_template = Proj.Transformation(
        "EPSG:3857", "EPSG:4326"; always_xy = true)
    finalize(finalized_template)
    @test finalized_template.pj == C_NULL
    @test_throws ArgumentError RasterGrid(raster;
        native_to_unit_sphere = finalized_template)

    # This package delegates every native operation to Proj.jl's Julia API.
    package_root = dirname(dirname(pathof(GlobalRegridding)))
    for subdir in ("src", "ext"), (root, _, files) in walkdir(joinpath(package_root, subdir))
        for file in files
            endswith(file, ".jl") || continue
            source = read(joinpath(root, file), String)
            @test !occursin(r"(?m)^\s*(?:@ccall|ccall)\b", source)
            @test !occursin("PROJ_jll", source)
        end
    end
end
