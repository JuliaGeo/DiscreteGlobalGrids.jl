# The write path: encoding choice, the coarse-ancestor chunk plan, the persisted
# manifest, and the dual-stamped attributes the store is read back through.
#
# Nothing here calls `dggread`. Every assertion reopens the store with Zarr
# directly and rebuilds the metadata snapshot out of plain dictionaries, so what
# is under test is the bytes on disk rather than a second helping of this
# package's own read code.

module DGGSIOWriteTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: IGeo7System, Z7Cell, CellVector, CellLookup, Cells,
    DGGSFormatError, levelgrid, descendants, ancestor, rawid, describe_store,
    dggread
import DimensionalData as DD

const ZARR_LOADED = try
    @eval using Zarr
    true
catch
    false
end

if !ZARR_LOADED
    @info "Zarr.jl is unavailable: the dggwrite suite is skipped."
else

    const SYS = IGeo7System()
    const LEVEL = 3

    # Three HEXAGON-rooted level-1 subtrees. A pentagon deletes one child digit,
    # so its subtree is short; these three hold exactly 7^2 level-3 cells each,
    # which is what makes the coarse-ancestor runs equal and the chunk plan's
    # alignment claim decidable.
    const ROOTS = [Z7Cell(DGG.idselect(levelgrid(SYS, 1), r)) for r in 1:3]
    const CELLS = sort!(reduce(vcat, [collect(descendants(SYS, c, LEVEL)) for c in ROOTS]))
    const NCELL = length(CELLS)
    const LOOKUP = CellLookup(CellVector(SYS, LEVEL, CELLS))
    const ELEV = Float32.(1:NCELL)
    const SLOPE = Float32.(0.5 .* (1:NCELL))

    # Named here rather than imported: the sidecar's array name and marker key
    # are the on-disk contract, and a test that reads them from the writer's own
    # constants would not notice one of them changing.
    const MANIFEST = "cell_chunk_manifest"
    const MARKER = "dggs_chunk_manifest"

    # The chunk plan's two fallbacks are properties of the plan itself, and one
    # of them needs a grid system no store spelling exists for, so that one is
    # asked of the planner directly.
    const WRITE = Base.get_extension(DGG, :DiscreteGlobalGridsZarrExt).DGGSZarrWrite

    # A downstream encoding that exists and does nothing else: it implements
    # none of the write verbs and is deliberately left out of the registry.
    struct Sketched <: DGG.CellEncoding end

    demostack() = DD.DimStack((elevation=copy(ELEV), slope=copy(SLOPE)), (Cells(LOOKUP),))
    dest(name) = joinpath(mktempdir(), name)

    # A source that remembers how it was read: how many times, and the largest
    # block any one read handed back. A chunk-at-a-time write asks once per
    # chunk and never for more than a chunk; materializing the layer first shows
    # up as one huge read, and copying it element by element as thousands.
    struct Recorder{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
        parent::A
        calls::Base.RefValue{Int}
        biggest::Base.RefValue{Int}
    end
    Recorder(p::AbstractArray) = Recorder(p, Ref(0), Ref(0))
    Base.size(r::Recorder) = size(r.parent)
    function Base.getindex(r::Recorder, I...)
        v = r.parent[I...]
        r.calls[] += 1
        r.biggest[] = max(r.biggest[], length(v))
        return v
    end

    # The two stages that own the cell ids, behind one function so `@allocated`
    # counts their work and nothing around it.
    function idbytes(A, grid)
        _, _, cells = WRITE._cellaxis(A)
        return WRITE._coordinate(DGG.DenseEncoding(), grid, cells, :rank)
    end

    # The `(n, 2)` arrays are written so the STORE's shape is `(n, 2)`; Zarr.jl
    # reverses between JSON and Julia, so what comes back here is `(2, n)`.
    rows(z) = permutedims(z[:, :])

    # A `StoreSnapshot` built from nothing but the reopened store's plain
    # dictionaries: the read path's own snapshot builder is a parallel task and
    # this suite must not depend on it.
    function handsnapshot(path)
        g = Zarr.zopen(path)
        arrays = DGG.ArrayEntry[]
        for name in sort!(collect(keys(g.arrays)))
            z = g.arrays[name]
            attrs = Dict{String,Any}(String(k) => v for (k, v) in z.attrs)
            dims = String[String(d) for d in get(attrs, "_ARRAY_DIMENSIONS", [])]
            push!(arrays, DGG.ArrayEntry(name=name, attrs=attrs,
                shape=reverse(size(z)), eltype=eltype(z), dims=dims))
        end
        return DGG.StoreSnapshot(identifier=path,
            attrs=Dict{String,Any}(String(k) => v for (k, v) in g.attrs),
            arrays=arrays)
    end

    described(enc, coordinate) = DGG.StoreDescription(gridname="igeo7",
        system=SYS, idscheme=:z7int, level=LEVEL, encoding=enc,
        coordinate=coordinate, spatial_dimension="cell_ids",
        variables=["elevation", "slope"])

    @testset "an eligible axis is written as ranges and reads back as itself" begin
        # The fixpoint: attrs -> description must return the description the
        # write path stamped. Kills a stamping mutant on either half of the dual
        # stamp -- a wrong `spatial_dimension`, a missing `compression`, a
        # coordinate name that does not name the array that was written.
        path = dest("ranges.zarr")
        @test DGG.dggwrite(path, demostack()) === path

        g = Zarr.zopen(path)
        @test sort!(collect(keys(g.arrays))) ==
              [MANIFEST, "cell_id_ranges", "elevation", "slope"]
        @test describe_store(handsnapshot(path)) ==
              described(DGG.RangesEncoding(), "cell_id_ranges")

        # Both stamps are present and each is readable on its own.
        @test g.attrs["dggs"]["compression"] == "ranges"
        @test g.attrs["dggs"]["spatial_dimension"] == "cell_ids"
        @test any(d -> d["uuid"] == DGG.ZARR_DGGS_UUID, g.attrs["zarr_conventions"])
        @test g["cell_id_ranges"].attrs["grid_name"] == "igeo7"
        @test g["cell_id_ranges"].attrs["level"] == LEVEL
    end

    @testset "the stored ranges and values are the cube that went in" begin
        # Kills a transposed range array, an off-by-one in the run scan, and any
        # reordering of the data against the axis.
        path = dest("values.zarr")
        DGG.dggwrite(path, demostack())
        g = Zarr.zopen(path)

        # Written under the default `:step` rule: a run breaks at every digit
        # rollover, so each of the 21 level-2 sibling sets in these three subtrees
        # is its own row. `merge = :rank` would make the same cells one row.
        @test size(g["cell_id_ranges"]) == (2, 21)   # store shape (21, 2)
        R = rows(g["cell_id_ranges"])
        grid = levelgrid(SYS, LEVEL)
        axis = DGG.cellaxis(DGG.RangesEncoding(), grid, R; declared_length=NCELL)
        @test collect(axis) == CELLS

        @test g["elevation"][:] == ELEV
        @test g["slope"][:] == SLOPE
        @test g["elevation"].attrs["_ARRAY_DIMENSIONS"] == ["cell_ids"]
    end

    @testset "the merge rule chooses what a run is, and the store shows it" begin
        # Ranks 200:220 of level 4 straddle two digit rollovers, so `:rank`
        # merges them into one interval that encloses ids naming no cell and
        # `:step` cannot. Kills a `merge` keyword that is dropped on the floor,
        # and a `:step` implementation that merges on rank anyway.
        grid = levelgrid(SYS, 4)
        ids = [DGG.idselect(grid, r) for r in 200:220]
        lk = CellLookup(CellVector(SYS, 4, Z7Cell.(ids)))
        A = DD.DimArray(Float32.(1:length(ids)), Cells(lk); name=:v)
        unit = DGG.idselect(grid, 1) - DGG.idselect(grid, 0)
        naive(a, b) = Int((b - a) ÷ unit) + 1

        prank, pstep = dest("rank.zarr"), dest("step.zarr")
        DGG.dggwrite(prank, A; merge=:rank)
        DGG.dggwrite(pstep, A; merge=:step)
        Rrank = rows(Zarr.zopen(prank)["cell_id_ranges"])
        Rstep = rows(Zarr.zopen(pstep)["cell_id_ranges"])

        @test size(Rrank, 1) < size(Rstep, 1)
        @test size(Rrank, 1) == 1
        # The defining property of each rule, not just the row count: a step run
        # holds every integer it spans, a rank run need not.
        @test all(i -> naive(Rstep[i, 1], Rstep[i, 2]) ==
                       DGG.idcount_between(grid, Rstep[i, 1], Rstep[i, 2]),
            1:size(Rstep, 1))
        @test naive(Rrank[1, 1], Rrank[1, 2]) !=
              DGG.idcount_between(grid, Rrank[1, 1], Rrank[1, 2])
        # And both are the same axis.
        for R in (Rrank, Rstep)
            @test collect(DGG.cellaxis(DGG.RangesEncoding(), grid, R;
                declared_length=length(ids))) == Z7Cell.(ids)
        end
    end

    @testset "encoding = :dense writes the ids themselves" begin
        # The interop escape. Kills an `encoding` keyword that `:auto` overrules.
        path = dest("dense.zarr")
        DGG.dggwrite(path, demostack(); encoding=:dense)
        g = Zarr.zopen(path)

        @test g["cell_ids"][:] == rawid.(CELLS)
        @test g["cell_ids"].attrs["_ARRAY_DIMENSIONS"] == ["cell_ids"]
        @test g.attrs["dggs"]["compression"] == "none"
        @test describe_store(handsnapshot(path)) ==
              described(DGG.DenseEncoding(), "cell_ids")
    end

    @testset "a cell axis that is not canonical is refused" begin
        # A reversed cube keeps its cell ids and loses their order, and a store
        # written that way could never be opened again -- the read path verifies
        # sortedness at open. Kills a writer that stores whatever it is handed.
        A = DD.DimArray(ELEV, Cells(LOOKUP); name=:elevation)
        reversed = reverse(A; dims=Cells)
        @test eltype(DD.val(DD.dims(reversed, Cells))) <: DGG.AbstractCellIndex

        @test_throws DGGSFormatError DGG.dggwrite(dest("rev1.zarr"), reversed)
        @test_throws DGGSFormatError DGG.dggwrite(dest("rev2.zarr"), reversed;
            encoding=:ranges)
        @test_throws DGGSFormatError DGG.dggwrite(dest("rev3.zarr"), reversed;
            encoding=:dense)
    end

    @testset "a mixed-level axis is refused by name" begin
        # No registered encoding writes mixed levels. Kills a writer that
        # routes a `MultiOrderLookup` into the unsorted-axis refusal, or one
        # that writes it as if it were single-level.
        mov = DGG.MultiOrderVector(SYS, [ROOTS[1]; collect(DGG.children(SYS, ROOTS[2]))])
        M = DD.DimArray(Float32.(1:length(mov)), Cells(DGG.MultiOrderLookup(mov));
            name=:elevation)
        err = try
            DGG.dggwrite(dest("moc.zarr"), M)
            nothing
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :mixed_level_axis
    end

    @testset "the auto chunk plan breaks on coarse-ancestor boundaries" begin
        # 49 level-3 cells per level-1 subtree, so two whole subtrees is the
        # largest whole number of runs under a 100-cell target. Kills a plan
        # that takes the target itself (100), which would split the third
        # subtree across both chunks.
        path = dest("chunked.zarr")
        DGG.dggwrite(path, demostack(); chunk_target=100)
        g = Zarr.zopen(path)
        @test g["elevation"].metadata.chunks == (98,)

        m = rows(g[MANIFEST])
        @test size(m, 1) == 2
        coarse(x) = ancestor(SYS, Z7Cell(x), 1)
        @test coarse(m[1, 2]) != coarse(m[2, 1])

        # And the manifest describes the chunks that were actually written.
        for c in 1:size(m, 1)
            lo, hi = (c - 1) * 98 + 1, min(c * 98, NCELL)
            @test m[c, 1] == rawid(CELLS[lo])
            @test m[c, 2] == rawid(CELLS[hi])
        end
        # The plan's claim, recorded where an aggregation can read it: these
        # chunks really are whole level-1 subtrees.
        marker = g[MANIFEST].attrs[MARKER]
        @test marker["chunk_length"] == 98
        @test marker["length"] == NCELL
        @test marker["ancestor_level"] == 1
        @test marker["ancestor_aligned"]
        # What the manifest was validated AGAINST, and not only how it is
        # shaped: a reader that skips the scan on this marker is reading the
        # ids at the level and in the id arithmetic named here, so a store
        # whose attributes were edited afterwards no longer matches its own
        # sidecar. Kills a marker that says only how big the chunks are.
        @test marker["level"] == LEVEL
        @test marker["grid"] == "igeo7"
    end

    @testset "partial coarse subtrees get the plan, not the guarantee" begin
        # Unequal ancestor runs. A uniform chunk length cannot land on all of
        # them -- Zarr has no other kind -- so `:auto` takes the largest whole
        # number of runs under the target, says so, and lets the manifest be the
        # commitment. Kills a writer that claims alignment it does not have, and
        # a manifest that is computed from the plan rather than from the axis.
        g1 = levelgrid(SYS, 1)
        # Five hexagons of base cell 0 and one of base cell 1; rank 0 and rank 6
        # are the pentagons, whose subtrees are short.
        roots = [Z7Cell(DGG.idselect(g1, r)) for r in (1, 2, 3, 4, 5, 7)]
        counts = [20, 30, 10, 20, 30, 10]
        cells = reduce(vcat, [collect(descendants(SYS, c, LEVEL))[1:k]
                              for (c, k) in zip(roots, counts)])
        @test issorted(cells) && length(cells) == 120

        lk = CellLookup(CellVector(SYS, LEVEL, cells))
        path = dest("partial.zarr")
        DGG.dggwrite(path, DD.DimArray(Float32.(1:120), Cells(lk); name=:v);
            chunk_target=55)
        g = Zarr.zopen(path)

        # 20 + 30 is the largest whole number of runs under 55; the boundary at
        # 100 then falls inside the fifth run.
        @test g["v"].metadata.chunks == (50,)
        marker = g[MANIFEST].attrs[MARKER]
        @test marker["ancestor_level"] == 1
        @test marker["ancestor_aligned"] == false

        m = rows(g[MANIFEST])
        @test size(m, 1) == 3
        for c in 1:3
            @test m[c, 1] == rawid(cells[(c-1)*50+1])
            @test m[c, 2] == rawid(cells[min(c*50, 120)])
        end
    end

    @testset "a small axis and an integer chunk size" begin
        # The default million-cell target degrades to a single chunk on a test
        # axis; an integer overrides the plan outright. Kills a `chunks` keyword
        # that only ever means `:auto`.
        one = dest("one.zarr")
        DGG.dggwrite(one, demostack())
        @test Zarr.zopen(one)["elevation"].metadata.chunks == (NCELL,)

        fixed = dest("fixed.zarr")
        DGG.dggwrite(fixed, demostack(); chunks=50)
        g = Zarr.zopen(fixed)
        @test g["elevation"].metadata.chunks == (50,)
        m = rows(g[MANIFEST])
        @test size(m, 1) == 3
        for c in 1:3
            @test m[c, 1] == rawid(CELLS[(c-1)*50+1])
            @test m[c, 2] == rawid(CELLS[min(c*50, NCELL)])
        end
        @test g[MANIFEST].attrs["_ARRAY_DIMENSIONS"] == ["chunks", "bounds"]
    end

    @testset "the store carries consolidated metadata" begin
        # Kills a writer that leaves `.zmetadata` off, which costs a cloud
        # reader one request per array.
        path = dest("consolidated.zarr")
        DGG.dggwrite(path, demostack())
        @test isfile(joinpath(path, ".zmetadata"))
        g = Zarr.zopen(path; consolidated=true)
        @test sort!(collect(keys(g.arrays))) ==
              [MANIFEST, "cell_id_ranges", "elevation", "slope"]
        @test g.attrs["dggs"]["refinement_level"] == LEVEL
    end

    @testset "a non-cell dimension becomes an ordinary Zarr dimension" begin
        # Cells are the fastest-varying axis on disk whichever way the cube is
        # laid out in Julia. Kills a writer that stamps `_ARRAY_DIMENSIONS` in
        # Julia order, or that writes the values without permuting to match.
        times = 1:4
        data = Float32[10i + t for i in 1:NCELL, t in times]
        A = DD.DimArray(data, (Cells(LOOKUP), DD.Dim{:time}(times)); name=:temperature)

        path = dest("time.zarr")
        DGG.dggwrite(path, A)
        g = Zarr.zopen(path)
        @test g["temperature"].attrs["_ARRAY_DIMENSIONS"] == ["time", "cell_ids"]
        @test size(g["temperature"]) == (NCELL, 4)
        @test g["temperature"][:, :] == data
        @test g["time"][:] == collect(times)
        @test describe_store(handsnapshot(path)).variables == ["temperature"]

        # The transposed cube is the same store.
        tpath = dest("time_t.zarr")
        DGG.dggwrite(tpath, permutedims(A))
        h = Zarr.zopen(tpath)
        @test h["temperature"].attrs["_ARRAY_DIMENSIONS"] == ["time", "cell_ids"]
        @test h["temperature"][:, :] == data
    end

    @testset "an already-open group is a destination too" begin
        # The same pipeline, one step further along: the group exists, so its
        # attributes are stamped rather than passed to its constructor. Kills a
        # writer that only ever stamps at creation time.
        dir = dest("group.zarr")
        g = Zarr.zgroup(dir)
        @test DGG.dggwrite(g, demostack()) === g
        @test describe_store(handsnapshot(dir)) ==
              described(DGG.RangesEncoding(), "cell_id_ranges")
    end

    @testset "the chunk plan falls back where no coarse level helps" begin
        # Both fallbacks the `ChunkPlan` docstring commits to, neither of which
        # any other fixture reaches. Kills a planner that always claims an
        # ancestor level.
        #
        # One: the level-2 runs are seven cells each, already past a five-cell
        # target, so the descent never starts and the target itself is the
        # chunk length.
        path = dest("toosmall.zarr")
        DGG.dggwrite(path, demostack(); chunk_target=5)
        g = Zarr.zopen(path)
        @test g["elevation"].metadata.chunks == (5,)
        @test !haskey(g[MANIFEST].attrs[MARKER], "ancestor_level")

        # Two: a system whose subtrees are not contiguous in canonical order has
        # no ancestor runs to group by at all.
        a5 = DGG.A5System()
        @test !DGG.has_sorted_subtrees(a5)
        grid = levelgrid(a5, 2)
        plan = WRITE._chunkplan(:auto, grid, [DGG.cellindex(grid, i) for i in 1:60], 20)
        @test plan.ancestor_level === nothing
        @test plan.chunklength == 20
    end

    @testset "chunk_target counts elements, not cells" begin
        # A forty-step time axis makes every cell forty values, so a chunk of a
        # million cells is a forty-million-element chunk. The target bounds the
        # elements: 400 ÷ 40 is a ten-cell target, and the largest whole number
        # of level-2 subtree runs under ten is one run of seven. Kills a plan
        # that measures a chunk in cells and lets the other dimensions span.
        steps = 40
        data = Float32[10i + t for i in 1:NCELL, t in 1:steps]
        A = DD.DimArray(data, (Cells(LOOKUP), DD.Dim{:time}(1:steps)); name=:temperature)
        path = dest("elements.zarr")
        DGG.dggwrite(path, A; chunk_target=400)
        g = Zarr.zopen(path)
        @test prod(g["temperature"].metadata.chunks) <= 400
        @test g["temperature"].metadata.chunks == (7, steps)
        @test g[MANIFEST].attrs[MARKER]["chunk_length"] == 7
    end

    @testset "a chunk size of zero is refused in both spellings" begin
        # Zero cells per chunk is one chunk per cell, which on the axis this is
        # built for is tens of millions of files. Kills an unvalidated
        # `chunk_target` next to a validated `chunks`.
        @test_throws ArgumentError DGG.dggwrite(dest("z1.zarr"), demostack(); chunks=0)
        @test_throws ArgumentError DGG.dggwrite(dest("z2.zarr"), demostack(); chunk_target=0)
        @test_throws ArgumentError DGG.dggwrite(dest("z3.zarr"), demostack(); chunk_target=-5)
    end

    @testset "an unknown merge rule is refused whatever the encoding" begin
        # `merge` only reaches the ranges coordinate, but a caller who misspells
        # it while writing dense has still asked for a rule that does not exist.
        # Kills validation that lives inside `idranges` alone.
        for enc in (:auto, :dense, :ranges)
            @test_throws ArgumentError DGG.dggwrite(dest("merge_$enc.zarr"),
                demostack(); encoding=enc, merge=:middle)
        end
    end

    @testset "the producer's attributes ride through to disk" begin
        # Units and long names are the producer's; `_ARRAY_DIMENSIONS` is this
        # writer's and wins the collision, because a stale one from another
        # layout would describe the array wrongly. Kills a layer plan that
        # builds its attributes from scratch.
        elev = DD.DimArray(copy(ELEV), (Cells(LOOKUP),); name=:elevation,
            metadata=Dict{String,Any}("units" => "m", "long_name" => "elevation",
                "_ARRAY_DIMENSIONS" => ["stale"]))
        slope = DD.DimArray(copy(SLOPE), (Cells(LOOKUP),); name=:slope)
        st = DD.DimStack((elevation=elev, slope=slope);
            metadata=Dict{String,Any}("attrs" =>
                Dict{String,Any}("title" => "Pori in miniature")))

        path = dest("attrs.zarr")
        DGG.dggwrite(path, st)
        g = Zarr.zopen(path)
        @test g["elevation"].attrs["units"] == "m"
        @test g["elevation"].attrs["long_name"] == "elevation"
        @test g["elevation"].attrs["_ARRAY_DIMENSIONS"] == ["cell_ids"]
        @test g.attrs["title"] == "Pori in miniature"
        # And the convention still won its own keys.
        @test g.attrs["dggs"]["refinement_level"] == LEVEL
    end

    @testset "a store rewritten from its own read is the same store" begin
        # The design's fixpoint, both encodings, more than one chunk each: the
        # values, the axis, the layer attributes and the group attributes all
        # survive dggread -> dggwrite -> dggread. Kills an encode/decode
        # asymmetry and every attribute the write path drops on the floor.
        #
        # The `:dense` arm reads both stores back through the TRUSTED path —
        # `dggwrite` persists a manifest and `dggread` believes it — so the ids
        # are not rescanned here. They were scanned before either store was
        # committed, by the write path's own `cellaxis` call, which is what the
        # marker's `validated = "strict"` records.
        for enc in (:dense, :ranges)
            elev = DD.DimArray(copy(ELEV), (Cells(LOOKUP),); name=:elevation,
                metadata=Dict{String,Any}("units" => "m", "long_name" => "elevation"))
            slope = DD.DimArray(copy(SLOPE), (Cells(LOOKUP),); name=:slope,
                metadata=Dict{String,Any}("units" => "m/m"))
            # Declared out of order on purpose: the store keeps one order.
            src = DD.DimStack((slope=slope, elevation=elev);
                metadata=Dict{String,Any}("attrs" =>
                    Dict{String,Any}("title" => "twin", "institution" => "DGG.jl")))

            once = dest("fix1_$enc.zarr")
            DGG.dggwrite(once, src; encoding=enc, chunks=50)
            r1 = dggread(once)
            twice = dest("fix2_$enc.zarr")
            DGG.dggwrite(twice, r1; encoding=enc, chunks=50)
            r2 = dggread(twice)

            # Layers come back alphabetical, whichever order went in.
            @test collect(keys(r1)) == [:elevation, :slope]
            @test collect(keys(r2)) == collect(keys(r1))

            @test collect(DD.lookup(r2[:elevation], Cells)) == CELLS
            @test DD.lookup(r2[:elevation], Cells) == DD.lookup(r1[:elevation], Cells)
            @test collect(parent(r2[:elevation])) == ELEV
            @test collect(parent(r2[:slope])) == SLOPE

            for st in (r1, r2)
                @test DD.metadata(st[:elevation])["units"] == "m"
                @test DD.metadata(st[:elevation])["long_name"] == "elevation"
                @test DD.metadata(st[:slope])["units"] == "m/m"
                # The one key the writer regenerates, and it rides back out.
                @test DD.metadata(st[:elevation])["_ARRAY_DIMENSIONS"] == ["cell_ids"]
                @test DD.metadata(st)["attrs"]["title"] == "twin"
                @test DD.metadata(st)["attrs"]["institution"] == "DGG.jl"
            end
            # A fixpoint from the first write onward: whatever normalization the
            # conventions apply has already happened by then.
            @test DD.metadata(r2)["attrs"] == DD.metadata(r1)["attrs"]
            @test DD.metadata(r2)["encoding"] == DD.metadata(r1)["encoding"]
            @test DD.metadata(r2)["description"] == DD.metadata(r1)["description"]
        end
    end

    @testset "a lazy source is copied one chunk at a time" begin
        # The premise is tens of millions of cells, where materializing a layer
        # to write it costs the whole axis in memory. Every read of the source
        # is counted, so a write that pulls the axis in one piece -- or that
        # falls back to copying it element by element -- is the wrong count.
        src = dest("lazysrc.zarr")
        DGG.dggwrite(src, demostack(); chunks=50)
        st = dggread(src)
        @test parent(st[:elevation]) isa Zarr.ZArray

        ks = Tuple(keys(st))
        recs = map(k -> Recorder(parent(st[k])), ks)
        wrapped = DD.DimStack(NamedTuple{ks}(map((k, r) -> DD.rebuild(st[k], r), ks, recs)))

        out = dest("lazyout.zarr")
        DGG.dggwrite(out, wrapped; chunks=50)
        nchunks = cld(NCELL, 50)
        for r in recs
            @test r.calls[] == nchunks
            @test r.biggest[] <= 50
        end

        # And the bytes are still the cube.
        g = Zarr.zopen(out)
        @test g["elevation"][:] == ELEV
        @test g["slope"][:] == SLOPE
    end

    @testset "the cell ids are materialized once, and they are the bytes on disk" begin
        # The axis of a cube on its way to disk used to be copied twice: once out
        # of the lookup as typed cells, once again as the raw ids the coordinate
        # writes. That is 160 MB at ten million cells, for a write whose premise
        # is that it never holds the axis twice. Kills the second copy: the ids
        # come out of the lookup ONCE, raw, and every stage after that passes the
        # same array along.
        n = 50_000
        grid = levelgrid(SYS, 6)
        lk = CellLookup(CellVector(SYS, 6, [Z7Cell(DGG.idselect(grid, r)) for r in 0:n-1]))
        A = DD.DimArray(Float32.(1:n), Cells(lk); name=:v)

        celldim, axisgrid, cells = WRITE._cellaxis(A)
        @test DD.name(celldim) === :Cells
        @test (DGG.system(axisgrid), DGG.level(axisgrid)) == (SYS, 6)
        # Raw ids, not typed cells: one vector, of what the store holds.
        @test cells isa Vector{UInt64}
        @test cells == DGG.rawid.(collect(lk))
        # Identity, not equality: the dense coordinate IS that vector, so nothing
        # between the lookup and the bytes copies it again.
        @test WRITE._coordinate(DGG.DenseEncoding(), grid, cells, :rank) === cells

        onecopy = 8 * n
        idbytes(A, grid)        # compile before measuring
        @test onecopy <= @allocated(idbytes(A, grid)) < 3 * onecopy ÷ 2
    end

    @testset "a group that already holds these arrays is refused, not half-stamped" begin
        # `zcreate` throws on the first name it cannot take, and by then the
        # group's attributes have been merged: the store would claim one
        # encoding and hold another, and would not open again. Kills a writer
        # that stamps before it knows the names are free.
        dir = dest("occupied.zarr")
        DGG.dggwrite(dir, demostack())      # ranges
        before = dggread(dir)

        g = Zarr.zopen(dir, "w")
        err = try
            DGG.dggwrite(g, demostack(); encoding=:dense)
        catch e
            e
        end
        @test err isa DGGSFormatError
        @test occursin("elevation", sprint(showerror, err))
        # And it reports the group by the path a person would type, not by the
        # `show` of whatever store object is behind it.
        @test occursin("occupied.zarr", err.store)
        @test !occursin("DirectoryStore", err.store)

        # Untouched: the same encoding, the same axis, the same values.
        after = dggread(dir)
        @test DD.metadata(before)["encoding"] == "ranges"
        @test DD.metadata(after)["encoding"] == "ranges"
        @test collect(DD.lookup(after[:elevation], Cells)) == CELLS
        @test collect(parent(after[:elevation])) == ELEV
    end

    @testset "an encoding with no write path is named, not a MethodError" begin
        # The write pipeline is a set of private verbs on the encoding, and a
        # downstream encoding that registers itself without implementing them
        # used to fall off the end of dispatch — a MethodError on a private name
        # that says nothing about what is missing. The read half already answers
        # this case by name (`storedaxis`), and this is the same answer on the
        # other side. Kills a write path with no fallback.
        err = try
            DGG.dggwrite(dest("noencoding.zarr"), demostack(); encoding=Sketched())
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :unsupported_encoding
        msg = sprint(showerror, err)
        @test occursin("Sketched", msg)
        # And what a downstream encoding has to do about it: the layouts that
        # are implemented are named, so the message is actionable.
        @test occursin("ranges", msg) && occursin("none", msg)
    end

    @testset "a layer that takes one of the writer's own array names is refused" begin
        # The cell coordinate, the range array and the manifest have fixed names,
        # and a layer called `cell_ids` collides with one of them. Refused before
        # anything is created: the alternative is `zcreate` failing partway
        # through and leaving a directory that holds some of a store. Kills a
        # writer that discovers the collision from Zarr.
        A = DD.DimArray(copy(ELEV), (Cells(LOOKUP),); name=:cell_ids)
        for enc in (:auto, :dense)
            path = dest("reserved_$enc.zarr")
            err = try
                DGG.dggwrite(path, A; encoding=enc)
            catch e
                e
            end
            @test err isa DGGSFormatError && err.check === :reserved_array_name
            msg = sprint(showerror, err)
            @test occursin("cell_ids", msg) && occursin(MANIFEST, msg)
            @test !isdir(path)
        end

        # The same check across the layers themselves: a `time` dimension writes
        # its own coordinate array, and a layer of that name would be a second
        # array with one name.
        B = DD.DimArray(Float32[10i + t for i in 1:NCELL, t in 1:3],
            (Cells(LOOKUP), DD.Dim{:time}(1:3)); name=:time)
        err = try
            DGG.dggwrite(dest("dupname.zarr"), B)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :duplicate_array_name
        @test occursin("time", sprint(showerror, err))
    end

    @testset "a remote destination is refused rather than half-written" begin
        # Kills a writer that hands `gs://bucket/x.zarr` to `DirectoryStore` and
        # silently creates a local directory called `gs:`.
        for url in ("gs://bucket/x.zarr", "s3://bucket/x.zarr",
            "https://example.com/x.zarr")
            @test_throws ArgumentError DGG.dggwrite(url, demostack())
        end
        @test !isdir("gs:")
    end

end # ZARR_LOADED

end # module DGGSIOWriteTests
