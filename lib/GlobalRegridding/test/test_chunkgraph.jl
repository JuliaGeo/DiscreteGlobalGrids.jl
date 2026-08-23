import Graphs
import DimensionalData as DD
import Serialization

# The four relations — truth, demand, cap join, and the graph's own rows — are
# defined once, in `graphoracles.jl`, and shared with the root suite and the G1
# harness. Do not re-spell any of them here.
include(joinpath(@__DIR__, "graphoracles.jl"))
using .ChunkGraphOracles: contributing_pairs, graph_pairs, demanded_pairs,
    capjoin_pairs, eager_pairs

# A raster whose chunk grid is deliberately out of step with the quadtree's
# power-of-two splits, so its leaves straddle chunk boundaries.
function misalignedraster(nx, ny, cx, cy)
    dlon, dlat = 360 / nx, 180 / ny
    array = DD.DimArray(zeros(nx, ny),
        (DD.X(range(-180 + dlon / 2, 180 - dlon / 2; length = nx)),
         DD.Y(range(-90 + dlat / 2, 90 - dlat / 2; length = ny))))
    return RasterGrid(array;
        chunks = ([lo:min(lo + cx - 1, nx) for lo in 1:cx:nx],
            [lo:min(lo + cy - 1, ny) for lo in 1:cy:ny]))
end

# Task G4 made a narrow phase an argument to plan construction and to nothing
# else, so every test below that needs one goes through a plan. A plan takes its
# radius from its method rather than from a keyword, so this method supplies the
# radius a test wants. It builds no weights and is never asked to.
struct G4RadiusMethod <: AbstractRegriddingMethod
    radius::Float64
end
GR.support_radius(m::G4RadiusMethod, ::RegridSpace) = m.radius

planned_dependencies(dst, src; radius = 0.0, kw...) =
    GR.dependencies(ChunkedPlan(G4RadiusMethod(radius), Weighted(0.5), dst, src;
        dependencies = true, kw...))

# A space that counts every question the relation machinery asks it.
# `chunkextents` is that question: the destination caps a relation is built by
# querying with, the caps the identity stamp hashes, and the vector the generic
# `chunkindex` packs all come through it, so this counter is nonzero exactly
# when a relation is being derived from this space and zero when nothing is.
#
# Until Task E2 the counter sat on `chunktree`, which `chunkextents` collected
# from. Removing that bridge moved the funnel up one level; it did not widen or
# narrow it, because `chunktree` had no caller of its own.
mutable struct G4ProbeSpace{S<:RegridSpace} <: RegridSpace
    space::S
    queries::Int
end
G4ProbeSpace(space::RegridSpace) = G4ProbeSpace(space, 0)
GR.chunkextents(p::G4ProbeSpace) = (p.queries += 1; GR.chunkextents(p.space))

# A source chunk index that counts the queries a builder makes of it. The
# Phase 4 gate is that ONE query implementation defines every edge, so the
# counter and the pairs together say what the builder did: how many questions it
# asked, and that the answers are the relation.
mutable struct E2QueryIndex{I}
    inner::I
    queries::Int
end

GR.candidatechunks!(out::Vector{Int}, ix::E2QueryIndex, dstcap::GR.Cap;
    radius::Real = 0.0) =
    (ix.queries += 1; GR.candidatechunks!(out, ix.inner, dstcap; radius))

struct E2QuerySpace{S<:RegridSpace,I} <: RegridSpace
    space::S
    index::E2QueryIndex{I}
end

E2QuerySpace(space::RegridSpace) =
    E2QuerySpace(space, E2QueryIndex(GR.chunkindex(space), 0))

GR.chunkindex(s::E2QuerySpace) = s.index
GR.chunkextents(s::E2QuerySpace) = GR.chunkextents(s.space)
ncells(s::E2QuerySpace) = ncells(s.space)
getcell(s::E2QuerySpace, i::Int) = getcell(s.space, i)
manifold(s::E2QuerySpace) = manifold(s.space)
hascellchart(s::E2QuerySpace) = hascellchart(s.space)
cellcentroid(s::E2QuerySpace, i::Int) = cellcentroid(s.space, i)
cellat(s::E2QuerySpace, x) = cellat(s.space, x)
nchunks(s::E2QuerySpace) = nchunks(s.space)
cellindices(s::E2QuerySpace, c::Int) = cellindices(s.space, c)
celltree(s::E2QuerySpace) = celltree(s.space)
ncells(p::G4ProbeSpace) = ncells(p.space)
getcell(p::G4ProbeSpace, i::Int) = getcell(p.space, i)
manifold(p::G4ProbeSpace) = manifold(p.space)
hascellchart(p::G4ProbeSpace) = hascellchart(p.space)
cellcentroid(p::G4ProbeSpace, i::Int) = cellcentroid(p.space, i)
cellat(p::G4ProbeSpace, x) = cellat(p.space, x)
nchunks(p::G4ProbeSpace) = nchunks(p.space)
cellindices(p::G4ProbeSpace, c::Int) = cellindices(p.space, c)
celltree(p::G4ProbeSpace) = celltree(p.space)

# A space that genuinely IS a sub-space of another: its chunks are `chunks` of
# the parent, in that order, and its cells are those chunks' cells renumbered
# `1:n`. Its geometry is the parent's, so its chunk caps ARE the parent's caps
# for those chunks — which is precisely the precondition
# `subspace_dependencies` checks before re-stamping a row view onto it. This is
# the shape `scripts/copdem_production.jl`'s per-column destination has.
struct E1Subspace{S<:RegridSpace} <: RegridSpace
    parent::S
    chunks::Vector{Int}
    cells::Vector{Int}                  # local cell -> parent cell position
    spans::Vector{UnitRange{Int}}       # local cells of each local chunk
end

function E1Subspace(parent::RegridSpace, chunks)
    cells = Int[]
    spans = UnitRange{Int}[]
    for c in chunks
        lo = length(cells) + 1
        append!(cells, Int.(cellindices(parent, Int(c))))
        push!(spans, lo:length(cells))
    end
    return E1Subspace(parent, collect(Int, chunks), cells, spans)
end

ncells(s::E1Subspace) = length(s.cells)
nchunks(s::E1Subspace) = length(s.chunks)
getcell(s::E1Subspace, i::Int) = getcell(s.parent, s.cells[i])
manifold(s::E1Subspace) = manifold(s.parent)
cellindices(s::E1Subspace, c::Int) = s.spans[c]
# `chunkextents` is public but not exported, so this must name it to extend it
# rather than define a second function that the generic fallback would shadow.
GR.chunkextents(s::E1Subspace) = GR.chunkextents(s.parent)[s.chunks]

@testset "chunk dependency graph" begin

    @testset "absolute support-radius dilation" begin
        # Regridding support is an angle, not an Extents growth factor. The
        # implementation constructs that absolute dilation once and delegates
        # the decision to GeometryOps' guarded cap-cap predicate.
        base = SphericalCap(toy_point(0.0, 0.0), 0.2)
        @test GR._dilatedcap(base, 0.1).radius == base.radius + 0.1
        @test GR._dilatedcap(base, 0.0) === base

        cases = (
            # Tiny caps whose centres are only picoradians apart.
            (SphericalCap(USPoint(1.0, 0.0, 0.0), 1e-12),
             SphericalCap(USPoint(cos(2.5e-12), sin(2.5e-12), 0.0), 1e-12),
             1e-12, true),
            # Just inside and just outside external contact.
            (SphericalCap(toy_point(0.0, 0.0), 0.1),
             SphericalCap(USPoint(cos(0.25 - 1e-10), sin(0.25 - 1e-10), 0.0), 0.1),
             0.05, true),
            (SphericalCap(toy_point(0.0, 0.0), 0.1),
             SphericalCap(USPoint(cos(0.25 + 1e-10), sin(0.25 + 1e-10), 0.0), 0.1),
             0.05, false),
            # Polar and antimeridian pairs exercise coordinate-independent
            # dilation; only the spherical centre distance matters.
            (SphericalCap(toy_point(0.0, 89.9), 0.001),
             SphericalCap(toy_point(180.0, 89.9), 0.001), 0.0016, true),
            (SphericalCap(toy_point(179.9, 0.0), 0.001),
             SphericalCap(toy_point(-179.9, 0.0), 0.001), 0.0014, false),
        )
        for (a, b, radius, expected) in cases
            pred = GR.DilatedIntersects(radius)
            public_result = GO.Extents.intersects(GR._dilatedcap(a, radius), b)
            legacy_result = US.spherical_distance(a.point, b.point) <=
                            a.radius + b.radius + radius
            @test pred(a, b) == public_result == legacy_result == expected
            @test pred(b, a) == expected
        end

        # Closed caps include exact tangency.
        tangent_a = SphericalCap(toy_point(0.0, 0.0), 0.1)
        tangent_b = SphericalCap(USPoint(cos(0.25), sin(0.25), 0.0), 0.1)
        @test GR.DilatedIntersects(0.05)(tangent_a, tangent_b)
    end

    @testset "matches brute force" begin
        # Offset grids so chunk caps genuinely straddle each other, and include
        # polar chunks, where caps are widest.
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))

        for radius in (0.0, 0.05, 0.4)
            # A toy space's index tests the caps `chunkextents` reports, so on
            # this pair the relation is exactly the brute-force cap join.
            g = GR.chunk_dependency_graph(dst, src; radius)
            @test graph_pairs(g) == capjoin_pairs(dst, src; radius)
        end

        # A radius wide enough to reach every pair: the relation must then be
        # complete, not silently truncated by the index's own pruning.
        full = GR.chunk_dependency_graph(dst, src; radius = 4.0)
        @test Graphs.ne(full) == nchunks(dst) * nchunks(src)
    end

    @testset "the graph dominates what a lazy read demands" begin
        # A refcount built from `consumersof` is only sound if every source
        # chunk a read can ask for is an edge. The graph is not free to use its
        # own idea of a source chunk's cap: a space's index tests the extents
        # its own hierarchy derives, and a cap derived from a node rectangle
        # neither contains nor is contained in the cap derived from the chunk's
        # own boundary. The two relations then CROSS, and the pairs on the
        # index's side of the difference are demands no refcount predicted.
        #
        # A `RasterGrid` source is the shipped case: its index is a quadtree
        # cursor over cells whose leaves straddle chunk boundaries, and it
        # answers with whole chunks that the chunk's own cap would reject.
        raster = misalignedraster(360, 180, 37, 23)
        dst = ToyLonLatSpace(24, 12; chunks = (2, 2))

        for radius in (0.0, 0.02)
            g = GR.chunk_dependency_graph(dst, raster; radius)
            demanded = demanded_pairs(dst, raster; radius)
            @test !isempty(demanded)
            @test demanded ⊆ graph_pairs(g)
        end

        # Toy spaces have no such divergence — their index tests the very caps
        # `chunkextents` reports — so the same assertion must also hold there,
        # and there it must hold with equality.
        toysrc = ToyLonLatSpace(8, 4; chunks = (2, 1))
        gt = GR.chunk_dependency_graph(dst, toysrc; radius = 0.05)
        @test graph_pairs(gt) == demanded_pairs(dst, toysrc; radius = 0.05)
    end

    @testset "no geometrically contributing pair is dropped" begin
        # The tests above compare one conservative relation against another;
        # both could be conservative about the wrong thing together. This one
        # builds REAL cell geometry and asks the only question that matters: is
        # every chunk pair a weight could be nonzero on an edge?
        #
        # `contributing_pairs` reaches that answer without a cap, a chunk extent
        # or a chunk index anywhere in its construction, so it is independent of
        # everything the builder does. A builder swap that keeps the candidate
        # tests green and breaks this one has lost real weight.
        cases = (
            # A generic (packed R-tree) source index, at zero and at a support
            # radius wide enough to connect chunks that do not touch.
            ("toy source", ToyLonLatSpace(12, 6; chunks = (3, 2)),
                ToyLonLatSpace(8, 4; chunks = (2, 1)), 0.0),
            ("toy source, support", ToyLonLatSpace(12, 6; chunks = (3, 2)),
                ToyLonLatSpace(8, 4; chunks = (2, 1)), 0.2),
            # The shipped raster path: a quadtree cursor whose leaves straddle
            # chunk boundaries, so its answers are not the chunk caps'.
            ("raster source", ToyLonLatSpace(24, 12; chunks = (2, 2)),
                misalignedraster(36, 18, 7, 5), 0.0),
            ("raster source, support", ToyLonLatSpace(24, 12; chunks = (2, 2)),
                misalignedraster(36, 18, 7, 5), 0.1),
            # Chunks that do not divide the raster evenly: the last chunk of
            # each axis is a different size from the rest.
            ("nonuniform raster chunks", ToyLonLatSpace(18, 9; chunks = (4, 3)),
                misalignedraster(37, 19, 11, 6), 0.0),
            # A polar source band: the widest caps, and where a cell's own cap
            # and its hierarchy's node extents diverge most.
            ("polar source", ToyLonLatSpace(24, 12; chunks = (3, 3)),
                ToyLonLatSpace(8, 2; lat = (60.0, 90.0), chunks = (3, 1)), 0.0),
            # A regional source against the antimeridian, so the destination has
            # isolated chunks and the source's own coverage wraps the cut.
            ("antimeridian source", ToyLonLatSpace(24, 12; chunks = (2, 2)),
                ToyLonLatSpace(6, 4; lon = (150.0, 180.0), lat = (-40.0, 40.0),
                    chunks = (2, 2)), 0.0),
        )
        for (name, dst, src, radius) in cases
            @testset "$name" begin
                truth = contributing_pairs(dst, src; radius)
                graph = graph_pairs(GR.chunk_dependency_graph(dst, src; radius))
                @test !isempty(truth)
                @test truth ⊆ graph
                # `truth ⊆ graph` alone would pass on a builder that returned
                # every pair, so pin the relation exactly as well: post-#69 the
                # rows ARE the `candidatechunks!` answers, with no `refine`, so
                # the graph and the demanded relation are equal and not merely
                # nested. That equality is what a builder swap must preserve.
                @test graph == demanded_pairs(dst, src; radius)
                # A chunk cap covers its own cells, so the cap join must hold
                # the same pairs. This is the `chunkextents` half of the same
                # obligation, and it fails if a chunk cap is too tight.
                @test truth ⊆ capjoin_pairs(dst, src; radius)
            end
        end
    end

    @testset "eager weights ⊆ graph rows ⊆ the broad cap relation" begin
        # Task E1's proof obligation, in the three terms it is stated in. The
        # middle one is what the lazy executor now reads: a tile takes its
        # source chunks from `sourcesof` and from nothing else.
        #
        # Left: the relation must hold every pair the EAGER path put a nonzero
        # weight on. If it did not, a lazy read would never load a source the
        # eager one used, and the two would return different numbers — which is
        # the failure the whole card is guarding against. `eager_pairs` reads
        # the weight builder's own answer, not a cap or an index.
        #
        # Right: the relation must stay inside the brute-force cap join, which
        # is what "conservative bound" means as opposed to "any superset at
        # all". That half holds where the source's chunk index is the generic
        # one over the caps `chunkextents` reports. On a native hierarchy the
        # two relations CROSS in both directions — pinned in the testset below —
        # so the raster arm asserts what is true there instead, and says so.
        generic = (
            ("toy source", ToyLonLatSpace(12, 6; chunks = (3, 2)),
                ToyLonLatSpace(8, 4; chunks = (2, 1))),
            ("polar source", ToyLonLatSpace(24, 12; chunks = (3, 3)),
                ToyLonLatSpace(8, 2; lat = (60.0, 90.0), chunks = (3, 1))),
            ("antimeridian source", ToyLonLatSpace(24, 12; chunks = (2, 2)),
                ToyLonLatSpace(6, 4; lon = (150.0, 180.0), lat = (-40.0, 40.0),
                    chunks = (2, 2))),
            ("coarser destination", ToyLonLatSpace(6, 4; chunks = (2, 2)),
                ToyLonLatSpace(18, 9; chunks = (3, 3))),
        )
        for (name, dst, src) in generic
            @testset "$name" begin
                # The rows the executor reads are the plan's own relation, taken
                # exactly as `LazyRegridArray` takes it.
                rows = graph_pairs(GR.dependencies(
                    ChunkedPlan(GR.Conservative(), Weighted(0.5), dst, src)))
                eager = eager_pairs(dst, src)
                @test !isempty(eager)
                @test eager ⊆ rows
                @test rows ⊆ capjoin_pairs(dst, src)
            end
        end

        @testset "raster source" begin
            dst = ToyLonLatSpace(24, 12; chunks = (2, 2))
            src = misalignedraster(36, 18, 7, 5)
            rows = graph_pairs(GR.dependencies(
                ChunkedPlan(GR.Conservative(), Weighted(0.5), dst, src)))
            eager = eager_pairs(dst, src)
            @test !isempty(eager)
            @test eager ⊆ rows
            # The right-hand containment does not hold on a native hierarchy and
            # is not meant to: the quadtree answers whole chunks for a
            # straddling leaf. Both supersets still dominate the eager weights,
            # which is the invariant that matters for correctness.
            @test eager ⊆ capjoin_pairs(dst, src)
        end
    end

    @testset "the cap join is an identity only on the generic index" begin
        # G1's second gate. A space with no chunk index of its own is answered
        # by a packed R-tree over the very caps `chunkextents` reports, so there
        # the relation must be EXACTLY the brute-force cap join — an equality,
        # not a containment.
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        for radius in (0.0, 0.05, 0.4)
            @test graph_pairs(GR.chunk_dependency_graph(dst, src; radius)) ==
                  capjoin_pairs(dst, src; radius)
        end

        # A native hierarchy is a different matter, and the difference is the
        # reason PR #69 exists. The raster quadtree answers whole chunks for a
        # straddling leaf without testing each chunk's own cap, so it holds
        # pairs the cap join rejects; and the cap join holds pairs the quadtree
        # never reaches. The two relations CROSS in both directions — measured
        # here, so the cap join cannot be used as a bound on the graph in
        # either direction, only as a second superset of the geometric truth.
        rdst = ToyLonLatSpace(24, 12; chunks = (2, 2))
        rsrc = misalignedraster(36, 18, 7, 5)
        native = graph_pairs(GR.chunk_dependency_graph(rdst, rsrc))
        capjoin = capjoin_pairs(rdst, rsrc)
        @test !isempty(setdiff(native, capjoin))
        @test !isempty(setdiff(capjoin, native))
        # Both still dominate the geometric truth; that is the invariant, and
        # it is the one a builder change must keep.
        truth = contributing_pairs(rdst, rsrc)
        @test truth ⊆ native
        @test truth ⊆ capjoin
    end

    @testset "regional and polar spaces" begin
        # A regional source touches only part of the destination: the graph must
        # have isolated destination chunks rather than dropping rows.
        dst = ToyLonLatSpace(8, 4; chunks = (2, 1))
        src = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (2, 1))
        g = GR.chunk_dependency_graph(dst, src)
        @test graph_pairs(g) == capjoin_pairs(dst, src)
        @test GR.ndestinationchunks(g) == nchunks(dst)
        @test GR.nsourcechunks(g) == nchunks(src)
        @test any(d -> GR.sourcedegree(g, d) == 0, 1:GR.ndestinationchunks(g))

        # Every source chunk of a global source is needed by someone.
        global_src = ToyLonLatSpace(4, 2; chunks = (2, 1))
        gg = GR.chunk_dependency_graph(dst, global_src)
        @test all(s -> GR.consumerdegree(gg, s) > 0, 1:GR.nsourcechunks(gg))
    end

    @testset "the two directions are transposes" begin
        dst = ToyLonLatSpace(10, 6; chunks = (2, 2))
        src = ToyLonLatSpace(6, 4; chunks = (3, 2))
        g = GR.chunk_dependency_graph(dst, src; radius = 0.02)

        forward = Set((d, Int(s)) for d in 1:GR.ndestinationchunks(g)
                      for s in GR.sourcesof(g, d))
        reverse = Set((Int(d), s) for s in 1:GR.nsourcechunks(g)
                      for d in GR.consumersof(g, s))
        @test forward == reverse
        @test Graphs.ne(g) == length(forward)

        # Rows are ascending and duplicate-free, which `insorted` lookups and
        # merge-based schedulers both rely on.
        @test all(issorted(GR.sourcesof(g, d); lt = <)
                  for d in 1:GR.ndestinationchunks(g))
        @test all(allunique(GR.consumersof(g, s)) for s in 1:GR.nsourcechunks(g))

        # Degrees agree with the rows and with Graphs.jl's own accessor.
        @test all(GR.sourcedegree(g, d) == length(GR.sourcesof(g, d))
                  for d in 1:GR.ndestinationchunks(g))
        @test sum(GR.consumerdegree(g, s) for s in 1:GR.nsourcechunks(g)) == Graphs.ne(g)
        @test all(Graphs.degree(g, GR.dstvertex(g, d)) == GR.sourcedegree(g, d)
                  for d in 1:GR.ndestinationchunks(g))
    end

    @testset "vertex numbering and role checks" begin
        dst = ToyLonLatSpace(8, 4; chunks = (2, 2))
        src = ToyLonLatSpace(4, 2; chunks = (2, 1))
        g = GR.chunk_dependency_graph(dst, src)
        nsrc, ndst = GR.nsourcechunks(g), GR.ndestinationchunks(g)

        @test GR.srcvertices(g) == 1:nsrc
        @test GR.dstvertices(g) == (nsrc+1):(nsrc+ndst)
        @test Graphs.nv(g) == nsrc + ndst
        @test all(GR.srcchunk(g, GR.srcvertex(g, s)) == s for s in 1:nsrc)
        @test all(GR.dstchunk(g, GR.dstvertex(g, d)) == d for d in 1:ndst)
        @test GR.issrcvertex(g, 1) && !GR.isdstvertex(g, 1)
        @test GR.isdstvertex(g, nsrc + 1) && !GR.issrcvertex(g, nsrc + 1)

        # Reading a vertex in the wrong role is an error, not a silent offset.
        @test_throws ArgumentError GR.srcchunk(g, nsrc + 1)
        @test_throws ArgumentError GR.dstchunk(g, 1)
        @test_throws BoundsError GR.sourcesof(g, ndst + 1)
        @test_throws BoundsError GR.consumersof(g, 0)

        # A destination's neighbours are source chunk numbers unchanged; a
        # source's neighbours are shifted destination vertices.
        d = findfirst(i -> GR.sourcedegree(g, i) > 0, 1:ndst)
        @test Graphs.outneighbors(g, GR.dstvertex(g, d)) == GR.sourcesof(g, d)
        s = Int(first(GR.sourcesof(g, d)))
        @test GR.dstvertex(g, d) in Graphs.outneighbors(g, GR.srcvertex(g, s))
        @test Graphs.inneighbors(g, GR.srcvertex(g, s)) ==
              Graphs.outneighbors(g, GR.srcvertex(g, s))
    end

    @testset "argument validation" begin
        dst = ToyLonLatSpace(4, 2; chunks = (2, 1))
        src = ToyLonLatSpace(4, 2; chunks = (2, 1))
        @test_throws ArgumentError GR.chunk_dependency_graph(dst, src; radius = -1.0)
        @test_throws ArgumentError GR.chunk_dependency_graph(dst, src; radius = NaN)
        @test_throws ArgumentError GR.chunk_dependency_graph(dst, src; radius = Inf)
    end

    @testset "refine narrows conservatively" begin
        dst = ToyLonLatSpace(8, 4; chunks = (2, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 2))
        full = GR.chunk_dependency_graph(dst, src)

        # A refinement that rejects everything must empty the graph; one that
        # accepts everything must leave it untouched.
        none = planned_dependencies(dst, src; refine = (d, s) -> false)
        @test Graphs.ne(none) == 0
        @test all(isempty(GR.sourcesof(none, d)) for d in 1:GR.ndestinationchunks(none))
        all_kept = planned_dependencies(dst, src; refine = (d, s) -> true)
        @test graph_pairs(all_kept) == graph_pairs(full)

        # A selective refinement drops exactly the rejected pairs, in both CSR
        # directions.
        odd = planned_dependencies(dst, src; refine = (d, s) -> isodd(s))
        @test graph_pairs(odd) == Set(p for p in graph_pairs(full) if isodd(p[2]))
        @test Graphs.ne(odd) == length(graph_pairs(odd))
        # The rejection must reach the reverse CSR too, not just the forward one.
        @test all(GR.consumerdegree(odd, s) == 0
                  for s in 1:GR.nsourcechunks(odd) if iseven(s))
    end

    @testset "Graphs.jl interface conformance" begin
        dst = ToyLonLatSpace(10, 6; chunks = (2, 2))
        src = ToyLonLatSpace(6, 4; chunks = (2, 2))
        g = GR.chunk_dependency_graph(dst, src)
        nsrc = GR.nsourcechunks(g)

        # Interface basics.
        @test g isa Graphs.AbstractGraph{Int32}
        @test eltype(g) == Int32
        @test !Graphs.is_directed(g)
        @test !Graphs.is_directed(typeof(g))
        @test Graphs.edgetype(g) == Graphs.SimpleEdge{Int32}
        @test Graphs.has_vertex(g, 1) && !Graphs.has_vertex(g, Graphs.nv(g) + 1)
        @test Graphs.nv(zero(typeof(g))) == 0 && Graphs.ne(zero(typeof(g))) == 0

        # `edges` yields each undirected edge exactly once, source side first.
        es = collect(Graphs.edges(g))
        @test length(es) == Graphs.ne(g)
        @test length(unique(es)) == Graphs.ne(g)
        @test all(e -> Graphs.src(e) <= nsrc < Graphs.dst(e), es)
        @test Set((Graphs.dst(e) - nsrc, Int(Graphs.src(e))) for e in es) == graph_pairs(g)

        # `has_edge` agrees with the edge list, is symmetric, and rejects
        # same-side pairs — the bipartite invariant.
        @test all(Graphs.has_edge(g, Graphs.src(e), Graphs.dst(e)) for e in es)
        @test all(Graphs.has_edge(g, Graphs.dst(e), Graphs.src(e)) for e in es)
        @test !Graphs.has_edge(g, 1, 2)
        @test !Graphs.has_edge(g, nsrc + 1, nsrc + 2)
        @test !Graphs.has_edge(g, 1, Graphs.nv(g) + 1)

        # Neighbour lists and the edge list are the same relation.
        @test sum(length(Graphs.outneighbors(g, v)) for v in Graphs.vertices(g)) ==
              2 * Graphs.ne(g)

        # Stock Graphs.jl algorithms run against the type. This is the proof of
        # conformance: nothing here knows about chunks.
        @test Graphs.is_bipartite(g)
        comps = Graphs.connected_components(g)
        @test sum(length, comps) == Graphs.nv(g)
        @test all(!isempty, comps)
        # A global source reaches everything, so the whole graph is one piece.
        @test length(comps) == 1

        parents = Graphs.bfs_parents(g, 1)
        @test length(parents) == Graphs.nv(g)
        @test parents[1] == 1
        @test Graphs.dijkstra_shortest_paths(g, 1).dists[1] == 0
        @test maximum(Graphs.degree(g)) ==
              max(maximum(GR.sourcedegree(g, d) for d in 1:GR.ndestinationchunks(g)),
                  maximum(GR.consumerdegree(g, s) for s in 1:nsrc))

        # Disconnected inputs really do split into components, which is what
        # makes `connected_components` a scheduling primitive here.
        far = ToyLonLatSpace(4, 2; lon = (-180.0, -140.0), lat = (-80.0, -60.0),
            chunks = (2, 1))
        wide = ToyLonLatSpace(8, 4; chunks = (2, 2))
        split = GR.chunk_dependency_graph(wide, far)
        @test length(Graphs.connected_components(split)) > 1
    end

    @testset "graph identity carries the radius" begin
        dst = ToyLonLatSpace(8, 4; chunks = (2, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 2))
        tight = GR.chunk_dependency_graph(dst, src; radius = 0.0)
        loose = GR.chunk_dependency_graph(dst, src; radius = 0.5)
        @test GR.dependency_radius(tight) == 0.0
        @test GR.dependency_radius(loose) == 0.5
        # A wider radius can only add edges: the relation is monotone, which is
        # what makes "valid for any method with a smaller support" true.
        @test graph_pairs(tight) ⊆ graph_pairs(loose)
        @test occursin("ChunkDependencyGraph", sprint(show, tight))
    end

    @testset "space stamps distinguish the inputs" begin
        # A stamp is a fingerprint, not a proof, so what is under test is that
        # every way of being a *different* pair of spaces moves it. The two
        # things it cannot catch — a digest collision, and equal caps over
        # different cell geometry — are documented on `spacestamp` and are not
        # testable here by construction.
        base = ToyLonLatSpace(8, 4; chunks = (2, 2))
        @test GR.spacestamp(base) == GR.spacestamp(ToyLonLatSpace(8, 4; chunks = (2, 2)))

        # A different resolution, a different chunking, and a different extent
        # each move it, and so does a different space type over the same globe.
        @test GR.spacestamp(base) != GR.spacestamp(ToyLonLatSpace(16, 8; chunks = (2, 2)))
        @test GR.spacestamp(base) != GR.spacestamp(ToyLonLatSpace(8, 4; chunks = (4, 2)))
        @test GR.spacestamp(base) != GR.spacestamp(
            ToyLonLatSpace(8, 4; lon = (-90.0, 90.0), chunks = (2, 2)))
        raster = misalignedraster(8, 4, 2, 2)
        @test GR.spacestamp(raster).tag != GR.spacestamp(base).tag

        # The components are readable, not just hashed, so a mismatch can be
        # reported in terms a caller recognizes.
        s = GR.spacestamp(base)
        @test s.ncells == ncells(base)
        @test s.nchunks == nchunks(base)
        @test occursin("ToyLonLatSpace", string(s.tag))

        # Stamps, identities and whole graphs survive serialization: the point
        # of a fingerprint over an object reference is that it can be written
        # down beside the relation and read back.
        roundtrip(x) = (io = IOBuffer(); Serialization.serialize(io, x);
                        seekstart(io); Serialization.deserialize(io))
        @test roundtrip(s) == s
        id = GR.dependency_identity(base, base; radius = 0.25, narrow = :box)
        @test roundtrip(id) == id
        @test id.dst == id.src == s
        @test id.radius == 0.25 && id.narrow == :box
    end

    @testset "a graph carries the identity of what built it" begin
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        g = GR.chunk_dependency_graph(dst, src; radius = 0.1)
        id = GR.dependency_identity(g)

        @test id == GR.dependency_identity(dst, src; radius = 0.1)
        @test id.dst == GR.spacestamp(dst)
        @test id.src == GR.spacestamp(src)
        @test id.radius == GR.dependency_radius(g) == 0.1
        @test GR.narrowphase(g) == :none
        # The stamps carry the chunk counts the graph was built over.
        @test id.dst.nchunks == GR.ndestinationchunks(g) == nchunks(dst)
        @test id.src.nchunks == GR.nsourcechunks(g) == nchunks(src)
        # Roles are not symmetric: swapping the spaces makes a different graph
        # with a different identity, even at the same counts.
        @test GR.dependency_identity(GR.chunk_dependency_graph(src, dst; radius = 0.1)) != id
    end

    @testset "the narrow-phase tag is a claim about the relation" begin
        dst = ToyLonLatSpace(8, 4; chunks = (2, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 2))
        odd = (d, s) -> isodd(s)

        # No refine and no tag is the only combination that means ":none".
        @test GR.narrowphase(GR.chunk_dependency_graph(dst, src)) == :none
        @test GR.narrowphase(planned_dependencies(dst, src)) == :none
        # A closure has no identity a stamp can carry, so an unnamed refine is
        # recorded as exactly that, and can never be validated for reuse.
        @test GR.narrowphase(planned_dependencies(dst, src; refine = odd)) ==
              GR.UNNAMED_NARROW
        @test GR.narrowphase(planned_dependencies(dst, src; refine = odd,
            narrow = :oddsources)) == :oddsources

        # Both halves of the claim must agree. Naming a phase that was not
        # applied, or applying one and claiming none, is an error at
        # construction rather than a lie in the record.
        @test_throws ArgumentError planned_dependencies(dst, src; narrow = :oddsources)
        @test_throws ArgumentError planned_dependencies(dst, src; refine = odd,
            narrow = :none)
    end

    @testset "an invalid reuse fails at construction" begin
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        g = GR.chunk_dependency_graph(dst, src; radius = 0.1)

        # The valid reuse: the same pair, at this radius or any narrower one.
        @test GR.validate_dependencies(g, dst, src; radius = 0.1) === g
        @test GR.validate_dependencies(g, dst, src; radius = 0.0) === g

        # A wider radius is NOT valid: the relation is monotone in the radius,
        # so a graph built tight is not a superset of a looser demand.
        @test_throws ArgumentError GR.validate_dependencies(g, dst, src; radius = 0.2)
        @test_throws ArgumentError GR.validate_dependencies(g, dst, src; radius = -1.0)

        # Either space being a different space fails, and so does swapping the
        # two: the stamps are compared in their own roles.
        other_dst = ToyLonLatSpace(12, 6; chunks = (2, 2))
        other_src = ToyLonLatSpace(16, 8; chunks = (2, 1))
        @test_throws ArgumentError GR.validate_dependencies(g, other_dst, src; radius = 0.1)
        @test_throws ArgumentError GR.validate_dependencies(g, dst, other_src; radius = 0.1)
        @test_throws ArgumentError GR.validate_dependencies(g, src, dst; radius = 0.1)
        @test_throws ArgumentError GR.validate_dependencies(g, dst,
            misalignedraster(8, 4, 2, 1); radius = 0.1)
        # The message names the component that moved, not just "mismatch".
        @test occursin("chunks against",
            sprint(showerror, try
                GR.validate_dependencies(g, other_dst, src; radius = 0.1)
            catch e
                e
            end))

        # A narrow phase in either direction. A caller expecting the full
        # candidate relation must not silently receive one somebody thinned,
        # and a caller expecting a thinning must not receive the full one.
        thinned = planned_dependencies(dst, src; radius = 0.1,
            refine = (d, s) -> isodd(s), narrow = :oddsources)
        @test GR.validate_dependencies(thinned, dst, src; radius = 0.1,
            narrow = :oddsources) === thinned
        @test_throws ArgumentError GR.validate_dependencies(thinned, dst, src; radius = 0.1)
        @test_throws ArgumentError GR.validate_dependencies(g, dst, src; radius = 0.1,
            narrow = :oddsources)

        # An unnamed refine is unusable in both directions.
        anon = planned_dependencies(dst, src; radius = 0.1, refine = (d, s) -> isodd(s))
        @test_throws ArgumentError GR.validate_dependencies(anon, dst, src; radius = 0.1)
        @test_throws ArgumentError GR.validate_dependencies(anon, dst, src; radius = 0.1,
            narrow = GR.UNNAMED_NARROW)
        @test_throws ArgumentError GR.validate_dependencies(g, dst, src; radius = 0.1,
            narrow = GR.UNNAMED_NARROW)
    end

    @testset "row views share the parent relation" begin
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        g = GR.chunk_dependency_graph(dst, src; radius = 0.05)
        sel = [2, 3, 7, 11]
        view = GR.restrict(g, sel)

        # It is the same type, over the same source side, with the selected
        # rows renumbered 1:k.
        @test view isa GR.ChunkDependencyGraph
        @test GR.isrestricted(view) && !GR.isrestricted(g)
        @test GR.ndestinationchunks(view) == length(sel)
        @test GR.nsourcechunks(view) == GR.nsourcechunks(g)

        # THE POINT: the destination-major direction is the parent's arrays, by
        # reference. Nothing about the relation is recomputed, and no space is
        # queried. If this ever becomes a copy, `restrict` has stopped being a
        # view and the per-column saving is gone.
        @test view.dstoff === g.dstoff
        @test view.srcof === g.srcof

        # Rows are the parent's rows, and they still know their global chunk.
        @test all(GR.sourcesof(view, i) == GR.sourcesof(g, sel[i]) for i in eachindex(sel))
        @test all(GR.sourcedegree(view, i) == GR.sourcedegree(g, sel[i])
                  for i in eachindex(sel))
        @test GR.globaldestinations(view) == sel
        @test all(GR.globaldestination(view, i) == sel[i] for i in eachindex(sel))
        @test all(GR.localdestination(view, sel[i]) == i for i in eachindex(sel))
        @test GR.localdestination(view, 1) === nothing
        @test GR.globaldestinations(g) == 1:GR.ndestinationchunks(g)
        @test GR.localdestination(g, 4) == 4

        # A view's refcounts are the view's own. `consumersof` must count only
        # the rows the view holds, or a refcount taken from it retires nothing.
        @test Graphs.ne(view) == sum(GR.sourcedegree(g, d) for d in sel)
        @test all(GR.consumerdegree(view, s) <= GR.consumerdegree(g, s)
                  for s in 1:GR.nsourcechunks(g))
        @test sum(GR.consumerdegree(view, s) for s in 1:GR.nsourcechunks(view)) ==
              Graphs.ne(view)
        @test all(GR.consumersof(view, s) ==
                  [i for i in eachindex(sel) if insorted(Int32(s), GR.sourcesof(g, sel[i]))]
                  for s in 1:GR.nsourcechunks(view))

        # Both CSR directions are still transposes of one another, and sorted.
        forward = Set((d, Int(s)) for d in 1:GR.ndestinationchunks(view)
                      for s in GR.sourcesof(view, d))
        reverse = Set((Int(d), s) for s in 1:GR.nsourcechunks(view)
                      for d in GR.consumersof(view, s))
        @test forward == reverse
        @test forward == Set((i, s) for (i, d) in pairs(sel)
                             for s in Int.(GR.sourcesof(g, d)))
        @test all(issorted(GR.consumersof(view, s); lt = <)
                  for s in 1:GR.nsourcechunks(view))

        # The identity still stamps the WHOLE destination space: a view knows
        # what it is a view of, which is what makes global refinement possible.
        @test GR.dependency_identity(view) == GR.dependency_identity(g)
        @test GR.dependency_identity(view).dst.nchunks == nchunks(dst)
        @test occursin("row view", sprint(show, view))

        # Graphs.jl still works over the same CSR — the view is not a second
        # graph type with a second interface.
        @test Graphs.nv(view) == GR.nsourcechunks(view) + length(sel)
        @test Graphs.ne(view) == length(collect(Graphs.edges(view)))
        @test Graphs.is_bipartite(view)
        @test sum(length(Graphs.outneighbors(view, v)) for v in Graphs.vertices(view)) ==
              2 * Graphs.ne(view)
        @test all(Graphs.has_edge(view, Graphs.src(e), Graphs.dst(e))
                  for e in Graphs.edges(view))
        @test Set((Graphs.dst(e) - GR.nsourcechunks(view), Int(Graphs.src(e)))
                  for e in Graphs.edges(view)) == graph_pairs(view)
        @test sum(length, Graphs.connected_components(view)) == Graphs.nv(view)

        # Views compose, and composing them composes the global numbering.
        inner = GR.restrict(view, [2, 4])
        @test GR.globaldestinations(inner) == [sel[2], sel[4]]
        @test GR.sourcesof(inner, 1) == GR.sourcesof(g, sel[2])
        @test inner.srcof === g.srcof
        @test GR.restrict(g, 1:GR.ndestinationchunks(g)) |> graph_pairs == graph_pairs(g)

        # Degenerate selections are legal and do not need special-casing.
        empty_view = GR.restrict(g, Int[])
        @test GR.ndestinationchunks(empty_view) == 0 && Graphs.ne(empty_view) == 0
        @test all(GR.consumerdegree(empty_view, s) == 0
                  for s in 1:GR.nsourcechunks(empty_view))

        # An ill-formed selection fails before anything is built.
        @test_throws ArgumentError GR.restrict(g, [3, 3])
        @test_throws ArgumentError GR.restrict(g, [5, 2])
        @test_throws ArgumentError GR.restrict(g, [0])
        @test_throws ArgumentError GR.restrict(g, [GR.ndestinationchunks(g) + 1])

        # And a view is validated against the rows its caller means to compute,
        # never against the whole space.
        @test GR.validate_dependencies(view, dst, src; radius = 0.05,
            destinations = sel) === view
        @test_throws ArgumentError GR.validate_dependencies(view, dst, src; radius = 0.05)
        @test_throws ArgumentError GR.validate_dependencies(view, dst, src;
            radius = 0.05, destinations = [2, 3])
        @test_throws ArgumentError GR.validate_dependencies(g, dst, src; radius = 0.05,
            destinations = sel)
    end

    @testset "a row view is the graph of those destinations" begin
        # The oracle question, asked of a view: restricting must not change
        # which pairs are held, only which rows are visible. This is the
        # property any per-column plan taking a row view depends on.
        dst = ToyLonLatSpace(24, 12; chunks = (2, 2))
        src = misalignedraster(36, 18, 7, 5)
        for radius in (0.0, 0.1)
            g = GR.chunk_dependency_graph(dst, src; radius)
            ndst = GR.ndestinationchunks(g)
            sel = collect(1:3:ndst)
            view = GR.restrict(g, sel)
            parent = graph_pairs(g)
            # Read the view's own rows back into GLOBAL destination numbers and
            # compare against the parent's rows for exactly those chunks.
            lifted = Set((GR.globaldestination(view, d), Int(s))
                         for d in 1:GR.ndestinationchunks(view)
                         for s in GR.sourcesof(view, d))
            @test lifted == Set(p for p in parent if p[1] in Set(sel))
            @test !isempty(lifted)
            # And every geometrically contributing pair over those chunks is
            # still there, through the independent cell-geometry oracle.
            truth = contributing_pairs(dst, src; radius)
            @test Set(p for p in truth if p[1] in Set(sel)) ⊆ lifted
        end
    end

    @testset "a row view re-stamped onto a sub-space is that space's relation" begin
        # Task E1. G4 proved that a plan over a *sub-space* of the destination
        # cannot adopt a row view, because the view stamps the whole space. This
        # is the mechanism that closes it, and what makes it sound: the
        # destination half of a relation is a function of the destination caps
        # alone, so a space reproducing those caps reproduces those rows.
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        method = G4RadiusMethod(0.1)
        g = planned_dependencies(dst, src; radius = 0.1)
        sel = [2, 5, 6, 9]
        sub = E1Subspace(dst, sel)
        @test nchunks(sub) == length(sel)
        @test GR.chunkextents(sub) == GR.chunkextents(dst)[sel]

        view = GR.subspace_dependencies(g, sub, sel)
        # Row for row, it is the parent's relation over those chunks.
        for (k, d) in enumerate(sel)
            @test collect(GR.sourcesof(view, k)) == collect(GR.sourcesof(g, d))
        end
        @test GR.globaldestinations(view) == sel     # provenance is retained
        @test GR.nsourcechunks(view) == GR.nsourcechunks(g)
        # ...and it is the relation the sub-space would have built for itself.
        own = GR.chunk_dependency_graph(sub, src; radius = 0.1)
        @test graph_pairs(view) == graph_pairs(own)
        @test GR.dependency_identity(view) == GR.dependency_identity(own)

        # Still a view: both shared arrays are the parent's, not copies.
        @test view.dstoff === g.dstoff
        @test view.srcof === g.srcof
        @test GR.destinationextents(view) === GR.destinationextents(g)
        @test GR.sourceextents(view) === GR.sourceextents(g)
        # Its refcounts are its own, over its own rows.
        @test sum(GR.consumerdegree(view, s) for s in 1:GR.nsourcechunks(view)) ==
              sum(GR.sourcedegree(view, d) for d in 1:GR.ndestinationchunks(view))

        # A plan over the SUB-space adopts it, by reference, with no
        # `destinations` argument: for that space there are no other rows.
        adopted = ChunkedPlan(method, Weighted(0.5), sub, src; dependencies = view)
        @test GR.dependencies(adopted) === view
        # The un-re-stamped view is still refused, which is the G4 behaviour
        # this mechanism exists beside rather than instead of.
        @test_throws ArgumentError ChunkedPlan(method, Weighted(0.5), sub, src;
            dependencies = GR.restrict(g, sel))
        @test_throws ArgumentError ChunkedPlan(method, Weighted(0.5), sub, src;
            dependencies = g)
        # And the re-stamped view is not a relation for the WHOLE destination.
        @test_throws ArgumentError ChunkedPlan(method, Weighted(0.5), dst, src;
            dependencies = view)

        # Refusals, at the re-stamp rather than later.
        @test_throws ArgumentError GR.subspace_dependencies(g, E1Subspace(dst, [1, 2]), sel)
        @test_throws ArgumentError GR.subspace_dependencies(g, sub, [2, 5, 6, 10])
        @test_throws ArgumentError GR.subspace_dependencies(g, sub, [9, 6, 5, 2])
        # A space whose caps are not the parent's caps for those chunks is the
        # one thing that makes the re-stamp unsound, so it is the one thing
        # checked exactly rather than approximately.
        other = E1Subspace(ToyLonLatSpace(12, 6; lat = (-80.0, 80.0), chunks = (3, 2)), sel)
        @test GR.chunkextents(other) != GR.chunkextents(dst)[sel]
        @test_throws ArgumentError GR.subspace_dependencies(g, other, sel)
        # A relation with no extents has nothing to compare against.
        bare = GR.ChunkDependencyGraph(GR.DependencyIdentity(), [1, 1], Int32[],
            [1, 1], Int32[])
        @test_throws ArgumentError GR.subspace_dependencies(bare, sub, [1])
    end

    # ----------------------------------------------------------------------
    # Task G4: the Phase 2 gate. One logical plan exposes exactly one
    # validated relation, and a narrow phase cannot be supplied after the plan
    # exists. Both halves are asserted here rather than argued in prose.
    # ----------------------------------------------------------------------

    @testset "a plan owns exactly one relation" begin
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        method = G4RadiusMethod(0.1)
        policy = Weighted(0.5)

        # The default builds one, at the plan's own radius, in every spelling of
        # the constructor — keyword and both positional forms. Task E1 made this
        # the default because a lazy read IS a read of these rows.
        bare = ChunkedPlan(method, policy, dst, src)
        @test GR.dependencies(bare) isa GR.ChunkDependencyGraph
        @test GR.dependency_radius(GR.dependencies(bare)) == 0.1
        for p in (ChunkedPlan(method, policy, dst, src, PerChunk(), 2^30, nothing),
            ChunkedPlan(method, policy, dst, src, PerChunk(), 2^30, nothing, nothing))
            @test graph_pairs(GR.dependencies(p)) == graph_pairs(GR.dependencies(bare))
        end

        # `dependencies = false` is the opt-out, and construction then does no
        # relation work at all: the probe is never asked a question. A plan that
        # takes it cannot back a `LazyRegridArray`, which is the point.
        quiet = G4ProbeSpace(src)
        @test GR.dependencies(
            ChunkedPlan(method, policy, dst, quiet; dependencies = false)) === nothing
        @test quiet.queries == 0

        # `dependencies = true` builds it once, at the plan's OWN radius — the
        # method's `support_radius`, never a keyword of its own.
        owned = ChunkedPlan(method, policy, dst, src; dependencies = true)
        g = GR.dependencies(owned)
        @test g isa GR.ChunkDependencyGraph
        @test GR.dependency_radius(g) == 0.1
        @test graph_pairs(g) ==
              graph_pairs(GR.chunk_dependency_graph(dst, src; radius = 0.1))
        @test GR.dependency_radius(
            GR.dependencies(ChunkedPlan(G4RadiusMethod(0.4), policy, dst, src;
                dependencies = true))) == 0.4

        # The accessor builds NOTHING: the same object every time, and not one
        # further question to either space after construction.
        probe = G4ProbeSpace(src)
        built = ChunkedPlan(method, policy, dst, probe; dependencies = true)
        asked = probe.queries
        @test asked > 0
        for _ in 1:5
            @test GR.dependencies(built) === GR.dependencies(built)
        end
        @test GR.dependencies(built) === built.dependencies
        @test probe.queries == asked

        # Exactly one: no plan method on the builder, so there is no second
        # relation a plan can be made to hand out.
        @test isempty(methods(GR.chunk_dependency_graph, Tuple{ChunkedPlan}))
        @test_throws MethodError GR.chunk_dependency_graph(owned)
        @test_throws MethodError GR.chunk_dependency_graph(owned;
            refine = (d, s) -> true)

        # And an eager plan has no chunk relation to expose.
        @test GR.dependencies(plan_regrid(zeros(8, 4); to = dst, from = src,
            lazy = false)) === nothing
    end

    @testset "a narrow phase cannot be supplied after the plan exists" begin
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        odd = (d, s) -> isodd(s)
        plan = ChunkedPlan(G4RadiusMethod(0.1), Weighted(0.5), dst, src;
            dependencies = true, refine = odd, narrow = :oddsources)

        # `refine`/`narrow` are keywords of plan construction — `plan_regrid`
        # and the `ChunkedPlan` constructor it forwards to — and of nothing
        # else. Not of the graph builder, not of the row view, not of the
        # one-shot API.
        haskw(f, kw) = any(m -> kw in Base.kwarg_decl(m), methods(f))
        for kw in (:refine, :narrow)
            @test !haskw(GR.chunk_dependency_graph, kw)
            @test !haskw(GR.restrict, kw)
            @test !haskw(GR.validate_dependencies, kw) || kw === :narrow
            @test !haskw(regrid, kw)
            @test !haskw(regrid!, kw)
            @test haskw(plan_regrid, kw)
        end
        @test_throws MethodError regrid(zeros(8, 4); to = dst, from = src,
            lazy = true, refine = odd)

        # The plan's relation is a field of an immutable struct, so it cannot
        # be swapped for another one after the fact either.
        @test !ismutabletype(typeof(plan))
        @test_throws ErrorException setfield!(plan, :dependencies, nothing)

        # The relation records which phase produced it, and is the thinned one.
        @test GR.narrowphase(GR.dependencies(plan)) == :oddsources
        @test graph_pairs(GR.dependencies(plan)) ==
              Set(p for p in graph_pairs(GR.chunk_dependency_graph(dst, src;
                  radius = 0.1)) if isodd(p[2]))
    end

    @testset "a plan validates a relation it did not build" begin
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        method = G4RadiusMethod(0.1)
        policy = Weighted(0.5)
        g = GR.chunk_dependency_graph(dst, src; radius = 0.1)
        adopt(; m = method, d = dst, s = src, kw...) =
            GR.dependencies(ChunkedPlan(m, policy, d, s; kw...))

        # Adopted by reference — this is the sharing that keeps a second plan
        # over the same pair from rebuilding the relation.
        @test adopt(; dependencies = g) === g
        # A narrower method may reuse a wider relation; a wider one may not.
        @test adopt(; dependencies = g, m = G4RadiusMethod(0.0)) === g
        @test_throws ArgumentError adopt(; dependencies = g, m = G4RadiusMethod(0.2))
        # Either space moving, or the pair being swapped, fails at construction.
        @test_throws ArgumentError adopt(; dependencies = g,
            d = ToyLonLatSpace(12, 6; chunks = (2, 2)))
        @test_throws ArgumentError adopt(; dependencies = g,
            s = ToyLonLatSpace(16, 8; chunks = (2, 1)))
        @test_throws ArgumentError adopt(; dependencies = g, d = src, s = dst)

        # A narrow phase applies while a relation is built, so it cannot be
        # applied to one that already exists. Name it instead.
        @test_throws ArgumentError adopt(; dependencies = g, refine = (d, s) -> isodd(s))
        thinned = planned_dependencies(dst, src; radius = 0.1,
            refine = (d, s) -> isodd(s), narrow = :oddsources)
        @test_throws ArgumentError adopt(; dependencies = thinned)
        @test adopt(; dependencies = thinned, narrow = :oddsources) === thinned
        @test_throws ArgumentError adopt(; dependencies = g, narrow = :oddsources)

        # A plan computes its WHOLE destination, so a row view is not a
        # relation any plan may adopt — and neither is a whole-space relation
        # over a DIFFERENT destination space that happens to cover those rows.
        # That pair of refusals is exactly why production's per-column plans
        # own no relation: their destination is a rooted one-chunk grid, not
        # the global space the graph and its views are stamped against.
        @test_throws ArgumentError adopt(; dependencies = GR.restrict(g, [1, 2]))
        onechunk = ToyLonLatSpace(12, 2; chunks = (12, 2))
        @test GR.nchunks(onechunk) == 1
        @test_throws ArgumentError adopt(; dependencies = GR.restrict(g, [1]),
            d = onechunk)
        @test_throws ArgumentError adopt(; dependencies = g, d = onechunk)

        # `false` is "hold none, and say so"; a stray narrow phase beside it is
        # an error rather than something silently dropped.
        @test adopt(; dependencies = false) === nothing
        @test_throws ArgumentError adopt(; dependencies = false, refine = (d, s) -> true)
        @test_throws ArgumentError adopt(; dependencies = false, narrow = :oddsources)
        @test_throws ArgumentError adopt(; dependencies = 7)
    end

    @testset "one query implementation defines every edge" begin
        # The Phase 4 gate, structurally and behaviourally.
        #
        # Structurally: the duplicate spellings are gone from the module — not
        # merely unexported — and `chunkextents` is a required hook rather than
        # a fallback that collects caps back out of a compatibility tree.
        for name in (:connectedchunks, :connectedchunks!, :connectedchunkpairs,
                     :chunktree, :RasterFlatTree, :_collectextents!)
            @test !isdefined(GR, name)
        end
        @test :chunktree ∉ names(GlobalRegridding)
        @test !Base.ispublic(GR, :chunktree)
        # No method for the abstract space type: nothing can fall back to a tree.
        @test !hasmethod(GR.chunkextents, Tuple{RegridSpace})
        @test hasmethod(GR.chunkextents, Tuple{ToyLonLatSpace})
        @test hasmethod(GR.chunkextents, Tuple{RasterGrid})
        # `chunkextent` survives with a cheap specialization, which is the whole
        # reason to keep the singular form at all.
        @test hasmethod(GR.chunkextent, Tuple{RegridSpace,Int})
        @test which(GR.chunkextent, Tuple{RasterGrid,Int}).sig !==
              which(GR.chunkextent, Tuple{RegridSpace,Int}).sig

        # Behaviourally: every edge comes through `candidatechunks!` on the
        # source space's own chunk index, one query per destination cap, and
        # from nothing else. The counter says how many questions were asked; the
        # oracle says the answers ARE the relation.
        src = ToyLonLatSpace(8, 4; chunks = (4, 2))
        dst = ToyLonLatSpace(8, 4; chunks = (8, 2))
        probe = E2QuerySpace(src)
        g = planned_dependencies(dst, probe)
        @test probe.index.queries == nchunks(dst)
        @test graph_pairs(g) == demanded_pairs(dst, src)
        @test GR.nsourcechunks(g) == nchunks(src)

        # A lazy read adds no query of its own: it takes the rows.
        A = LazyRegridArray(zeros(8, 4),
            ChunkedPlan(ToyDiagonalMethod(), Weighted(0.5), dst, probe))
        asked = probe.index.queries
        @test collect(A) == vec(zeros(8, 4))
        @test probe.index.queries == asked

        # And a cheap `chunkextent` agrees with the vector it is cheaper than.
        raster = misalignedraster(16, 8, 5, 3)
        @test [GR.chunkextent(raster, c) for c in 1:nchunks(raster)] ==
              GR.chunkextents(raster)
    end

    @testset "`refine` reaches the plan through `plan_regrid` only" begin
        src = ToyLonLatSpace(8, 4; chunks = (4, 2))
        dst = ToyLonLatSpace(8, 4; chunks = (8, 2))
        data = zeros(8, 4)
        odd = (d, s) -> isodd(s)

        lazyplan = plan_regrid(data; to = dst, from = src, lazy = true,
            refine = odd, narrow = :oddsources)
        @test GR.narrowphase(GR.dependencies(lazyplan)) == :oddsources
        @test GR.narrowphase(GR.dependencies(plan_regrid(data; to = dst, from = src,
            lazy = true))) == :none
        @test GR.dependencies(plan_regrid(data; to = dst, from = src,
            lazy = true, dependencies = true)) isa GR.ChunkDependencyGraph
        @test GR.dependencies(plan_regrid(data; to = dst, from = src,
            lazy = true, dependencies = false)) === nothing

        # An eager plan names each keyword it refuses rather than ignoring it.
        for (name, call) in (
            ("dependencies", () -> plan_regrid(data; to = dst, from = src,
                lazy = false, dependencies = true)),
            ("refine", () -> plan_regrid(data; to = dst, from = src,
                lazy = false, refine = odd)),
            ("narrow", () -> plan_regrid(data; to = dst, from = src,
                lazy = false, narrow = :oddsources)))
            err = try
                call()
                nothing
            catch e
                e
            end
            @test err isa ArgumentError && occursin(name, err.msg)
        end
    end
end
