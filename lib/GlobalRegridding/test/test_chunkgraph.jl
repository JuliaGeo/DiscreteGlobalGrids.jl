import Graphs
import DimensionalData as DD
import Serialization

# The four relations — truth, demand, cap join, and the graph's own rows — are
# defined once, in `graphoracles.jl`, and shared with the root suite and the G1
# harness. Do not re-spell any of them here.
include(joinpath(@__DIR__, "graphoracles.jl"))
using .ChunkGraphOracles: contributing_pairs, graph_pairs, demanded_pairs,
    capjoin_pairs

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
        none = GR.chunk_dependency_graph(dst, src; refine = (d, s) -> false)
        @test Graphs.ne(none) == 0
        @test all(isempty(GR.sourcesof(none, d)) for d in 1:GR.ndestinationchunks(none))
        all_kept = GR.chunk_dependency_graph(dst, src; refine = (d, s) -> true)
        @test graph_pairs(all_kept) == graph_pairs(full)

        # A selective refinement drops exactly the rejected pairs, in both CSR
        # directions.
        odd = GR.chunk_dependency_graph(dst, src; refine = (d, s) -> isodd(s))
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
        # A closure has no identity a stamp can carry, so an unnamed refine is
        # recorded as exactly that, and can never be validated for reuse.
        @test GR.narrowphase(GR.chunk_dependency_graph(dst, src; refine = odd)) ==
              GR.UNNAMED_NARROW
        @test GR.narrowphase(GR.chunk_dependency_graph(dst, src; refine = odd,
            narrow = :oddsources)) == :oddsources

        # Both halves of the claim must agree. Naming a phase that was not
        # applied, or applying one and claiming none, is an error at
        # construction rather than a lie in the record.
        @test_throws ArgumentError GR.chunk_dependency_graph(dst, src; narrow = :oddsources)
        @test_throws ArgumentError GR.chunk_dependency_graph(dst, src; refine = odd,
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
        thinned = GR.chunk_dependency_graph(dst, src; radius = 0.1,
            refine = (d, s) -> isodd(s), narrow = :oddsources)
        @test GR.validate_dependencies(thinned, dst, src; radius = 0.1,
            narrow = :oddsources) === thinned
        @test_throws ArgumentError GR.validate_dependencies(thinned, dst, src; radius = 0.1)
        @test_throws ArgumentError GR.validate_dependencies(g, dst, src; radius = 0.1,
            narrow = :oddsources)

        # An unnamed refine is unusable in both directions.
        anon = GR.chunk_dependency_graph(dst, src; radius = 0.1, refine = (d, s) -> isodd(s))
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
        # property G4's per-column plans depend on.
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
end
