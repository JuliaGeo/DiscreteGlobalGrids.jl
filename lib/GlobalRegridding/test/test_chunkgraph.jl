import Graphs
import DimensionalData as DD

# Brute-force reference: the definition of the relation, applied to every pair.
function reference_pairs(dst, src; radius = 0.0)
    dcaps, scaps = GR.chunkextents(dst), GR.chunkextents(src)
    intersects = GR.DilatedIntersects(Float64(radius))
    return Set((d, s) for d in eachindex(dcaps), s in eachindex(scaps)
               if intersects(dcaps[d], scaps[s]))
end

graph_pairs(g) = Set((d, Int(s)) for d in 1:GR.ndestinationchunks(g)
                     for s in GR.sourcesof(g, d))

# The relation a lazy read actually issues: one `candidatechunks!` per
# destination chunk, through the source space's own index.
function demanded_pairs(dst, src; radius = 0.0)
    index = GR.chunkindex(src)
    out = Set{Tuple{Int,Int}}()
    buf = Int[]
    for (d, cap) in pairs(GR.chunkextents(dst))
        GR.candidatechunks!(buf, index, cap; radius)
        for s in buf
            push!(out, (d, s))
        end
    end
    return out
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
            @test graph_pairs(g) == reference_pairs(dst, src; radius)
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
        raster = RasterGrid(DD.DimArray(zeros(360, 180),
            (DD.X(range(-179.5, 179.5; length = 360)),
             DD.Y(range(-89.5, 89.5; length = 180)))),
            chunks = ([lo:min(lo + 36, 360) for lo in 1:37:360],
                [lo:min(lo + 22, 180) for lo in 1:23:180]))
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

    @testset "regional and polar spaces" begin
        # A regional source touches only part of the destination: the graph must
        # have isolated destination chunks rather than dropping rows.
        dst = ToyLonLatSpace(8, 4; chunks = (2, 1))
        src = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (2, 1))
        g = GR.chunk_dependency_graph(dst, src)
        @test graph_pairs(g) == reference_pairs(dst, src)
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
end
