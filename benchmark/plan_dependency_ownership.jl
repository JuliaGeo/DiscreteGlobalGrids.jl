# Who owns the chunk dependency relation in the production run, and what a
# per-column plan would pay for one of its own (Task G4).
#
#     DGG_COPDEM_TILELIST=/path/to/tileList90.txt \
#         julia -t 8 --gcthreads=4 --project=benchmark \
#         benchmark/plan_dependency_ownership.jl
#
# The tile list is LOCAL and names tiles; no raster is opened and nothing is
# downloaded. Without `DGG_COPDEM_TILELIST` this script does nothing, exactly as
# `chunk_graph_gates.jl`'s production case does not.
#
# Two questions.
#
#   1. Does production's `dagplan` hand back the relation its own global
#      `GR.ChunkedPlan` owns, and is that relation the one the pair has? G4
#      moved the build inside the plan; this asserts the move changed nothing.
#
#   2. What would a per-column plan pay to own a relation? `regrid_chunk` builds
#      one plan per destination chunk — 66 175 of them — over a rooted one-chunk
#      subtree grid. Three arms over exactly that space:
#
#        A  `dependencies = nothing`  what production does: owns none
#        B  `dependencies = true`     a one-destination graph per column
#        C  `restrict(graph, [d])`    a one-row view of the global relation
#
#      B is what the card forbids. C is what the card suggested instead, and the
#      run also shows that a per-column plan cannot legally hold it: its
#      destination is a DIFFERENT space from the one the graph and its views are
#      stamped against, so `validate_dependencies` refuses both C and the whole
#      global graph. That refusal is the point, not a defect.
#
# Environment
#
#   DGG_COPDEM_TILELIST  local Copernicus tile list (required)
#   G4_NCOL              destination chunks sampled across the covering (25)
#   G4_SAMPLES           timed samples per arm per column (5)
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import Graphs
import Statistics
using Printf: @printf

const TILELIST = get(ENV, "DGG_COPDEM_TILELIST", "")
const NCOL = parse(Int, get(ENV, "G4_NCOL", "25"))
const NSAMP = parse(Int, get(ENV, "G4_SAMPLES", "5"))

if isempty(TILELIST)
    println("SKIPPED: set DGG_COPDEM_TILELIST to a local Copernicus tile list.")
    exit(0)
end
isfile(TILELIST) || error("DGG_COPDEM_TILELIST=$TILELIST is not a file")

# `listedtiles`, `TileIds`, `SubtreeIds`, `covering_chunks`, `dagplan`, `CONFIG`.
include(joinpath(@__DIR__, "..", "scripts", "copdem_production.jl"))

const SYS7 = DGG.IGeo7System()

med(f, n) = Statistics.median([(@timed f())[2] for _ in 1:n])
pairsof(g) = Set((d, Int(s)) for d in 1:GR.ndestinationchunks(g)
                 for s in GR.sourcesof(g, d))

function main_ownership()
    sys = DGG.CopernicusDEMSystem(90)
    tiles = listedtiles(sys, TILELIST, nothing)
    srcspace = DGG.DGGSpace(DGG.PartialGrid(sys, 1, TileIds(sys, tiles)); chunklevel = 0)
    chunks = covering_chunks(SYS7, sys, tiles, 5; nthreads = max(1, Threads.nthreads()))
    @printf("pair: %d destination chunks <- %d tiles, %d threads, julia %s\n",
        length(chunks), GR.nchunks(srcspace), Threads.nthreads(), VERSION)

    # ---- 1. production's own path ----------------------------------------
    dag = dagplan(sys, SYS7, tiles, chunks, srcspace, CONFIG)
    graph = dag.graph
    @printf("\ndagplan: graph === dependencies(globalplan): %s\n",
        graph === GR.dependencies(dag.globalplan))
    @printf("dagplan: %d edges, radius %g, narrow %s; graph %.3f s, order %.3f s\n",
        Graphs.ne(graph), dag.radius, GR.narrowphase(graph), dag.tgraph, dag.torder)

    direct = GR.chunk_dependency_graph(dag.dstspace, srcspace; radius = dag.radius)
    @printf("dagplan relation == chunk_dependency_graph(dst, src; radius): %s (%d vs %d edges)\n",
        pairsof(graph) == pairsof(direct), Graphs.ne(graph), Graphs.ne(direct))
    @printf("dagplan identity == direct identity: %s\n",
        GR.dependency_identity(graph) == GR.dependency_identity(direct))

    # The first call above pays compilation; a second one is the warm number
    # comparable with `chunk_graph_gates.jl`'s median.
    warm = dagplan(sys, SYS7, tiles, chunks, srcspace, CONFIG)
    @printf("dagplan warm: graph %.3f s (first call %.3f s, compilation included)\n",
        warm.tgraph, dag.tgraph)

    dagr = dagplan(sys, SYS7, tiles, chunks, srcspace,
        merge(CONFIG, (refinegraph = true,)))
    @printf("refinegraph = true: %d edges, narrow %s, graph %.3f s, subset of the cap relation: %s\n",
        Graphs.ne(dagr.graph), GR.narrowphase(dagr.graph), dagr.tgraph,
        issubset(pairsof(dagr.graph), pairsof(graph)))

    # ---- 2. the per-column arms ------------------------------------------
    g5 = DGG.levelgrid(SYS7, 5)
    # Exactly the space `regrid_chunk` builds.
    colspace(c) = DGG.DGGSpace(DGG.subtree(SYS7, DGG.cellindex(g5, c), CONFIG.level);
        chunklevel = CONFIG.ancestor)

    picks = chunks[round.(Int, range(1, length(chunks); length = NCOL))]
    ta, tb, tc = Float64[], Float64[], Float64[]
    aa, ab, ac = Int[], Int[], Int[]
    for c in picks
        cs = colspace(c)
        GR.nchunks(cs) == 1 || error("a column is not one chunk")
        d = findfirst(==(c), chunks)
        A() = GR.ChunkedPlan(DGG.Conservative(), DGG.Weighted(0.5), cs, srcspace;
            budget = CONFIG.budget)
        B() = GR.ChunkedPlan(DGG.Conservative(), DGG.Weighted(0.5), cs, srcspace;
            budget = CONFIG.budget, dependencies = true)
        C() = GR.restrict(graph, [d])
        A(); B(); C()
        GR.dependencies(A()) === nothing || error("arm A built a relation")
        push!(ta, med(A, NSAMP)); push!(aa, @allocated A())
        push!(tb, med(B, NSAMP)); push!(ab, @allocated B())
        push!(tc, med(C, NSAMP)); push!(ac, @allocated C())
    end
    m(v) = Statistics.median(v)
    n = length(chunks)
    @printf("\n%-34s %12s %14s %12s %10s\n",
        "per-column arm", "median s", "bytes", "x $n s", "x $n GB")
    for (label, t, a) in (("A  plan, no relation (today)", ta, aa),
                          ("B  plan, dependencies = true", tb, ab),
                          ("C  restrict(graph, [d])", tc, ac))
        @printf("%-34s %12.4g %14d %12.1f %10.1f\n",
            label, m(t), m(a), m(t) * n, m(a) * n / 1e9)
    end
    @printf("\nB/A = %.0fx time, %.0fx bytes;  C/A = %.0fx time, %.0fx bytes\n",
        m(tb) / m(ta), m(ab) / max(m(aa), 1), m(tc) / m(ta), m(ac) / max(m(aa), 1))

    # Neither B's product nor C is a relation a per-column plan may hold.
    cs = colspace(picks[1])
    d = findfirst(==(picks[1]), chunks)
    for (name, g) in (("a one-row view of the global graph", GR.restrict(graph, [d])),
                      ("the whole global graph", graph))
        outcome = try
            GR.ChunkedPlan(DGG.Conservative(), DGG.Weighted(0.5), cs, srcspace;
                dependencies = g)
            "ADOPTED — UNEXPECTED"
        catch e
            e isa ArgumentError ? "refused: " * first(split(e.msg, ". ")) : rethrow()
        end
        @printf("\nper-column plan offered %s:\n  %s\n", name, outcome)
    end
end

main_ownership()
