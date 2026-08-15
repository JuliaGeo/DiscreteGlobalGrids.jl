#!/usr/bin/env julia
# Exhaustive-ish deterministic sweep of the subtree halo API, driven by real
# terrain code (see SphericalTerrain.jl).
#
#   julia --project=test examples/geomorphometry/run_sweep.jl
#
# Every root at levels 0 and 1 of every system, every target depth 0..3 (capped
# by grid size), both connectivities, four elevation fields. Exhausting the
# roots is how pentagons, poles, face seams, cube corners and icosahedral
# vertices get covered without having to name them: the classifier only
# *reports* which of those each case touched.

include("Harness.jl")
using .Harness
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, Edge
using Printf

const MAXCELLS = 30_000     # cap on the level grid we will build a full ctx for
const FIELDS = (:noise, :harmonic, :step, :ramp)

function main()
    t0 = time()
    allfails = Failure[]
    ncases = 0
    nmetriccmp = 0
    haloed = 0
    coverage = Dict{Symbol,Int}()
    percase = Tuple{String,Int,Int}[]

    println("="^78)
    println("PART 1  adjacency symmetry (once per grid)")
    println("="^78)
    for sys in DGG.systems(), l in 0:3, conn in (Vertex(), Edge())
        DGG.ncells(sys, l) > MAXCELLS && continue
        f = symmetry_failures(sys, l, conn)
        isempty(f) || (append!(allfails, f);
            @printf("  %-16s L%d %-6s  %d asymmetric pairs\n",
                nameof(typeof(sys)), l, nameof(typeof(conn)), length(f)))
    end
    println("  checked ", length(Harness.SYM_CACHE), " grids")

    println()
    println("="^78)
    println("PART 2  every root, every depth, both connectivities, 4 fields")
    println("="^78)
    for sys in DGG.systems()
        sysfails = 0; syscases = 0
        for rootlevel in 0:1
            DGG.ncells(sys, rootlevel) > 200 && continue
            groot = DGG.levelgrid(sys, rootlevel)
            roots = [DGG.cellindex(groot, p) for p in 1:DGG.ncells(groot)]
            for depth in 0:3
                target = rootlevel + depth
                target > DGG.max_level(sys) && continue
                DGG.ncells(sys, target) > MAXCELLS && continue
                for conn in (Vertex(), Edge()), field in FIELDS
                    for root in roots
                        f, s = check_subtree_case(sys, root, target, conn, field;
                            outside_in = (depth <= 1))
                        ncases += 1; syscases += 1
                        nmetriccmp += s.nmetric
                        haloed += s.nhalo
                        if !isempty(f)
                            append!(allfails, f); sysfails += length(f)
                        end
                        if field === :noise && conn === Vertex()
                            for t in classify_root(sys, root, target, conn)
                                coverage[t] = get(coverage, t, 0) + 1
                            end
                        end
                    end
                end
            end
        end
        push!(percase, (string(nameof(typeof(sys))), syscases, sysfails))
        @printf("  %-16s %6d cases, %d failures  (%.1fs elapsed)\n",
            nameof(typeof(sys)), syscases, sysfails, time() - t0)
    end

    println()
    println("="^78)
    println("PART 3  subset containers: PartialGrid / CellVector / CellLookup")
    println("="^78)
    nsubset = 0
    for sys in DGG.systems(), conn in (Vertex(), Edge())
        level = 3
        DGG.ncells(sys, level) > MAXCELLS && continue
        k = gridctx(sys, level, conn)
        groot = DGG.levelgrid(sys, 1)
        for p in 1:min(DGG.ncells(groot), 24)
            root = DGG.cellindex(groot, p)
            r = chunk_range(sys, root, level)
            members = collect(r)
            # (a) the complete subtree, root forgotten
            append!(allfails, check_subset_case(sys, level, conn, members;
                label = "complete subtree of $root")); nsubset += 1
            # (b) the same subtree with an interior hole punched out
            if length(members) > 4
                holed = deleteat!(copy(members), cld(length(members), 2))
                append!(allfails, check_subset_case(sys, level, conn, holed;
                    label = "subtree of $root minus its middle cell")); nsubset += 1
            end
            # (c) two disjoint subtrees, so the halo has two components
            if p + 1 <= DGG.ncells(groot)
                r2 = chunk_range(sys, DGG.cellindex(groot, p + 1), level)
                append!(allfails, check_subset_case(sys, level, conn,
                    vcat(collect(r), collect(r2));
                    label = "two adjacent subtrees")); nsubset += 1
            end
            # (d) a scattered subset: every third cell of the subtree
            append!(allfails, check_subset_case(sys, level, conn, members[1:3:end];
                label = "every third cell of $root")); nsubset += 1
            p >= 6 && break
        end
    end
    println("  ", nsubset, " subset cases")

    println()
    println("="^78)
    println("PART 4  laziness: prefix cost must not scale with halo size")
    println("="^78)
    for sys in DGG.systems(), conn in (Vertex(), Edge())
        root = DGG.cellindex(DGG.levelgrid(sys, 0), 1)
        deep = min(DGG.max_level(sys), 9)
        f, s = laziness_failures(sys, root, 2, deep, conn)
        append!(allfails, f)
        @printf("  %-16s %-6s  L2: %6d cells %6d B ctor %6d B prefix | L%-2d: %8d cells %6d B ctor %6d B prefix\n",
            nameof(typeof(sys)), nameof(typeof(conn)), s.ns, s.cs, s.ps,
            deep, s.nd, s.cd, s.pd)
    end

    println()
    println("="^78)
    println("PART 5  argument validation")
    println("="^78)
    for sys in DGG.systems()
        root = DGG.cellindex(DGG.levelgrid(sys, 2), 1)
        for (what, l) in (("level below the cell", 1),
                          ("level above max_level", DGG.max_level(sys) + 1))
            try
                DGG.SubtreeHaloIterator(sys, root, l)
                push!(allfails, Failure(:no_error,
                    string(nameof(typeof(sys))),
                    "SubtreeHaloIterator accepted $what (l=$l) without throwing"))
            catch e
                e isa ArgumentError || push!(allfails, Failure(:wrong_error,
                    string(nameof(typeof(sys))),
                    "$what (l=$l) threw $(typeof(e)), expected ArgumentError"))
            end
        end
    end
    println("  6 systems x 2 invalid levels")

    println()
    println("="^78)
    @printf("TOTAL: %d subtree cases, %d metric comparisons, %d halo cells built, %.1fs\n",
        ncases, nmetriccmp, haloed, time() - t0)
    println("Boundary geometry touched (Vertex/:noise cases): ", coverage)
    println("FAILURES: ", length(allfails))
    for f in allfails[1:min(end, 40)]
        println(f)
    end
    return length(allfails)
end

exit(main() == 0 ? 0 : 1)
