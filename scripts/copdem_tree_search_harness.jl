# Minimal tree-search harness for Copernicus DEM GLO-30 -> IGeo7 level 13.
#
# No DEM values, files, or network access are needed. Conservative regridding
# searches only the source and destination cell trees before it reads values,
# so the exact production tree shapes are enough to reproduce that phase:
#
#   default source       one 1-degree GLO-30 tile, 12,960,000 pixels
#   default destination  one IGeo7 level-6 root refined to level 13
#   multi source         two adjacent GLO-30 tiles
#   multi destination    every level-6 root covering those tiles, refined to level 13
#   execution            serial, below production's outer worker wave
#
# Run the default steady-state measurement with a warm Julia daemon:
#
#   jld --project=benchmark run scripts/copdem_tree_search_harness.jl
#
# Extra modes:
#
#   jld --project=benchmark run scripts/copdem_tree_search_harness.jl multi
#   jld --project=benchmark run scripts/copdem_tree_search_harness.jl profile
#   jld --project=benchmark run scripts/copdem_tree_search_harness.jl allocs
#   jld --project=benchmark run scripts/copdem_tree_search_harness.jl sweep
#
# Pass `repeats=N` or set DGG_TREE_REPEATS to change the default timed search
# passes: three for the small modes, one for the full multi-tile covering. For
# a non-default mode in an existing daemon, set
# `ENV["DGG_TREE_MODE"]` in that daemon before running the script.

import ConservativeRegridding as CR
import DiscreteGlobalGrids as DGG
import Extents
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GlobalRegridding as GR
import Printf: @printf
import Profile

const STI = GO.SpatialTreeInterface

const CONFIG = (
    source_resolution = 30,
    longitude = 10.5,
    latitude = 46.5,
    destination_ancestor = 6,
    destination_level = 13,
    bucket_sizes = (0, 7, 49, 343),
    multi_source_longitudes = (10.5, 11.5),
)

"The exact grids and spaces used by the harness, without building their trees."
function build_spaces()
    source_system = DGG.CopernicusDEMSystem(CONFIG.source_resolution)
    source_tile = DGG.cellat(
        DGG.levelgrid(source_system, 0), CONFIG.longitude, CONFIG.latitude)
    source_grid = DGG.subtree(source_system, source_tile, 1)
    source_space = DGG.DGGSpace(source_grid; chunklevel = 0)

    destination_system = DGG.IGeo7System()
    destination_root = DGG.cellat(
        DGG.levelgrid(destination_system, CONFIG.destination_ancestor),
        CONFIG.longitude, CONFIG.latitude)
    destination_grid = DGG.subtree(
        destination_system, destination_root, CONFIG.destination_level)
    destination_space = DGG.DGGSpace(
        destination_grid; chunklevel = CONFIG.destination_ancestor)

    @assert GR.nchunks(source_space) == 1
    @assert GR.nchunks(destination_space) == 1
    @assert GR.ncells(source_space) == 12_960_000
    @assert GR.ncells(destination_space) == 823_543

    return (; source_tile, source_grid, source_space,
        destination_root, destination_grid, destination_space)
end

"Combine rooted grids without expanding their millions of cell identifiers."
function combine_grids(grids)
    cells = reduce(union, DGG.CellVector.(grids))
    return DGG.PartialGrid(cells)
end

"The nominal one-degree footprint used to select production destination chunks."
function source_tile_extent(source_system, tile)
    south, west = DGG.CopernicusDEM.tilecorner(source_system, tile)
    return Extents.Extent(
        X = (Float64(west), Float64(west) + 1),
        Y = (Float64(south), Float64(south) + 1),
    )
end

"Two adjacent source tiles and every IGeo7 level-6 root covering their extents."
function build_multi_spaces()
    source_system = DGG.CopernicusDEMSystem(CONFIG.source_resolution)
    source_tiles = map(CONFIG.multi_source_longitudes) do longitude
        DGG.cellat(DGG.levelgrid(source_system, 0), longitude, CONFIG.latitude)
    end
    @assert allunique(source_tiles)
    source_grids = map(source_tiles) do tile
        DGG.subtree(source_system, tile, 1)
    end
    source_grid = combine_grids(source_grids)
    source_space = DGG.DGGSpace(source_grid; chunklevel = 0)

    destination_system = DGG.IGeo7System()
    source_extents = map(source_tiles) do tile
        source_tile_extent(source_system, tile)
    end
    root_parts = map(source_extents) do extent
        coverage = DGG.query(destination_system, DGG.MultiOrderCoverage(extent);
            level = CONFIG.destination_ancestor)
        DGG.CellVector(coverage; level = CONFIG.destination_ancestor)
    end
    destination_roots = collect(reduce(union, root_parts))
    destination_grids = map(destination_roots) do root
        DGG.subtree(destination_system, root, CONFIG.destination_level)
    end
    destination_grid = combine_grids(destination_grids)
    destination_space = DGG.DGGSpace(
        destination_grid; chunklevel = CONFIG.destination_ancestor)

    @assert GR.nchunks(source_space) == 2
    @assert GR.ncells(source_space) == 2 * 12_960_000
    @assert length(destination_roots) == 55
    @assert GR.nchunks(destination_space) == length(destination_roots)
    @assert GR.ncells(destination_space) == length(destination_roots) * 7^7

    return (; source_tiles, source_extents, source_grids, source_grid, source_space,
        destination_roots, destination_grids, destination_grid, destination_space)
end

"Build the one source-chunk tree that each production block searches."
function build_source_tree(spaces)
    indices = GR.ownedindices(spaces.source_space, 1)
    return GR.subtree(spaces.source_space, indices)
end

"Build one tree whose packed top layer contains both source tiles."
build_multi_source_tree(spaces) = DGG.treeify(spaces.source_grid)

"Build the production destination tree; IGeo7 caps remain analytical and uncached."
function build_destination_tree(spaces)
    indices = 1:GR.ncells(spaces.destination_space)
    return GR.subtree(spaces.destination_space, indices)
end

destination_cursor(tree::DGG.CapCachedTree) = tree.node
destination_cursor(tree) = tree

# Keep compilation out of the reported tree-build time.
compile_destination_tree(spaces) = (build_destination_tree(spaces); nothing)

"Run the production serial candidate search into a reusable buffer."
function search!(pairs, source_tree, destination_tree)
    empty!(pairs)
    CR.get_all_candidate_pairs!(pairs, GOCore.False(), Extents.intersects,
        source_tree, destination_tree)
    return pairs
end

function first_leaf(tree)
    while !STI.isleaf(tree)
        tree = first(STI.getchild(tree))
    end
    return tree
end

function print_workload(spaces, source_tree, destination_tree, tree_build)
    println("workload")
    @printf("  GLO-%d tile:       %s\n", CONFIG.source_resolution, spaces.source_tile)
    @printf("  source pixels:     %d\n", GR.ncells(spaces.source_space))
    @printf("  IGeo7 root:        %s (level %d -> %d)\n", spaces.destination_root,
        CONFIG.destination_ancestor, CONFIG.destination_level)
    @printf("  destination cells: %d\n", GR.ncells(spaces.destination_space))
    @printf("  Julia threads:     %d (production search is serial per outer worker)\n",
        Threads.nthreads())
    println("  source tree:       ", nameof(typeof(source_tree)))
    cursor = destination_cursor(destination_tree)
    println("  destination tree:  ", nameof(typeof(destination_tree)),
        " (bucket ", cursor.bucket_size, ")")
    println("  destination caps:  ", destination_tree isa DGG.CapCachedTree ?
        "precomputed" : "analytical, constructed per fixed-leaf descent")
    sample_entries = STI.child_indices_extents(first_leaf(destination_tree))
    println("  leaf entries:      ", nameof(typeof(sample_entries)),
        " (sample leaf length ", length(sample_entries), ")")
    @printf("  leaf payload:      %.2f KiB traversal-owned, no retained cache\n",
        sizeof(sample_entries) / 2.0^10)
    @printf("  destination build: %.3f s, %.1f MiB\n",
        tree_build.time, tree_build.bytes / 2.0^20)
end

function print_multi_workload(spaces, source_tree, destination_tree, tree_build)
    println("workload (2 source tiles x their full IGeo7 covering)")
    println("  GLO-30 tiles:      ", join(spaces.source_tiles, ", "))
    println("  tile extents:      ", join(spaces.source_extents, ", "))
    @printf("  source pixels:     %d (%d per tile)\n",
        GR.ncells(spaces.source_space), GR.ncells(spaces.source_space) ÷ 2)
    @printf("  IGeo7 covering:    %d level-%d roots, each refined to level %d\n",
        length(spaces.destination_roots), CONFIG.destination_ancestor,
        CONFIG.destination_level)
    @printf("  destination cells: %d (%d per root)\n",
        GR.ncells(spaces.destination_space), 7^7)
    @printf("  destination area:  whole roots selected by tile coverage (boundary spill kept)\n")
    @printf("  destination levels: %d -> %d\n",
        CONFIG.destination_ancestor, CONFIG.destination_level)
    @printf("  Julia threads:     %d (serial inner search)\n", Threads.nthreads())
    println("  source tree:       ", nameof(typeof(source_tree)), " with ",
        length(source_tree.tree.tiles), " packed tiles")
    cursor = destination_cursor(destination_tree)
    println("  destination tree:  ", nameof(typeof(destination_tree)),
        " (bucket ", cursor.bucket_size, ")")
    println("  destination caps:  ", destination_tree isa DGG.CapCachedTree ?
        "precomputed" : "analytical, constructed per fixed-leaf descent")
    sample_entries = STI.child_indices_extents(first_leaf(destination_tree))
    @printf("  leaf payload:      %.2f KiB traversal-owned, no retained cache\n",
        sizeof(sample_entries) / 2.0^10)
    @printf("  destination build: %.3f s, %.1f MiB\n",
        tree_build.time, tree_build.bytes / 2.0^20)
end

"Summarise how candidate work is distributed over source tiles and covering roots."
function print_candidate_distribution(pairs, spaces)
    nsource = length(spaces.source_tiles)
    ndestination = length(spaces.destination_roots)
    source_stride = GR.ncells(spaces.source_space) ÷ nsource
    destination_stride = GR.ncells(spaces.destination_space) ÷ ndestination

    # `union` stores roots in canonical cell order, which need not match the
    # longitude order above. Match each contiguous index range back to its root.
    function component_order(combined_grid, component_grids, stride)
        combined_cells = DGG.CellVector(combined_grid)
        component_starts = first.(DGG.CellVector.(component_grids))
        return map(1:length(component_grids)) do part
            range_start = combined_cells[1 + (part - 1) * stride]
            something(findfirst(==(range_start), component_starts))
        end
    end
    source_order = component_order(
        spaces.source_grid, spaces.source_grids, source_stride)
    destination_order = component_order(
        spaces.destination_grid, spaces.destination_grids, destination_stride)

    counts = zeros(Int, nsource, ndestination)
    for (source, destination) in pairs
        source_part = cld(source, source_stride)
        destination_part = cld(destination, destination_stride)
        counts[source_part, destination_part] += 1
    end

    println("candidate distribution")
    for row in axes(counts, 1)
        active = count(!iszero, @view counts[row, :])
        @printf("  %-22s %10d candidates across %2d/%d destination roots\n",
            string(spaces.source_tiles[source_order[row]]), sum(@view counts[row, :]),
            active, ndestination)
    end

    root_totals = vec(sum(counts; dims = 1))
    busiest = partialsortperm(root_totals, 1:min(5, ndestination); rev = true)
    println("  busiest destination roots")
    for part in busiest
        root = spaces.destination_roots[destination_order[part]]
        @printf("    %-22s %10d candidates\n", string(root), root_totals[part])
    end
end

function timed_searches!(pairs, source_tree, destination_tree; repeats)
    # The first pass compiles the search and grows `pairs` to production size.
    # Every timed pass reuses that capacity, just like CR's task-local assembly
    # scratch. Its allocated bytes therefore describe traversal, not Vector growth.
    search!(pairs, source_tree, destination_tree)
    expected_count = length(pairs)
    expected_hash = hash(pairs, UInt(0))
    println("steady-state candidate search (candidate buffer retained)")
    @printf("  candidates:        %d\n", expected_count)
    @printf("  order hash:        0x%016x\n", expected_hash)
    for repetition in 1:repeats
        GC.gc()
        measured = @timed search!(pairs, source_tree, destination_tree)
        @assert length(pairs) == expected_count
        @assert hash(pairs, UInt(0)) == expected_hash
        @printf("  run %d:             %.3f s, %.3f MiB allocated, %.3f s GC\n",
            repetition, measured.time, measured.bytes / 2.0^20, measured.gctime)
    end
    return expected_count, expected_hash
end

"Change only the analytical destination leaf size; no cap vector is built."
function destination_tree_with_bucket(spaces, destination_tree, bucket_size)
    return DGG.HierarchicalGridCursor(
        spaces.destination_grid; bucket_size = bucket_size)
end

function sweep_buckets!(pairs, spaces, source_tree, destination_tree,
        expected_count, expected_hash)
    println("destination leaf-size sweep (analytical uncached caps)")
    println("  bucket 0 means descend to individual destination cells")
    for bucket_size in CONFIG.bucket_sizes
        tree = destination_tree_with_bucket(spaces, destination_tree, bucket_size)
        search!(pairs, source_tree, tree) # compile this tree shape and retain capacity
        bucket_count = length(pairs)
        bucket_hash = hash(pairs, UInt(0))
        GC.gc()
        measured = @timed search!(pairs, source_tree, tree)
        @assert length(pairs) == bucket_count
        @assert hash(pairs, UInt(0)) == bucket_hash
        same_default = bucket_count == expected_count && bucket_hash == expected_hash
        @printf("  bucket %3d:         %.3f s, %.3f MiB, %d candidates, default order: %s\n",
            bucket_size, measured.time, measured.bytes / 2.0^20,
            bucket_count, same_default ? "yes" : "no")
    end
end

function cpu_profile!(pairs, source_tree, destination_tree)
    println("CPU profile (flat, inclusive sample count)")
    Profile.clear()
    Profile.@profile search!(pairs, source_tree, destination_tree)
    Profile.print(format = :flat, sortedby = :count, mincount = 20)
end

function allocation_profile!(pairs, source_tree, destination_tree)
    println("allocation profile (0.1% sample, flat, inclusive allocated bytes)")
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=0.001 search!(
        pairs, source_tree, destination_tree)
    results = Profile.Allocs.fetch()
    if isempty(results.allocs)
        println("  no allocations were sampled")
    else
        Profile.print(results; format = :flat, sortedby = :count, mincount = 1)
    end
end

function main(mode::Symbol)
    mode in (:run, :multi, :profile, :allocs, :sweep) || error(
        "mode must be run, multi, profile, allocs, or sweep; got $(repr(mode))")
    multi = mode === :multi
    repeat_arg = findfirst(arg -> startswith(arg, "repeats="), ARGS)
    default_repeats = multi ? "1" : "3"
    repeat_text = repeat_arg === nothing ?
        get(ENV, "DGG_TREE_REPEATS", default_repeats) :
        split(ARGS[repeat_arg], '='; limit = 2)[2]
    repeats = parse(Int, repeat_text)
    repeats > 0 || error("DGG_TREE_REPEATS must be positive")

    spaces = multi ? build_multi_spaces() : build_spaces()
    source_tree = multi ? build_multi_source_tree(spaces) : build_source_tree(spaces)

    compile_destination_tree(spaces)
    GC.gc()
    tree_build = @timed build_destination_tree(spaces)
    destination_tree = tree_build.value
    multi ? print_multi_workload(spaces, source_tree, destination_tree, tree_build) :
            print_workload(spaces, source_tree, destination_tree, tree_build)

    pairs = Tuple{Int,Int}[]
    expected_count, expected_hash = timed_searches!(
        pairs, source_tree, destination_tree; repeats)
    multi && print_candidate_distribution(pairs, spaces)

    mode === :sweep && sweep_buckets!(pairs, spaces, source_tree, destination_tree,
        expected_count, expected_hash)
    mode === :profile && cpu_profile!(pairs, source_tree, destination_tree)
    mode === :allocs && allocation_profile!(pairs, source_tree, destination_tree)
    return nothing
end

# `jld run` adds its own bookkeeping flags to `ARGS`, so select only a mode
# name we recognize. DGG_TREE_MODE is also convenient in wrappers.
mode_names = ("run", "multi", "profile", "allocs", "sweep")
mode_arg = findfirst(in(mode_names), ARGS)
mode = Symbol(mode_arg === nothing ? get(ENV, "DGG_TREE_MODE", "run") : ARGS[mode_arg])
main(mode)
