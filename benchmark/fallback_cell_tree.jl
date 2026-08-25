# Production-shaped fallback-cell benchmark. It intentionally forces the
# unstructured destination adapter over the same IGeo7 level-5 / global 1°
# spaces used by conservative_block_baseline.jl.
#
#     julia -t 8 --project=benchmark benchmark/fallback_cell_tree.jl
#
# Run in a fresh process: `Sys.maxrss()` is a process-lifetime high-water mark.

import ConservativeRegridding as CR
import DimensionalData as DD
import DiscreteGlobalGrids as DGG
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GlobalRegridding as GR
import SparseArrays

const STI = GO.SpatialTreeInterface

function leafindices!(out, node)
    if STI.isleaf(node)
        append!(out, first(entry) for entry in STI.child_indices_extents(node))
    else
        for child in STI.getchild(node)
            leafindices!(out, child)
        end
    end
    return out
end

function candidate_count(node1, node2)
    count = 0
    STI.dual_depth_first_search(GO.Extents.intersects, node1, node2) do _, _
        count += 1
    end
    return count
end

function main()
    dst = DGG.DGGSpace(DGG.levelgrid(DGG.IGeo7System(), 5))
    raster = DD.DimArray(zeros(Float32, 360, 180),
        (DD.X(-179.5:1.0:179.5), DD.Y(-89.5:1.0:89.5)))
    src = GR.RasterGrid(raster)
    inds = 1:GR.ncells(dst)

    built = @timed GR.CellSpaceRTree(dst, inds)
    fallback = built.value
    native = GR.subtree(dst, inds)
    srctree = GR.celltree(src)

    candidates = @timed CR.get_all_candidate_pairs(
        GOCore.False(), GO.Extents.intersects, srctree, fallback)
    frontier = CR.MultithreadedDualDepthFirstSearch.frontier(
        GO.Extents.intersects, srctree, fallback; nchunks = 64)
    task_candidates = [candidate_count(pair[1], pair[3]) for pair in frontier]

    fallback_weights = @timed CR.intersection_areas(
        GR.manifold(dst), GOCore.True(), fallback, srctree; progress = false)
    native_weights = @timed CR.intersection_areas(
        GR.manifold(dst), GOCore.True(), native, srctree; progress = false)
    Wf, Wn = fallback_weights.value, native_weights.value

    GC.gc()
    maxrss = Int(Sys.maxrss())
    println((
        destination_cells = GR.ncells(dst),
        source_cells = GR.ncells(src),
        original_leaf_indices = sort!(leafindices!(Int[], fallback)) == collect(inds),
        tree_build_seconds = built.time,
        tree_build_bytes = built.bytes,
        candidate_count = length(candidates.value),
        candidate_seconds = candidates.time,
        candidate_bytes = candidates.bytes,
        frontier_tasks = length(frontier),
        task_candidate_min = minimum(task_candidates),
        task_candidate_median = sort(task_candidates)[cld(length(task_candidates), 2)],
        task_candidate_max = maximum(task_candidates),
        task_candidate_mean = sum(task_candidates) / length(task_candidates),
        fallback_weight_seconds = fallback_weights.time,
        fallback_weight_bytes = fallback_weights.bytes,
        native_weight_seconds = native_weights.time,
        native_weight_bytes = native_weights.bytes,
        nonzeros = SparseArrays.nnz(Wf),
        weight_sum = sum(Wf),
        native_identity = Wf == Wn,
        process_maxrss_bytes = maxrss,
        process_maxrss_mib = maxrss / 2^20,
    ))
end

main()
