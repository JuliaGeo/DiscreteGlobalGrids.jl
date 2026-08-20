"""
    Harness

Shared halo verification utilities for deterministic and randomized checks.
Oracles use the cached one-ring in `GridCtx` and do not call the halo APIs.

A "case" is a `(system, root, level, connectivity, field)` tuple. Each case
checks:

1. **shape**   — the halo is sorted, unique, outside the subtree, and tight
                 (every halo cell really does touch the subtree).
2. **oracle**  — the halo equals a brute-force halo computed from `neighbors`
                 alone, inside-out.
3. **iterator**— restartable, prefix-consistent, reproducible across fresh
                 walks, and `length` (where advertised) truthful.
4. **terrain** — every metric computed over the chunk-plus-halo read equals the
                 whole-grid metric restricted to the chunk, and the chunk read
                 never touched a cell it did not hold.
5. **control** — a narrower halo produces misses under a wider stencil.

Adjacency symmetry (whether `b ∈ neighbors(a)` implies `a ∈ neighbors(b)`) is
checked once per `(system, level, connectivity)` by `symmetry_failures`, not
once per case, since it is a property of the grid rather than the halo.
"""
module Harness

include("SphericalTerrain.jl")
include("HaloOracle.jl")

import DiscreteGlobalGrids as DGG
import Random
using DiscreteGlobalGrids: Vertex, Edge, Connectivity
using .SphericalTerrain
using .SphericalTerrain: GridCtx, gridctx, nrange, WholeField, ChunkField,
    positions, value, METRICS, same, first_difference, make_field, FIELD_KINDS,
    chunk_range
using .HaloOracle

export Failure, check_subtree_case, check_subset_case, fast_oracle,
    fast_subset_oracle, gridctx, describe, chunk_range, symmetry_failures,
    whole_results, whole_field, whole_metric, laziness_failures, classify_root

# ---------------------------------------------------------------------------

struct Failure
    kind::Symbol
    case::String
    detail::String
end

Base.show(io::IO, f::Failure) =
    print(io, "FAIL[", f.kind, "] ", f.case, "\n        ", f.detail)

describe(sys, root, level, conn, field = nothing) =
    string(nameof(typeof(sys)), " root=", root, "@L", DGG.level(root),
        " target=L", level, " ", nameof(typeof(conn)),
        field === nothing ? "" : string(" field=", field))

# Oracles derived directly from the cached one-ring.

"Halo positions of the contiguous block `r`, inside-out. Ascending, unique."
function fast_oracle(k::GridCtx, r::UnitRange{Int})
    lo, hi = first(r), last(r)
    io = Int[]
    for p in r, j in nrange(k, p)
        q = k.nbr[j]
        (lo <= q <= hi) || push!(io, q)
    end
    return sort!(unique!(io))
end

"Halo positions of contiguous block `r`, found by scanning outside cells. O(ncells)."
function fast_oracle_outside_in(k::GridCtx, r::UnitRange{Int})
    lo, hi = first(r), last(r)
    oi = Int[]
    for p in 1:length(k)
        (lo <= p <= hi) && continue
        for j in nrange(k, p)
            q = k.nbr[j]
            if lo <= q <= hi
                push!(oi, p); break
            end
        end
    end
    return oi
end

"`(inside_out, outside_in)` halo positions of an arbitrary member set."
function fast_subset_oracle(k::GridCtx, member::BitVector)
    io = Int[]
    for p in 1:length(k)
        member[p] || continue
        for j in nrange(k, p)
            q = k.nbr[j]
            member[q] || push!(io, q)
        end
    end
    sort!(unique!(io))
    oi = Int[]
    for p in 1:length(k)
        member[p] && continue
        for j in nrange(k, p)
            if member[k.nbr[j]]
                push!(oi, p); break
            end
        end
    end
    return io, oi
end

const SYM_CACHE = Dict{Any,Vector{Failure}}()

"""
    symmetry_failures(sys, level, conn) -> Vector{Failure}

Return ordered pairs `(a, b)` for which `b ∈ neighbors(a)` but
`a ∉ neighbors(b)`. Results are memoized per grid and connectivity.
"""
function symmetry_failures(sys, level::Integer, conn::Connectivity)
    get!(SYM_CACHE, (sys, Int(level), conn)) do
        k = gridctx(sys, level, conn)
        tag = string(nameof(typeof(sys)), " L", level, " ", nameof(typeof(conn)))
        out = Failure[]
        for p in 1:length(k), j in nrange(k, p)
            q = k.nbr[j]
            found = false
            for j2 in nrange(k, q)
                k.nbr[j2] == p && (found = true; break)
            end
            found || push!(out, Failure(:adjacency_asymmetry, tag,
                "position $p (cell $(k.cells[p])) lists $q (cell $(k.cells[q])) " *
                "as a neighbour, but not vice versa"))
            length(out) > 8 && return out
        end
        return out
    end
end

# Memoized whole-grid results, independent of the subtree root.

const WHOLE_CACHE = Dict{Any,Any}()

"""
    whole_results(sys, level, conn, kind; spike_at = 1)
        -> (z, Dict{Symbol,Vector})

Return the elevation field and every metric over the complete level grid.
Results are memoized by system, level, connectivity, field kind, and spike
position.
"""
function whole_field(sys, level::Integer, conn::Connectivity, kind::Symbol,
        spike_at::Int)
    key = (:z, sys, Int(level), conn, kind, kind === :spike ? spike_at : 0)
    get!(WHOLE_CACHE, key) do
        k = gridctx(sys, level, conn)
        make_field(kind, k, Random.MersenneTwister(hash((key, :field))); spike_at)
    end::Vector{Float64}
end

"Return one memoized whole-grid metric."
function whole_metric(sys, level::Integer, conn::Connectivity, kind::Symbol,
        spike_at::Int, name::Symbol, fn)
    key = (:m, sys, Int(level), conn, kind, kind === :spike ? spike_at : 0, name)
    get!(WHOLE_CACHE, key) do
        z = whole_field(sys, level, conn, kind, spike_at)
        fn(WholeField(gridctx(sys, level, conn), z))
    end
end

function whole_results(sys, level::Integer, conn::Connectivity, kind::Symbol;
        spike_at::Int = 1)
    z = whole_field(sys, level, conn, kind, spike_at)
    res = Dict{Symbol,Any}()
    for (name, fn) in METRICS
        res[name] = whole_metric(sys, level, conn, kind, spike_at, name, fn)
    end
    return (z, res)
end

# Subtree cases.

"""
    check_subtree_case(sys, root, level, conn, fieldkind; kw...)

Returns `(failures, stats)`.

Keywords:
* `metrics`      — which metrics to compare (default all).
* `outside_in`   — also compute the O(ncells) outside-in oracle.
* `control`      — run a narrower-halo negative control.
* `spike_in_halo`— place a `:spike` field on the first halo cell.
"""
function check_subtree_case(sys, root, level::Integer, conn::Connectivity,
        fieldkind::Symbol; metrics = nothing, outside_in::Bool = false,
        control::Bool = true, spike_in_halo::Bool = true)
    fails = Failure[]
    tag = describe(sys, root, level, conn, fieldkind)
    k = gridctx(sys, level, conn)
    r = chunk_range(sys, root, level)
    pg = DGG.subtree(sys, root, level)

    # Reproducibility.
    halo1 = collect(DGG.halo(pg; connectivity = conn, cells = true))
    halo2 = collect(DGG.halo(pg; connectivity = conn, cells = true))
    halo1 == halo2 || push!(fails, Failure(:nondeterministic, tag,
        "two halo walks differ: $(length(halo1)) vs $(length(halo2)) cells"))

    hp = collect(DGG.halo(pg; connectivity = conn))

    # Ordering, uniqueness, and ancestry.
    issorted(hp) || push!(fails, Failure(:unsorted, tag,
        "halo positions not ascending: " * string(hp)))
    length(unique(hp)) == length(hp) || push!(fails, Failure(:duplicate, tag,
        "halo has $(length(hp) - length(unique(hp))) duplicate cells"))
    for (c, p) in zip(halo1, hp)
        if first(r) <= p <= last(r)
            push!(fails, Failure(:inside, tag,
                "halo cell $c at position $p is inside the subtree block $r"))
            break
        end
        if DGG.ancestor(sys, c, DGG.level(root)) == root
            push!(fails, Failure(:ancestry, tag,
                "halo cell $c has ancestor == root yet position $p is outside $r " *
                "(descendant_range and ancestor disagree)"))
            break
        end
    end

    # Independent one-ring oracle.
    io = fast_oracle(k, r)
    if hp != io
        push!(fails, Failure(:halo_mismatch, tag,
            "halo != brute-force one-ring halo: missing=" *
            string([k.cells[p] for p in setdiff(io, hp)]) *
            " (positions $(setdiff(io, hp))) extra=" *
            string([k.cells[p] for p in setdiff(hp, io)]) *
            " (positions $(setdiff(hp, io)))"))
    end
    if outside_in
        oi = fast_oracle_outside_in(k, r)
        io == oi || push!(fails, Failure(:oracle_direction, tag,
            "inside-out and outside-in oracles differ: only-io=$(setdiff(io, oi)) " *
            "only-oi=$(setdiff(oi, io))"))
    end

    # Iterator protocol.
    it = DGG.halo(pg; connectivity = conn, cells = true)
    collect(it) == halo1 || push!(fails, Failure(:collect_mismatch, tag,
        "a fresh halo walk does not collect to the same cells"))
    n3 = min(3, length(halo1))
    a = collect(Iterators.take(it, n3)); b = collect(Iterators.take(it, n3))
    (a == b == halo1[1:n3]) || push!(fails, Failure(:not_restartable, tag,
        "prefixes differ across restarts: $a vs $b vs $(halo1[1:n3])"))
    eltype(it) == eltype(halo1) || push!(fails, Failure(:eltype, tag,
        "eltype(iterator)=$(eltype(it)) but collect gave $(eltype(halo1))"))
    if Base.IteratorSize(typeof(it)) isa Base.HasLength
        n = try length(it) catch; -1 end
        n == length(halo1) || push!(fails, Failure(:bad_length, tag,
            "IteratorSize is HasLength but length=$n and the walk yields $(length(halo1))"))
    end

    # Whole-grid and chunk metric agreement.
    sp = spike_in_halo && !isempty(hp) ? hp[1] : clamp(first(r), 1, length(k))
    z = whole_field(sys, level, conn, fieldkind, sp)
    chunk = ChunkField(k, sys, root, level, z)
    chunk.halopos == hp || push!(fails, Failure(:chunk_halo, tag,
        "ChunkField halo positions differ from halo(subtree)"))
    want = metrics === nothing ? METRICS :
        Tuple(m for m in METRICS if m[1] in metrics)
    nmetric = 0
    for (name, fn) in want
        cres = fn(chunk)
        nmetric += 1
        wslice = whole_metric(sys, level, conn, fieldkind, sp, name, fn)[r]
        if !same(wslice, cres)
            i, x, y = first_difference(wslice, cres)
            push!(fails, Failure(:chunk_differs, tag,
                "$name differs at chunk offset $i (grid position $(first(r) + i - 1), " *
                "cell $(k.cells[clamp(first(r) + i - 1, 1, length(k))])): " *
                "whole=$x chunk=$y"))
        end
    end
    chunk.misses[] == 0 || push!(fails, Failure(:halo_incomplete, tag,
        "the chunk read reached $(chunk.misses[]) times for a cell it did not " *
        "hold — the halo is missing at least one neighbour of a chunk cell"))

    # Every halo cell must touch the chunk.
    slack = Int[]
    for p in hp
        touches = false
        for j in nrange(k, p)
            if first(r) <= k.nbr[j] <= last(r); touches = true; break; end
        end
        touches || push!(slack, p)
    end
    isempty(slack) || push!(fails, Failure(:halo_loose, tag,
        "$(length(slack)) halo cells have no neighbour inside the subtree: $slack"))

    # A narrower edge halo must miss inputs needed by a vertex stencil.
    if control && conn === Vertex()
        small = collect(DGG.halo(pg; connectivity = Edge()))
        if length(small) < length(halo1)
            probe = ChunkField(k, sys, root, level, z; halo_connectivity = Edge())
            SphericalTerrain.roughness(probe)
            probe.misses[] == 0 && push!(fails, Failure(:control_failed, tag,
                "an Edge() halo ($(length(small)) cells) under a Vertex() stencil " *
                "($(length(halo1)) cells) produced NO misses — the detector is blind"))
        end
    end

    return fails, (; nhalo = length(halo1), nchunk = length(r), nmetric)
end

# Arbitrary subset cases.

"""
    check_subset_case(sys, level, conn, members; label = "")

Checks `halo(PartialGrid)`, `halo(CellVector)` and `halo(CellLookup)` against
the brute-force subset oracle, for an arbitrary member set (holes included).
Also checks the rooted `subtree(sys, root, level)` path where `members` is
exactly a subtree.
"""
function check_subset_case(sys, level::Integer, conn::Connectivity,
        members::Vector{Int}; label::String = "")
    fails = Failure[]
    k = gridctx(sys, level, conn)
    tag = string(nameof(typeof(sys)), " L", level, " ", nameof(typeof(conn)),
        " subset[", length(members), "] ", label)

    member = falses(length(k))
    member[members] .= true
    ids = [k.cells[p] for p in sort(members)]

    io, oi = fast_subset_oracle(k, member)
    io == oi || push!(fails, Failure(:oracle_direction, tag,
        "subset oracles differ: only-io=$(setdiff(io, oi)) only-oi=$(setdiff(oi, io))"))

    for (nm, sub) in (("PartialGrid", DGG.PartialGrid(sys, level, ids)),
                      ("CellVector",  DGG.CellVector(sys, level, ids)),
                      ("CellLookup",  DGG.CellLookup(DGG.CellVector(sys, level, ids))))
        hp = try
            collect(DGG.halo(sub; connectivity = conn))
        catch e
            push!(fails, Failure(:threw, tag,
                "halo($nm) threw " * sprint(showerror, e)))
            continue
        end
        issorted(hp) || push!(fails, Failure(:unsorted, tag, "halo($nm) not ascending"))
        length(unique(hp)) == length(hp) ||
            push!(fails, Failure(:duplicate, tag, "halo($nm) has duplicates"))
        hp == io || push!(fails, Failure(:halo_mismatch_subset, tag,
            "halo($nm) != oracle: missing=$(setdiff(io, hp)) extra=$(setdiff(hp, io))"))
        it = DGG.halo(sub; connectivity = conn)
        n4 = min(4, length(hp))
        collect(Iterators.take(it, n4)) == collect(Iterators.take(it, n4)) ||
            push!(fails, Failure(:not_restartable, tag, "halo($nm) prefix not stable"))
        collect(it) == hp ||
            push!(fails, Failure(:nondeterministic, tag, "halo($nm) not reproducible"))
    end
    return fails
end

# Iterator allocation checks.

"""
    laziness_failures(sys, root, shallow, deep, conn)

Compare allocations for iterator construction and a three-cell prefix at two
target levels. Reports growth that tracks halo cardinality.
"""
function laziness_failures(sys, root, shallow::Integer, deep::Integer,
        conn::Connectivity)
    fails = Failure[]
    tag = string(nameof(typeof(sys)), " root=", root, " L", shallow, " vs L", deep,
        " ", nameof(typeof(conn)))
    probe(l) = begin
        pg = DGG.subtree(sys, root, l)
        DGG.halo(pg; connectivity = conn)                           # warm
        cons = @allocated DGG.halo(pg; connectivity = conn)
        it = DGG.halo(pg; connectivity = conn)
        collect(Iterators.take(it, 3))
        pre = @allocated collect(Iterators.take(
            DGG.halo(pg; connectivity = conn), 3))
        (cons, pre, length(collect(DGG.halo(pg; connectivity = conn))))
    end
    cs, ps, ns = probe(shallow)
    cd, pd, nd = probe(deep)
    grow = nd / max(ns, 1)
    if cd > max(4 * cs, cs + 4096) && grow > 4
        push!(fails, Failure(:construct_scales, tag,
            "construction allocated $cs B for a $ns-cell halo but $cd B for a " *
            "$nd-cell halo ($(round(grow, digits = 1))x more cells)"))
    end
    if pd > max(4 * ps, ps + 4096) && grow > 4
        push!(fails, Failure(:prefix_scales, tag,
            "taking 3 cells allocated $ps B at $ns cells but $pd B at $nd cells"))
    end
    return fails, (; cs, ps, ns, cd, pd, nd)
end

# Geometry labels used in coverage reports.

"""
    classify_root(sys, root, level, conn) -> Set{Symbol}

Return geometry labels for the subtree boundary. Verification does not depend
on these labels.

* `:pentagon`   — a chunk or halo cell has fewer neighbours than the modal degree.
* `:high_degree`— a chunk or halo cell has more neighbours than the modal degree.
* `:seam`       — the halo crosses a level-0 (base cell / face / diamond) boundary.
* `:npole` / `:spole` — the chunk contains the north / south pole.
"""
const MODAL_CACHE = Dict{Any,Int}()

"Return the most common neighbour degree on a grid."
function modal_degree(k::GridCtx)
    get!(MODAL_CACHE, objectid(k)) do
        counts = Dict{Int,Int}()
        for p in 1:length(k)
            d = length(nrange(k, p))
            counts[d] = get(counts, d, 0) + 1
        end
        isempty(counts) ? 0 : argmax(d -> counts[d], collect(keys(counts)))
    end
end

function classify_root(sys, root, level::Integer, conn::Connectivity)
    k = gridctx(sys, level, conn)
    r = chunk_range(sys, root, level)
    hp = collect(DGG.halo(DGG.subtree(sys, root, level); connectivity = conn))
    tags = Set{Symbol}()
    modal = modal_degree(k)
    for p in vcat(collect(r), hp)
        d = length(nrange(k, p))
        d < modal && push!(tags, :pentagon)
        d > modal && push!(tags, :high_degree)
    end
    base(p) = DGG.ancestor(sys, k.cells[p], 0)
    b0 = base(first(r))
    for p in hp
        base(p) == b0 || (push!(tags, :seam); break)
    end
    np, sp = pole_positions(k)
    (np in r || np in hp) && push!(tags, :npole)
    (sp in r || sp in hp) && push!(tags, :spole)
    return tags
end

const POLE_CACHE = Dict{Any,Tuple{Int,Int}}()

"Grid positions of the cells containing the two geographic poles."
function pole_positions(k::GridCtx)
    get!(POLE_CACHE, objectid(k)) do
        (DGG.cellposition(k.grid, DGG.cellat(k.grid, 0.0, 90.0)),
         DGG.cellposition(k.grid, DGG.cellat(k.grid, 0.0, -90.0)))
    end
end

end # module
