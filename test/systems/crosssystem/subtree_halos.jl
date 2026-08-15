# ---------------------------------------------------------------------------
# The OUTSIDE face of a subtree boundary.
#
# `subtree_iterators.jl` walks the inside of that boundary — the descendants
# with a neighbour that is not one. This file walks the outside: the level-`l`
# cells that are NOT descendants but have a neighbour that is. Same boundary,
# opposite side, so the two files share a fixture vocabulary and nothing else.
#
# Written against the generic interface wherever it reaches, which is most of
# the file but not all of it. Two departures, both deliberate: `hex_ispentagon`
# names the twelve pentagons through `DGG.H3.ispentagon` and
# `DGG.IGeo7.z7_is_pentagon`, because they have to be named rather than found by
# degree; and `SQUARE_SYSTEMS` / `HEX_SYSTEMS` gate the sections that test a
# SPECIALIZATION, which has nothing to say to a system that does not have one.
# Every law that is system-generic — the defining law, the geometry oracle, the
# sweep bundle, the iterator contract and the allocation laws — runs over
# `systems()` and picks up a system registered later with no edit here.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# THE DESIGN'S VERIFICATION LIST, AND WHERE EACH ROW IS DISCHARGED
#
# The halo design closed with two checklists — "Verification" and "Iterator
# contract". That document is not in this repository, so the rows are RESTATED
# here rather than cited: a mapping pointing at a path a reader cannot open is
# worse than no mapping. Each row names the testsets that answer it and nothing
# more — the reasoning lives at the testset, once. Rows that are NOT covered are
# at the bottom, with the reason.
#
# FROM "Verification":
#
#   Oracles must be forced geometry, not `neighbors`/`subtree_border`, which
#   share indexed topology.  `forced_geometry_halo` on every specialization and
#   on the subset walk — whose `SubsetMembership` calls `neighbors`, so the row
#   binds it too; see `holed_halo_oracle`. `law_halo` is the second and stronger
#   oracle: an ascending position scan of the whole target level, sharing no
#   code with any walk in this package.
#
#   All bundled systems, both connectivities.  Every `for sys in systems()` arm;
#   "$(system) at level $base" runs the whole `check_halo_case` bundle under both
#   and pins `Edge() ⊆ Vertex()`.
#
#   Depth zero and deeper subtrees.  "depth zero is the cell's own one-ring",
#   "forced geometry at depth zero", the sweep's `l in base:base+2`, "$(system)
#   at depth", "the band walk at 64x64", "the seam walk at depth five", "the
#   band walk at max_level and at level 20" (thirty levels of descent on S2),
#   "the directed walk at depth four".
#
#   Ordinary, seam, pole, corner, pentagon, irregular-degree cells.
#   `irregular_cells` by DEGREE, `classify_roots` over whole generations of the
#   square systems, `hex_pentagons` by NAME, "the two HEALPix poles, pinned by
#   location rather than by degree", and the three-face corner block in "the band
#   walk at max_level and at level 20".
#
#   Sortedness, uniqueness, outside ancestry, both adjacency directions.
#   `check_halo_case`, from every sweep; `law_halo` pins SET and ORDER on top.
#
#   Rooted and root-forgotten `PartialGrid`, holes, `CellVector`, `CellLookup`.
#   "$(system): halo on subsets" — delegation asserted BY TYPE, the
#   root-forgotten form and three shapes of punched hole held to the GEOMETRY
#   oracle, and the rooted-but-incomplete form. The other subset rule — no
#   `halo` for a mixed-level `MultiOrderCellSet` — is "a mixed-level set has no
#   halo".
#
#   The ASSUMPTION the subset walk prunes by.  "the coarse-containment law the
#   subset prune rests on" — the one statement about a system that, if it were
#   false, would make `halo` drop a cell silently rather than loudly.
#
#   Equality through `AuthalicSystem`.  "AuthalicSystem forwards the halo walk",
#   over all seven engines and the subset verb on all three containers.
#
#   Restartability and prefix equality.  "the walk is resumable, not restarted",
#   whose closing set assertions cover all seven engines and both wrappers.
#
#   Truthful count guards where `length` exists.  "length is truthful where it
#   exists and absent where it does not" partitions the seven engines into the
#   two that count in closed form and the five that refuse, each of the five
#   pinned by the `MethodError` that IS the contract. TRUTHFULNESS ITSELF IS NOT
#   ASSERTED THERE AND CANNOT BE: `collect` routes through `collect_subtree`,
#   which `error`s on a claim-versus-walk mismatch, so `length(it) ==
#   length(collect(it))` is that call's post-condition — it can throw, it cannot
#   fail. The counts are pinned against a FORMULA in "the band count is closed
#   form", and the guard itself in "collect is the guarded path".
#
#   Construction and short-prefix allocation independent of halo cardinality;
#   scaling measured after warm-up, never by wall clock.  The four arms under
#   "The allocation laws", "a prefix of a deep halo costs O(depth), not
#   O(halo)", and `EagerHaloEngine` — a correct iterator with the same public
#   surface that all four refuse, which is what says they discriminate. No
#   `@elapsed`, no `@time`, no wall-clock threshold anywhere in this file.
#
# FROM "Iterator contract":
#
#   exact, unique, canonically ordered   `law_halo` + `check_halo_case`.
#   restartable and resumable            "the walk is resumable, not restarted".
#   type-stable in `eltype`              "eltype is the system's cell index
#                                        type, on every engine".
#   honest about `IteratorSize`/`length` "length is truthful where it exists
#                                        and absent where it does not".
#   constructible without output-sized   "construction does not allocate in
#     materialization                    proportion to the halo".
#   consumable incrementally, or in      "consumable incrementally and in
#     caller-selected batches            caller-chosen batches".
#
# WHAT IS *NOT* DISCHARGED HERE, AND WHY:
#
#   * The `O(depth)` MEMORY CLASS is not measured directly — nothing here reads
#     a stack depth. Allocation is the proxy: an engine whose state left the
#     inline `SmallList` would show it in the prefix and construction arms.
#
#   * THE GENERIC WALK'S PREFIX COST IS NOT DEPTH-INVARIANT, and no arm pretends
#     it is: `_admit` calls `node_extent` on every node it prunes, and that
#     computes a boundary, which allocates. Measured on IGeo7 from a level-0
#     root, a four-cell generic prefix is 3136 B at `l = 3` and 213216 B at
#     `l = 6` — proportional to the WORK, not the OUTPUT, and only the second is
#     what `SubtreeHaloIterator` promises. The ratio arms hold it to the promise
#     it actually made; the depth-invariance arm is restricted to the engines
#     that have the property.
#
#   * THE SUBSET WALK'S COST is not measured by a clock either, and does not
#     need to be: the walk asks its subset exactly two questions, so a wrapper
#     that counts them counts the work. "the subset walk's cost follows the halo,
#     not the subset" is that count at two sizes, and "the subset walk is lazy"
#     is the construction and prefix pair the other engines get.
#
#   * A SUBSET THAT IS NOT A SUBTREE HAS NO GEOMETRY ORACLE. `ForcedGeometry`
#     answers "does `x` touch a descendant of `c`", and an arbitrary set of ids
#     has no `c` to ask about. Every subset below is a subtree or a subtree with
#     pieces removed, the largest family the oracle reaches; a scattered subset
#     is covered by the container-agreement arms and by nothing stronger.
#
#   * THE BYTE FIGURES QUOTED IN COMMENTS ARE MEASUREMENTS, not bounds. Every
#     assertion below is a difference, a ratio or a `MethodError`; none pins a
#     raw byte count, a fact about one machine and one Julia version.
# ---------------------------------------------------------------------------

module SubtreeHaloTests

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGrids: systems, levelgrid, level, max_level, ncells,
    cellindex, cellposition, neighbors, ancestor, subtree_border, Vertex, Edge,
    SubtreeHaloIterator, subtree_halo
import DiscreteGlobalGrids as DGG

# ---------------------------------------------------------------------------
# The two fixtures that cannot live inside the outer testset
#
# Everything below runs inside ONE outer `@testset`, so a failure anywhere is
# recorded and every later section still runs: a TOP-LEVEL testset throws when
# it finishes with failures, while a nested one only reports upwards. A `struct`
# cannot be declared in that local scope, so the two engines the guard testsets
# need are declared out here.
# ---------------------------------------------------------------------------

# Claims three, yields one — the shape `collect_subtree` exists to catch.
# Without it `collect`'s own `HasLength` route sizes the vector from the claim
# and hands back two `undef` slots as cell ids; the square band walk claims a
# closed-form count, so this guard is load-bearing.
struct MiscountingEngine end
Base.iterate(::MiscountingEngine) = (DGG.LevelIndex(0, 0), 1)
Base.iterate(::MiscountingEngine, ::Int) = nothing
Base.eltype(::Type{MiscountingEngine}) = DGG.LevelIndex
Base.IteratorSize(::Type{MiscountingEngine}) = Base.HasLength()
Base.length(::MiscountingEngine) = 3

# The second fixture, out here for the same reason, playing the opposite role:
# `MiscountingEngine` is an engine every value assertion refuses, and this is
# one every value assertion ACCEPTS. `EagerHaloEngine` answers exactly, in
# canonical order, with a concrete `eltype` and an honest `length`, and is
# resumable, restartable and batchable — it calls `subtree_halo`. What it does
# not do is stay lazy: it materialises the whole halo in its constructor, so
# `length` is free, and again on every fresh walk, which is what an engine that
# built its answer on the first `iterate` would pay. The four allocation laws
# are the only things in the package that refuse it — see "an eager engine with
# the same surface fails every allocation law".
struct EagerHaloEngine{S,C,K}
    system::S
    root::C
    target::Int
    connectivity::K
    cells::Vector{C}
end

EagerHaloEngine(sys, c, l, conn) = EagerHaloEngine(sys, c, Int(l), conn,
    subtree_halo(sys, c, Int(l); connectivity = conn))

Base.iterate(e::EagerHaloEngine{S,C,K}) where {S,C,K} =
    iterate(e, (subtree_halo(e.system, e.root, e.target;
        connectivity = e.connectivity), 1))
Base.iterate(::EagerHaloEngine{S,C,K}, s::Tuple{Vector{C},Int}) where {S,C,K} =
    s[2] > length(s[1]) ? nothing : (@inbounds(s[1][s[2]]), (s[1], s[2] + 1))
Base.eltype(::Type{<:EagerHaloEngine{S,C,K}}) where {S,C,K} = C
Base.IteratorSize(::Type{<:EagerHaloEngine}) = Base.HasLength()
Base.length(e::EagerHaloEngine) = length(e.cells)

# THE THIRD FIXTURE, WHICH MAKES THE SUBSET WALK'S COST COUNTABLE WITHOUT A
# CLOCK. `SubsetMembership` asks the subset it was handed exactly two questions —
# `cellposition` for one cell and `subset_span` for a whole block — so a wrapper
# that forwards both and counts them counts the walk's work exactly. No timing,
# no allocation proxy, no threshold that a faster machine could move.
#
# It is a SUBSET, not an engine: `halo` decides its own container's engine, so
# the counting arm builds `subset_halo_engine` directly and hands it this. The
# answer is asserted equal to the container's on every use, which is what says
# the wrapper changed the measurement and not the walk.
mutable struct CountingSubset{S}
    inner::S
    calls::Int
end

CountingSubset(inner) = CountingSubset{typeof(inner)}(inner, 0)

DGG.cellposition(cs::CountingSubset, c::DGG.AbstractCellIndex) =
    (cs.calls += 1; DGG.cellposition(cs.inner, c))

DGG.Fallbacks.subset_span(cs::CountingSubset, lo::Int, hi::Int) =
    (cs.calls += 1; DGG.Fallbacks.subset_span(cs.inner, lo, hi))

counting_halo(sys, cs::CountingSubset, complete, l, conn) = DGG.collect_subtree(
    DGG.Fallbacks.SubsetHaloIterator(cs, conn,
        DGG.Fallbacks.subset_halo_engine(sys, cs, complete, Int(l), conn)))

# ---------------------------------------------------------------------------
# The allocation harness — also out here, and for a measurement reason
#
# A helper defined inside the testset would be a CLOSURE over that block's
# local scope, and one that captured a boxed local would put bytes on the heap
# belonging to the harness rather than the walk. Nothing below captures
# anything: every input travels as an argument.
#
# CONSTRUCTION IS INSIDE THE MEASUREMENT in every harness here but one: an
# engine that materialised its answer in its constructor would measure zero from
# a harness that built the iterator first and timed only the walk, which is the
# failure the design's last verification row exists to catch. THE ONE EXCEPTION
# is the prefix-against-collect arm of "the subset walk is lazy, and its
# construction is O(1)", which hoists the iterator out on purpose so that the
# WALK's laziness is what is measured — construction is measured separately in
# the same testset, by `subset_construct_bytes`, and what keeps the hoisted arm
# honest is `EagerHaloEngine`, which pays the halo on every walk as well as in
# its constructor and is refused there too.
#
# Every measurement is WARM: the untimed call before the timed one is the
# compile, not a courtesy.
# ---------------------------------------------------------------------------

take_n(it, n::Int) = (seen = 0; for _ in it
    seen += 1
    seen >= n && break
end; seen)

build_and_take(sys, c, l, n::Int) = take_n(SubtreeHaloIterator(sys, c, l), n)

lazy_bytes(sys, c, l, n::Int) =
    (build_and_take(sys, c, l, n); @allocated build_and_take(sys, c, l, n))
eager_bytes(sys, c, l) = (subtree_halo(sys, c, l); @allocated subtree_halo(sys, c, l))

# A sink, so `@allocated` cannot be handed a dead value and told the truth
# about a constructor that never ran.
const SINK = Ref{Any}(nothing)

construct!(sys, c, l) = (SINK[] = SubtreeHaloIterator(sys, c, l); nothing)
construct_bytes(sys, c, l) = (construct!(sys, c, l); @allocated construct!(sys, c, l))

# The same shape for the subset verb, whose argument is a container rather than
# a `(system, cell, level)` triple.
subset_construct!(sub) = (SINK[] = halo(sub); nothing)
subset_construct_bytes(sub) =
    (subset_construct!(sub); @allocated subset_construct!(sub))

# The same three, forced onto the generic outside-first walk. No system reaches
# it through the keyword constructor any more, so it has to be built.
generic_iterator(sys, c, l) = SubtreeHaloIterator(sys, c, Int(l), Vertex(),
    DGG.Fallbacks.generic_halo_engine(sys, c, Int(l), Vertex()))

generic_take(sys, c, l, n::Int) = take_n(generic_iterator(sys, c, l), n)
generic_collect(sys, c, l) = DGG.collect_subtree(generic_iterator(sys, c, l))
generic_construct!(sys, c, l) = (SINK[] = generic_iterator(sys, c, l); nothing)
generic_construct_bytes(sys, c, l) =
    (generic_construct!(sys, c, l); @allocated generic_construct!(sys, c, l))

# And the same harnesses pointed at the eager fixture, so it is measured by the
# same shapes the laws are. (`fixture_` rather than `eager_`: `eager_bytes`
# above is the eager VERB's cost, a different quantity.)
fixture_iterator(sys, c, l) = SubtreeHaloIterator(sys, c, Int(l), Vertex(),
    EagerHaloEngine(sys, c, Int(l), Vertex()))

fixture_collect(sys, c, l) = DGG.collect_subtree(fixture_iterator(sys, c, l))
fixture_construct!(sys, c, l) = (SINK[] = fixture_iterator(sys, c, l); nothing)
fixture_take(sys, c, l, n::Int) = take_n(fixture_iterator(sys, c, l), n)
fixture_construct_bytes(sys, c, l) =
    (fixture_construct!(sys, c, l); @allocated fixture_construct!(sys, c, l))
fixture_prefix_bytes(sys, c, l, n::Int) =
    (fixture_take(sys, c, l, n); @allocated fixture_take(sys, c, l, n))
fixture_collect_bytes(sys, c, l) =
    (fixture_collect(sys, c, l); @allocated fixture_collect(sys, c, l))

# Everything from here down is ONE testset, so a failure in an early section is
# recorded and the rest of the file still runs. See the fixture note above.
@testset "subtree halos" begin

    # -----------------------------------------------------------------------
    # Depth zero: a cell's own one-ring
    # -----------------------------------------------------------------------

    @testset "depth zero is the cell's own one-ring" begin
        for sys in systems()
            grid = levelgrid(sys, 1)
            c = cellindex(grid, 1)
            for conn in (Vertex(), Edge())
                expected = sort!(collect(neighbors(grid, c, 1; connectivity = conn)))
                it = SubtreeHaloIterator(sys, c, 1; connectivity = conn)
                @test collect(it) == expected
                @test subtree_halo(sys, c, 1; connectivity = conn) == expected
                @test eltype(it) == DGG.cellindextype(sys)
            end
        end
    end

    @testset "level validation" begin
        for sys in systems()
            grid = levelgrid(sys, 1)
            c = cellindex(grid, 1)
            @test_throws ArgumentError SubtreeHaloIterator(sys, c, 0)
            @test_throws ArgumentError SubtreeHaloIterator(sys, c, max_level(sys) + 1)
        end
    end

    # -----------------------------------------------------------------------
    # The defining law — the only oracle that tests the ENUMERATION
    # -----------------------------------------------------------------------

    # The law itself, computed the slow honest way: every level-l cell that is
    # not a descendant and has a descendant neighbour. `O(ncells)` and unusable
    # in anger, which is exactly why it is the oracle.
    #
    # WHAT IT BUYS OVER THE GEOMETRY ORACLE. It shares NO CODE with the halo
    # walk — it considers cells by scanning positions 1:ncells, where the walk
    # decides by a pruned depth-first descent — so it pins three things at once:
    # which cells the walk emits, in what order, and that the pruning threw
    # nothing away. Against the GENERIC engine the geometry oracle pins only the
    # third, since both sides run the same descent.
    function law_halo(sys, c, l; connectivity = Vertex())
        grid = levelgrid(sys, l)
        lc = level(c)
        out = DGG.cellindextype(sys)[]
        for p in 1:ncells(grid)
            x = cellindex(grid, p)
            ancestor(sys, x, lc) == c && continue
            any(nb -> ancestor(sys, nb, lc) == c,
                neighbors(grid, x, 1; connectivity)) && push!(out, x)
        end
        return out
    end

    check_law(sys, c, l, conn) =
        @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
              law_halo(sys, c, l; connectivity = conn)

    # Deterministic: no RNG, so a failure names the same cell every run.
    sample_cells(grid, n::Int) =
        [cellindex(grid, i) for i in 1:max(1, ncells(grid) ÷ n):ncells(grid)]

    # The cells whose one-ring is not the modal size — pentagons, face corners,
    # poles — found by DEGREE, so a sweep needs no system knowledge and a system
    # added later is covered without anyone remembering to list its oddities.
    function irregular_cells(grid, limit::Int)
        cells = [cellindex(grid, i) for i in 1:ncells(grid)]
        degrees = [length(neighbors(grid, c, 1)) for c in cells]
        counts = Dict{Int,Int}()
        for d in degrees
            counts[d] = get(counts, d, 0) + 1
        end
        modal = argmax(k -> counts[k], keys(counts))
        odd = [c for (c, d) in zip(cells, degrees) if d != modal]
        return odd[1:min(limit, length(odd))]
    end

    # The budget. `law_halo` is O(ncells) per call, so the sweep SAMPLES ROOTS
    # as it goes deeper rather than dropping the depth: depth is what exposes
    # the cap prune (an under-covering root cap only starts dropping cells once
    # the nodes it prunes are smaller than the subtree's own overhang), so a
    # shallow-only sweep would be the cheap half of the coverage and the useless
    # half.
    #
    # One caveat: on A5 the shipped engine is `ScanHaloEngine`, which enumerates
    # by the same ascending position scan the law does, so there the law checks
    # the adjacency test and the descendant skip but not the enumeration. A5's
    # independent check is the geometry oracle.

    @testset "$(nameof(typeof(sys))): the defining law" for sys in systems()
        grid0 = levelgrid(sys, 0)
        n0 = ncells(grid0)
        mx = max_level(sys)

        # EVERY level-0 root at depth 1, both connectivities. The full generation
        # rather than a sample because it is cheap and it contains every awkward
        # cell at once: all twelve IGeo7 and H3 pentagons, HEALPix's polar faces,
        # S2's six cube faces, ISEA4R's diamonds, A5's dodecahedral roots.
        for i in 1:n0, conn in (Vertex(), Edge())
            check_law(sys, cellindex(grid0, i), 1, conn)
        end

        # A spread of level-0 roots deeper. Depth 2 everywhere; depth 3 on the five
        # systems with sorted subtrees, where the walk is the pruned descent whose
        # prune this is testing (A5 has no `descendant_range`, so it scans and there
        # is no prune to break).
        if mx >= 2
            for c in sample_cells(grid0, 6), conn in (Vertex(), Edge())
                check_law(sys, c, 2, conn)
            end
        end
        if DGG.has_sorted_subtrees(sys) && mx >= 3
            for c in sample_cells(grid0, 3), conn in (Vertex(), Edge())
                check_law(sys, c, 3, conn)
            end
        end
        # Depth 4 from a level-0 root, on the two aperture-7 systems only: the
        # configuration where descendants overhang their parent's drawn polygon
        # by the largest margin, so the one that notices a root cap that has
        # stopped covering them. Not hypothetical — swapping the walk's
        # `rootcap` from `node_extent` to the under-covering `cell_cap` changes
        # no arithmetic, only the covering margin the prune's soundness rests
        # on, and fails 214 assertions of which 4 are here. (The other 210 are
        # in the two directed-walk arms: `forced_geometry_halo` IS the generic
        # walk, so once H3 and IGeo7 grew a directed walk that shares none of
        # it, every comparison of the two became sound-against-broken. The
        # narrow claim this arm still owns is the cap caught with the generic
        # walk on BOTH sides, the only shape that survives the specializations
        # being removed.) Two roots is all the runtime affords: H3's level-4
        # grid is 288k cells and the law visits every one.
        #
        # THE GENERIC WALK IS BUILT EXPLICITLY, not reached — no system's
        # keyword constructor returns `OutsideWalkEngine` any more, so this arm
        # would otherwise have stopped testing the cap prune silently, both
        # engines answering correctly.
        if (sys isa DGG.IGeo7System || sys isa DGG.H3System) && mx >= 4
            for c in sample_cells(grid0, 2), conn in (Vertex(), Edge())
                want = law_halo(sys, c, 4; connectivity = conn)
                @test collect(SubtreeHaloIterator(sys, c, 4; connectivity = conn)) ==
                      want
                @test collect(SubtreeHaloIterator(sys, c, 4, conn,
                    DGG.Fallbacks.generic_halo_engine(sys, c, 4, conn))) == want
            end
        end

        # Roots that are no longer whole faces or whole pentagon fans: sampled and
        # irregular-degree cells one and two levels down, each to depth 2.
        for base in 1:min(2, mx)
            gridb = levelgrid(sys, base)
            roots = unique(vcat(sample_cells(gridb, 4), irregular_cells(gridb, 2)))
            for c in roots, l in (base + 1):min(base + 2, mx), conn in (Vertex(), Edge())
                check_law(sys, c, l, conn)
            end
        end
    end

    # -----------------------------------------------------------------------
    # The independent oracle
    # -----------------------------------------------------------------------

    # The generic walk forced onto unit-sphere boundary comparison, reached through
    # the POSITIONAL constructor so the keyword one keeps choosing whatever engine
    # the system ships.
    forced_geometry_halo(sys, c, l, conn) = DGG.collect_subtree(
        DGG.SubtreeHaloIterator(sys, c, Int(l), conn,
            DGG.Fallbacks.geometry_halo_engine(sys, c, Int(l), conn)))

    # WHAT THE GEOMETRY ORACLE IS WORTH, AND WHERE. Stated once, here, because
    # four later sections lean on it.
    #
    # Against a SPECIALIZED engine it is END-TO-END: the square band walk, the
    # seam merge and the two hexagonal walks enumerate by their own arithmetic —
    # no `_admit`, no descendant range, no cap, and in the square case no
    # neighbour query at all — so the comparison pins ENUMERATION and adjacency
    # PREDICATE together. It is the only oracle here that can see a candidate a
    # specialization never proposed, which is exactly what `neighbors` and
    # `subtree_border` cannot: built from the same index arithmetic the walk is,
    # they agree with a gap for the reason it is there.
    #
    # Against the GENERIC engine it is weaker, and the two testsets below are
    # that case: both sides run the same `OutsideWalkEngine` and only the
    # adjacency provider differs, so what is under test is the predicate alone —
    # whether a drawn boundary and the hierarchy's adjacency agree, at
    # pentagons, poles and cube corners, under both connectivities.
    #
    # The sweep is EVERY level-0 root at depth 1: the only cheap sweep holding
    # every structurally awkward cell at once — twelve IGeo7 and H3 pentagons,
    # HEALPix's polar faces, S2's cube corners, ISEA4R's icosahedral-vertex
    # diamonds, A5's dodecahedral roots. NO EXCLUSION IS NEEDED: two were
    # anticipated, A5's connectivity split and the aperture-7 pentagons, and
    # neither materialised. A future disagreement belongs in this comment,
    # named, with the native indexed walk authoritative — never a silent `skip`.

    # Depth zero is the one case where the geometry provider has no subtree to
    # descend: `root` is its own only target-level descendant. The cursor cannot
    # express that — seeded at the target level it would descend past it to
    # `max_level` and throw on a cell with no children — so the provider answers it
    # directly against the root's own boundary. This pins that it does.
    @testset "forced geometry at depth zero" begin
        for sys in systems()
            for base in 0:1
                grid = levelgrid(sys, base)
                for c in (cellindex(grid, 1), cellindex(grid, ncells(grid))),
                    conn in (Vertex(), Edge())
                    @test forced_geometry_halo(sys, c, base, conn) ==
                          collect(SubtreeHaloIterator(sys, c, base; connectivity = conn))
                end
            end
        end
    end

    @testset "$(nameof(typeof(sys))): geometry agrees with topology" for sys in systems()
        grid0 = levelgrid(sys, 0)
        n0 = ncells(grid0)
        for i in 1:n0, conn in (Vertex(), Edge())
            c = cellindex(grid0, i)
            @test forced_geometry_halo(sys, c, 1, conn) ==
                  collect(SubtreeHaloIterator(sys, c, 1; connectivity = conn))
        end
        if max_level(sys) >= 2
            for i in 1:max(1, n0 ÷ 4):n0, conn in (Vertex(), Edge())
                c = cellindex(grid0, i)
                @test forced_geometry_halo(sys, c, 2, conn) ==
                      collect(SubtreeHaloIterator(sys, c, 2; connectivity = conn))
            end
            # And once more one level down, where a root is no longer a whole face
            # or a whole pentagon fan: sampled level-1 roots plus the ones whose
            # degree marks them as a seam, pole or pentagon child.
            grid1 = levelgrid(sys, 1)
            n1 = ncells(grid1)
            roots1 = unique(vcat([cellindex(grid1, i) for i in 1:max(1, n1 ÷ 4):n1],
                irregular_cells(grid1, 2)))
            for c in roots1, conn in (Vertex(), Edge())
                @test forced_geometry_halo(sys, c, 2, conn) ==
                      collect(SubtreeHaloIterator(sys, c, 2; connectivity = conn))
            end
        end
    end

    # -----------------------------------------------------------------------
    # The sweep harness
    #
    # `law_halo` is the strongest oracle here but `O(ncells)` per call, so it
    # can only be afforded shallow. The bundle below is the other half: a fixed
    # list of the iterator contract's laws, cheap enough to run at every root,
    # level and connectivity a sweep reaches, and the thing every specialization
    # added later is put through unchanged. It is not an oracle — it cannot tell
    # a walk that drops a cell from one that never should have emitted it — but
    # it pins uniqueness, sortedness, outside ancestry and both adjacency
    # directions.
    # -----------------------------------------------------------------------

    # Bases 0, 1 and 2, so a root is a whole face, then a quarter of one, then a
    # sixteenth — the last is the first that can be nowhere flush with its face
    # edge, which is the configuration the square band walk needs. A5 STOPS AT
    # BASE 1: it takes `ScanHaloEngine`, `O(ncells)` per halo with a `Set`-
    # allocating `neighbors`, so base 2 would cost more than the other five
    # systems together and pin nothing they do not. Bases 0 and 1 still run every
    # law at every level, and `law_halo` and the geometry oracle cover A5 at full
    # width.
    sweep_bases(sys) = filter(b -> b <= max_level(sys),
        sys isa DGG.A5System ? (0, 1) : (0, 1, 2))

    sweep_roots(sys, base::Int) = (grid = levelgrid(sys, base);
        unique(vcat(sample_cells(grid, 4), irregular_cells(grid, 2))))

    # The deepest target whose grid is within `budget` times the root generation's,
    # so "deep" means the same amount of WORK on an aperture-4 and an aperture-7
    # system rather than the same number of levels.
    function deep_depth(sys, base::Int, budget::Int = 70_000)
        d = 0
        while base + d + 1 <= max_level(sys) &&
            ncells(levelgrid(sys, base + d + 1)) ÷ ncells(levelgrid(sys, base)) <= budget
            d += 1
        end
        return d
    end

    # The per-case law bundle. Every engine, every specialization, goes through
    # it.
    #
    # THE `collect` IS CAUGHT, and that is not defensive tidiness: it is the
    # GUARDED path, so a count regression arrives as an EXCEPTION rather than a
    # false assertion, and uncaught it would abort the enclosing testset and
    # silently stop about 1500 of the file's assertions from running — a second,
    # unrelated regression masked by the first. Recorded, the sweep finishes and
    # both are reported.
    #
    # `length(it) == length(h)` IS DELIBERATELY ABSENT: it is the post-condition
    # `collect_subtree` enforced on the line above, so it can throw but cannot
    # fail. The counts that ARE pinned are pinned against formulas — see "the
    # band count is closed form".
    function check_halo_case(sys, c, l, conn)
        it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
        h = try
            collect(it)
        catch err
            @test err === nothing     # fails naming the throw, rather than aborting
            return eltype(it)[]
        end
        @test allunique(h)
        lc = level(c)
        @test all(x -> ancestor(sys, x, lc) != c, h)          # outside ancestry
        grid = levelgrid(sys, l)
        @test issorted([cellposition(grid, x) for x in h])    # canonical order
        # Both adjacency directions, under the same connectivity: every border cell
        # reaches the halo, and every halo cell reaches the border.
        border = subtree_border(sys, c, l; connectivity = conn)
        hs, bs = Set(h), Set(border)
        @test all(r -> any(in(hs), neighbors(grid, r, 1; connectivity = conn)), border)
        @test all(x -> any(in(bs), neighbors(grid, x, 1; connectivity = conn)), h)
        return h
    end

    @testset "$(nameof(typeof(sys))) at level $base" for sys in systems(),
            base in sweep_bases(sys)
        for c in sweep_roots(sys, base), l in base:min(base + 2, max_level(sys))
            hv = check_halo_case(sys, c, l, Vertex())
            he = check_halo_case(sys, c, l, Edge())
            # `Edge()` is `Vertex()` minus the cells that touch at a point only.
            @test issubset(Set(he), Set(hv))
        end
    end

    # One deep case per system, where `deep_depth` puts the target grid within a
    # fixed factor of the root generation. One root and one connectivity: the
    # point is that the laws still hold when the halo is hundreds of cells and
    # the descent is long, not to re-sweep at depth. A5 is excluded for the
    # reason `sweep_bases` gives — its `deep_depth` from a level-0 root is
    # level 7, a million-cell scan.
    @testset "$(nameof(typeof(sys))) at depth" for sys in
            filter(s -> !(s isa DGG.A5System), systems())
        d = deep_depth(sys, 0)
        d >= 1 || continue
        check_halo_case(sys, cellindex(levelgrid(sys, 0), 1), d, Vertex())
    end

    # -----------------------------------------------------------------------
    # The square band walk — HEALPix, S2 and ISEA4R away from a face edge
    # -----------------------------------------------------------------------

    # `SquareBandEngine` is the first specialization, so this is the first
    # section where the geometry oracle is the end-to-end one — see the note at
    # `forced_geometry_halo`. It descends the FACE's quadtree in curve order and
    # prunes by lattice overlap, never touching `_admit`, `descendant_range`, a
    # cap or a neighbour query, so the comparison pins the enumeration and the
    # predicate together. `law_halo` runs on one arm below as well, because it
    # is cheaper to be sure than to argue.

    # Which walk a root gets is a fact about the lattice, not something a test
    # should hard-code, so roots are CLASSIFIED by the engine the constructor
    # chose. Three classes, because there are three walks:
    #
    #   `inface`   — `SquareBandEngine` under `NoCheck`: the block is nowhere
    #                flush with its face edge, the band IS the halo, and the
    #                count is closed form.
    #   `seam`     — `SquareBandEngine` under `NativeCheck`: the block is flush
    #                somewhere, the rectangles are a conservative superset, and
    #                every candidate goes through the native one-ring.
    #   `fallback` — anything else, i.e. the generic outside-first walk.
    #
    # The two CLAIMED classes are checked against the oracle below, cell for
    # cell. `fallback` is only asserted EMPTY (`check_root_classes`) — on these
    # three systems it is supposed to have no members at all. The generic walk
    # it names is still oracled: see "the generic fallback still agrees with the
    # oracle".
    function classify_roots(sys, base::Int, l::Int, conn)
        grid = levelgrid(sys, base)
        C = DGG.cellindextype(sys)
        inface, seam, fallback = C[], C[], C[]
        for i in 1:ncells(grid)
            c = cellindex(grid, i)
            e = SubtreeHaloIterator(sys, c, l; connectivity = conn).engine
            if e isa DGG.Fallbacks.SquareBandEngine
                push!(e.check isa DGG.Fallbacks.NoCheck ? inface : seam, c)
            else
                push!(fallback, c)
            end
        end
        return inface, seam, fallback
    end

    # Evenly spaced picks, so a sample of a face-ordered list crosses faces
    # instead of staying on the first one. EVERY SAMPLE OF A CLASSIFIED LIST
    # GOES THROUGH THIS, for a measured reason: `classify_roots` returns its
    # lists in position order, which is face-major, so a raw prefix is one
    # face's cells — at base 3 there are 432 in-face roots on HEALPix, 216 on S2
    # and 360 on ISEA4R over 12, 6 and 10 faces, and `inface[1:6]` is six cells
    # of FACE 0 on all three. `BAND_BASES` includes base 3 to defeat the S2
    # orientation-seed bug, which is a PER-FACE seed, so a face-0-only sample
    # would have retired the base that was added to catch it.
    function spread(v, n::Int)
        isempty(v) && return v
        length(v) <= n && return v
        step = length(v) ÷ n
        return [v[1 + (i - 1) * step] for i in 1:n]
    end

    # HOW MANY ROOTS EACH WALK MUST CLAIM, in closed form, and why a COUNT is
    # needed on top of every oracle comparison. `classify_roots` reads the
    # classes OFF THE CODE UNDER TEST, so a guard that grows STRICTER is
    # invisible to every other assertion: the blocks it stops claiming fall to
    # the seam walk, which probes nothing on a non-flush block and reduces to
    # the same band box filtered by the native one-ring — right answer, slower
    # walk, every oracle comparison still passing element for element. Only a
    # count notices; narrowing HEALPix's interval test to `side <= x0 && x0 +
    # 2side <= n - 1` misroutes 16 blocks per base and fails the two counts.
    #
    # A depth-`d` block at base `b` sits at lattice origin `(ix, iy) · 2^d` with
    # `ix, iy ∈ [0, 2^b)`, and the width-one band fits inside `[0, n-1]²` exactly
    # when `1 <= ix <= 2^b - 2`, likewise `iy` — independent of `d`. So each face
    # contributes `max(0, 2^b - 2)²` in-face blocks, times the level-0
    # generation: 12 faces on HEALPix, 6 on S2, 10 on ISEA4R. At bases 0 and 1
    # that is ZERO, which is why those bases are entirely the seam walk and why
    # they are swept. Every remaining root is the seam walk and NOTHING falls
    # back — the third assertion, the one that fails if a seam configuration is
    # quietly handed to `generic_halo_engine`.
    #
    # (One mutation this CANNOT catch, because it is not one: `2 <= x0 && x0 +
    # side <= n - 2` admits exactly the same blocks, since `x0 = ix << d` with
    # `d >= 1` is even and `n` is a power of two.)
    inface_root_count(sys, base::Int) =
        ncells(levelgrid(sys, 0)) * max(0, (1 << base) - 2)^2

    function check_root_classes(sys, base::Int, inface, seam, fallback)
        total = ncells(levelgrid(sys, base))
        @test length(inface) == inface_root_count(sys, base)
        @test length(seam) == total - inface_root_count(sys, base)
        @test isempty(fallback)
    end

    SQUARE_SYSTEMS = (HEALPixSystem(), S2System(), ISEA4RSystem())

    # BASE 2 IS NOT ENOUGH, and this is measured, not defensive. At base 2
    # exactly four of the sixteen cells per face are non-flush, at lattice
    # `(1,1)`, `(1,2)`, `(2,1)`, `(2,2)` — every one mapped to itself by the
    # square symmetry a wrong S2 orientation seed induces (`SWAP` fixes the two
    # diagonal blocks, `SWAP|INVERT` the two anti-diagonal ones). Seeding the
    # descent with `_hilbert_orientation(c.index, ...)` instead of the face
    # root's `isodd(face) ? SWAP_MASK : 0x0` therefore passes every base-2 arm
    # of this file element for element, on all 144 band cases, while being
    # wrong. At base 3 the non-flush set is 36 of 64 per face and most are not
    # fixed by either symmetry: the same mutation fails 576 of 864 cases.
    BAND_BASES = (2, 3)

    @testset "$(nameof(typeof(sys))): the band walk against forced geometry" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, d in 1:2, conn in (Vertex(), Edge())
            l = base + d
            l <= max_level(sys) || continue
            inface, seam, fallback = classify_roots(sys, base, l, conn)
            # The specialization was reached, and reached on exactly the blocks the
            # lattice says it should be. See `inface_root_count`.
            check_root_classes(sys, base, inface, seam, fallback)
            # Six claimed roots, SPREAD rather than the first six: enough to
            # reach six different faces, but only if the picks are spaced.
            for c in spread(inface, 6)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.SquareBandEngine
                @test it.engine.check isa DGG.Fallbacks.NoCheck
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # One arm against the `O(ncells)` brute force as well: the one check a wrong
    # band and a wrong oracle cannot pass together. Depth 2 only, where a visit
    # to every cell of the target level is still cheap.
    @testset "$(nameof(typeof(sys))): the band walk against the law" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, conn in (Vertex(), Edge())
            l = base + 2
            l <= max_level(sys) || continue
            inface, _, _ = classify_roots(sys, base, l, conn)
            for c in spread(inface, 4)
                @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
                      law_halo(sys, c, l; connectivity = conn)
            end
        end
    end

    # The closed-form count, VERIFIED rather than declared: `4·side + 4` band
    # cells under `Vertex()` and `4·side` under `Edge()`, on every block size
    # the sweep can reach.
    #
    # THIS IS WHERE THE COUNT IS PINNED, and it is pinned against the FORMULA.
    # `collect_subtree` already holds the walk to whatever `length` claims — it
    # `error`s on a mismatch — so comparing a collect against `length` would be
    # comparing the guard against itself and could never fail. The open question
    # is whether the claim is the right NUMBER, and only a formula written out
    # independently of the engine can answer it.
    @testset "$(nameof(typeof(sys))): the band count is closed form" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, d in 1:4
            l = base + d
            l <= max_level(sys) || continue
            side = 1 << d
            for conn in (Vertex(), Edge())
                inface, _, _ = classify_roots(sys, base, l, conn)
                isempty(inface) && continue
                for c in spread(inface, 3)
                    it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                    @test Base.IteratorSize(typeof(it)) isa Base.HasLength
                    @test length(it) == (conn isa Vertex ? 4side + 4 : 4side)
                end
            end
        end
    end

    # A 64x64 block, whose band is 260 cells, element for element against the
    # oracle on both connectivities. The shallow cases above all have
    # `side <= 8`, where a wrong `_restore_code` on the way back up the face
    # descent can still land on the right cell by accident; at nine levels of
    # descent it cannot. Three roots from base 3 rather than one from base 2,
    # for `BAND_BASES`' reason: a base-2 block is symmetric under the very
    # transforms a wrong descent applies.
    @testset "the band walk at 64x64" begin
        sys = S2System()
        l = 3 + 6
        inface, _, _ = classify_roots(sys, 3, l, Vertex())
        @test !isempty(inface)
        for c in spread(inface, 3), conn in (Vertex(), Edge())
            it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.SquareBandEngine
            @test length(it) == (conn isa Vertex ? 260 : 256)
            @test collect(it) == forced_geometry_halo(sys, c, l, conn)
            # The contract bundle at nine levels of descent too: sortedness and
            # both adjacency directions are where a deep walk can go wrong
            # differently from a shallow one.
            check_halo_case(sys, c, l, conn)
        end
    end

    # The corner law. `Edge()` drops exactly the four cells that touch the block
    # at a vertex only, so the two halos differ by four and nothing else. Pinned
    # to a NON-FLUSH block deliberately: a flush block can lose a corner across
    # a cube corner, where three faces meet and the diagonal neighbour does not
    # exist, so the count would be three there.
    @testset "Edge drops exactly the four diagonal corners" begin
        sys = S2System()
        inface, _, _ = classify_roots(sys, 2, 4, Vertex())
        # Four blocks on four different faces, not the first four of face 0:
        # `_band_corner` works in lattice coordinates and the four in-face
        # blocks of one S2 face are exactly the set a wrong orientation seed
        # maps to itself.
        for c in spread(inface, 4)
            hv = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Vertex()))
            he = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Edge()))
            @test length(setdiff(Set(hv), Set(he))) == 4
            @test issubset(Set(he), Set(hv))
        end
    end

    # -----------------------------------------------------------------------
    # The seam walk — the same engine where the block touches a face edge
    # -----------------------------------------------------------------------

    # The in-face band is exact by construction. The seam band is derived only
    # as a COVERING superset filtered by the native one-ring, and that admits
    # THREE ways to be wrong where the in-face path has one:
    #
    #   * yielding a non-halo cell — loud, the filter catches it;
    #   * MISSING one — silent, and the reason every arm here goes through the
    #     geometry oracle: the rectangles never propose it, the filter never
    #     sees it, and `neighbors` and `subtree_border` agree with the gap
    #     because the missing cell is missing from the same index arithmetic
    #     they are built from;
    #   * being sound, covering and SLACK — answers correctly, costs more, and
    #     no oracle can see it. `band_candidate_count` below is for that one.
    #
    # BASES 0 AND 1 ARE THE POINT: a level-0 block is a whole face, flush on all
    # four sides; a level-1 block is flush on one side per axis with a face
    # corner of its own. Those are the configurations the in-face guard used to
    # send to the generic walk, so a quiet re-routing back fails here rather
    # than passing everywhere. DEPTH 3 because the monotonicity argument is
    # about the INTERIOR rim cells of a flush side and only two cells of each
    # side are ever probed: at depth 1 a side is two cells and both are probed,
    # at depth 3 it is two of eight and six cells of every flush side reach
    # faces no probe ever asked about.
    seam_roots(sys, base::Int, seam) = unique(vcat(spread(seam, 6),
        filter(in(Set(seam)), irregular_cells(levelgrid(sys, base), 4))))

    # THE CANDIDATE STREAM BEFORE THE CHECK: the engine's own rectangles under
    # `NoCheck` with the corner rule off, which is what "the derived band"
    # names. Counted by iterating and not by `collect`, because `NoCheck`
    # declares `HasLength` and `collect` would size its vector from the in-face
    # perimeter formula, which a seam rectangle list does not obey.
    function band_candidate_count(e)
        n = 0
        for _ in DGG.Fallbacks.SquareBandEngine(e.curve, DGG.Fallbacks.NoCheck(),
                e.level, e.faceside, e.homeface, e.x0, e.y0, e.side, true, e.rects)
            n += 1
        end
        return n
    end

    # THE SEAM RECTANGLES ARE TIGHT, AND NOTHING ELSE PINS THAT. Widening every
    # `FaceRect` by one cell is BEHAVIOUR-PRESERVING — the surplus candidates
    # are not halo cells, `NativeCheck` rejects them, and every oracle
    # comparison here still passes element for element while the walk visits 32%
    # more candidates. So the derivation could be replaced by a lazy bounding
    # box and the file would stay green through a silent performance regression.
    # It is the mirror of "the check filters a widened arc": that arm pins the
    # CHECK by widening the band, this one pins the BAND by counting it.
    #
    # THE LAW NEEDS NO THRESHOLD because the derived band is exact, not merely
    # covering: every rectangle is the bounding box of probe images, and a probe
    # image is a neighbour of a block rim cell on another face — a `Vertex()`
    # halo cell by definition — so the candidates ARE the `Vertex()` halo, cell
    # for cell. `_seam_probe` asks for `Vertex()` whatever was requested, so one
    # number serves both arms; under `Edge()` the stream is still the `Vertex()`
    # halo and `NativeCheck` removes the vertex-only contacts. Measured over 600
    # (system, base, depth, connectivity) cases, including a whole-face block at
    # depth 9 and a corner block at `max_level`: zero surplus anywhere.
    @testset "$(nameof(typeof(sys))): the seam walk against forced geometry" for sys in
            SQUARE_SYSTEMS
        for base in (0, 1, 2, 3), d in 1:3, conn in (Vertex(), Edge())
            l = base + d
            l <= max_level(sys) || continue
            inface, seam, fallback = classify_roots(sys, base, l, conn)
            check_root_classes(sys, base, inface, seam, fallback)
            @test !isempty(seam)                 # the seam path was reached
            for c in seam_roots(sys, base, seam)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.SquareBandEngine
                @test it.engine.check isa DGG.Fallbacks.NativeCheck
                h = collect(it)
                @test h == forced_geometry_halo(sys, c, l, conn)
                @test band_candidate_count(it.engine) ==
                      (conn isa Vertex ? length(h) : length(subtree_halo(sys, c, l)))
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # And one arm against the `O(ncells)` brute force. Bases 0 and 1 at depth 2,
    # where the target grids are at most 3072 cells on all three systems.
    @testset "$(nameof(typeof(sys))): the seam walk against the law" for sys in
            SQUARE_SYSTEMS
        for base in (0, 1), conn in (Vertex(), Edge())
            l = base + 2
            l <= max_level(sys) || continue
            _, seam, _ = classify_roots(sys, base, l, conn)
            for c in spread(seam, 4)
                @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
                      law_halo(sys, c, l; connectivity = conn)
            end
        end
    end

    # The count contract, in the negative. No perimeter formula survives a seam —
    # a cube corner is three cells where the in-face rule wants four, an ISEA4R
    # icosahedral vertex is five — so the seam engine declares `SizeUnknown()`
    # and defines NO `length`. The `MethodError` is the contract; a `length` that
    # silently walked the halo to answer would be the thing the design forbids.
    @testset "the seam walk declares no length" begin
        for sys in SQUARE_SYSTEMS, base in (0, 1)
            l = base + 2
            l <= max_level(sys) || continue
            _, seam, _ = classify_roots(sys, base, l, Vertex())
            isempty(seam) && continue
            it = SubtreeHaloIterator(sys, first(seam), l)
            @test Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
            @test_throws MethodError length(it)
        end
    end

    # Deeper than the sweep can afford at every root, on one root per system: a
    # 32x32 whole-face block from a level-0 root, where each of the four flush
    # sides is thirty-two cells long and only two were ever probed. If the
    # monotonicity argument were wrong anywhere, the gap would be widest here.
    @testset "the seam walk at depth five" begin
        for sys in SQUARE_SYSTEMS, conn in (Vertex(), Edge())
            c = cellindex(levelgrid(sys, 0), 1)
            it = SubtreeHaloIterator(sys, c, 5; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.SquareBandEngine
            @test collect(it) == forced_geometry_halo(sys, c, 5, conn)
            # And still tight thirty-two cells along a flush side, which is
            # where a bound taken lazily would have the most room to be slack.
            @test band_candidate_count(it.engine) == length(subtree_halo(sys, c, 5))
            check_halo_case(sys, c, 5, conn)
        end
    end

    # THE DEEP REGIME. Everything above tops out at level 9. Two things can only
    # go wrong deep: `FaceRect` stores its bounds as `Int32`, which holds a
    # lattice coordinate through LEVEL 31 and overflows at 32 — so S2's
    # `max_level` of 30 has exactly one level of headroom and a bump must be
    # evaluated against 31, not against 30 and not against `_SQUARE_CAP` (see
    # `FaceRect`'s docstring); and the face-quadtree descent is 30 levels long
    # here rather than nine, so a `code`/`x`/`y` restore off by a level has
    # thirty chances to show rather than nine.
    #
    # A MAX-LEVEL CORNER BLOCK is the sharpest cheap case: root at
    # `max_level - 1`, target `max_level`, so the block is 2x2, flush on two
    # sides, and its corner is a face corner — three faces meet there on S2 and
    # the diagonal candidate does not exist, which is why the counts differ by
    # system. The level-20 arm is the same shape one decade shallower, so a
    # failure about the DEPTH and not about the corner separates the two.
    @testset "the band walk at max_level and at level 20" begin
        for sys in SQUARE_SYSTEMS
            mx = max_level(sys)
            for base in (mx - 1, 19), conn in (Vertex(), Edge())
                l = base + 1
                # Position 1 is lattice (0, 0) of face 0 under both curves: the
                # Morton systems because min-code is min-corner, S2 because the
                # Hilbert curve enters face 0 at its origin.
                c = cellindex(levelgrid(sys, base), 1)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.SquareBandEngine
                @test it.engine.check isa DGG.Fallbacks.NativeCheck
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                # And tight at level 30, where the rectangles' `Int32` bounds
                # are one level from overflowing: a bound derived a level too
                # late would be wide here and nowhere else.
                @test band_candidate_count(it.engine) == length(subtree_halo(sys, c, l))
                check_halo_case(sys, c, l, conn)
            end
        end
        # The counts, pinned so a walk that agreed with a wrong oracle would
        # still have to explain itself. A 2x2 corner block has five in-face band
        # cells; the rest come across the two seams, and the diagonal one exists
        # on HEALPix but not at an S2 cube corner or an ISEA4R icosahedral
        # vertex. `Edge()` drops all four diagonal contacts everywhere.
        for (sys, nv) in ((HEALPixSystem(), 12), (S2System(), 11),
                          (ISEA4RSystem(), 11))
            mx = max_level(sys)
            c = cellindex(levelgrid(sys, mx - 1), 1)
            @test length(collect(SubtreeHaloIterator(sys, c, mx))) == nv
            @test length(collect(SubtreeHaloIterator(sys, c, mx;
                connectivity = Edge()))) == 8
        end
    end

    # -----------------------------------------------------------------------
    # The calibrated directed walk — H3 and IGeo7
    # -----------------------------------------------------------------------

    # `HexArcHaloEngine` seeds each NEIGHBOUR's rim automaton with a calibrated
    # arc and walks it down, sharing no step with the generic descent and never
    # asking about the root's own subtree, so the geometry oracle is again the
    # end-to-end one. What must NOT be the only oracle here is `subtree_border`:
    # the seeded automaton IS the border automaton, so the two would agree with
    # a wrong arc for precisely the reason it was wrong.
    #
    # THE TWO STRUCTURAL RISKS THIS SECTION EXISTS FOR.
    #
    #   * Pentagons. Twelve per level per system, and where the arc-3
    #     calibration lives: the only configuration in 52,182 measured pairs
    #     whose minimal covering arc is three directions wide is a pentagon
    #     neighbour whose deleted direction falls between the two touching
    #     children. Pinned BY NAME (`ispentagon`, `z7_is_pentagon`) rather than
    #     found by degree, so a change to `irregular_cells` cannot quietly stop
    #     covering them, and the count is asserted to be twelve at every base.
    #   * Parity. The two systems' automata have exchanged parity branches and
    #     test their `L < 6` guards in the opposite order, so a seed right on H3
    #     at an even level can be wrong on IGeo7 at the same level. Every arm
    #     therefore runs BOTH systems at consecutive root levels.

    HEX_SYSTEMS = (H3System(), IGeo7System())

    hex_ispentagon(sys, c) = sys isa DGG.H3System ?
        DGG.H3.ispentagon(c) : DGG.IGeo7.z7_is_pentagon(c.id)

    # The twelve pentagons of a level, by NAME. Deep levels have billions of
    # cells and must never be enumerated: the centre child of a pentagon is a
    # pentagon, so descending the first child from each level-0 pentagon names
    # the whole level-`base` pentagon set in twelve short walks.
    function hex_pentagons(sys, base::Int)
        grid0 = levelgrid(sys, 0)
        out = DGG.cellindextype(sys)[]
        for i in 1:ncells(grid0)
            c = cellindex(grid0, i)
            hex_ispentagon(sys, c) || continue
            for _ in 1:base
                c = first(DGG.children(sys, c))
            end
            hex_ispentagon(sys, c) && push!(out, c)
        end
        return out
    end

    # All twelve pentagons, the ring around two of them (the arc-3 neighbours),
    # and a spread of ordinary cells by position.
    function hex_roots(sys, base::Int, nhex::Int)
        grid = levelgrid(sys, base)
        pents = hex_pentagons(sys, base)
        @test length(pents) == 12
        n = ncells(grid)
        step = max(1, n ÷ nhex)
        rest = [cellindex(grid, i) for i in 1:step:n][1:min(nhex, end)]
        around = unique(vcat([collect(neighbors(grid, p, 1)) for p in pents[1:2]]...))
        return unique(vcat(pents, around, rest))
    end

    # Which walk a root gets is read OFF THE CODE UNDER TEST, exactly as
    # `classify_roots` does for the square systems. Three classes:
    #
    #   `child`    — `HexChildHaloEngine`: `target == level(root) + 1`, where the
    #                calibration is already the answer and no automaton runs.
    #   `arc`      — `HexArcHaloEngine`: one seeded rim automaton per neighbour.
    #   `fallback` — anything else, i.e. the generic outside-first walk.
    function classify_hex_roots(sys, roots, l::Int, conn)
        C = DGG.cellindextype(sys)
        child, arc, fallback = C[], C[], C[]
        for c in roots
            e = SubtreeHaloIterator(sys, c, l; connectivity = conn).engine
            if e isa DGG.Fallbacks.HexChildHaloEngine
                push!(child, c)
            elseif e isa DGG.Fallbacks.HexArcHaloEngine
                push!(arc, c)
            else
                push!(fallback, c)
            end
        end
        return child, arc, fallback
    end

    # HOW MANY ROOTS EACH WALK MUST CLAIM — the same argument
    # `inface_root_count` makes for the square systems. `hex_halo_engine` has
    # four guards that send a case back to the generic walk and none was ever
    # observed to fire, so a guard that grew stricter (a calibration that started
    # rejecting arc-3 pentagons, say) is INVISIBLE to every oracle arm: the
    # rejected roots fall to the generic walk, which answers correctly. Only a
    # count notices. The rule is not statistical — depth one is
    # `HexChildHaloEngine`, every deeper target is `HexArcHaloEngine`, for EVERY
    # root of both systems at every level — and it is asserted exhaustively over
    # the level-0 and level-1 generations, 122 + 842 H3 roots and 12 + 72 IGeo7.
    function check_hex_classes(sys, roots, base::Int, l::Int, conn)
        child, arc, fallback = classify_hex_roots(sys, roots, l, conn)
        if l == base + 1
            @test length(child) == length(roots)
            @test isempty(arc)
        else
            @test isempty(child)
            @test length(arc) == length(roots)
        end
        @test isempty(fallback)
    end

    @testset "$(nameof(typeof(sys))): every root takes the directed walk" for sys in
            HEX_SYSTEMS
        for base in (0, 1), conn in (Vertex(), Edge())
            grid = levelgrid(sys, base)
            roots = [cellindex(grid, i) for i in 1:ncells(grid)]
            for l in (base + 1):min(base + 4, max_level(sys))
                check_hex_classes(sys, roots, base, l, conn)
            end
        end
        # AND ONE DEEP BASE, because the guards are not all base-independent:
        # `_hex_validate` runs only from depth three, `_hex_calibrate` reads a
        # ring that is a whole base cell's at base 0 and an ordinary hexagon's
        # at base 8, and the seeded frames sit at the other parity. The
        # generation cannot be enumerated — H3's level-8 grid is 7e8 cells — so
        # the roots are the twelve pentagons BY NAME, the ring around two of
        # them, and a positional spread. Depth 4 is in every base because it is
        # the first at which a seeded arc has been through three transitions.
        base = 8
        if base + 1 <= max_level(sys)
            roots = hex_roots(sys, base, 4)
            for conn in (Vertex(), Edge()),
                    l in (base + 1):min(base + 4, max_level(sys))
                check_hex_classes(sys, roots, base, l, conn)
            end
        end
    end

    # The differential sweep. Bases 0 and 1 are whole base cells and their
    # children; bases 5 and 8 are ordinary cells deep in the hierarchy, where
    # the ring is six neighbours of an ordinary hexagon and the seeded frames
    # sit at both parities. Depths 1-3 everywhere: depth 1 is the automaton-free
    # path, depth 2 is the first seeded walk and runs WITHOUT the depth-two
    # validation, and depth 3 is the first depth that validates.
    @testset "$(nameof(typeof(sys))): the directed walk against forced geometry" for
            sys in HEX_SYSTEMS
        for base in (0, 1, 5, 8)
            base + 1 <= max_level(sys) || continue
            roots = hex_roots(sys, base, 4)
            for c in roots, d in 1:3, conn in (Vertex(), Edge())
                l = base + d
                l <= max_level(sys) || continue
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
            end
        end
    end

    # THE EXACTNESS CONTRACT — the hexagonal mirror of the seam band's tightness
    # arm, and the reverse of it. There the band is counted because the check
    # would hide a slack one; here the band is WIDENED because the check would
    # otherwise never be exercised at all: the calibrated walk is exact
    # (candidate-to-halo ratio 1.0000 at every depth), so deleting both
    # `_touches_subtree` calls in `HexArcHaloEngine` leaves this file at its full
    # pass count with nothing red.
    #
    # Widening every calibrated arc `(L, s)` to `(L + 1, s - 1)` keeps the
    # original arc inside the new one, so the widened walk is a strict SUPERSET —
    # asserted by counting the raw automaton output, an inequality and not a
    # pinned number — and the engine's answer must be unchanged, because the
    # check filters the surplus back out. Deliberately NOT a count of the
    # surplus: the invariant is "the emitted set is the halo whatever the band
    # proposes", and a number would pin the widening rather than the check.
    function widen_hex_arcs(e)
        ring = e.ring
        for i in 1:length(ring)
            h = ring[i]
            ring = DGG.Helpers.small_setindex(ring,
                DGG.Fallbacks.HexNeighbour(h.cell, h.lo, h.arclen + Int8(1),
                    Int8(mod(Int(h.start) - 1, 6))), i)
        end
        return DGG.Fallbacks.HexArcHaloEngine(e.system, e.grid, e.root,
            e.rootlevel, e.target, e.connectivity, ring)
    end

    # The engine's candidate stream BEFORE the check: the seeded automata alone,
    # which is what "conservative band" names.
    function hex_candidate_count(e)
        n = 0
        for i in 1:length(e.ring)
            nb = e.ring[i]
            for _ in DGG.seeded_rim_engine(e.system, nb.cell, e.target,
                    Int(nb.arclen), Int(nb.start))
                n += 1
            end
        end
        return n
    end

    @testset "$(nameof(typeof(sys))): the check filters a widened arc" for sys in
            HEX_SYSTEMS
        for base in (0, 1, 2), d in 2:3, conn in (Vertex(), Edge())
            l = base + d
            l <= max_level(sys) || continue
            for c in spread(hex_roots(sys, base, 3), 5)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.HexArcHaloEngine
                # Recorded and skipped rather than thrown: a calibration
                # regression sends the root to the generic walk, and
                # `widen_hex_arcs` would then throw a `FieldError` outside a
                # `@test` and abort this testset — measured, 314 of the file's
                # assertions stop running. The failure above is the report.
                it.engine isa DGG.Fallbacks.HexArcHaloEngine || continue
                wide = widen_hex_arcs(it.engine)
                # The band really did widen, so the equality below is the check
                # doing work and not the mutation being a no-op.
                @test hex_candidate_count(wide) > hex_candidate_count(it.engine)
                @test collect(SubtreeHaloIterator(sys, c, l, conn, wide)) ==
                      collect(it)
            end
        end
    end

    # The contract bundle on a smaller spread, because it builds two `Set`s and
    # a `subtree_border` per case. Pentagons first.
    @testset "$(nameof(typeof(sys))): the directed walk keeps the contract" for sys in
            HEX_SYSTEMS
        for base in (0, 2), conn in (Vertex(), Edge())
            base + 1 <= max_level(sys) || continue
            for c in spread(hex_roots(sys, base, 2), 6),
                    l in (base + 1):min(base + 3, max_level(sys))
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # Depth four, where the halo is 246 cells around a hexagon and 205 around a
    # pentagon, and a seeded arc has been through three transitions rather than
    # one. The oracle costs about 50 ms a call, so this is two roots — a
    # pentagon and a hexagon — per system rather than a sweep.
    @testset "$(nameof(typeof(sys))): the directed walk at depth four" for sys in
            HEX_SYSTEMS
        grid = levelgrid(sys, 2)
        pent = first(hex_pentagons(sys, 2))
        hex = first(filter(c -> !hex_ispentagon(sys, c),
            collect(neighbors(grid, pent, 1))))
        for c in (pent, hex), conn in (Vertex(), Edge())
            it = SubtreeHaloIterator(sys, c, 6; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.HexArcHaloEngine
            @test collect(it) == forced_geometry_halo(sys, c, 6, conn)
        end
    end

    # One arm against the `O(ncells)` brute force: level-0 roots at depths 1 and
    # 2, where the target grids are at most 5882 cells on H3 and 492 on IGeo7.
    @testset "$(nameof(typeof(sys))): the directed walk against the law" for sys in
            HEX_SYSTEMS
        for c in spread(hex_roots(sys, 0, 4), 6), d in 1:2, conn in (Vertex(), Edge())
            @test collect(SubtreeHaloIterator(sys, c, d; connectivity = conn)) ==
                  law_halo(sys, c, d; connectivity = conn)
        end
    end

    # THE COUNTS, PINNED BUT NOT PROMISED. `3^(d+1) + 3` around a hexagon and
    # `5(3^d + 1)/2` around a pentagon, from a per-neighbour census of
    # `(3^d + 1)/2` that holds for both arc lengths — but derived by ENUMERATION
    # rather than from the transition recurrence, so `IteratorSize` still says
    # `SizeUnknown()` with no `length` method. Pinning is not promising: this
    # arm turns "the census changed" into a failure instead of a silent drift.
    @testset "$(nameof(typeof(sys))): the directed walk's census" for sys in HEX_SYSTEMS
        for base in (0, 3)
            grid = levelgrid(sys, base)
            pents = hex_pentagons(sys, base)
            hexes = filter(c -> !hex_ispentagon(sys, c),
                collect(neighbors(grid, first(pents), 1)))
            # IGeo7's base tessellation is twelve pentagons and nothing else, so
            # at base 0 there is no hexagon to take; base 3 supplies both.
            cells = isempty(hexes) ? [first(pents)] : [first(pents), first(hexes)]
            for d in 1:4
                base + d <= max_level(sys) || continue
                for c in cells
                    n = length(collect(SubtreeHaloIterator(sys, c, base + d)))
                    @test n == (hex_ispentagon(sys, c) ? (5 * (3^d + 1)) ÷ 2 :
                                3^(d + 1) + 3)
                end
            end
        end
    end

    # The count contract in the negative, as for the seam walk: the census above
    # is evidence, not an API promise, so neither hex engine declares a length.
    @testset "the directed walk declares no length" begin
        for sys in HEX_SYSTEMS, d in 1:2
            c = cellindex(levelgrid(sys, 1), 1)
            it = SubtreeHaloIterator(sys, c, 1 + d)
            @test Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
            @test_throws MethodError length(it)
        end
    end

    # THE LAZINESS LAW. The generic walk descends from `rootcells` and prunes by
    # cap, so reaching the first halo cell of a deep target costs a traversal
    # that grows with the target level: ten cells of an IGeo7 L=5, d=7 halo cost
    # 42 ms and 779 KB. The directed walk holds one seeded automaton and its
    # frame stack, both isbits, so the same prefix is a flat 256 bytes at every
    # depth. The assertion that matters is INDEPENDENCE, not a threshold: the
    # same byte count at depths three, five and seven IS the `O(depth)` claim.
    @testset "a prefix of a deep halo costs O(depth), not O(halo)" begin
        prefix10(it) = take_n(it, 10)
        for sys in HEX_SYSTEMS
            base = 5
            grid = levelgrid(sys, base)
            root = cellindex(grid, ncells(grid) ÷ 2 + 1)
            depths = filter(d -> base + d <= max_level(sys), [3, 5, 7])
            allocs = map(depths) do d
                prefix10(SubtreeHaloIterator(sys, root, base + d))     # compile
                @allocated prefix10(SubtreeHaloIterator(sys, root, base + d))
            end
            @test all(a -> a <= 4096, allocs)
            @test length(unique(allocs)) == 1
            @test all(d -> prefix10(SubtreeHaloIterator(sys, root, base + d)) == 10,
                depths)
        end
        # And the walk it replaced, on the system where the violation was
        # measured: the generic prefix allocates hundreds of kilobytes where the
        # directed one allocates hundreds of bytes.
        sys = IGeo7System()
        grid = levelgrid(sys, 5)
        root = cellindex(grid, ncells(grid) ÷ 2 + 1)
        gen = () -> prefix10(SubtreeHaloIterator(sys, root, 12, Vertex(),
            DGG.Fallbacks.generic_halo_engine(sys, root, 12, Vertex())))
        gen()
        dir = () -> prefix10(SubtreeHaloIterator(sys, root, 12))
        dir()
        @test @allocated(gen()) > 100 * @allocated(dir())
    end

    # -----------------------------------------------------------------------
    # A5 — the one system with no specialization, and the assertion that says so
    # -----------------------------------------------------------------------

    # Five systems route every root to a specialization and the classifying arms
    # above assert exactly that. A5 is the sixth and is supposed to route
    # NOWHERE: no `descendant_range` means no range to skip the subject subtree
    # by and no ordering to make a pruned descent canonical, so the scan is the
    # walk. This is a POSITIVE assertion about the engine, and it fails if A5
    # ever acquires a fast path by analogy from its aperture or its Hilbert-like
    # indexing rather than from a proved boundary automaton. Because "took the
    # right path" is worthless without "and answered correctly", the same loop
    # runs the geometry oracle and the contract bundle. Targets stop at level 2
    # (192 cells): the oracle is a scan whose per-cell cost is a pruned geometry
    # walk. Depth is covered for A5 by `law_halo` and "geometry agrees with
    # topology", both at full width.
    @testset "A5 stays on the linear scan" begin
        sys = A5System()
        for base in (0, 1), conn in (Vertex(), Edge())
            grid = levelgrid(sys, base)
            for c in sample_cells(grid, 2), l in base:2
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                # Depth zero is the native one-ring on every system, A5 included;
                # everything deeper is the scan.
                @test it.engine isa (l == level(c) ? DGG.Fallbacks.RingHaloEngine :
                                     DGG.Fallbacks.ScanHaloEngine)
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # -----------------------------------------------------------------------
    # The generic walk, still oracled where it is still the walk
    # -----------------------------------------------------------------------

    # The classifying arms assert that no root of any system falls back, which
    # says the specializations were reached and nothing about the walk they
    # would have fallen back TO. That walk is not dead code: it is what every one
    # of `hex_halo_engine`'s and `_seam_band_engine`'s guards returns, and what a
    # system added later inherits until it writes an engine of its own. No system
    # reaches it through the keyword constructor, so it is BUILT here and held to
    # both oracles. Without this arm it would be exercised only as the oracle's
    # own carrier — `forced_geometry_halo` runs the same `OutsideWalkEngine`, so
    # a bug in the enumeration or the cap prune would move both sides together
    # and show up nowhere. (The one other place it is pinned is the depth-4
    # aperture-7 arm of "the defining law".)
    @testset "the generic fallback still agrees with the oracle" begin
        for sys in HEX_SYSTEMS, base in (0, 1), conn in (Vertex(), Edge())
            grid = levelgrid(sys, base)
            roots = unique(vcat(sample_cells(grid, 3), irregular_cells(grid, 2)))
            l = base + 1
            for c in roots
                it = SubtreeHaloIterator(sys, c, l, conn,
                    DGG.Fallbacks.generic_halo_engine(sys, c, l, conn))
                @test it.engine isa DGG.Fallbacks.OutsideWalkEngine
                h = collect(it)
                @test h == forced_geometry_halo(sys, c, l, conn)
                @test h == law_halo(sys, c, l; connectivity = conn)
            end
        end
    end

    # -----------------------------------------------------------------------
    # Subset halos — `halo` on the three containers
    # -----------------------------------------------------------------------

    # `halo` IS AN ITERATOR ON ALL THREE CONTAINERS, which is the point of this
    # section. PR #19 returned an iterator on one fast path and a `Vector`
    # everywhere else, so `for x in halo(sub)` got laziness or a materialised
    # ring depending on how the subset happened to be built. Both branches are
    # iterators now — `SubtreeHaloIterator` on the rooted-complete-subtree path,
    # `SubsetHaloIterator` everywhere else — decided once, at construction.
    #
    # The arms: the three containers must agree, since a `CellVector` forgets
    # its root and a `CellLookup` is a `CellVector` wearing a hat; the rooted
    # complete subtree must DELEGATE, which only a type assertion can see; every
    # element must be out of the set; and A PUNCHED HOLE JOINS THE HALO, the one
    # law that distinguishes `halo` from `subtree_halo`.
    #
    # AND THE ORACLE IS THE GEOMETRY WALK, which the design's oracle row demands
    # and this section did not have: `SubsetMembership` decides both membership
    # and contact by `neighbors`, so holding it to `subtree_halo` is holding one
    # indexed walk to another. What is NOT reached that way is a subset that is
    # not a subtree — `ForcedGeometry` asks about descendants of a root and a
    # scattered set of ids has none — so every subset below is a subtree or a
    # subtree with pieces removed.
    #
    # `holed_halo_oracle` is that oracle with a hole in it: a cell outside the
    # subtree is a halo cell of what remains exactly when something it touches
    # is still held, and a removed cell is one by the same rule. Geometry plus
    # `removed` is a COMPLETE candidate list, because a halo cell of a subset of
    # the subtree touches the subtree and so is a halo cell of the whole.
    function holed_halo_oracle(sys, c, l, removed, conn, whole)
        grid = levelgrid(sys, l)
        lc = level(c)
        gone = Set(removed)
        held(x) = ancestor(sys, x, lc) == c && !(x in gone)
        out = eltype(whole)[]
        for x in vcat(whole, collect(removed))
            any(held, neighbors(grid, x, 1; connectivity = conn)) && push!(out, x)
        end
        sort!(out; by = x -> cellposition(grid, x))
        return out
    end

    @testset "$(nameof(typeof(sys))): halo on subsets" for sys in systems()
        l = min(2, max_level(sys))
        c = cellindex(levelgrid(sys, 0), 1)
        pg = PartialGrid(sys, c, l)
        cv = CellVector(pg)
        expected = subtree_halo(sys, c, l)
        @test collect(halo(pg)) == expected
        @test collect(halo(cv)) == expected
        @test collect(halo(CellLookup(cv))) == expected
        # And against the geometry oracle, on both connectivities. `cv` is a
        # `CellVector`, so this reaches `SubsetHaloIterator` on every system
        # including the five whose rooted grid delegates.
        geom = Dict(conn => forced_geometry_halo(sys, c, l, conn)
                    for conn in (Vertex(), Edge()))
        for conn in (Vertex(), Edge())
            @test collect(halo(cv; connectivity = conn)) == geom[conn]
        end

        # The rooted complete subtree delegates rather than re-deriving. Asserted
        # by type, because every arm above passes either way.
        #
        # A5 IS THE EXCEPTION AND IT IS THE SAME ONE AS EVERYWHERE ELSE:
        # `_whole_subtree_range` needs `has_sorted_subtrees` to know that holding
        # `length(descendant_range)` cells means holding the whole subtree, and
        # A5 has no `descendant_range` at all. So A5's rooted grid takes the
        # subset walk — which is the scan, the same engine `subtree_halo` would
        # have used — and the two arms above are what say the answer is identical.
        @test halo(pg) isa (DGG.has_sorted_subtrees(sys) ? SubtreeHaloIterator :
                            DGG.Fallbacks.SubsetHaloIterator)
        @test halo(cv) isa DGG.Fallbacks.SubsetHaloIterator

        # The same cells with the root forgotten must give the same answer.
        loose = PartialGrid(sys, l, collect(pg.ids))
        @test halo(loose) isa DGG.Fallbacks.SubsetHaloIterator
        @test collect(halo(loose)) == expected

        @test collect(halo(pg; connectivity = Edge())) == geom[Edge()]
        @test all(x -> cellposition(pg, x) === nothing, collect(halo(pg)))
        @test all(x -> cellposition(loose, x) === nothing, collect(halo(loose)))

        # PUNCHED HOLES, IN THREE SHAPES AND UNDER BOTH CONNECTIVITIES. One
        # interior cell is the smallest hole there is and can only ever ADD to
        # the halo, so on its own it says nothing about the two ways a hole can
        # be interesting:
        #
        #   * THE WHOLE INTERIOR is a connected patch, so a removed cell's own
        #     neighbours are removed too and being outside the subset stops
        #     being sufficient — a cell of the patch is a halo cell only if
        #     something it touches is still held.
        #   * A BORDER PATCH straddles the seam. `c` is a level-0 root, so on
        #     the three square systems its subtree is a whole face and its
        #     border is flush with all four face edges; removing the first three
        #     border cells removes the corner of that face, and the cells that
        #     lose their last contact are on OTHER FACES. It is the only shape
        #     here where cells the whole subtree's halo held DROP OUT — measured
        #     two or three of them on every system, against none for either of
        #     the other two shapes.
        #
        # Each held to `holed_halo_oracle`, so the cells to be found are
        # enumerated by geometry rather than by the walk under test.
        interior = collect(DGG.subtree_interior(sys, c, l))
        # one interior cell, the whole interior, a border patch
        holes = (interior[1:1], interior, subtree_border(sys, c, l)[1:min(3, end)])
        for removed in holes
            isempty(removed) && continue
            ids = filter(!in(Set(removed)), collect(pg.ids))
            holed = PartialGrid(sys, l, ids)
            for conn in (Vertex(), Edge())
                want = holed_halo_oracle(sys, c, l, removed, conn, geom[conn])
                hh = collect(halo(holed; connectivity = conn))
                @test hh == want
                @test hh == collect(halo(CellVector(holed); connectivity = conn))
                @test hh == collect(halo(CellLookup(CellVector(holed));
                    connectivity = conn))
                @test allunique(hh)
                @test all(x -> cellposition(holed, x) === nothing, hh)
                # ROOTED BUT NOT COMPLETE, the branch neither the whole
                # subtree nor the root-forgotten form reaches:
                # `_whole_subtree_range` refuses it on the count, so it takes
                # the subset walk. THE ROOT MUST MAKE NO DIFFERENCE — the walk
                # prunes by the subset's own spans and never asks whether one was
                # declared — and this is the arm that says so. Same ids, same
                # answer.
                rooted = PartialGrid(sys, l, ids; root = c)
                @test halo(rooted) isa DGG.Fallbacks.SubsetHaloIterator
                @test collect(halo(rooted; connectivity = conn)) == hh
            end
        end
    end

    # THE ONE ASSUMPTION THE SUBSET WALK PRUNES BY, checked rather than argued.
    #
    # `SubsetMembership` has no root, so the subset walk has no root cap and no
    # covering law to prune with. What it prunes with instead is the
    # COARSE-CONTAINMENT LAW: for every pair of cells the system calls
    # vertex-adjacent at a level, the two parents are equal or vertex-adjacent
    # themselves. Compose it down the generations and it says that a neighbour of
    # a target-level cell has its level-`lc` ancestor inside the CLOSED one-ring
    # of that cell's own level-`lc` ancestor — which is exactly `_near_subset`'s
    # licence to skip a node no neighbour of which holds a member.
    #
    # WHY IT IS HERE AND NOT LEFT TO THE VALUE ARMS. A prune that is too eager
    # does not throw and does not mis-order: it DROPS a halo cell, and it drops
    # it only where the law is violated, which is a seam or a pentagon on some
    # level the sampled subsets above may never reach. Stating the law directly
    # and running it over EVERY adjacent pair is the difference between "no
    # sample found a violation" and "there is none".
    #
    # `Vertex()` only, because that is the connectivity the probe uses whatever
    # was requested, and the `Edge()` halo is a subset of the `Vertex()` one. The
    # sweep is exhaustive to level 6, and stops at the first level past 300000
    # cells — IGeo7 at 5 and H3 at 4 — which costs a couple of seconds in total
    # and is where the aperture-7 seams and pentagons already are.
    #
    # ONE assertion per system rather than one per level, because the levels are
    # one law and not six: a failure prints the `(level, count)` pairs, which is
    # everything a per-level arm would have said.
    @testset "the coarse-containment law the subset prune rests on" begin
        for sys in systems()
            escaped = Tuple{Int,Int}[]
            for l in 1:6
                l <= max_level(sys) || continue
                grid = levelgrid(sys, l)
                ncells(grid) <= 300_000 || continue
                coarse = levelgrid(sys, l - 1)
                out = 0
                for p in 1:ncells(grid)
                    x = cellindex(grid, p)
                    a = ancestor(sys, x, l - 1)
                    ring = neighbors(coarse, a, 1; connectivity = Vertex())
                    for y in neighbors(grid, x, 1; connectivity = Vertex())
                        b = ancestor(sys, y, l - 1)
                        (b == a || any(==(b), ring)) || (out += 1)
                    end
                end
                out == 0 || push!(escaped, (l, out))
            end
            @test escaped == Tuple{Int,Int}[]
        end
    end

    # A rooted subtree with one interior cell removed: the smallest departure
    # from a subtree there is, and the one that forces the subset walk — a
    # complete subtree would delegate and measure the wrong engine.
    function holed_subtree(sys, c, l)
        r = DGG.descendant_range(sys, c, l)
        grid = levelgrid(sys, l)
        ids = [cellindex(grid, p) for p in r]
        deleteat!(ids, length(ids) ÷ 2)
        return PartialGrid(sys, l, ids; root = c)
    end

    # AND THE COMPLEXITY CLASS THAT LAW BUYS, counted rather than timed.
    #
    # The walk asks its subset two questions and nothing else — `subset_span` for
    # a block, `cellposition` for a cell — so `CountingSubset` measures its work
    # exactly, with no clock and no allocation proxy. The fixture is a rooted
    # subtree with one interior cell punched out, which is the irregular chunk
    # `halo` exists for: five levels apart its MEMBER count grows by a thousand
    # and its HALO by thirty, so a walk sized by the input and one sized by the
    # answer cannot both pass.
    #
    # TWO STATEMENTS, because one of them alone is not enough. The RATIO arm says
    # the growth followed the answer, and it is what refuses a walk with no node
    # prune at all — but a walk that visits its members and prunes nothing else
    # can still slip under a ratio threshold, since both endpoints grow together.
    # The FRACTION arm is the one that refuses that: a walk asking fewer questions
    # than the subset has members has demonstrably not visited them. Measured
    # 0.237, 0.203, 0.205 (the square systems), 0.389 (IGeo7) and 0.373 (H3),
    # against ratio-to-halo figures of 0.94, 1.06, 0.93, 0.86 and 0.60 — so both
    # thresholds sit a factor of 1.5 clear, and a walk that descended into the
    # members misses the second by four.
    #
    # A5 IS EXCLUDED FOR ITS USUAL REASON: with no `descendant_range` there are
    # no spans to prune by, so `subset_halo_engine` hands back `ScanHaloEngine`
    # and the cost is `O(ncells)` by construction, not by regression.
    @testset "the subset walk's cost follows the halo, not the subset" begin
        for sys in systems()
            DGG.has_sorted_subtrees(sys) || continue
            c = cellindex(levelgrid(sys, 0), 1)
            # Five and three levels of descent, which are the same subtree size
            # either side of the aperture: about 260000 members on the aperture-4
            # systems and about 100000 on the aperture-7 ones. The whole arm
            # costs a tenth of a second.
            depths = any(h -> sys isa typeof(h), HEX_SYSTEMS) ? (3, 6) : (4, 9)
            all(l -> l <= max_level(sys), depths) || continue
            calls = Int[]; members = Int[]; halos = Int[]
            for l in depths
                cv = CellVector(holed_subtree(sys, c, l))
                counted = CountingSubset(cv)
                h = counting_halo(sys, counted, levelgrid(sys, l), l, Vertex())
                # The wrapper changed the measurement, not the walk.
                @test h == collect(halo(cv))
                push!(calls, counted.calls)
                push!(members, length(cv))
                push!(halos, length(h))
            end
            # The two hypotheses are far enough apart to tell apart ...
            @test members[2] / members[1] >= 4 * (halos[2] / halos[1])
            # ... the growth went with the answer ...
            @test calls[2] <= 1.6 * (halos[2] / halos[1]) * calls[1]
            # ... and the walk never asked about most of the members at all.
            @test calls[2] <= 0.6 * members[2]
        end
    end

    # A MIXED-LEVEL SET HAS NO `halo`, AND THAT ABSENCE IS A DECISION.
    # `MultiOrderCellSet`'s members sit at different levels, so "the cells just
    # outside it" has no one level to be answered at and the verb would have to
    # invent one. The `MethodError` is pinned here the way the missing `length`s
    # are. Mixed-level adjacency is `member_neighbors`, in `multiorder_*.jl`.
    @testset "a mixed-level set has no halo" begin
        for sys in systems()
            c = cellindex(levelgrid(sys, 1), 1)
            set = DGG.MultiOrderCellSet(sys, [c], [1], trues(1), level(c))
            @test_throws MethodError halo(set)
        end
    end

    # -----------------------------------------------------------------------
    # Resumability, on every engine this file can reach
    # -----------------------------------------------------------------------

    # THE PREFIX-EQUALITY LAW, mirroring `subtree_iterators.jl`'s: four cells
    # taken off the front must be the first four of the collected form, and a
    # second `collect` must reproduce the whole walk. That is the observable half
    # of the house rule that an engine is immutable and ALL walk state travels in
    # the value `iterate` threads — an engine that cached a cursor in a mutable
    # field passes every oracle here and fails on the second pass.
    #
    # Run against EVERY engine type rather than one per system, because
    # resumability is a property of the STATE and the seven walks thread seven
    # different ones. `SquareBandEngine` is two walks wearing one name, so
    # `engine_tag` reads the emit rule too; it is at file scope because three
    # testsets partition the engines by it. The tag sets asserted at the end
    # stop a list from quietly losing an engine: one that stopped being
    # reachable drops out of the set and fails the equality.
    engine_tag(e) = e isa DGG.Fallbacks.SquareBandEngine ?
        (e.check isa DGG.Fallbacks.NoCheck ? :SquareBandNoCheck :
         :SquareBandNativeCheck) : nameof(typeof(e))

    ALL_ENGINE_TAGS = Set((:RingHaloEngine, :OutsideWalkEngine, :ScanHaloEngine,
        :SquareBandNoCheck, :SquareBandNativeCheck, :HexChildHaloEngine,
        :HexArcHaloEngine))

    # An in-face square root, the only way to reach the counted emit rule: a
    # level-0 or level-1 block is flush with its face edge on all three systems.
    #
    # A SPREAD PICK, NOT `first`. The in-face class is face-major, so `first` is
    # face 0's lowest block on every system and at base 3 so are the first six —
    # every call site of this helper used to land on face 0. Taking the last of
    # a five-way spread lands on a high-numbered face instead, which matters
    # because the S2 orientation seed `BAND_BASES` exists to defend is per-face
    # and face 0 is the one face where a wrong seed is a no-op.
    inface_root(sys, base::Int, l::Int) =
        last(spread(first(classify_roots(sys, base, l, Vertex())), 5))

    @testset "the walk is resumable, not restarted" begin
        seen = Set{Symbol}()
        wrappers = Set{Symbol}()
        function check_prefix(it)
            push!(seen, engine_tag(it.engine))
            push!(wrappers, nameof(typeof(it)))
            full = collect(it)
            @test length(full) >= 4
            prefix = eltype(it)[]
            for x in it
                push!(prefix, x)
                length(prefix) >= 4 && break
            end
            @test prefix == full[1:4]
            @test collect(it) == full        # a second pass gives the same walk
        end
        for sys in systems()
            c = cellindex(levelgrid(sys, 1), 1)
            l = min(level(c) + 2, max_level(sys))
            check_prefix(SubtreeHaloIterator(sys, c, l))              # shipped
            check_prefix(SubtreeHaloIterator(sys, c, level(c)))       # one-ring
            # The generic walk is no longer reachable through the keyword
            # constructor on any system, so it is built explicitly — and on A5
            # that same call is the scan, which covers both fallbacks without
            # naming either system.
            check_prefix(SubtreeHaloIterator(sys, c, l, Vertex(),
                DGG.Fallbacks.generic_halo_engine(sys, c, l, Vertex())))
        end
        # Depth one on the aperture-7 systems is the automaton-free child walk,
        # which the `l = level(c) + 2` cases above never reach.
        for sys in HEX_SYSTEMS
            check_prefix(SubtreeHaloIterator(sys, cellindex(levelgrid(sys, 1), 1), 2))
        end
        # Both square emit rules, on the blocks the constructor actually claims.
        for sys in SQUARE_SYSTEMS
            inface, seam, _ = classify_roots(sys, 2, 4, Vertex())
            check_prefix(SubtreeHaloIterator(sys, last(spread(inface, 5)), 4))
            check_prefix(SubtreeHaloIterator(sys, last(spread(seam, 5)), 4))
        end
        # And the subset wrapper on both containers: `SubsetHaloIterator`
        # forwards the whole protocol itself, so resumability is a property of
        # the wrapper as much as of the engine inside it.
        for sys in systems()
            l = min(2, max_level(sys))
            c = cellindex(levelgrid(sys, 0), 1)
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c, l).ids))
            check_prefix(halo(loose))
            check_prefix(halo(CellVector(loose)))
        end
        @test seen == ALL_ENGINE_TAGS
        # And BOTH WRAPPERS, a separate claim: `SubsetHaloIterator` forwards
        # `iterate` itself, so a cursor cached in the wrapper would be invisible
        # to an engine-only tag set.
        @test wrappers == Set((:SubtreeHaloIterator, :SubsetHaloIterator))
    end

    # -----------------------------------------------------------------------
    # Type stability in `eltype`, on every engine and every system
    # -----------------------------------------------------------------------

    # WHY THE QUESTION IS ASKED OF THE TYPE AND NOT OF THE INSTANCE.
    # `eltype(typeof(it))` is what `collect`, `Iterators.partition` and every
    # `Vector{eltype(it)}` preallocation dispatch on, and an engine that only
    # knew its element type at run time would answer `Any` there while still
    # yielding perfectly good cells — a silent `Vector{Any}` in every caller's
    # hands, nothing red anywhere. The three assertions per engine are the
    # type-domain answer, its concreteness, and the run-time consequence.
    @testset "eltype is the system's cell index type, on every engine" begin
        seen = Set{Symbol}()
        wrappers = Set{Symbol}()
        function check_eltype(sys, it)
            push!(seen, engine_tag(it.engine))
            push!(wrappers, nameof(typeof(it)))
            C = DGG.cellindextype(sys)
            @test eltype(typeof(it)) === C
            @test isconcretetype(eltype(typeof(it)))
            @test collect(it) isa Vector{C}
        end
        for sys in systems()
            mx = max_level(sys)
            c0 = cellindex(levelgrid(sys, 0), 1)
            check_eltype(sys, SubtreeHaloIterator(sys, c0, 0))       # the one-ring
            for l in 1:min(2, mx)
                check_eltype(sys, SubtreeHaloIterator(sys, c0, l))   # what it ships
            end
            l = min(2, mx)
            check_eltype(sys, generic_iterator(sys, c0, l))          # walk, or scan
            pg = PartialGrid(sys, c0, l)
            loose = PartialGrid(sys, l, collect(pg.ids))
            check_eltype(sys, halo(pg))
            check_eltype(sys, halo(loose))
            check_eltype(sys, halo(CellVector(loose)))
            check_eltype(sys, halo(CellLookup(CellVector(loose))))
        end
        # The counted square emit rule, which no level-0 or level-1 root can
        # reach: those blocks are flush with their face edge on all three
        # systems.
        for sys in SQUARE_SYSTEMS
            check_eltype(sys, SubtreeHaloIterator(sys, inface_root(sys, 2, 4), 4))
        end
        @test seen == ALL_ENGINE_TAGS
        @test wrappers == Set((:SubtreeHaloIterator, :SubsetHaloIterator))
    end

    # -----------------------------------------------------------------------
    # The count contract, for ALL SEVEN ENGINES AT ONCE
    # -----------------------------------------------------------------------

    # Two engines count in closed form and five refuse; this pins that as a
    # PARTITION rather than engine by engine.
    #
    # THE CLAIMED COUNTS ARE NOT RE-CHECKED AGAINST A COLLECT, because they
    # cannot be: `collect` routes through `collect_subtree`, which `error`s when
    # a `HasLength` walk emits a different number of cells than it claims, so
    # `length(it) == length(collect(it))` is that call's post-condition. It can
    # throw. It cannot fail. Written as an assertion it would read as coverage
    # and be none. "the band count is closed form" holds the surviving count to
    # a FORMULA; "collect is the guarded path" pins the guard, on a lying engine.
    #
    # WHAT THE PARTITION BUYS is the other half of the contract: an engine that
    # cannot count declares `SizeUnknown()` and defines NO method, so the
    # `MethodError` is the honest answer and a `length` that quietly walked the
    # halo would satisfy every caller while violating the design. The two set
    # assertions at the end say an engine that moved sides leaves the set it was
    # in. `RingHaloEngine` counts because the ring is already in hand;
    # `SquareBandEngine`/`NoCheck` because `4·side + 4` is the band.
    # `SquareBandEngine`/`NativeCheck` does not, because no perimeter formula
    # survives a seam; neither hexagonal engine does, because its census is
    # validated by enumeration rather than derived; and neither generic walk
    # does, because counting would be walking.
    @testset "length is truthful where it exists and absent where it does not" begin
        counted = Set{Symbol}()
        refusing = Set{Symbol}()
        function check_count(it)
            tag = engine_tag(it.engine)
            h = collect(it)
            @test !isempty(h)
            if Base.IteratorSize(typeof(it)) isa Base.HasLength
                push!(counted, tag)
            else
                @test Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
                push!(refusing, tag)
                @test_throws MethodError length(it)
            end
        end
        for sys in systems()
            mx = max_level(sys)
            c0 = cellindex(levelgrid(sys, 0), 1)
            check_count(SubtreeHaloIterator(sys, c0, 0))
            for l in 1:min(2, mx)
                check_count(SubtreeHaloIterator(sys, c0, l))
            end
            l = min(2, mx)
            check_count(generic_iterator(sys, c0, l))
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c0, l).ids))
            check_count(halo(loose))
            check_count(halo(CellVector(loose)))
        end
        for sys in SQUARE_SYSTEMS
            for d in 1:3
                l = 2 + d
                l <= max_level(sys) || continue
                check_count(SubtreeHaloIterator(sys, inface_root(sys, 2, l), l))
            end
        end
        @test counted == Set((:RingHaloEngine, :SquareBandNoCheck))
        @test refusing == Set((:OutsideWalkEngine, :ScanHaloEngine,
            :SquareBandNativeCheck, :HexChildHaloEngine, :HexArcHaloEngine))
    end

    # -----------------------------------------------------------------------
    # "Consumable incrementally or in caller-selected batches"
    # -----------------------------------------------------------------------

    # The prefix law above says a partial walk is a prefix. This says the two
    # standard ways a caller cuts a lazy stream into pieces — `Iterators.take`
    # and `Iterators.partition` — put the pieces back together into exactly the
    # walk, at every chunk size including ones that do not divide it. Not free
    # given `SizeUnknown()`: `partition` takes a different route for a sized
    # iterator, and an engine whose `iterate` mutated shared state would
    # reassemble into something shorter than the collect.
    @testset "consumable incrementally and in caller-chosen batches" begin
        function check_batches(it)
            full = collect(it)
            @test length(full) >= 6
            @test collect(Iterators.take(it, 3)) == full[1:3]
            @test collect(Iterators.take(it, length(full) + 5)) == full
            @test isempty(collect(Iterators.take(it, 0)))
            for k in (1, 2, 5)
                parts = collect.(Iterators.partition(it, k))
                @test reduce(vcat, parts) == full
                @test all(p -> 1 <= length(p) <= k, parts)
            end
        end
        for sys in systems()
            l = min(2, max_level(sys))
            c0 = cellindex(levelgrid(sys, 0), 1)
            check_batches(SubtreeHaloIterator(sys, c0, l))
            check_batches(generic_iterator(sys, c0, l))
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c0, l).ids))
            check_batches(halo(loose))
        end
        for sys in SQUARE_SYSTEMS
            check_batches(SubtreeHaloIterator(sys, inface_root(sys, 2, 4), 4))
        end
    end

    # -----------------------------------------------------------------------
    # The one awkward cell that DEGREE does not find
    # -----------------------------------------------------------------------

    # `irregular_cells` catches pentagons, cube corners and icosahedral vertices
    # because their one-ring is the wrong size. A HEALPix POLE is not reliably
    # like that: it has degree eight, and from base 2 down that is the MODAL
    # degree (168 of 192 cells at base 2, 744 of 768 at base 3), so a
    # degree-based sweep passes it by and reaches it only by luck of sampling —
    # while it is a genuinely distinct configuration, the meeting point of the
    # four polar faces, where `nested_neighbors` runs its `-1` entries and the
    # seam walk's rectangles come from three faces rather than two. (BASE 1 IS
    # THE EXCEPTION, and this arm sweeps it: the degrees split 24 to 24,
    # `irregular_cells` resolves the tie to seven, and the north pole is the
    # first cell it returns. The argument is about bases 2 and 3.)
    #
    # So it is pinned BY LOCATION: `cellat` at ±89.999° names it at any level
    # without this file knowing a single HEALPix identifier. Both poles, three
    # bases, two depths, both connectivities, against the geometry oracle and
    # the whole contract bundle; and one arm against the `O(ncells)` law.
    @testset "the two HEALPix poles, pinned by location rather than by degree" begin
        sys = HEALPixSystem()
        for base in (1, 2, 3), lat in (89.999, -89.999)
            grid = levelgrid(sys, base)
            p = DGG.cellat(grid, 0.0, lat)
            @test p !== nothing
            for d in 1:2, conn in (Vertex(), Edge())
                l = base + d
                l <= max_level(sys) || continue
                @test collect(SubtreeHaloIterator(sys, p, l; connectivity = conn)) ==
                      forced_geometry_halo(sys, p, l, conn)
                check_halo_case(sys, p, l, conn)
            end
        end
        grid = levelgrid(sys, 1)
        for lat in (89.999, -89.999), conn in (Vertex(), Edge())
            p = DGG.cellat(grid, 0.0, lat)
            @test collect(SubtreeHaloIterator(sys, p, 3; connectivity = conn)) ==
                  law_halo(sys, p, 3; connectivity = conn)
        end
    end

    # -----------------------------------------------------------------------
    # THE ALLOCATION LAWS
    # -----------------------------------------------------------------------

    # WHAT IS ASSERTED HERE, AND WHAT DELIBERATELY IS NOT. The design's last
    # verification row is "construction and short-prefix allocation independent
    # of halo cardinality", measured "after warm-up, not by brittle wall-clock
    # thresholds". Four arms, in increasing order of what they claim:
    #
    #   1. CONSTRUCTION IS FLAT IN THE TARGET LEVEL, on every engine.
    #   2. THE SPECIALIZED PREFIX IS FLAT IN THE TARGET LEVEL — the strong form,
    #      asserted only where it is true: the square band walk and the two
    #      hexagonal walks reach their first cell in `O(depth)` node visits with
    #      no geometry at all.
    #   3. THE FRACTION LAW, on every engine including the generic walk: a short
    #      prefix costs a small and NON-GROWING share of a full collect while
    #      the halo grows by an order of magnitude. Unlike (2) it holds
    #      everywhere, and it is what separates a lazy iterator from an eager.
    #   4. The subset walk, weaker for the reason given at its own arm.
    #
    # And a fifth that is none of them: `EagerHaloEngine`, a correct iterator
    # every one of the four refuses. Without it the four are thresholds nobody
    # has watched fail.
    #
    # (2) DOES NOT HOLD FOR THE GENERIC WALK AND IS NOT ASSERTED OF IT —
    # `_admit` calls `node_extent` on every node it prunes and that allocates,
    # so a four-cell IGeo7 prefix is 3136 B at `l = 3` and 213216 B at `l = 6`:
    # proportional to the WORK, not the OUTPUT.
    #
    # NO ARM BELOW PINS A RAW BYTE COUNT. Every assertion is a difference
    # between two measurements of the same shape, or a ratio; the byte figures
    # in these comments are what this machine measured, and say how much
    # headroom a threshold has.

    # The five systems whose halo engine reaches its first cell in `O(depth)`
    # node visits. A5 IS EXCLUDED BY NAME rather than by a `filter` nobody can
    # read: `ScanHaloEngine` reaches its first cell by scanning positions upward
    # from 1, so the prefix cost is a fact about where that cell sits in the
    # target level's ordering and not about the target level at all — measured
    # 53776, 20864, 19504 and 56576 B at levels 1 to 4 from one root, neither
    # flat nor monotone. Arm 3 is where A5 is held to the law it does obey.
    DEPTH_FLAT_SYSTEMS = (HEALPixSystem(), S2System(), ISEA4RSystem(),
        H3System(), IGeo7System())

    @testset "construction does not allocate in proportion to the halo" begin
        for sys in systems()
            c = cellindex(levelgrid(sys, 0), 1)
            # A5's targets stop at 3: its `subtree_halo` at level 4 is a
            # 3840-cell scan whose per-cell cost is a `Set`-allocating
            # `neighbors`, and the law here needs only two comparable points.
            depths = filter(l -> l <= max_level(sys),
                sys isa DGG.A5System ? (1, 2, 3) : (3, 5, 7))
            ship = [construct_bytes(sys, c, l) for l in depths]
            gen = [generic_construct_bytes(sys, c, l) for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            @test last(sizes) >= 2 * first(sizes)
            # Measured EXACTLY flat on every system and both engines: 880 B
            # (HEALPix, ISEA4R band), 1392 (S2 band), 256 (both hex walks), 64
            # (A5 scan); 576/1280/1360/416/1328 for the generic walk. The 64 B
            # of slack is for a future engine that rounds an allocation
            # differently, not for a trend.
            @test maximum(ship) - minimum(ship) <= 64
            @test maximum(gen) - minimum(gen) <= 64
        end
    end

    @testset "the specialized prefix costs the same at every depth" begin
        for sys in DEPTH_FLAT_SYSTEMS
            c = cellindex(levelgrid(sys, 0), 1)
            depths = filter(l -> l <= max_level(sys), (3, 5, 7))
            allocs = [lazy_bytes(sys, c, l, 4) for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            # The halo grows 15x (HEALPix, ISEA4R: 34 -> 514), 16x (S2:
            # 32 -> 512) and 78x (H3 84 -> 6564, IGeo7 70 -> 5470) ...
            @test last(sizes) >= 8 * first(sizes)
            # ... and the four-cell prefix does not move at all: 880 B on
            # HEALPix and ISEA4R, 1904 on S2 (whose `neighbors` allocates a
            # `Vector`), 256 on both hexagonal systems, at every one of the
            # three depths.
            @test maximum(allocs) - minimum(allocs) <= 64
            @test all(>(0), allocs)
        end
        # And on an IN-FACE block, which is the square systems' OTHER emit
        # rule — a counted band with no native check between yields, where a
        # regression would look nothing like one in the seam walk. Flatness of
        # the in-face guard in the target level (`1 <= ix <= 2^b - 2`,
        # independent of the depth) is why one classification serves all three
        # targets. Measured 800 B at every depth from 1 to 6, against a halo
        # growing 12 -> 260.
        for sys in SQUARE_SYSTEMS
            c = inface_root(sys, 3, 4)
            depths = filter(l -> l <= max_level(sys), (4, 6, 9))
            allocs = [lazy_bytes(sys, c, l, 4) for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            @test last(sizes) >= 8 * first(sizes)
            @test maximum(allocs) - minimum(allocs) <= 64
        end
    end

    @testset "a short prefix is a small, non-growing fraction of the collect" begin
        for sys in systems()
            c = cellindex(levelgrid(sys, 0), 1)
            # A5 goes 1 -> 4 rather than 3 -> 7: its halo from a level-0 root
            # grows 15 -> 80 over those levels, and level 5 would be a
            # 15360-cell scan for no extra claim.
            depths = sys isa DGG.A5System ? (1, 4) : (3, 7)
            all(l -> l <= max_level(sys), depths) || continue
            fracs = [lazy_bytes(sys, c, l, 4) / eager_bytes(sys, c, l)
                     for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            # The halo grew: 15x to 78x on the five, 5.3x on A5.
            @test last(sizes) >= 5 * first(sizes)
            # The share the prefix costs did not.
            @test last(fracs) <= first(fracs)
            # And it is small. Measured at the deep end: 0.060 (HEALPix,
            # ISEA4R), 0.024 (S2), 0.0023 (IGeo7), 0.0011 (H3), 0.0032 (A5) —
            # against 0.50, 0.30, 0.13, 0.13 and 0.29 at the shallow end, which
            # is why the threshold is on the deep end only.
            @test last(fracs) < 0.15
        end
    end

    @testset "the generic walk obeys the same fraction law" begin
        for sys in systems()
            # A5's `halo_engine` IS `generic_halo_engine` — with no
            # `descendant_range` it returns the scan — so the arm above already
            # measured exactly this walk on it, and repeating it here would buy
            # a second copy of the same numbers and a 17 MB collect.
            sys isa DGG.A5System && continue
            c = cellindex(levelgrid(sys, 0), 1)
            depths = (3, 6)
            all(l -> l <= max_level(sys), depths) || continue
            fracs = Float64[]
            sizes = Int[]
            for l in depths
                h = generic_collect(sys, c, l)             # warm up, and count
                eb = @allocated generic_collect(sys, c, l)
                generic_take(sys, c, l, 4)                 # warm up
                lb = @allocated generic_take(sys, c, l, 4)
                push!(fracs, lb / eb)
                push!(sizes, length(h))
            end
            # 34 -> 258 on HEALPix and ISEA4R, 32 -> 256 on S2, 70 -> 1825 on
            # IGeo7, 84 -> 2190 on H3.
            @test last(sizes) >= 5 * first(sizes)
            # 0.067 -> 0.036 (IGeo7), 0.411 -> 0.024 (H3), 0.132 -> 0.069
            # (HEALPix), 0.111 -> 0.008 (S2), 0.060 -> 0.003 (ISEA4R). The
            # prefix cost grows with the descent — see the section comment —
            # but never as fast as the halo does, which is the whole claim.
            @test last(fracs) <= first(fracs)
            @test last(fracs) < 0.15
        end
    end

    # THE FIXTURE THE THREE LAWS ABOVE — AND THE ONE BELOW — EXIST TO REFUSE.
    # `EagerHaloEngine` is the engine they are written against: correct on every
    # system, same public surface, same cells in the same order, and eager. It
    # passes every other assertion in this file, which its first assertion here
    # says, so nothing else would notice it.
    #
    # Measured on IGeo7 from a level-0 root at levels 3, 5 and 7, where the halo
    # runs 70 -> 610 -> 5470 cells: its construction spans 111744 B where a
    # shipped engine's spans zero, its four-cell prefix spans 223488 B, and that
    # prefix is 0.82 of a full collect where the shipped walk is 0.0023. The
    # assertions are the LAWS NEGATED, not those numbers — each is the arm above
    # with its comparison turned round, so a law that was relaxed would stop
    # failing here and this arm would go red.
    @testset "an eager engine with the same surface fails every allocation law" begin
        for sys in systems()
            c = cellindex(levelgrid(sys, 0), 1)
            depths = filter(l -> l <= max_level(sys),
                sys isa DGG.A5System ? (1, 2, 3) : (3, 5, 7))
            # It is a correct iterator. That is the whole point of it.
            @test collect(fixture_iterator(sys, c, last(depths))) ==
                  subtree_halo(sys, c, last(depths))
            ctor = [fixture_construct_bytes(sys, c, l) for l in depths]
            pref = [fixture_prefix_bytes(sys, c, l, 4) for l in depths]
            @test maximum(ctor) - minimum(ctor) > 64      # arm 1 refuses it
            @test maximum(pref) - minimum(pref) > 64      # arm 2 refuses it
            @test fixture_prefix_bytes(sys, c, last(depths), 4) >=
                  0.15 * fixture_collect_bytes(sys, c, last(depths))   # arm 3
            # And arm 4's shape, which is the one that hoists construction out
            # of the measurement: an engine that pays the halo on every walk is
            # still refused there, which is why the fixture materialises per
            # walk and not only in its constructor.
            it = fixture_iterator(sys, c, first(depths))
            take_n(it, 4)
            lb = @allocated take_n(it, 4)
            DGG.collect_subtree(it)
            eb = @allocated DGG.collect_subtree(it)
            @test lb > 0.7 * eb
        end
    end

    # THE SUBSET WALK NOW OBEYS BOTH LAWS, AND IT DID NOT ALWAYS. `halo` used to
    # summarise the subset into a bounding cap at construction — one
    # `cell_boundary` per member, up to a fixed batch, and the whole sphere past
    # it — so its construction was INPUT-sized and the only law assertable here
    # was the walk's. Nothing is read from the subset before the first `iterate`
    # now, so construction is FLAT across two inputs whose member counts differ
    # by more than an order of magnitude, and any per-member construction cost
    # put back would fail here rather than in a benchmark nobody runs.
    #
    # The second arm is the one that was always assertable: a four-cell prefix
    # against a full collect, with the iterator hoisted out of the measurement so
    # that an engine which materialised in its constructor could not read zero.
    # `EagerHaloEngine` pays the halo on every walk and reads 0.78 to 1.00 on the
    # same six systems, which is why this threshold is 0.7 rather than a number
    # nobody has watched fail.
    @testset "the subset walk is lazy, and its construction is O(1)" begin
        for sys in systems()
            mx = max_level(sys)
            ctor = (grid = Int[], vector = Int[])
            sizes = Int[]
            for l in unique((1, min(4, mx)))
                c = cellindex(levelgrid(sys, 0), 1)
                loose = PartialGrid(sys, l, collect(PartialGrid(sys, c, l).ids))
                push!(sizes, ncells(loose))
                for (built, sub) in ((ctor.grid, loose),
                        (ctor.vector, CellVector(loose)))
                    push!(built, subset_construct_bytes(sub))
                    it = halo(sub)
                    take_n(it, 4)
                    lb = @allocated take_n(it, 4)
                    DGG.collect_subtree(it)
                    eb = @allocated DGG.collect_subtree(it)
                    @test eb > 0
                    @test lb <= 0.7 * eb
                end
            end
            # The input grew by more than an order of magnitude ...
            @test last(sizes) >= 10 * first(sizes)
            # ... and building the iterator did not notice, on either container.
            @test maximum(ctor.grid) - minimum(ctor.grid) <= 64
            @test maximum(ctor.vector) - minimum(ctor.vector) <= 64
        end
    end

    # -----------------------------------------------------------------------
    # The wrapper, and the guard on a lying count
    #
    # Last rather than first, because the wrapper test needs `classify_roots`.
    # -----------------------------------------------------------------------

    # The authalic transform moves where a cell is DRAWN, not which cells are
    # adjacent, so the halo through the wrapper must be the halo without it — the
    # same ids, in the same order. `halo_engine(::AuthalicSystem, ...)` is one
    # forwarding line, and this says the line is there.
    @testset "AuthalicSystem forwards the halo walk" begin
        seen = Set{Symbol}()
        for sys in systems()
            wrapped = DGG.AuthalicSystem(sys)
            grid0 = levelgrid(sys, 0)
            c = cellindex(grid0, 1)
            for l in level(c):min(level(c) + 2, max_level(sys))
                it = SubtreeHaloIterator(wrapped, c, l)
                push!(seen, engine_tag(it.engine))
                @test collect(it) == collect(SubtreeHaloIterator(sys, c, l))
            end
        end
        # A level-0 root is flush on all four sides of its face and is nobody's
        # ordinary cell, so the loop above reaches the one-ring, the SEAM band,
        # both hexagonal walks and A5's scan — but never the counted in-face
        # band. That one is picked up explicitly below, and the tag set is what
        # says which of the seven this testset actually forwarded.
        @test seen == Set((:RingHaloEngine, :SquareBandNativeCheck,
            :HexChildHaloEngine, :HexArcHaloEngine, :ScanHaloEngine))
        # And on a root the SPECIALIZATION claims. Forwarding that only ever ran
        # on a level-0 root would be forwarding that only ever reached the
        # generic walk — the wrapper would be free to lose the fast path.
        for sys in SQUARE_SYSTEMS
            wrapped = DGG.AuthalicSystem(sys)
            for base in BAND_BASES
                l = base + 2
                l <= max_level(sys) || continue
                inface, _, _ = classify_roots(sys, base, l, Vertex())
                isempty(inface) && continue
                c = last(spread(inface, 5))
                for conn in (Vertex(), Edge())
                    it = SubtreeHaloIterator(wrapped, c, l; connectivity = conn)
                    @test it.engine isa DGG.Fallbacks.SquareBandEngine
                    @test collect(it) ==
                          collect(SubtreeHaloIterator(sys, c, l; connectivity = conn))
                end
            end
        end
        # And on the aperture-7 specialization, for the same reason.
        for sys in HEX_SYSTEMS
            wrapped = DGG.AuthalicSystem(sys)
            c = cellindex(levelgrid(sys, 1), 1)
            for d in 1:2, conn in (Vertex(), Edge())
                it = SubtreeHaloIterator(wrapped, c, 1 + d; connectivity = conn)
                @test it.engine isa (d == 1 ? DGG.Fallbacks.HexChildHaloEngine :
                                     DGG.Fallbacks.HexArcHaloEngine)
                @test collect(it) ==
                      collect(SubtreeHaloIterator(sys, c, 1 + d; connectivity = conn))
            end
        end
        # AND THE SUBSET VERB, which reaches the wrapper by a different route:
        # `halo(pg)` reads `pg.system`, carrying it into `subset_halo_engine`
        # and into `_whole_subtree_range`'s `has_sorted_subtrees` question.
        # Every container, both branches.
        for sys in systems()
            wrapped = DGG.AuthalicSystem(sys)
            l = min(2, max_level(sys))
            c = cellindex(levelgrid(sys, 0), 1)
            pg = PartialGrid(wrapped, c, l)
            loose = PartialGrid(wrapped, l, collect(pg.ids))
            expected = subtree_halo(sys, c, l)
            @test halo(pg) isa (DGG.has_sorted_subtrees(sys) ? SubtreeHaloIterator :
                                DGG.Fallbacks.SubsetHaloIterator)
            @test collect(halo(pg)) == expected
            @test collect(halo(loose)) == expected
            @test collect(halo(CellVector(pg))) == expected
            @test collect(halo(CellLookup(CellVector(pg)))) == expected
        end
    end

    @testset "collect is the guarded path" begin
        sys = HEALPixSystem()
        c = cellindex(levelgrid(sys, 1), 1)
        lying = SubtreeHaloIterator(sys, c, 1, Vertex(), MiscountingEngine())
        @test_throws ErrorException collect(lying)
    end

end  # @testset "subtree halos"

end # module
