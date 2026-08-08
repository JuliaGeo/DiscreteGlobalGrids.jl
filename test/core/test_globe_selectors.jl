# test/core/test_globe_selectors.jl — the point selectors over a globe-complete
# lookup (src/core/lookups.jl; docs/design/full_globe_lookups.md §2.1, §4.1).
#
# `At(id)` and `Contains(point)` both reduce to one question — which position of
# the lookup holds this cell id — and on a globe that question is `cell_to_
# ordinal`, one call, rather than a binary search whose every probe of the id
# vector is an `ordinal_to_cell`. Two claims, and the suite is arranged so that
# neither can pass by accident:
#
#   * the answers do not change. At levels small enough to enumerate, `At` and
#     `Contains` on a globe lookup return, cell for cell, exactly what the
#     explicit all-ids lookup returns — the point-selector half of design §4.2
#     item 1. A fast path is worth having only if it is indistinguishable from
#     the path it replaces, and this is what makes that a claim rather than a
#     hope.
#   * an id naming no cell is rejected. The globe path never reads the ids it
#     "holds", so nothing in it would notice garbage but the round trip back
#     through `ordinal_to_cell` (§4.1); without that check a mistyped or
#     wrong-resolution id would resolve to whatever cell its bits decoded near,
#     silently and in range. What it raises instead is what the partial lookup
#     raises for the same id, message for message: the two lookups mean the same
#     thing by "not present", so they say it in the same words.
#
# Both claims are made at levels no explicit lookup could reach as well, since
# that is the case the fast path exists for and the only one where nothing else
# could answer at all.
#
# One module rather than four per-system additions, because the implementation
# is one generic method pair over the id vector. `num_cells` is deliberately
# unexported from the package (the system submodules own that name), so it is
# qualified here and never enters the shared test namespace.
module GlobeSelectorTests

using Test

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
import DimensionalData as DD
using DimensionalData.Lookups: At, Contains, selectindices

using DiscreteGlobalGrids.A5.A5Lookups: A5Lookup
using DiscreteGlobalGrids.H3.H3Lookups: H3Lookup
using DiscreteGlobalGrids.HEALPix.HealpixLookups: HealpixLookup
using DiscreteGlobalGrids.IGeo7.IGeo7Lookups: IGeo7Lookup

# One entry per kernel-wired system: its lookup type, the system, the levels
# whose every cell is compared one by one against the explicit lookup, and a
# level whose explicit id vector could not exist (the levels of
# `test/core/test_globe_ids.jl`). The small levels are the design's — H3 0-2,
# HEALPix 0-3 — with the aperture-7 and pentagon systems given the same span.
const SYSTEMS = (
    (H3Lookup, H3DGGS(), 0:2, 15),
    (A5Lookup, A5DGGS(), 0:2, 25),
    (IGeo7Lookup, IGEO7DGGS(), 0:2, 19),
    (HealpixLookup, HEALPixDGGS(), 0:3, 29),
)

# The two lookup constructors differ in one keyword name and nothing else.
explicit_lookup(::Type{HealpixLookup}, ids, level) = HealpixLookup(ids; level)
explicit_lookup(Lookup, ids, level) = Lookup(ids; resolution=level)

# `Contains` takes lon/lat degrees, and the kernel's `cell_center` is the one
# center every wired system has — on the unit sphere, hence the conversion.
function center_lonlat(system, level, id)
    point = DGG.cell_center(system, level, id)
    return (atand(point[2], point[1]), asind(clamp(point[3], -1, 1)))
end

# `Contains` answers with a one-element vector on the structural-id lookups and
# with the position itself on HEALPix's — a difference in those selectors that
# predates this path and that the equivalence testset below pins as it stands.
# Normalized only where the stronger claim is being made, that a cell's own
# center selects that cell.
contained_position(indices::Integer) = indices
contained_position(indices::AbstractVector) = only(indices)

# Ids that name no cell. `total` is `num_cells`, so it is the first ordinal id
# past the end of its level; the structural encodings need no such care, since
# 0 is not a cell in any of them and neither is an all-ones word.
garbage_ids(::Type{Int64}, total) = Int64[-1, total, total + 10^9]
garbage_ids(::Type{UInt64}, _) = UInt64[0x0, 0xdeadbeef, typemax(UInt64)]

raised(f) = try
    f()
    nothing
catch err
    err
end

@testset "globe point selectors answer what the explicit lookup answers" begin
    for (Lookup, system, levels, _) in SYSTEMS
        for level in levels
            globe = Lookup(DGGSGlobeIds(system, level))
            ids = collect(DGGSGlobeIds(system, level))
            stored = explicit_lookup(Lookup, ids, level)
            @test length(globe) == length(stored) == Int(DGG.num_cells(system, level))

            # Every cell of the level, in both selectors, against the lookup
            # that really does hold the ids. Nothing is sampled: the fast path
            # computes a position from an id rather than finding it, so an
            # off-by-one or a mis-ordered ordinal would be invisible anywhere
            # but at the cell it happens to strike.
            @test all(eachindex(ids)) do i
                selectindices(globe, At(ids[i])) == selectindices(stored, At(ids[i])) == i
            end
            points = [center_lonlat(system, level, id) for id in ids]
            @test all(eachindex(ids)) do i
                selectindices(globe, Contains(points[i])) ==
                    selectindices(stored, Contains(points[i]))
            end
            # ...and the answer both give is the identity, which is what makes
            # the equivalence above worth having rather than two matching
            # wrongs: a cell's own center falls in that cell.
            @test all(i -> contained_position(
                selectindices(globe, Contains(points[i]))) == i, eachindex(ids))
        end
    end
end

@testset "globe point selectors work where no explicit lookup could" begin
    for (Lookup, system, _, huge) in SYSTEMS
        globe = Lookup(DGGSGlobeIds(system, huge))
        total = length(globe)
        @test total == DGG.num_cells(system, huge)
        # Sampled, necessarily — the levels here run to 1.1e17 cells. Both ends
        # are included because they are where an ordinal that is one off, or
        # counted from zero, stops being in range at all.
        for i in (1, 2, total ÷ 7, total ÷ 2, total - 1, total)
            id = globe[i]
            @test selectindices(globe, At(id)) == i
            @test contained_position(
                selectindices(globe, Contains(center_lonlat(system, huge, id)))) == i
        end
    end
end

@testset "an id naming no cell is rejected, not resolved to a neighbour" begin
    for (Lookup, system, levels, huge) in SYSTEMS
        level = last(levels)
        globe = Lookup(DGGSGlobeIds(system, level))
        stored = explicit_lookup(Lookup, collect(DGGSGlobeIds(system, level)), level)
        huge_globe = Lookup(DGGSGlobeIds(system, huge))
        ID = DGG.cell_id_type(system)

        for id in garbage_ids(ID, Int(DGG.num_cells(system, level)))
            # An `ArgumentError`, never an index: an id whose bits decode to an
            # in-range ordinal of some other cell is exactly what the round trip
            # exists to catch, and it is the failure mode that would otherwise
            # be silent.
            @test_throws ArgumentError selectindices(globe, At(id))
            # ...and the sentence is the partial lookup's own, since "this
            # lookup does not hold that cell" is the same fact on both paths.
            @test sprint(showerror, raised(() -> selectindices(globe, At(id)))) ==
                  sprint(showerror, raised(() -> selectindices(stored, At(id))))
        end
        # The same at a level with 5.7e14 cells and up, where the alternative —
        # an `is_valid_cell` pass over the lookup — could not be run at all.
        for id in garbage_ids(ID, Int(DGG.num_cells(system, huge)))
            @test_throws ArgumentError selectindices(huge_globe, At(id))
        end

        # Ids of the wrong level, which are the garbage a caller actually
        # produces: a resolution mixed up between two lookups of the same
        # system. Only the structural encodings can be asked — an ordinal id
        # carries no level of its own, so a level-`n + 1` HEALPix id in range at
        # level `n` *is* a level-`n` id, on the globe path and on the partial
        # one alike (`has_ordinal_ids`, `src/core/kernel.jl`).
        if !DGG.has_ordinal_ids(system)
            for other in (level - 1, level + 1)
                wrong = DGG.ordinal_to_cell(system, other, 3)
                @test_throws ArgumentError selectindices(globe, At(wrong))
                @test sprint(showerror, raised(() -> selectindices(globe, At(wrong)))) ==
                      sprint(showerror, raised(() -> selectindices(stored, At(wrong))))
            end
            @test_throws ArgumentError selectindices(
                huge_globe, At(DGG.ordinal_to_cell(system, huge - 1, 3)))
        end
    end
end

# The mechanism the two testsets above are the behaviour of. Both branches of
# one generic function, dispatched on the id vector rather than on the lookup —
# the lookup supertype is parameterized by its *element* type, so the globe case
# cannot be named on it (`src/core/lookups.jl`).
@testset "cell_position is the ordinal on a globe and a search on a Vector" begin
    for (_, system, levels, huge) in SYSTEMS
        ids = DGGSGlobeIds(system, huge)
        for i in (1, 2, length(ids) ÷ 2, length(ids))
            @test DGG.cell_position(ids, ids[i]) ==
                  DGG.cell_to_ordinal(system, huge, ids[i]) == i
        end
        for id in garbage_ids(DGG.cell_id_type(system), Int(DGG.num_cells(system, huge)))
            @test DGG.cell_position(ids, id) === nothing
        end

        # The stored branch is the one this step did not change: a binary search
        # over whatever ids are there, absent reported as `nothing`.
        materialized = collect(DGGSGlobeIds(system, first(levels)))
        @test DGG.cell_position(materialized, materialized[3]) == 3
        @test DGG.cell_position(materialized[2:end], materialized[1]) === nothing
        # A slice of a globe is a `Vector`, so it takes that branch and its
        # positions are the slice's, not the globe's (§1.4).
        @test DGG.cell_position(materialized[3:end], materialized[5]) == 3
    end
end

end # module GlobeSelectorTests
