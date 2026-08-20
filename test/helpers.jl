# Selection and naming shared by the cross-system sweeps.
#
# Each suite is its own module and stays runnable on its own, so this file is
# `include`d by the suites that need it rather than by `runtests.jl`; the
# include defines `DGGTestHelpers` inside the including module.
#
# Which systems a law applies to is a question about traits, never about names:
# `AbstractQuadFaceGridSystem` is the radix-4 family, `hex_child_direction` is
# the hex directed-walk hooks, `has_sorted_subtrees` is the contiguous-subtree
# arithmetic that A5 alone lacks. Per-system *expectations* — budget tables,
# measured counts, oracle routes — stay keyed by name where they are stated.

module DGGTestHelpers

using Test
import DiscreteGlobalGrids as DGG

export syslabel, basesystem, isquadface, ishexwalk, hassortedsubtrees,
    forsystems, sweepcovers

"""
    syslabel(sys) -> String

Testset name for a grid system: its type name, with `AuthalicSystem` shown
wrapping its parent's label.
"""
syslabel(sys) = sys isa DGG.AuthalicSystem ?
                "Authalic($(syslabel(parent(sys))))" : string(nameof(typeof(sys)))

"""
    basesystem(sys)

The system whose hierarchy `sys` uses. `AuthalicSystem` remaps geometry and
leaves identity, parentage and neighbourhood alone, so every trait below asks
about this.
"""
basesystem(sys) = sys isa DGG.AuthalicSystem ? basesystem(parent(sys)) : sys

"""
    isquadface(sys) -> Bool

`sys` is one of the radix-4 quad-face family: four children that tile their
parent exactly, on an aligned per-face lattice.
"""
isquadface(sys) = basesystem(sys) isa DGG.AbstractQuadFaceGridSystem

"`sys` supplies the hex hooks the calibrated directed halo walk needs."
function ishexwalk(sys)
    b = basesystem(sys)
    return hasmethod(DGG.hex_child_direction,
                     Tuple{typeof(b),DGG.cellindextype(b)})
end

"`sys` names a subtree as a contiguous ascending id range."
hassortedsubtrees(sys) = DGG.has_sorted_subtrees(basesystem(sys))

"""
    forsystems(; quadface, hexwalk, sortedsubtrees) -> Tuple
    forsystems(f; kwargs...)

The registered systems matching every trait given, in registry order; each
keyword left out is not constrained. The one-argument form applies `f` to each.

    @testset "\$(syslabel(sys))" for sys in forsystems(quadface = true)
"""
function forsystems(; quadface = nothing, hexwalk = nothing,
                    sortedsubtrees = nothing)
    matches(sys) =
        (quadface === nothing || isquadface(sys) == quadface) &&
        (hexwalk === nothing || ishexwalk(sys) == hexwalk) &&
        (sortedsubtrees === nothing || hassortedsubtrees(sys) == sortedsubtrees)
    return filter(matches, DGG.systems())
end

forsystems(f; kwargs...) = foreach(f, forsystems(; kwargs...))

"""
    sweepcovers(sweep)

Assert that a suite's own sweep table — rows whose first element is the system —
reaches every registered system and at least one `AuthalicSystem` wrap, so that
a newly registered system fails here instead of going silently untested.
"""
function sweepcovers(sweep)
    swept = Set(typeof(first(row)) for row in sweep)
    for sys in DGG.systems()
        @test typeof(sys) in swept
    end
    @test any(row -> first(row) isa DGG.AuthalicSystem, sweep)
    return nothing
end

end # module DGGTestHelpers
