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

export syslabel, basesystem, isquadface, iscongruent, ishexwalk,
    hassortedsubtrees, forsystems, sweepcovers

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

"""
    iscongruent(sys) -> Bool

`sys` refines congruently: a cell's children tile its footprint exactly, so a
subtree and its root cover the same area.

The quad-face family has this, and so do the literature systems whose aperture
is a square on their own face shape — triangular ISEA4T, the rHEALPix
nonuple quads, and the IVEA/RTEA rhombi. It is a geometric fact per refinement
with no other trait standing for it, so it is declared here.
"""
iscongruent(sys) = isquadface(sys) || basesystem(sys) isa Union{
    DGG.ISEA4TSystem,DGG.RHEALPixSystem,DGG.AusPIXSystem,
    DGG.IVEA4RSystem,DGG.IVEA9RSystem,DGG.RTEA4RSystem,DGG.RTEA9RSystem}

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
    sweepcovers(sweep; except = ())

Assert that a suite's own sweep table — rows whose first element is the system —
reaches every registered system and at least one `AuthalicSystem` wrap, so that
a newly registered system fails here instead of going silently untested.

`except` names system types the suite's law does not apply to; each must be
absent from the sweep, so an exclusion cannot quietly shadow a swept row.
"""
function sweepcovers(sweep; except = ())
    swept = Set(typeof(first(row)) for row in sweep)
    for sys in DGG.systems()
        @test (typeof(sys) in swept) ⊻ (typeof(sys) in except)
    end
    @test any(row -> first(row) isa DGG.AuthalicSystem, sweep)
    return nothing
end

end # module DGGTestHelpers
