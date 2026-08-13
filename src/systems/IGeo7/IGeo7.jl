# ---------------------------------------------------------------------------
# T4 — IGeo7 (ISEA7H + Z7).
#
# The shared ISEA basis (icosahedron tables + Snyder charts) is included here
# rather than by the package module, because it belongs to the ISEA *family*
# and not to IGEO7: ISEA4R and friends will want the same module. Whichever
# system is included first defines it; the guard makes the second a no-op, so
# the include order of `src/systems/` never matters.
# ---------------------------------------------------------------------------

isdefined(@__MODULE__, :ISEA) || include("../ISEA/ISEA.jl")

"""
    IGeo7

The IGEO7 discrete global grid system: an aperture-7 hexagonal hierarchy on the
icosahedron (ISEA7H) with Z7 indexing.

The entry points are [`IGeo7System`](@ref), the singleton, and [`Z7Cell`](@ref),
its canonical cell id; everything else is reached through the package's generic
interface (`levelgrid`, `cellindex`, `cell_boundary`, `neighbors`, `parent`,
`children`, `descendant_range`, ...).

| file        | contents                                                      |
|:------------|:--------------------------------------------------------------|
| `z7.jl`     | the Z7 `UInt64` bit format, string/hex codecs, prefix ops     |
| `engine.jl` | Eisenstein integer arithmetic + the fitted digit tables       |
| `z7grid.jl` | encode/decode, areas, dense ordinals, subtree borders         |
| `system.jl` | `Z7Cell`, `IGeo7System`, `IGeo7Grid` — the interface wiring   |

The first three are ported verbatim from the verified clean-room
implementation, whose agreement with DGGRID is pinned by the sealed oracle
vectors: all 196,080 published cell centres at levels 1–5 reproduce to within
the dumps' own print noise and decode to their exact Z7 string. Only
`system.jl` is new, and it is adapter code — no projection maths is rederived.

`z7.jl` and `engine.jl` are pure integer layers; `z7grid.jl` composes them with
[`ISEA`](@ref)'s floating-point geometry.

# The subtree rim

The Z7 border automaton is a genuine fast path — `O(rim)` from the digits alone,
against the `O(subtree)` a neighbour sweep costs — and since T7 introduced the
generic hook it is reached as a method on [`subtree_border`](@ref) rather than
under a module-local name. `subtree_border_count` stays IGEO7's own: it answers
the size in closed form without a walk, which the generic interface has no
counterpart for.
"""
module IGeo7

import ..DiscreteGlobalGrids as DGG
# Only names that cannot collide with the native geometry layer are imported.
# `cell_boundary`, `cell_area` and `cell_centroid` are deliberately NOT among
# them: `z7grid.jl` defines native `(lon, lat)` functions of two of those names,
# and the interface generics are reached through `DGG.` so the two can never be
# conflated. See the header of `system.jl`.
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge, level, rawid,
    subtree_border
import ..Helpers
using ..ISEA

import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

const USPoint = GO.UnitSphericalPoint{Float64}

# Dependency order: integer id layer, integer lattice engine, then the geometry
# that composes them, then the interface wiring over all three.
include("z7.jl")
include("engine.jl")
include("z7grid.jl")
include("system.jl")

# The system's contract surface. The package-level exports happen in T7; these
# make the names reachable as `DiscreteGlobalGrids.IGeo7.<name>`.
export IGeo7System,
    IGeo7Grid,
    Z7Cell,
    InvalidZ7Error,
    MAX_RESOLUTION,
    equal_area_steradians,
    is_pentagon,
    is_valid_cell,
    subtree_border,
    subtree_border_count,
    z7_hex,
    z7_string

end # module IGeo7
