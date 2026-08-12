"""
    IGeo7

Clean-room implementation of the IGEO7 DGGS: an aperture-7 hexagonal
hierarchy on the icosahedron (ISEA7H) with Z7 indexing, following
`docs/IGeo7/spec/design.md`.

Provenance: every algorithmic choice is sourced from a published paper
(`spec/z7-paper-spec.md`, `spec/isea-projection-spec.md`,
`spec/aperture7-indexing-spec.md`), from first-principles mathematics, or is
**fitted** to recorded black-box oracle output; every such fit, and the
evidence behind it, is recorded in `spec/igeo7-geometry-diagnosis.md`.

The spherical icosahedron and the Snyder equal-area charts are not here: they
are shared with the rest of the ISEA family and live in [`ISEA`](@ref), from
which this module takes `VERTICES`, `Orientation`, `snyder_fwd`, `dev_to_xyz`
and friends. What remains is the IGEO7-specific layering (`include` order below
is the dependency order):

| file        | contents                                                     |
|:------------|:-------------------------------------------------------------|
| `z7.jl`     | Z7 `UInt64` bit format, string/hex, prefix ops (no geometry)  |
| `engine.jl` | Eisenstein integer arithmetic + fitted digit tables           |
| `grid.jl`   | encode/decode, areas, dense indexing, subtree borders          |

plus the `IGeo7Lookups` integration module (`DimensionalData`) and
`IGeo7Kernel.jl`, which wires the package's operations kernel — and through it
the generic `GeometryOps` / `ConservativeRegridding` tree family — for
`IGEO7DGGS`.

`z7.jl` and `engine.jl` are integer code; `grid.jl` composes them with `ISEA`'s
floating-point geometry. The grid is the IGEO7/ISEA7H grid — an integer
Eisenstein lattice in the per-face Snyder ISEA plane, standard ISEA placement
(conventions fitted in `spec/igeo7-geometry-diagnosis.md`); agreement with the
oracle is verified at 100% exact z7 decode on all 196,080 oracle cell centers
res 1–5. The native core is stdlib + [`Helpers`](@ref) only.
"""
module IGeo7

import ..Helpers
using ..ISEA

include("z7.jl")
include("engine.jl")
include("grid.jl")
include("IGeo7Lookups.jl")

# Public API (spec/interface-contract.md). Names are defined by the include
# files listed above; the list is the module's contract surface. The shared
# icosahedron/Snyder names reach users through `ISEA`, not from here.
export InvalidZ7Error,
    MAX_RESOLUTION,
    border_descendants,
    cell_area,
    cell_boundary,
    cell_boundary_cartesian,
    cell_center,
    cell_to_children,
    cell_to_index,
    cell_to_parent,
    cell_to_z7,
    get_resolution,
    index_to_cell,
    is_pentagon,
    is_valid_cell,
    is_valid_z7,
    lonlat_to_cell,
    lonlat_to_index,
    lonlat_to_z7,
    num_cells,
    res0_cells,
    z7_base_cell,
    z7_child,
    z7_children,
    z7_digit,
    z7_from_hex,
    z7_from_string,
    z7_is_descendant,
    z7_is_pentagon,
    z7_parent,
    z7_resolution,
    z7_to_cell,
    z7_to_hex,
    z7_to_string

# Operations-kernel wiring for the `IGEO7DGGS` singleton. It defines methods
# on the *package's* generics only — it adds no name to this module's contract
# surface above, which is why it is included after the export block.
include("IGeo7Kernel.jl")

end # module IGeo7
