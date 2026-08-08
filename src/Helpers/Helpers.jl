"""
    Helpers

Shared allocation-free primitives used by every grid system: the inline
[`SmallList`](@ref) vector, sorted-ID lookup helpers, cell-id parsing, and the
[`AuthalicTransform`](@ref) ellipsoid ↔ authalic-sphere latitude conversion.
"""
module Helpers

export AuthalicTransform,
    SmallList,
    WGS84_AUTHALIC,
    authalic_q,
    authalic_radius,
    authalic_to_geodetic,
    authalic_to_geodeticd,
    empty_small_list,
    geodetic_to_authalic,
    geodetic_to_authalicd,
    small_push,
    small_sort,
    sorted_index,
    strictly_increasing,
    to_uint64_id,
    tuple_set

include("small_list.jl")
include("ids.jl")
include("authalic.jl")

end # module Helpers
