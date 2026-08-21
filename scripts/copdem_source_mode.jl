# Source-mode contract shared by the CopDEM production driver and its focused
# unit test. Keep this file dependency-free so the guard can be exercised
# without constructing grids, opening a store, or touching a tile directory.

"""
    effective_realspec(source, spec)

Return the local-real-tile selection for a CopDEM run.

Synthetic source is absolute: it accepts only `real = :none`. In particular,
`real = :auto` is an error rather than permission to replace synthetic tiles
with whatever GeoTIFFs happen to be cached locally. Real source retains the
driver's existing `:auto`, `:none`, and explicit-stem behavior.
"""
function effective_realspec(source, spec)
    source in (:real, :synthetic) || throw(ArgumentError(
        "source must be :real or :synthetic, got $(repr(source))"))
    source === :real && return spec
    spec === :none || throw(ArgumentError(
        "source=:synthetic requires real=:none; got real=$(repr(spec)). " *
        "Synthetic runs must never read cached real GeoTIFF overrides."))
    return :none
end
