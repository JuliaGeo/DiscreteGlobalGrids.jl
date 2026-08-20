# The one error type of the store-IO layer, defined layer-neutrally so that
# every layer can throw it.
#
# INCLUDE THIS FIRST. `encodings.jl` and `chunked_lookup.jl` define submodules
# that `import ..DGGSFormatError` from whichever module includes them, so the
# type must already be bound there. Nothing in this file depends on anything
# else in `src/io/`.

"""
    DGGSFormatError(; check, store = nothing, conventions = nothing,
                    declared = nothing, observed = nothing, detail = "")

Raised when a store cannot be read as it says it can: conventions that
contradict each other, a grid name in no registry, an alias supplied twice, a
length that does not check out, an id that names no cell.

  - `check`: a symbol naming the check that failed, e.g. `:level_disagreement`
    or `:unknown_grid_name`. This is the field to test against.
  - `store`: the store identifier, so an error from a lazy read still says
    where it came from.
  - `conventions`: the conventions that fired, in the order they were tried.
  - `declared` / `observed`: the two values that failed to reconcile.
  - `detail`: one sentence saying what to do about it.

`store` and `conventions` are CONTEXT, and optional: the encoding and lookup
layers see ids and lengths, not stores, and throw without them. The layer that
does know the store adds it with [`with_store_context`](@ref); `showerror` omits
whatever is absent.

The policy this serves: vocabulary disagreement is refused rather than guessed.
"""
Base.@kwdef struct DGGSFormatError <: Exception
    check::Symbol
    store::Union{String,Nothing} = nothing
    conventions::Union{Vector{String},Nothing} = nothing
    declared::Any = nothing
    observed::Any = nothing
    detail::String = ""
end

_hascontext(::Nothing) = false
_hascontext(x::AbstractString) = !isempty(x)
_hascontext(x::AbstractVector) = !isempty(x)

# A reconciled value can be a whole id array. Show enough of it to recognise it.
function _brief(x)
    s = repr(x)
    return length(s) <= 120 ? s : first(s, 117) * "…"
end

function Base.showerror(io::IO, e::DGGSFormatError)
    print(io, "DGGSFormatError: ", e.check)
    _hascontext(e.store) && print(io, "\n  store: ", e.store)
    _hascontext(e.conventions) &&
        print(io, "\n  conventions fired: ", join(e.conventions, ", "))
    e.declared === nothing || print(io, "\n  declared: ", _brief(e.declared))
    e.observed === nothing || print(io, "\n  observed: ", _brief(e.observed))
    isempty(e.detail) || print(io, "\n  ", e.detail)
    return nothing
end

"""
    store_context(e::DGGSFormatError; store = nothing, conventions = nothing)

`e` with `store` and `conventions` filled in where it carries none. Context
already on the error wins: the throw site knew more than the boundary does.
"""
store_context(e::DGGSFormatError; store=nothing, conventions=nothing) =
    DGGSFormatError(check=e.check,
        store=_hascontext(e.store) ? e.store :
              store === nothing ? nothing : String(store),
        conventions=_hascontext(e.conventions) ? e.conventions : conventions,
        declared=e.declared, observed=e.observed, detail=e.detail)

"""
    with_store_context(f, store; conventions = nothing)

Run `f`, and rethrow any [`DGGSFormatError`](@ref) it raises with `store` and
`conventions` added.

This is how a store-aware caller — the Zarr extension around an encoding's
validation pass — turns a context-free error from a lower layer into one that
names the store it came from, without that layer ever learning what a store is.
"""
function with_store_context(f, store; conventions=nothing)
    try
        return f()
    catch e
        e isa DGGSFormatError || rethrow()
        throw(store_context(e; store, conventions))
    end
end
