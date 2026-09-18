# Tile sources. A source answers `loadtile(source, ordinal) -> Vector{Float32}`:
# one tile's posts in the tile's own index order, `NaN32` where there is no
# elevation. `CopernicusTiles` reads real COGs; `SyntheticTiles` (synthetic.jl)
# fabricates them.

"""
    loadtile(source, ordinal) -> Vector{Float32}

The posts of level-0 tile `ordinal`: raster row north first, column west to
east, column fastest. `NaN32` marks nodata.
"""
function loadtile end

"""
    CopernicusTiles(sys, listed; cachedir, download = true, baseurl, retries, backoff, timeout)

Real Copernicus DEM GeoTIFFs, kept in `cachedir` as `<stem>.tif`.

Construction touches neither disk nor network. With `download = true` the first
access to a listed tile fetches its COG from `baseurl` to `<stem>.tif.part` and
renames it into place, so only complete files carry the final name. One lock per
tile makes concurrent requests for the same tile a single GET. With
`download = false` a tile missing from `cachedir` is an error.

An unlisted tile is ocean: [`loadtile`](@ref) returns all-`NaN32` for it with no
filesystem or network access.

`ndownloads` counts successful GETs; `ncold` counts demands that had to start
one.
"""
struct CopernicusTiles{S<:DGG.CopernicusDEMSystem}
    sys::S
    listed::Set{Int}
    cachedir::String
    download::Bool
    baseurl::String
    retries::Int
    backoff::Float64
    timeout::Float64
    locks::Dict{Int,ReentrantLock}
    downloader::Downloads.Downloader
    ndownloads::Threads.Atomic{Int}
    ncold::Threads.Atomic{Int}
end

function CopernicusTiles(sys::DGG.CopernicusDEMSystem, listed;
        cachedir::AbstractString, download::Bool = true,
        baseurl::AbstractString = bucketurl(resolution(sys)),
        retries::Integer = 4, backoff::Real = 1.0, timeout::Real = 600.0)
    retries >= 1 || throw(ArgumentError("retries must be at least 1"))
    backoff >= 0 || throw(ArgumentError("backoff must be non-negative"))
    timeout > 0 || throw(ArgumentError("timeout must be positive"))
    tiles = Set(Int.(listed))
    return CopernicusTiles(sys, tiles, abspath(String(cachedir)), download,
        String(rstrip(String(baseurl), '/')), Int(retries), Float64(backoff),
        Float64(timeout), Dict(t => ReentrantLock() for t in tiles),
        Downloads.Downloader(), Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

"The local path of a tile; no filesystem or network access."
tilecachepath(source::CopernicusTiles, ordinal::Int) =
    joinpath(source.cachedir, tilestem(source.sys, DGG.LevelIndex(0, ordinal)) * ".tif")

"The public S3 URL of a tile."
function tileurl(source::CopernicusTiles, ordinal::Int)
    stem = tilestem(source.sys, DGG.LevelIndex(0, ordinal))
    return string(source.baseurl, "/", stem, "/", stem, ".tif")
end

_transientstatus(status::Integer) =
    status == 0 || status == 408 || status == 429 || 500 <= status < 600

"""
    tilepath!(source, ordinal; demand = true) -> Union{String,Nothing}

The local path of a listed tile, downloading it first if needed; `nothing` for
an unlisted (ocean) tile. `demand = false` marks a speculative fetch, which
`source.ncold` leaves uncounted.
"""
function tilepath!(source::CopernicusTiles, ordinal::Int; demand::Bool = true)
    ordinal in source.listed || return nothing
    path = tilecachepath(source, ordinal)
    isfile(path) && return path
    source.download || error("no GeoTIFF at $path and downloading is off; " *
                             "fetch it from $(tileurl(source, ordinal))")
    return lock(source.locks[ordinal]) do
        isfile(path) && return path
        demand && Threads.atomic_add!(source.ncold, 1)
        mkpath(source.cachedir)
        part = path * ".part"
        url = tileurl(source, ordinal)
        last_error = nothing
        for attempt in 1:source.retries
            try
                Downloads.download(url, part; timeout = source.timeout,
                    downloader = source.downloader)
                filesize(part) > 0 || error("downloaded an empty object from $url")
                Base.Filesystem.rename(part, path)
                Threads.atomic_add!(source.ndownloads, 1)
                return path
            catch err
                last_error = err
                status = err isa Downloads.RequestError ? err.response.status : 0
                rm(part; force = true)
                (status == 403 || status == 404) &&
                    error("listed Copernicus tile returned HTTP $status from $url")
                (!_transientstatus(status) || attempt == source.retries) && break
                delay = source.backoff * 2.0^(attempt - 1)
                @warn "Copernicus tile download failed; retrying" url attempt delay status
                sleep(delay)
            end
        end
        error("failed to download $url after $(source.retries) attempt(s): " *
              sprint(showerror, last_error))
    end
end

# ArchGDAL states no thread-safety guarantee and caches load outside their locks.
const GDALLOCK = ReentrantLock()

"The decoded band of the GeoTIFF at `path`, size-checked and flattened into index order."
function readtile(path, sys, tile::DGG.LevelIndex)
    band = lock(GDALLOCK) do
        ArchGDAL.read(ds -> ArchGDAL.read(ds, 1), path)
    end
    lat, _ = CD.tilecorner(sys, tile)
    expected = (Int(CD.ncols_at(sys, lat)), Int(CD.lat_intervals(sys)))
    size(band) == expected || error(
        "Copernicus tile at $path has raster size $(size(band)); expected $expected")
    return Float32.(vec(band))
end

function loadtile(source::CopernicusTiles, ordinal::Int)
    path = tilepath!(source, ordinal)
    tile = DGG.LevelIndex(0, ordinal)
    if path === nothing
        lat, _ = CD.tilecorner(source.sys, tile)
        return fill(NaN32, Int(CD.ncols_at(source.sys, lat) * CD.lat_intervals(source.sys)))
    end
    return readtile(path, source.sys, tile)
end
