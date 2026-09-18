# Tile names and the tile list: which 1x1-degree tiles exist at all.

"The public AWS Open Data bucket for GLO-`res` COGs."
bucketurl(res::Integer) = "https://copernicus-dem-$(res)m.s3.amazonaws.com"

"Nominal resolution in metres, 90 or 30."
resolution(sys::DGG.CopernicusDEMSystem) = CD.lat_intervals(sys) == 3600 ? 30 : 90

"The tile an AWS stem names, e.g. `Copernicus_DSM_COG_30_S90_00_E000_00_DEM`."
function stemtile(sys, stem::AbstractString)
    m = match(r"_([NS])(\d{2})_00_([EW])(\d{3})_00_DEM$", stem)
    m === nothing && throw(ArgumentError("$stem is not a Copernicus DEM stem"))
    lat = parse(Int, m[2]) * (m[1] == "S" ? -1 : 1)
    lon = parse(Int, m[4]) * (m[3] == "W" ? -1 : 1)
    return CD.tilecell(sys, lat, lon)
end

"The AWS object stem of `tile`. Inverse of [`stemtile`](@ref)."
function tilestem(sys, tile)
    lat, lon = CD.tilecorner(sys, tile)
    tag = resolution(sys) == 30 ? "10" : "30"
    return string("Copernicus_DSM_COG_", tag, "_", lat < 0 ? "S" : "N",
        lpad(abs(lat), 2, '0'), "_00_", lon < 0 ? "W" : "E",
        lpad(abs(lon), 3, '0'), "_00_DEM")
end

"""
    tilelist(datadir, res = 90; baseurl = bucketurl(res), timeout = 600) -> String

The local path of the bucket's `tileList.txt`, downloaded once into
`datadir/CopernicusDEM/` when absent.
"""
function tilelist(datadir::AbstractString, res::Integer = 90;
        baseurl::AbstractString = bucketurl(res), timeout::Real = 600)
    path = joinpath(datadir, "CopernicusDEM", "tileList-glo$(res).txt")
    isfile(path) && return path
    url = string(rstrip(String(baseurl), '/'), "/tileList.txt")
    mkpath(dirname(path))
    part, io = mktemp(dirname(path))
    close(io)
    @info "downloading the Copernicus tile list" url path
    try
        Downloads.download(url, part; timeout)
        filesize(part) > 0 || error("downloaded an empty object from $url")
        Base.Filesystem.rename(part, path)
    catch err
        isfile(path) && return path  # another process completed the same download
        throw(ArgumentError("no tile list at $path and downloading $url failed: " *
                            sprint(showerror, err)))
    finally
        rm(part; force = true)
    end
    return path
end

inregions(::Nothing, lon, lat) = true
inregions(regions, lon, lat) =
    any(r -> r[1] <= lon <= r[2] && r[3] <= lat <= r[4], regions)

"""
    listedtiles(sys, path, regions = nothing) -> Vector{Int}

The level-0 ordinals of every tile named in the tile list at `path`, ascending.
`regions` is `nothing` for the globe or a list of `(w, e, s, n)` boxes that a
tile's south-west corner must fall in.

An unlisted tile is open ocean: it has no object in the bucket, so a source grid
built over this list has no chunk there.
"""
function listedtiles(sys, path::AbstractString, regions = nothing)
    isfile(path) || throw(ArgumentError("no tile list at $path"))
    out = Int[]
    for line in eachline(path)
        stem = strip(line)
        isempty(stem) && continue
        t = stemtile(sys, stem)
        lat, lon = CD.tilecorner(sys, t)
        inregions(regions, lon, lat) || continue
        push!(out, Int(t.index))
    end
    return sort!(unique!(out))
end
