using DiscreteGlobalGrids
using DiscreteGlobalGrids.IGeo7
using ArgParse

const _helpers = DiscreteGlobalGrids.Helpers
const _authalic = _helpers.WGS84_AUTHALIC
const _isea = DiscreteGlobalGrids.ISEA

# `vert0_lon_orientation` — the `Orientation` placing icosahedron vertex 0 at a
# chosen longitude (DGGRID's `dggs_vert0_lon`) — used to be defined here. It is
# now part of `IGeo7` itself, because the Zarr DGGS reader needs the same thing
# to honour an archive's declared placement; see its docstring.

"""
    centroid(z7_string; vert0_lon = ISEA.ISEA_LON0) -> (authalic, wgs84)

Centroid of the cell named by the Z7 string `z7_string` (e.g. `"0800433"`),
with icosahedron vertex 0 placed at longitude `vert0_lon`
(see [`vert0_lon_orientation`](@ref)).

Returns and prints two `(lon, lat)` pairs in degrees:
  1. `cell_center` as returned — on the authalic sphere the grid tiles;
  2. the same point with latitude converted authalic -> WGS84 geodetic.
Longitude is unchanged by the authalic transform.
"""
function centroid(z7_string::AbstractString; vert0_lon::Real = _isea.ISEA_LON0)
    z7 = z7_from_string(z7_string)
    orientation = vert0_lon_orientation(vert0_lon)
    lon, lat = cell_center(z7; orientation = orientation)
    wgs_lat = _helpers.authalic_to_geodeticd(_authalic, lat)
    println("cell ", z7_string, " (0x", z7_to_hex(z7), ")")
    println("  vertex 0 longitude:                    ", vert0_lon,
        vert0_lon == _isea.ISEA_LON0 ? " (standard ISEA)" : " (rotated placement)")
    println("  centroid as returned (authalic sphere): ", (lon, lat))
    println("  centroid authalic -> WGS84:            ", (lon, wgs_lat))
    return (lon, lat), (lon, wgs_lat)
end

"""
    parse_cli_args(argv) -> (z7_string, vert0_lon)

Parse the CLI arguments: a positional Z7 string cell id and the optional
`--vert0-lon` placement longitude.
"""
function parse_cli_args(argv)
    settings = ArgParseSettings(
        prog = basename(PROGRAM_FILE),
        description = "Print the centroid of a Z7 cell (authalic sphere and WGS84).",
    )
    @add_arg_table! settings begin
        "z7_string"
            help = "Z7 string cell id (e.g. \"0800433\")"
            required = true
        "--vert0-lon"
            help = "longitude of icosahedron vertex 0 in degrees (e.g. 11.2)"
            arg_type = Float64
            default = _isea.ISEA_LON0
    end
    args = ArgParse.parse_args(argv, settings)
    return args["z7_string"], args["vert0-lon"]
end

if abspath(PROGRAM_FILE) == @__FILE__
    z7_string, vert0_lon = parse_cli_args(ARGS)
    centroid(z7_string; vert0_lon = vert0_lon)
end
