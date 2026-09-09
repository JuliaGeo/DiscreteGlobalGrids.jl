using DiscreteGlobalGrids
using Documenter
using DocumenterVitepress
using Bonito
using GeoMakie
using Makie
using WGLMakie
using Literate
using Zarr

WGLMakie.activate!()

# Make unqualified package names available to every doctest.
DocMeta.setdocmeta!(DiscreteGlobalGrids, :DocTestSetup,
                    :(using DiscreteGlobalGrids); recursive = true)

# Generate Markdown beside each Literate source without executing its examples.
for f in ("choosing_a_grid", "regridding", "stencils", "zonal", "multiorder",
          "between_grids", "hydrology", "store_io", "out_of_core",
          "healpix_astronomy")
    Literate.markdown(joinpath(@__DIR__, "src", "tutorials", f * ".jl"),
                      joinpath(@__DIR__, "src", "tutorials");
                      flavor = Literate.DocumenterFlavor(), execute = false)
end

# The regridding verbs and the DE9IM predicates are re-exported, so their
# docstrings belong to those packages; `GlobalRegridding` is not a direct
# dependency of this environment, so reach it through a name it owns.
const GlobalRegridding = parentmodule(DiscreteGlobalGrids.regrid)
const Trees = DiscreteGlobalGrids.Trees

makedocs(;
    # Register `Fallbacks` and `Engine` so Documenter can render their boundary
    # API docstrings, and `Encodings`/`ChunkedLookups` for the store-IO ones: the
    # main module re-binds those names, but the docstrings belong to the
    # submodules.
    modules = [DiscreteGlobalGrids, DiscreteGlobalGrids.Fallbacks,
               DiscreteGlobalGrids.Engine,
               DiscreteGlobalGrids.Encodings, DiscreteGlobalGrids.ChunkedLookups,
               # Same reason, for the grid-interface and selection pages:
               # `CellLookups` owns the DimensionalData layer and `Helpers` the
               # authalic transform.
               DiscreteGlobalGrids.CellLookups, DiscreteGlobalGrids.Helpers,
               # Documenter filters docstrings by module, so the Zarr
               # extension's dggread/dggwrite methods render only if it is
               # listed here.
               Base.get_extension(DiscreteGlobalGrids, :DiscreteGlobalGridsZarrExt),
               # Re-exported names whose docstrings this site renders:
               # `Intersects` and friends, `Weighted`, `ncells`/`treeify`.
               DiscreteGlobalGrids.DE9IM, GlobalRegridding, Trees],
    authors = "Anshul Singhvi and contributors",
    sitename = "DiscreteGlobalGrids.jl",
    repo = Documenter.Remotes.GitHub("JuliaGeo", "DiscreteGlobalGrids.jl"),
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "DGGS gallery" => "all_dggs.md",
        "Tutorials" => [
            "Choosing a grid" => "tutorials/choosing_a_grid.md",
            "Regridding: getting data onto a grid" => "tutorials/regridding.md",
            "Moving between DGGS" => "tutorials/between_grids.md",
            "Zonal statistics" => "tutorials/zonal.md",
            "Multi-order coverage" => "tutorials/multiorder.md",
            "Stencil operations" => "tutorials/stencils.md",
            "Hydrology: a DEM on an IGEO7 grid" => "tutorials/hydrology.md",
            "A round trip through a DGGS store" => "tutorials/store_io.md",
            "Out of core" => "tutorials/out_of_core.md",
            "The sky in HEALPix" => "tutorials/healpix_astronomy.md",
        ],
        "Writing a grid system" => "extending.md",
        "API" => [
            "The grid interface" => "api/grid-interface.md",
            "Selecting cells" => "api/selecting-cells.md",
            "Choosing a regridding method" => "api/regridding-methods.md",
            "Region boundaries" => "api/boundaries.md",
            "Neighbours and stencils" => "api/neighbors.md",
            "Reading and writing DGGS stores" => "api/store-io.md",
            "Sweeping a cube along its chunk lines" => "api/chunk-sweep.md",
            "Requesting neighbour fields" => "api/neighbor-fields.md",
            "The ancestor-subzone layout" => "api/subzone-layout.md",
        ],
        "Internals" => [
            "Boundary traversal engines" => "internals/boundary-engines.md",
        ],
    ],
    plugins = [DocumenterVitepress.BonitoPlugin()],
    checkdocs = :none,
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    # `deploydocs` expects host/user/repository syntax without a URL protocol.
    repo = "github.com/JuliaGeo/DiscreteGlobalGrids.jl.git",
    devbranch = "main",
    push_preview = true,
)
