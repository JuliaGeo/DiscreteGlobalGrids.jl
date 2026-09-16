const FAST_DOCS = lowercase(get(ENV, "DGG_DOCS_FAST", "false")) in ("1", "true", "yes")

using DiscreteGlobalGrids
using Documenter
using Literate
using Zarr

if !FAST_DOCS
    using DocumenterVitepress
    using Bonito
    using GeoMakie
    using Makie
    using WGLMakie

    WGLMakie.activate!()
end

# Make unqualified package names available to every doctest.
DocMeta.setdocmeta!(DiscreteGlobalGrids, :DocTestSetup,
                    :(using DiscreteGlobalGrids); recursive = true)

const TUTORIALS = (
    "choosing_a_grid",
    "regridding",
    "stencils",
    "zonal",
    "multiorder",
    "between_grids",
    "hydrology",
    "store_io",
    "out_of_core",
    "healpix_astronomy",
)

# A fast build stages the source in a temporary directory. Literate output and
# HTML therefore never modify the working tree.
const DOCS_ROOT = FAST_DOCS ? mktempdir(; prefix = "dgg-docs-fast-") : @__DIR__
const DOCS_SOURCE = joinpath(DOCS_ROOT, "src")

if FAST_DOCS
    cp(joinpath(@__DIR__, "src"), DOCS_SOURCE; force = true)
end

# Literate's execute=false controls conversion only. DocumenterFlavor still
# emits @example blocks, which makedocs executes later.
for tutorial in TUTORIALS
    Literate.markdown(joinpath(DOCS_SOURCE, "tutorials", tutorial * ".jl"),
                      joinpath(DOCS_SOURCE, "tutorials");
                      flavor = Literate.DocumenterFlavor(), execute = false)
end

if FAST_DOCS
    # Keep every page in the build so references and @docs blocks are checked.
    # These pages contain the examples that need graphics, downloads, optional
    # development packages, or large datasets. Convert their @example blocks
    # to ordinary Julia code fences in the staged copy only.
    skipped_example_pages = [
        joinpath(DOCS_SOURCE, "index.md"),
        joinpath(DOCS_SOURCE, "all_dggs.md"),
        joinpath(DOCS_SOURCE, "extending.md"),
        (joinpath(DOCS_SOURCE, "tutorials", tutorial * ".md") for tutorial in TUTORIALS)...,
    ]
    for page in skipped_example_pages
        source = read(page, String)
        write(page, replace(source, r"```@example[^\n]*" => "```julia"))
    end
    @info "Fast docs mode: skipped executable examples" pages = relpath.(skipped_example_pages, DOCS_SOURCE)
end

# The regridding verbs and the DE9IM predicates are re-exported, so their
# docstrings belong to those packages; `GlobalRegridding` is not a direct
# dependency of this environment, so reach it through a name it owns.
const GlobalRegridding = parentmodule(DiscreteGlobalGrids.regrid)
const Trees = DiscreteGlobalGrids.Trees

const DOCS_MODULES = [
    DiscreteGlobalGrids,
    DiscreteGlobalGrids.Fallbacks,
    DiscreteGlobalGrids.Engine,
    DiscreteGlobalGrids.Encodings,
    DiscreteGlobalGrids.ChunkedLookups,
    DiscreteGlobalGrids.CellLookups,
    DiscreteGlobalGrids.Helpers,
    Base.get_extension(DiscreteGlobalGrids, :DiscreteGlobalGridsZarrExt),
    DiscreteGlobalGrids.DE9IM,
    GlobalRegridding,
    Trees,
]

const DOCS_PAGES = [
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
    "Grids and cell indices" => "abstractions.md",
    "Writing a grid system" => "extending.md",
    "API" => [
        "The grid interface" => "api/grid-interface.md",
        "Selecting cells" => "api/selecting-cells.md",
        "Regridding calls and plans" => "api/regridding.md",
        "Choosing a regridding method" => "api/regridding-methods.md",
        "Region boundaries" => "api/boundaries.md",
        "Neighbours and stencils" => "api/neighbors.md",
        "Reading and writing DGGS stores" => "api/store-io.md",
        "Sweeping a cube along its chunk lines" => "api/chunk-sweep.md",
        "Assigning chunks to workers" => "api/partitioning.md",
        "Requesting neighbour fields" => "api/neighbor-fields.md",
        "The ancestor-subzone layout" => "api/subzone-layout.md",
    ],
    "Internals" => [
        "Grid extension reference" => "internals/grid-contracts.md",
        "Collection and traversal contracts" => "internals/collection-contracts.md",
        "System capabilities and traversal costs" => "internals/system-capabilities.md",
        "Workflow execution details" => "internals/workflow-contracts.md",
        "Boundary traversal engines" => "internals/boundary-engines.md",
    ],
]

const COMMON_MAKEDOCS = (
    root = DOCS_ROOT,
    source = "src",
    build = "build",
    modules = DOCS_MODULES,
    authors = "Anshul Singhvi and contributors",
    sitename = "DiscreteGlobalGrids.jl",
    pages = DOCS_PAGES,
    checkdocs = :none,
)

if FAST_DOCS
    makedocs(;
        COMMON_MAKEDOCS...,
        format = Documenter.HTML(; prettyurls = false),
        remotes = nothing,
        warnonly = false,
    )
else
    makedocs(;
        COMMON_MAKEDOCS...,
        repo = Documenter.Remotes.GitHub("JuliaGeo", "DiscreteGlobalGrids.jl"),
        format = DocumenterVitepress.MarkdownVitepress(;
            repo = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
            devbranch = "main",
            devurl = "dev",
        ),
        plugins = [DocumenterVitepress.BonitoPlugin()],
        warnonly = true,
    )

    DocumenterVitepress.deploydocs(;
        # `deploydocs` expects host/user/repository syntax without a URL protocol.
        repo = "github.com/JuliaGeo/DiscreteGlobalGrids.jl.git",
        devbranch = "main",
        push_preview = true,
    )
end
