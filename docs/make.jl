using DiscreteGlobalGrids
using Documenter
using DocumenterVitepress
using Bonito
using GeoMakie
using Makie
using WGLMakie
using Literate

WGLMakie.activate!()

# Make unqualified package names available to every doctest.
DocMeta.setdocmeta!(DiscreteGlobalGrids, :DocTestSetup,
                    :(using DiscreteGlobalGrids); recursive = true)

# Generate Markdown beside each Literate source without executing its examples.
for f in ("stencils", "zonal", "regridding", "multiorder", "hydrology",
          "healpix_astronomy")
    Literate.markdown(joinpath(@__DIR__, "src", "tutorials", f * ".jl"),
                      joinpath(@__DIR__, "src", "tutorials");
                      flavor = Literate.DocumenterFlavor(), execute = false)
end

makedocs(;
    # Register `Fallbacks` so Documenter can render its boundary API docstrings.
    modules = [DiscreteGlobalGrids, DiscreteGlobalGrids.Fallbacks],
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
            "Stencil operations" => "tutorials/stencils.md",
            "Zonal statistics" => "tutorials/zonal.md",
            "Regridding a time series" => "tutorials/regridding.md",
            "Multi-order coverage" => "tutorials/multiorder.md",
            "Hydrology: a DEM on an IGEO7 grid" => "tutorials/hydrology.md",
            "The sky in HEALPix" => "tutorials/healpix_astronomy.md",
        ],
        "API" => [
            "Subtree and subset boundaries" => "api/boundaries.md",
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
