using DiscreteGlobalGrids
using Documenter
using DocumenterVitepress
using Bonito
using GeoMakie
using Makie
using WGLMakie
using Literate

WGLMakie.activate!()
# Makie.inline!(true)

# The tutorials are Literate.jl scripts; generate their markdown next to the
# sources, where the `pages` list below expects it. Execution happens in
# Documenter's @example blocks, not here.
for f in ("stencils", "zonal", "regridding", "hydrology", "healpix_astronomy")
    Literate.markdown(joinpath(@__DIR__, "src", "tutorials", f * ".jl"),
                      joinpath(@__DIR__, "src", "tutorials");
                      flavor = Literate.DocumenterFlavor(), execute = false)
end

makedocs(;
    modules = [DiscreteGlobalGrids],
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
            "Hydrology: a DEM on an IGEO7 grid" => "tutorials/hydrology.md",
            "The sky in HEALPix" => "tutorials/healpix_astronomy.md",
        ],
    ],
    plugins = [DocumenterVitepress.BonitoPlugin()],
    checkdocs = :none,
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    # Unlike `makedocs` and `MarkdownVitepress` above, which want a full URL,
    # `deploydocs` parses this as host/user/repo and rejects a protocol.
    repo = "github.com/JuliaGeo/DiscreteGlobalGrids.jl.git",
    devbranch = "main",
    push_preview = true,
)
