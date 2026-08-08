using DiscreteGlobalGrids
using Documenter
using DocumenterVitepress
using Bonito
using GeoMakie
using Makie
using WGLMakie

WGLMakie.activate!()
# Makie.inline!(true)

makedocs(;
    modules = [DiscreteGlobalGrids],
    authors = "Anshul Singhvi and contributors",
    sitename = "DiscreteGlobalGrids.jl",
    repo = "https://github.com/asinghvi17/DiscreteGlobalGrids.jl",
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/asinghvi17/DiscreteGlobalGrids.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "DGGS gallery" => "all_dggs.md",
    ],
    plugins = [DocumenterVitepress.BonitoPlugin()],
    checkdocs = :none,
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "https://github.com/asinghvi17/DiscreteGlobalGrids.jl",
    devbranch = "main",
    push_preview = true,
)
