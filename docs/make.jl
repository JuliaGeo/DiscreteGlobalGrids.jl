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
    repo = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
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
    # Unlike `makedocs` and `MarkdownVitepress` above, which want a full URL,
    # `deploydocs` parses this as host/user/repo and rejects a protocol.
    repo = "github.com/JuliaGeo/DiscreteGlobalGrids.jl.git",
    devbranch = "main",
    push_preview = true,
)
