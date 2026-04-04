import Pkg

Pkg.develop(path = joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Documenter
using MeijerG

DocMeta.setdocmeta!(MeijerG, :DocTestSetup, :(using MeijerG); recursive=true)

makedocs(
    modules = [MeijerG],
    sitename = "MeijerG.jl",
    checkdocs = :exports,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true", assets = ["assets/theme.js"]),
    pages = [
        "Home" => "index.md",
        "Mathematical Notes" => "math.md",
        "Reductions" => "reductions.md",
        "Visual Examples" => "visual_examples.md",
        "API" => "api.md",
        "External Alternatives" => "external_alt.md",
    ],
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/TSGut/MeijerG.jl.git",
        devbranch = "main",
    )
else
    @info "Local docs build detected: skipping deploydocs()."
end
