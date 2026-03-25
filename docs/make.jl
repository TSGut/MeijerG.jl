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
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Reductions" => "reductions.md",
        "Mathematical Notes" => "math.md",
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
