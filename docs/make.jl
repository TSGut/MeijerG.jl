using Documenter
using MeijerG

DocMeta.setdocmeta!(MeijerG, :DocTestSetup, :(using MeijerG); recursive=true)

makedocs(
    modules = [MeijerG],
    sitename = "MeijerG.jl",
    checkdocs = :exports,
    # ConfluentJets is an internal submodule (its exports are not part of the
    # MeijerG API) so its docstrings are not part of the manual.
    checkdocs_ignored_modules = [MeijerG.ConfluentJets],
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
