# https://documenter.juliadocs.org/stable/man/guide/ 

using Documenter, FastGaussianProcesses

assetsdir = joinpath(@__DIR__,"src/assets")
if ~isdir(assetsdir) mkdir(assetsdir) end

println("DOCTEST")
makedocs(
    doctest = :only,
    strict = true
)

println("BUILDING DOCS")
makedocs(
    sitename = "FastGaussianProcesses.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true
        ), 
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
    ],
    doctest = false
)

println("DEPLOY DOCS")
deploydocs(
    repo = "github.com/alegresor/FastGaussianProcesses.jl.git",
    #devbranch = "main",
    #versions = ["stable" => "v^", "main" => "main"]
)
