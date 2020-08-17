using Documenter, ForneyLab

makedocs(modules = [ForneyLab],
    clean = true,
    sitename = "ForneyLab.jl",
    assets = [
        joinpath("assets", "favicon.ico"),
    ],
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "User guide" => "user-guide.md",
        "Library" => Any[
            "User API" => "library/user-api.md",
            "Developer API" => "library/developer-api.md"
        ],
        "Contributing" => "contributing.md",
        "Internals" => "internals.md"
    ],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/biaslab/ForneyLab.jl.git"
    )
end