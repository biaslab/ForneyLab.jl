using Documenter, ForneyLab

ENV["GKS_ENCODING"] = "utf-8"
ENV["GKSwstype"] = "100"

makedocs(modules = [ForneyLab],
    clean = true,
    sitename = "ForneyLab.jl",
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
    format = Documenter.HTML(
        assets = [
            joinpath("assets", "favicon.ico"),
        ],
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/biaslab/ForneyLab.jl.git"
    )
end