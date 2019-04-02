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
            "Public API" => "library/public-api.md",
            "Internal API" => "library/internal-api.md"
        ],
        "Contributing" => "contributing.md",
        "Internals" => "internals.md"
    ]
)

deploydocs(
    repo = "https://github.com/biaslab/ForneyLab.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
