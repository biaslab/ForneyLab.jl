using Documenter, ForneyLab

makedocs(modules = [ForneyLab],
    clean = true,
    sitename = "ForneyLab.jl",
    assets = [
        joinpath("assets", "favicon.ico"),
    ],
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started/introduction.md",
        "Library" => Any[
            "User API" => "library/user-api.md",
            "Internals API" => "library/internals-api.md"
        ],
        "Contributing" => "contributing.md",
    ]
)

deploydocs(
    repo = "https://github.com/biaslab/ForneyLab.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
