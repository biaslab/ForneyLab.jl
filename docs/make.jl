using Documenter, ForneyLab

makedocs(modules = [ForneyLab],
    clean = true,
    sitename = "ForneyLab.jl",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started/introduction.md",
        "API" => Any[
            "Algorithms" => "API/algorithms.md",
            "Factor nodes" => "API/factor_nodes.md",
            "Engines" => "API/engines.md",
            "General" => "API/general.md"
        ],
        "Internals" => "internals/internals.md",
    ]
)

# deploydocs(
#     repo = "github.com/biaslab/ForneyLab.jl.git",
#     target = "build",
#     julia = "1.0.0",
#     deps = nothing,
#     make = nothing,
# )
