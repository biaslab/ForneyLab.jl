using Documenter, ForneyLab

makedocs(modules = [ForneyLab],
    clean = true,
    sitename = "ForneyLab.jl",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started/introduction.md",
        "User API" => Any[
            "Construct a model" => "user-API/construct-a-model.md",
            "Schedule an algorithm" => "user-API/schedule-an-algorithm.md",
            "Construct algorithm code" => "user-API/construct-algorithm-code.md",
            "Execute an algorithm" => "user-API/execute-an-algorithm.md",
            "Helpers" => "user-API/helpers.md",
        ],
        "Developer API" => Any[
            "Extended factor nodes" => "developer-API/extended-factor-nodes.md",
            "Extended rules" => "developer-API/extended-rules.md",
            "Graph (low-level)" => "developer-API/graph.md",
            "Scheduler (low-level)" => "developer-API/scheduler.md",
        ]

    ]
)

# deploydocs(
#     repo = "github.com/biaslab/ForneyLab.jl.git",
#     target = "build",
#     julia = "1.0.0",
#     deps = nothing,
#     make = nothing,
# )
