using Documenter, ForneyLab

makedocs(modules = [ForneyLab],
    clean = true,
    sitename = "ForneyLab.jl",
    assets = [
        joinpath("assets", "favicon.ico"),
    ],
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

deploydocs(
    repo = "https://github.com/biaslab/ForneyLab.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
